### functions  ####
odgmm = function(X,n_max_comp,debug=FALSE){
    # X: data, a 1-dimensional vector
    # n_max_comp: maximum number of components
    
    # perform overdispersed Gaussian mixture modeling on input data X, with maximum K components
    
    # functions for transforming parameters
    alpha2mu = function(alpha){
        # global var mu0
        return ( cumsum(c(mu0,alpha)) )
    }
    mu2alpha = function(mus){
        return ( diff(mus) )
    }
    beta2sigma = function(beta){
        # beta is for the variance
        # global var sigma0
        return( cumprod( c(sigma0,sqrt(beta)) ) )
    }
    beta2var = function(beta){
        # global var sigma0
        return( cumprod(c(sigma0^2,beta)) )
    }
    sigma2beta = function(sigmas){
        return (exp(diff(log(sigmas)))^2)
    }
    var2beta = function(vars){
        return (exp(diff(log(vars))))
    }
    
    # functions for EM algorithms
    estep = function(X, alpha, beta, ws){
        mus = alpha2mu(alpha)
        sds = beta2sigma(beta)
        K = length(alpha)
        N = length(X)
        Z = matrix(0,nrow=N,ncol=K+2)
        for(i in seq(K+1))
        {
            Z[,i] = log(ws[i]) + dnorm(X, mean = mus[i], sd = sds[i], log = TRUE)
        }
        
        if(any(is.na(Z))){
            print('Z contains NA')
        }
        
        # if(!is.na(ws[K+2]) & ws[K+2]>0){
        #   Z[uni_valid_inds,K+2] = log(ws[K+2]) + uni_log_density
        # }else{
        #   print('ws[K+2] is 0. Quit.')
        # }
        # Z[!uni_valid_inds,K+2] = 0
        Z[,K+2] = log(ws[K+2]) + uni_log_density
        
        
        if(any(is.na(Z))){
            print('Z contains NA')
        }
        
        
        maxz = apply(Z,1,max)
        Z = exp(Z - maxz)
        return (Z/apply(Z,1,sum))
    }
    cal_alpha = function(X,Z,alpha,beta,j){
        mus = alpha2mu(alpha)
        vars = beta2var(beta)
        K = length(alpha)
        N = length(X)
        XX = rep(0,N)
        ZZ = rep(0,N)
        for(k in seq(j,K))
        {
            tmp = Z[,(k+1)]/vars[k+1]
            XX = XX + (X - mus[k+1] + alpha[j])*tmp
            ZZ = ZZ + tmp
        }
        tsum = sum(ZZ)  # vars can be 0
        if(is.na(tsum) || tsum==0){
            return(NA)
        }else{
            return (sum(XX)/tsum)
        }
    }
    cal_beta = function(X,Z,alpha,beta,j){
        mus = alpha2mu(alpha)
        vars = beta2var(beta)
        K = length(alpha)
        N = length(X)
        XX = rep(0,N)
        for(k in seq(j,K))
        {
            XX = XX + Z[,(k+1)] * (X - mus[k+1])^2 * beta[j] / vars[k+1]
        }
        tsum = sum(Z[,seq(j+1,K+1)])
        if(is.na(tsum) || tsum==0){
            return(NA)
        }else{
            return (sum(XX)/tsum)
        }
    }
    cal_ws = function(Z){
        N = dim(Z)[1]
        ws = apply(Z,2,sum)/N
        return (ws)
    }
    mstep = function(X,Z,alpha,beta,ws){
        K = length(alpha)
        new_ws = cal_ws(Z)
        new_alpha = alpha
        new_beta = beta
        for(j in seq(K,1,-1))
        {
            new_alpha[j] = cal_alpha(X,Z,alpha,beta,j)
            new_beta[j] = cal_beta(X,Z,alpha,beta,j)
        }
        return(list(alpha=new_alpha,beta=new_beta,ws=new_ws))
    }
    exp_log_like = function(X, Z, alpha, beta, ws){
        mus = alpha2mu(alpha)
        sds = beta2sigma(beta)
        K = length(alpha)
        N = length(X)
        XX = rep(0,N)
        for(i in seq(K+1))
        {
            if(ws[i]==0){
                next
            }
            inds = Z[,i]!=0
            XX[inds] = XX[inds] + Z[inds,i] * (log(ws[i])+dnorm(X[inds], mean = mus[i], sd = sds[i], log = TRUE))
        }
        XX = XX  + Z[,K+2] * (log(ws[K+2])+uni_log_density)
        #XX[uni_valid_inds] = XX[uni_valid_inds]  + Z[uni_valid_inds,K+2] * (log(ws[K+2])+uni_log_density)
        
        # if(sum(XX)>0){
        #   print('expected lik greater than 0')
        # }
        return (sum(XX))
    }
    elbo = function(X,Z,alpha,beta, ws){
        LZ = Z
        LZ[Z!=0] = log(Z[Z!=0])
        entropy = -1 * Z * LZ
        lb = exp_log_like(X,Z,alpha,beta, ws)+sum(entropy)
        return (lb)
    }
    bic = function(X,Z,alpha,beta,ws){
        N = length(X)
        K = length(alpha)
        res = -2*elbo(X,Z,alpha,beta,ws) + (3*K+2)*log(N)
    }
    # generate the result
    gen_res = function(alpha=NA,beta=NA,ws=NA,logml=NA,Z=NA,bic=NA){
        if(any(is.na(alpha))){
            mus = NA
        }else{
            mus = alpha2mu(alpha)
        }
        if(any(is.na(beta))){
            sds = NA
        }else{
            sds = beta2sigma(beta)
        }
        return(list(alpha=alpha,beta=beta,ws=ws,logml=logml,Z=Z,bic=bic,mus=mus,sigmas=sds))
    }
    em_optim0 = function(X, init_mus, init_sgs, init_ws, nround=200, debug=FALSE){
        alpha = mu2alpha(init_mus)
        beta = var2beta(init_sgs)
        ws = init_ws
        Z = estep(X, alpha, beta, ws)
        logml = -Inf
        logml_res = rep(NA,nround)
        for (i in seq(nround)){
            if(debug){
                print(paste0('iteration=',i,'  logml=',logml))
            }
            Z = estep(X, alpha, beta, ws)
            res = mstep(X,Z,alpha,beta,ws)
            
            if(any(is.na(res$alpha)) || any(is.na(res$beta)) || any(is.na(res$ws))){
                print(paste0('Inference failed after ',i,' iterations.'))
                return(gen_res())
            }
            logml_new = elbo(X, Z, res$alpha, res$beta, res$ws)
            logml_res[i] = logml_new
            alpha = res$alpha
            beta = res$beta
            ws = res$ws
            
            if(is.na(logml_new) || is.na(logml)){
                print(paste0('Inference failed after ',i,' iterations.'))
                return(gen_res())
            }
            if( abs(logml_new-logml) < abs(1e-6*logml) )
            {
                logml = logml_new
                break
            }
            logml = logml_new
        }
        if(i==nround){
            print(paste0('Run all ',i,' iterations.',' logml=',logml))
        }else{
            print(paste0('Converge in ',i,' iterations.',' logml=',logml))
        }
        if(debug){
            idx = seq(3,nround)
            plot(idx,logml_res[idx])
        }
        return(gen_res(alpha,beta,ws,logml,Z))
    }
    em_optim = function(X, init_mus, init_sgs, init_ws, nround=200,debug=FALSE){
        res = em_optim0(X, init_mus, init_sgs, init_ws, nround, debug=debug)
        
        #debug=T
        
        if(any(is.na(res$alpha))){
            if(debug){
                print('Get NA in the opmization.')
            }
            return(res)
        }
        alpha = res$alpha
        beta = res$beta
        ws = res$ws
        logml = res$logml
        Z = res$Z
        flag = FALSE
        # try making alpha greater than 0.1
        min_alpha = 3*sigma0
        if(any(alpha<min_alpha)){
            if(debug){
                cat('Some alpha is less than min_alpha',alpha,"\n",sep='')
                disp_paras(alpha,beta,ws,logml,'curr paras: ')
            }
            K = length(alpha)
            tmus = alpha2mu(alpha)
            idx = order(tmus)
            tmus = tmus[idx]
            tvars = beta2var(beta)[idx]
            
            tws = ws
            tws[1:(K+1)] = ws[idx]
            alpha = mu2alpha(tmus)
            beta = var2beta(tvars)
            ws = tws
            
            alpha[alpha<min_alpha] = min_alpha
            flag = TRUE
        }
        # update Z and logml if parameter updated
        if(flag){
            Z = estep(X, alpha, beta, ws)
            logml = elbo(X, Z, alpha, beta, ws)
            if(debug){
                disp_paras(alpha,beta,ws,logml,'updated paras: ')
            }
        }
        
        flag = FALSE
        min_beta = 0.8
        # try making beta greater than 1
        if(any(beta<min_beta)){
            beta[beta<min_beta]=min_beta
            flag = TRUE
            
            if(debug){
                print('Some beta is less than 1.')
                disp_paras(alpha,beta,ws,logml,'curr paras: ')
            }
        }
        # update Z and logml if parameter updated
        if(flag){
            Z = estep(X, alpha, beta, ws)
            logml = elbo(X, Z, alpha, beta, ws)
            if(debug){
                disp_paras(alpha,beta,ws,logml,'updated paras: ')
            }
        }
        
        return(gen_res(alpha,beta,ws,logml,Z))
    }
    disp_paras = function(alpha,beta,ws,logml, infostr=''){
        cat(paste( infostr,
        paste0('alpha=',sprintf("%.3f", alpha),collapse=" "),
        paste0('beta=',sprintf("%.3f", beta),collapse=" "),
        paste0('ws=',sprintf("%.3f", ws),collapse=" "),
        paste0('mus=',sprintf("%.3f", alpha2mu(alpha)),collapse=" "),
        paste0('vars=',sprintf("%.3f", beta2var(beta)),collapse=" "),
        paste0('logml=',logml),
        '','',
        sep="\n" ))
    }
    em_algo = function(X,K,debug=FALSE){
        # depends on global mu0, sigma0
        xqu = unname(quantile(X, probs = seq(0.05, 1, 0.05)))
        squ = seq(log(sigma0^2),log(var(X)),length.out=21)
        
        n_trial = 20
        res_list = list(n_trial)
        logml_arr = rep(NA,n_trial)
        for(i in seq(n_trial)){
            if(debug){
                print(paste0(i,'th trial. ntrial=',n_trial))
            }
            tmpw = runif(K+2)+1
            init_ws = tmpw/sum(tmpw)
            init_mus = c(mu0, sort( sample(xqu,K)  ))
            init_sgs = c(sigma0^2, exp( sample(squ[-1],K,replace = TRUE))  )
            #res_list[[i]] = em_optim(X, init_mus, init_sgs, init_ws, debug=debug)
            res_list[[i]] = em_optim(X, init_mus, init_sgs, init_ws, debug=FALSE)
            logml_arr[i] = res_list[[i]][4]
        }
        max_ind = which.max(logml_arr)
        # list(alpha,beta,ws,logml,Z)
        res = res_list[[max_ind]]
        if(!is.na(res$logml)){
            res$bic = bic(X,res$Z,res$alpha,res$beta,res$ws)
        }else{
            res$bic = NA
        }
        return(res)
    }
    
    # main code for odgmm
    # initialization
    # mu0,sigma0
    mu0 = 0
    sigma0 = 0.1
    K = n_max_comp - 1
    
    uni_valid_inds = X!=0
    uni_log_density = -log(max(X))
    
    kflag=n_max_comp>=1
    if(!kflag){
        print('There must be at least 1 components, the dropout component. n_max_comp>=1.')
        paste0('curr n_max_comp=',n_max_comp)
        stopifnot(kflag)
    }
    
    # auxiliary function
    norm_center = function(w) {x=w[2:(length(w)-1)]; return(x/sum(x))}
    
    # list(alpha,beta,ws,logml,Z)
    #curr_res = list(NA,NA,NA,NA,NA,NA)
    curr_res = gen_res()
    curr_bic = NA
    for(i in seq(K,1,-1)){
        cat("\n")
        print(paste0('Model estimation with ',i+1,' components'))
        flag = TRUE
        res = em_algo(X,i,debug=debug)
        res_bic = res$bic
        if(is.na(curr_bic)){
            flag=FALSE
        }else if(length(res$ws)>1 && any(norm_center(res$ws)<0.1)){
            # if some component has less than 10% weight over the "real" (not dropout, uniform) component
            flag=FALSE
        }else if(res_bic-curr_bic < 0){
            flag=FALSE
        }
        if(!flag){
            curr_bic = res_bic
            curr_res = res
        }else{
            break
        }
        cat("curr_bic = ",curr_bic,"\n")
    }
    
    alpha = curr_res$alpha
    beta = curr_res$beta
    ws = curr_res$ws
    logml = curr_res$logml
    infostr = paste0('final model: ',length(alpha)+1,' components, bic=', curr_bic, ' ')
    disp_paras(alpha,beta,ws,logml,infostr)
    
    return(curr_res)
}

gmmpdf = function(x, mus, sigmas, ws, log=FALSE){
    K = length(mus)
    N = length(x)
    y = rep(0,N)
    for(k in seq(K)){
        y = y + ws[k] * dnorm(x,mus[k],sigmas[k])
    }
    if(log){
        return(log(y))
    }else{
        return (y)
    }
}

####  main code #####
#set.seed(1)

## part 1: simulate data from Gaussian mixture
N = 2000
wghs = c(0.2,0.1,0.7)
components = sample(1:3,prob=wghs,size=N,replace=TRUE)
mus = c(0, 2,  5)
sds = c(0.1, 0.5, 2)
sgs = sds^2

#X = rnorm(n=N,mean=mus[components],sd=sds[components])
#X = X_tab_with_3_peaks_Tuba1a

for(i in seq(1,dim(gene_exprs_387)[1])){
#for(i in seq(104,104)){
    cat(paste('','','',paste0('i=',i),'',sep='\n'))
    #for(i in seq(46,46)){
    pdf_file = paste0("./img/res_",i,"_plot.pdf")
    pdf(file=pdf_file)
    X = unname(gene_exprs_387[i,])
    
    ## part 2: estimating the components
    res = odgmm(X,4,debug=T)
    
    alpha = res$alpha
    beta = res$beta
    ws = res$ws
    logml = res$logml
    Z = res$Z
    final_bic = res$bic
    
    ## part 3: display results
    # display fitting with original density
    emus = res$mus
    esgs = res$sigmas
    ews = ws
    
    x = seq(-2,max(X)+1,0.01)
    ey = gmmpdf(x, emus, esgs, ews) + tail(ews,n=1)/max(X)
    #ey[x>0] = ey[x>0] + tail(ews,n=1)/max(X)
    
    x1 = emus[-1]
    ey1 = gmmpdf(x1, emus, esgs, ews) + tail(ews,n=1)/max(X)
    
    K = length(alpha)
    tryCatch(    plot(density(X,bw='sj'),ylim=c(0,0.9),main=paste0("GMM vs kernel density, i=",i,', K=', K, ', ws=', paste(round(ews[2:(K+1)], digits = 2),collapse = ' '))) ,
     error = function(cond){ plot(density(X),ylim=c(0,0.9),main=paste0("GMM vs kernel density, i=",i,', K=', K, ', ws=', paste(round(ews[2:(K+1)], digits = 2),collapse = ' ')))})
    #plot(density(X,bw='sj'),ylim=c(0,0.9),main=paste0("GMM vs kernel density, i=",i,', K=', K, ', ws=', paste(round(ews[2:(K+1)], digits = 2),collapse = ' ')))
    lines(x,ey,col="red",lwd=2)
    lines(x1,ey1,col="red",type="h")
    legend("topright",c("Estimated GMM Density","Kernel Density"),col=c("red","black"),lwd=2)
    dev.off()
}



