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
    Z = matrix(0,nrow=N,ncol=K+1)
    for(i in seq(K+1))
    {
      Z[,i] = log(ws[i])+dnorm(X, mean = mus[i], sd = sds[i], log = TRUE)
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
    tsum = sum(ZZ)
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
    return(list(new_alpha,new_beta,new_ws))
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
    res = -2*elbo(X,Z,alpha,beta,ws) + (3*K+1)*log(N)
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
      
      if(any(is.na(res[[1]])) || any(is.na(res[[2]])) || any(is.na(res[[3]]))){
        return(list(NA,NA,NA,NA,NA))
      }
      logml_new = elbo(X, Z, res[[1]], res[[2]], res[[3]])
      logml_res[i] = logml_new
      alpha = res[[1]]
      beta = res[[2]]
      ws = res[[3]]
      #if(is.na(logml_new) || is.na(logml)){
      #  return(list(NA,NA,NA,NA,NA))
      #}
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
    return(list(alpha,beta,ws,logml,Z))
  }
  em_optim = function(X, init_mus, init_sgs, init_ws, nround=200,debug=FALSE){
    res = em_optim0(X, init_mus, init_sgs, init_ws, nround, debug=debug)
    if(any(is.na(res[[1]]))){
      if(debug){
        print('Get NA in the opmization.')
      }
      return(res)
    }
    alpha = res[[1]]
    beta = res[[2]]
    ws = res[[3]]
    logml = res[[4]]
    Z = res[[5]]
    flag = FALSE
    # try making alpha greater than 0
    if(any(alpha<0)){
      if(debug){
        print('Some alpha is less than 0.')
        disp_paras(alpha,beta,ws,logml,'curr paras: ')
      }
      tmus = alpha2mu(alpha)
      idx = order(tmus)
      tmus = tmus[idx]
      tvars = beta2var(beta)[idx]
      tws = ws[idx]
      alpha = mu2alpha(tmus)
      beta = var2beta(tvars)
      ws = tws
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
    # try making beta greater than 1
    if(any(beta<1)){
      beta[beta<1]=1
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
    
    return(list(alpha,beta,ws,logml,Z))
  }
  disp_paras = function(alpha,beta,ws,logml, infostr=''){
    cat(paste( infostr, 
                 paste0('alpha=',sprintf("%.3f", alpha),collapse=" "),
                 paste0('beta=',sprintf("%.3f", beta),collapse=" "),
                 paste0('ws=',sprintf("%.3f", ws),collapse=" "),
                 paste0('mus=',sprintf("%.3f", alpha2mu(alpha)),collapse=" "),
                 paste0('vars=',sprintf("%.3f", beta2var(beta)),collapse=" "),
                 paste0('logml=',logml), 
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
      tmpw = runif(K+1)+1
      init_ws = tmpw/sum(tmpw)
      init_mus = c(mu0, sort( sample(xqu,K)  ))
      init_sgs = c(sigma0^2, exp( sample(squ[-1],K,replace = TRUE))  )
      res_list[[i]] = em_optim(X, init_mus, init_sgs, init_ws, debug=debug)
      logml_arr[i] = res_list[[i]][4]
    }
    max_ind = which.max(logml_arr)
    # list(alpha,beta,ws,logml,Z)
    res = res_list[[max_ind]]
    if(!is.na(res[[4]])){
      res[[6]] = bic(X,res[[5]],res[[1]],res[[2]],res[[3]])
    }else{
      res[[6]] = NA
    }
    return(res)
  }
  
  # main code for odgmm
  # initialization
  # mu0,sigma0
  mu0 = 0
  sigma0 = 0.1
  K = n_max_comp - 1
  
  kflag=n_max_comp>=1
  if(!kflag){
    print('There must be at least 1 components, the dropout component. n_max_comp>=1.')
    paste0('curr n_max_comp=',n_max_comp)
    stopifnot(kflag)
  }
  
  # list(alpha,beta,ws,logml,Z)
  curr_res = list(NA,NA,NA,NA,NA,NA)
  curr_bic = NA
  for(i in seq(K,1,-1)){
    print('')
    print(paste0('Model estimation with ',i+1,' components'))
    flag = TRUE
    res = em_algo(X,i,debug=debug)
    res_bic = res[[6]]
    if(is.na(curr_bic)){
      flag=FALSE
    }else if(length(res[[3]])>1 && any(res[[3]][-1]<0.01)){
      # if some component has less than 1% weight
      flag=FALSE
    }else if(res_bic-curr_bic < -6){
      flag=FALSE
    }
    if(!flag){
      curr_bic = res_bic
      curr_res = res
    }else{
      break
    }
  }
  
  alpha = curr_res[[1]]
  beta = curr_res[[2]]
  ws = curr_res[[3]]
  logml = curr_res[[4]]
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
## part 1: simulate data from Gaussian mixture
N = 2000
wghs = c(0.2,0.1,0.7)
components = sample(1:3,prob=wghs,size=N,replace=TRUE)
mus = c(0, 2,  5)
sds = c(0.1, 0.5, 2)
sgs = sds^2

X = rnorm(n=N,mean=mus[components],sd=sds[components])

## part 2: estimating the components
res = odgmm(X,4,debug=F)

alpha = res[[1]]
beta = res[[2]]
ws = res[[3]]
logml = res[[4]]
Z = res[[5]]
final_bic = res[[6]]

## part 3: display results
# display fitting with original density
emus = alpha2mu(alpha)
esgs = beta2sigma(beta)
ews = ws

x = seq(-2,max(X)+1,0.01)
ey = gmmpdf(x, emus, esgs, ews)

plot(density(X,bw='sj'),ylim=c(0,0.9),main="GMM vs kernel density")
lines(x,ey,col="red",lwd=2)
legend("topright",c("Estimated GMM Density","Kernel Density"),col=c("red","black"),lwd=2)
