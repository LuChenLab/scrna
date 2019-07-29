library("tictoc")

isomix = function(X, Y, S, I, FX, FY,iso_theta_ind_arr, aux_reads_theta_map, theta_len,
                  mu0, sigma0, runmean_win_len, debug=FALSE){
  get_theta_sum = function(theta){
    K = length(iso_theta_ind_arr)
    theta_sum = rep(0,K)
    for(k in seq(K)){
      theta_sum[k] = sum(theta[iso_theta_ind_arr[[k]]])
    }
    return(theta_sum)
  }
  
  estep = function(ws, theta){
    N = n_fragment
    K = n_isoform
    # tic()
    
    Z = matrix(-Inf, nrow = N, ncol = K)
    tmpx = matrix(0,nrow=N, ncol=K)
    tmpy = matrix(0,nrow=N, ncol=K)
    ws_mat = matrix(rep(log(ws),each=N),nrow = N)
    
    tmpinds = FX>0
    tmpx[tmpinds] = theta[FX[tmpinds]]
    tmpinds = FY>0
    tmpy[tmpinds] = theta[FY[tmpinds]]
    Z = log(I) + ws_mat + log(tmpx) + log(tmpy) + S_log_density
    maxz = apply(Z,1,max)
    Z = exp(Z - maxz)
    
    # toc()
    if(debug){
      cat("e-step","\n")
    }
    
    return (Z/apply(Z,1,sum))
  }
  
  cal_ws = function(Z){
    N = n_fragment
    ws = apply(Z*I,2,sum)/N
    return(ws)
  }
  
  #iso_theta_ind_arr
  cal_theta = function(Z, theta_sum){
    N = n_fragment
    K = n_isoform
    stopifnot(length(theta_sum)==K)
    
    tic()
    
    Z = Z * I
    theta = rep(0,theta_len)
    theta_sum_mat = matrix(rep(theta_sum,each=N),nrow = N)
    
    for(i in seq(theta_len)){
      # cat("theta_",i,"\n")
      if(length(aux_reads_theta_map$numerator[[i]])==0) next
      
      #tmp_flag_mat = FX==i | FY==i
      tmpinds = aux_reads_theta_map$numerator[[i]]
      # tmp_numerator = sum(Z[tmp_flag_mat])
      tmp_numerator = sum(Z[tmpinds])
      # tmp_denominator = sum(2*tmp_flag_mat*Z/theta_sum_mat)
      tmpinds = aux_reads_theta_map$denominator[[i]]
      tmp_denominator = 2*Z[,tmpinds]/theta_sum_mat[,tmpinds]
      tmp_denominator = sum(tmp_denominator)
      theta[i] = tmp_numerator/tmp_denominator
    }
    
    toc()
    if(debug){
      cat("Finished calculating theta\n")
    }
    
    return(theta)
  }
  
  mstep = function(Z, theta){
    new_ws = cal_ws(Z)
    theta_sum = get_theta_sum(theta)
    new_theta = cal_theta(Z,theta_sum)
    # new_theta = runmean(new_theta, 50)
    new_theta = runmean(new_theta, runmean_win_len)
    new_theta = new_theta/sum(new_theta)
    return(list(ws=new_ws,theta=new_theta))
  }
  
  exp_log_like = function(Z,ws,theta){
    N = dim(Z)[1]
    K = dim(Z)[2]
    loglik = 0
    theta_sum = get_theta_sum(theta)
    
    for(n in seq(N)){
      for(k in seq(K)){
        if(I[n,k]>0){
          tmp = log(ws[k]) + log(theta[FX[n,k]]) + log(theta[FY[n,k]]) + S_log_density[n,k] - 2*log(theta_sum[k])
          loglik = loglik + tmp * Z[n,k]
        }
      }
    }
    return(loglik)
  }
  
  entropy = function(Z){
    LZ = I
    LZ[I>0] = Z[I>0] * log(Z[I>0])
    return(sum(LZ))
  }
  
  elbo = function(Z, ws, theta){
    lb = exp_log_like(Z,ws,theta) - entropy(Z)
    return(lb)
  }
  
  em_algo = function(ws, theta){
    Z = estep(ws, theta)
    lb = -Inf
    lb_arr = rep(NA,nround)
    for(i in seq(nround)){
      if(debug){
        cat('iteration=',i,'  lb=',lb, "\n")
      }
      Z = estep(ws, theta)
      res = mstep(Z, theta)
      theta = res$theta
      ws = res$ws
      lb_new = elbo(Z, ws, theta)
      lb_arr[i] = lb_new
      
      if( abs(lb_new-lb) < abs(1e-6*lb) )
      {
        break
      }else{
        lb = lb_new
      }
    }
    if(i==nround){
      cat('Run all ',i,' iterations.',' lb=',lb, "\n")
    }else{
      cat('Converge in ',i,' iterations.',' lb=',lb, "\n")
    }
    return(list(ws=ws, theta=theta, lb_arr=lb_arr))
  }
  
  # main code of isomix
  nround = 400
  
  n_fragment = dim(X)[1]
  n_isoform = dim(X)[2]
  
  S_log_density = log(S>0) + dnorm(S, mean = mu0, sd = sigma0, log = TRUE)  # log(0) is -Inf
  ws = runif(n_isoform)+1
  ws = ws/sum(ws)
  theta = rep(1,theta_len)/theta_len
  
  res = em_algo(ws, theta)
  if(debug){
    cat(paste0('ws=',sprintf("%.3f", res$ws),collapse=" "),"\n")
    idx = seq(3,nround)
    if(!is.na(res$lb_arr[idx[1]])){
      plot(idx,res$lb_arr[idx],type='o')
    }
  }
  
  return(res)
}

res = isomix(X, Y, S, I, FX, FY,iso_theta_ind_arr,aux_reads_theta_map,theta_len,
             mu0=fragment_avg_len, sigma0=fragment_len_sd, runmean_win_len = 30, T)
cat("estimated ws: ",res$ws,"\n")
cat("real ws: ",iso_weights,"\n")

est_pmf = rep(0, theta_len)
iso_flag_arr = vector("list",length=n_isoform)
for(i in seq(n_isoform)){
  tmpisowin = map_winset(all_iso_exon_on_gene, all_iso_exon_local, iso_arr[[i]])
  tmpflag = get_flag_vec(tmpisowin, theta_len)
  iso_flag_arr[[i]] = tmpflag
  est_pmf[tmpflag] = est_pmf[tmpflag] + iso_weights[i] * (res$theta[tmpflag]/sum(res$theta[tmpflag]))
}

plot(real_pmf[exon_ind_all], type="l",ylim=c(0,2*max(real_pmf[exon_ind_all])), main = "pmf of collapsed exons on gene")
lines(density(c(frag_st_gene,frag_en_gene),bw='sj'),col="red")
lines(est_pmf,col="green")
abline(v = get_start_pos(all_iso_exon_local), col="blue", lty=2)
legend("topleft", legend=c("true", "data", "estimated"),
       col=c("black","red", "green"), lty=1, cex=0.8)

# par(mfrow=c(2,1))
# par(mar=c(1,1,1,1))
# barplot(as.matrix(diff(get_end_pos(all_iso_exon_local))), horiz=T, col=greens(10)[10*c(1,1,1,1)], axes=F)
# # https://stackoverflow.com/questions/9100841/rectangle-bar-graph-filled-with-color-and-distance-using-r-base-r-or-ggplot2-or

# plot(local_pos_mass,type="l",col="red",ylim=c(0,max(local_pos_mass)*2))
# lines(res$theta,col="green")
# abline(v = get_start_pos(all_iso_exon_local), col="blue", lty=2)

# par(mfrow=c(n_isoform,1))
for(i in seq(n_isoform)){
  tmppmf = local_pos_mass[iso_flag_arr[[i]]]
  tmppmf = tmppmf/sum(tmppmf)
  plot(tmppmf,type="l",col="red",ylim=c(0,max(tmppmf)*2),main=paste0("isoform",i,"pmf"))
  tmppmf = res$theta[iso_flag_arr[[i]]]
  tmppmf = tmppmf/sum(tmppmf)
  lines(tmppmf,col="green")
  abline(v = get_start_pos(local_iso_arr[[i]]), col="blue", lty=2)
}

plot(local_pos_mass,type="l",main="theta pmf",ylim=c(0,2*max(local_pos_mass)), xlab = "index on full pos set")
lines(res$theta,col="red")
legend("topleft", legend=c("true", "estimated"),
       col=c("black","red"), lty=c(1,1), cex=0.8)


