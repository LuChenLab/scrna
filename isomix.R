library("tictoc")

isomix = function(X, Y, S, I, FX, FY,iso_theta_ind_arr,theta_len, debug=FALSE){
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
    
    ws_mat = matrix(rep(log(ws),each=N),nrow = N)
    Z = log(I) + ws_mat + log(FX) + log(FY) + S_log_density
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
      tmp_flag_mat = FX==i | FY==i
      tmp_numerator = sum(Z[tmp_flag_mat])
      tmp_denominator = sum(2*tmp_flag_mat*Z/theta_sum_mat)
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
    new_theta = runmean(new_theta, 10)
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
  nround = 200
  
  n_fragment = dim(X)[1]
  n_isoform = dim(X)[2]
  
  mu0 = 30
  sigma0 = 3
  S_log_density = log(S>0) + dnorm(S, mean = mu0, sd = sigma0, log = TRUE)  # log(0) is -Inf
  ws = rep(1,n_isoform)/n_isoform
  theta = rep(1,theta_len)/theta_len
  
  res = em_algo(ws, theta)
  if(debug){
    cat(paste0('ws=',sprintf("%.3f", res$ws),collapse=" "),"\n")
    idx = seq(3,nround)
    plot(idx,res$lb_arr[idx],type='o')
  }
  
  return(res)
}


res = isomix(X, Y, S, I, FX, FY,iso_theta_ind_arr,theta_len, T)
cat("estimated ws: ",res$ws,"\n")
cat("real ws: ",iso_weights,"\n")

