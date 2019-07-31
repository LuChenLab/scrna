library(e1071)
library(tictoc)

# progressbar
# From: https://stackoverflow.com/questions/51213293/is-it-possible-to-get-a-progress-bar-with-foreach-and-a-multicore-kind-of-back
progBar <- function(ii, N, per = 10) {
  if (ii %in% seq(1, N, per)) {
    x <- round(ii * 100 / N)
    message("[ ", 
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""), 
            " ] ", x, "%", "\r",
            appendLF = FALSE)
  }
}


# impute missing values into 
gmm_impute = function(X,X_mus_list,k=5, n.cores=10){
  n_gene = dim(X)[2]
  IX = X
  IX_val = X
  
  # @since 2019.07.31
  # using multi processes to accelerate calculation
  registerDoMC(n.cores)
  tic("Imputation")
  res = foreach(i = seq(n_gene), pb=icount()) %dopar% {
    gmm_impute0(X[,-i], X[,i], X_mus_list[[i]], k)
    progBar(pb, n_gene)
  }
  print("")
  toc()
  
  for(i in 1:length(res)) {
    IX[,i] = res[[i]][,1]
    IX_val[,i] = res[[i]][,2]
  }
  
  return(list(impute_label_mat=IX,impute_val_mat=IX_val))
}


# impute the missing values in y, based on similarity given by X
# X: n x d matrix, n cells, d genes/features
# y: n x 1 vector, labels (component id) for this gene
# y_mus: mean for each component
# k: k nearest neighbors
gmm_impute0 = function(X, y, y_mus, k){
  test_inds = which(y==0)
  train_inds = which(y!=0)
  y_impute_val = rep(NA,length(y))
  
  if(length(test_inds)==0){
    return(cbind(y,y_impute_val))
  }
  
  n_test = length(test_inds)
  n_train = length(train_inds)
  d = matrix(0, nrow = n_test, ncol = n_train)
  NN = matrix(0, nrow = n_test, ncol = n_train)
  for(i in seq(n_test)){
    for(j in seq(n_train)){
      d[i,j] = hamming.distance(          # @since 2019.07.31 force convert variable type to avoid C stack close to limit
        as.numeric(X[test_inds[i],]),
        as.numeric(X[train_inds[j],])
      )
    }
    NN[i,] = order(d[i,])
  }
  
  # d = hamming.distance(X)  # hamming distance between the rows
  # NN <- apply(d[test_inds, train_inds], 1, order)
  
  pred <- apply(NN[, 1:k, drop=FALSE], 1, function(nn){
    tmp = rle(sort(y[train_inds]))
    tmp$values[which.max(tmp$lengths)]
  })
  
  y[test_inds] = pred
  y_impute_val[test_inds] = y_mus[pred]
  return(cbind(y,y_impute_val))
}

# process the output of odgmm with uniform component
# extract the cell labels (component id) of this gene
proc_res = function(res){
  # Z = res[[5]]
  # mus = res[[1]]
  Z = res$Z  # res[[5]]#
  mus = res$mus # res[[1]]#
  
  L = dim(Z)[2]
  y = apply(Z,1,which.max)
  y = y-1
  
  new_mus = c(mus[-1],mean(mus[-1]))
  
  return(list(label=y, mus=new_mus))
}

####  NOTE: specify the GMM result and gene express matrix here  #######
GMM_res = GMM_output_in_K_4   # a list, 1 x n_gene
gene_expr_mat = t(gene_exprs_387)   # n_cell x n_gene

# collect results for all genes and select genes with high expression
n_gene = length(GMM_output_in_K_4)
n_cell = dim(gene_expr_mat)[1]
stopifnot(dim(gene_expr_mat)[2]==n_gene)

X = matrix(0, nrow = n_cell, ncol = n_gene)
xmus <- vector("list", length = n_gene)

##GMM_output is a 7879 long list contain res 
for(i in seq(n_gene)){
  cat("proc gene",i,"\n")
  tmpres = proc_res(GMM_res[[i]]$res)
  X[,i] = tmpres$label
  xmus[[i]] = tmpres$mus
}

# filter genes
keep_rate_cutoff = 0.8
keep_rate = apply(X>0,2,sum)/n_cell
gene_inds = which(keep_rate>keep_rate_cutoff)

# perform imputation
XX = X[,gene_inds]
XX_mus = xmus[gene_inds]

res = gmm_impute(XX, XX_mus)
X_impute = res$impute_label_mat
X_val_impute = res$impute_val_mat # NA means values that does not need to be imputed

orig_data = gene_expr_mat[,gene_inds]
valid_inds = is.na(X_val_impute)
X_val_impute[valid_inds] = orig_data[valid_inds]
