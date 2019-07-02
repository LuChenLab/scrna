library(e1071)

# impute missing values into 
gmm_impute = function(X,X_val,k=5){
  n_gene = dim(X)[2]
  IX = X
  IX_val = X_val
  for(i in seq(n_gene)){
    res = gmm_impute0(X[,-i], X[,i], X_val[,i], k)
    IX[,i] = res[,1]
    IX_val[,i] = res[,2]
  }
  return(list(IX,IX_val))
}

# impute the missing values in y, based on similarity given by X
# X: n x d matrix, n cells, d genes/features
# y: n x 1 vector, labels (component id) for this gene
# yval: mean for each component
# k: k nearest neighbors
gmm_impute0 = function(X, y, yval, k){
  test_inds = which(y==0)
  train_inds = which(y!=0)
  
  n_test = length(test_inds)
  n_train = length(train_inds)
  d = matrix(0, nrow = n_test, ncol = n_train)
  NN = matrix(0, nrow = n_test, ncol = n_train)
  for(i in seq(n_test)){
    for(j in seq(n_train)){
      d[i,j] = hamming.distance(X[i,],X[j,])
    }
    NN[i,] = order(d[i,])
  }
  
  #d = hamming.distance(X)  # hamming distance between the rows
  #NN <- apply(d[test_inds, train_inds], 1, order)
  
  pred <- apply(NN[, 1:k, drop=FALSE], 1, function(nn){
    tmp = rle(sort(y[train_inds]))
    tmp$values[which.max(tmp$lengths)]
  })
  
  y[test_inds] = pred
  return(cbind(y,yval[pred]))
  
}

alpha2mu = function(alpha){
  # global var mu0
  mu0 = 0 # warning: check the actual value in odgmm.R
  return ( cumsum(c(mu0,alpha)) )
}

# process the output of odgmm with uniform component
# extract the cell labels (component id) of this gene
# TODO: modify the output of odgmm.R
proc_res = function(res){
  Z = res[[5]]
  alpha = res[[1]]
  mus = alpha2mu(alpha)
  
  L = dim(Z)[2]
  y = apply(Z,1,which.max)
  y = y-1
  
  new_mus = c(mus[-1],mean(mus[-1]))
  
  return(list(y, new_mus))
}

# collect results for all genes and select genes with high expression
n_gene = 100
n_cell = 200

X = matrix(0, nrow = n_cell, ncol = n_gene)
xvals <- vector("list", length = n_gene)

for(i in seq(n_gene)){
  # load result of geie i
  #TO DO 
  
  tmpres = proc_res(res)
  X[,i] = tmpres[[1]]
  xvals[[i]] = tmpres[[2]]
}

# filter genes
keep_rate = apply(X>0,2,sum)/n_cell
gene_inds = which(keep_rate>0.8)

# perform imputation
XX = X[,gene_inds]
XX_val = xvals[gene_inds]

res = gmm_impute(XX, XX_val)
X_impute = res[[1]]
X_val_impute = res[[2]]
