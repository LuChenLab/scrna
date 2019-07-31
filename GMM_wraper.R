require("doMC")
require("tictoc")
source("odgmm.R")
source("gmm_impute.R")


# get qulified expression matrix, avoid repeat calculation of this matrix when selected different rate
# :param sce SingleCellExperiment object
# :param comp number of components
# :param seed
# :param n.cores cpu numbers
# :param verbose for debug information
# :return a list contains X, X_tab and xvals
get_matrix <- function(sce, comp=2, seed=42, n.cores=10, verbose=FALSE) {
    
    mu0 = 0
    sigma0 = 0.1
    set.seed(seed)
    comp = comp
    X_tab = as.matrix(t(assays(sce)[["counts"]]))
    
    # @since 2019.07.31 reduce memory usage by finish genes selection inside each process
    registerDoMC(n.cores)
    tic("Filter")
    res <- foreach(i=1:dim(X_tab)[1], pb=icount()) %dopar% {
        X = as.numeric(X_tab[i,])   #X_tab is the matrix need to be fitted by GMM
        
        res = tryCatch(
            odgmm(X,comp), 
            error=function(e) {
                if(verbose) print(e)
                return(NA)
            }
        )
        
        progBar(pb, dim(X_tab)[1])
        
        if (length(res) == 0 || sum(is.na(unlist(res))) != 0) {
            return(NULL)
        } else {
            return(proc_res(res)) 
        }
    }
    print("")
    toc()
    
    X = as.data.frame(matrix(NA, nrow = dim(X_tab)[2], ncol = dim(X_tab)[1]))
    
    colnames(X) = rownames(X_tab)
    rownames(X) = colnames(X_tab)
    
    xvals = list()
    
    for (i in 1:length(res)) {
        temp_res = res[[i]]
        if(!is.null(temp_res)) {
            X[, i] = temp_res[[1]]
            xvals[[i]] = temp_res[[2]]
        } 
    }
    
    X = na.omit(X)
    
    return(list(
        X=X,
        X_tab=X_tab,
        xvals=xvals
    ))
}


# mat = get_matrix(sce, n.cores=20)

# Wrapper to perform impute
# :param data output of get_matrix()
# :param rate: see gmm_impute
# :param n.cores: number of cpus
# :return imputed matrix
perform_GMM_impute <- function(data, rate=0.8, n.cores=10) {
    
    ## Second step:
    ## GMM output into imputing
    X = data[["X"]]
    X_tab = data[["X_tab"]]
    
    n_gene = dim(X_tab)[1]
    n_cell = dim(X_tab)[2]
    
    xvals = data[["xvals"]]
    
    keep_rate = apply(X>0,2,sum)/n_cell
    gene_inds = which(keep_rate > rate)
    
    print(summary(keep_rate))
    
    # perform imputation
    XX = X[,gene_inds]

    XX_vals = xvals[gene_inds]

    res_after_impute = gmm_impute(XX,XX_vals, n.cores=n.cores)
    
    X_impute = res_after_impute$impute_label_mat

    X_val_impute = res_after_impute$impute_val_mat # NA means values that does not need to be imputed
    
    orig_data = t(X_tab[gene_inds,])
    
    orig_data_verse = t(X_tab[-gene_inds,])
    
    valid_inds = is.na(X_val_impute)
    
    X_val_impute[valid_inds] = orig_data[valid_inds]
    
    GMM_imputed_X_tab = cbind(orig_data_verse,X_val_impute)
    
    return(t(GMM_imputed_X_tab))
}


