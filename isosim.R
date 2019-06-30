# simulate isofroms and reads
set.seed(5)

n_all_exon = 10
n_isoform = 4

gene_start_pos = 0
gene_length = 10000

exon_avg_len = 1000
exon_sd_len = 100
exon_start_pos_all = sort(sample(seq(gene_length-exon_avg_len), n_all_exon, replace = T))
exon_len_all = round(rnorm(n_all_exon, mean=exon_avg_len, sd=exon_sd_len))
exon_end_pos_all = exon_start_pos_all + exon_len_all - 1

# check if any exon out of boundary
tmpinds = exon_end_pos_all>gene_length
exon_end_pos_all[tmpinds] = gene_length
exon_len_all[tmpinds] = exon_end_pos_all[tmpinds] - exon_start_pos_all[tmpinds] + 1
if(any(tmpinds)){
  stopifnot(!any(diff(sort(exon_start_pos_all[tmpinds]))==0) )  # start postion cannot be the same
}

# generate isoforms
iso_mat = matrix(data=F, nrow=n_isoform, ncol=n_all_exon)
for(i in seq(n_isoform)){
  res = gen_isoform(n_all_exon, exon_start_pos_all, exon_end_pos_all)
  tmp_n_exon = res[1]
  tmp_exon_inds = res[-1]
  iso_mat[i, tmp_exon_inds] = T
}

# check no isoform is contained in another isoform
kep_iso_set = seq(n_isoform)
rm_iso_set = vector()
for(i in seq(n_isoform)){
  for(j in seq(n_isoform)){
    if(i==j){
      next
    }
    if(any(i == rm_iso_set)){
      next
    }
    if(any(j == rm_iso_set)){
      next
    }
    st_i = exon_start_pos_all[iso_mat[i,]]
    en_i = exon_end_pos_all[iso_mat[i,]]
    st_j = exon_start_pos_all[iso_mat[j,]]
    en_j = exon_end_pos_all[iso_mat[j,]]
    if(check_include(st_i, en_i, st_j, en_j)){  # iso i contain j
      rm_iso_set = c(rm_iso_set,j)
      next
    }
    if(check_include(st_j, en_j, st_i, en_i)){
      rm_iso_set = c(rm_iso_set,i)
      next
    }
  }
}
# if some isoforms needs to be removed
if(length(rm_iso_set)>0){
  paste('isoforms to be removed: ',paste0(rm_iso_set))
  kep_iso_set = setdiff(kep_iso_set,rm_iso_set)
  n_isoform = length(kep_iso_set)
  iso_mat = iso_mat[kep_iso_set,]
}

# generate start position probability mass distribution
exon_flag_all = iso_mat[1,]
for(i in seq(n_isoform)){
  exon_flag_all = exon_flag_all | iso_mat[i,]
}
gene_flag_vec = rep(F,gene_length)
for(iex in which(exon_flag_all)){
  tmp_st = exon_start_pos_all[iex]
  tmp_en = exon_end_pos_all[iex]
  gene_flag_vec[seq(tmp_st,tmp_en)] = T
}
exon_ind_all = which(gene_flag_vec)  # mapping from collapsed isoform to gene
plot(gene_flag_vec,type='l')

# get exon boundary
exon_boundary = get_boundary(gene_flag_vec)
exon_boundary_start = exon_boundary[,1]
exon_boundary_end = exon_boundary[,2]

pos_mass = runmean(runif(length(exon_ind_all)),window=20)  # mean over [p-window,p+window]
pos_pmf_all = rep(0,gene_length)
pos_pmf_all[exon_ind_all] = pos_mass/sum(pos_mass)

plot(pos_pmf_all,type='l')
lines(exon_boundary_start,pos_pmf_all[exon_boundary_start],col='red',type='h')
lines(exon_boundary_end,pos_pmf_all[exon_boundary_end],col='green',type='h')


#plot(seq(length(pos_pmf_all)),pos_pmf_all,type='l')
#plot(exon_ind_all,pos_pmf_all,type='l',ylim = c(-max(pos_pmf_all),max(pos_pmf_all)))
#lines(exon_boundary,pos_pmf_all,col='red')
#points(which(!gene_flag_vec),rep(0,sum(!gene_flag_vec)),col='red')

# get isoform length
iso_len = rep(0,n_isoform)
for(i in seq(n_isoform)){
  iso_len[i] = sum(exon_len_all[iso_mat[i,]])
}

# generate reads from isoforms
iso_weights = runif(n_isoform)+0.1
iso_weights = iso_weights/sum(iso_weights)
n_frag = sum(gene_flag_vec) * n_isoform * 20  # 20X coverage
tmpb = round(n_frag * iso_weights)
tmpb = cumsum(tmpb)
tmpb[length(tmpb)] = n_frag
tmpb = c(0,tmpb)
frag_label = rep(0,n_frag)
frag_direc = runif(n_frag)>=0.5   # d, 0 means start, 1 means end

fragment_avg_len = 300
fragment_len_sd = 30
fragment_len = gen_frag_len(n_frag, fragment_avg_len, fragment_len_sd)

read_avg_len = 100
read_len_sd = 10
read1_len = gen_frag_len(n_frag, read_avg_len, read_len_sd)
read2_len = gen_frag_len(n_frag, read_avg_len, read_len_sd)

# reads length should be shorter than fragment length
frag_label = ifelse(read1_len>fragment_len | read2_len>fragment_len, -1, 0)

# X: start Y: end S: iso length I: possible iso, all positions are relative to concatenated exons
X = matrix(data=0, nrow=n_frag, ncol=n_isoform)
Y = matrix(data=0, nrow=n_frag, ncol=n_isoform)
S = matrix(data=0, nrow=n_frag, ncol=n_isoform)
I = matrix(data=0, nrow=n_frag, ncol=n_isoform)

R1ST = matrix(data=0, nrow=n_frag, ncol=n_isoform)
R1EN = matrix(data=0, nrow=n_frag, ncol=n_isoform)
R2ST = matrix(data=0, nrow=n_frag, ncol=n_isoform)
R2EN = matrix(data=0, nrow=n_frag, ncol=n_isoform)

for(i in seq(n_isoform)){
  tmpinds = seq((tmpb[i]+1),tmpb[i+1])
  frag_label[tmpinds] = ifelse(frag_label[tmpinds]!=-1,i,-1)  # only assign labels to fragments with valid reads
  
  tmp_pmf = get_iso_pmf(pos_pmf_all, exon_start_pos_all, exon_end_pos_all, iso_mat[i,])
  # print(paste0('i=',i))
  # print(paste0('pmf_len=',length(tmp_pmf)))
  # print(paste0('iso_len=',iso_len[i]))
  
  tmppos = sample(seq(iso_len[i]),length(tmpinds), replace = TRUE, prob = tmp_pmf)
  tmp_direc = frag_direc[tmpinds]
  
  # generate reads from ith isoform
  for(jj in seq(length(tmpinds))){
    j = tmpinds[jj]   # index of fragment over all fragments (n_frag)
    tp = tmppos[jj]   # break position in the sampled fragment
    
    if(!frag_direc[j]){
      # tp is the start position of the fragment
      if(tp+fragment_len[j]>iso_len[i]){
        frag_label[j] = -1
        next
      }else{
        # set fragment position
        X[j,i] = tp
        Y[j,i] = tp+fragment_len[j]-1
        S[j,i] = fragment_len[j]
        I[j,i] = 1
        
        # set reads position
        R1ST[j,i] = tp
        R1EN[j,i] = tp + read1_len[j] - 1
        R2ST[j,i] = tp + fragment_len[j] - read2_len[j]
        R2EN[j,i] = tp + fragment_len[j] - 1
      }
    }else{
      # tp is the end position of the fragment
      if(tp-fragment_len[j]<0){
        frag_label[j] = -1
        next
      }else{
        # set fragment position
        X[j,i] = tp - fragment_len[j] + 1
        Y[j,i] = tp
        S[j,i] = fragment_len[j]
        I[j,i] = 1
        
        # set read position
        R1ST[j,i] = tp - fragment_len[j] + 1
        R1EN[j,i] = tp - fragment_len[j] + read1_len[j]
        R2ST[j,i] = tp - read2_len[j] + 1
        R2EN[j,i] = tp
      }
    }
  }
}

# get rid of invalid fragments
valid_inds = which(frag_label!=-1)
n_frag = length(valid_inds)
frag_label = frag_label[valid_inds]
frag_direc = frag_direc[valid_inds]
fragment_len = fragment_len[valid_inds]
read1_len = read1_len[valid_inds]
read2_len = read2_len[valid_inds]
X = X[valid_inds,]
Y = Y[valid_inds,]
I = I[valid_inds,]
S = S[valid_inds,]
R1ST = R1ST[valid_inds,]
R1EN = R1EN[valid_inds,]
R2ST = R2ST[valid_inds,]
R2EN = R2EN[valid_inds,]

junc_read1_label =  rep(F,n_frag)
junc_read2_label =  rep(F,n_frag)

read1_gene_pos = matrix(0,nrow=n_frag,ncol=2)
read2_gene_pos = matrix(0,nrow=n_frag,ncol=2)


# assign generated reads to other possible isoforms and update matrix X,Y,S,I
for(n in seq(n_frag)){
  if(n%%1000 == 0){
    cat(paste0('n=',n,"\n"))
  }
  
  real_iso_ind = frag_label[n]
  
  tmp1_st = R1ST[n, real_iso_ind]
  tmp1_en = R1EN[n, real_iso_ind]
  tmp2_st = R2ST[n, real_iso_ind]
  tmp2_en = R2EN[n, real_iso_ind]
  
  r1_gene_wins = map_isowin_to_gene(exon_start_pos_all, exon_end_pos_all, exon_len_all, iso_mat[real_iso_ind,], tmp1_st, tmp1_en)
  r2_gene_wins = map_isowin_to_gene(exon_start_pos_all, exon_end_pos_all, exon_len_all, iso_mat[real_iso_ind,], tmp2_st, tmp2_en)
  
  read1_gene_pos[n,1] = r1_gene_wins[1]
  read1_gene_pos[n,2] = r1_gene_wins[length(r1_gene_wins[1])]
  read2_gene_pos[n,1] = r2_gene_wins[1]
  read2_gene_pos[n,2] = r2_gene_wins[length(r2_gene_wins[1])]
  
  if(nrow(r1_gene_wins)>1){
    junc_read1_label[n] = T
  }
  if(nrow(r2_gene_wins)>1){
    junc_read2_label[n] = T
  }
  
  
  for(k in seq(n_isoform)){
    if(k==real_iso_ind){
      next
    }
    
    tmp_k_st = exon_start_pos_all[iso_mat[k,]]
    tmp_k_en = exon_end_pos_all[iso_mat[k,]]
    
    tmp_flag1 = check_include(tmp_k_st,tmp_k_en,r1_gene_wins[,1],r1_gene_wins[,2])
    tmp_flag2 = check_include(tmp_k_st,tmp_k_en,r2_gene_wins[,1],r2_gene_wins[,2])
    if(tmp_flag1 && tmp_flag2){
      R1ST[n,k] = map_gene_ind_to_iso(exon_start_pos_all, exon_end_pos_all, exon_len_all, iso_mat[k,], r1_gene_wins[1])
      R1EN[n,k] = map_gene_ind_to_iso(exon_start_pos_all, exon_end_pos_all, exon_len_all, iso_mat[k,], r1_gene_wins[length(r1_gene_wins)])
      R2ST[n,k] = map_gene_ind_to_iso(exon_start_pos_all, exon_end_pos_all, exon_len_all, iso_mat[k,], r2_gene_wins[1])
      R2EN[n,k] = map_gene_ind_to_iso(exon_start_pos_all, exon_end_pos_all, exon_len_all, iso_mat[k,], r2_gene_wins[length(r2_gene_wins)])
      
      X[n,k] = R1ST[n,k]
      Y[n,k] = R2EN[n,k]
      S[n,k] = Y[n,k] - X[n,k] + 1
      I[n,k] = 1
      stopifnot(S[n,k]>0)
    }
  }
}


# map an index on the gene to an index on the isoform
map_gene_ind_to_iso = function(exon_sta_all, exon_end_all, exon_len_all, iso_exon_flag, ind){
  iso_boundary = c(0,cumsum(exon_len_all[iso_exon_flag]))
  
  res_pos = get_win_index(exon_sta_all[iso_exon_flag], exon_end_all[iso_exon_flag], ind)  # exon_ind, ind_in_exon
  iex = res_pos[1]
  iex_ind = res_pos[2]
  iso_ind = iso_boundary[iex] + iex_ind
  
  return(iso_ind)
}

# map a window on isoform to windows on the gene, may split into multiple windows
map_isowin_to_gene = function(exon_sta_all, exon_end_all, exon_len_all, iso_exon_flag, isowin_sta, isowin_end){
  exon_inds = which(iso_exon_flag)
  iso_boundary = c(0,cumsum(exon_len_all[iso_exon_flag]))
  n_exon = length(exon_inds)
  iso_boundary_sta = iso_boundary[1:n_exon]+1
  iso_boundary_end = iso_boundary[2:(n_exon+1)]
  
  res_sta = get_win_index(iso_boundary_sta, iso_boundary_end, isowin_sta)  # exon_ind, ind_in_exon
  res_end = get_win_index(iso_boundary_sta, iso_boundary_end, isowin_end)
  
  stopifnot(length(res_sta)==2)
  stopifnot(length(res_end)==2)
  
  n_ret_win = res_end[1]-res_sta[1]+1
  ret_win = matrix(0, nrow=n_ret_win, ncol = 2)
  for(i in seq(res_sta[1], res_end[1])){
    if(i==res_sta[1]){
      ret_win[i-res_sta[1]+1,1] = exon_sta_all[exon_inds[i]] + res_sta[2] -1
    }else{
      ret_win[i-res_sta[1]+1,1] = exon_sta_all[exon_inds[i]]
    }
    
    if(i==res_end[1]){
      ret_win[i-res_sta[1]+1,2] = exon_sta_all[exon_inds[i]] + res_end[2] -1
    }else{
      ret_win[i-res_sta[1]+1,2] = exon_end_all[exon_inds[i]]
    }
  }
  
  return(ret_win)
}

# get the window and index in the window for a given index
get_win_index = function(win_sta_pos, win_end_pos, ind){
  res = c(-1,NA)
  for(i in seq(length(win_sta_pos))){
    if(ind>=win_sta_pos[i] && ind<=win_end_pos[i]){
      res = c(i, ind-win_sta_pos[i]+1)
      break
    }
  }
  return(res)
}

# get the density of break points on an isoform
# pos_pmf_all is the density over *the whole gene*, where unrelated parts are zero density
get_iso_pmf = function(pos_pmf_all, exon_start_pos_all, end_pos_all, iso_exon_flag){
  stopifnot(length(exon_start_pos_all)==length(exon_end_pos_all))
  iso_exon_start_pos = exon_start_pos_all[iso_exon_flag]
  iso_exon_end_pos = exon_end_pos_all[iso_exon_flag]
  
  tmpflag = rep(F,length(pos_pmf_all))
  for(i in seq(length(iso_exon_start_pos))){
    tmpflag[iso_exon_start_pos[i]:iso_exon_end_pos[i]] = T
  }
  
  iso_den = pos_pmf_all[tmpflag]
  iso_den = iso_den/sum(iso_den)
  return(iso_den)
}

# generate integers that are normally distributed
gen_frag_len = function(n, mu=0, sigma=1){
  mu = round(mu)
  sigma = round(sigma)
  ymin = mu - 3*sigma
  ymax = mu + 3*sigma
  
  y = round(rnorm(n, mu, sigma))
  y[y<ymin] = ymin
  y[y>ymax] = ymax
  
  return(y)
}



# simulate an isoform, return number of exon, exon indices
gen_isoform = function(n_all_exon, exon_st_pos, exon_en_pos){
  n_exon = round((0.2+runif(1)*0.4)*n_all_exon)
  exon_inds = sort(sample(seq(n_all_exon),n_exon))
  stopifnot(!is.unsorted(exon_st_pos))  # exon_st_pos need to be sorted
  tmp_st = exon_st_pos[exon_inds]
  tmp_en = exon_en_pos[exon_inds]
  
  # remove overlap exon by iteratively ramdom taking out one exon
  res = is_overlap(exon_st_pos[exon_inds], exon_en_pos[exon_inds])
  overlap_flag = res[[1]]
  overlap_inds = res[[2]]
  while(overlap_flag){
    tmp_rm_ind = sample(overlap_inds,1)
    exon_inds = exon_inds[-tmp_rm_ind]
    res = is_overlap(exon_st_pos[exon_inds], exon_en_pos[exon_inds])
    overlap_flag = res[[1]]
    overlap_inds = res[[2]]
  }
  n_exon = exon_inds
  
  return(c(n_exon,exon_inds))
}

# check if a given set of exon windows overlaps
is_overlap = function(start_pos, end_pos){
  stopifnot(length(start_pos)==length(end_pos))
  stopifnot(!is.unsorted(start_pos))
  stopifnot(all(start_pos<end_pos))
  
  nw = length(start_pos)
  res = list(FALSE,NA)
  
  if(nw>1){
    tmp1 = start_pos[-1]
    tmp2 = end_pos[1:(nw-1)]
    flag = any(tmp1<=tmp2)
    if(flag){
      tmpinds = which(tmp1<=tmp2)
      retinds = unique(c(tmpinds, tmpinds+1))
      res = list(T,retinds)
    }
  }
  return(res)
}

# check if isoform A include isoform B
check_include = function(st_A, en_A, st_B, en_B){
  n_ex_A = length(st_A)
  n_ex_B = length(st_B)
  ia = 1
  ib = 1
  flag_arr = rep(FALSE,n_ex_B)
  while(ia<=n_ex_A && ib<=n_ex_B){
    if(st_A[ia]<=st_B[ib] && en_A[ia]>=en_B[ib]){
      flag_arr[ib] = TRUE
      ib = ib + 1
    }else if(st_B[ib]>en_A[ia]){
      ia = ia + 1
    }else{
      return(FALSE)
    }
  }
  if(all(flag_arr)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


# get boundaries of consecutive 1s in a flag vector
get_boundary = function(bool_vec){
  stopifnot(any(bool_vec))  # at least one element need to be true
  tmpvec = c(F,bool_vec,F)
  tmpdiff = diff(tmpvec)
  st_pos = which(tmpdiff==1)
  en_pos = which(tmpdiff==-1)-1
  return(cbind(st_pos,en_pos))
}

# sliding window mean, average of the input data over "-window, +window" , jump by "step"
# http://coleoguy.blogspot.com/2014/04/sliding-window-analysis.html
runmean <- function(data, window, step=1){
  total <- length(data)
  spots <- seq(from=1, to=(total), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    st = spots[i]-window
    en = spots[i]+window
    if(st<0){
      st = 1
    }
    if(en>total){
      en = total
    }
    result[i] <- mean(data[st:en])
    #result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}