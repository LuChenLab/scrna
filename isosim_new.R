# this script simulate isofroms and reads
source('sim_class.R')  # import all necessary functions
set.seed(5)

############  part 1:  parameter configurations  ######################

n_all_exon = 10
n_isoform = 4

gene_start_pos = 0
gene_length = 10000

exon_avg_len = 1000
exon_sd_len = 100

fragment_avg_len = 300
fragment_len_sd = 30

read_avg_len = 100
read_len_sd = 10

read_coverage_depth = 20

############  part 2:  generate gene and isoforms  ######################

# step 2.1: generate exons on gene
exon_start_pos_all = sort(sample(seq(gene_length-exon_avg_len), n_all_exon, replace = T))
exon_len_all = round(rnorm(n_all_exon, mean=exon_avg_len, sd=exon_sd_len))
exon_end_pos_all = exon_start_pos_all + exon_len_all - 1

# check if any exon out of boundary
tmpinds = exon_end_pos_all>gene_length
exon_end_pos_all[tmpinds] = gene_length
exon_len_all[tmpinds] = exon_end_pos_all[tmpinds] - exon_start_pos_all[tmpinds] + 1

wins = win_list_factory(exon_start_pos_all, exon_end_pos_all)
all_exons = winset_factory(rm_duplicate(wins), allow_overlap=TRUE) # remove identical windows

exon_start_pos_all = get_start_pos(all_exons)
exon_end_pos_all = get_end_pos(all_exons)
exon_len_all = exon_end_pos_all - exon_start_pos_all + 1

# step 2.2: generate isoforms
iso_mat = matrix(data=F, nrow=n_isoform, ncol=n_all_exon)
iso_arr = vector("list",length=n_isoform)               # isoform on gene, all index w.r.t. gene start
local_iso_arr = vector("list",length=n_isoform)         # isoform on isoform, all index w.r.t. local start
i = 1
n_trial = 0
while(i<=n_isoform){
  tmp_exon_inds = gen_isoform(all_exons)
  tmp_iso = winset_factory(all_exons@windows[tmp_exon_inds])
  if(i>1){
    flag = FALSE
    # check if the current sampled isoform is a subset/superset of previous isoforms
    for(j in seq(i-1)){
      if(tmp_iso<=iso_arr[[j]] || tmp_iso>=iso_arr[[j]]){
        flag = TRUE
        break
      }
    }
    if(flag){
      n_trail = n_trail + 1
      if(n_trail>1000) stop("can not generate suitable isoform after 1000 trails.")
      next
    }
  }
  
  iso_mat[i, tmp_exon_inds] = T
  iso_arr[[i]] = tmp_iso
  local_iso_arr[[i]] = get_local_winset(tmp_iso)
  i = i + 1
}

# step 2.3: generate collapsed isoforms and break position probability mass distribution
exon_flag_all = iso_mat[1,]
for(i in seq(n_isoform)){
  exon_flag_all = exon_flag_all | iso_mat[i,]
}
all_iso_exons = winset_factory(all_exons@windows[exon_flag_all],allow_overlap = TRUE)    # all exons of all isoforms
all_iso_exon_on_gene = collapse(all_iso_exons)                      # all exons of all isoforms collapsed on gene
all_iso_exon_local = get_local_winset(all_iso_exon_on_gene)         # all exons of all isoforms collapsed local index system

gene_flag_vec = get_flag_vec(all_iso_exon_on_gene, gene_length)
plot(gene_flag_vec,type='l')

# generate break position probability mass distribution
exon_ind_all = which(gene_flag_vec)  # mapping from collapsed isoform to gene
n_pos = sum(gene_flag_vec)
sin_period = n_pos/3     # 3 periods in all data
pos_mass = sin(seq(n_pos)/sin_period*2*pi) + 1.5
local_pos_mass = pos_mass/sum(pos_mass)
pos_pmf_all = rep(0,gene_length)
pos_pmf_all[exon_ind_all] = local_pos_mass

pos_pmf_iso_arr = vector("list",length=n_isoform)
for(i in seq(n_isoform)){
  pos_pmf_iso_arr[[i]] = get_iso_pmf(pos_pmf_all, iso_arr[[i]])
}

exon_boundary_start = get_start_pos(all_iso_exon_on_gene)
exon_boundary_end = get_end_pos(all_iso_exon_on_gene)

plot(pos_pmf_all,type='l')
lines(exon_boundary_start,pos_pmf_all[exon_boundary_start],col='red',type='h')
lines(exon_boundary_end,pos_pmf_all[exon_boundary_end],col='green',type='h')


############  part 3:  generate fragments and pair-end reads  ######################

# step 3.1: generate isoform weights, fragment lengths, read lengths, directions

# get isoform length
iso_len = rep(0,n_isoform)
for(i in seq(n_isoform)){
  iso_len[i] = sum(exon_len_all[iso_mat[i,]])
}

# generate reads from isoforms
iso_weights = runif(n_isoform)+0.1
iso_weights = iso_weights/sum(iso_weights)
n_frag = sum(gene_flag_vec) * n_isoform * read_coverage_depth  # 20X coverage
reads_boarder = round(n_frag * iso_weights)
reads_boarder = cumsum(reads_boarder)
reads_boarder[length(reads_boarder)] = n_frag
reads_boarder = c(0,reads_boarder)
frag_label = rep(0,n_frag)
frag_direc = runif(n_frag)>=0.5   # d, 0 means start, 1 means end

# fragment_avg_len = 300
# fragment_len_sd = 30
fragment_len = gen_frag_len(n_frag, fragment_avg_len, fragment_len_sd)

# read_avg_len = 100
# read_len_sd = 10
read1_len = gen_frag_len(n_frag, read_avg_len, read_len_sd)
read2_len = gen_frag_len(n_frag, read_avg_len, read_len_sd)

# reads length should be shorter than fragment length
frag_label = ifelse(read1_len>fragment_len | read2_len>fragment_len, -1, 0)

# step 3.2: generate fragments and reads, choose the actual isoforms

# X: start Y: end S: iso length I: possible iso, all positions are relative to isoforms
X = matrix(data=0, nrow=n_frag, ncol=n_isoform)
Y = matrix(data=0, nrow=n_frag, ncol=n_isoform)
S = matrix(data=0, nrow=n_frag, ncol=n_isoform)
I = matrix(data=0, nrow=n_frag, ncol=n_isoform)

R1ST = matrix(data=0, nrow=n_frag, ncol=n_isoform)
R1EN = matrix(data=0, nrow=n_frag, ncol=n_isoform)
R2ST = matrix(data=0, nrow=n_frag, ncol=n_isoform)
R2EN = matrix(data=0, nrow=n_frag, ncol=n_isoform)

for(i in seq(n_isoform)){
  tmpinds = seq((reads_boarder[i]+1),reads_boarder[i+1])
  frag_label[tmpinds] = ifelse(frag_label[tmpinds]!=-1,i,-1)  # only assign labels to fragments with valid reads
  
  tmp_pmf = pos_pmf_iso_arr[[i]]
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

# step 3.3: assign generated reads to other possible isoforms and update matrix X,Y,S,I

for(n in seq(n_frag)){
  if(n%%1000 == 0){
    cat(paste0('n=',n,"\n"))
  }
  
  real_iso_ind = frag_label[n]
  
  tmp1_st = R1ST[n, real_iso_ind]
  tmp1_en = R1EN[n, real_iso_ind]
  tmp2_st = R2ST[n, real_iso_ind]
  tmp2_en = R2EN[n, real_iso_ind]
  
  tmp1_winset = intersect(local_iso_arr[[real_iso_ind]],window_factory(tmp1_st,tmp1_en))
  tmp2_winset = intersect(local_iso_arr[[real_iso_ind]],window_factory(tmp2_st,tmp2_en))
  tmp1_read_gene_winset = map_winset(local_iso_arr[[real_iso_ind]], iso_arr[[real_iso_ind]], tmp1_winset) # read on gene
  tmp2_read_gene_winset = map_winset(local_iso_arr[[real_iso_ind]], iso_arr[[real_iso_ind]], tmp2_winset) # read on gene
  
  if(length(tmp1_winset@windows)>1){
    junc_read1_label[n] = T
  }
  if(length(tmp2_winset@windows)>1){
    junc_read2_label[n] = T
  }
  
  
  for(k in seq(n_isoform)){
    if(k==real_iso_ind){
      next
    }
    
    if(tmp1_read_gene_winset<=iso_arr[[k]] && tmp2_read_gene_winset<=iso_arr[[k]]){
      tmp1_read_local_k_winset = map_winset(iso_arr[[k]], local_iso_arr[[k]], tmp1_read_gene_winset)
      tmp2_read_local_k_winset = map_winset(iso_arr[[k]], local_iso_arr[[k]], tmp2_read_gene_winset)
      if(is_consecutive(tmp1_read_local_k_winset) && is_consecutive(tmp2_read_local_k_winset)){
        tmpst1 = get_start_pos(tmp1_read_local_k_winset)
        tmpen1 = get_end_pos(tmp1_read_local_k_winset)
        tmpst2 = get_start_pos(tmp2_read_local_k_winset)
        tmpen2 = get_end_pos(tmp2_read_local_k_winset)
        R1ST[n,k] = tmpst1[1]
        R1EN[n,k] = tmpen1[length(tmpen1)]
        
        R2ST[n,k] = tmpst2[1]
        R2EN[n,k] = tmpen2[length(tmpen2)]
        
        X[n,k] = R1ST[n,k]
        Y[n,k] = R2EN[n,k]
        S[n,k] = Y[n,k] - X[n,k] + 1
        I[n,k] = 1
        stopifnot(S[n,k]>0)
      }
    }
  }
}

