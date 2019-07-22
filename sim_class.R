#####################   Window class start  ############################

setClass("Window", representation(start = "numeric", end = "numeric"),
         prototype(start = NA_integer_, end = NA_integer_),
         validity = function(object){
           st = object@start
           en = object@end
           if(st<=en){         # window cannot be empty
             return(TRUE)
           }else{
             return(paste0("start should be less than end, start=",st," end=",en))
           }
        })


setMethod("<","Window",function(e1,e2){    # e1 is on the left of e2, no overlap
  #if(e1@empty || e2@empty) return(FALSE)
  return(e1@end<e2@start)
})
setMethod(">","Window",function(e1,e2){    # e1 is on the right of e2, no overlap
  #if(e1@empty || e2@empty) return(FALSE)
  return(e1@start>e2@end)
})
setMethod("<=","Window",function(e1,e2){      # e1 is a subset of e2
  return(e1@start>=e2@start && e1@end<=e2@end)
})
setMethod(">=","Window",function(e1,e2){      # e1 contains e2
  return(e1@start<=e2@start && e1@end>=e2@end)
})
setMethod("==","Window",function(e1,e2){      # e1 is identical to e2
  return(e1@start==e2@start && e1@end==e2@end)
})
setMethod("!=","Window",function(e1,e2){      # e1 and e2 overlap
  #if(e1@empty || e2@empty) return(FALSE)
  return(!(e1>e2 || e1<e2))
})

# return the union of two windows, only valid if two windows overlap
setMethod("+","Window",function(e1,e2){
  stopifnot(e1!=e2)
  st1 = e1@start
  en1 = e1@end
  st2 = e2@start
  en2 = e2@end
  return(new("Window",start=min(st1,st2),end=max(en1,en2)))
})

# return the intersection of two overlapping windows
setMethod("*","Window",function(e1,e2){
  stopifnot(e1!=e2)
  st1 = e1@start
  en1 = e1@end
  st2 = e2@start
  en2 = e2@end
  return(new("Window",start=max(st1,st2),end=min(en1,en2)))
})

# only works if e1 overlaps e2 and is not a subset or superset of it
setMethod("-","Window",function(e1,e2){
  stopifnot(e1!=e2)
  stopifnot(!e1<=e2)
  stopifnot(!e1>=e2)

  st1 = e1@start
  en1 = e1@end
  st2 = e2@start
  en2 = e2@end

  if(st2<=en1 && st1<=st2){
    return(new("Window",start=st1,end=st2-1))
  }else{
    return(new("Window",start=en2+1,end=en1))
  }
})

setMethod(f = "show",
          signature = "Window",
          definition = function(object){
            cat("start =",object@start)
            cat("  end =",object@end, "\n")
          })

# function to generate a window
window_factory = function(start, end){
  return(new("Window",start=start, end=end))
}

# function to generate a window list
win_list_factory = function(start_pos_arr, end_pos_arr){
  idx = order(start_pos_arr,end_pos_arr)
  start_pos_arr = start_pos_arr[idx]
  end_pos_arr = end_pos_arr[idx]
  nw = length(idx)
  wins = vector("list",nw)
  for(i in seq(nw)){
    wins[[i]] = window_factory(start_pos_arr[i], end_pos_arr[i])
  }
  return(wins)
}

#####################   Window class end  ############################


#####################   WindowSet class start  ############################
check_window_set = function(object){
  wins = object@windows
  allow_overlap_flag = object@allow_overlap
  
  # check each element is a "Window" object
  if(length(wins)<=0){
    return(paste0("windows list is empty"))
  }
  for(i in seq(length(wins))){
    if(class(wins[[i]])[1]!="Window"){
      return(paste0("The ",i,"th object of windows is ",class(wins[[i]])[1], '; should be Window.'))
    }
  }
  
  # no need for further check
  if(length(wins)==1){
    return(TRUE)
  }
  
  # start position in ascending order
  st_pos = get_start_pos(wins)
  en_pos = get_end_pos(wins)
  win_idx = order(st_pos,en_pos)
  if(!all(win_idx==seq(length(win_idx)))){
    #if(is.unsorted(st_pos)){
    return("windows are not sorted first acccording to start position then according to end position")
  }
  
  # consecutive windows cannot be identical
  for(i in seq(length(wins)-1)){
    if(wins[[i]]==wins[[i+1]]){
      return(paste0("window ",i," and ",i+1, "are identical."))
    }
  }
  
  # consecutive windows cannot overlap
  if(!allow_overlap_flag){
    if(is_overlap(wins)){
      return("Overlaps detected in input windows")
    }
  }
  return(TRUE)
}

setClass("WindowSet", representation(windows ="list", allow_overlap = "logical"),
         prototype(windows=list(), allow_overlap=FALSE),
         validity = check_window_set)



setMethod(f = "show",
          signature = "WindowSet",
          definition = function(object){
            cat("allow_overlap = ",object@allow_overlap, "\n")
            cat("start: ",get_start_pos(object@windows), "\n")
            cat("  end: ",get_end_pos(object@windows), "\n")
          })

setMethod("<=","WindowSet",function(e1,e2){      # e1 is a subset of e2
  stopifnot(!e1@allow_overlap)
  stopifnot(!e2@allow_overlap)
  wins1 = e1@windows
  wins2 = e2@windows
  ne1 = length(wins1)
  ne2 = length(wins2)
  
  i1 = 1
  i2 = 1
  while(i1<=ne1 && i2<=ne2){
    if(wins1[[i1]] < wins2[[i2]]){
      return(FALSE)
    }else if(wins1[[i1]] > wins2[[i2]]){
      i2 = i2 + 1
      next
    }else if(wins1[[i1]]<=wins2[[i2]]){
      i1 = i1 + 1
      next
    }else{  # case of superset (>=) or overlap (!=)
      return(FALSE)
    }
  }
  if(i1==ne1+1){
    return(TRUE)
  }else{
    return(FALSE)
  }
})
setMethod(">=","WindowSet",function(e1,e2){      # e1 contains e2
  return(e2<=e1)
})
setMethod("==","WindowSet",function(e1,e2){      # e1 is identical to e2
  return(e1<=e2 && e2<=e1)
})
setMethod("!=","WindowSet",function(e1,e2){      # e1 and e2 overlap
  stopifnot(!e1@allow_overlap)
  stopifnot(!e2@allow_overlap)
  wins1 = e1@windows
  wins2 = e2@windows
  ne1 = length(wins1)
  ne2 = length(wins2)
  
  i1 = 1
  i2 = 1
  
  while(i1<=ne1 && i2<=ne2){
    if(wins1[[i1]] < wins2[[i2]]){
      i1 = i1 + 1
      next
    }else if(wins1[[i1]] > wins2[[i2]]){
      i2 = i2 + 1
      next
    }else{
      return(TRUE)
    }
  }
  return(FALSE)
})
# return the union of two window sets（sorted by start pos, remove duplicated windows）, return a window set that allows overlap
setMethod("+","WindowSet",function(e1,e2){
  wins1 = e1@windows
  wins2 = e2@windows
  joint_wins = c(wins1, wins2)
  joint_wins = rm_duplicate(sort_win_list(joint_wins))
  
  return(new("WindowSet",windows=joint_wins, allow_overlap=TRUE))
})

# generate a WindowSet
winset_factory = function(wins_list, allow_overlap=FALSE){
  return(new("WindowSet",windows=wins_list, allow_overlap=allow_overlap))
}

# get start positions
get_start_pos = function(wins){
  if(is(wins,"list")){
    return( sapply(wins, function(n){n@start}) )
  }else if(is(wins, "WindowSet")){
    return( sapply(wins@windows, function(n){n@start}) )
  }else{
    stop("unknown input type.")
  }
}

# get end positions
get_end_pos = function(wins){
  if(is(wins,"list")){
    return( sapply(wins, function(n){n@end}) )
  }else if(is(wins, "WindowSet")){
    return( sapply(wins@windows, function(n){n@end}) )
  }else{
    stop("unknown input type.")
  }
}

# get length of each window
get_win_len = function(wins){
  if(is(wins,"list")){
    return( sapply(wins, function(n){n@end-n@start+1}) )
  }else if(is(wins, "WindowSet")){
    return( sapply(wins@windows, function(n){n@end-n@start+1}) )
  }else{
    stop("unknown input type.")
  }
}

# sort windows list
sort_win_list = function(wins){
  st_pos = get_start_pos(wins)
  en_pos = get_end_pos(wins)
  res_idx = order(st_pos,en_pos)
  return(wins[res_idx])
}

# check if a given set of exon windows (sorted acc. start pos) overlaps
# return: 1st flag, 2nd overlapping window indexes
is_overlap = function(wins){
  stopifnot(is(wins,"list"))
  
  if(length(wins)==1){
    return(FALSE)
  }
  
  start_pos = get_start_pos(wins)
  end_pos = get_end_pos(wins)
  idx = order(start_pos,end_pos)
  start_pos = start_pos[idx]
  end_pos = end_pos[idx]
  
  nw = length(start_pos)
  tmp1 = start_pos[-1]
  tmp2 = end_pos[1:(nw-1)]
  flag = any(tmp1<=tmp2)
  
  return(flag)
}

# get index of overlapping windows
get_overlap_inds = function(wins){
  stopifnot(is(wins,"list"))
  stopifnot(is_overlap(wins))
  
  start_pos = get_start_pos(wins)
  end_pos = get_end_pos(wins)
  
  nw = length(start_pos)
  tmp1 = start_pos[-1]
  tmp2 = end_pos[1:(nw-1)]
  tmpinds = which(tmp1<=tmp2)
  retinds = sort(unique(c(tmpinds, tmpinds+1)))
  
  return(retinds)
}

# is_overlap = function(wins){
#   stopifnot(is(wins,"list"))
#   start_pos = get_start_pos(wins)
#   end_pos = get_end_pos(wins)
#   
#   nw = length(start_pos)
#   res = list(flag=FALSE,overlap_inds=NA)
#   
#   if(nw>1){
#     tmp1 = start_pos[-1]
#     tmp2 = end_pos[1:(nw-1)]
#     flag = any(tmp1<=tmp2)
#     if(flag){
#       tmpinds = which(tmp1<=tmp2)
#       retinds = sort(unique(c(tmpinds, tmpinds+1)))
#       res = list(flag=TRUE,overlap_inds=retinds)
#     }
#   }
#   return(res)
# }

# collapse all overlaping windows into disjoint windowset
collapse = function(winset){
  if(!winset@allow_overlap){
    return(winset)
  }
  
  wins = winset@windows
  flag = is_overlap(wins)
  #overlap_inds = res$overlap_inds
  
  if(!flag){
    winset@allow_overlap=FALSE
    return(winset)
  }
  
  newwins = vector("list", length(wins))
  curr_ind = 1
  curr_win = NA
  
  for(i in seq(length(wins))){
    if(!is(curr_win,"Window") && is.na(curr_win)){
      curr_win = wins[[i]]
      next
    }else if(!(curr_win!=wins[[i]])){
      # output previous collapsed window first
      newwins[[curr_ind]] = curr_win
      curr_ind = curr_ind + 1
      curr_win = wins[[i]]
      next
    }else if(curr_win!=wins[[i]]){
      curr_win = curr_win + wins[[i]]
      next
    }
  }
  if(is(curr_win,"Window")){
    newwins[[curr_ind]] = curr_win
  }
  
  newwins = newwins[lapply(newwins, length) > 0]
  
  return(new("WindowSet", windows=newwins))
}

# remove duplicate windows from windows list
rm_duplicate = function(wins){
  nw = length(wins)
  if(nw==1){
    return(wins)
  }
  flag_arr = rep(T,nw)
  for(i in seq(2,nw)){
    if(wins[[i]]==wins[[i-1]]){
      flag_arr[i] = F
    }
  }
  return(wins[flag_arr])
}

# check if the windows set contains the win
contains = function(winset, win){
  stopifnot(is(winset,"WindowSet"))
  stopifnot(is(win,"Window"))
  
  wins = winset@windows
  nw = length(wins)
  
  for(i in seq(nw)){
    if(win<wins[[i]]){
      return(FALSE)
    }else if(win<=wins[[i]]){
      return(TRUE)
    }
  }
  return(FALSE)
}

# check if a WindowSet is consecutive
is_consecutive = function(winset){
  stopifnot(is(winset,"WindowSet"))
  stopifnot(!winset@allow_overlap)
  wins = winset@windows
  nw = length(wins)
  
  if(nw==1) return(TRUE)
  
  for(i in seq(nw-1)){
    if((wins[[i]]@end+1) != wins[[i+1]]@start){
      return(FALSE)
    }
  }
  return(TRUE)
}

# convert reference index (a single value) to relative index (window id + window offset)
ref_to_winset_ind = function(winset, ref_ind){
  stopifnot(is(winset,"WindowSet"))
  stopifnot(!winset@allow_overlap)
  wins = winset@windows
  st_pos = get_start_pos(wins)
  en_pos = get_end_pos(wins)
  win_ind = which(st_pos<=ref_ind & ref_ind<=en_pos)
  
  if(length(win_ind)==0){
    warning("display winset")
    print(winset)
    warning(paste0("ref_ind=",ref_ind))
    stop("ref_ind does not lie in winset.")
  }
  
  win_offset = ref_ind - st_pos[win_ind] + 1
  return(list(id=win_ind, offset=win_offset))
}

# convert winset ind (window id + offset) to reference ind (a single value)
winset_to_ref_ind = function(winset, winset_ind){
  stopifnot(is(winset,"WindowSet"))
  stopifnot(!winset@allow_overlap)
  stopifnot(is(winset_ind, "list"))
  win_ind = winset_ind$id
  win_offset = winset_ind$offset
  win_start =  winset@windows[[win_ind]]@start
  win_end =  winset@windows[[win_ind]]@end
  ref_ind = win_start + win_offset - 1
  stopifnot(ref_ind<=win_end)
  return(ref_ind)
}

# get the overlaped window sets between winset and win, winset is for a isoform, win is for a read
intersect = function(winset, win){
  stopifnot(is(winset,"WindowSet"))
  stopifnot(!winset@allow_overlap)
  stopifnot(is(win,"Window"))
  stopifnot(is_consecutive(winset))
  st_pos = win@start
  en_pos = win@end
  st_win_pos = ref_to_winset_ind(winset, st_pos)
  en_win_pos = ref_to_winset_ind(winset, en_pos)
  
  nw = en_win_pos$id - st_win_pos$id + 1
  
  if(nw==1){
    return(winset_factory(list(win)))
  }
  
  winset_st_pos = get_start_pos(winset@windows)
  winset_en_pos = get_end_pos(winset@windows)
  
  st_arr = rep(st_pos,nw)
  en_arr = rep(en_pos,nw)
  tmpinds = seq(nw-1)
  st_arr[tmpinds+1] = winset_st_pos[st_win_pos$id+tmpinds]
  en_arr[tmpinds] = winset_en_pos[st_win_pos$id+tmpinds-1]
  
  wins = win_list_factory(st_arr, en_arr)
  return(winset_factory(wins))
}

# map a WindowSet (subWinSet) from a reference system (fromWinSet) to another (toWinSet)
map_winset = function(fromWinSet, toWinSet, subWinSet){
  stopifnot(subWinSet<=fromWinSet)
  stopifnot(length(fromWinSet@windows)==length(toWinSet@windows))
  st_pos = get_start_pos(subWinSet@windows)
  en_pos = get_end_pos(subWinSet@windows)
  
  nw = length(st_pos)
  new_st_pos = st_pos
  new_en_pos = en_pos
  
  for(i in seq(nw)){
    tmpwinind = ref_to_winset_ind(fromWinSet, st_pos[i])
    new_st_pos[i] = winset_to_ref_ind(toWinSet, tmpwinind)
    tmpwinind = ref_to_winset_ind(fromWinSet, en_pos[i])
    new_en_pos[i] = winset_to_ref_ind(toWinSet, tmpwinind)
  }
  
  wins = win_list_factory(new_st_pos, new_en_pos)
  return(winset_factory(wins))
}

# get the window set as if the windows are concatenated
# this is useful to derive the local index system for an isoform
get_local_winset = function(winset){
  stopifnot(is(winset,"WindowSet"))
  stopifnot(!winset@allow_overlap)
  
  win_len = get_win_len(winset@windows)
  
  tmp = cumsum(win_len)
  new_en_pos = tmp
  new_st_pos = c(1,tmp[1:(length(tmp)-1)]+1)
  wins = win_list_factory(new_st_pos, new_en_pos)
  return(winset_factory(wins))
}

#####################   WindowSet class end  ############################


############ test functions start ####################

window_tests = function(){
  # tests windows class
  #w0 = new('Window',start=12, end=10)  # should be wrong
  
  w0 = new('Window',start=1, end=5) 
  w1 = new('Window',start=1, end=10)
  w2 = new('Window',start=10, end=10)
  w3 = new('Window',start=5, end=15)
  w4 = new('Window',start=10, end=20)
  w5 = new('Window',start=30, end=40)
  w6 = new('Window',start=1, end=40)
  w7 = new('Window',start=5, end=15)
  w8 = new('Window',start=1, end=20)
  w9 = new('Window',start=1, end=15)
  w10 = new('Window',start=1, end=4)
  w11 = new('Window',start=6, end=15)
  
  stopifnot(w1<w5)
  stopifnot(w4>w0)
  stopifnot(w2<=w1)
  stopifnot(w1>=w2)
  stopifnot(w1<=w6)
  stopifnot(w6>=w3)
  stopifnot(w3==w7)
  stopifnot(w1!=w3)
  stopifnot(!w1!=w5)
  
  stopifnot((w0+w1+w2)==w1)
  stopifnot((w3+w0)==w9)
  stopifnot((w1+w3)==w9)
  stopifnot((w3-w0)==w11)
  stopifnot((w9*w10)==w10)
  
  # test windowsets
  lis = c(w0,w1,w2,w3)
  # tmp = new("WindowSet", windows=lis)       # should error not sorted
  # tmp = new("WindowSet", windows=sort_win_list(lis))  # should error overlap windows
  tmp = new("WindowSet", windows=sort_win_list(lis), allow_overlap=TRUE)
  print(tmp)
  
  print(winset_factory(tmp@windows,tmp@allow_overlap))
  #print(win_list_factory(c(15,30,15,25),c(14,31,23,30)))  # first window error
  wins = win_list_factory(c(15,30,15,25),c(23,31,16,30))
  print(winset_factory(wins,allow_overlap = TRUE))
  print(collapse(winset_factory(wins,allow_overlap = TRUE)))
}

winset_tests = function(){
  w0 = new('Window',start=1, end=5) 
  w1 = new('Window',start=6, end=10)
  w2 = new('Window',start=16, end=20)
  w3 = new('Window',start=26, end=30)
  w4 = new('Window',start=16, end=30)
  
  w5 = new('Window',start=2, end=12)
  w6 = new('Window',start=13, end=22)
  
  lis = c(w0,w1,w2,w3,w4,w5,w6)
  ws1 = new("WindowSet", windows=sort_win_list(lis[1:5]), allow_overlap=TRUE)  # collapse
  ws2 = collapse(ws1)
  ws2
  
  ws3 = new("WindowSet", windows=sort_win_list(lis[1:4]))  # collapse
  collapse(ws3)
  
  wins = rm_duplicate(c(w0,w0,w1,w1,w2,w3))     # remove duplicate windows
  ws4 = new("WindowSet",windows=wins,allow_overlap=TRUE)
  ws4
  stopifnot(ws3==collapse(ws4))
  
  stopifnot(ws3<=ws2)
  stopifnot(ws2>=ws3)
  stopifnot(ws2==ws2)
  
  stopifnot(ws3!=ws2)
  
  ws5 = new("WindowSet", windows=sort_win_list(lis[1:7]), allow_overlap=TRUE)  # collapse
  ws5 = collapse(ws5)
  stopifnot(ws2<=ws5)
  
  ws6 = new("WindowSet", windows=sort_win_list(lis[6:7]))
  stopifnot(ws2!=ws6)
  
  print(ws2)
  print(ws5)
  print(ws2+ws5)
  print(ws2+ws5+ws2)
  
  print(get_start_pos(ws5))
  
  flag = contains(ws2+ws5,w1)
  stopifnot(flag)
  
  flag = contains(ws2+ws5,window_factory(2,13))
  stopifnot(!flag)
  
  wins = win_list_factory(c(1,11,21,31),c(10,20,30,40))
  winset = winset_factory(wins,allow_overlap = FALSE)
  flag = is_consecutive(winset)
  stopifnot(flag)
  
  flag = is_consecutive(ws2)
  stopifnot(!flag)
  
  ref_ind = 35
  print(winset)
  winset_ind = ref_to_winset_ind(winset, ref_ind)
  print(winset_ind)
  ref_ind1 = winset_to_ref_ind(winset, winset_ind)
  stopifnot(ref_ind1==ref_ind)
  
  # print(ws2)
  # ref_to_winset_ind(ws2, 15)                      # should error since 15 is not in ws2
  #  winset_to_ref_ind(ws2, list(id=2,offset=10))   # should error since ref_ind is out of range
  
  gene_st_pos = c( 1, 21, 41, 61)
  gene_en_pos = c(10, 30, 50, 70)
  wins = win_list_factory(gene_st_pos, gene_en_pos)
  
  gene_winset = winset_factory(wins)
  iso_winset = get_local_winset(gene_winset)
  
  read_win = window_factory(5,25)
  read_winset = intersect(iso_winset, read_win)
  print(read_winset)
  
  read_gene_winset = map_winset(iso_winset, gene_winset, read_winset)
  print(read_gene_winset)
  
  gene_st_pos = c( 1, 41, 61)
  gene_en_pos = c(10, 50, 70)
  wins = win_list_factory(gene_st_pos, gene_en_pos)
  
  gene_winset1 = winset_factory(wins)
  iso_winset1 = get_local_winset(gene_winset1)
  
  stopifnot(!read_gene_winset<=gene_winset1)
  
  read_win1 = window_factory(5,15)
  read_winset1 = intersect(iso_winset1, read_win1)
  print(read_winset1)

  read_gene_winset1 = map_winset(iso_winset1, gene_winset1, read_winset1)
  print(read_gene_winset1)
  
  stopifnot(read_gene_winset1<=gene_winset)
  
  mapwin_iso = map_winset(gene_winset,iso_winset,read_gene_winset1)
  print(mapwin_iso)
  
  stopifnot(!is_consecutive(mapwin_iso))
  
}

############ test functions end ####################

############ simulation functions start ####################
# get flag vector on the gene
get_flag_vec = function(winset,gene_length){
  gene_flag_vec = rep(FALSE,gene_length)
  tmp_start_pos = get_start_pos(winset)
  tmp_end_pos = get_end_pos(winset)
  
  for(iw in seq(length(tmp_start_pos))){
    tmp_st = tmp_start_pos[iw]
    tmp_en = tmp_end_pos[iw]
    gene_flag_vec[tmp_st:tmp_en] = TRUE
  }
  return(gene_flag_vec)
}

# simulate an isoform, return number of exon, exon indices
gen_isoform = function(all_exons, exon_prob=list(a=0.2,b=0.6)){
  stopifnot(is(all_exons,"WindowSet"))
  stopifnot(exon_prob$a<=exon_prob$b)
  n_all_exon = length(all_exons@windows)
  
  n_exon = round((exon_prob$a+runif(1)*(exon_prob$b-exon_prob$a))*n_all_exon)
  exon_inds = sort(sample(seq(n_all_exon),n_exon))
  
  while(is_overlap(all_exons@windows[exon_inds])){
    overlap_inds = get_overlap_inds(all_exons@windows[exon_inds])
    tmp_rm_ind = sample(overlap_inds,1)
    exon_inds = exon_inds[-tmp_rm_ind]
  }
  
  # wins = all_exons@windows[exon_inds]
  # winset_factory(wins)
  
  return( exon_inds )
}

get_iso_pmf = function(pos_pmf_all, winset){
  stopifnot(is(winset,"WindowSet"))
  
  iso_exon_start_pos = get_start_pos(winset)
  iso_exon_end_pos = get_end_pos(winset)
  
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

############ simulation functions end ####################

# window_tests()
# winset_tests()

