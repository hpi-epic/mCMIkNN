getCounts = function(sr){
    counts = sapply(sr$local_index, length)
    return(counts)
}
getVols = function(sr){
  ll = (length(sr) - 2) / 2
  if(ll < 1){
    stop("Error while calculationg volumes: Number of dimensions is smaller than 1 or sr is corrupted.")
  }
  vols = sr[[2]] - sr[[1]]
  if(ll > 1){
    for(i in 2:ll){
      ii = i * 2
      vols = vols * (sr[[ii]] - sr[[ii-1]])
    }
  }
  return(vols)
}
getSCLog = function(sr){
    vols = getVols(sr)
    counts = getCounts(sr)
    notZero = counts != 0 ## only sum non-empty bins
    n = sum(counts)
    dims = length(counts)
    lh = sum( -counts[notZero]*log(counts[notZero] / n))
    reg = regret(n=n, k=dims)
    return(list(L=lh, R=reg))
}
getSC = function(sr, n=-1, nonEmpty=F, asymp=F, noVol=F){
    vols = getVols(sr)
    counts = getCounts(sr)
    notZero = counts != 0 ## only sum non-empty bins
    if(n == -1){
      n = sum(counts)
    }
    dims = length(counts)
    dims.ne = length(counts[notZero])
    lh = sum( counts[notZero] * (-log2(counts[notZero]) + log2(n * vols[notZero])) )
    if(noVol){
      lh = sum( -counts[notZero]*log2(counts[notZero] / n))
    }
    reg = 0
    if(nonEmpty){
      reg = regret(n=n, k=dims.ne) + log2nChoosek(max(dims, dims.ne), min(dims, dims.ne))
    }else{
      if(asymp){
        reg = (dims - 1) / 2.0 * log2(n)
      }else{
        reg = regret(n=n, k=dims)
      }
    }
    return(list(L=lh, R=reg))
}
srToVec = function(sr){
  inds = rep(-1, length(unlist(sr$local_index)))
  for(ind in 1:length(sr$local_index)){
    for(i in sr$local_index[[ind]]){
      inds[i] = ind
    }
  }
  if(sum(inds == -1) > 0){
    stop("Error transforming lists to vec.")
  }
  return(inds)
}
transformSRVec = function(sr, old_inds, breaks, inds){
  ### first transform sr
  ll = (length(sr) - 2) / 2
  cut_mat <- t(rep(NA, 2*(ll+1)))
  new_inds = list()
  indsI = 1
  if(length(breaks) < 2){
    stop("Error in getSCBreaks: Less than two break points")
  }
  count_bins = length(sr[[1]])
  for(i in 2:length(breaks)){
    # go over all breaks
    blow  = breaks[i-1]
    bhigh = breaks[i]
    curr_ind = inds[[i-1]]
    curr_ind_old = old_inds[curr_ind]
    for(j in 1:count_bins){
      # for each break point create a bin for each old bin
      new_inds[[indsI]] = curr_ind[curr_ind_old == j]
      indsI = indsI + 1
      new_row = c()
      for(l in 1:(2*ll)){
        new_row = c(new_row, sr[[l]][j])
      }
      new_row = c(new_row, blow, bhigh)
      cut_mat = rbind(cut_mat, new_row, deparse.level = 0)
    }
  }
  # cancel out first NA row
  dimCM = dim(cut_mat)[1]
  df_mat = as.data.frame(cut_mat)
  df_mat = df_mat[2:(dimCM),]
  #create new sr
  new_sr = df_mat
  new_sr$"local_index" = new_inds
  new_sr$"dim_split_history" = sr$dim_split_history
  return(new_sr)
}
transformCategoricalData = function(x, sr=NULL){
  n = length(x)
  categories = unique(sort(x))
  if(length(categories) < 1){
    stop("Error: transformCategoricalData -- less than one category.")
  }
  inds = list()
  breaks = c(0)
  for(i in 1:length(categories)){
    inds[[i]] = which(x == categories[i])
    breaks = c(breaks, i)
  }
  # fixing first dimension issues
  srWasNull = F
  if(is.null(sr)){
    srWasNull = T
    sr = list(V1=c(0), V2=c(1), dim_split_history = c(0), local_index=list())
    sr$local_index[[1]] = 1:n
  }
  new_sr = transformSR(sr=sr, breaks=breaks, inds=inds, srSkip=srWasNull)
  scp = getSC(new_sr)
  new_sr$"L" = scp$L
  new_sr$"R" = scp$R
  new_sr$"L1" = scp$L
  new_sr$"R1" = scp$R
  return(new_sr)
}
#### assumes that both are only 1D discretizations
mergeSR = function(sr1, sr2){
  breaks = c(sr2[[1]], sr2[[2]][dim(sr2)[1]])
  inds = sr2$local_index
  return(transformSR(sr1[1:4], breaks, inds))
}
transformSR = function(sr, breaks, inds, srSkip=F){
  ### first transform sr
  ll = (length(sr) - 2) / 2
  if(srSkip){
    ll = 0
  }
  cut_mat <- t(rep(NA, 2*(ll+1)))
  new_inds = list()
  indsI = 1
  if(length(breaks) < 2){
    stop("Error in transformSR: Less than two break points")
  }
  for(i in 2:length(breaks)){
    # go over all breaks
    blow  = breaks[i-1]
    bhigh = breaks[i]
    for(j in 1:length(sr[[1]])){
      # for each break point create a bin for each old bin
      dps = sr$local_index[[j]]
      new_inds[[indsI]] = dps[ dps %in% inds[[i-1]] ]
      indsI = indsI + 1
      new_row = c()
      if(ll > 0){
        for(l in 1:(2*ll)){
          new_row = c(new_row, sr[[l]][j])
        }
      }
      new_row = c(new_row, blow, bhigh)
      cut_mat = rbind(cut_mat, new_row, deparse.level = 0)
    }
  }
  # cancel out first NA row
  dimCM = dim(cut_mat)[1]
  df_mat = as.data.frame(cut_mat)
  df_mat = df_mat[2:(dimCM),]
  # create new sr
  new_sr = df_mat
  new_sr$"local_index" = new_inds
  new_sr$"dim_split_history" = sr$dim_split_history
  return(new_sr)
}
transformSRMixture = function(sr, v1, v2, inds, srSkip=F){
  ### first transform sr
  ll = (length(sr) - 2) / 2
  if(srSkip){
    ll = 0
  }
  cut_mat <- t(rep(NA, 2*(ll+1)))
  new_inds = list()
  indsI = 1
  if(length(v1) != length(v2) | length(v1) < 1){
    stop("Error in transformSRMixture: Less than two break points")
  }
  for(i in 1:length(v1)){
    # go over all breaks
    blow  = v1[i]
    bhigh = v2[i]
    for(j in 1:length(sr[[1]])){
      # for each break point create a bin for each old bin
      dps = sr$local_index[[j]]
      new_inds[[indsI]] = dps[ dps %in% inds[[i]] ]
      indsI = indsI + 1
      new_row = c()
      if(ll > 0){
        for(l in 1:(2*ll)){
          new_row = c(new_row, sr[[l]][j])
        }
      }
      new_row = c(new_row, blow, bhigh)
      cut_mat = rbind(cut_mat, new_row, deparse.level = 0)
    }
  }
  # cancel out first NA row
  dimCM = dim(cut_mat)[1]
  df_mat = as.data.frame(cut_mat)
  df_mat = df_mat[2:(dimCM),]
  #create new sr
  new_sr = df_mat
  new_sr$"local_index" = new_inds
  new_sr$"dim_split_history" = sr$dim_split_history
  return(new_sr)
}

summarize = function(sr){
  counts = sapply(sr$local_index, length)
  dummy = sr
  dummy$local_index = counts
  return(dummy)
}
logg = function(x){
    if(x == 0){
        return(0)
    }else{
        return(log2(x))
    }
}
log2fac = function(n){
    sum = 0
    for(i in 2:n){
        sum = sum + logg(i)
    }
    return(sum)
}
log2nChoosek = function(n, k){
    if(k > n | k == 0){
        return(0)
    }else{
        return(log2fac(n) - log2fac(k) - log2fac(n-k))
    }
}
getCountsX = function(sr, ind){
  df = data.frame(sr[ind])
  df.un = unique(df)
  new_inds = list()
  countsX = rep(0, dim(df.un)[1])
  for(i in 1:dim(df.un)[1]){
    all_cols = c()
    if(dim(df.un)[2] > 1){
      all_cols = which(apply(df, 1, function(x) identical(unlist(x), unlist(df.un[i,]))))
    }else{
      all_cols = which(df[,1] == df.un[i,1])
      countsX[i] = length(unlist(sr$local_index[all_cols]))
    }
  }
  return(countsX)
}
getLH = function(counts,n){
  notZero = counts != 0
  lh = sum(-counts[notZero] * log2(counts[notZero] / n))
  return(lh)
}
# expects sr for 2 dim, no L, R, L1, R1
# expects domain of X dX and of Y dY
subScore = function(sr, dX, dY){
  ## extract dim 1 (X)
  countsX = getCountsX(sr, 1)
  countsY = getCountsX(sr, 3)
  countsXY = getCounts(sr)
  n = sum(countsX)
  lhX = getLH(countsX,n)
  lhY = getLH(countsY,n)
  lhXY = getLH(countsXY,n)
  score = lhX + lhY - lhXY + regret(k=dX, n=n) + regret(k=dY, n=n) - regret(k=dX*dY, n=n)
  return(score)
}
extractConditionals = function(sr, indX, indY, indsZ){
  #### get domain X
  dX = length(unique(sort(unlist(sr[indX * 2]))))
  dY = length(unique(sort(unlist(sr[indY * 2]))))
  #### first extract conditional Z
  inds = c(indsZ)
  sec_inds = c(inds * 2)
  fir_inds = c(sec_inds - 1)
  df = data.frame(sr[fir_inds])
  df.un = unique(df)
  new_inds = list()
  score = 0
  for(i in 1:dim(df.un)[1]){
    all_cols = c()
    if(dim(df.un)[2] > 1){
      all_cols = which(apply(df, 1, function(x) identical(unlist(x), unlist(df.un[i,]))))
    }else{
      all_cols = which(df[,1] == df.un[i,1])
    }
    xysrcols = 1:(dim(sr)[2]-4) ### L, R, L1, R1
    xysrcols = xysrcols[!(xysrcols %in% c(fir_inds,sec_inds))]
    ## extract sub sr over x and y for Z = z
    sub_sr = (sr[xysrcols])[all_cols,]
    ## compute score of sub sr (SC(X) + SC(Y) - SC(X,Y))
    score = score + subScore(sub_sr, dX, dY)
  }
  return(score)
}
extractFNML = function(sr, indX, indY, indsZ, logE=F){
  #### get domain X
  dX = length(unique(sort(unlist(sr[indX * 2]))))
  dY = length(unique(sort(unlist(sr[indY * 2]))))
  #### first extract conditional Z
  inds = c(indsZ)
  sec_inds = c(inds * 2)
  fir_inds = c(sec_inds - 1)
  df = data.frame(sr[fir_inds])
  df.un = unique(df)
  new_inds = list()
  reg = 0
  for(i in 1:dim(df.un)[1]){
    all_cols = c()
    if(dim(df.un)[2] > 1){
      all_cols = which(apply(df, 1, function(x) identical(unlist(x), unlist(df.un[i,]))))
    }else{
      all_cols = which(df[,1] == df.un[i,1])
    }
    sum_inds = length(unlist(sr$local_index[all_cols]))
    if(sum_inds == 0){
      next
    }
    if(logE){
      reg = reg + (regret(k=dX,n=sum_inds)/log2(exp(1)) + regret(k=dY,n=sum_inds)/log2(exp(1)) - regret(k=dX*dY,n=sum_inds))/log2(exp(1))
    }else{
      reg = reg + (regret(k=dX,n=sum_inds) + regret(k=dY,n=sum_inds) - regret(k=dX*dY,n=sum_inds))
    }
  }
  return(reg)
}
extractScores = function(sr, indX, indY, indsZ, logE=F){
  #### get domain X
  dX = length(unique(sort(unlist(sr[indX * 2]))))
  dY = length(unique(sort(unlist(sr[indY * 2]))))
  #### first extract conditional Z
  reg = 0
  n = 0
  dZ = 1
  if(length(indsZ) == 0){
    n = length(unlist(sr$local_index))
    if(logE){
      reg = reg + (regret(k=dX,n=n)/log2(exp(1)) + regret(k=dY,n=n)/log2(exp(1)) - regret(k=dX*dY,n=n))/log2(exp(1))
    }else{
      reg = reg + (regret(k=dX,n=n) + regret(k=dY,n=n) - regret(k=dX*dY,n=n))
    }
  }else{
    inds = c(indsZ)
    sec_inds = c(inds * 2)
    fir_inds = c(sec_inds - 1)
    df = data.frame(sr[fir_inds])
    df.un = unique(df)
    ### get dZ
    dZ = dim(df.un)[1]
    new_inds = list()
    for(i in 1:dim(df.un)[1]){
      all_cols = c()
      if(dim(df.un)[2] > 1){
        all_cols = which(apply(df, 1, function(x) identical(unlist(x), unlist(df.un[i,]))))
      }else{
        all_cols = which(df[,1] == df.un[i,1])
      }
      sum_inds = length(unlist(sr$local_index[all_cols]))
      n = n + sum_inds
      if(sum_inds == 0){
        next
      }
      if(logE){
        reg = reg + (regret(k=dX,n=sum_inds)/log2(exp(1)) + regret(k=dY,n=sum_inds)/log2(exp(1)) - regret(k=dX*dY,n=sum_inds))/log2(exp(1))
      }else{
        reg = reg + (regret(k=dX,n=sum_inds) + regret(k=dY,n=sum_inds) - regret(k=dX*dY,n=sum_inds))
      }
    }
  } ### end if
  ### Determine which log to use
  regQ = 0
  if(length(indsZ) == 0){
    regQ = reg
  }else{
    if(logE){
      regQ = regret(k=dX*dZ,n=n)/log2(exp(1)) + regret(k=dY*dZ,n=n)/log2(exp(1)) - regret(k=dX*dY*dZ,n=n)/log2(exp(1)) - regret(k=dZ,n=n)/log2(exp(1))
    }else{
      regQ = regret(k=dX*dZ,n=n) + regret(k=dY*dZ,n=n) - regret(k=dX*dY*dZ,n=n) - regret(k=dZ,n=n)
    }
  }
  l = (dX-1)*(dY-1)*dZ
  myLog = log2
  if(logE){
    myLog = log
  }
  jic = -(l / 2) * myLog(n)
  chisq99 = - qchisq(df=l, p=0.99) / 2
  chisq95 = - qchisq(df=l, p=0.95) / 2
  return(list(qNML=regQ, fNML=reg, JIC=jic, Chisq99=chisq99, Chisq95=chisq95))
}
extractSC = function(sr, inds, noVol=F, ll=F){
  # select the first column of each bin related to the given dimensions
  sec_inds = c(inds * 2)
  fir_inds = c(sec_inds - 1)
  df = data.frame(sr[fir_inds])
  df.un = unique(df)
  new_inds = list()
  for(i in 1:dim(df.un)[1]){
    all_cols = c()
    if(dim(df.un)[2] > 1){
      all_cols = which(apply(df, 1, function(x) identical(unlist(x), unlist(df.un[i,]))))
    }else{
      all_cols = which(df[,1] == df.un[i,1])
    }
    sum_inds = unlist(sr$local_index[all_cols])
    new_inds[[i]] = sum_inds
  }
  all_inds = sort(c(fir_inds, sec_inds))
  all_rows = as.numeric(row.names(df.un))
  df2 = data.frame(sr[all_inds])
  all_rows_old = as.numeric(row.names(df2))
  new_sr = df2[all_rows_old %in% all_rows,]
  new_sr$"local_index" = new_inds
  new_sr$"dim_split_history" = -1
  if(ll){
    return(getSCLog(new_sr))
  }else{
    return(getSC(new_sr, noVol=noVol))
  }
}
duplicated_row_indices = function(a){
  # a is a matrix
  ind_base = which(!duplicated(a)) # rows indices that are unique
  ind_dup = setdiff(1:nrow(a), ind_base) # row indices that are equal to one of rows of ind_base
  if(length(ind_dup) == 0){
    return(list(ind_base = ind_base, ind_dup = ind_dup, corresponding_ind = numeric(0)))
  }
  corresponding_ind = rep(0, length(ind_dup))
  for(i in 1:length(ind_dup)){
    corresponding_ind[i] = which(apply(a[ind_base,,drop=F],1,function(b){isTRUE(all.equal(b,a[ind_dup,,drop=F][i,]))})  )
  }
  return(list(ind_base = ind_base, ind_dup = ind_dup, corresponding_ind = corresponding_ind))
}

# random smooth functions
random_fun = function(x){
  dice = sample.int(4,1)
  if(dice == 1){
    return(x)
  } else if(dice == 2){
    return(x^2)
  } else if(dice == 3){
    return(x^3)
  } else if(dice == 4){
    return(tanh(x))
  } else{
    stop("error: random_fun!")
  }
}

# choose the order of dimensions to be split by the compression rate
dim_order_by_compression_rate = function(data, isCat = c()){
  continuous_dimensions = setdiff(1:ncol(data), isCat)
  
  if(length(continuous_dimensions) == 0){
    stop("error: there is no continuous dimension!")
  }
  
  compression_rate = rep(0, ncol(data)) 
  
  for(i in 1:length(continuous_dimensions)){
    dim = continuous_dimensions[i]
    res = multi_hist_splitting_seed_based_simple(data[,dim,drop=F])
    compression_rate[dim] = (res$L[1] + res$R[1]) / (res$L1[1] + res$R1[1])
  }
  
  return(order(compression_rate))
}

res_after_drop_one_dim = function(res, jj, data){
  #jj = 1
  sub_cols_to_remove = c(2*jj-1, 2*jj) 
  sub_cols = setdiff(1:(2*ncol(data)), sub_cols_to_remove)
  duplicated_info_list = res[,sub_cols] %>% as.matrix() 
  duplicated_info_list = duplicated_row_indices(duplicated_info_list)
  
  if(length(duplicated_info_list$ind_dup) == 0){
    res_drop = res[,-sub_cols_to_remove,drop=F] # no duplicated rows to remove
  } else{
    res_drop = res[duplicated_info_list$ind_base,-sub_cols_to_remove, drop=F] # the res if we remove the split of one dimension
    for(i in 1:length(duplicated_info_list$ind_dup)){
      res_drop[duplicated_info_list$corresponding_ind[i],"local_index"][[1]] = 
        c(unlist(res_drop[duplicated_info_list$corresponding_ind[i],"local_index"]), 
          unlist(res[duplicated_info_list$ind_dup[i],"local_index"])) %>% list()
    }
  }
  
  sc = getSC(res_drop[,1:(ncol(res_drop) - 4)])
  res_drop$L = sc$L
  res_drop$R = sc$R
  return(res_drop)
}

