oned_hist_iterative = function(x, eps, Kmax, sr = NULL, useAsymptoticRegret=F){
  # This funciton constructs MDL-optimal one-dimensional histogram
  
  # x is the 1d data, continuous
  # eps is the precision used to pre-discreitize the data
  # eps can be 0.1,0.01, etc
  # Kmax is the maximum of bins of the segmentation
  require("SCCI") # used for regret calculation
  
  n = length(x)
  
  ## needed for mixtures
  x.copy = x
  max.x = round(max(c(x,1))) + 1
  disc.vals = c()
  v1d = c()
  v2d = c()
  disc.inds = list()
  remaining.inds = 1:n
  
  # preprocess to get all discrete points
  if(n > 1){
    minpoints = 10
    if(n / 100 < 10){
      minpoints = 5
    }
    x.sort = sort(x)
    val = x.sort[1]
    count = 1
    for(i in 2:n){
      if(x.sort[i] == val){
        count = count + 1
        if(i == n){
          if(count >= minpoints){
            disc.vals = c(disc.vals, val)
          }
        }
      }else{
        if(count >= minpoints){
          disc.vals = c(disc.vals, val)
        }
        count = 1
        val = x.sort[i]
      }
    }
    for(v in disc.vals){
      x = x[x != v]
    }
    ldv = length(disc.vals)
    if(ldv > 0){
      for(i in 1:ldv){
        v1d = c(v1d, max.x + i)
        v2d = c(v2d, max.x + i + 1)
        disc.inds[[i]] = which(x.copy == disc.vals[i])
      }
      all.inds = 1:n
      remaining.inds = all.inds[!(all.inds %in% unlist(disc.inds))]
      # recompute eps and Kmax
      if(eps == 0){
        x.not_disc = x
        x.range = max(x.not_disc) - min(x.not_disc)
        n.dp = length(x)
        eps = x.range / (20 * log(n.dp))
      }
      if(Kmax == 0){
        Kmax = ceiling(2 * log(n.dp))
      }
    }
  }

  C_til_complete = generate_candidate_cuts(x, NULL, eps)
  
  # in case that all the data points are actually the same
  if(length(C_til_complete) == 1){
    C_til_complete = c(C_til_complete, C_til_complete + eps)
  }
  C_til = refine_candidate_cuts(C_til_complete, x) ## cutpoints

  # init default for already discretized dimensions (the default is if it is the first disrectized dimension)
  srWasNull = F
  if(is.null(sr)){
    srWasNull = T
    sr = list(V1=c(0), V2=c(1), dim_split_history = c(0), local_index=list())
    sr$local_index[[1]] = 1:n
  }
  
  # No model selection if no candidate cut points
  if(length(C_til) == 2){
    v1 = c(C_til[1], v1d)
    v2 = c(C_til[2], v2d)
    inds = list()
    inds[[1]] = remaining.inds
    if(ldv > 0){
      for(i in 1:ldv){
        inds[[i+1]] = disc.inds[[i]]
      }
    }
    new_sr = transformSRMixture(sr=sr, v1=v1, v2=v2, inds=inds, srSkip=srWasNull)
    scp = getSC(new_sr, asymp=useAsymptoticRegret)
    new_sr$"L" = scp$L
    new_sr$"R" = scp$R
    new_sr$"L1" = scp$L
    new_sr$"R1" = scp$R
    return(new_sr)
  }
  
  # no model selection if Kmax = 1
  if(Kmax == 1){
    v1 = c(C_til[1], v1d)
    v2 = c(C_til[length(C_til)], v2d)
    inds = list()
    inds[[1]] = remaining.inds
    if(ldv > 0){
      for(i in 1:ldv){
        inds[[i+1]] = disc.inds[[i]]
      }
    }
    new_sr = transformSRMixture(sr=sr, v1=v1, v2=v2, inds=inds, srSkip=srWasNull)
    scp = getSC(new_sr, asymp=useAsymptoticRegret)
    new_sr$"L" = scp$L
    new_sr$"R" = scp$R
    new_sr$"L1" = scp$L
    new_sr$"R1" = scp$R
    return(new_sr)
  }
  
  
  C_0 = C_til[1]                        ##smallest
  C_E_plus_1 = C_til[length(C_til)]     ##largest
  C_til = C_til[-1] #exclude the (x_min - eps/2) but KEEP (x_max + eps/2)
  
  # ne in the paper, also make the last interval close on both sides
  C_til[length(C_til)] = C_til[length(C_til)] + 100
  volumes = getVols(sr) ### volumes from pre-discretized variables
  lv = length(volumes)
  nee = sapply(C_til, function(c){
    curr_inds = which(x < c)
    dummy = rep(0, lv)
    for(i in 1:lv){
      dummy[i] = sum(curr_inds %in% sr$local_index[[i]])
    }
    return(dummy)
  })
  if(lv > 1){
    nee = t(nee)
  }else{
    nee = as.matrix(nee)
  }
  C_til[length(C_til)] = C_til[length(C_til)] - 100
  
  ## Create volumes
  C_vol = matrix(0, nrow=length(C_til), ncol=lv)
  for(i in 1:lv){
    C_vol[,i] = C_til * volumes[i]
  }
  
  # initiate the dynamical programming
  # first row for the dynamic programming: lh until this cut point
  # we do not need to consider the discret points, because they are constant for the optimization
  C_00 = C_0 * volumes
  B_1e = rowSums( -nee * (modified_log(nee) - log2(sweep(C_vol, 2, C_00) * n) ) )
  
  E = length(C_til)
  
  Kmax = min(E,Kmax)
  
  ## accounts for previous bins
  R = sapply(1:Kmax, function(cur_k){
    return(regret(k=(cur_k + ldv) * lv, n=n)) #### ldv is number of discrete variables
  })
  
  B = matrix(rep(Inf, E * Kmax), ncol = Kmax)
  e_prime_min = matrix(rep(Inf, E * Kmax), ncol = Kmax)
  
  B[,1] = B_1e
  
  for(k in 2:Kmax){
    minss = lapply(k:E, function(e){
      eprime_candidates = (k-1):(e-1)
      if(lv == 1){ ## 1D
        counts = nee[e,] - nee[eprime_candidates,]
        candidate = B[eprime_candidates,k-1] - (counts) * (modified_log(counts) - log2((C_vol[eprime_candidates,] - C_vol[e,]) * (-n))) ## -n since we reversed the order of subtraction
        return(list(min(candidate), which.min(candidate)))
      }else if(length(eprime_candidates) < 2){ ## Cannot use sweep or row sum for single element
        counts = nee[e,] - nee[eprime_candidates,]
        candidate = B[eprime_candidates,k-1] - sum( (counts) * (modified_log(counts) - log2((C_vol[eprime_candidates,] - C_vol[e,]) * (-n))) ) ## -n since we reversed the order of subtraction
        return(list(min(candidate), which.min(candidate)))
      }else{
        counts = - sweep(nee[eprime_candidates,], 2, nee[e,]) ## neg cause reversed
        candidate = B[eprime_candidates,k-1] - rowSums((counts) * (modified_log(counts) - log2(sweep(C_vol[eprime_candidates,],2,C_vol[e,]) * (-n))) ) ## -n since we reversed the order of subtraction
        return(list(min(candidate), which.min(candidate)))
      }
    })
    B[k:E,k] = unlist(minss)[seq(1,2*length(minss)-1,2)]
    e_prime_min[k:E,k] = unlist(minss)[seq(2,2*length(minss),2)] + (k-2)
  }
  
  L_model = log2(choose(length(C_til_complete), 0:(Kmax-1)))
  total_length = B[E,] + R + L_model
  mm = which.min(total_length)
  
  breakss = rep(0, mm)
  EE = E
  
  v11 = c(C_0, v1d)
  v21 = c(C_til[length(C_til)], v2d)
  inds1 = list()
  inds1[[1]] = remaining.inds
  if(ldv > 0){
    for(i in 1:ldv){
      inds1[[i+1]] = disc.inds[[i]]
    }
  }
  new_sr1 = transformSRMixture(sr=sr, v1=v11, v2=v21, inds=inds1, srSkip=srWasNull)
  scp1 = getSC(new_sr1, asymp=useAsymptoticRegret)
  new_sr1$"L" = scp1$L
  new_sr1$"R" = scp1$R

  new_sr1$"L1" = scp1$L
  new_sr1$"R1" = scp1$R
  l1 = scp1$L
  r1 = scp1$R
  if(mm == 1){
    return(new_sr1)
  }
  
  # get breakpoints
  for(j in mm:2){
    breakss[j] = e_prime_min[EE,j]
    EE = breakss[j]
  }
  breaks = c(C_0, C_til[breakss],C_til[length(C_til)])
  local_data_index = list()
  # make the last interval closed on both sides.
  v1 =c()
  v2 = c()
  breaks[length(breaks)] = breaks[length(breaks)] + 100
  for(i in 2:length(breaks)){
    v1 = c(v1, breaks[i-1])
    v2 = c(v2, breaks[i])
    local_data_index[[i-1]] = which(x < breaks[i] & x >= breaks[i-1])
  }
  breaks[length(breaks)] = breaks[length(breaks)] - 100
  v2[length(v2)] = v2[length(v2)] - 100
  
  ## add discrete points
  if(ldv > 0){
    linds = length(local_data_index)
    v1 = c(v1, v1d)
    v2 = c(v2, v2d)
    for(i in 1:ldv){
      local_data_index[[i+linds]] = disc.inds[[i]]
    }
  }
  
  # compute score
  new_sr = transformSRMixture(sr=sr, v1=v1, v2=v2, inds=local_data_index, srSkip=srWasNull)
  scp = getSC(new_sr, asymp=useAsymptoticRegret)
  new_sr$"L" = scp$L
  new_sr$"R" = scp$R

  new_sr$"L1" = l1
  new_sr$"R1" = r1
  return(new_sr)
}
