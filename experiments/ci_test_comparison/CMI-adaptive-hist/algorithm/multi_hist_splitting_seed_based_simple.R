# @data is a data frame (one or more columns)
# @eps is the minimum bin width
# @Kmax it the maximum number of bins
# @res can be NULL (by default) or a already created result file outputted by this function; which is then used as a seed (e.g. X and Y are already discretized)
# @isCat is a vector of booleans indicating which variables are categorical (default is NULL, indicating so all variables are continuous)
multi_hist_splitting_seed_based_simple = function(data, eps=0, Kmax=0, res=NULL, isCat=NULL){
  # This function discretizes the data iteratively until the "local"
  #   code lengths reaches the minimum. 
  
  # init
  dimD = dim(data)[2]
  n.dp = dim(data)[1]
  remaining_dims = 1:dimD
  cat_dims = c()
  # handle categorical data
  if(!is.null(isCat)){
    if(length(isCat) != length(remaining_dims)){
      stop("Error: multi_hist_splitting_seed_based_simple.R -- isCat cannot contain a different number of elements than there are rows!")
    }else{
      cat_dims = remaining_dims[isCat]
    }
  }
  history = c()
  
  while(length(remaining_dims) > 0){
    # discretize next dimension
    new_res = NULL
    best_res = -1
    dummy_res = NULL
    if(!is.null(res)){
      dummy_res = res[1:(dim(res)[2]-4)] ## L,R,L1,R1 are not needed
    }
    for(d in remaining_dims){
      curr_res = NULL
      if(d %in% cat_dims){
        curr_res = transformCategoricalData(x=data[,d], sr=dummy_res)
      }else{
        # calculate current eps and Kmax
        x.not_disc = data[,d]
        x.range = max(x.not_disc) - min(x.not_disc)
        current_eps = x.range / (20 * log(n.dp))
        if(eps > 0){
          current_eps = eps
        }
        current_Kmax = ceiling(2 * log(n.dp))
        if(Kmax > 0){
          current_Kmax = Kmax
        }
        curr_res = oned_hist_iterative(x=x.not_disc, eps=current_eps, Kmax=current_Kmax, sr=dummy_res)
      }
      if(is.null(new_res)){
        new_res = curr_res
        best_res = d
      }else{
        sc1 = new_res$L[1] + new_res$R[1]
        sc2 = curr_res$L[1] + curr_res$R[1]
        if(sc2 < sc1){
          new_res = curr_res
          best_res = d
        }
      }
    }
    if(best_res == -1 | is.null(new_res)){
      stop("Error: res not updated... -1.")
    }else{
      res = new_res
      remaining_dims = remaining_dims[remaining_dims != best_res]
      history = c(history,best_res)
    }
  }
  
  if(is.null(res)){
    stop("Discretization ended with res being NULL!")
  }else{
    lh = length(history)
    sh = 1
    if(dim(res)[1] < lh){
      sh = 1 + lh - dim(res)[1]
    }
    if(sh == lh){
      res$dim_split_history = history[lh]
    }else if(sh > 1){
      res$dim_split_history = history[sh:lh]
    }else{
      res$dim_split_history[1:lh] = history
    }
    return(res)
  }
}
