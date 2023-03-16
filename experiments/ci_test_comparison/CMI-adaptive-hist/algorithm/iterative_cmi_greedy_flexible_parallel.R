iterative_cmi_greedy_flexible_parallel = function(data, eps = 0, Kmax = 0, v=F, isCat = c(), max_num_updating = 20, cores=1){
  require("parallel")
  cores = min(cores, detectCores()) ## check if number of available cores is sufficient
  # transform the categorical data into 1,2,3,...,(number of categories)
  if(length(isCat)>0){
    data[,isCat] = apply(data[,isCat,drop = F],2, function(x){as.integer(as.factor(x))})
  }
  
  
  # In case it is 1D data, then no iteration is needed. 
  if(ncol(data) == 1){
    if(length(isCat) == 1){
      res = multi_hist_splitting_seed_based_simple(data, isCat = T)
    } else {
      res = multi_hist_splitting_seed_based_simple(data)  
    }
    return(res)
  }
  
  
  # initialize the grid
  prev_res = NULL
  for(i in 1:ncol(data)){
    if(i %in% isCat){
      res = multi_hist_splitting_seed_based_simple(data[,i,drop=F], res = prev_res, isCat = T)
    } else {
      res = multi_hist_splitting_seed_based_simple(data[,i,drop=F], res = prev_res, Kmax = 1)
    }
    prev_res = res
  }
  l1 = res$L1
  r1 = res$R1
  if(length(res$L1) > 1){
    l1 = res$L1[1]
    r1 = res$R1[1]
  }
  
  
  # iteratively updating each dimension in a greedy manner
  min_sc = res$L[1] + res$R[1]
  min_res = NULL
  min_dim = 1:ncol(data)
  remaining_cat = isCat
  dims_in_order = 1:ncol(data)
  for(iter in 1:max_num_updating){
    min_dim_here = 0
    num_updates = 0
    result_list = mclapply(1:(ncol(data)), function(jj){
      # drop the current split of this dimension
      dim = dims_in_order[jj]
      if(dim %in% isCat){
        return(NULL)
      }
      sub_cols_to_remove = c(2*jj-1, 2*jj)
      sub_cols = setdiff(1:(2*ncol(data)), sub_cols_to_remove)
      
      duplicated_info_list = prev_res[,sub_cols] %>%
        as.matrix()
      duplicated_info_list = duplicated_row_indices(duplicated_info_list) # The function duplicated_row_indices() is in utils.R
      
      if(length(duplicated_info_list$ind_dup) == 0){
        prev_res_drop = prev_res[,-sub_cols_to_remove,drop=F] # no duplicated rows to remove
      } else{
        prev_res_drop = prev_res[duplicated_info_list$ind_base,-sub_cols_to_remove, drop=F] # the result if we remove the split of one dimension
        for(i in 1:length(duplicated_info_list$ind_dup)){
          prev_res_drop[duplicated_info_list$corresponding_ind[i],"local_index"][[1]] =
            c(unlist(prev_res_drop[duplicated_info_list$corresponding_ind[i],"local_index"]),
              unlist(prev_res[duplicated_info_list$ind_dup[i],"local_index"])) %>% list()
        }
      }
      # update
      res = multi_hist_splitting_seed_based_simple(data[,dim,drop=F], res = prev_res_drop, isCat = dim %in% isCat, Kmax = Kmax)
      return(res)
    }, mc.cores=cores)
    for(jj in 1:(ncol(data))){
      dim = dims_in_order[jj]
      if(dim %in% isCat){
        next
      }
      num_updates = num_updates + 1
      res = ((result_list[jj]))[[1]]
      if(res$L[1] + res$R[1] < min_sc){ # note that min_sc is not equal to min_sc if we are considering a categorical dimension
        min_res = res # update min_res
        min_sc = res$L[1] + res$R[1] # update min_sc
        min_dim_here = dim
      }
    }
    if(min_dim_here == 0){ # break if further split will not have lower SC 
      break
    }
    min_dim = c(min_dim, min_dim_here) # update the min_dim
    if(min_dim_here %in% isCat){
      remaining_cat = setdiff(remaining_cat, min_dim_here)
    }
    prev_res = min_res # update the prev_res
    dims_in_order = c(setdiff(dims_in_order, min_dim_here),min_dim_here)
    if(num_updates <= 1){
      break
    }
  }
  res = prev_res
  
  # make the res in the right order 
  min_dim = unique(min_dim, fromLast = T)
  correct_order = c(min_dim * 2 - 1, min_dim * 2) %>% matrix(byrow = T, nrow = 2) %>% as.numeric()
  res[,correct_order] = res[,1:(2*ncol(data))]
  res$L1 = l1
  res$R1 = r1
  return(res)
}
