generate_candidate_cuts = function(x, outer_bound=NULL, eps){
  # This function generate the search space of cut-points for
  #   constructing one-dimensional MDL histograms
  max_x = max(x)
  min_x = min(x)
  
  if(outer_bound %>% is.null()){
    C_til = seq(floor(min_x / eps), 
                ceiling(max_x / eps), 
                1) * eps
  } else {
    if(outer_bound[1] > min_x | outer_bound[2] < max_x){
      stop("outer_bound must lie outside the interval [min_x,max_x]. \n 
           Restrict the data to a smaller region before running this function if necessary.")
    }
    C_til = seq(floor(outer_bound[1]/eps), ceiling(outer_bound[2]/eps), by = 1) * eps
  }
  
  if(max_x == tail(C_til,1)){
    C_til = c(C_til, tail(C_til,1) + eps)
  }
  
  return(C_til)
}
