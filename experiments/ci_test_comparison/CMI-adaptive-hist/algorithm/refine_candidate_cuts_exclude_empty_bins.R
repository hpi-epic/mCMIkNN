refine_candidate_cuts = function(candidate_cuts, x){
  # This function removes those candidate_cuts that split an empty bin 
  # into more empty bins. 
  
  # we always consider interval [ ), except that the last interval is [ ].
  
  ne = sapply(candidate_cuts, function(y){sum(x<y)}) 

  
  counts_each_bin = diff(ne)
  
  index = which(counts_each_bin == 0)[
      which(
        diff(which(counts_each_bin == 0)) == 1
      )
    ]+1
  
  if(length(index)==0){ # in case that index is numeric(0)
    return(candidate_cuts) 
  } else {
    refined_cuts = candidate_cuts[-index]
    return(refined_cuts)
  }
}