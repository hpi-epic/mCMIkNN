modified_log = function(x, base = 2){
  # This is a modified version of $log$ function, in order to set 0log0=0
  
  lg = ifelse(x==0, -.Machine$double.xmax, log(base = base, x))

  return(lg)
}

