#### Each function CMIp.XX computes a pseudo p-value based on the output of CMI.XX. Such a pseudo p-value is > 0.01 iff our estimate is <= 0, i.e., indicates independence and <= 0.01 otherwise. With these wrapper functions, our CMI-based test can be applied in algorithms of the pcalg R package such as the PC algorithm to do causal discovery. For an example see "experiments/causal_graph_learning/syntheticNonLinear.R".
# @x index of x
# @y index of y
# @S vector of indicies of the separating set (can be empty)
# @suffStat a list that contains "dm", the data matrix (x is e.g. the index of x in the dm, ...), "type" a vector that idicates for each column in dm if the corresponding feature is discrete (0) or continuous (1) and 'n' the number of data points
CMIp.qNML = function(x,y,S,suffStat){
  if(!is.numeric(x) | !is.numeric(y)){
    stop("CMIp: x and y have to be numeric!")
  }
  if(is.null(suffStat$dm) | is.null(suffStat$type)){
    stop("CMIp: suffStat has to contain a data matrix 'dm' and a vector 'type' indicating the type of the variable for the corresponding column!")
  }
  if(length(S) == 0){
    S = vector()
  }
  if(!is.vector(S) | is.list(S)){
    stop("CMIp: S hast to be a vector!")
  }
  # extract dimensions form dm
  ncols = dim(suffStat$dm)[2]
  # check if indicies are in bounds
  outOfBounds = FALSE
  for(s in S){
    if(s < 1 | s > ncols){
      outofBounds = TRUE
      break
    }
  }
  if((x < 1 | x > ncols) | (y < 1 | y > ncols) | outOfBounds){
    stop("CMIp: The indicies x, y and the set of indicies S have to be between 1 and the number of columns of 'dm'!")
  }
  if( length(suffStat$type) != ncols | (sum(suffStat$type == 0) + sum(suffStat$type == 1) != ncols) ){
    stop("CMIp: 'type' has to be defined for each column in 'dm' and has to be either '0' for discrete or '1' for continuous!")
  }
  # extracting the categorical variables
  isCat = c(x,y,S)
  isCat = isCat[suffStat$type[isCat] == 0]
  val = CMI.qNML(data=suffStat$dm, xind=x, yind=y, zinds=S, isCat=isCat) / suffStat$n
  # shift value such that pv > 0.01 for indep and <= 0.01 for rejecting independence
  val = max(0, val)
  pv = 2^(-6.643855 - val)
  pv = min(pv, 1)
  pv = max(pv, 0)
  return(pv)
}

CMIp.fNML = function(x,y,S,suffStat){
  if(!is.numeric(x) | !is.numeric(y)){
    stop("CMIp: x and y have to be numeric!")
  }
  if(is.null(suffStat$dm) | is.null(suffStat$type)){
    stop("CMIp: suffStat has to contain a data matrix 'dm' and a vector 'type' indicating the type of the variable for the corresponding column!")
  }
  if(length(S) == 0){
    S = vector()
  }
  if(!is.vector(S) | is.list(S)){
    stop("CMIp: S hast to be a vector!")
  }
  # extract dimensions form dm
  ncols = dim(suffStat$dm)[2]
  # check if indicies are in bounds
  outOfBounds = FALSE
  for(s in S){
    if(s < 1 | s > ncols){
      outofBounds = TRUE
      break
    }
  }
  if((x < 1 | x > ncols) | (y < 1 | y > ncols) | outOfBounds){
    stop("CMIp: The indicies x, y and the set of indicies S have to be between 1 and the number of columns of 'dm'!")
  }
  if( length(suffStat$type) != ncols | (sum(suffStat$type == 0) + sum(suffStat$type == 1) != ncols) ){
    stop("CMIp: 'type' has to be defined for each column in 'dm' and has to be either '0' for discrete or '1' for continuous!")
  }
  # extracting the categorical variables
  isCat = c(x,y,S)
  isCat = isCat[suffStat$type[isCat] == 0]
  val = CMI.fNML(data=suffStat$dm, xind=x, yind=y, zinds=S, isCat=isCat) / suffStat$n
  # shift value such that pv > 0.01 for indep and <= 0.01 for rejecting independence
  val = max(0, val)
  pv = 2^(-6.643855 - val)
  pv = min(pv, 1)
  pv = max(pv, 0)
  return(pv)
}

CMIp.Chisq99 = function(x,y,S,suffStat){
  if(!is.numeric(x) | !is.numeric(y)){
    stop("CMIp: x and y have to be numeric!")
  }
  if(is.null(suffStat$dm) | is.null(suffStat$type)){
    stop("CMIp: suffStat has to contain a data matrix 'dm' and a vector 'type' indicating the type of the variable for the corresponding column!")
  }
  if(length(S) == 0){
    S = vector()
  }
  if(!is.vector(S) | is.list(S)){
    stop("CMIp: S hast to be a vector!")
  }
  # extract dimensions form dm
  ncols = dim(suffStat$dm)[2]
  # check if indicies are in bounds
  outOfBounds = FALSE
  for(s in S){
    if(s < 1 | s > ncols){
      outofBounds = TRUE
      break
    }
  }
  if((x < 1 | x > ncols) | (y < 1 | y > ncols) | outOfBounds){
    stop("CMIp: The indicies x, y and the set of indicies S have to be between 1 and the number of columns of 'dm'!")
  }
  if( length(suffStat$type) != ncols | (sum(suffStat$type == 0) + sum(suffStat$type == 1) != ncols) ){
    stop("CMIp: 'type' has to be defined for each column in 'dm' and has to be either '0' for discrete or '1' for continuous!")
  }
  # extracting the categorical variables
  isCat = c(x,y,S)
  isCat = isCat[suffStat$type[isCat] == 0]
  val = CMI.Chisq99(data=suffStat$dm, xind=x, yind=y, zinds=S, isCat=isCat) / suffStat$n
  # shift value such that pv > 0.01 for indep and <= 0.01 for rejecting independence
  val = max(0, val)
  pv = 2^(-6.643855 - val)
  pv = min(pv, 1)
  pv = max(pv, 0)
  return(pv)
}

CMIp.Chisq95 = function(x,y,S,suffStat){
  if(!is.numeric(x) | !is.numeric(y)){
    stop("CMIp: x and y have to be numeric!")
  }
  if(is.null(suffStat$dm) | is.null(suffStat$type)){
    stop("CMIp: suffStat has to contain a data matrix 'dm' and a vector 'type' indicating the type of the variable for the corresponding column!")
  }
  if(length(S) == 0){
    S = vector()
  }
  if(!is.vector(S) | is.list(S)){
    stop("CMIp: S hast to be a vector!")
  }
  # extract dimensions form dm
  ncols = dim(suffStat$dm)[2]
  # check if indicies are in bounds
  outOfBounds = FALSE
  for(s in S){
    if(s < 1 | s > ncols){
      outofBounds = TRUE
      break
    }
  }
  if((x < 1 | x > ncols) | (y < 1 | y > ncols) | outOfBounds){
    stop("CMIp: The indicies x, y and the set of indicies S have to be between 1 and the number of columns of 'dm'!")
  }
  if( length(suffStat$type) != ncols | (sum(suffStat$type == 0) + sum(suffStat$type == 1) != ncols) ){
    stop("CMIp: 'type' has to be defined for each column in 'dm' and has to be either '0' for discrete or '1' for continuous!")
  }
  # extracting the categorical variables
  isCat = c(x,y,S)
  isCat = isCat[suffStat$type[isCat] == 0]
  val = CMI.Chisq95(data=suffStat$dm, xind=x, yind=y, zinds=S, isCat=isCat) / suffStat$n
  # shift value such that pv > 0.01 for indep and <= 0.01 for rejecting independence
  val = max(0, val)
  pv = 2^(-6.643855 - val)
  pv = min(pv, 1)
  pv = max(pv, 0)
  return(pv)
}
