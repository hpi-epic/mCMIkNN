## This file contains the functions to compute our CMI estimate based on adaptive histograms learned from the data. See paper cited in README.
#!! Important: All estimates approximate n*CMI to get the true CMI estimate, the outputs have to be divided by the number of samples n.
## CMI.XX where XX is not equal to "estimates" outputs a CMI estimate w.r.t. a certain correction criterium, which are all explained in the paper. CMI.estimates outputs a list containing the empirical estimate based on the histogram and each of the correction terms; which has to be added to compute the corresponding corrected score. The remaining functions CMI.XX compute: 'LH' the likelihood estimate, 'JIC' this JIC correction (see Suzuki, 2016, Entropy, 'An estimator for mutual information and its application to independence testing') 'qNML' the quotient NML correction, 'fNML' the factorized NML correction, 'Chisq99' the chi-squared correction with alpha=0.01 and 'Chisq95' the chi-squared correction with alpha=0.05.

## Input parameters
# @data is a data frame containing the data
# @xind is the index of X in the data frame
# @yind is the index of X in the data frame
# @zinds is a vector that indicates the positions of the conditional; this can be empty, like "c()" or "c(2)", "c(3,2,1)"...
# @eps and @Kmax are two initialization parameters for the discretization as explained in the paper; using 0 for both initiates the default
# @v is for verbose output: if true it will output some intermediate results
# @isCat is a vector containing the indices of variables that are categorical
CMI.estimates = function(data, xind, yind, zinds, eps=0, Kmax=0, v=F, isCat=c(), noVol=T, logE=T, number_cores_use = 2){
    max_ups = 5
    if(length(xind) != 1 | length(yind) != 1){
      stop("Error testIndep: X and Y should only be one column!")
    }
    result = NULL
    if(length(zinds) > 0){
      inds = c(xind,yind,zinds)
      resXYZ = iterative_cmi_greedy_flexible_parallel(data=as.matrix(data[,inds]), isCat=c(1:length(inds))[c(inds %in% isCat)], eps = eps, Kmax = Kmax, max_num_updating=max_ups, cores=4)
      xsplit = 1
      ysplit = 2
      zsplits = 3:(length(inds))
      scXYZ = getSC(resXYZ[1:(dim(resXYZ)[2]-4)], noVol=noVol)
      scXZ = extractSC(resXYZ, c(xsplit,zsplits), noVol=noVol)
      scYZ = extractSC(resXYZ, c(ysplit,zsplits), noVol=noVol)
      scZ = extractSC(resXYZ, c(zsplits), noVol=noVol)
      lh = scXZ$L + scYZ$L - scXYZ$L - scZ$L
      if(logE){
        lh = scXZ$L/log2(exp(1)) + scYZ$L/log2(exp(1)) - scXYZ$L/log2(exp(1)) - scZ$L/log2(exp(1))
      }
      corr = extractScores(sr=resXYZ, indX = xsplit, indY = ysplit, indsZ=zsplits,logE=logE)
      if(v){
        print(paste(c("Z: l = ", scZ$L, " r = ", scZ$R), collapse = ""))
        print(paste(c("XZ: l = ", scXZ$L, " r = ", scXZ$R), collapse = ""))
        print(paste(c("YZ: l = ", scYZ$L, " r = ", scYZ$R), collapse = ""))
        print(paste(c("XYZ: l = ", scXYZ$L, " r = ", scXYZ$R), collapse = ""))
        print(paste(c("JIC = ", corr$JIC, " Rq = ", corr$qNML, " Rf = ", corr$fNML, " Chisq99 = ", corr$Chisq99, " Chisq95 = ", corr$Chisq95), collapse = ""))
        print(paste(c("l = ", lh), collapse = ""))
      }
      numCont = 2 + length(zinds) - length(c(1:length(inds))[c(inds %in% isCat)])
      result = list(LH=lh, JIC=corr$JIC, qNML=corr$qNML, fNML=corr$fNML, Chisq99=corr$Chisq99, Chisq95=corr$Chisq95)
    }else{
      inds = c(xind,yind)
      resXY = iterative_cmi_greedy_flexible_parallel(data=as.matrix(data[,inds]), eps=eps, Kmax=Kmax, isCat=c(1:length(inds))[c(inds %in% isCat)], 
                                                      max_num_updating=max_ups, cores=number_cores_use)
      xsplit = 1
      ysplit = 2
      scXY = getSC(resXY[1:(dim(resXY)[2]-4)], noVol=noVol)
      scX = extractSC(resXY, c(xsplit), noVol=noVol,ll=F)
      scY = extractSC(resXY, c(ysplit), noVol=noVol,ll=F)
      lh = scX$L + scY$L - scXY$L
      if(logE){
        lh = scX$L/log2(exp(1)) + scY$L/log2(exp(1)) - scXY$L/log2(exp(1))
      }
      corr = extractScores(sr=resXY, indX = xsplit, indY = ysplit, indsZ=c(),logE=logE)
      if(v){
        print(paste(c("X: l = ", scX$L, " r = ", scX$R), collapse = ""))
        print(paste(c("Y: l = ", scY$L, " r = ", scY$R), collapse = ""))
        print(paste(c("XY: l = ", scXY$L, " r = ", scXY$R), collapse = ""))
        print(paste(c("JIC = ", corr$JIC, " Rq = ", corr$qNML, " Rf = ", corr$fNML, " Chisq99 = ", corr$Chisq99, " Chisq95 = ", corr$Chisq95), collapse = ""))
          print(paste(c("l = ", lh), collapse = ""))
        }
        numCont = 2 - length(c(1:length(inds))[c(inds %in% isCat)])
        result = list(LH=lh, JIC=corr$JIC, qNML=corr$qNML, fNML=corr$fNML, Chisq99=corr$Chisq99, Chisq95=corr$Chisq95)
    }
    return(result)
  }

# The below estimates estimate a particular CMI correction.
CMI.LH = function(data, xind, yind, zinds, eps=0, Kmax=0, v=F, isCat=c(), noVol=T, logE=T, number_cores_use = 2){
  res = CMI.estimates(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  return(res$LH)
}
CMI.qNML = function(data, xind, yind, zinds, eps=0, Kmax=0, v=F, isCat=c(), noVol=T, logE=T, number_cores_use = 2){
  res = CMI.estimates(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  return(res$LH + res$qNML)
}
CMI.fNML = function(data, xind, yind, zinds, eps=0, Kmax=0, v=F, isCat=c(), noVol=T, logE=T, number_cores_use = 2){
  res = CMI.estimates(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  return(res$LH + res$fNML)
}
CMI.Chisq99 = function(data, xind, yind, zinds, eps=0, Kmax=0, v=F, isCat=c(), noVol=T, logE=T, number_cores_use = 2){
  res = CMI.estimates(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  return(res$LH + res$Chisq99)
}
CMI.Chisq95 = function(data, xind, yind, zinds, eps=0, Kmax=0, v=F, isCat=c(), noVol=T, logE=T, number_cores_use = 2){
  res = CMI.estimates(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  return(res$LH + res$Chisq95)
}
  
