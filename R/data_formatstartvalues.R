#'Function to create a Starting Values for the Parameters
#'to be Optimized.
#'
#'@name gen_st_val
#'@param n number of datasets you want to generate
#'@param samp Which sampling method to use. "NORM" for random normal samples or "LHS" for latin hypercube sampling.
#'@return rData file of starting values
#'@export
gen_st_val <-function(n=10, samp="LHS"){
  if(identical(samp,"LHS")==FALSE & identical(samp,"NORM")==FALSE) stop("'samp' must be either 'LHS' or 'NORM'")
  if(identical(samp,"NORM")==TRUE) StartVal <-sample.prior1(n)
  if(identical(samp,"LHS") ==TRUE) StartVal <- sample.prior2(n)
  colnames(StartVal)<-rownames(ParamInitZ)
  saveRDS(StartVal, file=paste("~/MITUS/inst/US/US_StartVal_", Sys.Date(),".rds", sep=""),version=2)

}

#'@name gen_st_val_st
#'@param n number of datasets you want to generate
#'@param samp Which sampling method to use. "NORM" for random normal samples or "LHS" for latin hypercube sampling.
#'@return rData file of starting values
#'@export
gen_st_val_st <-function(n=10, samp="LHS"){
  if(identical(samp,"LHS")==FALSE & identical(samp,"NORM")==FALSE) stop("'samp' must be either 'LHS' or 'NORM'")
  if(identical(samp,"NORM")==TRUE) StartVal <-sample.prior1(n)
  if(identical(samp,"LHS") ==TRUE) StartVal <- sample.prior2(n)
  colnames(StartVal)<-rownames(ParamInitZ)
  StartVal_st<-StartVal
  saveRDS(StartVal_st, file=paste("~/MITUS/inst/ST/ST_StartVal_", Sys.Date(),".rds", sep=""),version=2)

}


#makes a matrix of the Parameter Vectors generated from StartVals
#or optimized data opt_all
#'@name gen_par_matrix
#'@param startMat a matrix of values in the normal space
#'@return matrix of Params in their original distributions
gen_par_matrix<-function(startMat){
  parAll<-matrix(NA,nrow(startMat),nrow(ParamInit))
  colnames(parAll)<-rownames(ParamInit)
  for(i in 1:nrow(startMat)){
  if(min(dim(as.data.frame(startMat)))==1) {
    Par <- as.numeric(startMat);
    names(Par) <- names(startMat)
  } else {  Par <- as.numeric(startMat[i,]);
  names(Par) <- colnames(startMat) }  ##previously, the distribution of parameters were transformed to normal distribution in
  ##to facilitate comparisons. These first two steps convert these parameters back to their
  ##distributions
  # normal to uniform
  Par2 <- pnorm(Par,0,1)
  # uniform to true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  P <- P
  parAll[i,]<-P
  }
saveRDS(parAll,file=paste("parAll",nrow(startMat),"_",Sys.time(),".rds", sep=""), version=2)
  return(parAll)

}

