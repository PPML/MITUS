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
  save(StartVal, file=paste("StartVal_", Sys.Date(),".rData", sep=""))
}

