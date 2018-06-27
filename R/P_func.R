###this script is to move this code into a place that I can find
### it until I can decide the best place to move it into for final
### model run
#'@name define_P
#'@param ParMatrix
#'@param samp_i
#'@return P a vector of values
define_P <-function (ParMatrix,samp_i){
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par <- as.numeric(ParMatrix);
    names(Par) <- names(ParMatrix)
  } else {
    Par <- as.numeric(ParMatrix[samp_i,])
    names(Par) <- colnames(ParMatrix) }
}

