#'This script defines the functions to sample the prior.


#'This function samples the prior using a random normal generation.
#'@name sample.prior1
#'@param n number of sample sets to generate
#'@return n sample sets from the prior in a single object
sample.prior1 <- function(n) {
  rmnorm(n,rep(0,sum(ParamInit$Calib==1)),diag(sum(ParamInit$Calib==1))) }
#'This function samples the prior using Latin Hypercube sampling method.
#'@name sample.prior2
#'@param n number of sample sets to generate
#'@return n sample sets from the prior in a single object
sample.prior2 <- function(n) {
  qnorm(randomLHS(n,sum(ParamInit$Calib==1)),0,1)   }
#sample.prior  <- sample.prior2
