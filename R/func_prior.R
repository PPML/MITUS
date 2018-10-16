#'@name lPrior2 calculates the parameter prior with Jacobians
#'@param Par parameters translated to the normal space
#'@param Par3 parameters in their original distributions
#'@return ldensity3

lPrior2 <- function(Par,Par3) {
  if(dim(as.matrix(Par))[2]==1) Par <- t(as.matrix(Par))
  ldensity <- dmnorm(Par,rep(0,nrow(ParamInitZ)),diag(nrow(ParamInitZ)),log=T)
  ldensity2 <- ldensity-sum(dnorm(Par,0,1,log=T))
  lDensTrue <- rep(NA,nrow(ParamInitZ))
  lDensTrue[idZ0] <- dbeta( Par3[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7],log=T)
  lDensTrue[idZ1] <- dgamma(Par3[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7],log=T)
  lDensTrue[idZ2] <- dnorm( Par3[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7],log=T)
  ldensity3 <- ldensity2+sum(lDensTrue)
  ldensity3  }

#'@name lprior
#'@param ParMatrix matrix of Initial parameters
#'@return lPri

lprior <- function(ParMatrix = ParInit) { # Par = ParInit
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    ParMatrix <- t(as.data.frame(ParMatrix)) }
  lPri <- rep(0,nrow(ParMatrix))
  for(samp_i in 1:nrow(ParMatrix)) {
    # norm2unif
    Par <- ParMatrix[samp_i,]
    Par2 <- pnorm(Par,0,1)
    # unif2true
    Par3 <- Par2
    Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
    Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
    Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
    lPri[samp_i] <- lPrior2(Par,Par3)
    if(is.na(lPri[samp_i]))         {
      lPri[samp_i] <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
    if(lPri[samp_i]%in%c(-Inf,Inf)) {
      lPri[samp_i] <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  } # end inner loop
  return(lPri)
}
