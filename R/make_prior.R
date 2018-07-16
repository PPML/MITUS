################################################################################
##### THIS CODE WAS ORIGINALLY AUTHORED BY DR. NICK MENZIES AND HAS BEEN   #####
##### EDITED BY NICOLE A. B. SWARTWOOD FOR USE ON A PROJECT INVESTIGATING  #####
##### THE XPERT DIAGNOSTIC EFFECTIVENESS IN BOTSWANA IN THE YEAR 2017-2018.#####
################################################################################

################################################################################
#####						PACKAGES REQUIRED				               #####
#####          CODE USES logitnorm, mvtnorm, and lhs PACKAGES              #####
################################################################################
#library(logitnorm)
library(mvtnorm)
#library(lhs)


## Load files

logit <- function(x) log(x/(1-x))

gammapar <- function(tgt) {
  tgt <- as.numeric(tgt)
  mn <- tgt[1]; cir <- (tgt[3]-tgt[2])
  xopt <- function(b,mn,cir) {
    cir2 <- qgamma(c(1,39)/40,mn*b,b); cir2 <- cir2[2]-cir2[1]
    (cir2-cir)^2 }
  zz <- optimize(xopt,c(0.1,100000),mn=mn,cir=cir)$minimum
  c(zz*mn,zz)
}

# Function for calculating beta parameters
betapar <- function(tgt) {
  tgt <- as.numeric(tgt)
  mn <- tgt[1]; cir <- (tgt[3]-tgt[2])
  xopt <- function(xx,mn=mn,cir=cir) {
    cir2 <- qbeta(c(1,39)/40,xx*(mn/(1-mn)),xx); cir2 <- cir2[2]-cir2[1]
    sum((cir2-cir)^2)
  }
  zz <- optimize(xopt,c(0.2,100000),mn=mn,cir=cir)$minimum
  bp <-  c(zz*(mn/(1-mn)),zz)
  if(sum(abs(bp-1))<0.2) {
    c(1,1) } else { bp }
}


########################
### Calc Normal values

ParamInit<-as.data.frame(read.csv("~/MITUS/inst/extdata/ParamInitUS.csv")[,2:6])
rownames(ParamInit)<-read.csv("~/MITUS/inst/extdata/ParamInitUS.csv")[,1]

ParamInit <- cbind(ParamInit,NA,NA)
colnames(ParamInit) <- c(colnames(ParamInit)[1:5],"Par1","Par2")

for (i in 1:nrow(ParamInit)) {
  if(ParamInit[i,4]==0)   { ParamInit[i,6:7] <- betapar(ParamInit[i,1:3]) }
  if(ParamInit[i,4]==1)   { ParamInit[i,6:7] <- gammapar(ParamInit[i,1:3])  }
  if(ParamInit[i,4]==2) 	{ ParamInit[i,6:7] <- c(ParamInit[i,1],(ParamInit[i,3]-ParamInit[i,2])/3.92) } }

### (re)Calculate mean and bounds
ParamInit[ParamInit[,4]==0,1] <- ParamInit[ParamInit[,4]==0,6]/rowSums(ParamInit[ParamInit[,4]==0,6:7])
ParamInit[ParamInit[,4]==1,1] <- ParamInit[ParamInit[,4]==1,6]/ParamInit[ParamInit[,4]==1,7]
ParamInit[ParamInit[,4]==2,1] <- ParamInit[ParamInit[,4]==2,6]

ParamInit[ParamInit[,4]==0,2] <- qbeta( 0.025,ParamInit[ParamInit[,4]==0,6],ParamInit[ParamInit[,4]==0,7])
ParamInit[ParamInit[,4]==0,3] <- qbeta( 0.975,ParamInit[ParamInit[,4]==0,6],ParamInit[ParamInit[,4]==0,7])
ParamInit[ParamInit[,4]==1,2] <- qgamma(0.025,ParamInit[ParamInit[,4]==1,6],ParamInit[ParamInit[,4]==1,7])
ParamInit[ParamInit[,4]==1,3] <- qgamma(0.975,ParamInit[ParamInit[,4]==1,6],ParamInit[ParamInit[,4]==1,7])
ParamInit[ParamInit[,4]==2,2] <- qnorm( 0.025,ParamInit[ParamInit[,4]==2,6],ParamInit[ParamInit[,4]==2,7])
ParamInit[ParamInit[,4]==2,3] <- qnorm( 0.975,ParamInit[ParamInit[,4]==2,6],ParamInit[ParamInit[,4]==2,7])


write.csv(ParamInit,file=paste0("ParamInitUS_",Sys.Date(),"_final.csv"))
save(ParamInit,file=paste0("ParamInitUS_",Sys.Date(),"tab.rda"))
################################################################################
######################### DRAWS SOME STARTING VALUES ###########################
sample.prior2 <- function(n)  { # n=5
  jj <-ParamInit[ParamInit$Calib==1,]
  pari1 <- rmvnorm(n,rep(0,nrow(jj)),diag(nrow(jj)))
  pari1 }

set.seed(321); StartVal <- sample.prior2(n=20)
save(StartVal,file=paste0("StartVal_", Sys.Date(), ".rda"))
################################################################################



