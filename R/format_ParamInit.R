format_paramInit<-function(){
  library(mvtnorm)
#library(lhs)


## Load files

lgt <- function(x) log(x/(1-x))
#'takes a mean and interval and creates a distribution
#'fitting mean and the width of the interval
gammapar <- function(tgt) {
  tgt <- as.numeric(tgt)
  mn <- tgt[1]; cir <- (tgt[3]-tgt[2])
  xopt <- function(b,mn,cir) {
    cir2 <- qgamma(c(1,39)/40,mn*b,b); cir2 <- cir2[2]-cir2[1]
    (cir2-cir)^2 }
  zz <- optimize(xopt,c(0.1,100000),mn=mn,cir=cir)$minimum
  c(zz*mn,zz) }

# Function for calculating beta parameters
betapar <- function(tgt) {
  tgt <- as.numeric(tgt)
  mn <- tgt[1]; cir <- (tgt[3]-tgt[2])
  xopt <- function(xx,mn=mn,cir=cir) {
    cir2 <- qbeta(c(1,39)/40,xx*(mn/(1-mn)),xx); cir2 <- cir2[2]-cir2[1]
    sum((cir2-cir)^2) }
  zz <- optimize(xopt,c(0.2,100000),mn=mn,cir=cir)$minimum
  bp <-  c(zz*(mn/(1-mn)),zz)
  if(sum(abs(bp-1))<0.2) { c(1,1) } else { bp }  }


########################
### Calc Normal values
ParamInit <- as.data.frame(read.csv("~/MITUS/inst/extdata/US/ParamInitUS.csv")[,2:6]); rownames(ParamInit) <- read.csv("~/MITUS/inst/extdata/US/ParamInitUS.csv")[,1]

ParamInit <- cbind(ParamInit,NA,NA,NA)
colnames(ParamInit) <- c(colnames(ParamInit)[1:5],"Par1","Par2","TransMean")
#'calculates the distribution values
for (i in 1:nrow(ParamInit)) {
  if(ParamInit[i,4]==0)   { ParamInit[i,6:7] <- betapar(ParamInit[i,1:3]) }
  if(ParamInit[i,4]==1)   { ParamInit[i,6:7] <- gammapar(ParamInit[i,1:3])  }
  if(ParamInit[i,4]==2) 	{ ParamInit[i,6:7] <- c(ParamInit[i,1],(ParamInit[i,3]-ParamInit[i,2])/3.92) } }

### (re)Calculate mean and bounds
#'might be unnecessary
ParamInit[ParamInit[,4]==0,1] <- ParamInit[ParamInit[,4]==0,6]/rowSums(ParamInit[ParamInit[,4]==0,6:7])
ParamInit[ParamInit[,4]==1,1] <- ParamInit[ParamInit[,4]==1,6]/ParamInit[ParamInit[,4]==1,7]
ParamInit[ParamInit[,4]==2,1] <- ParamInit[ParamInit[,4]==2,6]

ParamInit[ParamInit[,4]==0,2] <- qbeta( 0.025,ParamInit[ParamInit[,4]==0,6],ParamInit[ParamInit[,4]==0,7])
ParamInit[ParamInit[,4]==0,3] <- qbeta( 0.975,ParamInit[ParamInit[,4]==0,6],ParamInit[ParamInit[,4]==0,7])
ParamInit[ParamInit[,4]==1,2] <- qgamma(0.025,ParamInit[ParamInit[,4]==1,6],ParamInit[ParamInit[,4]==1,7])
ParamInit[ParamInit[,4]==1,3] <- qgamma(0.975,ParamInit[ParamInit[,4]==1,6],ParamInit[ParamInit[,4]==1,7])
ParamInit[ParamInit[,4]==2,2] <- qnorm( 0.025,ParamInit[ParamInit[,4]==2,6],ParamInit[ParamInit[,4]==2,7])
ParamInit[ParamInit[,4]==2,3] <- qnorm( 0.975,ParamInit[ParamInit[,4]==2,6],ParamInit[ParamInit[,4]==2,7])

### Calculate transformed
#'calculating the transformed by the log or lgt
ParamInit[ParamInit[,4]==0,8] <- lgt(ParamInit[ParamInit[,4]==0,1])
ParamInit[ParamInit[,4]==1,8] <- log(  ParamInit[ParamInit[,4]==1,1])
ParamInit[ParamInit[,4]==2,8] <-       ParamInit[ParamInit[,4]==2,1]
ParamInit[ParamInit[,8]==-Inf,8] <- -100;  ParamInit[ParamInit[,8]==Inf,8] <- 100
# ParamInit_st<-ParamInit
#'save it all
write.csv(ParamInit,file="ParamInitUS_final.csv")
saveRDS(ParamInit,file=paste0("~/MITUS/inst/US/US_ParamInit_", Sys.Date(),".rds"))
}
