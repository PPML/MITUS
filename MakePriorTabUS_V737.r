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
ParamInit <- as.data.frame(read.csv("ParamInitUS_V737.csv")[,2:6]); rownames(ParamInit) <- read.csv("ParamInitUS_V737.csv")[,1]

ParamInit <- cbind(ParamInit,NA,NA,NA)
colnames(ParamInit) <- c(colnames(ParamInit)[1:5],"Par1","Par2","TransMean")

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

### Calculate transformed
  ParamInit[ParamInit[,4]==0,8] <- logit(ParamInit[ParamInit[,4]==0,1])
  ParamInit[ParamInit[,4]==1,8] <- log(  ParamInit[ParamInit[,4]==1,1])
  ParamInit[ParamInit[,4]==2,8] <-       ParamInit[ParamInit[,4]==2,1] 
  ParamInit[ParamInit[,8]==-Inf,8] <- -100;  ParamInit[ParamInit[,8]==Inf,8] <- 100

### Create VCV for MV unit normal
  VCVall <- diag(sum(ParamInit$Calib))

write.csv(ParamInit,file="ParamInitUS_V737_final.csv")
save(ParamInit,file="ParamInitUS_V737tab.rData")
save(VCVall,file="VCVall_V737.rData")

##################################################
####### DRAWS SOME STARTING VALUES

  sample.prior2 <- function(n)  { # n=5
    jj <-ParamInit[ParamInit$Calib==1,]
    pari1 <- rmvnorm(n,rep(0,nrow(jj)),VCVall)
    pari1 }

  set.seed(351); StartVal <- sample.prior2(n=20) 
  save(StartVal,file="StartValUS_9-5-2016.rData")

#########################################
## Parameter stuff
  ParamInitZ <- ParamInit[ParamInit$Calib==1,]
  idZ0 <- ParamInitZ[,4]==0
  idZ1 <- ParamInitZ[,4]==1
  idZ2 <- ParamInitZ[,4]==2

  library(mnormt)
  lPrior2 <- function(Par,Par3) {
    if(dim(as.matrix(Par))[2]==1) Par <- t(as.matrix(Par))
    ldensity  <- dmnorm(Par,rep(0,nrow(ParamInitZ)),diag(nrow(ParamInitZ)),log=T)
    ldensity2 <- ldensity-sum(dnorm(Par,0,1,log=T))
    lDensTrue <- rep(NA,nrow(ParamInitZ))
    lDensTrue[idZ0] <- dgamma(Par3[idZ0], shape   = ParamInitZ[idZ0,6], rate   = ParamInitZ[idZ0,7],log=T)
    lDensTrue[idZ1] <- dbeta( Par3[idZ1], shape1  = ParamInitZ[idZ1,6], shape2 = ParamInitZ[idZ1,7],log=T)
    lDensTrue[idZ2] <- dnorm( Par3[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7],log=T)  
    ldensity3 <- ldensity2+sum(lDensTrue)
    ldensity3  }

