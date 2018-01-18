########## Optimize  for each country x5

setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC Large Grant/Analysis Transmission")

ScriptMak <- function(b) {
  
  ################################################
  Script <- paste('
  
  b    <-  ',b,'   # batch number

#########################  SET-UP  #########################  

## Parameter stuff
  load("ParamInitUS_V738tab.rData") # ParamInit
  P  <- ParamInit[,1]; names(P) <- rownames(ParamInit)
  ii <-  ParamInit[,5]==1
  ParamInitZ <- ParamInit[ParamInit$Calib==1,]
  idZ0 <- ParamInitZ[,4]==0
  idZ1 <- ParamInitZ[,4]==1
  idZ2 <- ParamInitZ[,4]==2

## Scripts and functions
load("ModelInputs_9-2-16.rData")
  source("ParamUS_V504.r")
  source("PriorFuncUS_V22.r")
  source("CalibFunctionsUS_V22.r")
  source("TimeStepUS_V891.r")
  source("IMISfunctionsUS_V236.r")
  source("OptUniv_V3.r")
  load("StartValUS_9-6-2016.rData") # StartVal

  posterior = function(theta) { -lprior(theta) - llikelihood(theta,n_cores) }

##### OPTIMISE ###########
    o1  <- optim(StartVal[b,], posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
      save(o1,file=paste("Opt_US540_r1_",b,"_1-27-16.rData",sep=""))
    o2  <- optim(o1$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
      save(o2,file=paste("Opt_US540_r2_",b,"_1-27-16.rData",sep=""))
    o3  <- optim(o2$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
      save(o3,file=paste("Opt_US540_r3_",b,"_1-27-16.rData",sep=""))
    o4  <- optim(o3$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
      save(o4,file=paste("Opt_US540_r4_",b,"_1-27-16.rData",sep=""))
    o5  <- optim(o4$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
      save(o5,file=paste("Opt_US540_r5_",b,"_1-27-16.rData",sep=""))
    o6  <- optim(o5$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
      save(o6,file=paste("Opt_US540_r6_",b,"_1-27-16.rData",sep=""))
    o7  <- optim(o6$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
      save(o7,file=paste("Opt_US540_r7_",b,"_1-27-16.rData",sep=""))
    o8 <- UnivOptim(o7$par) ; o8$value
     save(o8,file=paste("Opt_US540_r8_",b,"_1-27-16.rData",sep=""))

############ Done!   #####################################

  quit("no")  # Kill
  
                    ',sep=""); return(Script)  }

####################################

for (b in 1:10) {  kkk <- ScriptMak(b); write(kkk,file=paste("OptUS540_",b,".r",sep=""))  } 




