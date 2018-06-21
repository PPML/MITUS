#'Function to create scripts for optimization
#'@name ScriptMak
#'@param b number of scripts to generate
#'@return b number of scripts for optimization
#'@export
ScriptMak <- function(b) {
for (i in 1:b){
Script <- paste('

  b    <-  ',i,'   # batch number

#setup the parameters for the model
load("data/ParamInitUS_V738tab.rData") # ParamInit
P  <- ParamInit[,1]
names(P) <- rownames(ParamInit)
#logical vector for whether the parameter is to be calibrated
ii <-  ParamInit[,5]==1
#create a list of the parameters (and values across the dataframe) to be calibrated
ParamInitZ <- ParamInit[ParamInit$Calib==1,]
#logical vectors dependent on the distribution of the parameters
idZ0 <- ParamInitZ[,4]==0
idZ1 <- ParamInitZ[,4]==1
idZ2 <- ParamInitZ[,4]==2

## load data
load("data/ModelInputs_9-2-16.rData")
load("data/StartValUS_9-5-2016.rData") # StartVal
 posterior = function(theta) { -lprior(theta) - llikelihood(theta,n_cores) }

##### OPTIMISE ###########
o1  <- optim(StartVal[b,], posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
save(o1,file=paste("Opt_US540_r1_",b,"_", Sys.Date(), ".rData",sep=""))
o2  <- optim(o1$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
save(o2,file=paste("Opt_US540_r2_",b,"_", Sys.Date(), ".rData",sep=""))
o3  <- optim(o2$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
save(o3,file=paste("Opt_US540_r3_",b,"_", Sys.Date(), ".rData",sep=""))
o4  <- optim(o3$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
save(o4,file=paste("Opt_US540_r4_",b,"_", Sys.Date(), ".rData",sep=""))
o5  <- optim(o4$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
save(o5,file=paste("Opt_US540_r5_",b,"_", Sys.Date(), ".rData",sep=""))
o6  <- optim(o5$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
save(o6,file=paste("Opt_US540_r6_",b,"_", Sys.Date(), ".rData",sep=""))
o7  <- optim(o6$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
save(o7,file=paste("Opt_US540_r7_",b,"_", Sys.Date(), ".rData",sep=""))
o8 <- UnivOptim(o7$par) ; o8$value
save(o8,file=paste("Opt_US540_r8_",b,"_", Sys.Date(), ".rData",sep=""))

############ Done!   #####################################

quit("no")  # Kill

',sep="")

write(Script,file=paste("OptUS540_",i,".r",sep=""))

} }









