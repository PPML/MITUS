#'This function is used to run a optimization on the input parameters.
#'The only input to the function is a "batch number" which will determine
#'which row of the starting values dataframe the optimization will use.
#'It will return 8 datasets of the optimized parameters --one from each
#'optimization step. The final (8th) optimization step is an univariate
#'optimization.
#'The function requires the use of the ParamInit Rdata file and the
#'StartValues Rdata file.
#'@name optim_b
#'@param b batch number; must be > 21; corresponds to a row of start vals
#'@return 8 datasets from optimization loop
#'@export

optim_b <- function(b){

load("~/MITUS/data/ParamInit_2018.rData")
P  <- ParamInit[,1]
names(P) <- rownames(ParamInit)
ii <-  ParamInit[,5]==1
ParamInitZ <- ParamInit[ParamInit$Calib==1,]
idZ0 <- ParamInitZ[,4]==0
idZ1 <- ParamInitZ[,4]==1
idZ2 <- ParamInitZ[,4]==2
load("data/StartVal_ 2018-06-28 .rData") # StartVal

posterior = function(theta) { -lprior(theta) - llikelihood(theta,n_cores) }


# for (i in min(b, nrow(StartVal))){
  o1  <- optim(StartVal[b,], posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
  save(o1,file=paste("Opt_US_r1_",b,"_new.rData",sep=""))
  o2  <- optim(o1$par,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
  save(o2,file=paste("Opt_US_r2_",b,"_new.rData",sep=""))
  o3  <- optim(o2$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
  save(o3,file=paste("Opt_US_r3_",b,"_new.rData",sep=""))
  o4  <- optim(o3$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
  save(o4,file=paste("Opt_US_r4_",b,"_new.rData",sep=""))
  o5  <- optim(o4$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
  save(o5,file=paste("Opt_US_r5_",b,"_new.rData",sep=""))
  o6  <- optim(o5$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
  save(o6,file=paste("Opt_US_r6_",b,"_new.rData",sep=""))
  o7  <- optim(o6$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
  save(o7,file=paste("Opt_US_r7_",b,"_new.rData",sep=""))
  o8 <- UnivOptim(o7$par) ; o8$value
  save(o8,file=paste("Opt_US_r8_",b,"_new.rData",sep=""))
# }
}
