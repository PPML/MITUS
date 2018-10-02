#'This function is used to run a optimization on the input parameters.
#'The only input to the function is a "batch number" which will determine
#'which row of the starting values dataframe the optimization will use.
#'It will return 8 datasets of the optimized parameters --one from each
#'optimization step. The final (8th) optimization step is an univariate
#'optimization.
#'The function requires the use of the ParamInit Rdata file and the
#'StartValues Rdata file.
#'@name optim_b
#'@param df dataframe or matrix of starting values data frame
#'@param samp_i which rows of the data frame to use
#'@return 8 datasets from optimization loop
#'@export

optim_b <- function(df, samp_i=1){
# data("StartVal_2018-08-06", package = "MITUS")


posterior = function(theta) {
   -lprior(theta) - llikelihood(theta,n_cores)
 }

if(min(dim(as.data.frame(df)))==1) {
  df1 <- as.numeric(df)
  names(df1) <- names(df)
} else{
  df1 <- as.numeric(df[samp_i,])
}

b<-samp_i
# for (i in min(b, nrow(StartVal))){
  o1  <- optim(df1, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
  save(o1,file=paste("Opt_US_r1_", b,"_", Sys.Date(),".rda",sep=""))
  o2  <- optim(o1$par,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
  save(o2,file=paste("Opt_US_r2_", b,"_", Sys.Date(),".rda",sep=""))
  o3  <- optim(o2$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
  save(o3,file=paste("Opt_US_r3_", b,"_", Sys.Date(),".rda",sep=""))
  o4  <- optim(o3$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
  save(o4,file=paste("Opt_US_r4_", b,"_", Sys.Date(),".rda",sep=""))
  o5  <- optim(o4$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
  save(o5,file=paste("Opt_US_r5_", b,"_", Sys.Date(),".rda",sep=""))
  o6  <- optim(o5$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
  save(o6,file=paste("Opt_US_r6_", b,"_", Sys.Date(),".rda",sep=""))
  o7  <- optim(o6$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
  save(o7,file=paste("Opt_US_r7_", b,"_", Sys.Date(),".rda",sep=""))
  # o8 <- UnivOptim(o7$par) ; o8$value
  # save(o8,file=paste("Opt_US_r8_", b,"_", Sys.Date(),".rda",sep=""))
# }
}

#'This function is used to run a optimization on the input parameters.
#'The only input to the function is a "batch number" which will determine
#'which row of the starting values dataframe the optimization will use.
#'It will return 8 datasets of the optimized parameters --one from each
#'optimization step. The final (8th) optimization step is an univariate
#'optimization.
#'The function requires the use of the ParamInit Rdata file and the
#'StartValues Rdata file.
#'@name optim_b_st
#'@param df dataframe or matrix of starting values data frame
#'@param samp_i which rows of the data frame to use
#'@return 8 datasets from optimization loop
#'@export

optim_b_st <- function(df, samp_i=1, st){
  # data("StartVal_2018-08-06", package = "MITUS")
data("stateID", package="MITUS")
  loc<-stateID[st,3]

  posterior_st = function(theta,st,n_cores) {
    -lprior(theta) - llikelihood_st(theta,st,n_cores)
  }

  if(min(dim(as.data.frame(df)))==1) {
    df1 <- as.numeric(df)
    names(df1) <- names(df)
  } else{
    df1 <- as.numeric(df[samp_i,])
  }

  b<-samp_i
  # for (i in min(b, nrow(StartVal))){
  o1  <- optim(df1, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
  save(o1,file=paste("Opt_",loc,"_r1_", b,"_", Sys.Date(),".rda",sep=""))
  o2  <- optim(o1$par,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
  save(o2,file=paste("Opt_",loc,"_r2_", b,"_", Sys.Date(),".rda",sep=""))
  o3  <- optim(o2$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
  save(o3,file=paste("Opt_",loc,"_r3_", b,"_", Sys.Date(),".rda",sep=""))
  o4  <- optim(o3$par ,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
  save(o4,file=paste("Opt_",loc,"_r4_", b,"_", Sys.Date(),".rda",sep=""))
  o5  <- optim(o4$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
  save(o5,file=paste("Opt_",loc,"_r5_", b,"_", Sys.Date(),".rda",sep=""))
  o6  <- optim(o5$par ,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
  save(o6,file=paste("Opt_",loc,"_r6_", b,"_", Sys.Date(),".rda",sep=""))
  o7  <- optim(o6$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
  save(o7,file=paste("Opt_",loc,"_r7_", b,"_", Sys.Date(),".rda",sep=""))
  # o8 <- UnivOptim(o7$par) ; o8$value
  # save(o8,file=paste("Opt_ST_r8_", b,"_", Sys.Date(),".rda",sep=""))
  # }
}
