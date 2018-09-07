#'load initial conditions
#'load necessary libraries
library(mnormt)
library(parallel)
library(lhs)
library(Rcpp)
library(MCMCpack)
library(MASS)

#' Calibration Target Environment
# targets <- new.env()

#'all other loads will be dependent on geography
model_load<-function(loc="US"){
#'lazy load necessary datasets
#'Model Input
model_inputs<-paste0(loc,"_ModelInputs_9-6-18")
par_init<-paste0(loc,"_ParamInit_09-7-2018")
start_val<-paste0(loc,"_StartVal_2018-09-07")
calib_dat<-paste0(loc,"_CalibDat_2018-09-01")

data(list=model_inputs, package = 'MITUS')
data(list=par_init, package = 'MITUS')
data(list=start_val, package = 'MITUS')
data(list=calib_dat, package='MITUS')# ParamInit


#'creation of background parameters
#'elements of P will be replaced from either the StartVals in the case
#'of optimization or the user inputted dataset
P  <<- ParamInit[,1];
names(P) <<- rownames(ParamInit)

ii <<-  ParamInit[,5]==1
ParamInitZ <<- ParamInit[ParamInit$Calib==1,]
idZ0 <<- ParamInitZ[,4]==0
idZ1 <<- ParamInitZ[,4]==1
idZ2 <<- ParamInitZ[,4]==2

#'format the calibration data
wts <<- CalibDat[["ImptWeights"]]
# targets[["ParamInitZ"]]<-ParamInitZ
# targets[["idZ0"]]<-idZ0
# targets[["idZ1"]]<-idZ1
# targets[["idZ2"]]<-idZ2
# targets[["P"]]<-P

return(invisible(NULL))
}
