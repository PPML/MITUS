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
  rm(Inputs)
  rm(CalibDat)
#'lazy load necessary datasets
#'Model Input
if (loc=="US")
{
  calib_dat<-paste0(loc,"_CalibDat_2018-09-17")

} else {
  calib_dat<-paste0("CalibDatState_7-2-18")

}
model_inputs<-paste0(loc,"_ModelInputs_9-6-18")
par_init<-paste0("US_ParamInit_2018-10-15")
start_val<-paste0("US_StartVal_2018-10-15")

data(list=model_inputs, package = 'MITUS')
data(list=par_init, package = 'MITUS')
data(list=start_val, package = 'MITUS')
data(list=calib_dat, package='MITUS')# ParamInit

if (loc=="US")
{
  wts <<- CalibDat[["ImptWeights"]]

} else {
  wts <<- CalibDatState[["ImptWeights"]]
  wtZ <- wts[44:67];  wtZ["2016"] <- 4
  wtZ<<-wtZ

}
#'creation of background parameters
#'elements of P will be replaced from either the StartVals in the case
#'of optimization or the user inputted dataset
#'
P  <<- ParamInit[,1];
names(P) <<- rownames(ParamInit)

ii <<-  ParamInit[,5]==1
ParamInitZ <<- ParamInit[ParamInit$Calib==1,]
idZ0 <<- ParamInitZ[,4]==0
idZ1 <<- ParamInitZ[,4]==1
idZ2 <<- ParamInitZ[,4]==2

#'format the calibration data
# targets[["ParamInitZ"]]<-ParamInitZ
# targets[["idZ0"]]<-idZ0
# targets[["idZ1"]]<-idZ1
# targets[["idZ2"]]<-idZ2
# targets[["P"]]<-P

return(invisible(NULL))
}
