#'loads data for the geography of interest
#'@name model_load
#'@param loc two letter postal abbreviation for states; US for national
#'@return void
model_load<-function(loc="US"){
#' add loc as a global variable
loc<<-loc
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)
#'load necessary datasets
#'Model Input
if (loc=="US"){
  CalibDat<<-readRDS(system.file("US/US_CalibDat_03-06-19.rds", package="MITUS"))
  ParamInit<<-readRDS(system.file("US/US_ParamInit_07-11-19.rds", package="MITUS"))
  StartVal<<-readRDS(system.file("US/US_StartVal_07-11-19.rds", package="MITUS"))
  Inputs<<-readRDS(system.file("US/US_Inputs_06-26-19.rds", package="MITUS"))

} else {
  CalibDat<<-CalibDatState<<-readRDS(system.file("ST/ST_CalibDat_07-15-19.rds", package="MITUS"))
  ParamInit_st<<-ParamInit<<-readRDS(system.file("ST/ST_ParamInit_2019-07-25.rds", package="MITUS"))
  StartVal_st<<-StartVal<<-readRDS(system.file("ST/ST_StartVal_2019-07-25.rds", package="MITUS"))
  Inputs<<-readRDS(system.file(paste0(loc,"/",loc,"_ModelInputs_07-03-19.rds"), package="MITUS"))
}

if (loc=="US"){
  wts <<- CalibDat[["ImptWeights"]]
  P  <<- ParamInit[,1]
  names(P) <<- rownames(ParamInit)

  ii <<-  ParamInit[,5]==1
  ParamInitZ <<- ParamInit[ParamInit$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
} else {
  wts <<- CalibDatState[["ImptWeights"]]
  W <- wts[44:69];  W["2016"] <- 4
  wtZ <<-W

  #creation of background parameters
  #elements of P will be replaced from either the StartVals in the case
  #of optimization or the user inputted dataset

  P  <<- ParamInit_st[,1]
  names(P) <<- rownames(ParamInit_st)
  ii <<-  ParamInit_st[,5]==1
  ParamInitZ <<- ParamInit_st[ParamInit_st$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
  ParamInit<-ParamInit_st
}
  prgchng<<-def_prgchng(P)

return(invisible(NULL))
}

#'loads data for the geography of interest
#'@name model_load_demo
#'@param loc two letter postal abbreviation for states; US for national
#'@return void
model_load_demo<-function(loc="US"){
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)
  #'lazy load necessary datasets
  #'Model Input
  if (loc=="US"){
    CalibDat<<-readRDS(system.file("US/US_CalibDat_03-06-19.rds", package="MITUS"))

    # CalibDat<<-readRDS(system.file("US/US_CalibDat_02-14-19.rds", package="MITUS"))
    ParamInit_demo<<-readRDS(system.file("US/US_ParamInitdemo_2019-02-21.rds", package="MITUS"))
    StartVal_demo<<-readRDS(system.file("US/US_StartValdemo_2019-02-21.rds", package="MITUS"))
    Inputs<<-readRDS(system.file("US/US_Inputs_01-24-19.rds", package="MITUS"))

  } else {
    CalibDat<<-CalibDatState<<-readRDS(system.file("ST/ST_CalibDat_01-24-19.rds", package="MITUS"))
    ParamInit_st<<-ParamInit<<-readRDS(system.file("ST/ST_ParamInit_02-09-19.rds", package="MITUS"))
    StartVal_st<<-StartVal<<-readRDS(system.file("ST/ST_StartVal_02-09-19.rds", package="MITUS"))
    Inputs<<-readRDS(system.file(paste0(loc,"/",loc,"_ModelInputs_01-24-19.rds"), package="MITUS"))
  }

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
  P  <<- ParamInit_demo[,1]
  names(P) <<- rownames(ParamInit_demo)

  ii <<-  ParamInit_demo[,5]==1
  ParamInitZ <<- ParamInit_demo[ParamInit_demo$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2

  #'format the calibration data
  # targets[["ParamInitZ"]]<-ParamInitZ
  # targets[["idZ0"]]<-idZ0
  # targets[["idZ1"]]<-idZ1
  # targets[["idZ2"]]<-idZ2
  # targets[["P"]]<-P
  ## ALSO LOAD IN THE BASELINE PROGRAM CHANGE VECTOR
  prgchng<<-def_prgchng(P)
  return(invisible(NULL))
}
