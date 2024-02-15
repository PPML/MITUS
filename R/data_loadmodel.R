#'loads data for the geography of interest
#'@name model_load
#'@param loc two letter postal abbreviation for states; US for national
#'@return void
model_load<-function(loc="US"){
#' add loc as a global variable
# loc<<-loc
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)
#'load necessary datasets
#'Model Input
if (loc=="US"){
  CalibDat<<-readRDS(system.file("US/US_CalibDat_2022-04-13.rds", package="MITUS"))
  ## ParamInit and StartVal were replaced on 8/1 to clean up, but 07-07 are last calibrated version
  ParamInit<<-as.data.frame(readRDS(system.file("US/US_ParamInit_2022-08-01.rds", package="MITUS")))
  StartVal<<-readRDS(system.file("US/US_StartVal_2022-08-01.rds", package="MITUS"))
  Inputs<<-readRDS(system.file("US/US_Inputs_08-31-20.rds", package="MITUS"))
  Opt <<- readRDS(system.file("US/US_Optim_all_25_0119.rds", package="MITUS"))
  Par <<- readRDS(system.file("US/US_Param_all_25_0119.rds", package="MITUS"))
} else {
  CalibDat<<-CalibDatState<<-readRDS(system.file("ST/ST_CalibDat_04-20-22.rds", package="MITUS"))
  ParamInit_st<<-ParamInit<<-readRDS(system.file("ST/ST_ParamInit_2022-07-08.rds", package="MITUS"))
  StartVal_st<<-StartVal<<-readRDS(system.file("ST/ST_StartVal_2022-07-08.rds", package="MITUS"))
  Inputs<<-readRDS(system.file(paste0(loc,"/",loc,"_ModelInputs_11-12-21.rds"), package="MITUS"))
  #last input change was to update the RR active TB by age in immigrants
}
if (loc=="US"){
  LgtCurveY2 <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
  zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2019-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
  zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }
  ImptWeights <- LgtCurveY2(2000,2019,0.95)+0.05
  names(ImptWeights) <- 1950:2019

  wts <<- ImptWeights
  P  <<- ParamInit[,1]
  names(P) <<- rownames(ParamInit)

  ii <<-  ParamInit[,5]==1
  ParamInitZ <<- ParamInit[ParamInit$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
} else {
  LgtCurveY2 <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
  zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2019-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
  zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }

  ImptWeights <- LgtCurveY2(2000,2019,0.95)+0.05
  names(ImptWeights) <- 1950:2019
  wts <<- ImptWeights
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
  ParamInit<<-ParamInit_st
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
