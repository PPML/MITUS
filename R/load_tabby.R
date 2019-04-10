#'loads data for the geography of interest
#'@name load_tabby
#'@param loc two letter postal abbreviation for states; US for national
#'@return env
load_tabby<-function(loc="US"){
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)
  #'create a new environment
  tabby_env<-new.env()
  #'load necessary datasets into the environment
  with(tabby_env,{
  if (loc=="US"){
    CalibDat<-readRDS(system.file("US/US_CalibDat_03-06-19.rds", package="MITUS"))
    ParamInit<-readRDS(system.file("US/US_ParamInit_01-24-19.rds", package="MITUS"))
    Inputs<-readRDS(system.file("US/US_Inputs_01-24-19.rds", package="MITUS"))
    # OptParam<-readRDS(system.file("US/US_Optim_all_01-24-19.rds"), package="MITUS")
    #add reshaped data
    #add calibplots
  } else {
    CalibDat<-CalibDatState<-readRDS(system.file("ST/ST_CalibDat_01-24-19.rds", package="MITUS"))
    ParamInit_st<-ParamInit<-readRDS(system.file("ST/ST_ParamInit_02-09-19.rds", package="MITUS"))
    Inputs<-readRDS(system.file(paste0(loc,"/",loc,"_ModelInputs_01-24-19.rds"), package="MITUS"))
    # OptParam<-readRDS(system.file(paste0(loc,"/",loc,"_Optim_all_01-24-19.rds"), package="MITUS"))
    #add reshaped data
    #add calibplots
  }

  if (loc=="US"){
    P  <- ParamInit[,1] #definition of P will change
    names(P) <- rownames(ParamInit)

    ii <-  ParamInit[,5]==1
    ParamInitZ <- ParamInit[ParamInit$Calib==1,]
    idZ0 <- ParamInitZ[,4]==0
    idZ1 <- ParamInitZ[,4]==1
    idZ2 <- ParamInitZ[,4]==2
  } else {

    #creation of background parameters
    #elements of P will be replaced from either the StartVals in the case
    #of optimization or the user inputted dataset

    P  <- ParamInit_st[,1] #definition of P will change
    names(P) <- rownames(ParamInit_st)
    ii <-  ParamInit_st[,5]==1
    ParamInitZ <- ParamInit_st[ParamInit_st$Calib==1,]
    idZ0 <- ParamInitZ[,4]==0
    idZ1 <- ParamInitZ[,4]==1
    idZ2 <- ParamInitZ[,4]==2
    ParamInit<-ParamInit_st
  }
})
  return(tabby_env)
}
