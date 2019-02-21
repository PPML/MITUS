#' Load data for the geography of interest
#' 
#' Load necessary datasets for the model; Export them to the global scope by default, 
#' or return them in an environment. 
#'
#'@name model_load_env
#'@param loc two letter postal abbreviation for states; US for national
#'@return void
model_load<-function(loc="US", return_env=F){
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)
  
  # Construct an environment where we will store the calibration data
  output_env <- new.env()

	# Attach to the output_env and load the necessary datasets 
	with(output_env, {
  if (loc=="US"){
		CalibDat<-readRDS(system.file("US/US_CalibDat_02-14-19.rds", package="MITUS"))
		ParamInit<-readRDS(system.file("US/US_ParamInit_02-14-19.rds", package="MITUS"))
		Inputs<-readRDS(system.file("US/US_Inputs_01-24-19.rds", package="MITUS"))
		StartVal<-readRDS(system.file("US/US_StartVal_02-14-19.rds", package="MITUS"))

    wts <- CalibDat[["ImptWeights"]]
    P  <- ParamInit[,1]
    names(P) <- rownames(ParamInit)
    
    ii <-  ParamInit[,5]==1
    ParamInitZ <- ParamInit[ParamInit$Calib==1,]
    idZ0 <- ParamInitZ[,4]==0
    idZ1 <- ParamInitZ[,4]==1
    idZ2 <- ParamInitZ[,4]==2

  } else {

		CalibDat<-CalibDatState<-readRDS(system.file("ST/ST_CalibDat_01-24-19.rds", package="MITUS"))
		ParamInit_st<-ParamInit<-readRDS(system.file("ST/ST_ParamInit_02-09-19.rds", package="MITUS"))
		StartVal_st<-StartVal<-readRDS(system.file("ST/ST_StartVal_02-09-19.rds", package="MITUS"))
		Inputs<-readRDS(system.file(paste0(loc,"/",loc,"_ModelInputs_01-24-19.rds"), package="MITUS"))

    wts <- CalibDatState[["ImptWeights"]]
    W <- wts[44:67];  W["2016"] <- 4
    wtZ <-W
    
    #creation of background parameters
    #elements of P will be replaced from either the StartVals in the case
    #of optimization or the user inputted dataset
    
    P  <- ParamInit_st[,1]
    names(P) <- rownames(ParamInit_st)
    ii <-  ParamInit_st[,5]==1
    ParamInitZ <- ParamInit_st[ParamInit_st$Calib==1,]
    idZ0 <- ParamInitZ[,4]==0
    idZ1 <- ParamInitZ[,4]==1
    idZ2 <- ParamInitZ[,4]==2
    ParamInit<-ParamInit_st
		}
	})
  
	# If the return_env argument is true, we return an environment containing the loaded datasets.
	# Otherwise we export each element in the output_env to the global environment.
	if (return_env) {
		return(output_env)
	} else {
		for (object_name in ls(output_env)) {
			assign(x = object_name, value = get(object_name, envir = output_env), envir = .GlobalEnv)
		}
		return(invisible(NULL))
	}
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
  if (loc=="US")
  {
    calib_dat<-paste0(loc,"_CalibDat_2018-09-17")

  } else {
    calib_dat<-paste0("CalibDatState_7-2-18")

  }
  model_inputs<-paste0(loc,"_ModelInputs_11-13-18")
  par_init<-paste0("US_ParamInit_demo")
  start_val<-paste0("US_StartVal_demo")

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

  return(invisible(NULL))
}
