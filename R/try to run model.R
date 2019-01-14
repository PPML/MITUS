#' library(mnormt)
#' library(mvtnorm)
#' library(parallel)
#' library(lhs)
#' library(Rcpp)
#' library(MCMCpack)
#' library(MASS)
#'
#' load("data/CalibDatState_7-2-18.rda")
#' load("data/MA_ModelInputs_11-13-18.rda")
#' load("data/ST_StartVal_2018-11-14.rda")
#' load("data/ParamInit_st_2018-11-14.rda")
#' load("data/wonder_pop.rda")
#' load("data/ST_tot_mort.rda")
#' load("data/mort_ag_16.rda")
#' load("data/stateID.rda")
#'
#' wts <<- CalibDatState[["ImptWeights"]]
#' wtZ <- wts[44:67];  wtZ["2016"] <- 4
#' wtZ<<-wtZ
#'
#' #'creation of background parameters
#' #'elements of P will be replaced from either the StartVals in the case
#' #'of optimization or the user inputted dataset
#'
#' P  <<- ParamInit_st[,1]
#' names(P) <<- rownames(ParamInit_st)
#'
#' ii <<-  ParamInit_st[,5]==1
#' ParamInitZ <<- ParamInit_st[ParamInit_st$Calib==1,]
#' idZ0 <<- ParamInitZ[,4]==0
#' idZ1 <<- ParamInitZ[,4]==1
#' idZ2 <<- ParamInitZ[,4]==2
#' ParamInit<-ParamInit_st
#'
#' Inputs<-MA_Inputs
#' rm(MA_Inputs)
#'
#' ParamInit<-ParamInit_st
#' rm(ParamInit_st)
#'
#' source("R/calib_functions_st.r")
#' source("R/IMIS_functions_st.r")
#'
#' source("R/basic_functions.R")
#' source("R/calib_functions.R")
#' source("R/define_dist.R")
#' source("R/func_prior.R")
#' source("R/func_optim_demo.R")
#' source("R/param_init.R")
#' source("R/param.R")
#'
#' sourceCpp("src/reblncd.cpp")
#' sourceCpp("src/tb_model.cpp")
