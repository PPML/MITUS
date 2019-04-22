#'loads data for the geography of interest
#'@name load_tabby
#'@param loc two letter postal abbreviation for states; US for national
#'@return tabby_env
load_tabby<-function(loc="US"){
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)
  #'create a new environment
  #be careful about defaults and parent environments
  tabby_env<-new.env()
  #'load necessary datasets into the environment
  with(tabby_env,{
    if (loc=="US"){
      CalibDat<<-readRDS(system.file("US/US_CalibDat_2018-09-17.rds", package="MITUS"))
      ParamInit<<-readRDS(system.file("US/US_ParamInit_2018-12-05.rds", package="MITUS"))
      StartVal<<-readRDS(system.file("US/US_StartVal_2018-12-05.rds", package="MITUS"))
      Inputs<<-readRDS(system.file("US/US_Inputs_12-05-2018.rds", package="MITUS"))

    } else {
      CalibDat<<-CalibDatState<<-readRDS(system.file("ST/CalibDatState_7-2-18.rds", package="MITUS"))
      ParamInit_st<<-ParamInit<<-readRDS(system.file("ST/ST_ParamInit_2018-12-10.rds", package="MITUS"))
      StartVal_st<<-StartVal<<-readRDS(system.file("ST/ST_StartVal_2018-12-10.rds", package="MITUS"))
      Inputs<<-readRDS(system.file(paste0(loc,"/",loc,"_ModelInputs_01-24-19.rds"), package="MITUS"))
      Opt <<- readRDS(system.file(paste0(loc,"/",loc,"_opt_all_2018-12-12.rds"), package="MITUS"))
      Par <<- readRDS(system.file(paste0(loc,"/parAll10_2018-12-12.rds"), package="MITUS"))
    }
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    # use getmode function to determine the mode of the optim dataset and then select one
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    # j<-0
    # Par<-vector()
    # for (i in 1:nrow(Opt)){
    #   while(j<1){
    #     if(Opt[i,ncol(Opt)]==getmode(Opt[i,ncol(Opt)])){
    #       Par<-Opt[i,1:(ncol(Opt)-1)]
    #       j<-j+1
    #     }}}
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    # make the calibration plots for the selected location
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    # tabby_calib_graphs(loc=loc,out_i = 1) #out_i needs to be updated; perhaps based off the reshaped data
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    # update the Parameter Vector for the optimized values
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    if (loc=="US"){
      P  <- ParamInit[,1]
      names(P) <- rownames(ParamInit)
      ii <-  ParamInit[,5]==1
      ParamInitZ <- ParamInit[ParamInit$Calib==1,]
      idZ0 <- ParamInitZ[,4]==0
      idZ1 <- ParamInitZ[,4]==1
      idZ2 <- ParamInitZ[,4]==2
      #add in the optimized values to P
      # Par2 <- pnorm(Par,0,1)
      # # uniform to true
      # Par3 <- Par2
      # Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
      # Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
      # Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
      # P[ii] <- Par3
      # P <- P

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
      # Par2 <- pnorm(Par,0,1)
      # # uniform to true
      # Par3 <- Par2
      # Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
      # Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
      # Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
      # P[ii] <- Par3
      # P <- P
    }
  })
  return(tabby_env)
}
