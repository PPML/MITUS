# # locs<-c("CA","FL","GA","IL","MA","NJ" ,"NY", "PA","TX", "VA","WA")
# # # #  make_results<-function(locs){
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#   CalibDat<-CalibDatState<-readRDS("~/Desktop/GoodStateOptims/CalibDat_2018-11-14.rds")
#   ParamInit_st<-ParamInit<-readRDS("~/Desktop/GoodStateOptims/ParamInit_2018-11-14.rds")
#   StartVal_st<-StartVal<-readRDS("~/Desktop/GoodStateOptims/StartVal_2018-11-14.rds")
#   Inputs<-readRDS(paste0("~/Desktop/GoodStateOptims/",loc,"/",loc,"_ModelInputs_01-24-19.rds"))
#
#    Opt <- readRDS(system.file(paste0(loc,"/",loc,"_opt_all_2018-12-12.rds"), package="MITUS"))
#   Par <- readRDS(paste0("~/Desktop/GoodStateOptims/parAll_",loc,"_10_2018_11_27.rds"))
#
#
  # for (intv in 0:8) {
  #   # intvs will be a vector of 0s except for (possibly) one activated intervention
  #   intvs <- rep(0, 8)
  #
  #   # If intv is 0, disable all interventions to run the basecase.
  #   # Otherwise, activate one of the 8 interventions.
  #   if (intv != 0) intvs[intv] <- 1
  #   #use the basecase program change
  #   defprg<-def_prgchng(Par[10,])
  #   # Simulate using the elements of intvs to control
  #   # whether or not each intervention is on.
  #   new_OutputsInt(loc,Par,n_cores=1,endyr=2050,
  #                     Int1 = intvs[[1]],
  #                     Int2 = intvs[[2]],
  #                     Int3 = intvs[[3]],
  #                     Int4 = intvs[[4]],
  #                     Int5 = intvs[[5]],
  #                     Scen1 = intvs[[6]],
  #                     Scen2 = intvs[[7]],
  #                     Scen3 = intvs[[8]],
  #                     prg_chng = defprg)
  #
  #   # Do something to save the simulation outcomes
  # }
# # format_as_restab(loc)
#   }
# }
# #s
#  # for (i in 1:length(locs)){
#  #   loc<-locs[i]
#  #   print(loc)
#  # tabby_calib_graphs(loc)
#  #   }
#
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#
#   model_load(loc)
#   format_as_restab(loc)
# }
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#
# # model_load(loc)
# make_results(loc)
#   }
# #
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#  #onesim(loc)
# make_results(loc)
#   }
