# locs<-c("CA","FL","GA","IL","NJ" ,"NY", "PA","TX", "VA","WA")
# # # #  make_results<-function(locs){
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#   CalibDat<-CalibDatState<-readRDS("~/Desktop/GoodStateOptims/CalibDat_2018-11-14.rds")
#   ParamInit_st<-ParamInit<-readRDS("~/Desktop/GoodStateOptims/ParamInit_2018-11-14.rds")
#   StartVal_st<-StartVal<-readRDS("~/Desktop/GoodStateOptims/StartVal_2018-11-14.rds")
#   Inputs<-readRDS(paste0("~/Desktop/GoodStateOptims/",loc,"/",loc,"_ModelInputs_01-24-19.rds"))
#
#   Opt <- readRDS(system.file(paste0(loc,"/",loc,"_opt_all_2018-12-12.rds"), package="MITUS"))
#   Par <- readRDS(paste0("~/Desktop/GoodStateOptims/parAll_",loc,"_10_2018_11_27.rds"))
#
#   OutputsInt(Par,loc,n_cores=1,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0)
#
#   }
#  # }
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
# model_load(loc)
# make_results(loc)
#   }
#
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#   reshaping(loc)
#   }
