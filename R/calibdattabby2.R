# ##format calibration data for tabby2 plots
# ##united states
# NCD<-readRDS(system.file("US/US_CalibDat_03-06-19.rds", package="MITUS"))
#
# #states
# SCD<-readRDS(system.file("ST/ST_CalibDat_07-03-19.rds", package="MITUS"))
#
# for(i in 1:length(NCD)){
#   x<-NCD[[i]]
#   saveRDS(x,file = paste0("~/MITUS/inst/US/calibration_targets/US_",names(NCD)[i],"_03-06-19.rds"))
# }
#
# for (i in 1:length(SCD)){
#   x<-SCD[[i]];
#   # if (length(x)==51){
#   for (st in 1:51){
#   y<-x[st,]
#   loc2<-stateID[st,3]
#   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_",names(SCD)[i],"_07-03-19.rds"))
#   }
#
#
#   else {
#   saveRDS(x,file = paste0("~/MITUS/inst/ST/calibration_targets/ST_",names(SCD)[i],"_07-03-19.rds"))
#   }
#   }
# }
#
# ###need to fix
# #hrcasessm
# for (st in 1:51){
# x<-SCD[["hr_cases_sm"]]
# hrsm<-matrix(NA,4,ncol(SCD[["hr_cases_sm"]]))
# loc2<-stateID[st,3]
# # for (i in 1:4){
#   y1<-x[((51*(0))+st),]
#   y2<-x[((51*(1))+st),]
#   y3<-x[((51*(2))+st),]
#   y4<-x[((51*(3))+st),]
#   y<-rbind(y1,y2,y3,y4)
#   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_hr_cases_sm_07-03-19.rds"))
#
# # }
# }
# #rt fb case
# for (i in 1:51){
#   x<-SCD[["rt_fb_cases"]]
#   rtfb<-matrix(NA,18,ncol(SCD[["rt_fb_cases"]]))
#   loc2<-stateID[i,3]
#   y<-dplyr::filter(x,st==loc2)
#   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_rt_fb_cases_07-03-19.rds"))
#
#   # }
# }
#
#
# #tbdeathsstyr
# for (i in 1:51){
#   x<-SCD[["tb_deaths_st_yr"]]
#   loc2<-stateID[i,3]
#   y<-dplyr::filter(x,State==stateID[i,1])
#   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_tb_deaths_st_yr_07-03-19.rds"))
#
#   # }
# }
#
# ##make new data directories
# # for (i in 1:length(stateID[,3])){
# #   st<-stateID[i,3]
# #   newdir<-paste0("~/MITUS/inst/",st,"/calibration_targets")
# #   dir.create(newdir)
# # }
