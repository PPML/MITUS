##format calibration data for tabby2 plots
##united states
# model_load("US")
# NCD<-CalibDat
#
# #states
# model_load("CA")
# SCD<-readRDS(system.file("ST/ST_CalibDat_07-03-19.rds", package="MITUS"))
#
# for(i in 1:length(NCD)){
#   x<-NCD[[i]]
#   saveRDS(x,file = paste0("~/MITUS/inst/US/calibration_targets/US_",names(NCD)[i],"_04-02-21.rds"))
# }
#
# for (i in 1:length(SCD)){
#   x<-SCD[[i]];
#   # if (length(x)==51){
#   for (st in 1:51){
#   y<-x[st,]
#   loc2<-stateID[st,3]
#   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_",names(SCD)[i],"_04-02-21.rds"))
#   }
#
#
#   else {
#   saveRDS(x,file = paste0("~/MITUS/inst/ST/calibration_targets/ST_",names(SCD)[i],"_04-02-21.rds"))
#   }
#   }
# }

# ##manual reformatting of state data to match national data
# #decennial population
# # for (st in 1:51){
# #   x<-SCD[["pop_50_10"]][[st]]
# #   loc2<-stateID[st,3]
# #   pop_yr_fb<-matrix(NA,7,4);
# #   colnames(pop_yr_fb)<-c("year","total","native","foreign.born")
# #   pop_yr_fb[,1]<-seq(1950,2010,10)
# #   pop_yr_fb[,2]<-colSums(x[,3:9])/1e6
# #   pop_yr_fb[,3]<-colSums(x[x$usb==1,3:9])/1e6
# #   pop_yr_fb[,4]<-colSums(x[x$usb==0,3:9])/1e6
# #   saveRDS(pop_yr_fb,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_tot_pop_yr_fb_07-03-19.rds"))
# # }
#
# #population by age
# # for (st in 1:51){
# #   x<-SCD[["pop_00_17"]][[st]]
# #   #subset this to only 2016
# #   x<-x[,c(1,2,19)]
# #   loc2<-stateID[st,3]
# #   tot_pop16_ag_fb<-matrix(NA,9,4);
# #   tot_pop16_ag_fb<-as.data.frame(tot_pop16_ag_fb)
# #   colnames(tot_pop16_ag_fb)<-c("age_grp","total","us","fb")
# #   tot_pop16_ag_fb[,1]<-as.factor(c("0_4", "5_24","25_44","45_54","55_64","65_74", "75_84", "85p", "all"))
# #   tot_pop16_ag_fb[1,2:4]<-cbind(sum(x[x$age_group==1,3]), x[(x$age_group==1 & x$usb==1),3], x[x$age_group==1 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[2,2:4]<-cbind(sum(x[(x$age_group==2 | x$age_group==3),3]), sum(x[((x$age_group==2 | x$age_group==3)  & x$usb==1),3]), sum(x[(x$age_group==2 | x$age_group==3) & x$usb==0,3]))/1e6
# #   tot_pop16_ag_fb[3,2:4]<-cbind(sum(x[(x$age_group==4 | x$age_group==5),3]), sum(x[((x$age_group==4 | x$age_group==5)  & x$usb==1),3]), sum(x[(x$age_group==4 | x$age_group==5) & x$usb==0,3]))/1e6
# #   tot_pop16_ag_fb[4,2:4]<-cbind(sum(x[x$age_group==6,3]), x[(x$age_group==6 & x$usb==1),3], x[x$age_group==6 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[5,2:4]<-cbind(sum(x[x$age_group==7,3]), x[(x$age_group==7 & x$usb==1),3], x[x$age_group==7 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[6,2:4]<-cbind(sum(x[x$age_group==8,3]), x[(x$age_group==8 & x$usb==1),3], x[x$age_group==8 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[7,2:4]<-cbind(sum(x[x$age_group==9,3]), x[(x$age_group==9 & x$usb==1),3], x[x$age_group==9 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[8,2:4]<-cbind(sum(x[(x$age_group==10 | x$age_group==11),3]), sum(x[((x$age_group==10 | x$age_group==11)  & x$usb==1),3]), sum(x[(x$age_group==10 | x$age_group==11) & x$usb==0,3]))/1e6
# #   tot_pop16_ag_fb[9,2:4]<-cbind(sum(x[,3]), sum(x[(x$usb==1),3]), sum(x[x$usb==0,3]))/1e6
# #
# #   saveRDS(tot_pop16_ag_fb,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_tot_pop16_ag_fb_07-03-19.rds"))
# # }
#
# # for (st in 1:51){
# #   x<-SCD[["pop_00_17"]][[st]]
# #   #subset this to only 2016
# #   x<-x[,c(1,2,19)]
# #   loc2<-stateID[st,3]
# #   tot_pop16_ag_fb<-matrix(NA,9,4);
# #   tot_pop16_ag_fb<-as.data.frame(tot_pop16_ag_fb)
# #   colnames(tot_pop16_ag_fb)<-c("age_grp","total","us","fb")
# #   tot_pop16_ag_fb[,1]<-as.factor(c("0_4", "5_24","25_44","45_54","55_64","65_74", "75_84", "85p", "all"))
# #   tot_pop16_ag_fb[1,2:4]<-cbind(sum(x[x$age_group==1,3]), x[(x$age_group==1 & x$usb==1),3], x[x$age_group==1 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[2,2:4]<-cbind(sum(x[(x$age_group==2 | x$age_group==3),3]), sum(x[((x$age_group==2 | x$age_group==3)  & x$usb==1),3]), sum(x[(x$age_group==2 | x$age_group==3) & x$usb==0,3]))/1e6
# #   tot_pop16_ag_fb[3,2:4]<-cbind(sum(x[(x$age_group==4 | x$age_group==5),3]), sum(x[((x$age_group==4 | x$age_group==5)  & x$usb==1),3]), sum(x[(x$age_group==4 | x$age_group==5) & x$usb==0,3]))/1e6
# #   tot_pop16_ag_fb[4,2:4]<-cbind(sum(x[x$age_group==6,3]), x[(x$age_group==6 & x$usb==1),3], x[x$age_group==6 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[5,2:4]<-cbind(sum(x[x$age_group==7,3]), x[(x$age_group==7 & x$usb==1),3], x[x$age_group==7 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[6,2:4]<-cbind(sum(x[x$age_group==8,3]), x[(x$age_group==8 & x$usb==1),3], x[x$age_group==8 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[7,2:4]<-cbind(sum(x[x$age_group==9,3]), x[(x$age_group==9 & x$usb==1),3], x[x$age_group==9 & x$usb==0,3])/1e6
# #   tot_pop16_ag_fb[8,2:4]<-cbind(sum(x[(x$age_group==10 | x$age_group==11),3]), sum(x[((x$age_group==10 | x$age_group==11)  & x$usb==1),3]), sum(x[(x$age_group==10 | x$age_group==11) & x$usb==0,3]))/1e6
# #   tot_pop16_ag_fb[9,2:4]<-cbind(sum(x[,3]), sum(x[(x$usb==1),3]), sum(x[x$usb==0,3]))/1e6
# #
# #   saveRDS(tot_pop16_ag_fb,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_tot_pop16_ag_fb_07-03-19.rds"))
# # }
#
# # fb case
# # for (st in 1:51){
# #   loc2<-stateID[st,3]
# #   x<-SCD[["cases_yr_ag_nat_st"]][[i]][,c(1,12),"nusb"]
# #   y<-SCD[["cases_yr_ag_nat_st"]][[i]][,12,"nusb"]+SCD[["cases_yr_ag_nat_st"]][[i]][,12,"usb"]
# #   yy<-x
# #   yy[,2]<-x[,2]/y
# #   yy<-cbind(yy,y)
# #   colnames(yy)<-c("year","pct_fb","sample_size")
# #   saveRDS(yy,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_fb_cases_07-03-19.rds"))
# # }
#
# ##total mortality
# # for (st in 1:51){
# # data("stateID",package="MITUS")
# # data("ST_tot_mort",package="MITUS")
# # loc2<-stateID[st,3]
# #
# # x  <- ST_tot_mort[which(ST_tot_mort$State==stateID[st,1]),]
# # y  <- x[,3:4]
# #   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_tot_mort_07-03-19.rds"))
# # }
#
#
#
# # #hrcasessm
# # for (st in 1:51){
# # x<-SCD[["hr_cases_sm"]]
# # hrsm<-matrix(NA,4,ncol(SCD[["hr_cases_sm"]]))
# # loc2<-stateID[st,3]
# # # for (i in 1:4){
# #   y1<-x[((51*(0))+st),]
# #   y2<-x[((51*(1))+st),]
# #   y3<-x[((51*(2))+st),]
# #   y4<-x[((51*(3))+st),]
# #   y<-rbind(y1,y2,y3,y4)
# #   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_hr_cases_sm_07-03-19.rds"))
# #
# # # }
# # }
# # #recent fb case
# # for (i in 1:51){
# #   x<-SCD[["rt_fb_cases"]]
# #   rtfb<-matrix(NA,18,ncol(SCD[["rt_fb_cases"]]))
# #   loc2<-stateID[i,3]
# #   y<-dplyr::filter(x,st==loc2)
# #   yy<-cbind(y[,2],y[,7],y[,3])
# #   colnames(yy)<-c("year","pct_fb_rec","sample_size")
# #   saveRDS(yy,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_fb_recent_cases_07-03-19.rds"))
# # }
#
# # #total tb cases
# # for (i in 1:51){
# #   loc2<-stateID[i,3]
# #   x<-SCD[["cases_yr_st"]][[i]]
# #   x[,2]<-x[,2]/1e6
# #   saveRDS(x,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_cases_yr_07-03-19.rds"))
# # }
#
# # #total age tb cases
# # for (i in 1:51){
# #   loc2<-stateID[i,3]
# #   x<-SCD[["cases_yr_ag_nat_st"]][[i]][,-c(1,12),"usb"]*SCD[["cases_yr_ag_nat_st"]][[i]][,12,"usb"]+
# #      SCD[["cases_yr_ag_nat_st"]][[i]][,-c(1,12),"nusb"]*SCD[["cases_yr_ag_nat_st"]][[i]][,12,"nusb"]
# #   y<-cbind(SCD[["cases_yr_ag_nat_st"]][[i]][,1,"usb"],x/rowSums(x),rowSums(x))
# #   saveRDS(y,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_age_cases_tot_07-03-19.rds"))
# # }
#
#
# #
# # #tbdeathsstyr
# # for (i in 1:51){
# #   x<-CalibDatState$tbdeaths[[i]]
# #   loc2<-stateID[i,3]
# #   y<-dplyr::filter(x,State==stateID[i,1])
# #   dist<-matrix(unlist(CalibDatState[["tbdeaths_age_yr"]]),18,11)
# #   z<-yy<-matrix(NA,18,10)
# #   for(i in 1:18){
# #   yy[i,]<-(dist[i,2:11]/sum(dist[i,2:11]))
# #   z[i,]<-x[i,3]*yy[i,]
# #   }
# #   z<-(cbind(x[,2],z))
# #   saveRDS(z,file = paste0("~/MITUS/inst/",loc2,"/calibration_targets/",loc2,"_tb_deaths_07-03-19.rds"))
# #
# #   # }
# # }
# #
# # ##make new data directories
# # # for (i in 1:length(stateID[,3])){
# # #   st<-stateID[i,3]
# # #   newdir<-paste0("~/MITUS/inst/",st,"/calibration_targets")
# # #   dir.create(newdir)
# # # }
