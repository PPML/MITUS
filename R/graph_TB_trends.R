loc_vec<-c("CA","FL","GA","IL","MA","NJ","NY","PA","TX","VA","WA")

state_fits<-function(vec){
library(MCMCpack)
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  pdfname<-paste("MITUS_results/state_TBtrend_graphs",Sys.time(),".pdf",sep="")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  for (i in 1:length(vec)) {
  loc<-vec[i]
  model_load(loc)
  if (loc %in% c("NJ","IL", "FL")){
  Opt<-readRDS(system.file(paste0(loc,"/",loc,"_Optim_all_10_923.rds"), package="MITUS"))
  } else {
  Opt<-readRDS(system.file(paste0(loc,"/",loc,"_Optim_all_10_918.rds"), package="MITUS"))
  }
  st<-which(StateID$USPS==loc)
  Par<-Opt[10,-ncol(Opt)]
  Par2 <- pnorm(Par,0,1)
  # uniform to true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  P <- P
#
  prg_chng<-def_prgchng(P)
  prms <-list()
  prms <- fin_param(P,loc,prg_chng)

  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  zz <- fin_cSim( nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]])  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                  Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
                  rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                  vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
                  HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
                  net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
                  mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                  TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                  LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                  EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
  M <- zz$Outputs
  colnames(M) <- prms[["ResNam"]]

  df<-as.data.frame(M)

#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#   ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(df[1:68,30], df[1:68,31]+df[1:68,32])
  #read in decade based stuff
  tot_pop<- CalibDatState[["pop_50_10"]][[st]]
  #get the FB pop from the decade
  tot_pop_yr_fb   <- tot_pop[tot_pop[,2]==0,]
  #get 2017 population
  pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
  #get 2017 fb population
  pop_ag_11_17nus <-sum(pop_ag_11_170[pop_ag_11_170[,2]==0,3][-11])
  #append the foreign born population
  tot_pop_yr_nus  <- c(colSums(tot_pop_yr_fb[,-c(1:2)]), pop_ag_11_17nus)/1e6

  #get the FB pop from the decade
  tot_pop_yr_us   <- tot_pop[tot_pop[,2]==1,]
  #get 2017 population
  pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
  #get 2017 fb population
  pop_ag_11_17us <-sum(pop_ag_11_170[pop_ag_11_170[,2]==1,3][-11])
  #append the us born population
  tot_pop_yr_us   <- c(colSums(tot_pop_yr_us[,-c(1:2)]), pop_ag_11_17us)/1e6

  tot_pop_yr<-(tot_pop_yr_nus+tot_pop_yr_us)
  time<-c(1950,1960,1970,1980,1990,2000,2010,2017)
  plot(1,1,ylim=c(min(V,tot_pop_yr)*.5,max(V,tot_pop_yr)*2),xlim=c(1950,2017),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(time,tot_pop_yr,pch=19,cex=0.6,col="grey50");  lines(time,tot_pop_yr,lty=3,col="grey50")
  points(time,tot_pop_yr_us,pch=19,cex=0.6,col="blue"); lines(time,tot_pop_yr_us,lty=3,col="blue")
  points(time,tot_pop_yr_nus,pch=19,cex=0.6,col="red3");lines(time,tot_pop_yr_nus,lty=3,col="red3")


  lines(1950:2017,V[,2],lwd=2,col="red3")
  lines(1950:2017,V[,1],lwd=2,col="blue")
  lines(1950:2017,rowSums(V),lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=1.2)
  mtext(paste(loc,"Population: Total, US, and Non-US Born (mil, log-scale)", sep=" "),3,.8,font=2,cex=1)
  legend("bottomleft",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
# # graph of total diagnosed cases
# by total population, US born population, and non-US born population
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
V0 <- df[44:67,"NOTIF_ALL"]+df[44:67,"NOTIF_MORT_ALL"] #total population
V1 <- df[44:67,"NOTIF_US"]+df[44:67,"NOTIF_MORT_US"]   #US born population
V2 <- df[44:67,"NOTIF_F1"]+df[44:67,"NOTIF_F2"]+df[44:67,"NOTIF_MORT_F1"]+df[44:67,"NOTIF_MORT_F2"]   #non-US born population

tot_cases<-(CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,1]+CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,2])
# tot_cases<-tot_cases/100;
#format the plot
plot(0,0,ylim=c(0,max(V0)*2*1e6),xlim=c(1993,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
#multiply raw output by 1,000 to convert from millions to hundredscali
lines(1993:2016,V0*1e6,lwd=3,col="white"); lines(1993:2016,V0*1e6,lwd=2,col=1) #total population
lines(1993:2016,V1*1e6,lwd=3,col="white"); lines(1993:2016,V1*1e6,lwd=2,col=4) #US born population
lines(1993:2016,V2*1e6,lwd=3,col="white"); lines(1993:2016,V2*1e6,lwd=2,col=3) #non-US born population

#reported data for comparison
points(1993:2016,tot_cases,pch=19,cex=0.3) #total population
lines(1993:2016,tot_cases,lty=3,col=1)

points(1993:2016,CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,"usb"],pch=19,cex=0.3,col=4) #US born population
lines(1993:2016,CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,"usb"],pch=19,lty=3,col=4)

points(1993:2016,CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,"nusb"],pch=19,cex=0.3,col=3) #non-US born population
lines(1993:2016,CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,"nusb"],lty=3,col=3)

#plot text
mtext("Year",1,2.5,cex=1.2)
mtext(paste0("Total TB Cases Identified in ", loc,", 1993-2016"),3,.8,font=2,cex=1.2)
legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                    "Model (all)","Model (US born)","Model (non-US born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)
#percentcases recent tb
V <- cbind(df[58:68,"NOTIF_F1"]+df[58:68,"NOTIF_MORT_F1"],df[58:68,"NOTIF_F2"]+df[58:68,"NOTIF_MORT_F2"])
V <- V[,1]/rowSums(V)*100
#reported data for comparison
notif_fb_per          <- CalibDatState[["rt_fb_cases"]][CalibDatState[["rt_fb_cases"]][,1]==loc,][8:18,c(2,7)]
notif_fb_per[,2]      <-notif_fb_per[,2]*100
#format the plot
plot(0,0,ylim=c(0,min(max(V,notif_fb_per[,2])*1.5,100)),xlim=c(2007,2017),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2007:2017,V,lwd=2,col=4)

#reported data for comparison
points(2007:2017,notif_fb_per[,2],pch=19,cex=0.6)
lines(2007:2017,notif_fb_per[,2],lty=3)

#plot text
mtext("Year",1,2.5,cex=1.2)
mtext(paste0("Percent of Non-US Born Cases Arrived in Past 2 Yrs in ",loc),3,.8,font=2,cex=1.2)
legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)
V   <- rowSums(df[57:67,227:237])*1e6
tb_death_tot<-CalibDatState$tbdeaths[[st]][8:18,2]
tb_death_tot[is.na(tb_death_tot)]<-0

#format the plot
plot(0,0,ylim=c(0,max(V,tb_death_tot)*1.5),xlim=c(2006,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2006:2016,V,lwd=2,col="blue")

#reported data for comparison
points(2006:2016,CalibDatState$tbdeaths[[st]][8:18,2],pch=19,cex=0.6,col="black")
lines (2006:2016,CalibDatState$tbdeaths[[st]][8:18,2],lty=3,col="black")

#plot text

mtext("Year",1,2.5,cex=1.2)
mtext(paste("Total TB Deaths in ", loc, " by Year 2006-2016", sep= " "),3,.8,font=2,cex=1.2)
legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
       col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))
print(i)
#   ################################################################################
#
}
  dev.off()
#
}
