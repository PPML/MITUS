#' vec<-c("CA","FL","GA","IL","MA","NJ","NY","PA","TX","VA","WA")
#' library(MCMCpack)
#'   data("stateID",package="MITUS")
#'   StateID<-as.data.frame(stateID)
#'   pdfname<-paste("MITUS_results/state_TBtrend_graphs",Sys.time(),".pdf",sep="")
#'   pdf(file=pdfname, width = 11, height = 8.5)
#'   par(mfrow=c(2,2),mar=c(4,4.5,3,1))
#'
#'   for (i in 1:length(vec)) {
#'   loc<-vec[i]
#'   model_load(loc)
#'   st_input<-paste0(loc,"_Inputs")
#'   Inputs<-get(st_input)
#'   rm(list = as.character(st_input),envir = globalenv())
#'   ParamInit<-ParamInit_st
#'   rm(ParamInit_st,envir = globalenv())
#'   opt_dat<-paste0(loc,"_Optim_all_10_2018-11-28")
#'   data(list=opt_dat, package = 'MITUS')
#'   opt<-paste0(loc,"_opt_all")
#'   opt_all<-get(opt)
#'   rm(list = as.character(opt),envir = globalenv())
#'   st<-which(StateID$USPS==loc)
#'   Par<-opt_all[10,-60]
#'   Par2 <- pnorm(Par,0,1)
#'   # uniform to true
#'   Par3 <- Par2
#'   Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
#'   Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
#'   Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
#'   P[ii] <- Par3
#'   P <- P
#'
#'   prms <-list()
#'   prms <- param(P)
#'   IP <- list()
#'   IP <- param_init(P)
#'   trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
#'   if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
#'   zz <- cSim(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
#'                Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
#'                rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
#'                vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
#'                HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
#'                mubt       = prms[["mubt"]]      , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
#'                TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
#'                LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
#'                EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
#'   #'if any output is missing or negative or if any model state population is negative
#'   #'set the likelihood to a hugely negative number (penalized)
#'
#'   M <- zz$Outputs
#'   colnames(M) <- prms[["ResNam"]]
#'
#'   df<-as.data.frame(M)
#'
#'   # graph of total diagnosed cases
#'   # by total population, US born population, and non-US born population
#'   V0 <- df[57:67,"NOTIF_ALL"]+df[57:67,"NOTIF_MORT_ALL"] #total population
#'   V1 <- df[57:67,"NOTIF_US"]+df[57:67,"NOTIF_MORT_US"]   #US born population
#'   V2 <- df[57:67,"NOTIF_F1"]+df[57:67,"NOTIF_F2"]+df[57:67,"NOTIF_MORT_F1"]+df[57:67,"NOTIF_MORT_F2"]   #non-US born population
#'
#'   tot_cases<-(CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,1]+CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,2])
#'   # tot_cases<-tot_cases/100;
#'   #format the plot
#'   plot(0,0,ylim=c(0,max(tot_cases)*1.25),xlim=c(2006,2016),xlab="",ylab="",axes=F)
#'   axis(1);axis(2,las=2);box()
#'   abline(h=axTicks(2),col="grey85")
#'
#'   #plot the model data
#'   #multiply raw output by 1,000 to convert from millions to hundredscali
#'   lines(2006:2016,V0*1e6,lwd=3,col="white"); lines(2006:2016,V0*1e6,lwd=2,col=1) #total population
#'   lines(2006:2016,V1*1e6,lwd=3,col="white"); lines(2006:2016,V1*1e6,lwd=2,col=4) #US born population
#'   lines(2006:2016,V2*1e6,lwd=3,col="white"); lines(2006:2016,V2*1e6,lwd=2,col=3) #non-US born population
#'
#'   #reported data for comparison
#'   points(2006:2016,tot_cases,pch=19,cex=0.3) #total population
#'   lines(2006:2016,tot_cases,lty=3,col=1)
#'
#'   points(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,1],pch=19,cex=0.3,col=4) #US born population
#'   lines(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,1],pch=19,lty=3,col=4)
#'
#'   points(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,2],pch=19,cex=0.3,col=3) #non-US born population
#'   lines(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,2],lty=3,col=3)
#'
#'   #plot text
#'   mtext("Year",1,2.5,cex=1.2)
#'   mtext(paste(loc,"Total TB Cases Identified, 2006-2016"),3,.8,font=2,cex=1.2)
#'   legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
#'                       "Model (all)","Model (US born)","Model (non-US born)"),
#'          pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)
#'   # total tb deaths over time 1999-2016
#'   V   <- rowSums(df[57:67,227:237])
#'   tb_death_tot<-CalibDatState$tbdeaths[[st]][8:18,c(-1,-2)]
#'   tb_death_tot[is.na(tb_death_tot)]<-0
#'
#'   #format the plot
#'   plot(0,0,ylim=c(0,max(tb_death_tot)*1.2),xlim=c(2006,2016),xlab="",ylab="",axes=F)
#'   axis(1);axis(2,las=2);box()
#'   abline(h=axTicks(2),col="grey85")
#'
#'   #plot the model data
#'   lines(2006:2016,V*1e6,lwd=2,col="blue")
#'
#'   #reported data for comparison
#'   points(2006:2016,CalibDatState$tbdeaths[[st]][8:18,c(-1,-2)],pch=19,cex=0.6,col="black")
#'   lines (2006:2016,CalibDatState$tbdeaths[[st]][8:18,c(-1,-2)],lty=3,col="black")
#'
#'   #plot text
#'
#'   mtext("Year",1,2.5,cex=1.2)
#'   mtext(paste(loc,"Total TB Deaths by Year 2006-2016"),3,.8,font=2,cex=1.2)
#'   legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
#'          col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))
#'   ################################################################################
#'
#' }
#'   dev.off()
#'
