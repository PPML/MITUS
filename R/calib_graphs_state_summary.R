calib_graphs_st_summary<-function(locvec,date){
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  # df<-as.data.frame(df)
  pdfname<-paste("MITUS_results/calib_graphs_state_summary_",date,".pdf",sep="")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))
  for (i in 1:length(locvec)){
    loc<-locvec[i] #assign loc
    st<-which(StateID$USPS==loc) #find numerical representation of this location
    model_load(loc) ##make sure to add our Opts into here
    Opt<-readRDS(system.file(paste0(loc,"/",loc,"_Optim_all_10_0713.rds"), package="MITUS"))

  #get the right run from the model
    posterior<-round(Opt[,ncol(Opt)],2); print(posterior)
    mode<-round(getmode(posterior), 2);print(mode)
    if (mode>1e11){ print(paste(loc, "did not optimize. Check optim manually", sep = " ")); next }
    run_no<-sample(which(posterior==mode), 1);print(run_no)

    Par<-Opt[run_no,-ncol(Opt)]
  #reformat the Par for model run

    Par2 <- pnorm(Par,0,1)
    # uniform to true
    Par3 <- Par2
    Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
    Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
    Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
    P[ii] <- Par3
    P <- P

    prms <-list()
    prms<-param_init(PV=P, loc=loc, prg_chng = def_prgchng(P), ttt_list = def_ttt())

    trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
    rownames(trans_mat_tot_ages) <-  paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
    colnames(trans_mat_tot_ages) <-  rep(paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4)),11)

    if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
    zz <- cSim( nYrs       = 2020-1950         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
                rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
                HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
                net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
                mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
                ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    M <- zz$Outputs
    colnames(M) <- prms[["ResNam"]]

    df<-as.data.frame(M)

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    V  <- cbind(df[1:68,30], df[1:68,31]+df[1:68,32])*1e6
    #read in decennial data
    tot_pop<- CalibDatState[["pop_50_10"]][[st]]
    #get the FB pop from the decade
    tot_pop_yr_fb   <- tot_pop[tot_pop[,2]==0,]
    #get 2017 population
    pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
    #get 2017 fb population
    pop_ag_11_17nus <-sum(pop_ag_11_170[pop_ag_11_170[,2]==0,3][-11])
    #append the foreign born population
    tot_pop_yr_nus  <- c(tot_pop_yr_fb[,-c(1:2)], pop_ag_11_17nus)

    #get the US pop from the decade
    tot_pop_yr_us   <- tot_pop[tot_pop[,2]==1,]
    #get 2017 population
    pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
    #get 2017 fb population
    pop_ag_11_17us <-sum(pop_ag_11_170[pop_ag_11_170[,2]==1,3][-11])
    #append the us born population
    tot_pop_yr_us   <- c(tot_pop_yr_us[,-c(1:2)], pop_ag_11_17us)

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
    mtext(paste("Population: Total, US, and Non-US Born (mil, log-scale) in", loc, sep=" "),3,.8,font=2,cex=1)
    legend("bottomleft",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
           pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    # graph of total diagnosed cases 2006-2016
    # by total population, US born population, and non-US born population
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    V0 <- df[59:69,"NOTIF_ALL"]+df[59:69,"NOTIF_MORT_ALL"] #total population
    V1 <- df[59:69,"NOTIF_US"]+df[59:69,"NOTIF_MORT_US"]   #US born population
    V2 <- df[59:69,"NOTIF_F1"]+df[59:69,"NOTIF_F2"]+df[59:69,"NOTIF_MORT_F1"]+df[59:69,"NOTIF_MORT_F2"]   #non-US born population

    tot_cases<-(CalibDatState$cases_yr_ag_nat_st[[st]][16:26,12,"usb"]+CalibDatState$cases_yr_ag_nat_st[[st]][16:26,12,"nusb"])
    # tot_cases<-tot_cases/100;
    #format the plot
    plot(0,0,ylim=c(min(V2,V1)*.5*1e6,max(V0)*2*1e6),xlim=c(2008,2018),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    #multiply raw output by 1,000 to convert from millions to hundredscali
    lines(2008:2018,V0*1e6,lwd=3,col="white"); lines(2008:2018,V0*1e6,lwd=2,col=1) #total population
    lines(2008:2018,V1*1e6,lwd=3,col="white"); lines(2008:2018,V1*1e6,lwd=2,col=4) #US born population
    lines(2008:2018,V2*1e6,lwd=3,col="white"); lines(2008:2018,V2*1e6,lwd=2,col=3) #non-US born population

    #reported data for comparison
    points(2008:2018,tot_cases,pch=19,cex=0.3) #total population
    lines(2008:2018,tot_cases,lty=3,col=1)

    points(2008:2018,CalibDatState$cases_yr_ag_nat_st[[st]][16:26,12,"usb"],pch=19,cex=0.3,col=4) #US born population
    lines(2008:2018,CalibDatState$cases_yr_ag_nat_st[[st]][16:26,12,"usb"],pch=19,lty=3,col=4)

    points(2008:2018,CalibDatState$cases_yr_ag_nat_st[[st]][16:26,12,"nusb"],pch=19,cex=0.3,col=3) #non-US born population
    lines(2008:2018,CalibDatState$cases_yr_ag_nat_st[[st]][16:26,12,"nusb"],lty=3,col=3)

    #plot text
    mtext("Year",1,2.5,cex=1.2)
    mtext(paste("Total TB Cases Identified, 2008-2018 in ", loc, sep=" "),3,.8,font=2,cex=1.2)
    legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                        "Model (all)","Model (US born)","Model (non-US born)"),
           pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

    ################################################################################
    #LTBI Prevalance by Age in 2011
    #specify our sensitivities and specificities (from Stout paper)
    Sens_IGRA <-c(.780,.780,.712,.789,.789)
    Spec_IGRA <-c(.979,.979,.989,.985,.985)
    names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
    #pull the model data
    Vus  <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
    Vfb  <- cbind(t(df[62,66:76]),t(df[62,44:54]-df[62,66:76]))
    #assign column names
    colnames(Vus) <-colnames(Vfb) <- c("LTBI", "No-LTBI")
    pIGRA<-1
    ###NUSB
    v1<-Vfb*pIGRA
    #under age 5
    v1b <- (v1[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v1[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
    #over age 5
    v1c <- outer(v1[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v1[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
    v1d<-rbind(v1b,v1c)
    V1nusb <- v1d[-11,]; V1nusb<-V1nusb[-10,]
    V1nusb[9,] <- V1nusb[9,]+v1d[10,]+v1d[11,]

    V2nusb <- rep(NA,8)
    V2nusb <- V1nusb[2:9,1]/rowSums(V1nusb[2:9,])*100
###USB
    v2<-Vus*pIGRA

    Va <- outer(v2[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v2[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))

    v2 <- Va[-11,]; v2<-v2[-10,]
    v2[9,] <- v2[9,]+Va[10,]+Va[11,]

    V2usb <- rep(NA,8)
    V2usb <- v2[2:9,1]/rowSums(v2[2:9,])*100

    #reported data for comparison
    ltbi_us_11      <- CalibDatState[["LTBI_prev_US_11_IGRA"]]
    ltbi_fb_11      <- CalibDatState[["LTBI_prev_FB_11_IGRA"]]

    ltbius<-ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100
    ltbifb<-ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100

    #format the plot
    plot(0,0,ylim=c(0,max(V2usb,V2nusb,ltbius, ltbifb)*1.5),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)

    axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
    axis(1,1:8-0.5,rep("",8))
    axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    for(i in 1:8) polygon(i+c(-.5,0,0,-.5),c(0,0,V2usb[i],V2usb[i]),border="white",col="lightblue")
    for(i in 1:8) polygon(i+c(0,.5,.5,0),c(0,0,V2nusb[i],V2nusb[i]),border="white",col="pink")

    #reported data for comparison
    #usb
    us_x<-(1:8)-.25
    points((1:8)-.25,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2, col="darkblue")
    for(i in 1:8) lines(us_x[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)
    #nusb
    nus_x<-(1:8)+.25
    points((1:8)+.25,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2, col="darkred")
    for(i in 1:8) lines(nus_x[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=.8)

    #plot text
    mtext("Age Group",1,2.5,cex=1.2)
    mtext(paste("IGRA+ LTBI % in 2011 by Age, Stratified by Nativity in", loc,"[NATIONAL DATA]"),3,.8,font=2,cex=.8)
    legend("topleft",c("Reported data","USB (model)", "NUSB (model)"),pch=c(19,15,15),lwd=c(0,NA,NA),
           pt.cex=c(1,2,2),col=c("black","lightblue", "pink"),bg="white")


    ################################################################################
    #tb deaths 2006-2016
    V   <- rowSums(df[57:67,227:237])*1e6
    tb_death_tot<-CalibDatState$tbdeaths[[st]][8:18,3]
    tb_death_tot[is.na(tb_death_tot)]<-0

    #format the plot
    plot(0,0,ylim=c(0,max(V,tb_death_tot)*1.5),xlim=c(2006,2016),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    lines(2006:2016,V,lwd=2,col="blue")

    #reported data for comparison
    points(2006:2016,tb_death_tot,pch=19,cex=0.6,col="black")
    lines (2006:2016,tb_death_tot,lty=3,col="black")

    #plot text

    mtext("Year",1,2.5,cex=1.2)
    mtext(paste("Total TB Deaths by Year 2006-2016 in", loc, sep=" "),3,.8,font=2,cex=1.2)
    legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
           col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))
  }
  dev.off()
}

calib_graphs_st_locs<-function(locvec,date){
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  # df<-as.data.frame(df)
  pdfname<-paste("MITUS_results/calib_graphs_states_",date,".pdf",sep="")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(3,3),mar=c(4,4,1.5,1.5))
  for (i in 1:length(locvec)){
    loc_i<-locvec[i]; print(loc_i)
    model_load(loc_i)
    Opt<-readRDS(system.file(paste0(loc_i,"/",loc_i,"_Optim_all_10_0713.rds"), package="MITUS"))

    #get the right run from the model
    posterior<-round(Opt[,ncol(Opt)],2); print(posterior)
    mode<-round(getmode(posterior), 2);print(mode)
    if (mode>1e11){ print(paste(loc, "did not optimize. Check optim manually", sep = " ")); next }
    run_no<-sample(which(posterior==mode), 1);print(run_no)

    #reformat the Par for model run

    # Par2 <- pnorm(Par,0,1)
    # # uniform to true
    # Par3 <- Par2
    # Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
    # Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
    # Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
    # P[ii] <- Par3
    # P <- P
    # #run the model
    # prms <-list()
    # prms<-param_init(PV=P, loc=loc, prg_chng = def_prgchng(P), ttt_list = def_ttt())
    #
    # trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
    # rownames(trans_mat_tot_ages) <-  paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
    # colnames(trans_mat_tot_ages) <-  rep(paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4)),11)
    #
    # if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
    # zz <- cSim( nYrs       = 2020-1950         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
    #             Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
    #             rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
    #             vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
    #             HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
    #             net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
    #             mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
    #             TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
    #             LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
    #             ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
    #             EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    # M <- zz$Outputs
    # colnames(M) <- prms[["ResNam"]]

    calib(run_no,Opt[,-40],loc_i, pdf=FALSE, cex.size=.5)
    plot.new()
    plot.new()

  }
  dev.off()
}
