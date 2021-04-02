calib_graphs_st_summary<-function(locvec,date){
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  # df<-as.data.frame(df)
  pdfname<-paste("MITUS_results/calib_graphs_state_summary_",date,".pdf",sep="")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))
  for (i in 1:length(locvec)){
    loc<-locvec[i] #assign loc
    print(loc)
    st<-which(StateID$USPS==loc) #find numerical representation of this location
    model_load(loc) ##make sure to add our Opts into here
    Opt<-readRDS(system.file(paste0(loc,"/",loc,"_Optim_all_10_0305.rds"), package="MITUS"))

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
                net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]], RRcrAG = prms[["RRcrAG"]],
                mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
                ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    M <- zz$Outputs
    colnames(M) <- prms[["ResNam"]]

    df<-as.data.frame(M)

    #set cex.size
    cex.size<-.75
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

    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Population: Total, US, and Non-US Born (mil, log-scale) in", loc, sep=" "),3,.8,font=2,cex=cex.size)
    legend("bottomleft",c("Total","US born","Non-US Born","Reported data","model"),cex=cex.size,
           pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    V  <- cbind(t(df[65,33:43]), t(df[65,44:54]))
    V3  <- V[-11,]
    V3[10,] <- V3[10,]+V[11,]
    pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
    pop_ag_11_17us <-pop_ag_11_170[pop_ag_11_170[,2]==1,3][-11]
    pop_ag_11_17nus <-pop_ag_11_170[pop_ag_11_170[,2]==0,3][-11]

    plot(1,1,ylim=c(min(pop_ag_11_17nus)*.5/1e6,max(pop_ag_11_17us*1.2/1e6)),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA,log="y" )
    axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
    axis(1,1:11-0.5,rep("",11))
    axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    for(i in 1:10) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
    for(i in 1:10) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

    points(1:10+0.2,pop_ag_11_17us/1e6,pch=19,cex=cex.size,col="blue")
    points(1:10-0.2,pop_ag_11_17nus/1e6,pch=19,cex=cex.size,col="red3")

    mtext("Age Group",1,2.5,cex=cex.size)
    box()
    mtext(paste("Total Population in by Age Group 2017 in",loc,"(mil,log-scale)", sep = " "),3,.3,font=2,cex=cex.size)
    legend("bottom",c("US born","Non-US Born","Reported data"),cex=cex.size,
           pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))

    # graph of total diagnosed cases 5 year bands
    # by total population, US born population, and non-US born population
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    Va <- df[44:70,"NOTIF_ALL"]+df[44:70,"NOTIF_MORT_ALL"] #total population

    # tot_cases<-rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14])+rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14])
    tot_cases<-CalibDatState[["cases_yr_st"]][[st]][,2]
    #format the plot
    plot(0,0,ylim=c(0,max(Va)*1.25*1e6),xlim=c(1993,2019),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    #multiply raw output by 1,000 to convert from millions to hundredscale
    lines(1993:2019,Va*1e6,lwd=3); lines(1993:2019,Va*1e6,lwd=2,col=1) #total population

    #reported data for comparison
    points(1993:2019,tot_cases,pch=19,cex=0.3) #total population
    lines(1993:2019,tot_cases,lty=3,col=1)

    #plot text
    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Total TB Cases Identified in", loc, ", 1993-2019", sep = " "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data (all)",
                        "Model (all)"),
           pch=c(19,NA),lwd=c(1,2),lty=c(3,1),col=c(1,1),bg="white",ncol=2,cex=cex.size,pt.cex=0.4)

    # graph of total diagnosed cases 5 year bands
    # by total population, US born population, and non-US born population
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    Vu <- df[46:70,"NOTIF_US"]+df[46:70,"NOTIF_MORT_US"]   #US born population
    Vn <- df[46:70,"NOTIF_F1"]+df[46:70,"NOTIF_F2"]+df[46:70,"NOTIF_MORT_F1"]+df[46:70,"NOTIF_MORT_F2"]   #non-US born population
    #remove create the 85+ age band
    Vn2<-Vu2<-rep(0,5)
    # Va2[1]<-sum(Va[1:5]);Va2[2]<-sum(Va[6:10]);Va2[3]<-sum(Va[11:15]); Va2[4]<-sum(Va[16:20]); Va2[5]<-sum(Va[21:25])
    Vu2[1]<-sum(Vu[1:5]);Vu2[2]<-sum(Vu[6:10]);Vu2[3]<-sum(Vu[11:15]); Vu2[4]<-sum(Vu[16:20]); Vu2[5]<-sum(Vu[21:25])
    Vn2[1]<-sum(Vn[1:5]);Vn2[2]<-sum(Vn[6:10]);Vn2[3]<-sum(Vn[11:15]); Vn2[4]<-sum(Vn[16:20]); Vn2[5]<-sum(Vn[21:25])

    # tot_cases<-rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14])+rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14])
    #format the plot
    plot(0,0,ylim=c(0,max(Vn2,Vu2)*1.75*1e6),xlim=c(1995,2020),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    #multiply raw output by 1,000 to convert from millions to hundredscale
    lines(c(1997,2002,2007,2012,2017),Vu2*1e6,lwd=3,col="white"); lines(c(1997,2002,2007,2012,2017),Vu2*1e6,lwd=2,col=4) #US born population
    lines(c(1997,2002,2007,2012,2017),Vn2*1e6,lwd=3,col="white"); lines(c(1997,2002,2007,2012,2017),Vn2*1e6,lwd=2,col=3) #non-US born population

    #reported data for comparison
    points(c(1997,2002,2007,2012,2017),rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14]),pch=19,cex=0.3,col=4) #US born population
    lines(c(1997,2002,2007,2012,2017),rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14]),pch=19,lty=3,col=4)

    points(c(1997,2002,2007,2012,2017),rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14]),pch=19,cex=0.3,col=3) #non-US born population
    lines(c(1997,2002,2007,2012,2017),rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14]),lty=3,col=3)

    #plot text
    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("TB Cases Identified in", loc," by Nativity, 1995-2019", sep = " "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data (US born)","Reported data (non-US born)",
                        "Model (US born)","Model (non-US born)"),
           pch=c(19,19,NA,NA),lwd=c(1,1,2,2),lty=c(3,3,1,1),col=c(4,3,4,3),bg="white",ncol=2,cex=cex.size,pt.cex=0.4)

    ################################################################################
    #Age Distribution of TB Cases in Percentages
    #updated for 5 year data

    V   <- (df[46:70,136:146]+df[46:70,189:199])*1e6
    V2  <- V[,-11]
    V2[,10] <- V2[,10]+V[,11]
    V2<-colSums(V2)
    V2<-(V2/sum(V2))*100
    #format the plot
    plot(0,0,ylim=c(0,max(range(V2))*1.5),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
    axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
         tick=F,cex.axis=0.6)
    axis(1,1:11-0.5,rep("",11))
    axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

    #reported data for comparison
    notif_age_us     <- matrix(unlist(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14]),5,10)
    notif_age_nus    <- matrix(unlist(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14]),5,10)
    notif_age<- notif_age_us+notif_age_nus
    notif_age_10<-colSums(notif_age[,])
    points(1:10,notif_age_10[]/sum(notif_age_10[])*100,pch=19,cex=cex.size)

    #plot text
    mtext("Age Group",1,2.5,cex=cex.size)
    mtext(paste("Age Distribution of TB Cases (%) in",loc,", 1995-19", sep = " "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
           pt.cex=c(1,2),col=c("black","lightblue"),bg="white",cex=cex.size)

    ################################################################################
    ### Recent infection
    #colnames(M)
    Vall <- (df[69,172:187]/df[69,156:171])
    plot(-1,0,ylim=c(0.02,1),xlim=c(0.5,16.5),xlab="",ylab="",axes=F)
    axis(2,las=2);box()

    axis(1,1:16,c("All",paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),
                              "yrs"),"US born","Foreign born","FB >2yrs","Homeless"),tick=F,cex.axis=0.7,las=2,
         mgp=c(3, 0.25, 0))

    abline(h=axTicks(2),col="grey85")
    mtext(paste("Fraction of Incident TB from Recent Infection (<2 years) in", loc, sep = " "),3,.3,font=2,cex=cex.size)
    #bring in the target data
    rct_trans_dist        <- CalibDat[["rct_cases_sm"]][st,5]
    for(i in 1:16) lines(rep(i,2),c(0,Vall[i]),col="forestgreen",lwd=10,lend="butt")
    points(1,rct_trans_dist,pch=19,cex=cex.size*1.5)
    text(1:16,Vall,format(round(Vall,2),nsmall=2),cex=cex.size*1.25,pos=3)
    legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
           pt.cex=c(1,2),col=c("black","forestgreen"),bg="white",cex=cex.size)
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
    mtext("Age Group",1,2.5,cex=cex.size)
    mtext(paste("IGRA+ LTBI % in 2011 by Age, Stratified by Nativity in", loc,"[NATIONAL DATA]"),3,.8,font=2,cex=cex.size)
    legend("topleft",c("Reported data","USB (model)", "NUSB (model)"),pch=c(19,15,15),lwd=c(0,NA,NA),
           pt.cex=c(1,2,2),col=c("black","lightblue", "pink"),bg="white",cex=cex.size)

    ################################################################################
    # total tb deaths over time 1999-2016
    #tb deaths 2006-2016
    V   <- rowSums(df[59:69,227:237])*1e6
    tb_death_tot<-as.numeric(CalibDatState$tbdeaths[[st]][10:20,3])
    tb_death_tot[is.na(tb_death_tot)]<-0

    #format the plot
    plot(0,0,ylim=c(0,max(V,tb_death_tot)*1.5),xlim=c(2008,2018),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    lines(2008:2018,V,lwd=2,col="blue")

    #reported data for comparison
    points(2008:2018,tb_death_tot,pch=19,cex=0.6,col="black")
    lines (2008:2018,tb_death_tot,lty=3,col="black")

    #plot text

    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Total TB Deaths by Year 2008-2018 in", loc, sep=" "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
           col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA),cex=cex.size)
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
    Opt<-readRDS(system.file(paste0(loc_i,"/",loc_i,"_Optim_all_10_0304.rds"), package="MITUS"))

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
