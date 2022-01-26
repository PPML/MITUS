
  library(MCMCpack)

  future_graphs <- function(df,endyr){
    df<-as.data.frame(df)
    pdf(file=paste("MITUS_results/future_graphs",Sys.time(),".pdf"), width = 11, height = 8.5)
    par(mfrow=c(2,2),mar=c(4,4.5,3,1))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    V  <- cbind(df[1:(endyr-1949),30], df[1:(endyr-1949),31]+df[1:(endyr-1949),32])
    plot(1,1,ylim=c(2,500),xlim=c(1950,endyr),xlab="",ylab="",axes=F,log="y")
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")
    points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
    points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
    points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
    lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")
    lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],lty=3,col="blue")
    lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],lty=3,col="red3")
    lines(1950:endyr,V[,2],lwd=2,col="red3")
    lines(1950:endyr,V[,1],lwd=2,col="blue")
    lines(1950:endyr,rowSums(V),lwd=2,col="grey50")

    mtext("Year",1,2.5,cex=0.9)
    mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=0.8)
    legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=0.9,
           pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    V  <- cbind(rowSums(df[1:(endyr-1949),255:265]), rowSums(df[1:(endyr-1949),266:276]))
    V1c <- rowSums(df[1:(endyr-1949),121:131])
    plot(1,1,ylim=c(0,max(range(V1c))),xlim=c(1950,endyr),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    lines(1950:endyr,V[,2],lwd=2,col="red3")
    lines(1950:endyr,V[,1],lwd=2,col="blue")
    lines(1950:endyr,V1c,lwd=2,col="grey50")
    points(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,pch=19,cex=0.6,col="grey50")
    lines(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,lty=3,col="grey50")

    mtext("Year",1,2.5,cex=0.9)
    mtext("Mortality: Total, US, and Non-US Born",3,.8,font=2,cex=0.8)
    legend("bottomleft",c("Total","US born","Non-US Born","Reported data","model"),cex=0.9,
           pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

    #' graph of total diagnosed cases
    #' by total population, US born population, and non-US born population
    V0 <- df[4:(endyr-1949),"NOTIF_ALL"]+df[4:(endyr-1949),"NOTIF_MORT_ALL"] #total population
    V1 <- df[44:(endyr-1949),"NOTIF_US"]+df[44:(endyr-1949),"NOTIF_MORT_US"]   #US born population
    V2 <- df[44:(endyr-1949),"NOTIF_F1"]+df[44:(endyr-1949),"NOTIF_F2"]+df[44:(endyr-1949),"NOTIF_MORT_F1"]+df[44:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

    #'format the plot
    plot(1,1,ylim=c(2,100),xlim=c(1954,endyr),xlab="",ylab="",axes=F, log = "y")
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #'plot the model data
    #'multiply raw output by 1,000 to convert from millions to thousands
    lines(1953:endyr,V0*1e3,lwd=3,col="white"); lines(1953:endyr,V0*1e3,lwd=2,col=1) #total population
    lines(1993:endyr,V1*1e3,lwd=3,col="white"); lines(1993:endyr,V1*1e3,lwd=2,col=4) #US born population
    lines(1993:endyr,V2*1e3,lwd=3,col="white"); lines(1993:endyr,V2*1e3,lwd=2,col=3) #non-US born population

    #'reported data for comparison
    points(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,pch=19,cex=0.3) #total population
    lines(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,lty=3,col=1)

    notif_fb      <- cbind(CalibDat[["fb_cases"]][,2],1-CalibDat[["fb_cases"]][,2])*CalibDat[["fb_cases"]][,3]
    notif_fb <-notif_fb/1000

    points(1993:2020,CalibDat[["age_cases_us"]][,12]/1e3,pch=19,cex=0.3,col=4) #US born population
    lines(1993:2020,CalibDat[["age_cases_us"]][,12]/1e3,pch=19,lty=3,col=4)

    points(1993:2020,CalibDat[["age_cases_fb"]][,12]/1e3,pch=19,cex=0.3,col=3) #non-US born population
    lines(1993:2020,CalibDat[["age_cases_fb"]][,12]/1e3,lty=3,col=3)

    #'plot text
    mtext("Year",1,2.5,cex=0.9)
    mtext("Total TB Cases Identified (000s)",3,.8,font=2,cex=0.8)
    legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                        "Model (all)","Model (US born)","Model (foreign born)"),
           pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

    #' graph of total diagnosed cases
    #' by total population, US born population, and non-US born population
    V0 <- df[51:(endyr-1949),"NOTIF_ALL"]+df[51:(endyr-1949),"NOTIF_MORT_ALL"] #total population
    V1 <- df[51:(endyr-1949),"NOTIF_US"]+df[51:(endyr-1949),"NOTIF_MORT_US"]   #US born population
    V2 <- df[51:(endyr-1949),"NOTIF_F1"]+df[51:(endyr-1949),"NOTIF_F2"]+df[51:(endyr-1949),"NOTIF_MORT_F1"]+df[51:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

    #'format the plot
    plot(0,0,ylim=c(0,20),xlim=c(2000,endyr),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #'plot the model data
    #'multiply raw output by 1,000 to convert from millions to thousands
    lines(2000:endyr,V0*1e3,lwd=3,col="white"); lines(2000:endyr,V0*1e3,lwd=2,col=1) #total population
    lines(2000:endyr,V1*1e3,lwd=3,col="white"); lines(2000:endyr,V1*1e3,lwd=2,col=4) #US born population
    lines(2000:endyr,V2*1e3,lwd=3,col="white"); lines(2000:endyr,V2*1e3,lwd=2,col=3) #non-US born population

    #'reported data for comparison
    points(CalibDat[["tot_cases"]][48:64,1],(CalibDat[["tot_cases"]][48:64,2])*1e3,pch=19,cex=0.3) #total population
    lines(CalibDat[["tot_cases"]][48:64,1],(CalibDat[["tot_cases"]][48:64,2])*1e3,lty=3,col=1)

    notif_fb      <- cbind(CalibDat[["fb_cases"]][8:28,2],1-CalibDat[["fb_cases"]][8:28,2])*CalibDat[["fb_cases"]][8:28,3]
    notif_fb <-notif_fb/1000

    points(2000:2020,CalibDat[["age_cases_us"]][8:28,12]/1e3,pch=19,cex=0.3,col=4) #US born population
    lines(2000:2020,CalibDat[["age_cases_us"]][8:28,12]/1e3,pch=19,lty=3,col=4)

    points(2000:2020,CalibDat[["age_cases_fb"]][8:28,12]/1e3,pch=19,cex=0.3,col=3) #non-US born population
    lines(2000:2020,CalibDat[["age_cases_fb"]][8:28,12]/1e3,lty=3,col=3)

    #'plot text
    mtext("Year",1,2.5,cex=0.9)
    mtext("Total TB Cases Identified (000s)",3,.8,font=2,cex=0.8)
    legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                        "Model (all)","Model (US born)","Model (foreign born)"),
           pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

    ################################################################################
    #'Percent of Total Cases Non-US Born Population

    V <- cbind(df[44:(endyr-1949),"NOTIF_US"]+df[44:(endyr-1949),"NOTIF_MORT_US"], #US born population
               df[44:(endyr-1949),"NOTIF_F1"]+df[44:(endyr-1949),"NOTIF_F2"]+  #non-US born population
                 df[44:(endyr-1949),"NOTIF_MORT_F1"]+df[44:(endyr-1949),"NOTIF_MORT_F2"])
    V <- V[,2]/rowSums(V)

    #'format the plot
    plot(0,0,ylim=c(2.5,97.5),xlim=c(2000,endyr),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #'plot the model data
    lines(1993:endyr,V*100,lwd=2,col=4)

    #reported data for comparison
    notif_fb      <- cbind(CalibDat[["age_cases_fb"]][18:28,12],CalibDat[["age_cases_us"]][18:28,12])
    points(2010:2020,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
    lines(2010:2020,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

    #'plot text
    mtext("Year",1,2.5,cex=0.9)
    mtext("Percent of TB Cases Non-US-Born",3,.8,font=2,cex=0.8)
    legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

    ################################################################################
    #'Percent of Non-US Born Cases from Recent Immigrant Population

    V <- cbind(df[44:(endyr-1949),"NOTIF_F1"]+df[44:(endyr-1949),"NOTIF_MORT_F1"],df[44:(endyr-1949),"NOTIF_F2"]+df[44:(endyr-1949),"NOTIF_MORT_F2"])
    V <- V[,1]/rowSums(V)

    #'format the plot
    plot(0,0,ylim=c(0,60),xlim=c(1993,endyr),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #'plot the model data
    lines(1993:endyr,V*100,lwd=2,col=4)

    #reported data for comparison
    notif_fb_rec2   <- cbind(CalibDat[["fb_recent_cases2"]][,2],1-CalibDat[["fb_recent_cases2"]][,2])*CalibDat[["fb_recent_cases2"]][,3]

    points(1993:2020,notif_fb_rec2[1:28,1]/rowSums(notif_fb_rec2[1:28,])*100,pch=19,cex=0.6, col="darkred")
    lines(1993:2020,notif_fb_rec2[1:28,1]/rowSums(notif_fb_rec2[1:28,])*100,lty=3, col="darkred")

    #'plot text
    mtext("Year",1,2.5,cex=0.9)
    mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=0.8)
    legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)


    ################################################################################
    #' Treatment Outcomes 1993-2014

    #' V   <- df[44:(endyr-1949),132:134]
    #' Vdisc <- V[,2]/rowSums(V)
    #' Vdead <- V[,3]/rowSums(V)
    #'
    #' #'format the plot
    #' plot(0,0,ylim=c(0,15),xlim=c(1993,endyr),xlab="",ylab="",axes=F)
    #' axis(1);axis(2,las=2);box()
    #' abline(h=axTicks(2),col="grey85")
    #'
    #' #'plot the model data
    #' lines(1993:endyr,Vdisc*100,lwd=2,col="red3")
    #' lines(1993:endyr,Vdead*100,lwd=2,col="blue")
    #'
    #' #'reported data for comparison
    #' tx_outcomes      <- cbind(1-rowSums(CalibDat[["tx_outcomes"]][,2:3]),CalibDat[["tx_outcomes"]][,2],CalibDat[["tx_outcomes"]][,3])*CalibDat[["tx_outcomes"]][,4]
    #'
    #' points(1993:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
    #' points(1993:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
    #' lines (1993:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
    #' lines (1993:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")
    #'
    #' #'plot text
    #'
    #' mtext("Year",1,2.5,cex=0.9)
    #' mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=0.8)
    #' legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
    #'        col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))


    # total tb deaths over time 2004-2014
    V   <- rowSums(df[59:(endyr-1949),227:237])*1e6
    tb_death_tot<-rowSums(CalibDat$tb_deaths[9:19,-1])

    #format the plot
    plot(0,0,ylim=c(0,max(V,tb_death_tot)*1.2),xlim=c(2006,2050),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    lines(2008:2049,V,lwd=2,col="blue")

    #reported data for comparison
    points(2008:2018,tb_death_tot,pch=19,cex=0.6,col="black")
    lines (2008:2018,tb_death_tot,lty=3,col="black")

    #plot text

    mtext("Year",1,2.5,cex=1.2)
    mtext("Total TB Deaths by Year 2008-2050",3,.8,font=2,cex=1.2)
    legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
           col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))################################################################################


    dev.off() # code for ma

  }

  future_graphs_st <-function(loc,df, cex.size = .7){
    data("stateID",package="MITUS")
    StateID<-as.data.frame(stateID)
    st<-which(StateID$USPS==loc)
    df<-as.data.frame(df)

    pdfname<-paste("MITUS_results/",loc,"_future_graphs_", Sys.Date(), ".pdf",sep="")
    pdf(file=pdfname, width = 11, height = 8.5)
    par(mfrow=c(2,2),mar=c(4,4.5,3,1))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    V  <- cbind(df[1:100,30], df[1:100,31]+df[1:100,32])*1e6
    #read in decade based stuff
    tot_pop<- CalibDatState[["pop_50_10"]][[st]]
    #get the FB pop from the decade
    tot_pop_yr_fb   <- tot_pop[tot_pop[,2]==0,]
    #get 2017 population
    pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c(1,2,22)]
    #get 2020 fb population
    pop_ag_11_19nus <-sum(pop_ag_11_190[pop_ag_11_190[,2]==0,3][-11])
    #append the foreign born population
    tot_pop_yr_nus  <- c(tot_pop_yr_fb[,-c(1:2)], pop_ag_11_19nus)

    #get the FB pop from the decade
    tot_pop_yr_us   <- tot_pop[tot_pop[,2]==1,]
    #get 2020 population
    pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c(1,2,20)]
    #get 2020 fb population
    pop_ag_11_19us <-sum(pop_ag_11_190[pop_ag_11_190[,2]==1,3][-11])
    #append the us born population
    tot_pop_yr_us   <- c(tot_pop_yr_us[,-c(1:2)], pop_ag_11_19us)

    tot_pop_yr<-(tot_pop_yr_nus+tot_pop_yr_us)
    time<-c(1950,1960,1970,1980,1990,2000,2010,2020)
    plot(1,1,ylim=c(min(V,tot_pop_yr)*.15,max(V,tot_pop_yr)*2),xlim=c(1950,2050),xlab="",ylab="",axes=F,log="y")
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")
    points(time,tot_pop_yr,pch=19,cex=0.6,col="grey50");  lines(time,tot_pop_yr,lty=3,col="grey50")
    points(time,tot_pop_yr_us,pch=19,cex=0.6,col="blue"); lines(time,tot_pop_yr_us,lty=3,col="blue")
    points(time,tot_pop_yr_nus,pch=19,cex=0.6,col="red3");lines(time,tot_pop_yr_nus,lty=3,col="red3")


    lines(1950:2049,V[,2],lwd=2,col="red3")
    lines(1950:2049,V[,1],lwd=2,col="blue")
    lines(1950:2049,rowSums(V),lwd=2,col="grey50")

    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Population in: Total, US, and Non-US Born in",loc,"mil, log-scale)", sep = " "),3,.3,font=2,cex=cex.size)

    legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=cex.size,
           pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),
           bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA), ncol=2)

    # graph of total diagnosed cases 5 year bands
    # by total population, US born population, and non-US born population
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    Va <- df[61:100,"NOTIF_ALL"]+df[61:100,"NOTIF_MORT_ALL"] #total population
    Vu <- df[61:100,"NOTIF_US"]+df[61:100,"NOTIF_MORT_US"]   #US born population
    Vn <- df[61:100,"NOTIF_F1"]+df[61:100,"NOTIF_F2"]+df[61:100,"NOTIF_MORT_F1"]+df[61:100,"NOTIF_MORT_F2"]   #non-US born population
    # tot_cases<-rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14])+rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14])
    tot_cases<-CalibDatState[["cases_yr_st"]][[st]][,2]
    #format the plot
    plot(0,0,ylim=c(0,max(Va*1e6,tot_cases[18:28])*1.25),xlim=c(2010,2049),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    #multiply raw output by 1,000 to convert from millions to hundredscale
    lines(2010:2049,Va*1e6,lwd=3); lines(2010:2049,Va*1e6,lwd=2,col=1) #total population
    lines(2010:2049,Vu*1e6,lwd=3); lines(2010:2049,Vu*1e6,lwd=2,col=4) #total population
    lines(2010:2049,Vn*1e6,lwd=3); lines(2010:2049,Vn*1e6,lwd=2,col=3) #total population

    #reported data for comparison
    points(2010:2020,tot_cases[18:28],pch=19,cex=0.3) #total population
    lines(2010:2020,tot_cases[18:28],lty=3,col=1)

    #plot text
    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Total TB Cases Identified in", loc, ", 1993-2050", sep = " "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data (all)",
                        "Model (all)"),
           pch=c(19,NA),lwd=c(1,2),lty=c(3,1),col=c(1,1),bg="white",ncol=2,cex=cex.size,pt.cex=0.4)


    #Percent of Non-US Born Cases from Recent Immigrant Population
    #updated for 5 year data
    #check this plot and otis
    V <- cbind(df[61:100,"NOTIF_F1"]+df[61:100,"NOTIF_MORT_F1"],df[61:100,"NOTIF_F2"]+df[61:100,"NOTIF_MORT_F2"])
    #create a five year band for this data

    V0<-V/rowSums(V)*100

    #format the plot
    # plot(0,0,ylim=c(min(V,notif_fb_rec[,1])*.25,min(max(V,notif_fb_rec[,1])*1.5,100)),xlim=c(2007,2017),xlab="",ylab="",axes=F)
    plot(0,0,ylim=c(0,min(max(V0[,1])*1.5,100)),xlim=c(2010,2049),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    lines(2010:2049,V0[,1],lwd=2,col=4)
    #plot text
    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Percent of Non-US Born Cases Arrived in Past 2 Yrs in",loc, sep = " "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),
           lty=c(3,1),bg="white",pt.cex=0.6,cex=cex.size)

    ################################################################################
    # total tb deaths over time 1999-2016
    #tb deaths 2006-2016
    V   <- rowSums(df[61:100,227:237])*1e6
    tb_death_tot<-as.numeric(CalibDatState$tbdeaths[[st]][12:20,3])
    tb_death_tot[is.na(tb_death_tot)]<-0

    #format the plot
    plot(0,0,ylim=c(0,max(V,tb_death_tot)*1.5),xlim=c(2010,2049),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    lines(2010:2049,V,lwd=2,col="blue")

    #reported data for comparison
    points(2010:2018,tb_death_tot,pch=19,cex=0.6,col="black")
    lines (2010:2018,tb_death_tot,lty=3,col="black")

    #plot text
    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Total TB Deaths by Year 2008-2020 in", loc, sep=" "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
           col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA),cex=cex.size)
    dev.off()
  }


