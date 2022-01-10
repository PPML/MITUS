#'
#'This script creates graphs that compare model simulations
#'versus calibration data.
#'@name graphs_calib
#'@param df dataframe of results
#'@return .pdf of calibration graphs
#'@export
calib_graphs <- function(df, Par_list){
  library(MCMCpack)
  df<-as.data.frame(df)
  pdfname<-paste("~/MITUS/MITUS_results/calib_graphs",Sys.time(),".pdf")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(df[1:69,30], df[1:69,31]+df[1:69,32])
  plot(1,1,ylim=c(2,600),xlim=c(1950,2018),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
  points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
  points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
  lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")
  lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],lty=3,col="blue")
  lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],lty=3,col="red3")
  lines(1950:2018,V[,2],lwd=2,col="red3")
  lines(1950:2018,V[,1],lwd=2,col="blue")
  lines(1950:2018,rowSums(V),lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=1.2)
  mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=1)
  legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[69,33:43]), t(df[69,44:54]))
  V1  <- V[-3,]
  V1[2,] <- V1[2,]+V[3,]
  V2 <- V1[-4,]
  V2[3,] <- V2[3,]+V1[4,]
  V3 <- V2[-9,]
  V3[8,] <- V3[8,]+V2[9,]

  plot(1,1,ylim=c(0.05,125),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y" )
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,5,10,25,50,75,125),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

  points(1:8+0.2,CalibDat[["tot_pop18_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  points(1:8-0.2,CalibDat[["tot_pop18_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")

  mtext("Age Group",1,2.5,cex=1.2)
  box()
  mtext("Total Population by Age Group 2018 (mil,log-scale)",3,.8,font=2,cex=1)
  legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
         pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))

  ### ### ### Population HR DISTRIBUTION 1993-2013  ### ### ### ### ### ###
  V   <- df[44:69,29]
  us_homeless<-CalibDat$homeless_pop
  plot(0,0,ylim=c(0,max(V)*1.1),xlim=c(1993,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(rep(2010,3),us_homeless,pch=19,cex=0.8,col=c(4,4,4))
  lines(1993:2018,V,lwd=2,col="blue")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Population Homeless in Past Yr",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(4,4),lty=c(NA,1),bg="white",pt.cex=0.6)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014 all ages ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  # V  <- cbind(t(df[65,33:43]),t(df[65,44:54]))
  #
  # plot(1,1,ylim=c(0.05,max(range(V))),xlim=c(0.6,11.4),xlab="",ylab="",axes=F,col=NA, log="y")
  # axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94", "95p"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  # axis(1,1:12-0.5,rep("",12))
  # axis(2,c(0,.1,1,10,20,30,50),las=2);box()
  # abline(h=axTicks(2),col="grey85")
  #
  # for(i in 1:11) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
  # for(i in 1:11) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")
  #
  # mtext("Age Group",1,2.5,cex=1.2)
  # mtext("Millions",2,2.5,cex=1.2)
  #
  # box()
  # mtext("Total Population by Age Group 2014 (mil,log-scale)",3,.8,font=2,cex=1.2)
  # legend("topright",c("US born","Non-US Born","Reported data"),cex=1.0,
  #        pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(rowSums(df[1:69,255:265]), rowSums(df[1:69,266:276]))
  V1c <- rowSums(df[1:69,121:131])
  mort_tt<-CalibDat[["US_tot_mort"]][,2]/1e6
  plot(1,1,ylim=c(0,3.5),xlim=c(1950,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2018,V[,2],lwd=2,col="red3")
  lines(1950:2018,V[,1],lwd=2,col="blue")
  lines(1950:2018,V1c,lwd=2,col="grey50")
  points(1950:2018,mort_tt,pch=19,cex=0.6,col="grey50")
  lines(1950:2018,mort_tt,lty=3,col="grey50")

  mtext("Year",1,2.5,cex=1.2)
  mtext("Mortality: Total, US, and Non-US Born (mil)",3,.8,font=2,cex=1.2)
  legend("topleft",c("Total","US born","Non-US Born","Reported data","model"),cex=1.0,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # pdfname<-paste("MITUS_results/mort_age",Sys.time(),".pdf")
  # pdf(file=pdfname, width = 11, height = 8.5)
  # par(mfrow=c(2,2),mar=c(4,4.5,3,1))
  V  <- cbind((df[68,255:265])+(df[68,266:276]))
  V1  <- V[,-3]
  V1[,2] <- V1[,2]+V[,3]
  V2 <- V1[,-4]
  V2[,3] <- V2[,3]+V1[,4]
  V3 <- V2[,-9]
  V3[,8] <- V3[,8]+V2[,9]
  V3<-V3/rowSums(V3)*100
  death_age2 <-readRDS(system.file("US/US_MortalityCountsByAge.rds", package="MITUS"))[,68]

  death_agegrp<-rep(NA,8)
  death_agegrp[1]<-sum(death_age2[1:5])
  death_agegrp[2]<-sum(death_age2[6:25])
  death_agegrp[3]<-sum(death_age2[26:45])
  death_agegrp[4]<-sum(death_age2[46:55])
  death_agegrp[5]<-sum(death_age2[56:65])
  death_agegrp[6]<-sum(death_age2[66:75])
  death_agegrp[7]<-sum(death_age2[76:85])
  death_agegrp[8]<-sum(death_age2[86:111])
  death_agegrp<-death_agegrp/sum(death_agegrp)
  # for (x in 1:67){
  plot(0,0,ylim=c(0.05,max(range(V3))*1.25),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,.2,.3,.4,.6,.8,1.0)*100,las=2);box()
  abline(h=axTicks(2),col="grey85")

  # for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[65,i],V3[65,i]),border=NA,col="gray")
  # for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[,i],V3[,i]),border=NA,col="gray")
  for(i in 1:8) points(i+.2,(death_agegrp[i])*100,pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=1.2)
  box()
  mtext(paste("Mortality by Age 2017 (%)"),3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("black","gray"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # pdfname<-paste("MITUS_results/mort_age",Sys.time(),".pdf")
  # pdf(file=pdfname, width = 11, height = 8.5)
  # par(mfrow=c(2,2),mar=c(4,4.5,3,1))
  V  <- cbind((df[68,255:265])+(df[68,266:276]))
  V1  <- V[,-3]
  V1[,2] <- V1[,2]+V[,3]
  V2 <- V1[,-4]
  V2[,3] <- V2[,3]+V1[,4]
  V3 <- V2[,-9]
  V3[,8] <- V3[,8]+V2[,9]
  death_age2 <-readRDS(system.file("US/US_MortalityCountsByAge.rds", package="MITUS"))[,69]/1e6

  death_agegrp<-rep(NA,8)
  death_agegrp[1]<-sum(death_age2[1:5])
  death_agegrp[2]<-sum(death_age2[6:25])
  death_agegrp[3]<-sum(death_age2[26:45])
  death_agegrp[4]<-sum(death_age2[46:55])
  death_agegrp[5]<-sum(death_age2[56:65])
  death_agegrp[6]<-sum(death_age2[66:75])
  death_agegrp[7]<-sum(death_age2[76:85])
  death_agegrp[8]<-sum(death_age2[86:111])
  # for (x in 1:67){
  plot(0,0,ylim=c(0.05,max(range(V3))*1.25),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,.2,.3,.4,.6,.8,1.0)*100,las=2);box()
  abline(h=axTicks(2),col="grey85")

  # for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[65,i],V3[65,i]),border=NA,col="gray")
  # for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[,i],V3[,i]),border=NA,col="gray")
  for(i in 1:8) points(i+.2,death_agegrp[i],pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=1.2)
  box()
  mtext(paste("Mortality by Age 2017"),3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("black","gray"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT RATE 1950-2013 ### ### ### ### ### ###
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # V   <- cbind(df[1:65,510:520])
  # x<-seq(6,781,12)
  # BgMort           <- Inputs[["BgMort"]]
  # mubt   <- matrix(NA,1801,11)
  # for(i in 1:11) {
  #   mubt[,i] <- SmoCurve(BgMort[,i+1])
  # }
  # V2  <-mubt[x,]
  #
  # col<-rainbow(11)
  #
  # plot(1,1,ylim=c(min(V,V2)*.5,max(V,V2)*2),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  # # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  # for (i in 1:11){
  #   lines(1950:2014,V[,i],lwd=3,col=col[i])
  #   lines(1950:2014,V2[,i],lty=3,lwd=2, col=col[i])
  #
  #   mtext("Year",1,2.5,cex=1.2)
  #   mtext("Age Specific Mortality Rates from 1950 to 2014",3,.8,font=2,cex=1.2)
  #   # legend("bottomleft",colnames(V),cex=1.0,
  #   #   pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))
  # }

  # graph of total diagnosed cases
  # by total population, US born population, and non-US born population
  V0 <- df[4:71,"NOTIF_ALL"]+df[4:71,"NOTIF_MORT_ALL"] #total population
  V1 <- df[44:71,"NOTIF_US"]+df[44:71,"NOTIF_MORT_US"]   #US born population
  V2 <- df[44:71,"NOTIF_F1"]+df[44:71,"NOTIF_F2"]+df[44:71,"NOTIF_MORT_F1"]+df[44:71,"NOTIF_MORT_F2"]   #non-US born population

  #format the plot
  plot(0,0,ylim=c(0,max(CalibDat[["tot_cases"]][,2])*1e3),xlim=c(1954,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  #multiply raw output by 1,000 to convert from millions to thousands
  lines(1953:2020,V0*1e3,lwd=3,col="white"); lines(1953:2020,V0*1e3,lwd=2,col=1) #total population
  lines(1993:2020,V1*1e3,lwd=3,col="white"); lines(1993:2020,V1*1e3,lwd=2,col=4) #US born population
  lines(1993:2020,V2*1e3,lwd=3,col="white"); lines(1993:2020,V2*1e3,lwd=2,col=3) #non-US born population

  #reported data for comparison
  points(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,pch=19,cex=0.3) #total population
  lines(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,lty=3,col=1)

  points(1993:2020,CalibDat[["age_cases_us"]][,12]/1e3,pch=19,cex=0.3,col=4) #US born population
  lines(1993:2020,CalibDat[["age_cases_us"]][,12]/1e3,pch=19,lty=3,col=4)

  points(1993:2020,CalibDat[["age_cases_fb"]][,12]/1e3,pch=19,cex=0.3,col=3) #non-US born population
  lines(1993:2020,CalibDat[["age_cases_fb"]][,12]/1e3,lty=3,col=3)


  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Cases Identified (000s), 1953-2020",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                        "Model (all)","Model (US born)","Model (non-US born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  # graph of total diagnosed cases
  # by total population, US born population, and non-US born population
  V0 <- df[61:71,"NOTIF_ALL"]+df[61:71,"NOTIF_MORT_ALL"] #total population
  V1 <- df[61:71,"NOTIF_US"]+df[61:71,"NOTIF_MORT_US"]   #US born population
  V2 <- df[61:71,"NOTIF_F1"]+df[61:71,"NOTIF_F2"]+df[61:71,"NOTIF_MORT_F1"]+df[61:71,"NOTIF_MORT_F2"]   #non-US born population

  #format the plot
  plot(0,0,ylim=c(0,max(CalibDat[["tot_cases"]][58:68,2])*1e3),xlim=c(2010,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  #multiply raw output by 1,000 to convert from millions to thousands
  lines(2010:2020,V0*1e3,lwd=3,col="white"); lines(2010:2020,V0*1e3,lwd=2,col=1) #total population
  lines(2010:2020,V1*1e3,lwd=3,col="white"); lines(2010:2020,V1*1e3,lwd=2,col=4) #US born population
  lines(2010:2020,V2*1e3,lwd=3,col="white"); lines(2010:2020,V2*1e3,lwd=2,col=3) #non-US born population

  #reported data for comparison
  points(CalibDat[["tot_cases"]][58:68,1],(CalibDat[["tot_cases"]][58:68,2])*1e3,pch=19,cex=0.3) #total population
  lines(CalibDat[["tot_cases"]][58:68,1],(CalibDat[["tot_cases"]][58:68,2])*1e3,lty=3,col=1)

  points(2010:2020,CalibDat[["age_cases_us"]][18:28,12]/1e3,pch=19,cex=0.3,col=4) #US born population
  lines(2010:2020,CalibDat[["age_cases_us"]][18:28,12]/1e3,pch=19,lty=3,col=4)

  points(2010:2020,CalibDat[["age_cases_fb"]][18:28,12]/1e3,pch=19,cex=0.3,col=3) #non-US born population
  lines(2010:2020,CalibDat[["age_cases_fb"]][18:28,12]/1e3,lty=3,col=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Cases Identified (000s), 2008-2020",3,.8,font=2,cex=1.2)
  legend("bottomleft",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                        "Model (all)","Model (US born)","Model (non-US born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  ################################################################################
  #Percent of Total Cases Non-US Born Population

  V <- cbind(df[61:71,"NOTIF_US"]+df[61:71,"NOTIF_MORT_US"], #US born population
             df[61:71,"NOTIF_F1"]+df[61:71,"NOTIF_F2"]+  #non-US born population
               df[61:71,"NOTIF_MORT_F1"]+df[61:71,"NOTIF_MORT_F2"])
  V <- V[,2]/rowSums(V)

  #format the plot
  plot(0,0,ylim=c(2.5,97.5),xlim=c(2010,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2010:2020,V*100,lwd=2,col=4)

  #reported data for comparison
  notif_fb      <- cbind(CalibDat[["age_cases_fb"]][18:28,12],CalibDat[["age_cases_us"]][18:28,12])
  points(2010:2020,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
  lines(2010:2020,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Percent of TB Cases Non-US Born, 2010-2020",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #Percent of Non-US Born Cases from Recent Immigrant Population

  V <- cbind(df[44:71,"NOTIF_F1"]+df[44:71,"NOTIF_MORT_F1"],df[44:71,"NOTIF_F2"]+df[44:71,"NOTIF_MORT_F2"])
  V <- V[,1]/rowSums(V)

  #format the plot
  plot(0,0,ylim=c(0,60),xlim=c(1993,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(1993:2020,V*100,lwd=2,col=4)

  #reported data for comparison
  notif_fb_rec2   <- cbind(CalibDat[["fb_recent_cases2"]][,2],1-CalibDat[["fb_recent_cases2"]][,2])*CalibDat[["fb_recent_cases2"]][,3]

  points(1993:2020,notif_fb_rec2[1:28,1]/rowSums(notif_fb_rec2[1:28,])*100,pch=19,cex=0.6, col="darkred")
  lines(1993:2020,notif_fb_rec2[1:28,1]/rowSums(notif_fb_rec2[1:28,])*100,lty=3, col="darkred")

  # notif_fb_rec   <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]
  #
  # points(1993:2014,notif_fb_rec[1:22,1]/rowSums(notif_fb_rec[1:22,])*100,pch=19,cex=0.6)
  # lines(1993:2014,notif_fb_rec[1:22,1]/rowSums(notif_fb_rec[1:22,])*100,lty=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=1.2)
  legend("bottomright",c("Reported data","Smoothed Data","Model"),pch=c(19,19,NA),lwd=c(1,1,2),col=c(1,"darkred",4),lty=c(3,3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #Age distribution of Cases
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

  V   <- (df[61:71,136:146]+df[61:71,189:199])
  V2  <- V[,-11]
  V2[,10] <- V2[,10]+V[,11]

  #format the plot
  cls <- colorRampPalette(c("blue", "red"))( 4 )
  plot(0,0,ylim=c(0,6),xlim=c(2010,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2010:2020,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])    #0-24 yrs
  lines(2010:2020,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])    #25-44 yrs
  lines(2010:2020,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])    #45-64 yrs
  lines(2010:2020,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])   #65+ yrs

  #reported data for comparison
  notif_age     <- CalibDat[["age_cases_us"]][18:28,-c(1,12)]*CalibDat[["age_cases_us"]][18:28,12]+
                   CalibDat[["age_cases_fb"]][18:28,-c(1,12)]*CalibDat[["age_cases_fb"]][18:28,12]
  points(2010:2020,rowSums(notif_age[,1:3])/1e3,pch=19,cex=0.6,col=cls[1]) #0-24 yrs
  lines(2010:2020,rowSums(notif_age[,1:3])/1e3,col=cls[1],lty=3)
  points(2010:2020,rowSums(notif_age[,4:5])/1e3,pch=19,cex=0.6,col=cls[2]) #25-44 yrs
  lines(2010:2020,rowSums(notif_age[,4:5])/1e3,col=cls[2],lty=3)
  points(2010:2020,rowSums(notif_age[,6:7])/1e3,pch=19,cex=0.6,col=cls[3]) #45-64 yrs
  lines(2010:2020,rowSums(notif_age[,6:7])/1e3,col=cls[3],lty=3)
  points(2010:2020,rowSums(notif_age[,8:10])/1e3,pch=19,cex=0.6,col=cls[4]) #65+ yrs
  lines(2010:2020,rowSums(notif_age[,8:10])/1e3,col=cls[4],lty=3)

  #plot text
  mtext("TB Cases By Age (000s), 2009-19",3,.8,font=2,cex=1.2)
  mtext("Year",1,2.5,cex=1.2)

  legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
         lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
         pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

  ################################################################################
  #Age Distribution of TB Cases in Percentages
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

  # V   <- (df[59:69,136:146]+df[59:69,189:199])
  # V2  <- V[,-11]
  # V2[,10] <- V2[,10]+V[,11]
  # V3<-colSums(V2)
  # V4<-rep(NA,10)
  # for (i in 1:length(V3)){
  #   V4[i]<-(V3[i]/sum(V3))*100
  # }

  Va   <- df[70:71,205:215]+df[70:71,216:226]
  Vb   <-(df[70:71,136:146]+df[70:71,189:199]) - (df[70:71,205:215]+df[70:71,216:226])
  Va2  <- Va[,-11]
  Vb2  <- Vb[,-11]

  Va2[,10] <- Va2[,10]+Va[,11]
  Vb2[,10] <- Vb2[,10]+Vb[,11]

  Va3<-colSums(Va2)
  Vb3<-colSums(Vb2)

  V4nusb<-V4usb<-rep(NA,10)

  for (i in 1:length(Va3)){
    V4usb[i]<-(Va3[i]/sum(Va3))*100
    V4nusb[i]<-(Vb3[i]/sum(Vb3))*100
  }
  #format the plot
  plot(0,0,ylim=c(0,max(range(V4usb,V4nusb))+5),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
       tick=F,cex.axis=0.6)
  axis(1,1:11-0.5,rep("",11))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:10) polygon(i+c(-.5,0,0,-.5),c(0,0,V4usb[i],V4usb[i]),border="white",col="lightblue")
  for(i in 1:10) polygon(i+c(0,.5,.5,0),c(0,0,V4nusb[i],V4nusb[i]),border="white",col="pink")

  # for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V4[i],V4[i]),border="white",col="lightblue")

  #reported data for comparison
  notif_age_us     <- colSums(CalibDat[["age_cases_us"]][27:28,-c(1,12)]*CalibDat[["age_cases"]][27:28,12])
  notif_age_fb     <- colSums(CalibDat[["age_cases_fb"]][27:28,-c(1,12)]*CalibDat[["age_cases"]][27:28,12])
  # notif_age        <- notif_age_us + notif_age_fb
  # # x<-notif_age/sum(notif_age)*100
  # points(1:10, notif_age/sum(notif_age)*100,pch=19,cex=1.2)



  #reported data for comparison
  #usb
  us_x<-(1:10)-.25
  points((1:10)-.25,(notif_age_us/sum(notif_age_us))*100,pch=19,cex=1.2, col="darkblue")
  #nusb
  nus_x<-(1:10)+.25
  points((1:10)+.25,(notif_age_fb/sum(notif_age_fb))*100,pch=19,cex=1.2, col="darkred")

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("Age Distribution of TB Cases (%), 2009-19",3,.8,font=2,cex=1.2)
  legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
         pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))


  ################################################################################
  #Average Age of TB Cases
  age_case   <- df[,136:146]+df[,205:215]
  ages<-c(2,9.5,19.5,29.5,39.5,49.5,
          59.5,69.5,79.5,89.5,99.5)
  avg_age<-rep(0,nrow(age_case))
  for (i in 1:nrow(age_case)){
    avg_age[i]<-sum(age_case[i,]*ages)/sum(age_case[i,])
  }

  #format the plot
  plot(0,0,ylim=c(0,max(avg_age)+5),xlim=c(1950,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(1950:2020,avg_age,lwd=2,col="blue")    #0-24 yrs

  #plot text
  mtext("Average Age of Notified TB Case, 1950-2020",3,.8,font=2,cex=1.2)
  mtext("Year",1,2.5,cex=1.2)

  ################################################################################
  ### ### ### CASES HR DISTRIBUTION 1993-2013  ### ### ### ### ### ###
  V <-cbind(df[44:71,685]/(df[44:71,685]+df[44:71,684]),
            df[44:71,687]/(df[44:71,687]+df[44:71,686]))
  V2<-(df[44:71,685]+df[44:71,687])/(df[44:71,687]+df[44:71,686]+df[44:71,685]+df[44:71,684])
  hr_dist_us<-CalibDat[["us_homeless_cases"]][,2]*CalibDat[["us_homeless_cases"]][,3]
  hr_dist_tot<-CalibDat[["homeless_cases"]][,2]*CalibDat[["homeless_cases"]][,3]
  #US, NUSB
  # notif_hr<-cbind(hr_dist_us,(hr_dist_tot-hr_dist_us))
  notif_hr_us<-CalibDat[["us_homeless_cases"]][,2]
  notif_hr_nus<-(hr_dist_tot-hr_dist_us)/(CalibDat[["homeless_cases"]][,3]-CalibDat[["us_homeless_cases"]][,3])
  # V   <- (df[44:69,151]+df[44:69,204])/(df[44:69,151]+df[44:69,150]+df[44:69,204]+df[44:69,203])
  notif_us_hr<-CalibDat$homeless_cases

  plot(0,0,ylim=c(0,max(V,CalibDat[["homeless_cases"]][,2])*110),xlim=c(1993,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(1993:2020,CalibDat[["homeless_cases"]][,2]*100,pch=19,cex=0.6, col="black")
  points(1993:2020,notif_hr_us*100,pch=19,cex=0.6, col="blue")
  points(1993:2020,notif_hr_nus*100,pch=19,cex=0.6, col="red")

  lines(1993:2020,CalibDat[["homeless_cases"]][,2]*100,lty=3, col="black")
  lines(1993:2020,notif_hr_us*100,lty=3, col="blue")
  lines(1993:2020,notif_hr_nus*100,lty=3, col="red")

  lines(1993:2020,V2*100,lwd=2, col="grey")
  lines(1993:2020,V[,1]*100,lwd=2, col="lightblue")
  lines(1993:2020,V[,2]*100,lwd=2,col="pink")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent of TB Cases Homeless in Past Yr",3,.8,font=2,cex=0.8)
  legend("bottomleft",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
  ################################################################################
  ### Recent infection
  #colnames(M)
  Vall <- (df[70,172:187]/df[70,156:171])
  plot(-1,0,ylim=c(0.02,1),xlim=c(0.5,16.5),xlab="",ylab="",axes=F)
  axis(2,las=2);box()
  axis(1,1:16,c("All",paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),
                            "yrs"),"US born","Foreign born","FB >2yrs","Homeless"),tick=F,cex.axis=0.7,las=2,
       mgp=c(3, 0.25, 0))

  abline(h=axTicks(2),col="grey85")
  mtext("Fraction of Incident TB from Recent Infection (<2 years) in 2020",3,.8,font=2,cex=0.8)

  for(i in 1:16) lines(rep(i,2),c(0,Vall[i]),col="forestgreen",lwd=10,lend="butt")
  rct_trans_dist<-CalibDat$recent_trans_cases[,1]/CalibDat$recent_trans_cases[,3]
  points(1,rct_trans_dist[3],pch=19,cex=1.5, col="limegreen")
  points(13,rct_trans_dist[1],pch=19,cex=1.5,col="limegreen")
  points(14,rct_trans_dist[2],pch=19,cex=1.5,col="limegreen")

  text(1:16,Vall,format(round(Vall,2),nsmall=2),cex=0.7,pos=3)
  text(c(13,14,1),rct_trans_dist,format(round(rct_trans_dist,2),nsmall=2),cex=0.7,col="black")

  legend("topright",c("2019 Reported data\n(Assumed no change)","Model"),pch=c(19,15),lwd=NA,
         pt.cex=c(1,2),col=c("limegreen","forestgreen"),bg="white")

  # ################################################################################
  # #Recent Infection 4 year time trend by nativity
  # Vus <- df[66:69,184]/df[66:69,168]
  # Vnus <- rowSums(df[66:69,185:186])/rowSums(df[66:69,169:170])
  # plot(0,0,ylim=c(0,max(Vnus,Vus)*110),xlim=c(2015,2018),xlab="",ylab="",axes=F)
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  # points(2016:2018,CalibDat[["homeless_cases"]][,2]*100,pch=19,cex=0.6, col="black")
  ################################################################################
  ### ### ### LTBI INITIATIONS 1993-2011 Distribution ### ### ### ### ### ###
  v13  <- df[43:65,153:154]/df[43:65,152]
  TLTBI_dist<-CalibDat$TLTBI_dist
  tltbi_vol<-CalibDat$TLTBI_volume
  plot(1,1,ylim=c(0.001,1)*100,xlim=c(1992,2015),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(rep(2002,2),TLTBI_dist*100,pch=19,cex=0.8,col=c("red3",4,6))
  for(i in 1:2) lines(1992:2014,v13[,i]*100,lwd=2,col=c("red3",4,6)[i])
  lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

  mtext("Year",1,2.5,cex=0.9)
  mtext("IPT Treatment Initiations By Risk Group (%)",3,.8,font=2,cex=0.8)
  legend("bottomleft",c("Foreign-born","Homeless","Reported data","Fitted model"),
         pch=c(15,15,19,NA),lwd=c(NA,NA,NA,2),col=c("red3",4,1,1),bg="white",pt.cex=c(1.8,1.8,0.8,NA))
  ################################################################################
  # Treatment Outcomes 1993-2014

  V   <- df[55:65,132:134]
  Vdisc <- V[,2]/rowSums(V)
  Vdead <- V[,3]/rowSums(V)

  #format the plot
  plot(0,0,ylim=c(0,10),xlim=c(2004,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2004:2014,Vdisc*100,lwd=2,col="red3")
  lines(2004:2014,Vdead*100,lwd=2,col="blue")

  #reported data for comparison
  points(2004:2014,CalibDat[["tx_outcomes"]][12:22,2]*100,pch=19,cex=0.6,col="red3")
  points(2004:2014,CalibDat[["tx_outcomes"]][12:22,3]*100,pch=19,cex=0.6,col="blue")
  lines (2004:2014,CalibDat[["tx_outcomes"]][12:22,2]*100,lty=3,col="red3")
  lines (2004:2014,CalibDat[["tx_outcomes"]][12:22,3]*100,lty=3,col="blue")

  #plot text

  mtext("Year",1,2.5,cex=1.2)
  mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=1.2)
  legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
         col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

  ################################################################################
  #LTBI Prevalance by Age in 2011, US born

  V <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
  pIGRA<-1
  v1<-V*pIGRA
  Sens_IGRA <-c(.780,.780,.712,.789,.789)
  Spec_IGRA <-c(.979,.979,.989,.985,.985)
  names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
  Va <- outer(v1[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v1[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))

  # Va <- outer(V[,1],c(0.74382,(1-0.74382)))+outer(V[,2],c((1-0.94014),0.94014))
  # Va <- outer(V[,1],c(Sens_IGRA,(1-Sens_IGRA)))+outer(V[,2],c((1-Spec_IGRA),Spec_IGRA))
  # Va<-V
  colnames(V) <- c("LTBI", "No-LTBI")

  V1 <- Va[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+Va[10,]+Va[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #format the plot
  plot(0,0,ylim=c(0,8.5),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #reported data for comparison
  ltbi_us_11      <- CalibDat[["LTBI_prev_US_11_IGRA"]]
  ltbi_fb_11      <- CalibDat[["LTBI_prev_FB_11_IGRA"]]

  points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("IGRA+ LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  #LTBI Prevalance by Age in 2011, non-US born
  V <- cbind(t(df[62,66:76]),t(df[62,44:54]-df[62,66:76]))
  Va <- outer(V[,1],c(0.74382,(1-0.74382)))+outer(V[,2],c((1-0.94014),0.94014))
 # make this IGRA positive
  pIGRA<-1
  v1<-V*pIGRA
  #under age 5
  v1b <- (v1[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v1[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
  #over age 5
  v1c <- outer(v1[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v1[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
  v1d<-rbind(v1b,v1c)
  colnames(v1d) <- c("LTBI", "No-LTBI")

  V1 <- v1d[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+v1d[10,]+v1d[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #format the plot
  plot(0,0,ylim=c(0,55),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #reported data for comparison
  points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=1.2)

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("IGRA+LTBI in Non-US Born Population 2011 by Age (%)",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  # Age Distribution of TB Deaths 1999-2014

  V  <- df[50:69,227:237]

  V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
  V3 <- colSums(V2)*1e6
  tb_deaths      <- CalibDat[["tb_deaths"]][,-1]

  #format the plot
  plot(0,0,ylim=c(0,max(range(V3),range(colSums(tb_deaths)))+1000),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

  #plot the model data
  for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")

  #reported data for comparison
  points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("Total TB Deaths by Age Group 1999-2018",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=NA,
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  #######################################
  # total tb deaths over time 2004-2014
  V   <- rowSums(df[59:69,227:237])
  tb_death_tot<-rowSums(CalibDat$tb_deaths[9:19,-1])

  #format the plot
  plot(0,0,ylim=c(0,max(tb_death_tot)*1.2),xlim=c(2008,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2008:2018,V*1e6,lwd=2,col="blue")

  #reported data for comparison
  points(2008:2018,tb_death_tot,pch=19,cex=0.6,col="black")
  lines (2008:2018,tb_death_tot,lty=3,col="black")

  #plot text

  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Deaths by Year 2008-2018",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
         col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))################################################################################

  ###rate of decline in TB deaths over time
  # V1<-rowSums(df[2:70,227:237])-rowSums(df[1:69,227:237])
  # V2<-(df[2:70,"NOTIF_ALL"]+df[2:70,"NOTIF_MORT_ALL"])-df[1:69,"NOTIF_ALL"]+df[1:69,"NOTIF_MORT_ALL"] #total population
  #
  # plot(0,0,ylim=c(-.04,.005),xlim=c(1951,2019),xlab="",ylab="",axes=F)
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  # #plot the model data
  # lines(1951:2019,V1,lwd=2,col="orange")
  # lines(1951:2019,V2,lwd=2,col="purple")

  ###rate of decline in TB deaths over time
  V1<-(rowSums(df[60:69,227:237])-rowSums(df[59:68,227:237]))*1e6
  V2<-((df[60:69,"NOTIF_ALL"]+df[60:69,"NOTIF_MORT_ALL"])-(df[59:68,"NOTIF_ALL"]+df[59:68,"NOTIF_MORT_ALL"]))*1e6 #total population
  case_decline<-(CalibDat$tot_cases[57:66,2]-CalibDat$tot_cases[56:65,2])*1e6
  death_decline<-rowSums(CalibDat$tb_deaths[,-1])[10:19]-rowSums(CalibDat$tb_deaths[,-1])[9:18]
  plot(0,0,ylim=c(-1800,200),xlim=c(2009,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  #plot the model data
  lines(2009:2018,V2,lwd=2,col="violet")
  lines(2009:2018,V1,lwd=2,col="orange")

  points(2009:2018,case_decline,pch=19,cex=0.6,col="purple")
  lines (2009:2018,case_decline,lty=3,col="purple")
  # abline(lm(case_decline~time), col="blue")
  points(2009:2018,death_decline,pch=19,cex=0.6,col="sienna")
  lines (2009:2018,death_decline,lty=3,col="sienna")

  mtext("Year",1,2.5,cex=1.2)
  mtext("One Year Change in TB Cases and Deaths 2008-2018",3,.8,font=2,cex=1.2)
  legend(ncol=2,"bottom",c("Reported deaths data", "Reported deaths cases","Model Deaths", "Model Cases"),pch=c(19,19,NA,NA),lwd=c(1,1,2,2),
         col=c("sienna","purple", "orange", "violet"),lty=c(3,3,1,1),bg="white",pt.cex=c(0.6,.6,NA,NA))################################################################################

  # V1<-rowSums(df[60:69,227:237])*1e6
  # V2<-((df[60:69,"NOTIF_ALL"]+df[60:69,"NOTIF_MORT_ALL"]))*1e6 #total population
  # case_decline<-(CalibDat$tot_cases[57:66,2])*1e6
  # death_decline<-rowSums(CalibDat$tb_deaths[,-1])[10:19]




  ################################################################################
  graphs_pub(Par_list=Par_list)
  # dev.off()

  # }
  dev.off()
  # system(paste("open", pdfname))# code for ma
  # pdfname<-paste("MITUS_results/pop_age",Sys.time(),".pdf")
  # pdf(file=pdfname, width = 11, height = 8.5)
  # par(mfrow=c(2,2),mar=c(4,4.5,3,1))
  # I<-as.matrix(readRDS(system.file("US/immig_pop_ag.rds", package="MITUS")))[,2:11]
  # I<-I/1e6
  # I1<-I[,-3]
  # I1[,2]<-I1[,2]+I[,3]
  # I2<-I1[,-4]
  # I2[,3]<-I2[,3]+I1[,4]
  # for (x in c(1,11,21,31,41)){
  # V  <- cbind(t(df[x,33:43]), t(df[x,44:54]))
  # V1  <- V[-3,]
  # V1[2,] <- V1[2,]+V[3,]
  # V2 <- V1[-4,]
  # V2[3,] <- V2[3,]+V1[4,]
  # V3 <- V2[-9,]
  # V3[8,] <- V3[8,]+V2[9,]
  #
  # plot(1,1,ylim=c(0.05,125),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y" )
  # axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  # axis(1,1:9-0.5,rep("",9))
  # axis(2,c(0,.5,1,3,5,10,25,50,75,125),las=2);box()
  # abline(h=axTicks(2),col="grey85")
  #
  # for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  # for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")
  #
  # # points(1:8+0.2,CalibDat[["tot_pop16_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  # points(1:8-0.2,I2[x%%9,],pch=19,cex=1.2,col="red3")
  #
  # mtext("Age Group",1,2.5,cex=1.2)
  # box()
  # mtext(paste("Total Population by Age Group,",x+1949, " (mil,log-scale)"),3,.8,font=2,cex=1)
  # legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
  #        pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))
  # }
  # dev.off()
}

calib_graphs_2020 <- function(df, Par_list){
  library(MCMCpack)
  df<-as.data.frame(df)
  pdfname<-paste("~/MITUS/MITUS_results/calib_graphs",Sys.time(),".pdf")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # V  <- cbind(df[1:69,30], df[1:69,31]+df[1:69,32])
  # plot(1,1,ylim=c(2,600),xlim=c(1950,2018),xlab="",ylab="",axes=F,log="y")
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  # points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
  # points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
  # points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
  # lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")
  # lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],lty=3,col="blue")
  # lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],lty=3,col="red3")
  # lines(1950:2018,V[,2],lwd=2,col="red3")
  # lines(1950:2018,V[,1],lwd=2,col="blue")
  # lines(1950:2018,rowSums(V),lwd=2,col="grey50")
  #
  # mtext("Year",1,2.5,cex=1.2)
  # mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=1)
  # legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
  #        pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[69,33:43]), t(df[69,44:54]))
  V1  <- V[-3,]
  V1[2,] <- V1[2,]+V[3,]
  V2 <- V1[-4,]
  V2[3,] <- V2[3,]+V1[4,]
  V3 <- V2[-9,]
  V3[8,] <- V3[8,]+V2[9,]

  plot(1,1,ylim=c(0.05,125),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y" )
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,5,10,25,50,75,125),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

  points(1:8+0.2,CalibDat[["tot_pop18_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  points(1:8-0.2,CalibDat[["tot_pop18_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")

  mtext("Age Group",1,2.5,cex=1.2)
  box()
  mtext("Total Population by Age Group 2018 (mil,log-scale)",3,.8,font=2,cex=1)
  legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
         pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))

  # graph of total diagnosed cases
  # by total population, US born population, and non-US born population
  # V0 <- df[4:71,"NOTIF_ALL"]+df[4:71,"NOTIF_MORT_ALL"] #total population
  # V1 <- df[44:71,"NOTIF_US"]+df[44:71,"NOTIF_MORT_US"]   #US born population
  # V2 <- df[44:71,"NOTIF_F1"]+df[44:71,"NOTIF_F2"]+df[44:71,"NOTIF_MORT_F1"]+df[44:71,"NOTIF_MORT_F2"]   #non-US born population
  #
  # #format the plot
  # plot(0,0,ylim=c(0,max(CalibDat[["tot_cases"]][,2])*1e3),xlim=c(1954,2020),xlab="",ylab="",axes=F)
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  #
  # #plot the model data
  # #multiply raw output by 1,000 to convert from millions to thousands
  # lines(1953:2020,V0*1e3,lwd=3,col="white"); lines(1953:2020,V0*1e3,lwd=2,col=1) #total population
  # lines(1993:2020,V1*1e3,lwd=3,col="white"); lines(1993:2020,V1*1e3,lwd=2,col=4) #US born population
  # lines(1993:2020,V2*1e3,lwd=3,col="white"); lines(1993:2020,V2*1e3,lwd=2,col=3) #non-US born population
  #
  # #reported data for comparison
  # points(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,pch=19,cex=0.3) #total population
  # lines(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,lty=3,col=1)
  #
  # points(1993:2020,CalibDat[["age_cases_us"]][,12]/1e3,pch=19,cex=0.3,col=4) #US born population
  # lines(1993:2020,CalibDat[["age_cases_us"]][,12]/1e3,pch=19,lty=3,col=4)
  #
  # points(1993:2020,CalibDat[["age_cases_fb"]][,12]/1e3,pch=19,cex=0.3,col=3) #non-US born population
  # lines(1993:2020,CalibDat[["age_cases_fb"]][,12]/1e3,lty=3,col=3)
  #
  #
  # #plot text
  # mtext("Year",1,2.5,cex=1.2)
  # mtext("Total TB Cases Identified (000s), 1953-2020",3,.8,font=2,cex=1.2)
  # legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
  #                     "Model (all)","Model (US born)","Model (non-US born)"),
  #        pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  # graph of total diagnosed cases
  # by total population, US born population, and non-US born population
  V0 <- df[61:71,"NOTIF_ALL"]+df[61:71,"NOTIF_MORT_ALL"] #total population
  V1 <- df[61:71,"NOTIF_US"]+df[61:71,"NOTIF_MORT_US"]   #US born population
  V2 <- df[61:71,"NOTIF_F1"]+df[61:71,"NOTIF_F2"]+df[61:71,"NOTIF_MORT_F1"]+df[61:71,"NOTIF_MORT_F2"]   #non-US born population

  #format the plot
  plot(0,0,ylim=c(0,max(CalibDat[["tot_cases"]][58:68,2])*1e3),xlim=c(2010,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  #multiply raw output by 1,000 to convert from millions to thousands
  lines(2010:2020,V0*1e3,lwd=3,col="white"); lines(2010:2020,V0*1e3,lwd=2,col=1) #total population
  lines(2010:2020,V1*1e3,lwd=3,col="white"); lines(2010:2020,V1*1e3,lwd=2,col=4) #US born population
  lines(2010:2020,V2*1e3,lwd=3,col="white"); lines(2010:2020,V2*1e3,lwd=2,col=3) #non-US born population

  #reported data for comparison
  points(CalibDat[["tot_cases"]][58:68,1],(CalibDat[["tot_cases"]][58:68,2])*1e3,pch=19,cex=0.3) #total population
  lines(CalibDat[["tot_cases"]][58:68,1],(CalibDat[["tot_cases"]][58:68,2])*1e3,lty=3,col=1)

  points(2010:2020,CalibDat[["age_cases_us"]][18:28,12]/1e3,pch=19,cex=0.3,col=4) #US born population
  lines(2010:2020,CalibDat[["age_cases_us"]][18:28,12]/1e3,pch=19,lty=3,col=4)

  points(2010:2020,CalibDat[["age_cases_fb"]][18:28,12]/1e3,pch=19,cex=0.3,col=3) #non-US born population
  lines(2010:2020,CalibDat[["age_cases_fb"]][18:28,12]/1e3,lty=3,col=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Cases Identified (000s), 2008-2020",3,.8,font=2,cex=1.2)
  legend("bottomleft",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                        "Model (all)","Model (US born)","Model (non-US born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  ################################################################################
  #Percent of Total Cases Non-US Born Population

  V <- cbind(df[61:71,"NOTIF_US"]+df[61:71,"NOTIF_MORT_US"], #US born population
             df[61:71,"NOTIF_F1"]+df[61:71,"NOTIF_F2"]+  #non-US born population
               df[61:71,"NOTIF_MORT_F1"]+df[61:71,"NOTIF_MORT_F2"])
  V <- V[,2]/rowSums(V)

  #format the plot
  plot(0,0,ylim=c(2.5,97.5),xlim=c(2010,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2010:2020,V*100,lwd=2,col=4)

  #reported data for comparison
  notif_fb      <- cbind(CalibDat[["age_cases_fb"]][18:28,12],CalibDat[["age_cases_us"]][18:28,12])
  points(2010:2020,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
  lines(2010:2020,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Percent of TB Cases Non-US Born, 2010-2020",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #Percent of Non-US Born Cases from Recent Immigrant Population

  # V <- cbind(df[44:71,"NOTIF_F1"]+df[44:71,"NOTIF_MORT_F1"],df[44:71,"NOTIF_F2"]+df[44:71,"NOTIF_MORT_F2"])
  # V <- V[,1]/rowSums(V)
  #
  # #format the plot
  # plot(0,0,ylim=c(0,60),xlim=c(1993,2020),xlab="",ylab="",axes=F)
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  #
  # #plot the model data
  # lines(1993:2020,V*100,lwd=2,col=4)
  #
  # #reported data for comparison
  # notif_fb_rec2   <- cbind(CalibDat[["fb_recent_cases2"]][,2],1-CalibDat[["fb_recent_cases2"]][,2])*CalibDat[["fb_recent_cases2"]][,3]
  #
  # points(1993:2020,notif_fb_rec2[1:28,1]/rowSums(notif_fb_rec2[1:28,])*100,pch=19,cex=0.6, col="darkred")
  # lines(1993:2020,notif_fb_rec2[1:28,1]/rowSums(notif_fb_rec2[1:28,])*100,lty=3, col="darkred")
  #
  # # notif_fb_rec   <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]
  # #
  # # points(1993:2014,notif_fb_rec[1:22,1]/rowSums(notif_fb_rec[1:22,])*100,pch=19,cex=0.6)
  # # lines(1993:2014,notif_fb_rec[1:22,1]/rowSums(notif_fb_rec[1:22,])*100,lty=3)
  #
  # #plot text
  # mtext("Year",1,2.5,cex=1.2)
  # mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=1.2)
  # legend("bottomright",c("Reported data","Smoothed Data","Model"),pch=c(19,19,NA),lwd=c(1,1,2),col=c(1,"darkred",4),lty=c(3,3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #Age distribution of Cases
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

  # V   <- (df[61:71,136:146]+df[61:71,189:199])
  # V2  <- V[,-11]
  # V2[,10] <- V2[,10]+V[,11]

  #format the plot
  cls <- colorRampPalette(c("blue", "red"))( 4 )
  # plot(0,0,ylim=c(0,6),xlim=c(2010,2020),xlab="",ylab="",axes=F)
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  #
  # #plot the model data
  # lines(2010:2020,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])    #0-24 yrs
  # lines(2010:2020,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])    #25-44 yrs
  # lines(2010:2020,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])    #45-64 yrs
  # lines(2010:2020,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])   #65+ yrs

  #reported data for comparison
  notif_age     <- CalibDat[["age_cases_us"]][18:28,-c(1,12)]*CalibDat[["age_cases_us"]][18:28,12]+
    CalibDat[["age_cases_fb"]][18:28,-c(1,12)]*CalibDat[["age_cases_fb"]][18:28,12]
  # points(2010:2020,rowSums(notif_age[,1:3])/1e3,pch=19,cex=0.6,col=cls[1]) #0-24 yrs
  # lines(2010:2020,rowSums(notif_age[,1:3])/1e3,col=cls[1],lty=3)
  # points(2010:2020,rowSums(notif_age[,4:5])/1e3,pch=19,cex=0.6,col=cls[2]) #25-44 yrs
  # lines(2010:2020,rowSums(notif_age[,4:5])/1e3,col=cls[2],lty=3)
  # points(2010:2020,rowSums(notif_age[,6:7])/1e3,pch=19,cex=0.6,col=cls[3]) #45-64 yrs
  # lines(2010:2020,rowSums(notif_age[,6:7])/1e3,col=cls[3],lty=3)
  # points(2010:2020,rowSums(notif_age[,8:10])/1e3,pch=19,cex=0.6,col=cls[4]) #65+ yrs
  # lines(2010:2020,rowSums(notif_age[,8:10])/1e3,col=cls[4],lty=3)
  #
  # #plot text
  # mtext("TB Cases By Age (000s), 2009-19",3,.8,font=2,cex=1.2)
  # mtext("Year",1,2.5,cex=1.2)
  #
  # legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
  #        lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
  #        pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

  ################################################################################
  #Age Distribution of TB Cases in Percentages
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

  # V   <- (df[59:69,136:146]+df[59:69,189:199])
  # V2  <- V[,-11]
  # V2[,10] <- V2[,10]+V[,11]
  # V3<-colSums(V2)
  # V4<-rep(NA,10)
  # for (i in 1:length(V3)){
  #   V4[i]<-(V3[i]/sum(V3))*100
  # }

  Va   <- df[70:71,205:215]+df[70:71,216:226]
  Vb   <-(df[70:71,136:146]+df[70:71,189:199]) - (df[70:71,205:215]+df[70:71,216:226])
  Va2  <- Va[,-11]
  Vb2  <- Vb[,-11]

  Va2[,10] <- Va2[,10]+Va[,11]
  Vb2[,10] <- Vb2[,10]+Vb[,11]

  Va3<-colSums(Va2)
  Vb3<-colSums(Vb2)

  V4nusb<-V4usb<-rep(NA,10)

  for (i in 1:length(Va3)){
    V4usb[i]<-(Va3[i]/sum(Va3))*100
    V4nusb[i]<-(Vb3[i]/sum(Vb3))*100
  }
  #format the plot
  plot(0,0,ylim=c(0,max(range(V4usb,V4nusb))+5),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
       tick=F,cex.axis=0.6)
  axis(1,1:11-0.5,rep("",11))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:10) polygon(i+c(-.5,0,0,-.5),c(0,0,V4usb[i],V4usb[i]),border="white",col="lightblue")
  for(i in 1:10) polygon(i+c(0,.5,.5,0),c(0,0,V4nusb[i],V4nusb[i]),border="white",col="pink")

  # for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V4[i],V4[i]),border="white",col="lightblue")

  #reported data for comparison
  notif_age_us     <- colSums(CalibDat[["age_cases_us"]][27:28,-c(1,12)]*CalibDat[["age_cases"]][27:28,12])
  notif_age_fb     <- colSums(CalibDat[["age_cases_fb"]][27:28,-c(1,12)]*CalibDat[["age_cases"]][27:28,12])
  # notif_age        <- notif_age_us + notif_age_fb
  # # x<-notif_age/sum(notif_age)*100
  # points(1:10, notif_age/sum(notif_age)*100,pch=19,cex=1.2)


  #reported data for comparison
  #usb
  us_x<-(1:10)-.25
  points((1:10)-.25,(notif_age_us/sum(notif_age_us))*100,pch=19,cex=1.2, col="darkblue")
  #nusb
  nus_x<-(1:10)+.25
  points((1:10)+.25,(notif_age_fb/sum(notif_age_fb))*100,pch=19,cex=1.2, col="darkred")

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("Age Distribution of TB Cases (%), 2009-19",3,.8,font=2,cex=1.2)
  legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
         pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))


  ################################################################################
  #LTBI Prevalance by Age in 2011, US born

  V <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
  pIGRA<-1
  v1<-V*pIGRA
  Sens_IGRA <-c(.780,.780,.712,.789,.789)
  Spec_IGRA <-c(.979,.979,.989,.985,.985)
  names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
  Va <- outer(v1[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v1[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))

  # Va <- outer(V[,1],c(0.74382,(1-0.74382)))+outer(V[,2],c((1-0.94014),0.94014))
  # Va <- outer(V[,1],c(Sens_IGRA,(1-Sens_IGRA)))+outer(V[,2],c((1-Spec_IGRA),Spec_IGRA))
  # Va<-V
  colnames(V) <- c("LTBI", "No-LTBI")

  V1 <- Va[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+Va[10,]+Va[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #format the plot
  plot(0,0,ylim=c(0,8.5),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #reported data for comparison
  ltbi_us_11      <- CalibDat[["LTBI_prev_US_11_IGRA"]]
  ltbi_fb_11      <- CalibDat[["LTBI_prev_FB_11_IGRA"]]

  points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("IGRA+ LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  #LTBI Prevalance by Age in 2011, non-US born
  V <- cbind(t(df[62,66:76]),t(df[62,44:54]-df[62,66:76]))
  Va <- outer(V[,1],c(0.74382,(1-0.74382)))+outer(V[,2],c((1-0.94014),0.94014))
  # make this IGRA positive
  pIGRA<-1
  v1<-V*pIGRA
  #under age 5
  v1b <- (v1[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v1[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
  #over age 5
  v1c <- outer(v1[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v1[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
  v1d<-rbind(v1b,v1c)
  colnames(v1d) <- c("LTBI", "No-LTBI")

  V1 <- v1d[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+v1d[10,]+v1d[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #format the plot
  plot(0,0,ylim=c(0,55),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #reported data for comparison
  points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=1.2)

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("IGRA+LTBI in Non-US Born Population 2011 by Age (%)",3,.8,font=2,cex=1.2)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")


  #######################################
  # total tb deaths over time 2004-2014
  V   <- rowSums(df[59:69,227:237])
  tb_death_tot<-rowSums(CalibDat$tb_deaths[9:19,-1])

  #format the plot
  plot(0,0,ylim=c(0,max(tb_death_tot)*1.2),xlim=c(2008,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2008:2018,V*1e6,lwd=2,col="blue")

  #reported data for comparison
  points(2008:2018,tb_death_tot,pch=19,cex=0.6,col="black")
  lines (2008:2018,tb_death_tot,lty=3,col="black")

  #plot text

  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Deaths by Year 2008-2018",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
         col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))################################################################################

dev.off()
}
