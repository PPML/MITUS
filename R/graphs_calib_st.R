#'
#'This script creates graphs that compare model simulations
#'versus calibration data.
#'@name graphs_calib_st
#'@param df dataframe of results
#'@return .pdf of calibration graphs
#'@export
calib_graphs_st <- function(df,loc, Par_list){

  library(MCMCpack)
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  st<-which(StateID$USPS==loc)
  df<-as.data.frame(df)
  pdfname<-paste("MITUS_results/",loc,"_calib_graphs_",Sys.time(),".pdf",sep="")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

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
  mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=1)
  legend("bottomleft",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
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

  points(1:10+0.2,pop_ag_11_17us/1e6,pch=19,cex=1.2,col="blue")
  points(1:10-0.2,pop_ag_11_17nus/1e6,pch=19,cex=1.2,col="red3")

  mtext("Age Group",1,2.5,cex=1.2)
  box()
  mtext("Total Population by Age Group 2017 (mil,log-scale)",3,.8,font=2,cex=1)
  legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
         pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL FB POP EACH DECADE, BY REC/LONG   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  # V  <- cbind(df[1:68,31],df[1:68,32])
#
#   plot(0,0,ylim=c(min((V*.5),0),max(rowSums(V))*1.25),xlim=c(1950,2017),xlab="",ylab="",axes=F)
#   axis(1);axis(2,las=2);box()
#   abline(h=axTicks(2),col="grey85")
#
#   lines(1950:2017,V[,2],lwd=2,col="red3")
#   lines(1950:2017,V[,1],lwd=2,col="blue")
#   lines(1950:2017,rowSums(V),lwd=2,col="grey50")
#
#   mtext("Year",1,2.5,cex=1.2)
#   mtext("Non-US Born Population: Total, Recent, Long-Term (mil)",3,.8,font=2,cex=1)
#   legend("topleft",c("Total","Recent","Long-Term","model"),cex=1,
#          pch=c(15,15,15,NA),lwd=c(NA,NA,NA,2),lty=c(NA,NA,NA,1),col=c("grey50",4,"red3",1),bg="white",pt.cex=c(1.8,1.8,1.8,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL FB RECENT POP EACH DECADE  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # V1<-V*1e3
  # if (loc =="CA"){
  # est<-c(623500,596598, 569657, 544240,534542, 544490, 561243)/1e3
  # } else if (loc=="TX"){
  #   est<-c(428585,
  #          410480,
  #          391566,
  #          375654,
  #          367044,
  #          364898,
  #          364670)/1e3
  # } else if (loc=="MA"){
  #   est<-c(133153,
  #          127403,
  #          124165,
  #          118684,
  #          114787,
  #          114683,
  #          113249)/1e3
  # } else if (loc=="NY"){
  #   est<-c(311010,
  #          311062,
  #          307475,
  #          302396,
  #          298544,
  #          295435,
  #          290282)/1e3
  # } else {est<-rep(0,7)}
  #
  # plot(0,0,ylim=c(min((V1[,1]*.5),0),max(V1[,1],est)*1.25),xlim=c(1950,2015),xlab="",ylab="",axes=F)
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  #
  # lines(1950:2017,V1[,1],lwd=2,col="blue")
  # if (loc =="CA" |loc == "NY"|loc == "MA"|loc =="TX"){
  #   points(2011:2017,est,pch=19,cex=0.6,col="blue")
  #   lines(2011:2017,est,lty=3,col="blue")
  # }
  # mtext("Year",1,2.5,cex=1.2)
  #
  # mtext("Recent Immigration Non-US Born Population (<2 yrs) (000s)",3,.8,font=2,cex=1)
  # legend("topleft",c("Total"),cex=1,
  #        pch=c(15),lwd=c(NA),lty=c(NA),col=c("blue"),bg="white",pt.cex=c(1.8))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL FB POP EACH DECADE, BY REC/LONG   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#   V  <- cbind(df[1:66,31],df[1:66,32])
#
#   plot(1,1,ylim=c(min((V*.5),1e-4),max(rowSums(V))*2),xlim=c(1950,2015),xlab="",ylab="",axes=F, log="y")
#   axis(1);axis(2,las=2);box()
#   abline(h=axTicks(2),col="grey85")
#
#   lines(1950:2015,V[,2],lwd=2,col="red3")
#   lines(1950:2015,V[,1],lwd=2,col="blue")
#   lines(1950:2015,rowSums(V),lwd=2,col="grey50")
#
#   mtext("Year",1,2.5,cex=1.2)
#   mtext("Non-US Born Population: Total, Recent, Long-Term (mil, log-scale)",3,.8,font=2,cex=1)
#   legend("bottomright",c("Total","Recent","Long-Term","model"),cex=1,
#          pch=c(15,15,15,NA),lwd=c(NA,NA,NA,2),lty=c(NA,NA,NA,1),col=c("grey50",4,"red3",1),bg="white",pt.cex=c(1.8,1.8,1.8,NA))
#
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#   ### ### ### ### ### ###   TOTAL FB POP EACH DECADE, BY REC/LONG % ### ### ### ### ### ###
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#   V  <- cbind(df[1:66,31],df[1:66,32])
#   V <-V/rowSums(V)*100
#
#   plot(0,0,ylim=c(max(min(V)*.5,0),min(max(V)*2,100)),xlim=c(1950,2015),xlab="",ylab="",axes=F)
#   axis(1);axis(2,las=2);box()
#   abline(h=axTicks(2),col="grey85")
#
#   lines(1950:2015,V[,2],lwd=2,col="red3")
#   lines(1950:2015,V[,1],lwd=2,col="blue")
#   # lines(1950:2015,rowSums(V),lwd=2,col="grey50")
#
#   mtext("Year",1,2.5,cex=1.2)
#   mtext("Non-US Born Population: Recent, Long-Term (%)",3,.8,font=2,cex=1)
#   legend("left",c("Recent","Long-Term","model"),cex=1,
#          pch=c(15,15,NA),lwd=c(NA,NA,2),lty=c(NA,NA,1),col=c(4,"red3",1),bg="white",pt.cex=c(1.8,1.8,NA))
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
  # V  <- cbind(rowSums(df[30:67,255:265]), rowSums(df[30:67,266:276]))*1e6
  V1c <- rowSums(df[1:67,121:131])
  #1979-2016 total deaths
  ST_deaths_tot <- readRDS(system.file("ST/STdeathbyAge.rds",package="MITUS"))[[st]][,c(1,12)]
  ST_deaths_tot[,2]<-ST_deaths_tot[,2]/1e6
  plot(1,1,ylim=c(min(V1c,ST_deaths_tot[,2])*.5,max(V1c,ST_deaths_tot[,2])*2),xlim=c(1950,2016),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  lines(1950:2016,V1c,lwd=2,col="grey50")
  points(ST_deaths_tot[,1],(ST_deaths_tot[,2]),pch=19,cex=0.6,col="grey50")
  lines(ST_deaths_tot[,1],(ST_deaths_tot[,2]),lty=3,col="grey50")

  mtext("Year",1,2.5,cex=1.2)
  mtext("Total Mortality",3,.8,font=2,cex=1.2)
  legend("topright",c("Total","Reported data","model"),cex=1.0,
         pch=c(15,19,NA),lwd=c(NA,1,2),lty=c(NA,3,1),col=c("grey50",1,1),bg="white",pt.cex=c(1.8,0.3,NA))
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2017  ### ### ### ### ### ###
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind((df[68,255:265])+(df[68,266:276]))
  V1<-V[-11]
  V1[10] <- V1[10]+V[11]
  V2<-V1/sum(V1)

  tda <- readRDS(system.file("ST/STdeathbyAge.rds",package="MITUS"))[[st]][48,-c(1,12)]
  tda<-tda/sum(tda)
  plot(0,0,ylim=c(min(range(V2,tda))*.5,max(range(V2,tda))*1.25),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:11-0.5,rep("",11))
  axis(2,c(0,.2,.4,.6,.8),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:10) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V2[i],V2[i]),border=NA,col="gray")
  for(i in 1:10) points(i+.2,(tda[i]),pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=1.2)
  box()
  mtext("Mortality Distribution by Age, 2016 (%) [NATIONAL]",3,.8,font=2,cex=1.2)
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
  #   mubt[,i] <- SmoCurve(BgMort[,i+1])*P["TunMubt"]/12
  # }
  # V2  <-mubt[x,]
  #
  # col<-rainbow(11)
  #
  # plot(1,1,ylim=c(.00001,.04),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
  # axis(1);axis(2,las=2);box()
  # abline(h=axTicks(2),col="grey85")
  # # points(1993:2014,CalibDatState$notif_us_hr[,1]/rowSums(CalibDatState$notif_us_hr)*100,pch=19,cex=0.6)
  # for (i in 1:11){
  #   lines(1950:2014,V[,i],lwd=3,col=col[i])
  #   lines(1950:2014,V2[,i],lty=3,lwd=2, col="black")
  #
  #   mtext("Year",1,2.5,cex=1.2)
  #   mtext("Age Specific Mortality Rates from 1950 to 2014",3,.8,font=2,cex=1.2)
  #   legend("bottomleft",colnames(V),cex=1.0,
  #          pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))
  # }
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # graph of total diagnosed cases
  # by total population, US born population, and non-US born population
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V0 <- df[44:67,"NOTIF_ALL"]+df[44:67,"NOTIF_MORT_ALL"] #total population
  V1 <- df[44:67,"NOTIF_US"]+df[44:67,"NOTIF_MORT_US"]   #US born population
  V2 <- df[44:67,"NOTIF_F1"]+df[44:67,"NOTIF_F2"]+df[44:67,"NOTIF_MORT_F1"]+df[44:67,"NOTIF_MORT_F2"]   #non-US born population

  tot_cases<-(CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,1]+CalibDatState$cases_yr_ag_nat_st[[st]][1:24,12,2])
  # tot_cases<-tot_cases/100;
  #format the plot
  plot(0,0,ylim=c(min(V2,V1)*.5*1e6,max(V0)*2*1e6),xlim=c(1993,2016),xlab="",ylab="",axes=F)
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
  mtext("Total TB Cases Identified, 1993-2016",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                      "Model (all)","Model (US born)","Model (non-US born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # graph of total diagnosed cases 2006-2016
  # by total population, US born population, and non-US born population
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V0 <- df[57:67,"NOTIF_ALL"]+df[57:67,"NOTIF_MORT_ALL"] #total population
  V1 <- df[57:67,"NOTIF_US"]+df[57:67,"NOTIF_MORT_US"]   #US born population
  V2 <- df[57:67,"NOTIF_F1"]+df[57:67,"NOTIF_F2"]+df[57:67,"NOTIF_MORT_F1"]+df[57:67,"NOTIF_MORT_F2"]   #non-US born population

  tot_cases<-(CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"usb"]+CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"nusb"])
  # tot_cases<-tot_cases/100;
  #format the plot
  plot(0,0,ylim=c(min(V2,V1)*.5*1e6,max(V0)*2*1e6),xlim=c(2006,2016),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  #multiply raw output by 1,000 to convert from millions to hundredscali
  lines(2006:2016,V0*1e6,lwd=3,col="white"); lines(2006:2016,V0*1e6,lwd=2,col=1) #total population
  lines(2006:2016,V1*1e6,lwd=3,col="white"); lines(2006:2016,V1*1e6,lwd=2,col=4) #US born population
  lines(2006:2016,V2*1e6,lwd=3,col="white"); lines(2006:2016,V2*1e6,lwd=2,col=3) #non-US born population

  #reported data for comparison
  points(2006:2016,tot_cases,pch=19,cex=0.3) #total population
  lines(2006:2016,tot_cases,lty=3,col=1)

  points(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"usb"],pch=19,cex=0.3,col=4) #US born population
  lines(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"usb"],pch=19,lty=3,col=4)

  points(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"nusb"],pch=19,cex=0.3,col=3) #non-US born population
  lines(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"nusb"],lty=3,col=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Cases Identified, 2006-2016",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                        "Model (all)","Model (US born)","Model (non-US born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  ################################################################################
  #Percent of Total Cases Non-US Born Population

  V <- cbind(df[57:67,"NOTIF_US"]+df[57:67,"NOTIF_MORT_US"], #US born population
             df[57:67,"NOTIF_F1"]+df[57:67,"NOTIF_F2"]+  #non-US born population
               df[57:67,"NOTIF_MORT_F1"]+df[57:67,"NOTIF_MORT_F2"])
  V <- V[,2]/rowSums(V)

  #format the plot
  plot(0,0,ylim=c(min(V)*.5*100,min((V*2*100),100)),xlim=c(2006,2016),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2006:2016,V*100,lwd=2,col=4)

  #reported data for comparison
  points(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"nusb"]/tot_cases*100,pch=19,cex=0.6)
  lines(2006:2016,CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,"nusb"]/tot_cases*100,lty=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Percent of TB Cases Non-US Born, 2006-2016",3,.8,font=2,cex=1.2)
  legend("bottomright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #Percent of Non-US Born Cases from Recent Immigrant Population

  V <- cbind(df[58:68,"NOTIF_F1"]+df[58:68,"NOTIF_MORT_F1"],df[58:68,"NOTIF_F2"]+df[58:68,"NOTIF_MORT_F2"])
  V <- V[,1]/rowSums(V)*100
  #reported data for comparison
  notif_fb_per          <- CalibDatState[["rt_fb_cases"]][CalibDatState[["rt_fb_cases"]][,1]==loc,][8:18,c(2,7)]
  notif_fb_per[,2]      <-notif_fb_per[,2]*100
  #format the plot
  plot(0,0,ylim=c(min(V,notif_fb_per[,2])*.25,min(max(V,notif_fb_per[,2])*1.5,100)),xlim=c(2007,2017),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2007:2017,V,lwd=2,col=4)

  #reported data for comparison
  points(2007:2017,notif_fb_per[,2],pch=19,cex=0.6)
  lines(2007:2017,notif_fb_per[,2],lty=3)

  #plot text
  mtext("Year",1,2.5,cex=1.2)
  mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #Age distribution of Cases
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs
  V   <- (df[57:67,136:146]+df[57:67,189:199])*1e6
  V2  <- V[,-11]
  V2[,10] <- V2[,10]+V[,11]
  notif_age_us     <- CalibDatState$cases_yr_ag_nat_st[[st]][14:24,-c(1,12),1]*CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,1]
  notif_age_nus    <- CalibDatState$cases_yr_ag_nat_st[[st]][14:24,-c(1,12),2]*CalibDatState$cases_yr_ag_nat_st[[st]][14:24,12,2]
  notif_age<- notif_age_us+notif_age_nus
  #format the plot
  cls <- colorRampPalette(c("blue", "red"))( 4 )
  plot(0,0,ylim=c(min(V2,notif_age)*.5,max(V2,notif_age)*2),xlim=c(2006,2016),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2006:2016,rowSums(V2[,1:3]),lwd=2,col=cls[1])    #0-24 yrs
  lines(2006:2016,rowSums(V2[,4:5]),lwd=2,col=cls[2])    #25-44 yrs
  lines(2006:2016,rowSums(V2[,6:7]),lwd=2,col=cls[3])    #45-64 yrs
  lines(2006:2016,rowSums(V2[,8:10]),lwd=2,col=cls[4])   #65+ yrs

  #reported data for comparison

  points(2006:2016,rowSums(notif_age[,1:3]),pch=19,cex=0.6,col=cls[1]) #0-24 yrs
  lines(2006:2016,rowSums(notif_age[,1:3]),col=cls[1],lty=3)
  points(2006:2016,rowSums(notif_age[,4:5]),pch=19,cex=0.6,col=cls[2]) #25-44 yrs
  lines(2006:2016,rowSums(notif_age[,4:5]),col=cls[2],lty=3)
  points(2006:2016,rowSums(notif_age[,6:7]),pch=19,cex=0.6,col=cls[3]) #45-64 yrs
  lines(2006:2016,rowSums(notif_age[,6:7]),col=cls[3],lty=3)
  points(2006:2016,rowSums(notif_age[,8:10]),pch=19,cex=0.6,col=cls[4]) #65+ yrs
  lines(2006:2016,rowSums(notif_age[,8:10]),col=cls[4],lty=3)

  #plot text
  mtext("TB Cases By Age, 2006-16",3,.8,font=2,cex=1.2)
  mtext("Year",1,2.5,cex=1.2)

  legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
         lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
         pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

  ################################################################################
  #Age Distribution of TB Cases in Percentages

    V   <- (df[57:67,136:146]+df[57:67,189:199])
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
    notif_age_10<-colSums(notif_age[,])
    points(1:10,notif_age_10[]/sum(notif_age_10[])*100,pch=19,cex=1.2)

    #plot text
    mtext("Age Group",1,2.5,cex=1.2)
    mtext("Age Distribution of TB Cases (%), 2006-16",3,.8,font=2,cex=1.2)
    legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
           pt.cex=c(1,2),col=c("black","lightblue"),bg="white")
  ###############################################################################
  ################################################################################
  # Treatment Outcomes 1993-2014

  V   <- df[53:63,132:134]
  Vdisc <- V[,2]/rowSums(V)
  Vdead <- V[,3]/rowSums(V)
  tx_outcomes      <- CalibDatState$tx_outcomes[10:20,2:3]*100

  #format the plot
  plot(0,0,ylim=c(0,10),xlim=c(2002,2012),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2002:2012,Vdisc*100,lwd=2,col="red3")
  lines(2002:2012,Vdead*100,lwd=2,col="blue")

  #reported data for comparison

  points(2002:2012,tx_outcomes[,1],pch=19,cex=0.6,col="red3")
  points(2002:2012,tx_outcomes[,2],pch=19,cex=0.6,col="blue")
  lines (2002:2012,tx_outcomes[,1],lty=3,col="red3")
  lines (2002:2012,tx_outcomes[,2],lty=3,col="blue")

  #plot text

  mtext("Year",1,2.5,cex=1.2)
  mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=1.2)
  legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
         col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

  ################################################################################
  #LTBI Prevalance by Age in 2011, US born

  V  <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
  colnames(V) <- c("LTBI", "No-LTBI")

  pIGRA<-1
  v1<-V*pIGRA
  Sens_IGRA <-c(.780,.675,.712,.789,.591)
  Spec_IGRA <-c(.979,.958,.989,.985,.931)
  names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
  Va <- outer(v1[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v1[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))

  # Va <- outer(V[,1],c(0.74382,(1-0.74382)))+outer(V[,2],c((1-0.94014),0.94014))
  # Va <- outer(V[,1],c(Sens_IGRA,(1-Sens_IGRA)))+outer(V[,2],c((1-Spec_IGRA),Spec_IGRA))
  # Va<-V

  V1 <- Va[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+Va[10,]+Va[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100
  #
  # V1 <- V[-11,]; V1<-V1[-10,]
  # V1[9,] <- V1[9,]+V[10,]+V[11,]
  #
  # V2 <- rep(NA,8)
  # V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100
  #reported data for comparison
  ltbi_us_11      <- CalibDatState[["LTBI_prev_US_11_IGRA"]]
  ltbi_fb_11      <- CalibDatState[["LTBI_prev_FB_11_IGRA"]]

  ltbius<-ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100
  #format the plot
  plot(0,0,ylim=c(0,max(V2,ltbius)*1.5),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #reported data for comparison

  points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("IGRA+ LTBI in US Born Population 2011 by Age (%) [NATIONAL]",3,.8,font=2,cex=.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  #LTBI Prevalance by Age in 2011, non-US born

  V  <- cbind(t(df[62,66:76]),t(df[62,44:54]-df[62,66:76]))
  colnames(V) <- c("LTBI", "No-LTBI")

  pIGRA<-1
  v1<-V*pIGRA
  #under age 5
  v1b <- (v1[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v1[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
  #over age 5
  v1c <- outer(v1[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v1[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
  v1d<-rbind(v1b,v1c)
  colnames(V) <- c("LTBI", "No-LTBI")

  V1 <- v1d[-11,]; V1<-V1[-10,]
  V1[9,] <- V[9,]+v1d[10,]+v1d[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  # V1 <- V[-11,]; V1<-V1[-10,]
  # V1[9,] <- V[9,]+V[10,]+V[11,]
  #
  # V2 <- rep(NA,8)
  # V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #format the plot
  plot(0,0,ylim=c(0,max(V2)*2),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #reported data for comparison
  points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=.8)

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("IGRA+ LTBI in Non-US Born Population 2011 by Age (%) [NATIONAL]",3,.8,font=2,cex=.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  # Age Distribution of TB Deaths 1999-2014
  V  <- df[50:67,227:237]
  tb_deaths_dist  <- CalibDatState$tbdeaths_age_yr[,-1]/rowSums(CalibDatState$tbdeaths_age_yr[,-1])
  tb_deaths      <- as.data.frame(CalibDatState$tbdeaths[[st]][,2]*tb_deaths_dist[,])
  tb_deaths[is.na(tb_deaths)]<-0
  V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
  V3 <- colSums(V2)*1e6

  # for (i in length(tb_deaths)){ if(tb_deaths[i]=="NA"){ tb_deaths[i]<-0}}
  #format the plot
  plot(0,0,ylim=c(0,max(colSums(tb_deaths),V3)*1.5),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

  #plot the model data
  for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")

  #reported data for comparison
  points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

  #plot text
  mtext("Age Group",1,2.5,cex=1.2)
  mtext("Total TB Deaths by Age Group 1999-2016 [NATIONAL]",3,.8,font=2,cex=.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=NA,
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  # total tb deaths over time 1999-2016
  V   <- rowSums(df[57:67,227:237])*1e6
  tb_death_tot<-CalibDatState$tbdeaths[[st]][8:18,2]
  tb_death_tot[is.na(tb_death_tot)]<-0

  #format the plot
  plot(0,0,ylim=c(min(V)*.5,max(V,tb_death_tot)*1.5),xlim=c(2006,2016),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2006:2016,V,lwd=2,col="blue")

  #reported data for comparison
  points(2006:2016,CalibDatState$tbdeaths[[st]][8:18,2],pch=19,cex=0.6,col="black")
  lines (2006:2016,CalibDatState$tbdeaths[[st]][8:18,2],lty=3,col="black")

  #plot text

  mtext("Year",1,2.5,cex=1.2)
  mtext("Total TB Deaths by Year 2006-2016",3,.8,font=2,cex=1.2)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
         col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))
  ################################################################################
  graphs_pub(Par_list=Par_list)

  dev.off()
  # system(paste("open", pdfname))# code for ma

}

