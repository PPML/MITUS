#'
#'This script creates graphs that compare model simulations
#'versus calibration data.
#'@name graphs_calib_st
#'@param df dataframe of results
#'@return .pdf of calibration graphs
#'@export
calib_graphs_st <- function(df,loc, Par_list,pdf=TRUE, cex.size=.7){

  library(MCMCpack)
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  st<-which(StateID$USPS==loc)
  df<-as.data.frame(df)
  if (pdf==TRUE){
  pdfname<-paste("MITUS_results/",loc,"_calib_graphs_", Sys.Date(), ".pdf",sep="")
  pdf(file=pdfname, width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))
}
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(df[1:68,30], df[1:68,31]+df[1:68,32])*1e6
  #read in decade based stuff
  tot_pop<- CalibDatState[["pop_50_10"]][[st]]
  #get the FB pop from the decade
  tot_pop_yr_fb   <- tot_pop[tot_pop[,2]==0,]
  #get 2017 population
  pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
  #get 2017 fb population
  pop_ag_11_17nus <-sum(pop_ag_11_170[pop_ag_11_170[,2]==0,3][-11])
  #append the foreign born population
  tot_pop_yr_nus  <- c(tot_pop_yr_fb[,-c(1:2)], pop_ag_11_17nus)

  #get the FB pop from the decade
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
  mtext(paste("Population in: Total, US, and Non-US Born in",loc,"mil, log-scale)", sep = " "),3,.3,font=2,cex=cex.size)

  legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=cex.size,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),
         bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA), ncol=2)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,33:43]), t(df[65,44:54]))
  V3  <- V[-11,]
  V3[10,] <- V3[10,]+V[11,]
  pop_ag_11_170  <- CalibDatState[["pop_00_17"]][[st]][,c(1,2,20)]
  pop_ag_11_17us <-pop_ag_11_170[pop_ag_11_170[,2]==1,3][-11]
  pop_ag_11_17nus <-pop_ag_11_170[pop_ag_11_170[,2]==0,3][-11]
  pop_ag_11_17nus <-pop_ag_11_17nus + 100
  plot(1,1,ylim=c(min(pop_ag_11_17nus/1e6,V3)*.75,max(pop_ag_11_17us/1e6,V3)*1.5),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA,log="y")
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

  ### ### ### Population HR DISTRIBUTION 1993-2013  ### ### ### ### ### ###
  V   <- df[44:69,29]
  us_homeless<-CalibDat$homeless_pop[[st]]
  plot(0,0,ylim=c(0,max(V,us_homeless)*1.2),xlim=c(1993,2018),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(rep(2010,3),us_homeless,pch=19,cex=cex.size,col=c(4,4,4))
  lines(1993:2018,V,lwd=2,col="blue")

  mtext("Year",1,2.5,cex=cex.size)
  mtext(paste("Population Homeless in Past Yr in",loc, sep = " "),3,.3,font=2,cex=cex.size)
  legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),
         cex=cex.size,lwd=c(1,2),col=c(4,4),lty=c(NA,1),bg="white",pt.cex=0.6)

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
#   mtext("Year",1,2.5,cex=cex.size)
#   mtext("Non-US Born Population: Total, Recent, Long-Term (mil)",3,.3,font=2,cex=cex.size)
#   legend("topleft",c("Total","Recent","Long-Term","model"),cex=cex.size,
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
  # mtext("Year",1,2.5,cex=cex.size)
  #
  # mtext("Recent Immigration Non-US Born Population (<2 yrs) (000s)",3,.3,font=2,cex=cex.size)
  # legend("topleft",c("Total"),cex=cex.size,
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
#   mtext("Year",1,2.5,cex=cex.size)
#   mtext("Non-US Born Population: Total, Recent, Long-Term (mil, log-scale)",3,.3,font=2,cex=cex.size)
#   legend("bottomright",c("Total","Recent","Long-Term","model"),cex=cex.size,
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
#   mtext("Year",1,2.5,cex=cex.size)
#   mtext("Non-US Born Population: Recent, Long-Term (%)",3,.3,font=2,cex=cex.size)
#   legend("left",c("Recent","Long-Term","model"),cex=cex.size,
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
  # mtext("Age Group",1,2.5,cex=cex.size)
  # mtext("Millions",2,2.5,cex=cex.size)
  #
  # box()
  # mtext("Total Population by Age Group 2014 (mil,log-scale)",3,.3,font=2,cex=cex.size)
  # legend("topright",c("US born","Non-US Born","Reported data"),cex=cex.size,
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

  mtext("Year",1,2.5,cex=cex.size)
  mtext(paste("Total Mortality in", loc, "(in millions)", sep = " "),3,.3,font=2,cex=cex.size)
  legend("topright",c("Total","Reported data","model"),cex=cex.size,
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
  axis(2,c(0,.2,.4,.6,.3),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:10) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V2[i],V2[i]),border=NA,col="gray")
  for(i in 1:10) points(i+.2,(tda[i]),pch=19,cex=cex.size,col="black")


  mtext("Age Group",1,2.5,cex=cex.size)
  box()
  mtext(paste("Mortality Distribution by Age in", loc, ", 2016 (%)", sep = " "),3,.3,font=2,cex=cex.size)
  legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("black","gray"),bg="white",cex=cex.size)


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
  #   mtext("Year",1,2.5,cex=cex.size)
  #   mtext("Age Specific Mortality Rates from 1950 to 2014",3,.3,font=2,cex=cex.size)
  #   legend("bottomleft",colnames(V),cex=cex.size,
  #          pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))
  # }

  # graph of total diagnosed cases 5 year bands
  # by total population, US born population, and non-US born population
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  Va <- df[44:70,"NOTIF_ALL"]+df[44:70,"NOTIF_MORT_ALL"] #total population

  # tot_cases<-rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14])+rowSums(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14])
  tot_cases<-CalibDatState[["cases_yr_st"]][[st]][,2]
  #format the plot
  plot(0,0,ylim=c(0,max(Va*1e6,tot_cases)*1.25),xlim=c(1993,2019),xlab="",ylab="",axes=F)
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
  usb_5yr_cases <- as.numeric(unlist(CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==1,4:8]))
  nusb_5yr_cases<-as.numeric(unlist(CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==0,4:8]))

  plot(0,0,ylim=c(0,max(Vn2*1e6,Vu2*1e6,nusb_5yr_cases, usb_5yr_cases)*1.25),xlim=c(1995,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  #multiply raw output by 1,000 to convert from millions to hundredscale
  lines(c(1997,2002,2007,2012,2017),Vu2*1e6,lwd=3,col="white"); lines(c(1997,2002,2007,2012,2017),Vu2*1e6,lwd=2,col=4) #US born population
  lines(c(1997,2002,2007,2012,2017),Vn2*1e6,lwd=3,col="white"); lines(c(1997,2002,2007,2012,2017),Vn2*1e6,lwd=2,col=3) #non-US born population

  #reported data for comparison
  points(c(1997,2002,2007,2012,2017),usb_5yr_cases,pch=19,cex=0.3,col=4) #US born population
  lines(c(1997,2002,2007,2012,2017),usb_5yr_cases,pch=19,lty=3,col=4)
  points(c(1997,2002,2007,2012,2017),nusb_5yr_cases,pch=19,cex=0.3,col=3) #non-US born population
  lines(c(1997,2002,2007,2012,2017),nusb_5yr_cases,lty=3,col=3)

  #plot text
  mtext("Year",1,2.5,cex=cex.size)
  mtext(paste("TB Cases Identified in", loc," by Nativity, 1995-2019", sep = " "),3,.3,font=2,cex=cex.size)
  legend("topright",c("Reported data (US born)","Reported data (non-US born)",
                      "Model (US born)","Model (non-US born)"),
         pch=c(19,19,NA,NA),lwd=c(1,1,2,2),lty=c(3,3,1,1),col=c(4,3,4,3),bg="white",ncol=2,cex=cex.size,pt.cex=0.4)

  ################################################################################
  #Percent of Total Cases Non-US Born Population
  #updated for 5 year data
  V<-Vn2/(Vn2+Vu2)
  #create 5 year totals
  tot_cases<- c(sum(CalibDatState[["cases_yr_st"]][[st]][3:7,2]),
                sum(CalibDatState[["cases_yr_st"]][[st]][8:12,2]),
                sum(CalibDatState[["cases_yr_st"]][[st]][13:17,2]),
                sum(CalibDatState[["cases_yr_st"]][[st]][18:22,2]),
                sum(CalibDatState[["cases_yr_st"]][[st]][23:27,2]))
  #format the plot
  plot(0,0,ylim=c(0,min((max(V,(nusb_5yr_cases/tot_cases))*1.2*100),100)),xlim=c(1995,2020),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(c(1997,2002,2007,2012,2017),V*100,lwd=2,col=4)

  #reported data for comparison

  points(c(1997,2002,2007,2012,2017),nusb_5yr_cases/tot_cases*100,pch=19,cex=0.6)
  lines(c(1997,2002,2007,2012,2017), nusb_5yr_cases/tot_cases*100,lty=3)

  #plot text
  mtext("Year",1,2.5,cex=cex.size)
  mtext(paste("Percent of TB Cases Non-US Born in",loc,", 1995-2019", sep = " "),3,.3,font=2,cex=cex.size)
  legend("bottomright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),
         lty=c(3,1),bg="white",pt.cex=0.6,cex=cex.size)

  ################################################################################
  #Percent of Non-US Born Cases from Recent Immigrant Population
  #updated for 5 year data
  #check this plot and otis
  V <- cbind(df[46:70,"NOTIF_F1"]+df[46:70,"NOTIF_MORT_F1"],df[46:70,"NOTIF_F2"]+df[46:70,"NOTIF_MORT_F2"])
  #create a five year band for this data
  V0<-rep(0,5)
  V0[1]<-sum(V[1:5,1])/sum(V[1:5,]);V0[2]<-sum(V[6:10,1])/sum(V[6:10,]);
  V0[3]<-sum(V[11:15,1])/sum(V[11:15,]);V0[4]<-sum(V[16:20,1])/sum(V[16:20,]);V0[5]<-sum(V[21:25,1])/sum(V[21:25,]);
  V0<-V0*100
  #reported data for comparison

  notif_rec<-CalibDatState[["rt_fb_cases_sm"]][which(CalibDatState[["rt_fb_cases_sm"]][,1]==stateID[st,1]),9]
  notif_time<-c(1997,2002,2007,2012,2017)
  notif_fb_rec<-cbind(notif_rec*100, notif_time)
  #format the plot
  # plot(0,0,ylim=c(min(V,notif_fb_rec[,1])*.25,min(max(V,notif_fb_rec[,1])*1.5,100)),xlim=c(2007,2017),xlab="",ylab="",axes=F)
  plot(0,0,ylim=c(min(V,notif_fb_rec[,1])*.25,min(max(V,notif_fb_rec[,1])*1.5,100)),xlim=c(1995,2019),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(c(1997,2002,2007,2012,2017),V0,lwd=2,col=4)

  #reported data for comparison
  points(notif_fb_rec[,2],notif_fb_rec[,1],pch=19,cex=0.6)
  lines(notif_fb_rec[,2],notif_fb_rec[,1],lty=3)

  #plot text
  mtext("Year",1,2.5,cex=cex.size)
  mtext(paste("Percent of Non-US Born Cases Arrived in Past 2 Yrs in",loc, sep = " "),3,.3,font=2,cex=cex.size)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),
         lty=c(3,1),bg="white",pt.cex=0.6,cex=cex.size)

  ################################################################################
  #Age distribution of Cases
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs
  #updated for 5 year data
  V   <- (df[46:70,136:146]+df[46:70,189:199])*1e6
  V2  <- V[,-11]
  V2[,10] <- V2[,10]+V[,11]
  #create 5 year bands
  V3<-matrix(0,5,ncol(V2))
  V3[1,]<-colSums(V2[1:5,]);V3[2,]<-colSums(V2[6:10,]);V3[3,]<-colSums(V2[11:15,]); V3[4,]<-colSums(V2[16:20,]); V3[5,]<-colSums(V2[21:25,])

  notif_age_us     <- matrix(unlist(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][1:5,5:14]),5,10)
  notif_age_nus    <- matrix(unlist(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][6:10,5:14]),5,10)
  notif_age<- notif_age_us+notif_age_nus
  #format the plot
  cls <- colorRampPalette(c("blue", "red"))( 4 )
  plot(0,0,ylim=c(min(V3,notif_age[!is.na(notif_age)])*.5,max(V3,notif_age[!is.na(notif_age)])*2),xlim=c(1995,2019),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(c(1997,2002,2007,2012,2017),rowSums(V3[,1:3]),lwd=2,col=cls[1])    #0-26 yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(V3[,4:5]),lwd=2,col=cls[2])    #25-44 yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(V3[,6:7]),lwd=2,col=cls[3])    #45-64 yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(V3[,8:10]),lwd=2,col=cls[4])   #65+ yrs

  #reported data for comparison

  points(c(1997,2002,2007,2012,2017),rowSums(notif_age[,1:3]),pch=19,cex=0.6,col=cls[1]) #0-26 yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(notif_age[,1:3]),col=cls[1],lty=3)
  points(c(1997,2002,2007,2012,2017),rowSums(notif_age[,4:5]),pch=19,cex=0.6,col=cls[2]) #25-44 yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(notif_age[,4:5]),col=cls[2],lty=3)
  points(c(1997,2002,2007,2012,2017),rowSums(notif_age[,6:7]),pch=19,cex=0.6,col=cls[3]) #45-64 yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(notif_age[,6:7]),col=cls[3],lty=3)
  points(c(1997,2002,2007,2012,2017),rowSums(notif_age[,8:10]),pch=19,cex=0.6,col=cls[4]) #65+ yrs
  lines(c(1997,2002,2007,2012,2017),rowSums(notif_age[,8:10]),col=cls[4],lty=3)

  #plot text
  mtext(paste("TB Cases By Age in",loc,", 1995-2019", sep = " "),3,.3,font=2,cex=cex.size)
  mtext("Year",1,2.5,cex=cex.size)

  legend("topright",c("0-26 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
         lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
         pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA),cex=cex.size,ncol=2)

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
    notif_age_10<-colSums(notif_age[,])
    points(1:10,notif_age_10[]/sum(notif_age_10[])*100,pch=19,cex=cex.size)

    #plot text
    mtext("Age Group",1,2.5,cex=cex.size)
    mtext(paste("Age Distribution of TB Cases (%) in",loc,", 1995-19", sep = " "),3,.3,font=2,cex=cex.size)
    legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
           pt.cex=c(1,2),col=c("black","lightblue"),bg="white",cex=cex.size)

    ################################################################################
    #Average Age of TB Cases
    age_case   <- df[,136:146]+df[,189:199]
    ages<-c(2,9.5,19.5,29.5,39.5,49.5,
            59.5,69.5,79.5,89.5,99.5)
    avg_age<-rep(0,nrow(age_case))
    for (i in 1:nrow(age_case)){
      if (sum(age_case[i,]) == 0){
        avg_age[i] = 0
      } else{
        avg_age[i]<-sum(age_case[i,]*ages)/sum(age_case[i,])
      }
    }

    #format the plot
    plot(0,0,ylim=c(0,max(avg_age[2:length(avg_age)])+5),xlim=c(1951,2019),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")

    #plot the model data
    lines(1951:2019,avg_age[2:length(avg_age)],lwd=2,col="blue")    #0-24 yrs

    #plot text
    mtext(paste("Average Age of Notified TB Case in",loc,", 1950-2019", sep = " "),3,.3,font=2,cex=cex.size)
    mtext("Year",1,2.5,cex=cex.size)

    ################################################################################
    ### ### ### CASES HR DISTRIBUTION 1993-2013  ### ### ### ### ### ###

    X   <- (df[46:70,151]+df[46:70,204])
    Xa  <- rbind(sum(X[1:5]),sum(X[6:10]), sum(X[11:15]),
                 sum(X[16:20]), sum(X[21:25]))
    W<-(df[46:70,151]+df[46:70,150]+df[46:70,204]+df[46:70,203])
    Wa  <- rbind(sum(W[1:5]),sum(W[6:10]), sum(W[11:15]),
                 sum(W[16:20]), sum(W[21:25]))

    V<-Xa/Wa

    notif_hr<-CalibDatState[["hr_cases_sm"]][which(CalibDatState[["hr_cases_sm"]][,1]==stateID[st,1]),7]
    notif_hr<-cbind(notif_time,notif_hr)
    plot(0,0,ylim=c(0,max(V,notif_hr[,2])*110),xlim=c(1995,2020),xlab="",ylab="",axes=F)
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")
    points(notif_time,notif_hr[,2]*100,pch=19,cex=0.6)
    lines(notif_time,notif_hr[,2]*100,lty=3)
    lines(notif_time,V*100,lwd=2,col=4)
    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("Percent of TB Cases Homeless in Past Yr in", loc, sep = " "),3,.3,font=2,cex=cex.size)
    legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),
           lty=c(3,1),bg="white",pt.cex=0.6,cex=cex.size)
  ###############################################################################
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
    points(1,rct_trans_dist,pch=19,cex=cex.size)
    text(1:16,Vall,format(round(Vall,2),nsmall=2),cex=cex.size*1.25,pos=3)
    legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
           pt.cex=c(1,2),col=c("black","forestgreen"),bg="white",cex=cex.size)
  ################################################################################
    ### ### ### LTBI INITIATIONS 1993-2011 Distribution ### ### ### ### ### ###
    v13  <- df[43:65,153:154]/df[43:65,152]
    TLTBI_dist<-CalibDat$TLTBI_dist[1:2]
    tltbi_vol<-CalibDat$TLTBI_volume[[st]]
    plot(1,1,ylim=c(0.001,1)*100,xlim=c(1992,2015),xlab="",ylab="",axes=F,log="y")
    axis(1);axis(2,las=2);box()
    abline(h=axTicks(2),col="grey85")
    points(rep(2002,2),TLTBI_dist*100,pch=19,cex=cex.size,col=c("red3",4,6))
    for(i in 1:2) lines(1992:2014,v13[,i]*100,lwd=2,col=c("red3",4,6)[i])
    # lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

    mtext("Year",1,2.5,cex=cex.size)
    mtext(paste("IPT Treatment Initiations By Risk Group (%) in", loc, sep = " "),3,.3,font=2,cex=cex.size)
    legend("bottomleft",c("Foreign-born","Homeless","Reported data","Fitted model"),
           pch=c(15,15,19,NA),lwd=c(NA,NA,NA,2),col=c("red3",4,1,1),bg="white",pt.cex=c(1.8,1.8,0.8,NA),cex=cex.size)
    ################################################################################
  # LTBI Outcomes 1993-2014
  V   <- df[56:66,132:134]
  Vdisc <- V[,2]/rowSums(V)
  Vdead <- V[,3]/rowSums(V)
  tx_outcomes      <- CalibDatState$tx_outcomes[13:23,2:3]*100

  #format the plot
  plot(0,0,ylim=c(0,10),xlim=c(2005,2015),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the model data
  lines(2005:2015,Vdisc*100,lwd=2,col="red3")
  lines(2005:2015,Vdead*100,lwd=2,col="blue")

  #reported data for comparison

  points(2005:2015,tx_outcomes[,1],pch=19,cex=0.6,col="red3")
  points(2005:2015,tx_outcomes[,2],pch=19,cex=0.6,col="blue")
  lines (2005:2015,tx_outcomes[,1],lty=3,col="red3")
  lines (2005:2015,tx_outcomes[,2],lty=3,col="blue")

  #plot text

  mtext("Year",1,2.5,cex=cex.size)
  mtext(paste("Treatment Outcomes: Discontinued and Died (%) in", loc, sep = " "),3,.3,font=2,cex=cex.size)
  legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
         col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA),cex=cex.size)

  ################################################################################
  #LTBI Prevalance by Age in 2011, US born

  V  <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
  colnames(V) <- c("LTBI", "No-LTBI")

  pIGRA<-1
  v1<-V*pIGRA
  Sens_IGRA <-c(.780,.780,.712,.789,.789)
  Spec_IGRA <-c(.979,.979,.989,.985,.985)
  names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
  Va <- outer(v1[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v1[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))

  # Va <- outer(V[,1],c(0.74382,(1-0.74382)))+outer(V[,2],c((1-0.94014),0.94014))
  # Va <- outer(V[,1],c(Sens_IGRA,(1-Sens_IGRA)))+outer(V[,2],c((1-Spec_IGRA),Spec_IGRA))
  # Va<-V

  V1 <- Va[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+Va[10,]+Va[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  # Vtot<-sum(V1[2:9,1])/sum(V1[2:9])*100
  #
  # V2<-c(V2,Vtot)
  #
  # V1 <- V[-11,]; V1<-V1[-10,]
  # V1[9,] <- V1[9,]+V[10,]+V[11,]
  #
  # V2 <- rep(NA,8)
  # V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100
  #reported data for comparison
  ltbi_us_11      <- CalibDatState[["LTBI_prev_US_11_IGRA"]]
  # #read in Maryam estimates
  # ltbi_us_11      <- rbind(ltbi_us_11)
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

  points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=cex.size)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=cex.size)

  #plot text
  mtext("Age Group",1,2.5,cex=cex.size)
  mtext(paste("IGRA+ LTBI in US Born Population 2011 by Age (%) in",loc,"[NATIONAL]", sep = " "),3,.3,font=2,cex=cex.size)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white",cex=cex.size)

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
  points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=cex.size)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=cex.size)

  #plot text
  mtext("Age Group",1,2.5,cex=cex.size)
  mtext(paste("IGRA+ LTBI in Non-US Born Population 2011 by Age (%) in",loc,"[NATIONAL]", sep = " "),3,.3,font=2,cex=cex.size)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white",cex=cex.size)

  ################################################################################
  # Age Distribution of TB Deaths 1999-2014
  V  <- df[50:69,227:237]
  V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
  V3 <- colSums(V2)*1e6
  V3 <- V3/sum(V3)*100

  tb_deaths_dist  <- colSums(CalibDatState$tbdeaths_age_yr[,-1])/sum(CalibDatState$tbdeaths_age_yr[,-1])*100
  # tb_deaths      <- as.data.frame(as.numeric(CalibDatState$tbdeaths[[st]][,3])*tb_deaths_dist[,])
  # tb_deaths[is.na(tb_deaths)]<-0

  # for (i in length(tb_deaths)){ if(tb_deaths[i]=="NA"){ tb_deaths[i]<-0}}
  #format the plot
  plot(0,0,ylim=c(0,max(tb_deaths_dist,V3)*1.5),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

  #plot the model data
  for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")

  #reported data for comparison
  points(1:10,tb_deaths_dist,pch=19,cex=cex.size,col="black")

  #plot text
  mtext("Age Group",1,2.5,cex=cex.size)
  mtext(paste("Total TB Deaths by Age Group 1999-2018 in",loc," (%) [NATIONAL]", sep = " "),3,.3,font=2,cex=cex.size)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=NA,
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white",cex=cex.size)

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
  ################################################################################
  graphs_pub(Par_list=Par_list)
if (pdf==TRUE) {dev.off()}
  # system(paste("open", pdfname))# code for ma

}

