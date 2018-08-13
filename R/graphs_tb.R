#'creates plots specific to the TB dynamics of the model
#'and writes these results to a .pdf file.
#'Use to check outputs of the model; not for distribution


#'Create a function to be run on a specific mdoel run output to
#'create simple graphs.
#'@param df dataframe of output for all years
#'@return .pdf of all graphs
#'@export

# pop_dist aross risk group and ltbi

tbdyn_graphs <-function(df){


  data("CalibDat_2018-07-12", package='MITUS')

  pdf(file=paste("MITUS_results/graphs_tbdyn",Sys.time(),".pdf"), width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  #' graph of total diagnosed cases
  #' by total population, US born population, and non-US born population
  V0 <- df[4:67,"NOTIF_ALL"]+df[4:67,"NOTIF_MORT_ALL"] #total population
  V1 <- df[44:67,"NOTIF_US"]+df[44:67,"NOTIF_MORT_US"]   #US born population
  V2 <- df[44:67,"NOTIF_F1"]+df[44:67,"NOTIF_F2"]+df[44:67,"NOTIF_MORT_F1"]+df[44:67,"NOTIF_MORT_F2"]   #non-US born population

  #'format the plot
  plot(1,1,ylim=c(0.01,10000),xlim=c(1954,2015),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  #'multiply raw output by 1,000 to convert from millions to thousands
  lines(1953:2016,V0*1e3,lwd=3,col="white"); lines(1953:2016,V0*1e3,lwd=2,col=1) #total population
  lines(1993:2016,V1*1e3,lwd=3,col="white"); lines(1993:2016,V1*1e3,lwd=2,col=4) #US born population
  lines(1993:2016,V2*1e3,lwd=3,col="white"); lines(1993:2016,V2*1e3,lwd=2,col=3) #non-US born population

  #'reported data for comparison
  points(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2],pch=19,cex=0.3) #total population
  lines(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2],lty=3,col=1)

  notif_fb      <- cbind(CalibDat[["fb_cases"]][,2],1-CalibDat[["fb_cases"]][,2])*CalibDat[["fb_cases"]][,3]

  points(1993:2016,notif_fb[,2],pch=19,cex=0.3,col=4) #US born population
  lines(1993:2016,notif_fb[,2],pch=19,lty=3,col=4)

  points(1993:2016,notif_fb[,1],pch=19,cex=0.3,col=3) #non-US born population
  lines(1993:2016,notif_fb[,1],lty=3,col=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Total TB Cases Identified (000s), 1953-2014",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                      "Model (all)","Model (US born)","Model (foreign born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

################################################################################

  #'Percent of Total Cases Non-US Born Population

  V <- cbind(df[44:67,"NOTIF_US"]+df[44:67,"NOTIF_MORT_US"], #US born population
             df[44:67,"NOTIF_F1"]+df[44:67,"NOTIF_F2"]+  #non-US born population
             df[44:67,"NOTIF_MORT_F1"]+df[44:67,"NOTIF_MORT_F2"])
  V <- V[,2]/rowSums(V)

  #'format the plot
  plot(0,0,ylim=c(2.5,97.5),xlim=c(2000,2016),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(1993:2016,V*100,lwd=2,col=4)

  #'reported data for comparison
  points(1993:2016,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
  lines(1993:2016,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent of TB Cases Non-US-Born, 2000-14",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
  #'Percent of Non-US Born Cases from Recent Immigrant Population

  V <- cbind(df[44:65,"NOTIF_F1"]+df[44:65,"NOTIF_MORT_F1"],df[44:65,"NOTIF_F2"]+df[44:65,"NOTIF_MORT_F2"])
  V <- V[,1]/rowSums(V)

  #'format the plot
  plot(0,0,ylim=c(0,60),xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(1993:2014,V*100,lwd=2,col=4)

  #'reported data for comparison
  notif_fb_rec   <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]

  points(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,pch=19,cex=0.6)
  lines(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,lty=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
  #'Age distribution of Cases
  #'0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

  V   <- (df[51:67,136:146]+df[51:67,189:199])
  V2  <- V[,-11]
  V2[,10] <- V2[,10]+V[,11]

  #'format the plot
  cls <- colorRampPalette(c("blue", "red"))( 4 )
  plot(0,0,ylim=c(0,40),xlim=c(2000,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(2000:2016,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])    #0-24 yrs
  lines(2000:2016,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])    #25-44 yrs
  lines(2000:2016,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])    #45-64 yrs
  lines(2000:2016,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])   #65+ yrs

  #'reported data for comparison
  notif_age     <- CalibDat[["age_cases"]][,-c(1,12)]*CalibDat[["age_cases"]][,12]

  points(1993:2016,rowSums(notif_age[,1:3])*1e3,pch=19,cex=0.6,col=cls[1]) #0-24 yrs
  lines(1993:2016,rowSums(notif_age[,1:3])*1e3,col=cls[1],lty=3)
  points(1993:2016,rowSums(notif_age[,4:5])*1e3,pch=19,cex=0.6,col=cls[2]) #25-44 yrs
  lines(1993:2016,rowSums(notif_age[,4:5])*1e3,col=cls[2],lty=3)
  points(1993:2016,rowSums(notif_age[,6:7])*1e3,pch=19,cex=0.6,col=cls[3]) #45-64 yrs
  lines(1993:2016,rowSums(notif_age[,6:7])*1e3,col=cls[3],lty=3)
  points(1993:2016,rowSums(notif_age[,8:10])*1e3,pch=19,cex=0.6,col=cls[4]) #65+ yrs
  lines(1993:2016,rowSums(notif_age[,8:10])*1e3,col=cls[4],lty=3)

  #'plot text
  mtext("TB Cases By Age (000s), 2000-16",3,.8,font=2,cex=0.8)
  mtext("Year",1,2.5,cex=0.9)

  legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
         lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
         pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

################################################################################
  #'Age Distribution of TB Cases in Percentages
  #'0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs
#'
#'   V   <- (df[51:65,136:146]+df[51:65,189:199])
#'   V2  <- V[,-11]
#'   V2[,10] <- V2[,10]+V[,11]
#'
#'   #'format the plot
#'   plot(0,0,ylim=c(0,max(range(V2))),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
#'   axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
#'        tick=F,cex.axis=0.6)
#'   axis(1,1:11-0.5,rep("",11))
#'   axis(2,las=2);box()
#'   abline(h=axTicks(2),col="grey85")
#'
#'   #'plot the model data
#'   for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")
#'
#'   #'reported data for comparison
#'   points(1:10,colSums(notif_age[7:21,])/sum(notif_age[7:21,])*100,pch=19,cex=1.2)
#'
#'   #'plot text
#'   mtext("Age Group",1,2.5,cex=0.9)
#'   mtext("Age Distribution of TB Cases (%), 2000-14",3,.8,font=2,cex=0.8)
#'   legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
#'          pt.cex=c(1,2),col=c("black","lightblue"),bg="white")
################################################################################
  #' Distribution of Cases across the TB Progression Categories


  #'format the plot
  #'plot the model data
  #'reported data for comparison
  #'plot text

################################################################################
  #' Treatment Outcomes 1993-2014

  V   <- df[44:65,132:134]
  Vdisc <- V[,2]/rowSums(V)
  Vdead <- V[,3]/rowSums(V)

  #'format the plot
  plot(0,0,ylim=c(0,15),xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(1993:2014,Vdisc*100,lwd=2,col="red3")
  lines(1993:2014,Vdead*100,lwd=2,col="blue")

  #'reported data for comparison
  tx_outcomes      <- cbind(1-rowSums(CalibDat[["tx_outcomes"]][,2:3]),CalibDat[["tx_outcomes"]][,2],CalibDat[["tx_outcomes"]][,3])*CalibDat[["tx_outcomes"]][,4]

  points(1993:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
  points(1993:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
  lines(1993:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
  lines(1993:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")

  #'plot text

  mtext("Year",1,2.5,cex=0.9)
  mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=0.8)
  legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
         col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

################################################################################
  #'LTBI Prevalance by Age in 2011, US born

  V  <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
  colnames(V) <- c("LTBI", "No-LTBI")

  V1 <- V[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+V[10,]+V[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #'format the plot
  plot(0,0,ylim=c(0,max(range(V2))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #'reported data for comparison
  ltbi_us_11      <- CalibDat[["LTBI_prev_US_11_IGRA"]]
  ltbi_fb_11      <- CalibDat[["LTBI_prev_FB_11_IGRA"]]

  points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)

  #'plot text
  mtext("Age Group",1,2.5,cex=0.9)
  mtext("LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  #'LTBI Prevalance by Age in 2011, non-US born

  V  <- cbind(t(df[62,66:76]),t(df[62,44:54]-df[62,66:76]))

  colnames(V) <- c("LTBI", "No-LTBI")

  V1 <- V[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+V[10,]+V[11,]

  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

  #'format the plot
  plot(0,0,ylim=c(0,max(range(V2))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
  axis(1,1:8-0.5,rep("",8))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

  #'reported data for comparison
  points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
  for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=1.2)

  #'plot text
  mtext("Age Group",1,2.5,cex=0.9)
  mtext("LTBI in Non-US Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  ################################################################################
  #' Age Distribution of TB Deaths 1999-2013

  V  <- df[50:65,88:98]+df[50:65,99:109]
  V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
  V3 <- colSums(V2)*1e6

  #'format the plot
  plot(0,0,ylim=c(0,max(range(V3))+500),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

  #'plot the model data
  for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")

  #'reported data for comparison
  tb_deaths      <- CalibDat[["tb_deaths"]][,-1]

  points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

  #'plot text
  mtext("Age Group",1,2.5,cex=0.9)
  mtext("Total TB Deaths by Age Group",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=NA,
         pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

  dev.off()

}

