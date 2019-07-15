
library(MCMCpack)

futurePC_graphs <- function(dfBc,dfP1,dfP2,dfP3,endyr){
  dfBc<-as.data.frame(dfBc)
  dfP1<-as.data.frame(dfP1)
  dfP2<-as.data.frame(dfP2)
  dfP3<-as.data.frame(dfP3)

  # pdf(file=paste("MITUS_results/future_graphs",Sys.time(),".pdf"), width = 11, height = 8.5)
  # par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  #' graph of total diagnosed cases
  #' by total population, US born population, and non-US born population
  V0Bc <- dfBc[4:(endyr-1949),"NOTIF_ALL"]+dfBc[4:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1Bc <- dfBc[44:(endyr-1949),"NOTIF_US"]+dfBc[44:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2Bc <- dfBc[44:(endyr-1949),"NOTIF_F1"]+dfBc[44:(endyr-1949),"NOTIF_F2"]+dfBc[44:(endyr-1949),"NOTIF_MORT_F1"]+dfBc[44:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  V0P1 <- dfP1[4:(endyr-1949),"NOTIF_ALL"]+dfP1[4:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1P1 <- dfP1[44:(endyr-1949),"NOTIF_US"]+dfP1[44:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2P1 <- dfP1[44:(endyr-1949),"NOTIF_F1"]+dfP1[44:(endyr-1949),"NOTIF_F2"]+dfP1[44:(endyr-1949),"NOTIF_MORT_F1"]+dfP1[44:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  V0P2 <- dfP2[4:(endyr-1949),"NOTIF_ALL"]+dfP2[4:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1P2 <- dfP2[44:(endyr-1949),"NOTIF_US"]+dfP2[44:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2P2 <- dfP2[44:(endyr-1949),"NOTIF_F1"]+dfP2[44:(endyr-1949),"NOTIF_F2"]+dfP2[44:(endyr-1949),"NOTIF_MORT_F1"]+dfP2[44:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  V0P3 <- dfP3[4:(endyr-1949),"NOTIF_ALL"]+dfP3[4:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1P3 <- dfP3[44:(endyr-1949),"NOTIF_US"]+dfP3[44:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2P3 <- dfP3[44:(endyr-1949),"NOTIF_F1"]+dfP3[44:(endyr-1949),"NOTIF_F2"]+dfP3[44:(endyr-1949),"NOTIF_MORT_F1"]+dfP3[44:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  #'format the plot
  plot(1,1,ylim=c(min(range(V1),range(V2)),max(range(V0)))*1e3,xlim=c(1954,endyr),xlab="",ylab="",axes=F)

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

  points(1993:2016,notif_fb[,2],pch=19,cex=0.3,col=4) #US born population
  lines(1993:2016,notif_fb[,2],pch=19,lty=3,col=4)

  points(1993:2016,notif_fb[,1],pch=19,cex=0.3,col=3) #non-US born population
  lines(1993:2016,notif_fb[,1],lty=3,col=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Total TB Cases Identified (000s)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                      "Model (all)","Model (US born)","Model (foreign born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  #' graph of total diagnosed cases
  #' by total population, US born population, and non-US born population
  V0Bc <- dfBc[51:(endyr-1949),"NOTIF_ALL"]+dfBc[51:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1Bc <- dfBc[51:(endyr-1949),"NOTIF_US"]+dfBc[51:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2Bc <- dfBc[51:(endyr-1949),"NOTIF_F1"]+dfBc[51:(endyr-1949),"NOTIF_F2"]+dfBc[51:(endyr-1949),"NOTIF_MORT_F1"]+dfBc[51:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  V0P1 <- dfP1[51:(endyr-1949),"NOTIF_ALL"]+dfP1[51:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1P1 <- dfP1[51:(endyr-1949),"NOTIF_US"]+dfP1[51:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2P1 <- dfP1[51:(endyr-1949),"NOTIF_F1"]+dfP1[51:(endyr-1949),"NOTIF_F2"]+dfP1[51:(endyr-1949),"NOTIF_MORT_F1"]+dfP1[51:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  V0P2 <- dfP2[51:(endyr-1949),"NOTIF_ALL"]+dfP2[51:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1P2 <- dfP2[51:(endyr-1949),"NOTIF_US"]+dfP2[51:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2P2 <- dfP2[51:(endyr-1949),"NOTIF_F1"]+dfP2[51:(endyr-1949),"NOTIF_F2"]+dfP2[51:(endyr-1949),"NOTIF_MORT_F1"]+dfP2[51:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  V0P3 <- dfP3[51:(endyr-1949),"NOTIF_ALL"]+dfP3[51:(endyr-1949),"NOTIF_MORT_ALL"] #total population
  V1P3 <- dfP3[51:(endyr-1949),"NOTIF_US"]+dfP3[51:(endyr-1949),"NOTIF_MORT_US"]   #US born population
  V2P3 <- dfP3[51:(endyr-1949),"NOTIF_F1"]+dfP3[51:(endyr-1949),"NOTIF_F2"]+dfP3[51:(endyr-1949),"NOTIF_MORT_F1"]+dfP3[51:(endyr-1949),"NOTIF_MORT_F2"]   #non-US born population

  #'format the plot
  plot(1,1,ylim=c(min(range(V1),range(V2)),max(range(V0)))*1e3,xlim=c(2000,endyr),xlab="",ylab="",axes=F)
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

  notif_fb      <- cbind(CalibDat[["fb_cases"]][,2],1-CalibDat[["fb_cases"]][,2])*CalibDat[["fb_cases"]][,3]
  notif_fb <-notif_fb/1000

  points(2000:2016,notif_fb[8:24,2],pch=19,cex=0.3,col=4) #US born population
  lines(2000:2016,notif_fb[8:24,2],pch=19,lty=3,col=4)

  points(2000:2016,notif_fb[8:24,1],pch=19,cex=0.3,col=3) #non-US born population
  lines(2000:2016,notif_fb[8:24,1],lty=3,col=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Total TB Cases Identified (000s)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                      "Model (all)","Model (US born)","Model (foreign born)"),
         pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

  ################################################################################
  #'Percent of Total Cases Non-US Born Population

  VBc <- cbind(dfBc[44:(endyr-1949),"NOTIF_US"]+dfBc[44:(endyr-1949),"NOTIF_MORT_US"], #US born population
             dfBc[44:(endyr-1949),"NOTIF_F1"]+dfBc[44:(endyr-1949),"NOTIF_F2"]+  #non-US born population
               dfBc[44:(endyr-1949),"NOTIF_MORT_F1"]+dfBc[44:(endyr-1949),"NOTIF_MORT_F2"])
  VBc <- VBc[,2]/rowSums(VBc)

  VP1 <- cbind(dfP1[44:(endyr-1949),"NOTIF_US"]+dfP1[44:(endyr-1949),"NOTIF_MORT_US"], #US born population
               dfP1[44:(endyr-1949),"NOTIF_F1"]+dfP1[44:(endyr-1949),"NOTIF_F2"]+  #non-US born population
                 dfP1[44:(endyr-1949),"NOTIF_MORT_F1"]+dfP1[44:(endyr-1949),"NOTIF_MORT_F2"])
  VP1 <- VP1[,2]/rowSums(VP1)

  VP2 <- cbind(dfP2[44:(endyr-1949),"NOTIF_US"]+dfP2[44:(endyr-1949),"NOTIF_MORT_US"], #US born population
               dfP2[44:(endyr-1949),"NOTIF_F1"]+dfP2[44:(endyr-1949),"NOTIF_F2"]+  #non-US born population
                 dfP2[44:(endyr-1949),"NOTIF_MORT_F1"]+dfP2[44:(endyr-1949),"NOTIF_MORT_F2"])
  VP2 <- VP2[,2]/rowSums(VP2)

  VP3 <- cbind(dfP3[44:(endyr-1949),"NOTIF_US"]+dfP3[44:(endyr-1949),"NOTIF_MORT_US"], #US born population
               dfP3[44:(endyr-1949),"NOTIF_F1"]+dfP3[44:(endyr-1949),"NOTIF_F2"]+  #non-US born population
                 dfP3[44:(endyr-1949),"NOTIF_MORT_F1"]+dfP3[44:(endyr-1949),"NOTIF_MORT_F2"])
  VP3 <- VP3[,2]/rowSums(VP3)

  #'format the plot
  plot(0,0,ylim=c(2.5,97.5),xlim=c(2000,endyr),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(1993:endyr,VBc*100,lwd=2,col="grey50")
  lines(1993:endyr,VP1*100,lwd=2,col="red")
  lines(1993:endyr,VP2*100,lwd=2,col="green")
  lines(1993:endyr,VP3*100,lwd=2,col="blue")


  #'reported data for comparison
  points(1993:2016,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
  lines(1993:2016,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent of TB Cases Non-US-Born",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

  ################################################################################
  #'Percent of Non-US Born Cases from Recent Immigrant Population

  V <- cbind(dfBc[44:(endyr-1949),"NOTIF_F1"]+dfBc[44:(endyr-1949),"NOTIF_MORT_F1"],dfBc[44:(endyr-1949),"NOTIF_F2"]+dfBc[44:(endyr-1949),"NOTIF_MORT_F2"])
  V <- V[,1]/rowSums(V)

  #'format the plot
  plot(0,0,ylim=c(0,60),xlim=c(1993,endyr),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(1993:endyr,V*100,lwd=2,col=4)

  #'reported data for comparison
  notif_fb_rec   <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]

  points(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,pch=19,cex=0.6)
  lines(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,lty=3)

  #'plot text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)


  ################################################################################
  #' Treatment Outcomes 1993-2014

  V   <- dfBc[44:(endyr-1949),132:134]
  Vdisc <- V[,2]/rowSums(V)
  Vdead <- V[,3]/rowSums(V)

  #'format the plot
  plot(0,0,ylim=c(0,15),xlim=c(1993,endyr),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #'plot the model data
  lines(1993:endyr,Vdisc*100,lwd=2,col="red3")
  lines(1993:endyr,Vdead*100,lwd=2,col="blue")

  #'reported data for comparison
  tx_outcomes      <- cbind(1-rowSums(CalibDat[["tx_outcomes"]][,2:3]),CalibDat[["tx_outcomes"]][,2],CalibDat[["tx_outcomes"]][,3])*CalibDat[["tx_outcomes"]][,4]

  points(1993:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
  points(1993:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
  lines (1993:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
  lines (1993:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")

  #'plot text

  mtext("Year",1,2.5,cex=0.9)
  mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=0.8)
  legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
         col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))




  # dev.off() # code for ma

}


