tb_graph_immig <- function(df){

  pdf(file=paste("MITUS_results/graphs_immig",Sys.time(),".pdf"), width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))


  V<-cbind(M[1:66,31],M[1:66,32])
  plot(0,0,ylim=c(0,40),xlim=c(1950,2015),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2015,V[,2],lwd=2,col="red3")
  lines(1950:2015,V[,1],lwd=2,col="blue")
  lines(1950:2015,rowSums(V),lwd=2,col="grey50")
  points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],pch=19,cex=0.6,col="grey50")
  lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],lty=3,col="grey50")



  mtext("Year",1,2.5,cex=0.9)
  mtext("Population: Total Non-US Born (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Total Non-US","recent Non-US","longterm Non-US","Reported data","model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  V<-cbind(M[1:66,31],M[1:66,32])
  plot(0,0,ylim=c(0,20),xlim=c(1950,2015),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2015,V[,2],lwd=2,col="red3")
  lines(1950:2015,V[,1],lwd=2,col="blue")
  lines(1950:2015,rowSums(V),lwd=2,col="grey50")


  mtext("Year",1,2.5,cex=0.9)
  mtext("Population: Total Non-US Born (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Total Non-US","recent Non-US","longterm Non-US","model"),cex=0.9,
         pch=c(15,15,15,NA),lwd=c(NA,NA,NA,2),lty=c(NA,NA,NA,1),col=c("grey50",4,"red3",1),bg="white",pt.cex=c(1.8,1.8,1.8,NA))



  #immigration over time with no tb
  V<-as.data.frame(IP$ImmNon[1:65,])
  col<-rainbow(11)

  #format the plot
  plot(1,1,ylim=c(.0000001,.01),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the lines
  for (i in 1:11){
    lines(1950:2014,V[,i],lwd=3,col=col[i])
  }

  #text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Age Specific ImmNon Rates from 1950 to 2014, log-scale",3,.8,font=2,cex=0.8)
  legend("topright",colnames(V),cex=0.9,
         pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))

  ##########  ##########  ##########  ##########  ##########  ##########  ##########

  #immigration over time with latent slow tb
  V<-as.data.frame(IP$ImmLat[1:65,])
  col<-rainbow(11)

  #format the plot
  plot(1,1,ylim=c(.0000001,.001),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the lines
  for (i in 1:11){
    lines(1950:2014,V[,i],lwd=3,col=col[i])
  }

  #text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Age Specific ImmLat Rates from 1950 to 2014, log-scale",3,.8,font=2,cex=0.8)
  legend("topright",colnames(V),cex=0.9,
         pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))

  ##########  ##########  ##########  ##########  ##########  ##########  ##########

  #immigration over time with latent fast tb
  V<-as.data.frame(IP$ImmFst[1:65,])
  col<-rainbow(11)

  #format the plot
  plot(1,1,ylim=c(.00000001,.0001),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the lines
  for (i in 1:11){
    lines(1950:2014,V[,i],lwd=3,col=col[i])
  }

  #text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Age Specific ImmFst Rates from 1950 to 2014, log-scale",3,.8,font=2,cex=0.8)
  legend("topright",colnames(V),cex=0.9,
         pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))

  ##########  ##########  ##########  ##########  ##########  ##########  ##########

  #immigration over time with active tb
  V<-as.data.frame(IP$ImmAct[1:65,])
  col<-rainbow(11)

  #format the plot
  plot(1,1,ylim=c(.00000004,.000002),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the lines
  for (i in 1:11){
    lines(1950:2014,V[,i],lwd=3,col=col[i])
  }

  #text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Age Specific ImmAct Rates from 1950 to 2014, log-scale",3,.8,font=2,cex=0.8)

  legend("topright",colnames(V),cex=0.9,
         pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))

  ##########  ##########  ##########  ##########  ##########  ##########  ##########

  #death in non us born population
  V<-as.data.frame(M[1:65,266:276])
  col<-rainbow(11)

  #format the plot
  plot(0,0,ylim=c(0, max(range(V))),xlim=c(1950,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  #plot the lines
  for (i in 1:11){
    lines(1950:2014,V[,i],lwd=3,col=col[i])
  }

  #text
  mtext("Year",1,2.5,cex=0.9)
  mtext("Age Specific Mortality in NUS from 1950 to 2014 (millions)",3,.8,font=2,cex=0.8)

  legend("topright",colnames(V),cex=0.9,
         pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))



dev.off()

}
