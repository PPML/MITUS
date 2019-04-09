
inputs_graphs <-function(prms){


  pdf(file=paste("MITUS_results/graphs_inputs",Sys.time(),".pdf"), width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ANNUAL MORT RATE FOR MORT GROUPS ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  plot(prms$Birthst,xlab="",ylab="")
  box()
  mtext("Birth Params 1950-2014",3,.8,font=2,cex=0.8)


  V<-prms$mubt
  col<-rainbow(11)

  plot(1,1,ylim=c(min(V)*.5,max(V)*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific Mortality Params, log scale",3,.8,font=2,cex=1.2)
  }

  V<-prms$mubt
  col<-rainbow(11)

  plot(0,0,ylim=c(min(V)*.5,max(V)*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific Mortality Params",3,.8,font=2,cex=1.2)
  }

  #plot BgMort Inputs
  V<-Inputs$BgMort
  col<-rainbow(11)

  plot(1,1,ylim=c(min(V[2:12])*.5,max(V[,2:12])*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 2:12){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific Mortality Inputs, log scale",3,.8,font=2,cex=1.2)
  }

  V<-Inputs$BgMort
  col<-rainbow(11)

  plot(0,0,ylim=c(min(V[2:12])*.5,max(V[,2:12])*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 2:12){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific Mortality Inputs",3,.8,font=2,cex=1.2)
  }

#######IMMIGRATION PARAMETERS
  V<-prms$ImmNon
  col<-rainbow(11)

  plot(1,1,ylim=c(min(V)*.5,max(V)*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific ImmigNon Params, log scale",3,.8,font=2,cex=1.2)
  }

  V<-prms$ImmLat
  col<-rainbow(11)

  plot(1,1,ylim=c(min(V)*.5,max(V)*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific ImmigLat Params, log scale",3,.8,font=2,cex=1.2)
  }

  V<-prms$ImmFst
  col<-rainbow(11)

  plot(1,1,ylim=c(min(V)*.5,max(V)*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific ImmigFst Params, log scale",3,.8,font=2,cex=1.2)
  }

  V<-prms$ImmAct
  col<-rainbow(11)

  plot(1,1,ylim=c(min(V)*.5,max(V)*1.25),xlim=c(1,nrow(V)),xlab="",ylab="",axes=F, log = "y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
    lines(V[,i],lwd=3,col=col[i])

    mtext("Year",1,2.5,cex=1.2)
    mtext("Age Specific ImmigAct Params, log scale",3,.8,font=2,cex=1.2)
  }


  #####TB STUFF
  plot(prms$rDxt[,1],xlab="",ylab="")
  box()
  mtext("rDxt Params 1950-2014",3,.8,font=2,cex=0.8)

  plot(prms$TxQualt,xlab="",ylab="")
  box()
  mtext("TxQualt Params 1950-2014",3,.8,font=2,cex=0.8)

  plot(prms$rDeft,xlab="",ylab="")
  box()
  mtext("rDeft Params 1950-2014",3,.8,font=2,cex=0.8)

  plot(prms$pReTx,xlab="",ylab="")
  box()
  mtext("pReTx Params 1950-2014",3,.8,font=2,cex=0.8)

  plot(prms$rLtScrt,xlab="",ylab="")
  box()
  mtext("rLtScrt Params 1950-2014",3,.8,font=2,cex=0.8)

  plot(prms$EarlyTrend,xlab="",ylab="")
  box()
  mtext("EarlyTrend Params 1950-2014",3,.8,font=2,cex=0.8)

  dev.off()

}
