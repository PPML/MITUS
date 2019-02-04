calib_demo_agna<-function(samp_i,ParMatrix){
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par <- as.numeric(ParMatrix);
    names(Par) <- names(ParMatrix)
  } else {  Par <- as.numeric(ParMatrix[samp_i,]);
  names(Par) <- colnames(ParMatrix) }  ##previously, the distribution of parameters were transformed to normal distribution in
  ##to facilitate comparisons. These first two steps convert these parameters back to their
  ##distributions
  # normal to uniform
  Par2 <- pnorm(Par,0,1)
  # uniform to true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  P <- P

  prms <-list()
  prms <- param(P)

  dd <- cSim_demo_agna(  nYrs       = 2018-1950         ,
                    nRes      = 27,
                    InitPop  = prms[["InitPop"]]    ,
                    Birthst  = prms[["Birthst"]]    ,
                    ImmNon    = prms[["ImmNon"]]    ,
                    ImmLat   = prms[["ImmLat" ]]    ,
                    ImmAct     = prms[["ImmAct"]]   ,
                    ImmFst   = prms[["ImmFst" ]]    ,
                    rEmmigFB =prms[["rEmmigFB"]]    ,
                    mubt       = prms[["mubt"]]
  )

  D <- dd$Outputs
  colnames(D) <- c(prms$ResNam[1:13],prms$ResNam[121:131], prms$ResNam[30:32])

  pdf(file=paste("MITUS_results/graphs_only_demo",Sys.time(),".pdf"), width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <-D[1:66,2]
  plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")

  lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")


  lines(1950:2015,V,lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Population: Total (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","Reported data","model"),cex=0.9,
         pch=c(15,19,NA),lwd=c(NA,1,2),lty=c(NA,3,1),col=c("grey50",1,1),bg="white",pt.cex=c(1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- t(D[65,3:13])
  V1  <- V[-3]
  V1[2] <- V1[2]+V[3]
  V2 <- V1[-4]
  V2[3] <- V2[3]+V1[4]
  V3 <- V2[-9]
  V3[8] <- V3[8]+V2[9]

  plot(1,1,ylim=c(1,80),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y" )
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,10,20,40,60,80),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i],V3[i]),border=NA,col="grey50")

  points(1:8+0.2,(CalibDat[["tot_pop16_ag_fb"]][-9,3]+CalibDat[["tot_pop16_ag_fb"]][-9,4]),pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Population by Age for FB (red) and US (blue), 2014 (mil,log-scale)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("grey30","grey80"),bg="white")


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014 all ages ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- t(D[65,3:13])

  plot(0,0,ylim=c(0.05,max(range(V))),xlim=c(0.6,11.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94", "95p"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:12-0.5,rep("",12))
  axis(2,c(0,5,10,15,20,25,30,35,40,50),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:11) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i],V[i]),border=NA,col="grey")

  mtext("Age Group",1,2.5,cex=0.9)
  mtext("Millions",2,2.5,cex=0.9)

  box()
  mtext("Total Population by Age Group 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Model"),pch=15,pt.cex=2,
         lwd=NA,col=c("gray"),bg="white")


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT EACH DECADE,  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V1 <- rowSums(D[1:66,14:24])
  plot(1,1,ylim=c(0,3.5),xlim=c(1950,2015),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2015,V1,lwd=2,col="grey50")
  points(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,pch=19,cex=0.6,col="grey50")
  lines(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,lty=3,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Mortality: Total, US, and Non-US Born",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","Reported data","model"),cex=0.9,
         pch=c(15,19,NA),lwd=c(NA,1,2),lty=c(NA,3,1),col=c("grey50",1,1),bg="white",pt.cex=c(1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- D[65,14:24]
  V1  <- V[-3]
  V1[2] <- V1[2]+V[3]
  V2 <- V1[-4]
  V2[3] <- V2[3]+V1[4]
  V3 <- V2[-9]
  V3[8] <- V3[8]+V2[9]


  plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i],V3[i]),border=NA,col="gray")
  for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality by Age, 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("black","gray"),bg="white")
  dev.off()


}
