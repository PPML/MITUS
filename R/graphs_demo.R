#' creates plots of the demographics of the model
#' and writes these results to a .pdf file.
#' Use to check both with and without tb in the model

US_tot_mort <- read.csv(file="data/US_total_mort.csv", header = FALSE)
US_mort_age <- read.csv(file="data/US_mort_age.csv", header = TRUE)

#'Create a function to be run on a specific model run output to
#'create simple graphs of all the output for a selected year range
#'@param df dataframe of output for all years
#'@return .pdf of the graphs
#'@export

tb_graph_demo <- function(df){

  load("data/CalibDat_9-14-16.rData")
  source("R/CalibFunctionsUS_V23.r")

  pdf(file=paste("MITUS_results/graphs_demo",Sys.time(),".pdf"), width = 11, height = 8.5)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(df[1:66,30], df[1:66,31]+df[1:66,32])
  plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
  points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
  points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
  lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],lty=3,col="grey50")
  lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],lty=3,col="blue")
  lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],lty=3,col="red3")
  lines(1950:2015,V[,2],lwd=2,col="red3")
  lines(1950:2015,V[,1],lwd=2,col="blue")
  lines(1950:2015,rowSums(V),lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Population: Total, US, and Foreign Born (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,33:43]), t(df[65,44:54]))
  V1  <- V[-3,]
  V1[2,] <- V1[2,]+V[3,]
  V2 <- V1[-4,]
  V2[3,] <- V2[3,]+V1[4,]
  V3 <- V2[-9,]
  V3[8,] <- V3[8,]+V2[9,]

  plot(0,1,ylim=c(0.05,135),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y")
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

  points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")

  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Population by Age for FB (red) and US (blue), 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","Fitted model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("grey30","grey80"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP PROG RG DISTRIBUTION 2014  ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- t(df[65,20:23])

  plot(0,1,ylim=c(0.05,135),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA,log="y")
  axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:5-0.5,rep("",5))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="gray")
  mtext("Risk Group",1,2.5,cex=0.9)
  box()
  mtext("Population by TB Progression Group, 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("topright",c("Fitted model"),pch=15,pt.cex=2,,
         lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP, MORT RG DISTRIBUTION 2014 ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,305:308]), t(df[65,309:312]))

  plot(0,1,ylim=c(0.05,1e4),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA,log="y")
  axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:5-0.5,rep("",5))
  axis(2,c(0.1,1,10,100,1000,10,000),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="gray")
  mtext("Risk Group",1,2.5,cex=0.9)
  box()
  mtext("Population by Mortality Group, 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("topright",c("Fitted model"),pch=15,pt.cex=2,,
         lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP HR DISTRIBUTION 1993-2013 ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V   <- cbind(df[44:65,313],df[44:65,314],df[44:65,315],df[44:65,316] )
  V[,1] <-V[,2]/(V[,1]+ V[,2]) #high risk percentage
  V <- V[,-2]
  V[,2] <-V[,3]/(V[,2]+ V[,3]) #high risk percentage
  V <- V[,-3]

  V4 <- cbind(df[44:65,313]+df[44:65,315], df[44:65,314]+df[44:65,316])
  V4 <- V4[,2]/(V4[,1]+V4[,2])

  plot(0,0,ylim=c(0,max(range(rowSums(V)))*100),xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
  # lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)
  lines(1993:2014,V[,2]*100,lwd=2,col="red3")
  lines(1993:2014,V[,1]*100,lwd=2,col="blue")
  lines(1993:2014,V4*100,lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent Homeless Population for FB (red) and US (blue) 1993-2014 ",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(rowSums(df[1:66,255:265]), rowSums(df[1:66,266:276]))
  V1 <- df[1:66,317]
  plot(1,1,ylim=c(0,max(range(V1))),xlim=c(1950,2015),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2015,V[,2],lwd=2,col="red3")
  lines(1950:2015,V[,1],lwd=2,col="blue")
  lines(1950:2015,V1,lwd=2,col="grey50")
  points(US_tot_mort[,1],US_tot_mort[,2]/1e6,pch=19,cex=0.6,col="grey50")
  lines(US_tot_mort[,1],US_tot_mort[,2]/1e6,lty=3,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Mortality: Total, US, and Foreign Born",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind((df[65,255:265])+(df[65,266:276]))
  V1  <- V[,-3]
  V1[,2] <- V1[,2]+V[,3]
  V2 <- V1[,-4]
  V2[,3] <- V2[,3]+V1[,4]
  V3 <- V2[,-9]
  V3[,8] <- V3[,8]+V2[,9]


  plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[1,i],V3[1,i]),border=NA,col="gray")
  for(i in 1:8) points(i+.2,(US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality by Age, 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("black","gray"),bg="white")


   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,255:265]), t(df[65,266:276]))
  V1  <- V[-3,]
  V1[2,] <- V1[2,]+V[3,]
  V2 <- V1[-4,]
  V2[3,] <- V2[3,]+V1[4,]
  V3 <- V2[-9,]
  V3[8,] <- V3[8,]+V2[9,]

  plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality by Age for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Fitted model"),pch=15,pt.cex=2,,
         lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION % 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#   V  <- cbind(t(df[65,255:265]), t(df[65,266:276]))
#   V1  <- V[-3,]
#   V1[2,] <- V1[2,]+V[3,]
#   V2 <- V1[-4,]
#   V2[3,] <- V2[3,]+V1[4,]
#   V3 <- V2[-9,]
#   V3[8,] <- V3[8,]+V2[9,]
#
#   V4 <- sum(V3)
#
#   V3<- V3/V4*100
#
#   plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
#   axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
#   axis(1,1:9-0.5,rep("",9))
#   axis(2,c(0.1,1,5,10,25,32),las=2);box()
#   abline(h=axTicks(2),col="grey85")
#
#   for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
#   for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")
#
#   mtext("Age Group",1,2.5,cex=0.9)
#   box()
#   mtext("Percent Mortality by Age for FB (red) and US (blue), 2014",3,.8,font=2,cex=0.8)
#   legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),pt.cex=c(1,2),
#          lwd=NA,col=c("grey30","grey80"),bg="white")

   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT PROG DISTRIBUTION 2014 ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,277:280]),t(df[65,281:284]))

  plot(0,0,ylim=c(0.05,max(range(V))*1.001),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:5-0.5,rep("",5))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
  for(i in 1:4) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")

  mtext("Risk Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality by TB Progression Group for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topright",c("Fitted model"),pch=15,pt.cex=2,,
         lwd=NA,col="gray",bg="white")
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT MORT GROUP DISTRIBUTION 2014 ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,285:288]),t(df[65,289:292]))
  V1  <- sum(V)

  plot(0,0,ylim=c(0.05,max(range(V/V1))*100),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:5-0.5,rep("",5))
  axis(2,c(0,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,(V[i,1]/V1[1])*100,(V[i,1]/V1[1])*100),border=NA,col="lightblue")
  for(i in 1:4) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,(V[i,2]/V1[2])*100,(V[i,2]/V1[2])*100),border=NA,col="pink")

  mtext("Risk Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality % by Mortality Group for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topright",c("Fitted model"),pch=15,pt.cex=2,,
         lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT % HR DISTRIBUTION 1993-2013 ### ### ### ### ### ###
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V   <- cbind(df[44:65,313],df[44:65,315],df[44:65,314],df[44:65,316] )
  V[,1] <-V[,2]/(V[,1]+ V[,2]) #high risk percentage
  V <- V[,-2]
  V[,2] <-V[,3]/(V[,2]+ V[,3]) #high risk percentage
  V <- V[,-3]

  V4 <- cbind(df[44:65,313]+df[44:65,315], df[44:65,314]+df[44:65,316])
  V4 <- V4[,2]/(V4[,1]+V4[,2])

  plot(0,0,ylim=c(0,max(range(V4)))*100,xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
  # lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)
  lines(1993:2014,V[,2]*100,lwd=2,col="red3")
  lines(1993:2014,V[,1]*100,lwd=2,col="blue")
  lines(1993:2014,V4*100,lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent Homeless Mortality for FB (red) and US (blue) 1993-2014 ",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

dev.off()
  }
