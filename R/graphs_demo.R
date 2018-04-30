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
  legend("bottomright",c("Total","US born","Foreign born","Reported data","model"),cex=0.9,
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

  plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0,10,25,40,50,60,75,78 ),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

  points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")

  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Population by Age for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("grey30","grey80"),bg="white")


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014 all ages ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,33:43])+t(df[65,44:54]))

  plot(0,0,ylim=c(0.05,max(range(V))),xlim=c(0.6,11.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94", "95p"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:12-0.5,rep("",12))
  axis(2,c(0,5,10,15,20,25,30,35,40,50),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:11) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="grey")

  mtext("Age Group",1,2.5,cex=0.9)
  mtext("Millions",2,2.5,cex=0.9)

  box()
  mtext("Total Population by Age Group 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Model"),pch=15,pt.cex=2,
         lwd=NA,col=c("gray"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP PROG RG DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- t(df[65,20:23])

  plot(0,0,ylim=c(0,max(range(V))),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:5-0.5,rep("",5))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="gray")
  mtext("Risk Group",1,2.5,cex=0.9)
  box()
  mtext("Population by TB Progression Group, 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topright",c("model"),pch=15,pt.cex=2,
         lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP, MORT RG DISTRIBUTION 2014 ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,309:312]), t(df[65,313:316]))

  plot(0,0,ylim=c(0,max(range(V))),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:5-0.5,rep("",5))
  axis(2,c(0,100,500,1000,1600),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="gray")
  mtext("Risk Group",1,2.5,cex=0.9)
  box()
  mtext("Population by Mortality Group, 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topright",c("model"),pch=15,pt.cex=2,lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP JOINT DISTRIBUTION 2014 ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  mTrsp <- function(cl,a)  { apply(col2rgb(cl), 2, function(x){ rgb(x[1],x[2],x[3],a,maxColorValue=255)}) }


  V<-cbind(t(df[65,318:333]))
  dist_mat <-matrix(V,4,4)
  rownames(dist_mat) <-c("M1","M2","M3","M4")
  colnames(dist_mat) <-c("P1","P2","P3","P4")
  dist_mat <- dist_mat/sum(dist_mat)

  plot(0:4,0:4,xlab="",ylab="",col=NA)

  for(i in 1:4) points(1:4-0.5,rep(i-0.5,4),cex=dist_mat[i,]*50,pch=16,col="grey40")
  for(i in 1:4) points(1:4-0.5,rep(i-0.5,4),cex=dist_goal[i,]*50,pch=16,col=mTrsp(2,75))


  mtext("TB Progression",1,2.5,cex=0.9)
  mtext("Mortality",2,2.5,cex=0.9)

  box()
  mtext("Population by Joint Risk Factor Distribution, 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topright",c("model", "Goal Distribution"),pch=19,pt.cex=1,
         lwd=NA,col=c("grey40","darkred"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP HR DISTRIBUTION 1993-2013 ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V   <- cbind(df[44:65,305],df[44:65,306],df[44:65,307],df[44:65,308] )
  V[,1] <-V[,2]/(V[,1]+ V[,2]) #high risk percentage
  V <- V[,-2]
  V[,2] <-V[,3]/(V[,2]+ V[,3]) #high risk percentage
  V <- V[,-3]

  V4 <- cbind(df[44:65,305]+df[44:65,307], df[44:65,306]+df[44:65,08])
  V4 <- V4[,2]/(V4[,1]+V4[,2])

  plot(0,0,ylim=c(0,max(range(V4))*100),xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
  # lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)
  lines(1993:2014,V[,2]*100,lwd=2,col="red3")
  lines(1993:2014,V[,1]*100,lwd=2,col="blue")
  lines(1993:2014,V4*100,lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent High Risk Population for FB (red) and US (blue) 1993-2014 ",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","model"),cex=0.9,
         pch=c(15,15,15,NA),lwd=c(NA,NA,NA,2),lty=c(NA,NA,NA,1),col=c("grey50",4,"red3",1),bg="white",pt.cex=c(1.8,1.8,1.8,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

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
  legend("bottomright",c("Total","US born","Foreign born","Reported data","model"),cex=0.9,
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
  axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[1,i],V3[1,i]),border=NA,col="gray")
  for(i in 1:8) points(i+.2,(US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality by Age, 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("black","gray"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT RATE BY AGE DISTRIBUTION 1993-2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind((df[44:65,255:265])+(df[44:65,266:276]))

  V1  <- cbind(t(df[44:65,33:43])+t(df[44:65,44:54]))
  V2  <-  V/V1

  plot(0,0,ylim=c(0,1),xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  color=c("red", "orange", "gold", "green", "blue", "purple", "darkred", "light blue", "pink", "darkblue", "light green")
  for(i in 1:11) lines(1993:2014,V[,i],lwd=2,col=color[i])
  for(i in 1:11) lines(mubt[528])
  # for(i in 1:8) points(i+.2,(US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


  mtext("Year",1,2.5,cex=0.9)
  box()
  mtext("Mortality Rate by Age over Time, 1993-2014",3,.8,font=2,cex=0.8)
  legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
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
  axis(2,c(0.15,.3,.45,.6,.75,.90,1.05,1.15),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

  mtext("Age Group",1,2.5,cex=0.9)
  mtext("Number of Deaths",2,2.5,cex=0.9)

  box()
  mtext("Mortality by Age for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Non US Born", "US Born"),pch=c(15,15),pt.cex=c(2,2),
         lwd=c(NA,NA),col=c("pink","lightblue"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE % DISTRIBUTION 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,121:131]))
  V1 <-sum(V)
  V2 <-V/V1*100

  plot(0,0,ylim=c(0.05,max(range(V/V1))*100),xlim=c(0.6,11.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94", "95p"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:12-0.5,rep("",12))
  axis(2,c(0,5,10,15,20,25,30),las=2);box()
  abline(h=axTicks(2),col="grey85")

  for(i in 1:11) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V2[i,1],V2[i,1]),border=NA,col="grey")

  mtext("Age Group",1,2.5,cex=0.9)
  mtext("Percent of Total Deaths",2,2.5,cex=0.9)

  box()
  mtext("Percent of Total Mortality by Age Group 2014 (mil)",3,.8,font=2,cex=0.8)
  legend("topleft",c("Model"),pch=15,pt.cex=2,
         lwd=NA,col=c("gray"),bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION % 2014  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

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
#   legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
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
  legend("topright",c("model"),pch=15,pt.cex=2,lwd=NA,col="gray",bg="white")
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT MORT GROUP DISTRIBUTION 2014 ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#   V  <- cbind(t(df[65,285:288]),t(df[65,289:292]))
#   V1  <- colSums(V)
#
#   plot(0,0,ylim=c(0.05,50),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA)
#   axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
#   axis(1,1:5-0.5,rep("",5))
#   axis(2,c(0,10,20,30,40,50),las=2);box()
#   abline(h=axTicks(2),col="grey85")
#   for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,(V[i,1]/V1[1])*100,(V[i,1]/V1[1])*100),border=NA,col="lightblue")
#   for(i in 1:4) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,(V[i,2]/V1[2])*100,(V[i,2]/V1[2])*100),border=NA,col="pink")
#
#   mtext("Risk Group",1,2.5,cex=0.9)
#   box()
#   mtext("Percent of Mortality by Mortality Group for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
#   legend("topright",c("model"),pch=15,pt.cex=2,lwd=NA,col="gray",bg="white")
#
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#   ### ### ### ### ### ### PERCENT OF EACH MORT GROUP KILLED IN 2014 ### ### ### ### ### ### ###
#   ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#                   #US
  V  <- cbind(t(df[2:66,285:288])+t(df[2:66,289:292])) #mortality by mortality group
  V1  <- cbind(t(df[1:65,24:27])) #total population of each mortality group
  V2 <- V[,]/V1[,] #fraction of each mortality group dies
  plot(0,0,ylim=c(0,max(range(V2))*100),xlim=c(1950,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2014,V2[4,]*100,lwd=2,col="red3")
  lines(1950:2014,V2[3,]*100,lwd=2,col="orange")
  lines(1950:2014,V2[2,]*100,lwd=2,col="gold")
  # lines(1950:2014,V2[1,]*100,lwd=2,col="green")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent",2,2.5,cex=0.9)

  box()
  mtext("Percent of Each Mortality Group Killed, 1950-2014",3,.8,font=2,cex=0.8)
  legend("bottomright",c("MG1","MG2","MG3","MG4","Model"),pch=c(15,15,15,15,NA),lwd=c(NA,NA,NA,NA,2),lty=c(NA,NA,NA,NA,1), col=c("green","gold","orange","red3","black"),bg="white")


  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # ### ### ### ### ### ### Percent of each mort group killed by Age Group ### ### ### ### ### ### ###
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  #
  # V  <- cbind(t(df[65,285:288]),t(df[65,289:292]))
  # V1  <- colSums(V)
  #
  # plot(0,0,ylim=c(0.05,50),xlim=c(0.6,4.4),xlab="",ylab="",axes=F,col=NA)
  # axis(1,1:4,paste(c("1st","2nd","3rd","4th"),"\ngroup",sep=""),tick=F,cex.axis=0.75)
  # axis(1,1:5-0.5,rep("",5))
  # axis(2,c(0,10,20,30,40,50),las=2);box()
  # abline(h=axTicks(2),col="grey85")
  # for(i in 1:4) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,(V[i,1]/V1[1])*100,(V[i,1]/V1[1])*100),border=NA,col="lightblue")
  # for(i in 1:4) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,(V[i,2]/V1[2])*100,(V[i,2]/V1[2])*100),border=NA,col="pink")
  #
  # mtext("Risk Group",1,2.5,cex=0.9)
  # box()
  # mtext("Percent of Mortality by Mortality Group for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
  # legend("topright",c("model"),pch=15,pt.cex=2,lwd=NA,col="gray",bg="white")

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT % HR DISTRIBUTION 1993-2013 ### ### ### ### ### ###
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V   <- cbind(df[44:65,293],df[44:65,294],df[44:65,295],df[44:65,296] )
  V[,1] <-V[,2]/(V[,1]+ V[,2]) #high risk percentage
  V <- V[,-2]
  V[,2] <-V[,3]/(V[,2]+ V[,3]) #high risk percentage
  V <- V[,-3]

  V4 <- cbind(df[44:65,293]+df[44:65,295], df[44:65,294]+df[44:65,296])
  V4 <- V4[,2]/(V4[,1]+V4[,2])

  plot(0,0,ylim=c(0,max(range(V))*100),xlim=c(1993,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
  # lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)
  lines(1993:2014,V[,2]*100,lwd=2,col="red3")
  lines(1993:2014,V[,1]*100,lwd=2,col="blue")
  lines(1993:2014,V4*100,lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent of Total Mortality from High Risk Populations for FB (red) and US (blue) 1993-2014 ",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","model"),cex=0.9,
         pch=c(15,15,15,NA),lwd=c(NA,NA,NA,2),lty=c(NA,NA,NA,1),col=c("grey50",4,"red3",1),bg="white",pt.cex=c(1.8,1.8,1.8,NA))


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL MORT RATE 1950-2013 ### ### ### ### ### ###
  # ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V   <- cbind(df[1:65,510:520])
  x<-seq(6,781,12)
  V2  <-mubt[x,]

  col<-rainbow(11)

  plot(0,0,ylim=c(0,max(range(V), range(V2))),xlim=c(1950,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
  for (i in 1:11){
  lines(1950:2014,V[,i],lwd=3,col=col[i])
  lines(1950:2014,V2[,i],lty=3,lwd=2, col="black")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Age Specific Mortality Rates from 1950 to 2014",3,.8,font=2,cex=0.8)
  legend("topright",colnames(V),cex=0.9,
         pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))
  }
dev.off()
  }
