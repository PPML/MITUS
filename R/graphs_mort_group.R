
mugrp_graphs <-function(df){


  pdf(file=paste("MITUS_results/graphs_mugrp",Sys.time(),".pdf"), width = 11, height = 8.5)
  par(mfrow=c(2,2),mar=c(4,4.5,3,1))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ANNUAL MORT RATE FOR MORT GROUPS ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

V  <- cbind(t(df[2:66,285:288])+t(df[2:66,289:292])) #mortality by mortality group
V1  <- cbind(t(df[1:65,24:27])) #total population of each mortality group
V2 <- V[,]/V1[,] #fraction of each mortality group dies
V3<-V1/colSums(V1)
plot(0,0,ylim=c(0,max(range(V3))),xlim=c(1950,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(1950:2014,V3[4,],lwd=2,col="red3")
lines(1950:2014,V3[3,],lwd=2,col="orange")
lines(1950:2014,V3[2,],lwd=2,col="gold")
lines(1950:2014,V3[1,],lwd=2,col="green")

mtext("Year",1,2.5,cex=0.9)
mtext("Percent",2,2.5,cex=0.9)

box()
mtext("Annual Mortality Rate by Mortality Group, 1950-2014",3,.8,font=2,cex=0.8)
legend("bottomleft",c("MG1","MG2","MG3","MG4","Model"),pch=c(15,15,15,15,NA),lwd=c(NA,NA,NA,NA,2),lty=c(NA,NA,NA,NA,1), col=c("green","gold","orange","red3","black"),bg="white")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##for age group 1
for (i in 1:11){
V  <- df[1:65,(521+(4*(i-1))):(524+(4*(i-1)))] #mortality by mortality group
plot(0,0,ylim=c(0,max(range(V))),xlim=c(1950,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(1950:2014,V[,4],lwd=2,col="red3")
lines(1950:2014,V[,3],lwd=2,col="orange")
lines(1950:2014,V[,2],lwd=2,col="gold")
lines(1950:2014,V[,1],lwd=2,col="green")

mtext("Year",1,2.5,cex=0.9)
mtext("Percent",2,2.5,cex=0.9)

box()
mtext("Annual Mortality Rate by Mortality Group & Age, 1950-2014",3,.8,font=2,cex=0.8)
legend("bottomleft",c("MG1","MG2","MG3","MG4","Model"),pch=c(15,15,15,15,NA),lwd=c(NA,NA,NA,NA,2),lty=c(NA,NA,NA,NA,1), col=c("green","gold","orange","red3","black"),bg="white")
}

for (i in 1:11){
  V  <- cbind(df[1:65,(334+(16*(i-1)))]+df[1:65,(335+(16*(i-1)))]+df[1:65,(336+(16*(i-1)))]+df[1:65,(337+(16*(i-1)))],
              df[1:65,(338+(16*(i-1)))]+df[1:65,(339+(16*(i-1)))]+df[1:65,(340+(16*(i-1)))]+df[1:65,(341+(16*(i-1)))],
              df[1:65,(342+(16*(i-1)))]+df[1:65,(343+(16*(i-1)))]+df[1:65,(344+(16*(i-1)))]+df[1:65,(345+(16*(i-1)))],
              df[1:65,(346+(16*(i-1)))]+df[1:65,(347+(16*(i-1)))]+df[1:65,(348+(16*(i-1)))]+df[1:65,(348+(16*(i-1)))]) #mortality by mortality group
  plot(0,0,ylim=c(0,max(range(V))),xlim=c(1950,2014),xlab="",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2014,V[,4],lwd=2,col="red3")
  lines(1950:2014,V[,3],lwd=2,col="orange")
  lines(1950:2014,V[,2],lwd=2,col="gold")
  lines(1950:2014,V[,1],lwd=2,col="green")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Percent",2,2.5,cex=0.9)

  box()
  mtext("Population by Mortality Group & Age, 1950-2014",3,.8,font=2,cex=0.8)
  legend("bottomleft",c("MG1","MG2","MG3","MG4","Model"),pch=c(15,15,15,15,NA),lwd=c(NA,NA,NA,NA,2),lty=c(NA,NA,NA,NA,1), col=c("green","gold","orange","red3","black"),bg="white")
}


dev.off()





}
