#'Create a function to be run on a specific model run output to
#'create simple graphs of all the output for a selected year range
#'@param df dataframe of output for all years
#'@param dist boolean, if you want generic risk group functions
#'@return .pdf of the graphs
#'@export

tb_graph_demomod <- function(df){
##total population
V  <- df[1:67,2]
plot(1,1,ylim=c(100,500),xlim=c(1950,2016),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")

lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")


lines(1950:2016,V,lwd=2,col="grey50")

mtext("Year",1,2.5,cex=0.9)
mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Total","Reported data","model"),cex=0.9,
       pch=c(15,19,NA),lwd=c(NA,1,2),lty=c(NA,3,1),col=c("grey50",1,1),bg="white",pt.cex=c(1.8,0.3,NA))

##population by age in 2016
V  <- df[67,3:13]
V1  <- V[-3]
V1[2] <- V1[2]+V[3]
V2 <- V1[-4]
V2[3] <- V2[3]+V1[4]
V3 <- V2[-9]
V3[8] <- V3[8]+V2[9]

plot(1,1,ylim=c(0.05,75),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y" )
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0,10,25,50,75 ),las=2);box()
abline(h=axTicks(2),col="grey85")

for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i],V3[i]),border=NA,col="grey")

points(1:8+0.2,(CalibDat[["tot_pop16_ag_fb"]][-9,3]+CalibDat[["tot_pop16_ag_fb"]][-9,4]),pch=19,cex=1.2,col="black")

mtext("Age Group",1,2.5,cex=0.9)
box()
mtext("Population by Age, 2016 (mil,log-scale)",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("grey30","grey80"),bg="white")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
V  <- rowSums(df[1:67,14:24])
plot(1,1,ylim=c(0,max(range(V))),xlim=c(1950,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(1950:2016,V,lwd=2,col="grey50")
points(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,pch=19,cex=0.6,col="grey50")
lines(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,lty=3,col="grey50")

mtext("Year",1,2.5,cex=0.9)
mtext("Total Mortality 1950-2016",3,.8,font=2,cex=0.8)
legend("bottomright",c("Total","Reported data","model"),cex=0.9,
       pch=c(15,19,NA),lwd=c(NA,1,2),lty=c(NA,3,1),col=c("grey50",1,1),bg="white",pt.cex=c(1.8,0.3,NA))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2016  ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

V  <- df[67,14:24]
V1  <- V[-3]
V1[2] <- V1[2]+V[3]
V2 <- V1[-4]
V2[3] <- V2[3]+V1[4]
V3 <- V2[-9]
V3[8] <- V3[8]+V2[9]
V3<-V3/sum(V3)

plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
abline(h=axTicks(2),col="grey85")

for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i],V3[i]),border=NA,col="gray")
for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[67,i+1])/rowSums(CalibDat$US_mort_age[67,2:9]),pch=19,cex=1.2,col="black")


mtext("Age Group",1,2.5,cex=0.9)
box()
mtext("Mortality by Age, 2016 (%)",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("black","gray"),bg="white")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2016  ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

V  <- df[1,14:24]
V1  <- V[-3]
V1[2] <- V1[2]+V[3]
V2 <- V1[-4]
V2[3] <- V2[3]+V1[4]
V3 <- V2[-9]
V3[8] <- V3[8]+V2[9]
V3<-V3/sum(V3)

plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
abline(h=axTicks(2),col="grey85")

for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i],V3[i]),border=NA,col="gray")
for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[1,i+1])/rowSums(CalibDat$US_mort_age[1,2:9]),pch=19,cex=1.2,col="black")


mtext("Age Group",1,2.5,cex=0.9)
box()
mtext("Mortality by Age, 1950 (%)",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("black","gray"),bg="white")
}
