#'
#'This script creates graphs that compare model simulations
#'versus calibration data.
#'@name graphs_calib
#'@param df dataframe of results
#'@return .pdf of calibration graphs
#'@export
calib_graphs <- function(df, Par_list){
  library(MCMCpack)
df<-as.data.frame(df)
pdfname<-paste("MITUS_results/calib_graphs",Sys.time(),".pdf")
pdf(file=pdfname, width = 11, height = 8.5)
par(mfrow=c(2,2),mar=c(4,4.5,3,1))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

V  <- cbind(df[1:66,30], df[1:66,31]+df[1:66,32])
plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")
lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],lty=3,col="blue")
lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],lty=3,col="red3")
lines(1950:2015,V[,2],lwd=2,col="red3")
lines(1950:2015,V[,1],lwd=2,col="blue")
lines(1950:2015,rowSums(V),lwd=2,col="grey50")

mtext("Year",1,2.5,cex=1.2)
mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=1)
legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
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

plot(1,1,ylim=c(0.05,125),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y" )
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0,5,10,25,50,75,125),las=2);box()
abline(h=axTicks(2),col="grey85")

for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

points(1:8+0.2,CalibDat[["tot_pop16_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
points(1:8-0.2,CalibDat[["tot_pop16_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")

mtext("Age Group",1,2.5,cex=1.2)
box()
mtext("Total Population by Age Group 2014 (mil,log-scale)",3,.8,font=2,cex=1)
legend("topright",c("US born","Non-US Born","Reported data"),cex=1,
       pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014 all ages ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

V  <- cbind(t(df[65,33:43]),t(df[65,44:54]))

plot(1,1,ylim=c(0.05,max(range(V))),xlim=c(0.6,11.4),xlab="",ylab="",axes=F,col=NA, log="y")
axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94", "95p"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:12-0.5,rep("",12))
axis(2,c(0,.1,1,10,20,30,50),las=2);box()
abline(h=axTicks(2),col="grey85")

for(i in 1:11) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
for(i in 1:11) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")

mtext("Age Group",1,2.5,cex=1.2)
mtext("Millions",2,2.5,cex=1.2)

box()
mtext("Total Population by Age Group 2014 (mil,log-scale)",3,.8,font=2,cex=1.2)
legend("topright",c("US born","Non-US Born","Reported data"),cex=1.0,
       pch=c(15,15,19),lwd=c(NA,NA,1),lty=c(NA,NA,3),col=c("lightblue","pink",1),bg="white",pt.cex=c(1.8,1.8,0.3))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
V  <- cbind(rowSums(df[1:66,255:265]), rowSums(df[1:66,266:276]))
V1c <- rowSums(df[1:66,121:131])
plot(1,1,ylim=c(0,3.5),xlim=c(1950,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(1950:2015,V[,2],lwd=2,col="red3")
lines(1950:2015,V[,1],lwd=2,col="blue")
lines(1950:2015,V1c,lwd=2,col="grey50")
points(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,pch=19,cex=0.6,col="grey50")
lines(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,lty=3,col="grey50")

mtext("Year",1,2.5,cex=1.2)
mtext("Mortality: Total, US, and Non-US Born (mil)",3,.8,font=2,cex=1.2)
legend("topleft",c("Total","US born","Non-US Born","Reported data","model"),cex=1.0,
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


plot(0,0,ylim=c(0.05,max(range(V3))+.2),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
abline(h=axTicks(2),col="grey85")

for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[1,i],V3[1,i]),border=NA,col="gray")
for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")


mtext("Age Group",1,2.5,cex=1.2)
box()
mtext("Mortality by Age, 2014 (mil)",3,.8,font=2,cex=1.2)
legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("black","gray"),bg="white")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### TOTAL MORT RATE 1950-2013 ### ### ### ### ### ###
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
V   <- cbind(df[1:65,510:520])
x<-seq(6,781,12)
BgMort           <- Inputs[["BgMort"]]
mubt   <- matrix(NA,1801,11)
for(i in 1:11) {
  mubt[,i] <- SmoCurve(BgMort[,i+1])*P["TunMubt"]/12
}
V2  <-mubt[x,]

col<-rainbow(11)

plot(1,1,ylim=c(.00001,.04),xlim=c(1950,2014),xlab="",ylab="",axes=F, log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
# points(1993:2014,CalibDat$notif_us_hr[,1]/rowSums(CalibDat$notif_us_hr)*100,pch=19,cex=0.6)
for (i in 1:11){
  lines(1950:2014,V[,i],lwd=3,col=col[i])
  lines(1950:2014,V2[,i],lty=3,lwd=2, col="black")

  mtext("Year",1,2.5,cex=1.2)
  mtext("Age Specific Mortality Rates from 1950 to 2014",3,.8,font=2,cex=1.2)
  legend("bottomleft",colnames(V),cex=1.0,
    pch=rep(15,i),lwd=rep(NA,i),lty=rep(NA,i),col=col,bg="white",pt.cex=rep(1.8,i))
}

# graph of total diagnosed cases
# by total population, US born population, and non-US born population
V0 <- df[4:67,"NOTIF_ALL"]+df[4:67,"NOTIF_MORT_ALL"] #total population
V1 <- df[44:67,"NOTIF_US"]+df[44:67,"NOTIF_MORT_US"]   #US born population
V2 <- df[44:67,"NOTIF_F1"]+df[44:67,"NOTIF_F2"]+df[44:67,"NOTIF_MORT_F1"]+df[44:67,"NOTIF_MORT_F2"]   #non-US born population

#format the plot
plot(0,0,ylim=c(0,max(CalibDat[["tot_cases"]][,2])*1e3),xlim=c(1954,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
#multiply raw output by 1,000 to convert from millions to thousands
lines(1953:2016,V0*1e3,lwd=3,col="white"); lines(1953:2016,V0*1e3,lwd=2,col=1) #total population
lines(1993:2016,V1*1e3,lwd=3,col="white"); lines(1993:2016,V1*1e3,lwd=2,col=4) #US born population
lines(1993:2016,V2*1e3,lwd=3,col="white"); lines(1993:2016,V2*1e3,lwd=2,col=3) #non-US born population

#reported data for comparison
points(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,pch=19,cex=0.3) #total population
lines(CalibDat[["tot_cases"]][,1],(CalibDat[["tot_cases"]][,2])*1e3,lty=3,col=1)

notif_fb      <- cbind(CalibDat[["fb_cases"]][,2],1-CalibDat[["fb_cases"]][,2])*CalibDat[["fb_cases"]][,3]
notif_fb <-notif_fb/1000

points(1993:2016,notif_fb[,2],pch=19,cex=0.3,col=4) #US born population
lines(1993:2016,notif_fb[,2],pch=19,lty=3,col=4)

points(1993:2016,notif_fb[,1],pch=19,cex=0.3,col=3) #non-US born population
lines(1993:2016,notif_fb[,1],lty=3,col=3)

#plot text
mtext("Year",1,2.5,cex=1.2)
mtext("Total TB Cases Identified (000s), 1953-2014",3,.8,font=2,cex=1.2)
legend("bottomleft",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                    "Model (all)","Model (US born)","Model (non-US born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

# graph of total diagnosed cases
# by total population, US born population, and non-US born population
V0 <- df[57:67,"NOTIF_ALL"]+df[57:67,"NOTIF_MORT_ALL"] #total population
V1 <- df[57:67,"NOTIF_US"]+df[57:67,"NOTIF_MORT_US"]   #US born population
V2 <- df[57:67,"NOTIF_F1"]+df[57:67,"NOTIF_F2"]+df[57:67,"NOTIF_MORT_F1"]+df[57:67,"NOTIF_MORT_F2"]   #non-US born population

#format the plot
plot(0,0,ylim=c(0,max(CalibDat[["tot_cases"]][48:64,2])*1e3),xlim=c(2006,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
#multiply raw output by 1,000 to convert from millions to thousands
lines(2006:2016,V0*1e3,lwd=3,col="white"); lines(2006:2016,V0*1e3,lwd=2,col=1) #total population
lines(2006:2016,V1*1e3,lwd=3,col="white"); lines(2006:2016,V1*1e3,lwd=2,col=4) #US born population
lines(2006:2016,V2*1e3,lwd=3,col="white"); lines(2006:2016,V2*1e3,lwd=2,col=3) #non-US born population

#reported data for comparison
points(CalibDat[["tot_cases"]][54:64,1],(CalibDat[["tot_cases"]][54:64,2])*1e3,pch=19,cex=0.3) #total population
lines(CalibDat[["tot_cases"]][54:64,1],(CalibDat[["tot_cases"]][54:64,2])*1e3,lty=3,col=1)

notif_fb      <- cbind(CalibDat[["fb_cases"]][14:24,2],1-CalibDat[["fb_cases"]][14:24,2])*CalibDat[["fb_cases"]][14:24,3]
notif_fb <-notif_fb/1000

points(2006:2016,notif_fb[,2],pch=19,cex=0.3,col=4) #US born population
lines(2006:2016,notif_fb[,2],pch=19,lty=3,col=4)

points(2006:2016,notif_fb[,1],pch=19,cex=0.3,col=3) #non-US born population
lines(2006:2016,notif_fb[,1],lty=3,col=3)

#plot text
mtext("Year",1,2.5,cex=1.2)
mtext("Total TB Cases Identified (000s), 2000-2014",3,.8,font=2,cex=1.2)
legend("bottomleft",c("Reported data (all)","Reported data (US born)","Reported data (non-US born)",
                    "Model (all)","Model (US born)","Model (non-US born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

################################################################################
#Percent of Total Cases Non-US Born Population

V <- cbind(df[57:67,"NOTIF_US"]+df[57:67,"NOTIF_MORT_US"], #US born population
           df[57:67,"NOTIF_F1"]+df[57:67,"NOTIF_F2"]+  #non-US born population
           df[57:67,"NOTIF_MORT_F1"]+df[57:67,"NOTIF_MORT_F2"])
V <- V[,2]/rowSums(V)

#format the plot
plot(0,0,ylim=c(2.5,97.5),xlim=c(2006,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2006:2016,V*100,lwd=2,col=4)

#reported data for comparison
points(2006:2016,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
lines(2006:2016,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

#plot text
mtext("Year",1,2.5,cex=1.2)
mtext("Percent of TB Cases Non-US Born, 2006-2016",3,.8,font=2,cex=1.2)
legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
#Percent of Non-US Born Cases from Recent Immigrant Population

V <- cbind(df[55:65,"NOTIF_F1"]+df[55:65,"NOTIF_MORT_F1"],df[55:65,"NOTIF_F2"]+df[55:65,"NOTIF_MORT_F2"])
V <- V[,1]/rowSums(V)

#format the plot
plot(0,0,ylim=c(0,60),xlim=c(2004,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2004:2014,V*100,lwd=2,col=4)

#reported data for comparison
notif_fb_rec   <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]

points(2004:2014,notif_fb_rec[12:22,1]/rowSums(notif_fb_rec[12:22,])*100,pch=19,cex=0.6)
lines(2004:2014,notif_fb_rec[12:22,1]/rowSums(notif_fb_rec[12:22,])*100,lty=3)

#plot text
mtext("Year",1,2.5,cex=1.2)
mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=1.2)
legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
#Age distribution of Cases
#0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

V   <- (df[57:67,136:146]+df[57:67,189:199])
V2  <- V[,-11]
V2[,10] <- V2[,10]+V[,11]

#format the plot
cls <- colorRampPalette(c("blue", "red"))( 4 )
plot(0,0,ylim=c(0,6),xlim=c(2006,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2006:2016,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])    #0-24 yrs
lines(2006:2016,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])    #25-44 yrs
lines(2006:2016,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])    #45-64 yrs
lines(2006:2016,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])   #65+ yrs

#reported data for comparison
notif_age     <- CalibDat[["age_cases"]][14:24,-c(1,12)]*CalibDat[["age_cases"]][14:24,12]

points(2006:2016,rowSums(notif_age[,1:3])/1e3,pch=19,cex=0.6,col=cls[1]) #0-24 yrs
lines(2006:2016,rowSums(notif_age[,1:3])/1e3,col=cls[1],lty=3)
points(2006:2016,rowSums(notif_age[,4:5])/1e3,pch=19,cex=0.6,col=cls[2]) #25-44 yrs
lines(2006:2016,rowSums(notif_age[,4:5])/1e3,col=cls[2],lty=3)
points(2006:2016,rowSums(notif_age[,6:7])/1e3,pch=19,cex=0.6,col=cls[3]) #45-64 yrs
lines(2006:2016,rowSums(notif_age[,6:7])/1e3,col=cls[3],lty=3)
points(2006:2016,rowSums(notif_age[,8:10])/1e3,pch=19,cex=0.6,col=cls[4]) #65+ yrs
lines(2006:2016,rowSums(notif_age[,8:10])/1e3,col=cls[4],lty=3)

#plot text
mtext("TB Cases By Age (000s), 2006-16",3,.8,font=2,cex=1.2)
mtext("Year",1,2.5,cex=1.2)

legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
       lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
       pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

################################################################################
#Age Distribution of TB Cases in Percentages
#0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs
#
#   V   <- (df[51:65,136:146]+df[51:65,189:199])
#   V2  <- V[,-11]
#   V2[,10] <- V2[,10]+V[,11]
#
#   #format the plot
#   plot(0,0,ylim=c(0,max(range(V2))),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
#   axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
#        tick=F,cex.axis=0.6)
#   axis(1,1:11-0.5,rep("",11))
#   axis(2,las=2);box()
#   abline(h=axTicks(2),col="grey85")
#
#   #plot the model data
#   for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")
#
#   #reported data for comparison
#   points(1:10,colSums(notif_age[7:21,])/sum(notif_age[7:21,])*100,pch=19,cex=1.2)
#
#   #plot text
#   mtext("Age Group",1,2.5,cex=1.2)
#   mtext("Age Distribution of TB Cases (%), 2000-14",3,.8,font=2,cex=1.2)
#   legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
#          pt.cex=c(1,2),col=c("black","lightblue"),bg="white")
################################################################################
# Distribution of Cases across the TB Progression Categories


#format the plot
#plot the model data
#reported data for comparison
#plot text

################################################################################
# Treatment Outcomes 1993-2014

V   <- df[55:65,132:134]
Vdisc <- V[,2]/rowSums(V)
Vdead <- V[,3]/rowSums(V)

#format the plot
plot(0,0,ylim=c(0,10),xlim=c(2004,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2004:2014,Vdisc*100,lwd=2,col="red3")
lines(2004:2014,Vdead*100,lwd=2,col="blue")

#reported data for comparison
tx_outcomes      <- cbind((1-rowSums(CalibDat[["tx_outcomes"]][12:22,2:3])),CalibDat[["tx_outcomes"]][12:22,2],CalibDat[["tx_outcomes"]][12:22,3])*CalibDat[["tx_outcomes"]][12:22,4]

points(2004:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
points(2004:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
lines (2004:2014,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
lines (2004:2014,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")

#plot text

mtext("Year",1,2.5,cex=1.2)
mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=1.2)
legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

################################################################################
#LTBI Prevalance by Age in 2011, US born

V  <- cbind(t(df[62,55:65]),t(df[62,33:43]-df[62,55:65]))
Va <- outer(V[,1],c(P[["SensLt"]],1-P[["SensLt"]]))+outer(V[,2],c(1-P[["SpecLt"]],P[["SpecLt"]]))
colnames(V) <- c("LTBI", "No-LTBI")

V1 <- Va[-11,]; V1<-V1[-10,]
V1[9,] <- V1[9,]+Va[10,]+Va[11,]

V2 <- rep(NA,8)
V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

#format the plot
plot(0,0,ylim=c(0,8.5),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:8-0.5,rep("",8))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

#reported data for comparison
ltbi_us_11      <- CalibDat[["LTBI_prev_US_11_IGRA"]]
ltbi_fb_11      <- CalibDat[["LTBI_prev_FB_11_IGRA"]]

points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)

#plot text
mtext("Age Group",1,2.5,cex=1.2)
mtext("IGRA+ LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=1.2)
legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

################################################################################
#LTBI Prevalance by Age in 2011, non-US born

V  <- cbind(t(df[62,66:76]),t(df[62,44:54]-df[62,66:76]))
Va <- outer(V[,1],c(P[["SensLt"]],1-P[["SensLt"]]))+outer(V[,2],c(1-P[["SpecLt"]],P[["SpecLt"]]))

colnames(V) <- c("LTBI", "No-LTBI")

V1 <- Va[-11,]; V1<-V1[-10,]
V1[9,] <- V[9,]+Va[10,]+Va[11,]

V2 <- rep(NA,8)
V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100

#format the plot
plot(0,0,ylim=c(0,55),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:8-0.5,rep("",8))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

#reported data for comparison
points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=1.2)

#plot text
mtext("Age Group",1,2.5,cex=1.2)
mtext("IGRA+ LTBI in Non-US Born Population 2011 by Age (%)",3,.8,font=2,cex=1.2)
legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

################################################################################
# Age Distribution of TB Deaths 1999-2014

V  <- df[50:65,227:237]

V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
V3 <- colSums(V2)*1e6
tb_deaths      <- CalibDat[["tb_deaths"]][,-1]

#format the plot
plot(0,0,ylim=c(0,max(range(V3),range(colSums(tb_deaths)))+1000),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

#plot the model data
for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")

#reported data for comparison
points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

#plot text
mtext("Age Group",1,2.5,cex=1.2)
mtext("Total TB Deaths by Age Group 1999-2014",3,.8,font=2,cex=1.2)
legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

#######################################
# total tb deaths over time 2004-2014
V   <- rowSums(df[55:65,227:237])
tb_death_tot<-rowSums(CalibDat$tb_deaths[6:16,-1])

#format the plot
plot(0,0,ylim=c(0,max(tb_death_tot)*1.2),xlim=c(2006,2016),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#plot the model data
lines(2006:2016,V*1e6,lwd=2,col="blue")

#reported data for comparison
points(2006:2016,tb_death_tot,pch=19,cex=0.6,col="black")
lines (2006:2016,tb_death_tot,lty=3,col="black")

#plot text

mtext("Year",1,2.5,cex=1.2)
mtext("Total TB Deaths by Year 2004-2014",3,.8,font=2,cex=1.2)
legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),
       col=c("black","blue"),lty=c(3,1),bg="white",pt.cex=c(0.6,NA))################################################################################

################################################################################
graphs_pub(Par_list=Par_list)


dev.off()
# system(paste("open", pdfname))# code for ma

}

