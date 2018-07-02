#'This script creates graphs that compare model simulations
#'versus calibration data.
#'
#'load the necessary library
library(MCMCpack)

calib_graphs <- function(dataset){
#'load the parameter data set & model inputs
data("ParamInit_2018", package='MITUS')# ParamInit
data("ModelInputs_9-2-16", package='MITUS')
data("CalibDat_2018-06-28", package='MITUS')


#'Define the Posterior Function
#'@param theta
#'@return posterior
posterior <- function(theta){
  -lprior(theta) - llikelihood(theta,n_cores)
}
#' load the optimized data set
#load("~MITUS/data/Opt_US5401_r8_1_1-27-16.rData")
dataset #this is the name of the data set above

par_1 <- dataset$par

llikelihood(par_1)
lprior(par_1)
posterior(par_1)
M1    <- OutputsZint(ParMatrix=par_1)
M <- M1
#save(M,file="M_6-2-16c.rData")
lPrior2

Par2 <- pnorm(par_1,0,1)
# unif2true
Par3 <- Par2
Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
P[ii] <- Par3
P <- P
reldens <- Par3
reldens[idZ0] <- dbeta( Par3[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7],log=T)-dbeta( ParamInit$BcValue[ii][idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7],log=T)
reldens[idZ1] <- dgamma(Par3[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7],log=T)-dgamma(ParamInit$BcValue[ii][idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7],log=T)
reldens[idZ2] <- dnorm( Par3[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7],log=T)-dnorm( ParamInit$BcValue[ii][idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7],log=T)
reldens2 <- P
reldens2[] <- 0
reldens2[ii] <- round(reldens,2)

fn = paste("Parameter_table_", Sys.time(),".csv", sep='')
write.csv(cbind(ParamInit,P,reldens2),file=fn)

source("R/param.R")

#'Name and Open a .pdf file to which the plots will be written.
pdfnam <-paste("Fig_Calib_Checks_",Sys.time(),".pdf", sep = "")
pdf(file=pdfnam,width=11, height=8.5)
par(mfrow=c(2,2),mar=c(4,4.5,3,1))

#' graph of total diagnosed cases
#' by total population, US born population, and non-US born population
#' 1953-2013
V0 <- M[4:66,"NOTIF_ALL"]+M[4:66,"NOTIF_MORT_ALL"]   #total population
V1 <- M[44:66,"NOTIF_US"]+M[44:66,"NOTIF_MORT_US"]   #US born population
V2 <- M[44:66,"NOTIF_F1"]+M[44:66,"NOTIF_F2"]+M[44:66,"NOTIF_MORT_F1"]+M[44:66,"NOTIF_MORT_F2"]   #non-US born population

#'format the plot
plot(0,0,ylim=c(0,max(range(V0*1e3))),xlim=c(1954,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
lines(1953:2015,V0*1e3,lwd=3,col="white"); lines(1953:2015,V0*1e3,lwd=2,col=1) #total population
lines(1993:2015,V1*1e3,lwd=3,col="white"); lines(1993:2015,V1*1e3,lwd=2,col=4) #US born population
lines(1993:2015,V2*1e3,lwd=3,col="white"); lines(1993:2015,V2*1e3,lwd=2,col=3) #non-US born population

#'reported data for comparison
points(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,pch=19,cex=0.3) #total population
lines(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,lty=3,col=1)

points(1993:2015,notif_fb[,2]/1000,pch=19,cex=0.3,col=4) #US born population
lines(1993:2015,notif_fb[,2]/1000,pch=19,lty=3,col=4)

points(1993:2015,notif_fb[,1]/1000,pch=19,cex=0.3,col=3) #non-US born population
lines(1993:2015,notif_fb[,1]/1000,lty=3,col=3)

#'plot text
mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Cases Identified (000s), 1953-2014",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                    "Model (all)","Model (US born)","Model (foreign born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

#' 2000-2014
#' format the plot
plot(0,0,ylim=c(0,19),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#' plot the model data
lines(1953:2015,V0*1e3,lwd=3,col="white"); lines(1953:2015,V0*1e3,lwd=2,col=1) #total population
lines(1993:2015,V1*1e3,lwd=3,col="white"); lines(1993:2015,V1*1e3,lwd=2,col=4) #US born population
lines(1993:2015,V2*1e3,lwd=3,col="white"); lines(1993:2015,V2*1e3,lwd=2,col=3) #non-US born population

#' reported data for comparison
points(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,pch=19,cex=0.6)
points(1993:2015,notif_fb[,1]/1000,pch=19,cex=0.6,col=3)
points(1993:2015,notif_fb[,2]/1000,pch=19,cex=0.6,col=4)
lines(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,lty=3,col=1)
lines(1993:2015,notif_fb[,1]/1000,lty=3,col=3)
lines(1993:2015,notif_fb[,2]/1000,pch=19,lty=3,col=4)

#' plot text
mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Cases Identified (000s), 2000-2014",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                    "Fitted model (all)","Fitted model (US born)","Fitted model (foreign born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.6)

########################################################################################

#'Percent of Total Cases Non-US Born Population

V <- cbind(M[44:66,"NOTIF_US"]+M[44:66,"NOTIF_MORT_US"], #US born population
           M[44:66,"NOTIF_F1"]+M[44:66,"NOTIF_F2"]+  #non-US born population
             M[44:66,"NOTIF_MORT_F1"]+M[44:66,"NOTIF_MORT_F2"])
V <- V[,2]/rowSums(V)

#'format the plot
plot(0,0,ylim=c(2.5,97.5),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
lines(1993:2015,V*100,lwd=2,col=4)

#'reported data for comparison
points(1993:2015,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
lines(1993:2015,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)

#'plot text
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Non-US-Born, 2000-14",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
#'Percent of Non-US Born Cases from Recent Immigrant Population

V <- cbind(M[44:65,"NOTIF_F1"]+M[44:65,"NOTIF_MORT_F1"],M[44:65,"NOTIF_F2"]+M[44:65,"NOTIF_MORT_F2"])
V <- V[,1]/rowSums(V)

#'format the plot
plot(0,0,ylim=c(0,60),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
lines(1993:2014,V*100,lwd=2,col=4)

#'reported data for comparison
points(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,pch=19,cex=0.6)
lines(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,lty=3)

#'plot text
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of Non-US Born Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
#'Age distribution of Cases
#'0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs

V   <- (M[51:65,136:146]+M[51:65,189:199])
V2  <- V[,-11]
V2[,10] <- V2[,10]+V[,11]

#'format the plot
cls <- colorRampPalette(c("blue", "red"))( 4 )
plot(0,0,ylim=c(0,40),xlim=c(2000,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
lines(2000:2014,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])    #0-24 yrs
lines(2000:2014,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])    #25-44 yrs
lines(2000:2014,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])    #45-64 yrs
lines(2000:2014,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])   #65+ yrs

#'reported data for comparison
notif_age     <- CalibDat[["age_cases"]][,-c(1,12)]*CalibDat[["age_cases"]][,12]

points(1993:2014,rowSums(notif_age[,1:3])/1e3,pch=19,cex=0.6,col=cls[1]) #0-24 yrs
lines(1993:2014,rowSums(notif_age[,1:3])/1e3,col=cls[1],lty=3)
points(1993:2014,rowSums(notif_age[,4:5])/1e3,pch=19,cex=0.6,col=cls[2]) #25-44 yrs
lines(1993:2014,rowSums(notif_age[,4:5])/1e3,col=cls[2],lty=3)
points(1993:2014,rowSums(notif_age[,6:7])/1e3,pch=19,cex=0.6,col=cls[3]) #45-64 yrs
lines(1993:2014,rowSums(notif_age[,6:7])/1e3,col=cls[3],lty=3)
points(1993:2014,rowSums(notif_age[,8:10])/1e3,pch=19,cex=0.6,col=cls[4]) #65+ yrs
lines(1993:2014,rowSums(notif_age[,8:10])/1e3,col=cls[4],lty=3)

#'plot text
mtext("TB Cases By Age (000s), 2000-14",3,.8,font=2,cex=0.8)
mtext("Year",1,2.5,cex=0.9)

legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Model"),
       lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
       pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

################################################################################
#'Age Distribution of TB Cases in Percentages
#'0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs
#'
V   <- (M[51:65,136:146]+M[51:65,189:199])
V2  <- V[,-11]
V2[,10] <- V2[,10]+V[,11]
V2  <- colSums(V2)/sum(V2)*100

#'format the plot
plot(0,0,ylim=c(0,max(range(V2)+.5)),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
     tick=F,cex.axis=0.6)
axis(1,1:11-0.5,rep("",11))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

#'reported data for comparison
points(1:10,colSums(notif_age[7:21,])/sum(notif_age[7:21,])*100,pch=19,cex=1.2)

#'plot text
mtext("Age Group",1,2.5,cex=0.9)
mtext("Age Distribution of TB Cases (%), 2000-14",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

################################################################################
#'Are we interested in the distribution of cases across latent tx history?
################################################################################

#'Distribution of TB Cases Across High/Low Risk Cateogries 1993-2015
V   <- cbind(M[44:65,204],M[44:65,205]) + cbind(M[44:65,151],M[44:65,150])
V <- (V[,1]/rowSums(V))*100

#' format the plot
plot(0,0,ylim=c(0,(max(range(V))+2.5)),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#' plot the model data
lines(1993:2014,V,lwd=2,col=4)

#' reported data for comparison
points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)

#' plot text
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Homeless in Past Yr",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

################################################################################
#' Treatment Outcomes 1993-2011

V   <- M[44:63,132:134]
Vdisc <- V[,2]/rowSums(V)
Vdead <- V[,3]/rowSums(V)

#'format the plot
plot(0,0,ylim=c(0,15),xlim=c(1993,2011),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
lines(1993:2012,Vdisc*100,lwd=2,col="red3")
lines(1993:2012,Vdead*100,lwd=2,col="blue")

#'reported data for comparison
points(1993:2012,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
points(1993:2012,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
lines(1993:2012,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
lines(1993:2012,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")

#'plot text

mtext("Year",1,2.5,cex=0.9)
mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=0.8)
legend("topright",c("Discontinued","Died","Reported data","Model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

################################################################################
#'LTBI Initiations 1992-2011 Distribution
v13  <- M[43:65,153:154]/M[43:65,152]

#'format the plot
plot(1,1,ylim=c(0.001,.74)*100,xlim=c(1992,2015),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
for(i in 1:2) lines(1992:2014,v13[,i]*100,lwd=2,col=c("red3",4,6)[i])

#'reported data for comparison
points(rep(2002,3),TLTBI_dist*100,pch=19,cex=0.8,col=c("red3",4,6))
lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

#'plot text
mtext("Year",1,2.5,cex=0.9)
mtext("IPT Treatment Initiations By Risk Group (%)",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Foreign-born","Homeless","Reported data","Fitted model"),
       pch=c(15,15,19,NA),lwd=c(NA,NA,NA,2),col=c("red3",4,1,1),bg="white",pt.cex=c(1.8,1.8,0.8,NA))
################################################################################

#'LTBI Prevalance by Age in 2011, US born

V  <- cbind(M[62,55:65],(M[62,33:43]-M[62,55:65]))

V1 <- V[,-22]; V1<-V1[,-21]
V1[,20] <- V1[,20]+V[,21]+V[,22]

V1 <- V1[,-11]; V1<-V1[,-10]
V1[,9] <- V1[,9]+V[,10]+V[,11]

V2 <- rep(NA,8)
for (i in 2:9) V2[i] <- V1[,i] / (V1[,i]+V1[,i+9])
V2 <- V2[-1]*100

#'format the plot
plot(0,0,ylim=c(0,max(range(V2))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:8-0.5,rep("",8))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#'plot the model data
for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")

#'reported data for comparison
points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)

#'plot text
mtext("Age Group",1,2.5,cex=0.9)
mtext("LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=c(0,NA),
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

################################################################################
#'LTBI Prevalance by Age in 2011, non-US born

V  <- cbind(M[62,66:76],(M[62,44:54]-M[62,66:76]))

V1 <- V[,-22]; V1<-V1[,-21]
V1[,20] <- V1[,20]+V[,21]+V[,22]

V1 <- V1[,-11]; V1<-V1[,-10]
V1[,9] <- V1[,9]+V[,10]+V[,11]

V2 <- rep(NA,8)
for (i in 2:9) V2[i] <- V1[,i] / (V1[,i]+V1[,i+9])
V2 <- V2[-1]*100

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

V  <- M[50:65,88:98]+M[50:65,99:109]
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
points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

#'plot text
mtext("Age Group",1,2.5,cex=0.9)
mtext("Total TB Deaths by Age Group",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")
################################################################################
################################################################################
################################################################################
#'Begin Demographic Plots

#'Total Population Each Decade
#'US and non-US born

V  <- cbind(M[1:66,30], M[1:66,31]+M[1:66,32])

#' format the plot
plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

#' plot the model data
lines(1950:2015,V[,2],lwd=2,col="red3")
lines(1950:2015,V[,1],lwd=2,col="blue")
lines(1950:2015,rowSums(V),lwd=2,col="grey50")

#' reported data for comparison
points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],lty=3,col="grey50")
lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],lty=3,col="blue")
lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],lty=3,col="red3")

#'plot text
mtext("Year",1,2.5,cex=0.9)
mtext("Population: Total, US, and Foreign Born (mil, log-scale)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Total","US born","Foreign born","Reported data","model"),cex=0.9,
       pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
################################################################################

#'Total Population Age Distribution 2014

V  <- cbind(t(M[65,33:43]), t(M[65,44:54]))
V1  <- V[-3,]
V1[2,] <- V1[2,]+V[3,]
V2 <- V1[-4,]
V2[3,] <- V2[3,]+V1[4,]
V3 <- V2[-9,]
V3[8,] <- V3[8,]+V2[9,]

#' format the plot
plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0,10,25,40,50,60,75,78 ),las=2);box()
abline(h=axTicks(2),col="grey85")

#' plot the model data
for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[i,1],V3[i,1]),border=NA,col="lightblue")
for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V3[i,2],V3[i,2]),border=NA,col="pink")

#' reported data for comparison
points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")

#'plot text
mtext("Age Group",1,2.5,cex=0.9)
box()
mtext("Population by Age for FB (red) and US (blue), 2014 (mil)",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("grey30","grey80"),bg="white")

################################################################################
#' Recent TB infection 2014 (indices have been updated. )
Vall <- (M[65,156:171]/M[65,172:187])
#'format the plot
plot(-1,0,ylim=c(0.02,0.8),xlim=c(0.5,17.5),xlab="",ylab="",axes=F)
axis(2,las=2);box()
axis(1,1:16,c("All",paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),
            "yrs"),"US born","Foreign born","FB >2yrs","Homeless"),tick=F,cex.axis=0.7,las=2, mgp=c(3, 0.25, 0))
abline(h=axTicks(2),col="grey85")

#'plot the model data
for(i in 1:16) lines(rep(i,2),c(0,Vall[i]),col="forestgreen",lwd=10,lend="butt")

#'plot text
text(1:16,Vall,format(round(Vall,2),nsmall=2),cex=0.7,pos=3)
mtext("Fraction of Incident TB from Recent Infection (<2 years)",3,.8,font=2,cex=0.8)
################################################################################
################################################################################
### ### ### VARIOUS DEATHS DATA ### ### ### ### ### ###
# DdatDx <- read.csv("DeadAtDx_2-5-16.csv")
# load("CalibDat_4-1-16.rData")
# DdRpt <- read.csv("tb_deaths_fromreport_2-6-16.csv")
# Ddx <- rowSums(DdatDx[7:20,-1]) # dead at dx
# Dtx <- CalibDat$tx_outcomes[7:20,3]* CalibDat$tx_outcomes[7:20,4]
#
# #######
# plot(1999:2012,rowSums(CalibDat$tb_deaths[1:14,-1]),
#      las=1,ylim=c(0,2300),type="l",lwd=2,col=6,xlab="",ylab="Annual deaths")
# abline(h=axTicks(2),col="grey85")
# lines(1999:2012,Ddx+Dtx,las=1,ylim=c(0,1500),type="l",lwd=2,col=4)
# lines(1999:2012,Dtx,las=1,ylim=c(0,1500),type="l",lwd=2,col=3)
# lines(1999:2012,rowSums(CalibDat$tb_deaths[1:14,-1]),las=1,ylim=c(0,1500),type="l",lwd=2,col=2)
# mtext("Year",1,2.5,cex=0.9)
# mtext("Various Data Sources on TB Deaths",3,.8,font=2,cex=0.8)
# legend("topright",c("Dead at TB diagnosis or while on TB treatment","Died while on TB treatment","ICD-10 MCD TB codes"),
#        col=c(4,3,2),lwd=3,bg="white")

################################################################################
#' Other TB things! comparisons to published estimates
#'Reactivation Rates vs. Ferebee and Sutherland estimates
#'source in param_init.R or param.R (will need to finalize)
# source("R/param_init.R")
n_yr_Z = 50
pfast_v <- rslow_v <- rep(NA,4)
pfast_v <- pfast * exp(0:3)/3*log(ORpfastRF)
rslow_v <- rslow * exp(0:3)/3*log(RRrslowRF)
s2 <- rep(NA,4)
#' calculate the joint probability of pfast/pslow & not progressing to active disease & not recovering
for (i in 1:4){
  p2 <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*(0:n_yr_F)))*(rfast/(rfast+rRecov))) +
    (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*(0:n_yr_F)))*(rslow_v[i]/(rslow_v[i]+rRecov)))
  r2 <- -log(1-diff(-p2)/p2[-(n_yr_F+1)])/1
  s2[i]<- sum((datF[,3]*log(r2) + (datF[,2]-datF[,3])*log(1-r2))/N_red)-adj_25/N_red
}
V<-(sum(s2*colSums(dist_gen)))

#' #' format the plot
#' plot(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,type="l",lwd=2,xlim=c(0,15),ylim=c(0,60),las=1,xlab="",ylab="")
#'
#' #' plot the model data
#' lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2, col=)
#' lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2, col=)
#'
#' #' reported data for comparison
#' lines(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,lwd=2)
#' lines(1:nrow(datS)-.5,datS[,3]/datS[,2]*1000,lwd=2,lty=2)
#'
#' #' plot text
#' mtext("Years since infection",1,2.5,cex=0.9)
#' mtext("Reactivation rate per 1000 PY",2,3,cex=0.8)
#' legend("topright",c("Ferebee 1970","Sutherland 1968","Fitted model"),
#'        col=c(1,1,4),lty=c(1,2,1),lwd=2,cex=0.9)
#' mtext("Reactivation rate vs. Ferebee, Sutherland estimates",3,0.8,cex=0.8,font=2)
#' ################################################################################
#' #'Reactivation Rates vs. Ferebee and Sutherland estimates on Log Scale
#'
#' #'format the plot
#' plot(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,type="l",lwd=2,xlim=c(0,20),ylim=c(0.2,55),las=1,xlab="",ylab="",log="y")
#'
#' #' plot the model data
#' lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2, col=)
#' lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2, col=)
#'
#' #' reported data for comparison
#' lines(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,lwd=2)
#' lines(1:nrow(datS)-.5,datS[,3]/datS[,2]*1000,lwd=2,lty=2)
#'
#' #' plot text
#' mtext("Years since infection",1,2.5,cex=0.9)
#' mtext("Reactivation rate per 1000 PY (log scale)",2,3,cex=0.8)
#' legend("topright",c("Ferebee 1970","Sutherland 1968","Fitted model"),
#'        col=c(1,1,4),lty=c(1,2,1),lwd=2,cex=0.9)
#' mtext("Reactivation rate vs Ferebee & Sutherland estimates, log-scale",3,0.8,cex=0.8,font=2)
#' ################################################################################
#' #'Fraction without active TB (of those progressing <15yrs)
pfast_v <- rslow_v <- rep(NA,4)
pfast_v <- pfast * exp(0:3)/3*log(ORpfastRF)
rslow_v <- rslow * exp(0:3)/3*log(RRrslowRF)

for (i in 1:4){
  p <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*datB[,1]))*(rfast/(rfast+rRecov))) +
    (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*datB[,1]))*(rslow_v[i]/(rslow_v[i]+rRecov)))
  p <- 1-(1-p)/(1-p)[nrow(datB)]
  s[i]<-sum(diff(-datB[,2])*ss_borgdorff*log(diff(-p)) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-p)))/N_red - adj_24/N_red
}
s*(colSums(dist_gen)) #weighted average of the likelihoods
#'
#' #'format the plot
#' plot(datB[,1:2],type="l",lwd=1,lty=3,las=1,xlab="",ylab="")
#'
#' #' plot the model data
#' points(datB[,1:2],cex=0.8,pch=19, col=)
#' lines(datB[,1],p,col=4,lwd=2, col=)
#' points(datB[,1:2],cex=0.8,pch=19, col=)
#' lines(datB[,1],p,col=4,lwd=2, col=)
#'
#'
#' #' plot text
#' mtext("Years since infection",1,2.5,cex=0.9)
#' mtext("Fraction without active TB",2,3,cex=0.8)
#' legend("topright",c("Borgdorff 2011","Fitted model"), col=c(1,4),lwd=3,cex=0.9)
#' mtext("Fraction without active TB (of those progressing <15yrs)",3,0.8,cex=0.8,font=2)
#'
#' ################################################################################
#' #'Cumulative incidence of active TB following M.tb infection at age 15 (%)
#' #'needs to be updated!!!
#' vec <- rep(0,3);
#' names(vec) <- c("safe","slow","fast")
#' survM15 <- NULL
#' for(i in 1:(50*12)) {
#'   if(i==1) {
#'     vec <- c(0,1-P[["pfast"]][3,1],P[["pfast"]][3,1]) #need to update these values
#'   } else {
#'     vec[1] <- vec[1]+vec[2]*P[["rRecov"]]
#'     vec[4] <- vec[4]+vec[3]*P[["rfast"]]+vec[2]*rep(P[["rslow"]][3:11,1],each=10*12)[i]
#'     vec[2] <- vec[2]-vec[2]*P[["rRecov"]]-vec[2]*rep(P[["rslow"]][3:11,1],each=10*12)[i]
#'     vec[3] <- vec[3]-vec[3]*P[["rfast"]]  }
#'   survM15[i] <- 1- vec[4]  }
#'
#' #'format the plot
#' plot(-1,-1,col=0,las=1,ylim=c(0,10),xlim=c(0,50),xlab="",ylab="")
#' abline(h=axTicks(2),col="grey90");box()
#' #'plot the model data
#' lines(seq(0,50*12)/12,c(0,(1-survM15)*100),col=4,lwd=2)
#' #' plot text
#' mtext("Years since infection",1,2.5,cex=0.9)
#' mtext("Cumulative incidence of active TB following M.tb infection at age 15 (%)",3,0.8,cex=0.8,font=2)

################################################################################


dev.off(); system(paste("open", pdfnam)) # code for ma

}
