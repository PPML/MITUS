setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC Large Grant/Analysis Transmission")

###### GRAPHS COMPARING SIMS VS CALIB DATA ############
## Parameter stuff
load("ParamInitUS_V732tab.rData") # ParamInit
P  <- ParamInit[,1]; names(P) <- rownames(ParamInit)
ii <-  ParamInit[,5]==1
ParamInitZ <- ParamInit[ParamInit$Calib==1,]
idZ0 <- ParamInitZ[,4]==0
idZ1 <- ParamInitZ[,4]==1
idZ2 <- ParamInitZ[,4]==2

## Scripts and functions
load("ModelInputs_4-15-16.rData")
source("ParamUS_V501.r")
source("PriorFuncUS_V22.r")
source("CalibFunctionsUS_V20.r")
source("TimeStepUS_V891.r")
source("IMISfunctionsUS_V231.r")

posterior = function(theta) { -lprior(theta) - llikelihood(theta,n_cores) }


source("ProcessSims_V231.r")

load("Opt_US522_r10_15_1-27-16.rData")
o10
par_1 <- o10$par

  llikelihood(par_1)
  lprior(par_1)
  posterior(par_1)
  M1    <- Outputs(par_1)
  M <- M1
#save(M,file="M_6-2-16c.rData")
  lPrior2

# norm2unif
Par2 <- pnorm(par_1,0,1)
# unif2true
Par3 <- Par2
Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])  
P[ii] <- Par3
P <<- P
reldens <- Par3
reldens[idZ0] <- dbeta( Par3[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7],log=T)-dbeta( ParamInit$BcValue[ii][idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7],log=T)
reldens[idZ1] <- dgamma(Par3[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7],log=T)-dgamma(ParamInit$BcValue[ii][idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7],log=T)
reldens[idZ2] <- dnorm( Par3[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7],log=T)-dnorm( ParamInit$BcValue[ii][idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7],log=T)
reldens2 <- P
reldens2[] <- 0
reldens2[ii] <- round(reldens,2)
write.csv(cbind(ParamInit,P,reldens2),file="Parameter_table_8-28-16.csv")

source("ParamUS_V501.r")

############
# colnames(M)
pdfnam <- "Fig_Calib_Checks_8-28-16.pdf"

### ### ### TOTAL DIAGNOSED CASES 1953-2013  ### ### ### ### ### ### 
pdf(file=pdfnam,width=11, height=8.5) 
par(mfrow=c(2,2),mar=c(4,4.5,3,1))
V0   <- M[4:66,"NOTIF_ALL"]+M[4:66,"NOTIF_MORT_ALL"]
V1   <- rowSums(M[44:66,235:245]+M[44:66,246:256])
V2   <- rowSums((M[44:66,137:147]+M[44:66,216:226]) - (M[44:66,235:245]+M[44:66,246:256]))

plot(0,0,ylim=c(0,85),xlim=c(1954,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,pch=19,cex=0.3)
points(1993:2015,notif_fb[,1]/1000,pch=19,cex=0.3,col=3)
points(1993:2015,notif_fb[,2]/1000,pch=19,cex=0.3,col=4)
lines(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,lty=3,col=1)
lines(1993:2015,notif_fb[,1]/1000,lty=3,col=3)
lines(1993:2015,notif_fb[,2]/1000,pch=19,lty=3,col=4)

lines(1953:2015,V0*1e3,lwd=3,col="white"); lines(1953:2015,V0*1e3,lwd=2,col=1)
lines(1993:2015,V1*1e3,lwd=3,col="white"); lines(1993:2015,V1*1e3,lwd=2,col=4)
lines(1993:2015,V2*1e3,lwd=3,col="white"); lines(1993:2015,V2*1e3,lwd=2,col=3)

mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Cases Identified (000s), 1953-2014",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                    "Fitted model (all)","Fitted model (US born)","Fitted model (foreign born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

#############

plot(0,0,ylim=c(0,19),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,pch=19,cex=0.6)
points(1993:2015,notif_fb[,1]/1000,pch=19,cex=0.6,col=3)
points(1993:2015,notif_fb[,2]/1000,pch=19,cex=0.6,col=4)
lines(CalibDat[["tot_cases"]][,1],CalibDat[["tot_cases"]][,2]/1000,lty=3,col=1)
lines(1993:2015,notif_fb[,1]/1000,lty=3,col=3)
lines(1993:2015,notif_fb[,2]/1000,pch=19,lty=3,col=4)

lines(1953:2015,V0*1e3,lwd=3,col="white"); lines(1953:2015,V0*1e3,lwd=2,col=1)
lines(1993:2015,V1*1e3,lwd=3,col="white"); lines(1993:2015,V1*1e3,lwd=2,col=4)
lines(1993:2015,V2*1e3,lwd=3,col="white"); lines(1993:2015,V2*1e3,lwd=2,col=3)

mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Cases Identified (000s), 2000-2014",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data (all)","Reported data (US born)","Reported data (foreign born)",
                    "Fitted model (all)","Fitted model (US born)","Fitted model (foreign born)"),
       pch=c(19,19,19,NA,NA,NA),lwd=c(1,1,1,2,2,2),lty=c(3,3,3,1,1,1),col=c(1,4,3,1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.6)

### ### ### CASES FB DISTRIBUTION 1993-2013  ### ### ### ### ### ###
V   <- cbind(M[44:66,154]+M[44:66,155]+(M[44:66,233]+M[44:66,234]),
              M[44:66,152]+M[44:66,153]+(M[44:66,231]+M[44:66,232]))
V <- V[,1]/rowSums(V )
plot(0,0,ylim=c(2.5,97.5),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2015,notif_fb[,1]/rowSums(notif_fb)*100,pch=19,cex=0.6)
lines(1993:2015,notif_fb[,1]/rowSums(notif_fb)*100,lty=3)
lines(1993:2015,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Foreign-Born, 2000-14",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ### 
V   <- cbind(M[44:65,154]+M[44:65,233],M[44:65,155]+M[44:65,234])
V <- V[,1]/rowSums(V)

plot(0,0,ylim=c(0,60),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,pch=19,cex=0.6)
lines(1993:2014,notif_fb_rec[,1]/rowSums(notif_fb_rec)*100,lty=3)
lines(1993:2014,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of FB Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### CASES AGE DISTRIBUTION 2000-2014  ### ### ### ### ### ### 
notif_age     <- CalibDat[["age_cases"]][,-c(1,12)]*CalibDat[["age_cases"]][,12]
V   <- (M[51:65,137:147]+M[51:65,216:226])
V2  <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
cls <- colorRampPalette(c("blue", "red"))( 4 ) 
plot(0,0,ylim=c(0,6.5),xlim=c(2000,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,rowSums(notif_age[,1:3])/1e3,pch=19,cex=0.6,col=cls[1])
lines(1993:2014,rowSums(notif_age[,1:3])/1e3,col=cls[1],lty=3)
lines(2000:2014,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])
points(1993:2014,rowSums(notif_age[,4:5])/1e3,pch=19,cex=0.6,col=cls[2])
lines(2000:2014,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])
lines(1993:2014,rowSums(notif_age[,4:5])/1e3,col=cls[2],lty=3)
points(1993:2014,rowSums(notif_age[,6:7])/1e3,pch=19,cex=0.6,col=cls[3])
lines(1993:2014,rowSums(notif_age[,6:7])/1e3,col=cls[3],lty=3)
lines(2000:2014,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])
points(1993:2014,rowSums(notif_age[,8:10])/1e3,pch=19,cex=0.6,col=cls[4])
lines(2000:2014,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])
lines(1993:2014,rowSums(notif_age[,8:10])/1e3,col=cls[4],lty=3)

mtext("TB Cases By Age (000s), 2000-14",3,.8,font=2,cex=0.8)
mtext("Year",1,2.5,cex=0.9)

legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years","Reported data","Fitted model"),
       lwd=c(NA,NA,NA,NA,1,2),lty=c(NA,NA,NA,NA,3,1),col=c(cls,1,1),bg="white",
       pt.cex=c(1.8,1.8,1.8,1.8,0.6,NA),pch=c(15,15,15,15,19,NA))

### ### ### CASES AGE DISTRIBUTION 2000-2014  ### ### ### ### ### ###
V   <- (M[51:65,137:147]+M[51:65,216:226])
V2  <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  V2  <- colSums(V2)/sum(V2)*100

  plot(0,0,ylim=c(0,20),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
  axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
       tick=F,cex.axis=0.6)
  axis(1,1:11-0.5,rep("",11))
  axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")
  points(1:10,colSums(notif_age[7:21,])/sum(notif_age[7:21,])*100,pch=19,cex=1.2)
  mtext("Age Group",1,2.5,cex=0.9)
  mtext("Age Distribution of TB Cases (%), 2000-14",3,.8,font=2,cex=0.8)
 legend("topright",c("Reported data","Fitted model"),pch=c(19,15),lwd=NA,
        pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### CASES TX HISTORY DISTRIBUTION 1993-2013  ### ### ### ### ### ### 
V  <- M[44:65,148:149]+M[44:65,227:228]
V <- V[,2]/rowSums(V )

plot(0,0,ylim=c(0,13),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_prev[,2]/rowSums(notif_prev)*100,pch=19,cex=0.6)
lines(1993:2014,notif_prev[,2]/rowSums(notif_prev)*100,lty=3)
lines(1993:2014,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Previously Treated",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### CASES HIV DISTRIBUTION 1993-2013  ### ### ### ### ### ### 
V   <- M[44:65,150:151]+M[44:65,229:230]
V <- V[,1]/rowSums(V )
plot(0,0,ylim=c(0,50),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

points(1993:2014,notif_hiv0[,2]/rowSums(notif_hiv0[,2:4])*100,pch=19,cex=0.6)
lines(1993:2014,notif_hiv0[,2]/rowSums(notif_hiv0[,2:4])*100,lty=3)

points(1993:2014,notif_hiv0[,2]/rowSums(notif_hiv0[,2:3])*100,pch=19,cex=0.6,col="grey50")
lines(1993:2014,notif_hiv0[,2]/rowSums(notif_hiv0[,2:3])*100,lty=3,col="grey50")


lines(1993:2014,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases HIV-Positive",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data (low estimate)","Reported data (high estimate)","Fitted model"),pch=c(19,19,NA),lwd=c(1,1,2),col=c(1,"grey50",4),lty=c(3,3,1),bg="white",pt.cex=0.6)

### ### ### CASES HR DISTRIBUTION 1993-2013  ### ### ### ### ### ### 
V   <- cbind(M[44:65,153],M[44:65,152]) + cbind(M[44:65,232],M[44:65,231])
V <- V[,1]/rowSums(V )

plot(0,0,ylim=c(0,15),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)
lines(1993:2014,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Homeless in Past Yr",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### CASES PCT INH & MDR RES US N 1993-2013  ### ### ### ### ### ### 
V   <- M[44:65,156:160]
Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)

plot(0,0,ylim=c(0,7.5),xlim=c(1998,2013),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_inh_us_n[,1]/rowSums(notif_inh_us_n)*100,pch=19,cex=0.6,col="red3")
points(1993:2014,notif_mdr_us_n[,1]/rowSums(notif_mdr_us_n)*100,pch=19,cex=0.6,col="blue")
lines(1993:2014,notif_inh_us_n[,1]/rowSums(notif_inh_us_n)*100,lty=3,col="red3")
lines(1993:2014,notif_mdr_us_n[,1]/rowSums(notif_mdr_us_n)*100,lty=3,col="blue")
lines(1993:2014,Vinh[,1]*100,lwd=2,col="red3")
lines(1993:2014,Vmdr[,1]*100,lwd=2,col="blue")
mtext("Year",1,2.5,cex=0.9)
mtext("INH Resistance and MDR-TB in Tx Naive US Cases (%)",3,1,font=2,cex=0.8)
legend("topleft",c("INH resistance","MDR-TB","Reported data","Fitted model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

### ### ### CASES PCT INH & MDR RES US E 1993-2013  ### ### ### ### ### ### 
V   <- M[44:65,161:165]
Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)
#plot(Vmdr[,1])
plot(0,0,ylim=c(0,13),xlim=c(1998,2013),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_inh_us_e[,1]/rowSums(notif_inh_us_e)*100,pch=19,cex=0.6,col="red3")
points(1993:2014,notif_mdr_us_e[,1]/rowSums(notif_mdr_us_e)*100,pch=19,cex=0.6,col="blue")
lines(1993:2014,notif_inh_us_e[,1]/rowSums(notif_inh_us_e)*100,lty=3,col="red3")
lines(1993:2014,notif_mdr_us_e[,1]/rowSums(notif_mdr_us_e)*100,lty=3,col="blue")
lines(1993:2014,Vinh[,1]*100,lwd=2,col="red3")
lines(1993:2014,Vmdr[,1]*100,lwd=2,col="blue")
mtext("Year",1,2.5,cex=0.9)
mtext("INH Resistance and MDR-TB in Tx Experienced US Cases (%)",3,1,font=2,cex=0.8)
legend("topright",c("INH resistance","MDR-TB","Reported data","Fitted model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

### ### ### CASES PCT INH & MDR RES FB N 1993-2013  ### ### ### ### ### ### 
V   <- M[44:65,166:170]
Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)

plot(0,0,ylim=c(0,16),xlim=c(1998,2013),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_inh_fb_n[,1]/rowSums(notif_inh_fb_n)*100,pch=19,cex=0.6,col="red3")
points(1993:2014,notif_mdr_fb_n[,1]/rowSums(notif_mdr_fb_n)*100,pch=19,cex=0.6,col="blue")
lines(1993:2014,notif_inh_fb_n[,1]/rowSums(notif_inh_fb_n)*100,lty=3,col="red3")
lines(1993:2014,notif_mdr_fb_n[,1]/rowSums(notif_mdr_fb_n)*100,lty=3,col="blue")
lines(1993:2014,Vinh[,1]*100,lwd=2,col="red3")
lines(1993:2014,Vmdr[,1]*100,lwd=2,col="blue")
mtext("Year",1,2.5,cex=0.9)
mtext("INH Resistance and MDR-TB in Tx Naive FB Cases (%)",3,1,font=2,cex=0.8)
legend("topright",c("INH resistance","MDR-TB","Reported data","Fitted model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

### ### ### CASES PCT INH & MDR RES FB E 1993-2013  ### ### ### ### ### ### 
V   <- M[44:65,171:175]
Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)

plot(0,0,ylim=c(0,33),xlim=c(1998,2013),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2014,notif_inh_fb_e[,1]/rowSums(notif_inh_fb_e)*100,pch=19,cex=0.6,col="red3")
points(1993:2014,notif_mdr_fb_e[,1]/rowSums(notif_mdr_fb_e)*100,pch=19,cex=0.6,col="blue")
lines(1993:2014,notif_inh_fb_e[,1]/rowSums(notif_inh_fb_e)*100,lty=3,col="red3")
lines(1993:2014,notif_mdr_fb_e[,1]/rowSums(notif_mdr_fb_e)*100,lty=3,col="blue")
lines(1993:2014,Vinh[,1]*100,lwd=2,col="red3")
lines(1993:2014,Vmdr[,1]*100,lwd=2,col="blue")
mtext("Year",1,2.5,cex=0.9)
mtext("INH Resistance and MDR-TB in Tx Experienced FB Cases (%)",3,1,font=2,cex=0.8)
legend("topleft",c("INH resistance","MDR-TB","Reported data","Fitted model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

### ### ### TREATMENT OUTCOMES 1993-2011  ### ### ### ### ### ### 
V   <- M[44:63,133:135]
Vdisc <- V[,2]/rowSums(V)
Vdead <- V[,3]/rowSums(V)

plot(0,0,ylim=c(0,15),xlim=c(1993,2011),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1993:2012,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
points(1993:2012,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
lines(1993:2012,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
lines(1993:2012,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")
lines(1993:2012,Vdisc*100,lwd=2,col="red3")
lines(1993:2012,Vdead*100,lwd=2,col="blue")
mtext("Year",1,2.5,cex=0.9)
mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=0.8)
legend("topright",c("Discontinued","Died","Reported data","Fitted model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

### ### ### TLTBI VOL 1993-2011  ### ### ### ### ### ### 
v12  <- M[43:65,176]

plot(0,0,ylim=c(0,500),xlim=c(1992,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(2002,tltbi_vol[1]/1e3,pch=19,cex=0.8,col="black")
lines(1992:2014,v12*1e3,lwd=2,col="blue")
lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

mtext("Year",1,2.5,cex=0.9)
mtext("IPT Treatment Initiations Per Year (000s)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(2,2),col=c(1,"blue"),bg="white",pt.cex=0.8)

### ### ### LTBI INITIATIONS 1993-2011 Distribution ### ### ### ### ### ### 
v13  <- M[43:65,177:179]/M[43:65,176]

plot(1,1,ylim=c(0.001,.74)*100,xlim=c(1992,2015),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(rep(2002,3),TLTBI_dist*100,pch=19,cex=0.8,col=c("red3",4,6))
for(i in 1:3) lines(1992:2014,v13[,i]*100,lwd=2,col=c("red3",4,6)[i])
lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

mtext("Year",1,2.5,cex=0.9)
mtext("IPT Treatment Initiations By Risk Group (%)",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Foreign-born","Homeless","HIV","Reported data","Fitted model"),
       pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,NA,2),col=c("red3",4,6,1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.8,NA))

### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ### 
V  <- cbind(M[62,56:66],M[62,34:44]-M[62,56:66])
V[9,] <- colSums(V[9:11,])
V <- V[2:9,1]/rowSums(V[2:9,])*100
plot(0,0,ylim=c(0,9.5),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:8-0.5,rep("",8))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V[i],V[i]),border="white",col="lightblue")
points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)
mtext("Age Group",1,2.5,cex=0.9)
mtext("LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=c(0,NA),
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ### 
V  <- cbind(M[62,67:77],M[62,45:55]-M[62,67:77])
V[9,] <- colSums(V[9:11,])
V <- V[2:9,1]/rowSums(V[2:9,])*100

plot(0,0,ylim=c(0,52),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:8-0.5,rep("",8))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V[i],V[i]),border="white",col="lightblue")
points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=1.2)
mtext("Age Group",1,2.5,cex=0.9)
mtext("LTBI in Foreign Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=c(0,NA),
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### TOTAL POP EACH DECADE, BY US/FB  ### ### ### ### ### ### 
V  <- M[,c(30,32)]+M[,c(31,33)]

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

### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ### 
V  <- cbind(M[65,34:44],M[65,45:55])
V <- rbind(V[1,],V[2,]+V[3,],V[4,]+V[5,],V[6,],V[7,],V[8,],V[9,],V[10,]+V[11,])

plot(0,1,ylim=c(0.05,135),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y")
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0.1,1,10,100),las=2);box()
abline(h=axTicks(2),col="grey85")
for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")
points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")
mtext("Age Group",1,2.5,cex=0.9)
box()
mtext("Population by Age for FB (red) and US (blue), 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("grey30","grey80"),bg="white")

### ### ### TB DEATHS 1999-2013 ### ### ### ### ### ### 
V  <- M[50:65,89:99]+M[50:65,100:110]
V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]

plot(0,0,ylim=c(0,2500),xlim=c(1998.5,2014.5),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(1999:2014,rowSums(tb_deaths),pch=19,cex=0.6,col=1)
lines(1999:2014,rowSums(tb_deaths),lty=3,col=1)
lines(1999:2014,rowSums(V2)*1e6,lwd=2,col="blue")

mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Deaths by Year",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### AGE DISTRIBUTION, TB DEATHS 1999-2013 ### ### ### ### ### ### 
V  <- M[50:65,89:99]+M[50:65,100:110]
V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
V3 <- colSums(V2)*1e6

plot(0,0,ylim=c(0,5500),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")
points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

mtext("Age Group",1,2.5,cex=0.9)
mtext("Total TB Deaths by Age Group",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### HIV DEATHS 1999-2013 ### ### ### ### ### ###
V  <-M[50:66,111:121]
V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]

plot(0,0,ylim=c(0,35000),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(2010:2013,hv_deaths_tot,pch=19,cex=0.6,col=1)
lines(2010:2013,hv_deaths_tot,lty=3,col=1)
for(i in 1:4) lines(rep(2009+i,2),hv_deaths_tot[i]*c(0.8,1.2),lty=1,col=1)

lines(1999:2015,rowSums(V)*1e6,lwd=2,col="blue")

mtext("Year",1,2.5,cex=0.9)
mtext("Total HIV Deaths by Year",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### AGE DISTRIBUTION, HIV DEATHS 1999-2013 ### ### ### ### ### ###
V  <- M[61:64,111:121]
V2 <- cbind(V[,1]+V[,2],V[,3:7],rowSums(V[,8:11]))
V3 <- colSums(V2)*1e6

plot(0,0,ylim=c(0,30000),xlim=c(0.6,7.4),xlab="",ylab="",axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:7,paste(c("0-15","15-24","25-34","35-44","45-54","55-64","65+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

for(i in 1:7) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")
points(1:7,colSums(hv_deaths_age),pch=19,cex=1.2,col="black")

mtext("Age Group",1,2.5,cex=0.9)
mtext("Total HIV Deaths by Age Group, 2010-2013",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### HIV Prevalence 2006-2012 ### ### ### ### ### ###
v21  <- rowSums(M[50:66,26:29])

plot(0,0,ylim=c(0,1.8),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(2006:2012,hiv_prev_by_year/1e6,pch=19,cex=0.6,col=1)
lines(2006:2012,hiv_prev_by_year/1e6,lty=3,col=1)
for(i in 1:7) lines(rep(2005+i,2),hiv_prev_by_year[i]/1e6*c(0.8,1.2),lty=1,col=1)
lines(1999:2015,v21,lwd=2,col="blue")

mtext("Year",1,2.5,cex=0.9)
mtext("Total HIV Population (mil)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### ART VOL 1993-2011  ### ### ### ### ### ### 
v23  <- M[43:66,27]+M[43:66,29]
plot(0,0,ylim=c(0,800),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
points(2010,art_vol_10/1e3,pch=19,cex=0.8,col="black")
lines(1992:2015,v23*1e3,lwd=2,col="blue")
lines(rep(2010,2),art_vol_10/1e3*c(.8,1.2),lwd=1,col="black")

mtext("Year",1,2.5,cex=0.9)
mtext("ART Volume (000s)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,"blue"),bg="white",pt.cex=0.8)

##### HIV survival without treatment
pp  <- list(vHMort=vHMort,vHxtoHy=vHxtoHy,mubt00=mubt[12*50,])
sv <- (1/(pp[[1]][,2]+pp[[2]][,1]+pp[[3]])+pp[[2]][,1]/(pp[[1]][,2]+pp[[2]][,1]+pp[[3]])*1/(pp[[1]][,4]+pp[[3]]))[3:8]/12

plot(0,0,xlim=c(0.5,6.5),ylim=c(0,15),axes=F,xlab="",ylab="");box()
axis(1,1:6,paste(c("15-24","25-34","35-44","45-54","55-64","65+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
points(1:6,hiv_surv[,2],pch=19,cex=0.5,col=1)
for(i in 1:6) lines(c(i,i),hiv_surv[i,2]+(hiv_surv[i,3:4]-hiv_surv[i,2])*2,lwd=1)
points(1:6,sv,pch=19,col=2,cex=.5)
text(1:6,sv,format(round(sv,1),nsmall=1),col=2,cex=.9,pos=4)
text(1:6,hiv_surv[,2],format(round(hiv_surv[,2],1),nsmall=1),col=1,cex=.9,pos=4)
axis(2,las=1)
mtext("Years",2,2,cex=0.9)
mtext("Age Group",1,2.5,cex=0.9)
mtext("HIV Untreated Survival from Seroconversion",3,0.8,cex=0.8,font=2)
legend("topright",c("Study data","Fitted model"),pch=19,lwd=c(1,NA),pt.cex=0.5,col=1:2)

### ### ### VARIOUS DEATHS DATA ### ### ### ### ### ### 
DdatDx <- read.csv("DeadAtDx_2-5-16.csv")
load("CalibDat_4-1-16.rData")
DdRpt <- read.csv("tb_deaths_fromreport_2-6-16.csv")
Ddx <- rowSums(DdatDx[7:20,-1]) # dead at dx
Dtx <- CalibDat$tx_outcomes[7:20,3]* CalibDat$tx_outcomes[7:20,4]

#######
plot(1999:2012,rowSums(CalibDat$tb_deaths[1:14,-1]),
     las=1,ylim=c(0,2300),type="l",lwd=2,col=6,xlab="",ylab="Annual deaths")
abline(h=axTicks(2),col="grey85")
lines(1999:2012,Ddx+Dtx,las=1,ylim=c(0,1500),type="l",lwd=2,col=4)
lines(1999:2012,Dtx,las=1,ylim=c(0,1500),type="l",lwd=2,col=3)
lines(1999:2012,rowSums(CalibDat$tb_deaths[1:14,-1]),las=1,ylim=c(0,1500),type="l",lwd=2,col=2)
mtext("Year",1,2.5,cex=0.9)
mtext("Various Data Sources on TB Deaths",3,.8,font=2,cex=0.8)
legend("topright",c("Dead at TB diagnosis or while on TB treatment","Died while on TB treatment","ICD-10 MCD TB codes"),
       col=c(4,3,2),lwd=3,bg="white")


####################################
## Borgdo
n_yr_Z = 50
p2 <- pfast*pimmed*c(1,rep(0,n_yr_Z)) + 
  pfast*(1-pimmed)*(1-(1-exp(-(rfast+rRecov)*12*(0:n_yr_Z)))*(rfast/(rfast+rRecov))) + 
  (1-pfast)*(1-(1-exp(-(rslow+rRecov)*12*(0:n_yr_Z)))*(rslow/(rslow+rRecov)))   
r2 <- -log(1-diff(-p2)/p2[-(n_yr_Z+1)])/1

##
plot(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,type="l",lwd=2,xlim=c(0,15),ylim=c(0,60),las=1,xlab="",ylab="")
lines(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,lwd=2)
lines(1:nrow(datS)-.5,datS[,3]/datS[,2]*1000,lwd=2,lty=2)
lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2)

mtext("Years since infection",1,2.5,cex=0.9)
mtext("Reactivation rate per 1000 PY",2,3,cex=0.8)
legend("topright",c("Ferebee 1970","Sutherland 1968","Fitted model"),
       col=c(1,1,4),lty=c(1,2,1),lwd=2,cex=0.9)
mtext("Reactivation rate vs. Ferebee, Sutherland estimates",3,0.8,cex=0.8,font=2)

## Plot 2  
plot(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,type="l",lwd=2,xlim=c(0,20),ylim=c(0.2,55),las=1,xlab="",ylab="",log="y")
lines(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,lwd=2)
lines(1:nrow(datS)-.5,datS[,3]/datS[,2]*1000,lwd=2,lty=2)
lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2)

mtext("Years since infection",1,2.5,cex=0.9)
mtext("Reactivation rate per 1000 PY (log scale)",2,3,cex=0.8)
legend("topright",c("Ferebee 1970","Sutherland 1968","Fitted model"),
       col=c(1,1,4),lty=c(1,2,1),lwd=2,cex=0.9)
mtext("Reactivation rate vs Ferebee & Sutherland estimates, log-scale",3,0.8,cex=0.8,font=2)

## Plot 3
p <- pfast*pimmed*c(1,rep(0,nrow(datB)-1)) + 
  pfast*(1-pimmed)*(1-(1-exp(-(rfast+rRecov)*12*datB[,1]))*(rfast/(rfast+rRecov))) + 
  (1-pfast)*(1-(1-exp(-(rslow+rRecov)*12*datB[,1]))*(rslow/(rslow+rRecov))) 
p <- 1-(1-p)/(1-p)[nrow(datB)]

plot(datB[,1:2],type="l",lwd=1,lty=3,las=1,xlab="",ylab="")
points(datB[,1:2],cex=0.8,pch=19)
lines(datB[,1],p,col=4,lwd=2)

mtext("Years since infection",1,2.5,cex=0.9)
mtext("Fraction without active TB",2,3,cex=0.8)
legend("topright",c("Borgdorff 2011","Fitted model"), col=c(1,4),lwd=3,cex=0.9)
mtext("Fraction without active TB (of those progressing <15yrs)",3,0.8,cex=0.8,font=2)

## Plot 4

vec <- rep(0,4); names(vec) <- c("safe","slow","fast","prog")
survM15 <- NULL
for(i in 1:(50*12)) {
  if(i==1) {  vec <- c(0,1-Mpfast[3,1],Mpfast[3,1]*(1-pimmed),Mpfast[3,1]*pimmed) } else {
    vec[1] <- vec[1]+vec[2]*rRecov
    vec[4] <- vec[4]+vec[3]*rfast+vec[2]*rep(Mrslow[3:11,1],each=10*12)[i]
    vec[2] <- vec[2]-vec[2]*rRecov-vec[2]*rep(Mrslow[3:11,1],each=10*12)[i]
    vec[3] <- vec[3]-vec[3]*rfast  }
  survM15[i] <- 1- vec[4]  }

  plot(-1,-1,col=0,las=1,ylim=c(0,10),xlim=c(0,50),xlab="",ylab="")
  abline(h=axTicks(2),col="grey90");box()
  lines(seq(0,50*12)/12,c(0,(1-survM15)*100),col=4,lwd=2)
  mtext("Years since infection",1,2.5,cex=0.9)
  mtext("Cumulative incidence of active TB following M.tb infection at age 15 (%)",3,0.8,cex=0.8,font=2)

##############################################
library(RColorBrewer)

cls <-colorRampPalette( brewer.pal(11,"Spectral")[-(4:7)])(10)

#cls <- colorRampPalette(c("red","blue"))(10)
#par(mfrow=c(2,2),mar=c(3,0,0,0),oma=c(.5,5,3,1))

plot(1,1,xlim=c(1990,2017),ylim=range(notif_age_us),log="y",axes=F,xlab="",ylab="")
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(1);axis(2,las=1); box()
for(i in 1:10) lines(1993:2014,notif_age_us[,i],col=cls[i],lwd=2)
for(i in 1:10) { text(c(1993,2014),notif_age_us[c(1,22),i],
                      paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"yrs")[i],col=cls[i],cex=0.7,pos=c(2,4)) }
mtext("Total US Cases, Age Groups Plotted By Year, REPORTED",3,0.8,font=2,cex=0.8)
####
num <- (M[,235:245]+M[,246:256])[44:65,]
num[,10] <- num[,10]+num[,11]; num <- num[,-11] 


plot(1,1,xlim=c(1990,2017),ylim=range(notif_age_us),log="y",axes=F,xlab="",ylab="")
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(1);axis(2,las=1); box()

for(i in 1:10) lines(1993:2014,(num*1e6)[,i],col=cls[i],lwd=2)
for(i in 1:10) { text(c(1993,2014),(num*1e6)[c(1,22),i],
                      paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"yrs")[i],col=cls[i],cex=0.7,pos=c(2,4)) }

mtext("Total US Cases, Age Groups Plotted By Year, MODEL",3,0.8,font=2,cex=0.8)

####################

cls <-colorRampPalette( brewer.pal(11,"Spectral")[-(4:7)])(22)

plot(1,1,xlim=c(0.2,10.8),ylim=range(notif_age_us),las=1,log="y",axes=F,xlab="",ylab="")
axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
     tick=F,cex.axis=0.7)
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(2,las=1); box()

for(i in 1:22) lines(1:10,notif_age_us[i,],col=cls[i],lwd=2)
for(i in 1:22) text(c(1,10),notif_age_us[i,c(1,10)],(1993:2014)[i],col=cls[i],cex=0.7,pos=c(2,4))

mtext("Total US Cases, Years Plotted By Age Group, REPORTED",3,0.8,font=2,cex=0.8)

####

plot(1,1,xlim=c(0.2,10.8),ylim=range(notif_age_us),las=1,log="y",axes=F,xlab="",ylab="")
axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
     tick=F,cex.axis=0.7)
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(2,las=1); box()

for(i in 1:22) lines(1:10,(num*1e6)[i,],col=cls[i],lwd=2)
for(i in 1:22) text(c(1,10),(num*1e6)[i,c(1,10)],(1993:2014)[i],col=cls[i],cex=0.7,pos=c(2,4))
mtext("Total US Cases, Years Plotted By Age Group, MODEL",3,0.8,font=2,cex=0.8)

##############################################
##############################################
##############################################

cls <-colorRampPalette( brewer.pal(11,"Spectral")[-(4:7)])(10)

#cls <- colorRampPalette(c("red","blue"))(10)
#par(mfrow=c(2,2),mar=c(3,0,0,0),oma=c(.5,5,3,1))

plot(1,1,xlim=c(1990,2017),ylim=range(notif_age_fb),log="y",axes=F,xlab="",ylab="")
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(1);axis(2,las=1); box()
for(i in 1:10) lines(1993:2014,notif_age_fb[,i],col=cls[i],lwd=2)
for(i in 1:10) { text(c(1993,2014),notif_age_fb[c(1,22),i],
                      paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"yrs")[i],col=cls[i],cex=0.7,pos=c(2,4)) }
mtext("Total FB Cases, Age Groups Plotted By Year, REPORTED",3,0.8,font=2,cex=0.8)
####
num <- (M[,137:147]-M[,235:245]+M[,216:226]-M[,246:256])[44:65,]
num[,10] <- num[,10]+num[,11]; num <- num[,-11]

plot(1,1,xlim=c(1990,2017),ylim=range(notif_age_fb),log="y",axes=F,xlab="",ylab="")
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(1);axis(2,las=1); box()

for(i in 1:10) lines(1993:2014,(num*1e6)[,i],col=cls[i],lwd=2)
for(i in 1:10) { text(c(1993,2014),(num*1e6)[c(1,22),i],
                      paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"yrs")[i],col=cls[i],cex=0.7,pos=c(2,4)) }

mtext("Total FB Cases, Age Groups Plotted By Year, MODEL",3,0.8,font=2,cex=0.8)

####################

cls <-colorRampPalette( brewer.pal(11,"Spectral")[-(4:7)])(22)

plot(1,1,xlim=c(0.2,10.8),ylim=range(notif_age_fb),las=1,log="y",axes=F,xlab="",ylab="")
axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
     tick=F,cex.axis=0.7)
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(2,las=1); box()

for(i in 1:22) lines(1:10,notif_age_fb[i,],col=cls[i],lwd=2)
for(i in 1:22) text(c(1,10),notif_age_fb[i,c(1,10)],(1993:2014)[i],col=cls[i],cex=0.7,pos=c(2,4))

mtext("Total FB Cases, Years Plotted By Age Group, REPORTED",3,0.8,font=2,cex=0.8)

####

plot(1,1,xlim=c(0.2,10.8),ylim=range(notif_age_fb),las=1,log="y",axes=F,xlab="",ylab="")
axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
     tick=F,cex.axis=0.7)
abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
axis(2,las=1); box()

for(i in 1:22) lines(1:10,(num*1e6)[i,],col=cls[i],lwd=2)
for(i in 1:22) text(c(1,10),(num*1e6)[i,c(1,10)],(1993:2014)[i],col=cls[i],cex=0.7,pos=c(2,4))
mtext("Total FB Cases, Years Plotted By Age Group, MODEL",3,0.8,font=2,cex=0.8)

####################################
#################################### 
### Recent infection
#colnames(M)
Vall <- (M[65,198:214]/M[65,181:197])
plot(-1,0,ylim=c(0.02,0.8),xlim=c(0.5,17.5),xlab="",ylab="",axes=F)
axis(2,las=2);box()
axis(1,1:17,c("All",paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),
                          "yrs"),"US born","Foreign born","FB >2yrs","Homeless","HIV positive"),tick=F,cex.axis=0.7,las=2,
     mgp=c(3, 0.25, 0))

abline(h=axTicks(2),col="grey85")
mtext("Fraction of Incident TB from Recent Infection (<2 years)",3,.8,font=2,cex=0.8)

for(i in 1:17) lines(rep(i,2),c(0,Vall[i]),col="forestgreen",lwd=10,lend="butt")
text(1:17,Vall,format(round(Vall,2),nsmall=2),cex=0.7,pos=3)

####################################
####################################

dev.off(); system(paste("open", pdfnam)) # code for ma

