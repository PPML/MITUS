
load("data/results 2018-04-12 09:17:22 .rData")
M <- results
pdf(file=paste("MITUS_results/orig_graphs",Sys.time(),".pdf"), width = 11, height = 8.5)

### ### ### ### ### ### TOTAL DIAGNOSED CASES 1953-2013  ### ### ### ### ### ###
V0   <- M[4:66,"NOTIF_ALL"]+M[4:66,"NOTIF_MORT_ALL"]
### ### ### ### ### US notif by age + US dead @ diag by age ### ### ### ### ###
V1   <- rowSums(M[44:66,205:215]+M[44:66,216:226])
### ### ### ### ### ### notif by age + dead @ diag by age    ### ### ### ### ###
V2   <- rowSums((M[44:66,136:146]+M[44:66,189:199]) - (M[44:66,205:215]+M[44:66,216:226]))

plot(0,0,ylim=c(0,max(range(V0)))*1e3,xlim=c(1954,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(1953:2015,V0*1e3,lwd=3,col="white"); lines(1953:2015,V0*1e3,lwd=2,col=1)
lines(1993:2015,V1*1e3,lwd=3,col="white"); lines(1993:2015,V1*1e3,lwd=2,col=4)
lines(1993:2015,V2*1e3,lwd=3,col="white"); lines(1993:2015,V2*1e3,lwd=2,col=3)

mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Cases Identified (000s), 1953-2014",3,.8,font=2,cex=0.8)

legend("topright",c("Fitted model (all)","Fitted model (US born)","Fitted model (foreign born)"),
       pch=c(NA,NA,NA),lwd=c(2,2,2),lty=c(1,1,1),col=c(1,4,3),bg="white",ncol=2,cex=.8,pt.cex=0.4)

### ### ### ### ### ###  CASES FB DISTRIBUTION 1993-2013  ### ### ### ### ### ###
### ### ### dx FB +  dead @ dx  FB, dx US + dx dead @ dx US
V   <- cbind(M[44:66,148]+M[44:66,149]+M[44:66,201]+M[44:66,202], M[44:66,147]+M[44:66,200])
V <- V[,1]/rowSums(V )
plot(0,0,ylim=c(2.5,97.5),xlim=c(2000,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(1993:2015,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Foreign-Born, 2000-14",3,.8,font=2,cex=0.8)
legend("topleft","Fitted model",pch=NA,lwd=2,col=4,lty=1,bg="white",pt.cex=0.6)

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ###
V   <- cbind(M[44:65,148]+M[44:65,201],M[44:65,149]+M[44:65,202])
V <- V[,1]/rowSums(V)

plot(0,0,ylim=c(0,60),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(1993:2014,V*100,lwd=2,col=4)

mtext("Year",1,2.5,cex=0.9)
mtext("Percent of FB Cases Arrived in Past 2 Yrs",3,.8,font=2,cex=0.8)
legend("topright","Fitted model",pch=19,lwd=2,col=4,lty=1,bg="white",pt.cex=0.6)

### ### ### ### ### ### CASES AGE DISTRIBUTION 2000-2014  ### ### ### ### ### ###
V   <- (M[51:65,136:146]+M[51:65,189:199])
V2  <- V[,-11]
V2[,10] <- V2[,10]+V[,11]
cls <- colorRampPalette(c("blue", "red"))( 4 )
plot(0,0,ylim=c(0,1100),xlim=c(2000,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")

lines(2000:2014,rowSums(V2[,1:3])*1e3,lwd=2,col=cls[1])
lines(2000:2014,rowSums(V2[,4:5])*1e3,lwd=2,col=cls[2])
lines(2000:2014,rowSums(V2[,6:7])*1e3,lwd=2,col=cls[3])
lines(2000:2014,rowSums(V2[,8:10])*1e3,lwd=2,col=cls[4])

mtext("TB Cases By Age (000s), 2000-14",3,.8,font=2,cex=0.8)
mtext("Year",1,2.5,cex=0.9)

legend("topright",c("0-24 years","25-44 years","45-64 years","65+ years"),
       lwd=NA,lty=NA,col=cls,bg="white",
       pt.cex=c(1.8,1.8,1.8,1.8),pch=c(15,15,15,15))

### ### ### ### ### ### CASES AGE DISTRIBUTION 2000-2014  ### ### ### ### ### ###
V   <- (M[51:65,136:146]+M[51:65,189:199])
V2  <- V[,-11]
V2[,10] <- V2[,10]+V[,11]
V2  <- colSums(V2)/sum(V2)*100

plot(0,0,ylim=c(0,20),xlim=c(0.6,10.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:10,paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"\nyears",sep=""),
     tick=F,cex.axis=0.6)
axis(1,1:11-0.5,rep("",11))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V2[i],V2[i]),border="white",col="lightblue")
#points(1:10,colSums(notif_age[7:21,])/sum(notif_age[7:21,])*100,pch=19,cex=1.2)
mtext("Age Group",1,2.5,cex=0.9)
mtext("Age Distribution of TB Cases (%), 2000-14",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### ### ### ### CASES HR DISTRIBUTION 1993-2013  ### ### ### ### ### ###
V   <- cbind(M[44:65,"NOTIF_HR"],M[44:65,"NOTIF_LR"]) + cbind(M[44:65,"NOTIF_MORT_HR"],M[44:65,"NOTIF_MORT_LR"])
V <- V[,1]/rowSums(V)
plot(0,0,ylim=c(0,15),xlim=c(1993,2014),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
# points(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,pch=19,cex=0.6)
# lines(1993:2014,notif_us_hr[,1]/rowSums(notif_us_hr)*100,lty=3)
lines(1993:2014,V*100,lwd=2,col=4)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of TB Cases Homeless in Past Yr",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### ### ### ### TREATMENT OUTCOMES 1993-2011  ### ### ### ### ### ###
V   <- M[44:63,132:134]
Vdisc <- V[,2]/rowSums(V)
Vdead <- V[,3]/rowSums(V)

plot(0,0,ylim=c(0,12.5),xlim=c(1993,2011),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
# points(1993:2012,tx_outcomes[,2]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="red3")
# points(1993:2012,tx_outcomes[,3]/rowSums(tx_outcomes)*100,pch=19,cex=0.6,col="blue")
# lines(1993:2012,tx_outcomes[,2]/rowSums(tx_outcomes)*100,lty=3,col="red3")
# lines(1993:2012,tx_outcomes[,3]/rowSums(tx_outcomes)*100,lty=3,col="blue")
lines(1993:2012,Vdisc*100,lwd=2,col="red3")
lines(1993:2012,Vdead*100,lwd=2,col="blue")
mtext("Year",1,2.5,cex=0.9)
mtext("Treatment Outcomes: Discontinued and Died (%)",3,.8,font=2,cex=0.8)
legend("topright",c("Discontinued","Died","Reported data","Fitted model"),pch=c(15,15,19,NA),lwd=c(NA,NA,1,2),
       col=c("red3",4,1,1),lty=c(NA,NA,3,1),bg="white",pt.cex=c(1.8,1.8,0.6,NA))

### ### ### ### ### ### TLTBI VOL 1993-2011  ### ### ### ### ### ###
v12  <- M[43:65,152]

plot(0,0,ylim=c(0,500),xlim=c(1992,2015),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
# points(2002,tltbi_vol[1]/1e3,pch=19,cex=0.8,col="black")
lines(1992:2014,v12*1e3,lwd=2,col="blue")
# lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

mtext("Year",1,2.5,cex=0.9)
mtext("IPT Treatment Initiations Per Year (000s)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(2,2),col=c(1,"blue"),bg="white",pt.cex=0.8)

### ### ### ### ### ### LTBI INITIATIONS 1993-2011 Distribution ### ### ### ### ### ###
v13  <- cbind(M[43:65,"TLTBI_INITS_FB"],M[43:65,"TLTBI_INITS_HR"])
v13  <-   v13/M[43:65,"TLTBI_INITS"]
plot(1,1,ylim=c(0.001,.5)*100,xlim=c(1992,2015),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
# points(rep(2002,3),TLTBI_dist*100,pch=19,cex=0.8,col=c("red3",4,6))
for(i in 1:3) lines(1992:2014,v13[,i]*100,lwd=2,col=c("red3",4,6)[i])
# lines(rep(2002,2),tltbi_vol[2:3]/1e3,lwd=2,col="black")

mtext("Year",1,2.5,cex=0.9)
mtext("IPT Treatment Initiations By Risk Group (%)",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Foreign-born","Homeless","Reported data","Fitted model"),
       pch=c(15,15,19,NA),lwd=c(NA,NA,NA,2),col=c("red3",4,1,1),bg="white",pt.cex=c(1.8,1.8,0.8,NA))

## ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ###
## plot not functioning
V  <- cbind(M[62,55:65], M[62,33:43]-M[62,55:65])
V[9,] <- colSums(V[9:11,])
V <- V[2:9,1]/rowSums(V[2:9,])*100
plot(0,0,ylim=c(0,100),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:8-0.5,rep("",8))
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V[i],V[i]),border="white",col="lightblue")
# points(1:8,ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3])*100,pch=19,cex=1.2)
# for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_us_11[i,2],ltbi_us_11[i,3])*100,pch=19,cex=1.2)
mtext("Age Group",1,2.5,cex=0.9)
mtext("LTBI in US Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=c(0,NA),
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ###
### plot not functioning NaN with division line
# V  <- cbind(M[62,65:75],M[62,43:53]-M[62,65:75])
# V
# V[9,] <- colSums(V[9:11,])
# V <- V[2:9,1]/rowSums(V[2:9,])*100
#
# plot(0,0,ylim=c(0,52),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
# axis(1,1:8,paste(c(paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
# axis(1,1:8-0.5,rep("",8))
# axis(2,las=2);box()
# abline(h=axTicks(2),col="grey85")
# for(i in 1:8) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V[i],V[i]),border="white",col="lightblue")
# # points(1:8,ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3])*100,pch=19,cex=1.2)
# # for(i in 1:8) lines((1:8)[c(i,i)],qbeta(c(1,39)/40,ltbi_fb_11[i,2],ltbi_fb_11[i,3])*100,pch=19,cex=1.2)
# mtext("Age Group",1,2.5,cex=0.9)
# mtext("LTBI in Foreign Born Population 2011 by Age (%)",3,.8,font=2,cex=0.8)
# legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=c(0,NA),
#        pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### ### ### ### TOTAL POP EACH DECADE, BY US/FB  ### ### ### ### ### ###
V  <- cbind(M[1:66,30], M[1:66,31]+M[1:66,32])
plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
# points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
# points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
# points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
# lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],lty=3,col="grey50")
# lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],lty=3,col="blue")
# lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],lty=3,col="red3")
lines(1950:2015,V[,2],lwd=2,col="red3")
lines(1950:2015,V[,1],lwd=2,col="blue")
lines(1950:2015,rowSums(V),lwd=2,col="grey50")

mtext("Year",1,2.5,cex=0.9)
mtext("Population: Total, US, and Foreign Born (mil, log-scale)",3,.8,font=2,cex=0.8)
legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
       pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
### need to update
V  <- cbind(t(M[65,33:43]), t(M[65,44:54]))

plot(0,1,ylim=c(0.05,135),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y")
axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
axis(1,1:9-0.5,rep("",9))
axis(2,c(0.1,1,10,100),las=2);box()
abline(h=axTicks(2),col="grey85")
for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")
# points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
# points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")
mtext("Age Group",1,2.5,cex=0.9)
box()
mtext("Population by Age for FB (red) and US (blue), 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
legend("topright",c("Reported data","Fitted model"),pch=c(19,15),pt.cex=c(1,2),
       lwd=NA,col=c("grey30","grey80"),bg="white")

### ### ### TB DEATHS 1999-2013 ### ### ### ### ### ###
###   This is a logically invlaid value
V  <- M[50:65,88:98]+M[50:65,99:109]
V2 <- V[,-11];
V2[,10] <- V[,10]+V[,11]

plot(0,0,ylim=c(0,250000),xlim=c(1998.5,2014.5),xlab="",ylab="",axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
#points(1999:2014,rowSums(tb_deaths),pch=19,cex=0.6,col=1)
#lines(1999:2014,rowSums(tb_deaths),lty=3,col=1)
lines(1999:2014,rowSums(V2)*1e6,lwd=2,col="blue")

mtext("Year",1,2.5,cex=0.9)
mtext("Total TB Deaths by Year",3,.8,font=2,cex=0.8)
legend("bottomleft",c("Reported data","Fitted model"),pch=c(19,NA),lwd=c(1,2),col=c(1,4),lty=c(3,1),bg="white",pt.cex=0.6)

### ### ### ### ### ### AGE DISTRIBUTION, TB DEATHS 1999-2013 ### ### ### ### ### ###
### ### ### values are off by two orders of magnitude
V  <- M[50:65,88:98]+M[50:65,99:109]
V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
V3 <- colSums(V2)*1e6

plot(0,0,ylim=c(0,500000),xlim=c(0.6,10.4),xlab="",ylab="",axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:10,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)

for(i in 1:10) polygon(i+c(-.5,.5,.5,-.5),c(0,0,V3[i],V3[i]),border="white",col="lightblue")
# points(1:10,colSums(tb_deaths),pch=19,cex=1.2,col="black")

mtext("Age Group",1,2.5,cex=0.9)
mtext("Total TB Deaths by Age Group 1999-2014",3,.8,font=2,cex=0.8)
legend("topleft",c("Reported data","Fitted model"),pch=c(19,15),lwd=NA,
       pt.cex=c(1,2),col=c("black","lightblue"),bg="white")

### ### ### ### ### ### Total US Cases, Age Groups Plotted By Year, MODEL
# num <- (M[,205:215]+M[,216:226])[44:65,]
# num[,10] <- num[,10]+num[,11]; num <- num[,-11]
#
#
# plot(1,1,xlim=c(1990,2017),ylim=range(num),log="y",axes=F,xlab="",ylab="")
# abline(h=c(1:19/2,1:19*5,1:19*50,1:19*500),col="grey95");abline(h=c(100,1000),col="grey80"); box()
# axis(1);axis(2,las=1); box()
#
# for(i in 1:10) lines(1993:2014,(num*1e6)[,i],col=cls[i],lwd=2)
# for(i in 1:10) { text(c(1993,2014),(num*1e6)[c(1,22),i],
#                       paste(c("0-4",paste(0:7*10+5,1:8*10+4,sep="-"),"85+"),"yrs")[i],col=cls[i],cex=0.7,pos=c(2,4)) }
#
# mtext("Total US Cases, Age Groups Plotted By Year, MODEL",3,0.8,font=2,cex=0.8)

dev.off()
