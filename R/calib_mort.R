# ###CHECK CALIBRATION TARGET GRAPHS
#
# #IF THE POPULATION TARGET IN 2014 IS RIGHT AND THE MORT RATE IS RIGHT WHAT DOES IT LOOK LIKE
#
# pop<-CalibDat$tot_pop_yr_fb[7,2]*1e6
# age_pop<-CalibDat$tot_pop16_ag_fb[,2]
# mort<-CalibDat$US_tot_mort[61,2]
# age_mort<-CalibDat$US_mort_age[12,2:9]
# mu<-Inputs$BgMort[61,2:12]
#
# V<-age_pop*mubt[793,]
# V3<-V/sum(V)
# V3[,8] <- V3[,8]+V2[,9]
# V3<-V3/rowSums(V3)
#
# plot(0,0,ylim=c(0.05,max(range(V3))+.2),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
# axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
# axis(1,1:9-0.5,rep("",9))
# axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
# abline(h=axTicks(2),col="grey85")
#
# for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[1,i],V3[1,i]),border=NA,col="gray")
# for(i in 1:8) points(i+.2,CalibDat$US_mort_age[16,]/rowSums(CalibDat$US_mort_age[16,]),pch=19,cex=1.2,col="black")
#
#
# mtext("Age Group",1,2.5,cex=1.2)
# box()
# mtext("Mortality by Age, 2014 (mil)",3,.8,font=2,cex=1.2)
# legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
#        lwd=NA,col=c("black","gray"),bg="white")
