# V  <- cbind(rowSums(df[1:11,255:265]), rowSums(df[1:11,266:276]))
# V1c <- rowSums(df[1:11,121:131])
# plot(1,1,ylim=c(0,2.5),xlim=c(1950,1960),xlab="",ylab="",axes=F)
# axis(1);axis(2,las=2);box()
# abline(h=axTicks(2),col="grey85")
#
# lines(1950:1960,V[,2],lwd=2,col="red3")
# lines(1950:1960,V[,1],lwd=2,col="blue")
# lines(1950:1960,V1c,lwd=2,col="grey50")
# points(CalibDat$US_tot_mort[1:11,1],(CalibDat$US_tot_mort[1:11,2])/1e6,pch=19,cex=0.6,col="grey50")
# lines(CalibDat$US_tot_mort[1:11,1],(CalibDat$US_tot_mort[1:11,2])/1e6,lty=3,col="grey50")
#
# mtext("Year",1,2.5,cex=1.2)
# mtext("Mortality: Total, US, and Non-US Born (mil)",3,.8,font=2,cex=1.2)
# legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=1.0,
#        pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
