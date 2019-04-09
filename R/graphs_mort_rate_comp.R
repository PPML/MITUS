# #MORTALITY RATE INPUT COMPARISON PLOTS
#
# #READ IN THE NCHS DATA
# NCHS_mort<-readRDS(system.file("US/US_NCHS_mort.rds", package="MITUS"))
#
# #SET THE COLORS
# col<-rainbow(11)
# pdfname<-paste("MITUS_results/mort_rate_comp",Sys.time(),".pdf")
# pdf(file=pdfname, width = 11, height = 8.5)
# par(mfrow=c(2,2),mar=c(4,4.5,3,1))
#
# for (i in 2:12){
#   plot(0,0,ylim=c(min(NCHS_mort[,i],Inputs$BgMort[1:68,i])*.5,max(NCHS_mort[,i],Inputs$BgMort[1:68,i])*1.5),xlim=c(1950,2017),xlab="",ylab="",axes=F)
#   axis(1);axis(2,las=2);box()
#   abline(h=axTicks(2),col="grey85")
#   lines(1950:2017, NCHS_mort[,i], lwd=3, col=col[i-1])
#   lines(1950:2017, Inputs$BgMort[1:68,i], lty=3, lwd=2, col=col[i-1])
#   mtext("Year",1,2.5,cex=1.2)
#   mtext(paste("Age Specific Mortality Rates from 1950 to 2017", colnames(Inputs$BgMort)[i]),3,.8,font=2,cex=1.2)
#   legend("topright",c("NCHS","lifetables"),cex=1.0,
#   pch=c(NA,NA),lwd=c(3,3),lty=c(1,3),col=col[i-1],bg="white",pt.cex=c(1.8,1.8))
# }
#
# dev.off()
#
