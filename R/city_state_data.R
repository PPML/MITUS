# #format states and city report tx data
# pdf(file=paste("MITUS_results/city_state_data",Sys.time(),".pdf"), width = 11, height = 8.5)
# par(mfrow=c(2,2),mar=c(4,4.5,3,1))
# x<-(0:150)
# #recent infections active treatment completion in 12 months or less
# cs_txcomp<-read.csv("~/MITUS/inst/extdata/city_state_txcomp.csv")
# cs_txcomp[,1]<-stateID[,3]
# x1<-(0:max(cs_txcomp[,4]))
#
# boxplot(cs_txcomp[,2], horizontal = TRUE, xlab="Percent Active Tx Completion")
#
# plot(cs_txcomp[,3]~cs_txcomp[,4],xlab="total ActTB diagnoses",ylab="total ActTX completions")
# lines(x=x1,y=x1,col="red")
# legend("bottomright", c("reported % Active TX initiation","100%  Active TX initiation"),cex=0.9,
#        pch=c(1,NA),lwd=c(NA,1),lty=c(NA,1),col=c(1,"red"),bg="white",pt.cex=c(1.8,NA))
#
# plot(cs_txcomp[,3]~cs_txcomp[,4],xlim=c(0,150),ylim=c(0,150),xlab="total ActTB diagnoses",ylab="total ActTX completions")
# lines(x=x,y=x,col="red")
# legend("bottomright", c("reported % Active TX initiation","100%  Active TX initiation"),cex=0.9,
#        pch=c(1,NA),lwd=c(NA,1),lty=c(NA,1),col=c(1,"red"),bg="white",pt.cex=c(1.8,NA))
#
# #recent LTBI diagnosed contacts begin TLTBI
# cs_tltbi_init<-read.csv("~/MITUS/inst/extdata/city_state_tltbi_init.csv")
# cs_tltbi_init[,1]<-stateID[,3]
# x1<-(0:max(cs_tltbi_init[,4]))
#
# boxplot(cs_tltbi_init[,2], horizontal = TRUE, xlab="Percent TLTBI Initiation")
# plot(cs_tltbi_init[,3]~cs_tltbi_init[,4],xlab="total LTBI diagnoses",ylab="total TLTBI initiates")
# lines(x=x1,y=x1,col="red")
# legend("bottomright", c("reported % TLTBI initiation","100% TLTBI initiation"),cex=0.9,
#        pch=c(1,NA),lwd=c(NA,1),lty=c(NA,1),col=c(1,"red"),bg="white",pt.cex=c(1.8,NA))
#
# plot(cs_tltbi_init[,3]~cs_tltbi_init[,4],xlim=c(0,150),ylim=c(0,150),xlab="total LTBI diagnoses",ylab="total TLTBI initiates")
# lines(x=x,y=x,col="red")
# legend("bottomright", c("reported % TLTBI initiation","100% TLTBI initiation"),cex=0.9,
#        pch=c(1,NA),lwd=c(NA,1),lty=c(NA,1),col=c(1,"red"),bg="white",pt.cex=c(1.8,NA))
# #recent LTBI diagnosed contacts complete TLTBI
# cs_tltbi_comp<-read.csv("~/MITUS/inst/extdata/city_state_tltbi_comp.csv")
# cs_tltbi_comp[,1]<-stateID[,3]
# x1<-(0:max(cs_tltbi_comp[,4]))
#
# boxplot(cs_tltbi_comp[,2], horizontal = TRUE, xlab="Percent TLTBI Completion")
# plot(cs_tltbi_comp[,3]~cs_tltbi_comp[,4],xlab="total TLTBI initiates",ylab="total TLTBI completions")
# lines(x=x1,y=x1,col="red")
# legend("bottomright", c("reported % TLTBI completion","100% TLTBI completion"),cex=0.9,
#        pch=c(1,NA),lwd=c(NA,1),lty=c(NA,1),col=c(1,"red"),bg="white",pt.cex=c(1.8,NA))
#
# plot(cs_tltbi_comp[,3]~cs_tltbi_comp[,4],xlim=c(0,150),ylim=c(0,150),xlab="total TLTBI initiates",ylab="total TLTBI completions")
# lines(x=x,y=x,col="red")
# legend("bottomright", c("reported % TLTBI completion","100% TLTBI completion"),cex=0.9,
#        pch=c(1,NA),lwd=c(NA,1),lty=c(NA,1),col=c(1,"red"),bg="white",pt.cex=c(1.8,NA))
# # with(cs_tltbi_comp[,], text(cs_tltbi_comp[,3]~cs_tltbi_comp[,4], labels = cs_tltbi_comp[,1]), pos = 1,)
#
# # ggplot(data = cs_txcomp, aes(cs_txcomp[,2]))+geom_violin()
#
# dev.off()
