graphs_pub<-function(Par_list, loc="US"){
  ##get data from the calibdat
  if (loc=="US"){
  datF         <- CalibDat[["ferebee_data"]]
  datB         <- CalibDat[["borgdorff_data"]]
  datS         <- CalibDat[["sutherland_data"]]
  } else {
    datF         <- CalibDatState[["ferebee_data"]]
    datB         <- CalibDatState[["borgdorff_data"]]
    datS         <- CalibDatState[["sutherland_data"]]
  }
  datSz        <- datS; datSz[datS[,3]==0,3] <- 0.01
####################################
## Borgdo
##updated based on the calculations for the llikelihood function 2/5/19

n_yr_Z = 50

  Mpfast     <- Par_list[[1]];
  # ORpfastRF <- Par[2];
  Mrslow     <- Par_list[[2]]*12;
  # RRrSlowRF <- Par[4];
  rfast     <- Par_list[[3]]*12;
  rRecov    <- Par_list[[4]]*12;
  pfast_v <- rslow_v <- rep(NA,4)
  pfast_v <- Mpfast[3,]
  rslow_v <- Mrslow[3,]
  p0 <- matrix(NA,4,n_yr_Z+1)  # this added
  for (i in 1:4){
    p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*(0:n_yr_Z)))*(rfast/(rfast+rRecov))) +
      (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*(0:n_yr_Z)))*(rslow_v[i]/(rslow_v[i]+rRecov)))
  }
  p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
  r2 <- -log(1-diff(-p1)/p1[-(n_yr_Z+1)])/1

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

###################################### Plot 2
n_yr_F       <- nrow(datF)
plot(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,type="l",lwd=2,xlim=c(0,20),ylim=c(0.2,55),las=1,xlab="",ylab="",log="y")
lines(1:nrow(datF)-.5,datF[,3]/datF[,2]*1000,lwd=2)
lines(1:nrow(datS)-.5,datS[,3]/datS[,2]*1000,lwd=2,lty=2)
lines(1:n_yr_Z-.5,r2*1000,col=4,lwd=2)

mtext("Years since infection",1,2.5,cex=0.9)
mtext("Reactivation rate per 1000 PY (log scale)",2,3,cex=0.8)
legend("topright",c("Ferebee 1970","Sutherland 1968","Fitted model"),
       col=c(1,1,4),lty=c(1,2,1),lwd=2,cex=0.9)
mtext("Reactivation rate vs Ferebee & Sutherland estimates, log-scale",3,0.8,cex=0.8,font=2)

###################################### Plot 3
##updated based on the calculations for the llikelihood function 2/5/19
Mpfast     <- Par_list[[1]];
# ORpfastRF <- Par[2];
Mrslow     <- Par_list[[2]]*12;
# RRrSlowRF <- Par[4];
rfast     <- Par_list[[3]]*12;
rRecov    <- Par_list[[4]]*12;
pfast_v <- rslow_v <- rep(NA,4)
pfast_v <- Mpfast[3,]
rslow_v <- Mrslow[3,]
p0 <- matrix(NA,4,10)  # this added
for (i in 1:4){
  p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*datB[,1]))*(rfast/(rfast+rRecov))) +
    (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*datB[,1]))*(rslow_v[i]/(rslow_v[i]+rRecov)))
}
# p1 <- as.numeric(t(p0)%*%v21a[17,1:4])

p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
p <- 1-(1-p1)/(1-p1)[nrow(datB)]

plot(datB[,1:2],type="l",lwd=1,lty=3,las=1,xlab="",ylab="")
points(datB[,1:2],cex=0.8,pch=19)
lines(datB[,1],p,col=4,lwd=2)

mtext("Years since infection",1,2.5,cex=0.9)
mtext("Fraction without active TB",2,3,cex=0.8)
legend("topright",c("Borgdorff 2011","Fitted model"), col=c(1,4),lwd=3,cex=0.9)
mtext("Fraction without active TB (of those progressing <15yrs)",3,0.8,cex=0.8,font=2)

#################################### Plot 4

# vec <- rep(0,4); names(vec) <- c("safe","slow","fast","prog")
# survM15 <- NULL
# for(i in 1:(50*12)) {
#   if(i==1) {  vec <- c(0,1-Mpfast[3,1],Mpfast[3,1]*(1-pimmed),Mpfast[3,1]*pimmed) } else {
#     vec[1] <- vec[1]+vec[2]*rRecov
#     vec[4] <- vec[4]+vec[3]*rfast+vec[2]*rep(Mrslow[3:11,1],each=10*12)[i]
#     vec[2] <- vec[2]-vec[2]*rRecov-vec[2]*rep(Mrslow[3:11,1],each=10*12)[i]
#     vec[3] <- vec[3]-vec[3]*rfast  }
#   survM15[i] <- 1- vec[4]  }
#
# plot(-1,-1,col=0,las=1,ylim=c(0,10),xlim=c(0,50),xlab="",ylab="")
# abline(h=axTicks(2),col="grey90");box()
# lines(seq(0,50*12)/12,c(0,(1-survM15)*100),col=4,lwd=2)
# mtext("Years since infection",1,2.5,cex=0.9)
# mtext("Cumulative incidence of active TB following M.tb infection at age 15 (%)",3,0.8,cex=0.8,font=2)
}
