###############
library(Hmisc)
library(abind)
setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC Large Grant/Analysis Transmission")

load("parAll_9-14-16.rData") # parAll
wt = parAll[,"wt"]

dim(parAll)
sum(wt)^2/sum(wt^2) # 3013.833
for(rr in 0:14) { load(paste("Cluster Results/MiAllex",rr,"_0_9-21-16.rData",sep="")); print(rr); flush.console() }
# load("AddnOutcomes_V1.rData") # ImpTb16","ImpAct16","ImpAct16Mdr
dim(MiAllex0)
load("AddnOutcomes_V2.rData")

LE <- read.csv("US_LE_2011.csv")[,2]

### Sample Quality weights
##### Sample and insert additional variables 
library(lhs)
set.seed(988)
qol_samp <- randomLHS(nrow(parAll),3)
colnames(qol_samp) <- c("tltbi","activeTB","TbTx")
meanQOL <- c(0.03,0.16,0.10) # GUO 2008 HUI3 medians also Holland 2009, for active disease = Tan 2008
for(i in 1:3)   qol_samp[,i] <- qgamma(qol_samp[,i], shape   = 2^2, rate   = 2^2/meanQOL[i])  

# Helper function to subset out columns and rows
gc <- function(yr,id) {
  out = NULL
  if(length(yr)==1 & length(id)==1) {
    for(bb in 0:14) out = c(out,get(paste("MiAllex",bb,sep=""))[,yr-1949,id])  }
  if( (length(yr)>1 & length(id)==1) | (length(yr)==1 & length(id)>1) ) {
    for(bb in 0:14) out = rbind(out,get(paste("MiAllex",bb,sep=""))[,yr-1949,id])  }
  if( length(yr)>1 & length(id)>1 ) {
    for(bb in 0:14) out = abind(out,get(paste("MiAllex",bb,sep=""))[,yr-1949,id],along=1)   }
  return(out)  }

# works
############################
####  RESULTS BY LOCATION IN MANUSCRIPT ############################
############################
# Foreign born as faction of total pop 2100
  Pop99      <- gc(2099,2) 
  Pop99us    <- gc(2099,30) + gc(2099,31) 
  Pop99fb    <- gc(2099,32) + gc(2099,33)
  wtd.mean(Pop99fb/Pop99,wt)*1e2; wtd.quantile(Pop99fb/Pop99,wt,c(1,39)/40)*1e2

# incidence in 2050, pct reduction 2016, us then fb
  Pop16      <- gc(2016,2) 
  Pop16us    <- rowSums(gc(2016,30:31))  
  Pop16fb    <- rowSums(gc(2016,32:33)) 
  Pop16y     <- rowSums(gc(2016,3:5))
  Pop16o     <- rowSums(gc(2016,10:13))
  
  Pop50      <- gc(2050,2)
  Pop50us    <- rowSums(gc(2050,30:31))  
  Pop50fb    <- rowSums(gc(2050,32:33)) 
  Pop50y     <- rowSums(gc(2050,3:5))
  Pop50o     <- rowSums(gc(2050,10:13))

  Cases16    <- rowSums(gc(2016,c(136,215)))
  Cases16us  <- rowSums(gc(2016,c(152,153,231,232))) 
  Cases16fb  <- rowSums(gc(2016,c(154,155,233,234)))
  Cases16y   <- rowSums(gc(2016,c(137:139,216:218)))
  Cases16o   <- rowSums(gc(2016,c(144:147,223:226)))
  
  Cases50    <- rowSums(gc(2050,c(136,215)))
  Cases50us  <- rowSums(gc(2050,c(152,153,231,232))) 
  Cases50fb  <- rowSums(gc(2050,c(154,155,233,234)))
  Cases50y   <- rowSums(gc(2050,c(137:139,216:218)))
  Cases50o   <- rowSums(gc(2050,c(144:147,223:226)))
  
# USB inc 2050
wtd.mean(Cases50us/Pop50us,wt)*1e6; wtd.quantile(Cases50us/Pop50us,wt,c(1,39)/40)*1e6
## Cases 50 US pct red vs 2016
wtd.mean(1-(Cases50us/Pop50us)/(Cases16us/Pop16us),wt)*1e2; wtd.quantile(1-(Cases50us/Pop50us)/(Cases16us/Pop16us),wt,c(1,39)/40)*1e2
# FB inc 2050
wtd.mean(Cases50fb/Pop50fb,wt)*1e6; wtd.quantile(Cases50fb/Pop50fb,wt,c(1,39)/40)*1e6
## Cases 50 FB pct red vs 2016
wtd.mean(1-(Cases50fb/Pop50fb)/(Cases16fb/Pop16fb),wt)*1e2; wtd.quantile(1-(Cases50fb/Pop50fb)/(Cases16fb/Pop16fb),wt,c(1,39)/40)*1e2
## Cases 50 pct FB
wtd.mean(Cases50fb/Cases50,wt)*1e2; wtd.quantile(Cases50fb/Cases50,wt,c(1,39)/40)*1e2

## Cases 50 young (<25) pct red vs 2016
wtd.mean(1-(Cases50y/Pop50y)/(Cases16y/Pop16y),wt)*1e2; wtd.quantile(1-(Cases50y/Pop50y)/(Cases16y/Pop16y),wt,c(1,39)/40)*1e2
## Cases 50 old (+65) pct red vs 2016
wtd.mean(1-(Cases50o/Pop50o)/(Cases16o/Pop16o),wt)*1e2; wtd.quantile(1-(Cases50o/Pop50o)/(Cases16o/Pop16o),wt,c(1,39)/40)*1e2

### Elimination 2100
Cases99us  <- rowSums(gc(2099,c(152,153,231,232))) 
Cases99fb  <- rowSums(gc(2099,c(154,155,233,234)))
Cases99    <- rowSums(gc(2099,c(136,215)))

Pop99us    <- rowSums(gc(2099,30:31))  
Pop99fb    <- rowSums(gc(2099,32:33)) 
Pop99      <- gc(2099,2)

# Prob US elimi 99
sum((1.0>Cases99us/Pop99us*1e6)*wt)/sum(wt)*100
# Prob US elimi by year
cbind(2050:2099,apply(as.matrix(2050:2099),1,function(x) sum((1.0>rowSums(gc(x,c(152,153,231,232)))/rowSums(gc(x,30:31))*1e6)*wt)/sum(wt) )*100)
# Cases 99
wtd.mean(Cases99/Pop99,wt)*1e6; wtd.quantile(Cases99/Pop99,wt,c(1,39)/40)*1e6
## Cases 99 pct FB
wtd.mean(Cases99fb/Cases99,wt)*1e2; wtd.quantile(Cases99fb/Cases99,wt,c(1,39)/40)*1e2

##### RT
Incid16    <- gc(2016,181)
Rt16       <- gc(2016,198) 

## RT
wtd.mean(Rt16/Incid16,wt)*100; wtd.quantile(Rt16/Incid16,wt,c(1,39)/40)*100

## SOURCES OF LTBI
Inf16      <- rowSums(gc(2016,304:307)) 
Inf16us    <- rowSums(gc(2016,304:305)) 
Inf16fb    <- rowSums(gc(2016,306:307)) 

## INF  fix super infection
wtd.mean(Inf16,wt)*1e3; wtd.quantile(Inf16,wt,c(1,39)/40)*1e3
## INF US
wtd.mean(Inf16us,wt)*1e3; wtd.quantile(Inf16us,wt,c(1,39)/40)*1e3
## INF FB
wtd.mean(Inf16fb,wt)*1e3; wtd.quantile(Inf16fb,wt,c(1,39)/40)*1e3
# IMP TB lf ls
wtd.mean(ImpTb16-ImpAct16,wt)*1e3; wtd.quantile(ImpTb16-ImpAct16,wt,c(1,39)/40)*1e3

### MDR-TB
Inc16Mdr   <- rowSums(gc(2016,c(159:160,164:165,169:170,174:175)))/rowSums(gc(2016,156:175))
Inc50Mdr   <- rowSums(gc(2050,c(159:160,164:165,169:170,174:175)))/rowSums(gc(2050,156:175))
Inc16MdrT  <- rowSums(gc(2016,c(159:160,164:165,169:170,174:175)))
Inc50MdrT  <- rowSums(gc(2050,c(159:160,164:165,169:170,174:175)))

## MDR TB per mil
wtd.mean(Inc16MdrT/Pop16,wt)*1e6; wtd.quantile(Inc16MdrT/Pop16,wt,c(1,39)/40)*1e6
wtd.mean(Inc50MdrT/Pop50,wt)*1e6; wtd.quantile(Inc50MdrT/Pop50,wt,c(1,39)/40)*1e6

## MDR TB %
wtd.mean(Inc16Mdr,wt)*100; wtd.quantile(Inc16Mdr,wt,c(1,39)/40)*100
wtd.mean(Inc50Mdr,wt)*100; wtd.quantile(Inc50Mdr,wt,c(1,39)/40)*100

## LY lost
TBdeath16   <- rowSums(gc(2016,c(89:99,100:110)))
TBdeath50   <- rowSums(gc(2050,c(89:99,100:110)))
LYlost16   <- (gc(2016,89:99)+gc(2016,100:110))%*%LE
LYlost50   <- (gc(2050,89:99)+gc(2050,100:110))%*%LE
Prev16     <- rowSums(gc(2016,18:19))
Prev50     <- rowSums(gc(2050,18:19))
Prev16tx   <- rowSums(gc(2016,c(18:19,21:24)))
Prev50tx   <- rowSums(gc(2050,c(18:19,21:24)))

## Burden LY lost
wtd.mean(TBdeath16,wt)*1e6; wtd.quantile(TBdeath16,wt,c(1,39)/40)*1e6
wtd.mean(LYlost16,wt)*1e3; wtd.quantile(LYlost16,wt,c(1,39)/40)*1e3
wtd.mean(Prev16tx,wt)*1e3; wtd.quantile(Prev16tx,wt,c(1,39)/40)*1e3

wtd.mean(TBdeath50,wt)*1e6; wtd.quantile(TBdeath50,wt,c(1,39)/40)*1e6
wtd.mean(LYlost50,wt)*1e3; wtd.quantile(LYlost50,wt,c(1,39)/40)*1e3
wtd.mean(Prev50tx,wt)*1e3; wtd.quantile(Prev50tx,wt,c(1,39)/40)*1e3


######################################

ImpTb16 = ImpAct16 = ImpAct16Mdr = NULL
load("ModelInputs_9-2-16.rData")
ImmigInputs      <- Inputs[["ImmigInputs"]]
library(MASS)
##########  SQ EXP COVARIANCE FUNCTION
# http://www.jameskeirstead.ca/blog/gaussian-process-regression-with-r/
calcSigma <- function(X1,X2,l=25,sd=0.1) {  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
for (i in 1:nrow(Sigma)) { for (j in 1:ncol(Sigma)) { Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)  }  }
return(Sigma*sd^2)  }
## l=10,sd=0.1
x <- c(0,0.1)
k.xx <- calcSigma(x,x);     k.xxs <- calcSigma(x,0:84)
k.xsx <- calcSigma(0:84,x); k.xsxs <- calcSigma(0:84,0:84)
vcv_gp_l10_sd0.1 <- k.xsxs - k.xsx%*%solve(k.xx)%*%k.xxs
SmoCurve <- function(vec) { jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2),x=seq(1,length(vec),1/12))$y; jj[jj<0] <- 0 ; jj }

for(ee in 1:nrow(parAll)) { # ee=1

  P <- parAll[ee,-(1:2)]
  # ######  IMMIGRATION  ######
  # overall
  # TotImmig1       <- c(ImmigInputs[[1]][1:65],(ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*P["ImmigVolFut"]+ImmigInputs[[1]][66])/12*P["ImmigVol"]
  TotImmig0       <- (c(ImmigInputs[[1]][1:151])+c(rep(0,65),cumsum(rep(P["ImmigVolFut"],86))))/12*P["ImmigVol"]
  set.seed( P["rand_seed"]+0)
  TotImmig1       <-   c(TotImmig0[1:66],exp(mvrnorm(1, log(TotImmig0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  #  TotImmig1       <-   TotImmig0
  TotImmig        <- SmoCurve(TotImmig1)
  TotImmAge       <- outer(TotImmig,ImmigInputs[["AgeDist"]])
  # latent 1
  PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:65]^P["TunLtbiTrend"]*ImmigInputs[["PrevTrend25_34"]][65]^(1-P["TunLtbiTrend"]),
                        ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutLat"]/0.99)^(1:86))
  set.seed( P["rand_seed"]+1)
  PrevTrend25_341l <-   c(PrevTrend25_340l[1:66],exp(mvrnorm(1, log(PrevTrend25_340l[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  #  PrevTrend25_341l <-   PrevTrend25_340l
  PrevTrend25_34l  <- SmoCurve(PrevTrend25_341l)
  PrevTrend25_34_ls <- (PrevTrend25_34l); PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
  ImmLat          <- matrix(NA,length(PrevTrend25_34_ls),11)
  for(i in 1:11) ImmLat[,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]*P["LtbiPar1"]-(c(2.5,1:9*10,100)/100)[i]^2*P["LtbiPar2"])*PrevTrend25_34_ls))*TotImmAge[,i]
  
  # active 
  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:65],ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutAct"]/0.99)^(1:86))
  set.seed( P["rand_seed"]+2)
  PrevTrend25_341a <-   c(PrevTrend25_340a[1:66],exp(mvrnorm(1, log(PrevTrend25_340a[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  # PrevTrend25_341a <-   PrevTrend25_340a
  PrevTrend25_34a  <- SmoCurve(PrevTrend25_341a)
  ImDxChngV      <- SmoCurve(c(rep(1,57),seq(1,P["ImDxChng"],length.out=6)[-1],rep(P["ImDxChng"],89)))
  ImmAct         <- outer(PrevTrend25_34a*P["RRtbprev"]*ImDxChngV,ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*P["pImAct"]
  ImmFst         <- outer(PrevTrend25_34a*P["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-P["pImAct"])
  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat
  
  ###
  # tx history, dr
  TxExpAge       <- ImmigInputs[["Fract_TxE_by_Age"]]
  DrN0            <- ImmigInputs[["DR_TB_by_year_N"]]
  DrE0            <- ImmigInputs[["DR_TB_by_year_E"]]
  # DR variability trend for Immigrants 
  DrTrend0       <- c(rep(P["TunImmDr"],65),P["TunImmDr"]*P["TunImmDrFut"]^(1:86))
  set.seed( P["rand_seed"]+3)
  DrTrend1       <- c(DrTrend0[1:66],exp(mvrnorm(1, log(DrTrend0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  #  DrTrend1       <- DrTrend0
  DrTrend2       <- SmoCurve(DrTrend1)
  
  OddsDrN        <- (1-DrN0[,1])/DrN0[,1] * DrTrend2;
  pDrN           <- OddsDrN/(1+OddsDrN)
  DrN            <- DrN0
  DrN[,-1]       <- DrN0[,-1]*pDrN/rowSums(DrN0[,-1])
  DrN[is.nan(DrN)] <- 0
  DrN[,1]        <- 1-rowSums(DrN[,-1])
  
  OddsDrE        <- (1-DrE0[,1])/DrE0[,1] * DrTrend2
  pDrE           <- OddsDrE/(1+OddsDrE)
  DrE            <- DrE0
  DrE[,-1]       <- DrE0[,-1]*pDrE/rowSums(DrE0[,-1])
  DrE[is.nan(DrE)] <- 0
  DrE[,1]        <- 1-rowSums(DrE[,-1])
  
  ImpTb16[ee]     <- sum(ImmFst[(2016-1950)*12+7,]+ImmAct[(2016-1950)*12+7,]+ImmLat[(2016-1950)*12+7,])*12
  ImpAct16[ee]    <- sum(ImmAct[(2016-1950)*12+7,])*12
  ImpAct16Mdr[ee] <- sum(ImmAct[(2016-1950)*12+7,]*(sum(DrN[(2016-1950)*12+7,4:5])*(1-TxExpAge)+sum(DrE[(2016-1950)*12+7,4:5])*TxExpAge))*12
  if(ee/100==round(ee/100)) { cat('\r',ee); flush.console() } }


###############

 save(list=c("ImpTb16","ImpAct16","ImpAct16Mdr"),file="AddnOutcomes_V2.rData")

# load("AddnOutcomes_V2.rData") # ImpTb16","ImpAct16","ImpAct16Mdr

###########################################################
################## SCENARIOS
###########################################################
# SCENARIOS
for(aa in 0:8) { # aa=8
  ResTabfb <- ResTabus <- ResTab <- array(NA,dim=c(nrow(parAll),90,6))
    for(bb in 0:14) { load(paste("cluster results/MiAllex",bb,"_",aa,"_9-21-16.rData",sep=""))
    nr <- if(bb==14) { nrow(parAll)-14000 } else { 1000 }
    o <- get(paste("MiAllex",bb,sep=""))
    ResTab[1:nr+bb*1000,,1] <- o[,-(1:60),1]
    ResTab[1:nr+bb*1000,,2] <- apply(o[,-(1:60),304:307],c(1,2),sum)/o[,-(1:60),2]*1e6
    ResTab[1:nr+bb*1000,,3] <- apply(o[,-(1:60),c(56:66,67:77)],c(1,2),sum)/o[,-(1:60),2]*1e2
    ResTab[1:nr+bb*1000,,4] <-  (o[,-(1:60),136]+o[,-(1:60),215])/o[,-(1:60),2]*1e6
    ResTab[1:nr+bb*1000,,5] <- apply(o[,-(1:60),c(159:160,164:165,169:170,174:175)],c(1,2),sum)/apply(o[,-(1:60),156:175],c(1,2),sum)*1e2
    ResTab[1:nr+bb*1000,,6] <- apply(o[,-(1:60),89:110],c(1,2),sum)/o[,-(1:60),2]*1e6
    
    ResTabus[1:nr+bb*1000,,1] <- o[,-(1:60),1]
    ResTabus[1:nr+bb*1000,,2] <- apply(o[,-(1:60),304:305],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
    ResTabus[1:nr+bb*1000,,3] <- apply(o[,-(1:60),56:66],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e2
    ResTabus[1:nr+bb*1000,,4] <-apply(o[,-(1:60),c(235:245,246:256)],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
    ResTabus[1:nr+bb*1000,,5] <- apply(o[,-(1:60),c(159:160,164:165)],c(1,2),sum)/apply(o[,-(1:60),156:165],c(1,2),sum)*1e2
    ResTabus[1:nr+bb*1000,,6] <- o[,-(1:60),308]/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
    
    ResTabfb[1:nr+bb*1000,,1] <- o[,-(1:60),1]
    ResTabfb[1:nr+bb*1000,,2] <- apply(o[,-(1:60),306:307],c(1,2),sum)/(o[,-(1:60),32]+o[,-(1:60),33])*1e6
    ResTabfb[1:nr+bb*1000,,3] <- apply(o[,-(1:60),67:77],c(1,2),sum)/(o[,-(1:60),32]+o[,-(1:60),33])*1e2
    ResTabfb[1:nr+bb*1000,,4] <- (apply(o[,-(1:60),c(137:147,216:226)],c(1,2),sum)-apply(o[,-(1:60),c(235:245,246:256)],c(1,2),sum))/
      (o[,-(1:60),32]+o[,-(1:60),33])*1e6
    ResTabfb[1:nr+bb*1000,,5] <- apply(o[,-(1:60),c(169:170,174:175)],c(1,2),sum)/apply(o[,-(1:60),166:175],c(1,2),sum)*1e2
    ResTabfb[1:nr+bb*1000,,6] <- o[,-(1:60),309]/(o[,-(1:60),32]+o[,-(1:60),33])*1e6
    
    print(paste(aa,"--",bb)); flush.console()  }
  assign(paste("ResTab",aa,sep=""),ResTab)
  assign(paste("ResTabus",aa,sep=""),ResTabus)
  assign(paste("ResTabfb",aa,sep=""),ResTabfb)  }

  ResTabAllfb <- ResTabAllus <- ResTabAll <- list()
  for(aa in 0:8)  ResTabAll[[aa+1]] <- get(paste("ResTab",aa,sep=""))
  for(aa in 0:8)  ResTabAllfb[[aa+1]] <- get(paste("ResTabfb",aa,sep=""))
  for(aa in 0:8)  ResTabAllus[[aa+1]] <- get(paste("ResTabus",aa,sep=""))
   save(ResTabAll  ,file="ScenResTab_9-21-16.rData")
   save(ResTabAllus,file="ScenResTabus_9-21-16.rData")
   save(ResTabAllfb,file="ScenResTabfb_9-21-16.rData")

  # load("ScenResTab_9-21-16.rData") #  ResTabAll
  # load("ScenResTabus_9-21-16.rData") #  ResTabAllus
  # load("ScenResTabfb_9-21-16.rData") #  ResTabAllfb

  MeanResTab <- MeanResTabUS <- MeanResTabFB <- list()
  for(aa in 0:8)  MeanResTab[[aa+1]] <- apply(ResTabAll[[aa+1]],2:3,function(x) weighted.mean(x,wt))
  for(aa in 0:8)  MeanResTabUS[[aa+1]] <- apply(ResTabAllus[[aa+1]],2:3,function(x) weighted.mean(x,wt))
  for(aa in 0:8)  MeanResTabFB[[aa+1]] <- apply(ResTabAllfb[[aa+1]],2:3,function(x) weighted.mean(x,wt))
  
    save(MeanResTab,file="ScenResTab2_9-21-16.rData")
    save(MeanResTabUS,file="ScenResTab2US_9-21-16.rData")
    save(MeanResTabFB,file="ScenResTab2FB_9-21-16.rData")

# Outcomes shown in graph at 2025, 2050, and 2100, + percent reduction compared to base case
ResTab <- array(NA,dim=c(5*2,9*3,3))
dimnames(ResTab)[[1]] <- c("Incident M. tb Infection (per Mil)","Percent of Base Case Value",
                            "LTBI Prevalence (%)","Percent of Base Case Value",
                            "New TB Cases (per Mil)","Percent of Base Case Value",
                            "MDR-TB in New TB Cases (%)","Percent of Base Case Value",
                            "TB-Related Deaths (per Mil)","Percent of Base Case Value")          
dimnames(ResTab)[[2]] <- c(paste("Int",rep(0:5,each=3),rep(c("_mu","_lo","_hi"),6),sep=""),paste("Scen",rep(1:3,each=3),rep(c("_mu","_lo","_hi"),3),sep=""))
dimnames(ResTab)[[3]] <- c(2025,2050,2100)
ResTabfb <- ResTabus <- ResTab

# func for getting summary stats
l <- function(x,w,q){if(q==1){wtd.mean(x,w)}else{if(q==2){wtd.quantile(x,w,1/40)}else{wtd.quantile(x,w,39/40)}}}
l(ResTabAll[[1]][,90,2],wt,3)

for(y0 in 1:3) {
  y = c(2025,2050,2099)[y0]-2009 
  for(oc in 1:5) {
    for(sc in 1:9) {
      ResTab[oc*2-1,sc*3-2,y0] <- l(ResTabAll[[sc]][,y,oc+1],wt,1)
      ResTab[oc*2-1,sc*3-1,y0] <- l(ResTabAll[[sc]][,y,oc+1],wt,2)
      ResTab[oc*2-1,sc*3-0,y0] <- l(ResTabAll[[sc]][,y,oc+1],wt,3)
      ResTab[oc*2,sc*3-2,y0] <- l(ResTabAll[[sc]][,y,oc+1]/ResTabAll[[1]][,y,oc+1],wt,1)*100
      ResTab[oc*2,sc*3-1,y0] <- l(ResTabAll[[sc]][,y,oc+1]/ResTabAll[[1]][,y,oc+1],wt,2)*100
      ResTab[oc*2,sc*3-0,y0] <- l(ResTabAll[[sc]][,y,oc+1]/ResTabAll[[1]][,y,oc+1],wt,3)*100
      
      ResTabus[oc*2-1,sc*3-2,y0] <- l(ResTabAllus[[sc]][,y,oc+1],wt,1)
      ResTabus[oc*2-1,sc*3-1,y0] <- l(ResTabAllus[[sc]][,y,oc+1],wt,2)
      ResTabus[oc*2-1,sc*3-0,y0] <- l(ResTabAllus[[sc]][,y,oc+1],wt,3)
      ResTabus[oc*2,sc*3-2,y0] <- l(ResTabAllus[[sc]][,y,oc+1]/ResTabAllus[[1]][,y,oc+1],wt,1)*100
      ResTabus[oc*2,sc*3-1,y0] <- l(ResTabAllus[[sc]][,y,oc+1]/ResTabAllus[[1]][,y,oc+1],wt,2)*100
      ResTabus[oc*2,sc*3-0,y0] <- l(ResTabAllus[[sc]][,y,oc+1]/ResTabAllus[[1]][,y,oc+1],wt,3)*100
      
      ResTabfb[oc*2-1,sc*3-2,y0] <- l(ResTabAllfb[[sc]][,y,oc+1],wt,1)
      ResTabfb[oc*2-1,sc*3-1,y0] <- l(ResTabAllfb[[sc]][,y,oc+1],wt,2)
      ResTabfb[oc*2-1,sc*3-0,y0] <- l(ResTabAllfb[[sc]][,y,oc+1],wt,3)
      ResTabfb[oc*2,sc*3-2,y0] <- l(ResTabAllfb[[sc]][,y,oc+1]/ResTabAllfb[[1]][,y,oc+1],wt,1)*100
      ResTabfb[oc*2,sc*3-1,y0] <- l(ResTabAllfb[[sc]][,y,oc+1]/ResTabAllfb[[1]][,y,oc+1],wt,2)*100
      ResTabfb[oc*2,sc*3-0,y0] <- l(ResTabAllfb[[sc]][,y,oc+1]/ResTabAllfb[[1]][,y,oc+1],wt,3)*100
    } } }

# Outcomes shown in graph at 2025, 2050, and 2100, + percent reduction compared to base case
ResTabb <- matrix(NA,5*2,9)
colnames(ResTabb) <- c("Base case",paste("Int",1:5,sep=""),paste("Scen",1:3,sep=""))
rownames(ResTabb) <- c("New infections per million","Percent of base case value",
                       "LTBI prevalence (%)","Percent of base case value",
                       "TB cases per million","Percent of base case value",
                       "MDR-TB in new TB cases (%)","Percent of base case value",
                       "TB deaths per million","Percent of base case value")    
ResTabb <- list(ResTab25=ResTabb,ResTab50=ResTabb,ResTab00=ResTabb)
ResTabbus <- ResTabbfb <- ResTabb

for(j in 1:3) {
  for(rw in 1:10) {
  ns <- c(1,0,2,0,1,0,2,0,2,0)[rw]
ResTabb[[j]][rw,] <-  paste(gsub(" ","",format(round(ResTab[rw,1+0:8*3,j],ns),nsmall=ns),fixed=T)," (",
                            gsub(" ","",format(round(ResTab[rw,2+0:8*3,j],ns),nsmall=ns),fixed=T),", ",
                            gsub(" ","",format(round(ResTab[rw,3+0:8*3,j],ns),nsmall=ns),fixed=T),")",sep="") 
ResTabbus[[j]][rw,] <-  paste(gsub(" ","",format(round(ResTabus[rw,1+0:8*3,j],ns),nsmall=ns),fixed=T)," (",
                            gsub(" ","",format(round(ResTabus[rw,2+0:8*3,j],ns),nsmall=ns),fixed=T),", ",
                            gsub(" ","",format(round(ResTabus[rw,3+0:8*3,j],ns),nsmall=ns),fixed=T),")",sep="") 
ResTabbfb[[j]][rw,] <-  paste(gsub(" ","",format(round(ResTabfb[rw,1+0:8*3,j],ns),nsmall=ns),fixed=T)," (",
                            gsub(" ","",format(round(ResTabfb[rw,2+0:8*3,j],ns),nsmall=ns),fixed=T),", ",
                            gsub(" ","",format(round(ResTabfb[rw,3+0:8*3,j],ns),nsmall=ns),fixed=T),")",sep="") 
  }
  ResTabb[[j]][1:5*2,1] <- ResTabbus[[j]][1:5*2,1] <-  ResTabbfb[[j]][1:5*2,1] <- "100%"
  }

write.csv(ResTabb[[1]],file="ResTab_2025_9-21-16.csv")
write.csv(ResTabb[[2]],file="ResTab_2050_9-21-16.csv")
write.csv(ResTabb[[3]],file="ResTab_2100_9-21-16.csv")

write.csv(ResTabbus[[1]],file="ResTabUS_2025_9-21-16.csv")
write.csv(ResTabbus[[2]],file="ResTabUS_2050_9-21-16.csv")
write.csv(ResTabbus[[3]],file="ResTabUS_2100_9-21-16.csv")

write.csv(ResTabbfb[[1]],file="ResTabFB_2025_9-21-16.csv")
write.csv(ResTabbfb[[2]],file="ResTabFB_2050_9-21-16.csv")
write.csv(ResTabbfb[[3]],file="ResTabFB_2100_9-21-16.csv")

#####################################
####### Results for base case at various points, by nativity
ResTabBC <- array(NA,dim=c(nrow(parAll),90,6,3))
# BASE CASE
for(bb in 0:14) { load(paste("cluster results/MiAllex",bb,"_0_9-21-16.rData",sep="")); print(bb); flush.console() }

dimnames(ResTabBC)[[2]] <- 2010:2099
dimnames(ResTabBC)[[3]] <- c("Year","New infections per million","LTBI prevalence (%)","TB cases per million",
                           "MDR-TB in new TB cases (%)","TB deaths per million")   
dimnames(ResTabBC)[[4]] <- c("All","US-born","Foreign-born")

for(bb in 0:14) { # bb=1
  nr <- if(bb==14) { nrow(parAll)-14000 } else { 1000 }
  o <- get(paste("MiAllex",bb,sep=""))
  ResTabBC[1:nr+bb*1000,,1,1] <- o[,-(1:60),1]
  ResTabBC[1:nr+bb*1000,,2,1] <- apply(o[,-(1:60),304:307],c(1,2),sum)/o[,-(1:60),2]*1e6
  ResTabBC[1:nr+bb*1000,,3,1] <- apply(o[,-(1:60),c(56:66,67:77)],c(1,2),sum)/o[,-(1:60),2]*1e2
  ResTabBC[1:nr+bb*1000,,4,1] <-  (o[,-(1:60),136]+o[,-(1:60),215])/o[,-(1:60),2]*1e6
  ResTabBC[1:nr+bb*1000,,5,1] <- apply(o[,-(1:60),c(159:160,164:165,169:170,174:175)],c(1,2),sum)/apply(o[,-(1:60),156:175],c(1,2),sum)*1e2
  ResTabBC[1:nr+bb*1000,,6,1] <- apply(o[,-(1:60),308:309],c(1,2),sum)/o[,-(1:60),2]*1e6

  ResTabBC[1:nr+bb*1000,,1,2] <- o[,-(1:60),1]
  ResTabBC[1:nr+bb*1000,,2,2] <- apply(o[,-(1:60),304:305],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
  ResTabBC[1:nr+bb*1000,,3,2] <- apply(o[,-(1:60),56:66],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e2
  ResTabBC[1:nr+bb*1000,,4,2] <-apply(o[,-(1:60),c(235:245,246:256)],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
  ResTabBC[1:nr+bb*1000,,5,2] <- apply(o[,-(1:60),c(159:160,164:165)],c(1,2),sum)/apply(o[,-(1:60),156:165],c(1,2),sum)*1e2
  ResTabBC[1:nr+bb*1000,,6,2] <- o[,-(1:60),308]/(o[,-(1:60),30]+o[,-(1:60),31])*1e6

  ResTabBC[1:nr+bb*1000,,1,3] <- o[,-(1:60),1]
  ResTabBC[1:nr+bb*1000,,2,3] <- apply(o[,-(1:60),306:307],c(1,2),sum)/(o[,-(1:60),32]+o[,-(1:60),33])*1e6
  ResTabBC[1:nr+bb*1000,,3,3] <- apply(o[,-(1:60),67:77],c(1,2),sum)/(o[,-(1:60),32]+o[,-(1:60),33])*1e2
  ResTabBC[1:nr+bb*1000,,4,3] <- (apply(o[,-(1:60),c(137:147,216:226)],c(1,2),sum)-apply(o[,-(1:60),c(235:245,246:256)],c(1,2),sum))/
                                  (o[,-(1:60),32]+o[,-(1:60),33])*1e6
  ResTabBC[1:nr+bb*1000,,5,3] <- apply(o[,-(1:60),c(169:170,174:175)],c(1,2),sum)/apply(o[,-(1:60),166:175],c(1,2),sum)*1e2
  ResTabBC[1:nr+bb*1000,,6,3] <- o[,-(1:60),309]/(o[,-(1:60),32]+o[,-(1:60),33])*1e6
  print(bb); flush.console()  }

 save("ResTabBC",file="ScenResTabBC_9-21-16.rData")
# load("ScenResTabBC_9-21-16.rData")  # ResTab0
#  SUMMARY STATS
ResTabBC2 <- array(NA,dim=c(5*2,4*3,3))
dimnames(ResTabBC2)[[1]] <- c("Incident M. tb Infection (per Mil)","Percent of Base Case Value",
                           "LTBI Prevalence (%)","Percent of Base Case Value",
                           "New TB Cases (per Mil)","Percent of Base Case Value",
                           "MDR-TB in New TB Cases (%)","Percent of Base Case Value",
                           "TB-Related Deaths (per Mil)","Percent of Base Case Value")          
dimnames(ResTabBC2)[[2]] <- paste("Year",rep(c(2016,2025,2050,2100),each=3),rep(c("_mu","_lo","_hi"),4),sep="")
dimnames(ResTabBC2)[[3]] <- c("All","US","FB")
# func for getting summary stats
l <- function(x,w,q){if(q==1){wtd.mean(x,w)}else{if(q==2){wtd.quantile(x,w,1/40)}else{wtd.quantile(x,w,39/40)}}}

for(pp in 1:3) {
  for(oc in 1:5) {
    for(y0 in 1:4) {
      y = c(2016,2025,2050,2099)[y0]-2009 
      ResTabBC2[oc*2-1,y0*3-2,pp] <- l(ResTabBC[,y,oc+1,pp],wt,1)
      ResTabBC2[oc*2-1,y0*3-1,pp] <- l(ResTabBC[,y,oc+1,pp],wt,2)
      ResTabBC2[oc*2-1,y0*3-0,pp] <- l(ResTabBC[,y,oc+1,pp],wt,3)
      
      ResTabBC2[oc*2,y0*3-2,pp] <- l(ResTabBC[,y,oc+1,pp]/ResTabBC[,2016-2009,oc+1,pp],wt,1)*100
      ResTabBC2[oc*2,y0*3-1,pp] <- l(ResTabBC[,y,oc+1,pp]/ResTabBC[,2016-2009,oc+1,pp],wt,2)*100
      ResTabBC2[oc*2,y0*3-0,pp] <- l(ResTabBC[,y,oc+1,pp]/ResTabBC[,2016-2009,oc+1,pp],wt,3)*100
    } } }

# Outcomes shown in graph at 2025, 2050, and 2100, + percent reduction compared to base case
ResTabb <- matrix(NA,5*2,4)
colnames(ResTabb) <- c(2016,2025,2050,2100)
rownames(ResTabb) <- c("New infections per million","Percent of 2016 value",
                       "LTBI prevalence (%)","Percent of 2016 value",
                       "TB cases per million","Percent of 2016 value",
                       "MDR-TB in new TB cases (%)","Percent of 2016 value",
                       "TB deaths per million","Percent of 2016 value")    

ResTabb <- list(ResTabAll=ResTabb,ResTabUS=ResTabb,ResTabFB=ResTabb)

for(j in 1:3) {
  for(rw in 1:10) {
    ns <- c(1,0,2,0,1,0,2,0,2,0)[rw]
    ResTabb[[j]][rw,] <-  paste(gsub(" ","",format(round(ResTabBC2[rw,1+0:3*3,j],ns),nsmall=ns),fixed=T)," (",
                                gsub(" ","",format(round(ResTabBC2[rw,2+0:3*3,j],ns),nsmall=ns),fixed=T),", ",
                                gsub(" ","",format(round(ResTabBC2[rw,3+0:3*3,j],ns),nsmall=ns),fixed=T),")",sep="") 
  }
  
  ResTabb[[j]][1:5*2,1] <- "100%"   }

write.csv(ResTabb[[1]],file="ResTabBC_All_9-21-2016.csv")
write.csv(ResTabb[[2]],file="ResTabBC_US_9-21-2016.csv")
write.csv(ResTabb[[3]],file="ResTabBC_FB_9-21-2016.csv")



###################################
## 'End TB' reach elim target?

dimnames(ResTabAll[[1]])

1-l(ResTabAll[[1]][,41,4]/ResTabAll[[1]][,7,4],wt,q=1)
1-l(ResTabAll[[1]][,41,4]/ResTabAll[[1]][,7,4],wt,q=2)
1-l(ResTabAll[[1]][,41,4]/ResTabAll[[1]][,7,4],wt,q=3)

1-l(ResTabAll[[6]][,41,4]/ResTabAll[[1]][,7,4],wt,q=1)
1-l(ResTabAll[[6]][,41,4]/ResTabAll[[1]][,7,4],wt,q=2)
1-l(ResTabAll[[6]][,41,4]/ResTabAll[[1]][,7,4],wt,q=3)
