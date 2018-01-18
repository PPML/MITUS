##########     TB TRANSMISSION MODEL 2015    ##########   
##########            PARAM FILE             ##########

if(Int5==1) {  Int1 = Int2 = Int3 = Int4 = 1   }

library(MASS)

## Some functions
  ORAdd <- function(val,OR) { x <- (val/(1-val))*OR; x/(1+x) }
  lgt <-  function(x) log(x/(1-x));  invlgt <- function(x) 1/(1+exp(-x))

##########  Logistic Curve function 
  LgtCurve <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
    zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2100-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr)/12)
    zz  <- as.numeric(EndVal)/(1+exp(-zz));  if(StYr>1950) { zz[1:((StYr-1950)*12)] <- 0 };    zz  }

  SmoCurve <- function(vec) { jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2),x=seq(1,length(vec),1/12))$y; jj[jj<0] <- 0 ; jj }

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

######### INPUTS ######### 
  BgMort           <- Inputs[["BgMort"]]
  InitPop          <- Inputs[["InitPop"]]
  Births           <- Inputs[["Births"]]
  ImmigInputs      <- Inputs[["ImmigInputs"]]
  HIV_incid_inputs <- Inputs[["HIV_incid_inputs"]]
  TxInputs         <- Inputs[["TxInputs"]]

######### PARAMETER DEFINITIONS ######### 

## BIRTHS 
  Birthst   <- SmoCurve(Births)*P["TunBirths"]/12		# indexed by t, abs number of new adult entrants over time 

######  MORTALITY RATES  ######  
## background 
  mubt      <- matrix(NA,1801,11)
  for(i in 1:11) { mubt[,i] <- SmoCurve(BgMort[,i+1])*P["TunMubt"]/12  } 		# indexed by t, background mortality rate over time 
## disease specific
  muIp  	  <- P["muIp"]/12
  RRmuIn	  <- P["RRmuIn"]
  muTbH     <- P["muTbH"]/12  
  TunmuTbAg <- P["TunmuTbAg"] # multiplier of mort rate above   
  TunmuHvAg <- P["TunmuTbAg"] # ffs.
  
  muH1      <- P["muH1"]/12
  muH2      <- P["muH2"]/12
  muT2      <- P["muT2"]/12
  RRmuHR    <- c(1,P["RRmuHR"],1,1)
  vHMort    <- matrix(0,11,5); 
  rownames(vHMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vHMort) <- c("N0","H1","T1","H2","T2")
  vHMort[,c(1,3)] <- 0; 
  vHMort[,2] <- muH1*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  vHMort[,4] <- muH2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  vHMort[,5] <- muT2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  
  vTMort   <- matrix(0,11,11); 
  rownames(vTMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vTMort) <- c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp")
  vTMort[,c(5,8,10)] <- muIp*RRmuIn;   vTMort[,c(6,9,11)] <- muIp
  RRmuTbAg <- exp(c(0,0,1:9)*TunmuTbAg)
  for(i in 1:ncol(vTMort)) {  vTMort[,i] <- vTMort[,i] * RRmuTbAg  }

# ######  IMMIGRATION  ######
# overall
 # TotImmig1       <- c(ImmigInputs[[1]][1:65],(ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*P["ImmigVolFut"]+ImmigInputs[[1]][66])/12*P["ImmigVol"]
  TotImmig0       <- (c(ImmigInputs[[1]][1:151])+c(rep(0,65),cumsum(rep(P["ImmigVolFut"],86))))/12*P["ImmigVol"]
  set.seed( P["rand_seed"]+0)
#  TotImmig1       <-   c(TotImmig0[1:66],exp(mvrnorm(1, log(TotImmig0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  TotImmig1       <-   TotImmig0
  TotImmig        <- SmoCurve(TotImmig1)
  TotImmAge       <- outer(TotImmig,ImmigInputs[["AgeDist"]])
# latent 1
  PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:65]^P["TunLtbiTrend"]*ImmigInputs[["PrevTrend25_34"]][65]^(1-P["TunLtbiTrend"]),
                        ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutLat"]/0.99)^(1:86))
  set.seed( P["rand_seed"]+1)
 # PrevTrend25_341l <-   c(PrevTrend25_340l[1:66],exp(mvrnorm(1, log(PrevTrend25_340l[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
    PrevTrend25_341l <-   PrevTrend25_340l
  PrevTrend25_34l  <- SmoCurve(PrevTrend25_341l)
  PrevTrend25_34_ls <- (PrevTrend25_34l); PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
  ImmLat          <- matrix(NA,length(PrevTrend25_34_ls),11)
  for(i in 1:11) ImmLat[,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]*P["LtbiPar1"]-(c(2.5,1:9*10,100)/100)[i]^2*P["LtbiPar2"])*PrevTrend25_34_ls))*TotImmAge[,i]
 
  # active 
  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:65],ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutAct"]/0.99)^(1:86))
  set.seed( P["rand_seed"]+2)
 # PrevTrend25_341a <-   c(PrevTrend25_340a[1:66],exp(mvrnorm(1, log(PrevTrend25_340a[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  PrevTrend25_341a <-   PrevTrend25_340a
  PrevTrend25_34a  <- SmoCurve(PrevTrend25_341a)
  ImDxChngV      <- SmoCurve(c(rep(1,57),seq(1,P["ImDxChng"],length.out=6)[-1],rep(P["ImDxChng"],89)))
  ImmAct         <- outer(PrevTrend25_34a*P["RRtbprev"]*ImDxChngV,ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*P["pImAct"]
  ImmFst         <- outer(PrevTrend25_34a*P["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-P["pImAct"])
  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat
  
  #### #### #### INT 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  pctDoc <- (1-0.28)
  if(Int1==1) {
    for(i in 1:11) ImmLat[,i] <- ImmLat[,i]*(1-LgtCurve(2016,2021,1)*P["SensLt"]*P["EffLt"]*(1-P["pDefLt"])*pctDoc)
    for(i in 1:11) ImmAct[,i] <- ImmAct[,i]*(1-LgtCurve(2016,2021,1)*P["SensLt"]*P["EffLt"]*(1-P["pDefLt"])*pctDoc)
    for(i in 1:11) ImmFst[,i] <- ImmFst[,i]*(1-LgtCurve(2016,2021,1)*P["SensLt"]*P["EffLt"]*(1-P["pDefLt"])*pctDoc)
    ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat   }
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  #### #### #### SCEN 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Scen2==1) {
    for(i in 1:11) ImmLat[,i] <- ImmLat[,i]*(1-LgtCurve(2016,2017,1))
    for(i in 1:11) ImmAct[,i] <- ImmAct[,i]*(1-LgtCurve(2016,2017,1))  
    for(i in 1:11) ImmFst[,i] <- ImmFst[,i]*(1-LgtCurve(2016,2017,1))  
    ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat   }
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  #### #### #### SCEN 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Scen3==1) {
    adjf <- (1-0.5)^(1/10/12);  adjfV <- adjf^(0:1008)
    for(i in 1:11) ImmLat[793:1801,i] <- ImmLat[793:1801,i]*adjfV
    for(i in 1:11) ImmAct[793:1801,i] <- ImmAct[793:1801,i]*adjfV  
    for(i in 1:11) ImmFst[793:1801,i] <- ImmFst[793:1801,i]*adjfV
    ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat   }
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
# tx history, dr
  TxExpAge       <- ImmigInputs[["Fract_TxE_by_Age"]]
  DrN0            <- ImmigInputs[["DR_TB_by_year_N"]]
  DrE0            <- ImmigInputs[["DR_TB_by_year_E"]]
# DR variability trend for Immigrants 
  DrTrend0       <- c(rep(P["TunImmDr"],65),P["TunImmDr"]*P["TunImmDrFut"]^(1:86))
  set.seed( P["rand_seed"]+3)
 # DrTrend1       <- c(DrTrend0[1:66],exp(mvrnorm(1, log(DrTrend0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  DrTrend1       <- DrTrend0
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

# p(smr-pos|active TB)
  p_Imm_SP       <- 0.2
# Exogenous infection risk  
  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  for(i in 1:5) ExogInf[,i] <- P["ExogInf"]*PrevTrend25_34a/PrevTrend25_341a["2013"]*(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])/12

######  EMIGRATION  ###### 
  rEmmigFB <- c(P["rEmmigF1"],P["rEmmigF2"])/12

######  HIGH-RISK ENTRY/EXIT  ######
  p_HR     <- P["pHR"]
  yr       <- c(2.5,1:10*10)
  r0_5     <- 1/3; r45_55 <- 1/20
  HR_exit  <- r0_5*((r45_55/r0_5)^(1/(50-2.5)))^(yr-2.5)
  HR_entry <- HR_exit*p_HR*1.3
  HrEntEx  <- cbind(HR_entry,HR_exit)/12

######  HIV NATURAL HISTORY  ######

## HIV incidence
  rHIVt0 <- c(HIV_incid_inputs[[1]][1:61],HIV_incid_inputs[[1]][61]*P["HivIncFut"]^(1:90))*P["HivInc"]/12
  set.seed( P["rand_seed"]+4)
 # rHIVt1 <-  c(rHIVt0[1:66],exp(mvrnorm(1, log(rHIVt0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  rHIVt1 <-  rHIVt0
  rHIVt2   <- SmoCurve(rHIVt1)
  rHIVt2[1:(20*12)] <- 0
  rHIVt    <- outer(rHIVt2,HIV_incid_inputs[[2]])
  HivHrPar <- P["HivHrPar"]

## HIV progression
  mHxtoHy      <- c(P["rH1toH2"],P["rT1toT2"])/12
  TunHIVProgAg <- P["TunHIVProgAg"]
  vHxtoHy   <- outer(exp(c(0,0,1:7,7,7)*TunHIVProgAg),mHxtoHy)

## Tuning TB natural history for HIV/ART status
  ArtTbEff1   <- P["ArtTbEff1"] 
  ArtTbEff2   <- P["ArtTbEff2"] 
  TbHivEarly  <- P["TbHivEarly"]  # par values for HIV>350 = Advanced HIV vals * TbHivEarly + HIV neg vals * (1-TbHivEarly)

## STRAIN FITNESS (with pansensitive=1.0) N.B. applied to contact rates
  RelFit     <- c(1,ORAdd(c(0.95,0.85,0.73,0.73),P["TunRelFit"]))  
                
######  TB TRANSMISSION  ######  
  CR           <- P["CR"]/12
  TrIn         <- P["TrIn"]	# Contact rate for In as a fraction of Ip
  RelInfRg     <- c(1.0,P["RelCrHr"],1.0)*CR
  RelInfHivt   <- LgtCurve(1990,1996,1) + (1-LgtCurve(1990,1996,1))*P["RelCrHr"]
  TunTbTransTx <- P["TunTbTransTx"]  # set to zero?
  Vmix         <- 1-c(P["sigmaHiv"],P["sigmaHr"],P["sigmaFb"])
  RelInf       <- matrix(0,11,5); rownames(RelInf) <- c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp"); colnames(RelInf) <- 1:5
  RelInf[c(5,8,10),] <- TrIn; RelInf[c(6,9,11),] <- 1
  RelInf[8:11,] <- RelInf[8:11,]*TunTbTransTx
  for(i in 1:5) RelInf[,i] <- RelInf[,i]*RelFit[i]

######  TB NATURAL HISTORY  ######  
### EARLY EPIDEMIC
  Early0 <- P["Early0"]
  EarlyTrend <- c(rep(1+Early0,200*12),seq(1+Early0,1.0,length.out=50*12+2))

## PROGRESSION TO DISEASE
  pfast     <- P["pfast"]
  pimmed    <- P["pimmed"]
  ORpfast1  <- P["ORpfast1"]
  ORpfast2  <- P["ORpfast2"]
  ORpfastH  <- P["ORpfastH"]
  ORpfastPI <- P["ORpfastPI"]
  rslow     <- P["rslow"]/12
  rslowH    <- P["rslowH"]/12
  rfast     <- P["rfast"]/12
  rrSlowFB0 <- P["rrSlowFB"]
  rrSlowFB  <- c(1,1,rrSlowFB0,rrSlowFB0)

  Mpfast       <- matrix(NA,11,5)
  Mpfast[,]    <- pfast/(1-pfast)
  Mpfast[1,]   <- Mpfast[1,]*ORpfast1
  Mpfast[2,]   <- Mpfast[2,]*ORpfast2
  Mpfast[,4]   <- Mpfast[,4]*ORpfastH
  
  MpfastPI     <- Mpfast
  MpfastPI[,1] <- Mpfast[,1]*ORpfastPI
  
  Mpfast[,]    <- Mpfast[,]  /(1+Mpfast[,]);
  MpfastPI[,]  <- MpfastPI[,]/(1+MpfastPI[,])
  Mpfast[,2]   <- Mpfast[,4]  *TbHivEarly   + Mpfast[,1]  *(1-TbHivEarly); 
  MpfastPI[,2] <- MpfastPI[,4]*TbHivEarly   + MpfastPI[,1]*(1-TbHivEarly);
  Mpfast[,3]   <- Mpfast[,2]  *(1-ArtTbEff1)+ Mpfast[,1]  *ArtTbEff1;     
  MpfastPI[,3] <- MpfastPI[,2]*(1-ArtTbEff1)+ MpfastPI[,1]*ArtTbEff1;
  Mpfast[,5]   <- Mpfast[,4]  *(1-ArtTbEff2)+ Mpfast[,1]  *ArtTbEff2;  
  MpfastPI[,5] <- MpfastPI[,4]*(1-ArtTbEff2)+ MpfastPI[,1]*ArtTbEff2;

  Vrslow     <- rep(rslow,5)
  Vrslow[4]  <- rslowH
  Vrslow[2]  <- Vrslow[4]*TbHivEarly   +Vrslow[1]*(1-TbHivEarly)
  Vrslow[3]  <- Vrslow[2]*(1-ArtTbEff1)+Vrslow[1]*ArtTbEff1
  Vrslow[5]  <- Vrslow[4]*(1-ArtTbEff2)+Vrslow[1]*ArtTbEff2

  TunrslowAge  <- P["TunrslowAge"]
  rrReactAg       <- exp(c(0,0,0,0,0,0,0.5,1:4)*P["TunrslowAge"]) 
  Mrslow <- outer(rrReactAg,Vrslow)
  
## RATE OF SELF CURE, by smear status and HIV status
  rRecov     <-  P["rRecov"]/12
  
## PROBABILITY OF GOING INTO Ip FOLLOWING BREAKDOWN, by HIV status
  pSmPos    <- P["pSmPos"] 
  ORpSmPos1 <- P["ORpSmPos1"]
  ORpSmPosH <- P["ORpSmPosH"]
  
  MpSmPos      <- matrix(NA,11,5)
  MpSmPos[,]  <- pSmPos/(1-pSmPos)
  MpSmPos[1,] <- MpSmPos[1,]*ORpSmPos1
  MpSmPos[2,] <- MpSmPos[2,]
  MpSmPos[,4] <- MpSmPos[,4]*ORpSmPosH
  MpSmPos[,]  <- MpSmPos[,]/(1+MpSmPos[,])
  
  MpSmPos[,2] <- MpSmPos[,4]*TbHivEarly+MpSmPos[,1]*(1-TbHivEarly)
  MpSmPos[,3] <- MpSmPos[,2]*(1-ArtTbEff1)+MpSmPos[,1]*ArtTbEff1
  MpSmPos[,5] <- MpSmPos[,4]*(1-ArtTbEff2)+MpSmPos[,1]*ArtTbEff2
  
## RATE OF SELF CURE, by smear status and HIV status
  rSlfCur      <- P["rSlfCur"]/12
  VrSlfCur     <- rep(rSlfCur,5)
  VrSlfCur[4]  <- 0
  VrSlfCur[2]  <- VrSlfCur[4]*TbHivEarly   +VrSlfCur[1]*(1-TbHivEarly)
  VrSlfCur[3]  <- VrSlfCur[2]*(1-ArtTbEff1)+VrSlfCur[1]*ArtTbEff1
  VrSlfCur[5]  <- VrSlfCur[4]*(1-ArtTbEff2)+VrSlfCur[1]*ArtTbEff2

## RATE OF CONVERSION FROM In TO Ip
  rSmConv  	   <- P["rSmConv"]/12
  VrSmConv     <- rep(rSmConv,5)
  VrSmConv[4]  <- 0
  VrSmConv[2]  <- VrSmConv[4]*TbHivEarly   +VrSmConv[1]*(1-TbHivEarly)
  VrSmConv[3]  <- VrSmConv[2]*(1-ArtTbEff1)+VrSmConv[1]*ArtTbEff1
  VrSmConv[5]  <- VrSmConv[4]*(1-ArtTbEff2)+VrSmConv[1]*ArtTbEff2

### LTBI DIAGNOSIS  
  rLtScrt       <- LgtCurve(1985,2015,P["rLtScr"])/12
  SensLt        <- P["SensLt"]    #  sens of test for latent TB infection (based on IGRA QFT-GIT)
  SensLtHiv     <- SensLt*P["rrSensLtHiv"]    #  sens of test for latent TB infection (based on IGRA QFT-GIT) with HIV infection
  SpecLt        <- P["SpecLt"]    #  spec of test for latent TB infection (based on IGRA QFT-GIT)
  SpecLtFb      <- SpecLt         #  spec of test for latent TB infection (based on IGRA QFT-GIT) in foreign-born (assumed BCG exposed)
  rrTestHiv     <- P["rrTestHiv"] # RR of LTBI screening for HIV and HR as cmpared to general
  rrTestHr      <- P["rrTestHr"] # RR of LTBI screening for HIV and HR as cmpared to general
  rrTestLrNoTb  <- P["rrTestLrNoTb"] # RR of LTBI screening for individuals with no risk factors
  dLt           <- 1/9
  rDefLt        <- dLt*P["pDefLt"]/(1-P["pDefLt"])  # based on 50% tx completion with 6 mo INH regimen 2.0 [1.0,3.0] from Menzies Ind J Med Res 2011  
  EffLt         <- P["EffLt"]     
  LtTxPar       <- c(dLt,rDefLt,EffLt)
  
  dLtt          <- (1-LgtCurve(2016,2021,0)) * 1/9
  
  EffLt0        <- LgtCurve(2016,2021,0)+1
  EffLtX        <- cbind(EffLt0,0,EffLt0,0,0)
  
  #### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Int2==1) { rLtScrt    <- rLtScrt + LgtCurve(2016,2021,1)*rLtScrt*1
                EffLtX[,2] <- LgtCurve(2016,2021,1)
                dLtt       <- (1-LgtCurve(2016,2021,1))*(1/9) + LgtCurve(2016,2021,1)*(1/3)  }
  
  #### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  LtDxPar <- matrix(NA,6,2);  colnames(LtDxPar) <- c("latent","no latent");  rownames(LtDxPar) <- c("LR","HR","FB","LRh","HRh","FBh")
  LtDxPar[,1] <- c(SensLt               , rrTestHr*SensLt   , SensLt,
                   rrTestHiv*SensLtHiv, max(rrTestHiv,rrTestHr)*SensLtHiv, rrTestHiv*SensLtHiv)
  LtDxPar[,2] <- c(rrTestLrNoTb*(1-SpecLt), rrTestHr*(1-SpecLt), (1-SpecLtFb), 
                   rrTestHiv*(1-SpecLt) , max(rrTestHiv,rrTestHr)*(1-SpecLt), rrTestHiv*(1-SpecLtFb))

  pImmScen    <- P["pImmScen"] # lack of reactivitiy to IGRA for Sp

######  TB DIAGNOSIS  ######  

## TEST CHARACTERISTICS
  SensSn    <- P["SensSn"]
  SensSp    <- P["SensSp"]

## PROVIDER DELAY
  DelaySn    <- P["DelaySn"]
  DelaySp    <- P["DelaySp"]

## RR Delay Elderly
  TunRRdxAge    <- P["TunRRdxAge"]
  RRdxAge       <- 1+(c(rep(1,6),cumprod(seq(1.05,1.3,length.out=5)))-1)*TunRRdxAge

## Attendance rate
  bspline <- function(x,k,i,m) {
      if (m==-1) { res <- as.numeric(x<k[i+1]&x>=k[i])  } else { 
      z0  <- (x-k[i]) / (k[i+m+1]-k[i]);  z1  <- (k[i+m+2]-x) / (k[i+m+2]-k[i+1])
      z0*bspline(x,k,i,m-1) + z1*bspline(x,k,i+1,m-1) }  }

  n_Spln   <- 5; n_Stps   <- 2010-1950+1; dif_pen   <- 1 # quadratic spline 
# Working...
  x1    <- seq(1,n_Stps);  k1  <- seq(min(x1),max(x1),length=n_Spln-dif_pen)   
  dk1   <- k1[2]-k1[1] ;   k1  <- c(k1[1]-dk1*((dif_pen+1):1),k1,k1[n_Spln-dif_pen]+dk1*(1:(dif_pen+1))) 
  SpMat <- matrix(NA,nrow=n_Stps,ncol=n_Spln); for(i in 1:n_Spln){ SpMat[,i] <- bspline(x=x1,k=k1,m=dif_pen,i=i)  }

  DxPri          <- (seq(0.5,4,length.out=5)-0.75)*3.5/2.626+0.25  # Prior from 0.5 t to 4.0
  rDx            <- c(SpMat%*%(DxPri+c(P["Dx1"],P["Dx2"],P["Dx3"],P["Dx4"],P["Dx5"])),62:151)
  rDx[62:151]    <- rDx[61] + (rDx[61]-rDx[60])*cumsum((0.75^(1:90)))
  rDxt0          <- SmoCurve(rDx)/12;  
  rDxt1          <- cbind(rDxt0,rDxt0)  

# Put it all together
  rDxt           <- cbind(1/(1/rDxt1+DelaySn)*SensSn,1/(1/rDxt1+DelaySp)*SensSp)
  rDxt[,2]       <- (rDxt[,1]-min(rDxt[,1]))/P["rrDxH"]+min(rDxt[,1])
  rDxt[,4]       <- (rDxt[,2]-min(rDxt[,2]))/P["rrDxH"]+min(rDxt[,2])
  colnames(rDxt) <- c("Sn","Sn_H","Sp","SP_H")
  
  #### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Int3==1) { for(i in 1:4) { rDxt[,i] <- rDxt[,i]+ rDxt[,i]*LgtCurve(2016,2021,1)     }   }
  
  #### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
## DST
  pDstt       <- LgtCurve(1985,2000,P["pDst"]) 

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Int4==1) {  pDstt <- 1-(1-pDstt)*(1-LgtCurve(2016,2021,0.5))      }
  
  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
######  TREATMENT OUTCOMES  ######  
  TunTxMort	<- P["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

# Regimen duration (no uncertainty?)
  d1st <- 1/9;  dInh <- 1/9;  dRif <- 1/21;  dMdr <- 1/21; dXdr <- 1/21

## Regimen efficacy
  pCurPs  <- P["pCurPs"]    # probability of cure with pansensitive TB, 1st line regimen (Menzies 09)
  TxE1		<- P["TxEf1"]     # RR cure given mono-resistance, 1st line regimen 0.90 
  TxE2		<- P["TxEf2"]     # RR cure given multiresistance, 1st line or 2nd line 0.50 
  TxE3		<- P["TxEf3"]     # RR cure given effective 2nd line 0.90 

## Default
  rDef0         <- rep(NA,151)
  rDef0[1:30]   <- P["TxDefEarly"]
  rDef0[44:63]  <- ORAdd(TxInputs[[1]][,2],P["TunTxDef"])
  rDef0[64:151] <- rDef0[63] 
  rDef1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=rDef0[-(31:43)],spar=0.4),x=1950:2100)$y
  rDeft         <- SmoCurve(rDef1)/12;
  rDeftH        <- rDeft*P["RRdefHR"] # second col is HR default rate
  
  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Int4==1) {  rDeftH <- rDeftH*(1-LgtCurve(2016,2021,0.5))      }
  if(Int4==1) {  rDeft  <- rDeft* (1-LgtCurve(2016,2021,0.5))      }
  
  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  ## Tx quality
  TxQual0         <- rep(NA,151)
  TxQual0[1:30]   <- P["TxQualEarly"]
  TxQual0[44:62]  <- ORAdd(TxInputs[[2]][,2],P["TunTxQual"])
  TxQual0[63:151] <- TxQual0[62] 
  TxQual1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=TxQual0[-(31:43)],spar=0.4),x=1950:2100)$y
  TxQualt         <- SmoCurve(TxQual1);
  RRcurDef        <- P["RRcurDef"]
  
  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  if(Int4==1) {  TxQualt <- 1-(1-TxQualt)*(1-LgtCurve(2016,2021,0.5))      }
  
  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  #### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  NixTrans <- 1-LgtCurve(2016,2017,1)
  if(Scen1==0) {  NixTrans[] <- 1      }
  
  #### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  ## Retreatment  
  pReTx   <- LgtCurve(1985,2000,P["pReTx"])   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)

  ## Aquired resistance
  pAR1  	<- P["pAR1"]
  pAR2		<- P["pAR2"]
  pAR3		<- P["pAR3"]
  pAR4		<- P["pAR4"]
  pAR5  	<- P["pAR5"]
  fAR  	  <- P["fAR"]

## TB treatment table 
  TxMat           <- matrix(NA,7,10)
  rownames(TxMat) <- c("TxCompRate","TxEff","PanSens","MonoINH","MonoRIF","MDRTB","XDRTB")
  colnames(TxMat) <- paste(rep(c("FirstLine","Mdr"),5),"_Str",rep(1:5,each=2),sep="")
# COLS                 1LS1 2LS1 1LS2 2LS2 1LS3 2LS3 1LS4 2LS4 1LS5 2LS5 
  TxMat[1,]       <- c(d1st,d1st,d1st,dInh,d1st,dRif,d1st,dMdr,d1st,dMdr)
  TxMat[2,]       <- c(1   ,1   ,TxE1,TxE3,TxE1,TxE3,TxE2,TxE3,TxE2,TxE2)*pCurPs
  TxMat[3,]       <- c(0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   )
  TxMat[4,]       <- c(pAR1,pAR1,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   )
  TxMat[5,]       <- c(pAR2,pAR2,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   )
  TxMat[6,]       <- c(pAR3,pAR3,pAR4,pAR2,pAR4,pAR1,0   ,0   ,0   ,0   )
  TxMat[7,]       <- c(0   ,0   ,0   ,0   ,0   ,0   ,0   ,pAR5,0   ,0   )

######  ART TREATMENT  ###### 
## ART initiation
  rArtInitH <- LgtCurve(1994,1997,1)*P["rArtInit15"]/12
  rArtInitL <- rArtInitH*LgtCurve(2005,2015,1)*P["rrArtInitL"]
  rArtInit  <- cbind(rArtInitL,rArtInitH,rArtInitL*P["rrArtInitHR"],rArtInitH*P["rrArtInitHR"])

## ART default rate
  rArtDef    <- P["r_ArtLtfu"]/12

StatList <- noquote(list(c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p"),
                         c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp"),
                         1:5,c("NT","TL","TX"),c("N0","H1","T1","H2","T2"),c("LR","HR","F1","F2")))

ResNam <- c("Year",                                         # year
            "N_ALL",                                        # total pop  
            paste("N",StatList[[1]],sep="_"),               # pop by ag cat
            paste("N",StatList[[2]],sep="_"),               # pop by tb cat
            paste("N",StatList[[5]],sep="_"),               # pop by hv cat
            paste("N",StatList[[6]],sep="_"),               # pop by rg cat
            paste("N_US",StatList[[1]],sep="_"),            # US pop by ag cat
            paste("N_FB",StatList[[1]],sep="_"),            # FB pop by ag cat
            paste("N_US_LTBI",StatList[[1]],sep="_"),       # US LTBI pop by ag cat 
            paste("N_FB_LTBI",StatList[[1]],sep="_"),       # FB LTBI pop by ag cat
            paste("N_HIV",StatList[[1]],sep="_"),           # HIV pop by ag cat
            paste("TBMORT_HIVNEG",StatList[[1]],sep="_"),   # TB mort, HIV neg, by ag cat
            paste("TBMORT_HIVPOS",StatList[[1]],sep="_"),   # TB mort, HIV pos, by ag cat
            paste("HIVMORT",StatList[[1]],sep="_"),         # HIV mort, by ag cat
            paste("TOTMORT",StatList[[1]],sep="_"),         # total mort, by ag cat
            "TBTX_COMPLT","TBTX_DISCONT","TBTX_DIED",       # TB treatment outcomes complete, discontinue, death
            "NOTIF_ALL",                                    # total notif
            paste("NOTIF",StatList[[1]],sep="_"),           # notif by ag cat
            paste("NOTIF",c("N","E"),sep="_"),              # notif by tx cat
            paste("NOTIF_HIV",c("POS","NEG"),sep="_"),      # notif by HIV pos/neg
            paste("NOTIF",StatList[[6]],sep="_"),           # notif by rg cat
            paste("NOTIF_US_N",StatList[[3]],sep="_"),      # notif, US, N, by dr cat
            paste("NOTIF_US_E",StatList[[3]],sep="_"),      # notif, US, E, by dr cat
            paste("NOTIF_FB_N",StatList[[3]],sep="_"),      # notif, FB, N, by dr cat
            paste("NOTIF_FB_E",StatList[[3]],sep="_"),      # notif, FB, E, by dr cat
            "TLTBI_INITS",                                  # Initiations on LTBI tx
            "TLTBI_INITS_FB",                               # Initiations on LTBI tx FB
            "TLTBI_INITS_HR",                               # Initiations on LTBI tx HR
            "TLTBI_INITS_HV",                               # Initiations on LTBI tx HV
            "TLTBI_INITS_TP",                               # Initiations on LTBI tx, with LTBI
            "INCID_ALL",                                    # Total incidence
            paste("INCID_ALL",StatList[[1]],sep="_"),       # Total incidence by ag cat
            "INCID_ALL_US",                                 # Total incidence, US born
            "INCID_ALL_FB",                                 # Total incidence, foreign born
            "INCID_ALL_FB2",                                # Total incidence, foreign born
            "INCID_ALL_HR",                                 # Total incidence, high risk
            "INCID_ALL_HV",                                 # Total incidence, HIV pos
            "INCID_REC",                                    # Total incidence, recent infection
            paste("INCID_REC",StatList[[1]],sep="_"),       # Total incidence by ag cat, recent infection
            "INCID_REC_US",                                 # Total incidence, US born, recent infection
            "INCID_REC_FB",                                 # Total incidence, foreign born, recent infection
            "INCID_REC_FB2",                                # Total incidence, foreign born, recent infection
            "INCID_REC_HR",                                 # Total incidence, high risk, recent infection
            "INCID_REC_HV",                                 # Total incidence, HIV pos, recent infection
            "NOTIF_MORT_ALL",                               # total notif, dead at diagnosis
            paste("NOTIF_MORT",StatList[[1]],sep="_"),      # notif by ag cat, dead at diagnosis
            paste("NOTIF_MORT",c("N","E"),sep="_"),         # notif by tx cat, dead at diagnosis
            paste("NOTIF_MORT_HIV",c("POS","NEG"),sep="_"), # notif by HIV pos/neg, dead at diagnosis
            paste("NOTIF_MORT",StatList[[6]],sep="_"),      # notif by rg cat, dead at diagnosis
            paste("NOTIF_US",StatList[[1]],sep="_"),        # notif by ag cat, US only
            paste("NOTIF_US_MORT",StatList[[1]],sep="_"),   # notif by ag cat, dead at diagnosis US only
            paste("NOTIF_MORT_HIV_Neg",StatList[[1]],sep="_"),    # notif by ag cat, dead at diagnosis HIV neg
            paste("TOTMORT_W_HIV",StatList[[1]],sep="_"),   # total mort, by ag cat, have HIV
            paste("TOTMORT_W_TB",StatList[[1]],sep="_"),     # total mort, by ag cat, have active TB
            c("N_Ls_US","N_Lf_US","N_In_US","N_Ip_US"),
            c("N_Ls_FB","N_Lf_FB","N_In_FB","N_Ip_FB"),
            c("ARTI_LR","ARTI_HR","ARTI_FB","ARTI_LR_H","ARTI_HR_H","ARTI_FB_H"), # Force of infection
            c("TB_INF_LR","TB_INF_HR","TB_INF_F1","TB_INF_F2"), # NEW TB INFECTIONS
            paste("TBMORT",c("US","FB"),sep="_") )
length(ResNam)

###################################################

