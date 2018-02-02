################################################################################
##### THE CODE BELOW WILL GATHER AND FORMAT INPUTS FOR THE TB_MODEL        #####
##### FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME; #####
##### VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.        #####
################################################################################
##### THIS FILE USES FUNCTIONS FOUND IN HELPER_FUNCTIONS.R                 #####
################################################################################
library(MASS)
################################################################################
###########################          INPUTS            #########################
################################################################################
  BgMort           <- Inputs[["BgMort"]]
  InitPop          <- Inputs[["InitPop"]]
  Births           <- Inputs[["Births"]]
  ImmigInputs      <- Inputs[["ImmigInputs"]]
####### WILL NEED TO ADD RF FACTOR INCIDENCE RATES IN PLACE OF HIV INCIDENCE
  TxInputs         <- Inputs[["TxInputs"]]

##########                PARAMETER DEFINITIONS                      ###########
#######################           BIRTHS                 #######################
####### INDEXED BY TIME, ABSOLUTE NUMBER OF NEW ADULT ENTRANTS OVER TIME #######

<<<<<<< HEAD
Birthst   <- SmoCurve(Births)*P["TunBirths"]/12
=======
  Birthst   <- SmoCurve(Births)*P["TunBirths"]/12
>>>>>>> 6a3d0294ae8f4e35c62d8b0a3fb6f50350f3d04f

##########################      MORTALITY RATES       ##########################
########################## BACKGROUND MORTALITY BY TIME ########################
  mubt      <- matrix(NA,1801,11)
  for(i in 1:11) {
  	mubt[,i] <- SmoCurve(BgMort[,i+1])*P["TunMubt"]/12
  }
<<<<<<< HEAD
## disease specific
  muIp  	  <- P["muIp"]/12
  RRmuIn	  <- P["RRmuIn"]
  muTbH     <- P["muTbH"]/12
  TunmuTbAg <- P["TunmuTbAg"] # multiplier of mort rate above
  TunmuHvAg <- P["TunmuTbAg"]

  muH1      <- P["muH1"]/12
  muH2      <- P["muH2"]/12
  muT2      <- P["muT2"]/12
  RRmuHR    <- c(1,P["RRmuHR"],1,1)
=======
#########################     DISEASE SPECIFIC       ###########################
  muIp  	  <- P["muIp"]/12
  RRmuIn	  <- P["RRmuIn"]
  muTbH     <- P["muTbH"]/12
  TunmuTbAg <- P["TunmuTbAg"] # multiplier of mort rate above
######################## MULTIPLER OF MORT RATE ABOVE ########################
  TunmuHvAg <- P["TunmuTbAg"] ##is this a bug??
############ CONVERT ANNUAL RATES OF HIV MORTALITY TO MONTHLY RATES ##########
  muH1      <- P["muH1"]/12
  muH2      <- P["muH2"]/12
#### IS THIS ART?
  muT2      <- P["muT2"]/12
###############  RATE RATIO OF MORTALITY INCREASE FOR HIGH RISK ###############
h RRmuHR    <- c(1,P["RRmuHR"],1,1)
############### CREATE A MATRIX OF HIV MORTALITIES BY AGE GROUP ###############
#ALL THIS NEEDS TO BE UPDATED
>>>>>>> 6a3d0294ae8f4e35c62d8b0a3fb6f50350f3d04f
  vHMort    <- matrix(0,11,5);
  rownames(vHMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vHMort) <- c("N0","H1","T1","H2","T2")
  vHMort[,c(1,3)] <- 0;
  vHMort[,2] <- muH1*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  vHMort[,4] <- muH2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  vHMort[,5] <- muT2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
############### CREATE A MATRIX OF TB MORTALITIES BY AGE GROUP ###############
  vTMort   <- matrix(0,11,11);
  rownames(vTMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vTMort) <- c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp")
  vTMort[,c(5,8,10)] <- muIp*RRmuIn;   vTMort[,c(6,9,11)] <- muIp
  RRmuTbAg <- exp(c(0,0,1:9)*TunmuTbAg)
  for(i in 1:ncol(vTMort)) {  vTMort[,i] <- vTMort[,i] * RRmuTbAg  }

# ######  IMMIGRATION  ######
# overall

######################         IMMIGRATION             ########################
######################         OVERALL IMM.            ########################
 # TotImmig1       <- c(ImmigInputs[[1]][1:65],(ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*P["ImmigVolFut"]+ImmigInputs[[1]][66])/12*P["ImmigVol"]
  TotImmig0       <- (c(ImmigInputs[[1]][1:151])+c(rep(0,65),cumsum(rep(P["ImmigVolFut"],86))))/12*P["ImmigVol"]
#  TotImmig1       <-   c(TotImmig0[1:66],exp(mvrnorm(1, log(TotImmig0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  TotImmig1       <-   TotImmig0
  TotImmig        <- SmoCurve(TotImmig1)
  TotImmAge       <- outer(TotImmig,ImmigInputs[["AgeDist"]])
  # latent 1

######################           LTBI IMM.             ########################
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
######################         ACTIVE TB IMM.           ########################
  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:65],ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutAct"]/0.99)^(1:86))
#  PrevTrend25_341a <-   c(PrevTrend25_340a[1:66],exp(mvrnorm(1, log(PrevTrend25_340a[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  PrevTrend25_341a <-   PrevTrend25_340a
  PrevTrend25_34a  <- SmoCurve(PrevTrend25_341a)
  ImDxChngV      <- SmoCurve(c(rep(1,57),seq(1,P["ImDxChng"],length.out=6)[-1],rep(P["ImDxChng"],89)))
  ImmAct         <- outer(PrevTrend25_34a*P["RRtbprev"]*ImDxChngV,ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*P["pImAct"]
  ImmFst         <- outer(PrevTrend25_34a*P["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-P["pImAct"])
  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat
# tx history, dr
  TxExpAge       <- ImmigInputs[["Fract_TxE_by_Age"]]
 ####REMOVED ALL THE DRUG RESISTANCE PARAMETERS #########

# p(smr-pos|active TB)
#  p_Imm_SP       <- 0.2
# Exogenous infection risk
  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  for(i in 1:5) ExogInf[,i] <- P["ExogInf"]*PrevTrend25_34a/PrevTrend25_341a["2013"]*(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])/12

######  EMIGRATION  ######
  rEmmigFB <- c(P["rEmmigF1"],P["rEmmigF2"])/12

######  HIGH-RISK ENTRY/EXIT  ######

######################        p(smr-pos|active TB)        ######################
 ## p_Imm_SP       <- 0.2
######################         EXOGENOUS INF. RISK        ######################
  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  for(i in 1:5)
      ExogInf[,i] <- P["ExogInf"]*PrevTrend25_34a/PrevTrend25_341a["2013"]*(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])/12

######################         EMMIGRATION             ########################
  rEmmigFB <- c(P["rEmmigF1"],P["rEmmigF2"])/12

######################       HIGH-RISK ENTRY/EXIT      ########################
  p_HR     <- P["pHR"]
  yr       <- c(2.5,1:10*10)
  r0_5     <- 1/3; r45_55 <- 1/20
  HR_exit  <- r0_5*((r45_55/r0_5)^(1/(50-2.5)))^(yr-2.5)
  HR_entry <- HR_exit*p_HR*1.3
  HrEntEx  <- cbind(HR_entry,HR_exit)/12

## STRAIN FITNESS (with pansensitive=1.0) N.B. applied to contact rates
  RelFit     <- c(1,ORAdd(c(0.95,0.85,0.73,0.73),P["TunRelFit"]))

######  TB TRANSMISSION  ######
######################       TB TRANSMISSION           #######################
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
#  pimmed    <- P["pimmed"]
  ORpfast1  <- P["ORpfast1"]
  ORpfast2  <- P["ORpfast2"]
  ORpfastH  <- P["ORpfastH"]
  ORpfastPI <- P["ORpfastPI"]
  rslow     <- P["rslow"]/12
#  rslowH    <- P["rslowH"]/12
  rfast     <- P["rfast"]/12
  rrSlowFB0 <- P["rrSlowFB"]
  rrSlowFB  <- c(1,1,rrSlowFB0,rrSlowFB0)

  Mpfast       <- matrix(NA,11,5)
######################      TB NATURAL HISTORY       ##########################
######################        EARLY EPIDEMIC         ##########################
  Early0 <- P["Early0"]
  EarlyTrend <- c(rep(1+Early0,200*12),seq(1+Early0,1.0,length.out=50*12+2))
######################     PROGRESSION TO DISEASE     ##########################
  pfast      <- P["pfast"]
  pimmed     <- P["pimmed"]
  ORpfast1   <- P["ORpfast1"]
  ORpfast2   <- P["ORpfast2"]
  ORpfastRF  <- P["ORpfastRF"]
  ORpfastPI  <- P["ORpfastPI"]
  rslow      <- P["rslow"]/12
  rslowRF    <- P["rslowRF"]/12
  rfast      <- P["rfast"]/12
  rrSlowFB0  <- P["rrSlowFB"]
  rrSlowFB   <- c(1,1,rrSlowFB0,rrSlowFB0)

  Mpfast       <- matrix(NA,11,5)
############## CREATE AN ODDS FROM THE PROB OF FAST PROGRESSION ##############
  Mpfast[,]    <- pfast/(1-pfast)
  Mpfast[1,]   <- Mpfast[1,]*ORpfast1
  Mpfast[2,]   <- Mpfast[2,]*ORpfast2
  Mpfast[,4]   <- Mpfast[,4]*ORpfastH

  MpfastPI     <- Mpfast
  MpfastPI[,1] <- Mpfast[,1]*ORpfastPI


#################       CREATE A NEW MATRIX PARTIAL. IMM.     #################
  MpfastPI     <- Mpfast
  MpfastPI[,1] <- Mpfast[,1]*ORpfastPI
#######NEED TO UPDATE THIS SECTION
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
#######################       RATE OF SELF CURE         ########################
#by smear status and HIV status
  rRecov     <-  P["rRecov"]/12

######## PROBABILITY OF GOING INTO Ip FOLLOWING BREAKDOWN, by HIV status   #####
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

#######################       RATE OF SELF CURE         ########################
# by smear status and HIV status
  rSlfCur      <- P["rSlfCur"]/12
  VrSlfCur     <- rep(rSlfCur,5)
  VrSlfCur[4]  <- 0
  VrSlfCur[2]  <- VrSlfCur[4]*TbHivEarly   +VrSlfCur[1]*(1-TbHivEarly)
  VrSlfCur[3]  <- VrSlfCur[2]*(1-ArtTbEff1)+VrSlfCur[1]*ArtTbEff1
  VrSlfCur[5]  <- VrSlfCur[4]*(1-ArtTbEff2)+VrSlfCur[1]*ArtTbEff2

## RATE OF CONVERSION FROM In TO Ip
######################     RATE OF CONVERSION FROM      ########################
##############       TB SMEAR NEGATIVE TO TB SMEAR POS        ##################
  rSmConv  	   <- P["rSmConv"]/12
  VrSmConv     <- rep(rSmConv,5)
  VrSmConv[4]  <- 0
  VrSmConv[2]  <- VrSmConv[4]*TbHivEarly   +VrSmConv[1]*(1-TbHivEarly)
  VrSmConv[3]  <- VrSmConv[2]*(1-ArtTbEff1)+VrSmConv[1]*ArtTbEff1
  VrSmConv[5]  <- VrSmConv[4]*(1-ArtTbEff2)+VrSmConv[1]*ArtTbEff2

### LTBI DIAGNOSIS
######################          LTBI DIAGNOSIS           ########################
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

################################################################################
#######################         TB DIAGNOSIS            ########################
#######################      TEST CHARACTERISTICS       ########################
  SensSn    <- P["SensSn"]
  SensSp    <- P["SensSp"]

#######################        PROVIDER DELAY           ########################

  DelaySn    <- P["DelaySn"]
  DelaySp    <- P["DelaySp"]

#######################     PROVIDER DELAY RR ELDERLY     ########################


  TunRRdxAge    <- P["TunRRdxAge"]
  RRdxAge       <- 1+(c(rep(1,6),cumprod(seq(1.05,1.3,length.out=5)))-1)*TunRRdxAge

#######################         ATTENDANCE RATE           ########################

#######################     B SPLINE HELPER FUNCTION      ########################
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
## DST
  pDstt       <- LgtCurve(1985,2000,P["pDst"])

######  TREATMENT OUTCOMES  ######
  TunTxMort	<- P["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

# Regimen duration (no uncertainty?)
  d1st <- 1/9;  dInh <- 1/9;  dRif <- 1/21;  dMdr <- 1/21; dXdr <- 1/21

## Regimen efficacy

################################################################################
###########################     TREATMENT OUTCOMES    ##########################
################################################################################

  TunTxMort	<- P["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

###########################      REGIMEN DURATION      ##########################
# (no uncertainty?)

  d1st <- 1/9;
  dInh <- 1/9;
  dRif <- 1/21;
  dMdr <- 1/21;
  dXdr <- 1/21

###########################       REGIMEN EFFICACY      ##########################

  pCurPs  <- P["pCurPs"]    # probability of cure with pansensitive TB, 1st line regimen (Menzies 09)
  TxE1		<- P["TxEf1"]     # RR cure given mono-resistance, 1st line regimen 0.90
  TxE2		<- P["TxEf2"]     # RR cure given multiresistance, 1st line or 2nd line 0.50
  TxE3		<- P["TxEf3"]     # RR cure given effective 2nd line 0.90

## Default

###########################        REGIMEN DEFAULT       ##########################

  rDef0         <- rep(NA,151)
  rDef0[1:30]   <- P["TxDefEarly"]
  rDef0[44:63]  <- ORAdd(TxInputs[[1]][,2],P["TunTxDef"])
  rDef0[64:151] <- rDef0[63]
  rDef1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=rDef0[-(31:43)],spar=0.4),x=1950:2100)$y
  rDeft         <- SmoCurve(rDef1)/12;
  rDeftH        <- rDeft*P["RRdefHR"] # second col is HR default rate
## Tx quality


#########################        REGIMEN QUALITY       ##########################

  TxQual0         <- rep(NA,151)
  TxQual0[1:30]   <- P["TxQualEarly"]
  TxQual0[44:62]  <- ORAdd(TxInputs[[2]][,2],P["TunTxQual"])
  TxQual0[63:151] <- TxQual0[62]
  TxQual1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=TxQual0[-(31:43)],spar=0.4),x=1950:2100)$y
  TxQualt         <- SmoCurve(TxQual1);
  RRcurDef        <- P["RRcurDef"]

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
#########################         RETREATMENT         ##########################
  pReTx   <- LgtCurve(1985,2000,P["pReTx"])   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)

################################################################################
#######    CREATED COMPLETION RATE THAT IS AN AVERAGE OF ALL TREATMENTS   ######
#######    RATE SHOULD BE UPDATED TO A WEIGHTED AVERAGE WHEN DATA AVAIL.  ######
################################################################################
  TxCompRate  <-(d1st+dInh+dRif+dMdr+dXdr)/5

################################################################################
#######      CREATED TX QUALITY THAT IS AN AVERAGE OF ALL TREATMENTS      ######
#######    QUAL SHOULD BE UPDATED TO A WEIGHTED AVERAGE WHEN DATA AVAIL.  ######
################################################################################

  TxEff       <- [(1+TxE1+TxE2+TxE3)/4]*pCurPs

################################################################################
##########################      ART TREATMENT     ##############################
################################################################################
##########################     ART INITIATION     ##############################

  rArtInitH <- LgtCurve(1994,1997,1)*P["rArtInit15"]/12
  rArtInitL <- rArtInitH*LgtCurve(2005,2015,1)*P["rrArtInitL"]
  rArtInit  <- cbind(rArtInitL,rArtInitH,rArtInitL*P["rrArtInitHR"],rArtInitH*P["rrArtInitHR"])

## ART default rate
  rArtDef    <- P["r_ArtLtfu"]/12

StatList <- noquote(list(c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p"),
                         c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp"),
                         1:5,c("NT","TL","TX"),c("N0","H1","T1","H2","T2"),c("LR","HR","F1","F2")))

##########################    ART DEFAULT RATE    ##############################

  rArtDef    <- P["r_ArtLtfu"]/12

################################################################################
##### CREATE A LIST TO HOLD THE VECTORS FOR AGE CATEGORIES, TB STATES,     #####
##### DRUG RESISTANCE, TREATMENT HISTORY, HIV STATUS, AND RISK CATEGORY.   #####
################################################################################
################################################################################
StatList <- noquote(list(
#####					                  AGE CATEGORIES                             #####
#####           AGE GROUPS 0-4, 5-14, 15-24, ... , 95-94, 95+              #####
  c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p"),
#####					     				    TUBERCULOSIS STATES		              			   #####
##### SUSCEPTIBLE; UNINFECTED & PARTIALLY IMMUNE; LATENT SLOW; LATENT FAST;#####
#####       ACTIVE TB SMEAR NEG.; ACTIVE TB SMEAR POS.; TB TREATMENT       #####
  c("Su","Sp","Ls","Lf","An", "Ap", "Tx"),
#####                          TREATMENT HISTORY                           #####
#####              NO TB TREATMENT HISTORY, PRIOR LTBI TREATMENT           #####
  c("NT","LT"),
#####                        TB REACTIVATION RISK                          #####
  c("I1","I2","I3","I4"),
#####                          NON-TB MORTALITY                            #####
  c("M1","M2","M3","M4"),
#####                  LIVING CONDITIONS/RISK OF INFECTION                 #####
  c("L1","L2","L3"),
#####                              NATIVITY                                #####
  c("US","F1","F2")))
################################################################################
################################################################################
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
            c("TB_INF_LR","TB_INF_HR","TB_INF_F1","TB_INF_F2") # NEW TB INFECTIONS
)
length(ResNam)

###################################################

########################### HOW MANY VARS IN ResNam ############################
length(ResNam)
################################################################################
################################################################################
################################################################################

