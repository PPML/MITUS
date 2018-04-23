################################################################################
##### THE CODE BELOW WILL GATHER AND FORMAT INPUTS FOR THE TB_MODEL        #####
##### FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME; #####
##### VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.        #####
################################################################################
##### THIS FILE USES FUNCTIONS FOUND IN HELPER_FUNCTIONS.R                 #####
################################################################################
library(MASS)
source("R/basic_functions.R")
source("R/define_P.R")
load("data/ModelInputs_9-2-16.rData")
################################################################################
###########################          INPUTS            #########################
################################################################################
  BgMort           <- Inputs[["BgMort"]]
  InitPop          <- Inputs[["InitPop"]]
  Births           <- Inputs[["Births"]]
  ImmigInputs      <- Inputs[["ImmigInputs"]]
  TxInputs         <- Inputs[["TxInputs"]]

##########                PARAMETER DEFINITIONS                      ###########
##########                RISK FACTOR DISTRIBUTIONS   ##########################

  #dist   <- dist

#######################           BIRTHS                 #######################
####### INDEXED BY TIME, ABSOLUTE NUMBER OF NEW ADULT ENTRANTS OVER TIME #######

  Birthst   <- SmoCurve(Births)*P["TunBirths"]/12

##########################      MORTALITY RATES       ##########################
########################## BACKGROUND MORTALITY BY TIME ########################
  mubt      <- matrix(NA,1801,11)
  for(i in 1:11) {
  	mubt[,i] <- SmoCurve(BgMort[,i+1])*P["TunMubt"]/12
  }

#########################     DISEASE SPECIFIC       ###########################
#############    ACTIVE TB RATES DEFAULT TO THE SMEAR POS LEVELS   #############

  muIp  	  <- P["muIp"]/12

######################## MULTIPLER OF MORT RATE ABOVE ########################

  TunmuTbAg <- P["TunmuTbAg"]
  TunmuHvAg <- P["TunmuTbAg"]

############ CONVERT ANNUAL RATES OF RF MORTALITY TO MONTHLY RATES ##########

############ THESE MUST BE UPDATED
  muRF1      <- P["muH1"]/12*.01
  muRF2      <- (P["muH2"]+P["muH1"])/24 *.01
  muRF3      <- P["muH2"]/12*.01
  muTbRF    <- P["muTbH"]/12

###############  RATE RATIO OF MORTALITY INCREASE FOR HIGH RISK ###############

  RRmuHR    <- c(1,P["RRmuHR"])

############### CREATE A MATRIX OF RF MORTALITIES BY AGE GROUP ###############

  # vRFMort    <- c(0,0,0,0);
  # names(vRFMort) <- c("RF1","RF2","RF3","RF4")
  # vRFMort[2] <- (1/3)*muRF
  # vRFMort[3] <- (2/3)*muRF
  # vRFMort[4] <- muRF

  vRFMort    <- matrix(0,11,4);
  rownames(vRFMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vRFMort) <- c("RF1","RF2","RF3","RF4")
  vRFMort[,1] <- 0;
  vRFMort[,2] <- muRF1*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  vRFMort[,3] <- muRF2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  vRFMort[,4] <- muRF3*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)


############### CREATE A MATRIX OF TB MORTALITIES BY AGE GROUP ###############

  vTMort   <- matrix(0,11,6);
  rownames(vTMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vTMort) <- c("Su","Sp","Ls","Lf","Ac","Tx")
  vTMort[,5] <- muIp #active disease rates default to smear positive
  RRmuTbAg <- exp(c(0,0,1:9)*TunmuTbAg)
  for(i in 1:ncol(vTMort)) {
    vTMort[,i] <- vTMort[,i] * RRmuTbAg
  }

######################         IMMIGRATION             ########################
######################         OVERALL IMM.            ########################

  TotImmig0       <- (c(ImmigInputs[[1]][1:151])+c(rep(0,65),cumsum(rep(P["ImmigVolFut"],86))))/12*P["ImmigVol"]
  TotImmig1       <- TotImmig0
  TotImmig        <- SmoCurve(TotImmig1)
  TotImmAge       <- outer(TotImmig,ImmigInputs[["AgeDist"]])

######################           LTBI IMM.             ########################
  PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:65]^P["TunLtbiTrend"]*ImmigInputs[["PrevTrend25_34"]][65]^(1-P["TunLtbiTrend"]),
                        ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutLat"]/0.99)^(1:86))
  set.seed(982378)
  #set.seed( P["rand_seed"]+1)
  # PrevTrend25_341l <-   c(PrevTrend25_340l[1:66],exp(mvrnorm(1, log(PrevTrend25_340l[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  PrevTrend25_341l <-   PrevTrend25_340l
  PrevTrend25_34l  <- SmoCurve(PrevTrend25_341l)
  PrevTrend25_34_ls <- (PrevTrend25_34l);
  PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
  ImmLat          <- matrix(NA,length(PrevTrend25_34_ls),11)
  for(i in 1:11) ImmLat[,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]*P["LtbiPar1"]-(c(2.5,1:9*10,100)/100)[i]^2*P["LtbiPar2"])*PrevTrend25_34_ls))*TotImmAge[,i]

######################         ACTIVE TB IMM.           ########################
  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:65],ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutAct"]/0.99)^(1:86))
#  PrevTrend25_341a <-   c(PrevTrend25_340a[1:66],exp(mvrnorm(1, log(PrevTrend25_340a[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  PrevTrend25_341a <-   PrevTrend25_340a
  PrevTrend25_34a  <- SmoCurve(PrevTrend25_341a)
  ImDxChngV      <- SmoCurve(c(rep(1,57),seq(1,P["ImDxChng"],length.out=6)[-1],rep(P["ImDxChng"],89)))
  ImmAct         <- outer(PrevTrend25_34a*P["RRtbprev"]*ImDxChngV,ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*P["pImAct"]
  ImmFst         <- outer(PrevTrend25_34a*P["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-P["pImAct"])
  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat

######################   EXOGENEOUS INFECTION RISK      ########################

  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  ExogInf        <- P["ExogInf"]*PrevTrend25_34a/PrevTrend25_341a["2013"]/12

#removed *(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])

######################             EMIGRATION          #########################

   rEmmigFB <- c(P["rEmmigF1"],P["rEmmigF2"])/12

######################       HIGH-RISK ENTRY/EXIT      ########################

  p_HR     <- P["pHR"]
  yr       <- c(2.5,1:10*10)
  r0_5     <- 1/3; r45_55 <- 1/20
  HR_exit  <- r0_5*((r45_55/r0_5)^(1/(50-2.5)))^(yr-2.5)
  HR_entry <- HR_exit*p_HR*1.3
  HrEntEx  <- cbind(HR_entry,HR_exit)/12

######################       TB TRANSMISSION           #######################

  CR           <- P["CR"]/12
  RelInfRg     <- c(1.0,P["RelCrHr"])*CR
  TunTbTransTx <- .1#P["TunTbTransTx"]  # set to zero?
  Vmix         <- 1-c(P["sigmaHr"],P["sigmaFb"])
  RelInf       <- rep(0,6)
  names(RelInf) <- c("Su","Sp","Ls","Lf","Ac", "Tx")
  RelInf[5] <- 1; #set to 1 as what it was in the old model
  RelInf[6] <- RelInf[5]*TunTbTransTx

######################      TB NATURAL HISTORY       ##########################
######################        EARLY EPIDEMIC         ##########################

  Early0 <- P["Early0"]
  EarlyTrend <- c(rep(1+Early0,200*12),seq(1+Early0,1.0,length.out=50*12+2))

######################     PROGRESSION TO DISEASE     ##########################

  pfast      <- P["pfast"]
  ORpfast1   <- P["ORpfast1"] ## age group 1
  ORpfast2   <- P["ORpfast2"] ## age group 2
  ORpfastRF  <- P["ORpfastH"] ##riskfactor
  ORpfastPI  <- P["ORpfastPI"]
  rslow      <- P["rslow"]/12
  rslowRF    <- P["rslowH"]/12
  RRrslowRF  <- rslowRF/rslow
  rfast      <- P["rfast"]/12
  rrSlowFB0  <- P["rrSlowFB"]
  rrSlowFB   <- c(1,rrSlowFB0,rrSlowFB0)
##############            ORIGINAL Mpfast[ag][hv]             ################
##############          CREATE NEW Mpfast[ag][im]               ##############
############## MIGHT WRITE A NEW SCRIPT FOR THIS PART
  Mpfast       <- matrix(NA,11,4)
############## CREATE AN ODDS FROM THE PROB OF FAST PROGRESSION ##############
  Mpfast[,]    <- pfast/(1-pfast)
  Mpfast[1,]   <- Mpfast[1,]*ORpfast1 # progression for age group 1
  Mpfast[2,]   <- Mpfast[2,]*ORpfast2 # progression for age group 2
############ UPDATE PROBS FOR LEVEL 2 OF REACTIVATION ###########
  Mpfast[,2]   <- (ORpfastRF*1/3)*Mpfast[,2]
############ UPDATE PROBS FOR LEVEL 3 OF REACTIVATION ###########
  Mpfast[,3]   <- (ORpfastRF*2/3)*Mpfast[,3]
############ UPDATE PROBS FOR LEVEL 4 OF REACTIVATION ###########
  Mpfast[,4]   <- Mpfast[,4]*ORpfastRF #progression for tb reactivation group 4

#################       CREATE A NEW MATRIX PARTIAL. IMM.     #################
  MpfastPI     <- Mpfast
  MpfastPI[,1] <- Mpfast[,1]*ORpfastPI
  MpfastPI[,2] <- (ORpfastRF*1/3)*MpfastPI[,1]
  MpfastPI[,3] <- (ORpfastRF*2/3)*MpfastPI[,1]

### ADD IN THE INV LOGIT FOR THE OTHER FACTOR LEVELS

##### UPDATE BOTH MATRICES WITH PROBABILITIES, NOT RATES
  Mpfast[,]    <- Mpfast[,]  /(1+Mpfast[,]);
  MpfastPI[,]  <- MpfastPI[,]/(1+MpfastPI[,]);

############# CREATE A VECTOR FOR RATE OF SLOW PROGRESSION THAT WILL
############# VARY BASED ON LEVELS OF TB REACTIVATION RATES
  Vrslow     <- rep(rslow,4)
############# UPDATE LEVEL FOUR OF THE RATE OF SLOW BASED ON CALCULATED RR FROM
############# USER INPUTTED RR FOR THE RISK FACTOR
  for (i in 2:4){
    Vrslow[i]=Vrslow[i]*RRrslowRF*((i-1)/3)
  }

  TunrslowAge  <- P["TunrslowAge"]
  rrReactAg       <- exp(c(0,0,0,0,0,0,0.5,1:4)*P["TunrslowAge"])
  Mrslow <- outer(rrReactAg,Vrslow)

#######################       RATE OF RECOVERY          ########################

  rRecov     <-  P["rRecov"]/12

#######################       RATE OF SELF CURE         ########################

  rSlfCur      <- P["rSlfCur"]/12

######################          LTBI DIAGNOSIS           ########################
#  rLtScrt       <- LgtCurve(1985,2015,P["rLtScr"])/12
  rLtScrt       <- c(rep(0,888),LgtCurve(1985,2015,P["rLtScr"])/12)
  SensLt        <- P["SensLt"]    #  sens of test for latent TB infection (based on IGRA QFT-GIT)
  SpecLt        <- P["SpecLt"]    #  spec of test for latent TB infection (based on IGRA QFT-GIT)
  SpecLtFb      <- SpecLt         #  spec of test for latent TB infection (based on IGRA QFT-GIT) in foreign-born (assumed BCG exposed)
###########   WILL THIS PARAMETER NEED TO BE REDUCED?
  rrTestHr      <- P["rrTestHr"] # RR of LTBI screening for HIV and HR as cmpared to general
  rrTestLrNoTb  <- P["rrTestLrNoTb"] # RR of LTBI screening for individuals with no risk factors
  dLt           <- 1/9

  rDefLt        <- dLt*P["pDefLt"]/(1-P["pDefLt"])  # based on 50% tx completion with 6 mo INH regimen 2.0 [1.0,3.0] from Menzies Ind J Med Res 2011
  EffLt         <- P["EffLt"]
  LtTxPar       <- c(dLt,rDefLt,EffLt)

  LtDxPar <- matrix(NA,3,2);
  colnames(LtDxPar) <- c("latent","no latent");
  rownames(LtDxPar) <- c("LR","HR","FB")
  LtDxPar[,1] <- c(SensLt                 , rrTestHr*SensLt    , SensLt)
  LtDxPar[,2] <- c(rrTestLrNoTb*(1-SpecLt), rrTestHr*(1-SpecLt), (1-SpecLtFb))

  pImmScen    <- P["pImmScen"] # lack of reactivitiy to IGRA for Sp

################################################################################
#######################         TB DIAGNOSIS            ########################
#######################      TEST CHARACTERISTICS       ########################

  SensSp    <- P["SensSp"]

######################           PROVIDER DELAY         ########################

  DelaySp    <- P["DelaySp"]

#######################     PROVIDER DELAY RR ELDERLY     ########################

  TunRRdxAge    <- P["TunRRdxAge"]
  RRdxAge       <- 1+(c(rep(1,6),cumprod(seq(1.05,1.3,length.out=5)))-1)*TunRRdxAge

#######################         ATTENDANCE RATE           ########################

  n_Spln   <- 5;
  n_Stps   <- 2010-1950+1;
  dif_pen   <- 1 # quadratic spline
# Working...
  x1    <- seq(1,n_Stps);
  k1  <- seq(min(x1),max(x1),length=n_Spln-dif_pen)
  dk1   <- k1[2]-k1[1] ;
  k1  <- c(k1[1]-dk1*((dif_pen+1):1),k1,k1[n_Spln-dif_pen]+dk1*(1:(dif_pen+1)))
  SpMat <- matrix(NA,nrow=n_Stps,ncol=n_Spln);

  for(i in 1:n_Spln){
    SpMat[,i] <- bspline(x=x1,k=k1,m=dif_pen,i=i)
  }

  DxPri          <- (seq(0.5,4,length.out=5)-0.75)*3.5/2.626+0.25  # Prior from 0.5 t to 4.0
  rDx            <- c(SpMat%*%(DxPri+c(P["Dx1"],P["Dx2"],P["Dx3"],P["Dx4"],P["Dx5"])),62:151)
  rDx[62:151]    <- rDx[61] + (rDx[61]-rDx[60])*cumsum((0.75^(1:90)))
  rDxt0          <- SmoCurve(rDx)/12;
  rDxt1          <- cbind(rDxt0,rDxt0)

# Put it all together
  rDxt           <- 1/(1/rDxt1+DelaySp)*SensSp
  rDxt[,2]       <- (rDxt[,1]-min(rDxt[,1]))/P["rrDxH"]+min(rDxt[,1]) #check this with Nick
  colnames(rDxt) <- c("Active","Active_HighRisk")

################################################################################
###########################     TREATMENT OUTCOMES    ##########################
################################################################################

  TunTxMort	<- P["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

###########################      REGIMEN DURATION      ##########################

  d1st <- 1/9

###########################       REGIMEN EFFICACY      ##########################

  pCurPs  <- P["pCurPs"]    # probability of cure with pansensitive TB, 1st line regimen (Menzies 09)

###########################        REGIMEN DEFAULT       ##########################

  rDef0         <- rep(NA,151)
  rDef0[1:30]   <- P["TxDefEarly"]
  rDef0[44:63]  <- ORAdd(TxInputs[[1]][,2],P["TunTxDef"])
  rDef0[64:151] <- rDef0[63]
  rDef1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=rDef0[-(31:43)],spar=0.4),x=1950:2100)$y
  rDeft         <- SmoCurve(rDef1)/12;
  rDeftH        <- rDeft*P["RRdefHR"] # second col is HR default rate

#########################        REGIMEN QUALITY       ##########################

  TxQual0         <- rep(NA,151)
  TxQual0[1:30]   <- P["TxQualEarly"]
  TxQual0[44:62]  <- ORAdd(TxInputs[[2]][,2],P["TunTxQual"])
  TxQual0[63:151] <- TxQual0[62]
  TxQual1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=TxQual0[-(31:43)],spar=0.4),x=1950:2100)$y
  TxQualt         <- SmoCurve(TxQual1);
  RRcurDef        <- P["RRcurDef"]

#########################         RETREATMENT         ##########################

   pReTx   <- c(rep(0,888),LgtCurve(1985,2000,P["pReTx"]))   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)

#####################         NEW TB TREATMENT VECTOR       ####################

  TxVec           <- rep(NA,2)
  names(TxVec) <- c("TxCompRate","TxEff")
  TxVec[1]       <-  d1st
  TxVec[2]       <-  pCurPs

################################################################################
##### CREATE A LIST TO HOLD THE VECTORS
################################################################################
################################################################################
StatList <- noquote(list(
#####					                  AGE CATEGORIES                             #####
#####           AGE GROUPS 0-4, 5-14, 15-24, ... , 95-94, 95+              #####
  c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p"),
#####					     				    TUBERCULOSIS STATES		              			   #####
##### SUSCEPTIBLE; UNINFECTED & PARTIALLY IMMUNE; LATENT SLOW; LATENT FAST;#####
#####       ACTIVE TB SMEAR NEG.; ACTIVE TB SMEAR POS.; TB TREATMENT       #####
  c("Su","Sp","Ls","Lf", "Ac", "Tx"),
#####                          TREATMENT HISTORY                           #####
#####              NO TB TREATMENT HISTORY, PRIOR LTBI TREATMENT           #####
  c("NT","LT"),
#####                        TB REACTIVATION RISK                          #####
  c("I1","I2","I3","I4"),
#####                          NON-TB MORTALITY                            #####
  c("M1","M2","M3","M4"),
#####                  LIVING CONDITIONS/RISK OF INFECTION                 #####
  c("LR","HR"),
#####                              NATIVITY                                #####
  c("US","F1","F2")))
################################################################################
################################################################################
ResNam <- c("Year",                                         # year
            ##############################    POPULATION    ################################
            "N_ALL",                                        # total pop
            paste("N",StatList[[1]],sep="_"),               # pop by ag cat
            paste("N",StatList[[2]],sep="_"),               # pop by tb cat
            paste("N",StatList[[4]],sep="_"),               # pop by im cat
            paste("N",StatList[[5]],sep="_"),               # pop by nm cat
            paste("N",StatList[[6]],sep="_"),               # pop by rg cat
            paste("N",StatList[[7]],sep="_"),               # pop by na cat
            ###################    POPULATION W/ NATIVITY  ################################
            paste("N_US",StatList[[1]],sep="_"),            # US pop by ag cat
            paste("N_FB",StatList[[1]],sep="_"),            # FB pop by ag cat
            paste("N_US_LTBI",StatList[[1]],sep="_"),       # US LTBI pop by ag cat
            paste("N_FB_LTBI",StatList[[1]],sep="_"),       # FB LTBI pop by ag cat
            paste("N_RF",StatList[[1]],sep="_"),            # RF pop by ag cat
            #############################    MORTALITY    ################################
            paste("TBMORT_RFNEG",StatList[[1]],sep="_"),   # TB mort, HIV neg, by ag cat
            paste("TBMORT_RFPOS",StatList[[1]],sep="_"),   # TB mort, HIV pos, by ag cat
            paste("RFMORT",StatList[[1]],sep="_"),         # RF mort, by ag cat
            paste("TOTMORT",StatList[[1]],sep="_"),         # total mort, by ag cat
            #### @object 131
            #############################   TX OUTCOMES   ################################

            "TBTX_COMPLT","TBTX_DISCONT","TBTX_DIED",       # TB treatment outcomes complete, discontinue, death
            "NOTIF_ALL",                                    # total notif
            paste("NOTIF",StatList[[1]],sep="_"),           # notif by ag cat
            paste("NOTIF",StatList[[7]],sep="_"),         # notif by nat cat
            #            paste("NOTIF_HIV",c("POS","NEG"),sep="_"),      # notif by HIV pos/neg
            paste("NOTIF",StatList[[6]],sep="_"),           # notif by rg cat
            #           paste("NOTIF_US_N",StatList[[3]],sep="_"),      # notif, US, N, by dr cat
            #            paste("NOTIF_US_E",StatList[[3]],sep="_"),      # notif, US, E, by dr cat
            #           paste("NOTIF_FB_N",StatList[[3]],sep="_"),      # notif, FB, N, by dr cat
            #           paste("NOTIF_FB_E",StatList[[3]],sep="_"),      # notif, FB, E, by dr cat
            ###########################   TLTBI INITIATION   ##############################
            "TLTBI_INITS",                                  # Initiations on LTBI tx
            "TLTBI_INITS_FB",                               # Initiations on LTBI tx FB
            "TLTBI_INITS_HR",                               # Initiations on LTBI tx HR
            #            "TLTBI_INITS_HV",                               # Initiations on LTBI tx HV
            "TLTBI_INITS_TP",                               # Initiations on LTBI tx, with LTBI
            ###########################     TB INCIDENCE     ##############################
            "INCID_ALL",                                    # Total incidence
            paste("INCID_ALL",StatList[[1]],sep="_"),       # Total incidence by ag cat
            "INCID_ALL_US",                                 # Total incidence, US born
            "INCID_ALL_FB",                                 # Total incidence, foreign born
            "INCID_ALL_FB2",                                # Total incidence, foreign born
            "INCID_ALL_HR",                                 # Total incidence, high risk
            #            "INCID_ALL_HV",                                 # Total incidence, HIV pos
            "INCID_REC",                                    # Total incidence, recent infection
            paste("INCID_REC",StatList[[1]],sep="_"),       # Total incidence by ag cat, recent infection
            "INCID_REC_US",                                 # Total incidence, US born, recent infection
            "INCID_REC_FB",                                 # Total incidence, foreign born, recent infection
            "INCID_REC_FB2",                                # Total incidence, foreign born, recent infection
            "INCID_REC_HR",                                 # Total incidence, high risk, recent infection
            #            "INCID_REC_HV",                                 # Total incidence, HIV pos, recent infection
            ###########################    NOTIFICATION DEAD      ##############################
            "NOTIF_MORT_ALL",                               # total notif, dead at diagnosis
            paste("NOTIF_MORT",StatList[[1]],sep="_"),      # notif by ag cat, dead at diagnosis
            paste("NOTIF_MORT",StatList[[7]],sep="_"),      # notif by nat cat, dead at diagnosis
            #            paste("NOTIF_MORT_HIV",c("POS","NEG"),sep="_"), # notif by HIV pos/neg, dead at diagnosis
            paste("NOTIF_MORT",StatList[[6]],sep="_"),      # notif by rg cat, dead at diagnosis
            paste("NOTIF_US",StatList[[1]],sep="_"),        # notif by ag cat, US only
            paste("NOTIF_US_MORT",StatList[[1]],sep="_"),   # notif by ag cat, dead at diagnosis US only
            #            paste("NOTIF_MORT_HIV_Neg",StatList[[1]],sep="_"),    # notif by ag cat, dead at diagnosis HIV neg
            #            paste("TOTMORT_W_HIV",StatList[[1]],sep="_"),   # total mort, by ag cat, have HIV
            paste("TOTMORT_W_TB",StatList[[1]],sep="_"),     # total mort, by ag cat, have active TB
            c("N_Ls_US","N_Lf_US","N_Act_US"),
            c("N_Ls_FB","N_Lf_FB","N_Act_FB"),
            c("FOI_LR_US","FOI_HR_US","FOI_LR_FB","FOI_HR_FB" ), #force of infection
            c("TB_INF_LR","TB_INF_HR","TB_INF_US","TB_INF_F1","TB_INF_F2"), # NEW TB INFECTIONS
            c("TBMORT_US","TBMORT_NUS") ,          # tbmortality by nativity
            paste("TOTMORT_US",StatList[[1]],sep="_"),       # mortality by nat cat
            paste("TOTMORT_NUS",StatList[[1]],sep="_"),       # mortality by nat cat

            paste("TOTMORT_US",StatList[[4]],sep="_"),       # mortality by nat & im cat
            paste("TOTMORT_NUS",StatList[[4]],sep="_"),       # mortality by nat & im cat

            paste("TOTMORT_US",StatList[[5]],sep="_"),       # mortality by nat & nm cat
            paste("TOTMORT_NUS",StatList[[5]],sep="_"),       # mortality by nat & nm cat

            paste("TOTMORT_US",StatList[[6]],sep="_"),
            paste("TOTMORT_NUS",StatList[[6]],sep="_"),


            paste("N_US",StatList[[4]],sep="_"),               # pop by nat and im cat
            paste("N_NUS",StatList[[4]],sep="_"),               # pop by nat and im cat

            paste("N_US",StatList[[5]],sep="_"),               # pop by nat and nm cat
            paste("N_NUS",StatList[[5]],sep="_"),               # pop by nat and nm cat

            paste("N_US",StatList[[6]],sep="_"),               # pop by nat and hr cat
            paste("N_NUS",StatList[[6]],sep="_"),              # pop by nat and hr cat
            paste("TOTMORT"),
            paste("N_NM1",StatList[[4]],sep="_"),
            paste("N_NM2",StatList[[4]],sep="_"),
            paste("N_NM3",StatList[[4]],sep="_"),
            paste("N_NM4",StatList[[4]],sep="_")

)
length(ResNam)

###################################################

########################### HOW MANY VARS IN ResNam ############################
length(ResNam)
################################################################################
################################################################################
################################################################################

