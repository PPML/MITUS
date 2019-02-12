#'THE CODE BELOW SOURCES THE MODEL INPUTS AND THEN FORMATS THEM FOR USE
#'IN THE OUTPUTSINTZ FUNCTION THAT CALLS CSIM FROM THE TB_MODEL.CPP
#'FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME
#'VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.

#'@name param_noRB
#' @param P vector of
#' @return Params list
#' @export
param_noRB <- function (PV){
  # load("~/MITUS/data/US_ModelInputs_9-6-18.rda")
  ################################################################################
  ###########################          INPUTS            #########################
  ################################################################################
  BgMort           <- Inputs[["BgMort"]]
  InitPop          <- Inputs[["InitPop"]]
  Births           <- Inputs[["Births"]]
  ImmigInputs      <- Inputs[["ImmigInputs"]]
  TxInputs         <- Inputs[["TxInputs"]]
  NetMig           <- Inputs[["NetMigrState"]]

  ########## DEFINE A VARIABLE THAT WILL DETERMINE HOW LONG THE TIME DEPENDENT
  ########## VARIABLES SHOULD BE (IN MONTHS)
  month<-1201;

  ##########                PARAMETER DEFINITIONS                      ###########
  ##########                RISK FACTOR DISTRIBUTIONS   ##########################

  #adj_fact<-exp(.0001*(10:0)/11 + .0011*(0:10)/11)
  adj_fact <- exp(PV[["adj_ag1"]]*(10:0)/11 + PV[["adj_ag11"]]*(0:10)/11)

  #######################           BIRTHS                 #######################
  ####### INDEXED BY TIME, ABSOLUTE NUMBER OF NEW ADULT ENTRANTS OVER TIME #######

  Birthst   <- SmoCurve(Births)*PV["TunBirths"]/12
  Birthst   <- Birthst[1:month]

  ##########################      MORTALITY RATES       ##########################
  ########################## BACKGROUND MORTALITY BY TIME ########################
  mubt      <- matrix(NA,1801,11)

  TunmuAg <- PV["TunmuAg"]
  RRmuAg <- exp((1:11)*TunmuAg)
  for(i in 1:11) {
    mubt[,i] <- SmoCurve(BgMort[,i+1])*PV["TunMubt"]/12
    mubt[,i] <- mubt[,i]*RRmuAg[i]
  }



  # for(i in 2:10) {
  #   mubt[,i]<-mubt[,i]*exp(TunmuAg)
  # }
  #########################     DISEASE SPECIFIC       ###########################
  #############    ACTIVE TB RATES DEFAULT TO THE SMEAR POS LEVELS   #############

  muIp  	  <- PV["muIp"]/12

  ######################## MULTIPLER OF MORT RATE ABOVE ########################

  TunmuTbAg <- PV["TunmuTbAg"]

  ###############  RATE RATIO OF MORTALITY INCREASE FOR HIGH RISK ###############

  RRmuHR    <- c(1,PV["RRmuHR"])

  ############### CREATE A MATRIX OF RF MORTALITIES BY AGE GROUP ###############
  ############### CREATE A MATRIX OF RF MORTALITIES BY AGE GROUP ###############

  RF_fact=20

  RRmuRF    <- rep(NA,4);
  names(RRmuRF) <- c("RF1","RF2","RF3","RF4")
  RRmuRF<-c(1,1,1,1)


  # vRFMort    <- matrix(0,11,4);
  # rownames(vRFMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  # colnames(vRFMort) <- c("RF1","RF2","RF3","RF4")
  # vRFMort[,1] <- 0;
  # vRFMort[,2] <- muRF1*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  # vRFMort[,3] <- muRF2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
  # vRFMort[,4] <- muRF3*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)

  ##### combine the two RRs to a single factor
  RRs_mu <-matrix(NA,2,4)
  RRs_mu[1,] <-RRmuRF
  RRs_mu[2,] <-RRmuRF*RRmuHR[2]
  ############### CREATE A MATRIX OF TB MORTALITIES BY AGE GROUP ###############

  vTMort   <- matrix(0,11,6);
  rownames(vTMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vTMort) <- c("Su","Sp","Ls","Lf","Ac","Tx")
  vTMort[,5:6] <- muIp #active disease rates default to smear positive
  RRmuTbAg <- exp(c(0,0,1:9)*TunmuTbAg)
  for(i in 1:ncol(vTMort)) {
    vTMort[,i] <- vTMort[,i] * RRmuTbAg
  }

  ######################         IMMIGRATION             ########################
  ######################         OVERALL IMM.            ########################

  TotImmig0       <- (c(ImmigInputs[[1]][1:151])+c(rep(0,65),cumsum(rep(PV["ImmigVolFut"],86))))/12*PV["ImmigVol"]
  TotImmig1       <- TotImmig0
  TotImmig        <- SmoCurve(TotImmig1)
  TotImmAge       <- outer(TotImmig,ImmigInputs[["AgeDist"]])

  ######################           LTBI IMM.             ########################
  PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:65]^PV["TunLtbiTrend"]*ImmigInputs[["PrevTrend25_34"]][65]^(1-PV["TunLtbiTrend"]),
                        ImmigInputs[["PrevTrend25_34"]][66:151]*(PV["ImmigPrevFutLat"]/0.99)^(1:86))
  set.seed(982378)
  #set.seed( PV["rand_seed"]+1)
  # PrevTrend25_341l <-   c(PrevTrend25_340l[1:66],exp(mvrnorm(1, log(PrevTrend25_340l[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  PrevTrend25_341l <-   PrevTrend25_340l
  PrevTrend25_34l  <- SmoCurve(PrevTrend25_341l)
  PrevTrend25_34_ls <- (PrevTrend25_34l);
  PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
  ImmLat          <- matrix(NA,length(PrevTrend25_34_ls),11)
  for(i in 1:11) ImmLat[,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]*PV["LtbiPar1"]-(c(2.5,1:9*10,100)/100)[i]^2*PV["LtbiPar2"])*PrevTrend25_34_ls))*TotImmAge[,i]

  ######################         ACTIVE TB IMM.           ########################
  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:65],ImmigInputs[["PrevTrend25_34"]][66:151]*(PV["ImmigPrevFutAct"]/0.99)^(1:86))
  #  PrevTrend25_341a <-   c(PrevTrend25_340a[1:66],exp(mvrnorm(1, log(PrevTrend25_340a[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
  PrevTrend25_341a <-   PrevTrend25_340a
  PrevTrend25_34a  <- SmoCurve(PrevTrend25_341a)
  # ImDxChngV      <- SmoCurve(c(rep(1,57),seq(1,PV["ImDxChng"],length.out=6)[-1],rep(PV["ImDxChng"],89)))
  ImmAct         <- outer(PrevTrend25_34a*PV["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*PV["pImAct"]
  ImmFst         <- outer(PrevTrend25_34a*PV["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-PV["pImAct"])
  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat

  ###################### TRUNCATE THESE VALS
  ImmAct<-ImmAct[1:month,];ImmFst<-ImmFst[1:month,]; ImmLat<-ImmLat[1:month,]; ImmNon<-ImmNon[1:month,]
  ######################   EXOGENEOUS INFECTION RISK      ########################

  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  ExogInf        <- PV["ExogInf"]*PrevTrend25_34a/PrevTrend25_341a["2013"]/12

  #removed *(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])

  ######################             EMIGRATION          #########################

  rEmmigFB <- c(PV["rEmmigF1"],PV["rEmmigF2"])/12

  ######################          NET MIGRATION          #########################

  net_mig_usb  <- (NetMig[,"usb" ]*PV["TunNetMig"])^(1/12)-1
  net_mig_nusb <- (NetMig[,"nusb"]*PV["TunNetMig"])^(1/12)-1

  ######################       HIGH-RISK ENTRY/EXIT      ########################

  p_HR     <- PV["pHR"]
  yr       <- c(2.5,1:10*10)
  r0_5     <- 1/3; r45_55 <- 1/20
  HR_exit  <- r0_5*((r45_55/r0_5)^(1/(50-2.5)))^(yr-2.5)
  HR_entry <- HR_exit*p_HR*1.3
  HrEntEx  <- cbind(HR_entry,HR_exit)/12

  ######################       TB TRANSMISSION           #######################

  CR           <- PV["CR"]/12
  RelInfRg     <- c(1.0,PV["RelCrHr"],1.0,PV["RelCrHr"])*CR
  TunTbTransTx <- 0#PV["TunTbTransTx"]  # set to zero?
  Vmix         <- 1-c(PV["sigmaHr"],PV["sigmaFb"])
  RelInf       <- rep(0,6)
  names(RelInf) <- c("Su","Sp","Ls","Lf","Ac","Tx")
  RelInf[5] <- 1; #set to 1 as what it was in the old model
  RelInf[6] <- RelInf[5]*TunTbTransTx

  ######################      TB NATURAL HISTORY       ##########################
  ######################        EARLY EPIDEMIC         ##########################

  Early0 <- PV["Early0"]
  EarlyTrend <- c(rep(1+Early0,200*12),seq(1+Early0,1.0,length.out=50*12+2))

  ######################     PROGRESSION TO DISEASE     ##########################

  pfast      <- PV["pfast"]
  ORpfast1   <- PV["ORpfast1"] ## age group 1
  ORpfast2   <- PV["ORpfast2"] ## age group 2
  ORpfastRF  <- PV["ORpfastH"] ##riskfactor
  ORpfastPI  <- PV["ORpfastPI"]

  ##############            ORIGINAL Mpfast[ag][hv]             ################
  ##############          CREATE NEW Mpfast[ag][im]               ##############
  ############## MIGHT WRITE A NEW SCRIPT FOR THIS PART
  Mpfast       <- matrix(NA,11,4)
  ############## CREATE AN ODDS FROM THE PROB OF FAST PROGRESSION ##############
  Mpfast[,]    <- pfast/(1-pfast)
  Mpfast[1,]   <- Mpfast[1,]*ORpfast1 # progression for age group 1
  Mpfast[2,]   <- Mpfast[2,]*ORpfast2 # progression for age group 2


  #################       CREATE A NEW MATRIX PARTIAL. IMM.     #################
  MpfastPI     <- Mpfast

  #vector of ORpfastRF
  vORpfastPIRF<-vORpfastRF  <-c(1,1,1,1)
  vORpfastRF  <-(exp((0:3)/3*log(ORpfastRF)))
  vORpfastPIRF  <- vORpfastRF*ORpfastPI

  ############ UPDATE PROBS FOR LEVEL 2 OF REACTIVATION ###########
  Mpfast[,2]   <- vORpfastRF[2]*Mpfast[,2]
  MpfastPI[,2]   <- vORpfastPIRF[2]*Mpfast[,2]
  ############ UPDATE PROBS FOR LEVEL 3 OF REACTIVATION ###########
  Mpfast[,3]   <- vORpfastRF[3]*Mpfast[,3]
  MpfastPI[,3]   <- vORpfastPIRF[3]*Mpfast[,3]
  ############ UPDATE PROBS FOR LEVEL 4 OF REACTIVATION ###########
  Mpfast[,4]   <- vORpfastRF[4]*Mpfast[,4]
  MpfastPI[,4]   <- vORpfastPIRF[4]*Mpfast[,4]

  ##### UPDATE BOTH MATRICES WITH PROBABILITIES, NOT RATES
  Mpfast[,]    <- Mpfast[,]  /(1+Mpfast[,]);
  MpfastPI[,]  <- MpfastPI[,]/(1+MpfastPI[,]);

  rslow      <- PV["rslow"]/12
  rslowRF    <- PV["rslowH"]/12
  RRrslowRF  <- rslowRF/rslow
  rfast      <- PV["rfast"]/12
  # rrSlowFB0  <- PV["rrSlowFB"]
  rrSlowFB   <- c(1,1,1)

  ############# CREATE A VECTOR FOR RATE OF SLOW PROGRESSION THAT WILL
  ############# VARY BASED ON LEVELS OF TB REACTIVATION RATES
  Vrslow     <- rep(1,4)
  ############# UPDATE LEVEL FOUR OF THE RATE OF SLOW BASED ON CALCULATED RR FROM
  ############# USER INPUTTED RR FOR THE RISK FACTOR
  Vrslow<-rslow*exp((0:3)/3*log(RRrslowRF))

  TunrslowAge  <- PV["TunrslowAge"]
  rrReactAg       <- exp(c(0,0,0,0,0,0,0.5,1:4)*PV["TunrslowAge"])
  Mrslow <- outer(rrReactAg,Vrslow)

  #######################       RATE OF RECOVERY          ########################

  rRecov     <-  PV["rRecov"]/12

  #######################       RATE OF SELF CURE         ########################

  rSlfCur      <- PV["rSlfCur"]/12

  ######################          LTBI DIAGNOSIS           ########################
  rLtScrt       <- LgtCurve(1985,2015,PV["rLtScr"])/12
  # rLtScrt       <- c(rep(0,888),LgtCurve(1985,2015,PV["rLtScr"])/12)
  SensLt        <- PV["SensLt"]    #  sens of test for latent TB infection (based on IGRA QFT-GIT)
  SpecLt        <- PV["SpecLt"]    #  spec of test for latent TB infection (based on IGRA QFT-GIT)
  SpecLtFb      <- SpecLt         #  spec of test for latent TB infection (based on IGRA QFT-GIT) in foreign-born (assumed BCG exposed)
  ###########   WILL THIS PARAMETER NEED TO BE REDUCED?
  rrTestHr      <- PV["rrTestHr"] # RR of LTBI screening for HIV and HR as cmpared to general
  rrTestLrNoTb  <- PV["rrTestLrNoTb"] # RR of LTBI screening for individuals with no risk factors
  #  dLt           <- 1/9

  rDefLt        <- PV["pDefLt"]/(1-PV["pDefLt"])  # based on 50% tx completion with 6 mo INH regimen 2.0 [1.0,3.0] from Menzies Ind J Med Res 2011
  EffLt         <- PV["EffLt"]
  ######NEW PARAMETER FOR MITUS MODEL
  pTlInt        <- .80
  LtTxPar       <- c(pTlInt,PV["pDefLt"],EffLt)

  LtDxPar <- matrix(NA,3,2);
  colnames(LtDxPar) <- c("latent","no latent");
  rownames(LtDxPar) <- c("LR","HR","FB")
  LtDxPar[,1] <- c(SensLt                 , rrTestHr*SensLt    , SensLt)
  LtDxPar[,2] <- c(rrTestLrNoTb*(1-SpecLt), rrTestHr*(1-SpecLt), (1-SpecLtFb))

  pImmScen    <- PV["pImmScen"] # lack of reactivitiy to IGRA for Sp

  ################################################################################
  #######################         TB DIAGNOSIS            ########################
  #######################      TEST CHARACTERISTICS       ########################

  SensSp    <- PV["SensSp"]

  ######################           PROVIDER DELAY         ########################

  DelaySp    <- PV["DelaySp"]

  #######################     PROVIDER DELAY RR ELDERLY     ########################

  TunRRdxAge    <- PV["TunRRdxAge"]
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
  rDx            <- c(SpMat%*%(DxPri+c(PV["Dx1"],PV["Dx2"],PV["Dx3"],PV["Dx4"],PV["Dx5"])),62:151)
  rDx[62:151]    <- rDx[61] + (rDx[61]-rDx[60])*cumsum((0.75^(1:90)))
  rDxt0          <- SmoCurve(rDx)/12;
  rDxt1          <- cbind(rDxt0,rDxt0)

  # Put it all together
  rDxt           <- 1/(1/rDxt1+DelaySp)*SensSp
  rDxt[,2]       <- (rDxt[,1]-min(rDxt[,1]))/PV["rrDxH"]+min(rDxt[,1]) #check this with Nick
  colnames(rDxt) <- c("Active","Active_HighRisk")
  rDxt<-rDxt[1:month,]
  ################################################################################
  ###########################     TREATMENT OUTCOMES    ##########################
  ################################################################################

  TunTxMort	<- PV["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

  ###########################      REGIMEN DURATION      ##########################

  d1st <- 1/9

  ###########################       REGIMEN EFFICACY      ##########################

  pCurPs  <- PV["pCurPs"]    # probability of cure with pansensitive TB, 1st line regimen (Menzies 09)

  ###########################        REGIMEN DEFAULT       ##########################

  rDef0         <- rep(NA,151)
  rDef0[1:30]   <- PV["TxDefEarly"]
  rDef0[44:63]  <- ORAdd(TxInputs[[1]][,2],PV["TunTxDef"])
  rDef0[64:151] <- rDef0[63]
  rDef1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=rDef0[-(31:43)],spar=0.4),x=1950:2100)$y
  rDeft         <- SmoCurve(rDef1)/12;
  rDeftH        <- rDeft*PV["RRdefHR"] # second col is HR default rate
  rDeft<-rDeft[1:month]
  rDeftH<-rDeftH[1:month]
  #########################        REGIMEN QUALITY       ##########################

  TxQual0         <- rep(NA,151)
  TxQual0[1:30]   <- PV["TxQualEarly"]
  TxQual0[44:62]  <- ORAdd(TxInputs[[2]][,2],PV["TunTxQual"])
  TxQual0[63:151] <- TxQual0[62]
  TxQual1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=TxQual0[-(31:43)],spar=0.4),x=1950:2100)$y
  TxQualt         <- SmoCurve(TxQual1);

  TxQualt<-TxQualt[1:month]

  RRcurDef        <- PV["RRcurDef"]

  #########################         RETREATMENT         ##########################

  pReTx   <- LgtCurve(1985,2000,PV["pReTx"])   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)
  pReTx   <- pReTx[1:month]

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
  ResNam <- c("Year",                                         # year            1
              ##############################    POPULATION    ################################
              "N_ALL",                                        # total pop       2
              paste("N",StatList[[1]],sep="_"),               # pop by ag cat    3-13
              paste("N",StatList[[2]],sep="_"),               # pop by tb cat    14-19
              paste("N",StatList[[4]],sep="_"),               # pop by im cat   20-23
              paste("N",StatList[[5]],sep="_"),               # pop by nm cat   24-27
              paste("N",StatList[[6]],sep="_"),               # pop by rg cat   28-29
              paste("N",StatList[[7]],sep="_"),               # pop by na cat   30-32
              ###################    POPULATION W/ NATIVITY  ################################
              paste("N_US",StatList[[1]],sep="_"),            # US pop by ag cat        33
              paste("N_FB",StatList[[1]],sep="_"),            # FB pop by ag cat        44
              paste("N_US_LTBI",StatList[[1]],sep="_"),       # US LTBI pop by ag cat   55
              paste("N_FB_LTBI",StatList[[1]],sep="_"),       # FB LTBI pop by ag cat   66
              paste("N_RF",StatList[[1]],sep="_"),            # RF pop by ag cat        77
              #############################    MORTALITY    ################################
              paste("TBMORT_US",StatList[[1]],sep="_"),   # TB mort, HIV neg, by ag cat   88
              paste("TBMORT_NUS",StatList[[1]],sep="_"),   # TB mort, HIV pos, by ag cat   99
              paste("RFMORT",StatList[[1]],sep="_"),         # RF mort, by ag cat           110
              paste("TOTMORT",StatList[[1]],sep="_"),         # total mort, by ag cat       121
              #### @object 131
              #############################   TX OUTCOMES   ################################

              "TBTX_COMPLT","TBTX_DISCONT","TBTX_DIED",       # TB treatment outcomes complete, discontinue, death  132-134
              "NOTIF_ALL",                                    # total notif 135
              paste("NOTIF",StatList[[1]],sep="_"),           # notif by ag cat 136-146
              paste("NOTIF",StatList[[7]],sep="_"),         # notif by nat cat  147-149
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
              "INCID_ALL",                                    # Total incidence  #156
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
              "NOTIF_MORT_ALL",                               # total notif, dead at diagnosis (188)
              paste("NOTIF_MORT",StatList[[1]],sep="_"),      # notif by ag cat, dead at diagnosis
              paste("NOTIF_MORT",StatList[[7]],sep="_"),      # notif by nat cat, dead at diagnosis
              #            paste("NOTIF_MORT_HIV",c("POS","NEG"),sep="_"), # notif by HIV pos/neg, dead at diagnosis
              paste("NOTIF_MORT",StatList[[6]],sep="_"),      # notif by rg cat, dead at diagnosis
              paste("NOTIF_US",StatList[[1]],sep="_"),        # notif by ag cat, US only
              paste("NOTIF_US_MORT",StatList[[1]],sep="_"),   # notif by ag cat, dead at diagnosis US only
              #            paste("NOTIF_MORT_HIV_Neg",StatList[[1]],sep="_"),    # notif by ag cat, dead at diagnosis HIV neg
              #            paste("TOTMORT_W_HIV",StatList[[1]],sep="_"),   # total mort, by ag cat, have HIV

              paste("TOTMORT_W_TB",StatList[[1]],sep="_"),     # total mort, by ag cat, have active TB ends at 237
              c("N_Ls_US","N_Lf_US","N_Act_US"),
              c("N_Ls_FB","N_Lf_FB","N_Act_FB"),
              c("FOI_LR_US","FOI_HR_US","FOI_LR_FB","FOI_HR_FB" ), #force of infection
              c("TB_INF_LR","TB_INF_HR","TB_INF_US","TB_INF_F1","TB_INF_F2"), # NEW TB INFECTIONS
              c("TBMORT_US","TBMORT_NUS") ,          # tbmortality by nativity

              paste("TOTMORT_US",StatList[[1]],sep="_"),       # mortality by nat cat ends at 265
              paste("TOTMORT_NUS",StatList[[1]],sep="_"),       # mortality by nat cat

              paste("TOTMORT_US",StatList[[4]],sep="_"),       # mortality by nat & im cat
              paste("TOTMORT_NUS",StatList[[4]],sep="_"),       # mortality by nat & im cat

              paste("TOTMORT_US",StatList[[5]],sep="_"),       # mortality by nat & nm cat
              paste("TOTMORT_NUS",StatList[[5]],sep="_"),       # mortality by nat & nm cat

              paste("TOTMORT_US",StatList[[6]],sep="_"),
              paste("TOTMORT_NUS",StatList[[6]],sep="_"),
              paste("N_US",StatList[[4]],sep="_"),               # pop by nat and im cat
              paste("N_NUS",StatList[[4]],sep="_"),               # pop by nat and im cat

              paste("N_US",StatList[[6]],sep="_"),               # pop by nat and nm cat
              paste("N_NUS",StatList[[6]],sep="_"),               # pop by nat and nm cat

              paste("N_US",StatList[[5]],sep="_"),               # pop by nat and hr cat
              paste("N_NUS",StatList[[5]],sep="_"),              # pop by nat and hr cat
              paste("TOTMORT"),

              paste("N_NM1",StatList[[4]],sep="_"),
              paste("N_NM2",StatList[[4]],sep="_"),
              paste("N_NM3",StatList[[4]],sep="_"),
              paste("N_NM4",StatList[[4]],sep="_"),

              paste("%_0-4","NM1",StatList[[4]],sep="_"),
              paste("%_0-4","NM2",StatList[[4]],sep="_"),
              paste("0-4","NM3",StatList[[4]],sep="_"),
              paste("0-4","NM4",StatList[[4]],sep="_"),

              paste("5-14","NM1",StatList[[4]],sep="_"),
              paste("5-14","NM2",StatList[[4]],sep="_"),
              paste("5-14","NM3",StatList[[4]],sep="_"),
              paste("5-14","NM4",StatList[[4]],sep="_"),

              paste("15-24","NM1",StatList[[4]],sep="_"),
              paste("15-24","NM2",StatList[[4]],sep="_"),
              paste("15-24","NM3",StatList[[4]],sep="_"),
              paste("15-24","NM4",StatList[[4]],sep="_"),

              paste("25-34","NM1",StatList[[4]],sep="_"),
              paste("25-34","NM2",StatList[[4]],sep="_"),
              paste("25-34","NM3",StatList[[4]],sep="_"),
              paste("25-34","NM4",StatList[[4]],sep="_"),

              paste("35-44","NM1",StatList[[4]],sep="_"),
              paste("35-44","NM2",StatList[[4]],sep="_"),
              paste("35-44","NM3",StatList[[4]],sep="_"),
              paste("35-44","NM4",StatList[[4]],sep="_"),

              paste("45-54","NM1",StatList[[4]],sep="_"),
              paste("45-54","NM2",StatList[[4]],sep="_"),
              paste("45-54","NM3",StatList[[4]],sep="_"),
              paste("45-54","NM4",StatList[[4]],sep="_"),

              paste("55-64","NM1",StatList[[4]],sep="_"),
              paste("55-64","NM2",StatList[[4]],sep="_"),
              paste("55-64","NM3",StatList[[4]],sep="_"),
              paste("55-64","NM4",StatList[[4]],sep="_"),

              paste("65-74","NM1",StatList[[4]],sep="_"),
              paste("65-74","NM2",StatList[[4]],sep="_"),
              paste("65-74","NM3",StatList[[4]],sep="_"),
              paste("65-74","NM4",StatList[[4]],sep="_"),

              paste("75-84","NM1",StatList[[4]],sep="_"),
              paste("75-84","NM2",StatList[[4]],sep="_"),
              paste("75-84","NM3",StatList[[4]],sep="_"),
              paste("75-84","NM4",StatList[[4]],sep="_"),

              paste("85-94","NM1",StatList[[4]],sep="_"),
              paste("85-94","NM2",StatList[[4]],sep="_"),
              paste("85-94","NM3",StatList[[4]],sep="_"),
              paste("85-94","NM4",StatList[[4]],sep="_"),

              paste("95p","NM1",StatList[[4]],sep="_"),
              paste("95p","NM2",StatList[[4]],sep="_"),
              paste("95p","NM3",StatList[[4]],sep="_"),
              paste("95p","NM4",StatList[[4]],sep="_"),

              paste("mort_rate",StatList[[1]],sep="_" ),

              paste("N_ag_1",StatList[[5]],sep="_" ),
              paste("N_ag_2",StatList[[5]],sep="_" ),
              paste("N_ag_3",StatList[[5]],sep="_" ),
              paste("N_ag_4",StatList[[5]],sep="_" ),
              paste("N_ag_5",StatList[[5]],sep="_" ),
              paste("N_ag_6",StatList[[5]],sep="_" ),
              paste("N_ag_7",StatList[[5]],sep="_" ),
              paste("N_ag_8",StatList[[5]],sep="_" ),
              paste("N_ag_9",StatList[[5]],sep="_" ),
              paste("N_ag_10",StatList[[5]],sep="_" ),
              paste("N_ag_11",StatList[[5]],sep="_" )
  )

  Params<-list()
  Params[["rDxt"]]      = rDxt
  Params[["TxQualt"]]   = TxQualt
  Params[["InitPop"]]   = InitPop
  Params[["Mpfast"]]    = Mpfast
  Params[["ExogInf"]]   = ExogInf
  Params[["MpfastPI"]]  = MpfastPI
  Params[["Mrslow"]]    = Mrslow
  Params[["rrSlowFB"]]  = rrSlowFB
  Params[["rfast"]]     = rfast
  Params[["RRcurDef"]]  = RRcurDef
  Params[["rSlfCur"]]   = rSlfCur
  Params[["p_HR"]]      = p_HR
  Params[["dist_gen"]]  = matrix(c(1,rep(0,15)),4,4)
  Params[["vTMort"]]    = vTMort
  Params[["RRmuRF"]]    = RRmuRF
  Params[["RRmuHR"]]    = RRmuHR
  Params[["Birthst"]]   = Birthst
  Params[["HrEntEx"]]   = HrEntEx
  Params[["ImmNon"]]    = ImmNon
  Params[["ImmLat"]]    = ImmLat
  Params[["ImmAct"]]    = ImmAct
  Params[["ImmFst"]]    = ImmFst
  Params[["mubt"]]      = mubt
  Params[["RelInf"]]    = RelInf
  Params[["RelInfRg"]]  = RelInfRg
  Params[["Vmix"]]      = Vmix
  Params[["rEmmigFB"]]  = rEmmigFB
  Params[["TxVec"]]     = TxVec
  Params[["TunTxMort"]] = TunTxMort
  Params[["rDeft"]]     = rDeft
  Params[["pReTx"]]     = pReTx
  Params[["LtTxPar"]]   = LtTxPar
  Params[["LtDxPar"]]   = LtDxPar
  Params[["rLtScrt"]]   = rLtScrt
  Params[["RRdxAge"]]   = RRdxAge
  Params[["rRecov"]]    = rRecov
  Params[["pImmScen"]]  = pImmScen
  Params[["EarlyTrend"]]= EarlyTrend
  Params[["net_mig_usb"]]  = net_mig_usb
  Params[["net_mig_nusb"]]= net_mig_nusb
  # Params[["NixTrans"]]  = NixTrans
  # Params[["can_go"]]    = can_go
  # Params[["dist_goal"]] = dist_goal
  # Params[["diff_i_v"]]  = diff_i_v
  # Params[["dist_orig_v"]]=dist_orig_v
  Params[["adj_fact"]] <- adj_fact;

  Params[["ResNam"]]    = ResNam
  return(Params)
}

