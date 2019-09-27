################################################################################
##### THE CODE BELOW WILL GATHER AND FORMAT INPUTS FOR THE TB_MODEL        #####
##### FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME; #####
##### VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.        #####
################################################################################

#' This function takes the same inputs as the Outputs
#'@name fin2_param_init
#'@param PV vector of Inputs to format
#'@param loc two letter abbreviation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param prg_chng vector of program change values
#'@param ttt_list list of ttt changes
#'@return InputParams list
#'@export
fin2_param_init <- function(PV,loc,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0,prg_chng, ttt_list){
  InputParams <-vector("list", 51)   #Create an empty list to hold the formatted intitial parameters
  names(InputParams) <- c("rDxt","TxQualt", "InitPop", "Mpfast", "ExogInf", "MpfastPI",
                          "Mrslow", "rrSlowFB", "rfast"    ,"RRcurDef"      , "rSlfCur"  ,
                          "p_HR"        , "dist_gen" , "vTMort"   ,"RRmuRF"          , "RRmuHR",
                          "Birthst" , "HrEntEx"  ,"ImmNon"       , "ImmLat"  ,
                          "ImmAct"      , "ImmFst" , "mubt"     ,"RelInf"        , "RelInfRg" ,
                          "Vmix"       , "rEmmigFB" , "TxVec"    , "TunTxMort"    , "rDeft"    ,
                          "pReTx"      , "LtTxPar"  , "LtDxPar_lt"  ,"LtDxPar_nolt"  , "rLtScrt"      ,
                          "RRdxAge"  ,
                          "rRecov"      , "pImmScen"  ,   "EarlyTrend","net_mig_usb", "net_mig_nusb",
                          "aging_denom", "adj_fact","NixTrans"   ,"ResNam","ttt_sampling_dist",
                          "ttt_ag", "ttt_na", "ttt_month" ,"ttt_pop_frc", "ttt_ltbi" )
  ########## DEFINE A VARIABLE THAT WILL DETERMINE HOW LONG THE TIME DEPENDENT
  ########## VARIABLES SHOULD BE (IN MONTHS)
  month<-1201;
  prg_yr <-prg_chng["start_yr"]
  prg_m  <-(prg_yr-1950)*12
  InputParams[["ttt_month"]] <-seq((ttt_list[["StartYr"]]-1950)*12,(ttt_list[["EndYr"]]-1949)*12,1)
  ################################################################################
  ##### INTERVENTION
  ################################################################################
  if(Int5==1) {
    Int1 = Int2 = Int3 = Int4 = 1
  }

  #params from general reblnc_pop:
  InputParams[["dist_gen"]] <- dist_gen;

  # ################################################################################
  ###########################          INPUTS            #########################
  ################################################################################
  ################################################################################
  BgMort              <- as.matrix(Inputs[["BgMort"]])
  if(loc=="US"){
    BgMort[1:68,2:12] <-weight_mort(loc)
  } else{
    BgMort[10:67,2:12]<-weight_mort(loc)
  }
  ###SSA projections for reductions in mortality going forward
  ### these should be calculated and stored as new input data
  for(j in 68:151){
    for (i in 1:2){
      BgMort[j,i]<-BgMort[j-1,i]*(1-.0159)
    }
    for (i in 3:5){
      BgMort[j,i]<-BgMort[j-1,i]*(1-.0090)
    }
    for (i in 6:7){
      BgMort[j,i]<-BgMort[j-1,i]*(1-.0107)
    }
    for (i in 8:9){
      BgMort[j,i]<-BgMort[j-1,i]*(1-.0083)
    }
    for (i in 10:11){
      BgMort[j,i]<-BgMort[j-1,i]*(1-.0069)
    }
  }
  InputParams[["InitPop"]] <- Inputs[["InitPop"]]
  Births           <- Inputs[["Births"]]
  ImmigInputs      <- Inputs[["ImmigInputs"]]
  ImmigInputs$PrevTrend25_34<-crude_rate(Inputs,loc)
  TxInputs         <- Inputs[["TxInputs"]]
  NetMig           <- Inputs[["NetMigrState"]]

  ##########                CALCULATION OF AGING DENOMINATORS           ##########
  InputParams[["aging_denom"]]<-age_denom(loc)

  ##########                PARAMETER DEFINITIONS                      ###########
  ##########                RISK FACTOR DISTRIBUTIONS   ##########################
  InputParams[["adj_fact"]] <- exp(PV[["adj_ag1"]]*(10:0)/11 + PV[["adj_ag11"]]*(0:10)/11)

  #######################           BIRTHS                 #######################
  ####### INDEXED BY TIME, ABSOLUTE NUMBER OF NEW ADULT ENTRANTS OVER TIME #######

  Birthst   <- SmoCurve(Births)*PV["TunBirths"]/12
  Birthst   <- Birthst[1:month]
  InputParams[["Birthst"]]<-Birthst
  ##########################      MORTALITY RATES       ##########################
  ########################## BACKGROUND MORTALITY BY TIME ########################
  mubt      <- matrix(NA,1801,11)
  for(i in 1:11){
    mubt[,i] <- SmoCurve(BgMort[,i+1])*PV[["TunMubt"]] /12
  }
  mubt<-mubt[1:month,]
  for(i in 1:11){
    mubt[,i] <- mubt[,i]*exp(PV[["TunmuAg"]]*i)
  }
  InputParams[["mubt"]] <-mubt

  #########################     DISEASE SPECIFIC       ###########################
  #############    ACTIVE TB RATES DEFAULT TO THE SMEAR POS LEVELS   #############
  muIp  	  <- PV["muIp"]/12
  ######################## MULTIPLER OF MORT RATE ABOVE ########################

  TunmuTbAg <- PV["TunmuTbAg"]

  ###############  RATE RATIO OF MORTALITY INCREASE FOR HIGH RISK ###############
  RRmuHR<-c(1,PV["RRmuHR"])
  InputParams[["RRmuHR"]]    <- RRmuHR

  ############### CREATE A MATRIX OF RF MORTALITIES BY rf GROUP ###############
  mort_dist<-rowSums(InputParams[["dist_gen"]])

  RF_fact=20
  #
  InputParams[["RRmuRF"]]   <- rep(NA,4);
  names(InputParams[["RRmuRF"]]) <- c("RF1","RF2","RF3","RF4")

  InputParams[["RRmuRF"]]<-exp((0:3)/3*log(RF_fact))
  InputParams[["RRmuRF"]]<-InputParams[["RRmuRF"]]/sum(InputParams[["RRmuRF"]]*mort_dist)

  ############### CREATE A MATRIX OF TB MORTALITIES BY AGE GROUP ###############
  TunmuTbAg <- PV["TunmuTbAg"]

  vTMort   <- matrix(0,11,6);
  rownames(vTMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
  colnames(vTMort) <- c("Su","Sp","Ls","Lf","Ac","Tx")
  vTMort[,5:6] <- muIp #active disease rates default to smear positive
  RRmuTbAg <- exp(c(0,0,1:9)*TunmuTbAg)
  for(i in 1:ncol(vTMort)) {
    vTMort[,i] <- vTMort[,i] * RRmuTbAg
  }
  InputParams[["vTMort"]]  <- vTMort

  ######################## MULTIPLER OF MORT RATE ABOVE ########################

  ######################         IMMIGRATION             ########################
  ######################         OVERALL IMM.            ########################
  TotImmig0       <- (c(Inputs$ImmigInputs[[1]][1:151])+c(rep(0,67),cumsum(rep(PV["ImmigVolFut"],84))))/12*PV["ImmigVol"]
  TotImmAge0      <-matrix(0,151,11)
  for (i in 1:151){
    for (j in 1:11){
      TotImmAge0[i,j]   <- TotImmig0[i]*as.matrix(ImmigInputs$AgeDist[j,i])
    } }
  TotImmAge <-matrix(0,1801,11)
  # for (i in 1:1801){
  for (j in 1:11){
    TotImmAge[,j]        <- SmoCurve(TotImmAge0[,j])
  }
  # }
  ######################           LTBI IMM.             ########################
  PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:67]^(exp(PV["TunLtbiTrend"]))*ImmigInputs[["PrevTrend25_34"]][67]^(1-exp(PV["TunLtbiTrend"])),
                        ImmigInputs[["PrevTrend25_34"]][68:151]*(PV["ImmigPrevFutLat"]/0.99)^(1:84))
  PrevTrend25_341l <-   PrevTrend25_340l
  PrevTrend25_34l  <- SmoCurve(PrevTrend25_341l)
  PrevTrend25_34_ls <- (PrevTrend25_34l);
  PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
  ImmLat          <- matrix(NA,length(PrevTrend25_34_ls),11)
  for(i in 1:11) ImmLat[,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]*PV["LtbiPar1"]-(c(2.5,1:9*10,100)/100)[i]^2*PV["LtbiPar2"])*PrevTrend25_34_ls))*TotImmAge[,i]

  ######################         ACTIVE TB IMM.           ########################
  # PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:69]^exp(PV["TunActTrend"])*ImmigInputs[["PrevTrend25_34"]][69]^exp((1-PV["TunActTrend"])),

  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:67]^(exp(PV["TunActTrend"]))*ImmigInputs[["PrevTrend25_34"]][67]^(1-exp(PV["TunActTrend"])),
                        ImmigInputs[["PrevTrend25_34"]][68:151]*(PV["ImmigPrevFutAct"]/0.99)^(1:84))

  PrevTrend25_34a  <- SmoCurve(PrevTrend25_340a)
  act_prob<-rep(0,1801)
  for (t in 1:1801){
    act_prob[t]<-((PV[["pImActSlp"]]*t)/1801)+PV[["pImActIntc"]]
  }

  ImmAct         <- outer(PrevTrend25_34a*PV["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*act_prob
  ImmFst         <- outer(PrevTrend25_34a*PV["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-act_prob)

  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat
  ###################### TRUNCATE THESE VALS
  InputParams[["ImmAct"]]<-ImmAct[1:month,];InputParams[["ImmFst"]]<-ImmFst[1:month,]; InputParams[["ImmLat"]]<-ImmLat[1:month,]; InputParams[["ImmNon"]]<-ImmNon[1:month,]
  # ImmNon[1:6,]<- ImmAct[1:6,]<-ImmFst[1:6,]<-ImmLat[1:6,]<-0

  #### #### #### SCEN 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Scen2==1) {
    for(i in 1:11) InputParams[["ImmLat"]][,i] <- InputParams[["ImmLat"]][,i]*(1-LgtCurve(2018,2019,1))
    for(i in 1:11) InputParams[["ImmAct"]][,i] <- InputParams[["ImmAct"]][,i]*(1-LgtCurve(2018,2019,1))
    for(i in 1:11) InputParams[["ImmFst"]][,i] <- InputParams[["ImmFst"]][,i]*(1-LgtCurve(2018,2019,1))
    InputParams[["ImmNon"]]         <- TotImmAge[1:1201,]-InputParams[["ImmAct"]]-InputParams[["ImmFst"]]-InputParams[["ImmLat"]]

  }


  #### #### #### SCEN 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Scen3==1) {
    adjf <- (1-0.5)^(1/10/12);  adjfV <- adjf^(0:384)
    for(i in 1:11) InputParams[["ImmLat"]][817:1201,i] <- InputParams[["ImmLat"]][817:1201,i]*adjfV
    for(i in 1:11) InputParams[["ImmAct"]][817:1201,i] <- InputParams[["ImmAct"]][817:1201,i]*adjfV
    for(i in 1:11) InputParams[["ImmFst"]][817:1201,i] <- InputParams[["ImmFst"]][817:1201,i]*adjfV
    InputParams[["ImmNon"]]         <- TotImmAge[1:1201,]-InputParams[["ImmAct"]]-InputParams[["ImmFst"]]-InputParams[["ImmLat"]]
  }

  #######################__EXOGENEOUS INFECTION RISK    #######################
  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  ExogInf        <- PV["ExogInf"]*PrevTrend25_34a/PrevTrend25_340a["2013"]/12
  InputParams[["ExogInf"]]        <- ExogInf[1:month]
  #removed *(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])
  ###############################    EMMIGRATION   ##############################

  InputParams[["rEmmigFB"]] <- c(PV["rEmmigF1"],PV["rEmmigF2"])/12

  ######################          NET MIGRATION          #########################

  InputParams[["net_mig_usb"]]  <- (NetMig[,"usb" ]*PV["TunNetMig"])^(1/12)-1
  InputParams[["net_mig_nusb"]] <- (NetMig[,"nusb"]*PV["TunNetMig"])^(1/12)-1

  ##########################  HIGH-RISK ENTRY/EXIT  #############################

  InputParams[["p_HR"]]    <- PV["pHR"]
  yr       <- c(2.5,1:10*10)
  r0_5     <- 1/3;
  r45_55 <- 1/20
  HR_exit  <- r0_5*((r45_55/r0_5)^(1/(50-2.5)))^(yr-2.5)
  HR_entry <- HR_exit*InputParams[["p_HR"]]*1.3
  InputParams[["HrEntEx"]]  <- cbind(HR_entry,HR_exit)/12

  ##########################   TB TRANSMISSION  #################################

  CR           <- PV["CR"]/12
  TrIn         <- PV["TrIn"]	# Contact rate for In as a fraction of Ip
  InputParams[["RelInfRg"]]    <- c(1.0,PV["RelCrHr"], 1.0, PV["RelCrHr"])*CR
  TunTbTransTx <- 0 #PV["TunTbTransTx"]  # set to zero?
  InputParams[["Vmix"]]         <- 1-c(PV["sigmaHr"],PV["sigmaFb"])
  RelInf <- rep(0,6)
  names(RelInf) <-c("Su","Sp","Ls","Lf","Ac","Tx")

  RelInf[5] <- 1;
  RelInf[6] <- RelInf[5]*TunTbTransTx
  InputParams[["RelInf"]] <-RelInf

  #########################   TB NATURAL HISTORY  ###############################
  #########################__   EARLY EPIDEMTIC   ###############################
  Early0 <- PV["Early0"]
  InputParams[["EarlyTrend"]] <- c(rep(1+Early0,200*12),seq(1+Early0,1.0,length.out=50*12+2))

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

  ############ UPDATE PROBS FOR LEVEL 1 OF REACTIVATION ###########
  Mpfast[,1]   <- vORpfastRF[1]*Mpfast[,1]
  MpfastPI[,1]   <- vORpfastPIRF[1]*MpfastPI[,1]
  ############ UPDATE PROBS FOR LEVEL 2 OF REACTIVATION ###########
  Mpfast[,2]   <- vORpfastRF[2]*Mpfast[,2]
  MpfastPI[,2]   <- vORpfastPIRF[2]*MpfastPI[,2]
  ############ UPDATE PROBS FOR LEVEL 3 OF REACTIVATION ###########
  Mpfast[,3]   <- vORpfastRF[3]*Mpfast[,3]
  MpfastPI[,3]   <- vORpfastPIRF[3]*MpfastPI[,3]
  ############ UPDATE PROBS FOR LEVEL 4 OF REACTIVATION ###########
  Mpfast[,4]   <- vORpfastRF[4]*Mpfast[,4]
  MpfastPI[,4]   <- vORpfastPIRF[4]*MpfastPI[,4]

  ##### UPDATE BOTH MATRICES WITH PROBABILITIES, NOT RATES
  Mpfast[,]    <- Mpfast[,]  /(1+Mpfast[,]);
  MpfastPI[,]  <- MpfastPI[,]/(1+MpfastPI[,]);
  InputParams[["Mpfast"]]    <- Mpfast
  InputParams[["MpfastPI"]]  <- MpfastPI


  rslow      <- PV["rslow"]/12
  rslowRF    <- PV["rslowH"]/12
  RRrslowRF  <- rslowRF/rslow
  InputParams[["rfast"]]    <- PV["rfast"]/12
  # rrSlowFB0  <- PV["rrSlowFB"]
  InputParams[["rrSlowFB"]]   <- c(1,1,1)

  ############# CREATE A VECTOR FOR RATE OF SLOW PROGRESSION THAT WILL
  ############# VARY BASED ON LEVELS OF TB REACTIVATION RATES
  Vrslow     <- rep(NA,4)
  ############# UPDATE LEVEL FOUR OF THE RATE OF SLOW BASED ON CALCULATED RR FROM
  ############# USER INPUTTED RR FOR THE RISK FACTOR
  Vrslow=rslow*exp((0:3)/3*log(RRrslowRF))

  TunrslowAge  <- PV["TunrslowAge"]
  rrReactAg       <- exp(c(0,0,0,0,0,0,0.5,1:4)*PV["TunrslowAge"])
  InputParams[["Mrslow"]] <- outer(rrReactAg,Vrslow)
  # InputParams[["Mrslow"]] <- InputParams["Mrslow"];


  #######################       RATE OF RECOVERY          ########################

  InputParams[["rRecov"]]     <-  PV["rRecov"]/12

  #######################       RATE OF SELF CURE         ########################

  InputParams[["rSlfCur"]]      <- PV["rSlfCur"]/12

  ######################          LTBI DIAGNOSIS           ########################

  ######################        TEST SPECIFICATIONS          ######################
  Sens_IGRA <-c(.780,.675,.712,.789,.591)
  Spec_IGRA <-c(.979,.958,.989,.985,.931)
  IGRA_frc<-.33
  Sens_TST <-c(.726,.540,.691,.807,.570)
  Spec_TST <-c(.921,.965,.739,.70,.885)
  names(Sens_TST)<- names(Spec_TST)<-names(Sens_IGRA)<- names(Spec_IGRA)<-c("US","hivUS","youngNUS","NUS","hivNUS")

  ###calculate the weighted mean of sensitivity and specificity
  SensLt<-matrix(NA,length(Sens_IGRA),month)
  SpecLt<-matrix(NA,length(Spec_IGRA),month)
  rownames(SensLt)<-rownames(SpecLt)<-names(Sens_IGRA)
  ###sensitivity and specificity must be time varying parameters in order to allow the user to change the
  ###percentage of IGRA used at a specific time range
  for (i in 1:nrow(SensLt)){
    SensLt[i,]          <-rep((Sens_IGRA[i]*IGRA_frc + (1-IGRA_frc)*Sens_TST[i]),month)
    SpecLt[i,]          <-rep((Spec_IGRA[i]*IGRA_frc + (1-IGRA_frc)*Spec_TST[i]),month)
  }

  ######################     IGRA FRACTION PROGRAM CHANGE    ########################
  if (prg_chng["IGRA_frc"] != IGRA_frc){
    for (i in 1:nrow(SensLt)){
      SensLt[i,prg_m:ncol(SensLt)]    <-rep((Sens_IGRA[i]*prg_chng["IGRA_frc"] + (1-prg_chng["IGRA_frc"])*Sens_TST[i]),1+1201-prg_m)
      SpecLt[i,prg_m:ncol(SpecLt)]    <-rep((Spec_IGRA[i]*prg_chng["IGRA_frc"] + (1-prg_chng["IGRA_frc"])*Spec_TST[i]),1+1201-prg_m)
    }}

  ### ADJUST THIS FOR THE FOREIGN BORN
  ######################        SCREENING RATES          ######################
  rLtScrt       <- LgtCurve(1985,2015,PV["rLtScr"])/12
  ######################  SCREENING RATE PROGRAM CHANGE ########################
  if (prg_chng["scrn_cov"] !=1) {
    rLtScrt[prg_m:length(rLtScrt)]<-rLtScrt[prg_m:length(rLtScrt)]*prg_chng["scrn_cov"];
  }
  ################################################################################
  #################### TTT ADDITIONAL SCREENING PROBABILITIES ####################
  ################################################################################
  ###this is dependent on a basecase run so load in that data
  if (ttt_list[[3]]!=0 & ttt_list[[4]]!=0){
    load(system.file("US/US_results_1.rda", package="MITUS"))
    x<-create_ttt_dist(ttt_list = ttt_list,
                       results = out[1,,],
                       PV = PV)
    ttt_sampling_dist<-x[[1]]/12
    ttt_pop_frc<-x[[2]]
    ttt_ag<-switch(ttt_list[["AgeGrp"]], "All"=0:10,
                   "0 to 24"=0:2,
                   "25 to 64"=3:6,
                   "65+"=7:10
    )
    ttt_na<-switch(ttt_list[["NativityGrp"]], "All"=0:2,
                   "USB"=0,
                   "NUSB"=1:2
    )
  } else{
    ttt_sampling_dist<-matrix(0,4,4)
    ttt_na<-99
    ttt_ag<-99
    ttt_pop_frc<-0
  }
  InputParams[["ttt_sampling_dist"]]<-ttt_sampling_dist
  InputParams[["ttt_na"]]<-ttt_na
  InputParams[["ttt_ag"]]<-ttt_ag
  InputParams[["ttt_pop_frc"]]<-ttt_pop_frc
  InputParams[["ttt_ltbi"]]<-ttt_list[["RRPrev"]]

  ###adjustments to the screening rates dependent on risk and TB status
  rrTestHr      <- PV["rrTestHr"] # RR of LTBI screening for HIV and HR as cmpared to general
  rrTestLrNoTb  <- PV["rrTestLrNoTb"] # RR of LTBI screening for individuals with no risk factors
  ### because of the introduction of new time varying parameters, we will create 2 matrices to
  ### hold the three different sensitivity and specificity measures; one will be for those whose
  ### true LTBI status is positive and the other is for those whose true TB status is negative.
  LtDxPar_nolt <- LtDxPar_lt <- matrix(NA,nrow(SensLt),month);
  rownames(LtDxPar_lt) <- rownames(LtDxPar_nolt) <- rownames(SensLt)

  ##adjust for no latent
  LtDxPar_lt   <-SensLt
  LtDxPar_nolt <- 1-SpecLt
  #adjust for High Risk Populations
  LtDxPar_lt[2,]     <-rrTestHr*LtDxPar_lt[1,]
  LtDxPar_nolt[2,]   <-rrTestHr*LtDxPar_nolt[1,]
  #High risk foriegn born
  LtDxPar_lt[5,]     <-rrTestHr*LtDxPar_lt[4,]
  LtDxPar_nolt[5,]   <-rrTestHr*LtDxPar_nolt[4,]
  #adjust for no latent
  LtDxPar_nolt[1,]   <-LtDxPar_nolt[1,]*rrTestLrNoTb

  InputParams[["LtDxPar_lt"]]<-LtDxPar_lt;   InputParams[["LtDxPar_nolt"]]<-LtDxPar_nolt;

  ######################          LTBI DIAGNOSIS           ########################
  ###################### LTBI TX EFFECTIVENESS PROGRAM CHANGE ########################
  if (prg_chng["ltbi_eff_frc"] != round(PV["EffLt"], 2)){
    EffLt         <- c(rep(PV["EffLt"],prg_m-1),rep(prg_chng["ltbi_eff_frc"],1+1201-prg_m))
  } else {
    EffLt         <- rep(PV["EffLt"],month)
  }

  ### PROBABILITY OF LATENT TREATMENT INTIATION
  pTlInt        <- rep(.775,month)
  ###################### LTBI TX INITIATION PROGRAM CHANGE ########################
  if (prg_chng["ltbi_init_frc"] !=pTlInt[prg_m]){
    pTlInt[prg_m:length(pTlInt)] <- prg_chng["ltbi_init_frc"];
  }
  ######################    LTBI TREATMENT DEFAULT          ######################
  pDefLt<-rep(PV["pDefLt"],month)   ### PROBABILITY OF LATENT TREATMENT DEFAULT
  ###################### LTBI TREATMENT COMPLETION PROGRAM CHANGE ######################
  if (prg_chng["ltbi_comp_frc"] != 1-(round(pDefLt[prg_m], 2))){
    pDefLt[prg_m:length(pDefLt)]<-(1-prg_chng["ltbi_comp_frc"])
  }
  ### because of the introduction of new time varying parameters, we will create a matrix to
  ### hold the latent treatment parameters
  LtTxPar       <- matrix(NA,3,month)
  InputParams[["LtTxPar"]]       <- cbind(pTlInt,pDefLt,EffLt)

  pImmScen    <- PV["pImmScen"] # lack of reactivitiy to IGRA for Sp

  #### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  ### HOW TO ADD PROGRAM CHANGE HERE?
  if(Int2==1) { rLtScrt     <- rLtScrt  + LgtCurve(2018,2023,1)*rLtScrt*1}
  InputParams[["rLtScrt"]]       <- rLtScrt

  InputParams[["pImmScen"]]   <- PV["pImmScen"] # lack of reactivitiy to IGRA for Sp

  ################################################################################
  #######################         TB DIAGNOSIS            ########################
  #######################      TEST CHARACTERISTICS       ########################

  SensSp    <- PV["SensSp"]

  ######################           PROVIDER DELAY         ########################

  DelaySp    <- rep(PV["DelaySp"],month)

  if (prg_chng["tb_tim2tx_frc"] !=100){
    DelaySp[prg_m:length(DelaySp)] <- PV["DelaySp"]*(prg_chng["tb_tim2tx_frc"])/100;
  }

  #######################     PROVIDER DELAY RR ELDERLY     ########################

  TunRRdxAge    <- PV["TunRRdxAge"]
  InputParams[["RRdxAge"]]       <- 1+(c(rep(1,6),cumprod(seq(1.05,1.3,length.out=5)))-1)*TunRRdxAge

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
  rDxt1<-rDxt1[1:1201,]

  # Put it all together
  rDxt           <- 1/(1/rDxt1+DelaySp)*SensSp
  rDxt[,2]       <- (rDxt[,1]-min(rDxt[,1]))/PV["rrDxH"]+min(rDxt[,1]) #check this with Nick
  colnames(rDxt) <- c("Active","Active_HighRisk")


  #### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Int3==1) { for(i in 1:2) { rDxt[,i] <- rDxt[,i]+ rDxt[,i]*LgtCurve(2018,2023,1)     }   }
  InputParams[["rDxt"]]<-rDxt

  #### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  ######  TREATMENT OUTCOMES  ######
  InputParams[["TunTxMort"]]	<- PV["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

  # Regimen duration (no uncertainty?)
  d1st <- 1/9

  ## Regimen efficacy
  pCurPs  <- PV["pCurPs"]    # probability of cure with pansensitive TB, 1st line regimen (Menzies 09)
  # TxE1		<- PV["TxEf1"]     # RR cure given mono-resistance, 1st line regimen 0.90
  # TxE2		<- PV["TxEf2"]     # RR cure given multiresistance, 1st line or 2nd line 0.50
  # TxE3		<- PV["TxEf3"]     # RR cure given effective 2nd line 0.90

  ###########################        REGIMEN DEFAULT       ##########################

  rDef0         <- rep(NA,151)
  rDef0[1:30]   <- PV["TxDefEarly"]
  rDef0[44:63]  <- ORAdd(TxInputs[[1]][,2],PV["TunTxDef"])
  rDef0[64:151] <- rDef0[63]
  rDef1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=rDef0[-(31:43)],spar=0.4),x=1950:2100)$y
  rDeft         <- SmoCurve(rDef1)/12;
  rDeft<-rDeft[1:month]


  if (round(prg_chng["tb_txdef_frc"],2) !=round(rDeft[prg_m], 2)){
    rDeft[prg_m:length(rDeft)] <- prg_chng["tb_txdef_frc"];
  }

  rDeftH        <- rDeft*PV["RRdefHR"] # second col is HR default rate
  rDeftH<-rDeftH[1:month]
  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Int4==1) {  rDeftH <- rDeftH*(1-LgtCurve(2018,2023,0.5))      }
  if(Int4==1) {  rDeft<-rDeft[1:1201]
  rDeft<- rDeft * (1-LgtCurve(2018,2023,0.5))      }

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  InputParams[["rDeft"]] <-rDeft
  ## Tx quality
  TxQual0         <- rep(NA,151)
  TxQual0[1:30]   <- PV["TxQualEarly"]
  TxQual0[44:62]  <- ORAdd(TxInputs[[2]][,2],PV["TunTxQual"])
  TxQual0[63:151] <- TxQual0[62]
  TxQual1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=TxQual0[-(31:43)],spar=0.4),x=1950:2100)$y
  InputParams[["TxQualt"]]         <- SmoCurve(TxQual1);
  InputParams[["TxQualt"]]         <-InputParams[["TxQualt"]][1:1201]
  InputParams[["RRcurDef"]]        <- PV["RRcurDef"]

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Int4==1) {  InputParams[["TxQualt"]] <- 1-(1-InputParams[["TxQualt"]])*(1-LgtCurve(2018,2023,0.5))      }

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  #### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  InputParams[["NixTrans"]] <- 1-LgtCurve(2018,2019,1)
  if(Scen1==0) {  InputParams[["NixTrans"]][] <- 1      }

  #### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  ## Retreatment
  InputParams[["pReTx"]]   <- LgtCurve(1985,2000,PV["pReTx"])   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)

  ### REMOVED ACQUIRED RESISTANCE PARAMETERS
  ### REPLACED TX MAT WITH VECTOR FOR SINGLE VALUES
  #####################   NEW TB TREATMENT TABLE ################
  InputParams[["TxVec"]]           <- rep(NA,2)
  names(InputParams[["TxVec"]]) <- c("TxCompRate","TxEff")
  InputParams[["TxVec"]][1]       <-  d1st
  InputParams[["TxVec"]][2]       <-  pCurPs
  ### REMOVED ART PARAMETERS


  InputParams[["ResNam"]]<-func_ResNam()
  return(InputParams)

  ###################################################
}
