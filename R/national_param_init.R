################################################################################
##### THE CODE BELOW WILL GATHER AND FORMAT INPUTS FOR THE TB_MODEL        #####
##### FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME; #####
##### VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.        #####
################################################################################

#' This function takes the same inputs as the Outputs
#'@name  national_param_init
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
#'@param Scen4 boolean for scenario 4 (additional for paper analysis)
#'@param Scen5 boolean for scenario 5 (additional for paper analysis)
#'@param Scen6 boolean for scenario 6 (additional for paper analysis)
#'@param prg_chng vector of program change values
#'@param ttt_list list of ttt changes
#'@return InputParams list
#'@export
national_param_init <- function(PV,loc,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,
                                Scen1=0,Scen2=0,Scen3=0,Scen4=0,Scen5=0, Scen6=0,
                                prg_chng, ttt_list, delay=0, immig=1, noiseParams=list(1,1,1)){
  ################################################################################
  ##### DEFINE VARIABLES THAT WILL DETERMINE HOW LONG THE TIME               #####
  ##### DEPENDENT VARIABLES SHOULD BE                                        #####
  ################################################################################
  month<-1801;
  intv_yr<-2022
  intv_m<-((intv_yr-1949)*12)+1
  prg_yr <-prg_chng["start_yr"]
  prg_m  <-((prg_yr-1949)*12)+1
  ttt_month <-seq(((ttt_list[[1]][["StartYr"]]-1950)*12) + 6,((ttt_list[[1]][["EndYr"]]-1950)*12) + 6,12)
  ################################################################################
  ##### CHECK IF ALL INTERVENTIONS OPTION IS SELECTED AND UPDATE INT VARS ########
  ################################################################################
  if(Int5==1) {
    Int1 = Int2 = Int3 = Int4 = 1
  }
  ################################################################################
  ###########################          INPUTS            #########################
  ################################################################################
  BgMort              <- as.matrix(Inputs[["BgMort"]])
  if(loc=="US"){
    BgMort[1:68,2:12] <-weight_mort(loc)
  } else{
    BgMort[10:67,2:12]<-weight_mort(loc)
  }
  ################################################################################
  ###SSA projections for reductions in mortality going forward
  ### these should be calculated and stored as new input data
  ################################################################################
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
  ################################################################################
  InitPop          <- Inputs[["InitPop"]]
  Births           <- Inputs[["Births"]]
  ImmigInputs      <- Inputs[["ImmigInputs"]]
  ################################################################################
  ##################### CHECK FOR SCENARIO 6             ########################
  ##################### CHANGES IMMIGRATION BURDEN       ########################
  if(Scen6==1){
    r_decline=.03
  } else if (Scen6==2){
    r_decline=0
  } else {
    r_decline=0.015
  }
  ImmigInputs$PrevTrend25_34<-crude_rate(Inputs,loc,r_decline)
  ImmigInputs$PrevTrend25_34[74:101] <- ImmigInputs$PrevTrend25_34[73]*noiseParams[[2]]
  ################################################################################
  TxInputs         <- Inputs[["TxInputs"]]
  NetMig           <- Inputs[["NetMigrState"]]
  ################################################################################
  ##########                CALCULATION OF AGING DENOMINATORS           ##########
  spl_den <-age_denom(loc, month)
  ################################################################################
  ##########                PARAMETER DEFINITIONS                      ###########
  ##########                RISK FACTOR DISTRIBUTIONS   ##########################
  #turned off because the rebalancing has been moved back into the model
  adj_fact <- exp(PV[["adj_ag1"]]*(10:0)/11 + PV[["adj_ag11"]]*(0:10)/11)

  ################################################################################
  #######################           BIRTHS                 #######################
  ####### INDEXED BY TIME, ABSOLUTE NUMBER OF NEW ADULT ENTRANTS OVER TIME #######

  Birthst   <- SmoCurve(Births)*PV["TunBirths"]/12
  Birthst   <- Birthst[1:month]

  ################################################################################
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

  #########################     DISEASE SPECIFIC       ###########################
  #############    ACTIVE TB RATES DEFAULT TO THE SMEAR POS LEVELS   #############

  muIp  	  <- PV["muIp"]/12

  ######################## MULTIPLER OF MORT RATE ABOVE ########################

  TunmuTbAg <- PV["TunmuTbAg"]

  ###############  RATE RATIO OF MORTALITY INCREASE FOR HIGH RISK ###############

  RRmuHR    <- c(1,PV["RRmuHR"])

  ############### CREATE A MATRIX OF RF MORTALITIES BY AGE GROUP ###############
  mort_dist<-rowSums(dist_gen)

  RF_fact=20

  RRmuRF    <- rep(NA,4);
  names(RRmuRF) <- c("RF1","RF2","RF3","RF4")

  RRmuRF<-exp((0:3)/3*log(RF_fact))
  RRmuRF<-RRmuRF/sum(RRmuRF*mort_dist)
  #check
  #RRmuRF%*%mort_dist
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

  ##################### CHECK FOR SCENARIO 5             ########################
  ##################### CHANGES IMMIGRATION VOLUME       ########################
  if(Scen5 == 0){
    Inputs$ImmigInputs[["TotByYear"]] <- Inputs$ImmigInputs[["TotByYear"]]
  } else if(Scen5 == 1){
    low_immig<-base_immig <- Inputs$ImmigInputs[["TotByYear"]]
    low_immig[69:151] <- exp(log(base_immig[69:151]) - (log(base_immig[69:151]*1.5) - log(base_immig[69:151])))
    Inputs$ImmigInputs[["TotByYear"]] <- low_immig
    #high immigration
  } else if (Scen5 == 2){
    high_immig <- Inputs$ImmigInputs[["TotByYear"]]
    high_immig[69:151] <- Inputs$ImmigInputs[["TotByYear"]][69:151] * 1.5
    Inputs$ImmigInputs[["TotByYear"]] <- high_immig
  }
  Inputs$ImmigInputs[["TotByYear"]][74:101]<-Inputs$ImmigInputs[["TotByYear"]][73]*noiseParams[[1]]
  ####################### calculate the age immigration #######################
  TotImmig0       <- (c(Inputs$ImmigInputs[[1]][1:151])+c(rep(0,71),cumsum(rep(PV["ImmigVolFut"],80))))/12*PV["ImmigVol"]
  # TotImmig0[74:151]       <- TotImmig0[73]*noiseParams[[1]]
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
  if(immig != 1){
    TotImmAge[843:855,]<-TotImmAge[843:855,]*immig;
  }
  # }
  ######################           LTBI IMM.             ########################
  PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:71]^(exp(PV["TunLtbiTrend"]))*ImmigInputs[["PrevTrend25_34"]][71]^(1-exp(PV["TunLtbiTrend"])),
                        ImmigInputs[["PrevTrend25_34"]][72:151]*(PV["ImmigPrevFutLat"]/0.99)^(1:80))
  PrevTrend25_341l <-   PrevTrend25_340l
  PrevTrend25_34_ls  <- SmoCurve(PrevTrend25_341l)
  PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
  # plot(PrevTrend25_34_ls[ (65*12):(75*12)], type="l")
  ImmLat          <- matrix(NA,length(PrevTrend25_34_ls),11)

  for(i in 1:11) ImmLat[,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]^PV["LtbiPar1"])*PrevTrend25_34_ls*PV["LtbiPar2"]))*TotImmAge[,i]
  ######################         ACTIVE TB IMM.           ########################
  PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:71]^(exp(PV["TunActTrend"]))*ImmigInputs[["PrevTrend25_34"]][71]^(1-exp(PV["TunActTrend"])),
                        ImmigInputs[["PrevTrend25_34"]][72:151]*(PV["ImmigPrevFutAct"]/0.99)^(1:80))

  PrevTrend25_34a  <- SmoCurve(PrevTrend25_340a)

  ImmAct         <- outer(PrevTrend25_34a*exp(PV["RRtbprev"]),ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*PV["pImmAct"]
  ImmFst         <- outer(PrevTrend25_34a*exp(PV["RRtbprev"]),ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-PV["pImmAct"])
  ImmNon         <- TotImmAge-ImmAct-ImmFst-ImmLat
  for (i in 1:length(ImmNon)) if (ImmNon[i]<0) ImmNon[i]<-0

  ###################### TRUNCATE THESE VALS
  ImmAct<-ImmAct[1:month,];ImmFst<-ImmFst[1:month,]; ImmLat<-ImmLat[1:month,]; ImmNon<-ImmNon[1:month,]
  SpImmNon       <- matrix(0,month, 11)

  #### #### #### SCEN 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  if(Scen2==1) {
    baseline_ImmNon<-ImmNon
    for(i in 1:11) ImmLat[,i] <- ImmLat[,i]*(1-LgtCurve(intv_yr,intv_yr+1,1))
    for(i in 1:11) ImmAct[,i] <- ImmAct[,i]*(1-LgtCurve(intv_yr,intv_yr+1,1))
    for(i in 1:11) ImmFst[,i] <- ImmFst[,i]*(1-LgtCurve(intv_yr,intv_yr+1,1))
    ImmNon        <- TotImmAge[1:month,]-ImmAct-ImmFst-ImmLat
    SpImmNon<-baseline_ImmNon-ImmNon
  }
#no immigration immediately at 2020
  if(Scen2==2) {
    baseline_ImmNon<-ImmNon
    vec.change<-c(rep(1,intv_m), rep(0,month-intv_m))
    for(i in 1:11) ImmLat[,i] <- ImmLat[,i]*vec.change
    for(i in 1:11) ImmAct[,i] <- ImmAct[,i]*vec.change
    for(i in 1:11) ImmFst[,i] <- ImmFst[,i]*vec.change
    ImmNon        <- TotImmAge[1:month,]-ImmAct-ImmFst-ImmLat
    SpImmNon<-baseline_ImmNon-ImmNon
  }

  #### #### #### SCEN 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Scen3==1) {
    adjf <- (1-0.5)^(1/10/12);  adjfV <- adjf^(0:(length(intv_m:month)-1))
    for(i in 1:11) ImmLat[intv_m:month,i] <- ImmLat[intv_m:month,i]*adjfV
    for(i in 1:11) ImmAct[intv_m:month,i] <- ImmAct[intv_m:month,i]*adjfV
    for(i in 1:11) ImmFst[intv_m:month,i] <- ImmFst[intv_m:month,i]*adjfV
    ImmNon         <- TotImmAge[1:month,]-ImmAct-ImmFst-ImmLat
  }

  ######################   EXOGENEOUS INFECTION RISK      ########################

  ExogInf        <- matrix(NA,length(PrevTrend25_34a),5)
  ExogInf        <- PV["ExogInf"]*PrevTrend25_34a/PrevTrend25_340a["2020"]/12
  ExogInf        <- ExogInf[1:month]
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
  RRcrAG       <- c(0.517233655,
                    0.77454377,
                    1,
                    0.765666833,
                    0.783155643,
                    0.754892753,
                    0.581666584,
                    0.321227236,
                    0.244200692,
                    0.244200692,
                    0.244200692)
  RelInfRg     <- c(1.0,PV["RelCrHr"],1.0,PV["RelCrHr"])*CR
  TunTbTransTx <- 0 #PV["TunTbTransTx"]  # set to zero?
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
  ORpfastRF  <- 40 ##riskfactor
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
  vORpfastRF  <-vORpfastRF/sum(vORpfastRF*mort_dist)

  vORpfastPIRF  <- vORpfastRF*ORpfastPI

  ############ UPDATE PROBS FOR LEVEL 2 OF REACTIVATION ###########
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

  rslow      <- PV["rslow"]/12
  rrslowRF    <- 40 #PV["rslowH"]/12

  RRrslowRF  <- exp((0:3)/3*log(rrslowRF))
  RRrslowRF<-RRrslowRF/sum(RRrslowRF*mort_dist)
  rfast      <- PV["rfast"]/12
  #rrSlowFB0  <- PV["rrSlowFB"] #removed
  rrSlowFB   <- c(1,1,1)

  ############# CREATE A VECTOR FOR RATE OF SLOW PROGRESSION THAT WILL
  ############# VARY BASED ON LEVELS OF TB REACTIVATION RATES
  Vrslow     <- rep(1,4)
  ############# UPDATE LEVEL FOUR OF THE RATE OF SLOW BASED ON CALCULATED RR FROM
  ############# USER INPUTTED RR FOR THE RISK FACTOR
  Vrslow<-rslow*RRrslowRF

  TunrslowAge  <- PV["TunrslowAge"]
  rrReactAg       <- exp(c(0,0,0,0,0,0,0.5,1:4)*PV["TunrslowAge"])
  Mrslow <- outer(rrReactAg,Vrslow)

  #######################       RATE OF RECOVERY          ########################

  rRecov     <-  PV["rRecov"]/12

  #######################       RATE OF SELF CURE         ########################

  rSlfCur      <- PV["rSlfCur"]/12

  ######################          LTBI DIAGNOSIS           ########################

  ######################        TEST SPECIFICATIONS          ######################
  ####numbers from Stout paper
  Sens_IGRA <-c(.780,.675,.712,.789,.591)
  Spec_IGRA <-c(.979,.958,.989,.985,.931)
  IGRA_frc<-.50 #updated to 50% on 6/29/2021 per Suzanne
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
  ### In intervention 2, we would like to assume a smooth transition to 100% igra as the
  ### test of choice. In order to take account of this, we will transition the test
  ### sensitivities and specificities to IGRA levels.
  if (Int2 == 1) {
    for (i in 1:nrow(SensLt)){
      SensLt[i,] <- LgtCurve(intv_yr,intv_yr + 5,Sens_IGRA[i] - SensLt[i,]) + SensLt[i,]

      SpecLt[i,] <- LgtCurve(intv_yr, intv_yr + 5,Spec_IGRA[i] - SpecLt[i,]) + SpecLt[i,]
    } }


  ######################     IGRA FRACTION PROGRAM CHANGE    ########################
  # for (i in 1:nrow(SensLt)){
  #   SensLt[i,913:ncol(SensLt)]    <-rep(Sens_IGRA[i],1+month-913)
  #   SpecLt[i,913:ncol(SpecLt)]    <-rep(Spec_IGRA[i],1+month-913)
  # }
  if (prg_chng["IGRA_frc"] != IGRA_frc){
    for (i in 1:nrow(SensLt)){
      SensLt[i,prg_m:ncol(SensLt)]    <-rep((Sens_IGRA[i]*prg_chng["IGRA_frc"] + (1-prg_chng["IGRA_frc"])*Sens_TST[i]),1+month-prg_m)
      SpecLt[i,prg_m:ncol(SpecLt)]    <-rep((Spec_IGRA[i]*prg_chng["IGRA_frc"] + (1-prg_chng["IGRA_frc"])*Spec_TST[i]),1+month-prg_m)
    }
  }

  ### ADJUST THIS FOR THE FOREIGN BORN
  ######################        SCREENING RATES          ######################
  rLtScrt       <- cbind(LgtCurve(1985,2015,(PV["rLtScr"]))/12, LgtCurve(1985,2015,(PV["rLtScr"]))/12)
  ######################  SCREENING RATE PROGRAM CHANGE ########################
  if (prg_chng["scrn_cov"] !=1) {
    rLtScrt[prg_m:nrow(rLtScrt),]<-rLtScrt[prg_m:nrow(rLtScrt),]*prg_chng["scrn_cov"];
  }

  ################################################################################
  #################### TTT ADDITIONAL SCREENING PROBABILITIES ####################
  ################################################################################
  ###    default sampling distribution is set to zero & ltbi prevalence to 1   ###
  ttt_sampling_dist<-matrix(0,length(ttt_list),22*16)
  ttt_ltbi<-as.vector(rep(1,length(ttt_list)))
  ###  if the population & fraction to be screened are non zero update values  ###
  if (ttt_list[[1]][[3]]!=0 & ttt_list[[1]][[4]]!=0){
    ###this is dependent on a basecase run so load in that data
    load(system.file(paste0(loc, "/", loc, "_results_1.rda"), package = "MITUS"))
    ttt_params<-create_ttt_mdist(ttt_list = ttt_list,
                                 results = out[3,,],
                                 PV = PV)
    for (i in 1:length(ttt_params)){
      ttt_sampling_dist[i,]<-as.vector(ttt_params[[i]][,-(ncol(ttt_params[[i]]))])
      ttt_ltbi[i]<-ttt_list[[i]][[9]]
    }
  }

################################################################################
################################################################################
################################################################################

  ###adjustments to the screening rates dependent on risk and TB status
  rrTestHr    <- PV["rrTestHr"] # RR of LTBI screening for HIV and HR as cmpared to general
  rrTestLrNoTb  <- PV["rrTestLrNoTb"] # RR of LTBI screening for individuals with no risk factors
  ### because of the introduction of new time varying parameters, we will create 2 matrices to
  ### hold the three different sensitivity and specificity measures; one will be for those whose
  ### true LTBI status is positive and the other is for those whose true TB status is negative.
  LtDxPar_nolt <- LtDxPar_lt <- matrix(NA,nrow(SensLt),month);

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
  rownames(LtDxPar_lt) <- rownames(LtDxPar_nolt) <- c("US","HR.US","youngNUS","NUS","HR.NUS")


  #### #### #### INT 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  #### #### #### INT 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  pctDoc <- (1-0.28)
  Int1Test <- Int1Init <- Int1Tx <- matrix(0,101,11)
  if(Int1==1) {
    ##### In order to properly cost this intervention, we need to record the number of
    ##### to individuals that transition TB states during this intervention
    ImmNon0 <- ImmNon; ImmLat0 <- ImmLat; ImmFst0 <- ImmFst; ImmAct0 <- ImmAct;
    ImmLatTx <- ImmActTx <- ImmFstTx <- ImmLatTest <- ImmActTest <- ImmFstTest <- ImmNonTest <- matrix(0, month, 11)
    TxPars <- PV["EffLt"]*(1-PV["pDefLt"])
    ##### Create the time varying, age-stratified change vectors for the number of tests
    TotImmAgeDoc  <- TotImmAge[1:month,] *(LgtCurve(intv_yr,intv_yr+5,1)*pctDoc)
    Int1Test <- rowsum(TotImmAgeDoc[-month,],rep(1:((month-1)/12),each=12))
    for(i in 1:11){
      ImmLat[,i] <- ImmLat0[,i]*(1-LgtCurve(intv_yr,intv_yr+5,1)*SensLt[4,1]*TxPars*pctDoc)
      # ImmLatTest[,i] <- ImmLat0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*pctDoc)
      ImmLatTx[,i] <- ImmLat0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*SensLt[4,1]*pctDoc)

      ImmAct[,i] <- ImmAct0[,i]*(1-LgtCurve(intv_yr,intv_yr+5,1)*SensLt[4,1]*TxPars*pctDoc)
      # ImmActTest[,i] <- ImmAct0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*pctDoc)
      ImmActTx[,i] <- ImmAct0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*SensLt[4,1]*pctDoc)

      ImmFst[,i] <- ImmFst0[,i]*(1-LgtCurve(intv_yr,intv_yr+5,1)*SensLt[4,1]*TxPars*pctDoc)
      # ImmFstTest[,i] <- ImmFst0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*pctDoc)
      ImmFstTx[,i] <- ImmFst0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*SensLt[4,1]*pctDoc)
      # ImmNonTest[,i] <- ImmNon0[,i]*(LgtCurve(intv_yr,intv_yr+5,1)*pctDoc)
    }
    ImmNon      <- TotImmAge[1:month,]-ImmAct-ImmFst-ImmLat
    for (i in 1:length(ImmNon)) if (ImmNon[i] < 0) ImmNon[i]<-0
    ##### Create the time varying, age-stratified change vectors for the number of inits and completions
    Int1Init <- rowsum(ImmLatTx[-month,],rep(1:((month-1)/12),each=12)) + rowsum(ImmActTx[-month,],rep(1:((month-1)/12),each=12)) +
      rowsum(ImmFstTx[-month,],rep(1:((month-1)/12),each=12))
    Int1Tx <- Int1Init * (1-PV["pDefLt"])
  }
  #### #### #### INT 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####


  ######################          LTBI DIAGNOSIS           ########################
  ##read in the default program change values
  default_pc<-def_prgchng(PV)
  ######################          LTBI DIAGNOSIS           ########################
  ### PROBABILITY OF LATENT TREATMENT INTIATION
  pTlInt        <- rep(.773,month)
  ###################### LTBI TX INITIATION PROGRAM CHANGE ########################
  if (prg_chng["ltbi_init_frc"] !=pTlInt[prg_m]){
    pTlInt[prg_m:length(pTlInt)] <- prg_chng["ltbi_init_frc"];
  }
  ########################################################################################
  ##### ##### LTBI TREATMENT COMPLETION RATE
  ########################################################################################
  bccomprate<-((default_pc[["frc_3hp"]]*default_pc[["comp_3hp"]])+(default_pc[["frc_3hr"]]*default_pc[["comp_3hr"]])+
                 (default_pc[["frc_4r"]]*default_pc[["comp_4r"]]))
  # print(bccomprate)
  if (prg_chng[["frc_3hp"]] != default_pc[["frc_3hp"]] | prg_chng[["comp_3hp"]] != default_pc[["comp_3hp"]] |
      prg_chng[["frc_3hr"]] != default_pc[["frc_3hr"]] | prg_chng[["comp_3hr"]] != default_pc[["comp_3hr"]] |
      prg_chng[["frc_4r"]] != default_pc[["frc_4r"]] | prg_chng[["comp_4r"]] != default_pc[["comp_4r"]]){
    #calculate the weighted treatment completion rate
    comprate<-(prg_chng[["frc_3hp"]]*prg_chng[["comp_3hp"]]+prg_chng[["frc_3hr"]]*prg_chng[["comp_3hr"]]+
                 prg_chng[["frc_4r"]]*prg_chng[["comp_4r"]])/1
    pDefLt         <- c(rep(1-bccomprate,prg_m-1),rep(1-comprate,1+month-prg_m))
  } else {
    pDefLt         <- rep(1-bccomprate,month)
  }

  ########################################################################################
  ##### ##### LTBI TREATMENT EFFECTIVENESS
  ########################################################################################

  EffLt         <- rep(PV["EffLt"],month)

  ### because of the introduction of new time varying parameters, we will create a matrix to
  ### hold the latent treatment parameters
  LtTxPar       <- matrix(NA,3,month)
  LtTxPar       <- cbind(pTlInt,pDefLt,EffLt)

  pImmScen    <- PV["pImmScen"] # lack of reactivitiy to IGRA for Sp



  #### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  if(Int2==1) { rLtScrt     <- rLtScrt  + LgtCurve(intv_yr,intv_yr+5,1)*rLtScrt*1}
  pImmScen   <- PV["pImmScen"] # lack of reactivitiy to IGRA for Sp

  #### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  ### HOW TO ADD PROGRAM CHANGE HERE?

  if(Int2==1) {
    ## INCREASE SCREENING AMONG THE HIGH RISK POPULATION
    rLtScrt[,2]     <- rLtScrt[,2]  + LgtCurve(intv_yr,intv_yr+5,1)*rLtScrt[,2]
    ## INCREASE THE EFFECTIVENESS OF TREATMENT
    LtTxPar[,3] <-    LtTxPar[,3]+LgtCurve(intv_yr,intv_yr+5, (1-EffLt[intv_m])/2)
    ## INCREASE TREATMENT COMPLETION RATE
    LtTxPar[,2] <-    LtTxPar[,2]-LgtCurve(intv_yr,intv_yr+5, pDefLt[intv_m]/2)

  }
  pImmScen   <- PV["pImmScen"] # lack of reactivitiy to IGRA for Sp

  ################################################################################
  #######################         TB DIAGNOSIS            ########################
  #######################      TEST CHARACTERISTICS       ########################

  SensSp    <- PV["SensSp"]

  ######################           PROVIDER DELAY         ########################
  DelaySp    <- rep(PV["DelaySp"],month)

  #add in a temporary change to provider delay
  if (delay == 1){
    DelaySp[841:859]<-8*DelaySp[841:859]
  }

  if (prg_chng["tb_tim2tx_frc"] !=100){
    DelaySp[prg_m:length(DelaySp)] <- PV["DelaySp"]*(prg_chng["tb_tim2tx_frc"])/100;
  }

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
  rDxt1<-rDxt1[1:month,]

  # Put it all together
  rDxt           <- 1/(1/rDxt1+DelaySp)*SensSp
  rDxt[,2]       <- (rDxt[,1]-min(rDxt[,1]))/PV["rrDxH"]+min(rDxt[,1]) #check this with Nick
  colnames(rDxt) <- c("Active","Active_HighRisk")

  #### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Int3==1) { for(i in 1:2) { rDxt[,i] <- rDxt[,i]+ rDxt[,i]*LgtCurve(intv_yr,intv_yr+5,.1)     }   }
  if(Int3==2) { for(i in 1:2) { rDxt[,i] <- rDxt[,i]- rDxt[,i]*LgtCurve(intv_yr,intv_yr+5,.1)     }   }

  #### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  ######  TREATMENT OUTCOMES  ######
  TunTxMort	<- PV["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

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
  NixTrans<- rep(1,month)
  NixTb<- rep(0,month)

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Int4==1) {  rDeftH <- rDeftH*(1-LgtCurve(intv_yr,intv_yr+5,0.5))
  for (i in 1:length(rDeftH)) if (rDeftH[i] > 1) rDeftH[i] <-1
  }
  if(Int4==1) {  rDeft<-rDeft[1:month]
  rDeft<- rDeft * (1-LgtCurve(intv_yr,intv_yr+5,0.5))
  for (i in 1:length(rDeft)) if (rDeft[i] > 1) rDeft[i] <-1
  }

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  ## Tx quality
  TxQual0         <- rep(NA,151)
  TxQual0[1:30]   <- PV["TxQualEarly"]
  TxQual0[44:62]  <- ORAdd(TxInputs[[2]][,2],PV["TunTxQual"])
  TxQual0[63:151] <- TxQual0[62]
  TxQual1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=TxQual0[-(31:43)],spar=0.4),x=1950:2100)$y
  TxQualt        <- SmoCurve(TxQual1);
  TxQualt         <-TxQualt[1:month]
  RRcurDef      <- PV["RRcurDef"]

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(Int4 == 1) {
    TxQualt<- 1 - (1 - TxQualt)*(1 - LgtCurve(intv_yr, intv_yr + 5, 0.5))
    for (i in 1:length(TxQualt)) if (TxQualt[i] > 1) TxQualt[i] <-1
  }

  #### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  #### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  NixTrans <- 1 - LgtCurve(intv_yr,intv_yr + 1, 1)
  if(Scen1 == 0) { NixTrans[] <- 1}
  if(Scen1 == 2) { Nixtrans <- c(rep(1, intv_m), rep(0, month - (intv_m)))}

  #### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  #### #### #### SCEN 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if (Scen4 == 1){ NixTb <- c(rep(0, intv_m), 1, rep(0, month - (intv_m - 1)))}

  #### #### #### SCEN 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  ## Retreatment
  pReTx   <- LgtCurve(1985,2000,PV["pReTx"])   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)

  ### REMOVED ACQUIRED RESISTANCE PARAMETERS
  ### REPLACED TX MAT WITH VECTOR FOR SINGLE VALUES
  #####################   NEW TB TREATMENT TABLE ################
  TxVec          <- rep(NA,2)
  names(TxVec) <- c("TxCompRate","TxEff")
  TxVec[1]       <-  d1st
  TxVec[2]       <-  pCurPs
  ### REMOVED ART PARAMETERS

  for(i in 1:2){
    for (j in 1:nrow(rDxt)){
      if (rDxt[j,i]>1) {rDxt[j,i]<-1
      }   }}


  InputParams<-list()
  InputParams[["rDxt"]]      = rDxt
  InputParams[["TxQualt"]]   = TxQualt
  InputParams[["InitPop"]]   = InitPop
  InputParams[["Mpfast"]]    = Mpfast
  InputParams[["ExogInf"]]   = ExogInf
  InputParams[["MpfastPI"]]  = MpfastPI
  InputParams[["Mrslow"]]    = Mrslow
  InputParams[["rrSlowFB"]]  = rrSlowFB
  InputParams[["rfast"]]     = rfast
  InputParams[["RRcurDef"]]  = RRcurDef
  InputParams[["rSlfCur"]]   = rSlfCur
  InputParams[["p_HR"]]      = p_HR
  InputParams[["dist_gen"]]  = dist_gen
  InputParams[["vTMort"]]    = vTMort
  InputParams[["RRmuRF"]]    = RRmuRF
  InputParams[["RRmuHR"]]    = RRmuHR
  InputParams[["Birthst"]]   = Birthst
  InputParams[["HrEntEx"]]   = HrEntEx
  InputParams[["ImmNon"]]    = ImmNon
  InputParams[["ImmLat"]]    = ImmLat
  InputParams[["ImmAct"]]    = ImmAct
  InputParams[["ImmFst"]]    = ImmFst
  InputParams[["Int1Test"]]  = Int1Test
  InputParams[["Int1Init"]]  = Int1Init
  InputParams[["Int1Tx"]]  = Int1Tx
  InputParams[["mubt"]]      = mubt
  InputParams[["RelInf"]]    = RelInf
  InputParams[["RelInfRg"]]  = RelInfRg
  InputParams[["RRcrAG"]]    = RRcrAG
  InputParams[["Vmix"]]      = Vmix
  InputParams[["rEmmigFB"]]  = rEmmigFB
  InputParams[["TxVec"]]     = TxVec
  InputParams[["TunTxMort"]] = TunTxMort
  InputParams[["rDeft"]]     = rDeft
  InputParams[["pReTx"]]     = pReTx
  InputParams[["LtTxPar"]]   = LtTxPar
  InputParams[["LtDxPar_lt"]]   = LtDxPar_lt
  InputParams[["LtDxPar_nolt"]]   = LtDxPar_nolt
  InputParams[["rrTestLrNoTb"]] = rrTestLrNoTb
  InputParams[["rrTestHr"]] = rrTestHr
  InputParams[["ttt_month"]]<-ttt_month
  InputParams[["ttt_sampling_dist"]]<-ttt_sampling_dist
  InputParams[["ttt_na"]]<-0
  InputParams[["ttt_ag"]]<-0
  InputParams[["ttt_pop_scrn"]]<-1
  InputParams[["ttt_ltbi"]]<-ttt_ltbi
  InputParams[["rLtScrt"]]   = rLtScrt
  InputParams[["RRdxAge"]]   = RRdxAge
  InputParams[["rRecov"]]    = rRecov
  InputParams[["pImmScen"]]  = pImmScen
  InputParams[["EarlyTrend"]]= EarlyTrend
  InputParams[["net_mig_usb"]]  = net_mig_usb
  InputParams[["net_mig_nusb"]]= net_mig_nusb
  InputParams[["aging_denom"]] <-spl_den
  InputParams[["adj_fact"]] <- adj_fact
  InputParams[["NixTrans"]] <- NixTrans
  InputParams[["NixTb"]] <- NixTb
  InputParams[["SpImmNon"]]  = SpImmNon

  InputParams[["ResNam"]]    <- func3_ResNam()
  return(InputParams)

  ###################################################
}
