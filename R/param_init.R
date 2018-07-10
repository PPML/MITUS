################################################################################
##### THE CODE BELOW WILL GATHER AND FORMAT INPUTS FOR THE TB_MODEL        #####
##### FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME; #####
##### VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.        #####
################################################################################

#' This function takes the same inputs as the Outputs
#'@name param_init
#'@param ParVector vector of Inputs to format
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@return InputParams list
#'@export

param_init <- function(ParVector,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0){
data("ModelInputs_9-2-16", package='MITUS')
#'Create an empty list to hold the formatted intitial parameters
InputParams <-vector("list", 45)
names(InputParams) <- c("rDxt","TxQualt", "InitPop", "Mpfast", "ExogInf", "MpfastPI",
                        "Mrslow", "rrSlowFB", "rfast"    ,"RRcurDef"      , "rSlfCur"  ,
                        "p_HR"        , "dist_gen" , "vTMort"   ,"RRmuRF"          , "RRmuHR",
                       "muTbRF"        , "Birthst" , "HrEntEx"  ,"ImmNon"       , "ImmLat"  ,
                        "ImmAct"      , "ImmFst" , "mubt"     ,"RelInf"        , "RelInfRg" ,
                        "Vmix"       , "rEmmigFB" , "TxVec"    , "TunTxMort"    , "rDeft"    ,
                       "pReTx"      , "LtTxPar"  , "LtDxPar"  , "rLtScrt"      , "RRdxAge"  ,
                        "rRecov"      , "pImmScen"  ,   "EarlyTrend", "NixTrans"    , "can_go"  ,
                        "dist_goal"      ,  "diff_i_v" ,"dist_orig_v", "ResNam")
################################################################################
##### INTERVENTION
################################################################################
if(Int5==1) {
  Int1 = Int2 = Int3 = Int4 = 1
}

#'params from general reblnc_pop:
InputParams[["dist_gen"]] <- dist_gen;
InputParams[["can_go"]] <- can_go;
InputParams[["dist_goal"]] <-dist_goal;
InputParams[["diff_i_v"]] <-diff_i_v;
InputParams[["dist_orig_v"]] <-dist_orig_v;
################################################################################
###########################          INPUTS            #########################
################################################################################
BgMort                   <- Inputs[["BgMort"]]
InputParams[["InitPop"]] <- Inputs[["InitPop"]]
Births                   <- Inputs[["Births"]]
ImmigInputs              <- Inputs[["ImmigInputs"]]
TxInputs                 <- Inputs[["TxInputs"]]

##########                PARAMETER DEFINITIONS                      ###########
#######################           BIRTHS                 #######################
####### INDEXED BY TIME, ABSOLUTE NUMBER OF NEW ADULT ENTRANTS OVER TIME #######

InputParams[["Birthst"]]  <- SmoCurve(Births)*P["TunBirths"]/12

##########################      MORTALITY RATES       ##########################
########################## BACKGROUND MORTALITY BY TIME ########################
InputParams[["mubt"]]    <- matrix(NA,1801,11)
for(i in 1:11) {
  InputParams[["mubt"]][,i] <- SmoCurve(BgMort[,i+1])*P["TunMubt"]/12
}
#########################     DISEASE SPECIFIC       ###########################
#############    ACTIVE TB RATES DEFAULT TO THE SMEAR POS LEVELS   #############
muIp  	  <- P["muIp"]/12
############ CONVERT ANNUAL RATES OF RF MORTALITY TO MONTHLY RATES ##########

############ THESE MUST BE UPDATED
muRF1      <- P["muH1"]/12*.01
muRF2      <- (P["muH2"]+P["muH1"])/24*.01
muRF3      <- P["muH2"]/12*.01
InputParams[["muTbRF"]]    <- P["muTbH"]/12
TunmuHvAg <- P["TunmuTbAg"] # ffs.

###############  RATE RATIO OF MORTALITY INCREASE FOR HIGH RISK ###############

InputParams[["RRmuHR"]]    <- c(1,P["RRmuHR"])

############### CREATE A MATRIX OF RF MORTALITIES BY AGE GROUP ###############
mort_dist<-rowSums(InputParams[["dist_goal"]])

RF_fact=2

InputParams[["RRmuRF"]]   <- rep(NA,4);
names(InputParams[["RRmuRF"]]) <- c("RF1","RF2","RF3","RF4")

InputParams[["RRmuRF"]]<-exp((0:3)/3*log(RF_fact))
# RRmuRF<-RRmuRF/sum(RRmuRF*mort_dist)

#check =1
# RRmuRF%*%mort_dist

# vRFMort    <- matrix(0,11,4);
# rownames(vRFMort) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
# colnames(vRFMort) <- c("RF1","RF2","RF3","RF4")
# vRFMort[,1] <- 0;
# vRFMort[,2] <- muRF1*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
# vRFMort[,3] <- muRF2*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)
# vRFMort[,4] <- muRF3*exp(c(0,0,1:6,6,6,6)*TunmuHvAg)

##### combine the two RRs to a single factor
# InputParams["RRs_mu"]<-matrix(NA,2,4)
# InputParams["RRs_mu"][1,] <-RRmuRF
# InputParams["RRs_mu"][2,] <-RRmuRF*RRmuHR[2]

######################## MULTIPLER OF MORT RATE ABOVE ########################

TunmuTbAg <- P["TunmuTbAg"]

############### CREATE A MATRIX OF TB MORTALITIES BY AGE GROUP ###############

InputParams[["vTMort"]]  <- matrix(0,11,6);
rownames(InputParams[["vTMort"]]) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
colnames(InputParams[["vTMort"]]) <- c("Su","Sp","Ls","Lf","Ac","Tx")
InputParams[["vTMort"]][,5:6] <- muIp #active disease rates default to smear positive
for(i in 1:11) {
  InputParams[["vTMort"]][i,] <-  InputParams[["vTMort"]][i,]*c(3,1+0:9*TunmuTbAg)[i]
}


######################## MULTIPLER OF MORT RATE ABOVE ########################

TunmuTbAg <- P["TunmuTbAg"]

#################                IMMIGRATION              #####################
 TotImmig1       <- c(ImmigInputs[[1]][1:65],(ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*P["ImmigVolFut"]+ImmigInputs[[1]][66])/12*P["ImmigVol"]
 TotImmig0       <- (c(ImmigInputs[[1]][1:151])+c(rep(0,65),cumsum(rep(P["ImmigVolFut"],86))))/12*P["ImmigVol"]
#set.seed( rand_seed)
# TotImmig1       <-   c(TotImmig0[1:66],exp(mvrnorm(1, log(TotImmig0[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
TotImmig1       <-   TotImmig0
TotImmig        <- SmoCurve(TotImmig1)
TotImmAge       <- outer(TotImmig,ImmigInputs[["AgeDist"]])
#######################   IMMIGRATION WITH LATENT TB   #######################

PrevTrend25_340l <- c(ImmigInputs[["PrevTrend25_34"]][1:65]^P["TunLtbiTrend"]*ImmigInputs[["PrevTrend25_34"]][65]^(1-P["TunLtbiTrend"]),
                      ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutLat"]/0.99)^(1:86))

#set.seed(rand_seed+1)
# PrevTrend25_341l <-   c(PrevTrend25_340l[1:66],exp(mvrnorm(1, log(PrevTrend25_340l[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
PrevTrend25_341l <-   PrevTrend25_340l
PrevTrend25_34l  <- SmoCurve(PrevTrend25_341l)
PrevTrend25_34_ls <- (PrevTrend25_34l); PrevTrend25_34_ls <- PrevTrend25_34_ls/PrevTrend25_34_ls[(2011-1950)*12+6]
InputParams[["ImmLat"]]        <- matrix(NA,length(PrevTrend25_34_ls),11)
for(i in 1:11) InputParams[["ImmLat"]][,i] <- (1-exp((-(c(2.5,1:9*10,100)/100)[i]*P["LtbiPar1"]-(c(2.5,1:9*10,100)/100)[i]^2*P["LtbiPar2"])*PrevTrend25_34_ls))*TotImmAge[,i]

#######################   IMMIGRATION WITH ACTIVE TB   #######################
PrevTrend25_340a <- c(ImmigInputs[["PrevTrend25_34"]][1:65],ImmigInputs[["PrevTrend25_34"]][66:151]*(P["ImmigPrevFutAct"]/0.99)^(1:86))
#set.seed( rand_seed+2)
# PrevTrend25_341a <-   c(PrevTrend25_340a[1:66],exp(mvrnorm(1, log(PrevTrend25_340a[-(1:66)]), vcv_gp_l10_sd0.1))) # Add gaussian process noise
PrevTrend25_341a <-   PrevTrend25_340a
PrevTrend25_34a  <- SmoCurve(PrevTrend25_341a)
ImDxChngV      <- SmoCurve(c(rep(1,57),seq(1,P["ImDxChng"],length.out=6)[-1],rep(P["ImDxChng"],89)))
InputParams[["ImmAct"]]        <- outer(PrevTrend25_34a*P["RRtbprev"]*ImDxChngV,ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*P["pImAct"]
InputParams[["ImmFst"]]        <- outer(PrevTrend25_34a*P["RRtbprev"],ImmigInputs[["RR_Active_TB_Age"]])*TotImmAge*(1-P["pImAct"])
InputParams[["ImmNon"]]        <- TotImmAge-InputParams[["ImmAct"]]-InputParams[["ImmFst"]]-InputParams[["ImmLat"]]

#### #### #### INT 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
pctDoc <- (1-0.28)
if(Int1==1) {
  for(i in 1:11) InputParams[["ImmLat"]][,i] <- InputParams[["ImmLat"]][,i]*(1-LgtCurve(2016,2021,1)*P["SensLt"]*P["EffLt"]*(1-P["pDefLt"])*pctDoc)
  for(i in 1:11) InputParams[["ImmAct"]][,i] <- InputParams[["ImmAct"]][,i]*(1-LgtCurve(2016,2021,1)*P["SensLt"]*P["EffLt"]*(1-P["pDefLt"])*pctDoc)
  for(i in 1:11) InputParams[["ImmFst"]][,i] <- InputParams[["ImmFst"]][,i]*(1-LgtCurve(2016,2021,1)*P["SensLt"]*P["EffLt"]*(1-P["pDefLt"])*pctDoc)
  InputParams["ImmNon"]         <- TotImmAge-InputParams[["ImmAct"]]-InputParams[["ImmFst"]]-InputParams[["ImmLat"]]   }


#### #### #### SCEN 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

if(Scen2==1) {
  for(i in 1:11) InputParams["ImmLat"][,i] <- InputParams[["ImmLat"]][,i]*(1-LgtCurve(2016,2017,1))
  for(i in 1:11) InputParams["ImmAct"][,i] <- InputParams[["ImmAct"]][,i]*(1-LgtCurve(2016,2017,1))
  for(i in 1:11) InputParams["ImmFst"][,i] <- InputParams[["ImmFst"]][,i]*(1-LgtCurve(2016,2017,1))
  InputParams[["ImmNon"]]         <- TotImmAge-InputParams[["ImmAct"]]-InputParams[["ImmFst"]]-InputParams[["ImmLat"]]    }


#### #### #### SCEN 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

if(Scen3==1) {
  adjf <- (1-0.5)^(1/10/12);  adjfV <- adjf^(0:1008)
  for(i in 1:11) InputParams[["ImmLat"]][793:1801,i] <- InputParams[["ImmLat"]][793:1801,i]*adjfV
  for(i in 1:11) InputParams[["ImmAct"]][793:1801,i] <- InputParams[["ImmAct"]][793:1801,i]*adjfV
  for(i in 1:11) InputParams[["ImmFst"]][793:1801,i] <- InputParams[["ImmFst"]][793:1801,i]*adjfV
  InputParams["ImmNon"]         <- TotImmAgeInputParams[["ImmAct"]]-InputParams[["ImmFst"]]-InputParams[["ImmLat"]] }

#######################__EXOGENEOUS INFECTION RISK    #######################

# NEED TO REMOVE DRUG RESISTANCE
InputParams[["ExogInf"]] <- rep(NA,length(PrevTrend25_34a))
InputParams[["ExogInf"]] <- P["ExogInf"]*PrevTrend25_34a/PrevTrend25_341a["2013"]/12
#removed *(ImmigInputs[[7]][4]*DrN[,i]+(1-ImmigInputs[[7]][4])*DrE[,i])
###############################    EMMIGRATION   ##############################

InputParams[["rEmmigFB"]] <- c(P["rEmmigF1"],P["rEmmigF2"])/12


##########################  HIGH-RISK ENTRY/EXIT  #############################

InputParams[["p_HR"]]    <- P["pHR"]
yr       <- c(2.5,1:10*10)
r0_5     <- 1/3;
r45_55 <- 1/20
HR_exit  <- r0_5*((r45_55/r0_5)^(1/(50-2.5)))^(yr-2.5)
HR_entry <- HR_exit*InputParams[["p_HR"]]*1.3
InputParams[["HrEntEx"]]  <- cbind(HR_entry,HR_exit)/12

##########################   TB TRANSMISSION  #################################

CR           <- P["CR"]/12
TrIn         <- P["TrIn"]	# Contact rate for In as a fraction of Ip
InputParams[["RelInfRg"]]    <- c(1.0,P["RelCrHr"], 1.0, P["RelCrHr"])*CR
TunTbTransTx <- 0 #P["TunTbTransTx"]  # set to zero?
InputParams[["Vmix"]]         <- 1-c(P["sigmaHr"],P["sigmaFb"])
RelInf <- rep(0,6)
names(RelInf) <-c("Su","Sp","Ls","Lf","Ac","Tx")

RelInf[5] <- 1;
RelInf[6] <- RelInf[5]*TunTbTransTx
InputParams[["RelInf"]] <-RelInf

#########################   TB NATURAL HISTORY  ###############################
#########################__   EARLY EPIDEMTIC   ###############################

Early0 <- P["Early0"]
InputParams[["EarlyTrend"]] <- c(rep(1+Early0,75*12),seq(1+Early0,1.0,length.out=25*12+2))

######################     PROGRESSION TO DISEASE     ##########################

pfast      <- P["pfast"]
pimmed    <- P["pimmed"]
ORpfast1   <- P["ORpfast1"] ## age group 1
ORpfast2   <- P["ORpfast2"] ## age group 2
ORpfastRF  <- P["ORpfastH"] ##riskfactor
ORpfastPI  <- P["ORpfastPI"]
rslow      <- P["rslow"]/12
rslowRF    <- P["rslowH"]/12
RRrslowRF  <- rslowRF/rslow
InputParams[["rfast"]]    <- P["rfast"]/12
# rfast      <- ((1-pimmed)*rfast)+pimmed
rrSlowFB0  <- P["rrSlowFB"]
InputParams[["rrSlowFB"]]   <- c(1,rrSlowFB0,rrSlowFB0)
##############            ORIGINAL Mpfast[ag][hv]             ################
##############          CREATE NEW Mpfast[ag][im]               ##############
############## MIGHT WRITE A NEW SCRIPT FOR THIS PART
InputParams[["Mpfast"]]      <- matrix(NA,11,4)
############## CREATE AN ODDS FROM THE PROB OF FAST PROGRESSION ##############
InputParams[["Mpfast"]][,]    <- pfast/(1-pfast)
InputParams[["Mpfast"]][1,]   <- InputParams[["Mpfast"]][1,]*ORpfast1 # progression for age group 1
InputParams[["Mpfast"]][2,]   <- InputParams[["Mpfast"]][2,]*ORpfast2 # progression for age group 2


#################       CREATE A NEW MATRIX PARTIAL. IMM.     #################
InputParams[["MpfastPI"]]     <- InputParams[["Mpfast"]]

#vector of ORpfastRF
vORpfastPIRF<-vORpfastRF  <-c(1,1,1,1)
vORpfastRF  <-  (exp((0:3)/3*log(ORpfastRF)))
vORpfastPIRF  <- vORpfastRF*ORpfastPI

############ UPDATE PROBS FOR LEVEL 2 OF REACTIVATION ###########
InputParams[["Mpfast"]][,2]   <- vORpfastRF[2]*InputParams[["Mpfast"]][,2]
InputParams[["MpfastPI"]][,2]   <- vORpfastPIRF[2]*InputParams[["Mpfast"]][,2]
############ UPDATE PROBS FOR LEVEL 3 OF REACTIVATION ###########
InputParams[["Mpfast"]][,3]   <- vORpfastRF[3]*InputParams[["Mpfast"]][,3]
InputParams[["MpfastPI"]][,3]   <- vORpfastPIRF[3]*InputParams[["Mpfast"]][,3]
############ UPDATE PROBS FOR LEVEL 4 OF REACTIVATION ###########
InputParams[["Mpfast"]][,4]   <- vORpfastRF[4]*InputParams[["Mpfast"]][,4]
InputParams[["MpfastPI"]][,4]   <- vORpfastPIRF[4]*InputParams[["Mpfast"]][,4]

##### UPDATE BOTH MATRICES WITH PROBABILITIES, NOT ODDS
InputParams[["Mpfast"]][,]    <- InputParams[["Mpfast"]][,]  /(1+InputParams[["Mpfast"]][,]);
InputParams[["MpfastPI"]][,]  <- InputParams[["MpfastPI"]][,]/(1+InputParams[["MpfastPI"]][,]);

############# CREATE A VECTOR FOR RATE OF SLOW PROGRESSION THAT WILL
############# VARY BASED ON LEVELS OF TB REACTIVATION RATES
Vrslow     <- rep(NA,4)
############# UPDATE LEVEL FOUR OF THE RATE OF SLOW BASED ON CALCULATED RR FROM
############# USER INPUTTED RR FOR THE RISK FACTOR
Vrslow=rslow*exp((0:3)/3*log(RRrslowRF))

TunrslowAge  <- P["TunrslowAge"]
rrReactAg       <- exp(c(0,0,0,0,0,0,0.5,1:4)*P["TunrslowAge"])
InputParams[["Mrslow"]] <- outer(rrReactAg,Vrslow)
# InputParams[["Mrslow"]] <- InputParams["Mrslow"];


#######################       RATE OF RECOVERY          ########################

InputParams[["rRecov"]]     <-  P["rRecov"]/12

#######################       RATE OF SELF CURE         ########################

InputParams[["rSlfCur"]]      <- P["rSlfCur"]/12

##################### __         LTBI DIAGNOSIS           #######################
InputParams[["rLtScrt"]]       <- LgtCurve(1985,2015,P["rLtScr"])/12
SensLt        <- P["SensLt"]    #  sens of test for latent TB infection (based on IGRA QFT-GIT)
###removed sens for HIV
SpecLt        <- P["SpecLt"]    #  spec of test for latent TB infection (based on IGRA QFT-GIT)
SpecLtFb      <- SpecLt         #  spec of test for latent TB infection (based on IGRA QFT-GIT) in foreign-born (assumed BCG exposed)
###removed rrTestHIV
rrTestHr      <- P["rrTestHr"] # RR of LTBI screening for HIV and HR as cmpared to general
rrTestLrNoTb  <- P["rrTestLrNoTb"] # RR of LTBI screening for individuals with no risk factors
dLt           <- 1/9
rDefLt        <- dLt*P["pDefLt"]/(1-P["pDefLt"])  # based on 50% tx completion with 6 mo INH regimen 2.0 [1.0,3.0] from Menzies Ind J Med Res 2011
EffLt         <- P["EffLt"]
######NEW PARAMETER FOR MITUS MODEL
pTlInt        <- .80
InputParams[["LtTxPar"]]       <- c(pTlInt,P["pDefLt"],EffLt)

#dLtt          <- (1-LgtCurve(2016,2021,0)) * 1/9

EffLt0        <- LgtCurve(2016,2021,0)+1
EffLtX        <- cbind(EffLt0,0,EffLt0,0,0)

#### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

if(Int2==1) { InputParams[["rLtScrt"]]    <- InputParams[["rLtScrt"]] + LgtCurve(2016,2021,1)*InputParams[["rLtScrt"]]*1
EffLtX[,2] <- LgtCurve(2016,2021,1)
dLtt       <- (1-LgtCurve(2016,2021,1))*(1/9) + LgtCurve(2016,2021,1)*(1/3)  }

#### #### #### INT 2 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

InputParams[["LtDxPar"]] <- matrix(NA,3,2);
colnames(InputParams[["LtDxPar"]]) <- c("latent","no latent");
rownames(InputParams[["LtDxPar"]]) <- c("LR","HR","FB")
InputParams[["LtDxPar"]][,1] <- c(SensLt                 , rrTestHr*SensLt    , SensLt)
InputParams[["LtDxPar"]][,2] <- c(rrTestLrNoTb*(1-SpecLt), rrTestHr*(1-SpecLt), (1-SpecLtFb))

InputParams[["pImmScen"]]   <- P["pImmScen"] # lack of reactivitiy to IGRA for Sp

################################################################################
#######################         TB DIAGNOSIS            ########################
#######################      TEST CHARACTERISTICS       ########################

SensSp    <- P["SensSp"]

######################           PROVIDER DELAY         ########################

DelaySp    <- P["DelaySp"]

#######################     PROVIDER DELAY RR ELDERLY     ########################

TunRRdxAge    <- P["TunRRdxAge"]
InputParams[["RRdxAge"]]       <- 1+(c(rep(1,6),cumprod(seq(1.05,1.3,length.out=5)))-1)*TunRRdxAge

#######################         ATTENDANCE RATE           ########################

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

InputParams[["rDxt"]] <- 1/(1/rDxt1+DelaySp)*SensSp
InputParams[["rDxt"]][,2]       <- (InputParams[["rDxt"]][,1]-min(InputParams[["rDxt"]][,1]))/P["rrDxH"]+min(InputParams[["rDxt"]][,1]) #check this with Nick
colnames(InputParams[["rDxt"]]) <- c("Active","Active_HighRisk")

#### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

if(Int3==1) { for(i in 1:4) { InputParams[["rDxt"]][,i] <- InputParams[["rDxt"]][,i]+ InputParams[["rDxt"]][,i]*LgtCurve(2016,2021,1)     }   }

#### #### #### INT 3 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

## DST
pDstt       <- LgtCurve(1985,2000,P["pDst"])

#### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

if(Int4==1) {  pDstt <- 1-(1-pDstt)*(1-LgtCurve(2016,2021,0.5))      }

#### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

######  TREATMENT OUTCOMES  ######
InputParams[["TunTxMort"]]	<- P["TunTxMort"]	# Multiplier to adjust mortality rates while on treatment into reasonable range (based on observed data) 0 = no TB mort on TX

# Regimen duration (no uncertainty?)
d1st <- 1/9

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
InputParams[["rDeft"]]         <- SmoCurve(rDef1)/12;
# rDeftH        <- rDeft*P["RRdefHR"] # second col is HR default rate

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
InputParams[["TxQualt"]]         <- SmoCurve(TxQual1);
InputParams[["RRcurDef"]]        <- P["RRcurDef"]

#### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

if(Int4==1) {  InputParams[["TxQualt"]] <- 1-(1-InputParams[["TxQualt"]])*(1-LgtCurve(2016,2021,0.5))      }

#### #### #### INT 4 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

InputParams[["NixTrans"]] <- 1-LgtCurve(2016,2017,1)
if(Scen1==0) {  InputParams[["NixTrans"]][] <- 1      }

#### #### #### SCEN 1 #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

## Retreatment
InputParams[["pReTx"]]   <- LgtCurve(1985,2000,P["pReTx"])   	# Probability Tx failure identified, patient initiated on tx experienced reg (may be same)

### REMOVED ACQUIRED RESISTANCE PARAMETERS
### REPLACED TX MAT WITH VECTOR FOR SINGLE VALUES
#####################   NEW TB TREATMENT TABLE ################
InputParams[["TxVec"]]           <- rep(NA,2)
names(InputParams[["TxVec"]]) <- c("TxCompRate","TxEff")
InputParams[["TxVec"]][1]       <-  d1st
InputParams[["TxVec"]][2]       <-  pCurPs
### REMOVED ART PARAMETERS
################################################################################
##### CREATE A LIST TO HOLD THE VECTORS FOR AGE CATEGORIES, TB STATES,     #####
##### DRUG RESISTANCE, TREATMENT HISTORY, HIV STATUS, AND RISK CATEGORY.   #####
################################################################################
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

            paste("mort_rate ag 1",StatList[[5]],sep="_" ),
            paste("mort_rate ag 2",StatList[[5]],sep="_" ),
            paste("mort_rate ag 3",StatList[[5]],sep="_" ),
            paste("mort_rate ag 4",StatList[[5]],sep="_" ),
            paste("mort_rate ag 5",StatList[[5]],sep="_" ),
            paste("mort_rate ag 6",StatList[[5]],sep="_" ),
            paste("mort_rate ag 7",StatList[[5]],sep="_" ),
            paste("mort_rate ag 8",StatList[[5]],sep="_" ),
            paste("mort_rate ag 9",StatList[[5]],sep="_" ),
            paste("mort_rate ag 10",StatList[[5]],sep="_" ),
            paste("mort_rate ag 11",StatList[[5]],sep="_" )
)

InputParams[["ResNam"]]<-ResNam
return(InputParams)

###################################################
}
