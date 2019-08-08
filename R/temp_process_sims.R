#' THIS FUNCTION INPUTS A TABLE OF PARAMETERS AND RUNS THE TB MODEL
#' AND GENERATES AN ARRAY OF OUTPUTS



#load("~/MITUS/data/parAll200_9-14-16.rData")
# load("data/ParamInit_2018.rData")
#'function to run the model
#'@name new_OutputsZint
#'@param samp_i how many samples
#'@param ParMatrix parameters to use in the simulation
#'@param loc
#'@param endyr year to end the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param prg_chng vector of program change values
#'@param ttt_list
#'@return results data frame of output
#'@export
new_OutputsZint <-  function(samp_i=1,ParMatrix,loc, startyr=1950, endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0,prg_chng, ttt_list) {
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par <- as.numeric(ParMatrix);
    names(Par) <- names(ParMatrix)
  } else {  Par <- as.numeric(ParMatrix[samp_i,]);
  names(Par) <- colnames(ParMatrix) }

  P <- Par
  Int1 <<- Int1;
  Int2 <<- Int2;
  Int3 <<- Int3;
  Int4 <<- Int4;
  Int5 <<- Int5;
  Scen1 <<- Scen1;
  Scen2 <<- Scen2;
  Scen3 <<- Scen3

  # P <- Par
  prms <-list()
  prms <- fin2_param(P,loc,prg_chng, ttt_list)
  IP <- list()
  IP <- fin2_param_init(P,loc,Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3,prg_chng, ttt_list)
  # data("trans_mat_nat",package="MITUS")
  # trans_mat_tot_ages<-trans_mat_tot_ages_nat
  # tm<-matrix(0,16,16)
  # diag(tm)<-1
  # trans_mat_tot_ages<<-matrix(tm,16,176)
  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  m <- fin2_cSim( nYrs       = 2050-1950         , nRes      = length(IP[["ResNam"]])  , rDxt     = IP[["rDxt"]]  , TxQualt    = IP[["TxQualt"]]   , InitPop  = IP[["InitPop"]]    ,
                  Mpfast     = IP[["Mpfast"]]    , ExogInf   = IP[["ExogInf"]]       , MpfastPI = IP[["MpfastPI"]], Mrslow     = IP[["Mrslow"]]    , rrSlowFB = IP[["rrSlowFB"]]  ,
                  rfast      = IP[["rfast"]]     , RRcurDef  = IP[["RRcurDef"]]      , rSlfCur  = IP[["rSlfCur"]] , p_HR       = IP[["p_HR"]]      , dist_gen = IP[["dist_gen"]]    ,
                  vTMort     = IP[["vTMort"]]    , RRmuRF    = IP[["RRmuRF"]]        , RRmuHR   = IP[["RRmuHR"]]  , Birthst  = IP[["Birthst"]]    ,
                  HrEntEx    = IP[["HrEntEx"]]   , ImmNon    = IP[["ImmNon"]]        , ImmLat   = IP[["ImmLat"]] , ImmAct     = IP[["ImmAct"]]    , ImmFst   = IP[["ImmFst"]]    ,
                  net_mig_usb = IP[["net_mig_usb"]], net_mig_nusb = IP[["net_mig_nusb"]],
                  mubt       = IP[["mubt"]]    , RelInf    = IP[["RelInf"]]        , RelInfRg = IP[["RelInfRg"]], Vmix       = IP[["Vmix"]]      , rEmmigFB = IP [["rEmmigFB"]]  ,
                  TxVec      = IP[["TxVec"]]     , TunTxMort = IP[["TunTxMort"]]     , rDeft    = IP[["rDeft"]]   , pReTx      = IP[["pReTx"]]     , LtTxPar  = IP[["LtTxPar"]]    ,
                  LtDxPar_lt    = IP[["LtDxPar_lt"]]   , LtDxPar_nolt    = IP[["LtDxPar_nolt"]]   , rLtScrt   = IP[["rLtScrt"]]       , ttt_samp_dist   = IP[["ttt_sampling_dist"]] ,
                  ttt_ag = IP[["ttt_ag"]], ttt_na = IP[["ttt_na"]], ttt_month = IP[["ttt_month"]], RRdxAge  = IP[["RRdxAge"]] , rRecov     = IP[["rRecov"]]    , pImmScen = IP[["pImmScen"]]   ,
                  EarlyTrend = IP[["EarlyTrend"]], ag_den=IP[["aging_denom"]],  NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)$Outputs
  colnames(m) <- func2_ResNam();
  results<<-as.matrix(m)

  return(results)
}


#'wrapper function for the above function
#'@name new_OutputsInt
#'@param ParMatrix parameters to use in the simulation
#'@param n_cores how many cores to use
#'@param endyr year to end the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param prg_chng vector of program change values
#'@return out outputs
#'@export
new_OutputsInt <- function(loc,ParMatrix,n_cores=1,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0, prg_chng, ttt_list) {
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    out <- new_OutputsZint(samp_i=1,ParMatrix=ParMatrix,endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3, prg_chng=prg_chng, ttt_list)
  } else {
    out0 <- mclapply(X=1:nrow(ParMatrix),FUN=new_OutputsZint,mc.cores=n_cores,
                     ParMatrix=ParMatrix,endyr=2050,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3,prg_chng=prg_chng, ttt_list)
    out <- array(NA,dim=c(length(out0),100,586))
    for(i in 1:length(out0)) out[i,,] <- as.matrix(out0[[i]])
  }
  if (sum(Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3)==0) intv<-1;
  if(Int1==1) intv<-2;if(Int2==1) intv<-3; if(Int3==1) intv<-4;
  if(Int4==1) intv<-5; if(Int5==1) intv<-6; if(Scen1==1) intv<-7;
  if(Scen2==1) intv<-8;if(Scen3==1) intv<-9;
  # saveRDS(out, file=paste("~MITUS/",loc,"_results_",intv,".rds",sep=""))
  # save(out,file=paste("/Users/nis100/MITUS/",loc,"_results_",intv,".rda",sep=""))
  return(out)
}

################################################################################
