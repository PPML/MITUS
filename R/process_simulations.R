#' THIS FUNCTION INPUTS A TABLE OF PARAMETERS AND RUNS THE TB MODEL
#' AND GENERATES AN ARRAY OF OUTPUTS
#'@name  OutputsZint
#'@param samp_i which row of the parameter matrix to use
#'@param ParMatrix parameter matrix to use in the simulation
#'@param loc two-digit abbreviation for location
#'@param startyr year to start the simulation
#'@param endyr year to end the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param Scen4 boolean for scenario 4
#'@param Scen5 boolean for scenario 3
#'@param Scen6 boolean for scenario 4
#'@param prg_chng vector of program change values
#'@param ttt_list list of targeted testing and treatment values
#'@return results data frame of output
#'@export
OutputsZint <-  function(samp_i=1,ParMatrix,loc, startyr=1950, endyr=2050,
                                  Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0,Scen4=0,Scen5=0,Scen6=0,
                                  prg_chng=def_prgchng(), ttt_list) {
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par1 <- as.numeric(ParMatrix);
    names(Par1) <- names(ParMatrix)
  } else {  Par1 <- as.numeric(ParMatrix[samp_i,]);
  names(Par1) <- colnames(ParMatrix) }

  P <- Par1

  Int1 <<- Int1;
  Int2 <<- Int2;
  Int3 <<- Int3;
  Int4 <<- Int4;
  Int5 <<- Int5;
  Scen1 <<- Scen1;
  Scen2 <<- Scen2;
  Scen3 <<- Scen3;
  Scen4 <<- Scen4;
  Scen5 <<- Scen5;
  Scen6 <<- Scen6;



  prms <- list()
  prms <- param_init(P,loc,Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3,Scen4,Scen5,Scen6,prg_chng,ttt_list)
  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  m <- cSim( nYrs       = endyr-(startyr-1)         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                      Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
                      rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                      vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
                      HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
                      net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]], SpImmNon    = prms[["SpImmNon"]]        ,
                      mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], RRcrAG = prms[["RRcrAG"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                      TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                      LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
                      ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]],
                      RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                      EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]], NixTb = prms[["NixTb"]],   trans_mat_tot_ages = trans_mat_tot_ages)$Outputs
  colnames(m) <- func_ResNam();
  results<<-as.matrix(m)

  return(results)
}


#'wrapper function for the above function
#'@name  OutputsInt
#'@param loc two-digit abbreviation for location
#'@param ParMatrix parameters to use in the simulation
#'@param n_cores how many cores to use
#'@param startyr year to end the simulation
#'@param endyr year to end the simulation
#'@param startyr year to start the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param prg_chng vector of program change values
#'@param ttt_list list of targeted testing and treatment values
#'@return out outputs
#'@export
OutputsInt <- function(loc,ParMatrix,n_cores=1,startyr=1950,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0,Scen4=0,Scen5=0,Scen6=0,prg_chng, ttt_list) {
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    out <- OutputsZint(samp_i=1,ParMatrix=ParMatrix,loc=loc,endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3,Scen4=Scen4,Scen5=Scen5,Scen6=Scen6, prg_chng=prg_chng, ttt_list=ttt_list)
  } else {
    out0 <- mclapply(X=1:nrow(ParMatrix),FUN=OutputsZint,mc.cores=n_cores,
                     ParMatrix=ParMatrix, loc=loc,endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3,Scen4=Scen4,Scen5=Scen5,Scen6=Scen6,prg_chng=prg_chng,ttt_list= ttt_list)
    out <- array(NA,dim=c(length(out0),endyr-(startyr-1),length(func_ResNam())))

    for(i in 1:length(out0)){
      out[i,,] <- as.matrix(out0[[i]])
    }
    dimnames(out)[[3]]<-func_ResNam()
  }
  if (sum(Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3)==0) intv<-1;
  if(Int1==1) intv<-2;if(Int2==1) intv<-3; if(Int3==1) intv<-4;
  if(Int4==1) intv<-5; if(Int5==1) intv<-6; if(Scen1==1) intv<-7;
  if(Scen2==1) intv<-8;if(Scen3==1) intv<-9;
  # saveRDS(out, file=paste("~MITUS/",loc,"_results_",intv,".rds",sep=""))
  # save(out,file=paste("/Users/nis100/MITUS/",loc,"_results_",intv,".rda",sep=""))
  return(out)
}
