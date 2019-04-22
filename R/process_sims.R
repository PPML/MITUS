#' THIS FUNCTION INPUTS A TABLE OF PARAMETERS AND RUNS THE TB MODEL
#' AND GENERATES AN ARRAY OF OUTPUTS


#load("~/MITUS/data/parAll200_9-14-16.rData")
# load("data/ParamInit_2018.rData")
#'function to run the model
#'@name OutputsZint
#'@param samp_i how many samples
#'@param ParMatrix parameters to use in the simulation
#'@param loc location
#'@param endyr year to end the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@return results data frame of output
#'@export
OutputsZint <-  function(samp_i=1,ParMatrix,loc="US",startyr=1950, endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0) {
   if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par <- as.numeric(ParMatrix);
    names(Par) <- names(ParMatrix)
} else {  Par <- as.numeric(ParMatrix[samp_i,]);
    names(Par) <- colnames(ParMatrix) }
    Int1 <<- Int1;
    Int2 <<- Int2;
    Int3 <<- Int3;
    Int4 <<- Int4;
    Int5 <<- Int5;
    Scen1 <<- Scen1;
    Scen2 <<- Scen2;
    Scen3 <<- Scen3

   P <- Par
   prms <-list()
   prms <- param(P)
   IP <- list()
   IP <- param_init(P)
   NY<-endyr-startyr
   trans_mat_tot_ages<<-reblncd(mubt = IP$mubt,can_go = can_go,RRmuHR = IP$RRmuHR[2], RRmuRF = IP$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v,adj_fact= IP[["adj_fact"]])
   if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
   m <- cSim( nYrs       = NY         , nRes      = length(prms[["ResNam"]])  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
               Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
               rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
               vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
               HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
               mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
               TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
               LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
               EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)$Outputs
   colnames(m) <- IP[["ResNam"]];
   results<<-as.data.frame(m)

   return(results)
}


#'wrapper function for the above function
#'@name OutputsInt
#'@param ParMatrix parameters to use in the simulation
#'@param loc location of simulation
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
#'@return out outputs
#'@export
OutputsInt <- function(ParMatrix, loc="US",n_cores=1,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0) {
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    out <- OutputsZint(samp_i=1,ParMatrix=ParMatrix,endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3)
} else {
    out0 <- mclapply(X=1:nrow(ParMatrix),FUN=OutputsZint,mc.cores=n_cores,
                       ParMatrix=ParMatrix,endyr=2050,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3)
    out <- array(NA,dim=c(length(out0),100,564))
    for(i in 1:length(out0)) out[i,,] <- as.matrix(out0[[i]])
}
  if (sum(Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3)==0) intv<-1;
  if(Int1==1) intv<-2;if(Int2==1) intv<-3;if(Int3==1) intv<-4;
  if(Int4==1) intv<-5; if(Int5==1) intv<-6; if(Scen1==1) intv<-7;
  if(Scen2==1) intv<-8;if(Scen3==1) intv<-9;
  saveRDS(out,file=paste0("/Users/nis100/MITUS/inst/",loc,"/results_",intv,"_",Sys.Date(),".rds"))
  return(out)
}

################################################################################
