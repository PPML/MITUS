library(Rcpp)
library(RcppArmadillo)

#source("R/start_model.R")
source("R/param.R")
source("R/add_params.R")
source("R/param_init.R")
source("R/gen_reblnc_pop.R")

sourceCpp("src/tb_model.cpp")
source("R/graphs.R")

#system.time(

m <- cSim( nYrs     =   2050-1950,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
           Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
           rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist_gen = dist_gen    ,
           vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
           HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
           mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
           TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
           LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
           EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans    ,
           can_go   = can_go      , did_go=did_go            , dist_goal=dist_goal, dist_goal_v=dist_goal_v  , diff_i_v = diff_i_v,
           dist_orig_v=dist_orig_v, dist_new = dist_new )$Outputs

colnames(m) <-ResNam

#)

results=data.frame(m)
save(results,file = paste("data/results", Sys.time(),".rData"))
write.csv(results, file = paste("MITUS_results/results", Sys.time(),".csv"))

tb_graph_all(1950,2050,results)
tb_graph_specific(1993,2005,results, "mort")


v1<-cSim(      nYrs     =   2100-1950 ,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
           Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
           rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist = dist            ,
           vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
           HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
           mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
           TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
           LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
           EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans )$V1

v0<-cSim( nYrs     =   2100-1950,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
      Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
      rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist = dist            ,
      vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
      HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
      mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
      TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
      LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
      EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans )$V0

m[,2]
m[3,]
v1
plot(v1)
v0
plot(v0)

