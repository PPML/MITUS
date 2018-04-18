library(Rcpp)
source("R/process_sims.R")
sourceCpp("src/tb_model.cpp")
#source("R/start_model.R")
# source("R/param.R")
# source("R/add_params.R")
 #source("R/param_init.R")
#source("R/gen_reblnc_pop.R")

# source("R/graphs.R")
# source("R/graphs_specific.R")

  OutputsZint(samp_i = 1, parAll200,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0)

save(results,file = paste("data/results", Sys.time(),".rData"))
write.csv(results, file = paste("MITUS_results/results", Sys.time(),".csv"))
source("R/graphs_demo.R")
tb_graph_demo(1950,2050,results)

# tb_graph_orig(results)

# source("R/graphs_specific.R")
# tb_graph_all(1950,2050,results)
# tb_graph_vital(1950,2049, results)
# tb_graph_specific(1993,2005,results, "mort")

##transmat
#system.time(


#
# colnames(m) <-ResNam
#
#
#
# results<-data.frame(m)

# v1<-cSim(      nYrs     =   2100-1950 ,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
#            Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
#            rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist = dist            ,
#            vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
#            HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
#            mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
#            TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
#            LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
#            EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans )$V1
#
# v0<-cSim( nYrs     =   2100-1950,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
#       Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
#       rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist = dist            ,
#       vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
#       HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
#       mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
#       TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
#       LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
#       EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans )$V0

#
g <- cSim( nYrs     =   2050-1950,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
           Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
           rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist_gen = dist_gen    ,
           vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
           HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
           mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
           TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
           LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
           EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans    ,
           can_go   = can_go      , dist_goal=dist_goal  , diff_i_v = diff_i_v,
           dist_orig_v=dist_orig_v )$trans_mat
#
# g
n <- cSim( nYrs     =   2050-1950,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
           Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
           rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist_gen = dist_gen    ,
           vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
           HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
           mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
           TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
           LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
           EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans    ,
           can_go   = can_go      , dist_goal=dist_goal  , diff_i_v = diff_i_v,
           dist_orig_v=dist_orig_v)$dist_i_v
 n

# m[,2]
# m[3,]
# v1
# plot(v1)
# v0
# plot(v0)

