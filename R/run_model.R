library(Rcpp)

source("R/param.R")
source("R/add_params.R")
sourceCpp("src/tb_model.cpp")


m <- cSim( nYrs     =   2100-1950, nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
      Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
      rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist = dist            ,
      vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
      HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
      mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
      TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
      LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
      EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLtX    = EffLtX   , dLtt     = dLtt        , NixTrans = NixTrans )$V0




source("R/process_sims.R")


OutputsInt(P)
#ParMatrix,n_cores=1,endyr=2100,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0)

OutputsInt(P, 1, 1990, 0, 0, 0, 0, 0, 0,0,0 )
