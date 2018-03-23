library(Rcpp)

source("R/param.R")
source("R/add_params.R")
sourceCpp("src/tb_model.cpp")

m <- cSim( nYrs     =   2100-1950,  nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt    , InitPop   = InitPop    ,
           Mpfast   = Mpfast      , ExogInf   = ExogInf      , MpfastPI  = MpfastPI , Mrslow    = Mrslow     , rrSlowFB = rrSlowFB    ,
           rfast    = rfast       , RRcurDef = RRcurDef      , rSlfCur  = rSlfCur   , p_HR     = p_HR        , dist = dist            ,
      vTMort   = vTMort      , vRFMort = vRFMort        , RRmuHR    = RRmuHR   , muTbRF = muTbRF        , Birthst   = Birthst    ,
      HrEntEx   = HrEntEx    , ImmNon    = ImmNon       , ImmLat    = ImmLat   , ImmAct   = ImmAct      , ImmFst    = ImmFst     ,
      mubt     = mubt        , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix       , rEmmigFB  = rEmmigFB   ,
      TxVec    = TxVec       , TunTxMort = TunTxMort    , rDeft     = rDeft    , pReTx     = pReTx      , LtTxPar  = LtTxPar     ,
      LtDxPar  = LtDxPar     , rLtScrt   = rLtScrt      , RRdxAge  = RRdxAge   , rRecov   = rRecov      , pImmScen  = pImmScen   ,
      EarlyTrend = EarlyTrend, EffLt    = EffLt         , EffLt0    = EffLt0   , dLtt     = dLtt        , NixTrans = NixTrans )$Outputs


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
v1
plot(v1)
v0


