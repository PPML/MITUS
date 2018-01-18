

######### TAKES TABLE OF PARAMETERS PRODUCES ARRAY OF OUTPUTS
library(parallel)
########  FUNCTION  ########

  OutputsZint <-  function(samp_i=1,ParMatrix,endyr=2100,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0) {
    if(min(dim(as.data.frame(ParMatrix)))==1) { ; Par <- as.numeric(ParMatrix); names(Par) <- names(ParMatrix)
    } else {  Par <- as.numeric(ParMatrix[samp_i,]); names(Par) <- colnames(ParMatrix) }
    Int1 <<- Int1; Int2 <<- Int2; Int3 <<- Int3; Int4 <<- Int4; Int5 <<- Int5;
    Scen1 <<- Scen1; Scen2 <<- Scen2; Scen3 <<- Scen3
    P <<- Par
    source("ParamUS_V5070int.r")
    M <-         cSim( nYrs     = 2100-1950, nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt   , InitPop   = InitPop, 
                       p_HR     = p_HR     , Mpfast   = Mpfast        , ExogInf   = ExogInf  , MpfastPI  = MpfastPI  , RelInfHivt = RelInfHivt,
                       pimmed   = pimmed   , MpSmPos  = MpSmPos       , HivHrPar  = HivHrPar , RRmuHR    = RRmuHR    , Mrslow    = Mrslow, 
                       rfast    = rfast    , RRcurDef = RRcurDef      , fAR       = fAR      , VrSlfCur  = VrSlfCur  , VrSmConv  = VrSmConv, 
                       vTMort   = vTMort   , vHMort   = vHMort        , Birthst   = Birthst  , ImmNon    = ImmNon    , ImmLat    = ImmLat, ImmFst    = ImmFst,
                       ImmAct   = ImmAct   , DrN      = DrN           , DrE       = DrE      , TxExpAge  = TxExpAge  , p_Imm_SP  = p_Imm_SP, 
                       mubt     = mubt     , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix      , rArtInit  = rArtInit , 
                       vHxtoHy  = vHxtoHy  , rArtDef  = rArtDef       , rHIVt     = rHIVt    , rEmmigFB  = rEmmigFB  , pDstt     = pDstt   , 
                       TxMat    = TxMat    , TunTxMort = TunTxMort    , rDeft     = rDeft    , rDeftH    = rDeftH    , pReTx     = pReTx    ,
                       LtTxPar  = LtTxPar  , LtDxPar  = LtDxPar       , rLtScrt   = rLtScrt  , HrEntEx   = HrEntEx   , muTbH     = muTbH,
                       RRdxAge  = RRdxAge  , rRecov   = rRecov        , pImmScen  = pImmScen, EarlyTrend = EarlyTrend, rrSlowFB = rrSlowFB,
                       EffLt    = EffLt    ,EffLtX    = EffLtX     , dLtt     = dLtt          , NixTrans = NixTrans )$Outputs
    colnames(M) <- ResNam;
    return(M)  }
  
  ######## WRAPPER FUNCTION  ########
  
  OutputsInt <- function(ParMatrix,n_cores=1,endyr=2100,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0) {
    if(min(dim(as.data.frame(ParMatrix)))==1) {
      out <- OutputsZint(samp_i=1,ParMatrix=ParMatrix,endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3)  } else {
      out0 <- mclapply(X=1:nrow(ParMatrix),FUN=OutputsZint,mc.cores=n_cores,
                       ParMatrix=ParMatrix,endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3)
      out <- array(NA,dim=c(dim(out0[[1]]),length(out0)))
      for(i in 1:length(out0)) out[,,i] <- out0[[i]]
      }
    return(out) }

###############################################
 