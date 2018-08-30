#'This scrprmst creates a function that loops over the log-likelihood
#'functions found in calib_functions.R and updates the
#'this is the function that goes into the optimizer
#'load the necessary libraries
library(mnormt)
library(parallel)
library(lhs)
#source("R/define_P.R")
#'@name llikelihoodZ
#'@param samp_i sample id
#'@param ParMatrix matrix of parameters  # Par = par_1
#'@return lLik
llikelihoodZ_noTB <-  function(samp_i,ParMatrix) {
  #'load the necessary calibration data
  data("CalibDat_2018-07-12", package='MITUS') # CalibDat
  #'Log-likelihood functions
  #'Assign the calibration importance weights from CalibDat
  #'These weights are based on year of the simulation.
  wts <- CalibDat[["ImptWeights"]]
  #'format P
  data("ParamInitUS_2018-08-06_final", package='MITUS')# ParamInit
  P  <- ParamInit[,1];
  names(P) <- rownames(ParamInit)
  ii <-  ParamInit[,5]==1
  ParamInitZ <- ParamInit[ParamInit$Calib==1,]
  idZ0 <- ParamInitZ[,4]==0
  idZ1 <- ParamInitZ[,4]==1
  idZ2 <- ParamInitZ[,4]==2

  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par <- as.numeric(ParMatrix);
    names(Par) <- names(ParMatrix)
  } else {  Par <- as.numeric(ParMatrix[samp_i,]);
  names(Par) <- colnames(ParMatrix) }  ##previously, the distribution of parameters were transformed to normal distribution in
  ##to facilitate comparisons. These first two steps convert these parameters back to their
  ##distributions
  # normal to uniform
  Par2 <- pnorm(Par,0,1)
  # uniform to true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  P <- P

  jj <- tryCatch({
    data("CalibDat_2018-07-12", package='MITUS') # CalibDat
    #'Log-likelihood functions
    #'Assign the calibration importance weights from CalibDat
    #'These weights are based on year of the simulation.
    wts <- CalibDat[["ImptWeights"]]
    #'format P

    prms <-list()
    prms <- param(P)
    IP <- list()
    IP <- param_init(P)
    trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, prms$adj_fact)

    zz <- cSim_noTB(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                 Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
                 rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                 vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , muTbRF     = prms[["muTbRF"]]    , Birthst  = prms[["Birthst"]]    ,
                 HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
                 mubt       = prms[["mubt"]]      , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                 TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                 LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                 EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    #'if any output is missing or negative or if any model state population is negative
    #'set the likelihood to a hugely negative number (penalized)
    if(sum(is.na(zz$Outputs[65,]))>0 | min(zz$Outputs[65,])<0 | min(zz$V1)<0 ) {
      lLik <- -10^12
    } else {
      data("ParamInitUS_2018-08-06_final", package='MITUS')# ParamInit
      P  <- ParamInit[,1];
      names(P) <- rownames(ParamInit)
      ii <-  ParamInit[,5]==1
      ParamInitZ <- ParamInit[ParamInit$Calib==1,]
      idZ0 <- ParamInitZ[,4]==0
      idZ1 <- ParamInitZ[,4]==1
      idZ2 <- ParamInitZ[,4]==2
      M <- zz$Outputs
      colnames(M) <- prms[["ResNam"]]
      lLik <- 0

      #' TOTAL POP EACH DECADE, BY US/FB - index updated (maybe)
      v17  <- M[,30]+(M[,31]+M[,32])
      addlik <- tot_pop_yr_fb_lLik(V=v17); addlik
      lLik <- lLik + addlik
      #' TOTAL POP AGE DISTRIBUTION 2016 index updated
      v18  <- cbind(M[67,33:43],M[67,44:54])
      addlik <- tot_pop16_ag_fb_lLik(V=v18); addlik
      lLik <- lLik + addlik

      #' Total DEATHS 1999-2016 BY AGE
      v20  <- M[50:67,255:265]+M[50:67,266:276]
      addlik <- tot_dth_age_lLik(V=v20); addlik
      lLik <- lLik + addlik
      #' HOMELESS POP 2010 - index updated
      v23b  <- M[61,29]
      addlik <- homeless_10_lLik(V=v23b); addlik
      lLik <- lLik + addlik


    } }, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

  return((lLik))  }



