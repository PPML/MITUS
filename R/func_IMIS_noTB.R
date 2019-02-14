#'This scrprmst creates a function that loops over the log-likelihood
#'functions found in calib_functions.R and updates the
#'this is the function that goes into the optimizer
#'@name llikelihoodZ_noTB
#'@param samp_i sample id
#'@param opt_mat matrix of parameters  # Par = par_1
#'@return lLik
llikelihoodZ_noTB <-  function(samp_i,opt_mat) {

  if(min(dim(as.data.frame(opt_mat)))==1) {
    Par <- as.numeric(opt_mat);
    names(Par) <- names(opt_mat)
  } else {  Par <- as.numeric(opt_mat[samp_i,]);
  names(Par) <- colnames(opt_mat) }  ##previously, the distribution of parameters were transformed to normal distribution in
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
    #'format P

    prms <-list()
    prms <- param(P)
    IP <- list()
    IP <- param_init(P)
    trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, prms$adj_fact)

    zz <- cSim_noTB(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                 Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
                 rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                 vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  ,  Birthst  = prms[["Birthst"]]    ,
                 HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
                 net_mig_usb = prms[["net_mig_usb"]]      , net_mig_nusb    = prms[["net_mig_nusb"]]        ,
                 mubt       = prms[["mubt"]]      , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                 TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                 LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                 EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    #'if any output is missing or negative or if any model state population is negative
    #'set the likelihood to a hugely negative number (penalized)
    if(sum(is.na(zz$Outputs[65,]))>0 | min(zz$Outputs[65,])<0 | min(zz$V1)<0 ) {
      lLik <- -10^12
    } else {

      M <- zz$Outputs
      colnames(M) <- prms[["ResNam"]]
      lLik <- 0

      #' TOTAL POP EACH DECADE, BY US/FB - index updated (maybe)
      v17  <- M[,31]+M[,32]
      addlik <- tot_pop_yr_fb_lLik(V=v17); addlik
      lLik <- lLik + addlik
      #' TOTAL POP AGE DISTRIBUTION 2016 index updated
      v18  <- cbind(M[67,33:43],M[67,44:54])
      addlik <- tot_pop16_ag_fb_lLik(V=v18); addlik
      lLik <- lLik + addlik
      # Total DEATHS 1979-2016
      v20a  <- rowSums(M[30:67,121:131])
      addlik <- US_dth_tot_lLik(V=v20a); addlik
      lLik <- lLik + addlik
      #' Total DEATHS 1999-2016 BY AGE
      v20b  <- M[66:67,121:131]
      addlik <- tot_dth_age_lLik(V=v20b); addlik
      lLik <- lLik + addlik
      #' #' Mort_dist 2016
      v21a<- v21  <- M[66:67,521:564]
      for (i in 1:11){
        denom<-M[66:67,2+i]
        for (j in 1:ncol(v21)){
          v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
        } }
      addlik <- mort_dist_lLik(V=v21a); addlik
      lLik <- lLik + addlik
      #' HOMELESS POP 2010 - index updated
      v23b  <- M[61,29]
      addlik <- homeless_10_lLik(V=v23b); addlik
      lLik <- lLik + addlik


    } }, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

  return((lLik))  }

#'This scrprmst creates a function that loops over the log-likelihood
#'functions found in calib_functions.R and updates the
#'this is the function that goes into the optimizer
#'@name llikelihoodZ_noTB_st
#'@param samp_i sample id
#'@param opt_mat matrix of parameters  # Par = par_1
#'@return lLik
llikelihoodZ_noTB_st <-  function(samp_i,opt_mat,loc) {
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)

  if(min(dim(as.data.frame(opt_mat)))==1) {
    Par <- as.numeric(opt_mat);
    names(Par) <- names(opt_mat)
  } else {  Par <- as.numeric(opt_mat[samp_i,]);
  names(Par) <- colnames(opt_mat) }  ##previously, the distribution of parameters were transformed to normal distribution in
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
    #'format P

    prms <-list()
    prms <- param(P)
    IP <- list()
    IP <- param_init(P)
    trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, prms$adj_fact)

    zz <- cSim_noTB(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                      Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
                      rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                      vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  ,  Birthst  = prms[["Birthst"]]    ,
                      HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
                      net_mig_usb = prms[["net_mig_usb"]]      , net_mig_nusb    = prms[["net_mig_nusb"]]        ,
                      mubt       = prms[["mubt"]]      , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                      TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                      LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                      EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    #'if any output is missing or negative or if any model state population is negative
    #'set the likelihood to a hugely negative number (penalized)
    if(sum(is.na(zz$Outputs[65,]))>0 | min(zz$Outputs[65,])<0 | min(zz$V1)<0 ) {
      lLik <- -10^12
    } else {

      M <- zz$Outputs
      colnames(M) <- prms[["ResNam"]]
      lLik <- 0
      st<-which(StateID$USPS==loc)


      ### ### ### TOTAL FB POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
      v17  <- M[,31]+M[,32]
      addlik <- tot_pop_yr_fb_lLik_st(V=v17,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### TOTAL US POP EACH DECADE, FOR US  ### ### ### ### ### ###  D ResNam
      v17b  <- M[,30]
      addlik <- tot_pop_yr_us_lLik_st_00_10(V=v17b,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ### D
      v18  <- cbind(M[65,33:43],M[65,44:54])
      addlik <- tot_pop14_ag_fb_lLik_st(V=v18,st=st); addlik
      lLik <- lLik + addlik

      #' Total DEATHS 1979-2016
      # v20a  <- rowSums(M[30:67,121:131])
      # v20a<-v20a*1e6
      # addlik <- dth_tot_lLik_st(V=v20a,st=st); addlik
      # lLik <- lLik + addlik
      #' Total DEATHS by Decade
      v20a  <- rowSums(M[1:67,121:131])
      addlik <- US_dth_10_tot_lLik(V=v20a); addlik
      lLik <- lLik + addlik
      #' #' Total DEATHS 2015-2016 BY AGE
      v20b  <- M[66:67,121:131]
      v20b  <- v20b/rowSums(v20b)
      #V is scaled in the calib script
      addlik <- tot_dth_age_lLik_st(V=v20b); addlik
      lLik <- lLik + addlik
      #' #' Mort_dist 2016
      v21a<- v21  <- M[51:67,521:564]
      for (i in 1:11){
        denom<-M[51:67,2+i]
        for (j in 1:ncol(v21)){
          v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
        } }
      # v21a<-v21a[,c(1:4,41:44)]
      addlik <- mort_dist_lLik_st(V=v21a); addlik
      lLik <- lLik + addlik


    } }, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

  return((lLik))  }

#'Local parallelization via multicore
#'@name llikelihood_noTB
#'@param start_mat matrix of parameters
#'@param loc two digit postal code of location
#'@param n_cores number of cores to use on the cluster
#'@return lLik
#'@export

llikelihood_noTB <- function(start_mat,loc="US",n_cores=1) {
  if(loc=="US"){
    if(dim(as.data.frame(start_mat))[2]==1) {
      lLik <- llikelihoodZ_noTB(1,t(as.data.frame(start_mat)))
    } else {
      lLik <- unlist(mclapply(1:nrow(start_mat),llikelihoodZ_noTB,opt_mat=start_mat,mc.cores=n_cores))
    }
  } else {
  if(dim(as.data.frame(start_mat))[2]==1) {
    lLik <- llikelihoodZ_noTB_st(1,t(as.data.frame(start_mat)),loc=loc)
    } else {
      lLik <- unlist(mclapply(1:nrow(start_mat),llikelihoodZ_noTB_st,opt_mat=start_mat,loc=loc,mc.cores=n_cores))
    }}
  return((lLik)) }


