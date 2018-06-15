#'This script creates a function that loops over the log-likelihood
#'functions found in calib_functions.R and updates the
#'load the necessary libraries
library(mnormt)
library(parallel)
library(lhs)

#'@name llikelihoodZ
#'@param samp_i
#'@param ParMatrix ParInit
#'@return lLik
#'
llikelihoodZ <-  function(samp_i,ParMatrix) {
  Par <- ParMatrix[samp_i,]
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
    source("R/param.R")
    IP <- param(P)
    zz <- cSim(  nYrs     =   2050-1950         , nRes      = length(IP[["ResNam"]]), rDxt     = IP[["rDxt"]]    , TxQualt    = IP[["TxQualt"]]   , InitPop  = IP[["InitPop"]]    ,
                 Mpfast     = IP[["Mpfast"]]    , ExogInf   = IP[["ExogInf"]]       , MpfastPI = IP[["MpfastPI"]], Mrslow     = IP[["Mrslow"]]    , rrSlowFB = IP[["rrSlowFB"]]    ,
                 rfast      = IP[["rfast"]]     , RRcurDef  = IP[["RRcurDef"]]      , rSlfCur  = IP[["rSlfCur"]] , p_HR       = IP[["p_HR"]]      , dist_gen = IP[["dist_gen"]]    ,
                 vTMort     = IP[["vTMort"]]    , RRmuRF    = IP[["RRmuRF"]]        , RRmuHR   = IP[["RRmuHR"]]  , muTbRF     = IP[["muTbRF"]]    , Birthst  = IP[["Birthst"]]    ,
                 HrEntEx    = IP[["HrEntEx"]]   , ImmNon    = IP[["ImmNon"]]        , ImmLat   = IP[["ImmLat" ]] , ImmAct     = IP[["ImmAct"]]    , ImmFst   = IP[["ImmFst" ]]    ,
                 mubt       = IP[["mubt"]]      , RelInf    = IP[["RelInf"]]        , RelInfRg = IP[["RelInfRg"]], Vmix       = IP[["Vmix"]]      , rEmmigFB = IP [["rEmmigFB"]]  ,
                 TxVec      = IP[["TxVec"]]     , TunTxMort = IP[["TunTxMort"]]     , rDeft    = IP[["rDeft"]]   , pReTx      = IP[["pReTx"]]     , LtTxPar  = IP[["LtTxPar"]]    ,
                 LtDxPar    = IP[["LtDxPar"]]   , rLtScrt   = IP[["rLtScrt"]]       , RRdxAge  = IP[["RRdxAge"]] , rRecov     = IP[["rRecov"]]    , pImmScen = IP[["pImmScen"]]   ,
                 EarlyTrend = IP[["EarlyTrend"]], NixTrans  = IP[["NixTrans"]]      , can_go   = IP[["can_go"]]  , dist_goal  = IP[["dist_goal"]] , diff_i_v = IP[["diff_i_v"]]   ,
                 dist_orig_v=IP[["dist_orig_v"]])
#'if any output is missing or negative or if any model state population is negative
#'set the likelihood to a hugely negative number (penalized)
    if(sum(is.na(zz$Outputs[65,]))>0 | min(zz$Outputs[65,])<0 | min(zz$V1)<0 ) {
      lLik <- -10^12
  } else {
      M <- zz$Outputs
      colnames(M) <- ResNam
      lLik <- 0
      #' TOTAL DIAGNOSED CASES 1953-2014 - index is same
      v1   <- M[4:66,"NOTIF_ALL"]+M[4:66,"NOTIF_MORT_ALL"]
      addlik <- notif_tot_lik(V=v1); addlik
      lLik <- lLik + addlik
      #' US CASES AGE DISTRIBUTION 1993-2013 - index updated
      v2a   <- M[44:65,205:215]+M[44:65,216:226]
      addlik <- notif_age_us_lLik(V=v2a); addlik
      lLik <- lLik + addlik
      #' FB CASES AGE DISTRIBUTION 1993-2013 - index updated
      v2b   <- (M[44:65,136:146]+M[44:65,189:199]) - (M[44:65,205:215]+M[44:65,216:226])
      addlik <- notif_age_fb_lLik(V=v2b); addlik
      lLik <- lLik + addlik
      #' CASES FB DISTRIBUTION 1993-2014 - index updated
      v3   <- cbind(M[44:66,148]+M[44:66,149]+(M[44:66,201]+M[44:66,202]),
                    M[44:66,147]+M[44:66,200])
      addlik <- notif_fb_lLik(V=v3); addlik
      lLik <- lLik + addlik
      #' CASES FB, US 2010-2014  SLOPE - index updated
      v3   <- cbind(M[44:66,148]+M[44:66,149]+(M[44:66,201]+M[44:66,202]),
                    M[44:66,147]+M[44:66,200])
      addlik <- notif_fbus_slp_lLik(V=v3); addlik
      lLik <- lLik + addlik
      #' CASES HR DISTRIBUTION 1993-2014 - index updated
      v5b   <- cbind(M[44:65,151],M[44:65,150]) + cbind(M[44:65,204],M[44:65,203])
      addlik <- notif_us_hr_lLik(V=v5b); addlik
      lLik <- lLik + addlik
      #' CASES FB RECENT ENTRY DISTRIBUTION 1993-2014 index updated
      v6   <- M[44:65,148:149]+M[44:65,201:202]
      addlik <- notif_fb_rec_lLik(V=v6); addlik
      lLik <- lLik + addlik
      #' TREATMENT OUTCOMES 1993-2012 - index updated
      v11  <- M[44:63,132:134]
      addlik <- tx_outcomes_lLik(V=v11); addlik
      lLik <- lLik + addlik
      #' TOTAL LTBI TREATMENT INITS 2002 - index updated
      v12  <- M[53,152]
      addlik <- tltbi_tot_lLik(V=v12); addlik
      lLik <- lLik + addlik
      #' DIST LTBI TREATMENT INITS 2002 - index updated
      v13  <- M[53,153:155]/M[53,152]
      addlik <- tltbi_dist_lLik(V=v13)*2; addlik
      lLik <- lLik + addlik
      #' LTBI PREVALENCE BY AGE 2011, US - index updated
      v15  <- cbind(M[62,55:65],M[62,33:43]-M[62,55:65])
      v15a <- outer(v15[,1],c(SensLt,1-SensLt))+outer(v15[,2],c(1-SpecLt,SpecLt))
      addlik <- ltbi_us_11_lLik(V=v15a)*2; addlik
      lLik <- lLik + addlik
      #' LTBI PREVALENCE BY AGE 2011, FB - index updated
      v16  <- cbind(M[62,66:76],M[62,44:54]-M[62,66:76])
      v16a <- outer(v16[,1],c(SensLt,1-SensLt))+outer(v16[,2],c(1-SpecLt,SpecLt))
      addlik <- ltbi_fb_11_lLik(V=v16a)*2; addlik
      lLik <- lLik + addlik
      #' TOTAL POP EACH DECADE, BY US/FB - index updated (maybe)
      v17  <- M[,30]+(M[,31]+M[,32])
      addlik <- tot_pop_yr_fb_lLik(V=v17); addlik
      lLik <- lLik + addlik
      #' TOTAL POP AGE DISTRIBUTION 2014 index updated
      v18  <- cbind(M[65,33:43],M[65,44:54])
      addlik <- tot_pop14_ag_fb_lLik(V=v18); addlik
      lLik <- lLik + addlik
      #' TOTAL DEATHS WITH TB 1999-2014 - index updated
      v19  <- M[50:65,227:237]
      addlik <- tb_dth_tot_lLik(V=v19); addlik
      lLik <- lLik + addlik
      #' TB DEATHS 1999-2014 BY AGE - index updated above
      addlik <- tb_dth_age_lLik(V=v19); addlik
      lLik <- lLik + addlik

      #' HOMELESS POP 2010 - index updated
      v23b  <- M[61,29]
      addlik <- homeless_10_lLik(V=v23b); addlik
      lLik <- lLik + addlik

      #' LIKELIHOOD FOR BORGDORFF, FEREBEE & SUTHERLAND ESTIMATES
      v2456  <- c(pfast,pimmed,rslow,rfast,rRecov)
      addlik <- borgdorff_lLik( Par=v2456); addlik
      lLik <- lLik + addlik
      addlik <- ferebee_lLik(   Par=v2456); addlik
      lLik <- lLik + addlik
      addlik <- sutherland_lLik(Par=v2456); addlik
      lLik <- lLik + addlik


      ### ### ### FB RT LIKELIHOOD
      #   v30  <-  (M[66,212]/M[66,195])
      #   addlik <- dnorm(v30,0.075,0.0125,log=T)-dnorm(0.075,0.075,0.0125,log=T); addlik
      #   lLik <- lLik + addlik
      ### ### ###  ALL LIKELIHOODS DONE !!

    } }, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

  return((lLik))  }
#'Local parallelization via multicore
#'@param ParMatrix
#'@param n_cores
#'@return lLik
llikelihood <- function(ParMatrix,n_cores=1) {
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    lLik <- llikelihoodZ(1,t(as.data.frame(ParMatrix)))
  } else {
    lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ,ParMatrix=ParMatrix,mc.cores=n_cores))
  }
return((lLik))
}

#'@name lprior
#'@param ParMatrix
#'@return lPri

lprior <- function(ParMatrix) { # Par = ParInit
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    ParMatrix <- t(as.data.frame(ParMatrix)) }
  lPri <- rep(0,nrow(ParMatrix))
  for(samp_i in 1:nrow(ParMatrix)) {
    # norm2unif
    Par <- ParMatrix[samp_i,]
    Par2 <- pnorm(Par,0,1)
    # unif2true
    Par3 <- Par2
    Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
    Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
    Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
    lPri[samp_i] <- lPrior2(Par,Par3)
    if(is.na(lPri[samp_i]))         {
      lPri[samp_i] <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
    if(lPri[samp_i]%in%c(-Inf,Inf)) {
      lPri[samp_i] <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  } # end inner loop
  return(lPri)
}

#'sample the prior function
#'@name sample.prior1
#'@param n
#'@return random sample of the prior
sample.prior1 <- function(n) {
  rmnorm(n,rep(0,sum(ParamInit$Calib==1)),diag(sum(ParamInit$Calib==1))) }
sample.prior2 <- function(n) {
  qnorm(randomLHS(n,sum(ParamInit$Calib==1)),0,1)   }
sample.prior  <- sample.prior2
