#'This scrprmst creates a function that loops over the log-likelihood
#'functions found in calib_functions.R and updates the
#'this is the function that goes into the optimizer
#'load the necessary libraries
library(mnormt)
library(parallel)
library(lhs)
#source("R/define_P.R")
#'@name llikelihoodZ_noRB
#'@param samp_i sample id
#'@param ParMatrix matrix of parameters  # Par = par_1
#'@return lLik
llikelihoodZ_noRB <-  function(samp_i,ParMatrix) {
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
    prms <-list()
    prms <- param(P)
    IP <- list()
    IP <- param_init(P)
    tm<-matrix(0,16,16)
    diag(tm)<-1
    trans_mat_tot_ages<<-matrix(tm,16,176)

    zz <- cSim(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                 Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
                 rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = matrix(c(1,rep(0,15)),4,4)    ,
                 vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
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
      M <- zz$Outputs
      colnames(M) <- prms[["ResNam"]]
      lLik <- 0
      #' TOTAL DIAGNOSED CASES 1953-2016 - index is same
      v1   <- M[4:67,"NOTIF_ALL"]+M[4:67,"NOTIF_MORT_ALL"]
      addlik <- notif_tot_lik(V=v1); addlik
      lLik <- lLik + addlik
      #' US CASES AGE DISTRIBUTION 1993-2016 - index updated
      v2   <- M[44:67,205:215]+M[44:67,216:226]
      v2a <- v2[,-11]; v2a[,10] <- v2a[,10]+v2[,11]
      addlik <- notif_age_us_lLik(V=v2a); addlik
      lLik <- lLik + addlik
      #' FB CASES AGE DISTRIBUTION 1993-2016 - index updated
      v2   <- (M[44:67,136:146]+M[44:67,189:199]) - (M[44:67,205:215]+M[44:67,216:226])
      v2b <- v2[,-11]; v2b[,10] <- v2b[,10]+v2[,11]
      addlik <- notif_age_fb_lLik(V=v2b); addlik
      lLik <- lLik + addlik
      #' CASES FB DISTRIBUTION 1993-2016 - index updated
      v3   <- cbind(M[44:67,148]+M[44:67,149]+(M[44:67,201]+M[44:67,202]),
                    M[44:67,147]+M[44:67,200])
      addlik <- notif_fb_lLik(V=v3); addlik
      lLik <- lLik + addlik
      #' CASES FB, US 2010-2014  SLOPE - index updated
      v3   <- cbind(M[44:66,148]+M[44:66,149]+(M[44:66,201]+M[44:66,202]),
                    M[44:66,147]+M[44:66,200])
      addlik <- notif_fbus_slp_lLik(V=v3); addlik
      lLik <- lLik + addlik
      #' CASES HR DISTRIBUTION 1993-2016 - index updated
      v5b   <- cbind(M[44:67,151],M[44:67,150]) + cbind(M[44:67,204],M[44:67,203])
      addlik <- notif_us_hr_lLik(V=v5b); addlik
      lLik <- lLik + addlik
      #' CASES FB RECENT ENTRY DISTRIBUTION 1993-2014 index updated
      v6   <- M[44:65,148:149]+M[44:65,201:202]
      addlik <- notif_fb_rec_lLik(V=v6); addlik
      lLik <- lLik + addlik
      #' TREATMENT OUTCOMES 1993-2012 - index updated
      v11  <- M[44:65,132:134]
      addlik <- tx_outcomes_lLik(V=v11); addlik
      lLik <- lLik + addlik
      #' TOTAL LTBI TREATMENT INITS 2002 - index updated
      v12  <- M[53,152]
      addlik <- tltbi_tot_lLik(V=v12); addlik
      lLik <- lLik + addlik
      #' DIST LTBI TREATMENT INITS 2002 - index updated
      v13  <- M[53,153:154]/M[53,152]
      addlik <- (tltbi_dist_lLik(V=v13))*2; addlik
      lLik <- lLik + addlik
      #' LTBI PREVALENCE BY AGE 2011, US - index updated
      v15  <- cbind(M[62,55:65],M[62,33:43]-M[62,55:65])
      v15a <- outer(v15[,1],c(P[["SensLt"]],1-P[["SensLt"]]))+outer(v15[,2],c(1-P[["SpecLt"]],P[["SpecLt"]]))
      addlik <- ltbi_us_11_lLik(V=v15a)*2; addlik
      lLik <- lLik + addlik
      #' LTBI PREVALENCE BY AGE 2011, FB - index updated
      v16  <- cbind(M[62,66:76],M[62,44:54]-M[62,66:76])
      v16a <- outer(v16[,1],c(P[["SensLt"]],1-P[["SensLt"]]))+outer(v16[,2],c(1-P[["SpecLt"]],P[["SpecLt"]]))
      addlik <- ltbi_fb_11_lLik(V=v16a)*2; addlik
      lLik <- lLik + addlik
      #' TOTAL POP EACH DECADE, BY US/FB - index updated (maybe)
       v17  <- M[,31]+M[,32]
       addlik <- tot_pop_yr_fb_lLik(V=v17); addlik
       lLik <- lLik + addlik
      #' TOTAL POP AGE DISTRIBUTION 2016 index updated
      v18  <- cbind(M[67,33:43],M[67,44:54])
      addlik <- tot_pop16_ag_fb_lLik(V=v18); addlik
      lLik <- lLik + addlik
      #' TOTAL DEATHS WITH TB 1999-2014 - index updated
      v19  <- M[50:65,227:237]
      addlik <- tb_dth_tot_lLik(V=rowSums(v19)); addlik
      lLik <- lLik + addlik
      #' TB DEATHS 1999-2014 BY AGE - index updated above
      addlik <- tb_dth_age_lLik(V=v19); addlik
      lLik <- lLik + addlik
      #' Total DEATHS 1979-2016
      v20a  <- rowSums(M[30:67,121:131])
      addlik <- US_dth_tot_lLik(V=v20a); addlik
      lLik <- lLik + addlik
      #' Total DEATHS 1999-2016 BY AGE
      v20b  <- M[50:67,121:131]
      addlik <- tot_dth_age_lLik(V=v20b); addlik
      lLik <- lLik + addlik

      #' HOMELESS POP 2010 - index updated
      v23b  <- M[61,29]
      addlik <- homeless_10_lLik(V=v23b); addlik
      lLik <- lLik + addlik

      #' LIKELIHOOD FOR BORGDORFF, FEREBEE & SUTHERLAND ESTIMATES

      v2456  <- list(prms[["Mpfast"]],prms[["Mrslow"]], prms[["rfast"]],prms[["rRecov"]])
      addlik <- borgdorff_lLik( Par=v2456); addlik
      lLik <- lLik + addlik
      addlik <- ferebee_lLik(Par=v2456); addlik
      lLik <- lLik + addlik
      addlik <- sutherland_lLik(Par=v2456); addlik
      lLik <- lLik + addlik

      # ### ### ### LIKELIHOOD FOR TIEMERSMA ESTS ### ### ### ### ### ### ~~~
      v35   <- c(P["rSlfCur"],P["muIp"])
      addlik <- tiemersma_lLik(Par=v35); addlik
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
#'@name llikelihood_noRB
#'@param ParMatrix matrix of parameters
#'@param n_cores number of cores to use on the cluster
#'@return lLik
#'@export
llikelihood_noRB <- function(ParMatrix,n_cores=1) {
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    lLik <- llikelihoodZ_noRB(1,t(as.data.frame(ParMatrix)))
  } else {
    lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ_noRB,ParMatrix=ParMatrix,mc.cores=n_cores))
  }
return((lLik))
}



