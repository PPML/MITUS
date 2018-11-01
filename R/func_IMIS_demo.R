#'This scrprmst creates a function that loops over the log-likelihood
#'functions found in calib_functions.R and updates the
#'this is the function that goes into the optimizer
#'@name llikelihoodZ_demo
#'@param samp_i sample id
#'@param ParMatrix matrix of parameters  # Par = par_1
#'@return lLik
llikelihoodZ_demo <-  function(samp_i,ParMatrix) {
  library(mnormt)
  library(parallel)
  library(lhs)
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

    dd <- cSim_demo(  nYrs       = 2018-1950         ,
                      nRes      = 24,
                      InitPop  = prms[["InitPop"]]    ,
                      Birthst  = prms[["Birthst"]]    ,
                      ImmNon    = prms[["ImmNon"]]    ,
                      ImmLat   = prms[["ImmLat" ]]    ,
                      ImmAct     = prms[["ImmAct"]]   ,
                      ImmFst   = prms[["ImmFst" ]]    ,
                      mubt       = prms[["mubt"]]
    )
    #'if any output is missing or negative or if any model state population is negative
    #'set the likelihood to a hugely negative number (penalized)
    if(sum(is.na(dd$Outputs[65,]))>0 | min(dd$Outputs[65,])<0 | min(dd$V1)<0 ) {
      lLik <- -10^12
    } else {
      # data("ParamInitUS_2018-08-06_final", package='MITUS')# ParamInit
      # P  <- ParamInit[,1];
      # names(P) <- rownames(ParamInit)
      # ii <-  ParamInit[,5]==1
      # ParamInitZ <- ParamInit[ParamInit$Calib==1,]
      # idZ0 <- ParamInitZ[,4]==0
      # idZ1 <- ParamInitZ[,4]==1
      # idZ2 <- ParamInitZ[,4]==2
      D <- dd$Outputs
      colnames(D) <- c(prms$ResNam[1:13],prms$ResNam[121:131])
      lLik <- 0

      #' TOTAL POP EACH DECADE, BY US/FB
      v17  <- D[,2]
      addlik <-US_pop_tot_lLik(V=v17); addlik
      lLik <- lLik + addlik
      #' TOTAL POP AGE DISTRIBUTION 2016 index updated
      v18  <- D[65,3:13]
      addlik <- tot_pop_age_lLik(V=v18); addlik
      lLik <- lLik + addlik

      #' Total DEATHS 1999-2016 BY AGE
      v20  <- D[50:67,14:24]
      addlik <- tot_dth_age_lLik(V=v20); addlik
      lLik <- lLik + addlik
      #' Total DEATHS 1999-2016 BY AGE
      v20b  <- D[50:67,14:24]
      addlik <- tot_dth_age_lLik(V=v20b); addlik
      lLik <- lLik + addlik


    } }, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

  return((lLik))  }

#'Local parallelization via multicore
#'@name llikelihood_demo
#'@param ParMatrix matrix of parameters
#'@param n_cores number of cores to use on the cluster
#'@return lLik
#'@export
llikelihood_demo <- function(ParMatrix,n_cores=1) {
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    lLik <- llikelihoodZ_demo(1,t(as.data.frame(ParMatrix)))
  } else {
    lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ_demo,ParMatrix=ParMatrix,mc.cores=n_cores))
  }
  return((lLik))
}


