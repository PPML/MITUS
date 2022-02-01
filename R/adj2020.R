llikelihood2020 <- function(samp_i, start_mat, TB=1){

  ### Create the parameter vector
  ### par2020 will hold the three parameter adjustment factors to be fed into the
  ### param_init function for the model run.

  if(min(dim(as.data.frame(start_mat)))==1) {
    par2020 <- as.numeric(start_mat)
  } else{
    par2020 <- as.numeric(start_mat[samp_i,])
  }
  # par2020 <- c(0.088, 0.43, 0.28)
  # par2020<-c(0.2046272,0.4325273,0.2712293)
  # par2020 <- c (0.1797244, 0.4360917, 0.2754278)
  # par2020 <- c(0.3952822,0.4029466,0.2458458)
  # par2020 <- c(0.2049939,0.4322080, 0.2785897)
  # par2020 <- c(0.2046310,0.4325234, 0.2712266)
  # par2020 <- c(0.2859564,0.4206157,0.2653754)
  # par2020 <- c(0.3205440,0.4151213,0.2593764)
  # par2020 <- c(0.4064649,0.4010572,0.2438105)
  # par2020 <- c(0.5003912,0.3857184,0.2265612)
  # par2020 <- c(0.4508752,0.3942389,0.2361109)
  # par2020 <- c(0.3865738,0.4044850,0.2477164)
  # par2020 <- c(0.3638802,0.4045582,0.2426989) #as of 1/29/22
  # par2020 <- c(0.3957942, 0.3992307, 0.2365767) #as of 1/31/22

  names(par2020) <- c("Immig", "Dxt", "Trans")
  # print(colnames(start_mat))
  # print(par2020)
  ### Set P to the optimal parameter set which is loaded in with load_model()
  P <- Par[1,]
  ### Setup the model run
  prg_chng<-def_prgchng(P)
  prms <-list()
  prms <- param_init(P,"US",prg_chng=prg_chng, ttt_list=def_ttt(), immig = par2020["Immig"])
  prms$rDxt[843:864,]<-prms$rDxt[843:864,] - (prms$rDxt[843:864,]*par2020["Dxt"])
  prms$NixTrans[843:864]<- (1-par2020["Trans"])
  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  # jj <- tryCatch({
  zz <- cSim( nYrs       = 2021-1950         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
              Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
              rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
              vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
              HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
              net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
              mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], RRcrAG = prms[["RRcrAG"]],
              Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
              TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
              LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
              ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
              EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)

  if(sum(is.na(zz$Outputs[68,]))>0 | min(zz$Outputs[68,])<0 | min(zz$V1)<0 ) {
    lLik <- -10^12
  } else {
    ### read in the basecase for creating comparison
    bcRes <- readRDS(system.file("US/US_basecase_0719.rds", package="MITUS"))
    M <- zz$Outputs
    colnames(M) <- func_ResNam()
    # saveRDS(M,"~/MITUS/inst/US/US_2020basecase_010621.rds")
    lLik <- 0
    # TOTAL DIAGNOSED CASES 2020
    v1bc   <- bcRes[70,"NOTIF_ALL"]+bcRes[70,"NOTIF_MORT_ALL"]
    v1M     <- M[71,"NOTIF_ALL"]+M[71,"NOTIF_MORT_ALL"]
    v1 <- 1 - (v1M / v1bc)
    addlik <- notif_tot_20_lik(V=v1); addlik
    lLik <- lLik + addlik

    # CASES FB RECENT ENTRY DISTRIBUTION 2020
    v2bc   <- (bcRes[70,148] + bcRes[70,201])/ sum(bcRes[70,201:202],bcRes[70,148:149])
    v2M    <- (M[71,148] + M[71,201])/ sum(M[71,201:202],M[71,148:149])
    v2 <- 1 - (v2M / v2bc)
    addlik <- notif_NUSBrec_20_lik(V=v2); addlik
    lLik <- lLik + addlik

    # CASES RECENT TRANSMISSION DISTRIBUTION 2020
    v3abc  <- sum(bcRes[70,184:185])/sum(bcRes[70,168:169])
    v3aM   <- sum(M[71,172])/sum(M[71,156])
    v3     <- 1- (v3aM / v3abc)
    addlik <- notif_RT_20_lik(V=v3); addlik
    lLik   <- lLik + addlik

  }
  return((lLik))
}

### Use three likelihoods to evaluate the differential changes in our mechanisms:
### 1. Rate of Case Finding
### 2. Immigration Volume
### 3. Contact Rates
### The functions ###
###  Measure the change in TB cases in 2020 from 2019 (basecase)
notif_tot_20_lik <- function(V) {
  ### We are basing this off of the preliminary data that suggests a 20% decrease
  case_diff_tot <-  0.1953791
  adj_1         <- dnorm(case_diff_tot,case_diff_tot,0.01/1.96,log=T)
  dnorm(case_diff_tot,V,0.01/1.96,log=T) - adj_1
}

### Measure the % change in NUSB recent entry in 2020 from 2019 (basecase)
notif_NUSBrec_20_lik <- function(V) {
  ### We are basing this off of the preliminary data that suggests a 6% decrease
  case_diff_NUSB <- 0.174
  adj_2         <- dnorm(case_diff_NUSB,case_diff_NUSB,0.025/1.96,log=T)
  dnorm(case_diff_NUSB,V,0.025/1.96,log=T) - adj_2
}

### Measure the % change in recent transmission cases
notif_RT_20_lik <- function(V) {
  ### We are basing this off of the preliminary data that suggests a 0% decrease
  case_diff_RT <- 0.001
  adj_3         <- dnorm(case_diff_RT,case_diff_RT,0.025/1.96,log=T)
  dnorm(case_diff_RT,V,0.025/1.96,log=T) - adj_3
}

### Create the starting parameter matrix

### First make a matrix of each value and its prior boundaries
# paraminit2020 <- matrix(0,3,5)
# rownames(paraminit2020) <- c("Immig", "Dxt", "Trans")
# paraminit2020[,1:3]<- cbind(c(.17, .5, .5), c(.07, .1, .3), c(.27, .9, .7))
# paraminit2020[,4:5] <- c(paraminit2020[,1],(paraminit2020[,3]-paraminit2020[,2])/3.92)
#
# startval2020 <-randomLHS(10,nrow(paraminit2020))
# colnames(startval2020) <- c("Immig", "Dxt", "Trans")
### Optimizing Function
optim_2020 <- function(df, samp_i=1,n_cores=2, TB=1){

  ### Define the posterior function
  posterior = function(theta) {
    -lprior2020(theta) - llikelihood2020(start_mat=theta)
  }

  ### Edge case for single parameter set
  if(min(dim(as.data.frame(df)))==1) {
    df1 <- as.numeric(df)
    names(df1) <- names(df)
  } else{
    df1 <- as.numeric(df[samp_i,])
    names(df1) <- names(df)
  }

  b<-samp_i
  # for (i in min(b, nrow(StartVal))){
  o1  <- optim(df1, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
  save(o1,file=paste("Opt_US2020_r1_", b,"_", Sys.Date(),".rda",sep=""))
  o2  <- optim(o1$par,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
  save(o2,file=paste("Opt_US2020_r2_", b,"_", Sys.Date(),".rda",sep=""))
  o3  <- optim(o2$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
  save(o3,file=paste("Opt_US2020_r3_", b,"_", Sys.Date(),".rda",sep=""))
}

lprior2020 <- function(ParMatrix = startval2020, InitPar = paraminit2020) { # Par = ParInit
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    ParMatrix <- t(as.data.frame(ParMatrix)) }
  lPri <- rep(0,nrow(ParMatrix))
  for(samp_i in 1:nrow(ParMatrix)) {
    # norm2unif
    Par <- ParMatrix[samp_i,]
    # unif2true
    Par3 <- Par
    Par3 <- qnorm(Par, mean    = InitPar[,4], sd = InitPar[,5])

    if(dim(as.matrix(Par))[2]==1) Par <- t(as.matrix(Par))
    ldensity <- dmnorm(Par,rep(0,nrow(InitPar)),diag(nrow(InitPar)),log=T)
    lDensTrue <- rep(NA,nrow(InitPar))
    lDensTrue <- dnorm(Par, mean  = InitPar[,4], sd = InitPar[,5], log=T)
    ldensity3 <- ldensity
    ldensity3[samp_i] <- ldensity+sum(lDensTrue)

    if(is.na(ldensity3[samp_i]))         {
      ldensity3[samp_i] <- -10^12 }
    if(ldensity3[samp_i]%in%c(-Inf,Inf)) {
      ldensity3[samp_i] <- -10^12  }
  }
  return(ldensity3)
}

