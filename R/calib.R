calib<-function(samp_i,optim_mat, loc, pdf=TRUE, cex.size=1){
  if(min(dim(as.data.frame(optim_mat)))==1) {
  Par <- as.numeric(optim_mat);
  names(Par) <- names(optim_mat)
} else {  Par <- as.numeric(optim_mat[samp_i,]);
names(Par) <- colnames(optim_mat) }  ##previously, the distribution of parameters were transformed to normal distribution in
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
# P <- Par[8,]
# P["TunmuAg"]<-P["adj_ag1"]<-P["adj_ag11"]<-0
# P[["pfast"]]<-.08
# P[["rslowH"]]<-oldpar["rslowH"]*2
# P[["rslow"]]<-oldpar["rslow"]*2

# P["TunrslowAge"]<-oldpar["TunrslowAge"]
# P[["rSlfCur"]]<-1
# P[["ImmigVol"]]<-1;
# P["RRtbprev"]<-2;
# P["TunLtbiTrend"]<-P["TunLtbiTrend"]*4
# P[["TunNetMig"]]<-1
# P[["TunMubt"]]<-.1
# P[["pImAct"]]<-.1
# P[["LtbiPar1"]]<-1
# P[["LtbiPar2"]]<-1
# P[["sigmaFb"]]<-.5 #.9#P[["sigmaFb"]]*10
# P[["sigmaHr"]]<-.5# .9#P[["sigmaHr"]]/10
#  P[["RelCrHr"]]<-  P[["RelCrHr"]]/2


prms <-list()
prms<-param_init(PV=P, loc=loc, prg_chng = def_prgchng(P), ttt_list = def_ttt())
# data("trans_mat_nat",package="MITUS")
# trans_mat_tot_ages<-trans_mat_tot_ages_nat
# tm<-matrix(0,16,16)
# diag(tm)<-1
# trans_mat_tot_ages<<-matrix(tm,16,176)
trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
rownames(trans_mat_tot_ages) <-  paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
colnames(trans_mat_tot_ages) <-  rep(paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4)),11)

if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
zz <- cSim( nYrs       = 2020-1950         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
            Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
            rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]] ,  dist_gen = prms[["dist_gen"]]    ,
            vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
            HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
            net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
            mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], RRcrAG = prms[["RRcrAG"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
            TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
            LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
            ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
            EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
#'if any output is missing or negative or if any model state population is negative
  M <- zz$Outputs
    colnames(M) <- prms[["ResNam"]]
    v21a<- v21  <- M[51:52,521:564]
    for (i in 1:11){
      denom<-M[51:52,2+i]
      for (j in 1:ncol(v21)){
        v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
      } }
   print(v21a)
    pub_list<-list(prms[["Mpfast"]],prms[["Mrslow"]], prms[["rfast"]],prms[["rRecov"]])
    # print(pub_list)
    if (loc=="US"){
    calib_graphs(M, pub_list)
    } else calib_graphs_st(M,loc, pub_list, pdf=TRUE, cex.size = .7)
    return(M)
}
calib_2020<-function(samp_i,optim_mat, loc, pdf=TRUE, cex.size=1){

  if(min(dim(as.data.frame(optim_mat)))==1) {
    Par <- as.numeric(optim_mat);
    names(Par) <- names(optim_mat)
  } else {  Par <- as.numeric(optim_mat[samp_i,]);
  names(Par) <- colnames(optim_mat) }  ##previously, the distribution of parameters were transformed to normal distribution in
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
  # P <- Par[1,]
  # par2020 <- c(0.088, 0.43, 0.28)
  # par2020 <- c (0.1797244, 0.4360917, 0.2754278)
  # par2020 <- c(0.1800718,.4351458,0.3015594)
  # par2020 <- c(0.3114697,0.4184191,0.2535374)
  # *.1, *.1, .1
  # par2020 <- c(0.3956915,0.4018355, 0.2696539)
  # par2020 <- c(0.2059258,0.4313880, 0.29800063)
  #as of 1/14/22
  # par2020 <- c(0.4064649,0.4010572,0.2438105)
  #as of 1/24/22
  par2020 <- c(0.5168886,0.3827480,0.2232671)
  names(par2020) <- c("Immig", "Dxt", "Trans")

  prg_chng<-def_prgchng(P)

  prms <-list()
  prms <- param_init(P,loc,prg_chng=prg_chng, ttt_list=def_ttt(), immig = par2020["Immig"])
  prms$rDxt[843:864,]<-prms$rDxt[843:864,] - (prms$rDxt[843:864,]*par2020["Dxt"])
  prms$NixTrans[843:864]<- (1-par2020["Trans"])  # data("trans_mat_nat",package="MITUS")
  # Bring up params to 50% by end of 2022 (smoothly)
  for (riskgrp in 1:ncol(prms$rDxt)){
    prms$rDxt[865:888,riskgrp] <- seq(prms$rDxt[864,riskgrp],prms$rDxt[842,riskgrp], length.out=24)
  }
    prms$NixTrans[865:888] <- seq(prms$NixTrans[864],prms$NixTrans[842], length.out=24)

  # trans_mat_tot_ages<-trans_mat_tot_ages_nat
  # tm<-matrix(0,16,16)
  # diag(tm)<-1
  # trans_mat_tot_ages<<-matrix(tm,16,176)
  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  rownames(trans_mat_tot_ages) <-  paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
  colnames(trans_mat_tot_ages) <-  rep(paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4)),11)

  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  zz <- cSim( nYrs       = 2050-1950         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
              Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
              rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]] ,  dist_gen = prms[["dist_gen"]]    ,
              vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
              HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
              net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
              mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], RRcrAG = prms[["RRcrAG"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
              TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
              LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
              ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
              EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
  #'if any output is missing or negative or if any model state population is negative
  M <- zz$Outputs
  colnames(M) <- prms[["ResNam"]]
  v21a<- v21  <- M[51:52,521:564]
  for (i in 1:11){
    denom<-M[51:52,2+i]
    for (j in 1:ncol(v21)){
      v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
    } }
  print(v21a)
  pub_list<-list(prms[["Mpfast"]],prms[["Mrslow"]], prms[["rfast"]],prms[["rRecov"]])
  # print(pub_list)
  if (loc=="US"){
    calib_graphs_2020(M, pub_list)
  } else calib_graphs_st_2020(M,loc, pdf=TRUE, cex.size = .7)
  future_graphs_st(loc=loc, df=M, cex.size = .7)
   return(M)
}

# calib_par<-function(samp_i,par_vec, loc){
#
#   P <-oldpar
#
#   # P["TunMubt"]<-1;
#   # P["TunmuAg"]<-P["adj_ag1"]<-0
#   # P[["rslow"]]<-oldpar["rslow"]
#   # P[["rslowH"]]<-oldpar["rslowH"]
#   # P["TunrslowAge"]<-oldpar["TunrslowAge"]
#   # # p[[""]]
#   # P[["ImmigVol"]]<-1.1
#
#   prms <-list()
#   prms <- param(P)
#   IP <- list()
#   IP <- param_init(P)
#   # data("trans_mat_nat",package="MITUS")
#   # trans_mat_tot_ages<-trans_mat_tot_ages_nat
#
#   trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
#   if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
#   zz <- cSim( nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]])  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
#               Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
#               rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
#               vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
#               HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
#               net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
#               mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
#               TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
#               LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
#               EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
#   M <- zz$Outputs
#   colnames(M) <- prms[["ResNam"]]
#   v21a<- v21  <- M[51:67,521:564]
#   for (i in 1:11){
#     denom<-M[51:67,2+i]
#     for (j in 1:ncol(v21)){
#       v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
#     } }
#   print(v21a)
#   if (loc=="US"){
#     calib_graphs(M)
#   } else calib_graphs_st(M,loc)
#   return(M)
# }
