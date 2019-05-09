calib<-function(samp_i,optim_mat, loc){
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

# P["TunMubt"]<-1;
# P["TunmuAg"]<-P["adj_ag1"]<-P["adj_ag11"]<-0
# # P[["pfast"]]<-.3
# P[["rslowH"]]<-oldpar["rslowH"]*2
# P[["rslow"]]<-oldpar["rslow"]*2

# P["TunrslowAge"]<-oldpar["TunrslowAge"]
# # p[[""]]
# P[["ImmigVol"]]<-1.1

  prms <-list()
  prms <- param(P)
  IP <- list()
  IP <- param_init(P)
  # data("trans_mat_nat",package="MITUS")
  # trans_mat_tot_ages<-trans_mat_tot_ages_nat
  # tm<-matrix(0,16,16)
  # diag(tm)<-1
  # trans_mat_tot_ages<<-matrix(tm,16,176)
  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  zz <- cSim_flow( nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]])  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
               Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
               rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
               vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
               HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
               net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
               mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
               TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
               LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
               EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    M <- zz$Outputs
    colnames(M) <- prms[["ResNam"]]
    v21a<- v21  <- M[1:2,521:564]
    for (i in 1:11){
      denom<-M[1:2,2+i]
      for (j in 1:ncol(v21)){
        v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
      } }
   print(v21a)
    pub_list<-list(prms[["Mpfast"]],prms[["Mrslow"]], prms[["rfast"]],prms[["rRecov"]])
    # print(pub_list)
    if (loc=="US"){
    calib_graphs(M, pub_list)
    } else calib_graphs_st(M,loc, pub_list)
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
