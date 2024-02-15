#' THIS FUNCTION INPUTS A TABLE OF PARAMETERS AND RUNS THE TB MODEL
#' AND GENERATES AN ARRAY OF OUTPUTS
#'@name  OutputsZint
#'@param samp_i which row of the parameter matrix to use
#'@param ParMatrix parameter matrix to use in the simulation
#'@param loc two-digit abbreviation for location
#'@param output_month when to output results
#'@param startyr year to start the simulation
#'@param endyr year to end the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param prg_chng vector of program change values
#'@param ttt_list list of targeted testing and treatment values
#'@param par2020 vector of 2020 adjustment parameters
#'@param return_params
#'@return results data frame of output
#'@export
OutputsZint <-  function(samp_i=1,ParMatrix,loc, output_month = 11, startyr=1950, endyr=2050,
                         Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0,
                         prg_chng=def_prgchng(Par[4,]), ttt_list=def_ttt(), care_cascade = def_care_cascade(),
                         par2020 = c(99,rep(0,12),rep(1,6)),
                         return_params = def_returnScenario())
{
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par1 <- as.numeric(ParMatrix);
    names(Par1) <- names(ParMatrix)
  } else {  Par1 <- as.numeric(ParMatrix[samp_i,]);
  names(Par1) <- colnames(ParMatrix) }

  P <- Par1

  Int1 <<- Int1;
  Int2 <<- Int2;
  Int3 <<- Int3;
  Int4 <<- Int4;
  Int5 <<- Int5;
  Scen1 <<- Scen1;
  Scen2 <<- Scen2;
  Scen3 <<- Scen3;

  ### add in the 2020 parameter adjustments
  ### basecase values are c(99, 0, 0, 1)
  names(par2020) <-  c("Immig",
                       "DxtKnot1", "DxtKnot2", "DxtKnot3", "DxtKnot4", "DxtKnot5", "DxtKnot6",
                       "TransKnot1", "TransKnot2", "TransKnot3", "TransKnot4", "TransKnot5", "TransKnot6",
                       "CaseFatKnot1", "CaseFatKnot2", "CaseFatKnot3", "CaseFatKnot4", "CaseFatKnot5", "CaseFatKnot6")

  prms <- list()
  prms <- param_init(P,loc,Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3,prg_chng,ttt_list,
                     immig = par2020["Immig"])
                     # return_months = return_params[["Immig"]][["return_months"]],
                     # multiplier = return_params[["Immig"]][["multiplier"]])

  ### adjust parameters for 2020 ###

  prms2020 <- adj_param_2020(rDxt = prms$rDxt,
                              NixTrans = prms$NixTrans,
                              par2020 = par2020,
                              return_params = return_params)

  # prms$rDxt[843:864,]<-prms$rDxt[843:864,] - (prms$rDxt[843:864,]*par2020["Dxt"])
  # prms$NixTrans[843:864]<- (1-par2020["Trans"])
  # # Bring up params to 50% by end of 2022 (smoothly)
  # for (riskgrp in 1:ncol(prms$rDxt)){
  #   prms$rDxt[865:888,riskgrp] <- seq(prms$rDxt[864,riskgrp],prms$rDxt[842,riskgrp], length.out=24)
  # }
  # prms$NixTrans[865:888] <- seq(prms$NixTrans[864],prms$NixTrans[842], length.out=24)
  #
  # RRmuTBPand <- rep(1,1812)
  # RRmuTBPand[843:888] <-c(rep(par2020["CaseFat"], 22), seq(par2020["CaseFat"], 1, length.out = 24))
  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")

  setup <- c(endyr-(startyr-1), length(func_ResNam()), output_month)

  m <- cSim(setup_pars = setup, rDxt               = prms2020[["rDxt"]]        , TxQualt      = prms[["TxQualt"]]     , InitPop       = prms[["InitPop"]],
            Mpfast   = prms[["Mpfast"]]   , ExogInf    = prms[["ExogInf"]]    , MpfastPI           = prms[["MpfastPI"]]    , Mrslow       = prms[["Mrslow"]]      , rrSlowFB      = prms[["rrSlowFB"]],
            rfast    = prms[["rfast"]]    , RRcurDef   = prms[["RRcurDef"]]   , rSlfCur            = prms[["rSlfCur"]]     , p_HR         = prms[["p_HR"]]        , vTMort        = prms[["vTMort"]],
            RRmuRF   = prms[["RRmuRF"]]   , RRmuHR     = prms[["RRmuHR"]]     , Birthst            = prms[["Birthst"]]     , HrEntEx      = prms[["HrEntEx"]]     , ImmNon        = prms[["ImmNon"]],
            ImmLat   = prms[["ImmLat"]]   , ImmAct     = prms[["ImmAct"]]     , ImmFst             = prms[["ImmFst"]]      , Int1Test     = prms[['Int1Test']]    , Int1Init     = prms[["Int1Init"]],
            Int1Tx     = prms[['Int1Tx']] , net_mig_usb  = prms[["net_mig_usb"]] , net_mig_nusb  = prms[["net_mig_nusb"]]  , RRmuTBPand   = prms2020[["RRmuTBPand"]]            ,
            mubt     = prms[["mubt"]]     , RelInf     = prms[["RelInf"]]     , RelInfRg           = prms[["RelInfRg"]]    , RRcrAG       = prms[["RRcrAG"]]      , Vmix          = prms[["Vmix"]],
            rEmmigFB = prms [["rEmmigFB"]], TxVec      = prms[["TxVec"]]      , TunTxMort          = prms[["TunTxMort"]]   , rDeft        = prms[["rDeft"]]       , ttt_samp_dist = prms[["ttt_sampling_dist"]],
            ttt_ag   = prms[["ttt_ag"]]   , ttt_na     = prms[["ttt_na"]]     , ttt_month          = prms[["ttt_month"]]   , ttt_pop_scrn = prms[["ttt_pop_scrn"]], ttt_ltbi      = prms[["ttt_ltbi"]],
            LtTxPar  = prms[["LtTxPar"]]  , LtDxPar_lt = prms[["LtDxPar_lt"]] , LtDxPar_nolt       = prms[["LtDxPar_nolt"]], rrTestLrNoTb = prms[["rrTestLrNoTb"]], rrTestHr = prms[["rrTestHr"]], rLtScrt      = prms[["rLtScrt"]]     , RRdxAge       = prms[["RRdxAge"]],
            ttt_ltbi_init=care_cascade[1], ttt_ltbi_comp=care_cascade[2], ttt_ltbi_eff=care_cascade[3], ttt_ltbi_sens=care_cascade[4], ttt_ltbi_spec=care_cascade[5], ttt_ltbi_accept=care_cascade[6],
            rRecov   = prms[["rRecov"]]   , pImmScen   = prms[["pImmScen"]]   , EarlyTrend         = prms[["EarlyTrend"]]  , pReTx        = prms[["pReTx"]]       , ag_den        = prms[["aging_denom"]],
            NixTrans = prms2020[["NixTrans"]] ,  dist_gen  = prms[["dist_gen"]]   , trans_mat_tot_ages = trans_mat_tot_ages)

  colnames(m$Outputs) <- func_ResNam();

  results <- as.matrix(m$Outputs)

  population <- as.matrix(m$OutputsI)

  colnames(m$OutputsI) <-colnames(m$Outputs) <- func_ResNam();

  # saveRDS(population, paste0("~/Documents/Tabby2 Paper/", loc, "_intv_pop_info_", Sys.time(), ".rds"))
  return(results)
}


#'wrapper function for the above function
#'@name  OutputsInt
#'@param loc two-digit abbreviation for location
#'@param ParMatrix parameters to use in the simulation
#'@param n_cores how many cores to use
#'@param output_month when to output results
#'@param startyr year to end the simulation
#'@param endyr year to end the simulation
#'@param startyr year to start the simulation
#'@param Int1 boolean for intervention 1
#'@param Int2 boolean for intervention 2
#'@param Int3 boolean for intervention 3
#'@param Int4 boolean for intervention 4
#'@param Int5 boolean for intervention 5
#'@param Scen1 boolean for scenario 1
#'@param Scen2 boolean for scenario 2
#'@param Scen3 boolean for scenario 3
#'@param prg_chng vector of program change values
#'@param ttt_list list of targeted testing and treatment values
#'@param par2020 vector of 2020 adjustment parameters
#'@param return_months
#'@param multiplier
#'@return out outputs
#'@export
OutputsInt <- function(loc,ParMatrix, n_cores=1, output_month = 11, startyr=1950,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0,prg_chng, ttt_list, par2020=c(99,rep(0,12),rep(1,6)), care_cascade = def_care_cascade(),
                       return_months =  865:888,
                       multiplier = 1) {
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    out <- OutputsZint(samp_i=1,ParMatrix=ParMatrix,loc=loc, output_month = output_month, endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3,
                       prg_chng=prg_chng, ttt_list=ttt_list, par2020 = par2020, care_cascade = care_cascade, return_months = return_months, multiplier = multiplier)
  } else {
    out0 <- mclapply(X=1:nrow(ParMatrix),FUN=OutputsZint,mc.cores=n_cores,
                     ParMatrix=ParMatrix, loc=loc, output_month = output_month, endyr=endyr,Int1=Int1,Int2=Int2,Int3=Int3,Int4=Int4,Int5=Int5,Scen1=Scen1,Scen2=Scen2,Scen3=Scen3,
                     prg_chng=prg_chng,ttt_list= ttt_list, par2020 = par2020, care_cascade = care_cascade, return_months = return_months, multiplier = multiplier)
    out <- array(NA,dim=c(length(out0),endyr-(startyr-1),length(func_ResNam())))

    for(i in 1:length(out0)){
      out[i,,] <- as.matrix(out0[[i]])
    }
    dimnames(out)[[3]]<-func_ResNam()
  }
  # if (sum(Int1,Int2,Int3,Int4,Int5,Scen1,Scen2,Scen3)==0) intv<-1;
  # if(Int1==1) intv<-2;if(Int2==1) intv<-3; if(Int3==1) intv<-4;
  # if(Int4==1) intv<-5; if(Int5==1) intv<-6; if(Scen1==1) intv<-7;
  # if(Scen2==1) intv<-8;if(Scen3==1) intv<-9;
  # saveRDS(out, file=paste("~MITUS/",loc,"_results_",intv,".rds",sep=""))
  # save(out,file=paste("/Users/nis100/MITUS/",loc,"_results_",intv,".rda",sep=""))
  return(out)
}
