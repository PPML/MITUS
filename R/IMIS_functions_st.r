###################### Par = par_1
# Function for calculating likelihood

llikelihoodZ_st <-  function(samp_i,ParMatrix,loc, TB=1) { # ParMatrix = ParInit
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  if(min(dim(as.data.frame(ParMatrix)))==1) {
    Par <- as.numeric(ParMatrix);
    names(Par) <- names(ParMatrix)
  } else {  Par <- as.numeric(ParMatrix[samp_i,]);
  names(Par) <- colnames(ParMatrix) }   # norm2unif

  Par2 <- pnorm(Par,0,1)
  # unif2true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  P <- P
  # print(P)

  ### add in the 2020 parameter adjustments
  par2020 = c(0.4232265, 0.3707595, 0.1984619, 1.1158255)
  names(par2020) <- c("Immig", "Dxt", "Trans", "CaseFat")
  ### CALL THE DEFAULT CARE CASCADE FOR BELOW
  care_cascade <- def_care_cascade()
  CalibDatCases <<-readRDS(system.file("ST/ST_CalibDat_10-29-21.rds", package="MITUS"))

  jj <- tryCatch({
    prg_chng<-def_prgchng(P)
    prms <-list()
    prms <- param_init(P,loc,prg_chng=prg_chng, ttt_list=def_ttt())

    ### ADD IN THE 2020 ADJUSTMENT PARAMETERS; FOR OPTIMS THROUGH 2021 THESE
    ### WILL NOT HAVE ANY EFFECT, BUT ARE REQUIRED FOR MODEL RUNS
    prms$rDxt[843:864,]<-prms$rDxt[843:864,] - (prms$rDxt[843:864,]*par2020["Dxt"])
    prms$NixTrans[843:864]<- (1-par2020["Trans"])
    # Bring up params to 50% by end of 2022 (smoothly)
    for (riskgrp in 1:ncol(prms$rDxt)){
      prms$rDxt[865:888,riskgrp] <- seq(prms$rDxt[864,riskgrp],prms$rDxt[842,riskgrp], length.out=24)
    }
    prms$NixTrans[865:888] <- seq(prms$NixTrans[864],prms$NixTrans[842], length.out=24)

    RRmuTBPand <- rep(1,1812)
    RRmuTBPand[843:888] <-c(rep(par2020["CaseFat"], 22), seq(par2020["CaseFat"], 1, length.out = 24))
    ### END OF 2020 ADJUSTMENT
    ### POPULATION DISTRIBUTION REBLANCING FUNCTION
    trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
    if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
    ### END POPULATION REBLANCING
    ### DUE TO PARAMETER LIMITS WE CREATE A VECTOR OF VARIOUS SETUP VALUES
    ### THESE CORRESPOND TO NYRS, NRES, AND OUTCOME MONTH
    setup <- c(2021-(1950-1), length(func_ResNam()), 11)
    ### END SETUP PARAMETERS
    ### CALL THE MODEL FUNCTION
    zz <- cSim(setup_pars = setup, rDxt               = prms[["rDxt"]]        , TxQualt      = prms[["TxQualt"]]     , InitPop       = prms[["InitPop"]],
              Mpfast   = prms[["Mpfast"]]   , ExogInf    = prms[["ExogInf"]]    , MpfastPI           = prms[["MpfastPI"]]    , Mrslow       = prms[["Mrslow"]]      , rrSlowFB      = prms[["rrSlowFB"]],
              rfast    = prms[["rfast"]]    , RRcurDef   = prms[["RRcurDef"]]   , rSlfCur            = prms[["rSlfCur"]]     , p_HR         = prms[["p_HR"]]        , vTMort        = prms[["vTMort"]],
              RRmuRF   = prms[["RRmuRF"]]   , RRmuHR     = prms[["RRmuHR"]]     , Birthst            = prms[["Birthst"]]     , HrEntEx      = prms[["HrEntEx"]]     , ImmNon        = prms[["ImmNon"]],
              ImmLat   = prms[["ImmLat"]]   , ImmAct     = prms[["ImmAct"]]     , ImmFst             = prms[["ImmFst"]]      , Int1Test     = prms[['Int1Test']]    , Int1Init     = prms[["Int1Init"]],
              Int1Tx     = prms[['Int1Tx']] , net_mig_usb  = prms[["net_mig_usb"]] , net_mig_nusb  = prms[["net_mig_nusb"]]  , RRmuTBPand   = RRmuTBPand            ,
              mubt     = prms[["mubt"]]     , RelInf     = prms[["RelInf"]]     , RelInfRg           = prms[["RelInfRg"]]    , RRcrAG       = prms[["RRcrAG"]]      , Vmix          = prms[["Vmix"]],
              rEmmigFB = prms [["rEmmigFB"]], TxVec      = prms[["TxVec"]]      , TunTxMort          = prms[["TunTxMort"]]   , rDeft        = prms[["rDeft"]]       , ttt_samp_dist = prms[["ttt_sampling_dist"]],
              ttt_ag   = prms[["ttt_ag"]]   , ttt_na     = prms[["ttt_na"]]     , ttt_month          = prms[["ttt_month"]]   , ttt_pop_scrn = prms[["ttt_pop_scrn"]], ttt_ltbi      = prms[["ttt_ltbi"]],
              LtTxPar  = prms[["LtTxPar"]]  , LtDxPar_lt = prms[["LtDxPar_lt"]] , LtDxPar_nolt       = prms[["LtDxPar_nolt"]], rrTestLrNoTb = prms[["rrTestLrNoTb"]], rrTestHr = prms[["rrTestHr"]], rLtScrt      = prms[["rLtScrt"]]     , RRdxAge       = prms[["RRdxAge"]],
              ttt_ltbi_init=care_cascade[1], ttt_ltbi_comp=care_cascade[2], ttt_ltbi_eff=care_cascade[3], ttt_ltbi_sens=care_cascade[4], ttt_ltbi_spec=care_cascade[5], ttt_ltbi_accept=care_cascade[6],
              rRecov   = prms[["rRecov"]]   , pImmScen   = prms[["pImmScen"]]   , EarlyTrend         = prms[["EarlyTrend"]]  , pReTx        = prms[["pReTx"]]       , ag_den        = prms[["aging_denom"]],
              NixTrans = prms[["NixTrans"]] ,  dist_gen  = prms[["dist_gen"]]   , trans_mat_tot_ages = trans_mat_tot_ages)

    if(sum(is.na(zz$Outputs[68,]))>0 | min(zz$Outputs[68,])<0 | min(zz$V1)<0 ) {
      lLik <- -10^12
      print("condition 0")
    } else {
      ######  ####  ######  ######  ####  ######  ######  ####  ######
      M <- zz$Outputs
      colnames(M) <- prms[["ResNam"]]
      # print(M[69,])
      lLik <- 0
      st<-which(StateID$USPS==loc)
      # print(st)
      ### ### ### TB SPECIFIC LIKELIHOODS
      if(TB==1){
        ### ### ### TOTAL DIAGNOSED CASES 1993-2018  ### ### ### ### ### ### D
        v1   <- M[44:70,"NOTIF_ALL"]+M[44:70,"NOTIF_MORT_ALL"]
        addlik <- notif_tot_lLik_st(V=v1,st=st); addlik
        lLik <- lLik + addlik
        # print(paste("1:", lLik))
        ### ### ### ANN DECLINE IN CASES 1953-2015  ### ### ### ### ### ### D
        v1b   <- M[4:44,"NOTIF_ALL"]+M[4:44,"NOTIF_MORT_ALL"]
        addlik <- notif_decline_lLik_st(V=v1b,st=st); addlik
        lLik <- lLik + addlik
        # print(paste("2:", lLik))
        ### ### ### US CASES AGE DISTRIBUTION 5year 1994-2016  ### ### ### ### ### ### D
        v2a   <- M[46:70,205:215]+M[46:70,216:226]
        addlik <- notif_age_us_5yr_lLik_st(V=v2a,st=st); addlik
        lLik <- lLik + addlik
        # print(paste("5:", lLik))
        ### ### ### NUSB CASES AGE DISTRIBUTION 5year 1994-2013  ### ### ### ### ### ### D
        v2b   <- (M[46:70,136:146]+M[46:70,189:199]) - (M[46:70,205:215]+M[46:70,216:226])
        addlik <- notif_age_nus_5yr_lLik_st(V=v2b,st=st); addlik
        lLik <- lLik + addlik
        # print(addlik)
        # print(paste("6:", lLik))
        ### ### ### NUSB CASE DISTRIBUTION 1995-2019 ### ### ###
        # v3   <-  cbind(M[46:70,148]+M[46:70,149]+(M[46:70,201]+M[46:70,202]),
        #                M[46:70,147]+M[46:70,200])
        # addlik <- notif_fb_5yr_lLik_st(V=v3,st=st); addlik
        # lLik <- lLik + addlik
        # print(paste("7:", lLik))
        ### ### ### NUSB CASE TOTALS 1995-2019 ### ### ###
        v3a   <-  M[46:70,148]+M[46:70,149]+(M[46:70,201]+M[46:70,202])
        addlik <- notif_fb_5yr_lik(V=v3a,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### USB CASE TOTALS 1995-2019 ### ### ###
        v3b   <- M[46:70,147]+M[46:70,200]
        addlik <- notif_us_5yr_lik(V=v3b,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### CASES NUSB, US 2010-2014  SLOPE ### ### ### ### ### ### D
        ### ### ### Removed as we now use the five year totals
        if (loc %in% c("MT")){
          lLik <-lLik;
        } else {
        v4   <- cbind(M[61:70,148]+M[61:70,149]+(M[61:70,201]+M[61:70,202]),
                      M[61:70,147]+M[61:70,200])
        addlik <- notif_fbus_slp_lLik_st(V=v4,st=st); addlik
        lLik <- lLik + addlik
        }
        # print(paste("8:", lLik))
        ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
        ### ### ### CHECK THIS BECAUSE IT'S A HUGE VALUE
        if (loc %in% c("RI", "ID", "NH", "VT", "WY")){
          lLik <-lLik;
        } else {
        v5   <- M[46:70,151] + M[46:70,204]
        v5b  <- rbind(sum(v5[1:5]),sum(v5[6:10]), sum(v5[11:15]),
                      sum(v5[16:20]), sum(v5[21:25]))
        addlik <- notif_hr_lLik_st(V=v5b,st=st); addlik
        lLik <- lLik + addlik; addlik
        }
        # print(paste("9:", lLik))
        ### ### ### CASES RCT TRANS DISTRIBUTION 2015-2018  ### ### ### ### ### ### D
        v4a <- (sum(M[66:69,184:185])/sum(M[66:69,168:169]))
        v4 <- c(v4a,1-v4a)
        addlik <- recent_trans_dist_lLik_st(V=v4,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### CASES NUSB RECENT ENTRY DISTRIBUTION 5YR 1994-2018 ### ### ### ### ### ### D
        ### ### ### recent, not recent
        v6   <-M[45:69,148]+M[45:69,201]
        v6b  <- rbind(sum(v6[1:5]),sum(v6[6:10]), sum(v6[11:15]),
                      sum(v6[16:20]), sum(v6[21:25]))
        addlik <- notif_fb_rec_lLik_st(V=v6b, st=st); addlik
        lLik <- lLik + addlik
        # print(paste("10:", lLik))
        ### ### ### TREATMENT OUTCOMES 1993-2015  ### ### ### ### ### ### D
        v11  <- M[44:66,132:134]
        addlik <- tx_outcomes_lLik_st(V=v11); addlik
        lLik <- lLik + addlik
        # print(paste("11:", lLik))
        ### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
        v12  <- M[53,152]
        addlik <- tltbi_tot_lLik_st(V=v12,st=st); addlik
        lLik <- lLik + addlik
        # print(paste("12:", lLik))
        ### ### ### DIST LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
        v13  <- M[53,153:154]/M[53,152]
        addlik <- tltbi_dist_lLik_st(V=v13); addlik
        lLik <- lLik + addlik
        # print(paste("13:", lLik))
        ### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ###  D
        v15  <- cbind(M[62,55:65],M[62,33:43]-M[62,55:65])
        #make this IGRA positive
        pIGRA<-1
        v15a<-v15*pIGRA
        Sens_IGRA <-c(.780,.780,.712,.789,.789)
        Spec_IGRA <-c(.979,.979,.989,.985,.985)
        names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
        v15b <- outer(v15a[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v15a[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))
        addlik <- ltbi_us_11_lLik_st(V=v15b)*2; addlik
        lLik <- lLik + addlik
        # print(paste("14:", lLik))
        #' LTBI PREVALENCE BY AGE 2011, NUSB - index updated
        v16  <- cbind(M[62,66:76],M[62,44:54]-M[62,66:76])
        v16a <- v16*pIGRA
        #under age 5
        v16b <- (v16a[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v16a[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
        #over age 5
        v16c <- outer(v16a[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v16a[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
        v16d<-rbind(v16b,v16c)
        addlik <- ltbi_fb_11_lLik_st(V=v16d)*2; addlik
        lLik <- lLik + addlik
        # # print(paste("15:", lLik))
        # ### ### ### TOTAL DEATHS WITH TB 1999-2018 ### ### ### ### ### ###  D
        v19  <- M[50:70,227:237]   ### THIS NOW ALL TB DEATHS
        addlik <- tbdeaths_lLik_st(V=v19,st=st); addlik
        lLik <- lLik + addlik
        # # print(paste("16:", lLik))
        # ### ### ### TB DEATHS 1999-2016 BY AGE ### ### ### ### ### ###  D
        addlik <- tb_dth_age_lLik_st(V=v19); addlik
        lLik <- lLik + addlik
        # print(paste("17:", lLik))
        ### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D
        ###not working
        # v19b  <- rowSums(M[19:66,227:237])   ### THIS NOW ALL TB DEATHS
        # addlik <- tbdeaths_decline_lLik_st(V=v19b); addlik
        # lLik <- lLik + addlik

        #' LIKELIHOOD FOR BORGDORFF, FEREBEE & SUTHERLAND ESTIMATES
        v2456  <- list(prms[["Mpfast"]],prms[["Mrslow"]], prms[["rfast"]],prms[["rRecov"]])
        addlik <- borgdorff_lLik_st(Par=v2456); addlik
        lLik <- lLik + addlik
        # print(paste("18:", lLik))
        addlik <- ferebee_lLik_st(Par=v2456); addlik

         lLik <- lLik + addlik
        # print(paste("19:", lLik))
        addlik <- sutherland_lLik_st(Par=v2456); addlik
        lLik <- lLik + addlik
        # print(paste("20:", lLik))

        # ### ### ### LIKELIHOOD FOR TIEMERSMA ESTS ### ### ### ### ### ### ~~~
        v35   <- c(P["rSlfCur"],P["muIp"])
        addlik <- tiemersma_lLik_st(Par=v35); addlik
        lLik <- lLik + addlik
        # print(paste("21:", lLik))
      } ### END OF TB LIKELIHOODS

    ### ### ### DEMOGRAPHIC LIKELIHOOD FUNCTIONS  ### ### ### ### ### ###
    ### ### ### TOTAL NUSB POP EACH DECADE, FOR NUSB  ### ### ### ### ### ###  D
    v17  <- M[,31]+M[,32]
    addlik <- tot_pop_yr_fb_lLik_st(V=v17,st=st); addlik
    lLik <- lLik + addlik
    ### ### ### TOTAL US POP EACH DECADE, FOR US  ### ### ### ### ### ###  D ResNam
    v17b  <- M[,30]
    addlik <- tot_pop_yr_us_lLik_st(V=v17b,st=st); addlik
    lLik <- lLik + addlik
    ### ### ### TOTAL POP AGE DISTRIBUTION 2019  ### ### ### ### ### ### D
    # v18a  <- cbind(M[70,33:43],M[70,44:54])
    # addlik <- tot_pop19_ag_fb_lLik_st(V=v18,st=st); addlik
    # lLik <- lLik + addlik
    ### ### ### TOTAL POP AGE DISTRIBUTION 2017-2019  ### ### ### ### ### ### D
    v18  <- cbind(colSums(M[68:70,33:43]),colSums(M[68:70,44:54]))
    addlik <- tot_pop1719_ag_fb_lLik_st(V=v18,st=st); addlik
    lLik <- lLik + addlik
    #' ### ### ### HOMELESS POP 2010  ### ### ### ### ### ###
    v23b  <- M[61,29]
    addlik <- homeless_10_lLik_st(V=v23b,st=st); addlik
    lLik <- lLik + addlik
    #' #' Total DEATHS 2016
    v20a  <- sum(M[67,121:131])
    addlik <- dth_tot_lLik_st(V=v20a,st=st); addlik
    lLik <- lLik + addlik
    #' #' #' #' Total DEATHS 2015-2016 BY AGE
    v20b  <- M[66:67,121:131]
    addlik <- tot_dth_age_lLik_st(V=v20b,st=st); addlik
    lLik <- lLik + addlik
    #' #' #' Mort_dist 2016 dirchlet
    v21a<- M[66:67,521:564]
    addlik <- mort_dist_lLik_st(V=v21a); addlik
    lLik <- lLik + addlik
    ### ### ###  ALL LIKELIHOODS DONE !!  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
}
}, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2); print("condition 1") }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2); print("condition 2") }
# print(paste("log likelihood is ", lLik))
  return((lLik))  }


###################### local parallelization via multicore
llikelihood_st <- function(ParMatrix,loc,n_cores=1, TB=1) {
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    lLik <- llikelihoodZ_st(1,t(as.data.frame(ParMatrix)),loc=loc, TB=TB)
    } else {
      lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ_st,ParMatrix=ParMatrix,loc=loc,mc.cores=n_cores, TB=TB))
    }
  return((lLik)) }
