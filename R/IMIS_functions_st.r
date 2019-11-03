###################### Par = par_1
# Function for calculating likelihood

llikelihoodZ_st <-  function(samp_i,ParMatrix,loc) { # ParMatrix = ParInit
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  Par <- ParMatrix[samp_i,]
  # norm2unif
  Par2 <- pnorm(Par,0,1)
  # unif2true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  P <- P

  jj <- tryCatch({
    defpc<-def_prgchng(P)
    prms <-list()
    prms <- fin_param(P,loc,defpc)
    trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
    if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
    zz <- fin_cSim( nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]])  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                   Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
                   rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                   vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
                   HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
                   net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
                   mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                   TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                   LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                   EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
    #'if any output is missing or negative or if any model state population is negative
    if(sum(is.na(zz$Outputs[68,]))>0 | min(zz$Outputs[68,])<0 | min(zz$V1)<0 ) {
      lLik <- -10^12
  }
    else {
      ######  ####  ######  ######  ####  ######  ######  ####  ######
      M <- zz$Outputs
      colnames(M) <- prms[["ResNam"]]
      lLik <- 0
      st<-which(StateID$USPS==loc)

      ### ### ### TB SPECIFIC LIKELIHOODS
      # if(TB==1){
        ### ### ### TOTAL DIAGNOSED CASES 1993-2017  ### ### ### ### ### ### D
        v1   <- M[44:68,"NOTIF_ALL"]+M[44:68,"NOTIF_MORT_ALL"]
        addlik <- notif_tot_lLik_st(V=v1,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### ANN DECLINE IN CASES 1953-2015  ### ### ### ### ### ### D
        # v1b   <- M[4:44,"NOTIF_ALL"]+M[4:44,"NOTIF_MORT_ALL"]
        # addlik <- notif_decline_lLik_st(V=v1b,st=st); addlik
        # lLik <- lLik + addlik
        ### ### ### US CASES 1993-2016  ### ### ### ### ### ### D
        # v1a   <- rowSums(M[44:68,205:215]+M[44:68,216:226])
        # addlik <- US_notif_tot_lLik_st(V=v1a,st=st); addlik
        # lLik <- lLik + addlik
        # ### ### ### NUS CASES 1993-2016  ### ### ### ### ### ### D
        # v1b   <- rowSums(M[44:68,136:146]+M[44:68,189:199]) - (M[44:68,205:215]+M[44:68,216:226])
        # addlik <- NUS_notif_tot_lLik_st(V=v1b,st=st); addlik
        # lLik <- lLik + addlik
        ### ### ### US CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
        v2a   <- M[44:68,205:215]+M[44:68,216:226]
        addlik <- notif_age_us_lLik_st(V=v2a,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### FB CASES AGE DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
        v2b   <- (M[44:68,136:146]+M[44:68,189:199]) - (M[44:68,205:215]+M[44:68,216:226])
        addlik <- notif_age_fb_lLik_st(V=v2b,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### CASES FB DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
        v3   <-  cbind(M[44:68,148]+M[44:68,149]+(M[44:68,201]+M[44:68,202]),
                       M[44:68,147]+M[44:68,200])
        addlik <- notif_fb_lLik_st(V=v3,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### CASES FB, US 2010-2014  SLOPE ### ### ### ### ### ### D
        v3   <- cbind(M[63:68,148]+M[63:68,149]+(M[63:68,201]+M[63:68,202]),
                      M[63:68,147]+M[63:68,200])
        addlik <- notif_fbus_slp_lLik_st(V=v3,st=st); addlik
        lLik <- lLik + addlik

        ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
        v5b   <- cbind(M[44:67,151],M[44:67,150]) + cbind(M[44:67,204],M[44:67,203])
        addlik <- notif_us_hr_lLik_st(V=v5b,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
                v6   <- M[51:68,148:149]+M[51:68,201:202]
           addlik <- notif_fb_rec_lLik_st(V=v6, loc=loc); addlik
           lLik <- lLik + addlik
        ### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
        v11  <- M[44:66,132:134]
        addlik <- tx_outcomes_lLik_st(V=v11); addlik
        lLik <- lLik + addlik
        ### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
        v12  <- M[53,152]
        addlik <- tltbi_tot_lLik_st(V=v12,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### DIST LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
        v13  <- M[53,153:154]/M[53,152]
        addlik <- tltbi_dist_lLik_st(V=v13); addlik
        lLik <- lLik + addlik
        ### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ###  D
        v15  <- cbind(M[62,55:65],M[62,33:43]-M[62,55:65])
        #make this IGRA positive
        pIGRA<-1
        v15a<-v15*pIGRA
        Sens_IGRA <-c(.780,.675,.712,.789,.591)
        Spec_IGRA <-c(.979,.958,.989,.985,.931)
        names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
        v15b <- outer(v15a[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v15a[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))
        addlik <- ltbi_us_11_lLik(V=v15b)*2; addlik
        lLik <- lLik + addlik
        #' LTBI PREVALENCE BY AGE 2011, FB - index updated
        v16  <- cbind(M[62,66:76],M[62,44:54]-M[62,66:76])
        v16a <- v16*pIGRA
        #under age 5
        v16b <- (v16a[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v16a[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
        #over age 5
        v16c <- outer(v16a[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v16a[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
        v16d<-rbind(v16b,v16c)
        addlik <- ltbi_fb_11_lLik_st(V=v16d)*2; addlik
        lLik <- lLik + addlik
        ### ### ### TOTAL DEATHS WITH TB 1999-2016 ### ### ### ### ### ###  D
        v19  <- M[50:67,227:237]   ### THIS NOW ALL TB DEATHS
        addlik <- tbdeaths_lLik_st(V=v19,st=st); addlik
        lLik <- lLik + addlik
        ### ### ### TB DEATHS 1999-2016 BY AGE ### ### ### ### ### ###  D
        addlik <- tb_dth_age_lLik_st(V=v19); addlik
        lLik <- lLik + addlik
        ### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D
        v19b  <- rowSums(M[19:66,227:237])   ### THIS NOW ALL TB DEATHS
        addlik <- tbdeaths_decline_lLik_st(V=v19b); addlik
        lLik <- lLik + addlik

        #' LIKELIHOOD FOR BORGDORFF, FEREBEE & SUTHERLAND ESTIMATES
        v2456  <- list(prms[["Mpfast"]],prms[["Mrslow"]], prms[["rfast"]],prms[["rRecov"]])
        addlik <- borgdorff_lLik_st( Par=v2456); addlik
        lLik <- lLik + addlik
        addlik <- ferebee_lLik_st(Par=v2456); addlik
        lLik <- lLik + addlik
        addlik <- sutherland_lLik_st(Par=v2456); addlik
        lLik <- lLik + addlik

        # ### ### ### LIKELIHOOD FOR TIEMERSMA ESTS ### ### ### ### ### ### ~~~
        v35   <- c(P["rSlfCur"],P["muIp"])
        addlik <- tiemersma_lLik_st(Par=v35); addlik
        lLik <- lLik + addlik
      # } ### END OF TB LIKELIHOODS


    ### ### ### DEMOGRAPHIC LIKELIHOOD FUNCTIONS  ### ### ### ### ### ###
    ### ### ### TOTAL FB POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
    v17  <- M[,31]+M[,32]
    addlik <- tot_pop_yr_fb_lLik_st(V=v17,st=st); addlik
    lLik <- lLik + addlik
    ### ### ### TOTAL US POP EACH DECADE, FOR US  ### ### ### ### ### ###  D ResNam
    # v17b  <- M[,30]
    # addlik <- tot_pop_yr_us_lLik_st_00_10(V=v17b,st=st); addlik
    # lLik <- lLik + addlik
    ### ### ### TOTAL POP AGE DISTRIBUTION 2017  ### ### ### ### ### ### D
    v18  <- cbind(M[68,33:43],M[68,44:54])
    addlik <- tot_pop17_ag_fb_lLik_st(V=v18,st=st); addlik
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
    v21a<- v21  <- M[66:67,521:564]
    addlik <- mort_dist_lLik_st(V=v21a); addlik
    lLik <- lLik + addlik
    ### ### ###  ALL LIKELIHOODS DONE !!  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
} }, error = function(e) NA)
  if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
  if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

  return((lLik))  }


###################### local parallelization via multicore
llikelihood_st <- function(ParMatrix,loc,n_cores=1) {
  if(dim(as.data.frame(ParMatrix))[2]==1) {
    lLik <- llikelihoodZ_st(1,t(as.data.frame(ParMatrix)),loc=loc)
    } else {
      lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ_st,ParMatrix=ParMatrix,loc=loc,mc.cores=n_cores))
    }
  return((lLik)) }
