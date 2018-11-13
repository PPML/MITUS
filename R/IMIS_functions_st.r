library(mnormt)
library(parallel)
library(lhs)

###################### Par = par_1
# Function for calculating likelihood

  llikelihoodZ_st <-  function(samp_i,ParMatrix,loc="CA") { # ParMatrix = ParInit
    StateID<-as.data.frame(read.csv(file="inst/extdata/state_ID.csv", header = TRUE))

    # model_inputs<-paste0(loc,"_ModelInputs_11-13-18")
    # data(list=model_inputs, package = 'MITUS')

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
      prms <-list()
      prms <- param(P)
      IP <- list()
      IP <- param_init(P)
      # tm<-matrix(0,16,16)
      # diag(tm)<-1
      # trans_mat_tot_ages<<-matrix(tm,16,176)
      trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])

      zz <- cSim(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
                   Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
                   rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
                   vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
                   HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
                   mubt       = prms[["mubt"]]      , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
                   TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
                   LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
                   EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
      if(sum(is.na(zz$Outputs[66,]))>0 | min(zz$Outputs[66,])<0 | min(zz$V1)<0 ) { lLik <- -10^12  } else {

      ######  ####  ######  ######  ####  ######  ######  ####  ######
        M <- zz$Outputs
        colnames(M) <- prms[["ResNam"]]
        lLik <- 0
        st<-which(StateID$USPS==loc)
        ### ### ### TOTAL DIAGNOSED CASES 1993-2016  ### ### ### ### ### ### D
      v1   <- M[44:67,"NOTIF_ALL"]+M[44:67,"NOTIF_MORT_ALL"]
      addlik <- notif_tot_lLik_st(V=v1,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### ANN DECLINE IN CASES 1953-2015  ### ### ### ### ### ### D
      v1b   <- M[4:44,"NOTIF_ALL"]+M[4:44,"NOTIF_MORT_ALL"]
      addlik <- notif_decline_lLik_st(V=v1b); addlik
      lLik <- lLik + addlik
      ### ### ### US CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
      v2a   <- M[44:67,205:215]+M[44:67,216:226]
      addlik <- notif_age_us_lLik_st(V=v2a,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### FB CASES AGE DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
      v2b   <- (M[44:67,136:146]+M[44:67,189:199]) - (M[44:67,205:215]+M[44:67,216:226])
      addlik <- notif_age_fb_lLik_st(V=v2b,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### CASES FB DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
      v3   <-  cbind(M[44:67,148]+M[44:67,149]+(M[44:67,201]+M[44:67,202]),
                     M[44:67,147]+M[44:67,200])
      addlik <- notif_fb_lLik_st(V=v3,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### CASES FB, US 2010-2014  SLOPE ### ### ### ### ### ### D
      v3   <- cbind(M[62:67,148]+M[62:67,149]+(M[62:67,201]+M[62:67,202]),
                    M[62:67,147]+M[62:67,200])
      addlik <- notif_fbus_slp_lLik_st(V=v3,st=st); addlik
      lLik <- lLik + addlik

      ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
      v5b   <- cbind(M[44:67,151],M[44:67,150]) + cbind(M[44:67,204],M[44:67,203])
      addlik <- notif_us_hr_lLik_st(V=v5b,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
   #         v6   <- M[44:65,148:149]+M[44:65,201:202]

  #    addlik <- notif_fb_rec_lLik_st(V=v6); addlik
  #    lLik <- lLik + addlik
   ### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
      v11  <- M[44:63,132:134]
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
      v15a <- outer(v15[,1],c(P[["SensLt"]],1-P[["SensLt"]]))+outer(v15[,2],c(1-P[["SpecLt"]],P[["SpecLt"]]))
      addlik <- ltbi_us_11_lLik_st(V=v15a)*2; addlik
      lLik <- lLik + addlik
      ### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ### D
      v16  <- cbind(M[62,66:76],M[62,44:54]-M[62,66:76])
      v16a <- outer(v16[,1],c(P[["SensLt"]],1-P[["SensLt"]]))+outer(v16[,2],c(1-P[["SpecLt"]],P[["SpecLt"]]))
      addlik <- ltbi_fb_11_lLik_st(V=v16a)*2; addlik
      lLik <- lLik + addlik
      ### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
      v17  <- M[,31]+M[,32]
      addlik <- tot_pop_yr_fb_lLik_st(V=v17,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### TOTAL POP EACH DECADE, FOR US  ### ### ### ### ### ###  D ResNam
      v17b  <- M[,30]
      addlik <- tot_pop_yr_us_lLik_st_00_10(V=v17b,st=st); addlik
      lLik <- lLik + addlik
      ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ### D
      v18  <- cbind(M[65,33:43],M[65,44:54])
      addlik <- tot_pop14_ag_fb_lLik_st(V=v18,st=st); addlik
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

      ### ### ### HOMELESS POP 2010  ### ### ### ### ### ###
      v23b  <- M[61,29]
      addlik <- homeless_10_lLik_st(V=v23b,st=st); addlik
      lLik <- lLik + addlik

      #' Total DEATHS 1979-2016
      v20a  <- rowSums(M[30:67,121:131])
      v20a<-v20a*1e6
      addlik <- dth_tot_lLik_st(V=v20a,st=st); addlik
      lLik <- lLik + addlik
      #' #' Total DEATHS 1999-2016 BY AGE
      #' v20b  <- M[50:67,121:131]
      #' addlik <- tot_dth_age_lLik(V=v20b); addlik
      #' lLik <- lLik + addlik
      #' #' Mort_dist 2016
      v21a<- v21  <- M[51:67,521:564]
      for (i in 1:11){
        denom<-M[51:67,2+i]
        for (j in 1:ncol(v21)){
          v21a[,(1:4)+4*(i-1)]<-v21[,(1:4)+4*(i-1)]/denom
        } }
      addlik <- mort_dist_lLik(V=v21a); addlik
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
      ### ### ###  ALL LIKELIHOODS DONE !!  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

      } }, error = function(e) NA)
      if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
      if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }

      return((lLik))  }


  ###################### local parallelization via multicore
  llikelihood_st <- function(ParMatrix,st,n_cores=1) {
    if(dim(as.data.frame(ParMatrix))[2]==1) {
      lLik <- llikelihoodZ_st(1,t(as.data.frame(ParMatrix)),st=st)  } else {
      lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ_st,ParMatrix=ParMatrix,st=st,mc.cores=n_cores))
      }
    return((lLik)) }






