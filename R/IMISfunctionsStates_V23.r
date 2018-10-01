# library(mnormt)
# library(parallel)
# library(lhs)
#
# ###################### Par = par_1
# # Function for calculating likelihood
#
#
#   llikelihoodZ <-  function(samp_i,ParMatrix) { # ParMatrix = ParInit
#     Par <- ParMatrix[samp_i,]
#     # norm2unif
#     Par2 <- pnorm(Par,0,1)
#     # unif2true
#     Par3 <- Par2
#     Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
#     Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
#     Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
#     P[ii] <- Par3
#     P <<- P
#
#     jj <- tryCatch({
#       source("ParamStates_V23.r")
#       zz <- cSim( nYrs     = 2017-1950, nRes     = length(ResNam), rDxt      = rDxt     , TxQualt   = TxQualt   , InitPop   = InitPop,
#                   p_HR     = p_HR     , Mpfast   = Mpfast        , ExogInf   = ExogInf  , MpfastPI  = MpfastPI  , RelInfHivt = RelInfHivt,
#                   pimmed   = pimmed   , MpSmPos  = MpSmPos       , HivHrPar  = HivHrPar , RRmuHR    = RRmuHR    , Mrslow    = Mrslow,
#                   rfast    = rfast    , RRcurDef = RRcurDef      , fAR       = fAR      , VrSlfCur  = VrSlfCur  , VrSmConv  = VrSmConv,
#                   vTMort   = vTMort   , vHMort   = vHMort        , Birthst   = Birthst  , ImmNon    = ImmNon    , ImmLat    = ImmLat, ImmFst    = ImmFst,
#                   ImmAct   = ImmAct   , DrN      = DrN           , DrE       = DrE      , TxExpAge  = TxExpAge  , p_Imm_SP  = p_Imm_SP,
#                   mubt     = mubt     , RelInf   = RelInf        , RelInfRg  = RelInfRg , Vmix      = Vmix      , rArtInit  = rArtInit ,
#                   vHxtoHy  = vHxtoHy  , rArtDef  = rArtDef       , rHIVt     = rHIVt    , rEmmigFB  = rEmmigFB  , pDstt     = pDstt   ,
#                   TxMat    = TxMat    , TunTxMort = TunTxMort    , rDeft     = rDeft    , rDeftH    = rDeftH    , pReTx     = pReTx    ,
#                   LtTxPar  = LtTxPar  , LtDxPar  = LtDxPar       , rLtScrt   = rLtScrt  , HrEntEx   = HrEntEx   , muTbH     = muTbH,
#                   RRdxAge  = RRdxAge  , rRecov   = rRecov        , pImmScen  = pImmScen, EarlyTrend = EarlyTrend, rrSlowFB = rrSlowFB,
#                   net_mig_usb = net_mig_usb,  net_mig_nusb = net_mig_nusb)
#       if(sum(is.na(zz$Outputs[66,]))>0 | min(zz$Outputs[66,])<0 | min(zz$V1)<0 ) { lLik <- -10^12  } else {
#
#       ######  ####  ######  ######  ####  ######  ######  ####  ######
#       M <- zz$Outputs; colnames(M) <- ResNam;   lLik <- 0
#       ### ### ### TOTAL DIAGNOSED CASES 1993-2016  ### ### ### ### ### ### D
#       v1   <- M[44:67,"NOTIF_ALL"]+M[44:67,"NOTIF_MORT_ALL"]
#       addlik <- notif_tot_lik(V=v1); addlik
#       lLik <- lLik + addlik
#       ### ### ### ANN DECLINE IN CASES 1953-2015  ### ### ### ### ### ### D
#       v1b   <- M[4:44,"NOTIF_ALL"]+M[4:44,"NOTIF_MORT_ALL"]
#       addlik <- notif_decline_lLik(V=v1b); addlik
#       lLik <- lLik + addlik
#       ### ### ### US CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
#       v2a   <- M[44:67,235:245]+M[44:67,246:256]
#       addlik <- notif_age_us_lLik(V=v2a); addlik
#       lLik <- lLik + addlik
#       ### ### ### FB CASES AGE DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
#       v2b   <- (M[44:67,137:147]+M[44:67,216:226]) - (M[44:67,235:245]+M[44:67,246:256])
#       addlik <- notif_age_fb_lLik(V=v2b); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES FB DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
#       v3   <- cbind(M[44:67,154]+M[44:67,155]+(M[44:67,233]+M[44:67,234]),
#                   M[44:67,152]+M[44:67,153]+(M[44:67,231]+M[44:67,232]))
#       addlik <- notif_fb_lLik(V=v3); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES FB, US 2010-2014  SLOPE ### ### ### ### ### ### D
#       v3   <- cbind(M[62:67,154]+M[62:67,155]+(M[62:67,233]+M[62:67,234]),
#                   M[62:67,152]+M[62:67,153]+(M[62:67,231]+M[62:67,232]))
#       addlik <- notif_fbus_slp_lLik(V=v3); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES TX HISTORY DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
#       v4   <- M[48:67,148:149]+M[48:67,227:228]
#       addlik <- notif_prev_lLik(V=v4); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES HIV DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
#       v5   <- M[48:67,150:151]+M[48:67,229:230]
#       addlik <- notif_hiv_lLik(V=v5); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
#       v5b   <- cbind(M[48:67,153],M[48:67,152]) + cbind(M[48:67,232],M[48:67,231])
#       addlik <- notif_us_hr_lLik(V=v5b); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
#    #   v6   <- M[44:65,154:155]+M[44:65,233:234]
#   #    addlik <- notif_fb_rec_lLik(V=v6); addlik
#   #    lLik <- lLik + addlik
#       ### ### ### CASES PCT MDR RES US N 1997-2016  ### ### ### ### ### ### D
#       v7   <- M[48:67,156:160]
#       addlik <- notif_mdr_us_n_lLik(V=v7); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES PCT MDR RES US E 1997-2016  ### ### ### ### ### ### D
#       v8   <- M[48:67,161:165]
#       addlik <- notif_mdr_us_e_lLik(V=v8); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES PCT MDR RES FB N 1997-2016  ### ### ### ### ### ### D
#       v9   <- M[48:67,166:170]
#       addlik <- notif_mdr_fb_n_lLik(V=v9); addlik
#       lLik <- lLik + addlik
#       ### ### ### CASES PCT MDR RES FB E 1997-2016  ### ### ### ### ### ### D
#       v10  <- M[48:67,171:175]
#       addlik <- notif_mdr_fb_e_lLik(V=v10); addlik
#       lLik <- lLik + addlik
#       ### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
#       v11  <- M[44:63,133:135]
#       addlik <- tx_outcomes_lLik(V=v11); addlik
#       lLik <- lLik + addlik
#       ### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
#       v12  <- M[53,176]
#       addlik <- tltbi_tot_lLik(V=v12); addlik
#       lLik <- lLik + addlik
#       ### ### ### DIST LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
#       v13  <- M[53,177:179]/M[53,176]
#       addlik <- tltbi_dist_lLik(V=v13); addlik
#       lLik <- lLik + addlik
#       ### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ###  D
#       v15  <- cbind(M[62,56:66],M[62,34:44]-M[62,56:66])
#       v15a <- outer(v15[,1],c(SensLt,1-SensLt))+outer(v15[,2],c(1-SpecLt,SpecLt))
#       addlik <- ltbi_us_11_lLik(V=v15a)*2; addlik
#       lLik <- lLik + addlik
#       ### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ### D
#       v16  <- cbind(M[62,67:77],M[62,45:55]-M[62,67:77])
#       v16a <- outer(v16[,1],c(SensLt,1-SensLt))+outer(v16[,2],c(1-SpecLt,SpecLt))
#       addlik <- ltbi_fb_11_lLik(V=v16a)*2; addlik
#       lLik <- lLik + addlik
#       ### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
#       v17  <- M[,32]+M[,33]
#       addlik <- tot_pop_yr_fb_lLik(V=v17); addlik
#       lLik <- lLik + addlik
#       ### ### ### TOTAL POP EACH DECADE, FOR US  ### ### ### ### ### ###  D ResNam
#       v17b  <- M[,30]+M[,31]
#       addlik <- tot_pop_yr_us_lLik_00_10(V=v17b); addlik
#       lLik <- lLik + addlik
#       ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ### D
#       v18  <- cbind(M[65,34:44],M[65,45:55])
#       addlik <- tot_pop14_ag_fb_lLik(V=v18); addlik
#       lLik <- lLik + addlik
#       ### ### ### TOTAL DEATHS WITH TB 1999-2016 ### ### ### ### ### ###  D
#       v19  <- M[50:67,279:289]   ### THIS NOW ALL TB DEATHS
#       addlik <- tbdeaths_lik(V=v19); addlik
#       lLik <- lLik + addlik
#       ### ### ### TB DEATHS 1999-2016 BY AGE ### ### ### ### ### ###  D
#       addlik <- tb_dth_age_lLik(V=v19); addlik
#       lLik <- lLik + addlik
#       ### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D
#       v19b  <- rowSums(M[19:66,279:289])   ### THIS NOW ALL TB DEATHS
#       addlik <- tbdeaths_decline_lLik(V=v19b); addlik
#       lLik <- lLik + addlik
#       ### ### ### TOTAL HIV DEATHS 2008-2015  ### ### ### ### ### ### D
#       v20  <- rowSums(M[59:66,111:121])
#       addlik <- hivdeaths_lik(V=v20); addlik
#       lLik <- lLik + addlik
#       ### ### ### HIV PREVALENCE 2010-2015  ### ### ### ### ### ### D
#       v21  <- rowSums(M[61:66,80:88])/rowSums(M[61:66,5:13])
#       addlik <- hiv_prev_by_year_lLik(V=v21); addlik
#       lLik <- lLik + addlik
#       ### ### ###  HIV AGE DISTRIBUTION 2011  ### ### ### ### ### ### D
#      # v22  <- M[62,78:88]
#     #  addlik <- hiv_prev_by_age_11_lLik(V=v22); addlik
#     #  lLik <- lLik + addlik
#       ### ### ### ART COVERAGE 2010-15  ### ### ### ### ### ###  M[61:66,"Year"];ResNam[26:29]
#       v23  <- (M[61:65,27]+M[61:65,29])/rowSums(M[61:65,26:29])
#       addlik <- art_vol_10_lLik(V=v23); addlik
#       lLik <- lLik + addlik
#       ### ### ### HOMELESS POP 2010  ### ### ### ### ### ###
#       v23b  <- M[61,31]
#       addlik <- homeless_10_lLik(V=v23b); addlik
#       lLik <- lLik + addlik
#       ### ### ### LIKELIHOOD FOR BORGDORFF, FEREBEE & SUTHERLAND ESTIMATES  ### ### ### ### ### ###
#       v2456  <- c(pfast,pimmed,rslow,rfast,rRecov)
#       addlik <- borgdorff_lLik( Par=v2456); addlik
#       lLik <- lLik + addlik
#       addlik <- ferebee_lLik(   Par=v2456); addlik
#       lLik <- lLik + addlik
#       addlik <- sutherland_lLik(Par=v2456); addlik
#       lLik <- lLik + addlik
#       ### ### ### PERCENT SMEAR-POS ### ### ### ### ### ###
#       v28  <- M[1,19]/sum(M[1,18:19])
#       addlik <- smr_pos_lLik(V=v28); addlik
#       lLik <- lLik + addlik
#       ### ### ### HIV SURVIVAL ### ### ### ### ### ###
#       v29  <- list(vHMort=vHMort,vHxtoHy=vHxtoHy,mubt00=mubt[12*50,])
#       addlik <- hiv_surv_lLik(Par=v29); addlik
#       lLik <- lLik + addlik
#
#       ### ### ###  ALL LIKELIHOODS DONE !!  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#       } }, error = function(e) NA)
#       if(is.na(jj))         { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
#       if(jj%in%c(-Inf,Inf)) { lLik <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
#
#       return((lLik))  }
#
#
#   ###################### local parallelization via multicore
#   llikelihood <- function(ParMatrix,n_cores=1) {
#     if(dim(as.data.frame(ParMatrix))[2]==1) {
#       lLik <- llikelihoodZ(1,t(as.data.frame(ParMatrix)))  } else {
#       lLik <- unlist(mclapply(1:nrow(ParMatrix),llikelihoodZ,ParMatrix=ParMatrix,mc.cores=n_cores)) }
#     return((lLik)) }
#
#
# #########################  FUNCTION2  #########################
#
#   lprior <- function(ParMatrix) { # Par = ParInit
#     if(dim(as.data.frame(ParMatrix))[2]==1) { ParMatrix <- t(as.data.frame(ParMatrix)) }
#     lPri <- rep(0,nrow(ParMatrix))
#     for(samp_i in 1:nrow(ParMatrix)) {
#       # norm2unif
#       Par <- ParMatrix[samp_i,]
#       Par2 <- pnorm(Par,0,1)
#       # unif2true
#       Par3 <- Par2
#       Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
#       Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
#       Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
#       lPri[samp_i] <- lPrior2(Par,Par3)
#       if(is.na(lPri[samp_i]))         { lPri[samp_i] <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
#       if(lPri[samp_i]%in%c(-Inf,Inf)) { lPri[samp_i] <- -10^12 - sum((ParamInitZ[,8]-Par)^2) }
#     } # end inner loop
#     return(lPri) }
#
# #########################  FUNCTION3  #########################
#
#   sample.prior1 <- function(n) { rmnorm(n,rep(0,sum(ParamInit$Calib==1)),diag(sum(ParamInit$Calib==1))) }
#   sample.prior2 <- function(n) { qnorm(randomLHS(n,sum(ParamInit$Calib==1)),0,1)   }
#   sample.prior  <- sample.prior2
#
#
#
