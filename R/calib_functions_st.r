#'Calibration Functions for use with the tb_model.cpp model for state level
#'This script creates several individual log likelihood functions
#'for the calibration of the State Level TB model in tb_model.cpp
#'These llikelihood functions are called in IMIS_functions.R
#'takes in the outputs and calibration data and creates likelihood functions

  library(MCMCpack)

  #'Total Diagnosed Cases 1953-2016
  #'Motivation: Normal, mean centered with CI = +/- 5% of the mean
  #'@name notif_tot_lLik_st
  #'@param V vector of total notifications 1953-2014
  #'@return likelihood
  notif_tot_lLik_st <- function(V,st) { # V = vector of total notifications 1993-2016
    notif_tot     <- CalibDatState[["cases_yr_st"]][[st]][,2];
    adj_1         <- sum(dnorm(notif_tot,notif_tot,notif_tot*0.1/1.96,log=T)*wtZ)
    sum(dnorm(notif_tot,V*1e6,notif_tot*0.1/1.96,log=T)*wtZ) - adj_1  }

  ### ### ### ANN DECLINE IN CASES 1953-1994  ### ### ### ### ### ### D
  # notif_decline      <- CalibDatState[["cases_prop_change_53_94"]]
  # adj_1a             <- sum(dnorm(rep(notif_decline,40),notif_decline,0.1,log=T))
  # notif_decline_lLik_st <- function(V) { # V = vector of notifications by fb 1953-2015
  #  sum(dnorm(exp(diff(log(V))),notif_decline,0.1,log=T)) - adj_1a  }

  ### ### ### TOTAL DIAGNOSED CASES 1953-1993  ### ### ### ### ### ### D

  notif_decline_lLik_st <- function(V, st=st) {
    notif_tot     <- CalibDatState[["cases_yr_st"]][[st]][,2];
    # V = vector of total notifications 1953-1993
    notif_decline  <- CalibDatState[["cases_prop_change_53_94"]]
    notif_tot2     <- cumprod(notif_decline)/prod(notif_decline)*notif_tot[1]
    adj_1b         <- sum(dnorm(notif_tot2,notif_tot2,notif_tot2*0.2/1.96,log=T)*wts[4:44])
    sum(dnorm(notif_tot2,V*1e6,notif_tot2*0.2/1.96,log=T)*wts[4:44]) - adj_1b  }

  ### ### ### US CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  notif_age_us_lLik_st <- function(V,st,rho=0.015) { # V = table of us notifications by age 1993-2016 (row=24 years, col=11 ages)
    notif_age_us0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"usb"]
    notif_age_us      <- notif_age_us0[,-c(1,12)]*notif_age_us0[,12]
    adj_2a            <- sum(dDirMult(M=notif_age_us+0.01,n=notif_age_us,Rho=0.015)*wts[44:67])
     V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=notif_age_us,Rho=rho)*wts[44:67]) - adj_2a  }

### ### ### FB CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  notif_age_fb_lLik_st <- function(V,st,rho=0.015) { # V = table of fb notifications by age 1993-2013 (row=21 years, col=11 ages)
    notif_age_fb0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"nusb"]
    notif_age_fb      <- notif_age_fb0[,-c(1,12)]*notif_age_fb0[,12]
    adj_2b            <- sum(dDirMult(M=notif_age_fb+0.01,n=notif_age_fb,Rho=0.015)*wts[44:67])
     V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=notif_age_fb,Rho=rho)*wts[44:67]) - adj_2b  }

### ### ### CASES FB DISTRIBUTION 1993-2016  ### ### ### ### ### ###  D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  notif_fb_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by fb 1993-2016 (row=24 years, col=fb then us)
    notif_age_fb0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"nusb"]
    notif_age_us0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"usb"]

    notif_fb      <- cbind(notif_age_fb0[,12],notif_age_us0[,12])
    adj_3         <- sum(dDirMult(M=notif_fb+0.01,n=notif_fb,Rho=0.005)*wts[44:67])
       sum(dDirMult(M=V,n=notif_fb,Rho=rho)*wts[44:67]) - adj_3  }

### ### ### CASES FB DISTRIBUTION SLOPES OVER PAST 5 year  ### ### ### ### ### ### D

  notif_fbus_slp_lLik_st <- function(V,st) {
    notif_fbus_slp5     <- CalibDatState[["case_change_5"]][[st]];
    adj_3a              <- sum(dnorm(notif_fbus_slp5,notif_fbus_slp5,0.005,log=T))# V = table of notifications by fb 2011-2016 (row=6 years, col=fb then us)
    V2 <- c(mean(diff(log(V[,2]))),mean(diff(log(V[,1]))))
    sum(dnorm(V2,notif_fbus_slp5,0.005,log=T)) - adj_3a  }

  ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  notif_us_hr_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by tx history (row=97:16, col=n then e)
    notif_us_hr0     <- CalibDatState[["hr_cases"]][[st]]
    notif_us_hr      <- cbind(notif_us_hr0$pct_hr,1-notif_us_hr0$pct_hr)*notif_us_hr0$sample_size
    adj_5b           <- sum(dDirMult(M=notif_us_hr+0.01,n=notif_us_hr,Rho=0.005)*wts[c(50,55,60,65)])
    V2 <- rbind(colSums(V[1:5,]),colSums(V[6:10,]),colSums(V[11:15,]),colSums(V[16:20,]))
    sum(dDirMult(M=V2,n=notif_us_hr,Rho=rho)*wts[c(50,55,60,65)]) - adj_5b  }

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  notif_fb_rec_lLik_st <- function(V,rho=0.02) { # V = table of notifications by rec 1993-2014 (row=22 years, col=pos then neg)
    notif_fb_rec      <- cbind(CalibDatState[["fb_recent_cases"]][,2],1-CalibDatState[["fb_recent_cases"]][,2])*CalibDatState[["fb_recent_cases"]][,3]
    adj_6             <- sum(dDirMult(M=notif_fb_rec+0.01,n=notif_fb_rec,Rho=0.02)*wts[44:65])
     sum(dDirMult(M=V,n=notif_fb_rec,Rho=rho)*wts[44:65]) - adj_6  }

### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  tx_outcomes_lLik_st <- function(V,rho=0.01) {
    tx_outcomes      <- cbind(1-rowSums(CalibDatState[["tx_outcomes"]][,2:3]),CalibDatState[["tx_outcomes"]][,2],CalibDatState[["tx_outcomes"]][,3])*CalibDatState[["tx_outcomes"]][,4]
    adj_11           <- sum(dDirMult(M=tx_outcomes+0.01,n=tx_outcomes,Rho=0.01)*wts[44:63])# V = table of treatment outcomes 1993-2012 (row=20 years, col= complete, discontinue, dead)
    sum(dDirMult(M=V,n=tx_outcomes,Rho=rho)*wts[44:63]) - adj_11  }

### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D

  tltbi_tot_lLik_st   <- function(V,st) { # V = total TLTBI inits in 2002 (scalar)
    # Motivation: norm, mean centered with CI = +/- 10% of mean
    tltbi_vol        <- CalibDatState[["TLTBI_volume_state"]][[st]]
    adj_12           <- dnorm(tltbi_vol[1],tltbi_vol[1],diff(tltbi_vol[2:3])/1.96,log=T)
        dnorm(tltbi_vol[1],V*1e6,diff(tltbi_vol[2:3])/1.96,log=T) - adj_12  }

### ### ### DISTRIBUTION OF LTBI TREATMENT INITS 2002  ### ### ### ### ### ###  D

  tltbi_dist_lLik_st  <- function(V) {
    TLTBI_dist       <- CalibDatState[["TLTBI_dist"]]
    adj_13           <- sum( dbeta(TLTBI_dist,TLTBI_dist*100,(1-TLTBI_dist)*100,log=T) )# V = dist TLTBI inits in 2002 (vector fraction FB, HR, HV in 2002)
    sum( dbeta(TLTBI_dist,V*100,(1-V)*100,log=T) ) - adj_13  }

### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ### D
# Motivation: additional prior on LTBI, using beta densities parameterized to Miramontes/Hill results

  ltbi_us_11_lLik_st <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    ltbi_us_11      <- CalibDatState[["LTBI_prev_US_11_IGRA"]]
    adj_15          <- sum( dbeta(ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) )
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) ) - adj_15)*2  }

  ltbi_us_11_dp_lLik_st <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    ltbi_us_11_dp      <- CalibDatState[["LTBI_prev_US_11_DoubPos"]]
    adj_15dp           <- sum( dbeta(ltbi_us_11_dp[,2]/rowSums(ltbi_us_11_dp[,2:3]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) )

        V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) ) - adj_15dp)*2  }

### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ### D
# Motivation: multinomial adjusted to match effective sample size due to survey weighting
  ltbi_fb_11_lLik_st <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)

    ltbi_fb_11      <- CalibDatState[["LTBI_prev_FB_11_IGRA"]]
    adj_16          <- sum( dbeta(ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) )
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) ) - adj_16)*2  }


  ltbi_fb_11_dp_lLik_st <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    ltbi_fb_11_dp      <- CalibDatState[["LTBI_prev_FB_11_DoubPos"]]
    adj_16dp           <- sum( dbeta(ltbi_fb_11_dp[,2]/rowSums(ltbi_fb_11_dp[,2:3]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) )
      V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) ) - adj_16dp)*2  }

### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
# Motivation: norm, mean centered with CI = +/- 2 million wts[1+0:6*10]

  tot_pop_yr_fb_lLik_st <- function(V,st) { # V = total pop (rows=year, cols=us, fb)
    tot_pop_yr      <- CalibDatState[["tot_pop_yr_fb"]][[st]]
    tot_pop_yr_fb   <- tot_pop_yr[tot_pop_yr$usb==0,]
    adj_17          <- sum(dnorm(tot_pop_yr_fb[-1,3],tot_pop_yr_fb[-1,3],tot_pop_yr_fb[7,3]*0.05/1.96,log=T)*wts[1+1:6*10])
        sum(dnorm(tot_pop_yr_fb[-1,3],V[c(11,21,31,41,51,61)]*1e6,tot_pop_yr_fb[7,3]*0.1/1.96,log=T)*wts[1+1:6*10]) - adj_17  } # CI = +/- 2mil

  ### D
  tot_pop_yr_us_lLik_st_00_10 <- function(V,st) {
    tot_pop_yr      <- CalibDatState[["tot_pop_yr_fb"]][[st]]
    # V = total pop (rows=year, cols=us, fb)
    tot_pop_yr_us  <- tot_pop_yr[tot_pop_yr$usb==1,]
    adj_17b        <- sum(dnorm(tot_pop_yr_us[6:7,3],tot_pop_yr_us[6:7,3],tot_pop_yr_us[7,3]*0.05/1.96,log=T)*wts[1+5:6*10])
    sum(dnorm(tot_pop_yr_us[6:7,3],V[c(51,61)]*1e6,tot_pop_yr_us[7,3]*0.05/1.96,log=T)*wts[1+5:6*10]) - adj_17b  } # CI = +/- 2mil

### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ### D
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500

  tot_pop14_ag_fb_lLik_st <- function(V,st,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
    pop_ag_11_160  <- CalibDatState[["tot_pop_ag_fb_11_16"]][[st]]
    pop_ag_11_16   <- cbind(pop_ag_11_160$pop[pop_ag_11_160$usb==1]/sum(pop_ag_11_160$pop[pop_ag_11_160$usb==1]),
                            pop_ag_11_160$pop[pop_ag_11_160$usb==0]/sum(pop_ag_11_160$pop[pop_ag_11_160$usb==0]))
    adj_18         <- sum(log(pop_ag_11_16[,1])*pop_ag_11_16[,1])+sum(log(pop_ag_11_16[,2])*pop_ag_11_16[,2])
        V1 <- rbind(V[1:9,],V[10,]+V[11,])
    (sum(log(V1[,1]/sum(V1[,1]))*pop_ag_11_16[,1])+sum(log(V1[,2]/sum(V1[,2]))*pop_ag_11_16[,2]))*ESS - adj_18*ESS  }

  ### ### ### Total TB DEATHS 1999-2016 ### ### ### ### ### ### D

  tbdeaths_lLik_st <- function(V,st) { # V = vector of total notifications 1999-2016
    tb_deaths <- CalibDatState[["tbdeaths"]][[st]]$Deaths
    adj_19    <- sum((dnorm(tb_deaths,tb_deaths,tb_deaths*0.1/1.96,log=T)*wts[50:67])[is.na(tb_deaths)==F])
    sum((dnorm(tb_deaths,rowSums(V)*1e6,tb_deaths*0.2/1.96,log=T)*wts[50:67])[is.na(tb_deaths)==F]) - adj_19 }
  ### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D

  tbdeaths_decline_lLik_st <- function(V) { # V = vector of tb deaths 1968-2015
    tbdeaths_decline      <- CalibDatState[["deaths_ann_decline_68_15"]]
    adj_19a            <- sum(dnorm(tbdeaths_decline,tbdeaths_decline,0.005,log=T))
    V2 <- (1-(V[48]/V[1])^(1/47))
    sum(dnorm(V2,tbdeaths_decline,0.015/1.96,log=T)) - adj_19a  }
  ### ### ### TB DEATHS AGE DISTRIBUTION 1999-2016  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

  tb_dth_age_lLik_st <- function(V,rho=0.01) { # V = table of deaths by age 1999-2016 (row=18 years, col=11 ages)
    tb_deaths_age  <- CalibDatState[["tbdeaths_age_yr"]][,-1]
    adj_19b        <- sum(dDirMult(M=tb_deaths_age+0.01,n=tb_deaths_age,Rho=0.01)*wts[50:67])
     V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=tb_deaths_age,Rho=rho)*wts[50:67]) - adj_19b  }

  #' TOTAL US DEATHS
  #' 1970,1975,1980,1985,1990-2007
  #' Motivation: norm, mean centered with CI = +/- 5% of mean
  #'@name dth_tot_lLik_st
  #'@param V
  #'@return likelihood

  dth_tot_lLik_st <- function(V,st) {
    data("ST_tot_mort",package="MITUS")
    ST_deaths_tot   <- ST_tot_mort[,-1]
    ST_deaths_tot   <- ST_deaths_tot[((st-1)*38)+(1:38),3]
    adj_20a         <- sum(dnorm(ST_deaths_tot,ST_deaths_tot,ST_deaths_tot*0.1/1.96,log=T)*wts[30:67])
    sum(dnorm(ST_deaths_tot,V,ST_deaths_tot*0.1/1.96,log=T)*wts[30:67]) - adj_20a
  }

  #'  #' TOTAL DEATHS AGE DISTRIBUTION 1999-2014
  #' #' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  #' #'@param V table of deaths by age 1999-2014 (row=16 years, col=11 ages)
  #' #'@param rho correlation parameter
  #' #'@return likelihood
  #' tot_dth_age_lLik_st <- function(V,rho=0.01) {
  #'   # CalibDat$US_mort_age <- read.csv(system.file("extdata","US_mort_age.csv", package="MITUS"))
  #'   tot_deaths_age  <- CalibDat[["ST_mort_age"]][,-1]
  #'   adj_20b        <- sum(dDirMult(M=tot_deaths_age+0.1,n=tot_deaths_age+0.1,Rho=0.01)*wts[50:67])
  #'   V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  #'   V2 <- V2[,-5]; V2[,4]  <- V2[,4]+V[,5]
  #'   V2 <- V2[,-3]; V2[,2]  <- V2[,2]+V[,3]
  #'
  #'   sum(dDirMult(M=(V2*1e6),n=tot_deaths_age+.1,Rho=rho)*wts[50:67]) - adj_20b  }


  #' Mortality Risk Group Distribution 1999-2014
  #' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  #'@name mort_dist_lLik_st
  #'@param V table of mort_dist 1999-2014 (row=16 years, col=11 ages)
  #'@param rho correlation parameter
  #'@return likelihood
  mort_dist_lLik_st <- function(V,rho=0.01) {
    md     <- rowSums(dist_gen)
    mort_dist     <-matrix(md,17,4, byrow = TRUE)
    adj_21        <- sum(dDirMult(M=mort_dist,n=mort_dist,Rho=0.01)*wts[51:67])
    tot_lik<-0
    for(ag in 1:11){
      V1<-V[,(1:4)+4*(ag-1)]
      x<-sum(dDirMult(M=V1,n=mort_dist,Rho=rho)*wts[51:67]) - adj_21
      tot_lik<-tot_lik+x
    }
    return(tot_lik)
  }
  ### ### ### HOMELESS POP 2010  ### ### ### ### ### ### names(CalibDatState)
  # Motivation: norm, mean centered with CI = +/- 25% of mean

  homeless_10_lLik_st <- function(V,st) { # V = homeless pop in 2010 (scalar)
    homeless_pop      <- CalibDatState[["homeless_pop"]][[st]][1]
    adj_23b          <- dnorm(homeless_pop,homeless_pop,homeless_pop*0.2/1.96,log=T)
        dnorm(homeless_pop,V,homeless_pop*0.2/1.96,log=T) - adj_23b   }

### ### ### LIKELIHOOD FOR BORGDORFF ESTIMATES  ### ### ### ### ### ###
  borgdorff_lLik_st <- function(Par_list,N_red=1) {  # Par_list = list(Mpfast[,c(1,3,2,4)], Mrslow[,c(1,3,2,4)], rfast, rRecov)
    ss_borgdorff   <- 854.5463/4
    datB           <- CalibDatState[["borgdorff_data"]]
    adj_24         <- sum(diff(-datB[,2])*ss_borgdorff*log(diff(-datB[,2])) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-datB[,2])))
    zz <- tryCatch({
      Mpfast     <- Par_list[[1]];
      # ORpfastRF <- Par[2];
      Mrslow     <- Par_list[[2]]*12;
      # RRrSlowRF <- Par[4];
      rfast     <- Par_list[[3]]*12;
      rRecov    <- Par_list[[4]]*12;
      pfast_v <- rslow_v <- rep(NA,4)
      pfast_v <- Mpfast[3,]
      rslow_v <- Mrslow[3,]
      p0 <- matrix(NA,4,10)  # this added
      for (i in 1:4){
        p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*datB[,1]))*(rfast/(rfast+rRecov))) +
          (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*datB[,1]))*(rslow_v[i]/(rslow_v[i]+rRecov)))
      }
      p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
      p <- 1-(1-p1)/(1-p1)[nrow(datB)]
      sum(diff(-datB[,2])*ss_borgdorff*log(diff(-p)) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-p)))/N_red - adj_24/N_red
    },error=function(e) -Inf )
    if(is.nan(zz)) { zz = -10^4
    } else {
      if(zz== -Inf) zz = -10^4 }
    zz }

### ### ### LIKELIHOOD FOR FEREBEE ESTIMATES  ### ### ### ### ### ###
  ferebee_lLik_st <- function(Par_list,N_red=4) {
    datF         <- CalibDatState[["ferebee_data"]]
    adj_25       <- sum(datF[,3]*log(datF[,3]/datF[,2]) + (datF[,2]-datF[,3])*log(1-datF[,3]/datF[,2]))
    n_yr_F       <- nrow(datF)
    zz <- tryCatch({

      Mpfast     <- Par_list[[1]];
      # ORpfastRF <- Par[2];
      Mrslow     <- Par_list[[2]]*12;
      # RRrSlowRF <- Par[4];
      rfast     <- Par_list[[3]]*12;
      rRecov    <- Par_list[[4]]*12;
      pfast_v <- rslow_v <- rep(NA,4)
      pfast_v <- Mpfast[3,]
      rslow_v <- Mrslow[3,]
      p0 <- matrix(NA,4,11)  # this added
      for (i in 1:4){
        p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*(0:n_yr_F)))*(rfast/(rfast+rRecov))) +
          (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*(0:n_yr_F)))*(rslow_v[i]/(rslow_v[i]+rRecov)))
      }
      p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
      r2 <- -log(1-diff(-p1)/p1[-(n_yr_F+1)])/1
      sum((datF[,3]*log(r2) + (datF[,2]-datF[,3])*log(1-r2)))/N_red - adj_25/N_red
    },error=function(e) -Inf )
    if(is.nan(zz)) {
      zz = -10^4
    } else {
      if(zz== -Inf) zz = -10^4 }
    zz  }


### ### ### LIKELIHOOD FOR SUTHERLAND ESTIMATES  ### ### ### ### ### ###
  sutherland_lLik_st <- function(Par_list,N_red=4) {
    datS            <- CalibDatState[["sutherland_data"]]
    datSz           <- datS; datSz[datS[,3]==0,3] <- 0.01
    adj_26          <- sum(datSz[,3]*log(datSz[,3]/datSz[,2]) + (datSz[,2]-datSz[,3])*log(1-datSz[,3]/datSz[,2]))
    n_yr_S          <- nrow(datS)
    zz <- tryCatch({
      Mpfast     <- Par_list[[1]];
      # ORpfastRF <- Par[2];
      Mrslow     <- Par_list[[2]]*12;
      # RRrSlowRF <- Par[4];
      rfast     <- Par_list[[3]]*12;
      rRecov    <- Par_list[[4]]*12;
      s<-pfast_v <- rslow_v <- rep(NA,4)
      pfast_v <- Mpfast[3,]
      rslow_v <- Mrslow[3,]
      p0 <- matrix(NA,4,16)  # this added
      for (i in 1:4){
        p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*(0:n_yr_S)))*(rfast/(rfast+rRecov))) +
          (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*(0:n_yr_S)))*(rslow_v[i]/(rslow_v[i]+rRecov)))
      }
      p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
      r2 <- -log(1-diff(-p1)/p1[-(n_yr_S+1)])/1
      sum((datS[,3]*log(r2) + (datS[,2]-datS[,3])*log(1-r2)))/N_red - adj_26/N_red
    },error=function(e) -Inf )
    if(is.nan(zz)) {
      zz = -10^4
    } else {
      if(zz== -Inf) zz = -10^4 }
    zz  }


### ### ### LIKELIHOOD FOR TIEMERSMA ESTS ### ### ### ### ### ###
  tiemersma_lLik_st <- function(Par) { # Par= c(rSlfCur,muIp)
    adj_27         <- dnorm(3.0,3.0,0.5/1.96,log=T)+dbeta(0.9,0.9*50,(1-0.9)*50,log=T)
    rSlfCur <- Par[1]*12;
    muIp    <- Par[2]*12 #muIp is the value for all TB active disease (not just smear pos)

    dur <- (1/(rSlfCur+muIp))
    cf_sp <- muIp/(rSlfCur+muIp);

    l1 <- dnorm(3.0,dur,0.5/1.96,log=T)
    l2 <- dbeta(0.9,cf_sp*50,(1-cf_sp)*50,log=T)
    l1+l2 - adj_27
  }


