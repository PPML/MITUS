###  Calibration functions
  library(MCMCpack)

### WHICH STATE???
  # st = 5

########################  ########################  ########################
###########################  CALIBRATION DATA  ########################
########################  ########################  ########################
  # load("CalibDat_9-14-16.rData") # CalibDat
data("CalibDatState_7-2-18", package="MITUS")  # names(CalibDatState)

########################  ########################  ########################
###########################  HELPFUL FUNCTIONS  ########################
########################  ########################  ########################

### ### ### Dirichlet multinomial density ### ### ### ### ### ### ### ### ###
  dDirMult <- function(M,n,Rho) {
    # M= params of the dirichlet, i.e. model strain distribution
    # Rho = correlation parameter = 1/ sample size
    # n = category counts from the survey
    if(dim(as.matrix(M))[2]==1) { M <- M/sum(M) } else {  M <- M/rowSums(M)  }
    rowSums(lgamma(n+M/Rho))-rowSums(lgamma(M/Rho)) }

########################  ########################  ########################
########################  LOG_LIKELIHOOD FUNCTIONS  ########################
########################  ########################  ########################

### ### ### calibration importance weights  ### ### ### ### ### ### D
  wts <- CalibDatState[["ImptWeights"]]

### ### ### TOTAL DIAGNOSED CASES 1993-2016  ### ### ### ### ### ### D
# Motivation: norm, mean centered with CI = +/- 5% of mean
  wtZ <- wts[44:67];  wtZ["2016"] <- 4
  notif_tot     <- CalibDatState[["cases_yr_st"]][[st]][,2];
  adj_1         <- sum(dnorm(notif_tot,notif_tot,notif_tot*0.1/1.96,log=T)*wtZ)
  notif_tot_lik <- function(V) { # V = vector of total notifications 1993-2016
    sum(dnorm(notif_tot,V*1e6,notif_tot*0.1/1.96,log=T)*wtZ) - adj_1  }

  ### ### ### ANN DECLINE IN CASES 1953-1994  ### ### ### ### ### ### D
  # notif_decline      <- CalibDatState[["cases_prop_change_53_94"]]
  # adj_1a             <- sum(dnorm(rep(notif_decline,40),notif_decline,0.1,log=T))
  # notif_decline_lLik <- function(V) { # V = vector of notifications by fb 1953-2015
  #  sum(dnorm(exp(diff(log(V))),notif_decline,0.1,log=T)) - adj_1a  }

  ### ### ### TOTAL DIAGNOSED CASES 1953-1993  ### ### ### ### ### ### D
  notif_decline  <- CalibDatState[["cases_prop_change_53_94"]]
  notif_tot2     <- cumprod(notif_decline)/prod(notif_decline)*notif_tot[1]
  adj_1b         <- sum(dnorm(notif_tot2,notif_tot2,notif_tot2*0.2/1.96,log=T)*wts[4:44])
  notif_decline_lLik <- function(V) { # V = vector of total notifications 1953-1993
    sum(dnorm(notif_tot2,V*1e6,notif_tot2*0.2/1.96,log=T)*wts[4:44]) - adj_1b  }

  ### ### ### US CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_age_us0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"usb"]
  notif_age_us      <- notif_age_us0[,-c(1,12)]*notif_age_us0[,12]
  adj_2a            <- sum(dDirMult(M=notif_age_us+0.01,n=notif_age_us,Rho=0.015)*wts[44:67])
  notif_age_us_lLik <- function(V,rho=0.015) { # V = table of us notifications by age 1993-2016 (row=24 years, col=11 ages)
    V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=notif_age_us,Rho=rho)*wts[44:67]) - adj_2a  }

### ### ### FB CASES AGE DISTRIBUTION 1993-2016  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_age_fb0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"nusb"]
  notif_age_fb      <- notif_age_fb0[,-c(1,12)]*notif_age_fb0[,12]
  adj_2b            <- sum(dDirMult(M=notif_age_fb+0.01,n=notif_age_fb,Rho=0.015)*wts[44:67])
  notif_age_fb_lLik <- function(V,rho=0.015) { # V = table of fb notifications by age 1993-2013 (row=21 years, col=11 ages)
    V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=notif_age_fb,Rho=rho)*wts[44:67]) - adj_2b  }

### ### ### CASES FB DISTRIBUTION 1993-2016  ### ### ### ### ### ###  D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_fb      <- cbind(notif_age_fb0[,12],notif_age_us0[,12])
  adj_3         <- sum(dDirMult(M=notif_fb+0.01,n=notif_fb,Rho=0.005)*wts[44:67])
  notif_fb_lLik <- function(V,rho=0.005) { # V = table of notifications by fb 1993-2016 (row=24 years, col=fb then us)
    sum(dDirMult(M=V,n=notif_fb,Rho=rho)*wts[44:67]) - adj_3  }

### ### ### CASES FB DISTRIBUTION SLOPES OVER PAST 5 year  ### ### ### ### ### ### D
  notif_fbus_slp5     <- CalibDatState[["case_change_5"]][[st]];
  adj_3a              <- sum(dnorm(notif_fbus_slp5,notif_fbus_slp5,0.005,log=T))
  notif_fbus_slp_lLik <- function(V) { # V = table of notifications by fb 2011-2016 (row=6 years, col=fb then us)
    V2 <- c(mean(diff(log(V[,2]))),mean(diff(log(V[,1]))))
    sum(dnorm(V2,notif_fbus_slp5,0.005,log=T)) - adj_3a  }

### ### ### CASES TX HISTORY DISTRIBUTION IN 5-YEAR BANDS  ### ### ### ### ### ###  D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_prev0     <- CalibDatState[["prev_cases"]][[st]]
  notif_prev      <- cbind(1-notif_prev0$pct_prev,notif_prev0$pct_prev)*notif_prev0$sample_size
  adj_4           <- sum(dDirMult(M=notif_prev+0.01,n=notif_prev,Rho=0.005)*wts[c(50,55,60,65)])
  notif_prev_lLik <- function(V,rho=0.005) { # V = table of notifications by tx history (row=1997:2016, col=n then e)
    V2 <- rbind(colSums(V[1:5,]),colSums(V[6:10,]),colSums(V[11:15,]),colSums(V[16:20,]))
    sum(dDirMult(M=V2,n=notif_prev,Rho=rho)*wts[c(50,55,60,65)]) - adj_4  }

  ### ### ### CASES HIV DISTRIBUTION 1993-2014  ### ### ### ### ### ###  D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_hiv0     <- CalibDatState[["hiv_cases"]][[st]]
  rnotif_hiv1    <- rowMeans(notif_hiv0[,c("pct_hiv_hi","pct_hiv_lo")])
  notif_hiv      <- cbind(rnotif_hiv1,1-rnotif_hiv1)*notif_hiv0$sample_size
  adj_5          <- sum((dDirMult(M=notif_hiv+0.01,n=notif_hiv,Rho=0.005)*wts[c(50,55,60,65)])[is.na(notif_hiv[,1])==F])
  notif_hiv_lLik <- function(V,rho=0.005) { # V = table of notifications by tx history (row=97:16, col=pos then neg)
    V2 <- rbind(colSums(V[1:5,]),colSums(V[6:10,]),colSums(V[11:15,]),colSums(V[16:20,]))
    sum((dDirMult(M=V2,n=notif_hiv,Rho=rho)*wts[c(50,55,60,65)])[is.na(notif_hiv[,1])==F]) - adj_5  }

  ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_us_hr0     <- CalibDatState[["hr_cases"]][[st]]
  notif_us_hr      <- cbind(notif_us_hr0$pct_hr,1-notif_us_hr0$pct_hr)*notif_us_hr0$sample_size
  adj_5b           <- sum(dDirMult(M=notif_us_hr+0.01,n=notif_us_hr,Rho=0.005)*wts[c(50,55,60,65)])
  notif_us_hr_lLik <- function(V,rho=0.005) { # V = table of notifications by tx history (row=97:16, col=n then e)
    V2 <- rbind(colSums(V[1:5,]),colSums(V[6:10,]),colSums(V[11:15,]),colSums(V[16:20,]))
    sum(dDirMult(M=V2,n=notif_us_hr,Rho=rho)*wts[c(50,55,60,65)]) - adj_5b  }

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_fb_rec      <- cbind(CalibDatState[["fb_recent_cases"]][,2],1-CalibDatState[["fb_recent_cases"]][,2])*CalibDatState[["fb_recent_cases"]][,3]
  adj_6             <- sum(dDirMult(M=notif_fb_rec+0.01,n=notif_fb_rec,Rho=0.02)*wts[44:65])
  notif_fb_rec_lLik <- function(V,rho=0.02) { # V = table of notifications by rec 1993-2014 (row=22 years, col=pos then neg)
    sum(dDirMult(M=V,n=notif_fb_rec,Rho=rho)*wts[44:65]) - adj_6  }

### ### ### CASES PCT MDR RES US N 1997-2016  ### ### ### ### ### ###  D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_mdr_us_n0     <- CalibDatState[["cases_mdr_usb_npr"]][[st]]
  notif_mdr_us_n      <- cbind(notif_mdr_us_n0[,2],1-notif_mdr_us_n0[,2])*notif_mdr_us_n0[,3]
  adj_7               <- sum((dDirMult(M=notif_mdr_us_n+0.01,n=notif_mdr_us_n,Rho=0.005)*wts[c(50,55,60,65)])[notif_mdr_us_n0$sample_size>0])
  notif_mdr_us_n_lLik <- function(V,rho=0.005) { # V = table of US N notifications by dr (row=97:16, col= 5 dr cats )
    V1 <- cbind(rowSums(V[,c(4,5)]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    V2 <- rbind(colSums(V1[1:5,]),colSums(V1[6:10,]),colSums(V1[11:15,]),colSums(V1[16:20,]))
    sum((dDirMult(M=V2,n=notif_mdr_us_n,Rho=rho)*wts[c(50,55,60,65)])[notif_mdr_us_n0$sample_size>0]) - adj_7  }
  #### D
  notif_mdr_us_e0     <- CalibDatState[["cases_mdr_usb_pr"]][[st]]
  notif_mdr_us_e      <- cbind(notif_mdr_us_e0[,2],1-notif_mdr_us_e0[,2])*notif_mdr_us_e0[,3]
  adj_8               <- sum((dDirMult(M=notif_mdr_us_e+0.01,n=notif_mdr_us_e,Rho=0.005)*wts[c(50,55,60,65)])[notif_mdr_us_e0$sample_size>0])
  notif_mdr_us_e_lLik <- function(V,rho=0.005) { # V = table of US E notifications by dr (row=97_01,02_06,07_11,12_16, col= 5 dr cats )
    V1 <- cbind(rowSums(V[,c(4,5)]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    V2 <- rbind(colSums(V1[1:5,]),colSums(V1[6:10,]),colSums(V1[11:15,]),colSums(V1[16:20,]))
    sum((dDirMult(M=V2,n=notif_mdr_us_e,Rho=rho)*wts[c(50,55,60,65)])[notif_mdr_us_e0$sample_size>0]) - adj_8  }
  #### D
  notif_mdr_fb_n0     <- CalibDatState[["cases_mdr_nusb_npr"]][[st]]
  notif_mdr_fb_n      <- cbind(notif_mdr_fb_n0[,2],1-notif_mdr_fb_n0[,2])*notif_mdr_fb_n0[,3]
  adj_9               <- sum((dDirMult(M=notif_mdr_fb_n+0.01,n=notif_mdr_fb_n,Rho=0.005)*wts[c(50,55,60,65)])[notif_mdr_fb_n0$sample_size>0])
  notif_mdr_fb_n_lLik <- function(V,rho=0.005) { # V = table of US N notifications by dr (row=97_01,02_06,07_11,12_16, col= 5 dr cats )
    V1 <- cbind(rowSums(V[,c(4,5)]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    V2 <- rbind(colSums(V1[1:5,]),colSums(V1[6:10,]),colSums(V1[11:15,]),colSums(V1[16:20,]))
    sum((dDirMult(M=V2,n=notif_mdr_fb_n,Rho=rho)*wts[c(50,55,60,65)])[notif_mdr_fb_n0$sample_size>0]) - adj_9  }
  #### D
  notif_mdr_fb_e0     <- CalibDatState[["cases_mdr_nusb_pr"]][[st]]
  notif_mdr_fb_e      <- cbind(notif_mdr_fb_e0[,2],1-notif_mdr_fb_e0[,2])*notif_mdr_fb_e0[,3]
  adj_10              <- sum((dDirMult(M=notif_mdr_fb_e+0.01,n=notif_mdr_fb_e,Rho=0.005)*wts[c(50,55,60,65)])[notif_mdr_fb_e0$sample_size>0])
  notif_mdr_fb_e_lLik <- function(V,rho=0.005) { # V = table of US N notifications by dr (row=97_01,02_06,07_11,12_16, col= 5 dr cats )
    V1 <- cbind(rowSums(V[,c(4,5)]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    V2 <- rbind(colSums(V1[1:5,]),colSums(V1[6:10,]),colSums(V1[11:15,]),colSums(V1[16:20,]))
    sum((dDirMult(M=V2,n=notif_mdr_fb_e,Rho=rho)*wts[c(50,55,60,65)])[notif_mdr_fb_e0$sample_size>0]) - adj_10  }

### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  tx_outcomes      <- cbind(1-rowSums(CalibDatState[["tx_outcomes"]][,2:3]),CalibDatState[["tx_outcomes"]][,2],CalibDatState[["tx_outcomes"]][,3])*CalibDatState[["tx_outcomes"]][,4]
  adj_11           <- sum(dDirMult(M=tx_outcomes+0.01,n=tx_outcomes,Rho=0.01)*wts[44:63])
  tx_outcomes_lLik <- function(V,rho=0.01) { # V = table of treatment outcomes 1993-2012 (row=20 years, col= complete, discontinue, dead)
    sum(dDirMult(M=V,n=tx_outcomes,Rho=rho)*wts[44:63]) - adj_11  }

### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D
# Motivation: norm, mean centered with CI = +/- 10% of mean
  tltbi_vol        <- CalibDatState[["TLTBI_volume_state"]][[st]]
  adj_12           <- dnorm(tltbi_vol[1],tltbi_vol[1],diff(tltbi_vol[2:3])/1.96,log=T)
  tltbi_tot_lLik   <- function(V) { # V = total TLTBI inits in 2002 (scalar)
    dnorm(tltbi_vol[1],V*1e6,diff(tltbi_vol[2:3])/1.96,log=T) - adj_12  }

### ### ### DISTRIBUTION OF LTBI TREATMENT INITS 2002  ### ### ### ### ### ###  D
  TLTBI_dist       <- CalibDatState[["TLTBI_dist"]]
  adj_13           <- sum( dbeta(TLTBI_dist,TLTBI_dist*100,(1-TLTBI_dist)*100,log=T) )
  tltbi_dist_lLik  <- function(V) { # V = dist TLTBI inits in 2002 (vector fraction FB, HR, HV in 2002)
    sum( dbeta(TLTBI_dist,V*100,(1-V)*100,log=T) ) - adj_13  }

### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ### D
# Motivation: additional prior on LTBI, using beta densities parameterized to Miramontes/Hill results
  ltbi_us_11      <- CalibDatState[["LTBI_prev_US_11_IGRA"]]
  adj_15          <- sum( dbeta(ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) )
  ltbi_us_11_lLik <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) ) - adj_15)*2  }

  ltbi_us_11_dp      <- CalibDatState[["LTBI_prev_US_11_DoubPos"]]
  adj_15dp           <- sum( dbeta(ltbi_us_11_dp[,2]/rowSums(ltbi_us_11_dp[,2:3]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) )
  ltbi_us_11_dp_lLik <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) ) - adj_15dp)*2  }

### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ### D
# Motivation: multinomial adjusted to match effective sample size due to survey weighting
  ltbi_fb_11      <- CalibDatState[["LTBI_prev_FB_11_IGRA"]]
  adj_16          <- sum( dbeta(ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) )
  ltbi_fb_11_lLik <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) ) - adj_16)*2  }

  ltbi_fb_11_dp      <- CalibDatState[["LTBI_prev_FB_11_DoubPos"]]
  adj_16dp           <- sum( dbeta(ltbi_fb_11_dp[,2]/rowSums(ltbi_fb_11_dp[,2:3]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) )
  ltbi_fb_11_dp_lLik <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) ) - adj_16dp)*2  }

### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
# Motivation: norm, mean centered with CI = +/- 2 million wts[1+0:6*10]
  tot_pop_yr      <- CalibDatState[["tot_pop_yr_fb"]][[st]]
  tot_pop_yr_fb   <- tot_pop_yr[tot_pop_yr$usb==0,]
  adj_17          <- sum(dnorm(tot_pop_yr_fb[-1,3],tot_pop_yr_fb[-1,3],tot_pop_yr_fb[7,3]*0.05/1.96,log=T)*wts[1+1:6*10])
  tot_pop_yr_fb_lLik <- function(V) { # V = total pop (rows=year, cols=us, fb)
    sum(dnorm(tot_pop_yr_fb[-1,3],V[c(11,21,31,41,51,61)]*1e6,tot_pop_yr_fb[7,3]*0.1/1.96,log=T)*wts[1+1:6*10]) - adj_17  } # CI = +/- 2mil

  ### D
  tot_pop_yr_us  <- tot_pop_yr[tot_pop_yr$usb==1,]
  adj_17b        <- sum(dnorm(tot_pop_yr_us[6:7,3],tot_pop_yr_us[6:7,3],tot_pop_yr_us[7,3]*0.05/1.96,log=T)*wts[1+5:6*10])
  tot_pop_yr_us_lLik_00_10 <- function(V) { # V = total pop (rows=year, cols=us, fb)
    sum(dnorm(tot_pop_yr_us[6:7,3],V[c(51,61)]*1e6,tot_pop_yr_us[7,3]*0.05/1.96,log=T)*wts[1+5:6*10]) - adj_17b  } # CI = +/- 2mil

### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ### D
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
  pop_ag_11_160  <- CalibDatState[["tot_pop_ag_fb_11_16"]][[st]]
  pop_ag_11_16   <- cbind(pop_ag_11_160$pop[pop_ag_11_160$usb==1]/sum(pop_ag_11_160$pop[pop_ag_11_160$usb==1]),
                             pop_ag_11_160$pop[pop_ag_11_160$usb==0]/sum(pop_ag_11_160$pop[pop_ag_11_160$usb==0]))
  adj_18         <- sum(log(pop_ag_11_16[,1])*pop_ag_11_16[,1])+sum(log(pop_ag_11_16[,2])*pop_ag_11_16[,2])
  tot_pop14_ag_fb_lLik <- function(V,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
    V1 <- rbind(V[1:9,],V[10,]+V[11,])
    (sum(log(V1[,1]/sum(V1[,1]))*pop_ag_11_16[,1])+sum(log(V1[,2]/sum(V1[,2]))*pop_ag_11_16[,2]))*ESS - adj_18*ESS  }

  ### ### ### Total TB DEATHS 1999-2016 ### ### ### ### ### ### D
  tb_deaths <- CalibDatState[["tbdeaths"]][[st]]$Deaths
  adj_19    <- sum((dnorm(tb_deaths,tb_deaths,tb_deaths*0.1/1.96,log=T)*wts[50:67])[is.na(tb_deaths)==F])
  tbdeaths_lik <- function(V) { # V = vector of total notifications 1999-2016
    sum((dnorm(tb_deaths,rowSums(V)*1e6,tb_deaths*0.2/1.96,log=T)*wts[50:67])[is.na(tb_deaths)==F]) - adj_19 }
  ### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D
  tbdeaths_decline      <- CalibDatState[["deaths_ann_decline_68_15"]]
  adj_19a            <- sum(dnorm(tbdeaths_decline,tbdeaths_decline,0.005,log=T))
  tbdeaths_decline_lLik <- function(V) { # V = vector of tb deaths 1968-2015
    V2 <- (1-(V[48]/V[1])^(1/47))
    sum(dnorm(V2,tbdeaths_decline,0.015/1.96,log=T)) - adj_19a  }
  ### ### ### TB DEATHS AGE DISTRIBUTION 1999-2016  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  tb_deaths_age  <- CalibDatState[["tbdeaths_age_yr"]][,-1]
  adj_19b        <- sum(dDirMult(M=tb_deaths_age+0.01,n=tb_deaths_age,Rho=0.01)*wts[50:67])
  tb_dth_age_lLik <- function(V,rho=0.01) { # V = table of deaths by age 1999-2016 (row=18 years, col=11 ages)
    V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=tb_deaths_age,Rho=rho)*wts[50:67]) - adj_19b  }

  ### ### ### TOTAL HIV DEATHS 2008-2015  ### ### ### ### ### ### D
  # Motivation: norm, mean centered with CI = +/- 20% of mean
  hiv_deaths0      <- CalibDatState[["hiv_deaths"]][[st]]
  hiv_deaths <- aggregate(deaths~year,hiv_deaths0,sum)$deaths
  adj_20    <- sum(dnorm(hiv_deaths,hiv_deaths,hiv_deaths*0.25/1.96,log=T)*wts[59:66])
  hivdeaths_lik <- function(V) { # V = vector of HIV DEATHS 2008-2015
    sum(dnorm(hiv_deaths,V*1e6,hiv_deaths*0.25/1.96,log=T)*wts[59:66]) - adj_20 }

  ### ### ### HIV DEATHS AGE DISTRIBUTION 1999-2015  ### ### ### ### ### ### D
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  hv_deaths_age  <- CalibDatState[["hiv_deaths_age"]][[st]][,-1]
  adj_20b        <- sum(dDirMult(M=hv_deaths_age+0.01,n=hv_deaths_age,Rho=0.02)*wts[59:66])
  hv_dth_age_lLik <- function(V,rho=0.01) { # V = table of deaths by age 2008-2015 (row=16 years, col=11 ages)
    V2 <- cbind(V[,3:6],rowSums(V[,7:11]))
    sum(dDirMult(M=V2,n=hv_deaths_age,Rho=rho)*wts[59:66]) - adj_20b  }

### ### ### HIV PREVALENCE 2010-2015 ### ### ### ### ### ### D
# Motivation: norm, mean centered with CI = 5% million
  hiv_prev     <- CalibDatState[["hiv_prev"]][[st]]$prev/1e5
  adj_21                <-  sum(dnorm(hiv_prev,hiv_prev,hiv_prev*0.1/1.96,log=T)*wts[61:66])
  hiv_prev_by_year_lLik <- function(V) { # V = HIV prev >15 yo, 2010-2015
    sum(dnorm(hiv_prev,V,hiv_prev*0.1/1.96,log=T)*wts[61:66]) - adj_21  }

### ### ### ART VOLUME 2010  ### ### ### ### ### ### D
# Motivation: norm, mean centered with CI = +/- 5% of mean
  art_cov      <- CalibDatState[["hiv_art"]][[st]]$Percent/100
  adj_23          <- sum((dnorm(art_cov,art_cov,art_cov*0.10/1.96,log=T))[is.na(art_cov)==F])
  art_vol_10_lLik <- function(V) { # V = ART coverage as fraction, 2010-2015
    sum((dnorm(art_cov,V,art_cov*0.1/1.96,log=T))[is.na(art_cov)==F]) - adj_23   }

  ### ### ### HOMELESS POP 2010  ### ### ### ### ### ### names(CalibDatState)
  # Motivation: norm, mean centered with CI = +/- 25% of mean
  homeless_pop      <- CalibDatState[["homeless_pop"]][[st]][1]
  adj_23b          <- dnorm(homeless_pop,homeless_pop,homeless_pop*0.2/1.96,log=T)
  homeless_10_lLik <- function(V) { # V = homeless pop in 2010 (scalar)
    dnorm(homeless_pop,V,homeless_pop*0.2/1.96,log=T) - adj_23b   }

### ### ### LIKELIHOOD FOR BORGDORFF ESTIMATES  ### ### ### ### ### ###
  ss_borgdorff   <- 854.5463/4
  datB           <- CalibDatState[["borgdorff_data"]]
  adj_24         <- sum(diff(-datB[,2])*ss_borgdorff*log(diff(-datB[,2])) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-datB[,2])))
  borgdorff_lLik <- function(Par,N_red=1) { # Par = c(pfast,pimmed,rslow,rfast)
    zz <- tryCatch({
      pfast <- Par[1]; pimmed <- Par[2]; rslow <- Par[3]*12;  rfast <- Par[4]*12;  rRecov <- Par[5]*12;
      p <- pfast*pimmed*c(1,rep(0,nrow(datB)-1)) +
           pfast*(1-pimmed)*(1-(1-exp(-(rfast+rRecov)*datB[,1]))*(rfast/(rfast+rRecov))) +
           (1-pfast)*(1-(1-exp(-(rslow+rRecov)*datB[,1]))*(rslow/(rslow+rRecov)))
      p <- 1-(1-p)/(1-p)[nrow(datB)]
      sum(diff(-datB[,2])*ss_borgdorff*log(diff(-p)) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-p)))/N_red - adj_24/N_red
    },error=function(e) -Inf )
    if(is.nan(zz)) { zz = -10^4 } else { if(zz== -Inf) zz = -10^4 }
    zz }

### ### ### LIKELIHOOD FOR FEREBEE ESTIMATES  ### ### ### ### ### ###
  datF         <- CalibDatState[["ferebee_data"]]
  adj_25       <- sum(datF[,3]*log(datF[,3]/datF[,2]) + (datF[,2]-datF[,3])*log(1-datF[,3]/datF[,2]))
  n_yr_F       <- nrow(datF)
  ferebee_lLik <- function(Par,N_red=4) { # Par = c(pfast,pimmed,rslow,rfast)
    zz <- tryCatch({
      pfast <- Par[1]; pimmed <- Par[2]; rslow <- Par[3]*12;  rfast <- Par[4]*12;  rRecov <- Par[5]*12;
      p2 <- pfast*pimmed*c(1,rep(0,n_yr_F)) +
            pfast*(1-pimmed)*(1-(1-exp(-(rfast+rRecov)*(0:n_yr_F)))*(rfast/(rfast+rRecov))) +
            (1-pfast)*(1-(1-exp(-(rslow+rRecov)*(0:n_yr_F)))*(rslow/(rslow+rRecov)))
      r2 <- -log(1-diff(-p2)/p2[-(n_yr_F+1)])/1
      sum((datF[,3]*log(r2) + (datF[,2]-datF[,3])*log(1-r2))/N_red)-adj_25/N_red
    },error=function(e) -Inf )
    if(is.nan(zz)) { zz = -10^4 } else { if(zz== -Inf) zz = -10^4 }
    zz  }

### ### ### LIKELIHOOD FOR SUTHERLAND ESTIMATES  ### ### ### ### ### ###
  datS            <- CalibDatState[["sutherland_data"]]
  datSz           <- datS; datSz[datS[,3]==0,3] <- 0.01
  adj_26          <- sum(datSz[,3]*log(datSz[,3]/datSz[,2]) + (datSz[,2]-datSz[,3])*log(1-datSz[,3]/datSz[,2]))
  n_yr_S          <- nrow(datS)
  sutherland_lLik <- function(Par,N_red=4) {
    zz <- tryCatch({
      pfast <- Par[1]; pimmed <- Par[2]; rslow <- Par[3]*12;  rfast <- Par[4]*12;  rRecov <- Par[5]*12;
      p2 <- pfast*pimmed*c(1,rep(0,n_yr_S)) +
        pfast*(1-pimmed)*(1-(1-exp(-(rfast+rRecov)*(0:n_yr_S)))*(rfast/(rfast+rRecov))) +
        (1-pfast)*(1-(1-exp(-(rslow+rRecov)*(0:n_yr_S)))*(rslow/(rslow+rRecov)))
      r2 <- -log(1-diff(-p2)/p2[-(n_yr_S+1)])/1
      sum((datS[,3]*log(r2) + (datS[,2]-datS[,3])*log(1-r2))/N_red) - adj_26/N_red
    },error=function(e) -Inf )
    if(is.nan(zz)) { zz = -10^4 } else { if(zz== -Inf) zz = -10^4 }
    zz  }

### ### ### LIKELIHOOD FOR TIEMERSMA ESTS ### ### ### ### ### ###
  adj_27         <- dnorm(3.0,3.0,0.5/1.96,log=T)+dbeta(0.2,0.2*50,(1-0.2)*50,log=T)+dbeta(0.7,0.7*50,(1-0.7)*50,log=T)
  tiemersma_lLik <- function(Par) { # Par= c(pSmPos,rSlfCur,muIn,muIp)
    pSmPos <- Par[1]; rSlfCur <- Par[2]*12; muIn <- Par[3]*12;  muIp <- Par[4]*12;
    dur <- pSmPos*(1/(rSlfCur+muIp)) + (1-pSmPos)*(1/(rSlfCur+muIn))
    cf_sp <- muIp/(rSlfCur+muIp);  cf_sn <- muIn/(rSlfCur+muIn)
    l1 <- dnorm(3.0,dur,0.5/1.96,log=T)
    l2 <- dbeta(0.2,cf_sn*50,(1-cf_sn)*50,log=T)
    l3 <- dbeta(0.7,cf_sp*50,(1-cf_sp)*50,log=T)
    l1+l2+l3 - adj_27  }

### ### ### PERCENT SMEAR-POS ### ### ### ### ### ###
  adj_28       <- dbeta(CalibDatState[["pct_sp_ss"]][1],CalibDatState[["pct_sp_ss"]][1]*100,(1-CalibDatState[["pct_sp_ss"]][1])*100,log=T)
  smr_pos_lLik <- function(V) { # V = fraction smr-pos in 1950 (scalar)
    dbeta(CalibDatState[["pct_sp_ss"]][1],V*100,(1-V)*100,log=T) - adj_28  }

### ### ### HIV SURVIVAL ### ### ### ### ### ###
  hiv_surv <- CalibDatState[["hiv_survival"]]
  adj_29  <- sum(dnorm(hiv_surv[,2],hiv_surv[,2],(hiv_surv[,4]-hiv_surv[,3])/1.96*2,log=T))
  hiv_surv_lLik <- function(Par) { # Par= c(muH2,mHxtoHy[1])
    sv <- (1/(Par[[1]][,2]+Par[[2]][,1]+Par[[3]])+Par[[2]][,1]/(Par[[1]][,2]+Par[[2]][,1]+Par[[3]])*1/(Par[[1]][,4]+Par[[3]]))[3:8]/12
    sum(dnorm(sv,hiv_surv[,2],(hiv_surv[,4]-hiv_surv[,3])/1.96*2,log=T)) - adj_29  }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


