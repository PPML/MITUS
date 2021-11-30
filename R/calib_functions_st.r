#'Calibration Functions for use with the tb_model.cpp model for state level
#'This script creates several individual log likelihood functions
#'for the calibration of the State Level TB model in tb_model.cpp
#'These llikelihood functions are called in IMIS_functions.R
#'takes in the outputs and calibration data and creates likelihood functions

#'Total Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@name notif_tot_lLik_st
#'@param V vector of total notifi
#'cations 1953-2014
#'@return likelihood
notif_tot_lLik_st <- function(V,st) {
  notif_tot     <- CalibDatState[["cases_yr_st"]][[st]][,2];
  wts2 <- wts; wts2[length(wts2)] <- wts2[length(wts2)]*2;
  adj_1         <- sum((dnorm(notif_tot,notif_tot,notif_tot*0.1/1.96,log=T)*wts2[44:70])[is.na(notif_tot)==F  & notif_tot != 0])
  #notif tot is in real scale must scale outputs up
  sum((dnorm(notif_tot,V*1e6,notif_tot*0.1/1.96,log=T)*wts2[44:70])[is.na(notif_tot)==F  & notif_tot != 0]) - adj_1
}

### ### ### TOTAL DIAGNOSED CASES 1953-1993  ### ### ### ### ### ### D
notif_decline_lLik_st <- function(V, st=st) {
  notif_tot     <- CalibDatState[["cases_yr_st"]][[st]][,2];
  # V = vector of total notifications 1953-1993
  notif_decline  <- CalibDatState[["cases_prop_change_53_94"]]
  notif_tot2     <- cumprod(notif_decline)/prod(notif_decline)*notif_tot[1]
  adj_1b         <- sum(dnorm(notif_tot2,notif_tot2,notif_tot2*0.2/1.96,log=T)*wts[4:44])
  #notif tot is in real scale must scale outputs up
  sum(dnorm(notif_tot2,V*1e6,notif_tot2*0.2/1.96,log=T)*wts[4:44]) - adj_1b
}

### ### ### CASES FB DISTRIBUTION 1993-2016  ### ### ### ### ### ###  D
# Motivatioµn: dirichlet-multinomial, multinomial data with additional non-sampling biases
notif_fb_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by fb 1993-2016 (row=24 years, col=fb then us)
  notif_age_fb0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"nusb"]
  notif_age_us0     <- CalibDatState[["cases_yr_ag_nat_st"]][[st]][,,"usb"]
  notif_fb      <- cbind(notif_age_fb0[,12],notif_age_us0[,12])
  adj_3         <- sum(dDirMult(M=notif_fb+0.01,n=notif_fb,Rho=rho)*wts[44:69])
  #scale does not matter for dirichlet llikelihood
  (sum(dDirMult(M=V,n=notif_fb,Rho=rho)*wts[44:69]) - adj_3)*2
}
###############################################################################################

notif_fb_5yr_lLik_st <- function(V,st,rho=0.0005) { # V = table of notifications by fb 1993-2016 (row=24 years, col=fb then us)
  ### Read in the nativity stratified data
  notif_age_fb0   <- CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==0,4:8]
  notif_age_us0   <- CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==1,4:8]
  notif_fb      <- cbind(t(notif_age_fb0),t(notif_age_us0))
  adj_3         <- sum(dDirMult(M=notif_fb,n=notif_fb,Rho=rho)*wts[c(49,54,59,64,70)])
  V2<-matrix(0,5,2)
  V2[1,]<-colSums(V[1:5,]);V2[2,]<-colSums(V[6:10,]);V2[3,]<-colSums(V[11:15,]); V2[4,]<-colSums(V[16:20,]); V2[5,]<-colSums(V[21:25,])
  #scale does not matter for dirichlet llikelihood
  (sum(dDirMult(M=V2,n=notif_fb,Rho=rho)*wts[c(49,54,59,64,70)]) - adj_3)*2
}

#'FB Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@param V vector of total notifications 1953-2014
#'@return likelihood
notif_fb_5yr_lik <- function(V,st=st) {
  notif_fb   <- as.numeric(unlist(CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==0,4:8]))
  adj_1         <- sum(dnorm(notif_fb,notif_fb,notif_fb*0.05/1.96,log=T)*wts[c(49,54,59,64,70)])
  V2<-rep(0,5)
  V2[1]<-sum(V[1:5]);V2[2]<-sum(V[6:10]);V2[3]<-sum(V[11:15]); V2[4]<-sum(V[16:20]); V2[5]<-sum(V[21:25])
  (sum(dnorm(notif_fb,V2*1e6,notif_fb*0.05/1.96,log=T)*wts[c(49,54,59,64,70)]) - adj_1)
}

#'US Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@param V vector of total notifications 1953-2014
#'@return likelihood

notif_us_5yr_lik <- function(V,st=st) {
  notif_us   <- as.numeric(unlist(CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==1,4:8]))
  adj_1         <- sum(dnorm(notif_us,notif_us,notif_us*0.05/1.96,log=T)*wts[c(49,54,59,64,70)])
  V2<-rep(0,5)
  V2[1]<-sum(V[1:5]);V2[2]<-sum(V[6:10]);V2[3]<-sum(V[11:15]); V2[4]<-sum(V[16:20]); V2[5]<-sum(V[21:25])
  (sum(dnorm(notif_us,V2*1e6,notif_us*0.05/1.96,log=T)*wts[c(49,54,59,64,70)]) - adj_1)
}


### ### ### CASES FB DISTRIBUTION SLOPES OVER PAST 5 year  ### ### ### ### ### ###
notif_fbus_slp_lLik_st <- function(V,st) {
  # notif_fbus_slp5     <- as.numeric(CalibDatState[["case_change_5"]][st,2:3])
  notif_age_fb0   <- CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==0,7:8]
  notif_age_us0   <- CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==1,7:8]
  tot_case_nat<-cbind(t(notif_age_fb0),t(notif_age_us0))
  #calculate the slopes
  notif_fbus_slp5<-apply(log(tot_case_nat),2,function(x) lm(x~I(1:2))$coef[2])
  adj_3a              <- sum(dnorm(notif_fbus_slp5,notif_fbus_slp5,0.005,log=T))# V = table of notifications by fb 2011-2016 (row=6 years, col=fb then us)
  V1a <- matrix(c(sum(V[1:5,1]),sum(V[1:5,2]), sum(V[6:10,1]), sum(V[6:10,2])),2,2, byrow = TRUE)
  V2 <- apply(log(V1a),2,function(x) lm(x~I(1:2))$coef[2])
  sum(dnorm(notif_fbus_slp5,V2,0.005,log=T)) - adj_3a
}

### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

# notif_hr_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by tx history (row=97:16, col=n then e)
#   notif_hr0     <- CalibDatState[["hr_cases"]][[st]]
#   notif_hr      <- cbind(notif_hr0[,1],1-notif_hr0[,1])#*notif_us_hr0[,2]
#   adj_5b           <- sum(dDirMult(M=notif_hr+0.01,n=notif_hr,Rho=rho)*wts[c(45,50,55,60,65)])
#   V2 <- rbind(colSums(V[1:5,]),colSums(V[6:10,]),colSums(V[11:15,]),colSums(V[16:20,]), colSums(V[21:25,]))
#   #scale does not matter for dirichlet llikelihood
#   sum(dDirMult(M=V2,n=notif_hr,Rho=rho)*wts[c(45,50,55,60,65)]) - adj_5b
# }

###############################################################################################
### ### ### US CASES AGE DISTRIBUTION 5yrs 1993-2016  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
notif_age_us_5yr_lLik_st <- function(V,st,rho=0.1) { # V = table of us notifications by age 1993-2016 (row=24 years, col=11 ages)
  ### Read in the age and nativity stratified data
  notif_age_us_5yr     <- as.data.frame(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][CalibDat$cases_yr_ag_nat_st_5yr[[st]][,4]==1,5:14])
  ### Format the model estimates to match the calibration data
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  V3<-matrix(0,5,10)
  V3[1,]<-colSums(V2[1:5,]);V3[2,]<-colSums(V2[6:10,]);V3[3,]<-colSums(V2[11:15,]); V3[4,]<-colSums(V2[16:20,]); V3[5,]<-colSums(V2[21:25,])
  ### Check for missing calibration data
  if (sum(is.na(notif_age_us_5yr > 0))) {
    ### Initialize the likelihood value to zero
    tot_lik <- 0; adj_2a <-0
    ### If data is missing we need to create a single new bucket for these data and estimates
    ### Read in the total US cases to use as a total
    notif_age_us0   <- CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==1,4:8]
    for (i in 1:nrow(notif_age_us_5yr)){
      print(i)
      ### Check that this row has an NA
      if (sum(is.na(notif_age_us_5yr[i,])) > 0){
        ### Calibration data
        ### Which of the age groups are NA
        index <- which(is.na(notif_age_us_5yr[i,]))
        ### Remove these age groups from the dataframe
        known_cases <- notif_age_us_5yr[i,-index]
        ### Set the total of missing cases equal to the difference of total cases and known cases
        missing_cases <- notif_age_us0[i] - sum(known_cases)
        ### Add missing cases to the dataframe
        tot_cases <- c(unlist(known_cases), unlist(missing_cases))
        ### Model estimates (match to estimates)
        ### Create a sum of the estimates that correspond to missing age groups
        known_est   <- V3[i,-index]
        missing_est <- sum(V3) - sum(known_est)
        tot_est <- c(known_est, missing_est)
      } else {
        tot_est <- V3[i,]
        tot_cases <- notif_age_us_5yr[i,]
      }
      adj_2a  <- adj_2a +  dDirMult(M=tot_cases+0.01,n=tot_cases,Rho=rho)*wts[c(49,54,59,64,70)][i]
      print(adj_2a)
      tot_lik <- tot_lik + dDirMult(M=tot_est*1e6,n=tot_cases,Rho=rho)*wts[c(49,54,59,64,70)][i]
      }
  } else {
    adj_2a            <- sum(dDirMult(M=notif_age_us_5yr+0.01,n=notif_age_us_5yr,Rho=rho)*wts[c(49,54,59,64,70)])
    #scale does not matter for dirichlet llikelihood
    tot_lik <- sum(dDirMult(M=V3*1e6,n=notif_age_us_5yr,Rho=rho)*wts[c(49,54,59,64,70)]) - adj_2a
  }
  return(tot_lik)
}

### ### ### NUS CASES AGE DISTRIBUTION 5yrs 1993-2016  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
notif_age_nus_5yr_lLik_st <- function(V,st,rho=0.1) { # V = table of us notifications by age 1993-2016 (row=24 years, col=11 ages)
  ### Read in the age and nativity stratified data
  notif_age_nus_5yr     <- as.data.frame(CalibDatState$cases_yr_ag_nat_st_5yr[[st]][CalibDat$cases_yr_ag_nat_st_5yr[[st]][,4]==0,5:14])
  ### Format the model estimates to match the calibration data
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  V3<-matrix(0,5,10)
  V3[1,]<-colSums(V2[1:5,]);V3[2,]<-colSums(V2[6:10,]);V3[3,]<-colSums(V2[11:15,]); V3[4,]<-colSums(V2[16:20,]); V3[5,]<-colSums(V2[21:25,])
  ### Check for missing calibration data
  if (sum(is.na(notif_age_nus_5yr > 0))) {
    ### Initialize the likelihood value to zero
    tot_lik <- 0; adj_2a <-0
    ### If data is missing we need to create a single new bucket for these data and estimates
    ### Read in the total US cases to use as a total
    notif_age_nus0   <- CalibDatState$cases_nat_st_5yr[CalibDatState$cases_nat_st_5yr$State.Code==st & CalibDatState$cases_nat_st_5yr$usb==0,4:8]
    for (i in 1:nrow(notif_age_nus_5yr)){
      print(i)
      ### Check that this row has an NA
      if (sum(is.na(notif_age_nus_5yr[i,])) > 0){
        ### Calibration data
        ### Which of the age groups are NA
        index <- which(is.na(notif_age_nus_5yr[i,]))
        ### Remove these age groups from the dataframe
        known_cases <- notif_age_nus_5yr[i,-index]
        ### Set the total of missing cases equal to the difference of total cases and known cases
        missing_cases <- notif_age_nus0[i] - sum(known_cases)
        ### Add missing cases to the dataframe
        tot_cases <- c(unlist(known_cases), unlist(missing_cases))
        ### Model estimates (match to estimates)
        ### Create a sum of the estimates that correspond to missing age groups
        known_est   <- V3[i,-index]
        missing_est <- sum(V3) - sum(known_est)
        tot_est <- c(known_est, missing_est)
      } else {
        tot_est <- V3[i,]
        tot_cases <- notif_age_nus_5yr[i,]
      }
      adj_2a  <- adj_2a +  dDirMult(M=tot_cases+0.01,n=tot_cases,Rho=rho)*wts[c(49,54,59,64,70)][i]
      tot_lik <- tot_lik + dDirMult(M=tot_est*1e6,n=tot_cases,Rho=rho)*wts[c(49,54,59,64,70)][i]
    }
  } else {
    adj_2a            <- sum(dDirMult(M=notif_age_nus_5yr+0.01,n=notif_age_nus_5yr,Rho=rho)*wts[c(49,54,59,64,70)])
    #scale does not matter for dirichlet llikelihood
    tot_lik <- sum(dDirMult(M=V3*1e6,n=notif_age_nus_5yr,Rho=rho)*wts[c(49,54,59,64,70)]) - adj_2a
  }
  return(tot_lik)
}

##smoothed estimates
notif_hr_lLik_st <- function(V,st) { # V = table of notifications by tx history (row=97:16, col=n then e)
  notif_hr_5yr<-CalibDatState[["hr_cases_sm"]][which(CalibDatState[["hr_cases_sm"]][,1]==stateID[st,1]),7]*
                CalibDatState[["hr_cases_sm"]][which(CalibDatState[["hr_cases_sm"]][,1]==stateID[st,1]),3]
  # notif_5yr          <-    c(sum(CalibDatState[["cases_yr_st"]][[st]][2:6,2]),
  #                             sum(CalibDatState[["cases_yr_st"]][[st]][7:11,2]),
  #                             sum(CalibDatState[["cases_yr_st"]][[st]][12:16,2]),
  #                             sum(CalibDatState[["cases_yr_st"]][[st]][17:21,2]),
  #                             sum(CalibDatState[["cases_yr_st"]][[st]][22:26,2]))
  # notif_hr_5yr      <- notif_hr*notif_5yr
  adj_6           <- sum(dnorm(notif_hr_5yr,notif_hr_5yr,notif_hr_5yr*0.1/1.96,log=T)*wts[c(45,50,55,60,65)])
  sum(dnorm(notif_hr_5yr,V*1e6,notif_hr_5yr*0.1/1.96,log=T)*wts[c(45,50,55,60,65)]) - adj_6
}

#' DISTRIBUTION OF CASES RECENT TRANSMISSION VS NO RECENT TRANSMISSION
#'@param V distribution of recent and not recent transmission
#'@return likelihood
recent_trans_dist_lLik_st  <- function(V,st) {
  rct_trans_dist        <- CalibDat[["rct_cases_sm"]][st,5]
  adj_13          <- dbeta(rct_trans_dist,rct_trans_dist*100,(1-rct_trans_dist)*100,log=T)
  dbeta(rct_trans_dist,V[1]*100,V[2]*100,log=T)  - adj_13
}

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
# Motivation: should be a normal distribution because it is based on a model result
notif_fb_rec_lLik_st<-function(V,st){
  notif_rec<-CalibDatState[["rt_fb_cases_sm"]][which(CalibDatState[["rt_fb_cases_sm"]][,1]==stateID[st,1]),9]
  notif_fb          <-  rbind(sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][2:6,12,"nusb"]),
                              sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][7:11,12,"nusb"]),
                              sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][12:16,12,"nusb"]),
                              sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][17:21,12,"nusb"]),
                              sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][22:26,12,"nusb"]))
  notif_fb_rec      <- notif_rec*notif_fb
  adj_6           <- sum(dnorm(notif_fb_rec,notif_fb_rec,notif_fb_rec*0.1/1.96,log=T)*wts[c(45,50,55,60,65)])
  sum(dnorm(notif_fb_rec,V*1e6,notif_fb_rec*0.1/1.96,log=T)*wts[c(45,50,55,60,65)]) - adj_6
}

# notif_fb_rec_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by rec 1993-2014 (row=22 years, col=pos then neg)
#   notif_rec<-CalibDatState[["rt_fb_cases_sm"]][which(CalibDatState[["rt_fb_cases_sm"]][,1]==stateID[st,1]),9]
#   notif_fb          <-  rbind(sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][2:6,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][7:11,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][12:16,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][17:21,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][22:26,12,"nusb"]))
#   notif_fb_rec      <- cbind(notif_rec*notif_fb, (1-notif_rec)*notif_fb)
#   adj_6             <- sum(dDirMult(M=notif_fb_rec,n=notif_fb_rec,Rho=rho)*wts[c(45,50,55,60,65)])
#   #scale does not matter for dirichlet llikelihood
#   (sum(dDirMult(M=V,n=notif_fb_rec,Rho=rho)*wts[c(45,50,55,60,65)]) - adj_6)*5
#   }

### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

tx_outcomes_lLik_st <- function(V,rho=0.01) {
  tx_outcomes      <- (cbind(1-rowSums(CalibDatState[["tx_outcomes"]][,2:3]),CalibDatState[["tx_outcomes"]][,2],CalibDatState[["tx_outcomes"]][,3])*CalibDatState[["tx_outcomes"]][,4])
  adj_11           <- sum(dDirMult(M=tx_outcomes+0.01,n=tx_outcomes,Rho=0.01)*wts[44:66])# V = table of treatment outcomes 1993-2012 (row=20 years, col= complete, discontinue, dead)
  #scale does not matter for dirichlet llikelihood
  sum(dDirMult(M=V,n=tx_outcomes,Rho=rho)*wts[44:66]) - adj_11
  }

### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D

tltbi_tot_lLik_st   <- function(V,st) { # V = total TLTBI inits in 2002 (scalar)
  # Motivation: norm, mean centered with CI = +/- 10% of mean
  tltbi_vol        <- CalibDatState[["TLTBI_volume_state"]][[st]]
  adj_12           <- dnorm(tltbi_vol[1],tltbi_vol[1],diff(tltbi_vol[2:3])/1.96,log=T)
  dnorm(tltbi_vol[1],V*1e6,diff(tltbi_vol[2:3])/1.96,log=T) - adj_12
  }

### ### ### DISTRIBUTION OF LTBI TREATMENT INITS 2002  ### ### ### ### ### ###  D
tltbi_dist_lLik_st  <- function(V) {
  TLTBI_dist       <- CalibDatState[["TLTBI_dist"]][1:2]
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
### ### ### Total TB DEATHS 1999-2016 ### ### ### ### ### ### D

tbdeaths_lLik_st <- function(V,st) { # V = vector of total notifications 1999-2016
  tb_deaths <- as.numeric(CalibDatState[["tbdeaths"]][[st]][,3])
  V2<-rowSums(V)*1e6
  adj_19    <- sum((dnorm(tb_deaths,tb_deaths,tb_deaths*0.1/1.96,log=T)*wts[50:69])[is.na(tb_deaths)==F])
  sum((dnorm(tb_deaths,V2,tb_deaths*0.1/1.96,log=T)*wts[50:69])[is.na(tb_deaths)==F]) - adj_19
}
### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D

tbdeaths_decline_lLik_st <- function(V) { # V = vector of tb deaths 1968-2015
  tbdeaths_decline      <- CalibDatState[["deaths_ann_decline_68_15"]]
  adj_19a            <- sum(dnorm(tbdeaths_decline,tbdeaths_decline,0.015/1.96,log=T))
  V2 <- (1-(V[48]/V[1])^(1/47))
  sum(dnorm(tbdeaths_decline,V2,0.015/1.96,log=T)) - adj_19a
}
### ### ### TB DEATHS AGE DISTRIBUTION 1999-2016  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

tb_dth_age_lLik_st <- function(V,rho=0.005) { # V = table of deaths by age 1999-2016 (row=18 years, col=11 ages)
  tb_deaths_age  <- CalibDatState[["tbdeaths_age_yr"]][,-1]
  adj_19b        <- sum(dDirMult(M=tb_deaths_age+0.005,n=tb_deaths_age,Rho=rho)*wts[50:69])
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  sum(dDirMult(M=V2,n=tb_deaths_age,Rho=rho)*wts[50:69]) - adj_19b
}
### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
# Motivation: norm, mean centered with CI = +/- 2 million wts[1+0:6*10]

tot_pop_yr_fb_lLik_st <- function(V,st) { # V = total pop (rows=year, cols=us, fb)
  tot_pop_yr      <- CalibDatState[["pop_50_10"]][[st]]
  tot_pop_yr_fb   <- tot_pop_yr[tot_pop_yr[,2]==0,]
  #get 2019 population
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c(1,2,22)]
  #get 2019 fb population
  pop_ag_11_19nus <-sum(pop_ag_11_190[pop_ag_11_190[,2]==0,3][-11])
  #append the foreign born population
  tot_pop_yr_fb   <- c(tot_pop_yr_fb[,-c(1:2)], pop_ag_11_19nus)
  # if (loc != "HI" & loc != "AK"){
  adj_17          <- sum(dnorm(tot_pop_yr_fb[-1],tot_pop_yr_fb[-1],tot_pop_yr_fb[7]*0.05/1.96,log=T)*wts[c(1+1:6*10,70)])
  #total population is in real numbers so we need to scale up output
  sum(dnorm(tot_pop_yr_fb[-1],V[c(11,21,31,41,51,61,70)]*1e6,tot_pop_yr_fb[7]*0.05/1.96,log=T)*wts[c(1+1:6*10,70)]) - adj_17}
  # else{
  #   adj_17          <- sum(dnorm(tot_pop_yr_fb,tot_pop_yr_fb,tot_pop_yr_fb[7]*0.1/1.96,log=T)*wts[c(1+1:6*10,68)])
  #   #total population is in real numbers so we need to scale up output
  #   sum(dnorm(tot_pop_yr_fb[-1],V[c(11,21,31,41,51,61,68)]*1e6,tot_pop_yr_fb[7]*0.1/1.96,log=T)*wts[c(1+1:6*10,68)]) - adj_17
  # }
# }
### D
tot_pop_yr_us_lLik_st_00_10 <- function(V,st) {
  tot_pop_yr      <- CalibDatState[["pop_50_10"]][[st]]
  # V = total pop (rows=year, cols=us, fb)
  tot_pop_yr_us  <- tot_pop_yr[tot_pop_yr[,2]==1,3]
  # tot_pop_yr_us<-as.matrix(tot_pop_yr_us)
  # tot_pop_yr_us<-colSums(as.matrix(tot_pop_yr_us)[,-c(1:2)])
  adj_17b        <- sum(dnorm(tot_pop_yr_us[6:7],tot_pop_yr_us[6:7],tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[1+5:6*10])
  sum(dnorm(tot_pop_yr_us[6:7],V[c(51,61)]*1e6,tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[1+5:6*10]) - adj_17b  } # CI = +/- 2mil

tot_pop_yr_us_lLik_st <- function(V,st) {

  tot_pop_yr      <- CalibDatState[["pop_50_10"]][[st]]
  tot_pop_yr_us   <- tot_pop_yr[tot_pop_yr[,2]==1,]
  #get 2019 population
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c(1,2,22)]
  #get 2019 fb population
  pop_ag_11_19us <-sum(pop_ag_11_190[pop_ag_11_190[,2]==1,3][-11])
  #append the foreign born population
  tot_pop_yr_us   <- c(tot_pop_yr_us[,-c(1:2)], pop_ag_11_19us)
  adj_17          <- sum(dnorm(tot_pop_yr_us[-1],tot_pop_yr_us[-1],tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[c(1+1:6*10,70)])
  #total population is in real numbers so we need to scale up output
  sum(dnorm(tot_pop_yr_us[-1],V[c(11,21,31,41,51,61,70)]*1e6,tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[c(1+1:6*10,70)]) - adj_17}

### ### ### TOTAL POP AGE DISTRIBUTION 2017  ### ### ### ### ### ### D
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
tot_pop19_ag_fb_lLik_st <- function(V,st,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c(1,2,22)]
  pop_ag_11_19us <-pop_ag_11_190[pop_ag_11_190[,2]==1,3][-11]
  pop_ag_11_19nus <-pop_ag_11_190[pop_ag_11_190[,2]==0,3][-11]

  pop_ag_11_19   <- cbind(pop_ag_11_19us/sum(pop_ag_11_19us)+.001, pop_ag_11_19nus/sum(pop_ag_11_19nus)+.001)
  adj_18         <- (sum(log(pop_ag_11_19[,1])*pop_ag_11_19[,1])+sum(log(pop_ag_11_19[,2])*pop_ag_11_19[,2]))*ESS
  V1 <- rbind(V[1:9,],V[10,]+V[11,])
  V2<-cbind(V1[,1]/sum(V1[,1]),V1[,2]/sum(V1[,2]))
  (sum(log(V2[,1])*pop_ag_11_19[,1])+sum(log(V2[,2])*pop_ag_11_19[,2]))*ESS - adj_18
}

### ### ### TOTAL POP AGE DISTRIBUTION 2017  ### ### ### ### ### ### D
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
tot_pop1719_ag_fb_lLik_st <- function(V,st,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c(1,2,20:22)]
  pop_ag_11_19us <- rowSums(pop_ag_11_190[pop_ag_11_190[,2]==1,3:5])[-11]
  #top code at 75p
  # pop_ag_11_19us[9] <- pop_ag_11_19us[9] + pop_ag_11_19us[10]
  # pop_ag_11_19us <- pop_ag_11_19us[-10]
  pop_ag_11_19nus <-rowSums(pop_ag_11_190[pop_ag_11_190[,2]==0,3:5])[-11]
  #top code at 75p
  pop_ag_11_19   <- cbind(pop_ag_11_19us/sum(pop_ag_11_19us), pop_ag_11_19nus/sum(pop_ag_11_19nus))
  adj_18         <- (sum(log(pop_ag_11_19[,1])*pop_ag_11_19[,1])+sum(log(pop_ag_11_19[,2])*pop_ag_11_19[,2]))*ESS
  V1 <- rbind(V[1:9,],V[10,]+V[11,])
  V2<-cbind(V1[,1]/sum(V1[,1]),V1[,2]/sum(V1[,2]))
  (sum(log(V2[,1])*pop_ag_11_19[,1])+sum(log(V2[,2])*pop_ag_11_19[,2]))*ESS - adj_18
}

#' TOTAL US DEATHS
#' 1970,1975,1980,1985,1990-2007
#' Motivation: norm, mean centered with CI = +/- 5% of mean
#'@name dth_tot_lLik_st
#'@param V
#'@return likelihood
dth_tot_lLik_st <- function(V,st) {
  ST_deaths_tot <- readRDS(system.file("ST/STdeathbyAge.rds",package="MITUS"))[[st]][48,12]
  adj_20a         <- sum(dnorm(ST_deaths_tot,ST_deaths_tot,ST_deaths_tot*0.1/1.96,log=T)*wts[67])
  sum(dnorm(ST_deaths_tot,V*1e6,ST_deaths_tot*0.1/1.96,log=T)*wts[67]) - adj_20a
}

#'  #' TOTAL DEATHS AGE DISTRIBUTION 1999-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of deaths by age 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
tot_dth_age_lLik_st <- function(V,st,rho=0.01) {
  tda <- readRDS(system.file("ST/STdeathbyAge.rds",package="MITUS"))[[st]][47:48,-c(1,12)]
  adj_20b        <- sum(dDirMult(M=tda+0.1,n=tda,Rho=rho)*wts[66:67])
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  # V2<-V2*1e6
  sum(dDirMult(M=V2,n=tda,Rho=rho)*wts[66:67]) - adj_20b
  }

#' Mortality Risk Group Distribution 1999-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@name mort_dist_lLik_st
#'@param V table of mort_dist 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
mort_dist_lLik_st <- function(V,rho=0.1) {
  md     <- rowSums(dist_gen)
  mort_dist     <-matrix(md,length(66:67),4, byrow = TRUE)
  adj_21        <- sum(dDirMult(M=mort_dist+0.01,n=mort_dist,Rho=0.1)*wts[66:67])
  tot_lik<-0
  for(ag in 1:11){
    V1<-V[,(1:4)+4*(ag-1)]
    x<-sum(dDirMult(M=(V1*1e6),n=mort_dist,Rho=rho)*wts[66:67]) - adj_21
    tot_lik<-tot_lik+x
    # print(x)
  }
  tot_lik<-tot_lik
  return(tot_lik)
}

#' Mortality Risk Group Distribution 1999-2014
#' Motivation: dnorm
#'@name mort_dist_lLik_norm_st
#'@param V table of mort_dist 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
mort_dist_lLik_norm_st <- function(V) {
  md     <- rowSums(dist_gen)
  mort_dist     <-matrix(md,length(51:67),4, byrow = TRUE)
  adj_21b        <- sum(dnorm(mort_dist,mort_dist,mort_dist*0.1/1.96, log=T)*wts[51:67])
  tot_lik<-0
  for(ag in 1:11){
    V1<-V[,(1:4)+4*(ag-1)]
    x<-sum(dnorm(mort_dist, V1,mort_dist*0.1/1.96, log=T)*wts[51:67]) - adj_21b
    tot_lik<-tot_lik+x
    # print(x)
  }
  tot_lik<-tot_lik
  return(tot_lik)
}
### ### ### HOMELESS POP 2010  ### ### ### ### ### ### names(CalibDatState)
# Motivation: norm, mean centered with CI = +/- 25% of mean

homeless_10_lLik_st <- function(V,st) { # V = homeless pop in 2010 (scalar)
  homeless_pop      <- CalibDatState[["homeless_pop"]][[st]]
  adj_23b          <- dnorm(homeless_pop[1],homeless_pop[1],diff(homeless_pop[2:3])/2/1.96,log=T)
  dnorm(homeless_pop[1],V,diff(homeless_pop[2:3])/2/1.96,log=T) - adj_23b   }

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
    # p1 <- as.numeric(t(p0)%*%v21a[17,41:44])

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
