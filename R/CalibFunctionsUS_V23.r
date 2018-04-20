###  Calibration functions
  library(MCMCpack)

########################  ########################  ########################
###########################  CALIBRATION DATA  ########################
########################  ########################  ########################
  load("data/CalibDat_9-14-16.rData") # CalibDat
  # names(CalibDat)

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

### ### ### calibration importance weights  ### ### ### ### ### ### ~~~*
  wts <- CalibDat[["ImptWeights"]]

### ### ### TOTAL DIAGNOSED CASES 1953-2014  ### ### ### ### ### ### ~~~*
# Motivation: norm, mean centered with CI = +/- 5% of mean
  notif_tot     <- CalibDat[["tot_cases"]][,2];
  adj_1         <- sum(dnorm(notif_tot,notif_tot,notif_tot*0.1/1.96,log=T)*wts[4:66])
  notif_tot_lik <- function(V) { # V = vector of total notifications 1953-2014
    sum(dnorm(notif_tot,V*1e6,notif_tot*0.1/1.96,log=T)*wts[4:66]) - adj_1  }

### ### ### US CASES AGE DISTRIBUTION 1993-2013  ### ### ### ### ### ### names(CalibDat)
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_age     <- CalibDat[["age_cases"]][,-c(1,12)]*CalibDat[["age_cases"]][,12]
  notif_age_us     <- CalibDat[["age_cases_us"]][,-c(1,12)]*CalibDat[["age_cases_us"]][,12]
  adj_2a           <- sum(dDirMult(M=notif_age_us,n=notif_age_us,Rho=0.005)*wts[44:65])
  notif_age_us_lLik <- function(V,rho=0.005) { # V = table of us notifications by age 1993-2013 (row=21 years, col=11 ages)
    V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=notif_age_us,Rho=rho)*wts[44:65]) - adj_2a  }

  ### ### ### FB CASES AGE DISTRIBUTION 1993-2013  ### ### ### ### ### ###
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_age_fb     <- CalibDat[["age_cases_fb"]][,-c(1,12)]*CalibDat[["age_cases_fb"]][,12]
  adj_2b           <- sum(dDirMult(M=notif_age_fb,n=notif_age_fb,Rho=0.005)*wts[44:65])
  notif_age_fb_lLik <- function(V,rho=0.005) { # V = table of fb notifications by age 1993-2013 (row=21 years, col=11 ages)
    V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=notif_age_fb,Rho=rho)*wts[44:65]) - adj_2b  }

### ### ### CASES FB DISTRIBUTION 1993-2014  ### ### ### ### ### ###  ~~~*
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_fb      <- cbind(CalibDat[["fb_cases"]][,2],1-CalibDat[["fb_cases"]][,2])*CalibDat[["fb_cases"]][,3]
  adj_3         <- sum(dDirMult(M=notif_fb,n=notif_fb,Rho=0.005)*wts[44:66])
  notif_fb_lLik <- function(V,rho=0.005) { # V = table of notifications by fb 1993-2014 (row=22 years, col=fb then us)
    sum(dDirMult(M=V,n=notif_fb,Rho=rho)*wts[44:66]) - adj_3  }

### ### ### CASES FB DISTRIBUTION SLOPES OVER PAST 4 year  ### ### ### ### ### ###  ~~~*
  notif_fbus_slp5     <- CalibDat[["fbus_cases_slope5"]];
  adj_3a         <- sum(dnorm(notif_fbus_slp5,notif_fbus_slp5,0.01/2,log=T))
  notif_fbus_slp_lLik <- function(V) { # V = table of notifications by fb 1993-2014 (row=22 years, col=fb then us)
    V2 <- apply(log(V[20:23,]),2,function(x) lm(x~I(1:4))$coef[2])
    sum(dnorm(V2,notif_fbus_slp5,0.01/2,log=T)) - adj_3a  }

### ### ### CASES TX HISTORY DISTRIBUTION 1993-2013  ### ### ### ### ### ###
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_prev      <- cbind(1-CalibDat[["prev_cases"]][,2],CalibDat[["prev_cases"]][,2])*CalibDat[["prev_cases"]][,3]
  adj_4           <- sum(dDirMult(M=notif_prev,n=notif_prev,Rho=0.005)*wts[44:65])
  notif_prev_lLik <- function(V,rho=0.005) { # V = table of notifications by tx history 1993-2013 (row=21 years, col=n then e)
    sum(dDirMult(M=V,n=notif_prev,Rho=rho)*wts[44:65]) - adj_4  }

  ### ### ### CASES HIV DISTRIBUTION 1993-2014  ### ### ### ### ### ###  ~~~*
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_hiv0      <- CalibDat[["hiv_cases2"]]
  ppp <- notif_hiv0[,2]/rowSums(notif_hiv0[,2:3])/2
  notif_hiv1 <- notif_hiv0[,2] + notif_hiv0[,4]*ppp
  notif_hiv <- cbind(notif_hiv1,notif_hiv0[,4]-notif_hiv1+notif_hiv0[,3])
  adj_5          <- sum(dDirMult(M=notif_hiv,n=notif_hiv,Rho=0.005)*wts[44:65])
  notif_hiv_lLik <- function(V,rho=0.005) { # V = table of notifications by hiv 1993-2014 (row=22 years, col=pos then neg)
    sum(dDirMult(M=V,n=notif_hiv,Rho=rho)*wts[44:65]) - adj_5  }

  ### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ###  ~~~*
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_us_hr      <- cbind(CalibDat[["us_homeless_cases"]][,2],1-CalibDat[["us_homeless_cases"]][,2])*CalibDat[["us_homeless_cases"]][,3]
  adj_5b          <- sum(dDirMult(M=notif_us_hr,n=notif_us_hr,Rho=0.005)*wts[44:65])
  notif_us_hr_lLik <- function(V,rho=0.005) { # V = table of notifications by hr 1993-2014 (row=22 years, col=pos then neg)
    sum(dDirMult(M=V,n=notif_us_hr,Rho=rho)*wts[44:65]) - adj_5b  }

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ###
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_fb_rec      <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]
  adj_6          <- sum(dDirMult(M=notif_fb_rec,n=notif_fb_rec,Rho=0.005)*wts[44:65])
  notif_fb_rec_lLik <- function(V,rho=0.005) { # V = table of notifications by hiv 1993-2014 (row=22 years, col=pos then neg)
    sum(dDirMult(M=V,n=notif_fb_rec,Rho=rho)*wts[44:65]) - adj_6  }

### ### ### CASES PCT INH & MDR RES US N 1993-2014  ### ### ### ### ### ###  ~~~*
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_inh_us_n     <- cbind(CalibDat[["inh_res_us_n"]][,2],1-CalibDat[["inh_res_us_n"]][,2])*CalibDat[["inh_res_us_n"]][,3]
  notif_mdr_us_n     <- cbind(CalibDat[["mdr_res_us_n"]][,2],1-CalibDat[["mdr_res_us_n"]][,2])*CalibDat[["mdr_res_us_n"]][,3]
  adj_7              <- sum(dDirMult(M=notif_inh_us_n,n=notif_inh_us_n,Rho=0.005)*wts[44:65])+sum(dDirMult(M=notif_mdr_us_n,n=notif_mdr_us_n,Rho=0.005)*wts[44:65])
  notif_dr_us_n_lLik <- function(V,rho=0.005) { # V = table of US N notifications by dr 1993-2014 (row=22 years, col= 5 dr cats )
    Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
    Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    sum(dDirMult(M=Vinh,n=notif_inh_us_n,Rho=rho)*wts[44:65]) + sum(dDirMult(M=Vmdr,n=notif_mdr_us_n,Rho=rho)*wts[44:65]) - adj_7  }

### ### ### CASES PCT INH & MDR RES US E 1993-2014  ### ### ### ### ### ###  ~~~*
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_inh_us_e     <- cbind(CalibDat[["inh_res_us_e"]][,2],1-CalibDat[["inh_res_us_e"]][,2])*CalibDat[["inh_res_us_e"]][,3]
  notif_mdr_us_e     <- cbind(CalibDat[["mdr_res_us_e"]][,2],1-CalibDat[["mdr_res_us_e"]][,2])*CalibDat[["mdr_res_us_e"]][,3]
  adj_8              <- sum(dDirMult(M=notif_inh_us_e,n=notif_inh_us_e,Rho=0.005)*wts[44:65])+sum(dDirMult(M=notif_mdr_us_e+.1,n=notif_mdr_us_e,Rho=0.005)*wts[44:65])
  notif_dr_us_e_lLik <- function(V,rho=0.005) { # V = table of US N notifications by dr 1993-2014 (row=22 years, col= 5 dr cats )
    Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
    Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    sum(dDirMult(M=Vinh,n=notif_inh_us_e,Rho=rho)*wts[44:65]) + sum(dDirMult(M=Vmdr,n=notif_mdr_us_e,Rho=rho)*wts[44:65]) - adj_8  }

### ### ### CASES PCT INH & MDR RES FB N 1993-2014  ### ### ### ### ### ###  ~~~*
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_inh_fb_n     <- cbind(CalibDat[["inh_res_fb_n"]][,2],1-CalibDat[["inh_res_fb_n"]][,2])*CalibDat[["inh_res_fb_n"]][,3]
  notif_mdr_fb_n     <- cbind(CalibDat[["mdr_res_fb_n"]][,2],1-CalibDat[["mdr_res_fb_n"]][,2])*CalibDat[["mdr_res_fb_n"]][,3]
  adj_9              <- sum(dDirMult(M=notif_inh_fb_n,n=notif_inh_fb_n,Rho=0.005)*wts[44:65])+sum(dDirMult(M=notif_mdr_fb_n,n=notif_mdr_fb_n,Rho=0.005)*wts[44:65])
  notif_dr_fb_n_lLik <- function(V,rho=0.005) { # V = table of FB N notifications by dr 1993-2014 (row=22 years, col= 5 dr cats )
    Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
    Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    sum(dDirMult(M=Vinh,n=notif_inh_fb_n,Rho=rho)*wts[44:65]) + sum(dDirMult(M=Vmdr,n=notif_mdr_fb_n,Rho=rho)*wts[44:65]) - adj_9  }

### ### ### CASES PCT INH & MDR RES FB E 1993-2014  ### ### ### ### ### ###  ~~~*
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  notif_inh_fb_e     <- cbind(CalibDat[["inh_res_fb_e"]][,2],1-CalibDat[["inh_res_fb_e"]][,2])*CalibDat[["inh_res_fb_e"]][,3]
  notif_mdr_fb_e     <- cbind(CalibDat[["mdr_res_fb_e"]][,2],1-CalibDat[["mdr_res_fb_e"]][,2])*CalibDat[["mdr_res_fb_e"]][,3]
  adj_10             <- sum(dDirMult(M=notif_inh_fb_e,n=notif_inh_fb_e,Rho=0.005)*wts[44:65])+sum(dDirMult(M=notif_mdr_fb_e,n=notif_mdr_fb_e,Rho=0.005)*wts[44:65])
  notif_dr_fb_e_lLik <- function(V,rho=0.005) { # V = table of FB N notifications by dr 1993-2014 (row=22 years, col= 5 dr cats )
    Vinh <- cbind(rowSums(V[,c(2,4,5)]),rowSums(V[,c(1,3  )]))/rowSums(V)
    Vmdr <- cbind(rowSums(V[,c(4,5  )]),rowSums(V[,c(1,2,3)]))/rowSums(V)
    sum(dDirMult(M=Vinh,n=notif_inh_fb_e,Rho=rho)*wts[44:65]) + sum(dDirMult(M=Vmdr,n=notif_mdr_fb_e,Rho=rho)*wts[44:65]) - adj_10  }

### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ###  ~~~*
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  tx_outcomes      <- cbind(1-rowSums(CalibDat[["tx_outcomes"]][,2:3]),CalibDat[["tx_outcomes"]][,2],CalibDat[["tx_outcomes"]][,3])*CalibDat[["tx_outcomes"]][,4]
  adj_11           <- sum(dDirMult(M=tx_outcomes,n=tx_outcomes,Rho=0.005)*wts[44:63])
  tx_outcomes_lLik <- function(V,rho=0.005) { # V = table of treatment outcomes 1993-2012 (row=20 years, col= complete, discontinue, dead)
    sum(dDirMult(M=V,n=tx_outcomes,Rho=rho)*wts[44:63]) - adj_11  }

### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ###
# Motivation: norm, mean centered with CI = +/- 10% of mean
  tltbi_vol        <- CalibDat[["TLTBI_volume"]]
  adj_12           <- dnorm(tltbi_vol[1],tltbi_vol[1],diff(tltbi_vol[2:3])/2/1.96,log=T)
  tltbi_tot_lLik   <- function(V) { # V = total TLTBI inits in 2002 (scalar)
    dnorm(tltbi_vol[1],V*1e6,diff(tltbi_vol[2:3])/2/1.96,log=T) - adj_12  }

### ### ### DISTRIBUTION OF LTBI TREATMENT INITS 2002  ### ### ### ### ### ###
  TLTBI_dist       <- CalibDat[["TLTBI_dist"]]
  adj_13          <- sum( dbeta(TLTBI_dist,TLTBI_dist*100,(1-TLTBI_dist)*100,log=T) )
  tltbi_dist_lLik  <- function(V) { # V = dist TLTBI inits in 2002 (vector fraction FB, HR, HV in 2002)
    sum( dbeta(TLTBI_dist,V*100,(1-V)*100,log=T) ) - adj_13  }

### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ###
# Motivation: additional prior on LTBI, using beta densities parameterized to Miramontes/Hill results
  ltbi_us_11      <- CalibDat[["LTBI_prev_US_11_IGRA"]]
  adj_15          <- sum( dbeta(ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) )
  ltbi_us_11_lLik <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) ) - adj_15)*2  }

  ltbi_us_11_dp      <- CalibDat[["LTBI_prev_US_11_DoubPos"]]
  adj_15dp           <- sum( dbeta(ltbi_us_11_dp[,2]/rowSums(ltbi_us_11_dp[,2:3]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) )
  ltbi_us_11_dp_lLik <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) ) - adj_15dp)*2  }

### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ###
# Motivation: multinomial adjusted to match effective sample size due to survey weighting
  ltbi_fb_11      <- CalibDat[["LTBI_prev_FB_11_IGRA"]]
  adj_16          <- sum( dbeta(ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) )
  ltbi_fb_11_lLik <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) ) - adj_16)*2  }

  ltbi_fb_11_dp      <- CalibDat[["LTBI_prev_FB_11_DoubPos"]]
  adj_16dp           <- sum( dbeta(ltbi_fb_11_dp[,2]/rowSums(ltbi_fb_11_dp[,2:3]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) )
  ltbi_fb_11_dp_lLik <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
    V[9,] <- colSums(V[9:11,])
    (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) ) - adj_16dp)*2  }

### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###
# Motivation: norm, mean centered with CI = +/- 2 million wts[1+0:6*10]
  tot_pop_yr_fb      <- CalibDat[["tot_pop_yr_fb"]]
  adj_17             <- sum(dnorm(tot_pop_yr_fb[-1,4],tot_pop_yr_fb[-1,4],tot_pop_yr_fb[7,4]*0.1/1.96,log=T)*wts[1+1:6*10])
  tot_pop_yr_fb_lLik <- function(V) { # V = total pop (rows=year, cols=us, fb)
    sum(dnorm(tot_pop_yr_fb[-1,4],V[c(11,21,31,41,51,61)],tot_pop_yr_fb[7,4]*0.1/1.96,log=T)*wts[1+1:6*10]) - adj_17  } # CI = +/- 2mil

### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
  tot_pop14_ag_fb      <- cbind(CalibDat[["tot_pop14_ag_fb"]][-9,3]/sum(CalibDat[["tot_pop14_ag_fb"]][-9,3]),
                           CalibDat[["tot_pop14_ag_fb"]][-9,4]/sum(CalibDat[["tot_pop14_ag_fb"]][-9,4]))
  adj_18               <- sum(log(tot_pop14_ag_fb[,1])*tot_pop14_ag_fb[,1])+sum(log(tot_pop14_ag_fb[,2])*tot_pop14_ag_fb[,2])
  tot_pop14_ag_fb_lLik <- function(V,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
    V1 <- rbind(V[1,],V[2,]+V[3,],V[4,]+V[5,],V[6,],V[7,],V[8,],V[9,],V[10,]+V[11,])
    (sum(log(V1[,1]/sum(V1[,1]))*tot_pop14_ag_fb[,1])+sum(log(V1[,2]/sum(V1[,2]))*tot_pop14_ag_fb[,2]))*ESS - adj_18*ESS  }

### ### ### Total TB DEATHS 1999-2014 ### ### ### ### ### ###
# Motivation: overdispersed poisson, modelled with negbin with overdispersion param = 100 *wts[50:65]
  tb_deaths      <- CalibDat[["tb_deaths"]][,-1]
  adj_19         <- 0
  for(i in 1:16) adj_19 <- adj_19 + sum(dnbinom(as.numeric(tb_deaths[i,]),mu=as.numeric(tb_deaths[i,]),size=50,log=T))*wts[50:65][i]
  tb_deaths_lLik <- function(V,sgsq=50) { # V =  TB deaths by age 1999-2013  (row=15 years, col= 11 ages)
    V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
    l1 <- 0
    for(i in 1:16) l1 <- l1 + sum(dnbinom(as.numeric(tb_deaths[i,]),mu=as.numeric(V2[i,])*1e6,size=sgsq,log=T))*wts[50:65][i]
    l1 - adj_19  }

  ### ### ### TOTAL TB DEATHS 1999-2014  ### ### ### ### ### ### ~~~*
  # Motivation: norm, mean centered with CI = +/- 5% of mean
  tb_deaths_tot   <- rowSums(CalibDat[["tb_deaths"]][,-1])
  adj_19a         <- sum(dnorm(tb_deaths_tot,tb_deaths_tot,tb_deaths_tot*0.25/1.96,log=T)*wts[50:65])
  tb_dth_tot_lLik <- function(V) {
    sum(dnorm(tb_deaths_tot,rowSums(V)*1e6,tb_deaths_tot*0.2/1.96,log=T)*wts[50:65]) - adj_19a  }

  ### ### ### TB DEATHS AGE DISTRIBUTION 1999-2014  ### ### ### ### ### ###
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  tb_deaths_age  <- CalibDat[["tb_deaths"]][,-1]
  adj_19b        <- sum(dDirMult(M=tb_deaths_age+0.1,n=tb_deaths_age+0.1,Rho=0.01)*wts[50:65])
  tb_dth_age_lLik <- function(V,rho=0.01) { # V = table of deaths by age 1999-2014 (row=16 years, col=11 ages)
    V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
    sum(dDirMult(M=V2,n=tb_deaths_age,Rho=rho)*wts[50:65]) - adj_19b  }

### ### ### HIV DEATHS 2008-2012 ### ### ### ### ### ###
# Motivation: overdispersed poisson, modelled with negbin with overdispersion param = 100
  hiv_deaths      <- CalibDat[["hiv_deaths"]][,-1]
  adj_20          <- 0
  for(i in 1:16) adj_20 <- adj_20 + sum(dnbinom(as.numeric(hiv_deaths[i,]),mu=as.numeric(hiv_deaths[i,]),size=50,log=T))*wts[50:65][i]
  hiv_deaths_lLik <- function(V,sgsq=50) { # V =  HIV deaths by age 2008-2012  (row=5 years, col= 7 ages)
    V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
    l1 <- 0
    for(i in 1:16) l1 <- l1 + sum(dnbinom(as.numeric(hiv_deaths[i,]),mu=as.numeric(V2[i,])*1e6,size=sgsq,log=T))*wts[50:65][i]
    l1 - adj_20  }

  ### ### ### TOTAL HIV DEATHS 1999-2014  ### ### ### ### ### ### ~~~*
  # Motivation: norm, mean centered with CI = +/- 5% of mean
  hv_deaths_tot   <- rowSums(CalibDat[["hiv_deaths_all"]][,-1])
  adj_20a         <- sum(dnorm(hv_deaths_tot,hv_deaths_tot,hv_deaths_tot*0.25/1.96,log=T)*wts[61:64])
  hv_dth_tot_lLik <- function(V) {
    sum(dnorm(hv_deaths_tot,rowSums(V)*1e6,hv_deaths_tot*0.2/1.96,log=T)*wts[61:64]) - adj_20a  }

  ### ### ### HIV DEATHS AGE DISTRIBUTION 1999-2014  ### ### ### ### ### ###
  # Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
  hv_deaths_age  <- CalibDat[["hiv_deaths_all"]][,-1]
  adj_20b        <- sum(dDirMult(M=hv_deaths_age,n=hv_deaths_age,Rho=0.01)*wts[61:64])
  hv_dth_age_lLik <- function(V,rho=0.01) { # V = table of deaths by age 1999-2014 (row=16 years, col=11 ages)
    V2 <- cbind(V[,1]+V[,2],V[,3:7],rowSums(V[,8:11]))
    sum(dDirMult(M=V2+0.5,n=hv_deaths_age+0.5,Rho=rho)*wts[61:64]) - adj_20b  }

### ### ### HIV PREVALENCE 2006-2012  ### ### ### ### ### ###
# Motivation: norm, mean centered with CI = 5% million
  hiv_prev_by_year      <- CalibDat[["hiv_prev_by_year"]][,2]
  adj_21                <-  sum(dnorm(hiv_prev_by_year,hiv_prev_by_year,hiv_prev_by_year*0.1/1.96,log=T)*wts[57:63])
  hiv_prev_by_year_lLik <- function(V) { # V = total HIV pop 2006-2012
    sum(dnorm(hiv_prev_by_year,V*1e6,hiv_prev_by_year*0.1/1.96,log=T)*wts[57:63]) - adj_21  }

### ### ###  HIV AGE DISTRIBUITION 2011  ### ### ### ### ### ###
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 100
  hiv_prev_by_age_11      <- CalibDat[["hiv_prev_by_age_11"]][,2]/sum(CalibDat[["hiv_prev_by_age_11"]][,2])
  adj_22                  <-  sum(log(hiv_prev_by_age_11)*hiv_prev_by_age_11)*100
  hiv_prev_by_age_11_lLik <- function(V,ESS=100) { # V =  HIV pop in 2011 (11 ages)
    sum(log(V[2:7]/sum(V[2:7]))*hiv_prev_by_age_11)*ESS - adj_22  }

### ### ### ART VOLUME 2010  ### ### ### ### ### ###
# Motivation: norm, mean centered with CI = +/- 5% of mean
  art_vol_10      <- CalibDat[["art_vol_10"]]
  adj_23          <- dnorm(art_vol_10,art_vol_10,art_vol_10*0.10/1.96,log=T)
  art_vol_10_lLik <- function(V) { # V = ART vol in 2010 (scalar)
    dnorm(art_vol_10,V*1e6,art_vol_10*0.1/1.96,log=T) - adj_23   }

  ### ### ### HOMELESS POP 2010  ### ### ### ### ### ###
  # Motivation: norm, mean centered with CI = +/- 25% of mean
  homeless_pop      <- CalibDat[["homeless_pop"]][1]
  adj_23b          <- dnorm(homeless_pop,homeless_pop,homeless_pop*0.25/1.96,log=T)
  homeless_10_lLik <- function(V) { # V = homeless pop in 2010 (scalar)
    dnorm(homeless_pop,V,homeless_pop*0.25/1.96,log=T) - adj_23b   }

### ### ### LIKELIHOOD FOR BORGDORFF ESTIMATES  ### ### ### ### ### ###
  ss_borgdorff   <- 854.5463/4
  datB           <- CalibDat[["borgdorff_data"]]
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
  datF         <- CalibDat[["ferebee_data"]]
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
  datS            <- CalibDat[["sutherland_data"]]
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
  adj_28       <- dbeta(CalibDat[["pct_sp_ss"]][1],CalibDat[["pct_sp_ss"]][1]*100,(1-CalibDat[["pct_sp_ss"]][1])*100,log=T)
  smr_pos_lLik <- function(V) { # V = fraction smr-pos in 1950 (scalar)
    dbeta(CalibDat[["pct_sp_ss"]][1],V*100,(1-V)*100,log=T) - adj_28  }

### ### ### HIV SURVIVAL ### ### ### ### ### ###
  hiv_surv <- CalibDat[["hiv_survival"]]
  adj_29  <- sum(dnorm(hiv_surv[,2],hiv_surv[,2],(hiv_surv[,4]-hiv_surv[,3])/1.96*2,log=T))
  hiv_surv_lLik <- function(Par) { # Par= c(muH2,mHxtoHy[1])
    sv <- (1/(Par[[1]][,2]+Par[[2]][,1]+Par[[3]])+Par[[2]][,1]/(Par[[1]][,2]+Par[[2]][,1]+Par[[3]])*1/(Par[[1]][,4]+Par[[3]]))[3:8]/12
    sum(dnorm(sv,hiv_surv[,2],(hiv_surv[,4]-hiv_surv[,3])/1.96*2,log=T)) - adj_29  }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


