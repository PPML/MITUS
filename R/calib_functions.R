#'Calibration Functions for use with the tb_model.cpp model
#'This script creates several individual log likelihood functions
#'for the calibration of the State Level TB model in tb_model.cpp
#'These llikelihood functions are called in IMIS_functions.R
#'takes in the outputs and calibration data and creates likelihood functions

#'Total Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@param V vector of total notifications 1953-2014
#'@return likelihood

notif_tot_lik <- function(V) {
  notif_tot     <- CalibDat[["tot_cases"]][,2]
  adj_1         <- sum(dnorm(notif_tot,notif_tot,notif_tot*0.1/1.96,log=T)*wts[4:67])
  sum(dnorm(notif_tot,V,notif_tot*0.1/1.96,log=T)*wts[4:67]) - adj_1  }

#'US Cases Age Distribution 1993-2013
#'Motivation: dirichlet-multinomial data with additional non-sampling biases
#'@param V table of us notifications by age 1993-2013 (row=21 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
notif_age_us_lLik <- function(V,rho=0.005) {
  notif_age         <- CalibDat[["age_cases"]][,-c(1,12)]*CalibDat[["age_cases"]][,12]
  notif_age_us      <- CalibDat[["age_cases_us"]][,-c(1,12)]*CalibDat[["age_cases_us"]][,12]
  #weighted sum across the years
  adj_2a            <- sum(dDirMult(M=notif_age_us,n=notif_age_us,Rho=0.005)*wts[44:67])
  sum(dDirMult(M=V,n=notif_age_us,Rho=rho)*wts[44:67]) - adj_2a
  }

#' FB CASES AGE DISTRIBUTION 1993-2013
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of fb notifications by age 1993-2013 (row=21 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
notif_age_fb_lLik <- function(V,rho=0.005) {
  notif_age_fb     <- CalibDat[["age_cases_fb"]][,-c(1,12)]*CalibDat[["age_cases_fb"]][,12]
  adj_2b           <- sum(dDirMult(M=notif_age_fb,n=notif_age_fb,Rho=0.005)*wts[44:67])
  sum(dDirMult(M=V,n=notif_age_fb,Rho=rho)*wts[44:67]) - adj_2b
}

#' CASES FB DISTRIBUTION 1993-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of notifications by fb 1993-2014 (row=22 years, col=fb then us)
#'@param rho correlation parameter
#'@return likelihood
notif_fb_lLik <- function(V,rho=0.005) {
  notif_fb      <- cbind(CalibDat[["fb_cases"]][,2],1-CalibDat[["fb_cases"]][,2])*CalibDat[["fb_cases"]][,3]
  adj_3         <- sum(dDirMult(M=notif_fb,n=notif_fb,Rho=0.005)*wts[44:67])
  sum(dDirMult(M=V,n=notif_fb,Rho=rho)*wts[44:67]) - adj_3  }

#' CASES FB DISTRIBUTION SLOPES OVER PAST 4 year
#'@param V table of notifications by fb 1993-2014 (row=22 years, col=fb then us)
#'@return likelihood
notif_fbus_slp_lLik <- function(V) {
  notif_fbus_slp5     <- CalibDat[["fbus_cases_slope5"]];
  adj_3a         <- sum(dnorm(notif_fbus_slp5,notif_fbus_slp5,0.01/2,log=T))
  V2 <- apply(log(V[20:23,]),2,function(x) lm(x~I(1:4))$coef[2])
  sum(dnorm(V2,notif_fbus_slp5,0.01/2,log=T)) - adj_3a  }

#' CASES HR DISTRIBUTION 1993-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of notifications by hr 1993-2014 (row=22 years, col=pos then neg)
#'@param rho correlation parameter
#'@return likelihood
notif_us_hr_lLik <- function(V,rho=0.005) {
  notif_us_hr      <- cbind(CalibDat[["us_homeless_cases"]][,2],1-CalibDat[["us_homeless_cases"]][,2])*CalibDat[["us_homeless_cases"]][,3]
  adj_5b           <- sum(dDirMult(M=notif_us_hr,n=notif_us_hr,Rho=0.005)*wts[44:67])
  sum(dDirMult(M=V,n=notif_us_hr,Rho=rho)*wts[44:67]) - adj_5b  }

#' CASES FB RECENT ENTRY DISTRIBUTION 1993-2013
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of notifications by FB 1993-2014 (row=22 years, col=pos then neg)
#'@param rho correlation parameter
#'@return likelihood
notif_fb_rec_lLik <- function(V,rho=0.005) {
  notif_fb_rec   <- cbind(CalibDat[["fb_recent_cases"]][,2],1-CalibDat[["fb_recent_cases"]][,2])*CalibDat[["fb_recent_cases"]][,3]
  adj_6          <- sum(dDirMult(M=notif_fb_rec,n=notif_fb_rec,Rho=0.005)*wts[44:65])
  sum(dDirMult(M=V,n=notif_fb_rec,Rho=rho)*wts[44:65]) - adj_6  }

#' TREATMENT OUTCOMES 1993-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of treatment outcomes 1993-2012 (row=20 years, col= complete, discontinue, dead)
#'@param rho correlation parameter
#'@return likelihood
tx_outcomes_lLik <- function(V,rho=0.005) {
  tx_outcomes      <- cbind(1-rowSums(CalibDat[["tx_outcomes"]][,2:3]),CalibDat[["tx_outcomes"]][,2],CalibDat[["tx_outcomes"]][,3])*CalibDat[["tx_outcomes"]][,4]
  adj_11           <- sum(dDirMult(M=tx_outcomes,n=tx_outcomes,Rho=0.005)*wts[44:65])
  sum(dDirMult(M=V,n=tx_outcomes,Rho=rho)*wts[44:65]) - adj_11
  }

#' TOTAL LTBI TREATMENT INITS 2002
#' Motivation: norm, mean centered with CI = +/- 10% of mean
#'@param V scalar total TLTBI inits in 2002
#'@return likelihood
tltbi_tot_lLik   <- function(V) {
  tltbi_vol        <- CalibDat[["TLTBI_volume"]]
  adj_12           <- dnorm(tltbi_vol[1],tltbi_vol[1],diff(tltbi_vol[2:3])/2/1.96,log=T)
  dnorm(tltbi_vol[1],V*1e6,diff(tltbi_vol[2:3])/2/1.96,log=T) - adj_12
  }

#' DISTRIBUTION OF LTBI TREATMENT INITS 2002
#'@param V dist TLTBI inits in 2002 (vector fraction FB, HR in 2002)
#'@return likelihood
tltbi_dist_lLik  <- function(V) {
  TLTBI_dist       <- CalibDat[["TLTBI_dist"]]
  TLTBI_dist       <-TLTBI_dist[-3]
  adj_13          <- sum( dbeta(TLTBI_dist,TLTBI_dist*100,(1-TLTBI_dist)*100,log=T) )
  sum( dbeta(TLTBI_dist,V*100,(1-V)*100,log=T) ) - adj_13  }

#' LTBI PREVALENCE BY AGE 2011, US
#' Motivation: additional prior on LTBI, using beta densities parameterized to Miramontes/Hill results
#'@name ltbi_us_11_lLik
#'@param V LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
#'@return likelihood
ltbi_us_11_lLik <- function(V) {
  ltbi_us_11      <- CalibDat[["LTBI_prev_US_11_IGRA"]]
  adj_15          <- sum( dbeta(ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) ) - adj_15)*2  }

#'Likelihood function for double positive LTBI in the US born based on 2011 NHANES

#'@name ltbi_us_11_dp_lLik
#'@param V LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
#'@return likelihood
ltbi_us_11_dp_lLik <- function(V) {
  ltbi_us_11_dp      <- CalibDat[["LTBI_prev_US_11_DoubPos"]]
  adj_15dp           <- sum(dbeta(ltbi_us_11_dp[,2]/rowSums(ltbi_us_11_dp[,2:3]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) ) - adj_15dp)*2  }

#' LTBI PREVALENCE BY AGE 2011, FB
#' Motivation: multinomial adjusted to match effective sample size due to survey weighting
#'@param V LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
#'@return likelihood
ltbi_fb_11_lLik <- function(V) {
  ltbi_fb_11      <- CalibDat[["LTBI_prev_FB_11_IGRA"]]
  adj_16          <- sum( dbeta(ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) ) - adj_16)*2  }

#'Likelihood function for double positive LTBI in the NUS born based on 2011 NHANES

#'@name ltbi_fb_11_dp_lLik
#'@param V LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
#'@return likelihood
ltbi_fb_11_dp_lLik <- function(V) {
  ltbi_fb_11_dp      <- CalibDat[["LTBI_prev_FB_11_DoubPos"]]
  adj_16dp           <- sum( dbeta(ltbi_fb_11_dp[,2]/rowSums(ltbi_fb_11_dp[,2:3]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) ) - adj_16dp)*2  }

#' TOTAL POP EACH DECADE, FOR FB
#' Motivation: norm, mean centered with CI = +/- 2 million wts[1+0:6*10]
#'@param V total pop (rows=year, cols=us, fb)
#'@return likelihood
tot_pop_yr_fb_lLik <- function(V) {
  tot_pop_yr_fb      <- CalibDat[["tot_pop_yr_fb"]]
  adj_17             <- sum(dnorm(tot_pop_yr_fb[-1,4],tot_pop_yr_fb[-1,4],tot_pop_yr_fb[7,4]*0.1/1.96,log=T)*wts[1+1:6*10])
  sum(dnorm(tot_pop_yr_fb[-1,4],V[c(11,21,31,41,51,61)],tot_pop_yr_fb[7,4]*0.1/1.96,log=T)*wts[1+1:6*10]) - adj_17  } # CI = +/- 2mil

#' #' TOTAL POP AGE DISTRIBUTION 2014 by nativity
#' #' Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
#' #'@param V US pop in 2016 (row=11 ages, col= rec fb, fb)
#' #'@param ESS explained sum of squares
#' #'@return likelihood
tot_pop16_ag_fb_lLik <- function(V,ESS=500) {
  tot_pop16_ag_fb      <- cbind(CalibDat[["tot_pop16_ag_fb"]][-9,3]/sum(CalibDat[["tot_pop16_ag_fb"]][-9,3]),
                                CalibDat[["tot_pop16_ag_fb"]][-9,4]/sum(CalibDat[["tot_pop16_ag_fb"]][-9,4]))
  adj_18               <- sum(log(tot_pop16_ag_fb[,1])*tot_pop16_ag_fb[,1])+sum(log(tot_pop16_ag_fb[,2])*tot_pop16_ag_fb[,2])
  V1 <- rbind(V[1,],V[2,]+V[3,],V[4,]+V[5,],V[6,],V[7,],V[8,],V[9,],V[10,]+V[11,])
  (sum(log(V1[,1]/sum(V1[,1]))*tot_pop16_ag_fb[,1])+sum(log(V1[,2]/sum(V1[,2]))*tot_pop16_ag_fb[,2]))*ESS - adj_18*ESS  }

# US_pop_tot_lLik_2<-function(X,ESS=500){
#   #target data
#   US_pop_tot  <- CalibDat[["tot_pop_yr_fb"]][-9,2]
#   #calculate an adjustment factor
#   adj_20s             <- sum(log(US_pop_tot[-1])*US_pop_tot[-1])
#   sum(log(X*US_pop_tot[-1])*ESS - adj_20s*ESS)
# }



#' TOTAL US POP No nativity (ONLY FOR USE IN DEMO MODEL)
#' 1970,1975,1980,1985,1990-2007
#' Motivation: norm, mean centered with CI = +/- 5% of mean
#'@param X vector of total deaths in US from 1971-2016, fraction of millions
#'@return j likelihood
US_pop_tot_lLik <- function(X) {
  # CalibDat$US_tot_mort <- read.csv(file="inst/extdata/US_total_mort.csv", header = FALSE)
  US_pop_tot  <- CalibDat[["tot_pop_yr_fb"]][-9,2];
  # X1<-X[c(11,21,31,41,51,61)]
  adj_20p             <- sum(dnorm(US_pop_tot[-1],US_pop_tot[-1],US_pop_tot[7]*0.01/1.96,log=T)*wts[1+1:6*10]);
 #this is the problematic line
  # print(X)
   j<-(sum(dnorm(US_pop_tot[-1], X ,US_pop_tot[7]*0.01/1.96,log=T)*wts[1+1:6*10]) - adj_20p);
   return(j)
   # } # CI = +/- 2mil *wts[1+1:6*10]) - adj_20p
   #dnorm(US_pop_tot[-1], X ,US_pop_tot[7]*0.1/1.96,log=T)*
}
#' #' TOTAL POP AGE DISTRIBUTION  NO NATIVITY (ONLY FOR USE IN DEMO MODEL) 1999-2014
#' #'@param V table of POP by age 2016 (row=16 years, col=11 ages)
#' #'@param ESS correlation parameter
#' #'@return likelihood
tot_pop_age_lLik <- function(V,ESS=500) {
  # CalibDat$US_mort_age <- read.csv(system.file("extdata","US_mort_age.csv", package="MITUS"))
  # tot_pop16_ag  <- CalibDat[["tot_pop16_ag_fb"]][-9,2]/sum(CalibDat[["tot_pop16_ag_fb"]][-9,2]);
  data("wonder_pop",package="MITUS")
  tot_pop16_ag<-wonder_pop/sum(wonder_pop)
  adj_18               <- sum(log(tot_pop16_ag)*tot_pop16_ag);
  sum(log(V[]/sum(V))*tot_pop16_ag[])*ESS - adj_18*ESS
   }


#' Total TB DEATHS 1999-2014
#' Motivation: overdispersed poisson, modelled with negbin with overdispersion param = 100 *wts[50:65]
#'@param V TB deaths by age 1999-2013  (row=15 years, col= 11 ages)
#'@return likelihood
tb_deaths_lLik <- function(V,sgsq=50) {
  tb_deaths      <- CalibDat[["tb_deaths"]][,-1]
  adj_19         <- 0
  for(i in 1:16) adj_19 <- adj_19 + sum(dnbinom(as.numeric(tb_deaths[i,]),mu=as.numeric(tb_deaths[i,]),size=50,log=T))*wts[50:65][i]
  V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
  l1 <- 0
  for(i in 1:16) l1 <- l1 + sum(dnbinom(as.numeric(tb_deaths[i,]),mu=as.numeric(V2[i,])*1e6,size=sgsq,log=T))*wts[50:65][i]
  l1 - adj_19  }

#' TOTAL TB DEATHS 1999-2014
#' Motivation: norm, mean centered with CI = +/- 5% of mean
#'@param V vector of TB deaths in fraction of millions
#'@return likelihood
tb_dth_tot_lLik <- function(V) {
  tb_deaths_tot   <- rowSums(CalibDat[["tb_deaths"]][,-1])
  adj_19a         <- sum(dnorm(tb_deaths_tot,tb_deaths_tot,tb_deaths_tot*0.2/1.96,log=T)*wts[50:65])
  sum(dnorm(tb_deaths_tot,V*1e6,tb_deaths_tot*0.2/1.96,log=T)*wts[50:65]) - adj_19a  }

#' TB DEATHS AGE DISTRIBUTION 1999-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of deaths by age 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
tb_dth_age_lLik <- function(V,rho=0.01) {
  tb_deaths_age  <- CalibDat[["tb_deaths"]][,-1]
  adj_19b        <- sum(dDirMult(M=tb_deaths_age+0.1,n=tb_deaths_age+0.1,Rho=0.01)*wts[50:65])
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  sum(dDirMult(M=V2,n=tb_deaths_age+.1,Rho=rho)*wts[50:65]) - adj_19b  }

#' TOTAL US DEATHS
#' 1970,1975,1980,1985,1990-2007
#' Motivation: norm, mean centered with CI = +/- 5% of mean
#'@param V vector of total deaths in US from 1971-2016, fraction of millions
#'@return likelihood

US_dth_tot_lLik <- function(V) {
  # CalibDat$US_tot_mort <- read.csv(file="inst/extdata/US_total_mort.csv", header = FALSE)
  US_deaths_tot   <- CalibDat[["US_tot_mort"]][22:38,-1]
  adj_20a         <- sum(dnorm(US_deaths_tot,US_deaths_tot,US_deaths_tot*0.01/1.96,log=T)*wts[51:67])
  sum(dnorm(US_deaths_tot,V,US_deaths_tot*0.01/1.96,log=T)*wts[51:67]) - adj_20a
  }

#' TOTAL DEATHS AGE DISTRIBUTION 1999-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of deaths by age 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
tot_dth_age_lLik <- function(V,ESS=500) {
  # CalibDat$US_mort_age <- read.csv(system.file("extdata","US_mort_age.csv", package="MITUS"))
  data("mort_ag_16", package = 'MITUS');
  mort_ag_16_d <-(mort_ag_16[]/sum(mort_ag_16[]))
  adj_20b               <- sum(log(mort_ag_16_d)*mort_ag_16_d);
  sum(log(V[]/sum(V))*mort_ag_16_d[])*ESS - adj_20b*ESS  }


#' Mortality Risk Group Distribution 1999-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V table of mort_dist 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
mort_dist_lLik <- function(V,rho=0.01) {
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


#'Homeless Population in 2010
#'Motivation: normally distributed, mean centered with CI = +/- 25% of mean
#'@param V scalar value of homeless pop in 2010
#'@return likelihood
homeless_10_lLik <- function(V) {
  homeless_pop      <- CalibDat[["homeless_pop"]][1]
  adj_23b          <- dnorm(homeless_pop,homeless_pop,homeless_pop*0.25/1.96,log=T)
  dnorm(homeless_pop,V,homeless_pop*0.25/1.96,log=T) - adj_23b   }
########################################################################################
#'Functions for likelihood of different published estimates
########################################################################################
# dist_gen <- matrix(1:16,4,4); dist_gen <- dist_gen/sum(dist_gen)
#' #'Likelihood of Borgdorff Estimates
#' #'across the tb_progression groups
#' #'average at the end
#' #'@param Par_list list of parameters Mpfast, Mrslow, rfast, rRecov
#' #'@param n_red crude reductions the likelihood (allows non-sampling bias)
#' #'@return likelihood
borgdorff_lLik <- function(Par_list,N_red=1) {  # Par_list = list(Mpfast[,c(1,3,2,4)], Mrslow[,c(1,3,2,4)], rfast, rRecov)
  ss_borgdorff   <- 854.5463/4
  datB           <- CalibDat[["borgdorff_data"]]
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

#'Likelihood of Ferebee Estimates
#'across the tb_progression groups
#'average at the end
#'@param Par_list list (pfast,rslow,rfast,rRecov)
#'@param n_red crude reductions the likelihood (allows non-sampling bias)
#'@return likelihood
ferebee_lLik <- function(Par_list,N_red=4) {
  datF         <- CalibDat[["ferebee_data"]]
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

#' #'Likelihood of Sutherland Estimates
#' #'across the tb_progression groups
#' #'average at the end
#' #'@param Par_list list (Mpfast,Mrslow,rfast,rRecov)
#' #'@param n_red crude reductions the likelihood (allows non-sampling bias)
#' #'@return likelihood
sutherland_lLik <- function(Par_list,N_red=4) {
  datS            <- CalibDat[["sutherland_data"]]
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

#'Likelihood of Tiemersa Estimates
#we decided on friday that there will be no comorbidity factor
#so this likelihood should be independent of the four risk groups
#since neither rSlfCur and muIp are independent of both the mort
#and the tb progression
#'@param Par vector (pfast,rslow,rfast,rRecov)
#'@return likelihood

tiemersma_lLik <- function(Par) { # Par= c(rSlfCur,muIp)
  adj_27         <- dnorm(3.0,3.0,0.5/1.96,log=T)+dbeta(0.9,0.9*50,(1-0.9)*50,log=T)
  rSlfCur <- Par[1]*12;
  muIp    <- Par[2]*12 #muIp is the value for all TB active disease (not just smear pos)

  dur <- (1/(rSlfCur+muIp))
  cf_sp <- muIp/(rSlfCur+muIp);

  l1 <- dnorm(3.0,dur,0.5/1.96,log=T)
  l2 <- dbeta(0.9,cf_sp*50,(1-cf_sp)*50,log=T)
  l1+l2 - adj_27
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

