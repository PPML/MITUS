#######  ORGANIZE DATA FOR CALIBRATION  #######

# homelessness prevalence by age

#XXX  FROM OTIS: Routine reporting of new TB cases, subdivided by age, foreign-born status, prior treatment, and HIV status (though some biases in the HIV values).  
#XXX Separately, the fraction of all foreign-born TB occurring within 2 years, from Suzanne’s results (single number).  
#XXX	Distribution of reported TB cases across drug resistance categories (INH resistance, MDR-TB), by foreign-born status and treatment history, from the 2013 TB report. 
#XXX	Mortality while on treatment, from the 2013 TB report. 
#XXX	A crude likelihood for the total volume of TLTBI in 2015 – estimated at 291,000–433,000 (Sterling 2006)
#XXX	Cross-sectional data from NHANES 2011 (Miramontes 2015) describing the distribution of latent infection across population groups. Am just using the IGRA data
#XXX	Time series estimates of total population size from decennial censuses.
#XXX	Estimate of total immigrant population from census and CPS.
#XXX	Death registry data on HIV and TB deaths among HIV positive individuals deaths (ICD-10 B20-24) and TB deaths among HIV negative individuals (ICD-10 pulmonary A16… or all TB? A16-19). 
#XXX     (National Center for Health Statistics 2014) http://wonder.cdc.gov/mortSQL.html 
#XXX	Recent HIV prevalence estimates (to ensure the HIV incidence trends, in combination with HIV mortality estimates, are producing appropriate results)
#XXX	Crude estimate of 2010 ART volume = 480,395 (S. Cohen et al. 2011). This value based on synthesizing multiple datasets, and allowing for imperfect performance in diagnosis, 
#XXX     linkage to care, and retention. uncertainty not given. 
#XXX	Information on progression risk as a function of time since infection, from Ferebee, Sutherland and Borgdorff articules
#XXX	Estimates of TB case fatality and duration of disease from Tiemersma.
#XXX	Estimates of smear-positivity in early years from Styblo
#XXX  Estimates of HIV survival in the absence of treatment from CASCADE Collaboration 2000
#XXX  Estimates of size of homeless pop  



CalibDat <- list()
############
## TB NOTIFICATIONS DATA BY YEAR, AGE, NATIVITY, PRIOR TX, AND HIV STATUS
############

# 1. total by year up to 2015 *************
  tot_cases  <- read.csv("calib files/Total_cases_4-1-16.csv")
  tot_cases2 <- tot_cases
  shift75    <- tot_cases[23,2]*2-tot_cases[22,2]-tot_cases[24,2]
  tot_cases2[1:22,2] <- tot_cases[1:22,2]+shift75
  CalibDat[["tot_cases"]] <- tot_cases2

# 2. age distribution by year up to 2014
  age_cases0         <- read.csv("calib files/age_cases_10-21-15.csv")
  age_cases_us0      <- read.csv("calib files/age_cases_us_3-1-16.csv")[,-1]
  age_cases_fb0      <- age_cases_us0
  age_cases_fb0[,-1] <- age_cases0[,-1] - age_cases_us0[,-1]

  CalibDat[["age_cases"]]    <- data.frame(year = age_cases0[,1],
                                           age_cases0[,-c(1,12)]/age_cases0[,12],
                                           sample_size= age_cases0[,12])

  CalibDat[["age_cases_us"]] <- data.frame(year = age_cases_us0[,1],
                                           age_cases_us0[,-c(1,12)]/age_cases_us0[,12],
                                           sample_size= age_cases_us0[,12])

  CalibDat[["age_cases_fb"]] <- data.frame(year = age_cases_fb0[,1],
                                           age_cases_fb0[,-c(1,12)]/age_cases_fb0[,12],
                                           sample_size= age_cases_fb0[,12])

# 3. nativity by year up to 2015   *************
  fb_cases0 <- read.csv("calib files/fb_cases_4-1-16.csv")
  CalibDat[["fb_cases"]] <- data.frame(year        = fb_cases0[,1],
                                       pct_fb      = fb_cases0[,3]/fb_cases0[,4],
                                       sample_size = fb_cases0[,4])
  
  CalibDat[["fbus_cases_slope5"]] <- apply(log(fb_cases0[20:23,3:2]),2,function(x) lm(x~I(1:4))$coef[2])

# 3. prior tx by year up to 2014
  pt_cases0 <- read.csv("calib files/prevTB_cases_10-21-15.csv")
  CalibDat[["prev_cases"]] <- data.frame(year        = pt_cases0[,1],
                                         pct_prev    = pt_cases0[,3]/(pt_cases0[,2]+pt_cases0[,3]),
                                         sample_size = pt_cases0[,2]+pt_cases0[,3])

# 4.HIV by year up to 2014
  hiv_cases0 <- read.csv("calib files/HIV_cases_10-21-15.csv")
  CalibDat[["hiv_cases"]] <- data.frame(year        = hiv_cases0$year,
                                        pct_hiv     = hiv_cases0[,3]/hiv_cases0[,2],
                                        sample_size = hiv_cases0[,2])

  CalibDat[["hiv_cases2"]] <- data.frame(year        = hiv_cases0$year,
                                         known_pos   = hiv_cases0[,3],
                                         known_neg   = hiv_cases0[,2]-hiv_cases0[,3],
                                         unknown     = tot_cases2[41:62,2]-hiv_cases0[,2])

# Homeless by year up to 2014
  us_hr_cases0 <- read.csv("calib files/us_hr_cases_3-9-16.csv")
  CalibDat[["us_homeless_cases"]] <- data.frame(year        = us_hr_cases0$year,
                                                   pct_hr      = us_hr_cases0[,2]/us_hr_cases0[,3],
                                                   sample_size = us_hr_cases0[,3])

# Homeless by year up to 2014
  hr_cases0 <- read.csv("calib files/Homeless_cases_10-21-15.csv")
  CalibDat[["homeless_cases"]] <- data.frame(year        = hr_cases0$year,
                                             pct_hr      = hr_cases0[,2]/hr_cases0[,3],
                                             sample_size = hr_cases0[,3])

# 5. Fb >2 yrs overall
  fb_rec_cases0 <- read.csv("calib files/FB_recent_cases_10-21-15.csv")
  CalibDat[["fb_recent_cases"]] <- data.frame(year        = fb_rec_cases0$year,
                                              pct_fb_rec  = fb_rec_cases0[,3]/fb_rec_cases0[,2],
                                              sample_size = fb_rec_cases0[,2])

# 6. Pct INH by year, fb, and tx history up to 2014
  inh_cases0 <- read.csv("calib files/inh_res_cases_10-21-15.csv")
  inh_cases0_us_n <- inh_cases0[inh_cases0$prior_tx==0 & inh_cases0$foreign_born==0, ]
  CalibDat[["inh_res_us_n"]] <- data.frame(year        = inh_cases0_us_n$year,
                                           pct_inh     = inh_cases0_us_n[,5]/inh_cases0_us_n[,4],
                                           sample_size = inh_cases0_us_n[,4])
  inh_cases0_us_e <- inh_cases0[inh_cases0$prior_tx==1 & inh_cases0$foreign_born==0, ]
  CalibDat[["inh_res_us_e"]] <- data.frame(year        = inh_cases0_us_e$year,
                                           pct_inh     = inh_cases0_us_e[,5]/inh_cases0_us_e[,4],
                                           sample_size = inh_cases0_us_e[,4])
  inh_cases0_fb_n <- inh_cases0[inh_cases0$prior_tx==0 & inh_cases0$foreign_born==1, ]
  CalibDat[["inh_res_fb_n"]] <- data.frame(year        = inh_cases0_fb_n$year,
                                           pct_inh     = inh_cases0_fb_n[,5]/inh_cases0_fb_n[,4],
                                           sample_size = inh_cases0_fb_n[,4])
  inh_cases0_fb_e <- inh_cases0[inh_cases0$prior_tx==1 & inh_cases0$foreign_born==1, ]
  CalibDat[["inh_res_fb_e"]] <- data.frame(year        = inh_cases0_fb_e$year,
                                           pct_inh     = inh_cases0_fb_e[,5]/inh_cases0_fb_e[,4],
                                           sample_size = inh_cases0_fb_e[,4])
# 7. Pct MDR by year and tx history up to 2014
  mdr_cases0 <- read.csv("calib files/mdr_res_cases_10-21-15.csv")
  mdr_cases0_us_n <- mdr_cases0[mdr_cases0$prior_tx==0 & mdr_cases0$foreign_born==0, ]
  CalibDat[["mdr_res_us_n"]] <- data.frame(year        = mdr_cases0_us_n$year,
                                           pct_mdr     = mdr_cases0_us_n[,5]/mdr_cases0_us_n[,4],
                                           sample_size = mdr_cases0_us_n[,4])
  mdr_cases0_us_e <- mdr_cases0[mdr_cases0$prior_tx==1 & mdr_cases0$foreign_born==0, ]
  CalibDat[["mdr_res_us_e"]] <- data.frame(year        = mdr_cases0_us_e$year,
                                           pct_mdr     = mdr_cases0_us_e[,5]/mdr_cases0_us_e[,4],
                                           sample_size = mdr_cases0_us_e[,4])
  mdr_cases0_fb_n <- mdr_cases0[mdr_cases0$prior_tx==0 & mdr_cases0$foreign_born==1, ]
  CalibDat[["mdr_res_fb_n"]] <- data.frame(year        = mdr_cases0_fb_n$year,
                                           pct_mdr     = mdr_cases0_fb_n[,5]/mdr_cases0_fb_n[,4],
                                           sample_size = mdr_cases0_fb_n[,4])
  mdr_cases0_fb_e <- mdr_cases0[mdr_cases0$prior_tx==1 & mdr_cases0$foreign_born==1, ]
  CalibDat[["mdr_res_fb_e"]] <- data.frame(year        = mdr_cases0_fb_e$year,
                                           pct_mdr     = mdr_cases0_fb_e[,5]/mdr_cases0_fb_e[,4],
                                           sample_size = mdr_cases0_fb_e[,4])
############
## TB OUTCOMES BY YEAR up to 2012
############
  tx_outcomes0 <- read.csv("calib files/tx_outcomes_10-21-15.csv")
  CalibDat[["tx_outcomes"]] <- data.frame(year        = tx_outcomes0$year,
                                        pct_disc     = tx_outcomes0[,4]/tx_outcomes0[,2],
                                        pct_dead     = tx_outcomes0[,5]/tx_outcomes0[,2],
                                        sample_size = tx_outcomes0$all)

############
## TLTBI VOLUME AND DISTRIBUTION
############
  CalibDat[["TLTBI_volume"]] <- c(sqrt(291000*433000),c(291000,433000))

  TLTBI_dist <- c(0.562,0.048,0.021); names(TLTBI_dist) <- c("FB","HR","HV")
  CalibDat[["TLTBI_dist"]] <- TLTBI_dist

############
## LTBI PREVALENCE, NHANES 1999+2011
############
 # load("calib files/nhanes_data_10-21-15.rData") # nhanes_dat

 # CalibDat[["LTBI_prev_US_99"]] <- nhanes_dat[["US_99"]][,c(1,3,5)]
 # CalibDat[["LTBI_prev_FB_99"]] <- nhanes_dat[["FB_99"]][,c(1,3,5)]
 # CalibDat[["LTBI_prev_US_11"]] <- nhanes_dat[["US_11"]][,c(1,3,5)]
 # CalibDat[["LTBI_prev_FB_11"]] <- nhanes_dat[["FB_11"]][,c(1,3,5)]

  load("calib files/IgraDat_1-13-2016.rData") # IgraDat

  CalibDat[["LTBI_prev_US_11_IGRA"]] <- IgraDat[IgraDat$nativity=="us",c(3,9:10)]
  CalibDat[["LTBI_prev_FB_11_IGRA"]] <- IgraDat[IgraDat$nativity=="fb",c(3,9:10)]
  
  load("calib files/DoubPosDat_9-14-2016.rData") # DoubPosDat
  
  CalibDat[["LTBI_prev_US_11_DoubPos"]] <- DoubPosDat[DoubPosDat$nativity=="us",c(3,9:10)]
  CalibDat[["LTBI_prev_FB_11_DoubPos"]] <- DoubPosDat[DoubPosDat$nativity=="fb",c(3,9:10)]

############
## TOTAL POP
############
  pop_year_fb0 <- read.csv("calib files/tot_by_year_fb_10-21-15.csv")
  pop_age_fb0 <- read.csv("calib files/tot2014_by_age_fb_10-21-15.csv")
  CalibDat[["tot_pop_yr_fb"]]   <- pop_year_fb0
  CalibDat[["tot_pop14_ag_fb"]] <- pop_age_fb0

############
## TB DEATHS
############
  deaths0 <- read.csv("calib files/icd10_mcd_mort_4-1-16.csv")
  deaths0$TB_code <- deaths0$Cause.of.death.Code=="tb"
  deaths0$HIV_code <- deaths0$Cause.of.death.Code=="hiv"

  tb_deaths0 <- matrix(NA,length(1999:2014),length(unique(deaths0$Age.Group)))
  colnames(tb_deaths0) <- unique(deaths0$Age.Group)
  rownames(tb_deaths0) <- 1999:2014
  hiv_deaths0 <- tb_deaths0

  for(i in 1:nrow(tb_deaths0)) { for(j in 1:ncol(tb_deaths0)) {
    idxtb <- deaths0$Year==(1999:2014)[i] & deaths0$Age.Group==unique(deaths0$Age.Group)[j] & deaths0$TB_code
    idxhiv <- deaths0$Year==(1999:2014)[i] & deaths0$Age.Group==unique(deaths0$Age.Group)[j] & deaths0$HIV_code
    tb_deaths0[i,j] <- sum(deaths0[idxtb,3]); hiv_deaths0[i,j] <- sum(deaths0[idxhiv,3]);  }  }
  
  tb_deaths1 <- data.frame(year   = 1999:2014,
                           "0_4"   = rowSums(tb_deaths0[,1:2]),
                           "5_14"  = tb_deaths0[,3],
                           "15_24" = tb_deaths0[,4],
                           "25_34" = tb_deaths0[,5],
                           "35_44" = tb_deaths0[,6],
                           "45_54" = tb_deaths0[,7],
                           "55_64" = tb_deaths0[,8],
                           "65_74" = tb_deaths0[,9],
                           "75_84" = tb_deaths0[,10],
                           "85p"   = tb_deaths0[,11] )

hiv_deaths1 <- data.frame(year   = 1999:2014,
                         "0_4"   = rowSums(hiv_deaths0[,1:2]),
                         "5_14"  = hiv_deaths0[,3],
                         "15_24" = hiv_deaths0[,4],
                         "25_34" = hiv_deaths0[,5],
                         "35_44" = hiv_deaths0[,6],
                         "45_54" = hiv_deaths0[,7],
                         "55_64" = hiv_deaths0[,8],
                         "65_74" = hiv_deaths0[,9],
                         "75_84" = hiv_deaths0[,10],
                         "85p"   = hiv_deaths0[,11] )

  CalibDat[["tb_deaths"]]  <- tb_deaths1
  CalibDat[["hiv_deaths"]] <- hiv_deaths1

############
## HIV DEATHS
############
  deathsHiv <- read.csv("calib files/hiv_mort_4-15-16.csv")[,-9]
  CalibDat[["hiv_deaths_all"]]  <- deathsHiv

############
## HIV PREVALENCE
############
  hiv_prev_year0 <- read.csv("calib files/Calib_HIV_by_year.csv")
  hiv_prev_age0 <- read.csv("calib files/Calib_HIV_by_age_2011.csv")

  CalibDat[["hiv_prev_by_year"]]   <- hiv_prev_year0
  CalibDat[["hiv_prev_by_age_11"]] <- hiv_prev_age0

############
## ART Volume
############
  CalibDat[["art_vol_10"]]   <- 480395

############
## TB progression with time
############
  borgdorff_data  <- read.csv("Borgorff 2011 data.csv")
  ferebee_data <- read.csv("Ferebee 1970 data_rev2.csv")
  sutherland_data <- read.csv("Sutherland 1968 data_rev2.csv")[,-3]
  CalibDat[["borgdorff_data"]]   <- borgdorff_data
  CalibDat[["ferebee_data"]]     <- ferebee_data
  CalibDat[["sutherland_data"]]  <- sutherland_data

############
## DURATION OF DISEASE AND CASE FATALITY
############
  CalibDat[["duration_of_dis"]]    <- 3.0
  CalibDat[["case_fatality_sn"]]   <- 0.2
  CalibDat[["case_fatality_sp"]]   <- 0.7

############
## SMEAR-POSITIVITY
############
  CalibDat[["pct_sp_ss"]]   <- c(148/383,383)

############
## HIV-SURVIVAL
############
  hiv_survival  <- read.csv("calib files/hiv_survival_age_10-21-15.csv")
  CalibDat[["hiv_survival"]]   <- hiv_survival

############
## HIV-SURVIVAL
############
CalibDat[["homeless_pop"]]   <- c(1.6,1.1,2.1) # Per 

############
## RECENCY WEIGHT
############
LgtCurveY <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
                                          zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2015-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
                                          zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }
ImptWeights <- LgtCurveY(1990,2015,0.90)+0.10
names(ImptWeights) <- 1950:2015
ImptWeights <- ImptWeights/max(ImptWeights)

CalibDat[["ImptWeights"]]   <- ImptWeights

############
## SAVE IT!
############

save(CalibDat,file="CalibDat_9-14-16.rData")


# 
# pdfnam <- "Fig_Importance_Weights.pdf"
# pdf(file=pdfnam,width=8, height=4) 
# par(mfrow=c(1,1),mar=c(4.0,4,0.5,0.5))
#  plot(1950:2015,ImptWeights,type="l",ylim=0:1,col=4,lwd=2,xlab="Year",ylab="Importance Weights",las=1)
#  abline(h=0:10/10,col="grey90")
# lines(1950:2015,ImptWeights,type="l",ylim=0:1,col=4,lwd=2)
# ###
#  dev.off(); system(paste("open", pdfnam)) # code for ma
