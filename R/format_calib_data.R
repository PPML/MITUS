#' #'This script converts the calibration target data from
#' #'.csv files to a single .rData file
#' #'
#' #'CalibDat is a list
#' CalibDat <- list()
#'
#' #' TB Notifications
#' #' Total by year up until 2015
#' tot_cases  <- read.csv("inst/extdata/6.15.18 Final Data/Total_cases_6-15-18.csv")
#' tot_cases2 <- tot_cases
#' shift75    <- tot_cases[23,2]*2-tot_cases[22,2]-tot_cases[24,2]
#' tot_cases2[1:22,2] <- tot_cases[1:22,2]+shift75
#' CalibDat[["tot_cases"]] <- tot_cases2
#'
#' #'Age distribution by year up until 2014
#' age_cases0         <- read.csv("inst/extdata/6.15.18 Final Data/age_cases_6-15-18.csv")
#' age_cases_us0      <- read.csv("inst/extdata/6.15.18 Final Data/age_cases_us_6-15-18.csv")
#' age_cases_fb0      <- age_cases_us0
#' age_cases_fb0[,-1] <- age_cases0[,-1] - age_cases_us0[,-1]
#'
#' CalibDat[["age_cases"]]    <- data.frame(year = age_cases0[,1],
#'                                          age_cases0[,-c(1,12)]/age_cases0[,12],
#'                                          sample_size= age_cases0[,12])
#'
#' CalibDat[["age_cases_us"]] <- data.frame(year = age_cases_us0[,1],
#'                                          age_cases_us0[,-c(1,12)]/age_cases_us0[,12],
#'                                          sample_size= age_cases_us0[,12])
#'
#' CalibDat[["age_cases_fb"]] <- data.frame(year = age_cases_fb0[,1],
#'                                          age_cases_fb0[,-c(1,12)]/age_cases_fb0[,12],
#'                                          sample_size= age_cases_fb0[,12])
#'
#' #'Nativity by year up to 2015
#' fb_cases0 <- read.csv("inst/extdata/6.15.18 Final Data/fb_cases_6-15-18.csv")
#' CalibDat[["fb_cases"]] <- data.frame(year        = fb_cases0[,1],
#'                                      pct_fb      = fb_cases0[,3]/fb_cases0[,4],
#'                                      sample_size = fb_cases0[,4])
#'
#' CalibDat[["fbus_cases_slope5"]] <- apply(log(fb_cases0[20:23,3:2]),2,function(x) lm(x~I(1:4))$coef[2])
#'
#' #'US Homeless by year up to 2014
#' us_hr_cases0 <- read.csv("inst/extdata/6.15.18 Final Data/us_hr_cases_6-15-18.csv")
#' CalibDat[["us_homeless_cases"]] <- data.frame(year        = us_hr_cases0$year,
#'                                               pct_hr      = us_hr_cases0[,2]/us_hr_cases0[,3],
#'                                               sample_size = us_hr_cases0[,3])
#'
#' #'Homeless by year up to 2014
#' hr_cases0 <- read.csv("inst/extdata/6.15.18 Final Data/Homeless_cases_6-15-18.csv")
#' CalibDat[["homeless_cases"]] <- data.frame(year        = hr_cases0$year,
#'                                            pct_hr      = hr_cases0[,2]/hr_cases0[,3],
#'                                            sample_size = hr_cases0[,3])
#'
#' #'Long Term Resident (Non-US born >2 yrs) overall
#' # fb_rec_cases0 <- read.csv("calib files/FB_recent_cases_10-21-15.csv")
#' # CalibDat[["fb_recent_cases"]] <- data.frame(year        = fb_rec_cases0$year,
#' #                                             pct_fb_rec  = fb_rec_cases0[,3]/fb_rec_cases0[,2],
#' #                                             sample_size = fb_rec_cases0[,2])
#'
#' #'TB Treatment Outcomes
#' #'Overall, up to 2012
#' tx_outcomes0 <- read.csv("inst/extdata/6.15.18 Final Data/tx_outcomes_6-15-18.csv")
#' CalibDat[["tx_outcomes"]] <- data.frame(year        = tx_outcomes0$year,
#'                                         pct_disc     = tx_outcomes0[,4]/tx_outcomes0[,2],
#'                                         pct_dead     = tx_outcomes0[,5]/tx_outcomes0[,2],
#'                                         sample_size = tx_outcomes0$all)
#' #'TLTBI Volume and Distribution
#' CalibDat[["TLTBI_volume"]] <- c(sqrt(291000*433000),c(291000,433000))
#'
#' TLTBI_dist <- c(0.562,0.048,0.021); names(TLTBI_dist) <- c("FB","HR","HV")
#' CalibDat[["TLTBI_dist"]] <- TLTBI_dist
#'
#' #'LTBI Prevalence, NHANES 1999-2011
#' # load("calib files/IgraDat_1-13-2016.rData") # IgraDat
#' #
#' # CalibDat[["LTBI_prev_US_11_IGRA"]] <- IgraDat[IgraDat$nativity=="us",c(3,9:10)]
#' # CalibDat[["LTBI_prev_FB_11_IGRA"]] <- IgraDat[IgraDat$nativity=="fb",c(3,9:10)]
#' #
#' # load("calib files/DoubPosDat_9-14-2016.rData") # DoubPosDat
#' #
#' # CalibDat[["LTBI_prev_US_11_DoubPos"]] <- DoubPosDat[DoubPosDat$nativity=="us",c(3,9:10)]
#' # CalibDat[["LTBI_prev_FB_11_DoubPos"]] <- DoubPosDat[DoubPosDat$nativity=="fb",c(3,9:10)]
#'
#' #'Total Population
#' #'by Year
#' pop_year_fb0 <- read.csv("inst/extdata/6.15.18 Final Data/tot_by_year_fb_6-15-18.csv")
#' #'age distribution
#'
#' pop_age_fb0 <- read.csv("inst/extdata/6.15.18 Final Data/total2016_by_age_fb-6-15-18.csv")
#' CalibDat[["tot_pop_yr_fb"]]   <- pop_year_fb0
#' CalibDat[["tot_pop14_ag_fb"]] <- pop_age_fb0
#'
#' #'TB Deaths
#' deaths0 <- read.csv("inst/extdata/6.15.18 Final Data/icd10_mcd_mort_6-15-18.csv")
#' deaths0$TB_code <- deaths0$Cause.of.death.Code=="tb"
#' deaths0$HIV_code <- deaths0$Cause.of.death.Code=="hiv"
#'
#' tb_deaths0 <- matrix(NA,length(1999:2014),length(unique(deaths0$Age.Group)))
#' colnames(tb_deaths0) <- unique(deaths0$Age.Group)
#' rownames(tb_deaths0) <- 1999:2014
#' hiv_deaths0 <- tb_deaths0
#'
#' for(i in 1:nrow(tb_deaths0)) {
#'   for(j in 1:ncol(tb_deaths0)) {
#'   idxtb <- deaths0$Year==(1999:2014)[i] & deaths0$Age.Group==unique(deaths0$Age.Group)[j] & deaths0$TB_code
#'   idxhiv <- deaths0$Year==(1999:2014)[i] & deaths0$Age.Group==unique(deaths0$Age.Group)[j] & deaths0$HIV_code
#'
#'   tb_deaths0[i,j] <- sum(deaths0[idxtb,3]); hiv_deaths0[i,j] <- sum(deaths0[idxhiv,3]);
#' }  }
#'
#'   tb_deaths1 <- data.frame(year   = 1999:2014,
#'                            "0_4"   = rowSums(tb_deaths0[,1:2]),
#'                            "5_14"  = tb_deaths0[,3],
#'                            "15_24" = tb_deaths0[,4],
#'                            "25_34" = tb_deaths0[,5],
#'                            "35_44" = tb_deaths0[,6],
#'                            "45_54" = tb_deaths0[,7],
#'                            "55_64" = tb_deaths0[,8],
#'                            "65_74" = tb_deaths0[,9],
#'                            "75_84" = tb_deaths0[,10],
#'                            "85p"   = tb_deaths0[,11] )
#'
#'   CalibDat[["tb_deaths"]]  <- tb_deaths1
#'
#' #'TB progression with time based off of several studies
#' borgdorff_data  <- read.csv("inst/extdata/Borgorff 2011 data.csv")
#' ferebee_data <- read.csv("inst/extdata/Ferebee 1970 data_rev2.csv")
#' sutherland_data <- read.csv("inst/extdata/Sutherland 1968 data_rev2.csv")[,-3]
#' CalibDat[["borgdorff_data"]]   <- borgdorff_data
#' CalibDat[["ferebee_data"]]     <- ferebee_data
#' CalibDat[["sutherland_data"]]  <- sutherland_data
#'
#' #'Duration of Disease and Case Fatality
#' CalibDat[["duration_of_dis"]]    <- 3.0
#' CalibDat[["case_fatality_sn"]]   <- 0.2
#' CalibDat[["case_fatality_sp"]]   <- 0.7
#'
#' #'RECENCY WEIGHT
#' LgtCurveY <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
#' zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2015-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
#' zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }
#' ImptWeights <- LgtCurveY(1990,2015,0.90)+0.10
#' names(ImptWeights) <- 1950:2015
#' ImptWeights <- ImptWeights/max(ImptWeights)
#'
#' CalibDat[["ImptWeights"]]   <- ImptWeights
#'
#' #'Temporarily add in 5 data points from the CalibDat_9-14-2016.rData file
#' load("data/OldCalibDat_2018-06-28.rData")
#'
#' CalibDat[["LTBI_prev_US_11_IGRA"]] <- OldCalibDat[["LTBI_prev_US_11_IGRA"]]
#' CalibDat[["LTBI_prev_FB_11_IGRA"]] <- OldCalibDat[["LTBI_prev_FB_11_IGRA"]]
#' CalibDat[["LTBI_prev_US_11_DoubPos"]] <- OldCalibDat[["LTBI_prev_US_11_DoubPos"]]
#' CalibDat[["LTBI_prev_FB_11_DoubPos"]] <- OldCalibDat[["LTBI_prev_FB_11_DoubPos"]]
#' CalibDat[["fb_recent_cases"]]  <- OldCalibDat[["fb_recent_cases"]]
#'
#'
#' #'Save CalibDat to an .rData file
#' save(CalibDat,file=paste("CalibDat_", Sys.Date(),".rData",sep=""))
#'
