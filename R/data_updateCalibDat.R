####yearly update to the CalibDat
# year<-2020
# # ##### OTIS DATA
# total_cases<-7174
# us_cases<-2018
# percent_nus_cases<-(5127/7174)*100
# # make sure that each of these sum to 100
# age_cases_us<- c(2.03+5.50, 5.05, 10.60, 13.78, 10.21, 11.15, 18.83, 12.69, 6.89, 3.22)
# age_cases_nus<-c(0.04+0.21, 0.98, 9.26, 16.34, 15.17, 14.28, 16.21, 14.14, 9.50, 3.86)
# age_cases<-((age_cases_us/100*(1-percent_nus_cases/100)*total_cases)+(age_cases_nus/100*percent_nus_cases/100*total_cases))/total_cases
#
# # percent_recent_nus_cases
# ## a smoothed function
# p_0_1 = 9.73
# p_1_4 = 17.81
# p_0_2 = -0.05746 + p_0_1*1.25768 + p_1_4*0.63609
#
# percent_recent_nus<-p_0_2
# percent_homeless_cases<-4.07
# percent_us_homeless_cases<-8.67
#
# #fb recent entry cases
# # percent_nus_recent_cases<-
# ##### ACS Table S0501
# total_population<-328239523/1e6
# percent_nus_pop<-44932901/328239523
#
# #need to also update the following CSV in inst/extdata/
# #tb.deaths.us.2015.2018.csv CDC Wonder Multiple Cause of Death
# #ag.nat.pop.us.2018.csv ACS table
# #ag.deaths.us.2017.csv mortality.org
#
# ##### load in the current calibration data
  # model_load()
# newCalibDat<-CalibDat
#   #replace with new values
#
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### TB CASES   ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   newCalibDat$tot_cases<-rbind(newCalibDat$tot_cases, c(year, total_cases/1e6))
#
#   newCalibDat$fb_cases<-rbind(newCalibDat$fb_cases, c(year,percent_nus_cases/100, total_cases))
#   newCalibDat$fb_cases[25,1] <- 2017
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### TB CASES BY AGE   ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   newCalibDat$age_cases<-rbind(newCalibDat$age_cases,
#                                c(year,
#                                  age_cases,
#                                  total_cases))
#
# newCalibDat$age_cases_fb<-rbind(newCalibDat$age_cases_fb,
#                                 c(year,
#                                 age_cases_nus/100,
#                                 total_cases*(percent_nus_cases/100)))
# newCalibDat$age_cases_us<-rbind(newCalibDat$age_cases_us,
#                                 c(year,
#                                 age_cases_us/100,
#                                 total_cases*(1-(percent_nus_cases/100))))
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### NUS RECENT CASES    ##### ##### ##### ##### ##### ##### ##### ##### ###
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   newCalibDat$fb_recent_cases2 <- rbind(newCalibDat$fb_recent_cases2,
#                                       c(year,
#                                         percent_recent_nus/100,
#                                         total_cases*(percent_nus_cases/100)))
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### HOMELESS CASES    ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#
#   newCalibDat$homeless_cases<-rbind(newCalibDat$homeless_cases,c(year,percent_homeless_cases/100,total_cases))
#
#   newCalibDat$us_homeless_cases<-rbind(newCalibDat$us_homeless_cases,
#                                      c(year,percent_us_homeless_cases/100,
#                                        us_cases
#                                        ))

#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### TB DEATHS   ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#
#   newCalibDat$tb_deaths<-newCalibDat$tb_deaths[1:16,]
#   tbdeaths<-read.csv(system.file("extdata", "tb.deaths.us.2015.2018.csv", package="MITUS"))
#   colnames(tbdeaths)<-colnames(newCalibDat$tb_deaths)
#   rownames(tbdeaths)<-c("2015","2016","2017", "2018")
#   newCalibDat$tb_deaths<-rbind(newCalibDat$tb_deaths,
#                                tbdeaths)
#
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### TOTAL POPULATION YR FB  ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   newCalibDat$tot_pop_yr_fb<-rbind(newCalibDat$tot_pop_yr_fb,
#                                    c(year, total_population,
#                                      total_population*(1-(percent_nus_pop/100)),
#                                      total_population*(percent_nus_pop/100)))
#   #replace the population age structure
#   newCalibDat$tot_pop19_ag_fb<-newCalibDat$tot_pop18_ag_fb
#   newCalibDat$tot_pop19_ag_fb[,2:4]<-cbind(   c(.061,.262,.264,.132,.128,.088,.045,.019)*newCalibDat$tot_pop_yr_fb[9,1],
#                                               c(.07,.284,.245,.123,.126,.088,.045,.019)*newCalibDat$tot_pop_yr_fb[9,2],
#                                               c(.007,.121,.386,.192,.142,.088, .047,.018)*newCalibDat$tot_pop_yr_fb[9,3]
# )
#
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### TOTAL MORTALITY COUNTS  ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   mortality<-rep(0,12)
#   mortality<-c(12649,5450,30156,58845,80386,164850,374856,543800, 675243, 705213,175102)
#   mortality_yr<-2018
#   newCalibDat$US_tot_mort<-rbind(newCalibDat$US_tot_mort,
#                                  c(mortality_yr,sum(mortality)))
#   #age mortality is condensed
#   newCalibDat$US_mort_age<-rbind(newCalibDat$US_mort_age,
#                                  c(mortality_yr,
#                                    mortality[1],
#                                    sum(mortality[2:3]),
#                                    sum(mortality[4:5]),
#                                    mortality[6:9],
#                                    sum(mortality[10:11])))
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### IMPORTANCE WEIGHTS      ##### ##### ##### ##### ##### ##### ##### #####
#   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#   LgtCurveY <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
#   zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(year-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
#   zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }
#
#   ImptWeights <- LgtCurveY(1950,year,0.90)+0.10
#   names(ImptWeights) <- 1950:year
#   ImptWeights <- ImptWeights/max(ImptWeights)
#
#   newCalibDat$ImptWeights<-ImptWeights
# saveRDS(newCalibDat, file = paste0("~/MITUS/inst/US/US_CalibDat_", Sys.Date(), "_.rds"), version=2)
