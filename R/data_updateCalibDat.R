# ####yearly update to the CalibDat
# year<-2018
# ##### OTIS DATA
# total_cases<-9025
# us_cases<-2666
# percent_nus_cases<-70.19
# percent_recent_nus_cases
# age_cases_us<- c(3.16,2.26,4.84,5.49,5.4,6.67,9.47,6.05,4.45,2.0)
# age_cases_nus<-c(.27,1.03,9.69,18.18,15.53,13.91,14.93,12.83,9.88,3.71)
# age_cases<-((age_cases_us/100*(1-percent_nus_cases/100)*total_cases)+(age_cases_nus/100*percent_nus_cases/100*total_cases))/total_cases
#
# percent_homeless_cases<-4.11
# percent_us_homeless_cases<-8.44
# ##### ACS Table S0501
# total_population<-327167439/1e6
# percent_nus_pop<-13.67
#
# #need to also update the following CSV in inst/extdata/
# #tb.deaths.us.2015.2018.csv CDC Wonder Multiple Cause of Death
# #ag.nat.pop.us.2018.csv ACS table
# #ag.deaths.us.2017.csv mortality.org
#
# ##### load in the current calibration data
#   # model_load()
#   newCalibDat<-CalibDat
#   #### 2017 updates
#   newCalibDat$tot_cases<-rbind(newCalibDat$tot_cases, c(2017, 9088/1e6))
#   newCalibDat$fb_cases<-rbind(newCalibDat$fb_cases, c(2-17,70.3/100, 9088))
#
#
#   newCalibDat$age_cases<-rbind(newCalibDat$age_cases,
#                                c(2017,
#                                  (((c(.33,1.14,9.24,18.53,14.92,13.93,15.87,12.81,9.51,3.67)/100* 9088*.703)
#                                + (c(3.85,2.41,4.77,6.41,5.12,6.84,9.03,5.92,3.87,1.96)/100*9088*(1-.703)))/9088),
#                                 9088))
#
#   newCalibDat$age_cases_fb<-rbind(newCalibDat$age_cases_fb,
#                                   c(2017,
#                                     c(.33,1.14,9.24,18.53,14.92,13.93,15.87,12.81,9.51,3.67)/100,
#                                     9088*.703))
#   newCalibDat$age_cases_us<-rbind(newCalibDat$age_cases_us,
#                                  c(2017,
#                                  (c(3.85,2.41,4.77,6.41,5.12,6.84,9.03,5.92,3.87,1.96)/100),
#                                  9088*(1-.703)))
#   newCalibDat$homeless_cases<-rbind(newCalibDat$homeless_cases,
#                                     c(2017,
#                                       0.45, 9088))
#   newCalibDat$us_homeless_cases<-rbind(newCalibDat$us_homeless_cases,
#                                        c(2017, 0.09155455904, 2676))
#   ####start new year ypdates
#   newCalibDat$tot_cases<-rbind(newCalibDat$tot_cases, c(year, total_cases/1e6))
#
#   newCalibDat$fb_cases<-rbind(newCalibDat$fb_cases, c(year,percent_nus_cases/100, total_cases))
#
#   newCalibDat$age_cases<-rbind(newCalibDat$age_cases,
#                                c(year,
#                                  age_cases,
#                                  total_cases))
#
#   newCalibDat$age_cases_fb<-rbind(newCalibDat$age_cases_fb,
#                                   c(year,
#                                   age_cases_nus/100,
#                                   total_cases*(percent_nus_cases/100)))
#   newCalibDat$age_cases_us<-rbind(newCalibDat$age_cases_us,
#                                   c(year,
#                                   age_cases_us/100,
#                                   total_cases*(1-(percent_nus_cases/100))))
#
#   newCalibDat$homeless_cases<-rbind(newCalibDat$homeless_cases,c(year,percent_homeless_cases))
#
#   newCalibDat$us_homeless_cases<-rbind(newCalibDat$us_homeless_cases,
#                                        c(year,percent_us_homeless_cases/100,
#                                          us_cases
#                                          ))
#   tbdeaths<-read.csv(system.file("extdata", "tb.deaths.us.2015.2018.csv", package="MITUS"))[2:4,]
#   colnames(tbdeaths)<-colnames(newCalibDat$tb_deaths)
#   rownames(tbdeaths)<-c("2015","2016","2017")
#   newCalibDat$tb_deaths<-rbind(newCalibDat$tb_deaths,
#                                tbdeaths)
#   newCalibDat$tot_pop_yr_fb<-rbind(newCalibDat$tot_pop_yr_fb,
#                                    c(year, total_population,
#                                      total_population*(1-(percent_nus_pop/100)),
#                                      total_population*(percent_nus_pop/100)))
#   #replace the population age structure
#   newCalibDat$tot_pop18_ag_fb<-newCalibDat$tot_pop16_ag_fb
#   newCalibDat$tot_pop18_ag_fb[,2:4]<-read.csv(file=system.file("extdata", "ag.nat.pop.us.2018.csv", package = "MITUS"), header <- FALSE)[1:9,1:3]
#
#   #update total mortality
#   mortality<-read.csv(file=system.file("extdata", "ag.deaths.us.2017.csv", package = "MITUS"), header <- FALSE)
# #age mortality is condensed
#   newCalibDat$US_tot_mort<-rbind(newCalibDat$US_tot_mort,
#                                  c(2017,mortality[12,2]))
#   #age mortality is condensed
#   newCalibDat$US_mort_age<-rbind(newCalibDat$US_mort_age,
#                                  c(2017,
#                                    mortality[1,2],
#                                    sum(mortality[2:3,2]),
#                                    sum(mortality[4:5,2]),
#                                    mortality[6:9,2],
#                                    sum(mortality[10:11,2])))
# ##finally update the importance weights
#   LgtCurveY <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
#   zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2018-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
#   zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }
#   ImptWeights <- LgtCurveY(1990,year,0.90)+0.10
#   names(ImptWeights) <- 1950:year
#   ImptWeights <- ImptWeights/max(ImptWeights)
#
#   newCalibDat$ImptWeights<-ImptWeights
#
# saveRDS(newCalibDat, file = "~/MITUS/inst/US/US_CalibDat_04-13-20.rds", version=2)
