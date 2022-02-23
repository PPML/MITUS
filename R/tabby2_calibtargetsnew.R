model_calib_targets <- function(){

#### united states ####
model_load("US")
loc <-"US"
#### total population ####

#### total tb cases ####
saveRDS(CalibDat$tot_cases*1e6, file=paste0("~/MITUS/inst/", loc, "/calibration_targets/",loc, "_cases_yr_02-10-22.rds"), version=2)

#### tb cases by nativity ####
saveRDS(CalibDat$age_cases, file=paste0("~/MITUS/inst/", loc, "/calibration_targets/",loc, "_age_cases_tot_02-10-22.rds"), version=2)

#### tb cases by age ####
saveRDS(CalibDat$fb_cases, file=paste0("~/MITUS/inst/", loc, "/calibration_targets/",loc, "_fb_cases_02-10-22.rds"), version=2)

#### tb cases by age and nativity ####
temp <- read_csv("~/Desktop/MITUS 2020 TB Data Update/us_age_nat_cases_5yr_2102022.csv")
saveRDS(temp, file="~/MITUS/inst/US/calibration_targets/US_ag_nat_cases_5yr_02-10-22.rds", version=2)
#### #### #### #### #### #### #### ####
#### states ####
model_load("CA")

#### total population ####

#### total tb cases ####
for (i in 1:51){
  loc<-stateID[i,3]
  saveRDS(CalibDat$cases_yr_st[[i]], file=paste0("~/MITUS/inst/", loc, "/calibration_targets/",loc, "_cases_yr_02-10-21.rds"), version=2)
}
#### 5 yr tb cases by nativity ####
for (i in 1:51){
  loc<-stateID[i,3]
  x <- CalibDat$cases_nat_st_5yr[((3*i)-2):(3*i),]
  saveRDS(x, file=paste0("~/MITUS/inst/", loc, "/calibration_targets/",loc, "_nat_cases_5yr_02-10-21.rds"), version=2)
}
#### 5 yr tb cases by age ####
for (i in 1:51){
  loc<-stateID[i,3]
  x<- (CalibDat[["cases_yr_ag_nat_st_5yr"]][[i]][1:5,-c(1,2,3,4)]+
         CalibDat[["cases_yr_ag_nat_st_5yr"]][[i]][6:10,-c(1,2,3,4)])
  y<-data.frame("Years"=CalibDat[["cases_yr_ag_nat_st_5yr"]][[i]][1:5,2],x/rowSums(x),"Total"=rowSums(x))
  saveRDS(y,file = paste0("~/MITUS/inst/",loc,"/calibration_targets/",loc,"_age_cases_tot_02-10-22.rds"), version=2)
}
#### ltbi in US ####

#### ltbi in NUS ####

#### tb deaths ####
}
