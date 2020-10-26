#read in the data from the spreadsheet
#this data is in 5 year bands for each state due to data restrictions
update_CalibDat_st<-function(){
age_state_cases<-read.csv("~/MITUS/inst/extdata/age_cases_states_07-22-20.csv",
                          stringsAsFactors = FALSE)
age_state_cases<-age_state_cases[,-15]
# age_state_cases[,5:14]<-age_state_cases[,5:14]/rowSums(age_state_cases[,5:14])
stateID<-as.matrix(read.csv(file="inst/extdata/state_ID.csv", header = TRUE))
cases_yr_ag_nat_st_5yr<-list()
for (i in 1:51){
  loc<-stateID[i,3]
  model_load(loc)
  NewCalibDat<-CalibDat
  library(dplyr)
  cases_yr_ag_nat_st_5yr[[i]]<-age_state_cases %>%  dplyr::filter(state==stateID[i,1])
}
NewCalibDat[["cases_yr_ag_nat_st_5yr"]]<- cases_yr_ag_nat_st_5yr

saveRDS(NewCalibDat,file="~/MITUS/inst/ST/ST_CalibDat_10-26-20.rds", version=2)


}
