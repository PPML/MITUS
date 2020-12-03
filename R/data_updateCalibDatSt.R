#all data is to be read in via spreadsheet; some of the examples add multiple
#years so it is unlikely it can be used without revision


update_CalibDat_st<-function(){
  #load in the model (use CA as example)
  model_load("CA")
  #assign old calibdat to a new version
  newCalibDat<-CalibDat
  ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
  #### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
  #Total TB cases
  #found in CDC OTIS
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
  TBcases<-read.csv("~/Desktop/TBCasesStateTotal.csv", header=TRUE)
  for (i in 1:51){
    newCalibDat$cases_yr_st[[i]]<-rbind(CalibDat$cases_yr_st[[i]],
                                        c(2019,TBcases[i,3]))
  }
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB cases by age and nativity in 5 year bands
#found in CDC OTIS
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# age_state_cases<-read.csv("~/MITUS/inst/extdata/age_cases_states_07-22-20.csv",
#                           stringsAsFactors = FALSE)
age_state_cases<-read.csv("~/Desktop/age_cases_states_12-2-20.csv")
age_state_cases<-age_state_cases[,-15]
# age_state_cases[,5:14]<-age_state_cases[,5:14]/rowSums(age_state_cases[,5:14])
stateID<-as.matrix(read.csv(file="inst/extdata/state_ID.csv", header = TRUE))
cases_yr_ag_nat_st_5yr<-list()
for (i in 1:51){
  loc<-stateID[i,3]
  library(dplyr)
  cases_yr_ag_nat_st_5yr[[i]]<-age_state_cases %>%  dplyr::filter(state==stateID[i,1])
}
newCalibDat[["cases_yr_ag_nat_st_5yr"]]<- cases_yr_ag_nat_st_5yr

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB cases homeless
#found in CDC OTIS
#needs to be smoothed
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

cases_hr_95_99 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_95-99.txt")
cases_hr_99_03 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_00-04.txt")
cases_hr_05_09 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_05-09.txt")
cases_hr_10_14 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_10-14.txt")
cases_hr_15_19 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_15-19.txt")

cases_hr_0  <- data.frame(years=c("95_99","99_03","05_09","10_14","15_19"),pct_hr=NA,sample_size=NA)
cases_hr <- list()

datt <- as.data.frame(expand.grid(stateID[,1],1:5))
names(datt) <- c("st","yr")
datt$hr <- datt$ss <- NA

for(st in 1:51){ # st=5
  cases_hr[[st]] <- cases_hr_0
  id <- cases_hr_95_99$State==unique(cases_hr_95_99$State)[st]
  cases_hr[[st]]$sample_size[1] <- sum(cases_hr_95_99$Cases[id])
  cases_hr[[st]]$pct_hr[1]      <- sum(cases_hr_95_99$Cases[id & cases_hr_95_99$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_95_99$Cases[id & cases_hr_95_99$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==1] <- sum(cases_hr_95_99$Cases[id & cases_hr_95_99$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==1] <- sum(cases_hr_95_99$Cases[id & cases_hr_95_99$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_99_03$State==unique(cases_hr_99_03$State)[st]
  cases_hr[[st]]$sample_size[2] <- sum(cases_hr_99_03$Cases[id])
  cases_hr[[st]]$pct_hr[2]      <- sum(cases_hr_99_03$Cases[id & cases_hr_99_03$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_99_03$Cases[id & cases_hr_99_03$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==2] <- sum(cases_hr_99_03$Cases[id & cases_hr_99_03$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==2] <- sum(cases_hr_99_03$Cases[id & cases_hr_99_03$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_05_09$State==unique(cases_hr_05_09$State)[st]
  cases_hr[[st]]$sample_size[3] <- sum(cases_hr_05_09$Cases[id])
  cases_hr[[st]]$pct_hr[3]      <- sum(cases_hr_05_09$Cases[id & cases_hr_05_09$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_05_09$Cases[id & cases_hr_05_09$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==3] <- sum(cases_hr_05_09$Cases[id & cases_hr_05_09$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==3] <- sum(cases_hr_05_09$Cases[id & cases_hr_05_09$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_10_14$State==unique(cases_hr_10_14$State)[st]
  cases_hr[[st]]$sample_size[4] <- sum(cases_hr_10_14$Cases[id])
  cases_hr[[st]]$pct_hr[4]      <- sum(cases_hr_10_14$Cases[id & cases_hr_10_14$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_10_14$Cases[id & cases_hr_10_14$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==4] <- sum(cases_hr_10_14$Cases[id & cases_hr_10_14$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==4] <- sum(cases_hr_10_14$Cases[id & cases_hr_10_14$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_15_19$State==unique(cases_hr_15_19$State)[st]
  cases_hr[[st]]$sample_size[5] <- sum(cases_hr_15_19$Cases[id])
  cases_hr[[st]]$pct_hr[5]      <- sum(cases_hr_15_19$Cases[id & cases_hr_15_19$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_15_19$Cases[id & cases_hr_15_19$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==5] <- sum(cases_hr_15_19$Cases[id & cases_hr_15_19$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==5] <- sum(cases_hr_15_19$Cases[id & cases_hr_15_19$Homeless.in.Past.Year=="Yes"])
}

newCalibDat[["hr_cases"]] <- cases_hr

# smooth
library(mgcv)

datt$yr_sq <- datt$yr^2
ft_hr <- gam(hr~ yr + yr_sq + ti(st,bs="re") + s(yr,st,bs="re") + s(yr_sq,st,bs="re"),gamma=2,offset=log(ss),family=poisson(link = "log"),data=datt)

summary(ft_hr)
hist( exp(predict(ft_hr)) )
datt$fract    <- datt$hr/datt$ss
datt$fract_sm <- exp(predict(ft_hr))
datt$yr_txt <- NA
datt$yr_txt[datt$yr==1] <- "95_99"
datt$yr_txt[datt$yr==2] <- "00_04"
datt$yr_txt[datt$yr==3] <- "05_09"
datt$yr_txt[datt$yr==4] <- "10_14"
datt$yr_txt[datt$yr==5] <- "15_19"

newCalibDat[["hr_cases_sm"]] <- datt

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB deaths
#multiple cause of death in CDC wonder
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
TBdeaths<-read.csv("~/Desktop/TBDeathsState1718.csv", header=TRUE)
for (i in 1:51){
  newCalibDat$tbdeaths[[i]] <-rbind(CalibDat$tbdeaths[[i]],
                                      c(CalibDat$tbdeaths[[i]][1,1],2017,TBdeaths[i,3]),
                                      c(CalibDat$tbdeaths[[i]][1,1],2018,TBdeaths[i+51,3] ))
}


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB deaths by age (use the national level)
#multiple cause of death in CDC wonder
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
model_load("US")
NatCalibDat<-CalibDat
newCalibDat$tbdeaths_age_yr<-NatCalibDat$tb_deaths

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#update the weights
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

newCalibDat$ImptWeights<-NatCalibDat$ImptWeights


saveRDS(newCalibDat,file=paste0("~/MITUS/inst/ST/ST_CalibDat_", Sys.Date(), ".rds"), version=2)



}
