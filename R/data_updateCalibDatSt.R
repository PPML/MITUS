#all data is to be read in via spreadsheet; some of the examples add multiple
#years so it is unlikely it can be used without revision


update_CalibDat_st<-function(){
  #load in the model (use CA as example)
  model_load("ND")
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
                                        c(2020,TBcases[i,3]))
  }
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB cases by age and nativity in 5 year bands
#found in CDC OTIS
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# age_state_cases<-read.csv("~/MITUS/inst/extdata/age_cases_states_07-22-20.csv",
#                           stringsAsFactors = FALSE)
age_state_cases<-read.csv("~/Desktop/TBcases_age_nat_st_5yr.csv")
age_state_cases<-age_state_cases[,-c(15,16)]
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
#TB cases by nativity in 5 year bands
#found in CDC OTIS
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# age_state_cases<-read.csv("~/MITUS/inst/extdata/age_cases_states_07-22-20.csv",
#                           stringsAsFactors = FALSE)
nat_state_cases<-read.csv("~/Desktop/TBcases_st_nat_5yrdata.csv")
newCalibDat[["cases_nat_st_5yr"]] <- nat_state_cases

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB cases homeless
#found in CDC OTIS
#needs to be smoothed
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

# cases_hr_95_99 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_95-99.txt")
# cases_hr_01_05 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_00-04.txt")
# cases_hr_05_09 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_05-09.txt")
# cases_hr_10_14 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_10-14.txt")
# cases_hr_15_19 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_15-19.txt")

cases_hr_96_00 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_96-00.txt")#, stringsAsFactors = FALSE)
cases_hr_01_05 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_01-05.txt")
cases_hr_06_10 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_06-10.txt")
cases_hr_06_10 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_11-15.txt")
cases_hr_16_20 <- read.delim("~/Downloads/OTIS TB Data 1993-2019_hr_16-20.txt")

cases_hr_0  <- data.frame(years=c("96_00","01_05","06_10","06_10","16_20"),pct_hr=NA,sample_size=NA)
cases_hr <- list()

datt <- as.data.frame(expand.grid(stateID[,1],1:5))
names(datt) <- c("st","yr")
datt$hr <- datt$ss <- NA

for(st in 1:51){ # st=5
  cases_hr[[st]] <- cases_hr_0
  id <- cases_hr_96_00$State==unique(cases_hr_96_00$State)[st]
  cases_hr[[st]]$sample_size[1] <- sum(cases_hr_96_00$Cases[id])
  cases_hr[[st]]$pct_hr[1]      <- sum(cases_hr_96_00$Cases[id & cases_hr_96_00$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_96_00$Cases[id & cases_hr_96_00$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==1] <- sum(cases_hr_96_00$Cases[id & cases_hr_96_00$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==1] <- sum(cases_hr_96_00$Cases[id & cases_hr_96_00$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_01_05$State==unique(cases_hr_01_05$State)[st]
  cases_hr[[st]]$sample_size[2] <- sum(cases_hr_01_05$Cases[id])
  cases_hr[[st]]$pct_hr[2]      <- sum(cases_hr_01_05$Cases[id & cases_hr_01_05$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_01_05$Cases[id & cases_hr_01_05$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==2] <- sum(cases_hr_01_05$Cases[id & cases_hr_01_05$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==2] <- sum(cases_hr_01_05$Cases[id & cases_hr_01_05$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_06_10$State==unique(cases_hr_06_10$State)[st]
  cases_hr[[st]]$sample_size[3] <- sum(cases_hr_06_10$Cases[id])
  cases_hr[[st]]$pct_hr[3]      <- sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==3] <- sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==3] <- sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_06_10$State==unique(cases_hr_06_10$State)[st]
  cases_hr[[st]]$sample_size[4] <- sum(cases_hr_06_10$Cases[id])
  cases_hr[[st]]$pct_hr[4]      <- sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==4] <- sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==4] <- sum(cases_hr_06_10$Cases[id & cases_hr_06_10$Homeless.in.Past.Year=="Yes"])

  id <- cases_hr_16_20$State==unique(cases_hr_16_20$State)[st]
  cases_hr[[st]]$sample_size[5] <- sum(cases_hr_16_20$Cases[id])
  cases_hr[[st]]$pct_hr[5]      <- sum(cases_hr_16_20$Cases[id & cases_hr_16_20$Homeless.in.Past.Year=="Yes"])/
    sum(cases_hr_16_20$Cases[id & cases_hr_16_20$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$ss[datt$st==stateID[,1][st] & datt$yr==5] <- sum(cases_hr_16_20$Cases[id & cases_hr_16_20$Homeless.in.Past.Year%in%c("No","Yes")])
  datt$hr[datt$st==stateID[,1][st] & datt$yr==5] <- sum(cases_hr_16_20$Cases[id & cases_hr_16_20$Homeless.in.Past.Year=="Yes"])
}

newCalibDat[["hr_cases"]] <- cases_hr

# smooth
library(mgcv)

datt$yr_sq <- datt$yr^2
ft_hr <- gam(hr~ yr + yr_sq + ti(st,bs="re") + s(yr,st,bs="re") + s(yr_sq,st,bs="re"),gamma=2,offset=log(ss),family=poisson(link = "log"),data=datt)

summary(ft_hr)
hist(exp(predict(ft_hr)) )
datt$fract    <- datt$hr/datt$ss
datt$fract_sm<-NA
datt$fract_sm[which(is.na(datt$fract)==FALSE)] <- exp(predict(ft_hr))
datt$yr_txt <- NA
datt$yr_txt[datt$yr==1] <- "96_00"
datt$yr_txt[datt$yr==2] <- "00_04"
datt$yr_txt[datt$yr==3] <- "06_10"
datt$yr_txt[datt$yr==4] <- "06_10"
datt$yr_txt[datt$yr==5] <- "16_20"

newCalibDat[["hr_cases_sm"]] <- datt

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
#TB cases among new migrants
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

cases_yse_96_00 <- read.delim("~/Downloads/OTIS TB Data 1993-2020_yse_96-00.txt")
cases_yse_01_05 <- read.delim("~/Downloads/OTIS TB Data 1993-2020_yse_01-05.txt")
cases_yse_06_10 <- read.delim("~/Downloads/OTIS TB Data 1993-2020_yse_06-10.txt")
cases_yse_11_15 <- read.delim("~/Downloads/OTIS TB Data 1993-2020_yse_11-15.txt")
cases_yse_16_20 <- read.delim("~/Downloads/OTIS TB Data 1993-2020_yse_16-20.txt")

#cases_yse_0  <- data.frame(years=c("96_00","01_05","06_10","11_15","16_20"),pct_0_1=NA,pct_1_4=NA,sample_size=NA)
#cases_yse <- list()

datt_0_1 <- as.data.frame(expand.grid(unique(StateID[,1]),1:5))
names(datt_0_1) <- c("st","yr")
datt_0_1$yse <- datt_0_1$ss <- NA
datt_1_4 <- datt_0_1

for(st in 1:51){ # st=1

  zz <- cases_yse_96_00
  id <- zz$State==unique(zz$State)[st]
  datt_0_1$ss[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==1] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_0_1$yse[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==1] <- sum(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"])==FALSE)])
  datt_1_4$ss[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==1] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_1_4$yse[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==1] <- sum(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"])==FALSE)])

  zz <- cases_yse_01_05
  id <- zz$State==unique(zz$State)[st]
  datt_0_1$ss[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==2] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_0_1$yse[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==2] <- sum(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"])==FALSE)])
  datt_1_4$ss[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==2] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_1_4$yse[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==2] <- sum(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"])==FALSE)])

  zz <- cases_yse_06_10
  id <- zz$State==unique(zz$State)[st]
  datt_0_1$ss[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==3] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_0_1$yse[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==3] <- sum(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"])==FALSE)])
  datt_1_4$ss[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==3] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_1_4$yse[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==3] <- sum(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"])==FALSE)])

  zz <- cases_yse_11_15
  id <- zz$State==unique(zz$State)[st]
  datt_0_1$ss[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==4] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_0_1$yse[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==4] <- sum(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"])==FALSE)])
  datt_1_4$ss[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==4] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_1_4$yse[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==4] <- sum(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"])==FALSE)])

  zz <- cases_yse_16_20
  id <- zz$State==unique(zz$State)[st]
  datt_0_1$ss[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==5] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_0_1$yse[datt_0_1$st==unique(StateID[,1])[st] & datt_0_1$yr==5] <- sum(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="< 1 year"])==FALSE)])
  datt_1_4$ss[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==5] <- sum(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")][which(is.na(zz$Cases[id & zz$Years.in.U.S.%in%c("< 1 year","1-4 years","5-14 years","15+ years")])==FALSE)])
  datt_1_4$yse[datt_1_4$st==unique(StateID[,1])[st] & datt_1_4$yr==5] <- sum(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"][which(is.na(zz$Cases[id & zz$Years.in.U.S.=="1-4 years"])==FALSE)])
}

# smooth
datt_0_1$yr_sq <- datt_1_4$yr_sq <- datt_1_4$yr^2

ft_yse01 <- gam(yse~ yr + yr_sq + s(st,bs="re") + s(yr,st,bs="re") + s(yr_sq,st,bs="re"),gamma=2,offset=log(ss),family=poisson(link = "log"),data=datt_0_1[datt_0_1$ss>0,])
summary(ft_yse01)
ft_yse14 <- gam(yse~ yr + yr_sq + s(st,bs="re") + s(yr,st,bs="re") + s(yr_sq,st,bs="re"),gamma=2,offset=log(ss),family=poisson(link = "log"),data=datt_1_4[datt_1_4$ss>0,])
summary(ft_yse14)

datt <- datt_1_4[,-4]
datt$yse01 <- datt_0_1[,4]
datt$yse14 <- datt_1_4[,4]
datt$p_0_1 <- NA
datt$p_0_1[which(is.na(datt$yse01)==FALSE)] <- exp(predict(ft_yse01,newdata=datt_0_1))
datt$p_1_4 <- NA
datt$p_1_4[which(is.na(datt$yse14)==FALSE)] <- exp(predict(ft_yse14,newdata=datt_1_4))
datt$p_recent <- -0.05746 + datt$p_0_1*1.25768 + datt$p_1_4*0.63609

hist(datt$p_recent)

datt$yr_txt <- NA
datt$yr_txt[datt$yr==1] <- "96_00"
datt$yr_txt[datt$yr==2] <- "01_05"
datt$yr_txt[datt$yr==3] <- "06_10"
datt$yr_txt[datt$yr==4] <- "11_15"
datt$yr_txt[datt$yr==5] <- "16_20"

newCalibDat[["rt_fb_cases_sm"]] <- datt
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
