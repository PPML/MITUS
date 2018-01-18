setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC large grant/Analysis Transmission")

###### BRING IN WHO TB DATA

# Load Data
Burden <- read.csv("WHO tables/TB_burden_countries_2015-07-31.csv")
Notifications <- read.csv("WHO tables/TB_notifications_2015-07-31.csv")
MdrBurden <- read.csv("WHO tables/MDR-TB_burden_estimates_2015-07-31.csv")

# Merge
DatAll <- merge(Burden,MdrBurden,by=c("iso2","iso3","iso_numeric","year"),all.x=T,all.y=T)
DatAll <- merge(DatAll,Notifications,by=c("iso2","iso3","iso_numeric","year"),all.x=T,all.y=T)
Remove <- which(colnames(DatAll)=="g_whoregion.y");DatAll <- DatAll[,-Remove]
Remove <- which(colnames(DatAll)=="country.y"); DatAll <- DatAll[,-Remove]
colnames(DatAll)
# Some data crunching

DatAll$country <- factor(DatAll$country); DatAll$iso2 <- factor(DatAll$iso2); 
DatAll$iso3 <- factor(DatAll$iso3); DatAll$iso_numeric <- factor(DatAll$iso_numeric)
DatAll$g_whoregion.x <- factor(DatAll$g_whoregion.x);

DatAll <- DatAll[DatAll[,"year"]>=1990,]

#### Tables by country and year   
 TotPrevNum <-  matrix(0,length(unique(DatAll$iso3)),length(1950:2100)); rownames(TotPrevNum) <- unique(DatAll$iso3);colnames(TotPrevNum) <- 1950:2100
TotPrev <- TotPop <- TotPrevNum
TotPrevMdrE <- rep(NA,length(unique(DatAll$iso3)));names(TotPrevMdrE) <- unique(DatAll$iso3)
TotPrevMdrN <- TotPrevMdrE
for(i in 1:length(unique(DatAll$iso3))) {
  TotPrevNum[i,colnames(TotPrevNum)%in%DatAll[DatAll$iso3==unique(DatAll$iso3)[i],"year"]] <- DatAll[DatAll$iso3==unique(DatAll$iso3)[i],"e_prev_num"]
  TotPop[i,colnames(TotPop)%in%DatAll[DatAll$iso3==unique(DatAll$iso3)[i],"year"]]         <- DatAll[DatAll$iso3==unique(DatAll$iso3)[i],"e_pop_num"]
  TotPrev[i,colnames(TotPrev)%in%DatAll[DatAll$iso3==unique(DatAll$iso3)[i],"year"]]       <- DatAll[DatAll$iso3==unique(DatAll$iso3)[i],"e_prev_100k"]
  if(length(DatAll[DatAll$iso3==unique(DatAll$iso3)[i] & DatAll$year==2013,"e_new_mdr_pcnt"])>0) {
    TotPrevMdrN[i]       <- DatAll[DatAll$iso3==unique(DatAll$iso3)[i] & DatAll$year==2013,"e_new_mdr_pcnt"] 
    TotPrevMdrE[i]       <- DatAll[DatAll$iso3==unique(DatAll$iso3)[i] & DatAll$year==2013,"e_ret_mdr_pcnt"] } }

for(i in 1:nrow(TotPrev)) TotPrev[i,TotPrev[i,]==0] <- min(TotPrev[i,TotPrev[i,]>0])

#### Regional prevalence 
RegPrevNum <- matrix(NA,length(unique(DatAll$g_whoregion.x))+1,length(1950:2100)); 
rownames(RegPrevNum) <- c(as.character(unique(DatAll$g_whoregion.x)),"ALL"); colnames(RegPrevNum) <- 1950:2100
RegPrev <- RegPop <- RegPrevNum
RegPrevMdrE <- rep(NA,length(unique(DatAll$g_whoregion.x))); names(RegPrevMdrE) <- unique(DatAll$g_whoregion.x)
RegPrevMdrN <- RegPrevMdrE
i=1
for(i in 1:length(unique(DatAll$g_whoregion.x))) {
  ctyID <- as.character(unique(DatAll$iso3[DatAll$g_whoregion.x%in%unique(DatAll$g_whoregion.x)[i]]))
  RegPop[i,] <- colSums(TotPop[ctyID,]); RegPrevNum[i,] <- colSums(TotPrevNum[ctyID,]);
  RegPrev[i,] <- RegPrevNum[i,]/RegPop[i,]*100000
  RegPrevMdrN[i] <- sum(TotPrevNum[ctyID,"2013"]*TotPrevMdrN[ctyID],na.rm=T)/sum(TotPrevNum[ctyID,"2013"][is.na(TotPrevMdrN[ctyID])==F],na.rm=T)
  RegPrevMdrE[i] <- sum(TotPrevNum[ctyID,"2013"]*TotPrevMdrE[ctyID],na.rm=T)/sum(TotPrevNum[ctyID,"2013"][is.na(TotPrevMdrE[ctyID])==F],na.rm=T)
}
RegPop["ALL",] <- colSums(TotPop); RegPrevNum["ALL",] <- colSums(TotPrevNum);
RegPrev["ALL",] <- RegPrevNum["ALL",]/RegPop["ALL",]*100000
RegPrevMdrN["ALL"] <- sum(TotPrevNum[,"2013"]*TotPrevMdrN,na.rm=T)/sum(TotPrevNum[,"2013"][is.na(TotPrevMdrN[ctyID])==F],na.rm=T)
RegPrevMdrE["ALL"] <- sum(TotPrevNum[,"2013"]*TotPrevMdrE,na.rm=T)/sum(TotPrevNum[,"2013"][is.na(TotPrevMdrE[ctyID])==F],na.rm=T)


## PLOT
library(RColorBrewer)
colZ=brewer.pal(9,"Set1")[-6]

plot(0,0,col=NA,xlim=c(1990,2013),ylim=c(0,500),las=1,ylab="Prevalence per 100K",xlab="Year")
for(i in 1:nrow(RegPrev)) lines(1990:2013,RegPrev[i,41:64],col=colZ[i],lwd=2)
lines(1990:2013,RegPrev["ALL",41:64],col=1,lwd=4)
legend("topright",rownames(RegPrev),lwd=3,col=c(colZ[1:6],1))
mtext("Prevalence 1990-2013 by major world region",3,1.5,font=2,cex=1.3)

#### Project prevalence backwards
PriorYrFactor <- (RegPrev["ALL",41]/RegPrev["ALL",51])^0.1
(1-1/PriorYrFactor)*100

# Prior years
for(i in 1:nrow(RegPrev)) RegPrev[i,1:40] <- RegPrev[i,41]*cumprod(c(seq(RegPrev[i,41]/RegPrev[i,42],PriorYrFactor,length.out=6)[-1],rep(PriorYrFactor,35)))[40:1]
for(i in 1:nrow(TotPrev)) TotPrev[i,1:40] <- TotPrev[i,41]*cumprod(c(seq(TotPrev[i,41]/TotPrev[i,42],PriorYrFactor,length.out=6)[-1],rep(PriorYrFactor,35)))[40:1]
dim(RegPrev)
for(i in 1:nrow(RegPrev)) RegPrev[i,65:151] <- RegPrev[i,64]*cumprod(c(seq(RegPrev[i,64]/RegPrev[i,63],0.99,length.out=6)[-1],rep(0.99,82)))
for(i in 1:nrow(TotPrev)) TotPrev[i,65:151] <- TotPrev[i,64]*cumprod(c(seq(TotPrev[i,64]/TotPrev[i,63],0.99,length.out=6)[-1],rep(0.99,82)))

#
  plot(0,0,col=NA,xlim=c(1950,2100),ylim=c(0,650),las=1,ylab="Prevalence per 100K",xlab="Year")
for(i in 1:nrow(RegPrev)) lines(1950:2100,RegPrev[i,],col=colZ[i],lwd=2,lty=3)
for(i in 1:nrow(RegPrev)) lines(1990:2013,RegPrev[i,41:64],col=colZ[i],lwd=2)
lines(1950:2100,RegPrev["ALL",],col=1,lwd=4,lty=3)
lines(1990:2013,RegPrev["ALL",41:64],col=1,lwd=4)
  abline(v=1990,col="grey80");text(1990,600,"Estimates 1950-1989 obtained by projecting 1990-2000 trends
                                 in WHO global estimate backwards",cex=0.8,pos=)
  legend("topright",rownames(RegPrev),lwd=3,col=c(colZ[1:6],1))
  mtext("Prevalence 1990-2013 by major world region",3,1.5,font=2,cex=1.3)

########################################
########################################
## Bring in immigration data
  ImmigDat <- read.csv("LPR_by_year_country3.csv")
  colnames(ImmigDat) <- c("iso3","cty",1945+0:5*10,2000:2013,2000+3:20*5)
  ImmigDat2 <- as.data.frame(matrix(NA,nrow(ImmigDat),length(1945:2100)))
  colnames(ImmigDat2) <- 1945:2100; rownames(ImmigDat2) <- ImmigDat$iso3
  for(i in 1:nrow(ImmigDat))   ImmigDat2[i,] <- predict(smooth.spline(c(1945+0:5*10,2000:2013,2000+3:20*5),ImmigDat[i,-(1:2)],spar=0.5),1945:2100)$y 
  ImmigDat2 <- ImmigDat2[,-(1:5)]
  ImmigDat2[ImmigDat2<0] <- 0
  colnames(ImmigDat2)
  rnk <- order(-rowSums((ImmigDat2)[,1:64]))

  rowSums((ImmigDat2/colSums(ImmigDat2))[,1:64])[rnk[c(1,3:13,15:21,24)]]


plot(1,1,xlim=c(1950,2100),ylim=c(0.001,10),log="y")
for(i in 1:20) lines(1950:2100,ImmigDat2[rnk[c(1,3:13,15:21,24)][i],]/1e6)

ImmigDat2prev <- ImmigDat2; ImmigDat2prev[,] <- NA
  ImmigDat2prevMDRe <- rep(NA,nrow(ImmigDat2prev)); names(ImmigDat2prevMDRe) <- rownames(ImmigDat2prev)
  ImmigDat2prevMDRn <- ImmigDat2prevMDRe
  for(i in 1:nrow(ImmigDat2))  { if(rownames(ImmigDat2)[i]%in%rownames(TotPrev)) {
    ImmigDat2prev[i,] <- TotPrev[rownames(ImmigDat2)[i],]
    ImmigDat2prevMDRn[i] <- TotPrevMdrN[rownames(ImmigDat2)[i]]
    ImmigDat2prevMDRe[i] <- TotPrevMdrE[rownames(ImmigDat2)[i]]  }} 
  ImmigDat2prev["TWN",]                  <- TotPrev["CHN",]
  ImmigDat2prev["YUG",]                  <- RegPrev["EUR",]
  ImmigDat2prev["Other_Africa",]         <- RegPrev["AFR",]
  ImmigDat2prev["Other_America ",]       <- RegPrev["AMR",]
  ImmigDat2prev["Other_Asia",]           <- RegPrev["SEA",]
  ImmigDat2prev["Other_Caribbean",]      <- RegPrev["AMR",]
  ImmigDat2prev["Other_Europe",]         <- RegPrev["EUR",]
  ImmigDat2prev["Other_Oceania",]        <- RegPrev["WPR",]
  ImmigDat2prev["Other_South_America ",] <- RegPrev["AMR",]
  ImmigDat2prev["Not_Specified",]        <- RegPrev["ALL",]

ImmigDat2prevMDRn["TWN"]                  <- TotPrevMdrN["CHN"]
ImmigDat2prevMDRn["YUG"]                  <- RegPrevMdrN["EUR"]
ImmigDat2prevMDRn["Other_Africa"]         <- RegPrevMdrN["AFR"]
ImmigDat2prevMDRn["Other_America "]       <- RegPrevMdrN["AMR"]
ImmigDat2prevMDRn["Other_Asia"]           <- RegPrevMdrN["SEA"]
ImmigDat2prevMDRn["Other_Caribbean"]      <- RegPrevMdrN["AMR"]
ImmigDat2prevMDRn["Other_Europe"]         <- RegPrevMdrN["EUR"]
ImmigDat2prevMDRn["Other_Oceania"]        <- RegPrevMdrN["WPR"]
ImmigDat2prevMDRn["Other_South_America "] <- RegPrevMdrN["AMR"]
ImmigDat2prevMDRn["Not_Specified"]        <- RegPrevMdrN["ALL"]

ImmigDat2prevMDRe["TWN"]                  <- TotPrevMdrE["CHN"]
ImmigDat2prevMDRe["YUG"]                  <- RegPrevMdrE["EUR"]
ImmigDat2prevMDRe["Other_Africa"]         <- RegPrevMdrE["AFR"]
ImmigDat2prevMDRe["Other_America "]       <- RegPrevMdrE["AMR"]
ImmigDat2prevMDRe["Other_Asia"]           <- RegPrevMdrE["SEA"]
ImmigDat2prevMDRe["Other_Caribbean"]      <- RegPrevMdrE["AMR"]
ImmigDat2prevMDRe["Other_Europe"]         <- RegPrevMdrE["EUR"]
ImmigDat2prevMDRe["Other_Oceania"]        <- RegPrevMdrE["WPR"]
ImmigDat2prevMDRe["Other_South_America "] <- RegPrevMdrE["AMR"]
ImmigDat2prevMDRe["Not_Specified"]        <- RegPrevMdrE["ALL"]

PrevalenceInImmigrants <- colSums(ImmigDat2*ImmigDat2prev)/colSums(ImmigDat2)
TotalImmigrants <- colSums(ImmigDat2)
EntrantsWithActiveTB <- PrevalenceInImmigrants*TotalImmigrants
save(PrevalenceInImmigrants,file="PrevalenceInImmigrants_8-12-16.rData")
################
LgtCurve <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
 zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2100-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
 zz  <- as.numeric(EndVal)/(1+exp(-zz));  if(StYr>1950) { zz[1:((StYr-1950))] <- 0 };    zz  }
plot(1950:2013,LgtCurve(1970,2010,1)[1:64],type="l")
##################

PrevalenceInImmigrantsMDRe <- colSums(ImmigDat2*ImmigDat2prev*ImmigDat2prevMDRe)/colSums(ImmigDat2*ImmigDat2prev)*LgtCurve(1970,2010,1)[1:64]
PrevalenceInImmigrantsMDRn <- colSums(ImmigDat2*ImmigDat2prev*ImmigDat2prevMDRn)/colSums(ImmigDat2*ImmigDat2prev)*LgtCurve(1970,2010,1)[1:64]
plot(1950:2013,PrevalenceInImmigrantsMDRe,type="l",ylim=c(0,20))
plot(1950:2013,PrevalenceInImmigrantsMDRn,type="l",ylim=c(0,4))



TotalImmigrantsEur <- colSums(ImmigDat2[c("AUS","AUT","BEL","CAN","DNK",
                                          "FIN","FRA","DEU","GRC","IRL","ITA",
                                          "NLD","NZL","NOR","Other_Oceania","Other_Europe",
                                          "PRT","ESP","SWE","CHE","GBR"),])


crudeLPR <- read.csv("TotalLPR_by_year.csv")
###############
PriorYrFactorZ <- c((1/.98)^(40:1)/PriorYrFactor,rep(1,111))



pdfnam <- "Fig_Immigrants.pdf"

pdf(file=pdfnam,width=7, height=8) 
par(mfrow=c(3,1),mar=c(4.5,4,3,4))
  plot(1950:2100,TotalImmigrants/10^6,type="l",lwd=2,las=1,col=4,
       ylab="Total immigrants per year (mil)",xlab="Year",ylim=c(0,2.3))
mtext("Total immigrants per year",3,0.8,font=2,cex=0.85)
mtext("Pct of Total",4,2.5,cex=0.7)

lines(1950:2100,TotalImmigrantsEur/TotalImmigrants*1.6,lty=2,col=4)
axis(4,(1:5)/5*1.6,(1:5)/5*100,las=1)
legend(2040,1.2,c("Total immigrants","% from Western Europe/Oceania"),
       col=4,lwd=c(2,1),lty=c(1,2),cex=0.9)
 points(crudeLPR[,1],as.numeric(crudeLPR[,2])/10^6)
###
plot(1950:2100,PrevalenceInImmigrants,type="l",lwd=2,las=1,col=2,
     ylab="Prevalence per 100K",xlab="Year",ylim=c(0,280))
lines(1950:2100,PrevalenceInImmigrants*PriorYrFactorZ,lwd=2,col=2,lty=2)
mtext("Average prevalence in sender countries",3,0.8,font=2,cex=0.85)
legend("topright",c("Projected backwards from 1990 with 0.55% incr per year",
                       "Projected backwards from 1990 with 1.5% incr per year"),
       col=2,lwd=2,lty=c(1,2),cex=0.9)
##
plot(1950:2100,EntrantsWithActiveTB/max(EntrantsWithActiveTB),type="l",lwd=2,las=1,
     ylab="",xlab="Year",ylim=c(0,1.1))
lines(1950:2100,EntrantsWithActiveTB*PriorYrFactorZ/max(EntrantsWithActiveTB),lwd=2,lty=2)
mtext("Total volume of immigrant TB cases, scaled to max = 1.0",3,0.8,font=2,cex=0.85)
legend("bottomright",c("Prev. projected backwards from 1990 with 0.55% incr per yr",
                       "Prev. projected backwards from 1990 with 1.5% incr per yr"),
       lwd=2,lty=c(1,2),cex=0.9)
##
dev.off(); system(paste("open", pdfnam)) # code for ma

###########

################################################
################################################
#### PREVALENCE BY AGE
age_dist <- read.csv("Age_dist_by_country.csv")

colnames(DatAll)

DatAll[,"newrel_04"]   <- DatAll[,"newrel_m04"]+DatAll[,"newrel_f04"]
DatAll[,"newrel_514"]  <- DatAll[,"newrel_m514"]+DatAll[,"newrel_f514"]
DatAll[,"newrel_1524"] <- DatAll[,"newrel_m1524"]+DatAll[,"newrel_f1524"]
DatAll[,"newrel_2534"] <- DatAll[,"newrel_m2534"]+DatAll[,"newrel_f2534"]
DatAll[,"newrel_3544"] <- DatAll[,"newrel_m3544"]+DatAll[,"newrel_f3544"]
DatAll[,"newrel_4554"] <- DatAll[,"newrel_m4554"]+DatAll[,"newrel_f4554"]
DatAll[,"newrel_5564"] <- DatAll[,"newrel_m5564"]+DatAll[,"newrel_f5564"]
DatAll[,"newrel_65p"]  <- DatAll[,"newrel_m65"]+DatAll[,"newrel_f65"]

#### Tables by country and age group, 2013    
TotNotifAge <-  matrix(0,length(unique(DatAll$iso3)),8); 
rownames(TotNotifAge) <- unique(DatAll$iso3);
colnames(TotNotifAge) <- c("newrel_04","newrel_514","newrel_1524","newrel_2534","newrel_3544","newrel_4554","newrel_5564","newrel_65p")
NotifFracAge <- TotNotifAge

for(i in 1:length(unique(DatAll$iso3))) { for(j in 1:8) { # i = j = 1
  zz <- DatAll[DatAll$iso3==unique(DatAll$iso3)[i] & DatAll$year==2013,colnames(TotNotifAge)[j]]
  TotNotifAge[i,j] <- if(length(zz)>0) { as.numeric(zz) } else { 0 }   } }

NotifFracAge <- t(apply(TotNotifAge,1,function(x) x/sum(x)))
for(i in 1:ncol(NotifFracAge)) NotifFracAge[is.na(NotifFracAge[,i]),i] <- mean(NotifFracAge[is.na(NotifFracAge[,i])==F,i])  ## filling in missing with global average

###
RelNotRateAge <- PopFractAge <- NotifFracAge

for(i in 1:length(unique(DatAll$iso3))) { #i = j = 2
  zz <- age_dist[as.character(age_dist$iso3)==unique(DatAll$iso3)[i],3:10]
  PopFractAge[i,] <- if(nrow(zz)>0) { as.numeric(zz) } else { NA }   } 

for(i in 1:ncol(PopFractAge)) PopFractAge[is.na(PopFractAge[,i]),i] <- age_dist[1,(3:10)[i]]  ## filling in missing with global average

RelNotRateAge <- TotNotifAge/PopFractAge
RelNotRateAge <- t(apply(RelNotRateAge,1,function(x) x/max(x)))

########
# REGIONAL VALUES
# sum_across_countries(total tb * fract_notif_age_grps)/sum_across_countries(total_pop*pop_fract_age_grp)*10^5


RegRelNotRateAge <- matrix(NA,length(unique(DatAll$g_whoregion.x))+1,8); 
rownames(RegRelNotRateAge) <- c(as.character(unique(DatAll$g_whoregion.x)),"ALL"); colnames(RegRelNotRateAge) <- colnames(TotNotifAge)

for(i in 1:length(unique(DatAll$g_whoregion.x))) { for(j in 1:8) { #i = j = 1 
  ctyID <- as.character(unique(DatAll$iso3[DatAll$g_whoregion.x%in%unique(DatAll$g_whoregion.x)[i]]))
  RegRelNotRateAge[i,j] <- sum(TotPrevNum[ctyID,"2013"]*NotifFracAge[ctyID,j],na.rm=T)/
    sum(TotPop[ctyID,"2013"]*PopFractAge[ctyID,j],na.rm=T)*10^5; } }

for(j in 1:8) { RegRelNotRateAge["ALL",j] <- sum(TotPrevNum[,"2013"]*NotifFracAge[,j],na.rm=T)/
  sum(TotPop[,"2013"]*PopFractAge[,j],na.rm=T)*10^5; }

plot(RegRelNotRateAge["ALL",],type="l",lwd=2,las=1,ylim=c(0,max(RegRelNotRateAge["ALL",])))
for(i in 1:6) lines(RegRelNotRateAge[i,],col=i+1,lwd=2)

## PLOT

library(RColorBrewer)
colZ=brewer.pal(9,"Set1")[-6]

pdfnam <- "Fig_Immigrants_Age.pdf"

pdf(file=pdfnam,width=8.5, height=13) 
par(mfrow=c(3,1),mar=c(5.5,4,3,4))

plot(0,0,col=NA,xlim=c(0.9,8.1),ylim=c(0,max(RegRelNotRateAge)),las=1,ylab="Prevalence per 100K",xlab="Age group",axes=F)
axis(2,las=1);axis(1,1:8,c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65+")); box()
for(i in 1:6) lines(RegRelNotRateAge[i,],col=colZ[i],lwd=2)
lines(RegRelNotRateAge[7,],col=1,lwd=4)

legend("topleft",rownames(RegPrev),lwd=3,col=c(colZ[1:6],1))
mtext("Prevalence 2013 by age group and major world region",3,.7,font=2,cex=1.0)

## PLOT
plot(0,0,col=NA,xlim=c(0.9,8.1),ylim=c(0,1),las=1,ylab="Relative prevalence",xlab="Age group",axes=F)
axis(2,las=1);axis(1,1:8,c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65+")); box()
for(i in 1:6) lines(RegRelNotRateAge[i,]/max(RegRelNotRateAge[i,]),col=colZ[i],lwd=2)
lines(RegRelNotRateAge[7,]/max(RegRelNotRateAge[7,]),col=1,lwd=4)

legend("topleft",rownames(RegPrev),lwd=3,col=c(colZ[1:6],1))
mtext("Prevalence 2013 by age group and major world region,
      relative to max(age group)",3,.7,font=2,cex=1.0)


## 
### Age distribution
Lpr_age <- read.csv("LPR_by_age_year.csv",row.names=1)
ImmigAgeDist <- rowSums(Lpr_age)/sum(Lpr_age)
ImmigAgeDist2 <- ImmigAgeDist[1]
ImmigAgeDist2[2] <- sum(ImmigAgeDist[2:3])
ImmigAgeDist2[3] <- sum(ImmigAgeDist[4:5])
ImmigAgeDist2[4] <- sum(ImmigAgeDist[6:7])
ImmigAgeDist2[5] <- sum(ImmigAgeDist[8:9])
ImmigAgeDist2[6] <- sum(ImmigAgeDist[10:11])
ImmigAgeDist2[7] <- sum(ImmigAgeDist[12:13])
ImmigAgeDist2[8] <- ImmigAgeDist[14]
ImmigAgeDist2[9:11] <- ImmigAgeDist[15]*c(16,4,1)/21
names(ImmigAgeDist2) <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="-"),"95+")

## Adjust older agegroups up
plot(ImmigAgeDist2)
ImmigAgeDist2[5:11] <- ImmigAgeDist2[5:11]*c(0.9,0.85,0.75,0.65,0.45,0.3,0.1)
ImmigAgeDist2 <- ImmigAgeDist2/sum(ImmigAgeDist2)
lines(ImmigAgeDist2)

Imig_TB_Num <-  ImmigAgeDist2*RegRelNotRateAge[7,c(1:8,8,8,8)]
Imig_TB_Num <-  Imig_TB_Num/sum(Imig_TB_Num)
Imig_TB_Cases <- c(20+198,398,2851,3685,1705,1058,836,607,319,63)
Imig_TB_Cases <- Imig_TB_Cases/sum(Imig_TB_Cases)
Imig_TB_Num[10] <- Imig_TB_Num[10]+Imig_TB_Num[11];Imig_TB_Num <- Imig_TB_Num[-11]
length(Imig_TB_Num)

## PLOT
plot(0,0,col=NA,xlim=c(0.9,10.1),ylim=c(0,.35),las=1,ylab="Fraction of total cases",xlab="Age group",axes=F)
axis(2,las=1);axis(1,1:10,c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+")); box()
lines(Imig_TB_Cases,col=4,lwd=2); lines(Imig_TB_Num,col=4,lwd=2,lty=2)

legend("topright",c("Reported TB cases in recent immigrants (arrival <5 yrs), 2008-2012",
       "Distribution of immigrants with active disease (model assumptions)"),lwd=2,col=4,lty=c(1,2))
mtext("Age distribution of TB cases in recent immigrants (empirical) vs. 
      fraction of all immigrants with active TB based on model assumptions",3,.7,font=2,cex=1.0)

dev.off(); system(paste("open", pdfnam)) # 

#####
# sum_across_countries(total tb * fract_notif_age_grps)/sum_across_countries(total_pop*pop_fract_age_grp)*10^5

## intermediate = regional values 
## final goal = rate ratio of TB by age group relative to highest age group

### LTBI Prevalence in immigrants (BENNET FOR DATA)
### Proxy with hyperbolic distribution  A(t) = 1 - (b/(b+t))^a
###############  1999  ###########################
## DATA IN
load("calib files/IgraDat_1-13-2016.rData") # IgraDat

LTBI_prev_FB_11_IGRA <- IgraDat[IgraDat$nativity=="fb",c(3,9:10)]


alimZ=60
xd <- optim(c(.5,.5),F2,method="L-BFGS-B",lower=0.0001,control=list(fnscale=-1),alim=alimZ,w=1)
p= xd$par
LtbiPars <- xd$par; names(LtbiPars) <- c("a","b")
t=0:90+0.5; 
b = 50 # 3.443236
a = 0.4# 0.07848647
InfByAge <- 1-(b/(b+t))^a;InfByAge
LtbiPars <- c(a,b)
InfByAge2 <- 1-exp(-t*0.0045);InfByAge2
#################################################

nhanes <- read.csv("calib files/Hill_NHANES_1-13-16.csv")
nhanes <- nhanes[nhanes$test=="igra",]
nhanes_fb <- nhanes[nhanes$nativity=="fb",]

plot((0:8)*10,c(0,nhanes_fb[,6]),type="l",xlim=c(0,90),ylim=c(0,40))
pts <- (0:8)/10
p <- c(0.01,0.005)
1-exp(-pts*p[1]-pts^2*p[2])

Func <- function(p) {
  pred <- 1-exp(-pts*p[1]-pts^2*p[2])
  sum((c(0,nhanes_fb[,6])/100-pred)^2) }

opt <- optim(c(0.01,0.005),Func,method="BFGS")
opt <- optim(c(0.01,0.005),Func,method="Nelder-Mead")
p <- opt$par;p
#p <- opt$par+c(0.0,-0.1)
#lines((0:8)*10,1-exp(-pts*p[1]-pts*p[2]^2),col=2)
p <- c(0.25,0.41)
plot((0:8)*10,c(0,nhanes_fb[,6])/100,type="l",xlim=c(0,90),ylim=c(0,.4))
lines((0:8)*10,1-exp(-pts*p[1]-pts^2*p[2]),col=2)

####
pdfnam <- "Fig_Immigrants_LTBI.pdf"

pdf(file=pdfnam,width=9, height=5) 
par(mfrow=c(1,1),mar=c(3.5,4,.5,.5))

plot(0:90+0.5,InfByAge2,type="l",ylim=c(0,.55),xlim=c(2,90),las=1,
     xlab="",ylab="")
for(i in 1:25) { lines(0:90+0.5,1-exp(-t*0.0045*seq(0.01,8,length.out=25)[i]),col="grey90") }
lines(0:90+0.5,InfByAge2,lwd=2)
mtext("Age",1,2.5);mtext("Fraction with LTBI",2,3)
ageNHANES <- c(10,20,30,40,50,60,70,85)
points(ageNHANES,LTBI_prev_FB_11_IGRA[,2]/rowSums(LTBI_prev_FB_11_IGRA[,2:3]),col=2,pch=16)
for(i in 1:14) {  lines(rep(ageNHANES[i],2),qbeta(c(1,39)/40,LTBI_prev_FB_11_IGRA[i,2],LTBI_prev_FB_11_IGRA[i,3]),col=2) }
legend("topleft",c("Survey estimates, NHANES 1999-2000","Fitted exponential: r=0.0045",
       "Adjustment in curve for time-changes in LTBI prevalence"),lwd=c(1,2,1),col=c(2,1,"grey75"),
       pch=c(16,NA,NA),cex=0.8,bg="white")

dev.off(); system(paste("open", pdfnam)) # 

###########
  # FB Data
FbMdr <- read.csv("Otis_FbMdrPrior2.csv")
colnames(FbMdr)
FbMdr2 <- matrix(NA,0,5); colnames(FbMdr2) <- colnames(FbMdr)[-6]

for(i in 1:nrow(FbMdr)) {  FbMdr2 <- rbind(FbMdr2,cbind(FbMdr[rep(i,FbMdr[i,6]),-(5:6)],ifelse(FbMdr[i,5]=="No",0,1)))   }
colnames(FbMdr2) <- colnames(FbMdr)[-6] 
FbMdr2$Age_groupNum <- as.numeric(FbMdr2$Age_group)
colnames(FbMdr2)
library(mgcv)
unique(FbMdr2$Age_group)
unique(FbMdr2$Years_in_US)
unique(FbMdr2$Prev_TB)
unique(FbMdr$MDR)


##########

chk <- read.csv("OTIS MDR Check.csv")
colnames(chk)
sum(chk$n[chk$mdr%in%unique(chk$mdr)[1:2]])/sum(chk$n[chk$mdr%in%unique(chk$mdr)[1:3]])

plot(sum(chk$n[chk$mdr%in%unique(chk$mdr)[1:2]])/sum(chk$n[chk$mdr%in%unique(chk$mdr)[1:3]])


#######

#ft <- gam(MDR~Age_group+s(Year)+Years_in_US,data=FbMdr2,family=binomial(link = "logit"))
#ft <- gam(MDR~Age_group+s(Year)+Years_in_US,data=FbMdr2)
ft <- gam(MDR~s(Year)+Years_in_US+Prev_TB,data=FbMdr2)
summary(ft)
plot(ft)
byYear <- predict(ft,data.frame(Year       =1993:2013,
                                Years_in_US=unique(FbMdr2$Years_in_US)[2],
                                Prev_TB    =unique(FbMdr2$Prev_TB)[1]))
byYrInUs <- predict(ft,data.frame(Year=2013,
                               Years_in_US=unique(FbMdr2$Years_in_US),
                                Prev_TB    =unique(FbMdr2$Prev_TB)[1]))

byPriorTb <- predict(ft,data.frame(Year=2013,
                                  Years_in_US=unique(FbMdr2$Years_in_US)[2],
                                  Prev_TB    =unique(FbMdr2$Prev_TB)))

plot(1993:2013,byYear,ylim=c(0,.10),type="l",las=1)
plot(unique(FbMdr2$Years_in_US),byYrInUs,ylim=c(0,.10),type="l",las=1)
plot(unique(FbMdr2$Prev_TB),byPriorTb,ylim=c(0,.10),type="l",las=1)
################

zz <- FbMdr2[FbMdr2$Prev_TB=="Yes" & FbMdr2$Years_in_US%in%unique(FbMdr2$Years_in_US)[1:2],"MDR"]
MdrE <- sum(zz)/length(zz);MdrE # 0.15 MDR TB in treatment experienced arrivals
zz <- FbMdr2[FbMdr2$Prev_TB=="No" & FbMdr2$Years_in_US%in%unique(FbMdr2$Years_in_US)[1:2],"MDR"]
MdrN <- sum(zz)/length(zz); MdrN # 0.022 MDR TB in treatment naive arrivals

pctXdr <- 0.094 # Zignol 2012
InhRifMdrN <- c(6.2-1.0,1.2-1.0,1.0)
InhRifMdrN <- InhRifMdrN/InhRifMdrN[3]*MdrN
InhRifMdrE <- c(19.6-9.3,12.0-9.3,9.3)
InhRifMdrE <- InhRifMdrE/InhRifMdrE[3]*MdrE
InhRifMdrN*100;InhRifMdrE*100
#################
################
LgtCurve <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
   zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2100-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
   zz  <- as.numeric(EndVal)/(1+exp(-zz));  if(StYr>1950) { zz[1:((StYr-1950))] <- 0 };    zz  }
plot(1950:2013,LgtCurve(1970,2010,1)[1:64],type="l")
##################
LgtCurve2 <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
 zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2100-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr)/12)
 zz  <- as.numeric(EndVal)/(1+exp(-zz));  if(StYr>1950) { zz[1:((StYr-1950)*12)] <- 0 };    zz  }
######

pDRtN <- outer(LgtCurve(1970,2010,1)[1:64],InhRifMdrN)
pDRtN <- cbind(1-rowSums(pDRtN),pDRtN)
pDRtN <- cbind(pDRtN[,1:3],pDRtN[,4]*(1-LgtCurve(1995,2010,pctXdr)[1:64]),pDRtN[,4]*LgtCurve(1995,2010,pctXdr)[1:64])

pDRtN2 <- outer(LgtCurve2(1970,2010,1),InhRifMdrN)
pDRtN2 <- cbind(1-rowSums(pDRtN2),pDRtN2)
pDRtN2 <- cbind(pDRtN2[,1:3],pDRtN2[,4]*(1-LgtCurve2(1995,2010,pctXdr)),pDRtN2[,4]*LgtCurve2(1995,2010,pctXdr))

pDRtE <- outer(LgtCurve(1970,2010,1)[1:64],InhRifMdrE)
pDRtE <- cbind(1-rowSums(pDRtE),pDRtE)
pDRtE <- cbind(pDRtE[,1:3],pDRtE[,4]*(1-LgtCurve(1995,2010,pctXdr)[1:64]),pDRtE[,4]*LgtCurve(1995,2010,pctXdr)[1:64])

pDRtE2 <- outer(LgtCurve2(1970,2010,1),InhRifMdrE)
pDRtE2 <- cbind(1-rowSums(pDRtE2),pDRtE2)
pDRtE2 <- cbind(pDRtE2[,1:3],pDRtE2[,4]*(1-LgtCurve2(1995,2010,pctXdr)),pDRtE2[,4]*LgtCurve2(1995,2010,pctXdr))

rownames(pDRtN) <- rownames(pDRtE) <- 1950:2013
colnames(pDRtN) <- colnames(pDRtE) <- colnames(pDRtN2) <- colnames(pDRtE2) <-c("PanSens","MonoINH","MonoRIF","MDR","XDR")
FbDrFrac <- list(pDRtN=pDRtN,pDRtE=pDRtE)

############################################
#########

# FB Data
Fb <- read.csv("Otis_FbPrior2.csv")
colnames(Fb)
head(Fb)

Fb2 <- Fb[Fb$Years_in_US%in%unique(Fb$Years_in_US)[c(1,3)],]
unique(Fb$Year)
FbAgeYear <- matrix(NA,length(unique(Fb2$Year)),length(unique(Fb2$Age_group))-1)
for(i in 1:21) {
  FbAgeYear[i,1] <- sum(Fb2$Count[Fb2$Year==unique(Fb2$Year)[i] & Fb2$Age_group%in%unique(Fb$Age_group)[1:2]])
  for(j in 1:10)  {
    FbAgeYear[i,j+1] <- sum(Fb2$Count[Fb2$Year==unique(Fb2$Year)[i] & Fb2$Age_group==unique(Fb$Age_group)[j+2]]) }}
  FbAgeYear <- FbAgeYear[,-11]
rownames(FbAgeYear) <- 1993:2013; colnames(FbAgeYear) <- c("0_4","5_14","15_24","25_34","35_44","45_54","55_64","65_74","75_84","85+")
    
FbAgeYearPrev <- array(NA,c(length(unique(Fb2$Year)),length(unique(Fb2$Age_group))-1,2))
for(k in 1:2) { for(i in 1:21) { 
  FbAgeYearPrev[i,1,k] <- sum(Fb2$Count[Fb2$Year==unique(Fb2$Year)[i] & Fb2$Age_group%in%unique(Fb$Age_group)[1:2] & Fb2$Prev_TB==unique(Fb$Prev_TB)[c(1,3)[k]]])
  for(j in 1:10)  {
    FbAgeYearPrev[i,j+1,k] <- sum(Fb2$Count[Fb2$Year==unique(Fb2$Year)[i] & Fb2$Age_group==unique(Fb$Age_group)[j+2] &  Fb2$Prev_TB%in%unique(Fb$Prev_TB)[c(1,3)[k]]]) } } }
FbAgeYearPrev <- FbAgeYearPrev[,-11,]
rownames(FbAgeYearPrev) <- 1993:2013; colnames(FbAgeYearPrev) <- c("0_4","5_14","15_24","25_34","35_44","45_54","55_64","65_74","75_84","85+")
dimnames(FbAgeYearPrev)[[3]] <- c("TxNaive","TxExperienced")


FbAgeYearPrevz <- FbAgeYearPrev[,-10,]
FbAgeYearPrevz[,9,] <- FbAgeYearPrev[,9,]+FbAgeYearPrev[,10,]
FractTxEbyAge <- colSums(FbAgeYearPrevz[,,2])/colSums(FbAgeYearPrevz[,,1]+FbAgeYearPrevz[,,2])
barplot(FractTxEbyAge)


plot(colSums(FbAgeYear/rowSums(FbAgeYear))/21,type="l",lwd=3,lty=3,axes=F)
   ff <- colorRampPalette(c(2,4))  
for(i in 1:21) lines((FbAgeYear/rowSums(FbAgeYear))[i,],col=ff(21)[i],lwd=2)
  
  NotifAgeDist <- colSums(FbAgeYear/rowSums(FbAgeYear))/21
ImmigAgeDist2z <- ImmigAgeDist2[-(10:11)]
ImmigAgeDist2z[9] <- sum(ImmigAgeDist2[9:11])
NotifAgeDistz <- NotifAgeDist[-10]
NotifAgeDistz[9] <- sum(NotifAgeDist[9:10])
RelActDisAge <- (NotifAgeDistz/ImmigAgeDist2z)/(NotifAgeDistz/ImmigAgeDist2z)[4]
barplot(RelActDisAge,ylim=c(0,3.5),las=1)



Panel: 
1 total notifications by year, with LPR overplotted
2a Age distribution of immigrants
2b Average prevalence in sender countries
3a Age distribution of latent infection
3b Age distribution of active disease
4 Implied time trend in active/latent TB cases entering the US

################################
### Panel plot of all inputs

LprCrude <- read.csv("Lpr_by_Year.csv")
library(RColorBrewer)
colZ=brewer.pal(9,"Set1")[-6]

pdfnam <- "Fig_ImmigInputs.pdf"

pdf(file=pdfnam,width=9, height=12)
layout(rbind(c(1,1),2:3,4:5))

par(mar=c(4,3,3,.5))

plot(0,0,type="l",lwd=2,las=1,ylab="",xlab="",ylim=c(0,1.82),xlim=c(1951.5,2011.5))
abline(h=0:10*0.5,col="grey85")
for(i in 1:64) lines(rep((1950:2013)[i],2),c(0,TotalImmigrants[i]/1e6),col=colZ[2],lwd=5,lend="butt")
lines(1950:2013,LprCrude[,2]/1e6,col=colZ[1],lwd=2)
mtext("Year",1,2.5,cex=0.75)
mtext("A. Total immigrants per year (millions)",3,0.8,font=2,cex=0.9)
legend("topleft",c("Raw count of LPR per year","Total annual immigrants used for model"),col=colZ[1:2],,lwd=c(3,NA),pch=c(NA,15),pt.cex=1.5)
###
plot(0,0,xlim=c(0.8,11.2),ylim=c(1,28),ylab="",xlab="",axes=F)
axis(2,las=1);box()
axis(1,1:11,paste(c("0-4",paste(0:8*10+5,1:9*10+4,sep="-"),"95+"),"\nyears",sep=""),tick=F,cex.axis=0.85)
axis(1,1:12-0.5,rep("",12))
for(i in 1:11) polygon(c(0,1,1,0)-0.5+i,c(0,0,rep(ImmigAgeDist2[i]*100+0.2,2)),col="goldenrod",border="white",lwd=2)
text(1:11,ImmigAgeDist2*100,paste(round(ImmigAgeDist2*100,2),"%",sep=""),pos=3,cex=0.9)
mtext("Age group",1,2.5,cex=0.75)
mtext("B. Age distribution of immigrants (percent of total)",3,0.8,font=2,cex=.9)
box()
###
plot(1950:2100,PrevalenceInImmigrants,type="l",lwd=2,las=1,col=2,
     ylab="",xlab="",ylim=c(2,280),xlim=c(1951,2012))
abline(h=0:10*50,col="grey85")
mtext("Year",1,2.5,cex=0.75)
lines(1950:2100,PrevalenceInImmigrants,col=colZ[3],lwd=2)
mtext("C. Average TB prevalence in sender countries (per 100K)",3,0.8,font=2,cex=0.9)
box()

###
plot(0,0,xlim=c(0.75,9.25),ylim=c(0.1,3),ylab="",xlab="",axes=F)
axis(2,las=1);box()
axis(1,1:9,paste(c("0-4",paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears"),tick=F,cex.axis=0.85)
axis(1,1:10-0.5,rep("",10))
for(i in 1:9) polygon(c(0,1,1,0)-0.5+i,c(0,0,rep(RelActDisAge[i],2)),col=colZ[4],border="white",lwd=2)
text(1:9,RelActDisAge,c(round(RelActDisAge,2)[1:3],"1.00\n(ref)",round(RelActDisAge,2)[5:9]),pos=3,cex=0.9)
mtext("Age group",1,2.5,cex=0.75)
mtext("D. Risk-ratio of active TB relative to 25-34 year olds",3,0.8,font=2,cex=.9)
#####
plot(0:90+0.5,InfByAge*100,type="l",ylim=c(1.4,60),xlim=c(2,70),las=1,
     xlab="",ylab="")
for(i in 1:30) { ai = a+seq(-a,a*5,length.out=30)[i]; lines(0:90+0.5,(1-(b/(b+t))^ai)*100,col="grey90") }
lines(0:90+0.5,InfByAge*100,lwd=2)
mtext("Age",1,2.5,cex=0.75)

for(i in 1:14) { ag = DatNH99f$AGE>=(0:13)[i]*5 & DatNH99f$AGE<(1:14)[i]*5
                 points(2.5+(0:13)[i]*5,sum(DatNH99f$LTBI[ag]*DatNH99f$WEIGHT[ag])/sum(DatNH99f$WEIGHT[ag])*100,col=2,pch=16)
                 lines(rep(2.5+(0:13)[i]*5,2),
                       qbeta(c(1,39)/40,sum(DatNH99f$LTBI[ag]*DatNH99f$WEIGHT[ag])*sum(DatNH99f$WEIGHT[ag])^2/sum(DatNH99f$WEIGHT[ag]^2)^2,
                             (sum(DatNH99f$WEIGHT[ag])-sum(DatNH99f$LTBI[ag]*DatNH99f$WEIGHT[ag]))*sum(DatNH99f$WEIGHT[ag])^2/sum(DatNH99f$WEIGHT[ag]^2)^2)*100,col=2) }
legend("topleft",c("Survey estimates, NHANES 1999-2000","Fitted hyperbolic function: a=0.079, b=3.44",
                   "Variation in curve for uncertainty in LTBI prevalence"),lwd=c(1,2,1),col=c(2,1,"grey75"),
       pch=c(16,NA,NA),cex=1,bg="white")
mtext("E. Percent with latent TB infection, by age group",3,0.8,font=2,cex=0.9)
box()
###

dev.off(); system(paste("open", pdfnam)) # code for ma
#########################

plot(1950:2100,EntrantsWithActiveTB/max(EntrantsWithActiveTB),type="l",lwd=2,las=1,
     ylab="",xlab="Year",ylim=c(0,1.03))
abline(h=0:10*0.2,col="grey85")
lines(1950:2100,EntrantsWithActiveTB/max(EntrantsWithActiveTB),lwd=2)
mtext("Implied time trend in total TB cases\nentering the US, scaled to max = 1.0",3,0.8,font=2,cex=0.85)

######### DR_TB
pdfnam <- "Fig_ImmigDrInputs.pdf"

pdf(file=pdfnam,width=9, height=9)
layout(rbind(c(1,1),c(1,1),c(1,1),2:3,2:3,2:3,2:3,2:3,c(4,4)))

par(mar=c(5.5,2.5,1.5,1),oma=c(0,2.5,1,0))
##
plot(0,0,xlim=c(0.75,9.25),ylim=c(0.0,11.5),ylab="",xlab="",axes=F)
axis(2,las=1,cex.axis=1.1);box()
axis(1,1:9,paste(c("0-4",paste(0:6*10+5,1:7*10+4,sep="-"),"75+"),"\nyears"),tick=F,cex.axis=1.1,line=0.5)
axis(1,1:10-0.5,rep("",10))
for(i in 1:9) polygon(c(0,1,1,0)-0.5+i,c(0,0,rep(FractTxEbyAge[i]*100,2)),col="forestgreen",border="white",lwd=2)
text(1:9,FractTxEbyAge*100,paste(format(round(FractTxEbyAge*100,2),nsmall=2),"%",sep=""),pos=3,cex=1.1)
mtext("Age group",1,3,cex=0.9)
mtext("A. Percent with prior TB treatment history, by age group",3,0.8,font=2,cex=.9)
mtext("Percent",2,3.2,cex=1)

##
par(mar=c(3,2.5,2.5,1))

plot(1,1,xlim=c(1945,2019),ylim=c(0.014,100),log="y",axes=F,xlab="",ylab="")

abline(h=rep(1:10,4)*10^rep(-2:1,each=10),col="grey95")
abline(h=10^rep(-2:2),col="grey85")

for(i in 1:5) lines(1950:2013,pDRtN[,i]*100,lwd=3,col=colZ[i])
axis(1,195:201*10,195:201*10,cex.axis=1.2)
axis(2,10^(-2:2),10^(-2:2),las=1,cex.axis=1.2);box()
text(rep(2012.5,5),pDRtN[64,]*100,paste(format(round(pDRtN[64,]*100,1),nsmall=1),"%",sep=""),pos=4,cex=1.1)
text(1951,100,"100%",pos=2,cex=1.1)
zz <- NULL
for(i in 1:21) zz[i] <- mean(FbMdr2[FbMdr2$Year==i+1992 & FbMdr2$Prev_TB=="No" & FbMdr2$Years_in_US%in%unique(FbMdr2$Years_in_US)[1:2],"MDR"])
points(1993:2013,zz*100,pch=23,col=colZ[4],bg=colZ[5],lwd=1,cex=0.9)
mtext("B. Strain distribution among treatment-naive",3,1,font=2,cex=1)
mtext("Year",1,2.5,cex=0.9)
mtext("Percent of total (log-scale)",2,3.2,cex=1)

###
plot(1,1,xlim=c(1945,2019),ylim=c(0.014,100),log="y",axes=F,xlab="",ylab="")
abline(h=rep(1:10,4)*10^rep(-2:1,each=10),col="grey95")
abline(h=10^rep(-2:2),col="grey85")

for(i in 1:5) lines(1950:2013,pDRtE[,i]*100,lwd=3,col=colZ[i])
axis(1,195:201*10,195:201*10,cex.axis=1.2)
axis(2,10^(-2:2),10^(-2:2),las=1,cex.axis=1.2);box()
text(rep(2012.5,5),pDRtE[64,]*100,paste(format(round(pDRtE[64,]*100,1),nsmall=1),"%",sep=""),pos=4,cex=1.1)
text(1951,100,"100%",pos=2,cex=1.1)
zz <- NULL
for(i in 1:21) zz[i] <- mean(FbMdr2[FbMdr2$Year==i+1992 & FbMdr2$Prev_TB=="Yes" & FbMdr2$Years_in_US%in%unique(FbMdr2$Years_in_US)[1:2],"MDR"])
points(1993:2013,zz*100,pch=23,col=colZ[4],bg=colZ[5],lwd=1,cex=0.9)

mtext("C. Strain distribution among treatment-experienced",3,1,font=2,cex=1)
mtext("Year",1,2.5,cex=0.9)

plot(0,0,xlim=c(1945,2019),ylim=c(0.014,100),axes=F,xlab="",ylab="")
legend("center",c("INH sensitive, RIF sensitive","INH resistant, RIF sensitive",
                  "INH sensitive, RIF resistant","MDR-TB (non-XDR)","XDR-TB","All MDR-TB, TB case reports")[c(1,4,2,5,3,6)],
       col=colZ[c(1,4,2,5,3,4)],pt.bg=colZ[5],lwd=c(rep(3,5),1),lty=c(rep(1,5),0),pch=c(rep(NA,5),23),ncol=3,cex=1.2,xpd=T)

dev.off(); system(paste("open", pdfnam)) # code for ma

#### SAVE ALL INPUTS
######
# assume no HIV, and fb1 by construction

# Dims: ag, tx, tb, dr, and year
#################  TOTAL IMMIGRANTS   #################
# LPR data by country and year for 1950-2013. 
# Early years are only reported as decades, spline interpolation used to fill gaps
# 


ImmigInputs <- list()
ImmigInputs[["TotByYear"]] <- TotalImmigrants/1e6
# age distribution
ImmigInputs[["AgeDist"]] <- ImmigAgeDist2 

#################  LATENT TB (all tx naive, no dr)   #################
Func_LtbiAgeDist <- function(a,b) { 1-(b/(b+c(2.5,1:9*10,100)))^a }
ImmigInputs[["Func_LtbiAgeDist"]] <- Func_LtbiAgeDist 
ImmigInputs[["LtbiPars"]] <- LtbiPars 

#################  ACTIVE TB  #################
# Time trend with prevalence fixed in 2013
PrevalenceInImmigrants <- colSums(ImmigDat2*ImmigDat2prev)/colSums(ImmigDat2)
Prev25_34 <- colSums(FbAgeYear)[4]/(ImmigAgeDist2[4]*sum(TotalImmigrants[-(1:43)]))
PrevTrend25_34 <- PrevalenceInImmigrants/PrevalenceInImmigrants[64]*Prev25_34
ImmigInputs[["PrevTrend25_34"]] <- PrevTrend25_34 

# Age distribution of active disease, rel to 25-34
ImmigInputs[["RR_Active_TB_Age"]] <- RelActDisAge[c(1:9,9,9)]
# fraction tx experienced, by age
ImmigInputs[["Fract_TxE_by_Age"]] <- FractTxEbyAge[c(1:9,9,9)]


#################  DR-TB  #################
## By year, tx naive'

pDRtN3 <- pDRtN2
pDRtN3[,2] <- pDRtN2[,2]*(3/2)
pDRtN3[,1] <- pDRtN3[,1]-(pDRtN3[,2]-pDRtN2[,2])

pDRtN3[,4] <- pDRtN2[,4]*(16/19)
pDRtN3[,1] <- pDRtN3[,1]-(pDRtN3[,4]-pDRtN2[,4])

pDRtE3 <- pDRtE2
pDRtE3[,2] <- pDRtE2[,2]*(5/6)
pDRtE3[,1] <- pDRtE3[,1]-(pDRtE3[,2]-pDRtE2[,2])

pDRtE3[,4] <- pDRtE2[,4]*(8/9)
pDRtE3[,1] <- pDRtE3[,1]-(pDRtE3[,4]-pDRtE2[,4])


ImmigInputs[["DR_TB_by_year_N"]] <- pDRtN3
ImmigInputs[["DR_TB_by_year_E"]] <- pDRtE3

save(ImmigInputs,file="ImmigInputs_4-15-16.rData")
names(ImmigInputs)
##########
load("ImmigInputs_8-26-15.rData") # ImmigInputs

plot(ImmigInputs[[1]][66:151],ylim=c(0,3.5))
xs <- (ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*2.5 + ImmigInputs[[1]][66]
lines(xs)
xs <- (ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*-0.5 + ImmigInputs[[1]][66]
lines(xs)

ImmigInputs[[1]]["2050"]
((ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*-.5 + ImmigInputs[[1]][66])["2050"]
((ImmigInputs[[1]][66:151]-ImmigInputs[[1]][66])*2.5 + ImmigInputs[[1]][66])["2050"]

############

plot(0:90+0.5,InfByAge*100,type="l",ylim=c(1.4,60),xlim=c(2,70),las=1,
     xlab="",ylab="")
for(i in 1:30) { ai = a+seq(-a,a*5,length.out=30)[i]; lines(0:90+0.5,(1-(b/(b+t))^ai)*100,col="grey90") }
lines(0:90+0.5,InfByAge*100,lwd=2)

#####
1-0.99^35

1-0.975^35
rgamma

