# ###Taking the TB burden data to PrevTrend form
crude_rate<-function(Inputs){
totcase<-Inputs$ImmigInputs$TBBurdenImmig*(90/1e5)*(Inputs$ImmigInputs$TotByYear[1:69]*1e6)

#RR of TB prevalence across age groups
RR<-Inputs$ImmigInputs[["RR_Active_TB_Age"]]


AgeDist<-Inputs$ImmigInputs$AgeDist

AgeDist[11,]<-.005690661*AgeDist[10,]
AgeDist[10,]<-(1-.005690661)*AgeDist[10,]

TotImmigAge<-matrix(NA,11,69)
for (i in 1:69){
  for (j in 1:11){
    TotImmigAge[j,i]   <- Inputs$ImmigInputs$TotByYear[i]*AgeDist[j,i]
  }}

newprev<-rep(NA,69)
for (j in 1:ncol(Inputs$ImmigInputs$AgeDist[,1:69])){
for (i in 1:length(RR)){
    newprev[j]<-totcase[j]/sum(RR[i]*TotImmigAge[i,j]*1e6)
  }
}

cruderate<-newprev/1e5
  return(cruderate)}
#
# ##what is the total of the old way
#
# # # check that this gives us the right amount
# # x<-rep(NA,69)
# # for (j in 1:ncol(TotImmigAge[,1:69])){
# #   for (i in 1:length(RR)){
# #     x[j]<-sum(RR[i]*newprev[j]*TotImmAge[i,j]*1e6)
# #   }
# # }
#
