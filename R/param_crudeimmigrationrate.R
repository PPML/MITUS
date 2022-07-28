# ###Taking the TB burden data to PrevTrend form
crude_rate<-function(Inputs, loc,r_decline=0.015){
# totcase<-Inputs$ImmigInputs$TBBurdenImmig*(90/1e5)*(Inputs$ImmigInputs$TotByYear[1:69]*1e6)

#RR of TB prevalence across age groups
# RR<-Inputs$ImmigInputs[["RR_Active_TB_Age"]]
#
# AgeDist<-as.matrix(Inputs$ImmigInputs$AgeDist)
# nonzeroages<-which(rowSums(AgeDist)!=0)
# #calculate the age specific population
# TotImmigAge<-matrix(NA,11,69)
# for (i in 1:69){
#   for (j in 1:11){
#     TotImmigAge[j,i]   <- Inputs$ImmigInputs$TotByYear[i]*AgeDist[j,i]
#   }}
#
# newprev<-rep(NA,69)
# for (j in 1:ncol(Inputs$ImmigInputs$AgeDist[,1:69])){
# for (i in nonzeroages){
#     newprev[j]<-totcase[j]/sum(RR[i]*TotImmigAge[i,j]*1e6)
#   }
# }
#
# cruderatepast<-newprev/1e5
#
# for (i in 1:5){
#   if((cruderatepast[i+64]/cruderatepast[(i+64)-1]) < .95){
#     cruderatepast[i+64]<-.95*cruderatepast[(i+64)-1]
#   }
# }
cruderatepast<-Inputs$ImmigInputs$TBBurdenImmig*(90/1e5)
cruderatefuture<-rep(NA,151-69)
cruderatefuture[1]<-cruderatepast[69]*(1-r_decline)
for (i in 2:length(cruderatefuture)){
  cruderatefuture[i]<-cruderatefuture[i-1]*(1-r_decline)
}
#combine past and future
cruderate<-c(cruderatepast,cruderatefuture)
names(cruderate)<-as.character(seq(1950,2100,1))


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
