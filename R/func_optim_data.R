#'The function requires the use of the ParamInit Rdata file and the
#'StartValues Rdata file.
#'@name optim_data
#'@param batches dataframe or matrix of starting values data frame
#'@return datafile of all the optimized data parameters
#'@export


optim_data <- function(batches){

opt_all<-matrix(NA,batches,(nrow(ParamInitZ)+1))
cnames<-c(rownames(ParamInitZ), "post_val")
colnames(opt_all)<-cnames
rnames<-rep(NA,batches)
for (j in 1:batches){
rnames[j]<-paste("b_no_", j, sep="")
}
rownames(opt_all)<-rnames

for (i in 1:batches){
  load(paste("/Users/nis100/Desktop/MA_703/Opt_MA_r7_",i,"_2019-07-03.rda", sep=""))
  opt_all[i,1:nrow(ParamInitZ)] <- o7$par
  opt_all[i,nrow(ParamInitZ)+1]<- o7$value
}
US_opt_all<-opt_all
saveRDS(US_opt_all, file=paste("~/MITUS/inst/MA/MA_Optim_all_", batches,"_703.rds", sep = ""))

}

