#'@name fix_vals
#'@description function generates new param_init files for a state with values fixed to national levels
#'@param samp_i which row of the optim_data frame to use
#'@return dataset
#'@export
fixed_vals<-function(samp_i, US_opt_all){
  library(dplyr)
  model_load("US")
  fixed_prior<-c("EffLt","SensSp","pCurPs", "pImmScen")
  pr_x<-filter(ParamInit, rownames(ParamInit) %in% fixed_prior)
  pr_x<-as.vector(unlist(pr_x[,1]))
#load the most up to date national optimized data set
# and format the data back to their original distributions
# US_opt_all<-readRDS(system.file("US/US_Optim_all_10_1031.rds", package="MITUS"))
# Par<-US_opt_all[samp_i,-(ncol(US_opt_all))]
  # Par<-o7$par
Opt <- readRDS("~/MITUS/inst/US/US_Optim_all_10_0417.rds")
Par <-Opt[5,-ncol(Opt)]
Par2 <- pnorm(Par,0,1)
# uniform to true
Par3 <- Par2
Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
P[ii] <- Par3
P<-P
US_opt <- as.data.frame(P)

###CREATE A VECTOR OF THE VALUES THAT WILL BE REPLACED AND
###FILL IT WITH THE VALUES SET TO EITHER THE PRIOR OR THE
###NATIONAL CALIBRATIONS
###THIS CALL IS NOT STATE DEPENDENT --USING CA FOR CONVENIENCE
fixed_national<-c("muIp","TunmuTbAg","RRmuHR","pfast","ORpfast1","ORpfast2",
                  "ORpfastH","ORpfastPI","rslow","rslowH","TunrslowAge", "rfast",
                  "rRecov","rSlfCur","TxQualEarly","TunTxQual","RRcurDef")

nat_x<-filter(US_opt, rownames(US_opt) %in% fixed_national)
nat_x<-as.vector(unlist(nat_x))
fixed_vals<-rep(NA,length(c(fixed_prior,fixed_national)))
fixed_vals<-as.vector(c(pr_x,nat_x))
names(fixed_vals)<-c(fixed_prior,fixed_national)

###LOAD IN THE STATE DATA TO CREATE A NEW PARAM_INIT,
###START_VAL DATA SETS FOR A NEW OPTIMIZATION RUN
# model_load("CA")
#read in the baseline state data from the csv in extdata
ParamInit_st<-read.csv(system.file("extdata", "ParamInitST_final.csv", package="MITUS"))[,2:9]
rownames(ParamInit_st)<-read.csv(system.file("extdata", "ParamInitST.csv", package="MITUS"))[,1]
newparaminitst<-ParamInit_st[,]
#IF THE COLNAME OF PARAM INIT IS IN THE fixed_vals
#VECTOR THEN MAKE SURE THAT ParamInit_st$Calib <-0
bool<-rownames(ParamInit_st) %in% names(fixed_vals)
newparaminitst[bool,"Calib"]<-0

#THEN REPLACE THE FIRST COLUMN WITH THE VALUES IN
#THE fixed_val VECTOR

#IF ROWNAME IN PARAMINIT_ST MATCHES NAME FIXED_VAL
#THEN REPLACE THE FIRST COL VALUE WITH THE VALUE IN FIXED_VAL
for (i in seq_along(fixed_vals)) {
  j <- which(rownames(ParamInit_st) == names(fixed_vals)[[i]])
  newparaminitst[j,1] <- fixed_vals[[i]]
}

#check that everything matches
y<-rownames(US_opt) %in% fixed_national
x<-rownames(newparaminitst) %in% fixed_national
identical(US_opt[y,1],newparaminitst[x,1])

#set TunNetMigration to be Calibrated
newparaminitst["TunNetMig","Calib"]<-1
ParamInit_st[,]<-newparaminitst[,]
saveRDS(ParamInit_st,file=paste0("~/MITUS/inst/ST/ST_ParamInit_", Sys.Date() ,".rds"),version=2)

}
