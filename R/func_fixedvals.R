#'@name fix_vals
#'@description function generates new param_init files for a state with values fixed to national levels
#'@param samp_i which row of the optim_data frame to use
#'@return dataset
#'@export
fixed_vals<-function(samp_i){
#load the most up to date national optimized data set
# and format the data back to their original distributions
load("~/MITUS/data/US_opt_all_2018-12-10.rda")
Par<-US_opt_all[samp_i,-(ncol(US_opt_all)-1)]
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
model_load("CA")
Inputs<-CA_Inputs
rm(CA_Inputs)
ParamInit<-ParamInit_st
fixed_prior<-c("SensLt","SpecLt","EffLt","SensSp","pCurPs")

fixed_national<-c("muIp","TunmuTbAg","RRmuHR","pfast","ORpfast1","ORpfast2",
                  "ORpfastH","ORpfastPI","rfast","rslow","rslowH","TunrslowAge",
                  "rRecov","rSlfCur","TxQualEarly","TunTxQual","RRcurDef")

pr_x<-filter(ParamInit_st, rownames(ParamInit_st) %in% fixed_prior)
pr_x<-as.vector(unlist(pr_x[,1]))

# US_opt<-as.data.frame(US_opt_all[10,])

nat_x<-filter(US_opt, rownames(US_opt) %in% fixed_national)
nat_x<-as.vector(unlist(nat_x))
fixed_vals<-rep(NA,length(c(fixed_prior,fixed_national)))
fixed_vals<-as.vector(c(pr_x,nat_x))
names(fixed_vals)<-c(fixed_prior,fixed_national)

###LOAD IN THE STATE DATA TO CREATE A NEW PARAM_INIT,
###START_VAL DATA SETS FOR A NEW OPTIMIZATION RUN
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

ParamInit_st[,]<-newparaminitst[,]
save(ParamInit_st,file=paste0("~/MITUS/data/ST_ParamInit_", Sys.Date() ,"rda"))

}
