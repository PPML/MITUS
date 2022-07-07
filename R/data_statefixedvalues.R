#'@name fixed_vals
#'@description function generates new param_init files for a state with values fixed to national levels
#'@param samp_i which row of the optim_data frame to use
#'@return dataset
#'@export
fixed_vals<-function(samp_i, US_opt_all){
  library(dplyr)
  model_load("US")
  # fixed prior values are those which we have large evidence informing
  # the prior value and we are confident using this for all states
  fixed_prior<-c("EffLt","SensSp","pCurPs")
  pr_x<-filter(ParamInit, rownames(ParamInit) %in% fixed_prior)
  pr_x<-as.vector(unlist(pr_x[,1]))
  # model_load loads in the most recent optimized parameters
  P<-Par[1,]
  US_opt <- as.data.frame(P)

  ###CREATE A VECTOR OF THE VALUES THAT WILL BE REPLACED AND
  ###FILL IT WITH THE VALUES SET TO EITHER THE PRIOR OR THE
  ###NATIONAL CALIBRATIONS
  ###THIS CALL IS NOT STATE DEPENDENT --USING CA FOR CONVENIENCE
  # fixed national level parameters are those that we have less data
  # for an would trust the model at the national level to inform
  fixed_national<-c("muIp","TunmuTbAg","RRmuHR","pfast","ORpfast1","ORpfast2",
                    "ORpfastH","ORpfastPI","rslow","rslowH","TunrslowAge", "rfast",
                    "rRecov","rSlfCur", "rrTestHr", "rrTestLrNoTb",  "pImmScen",
                    "TxQualEarly","TunTxQual","RRcurDef")

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
