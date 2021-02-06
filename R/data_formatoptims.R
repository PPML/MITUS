#'The function requires the use of the ParamInit Rdata file and the
#'StartValues Rdata file.
#'@name optim_data
#'@param batches dataframe or matrix of starting values data frame
#'@param loc two digit location code
#'@param date MM-DD date on optim files
#'@return datafile of all the optimized data parameters
#'@export

optim_data <- function(batches, loc, date){
opt_all<-matrix(NA,max(batches),(nrow(ParamInitZ)+1))
cnames<-c(rownames(ParamInitZ), "post_val")
colnames(opt_all)<-cnames
rnames<-rep(NA,max(batches))
for (j in 1:max(batches)){
rnames[j]<-paste("b_no_", j, sep="")
}
rownames(opt_all)<-rnames
month<-strsplit(date, "-")[[1]][1]
day<-strsplit(date, "-")[[1]][2]
for (i in batches){
  load(paste("/Users/nis100/Desktop/US_020421/","Opt_", loc, "_r7_",i,"_2021-", date, ".rda", sep=""))
  opt_all[i,1:nrow(ParamInitZ)] <- o7$par
  opt_all[i,nrow(ParamInitZ)+1]<- o7$value
}
saveRDS(opt_all, file=paste("~/MITUS/inst/", loc,"/", loc, "_Optim_all_", length(batches),"_", month, day,".rds", sep = ""))
}

optim_data_locs<-function(locs, date){
  for (loc in locs){
    model_load(loc)
    optim_data(1:10,loc, date)
  }
}

##check optim plots for locs
calib_plots_locs<-function(locs, simp.date, batches=10){
  for (loc in locs){
    model_load(loc)
    Opt<-readRDS(system.file(paste0(loc,"/", loc, "_Optim_all_", batches,"_", simp.date,".rds"), package="MITUS"))
    posterior<-round(Opt[,ncol(Opt)],2); print(posterior)
    mode<-round(getmode(posterior), 2);print(mode)
    if (mode>1e11){ print(paste(loc, "did not optimize. Check optim manually", sep = " ")); next }
    samp.i<-sample(which(posterior==mode), 1);print(samp.i)
    Opt<-Opt[,-ncol(Opt)]
    calib(samp.i, Opt, loc, cex.size = .65)
  }
}
