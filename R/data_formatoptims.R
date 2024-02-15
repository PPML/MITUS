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
    load(paste("/Users/nis100/Desktop/ModelOpts_011924/","Opt_", loc, "_r7_",i,"_2024-", date, ".rda", sep=""))
    opt_all[i,1:nrow(ParamInitZ)] <- o7$par
    opt_all[i,nrow(ParamInitZ)+1]<- o7$value
  }
  saveRDS(opt_all, file=paste("~/MITUS/inst/", loc,"/", loc, "_Optim_all_", length(batches),"_", month, day,".rds", sep = ""),version=2)
}

optim_data_locs<-function(locs, date){
  for (loc in locs){
    model_load(loc)
    tryCatch(optim_data(1:15,loc, date))->x
  }
}

##check optim plots for locs
# all_optim <- data.frame("State" = ordered_locs, "Sample" = rep(0,51))

calib_plots_locs<-function(locs, simp.date, batches=15){
  for (loc in locs){
    print(loc)
    model_load(loc)
    Opt<-readRDS(system.file(paste0(loc,"/", loc, "_Optim_all_", batches,"_", simp.date,".rds"), package="MITUS"))
    posterior<-rep(0,15)

    for(i in 1:15){ posterior[i]<-trunc_number_n_decimals(Opt[i,ncol(Opt)],1);  print(posterior[i]) }
    # mode<-round(getmode(posterior), 2); print(mode)
    # if (mode>1e11){ print(paste(loc, "did not optimize. Check optim manually", sep = " ")); next }
    samp.i<-which(posterior==min(posterior));print(samp.i)
    Opt<-Opt[,-ncol(Opt)]
    results <-calib(samp.i[1], Opt, loc, cex.size = .8)
    ### Create Calibration Targets
    model_calib_outputs(loc, results, 1, simp.date)
    # model_calib_outputs_2020(loc, results, 1, simp.date)
  }}

calib_plots_locs_2020<-function(locs, simp.date, batches=15){
  for (loc in locs){
    print(loc)
    model_load(loc)
    Opt<-readRDS(system.file(paste0(loc,"/", loc, "_Optim_all_", batches,"_", simp.date,".rds"), package="MITUS"))
    posterior<-rep(0,15)
    for(i in 1:15){ posterior[i]<-trunc_number_n_decimals(Opt[i,ncol(Opt)],1);  print(posterior[i]) }
    # mode<-round(getmode(posterior), 2); print(mode)
    # if (mode>1e11){ print(paste(loc, "did not optimize. Check optim manually", sep = " ")); next }
    samp.i<-which(posterior==min(posterior));print(samp.i)
    Opt<-Opt[,-ncol(Opt)]
    results <-calib_2020(samp.i[1], Opt, loc, cex.size = .8)
    ### Create Calibration Targets
    # model_calib_outputs_2020(loc, results, 1, simp.date)
  }}

#### Create the dataFrame of formatted parameters for model run
format_Opt_to_Par <- function(loc,
                              date,
                              batches = 25){
  ### format the date
  month<-strsplit(date, "-")[[1]][1]
  day<-strsplit(date, "-")[[1]][2]
  ### read in the Opt file
  Opt <- readRDS(system.file(paste0(loc,"/", loc, "_Optim_all_", batches,"_", month, day,".rds"), package="MITUS"))
  ### determine the posterior
  posterior<-rep(0,nrow(Opt))
  for(i in 1:nrow(Opt)){ posterior[i]<-trunc_number_n_decimals(Opt[i,ncol(Opt)],1);  print(posterior[i]) }
  samp.i<-which(posterior==min(posterior));print(samp.i)
  Par0 <- Opt[samp.i[1], -ncol(Opt)]
  ##to facilitate comparisons. These first two steps convert these parameters back to their
  ##distributions
  # normal to uniform
  Par2 <- pnorm(Par0,0,1)
  # uniform to true
  Par3 <- Par2
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii] <- Par3
  ### Both MITUS and Tabby2 have been designed to allow for the future calculation
  ### of uncertainty intervals. However, because of this design, we must create
  ### a matrix that repeats the same vector until the design of those uncertainty
  ### intervals have been completed.
  Par <- t(replicate(batches, P))
  ### Save the matrix to a datafile
  saveRDS(Par, file=paste("~/MITUS/inst/", loc,"/", loc, "_Param_all_", nrow(Par),"_", month, day,".rds", sep = ""),version=2)
}
