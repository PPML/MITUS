#' #how to design a multiple ttt function
#' ttt_functions<-function(num_intv,list_intv){
#' try(if(num_intv != length(list_intv)) stop("Number of Interventions and Intervention Parameters do not match"))
#' ttt_param_list<-vector("list",length(list_intv))
#' names(ttt_param_list)<-names(list_intv)
#'   for(i in 1:length(list_intv)){
#'   ttt_list<-list_intv[[i]]
#' ################################################################################
#' #################### TTT ADDITIONAL SCREENING PROBABILITIES ####################
#' ################################################################################
#' ttt_month <-seq((ttt_list[["StartYr"]]-1950)*12,(ttt_list[["EndYr"]]-1949)*12,1)
#' ttt_sampling_dist<-matrix(0,4,4)
#' ttt_na<-99
#' ttt_ag<-99
#' ttt_pop_frc<-0
#'
#' if (ttt_list[[3]]!=0 & ttt_list[[4]]!=0){
#'   load(system.file("US/US_results_1.rda", package="MITUS"))
#'   # res<-out[1,,]
#'   # colnames(res)<-func_ResNam()
#'   x<-create_mult_ttt_dist(ttt_param_list = ttt_input,
#'                      results = out[1,,],
#'                      PV = PV)
#'   # if (ttt_list[[7]]!=1 | ttt_list[[8]]!=1){
#'   ttt_sampling_dist<-x[[1]]/12
#'   # }
#'   # ttt_pop_frc<-x[[2]]
#'   ttt_ag<-switch(ttt_list[["AgeGrp"]], "All"=0:10,
#'                  "0 to 24"=0:2,
#'                  "25 to 64"=3:6,
#'                  "65+"=7:10
#'   )
#'   ttt_na<-switch(ttt_list[["NativityGrp"]], "All"=0:2,
#'                  "USB"=0,
#'                  "NUSB"=1:2
#'   )
#'
#' }
#'   ttt_ltbi<-ttt_list[[9]]
#'   ttt_param_list[[i]]<-list(ttt_na,
#'                             ttt_ag,
#'                             ttt_sampling_dist,
#'                             ttt_pop_frc,
#'                             ttt_month,
#'                             ttt_ltbi)
#'   names(ttt_param_list[[i]])<-c(   "ttt_na",
#'                                    "ttt_ag",
#'                                    "ttt_sampling_dist",
#'                                    "ttt_pop_frc",
#'                                    "ttt_month",
#'                                    "ttt_ltbi")
#'   }
#' return(ttt_param_list)
#' }
#' #this function returns default values for ttt which would be total population
#' #'@name create_mult_ttt_dist
#' #'@param ttt_param_list vector of values from interface that define the intervention
#' #'@param results matrix of results (can be 1 entry from basecase)
#' #'@param PV formatted parameter vector
#' #'@export
#' create_mult_ttt_dist<-function(ttt_param_list,results,PV){
#'   #create some empty objects to fill
#'   ttt_sampling_dists<-matrix(1,length(ttt_param_list),16)
#'   ttt_ltbi<-frc_of_pop<-rep(1,length(ttt_param_list))
#'   #loop through all the entries of the list
#'   for (intv in 1:length(ttt_param_list)){
#'   ttt_list<-ttt_param_list[[intv]]
#'   #get the appropriate distribution outputs from the MITUS simulation in the start year
#'   na<-switch(ttt_list[[1]], "All" =c("US","NUS"), "USB"="US","NUSB"="NUS")
#'
#'   ag<-switch(ttt_list[[2]], "All" =c("0-24","25-64","65p"), "0 to 24"="0-24","25 to 64"="25-64", "65+"="65p")
#'
#'   x<-list()
#'   for (i in 1:length(ag)){
#'     for(j in 1:length(na)){
#'       x[[i+(3*(j-1))]]<-grep(paste(ag[i], na[j], sep = "_"), colnames(results))
#'     }}
#'
#'   y<-x[unlist(lapply(x, length) != 0)]
#'
#'   dist<-rep(0,16)
#'   #need to format the start year
#'   start_yr<-as.numeric(ttt_list[[5]])-1950
#'   for (i in 1:length(y)){
#'     dist<-dist+results[start_yr,y[[i]]]
#'   }
#'   dist<-matrix(dist,4,4)
#'   colnames(dist) <- paste0("p",0:3) # progresison
#'   rownames(dist) <- paste0("m",0:3) # mortality
#'
#'   ##rate ratio for mortality
#'   mort_dist<-rowSums(dist_gen)
#'   RF_fact=20
#'   RRmuRF    <- rep(NA,4);
#'   names(RRmuRF) <- c("RF1","RF2","RF3","RF4")
#'   RRmuRF<-exp((0:3)/3*log(RF_fact))
#'   rrmort<-RRmuRF/sum(RRmuRF*mort_dist)
#'   ##rate ratio of TB Progression
#'   #might need to check this bc current applied to odds then converted to probability
#'   ORpfastRF  <- PV["ORpfastH"] ##riskfactor
#'   vORpfastRF  <-c(1,1,1,1)
#'   rrprog  <-(exp((0:3)/3*log(ORpfastRF)))
#'   #desired RR for the screening groups
#'   rrprog_i <- ttt_list[[7]]
#'   rrmort_i <- ttt_list[[8]]
#'   #create function for reweighting
#'   funcB <- function(par,rrmort_i,rrprog_i,rrprog,rrmort,dist){ # par = 2:3
#'     rr_samp <- (exp(par[1])^(0:3)) %*% t(exp(par[2])^(0:3))
#'     dist_i  <- dist * rr_samp / sum(dist * rr_samp) * sum(dist)
#'     (sum(colSums(dist_i)*rrprog)/sum(colSums(dist)*rrprog)-rrprog_i)^2 + (sum(rowSums(dist_i)*rrmort)/sum(rowSums(dist)*rrmort)-rrmort_i)^2 + diff(par)^2/100
#'   }
#'   #apply
#'   fit <- optim(c(1,1),funcB,rrmort_i=rrmort_i,rrprog_i=rrprog_i,rrprog=rrprog,rrmort=rrmort,dist=dist)
#'   par = fit$par
#'
#'   #5 calc transition rates for TTT
#'   ttt_pop_yr =ttt_list[[3]]*ttt_list[[4]] # divide by 1e6 since model in millions
#'   rr_samp <- (exp(par[1])^(0:3)) %*% t(exp(par[2])^(0:3))
#'   an_samp_rate <- rr_samp * ttt_pop_yr / sum(rr_samp*dist)
#'
#'   ttt_sampling_dists[intv,]<-as.vector(an_samp_rate)
#'   frc_of_pop[intv]<-(ttt_list[["NRiskGrp"]]*ttt_list[["FrcScrn"]])/results[start_yr,2]
#'   ttt_ltbi[intv]<-ttt_list[["RRPrev"]]
#'   }
#'   ttt_params<-list()
#'   ttt_params[['an_samp_rate']]<-ttt_sampling_dists
#'   ttt_params[['frc_of_totpop']]<-frc_of_pop
#'   ttt_params[['ttt_ltbi']]<-ttt_ltbi
#'   return(ttt_params)
#' }
