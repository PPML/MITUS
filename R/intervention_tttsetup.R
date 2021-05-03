#'this function returns default values for ttt which would be total population

#'@name def_ttt
#'@export
def_ttt<-function(){
  ttt_list<-vector("list", 9)
  names(ttt_list)<-c("NativityGrp", "AgeGrp", "NRiskGrp", "FrcScrn",
                    "StartYr", "EndYr", "RRprg", "RRmu", "RRPrev")
  ttt_list[[1]]<-"All" #or "USB" or "NUSB"
  ttt_list[[2]]<-"All" #or "0 to 24" or "25 to 64" or "65+"
  ttt_list[[3]]<-0 #size of population IN MILLIONS
  ttt_list[[4]]<-0 #fraction screened
  ttt_list[[5]]<-2020
  ttt_list[[6]]<-2050
  ttt_list[[7]]<-1
  ttt_list[[8]]<-1
  ttt_list[[9]]<-1
  return(ttt_list)
}

#'this function returns default values for ttt which would be total population

#'@name create_ttt_dist
#'@param ttt_list vector of values from interface that define the intervention
#'@param results matrix of results (can be 1 entry from basecase)
#'@param PV formatted parameter vector
#'@return ttt_params list of params for ttt interventions
#'@export
create_ttt_dist<-function(ttt_list,results,PV){
#get the appropriate distribution outputs from the MITUS simulation in the start year
#this returns not the total population size in each of the 16 strata; not the frc in each
  na<-switch(ttt_list[[1]], "All" =c("US","NUS"), "USB"="US","NUSB"="NUS")

  ag<-switch(ttt_list[[2]], "All" =c("0-24","25-64","65\\+"), "0 to 24"="0-24","25 to 64"="25-64", "65+"="65\\+")

  x<-list()
  for (i in 1:length(ag)){
    for(j in 1:length(na)){
  x[[i+(3*(j-1))]]<-grep(paste(ag[i], na[j], sep = "_"), colnames(results))
}}

y<-x[unlist(lapply(x, length) != 0)]

dist<-rep(0,16)
#need to format the start year
start_yr<-as.numeric(ttt_list[[5]])-1949
#distribution in millions
for (i in 1:length(y)){
  dist<-dist+results[start_yr,y[[i]]]
}
dist<-matrix(dist,4,4)
colnames(dist) <- paste0("p",0:3) # progresison
rownames(dist) <- paste0("m",0:3) # mortality

##rate ratio for mortality
  mort_dist<-rowSums(dist_gen)
  RF_fact=20
  RRmuRF    <- rep(NA,4);
  names(RRmuRF) <- c("RF1","RF2","RF3","RF4")
  RRmuRF<-exp((0:3)/3*log(RF_fact))
  rrmort<-RRmuRF/sum(RRmuRF*mort_dist)
##rate ratio of TB Progression
  #might need to check this bc current applied to odds then converted to probability
  ORpfastRF  <- 40 ##riskfactor
  # vORpfastRF  <-c(1,1,1,1)
  rrprog  <-(exp((0:3)/3*log(ORpfastRF)))
#desired RR for the screening groups
  rrprog_i <- ttt_list[[7]]
  rrmort_i <- ttt_list[[8]]
#create function for reweighting
  funcB <- function(par,rrmort_i,rrprog_i,rrprog,rrmort,dist){ # par = 2:3
    rr_samp <- (exp(par[1])^(0:3)) %*% t(exp(par[2])^(0:3))
    dist_i  <- dist * rr_samp / sum(dist * rr_samp) * sum(dist)
    (sum(colSums(dist_i)*rrprog)/sum(colSums(dist)*rrprog)-rrprog_i)^2 + (sum(rowSums(dist_i)*rrmort)/sum(rowSums(dist)*rrmort)-rrmort_i)^2 + diff(par)^2/100
  }
#apply
  fit <- optim(c(1,1),funcB,rrmort_i=rrmort_i,rrprog_i=rrprog_i,rrprog=rrprog,rrmort=rrmort,dist=dist)
  par = fit$par

  #5 calc transition rates for TTT
  ttt_pop_yr =ttt_list[[3]]*ttt_list[[4]] # divide by 1e6 since model in millions
  rr_samp <- (exp(par[1])^(0:3)) %*% t(exp(par[2])^(0:3))
  an_samp_rate <- rr_samp * ttt_pop_yr / sum(rr_samp*dist)
  ttt_params<-list()
  ttt_params[['an_samp_rate']]<-pmin(an_samp_rate,12)
  #what is the fraction of the total treatment native population?
  ttt_params[['frc_of_totpop']]<-(ttt_list[["NRiskGrp"]]*ttt_list[["FrcScrn"]])/results[start_yr,683]

  return(ttt_params)
}

