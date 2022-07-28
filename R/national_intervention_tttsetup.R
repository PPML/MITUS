#a secondary version of all TTT functions that allow for the
#age and nativity distribution of each risk population to be
#included when calculating differential screening rates.

#this function returns default values for ttt which would be total population
#'@name def_ttt_nat_ag
#'@export
def_ttt_nat_ag<-function(){
  ttt_list<-vector("list", 9)
  names(ttt_list)<-c("AgeDistUS", "AgeDistNUS", "NRiskGrp", "FrcScrn",
                     "StartYr", "EndYr", "RRprg", "RRmu", "RRPrev")
  #load in the population age distributions from the model
  load(system.file("US/US_results_1.rda", package="MITUS"))
  results<-out[1,,]
  ttt_list[[1]]<-results[71,33:43]/sum(results[71,33:54])
  ttt_list[[2]]<-results[71,44:54]/sum(results[71,33:54])
  ttt_list[[3]]<-0 #size of population
  ttt_list[[4]]<-0 #fraction screened
  ttt_list[[5]]<-2022
  ttt_list[[6]]<-2050
  ttt_list[[7]]<-1
  ttt_list[[8]]<-1
  ttt_list[[9]]<-1
  return(ttt_list)
}

#this function returns the sampling distribution for each age and nativity group
create_ttt_mdist<-function(ttt_input,results,PV){
  all_samp_rates<-list()
  frc_of_pop<-rep(1,22)
  samp_dist<-matrix(1,22,17)
  x<-matrix(0,22,16)
  #for each of the ttt populations
for (intv in 1:length(ttt_input)){
  # print(paste("intv # =",intv))
  ttt_list<-ttt_input[[intv]]
  start_yr<-as.numeric(ttt_list[[5]])-1949
  # US_dist<-results[start_yr,33:43]
  # NUS_dist<-results[start_yr,44:54]
  #get the appropriate distribution outputs from the MITUS simulation in the start year
ag<-c("0-4", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-94", "95p")
na<-c("US", "NUS")
yo<-0
for (n in 1:2){
  for (a in 1:11){
  y<-grep(paste(ag[a], na[n], sep = "_"), colnames(results))
  # print(paste("age is", a))
  # print(paste("nat is", n))

  dist<-rep(0,16)
  #need to format the start year
  # for (i in 1:length(y)){
  dist<-results[start_yr,y]
  dist<-matrix(dist,4,4)
  colnames(dist) <- paste0("m",0:3) # mortality
  rownames(dist) <- paste0("p",0:3) # progresison
  if (intv>1){
  dist<-dist-matrix(x[((n-1)*11)+a,],4,4)}
  # print(paste("dist is", dist))
  # print(paste0("dist=",dist))
  ##rate ratio for mortality
  mort_dist<-rowSums(dist_gen)
  RF_fact=20
  RRmuRF    <- rep(NA,4);
  names(RRmuRF) <- c("RF1","RF2","RF3","RF4")
  RRmuRF<-exp((0:3)/3*log(RF_fact))
  rrmort<-RRmuRF/sum(RRmuRF*mort_dist)
  ##rate ratio of TB Progression
  #might need to check this bc current applied to odds then converted to probability
  ORpfastRF  <- 40#PV["ORpfastH"] ##riskfactor
  vORpfastRF  <-c(1,1,1,1)
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
  if( sum(dist)!=0){
  fit <- optim(c(1,1),funcB,rrmort_i=rrmort_i,rrprog_i=rrprog_i,rrprog=rrprog,rrmort=rrmort,dist=dist)
  par = fit$par

  #5 calc transition rates for TTT
  # if(length(ttt_input)==1){
  ttt_pop_yr =ttt_list[[3]]*ttt_list[[4]]*ifelse(n==1,ttt_list[[1]][a], ttt_list[[2]][a])
  # } else{
  # ttt_pop_yr = ttt_list[[10]]*ttt_list[[4]]*ifelse(n==1,ttt_list[[1]][a], ttt_list[[2]][a])
  # } # divide by 1e6 since model in millions
  rr_samp <- (exp(par[1])^(0:3)) %*% t(exp(par[2])^(0:3))
  an_samp_rate <- rr_samp * ttt_pop_yr / sum(rr_samp*dist)
  } else {
    an_samp_rate<- matrix(0,4,4)
  }
  for (i in 1:length(an_samp_rate)) an_samp_rate[i]<-min(an_samp_rate[i],1)

  #calculate the # of people sampled from each of the 22 subgroups
  #when we have more than one population, this will need to be removed before
  #entering the next iteration of the loop

  x[((n-1)*11)+a,]<-x[((n-1)*11)+a,]+as.vector(dist *  an_samp_rate)
  # print(paste("samp rate is ", an_samp_rate))
  yo<-(sum(dist * an_samp_rate))
  # print(intv)
  # print(paste("sum is ",yo))
  samp_dist[((n-1)*11)+a,1:16]<-as.vector(an_samp_rate)
  # print(paste(a,n,(ttt_list[["NRiskGrp"]]*ttt_list[["FrcScrn"]]*ifelse(n==1,US_dist[a], NUS_dist[a]))/results[start_yr,(((n-1)*11)+a)+32]))
  samp_dist[((n-1)*11)+a,17]<-ttt_list[["RRPrev"]]

    ########(ttt_list[["NRiskGrp"]]*ttt_list[["FrcScrn"]]*ifelse(n==1,US_dist[a], NUS_dist[a]))/results[start_yr,(((n-1)*11)+a)+32]
#this will be the real implementation after the risk group age/nativity distributions are finalized
  # frc_of_pop<-(ttt_list[["NRiskGrp"]]*ttt_list[["FrcScrn"]]*ifelse(n==1,ttt_list[[1]][a], ttt_list[[2]][a]))/results[start_yr,(((n-1)*11)+a)+32]
  # samp_dist[((n-1)*11)+a,17]<-frc_of_pop
  } }##end of age and nativity loop
#  print(y)
  all_samp_rates[[intv]]<-samp_dist
}##end of population loop
  return(all_samp_rates)
}
