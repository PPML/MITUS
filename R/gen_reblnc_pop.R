#'Create parameters for population distribution
#'This script generates the necessary parameters to distribute
#'the model population across two multivariate normal distributions
#'

library(mvtnorm)
source("R/basic_functions.R")

#'Define the number of mortality and progression risk groups
num_mRF=4
num_pRF=4
#' Define the size of the target population;
# targ_pop=7.61; #set to approximately 5% of initialpop
# total_pop=sum(InitPop)
#' Define the distribution of the targeted population
params= c(0,1,0.8);
#'Define the number of cut points in the model; standard is 4

cuts <- lgt(0:num_pRF/num_pRF)

#'Use the pmvnorm function to computer the distribution function of the multivariate normal
#'distribution.
pars = c(-2.5,2,0.6); # parameters defining distribution
dist_gen <- matrix(NA,num_mRF,num_pRF)
for(i in 1:num_mRF) {
  for(j in 1:num_pRF) {
    dist_gen[i,j] <- pmvnorm(lower = cuts[c(i,j)],
                             upper = cuts[c(i,j)+1],
                             mean  = pars[c(1,1)],
                             sigma = matrix(pars[2]*c(1,pars[3],pars[3],1),2,2) )[[1]]
  }
}

dist_gen      <- dist_gen/sum(dist_gen)
colnames(dist_gen) <- paste0("p",0:3) # progression
rownames(dist_gen) <- paste0("m",0:3) # mortality

dist_targ <- dist0 <- dist_gen

#'Change distribution
#'Parameters decided by function above;

pars1=params
dist_targ <- matrix(NA,num_mRF,num_pRF)

for(i in 1:num_mRF) {
  for(j in 1:num_pRF) {
    dist_targ[i,j]     <- pmvnorm(lower = cuts[c(i,j)],
                                  upper = cuts[c(i,j)+1],
                                  mean  = pars1[c(1,1)],
                                  sigma = matrix(pars1[2]*c(1,pars1[3],pars1[3],1),2,2) ) [[1]]
  }
}

#'Define weights of population in the target group and gen pop
# wgt_targ=targ_pop/total_pop
# wgt_gen=1-wgt_targ

wgt_targ =0.05
wgt_gen = 0.95
dist_new=wgt_gen*dist_gen + wgt_targ*dist_targ

#'Rebalance to a New Distribution
#'First, setup a matrix of all allowable transitions between states.
#'For this model, individuals can only move a single step in either
#'direction, that is an individual in mortality level 2 can only move
#'to levels 1 or 3, but not to level 4, in a single time step.
#'The matrix can_go below captures these transition limitations. We also
#'create a matrix to hold the sum of transitions, did_go.
#'
can_go <- matrix(0,num_mRF^2,num_pRF^2)
rownames(can_go) <- colnames(can_go) <- paste0(rep(paste0("p",0:(num_pRF-1)),each=num_pRF),"_",rep(paste0("m",0:(num_mRF-1)), num_mRF))
did_go <- can_go
did_go[,] <- 0


for(i in 1:nrow(can_go)){ # i=1
  pi <- rep(0:(num_pRF-1),each=4)[i]
  mi <- rep(0:(num_mRF-1), 4)[i]
  can_go[i, rep(0:(num_pRF-1),each=num_pRF)%in%(pi + -1:1)  & rep(0:(num_mRF-1),num_mRF)%in%(mi + -1:1)] <- 1
}

#'Setup two new matrices for changes
#'The first matrix, dist_goal is the distribution we want, fed in as model input, a 4x1 vector
#'The second matrix, dist_orig is the starting distribution pulled from inside the timestep


dist_goal <- dist_new

#moved to timestep
#dist_orig <- dist_orig_t

#'Create several vectors populated with NA to input into the model;

diff_i_v <- rep(NA,16); names(diff_i_v) <- colnames(can_go)
dist_orig_v <- dist_goal_v <- diff_i_v

dist_goal_v <- rep(0,16)
#'Create a vector from dist_goal, in order to account for age & nativity, we will repeat this 33 times
#'because we haven't a separate distribution yet
for (m in 0:(num_mRF-1)){
  for(p in 0:(num_pRF-1)){
    dist_goal_v[1+m+p*4] <- dist_gen[m+1,p+1]
} }

# moved to time step
# for(m in 0:num_pRF-1) {
#   for(p in 0:num_pRF-1) {
#     dist_goal_v[1+m+p*4] <- dist_goal[m+1,p+1]; ###need to update these 4's for flexibility
#     dist_orig_v[1+m+p*4] <- dist_orig[m+1,p+1];
#   }
# }

# dist_gen[,]<-.0625
# dist_goal_v[]<-.0625
