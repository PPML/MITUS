#'an attempt to create a rebalancing matrix that can be used throughout the model
#'this matrix will be dependent on the mortality risk dimension, but not the tb
#'progression risk dimension. This is appropriate because it allows for the appropriate
#'drainage of the tb progression risk groups if there is an external factor that targets
#'those with that certain risk factor.
#'Example: when the risk factor of interest is HIV, if there were a disease that caused the
#'this code is dependent on the MV
#'@param IP
#'@return trans_mat_tot_ages
reblnc<-function(IP){
library(mvtnorm)
lgt <-  function(x) log(x/(1-x));
#setup the variables to determine the number of risk groups to consider
t=1

dist_gen_v<-rep(NA,16)
for(m in 0:3) {
  for(p in 0:3) {
    dist_t0_v[1+m+p*4] <- IP$dist_gen[m+1,p+1];
  } }
names(dist_t0_v)<-paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))

#this distribution will remain the target or goal distribution throughout the whole model
dist_t1_v<-dist_targ<-dist0<-dist_t0_v
trans_mat_tot_ages <-matrix(NA,16,176)
colnames(trans_mat_tot_ages)<-paste0(rep(paste0("ag",1:11), each=16),"_",paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4)))
rownames(trans_mat_tot_ages)<-names(dist_t0_v)
#create an empty vector for the current distribution of the model
# for (t in (1:1201)){
# if(((t-1)%%120)==0){
dist_t1_v[]<-0
#empty logical vectors for the transitions
can_go <- matrix(0,16,16)
rownames(can_go) <- colnames(can_go) <- paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
did_go <- can_go; did_go[,] <- 0
#define can go matrix
for(i in 1:nrow(can_go)){ # i=1
  pi <- rep(0:3,each=4)[i]; mi <- rep(0:3,4)[i]
  can_go[i, rep(0:3,each=4)%in%(pi + -1:1)  & rep(0:3,4)%in%(mi + -1:1)] <- 1
}


for (ag in 1:11){
  did_go <- can_go; did_go[] <- 0
  # [ag+(m+1)*11+(p+1)*44]
  #calculate the dist at time
   x<-c(.0113957,.0123302,.012528,.0124515,.0120711,.0111045,.00928025,.00676634,.00428208,.00267804,.00154398)
  for(m in 0:3) {
    for(p in 0:3) {
       if ((IP$RRmuRF[m+1]*IP$RRmuHR[2])<5){
         if ((4<ag) & (ag>10)){
      dist_t1_v[(1+m)+p*4] = dist_t0_v[(1+m)+p*4]*(1-((IP$mubt[t,ag]*IP$RRmuRF[m+1]*IP$RRmuHR[2])*x[ag])+.005)}
         else if(ag>9){
           dist_t1_v[(1+m)+p*4] = dist_t0_v[(1+m)+p*4]*(1-((IP$mubt[t,ag]*IP$RRmuRF[m+1]*IP$RRmuHR[2])*x[ag])+.01)
         }
         else{
           dist_t1_v[(1+m)+p*4] = dist_t0_v[(1+m)+p*4]*(1-((IP$mubt[t,ag]*IP$RRmuRF[m+1]*IP$RRmuHR[2])*x[ag]))
         }
} else {
      dist_t1_v[(1+m)+p*4] = dist_t0_v[(1+m)+p*4]*(1-(IP$mubt[t,ag]*(5+(.5*m))))}
    } }
#setup other empty vectors for the difference between distributions : diff_i_v
#dist_t0 and dist_t1
diff_i_v <- rep(NA,16); names(diff_i_v) <- names(can_go)

#create an empty 16x16 matrix for the transitions
trans_mat <- matrix(0,16,16)
rownames(trans_mat) <- colnames(trans_mat)<-names(diff_i_v)
dist_i_v <- dist_t1_v
#open a loop to iterate over the process
N=100
# if (ag>7){ N=50} else {N=100}
for (n in 1:N){
  frc<-.1      #this is a tuning parameter that limits the proportion of the population to transfer
  #calculate the difference between current distribution and the goal distribution
  diff_i_v<- dist_i_v-dist_t0_v
  #calculate the transition matrix needed to
  trans_mat[]<-0
  names(trans_mat) <- names(can_go)
  #this loop is probs the most complicated
  for(r in 1:16) {
    for(c in 1:16)  {
      trans_mat[r,c] <- can_go[r,c]*max(0,(diff_i_v[r]-diff_i_v[c]))
    }
  }


  # Adjust transition matrix, 1st scale up rates, 2nd make sure does not sum to over 1
  # approach seems quite sensitive to this value, = fraction of change to
  for(i in 1:16) {
    trans_mat[i,] <- trans_mat[i,] / (dist_i_v[i]+1e-200)*frc #
    trans_mat[i,] <- trans_mat[i,] / max(1.0,sum(trans_mat[i,])) # should be dist0_v?
  }

  # 5 finalize trans mat
  diag(trans_mat) <- 1-rowSums(trans_mat)

  # 6 record absolute transitions, update did_go
  for(i in 1:16) {
    did_go[i,] <- did_go[i,] + dist_i_v[i]*trans_mat[i,]
  }
  diag(did_go) <- 0

  # 7 update vector
  dist_i_v <- dist_i_v%*%trans_mat
}

#update in one step

trans_mat_tot <- did_go
# print(trans_mat_tot)
for(i in 1:16) {
  trans_mat_tot[i,] <- did_go[i,] / (dist_t0_v[i]+1e-200)
}

diag(trans_mat_tot) <- 1-rowSums(trans_mat_tot)

####### create a block matrix to hold the 11 trans_mat matrices

for (i in 1:16){
  for (j in 1:16){
trans_mat_tot_ages[((16*(((t-1)/120)+1))-(16-i)),((16*ag)-(16-j))]<-trans_mat_tot[i,j]
  }
}

}#end of age loop

return(trans_mat_tot_ages)
} #end of function loop
