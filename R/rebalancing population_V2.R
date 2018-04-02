mTrsp <- function(cl,a)  { apply(col2rgb(cl), 2, function(x){ rgb(x[1],x[2],x[3],a,maxColorValue=255)}) }


library(mvtnorm)
lgt <-  function(x) log(x/(1-x));  invlgt <- function(x) 1/(1+exp(-x))
cuts <- lgt(0:4/4)

## Make 4x4 distribution based on invLogit 2D Normal
pars = c(-1,2,0.6); # parameters defining distribution
dist <- matrix(NA,4,4)
for(i in 1:4) { 
  for(j in 1:4) {
   dist[i,j] <- pmvnorm(lower = cuts[c(i,j)],
                        upper = cuts[c(i,j)+1], 
                        mean  = pars[c(1,1)],
                        sigma = matrix(pars[2]*c(1,pars[3],pars[3],1),2,2) )[[1]]
  }
}
dist <- dist/sum(dist)

####  CHANGES
colnames(dist) <- paste0("p",0:3) # progresison
rownames(dist) <- paste0("m",0:3) # mortality
dist1 <- dist0 <- dist

# Example change to distribution no. 1 = Mortality, weighted towards the high mort group
  #  for(i in 1:4) dist1[i,] <- dist0[i,]*exp(-exp(0:3-4))
  #  dist1 <- dist1/sum(dist1)

# Example change to distribution no. 2 = "HIV" added as a mixture distribution of two Normals 
pars1 = c(0.5,1,0.8); #  parameters for HIV dist
dist1 <- matrix(NA,4,4)
for(i in 1:4) { 
  for(j in 1:4) {
    dist1[i,j] <- pmvnorm(lower = cuts[c(i,j)],
                          upper = cuts[c(i,j)+1], 
                          mean  = pars1[c(1,1)],
                          sigma = matrix(pars1[2]*c(1,pars1[3],pars1[3],1),2,2) )[[1]]
  }
}
dist1 <- dist1/sum(dist1)

# look at HIV vs general
 plot(0:4,0:4,col=NA)
 for(i in 1:4) points(1:4-0.5,rep(i-0.5,4),cex=dist[i,]*20,pch=16,col="grey40")
 for(i in 1:4) points(1:4-0.5,rep(i-0.5,4),cex=dist1[i,]*20,pch=16,col=mTrsp(2,100))
 
 # make the target distribution
 dist_new <- 0.95*dist + 0.05*dist1
 
######### CODE FOR REBALANCEING TO A NEW DISTRIBUTION
 
# 1 Set up matrix of allowable transitions (can_go), 
 #  and matrix for holding sum of transitions (did_go)
can_go <- matrix(0,16,16)
rownames(can_go) <- colnames(can_go) <- paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
did_go <- can_go; did_go[,] <- 0

for(i in 1:nrow(can_go)){ # i=1
  pi <- rep(0:3,each=4)[i]; mi <- rep(0:3,4)[i]
  can_go[i, rep(0:3,each=4)%in%(pi + -1:1)  & rep(0:3,4)%in%(mi + -1:1)] <- 1 
}

# Set up for changes
dist_goal <- dist_new # dist_goal = distribution we want, fed in as model input (probably as a 4x1 vector)
dist_orig <- dist     # dist_orig = starting distribution, inside timestep

# Turn 2X2 matrices as 1X4 vectors
#  dist_i_v is an intermediate distribution, iteratively updated as 
#  we head towards dist_goal_v
diff_i_v <- rep(NA,16); names(diff_i_v) <- colnames(can_go)
dist_orig_v <- dist_goal_v <- diff_i_v
for(m in 0:3) {
  for(p in 0:3) { 
    dist_goal_v[1+m+p*4] <- dist_goal[m+1,p+1]; 
    dist_orig_v[1+m+p*4] <- dist_orig[m+1,p+1];
  }
}

# Everything below will need to be in timestep (ie C++)
trans_mat <- matrix(0,16,16)
dist_i_v <- dist_orig_v

# code to plot results for each update (just for checking things during development)
plot(dist_orig_v-dist_goal_v,type="l")
abline(h=0,col="grey80",lty="11")

# Start loop
N <- 30 # no iterations 
for(n in 1:N){ # n = 1

  # calc distance from dist_i to dist_goal
  diff_i_v <- dist_i_v - dist_goal_v

  # Create transition matrix
  trans_mat[,] <- 0
  colnames(trans_mat) <- rownames(trans_mat) <- rownames(can_go)
  for(r in 1:16) {
    for(c in 1:16)  {
      trans_mat[r,c] <- can_go[r,c]*max(0,(diff_i_v[r]-diff_i_v[c]))
    }
  }

  # Adjust transition matrix, 1st scale up rates, 2nd make sure does not sum to over 1
  frc <- 0.1  # approach seems quite sensitive to this value, = fraction of change to 
  for(i in 1:16) {
    trans_mat[i,] <- trans_mat[i,] / dist_i_v[i]*frc  # 
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
  
  lines(as.numeric(dist_i_v)-dist_goal_v,col=2)
  
}

# final
lines(as.numeric(dist_i_v)-dist_goal_v,col=4)

### Now, update in one step
trans_mat_tot <- did_go
for(i in 1:16) {
  trans_mat_tot[i,] <- did_go[i,] / dist_orig_v[i]
}
diag(trans_mat_tot) <- 1-rowSums(trans_mat_tot)

# check  
  points(as.numeric(dist_orig_v%*%trans_mat_tot)-dist_goal_v,col=4)
  # win!

# Now finally update distribution
  dist_new <- dist_orig
  dist_new[,] <- 0
  for(m in 0:3) {
    for(p in 0:3) { 
      for(m2 in 0:3) {
        for(p2 in 0:3) {
          dist_new[m+1,p+1] <- dist_new[m+1,p+1]+dist_orig[m2+1,p2+1]*trans_mat_tot[1+m2+p2*4,1+m+p*4]
        }
      }
    }
  }
  
  sum((dist_goal-dist_new)^2)/
  sum((dist_goal-dist_orig)^2)
  # SSE in new dist as fraction of SSE in original
  
  