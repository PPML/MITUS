
#######   LOAD THE REQUIRED LIBRARY
library(mvtnorm)

source("R/basic_functions.R")

p = c(-0.75,1.5,0.8); # change parameters to change dist (p[3] is a correlation parameter)

##### CREATE A MATRIX OF 1,000,000 RANDOM OBS APPROX BI NORMAL DISTRIBUTED
samp <- invlgt(rmvnorm(1e6,p[c(1,1)],matrix(p[2]*c(1,p[3],p[3],1),2,2)))

dist <- matrix(NA,4,4)
##### SUM UP THESE 1,000,000 INTO THE 4 DISCRETE GROUPS
for(i in 1:4) {
  for(j in 1:4) {
    dist[i,j] <- sum((samp[,1] > (i-1)/4 & samp[,1] <= i/4) & (samp[,2] > (j-1)/4 & samp[ ,2 ]<=j/4))/nrow(samp)
  }
}

barplot(rowSums(dist))
plot(0:4,0:4,col=NA)
for(i in 1:4) points(1:4-0.5,rep(i-0.5,4),cex=dist[i,]*20,pch=16)


