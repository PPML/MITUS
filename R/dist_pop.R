library(mvtnorm)

source("R/basic_functions.R")

p = c(-0.75,1.5,0.8); # change parameters to change dist (p[3] is a correlation parameter)

samp <- invlgt(rmvnorm(1e6,p[c(1,1)],matrix(p[2]*c(1,p[3],p[3],1),2,2)))

dist <- matrix(NA,4,4)

for(i in 1:4) {
  for(j in 1:4) {
    dist[i,j] <- sum((samp[,1]>(i-1)/4 & samp[,1]<=i/4) &
                     (samp[,2]>(j-1)/4 & samp[,2]<=j/4))/nrow(samp)
  }
}
