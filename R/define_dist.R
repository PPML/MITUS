#define the population

num_mRF = 4
num_pRF = 4
pars = c(-2.5,2,0.6)
cuts <- lgt(0:4/4)
dist_gen<-matrix(NA,4,4)
dist_gen_v<-rep(NA,16)
for(i in 1:num_mRF) {
  for(j in 1:num_pRF) {
    dist_gen[i,j] <- pmvnorm(lower = cuts[c(i,j)],
                             upper = cuts[c(i,j)+1],
                             mean  = pars[c(1,1)],
                             sigma = matrix(pars[2]*c(1,pars[3],pars[3],1),2,2) )[[1]]
    dist_gen_v[(i)+((j-1)*4)]=dist_gen[i,j]
  }
}
names(dist_gen_v) <-paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))



can_go <- matrix(0,16,16)
rownames(can_go) <- colnames(can_go) <- paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))

#define can go matrix
for(i in 1:nrow(can_go)){ # i=1
  pi <- rep(0:3,each=4)[i]; mi <- rep(0:3,4)[i]
  can_go[i, rep(0:3,each=4)%in%(pi + -1:1)  & rep(0:3,4)%in%(mi + -1:1)] <- 1
}

# can_go[1:4,5:16]<-0
# can_go[5:8,1:4]<-0; can_go[5:8,9:16]<-0;
# can_go[9:12,1:8]<-0; can_go[9:12,13:16]<-0;
# can_go[13:16,1:12]<-0

HRdist<-c(.0113957,.0123302,.012528,.0124515,.0120711,.0111045,.00928025,.00676634,.00428208,.00267804,.00154398)
