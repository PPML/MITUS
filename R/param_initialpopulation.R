init_pop<-function(){
  InitPopN<-readRDS(system.file("US/NCHS_initpop.rds", package="MITUS"))
  InitPopN<-as.matrix(InitPopN)
  return(InitPopN)
}

# NCHS_mort<-c(754.593,60.1,128.1,178.7,358.7,853.9,1901,4104.4,9331.1,20196.9,20196.9)/1e5

# NCHS_mort<-as.matrix(read.csv("~/Desktop/NCHS_BgMort.csv",header = TRUE))[,2:12]
# row.names(NCHS_mort)=read.csv("~/Desktop/NCHS_BgMort.csv")[,1]
#
# n_mort<-matrix(NA,1801,11)
# for (i in 1:11){
#   n_mort[,i]<-SmoCurve(NCHS_mort[,i])
# }
