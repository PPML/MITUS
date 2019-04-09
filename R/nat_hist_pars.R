# ###
# #state paramname value
# ##load in all the data values and then reassign their names
# library(dplyr)
# nat_hist<-c("muIp","TunmuTbAg","pfast","ORpfast1","ORpfast2",
#            "ORpfastPI","ORpfastH","rfast","rslow","rslowH",
#            "TunrslowAge","rRecov", "rSlfCur", "SensLt","SpecLt",
#            "SensSn","SensSp","EffLt")
#
# load("/Users/nis100/Desktop/removed/Optim_all_10_2018-10-31.rda")
# US_opt<-as.data.frame(colMeans(par[c(-1,-5),]))
# nat_hist_val<-filter(US_opt, rownames(US_opt) %in% nat_hist )
# row.names(nat_hist_val)<-nat_hist
# col.names(nat_hist_val)<-"Oct31_means"
#
# save(nat_hist_val, file ="~/MITUS/data/nat_hist_val_11-30-2018.rda")
