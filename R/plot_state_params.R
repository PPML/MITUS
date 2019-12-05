# #plot all attempted optimized states parameters
#
# #read in all the optimized parameters
# for (loc in locs){
#   model_load(loc)
#   CA_opt_all<-readRDS(system.file(paste0(loc, "/", loc, "_Optim_all_10_1104.rds"), package="MITUS"))
#   CA_opt_all[,ncol(CA_opt_all)]<-round(CA_opt_all[,ncol(CA_opt_all)],3)
#   m<-getmode(CA_opt_all[,ncol(CA_opt_all)])
#   mode_opt<-  as.matrix(CA_opt_all[CA_opt_all[,ncol(CA_opt_all)]==m,][,-(ncol(CA_opt_all))])
#   CA_opt_all<-melt(pnorm(mode_opt))
#
#
#   }
