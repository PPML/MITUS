# dist_gen_v
# dist_new_norm<-dist_new<-matrix(NA,11,16)
# for (ag in 1:11){
#   for (nm in 1:4){
#     for (im in 1:4){
#      # print(nm+(im-1)*4)
#     dist_new[ag,(nm+(im-1)*4)] <-dist_gen_v[(nm+(im-1)*4)]*(1-(prms$mubt[6,ag]*prms$RRmuRF[nm]))
#   }}}
# for (i in 1:11){
#   dist_new_norm[i,]<-dist_new[i,]/rowSums(dist_new)[i]
# }
#
# ## apply the transmat tot calculated in Rcpp
# rblnc_dist<-matrix(NA,11,16)
# for (ag in 1:11){
#   for (nm in 1:4){
#     for (im in 1:4){
#       for (m2 in 1:4){
#         for (p2 in 1:4){
#       # print(nm+(im-1)*4)
#     # print((16*(ag))-(16-(nm+(im-1)*4)))
#      rblnc_dist[ag,nm+(im-1)*4]<-dist_new[ag,(m2+(p2-1)*4)]*trans_mat_tot_ages[(m2+(p2-1)*4),(16*(ag))-(16-(nm+(im-1)*4))]
#   }
#   }
#     }  }
# }
