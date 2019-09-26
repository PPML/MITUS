# ttt_list<-def_ttt()
# # # for (i in 1:length(loc_vec)){
# # # loc<-loc_vec[i]
# model_load(loc)
# results_ttt<-new2_OutputsInt(loc=loc,ParMatrix=Par[1,],n_cores = 1, endyr=2050,
#                               Int1 = 0, Int2 = 0,Int3 = 0,Int4 = 0,Int5 = 0,
#                               Scen1 = 0, Scen2 = 0, Scen3 = 0,
#                               prg_chng = def_prgchng(Par[1,]), ttt_list = ttt_list)
# results_pc<-new_OutputsInt(loc=loc,ParMatrix=Par[1,], n_cores = 1, endyr=2050,
#                            Int1 = 0, Int2 = 0,Int3 = 0,Int4 = 0,Int5 = 0,
#                            Scen1 = 0, Scen2 = 0, Scen3 = 0,prg_chng = def_prgchng(Par[1,]))
# print(loc)
# data_name <-load(system.file(paste0(loc, "/", loc, "_results_",1,".rda"), package='MITUS'))
# results_pre<-(get(data_name))
# results_pre<-results_pre[1:2,,]
# if(identical(results_ttt, results_pc)==FALSE){
#   print("ttt and pc are different")
# }
# if(identical(results_ttt, results_pre)==FALSE){
#   print("ttt and presim are different")
# }
# if(identical(results_pc, results_pre)==FALSE){
#   print("pc and presim are different")
# }
# # }
# #
# # ##is this a reshaper problem
# # # results_ttt<-list(results_ttt)
# # # results_pc<-list(results_pc)
# # # results_pre<-list(results_pre)
# # #
# # # rr_ttt<-reshape_results(loc,results_ttt)
# # # rr_pc<-reshape_results(loc,results_pc)
# # # rr_pre<-reshape_results(loc,results_pre)
# # #
# # #
# # # if(identical(rr_ttt, rr_pc)==FALSE){
# # #   print("ttt and pc are different")
# # # }
# # # if(identical(rr_ttt, rr_pre)==FALSE){
# # #   print("ttt and presim are different")
# # # }
# # # if(identical(rr_pc, rr_pre)==FALSE){
# # #   print("pc and presim are different")
# # # }
