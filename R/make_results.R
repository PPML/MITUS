# # locs<-c("CA","FL","GA","IL","MA","NJ" ,"NY", "PA","TX", "VA","WA")
# make_results<-function(locs){
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#   CalibDat<-CalibDatState<-readRDS("~/Desktop/GoodStateOptims/CalibDat_2018-11-14.rds")
#   ParamInit_st<-ParamInit<-readRDS("~/Desktop/GoodStateOptims/ParamInit_2018-11-14.rds")
#   StartVal_st<-StartVal<-readRDS("~/Desktop/GoodStateOptims/StartVal_2018-11-14.rds")
#   Inputs<-readRDS(paste0("~/Desktop/GoodStateOptims/",loc,"/",loc,"_ModelInputs_01-24-19.rds"))
#
#    Opt <- readRDS(system.file(paste0(loc,"/",loc,"_opt_all_2018-12-12.rds"), package="MITUS"))
#   Par <- readRDS(paste0("~/Desktop/GoodStateOptims/parAll_",loc,"_10_2018_11_27.rds"))
#
#   #define the cpp reshaper
#   library(tabus)
#   library(inline)
#   cpp_reshaper <- cxxfunction(
#     signature(ResTab='numeric', ResTabus='numeric', ResTabfb='numeric', res_tab2 = 'numeric'),
#     plugin='Rcpp',
#     body=readr::read_file(
#       system.file('inline_cpp/format_restab2.cpp', package='tabus')))
#
# for (intv in 0:8) {
#   # intvs will be a vector of 0s except for (possibly) one activated intervention
#   intvs <- rep(0, 8)
#
#   # If intv is 0, disable all interventions to run the basecase.
#   # Otherwise, activate one of the 8 interventions.
#   if (intv != 0) intvs[intv] <- 1
#   #use the basecase program change
#   defprg<-def_prgchng(Par[10,])
#   # Simulate using the elements of intvs to control
#   # whether or not each intervention is on.
#   new_OutputsInt(loc,Par,n_cores=1,endyr=2050,
#                     Int1 = intvs[[1]],
#                     Int2 = intvs[[2]],
#                     Int3 = intvs[[3]],
#                     Int4 = intvs[[4]],
#                     Int5 = intvs[[5]],
#                     Scen1 = intvs[[6]],
#                     Scen2 = intvs[[7]],
#                     Scen3 = intvs[[8]],
#                     prg_chng = defprg)
#
#   # Do something to save the simulation outcomes
# }
#
#   #format these results as restabs
# ResTabC<-format_as_restab_sm(loc)
#
# #make the empty restabs
# sm_restab_ints <- sm_restab <- make_empty_res_tab2sm()
# sm_restab_ints %<>% mutate_if(is.factor, as.integer) %>% as.matrix
# #run the reshaper
# #small
# sm_restab_ints<- cpp_reshaper(ResTabC[[1]], ResTabC[[2]], ResTabC[[3]], sm_restab_ints)
# sm_restab[,ncol(sm_restab)] <- sm_restab_ints[,ncol(sm_restab_ints)]
#
# saveRDS(sm_restab2,file=paste0("~/MITUS/inst/",loc, "/sm_restab2.rds"))
# #######
# ResTabC<-format_as_restab_bg(loc)
#
# bg_restab_ints <- bg_restab <- make_empty_res_tab2bg()
# bg_restab_ints %<>% mutate_if(is.factor, as.integer) %>% as.matrix
#
# #big
# bg_restab_ints<- cpp_reshaper(ResTabC[[1]], ResTabC[[2]], ResTabC[[3]], bg_restab_ints)
# bg_restab[,ncol(bg_restab)] <- bg_restab_ints[,ncol(bg_restab_ints)]
#
# saveRDS(bg_restab,file=paste0("~/MITUS/inst/",loc, "/bg_restab2.rds"))
# } }
# ###Make the CalibGraphs for Tabby2
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
# tabby_calib_graphs(loc)
#   }
#
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#
#   model_load(loc)
#   format_as_restab(loc)
# }
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#
# # model_load(loc)
# make_results(loc)
#   }
# #
# for (i in 1:length(locs)){
#   loc<-locs[i]
#   print(loc)
#  #onesim(loc)
# make_results(loc)
#   }
