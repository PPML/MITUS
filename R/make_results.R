# # locs<-c("CA","FL","GA","IL","MA","NJ" ,"NY", "PA","TX", "VA","WA")
make_results<-function(locs){
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
for (intv in 0:8) {
  # intvs will be a vector of 0s except for (possibly) one activated intervention
  intvs <- rep(0, 8)

  # If intv is 0, disable all interventions to run the basecase.
  # Otherwise, activate one of the 8 interventions.
  if (intv != 0) intvs[intv] <- 1
  #use the basecase program change
  defprg<-def_prgchng(Par[10,])
  # Simulate using the elements of intvs to control
  # whether or not each intervention is on.
  new_OutputsInt(loc,Par,n_cores=1,endyr=2050,
                    Int1 = intvs[[1]],
                    Int2 = intvs[[2]],
                    Int3 = intvs[[3]],
                    Int4 = intvs[[4]],
                    Int5 = intvs[[5]],
                    Scen1 = intvs[[6]],
                    Scen2 = intvs[[7]],
                    Scen3 = intvs[[8]],
                    prg_chng = defprg)

  # Do something to save the simulation outcomes
}
#load tabus
library(tabus)

#load mitus sims into data structure
load_US_data <- function(i) {
  data_name <-
    load(system.file(paste0("US/US_results_",i,".rda"), package='MITUS'))
  return(get(data_name))
}

US_results <- lapply(1:9, load_US_data)

#make some lists of reformatted data
ResTabC <- list()
ResTabC[['small_results']] <- format_as_restab_small_ages(US_results)
ResTabC[['big_results']] <- format_as_restab_big_ages(US_results)

#Average the results
ResTabC[['small_results']] <- mean_small_restabs(ResTabC, nr = 10, nints = 9)
ResTabC[['big_results']] <- mean_big_restabs(ResTabC, nr = 10, nints = 9)
#make empty restabs
intvs <- c('base_case', paste0('intervention_', 1:5), paste0('scenario_', 1:3))
restab_sm <- restab_ints_sm <- make_empty_res_tab2sm(intvs)
restab_ints_sm %<>% mutate_if(is.factor, as.integer) %>% as.matrix
restab_bg <- restab_ints_bg <- make_empty_res_tab2bg(intvs)
restab_ints_bg %<>% mutate_if(is.factor, as.integer) %>% as.matrix

#define the reshaper
library(inline)
cpp_reshaper <- cxxfunction(
  signature(ResTab='numeric', ResTabus='numeric', ResTabfb='numeric', res_tab2 = 'numeric'),
  plugin='Rcpp',
  body=readr::read_file(
    system.file('inline_cpp/format_restab2.cpp', package='tabus')))

restab_ints_sm <- cpp_reshaper(ResTabC[['small_results']][[1]],ResTabC[['small_results']][[2]],
                               ResTabC[['small_results']][[3]], restab_ints_sm)
restab_ints_bg <- cpp_reshaper(ResTabC[['big_results']][[1]],ResTabC[['big_results']][[2]],
                               ResTabC[['big_results']][[3]], restab_ints_bg)
#use factored dataframe
restab_sm[,ncol(restab_sm)] <- restab_ints_sm[,ncol(restab_ints_sm)]
restab_bg[,ncol(restab_bg)] <- restab_ints_bg[,ncol(restab_ints_bg)]

saveRDS(restab_sm,file=paste0("~/MITUS/inst/",loc, "/sm_restab2.rds"))
saveRDS(restab_bg,file=paste0("~/MITUS/inst/",loc, "/bg_restab2.rds"))
}
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
