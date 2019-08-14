#' THIS FUNCTION INPUTS A TABLE OF PARAMETERS AND RUNS THE TB MODEL
#' AND GENERATES AN ARRAY OF OUTPUTS FOR EACH POSSIBLE PREDEFINED
#' SCENARIO; THIS ARRAY IS THEN SAVED TO THE DATA FOLDER FOR THAT
#' LOCATION

#'@name reshape_results
#'@param loc two character location code for the model location

reshape_results<-function(loc){
#load tabus
library(tabus)

#load mitus sims into data structure
load_data <- function(i) {
  data_name <-
    load(system.file(paste0(loc, "/", loc, "_results_",i,".rda"), package='MITUS'))
  return(get(data_name))
}

results <- lapply(1:9, load_data)

#make some lists of reformatted data
ResTabC <- list()
ResTabC[['small_results']] <- format_as_restab_small_ages_indices(results)
ResTabC[['big_results']] <- format_as_restab_big_ages_indices(results)

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
