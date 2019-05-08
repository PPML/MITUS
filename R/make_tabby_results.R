#make results
make_results<-function(loc){
#load library inline
  library(inline)
  library(dplyr)
#####################################################################################

#define the cpp reshaper
  cpp_reshaper <- cxxfunction(
    signature(ResTab='numeric', ResTabus='numeric', ResTabfb='numeric', res_tab2 = 'numeric'),
    plugin='Rcpp',
    body=readr::read_file(
      system.file('inline_cpp/format_restab2.cpp', package='MITUS')))

#####################################################################################

#load in the little results
results<-readRDS(system.file(paste0(loc,"/sm_resTab1_2019-05-08.rds"), package="MITUS"))
#make an empty restab
sm_restab<-make_empty_res_tab2sm()
sm_restab %<>% mutate_if(is.factor, as.integer) %>% as.matrix

#reshape it baby
res_tab2<-cpp_reshaper(results[[1]],results[[2]],results[[3]],sm_restab)
#save the results
saveRDS(res_tab2,file=paste0("~/MITUS/inst/", loc, "/sm_restab2.rds"))

#####################################################################################
rm(results)
results<-list()
#load in the big results
results<-readRDS(system.file(paste0(loc,"/bg_resTab1_2019-05-08.rds"), package="MITUS"))
#make an empty restab
bg_restab<-make_empty_res_tab2bg()
bg_restab %<>% mutate_if(is.factor, as.integer) %>% as.matrix

#reshape it baby
res_tab2<-cpp_reshaper(results[[1]],results[[2]],results[[3]],bg_restab)
#save the results
saveRDS(res_tab2,file=paste0("~/MITUS/inst/", loc, "/bg_restab2.rds"))

}
