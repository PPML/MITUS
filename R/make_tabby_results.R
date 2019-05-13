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
results<-readRDS(system.file(paste0(loc,"/sm_resTab1_2019-05-13.rds"), package="MITUS"))
#make an empty restab
sm_restab<-make_empty_res_tab2sm()
sm_restab %<>% mutate_if(is.factor, as.integer) %>% as.matrix
#reshape it baby
res_tab2<-cpp_reshaper(results[[1]],results[[2]],results[[3]],sm_restab)
# Specify the levels of each dimension to the data
CatList <- list()
CatList[[1]] <- c(
  "ltbi_000s",
  "pct_ltbi",
  "tb_incidence_000s",
  "tb_incidence_per_mil",
  "tb_mortality_000s",
  "tb_mortality_per_mil")
CatList[[2]] <- c("base_case",paste("intervention_",1:5,sep=""),paste("scenario_",1:3,sep=""))
CatList[[3]] <- c("all_populations","usb_population","fb_population")
CatList[[4]] <- c("0-4",paste(0:8*10+5,1:9*10+4,sep="-"),"95+")
CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
CatList[[6]] <- 2018:2049

# Re-Factorize each column
for (i in 1:6) {
  res_tab2[,i] <- factor(res_tab2[,i], labels = CatList[[i]])
}
#save the results
saveRDS(res_tab2,file=paste0("~/MITUS/inst/", loc, "/sm_restab2.rds"))

#####################################################################################
rm(results)
results<-list()
#load in the big results
results<-readRDS(system.file(paste0(loc,"/bg_resTab1_2019-05-13.rds"), package="MITUS"))
#make an empty restab
bg_restab<-make_empty_res_tab2bg()
bg_restab %<>% mutate_if(is.factor, as.integer) %>% as.matrix
#reshape it baby
res_tab2<-cpp_reshaper(results[[1]],results[[2]],results[[3]],bg_restab)
# Specify the levels of each dimension to the data
CatList <- list()
CatList[[1]] <- c(
  "pct_ltbi",
  "tb_infection_per_mil",
  "tb_incidence_per_mil",
  "tb_deaths_per_mil"
)
CatList[[2]] <- c("base_case",paste("intervention_",1:5,sep=""),paste("scenario_",1:3,sep=""))
CatList[[3]] <- c("all_populations","usb_population","fb_population")
CatList[[4]] <- c("all_ages", "age_0_24","age_25_64","age_65p")
CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
CatList[[6]] <- 2018:2049

# Re-Factorize each column
for (i in 1:6) {
  res_tab2[,i] <- factor(res_tab2[,i], labels = CatList[[i]])
}
#save the results
saveRDS(res_tab2,file=paste0("~/MITUS/inst/", loc, "/bg_restab2.rds"))

}
