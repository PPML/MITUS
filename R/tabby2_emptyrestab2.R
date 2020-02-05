# make_empty_res_tab2sm <- function() {
#   # Specify the levels of each dimension to the data
#   CatList <- list()
#   CatList[[1]] <- c(
#     "ltbi_000s",
#     "pct_ltbi",
#     "tb_incidence_000s",
#     "tb_incidence_per_mil",
#     "tb_mortality_000s",
#     "tb_mortality_per_mil")
#   CatList[[2]] <- c('base_case', paste0('intervention_', 1:5), paste0('scenario_', 1:3))
#   CatList[[3]] <- c("all_populations","usb_population","fb_population")
#   CatList[[4]] <- c("0-4",paste(0:8*10+5,1:9*10+4,sep="-"),"95+")
#   CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
#   CatList[[6]] <- 2018:2049
#   # CatList[[7]] <- c("mean","ci_low","ci_high")
#
#   # Make the specified levels integer-leveled factors
#   CatList_factors <- lapply(CatList, function(x) {
#     factor(x=1:length(x), levels=1:length(x), labels=x) })
#
#   # Turn the integer-factors into all possible combinations in a dataframe
#   # with an extra column of NA values for a 'value' column
#   res_tab2 <-  cbind(expand.grid(CatList_factors),NA)
#
#   # Name the columns
#   colnames(res_tab2) <- c(
#     "outcome",
#     "scenario",
#     "population",
#     "age_group",
#     "comparator",
#     "year",
#     "value")
#
#   return(res_tab2)
# }
#
#
# make_empty_res_tab2bg <- function() {
#   # Specify the levels of each dimension to the data
#   CatList <- list()
#   CatList[[1]] <- c(
#     "pct_ltbi",
#     "tb_infection_per_mil",
#     "tb_incidence_per_mil",
#     "tb_deaths_per_mil")
#   CatList[[2]] <- c('base_case', paste0('intervention_', 1:5), paste0('scenario_', 1:3))
#   CatList[[3]] <- c("all_populations","usb_population","fb_population")
#   CatList[[4]] <- c("all", "0-24","25-64","65+")
#   CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
#   CatList[[6]] <- 2018:2049
#   # CatList[[7]] <- c("mean","ci_low","ci_high")
#
#   # Make the specified levels integer-leveled factors
#   CatList_factors <- lapply(CatList, function(x) {
#     factor(x=1:length(x), levels=1:length(x), labels=x) })
#
#   # Turn the integer-factors into all possible combinations in a dataframe
#   # with an extra column of NA values for a 'value' column
#   res_tab2 <-  cbind(expand.grid(CatList_factors),NA)
#
#   # Name the columns
#   colnames(res_tab2) <- c(
#     "outcome",
#     "scenario",
#     "population",
#     "age_group",
#     "comparator",
#     "year",
#     "value")
#
#   return(res_tab2)
# }
