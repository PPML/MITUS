# tabby_results <- function(loc,startyr=1950,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0) {
#
#   wt = parAll200[,"wt"]
#   ## HELPER FUNCTIONS
#   # l takes a vector of results, and weights, and calculates weighted mean, lower bound, and upper bound for q=1,2,3 respectively
#   l <- function(x,w,q) {
#     if(q==1) {
#       wtd.mean(x,w)
#     } else {
#       if(q==2) {
#         wtd.quantile(x,w,1/40)
#       } else {
#         wtd.quantile(x,w,39/40)
#       }
#     }
#   }
#
#   ## PULL OUT SPECIFIC RESULTS (5) FOR EACH SCENARIO (8), FOR DIFFERENT POPS (3)
#   # OUTCOMES = 1. Incident M. tb Infection (per Mil)
#   #           2. LTBI Prevalence (%)
#   #           3. New TB Cases (per Mil)
#   #           4. MDR-TB in New TB Cases (%)
#   #           5. TB-Related Deaths (per Mil)
#   #           Z. ADD MORE AS DESIRED!
#
#   # SCENARIOS = 1. Base case
#   #             2. Intervention scenario 1 --TLTBI for new immigrants
#   #             3. Intervention scenario 2 --Improved TLTBI in US
#   #             4. Intervention scenario 3 --Better case detection
#   #             5. Intervention scenario 4 --Better TB treatment
#   #             6. Intervention scenario 5 --All improvements together
#   #             7. Sensitivity analysis  1 --no transmission within the United States after 2016
#   #             8. Sensitivity analysis  2 --no LTBI among all immigrants after 2016
#
#   # POPULATIONS = 1. ALL
#   #               2. US-BORN
#   #               3. FOREIGN-BORN
#
#
#
#   age_id = (2016:endyr)-1949
#   #create 3 lists to hold output
#   ResTabAllfb <- ResTabAllus <- ResTabAll <- list()
#
#
#   #concatenate all lists
#   ResTabC <- list(ResTabAll,ResTabAllus,ResTabAllfb)
#
#
#
#
#
#   ### FOR EACH OUTCOME, EST MEAN AND POSTERIOR INTERVAL FOR SPECIFIC YEARS (2025, 2050, 2100), PLUS PERCENT REDUCTION COMPARED TO BASE CASE
#   ### long format
#   #1 outcome  5
#   #2 scenario  8
#   #3 population 3
#   #4 age 4  **** nothing so far
#   #5 comparator -- 3 (abs, 2016 BC, current BC)
#   #6 time points 4
#   #7 statistic
#
#   CatList <- list()
#   CatList[[1]] <- c("ltbi_000s","pct_ltbi","tb_incidence_000s","tb_incidence_per_mil","tb_mortality_000s","tb_mortality_per_mil")
#   CatList[[2]] <- c("base_case",paste("intervention_",1:5,sep=""),paste("scenario_",1:7,sep=""))[wch+1]
#   CatList[[3]] <- c("all_populations","usb_population","fb_population")
#   CatList[[4]] <- c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
#   CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
#   CatList[[6]] <- 2016:2099
#   CatList[[7]] <- c("mean","ci_low","ci_high")
#
#   res_tab2 <-  cbind(expand.grid(CatList),NA)
#   colnames(res_tab2) <- c("outcome","scenario","population","age_group","comparator","year","statistic","value")
#
#   for(i1 in 1:length(CatList[[1]])) { # i1=i3=i4=i5=i6=i7=1
#     for(i3 in 1:length(CatList[[3]])) { #
#       for(i4 in 1:length(CatList[[4]])) { #
#         for(i6 in 1:length(CatList[[6]])) { #
#           for(i7 in 1:length(CatList[[7]])) { #
#              rw = res_tab2[,1]==CatList[[1]][i1] & res_tab2[,2]==CatList[[2]][1] & res_tab2[,3]==CatList[[3]][i3] &
#              res_tab2[,4]==CatList[[4]][i4] &  res_tab2[,6]==CatList[[6]][i6] & res_tab2[,7]==CatList[[7]][i7]
#              res_tab2[rw & res_tab2[,5]==CatList[[5]][1],"value"] = l(ResTabC[[i3]][[2]][,i6,i1+1,i4],wt,i7)
#              res_tab2[rw & res_tab2[,5]==CatList[[5]][2],"value"] = l(ResTabC[[i3]][[2]][,i6,i1+1,i4]/
#                                                              ResTabC[[i3]][[1 ]][,i6,i1+1,i4],wt,i7)*100
#              res_tab2[rw & res_tab2[,5]==CatList[[5]][3],"value"] = l(ResTabC[[i3]][[2]][,i6,i1+1,i4]/
#                                                              ResTabC[[i3]][[1 ]][,i6,i1+1,i4],wt,i7)*100
# } } } } }
#
#   }
