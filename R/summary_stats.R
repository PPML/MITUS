#'@name sum_stats
#'@description takes arguments from Tabby2 user interface and returns summary statistics for different age and nativity groups
#'@param loc two letter postal abbreviation for location
#'@param age_grp which age group of interest
#'@param nat_grp which nativity group of interest
#'@return summary statistics of
#'@export

sum_stat<-function(loc,age_grp,nat_grp){
  library(dplyr)

  if (nat_grp=="All Nativity Groups") {ng=1}
  if (nat_grp=="U.S. Born") {ng=2}
  if (nat_grp=="Non-U.S. Born") {ng=3}
  nat_grp <-switch(ng, "all_populations", "usb_population","fb_population")

  bc_results<-paste0(loc,"_restab2")
  data(list=bc_results, package = 'MITUS')
  bc_results<-filter(res_tab2, scenario=="base_case" & year=="2018" & comparator=="absolute_value" & statistic=="mean")
  # bc_results[,"age_group"]<-as.character(bc_results[,"age_group"])
  bc_stats<-rep(NA,3)
  names(bc_stats)<-c("tb_incid","ltbi_prev","pop_size")
#incidence percentage
  bc_result_incid<-filter(bc_results, outcome=="tb_incidence_per_mil"& population==nat_grp)

  bc_stats["tb_incid"]<-switch(age_grp,
                               "All Ages"={sum(bc_result_incid[,"value"])},
                               "0 to 24"={sum(filter(bc_result_incid, age_group %in% c("0_4","5_14","15_24"))["value"])},
                               "25 to 64"={sum(filter(bc_result_incid, age_group %in% c("25_34","35_44","45_54","55_64"))["value"])},
                               "65+"={sum(filter(bc_result_incid, age_group %in% c("65_74","75_84","95p"))["value"])}
  )
  bc_stats["tb_incid"]<-bc_stats["tb_incid"]/10e4
#ltbi prevalence percentage
  bc_result_ltbi<-filter(bc_results, outcome=="pct_ltbi"& population==nat_grp)
  bc_stats["ltbi_prev"]<-switch(age_grp,
                                "All Ages"={sum(bc_result_ltbi[,"value"])},
                                "0 to 24"={sum(filter(bc_result_ltbi, age_group %in% c("0_4","5_14","15_24"))["value"])},
                                "25 to 64"={sum(filter(bc_result_ltbi, age_group %in% c("25_34","35_44","45_54","55_64"))["value"])},
                                "65+"={sum(filter(bc_result_ltbi, age_group %in% c("65_74","75_84","95p"))["value"])}
  )
#population size
  #this is problematic -- update and check
  bc_result_nltbi<-filter(bc_results, outcome=="ltbi_000s" & population==nat_grp)
  bc_result_nltbi[,"value"]<-  bc_result_nltbi[,"value"]*1000
  bc_result_pop<-(bc_result_nltbi[,"value"])/bc_result_ltbi[,"value"]
  names(bc_result_pop)<-levels(bc_result_ltbi[,"age_group"])
  bc_stats["pop_size"]<-switch(age_grp,
                              "All Ages"={sum(bc_result_pop[])},
                              "0 to 24"={sum(bc_result_pop[1:3])},
                              "25 to 64"={sum(bc_result_pop[4:7])},
                              "65+"={sum(bc_result_pop[8:11])}
  )
  return(bc_stats)
  }
