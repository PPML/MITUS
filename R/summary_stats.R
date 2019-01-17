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

  bc_results<-paste0(loc,"_restab_b")
  data(list=bc_results, package = 'MITUS')
  restab<-get(paste0(loc,"_restab_b"))
  bc_results<-filter(restab, scenario=="base_case" & year=="2018" & comparator=="absolute_value" & statistic=="mean")
  # bc_results[,"age_group"]<-as.character(bc_results[,"age_group"])
  bc_stats<-rep(NA,3)
  names(bc_stats)<-c("tb_incid","ltbi_prev","pop_size")
#incidence percentage
  bc_result_incid<-filter(bc_results, outcome=="tb_incidence_per_mil"& population==nat_grp)

  x<-switch(age_grp,
                               "All Ages"= {filter(bc_result_incid, age_group=="all")["value"]},
                               "0 to 24"={filter(bc_result_incid, age_group=="0-24")["value"]},
                               "25 to 64"={filter(bc_result_incid, age_group=="25-64")["value"]},
                               "65+"={filter(bc_result_incid, age_group=="65+")["value"]}
  )
  bc_stats["tb_incid"]<-unlist(x)/1e5
#ltbi prevalence percentage
  bc_result_ltbi<-filter(bc_results, outcome=="pct_ltbi"& population==nat_grp)
  bc_stats["ltbi_prev"]<-switch(age_grp,
                                "All Ages"={filter(bc_result_ltbi, age_group=="all")["value"]},
                                "0 to 24"={filter(bc_result_ltbi, age_group=="0-24")["value"]},
                                "25 to 64"={filter(bc_result_ltbi, age_group=="25-64")["value"]},
                                "65+"={filter(bc_result_ltbi, age_group=="65+")["value"]}
  )
#population size
  bc_result_nltbi<-filter(bc_results, outcome=="ltbi_000s" & population==nat_grp)
  bc_result_nltbi[,"value"]<-  bc_result_nltbi[,"value"]*1000
  bc_result_pop<-(bc_result_nltbi[,"value"])/bc_result_ltbi[,"value"]
  names(bc_result_pop)<-levels(bc_result_ltbi[,"age_group"])
  bc_stats["pop_size"]<-switch(age_grp,
                              "All Ages"={sum(bc_result_pop[])},
                              "0 to 24"={sum(bc_result_pop[1])},
                              "25 to 64"={sum(bc_result_pop[2])},
                              "65+"={sum(bc_result_pop[3])}
  )
  return(bc_stats)
  }
