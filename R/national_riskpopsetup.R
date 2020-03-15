national_risk_pop_setup<-function(){
  #population size, RRprogression, RRmu, RRLTBI
  HIV<-c(1.1,10,2.3,1)
  Diab<-c(30.3,1.5,1.9,1)
  Silica<-c(1800/1e6,1.26,5,1)#values from linas paper, mort from ?
  ESRD<-c(20,18,15.8,1)#match tabby2
  Child5yr<-c(19.9,2,1.2,1.3)
  PWID<-c(6.6,1,3,8.2)#values from linas paper, mort from https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002964
  Immunosup<-c(8.9,2,1,1.89) #value from linas paper, mort from Overall and cancer related mortality among patients with ocular inflammation treated with immunosuppressive drugs: retrospective cohort study
  Prisoners<-c(2.3, 1, 1, 1)
  Homeless<-c(96100/1e6, 1, 3.4, 5.3)
  load(system.file("US/US_results_1.rda", package="MITUS"))
  results<-out[1,,]
  Immigrants<-c(sum(results[71,635:682]),1,1,1) #match tabby2
  Healthcare<-c(18,1,1,5)
  # Congregate<-c(4,1,3.4,5.3) #values are for homeless rn
  ################################################################################
  ####################              HIV POPULATION            ####################
  ################################################################################
  HIV_ttt_vec<-def_ttt_nat_ag()
  HIV_ttt_vec[[3]]<-HIV[1]
  HIV_ttt_vec[[4]]<-1
  HIV_ttt_vec[[5]]<-2020
  HIV_ttt_vec[[6]]<-2020
  HIV_ttt_vec[[7]]<-HIV[2]
  HIV_ttt_vec[[8]]<-HIV[3]
  HIV_ttt_vec[[9]]<-HIV[4]
  ################################################################################
  ####################              Diab POPULATION            ####################
  ################################################################################
  Diab_ttt_vec<-def_ttt_nat_ag()
  Diab_ttt_vec[[3]]<-Diab[1]
  Diab_ttt_vec[[4]]<-1
  Diab_ttt_vec[[5]]<-2020
  Diab_ttt_vec[[6]]<-2020
  Diab_ttt_vec[[7]]<-Diab[2]
  Diab_ttt_vec[[8]]<-Diab[3]
  Diab_ttt_vec[[9]]<-Diab[4]
  ################################################################################
  ####################              Silicosis            ####################
  ################################################################################
  Silica_ttt_vec<-def_ttt_nat_ag()
  Silica_ttt_vec[[3]]<-Silica[1]
  Silica_ttt_vec[[4]]<-1
  Silica_ttt_vec[[5]]<-2020
  Silica_ttt_vec[[6]]<-2020
  Silica_ttt_vec[[7]]<-Silica[2]
  Silica_ttt_vec[[8]]<-Silica[3]
  Silica_ttt_vec[[9]]<-Silica[4]
  ################################################################################
  ####################              ESRD            ####################
  ################################################################################
  ESRD_ttt_vec<-def_ttt_nat_ag()
  ESRD_ttt_vec[[3]]<-ESRD[1]
  ESRD_ttt_vec[[4]]<-1
  ESRD_ttt_vec[[5]]<-2020
  ESRD_ttt_vec[[6]]<-2020
  ESRD_ttt_vec[[7]]<-ESRD[2]
  ESRD_ttt_vec[[8]]<-ESRD[3]
  ESRD_ttt_vec[[9]]<-ESRD[4]
  ################################################################################
  ####################              Child5yr            ####################
  ################################################################################
  Child5yr_ttt_vec<-def_ttt_nat_ag()
  Child5yr_ttt_vec[[1]]<- c(1, rep(0,10))
  Child5yr_ttt_vec[[2]]<- c(1, rep(0,10))
  Child5yr_ttt_vec[[3]]<-Child5yr[1]
  Child5yr_ttt_vec[[4]]<-1
  Child5yr_ttt_vec[[5]]<-2020
  Child5yr_ttt_vec[[6]]<-2020
  Child5yr_ttt_vec[[7]]<-Child5yr[2]
  Child5yr_ttt_vec[[8]]<-Child5yr[3]
  Child5yr_ttt_vec[[9]]<-Child5yr[4]

  ################################################################################
  ####################              PWID            ####################
  ################################################################################
  PWID_ttt_vec<-def_ttt_nat_ag()
  PWID_ttt_vec[[3]]<-PWID[1]
  PWID_ttt_vec[[4]]<-1
  PWID_ttt_vec[[5]]<-2020
  PWID_ttt_vec[[6]]<-2020
  PWID_ttt_vec[[7]]<-PWID[2]
  PWID_ttt_vec[[8]]<-PWID[3]
  PWID_ttt_vec[[9]]<-PWID[4]

  ################################################################################
  ####################              Immunosup            ####################
  ################################################################################
  Immunosup_ttt_vec<-def_ttt_nat_ag()
  Immunosup_ttt_vec[[3]]<-Immunosup[1]
  Immunosup_ttt_vec[[4]]<-1
  Immunosup_ttt_vec[[5]]<-2020
  Immunosup_ttt_vec[[6]]<-2020
  Immunosup_ttt_vec[[7]]<-Immunosup[2]
  Immunosup_ttt_vec[[8]]<-Immunosup[3]
  Immunosup_ttt_vec[[9]]<-Immunosup[4]

  ################################################################################
  ####################              Prisoners            ####################
  ################################################################################
  Prisoners_ttt_vec<-def_ttt_nat_ag()
  Prisoners_ttt_vec[[3]]<-Prisoners[1]
  Prisoners_ttt_vec[[4]]<-1
  Prisoners_ttt_vec[[5]]<-2020
  Prisoners_ttt_vec[[6]]<-2020
  Prisoners_ttt_vec[[7]]<-Prisoners[2]
  Prisoners_ttt_vec[[8]]<-Prisoners[3]
  Prisoners_ttt_vec[[9]]<-Prisoners[4]

  ################################################################################
  ####################              Homeless            ####################
  ################################################################################
  Homeless_ttt_vec<-def_ttt_nat_ag()
  Homeless_ttt_vec[[3]]<-Homeless[1]
  Homeless_ttt_vec[[4]]<-1
  Homeless_ttt_vec[[5]]<-2020
  Homeless_ttt_vec[[6]]<-2020
  Homeless_ttt_vec[[7]]<-Homeless[2]
  Homeless_ttt_vec[[8]]<-Homeless[3]
  Homeless_ttt_vec[[9]]<-Homeless[4]

  ################################################################################
  ####################              Immigrants            ####################
  ################################################################################
  Immigrants_ttt_vec<-def_ttt_nat_ag()
  Immigrants_ttt_vec[[1]]<-rep(0,11)
  Immigrants_ttt_vec[[3]]<-Immigrants[1]
  Immigrants_ttt_vec[[4]]<-1
  Immigrants_ttt_vec[[5]]<-2020
  Immigrants_ttt_vec[[6]]<-2020
  Immigrants_ttt_vec[[7]]<-Immigrants[2]
  Immigrants_ttt_vec[[8]]<-Immigrants[3]
  Immigrants_ttt_vec[[9]]<-Immigrants[4]

  ################################################################################
  ####################              Healthcare            ####################
  ################################################################################
  Healthcare_ttt_vec<-def_ttt_nat_ag()
  Healthcare_ttt_vec[[3]]<-Healthcare[1]
  Healthcare_ttt_vec[[4]]<-1
  Healthcare_ttt_vec[[5]]<-2020
  Healthcare_ttt_vec[[6]]<-2020
  Healthcare_ttt_vec[[7]]<-Healthcare[2]
  Healthcare_ttt_vec[[8]]<-Healthcare[3]
  Healthcare_ttt_vec[[9]]<-Healthcare[4]

  ################################################################################
  ####################              Congregate            ####################
  ################################################################################
  # Congregate_ttt_vec<-def_ttt_nat_ag()
  # Congregate_ttt_vec[[3]]<-Congregate[1]
  # Congregate_ttt_vec[[4]]<-1
  # Congregate_ttt_vec[[5]]<-2020
  # Congregate_ttt_vec[[6]]<-2020
  # Congregate_ttt_vec[[7]]<-Congregate[2]
  # Congregate_ttt_vec[[8]]<-Congregate[3]
  # Congregate_ttt_vec[[9]]<-Congregate[4]

  all_pop_ttt<-list(HIV_ttt_vec,
                    Diab_ttt_vec,
                    Silica_ttt_vec,
                    ESRD_ttt_vec,
                    Child5yr_ttt_vec,
                    PWID_ttt_vec,
                    Immunosup_ttt_vec,
                    Prisoners_ttt_vec,
                    Homeless_ttt_vec,
                    Immigrants_ttt_vec,
                    Healthcare_ttt_vec
                   )
  names(all_pop_ttt)<-c("HIV", "Diabetes","Silicosis","ESRD","Child5yr",
                        "PWID","ImmunosupTherapy","Prisoners","Homeless",
                        "Immigrants","Healthcare")
  return(all_pop_ttt)
}
