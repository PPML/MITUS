risk_pop_setup<-function(){
  HIV<-c(1.1,24.9,2.19,1)
  Diab<-c(30.3,3,1.89,2.12)
  Silica<-c(1,30,5,1)
  CKD<-c(20,15,12.8,1) #update these
  Child5yr<-c(19.9,2,1.2,1.3)
  PWID<-c(6.6,3,2,1.5)
  Immunosup<-c(1,10,3,1)
  abCXR<-c(.15,2,1,2)
  Contact<-c(.2,1,1,8)
  Immigrants<-c(21,1,1,10)
  Healthcare<-c(18,1,1,5)
  Congregate<-c(4,1,1,8)
################################################################################
####################              HIV POPULATION            ####################
################################################################################
HIV_ttt_vec<-def_ttt()
HIV_ttt_vec[[3]]<-HIV[1]
HIV_ttt_vec[[4]]<-1
HIV_ttt_vec[[5]]<-2020
HIV_ttt_vec[[6]]<-2025
HIV_ttt_vec[[7]]<-HIV[2]
HIV_ttt_vec[[8]]<-HIV[3]
HIV_ttt_vec[[9]]<-HIV[4]
################################################################################
####################              Diab POPULATION            ####################
################################################################################
Diab_ttt_vec<-def_ttt()
Diab_ttt_vec[[3]]<-Diab[1]
Diab_ttt_vec[[4]]<-1
Diab_ttt_vec[[5]]<-2020
Diab_ttt_vec[[6]]<-2025
Diab_ttt_vec[[7]]<-Diab[2]
Diab_ttt_vec[[8]]<-Diab[3]
Diab_ttt_vec[[9]]<-Diab[4]
################################################################################
####################              Silicosis            ####################
################################################################################
Silica_ttt_vec<-def_ttt()
Silica_ttt_vec[[3]]<-Silica[1]
Silica_ttt_vec[[4]]<-1
Silica_ttt_vec[[5]]<-2020
Silica_ttt_vec[[6]]<-2025
Silica_ttt_vec[[7]]<-Silica[2]
Silica_ttt_vec[[8]]<-Silica[3]
Silica_ttt_vec[[9]]<-Silica[4]
################################################################################
####################              CKD            ####################
################################################################################
CKD_ttt_vec<-def_ttt()
CKD_ttt_vec[[3]]<-CKD[1]
CKD_ttt_vec[[4]]<-1
CKD_ttt_vec[[5]]<-2020
CKD_ttt_vec[[6]]<-2025
CKD_ttt_vec[[7]]<-CKD[2]
CKD_ttt_vec[[8]]<-CKD[3]
CKD_ttt_vec[[9]]<-CKD[4]
################################################################################
####################              Child5yr            ####################
################################################################################
Child5yr_ttt_vec<-def_ttt()
Child5yr_ttt_vec[[3]]<-Child5yr[1]
Child5yr_ttt_vec[[4]]<-1
Child5yr_ttt_vec[[5]]<-2020
Child5yr_ttt_vec[[6]]<-2025
Child5yr_ttt_vec[[7]]<-Child5yr[2]
Child5yr_ttt_vec[[8]]<-Child5yr[3]
Child5yr_ttt_vec[[9]]<-Child5yr[4]

################################################################################
####################              PWID            ####################
################################################################################
PWID_ttt_vec<-def_ttt()
PWID_ttt_vec[[3]]<-PWID[1]
PWID_ttt_vec[[4]]<-1
PWID_ttt_vec[[5]]<-2020
PWID_ttt_vec[[6]]<-2025
PWID_ttt_vec[[7]]<-PWID[2]
PWID_ttt_vec[[8]]<-PWID[3]
PWID_ttt_vec[[9]]<-PWID[4]

################################################################################
####################              Immunosup            ####################
################################################################################
Immunosup_ttt_vec<-def_ttt()
Immunosup_ttt_vec[[3]]<-Immunosup[1]
Immunosup_ttt_vec[[4]]<-1
Immunosup_ttt_vec[[5]]<-2020
Immunosup_ttt_vec[[6]]<-2025
Immunosup_ttt_vec[[7]]<-Immunosup[2]
Immunosup_ttt_vec[[8]]<-Immunosup[3]
Immunosup_ttt_vec[[9]]<-Immunosup[4]

################################################################################
####################              abCXR            ####################
################################################################################
abCXR_ttt_vec<-def_ttt()
abCXR_ttt_vec[[3]]<-abCXR[1]
abCXR_ttt_vec[[4]]<-1
abCXR_ttt_vec[[5]]<-2020
abCXR_ttt_vec[[6]]<-2025
abCXR_ttt_vec[[7]]<-abCXR[2]
abCXR_ttt_vec[[8]]<-abCXR[3]
abCXR_ttt_vec[[9]]<-abCXR[4]

################################################################################
####################              Contact            ####################
################################################################################
Contact_ttt_vec<-def_ttt()
Contact_ttt_vec[[3]]<-Contact[1]
Contact_ttt_vec[[4]]<-1
Contact_ttt_vec[[5]]<-2020
Contact_ttt_vec[[6]]<-2025
Contact_ttt_vec[[7]]<-Contact[2]
Contact_ttt_vec[[8]]<-Contact[3]
Contact_ttt_vec[[9]]<-Contact[4]

################################################################################
####################              Immigrants            ####################
################################################################################
Immigrants_ttt_vec<-def_ttt()
Immigrants_ttt_vec[[3]]<-Immigrants[1]
Immigrants_ttt_vec[[4]]<-1
Immigrants_ttt_vec[[5]]<-2020
Immigrants_ttt_vec[[6]]<-2025
Immigrants_ttt_vec[[7]]<-Immigrants[2]
Immigrants_ttt_vec[[8]]<-Immigrants[3]
Immigrants_ttt_vec[[9]]<-Immigrants[4]

################################################################################
####################              Healthcare            ####################
################################################################################
Healthcare_ttt_vec<-def_ttt()
Healthcare_ttt_vec[[3]]<-Healthcare[1]
Healthcare_ttt_vec[[4]]<-1
Healthcare_ttt_vec[[5]]<-2020
Healthcare_ttt_vec[[6]]<-2025
Healthcare_ttt_vec[[7]]<-Healthcare[2]
Healthcare_ttt_vec[[8]]<-Healthcare[3]
Healthcare_ttt_vec[[9]]<-Healthcare[4]

################################################################################
####################              Congregate            ####################
################################################################################
Congregate_ttt_vec<-def_ttt()
Congregate_ttt_vec[[3]]<-Congregate[1]
Congregate_ttt_vec[[4]]<-1
Congregate_ttt_vec[[5]]<-2020
Congregate_ttt_vec[[6]]<-2025
Congregate_ttt_vec[[7]]<-Congregate[2]
Congregate_ttt_vec[[8]]<-Congregate[3]
Congregate_ttt_vec[[9]]<-Congregate[4]

all_pop_ttt<-list(HIV_ttt_vec,
                  Diab_ttt_vec,
                  Silica_ttt_vec,
                  CKD_ttt_vec,
                  Child5yr_ttt_vec,
                  PWID_ttt_vec,
                  Immunosup_ttt_vec,
                  abCXR_ttt_vec,
                  Contact_ttt_vec,
                  Immigrants_ttt_vec,
                  Healthcare_ttt_vec,
                  Congregate_ttt_vec)
names(all_pop_ttt)<-c("HIV", "Diabetes","Silicosis","CKD","Child5yr",
                                         "PWID","ImmunosupTherapy","abCXR","Contact",
                                         "Immigrants","Healthcare","Congregate")
return(all_pop_ttt)
}
