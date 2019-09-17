#' THIS FUNCTION READS IN THE BASE CASE SCENARIO AND THEN
#' CREATES INDIVIDUAL RDS FILES FOR OUTPUTS NECESSARY FOR
#' THE TABBY2 COMPARISON TO RECENT DATA PLOTS

#'@name make_all_scenarios
#'@param loc two character location code for the model location
#'@param bc.array base case results array
#'@param samp.i which simulation to use (1:10)
#'@param simp.date date to append to the files; should match the optim date
#'@export

model_calib_outputs<-function(loc="US",bc.array, samp_i=1,simp.date){
  res<-as.data.frame(bc.array[samp_i,,])
  colnames(res)<-func_ResNam()
  ############                   demographic targets                     ############
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  V  <- cbind(res[1:68,30], res[1:68,31]+res[1:68,32])
  saveRDS(V,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_pop_yr_nat_",simp.date,".rds"))
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  V  <- cbind(t(res[68,33:43]), t(res[68,44:54]))
  V1  <- V[-3,]
  V1[2,] <- V1[2,]+V[3,]
  V2 <- V1[-4,]
  V2[3,] <- V2[3,]+V1[4,]
  V3 <- V2[-9,]
  V3[8,] <- V3[8,]+V2[9,]
  saveRDS(V3,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_pop_ag_nat_",simp.date,".rds"))

  ### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
  V  <- cbind(rowSums(res[1:67,255:265]), rowSums(res[1:67,266:276]))
  V1c <- rowSums(res[1:67,121:131])
  saveRDS(V1c,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_mort_yr_nat_",simp.date,".rds"))

  ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  V  <- cbind((res[67,255:265])+(res[67,266:276]))
  V1  <- V[,-3]
  V1[,2] <- V1[,2]+V[,3]
  V2 <- V1[,-4]
  V2[,3] <- V2[,3]+V1[,4]
  V3 <- V2[,-9]
  V3[,8] <- V3[,8]+V2[,9]
  saveRDS(V3,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_mort_ag_nat_",simp.date,".rds"))

  ############                   tb specific targets                     ############
  # graph of total diagnosed cases
  # by total population, US born population, and non-US born population
  V0 <- res[4:68,"NOTIF_ALL"]+res[4:68,"NOTIF_MORT_ALL"] #total population
  V1 <- res[44:68,"NOTIF_US"]+res[44:68,"NOTIF_MORT_US"]   #US born population
  V2 <- res[44:68,"NOTIF_F1"]+res[44:68,"NOTIF_F2"]+res[44:68,"NOTIF_MORT_F1"]+res[44:68,"NOTIF_MORT_F2"]   #non-US born population
  tot_cases<-list()
  tot_cases[["allpop"]]<-V0
  tot_cases[["USBpop"]]<-V1
  tot_cases[["NUSBpop"]]<-V2
  saveRDS(tot_cases,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_TBcases_",simp.date,".rds"))

  #Percent of Total Cases Non-US Born Population
  V <- cbind(res[44:67,"NOTIF_US"]+res[44:67,"NOTIF_MORT_US"], #US born population
             res[44:67,"NOTIF_F1"]+res[44:67,"NOTIF_F2"]+  #non-US born population
               res[44:67,"NOTIF_MORT_F1"]+res[44:67,"NOTIF_MORT_F2"])
  V <- V[,2]/rowSums(V)
  #Percent of Non-US Born Cases from Recent Immigrant Population
  V <- cbind(res[44:65,"NOTIF_F1"]+res[44:65,"NOTIF_MORT_F1"],res[44:65,"NOTIF_F2"]+res[44:65,"NOTIF_MORT_F2"])
  V <- V[,1]/rowSums(V)*100
  saveRDS(V,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_percentRecentFBcases_",simp.date,".rds"))

  #Age distribution of Cases
  #0-24 yrs, 25-44 yrs, 45-64 yrs, 65+ yrs
  #ends in 2016
  V   <- (res[51:67,136:146]+res[51:67,189:199])
  V2<-matrix(NA,length(51:67),4)
  V2[,1]<-rowSums(V[,1:3]); V2[,2]<-rowSums(V[,4:5])
  V2[,3]<-rowSums(V[,6:8]); V2[,4]<-rowSums(V[,9:11])
  saveRDS(V2*1e6,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_age_cases_4grps_",simp.date,".rds"))
  #Age distribution of Cases
  #all age bands
  V   <- (res[51:67,136:146]+res[51:67,189:199])
  V2  <- cbind(2000:2016,V)
  saveRDS(V2,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_age_cases_tot_",simp.date,".rds"))

  # Treatment Outcomes 1993-2014
  V   <- res[44:65,132:134]
  Vdisc <- V[,2]/rowSums(V)
  Vdead <- V[,3]/rowSums(V)

  txoutcomes<-list()
  txoutcomes[["discontinued tx"]]<-Vdisc
  txoutcomes[["died on tx"]]<-Vdead

  saveRDS(txoutcomes,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_txOutcomes_",simp.date,".rds"))


  #LTBI Prevalance by Age in 2011, US born
  V  <- cbind(t(res[62,55:65]),t(res[62,33:43]-res[62,55:65]))
  pIGRA<-1
  v1<-V*pIGRA
  Sens_IGRA <-c(.780,.675,.712,.789,.591)
  Spec_IGRA <-c(.979,.958,.989,.985,.931)
  names(Sens_IGRA)<- names(Spec_IGRA)<-c("lrUS","hrUS","youngNUS","NUS","hrNUS")
  Va <- outer(v1[,1],c(Sens_IGRA[1],(1-Sens_IGRA[1])))+outer(v1[,2],c((1-Spec_IGRA[1]),Spec_IGRA[1]))
  V1 <- Va[-11,]; V1<-V1[-10,]
  V1[9,] <- V1[9,]+Va[10,]+Va[11,]
  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100
  # colnames(V2) <- c("LTBI", "No-LTBI")

  saveRDS(V2,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_USB_LTBI_pct_",simp.date,".rds"))

  #LTBI Prevalance by Age in 2011, non-US born
  V  <- cbind(t(res[62,66:76]),t(res[62,44:54]-res[62,66:76]))
  pIGRA<-1
  v1<-V*pIGRA
  #under age 5
  v1b <- (v1[1,1]*c(Sens_IGRA[3],(1-Sens_IGRA[3])))+(v1[1,2]*c((1-Spec_IGRA[3]),Spec_IGRA[3]))
  #over age 5
  v1c <- outer(v1[2:11,1],c(Sens_IGRA[4],(1-Sens_IGRA[4])))+outer(v1[2:11,2],c((1-Spec_IGRA[4]),Spec_IGRA[4]))
  v1d<-rbind(v1b,v1c)
  colnames(v1d) <- c("LTBI", "No-LTBI")
  V1 <- v1d[-11,]; V1<-V1[-10,]
  V1[9,] <- V[9,]+v1d[10,]+v1d[11,]
  V2 <- rep(NA,8)
  V2 <- V1[2:9,1]/rowSums(V1[2:9,])*100
  saveRDS(V2,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_NUSB_LTBI_pct_",simp.date,".rds"))

  # Age Distribution of TB Deaths 1999-2014

  V  <- res[50:65,227:237]
  V2 <- V[,-11]; V2[,10] <- V[,10]+V[,11]
  V3 <- colSums(V2)*1e6

  saveRDS(V3,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_TBdeathsAge_",simp.date,".rds"))

  # total tb deaths over time 2004-2014
  V   <- rowSums(res[55:65,227:237])
  saveRDS(V,file = paste0("~/MITUS/inst/",loc,"/calibration_outputs/",loc,"_TBdeaths_",simp.date,".rds"))


}
