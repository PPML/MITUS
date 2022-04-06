func_ResNam<-function(){
  ################################################################################
  ##### CREATE A LIST TO HOLD THE VECTORS FOR AGE CATEGORIES, TB STATES,     #####
  ##### DRUG RESISTANCE, TREATMENT HISTORY, HIV STATUS, AND RISK CATEGORY.   #####
  ################################################################################
  ################################################################################

  ################################################################################
  StatList <- noquote(list(
    #####					                  AGE CATEGORIES                             #####
    #####           AGE GROUPS 0-4, 5-14, 15-24, ... , 95-94, 95+              #####
    c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p"),
    #####					     				    TUBERCULOSIS STATES		              			   #####
    ##### SUSCEPTIBLE; UNINFECTED & PARTIALLY IMMUNE; LATENT SLOW; LATENT FAST;#####
    #####       ACTIVE TB SMEAR NEG.; ACTIVE TB SMEAR POS.; TB TREATMENT       #####
    c("Su","Sp","Ls","Lf", "Ac", "Tx"),
    #####                          TREATMENT HISTORY                           #####
    #####              NO TB TREATMENT HISTORY, PRIOR LTBI TREATMENT           #####
    c("NT","LT"),
    #####                        TB REACTIVATION RISK                          #####
    c("I1","I2","I3","I4"),
    #####                          NON-TB MORTALITY                            #####
    c("M1","M2","M3","M4"),
    #####                  LIVING CONDITIONS/RISK OF INFECTION                 #####
    c("LR","HR"),
    #####                              NATIVITY                                #####
    c("US","F1","F2")))
  ################################################################################
  ################################################################################
  ResNam <- c("Year",                                         # year
              ##############################    POPULATION    ################################
              "N_ALL",                                        # total pop
              paste("N",StatList[[1]],sep="_"),               # pop by ag cat
              paste("N",StatList[[2]],sep="_"),               # pop by tb cat
              paste("N",StatList[[4]],sep="_"),               # pop by im cat
              paste("N",StatList[[5]],sep="_"),               # pop by nm cat
              paste("N",StatList[[6]],sep="_"),               # pop by rg cat
              paste("N",StatList[[7]],sep="_"),               # pop by na cat
              ###################    POPULATION W/ NATIVITY  ################################
              paste("N_US",StatList[[1]],sep="_"),            # US pop by ag cat
              paste("N_FB",StatList[[1]],sep="_"),            # FB pop by ag cat
              paste("N_US_LTBI",StatList[[1]],sep="_"),       # US LTBI pop by ag cat
              paste("N_FB_LTBI",StatList[[1]],sep="_"),       # FB LTBI pop by ag cat
              paste("N_RF",StatList[[1]],sep="_"),            # RF pop by ag cat
              #############################    MORTALITY    ################################
              paste("TBMORT_US",StatList[[1]],sep="_"),   # TB mort, HIV neg, by ag cat
              paste("TBMORT_NUS",StatList[[1]],sep="_"),   # TB mort, HIV pos, by ag cat
              paste("RFMORT",StatList[[1]],sep="_"),         # RF mort, by ag cat
              paste("TOTMORT",StatList[[1]],sep="_"),         # total mort, by ag cat
              #### @object 131
              #############################   TX OUTCOMES   ################################

              "TBTX_COMPLT","TBTX_DISCONT","TBTX_DIED",       # TB treatment outcomes complete, discontinue, death
              "NOTIF_ALL",                                    # total notif
              paste("NOTIF",StatList[[1]],sep="_"),           # notif by ag cat
              paste("NOTIF",StatList[[7]],sep="_"),         # notif by nat cat
              #            paste("NOTIF_HIV",c("POS","NEG"),sep="_"),      # notif by HIV pos/neg
              paste("NOTIF",StatList[[6]],sep="_"),           # notif by rg cat
              #           paste("NOTIF_US_N",StatList[[3]],sep="_"),      # notif, US, N, by dr cat
              #            paste("NOTIF_US_E",StatList[[3]],sep="_"),      # notif, US, E, by dr cat
              #           paste("NOTIF_FB_N",StatList[[3]],sep="_"),      # notif, FB, N, by dr cat
              #           paste("NOTIF_FB_E",StatList[[3]],sep="_"),      # notif, FB, E, by dr cat
              ###########################   TLTBI INITIATION   ##############################
              "TLTBI_INITS",                                  # Initiations on LTBI tx
              "TLTBI_INITS_FB",                               # Initiations on LTBI tx FB
              "TLTBI_INITS_HR",                               # Initiations on LTBI tx HR
              #            "TLTBI_INITS_HV",                               # Initiations on LTBI tx HV
              "TLTBI_INITS_TP",                               # Initiations on LTBI tx, with LTBI
              ###########################     TB INCIDENCE     ##############################
              "INCID_ALL",                                    # Total incidence
              paste("INCID_ALL",StatList[[1]],sep="_"),       # Total incidence by ag cat
              "INCID_ALL_US",                                 # Total incidence, US born
              "INCID_ALL_FB",                                 # Total incidence, foreign born
              "INCID_ALL_FB2",                                # Total incidence, foreign born
              "INCID_ALL_HR",                                 # Total incidence, high risk
              #            "INCID_ALL_HV",                                 # Total incidence, HIV pos
              "INCID_REC",                                    # Total incidence, recent infection
              paste("INCID_REC",StatList[[1]],sep="_"),       # Total incidence by ag cat, recent infection
              "INCID_REC_US",                                 # Total incidence, US born, recent infection
              "INCID_REC_FB",                                 # Total incidence, foreign born, recent infection
              "INCID_REC_FB2",                                # Total incidence, foreign born, recent infection
              "INCID_REC_HR",                                 # Total incidence, high risk, recent infection
              #            "INCID_REC_HV",                                 # Total incidence, HIV pos, recent infection
              ###########################    NOTIFICATION DEAD      ##############################
              "NOTIF_MORT_ALL",                               # total notif, dead at diagnosis
              paste("NOTIF_MORT",StatList[[1]],sep="_"),      # notif by ag cat, dead at diagnosis
              paste("NOTIF_MORT",StatList[[7]],sep="_"),      # notif by nat cat, dead at diagnosis
              #            paste("NOTIF_MORT_HIV",c("POS","NEG"),sep="_"), # notif by HIV pos/neg, dead at diagnosis
              paste("NOTIF_MORT",StatList[[6]],sep="_"),      # notif by rg cat, dead at diagnosis
              paste("NOTIF_US",StatList[[1]],sep="_"),        # notif by ag cat, US only
              paste("NOTIF_US_MORT",StatList[[1]],sep="_"),   # notif by ag cat, dead at diagnosis US only
              #            paste("NOTIF_MORT_HIV_Neg",StatList[[1]],sep="_"),    # notif by ag cat, dead at diagnosis HIV neg
              #            paste("TOTMORT_W_HIV",StatList[[1]],sep="_"),   # total mort, by ag cat, have HIV
              paste("TOTMORT_W_TB",StatList[[1]],sep="_"),     # total mort, by ag cat, have active TB
              c("N_Ls_US","N_Lf_US","N_Act_US"),
              c("N_Ls_FB","N_Lf_FB","N_Act_FB"),
              c("FOI_LR_US","FOI_HR_US","FOI_LR_FB","FOI_HR_FB" ), #force of infection
              c("TB_INF_LR","TB_INF_HR","TB_INF_US","TB_INF_F1","TB_INF_F2"), # NEW TB INFECTIONS
              c("TBMORT_US","TBMORT_NUS") ,          # tbmortality by nativity
              paste("TOTMORT_US",StatList[[1]],sep="_"),       # mortality by nat cat
              paste("TOTMORT_NUS",StatList[[1]],sep="_"),       # mortality by nat cat

              paste("TOTMORT_US",StatList[[4]],sep="_"),       # mortality by nat & im cat
              paste("TOTMORT_NUS",StatList[[4]],sep="_"),       # mortality by nat & im cat

              paste("TOTMORT_US",StatList[[5]],sep="_"),       # mortality by nat & nm cat
              paste("TOTMORT_NUS",StatList[[5]],sep="_"),       # mortality by nat & nm cat

              paste("TOTMORT_US",StatList[[6]],sep="_"),
              paste("TOTMORT_NUS",StatList[[6]],sep="_"),

              paste("N_US",StatList[[4]],sep="_"),               # pop by nat and im cat
              paste("N_NUS",StatList[[4]],sep="_"),               # pop by nat and im cat

              paste("N_US",StatList[[6]],sep="_"),               # pop by nat and hr cat
              paste("N_NUS",StatList[[6]],sep="_"),               # pop by nat and hr cat

              paste("N_US",StatList[[5]],sep="_"),               # pop by nat and nm cat
              paste("N_NUS",StatList[[5]],sep="_"),              # pop by nat and nm cat
              paste("TOTMORT"),

              paste("N_NM1",StatList[[4]],sep="_"),
              paste("N_NM2",StatList[[4]],sep="_"),
              paste("N_NM3",StatList[[4]],sep="_"),
              paste("N_NM4",StatList[[4]],sep="_"),

              paste("%_0-4","NM1",StatList[[4]],sep="_"),
              paste("%_0-4","NM2",StatList[[4]],sep="_"),
              paste("0-4","NM3",StatList[[4]],sep="_"),
              paste("0-4","NM4",StatList[[4]],sep="_"),

              paste("5-14","NM1",StatList[[4]],sep="_"),
              paste("5-14","NM2",StatList[[4]],sep="_"),
              paste("5-14","NM3",StatList[[4]],sep="_"),
              paste("5-14","NM4",StatList[[4]],sep="_"),

              paste("15-24","NM1",StatList[[4]],sep="_"),
              paste("15-24","NM2",StatList[[4]],sep="_"),
              paste("15-24","NM3",StatList[[4]],sep="_"),
              paste("15-24","NM4",StatList[[4]],sep="_"),

              paste("25-34","NM1",StatList[[4]],sep="_"),
              paste("25-34","NM2",StatList[[4]],sep="_"),
              paste("25-34","NM3",StatList[[4]],sep="_"),
              paste("25-34","NM4",StatList[[4]],sep="_"),

              paste("35-44","NM1",StatList[[4]],sep="_"),
              paste("35-44","NM2",StatList[[4]],sep="_"),
              paste("35-44","NM3",StatList[[4]],sep="_"),
              paste("35-44","NM4",StatList[[4]],sep="_"),

              paste("45-54","NM1",StatList[[4]],sep="_"),
              paste("45-54","NM2",StatList[[4]],sep="_"),
              paste("45-54","NM3",StatList[[4]],sep="_"),
              paste("45-54","NM4",StatList[[4]],sep="_"),

              paste("55-64","NM1",StatList[[4]],sep="_"),
              paste("55-64","NM2",StatList[[4]],sep="_"),
              paste("55-64","NM3",StatList[[4]],sep="_"),
              paste("55-64","NM4",StatList[[4]],sep="_"),

              paste("65-74","NM1",StatList[[4]],sep="_"),
              paste("65-74","NM2",StatList[[4]],sep="_"),
              paste("65-74","NM3",StatList[[4]],sep="_"),
              paste("65-74","NM4",StatList[[4]],sep="_"),

              paste("75-84","NM1",StatList[[4]],sep="_"),
              paste("75-84","NM2",StatList[[4]],sep="_"),
              paste("75-84","NM3",StatList[[4]],sep="_"),
              paste("75-84","NM4",StatList[[4]],sep="_"),

              paste("85-94","NM1",StatList[[4]],sep="_"),
              paste("85-94","NM2",StatList[[4]],sep="_"),
              paste("85-94","NM3",StatList[[4]],sep="_"),
              paste("85-94","NM4",StatList[[4]],sep="_"),

              paste("95p","NM1",StatList[[4]],sep="_"),
              paste("95p","NM2",StatList[[4]],sep="_"),
              paste("95p","NM3",StatList[[4]],sep="_"),
              paste("95p","NM4",StatList[[4]],sep="_"),

              paste("mort_rate",StatList[[1]],sep="_" ),

              paste("N_ag_1",StatList[[5]],sep="_" ),
              paste("N_ag_2",StatList[[5]],sep="_" ),
              paste("N_ag_3",StatList[[5]],sep="_" ),
              paste("N_ag_4",StatList[[5]],sep="_" ),
              paste("N_ag_5",StatList[[5]],sep="_" ),
              paste("N_ag_6",StatList[[5]],sep="_" ),
              paste("N_ag_7",StatList[[5]],sep="_" ),
              paste("N_ag_8",StatList[[5]],sep="_" ),
              paste("N_ag_9",StatList[[5]],sep="_" ),
              paste("N_ag_10",StatList[[5]],sep="_" ),
              paste("N_ag_11",StatList[[5]],sep="_" ),
              ### new infections
              paste("N_newinf_USB",StatList[[1]],sep="_" ),
              paste("N_newinf_NUSB",StatList[[1]],sep="_" ),

              paste("0-24","US", "NM1",StatList[[4]],sep="_"),
              paste("0-24","US", "NM2",StatList[[4]],sep="_"),
              paste("0-24","US", "NM3",StatList[[4]],sep="_"),
              paste("0-24","US", "NM4",StatList[[4]],sep="_"),


              paste("25-64","US","NM1",StatList[[4]],sep="_"),
              paste("25-64","US","NM2",StatList[[4]],sep="_"),
              paste("25-64","US","NM3",StatList[[4]],sep="_"),
              paste("25-64","US","NM4",StatList[[4]],sep="_"),


              paste("65+","US","NM1",StatList[[4]],sep="_"),
              paste("65+","US","NM2",StatList[[4]],sep="_"),
              paste("65+","US","NM3",StatList[[4]],sep="_"),
              paste("65+","US","NM4",StatList[[4]],sep="_"),


              paste("0-24","NUS", "NM1",StatList[[4]],sep="_"),
              paste("0-24","NUS", "NM2",StatList[[4]],sep="_"),
              paste("0-24","NUS", "NM3",StatList[[4]],sep="_"),
              paste("0-24","NUS", "NM4",StatList[[4]],sep="_"),


              paste("25-64","NUS","NM1",StatList[[4]],sep="_"),
              paste("25-64","NUS","NM2",StatList[[4]],sep="_"),
              paste("25-64","NUS","NM3",StatList[[4]],sep="_"),
              paste("25-64","NUS","NM4",StatList[[4]],sep="_"),


              paste("65+","NUS","NM1",StatList[[4]],sep="_"),
              paste("65+","NUS","NM2",StatList[[4]],sep="_"),
              paste("65+","NUS","NM3",StatList[[4]],sep="_"),
              paste("65+","NUS","NM4",StatList[[4]],sep="_"),

              paste("N_LtTxNaive"),

              #counts of various services
              paste("N_LtbiTests_USB",StatList[[1]],sep="_"),
              paste("N_LtbiTests_NUSB",StatList[[1]],sep="_"),

              paste("N_LtbiTxInits_USB",StatList[[1]],sep="_"),
              paste("N_LtbiTxInits_NUSB",StatList[[1]],sep="_"),

              paste("N_LtbiTxComps_USB",StatList[[1]],sep="_"),
              paste("N_LtbiTxComps_NUSB",StatList[[1]],sep="_"),

              paste("N_TBTxInits_USB",StatList[[1]],sep="_"),
              paste("N_TBTxInits_NUSB",StatList[[1]],sep="_"),
              paste("N_TBTxComps_USB",StatList[[1]],sep="_"),
              paste("N_TBTxComps_NUSB",StatList[[1]],sep="_"),

              paste("N_LtbiTests_USB_TP",StatList[[1]],sep="_"),
              paste("N_LtbiTests_NUSB_TP",StatList[[1]],sep="_")
  )
}
