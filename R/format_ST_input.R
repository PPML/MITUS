#format state level input data
state_input <- function() {
  load("/Users/nis100/Desktop/desktop3/untitled folder/ModelInputsState_6-23-18.rData")
  stateID<-as.matrix(read.csv(file="inst/extdata/state_ID.csv", header = TRUE))
  for (i in 1:51){
    loc<-stateID[i,3]
    Inputs<-list()
    Inputs$BgMort<-InputsState$BgMort
    ##Create the Initial Pop from the Decennial Census data
    pop<-readRDS(system.file("ST/ST_decennialpop.rds", package="MITUS"))[[i]][,1:3]
    pop2<-reshape2::dcast(pop,age_group~usb,value.var = "1950")
    rownames(pop2)<-pop2[,1]
    pop2<-pop2[,-1]
    colnames(pop2)<-c("US","NUS")
    Inputs$InitPop<-as.matrix(pop2)/1e6
    Inputs$Births<-as.matrix(InputsState$BirthsState[[i]])
    Inputs$Births<-as.vector(Inputs$Births[,2])
    Inputs$TxInputs<-InputsState$TxInputs
    #Create the Net Migration from the ST_netmigration.rds file
    netmig<-readRDS(system.file("ST/ST_netmigration.rds", package="MITUS"))[[i]]
    netmig2<-reshape2::dcast(netmig,age_group~usb,value.var = "net_mig_m")
    rownames(netmig2)<-netmig2[,1]
    netmig2<-netmig2[,-1]
    colnames(netmig2)<-c("usb","nusb")
    Inputs$NetMigrState<-netmig2

    ImmigInputs<-list()
    #load in the data
    immig<-readRDS(system.file("ST/ST_immigration.rds", package="MITUS"))[[i]]/1e6
    ImmigInputs[["TotByYear"]]<-as.numeric(colSums(immig))
    ImmigInputs[["AgeDist"]]<-immig
    ImmigInputs[["Func_LtbiAgeDist"]]<-InputsState$ImmigInputsState[["Func_LtbiAgeDist"]]
    ImmigInputs[["LtbiPars"]]<-InputsState$ImmigInputsState[["LtbiPars"]]
    ImmigInputs[["PrevTrend25_34"]]<-InputsState$ImmigInputsState[["PrevTrend25_34"]]
    #Immigration Burden Trend (needs to be multiplied by something)
    burden<-readRDS(system.file("ST/ST_TBburdenimmig.rds", package="MITUS"))
    ImmigInputs[["TBBurdenImmig"]]<-burden
    ImmigInputs[["RR_Active_TB_Age"]]<-InputsState$ImmigInputsState[["RR_Active_TB_Age"]]
    ImmigInputs[["Frac_TxE_by_Age"]]<-InputsState$ImmigInputsState[["Frac_TxE_by_Age"]]
    ImmigInputs[["DR_TB_by_year_N"]]<-InputsState$ImmigInputsState[["DR_TB_by_year_N"]]
    ImmigInputs[["DR_TB_by_year_E"]]<-InputsState$ImmigInputsState[["DR_TB_by_year_E"]]
    Inputs$ImmigInputs<-ImmigInputs

    # assign(paste0(loc,"_Inputs"),Inputs)
    saveRDS(Inputs, file=paste0("~/MITUS/inst/",loc,"/",loc,"_ModelInputs_07-03-19.rds"))
  }

}

