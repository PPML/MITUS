#format state level input data
state_input <- function() {
  load("/Users/nis100/Desktop/untitled folder/ModelInputsState_6-23-18.rData")
  stateID<-as.matrix(read.csv(file="inst/extdata/state_ID.csv", header = TRUE))
  for (i in 1:51){
    loc<-stateID[i,3]
    Inputs<-list()
    Inputs$BgMort<-InputsState$BgMort
    Inputs$InitPop<-InputsState$InitPopState[[i]]
    Inputs$Births<-as.matrix(InputsState$BirthsState[[i]])
    Inputs$Births<-as.vector(Inputs$Births[,2])
    Inputs$TxInputs<-InputsState$TxInputs
    Inputs$NetMigrState<-InputsState$NetMigrState[[i]]

    ImmigInputs<-list()
    ImmigInputs[["TotByYear"]]<-as.numeric(InputsState$ImmigInputsState$TotByYearState[[i]][,-1])
    ImmigInputs[["AgeDist"]]<-as.numeric(InputsState$ImmigInputsState[["AgeDistState"]][[i]])
    ImmigInputs[["Func_LtbiAgeDist"]]<-InputsState$ImmigInputsState[["Func_LtbiAgeDist"]]
    ImmigInputs[["LtbiPars"]]<-InputsState$ImmigInputsState[["LtbiPars"]]
    ImmigInputs[["PrevTrend25_34"]]<-InputsState$ImmigInputsState[["PrevTrend25_34"]]
    ImmigInputs[["RR_Active_TB_Age"]]<-InputsState$ImmigInputsState[["RR_Active_TB_Age"]]
    ImmigInputs[["Frac_TxE_by_Age"]]<-InputsState$ImmigInputsState[["Frac_TxE_by_Age"]]
    ImmigInputs[["DR_TB_by_year_N"]]<-InputsState$ImmigInputsState[["DR_TB_by_year_N"]]
    ImmigInputs[["DR_TB_by_year_E"]]<-InputsState$ImmigInputsState[["DR_TB_by_year_E"]]
    Inputs$ImmigInputs<-ImmigInputs

    assign(paste0(loc,"_Inputs"),Inputs)
    saveRDS(get(paste0(loc,"_Inputs")), file=paste0("~/MITUS/inst/",loc,"/",loc,"_ModelInputs_01-24-19.rds"))
  }

}

