#'Function that will return default values of the program change values.
#'
#'@name def_prgchng
#'@param ParVec parameter vector
#'@return vector of the default values
#'@export

def_prgchng<-function(ParVec){
  #create an empty vector to hold the values
  DefPrgChngVec<-rep(NA,8)
  names(DefPrgChngVec)<-c("start_yr", #year in which the program change starts (discontinuous step up to the values below at this year)
                          "scrn_cov", #Screening Coverage Rate as a Multiple of the Current Rate
                          "IGRA_frc", #Fraction of Individuals Receiving IGRA
                          "ltbi_init_frc", #Fraction of Individuals Testing Positive who Accept Treatment
                          "ltbi_comp_frc", #Fraction of Individuals Initiating Treatment Who Complete Treatment
                          "ltbi_eff_frc",
                          "tb_tim2tx_frc", #Duration of Infectiousness
                          "tb_txdef_frc" #Fraction Discontinuing/Defaulting from Treatment
  )
  #default start year will always be 2020
  DefPrgChngVec[1]<-2020
  #default screening coverage multiplier will always default to 1
  DefPrgChngVec[2]<-1
  #the default IGRA fraction is a constant that is not calibrated or calculated
  #once this is set it should not be changed
  DefPrgChngVec[3]<-.50
  #Treatment Initiation Fraction is a constant that is not calibrated or calculated
  #once this is set it should not be changed
  DefPrgChngVec[4]<-.773
  #LTBI treatment completion fraction
  DefPrgChngVec[5]<- 1-ParVec["pDefLt"]
  #LTBI Efficacy
  DefPrgChngVec[6]<-ParVec["EffLt"]
  #Time to Treatment //Duration of Infectiousness Percent of Current Value
  DefPrgChngVec[7]<-100
  #Fraction Discontinuing/Defaulting from Treatment
  TxInputs         <- Inputs[["TxInputs"]]
  rDef0         <- rep(NA,151)
  rDef0[1:30]   <- ParVec["TxDefEarly"]
  rDef0[44:63]  <- ORAdd(TxInputs[[1]][,2],ParVec["TunTxDef"])
  rDef0[64:151] <- rDef0[63]
  rDef1         <- predict(smooth.spline(x=c(1950:1979,1993:2100),y=rDef0[-(31:43)],spar=0.4),x=1950:2100)$y
  rDeft         <- SmoCurve(rDef1)/12;
  rDef<-rDeft[(2020-1950)+1]
  DefPrgChngVec[8]<-rDef

  return(DefPrgChngVec)
}
