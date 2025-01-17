#'Function that will return default values of the program change values.
#'
#'@name def_costinputs
#'@return vector of the default values
#'@export

def_costinputs<-function(){
  #create an empty vector to hold the values
  DefcostinputsVec<-rep(NA,12)
  names(DefcostinputsVec)<-c('LTBIIdCost','TSTCost','IGRACost','NoTBCost',
                  '3HPCost','4RCost','3HRCost','TBIdCost', 'TBtest',
                  'TBtx', 'Discount', 'EndYear')

  DefcostinputsVec['LTBIIdCost']<-0
  DefcostinputsVec['TSTCost']<-9.38
  DefcostinputsVec['IGRACost']<-61.98
  DefcostinputsVec['NoTBCost']<-33.20

  DefcostinputsVec['3HPCost']<-411.87
  DefcostinputsVec['4RCost']<-354.87
  DefcostinputsVec['3HRCost']<-353.61

  DefcostinputsVec['TBIdCost']<-0
  #TST,XRAY,Sputum culture, & susceptibility
  #9.38+33.20+20.05+7
  DefcostinputsVec['TBtest']<-0
  DefcostinputsVec['TBtx']<-20211

  DefcostinputsVec['Discount']<-3
  DefcostinputsVec["EndYear"]<-2050

  return(DefcostinputsVec)
}
