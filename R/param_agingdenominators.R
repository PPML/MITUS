#'THE CODE BELOW SOURCES THE MODEL INPUTS AND THEN FORMATS THEM FOR USE
#'IN THE OUTPUTSINTZ FUNCTION THAT CALLS CSIM FROM THE TB_MODEL.CPP
#'FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME
#'VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.

#'@name age_denom
#' @param loc two digit mailing abbreviation of state
#' @return age_den matrix of denominators 10x1213
#' @export

age_denom<-function(loc){
  if (loc=="US"){
    #read in the population data that is 1 year age bands by each year 1950-2017
    popdist<- as.matrix(readRDS(system.file("US/US_PopCountsByAge.rds", package="MITUS")))
    popdist<-popdist[,-1]
    rownames(popdist)<-as.matrix(readRDS(system.file("US/US_PopCountsByAge.rds", package="MITUS"))[,1])

    #calculate the percentage of the total age band in each single year age
    ltd<-matrix(NA,10,69)
    ltd[1,]<-popdist[5,]/colSums(popdist[1:5,])
    ltd[2,]<-popdist[15,]/colSums(popdist[6:15,])
    ltd[3,]<-popdist[25,]/colSums(popdist[16:25,])
    ltd[4,]<-popdist[35,]/colSums(popdist[26:35,])
    ltd[5,]<-popdist[45,]/colSums(popdist[36:45,])
    ltd[6,]<-popdist[55,]/colSums(popdist[46:55,])
    ltd[7,]<-popdist[65,]/colSums(popdist[56:65,])
    ltd[8,]<-popdist[75,]/colSums(popdist[66:75,])
    ltd[9,]<-popdist[85,]/colSums(popdist[76:85,])
    ltd[10,]<-popdist[95,]/colSums(popdist[86:95,])

    #invert this for the aging rate
    ltd<-1/ltd

    td<-matrix(NA,10,1213)
    for (i in 1:10){
      td[i,1:817]<-SmoCurve(as.numeric(ltd[i,]))
      td[i,818:1213]<-td[i,817]
    }

    age_den<-t(td)*12
  } else {
    #find the state ID number
    #this indexing is based off of the fips data
    data("stateID",package="MITUS")
    StateID<-as.data.frame(stateID)
    index<-which(StateID$USPS==loc)
    st<-as.numeric(as.matrix(StateID[index,2]))
    #read in the population distribution data
    #this data is only decades + 2017
    pop<-readRDS(system.file("ST/ST_PopCountsByAge.rds", package = "MITUS"))
    pop<-pop[pop$fips==st,]
    popdist<-pop[,-(1:3)]
    rownames(popdist)<-pop[,2]
    popdist<-as.matrix(popdist)

    #calculate the percentage of the total age band in each single year age
    ltd<-matrix(NA,10,8)
    ltd[1,]<-popdist[5,]/colSums(popdist[1:5,])
    ltd[2,]<-popdist[15,]/colSums(popdist[6:15,])
    ltd[3,]<-popdist[25,]/colSums(popdist[16:25,])
    ltd[4,]<-popdist[35,]/colSums(popdist[26:35,])
    ltd[5,]<-popdist[45,]/colSums(popdist[36:45,])
    ltd[6,]<-popdist[55,]/colSums(popdist[46:55,])
    ltd[7,]<-popdist[65,]/colSums(popdist[56:65,])
    ltd[8,]<-popdist[75,]/colSums(popdist[66:75,])
    ltd[9,]<-popdist[85,]/colSums(popdist[76:85,])
    ltd[10,]<-popdist[95,]/colSums(popdist[86:95,])


    for (i in 1:length(ltd)){
    if (ltd[i]==0){
      ltd[i]<-0.01980337
    }}
    #invert this for the aging rate

    ltd<-1/ltd
    td<-matrix(NA,10,1213)
    for (i in 1:10){
      td[i,1:841]<-SmoCurve_decade_month(as.numeric(ltd[i,]))
      td[i,842:1213]<-td[i,841]
    }
  age_den<-t(td)*12
  } #end of else statement
  return(age_den)
}
