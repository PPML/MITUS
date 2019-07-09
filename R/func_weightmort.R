#'THE CODE BELOW SOURCES THE MODEL INPUTS AND THEN FORMATS THEM FOR USE
#'IN THE OUTPUTSINTZ FUNCTION THAT CALLS CSIM FROM THE TB_MODEL.CPP
#'FUNCTION FILE. ALL VARIABLE NAMES THAT END IN t ARE INDEXED BY TIME
#'VARIABLE NAMES BEGINNING WITH m ARE MATRICES & V ARE VECTORS.

#'@name weight_mort
#' @param loc two digit mailing abbreviation of state
#' @return  matrix of denominators 10x1201
#' @export

weight_mort<-function(loc){
  if (loc=="US"){
    death_age <-readRDS(system.file("US/US_MortalityCountsByAge.rds", package="MITUS"))[,2:69]
    rownames(death_age)<-readRDS(system.file("US/US_MortalityCountsByAge.rds", package="MITUS"))[,1]
    death_age<-as.matrix(death_age)

    popdist<- as.matrix(readRDS(system.file("US/US_PopCountsByAge.rds", package="MITUS")))
    popdist<-popdist[,-1]
    popdist<-popdist[,-69]
    rownames(popdist)<-as.matrix(readRDS(system.file("US/US_PopCountsByAge.rds", package="MITUS"))[,1])
    new_mort<-popdist
    new_mort<-death_age/popdist
    ##weighted mort for the age groups
    weight_mort<-matrix(NA,11,68)
    colnames(weight_mort)<-colnames(new_mort)
    weight_mort[1,]<-colSums(new_mort[1:5,]*popdist[1:5,])/colSums(popdist[1:5,])
    weight_mort[2,]<-colSums(new_mort[6:15,]*popdist[6:15,])/colSums(popdist[6:15,])
    weight_mort[3,]<-colSums(new_mort[16:25,]*popdist[16:25,])/colSums(popdist[16:25,])
    weight_mort[4,]<-colSums(new_mort[26:35,]*popdist[26:35,])/colSums(popdist[26:35,])
    weight_mort[5,]<-colSums(new_mort[36:45,]*popdist[36:45,])/colSums(popdist[36:45,])
    weight_mort[6,]<-colSums(new_mort[46:55,]*popdist[46:55,])/colSums(popdist[46:55,])
    weight_mort[7,]<-colSums(new_mort[56:65,]*popdist[56:65,])/colSums(popdist[56:65,])
    weight_mort[8,]<-colSums(new_mort[66:75,]*popdist[66:75,])/colSums(popdist[66:75,])
    weight_mort[9,]<-colSums(new_mort[76:85,]*popdist[76:85,])/colSums(popdist[76:85,])
    weight_mort[10,]<-colSums(new_mort[86:95,]*popdist[86:95,])/colSums(popdist[86:95,])
    weight_mort[11,]<-colSums(new_mort[96:111,]*popdist[96:111,])/colSums(popdist[96:111,])

    weight_mort<-t(weight_mort)
  } else {
    #find the state ID number
    #this indexing is based off of the fips data
    data("stateID",package="MITUS")
    StateID<-as.data.frame(stateID)
    index<-which(StateID$USPS==loc)
    #read in the lifetable data
    ST_mortrate<-readRDS(system.file("ST/ST_SingleYearLifeTable.rds",package="MITUS"))[[index]][,-1]

    #calculate the crude mortality rates for the US national
    death_age <-readRDS(system.file("US/US_MortalityCountsByAge.rds", package="MITUS"))[,2:69]
    rownames(death_age)<-readRDS(system.file("US/US_MortalityCountsByAge.rds", package="MITUS"))[,1]
    death_age<-as.matrix(death_age)

    popdist<- as.matrix(readRDS(system.file("US/US_PopCountsByAge.rds", package="MITUS")))
    popdist<-popdist[,-1]
    popdist<-popdist[,-69]
    rownames(popdist)<-as.matrix(readRDS(system.file("US/US_PopCountsByAge.rds", package="MITUS"))[,1])
    nat_mort<-popdist
    nat_mort<-(death_age/popdist)[,10:67]

    #ratio between state/national
    ratio<-ST_mortrate/nat_mort
    #calculate national weighted mort
    weight_mort<-matrix(NA,11,58)
    colnames(weight_mort)<-colnames(nat_mort)
    weight_mort[1,]<-colSums(nat_mort[1:5,]*popdist[1:5,10:67]*ratio[1:5,])/colSums(popdist[1:5,10:67])
    weight_mort[2,]<-colSums(nat_mort[6:15,]*popdist[6:15,10:67]*ratio[5:15,])/colSums(popdist[6:15,10:67])
    weight_mort[3,]<-colSums(nat_mort[16:25,]*popdist[16:25,10:67]*ratio[16:25,])/colSums(popdist[16:25,10:67])
    weight_mort[4,]<-colSums(nat_mort[26:35,]*popdist[26:35,10:67]*ratio[26:35,])/colSums(popdist[26:35,10:67])
    weight_mort[5,]<-colSums(nat_mort[36:45,]*popdist[36:45,10:67]*ratio[36:45,])/colSums(popdist[36:45,10:67])
    weight_mort[6,]<-colSums(nat_mort[46:55,]*popdist[46:55,10:67]*ratio[46:55,])/colSums(popdist[46:55,10:67])
    weight_mort[7,]<-colSums(nat_mort[56:65,]*popdist[56:65,10:67]*ratio[56:65,])/colSums(popdist[56:65,10:67])
    weight_mort[8,]<-colSums(nat_mort[66:75,]*popdist[66:75,10:67]*ratio[66:75,])/colSums(popdist[66:75,10:67])
    weight_mort[9,]<-colSums(nat_mort[76:85,]*popdist[76:85,10:67]*ratio[76:85,])/colSums(popdist[76:85,10:67])
    weight_mort[10,]<-colSums(nat_mort[86:95,]*popdist[86:95,10:67]*ratio[86:95,])/colSums(popdist[86:95,10:67])
    weight_mort[11,]<-colSums(nat_mort[96:111,]*popdist[96:111,10:67]*ratio[96:111,])/colSums(popdist[96:111,10:67])

weight_mort<-t(weight_mort)
}
  return (weight_mort)
}
