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
    fips<-as.numeric(as.matrix(StateID[index,2]))
    #read in the population distribution data
    #this data is only decades + 2017
    pop<-readRDS(system.file("ST/ST_PopCountsByAge.rds", package = "MITUS"))
    pop<-pop[pop$fips==fips,]
    popdist<-pop[,-(1:3)]
    rownames(popdist)<-pop[,2]
    popdist<-as.matrix(popdist)
    ###load in the mortality rate from the life tables
    #find the state ID number (this is ordinal, not fips)
    st<-which(StateID$USPS==loc)
    #load in the state life table
    ST_lifetable<-readRDS(system.file("ST/ST_lifetables.rds", package="MITUS"))
    lt<-ST_lifetable[[st]]
    lt_d<-reshape2::dcast(lt,Age~Year,value.var = "mx")
    lt_d$Age<-as.character(lt_d$Age)
    lt_d$Age[lt_d$Age == '110+'] <- '110'
    mort_rate<-dplyr::arrange(lt_d,as.integer(Age))
    mort_rate<-matrix(as.numeric(unlist(mort_rate)),111,7)
   rownames(mort_rate)<-mort_rate[,1]
     mort_rate<-mort_rate[,-1]
     colnames(mort_rate)<-c("1960","1970","1980","1990","2000","2010")
     mort_rate<-mort_rate[1:101,]
    dec_weight_mort<-matrix(NA,11,6)
    popdist<-popdist[,2:7]
    dec_weight_mort[1,]<-colSums(mort_rate[1:5,]*popdist[1:5,])/colSums(popdist[1:5,])
    dec_weight_mort[2,]<-colSums(mort_rate[6:15,]*popdist[6:15,])/colSums(popdist[6:15,])
    dec_weight_mort[3,]<-colSums(mort_rate[16:25,]*popdist[16:25,])/colSums(popdist[16:25,])
    dec_weight_mort[4,]<-colSums(mort_rate[26:35,]*popdist[26:35,])/colSums(popdist[26:35,])
    dec_weight_mort[5,]<-colSums(mort_rate[36:45,]*popdist[36:45,])/colSums(popdist[36:45,])
    dec_weight_mort[6,]<-colSums(mort_rate[46:55,]*popdist[46:55,])/colSums(popdist[46:55,])
    dec_weight_mort[7,]<-colSums(mort_rate[56:65,]*popdist[56:65,])/colSums(popdist[56:65,])
    dec_weight_mort[8,]<-colSums(mort_rate[66:75,]*popdist[66:75,])/colSums(popdist[66:75,])
    dec_weight_mort[9,]<-colSums(mort_rate[76:85,]*popdist[76:85,])/colSums(popdist[76:85,])
    dec_weight_mort[10,]<-colSums(mort_rate[86:95,]*popdist[86:95,])/colSums(popdist[86:95,])
    dec_weight_mort[11,]<-colSums(mort_rate[96:101,]*popdist[96:101,])/colSums(popdist[96:101,])

    for (i in 1:length(dec_weight_mort)){
    if(is.nan(dec_weight_mort[i])==TRUE) dec_weight_mort[i]<-0
}
spl_mort<-matrix(NA,11,51)
    for (i in 1:nrow(dec_weight_mort)){
      spl_mort[i,]<-SmoCurve_decade_year(dec_weight_mort[i,])
    }
# for (i in 1:11){
#   plot(spl_mort[i,],type="l")
# }
weight_mort<-t(spl_mort)
}
  return (weight_mort)
}
