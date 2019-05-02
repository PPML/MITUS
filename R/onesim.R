onesim<-function(loc){
  nsres<-list()
  nbres<-list()

  if(loc %in% c("CA","FL","GA","IL","NY","WA")){
    i<-1
  } else if (loc %in% c("TX","VA")){
    i<-10
  } else if (loc == "NJ"){
    i <-9
  } else if (loc == "PA"){
    i <-7
  }
  # bg_resTab_2019-05-02
    o<-readRDS(system.file(paste0(loc,"/bg_resTab_2019-04-23.rds"), package="MITUS"))
    nbres$ResTab<-o$ResTab[((10*(0:8))+i),,,]
    nbres$ResTabus<-o$ResTabus[((10*(0:8))+i),,,]
    nbres$ResTabfb<-o$ResTabfb[((10*(0:8))+i),,,]
    saveRDS(nbres, file=paste0("~/MITUS/inst/",loc,"/bg_ResTab_",Sys.Date(),".rds"))

    o<-readRDS(system.file(paste0(loc,"/sm_resTab_2019-04-23.rds"), package="MITUS"))
    nsres$ResTab<-o$ResTab[((10*(0:8))+i),,,]
    nsres$ResTabus<-o$ResTabus[((10*(0:8))+i),,,]
    nsres$ResTabfb<-o$ResTabfb[((10*(0:8))+i),,,]
    saveRDS(nsres, file=paste0("~/MITUS/inst/",loc,"/sm_ResTab_",Sys.Date(),".rds"))
}
