format_as_restab <- function(loc) {

  age_id = (2018:2049)-1949
  tt<-0
  #create 3 lists to hold output
  ResTabfb <- ResTabus <- ResTab <- array(NA,dim=c(90,length(age_id),7,11))

  # For now just set intv = 1 for a base-case scenario, but later we need to change this
  # to reflect that we're working on custom scenarios
  for (intv in 1:9){
    #load the results for all the runs (need to make this dataset)
  o<-readRDS(system.file(paste0(loc,"/results_", intv, "_2019-05-08.rds"), package="MITUS"))
  # for (i in 1:564){
  #   o[,,i]<-o[,,i]*1e3
  # }
  #dimensions of restab are:
  #number of datasets, number of ages, number of interventions
  nr<-10
  #gather the outputs from that model run
  for (ag in 1:11){
    ################################################################################
    ################################################################################
    #total population
    #dimensions are
    #scenarios; length age id;
    # print(1:nr+((intv-1)*10))
    ################################################################################
    #get the years
    ResTab[1:nr+((intv-1)*10),,1,ag]<-o[,age_id,1]
    #number of ltbi prevalence in thousands
    ResTab[1:nr+((intv-1)*10),,2,ag]<-apply(o[,age_id,c(54,65)+ag],c(1,2),sum)*1e3
    #percentage of ltbi prevalence
    ResTab[1:nr+((intv-1)*10),,3,ag]<-apply(o[,age_id,c(54,65)+ag],c(1,2),sum)/o[, age_id,32+ag]*1e2
    #TB notifications (incidence) (alive+dead at diagnosis) in thousands
    ResTab[1:nr+((intv-1)*10),,4,ag]<-apply(o[,age_id,c(135,188)+ag],c(1,2),sum)*1e3
    #TB notifications (incidence) (alive+dead at diagnosis) per million
    ResTab[1:nr+((intv-1)*10),,5,ag]<-apply(o[,age_id,c(135,188)+ag],c(1,2),sum)/o[, age_id,32+ag]*1e6
    #tb attributable deaths in thousands
    ResTab[1:nr+((intv-1)*10),,6,ag]<-apply(o[,age_id,c(87,98)+ag],c(1,2),sum)*1e3
    #  tb attributable deaths per million
    ResTab[1:nr+((intv-1)*10),,7,ag]<-apply(o[,age_id,c(87,98)+ag],c(1,2),sum)/o[, age_id,32+ag]*1e6

    ################################################################################
    #US Born population
    ################################################################################
    #get da years
    ResTabus[1:nr+((intv-1)*10),,1,ag]<-o[,age_id,1]
    #number of ltbi prevalence in thousands
    ResTabus[1:nr+((intv-1)*10),,2,ag]<-o[, age_id,54+ag]*1e3
    #percentage of ltbi prevalence
    ResTabus[1:nr+((intv-1)*10),,3,ag]<-o[, age_id,54+ag]/o[, age_id,43+ag]*1e2
    #TB notifications (alive+dead at diagnosis) in thousands
    ResTabus[1:nr+((intv-1)*10),,4,ag]<-apply(o[, age_id,c(204,215)+ag],c(1,2),sum)*1e3
    #TB notifications (alive+dead at diagnosis) per million
    ResTabus[1:nr+((intv-1)*10),,5,ag]<-apply(o[, age_id,c(204,215)+ag],c(1,2),sum)/o[, age_id,43+ag]*1e6
    #tb attributable deaths in thousands
    ResTabus[1:nr+((intv-1)*10),,6,ag]<-o[, age_id,87+ag]*1e3
    #  tb attributable deaths per million
    ResTabus[1:nr+((intv-1)*10),,7,ag]<-o[, age_id,87+ag]/o[, age_id,43+ag]*1e6

    ################################################################################
    #non-US Born population
    ################################################################################
    ResTabfb[1:nr+((intv-1)*10),,1,ag]<-o[,age_id,1]
    #number of ltbi prevalence in thousands
    ResTabfb[1:nr+((intv-1)*10),,2,ag]<-o[, age_id,65+ag]*1e3
    #percentage of ltbi prevalence
    ResTabfb[1:nr+((intv-1)*10),,3,ag]<-o[, age_id,65+ag]/o[, age_id,2+ag]*1e2
    #TB notifications (alive+dead at diagnosis)
    #calculated as total notifications minus US only notifications in thousands
    ResTabfb[1:nr+((intv-1)*10),,4,ag]<-(apply(o[, age_id,c(135,188)+ag],c(1,2),sum)-apply(o[, age_id,c(204,215)+ag],c(1,2),sum))*1e3
    #percentage TB notifications (alive+dead at diagnosis) per million
    ResTabfb[1:nr+((intv-1)*10),,5,ag]<-(apply(o[, age_id,c(135,188)+ag],c(1,2),sum)-apply(o[, age_id,c(204,215)+ag],c(1,2),sum))/o[, age_id,43+ag]*1e6
    #tb attributable deaths
    ResTabfb[1:nr+((intv-1)*10),,6,ag]<-o[, age_id,98+ag]*1e3
    # percentage tb attributable deaths
    ResTabfb[1:nr+((intv-1)*10),,7,ag]<-o[, age_id,98+ag]/o[, age_id,2+ag]*1e6
  }
  }
  results<-list(ResTab = ResTab,
    ResTabus = ResTabus,
    ResTabfb = ResTabfb)
  saveRDS(results, file=paste0("~/MITUS/inst/", loc, "/sm_resTab_", Sys.Date(), ".rds"))


################ BIG RESULTS TAB #############################################

  age_id = (2018:2049)-1949
  tt<-0
  #create 3 lists to hold output
  mResTabfb <- mResTabus <- mResTab <- array(NA,dim=c(9,length(age_id),5,4))

  ResTabfb <- ResTabus <- ResTab <- array(NA,dim=c(90,length(age_id),5,4))
  for (intv in 1:9){
    #load the results for all the runs (need to make this dataset)
    o<-readRDS(system.file(paste0(loc,"/results_", intv, "_2019-05-08.rds"), package="MITUS"))
    # for (i in 1:564){
    #   o[,,i]<-o[,,i]*1e3
    # }
    #dimensions of restab are:
    #number of datasets, number of ages, number of interventions
    nr<-10
    #gather the outputs from that model run
    for (b_ag in 1:4){
      ################################################################################
      #total population
      #dimensions are
      #scenarios; length age id;
      # print(1:nr+((intv-1)*10))
      ag_vec<- switch(b_ag,1:11, c(1,2,3),c(4,5,6,7),c(8,9,10,11))
      ################################################################################
      # for (ag_vec in ag_vector){
        ResTab[1:nr+((intv-1)*10),,1,b_ag]<-o[,age_id,1]
        #percentage of ltbi prevalence
        ResTab[1:nr+((intv-1)*10),,2,b_ag]<-apply(o[,age_id,c(54+ag_vec,65+ag_vec)],c(1,2),sum)/rowSums(o[, age_id,2+ag_vec], dims = 2)*1e2
        #incident tb infections per million
        ResTab[1:nr+((intv-1)*10),,3,b_ag]<-apply(o[,age_id,c(564+ag_vec,575+ag_vec)],c(1,2),sum)/rowSums(o[, age_id,2+ag_vec], dims = 2)*1e6
        # TB notifications (alive+dead at diagnosis) per million
        ResTab[1:nr+((intv-1)*10),,4,b_ag]<-apply(o[,age_id,c(135+ag_vec,188+ag_vec)],c(1,2),sum)/rowSums(o[, age_id,2+ag_vec], dims = 2)*1e6
        # tb attributable deaths per million
        ResTab[1:nr+((intv-1)*10),,5,b_ag]<-apply(o[,age_id,c(87+ag_vec,98+ag_vec)],c(1,2),sum)/rowSums(o[, age_id,2+ag_vec], dims = 2)*1e6

        ################################################################################
        #US Born population
        ################################################################################
        ResTabus[1:nr+((intv-1)*10),,1,b_ag]<-o[,age_id,1]/1e6
        # #percentage of ltbi prevalence
        ResTabus[1:nr+((intv-1)*10),,2,b_ag]<-(rowSums(o[, age_id,54+ag_vec], dims = 2)/rowSums(o[, age_id,32+ag_vec], dims = 2))*1e2
        #incident tb infections per million
        ResTabus[1:nr+((intv-1)*10),,3,b_ag]<-rowSums(o[, age_id,564+ag_vec], dims=2)/rowSums(o[, age_id,32+ag_vec], dims = 2)*1e6
        # # TB notifications (alive+dead at diagnosis) per million
        ResTabus[1:nr+((intv-1)*10),,4,b_ag]<-apply(o[,age_id,c(204+ag_vec,215+ag_vec)],c(1,2),sum)/rowSums(o[, age_id,32+ag_vec], dims = 2)*1e6
        # #  attributable deaths per million
        ResTabus[1:nr+((intv-1)*10),,5,b_ag]<-rowSums(o[, age_id,87+ag_vec], dims = 2)/rowSums(o[, age_id,32+ag_vec], dims = 2)*1e6
        #
        # ################################################################################
        # #non-US Born population
        # ################################################################################
        ResTabfb[1:nr+((intv-1)*10),,1,b_ag]<-o[,age_id,1]/1e6
        # #percentage of ltbi prevalence
        ResTabfb[1:nr+((intv-1)*10),,2,b_ag]<-rowSums(o[, age_id,65+ag_vec], dims=2)/rowSums(o[, age_id,43+ag_vec], dims = 2)*1e2
        #incident tb infections per million
        ResTabfb[1:nr+((intv-1)*10),,3,b_ag]<-rowSums(o[, age_id,575+ag_vec], dims=2)/rowSums(o[, age_id,43+ag_vec], dims = 2)*1e6
        # TB notifications (alive+dead at diagnosis) per million
        ResTabfb[1:nr+((intv-1)*10),,4,b_ag]<-(apply(o[,age_id,c(135+ag_vec,188+ag_vec)],c(1,2),sum)-apply(o[,age_id,c(204+ag_vec,215+ag_vec)],c(1,2),sum))/rowSums(o[, age_id,43+ag_vec], dims = 2)*1e6
        # tb attributable deaths per million
        ResTabfb[1:nr+((intv-1)*10),,5,b_ag]<-rowSums(o[, age_id,98+ag_vec],dims = 2)/rowSums(o[, age_id,43+ag_vec], dims = 2)*1e6
      # }
    }
    # print(paste(aa,"--",bb)); flush.console()
  } #end batch loop
  # ResTabAll[[tt]]   <- ResTab
  # ResTabAllfb[[tt]] <- ResTabfb
  # ResTabAllus[[tt]] <- ResTabus
  # }
  #concatenate all lists

  results<-list(ResTab = ResTab,
                ResTabus = ResTabus,
                ResTabfb = ResTabfb)
  saveRDS(results, file=paste0("~/MITUS/inst/", loc, "/bg_resTab_", Sys.Date(), ".rds"))
}
