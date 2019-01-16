#'@name tabby_results_s
#'@description reshapes output from OutputsInt to shape needed for Tabby2
#'@param loc two digit postal abbreviation for location of interest
#'@return res_tab2 large dataframe of results
#'@export
tabby_results_s<-function(loc){
  age_id = (2018:2049)-1949
  tt<-0
  #create 3 lists to hold output
  mResTabfb <- mResTabus <- mResTab <- array(NA,dim=c(9,length(age_id),7,11))

  ResTabfb <- ResTabus <- ResTab <- array(NA,dim=c(90,length(age_id),7,11))
  for (intv in 1:9){
    #load the results for all the runs (need to make this dataset)
    load(paste0("~/MITUS/",loc,"_results_",intv,".rda"))
      #dimensions of restab are:
      #number of datasets, number of ages, number of interventions
    nr<-10
    #define o as the results
    o<-out
#gather the outputs from that model run
    for (ag in 1:11){
################################################################################
      #total population
      #dimensions are
      #scenarios; length age id;
      # print(1:nr+((intv-1)*10))
################################################################################
      ResTab[1:nr+((intv-1)*10),,1,ag]<-o[,age_id,1]
      #number of ltbi prevalence
      ResTab[1:nr+((intv-1)*10),,2,ag]<-apply(o[,age_id,c(54,65)+ag],c(1,2),sum)*1e3
      #percentage of ltbi prevalence
      ResTab[1:nr+((intv-1)*10),,3,ag]<-apply(o[,age_id,c(54,65)+ag],c(1,2),sum)/o[, age_id,32+ag]*1e2
      #TB notifications (alive+dead at diagnosis)
      ResTab[1:nr+((intv-1)*10),,4,ag]<-apply(o[,age_id,c(135,188)+ag],c(1,2),sum)*1e3
      #percentage TB notifications (alive+dead at diagnosis)
      ResTab[1:nr+((intv-1)*10),,5,ag]<-apply(o[,age_id,c(135,188)+ag],c(1,2),sum)/o[, age_id,32+ag]*1e2
      #tb attributable deaths
      ResTab[1:nr+((intv-1)*10),,6,ag]<-apply(o[,age_id,c(87,98)+ag],c(1,2),sum)*1e3
      # percentage tb attributable deaths
      ResTab[1:nr+((intv-1)*10),,7,ag]<-apply(o[,age_id,c(87,98)+ag],c(1,2),sum)/o[, age_id,32+ag]*1e2

      ################################################################################
      #US Born population
      ################################################################################
      ResTabus[1:nr+((intv-1)*10),,1,ag]<-o[,age_id,1]
      #number of ltbi prevalence
      ResTabus[1:nr+((intv-1)*10),,2,ag]<-o[, age_id,54+ag]*1e3
      #percentage of ltbi prevalence
      ResTabus[1:nr+((intv-1)*10),,3,ag]<-o[, age_id,54+ag]/o[, age_id,43+ag]*1e2
      #TB notifications (alive+dead at diagnosis)
      ResTabus[1:nr+((intv-1)*10),,4,ag]<-apply(o[, age_id,c(204,215)+ag],c(1,2),sum)*1e3
      #percentage TB notifications (alive+dead at diagnosis)
      ResTabus[1:nr+((intv-1)*10),,5,ag]<-apply(o[, age_id,c(204,215)+ag],c(1,2),sum)/o[, age_id,43+ag]*1e2
      #tb attributable deaths
      ResTabus[1:nr+((intv-1)*10),,6,ag]<-o[, age_id,87+ag]*1e3
      # percentage tb attributable deaths
      ResTabus[1:nr+((intv-1)*10),,7,ag]<-o[, age_id,87+ag]/o[, age_id,43+ag]*1e2

      ################################################################################
      #non-US Born population
      ################################################################################
      ResTabfb[1:nr+((intv-1)*10),,1,ag]<-o[,age_id,1]
      #number of ltbi prevalence
      ResTabfb[1:nr+((intv-1)*10),,2,ag]<-o[, age_id,65+ag]*1e3
      #percentage of ltbi prevalence
      ResTabfb[1:nr+((intv-1)*10),,3,ag]<-o[, age_id,65+ag]/o[, age_id,2+ag]*1e2
      #TB notifications (alive+dead at diagnosis)
      ResTabfb[1:nr+((intv-1)*10),,4,ag]<-(apply(o[, age_id,c(135,188)+ag],c(1,2),sum)-apply(o[, age_id,c(204,215)+ag],c(1,2),sum))*1e3
      #percentage TB notifications (alive+dead at diagnosis)
      ResTabfb[1:nr+((intv-1)*10),,5,ag]<-(apply(o[, age_id,c(135,188)+ag],c(1,2),sum)-apply(o[, age_id,c(204,215)+ag],c(1,2),sum))/o[, age_id,43+ag]*1e2
      #tb attributable deaths
      ResTabfb[1:nr+((intv-1)*10),,6,ag]<-o[, age_id,98+ag]*1e3
      # percentage tb attributable deaths
      ResTabfb[1:nr+((intv-1)*10),,7,ag]<-o[, age_id,98+ag]/o[, age_id,2+ag]*1e2
    }
    # print(paste(aa,"--",bb)); flush.console()
    } #end batch loop
  # ResTabAll[[tt]]   <- ResTab
  # ResTabAllfb[[tt]] <- ResTabfb
  # ResTabAllus[[tt]] <- ResTabus
  # }
  #concatenate all lists

  ################################################################################
  #take the mean for each of these
  for (i in 1:length(age_id)){
    for (j in 1:7){
      for (k in 1:11){
        mResTab[,i,j,k]<-mean(na.omit(ResTab[,i,j,k]))
        mResTabus[,i,j,k]<-mean(na.omit(ResTabus[,i,j,k]))
        mResTabfb[,i,j,k]<-mean(na.omit(ResTabfb[,i,j,k]))
      } } }


  ResTabC <- list(mResTab,mResTabus,mResTabfb)


  ### FOR EACH OUTCOME, EST MEAN AND POSTERIOR INTERVAL FOR SPECIFIC YEARS (2025, 2050, 2100), PLUS PERCENT REDUCTION COMPARED TO BASE CASE
  ### long format
  #1 outcome  5
  #2 scenario  8
  #3 population 3
  #4 age 4  **** nothing so far
  #5 comparator -- 3 (abs, 2016 BC, current BC)
  #6 time points 4
  #7 statistic


  CatList <- list()
  CatList[[1]] <- c("ltbi_000s","pct_ltbi","tb_incidence_000s","tb_incidence_per_mil","tb_mortality_000s","tb_mortality_per_mil")
  CatList[[2]] <- c("base_case",paste("intervention_",1:5,sep=""),paste("scenario_",1:3,sep=""))
  CatList[[3]] <- c("all_populations","usb_population","fb_population")
  CatList[[4]] <- c("0-4",paste(0:8*10+5,1:9*10+4,sep="-"),"95+")
  CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
  CatList[[6]] <- 2018:2049
  CatList[[7]] <- c("mean","ci_low","ci_high")

  res_tab2 <-  cbind(expand.grid(CatList),NA)
  colnames(res_tab2) <- c("outcome","scenario","population","age_group","comparator","year","statistic","value")


  for (it in 1:nrow(res_tab2)) {
    i1 <- which(CatList[[1]] == res_tab2[it,'outcome'])
    i2 <- which(CatList[[2]] == res_tab2[it, 'scenario'])
    i3 <- which(CatList[[3]] == res_tab2[it, 'population'])
    i4 <- which(CatList[[4]] == res_tab2[it, 'age_group'])
    i5 <- which(CatList[[5]] == res_tab2[it, 'comparator'])
    i6 <- which(CatList[[6]] == res_tab2[it, 'year'])
    i7 <- which(CatList[[7]] == res_tab2[it, 'statistic'])

    res_tab2[it, 'value'] <-
      switch(
        as.character(res_tab2[it, 'comparator']),
        'absolute_value' = {
          ResTabC[[i3]][i2, i6, i1 + 1, i4]
        },
        'pct_basecase_same_year' = {
          ResTabC[[i3]][i2, i6, i1 + 1, i4] /
            ResTabC[[i3]][1, i6, i1 + 1, i4] * 100
        },
        'pct_basecase_2016' = {
          ResTabC[[i3]][i2, i6, i1 + 1, i4] /
            ResTabC[[i3]][1, i6, 1, i4] * 100
        }
      )
  }


  # ii5_1 <- res_tab2[,5]==CatList[[5]][1]
  # ii5_2 <- res_tab2[,5]==CatList[[5]][2]
  # ii5_3 <- res_tab2[,5]==CatList[[5]][3]
  #
  # for(i1 in 1:length(CatList[[1]])) { # i1=i3=i4=i5=i6=i7=1
  #  ii1<- res_tab2[,1]==CatList[[1]][i1]
  #   for(i2 in 1:length(CatList[[2]])) { # i1=i3=i4=i5=i6=i7=1
  #     ii2<- res_tab2[,2]==CatList[[2]][i2]
  #   for(i3 in 1:length(CatList[[3]])) { #
  #    ii3<- res_tab2[,3]==CatList[[3]][i3]
  #     for(i4 in 1:length(CatList[[4]])) { #
  #       ii4<-res_tab2[,4]==CatList[[4]][i4]
  #       for(i6 in 1:length(CatList[[6]])) { #
  #         ii6<-res_tab2[,6]==CatList[[6]][i6]
  #         for(i7 in 1:length(CatList[[7]])) { #
  #          ii7<- res_tab2[,7]==CatList[[7]][i7]
  # ####Okay so indexing of outputs
  # rw = ii1&ii2&ii3&ii4&ii6&ii7
  #
  # res_tab2[rw & ii5_1,"value"] = ResTabC[[i3]][i2,i6,i1+1,i4]
  #
  # res_tab2[rw & ii5_2,"value"] = ResTabC[[i3]][i2,i6,i1+1,i4]/
  #                                                            ResTabC[[i3]][1,i6,i1+1,i4]*100
  #
  # res_tab2[rw & ii5_3,"value"] = ResTabC[[i3]][i2,i6,i1+1,i4]/
  #                                                            ResTabC[[i3]][1,i6,1,i4]*100
  #
  #
  #         }}}}}}

  save(res_tab2, file=paste0("~/MITUS/",loc,"_sm_restab2.rda"))
}


#'@name tabby_results_b
#'@description reshapes output from OutputsInt to shape needed for Tabby2
#'@param loc two digit postal abbreviation for location of interest
#'@return res_tab2 large dataframe of results
#'@export
tabby_results_b<-function(loc){
  age_id = (2018:2049)-1949
  tt<-0
  #create 3 lists to hold output
  mResTabfb <- mResTabus <- mResTab <- array(NA,dim=c(9,length(age_id),7,3))

  ResTabfb <- ResTabus <- ResTab <- array(NA,dim=c(90,length(age_id),7,3))
  for (intv in 1:9){
    #load the results for all the runs (need to make this dataset)
    load(paste0("~/MITUS/",loc,"_results_",intv,".rda"))
    #dimensions of restab are:
    #number of datasets, number of ages, number of interventions
    nr<-10
    #define o as the results
    o<-out
    #gather the outputs from that model run
    for (b_ag in 1:3){
      ################################################################################
      #total population
      #dimensions are
      #scenarios; length age id;
      # print(1:nr+((intv-1)*10))
      ag_vec<- switch(b_ag, c(1,2,3),c(4,5,6,7),c(8,9,10,11))
      ################################################################################
      for (ag in ag_vec){
        ResTab[1:nr+((intv-1)*10),,1,b_ag]<-o[,age_id,1]
        #number of ltbi prevalence
        ResTab[1:nr+((intv-1)*10),,2,b_ag]<-sum(apply(o[,age_id,rowSums(expand.grid(c(54,65), ag_vec))],c(1,2),sum)*1e3)
        #percentage of ltbi prevalence
        ResTab[1:nr+((intv-1)*10),,3,b_ag]<-sum(apply(o[,age_id,rowSums(expand.grid(c(54,65), ag_vec))],c(1,2),sum))/sum(o[, age_id,rowSums(expand.grid(2, ag_vec))])*1e2
        #TB notifications (alive+dead at diagnosis)
        ResTab[1:nr+((intv-1)*10),,4,b_ag]<-sum(apply(o[,age_id,rowSums(expand.grid(c(135,188), ag_vec))],c(1,2),sum)*1e3)
        #percentage TB notifications (alive+dead at diagnosis)
        ResTab[1:nr+((intv-1)*10),,5,b_ag]<-sum(apply(o[,age_id,rowSums(expand.grid(c(135,188), ag_vec))],c(1,2),sum))/sum(o[, age_id,rowSums(expand.grid(2, ag_vec))])*1e2
        #tb attributable deaths
        ResTab[1:nr+((intv-1)*10),,6,b_ag]<-sum(apply(o[,age_id,rowSums(expand.grid(c(87,98), ag_vec))],c(1,2),sum)*1e3)
        # percentage tb attributable deaths
        ResTab[1:nr+((intv-1)*10),,7,b_ag]<-sum(apply(o[,age_id,rowSums(expand.grid(c(87,98), ag_vec))],c(1,2),sum))/sum(o[, age_id,rowSums(expand.grid(2, ag_vec))])*1e2

        ################################################################################
        #US Born population
        ################################################################################
        ResTabus[1:nr+((intv-1)*10),,1,b_ag]<-o[,age_id,1]
        #number of ltbi prevalence
        ResTabus[1:nr+((intv-1)*10),,2,b_ag]<-sum(o[, age_id,54+ag_vec]*1e3)
        #percentage of ltbi prevalence
        ResTabus[1:nr+((intv-1)*10),,3,b_ag]<-sum(o[, age_id,54+ag_vec])/sum(o[, age_id,32+ag_vec])*1e2
        #TB notifications (alive+dead at diagnosis)
        ResTabus[1:nr+((intv-1)*10),,4,b_ag]<-sum(apply(o[, age_id,rowSums(expand.grid(c(204,215), ag_vec))],c(1,2),sum)*1e3)
        #percentage TB notifications (alive+dead at diagnosis)
        ResTabus[1:nr+((intv-1)*10),,5,b_ag]<-sum(apply(o[, age_id,rowSums(expand.grid(c(204,215), ag_vec))],c(1,2),sum))/sum(o[, age_id,32+ag_vec])*1e2
        #tb attributable deaths
        ResTabus[1:nr+((intv-1)*10),,6,b_ag]<-sum(o[, age_id,87+ag_vec])*1e3
        # percentage tb attributable deaths
        ResTabus[1:nr+((intv-1)*10),,7,b_ag]<-sum(o[, age_id,87+ag_vec])/sum(o[, age_id,32+ag_vec])*1e2

        ################################################################################
        #non-US Born population
        ################################################################################
        ResTabfb[1:nr+((intv-1)*10),,1,b_ag]<-o[,age_id,1]
        #number of ltbi prevalence
        ResTabfb[1:nr+((intv-1)*10),,2,b_ag]<-sum(o[, age_id,65+ag_vec]*1e3)
        #percentage of ltbi prevalence
        ResTabfb[1:nr+((intv-1)*10),,3,b_ag]<-sum(o[, age_id,65+ag_vec])/sum(o[, age_id, 43+ag_vec])*1e2
        #TB notifications (alive+dead at diagnosis)
        ResTabfb[1:nr+((intv-1)*10),,4,b_ag]<-sum(apply(o[, age_id,rowSums(expand.grid(c(135,188), ag_vec))],c(1,2),sum)*1e3)
        #percentage TB notifications (alive+dead at diagnosis)
        ResTabfb[1:nr+((intv-1)*10),,5,b_ag]<-sum(apply(o[, age_id,rowSums(expand.grid(c(135,188), ag_vec))],c(1,2),sum))/sum(o[, age_id,43+ag_vec])*1e2
        #tb attributable deaths
        ResTabfb[1:nr+((intv-1)*10),,6,b_ag]<-sum(o[, age_id,98+ag_vec])*1e3
        # percentage tb attributable deaths
        ResTabfb[1:nr+((intv-1)*10),,7,b_ag]<-sum(o[, age_id,98+ag_vec])/sum(o[, age_id, 43+ag_vec])*1e2

      }
    }
    # print(paste(aa,"--",bb)); flush.console()
  } #end batch loop
  # ResTabAll[[tt]]   <- ResTab
  # ResTabAllfb[[tt]] <- ResTabfb
  # ResTabAllus[[tt]] <- ResTabus
  # }
  #concatenate all lists

  ################################################################################
  #take the mean for each of these
  for (i in 1:length(age_id)){
    for (j in 1:7){
      for (k in 1:11){
        mResTab[,i,j,k]<-mean(na.omit(ResTab[,i,j,k]))
        mResTabus[,i,j,k]<-mean(na.omit(ResTabus[,i,j,k]))
        mResTabfb[,i,j,k]<-mean(na.omit(ResTabfb[,i,j,k]))
      } } }


  ResTabC <- list(mResTab,mResTabus,mResTabfb)


  ### FOR EACH OUTCOME, EST MEAN AND POSTERIOR INTERVAL FOR SPECIFIC YEARS (2025, 2050, 2100), PLUS PERCENT REDUCTION COMPARED TO BASE CASE
  ### long format
  #1 outcome  5
  #2 scenario  8
  #3 population 3
  #4 age 4  **** nothing so far
  #5 comparator -- 3 (abs, 2016 BC, current BC)
  #6 time points 4
  #7 statistic


  CatList <- list()
  CatList[[1]] <- c("ltbi_000s","pct_ltbi","tb_incidence_000s","tb_incidence_per_mil","tb_mortality_000s","tb_mortality_per_mil")
  CatList[[2]] <- c("base_case",paste("intervention_",1:5,sep=""),paste("scenario_",1:3,sep=""))
  CatList[[3]] <- c("all_populations","usb_population","fb_population")
  CatList[[4]] <- c("0-24","25-64","65+")
  CatList[[5]] <- c("absolute_value","pct_basecase_same_year","pct_basecase_2016")
  CatList[[6]] <- 2018:2049
  CatList[[7]] <- c("mean","ci_low","ci_high")

  res_tab2 <-  cbind(expand.grid(CatList),NA)
  colnames(res_tab2) <- c("outcome","scenario","population","age_group","comparator","year","statistic","value")


  for (it in 1:nrow(res_tab2)) {
    i1 <- which(CatList[[1]] == res_tab2[it,'outcome'])
    i2 <- which(CatList[[2]] == res_tab2[it, 'scenario'])
    i3 <- which(CatList[[3]] == res_tab2[it, 'population'])
    i4 <- which(CatList[[4]] == res_tab2[it, 'age_group'])
    i5 <- which(CatList[[5]] == res_tab2[it, 'comparator'])
    i6 <- which(CatList[[6]] == res_tab2[it, 'year'])
    i7 <- which(CatList[[7]] == res_tab2[it, 'statistic'])

    res_tab2[it, 'value'] <-
      switch(
        as.character(res_tab2[it, 'comparator']),
        'absolute_value' = {
          ResTabC[[i3]][i2, i6, i1 + 1, i4]
        },
        'pct_basecase_same_year' = {
          ResTabC[[i3]][i2, i6, i1 + 1, i4] /
            ResTabC[[i3]][1, i6, i1 + 1, i4] * 100
        },
        'pct_basecase_2016' = {
          ResTabC[[i3]][i2, i6, i1 + 1, i4] /
            ResTabC[[i3]][1, i6, 1, i4] * 100
        }
      )
  }


  # ii5_1 <- res_tab2[,5]==CatList[[5]][1]
  # ii5_2 <- res_tab2[,5]==CatList[[5]][2]
  # ii5_3 <- res_tab2[,5]==CatList[[5]][3]
  #
  # for(i1 in 1:length(CatList[[1]])) { # i1=i3=i4=i5=i6=i7=1
  #  ii1<- res_tab2[,1]==CatList[[1]][i1]
  #   for(i2 in 1:length(CatList[[2]])) { # i1=i3=i4=i5=i6=i7=1
  #     ii2<- res_tab2[,2]==CatList[[2]][i2]
  #   for(i3 in 1:length(CatList[[3]])) { #
  #    ii3<- res_tab2[,3]==CatList[[3]][i3]
  #     for(i4 in 1:length(CatList[[4]])) { #
  #       ii4<-res_tab2[,4]==CatList[[4]][i4]
  #       for(i6 in 1:length(CatList[[6]])) { #
  #         ii6<-res_tab2[,6]==CatList[[6]][i6]
  #         for(i7 in 1:length(CatList[[7]])) { #
  #          ii7<- res_tab2[,7]==CatList[[7]][i7]
  # ####Okay so indexing of outputs
  # rw = ii1&ii2&ii3&ii4&ii6&ii7
  #
  # res_tab2[rw & ii5_1,"value"] = ResTabC[[i3]][i2,i6,i1+1,i4]
  #
  # res_tab2[rw & ii5_2,"value"] = ResTabC[[i3]][i2,i6,i1+1,i4]/
  #                                                            ResTabC[[i3]][1,i6,i1+1,i4]*100
  #
  # res_tab2[rw & ii5_3,"value"] = ResTabC[[i3]][i2,i6,i1+1,i4]/
  #                                                            ResTabC[[i3]][1,i6,1,i4]*100
  #
  #
  #         }}}}}}

  save(res_tab2, file=paste0("~/MITUS/",loc,"_bg_restab2.rda"))
}
