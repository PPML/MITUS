# ### dont use this one use optim_data
# #state paramname value
# ##load in all the data values and then reassign their names
#
# nat_hist<-c("pHR","RRmuHR","muIp","TunmuTbAg","pfast","ORpfast1","ORpfast2",
#             "ORpfastPI","ORpfastH","rfast","rslow","rslowH",
#             "TunrslowAge","rRecov", "rSlfCur", "SensLt","SpecLt",
#             "SensSp","EffLt")
#
# load("CA_parAll_2018-12-03.rda")
# CA_par_all<-melt(ParMatrix[c(1,3,5,9,10),])
#
# load("FL_parAll_2018-12-03.rda")
# FL_par_all<-melt(ParMatrix[c(-7,-8),])
#
# load("GA_parAll_2018-12-03.rda")
# GA_par_all<-melt(ParMatrix[c(1,3,5,8),])
#
# load("IL_parAll_2018-12-03.rda")
# IL_par_all<-melt(ParMatrix[c(1,3,4,5,6,8,9,10),])
#
# load("MA_parAll_2018-12-03.rda")
# MA_par_all<-melt(ParMatrix[c(1,3,4,7,8,10),])
#
# load("NY_parAll_2018-12-03.rda")
# NY_par_all<-melt(ParMatrix[c(1,2,5,6,7,8,9),])
#
# load("NJ_parAll_2018-12-03.rda")
# NJ_par_all<-melt(ParMatrix[c(-6,-8),])
#
# load("PA_parAll_2018-12-03.rda")
# PA_par_all<-melt(ParMatrix[c(4,5,6,7,10),])
#
# load("VA_parAll_2018-12-03.rda")
# VA_par_all<-melt(ParMatrix[c(1,2,3,5,9,10),])
#
# load("WA_parAll_2018-12-03.rda")
# WA_par_all<-melt(ParMatrix[c(1,4,5,9),])
#
# load("TX_parAll_2018-12-03.rda")
# TX_par_all<-melt(ParMatrix[c(-5,-8),])
# ST_colnam<-colnames(ParMatrix[,])
#
# states<-c(rep("CA",nrow(CA_par_all)),rep("FL",nrow(FL_par_all)),rep("GA",nrow(GA_par_all)),rep("IL",nrow(IL_par_all)),rep("MA",nrow(MA_par_all)),rep("NJ",nrow(NJ_par_all)),rep("NY",nrow(NY_par_all)),rep("PA",nrow(PA_par_all)),rep("TX",nrow(TX_par_all)),rep("VA",nrow(VA_par_all)),rep("WA",nrow(WA_par_all)))
# par_all_11st<-as.data.frame(cbind(states,rbind(CA_par_all,FL_par_all,GA_par_all,IL_par_all,MA_par_all,NJ_par_all,NY_par_all,PA_par_all,TX_par_all,VA_par_all,WA_par_all)))
# load("US_parAll_2018-12-03.rda")
# US_colnam<-colnames(ParMatrix[,])
# US_par_all<-melt(ParMatrix[c(-1,-5),])
# library(dplyr)
# param_plot1<-ggplot((data=filter(par_all_11st, Var2 %in% unique(nat_hist)[1:5])),aes(Var2,value)) +geom_boxplot(aes(colour=states)) +geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_par_all, Var2 %in% unique(nat_hist)[1:5])), aes(Var2,value))
# param_plot2<-ggplot(data=filter(par_all_11st, Var2 %in% unique(nat_hist)[6:10]),aes(Var2,value)) +geom_boxplot(aes(colour=states)) +geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_par_all, Var2 %in% unique(nat_hist)[6:10])), aes(Var2,value))
# param_plot3<-ggplot(data=filter(par_all_11st, Var2 %in% unique(nat_hist)[11:15]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_par_all, Var2 %in% unique(nat_hist)[11:15])), aes(Var2,value))
# param_plot4<-ggplot(data=filter(par_all_11st, Var2 %in% unique(nat_hist)[16:19]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5))+geom_point(data=(filter(US_par_all, Var2 %in% unique(nat_hist)[16:19])), aes(Var2,value))
# #non-natural history params
# com_col_nam2<-intersect(US_colnam,ST_colnam)
# com_col_nam<-setdiff(com_col_nam2,"SensSn")
# param_plot5<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[1:5]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[1:5]),aes(Var2,value))
# param_plot6<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[6:10]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[6:10]),aes(Var2,value))
# param_plot7<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[11:15]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[11:15]),aes(Var2,value))
# param_plot8<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[16:20]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[16:20]),aes(Var2,value))
# param_plot9<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[21:25]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[21:25]),aes(Var2,value))
# param_plot10<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[26:30]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[26:30]),aes(Var2,value))
# param_plot9<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[31:35]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[31:35]),aes(Var2,value))
# param_plot10<-ggplot(data=filter(par_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[36:39]),aes(Var2,value)) +geom_boxplot(aes(colour=states))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_par_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[36:39]),aes(Var2,value))
# pdf("~/MITUS/MITUS_results/state_par_test", width=10.5,height=7)
# param_plot1
# param_plot2
# param_plot3
# param_plot4
# param_plot5
# param_plot6
# param_plot7
# param_plot8
# param_plot9
# param_plot10
# dev.off()
#
