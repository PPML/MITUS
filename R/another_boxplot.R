# ###
# #state paramname value
# ##load in all the data values and then reassign their names
#
# nat_hist<-c("pHR","RRmuHR","muIp","TunmuTbAg","pfast","ORpfast1","ORpfast2",
#             "ORpfastPI","ORpfastH","rfast","rslow","rslowH",
#             "TunrslowAge","rRecov", "rSlfCur", "SensLt","SpecLt",
#             "SensSp","EffLt")
# load("~/MITUS/Optim_all_10_2018-12-03.rda")
# new_opt_all<-melt(pnorm(opt_all[,-61]))
# ST_colnam<-colnames(opt_all[,-60])
#
# opt_all_11st<-new_opt_all
# load("~/Desktop/removed/Optim_all_10_2018-10-31.rda")
# US_colnam<-colnames(opt_all[,-64])
# US_opt_all<-melt(pnorm(opt_all[c(-1,-5),-64]))
# library(dplyr)
# param_plot1<-ggplot((data=filter(opt_all_11st, Var2 %in% unique(nat_hist)[1:5])),aes(Var2,value)) +geom_boxplot(aes(colour=Var1)) +geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[1:5])), aes(Var2,value))
# param_plot2<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(nat_hist)[6:10]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1)) +geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[6:10])), aes(Var2,value))
# param_plot3<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(nat_hist)[11:15]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[11:15])), aes(Var2,value))
# param_plot4<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(nat_hist)[16:19]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[16:19])), aes(Var2,value))
# #non-natural history params
# com_col_nam2<-intersect(US_colnam,ST_colnam)
# com_col_nam<-setdiff(com_col_nam2,"SensSn")
# param_plot5<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[1:5]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[1:5]),aes(Var2,value))
# param_plot6<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[6:10]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[6:10]),aes(Var2,value))
# param_plot7<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[11:15]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[11:15]),aes(Var2,value))
# param_plot8<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[16:20]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[16:20]),aes(Var2,value))
# param_plot9<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[21:25]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[21:25]),aes(Var2,value))
# param_plot10<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[26:30]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[26:30]),aes(Var2,value))
# param_plot9<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[31:35]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[31:35]),aes(Var2,value))
# param_plot10<-ggplot(data=filter(opt_all_11st, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[36:39]),aes(Var2,value)) +geom_boxplot(aes(colour=Var1))+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[36:39]),aes(Var2,value))
# pdf("~/MITUS/MITUS_results/new_US_opt_test.pdf", width=10.5,height=7)
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
