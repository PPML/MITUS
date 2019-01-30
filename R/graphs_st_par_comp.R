#'@name stcomp
#'@return pdf of comparisons of the state parameters and national; normalized.
#'@export
st_comp<-function(){
library(reshape2)
#state paramname value
##load in all the data values and then reassign their names

nat_hist<-c("RRmuHR","muIp","TunmuTbAg","pfast","ORpfast1","ORpfast2",
            "ORpfastPI","ORpfastH","rfast","rslow","rslowH",
            "TunrslowAge","rRecov", "rSlfCur", "SensLt","SpecLt",
            "SensSp","EffLt")

load("~/Desktop/13_thurs/CA_Optim_all_10_2018-12-14.rda")
CA_opt_all<-as.data.frame(CA_opt_all)
CA_opt_all[,ncol(CA_opt_all)]<-round(CA_opt_all[,ncol(CA_opt_all)],3)
m<-getmode(CA_opt_all[,ncol(CA_opt_all)])
mode_opt<-  as.matrix(CA_opt_all[CA_opt_all$post_val==m,][,-(ncol(CA_opt_all))])
CA_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/FL_Optim_all_10_2018-12-14.rda")
FL_opt_all<-as.data.frame(FL_opt_all)
FL_opt_all[,ncol(FL_opt_all)]<-round(FL_opt_all[,ncol(FL_opt_all)],3)
m<-getmode(FL_opt_all[,ncol(FL_opt_all)])
mode_opt<-  as.matrix(FL_opt_all[FL_opt_all$post_val==m,][,-(ncol(FL_opt_all))])
FL_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/GA_Optim_all_10_2018-12-14.rda")
GA_opt_all<-as.data.frame(GA_opt_all)
GA_opt_all[,ncol(GA_opt_all)]<-round(GA_opt_all[,ncol(GA_opt_all)],3)
m<-getmode(GA_opt_all[,ncol(GA_opt_all)])
mode_opt<-  as.matrix(GA_opt_all[GA_opt_all$post_val==m,][,-(ncol(GA_opt_all))])
GA_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/IL_Optim_all_10_2018-12-14.rda")
IL_opt_all<-as.data.frame(IL_opt_all)
IL_opt_all[,ncol(IL_opt_all)]<-round(IL_opt_all[,ncol(IL_opt_all)],3)
m<-getmode(IL_opt_all[,ncol(IL_opt_all)])
mode_opt<-  as.matrix(IL_opt_all[IL_opt_all$post_val==m,][,-(ncol(IL_opt_all))])
IL_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/MA_Optim_all_10_2018-12-14.rda")
MA_opt_all<-as.data.frame(MA_opt_all)
MA_opt_all[,ncol(MA_opt_all)]<-round(MA_opt_all[,ncol(MA_opt_all)],3)
m<-getmode(MA_opt_all[,ncol(MA_opt_all)])
mode_opt<-  as.matrix(MA_opt_all[MA_opt_all$post_val==m,][,-(ncol(MA_opt_all))])
MA_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/NY_Optim_all_10_2018-12-14.rda")
NY_opt_all<-as.data.frame(NY_opt_all)
NY_opt_all[,ncol(NY_opt_all)]<-round(NY_opt_all[,ncol(NY_opt_all)],3)
m<-getmode(NY_opt_all[,ncol(NY_opt_all)])
mode_opt<-  as.matrix(NY_opt_all[NY_opt_all$post_val==m,][,-(ncol(NY_opt_all))])
NY_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/NJ_Optim_all_10_2018-12-14.rda")
NJ_opt_all<-as.data.frame(NJ_opt_all)
NJ_opt_all[,ncol(NJ_opt_all)]<-round(NJ_opt_all[,ncol(NJ_opt_all)],3)
m<-getmode(NJ_opt_all[,ncol(NJ_opt_all)])
mode_opt<-  as.matrix(NJ_opt_all[NJ_opt_all$post_val==m,][,-(ncol(NJ_opt_all))])
NJ_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/PA_Optim_all_10_2018-12-14.rda")
PA_opt_all<-as.data.frame(PA_opt_all)
PA_opt_all[,ncol(PA_opt_all)]<-round(PA_opt_all[,ncol(PA_opt_all)],3)
m<-getmode(PA_opt_all[,ncol(PA_opt_all)])
mode_opt<-  as.matrix(PA_opt_all[PA_opt_all$post_val==m,][,-(ncol(PA_opt_all))])
PA_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/VA_Optim_all_10_2018-12-14.rda")
VA_opt_all<-as.data.frame(VA_opt_all)
VA_opt_all[,ncol(VA_opt_all)]<-round(VA_opt_all[,ncol(VA_opt_all)],3)
m<-getmode(VA_opt_all[,ncol(VA_opt_all)])
mode_opt<-  as.matrix(VA_opt_all[VA_opt_all$post_val==m,][,-(ncol(VA_opt_all))])
VA_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/WA_Optim_all_10_2018-12-14.rda")
WA_opt_all<-as.data.frame(WA_opt_all)
WA_opt_all[,ncol(WA_opt_all)]<-round(WA_opt_all[,ncol(WA_opt_all)],3)
m<-getmode(WA_opt_all[,ncol(WA_opt_all)])
mode_opt<-  as.matrix(WA_opt_all[WA_opt_all$post_val==m,][,-(ncol(WA_opt_all))])
WA_opt_all<-melt(pnorm(mode_opt))

load("~/Desktop/13_thurs/TX_Optim_all_10_2018-12-14.rda")
ST_colnam<-colnames(TX_opt_all[,-(ncol(TX_opt_all))])
TX_opt_all<-as.data.frame(TX_opt_all)
TX_opt_all[,ncol(TX_opt_all)]<-round(TX_opt_all[,ncol(TX_opt_all)],3)
m<-getmode(TX_opt_all[,ncol(TX_opt_all)])
mode_opt<-  as.matrix(TX_opt_all[TX_opt_all$post_val==m,][,-(ncol(TX_opt_all))])
TX_opt_all<-melt(pnorm(mode_opt))

states<-c(rep("CA",nrow(CA_opt_all)),rep("FL",nrow(FL_opt_all)),rep("GA",nrow(GA_opt_all)),rep("IL",nrow(IL_opt_all)),rep("MA",nrow(MA_opt_all)),rep("NJ",nrow(NJ_opt_all)),rep("NY",nrow(NY_opt_all)),rep("PA",nrow(PA_opt_all)),rep("TX",nrow(TX_opt_all)),rep("VA",nrow(VA_opt_all)),rep("WA",nrow(WA_opt_all)))
opt_all_11st<-as.data.frame(cbind(states,rbind(CA_opt_all,FL_opt_all,GA_opt_all,IL_opt_all,MA_opt_all,NJ_opt_all,NY_opt_all,PA_opt_all,TX_opt_all,VA_opt_all,WA_opt_all)))
US_opt_all<-readRDS("~/MITUS/inst/US/US_Optim_all_10_2019-01-25.rds")
US_colnam<-colnames(opt_all[,-(ncol(US_opt_all))])

US_opt_all<-as.data.frame(US_opt_all)
US_opt_all[,ncol(US_opt_all)]<-round(US_opt_all[,ncol(US_opt_all)],3)
m<-getmode(US_opt_all[,ncol(US_opt_all)])
mode_opt<-  as.matrix(US_opt_all[US_opt_all$post_val==m,][,-(ncol(US_opt_all))])
US_opt_all<-melt(pnorm(mode_opt))
library(dplyr)
library(ggplot2)
param_plot1<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(nat_hist)[1:5] & Var1 == "b_no_1")), aes(Var2,value)) +geom_point(aes(colour=states),position="jitter") +geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[1:5])), aes(Var2,value))
param_plot1<-param_plot1+annotate("text", x = 1, y = -.1, label = "0.05 [0.025,0.15]")+annotate("text", x = 2, y = -.1, label = "0.5 [0.2,0.8]")+annotate("text", x = 3, y = -.1, label = "3.4 [2.3,4.8]")
param_plot1<-param_plot1+annotate("text", x = 4, y = -.1, label = "0.005 [0.0025,0.01]")+annotate("text", x = 5, y = -.1, label = "0.064 [0.003,0.2]")+expand_limits(y=c(0,1))

param_plot2<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(nat_hist)[6:10]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter") +geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[6:10])), aes(Var2,value))
param_plot2<-param_plot2+annotate("text", x = 1, y = -.1, label = "2.5 [0.8, 5.1]")+annotate("text", x = 2, y = -.1, label = "0.5 [0.2, 1]")+annotate("text", x = 3, y = -.1, label = " 10 [4,19]")
param_plot2<-param_plot2+annotate("text", x = 4, y = -.1, label = "0.21 [0.08,0.4]")+annotate("text", x = 5, y = -.1, label = "0.55 [0.19,1.1]")+expand_limits(y=c(0,1))

param_plot3<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(nat_hist)[11:15] & Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[11:15])), aes(Var2,value))
param_plot3<-param_plot3+annotate("text", x = 1, y = -.1, label = "0.0007 [0.00004,0.002]")+annotate("text", x = 2, y = -.1, label = "0.042 [0.001,0.1]")+annotate("text", x = 3, y = -.1, label = "0.2 [0.1,0.35]")
param_plot3<-param_plot3+annotate("text", x = 4, y = -.1, label = "0.03 [0.01,0.06]")+annotate("text", x = 5, y = -.1, label = "0.2 [0.15,0.35]")+expand_limits(y=c(0,1))

param_plot4<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(nat_hist)[16:19]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5))+geom_point(data=(filter(US_opt_all, Var2 %in% unique(nat_hist)[16:19])), aes(Var2,value))
param_plot4<-param_plot4+annotate("text", x = 1, y = -.1, label = "0.8 [0.75,0.84]")+annotate("text", x = 2, y = -.1, label = "0.975 [.965,0.985]")+annotate("text", x = 3, y = -.1, label = "0.9 [0.8,0.95]")
param_plot4<-param_plot4+annotate("text", x = 4, y = -.1, label = "0.995 [0.99,0.999]")+expand_limits(y=c(0,1))


#non-natural history params
com_col_nam2<-intersect(US_colnam,ST_colnam)
com_col_nam<-setdiff(com_col_nam2,"SensSn")
param_plot5<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[1:5] & Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states), position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[1:5]),aes(Var2,value))
param_plot5<-param_plot5+annotate("text", x = 1, y = -.1, label = "1 [0.5,1.7]")+annotate("text", x = 2, y = -.1, label = "1 [0.4,1.9]")+annotate("text", x = 3, y = -.1, label = "1 [0.8,3.5]")
param_plot5<-param_plot5+annotate("text", x = 4, y = -.1, label = "0.2 [0.1,0.4]")+annotate("text", x = 5, y = -.1, label = "0.5 [0.3,0.9]")+expand_limits(y=c(0,1))

param_plot6<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[6:10] & Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[6:10]),aes(Var2,value))
param_plot6<-param_plot6+annotate("text", x = 1, y = -.1, label = "0.25 [0.125,0.5]")+annotate("text", x = 2, y = -.1, label = "0.6 [0.3,1.5]")+annotate("text", x = 3, y = -.1, label = "0.0001 [0.00001,0.0001]")
param_plot6<-param_plot6+annotate("text", x = 4, y = -.1, label = "0.05 [0.03,0.07]")+annotate("text", x = 5, y = -.1, label = "0.003 [0.002,0.004]")+expand_limits(y=c(0,1))

param_plot7<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[11:15]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[11:15]),aes(Var2,value))
param_plot7<-param_plot7+annotate("text", x = 1, y = -.1, label = "10 [5,15]")+annotate("text", x = 2, y = -.1, label = "0.5 [0.1,1]")+annotate("text", x = 3, y = -.1, label = "3 [2,5]")
param_plot7<-param_plot7+annotate("text", x = 4, y = -.1, label = "0.5 [0.1,0.9]")+annotate("text", x = 5, y = -.1, label = "0.5 [0.1,0.9]")+expand_limits(y=c(0,1))

param_plot8<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[16:20]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[16:20]),aes(Var2,value))
param_plot8<-param_plot8+annotate("text", x = 1, y = -.1, label = "1.2 [1,1.5]")+annotate("text", x = 2, y = -.1, label = "0.025 [0.01,0.1]")+annotate("text", x = 3, y = -.1, label = "3 [1.5,6.0]")
param_plot8<-param_plot8+annotate("text", x = 4, y = -.1, label = "0.2 [0.1,0.4]")+annotate("text", x = 5, y = -.1, label = "0.33 [0.29,0.37]")+expand_limits(y=c(0,1))

param_plot9<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[21:25]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[21:25]),aes(Var2,value))
param_plot9<-param_plot9+annotate("text", x = 1, y = -.1, label = "0.5 [0.25,0.75]")+annotate("text", x = 2, y = -.1, label = "0.6 [0.4,0.9]")+annotate("text", x = 3, y = -.1, label = "0 [-1,1]")
param_plot9<-param_plot9+annotate("text", x = 4, y = -.1, label = "0 [-1,1]")+annotate("text", x = 5, y = -.1, label = "0 [-1,1]")+expand_limits(y=c(0,1))

param_plot10<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[26:30]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[26:30]),aes(Var2,value))
param_plot10<-param_plot10+annotate("text", x = 1, y = -.1, label = "0 [-1,1]")+annotate("text", x = 2, y = -.1, label = "0 [-1,1]")+annotate("text", x = 3, y = -.1, label = "3 [2,4]")
param_plot10<-param_plot10+annotate("text", x = 4, y = -.1, label = "1 [0.5,2]")+annotate("text", x = 5, y = -.1, label = "0.8 [0.7,0.9]")+expand_limits(y=c(0,1))

param_plot11<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[31:35]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[31:35]),aes(Var2,value))
param_plot11<-param_plot11+annotate("text", x = 1, y = -.1, label = "0.5 [0.2,0.8]")+annotate("text", x = 2, y = -.1, label = "0.96 [0.94,0.98]")+annotate("text", x = 3, y = -.1, label = "0.1 [0.05,0.3]")
param_plot11<-param_plot11+annotate("text", x = 4, y = -.1, label = "1 [0.75,1.3]")+annotate("text", x = 5, y = -.1, label = "2 [1.2,3]")+expand_limits(y=c(0,1))

param_plot12<-ggplot(data=filter(opt_all_11st, (Var2 %in% unique(setdiff(com_col_nam,nat_hist))[36:39]& Var1 == "b_no_1")),aes(Var2,value)) +geom_point(aes(colour=states),position="jitter")+geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+geom_point(data=filter(US_opt_all, Var2 %in% unique(setdiff(com_col_nam,nat_hist))[36:39]),aes(Var2,value))
param_plot12<-param_plot12+annotate("text", x = 1, y = -.1, label = "0.85 [0.7,0.9]")+annotate("text", x = 2, y = -.1, label = "1 [0.75,1.3]")+annotate("text", x = 3, y = -.1, label = "0.25 [0.1,0.5]")
param_plot12<-param_plot12+annotate("text", x = 4, y = -.1, label = "0.5 [0.2,0.8]")+expand_limits(y=c(0,1))

pdf("~/MITUS/MITUS_results/st_new_opt_1_30_2019.pdf", width=10.5,height=7)
param_plot1
param_plot2
param_plot3
param_plot4
param_plot5
param_plot6
param_plot7
param_plot8
param_plot9
param_plot10
param_plot11
param_plot12
dev.off()

}
