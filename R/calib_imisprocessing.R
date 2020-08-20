imis_processing<-function(){
###### GRAPHS COMPARING SIMS VS CALIB DATA ############
# a helper function for making transparent colors
mTrsp <- function(cl,a)  { apply(col2rgb(cl), 2, function(x){ rgb(x[1],x[2],x[3],a,maxColorValue=255)}) }
library(Hmisc)
# ###################### ###################### ######################
model_load("US")
# # SAMP ADDN PARS
iiB <-  rownames(ParamInit)%in%c("ImmigVolFut","pDefLt","EffLt")
ParamInitZB <- ParamInit[iiB,]
idZ0B <- ParamInitZB[,4]==0
idZ1B <- ParamInitZB[,4]==1
idZ2B <- ParamInitZB[,4]==2
par_sampB <- matrix(NA,1000,nrow(ParamInitZB))
for(y in 1:nrow(ParamInitZB)){ set.seed(123+y); par_sampB[,y] <- rnorm(1000)  }
#
# ## Scripts and functions
# load("ModelInputs_7-5-2019.rData")
# source("ParamUS2019_V45.r")
# source("PriorFuncUS_V22.r")
# source("CalibFunctions2019_V3.r")
# source("TimeStepUS_V94.r")
# source("IMISfunctions2019_V61.r")
# load("CalibDat_11-26-19.rData") # CalibDat
#
# ########  Load par set

load("/Users/nis100/Desktop/MetaTabby/imis_result_temp_061720.rData")
par_samp <- imis_res_tmp$resample[1:1000,]
###what are we doing here???
for(i in 2:5) { # i=2
  load(paste0("imis_result_US2019b_1_",i,"_11-26-2019.rData"))
  par_samp <- rbind(par_samp,imis_result$resample)
}
#
#
# ############ FUNCTION
#
pred_al_res <- function(Par,ParB) {
  Par <<- Par
  Par3 <- Par2 <- pnorm(Par,0,1) # norm2unif
  Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7]) # unif2true
  Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
  Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
  P[ii]      <- Par3

  ParB <<- ParB
  Par3B <- Par2B <- pnorm(ParB,0,1) # norm2unif
  Par3B[idZ0B] <- qbeta( Par2B[idZ0B], shape1  = ParamInitZB[idZ0B,6], shape2 = ParamInitZB[idZ0B,7]) # unif2true
  Par3B[idZ1B] <- qgamma(Par2B[idZ1B], shape   = ParamInitZB[idZ1B,6], rate   = ParamInitZB[idZ1B,7])
  Par3B[idZ2B] <- qnorm( Par2B[idZ2B], mean    = ParamInitZB[idZ2B,6], sd     = ParamInitZB[idZ2B,7])
  P[iiB]      <- Par3B

  P <<- P
  prms <-list()
  prms<-param_init(PV=P, loc=loc, prg_chng = def_prgchng(P), ttt_list = def_ttt())

  trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v, adj_fact=prms[["adj_fact"]])
  rownames(trans_mat_tot_ages) <-  paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4))
  colnames(trans_mat_tot_ages) <-  rep(paste0(rep(paste0("p",0:3),each=4),"_",rep(paste0("m",0:3),4)),11)

  if(any(trans_mat_tot_ages>1)) print("transition probabilities are too high")
  zz <- cSim( nYrs       = 2019-1950         , nRes      = length(func_ResNam())  , rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
              Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]  ,
              rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]] ,  dist_gen = prms[["dist_gen"]]    ,
              vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
              HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat"]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst"]]    ,
              net_mig_usb = prms[["net_mig_usb"]], net_mig_nusb = prms[["net_mig_nusb"]],
              mubt       = prms[["mubt"]]    , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
              TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
              LtDxPar_lt    = prms[["LtDxPar_lt"]]   , LtDxPar_nolt    = prms[["LtDxPar_nolt"]]   , rLtScrt   = prms[["rLtScrt"]]       , ttt_samp_dist   = prms[["ttt_sampling_dist"]] ,
              ttt_ag = prms[["ttt_ag"]], ttt_na = prms[["ttt_na"]], ttt_month = prms[["ttt_month"]], ttt_ltbi = prms[["ttt_ltbi"]], ttt_pop_scrn = prms[["ttt_pop_scrn"]], RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
              EarlyTrend = prms[["EarlyTrend"]], ag_den=prms[["aging_denom"]],  NixTrans = prms[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)

  M <- zz$Outputs
  colnames(M) <- prms[["ResNam"]]
  return(M)
}
#
# ######### SIM RESULTS
loc<-"US"
mean_res <-   pred_al_res(Par = Opt[8,-58], ParB = rep(0,3))
#
par_samp_unique <- unique(par_samp) # 533 long
#lets save the par_samp_unique for the metatabby runs
saveRDS(par_samp_unique,"~/Desktop/MetaTabby/UniqueIMISParams_061720.rds", version=2)
#
n_reps <- NULL
for(i in 1:dim(par_samp_unique)[1]) n_reps[i] <- sum(par_samp[,1]==par_samp_unique[i,1])
#what are the dimensions here??
sim_res<- array(NA,dim=c(70,length(func_ResNam()),dim(par_samp_unique)[1]))
dimnames(sim_res)[[2]] <- func_ResNam()
for(hh in 1:dim(par_samp_unique)[1]){
  sim_res[,,hh] <- pred_al_res(Par = par_samp_unique[hh,], ParB = par_sampB[hh,])
  cat('\r',hh,"     "); flush.console() }

save(mean_res,sim_res,n_reps,file="metamod_intervals_draws_Int1_6-23-2020.rData")


###############  Dataset with 2000 model draws, plus a mean vector
load("metamod_intervals_draws_6-23-2020_Int1.rData")
# sim_res = array of model simulations
#     dim 1: calendar year, from 1950
#     dim 2: outcome
#     dim 3: parameter set
# mean_res = matrix of model results from best-fitting parameer set
#     dim 1: calendar year, from 1950
#     dim 2: outcome

# n_reps = a vector of number of times each parameter set is included in the posterior sample

## Calculate some outcomes and intervals
# total cases 1950-2019, then by usb and nusb

# point estimate results
notif_tot_m   <- rowSums(mean_res[,c("NOTIF_ALL","NOTIF_MORT_ALL")])*1e6
notif_usb_m   <- rowSums(mean_res[,c("NOTIF_US","NOTIF_MORT_US")])*1e6
notif_nusb_m   <- rowSums(mean_res[,c("NOTIF_F1","NOTIF_F2","NOTIF_MORT_F1", "NOTIF_MORT_F2")])*1e6

# intervals
notif_tot_i   <- apply(apply(sim_res[,c("NOTIF_ALL","NOTIF_MORT_ALL"),],c(1,3),sum),1,function(x) wtd.quantile(x,n_reps,c(1,39)/40))*1e6
notif_usb_i   <- apply(apply(sim_res[,c("NOTIF_US","NOTIF_MORT_US"),],c(1,3),sum),1,function(x) wtd.quantile(x,n_reps,c(1,39)/40))*1e6
notif_nusb_i   <- apply(apply(sim_res[,c("NOTIF_F1","NOTIF_F2","NOTIF_MORT_F1", "NOTIF_MORT_F2"),],c(1,3),sum),1,function(x) wtd.quantile(x,n_reps,c(1,39)/40))*1e6

warnings()
dim(sim_res)

wtd.quantile(x, weights=NULL, probs=c(0, .25, .5, .75, 1),
             type=c('quantile','(i-1)/(n-1)','i/(n+1)','i/n'),
             normwt=FALSE, na.rm=TRUE)
## Look at them

plot(1,1,ylim=c(1,150),xlim=c(1950,2019),xlab="",ylab="",axes=F,log="y")

polygon(c(1950:2019,2019:1950),c(notif_tot_i[1,],notif_tot_i[2,70:1])/1e3,border=F,col=mTrsp(1,100))
polygon(c(1950:2019,2019:1950),c(notif_usb_i[1,],notif_usb_i[2,70:1])/1e3,border=F,col=mTrsp(4,100))
polygon(c(1950:2019,2019:1950),c(notif_nusb_i[1,],notif_nusb_i[2,70:1])/1e3,border=F,col=mTrsp("forestgreen",100))

lines(1950:2019,notif_tot_m/1e3,pch=21,lwd=2,col=1)
lines(1950:2019,notif_usb_m/1e3,pch=21,lwd=2,col=4)
lines(1950:2019,notif_nusb_m/1e3,pch=21,lwd=2,col="forestgreen")

axis(2,las=1,tcl=-.2,mgp=c(3, 0.5, 0),col="grey50")
axis(1,las=1,tcl=-.2,mgp=c(3, 0.4, 0),col="grey50");
mtext("Year",1,1.8,cex=0.9)
mtext("Annual TB cases (000s)",2,2,cex=.9)
box(col="grey50")


# Fit some initial functions to predict intervals
#  Below assumes the upper interval is calculated as exp(log(y)+a) and the
#  lower interval is exp(log(y)-a), where a is something we fit.

notif_tot_a1 <- log(notif_tot_i[2,]) - log(notif_tot_m)
notif_tot_a2 <- log(notif_tot_m) - log(notif_tot_i[1,])

# organize values for fitting model
obs_a <- c(notif_tot_a1,notif_tot_a2) # observed values of a
upper <- rep(1:0,each=70) # predictor = 1 if upper interval, 0 if lower interval
year  <- rep(1950:2019,2) # calendar year of estimate
point_est <- rep(notif_tot_m,2) # point estimate
dat1 <- data.frame(obs_a=obs_a, upper=upper, year=year, point_est=point_est)

# fit some models, check AIC
fit_11 <- lm(obs_a~1,data=dat1);
summary(fit_11); AIC(fit_11) # a = 0.139

fit_12 <- lm(obs_a~upper,data=dat1);
summary(fit_12) ; AIC(fit_12)

fit_13 <- lm(obs_a~upper+year,data=dat1)
summary(fit_13) ; AIC(fit_13)

fit_14 <- lm(obs_a~upper+year+I(year^2),data=dat1)
summary(fit_14) ; AIC(fit_14)

fit_15 <- lm(obs_a~upper+year+I(year^2)+point_est,data=dat1)
summary(fit_15) ; AIC(fit_15)

fit_16 <- lm(obs_a~upper+year+I(year^2)+I(year^3)+point_est,data=dat1)
summary(fit_16) ; AIC(fit_16)

# Look at first model --> not adequate!
preds <- predict(fit_11)
pred_intervals0 <- exp(log(point_est)+preds*rep(c(1,-1),each=70))
pred_intervals  <- rbind(pred_intervals0[-(1:70)],pred_intervals0[1:70])

plot(1,1,ylim=c(5,85),xlim=c(1950,2019),xlab="",ylab="",axes=F,log="y")

polygon(c(1950:2019,2019:1950),c(notif_tot_i[1,],notif_tot_i[2,70:1])/1e3,border=F,col=mTrsp(1,80))
polygon(c(1950:2019,2019:1950),c(pred_intervals[1,],pred_intervals[2,70:1])/1e3,border=F,col=mTrsp(4,80))

axis(2,las=1,tcl=-.2,mgp=c(3, 0.5, 0),col="grey50")
axis(1,las=1,tcl=-.2,mgp=c(3, 0.4, 0),col="grey50");
mtext("Year",1,1.8,cex=0.9)
mtext("Annual TB cases (000s)",2,2,cex=.9)
box(col="grey50")

# Look at last model --> better, though not so good for future years
preds <- predict(fit_16)
pred_intervals0 <- exp(log(point_est)+preds*rep(c(1,-1),each=70))
pred_intervals  <- rbind(pred_intervals0[-(1:70)],pred_intervals0[1:70])

plot(1,1,ylim=c(5,85),xlim=c(1950,2019),xlab="",ylab="",axes=F,log="y")

polygon(c(1950:2019,2019:1950),c(notif_tot_i[1,],notif_tot_i[2,70:1])/1e3,border=F,col=mTrsp(1,80))
polygon(c(1950:2019,2019:1950),c(pred_intervals[1,],pred_intervals[2,70:1])/1e3,border=F,col=mTrsp(4,80))

axis(2,las=1,tcl=-.2,mgp=c(3, 0.5, 0),col="grey50")
axis(1,las=1,tcl=-.2,mgp=c(3, 0.4, 0),col="grey50");
mtext("Year",1,1.8,cex=0.9)
mtext("Annual TB cases (000s)",2,2,cex=.9)
box(col="grey50")

#####################################################
}
