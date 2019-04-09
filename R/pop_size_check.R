#'
#'
#'
#' for (i in 1:nrow(df)){
#'   pdf(file=paste("MITUS_results/pop_size_",i,".pdf"), width = 11, height = 8.5)
#'   par(mfrow=c(2,2),mar=c(4,4.5,3,1))
#' # data("CalibDat_2018-07-12", package='MITUS') # CalibDat
#' #'Log-likelihood functions
#' #'Assign the calibration importance weights from CalibDat
#' #'These weights are based on year of the simulation.
#' # wts <- CalibDat[["ImptWeights"]]
#' #'format P
#' # data("ParamInitUS_2018-08-06_final", package='MITUS')# ParamInit
#' # P  <- ParamInit[,1];
#' # names(P) <- rownames(ParamInit)
#' # ii <-  ParamInit[,5]==1
#' # ParamInitZ <- ParamInit[ParamInit$Calib==1,]
#' # idZ0 <- ParamInitZ[,4]==0
#' # idZ1 <- ParamInitZ[,4]==1
#' # idZ2 <- ParamInitZ[,4]==2
#'
#'
#' names(Par) <- colnames(df[,-62])   ##previously, the distribution of parameters were transformed to normal distribution in
#' ##to facilitate comparisons. These first two steps convert these parameters back to their
#' ##distributions
#' # normal to uniform
#' Par2 <- pnorm(Par,0,1)
#' # uniform to true
#' Par3 <- Par2
#' Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
#' Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
#' Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
#' P[ii] <- Par3
#' P <- P
#'
#'
#'   prms <-list()
#'   prms <- param(P)
#'   IP <- list()
#'   IP <- param_init(P)
#'   trans_mat_tot_ages<<-reblncd(mubt = prms$mubt,can_go = can_go,RRmuHR = prms$RRmuHR[2], RRmuRF = prms$RRmuRF, HRdist = HRdist, dist_gen_v=dist_gen_v)
#'
#'   zz <- cSim(  nYrs       = 2018-1950         , nRes      = length(prms[["ResNam"]]), rDxt     = prms[["rDxt"]]  , TxQualt    = prms[["TxQualt"]]   , InitPop  = prms[["InitPop"]]    ,
#'                Mpfast     = prms[["Mpfast"]]    , ExogInf   = prms[["ExogInf"]]       , MpfastPI = prms[["MpfastPI"]], Mrslow     = prms[["Mrslow"]]    , rrSlowFB = prms[["rrSlowFB"]]    ,
#'                rfast      = prms[["rfast"]]     , RRcurDef  = prms[["RRcurDef"]]      , rSlfCur  = prms[["rSlfCur"]] , p_HR       = prms[["p_HR"]]      , dist_gen = prms[["dist_gen"]]    ,
#'                vTMort     = prms[["vTMort"]]    , RRmuRF    = prms[["RRmuRF"]]        , RRmuHR   = prms[["RRmuHR"]]  , Birthst  = prms[["Birthst"]]    ,
#'                HrEntEx    = prms[["HrEntEx"]]   , ImmNon    = prms[["ImmNon"]]        , ImmLat   = prms[["ImmLat" ]] , ImmAct     = prms[["ImmAct"]]    , ImmFst   = prms[["ImmFst" ]]    ,
#'                mubt       = prms[["mubt"]]      , RelInf    = prms[["RelInf"]]        , RelInfRg = prms[["RelInfRg"]], Vmix       = prms[["Vmix"]]      , rEmmigFB = prms [["rEmmigFB"]]  ,
#'                TxVec      = prms[["TxVec"]]     , TunTxMort = prms[["TunTxMort"]]     , rDeft    = prms[["rDeft"]]   , pReTx      = prms[["pReTx"]]     , LtTxPar  = prms[["LtTxPar"]]    ,
#'                LtDxPar    = prms[["LtDxPar"]]   , rLtScrt   = prms[["rLtScrt"]]       , RRdxAge  = prms[["RRdxAge"]] , rRecov     = prms[["rRecov"]]    , pImmScen = prms[["pImmScen"]]   ,
#'                EarlyTrend = prms[["EarlyTrend"]], NixTrans = IP[["NixTrans"]],   trans_mat_tot_ages = trans_mat_tot_ages)
#'
#'     M <- zz$Outputs
#'     colnames(M) <- prms[["ResNam"]]
#'     results<-as.data.frame(M)
#'
#'     V  <- cbind(results[1:66,30], results[1:66,31]+results[1:66,32])
#'     plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
#'     axis(1);axis(2,las=2);box()
#'     abline(h=axTicks(2),col="grey85")
#'     points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
#'     points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
#'     points(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
#'     lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,2],lty=3,col="grey50")
#'     lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,3],lty=3,col="blue")
#'     lines(CalibDat$tot_pop_yr_fb[,1],CalibDat$tot_pop_yr_fb[,4],lty=3,col="red3")
#'     lines(1950:2015,V[,2],lwd=2,col="red3")
#'     lines(1950:2015,V[,1],lwd=2,col="blue")
#'     lines(1950:2015,rowSums(V),lwd=2,col="grey50")
#'
#'     mtext("Year",1,2.5,cex=0.9)
#'     mtext(paste("Population: Total, US, and Non-US Born (mil, log-scale)",i, sep="_"),3,.8,font=2,cex=0.8)
#'     legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=0.9,
#'            pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
#'
#'     ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#'     ### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
#'     ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#'     V  <- cbind(rowSums(results[1:66,255:265]), rowSums(results[1:66,266:276]))
#'     V1 <- results[1:66,317]
#'     plot(1,1,ylim=c(0,max(range(V1))),xlim=c(1950,2015),xlab="",ylab="",axes=F)
#'     axis(1);axis(2,las=2);box()
#'     abline(h=axTicks(2),col="grey85")
#'
#'     lines(1950:2015,V[,2],lwd=2,col="red3")
#'     lines(1950:2015,V[,1],lwd=2,col="blue")
#'     lines(1950:2015,V1,lwd=2,col="grey50")
#'     points(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,pch=19,cex=0.6,col="grey50")
#'     lines(CalibDat$US_tot_mort[,1],(CalibDat$US_tot_mort[,2])/1e6,lty=3,col="grey50")
#'
#'     mtext("Year",1,2.5,cex=0.9)
#'     mtext("Mortality: Total, US, and Non-US Born",3,.8,font=2,cex=0.8)
#'     legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=0.9,
#'            pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
#'
#'     ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#'     ### ### ### ### ### ###   TOTAL MORT AGE DISTRIBUTION 2014  ### ### ### ### ### ###
#'     ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#' for (year in 1:65){
#'     V  <- cbind((results[year,255:265])+(results[year,266:276]))
#'     V1  <- V[,-3]
#'     V1[,2] <- V1[,2]+V[,3]
#'     V2 <- V1[,-4]
#'     V2[,3] <- V2[,3]+V1[,4]
#'     V3 <- V2[,-9]
#'     V3[,8] <- V3[,8]+V2[,9]
#'
#'
#'     plot(0,0,ylim=c(0.05,max(range(V3))),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA)
#'     axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
#'     axis(1,1:9-0.5,rep("",9))
#'     axis(2,c(0,.2,.4,.6,.8,1.0,1.2),las=2);box()
#'     abline(h=axTicks(2),col="grey85")
#'
#'     for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V3[1,i],V3[1,i]),border=NA,col="gray")
#'     # for(i in 1:8) points(i+.2,(CalibDat$US_mort_age[16,i+1])/1e6,pch=19,cex=1.2,col="black")
#'
#'
#'     mtext("Age Group",1,2.5,cex=0.9)
#'     box()
#'     mtext(paste("Mortality by Age,",year, "(mil)", sep=''),3,.8,font=2,cex=0.8)
#'     legend("topleft",c("Reported data","model"),pch=c(19,15),pt.cex=c(1,2),
#'            lwd=NA,col=c("black","gray"),bg="white")
#'     }
#'     dev.off()
#'
#'     }
#'
