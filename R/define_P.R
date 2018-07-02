###this script is to move this code into a place that I can find
### it until I can decide the best place to move it into for final
### model run

load("~/MITUS/data/ParamInit_2018.rData") # ParamInit
P  <- ParamInit[,1];
names(P) <- rownames(ParamInit)
ii <-  ParamInit[,5]==1
ParamInitZ <- ParamInit[ParamInit$Calib==1,]
idZ0 <- ParamInitZ[,4]==0
idZ1 <- ParamInitZ[,4]==1
idZ2 <- ParamInitZ[,4]==2

