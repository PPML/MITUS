 ########## Optimize  for each country x5

setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC Large Grant/Analysis Transmission")

ScriptMak <- function(b) {
  
  ################################################
  Script <- paste('
  
  b    <-  ',b,'   # batch number

#########################  SET-UP  #########################  
## Parameter stuff
load("ParamInitUS_V738tab.rData") # ParamInit
P  <- ParamInit[,1]; names(P) <- rownames(ParamInit)
ii <-  ParamInit[,5]==1
ParamInitZ <- ParamInit[ParamInit$Calib==1,]
idZ0 <- ParamInitZ[,4]==0
idZ1 <- ParamInitZ[,4]==1
idZ2 <- ParamInitZ[,4]==2

## Scripts and functions
load("ModelInputs_9-2-16.rData")
source("ParamUS_V504.r")
source("PriorFuncUS_V22.r")
source("CalibFunctionsUS_V22.r")
source("TimeStepUS_V891.r")
source("IMISfunctionsUS_V236.r")
 source("BetterImisUS_V7.r")


load("Opt_US540_r8_4_1-27-16.rData")


  imis_result <- Imis2(B=2000,B.re=5000,number_k = 200,seed_val=b,SaveNam=paste(b,"_8",sep=""),n_cores=1,OptInit=o8$par) 

  save(imis_result,file=paste("imis_result_US8_",b,"_9-7-2016.rData",sep=""))

############ Done!   #####################################
  warnings()

  quit("no")  # Kill
  
                    ',sep=""); return(Script)  }

####################################

for (b in 1:200) {  kkk <- ScriptMak(b); write(kkk,file=paste("ImisSampUS8_",b,".r",sep=""))  } 


