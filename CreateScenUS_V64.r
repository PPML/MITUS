########## Optimize  for each country x5

setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC Large Grant/Analysis Transmission")

ScriptMak <- function(b,bb) {
  
  ################################################
  Script <- paste('
  
  b    <-  ',b,'   # batch number
  bb    <-  ',bb,'   # batch number 2

############# RUN All
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
                  
 source("TimeStepUS_V892intExtra3.r")
 source("ProcessSims_V33intExtra3.r")
 load("parAll_9-14-16.rData")  # parAll
     
  ww <- rep(0,8); 
  if(b>0) { ww[b] <- 1 } ##  Which scenario

  nr <- if(bb==15) { nrow(parAll)-14000 } else { 1000 }
    MiAllex <- array(NA,dim=c(nr,150,390))
           
  for(iy in 1:nr) {
     MiAllex[iy,,] <- OutputsInt(ParMatrix=parAll[iy+(bb-1)*1000,-(1:2)],endyr=2100,Int1=ww[1],Int2=ww[2],Int3=ww[3],Int4=ww[4],Int5=ww[5],Scen1=ww[6],Scen2=ww[7],Scen3=ww[8])
     print(iy); flush.console() } 
  assign(paste("MiAllex",bb-1,sep=""),MiAllex)
  save(list=paste("MiAllex",bb-1,sep=""),file=paste("MiAllex",bb-1,"_",b,"_3-13-17.rData",sep=""))
          
############ Done!   #####################################
  warnings()

  quit("no")  # Kill
  
                    ',sep=""); return(Script)  }

####################################

for (b in 0:8) { for (bb in 1:15) {  kkk <- ScriptMak(b,bb); write(kkk,file=paste("RunScenUS64_",b,"_",bb,".r",sep=""))  }  }


