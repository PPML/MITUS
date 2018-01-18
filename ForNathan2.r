
# setwd("/Users/nicolasmenzie/Google Drive/Harvard/CDC Large Grant/Analysis Transmission")


library(Hmisc)

## LOAD WEIGHTS
  load("parAll_9-14-16.rData") # parAll
  wt = parAll[,"wt"]

### LOAD INDICATOR LABELS, for reference
  load("ResNam.rData") # ResNam

## HELPER FUNCTIONS
# l takes a vector of results, and weights, and calculates weighted mean, lower bound, and upper bound for q=1,2,3 respectively
  l <- function(x,w,q) {
    if(q==1) {
      wtd.mean(x,w) 
    } else {
      if(q==2) { 
        wtd.quantile(x,w,1/40) 
      } else { 
        wtd.quantile(x,w,39/40) 
      }
    }
  }
  
## PULL OUT SPECIFIC RESULTS (5) FOR EACH SCENARIO (8), FOR DIFFERENT POPS (3)
  # OUTCOMES = 1. Incident M. tb Infection (per Mil)
   #           2. LTBI Prevalence (%)          
   #           3. New TB Cases (per Mil)        
   #           4. MDR-TB in New TB Cases (%)          
   #           5. TB-Related Deaths (per Mil)
   #           Z. ADD MORE AS DESIRED!
  
  # SCENARIOS = 1. Base case
  #             2. Intervention scenario 1 --TLTBI for new immigrants        
  #             3. Intervention scenario 2 --Improved TLTBI in US        
  #             4. Intervention scenario 3 --Better case detection        
  #             5. Intervention scenario 4 --Better TB treatment        
  #             6. Intervention scenario 5 --All improvements together        
  #             7. Sensitivity analysis  1 --no transmission within the United States after 2016        
  #             8. Sensitivity analysis  2 --no LTBI among all immigrants after 2016        

  # POPULATIONS = 1. ALL
  #               2. US-BORN     
  #               3. FOREIGN-BORN
  
    for(aa in 0:7) { # aa=7
    ResTabfb <- ResTabus <- ResTab <- array(NA,dim=c(nrow(parAll),90,6))
    for(bb in 0:14) { load(paste("cluster results/MiAllex",bb,"_",aa,"_9-21-16.rData",sep=""))
      nr <- if(bb==14) { nrow(parAll)-14000 } else { 1000 }
      o <- get(paste("MiAllex",bb,sep=""))
      ResTab[1:nr+bb*1000,,1] <- o[,-(1:60),1]
      ResTab[1:nr+bb*1000,,2] <- apply(o[,-(1:60),304:307],c(1,2),sum)/o[,-(1:60),2]*1e6
      ResTab[1:nr+bb*1000,,3] <- apply(o[,-(1:60),c(56:66,67:77)],c(1,2),sum)/o[,-(1:60),2]*1e2
      ResTab[1:nr+bb*1000,,4] <-  (o[,-(1:60),136]+o[,-(1:60),215])/o[,-(1:60),2]*1e6
      ResTab[1:nr+bb*1000,,5] <- apply(o[,-(1:60),c(159:160,164:165,169:170,174:175)],c(1,2),sum)/apply(o[,-(1:60),156:175],c(1,2),sum)*1e2
      ResTab[1:nr+bb*1000,,6] <- apply(o[,-(1:60),89:110],c(1,2),sum)/o[,-(1:60),2]*1e6
      
      ResTabus[1:nr+bb*1000,,1] <- o[,-(1:60),1]
      ResTabus[1:nr+bb*1000,,2] <- apply(o[,-(1:60),304:305],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
      ResTabus[1:nr+bb*1000,,3] <- apply(o[,-(1:60),56:66],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e2
      ResTabus[1:nr+bb*1000,,4] <-apply(o[,-(1:60),c(235:245,246:256)],c(1,2),sum)/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
      ResTabus[1:nr+bb*1000,,5] <- apply(o[,-(1:60),c(159:160,164:165)],c(1,2),sum)/apply(o[,-(1:60),156:165],c(1,2),sum)*1e2
      ResTabus[1:nr+bb*1000,,6] <- o[,-(1:60),308]/(o[,-(1:60),30]+o[,-(1:60),31])*1e6
      
      ResTabfb[1:nr+bb*1000,,1] <- o[,-(1:60),1]
      ResTabfb[1:nr+bb*1000,,2] <- apply(o[,-(1:60),306:307],c(1,2),sum)/(o[,-(1:60),32]+o[,-(1:60),33])*1e6
      ResTabfb[1:nr+bb*1000,,3] <- apply(o[,-(1:60),67:77],c(1,2),sum)/(o[,-(1:60),32]+o[,-(1:60),33])*1e2
      ResTabfb[1:nr+bb*1000,,4] <- (apply(o[,-(1:60),c(137:147,216:226)],c(1,2),sum)-apply(o[,-(1:60),c(235:245,246:256)],c(1,2),sum))/
        (o[,-(1:60),32]+o[,-(1:60),33])*1e6
      ResTabfb[1:nr+bb*1000,,5] <- apply(o[,-(1:60),c(169:170,174:175)],c(1,2),sum)/apply(o[,-(1:60),166:175],c(1,2),sum)*1e2
      ResTabfb[1:nr+bb*1000,,6] <- o[,-(1:60),309]/(o[,-(1:60),32]+o[,-(1:60),33])*1e6
      
      print(paste(aa,"--",bb)); flush.console()  }
    assign(paste("ResTab",aa,sep=""),ResTab)
    assign(paste("ResTabus",aa,sep=""),ResTabus)
    assign(paste("ResTabfb",aa,sep=""),ResTabfb)  }
  
  ## Put together as a list
  ResTabAllfb <- ResTabAllus <- ResTabAll <- list()
  for(aa in 0:7)  ResTabAll[[aa+1]] <- get(paste("ResTab",aa,sep=""))
  for(aa in 0:7)  ResTabAllfb[[aa+1]] <- get(paste("ResTabfb",aa,sep=""))
  for(aa in 0:7)  ResTabAllus[[aa+1]] <- get(paste("ResTabus",aa,sep=""))
  
  
  ### FOR EACH OUTCOME, EST MEAN AND POSTERIOR INTERVAL FOR SPECIFIC YEARS (2025, 2050, 2100), PLUS PERCENT REDUCTION COMPARED TO BASE CASE
  ResTab <- array(NA,dim=c(5*2,8,3,3))
  dimnames(ResTab)[[1]] <- c("Incident M. tb Infection (per Mil)","Percent of Base Case Value",
                             "LTBI Prevalence (%)","Percent of Base Case Value",
                             "New TB Cases (per Mil)","Percent of Base Case Value",
                             "MDR-TB in New TB Cases (%)","Percent of Base Case Value",
                             "TB-Related Deaths (per Mil)","Percent of Base Case Value")          
  dimnames(ResTab)[[2]] <- c("BaseCase",paste("Int",1:5,sep=""),paste("Scen",1:2,sep=""))
  dimnames(ResTab)[[3]] <- c(2025,2050,2100)
  dimnames(ResTab)[[4]] <- c("mean","ci_low","ci_high")
  ResTabfb <- ResTabus <- ResTab
  
  for(y0 in 1:3) {
    y = c(2025,2050,2099)[y0]-2009 
    for(oc in 1:5) {
      for(sc in 1:8) {
        for(st in 1:3) {
        ResTab[oc*2-1,sc,y0,st] <- l(ResTabAll[[sc]][,y,oc+1],wt,st)
        ResTab[oc*2  ,sc,y0,st] <- l(ResTabAll[[sc]][,y,oc+1]/ResTabAll[[1]][,y,oc+1],wt,st)*100
        ResTabus[oc*2-1,sc,y0,st] <- l(ResTabAllus[[sc]][,y,oc+1],wt,st)
        ResTabus[oc*2  ,sc,y0,st] <- l(ResTabAllus[[sc]][,y,oc+1]/ResTabAllus[[1]][,y,oc+1],wt,st)*100
        ResTabfb[oc*2-1,sc,y0,st] <- l(ResTabAllfb[[sc]][,y,oc+1],wt,st)
        ResTabfb[oc*2  ,sc,y0,st] <- l(ResTabAllfb[[sc]][,y,oc+1]/ResTabAllfb[[1]][,y,oc+1],wt,st)*100
      } } } }

    save(ResTab,ResTabus,ResTabfb,file="ResultsForNathan_1-23-17.rData")
  
## DONE