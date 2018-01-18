
library(Hmisc)

## LOAD WEIGHTS
  load("parAll_9-14-16.rData") # parAll
  wt = parAll[,"wt"]

### LOAD INDICATOR LABELS, for reference
  load("ResNam.rData") # ResNam

## LOAD SIMULATION RESULTS
  for(rr in 0:14) { load(paste("MiAllex",rr,"_0_9-21-16.rData",sep="")); print(rr); flush.console() }
# MiAllexfoo (for foo in 0:14) has dim1 = sim (14245), dim2 = year (1950:2099), dim3 = indicator (ResNam = 309 long)

## HELPER FUNCTIONS
# gc pulls a vector or matrix from the simulation results
  gc <- function(yr,id) {
    out = NULL
    if(length(yr)==1 & length(id)==1) {
      for(bb in 0:14) out = c(out,get(paste("MiAllex",bb,sep=""))[,yr-1949,id])  }
    if( (length(yr)>1 & length(id)==1) | (length(yr)==1 & length(id)>1) ) {
      for(bb in 0:14) out = rbind(out,get(paste("MiAllex",bb,sep=""))[,yr-1949,id])  }
    if( length(yr)>1 & length(id)>1 ) {
      for(bb in 0:14) out = abind(out,get(paste("MiAllex",bb,sep=""))[,yr-1949,id],along=1)   }
    return(out)  }

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
  
## PULL OUT RESULTS
# Pull out total pop by age/year/nativity
  # yr=2000:2016, then 2016:2050
  # id = 34:44 for US born age groups from youngest to oldest (see ResNam[34:44])
  # id = 45:55 for foreign born age groups from youngest to oldest

  pop_0_4_US <-   gc(yr=2016,id=34)  # this is a vector of 14245 
  pop_0_4_FB <-   gc(yr=2016,id=45)  # this is a vector of 14245 
  
# Pull out LTBI pop by age/year/nativity
  # yr=2000:2016, then 2016:2050
  # id = 56:66 for US born age groups from youngest to oldest (see ResNam[56:66])
  # id = 67:77 for foreign born age groups from youngest to oldest
  
  ltbi_0_4_US <-   gc(yr=2016,id=56)  # this is a vector of 14245 
  ltbi_0_4_FB <-   gc(yr=2016,id=67)  # this is a vector of 14245 
  
# Pull out TB incidence by age/year/nativity
  # yr=2000:2016, then 2016:2050
  # id = 137:147 for USB+FB age groups from youngest to oldest ALIVE AT DIAGNOSIS
  # id = 216:226 for USB+FB age groups from youngest to oldest DEAD AT DIAGNOSIS
  # id = 235:245 for US born age groups from youngest to oldest ALIVE AT DIAGNOSIS
  # id = 246:256 for US born age groups from youngest to oldest DEAD AT DIAGNOSIS
  
  tb_0_4_US <-   gc(yr=2016,id=235)+gc(yr=2016,id=246)  # this is a vector of 14245 
  tb_0_4_FB <-   gc(yr=2016,id=137)+gc(yr=2016,id=216)-tb_0_4_US  # this is a vector of 14245 
  
### TO PLOT:
# histogram with US born on top, FB on bottom, agegroups on the horizontal
# seperate histogram for (1) total pop (2) total ltbi, (3) total TB cases, (4) LTBI prevalence (in %), (5) TB incidence (in cases per million)
  
  
  
  
