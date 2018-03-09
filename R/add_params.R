#random variables from initial params that need to be sourced;
LgtCurve <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2100-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr)/12)
zz  <- as.numeric(EndVal)/(1+exp(-zz));  if(StYr>1950) { zz[1:((StYr-1950)*12)] <- 0 };    zz  }

NixTrans <- rep(1,1801)
#(c(1-LgtCurve(2016,2017,1)), rep(1,888))
#if(Scen1==0) {  NixTrans[] <- 1      }



dLtt          <-  rep(1/9,1801)
#(1-LgtCurve(2016,2021,0))

EffLt0        <- rep(1,1801)

#  LgtCurve(2016,2021,0)+1

#EffLtX        <- cbind(EffLt0,0,EffLt0,0,0)
