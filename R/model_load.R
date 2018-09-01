#'load initial conditions
#'load necessary libraries
library(mnormt)
library(parallel)
library(lhs)
library(Rcpp)
library(MCMCpack)
library(MASS)

#'lazy load necessary datasets
#'Model Input
data("ModelInputs_9-2-16", package = 'MITUS')
data("ParamInitUS_2018-08-06_final", package='MITUS')# ParamInit
data("StartVal_2018-08-30", package='MITUS')# ParamInit
data("CalibDat_2018-09-01", package='MITUS')# ParamInit

#'creation of background parameters
#'elements of P will be replaced from either the StartVals in the case
#'of optimization or the user inputted dataset
P  <- ParamInit[,1];
names(P) <- rownames(ParamInit)
ii <-  ParamInit[,5]==1
ParamInitZ <- ParamInit[ParamInit$Calib==1,]
idZ0 <- ParamInitZ[,4]==0
idZ1 <- ParamInitZ[,4]==1
idZ2 <- ParamInitZ[,4]==2

#'format the calibration data
wts <- CalibDat[["ImptWeights"]]

