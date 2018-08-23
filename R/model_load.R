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
data("StartVal_2018-08-23", package='MITUS')# ParamInit
data("CalibDat_2018-07-12", package='MITUS')# ParamInit
