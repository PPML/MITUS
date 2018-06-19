#'load initial conditions
#'load necessary libraries
library(mnormt)
library(parallel)
library(lhs)
library(Rcpp)
library(MCMCpack)
library(MASS)

#'lazy load necessary datasets
#'Model Input style
load("data/ModelInputs_9-2-16.rData")
