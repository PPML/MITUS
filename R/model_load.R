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
load("~/MITUS/data/ModelInputs_9-2-16.rData")

