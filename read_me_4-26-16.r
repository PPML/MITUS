#####  READ ME 4-26-16
# Description of files, functions etc
setwd("~/Google Drive/files_4-26-16")

#### Packages needed:
# mvtnorm -- standard MV Normal tools (density, quantile etc)
# mnormt -- this has a function for evaluating a MV normal density without any 
#    monte carlo error (this noise is poison for non-stochastic optimization)
# Rcpp -- create function in C++ 
# lhs -- for doing a latin hypercube sample
# parallel -- for creating a simple function for 
#    parallelization.. though I think this is now part of the base package
# MCMCpack -- used for some density used in the likelihood, can't remember which one right now
# RColorBrewer -- probably used in there somewhere

#### Scripts
#### GraphsVsCalibDat -- runs model, creates graphs to compare calib data vs post. mode. Code may need to be modified if non-Mac (last line)
####
#### Timestep -- this script sets up the model function (takes parameters, produces projections)
#     -- the "Extra" version has additional outputs produced, but otherwsie the same model
####
# MakePriorTab -- this produces several files used during optimization
#     a. estimates prior distributions to match a given mean and CI width
#     b. draws some starting values for optimization. These drawn from standard MV normal
####
# Param -- this loads up all data needed, takes vector of parameters and pre-processes then to use in the model
#     and creates some other objects which are used (e.g. vector of output names)
####
# IMIS functions -- creates functions used in optimization and IMIS (parameter vector is the subset which are included in calibration, not all)
#     a. llikelihood = parallelized log-likelihood of a given parameter vector (actually log-likelihood minus max log-likelihood)
#     b. lprior =  log-prior of a given parameter vector
#     c. sample.prior =  lhs sample from prior
# ProcessSims: function for taking vector/table of parameters and producing matrix/array of results (used after calibration)
# Calib functions - individual functions that make up the likelihood
# Better IMIS -- a version of the IMIS function from Raftery and Bao (IMIS package) that I have made small adjustments to 
# OptUniv -- small routine for taking the results from Optim and optimizing along each parameter individualy, 
#      just to make sure those parameters which have little effect are at their opt (or at least a quick hack towards it)
# PriorFunc -- calculates prior density (note -- parameters fitted in an unbounded space, MVN copula used to convert to parameter space)
# MainResults -- messy script I have been using to load and produce results
# CalibData   -- file where I have be preprocessing some inputs used for calibration
# ImmigrantInputs  -- preprocessing immigration data
# CreateOptimSimple -- creates scripts to send to cluster to identify posterior mode via optimization
# CreateOptim -- creates scripts to send to cluster to sample from posterior via IMIS
# CreateScen -- creates scripts to send to cluster to produce results for multiple scenarios


#### Files
# CalibDat -- list of data used for calibration
# M0i -- table of results from 100 draws from posterior, simulated to 2100 under status quo
# Mi -- table of results from 1000 draws from posterior, simulated to 2015 
# Mi_mean -- mean across Mi (post mean)
# ModelInputs -- as stated, inputs for model
# Opt_US48: results of optimization (i.e. posterior mode)
# ParamInit: parameter table for detals of prior
# parUnique1000 -- 1000 parameter sets draws from posterior


