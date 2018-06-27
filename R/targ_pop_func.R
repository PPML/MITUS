#' #' creates target population based on user input into the model
#' #' This script takes input from the TABBY2 user interface and
#' #' creates the parameters needed to fit an invlogit 2 dimensional
#' #'normal distribution for the target populations.
#'
#' #'This code creates a vector of weights for each age group based
#' #'on user entered percentages of age groups targeted.
#' #'do i need this or will tabby generate
#'
#'
#'
#' #'This code creates a vector of parameters to feed into an
#' #'inverse logit 2 dimensional normal distribution that will
#' #'be used to subset a portion of the population into an altered
#' #'risk group profile.
#' #'@param pop_size integer value bounded by
#' #'@param nat integer value; 0=US, 1=nonUS, 2=all
#' #'@param age_dist
#' #'@param RR_mort
#' #'@param RR_prog
#' #'@param age_wgts
#'
#' #'@param
#' #'@return vector of parameters for the distribution
#' #'@return targ_perc numeric value of target population's percent of the total population
#' #'@export
#'
#' targ_dist <- function(pop_size, age_dist, nat, RR_mort, RR_prog, age_wgts) {
#'
#'
#'
#'   return(targ_parm)
#'   return(targ_perc)
#' }
