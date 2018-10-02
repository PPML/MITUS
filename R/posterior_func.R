#' #'This function calculates the posterior using the
#' #'lprior and llikelihood functions.
#' #'@name posterior
#' #'@param theta vector of values
#' #'@param n_cores number of cores to use for the optimization
#' #'@return posterior
#' #'@export
#' posterior <- function(theta) {
#'   -lprior(theta) - llikelihood(theta,n_cores)
#' }
