#' Function to estimate the defining parameters of a beta distribution
#' when given the distribution mean and variance
#' useful for distributing the initial population among the generic
#' risk groups
#'
#'if user can provide a risk or even risk ratio of the reactivation
#'risk and it's variance, we can distribute the population based on
#'the calculated alpha and beta below
beta_params <- function(mu, var) {
  a <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  b <- a * (1 / mu - 1)
  return(params = list(alpha = a, beta = b))
}

beta_dist  <-function(N, num_group, alpha, beta){
  pbeta((1:num_group)/num_group, alpha, beta)-pbeta((0:num_group-1)/num_group, alpha, beta)
  return
}
