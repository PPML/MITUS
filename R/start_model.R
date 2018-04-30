#'Prompts the user to input the necessary parameters for a target population
#'
#'This function interacts with the user and prompts them to enter the following
#'parameters:
#'
#'@param pop_size integer value bounded by
#'@param nat integer value; 0=US, 1=nonUS, 2=all
#'@param age_dist vector length 11 of age distribution of risk factors
#'                we will provide a default distribution if none is provided.
#'@param RR_mort double value of non-TB mortality risk ratio compared to general population levels
#'               conditional on age, nativity, RR of TB progression, living conditions, treatment and TB history.
#'@param RR_prog double value of TB progression risk ratio compared to general population levels conditional on
#'               age, nativity, RR of non-TB mortality, living conditions, treatment and TB history.
#'@param age_wgts vector length 11 of the percent of age group should make up the pop_size
#'@return target_params
#'@export

name <- function(variables) {
  if (pop_size > sum(v1)) stop("Your target population is greater than total population")

  if (pop_size > ) stop()

}




