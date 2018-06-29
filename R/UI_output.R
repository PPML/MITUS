#' #'This is a script that will return the values of key calibration
#' #'targets that need to be shown to the user at the beginning of the
#' #'online Tabby II interface experience.
#' #'
#' #'It will grab values from the calibration target input file and
#' #'save them to an easily displayed value.
#' #'
#'
#' #'@param loc state for the model run
#' #'@param yrs years for which to express the historic screening estimates
#' #'@return one dimensional vector of historic screening estimates
#' #'@export
#'
#' ltscrt <- function (loc, yrs){
#'   file_name <- paste(loc,"_parAll_fin",".rData", sep = "")
#'   load(file_name)
#'   ltscrt <- rep(NA,(length(yrs)))
#' #'Because rLtScrt is on a monthly timescale, this function returns the
#' #'value at the six month of each year.
#'   for (i in length(yrs)){
#'     ltscrt[i] <- rLtScrt[(6*(yrs[i]-1949))]
#'   }
#'   return(ltscrt)
#' }
#'
#' #'This function will return the population size for the general
#' #'age and nativity group that the user will select their target
#' #'population from.
#' #'@param age_grp
#' #'@param nat_grp
#' #'@param yrs years for which to express the historic screening estimates
#' #'@return numeric vector of length yrs
#' #'@export
#'
#' pop_size <- function(age_grp, nat_grp, yrs) {
#'
#' }
#'
#'
#' #'This function will return TB incidence for a certain custom
#' #'intervention group specification over a certain number of years.
#' #'@param RR_tbprog
#' #'@param RR_ltbiprev
#' #'@param RR_mu
#' #'@param yrs years for which to express the historic screening estimates
#' #'@return numeric vector of length yrs
#' #'@export
#'
#' tb_incid <- function(RR_tbprog,RR_ltbiprev,RR_mu, yrs){
#'
#' }
#'
#' #'This function will return LTBI prevalence for a certain custom
#' #'intervention group specification over a certain number of years.
#' #'have the prevalance across each.
#' #'Transformation of the model runs that have been.
#' #'do we want an end yr, start year, rate per year for the intervention
#' #'implications of this, add a dashed line where intervention starts and stops.
#' #'@param RR_tbprog
#' #'@param RR_ltbiprev
#' #'@param RR_mu
#' #'@param yrs years for which to express the historic screening estimates
#' #'@return numeric vector of length yrs
#' #'@export
#' ltbi_prev <- function(RR_tbprog,RR_ltbiprev,RR_mu, yrs){
#'
#' }
