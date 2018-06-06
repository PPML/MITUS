#'This is a script that will return the values of key calibration
#'targets that need to be shown to the user at the beginning of the
#'online Tabby II interface experience.
#'
#'It will grab values from the calibration target input file and
#'save them to an easily displayed value.
#'



load_ltscrt <- function (loc){
  if (grep(loc)){
    load("data/CalibDat_9-14-16.rData")
    screening(length)
  }

  return(rLtScrt)
}


