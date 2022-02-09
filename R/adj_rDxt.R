adj_rDxt <- function(rDxt){
  ### calculate the average rDxt from 2000 to 2020
  st_month <- (2010-1950)*12+1
  end_month <- (2021-1950)*12
  for (i in 1:ncol(rDxt)){
  av_rDxt <- mean(rDxt[(st_month:end_month),i])

  ### set 2020:2021
  rDxt[(end_month-11):(end_month+12),i] <- seq(rDxt[(end_month-11),i], av_rDxt, length.out=24)
  ### set 2022:2050
  rDxt[(end_month+13):nrow(rDxt),i] <- av_rDxt
  }
  return(rDxt)
}
