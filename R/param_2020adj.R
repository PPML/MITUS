#'@name  adj_immig_2020
#'@param TotImmAge
#'@param immig two letter abbreviation
#'@param return_months boolean for intervention 1
#'@param multiplier boolean for intervention 2
#'@return TotImmAge
#'@export
adj_immig_2020 <- function(TotImmAge,
                           immig = 99,
                           return_months =  865:888,
                           multiplier = 1
){
  #Hold reduction through December 2021
  TotImmAge[843:864,]<-TotImmAge[843:864,]-(TotImmAge[843:864,]*immig);
  # Bring up immigration to 50% by end of 2022 (smoothly)
  for (agegrp in 1:ncol(TotImmAge)){
    TotImmAge[return_months,agegrp] <- seq(TotImmAge[864,agegrp],TotImmAge[842,agegrp]*multiplier, length.out=length(return_months))
  }
  return(TotImmAge)
}

#'@name  adj_param_2020
#'@param rDxt
#'@param NixTrans
#'@param return_months boolean for intervention 1
#'@param multiplier boolean for intervention 2
#'@return TotImmAge
#'@export
adj_param_2020 <- function(rDxt,
                           NixTrans,
                           par2020,
                           return_months =  865:888,
                           multiplier = 1
){
  rDxt[843:864,]<-rDxt[843:864,] - (rDxt[843:864,]*par2020["Dxt"])
  NixTrans[843:864]<- (1-par2020["Trans"])
  # Bring up params to 50% by end of 2022 (smoothly)
  for (riskgrp in 1:ncol(rDxt)){
    rDxt[return_months,riskgrp] <- seq(rDxt[864,riskgrp],rDxt[842,riskgrp]*multiplier, length.out=length(return_months))
  }
  NixTrans[return_months] <- seq(NixTrans[864],NixTrans[842]*multiplier, length.out=length(return_months))

  RRmuTBPand <- rep(1,1812)
  time_horizon <- 843:return_months[length(return_months)]
  RRmuTBPand[time_horizon] <-c(rep(par2020["CaseFat"], length(time_horizon) - length(return_months)),
                               seq(par2020["CaseFat"], 1*multiplier, length.out = length(return_months)))

  param_list <- list("rDxt" = rDxt,
                     "NixTrans" = NixTrans,
                     "RRmuTBPand" = RRmuTBPand)
  return(param_list)
}
