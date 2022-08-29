#'@name  adj_immig_2020
#'@param TotImmAge
#'@param immig
#'@param return_months
#'@param multiplier
#'@return TotImmAge
#'@export
adj_immig_2020 <- function(TotImmAge,
                           immig = 99,
                           return_months =  865:888,
                           multiplier = 1
){
  if (immig != 99){
  #Calculate first month after return_months
  postMonth <- return_months[length(return_months)] + 1

  #Hold reduction through start of return
  TotImmAge[843:864,]<-TotImmAge[843:864,]-(TotImmAge[843:864,]*immig);
  for (agegrp in 1:ncol(TotImmAge)){
    # Bring up immigration to 50% by end of return (smoothly)
    TotImmAge[return_months,agegrp] <- seq(TotImmAge[864,agegrp],TotImmAge[842,agegrp]*multiplier, length.out=length(return_months))
    # Continue this trend from here to end of model run
    TotImmAge[postMonth:nrow(TotImmAge),agegrp] <- TotImmAge[postMonth:nrow(TotImmAge),agegrp] * multiplier
  }
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
                           return_params = def_returnScenario()
){

  #############################################################################
  ##### rDxt #####
  #############################################################################
  # Setup params
  rDxt_RM <- return_params[["rDxt"]][["return_months"]]
  rDxt_postMonth <- rDxt_RM[length(rDxt_RM)] + 1
  rDxt_mult <- return_params[["rDxt"]][["multiplier"]]

  # Initial covid effects
  rDxt[843:864,] <- rDxt[843:864,] - (rDxt[843:864,]*par2020["Dxt"])

  # Bring up params to 50% by end of return months (smoothly)
  for (riskgrp in 1:ncol(rDxt)){
    rDxt[rDxt_RM,riskgrp] <- seq(rDxt[864,riskgrp],rDxt[842,riskgrp]*rDxt_mult, length.out=length(rDxt_RM))
    rDxt[rDxt_postMonth:nrow(rDxt), riskgrp] <-  rDxt[rDxt_postMonth:nrow(rDxt), riskgrp]*rDxt_mult
  }

  #############################################################################
  ##### NixTrans #####
  #############################################################################

  # Setup params
  NixTrans_RM <- return_params[["Trans"]][["return_months"]]
  NixTrans_postMonth <- NixTrans_RM[length(NixTrans_RM)] + 1
  NixTrans_mult <- return_params[["Trans"]][["multiplier"]]


  NixTrans[843:864]<- (1-par2020["Trans"])
  NixTrans[NixTrans_RM] <- seq(NixTrans[864],NixTrans[842]*NixTrans_mult, length.out=length(NixTrans_RM))
  NixTrans[NixTrans_postMonth:length(NixTrans)] <- NixTrans[NixTrans_postMonth:length(NixTrans)] * NixTrans_mult

  #############################################################################
  ##### Changes to mortality with TB #####
  #############################################################################

  # Setup params
  CaseFat_RM <- return_params[["CaseFat"]][["return_months"]]
  CaseFat_postMonth <- CaseFat_RM[length(CaseFat_RM)] + 1
  CaseFat_mult <- return_params[["CaseFat"]][["multiplier"]]


  RRmuTBPand <- rep(1,1812)
  time_horizon <- 843:CaseFat_RM[length(CaseFat_RM)]
  RRmuTBPand[time_horizon[1]:1812] <-c(rep(par2020["CaseFat"], length(time_horizon) - length(CaseFat_RM)),
                                       seq(par2020["CaseFat"], 1*CaseFat_mult, length.out = length(CaseFat_RM)),
                                       rep(1*CaseFat_mult, length.out = length(CaseFat_postMonth:1812)))


  param_list <- list("rDxt" = rDxt,
                     "NixTrans" = NixTrans,
                     "RRmuTBPand" = RRmuTBPand)
  return(param_list)
}
