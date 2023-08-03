#'@name  adj_immig_2020
#'@param TotImmAge
#'@param immig
#'@param return_months
#'@param multiplier
#'@return TotImmAge
#'@export
adj_immig_2020 <- function(TotImmAge,
                           immig = 99,
                           return_months =  846:866,
                           multiplier = 1.33 # Default is from DHS OIS LPR Annual flow
                                             # https://www.dhs.gov/sites/default/files/2023-02/2022_0405_plcy_lawful_permanent_residents_fy2021v2.pdf
){
  if (immig != 99){
  #Calculate first month after return_months
  postMonth <- return_months[length(return_months)]
  # Calculate the last month of the initial pandemic effect
  lastMonth <- return_months[1] - 1
  #Hold reduction through start of return
  TotImmAge[843:lastMonth,] <- TotImmAge[843:lastMonth,] - (TotImmAge[843:lastMonth, ] * immig);
  for (agegrp in 1:ncol(TotImmAge)){
    # Bring up immigration to 50% by end of return (smoothly)
    # Month 888 = last month of 2022
    TotImmAge[return_months,agegrp] <- seq(TotImmAge[lastMonth,agegrp],TotImmAge[842,agegrp]*multiplier, length.out=length(return_months))
    # Continue this trend from here to end of model run
    TotImmAge[postMonth:nrow(TotImmAge),agegrp] <- TotImmAge[postMonth:nrow(TotImmAge),agegrp] * multiplier
  }
  }
  return(TotImmAge)
}

#'@name  adj_param_2020
#'@param rDxt
#'@param NixTrans
#'@param return_params
#'@return param_list
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
  rDxt_lastMonth <- rDxt_RM[1] - 1
  rDxt_mult <- return_params[["rDxt"]][["multiplier"]]

  # Initial covid effects
  rDxt[843:rDxt_lastMonth,] <- rDxt[843:rDxt_lastMonth,] - (rDxt[843:rDxt_lastMonth,] * par2020["Dxt"])

  # Bring up params to X% by end of return months (smoothly)
  for (riskgrp in 1:ncol(rDxt)){
    rDxt[rDxt_RM,riskgrp] <- seq(rDxt[843,riskgrp], rDxt[842,riskgrp]*rDxt_mult, length.out=length(rDxt_RM))
    rDxt[rDxt_postMonth:nrow(rDxt), riskgrp] <-  rDxt[rDxt_postMonth:nrow(rDxt), riskgrp]*rDxt_mult
  }

  #############################################################################
  ##### NixTrans #####
  #############################################################################

  # Setup params
  NixTrans_RM <- return_params[["Trans"]][["return_months"]]
  NixTrans_postMonth <- NixTrans_RM[length(NixTrans_RM)] + 1
  NixTrans_mult <- return_params[["Trans"]][["multiplier"]]
  NixTrans_lastMonth <- NixTrans_RM[1] - 1

  NixTrans[843:NixTrans_lastMonth]<- (1-par2020["Trans"])
  NixTrans[NixTrans_RM] <- seq(NixTrans[NixTrans_lastMonth],NixTrans[842]*NixTrans_mult, length.out=length(NixTrans_RM))
  NixTrans[NixTrans_postMonth:length(NixTrans)] <- NixTrans[NixTrans_postMonth:length(NixTrans)] * NixTrans_mult

  #############################################################################
  ##### Changes to mortality with TB #####
  #############################################################################

  # Setup params
  CaseFat_RM <- return_params[["CaseFat"]][["return_months"]]
  CaseFat_postMonth <- CaseFat_RM[length(CaseFat_RM)] + 1
  CaseFat_mult <- return_params[["CaseFat"]][["multiplier"]]
  CaseFat_lastMonth <- CaseFat_RM[1] - 1

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
