#'@name  adj_immig_2020
#'@param TotImmAge
#'@param immig
#'@param return_months
#'@param multiplier
#'@return TotImmAge
#'@export
adj_immig_2020 <- function(TotImmAge,
                           immig = 99,
                           effect_months = 843:861,
                           return_months =  862:888,
                           returnTrendTunParam = 1,
                           multiplier = 1 # Default is from DHS OIS LPR Annual flow
                                             # https://www.dhs.gov/sites/default/files/2023-02/2022_0405_plcy_lawful_permanent_residents_fy2021v2.pdf
){
  if (immig != 99){
  #Calculate first month after return_months
  postMonth <- return_months[length(return_months)] + 1
  # Calculate the last month of the initial pandemic effect
  lastMonth <- return_months[1] - 1
  # Calculate the during pandemic effect
  # This period represents the period for which we have data
  # Data shows that immigration decreased through May to ~ 1/3rd of level in 2019
  # then rebounded to 1.33 time 2019 level by September 2021
  # TotImmAge[843:846,] <- TotImmAge[843:846,] - (TotImmAge[843:846, ] * min(max(immig,1e-12), .99));
  immigrationData <- c(74,27,21,28,38,48,47,52,39, 40,40,40,60,59, 59, 66, 80, 96, 111)
  immigrationTrend <- immigrationData/immigrationData[1]
  TotImmAge[effect_months,] <- TotImmAge[effect_months,] - (TotImmAge[effect_months, ] * immig*immigrationTrend);

  # Return to pre-pandemic rate (1 * 843) by end of 2023
  # Month of final return rate is 888 (december 2023)
  for (agegrp in 1:ncol(TotImmAge)){
    # TotImmAge[846:lastMonth,agegrp] <- seq(TotImmAge[846,agegrp],(TotImmAge[843,agegrp] * 1.3), length.out=length(846:lastMonth));

    # TotImmAge[846:last(effect_months), agegrp] <- logseq(TotImmAge[846,agegrp], (TotImmAge[843,agegrp] * 1.3), n=length(846:last(effect_months)));

    # TotImmAge[return_months,agegrp] <- seq(TotImmAge[lastMonth,agegrp],TotImmAge[843,agegrp]*multiplier, length.out=length(return_months))
    TotImmAge[return_months,agegrp] <- logseq(TotImmAge[lastMonth,agegrp], TotImmAge[postMonth-1,agegrp] * multiplier, n=length(return_months))
    # Continue this trend from here to end of model run
    TotImmAge[postMonth:nrow(TotImmAge),agegrp] <- TotImmAge[postMonth:nrow(TotImmAge),agegrp] * multiplier
  }
  }
  # print("Immig: GOOD")
  return(TotImmAge)
}

# plot(1, xlim = c(0 , 1), ylim = c(0, 21))
# plot((TotImmAge[return_months,3])^exp(0)*(immig)^(1-exp(0)), type="l")
# lines((TotImmAge[return_months,3])^exp(.1)*(TotImmAge[postMonth,3])^(1-exp(.1)), col="red")
# lines((TotImmAge[return_months,3])^exp(.25)*(TotImmAge[postMonth,3])^(1-exp(.25)), col="orange")
# lines((TotImmAge[return_months,3])^exp(.5)*(TotImmAge[postMonth,3])^(1-exp(.5)), col="yellow")
# lines((TotImmAge[return_months,3])^exp(.75)*(TotImmAge[postMonth,3])^(1-exp(.75)), col="green")
# lines((TotImmAge[return_months,3])^exp(.9)*(TotImmAge[postMonth,3])^(1-exp(.9)), col="blue")
# lines((TotImmAge[return_months,3])^exp(1)*(TotImmAge[postMonth,3])^(1-exp(1)), col="purple")

# plot((TotImmAge[return_months,1])^exp(.1)*(TotImmAge[postMonth,1])^(1-exp(.1)), col="red")
# points((TotImmAge[return_months,1])^exp(.25)*(TotImmAge[postMonth,1])^(1-exp(.25)), col="orange")
# points((TotImmAge[return_months,1])^exp(.5)*(TotImmAge[postMonth,1])^(1-exp(.5)), col="yellow")
# points((TotImmAge[return_months,1])^exp(.75)*(TotImmAge[postMonth,1])^(1-exp(.75)), col="green")
# points((TotImmAge[return_months,1])^exp(.9)*(TotImmAge[postMonth,1])^(1-exp(.9)), col="blue")
# points((TotImmAge[return_months,1])^exp(1)*(TotImmAge[postMonth,1])^(1-exp(1)), col="purple")


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

  ###
  #############################################################################
  ##### rDxt #####
  #############################################################################
  ### We would like this trend to be piecewise linear splines

  # Setup params
  rDxt_RM <- return_params[["rDxt"]][["return_months"]]
  # rDxt_postMonth <- rDxt_RM[length(rDxt_RM)] + 1
  # rDxt_lastMonth <- rDxt_RM[1] - 1
  rDxt_mult <- return_params[["rDxt"]][["multiplier"]]

  rDxt_trend <- readRDS(file = "~/Documents/COVIDTB Paper/Data/careSeekingTrend.rds")
  rDxt_trendNorm <- rDxt_trend / rDxt_trend[1]
  rDxt_lastMonth <- 843+length(rDxt_trend) - 1
  rDxt_postMonth <- rDxt_lastMonth + 1


  # Initial covid effects
  rDxt[843:rDxt_lastMonth,] <- rDxt[843:rDxt_lastMonth,] - (rDxt[843:rDxt_lastMonth,] * par2020["Dxt"] * rDxt_trendNorm)

  # Bring up params to X% by end of return months (smoothly)
  for (riskgrp in 1:ncol(rDxt)){
    # rDxt[rDxt_RM,riskgrp] <- seq(rDxt[843,riskgrp], rDxt[843,riskgrp]*rDxt_mult, length.out=length(rDxt_RM))

    # rDxt[rDxt_RM,riskgrp] <- logseq(rDxt[843,riskgrp], rDxt[843,riskgrp]*rDxt_mult, n=length(rDxt_RM))
    rDxt[rDxt_postMonth:nrow(rDxt), riskgrp] <-  rDxt[rDxt_postMonth:nrow(rDxt), riskgrp]*rDxt_mult
  }

  # plot(rDxt[,1])
  # print("Dxt: GOOD")

  #############################################################################
  ##### NixTrans #####
  #############################################################################

  # Setup params
  # NixTrans_RM <- return_params[["Trans"]][["return_months"]]
  # NixTrans_postMonth <- NixTrans_RM[length(NixTrans_RM)] + 1
  NixTrans_mult <- return_params[["Trans"]][["multiplier"]]
  # NixTrans_lastMonth <- NixTrans_RM[1] - 1

  NixTrans_trend <- readRDS(file = "~/Documents/COVIDTB Paper/Data/contactRateTrend.rds")
  NixTrans_trendNorm <- NixTrans_trend/NixTrans_trend[1]
  NixTrans_lastMonth <- 843+length(NixTrans_trend)-1
  NixTrans_postMonth <- NixTrans_lastMonth + 1

  NixTrans[843:NixTrans_lastMonth]<- 1-par2020["Trans"]*NixTrans_trendNorm
  # NixTrans[NixTrans_RM] <- seq(NixTrans[NixTrans_lastMonth], NixTrans[843]*NixTrans_mult, length.out=length(NixTrans_RM))
  # NixTrans[NixTrans_RM] <- logseq(NixTrans[NixTrans_lastMonth], NixTrans[843]*NixTrans_mult, n=length(NixTrans_RM))
  NixTrans[NixTrans_postMonth:length(NixTrans)] <- NixTrans[NixTrans_postMonth:length(NixTrans)] * NixTrans_mult
  # plot(NixTrans)
  # print("Trans: GOOD")

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
                                       # seq(par2020["CaseFat"], 1*CaseFat_mult, length.out = length(CaseFat_RM)),

                                       logseq(par2020["CaseFat"], 1*CaseFat_mult, n = length(CaseFat_RM)),
                                       rep(1*CaseFat_mult, length.out = length(CaseFat_postMonth:1812)))

  # print("CaseFat: GOOD")
  # plot(RRmuTBPand)

  param_list <- list("rDxt" = rDxt,
                     "NixTrans" = NixTrans,
                     "RRmuTBPand" = RRmuTBPand)

  return(param_list)
}
