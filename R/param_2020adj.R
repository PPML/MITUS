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
                           multiplier = 1 # Default is from DHS OIS LPR Annual flow
                                             # https://www.dhs.gov/sites/default/files/2023-02/2022_0405_plcy_lawful_permanent_residents_fy2021v2.pdf
){
  if (immig != 99){
  # newImmig <- 1 - immig
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
  immigrationTrend <- immigrationData/85

  ### model trend
  # modelTrend0 <- immig*immigrationData
  # modelTrend <- modelTrend0/modelTrend0[1]; plot(modelTrend)
  # TotImmAge[effect_months,] <- TotImmAge[effect_months,] - (TotImmAge[effect_months, ] * immig*immigrationTrend);


  # Return to pre-pandemic rate (1 * 843) by end of 2023
  # Month of final return rate is 888 (december 2023)
  for (agegrp in 1:ncol(TotImmAge)){
    # TotImmAge[effect_months,agegrp] <- TotImmAge[843,agegrp] * immig*immigrationTrend;
    TotImmAge[effect_months,agegrp] <- TotImmAge[843,agegrp] + (TotImmAge[843,agegrp] * (immigrationTrend-1)*immig);


    # TotImmAge[846:lastMonth,agegrp] <- seq(TotImmAge[846,agegrp],(TotImmAge[843,agegrp] * 1.3), length.out=length(846:lastMonth));

    # TotImmAge[846:last(effect_months), agegrp] <- logseq(TotImmAge[846,agegrp], (TotImmAge[843,agegrp] * 1.3), n=length(846:last(effect_months)));

    TotImmAge[return_months,agegrp] <- seq(TotImmAge[lastMonth,agegrp],TotImmAge[postMonth,agegrp]*multiplier, length.out=length(return_months))
    # TotImmAge[return_months,agegrp] <- logseq(TotImmAge[lastMonth,agegrp], TotImmAge[postMonth-1,agegrp] * multiplier, n=length(return_months))
    # Continue this trend from here to end of model run
    TotImmAge[postMonth:nrow(TotImmAge),agegrp] <- TotImmAge[postMonth:nrow(TotImmAge),agegrp] * multiplier
  }
  }
  # print("Immig: GOOD")
  return(TotImmAge)
}


# plot(rowSums(TotImmAge[840:890,]), type = "l")
# abline(v=which(840:890 == 862))

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

  if (par2020["DxtKnot1"] != 0){
  # Setup params
  # rDxt_RM <- return_params[["rDxt"]][["return_months"]]
  # rDxt_postMonth <- rDxt_RM[length(rDxt_RM)] + 1
  # rDxt_lastMonth <- rDxt_RM[1] - 1
  rDxt_mult <- return_params[["rDxt"]][["multiplier"]]

  # rDxt_trend <- readRDS(file = "~/Documents/COVIDTB Paper/Data/careSeekingTrend.rds")
  # rDxt_trendNorm <- rDxt_trend / rDxt_trend[1]
  rDxt_lastMonth <- 843+23
  rDxt_RM <- (rDxt_lastMonth+1):(rDxt_lastMonth + 24)
  rDxt_postMonth <- last(rDxt_RM) + 1
  # rDxt_mult <- 1.1
  ### Calculate the rate of diagnosis

  # rDxt_trend <- c(seq(par2020["DxtKnot1"],par2020["DxtKnot2"],length.out = 6),
  #                  seq(par2020["DxtKnot2"],par2020["DxtKnot3"],length.out = 6),
  #                  seq(par2020["DxtKnot3"],par2020["DxtKnot4"],length.out = 6),
  #                  seq(par2020["DxtKnot4"],par2020["DxtKnot5"],length.out = 6))

  rDxtDF<-data.frame(Month = c(0,3,9,15,21,27),
                     Value = c(0,par2020["DxtKnot1"], par2020["DxtKnot2"], par2020["DxtKnot2"],
                               par2020["DxtKnot4"], par2020["DxtKnot5"]))


  rDxtfit = lm(Value~bs(Month, knots = c(0,3,9,15,21,27), degree=1), data = rDxtDF)

  rDxt_trend <- predict(rDxtfit, newdata=list(Month=1:24))
  # Initial covid effects
  rDxt[843:rDxt_lastMonth,] <- rDxt[843:rDxt_lastMonth,] - (rDxt[843:rDxt_lastMonth,] * rDxt_trend)

  # Bring up params to X% by end of return months (smoothly)
  for (riskgrp in 1:ncol(rDxt)){
    rDxt[rDxt_RM,riskgrp] <- seq(rDxt[rDxt_lastMonth,riskgrp], rDxt[842,riskgrp]*rDxt_mult, length.out=length(rDxt_RM))

    # rDxt[rDxt_RM,riskgrp] <- logseq(rDxt[843,riskgrp], rDxt[843,riskgrp]*rDxt_mult, n=length(rDxt_RM))
    rDxt[rDxt_postMonth:nrow(rDxt), riskgrp] <-  rDxt[rDxt_postMonth:nrow(rDxt), riskgrp]*rDxt_mult
  }
  } else {
    rDxt <- rDxt
  }
  # plot(rDxt[840:890,1])
  # print("Dxt: GOOD")

  #############################################################################
  ##### NixTrans #####
  #############################################################################
  if (par2020["TransKnot1"] != 0){
  # Setup params
  # NixTrans_RM <- return_params[["Trans"]][["return_months"]]
  # NixTrans_postMonth <- NixTrans_RM[length(NixTrans_RM)] + 1
  NixTrans_mult <- return_params[["Trans"]][["multiplier"]]
  # NixTrans_lastMonth <- NixTrans_RM[1] - 1

  # NixTrans_trend <- readRDS(file = "~/Documents/COVIDTB Paper/Data/contactRateTrend.rds")
  # NixTrans_trendNorm <- NixTrans_trend/NixTrans_trend[1]
  # NixTrans_lastMonth <- 843+length(NixTrans_trend)-1
  # NixTrans_postMonth <- NixTrans_lastMonth + 1

  NixTrans_lastMonth <- 843+23
  NixTrans_RM <- (NixTrans_lastMonth+1):(NixTrans_lastMonth + 24)
  NixTrans_postMonth <- last(NixTrans_RM) + 1


  ### Calculate the trend
  # NixTrans_trend <- c(seq(par2020["TransKnot1"],par2020["TransKnot2"],length.out = 6),
  #                     seq(par2020["TransKnot2"],par2020["TransKnot3"],length.out = 6),
  #                     seq(par2020["TransKnot3"],par2020["TransKnot4"],length.out = 6),
  #                     seq(par2020["TransKnot4"],par2020["TransKnot5"],length.out = 6))

  NixTransDF<-data.frame(Month = c(0,3,9,15,21,27),
                         Value = c(0,par2020["TransKnot1"], par2020["TransKnot2"], par2020["TransKnot2"],
                                   par2020["TransKnot4"], par2020["TransKnot5"]))


  NixTrans_fit = lm(Value~bs(Month, knots = c(0,3,9,15,21,27), degree=1), data = NixTransDF)

  NixTrans_trend <- predict(NixTrans_fit, newdata=list(Month=1:24))
  # Initial covid effects
  NixTrans[843:NixTrans_lastMonth] <- NixTrans[843:NixTrans_lastMonth] - (NixTrans[843:NixTrans_lastMonth] * NixTrans_trend)
  NixTrans[NixTrans_RM] <- seq(NixTrans[NixTrans_lastMonth], NixTrans[842]*NixTrans_mult, length.out=length(NixTrans_RM))



  # NixTrans[843:NixTrans_lastMonth]<- 1-par2020["Trans"]*NixTrans_trendNorm
  # NixTrans[NixTrans_RM] <- seq(NixTrans[NixTrans_lastMonth], NixTrans[843]*NixTrans_mult, length.out=length(NixTrans_RM))
  # NixTrans[NixTrans_RM] <- logseq(NixTrans[NixTrans_lastMonth], NixTrans[843]*NixTrans_mult, n=length(NixTrans_RM))
  NixTrans[NixTrans_postMonth:length(NixTrans)] <- NixTrans[NixTrans_postMonth:length(NixTrans)] * NixTrans_mult
  # plot(NixTrans[840:890])
  # print("Trans: GOOD")
  } else {
  NixTrans<-NixTrans
}
  #############################################################################
  ##### Changes to mortality with TB #####
  #############################################################################
  if (par2020["CaseFatKnot1"] != 1){

  # Setup params
  # CaseFat_RM <- return_params[["CaseFat"]][["return_months"]]
  # CaseFat_postMonth <- CaseFat_RM[length(CaseFat_RM)] + 1
  CaseFat_mult <- return_params[["CaseFat"]][["multiplier"]]
  # CaseFat_lastMonth <- CaseFat_RM[1] - 1

  CaseFat_lastMonth <- 843+23
  CaseFat_RM <- (CaseFat_lastMonth+1):(CaseFat_lastMonth + 24)
  CaseFat_postMonth <- last(CaseFat_RM) + 1

  ## Calculate the trend
  # CaseFat_trend <- c(seq(par2020["CaseFatKnot1"],par2020["CaseFatKnot2"],length.out = 6),
  #                    seq(par2020["CaseFatKnot2"],par2020["CaseFatKnot3"],length.out = 6),
  #                    seq(par2020["CaseFatKnot3"],par2020["CaseFatKnot4"],length.out = 6),
  #                    seq(par2020["CaseFatKnot4"],par2020["CaseFatKnot5"],length.out = 6))

  CaseFatDF<-data.frame(Month = c(0,3,9,15,21,27),
                         Value = c(0,par2020["CaseFatKnot1"], par2020["CaseFatKnot2"],
                                   par2020["CaseFatKnot2"],
                                   par2020["CaseFatKnot4"], par2020["CaseFatKnot5"]))


  CaseFat_fit = lm(Value~bs(Month, knots = c(0,3,9,15,21,27), degree=1), data = CaseFatDF)

  CaseFat_trend <- predict(CaseFat_fit, newdata=list(Month=1:24))

  RRmuTBPand <- rep(1,1812)
  time_horizon <- 843:888
  RRmuTBPand[time_horizon[1]:1812] <-c(CaseFat_trend,
                                       seq(last(CaseFat_trend), 1*CaseFat_mult, length.out = length(CaseFat_RM)),
                                       # logseq(par2020["CaseFat"], 1*CaseFat_mult, n = length(CaseFat_RM)),
                                       rep(1*CaseFat_mult, length.out = length(CaseFat_postMonth:1812)))
  }
  else {RRmuTBPand <- rep(1,1812) }
  # print("CaseFat: GOOD")
  # plot(RRmuTBPand[840:890])

  param_list <- list("rDxt" = rDxt,
                     "NixTrans" = NixTrans,
                     "RRmuTBPand" = RRmuTBPand)

  return(param_list)
}

# plot(y=RRmuTBPand[843:888], x=843:888, type = "l")
# abline(v=876);abline(v=877);
# abline(v=which(843:888 == 862))
