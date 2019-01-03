#' the programatic changes provided by tabbyII are:
#' LTBI treatment cascade:
#' 1.Screening coverage for currently defined risk groups (Multiple of current coverage, from 1.0 to 5.0 times)
#' 2.Fraction receiving IGRA (% from 0 to 100) (rest assumed to get TST)
#' 3.Fraction of individuals testing positive who accept treatment (0% to 100%)
#' 4.Fraction of individuals initiating treatment who complete treatment (0% to 100%)
#' 5.Treatment efficacy for individuals completing treatment (0% to 100%)
#'
#' TB Treatment Cascade
#' 1.Average time to treatment for incident TB cases (0-100% of current value)
#' 2.Fraction discontinuing/defaulting from treatment (0% to 100%)
#' @name prg_chg
#' @param P vector of baseline Parameter Values
#' @param UI_prgm_chng vector of UI values for program changes
#' @return P vector of updated baseline Parameter Values
prm_chg <- function(P,UI_prgm_chng) {
if (prgm_chng==TRUE){
  #load the User inputs
  #rate of screening coverage as a multiple of the current rate
  #baseline rate is 0.025 so we must limit how high this can go (output to user)
  if(UI_prgm_chng$scrn_cov!=1){
    P["rLtScr"]<-P["rLtScr"]*UI_prgm_chng$scrn_cov
  }
  if(UI_prgm_chng$IGRA_frc!=b_IGRA_frc){
    #changes the sensitivity and specificity params
    #SensLt and SpecLt
    #also add these to P and allow them to change
    #however the default values are based on IGRA
    #need to update these fields but for now
    P["b_IGRA_frc"]<-UI_prgm_chng$IGRA_frc
  }
  if(UI_prgm_chng$ltbi_init_frc){
    #this parameter is not something fed into the model
    #needs to be added to P and held fixed so I can manipulate it easier
    #baseline level is 80% for national level (will differ based off of geography?)
    P["pTlInt"] <- UI_prgm_chng$ltbi_init_frc
  }
  if(UI_prgm_chng$ltbi_comp_frc){
    P["pDefLt"]<-1-(UI_prgm_chng$ltbi_comp_frc)
  }
  if(UI_prgm_chng$ltbi_eff_frc){
    P["EffLt"]<-1-(UI_prgm_chng$ltbi_eff_frc)
  }
  if(UI_prgm_chng$tb_tim2tx_frc){
    P["DelaySp"]<-P["DelaySp"]*UI_prgm_chng$tb_tim2tx_frc
  }
  if(UI_prgm_chng$tb_txdef_frc){
    rDeft[yrs,]<-1-(UI_prgm_chng$ltbi_txdef_frc)
  }
} else{P<-P}
  return(P)
} #end of function
