#' THIS FUNCTION INPUTS A TABLE OF PARAMETERS AND RUNS THE TB MODEL
#' AND GENERATES AN ARRAY OF OUTPUTS FOR EACH POSSIBLE PREDEFINED
#' SCENARIO; THIS ARRAY IS THEN SAVED TO THE DATA FOLDER FOR THAT
#' LOCATION

#'@name make_all_scenarios
#'@param loc two character location code for the model location
#'@param ParMatrix parameters to use in the simulation
#'@export
make_all_scenarios<-function(loc, ParMatrix){
  sim_list<-list()
  for (intv in 0:8) {
    # intvs will be a vector of 0s except for (possibly) one activated intervention
    intvs <- rep(0, 8)
    # If intv is 0, disable all interventions to run the basecase.
    # Otherwise, activate one of the 8 interventions.
    if (intv != 0) intvs[intv] <- 1
    #use the basecase program change
    # defprg<-def_prgchng(ParMatrix[1,])
    # Simulate using the elements of intvs to control
    # whether or not each intervention is on.
    out<-new_OutputsInt(loc=loc,ParMatrix=Par,n_cores = 1, endyr=2050,
                        Int1 = intvs[[1]],
                        Int2 = intvs[[2]],
                        Int3 = intvs[[3]],
                        Int4 = intvs[[4]],
                        Int5 = intvs[[5]],
                        Scen1 = intvs[[6]],
                        Scen2 = intvs[[7]],
                        Scen3 = intvs[[8]],
                        prg_chng =def_prgchng(Par[1,]))
    # print(intvs)
    # Do something to save the simulation outcomes
    # save(out,file=paste0("/Users/nis100/MITUS/inst/",loc,"/",loc,"_results_",intv+1,".rda"))
    #create a list of the simulations
    sim_list[[intv+1]]<-out
  }
  return(sim_list)
}
