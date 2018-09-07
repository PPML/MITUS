#'@name run_MITUS
#'@param loc character string length 2 that determines the geography of the simulation
#'@return dataframe of results

#enter location:
run_MITUS <-function(loc="US"){
model_load(loc)
OutputsZint(1,P,startyr=1950,endyr=2050,Int1=0,Int2=0,Int3=0,Int4=0,Int5=0,Scen1=0,Scen2=0,Scen3=0)
tabby_results(results)
return(t_results)
}
