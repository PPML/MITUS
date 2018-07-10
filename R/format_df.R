#' Formats Outputs to dataframe format for the TABBY2 interface

#'@name format_df
#'@param years vector of years
#'@param loc location of results
#'@param nat nationality of results
#'@param scen scenario
#'@param int  intervention
#'@param df dataframe of results
#'@return csv of results formatted
#'@export


format_df<-function(year, loc, nat, scen, int, df){

Date <- Sys.Date();

#' pull the year from df data frame
Year <- rep(NA,nrow(df)*(ncol(df)-1))
Year <-rep(df[,1], (ncol(df)-1))

#'location will need to be inputted from user
#'for now we will input the national example
Location <- "Nat"

#'Nativity Selection will need to be input by the User
#'for now, will input all;
#update this so that it is just the nativity group of the outputs
Nativity <- "all"

#' fill in the scenario number
Scenario <-0

if (Scen1==1)
 Scenario <- 1
if (Scen2==1)
  Scenario <- 2
if (Scen3==1)
  Scenario <- 3

#' fill in the Intervention Number
Intervention <-0

if (Int1==1)
Intervention <- 1
if (Int2==1)
Intervention <- 2
if (Int3==1)
Intervention <- 3
if (Int4==1)
Intervention <- 4
if (Int5==1)
Intervention <- 5



#' Input the Output Name from InputParams[["ResNam"]]
Output <- rep(NA,(length(InputParams[["ResNam"]])-1)*nrow(df))
  for (i in 1:(length(InputParams[["ResNam"]])-1)){
      Output[(((i-1)*100)+1):(i*100)] <-InputParams[["ResNam"]][i+1]
 }

#' Input the values from df dataframe
Value <- rep(NA,(nrow(df)*(ncol(df)-1)))

for (i in 1:nrow(df)){
  for (j in 1:(ncol(df)-1)){

    Value[i+(j-1)*100] <- df[i,j+1]

} }

results_form <- data.frame(Date,Location,Year,Nativity,Scenario,Intervention,Output,Value) #8

save(results_form,file = paste("data/results_form", Sys.time(),".rData"))
write.csv(results_form, file = paste("MITUS_results/results_form", Sys.time(),".csv"))

}
