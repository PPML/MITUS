#' Formats Outputs to dataframe format for the TABBY2 interface

Date <- Sys.Date();

#' pull the year from results data frame
Year <- rep(NA,nrow(results)*(ncol(results)-1))
Year <-rep(results[,1], (ncol(results)-1))

#'location will need to be inputted from user
#'for now we will input the national example
Location <- "Nat"

#'Nativity Selection will need to be input by the User
#'for now, will input all;
#'
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



#' Input the Output Name from ResNam
Output <- rep(NA,(length(ResNam)-1)*nrow(results))
  for (i in 1:(length(ResNam)-1)){
      Output[(((i-1)*100)+1):(i*100)] <-ResNam[i+1]
 }

#' Input the values from results dataframe
Value <- rep(NA,(nrow(results)*(ncol(results)-1)))

for (i in 1:nrow(results)){
  for (j in 1:(ncol(results)-1)){

    Value[i+(j-1)*100] <- results[i,j+1]

} }

results_form <- data.frame(Date,Location,Year,Nativity,Scenario,Intervention,Output,Value) #8

save(results_form,file = paste("data/results_form", Sys.time(),".rData"))
write.csv(results_form, file = paste("MITUS_results/results_form", Sys.time(),".csv"))

