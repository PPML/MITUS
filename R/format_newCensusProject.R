#'This function reads in the 2023 census projections and appends
#'the revised non-U.S.-born projections to the ImmigInputs list
#'contained in the Inputs list, which is loaded in at model start.
#'When running a simulation, the user will need to designate
#'whether they want to use these projections.
#'


updatePopProjections <- function(){
  ### Read in the revised census projections
  ### Note: race and sex 0 stratification represent totals
  revisedProjections <- as.data.frame(read.csv("~/MITUS/inst/extdata/US/netInternalMigrationStratifiedProjections.csv")) %>% filter(SEX == 0, RACE_HISP == 0)

  ### Check with a plot
  plot(x = revisedProjections$YEAR, y = revisedProjections$TOTAL_NIM / 1e6,
       ylab = "Total net migration", xlab = "Year",
       main = "Revised Census projections of net international migration (Total in millions)",
       type = "l")

  ### Calculate the annual rate of change in international migration
  revisedProjections$ANNUAL_PC <- ((revisedProjections$TOTAL_NIM - lag(revisedProjections$TOTAL_NIM)) / revisedProjections$TOTAL_NIM)

  ### Check with a plot
  plot(x = revisedProjections$YEAR, y = revisedProjections$ANNUAL_PC*100,
       ylab = "Annual percentage change in migration", xlab = "Year",
       main = "Revised Census projections of net international migration (Annual percentage change)",
       type = "l")

  ### Create an object to hold the new immigration estimates
  newImmigPop <- oldImmigPop <- Inputs[["ImmigInputs"]]$TotByYear

  newImmigPop[(which(names(newImmigPop) == 2024)):
              (which(names(newImmigPop) == 2100))] <- newImmigPop[(which(names(newImmigPop) == 2023)):
                                                                  (which(names(newImmigPop) == 2099))] * (1 + revisedProjections$ANNUAL_PC[-1])

  ### Check with plot
  plot(x=2023:2100, y=newImmigPop[(which(names(newImmigPop) == 2023)):
                                  (which(names(newImmigPop) == 2100))],
       xlab = "Year", ylab = "Net immigration (in millions)",
       main = "Updated net immigration (in millions)",
       type = "l")
  lines(x=2023:2100, y=oldImmigPop[(which(names(oldImmigPop) == 2023)):
                                   (which(names(oldImmigPop) == 2100))], col = "red")
}

# updatePopProjections <- function(){
#   revisedProjections <- as.data.frame(read.csv("~/MITUS/inst/extdata/US/netInter.csv"))
#
#   nusbProjectedPop <- (revisedProjections$NUSB.Population.Difference / 1e6) [-1]
#   newImmigPop <- oldImmigPop <- Inputs[["ImmigInputs"]]$TotByYear
#   range(oldImmigPop)
#
#
#   newImmigPop[(which(names(newImmigPop) == 2023)):
#               (which(names(newImmigPop) == 2100))] <- nusbProjectedPop
#
#   ### Check with plot
#
#   plot(x=names(newImmigPop), y=newImmigPop,
#        xlab = "Year", ylab = "New NUSB population (in millions)",
#        main = "New NUSB population over time")
#   abline(v=2023)
# }
