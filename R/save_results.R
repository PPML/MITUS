#'function to save results
#'@title save.results
#'@param df dataframe of output
#'@return .csv of results (unformatted)
#'@export

save.results <- function(df)
{
save(df,file = paste("data/results", Sys.time(),".rData"))
write.csv(df, file = paste("MITUS_results/results", Sys.time(),".csv"))
}
