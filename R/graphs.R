#' creates plots of various outputs of the model
#' This script plots a subset of the outputs of the model
#' and writes these results to a .pdf file.
#' This script is primarily useful for evaluating the model
#' output for logic and for identifying potential bugs.
#'

#'Create a function to be run on a specific model run output to
#'create easy to run graphs of the output
#'@param start_yr year to start the graphs
#'@param end_yr year to end the graphs
#'@param df dataframe of output for all years
#'@return .pdf of the graphs
#'@export


tb_graphs <- function(start_yr, end_yr, df){

#'Check that the column names of the selected dataframe are equal to those
#'defined in ResNam. The column names are used to identify the specific
#'results plotted below. The code will produce an error if they are not equal.
  if (all.equal(colnames(df),ResNam)==FALSE)
    stop("Column Names of 'df' do not match those defined in ResNam")

  #pdfname <- paste("MITUS_results/graphs",Sys.time(),".pdf")
  pdf(file=paste("MITUS_results/graphs",Sys.time(),".pdf"), width = 11, height = 8.5)

  plot(0,0,ylim=c(0,85),xlim=c(start_yr,end_yr),xlab="Year",ylab="",axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

#'Total diagnosed cases 1953-2013
 notif_tot=df$NOTIF_ALL+df$NOTIF_MORT_ALL
 lines(start_yr:end_yr, (notif_tot[(start_yr-1950):(end_yr-1950)])*10)


 dev.off();
 #system(paste("open", pdfname))
}


