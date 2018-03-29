#' creates plots of various outputs of the model
#' This script plots a subset of the outputs of the model
#' and writes these results to a .pdf file.
#' This script is primarily useful for evaluating the model
#' output for logic and for identifying potential bugs.
#'

#'Create a function to be run on a specific model run output to
#'create simple graphs of all the output for a selected year range
#'@param start_yr year to start the graphs
#'@param end_yr year to end the graphs
#'@param df dataframe of output for all years
#'@return .pdf of the graphs
#'@export
tb_graph_all <- function(start_yr, end_yr, df){

  pdf(file=paste("MITUS_results/graphs_all",Sys.time(),".pdf"), width = 11, height = 8.5)
  i<-1
  while (i <= ncol(df)){

  plot(0,0,ylim=range(df[,i]),xlim=c(start_yr,end_yr),xlab="Year",ylab=colnames(df)[i],axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  lines(start_yr:end_yr,results[(start_yr-1949):(end_yr-1949),i])

  print(i)
  i<- i+1
}

  dev.off();
  #system(paste("open", pdfname))
}


#'Create a function to be run on a specific model run output to
#'create easy to run graphs of the output
#'@param start_yr year to start the graphs
#'@param end_yr year to end the graphs
#'@param df dataframe of output for all years
#'@param output character string of output desired
#'@return .pdf of the graphs
#'@export


tb_graph_specific <- function(start_yr, end_yr, df, output){

  #'Check that the column names of the selected dataframe are equal to those
  #'defined in ResNam. The column names are used to identify the specific
  #'results plotted below. The code will produce an error if they are not equal.
  if (all.equal(colnames(df),ResNam)==FALSE)
    stop("Column Names of 'df' do not match those defined in ResNam")

  #'perform character pattern matching to create a vector of the results to graph
  if(length(i <- grep(output, colnames(df),ignore.case = TRUE)))

     plot_these <- df[,i]

  pdf(file=paste("MITUS_results/graphs_",output,Sys.time(),".pdf"), width = 11, height = 8.5)

  j<-1
  while(j <= ncol(plot_these)){

  plot(0,0,ylim=range(plot_these[,j]),xlim=c(start_yr,end_yr),xlab="Year",ylab=colnames(plot_these)[j],axes=F)
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(start_yr:end_yr, plot_these[(start_yr-1949):(end_yr-1949),j])
  print(j)
  j<-j+1
}
  dev.off()
}



#' start_yr=1950
#' end_yr=2050
#' #'Total diagnosed cases 1953-2013
#' notif_tot=results$NOTIF_ALL+results$NOTIF_MORT_ALL
#' lines(start_yr:end_yr, (notif_tot[(start_yr-1949):(end_yr-1949)])*1e1)
#'
#' dev.off();

