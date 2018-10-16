#' creates plots of various outputs of the model
#' This script plots a subset of the outputs of the model
#' and writes these results to a .pdf file.
#' This script is primarily useful for evaluating the model
#' output for logic and for identifying potential bugs.
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
#'@name tb_graph_specific
#'@param start_yr year to start the graphs
#'@param end_yr year to end the graphs
#'@param df dataframe of output for all years
#'@param output character string of output desired
#'@return .pdf of the graphs
#'@export
tb_graph_specific <- function(start_yr, end_yr, df, output){
  if (all.equal(colnames(df),ResNam)==FALSE)
    #Check that the column names of the selected dataframe are equal to those
    #defined in ResNam. The column names are used to identify the specific
    #results plotted below. The code will produce an error if they are not equal.
    stop("Column Names of 'df' do not match those defined in ResNam")

  #perform character pattern matching to create a vector of the results to graph
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

#'Create a function to be run that will graph the most important outputs
#'often sums for a designated year range.
#'@name tb_graph_vital
#'@param start_yr year to start the graphs
#'@param end_yr year to end the graphs
#'@param df dataframe of output for all years
#'@return .pdf of the graphs
#'@export
tb_graph_vital <- function(start_yr, end_yr, df){
  if (all.equal(colnames(df),ResNam)==FALSE)  #Create an empty list to hold the formatted intitial parameters
    stop("Column Names of 'df' do not match those defined in ResNam")
# open pdf file for the graphs
pdf(file=paste("MITUS_results/graphs_vital",Sys.time(),".pdf"), width = 11, height = 8.5)
# total diagnosed cases
plot_this <- results[(start_yr-1949):(end_yr-1949),"NOTIF_MORT_ALL"]+results[(start_yr-1949):(end_yr-1949),"NOTIF_ALL"]
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="All TB Cases", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)
# cases by nativity
plot_this <- results[(start_yr-1949):(end_yr-1949),"INCID_ALL_US"]
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="US Incident TB Cases", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)

plot_this<-results[(start_yr-1949):(end_yr-1949),"INCID_ALL_FB"] + results[(start_yr-1949):(end_yr-1949),"INCID_ALL_FB2"]
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="All FB Incident TB Cases", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)

plot_this<- results[(start_yr-1949):(end_yr-1949),"INCID_ALL_FB"]
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="All Recent FB Incident TB Cases", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)

# incident cases by age distribution in percents
plot_this <- cbind(results[(start_yr-1949):(end_yr-1949),"INCID_ALL_0_4"]   ,results[(start_yr-1949):(end_yr-1949),"INCID_ALL_5_14"] ,
                   results[(start_yr-1949):(end_yr-1949),"INCID_ALL_15_24"] ,results[(start_yr-1949):(end_yr-1949),"INCID_ALL_25_34"] ,
                   results[(start_yr-1949):(end_yr-1949),"INCID_ALL_35_44"] ,results[(start_yr-1949):(end_yr-1949),"INCID_ALL_45_54"] ,
                   results[(start_yr-1949):(end_yr-1949),"INCID_ALL_55_64"] ,results[(start_yr-1949):(end_yr-1949),"INCID_ALL_65_74"] ,
                   results[(start_yr-1949):(end_yr-1949),"INCID_ALL_75_84"] ,results[(start_yr-1949):(end_yr-1949),"INCID_ALL_85_94"] ,
                   results[(start_yr-1949):(end_yr-1949),"INCID_ALL_95p"]   )
plot_this <- colSums(plot_this)/sum(plot_this)*100
plot(0,0,ylim=range(plot_this),xlim=c(0.6,11.4),xlab="Age Group",ylab="Total TB Cases by Age Group", axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94", "95+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
lines(plot_this)

# cases by HR distribution
plot_this <- (results[(start_yr-1949):(end_yr-1949), "INCID_ALL_HR"]) / (results[(start_yr-1949):(end_yr-1949), "INCID_ALL"])
plot_this <- plot_this*100
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="Incident High Risk TB Cases (%)", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
plot(start_yr:end_yr,plot_this)

# treatment outcomes
plot_this <- results[(start_yr-1949):(end_yr-1949),"TBTX_COMPLT"]
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="TX Completion", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)

plot_this <- cbind(results[(start_yr-1949):(end_yr-1949),"TBTX_DIED"] +results[(start_yr-1949):(end_yr-1949),"TBTX_DISCONT"])
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="Treatment Outcomes: Discontinued and Died", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)
#TLTBI volume (very, very low!!)
plot_this <- results[(start_yr-1949):(end_yr-1949),"TLTBI_INITS"]
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="IPT Treatment Initiations Per Year (000s)", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(start_yr:end_yr,plot_this)

# LTBI initiations by risk group
# end_yr=1990
if (end_yr > 1985){
plot_this <- (results[(1985-1949):(end_yr-1949),"TLTBI_INITS_HR"]) / (results[(1985-1949):(end_yr-1949),"TLTBI_INITS"])
plot(0,0,ylim=range(plot_this),xlim=c(1985,end_yr),xlab="Year",ylab="IPT Treatment Initiations, High Risk (%)", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(1985:end_yr,plot_this)

plot_this <- ((results[(1985-1949):(end_yr-1949),"TLTBI_INITS_FB"]) / (results[(1985-1949):(end_yr-1949),"TLTBI_INITS"]))*100
plot(0,0,ylim=range(plot_this),xlim=c(1985,end_yr),xlab="Year",ylab="IPT Treatment Initiations, Foreign Born (%)", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
lines(1985:end_yr,plot_this)
}
# LTBI prevalence by age US
plot_this <-cbind(results[(start_yr-1949):(end_yr-1949),55:65])
plot_this <-colSums(plot_this)/sum(plot_this)*100
plot(0,0,ylim=range(plot_this),xlim=c(0.6,11.4),xlab="Age Group",ylab="US-Born LTBI Prevalance by Age", axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
lines(plot_this)
#LTBI prevalence by age non-US
plot_this <-cbind(results[(start_yr-1949):(end_yr-1949),66:76])
plot_this <-colSums(plot_this)/sum(plot_this)*100
plot(0,0,ylim=range(plot_this),xlim=c(0.6,11.4),xlab="Age Group",ylab="Non-US-Born LTBI Prevalance by Age", axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
lines(plot_this)

# TB deaths
plot_this <- cbind(results[(start_yr-1949):(end_yr-1949),"TBMORT_US"]+results[(start_yr-1949):(end_yr-1949),"TBMORT_NUS"])
plot(0,0,ylim=range(plot_this),xlim=c(start_yr,end_yr),xlab="Year",ylab="Total TB Deaths by Year", axes=F)
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
plot(start_yr:end_yr,plot_this)

# TB deaths by age
plot_this <-cbind(results[(start_yr-1949):(end_yr-1949),88:98]+results[(start_yr-1949):(end_yr-1949),99:109])
plot_this <- colSums(plot_this)
plot(0,0,ylim=range(plot_this),xlim=c(0.6,11.4),xlab="Age Group",ylab="Total TB Deaths by Age Group", axes=F)
axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
axis(1,1:11,paste(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-94","95+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
lines(plot_this)


  dev.off();
  #system(paste("open", pdfname))
}

