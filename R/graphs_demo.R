#' creates plots of the demographics of the model
#' and writes these results to a .pdf file.
#' Use to check both with and without tb in the model


#'Create a function to be run on a specific model run output to
#'create simple graphs of all the output for a selected year range
#'@param start_yr year to start the graphs
#'@param end_yr year to end the graphs
#'@param df dataframe of output for all years
#'@return .pdf of the graphs
#'@export

tb_graph_demo <- function(start_yr, end_yr, df){

  pdf(file=paste("MITUS_results/graphs_all",Sys.time(),".pdf"), width = 11, height = 8.5)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(df[1:66,30], df[1:66,31]+df[1:66,32])
  plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")
  # points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],pch=19,cex=0.6,col="grey50")
  # points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],pch=19,cex=0.6,col="blue")
  # points(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],pch=19,cex=0.6,col="red3")
  # lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,2],lty=3,col="grey50")
  # lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,3],lty=3,col="blue")
  # lines(tot_pop_yr_fb[,1],tot_pop_yr_fb[,4],lty=3,col="red3")
  lines(1950:2015,V[,2],lwd=2,col="red3")
  lines(1950:2015,V[,1],lwd=2,col="blue")
  lines(1950:2015,rowSums(V),lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Population: Total, US, and Foreign Born (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### TOTAL POP AGE DISTRIBUTION 2014  ### ### ### ### ### ###
  ### need to update
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  V  <- cbind(t(df[65,33:43]), t(df[65,44:54]))

  plot(0,1,ylim=c(0.05,135),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y")
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")
  # points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  # points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")
  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Population by Age for FB (red) and US (blue), 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","Fitted model"),pch=c(19,15),pt.cex=c(1,2),
         lwd=NA,col=c("grey30","grey80"),bg="white")
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###   TOTAL MORT EACH DECADE, BY US/FB  ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(rowsums(df[1:66,154:164]), rowsums(df[1:66,165:175]))
  plot(1,1,ylim=c(2,500),xlim=c(1950,2015),xlab="",ylab="",axes=F,log="y")
  axis(1);axis(2,las=2);box()
  abline(h=axTicks(2),col="grey85")

  lines(1950:2015,V[,2],lwd=2,col="red3")
  lines(1950:2015,V[,1],lwd=2,col="blue")
  lines(1950:2015,rowSums(V),lwd=2,col="grey50")

  mtext("Year",1,2.5,cex=0.9)
  mtext("Mortality: Total, US, and Foreign Born (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("bottomright",c("Total","US born","Foreign born","Reported data","Fitted model"),cex=0.9,
         pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ###  TOTAL MORT AGE DISTRIBUTION 2014   ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  V  <- cbind(t(df[65,154:164]), t(df[65,165:175]))

  plot(0,1,ylim=c(0.05,135),xlim=c(0.6,8.4),xlab="",ylab="",axes=F,col=NA,log="y")
  axis(1,1:8,paste(c("0-4","5-24","25-44","45-54","55-64","65-74","75-84","85+"),"\nyears",sep=""),tick=F,cex.axis=0.75)
  axis(1,1:9-0.5,rep("",9))
  axis(2,c(0.1,1,10,100),las=2);box()
  abline(h=axTicks(2),col="grey85")
  for(i in 1:8) polygon(i+c(.4,0,0,.4),c(0.0001,0.0001,V[i,1],V[i,1]),border=NA,col="lightblue")
  for(i in 1:8) polygon(i+c(-.4,0,0,-.4),c(0.0001,0.0001,V[i,2],V[i,2]),border=NA,col="pink")
  # points(1:8+0.2,CalibDat[["tot_pop14_ag_fb"]][-9,3],pch=19,cex=1.2,col="blue")
  # points(1:8-0.2,CalibDat[["tot_pop14_ag_fb"]][-9,4],pch=19,cex=1.2,col="red3")
  mtext("Age Group",1,2.5,cex=0.9)
  box()
  mtext("Mortality by Age for FB (red) and US (blue), 2014 (mil, log-scale)",3,.8,font=2,cex=0.8)
  legend("topright",c("Reported data","Fitted model"))

  }
