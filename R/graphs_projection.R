####Generalized Projection plots
#'@name projection_plots
#'@param resultsList list of results from n model runs (names included)
#'@param loc two digit location code of results
#'@export
projection_plots<-function(resultsList, loc){
index<-length(resultsList)
endyr_index<-nrow(resultsList[[1]])
endyr<-1949+endyr_index
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ###   TOTAL POP EACH DECADE, BY US/FB   ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#set up the plot
plot(1,1,ylim=c(2,500),xlim=c(1950,endyr),xlab="",ylab="",axes=F,log="y")
axis(1);axis(2,las=2);box()
abline(h=axTicks(2),col="grey85")
#plot the results for each run
for (i in 1:length(resultsList)){
V  <- cbind(df[1:endyr_index,30], df[1:endyr_index,31]+df[1:endyr_index,32])


lines(1950:endyr,V[,2],lwd=2,col="red3")
lines(1950:endyr,V[,1],lwd=2,col="blue")
lines(1950:endyr,rowSums(V),lwd=2,col="grey50")
}
mtext("Year",1,2.5,cex=1.2)
mtext("Population: Total, US, and Non-US Born (mil, log-scale)",3,.8,font=2,cex=1)
legend("bottomright",c("Total","US born","Non-US Born","Reported data","model"),cex=1,
       pch=c(15,15,15,19,NA),lwd=c(NA,NA,NA,1,2),lty=c(NA,NA,NA,3,1),col=c("grey50",4,"red3",1,1),bg="white",pt.cex=c(1.8,1.8,1.8,0.3,NA))

}
