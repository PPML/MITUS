#'@name format_results
#'@param startyr year results start
#'@param endyr year results end
#'@param ParMatrix this is a matrix of P's corresponding to optim
#'@return an array of the results
#'@export
format_results <- function(ParMatrix, startyr,endyr) {
  array_results<-array(NA,dim=c(nrow(ParMatrix),endyr-startyr,length(IP[["ResNam"]])))
  for(i in 1:nrow(ParMatrix)){
    array_results[i,,]<-OutputsZint(ParMatrix = ParMatrix)
  }
  assign(paste("array_results_",batch,sep = ""),array_results)
  save(list=paste("array_results_",batch,sep = ""),file=paste("array_results_",batch,"_",Sys.Date(),sep = ""))
  }
