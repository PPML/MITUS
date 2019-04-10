#'loads data for the geography of interest
#'@name optim_load
#'@return env
optim_load<-function(){
# create a logical vector for all 51 geographies
  # optim_vec<-rep("F",nrow(stateID))
#create an empty vector for the locations
  loc_vec<-vector()
#loop through all the geographies to check for optim files
for (i in 1:nrow(stateID)){
  loc<-stateID[i,3]
  if (any(grepl(paste0(loc,"_Optim_all"),list.files(paste0("inst/",loc,"/"), full.names = TRUE)))==TRUE){
    # optim_vec[i]<-"T"
    loc_vec<-c(loc_vec,loc)
  }
}
if (length(loc_vec)>1){
  loc_vec<-c("US",loc_vec)
}
return(loc_vec)
}
