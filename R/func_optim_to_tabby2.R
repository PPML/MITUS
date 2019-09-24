#' THIS FUNCTION READS IN THE OPTIM FILE THAT MATCHES THE SIMP.DATE
#' INPUT AND CREATES ALL NECESSARY FILES FOR TABBY2 FOR THAT LOCATION.
#' ALL PREDEFINED SCENARIOS (RESHAPED AND ORIGINAL FILES) AND MODEL
#' CALIBRATION OUTPUTS. IT ALSO REMOVES ALL OLD CALIBRATION OUTPUTS.

#'@name optim_to_tabby2
#'@param loc two character location code for the model location
#'@param simp.date date to append to the files; should match the optim date
#'@export

optim_to_tabby2<-function(loc, simp.date="724"){
  #load the optim matrix from the file
  params<-readRDS(system.file(paste0(loc,"/", loc, "_Param_all_10_",simp.date,".rds"), package="MITUS"))
  #remove the posterior value
  # optims<-optim[,-ncol(optim)]
  #create a transformed parameter vector from this matrix
  #generate the results vectors that we need
  results.list<-make_all_scenarios(loc,params)
  #reshape those results
  reshape_results(loc,results.list)

  #if you are reshaping old saved results use this workflow
  # load_data <- function(i) {
  #   data_name <-
  #     load(system.file(paste0(loc, "/", loc, "_results_",i,".rda"), package='MITUS'))
  #   return(get(data_name))
  # }
  # results.list <- lapply(1:9, load_data)
  # reshape_results(loc,results.list)

  #remove old calibration plot data
  #get the files in the directory
  file.list<-list.files(system.file(paste0(loc,"/calibration_outputs/"),package="MITUS"))
  if (length(file.list)>0){
  for (i in 1:length(file.list)){
    file.remove(paste0("~/MITUS/inst/",loc,"/calibration_outputs/", file.list[i]))
  }}
  #create new outputs for calibration plots
    #get the data array for the base case
    base.case<-results.list[[1]]
    #if you are reshaping old saved results use this workflow
    # load(system.file(paste0(loc,"/",loc,"_results_1.rda"),package="MITUS"))
    # model_calib_outputs(loc,out,samp_i=1,simp.date = "724")
    model_calib_outputs(loc,base.case,samp_i=1,simp.date =simp.date)
}
