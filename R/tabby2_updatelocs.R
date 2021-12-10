vector_all_locs<-c("CA", "MA", "FL", "GA", "IL", "TX", "NY", "NJ", "PA", "WA", "VA" )

vector_all_locs<-ordered_locs[-37]
vector_all_locs<- vector_all_locs[-48]; vector_all_locs <- vector_all_locs[-44]
#the simp.date here is to append to the calibration plots;
#the loaded optim data is not dependent on the simp.date
#it is pulled directly from the model_load function.
update_all_states<-function(loc_vec, simp.date){
  for (i in 1:length(loc_vec)){
    loc<-loc_vec[i]
    print(loc)
    model_load(loc)
    optim_to_tabby2(loc, simp.date)
  }
}

# update_all_states(vector_all_locs, "1201")

# for (loc in vector_all_locs[26:48]){
#   dir.create(paste0('~/MITUS/inst/', loc, "/calibration_outputs"))
#
# }

