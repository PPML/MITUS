vector_all_locs<-c("CA", "MA", "FL", "GA", "IL", "TX", "NY", "NJ", "PA", "WA", "VA")

update_all_states<-function(loc_vec, simp.date){
  for (i in 1:length(loc_vec)){
    loc<-loc_vec[i]
    model_load(loc)
    optim_to_tabby2(loc, simp.date)
  }
}
