def_care_cascade<-function(){
  care_cascade<-list("ttt_ltbi_init" = 0,
                     "ttt_ltbi_comp" = 0,
                     "ttt_ltbi_eff" = 0,
                     "ttt_ltbi_sens" = rep(0,2),
                     "ttt_ltbi_spec" = rep(0,2),
                     "ttt_ltbi_accept" =0)
  return(care_cascade)
  }
