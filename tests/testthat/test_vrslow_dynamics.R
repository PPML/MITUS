test_that("probability of slow progression decreases as RF level increases", {
  model_load()
  prg_chng<-def_prgchng(P)
  prms <-list()
  prms <- param_init(P,"US",prg_chng=prg_chng, ttt_list=def_ttt())
  for (i in 1:11){
    for (j in 2:4){
      expect_true(prms$Mrslow[i,j] > prms$Mrslow[i,j-1] )
    } }
})
