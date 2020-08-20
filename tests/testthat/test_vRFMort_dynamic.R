test_that("mortality increases over RF level", {
  prg_chng<-def_prgchng(P)
  prms <-list()
  prms <- param_init(P,"US",prg_chng=prg_chng, ttt_list=def_ttt())
  for (i in 2:4){
    expect_true(prms$RRmuRF[i-1] < prms$RRmuRF[i])
  }
})
