test_that("state simulations run", {
  for (i in 1:length(ordered_locs)){
  loc<-ordered_locs[i]
  model_load(loc)

  testthat::expect_error(
    OutputsZint(samp_i = 1, ParMatrix=P,
              loc = loc, startyr = 1950,
              endyr = 2020, prg_chng = def_prgchng(P),
              ttt_list = def_ttt())
    , NA)

  }
})
