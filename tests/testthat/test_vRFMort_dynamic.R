test_that("mortality increases over RF level", {
  for (i in 1:4){
    expect_true(vRFMort[i-1] < vRFMort[i])
  }
})
