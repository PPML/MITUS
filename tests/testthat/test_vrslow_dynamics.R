test_that("probability of slow progression decreases as RF level increases", {
  for (i in 2:4){
    expect_true(Vrslow[i] > Vrslow[i-1])
  }
})
