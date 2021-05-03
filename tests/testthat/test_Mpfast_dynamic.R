test_that("probability of fast progression increases as RF level increases", {
  for (i in 1:11){
  for (j in 2:4){
    expect_true(Mpfast[i,j] > Mpfast[i,j-1] )
  } }
})
