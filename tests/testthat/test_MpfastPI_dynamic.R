test_that("rate of fast progression w/PI increases as RF level increases", {
  for (i in 1:11){
    for (j in 2:4){
      expect_true(MpfastPI[i,j] > MpfastPI[i,j-1] )
    } }
})
