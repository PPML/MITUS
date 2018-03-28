test_that("no output is non-positive", {
  for (i in 1:nYrs){
    for (j in 1:nRes){
    expect_true(m[i,j]< 0 )
  }}
})
