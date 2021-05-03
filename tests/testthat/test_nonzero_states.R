test_that("no model state has a negative population", {
  for (i in 1:12672){
  expect_true(v1[i] != 0 )
  }
})
