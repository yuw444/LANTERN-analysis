test_that("count_ancestry_codes works", {
  mat <- matrix(c(3, 3, 1, 2, 1, 3, 2, 1, 1), nrow = 3, ncol = 3)
  result <- count_ancestry_codes(mat, 3)
  expect_equal(result[1], 2)
  expect_equal(result[2], 0)
  expect_equal(result[3], 1)
})

test_that("split_by_ancestry works", {
  gt <- matrix(c(1, 2, 0, 1, 0, 2), nrow = 2, ncol = 3)
  an <- matrix(c(3, 1, 2, 2, 1, 3), nrow = 2, ncol = 3)
  
  result <- split_by_ancestry(gt, an)
  
  expect_type(result, "list")
  expect_named(result, c("african", "european"))
})
