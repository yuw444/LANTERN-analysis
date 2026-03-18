context("Core C functions")

test_that("count_ancestry_codes works", {
  mat <- matrix(c(3, 3, 1, 2, 1, 3, 2, 1, 1), nrow = 3, ncol = 3, byrow = TRUE)
  result <- count_ancestry_codes(mat, 3)
  
  expect_type(result, "integer")
  expect_length(result, 3)
  expect_gt(result[1], 0)
})

test_that("count_ancestry_codes with different codes", {
  mat <- matrix(c(3, 3, 1, 2, 1, 3), nrow = 3, ncol = 2, byrow = TRUE)
  
  count_3 <- count_ancestry_codes(mat, 3)
  count_1 <- count_ancestry_codes(mat, 1)
  
  expect_gt(count_3[1], count_1[1])
  expect_equal(count_3[1] + count_1[1], 2)
})

test_that("split_by_ancestry returns correct structure", {
  gt <- matrix(c(1, 2, 0, 1, 0, 2), nrow = 2, ncol = 3, byrow = TRUE)
  an <- matrix(c(3, 1, 2, 2, 1, 3), nrow = 2, ncol = 3, byrow = TRUE)
  
  result <- split_by_ancestry(gt, an)
  
  expect_type(result, "list")
  expect_named(result, c("african", "european"))
  expect_equal(dim(result$african), c(2, 3))
  expect_equal(dim(result$european), c(2, 3))
})

test_that("split_by_ancestry handles pure ancestries", {
  gt <- matrix(c(1, 0, 2, 0), nrow = 2, ncol = 2)
  an <- matrix(c(3, 3, 1, 1), nrow = 2, ncol = 2)
  
  result <- split_by_ancestry(gt, an)
  
  expect_equal(result$african[1, 1], 2)
  expect_equal(result$european[1, 1], 0)
  expect_equal(result$african[2, 1], 0)
  expect_equal(result$european[2, 1], 0)
})

test_that("split_by_ancestry handles mixed ancestry", {
  gt <- matrix(c(1, 2, 0), nrow = 1, ncol = 3)
  an <- matrix(c(2, 2, 2), nrow = 1, ncol = 3)
  
  result <- split_by_ancestry(gt, an)
  
  expect_equal(result$african[1, 1], result$european[1, 1])
})

test_that("NA and invalid codes handled", {
  gt <- matrix(c(0, 1, 2, 0), nrow = 2, ncol = 2)
  an <- matrix(c(3, 0, 1, 4), nrow = 2, ncol = 2)
  
  result <- split_by_ancestry(gt, an)
  
  expect_equal(result$african[1, 1], 0)
  expect_equal(result$european[2, 1], 0)
})
