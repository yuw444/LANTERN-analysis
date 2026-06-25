context("Phased C functions")

test_that("split_phased deterministically splits phased haplotypes", {
  gt0 <- matrix(c(1L, 0L, 0L, 1L), 2, 2)
  gt1 <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
  a0 <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
  a1 <- matrix(c(1L, 0L, 1L, 0L), 2, 2)

  result <- split_phased(gt0, gt1, a0, a1)

  expect_equal(result$african, matrix(c(1, 1, 0, 0), 2, 2))
  expect_equal(result$european, matrix(c(0, 0, 1, 1), 2, 2))
  expect_equal(dim(result$african), c(2, 2))
  expect_equal(dim(result$european), c(2, 2))
})

test_that("split_phased handles homozygous alt and mixed haplotypes", {
  afr <- split_phased(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(0L, 1, 1), matrix(0L, 1, 1))
  eur <- split_phased(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(1L, 1, 1))
  mixed1 <- split_phased(matrix(1L, 1, 1), matrix(0L, 1, 1), matrix(0L, 1, 1), matrix(1L, 1, 1))
  mixed2 <- split_phased(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(0L, 1, 1))

  expect_equal(afr$african, matrix(2, 1, 1))
  expect_equal(afr$european, matrix(0, 1, 1))
  expect_equal(eur$african, matrix(0, 1, 1))
  expect_equal(eur$european, matrix(2, 1, 1))
  expect_equal(mixed1$african, matrix(1, 1, 1))
  expect_equal(mixed1$european, matrix(0, 1, 1))
  expect_equal(mixed2$african, matrix(1, 1, 1))
  expect_equal(mixed2$european, matrix(1, 1, 1))
})

test_that("split_phased treats NA gt and invalid ancestry as zero", {
  gt0 <- matrix(c(NA, 1), 1, 2)
  gt1 <- matrix(c(1, 0), 1, 2)
  a0 <- matrix(c(0L, 1L), 1, 2)
  a1 <- matrix(c(1L, 1L), 1, 2)

  result <- split_phased(gt0, gt1, a0, a1)

  expect_equal(result$african, matrix(c(0, 0), 1, 2))
  expect_equal(result$european, matrix(c(1, 1), 1, 2))

  invalid <- split_phased(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(9L, 1, 1), matrix(8L, 1, 1))
  expect_equal(invalid$african, matrix(0, 1, 1))
  expect_equal(invalid$european, matrix(0, 1, 1))
})

test_that("split_phased supports custom pop codes and default pop codes", {
  gt0 <- matrix(1L, 1, 1)
  gt1 <- matrix(1L, 1, 1)
  a0 <- matrix(2L, 1, 1)
  a1 <- matrix(5L, 1, 1)

  custom <- split_phased(gt0, gt1, a0, a1, pop_codes = c(AFR = 2L, EUR = 5L))
  default <- split_phased(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(0L, 1, 1), matrix(1L, 1, 1))

  expect_equal(custom$african, matrix(1, 1, 1))
  expect_equal(custom$european, matrix(1, 1, 1))
  expect_equal(default$african, matrix(1, 1, 1))
  expect_equal(default$european, matrix(1, 1, 1))
})
