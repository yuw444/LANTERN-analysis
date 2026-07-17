test_that("cauchy_combine returns a valid p-value in (0, 1]", {
  p <- cauchy_combine(c(0.01, 0.04))
  expect_type(p, "double")
  expect_true(p > 0 && p <= 1)
})

test_that("cauchy_combine is symmetric in input order", {
  p1 <- cauchy_combine(c(0.01, 0.5, 0.2))
  p2 <- cauchy_combine(c(0.2, 0.01, 0.5))
  expect_equal(p1, p2, tolerance = 1e-9)
})

test_that("cauchy_combine of identical p-values returns approximately that p-value", {
  p <- cauchy_combine(c(0.05, 0.05, 0.05))
  expect_equal(p, 0.05, tolerance = 1e-6)
})

test_that("cauchy_combine supports weights and matches equal-weight default", {
  p_equal    <- cauchy_combine(c(0.01, 0.04))
  p_weighted <- cauchy_combine(c(0.01, 0.04), weights = c(1, 1))
  expect_equal(p_equal, p_weighted, tolerance = 1e-9)

  # Heavier weight on the smaller p-value should pull the combined p-value down
  p_skewed <- cauchy_combine(c(0.01, 0.5), weights = c(10, 1))
  p_flat   <- cauchy_combine(c(0.01, 0.5), weights = c(1, 1))
  expect_true(p_skewed < p_flat)
})

test_that("cauchy_combine validates inputs", {
  expect_error(cauchy_combine(c(0.5)), "length >= 2")
  expect_error(cauchy_combine(c(0.5, NA)), "contains NA")
  expect_error(cauchy_combine(c(0.5, 1.5)), "must be in \\(0, 1\\]")
  expect_error(cauchy_combine(c(0.5, -0.1)), "must be in \\(0, 1\\]")
  expect_error(cauchy_combine(c(0.1, 0.2), weights = c(1, 1, 1)),
               "same length")
  expect_error(cauchy_combine(c(0.1, 0.2), weights = c(-1, 1)),
               "non-negative")
})

test_that("cauchy_combine clips exact zeros instead of erroring", {
  p <- cauchy_combine(c(0, 0.5))
  expect_true(is.finite(p))
  expect_true(p > 0 && p <= 1)
})
