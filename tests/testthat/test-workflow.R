context("Workflow functions - overlap handling")

test_that("run_ancestry_pipeline handles sample overlap", {
  # GT: 3 samples (A, B, C)
  # PT: 3 samples (A, B, D)
  # Expected: only A and B are common
  gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1), nrow = 3, ncol = 3,
               dimnames = list(c("var1", "var2", "var3"),
                               c("sample_A", "sample_B", "sample_C")))

  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
               dimnames = list(c("sample_A", "sample_B", "sample_D"),
                               c("var1", "var2")))

  result <- run_ancestry_pipeline(gt, pt, verbose = FALSE)

  # Check overlap info
  expect_type(result, "list")
  expect_true("overlap" %in% names(result))
  expect_equal(result$overlap$n_samples_kept, 2)
  expect_true("sample_C" %in% result$overlap$dropped_samples)
  expect_true("sample_D" %in% result$overlap$dropped_samples)

  # Check dimensions after filtering
  expect_equal(ncol(result$african), 2)  # 2 common samples
  expect_equal(nrow(result$african), 2)  # 2 common variants
})

test_that("run_ancestry_pipeline handles coordinate-based variant matching", {
  # GT variants: chr22:100, chr22:200, chr22:300
  # PT regions: chr22:50-150, chr22:150-250
  # Expected: chr22:100 matches chr22:50-150, chr22:200 matches chr22:150-250

  gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1), nrow = 3, ncol = 3,
               dimnames = list(c("chr22:100", "chr22:200", "chr22:300"),
                               c("s1", "s2", "s3")))

  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
               dimnames = list(c("s1", "s2", "s3"),
                               c("chr22:50-150", "chr22:150-250")))

  result <- run_ancestry_pipeline(gt, pt, verbose = FALSE)

  # Only 2 variants should match
  expect_equal(result$overlap$n_variants_kept, 2)
  expect_equal(nrow(result$african), 2)
  expect_equal(ncol(result$african), 3)
})

test_that("run_ancestry_pipeline returns correct counts", {
  gt <- matrix(c(2, 1, 0, 1), nrow = 2, ncol = 2,
               dimnames = list(c("var1", "var2"),
                               c("s1", "s2")))

  pt <- matrix(c(3, 1, 2, 2), nrow = 2, ncol = 2,
               dimnames = list(c("s1", "s2"),
                               c("var1", "var2")))

  result <- run_ancestry_pipeline(gt, pt, verbose = FALSE)

  expect_true("counts" %in% names(result))
  expect_true("african" %in% names(result$counts))
  expect_true("european" %in% names(result$counts))
  expect_true("mixed" %in% names(result$counts))
})

test_that("run_ancestry_pipeline handles unnamed matrices", {
  # Unnamed matrices - should match by order
  gt <- matrix(c(2, 1, 0, 1, 2, 1), nrow = 2, ncol = 3)
  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2)

  expect_warning(result <- run_ancestry_pipeline(gt, pt, verbose = FALSE),
                 "different dimensions")

  expect_equal(result$overlap$n_samples_kept, 3)
})

test_that("run_ancestry_pipeline stops on no common samples", {
  gt <- matrix(0, nrow = 2, ncol = 2,
               dimnames = list(c("v1", "v2"), c("s1", "s2")))

  pt <- matrix(1, nrow = 2, ncol = 2,
               dimnames = list(c("s3", "s4"), c("v1", "v2")))

  expect_error(run_ancestry_pipeline(gt, pt, verbose = FALSE),
               "No common samples found")
})

test_that("result contains all expected elements", {
  gt <- matrix(c(2, 1, 0, 1), nrow = 2, ncol = 2,
               dimnames = list(c("v1", "v2"), c("s1", "s2")))

  pt <- matrix(c(3, 1, 2, 2), nrow = 2, ncol = 2,
               dimnames = list(c("s1", "s2"), c("v1", "v2")))

  result <- run_ancestry_pipeline(gt, pt, verbose = FALSE)

  expect_named(result, c("african", "european", "counts", "overlap"))
  expect_equal(dim(result$african), c(2, 2))
  expect_equal(dim(result$european), c(2, 2))
})

test_that("run_ancestry_pipeline filters monomorphic variants", {
  # GT: var1 has alt alleles, var2 is monomorphic (all 0)
  gt <- matrix(c(2, 1, 0, 0, 0, 0), nrow = 2, ncol = 3,
               dimnames = list(c("var1", "var2"),
                               c("s1", "s2", "s3")))

  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
               dimnames = list(c("s1", "s2", "s3"),
                               c("var1", "var2")))

  result <- run_ancestry_pipeline(gt, pt, verbose = FALSE)

  # var2 should be filtered out
  expect_equal(result$overlap$n_monomorphic_filtered, 1)
  expect_equal(result$overlap$n_variants_kept, 1)
  expect_equal(nrow(result$african), 1)
})

test_that("run_ancestry_pipeline returns empty result if all variants are monomorphic", {
  gt <- matrix(0, nrow = 3, ncol = 2,
               dimnames = list(c("v1", "v2", "v3"),
                               c("s1", "s2")))

  pt <- matrix(1, nrow = 2, ncol = 3,
               dimnames = list(c("s1", "s2"),
                               c("v1", "v2", "v3")))

  expect_warning(result <- run_ancestry_pipeline(gt, pt, verbose = FALSE),
                 "All variants are monomorphic")

  expect_equal(result$overlap$n_variants_kept, 0)
  expect_equal(nrow(result$african), 0)
})
