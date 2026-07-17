test_that("write_ancestry_gds writes one GDS per population", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  split_result <- list(
    AFR = matrix(c(1, 0, 2, 1), nrow = 2, ncol = 2),
    EUR = matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2),
    variant_info = data.frame(chrom = c("chr19", "chr19"), pos = c(100L, 200L),
                               ref = c("A", "G"), alt = c("T", "C"),
                               stringsAsFactors = FALSE),
    sample_ids = c("S1", "S2"),
    mode = "dosage",
    ancestry_counts = matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                              dimnames = list(NULL, c("AFR", "EUR"))),
    overlap = list()
  )

  gds_paths <- write_ancestry_gds(split_result, td, verbose = FALSE)

  expect_named(gds_paths, c("AFR", "EUR"))
  expect_true(file.exists(gds_paths$AFR))
  expect_true(file.exists(gds_paths$EUR))

  g <- SeqArray::seqOpen(gds_paths$AFR)
  on.exit(SeqArray::seqClose(g), add = TRUE)
  expect_equal(SeqArray::seqGetData(g, "sample.id"), c("S1", "S2"))
  expect_equal(length(SeqArray::seqGetData(g, "position")), 2)
})

test_that("write_ancestry_gds creates out_path if missing", {
  td <- file.path(tempfile(), "nested", "out")
  on.exit(unlink(dirname(dirname(td)), recursive = TRUE), add = TRUE)
  expect_false(dir.exists(td))

  split_result <- list(
    POP1 = matrix(c(1, 2), nrow = 1, ncol = 2),
    variant_info = data.frame(chrom = "chr1", pos = 500L, ref = "A", alt = "G",
                               stringsAsFactors = FALSE),
    sample_ids = c("S1", "S2")
  )

  gds_paths <- write_ancestry_gds(split_result, td, verbose = FALSE)
  expect_true(dir.exists(td))
  expect_true(file.exists(gds_paths$POP1))
})

test_that("write_ancestry_gds writes one GDS per population for K > 2", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  pop_names <- c("AFR", "EUR", "NAT", "EAS", "SAS")   # K = 5
  split_result <- setNames(
    lapply(seq_along(pop_names), function(i) matrix(i, nrow = 2, ncol = 2)),
    pop_names
  )
  split_result$variant_info <- data.frame(
    chrom = c("chr1", "chr1"), pos = c(100L, 200L),
    ref = c("A", "G"), alt = c("T", "C"), stringsAsFactors = FALSE
  )
  split_result$sample_ids <- c("S1", "S2")

  gds_paths <- write_ancestry_gds(split_result, td, verbose = FALSE)

  expect_named(gds_paths, pop_names)
  for (pop in pop_names) expect_true(file.exists(gds_paths[[pop]]))
})
