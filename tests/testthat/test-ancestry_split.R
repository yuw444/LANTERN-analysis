# ============================================================================
# Core C primitives: count_ancestry_codes, split_diploid, split_diploid_multi
# ============================================================================

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

test_that("split_diploid returns correct structure", {
  gt <- matrix(c(1, 2, 0, 1, 0, 2), nrow = 2, ncol = 3, byrow = TRUE)
  an <- matrix(c(3, 1, 2, 2, 1, 3), nrow = 2, ncol = 3, byrow = TRUE)

  result <- split_diploid(gt, an)

  expect_type(result, "list")
  expect_named(result, c("african", "european"))
  expect_equal(dim(result$african), c(2, 3))
  expect_equal(dim(result$european), c(2, 3))
})

test_that("split_diploid handles pure ancestries", {
  gt <- matrix(c(1, 0, 2, 0), nrow = 2, ncol = 2)
  an <- matrix(c(3, 3, 1, 1), nrow = 2, ncol = 2)

  result <- split_diploid(gt, an)

  # var1, s1: gt=1, anc=3 (pure AFR) → african=1, european=0
  expect_equal(result$african[1, 1], 1)
  expect_equal(result$european[1, 1], 0)
  # var1, s2: gt=2, anc=1 (pure EUR) → african=0, european=2
  expect_equal(result$african[1, 2], 0)
  expect_equal(result$european[1, 2], 2)
  expect_equal(result$african[2, 1], 0)
  expect_equal(result$european[2, 1], 0)
})

test_that("split_diploid handles homozygous alt mixed ancestry", {
  gt <- matrix(c(2, 2, 2), nrow = 1, ncol = 3)
  an <- matrix(c(2, 2, 2), nrow = 1, ncol = 3)

  result <- split_diploid(gt, an)

  expect_equal(result$african[1, 1], 1.0)
  expect_equal(result$european[1, 1], 1.0)
})

test_that("split_diploid p1/p2 calculation with known values", {
  # For a variant with:
  # - 1 pure African het (pt=3, gt=1) -> N1=0, N2=1
  # - 1 mixed het (pt=2, gt=1) -> N4=0, N5=1
  # - 1 pure European het (pt=1, gt=1) -> N7=0, N8=1
  #
  # total_alt = 0*2 + 1 + 0 + 1 + 0*2 + 1 = 3
  # p1 = (0 + 1 + 0) / (3 - 1) = 0.5
  # p2 = (0 + 0 + 1) / (3 - 1) = 0.5

  gt <- matrix(c(1, 1, 1), nrow = 1, ncol = 3)
  an <- matrix(c(3, 2, 1), nrow = 1, ncol = 3)

  result <- split_diploid(gt, an)

  expect_equal(result$african[1, 1], 1.0)  # Pure African het
  expect_equal(result$european[1, 3], 1.0)  # Pure European het
  expect_equal(result$african[1, 2], 0.5)   # Mixed het -> p1 = 0.5
  expect_equal(result$european[1, 2], 0.5)  # Mixed het -> p2 = 0.5
})

test_that("split_diploid singleton case", {
  # Singleton: only one het with mixed ancestry
  gt <- matrix(c(1), nrow = 1, ncol = 1)
  an <- matrix(c(2), nrow = 1, ncol = 1)

  result <- split_diploid(gt, an)

  expect_equal(result$african[1, 1], 0.5)
  expect_equal(result$european[1, 1], 0.5)
})

test_that("split_diploid non-singleton mixed p1 calculation", {
  # For a variant with:
  # - 2 pure African het (pt=3, gt=1) -> N1=0, N2=2
  # - 1 mixed het (pt=2, gt=1) -> N4=0, N5=1
  # - 0 pure European het (pt=1, gt=1) -> N7=0, N8=0
  #
  # total_alt = 0*2 + 2 + 0 + 1 + 0*2 + 0 = 3
  # p1 = (0 + 2 + 0) / (3 - 1) = 1.0
  # p2 = (0 + 0 + 0) / (3 - 1) = 0.0

  gt <- matrix(c(1, 1, 1), nrow = 1, ncol = 3)
  an <- matrix(c(3, 3, 2), nrow = 1, ncol = 3)

  result <- split_diploid(gt, an)

  expect_equal(result$african[1, 1], 1.0)
  expect_equal(result$african[1, 2], 1.0)
  expect_equal(result$african[1, 3], 1.0)   # Mixed het -> p1 = 1.0
  expect_equal(result$european[1, 3], 0.0)   # Mixed het -> p2 = 0.0
})

test_that("NA and invalid codes handled", {
  gt <- matrix(c(0, 1, 2, 0), nrow = 2, ncol = 2)
  an <- matrix(c(3, 0, 1, 4), nrow = 2, ncol = 2)

  result <- split_diploid(gt, an)

  expect_equal(result$african[1, 1], 0)
  expect_equal(result$european[2, 1], 0)
})

test_that("split_diploid homozygous alt pure ancestries", {
  gt <- matrix(c(2, 2), nrow = 1, ncol = 2)
  an <- matrix(c(3, 1), nrow = 1, ncol = 2)

  result <- split_diploid(gt, an)

  expect_equal(result$african[1, 1], 2.0)  # Homozygous alt with African
  expect_equal(result$european[1, 2], 2.0)  # Homozygous alt with European
})

# ============================================================================
# split_haplotype / split_haplotype_multi
# ============================================================================

test_that("split_haplotype deterministically splits phased haplotypes", {
  gt0 <- matrix(c(1L, 0L, 0L, 1L), 2, 2)
  gt1 <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
  a0 <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
  a1 <- matrix(c(1L, 0L, 1L, 0L), 2, 2)

  result <- split_haplotype(gt0, gt1, a0, a1)

  expect_equal(result$african, matrix(c(1, 1, 0, 0), 2, 2))
  expect_equal(result$european, matrix(c(0, 0, 1, 1), 2, 2))
  expect_equal(dim(result$african), c(2, 2))
  expect_equal(dim(result$european), c(2, 2))
})

test_that("split_haplotype handles homozygous alt and mixed haplotypes", {
  afr <- split_haplotype(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(0L, 1, 1), matrix(0L, 1, 1))
  eur <- split_haplotype(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(1L, 1, 1))
  mixed1 <- split_haplotype(matrix(1L, 1, 1), matrix(0L, 1, 1), matrix(0L, 1, 1), matrix(1L, 1, 1))
  mixed2 <- split_haplotype(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(0L, 1, 1))

  expect_equal(afr$african, matrix(2, 1, 1))
  expect_equal(afr$european, matrix(0, 1, 1))
  expect_equal(eur$african, matrix(0, 1, 1))
  expect_equal(eur$european, matrix(2, 1, 1))
  expect_equal(mixed1$african, matrix(1, 1, 1))
  expect_equal(mixed1$european, matrix(0, 1, 1))
  expect_equal(mixed2$african, matrix(1, 1, 1))
  expect_equal(mixed2$european, matrix(1, 1, 1))
})

test_that("split_haplotype treats NA gt and invalid ancestry as zero", {
  gt0 <- matrix(c(NA, 1), 1, 2)
  gt1 <- matrix(c(1, 0), 1, 2)
  a0 <- matrix(c(0L, 1L), 1, 2)
  a1 <- matrix(c(1L, 1L), 1, 2)

  result <- split_haplotype(gt0, gt1, a0, a1)

  expect_equal(result$african, matrix(c(0, 0), 1, 2))
  expect_equal(result$european, matrix(c(1, 1), 1, 2))

  invalid <- split_haplotype(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(9L, 1, 1), matrix(8L, 1, 1))
  expect_equal(invalid$african, matrix(0, 1, 1))
  expect_equal(invalid$european, matrix(0, 1, 1))
})

test_that("split_haplotype supports custom pop codes and default pop codes", {
  gt0 <- matrix(1L, 1, 1)
  gt1 <- matrix(1L, 1, 1)
  a0 <- matrix(2L, 1, 1)
  a1 <- matrix(5L, 1, 1)

  custom <- split_haplotype(gt0, gt1, a0, a1, pop_codes = c(AFR = 2L, EUR = 5L))
  default <- split_haplotype(matrix(1L, 1, 1), matrix(1L, 1, 1), matrix(0L, 1, 1), matrix(1L, 1, 1))

  expect_equal(custom$african, matrix(1, 1, 1))
  expect_equal(custom$european, matrix(1, 1, 1))
  expect_equal(default$african, matrix(1, 1, 1))
  expect_equal(default$european, matrix(1, 1, 1))
})

# ============================================================================
# ancestry_split() (Step 1 entry point)
# ============================================================================

test_that("ancestry_split validates mode argument", {
  expect_error(
    ancestry_split("x.vcf", "x.msp", mode = "bogus"),
    "should be one of"
  )
})

test_that("ancestry_split mode = 'haplotype' matches ancestry_split_phased", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # MSP: 2 samples, 1 tract at chr19:100-1000
  # S1: hap0=AFR(0), hap1=AFR(0); S2: hap0=EUR(1), hap1=EUR(1)
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0\t1\t1"
  ), con)
  close(con)

  # VCF: 1 variant at pos 150, S1=0|1, S2=1|0
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|0"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "haplotype", chrom = "chr19",
                            verbose = FALSE)

  expect_true(all(c("AFR", "EUR") %in% names(result)))
  expect_equal(result$AFR[1, 1], 1, tolerance = 1e-9)
  expect_equal(result$EUR[1, 1], 0, tolerance = 1e-9)
  expect_equal(result$AFR[1, 2], 0, tolerance = 1e-9)
  expect_equal(result$EUR[1, 2], 1, tolerance = 1e-9)
  expect_equal(result$mode, "haplotype")
  expect_equal(nrow(result$variant_info), 1)
  expect_equal(result$sample_ids, c("S1", "S2"))
})

test_that("ancestry_split mode = 'dosage' splits mixed hets proportionally", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # S1 is mixed ancestry (hap0=AFR, hap1=EUR); S2 pure EUR
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1"
  ), con)
  writeLines(c(
    "chr19\t100\t500\t0.1\t0.2\t3\t0\t1\t1\t1"
  ), con)
  close(con)

  # S1 is the only carrier of the alt allele (singleton) -> p1 = p2 = 0.5
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr19\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0|1\t0|0"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr19",
                            verbose = FALSE)

  expect_equal(result$mode, "dosage")
  # Singleton mixed het: proportional split defaults to 0.5 / 0.5
  expect_equal(result$AFR[1, 1] + result$EUR[1, 1], 1, tolerance = 1e-6)
  expect_equal(result$AFR[1, 1], 0.5, tolerance = 1e-6)
  expect_equal(result$EUR[1, 1], 0.5, tolerance = 1e-6)
})

test_that("ancestry_split returns ancestry_counts aligned with variant_info", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # 3 samples: S1, S2 pure AFR; S3 pure EUR
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1\tS3.0\tS3.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0\t0\t0\t1\t1"
  ), con)
  close(con)

  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr19",
                            verbose = FALSE)

  expect_true("ancestry_counts" %in% names(result))
  expect_equal(dim(result$ancestry_counts), c(1, 2))
  expect_equal(colnames(result$ancestry_counts), c("AFR", "EUR"))
  # 2 pure-AFR samples (S1, S2), 1 pure-EUR sample (S3)
  expect_equal(unname(result$ancestry_counts[1, "AFR"]), 2)
  expect_equal(unname(result$ancestry_counts[1, "EUR"]), 1)
})

test_that("ancestry_split supports K = 3 populations", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1\tNAT=2",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1\tS3.0\tS3.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0\t1\t1\t2\t2"
  ), con)
  close(con)

  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "haplotype", chrom = "chr19",
                            verbose = FALSE)

  expect_true(all(c("AFR", "EUR", "NAT") %in% names(result)))
  expect_equal(ncol(result$ancestry_counts), 3)
  expect_equal(result$AFR[1, 1], 1, tolerance = 1e-9)  # S1 hap1 is AFR w/ alt
  expect_equal(result$EUR[1, 2], 2, tolerance = 1e-9)  # S2 pure EUR hom alt
  expect_equal(result$NAT[1, 3], 0, tolerance = 1e-9)  # S3 pure NAT, no alt
})

# ============================================================================
# ancestry_split_dosage() (backward-compatible matrix-input Step 1 wrapper)
# ============================================================================

test_that("ancestry_split_dosage handles sample overlap", {
  # GT: 3 samples (A, B, C)
  # PT: 3 samples (A, B, D)
  # Expected: only A and B are common
  gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1), nrow = 3, ncol = 3,
               dimnames = list(c("var1", "var2", "var3"),
                               c("sample_A", "sample_B", "sample_C")))

  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
               dimnames = list(c("sample_A", "sample_B", "sample_D"),
                               c("var1", "var2")))

  result <- ancestry_split_dosage(gt, pt, verbose = FALSE)

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

test_that("ancestry_split_dosage handles coordinate-based variant matching", {
  # GT variants: chr22:100, chr22:200, chr22:300
  # PT regions: chr22:50-150, chr22:150-250
  # Expected: chr22:100 matches chr22:50-150, chr22:200 matches chr22:150-250

  gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1), nrow = 3, ncol = 3,
               dimnames = list(c("chr22:100", "chr22:200", "chr22:300"),
                               c("s1", "s2", "s3")))

  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
               dimnames = list(c("s1", "s2", "s3"),
                               c("chr22:50-150", "chr22:150-250")))

  result <- ancestry_split_dosage(gt, pt, verbose = FALSE)

  # Only 2 variants should match
  expect_equal(result$overlap$n_variants_kept, 2)
  expect_equal(nrow(result$african), 2)
  expect_equal(ncol(result$african), 3)
})

test_that("ancestry_split_dosage returns correct counts", {
  gt <- matrix(c(2, 1, 0, 1), nrow = 2, ncol = 2,
               dimnames = list(c("var1", "var2"),
                               c("s1", "s2")))

  pt <- matrix(c(3, 1, 2, 2), nrow = 2, ncol = 2,
               dimnames = list(c("s1", "s2"),
                               c("var1", "var2")))

  result <- ancestry_split_dosage(gt, pt, verbose = FALSE)

  expect_true("counts" %in% names(result))
  expect_true("african" %in% names(result$counts))
  expect_true("european" %in% names(result$counts))
  expect_true("mixed" %in% names(result$counts))
})

test_that("ancestry_split_dosage handles unnamed matrices", {
  # Unnamed matrices - positional matching by order
  gt <- matrix(c(2, 1, 0, 1, 2, 1), nrow = 2, ncol = 3)
  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2)

  result <- ancestry_split_dosage(gt, pt, verbose = FALSE)

  expect_equal(result$overlap$n_samples_kept, 3)
  expect_equal(result$overlap$n_variants_kept, 2)
})

test_that("ancestry_split_dosage stops on no common samples", {
  gt <- matrix(0, nrow = 2, ncol = 2,
               dimnames = list(c("v1", "v2"), c("s1", "s2")))

  pt <- matrix(1, nrow = 2, ncol = 2,
               dimnames = list(c("s3", "s4"), c("v1", "v2")))

  expect_error(ancestry_split_dosage(gt, pt, verbose = FALSE),
               "No common samples found")
})

test_that("result contains all expected elements", {
  gt <- matrix(c(2, 1, 0, 1), nrow = 2, ncol = 2,
               dimnames = list(c("v1", "v2"), c("s1", "s2")))

  pt <- matrix(c(3, 1, 2, 2), nrow = 2, ncol = 2,
               dimnames = list(c("s1", "s2"), c("v1", "v2")))

  result <- ancestry_split_dosage(gt, pt, verbose = FALSE)

  expect_named(result, c("african", "european", "counts", "overlap"))
  expect_equal(dim(result$african), c(2, 2))
  expect_equal(dim(result$european), c(2, 2))
})

test_that("ancestry_split_dosage filters monomorphic variants", {
  # GT: var1 has alt alleles, var2 is monomorphic (all 0)
  gt <- matrix(c(2, 0, 0, 0, 0, 0), nrow = 2, ncol = 3,
               dimnames = list(c("var1", "var2"),
                               c("s1", "s2", "s3")))

  pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
               dimnames = list(c("s1", "s2", "s3"),
                               c("var1", "var2")))

  result <- ancestry_split_dosage(gt, pt, verbose = FALSE)

  # var2 should be filtered out
  expect_equal(result$overlap$n_monomorphic_filtered, 1)
  expect_equal(result$overlap$n_variants_kept, 1)
  expect_equal(nrow(result$african), 1)
})

test_that("ancestry_split_dosage errors if all variants are monomorphic", {
  gt <- matrix(0, nrow = 3, ncol = 2,
               dimnames = list(c("v1", "v2", "v3"),
                               c("s1", "s2")))

  pt <- matrix(1, nrow = 2, ncol = 3,
               dimnames = list(c("s1", "s2"),
                               c("v1", "v2", "v3")))

  expect_error(ancestry_split_dosage(gt, pt, verbose = FALSE),
               "All variants are monomorphic")
})

# ============================================================================
# ancestry_split_phased() (backward-compatible phased Step 1 + Step 2 wrapper)
# ============================================================================

test_that("ancestry_split_phased works with synthetic data", {
  skip_if_not(file.exists(file.path(.libPaths()[1], "lantern/libs/lantern.so")),
              "lantern package not installed")

  # Create temp directory
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # Synthetic MSP: 2 samples, 3 tracts on chr19
  # Line 1: pop codes, Line 2: column headers
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0\t0\t1",
    "chr19\t1000\t5000\t0.2\t0.3\t10\t1\t0\t0\t1",
    "chr19\t5000\t10000\t0.3\t0.4\t8\t0\t1\t1\t0"
  ), con)
  close(con)

  # Synthetic VCF: 3 variants, 2 samples, phased GT
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|0",
    "chr19\t2000\t.\tG\tC\t.\tPASS\t.\tGT\t1|1\t0|1",
    "chr19\t6000\t.\tT\tG\t.\tPASS\t.\tGT\t0|0\t1|0"
  ), vcf)

  # Run pipeline
  result <- ancestry_split_phased(
    vcf_path = vcf,
    msp_path = msp,
    out_path = td,
    chrom = "chr19",
    write_vcf = FALSE,
    verbose = FALSE
  )

  # Verify structure
  expect_type(result, "list")
  expect_true("african" %in% names(result))
  expect_true("european" %in% names(result))
  expect_true("variant_info" %in% names(result))
  expect_true("sample_ids" %in% names(result))
  expect_true("tract_info" %in% names(result))
  expect_true("overlap" %in% names(result))

  # Verify dimensions
  expect_equal(dim(result$african), c(3, 2))
  expect_equal(dim(result$european), c(3, 2))
  expect_equal(length(result$sample_ids), 2)
  expect_equal(result$sample_ids, c("S1", "S2"))

  # Verify overlap stats
  expect_equal(result$overlap$n_common, 2)
  expect_equal(result$overlap$n_variants_kept, 3)

  # Verify tract assignment
  expect_equal(nrow(result$tract_info), 3)

  # Verify variant_info
  expect_equal(nrow(result$variant_info), 3)
  expect_true("chrom" %in% names(result$variant_info))
  expect_true("pos" %in% names(result$variant_info))
})

test_that("ancestry_split_phased computes correct dosages", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # MSP: 2 samples, 1 tract at chr19:100-1000
  # S1: hap0=AFR(0), hap1=AFR(0)
  # S2: hap0=EUR(1), hap1=EUR(1)
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0\t1\t1"
  ), con)
  close(con)

  # VCF: 1 variant at pos 150, S1=0|1, S2=1|0
  # S1: hap0=0 (AFR), hap1=1 (AFR) -> african=0+1=1, european=0
  # S2: hap0=1 (EUR), hap1=0 (EUR) -> african=0, european=0+1=1
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|0"
  ), vcf)

  result <- ancestry_split_phased(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  # Check dosage values
  expect_equal(result$african[1, 1], 1, tolerance = 1e-9)
  expect_equal(result$european[1, 1], 0, tolerance = 1e-9)
  expect_equal(result$african[1, 2], 0, tolerance = 1e-9)
  expect_equal(result$european[1, 2], 1, tolerance = 1e-9)
})

test_that("ancestry_split_phased handles mixed ancestry", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # MSP: S1 has mixed ancestry (hap0=AFR, hap1=EUR), S2 pure EUR
  # Tract at chr19:100-500
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1"
  ), con)
  writeLines(c(
    "chr19\t100\t500\t0.1\t0.2\t3\t0\t1\t1\t1"
  ), con)
  close(con)

  # VCF: 1 variant, S1=1|1 (hom alt), S2=0|0 (hom ref)
  # S1: hap0=1 (AFR, code 0), hap1=1 (EUR, code 1) -> african=1, european=1
  # S2: hap0=0 (EUR), hap1=0 (EUR) -> african=0, european=0
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr19\t200\t.\tG\tC\t.\tPASS\t.\tGT\t1|1\t0|0"
  ), vcf)

  result <- ancestry_split_phased(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  expect_equal(result$african[1, 1], 1, tolerance = 1e-9)  # S1: AFR hap0 has alt
  expect_equal(result$european[1, 1], 1, tolerance = 1e-9)  # S1: EUR hap1 has alt
  expect_equal(result$african[1, 2], 0, tolerance = 1e-9)  # S2: no alt
  expect_equal(result$european[1, 2], 0, tolerance = 1e-9)  # S2: no alt
})

test_that("ancestry_split_phased filters monomorphic variants", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0"
  ), con)
  close(con)

  # VCF: 2 variants, second is monomorphic (0|0)
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1",
    "chr19\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0|0"
  ), vcf)

  result <- ancestry_split_phased(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  expect_equal(result$overlap$n_monomorphic_filtered, 1)
  expect_equal(nrow(result$african), 1)
  expect_equal(result$overlap$n_variants_kept, 1)
})

test_that("ancestry_split_phased handles sample mismatch", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # MSP: samples S1, S2
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0\t1\t0"
  ), con)
  close(con)

  # VCF: samples S1, S3 (S2 not in VCF, S3 not in MSP)
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS3",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t0|0"
  ), vcf)

  result <- ancestry_split_phased(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  # Only S1 is common
  expect_equal(length(result$sample_ids), 1)
  expect_equal(result$sample_ids, "S1")
  expect_true("S3" %in% result$overlap$dropped_samples_vcf)
  expect_true("S2" %in% result$overlap$dropped_samples_msp)
})

test_that("ancestry_split_phased errors on no common samples", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tX1.0\tX1.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0"
  ), con)
  close(con)

  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tY1\tY2",
    "chr19\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t0|0"
  ), vcf)

  expect_error(
    ancestry_split_phased(vcf, msp, out_path = td, chrom = "chr19",
                        write_vcf = FALSE, verbose = FALSE),
    "No common samples"
  )
})

test_that("ancestry_split_phased handles variants outside tracts", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # MSP: tract 100-1000
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1"
  ), con)
  writeLines(c(
    "chr19\t100\t1000\t0.1\t0.2\t5\t0\t0"
  ), con)
  close(con)

  # VCF: 2 variants, one inside tract (pos=500), one outside (pos=2000)
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr19\t500\t.\tA\tT\t.\tPASS\t.\tGT\t0|1",
    "chr19\t2000\t.\tG\tC\t.\tPASS\t.\tGT\t1|1"
  ), vcf)

  result <- ancestry_split_phased(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  # Only variant at pos=500 should be kept
  expect_equal(result$overlap$n_no_tract, 1)
  expect_equal(nrow(result$african), 1)
  expect_equal(result$variant_info$pos, 500)
})

# ============================================================================
# read_bed_file: PLINK .bed reader for BED-encoded local ancestry
# ============================================================================

test_that("read_bed_file matches snpStats::read.plink()'s raw numeric convention", {
  td <- tempdir()

  # 1 variant, 4 samples. True PLINK 2-bit codes (low bits = first sample
  # in the byte): 00 = hom. first allele, 01 = missing, 10 = het,
  # 11 = hom. second allele.
  #   S0 = 00, S1 = 01, S2 = 10, S3 = 11
  byte_val <- 0L + bitwShiftL(1L, 2) + bitwShiftL(2L, 4) + bitwShiftL(3L, 6)

  bed <- file.path(td, "rbf_test.bed")
  bim <- file.path(td, "rbf_test.bim")
  fam <- file.path(td, "rbf_test.fam")

  con <- file(bed, "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  writeBin(as.raw(byte_val), con)
  close(con)

  writeLines("22\trs1\t0\t1000\tA\tG", bim)
  writeLines(c("FAM\tS0\t0\t0\t0\t-9", "FAM\tS1\t0\t0\t0\t-9",
               "FAM\tS2\t0\t0\t0\t-9", "FAM\tS3\t0\t0\t0\t-9"), fam)

  mat <- read_bed_file(bed, bim, fam)

  expect_equal(dim(mat), c(1L, 4L))
  expect_equal(rownames(mat), "rs1")
  expect_equal(colnames(mat), c("S0", "S1", "S2", "S3"))
  # snpStats convention: 0 = missing, 1 = hom. first allele, 2 = het,
  # 3 = hom. second allele -- i.e. NOT the raw on-disk bit order (0,1,2,3).
  expect_equal(as.integer(mat["rs1", ]), c(1L, 0L, 2L, 3L))
})

# ============================================================================
# ancestry_split: chrom filtering with a multi-chromosome MSP/VCF
# ============================================================================

test_that("ancestry_split(chrom=) does not confuse tracts across chromosomes with overlapping positions", {
  td <- tempdir()

  # chr1 tract 100-200: S0 = pure EUR. chr2 tract 100-200 (same range!):
  # S0 = pure AFR. Without per-chromosome tract filtering, findInterval()
  # on a chromosome-agnostic tract_df could match a chr1 variant to chr2's
  # tract purely by position overlap.
  msp <- tempfile(fileext = ".tsv", tmpdir = td)
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS0.0\tS0.1",
    "chr1\t100\t200\t0.0\t0.1\t5\t1\t1",
    "chr2\t100\t200\t0.0\t0.1\t5\t0\t0"
  ), msp)

  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=248956422>",
    "##contig=<ID=chr2,length=242193529>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0",
    "chr1\t150\t.\tA\tT\t.\tPASS\t.\tGT\t0|1",
    "chr2\t150\t.\tC\tG\t.\tPASS\t.\tGT\t0|1"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr1", verbose = FALSE)

  expect_equal(nrow(result$variant_info), 1)
  expect_equal(result$variant_info$pos, 150)
  # S0 is pure EUR on chr1's tract -- must NOT pick up chr2's pure-AFR tract.
  expect_equal(result$EUR[1, 1], 1)
  expect_equal(result$AFR[1, 1], 0)
})
