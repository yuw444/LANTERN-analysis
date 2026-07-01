context("Phased pipeline")

test_that("run_phased_pipeline works with synthetic data", {
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
  result <- run_phased_pipeline(
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

test_that("run_phased_pipeline computes correct dosages", {
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

  result <- run_phased_pipeline(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  # Check dosage values
  expect_equal(result$african[1, 1], 1, tolerance = 1e-9)
  expect_equal(result$european[1, 1], 0, tolerance = 1e-9)
  expect_equal(result$african[1, 2], 0, tolerance = 1e-9)
  expect_equal(result$european[1, 2], 1, tolerance = 1e-9)
})

test_that("run_phased_pipeline handles mixed ancestry", {
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

  result <- run_phased_pipeline(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  expect_equal(result$african[1, 1], 1, tolerance = 1e-9)  # S1: AFR hap0 has alt
  expect_equal(result$european[1, 1], 1, tolerance = 1e-9)  # S1: EUR hap1 has alt
  expect_equal(result$african[1, 2], 0, tolerance = 1e-9)  # S2: no alt
  expect_equal(result$european[1, 2], 0, tolerance = 1e-9)  # S2: no alt
})

test_that("run_phased_pipeline filters monomorphic variants", {
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

  result <- run_phased_pipeline(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  expect_equal(result$overlap$n_monomorphic_filtered, 1)
  expect_equal(nrow(result$african), 1)
  expect_equal(result$overlap$n_variants_kept, 1)
})

test_that("run_phased_pipeline handles sample mismatch", {
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

  result <- run_phased_pipeline(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  # Only S1 is common
  expect_equal(length(result$sample_ids), 1)
  expect_equal(result$sample_ids, "S1")
  expect_true("S3" %in% result$overlap$dropped_samples_vcf)
  expect_true("S2" %in% result$overlap$dropped_samples_msp)
})

test_that("run_phased_pipeline errors on no common samples", {
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
    run_phased_pipeline(vcf, msp, out_path = td, chrom = "chr19",
                        write_vcf = FALSE, verbose = FALSE),
    "No common samples"
  )
})

test_that("run_phased_pipeline handles variants outside tracts", {
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

  result <- run_phased_pipeline(vcf, msp, out_path = td, chrom = "chr19",
                                 write_vcf = FALSE, verbose = FALSE)

  # Only variant at pos=500 should be kept
  expect_equal(result$overlap$n_no_tract, 1)
  expect_equal(nrow(result$african), 1)
  expect_equal(result$variant_info$pos, 500)
})
