# ============================================================================
# VCF-reading backends: bcftools (when on PATH) vs SeqArray (fallback)
# ============================================================================

# Fixture exercising every GT-string edge case the two backends must agree
# on: phased/unphased, fully and partially missing, multiallelic (dropped),
# a site with no ALT, and a second chromosome (filtered out by `chrom`).
make_reader_fixture <- function(td) {
  samples <- paste0("S", 1:4)
  hdr_samples <- paste(vapply(samples, function(s) paste0(s, c(".0", ".1")),
                              character(2)), collapse = "\t")
  msp <- file.path(td, "toy.msp.tsv")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    paste0("#chm\tspos\tepos\tsgpos\tegpos\tn snps\t", hdr_samples),
    "chr19\t100\t1000\t0.1\t0.2\t10\t0\t0\t1\t1\t0\t1\t1\t0",
    "chr19\t1001\t5000\t0.2\t0.5\t10\t0\t1\t1\t1\t0\t0\t1\t0",
    "chr20\t100\t5000\t0.1\t0.5\t10\t0\t0\t0\t0\t0\t0\t0\t0"
  ), msp)

  site <- function(chr, pos, ref, alt, gts)
    paste(c(chr, pos, ".", ref, alt, ".", "PASS", ".", "GT", gts),
          collapse = "\t")
  vcf <- file.path(td, "toy.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr19>",
    "##contig=<ID=chr20>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", samples), collapse = "\t"),
    site("chr19", 150,  "A", "G",   c("0|1", "1|1", "0|0", "1|0")),
    site("chr19", 300,  "C", "T",   c("0/1", "0|0", "./.", "0|.")),
    site("chr19", 600,  "G", "A,T", c("0|2", "0|1", "0|0", "0|0")),
    site("chr19", 900,  "T", ".",   c("0|0", "0|0", "0|0", "0|0")),
    site("chr19", 2000, "T", "C",   c("0|0", "0|1", "1|1", "0|0")),
    site("chr19", 4000, "A", "C",   c("1|0", "0|0", "0|0", "0|0")),
    site("chr20", 500,  "G", "C",   c("0|1", "0|0", "0|0", "0|0"))
  ), vcf)
  list(vcf = vcf, msp = msp)
}

split_with_reader <- function(reader, fx, ...) {
  old <- options(lantern.vcf_reader = reader)
  on.exit(options(old), add = TRUE)
  suppressWarnings(ancestry_split(fx$vcf, fx$msp, verbose = FALSE, ...))
}

test_that(".vcf_reader() honours the lantern.vcf_reader option", {
  old <- options(lantern.vcf_reader = "seqarray")
  on.exit(options(old), add = TRUE)
  expect_equal(lantern:::.vcf_reader(), "seqarray")
  options(lantern.vcf_reader = "bcftools")
  expect_equal(lantern:::.vcf_reader(), "bcftools")
  options(lantern.vcf_reader = "auto")
  expect_equal(lantern:::.vcf_reader(),
               if (nzchar(Sys.which("bcftools"))) "bcftools" else "seqarray")
  options(lantern.vcf_reader = "bogus")
  expect_error(lantern:::.vcf_reader())
})

test_that("SeqArray reader produces bcftools-format GT strings", {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  fx <- make_reader_fixture(td)

  res <- lantern:::.vcf_read_gt(fx$vcf, chrom = "chr19", reader = "seqarray",
                                verbose = FALSE)
  expect_equal(res$chrom, rep("chr19", 6))
  expect_equal(res$pos, c(150L, 300L, 600L, 900L, 2000L, 4000L))
  expect_equal(res$alt, c("G", "T", "A,T", ".", "C", "C"))
  expect_equal(unname(res$gt[2, ]), c("0/1", "0|0", "./.", "0|."))
  expect_equal(unname(res$gt[3, ]), c("0|2", "0|1", "0|0", "0|0"))
  expect_equal(colnames(res$gt), paste0("S", 1:4))

  expect_error(lantern:::.vcf_read_gt(fx$vcf, chrom = "chr7",
                                      reader = "seqarray", verbose = FALSE),
               "No variants on chromosome chr7")
})

test_that("ancestry_split() runs without bcftools via SeqArray", {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  fx <- make_reader_fixture(td)

  res <- split_with_reader("seqarray", fx, mode = "haplotype", chrom = "chr19")
  # multiallelic (600) dropped; monomorphic (900) dropped by the split
  expect_equal(res$variant_info$pos, c(150, 300, 2000, 4000))
  expect_equal(res$overlap$n_multiallelic_filtered, 1)
})

test_that("bcftools and SeqArray readers give identical ancestry_split() output", {
  skip_if(!nzchar(Sys.which("bcftools")), "bcftools not on PATH")
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  fx <- make_reader_fixture(td)

  for (mode in c("dosage", "haplotype")) {
    expect_identical(
      split_with_reader("seqarray", fx, mode = mode, chrom = "chr19"),
      split_with_reader("bcftools", fx, mode = mode, chrom = "chr19"),
      info = mode)
  }
})
