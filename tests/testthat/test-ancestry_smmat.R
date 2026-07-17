test_that(".compute_gene_weights computes per-gene median ancestry counts", {
  variant_info <- data.frame(chrom = rep("chr19", 3), pos = c(100L, 200L, 300L),
                              ref = "A", alt = "T", stringsAsFactors = FALSE)
  ancestry_counts <- matrix(c(2, 5, 9, 8, 5, 1), nrow = 3, ncol = 2,
                             dimnames = list(NULL, c("AFR", "EUR")))

  gg_file <- tempfile(fileext = ".tsv")
  gg <- data.frame(
    gene = c("GENE1", "GENE1", "GENE2"),
    chr  = c("19", "19", "19"),
    pos  = c(100L, 200L, 300L),
    ref  = "A", alt = "T", weight = 1
  )
  write.table(gg, gg_file, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)

  w <- .compute_gene_weights(ancestry_counts, variant_info, gg_file,
                              pop_names = c("AFR", "EUR"))

  expect_equal(dim(w), c(2, 2))
  expect_equal(unname(w["GENE1", "AFR"]), median(c(2, 5)))
  expect_equal(unname(w["GENE1", "EUR"]), median(c(8, 5)))
  expect_equal(unname(w["GENE2", "AFR"]), 9)
  expect_equal(unname(w["GENE2", "EUR"]), 1)
})

test_that(".compute_gene_weights falls back to NULL when no population names match", {
  variant_info <- data.frame(chrom = "chr19", pos = 100L, ref = "A", alt = "T",
                              stringsAsFactors = FALSE)
  ancestry_counts <- matrix(1, nrow = 1, ncol = 2,
                             dimnames = list(NULL, c("AFR", "EUR")))
  gg_file <- tempfile(fileext = ".tsv")
  write.table(data.frame(gene = "GENE1", chr = "19", pos = 100L, ref = "A",
                          alt = "T", weight = 1),
              gg_file, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)

  expect_warning(
    w <- .compute_gene_weights(ancestry_counts, variant_info, gg_file,
                                pop_names = c("FOO", "BAR")),
    "falling back to equal weights"
  )
  expect_null(w)
})

test_that("ancestry_smmat() runs SMMAT per population and Cauchy-combines p-values", {
  skip_if_not_installed("GMMAT")
  skip_if_not_installed("SeqArray")

  set.seed(42)
  n     <- 30
  ids   <- paste0("S", seq_len(n))
  n_var <- 4

  afr_mat <- matrix(rbinom(n_var * n, 2, 0.2), nrow = n_var, ncol = n,
                     dimnames = list(NULL, ids))
  eur_mat <- matrix(rbinom(n_var * n, 2, 0.2), nrow = n_var, ncol = n,
                     dimnames = list(NULL, ids))

  variant_info <- data.frame(chrom = rep("chr19", n_var),
                              pos = seq(100L, by = 100L, length.out = n_var),
                              ref = "A", alt = "T", stringsAsFactors = FALSE)
  ancestry_counts <- matrix(sample(5:15, n_var * 2, replace = TRUE),
                             nrow = n_var, ncol = 2,
                             dimnames = list(NULL, c("AFR", "EUR")))

  split_result <- list(AFR = afr_mat, EUR = eur_mat, variant_info = variant_info,
                        sample_ids = ids, mode = "dosage",
                        ancestry_counts = ancestry_counts, overlap = list())

  out_dir <- tempfile()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  gds_paths <- write_ancestry_gds(split_result, out_dir, verbose = FALSE)

  # gene_group_file passed as a data.frame directly -- no tempfile/unlink
  # management needed by the caller.
  gg <- data.frame(gene = "GENE1", chr = "19", pos = variant_info$pos,
                    ref = "A", alt = "T", weight = 1)

  pheno <- data.frame(id = ids, y = rnorm(n), age = rnorm(n, 50, 10),
                       stringsAsFactors = FALSE)
  kinship <- diag(n)
  rownames(kinship) <- colnames(kinship) <- ids

  out <- suppressWarnings(ancestry_smmat(
    gds_paths, pheno, y ~ age, kinship, gg,
    ancestry_counts = split_result$ancestry_counts,
    variant_info    = split_result$variant_info,
    verbose = FALSE
  ))

  expect_type(out, "list")
  expect_named(out, c("results", "smmat_results"))

  res <- out$results
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(res$gene, "GENE1")
  expect_true(all(c("p_AFR", "p_EUR", "w_AFR", "w_EUR", "p_cauchy") %in% names(res)))
  expect_true(res$p_cauchy > 0 && res$p_cauchy <= 1)
  expect_true(!is.na(res$w_AFR) && !is.na(res$w_EUR))

  smmat_results <- out$smmat_results
  expect_named(smmat_results, c("AFR", "EUR"))
  expect_true("E.pval" %in% names(smmat_results$AFR))
})

test_that("ancestry_smmat() still accepts gene_group_file as a plain file path", {
  skip_if_not_installed("GMMAT")
  skip_if_not_installed("SeqArray")

  set.seed(7)
  n     <- 20
  ids   <- paste0("S", seq_len(n))
  n_var <- 3

  pop_mat <- matrix(rbinom(n_var * n, 2, 0.2), nrow = n_var, ncol = n,
                     dimnames = list(NULL, ids))
  variant_info <- data.frame(chrom = rep("chr2", n_var),
                              pos = seq(50L, by = 50L, length.out = n_var),
                              ref = "A", alt = "T", stringsAsFactors = FALSE)

  out_dir <- tempfile()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  gds_paths <- write_ancestry_gds(
    list(POP1 = pop_mat, variant_info = variant_info, sample_ids = ids),
    out_dir, verbose = FALSE)

  gg_file <- tempfile(fileext = ".tsv")
  write.table(data.frame(gene = "GENEY", chr = "2", pos = variant_info$pos,
                          ref = "A", alt = "T", weight = 1),
              gg_file, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)

  pheno <- data.frame(id = ids, y = rnorm(n), stringsAsFactors = FALSE)
  kinship <- diag(n)
  rownames(kinship) <- colnames(kinship) <- ids

  out <- suppressWarnings(ancestry_smmat(gds_paths, pheno, y ~ 1, kinship, gg_file,
                                          verbose = FALSE))

  expect_false(file.exists(gg_file))  # ancestry_smmat() unlinks it, even a caller-supplied path
  expect_equal(out$results$gene, "GENEY")
})

test_that("ancestry_smmat() falls back to equal weights when ancestry_counts is omitted", {
  skip_if_not_installed("GMMAT")
  skip_if_not_installed("SeqArray")

  set.seed(1)
  n     <- 20
  ids   <- paste0("S", seq_len(n))
  n_var <- 3

  pop_mat <- matrix(rbinom(n_var * n, 2, 0.2), nrow = n_var, ncol = n,
                     dimnames = list(NULL, ids))
  variant_info <- data.frame(chrom = rep("chr1", n_var),
                              pos = seq(10L, by = 10L, length.out = n_var),
                              ref = "A", alt = "T", stringsAsFactors = FALSE)

  out_dir <- tempfile()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  gds_paths <- write_ancestry_gds(
    list(POP1 = pop_mat, variant_info = variant_info, sample_ids = ids),
    out_dir, verbose = FALSE)

  gg <- data.frame(gene = "GENEX", chr = "1", pos = variant_info$pos,
                    ref = "A", alt = "T", weight = 1)

  pheno <- data.frame(id = ids, y = rnorm(n), stringsAsFactors = FALSE)
  kinship <- diag(n)
  rownames(kinship) <- colnames(kinship) <- ids

  out <- suppressWarnings(ancestry_smmat(gds_paths, pheno, y ~ 1, kinship, gg,
                                          verbose = FALSE))

  expect_false(any(grepl("^w_", names(out$results))))
  expect_true("p_cauchy" %in% names(out$results))
})

test_that("full Step1->2->3 pipeline works for K = 4 populations", {
  skip_if_not_installed("GMMAT")
  skip_if_not_installed("SeqArray")

  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  set.seed(4321)
  n_samp <- 20
  samp   <- paste0("S", seq_len(n_samp))

  # MSP: 4 populations, one tract spanning all variants, random hap ancestry
  pop_codes <- c(AFR = 0L, EUR = 1L, NAT = 2L, EAS = 3L)
  hap_cols  <- unlist(lapply(samp, function(s) paste0(s, c(".0", ".1"))))
  hap_anc   <- sample(pop_codes, length(hap_cols), replace = TRUE)
  msp <- tempfile(fileext = ".tsv", tmpdir = td)
  writeLines(c(
    paste(c("#Subpopulation order/codes: AFR=0", "EUR=1", "NAT=2", "EAS=3"),
          collapse = "\t"),
    paste(c("#chm","spos","epos","sgpos","egpos","n snps", hap_cols), collapse = "\t"),
    paste(c("chr1", "1000", "3000", "0.0", "0.5", "5", hap_anc), collapse = "\t")
  ), msp)

  # VCF: 5 variants, phased, mostly rare
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  gt_lines <- vapply(seq_len(5), function(i) {
    gts <- sample(c("0|0","0|1","1|0","1|1"), n_samp, replace = TRUE,
                  prob = c(0.85, 0.06, 0.06, 0.03))
    paste(c("chr1", 1000L + i * 100L, ".", "A", "T", "100", "PASS", ".", "GT", gts),
          collapse = "\t")
  }, character(1))
  writeLines(c(
    "##fileformat=VCFv4.1",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", samp),
          collapse = "\t"),
    gt_lines
  ), vcf)

  # Step 1
  split4 <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr1", verbose = FALSE)
  pop4 <- setdiff(names(split4),
                  c("variant_info", "sample_ids", "mode",
                    "ancestry_counts", "tract_info", "overlap"))
  expect_setequal(pop4, c("AFR", "EUR", "NAT", "EAS"))
  expect_equal(ncol(split4$ancestry_counts), 4)

  # Step 2
  gds_paths <- write_ancestry_gds(split4, file.path(td, "gds"), verbose = FALSE)
  expect_named(gds_paths, pop4)

  # Step 3 -- gene_group_file passed as a data.frame directly
  gg <- data.frame(gene = "GENE_K4", chr = "1", pos = split4$variant_info$pos,
                    ref = split4$variant_info$ref, alt = split4$variant_info$alt,
                    weight = 1)

  pheno   <- data.frame(id = split4$sample_ids, y = rnorm(length(split4$sample_ids)))
  kinship <- diag(length(split4$sample_ids))
  rownames(kinship) <- colnames(kinship) <- split4$sample_ids

  out <- suppressWarnings(ancestry_smmat(
    gds_paths, pheno, y ~ 1, kinship, gg,
    ancestry_counts = split4$ancestry_counts,
    variant_info    = split4$variant_info,
    verbose = FALSE
  ))

  res <- out$results
  expect_equal(nrow(res), 1)
  for (pop in pop4) {
    expect_true(paste0("p_", pop) %in% names(res))
    expect_true(paste0("w_", pop) %in% names(res))
  }
  expect_true(res$p_cauchy > 0 && res$p_cauchy <= 1)
  expect_named(out$smmat_results, pop4)
})
