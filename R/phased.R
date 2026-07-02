# ============================================================================
# Phased ancestry splitting
# ============================================================================

.split_phased <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                           pop_codes = c(AFR = 0L, EUR = 1L)) {
    storage.mode(gt_hap0)  <- "integer"
    storage.mode(gt_hap1)  <- "integer"
    storage.mode(anc_hap0) <- "integer"
    storage.mode(anc_hap1) <- "integer"
    pop_codes <- as.integer(pop_codes)
    .Call("split_phased_by_ancestry_C", gt_hap0, gt_hap1, anc_hap0, anc_hap1,
          pop_codes, PACKAGE = "lantern")
}

#' Split phased haplotypes by ancestry
#'
#' Split phased haplotypes deterministically into African and European
#' ancestry-specific dosage matrices.  Each haplotype's allele is added
#' to whichever ancestry pool matches its local ancestry call; haplotypes
#' with unrecognised ancestry codes contribute nothing to either pool.
#'
#' @param gt_hap0 Integer matrix (variants × samples) of haplotype-0 alleles (0/1).
#' @param gt_hap1 Integer matrix (variants × samples) of haplotype-1 alleles (0/1).
#' @param anc_hap0 Integer matrix (variants × samples) of haplotype-0 ancestry codes.
#' @param anc_hap1 Integer matrix (variants × samples) of haplotype-1 ancestry codes.
#' @param pop_codes Named integer vector of length 2 giving the AFR and EUR
#'   ancestry codes as used in the MSP file.  Defaults to
#'   \code{c(AFR = 0L, EUR = 1L)} (RFMix convention).
#'
#' @return List with two numeric matrices:
#'   \item{african}{African ancestry-specific dosage (variants × samples)}
#'   \item{european}{European ancestry-specific dosage (variants × samples)}
#'
#' @examples
#' gt0 <- matrix(c(1L, 0L, 0L, 1L), 2, 2)
#' gt1 <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
#' a0  <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
#' a1  <- matrix(c(1L, 0L, 1L, 0L), 2, 2)
#' split_phased(gt0, gt1, a0, a1)
#'
#' @export
split_phased <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                         pop_codes = c(AFR = 0L, EUR = 1L)) {
    .split_phased(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes)
}

#' Split phased haplotypes into K population-specific dosage matrices
#'
#' Generalises \code{\link{split_phased}} to an arbitrary number of ancestries.
#' Each haplotype's allele is routed to the dosage pool whose code matches its
#' local ancestry call.  Haplotypes with unrecognised codes contribute 0 to
#' all pools.  NA genotypes are treated as 0 (reference allele).
#'
#' For the two-population case \code{split_phased_multi} and \code{split_phased}
#' are equivalent; prefer \code{split_phased_multi} for new code that may later
#' be extended to three or more ancestries.
#'
#' @param gt_hap0 Integer matrix (variants x samples) of haplotype-0 alleles (0/1).
#' @param gt_hap1 Integer matrix (variants x samples) of haplotype-1 alleles (0/1).
#' @param anc_hap0 Integer matrix (variants x samples) of haplotype-0 ancestry codes.
#' @param anc_hap1 Integer matrix (variants x samples) of haplotype-1 ancestry codes.
#' @param pop_codes Named integer vector mapping population label to ancestry code,
#'   e.g. \code{c(AFR = 0L, EUR = 1L, NAT = 2L)}.  Codes must match those used
#'   in the MSP file (RFMix default: 0-based integers).
#'
#' @return Named list of K numeric matrices (variants x samples), one per
#'   entry in \code{pop_codes}, in the same order.  List names equal the
#'   names of \code{pop_codes}.
#'
#' @seealso \code{\link{split_phased}} for the two-population convenience wrapper,
#'   \code{\link{split_by_ancestry}} for unphased (diploid-code) splitting.
#'
#' @examples
#' # Three-population panel: AFR=0, EUR=1, NAT=2
#' gt0 <- matrix(c(1L, 0L, 0L, 1L, 1L, 0L), nrow = 3, ncol = 2)
#' gt1 <- matrix(c(0L, 1L, 1L, 0L, 0L, 1L), nrow = 3, ncol = 2)
#' a0  <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 3, ncol = 2)
#' a1  <- matrix(c(1L, 2L, 0L, 2L, 0L, 1L), nrow = 3, ncol = 2)
#' out <- split_phased_multi(gt0, gt1, a0, a1, c(AFR = 0L, EUR = 1L, NAT = 2L))
#' names(out)   # "AFR" "EUR" "NAT"
#' out$NAT      # native-ancestry dosage matrix
#'
#' @export
split_phased_multi <- function(gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                                pop_codes = c(AFR = 0L, EUR = 1L)) {
    if (is.null(names(pop_codes)))
        stop("pop_codes must be a named integer vector, e.g. c(AFR=0L, EUR=1L, NAT=2L)")
    if (length(pop_codes) < 2L)
        stop("pop_codes must contain at least 2 populations")
    pop_names <- names(pop_codes)
    storage.mode(gt_hap0)  <- "integer"
    storage.mode(gt_hap1)  <- "integer"
    storage.mode(anc_hap0) <- "integer"
    storage.mode(anc_hap1) <- "integer"
    storage.mode(pop_codes) <- "integer"   # preserves names unlike as.integer()
    result <- .Call("split_phased_multi_C", gt_hap0, gt_hap1, anc_hap0, anc_hap1,
                    pop_codes, PACKAGE = "lantern")
    names(result) <- pop_names             # guarantee names even if C skips them
    result
}

# ============================================================================
# MSP / VCF parsing helpers
# ============================================================================

#' Parse RFMix MSP file
#'
#' Parse the RFMix MSP (local ancestry tract) file into ancestry matrices.
#' MSP format: two header lines followed by data rows.
#' Line 1: population codes (e.g., "#Subpopulation order/codes: AFR=0  EUR=1").
#' Line 2: column headers.
#' Data rows: ancestry calls per haplotype per tract.
#'
#' @param msp_path Path to MSP file (plain text or gzipped).
#' @param verbose Print progress messages.
#' @return List with: \code{pop_codes} (named integer vector),
#'   \code{sample_ids} (character), \code{tract_df} (data.frame sorted by
#'   spos), \code{anc_hap0} (integer matrix: samples × tracts),
#'   \code{anc_hap1} (integer matrix: samples × tracts).
#' @keywords internal
.parse_msp <- function(msp_path, verbose = TRUE) {
  if (verbose) message("  Parsing MSP file: ", msp_path)

  con <- if (grepl("\\.gz$", msp_path)) gzfile(msp_path, "rt") else file(msp_path, "r")
  on.exit(close(con))

  line1 <- readLines(con, n = 1)
  line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
  code_parts <- strsplit(line1, if (grepl("\t", line1)) "\t" else "\\s+")[[1]]
  pop_codes <- integer()
  for (part in code_parts) {
    m <- regmatches(part, regexec("([A-Za-z0-9_]+)\\s*=\\s*([0-9]+)", part))[[1]]
    if (length(m) == 3) pop_codes[m[2]] <- as.integer(m[3])
  }
  if (length(pop_codes) == 0) {
    warning("Could not parse pop codes from MSP; using default AFR=0, EUR=1")
    pop_codes <- c(AFR = 0L, EUR = 1L)
  }

  line2      <- readLines(con, n = 1)
  headers    <- strsplit(line2, "\t")[[1]]
  sample_hap_cols <- headers[7:length(headers)]
  sample_ids <- unique(sub("\\.[01]$", "", sample_hap_cols))

  if (verbose) {
    message("    Found ", length(sample_ids), " samples, ", length(pop_codes), " populations")
    message("    Pop codes: ", paste(names(pop_codes), "=", pop_codes, collapse = ", "))
  }

  data_lines <- readLines(con)
  if (length(data_lines) == 0) stop("MSP file has no data rows")

  data_mat <- do.call(rbind, lapply(data_lines, function(l) strsplit(l, "\t")[[1]]))

  tract_df <- data.frame(
    chrom = data_mat[, 1],
    spos  = as.integer(data_mat[, 2]),
    epos  = as.integer(data_mat[, 3]),
    stringsAsFactors = FALSE
  )
  tract_order <- order(tract_df$spos)
  tract_df    <- tract_df[tract_order, ]
  rownames(tract_df) <- NULL

  # Ancestry columns start at column 7 in data rows (6 metadata fields: chm spos epos sgpos egpos "n snps")
  anc_values   <- data_mat[, 7:ncol(data_mat), drop = FALSE]
  n_tracts_orig <- nrow(tract_df)
  n_samples     <- length(sample_ids)

  anc_hap0_mat <- matrix(0L, nrow = n_samples, ncol = n_tracts_orig)
  anc_hap1_mat <- matrix(0L, nrow = n_samples, ncol = n_tracts_orig)
  rownames(anc_hap0_mat) <- rownames(anc_hap1_mat) <- sample_ids

  for (j in seq_len(n_tracts_orig)) {
    vals <- anc_values[tract_order[j], ]
    for (s in seq_len(n_samples)) {
      anc_hap0_mat[s, j] <- as.integer(vals[2 * s - 1])
      anc_hap1_mat[s, j] <- as.integer(vals[2 * s])
    }
  }

  list(pop_codes  = pop_codes,
       sample_ids = sample_ids,
       tract_df   = tract_df,
       anc_hap0   = anc_hap0_mat,
       anc_hap1   = anc_hap1_mat)
}

#' Parse phased GT matrix into haplotype integer matrices (vectorised)
#'
#' @param gt_mat Character matrix (variants × samples).
#' @param n_variants Number of variants.
#' @param n_samples Number of samples.
#' @param sample_order Column names for output matrices.
#' @param verbose Print progress.
#' @return List: \code{hap0}, \code{hap1} (integer matrices), \code{missing}
#'   (logical matrix).
#' @keywords internal
.parse_phased_gt_matrix <- function(gt_mat, n_variants, n_samples,
                                    sample_order, verbose = FALSE) {
  gt_vec  <- as.character(gt_mat)
  missing_vec <- is.na(gt_vec) | gt_vec %in% c(".", "./.", ".|.")

  unphased <- grepl("/", gt_vec, fixed = TRUE) & !grepl("|", gt_vec, fixed = TRUE)
  if (any(unphased))
    warning("Unphased GT (\"/\") encountered for ", sum(unphased),
            " cells. Treating first allele as hap0.")

  gt_norm <- gt_vec
  gt_norm[unphased] <- sub("/", "|", gt_norm[unphased], fixed = TRUE)

  nch      <- nchar(gt_norm)
  hap0_str <- ifelse(nch >= 1, substr(gt_norm, 1, 1), "0")
  hap1_str <- ifelse(nch >= 3, substr(gt_norm, 3, 3),
                     ifelse(nch >= 1, substr(gt_norm, 1, 1), "0"))
  hap0_str[missing_vec] <- "0"
  hap1_str[missing_vec] <- "0"

  hap0_vec <- suppressWarnings(as.integer(hap0_str))
  hap1_vec <- suppressWarnings(as.integer(hap1_str))
  hap0_vec[is.na(hap0_vec)] <- 0L
  hap1_vec[is.na(hap1_vec)] <- 0L

  multi <- hap0_vec > 1L | hap1_vec > 1L
  if (any(multi)) {
    warning("GT with allele >1 (multiallelic?) for ", sum(multi),
            " cells. Treating as missing (0/0).")
    hap0_vec[multi] <- 0L; hap1_vec[multi] <- 0L; missing_vec[multi] <- TRUE
  }

  hap0_mat <- matrix(hap0_vec, nrow = n_variants, ncol = n_samples)
  hap1_mat <- matrix(hap1_vec, nrow = n_variants, ncol = n_samples)
  missing_mat <- matrix(missing_vec, nrow = n_variants, ncol = n_samples)
  colnames(hap0_mat) <- colnames(hap1_mat) <- sample_order

  list(hap0 = hap0_mat, hap1 = hap1_mat, missing = missing_mat)
}

# ============================================================================
# run_phased_pipeline
# ============================================================================

#' Run phased ancestry splitting pipeline
#'
#' Complete pipeline for phased data: parse an RFMix MSP file and a phased
#' VCF, intersect samples, map each variant to its ancestry tract, call
#' \code{\link{split_phased}} via the C backend, and optionally write
#' ancestry-specific VCFs with a DS field that can be converted to GDS.
#'
#' @param vcf_path Path to phased VCF/BCF file (plain or gzipped).
#'   \code{bcftools} must be in \code{PATH}.
#' @param msp_path Path to RFMix MSP file (plain text or gzipped TSV).
#' @param out_path Output directory for VCFs, GDS, and cache files.
#' @param chrom Chromosome to process (e.g., \code{"chr19"} or \code{"19"}).
#'   If \code{NULL}, all chromosomes present in the VCF are used.
#' @param write_vcf Logical; if \code{TRUE}, write ancestry-specific VCFs
#'   and convert them to GDS format using \code{SeqArray::seqVCF2GDS()}.
#' @param verbose Print step-by-step progress messages.
#'
#' @return Invisibly, a list with elements:
#'   \item{african}{Numeric matrix (variants × samples) of African dosages.}
#'   \item{european}{Numeric matrix (variants × samples) of European dosages.}
#'   \item{variant_info}{data.frame with columns chrom, pos, ref, alt.}
#'   \item{sample_ids}{Character vector of common sample IDs.}
#'   \item{tract_info}{data.frame of ancestry tracts from the MSP file.}
#'   \item{overlap}{List of intersection statistics.}
#'   \item{vcf_paths}{(when \code{write_vcf = TRUE}) Named list of output
#'     file paths.}
#'
#' @examples
#' \dontrun{
#' result <- run_phased_pipeline(
#'   vcf_path = "data/chr19.phased.bcf",
#'   msp_path = "data/chr19.msp.tsv.gz",
#'   out_path = "output/",
#'   chrom    = "chr19",
#'   write_vcf = TRUE
#' )
#' head(result$variant_info)
#' }
#'
#' @export
run_phased_pipeline <- function(vcf_path, msp_path, out_path,
                                chrom = NULL, write_vcf = TRUE,
                                verbose = TRUE) {

  if (verbose) message("=== LANTERN Phased Ancestry Pipeline ===\n")

  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  cache_dir <- file.path(out_path, "cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  # ---- Step 1: Parse MSP ----
  if (verbose) message("Step 1: Parsing RFMix MSP file...")
  msp_data      <- .parse_msp(msp_path, verbose = verbose)
  msp_samples   <- msp_data$sample_ids
  tract_df      <- msp_data$tract_df
  anc_hap0_mat  <- msp_data$anc_hap0
  anc_hap1_mat  <- msp_data$anc_hap1
  pop_codes     <- msp_data$pop_codes
  n_tracts      <- nrow(tract_df)
  if (verbose) message("  MSP samples: ", length(msp_samples),
                       "  Tracts: ", n_tracts)

  # ---- Step 2: VCF sample IDs ----
  if (verbose) message("\nStep 2: Querying VCF sample IDs...")
  vcf_samples <- tryCatch({
    res <- system(paste0("bcftools query -l ", shQuote(vcf_path)),
                  intern = TRUE, ignore.stderr = TRUE)
    if (length(res) == 0) stop("bcftools query -l returned nothing")
    res
  }, error = function(e) stop("Could not read VCF sample IDs: ", e$message))
  if (verbose) message("  VCF samples: ", length(vcf_samples))

  # ---- Step 3: Intersect samples ----
  if (verbose) message("\nStep 3: Intersecting samples...")
  common_samples <- intersect(vcf_samples, msp_samples)
  if (length(common_samples) == 0)
    stop("No common samples between VCF and MSP")
  dropped_vcf <- setdiff(vcf_samples, common_samples)
  dropped_msp <- setdiff(msp_samples, common_samples)
  if (verbose) message("  Common: ", length(common_samples))

  anc_hap0_common <- anc_hap0_mat[common_samples, , drop = FALSE]
  anc_hap1_common <- anc_hap1_mat[common_samples, , drop = FALSE]

  # ---- Step 4: Parse VCF GT ----
  if (verbose) message("\nStep 4: Parsing VCF genotypes...")
  cmd_gt <- paste0(
    "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ",
    shQuote(vcf_path))
  gt_output <- tryCatch(
    system(cmd_gt, intern = TRUE, ignore.stderr = TRUE),
    error = function(e) stop("bcftools query failed: ", e$message))
  if (length(gt_output) == 0) stop("bcftools query returned no variants")
  if (verbose) message("  ", length(gt_output), " variant lines read")

  tmp_gt <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_gt), add = TRUE)
  writeLines(gt_output, tmp_gt)
  gt_dt <- data.table::fread(tmp_gt, header = FALSE, sep = "\t")

  vcf_chrom <- gt_dt[[1]]
  vcf_pos   <- as.integer(gt_dt[[2]])
  vcf_ref   <- as.character(gt_dt[[3]])
  vcf_alt   <- as.character(gt_dt[[4]])

  is_biallelic <- !grepl(",", vcf_alt, fixed = TRUE)
  n_multi <- sum(!is_biallelic)
  if (n_multi > 0) {
    if (verbose) message("  Skipping ", n_multi, " multiallelic variants")
    vcf_chrom <- vcf_chrom[is_biallelic]; vcf_pos <- vcf_pos[is_biallelic]
    vcf_ref   <- vcf_ref[is_biallelic];   vcf_alt <- vcf_alt[is_biallelic]
    gt_dt     <- gt_dt[is_biallelic]
  }
  n_variants <- length(vcf_pos)

  sample_idx <- match(common_samples, vcf_samples)
  gt_mat_raw <- as.matrix(gt_dt[, 5:ncol(gt_dt), drop = FALSE])
  gt_mat     <- gt_mat_raw[, sample_idx, drop = FALSE]

  gt_parsed  <- .parse_phased_gt_matrix(gt_mat, n_variants,
                                         length(common_samples),
                                         common_samples, verbose = verbose)
  gt_hap0_mat <- gt_parsed$hap0
  gt_hap1_mat <- gt_parsed$hap1

  n_missing <- sum(gt_parsed$missing)
  if (n_missing > 0 && verbose)
    message("  Missing GT: ", n_missing, " (treated as 0/0)")

  # ---- Step 5: Map variants to tracts ----
  if (verbose) message("\nStep 5: Mapping variants to ancestry tracts...")
  vcf_chrom_clean   <- sub("^chr", "", vcf_chrom)
  tract_chrom_clean <- sub("^chr", "", tract_df$chrom)

  if (!is.null(chrom)) {
    chrom_clean <- sub("^chr", "", chrom)
    keep <- vcf_chrom_clean == chrom_clean
    if (!any(keep)) stop("No variants on chromosome ", chrom)
    vcf_chrom <- vcf_chrom[keep]; vcf_pos <- vcf_pos[keep]
    vcf_ref   <- vcf_ref[keep];   vcf_alt <- vcf_alt[keep]
    gt_hap0_mat <- gt_hap0_mat[keep, , drop = FALSE]
    gt_hap1_mat <- gt_hap1_mat[keep, , drop = FALSE]
    vcf_chrom_clean <- vcf_chrom_clean[keep]
    n_variants <- sum(keep)
    if (verbose) message("  Filtered to chrom ", chrom, ": ", n_variants, " variants")
  }

  tract_idx   <- findInterval(vcf_pos, tract_df$spos)
  valid_tract <- tract_idx > 0
  if (any(valid_tract))
    valid_tract[valid_tract] <-
      vcf_pos[valid_tract] <= tract_df$epos[tract_idx[valid_tract]]
  valid_chrom <- if (!is.null(chrom)) {
    rep(TRUE, length(tract_idx))
  } else {
    vapply(seq_along(tract_idx), function(i) {
      tract_idx[i] > 0 &&
        vcf_chrom_clean[i] == tract_chrom_clean[tract_idx[i]]
    }, logical(1))
  }
  keep_var <- valid_tract & valid_chrom
  n_no_tract <- sum(!keep_var)
  if (n_no_tract > 0 && verbose)
    message("  Variants with no tract: ", n_no_tract, " (dropped)")

  if (n_no_tract > 0) {
    vcf_chrom   <- vcf_chrom[keep_var]; vcf_pos <- vcf_pos[keep_var]
    vcf_ref     <- vcf_ref[keep_var];   vcf_alt <- vcf_alt[keep_var]
    gt_hap0_mat <- gt_hap0_mat[keep_var, , drop = FALSE]
    gt_hap1_mat <- gt_hap1_mat[keep_var, , drop = FALSE]
    tract_idx   <- tract_idx[keep_var]
    n_variants  <- sum(keep_var)
  }
  if (n_variants == 0) stop("No variants with valid ancestry tracts")

  # ---- Step 6: Broadcast tract ancestry to variants ----
  if (verbose) message("\nStep 6: Broadcasting ancestry to variants...")
  anc_hap0_var <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  anc_hap1_var <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  for (i in seq_len(n_variants)) {
    t <- tract_idx[i]
    anc_hap0_var[i, ] <- anc_hap0_common[, t]
    anc_hap1_var[i, ] <- anc_hap1_common[, t]
  }

  # ---- Step 7: Split phased ----
  if (verbose) message("\nStep 7: Splitting haplotypes by ancestry...")
  res <- split_phased(gt_hap0_mat, gt_hap1_mat,
                      anc_hap0_var, anc_hap1_var,
                      pop_codes = pop_codes)
  african_mat  <- res$african
  european_mat <- res$european
  if (verbose) message("  African: ", nrow(african_mat), " x ", ncol(african_mat))

  # ---- Step 8: Filter monomorphic ----
  if (verbose) message("\nStep 8: Filtering monomorphic variants...")
  has_alt       <- rowSums(african_mat + european_mat) > 0
  n_monomorphic <- sum(!has_alt)
  if (n_monomorphic > 0) {
    if (verbose) message("  Removed ", n_monomorphic, " monomorphic variants")
    african_mat  <- african_mat[has_alt, , drop = FALSE]
    european_mat <- european_mat[has_alt, , drop = FALSE]
    vcf_chrom <- vcf_chrom[has_alt]; vcf_pos <- vcf_pos[has_alt]
    vcf_ref   <- vcf_ref[has_alt];   vcf_alt <- vcf_alt[has_alt]
    tract_idx <- tract_idx[has_alt]
  }
  n_final <- length(vcf_pos)
  if (n_final == 0) stop("All variants are monomorphic after filtering")

  variant_info <- data.frame(chrom = vcf_chrom, pos = vcf_pos,
                              ref = vcf_ref, alt = vcf_alt,
                              stringsAsFactors = FALSE)

  # ---- Step 9: Write VCFs / GDS ----
  vcf_paths <- NULL
  if (write_vcf) {
    if (verbose) message("\nStep 9: Writing ancestry-specific VCFs...")

    sample_file <- file.path(cache_dir, "common_samples.txt")
    writeLines(common_samples, sample_file)

    subset_vcf <- file.path(cache_dir, "subset.vcf")
    status <- system2("bcftools",
                      c("view", "-S", shQuote(sample_file), shQuote(vcf_path)),
                      stdout = subset_vcf)
    if (status != 0) stop("bcftools view -S failed")

    subset_header <- readLines(subset_vcf, warn = FALSE)
    meta_lines  <- subset_header[startsWith(subset_header, "##")]
    meta_lines  <- meta_lines[!grepl("ID=DS", meta_lines, fixed = TRUE)]
    chrom_line  <- subset_header[grepl("^#CHROM\t", subset_header)]
    ds_header   <- paste0('##FORMAT=<ID=DS,Number=1,Type=Float,',
                          'Description="Ancestry-specific dosage">')

    write_vcf_ds <- function(out_file, dosage_mat) {
      gt_map <- c("0/0", "0/1", "1/1")
      con <- file(out_file, open = "wt")
      on.exit(close(con), add = TRUE)
      writeLines(c(meta_lines, ds_header, chrom_line), con)
      for (i in seq_len(nrow(dosage_mat))) {
        ds     <- dosage_mat[i, ]
        gt_idx <- pmin(2L, pmax(0L, as.integer(round(ds))))
        gt     <- gt_map[gt_idx + 1L]
        sfields <- ifelse(is.na(ds), "./.:.", paste0(gt, ":", sprintf("%.6f", ds)))
        writeLines(paste(c(vcf_chrom[i], as.character(vcf_pos[i]), ".",
                           vcf_ref[i], vcf_alt[i], ".", "PASS", ".", "GT:DS",
                           sfields), collapse = "\t"), con)
      }
    }

    african_vcf  <- file.path(out_path, "african_ds.vcf")
    european_vcf <- file.path(out_path, "european_ds.vcf")
    african_gds  <- file.path(out_path, "african_ds.gds")
    european_gds <- file.path(out_path, "european_ds.gds")

    write_vcf_ds(african_vcf,  african_mat)
    write_vcf_ds(european_vcf, european_mat)

    if (verbose) message("  Converting to GDS...")
    tryCatch({
      SeqArray::seqVCF2GDS(african_vcf,  african_gds,  verbose = FALSE)
      SeqArray::seqVCF2GDS(european_vcf, european_gds, verbose = FALSE)
    }, error = function(e)
      warning("GDS conversion failed: ", e$message, ". VCFs written, GDS skipped."))

    vcf_paths <- list(african_vcf  = african_vcf,
                      european_vcf = european_vcf,
                      african_gds  = african_gds,
                      european_gds = european_gds)
    if (verbose) message("  Done: ", african_vcf, " / ", european_vcf)
  }

  if (verbose) message("\n=== Pipeline Complete ===\n")

  overlap <- list(
    n_vcf_samples          = length(vcf_samples),
    n_msp_samples          = length(msp_samples),
    n_common               = length(common_samples),
    n_variants_kept        = n_final,
    n_multiallelic_filtered = n_multi,
    n_monomorphic_filtered = n_monomorphic,
    n_no_tract             = n_no_tract,
    dropped_samples_vcf    = dropped_vcf,
    dropped_samples_msp    = dropped_msp
  )

  invisible(list(african      = african_mat,
                 european     = european_mat,
                 variant_info = variant_info,
                 sample_ids   = common_samples,
                 tract_info   = tract_df,
                 overlap      = overlap,
                 vcf_paths    = vcf_paths))
}

# ============================================================================
# run_combined_pipeline
# ============================================================================

#' Run phased and unphased ancestry splitting on the same BCF + MSP
#'
#' Parses the VCF/BCF and RFMix MSP file once, then runs both
#' deterministic per-haplotype splitting (\code{\link{split_phased}}) and
#' proportional diploid splitting (\code{\link{split_by_ancestry}}).
#' Returns both results so that they can be compared directly.
#'
#' @param vcf_path Path to a phased VCF or BCF file.
#'   \code{bcftools} must be in \code{PATH}.
#' @param msp_path Path to an RFMix MSP file (plain text or gzipped TSV).
#' @param chrom Chromosome to process (e.g. \code{"chr19"} or \code{"19"}).
#'   If \code{NULL}, all chromosomes present in the VCF are used.
#' @param verbose Print step-by-step progress messages.
#'
#' @return Invisibly, a list with two named elements:
#'   \describe{
#'     \item{\code{phased}}{Deterministic per-haplotype split via
#'       \code{\link{split_phased}}: a list with \code{african},
#'       \code{european} (variants \eqn{\times} samples integer matrices),
#'       \code{variant_info}, \code{sample_ids}, and \code{overlap}.}
#'     \item{\code{unphased}}{Proportional diploid split via
#'       \code{\link{split_by_ancestry}} using MSP-derived ancestry codes:
#'       same structure as \code{phased}.}
#'   }
#'
#' @details
#' Both splits share the same parsed MSP tracts and VCF genotypes.
#' For the \strong{unphased} path the two haplotype alleles are summed to a
#' diploid dosage (0/1/2) and per-haplotype population codes from the MSP are
#' combined into diploid codes: AFR/AFR \eqn{\to} 3, EUR/EUR \eqn{\to} 1,
#' mixed \eqn{\to} 2.  Only two-population MSP files are supported.
#'
#' Monomorphic variants are filtered independently in each mode, so the two
#' \code{variant_info} frames may differ slightly in row count.
#'
#' @examples
#' \dontrun{
#' res <- run_combined_pipeline(
#'   vcf_path = "data/chr19.phased.bcf",
#'   msp_path = "data/chr19.msp.tsv.gz",
#'   chrom    = "chr19"
#' )
#' res$phased$african      # deterministic AFR dosage matrix
#' res$unphased$african    # proportional AFR dosage matrix
#' }
#'
#' @export
run_combined_pipeline <- function(vcf_path, msp_path, chrom = NULL,
                                   verbose = TRUE) {

  if (verbose) message("=== LANTERN Combined Ancestry Pipeline ===\n")

  # ---- Step 1: Parse MSP ----
  if (verbose) message("Step 1: Parsing RFMix MSP file...")
  msp_data     <- .parse_msp(msp_path, verbose = verbose)
  msp_samples  <- msp_data$sample_ids
  tract_df     <- msp_data$tract_df
  anc_hap0_mat <- msp_data$anc_hap0
  anc_hap1_mat <- msp_data$anc_hap1
  pop_codes    <- msp_data$pop_codes
  n_tracts     <- nrow(tract_df)
  if (length(pop_codes) != 2)
    stop("run_combined_pipeline requires exactly 2 populations in the MSP; ",
         "use run_phased_pipeline() for K > 2")
  if (verbose) message("  MSP samples: ", length(msp_samples),
                       "  Tracts: ", n_tracts)

  # ---- Step 2: VCF sample IDs ----
  if (verbose) message("\nStep 2: Querying VCF sample IDs...")
  vcf_samples <- tryCatch({
    res <- system(paste0("bcftools query -l ", shQuote(vcf_path)),
                  intern = TRUE, ignore.stderr = TRUE)
    if (length(res) == 0) stop("bcftools query -l returned nothing")
    res
  }, error = function(e) stop("Could not read VCF sample IDs: ", e$message))
  if (verbose) message("  VCF samples: ", length(vcf_samples))

  # ---- Step 3: Intersect samples ----
  if (verbose) message("\nStep 3: Intersecting samples...")
  common_samples <- intersect(vcf_samples, msp_samples)
  if (length(common_samples) == 0)
    stop("No common samples between VCF and MSP")
  dropped_vcf <- setdiff(vcf_samples, common_samples)
  dropped_msp <- setdiff(msp_samples, common_samples)
  if (verbose) message("  Common: ", length(common_samples))

  anc_hap0_common <- anc_hap0_mat[common_samples, , drop = FALSE]
  anc_hap1_common <- anc_hap1_mat[common_samples, , drop = FALSE]

  # ---- Step 4: Parse VCF GT ----
  if (verbose) message("\nStep 4: Parsing VCF genotypes...")
  cmd_gt <- paste0(
    "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ",
    shQuote(vcf_path))
  gt_output <- tryCatch(
    system(cmd_gt, intern = TRUE, ignore.stderr = TRUE),
    error = function(e) stop("bcftools query failed: ", e$message))
  if (length(gt_output) == 0) stop("bcftools query returned no variants")
  if (verbose) message("  ", length(gt_output), " variant lines read")

  tmp_gt <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_gt), add = TRUE)
  writeLines(gt_output, tmp_gt)
  gt_dt <- data.table::fread(tmp_gt, header = FALSE, sep = "\t")

  vcf_chrom <- gt_dt[[1]]
  vcf_pos   <- as.integer(gt_dt[[2]])
  vcf_ref   <- as.character(gt_dt[[3]])
  vcf_alt   <- as.character(gt_dt[[4]])

  is_biallelic <- !grepl(",", vcf_alt, fixed = TRUE)
  n_multi <- sum(!is_biallelic)
  if (n_multi > 0) {
    if (verbose) message("  Skipping ", n_multi, " multiallelic variants")
    vcf_chrom <- vcf_chrom[is_biallelic]; vcf_pos <- vcf_pos[is_biallelic]
    vcf_ref   <- vcf_ref[is_biallelic];   vcf_alt <- vcf_alt[is_biallelic]
    gt_dt     <- gt_dt[is_biallelic]
  }
  n_variants <- length(vcf_pos)

  sample_idx <- match(common_samples, vcf_samples)
  gt_mat_raw <- as.matrix(gt_dt[, 5:ncol(gt_dt), drop = FALSE])
  gt_mat     <- gt_mat_raw[, sample_idx, drop = FALSE]

  gt_parsed   <- .parse_phased_gt_matrix(gt_mat, n_variants,
                                          length(common_samples),
                                          common_samples, verbose = verbose)
  gt_hap0_mat <- gt_parsed$hap0
  gt_hap1_mat <- gt_parsed$hap1

  n_missing <- sum(gt_parsed$missing)
  if (n_missing > 0 && verbose)
    message("  Missing GT: ", n_missing, " (treated as 0/0)")

  # ---- Step 5: Map variants to tracts ----
  if (verbose) message("\nStep 5: Mapping variants to ancestry tracts...")
  vcf_chrom_clean   <- sub("^chr", "", vcf_chrom)
  tract_chrom_clean <- sub("^chr", "", tract_df$chrom)

  if (!is.null(chrom)) {
    chrom_clean <- sub("^chr", "", chrom)
    keep <- vcf_chrom_clean == chrom_clean
    if (!any(keep)) stop("No variants on chromosome ", chrom)
    vcf_chrom <- vcf_chrom[keep]; vcf_pos <- vcf_pos[keep]
    vcf_ref   <- vcf_ref[keep];   vcf_alt <- vcf_alt[keep]
    gt_hap0_mat <- gt_hap0_mat[keep, , drop = FALSE]
    gt_hap1_mat <- gt_hap1_mat[keep, , drop = FALSE]
    vcf_chrom_clean <- vcf_chrom_clean[keep]
    n_variants <- sum(keep)
    if (verbose) message("  Filtered to chrom ", chrom, ": ", n_variants,
                         " variants")
  }

  tract_idx   <- findInterval(vcf_pos, tract_df$spos)
  valid_tract <- tract_idx > 0
  if (any(valid_tract))
    valid_tract[valid_tract] <-
      vcf_pos[valid_tract] <= tract_df$epos[tract_idx[valid_tract]]
  valid_chrom <- if (!is.null(chrom)) {
    rep(TRUE, length(tract_idx))
  } else {
    vapply(seq_along(tract_idx), function(i) {
      tract_idx[i] > 0 &&
        vcf_chrom_clean[i] == tract_chrom_clean[tract_idx[i]]
    }, logical(1))
  }
  keep_var   <- valid_tract & valid_chrom
  n_no_tract <- sum(!keep_var)
  if (n_no_tract > 0 && verbose)
    message("  Variants with no tract: ", n_no_tract, " (dropped)")

  if (n_no_tract > 0) {
    vcf_chrom   <- vcf_chrom[keep_var]; vcf_pos <- vcf_pos[keep_var]
    vcf_ref     <- vcf_ref[keep_var];   vcf_alt <- vcf_alt[keep_var]
    gt_hap0_mat <- gt_hap0_mat[keep_var, , drop = FALSE]
    gt_hap1_mat <- gt_hap1_mat[keep_var, , drop = FALSE]
    tract_idx   <- tract_idx[keep_var]
    n_variants  <- sum(keep_var)
  }
  if (n_variants == 0) stop("No variants with valid ancestry tracts")

  # ---- Step 6: Broadcast tract ancestry to variants ----
  if (verbose) message("\nStep 6: Broadcasting ancestry to variants...")
  anc_hap0_var <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  anc_hap1_var <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  for (i in seq_len(n_variants)) {
    t <- tract_idx[i]
    anc_hap0_var[i, ] <- anc_hap0_common[, t]
    anc_hap1_var[i, ] <- anc_hap1_common[, t]
  }

  # shared overlap stats (independent of split mode)
  overlap_base <- list(
    n_vcf_samples           = length(vcf_samples),
    n_msp_samples           = length(msp_samples),
    n_common                = length(common_samples),
    n_multiallelic_filtered = n_multi,
    n_no_tract              = n_no_tract,
    dropped_samples_vcf     = dropped_vcf,
    dropped_samples_msp     = dropped_msp
  )

  # ---- Step 7a: Phased split ----
  if (verbose) message("\nStep 7a: Splitting haplotypes by ancestry (phased)...")
  res_ph <- split_phased(gt_hap0_mat, gt_hap1_mat,
                          anc_hap0_var, anc_hap1_var,
                          pop_codes = pop_codes)
  afr_ph <- res_ph$african
  eur_ph <- res_ph$european

  has_alt_ph <- rowSums(afr_ph + eur_ph) > 0
  n_mono_ph  <- sum(!has_alt_ph)
  afr_ph <- afr_ph[has_alt_ph, , drop = FALSE]
  eur_ph <- eur_ph[has_alt_ph, , drop = FALSE]
  vi_ph  <- data.frame(chrom = vcf_chrom[has_alt_ph],
                       pos   = vcf_pos[has_alt_ph],
                       ref   = vcf_ref[has_alt_ph],
                       alt   = vcf_alt[has_alt_ph],
                       stringsAsFactors = FALSE)
  if (verbose) message("  Phased: ", nrow(afr_ph), " variants kept (",
                       n_mono_ph, " monomorphic filtered)")

  # ---- Step 7b: Unphased split ----
  if (verbose) message("\nStep 7b: Splitting diploid dosages by ancestry (unphased)...")

  gt_diploid <- gt_hap0_mat + gt_hap1_mat

  # Derive diploid ancestry codes (1=EUR/EUR, 2=AFR/EUR, 3=AFR/AFR)
  # from per-haplotype pop codes (pop_codes[1]=AFR code, pop_codes[2]=EUR code)
  afr_hap_code <- pop_codes[[1]]
  eur_hap_code <- pop_codes[[2]]
  anc_diploid <- ifelse(
    anc_hap0_var == afr_hap_code & anc_hap1_var == afr_hap_code, 3L,
    ifelse(
      anc_hap0_var == eur_hap_code & anc_hap1_var == eur_hap_code, 1L,
      2L
    )
  )

  res_un <- split_by_ancestry(gt_diploid, anc_diploid)
  afr_un <- res_un$african
  eur_un <- res_un$european

  has_alt_un <- rowSums(afr_un + eur_un) > 0
  n_mono_un  <- sum(!has_alt_un)
  afr_un <- afr_un[has_alt_un, , drop = FALSE]
  eur_un <- eur_un[has_alt_un, , drop = FALSE]
  vi_un  <- data.frame(chrom = vcf_chrom[has_alt_un],
                       pos   = vcf_pos[has_alt_un],
                       ref   = vcf_ref[has_alt_un],
                       alt   = vcf_alt[has_alt_un],
                       stringsAsFactors = FALSE)
  if (verbose) message("  Unphased: ", nrow(afr_un), " variants kept (",
                       n_mono_un, " monomorphic filtered)")

  if (verbose) message("\n=== Pipeline Complete ===\n")

  invisible(list(
    phased = list(
      african      = afr_ph,
      european     = eur_ph,
      variant_info = vi_ph,
      sample_ids   = common_samples,
      overlap      = c(overlap_base,
                       n_variants_kept        = nrow(afr_ph),
                       n_monomorphic_filtered = n_mono_ph)
    ),
    unphased = list(
      african      = afr_un,
      european     = eur_un,
      variant_info = vi_un,
      sample_ids   = common_samples,
      overlap      = c(overlap_base,
                       n_variants_kept        = nrow(afr_un),
                       n_monomorphic_filtered = n_mono_un)
    )
  ))
}
