#' lantern: High-Level Workflow Functions
#'
#' Convenience functions for common ancestry-specific analysis workflows.

# ============================================================================
# Main pipeline function
# ============================================================================

#' Run ancestry splitting pipeline
#'
#' Complete pipeline: split genotype matrix by ancestry, return dosage matrices.
#' Handles sample and variant mismatches by using only overlapping data.
#' Memory efficient: single-pass subsetting using index vectors.
#'
#' @param gt_matrix Integer matrix of genotypes (rows=variants, cols=samples).
#'   Values: 0=homozygous ref, 1=heterozygous, 2=homozygous alt.
#'   rownames: variant IDs (e.g., \code{"chr:pos"}). colnames: sample IDs.
#'   May be \code{NULL} when \code{vcf_path} and \code{msp_path} are supplied.
#' @param pt_matrix Integer matrix of parent-of-origin ancestry codes
#'   (rows=samples, cols=regions/variants).
#'   Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR.
#'   rownames: sample IDs (must match colnames of \code{gt_matrix}).
#'   colnames: region IDs (e.g., \code{"chr:start-end"}).
#'   May be \code{NULL} when \code{vcf_path} and \code{msp_path} are supplied.
#' @param vcf_path Path to a phased VCF or BCF file.  When supplied together
#'   with \code{msp_path}, the function builds the genotype and ancestry
#'   matrices directly from these files (bcftools must be in PATH), bypassing
#'   the \code{gt_matrix}/\code{pt_matrix} arguments.
#' @param msp_path Path to an RFMix MSP file (plain text or gzipped TSV).
#'   Used together with \code{vcf_path} as a shortcut input.
#' @param chrom Chromosome name to restrict to (e.g. \code{"chr19"}).
#'   Only used when \code{vcf_path}/\code{msp_path} are supplied.
#' @param verbose Print progress messages (default TRUE)
#'
#' @return List with elements:
#'   \item{african}{African ancestry-specific dosage matrix}
#'   \item{european}{European ancestry-specific dosage matrix}
#'   \item{counts}{List of ancestry counts per region/sample}
#'   \item{overlap}{Overlap information}
#'
#' @examples
#' # Example with mismatched dimensions
#' gt <- matrix(c(2, 1, 0, 1, 2, 1), nrow = 2, ncol = 3,
#'              dimnames = list(c("chr1:100", "chr1:200"),
#'                             c("sample_A", "sample_B", "sample_C")))
#'
#' pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
#'              dimnames = list(c("sample_A", "sample_B", "sample_D"),
#'                             c("chr1:50-150", "chr1:150-250")))
#'
#' result <- run_ancestry_pipeline(gt, pt)
#' # Only sample_A and sample_B are kept (common to both)
#' # Only chr1:100 and chr1:200 that overlap regions are kept
#'
#' \dontrun{
#' # Shortcut: build matrices directly from VCF + MSP
#' result <- run_ancestry_pipeline(vcf_path = "data/chr19.phased.bcf",
#'                                  msp_path = "data/chr19.msp.tsv.gz",
#'                                  chrom    = "chr19")
#' }
#'
#' @export
run_ancestry_pipeline <- function(gt_matrix = NULL, pt_matrix = NULL,
                                   vcf_path = NULL, msp_path = NULL,
                                   chrom = NULL, verbose = TRUE) {

  # ---- VCF/MSP shortcut path ----
  if (!is.null(vcf_path) && !is.null(msp_path)) {
    return(invisible(.run_ancestry_from_vcf_msp(vcf_path, msp_path,
                                                 chrom, verbose)))
  }
  if (is.null(gt_matrix) || is.null(pt_matrix))
    stop("Provide either (gt_matrix, pt_matrix) or (vcf_path, msp_path)")

  if (verbose) message("=== LANTERN Ancestry Pipeline ===\n")

  # Record original dimensions
  n_gt_samples_orig <- ncol(gt_matrix)
  n_gt_variants_orig <- nrow(gt_matrix)
  n_pt_samples_orig <- nrow(pt_matrix)
  n_pt_regions_orig <- ncol(pt_matrix)

  # ========================================================================
  # Step 1: Find sample overlap (memory efficient - use indices, not copies)
  # ========================================================================
  if (verbose) message("Step 1: Finding sample overlap...")

  gt_sample_names <- colnames(gt_matrix)
  pt_sample_names <- rownames(pt_matrix)

  # Handle unnamed matrices
  if (is.null(gt_sample_names) || is.null(pt_sample_names)) {
    n_gt <- ncol(gt_matrix)
    n_pt <- nrow(pt_matrix)
    min_n <- min(n_gt, n_pt)
    if (n_gt != n_pt) {
      warning("Matrices have different dimensions (GT: ", n_gt, " cols, PT: ", n_pt,
              " rows). Using first ", min_n, " samples.")
    }
    gt_sample_names <- paste0("sample_", 1:min_n)
    pt_sample_names <- paste0("sample_", 1:min_n)
    colnames(gt_matrix) <- gt_sample_names
    rownames(pt_matrix) <- pt_sample_names
  }

  # Find common samples by name
  common_samples <- intersect(gt_sample_names, pt_sample_names)

  if (length(common_samples) == 0) {
    stop("No common samples found between GT and PT matrices")
  }

  # Calculate indices for subsetting (avoid copying until needed)
  gt_sample_idx <- match(common_samples, gt_sample_names)
  pt_sample_idx <- match(common_samples, pt_sample_names)

  n_common_samples <- length(common_samples)
  dropped_samples <- c(
    setdiff(gt_sample_names, common_samples),
    setdiff(pt_sample_names, common_samples)
  )

  if (verbose) {
    message("  Samples: ", n_gt_samples_orig, " in GT, ", n_pt_samples_orig, " in PT")
    message("  Common samples: ", n_common_samples)
    if (length(dropped_samples) > 0) {
      message("  Dropped: ", paste(head(dropped_samples, 3), collapse = ", "),
              if (length(dropped_samples) > 3) "..." else "")
    }
  }

  # ========================================================================
  # Step 2: Find variant/region overlap
  # ========================================================================
  if (verbose) message("\nStep 2: Finding variant/region overlap...")

  gt_var_names <- rownames(gt_matrix)
  pt_region_names <- colnames(pt_matrix)

  # Handle unnamed matrices
  if (is.null(gt_var_names) || is.null(pt_region_names)) {
    n_gt <- nrow(gt_matrix)
    n_pt <- ncol(pt_matrix)
    min_n <- min(n_gt, n_pt)
    if (n_gt != n_pt) {
      warning("Matrices have different dimensions (GT: ", n_gt, " rows, PT: ",
              n_pt, " cols). Using first ", min_n, " variants.")
    }
    gt_var_names <- paste0("var_", 1:min_n)
    pt_region_names <- paste0("var_", 1:min_n)
    rownames(gt_matrix) <- gt_var_names
    colnames(pt_matrix) <- pt_region_names
  }

  # Try exact matching first (fast)
  common_vars <- intersect(gt_var_names, pt_region_names)

  # If no exact match, try coordinate-based matching
  if (length(common_vars) == 0) {
    if (verbose) message("  No exact matches. Trying coordinate-based matching...")

    # Pre-parse GT coordinates: chr:pos -> list(chr, pos)
    parse_gt_coord <- function(name) {
      parts <- strsplit(name, ":")[[1]]
      if (length(parts) >= 2) {
        pos <- suppressWarnings(as.numeric(parts[2]))
        if (is.na(pos)) pos <- NULL
        list(chr = parts[1], pos = pos)
      } else {
        NULL
      }
    }

    # Pre-parse PT coordinates: chr:start-end -> list(chr, start, end)
    parse_pt_coord <- function(name) {
      chr_part <- strsplit(name, ":")[[1]]
      if (length(chr_part) >= 2) {
        coords <- strsplit(chr_part[2], "-")[[1]]
        if (length(coords) == 2) {
          start <- suppressWarnings(as.numeric(coords[1]))
          end <- suppressWarnings(as.numeric(coords[2]))
          if (!is.na(start) && !is.na(end)) {
            return(list(chr = chr_part[1], start = start, end = end))
          }
        }
      }
      NULL
    }

    # Parse all coordinates upfront (vectorized approach)
    gt_coords <- lapply(gt_var_names, parse_gt_coord)
    pt_coords <- lapply(pt_region_names, parse_pt_coord)

    # Find valid matches using matrix operations where possible
    matched_gt_idx <- integer()
    matched_pt_idx <- integer()

    # Build match matrix for O(n*m) scan
    for (i in seq_along(gt_coords)) {
      if (is.null(gt_coords[[i]]$pos)) next
      for (j in seq_along(pt_coords)) {
        pt <- pt_coords[[j]]
        if (is.null(pt$start)) next
        if (gt_coords[[i]]$chr == pt$chr &&
            gt_coords[[i]]$pos >= pt$start &&
            gt_coords[[i]]$pos <= pt$end) {
          matched_gt_idx <- c(matched_gt_idx, i)
          matched_pt_idx <- c(matched_pt_idx, j)
          break  # Take first matching region
        }
      }
    }

    # Clean up intermediate objects
    rm(gt_coords, pt_coords)
    gc()

    if (length(matched_gt_idx) > 0) {
      common_vars <- unique(c(gt_var_names[matched_gt_idx], pt_region_names[matched_pt_idx]))
      if (verbose) {
        message("  Found ", length(matched_gt_idx), " variants in ",
                length(unique(matched_pt_idx)), " regions")
      }
    } else {
      # Fallback to positional match by order
      warning("No variant/region name or coordinate overlap. Matching by order.")
      min_n <- min(length(gt_var_names), length(pt_region_names))
      rownames(gt_matrix) <- paste0("var_", seq_len(nrow(gt_matrix)))
      colnames(pt_matrix) <- paste0("region_", seq_len(ncol(pt_matrix)))
      common_vars <- paste0("shared_", 1:min_n)
    }
  }

  # Calculate indices for variant/region subsetting
  gt_var_idx <- match(intersect(gt_var_names, common_vars), gt_var_names)
  pt_region_idx <- match(intersect(pt_region_names, common_vars), pt_region_names)

  n_common_vars <- length(common_vars)
  dropped_variants <- c(
    setdiff(gt_var_names, common_vars),
    setdiff(pt_region_names, common_vars)
  )

  if (verbose) {
    message("  Variants in GT: ", n_gt_variants_orig)
    message("  Regions in PT: ", n_pt_regions_orig)
    message("  Common variants/regions: ", n_common_vars)
    if (length(dropped_variants) > 0) {
      message("  Dropped: ", paste(head(dropped_variants, 3), collapse = ", "),
              if (length(dropped_variants) > 3) "..." else "")
    }
  }

  # ========================================================================
  # Step 3: Single-pass subsetting (memory efficient)
  # ========================================================================
  if (verbose) message("\nStep 3: Subsetting matrices...")

  # Subset both matrices in single operation using calculated indices
  # GT: rows = variants we want, cols = samples we want
  # PT: rows = samples we want, cols = regions we want
  gt_subset <- gt_matrix[gt_var_idx, gt_sample_idx, drop = FALSE]
  pt_subset <- pt_matrix[pt_sample_idx, pt_region_idx, drop = FALSE]

  # Clean up intermediates immediately
  rm(gt_matrix, pt_matrix)
  gc()

  # ========================================================================
  # Step 4: Validate dimensions
  # ========================================================================
  if (ncol(gt_subset) != nrow(pt_subset)) {
    stop("After overlap filtering: GT cols (", ncol(gt_subset),
         ") != PT rows (", nrow(pt_subset), ")")
  }
  if (nrow(gt_subset) != ncol(pt_subset)) {
    stop("After overlap filtering: GT rows (", nrow(gt_subset),
         ") != PT cols (", ncol(pt_subset), "). ",
         "Each GT variant should match to exactly one PT region.")
  }

  # ========================================================================
  # Step 4: Filter monomorphic variants (no alt alleles)
  # ========================================================================
  if (verbose) message("\nStep 3b: Filtering monomorphic variants...")

  variant_has_alt <- rowSums(gt_subset > 0) > 0
  n_monomorphic <- sum(!variant_has_alt)

  if (n_monomorphic > 0) {
    if (verbose) message("  Removed ", n_monomorphic, " monomorphic variants (no alt alleles)")

    dropped_variants <- c(dropped_variants, rownames(gt_subset)[!variant_has_alt])

    gt_subset <- gt_subset[variant_has_alt, , drop = FALSE]
    pt_subset <- pt_subset[, variant_has_alt, drop = FALSE]
  }

  if (nrow(gt_subset) == 0) {
    stop("All variants are monomorphic (no alt alleles)")
  }

  # ========================================================================
  # Step 5: Split genotypes by ancestry (C backend)
  # ========================================================================
  if (verbose) message("\nStep 4: Splitting genotypes by ancestry...")

  result <- split_by_ancestry(gt_subset, pt_subset)

  # Capture post-filter variant count BEFORE removing gt_subset
  n_variants_kept_final <- nrow(gt_subset)

  # Clean up gt_subset after C call (keep pt_subset for counts - it's small)
  rm(gt_subset)
  gc()

  if (verbose) {
    message("  -> African dosage matrix: ", nrow(result$african), " x ", ncol(result$african))
    message("  -> European dosage matrix: ", nrow(result$european), " x ", ncol(result$european))
  }

  # ========================================================================
  # Step 6: Count ancestries per region
  # ========================================================================
  if (verbose) message("\nStep 5: Counting ancestries per region...")

  counts <- list(
    african = count_ancestry_codes(pt_subset, 3),
    european = count_ancestry_codes(pt_subset, 1),
    mixed = count_ancestry_codes(pt_subset, 2)
  )

  if (verbose) {
    message("  -> African (3): median = ", median(counts$african))
    message("  -> European (1): median = ", median(counts$european))
    message("  -> Mixed (2): median = ", median(counts$mixed))
  }

  if (verbose) message("\n=== Pipeline Complete ===\n")

  invisible(list(
    african = result$african,
    european = result$european,
    counts = counts,
    overlap = list(
      n_samples_total = n_gt_samples_orig + n_pt_samples_orig - n_common_samples,
      n_samples_kept = n_common_samples,
      n_variants_total = n_gt_variants_orig,
      n_variants_kept = n_variants_kept_final,
      n_monomorphic_filtered = n_monomorphic,
      dropped_samples = unique(dropped_samples),
      dropped_variants = unique(dropped_variants)
    )
  ))
}

# ============================================================================
# Utility functions
# ============================================================================

#' Calculate ancestry weights per gene
#'
#' Calculate median ancestry counts for samples within each gene's region.
#'
#' @param pt_matrix PT (ancestry) matrix (rows=samples, cols=regions)
#' @param sample_ids Character vector of sample IDs
#' @param gene_data Data frame with columns: gene, chr, start, end
#'
#' @return Data frame with gene_id, n_afr, n_eur, n_mixed columns
#' @export
calc_gene_ancestry_weights <- function(pt_matrix, sample_ids, gene_data) {

  if (!is.null(rownames(pt_matrix)) && !all(sample_ids %in% rownames(pt_matrix))) {
    stop("sample_ids must match rownames of pt_matrix")
  }

  results <- vector("list", nrow(gene_data))

  for (i in seq_len(nrow(gene_data))) {
    gene <- gene_data[i, ]

    results[[i]] <- data.frame(
      gene_id = gene$gene,
      n_afr = sum(pt_matrix[sample_ids == gene$sample, ] == 3, na.rm = TRUE),
      n_eur = sum(pt_matrix[sample_ids == gene$sample, ] == 1, na.rm = TRUE),
      n_mixed = sum(pt_matrix[sample_ids == gene$sample, ] == 2, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Create ancestry-specific VCF files
#'
#' @param vcf_path Input VCF path
#' @param gt_matrix Genotype matrix
#' @param pt_matrix PT (ancestry) matrix
#' @param out_african Output path for African-specific VCF
#' @param out_european Output path for European-specific VCF
#'
#' @export
create_ancestry_vcfs <- function(vcf_path, gt_matrix, pt_matrix,
                                  out_african, out_european) {

  message("Creating ancestry-specific VCFs...")

  .write_vcf_with_ancestry(
    vcf_path, gt_matrix, pt_matrix,
    out_african, out_european
  )

  message("  African: ", out_african)
  message("  European: ", out_european)

  invisible(c(african = out_african, european = out_european))
}

# ============================================================================
# Internal: unphased pipeline from VCF + MSP
# ============================================================================

.run_ancestry_from_vcf_msp <- function(vcf_path, msp_path, chrom, verbose) {

  if (verbose) message("=== LANTERN Ancestry Pipeline (VCF/MSP input) ===\n")

  # ---- Step 1: Parse MSP ----
  if (verbose) message("Step 1: Parsing RFMix MSP file...")
  msp_data     <- .parse_msp(msp_path, verbose = verbose)
  msp_samples  <- msp_data$sample_ids
  tract_df     <- msp_data$tract_df
  anc_hap0_mat <- msp_data$anc_hap0
  anc_hap1_mat <- msp_data$anc_hap1
  pop_codes    <- msp_data$pop_codes
  if (length(pop_codes) != 2)
    stop(".run_ancestry_from_vcf_msp requires exactly 2 populations in the MSP")
  if (verbose) message("  MSP samples: ", length(msp_samples),
                       "  Tracts: ", nrow(tract_df))

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

  # ---- Step 6: Build diploid GT and ancestry matrices ----
  if (verbose) message("\nStep 6: Building diploid genotype and ancestry matrices...")

  gt_diploid <- gt_hap0_mat + gt_hap1_mat  # (n_variants x n_samples), values 0/1/2

  # Broadcast tract ancestry to per-variant ancestry
  afr_hap_code <- pop_codes[[1]]
  eur_hap_code <- pop_codes[[2]]
  anc_diploid <- matrix(0L, nrow = n_variants, ncol = length(common_samples))
  for (i in seq_len(n_variants)) {
    h0 <- anc_hap0_common[, tract_idx[i]]
    h1 <- anc_hap1_common[, tract_idx[i]]
    anc_diploid[i, ] <- ifelse(h0 == afr_hap_code & h1 == afr_hap_code, 3L,
                        ifelse(h0 == eur_hap_code & h1 == eur_hap_code, 1L, 2L))
  }
  colnames(gt_diploid)  <- common_samples
  colnames(anc_diploid) <- common_samples

  # ---- Step 7: Filter monomorphic variants ----
  if (verbose) message("\nStep 7: Filtering monomorphic variants...")
  has_alt <- rowSums(gt_diploid) > 0
  n_mono  <- sum(!has_alt)
  if (n_mono > 0) {
    if (verbose) message("  Removed ", n_mono, " monomorphic variants")
    gt_diploid  <- gt_diploid[has_alt, , drop = FALSE]
    anc_diploid <- anc_diploid[has_alt, , drop = FALSE]
    vcf_chrom <- vcf_chrom[has_alt]; vcf_pos <- vcf_pos[has_alt]
    vcf_ref   <- vcf_ref[has_alt];   vcf_alt <- vcf_alt[has_alt]
  }
  n_final <- length(vcf_pos)
  if (n_final == 0) stop("All variants are monomorphic after filtering")
  rownames(gt_diploid)  <- paste0(vcf_chrom, ":", vcf_pos)
  rownames(anc_diploid) <- rownames(gt_diploid)

  # ---- Step 8: Split by ancestry ----
  if (verbose) message("\nStep 8: Splitting genotypes by ancestry...")
  res <- split_by_ancestry(gt_diploid, anc_diploid)
  if (verbose) {
    message("  -> African dosage matrix: ", nrow(res$african), " x ",
            ncol(res$african))
    message("  -> European dosage matrix: ", nrow(res$european), " x ",
            ncol(res$european))
  }

  # ---- Step 9: Count ancestries per sample ----
  if (verbose) message("\nStep 9: Counting ancestries per sample...")
  counts <- list(
    african  = count_ancestry_codes(t(anc_diploid), 3L),
    european = count_ancestry_codes(t(anc_diploid), 1L),
    mixed    = count_ancestry_codes(t(anc_diploid), 2L)
  )
  if (verbose) {
    message("  -> African (3): median = ", median(counts$african))
    message("  -> European (1): median = ", median(counts$european))
    message("  -> Mixed (2): median = ", median(counts$mixed))
  }

  if (verbose) message("\n=== Pipeline Complete ===\n")

  list(
    african  = res$african,
    european = res$european,
    counts   = counts,
    overlap  = list(
      n_samples_total          = length(vcf_samples) + length(msp_samples) -
                                   length(common_samples),
      n_samples_kept           = length(common_samples),
      n_variants_total         = n_variants + n_mono,
      n_variants_kept          = n_final,
      n_multiallelic_filtered  = n_multi,
      n_monomorphic_filtered   = n_mono,
      n_no_tract               = n_no_tract,
      dropped_samples          = unique(c(dropped_vcf, dropped_msp))
    )
  )
}
