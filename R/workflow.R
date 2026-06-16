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
#' @param gt_matrix Integer matrix of genotypes (rows=variants, cols=samples)
#'   Values: 0=homozygous ref, 1=heterozygous, 2=homozygous alt
#'   rownames: variant IDs (e.g., "chr:pos" or any identifier)
#'   colnames: sample IDs
#' @param pt_matrix Integer matrix of parent-of-origin ancestry codes
#'   (rows=samples, cols=regions)
#'   Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR
#'   rownames: sample IDs (must match colnames of gt_matrix)
#'   colnames: region IDs (e.g., "chr:start-end")
#' @param vcf_path Input VCF path (optional, for coordinate-based variant matching)
#' @param out_prefix Output file prefix for VCFs (optional)
#' @param verbose Print progress messages (default TRUE)
#'
#' @return List with elements:
#'   \item{african}{African ancestry-specific dosage matrix}
#'   \item{european}{European ancestry-specific dosage matrix}
#'   \item{counts}{List of ancestry counts per region}
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
#' @export
run_ancestry_pipeline <- function(gt_matrix, pt_matrix,
                                  vcf_path = NULL,
                                  out_prefix = NULL,
                                  verbose = TRUE) {

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
    pt_region_names <- paste0("region_", 1:min_n)
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

  # ========================================================================
  # Step 7: Write VCFs if requested
  # ========================================================================
  if (!is.null(vcf_path) && !is.null(out_prefix)) {
    if (verbose) message("\nStep 6: Writing ancestry-specific VCFs...")

    .write_vcf_with_ancestry(
      vcf_path, result$african, result$european, pt_subset,
      paste0(out_prefix, "_african.vcf"),
      paste0(out_prefix, "_european.vcf")
    )

    if (verbose) {
      message("  -> African: ", out_prefix, "_african.vcf")
      message("  -> European: ", out_prefix, "_european.vcf")
    }
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
