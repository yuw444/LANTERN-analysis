#' lantern: High-Level Workflow Functions
#'
#' Convenience functions for common ancestry-specific analysis workflows.

# ============================================================================
# Internal helper functions
# ============================================================================

#' Find sample overlap between GT and PT matrices
#'
#' @param gt_matrix Genotype matrix (cols = samples)
#' @param pt_matrix PT matrix (rows = samples)
#' @return List with common_samples, dropped_from_gt, dropped_from_pt
find_sample_overlap <- function(gt_matrix, pt_matrix) {
  gt_samples <- colnames(gt_matrix)
  pt_samples <- rownames(pt_matrix)

  # If no names, assume same order
  if (is.null(gt_samples) || is.null(pt_samples)) {
    n_gt <- ncol(gt_matrix)
    n_pt <- nrow(pt_matrix)
    if (n_gt != n_pt) {
      min_n <- min(n_gt, n_pt)
      warning("Matrices have different dimensions (GT: ", n_gt, " cols, PT: ", n_pt,
              " rows). Using first ", min_n, " samples.")
      gt_samples <- paste0("sample_", 1:min_n)
      pt_samples <- paste0("sample_", 1:min_n)
      colnames(gt_matrix) <- gt_samples
      rownames(pt_matrix) <- pt_samples
    } else {
      gt_samples <- paste0("sample_", 1:ncol(gt_matrix))
      pt_samples <- paste0("sample_", 1:nrow(pt_matrix))
      colnames(gt_matrix) <- gt_samples
      rownames(pt_matrix) <- pt_samples
    }
  }

  common <- intersect(gt_samples, pt_samples)

  if (length(common) == 0) {
    stop("No common samples found between GT and PT matrices")
  }

  list(
    gt_matrix = gt_matrix[, common, drop = FALSE],
    pt_matrix = pt_matrix[common, , drop = FALSE],
    common_samples = common,
    n_common = length(common),
    dropped_from_gt = setdiff(gt_samples, common),
    dropped_from_pt = setdiff(pt_samples, common)
  )
}

#' Find variant/region overlap between GT and PT matrices
#'
#' @param gt_matrix Genotype matrix (rows = variants, cols = samples)
#' @param pt_matrix PT matrix (rows = samples, cols = regions)
#' @param vcf_path Optional VCF path for coordinate-based matching
#' @return List with aligned matrices and overlap info
find_variant_overlap <- function(gt_matrix, pt_matrix, vcf_path = NULL) {
  gt_variants <- rownames(gt_matrix)
  pt_regions <- colnames(pt_matrix)

  # If no names, assume same order
  if (is.null(gt_variants) || is.null(pt_regions)) {
    n_gt <- nrow(gt_matrix)
    n_pt <- ncol(pt_matrix)
    if (n_gt != n_pt) {
      min_n <- min(n_gt, n_pt)
      warning("Matrices have different dimensions (GT: ", n_gt, " rows, PT: ",
              n_pt, " cols). Using first ", min_n, " variants.")
      gt_variants <- paste0("var_", 1:min_n)
      pt_regions <- paste0("region_", 1:min_n)
      rownames(gt_matrix) <- gt_variants
      colnames(pt_matrix) <- pt_regions
    } else {
      gt_variants <- paste0("var_", 1:nrow(gt_matrix))
      pt_regions <- paste0("region_", 1:ncol(pt_matrix))
      rownames(gt_matrix) <- gt_variants
      colnames(pt_matrix) <- pt_regions
    }
  }

  # Try exact matching first
  common <- intersect(gt_variants, pt_regions)

  # If no exact match, try coordinate overlap
  if (length(common) == 0) {
    message("No exact variant/region name matches. Attempting coordinate-based matching...")

    # Try to parse coordinates
    # GT format: "chr:pos" (e.g., "22:123456")
    # PT format: "chr:start-end" (e.g., "22:100000-200000")

    matched_gt <- character()
    matched_pt <- character()

    for (i in seq_along(gt_variants)) {
      gt_name <- gt_variants[i]
      parts <- strsplit(gt_name, ":")[[1]]
      if (length(parts) >= 2) {
        gt_chr <- parts[1]
        gt_pos <- as.numeric(parts[2])

        for (j in seq_along(pt_regions)) {
          pt_name <- pt_regions[j]
          pt_parts <- strsplit(pt_name, ":")[[1]]
          pt_chr <- pt_parts[1]
          coords <- strsplit(pt_parts[2], "-")[[1]]
          pt_start <- as.numeric(coords[1])
          pt_end <- as.numeric(coords[2])

          if (!is.na(gt_pos) && !is.na(pt_start) && !is.na(pt_end)) {
            if (gt_chr == pt_chr && gt_pos >= pt_start && gt_pos <= pt_end) {
              matched_gt <- c(matched_gt, gt_variants[i])
              matched_pt <- c(matched_pt, pt_regions[j])
              break  # Take first matching region
            }
          }
        }
      }
    }

    if (length(matched_gt) > 0) {
      # Create unique mapping
      common <- unique(c(matched_gt, matched_pt))
      message("Found ", length(matched_gt), " variants within ", length(unique(matched_pt)), " regions.")
    }
  }

  if (length(common) == 0) {
    # Fall back to positional matching by order
    warning("No variant/region name or coordinate overlap. Matching by order.")
    min_n <- min(nrow(gt_matrix), ncol(pt_matrix))
    rownames(gt_matrix) <- paste0("var_", 1:nrow(gt_matrix))
    colnames(pt_matrix) <- paste0("region_", 1:ncol(pt_matrix))
    common <- paste0("shared_", 1:min_n)
  }

  list(
    gt_matrix = gt_matrix[intersect(rownames(gt_matrix), common), , drop = FALSE],
    pt_matrix = pt_matrix[, intersect(colnames(pt_matrix), common), drop = FALSE],
    common_variants = common,
    n_common = length(common),
    dropped_from_gt = setdiff(rownames(gt_matrix), common),
    dropped_from_pt = setdiff(colnames(pt_matrix), common)
  )
}

# ============================================================================
# Main pipeline function
# ============================================================================

#' Run ancestry splitting pipeline
#'
#' Complete pipeline: split genotype matrix by ancestry, return dosage matrices.
#' Handles sample and variant mismatches by using only overlapping data.
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
#'     \itemize{
#'       \item{n_samples_total}{Total samples before filtering}
#'       \item{n_samples_kept}{Samples after overlap filtering}
#'       \item{n_variants_total}{Total variants before filtering}
#'       \item{n_variants_kept}{Variants after overlap filtering}
#'       \item{dropped_samples}{Samples removed due to mismatch}
#'       \item{dropped_variants}{Variants removed due to mismatch}
#'     }
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

  # Step 1: Find sample overlap
  if (verbose) message("Step 1: Finding sample overlap...")

  overlap <- find_sample_overlap(gt_matrix, pt_matrix)
  gt_matrix <- overlap$gt_matrix
  pt_matrix <- overlap$pt_matrix

  if (verbose) {
    message("  Samples: ", n_gt_samples_orig, " in GT, ", n_pt_samples_orig, " in PT")
    message("  Common samples: ", overlap$n_common)
    if (length(overlap$dropped_from_gt) > 0) {
      message("  Dropped from GT: ", paste(head(overlap$dropped_from_gt, 3), collapse = ", "),
              if (length(overlap$dropped_from_gt) > 3) "..." else "")
    }
    if (length(overlap$dropped_from_pt) > 0) {
      message("  Dropped from PT: ", paste(head(overlap$dropped_from_pt, 3), collapse = ", "),
              if (length(overlap$dropped_from_pt) > 3) "..." else "")
    }
  }

  # Step 2: Find variant/region overlap
  if (verbose) message("\nStep 2: Finding variant/region overlap...")

  var_overlap <- find_variant_overlap(gt_matrix, pt_matrix, vcf_path)
  gt_matrix <- var_overlap$gt_matrix
  pt_matrix <- var_overlap$pt_matrix

  if (verbose) {
    message("  Variants in GT: ", n_gt_variants_orig)
    message("  Regions in PT: ", n_pt_regions_orig)
    message("  Common variants/regions: ", var_overlap$n_common)
    if (length(var_overlap$dropped_from_gt) > 0) {
      message("  Dropped variants: ", paste(head(var_overlap$dropped_from_gt, 3), collapse = ", "),
              if (length(var_overlap$dropped_from_gt) > 3) "..." else "")
    }
  }

  # Validate dimensions match after overlap
  if (ncol(gt_matrix) != nrow(pt_matrix)) {
    stop("After overlap filtering: GT cols (", ncol(gt_matrix),
         ") != PT rows (", nrow(pt_matrix), ")")
  }
  if (nrow(gt_matrix) != ncol(pt_matrix)) {
    stop("After overlap filtering: GT rows (", nrow(gt_matrix),
         ") != PT cols (", ncol(pt_matrix), "). ",
         "Each GT variant should match to exactly one PT region.")
  }

  # Step 3: Split genotypes by ancestry
  if (verbose) message("\nStep 3: Splitting genotypes by ancestry...")

  result <- split_by_ancestry(gt_matrix, pt_matrix)

  if (verbose) {
    message("  -> African dosage matrix: ", nrow(result$african), " x ", ncol(result$african))
    message("  -> European dosage matrix: ", nrow(result$european), " x ", ncol(result$european))
  }

  # Step 4: Count ancestries per region
  if (verbose) message("\nStep 4: Counting ancestries per region...")

  counts <- list(
    african = count_ancestry_codes(pt_matrix, 3),
    european = count_ancestry_codes(pt_matrix, 1),
    mixed = count_ancestry_codes(pt_matrix, 2)
  )

  if (verbose) {
    message("  -> African (3): median = ", median(counts$african))
    message("  -> European (1): median = ", median(counts$european))
    message("  -> Mixed (2): median = ", median(counts$mixed))
  }

  # Step 5: Write VCFs if requested
  if (!is.null(vcf_path) && !is.null(out_prefix)) {
    if (verbose) message("\nStep 5: Writing ancestry-specific VCFs...")

    .write_vcf_with_ancestry(
      vcf_path, gt_matrix, pt_matrix,
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
      n_samples_total = n_gt_samples_orig + n_pt_samples_orig - overlap$n_common,
      n_samples_kept = overlap$n_common,
      n_variants_total = n_gt_variants_orig + n_pt_regions_orig - var_overlap$n_common,
      n_variants_kept = var_overlap$n_common,
      dropped_samples = unique(c(overlap$dropped_from_gt, overlap$dropped_from_pt)),
      dropped_variants = unique(c(var_overlap$dropped_from_gt, var_overlap$dropped_from_pt))
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
