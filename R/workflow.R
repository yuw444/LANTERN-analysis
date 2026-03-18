#' lantern: High-Level Workflow Functions
#'
#' Convenience functions for common ancestry-specific analysis workflows.

#' Run ancestry splitting pipeline
#'
#' Complete pipeline: split genotype matrix by ancestry, return dosage matrices.
#'
#' @param gt_matrix Integer matrix of genotypes (rows=variants, cols=samples)
#'   Values: 0=homozygous ref, 1=heterozygous, 2=homozygous alt
#' @param pt_matrix Integer matrix of parent-of-origin ancestry codes
#'   (rows=samples, cols=variants/regions)
#'   Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR
#' @param out_prefix Output file prefix for VCFs (optional)
#' @param vcf_path Input VCF path for VCF output (optional)
#'
#' @return List with elements:
#'   \item{african}{African ancestry-specific dosage matrix}
#'   \item{european}{European ancestry-specific dosage matrix}
#'   \item{counts}{List of ancestry counts per region}
#'
#' @examples
#' # Example with 3 variants and 4 samples
#' gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1, 0, 1, 1), nrow = 3, ncol = 4)
#' pt <- matrix(c(3, 2, 1, 3, 2, 1, 2, 2, 1, 1, 3, 1), nrow = 4, ncol = 3)
#'
#' result <- run_ancestry_pipeline(gt, pt)
#' result$african   # African ancestry dosages
#' result$european  # European ancestry dosages
#'
#' @export
run_ancestry_pipeline <- function(gt_matrix, pt_matrix,
                                  out_prefix = NULL, vcf_path = NULL) {

  message("=== LANTERN Ancestry Pipeline ===\n")

  # Validate dimensions
  if (ncol(gt_matrix) != nrow(pt_matrix)) {
    stop("GT matrix cols (", ncol(gt_matrix),
         ") != PT matrix rows (", nrow(pt_matrix), ")")
  }

  # Step 1: Split genotypes by ancestry
  message("Step 1: Splitting genotypes by ancestry...")
  result <- split_by_ancestry(gt_matrix, pt_matrix)
  message("  -> African dosage matrix: ", nrow(result$african), " x ", ncol(result$african))
  message("  -> European dosage matrix: ", nrow(result$european), " x ", ncol(result$european))

  # Step 2: Count ancestries per region
  message("\nStep 2: Counting ancestries per region...")
  counts <- list(
    african = count_ancestry_codes(pt_matrix, 3),
    european = count_ancestry_codes(pt_matrix, 1),
    mixed = count_ancestry_codes(pt_matrix, 2)
  )
  message("  -> African (3): median = ", median(counts$african))
  message("  -> European (1): median = ", median(counts$european))
  message("  -> Mixed (2): median = ", median(counts$mixed))

  # Step 3: Write VCFs if requested
  if (!is.null(vcf_path) && !is.null(out_prefix)) {
    message("\nStep 3: Writing ancestry-specific VCFs...")
    .write_vcf_with_ancestry(
      vcf_path, gt_matrix, pt_matrix,
      paste0(out_prefix, "_african.vcf"),
      paste0(out_prefix, "_european.vcf")
    )
    message("  -> African: ", out_prefix, "_african.vcf")
    message("  -> European: ", out_prefix, "_european.vcf")
  }

  message("\n=== Pipeline Complete ===\n")

  invisible(list(
    african = result$african,
    european = result$european,
    counts = counts
  ))
}

#' Calculate ancestry weights per gene
#'
#' Calculate median ancestry counts for samples within each gene's region.
#'
#' @param pt_matrix PT (ancestry) matrix (rows=samples, cols=regions)
#' @param sample_ids Character vector of sample IDs (must match rownames of pt_matrix)
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
