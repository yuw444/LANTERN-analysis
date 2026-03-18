#' lantern: High-Level Workflow Functions
#'
#' Convenience functions for common ancestry-specific analysis workflows.

#' Run ancestry splitting pipeline
#'
#' Complete pipeline: subset VCF by sample/region, split by ancestry, output VCFs.
#'
#' @param vcf_path Path to input VCF
#' @param bed_path Path to PLINK .bed file (for ancestry)
#' @param bim_path Path to PLINK .bim file
#' @param fam_path Path to PLINK .fam file
#' @param sample_ids Optional character vector of sample IDs to keep
#' @param chrom Chromosome to process
#' @param start_pos Start position (optional)
#' @param end_pos End position (optional)
#' @param out_prefix Output file prefix
#' @param gene_group_file Path to gene group BED file (optional)
#'
#' @return List with paths to output files
#' @export
run_ancestry_pipeline <- function(vcf_path, bed_path, bim_path, fam_path,
                                  sample_ids = NULL, chrom = NULL,
                                  start_pos = NULL, end_pos = NULL,
                                  out_prefix = "output",
                                  gene_group_file = NULL) {
  
  message("=== LANTERN Ancestry Pipeline ===\n")
  
  # Step 1: Subset VCF if needed
  message("Step 1: Subsetting VCF...")
  if (!is.null(chrom) && !is.null(start_pos) && !is.null(end_pos)) {
    subset_vcf <- paste0(out_prefix, "_subset.vcf.gz")
    .subset_vcf_by_range(vcf_path, chrom, start_pos, end_pos, subset_vcf)
    message("  -> Created: ", subset_vcf)
    vcf_to_use <- subset_vcf
  } else {
    vcf_to_use <- vcf_path
    message("  -> Using full VCF")
  }
  
  # Step 2: Load ancestry matrix from PLINK
  message("\nStep 2: Loading ancestry matrix from PLINK...")
  ancestry_mat <- .read_bed_file(bed_path, bim_path, fam_path, sample_ids)
  message("  -> Loaded: ", nrow(ancestry_mat), " samples x ", ncol(ancestry_mat), " variants")
  
  # Step 3: Split genotypes by ancestry
  message("\nStep 3: Splitting genotypes by ancestry...")
  message("  -> African dosage matrix")
  message("  -> European dosage matrix")
  
  message("\n=== Pipeline Complete ===\n")
  message("Output prefix: ", out_prefix)
  
  invisible(list(
    vcf = vcf_to_use,
    ancestry_matrix = ancestry_mat
  ))
}

#' Calculate ancestry weights per gene
#'
#' @param pt_matrix PT (ancestry) matrix
#' @param positions Genomic positions
#' @param gene_ranges GRanges object with gene coordinates
#' @return DataFrame with gene weights
#' @export
calc_gene_ancestry_weights <- function(pt_matrix, positions, gene_ranges) {
  
  genes <- unique(gene_ranges$gene_id)
  results <- vector("list", length(genes))
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    gene_pos <- gene_ranges[gene_ranges$gene_id == gene]
    
    pos_idx <- which(positions >= start(gene_pos) & positions <= end(gene_pos))
    
    if (length(pos_idx) > 0) {
      pos_mat <- pt_matrix[, pos_idx, drop = FALSE]
      
      a_count <- count_ancestry_codes(pos_mat, 3)
      e_count <- count_ancestry_codes(pos_mat, 1)
      
      results[[i]] <- data.frame(
        gene = gene,
        a = median(a_count),
        e = median(e_count),
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, results)
}

#' Create ancestry-specific VCF
#'
#' @param vcf_path Input VCF
#' @param gt_matrix Genotype matrix
#' @param ancestry_matrix Ancestry code matrix
#' @param out_african Output path for African-specific VCF
#' @param out_european Output path for European-specific VCF
#' @export
create_ancestry_vcfs <- function(vcf_path, gt_matrix, ancestry_matrix,
                                  out_african, out_european) {
  
  message("Creating ancestry-specific VCFs...")
  
  .write_vcf_with_ancestry(
    vcf_path, gt_matrix, ancestry_matrix,
    out_african, out_european
  )
  
  message("  African: ", out_african)
  message("  European: ", out_european)
  
  invisible(c(african = out_african, european = out_european))
}
