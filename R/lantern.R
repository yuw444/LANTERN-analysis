#' lantern: Ancestry-Specific Rare Variant Analysis
#'
#' R package for ancestry-specific rare variant association analysis
#' with pure C backend for performance-critical operations.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{count_ancestry_codes()}}{Count ancestry codes in matrix}
#'   \item{\code{split_by_ancestry()}}{Split genotype by ancestry}
#'   \item{\code{write_vcf_with_ancestry()}}{Write VCF with ancestry annotations}
#' }
#'
#' @docType package
#' @name lantern-package
#' @useDynLib lantern, .registration = TRUE
#' @importFrom utils head
#' @importFrom stats median
NULL

# ============================================================================
# Internal C wrappers
# ============================================================================

.count_ancestry_codes <- function(mat, code) {
    storage.mode(mat) <- "integer"
    code <- as.integer(code)
    .Call("count_ancestry_codes_C", mat, code, PACKAGE = "lantern")
}

.split_by_ancestry <- function(gt_genotype, ancestry) {
    storage.mode(gt_genotype) <- "integer"
    storage.mode(ancestry) <- "integer"
    .Call("split_by_ancestry_C", gt_genotype, ancestry, PACKAGE = "lantern")
}

.write_vcf_with_ancestry <- function(vcf_path, gt_matrix, ancestry_matrix,
                                      output_african, output_european) {
    .Call("write_vcf_with_ancestry_C", vcf_path, gt_matrix, ancestry_matrix,
          output_african, output_european, PACKAGE = "lantern")
}

.subset_vcf_by_range <- function(vcf_path, chrom, start, end, output_path) {
    .Call("subset_vcf_by_range_C", vcf_path, chrom, start, end, output_path, PACKAGE = "lantern")
}

# ============================================================================
# Exported Functions
# ============================================================================

#' Count ancestry codes in matrix
#'
#' Count occurrences of a specific ancestry code in each row of a matrix.
#'
#' @param mat Integer matrix with ancestry codes (1, 2, 3)
#'   - Rows represent genomic regions/windows
#'   - Columns represent samples
#' @param code Ancestry code to count:
#'   - 1 = EUR/EUR (Pure European)
#'   - 2 = AFR/EUR (Mixed)
#'   - 3 = AFR/AFR (Pure African)
#'
#' @return Integer vector with count of \code{code} per row
#'
#' @examples
#' # PT matrix: 3 regions x 4 samples
#' pt <- matrix(c(3, 2, 1, 3, 2, 2, 1, 1, 3, 3, 2, 1), nrow = 3, ncol = 4)
#' count_ancestry_codes(pt, code = 3)  # Count African ancestry
#'
#' @export
count_ancestry_codes <- function(mat, code) {
    .count_ancestry_codes(mat, code)
}

#' Split genotype matrix by ancestry
#'
#' Split a genotype matrix into African and European ancestry-specific
#' dosage matrices based on a parent-of-origin ancestry matrix.
#'
#' @param gt_genotype Integer matrix of genotypes (0, 1, 2)
#'   - Rows represent variants
#'   - Columns represent samples
#' @param ancestry Integer matrix of ancestry codes (1, 2, 3)
#'   - Same dimensions as gt_genotype
#'   - Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR
#'
#' @return List with two elements:
#'   \item{african}{Matrix of African ancestry-specific dosages}
#'   \item{european}{Matrix of European ancestry-specific dosages}
#'
#' @examples
#' # Genotype matrix: 5 variants x 4 samples
#' gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1, 0,
#'                1, 1, 1, 0, 2, 1, 0, 1, 0, 1), nrow = 5, ncol = 4)
#'
#' # Ancestry matrix: same dimensions
#' ancestry <- matrix(c(3, 2, 1, 3, 2, 1, 2, 2, 1, 1,
#'                      3, 1, 2, 1, 3, 2, 1, 2, 1, 3), nrow = 5, ncol = 4)
#'
#' result <- split_by_ancestry(gt, ancestry)
#' result$african
#' result$european
#'
#' @export
split_by_ancestry <- function(gt_genotype, ancestry) {
    .split_by_ancestry(gt_genotype, ancestry)
}

#' Write VCF files with ancestry-specific dosages
#'
#' Write separate VCF files for African and European ancestry-specific
#' variant dosages.
#'
#' @param vcf_path Input VCF path (for header template)
#' @param gt_matrix Genotype matrix (rows=variants, cols=samples)
#' @param ancestry_matrix Ancestry code matrix (same dimensions)
#' @param output_african Output VCF path for African ancestry
#' @param output_european Output VCF path for European ancestry
#'
#' @export
write_vcf_with_ancestry <- function(vcf_path, gt_matrix, ancestry_matrix,
                                     output_african, output_european) {
    .write_vcf_with_ancestry(vcf_path, gt_matrix, ancestry_matrix,
                             output_african, output_european)
}

#' Subset VCF by genomic range
#'
#' Extract variants within a genomic region from a VCF file.
#'
#' @param vcf_path Input VCF path
#' @param chrom Chromosome (e.g., "22" or "chr22")
#' @param start Start position
#' @param end End position
#' @param output_path Output VCF path
#'
#' @export
subset_vcf_by_range <- function(vcf_path, chrom, start, end, output_path) {
    .subset_vcf_by_range(vcf_path, chrom, start, end, output_path)
}
