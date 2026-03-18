#' lantern: Ancestry-Specific Rare Variant Analysis
#'
#' R package for ancestry-specific rare variant association analysis
#' with pure C backend for performance-critical operations.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{count_ancestry_codes()}}{Count ancestry codes in matrix}
#'   \item{\code{split_by_ancestry()}}{Split genotype by ancestry}
#'   \item{\code{read_bed_file()}}{Read PLINK bed file}
#'   \item{\code{write_vcf_with_ancestry()}}{Write VCF with ancestry annotations}
#'   \item{\code{subset_vcf_by_range()}}{Subset VCF by genomic range}
#' }
#'
#' @docType package
#' @name lantern-package
NULL

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

.read_bed_file <- function(bed_path, bim_path, fam_path, sample_indices = NULL) {
    .Call("read_bed_file_C", bed_path, bim_path, fam_path, sample_indices, PACKAGE = "lantern")
}

.write_vcf_with_ancestry <- function(vcf_path, gt_matrix, ancestry_matrix, 
                                      output_african, output_european) {
    .Call("write_vcf_with_ancestry_C", vcf_path, gt_matrix, ancestry_matrix, 
          output_african, output_european, PACKAGE = "lantern")
}

.subset_vcf_by_range <- function(vcf_path, chrom, start, end, output_path) {
    .Call("subset_vcf_by_range_C", vcf_path, chrom, start, end, output_path, PACKAGE = "lantern")
}

#' Count ancestry codes in matrix
#' @param mat Integer matrix with ancestry codes (1, 2, 3)
#' @param code Code to count
#' @export
count_ancestry_codes <- function(mat, code) {
    .count_ancestry_codes(mat, code)
}

#' Split genotype matrix by ancestry
#' @param gt_genotype Integer matrix of genotypes (0, 1, 2)
#' @param ancestry Integer matrix of ancestry codes (1, 2, 3)
#' @export
split_by_ancestry <- function(gt_genotype, ancestry) {
    .split_by_ancestry(gt_genotype, ancestry)
}

#' Read PLINK bed file into matrix
#' @param bed_path Path to .bed file
#' @param bim_path Path to .bim file
#' @param fam_path Path to .fam file
#' @param sample_indices Optional vector of sample indices to keep
#' @export
read_bed_file <- function(bed_path, bim_path, fam_path, sample_indices = NULL) {
    .read_bed_file(bed_path, bim_path, fam_path, sample_indices)
}

#' Write VCF files with ancestry-specific dosages
#' @param vcf_path Input VCF path
#' @param gt_matrix Genotype matrix
#' @param ancestry_matrix Ancestry code matrix
#' @param output_african Output VCF for African ancestry
#' @param output_european Output VCF for European ancestry
#' @export
write_vcf_with_ancestry <- function(vcf_path, gt_matrix, ancestry_matrix, 
                                     output_african, output_european) {
    .write_vcf_with_ancestry(vcf_path, gt_matrix, ancestry_matrix, 
                             output_african, output_european)
}

#' Subset VCF by genomic range
#' @param vcf_path Input VCF path
#' @param chrom Chromosome
#' @param start Start position
#' @param end End position
#' @param output_path Output VCF path
#' @export
subset_vcf_by_range <- function(vcf_path, chrom, start, end, output_path) {
    .subset_vcf_by_range(vcf_path, chrom, start, end, output_path)
}
