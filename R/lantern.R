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

#' Split unphased genotype matrix by ancestry for K populations
#'
#' Generalises \code{\link{split_by_ancestry}} to an arbitrary number of
#' ancestries.  For each variant, population allele proportions
#' \eqn{p_1, \ldots, p_K} (summing to 1) are estimated from samples with
#' unambiguous ancestry (pure-ancestry homozygotes and heterozygotes, plus
#' mixed hom-alts which each contribute one allele to each parent pool).
#' Ambiguous heterozygotes in each mixed pair are then split by the
#' conditional pairwise ratio \eqn{p_i / (p_i + p_j)}.
#' When no unambiguous alt alleles exist for a pair (singleton edge case),
#' the split defaults to 0.5 / 0.5.
#'
#' @param gt_genotype Integer matrix (variants x samples) of genotypes (0/1/2).
#' @param ancestry Integer matrix (variants x samples) of ancestry codes,
#'   same dimensions as \code{gt_genotype}.
#' @param pure_codes Named integer vector mapping each population label to its
#'   pure-ancestry diploid code.
#'   Example: \code{c(AFR = 3L, EUR = 1L, NAT = 4L)}.
#' @param mixed_codes Data frame with three columns:
#'   \describe{
#'     \item{code}{Integer ancestry code for this mixed pair.}
#'     \item{pop1}{Character name of the first parent population
#'       (must be a name in \code{pure_codes}).}
#'     \item{pop2}{Character name of the second parent population
#'       (must be a name in \code{pure_codes}).}
#'   }
#'   One row per mixed-ancestry diploid type.
#'   Example for 3 populations:
#'   \code{data.frame(code=c(2L,5L,6L), pop1=c("AFR","AFR","EUR"), pop2=c("EUR","NAT","NAT"))}
#'
#' @return Named list of K numeric matrices (variants x samples), one per
#'   population in \code{pure_codes}, in the same order.  Each element is
#'   named after the corresponding entry in \code{pure_codes}.
#'
#' @seealso \code{\link{split_by_ancestry}} for the hardcoded 2-population
#'   (AFR + EUR) wrapper; \code{\link{split_phased_multi}} for the phased
#'   K-population analogue.
#'
#' @examples
#' # 3-population example: AFR=3, EUR=1, NAT=4; mixed codes 2/5/6
#' set.seed(1)
#' gt  <- matrix(sample(0:2, 20, replace = TRUE, prob = c(.8,.15,.05)),
#'               nrow = 4, ncol = 5)
#' # ancestry codes: 3=AFR/AFR, 1=EUR/EUR, 4=NAT/NAT,
#' #                 2=AFR/EUR, 5=AFR/NAT, 6=EUR/NAT
#' anc <- matrix(sample(c(1L,2L,3L,4L,5L,6L), 20, replace = TRUE),
#'               nrow = 4, ncol = 5)
#' pure  <- c(AFR = 3L, EUR = 1L, NAT = 4L)
#' mixed <- data.frame(code = c(2L, 5L, 6L),
#'                     pop1 = c("AFR", "AFR", "EUR"),
#'                     pop2 = c("EUR", "NAT", "NAT"))
#' out <- split_by_ancestry_multi(gt, anc, pure, mixed)
#' names(out)   # "AFR" "EUR" "NAT"
#'
#' @export
split_by_ancestry_multi <- function(gt_genotype, ancestry,
                                     pure_codes, mixed_codes) {
    if (is.null(names(pure_codes)))
        stop("pure_codes must be a named integer vector, e.g. c(AFR=3L, EUR=1L)")
    if (!is.data.frame(mixed_codes) ||
        !all(c("code", "pop1", "pop2") %in% names(mixed_codes)))
        stop("mixed_codes must be a data.frame with columns: code, pop1, pop2")
    pop_names <- names(pure_codes)
    bad <- setdiff(c(mixed_codes$pop1, mixed_codes$pop2), pop_names)
    if (length(bad))
        stop("mixed_codes references populations not in pure_codes: ",
             paste(bad, collapse = ", "))

    storage.mode(gt_genotype) <- "integer"
    storage.mode(ancestry)    <- "integer"
    storage.mode(pure_codes)  <- "integer"

    m_code <- as.integer(mixed_codes$code)
    m_pop1 <- match(mixed_codes$pop1, pop_names) - 1L   # 0-based index
    m_pop2 <- match(mixed_codes$pop2, pop_names) - 1L

    result <- .Call("split_by_ancestry_multi_C",
                    gt_genotype, ancestry,
                    pure_codes, m_code, m_pop1, m_pop2,
                    PACKAGE = "lantern")
    names(result) <- pop_names
    result
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
