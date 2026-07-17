# ============================================================================
# write_ancestry_gds (Step 2): materialize a Step 1 split result to GDS
# ============================================================================
#
# write_dosage_gds()/write_ancestry_gds(): the GDS writers used by Steps 2-3.

# ============================================================================
# write_dosage_gds
# ============================================================================

#' Write an ancestry-specific dosage matrix to a SeqArray GDS file
#'
#' Converts a dosage matrix (rows = variants, columns = samples) to a
#' SeqArray GDS file via a temporary VCF with a \code{DS} FORMAT field.
#' The resulting GDS can be passed directly to \code{GMMAT::SMMAT()} via the
#' \code{gds.fn} argument with \code{is.dosage = TRUE}.
#'
#' @param dosage_mat Numeric matrix; rows = variants, columns = samples.
#'   Values should be in \eqn{[0, 2]}.
#' @param variant_info Data frame with columns \code{chrom}, \code{pos},
#'   \code{ref}, \code{alt} (one row per variant in the same order as rows of
#'   \code{dosage_mat}).  Returned directly by \code{\link{ancestry_split}}.
#' @param sample_ids Character vector of sample identifiers, length equal to
#'   \code{ncol(dosage_mat)}.
#' @param gds_path Output file path for the GDS file (e.g. \code{"afr.gds"}).
#'
#' @return Invisibly returns \code{gds_path}.
#'
#' @seealso \code{\link{ancestry_split}}, \code{\link{write_ancestry_gds}},
#'   \code{\link{cauchy_combine}}
#'
#' @examples
#' \dontrun{
#' split <- ancestry_split("cohort.bcf", "cohort.msp.tsv", mode = "haplotype")
#' write_dosage_gds(split$AFR, split$variant_info, split$sample_ids, "afr_phased.gds")
#' write_dosage_gds(split$EUR, split$variant_info, split$sample_ids, "eur_phased.gds")
#' }
#'
#' @export
write_dosage_gds <- function(dosage_mat, variant_info, sample_ids, gds_path) {
  if (!requireNamespace("SeqArray", quietly = TRUE))
    stop("Package 'SeqArray' is required. Install with: BiocManager::install('SeqArray')")

  n_var  <- nrow(dosage_mat)
  n_samp <- ncol(dosage_mat)
  if (n_var != nrow(variant_info))
    stop("nrow(dosage_mat) must equal nrow(variant_info)")
  if (n_samp != length(sample_ids))
    stop("ncol(dosage_mat) must equal length(sample_ids)")

  vcf_path <- paste0(tools::file_path_sans_ext(gds_path), "_tmp_ds.vcf")
  on.exit(unlink(vcf_path), add = TRUE)

  header <- c(
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste0('##FORMAT=<ID=DS,Number=1,Type=Float,',
           'Description="Ancestry-specific dosage of the alternate allele">'),
    paste(c("#CHROM", "POS", "ID", "REF", "ALT",
            "QUAL", "FILTER", "INFO", "FORMAT", sample_ids),
          collapse = "\t")
  )

  # Vectorised field construction (avoids per-row sprintf loop)
  gt_map <- c("0/0", "0/1", "1/1")
  gt_idx <- matrix(pmin(2L, pmax(0L, as.integer(round(dosage_mat)))),
                   nrow = n_var, ncol = n_samp)
  gt_mat <- matrix(gt_map[gt_idx + 1L], nrow = n_var, ncol = n_samp)

  # Use compact integer format when all values are 0/1/2; otherwise 4 d.p.
  if (all(dosage_mat == floor(dosage_mat), na.rm = TRUE)) {
    ds_mat <- matrix(as.character(as.integer(dosage_mat)), nrow = n_var, ncol = n_samp)
  } else {
    ds_mat <- matrix(sprintf("%.4f", dosage_mat), nrow = n_var, ncol = n_samp)
  }

  sfields <- matrix(paste0(gt_mat, ":", ds_mat), nrow = n_var, ncol = n_samp)
  if (anyNA(dosage_mat)) sfields[is.na(dosage_mat)] <- "./.:."

  # Prefix columns (CHROM POS ID REF ALT QUAL FILTER INFO FORMAT)
  pre <- cbind(variant_info$chrom,
               as.character(variant_info$pos),
               ".", variant_info$ref, variant_info$alt,
               ".", "PASS", ".", "GT:DS")

  # Combine and write in one shot
  data_lines <- apply(cbind(pre, sfields), 1L, paste, collapse = "\t")
  writeLines(c(header, data_lines), vcf_path)

  SeqArray::seqVCF2GDS(vcf_path, gds_path, verbose = FALSE)
  invisible(gds_path)
}

# ============================================================================
# write_ancestry_gds (Step 2 entry point)
# ============================================================================

#' Write all populations of an ancestry split to GDS files (Step 2)
#'
#' Takes the list returned by \code{\link{ancestry_split}} (or the
#' backward-compatible \code{\link{ancestry_split_dosage}} /
#' \code{\link{ancestry_split_phased}} output, which has the same
#' per-population matrix + \code{variant_info} + \code{sample_ids} shape)
#' and writes one GDS file per population via \code{\link{write_dosage_gds}},
#' printing progress as it goes. This is Step 2 of the package's
#' Step1/Step2/Step3 workflow: most users only need to call it so they have
#' GDS paths to pass to \code{\link{ancestry_smmat}} (Step 3); it's exported
#' mainly for standalone use (e.g. caching GDS files between SLURM jobs).
#'
#' @param split_result List as returned by \code{\link{ancestry_split}}:
#'   named per-population dosage matrices plus \code{variant_info} and
#'   \code{sample_ids}.
#' @param out_path Output directory (created if it doesn't exist). One file
#'   \code{<population>.gds} is written per population.
#' @param verbose Print per-population progress messages.
#'
#' @return Named list of GDS file paths, one per population, e.g.
#'   \code{list(AFR = "out/AFR.gds", EUR = "out/EUR.gds")}. Pass this
#'   directly as \code{gds_paths} to \code{\link{ancestry_smmat}}.
#'
#' @seealso \code{\link{ancestry_split}}, \code{\link{ancestry_smmat}}
#'
#' @examples
#' \dontrun{
#' split <- ancestry_split("cohort.bcf", "cohort.msp.tsv.gz", mode = "dosage")
#' gds   <- write_ancestry_gds(split, "out/")
#' }
#'
#' @export
write_ancestry_gds <- function(split_result, out_path, verbose = TRUE) {
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

  pop_names <- setdiff(names(split_result),
                        c("variant_info", "sample_ids", "mode",
                          "ancestry_counts", "tract_info", "overlap"))
  if (length(pop_names) == 0)
    stop("split_result has no per-population dosage matrices to write")

  variant_info <- split_result$variant_info
  sample_ids   <- split_result$sample_ids

  if (verbose) message("=== LANTERN Write Ancestry GDS (Step 2) ===\n")

  gds_paths <- vector("list", length(pop_names))
  names(gds_paths) <- pop_names
  for (i in seq_along(pop_names)) {
    pop <- pop_names[i]
    if (verbose)
      message("Writing GDS for ", pop, " (", i, "/", length(pop_names), ")...")
    gds_path <- file.path(out_path, paste0(pop, ".gds"))
    write_dosage_gds(split_result[[pop]], variant_info, sample_ids, gds_path)
    gds_paths[[pop]] <- gds_path
  }

  if (verbose)
    message("\n=== Step 2 Complete: ", length(pop_names),
            " GDS file(s) written to ", out_path, " ===\n")

  gds_paths
}
