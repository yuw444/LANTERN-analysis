# ============================================================================
# ancestry_smmat (Step 3): SMMAT association testing + Cauchy combination
# ============================================================================

#' Per-gene ancestry weights for the Cauchy combination (Step 3 helper)
#'
#' Generalises the retired \code{src/step3_weight_finding.R} script's
#' median per-gene ancestry count logic to K populations (that
#' functionality now lives here, fed by \code{ancestry_counts}/
#' \code{variant_info} from Step 1): for each gene in \code{gene_group_file},
#' takes the median (across the gene's variants) of the per-variant pure-
#' ancestry sample counts in \code{ancestry_counts}.
#'
#' @param ancestry_counts Numeric matrix (variants x populations), as
#'   returned by \code{\link{ancestry_split}}.
#' @param variant_info data.frame with columns chrom/pos, row-aligned to
#'   \code{ancestry_counts}, as returned by \code{\link{ancestry_split}}.
#' @param gene_group_file Path to the gene_group file (no header; columns
#'   gene, chr, pos, ref, alt, weight) — the same file passed to
#'   \code{GMMAT::SMMAT()}.
#' @param pop_names Character vector of population names to produce weight
#'   columns for (typically \code{names(gds_paths)} from \code{\link{ancestry_smmat}}).
#'
#' @return Numeric matrix (genes x \code{pop_names}), or \code{NULL} if
#'   none of \code{pop_names} match \code{colnames(ancestry_counts)}.
#'
#' @keywords internal
.compute_gene_weights <- function(ancestry_counts, variant_info,
                                   gene_group_file, pop_names) {
  gg <- data.table::fread(gene_group_file, header = FALSE,
                           col.names = c("gene", "chr", "pos", "ref", "alt", "weight"))

  avail_pops <- intersect(pop_names, colnames(ancestry_counts))
  if (length(avail_pops) == 0) {
    warning("None of the population names (", paste(pop_names, collapse = ", "),
            ") match ancestry_counts columns (",
            paste(colnames(ancestry_counts), collapse = ", "),
            "); falling back to equal weights.")
    return(NULL)
  }

  vi_key <- paste0(sub("^chr", "", as.character(variant_info$chrom)), ":",
                    variant_info$pos)
  gg_key <- paste0(sub("^chr", "", as.character(gg$chr)), ":", gg$pos)
  row_idx <- match(gg_key, vi_key)

  genes_unique <- unique(gg$gene)
  weights_mat  <- matrix(1, nrow = length(genes_unique), ncol = length(pop_names),
                          dimnames = list(genes_unique, pop_names))

  for (gene in genes_unique) {
    idx <- row_idx[gg$gene == gene]
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) next
    for (pop in avail_pops) {
      weights_mat[gene, pop] <- stats::median(ancestry_counts[idx, pop])
    }
  }

  weights_mat
}

#' Run ancestry-stratified SMMAT + Cauchy combination (Step 3)
#'
#' Wraps \code{GMMAT::SMMAT()} across all ancestry-split GDS files (as
#' produced by \code{\link{write_ancestry_gds}}), computes per-gene
#' ancestry weights from Step 1's \code{ancestry_counts} (this replaces the
#' retired \code{src/step3_weight_finding.R} script), and combines the per-gene p-values
#' across populations with \code{\link{cauchy_combine}}. This is the
#' package's Step 3 entry point; it does not run Step 1/2 for you (SMMAT
#' runs can take hours, so this is always an explicit, separate call).
#'
#' @param gds_paths Named list/vector of GDS file paths, one per ancestry
#'   population, e.g. \code{c(AFR = "afr.gds", EUR = "eur.gds")} — as
#'   returned by \code{\link{write_ancestry_gds}}. Any K >= 1 populations
#'   are supported; add an extra named entry (e.g. \code{OBSERVED = ...})
#'   to include an additional (e.g. non-ancestry-split) GDS in the same
#'   SMMAT + Cauchy combination run.
#' @param pheno data.frame with an \code{id_col} column plus the response
#'   and covariates referenced in \code{formula}.
#' @param formula Model formula for the null model, e.g. \code{y ~ age + sex}.
#' @param kinship Numeric matrix (or object coercible via \code{as.matrix})
#'   with row/column names equal to sample IDs.
#' @param gene_group_file Either a path to the gene_group file passed to
#'   \code{GMMAT::SMMAT()} (no header; columns gene, chr, pos, ref, alt,
#'   weight), or a data.frame/matrix with those same six columns (written to
#'   a temporary file internally). Either way, the file is deleted with
#'   \code{unlink()} when \code{ancestry_smmat()} returns — no
#'   \code{tempfile()}/\code{write.table()}/\code{unlink()} needed by the
#'   caller, but a \emph{path you want to keep must not be passed here}; pass
#'   a fresh copy (or the data.frame) each time, including across repeated
#'   calls with the same gene set (e.g. one call per phenotype).
#' @param ancestry_counts,variant_info Optional: pass
#'   \code{ancestry_split()$ancestry_counts} / \code{$variant_info}
#'   straight through to compute per-gene ancestry weights for the Cauchy
#'   combination. If omitted, equal weights are used.
#' @param id_col Name of the sample ID column in \code{pheno} and
#'   \code{formula}'s \code{id} argument to \code{glmmkin} (default \code{"id"}).
#' @param family A family object for \code{GMMAT::glmmkin()} (default
#'   \code{gaussian()}).
#' @param gene_col Name of the gene/group ID column in \code{SMMAT()}'s
#'   output (default \code{"group"}, matching \code{GMMAT::SMMAT()}).
#' @param p_col Name of the p-value column in \code{SMMAT()}'s output to
#'   combine (default \code{"E.pval"}, matching \code{GMMAT::SMMAT()}'s
#'   default \code{tests = "E"}; use \code{"O.pval"}/\code{"S.pval"}/
#'   \code{"B.pval"} if you pass a different \code{tests} via \code{...}).
#' @param MAF.range,miss.cutoff,method,ncores Passed straight through to
#'   \code{GMMAT::SMMAT()} — same names, same meaning as in a direct
#'   \code{SMMAT()} call.
#' @param verbose Print this function's own step-by-step progress messages
#'   (Steps 1-6 below). Distinct from \code{GMMAT::SMMAT()}'s own internal
#'   \code{verbose} (its per-call progress bar), which is not exposed
#'   separately and stays at SMMAT's default (\code{FALSE}); every other
#'   \code{SMMAT()} argument is available via \code{...}.
#' @param ... Any other \code{GMMAT::SMMAT()} argument, forwarded as-is —
#'   for example \code{tests} (\code{"B"}, \code{"S"}, \code{"O"}, or the
#'   default \code{"E"}), \code{rho}, \code{MAF.weights.beta},
#'   \code{missing.method}, \code{use.minor.allele}, \code{auto.flip},
#'   \code{Garbage.Collection}, \code{group.file.sep}, or
#'   \code{meta.file.prefix}. \code{null.obj}, \code{geno.file},
#'   \code{group.file}, and \code{is.dosage} are always supplied internally
#'   by \code{ancestry_smmat()} and cannot be overridden this way. See
#'   \code{?GMMAT::SMMAT} for the full list and defaults.
#'
#' @return A list with two elements:
#'   \item{results}{data.frame with one row per gene: \code{gene}, one
#'     \code{p_<POP>} column per population, one \code{w_<POP>} weight
#'     column per population (only if \code{ancestry_counts}/
#'     \code{variant_info} were supplied), and \code{p_cauchy}.}
#'   \item{smmat_results}{Named list (one entry per population, same names
#'     as \code{gds_paths}) of the raw \code{GMMAT::SMMAT()} result
#'     data.frames.}
#'
#' @seealso \code{\link{ancestry_split}}, \code{\link{write_ancestry_gds}},
#'   \code{\link{cauchy_combine}}
#'
#' @examples
#' \dontrun{
#' split <- ancestry_split("cohort.bcf", "cohort.msp.tsv.gz", mode = "dosage")
#' gds   <- write_ancestry_gds(split, "out/")
#'
#' # gene_group_file can be a data.frame -- no tempfile/unlink needed
#' genes <- data.frame(gene = "MYGENE", chr = "19", pos = split$variant_info$pos,
#'                      ref = split$variant_info$ref, alt = split$variant_info$alt,
#'                      weight = 1)
#' out <- ancestry_smmat(gds, pheno, y ~ age + sex, kinship, genes,
#'                        ancestry_counts = split$ancestry_counts,
#'                        variant_info    = split$variant_info)
#' head(out$results)         # per-gene p-values + Cauchy combination
#' out$smmat_results$AFR     # raw GMMAT::SMMAT() output for AFR
#' }
#'
#' @export
ancestry_smmat <- function(gds_paths, pheno, formula, kinship, gene_group_file,
                            ancestry_counts = NULL, variant_info = NULL,
                            id_col = "id", family = stats::gaussian(),
                            gene_col = "group", p_col = "E.pval",
                            MAF.range = c(0, 0.5), miss.cutoff = 1,
                            method = "davies", ncores = 1, verbose = TRUE, ...) {

  if (!requireNamespace("GMMAT", quietly = TRUE))
    stop("Package 'GMMAT' is required for ancestry_smmat(). Install with: ",
         "install.packages('GMMAT') (not available via conda/Bioconductor).")

  if (!is.character(gene_group_file)) {
    if (!is.data.frame(gene_group_file) && !is.matrix(gene_group_file))
      stop("gene_group_file must be a file path (character) or a ",
           "data.frame/matrix with columns gene, chr, pos, ref, alt, weight")
    gg_tmp <- tempfile(fileext = ".tsv")
    utils::write.table(gene_group_file, gg_tmp, sep = "\t", row.names = FALSE,
                        col.names = FALSE, quote = FALSE)
    gene_group_file <- gg_tmp
  }
  # gene_group_file is unlinked when this call returns, whether it was
  # supplied as a path or written internally from a data.frame/matrix above.
  on.exit(unlink(gene_group_file), add = TRUE)

  gds_paths <- as.list(gds_paths)
  pop_names <- names(gds_paths)
  if (is.null(pop_names) || any(pop_names == ""))
    stop("gds_paths must be a named list/vector, e.g. c(AFR = ..., EUR = ...)")
  K <- length(pop_names)

  if (verbose) message("=== LANTERN Association Testing (Step 3) ===\n")

  # ---- Step 1: Sample overlap across pheno, kinship, and all GDS files ----
  if (verbose) message("Step 1: Intersecting sample IDs...")
  if (!id_col %in% names(pheno))
    stop("pheno must contain an id column named '", id_col, "'")
  if (is.null(rownames(kinship)) || is.null(colnames(kinship)))
    stop("kinship must have row and column names corresponding to sample IDs")

  gds_sample_ids <- lapply(gds_paths, function(p) {
    g <- SeqArray::seqOpen(p)
    on.exit(SeqArray::seqClose(g), add = TRUE)
    SeqArray::seqGetData(g, "sample.id")
  })

  ids_common <- Reduce(intersect, c(list(pheno[[id_col]], rownames(kinship)),
                                     gds_sample_ids))
  if (length(ids_common) == 0)
    stop("No overlapping sample IDs found among phenotype, kinship, and GDS files")
  if (verbose) message("  Common samples: ", length(ids_common))

  pheno_sub <- pheno[match(ids_common, pheno[[id_col]]), , drop = FALSE]
  kin_idx   <- match(ids_common, rownames(kinship))
  kin_sub   <- as.matrix(kinship[kin_idx, kin_idx, drop = FALSE])

  # ---- Step 2: Fit the null model once ----
  if (verbose) message("\nStep 2: Fitting null model (glmmkin)...")
  model0 <- GMMAT::glmmkin(formula, data = pheno_sub, kins = kin_sub,
                            id = id_col, family = family)

  # ---- Step 3: Run SMMAT for each population ----
  smmat_results <- vector("list", K)
  names(smmat_results) <- pop_names
  for (i in seq_len(K)) {
    pop <- pop_names[i]
    if (verbose) message("\nStep 3: Running SMMAT for ", pop,
                         " (", i, "/", K, ")...")
    smmat_results[[pop]] <- GMMAT::SMMAT(
      model0, gds_paths[[pop]], gene_group_file,
      MAF.range = MAF.range, miss.cutoff = miss.cutoff,
      method = method, is.dosage = TRUE, ncores = ncores, ...)
  }

  # ---- Step 4: Merge per-gene p-values ----
  if (verbose) message("\nStep 4: Merging per-gene p-values...")
  for (pop in pop_names) {
    if (!gene_col %in% names(smmat_results[[pop]]))
      stop("Column '", gene_col, "' not found in SMMAT() output for ", pop)
    if (!p_col %in% names(smmat_results[[pop]]))
      stop("Column '", p_col, "' not found in SMMAT() output for ", pop,
           "; available columns: ",
           paste(names(smmat_results[[pop]]), collapse = ", "))
  }
  gene_ids <- Reduce(union, lapply(smmat_results, function(r) as.character(r[[gene_col]])))
  p_mat <- matrix(NA_real_, nrow = length(gene_ids), ncol = K,
                   dimnames = list(gene_ids, pop_names))
  for (pop in pop_names) {
    r <- smmat_results[[pop]]
    p_mat[as.character(r[[gene_col]]), pop] <- r[[p_col]]
  }

  # ---- Step 5: Per-gene ancestry weights ----
  weights_mat <- NULL
  if (!is.null(ancestry_counts) && !is.null(variant_info)) {
    if (verbose) message("\nStep 5: Computing per-gene ancestry weights...")
    weights_mat <- .compute_gene_weights(ancestry_counts, variant_info,
                                          gene_group_file, pop_names)
  } else if (verbose) {
    message("\nStep 5: ancestry_counts/variant_info not supplied; ",
            "Cauchy combination will use equal weights.")
  }

  # ---- Step 6: Cauchy-combine p-values per gene ----
  if (verbose) message("\nStep 6: Cauchy-combining p-values per gene...")
  p_cauchy <- vapply(seq_along(gene_ids), function(i) {
    p_row <- p_mat[i, ]
    keep  <- !is.na(p_row)
    if (sum(keep) < 2) return(NA_real_)
    w <- if (!is.null(weights_mat) && gene_ids[i] %in% rownames(weights_mat)) {
      weights_mat[gene_ids[i], pop_names][keep]
    } else NULL
    cauchy_combine(p_row[keep], weights = w)
  }, numeric(1))

  # ---- Step 7: Assemble result ----
  results <- data.frame(gene = gene_ids, stringsAsFactors = FALSE)
  for (pop in pop_names) results[[paste0("p_", pop)]] <- p_mat[, pop]
  if (!is.null(weights_mat)) {
    row_match <- match(gene_ids, rownames(weights_mat))
    for (pop in pop_names) {
      results[[paste0("w_", pop)]] <- unname(weights_mat[row_match, pop])
    }
  }
  results$p_cauchy <- p_cauchy
  rownames(results) <- NULL

  if (verbose) message("\n=== Step 3 Complete: ", nrow(results), " genes ===\n")

  list(results = results, smmat_results = smmat_results)
}
