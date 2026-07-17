#' lantern: Ancestry-Specific Rare Variant Analysis
#'
#' R package for ancestry-specific rare variant association analysis
#' with pure C backend for performance-critical operations.
#'
#' @section Step 1 — split: \code{\link{ancestry_split}}
#' @section Step 2 — write GDS: \code{\link{write_ancestry_gds}}
#' @section Step 3 — associate: \code{\link{ancestry_smmat}}
#'
#' @docType package
#' @name lantern-package
#' @useDynLib lantern, .registration = TRUE
#' @importFrom utils head
#' @importFrom stats median
NULL
