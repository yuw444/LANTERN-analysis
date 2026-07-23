# ============================================================================
# Generate lantern:::centromeres_hg38 (internal package data, R/sysdata.rda)
# ============================================================================
#
# Source: UCSC REST API, hg38 "centromeres" track (modeled centromere
# sequences; each chromosome may have several contig segments spanning the
# centromeric gap). We collapse each chromosome to a single p-arm/q-arm
# boundary: p_end = min(chromStart) across all segments, q_start =
# max(chromEnd) across all segments. Autosomes only (chr1-chr22) -- LANTERN's
# rare-variant burden pipeline is autosomal only.
#
# Re-run this script (`pixi run Rscript lantern/data-raw/centromeres_hg38.R`
# from the repo root) to refresh R/sysdata.rda if UCSC revises the hg38
# centromere models.

url <- "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=centromeres"
resp <- jsonlite::fromJSON(url, simplifyVector = FALSE)
cents <- resp$centromeres

chroms <- paste0("chr", 1:22)
centromeres_hg38 <- do.call(rbind, lapply(chroms, function(chrom) {
  segs <- cents[[chrom]]
  if (is.null(segs) || length(segs) == 0)
    stop("No centromere data for ", chrom, " -- UCSC track may have changed")
  starts <- vapply(segs, function(s) s$chromStart, numeric(1))
  ends   <- vapply(segs, function(s) s$chromEnd,   numeric(1))
  data.frame(chrom = chrom, p_end = min(starts), q_start = max(ends))
}))
rownames(centromeres_hg38) <- NULL

# Sanity check: p_end < q_start for every autosome
stopifnot(all(centromeres_hg38$p_end < centromeres_hg38$q_start))

usethis_dir <- file.path("lantern", "R")
if (!dir.exists(usethis_dir)) usethis_dir <- "R"  # allow running from lantern/ directly
save(centromeres_hg38, file = file.path(usethis_dir, "sysdata.rda"), version = 2)

message("Wrote ", nrow(centromeres_hg38), " autosomes to ",
        file.path(usethis_dir, "sysdata.rda"))
