library(optparse)

option_list <- list(
    make_option(
        c("--msp_path"),
        type = "character",
        default = NULL,
        help = "Path to RFMix MSP file (plain text or gzipped TSV; see vignette('generating-local-ancestry'))"
    ),
    make_option(
        c("-i", "--vcf_path"),
        type = "character",
        default = NULL,
        help = "Path to phased VCF/BCF (bgzipped) with genotypes to split by ancestry"
    ),
    make_option(
        c("-o", "--out_path"),
        type = "character",
        default = NULL,
        help = "Output directory (created if missing)"
    ),
    make_option(
        c("-c", "--chr_id"),
        type = "character",
        default = NULL,
        help = "Chromosome to process (must match both the MSP file and the VCF CHROM column, e.g. '22')"
    ),
    make_option(
        c("--mode"),
        type = "character",
        default = "dosage",
        help = "Split algorithm: 'dosage' (proportional p1/p2 split, unphased; ambiguous mixed-ancestry hets shrunk toward 1/2) or 'haplotype' (deterministic per-haplotype split, requires phased VCF) [default %default]"
    )
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("msp_path", "vcf_path", "out_path", "chr_id")
missing_opts <- required[vapply(opt[required], is.null, logical(1))]
if (length(missing_opts) > 0) {
    stop(
        "Missing required option(s): --", paste(missing_opts, collapse = ", --"),
        "\nSee: Rscript src/step1_vcf_split_by_ancestry.R --help"
    )
}

suppressPackageStartupMessages(library(lantern))

dir.create(opt$out_path, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Step 1: Split the VCF by local ancestry (lantern::ancestry_split())
# ============================================================================
split_result <- ancestry_split(
    vcf_path = opt$vcf_path,
    msp_path = opt$msp_path,
    mode     = opt$mode,
    chrom    = opt$chr_id,
    verbose  = TRUE
)

# ============================================================================
# Step 2: Write ancestry-specific GDS files (lantern::write_ancestry_gds())
# ============================================================================
gds_paths <- write_ancestry_gds(split_result, opt$out_path, verbose = TRUE)

meta_path <- file.path(opt$out_path, paste0("split_meta_chr", opt$chr_id, ".rds"))
saveRDS(
    list(
        gds_paths       = gds_paths,
        variant_info    = split_result$variant_info,
        ancestry_counts = split_result$ancestry_counts,
        sample_ids      = split_result$sample_ids,
        mode            = split_result$mode,
        shrink_target   = if (opt$mode == "dosage") "1/2" else NA_character_,
        lantern_version = as.character(utils::packageVersion("lantern")),
        chr_id          = opt$chr_id
    ),
    meta_path
)

message("Wrote: ", paste(unlist(gds_paths), collapse = ", "))
message("Wrote: ", meta_path,
        " (gds_paths + variant_info + ancestry_counts, feeds step2)")
