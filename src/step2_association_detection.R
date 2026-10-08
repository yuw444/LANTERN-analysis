library(optparse)
option_list <- list(
  make_option(
    c("--split_meta"),
    type = "character",
    default = NULL,
    help = "Path to split_meta_chr*.rds from step1 (supplies gds_paths for every population, plus per-variant ancestry counts used to weight the Cauchy combination)",
    metavar = "file"
  ),
  make_option(
    c("--data_file"),
    type = "character",
    default = NULL,
    help = "Path to phenotype file (csv, tsv, or rds). Must contain an 'id' column; column 2 is the response, remaining columns are covariates.",
    metavar = "file"
  ),
  make_option(
    c("--gene_group_file"),
    type = "character",
    default = NULL,
    help = "Path to gene group file (no header: gene, chr, pos, ref, alt, weight)",
    metavar = "file"
  ),
  make_option(
    c("--response_type"),
    type = "character",
    default = "continuous",
    help = "Response type: continuous, binary, or count [default: %default]",
    metavar = "type"
  ),
  make_option(
    c("--kinship_rds"),
    type = "character",
    default = NULL,
    help = "Path to kinship RDS file (square matrix, row/col names = sample IDs)",
    metavar = "file"
  ),
  make_option(
    c("--out_file"),
    type = "character",
    default = NULL,
    help = "Output TSV path for the per-gene results table (gene, p_<POP>, w_<POP>, p_cauchy columns). A sibling RDS with the same basename (extension swapped to .rds) is also written alongside it, containing ancestry_smmat()'s full return value (list(results, smmat_results)) for anyone who needs the raw per-population GMMAT::SMMAT() output too.",
    metavar = "file"
  ),
  make_option(
    c("--ncores"),
    type = "integer",
    default = NA_integer_,
    help = "Cores for GMMAT::SMMAT() (passed through as ancestry_smmat()'s ncores). Defaults to $SLURM_CPUS_PER_TASK when running under SLURM, else 1.",
    metavar = "n"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("split_meta", "data_file", "gene_group_file", "kinship_rds", "out_file")
missing_opts <- required[vapply(opt[required], is.null, logical(1))]
if (length(missing_opts) > 0) {
  stop(
    "Missing required option(s): --", paste(missing_opts, collapse = ", --"),
    "\nSee: Rscript src/step2_association_detection.R --help"
  )
}

suppressPackageStartupMessages({
  library(lantern)
  library(data.table)
})

# read phenotype data (supports .rds or delimited text)
if (grepl("\\.rds$", opt$data_file, ignore.case = TRUE)) {
  df_pheno <- readRDS(opt$data_file)
} else {
  df_pheno <- as.data.frame(fread(opt$data_file))
}
if (!"id" %in% colnames(df_pheno)) {
  stop("Phenotype data must contain an 'id' column.")
}

name_res <- colnames(df_pheno)[2]
name_covs <- colnames(df_pheno)[3:ncol(df_pheno)]
formula_to_fit <- as.formula(paste0(name_res, " ~ ", paste(name_covs, collapse = " + ")))

family_to_use <- switch(
  opt$response_type,
  continuous = gaussian(link = "identity"),
  binary = binomial(link = "logit"),
  count = poisson(link = "log"),
  stop("Unsupported --response_type: ", opt$response_type,
       " (use continuous, binary, or count)")
)

kinship <- readRDS(opt$kinship_rds)

# split_meta (from step1) supplies gds_paths for every population, plus
# per-variant ancestry counts, enabling ancestry_smmat() to weight the
# Cauchy combination per gene.
meta <- readRDS(opt$split_meta)
if (is.null(meta$gds_paths))
  stop("--split_meta has no $gds_paths; regenerate it with the current step1_vcf_split_by_ancestry.R")

ncores <- opt$ncores
if (is.na(ncores)) {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
  ncores <- if (nzchar(slurm_cpus)) as.integer(slurm_cpus) else 1L
}
message("Using ncores = ", ncores, " for GMMAT::SMMAT()")

result <- ancestry_smmat(
  gds_paths       = meta$gds_paths,
  pheno           = df_pheno,
  formula         = formula_to_fit,
  kinship         = kinship,
  gene_group_file = opt$gene_group_file,
  ancestry_counts = meta$ancestry_counts,
  variant_info    = meta$variant_info,
  family          = family_to_use,
  ncores          = ncores
)

dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
fwrite(result$results, opt$out_file, sep = "\t")
message(
  "Wrote: ", opt$out_file,
  " (per-gene p_<POP>/w_<POP>/p_cauchy table, tab-separated)"
)

rds_file <- paste0(tools::file_path_sans_ext(opt$out_file), ".rds")
saveRDS(result, file = rds_file)
message(
  "Wrote: ", rds_file,
  " (list(results = <same table as the TSV above>,",
  " smmat_results = <raw GMMAT::SMMAT() output per population>))"
)
