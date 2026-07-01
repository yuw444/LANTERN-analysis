# LANTERN

**L**everaging Local **AN**cestry **T**racts to **E**nhance **R**are-Varia**N**t Aggregate Association Testing

[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://yuw444.github.io/LANTERN/)
[![GitHub](https://img.shields.io/badge/source-GitHub-lightgrey)](https://github.com/yuw444/LANTERN)

Full documentation, vignettes, and function reference: **<https://yuw444.github.io/LANTERN/>**

## Features

- **Pure C backend** for performance-critical operations
- Efficient matrix operations for ancestry code counting
- Genotype splitting by local ancestry (African/European)
- Supporting multiple mixed ancestry (Up to 5)
- Direct data frame/matrix input (no PLINK dependency)
- **Automatic overlap handling** for sample and variant mismatches
- **Monomorphic filtering** - removes variants with no alt alleles

## Installation

SeqArray and SeqVarTools are Bioconductor packages and are not found by
`devtools::install_github()` automatically. Use BiocManager to handle all
dependencies in one step:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("yuw444/LANTERN/lantern")
```

Or pre-install the Bioconductor packages first, then use devtools:

```r
BiocManager::install(c("SeqArray", "SeqVarTools"))
devtools::install_github("yuw444/LANTERN", subdir = "lantern")
```

## Quick Start

```r
library(lantern)

# GT matrix: 5 variants x 4 samples
# 0 = homozygous ref, 1 = heterozygous, 2 = homozygous alt
gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1, 0,
               1, 1, 1, 0, 2, 1, 0, 1, 0, 1), 
             nrow = 5, ncol = 4)

# PT matrix: parent-of-origin ancestry codes (rows=samples, cols=variants)
# 1 = EUR/EUR, 2 = AFR/EUR (mixed), 3 = AFR/AFR
pt <- matrix(c(3, 2, 1, 3, 2, 1, 2, 2, 1, 1,
               3, 1, 2, 1, 3, 2, 1, 2, 1, 3), 
             nrow = 4, ncol = 5)

# Run full pipeline
result <- run_ancestry_pipeline(gt, pt)
result$african    # African ancestry-specific dosages
result$european   # European ancestry-specific dosages
result$counts     # Ancestry counts per region
```

## Input Format

### PT (Parent-of-Origin) Matrix

| sample_id | chr1:1000-2000 | chr1:2000-3000 | chr2:5000-6000 |
|-----------|:--------------:|:--------------:|:--------------:|
| S1        | 03             | 02             | 03             |
| S2        | 02             | 02             | 03             |
| S3        | 01             | 01             | 02             |

- **Rows**: Samples (must have sample_id column or rownames)
- **Columns**: Genomic regions/windows
- **Values**: Ancestry codes
  - `1` = EUR/EUR (Pure European)
  - `2` = AFR/EUR (Mixed)
  - `3` = AFR/AFR (Pure African)

### GT (Genotype) Matrix

|            | S1 | S2 | S3 |
|------------|:--:|:--:|:--:|
| chr1:1234  | 0/1| 1/1| 0/0|
| chr1:2345  | 1/0| 0/1| 0/0|

- **Rows**: Variants
- **Columns**: Samples (must match PT matrix columns)
- **Values**: 0, 1, 2 (dosage of alternate allele)

## Automatic Overlap Handling

The pipeline automatically handles mismatches between GT and PT matrices:

### Sample Mismatches
```r
# GT has samples A, B, C
# PT has samples A, B, D
# -> Only A and B are used

gt <- matrix(0, nrow = 2, ncol = 3,
             dimnames = list(c("v1", "v2"), c("A", "B", "C")))
pt <- matrix(1, nrow = 3, ncol = 2,
             dimnames = list(c("A", "B", "D"), c("v1", "v2")))

result <- run_ancestry_pipeline(gt, pt)
# result$overlap$n_samples_kept = 2
# result$overlap$dropped_samples = c("C", "D")
```

### Variant/Region Mismatches
```r
# GT variants: chr22:100, chr22:200, chr22:300
# PT regions: chr22:50-150, chr22:150-250
# -> Only chr22:100 and chr22:200 are used (matched by coordinate)

gt <- matrix(0, nrow = 3, ncol = 2,
             dimnames = list(c("chr22:100", "chr22:200", "chr22:300"),
                             c("s1", "s2")))
pt <- matrix(1, nrow = 2, ncol = 2,
             dimnames = list(c("s1", "s2"),
                             c("chr22:50-150", "chr22:150-250")))

result <- run_ancestry_pipeline(gt, pt)
# result$overlap$n_variants_kept = 2
```

## Core Functions

| Function | Description |
|----------|-------------|
| `count_ancestry_codes(mat, code)` | Count ancestry codes per row |
| `split_by_ancestry(gt, pt)` | Split genotype matrix by ancestry |
| `run_ancestry_pipeline(gt, pt)` | Full pipeline with counts |
| `write_vcf_with_ancestry(...)` | Write ancestry-specific VCFs |
| `subset_vcf_by_range(...)` | Filter VCF by genomic region |

## Ancestry Splitting Algorithm

For heterozygous genotypes (gt=1) with mixed ancestry (pt=2), the algorithm uses
population-based proportions (p1, p2):

### Formulas

```
p1 = (2*N1 + N2 + N4) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
p2 = (N4 + 2*N7 + N8) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
```

Where N1-N8 are counts per variant:

| Code | GT | Description |
|------|-----|-------------|
| N1 | pt=3, gt=2 | Pure African, homozygous alt |
| N2 | pt=3, gt=1 | Pure African, heterozygous |
| N4 | pt=2, gt=2 | Mixed, homozygous alt |
| N5 | pt=2, gt=1 | Mixed, heterozygous |
| N7 | pt=1, gt=2 | Pure European, homozygous alt |
| N8 | pt=1, gt=1 | Pure European, heterozygous |

### Special Cases

- **Singleton**: When all alt alleles come from mixed het individuals (N5 = sum(gt)),
  p1 = p2 = 0.5
- **Homozygous alt (gt=2)**: 1 allele to each ancestry regardless of pt
- **Pure ancestry (pt=1 or pt=3)**: All alt alleles to that ancestry

## Dependencies

### R packages

| Package | Type | Source | Notes |
|---------|------|--------|-------|
| [data.table](https://cran.r-project.org/package=data.table) | Required | CRAN | — |
| [SeqArray](https://bioconductor.org/packages/SeqArray/) | Required | Bioconductor | GDS file I/O |
| [SeqVarTools](https://bioconductor.org/packages/SeqVarTools/) | Required | Bioconductor | Variant iteration |
| [GMMAT](https://cran.r-project.org/package=GMMAT) | Suggested | CRAN | SMMAT gene-level tests |
| [dplyr](https://cran.r-project.org/package=dplyr) | Suggested | CRAN | Vignettes only |

**Bioconductor packages are not installed automatically by `devtools::install_github()`.**
Use the `BiocManager` installation instructions above.

### System tools

| Tool | Version | Required for |
|------|---------|-------------|
| [`bcftools`](https://samtools.github.io/bcftools/) | ≥ 1.10 | Phased-mode functions only (`read_phased_vcf`, `run_phased_pipeline`, `filter_phased_vcf_samples`) |

Install bcftools via your system package manager:

```bash
# macOS
brew install bcftools

# Debian / Ubuntu
sudo apt install bcftools

# Conda / pixi
conda install -c bioconda bcftools
```

The unphased splitting functions (`split_by_ancestry`, `split_by_ancestry_multi`,
`run_ancestry_pipeline`) do **not** require bcftools.

### Compiler

A C compiler (gcc or clang) is required to build the package from source. No
external C libraries are needed — the C backend uses only standard R headers
and plain C file I/O.

## License

MIT
