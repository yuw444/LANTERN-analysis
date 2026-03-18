# LANTERN

**L**ocal **A**ncestry-**N**oted **T**ransmission **R**egion **N**ucleotides

R package for **Ancestry-Specific Rare Variant Association Analysis** with pure C backend.

## Features

- **Pure C backend** for performance-critical operations
- Efficient matrix operations for ancestry code counting
- Genotype splitting by ancestry (African/European)
- Direct data frame/matrix input (no PLINK dependency)

## Installation

```r
# Install from GitHub
devtools::install_github("yuw444/LANTERN", subdir = "lantern")

# Or from local source
R CMD INSTALL lantern/
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

## Core Functions

| Function | Description |
|----------|-------------|
| `count_ancestry_codes(mat, code)` | Count ancestry codes per row |
| `split_by_ancestry(gt, pt)` | Split genotype matrix by ancestry |
| `run_ancestry_pipeline(gt, pt)` | Full pipeline with counts |
| `write_vcf_with_ancestry(...)` | Write ancestry-specific VCFs |
| `subset_vcf_by_range(...)` | Filter VCF by genomic region |

## Dependencies

- R >= 4.0
- data.table
- dplyr

## License

MIT
