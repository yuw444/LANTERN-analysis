# LANTERN

R package for **Ancestry-Specific Rare Variant Association Analysis** with pure C backend.

## Features

- **Pure C backend** for performance-critical operations
- Efficient matrix operations for ancestry code counting
- Genotype splitting by ancestry (African/European)
- VCF manipulation utilities

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

# Create test matrices
gt <- matrix(c(1, 2, 0, 1, 0, 2), nrow = 2, ncol = 3)
ancestry <- matrix(c(3, 1, 2, 2, 1, 3), nrow = 2, ncol = 3)

# Count ancestry codes per sample
count_ancestry_codes(ancestry, 3)  # African (code 3)

# Split genotypes by ancestry
result <- split_by_ancestry(gt, ancestry)
result$african   # African-specific dosage
result$european  # European-specific dosage
```

## Core Functions

| Function | Description |
|----------|-------------|
| `count_ancestry_codes(mat, code)` | Count ancestry codes per row |
| `split_by_ancestry(gt, ancestry)` | Split genotype matrix by ancestry |
| `read_bed_file(bed, bim, fam)` | Read PLINK files |
| `write_vcf_with_ancestry(...)` | Write ancestry-specific VCFs |
| `subset_vcf_by_range(...)` | Filter VCF by genomic region |

## Ancestry Codes

| Code | Meaning |
|------|---------|
| 1 | EUR/EUR (Pure European) |
| 2 | AFR/EUR (Mixed) |
| 3 | AFR/AFR (Pure African) |

## Testing

```bash
# Download test data
./download_1000g.sh

# Prepare data
Rscript prepare_1000g.R

# Run tests
Rscript -e 'devtools::test()'

# Run integration test
Rscript test_lantern_full.R
```

## Dependencies

- R >= 4.0
- data.table
- dplyr
- SeqArray (for GDS operations)
- SeqVarTools

## License

MIT
