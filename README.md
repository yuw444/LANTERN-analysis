# LANTERN 

**L**everaging Local **AN**cestry **T**racts to **E**nhance **R**are-Varia**N**t Aggregate Association Testing

[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://yuw444.github.io/LANTERN/)
[![GitHub](https://img.shields.io/badge/source-GitHub-lightgrey)](https://github.com/yuw444/LANTERN)

Full documentation, vignettes, and function reference: **[https://yuw444.github.io/LANTERN/](https://yuw444.github.io/LANTERN/)**. MedRxiv: [2026.04. 24.26351693](https://www.medrxiv.org/content/10.64898/2026.04.24.26351693v1.full.pdf)

## Features

- **Pure C backend** for performance-critical operations
- Efficient matrix operations for ancestry code counting
- Genotype splitting by local ancestry (African/European)
- Supporting multiple mixed ancestry (Up to 5)
- Direct data frame/matrix input (no PLINK dependency)
- **Automatic overlap handling** for sample and variant mismatches
- **Monomorphic filtering** - removes variants with no alt alleles

## Installation

SeqArray is a Bioconductor package and is not found by
`devtools::install_github()` automatically. Use BiocManager to handle
dependencies in one step:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("yuw444/LANTERN/lantern")
```

Or pre-install the Bioconductor dependency first, then use devtools:

```r
BiocManager::install("SeqArray")
devtools::install_github("yuw444/LANTERN", subdir = "lantern")
```

## Quick Start

LANTERN supports two perspectives on ancestry splitting, depending on what
input data you have:

- **Dosage (proportional) split** — unphased genotype dosages (0/1/2) +
  population-level parent-of-origin ancestry codes. Heterozygous mixed-ancestry
  genotypes are split by an estimated population allele proportion. Entry
  point: `ancestry_split_dosage()` / `split_diploid()`.
- **Phased split** — phased haplotypes (`0|1`-style genotypes) + per-haplotype
  local ancestry tracts (e.g. from RFMix). Each haplotype's allele is
  deterministically assigned to the ancestry pool it was inferred to
  originate from — no proportional estimation needed. Entry point:
  `ancestry_split_phased()` / `split_haplotype()`.

```r
library(lantern)
```

### Dosage split

```r
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
result <- ancestry_split_dosage(gt, pt)
result$african    # African ancestry-specific dosages
result$european   # European ancestry-specific dosages
result$counts     # Ancestry counts per region
```

### Phased split

```r
# Haplotype allele matrices: 2 variants x 2 samples, alleles (0/1)
gt_hap0 <- matrix(c(1L, 0L, 0L, 1L), nrow = 2, ncol = 2)
gt_hap1 <- matrix(c(0L, 1L, 1L, 0L), nrow = 2, ncol = 2)

# Local ancestry per haplotype (RFMix convention: AFR=0, EUR=1)
anc_hap0 <- matrix(c(0L, 1L, 0L, 1L), nrow = 2, ncol = 2)
anc_hap1 <- matrix(c(1L, 0L, 1L, 0L), nrow = 2, ncol = 2)

result <- split_haplotype(gt_hap0, gt_hap1, anc_hap0, anc_hap1)
result$african    # African ancestry-specific dosages
result$european   # European ancestry-specific dosages
```

Or run the full pipeline directly from a phased VCF/BCF + RFMix MSP file:

```r
result <- ancestry_split_phased(
  vcf_path = "data/chr19.phased.bcf",
  msp_path = "data/chr19.msp.tsv.gz",
  out_path = "output/"
)
```

## Input Format

### Dosage split

#### PT (Parent-of-Origin) Matrix

| sample_id | chr1:1000-2000 | chr1:2000-3000 | chr2:5000-6000 |
| --------- | :------------: | :------------: | :------------: |
| S1        |       3       |       2       |       3       |
| S2        |       2       |       2       |       3       |
| S3        |       1       |       1       |       2       |

- **Rows**: Samples (must have sample_id column or rownames)
- **Columns**: Genomic regions/windows
- **Values**: Ancestry codes
  - `1` = EUR/EUR (Pure European)
  - `2` = AFR/EUR (Mixed)
  - `3` = AFR/AFR (Pure African)

#### GT (Genotype) Matrix

|           | S1 | S2 | S3 |
| --------- | :-: | :-: | :-: |
| chr1:1234 | 0 | 2 | 0 |
| chr1:2345 | 1 | 1 | 0 |

- **Rows**: Variants
- **Columns**: Samples (must match PT matrix columns)
- **Values**: 0, 1, 2 (dosage of alternate allele)

### Phased split

#### Ancestry tract matrices (`anc_hap0` / `anc_hap1`)

Produced by parsing an RFMix `.msp` file, or supplied directly.

- **Rows**: Variants (tract calls broadcast to each variant they cover)
- **Columns**: Samples (must match haplotype matrix columns)
- **Values**: Population code, per `pop_codes` (RFMix default: `AFR = 0`, `EUR = 1`) — **not** the same 1/2/3 diploid codes used by the dosage split's PT matrix

##### Where the MSP file comes from

The `.msp.tsv[.gz]` file is one of several outputs written by
[RFMix2](https://github.com/slowkoni/rfmix), a local-ancestry inference tool.
A typical run that produces it looks like:

```bash
rfmix -f query.phased.vcf.gz \
      -r reference_panel.phased.vcf.gz \
      -m sample_map.tsv \
      -g genetic_map.txt \
      -o LANTERN_chr19 \
      --chromosome=chr19
```

This writes several sibling files sharing the `LANTERN_chr19` prefix; LANTERN
only reads the `.msp.tsv.gz`:

| File             | Content                                                                           | Used by LANTERN?                                                                                |
| ---------------- | --------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `*.msp.tsv.gz` | Viterbi (most-likely) local-ancestry call per haplotype per tract                 | **Yes** — `.parse_msp()` / `ancestry_split_phased()` / `ancestry_split_combined()` |
| `*.fb.tsv.gz`  | Forward-backward posterior probability per population per haplotype per marker    | No                                                                                              |
| `*.rfmix.Q.gz` | Global (genome-wide) ancestry proportion per sample, ADMIXTURE-style`.Q` format | No                                                                                              |
| `*.sis.tsv.gz` | Per-SNP interpolated ancestry probability                                         | No                                                                                              |

##### MSP file structure

```
#Subpopulation order/codes: AFR=0	EUR=1
#chm	spos	epos	sgpos	egpos	n snps	SAMPLE1.0	SAMPLE1.1	SAMPLE2.0	SAMPLE2.1	...
chr19	226776	518686	0.00	1.17	457	0	0	0	1	...
chr19	518686	554919	1.17	1.35	85	0	0	0	1	...
```

- **Line 1** (`#`-prefixed comment): population name → code mapping, e.g.
  `AFR=0  EUR=1`. Parsed into `pop_codes`.
- **Line 2** (`#`-prefixed comment): column headers.
- **Data rows**: one row per contiguous local-ancestry **tract** — a genomic
  segment RFMix2 called as a single ancestry block, *not* one row per
  variant:
  - `chm` — chromosome
  - `spos` / `epos` — tract start/end, physical position (bp)
  - `sgpos` / `egpos` — tract start/end, genetic position (cM)
  - `n snps` — number of markers RFMix2 used to call this tract
  - remaining columns — one per **haplotype** (named `<sample_id>.0`,
    `<sample_id>.1`), holding the ancestry code from line 1 for that
    haplotype over that tract

`.parse_msp()` broadcasts each tract's call out to every variant whose
position falls within `[spos, epos]`, producing the `anc_hap0`/`anc_hap1`
matrices at variant resolution.

#### Haplotype genotype matrices (`gt_hap0` / `gt_hap1`)

Produced by splitting a phased VCF's `0|1`-style genotype calls into one
single-allele matrix per haplotype, or supplied directly.

- **Rows**: Variants
- **Columns**: Samples (must match ancestry tract matrix columns)
- **Values**: 0 or 1 (the allele carried on that haplotype)

## Automatic Overlap Handling

`ancestry_split_dosage()` automatically handles mismatches between GT and PT
matrices (examples below use the dosage split; `ancestry_split_phased()` and
`ancestry_split_combined()` perform the equivalent sample intersection
between the VCF and MSP file automatically):

### Sample Mismatches

```r
# GT has samples A, B, C
# PT has samples A, B, D
# -> Only A and B are used

gt <- matrix(0, nrow = 2, ncol = 3,
             dimnames = list(c("v1", "v2"), c("A", "B", "C")))
pt <- matrix(1, nrow = 3, ncol = 2,
             dimnames = list(c("A", "B", "D"), c("v1", "v2")))

result <- ancestry_split_dosage(gt, pt)
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

result <- ancestry_split_dosage(gt, pt)
# result$overlap$n_variants_kept = 2
```

## Core Functions

### High-level pipelines

Parse VCF/MSP files (or accept pre-built matrices) and run a full ancestry split.

| Function                                                     | Description                                                                                                                                                                                                    |
| ------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ancestry_split_dosage(gt, pt, ...)`                       | Full pipeline: proportional dosage split of a genotype matrix by parent-of-origin ancestry, with automatic sample/variant overlap handling. Accepts matrices directly, or a`vcf_path`/`msp_path` shortcut. |
| `ancestry_split_phased(vcf_path, msp_path, out_path, ...)` | Full pipeline for phased data: parse a phased VCF + RFMix MSP file, deterministically split each haplotype by its local ancestry tract, and optionally write ancestry-specific VCFs/GDS.                       |
| `ancestry_split_combined(vcf_path, msp_path, ...)`         | Parses the VCF/BCF + MSP file once and runs both the phased and proportional splits together, for any number of populations K ≥ 2.                                                                            |

### Low-level splitters

Matrix-in, matrix-out primitives (C backend), used internally by the pipelines above but also usable directly.

| Function                                                                   | Description                                                                                                                       |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `split_diploid(gt, ancestry)`                                            | Split a genotype matrix into African/European dosage matrices using parent-of-origin ancestry codes (2-population, proportional). |
| `split_diploid_multi(gt, ancestry, pure_codes, mixed_codes)`             | K-population generalisation of`split_diploid`.                                                                                  |
| `split_haplotype(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes)`       | Deterministic per-haplotype ancestry split (2-population) from phased genotype + ancestry-tract matrices.                         |
| `split_haplotype_multi(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes)` | K-population generalisation of`split_haplotype`.                                                                                |

### I/O and utilities

| Function                                                             | Description                                                                                |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `count_ancestry_codes(mat, code)`                                  | Count occurrences of an ancestry code in each row of a PT matrix.                          |
| `write_dosage_gds(dosage_mat, variant_info, sample_ids, gds_path)` | Convert an ancestry-specific dosage matrix to a SeqArray GDS file (for`GMMAT::SMMAT()`). |

### Statistics

| Function                                     | Description                                                                                             |
| -------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `cauchy_combine(p_values, weights = NULL)` | Cauchy combination test — merges K (possibly correlated) p-values into a single meta-analysis p-value. |

## Splitting Algorithms

The two perspectives differ in one fundamental way: **the dosage split has to
estimate**, because an unphased genotype only tells you *how many* alt alleles
a sample carries, not *which parental haplotype* each one sits on. The
**phased split never has to estimate**, because phasing plus a local-ancestry
call already answers that question directly, allele by allele.

### Dosage (Proportional) Split

For heterozygous genotypes (gt=1) with mixed ancestry (pt=2), the algorithm uses
population-based proportions (p1, p2):

#### Formulas

```
p1 = (2*N1 + N2 + N4) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
p2 = (N4 + 2*N7 + N8) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
```

Where N1-N8 are counts per variant, over PT and GT matrix entries for the same sample:

| Code | PT, GT     | Description                   |
| ---- | ---------- | ----------------------------- |
| N1   | pt=3, gt=2 | Pure African, homozygous alt  |
| N2   | pt=3, gt=1 | Pure African, heterozygous    |
| N4   | pt=2, gt=2 | Mixed, homozygous alt         |
| N5   | pt=2, gt=1 | Mixed, heterozygous           |
| N7   | pt=1, gt=2 | Pure European, homozygous alt |
| N8   | pt=1, gt=1 | Pure European, heterozygous   |

#### Special Cases

- **Homozygous alt (gt=2)**: 1 allele to each ancestry regardless of pt
- **Pure ancestry (pt=1 or pt=3)**: All alt alleles to that ancestry
- **N5 dominates**: When ambiguous mixed-ancestry heterozygotes (N5) make up
  most or all of a variant's carriers — up to and including the singleton
  case, N5 = every carrier, where the formula is literally undefined
  (denominator = 0) — the raw p1/p2 ratio becomes unreliable. `ancestry_split()`
  applies **GLA shrinkage** by default to fix this; see below.

#### GLA shrinkage: why an ambiguous-dominated variant needs outside evidence

The raw formula estimates a variant's ancestry split from its own
unambiguous carriers, then applies that ratio to the ambiguous ones. That
works fine when unambiguous carriers are plentiful — but consider N5 = 10
(ten ambiguous mixed-ancestry heterozygotes) with everything else zero
except N7 = 1 (one pure-European homozygous-alt carrier):

```
total_alt   = 2*N1 + N2 + 2*N4 + N5 + 2*N7 + N8 = 2*0+0+2*0+10+2*1+0 = 12
denominator = total_alt - N5                    = 12 - 10             = 2
p1 = (2*N1 + N2 + N4) / denominator = (0+0+0) / 2 = 0
p2 = (N4 + 2*N7 + N8) / denominator = (0+2+0) / 2 = 1
```

The formula confidently assigns **all ten** ambiguous carriers 100% European
ancestry (p1=0, p2=1), extrapolated from a single unambiguous data point.
Unlike the pure singleton, this isn't a "no answer" case — the denominator is
nonzero, so nothing flags it as unreliable. It's a *confidently
wrong-looking* answer: the raw formula can't distinguish "1 unambiguous
carrier informing 2 total alt alleles" from "1 unambiguous carrier informing
100 total" — it takes whatever ratio it computes at face value, regardless of
how little evidence backs it.

**GLA (global local ancestry) shrinkage** fixes this — and the pure
singleton — with one mechanism: discount the raw ratio in proportion to how
much of a variant's evidence is actually ambiguous, and fill in the rest from
the chromosome arm's independently-estimated ancestry mixture (from RFMix's
tracts, not this variant's genotypes):

```
w  = N5 / (N1 + N2 + N4 + N5 + N7 + N8)   # fraction of carriers that are ambiguous
p1 = (1 - w) * p1_raw + w * GLA_AFR[arm]
```

For the example above, `w = 10/11 ≈ 0.91` — the near-total ambiguity is
almost fully discounted. If this cohort's arm-level ancestry runs 60% AFR /
40% EUR, GLA shrinkage gives `p1 ≈ 0.55, p2 ≈ 0.45` instead of the raw
formula's `p1=0, p2=1` — a moderate, defensible estimate instead of a hard
call backed by one data point. `w = 0` (unambiguous carriers dominate) leaves
the raw formula untouched; `w = 1` (the pure singleton) falls back entirely
to the arm's GLA proportion instead of a coin flip — the old special case is
just one end of this same continuum, not a separate mechanism.

`ancestry_split(mode = "dosage")` computes and applies this automatically.
Pass `use_gla = FALSE` to disable it and reproduce the original raw-ratio /
flat-0.5-singleton estimator exactly (or call
`split_diploid()`/`split_diploid_multi()` directly with `gla = NULL`, the
default, for the same effect at the matrix level). See
`vignette("split-intuition")` for a worked example and the K > 2
generalisation.

#### Intuition for K > 2 populations (`split_diploid_multi`)

With more than two source populations, each mixed ancestry code names an
*unordered pair* of parent populations (e.g. code 5 = AFR/NAT). The algorithm
generalises p1/p2 to K populations in two passes per variant:

1. **Build a per-population weight from unambiguous carriers.** A pure
   homozygous-alt sample contributes 2 alleles to its own population's
   weight; a pure heterozygote contributes 1; and — the key unambiguous case —
   a **mixed-pair homozygous-alt** sample contributes exactly 1 allele to
   *each* of its two parent populations, because carrying two alt copies with
   ancestry split between two populations means one copy must have come from
   each side.
2. **Distribute the ambiguous alleles proportionally.** A mixed-pair
   heterozygote carries exactly one alt allele but its parental origin is
   unknown, so it's split between its two parent populations `i` and `j`
   using the conditional ratio `w_i / (w_i + w_j)` — i.e. in proportion to how
   often each population's allele shows up unambiguously elsewhere at this
   variant. If neither population has unambiguous evidence at this variant
   (`w_i + w_j == 0`, the singleton case), GLA shrinkage (see above) supplies
   the arm's conditional ancestry ratio for that pair instead of a flat
   0.5/0.5 — same idea as the two-population case, generalised per pair.

With K=2 there's only one mixed code and one pair, so this reduces exactly to
the p1/p2 formulas above.

### Phased Split

Deterministic, per-haplotype, no estimation:

```
for each sample, each variant, each haplotype h in {hap0, hap1}:
    pop <- lookup(anc_hap_h[variant, sample], pop_codes)
    if pop is known:
        dosage[pop][variant, sample] += gt_hap_h[variant, sample]
```

Each haplotype carries exactly one allele (0 or 1), and its local-ancestry
tract call says exactly which population pool that allele belongs to — so a
mixed-ancestry individual's het genotype is never actually ambiguous once
phased: hap0's allele goes to whichever population hap0's tract says, hap1's
allele goes to whichever population hap1's tract says. Haplotypes whose
ancestry code isn't in `pop_codes` (e.g. RFMix's "unassigned" call) contribute
0 to every population.

#### Intuition for K > 2 populations (`split_haplotype_multi`)

Nothing changes conceptually — `pop_codes` just grows from 2 entries to K, and
the per-haplotype lookup routes each allele into one of K dosage matrices
instead of 2. There's no proportional step to generalise, because the whole
reason the dosage split needs proportions (unresolved parental origin within
a genotype) doesn't exist once haplotypes are individually ancestry-labeled.

## Association Testing

`ancestry_smmat()` (Step 3) runs `GMMAT::SMMAT()` once per population's GDS
file (plus, typically, once more on the unsplit "observed" GT GDS), then
combines the per-population p-values into a single gene-level result with
`cauchy_combine()`.

### Why per-population p-values aren't combined equally

Not every population contributes equally reliable evidence for every gene.
A gene sitting in a genomic region where only a handful of cohort samples
have *pure* ancestry $k$ gives that population's SMMAT p-value little to
work with, and shouldn't count as much as a population with abundant
pure-ancestry representation there. So instead of feeding `cauchy_combine()`
equal weights, `ancestry_smmat()` weights each population's p-value by how
much pure-ancestry evidence backs it *at that specific gene*:

```
w_k(gene) = median over the gene's variants of:
              (count of cohort samples with pure ancestry k at that variant)
```

using the **median** across the gene's variants — from `ancestry_split()`'s
`ancestry_counts` output — rather than, say, the mean, so a handful of
unusually ancestry-rich or ancestry-poor variants at the gene's edges don't
dominate the estimate. A population with little pure-ancestry
representation at this gene is discounted, not ignored outright — a gene
well-represented in one ancestry but barely in another naturally leans
toward the better-supported population's signal in the combined p-value.
If none of the populations' names match `ancestry_counts`'s columns (e.g. it
wasn't supplied), `ancestry_smmat()` falls back to equal weighting across
all K populations instead.

This weight computation happens automatically inside `ancestry_smmat()` —
there's no separate weight-finding step to run. See `?ancestry_smmat` for
the full argument/return-value reference.

## Dependencies

### R packages

| Package                                                      | Type      | Source       | Notes                  |
| ------------------------------------------------------------ | --------- | ------------ | ---------------------- |
| [data.table](https://cran.r-project.org/package=data.table)   | Required  | CRAN         | —                     |
| [SeqArray](https://bioconductor.org/packages/SeqArray/)       | Required  | Bioconductor | GDS file I/O           |
| [GMMAT](https://cran.r-project.org/package=GMMAT)             | Suggested | CRAN         | SMMAT gene-level tests |
| [dplyr](https://cran.r-project.org/package=dplyr)             | Suggested | CRAN         | Vignettes only         |

**Bioconductor packages (SeqArray) are not installed automatically by `devtools::install_github()`.**
Use the `BiocManager` installation instructions above.

### System tools

| Tool                                                | Version | Required for                                                                                               |
| --------------------------------------------------- | ------- | ---------------------------------------------------------------------------------------------------------- |
| [`bcftools`](https://samtools.github.io/bcftools/) | ≥ 1.10 | VCF/BCF-reading functions only (`ancestry_split()`, `ancestry_split_phased()`) — not needed for the matrix-in/matrix-out functions (`split_diploid()`, `ancestry_split_dosage()`, etc.) |

Install bcftools via your system package manager:

```bash
# macOS
brew install bcftools

# Debian / Ubuntu
sudo apt install bcftools

# Conda / pixi
conda install -c bioconda bcftools
```

The unphased splitting functions (`split_diploid`, `split_diploid_multi`,
`ancestry_split_dosage`) do **not** require bcftools.

### Compiler

A C compiler (gcc or clang) is required to build the package from source. No
external C libraries are needed — the C backend uses only standard R headers
and plain C file I/O.

## License

MIT
