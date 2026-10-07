# LANTERN <img src="lantern/man/figures/logo.png" align="right" height="139" alt="LANTERN logo" />

**L**everaging Local **AN**cestry **T**racts to **E**nhance **R**are-Varia**N**t Aggregate Association Testing

## 1. Background

* LANTERN is a method for conducting aggregate rare-variant association tests using inferred local ancestry as additional information. 
* For simplicity, in the following we demonstrate LANTERN on two ancestries, **african(AFR)** and **european(EUR)**.
* LANTERN can be implemented on up to 5 ancestries. 



## 2. Data

### Ancestry Matrix

| Genome track        | Subject_S1(AN) | Subject_S2(AN) | Subject_S3(AN) |
|---------------------|:------------:|:------------:|:------------:|
| chr1:1,000-2,000    | 03         | 02         | 01         |
| chr1:2,000-3,000    | 02         | 02         | 01         |
| chr2:5,000-6,000    | 03         | 03         | 02         |

Explanation:
- Track: Half-open [start, end) (end excluded)
- Rows: contiguous genome tracks or windows used in analysis
- Columns: individual subjects in the cohort.
- Entries: Inferred ancestry (03 = 2 AFR chromosomes, 02 = 1 AFR and 1 EUR chromosome , 01 = 2 EUR chromosomes).

### Variant Matrix

| Variant (CHROM:POS REF>ALT) | Subject_S1 (GT) | Subject_S2 (GT) | Subject_S3 (GT) |
|-----------------------------|:-----------------:|:-----------------:|:-----------------:|
| chr1:1,234 A>G              | 0/1             | 1/1             | 0/0             |
| chr1:2,345 T>C              | 1/0             | 0/1             | 0/0             |
| chr2:5,678 G>A              | 0/0             | 0/1             | 1/1             |

Notes:
- Rows: variant loci (chromosome:position, REF>ALT).
- Columns: per-sample genotype (GT) field from VCF.
- GT meanings: 1/1 = homozygous alternate, 0/1 or 1/0 = heterozygous (phase unknown), 0/0 = homozygous reference.

### Variant Matrix with Ancestry Annotation

Combine GT with the ancestry matrix to map alternate alleles to ancestry-of-origin where possible.

| Variant (CHROM:POS REF>ALT) | Subject_S1 (GT:AN) | Subject_S2 (GT:AN) | Subject_S3 (GT:AN) |
|-----------------------------|:--------------------:|:--------------------:|:--------------------:|
| chr1:1,234 A>G              | 0/1:01             | 1/1:02             | 0/0:02             |
| chr1:2,345 T>C              | 1/0:03             | 0/1:03             | 0/0:03             |
| chr2:5,678 G>A              | 0/0:02             | 0/1:01             | 1/1:01             |

## 3. Method

### Ancestry Specific Variant Matrix Split

#### 2-Ancestry Case (AFR + EUR)

All counts are per variant.

| Ancestry         | Genotype         | Freq  | Condition           | $\mathbf{x}_{AFR}$ | $\mathbf{x}_{EUR}$ |
|------------------|------------------|-------|---------------------|:------------------:|:------------------:|
| **AFR/AFR (03)** | 1/1 (2)          | $N_1$ | pt = 3 & gt = 2     | 2                  | 0                  |
|                  | 0/1 (1)          | $N_2$ | pt = 3 & gt = 1     | 1                  | 0                  |
|                  | 0/0 (0)          | $N_3$ | pt = 3 & gt = 0     | 0                  | 0                  |
| **AFR/EUR (02)** | 1/1 (2)          | $N_4$ | pt = 2 & gt = 2     | 1                  | 1                  |
|                  | 0/1 (1)          | $N_5$ | pt = 2 & gt = 1     | $p_1$              | $p_2$              |
|                  | 0/0 (0)          | $N_6$ | pt = 2 & gt = 0     | 0                  | 0                  |
| **EUR/EUR (01)** | 1/1 (2)          | $N_7$ | pt = 1 & gt = 2     | 0                  | 2                  |
|                  | 0/1 (1)          | $N_8$ | pt = 1 & gt = 1     | 0                  | 1                  |
|                  | 0/0 (0)          | $N_9$ | pt = 1 & gt = 0     | 0                  | 0                  |

$$p_1 = \frac{2N_1 + N_2 + N_4}{D}, \qquad p_2 = \frac{N_4 + 2N_7 + N_8}{D}, \qquad D = 2N_1 + N_2 + 2N_4 + 2N_7 + N_8$$

Special case: if the only alt carriers are AFR/EUR heterozygotes ($N_5 \ge 1$, all other $N = 0$), the formula above is undefined ($D = 0$).

**Shrinkage toward 1/2.** The formula above only uses evidence *at this
variant*: it estimates $p_1, p_2$ from unambiguous carriers, then applies
that ratio to the ambiguous AFR/EUR hets ($N_5$). That's unreliable whenever
$N_5$ dominates the carrier count — not only in the pure singleton above
($D=0$), but also just short of it. For example, $N_5 = 10$ with only
$N_7 = 1$ (one pure-EUR hom-alt) and everything else 0 gives $D = 2$,
$p_1 = 0$, $p_2 = 1$ — a fully confident call extrapolated from a single
unambiguous carrier.

`lantern::ancestry_split()` shrinks $p_1$ toward $1/2$, weighted by how much
of the variant's own evidence is ambiguous:

$$w = \frac{N_5}{D+N_5}, \qquad p_1 \leftarrow (1-w) p_1 + \frac{w}{2}, \qquad p_2 = 1 - p_1$$

using $D$ from the raw formula above. Hom-alt carriers ($N_1, N_4, N_7$)
count 2x toward this denominator, matching their 2x weight in $p_1$'s
numerator, since each carries two alt alleles' worth of ancestry evidence
versus one for a het or an ambiguous mixed het ($N_5$, which always carries
exactly one, unresolved, allele). In the example above, $w = 10/12$ and each
ambiguous het gets $p_1 = 5/12$ instead of 0.

$w = 0$ (unambiguous carriers dominate) reduces to the raw formula above;
$w = 1$ (the pure singleton, $D=0$) gives exactly $0.5/0.5$ — the old
special case is just one end of this continuum.

*Why $1/2$?* A mixed-ancestry heterozygote carries exactly one AFR and one
EUR haplotype at the site, so which one holds the allele depends on how
common the allele is on AFR versus EUR haplotypes, not on how many AFR
haplotypes the cohort has. An earlier version shrank toward the chromosome
arm's global local ancestry (GLA, ~0.82 AFR in the chr19 cohort); in chr19
simulations the true AFR share of mixed-het alleles was 0.50 at every allele
count, and the $1/2$ target both assigned alleles more accurately and gave
more power when the causal effect was on EUR alleles. See the `lantern`
package [README](lantern/README.md) and `vignette("split-intuition")` for
worked examples.

---

#### 3-Ancestry Case (AFR + EUR + NAT)

With a third ancestry (NAT = Native/other), there are 3 pure types and 3 mixed pairs.  
All counts are per variant; 0/0 rows (gt = 0) contribute nothing and are omitted for brevity.

| Ancestry         | Genotype | Freq       | $\mathbf{x}_{AFR}$               | $\mathbf{x}_{EUR}$               | $\mathbf{x}_{NAT}$               |
|------------------|----------|------------|:--------------------------------:|:--------------------------------:|:--------------------------------:|
| **AFR/AFR**      | 1/1 (2)  | $N_1$      | 2                                | 0                                | 0                                |
|                  | 0/1 (1)  | $N_2$      | 1                                | 0                                | 0                                |
| **EUR/EUR**      | 1/1 (2)  | $N_3$      | 0                                | 2                                | 0                                |
|                  | 0/1 (1)  | $N_4$      | 0                                | 1                                | 0                                |
| **NAT/NAT**      | 1/1 (2)  | $N_5$      | 0                                | 0                                | 2                                |
|                  | 0/1 (1)  | $N_6$      | 0                                | 0                                | 1                                |
| **AFR/EUR**      | 1/1 (2)  | $N_7$      | 1                                | 1                                | 0                                |
|                  | 0/1 (1)  | $N_8$      | $\dfrac{p_1}{p_1+p_2}$           | $\dfrac{p_2}{p_1+p_2}$           | 0                                |
| **AFR/NAT**      | 1/1 (2)  | $N_9$      | 1                                | 0                                | 1                                |
|                  | 0/1 (1)  | $N_{10}$   | $\dfrac{p_1}{p_1+p_3}$           | 0                                | $\dfrac{p_3}{p_1+p_3}$           |
| **EUR/NAT**      | 1/1 (2)  | $N_{11}$   | 0                                | 1                                | 1                                |
|                  | 0/1 (1)  | $N_{12}$   | 0                                | $\dfrac{p_2}{p_2+p_3}$           | $\dfrac{p_3}{p_2+p_3}$           |

$$D = 2N_1 + N_2 + 2N_3 + N_4 + 2N_5 + N_6 + 2N_7 + 2N_9 + 2N_{11}$$

$$p_1 = \frac{2N_1 + N_2 + N_7 + N_9}{D}, \quad p_2 = \frac{2N_3 + N_4 + N_7 + N_{11}}{D}, \quad p_3 = \frac{2N_5 + N_6 + N_9 + N_{11}}{D}$$

By construction $p_1 + p_2 + p_3 = 1$.  Each mixed hom-alt individual (e.g. $N_7$, AFR/EUR 1/1) contributes one alt allele to each parent population's numerator and two to the shared denominator $D$.  Ambiguous hets in each mixed pair are then split by the **pairwise ratio** of the two parent populations' proportions.

**Shrinkage toward 1/2.** Each mixed pair's ambiguous heterozygote count ($N_8$ for AFR/EUR, $N_{10}$ for AFR/NAT, $N_{12}$ for EUR/NAT) is shrunk the same way as the 2-ancestry case, toward an even split within that pair. For the AFR/EUR pair:

$$T = D + N_8 + N_{10} + N_{12}, \qquad w_{12} = \frac{N_8}{T}, \qquad \frac{p_1}{p_1+p_2} \leftarrow (1-w_{12})\frac{p_1}{p_1+p_2} + \frac{w_{12}}{2}$$

using $D$ from the raw formula above; $T$ (this variant's total allele-weighted confidence) adds every mixed pair's own ambiguous-het count on top — not just pair $12$'s — matching $D$'s existing 2x weighting of hom-alt carriers. $N_{10}$ (AFR/NAT) and $N_{12}$ (EUR/NAT) shrink identically off the same shared $T$ ($w_{9,10}=N_{10}/T$, $w_{11,12}=N_{12}/T$), each toward $1/2$ within its own pair, so an ambiguous AFR/EUR het never receives NAT dosage. $T$ is shared across all three pairs, so $w_{12}=1$ requires not just $D=0$ but also $N_{10}=N_{12}=0$ — i.e. this variant's *only* alt-carrying evidence, of any pair or population, is AFR/EUR ambiguous hets; that reduces exactly to the 2-ancestry pure-singleton case.

---

#### General K-Ancestry Case

With $K$ ancestries labelled $1, \ldots, K$ there are $K$ pure types and $\binom{K}{2}$ mixed pairs.

**Counts per variant:**

| Symbol | Meaning |
|--------|---------|
| $N_{kk}^{(2)}$ | count of pure pop-$k$ samples with gt = 2 |
| $N_{kk}^{(1)}$ | count of pure pop-$k$ samples with gt = 1 |
| $N_{ij}^{(2)}$ | count of pop-$i$/pop-$j$ mixed samples with gt = 2, $i < j$ |
| $N_{ij}^{(1)}$ | count of pop-$i$/pop-$j$ mixed samples with gt = 1 (ambiguous het), $i < j$ |

**Splitting rules:**

| Ancestry type | gt | $\mathbf{x}_k$ |
|---|---|---|
| pure pop-$k$ | 2 | 2 for $k$, 0 for all others |
| pure pop-$k$ | 1 | 1 for $k$, 0 for all others |
| mixed pop-$i$/pop-$j$ | 2 | 1 for $i$, 1 for $j$, 0 for all others |
| mixed pop-$i$/pop-$j$ | 1 | $\dfrac{p_i}{p_i + p_j}$ for $i$, $\dfrac{p_j}{p_i + p_j}$ for $j$, 0 for others |
| any | 0 | 0 for all |

**Denominator** (total unambiguous alt alleles):

$$D = \sum_{k=1}^{K} \left(2N_{kk}^{(2)} + N_{kk}^{(1)}\right) + 2\sum_{i < j} N_{ij}^{(2)}$$

**Population proportions** ($k = 1, \ldots, K$):

$$p_k = \frac{2N_{kk}^{(2)} + N_{kk}^{(1)} + \displaystyle\sum_{j \neq k} N_{\langle kj \rangle}^{(2)}}{D}$$

where $N_{\langle kj \rangle}^{(2)} = N_{\min(k,j),\max(k,j)}^{(2)}$.

It follows that $\sum_{k=1}^{K} p_k = 1$.

**Interpretation:** $p_k$ estimates the fraction of population-level alt alleles attributable to ancestry $k$, using only observations with unambiguous ancestry.  The mixed-pair hom-alts ($N_{ij}^{(2)}$) contribute one allele to each parent population.  Ambiguous hets in pair $(i,j)$ are then split by the conditional pairwise ratio $p_i / (p_i + p_j)$, i.e. the AFR fraction among AFR and EUR only for an AFR/EUR het, regardless of how many other populations exist.

**Singleton special case:** if the only alt carriers in a mixed pair $(i,j)$ are heterozygotes and all pure-ancestry and hom-alt counts are zero, $D=0$ for that pair and the raw ratio is undefined.

**Shrinkage toward 1/2 (general $K$).** The formula above generalises **per pair**, not per population: each mixed pair $(i,j)$ gets its own ambiguous-fraction weight, and is shrunk toward an even split within the pair.

$$T = D + \sum_{i' \lt j'} N_{i'j'}^{(1)}, \qquad w_{ij} = \frac{N_{ij}^{(1)}}{T}, \qquad \frac{p_i}{p_i+p_j} \leftarrow (1-w_{ij})\frac{p_i}{p_i+p_j} + \frac{w_{ij}}{2}$$

using $D$ from the denominator formula above and the same $p_i$ from the population-proportion formula (the full-variant estimate, not a pair-restricted recount) — so the raw ratio is defined even when pair $(i,j)$ itself has no unambiguous hom-alt carriers, as long as $i$ and $j$ each have unambiguous evidence elsewhere in the variant; if $p_i+p_j=0$ too, the pair's split is $1/2$. $T$ (this variant's total allele-weighted confidence) is shared across *every* pair, not recomputed per pair: it sums $D$ (which already counts hom-alt carriers, pure or mixed, 2x) plus *all* pairs' ambiguous-het counts, not just $(i,j)$'s. $w_{ij}=0$ recovers the raw pairwise ratio; $w_{ij}=1$ (this variant's only alt-carrying evidence, of any pair or population, is pair $(i,j)$'s own ambiguous hets — the singleton case above) gives exactly $1/2 : 1/2$.

**Number of ancestry types by K:**

| K | Pure types | Mixed pairs | Total types |
|---|-----------|-------------|-------------|
| 2 | 2 | 1 | 3 |
| 3 | 3 | 3 | 6 |
| 4 | 4 | 6 | 10 |
| 5 | 5 | 10 | 15 |

---

### Association Detection

* Per-ancestry association: $\mathbf{x}_k \sim \text{Response} \to p_k^{\text{assoc}}$ for each $k = 1,\ldots,K$
* Observed association: $\mathbf{x}_{obs} \sim \text{Response} \to p_{obs}$
* Cauchy combination of per-ancestry p-values $\to p_c$

For the 2-ancestry case: pure African $\to p_a$, pure European $\to p_e$, Cauchy combination $\to p_c$.

#### Per-gene ancestry weights for the Cauchy combination

Not every population contributes equally reliable evidence for every gene:
a gene sitting in a genomic region where only a handful of cohort samples
have pure ancestry $k$ gives $p_k^{\text{assoc}}$ little to work with, and
shouldn't count as much as a population with abundant pure-ancestry
representation there. LANTERN weights each population's p-value by how much
pure-ancestry evidence backs it *at that specific gene*, rather than
combining all $p_k$ equally:

$$w_k(\text{gene}) = \text{median}_{v \in \text{gene}}\bigl(\text{count of cohort samples with pure ancestry } k \text{ at variant } v\bigr)$$

taking the median across the gene's variants (from `ancestry_split()`'s
`ancestry_counts` output) so a handful of unusually ancestry-rich or
ancestry-poor variants at the gene's edges don't dominate the estimate.
These weights feed `cauchy_combine()`:

$$p_c = \mathrm{CCT}\bigl(p_1^{\text{assoc}}, \ldots, p_K^{\text{assoc}} ; w_1, \ldots, w_K\bigr)$$

so a population with little pure-ancestry representation at this gene is
discounted, not ignored outright — a gene well-represented in one ancestry
but barely in another naturally leans toward the better-supported
population's signal. If none of the populations' names match
`ancestry_counts`'s columns, `ancestry_smmat()` falls back to equal
weighting across all $K$ populations instead.

### 4. R Package Installation

The `lantern` R package is on GitHub. It has two categories of dependencies:

- **Bioconductor R package** — [SeqArray](https://bioconductor.org/packages/SeqArray/) is not on CRAN; `devtools::install_github()` will not find it automatically.
- **Optional system binary** — when [`bcftools`](https://samtools.github.io/bcftools/) (≥ 1.10) is on `PATH`, VCF reading uses it (faster, and required for BCF input). Otherwise lantern reads VCFs through SeqArray. Install it via your system package manager (`brew install bcftools` on macOS, `apt install bcftools` on Debian/Ubuntu).
- **`GMMAT`** (CRAN) — only needed for Step 2/3 (`ancestry_smmat()`), not for Step 1 splitting alone: `install.packages("GMMAT")`.

**Recommended installation** (handles the Bioconductor dependency):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("yuw444/LANTERN/lantern")
```

**Alternative**: pre-install the Bioconductor dependency first, then use devtools:

```r
BiocManager::install("SeqArray")
devtools::install_github("yuw444/LANTERN", subdir = "lantern")
```

### 5. Pipeline (legacy `./src`)

Three SLURM-friendly CLI scripts (`step0_observed_association.R`,
`step1_vcf_split_by_ancestry.R`, `step2_association_detection.R`). Step 1
and Step 2 wrap the `lantern` R package — `optparse` argument parsing and
file I/O around `lantern::ancestry_split()` / `write_ancestry_gds()` /
`ancestry_smmat()`. Step 0 calls `GMMAT::glmmkin()`/`GMMAT::SMMAT()`
directly (no `lantern` calls) to produce $p_{obs}$, the plain,
non-ancestry-split association p-value from the "Association Detection"
section above — it is entirely independent of Steps 1/2 and of each
other's outputs; run any subset, in any order. Use these if you want each
step as its own SLURM job; if you're scripting in R directly, just call
the package functions yourself (see `vignette("split-intuition")`,
`vignette("toy-vcf-msp-example")`, `vignette("real-data-chr19")`). There
is no separate weight-finding script — per-gene ancestry weights for the
Cauchy combination are computed automatically inside Step 2.

#### Prerequisites

* The `lantern` R package installed (Section 4 above)
* `GMMAT`: `install.packages("GMMAT")` (CRAN). If you're running from a clone
  of this repo via `pixi`, its conda-forge R environment doesn't have a CRAN
  mirror configured by default — use `pixi run install-cran-gmmat` instead.
* Optional: `bcftools` ≥ 1.10 on `PATH` for faster VCF reading (required for BCF input)
* An RFMix2-produced `.msp.tsv` file. There's no vignette covering RFMix2
  itself yet — see [RFMix2's own docs](https://github.com/slowkoni/rfmix)
  for producing one from a phased VCF + reference panel.

No installation of `./src` itself is needed — just run the scripts with `Rscript`.

#### Step 0 — plain association test on observed (unsplit) genotypes

```bash
Rscript src/step0_observed_association.R \
  --vcf_path        cohort.phased.vcf.gz \
  --chr_id          22 \
  --out_path        out/ \
  --data_file       pheno.tsv \
  --gene_group_file genes.tsv \
  --kinship_rds     kinship.rds \
  --response_type   continuous \
  --out_file        out/observed_results.tsv \
  --ncores          8
```

Converts `--vcf_path` (filtered to `--chr_id`) directly to
`out/OBSERVED.gds` via `SeqArray::seqVCF2GDS()` — true hard-call
genotypes read straight from the input file, no ancestry-split dosage
reconstruction — then fits its own null model (`GMMAT::glmmkin()`) and
runs a plain `GMMAT::SMMAT(..., is.dosage = FALSE)` on it. This is
$p_{obs}$ from the "Association Detection" section above: the ordinary,
non-ancestry-aware rare-variant test, for comparison against the
ancestry-stratified Step 1/Step 2 results.

| Flag | Required | Meaning |
|------|----------|---------|
| `--vcf_path` / `-i` | yes | Same phased VCF/BCF you'd pass to Step 1's `--vcf_path`. |
| `--chr_id` / `-c` | yes | Chromosome to process, must match `--vcf_path`'s `CHROM` column (e.g. `22` or `chr22`, tried both ways). |
| `--out_path` / `-o` | yes | Output directory (created if missing); `OBSERVED.gds` is written here. |
| `--data_file` | yes | Phenotype file — same format as Step 2's. |
| `--gene_group_file` | yes | Gene-group file passed straight to `GMMAT::SMMAT()` — same format as Step 2's. |
| `--kinship_rds` | yes | Kinship matrix RDS — same format as Step 2's. |
| `--response_type` | no (default `continuous`) | `continuous` (`gaussian`), `binary` (`binomial`), or `count` (`poisson`). |
| `--out_file` | yes | Output TSV path — see Output below. |
| `--ncores` | no | Cores for `GMMAT::SMMAT()`. Defaults to `$SLURM_CPUS_PER_TASK` when set, else `1`. |

**Output**: `out/OBSERVED.gds` plus two files derived from `--out_file`:
* `--out_file` itself — a two-column tab-separated table, `gene` and
  `p_OBSERVED`.
* `<out_file, .tsv swapped for .rds>` — an RDS of
  `list(results, smmat_results)`, where `results` is the same table as the
  TSV above and `smmat_results` is the raw `GMMAT::SMMAT()` output
  data.frame.

#### Step 1 — split a VCF by local ancestry

```bash
Rscript src/step1_vcf_split_by_ancestry.R \
  --vcf_path  cohort.phased.vcf.gz \
  --msp_path  cohort.msp.tsv \
  --out_path  out/ \
  --chr_id    22 \
  --mode      dosage
```

| Flag | Required | Meaning |
|------|----------|---------|
| `--vcf_path` / `-i` | yes | Phased VCF/BCF (bgzipped). May contain other chromosomes besides `--chr_id` — see "Multi-chromosome input" below. |
| `--msp_path` | yes | RFMix `.msp.tsv` (plain text or gzipped). May also contain other chromosomes' tracts; only `--chr_id`'s are used. |
| `--out_path` | yes | Output directory, created if missing. |
| `--chr_id` | yes | Chromosome to process, e.g. `22` or `chr22` — must identify the same chromosome in both `--vcf_path` and `--msp_path` (a `chr` prefix mismatch between the two is handled automatically). |
| `--mode` | no (default `dosage`) | `dosage` = proportional p1/p2 split (unphased-friendly) or `haplotype` = deterministic per-haplotype split (needs a truly phased VCF). See `vignette("split-intuition")` for the difference. |

**Multi-chromosome input**: `--vcf_path`/`--msp_path` don't need to be
pre-split per chromosome — `--chr_id` selects one chromosome out of a
genome-wide file. When possible (the VCF's contigs can be resolved),
`bcftools query` itself is restricted to `--chr_id`, so the other
chromosomes' genotypes are never read into R — this keeps memory bounded
to one chromosome's variant count even if you hand it a whole-genome VCF.
Run one Step 1 job per chromosome (e.g. a SLURM array over `--chr_id`)
rather than combining chromosomes into a single call.

**Output** (in `--out_path`):

| File | Contents |
|------|----------|
| `<POP>.gds` (one per population named in the MSP header, e.g. `AFR.gds`, `EUR.gds`) | Ancestry-specific dosage GDS, ready for `GMMAT::SMMAT()` (`is.dosage = TRUE`) |
| `split_meta_chr<chr_id>.rds` | `list(gds_paths, variant_info, ancestry_counts, sample_ids, mode, shrink_target, lantern_version, chr_id)` — bundles the GDS paths above plus per-variant ancestry counts. This whole file is Step 2's `--split_meta` input. |

#### Step 2 — ancestry-stratified association testing

```bash
Rscript src/step2_association_detection.R \
  --split_meta      out/split_meta_chr22.rds \
  --data_file       pheno.tsv \
  --gene_group_file genes.tsv \
  --kinship_rds     kinship.rds \
  --response_type   continuous \
  --out_file        results.tsv \
  --ncores          8
```

| Flag | Required | Meaning |
|------|----------|---------|
| `--split_meta` | yes | `split_meta_chr<chr_id>.rds` from Step 1 — supplies the GDS paths for every population plus per-variant ancestry counts (used to weight the Cauchy combination; equal weights if missing). |
| `--data_file` | yes | Phenotype file, csv/tsv/rds — see format below. |
| `--gene_group_file` | yes | Gene-group file passed straight to `GMMAT::SMMAT()` — see format below. |
| `--kinship_rds` | yes | Kinship matrix RDS — see format below. |
| `--response_type` | no (default `continuous`) | `continuous` (`gaussian`), `binary` (`binomial`), or `count` (`poisson`). |
| `--out_file` | yes | Output TSV path — see Output below. |
| `--ncores` | no | Cores for `GMMAT::SMMAT()`. Defaults to `$SLURM_CPUS_PER_TASK` when set (i.e. automatically picks up `--cpus-per-task` in a SLURM job), else `1`. Pass explicitly to override. |

**Input file formats**:

* **Phenotype** (`--data_file`, header row required): column 1 must be
  `id`, column 2 is the response, remaining columns are covariates (all
  used — the formula is built as `<col2> ~ <col3> + <col4> + ...`).
  ```
  id	response	age	sex	PC1	PC2
  sample_001	1	63	M	-0.012	0.034
  sample_002	0	57	F	0.104	-0.021
  ```
* **Gene group** (`--gene_group_file`, no header): `gene, chr, pos, ref,
  alt, weight`, one row per variant per gene.
  ```
  GOLGA6L22	15	22460882	G	T	1
  GOLGA6L22	15	22462401	G	C	1
  HERC2P2	15	22554572	G	A	1
  ```
  **The `chr` column must be a bare chromosome number/name with no `chr`
  prefix** (`15`, not `chr15`) — regardless of whether your original VCF
  uses `chr15` or `15` internally. Step 1's GDS files are written via
  `SeqArray::seqVCF2GDS()`, which always strips any `chr` prefix when it
  stores the chromosome field; `GMMAT::SMMAT()` then matches variants to
  genes by plain string equality (`chr:pos:ref:alt`) against exactly that
  stripped value. A mismatch doesn't error: every gene whose variants fail
  to match silently drops out of the result entirely.
* **Kinship** (`--kinship_rds`): an RDS of a square numeric matrix with
  row and column names equal to sample IDs.

**Output**: two files, both derived from `--out_file`:
* `--out_file` itself — the per-gene results table as tab-separated text
  (via `data.table::fwrite()`): `gene`, one `p_<POP>` and one `w_<POP>`
  column per population (e.g. `p_AFR`, `w_AFR`, `p_EUR`, `w_EUR`), and
  `p_cauchy` (the combined p-value). This is `ancestry_smmat()`'s
  `results` data.frame, and the one most people want.
* `<out_file, .tsv swapped for .rds>` — an RDS of `ancestry_smmat()`'s
  full return value, `list(results, smmat_results)`, where `results` is
  the same table as the TSV above and `smmat_results` is a named list (one
  entry per population) of the raw `GMMAT::SMMAT()` output data.frames —
  for anyone who needs the full per-gene detail (variant counts,
  missingness, MAF range, etc.) beyond just the p-values.
