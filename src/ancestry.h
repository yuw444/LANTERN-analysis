#ifndef ANCESTRY_H
#define ANCESTRY_H

#include <R.h>
#include <Rinternals.h>

/* Count ancestry codes in PT matrix */
SEXP count_ancestry_codes(SEXP mat, SEXP code);

/* Split genotype matrix by ancestry (gla/arm_id: GLA-shrinkage; pass a
   zero-length gla to disable and reproduce pre-shrinkage behavior) */
SEXP split_by_ancestry(SEXP gt_genotype, SEXP ancestry, SEXP gla, SEXP arm_id);

/* Read PLINK .bed/.bim/.fam ancestry file */
SEXP read_bed_file(SEXP bed_path, SEXP bim_path, SEXP fam_path);

/* Split phased haplotypes by per-haplotype local ancestry (2-population) */
SEXP split_phased_by_ancestry(SEXP gt_hap0, SEXP gt_hap1,
                               SEXP anc_hap0, SEXP anc_hap1,
                               SEXP pop_codes);

/* Split phased haplotypes into K population-specific dosage matrices */
SEXP split_phased_multi(SEXP gt_hap0, SEXP gt_hap1,
                         SEXP anc_hap0, SEXP anc_hap1,
                         SEXP pop_codes);

/* Split unphased genotype matrix into K population-specific dosage matrices
   (gla/arm_id: GLA-shrinkage; pass a zero-length gla to disable) */
SEXP split_by_ancestry_multi(SEXP gt, SEXP ancestry, SEXP pure_codes,
                               SEXP m_code, SEXP m_pop1, SEXP m_pop2,
                               SEXP gla, SEXP arm_id);

#endif
