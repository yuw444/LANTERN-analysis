#ifndef ANCESTRY_H
#define ANCESTRY_H

#include <R.h>
#include <Rinternals.h>

/* Count ancestry codes in PT matrix */
SEXP count_ancestry_codes(SEXP mat, SEXP code);

/* Split genotype matrix by ancestry */
SEXP split_by_ancestry(SEXP gt_genotype, SEXP ancestry);

/* Read PLINK .bed/.bim/.fam ancestry file */
SEXP read_bed_file(SEXP bed_path, SEXP bim_path, SEXP fam_path, SEXP sample_indices);

/* Write ancestry-specific VCFs */
SEXP write_vcf_with_ancestry(SEXP vcf_path, SEXP gt_matrix, SEXP ancestry_matrix,
                              SEXP output_african, SEXP output_european);

/* Subset VCF by genomic range */
SEXP subset_vcf_by_range(SEXP vcf_path, SEXP chrom, SEXP start, SEXP end, SEXP output_path);

/* Split phased haplotypes by per-haplotype local ancestry */
SEXP split_phased_by_ancestry(SEXP gt_hap0, SEXP gt_hap1,
                               SEXP anc_hap0, SEXP anc_hap1,
                               SEXP pop_codes);

#endif
