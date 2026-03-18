#ifndef ANCESTRY_H
#define ANCESTRY_H

#include <R.h>
#include <Rinternals.h>

/* Count ancestry codes in PT matrix */
SEXP count_ancestry_codes(SEXP mat, SEXP code);

/* Split genotype matrix by ancestry */
SEXP split_by_ancestry(SEXP gt_genotype, SEXP ancestry);

/* Read RFMix/PLINK ancestry file (.bed/.bim/.fam) */
SEXP read_ancestry_plink(SEXP bed_path, SEXP bim_path, SEXP fam_path);

/* Read PT matrix from TSV file */
SEXP read_pt_matrix(SEXP path);

/* Write ancestry-specific VCFs */
SEXP write_ancestry_vcf(SEXP vcf_path, SEXP gt_matrix, SEXP ancestry_matrix,
                         SEXP out_african, SEXP out_european);

/* Subset VCF by genomic range */
SEXP subset_vcf_by_range(SEXP vcf_path, SEXP chrom, SEXP start, SEXP end, SEXP output_path);

#endif
