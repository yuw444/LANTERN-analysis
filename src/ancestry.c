#include "ancestry.h"

// ============================================================================
// Helper: Count ancestry-genotype combinations for a variant
// ============================================================================
// Returns counts for a single variant (row):
// N1 = pt==3 & gt==2, N2 = pt==3 & gt==1
// N4 = pt==2 & gt==2, N5 = pt==2 & gt==1 (singleton case)
// N7 = pt==1 & gt==2, N8 = pt==1 & gt==1
typedef struct {
    int N1, N2, N4, N5, N7, N8;
} VariantCounts;

static VariantCounts count_variant_combinations(int *gt_ptr, int *an_ptr, int nrow, int ncol, int row) {
    VariantCounts cnt = {0, 0, 0, 0, 0, 0};
    
    for (int j = 0; j < ncol; j++) {
        int gt = gt_ptr[row + nrow * j];
        int an = an_ptr[row + nrow * j];
        
        if (an == 3 && gt == 2) cnt.N1++;
        else if (an == 3 && gt == 1) cnt.N2++;
        else if (an == 2 && gt == 2) cnt.N4++;
        else if (an == 2 && gt == 1) cnt.N5++;
        else if (an == 1 && gt == 2) cnt.N7++;
        else if (an == 1 && gt == 1) cnt.N8++;
    }
    
    return cnt;
}

// ============================================================================
// count_ancestry_codes: Count occurrences of a code per row
// ============================================================================
static SEXP count_ancestry_codes_c(SEXP mat, SEXP code) {
    SEXP mat_int = PROTECT(coerceVector(mat, INTSXP));
    SEXP dim = PROTECT(getAttrib(mat_int, R_DimSymbol));
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    int target_code = INTEGER(code)[0];
    
    SEXP result = PROTECT(allocVector(INTSXP, nrow));
    int *res_ptr = INTEGER(result);
    int *mat_ptr = INTEGER(mat_int);
    
    for (int i = 0; i < nrow; i++) {
        res_ptr[i] = 0;
        for (int j = 0; j < ncol; j++) {
            if (mat_ptr[i + nrow * j] == target_code) {
                res_ptr[i]++;
            }
        }
    }
    
    UNPROTECT(3);
    return result;
}

// ============================================================================
// split_by_ancestry: Split genotype matrix by ancestry with proper p1/p2
// ============================================================================
static SEXP split_by_ancestry_c(SEXP gt_genotype, SEXP ancestry) {
    SEXP gt_int = PROTECT(coerceVector(gt_genotype, INTSXP));
    SEXP an_int = PROTECT(coerceVector(ancestry, INTSXP));
    SEXP dim = PROTECT(getAttrib(gt_int, R_DimSymbol));
    
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    
    // Allocate output matrices
    SEXP african = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    SEXP european = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    
    int *gt_ptr = INTEGER(gt_int);
    int *an_ptr = INTEGER(an_int);
    double *afr_ptr = REAL(african);
    double *eur_ptr = REAL(european);
    
    // Process each variant (row)
    for (int i = 0; i < nrow; i++) {
        // Count ancestry-genotype combinations for this variant
        VariantCounts cnt = count_variant_combinations(gt_ptr, an_ptr, nrow, ncol, i);
        
        // Calculate p1 and p2
        // p1 = (2*N1 + N2 + N4) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
        // p2 = (N4 + 2*N7 + N8) / (2*N1 + N2 + 2*N4 + 2*N7 + N8)
        // 
        // Singleton case: if N5 == sum(gt), p1 = p2 = 0.5
        
        double total_alt = 2 * cnt.N1 + cnt.N2 + cnt.N4 + cnt.N5 + 2 * cnt.N7 + cnt.N8;
        double p1, p2;
        
        if (cnt.N5 > 0 && total_alt == (double)cnt.N5) {
            // Singleton case: all alt alleles are from mixed het individuals
            p1 = 0.5;
            p2 = 0.5;
        } else {
            // Denominator excludes heterozygous mixed (N5) - those are split
            double denominator = total_alt - cnt.N5;
            
            if (denominator > 0) {
                double numerator_p1 = 2.0 * cnt.N1 + cnt.N2 + cnt.N4;
                double numerator_p2 = cnt.N4 + 2.0 * cnt.N7 + cnt.N8;
                p1 = numerator_p1 / denominator;
                p2 = numerator_p2 / denominator;
            } else {
                // Edge case: no non-singleton data
                p1 = 0.5;
                p2 = 0.5;
            }
        }
        
        // Apply splitting to each sample
        for (int j = 0; j < ncol; j++) {
            int gt = gt_ptr[i + nrow * j];
            int an = an_ptr[i + nrow * j];
            int idx = i + nrow * j;
            
            if (an == 3) {
                // Pure African (EUR/EUR = 0)
                afr_ptr[idx] = gt;  // All alt alleles from African
                eur_ptr[idx] = 0.0;
            } else if (an == 1) {
                // Pure European
                afr_ptr[idx] = 0.0;
                eur_ptr[idx] = gt;  // All alt alleles from European
            } else if (an == 2) {
                // Mixed ancestry: split based on p1/p2
                if (gt == 2) {
                    // Homozygous alt: 1 allele to each ancestry
                    afr_ptr[idx] = 1.0;
                    eur_ptr[idx] = 1.0;
                } else if (gt == 1) {
                    // Heterozygous: split by p1/p2
                    afr_ptr[idx] = p1;
                    eur_ptr[idx] = p2;
                } else {
                    // Homozygous ref
                    afr_ptr[idx] = 0.0;
                    eur_ptr[idx] = 0.0;
                }
            } else {
                // Invalid/missing ancestry code
                afr_ptr[idx] = 0.0;
                eur_ptr[idx] = 0.0;
            }
        }
    }
    
    // Create named list result
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, african);
    SET_VECTOR_ELT(result, 1, european);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("african"));
    SET_STRING_ELT(names, 1, mkChar("european"));
    setAttrib(result, R_NamesSymbol, names);
    
    UNPROTECT(7);
    return result;
}

// ============================================================================
// read_bed_file: Read PLINK binary files
// ============================================================================
static SEXP read_bed_file_c(SEXP bed_path, SEXP bim_path, SEXP fam_path, SEXP sample_indices) {
    const char *bed = CHAR(STRING_ELT(bed_path, 0));
    const char *bim = CHAR(STRING_ELT(bim_path, 0));
    const char *fam = CHAR(STRING_ELT(fam_path, 0));
    
    FILE *fp = fopen(bed, "rb");
    if (!fp) error("Cannot open bed file: %s", bed);
    
    unsigned char magic[3];
    if (fread(magic, 1, 3, fp) != 3) {
        fclose(fp);
        error("Cannot read bed file header");
    }
    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
        fclose(fp);
        error("Invalid bed file format");
    }
    
    int n_samples = 0;
    fp = fopen(fam, "r");
    if (!fp) error("Cannot open fam file: %s", fam);
    char line[1000];
    while (fgets(line, sizeof(line), fp)) n_samples++;
    fclose(fp);
    
    int n_variants = 0;
    fp = fopen(bim, "r");
    if (!fp) error("Cannot open bim file: %s", bim);
    while (fgets(line, sizeof(line), fp)) n_variants++;
    fclose(fp);
    
    int store_size = (n_samples + 3) / 4;
    
    SEXP result = PROTECT(allocVector(INTSXP, n_samples * n_variants));
    int *gt = INTEGER(result);
    memset(gt, 0, n_samples * n_variants * sizeof(int));
    
    fp = fopen(bed, "rb");
    fseek(fp, 3, SEEK_SET);
    
    for (int v = 0; v < n_variants; v++) {
        unsigned char buffer[store_size];
        if (fread(buffer, 1, store_size, fp) != (size_t)store_size) {
            fclose(fp);
            error("Truncated bed file at variant %d", v);
        }
        
        for (int s = 0; s < n_samples; s++) {
            int byte_idx = s / 4;
            int bit_idx = (s % 4) * 2;
            int genotype = (buffer[byte_idx] >> bit_idx) & 0x03;
            
            if (genotype == 3) gt[s + n_samples * v] = 3;  // Missing
            else if (genotype == 2) gt[s + n_samples * v] = 2;  // Homozygous alt
            else if (genotype == 1) gt[s + n_samples * v] = 1;  // Heterozygous
            else gt[s + n_samples * v] = 0;  // Homozygous ref
        }
    }
    fclose(fp);
    
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = n_samples;
    INTEGER(dim)[1] = n_variants;
    setAttrib(result, R_DimSymbol, dim);
    
    UNPROTECT(2);
    return result;
}

// ============================================================================
// write_vcf_with_ancestry: Write VCF with ancestry-specific dosages
// ============================================================================
static SEXP write_vcf_with_ancestry_c(SEXP vcf_path, SEXP gt_matrix, SEXP ancestry_matrix, 
                                        SEXP output_african, SEXP output_european) {
    const char *vcf_in = CHAR(STRING_ELT(vcf_path, 0));
    const char *vcf_afr = CHAR(STRING_ELT(output_african, 0));
    const char *vcf_eur = CHAR(STRING_ELT(output_european, 0));
    
    int nrow = INTEGER(gt_matrix)[0];
    
    FILE *in = fopen(vcf_in, "r");
    if (!in) error("Cannot open input VCF: %s", vcf_in);
    
    FILE *out_afr = fopen(vcf_afr, "w");
    if (!out_afr) { fclose(in); error("Cannot open output VCF: %s", vcf_afr); }
    
    FILE *out_eur = fopen(vcf_eur, "w");
    if (!out_eur) { fclose(in); fclose(out_afr); error("Cannot open output VCF: %s", vcf_eur); }
    
    char line[10000];
    while (fgets(line, sizeof(line), in)) {
        if (line[0] == '#') {
            fprintf(out_afr, "%s", line);
            fprintf(out_eur, "%s", line);
            continue;
        }
        
        char *fields[10];
        int n = 0;
        char *ptr = line;
        while ((fields[n] = strtok(ptr, "\t\n")) != NULL) {
            n++;
            ptr = NULL;
        }
        
        fprintf(out_afr, "%s\t%s\t%s\t%s\t%s\t", fields[0], fields[1], fields[2], fields[3], fields[4]);
        fprintf(out_eur, "%s\t%s\t%s\t%s\t%s\t", fields[0], fields[1], fields[2], fields[3], fields[4]);
        
        for (int i = 9; i < n; i++) {
            if (i > 9) {
                fprintf(out_afr, "\t");
                fprintf(out_eur, "\t");
            }
            
            int sample_idx = i - 9;
            int gt = INTEGER(gt_matrix)[sample_idx + nrow * (i - 9)];
            int an = INTEGER(ancestry_matrix)[sample_idx + nrow * (i - 9)];
            
            char gt_afr[20], gt_eur[20];
            double ds_afr, ds_eur;
            
            if (an == 3) {
                snprintf(gt_afr, 20, "%d/%d", gt > 0 ? 1 : 0, gt > 1 ? 1 : 0);
                snprintf(gt_eur, 20, "0/0");
                ds_afr = gt > 0 ? (gt > 1 ? 2.0 : 1.0) : 0.0;
                ds_eur = 0.0;
            } else if (an == 1) {
                snprintf(gt_afr, 20, "0/0");
                snprintf(gt_eur, 20, "%d/%d", gt > 0 ? 1 : 0, gt > 1 ? 1 : 0);
                ds_afr = 0.0;
                ds_eur = gt > 0 ? (gt > 1 ? 2.0 : 1.0) : 0.0;
            } else if (an == 2) {
                if (gt == 2) {
                    snprintf(gt_afr, 20, "1/0");
                    snprintf(gt_eur, 20, "0/1");
                    ds_afr = 1.0;
                    ds_eur = 1.0;
                } else if (gt == 1) {
                    snprintf(gt_afr, 20, "0/1");
                    snprintf(gt_eur, 20, "0/1");
                    ds_afr = 0.5;
                    ds_eur = 0.5;
                } else {
                    snprintf(gt_afr, 20, "0/0");
                    snprintf(gt_eur, 20, "0/0");
                    ds_afr = 0.0;
                    ds_eur = 0.0;
                }
            } else {
                snprintf(gt_afr, 20, "0/0");
                snprintf(gt_eur, 20, "0/0");
                ds_afr = 0.0;
                ds_eur = 0.0;
            }
            
            if (i == 9) {
                fprintf(out_afr, "GT:DS");
                fprintf(out_eur, "GT:DS");
            }
            
            fprintf(out_afr, "\t%s:%.2f", gt_afr, ds_afr);
            fprintf(out_eur, "\t%s:%.2f", gt_eur, ds_eur);
        }
        fprintf(out_afr, "\n");
        fprintf(out_eur, "\n");
    }
    
    fclose(in);
    fclose(out_afr);
    fclose(out_eur);
    
    return R_NilValue;
}

// ============================================================================
// subset_vcf_by_range: Extract VCF region
// ============================================================================
static SEXP subset_vcf_by_range_c(SEXP vcf_path, SEXP chrom, SEXP start, SEXP end, SEXP output_path) {
    const char *in_vcf = CHAR(STRING_ELT(vcf_path, 0));
    const char *out = CHAR(STRING_ELT(output_path, 0));
    const char *chr = CHAR(STRING_ELT(chrom, 0));
    int pos_start = INTEGER(start)[0];
    int pos_end = INTEGER(end)[0];
    
    FILE *in = fopen(in_vcf, "r");
    if (!in) error("Cannot open VCF: %s", in_vcf);
    
    FILE *outf = fopen(out, "w");
    if (!outf) { fclose(in); error("Cannot open output: %s", out); }
    
    char line[10000];
    while (fgets(line, sizeof(line), in)) {
        if (line[0] == '#') {
            fprintf(outf, "%s", line);
            continue;
        }
        
        char *fields[10];
        int n = 0;
        char *ptr = line;
        while ((fields[n] = strtok(ptr, "\t\n")) != NULL) {
            n++;
            ptr = NULL;
        }
        
        if (n < 2) continue;
        
        if (strcmp(fields[0], chr) == 0) {
            int pos = atoi(fields[1]);
            if (pos >= pos_start && pos <= pos_end) {
                fprintf(outf, "%s", line);
            }
        }
    }
    
    fclose(in);
    fclose(outf);
    
    return R_NilValue;
}

// ============================================================================
// split_phased_by_ancestry: Deterministic split using per-haplotype ancestry
// ============================================================================
// Each haplotype contributes its allele to whichever ancestry population it
// was called as.  NA genotypes are treated as 0 (reference).  Any ancestry
// code that is neither AFR nor EUR contributes nothing to either output.
//
// pop_codes: integer vector of length 2 — pop_codes[0] = AFR code,
//            pop_codes[1] = EUR code (RFMix default: 0 and 1).
static SEXP split_phased_by_ancestry_c(SEXP gt_hap0, SEXP gt_hap1,
                                        SEXP anc_hap0, SEXP anc_hap1,
                                        SEXP pop_codes) {
    SEXP dim = PROTECT(getAttrib(gt_hap0, R_DimSymbol));
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    int afr_code = INTEGER(pop_codes)[0];
    int eur_code = INTEGER(pop_codes)[1];

    SEXP gh0 = PROTECT(coerceVector(gt_hap0, INTSXP));
    SEXP gh1 = PROTECT(coerceVector(gt_hap1, INTSXP));
    SEXP ah0 = PROTECT(coerceVector(anc_hap0, INTSXP));
    SEXP ah1 = PROTECT(coerceVector(anc_hap1, INTSXP));

    SEXP african  = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    SEXP european = PROTECT(allocMatrix(REALSXP, nrow, ncol));

    int *gh0_ptr = INTEGER(gh0);
    int *gh1_ptr = INTEGER(gh1);
    int *ah0_ptr = INTEGER(ah0);
    int *ah1_ptr = INTEGER(ah1);
    double *afr_ptr = REAL(african);
    double *eur_ptr = REAL(european);

    int n = nrow * ncol;
    for (int k = 0; k < n; k++) {
        int g0 = (gh0_ptr[k] == NA_INTEGER) ? 0 : gh0_ptr[k];
        int g1 = (gh1_ptr[k] == NA_INTEGER) ? 0 : gh1_ptr[k];
        int a0 = ah0_ptr[k];
        int a1 = ah1_ptr[k];

        double afr = 0.0, eur = 0.0;
        if (a0 == afr_code)      afr += g0;
        else if (a0 == eur_code) eur += g0;

        if (a1 == afr_code)      afr += g1;
        else if (a1 == eur_code) eur += g1;

        afr_ptr[k] = afr;
        eur_ptr[k] = eur;
    }

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, african);
    SET_VECTOR_ELT(result, 1, european);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("african"));
    SET_STRING_ELT(names, 1, mkChar("european"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(9);
    return result;
}

// ============================================================================
// Exported functions (registered in init.c)
// ============================================================================
SEXP count_ancestry_codes(SEXP mat, SEXP code) {
    return count_ancestry_codes_c(mat, code);
}

SEXP split_by_ancestry(SEXP gt_genotype, SEXP ancestry) {
    return split_by_ancestry_c(gt_genotype, ancestry);
}

SEXP read_bed_file(SEXP bed_path, SEXP bim_path, SEXP fam_path, SEXP sample_indices) {
    return read_bed_file_c(bed_path, bim_path, fam_path, sample_indices);
}

SEXP write_vcf_with_ancestry(SEXP vcf_path, SEXP gt_matrix, SEXP ancestry_matrix, 
                              SEXP output_african, SEXP output_european) {
    return write_vcf_with_ancestry_c(vcf_path, gt_matrix, ancestry_matrix, output_african, output_european);
}

SEXP subset_vcf_by_range(SEXP vcf_path, SEXP chrom, SEXP start, SEXP end, SEXP output_path) {
    return subset_vcf_by_range_c(vcf_path, chrom, start, end, output_path);
}

SEXP split_phased_by_ancestry(SEXP gt_hap0, SEXP gt_hap1,
                               SEXP anc_hap0, SEXP anc_hap1,
                               SEXP pop_codes) {
    return split_phased_by_ancestry_c(gt_hap0, gt_hap1, anc_hap0, anc_hap1, pop_codes);
}
