#include "ancestry_code.h"

static SEXP count_ancestry_codes_c(SEXP mat, SEXP code) {
    SEXP mat_int = PROTECT(coerceVector(mat, INTSXP));
    SEXP dim = PROTECT(getAttrib(mat_int, R_DimSymbol));
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    int target_code = INTEGER(code)[0];
    
    SEXP result = PROTECT(allocVector(INTSXP, nrow));
    int *res_ptr = INTEGER(result);
    
    for (int i = 0; i < nrow; i++) {
        res_ptr[i] = 0;
        for (int j = 0; j < ncol; j++) {
            if (INTEGER(mat_int)[i + nrow * j] == target_code) {
                res_ptr[i]++;
            }
        }
    }
    
    UNPROTECT(3);
    return result;
}

static SEXP split_by_ancestry_c(SEXP gt_genotype, SEXP ancestry) {
    SEXP gt_int = PROTECT(coerceVector(gt_genotype, INTSXP));
    SEXP an_int = PROTECT(coerceVector(ancestry, INTSXP));
    SEXP dim = PROTECT(getAttrib(gt_int, R_DimSymbol));
    
    int nrow = INTEGER(dim)[0];
    int ncol = INTEGER(dim)[1];
    
    SEXP african = PROTECT(allocMatrix(INTSXP, nrow, ncol));
    SEXP european = PROTECT(allocMatrix(INTSXP, nrow, ncol));
    
    int *gt_ptr = INTEGER(gt_int);
    int *an_ptr = INTEGER(an_int);
    int *afr_ptr = INTEGER(african);
    int *eur_ptr = INTEGER(european);
    
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            int gt = gt_ptr[i + nrow * j];
            int an = an_ptr[i + nrow * j];
            
            if (an == 3) {
                afr_ptr[i + nrow * j] = gt * 2;
                eur_ptr[i + nrow * j] = 0;
            } else if (an == 1) {
                afr_ptr[i + nrow * j] = 0;
                eur_ptr[i + nrow * j] = gt * 2;
            } else if (an == 2) {
                afr_ptr[i + nrow * j] = gt;
                eur_ptr[i + nrow * j] = gt;
            } else {
                afr_ptr[i + nrow * j] = 0;
                eur_ptr[i + nrow * j] = 0;
            }
        }
    }
    
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

static SEXP read_bed_file_c(SEXP bed_path, SEXP bim_path, SEXP fam_path, SEXP sample_indices) {
    const char *bed = CHAR(STRING_ELT(bed_path, 0));
    const char *bim = CHAR(STRING_ELT(bim_path, 0));
    const char *fam = CHAR(STRING_ELT(fam_path, 0));
    
    FILE *fp = fopen(bed, "rb");
    if (!fp) error("Cannot open bed file: %s", bed);
    
    unsigned char magic[3];
    fread(magic, 1, 3, fp);
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
        fread(buffer, 1, store_size, fp);
        
        for (int s = 0; s < n_samples; s++) {
            int byte_idx = s / 4;
            int bit_idx = (s % 4) * 2;
            int genotype = (buffer[byte_idx] >> bit_idx) & 0x03;
            
            if (genotype == 3) gt[s + n_samples * v] = 3;
            else if (genotype == 2) gt[s + n_samples * v] = 2;
            else if (genotype == 1) gt[s + n_samples * v] = 1;
            else gt[s + n_samples * v] = 0;
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
            
            if (an == 3) {
                sprintf(gt_afr, "%d/%d", gt, gt);
                sprintf(gt_eur, "0/0");
            } else if (an == 1) {
                sprintf(gt_afr, "0/0");
                sprintf(gt_eur, "%d/%d", gt, gt);
            } else if (an == 2) {
                int a1 = gt > 0 ? 1 : 0;
                int a2 = gt > 1 ? 1 : 0;
                sprintf(gt_afr, "%d/%d", a1, a2);
                sprintf(gt_eur, "%d/%d", gt - a1, gt - a2);
            } else {
                sprintf(gt_afr, "0/0");
                sprintf(gt_eur, "0/0");
            }
            
            double ds_afr = (an == 3) ? gt : ((an == 2) ? gt * 0.5 : 0);
            double ds_eur = (an == 1) ? gt : ((an == 2) ? gt * 0.5 : 0);
            
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

/* Exported functions - called from init.c registration */
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
