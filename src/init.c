#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "ancestry.h"

static const R_CallMethodDef CallEntries[] = {
    {"count_ancestry_codes_C", (DL_FUNC) &count_ancestry_codes, 2},
    {"split_by_ancestry_C", (DL_FUNC) &split_by_ancestry, 2},
    {"read_bed_file_C", (DL_FUNC) &read_bed_file, 4},
    {"write_vcf_with_ancestry_C", (DL_FUNC) &write_vcf_with_ancestry, 5},
    {"subset_vcf_by_range_C", (DL_FUNC) &subset_vcf_by_range, 5},
    {NULL, NULL, 0}
};

void R_init_lantern(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
