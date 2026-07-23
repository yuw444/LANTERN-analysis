#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "ancestry.h"

static const R_CallMethodDef CallEntries[] = {
    {"count_ancestry_codes_C",      (DL_FUNC) &count_ancestry_codes,      2},
    {"split_by_ancestry_C",         (DL_FUNC) &split_by_ancestry,         4},
    {"split_phased_by_ancestry_C",  (DL_FUNC) &split_phased_by_ancestry,  5},
    {"split_phased_multi_C",        (DL_FUNC) &split_phased_multi,        5},
    {"split_by_ancestry_multi_C",  (DL_FUNC) &split_by_ancestry_multi,   8},
    {"read_bed_file_C",            (DL_FUNC) &read_bed_file,              3},
    {NULL, NULL, 0}
};

void R_init_lantern(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
