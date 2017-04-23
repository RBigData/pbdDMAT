/* Automatically generated. Do not edit by hand. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP R_convert_csr_to_dense(SEXP dim, SEXP data, SEXP row_ptr, SEXP col_ind);
extern SEXP R_convert_dense_to_csr(SEXP x);
extern SEXP R_int_sparse_count_zeros(SEXP x);
extern SEXP R_sparse_count_zeros(SEXP x, SEXP tol);

static const R_CallMethodDef CallEntries[] = {
  {"R_convert_csr_to_dense", (DL_FUNC) &R_convert_csr_to_dense, 4},
  {"R_convert_dense_to_csr", (DL_FUNC) &R_convert_dense_to_csr, 1},
  {"R_int_sparse_count_zeros", (DL_FUNC) &R_int_sparse_count_zeros, 1},
  {"R_sparse_count_zeros", (DL_FUNC) &R_sparse_count_zeros, 2},
  {NULL, NULL, 0}
};
void R_init_pbdDMAT(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
