#ifndef __DMAT_PACKAGE__
#define __DMAT_PACKAGE__


#include <R.h>
#include <Rinternals.h>
#include <RNACI.h>


// converters.c


// sparse_utils.c
int int_sparse_count_zeros(int m, int n, int *x);
SEXP R_int_sparse_count_zeros(SEXP x);
int sparse_count_zeros(int m, int n, double *x, double tol);
SEXP R_sparse_count_zeros(SEXP x, SEXP tol);
int sparse_count_zeros_withrows(int m, int n, int *rows, double *x);


#endif

