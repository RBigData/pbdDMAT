#ifndef __DMAT_PACKAGE__
#define __DMAT_PACKAGE__


#include <R.h>
#include <Rinternals.h>

#define false 0
#define true 1

#define fequals(x,y,tol) (fabs(x-y)<tol?true:false)


// converters.c


// sparse_utils.c
int int_sparse_count_zeros(int m, int n, int *x);
SEXP R_int_sparse_count_zeros(SEXP x);
int sparse_count_zeros(int m, int n, double *x, double tol);
SEXP R_sparse_count_zeros(SEXP x, SEXP tol);


#endif

