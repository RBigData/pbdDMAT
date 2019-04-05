#ifndef __DMAT_PACKAGE__
#define __DMAT_PACKAGE__


#include <R.h>
#include <Rinternals.h>
#include <RNACI.h>

int sparse_count_zeros_withrows(const int m, const int n, int *const restrict rows, const double *const restrict x);


#endif
