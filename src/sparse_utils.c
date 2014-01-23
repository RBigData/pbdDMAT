#include <math.h>
#include "dmat.h"


int int_sparse_count_zeros(int m, int n, int *x)
{
  int count = 0;
  int i, j;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (x[i + m*j] == 0)
        count++;
    }
  }
  
  return count;
}



SEXP R_int_sparse_count_zeros(SEXP x)
{
  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 1));
  
  INTEGER(ret)[0] = int_sparse_count_zeros(nrows(x), ncols(x), INTEGER(x));
  
  UNPROTECT(1);
  return ret;
}



int sparse_count_zeros(int m, int n, double *x, double tol)
{
  int count = 0;
  int i, j;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (fequals(x[i + m*j], 0.0, tol))
        count++;
    }
  }
  
  return count;
}



SEXP R_sparse_count_zeros(SEXP x, SEXP tol)
{
  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 1));
  
  INTEGER(ret)[0] = sparse_count_zeros(nrows(x), ncols(x), REAL(x), REAL(tol)[0]);
  
  UNPROTECT(1);
  return ret;
}
