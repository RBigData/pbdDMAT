#include <math.h>
#include <stdbool.h>

#include <R.h>
#include <Rinternals.h>


#define IS_ZERO(x, tol) (fabs(x) < tol)
#define TOL 1e-10


static int int_sparse_count_zeros(int m, int n, int *x)
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



static int sparse_count_zeros(int m, int n, double *x, double tol)
{
  int count = 0;
  int i, j;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (IS_ZERO(x[i + m*j], tol))
        count++;
    }
  }
  
  return count;
}



int sparse_count_zeros_withrows(const int m, const int n, int *const restrict rows, const double *const restrict x)
{
  int count = 0;
  int i, j;
  int first;
  
  *rows = 0;
  
  for (i=0; i<m; i++)
  {
    first = true;
    for (j=0; j<n; j++)
    {
      if (IS_ZERO(x[i + m*j], TOL))
      {
        count++;
      
        if (first == true)
        {
          (*rows)++;
          first = false;
        }
      }
    }
  }
  
  return count;
}


//FIXME remove tol
SEXP R_sparse_count_zeros(SEXP x, SEXP tol)
{
  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 1));
  
  INTEGER(ret)[0] = sparse_count_zeros(nrows(x), ncols(x), REAL(x), REAL(tol)[0]);
  
  UNPROTECT(1);
  return ret;
}
