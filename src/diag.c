#include <R.h>
#include <Rinternals.h>


static inline int indxl2g(const int indxloc, const int nb, const int iproc, const int isrcproc, const int nprocs)
{
  return nprocs*nb*((indxloc - 1)/nb) + ((indxloc - 1) % nb) + ((nprocs + iproc - isrcproc) % nprocs)*nb + 1;
}

// grid: [nprow, npcol, ictxt, myrow, mycol]
SEXP R_diag_set(SEXP x_, SEXP val_, SEXP bldim, SEXP grid)
{
  const int m = nrows(x_);
  const int n = ncols(x_);
  const double *x = REAL(x_);
  
  const double val = REAL(val_)[0];
  
  const int mb = INTEGER(bldim)[0];
  const int nb = INTEGER(bldim)[1];
  
  const int nprow = INTEGER(grid)[0];
  const int npcol = INTEGER(grid)[1];
  const int myrow = INTEGER(grid)[3];
  const int mycol = INTEGER(grid)[4];
  
  SEXP d_;
  PROTECT(d_ = allocMatrix(REALSXP, m, n));
  double *d = REAL(d_);
  
  int gi, gj;
  
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<m; i++)
    {
      gi = indxl2g(i, mb, myrow, 0, nprow);
      gj = indxl2g(j, nb, mycol, 0, npcol);
      
      if (gi == gj)
        d[i + m*j] = val;
      else
        d[i + m*j] = x[i + m*j];
    }
  }
  
  UNPROTECT(1);
  return d_;
}
