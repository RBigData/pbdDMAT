// Copyright 2013, Schmidt


#include <SEXPtools.h>
#include "dmat.h"


SEXP sbase_convert_csr_to_dense(SEXP dim, SEXP data, SEXP row_ptr, SEXP col_ind)
{
  R_INIT;
  int i, j;
  const int m = INT(dim, 0), n = INT(dim, 1);
  SEXP dense_mat;
  
  newRmat(dense_mat, m, n, "dbl");
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
      MatDBL(dense_mat, i, j) = 0.0;
  }
  
  
  
  
  
  R_END;
  return dense_mat;
}




SEXP sbase_convert_dense_to_csr(SEXP x)
{
  R_INIT;
  SEXP data, row_ptr, col_ind;
  const int m = nrows(x), n = ncols(x);
  int i, j;
  int sparsity, density;
  int ct = 0, first;
  
/*  sparsity = sparse_count_zeros(m, n, REAL(x), tol);*/ //FIXME
  density = m*n - sparsity;
  
  newRvec(data, density, "dbl");
  newRvec(col_ind, density, "int");
  newRvec(row_ptr, m+1, "int");
  
  for (j=0; j<n; j++)
  {
    first = 0;
    
    for (i=0; i<m; i++)
    {
      if (MatDBL(x, i, j) > 0.0)
      {
        DBL(data, ct) = MatDBL(x, i, j);
        DBL(col_ind, ct) = j+1;
        ct++;
        
        if (first == 0)
        {
          INT(row_ptr, i) = i + n*j; //FIXME
        }
      }
    }
  }
  
  
  R_END;
  return RNULL; // FIXME
}




