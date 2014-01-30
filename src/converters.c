// Copyright 2013, Schmidt


#include "dmat.h"


SEXP sbase_convert_csr_to_dense(SEXP dim, SEXP data, SEXP row_ptr, SEXP col_ind)
{
  R_INIT;
  int i, j;
  int c = 0, r = 0, diff;
  const int m = INT(dim, 0), n = INT(dim, 1);
  SEXP dense_mat;
  
  
  // Initialize
  newRmat(dense_mat, m, n, "dbl");
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
      MatDBL(dense_mat, i, j) = 0.0;
  }
  
  
  // This is disgusting
  i = 0;
  while (r < m && INT(row_ptr, r) < INT(row_ptr, m))
  {
    diff = INT(row_ptr, r+1) - INT(row_ptr, r);
    
    if (diff == 0)
      goto increment;
    else
    {
      while (diff)
      {
        j = INT(col_ind, c)-1;
        MatDBL(dense_mat, i, j) = DBL(data, c);
        
        c++; // hehehe
        diff--;
      }
    }
    
    increment:
      r++;
      i++;
  }
  
  
  R_END;
  return dense_mat;
}



SEXP sbase_convert_dense_to_csr(SEXP x)
{
  R_INIT;
  SEXP data, row_ptr, col_ind;
  SEXP R_list, R_list_names;
  const int m = nrows(x), n = ncols(x);
  int i, j;
  int row_ptr_len;
  int sparsity, density;
  int ct = 0, rct = 0, first;
  
  
  sparsity = sparse_count_zeros_withrows(m, n, &row_ptr_len, REAL(x));
  density = m*n - sparsity;
  
  newRvec(data, density, "dbl");
  newRvec(col_ind, density, "int");
  newRvec(row_ptr, m+1, "int");
  
  
  for (i=0; i<m; i++)
  {
    first = true;
    
    for (j=0; j<n; j++)
    {
      if (MatDBL(x, i, j) > 0.0)
      {
        DBL(data, ct) = MatDBL(x, i, j);
        INT(col_ind, ct) = j+1;
        ct++;
        
        if (first == true)
        {
          INT(row_ptr, rct) = ct;
          first = false;
          rct++;
        }
      }
    }
    
    if (first == true)
    {
      INT(row_ptr, rct) = ct+1;
      rct++;
    }
  }
  
  INT(row_ptr, m) = ct+1;
  
  R_list_names = make_list_names(3, "Data", "row_ptr", "col_ind");
  R_list = make_list(R_list_names, data, row_ptr, col_ind);
  
  R_END;
  return R_list;
}




