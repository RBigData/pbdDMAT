#include "dmat.h"
#include "matexp.h"

SEXP R_matexp_pade(SEXP n, SEXP A)
{
  SEXP N, D;
  SEXP RET, RET_NAMES;
  
  // Allocate N and D
  PROTECT(N = allocMatrix(REALSXP, INT(n,0), INT(n,0)));
  PROTECT(D = allocMatrix(REALSXP, INT(n,0), INT(n,0)));
  
  // Compute N and D
  matexp_pade(INT(n,0), REAL(A), REAL(N), REAL(D));
  
  // Wrangle the return
  PROTECT(RET = allocVector(VECSXP, 2));
  PROTECT(RET_NAMES = allocVector(STRSXP, 2));
  
  SET_VECTOR_ELT(RET, 0, N);
  SET_VECTOR_ELT(RET, 1, D);
  
  SET_STRING_ELT(RET_NAMES, 0, mkChar("N")); 
  SET_STRING_ELT(RET_NAMES, 1, mkChar("D")); 
  
  setAttrib(RET, R_NamesSymbol, RET_NAMES);
  
  
  UNPROTECT(4);
  return(RET);
}


SEXP R_matpow_by_squaring(SEXP n, SEXP A, SEXP b)
{
  double *cpA;
  const int N = INT(n,0);
  
  SEXP P;
  PROTECT(P = allocMatrix(REALSXP, N, N));
  
/*  cpA = R_alloc(N*N, sizeof(double));*/
  cpA = malloc(N*N*sizeof(double));
  memcpy(cpA, REAL(A), N*N*sizeof(double));
  
  matpow_by_squaring(cpA, INT(n,0), INT(b,0), REAL(P));
  
  free(cpA);
  
  UNPROTECT(1);
  return(P);
}

