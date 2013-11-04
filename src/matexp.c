  /* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <stdio.h>
#include <stdlib.h>

#define EXPSGN(x,pow) (x==0?(pow==0?1:0):(x>0?1:(pow%2==0?1:(-1))))


// C = A * B
static inline void matprod(const unsigned int n, double *A, double *B, double *C)
{
  char trans = 'N';
  double one = 1.0, zero = 0.0;
  
  dgemm_(&trans, &trans, &n, &n, &n, &one, B, &n, A, &n, &zero, C, &n);
}

// Copy A ONTO B, i.e. B = A
static inline void matcopy(const unsigned int n, double *A, double *B)
{
  int i;
  
  #pragma omp for simd
  {
    for (i=0; i<n*n; i++)
      B[i] = A[i];
  }
}

// Identity matrix
static inline void mateye(const unsigned int n, double *A)
{
  int i, j;
  
  #pragma omp for simd
  {
    for (i=0; i<n*n; i++)
      A[i] = 0.0;
    
    // Fill diagonal with 1
    j = 0;
    for (i=0; i<n*n; i+=n)
    {
      i += j;
      
      A[i] = 1;
      
      j = 1;
    }
  }
}


/* r_m(x) = p_m(x) / q_m(x), where
   p_m(x) = sum_{j=0}^m (2m-j)!m!/(2m)!/(m-j)!/j! * x^j
   
   and q_m(x) = p_m(-x)
*/


// Generated from bc via:
// scale=50
// define f(x) {
// if (x <= 1) return (1);
// return (x*f(x-1));
// }
// m=13
// for(j=0; j<=m; j++) {f(2*m-j)*f(m)/f(2*m)/f(m-j)/f(j);}
const double dmat_pade_coefs[14] = 
{
  1.0000000000000000000,
  0.5,
  0.12,
  1.833333333333333333333e-2,
  1.992753623188405797101e-3,
  1.630434782608695652174e-4,
  1.035196687370600414079e-5,
  5.175983436853002070393e-7,
  2.043151356652500817261e-8,
  6.306022705717595115002e-10,
  1.483770048404140027059e-11,
  2.529153491597965955215e-13,
  2.810170546219962172461e-15,
  1.544049750670308885967e-17
};



// Matrix exponentiation using Pade' approximations
// p==q==13
void matexp_pade(const unsigned int n, double *A, double *N, double *D)
{
  int i, j;
  int itmp;
  double tmp, tmpj;
  double *B, *C;
  
  B = malloc(n*n*sizeof(double));
  C = malloc(n*n*sizeof(double));
  
  // Initialize
  #pragma omp for simd
  {
    for (i=0; i<n*n; i++)
    {
      N[i] = 0.0;
      D[i] = 0.0;
    }
    
    // Fill diagonal with 1
    j = 0;
    for (i=0; i<n*n; i+=n)
    {
      i += j;
      
      N[i] = 1;
      D[i] = 1;
      
      j = 1;
    }
  }
  
  // Fill N and D
  for (i=1; i<=13; i++)
  {
    // C = A*B
    if (i>1)
      matprod(n, A, B, C);
    else
    {
      #pragma omp for simd
      {
        for (j=0; j<n*n; j++)
          C[j] = A[j];
      }
    }
    
    #pragma omp for simd
    {
      // B = C
      for (j=0; j<n*n; j++)
        B[j] = C[j];
      
      // N = pade_coef[i] * C
      // D = (-1)^j * pade_coef[i] * C
      tmp = dmat_pade_coefs[i];
      itmp = EXPSGN(-1, i);
      
      if (itmp == 1)
      {
        for (j=0; j<n*n; j++)
        {
          tmpj = tmp * C[j];
          N[j] += tmpj;
          D[j] += tmpj;
        }
      }
      else
      {
        for (j=0; j<n*n; j++)
        {
          tmpj = tmp * C[j];
          N[j] += tmpj;
          D[j] -= tmpj;
        }
      }
    }
  }
  
  free(B);
  free(C);
}



// Exponentiation by squaring
// P = A^b
void matpow_by_squaring(double *A, int n, int b, double *P)
{
  int i, j;
  int itmp;
  double tmp, tmpj;
  double *TMP;
  
  mateye(n, P);
  
  
  // Trivial cases
  if (b == 0)
    return;
  
  if (b == 1)
  {
    matcopy(n, A, P);
    return;
  }
  
  
  // General case
  TMP = malloc(n*n*sizeof(double));
  
  while (b)
  {
    if (b&1)
    {
      matprod(n, P, A, TMP);
      matcopy(n, TMP, P);
    }
    
    b >>=1;
    matprod(n, A, A, TMP);
    matcopy(n, TMP, A);
  }
  
  free(TMP);
}

