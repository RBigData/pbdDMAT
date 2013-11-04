#ifndef DMAT_MATEXP_H
#define DMAT_MATEXP_H


void dmat_matexp_pade(const unsigned int n, double *A, double *N, double *D);
void matpow_by_squaring(double *A, int n, int b, double *P);


#endif
