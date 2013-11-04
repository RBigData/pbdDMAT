#ifndef R_DMAT_H
#define R_DMAT_H


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

// character pointers
#define CHARPT(x,i) ((char*)CHAR(STRING_ELT(x,i)))

#define INT(x,i) (INTEGER(x)[i])
#define RL(x,i) (REAL(x)[i])


#endif
