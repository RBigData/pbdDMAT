// Copyright 2013, Schmidt

#include "dmat.h"

double get_machine_eps()
{
  R_INIT;
  SEXP dmatPackage;
  SEXP tmp;
  
  PT(
    dmatPackage = eval( lang2( install("getNamespace"), ScalarString(mkChar("pbdDMAT")) ), R_GlobalEnv )
  );
  
  tmp = eval( lang1( install("get_machine_eps")), dmatPackage);
  
  R_END;
  return DBL(tmp,0);
}

