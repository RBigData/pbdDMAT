dnl Process this file with autoconf to produce a configure script
AC_INIT(DESCRIPTION)

R_SCMD="${R_HOME}/bin/Rscript -e"

dnl Rcpp
RCPP_CXXFLAGS=`${R_SCMD} "Rcpp:::CxxFlags()"`
RCPP_LDFLAGS=`${R_SCMD} "Rcpp:::LdFlags()"`

dnl Start to substitute templates
AC_SUBST(RCPP_CXXFLAGS)
AC_SUBST(RCPP_LDFLAGS)
AC_OUTPUT(src/Makevars)
