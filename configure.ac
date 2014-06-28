## autoconf configure.ac && autoreconf -if

AC_PREREQ([2.69])
AC_INIT([pbdDMAT], [0.2-4])
AC_CONFIG_SRCDIR([src/sparse_utils.c])

# Checks for programs.
AC_PROG_CC

R_SCMD="${R_HOME}/bin/Rscript -e"

# RNACI
RNACI_CPPFLAGS=`${R_SCMD} "RNACI:::cppflags()"`
RNACI_LDFLAGS=`${R_SCMD} "RNACI:::ldflags()"`


# Substitute templates
AC_SUBST(RNACI_CPPFLAGS)
AC_SUBST(RNACI_LDFLAGS)

AC_OUTPUT(src/Makevars)
