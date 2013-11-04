## autoconf configure.ac && autoreconf -if

AC_PREREQ([2.69])
AC_INIT([pbdDMAT], [0.2-2])
AC_CONFIG_SRCDIR([src/matexp.c])

# Checks for programs.
AC_PROG_CC

AC_CHECK_FUNCS([sqrt])


# Determine R_HOME.
: ${R_HOME=‘R RHOME‘}
echo ${R_HOME}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

# Obtain from R.
R_SCMD="${R_HOME}/bin/R CMD config"
FLIBS=`${R_SCMD} FLIBS`
BLAS_LIBS=`${R_SCMD} BLAS_LIBS`
LAPACK_LIBS=`${R_SCMD} LAPACK_LIBS`

# Check tools.
AC_CHECK_PROG([ls_ok], [ls], [yes], [no], ,)
AC_CHECK_PROG([echo_ok], [echo], [yes], [no], ,)
AC_CHECK_PROG([sed_ok], [sed], [yes], [no], ,)
AC_CHECK_PROG([grep_ok], [grep], [yes], [no], ,)

# Check Rlapack.
if test "X$echo_ok" = "Xyes" -a "X$grep_ok" = "Xyes"; then
  Rlapack_ok=`echo ${LAPACK_LIBS} | grep "\-lRlapack"`
  if test "X$Rlapack_ok" = "X"; then
    # Rlapack is not found.
    Rlapack_ok="NA"
  else
    echo "Found Rlpack from ${Rlapack_ok}"
  fi
else
  # echo and grep are not found.
  Rlapack_ok="NA"
fi

if test "X$Rlapack_ok" = "XNA"; then
  # An external LAPACK is used in run time.
  QZ_OBJS=""
else
  if test "X$ls_ok" = "Xyes" -a "X$sed_ok" = "Xyes" -a "X$grep_ok" = "Xyes"; then
    ### Check R_ext/Lapack.h
    LAPACK_HEADER="${R_HOME}/include/R_ext/Lapack.h"

    funcs=`ls src/qz/*.f | sed -e "s/.*\/\(.*\)\.f/\\1/"`
    QZ_OBJS=""
    for func in ${funcs}
    do
      func_ok=`grep "F77_NAME($func)" ${LAPACK_HEADER}`
      echo "Found $func from ${func_ok}"
      if test "X${func_ok}" = "X"; then
        QZ_OBJS="${QZ_OBJS} qz/$func.o"
      fi
    done
  else
    ### Use default if all above fail.
    QZ_OBJS="\$(QZ_OBJS)"
  fi
fi




AC_OUTPUT(src/Makevars)


