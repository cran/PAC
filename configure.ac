AC_INIT([PAC], [0.0.1])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "PAC - configure: could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

AC_PROG_CC
AC_PROG_CXX

AC_LANG(C++)


PKG_CPPFLAGS="-g"
if test "x$GCC" = xyes; then
  PKG_CFLAGS="-ffloat-store -g"
fi

AC_SUBST(PKG_CFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
