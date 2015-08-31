#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([1.0])
AC_INIT([libAMPI], [1.0], [schanen@stce.rwth-aachen.de])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror])
m4_pattern_allow([AM_PROG_AR])
#AM_PROG_AR([act-if-fail])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP
AC_PROG_RANLIB
AC_LANG([C])
#AM_PROG_AR


AC_CHECK_LIB([m], [sqrt],,AC_MSG_ERROR(Could not find -lm libary.))
AC_LANG([C])
AC_LANG([C++])

#if test "X${MPI_ENABLED}" = "Xyes"; then
    echo ""
    echo "MPI-C Settings"
    echo "--------------"
    echo ""
    AMPI_SET_MPICC
#fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_SUBST(MPICC)
AC_SUBST(MPI_INC_DIR)
AC_SUBST(MPI_LIB_DIR)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_FLAGS)

AC_CANONICAL_BUILD
AC_CANONICAL_HOST

AC_OUTPUT
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN],
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN],
[AC_CONFIG_FILES([doc/Doxyfile])])
AMPI_REPORT
