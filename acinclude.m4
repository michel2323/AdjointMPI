
#=================================================================#
#                                                                 #
#                                                                 #
#           P A R A L L E L     S U P P O R T                     #
#                                                                 #
#                                                                 #
#==================================================================


#------------------------------------------------------------------
# CHECK MPI-C COMPILER
#------------------------------------------------------------------

AC_DEFUN([AMPI_SET_MPICC],
[

AC_ARG_WITH([],[    ],[])

# MPI root directory
AC_ARG_WITH(mpi-root,
[AC_HELP_STRING([--with-mpi-root=MPIROOT],[use MPI root directory])],
[
MPI_ROOT_DIR="${withval}"
],
[
MPI_ROOT_DIR=""
])

# MPI include directory
AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
[
MPI_INC_DIR="${withval}"
],
[
MPI_INC_DIR=""
])

# MPI library directory
AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
[
MPI_LIB_DIR="${withval}"
],
[
MPI_LIB_DIR=""
])

# MPI libraries
AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs=ARG],[MPI libraries])],
[
MPI_LIBS="${withval}"
],
[
MPI_LIBS=""
])

# MPI flags
AC_ARG_WITH(mpi-flags,
[AC_HELP_STRING([--with-mpi-flags=ARG],[MPI-specific flags])],
[
MPI_FLAGS="${withval}"
MPI_FLAGS_OK="yes"
],
[
MPI_FLAGS=""
MPI_FLAGS_OK="no"
])

# MPI-C compiler
MPICC_COMP_GIVEN="yes"
AC_MSG_CHECKING([if using MPI-C script])
AC_ARG_WITH(mpicc,
[AC_HELP_STRING([--with-mpicc[[[[=ARG]]]]],[specify MPI-C compiler to use @<:@mpicc@:>@])],
[
if test "X${withval}" = "Xno"; then
  USE_MPICC_SCRIPT="no"
else
  USE_MPICC_SCRIPT="yes"
  MPICC_COMP="${withval}"
fi
],
[
  USE_MPICC_SCRIPT="yes"
  MPICC_COMP="mpicc"
  MPICC_COMP_GIVEN="no"
])
AC_MSG_RESULT([${USE_MPICC_SCRIPT}])

# If CC is a C++ compiler, then we certainly do NOT want to use an MPI-C script
# Note: USING_CPLUSPLUS_COMP was defined by a call to AMPI_CPLUSPLUS_CHECK
# in AMPI_SET_CC
# Note: If the user specified an MPI-C script, then we will NOT do anything for now
if test "X${MPICC_COMP_GIVEN}" = "Xno" && test "X${USING_CPLUSPLUS_COMP}" = "Xyes"; then
  MPICC_COMP="mpiCC"
fi

# Check MPI-C compiler (either MPI compiler script or regular C compiler)
if test "X${USE_MPICC_SCRIPT}" = "Xyes"; then
  AMPI_CHECK_MPICC
else
  MPICC_COMP="${CC}"
  MPICC="${CC}"
  AMPI_CC_WITH_MPI_CHECK
fi

]) dnl END AMPI_SET_MPICC

#------------------------------------------------------------------
# TEST MPI-C COMPILER
#------------------------------------------------------------------

AC_DEFUN([AMPI_CHECK_MPICC],
[

# Test MPI-C compiler (meaning test MPICC_COMP)
# Check if MPI-C compiler can be found

AC_MSG_CHECKING([if absolute path to ${MPICC_COMP} was given])

# CASE 1: MPICC_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPICC_COMP} ; then

  AC_MSG_RESULT([yes])
  MPICC_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPICC_COMP}"])`
  TMP_MPI_INC_DIR="${MPI_BASE_DIR}/../include"
  TMP_MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"

# CASE 2: MPICC_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else

  AC_MSG_RESULT([no])

  if test "X${MPI_ROOT_DIR}" = "X"; then
    # Try to find location of executable (perhaps directory was entered
    # incorrectly)
    TEMP_MPICC_COMP=`basename "${MPICC_COMP}"`
    AC_PATH_PROG([MPICC_COMP],[${TEMP_MPICC_COMP}],[none])
    # Cannot find executable in PATH
    if test "X${MPICC_COMP}" = "Xnone"; then
      MPICC_COMP_EXISTS="no"
      MPICC_COMP=""
    # Found executable and set MPICC_COMP to absolute pathname
    else
      MPICC_COMP_EXISTS="yes"
      MPI_BASE_DIR=`AS_DIRNAME(["${MPICC_COMP}"])`
      TMP_MPI_INC_DIR="${MPI_BASE_DIR}/../include"
      TMP_MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
    fi

  # CASE 3: MPICC_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else

    AC_MSG_CHECKING([if ${MPICC_COMP} exists in ${MPI_ROOT_DIR}/bin])
    # MPICC_COMP should really only contain an executable name
    # Found location of MPICC_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPICC_COMP} ; then
      AC_MSG_RESULT([yes])
      MPICC_COMP_EXISTS="yes"
      MPICC_COMP="${MPI_ROOT_DIR}/bin/${MPICC_COMP}"
      TMP_MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      TMP_MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    # Could NOT find MPICC_COMP anywhere
    else
      AC_MSG_RESULT([no])
      MPICC_COMP_EXISTS="no"
      MPICC_COMP=""
    fi

  fi

fi

# If MPICC_COMP exists, set MPICC and (conditionally) set MPI_INC_DIR
# and MPI_LIB_DIR so that we do not end up with empty -I options.
# Otherwise, issue warning message
if test "X${MPICC_COMP_EXISTS}" = "Xyes"; then

  MPICC="${MPICC_COMP}"
  MPI_C_COMP_OK="yes"

  # If MPI_INC_DIR is empty, set it to TMP_MPI_INC_DIR
  if test "X${MPI_INC_DIR}" = "X"; then
    MPI_INC_DIR="$TMP_MPI_INC_DIR"
  fi

  # If MPI_LIB_DIR is empty, set it to TMP_MPI_LIB_DIR
  if test "X${MPI_LIB_DIR}" = "X"; then
    MPI_LIB_DIR="$TMP_MPI_LIB_DIR"
  fi

else

  AC_MSG_WARN([cannot find MPI-C compiler])
  echo ""
  echo "   Unable to find a functional MPI-C compiler."
  echo ""
  echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
  echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
  echo "   to specify the locations of all relevant MPI files, or"
  echo "   --with-mpi-root to specify the base installation directory"
  echo "   of the MPI implementation to be used."
  echo ""
  echo "   Disabling the parallel NVECTOR module and all parallel examples..."
  echo ""
  MPICC=""
  MPI_C_COMP_OK="no"
  AMPI_WARN_FLAG="yes"

fi

]) dnl END AMPI_CHECK_MPICC

#------------------------------------------------------------------
# TEST C COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([AMPI_CC_WITH_MPI_CHECK],
[

# Test if we can compile MPI programs using the CC compiler
# and current MPI settings

AC_MSG_NOTICE([Testing CC with MPI settings])

# Save copies of CPPFLAGS, LDFLAGS and LIBS (preserve information)
# Temporarily overwritten so we can test MPI implementation
SAVED_CPPFLAGS="${CPPFLAGS}"
SAVED_LDFLAGS="${LDFLAGS}"
SAVED_LIBS="${LIBS}"

# Determine location of MPI header files (find MPI include directory)
MPI_EXISTS="yes"

AC_MSG_CHECKING([for location of MPI implementation])

# If MPI include directory was NOT explicitly specified, check if MPI root
# directory was given by user
if test "X${MPI_INC_DIR}" = "X"; then
  # If MPI root directory was NOT given so issue a warning message
  if test "X${MPI_ROOT_DIR}" = "X"; then
    AC_MSG_RESULT([not found])
    MPI_EXISTS="no"
    AC_MSG_WARN([cannot find MPI implementation files])
    echo ""
    echo "   Unable to find MPI implementation files."
    echo ""
    echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
    echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
    echo "   to specify the locations of all relevant MPI files, or"
    echo "   --with-mpi-root to specify the base installation directory"
    echo "   of the MPI implementation to be used."
    echo ""
    echo "   Disabling the parallel NVECTOR module and all parallel examples..."
    echo ""
    AMPI_WARN_FLAG="yes"
  # MPI root directory was given so set MPI_INC_DIR accordingly
  # Update CPPFLAGS
  else
    MPI_INC_DIR="${MPI_ROOT_DIR}/include"
    AC_MSG_RESULT([${MPI_INC_DIR}])
    if test "X${CPPFLAGS}" = "X"; then
      CPPFLAGS="-I${MPI_INC_DIR}"
    else
      CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
    fi
    # Add MPI_FLAGS if non-empty
    if test "X${MPI_FLAGS}" = "X"; then
      CPPFLAGS="${CPPFLAGS}"
    else
      CPPFLAGS="${CPPFLAGS} ${MPI_FLAGS}"
    fi
  fi
# MPI include directory was specified so update CPPFLAGS
else
  AC_MSG_RESULT([${MPI_INC_DIR}])
  if test "X${CPPFLAGS}" = "X"; then
    CPPFLAGS="-I${MPI_INC_DIR}"
  else
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
  fi
  # Add MPI_FLAGS if non-empty
  if test "X${MPI_FLAGS}" = "X"; then
    CPPFLAGS="${CPPFLAGS}"
  else
    CPPFLAGS="${CPPFLAGS} ${MPI_FLAGS}"
  fi
fi

# Only continue if found an MPI implementation
if test "X${MPI_EXISTS}" = "Xyes"; then

  AC_MSG_CHECKING([for location of MPI libraries])

  # Determine location of MPI libraries
  # MPI library directory was NOT specified by user so set based upon MPI_ROOT_DIR
  # Update LDFLAGS
  if test "X${MPI_LIB_DIR}" = "X"; then
    MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    AC_MSG_RESULT([${MPI_LIB_DIR}])
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  # MPI library directory was specified so update LDFLAGS
  else
    AC_MSG_RESULT([${MPI_LIB_DIR}])
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  fi

  # Check if user specified which MPI libraries must be included
  # If no libraries are given, then issue a warning message
  AC_MSG_CHECKING([for MPI libraries])
  if test "X${MPI_LIBS}" = "X"; then
    AC_MSG_RESULT([none])
    AC_MSG_WARN([no MPI libraries were given])
    echo ""
    echo "   Unable to compile MPI program using C compiler because"
    echo "   MPI libraries were not specified."
    echo ""
    echo "   Try using --with-mpi-libdir and --with-mpi-libs to"
    echo "   specify the location and names of the MPI libraries."
    echo ""
    echo "   Disabling the parallel NVECTOR module and all parallel examples..."
    echo ""
    MPI_C_COMP_OK="no"
    AMPI_WARN_FLAG="yes"
  # MPI libraries were specified so update LIBS
  else
    AC_MSG_RESULT([${MPI_LIBS}])
    if test "X${LIBS}" = "X"; then
      LIBS="${MPI_LIBS}"
    else
      LIBS="${LIBS} ${MPI_LIBS}"
    fi
    # Set the MPI_C_COMP_OK variable to NULL so we can conditionally execute
    # the next test
    MPI_C_COMP_OK=""
  fi

  if test "X${MPI_C_COMP_OK}" = "X"; then
    AC_MSG_CHECKING([if C compiler can compile MPI programs])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include "mpi.h"]],[[int c; char **v; MPI_Init(&c,&v);]])],
    [AC_MSG_RESULT([yes])
     MPI_C_COMP_OK="yes"],
    [AC_MSG_RESULT([no])
     AC_MSG_WARN([C compiler cannot compile MPI programs])
     echo ""
     echo "   Unable to compile MPI program using C compiler."
     echo ""
     echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
     echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
     echo "   to specify the locations of all relevant MPI files, or"
     echo "   --with-mpi-root to specify the base installation directory"
     echo "   of the MPI implementation to be used."
     echo ""
     echo "   Disabling the parallel NVECTOR module and all parallel examples..."
     echo ""
     MPI_C_COMP_OK="no"
     AMPI_WARN_FLAG="yes"])
  fi
else
  MPI_C_COMP_OK="no"
fi
  
# Restore CPPFLAGS, LDFLAGS and LIBS
CPPFLAGS="${SAVED_CPPFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

]) dnl END AMPI_CC_WITH_MPI_CHECK

#------------------------------------------------------------------
# PRINT STATUS REPORT
#------------------------------------------------------------------

AC_DEFUN([AMPI_REPORT],
[

if test "X${AMPI_WARN_FLAG}" = "Xyes"; then
echo "
***************
*   WARNING   *
***************

At least one warning was issued. Some features were disabled.

Review the configure output and/or the contents of config.log 
before proceeding with the build.
"
fi

echo "
------------------------------
AMPI Configuration Summary
------------------------------"

echo "
Configuration
-------------

  Host System:               ${host}
  Build System:              ${build}

  C Preprocessor:            ${CPP} 
  C Preprocessor Flags:      ${CPPFLAGS}
  C Compiler:	             ${CC}
  C Compiler Flags           ${CFLAGS}
  C Linker:                  ${CC}
  Linker Flags:              ${LDFLAGS}
  Libraries:                 ${LIBS}"


echo "
  MPI Root Directory:        ${MPI_ROOT_DIR}
  MPI Include Directory:     ${MPI_INC_DIR}
  MPI Library Directory:     ${MPI_LIB_DIR}
  MPI Flags:                 ${MPI_FLAGS}
  Extra MPI Libraries:       ${MPI_LIBS}

  Using MPI-C script?        ${USE_MPICC_SCRIPT}
  MPI-C:                     ${MPICC}"


# Determine SOURCE, BUILD, and EXEC_PREFIX directories
cv_srcdir=`( cd ${srcdir} ; pwd )`
cv_builddir=`pwd`
if test "X${exec_prefix}" = "XNONE"; then
  cv_exec_prefix="${prefix}"
else
  cv_exec_prefix="${exec_prefix}"
fi

echo "
  srcdir:                    ${cv_srcdir}
  builddir:                  ${cv_builddir}
  prefix:                    ${prefix}
  exec_prefix:               ${cv_exec_prefix}
  includedir:                ${includedir}
  libdir:                    ${libdir}"


echo "  
  Type 'make' and then 'make install' to build and install ${PACKAGE_STRING}."



echo "
----------------------------------
Finished AMPI Configure Script
----------------------------------
"

]) dnl END AMPI_REPORT
