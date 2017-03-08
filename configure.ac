AC_INIT(myconfig, version-0.1)
AC_LANG(C++)
AC_PROG_FC

dnl include m4 macro for netcdf4 check
m4_include([m4/ax_lib_netcdf4.m4])

dnl list of all possible optional components
all_options="netcdf4"

dnl When a component is found, append it here
options=

dnl Force the compiler to run once and do all the basic checks
dnl if you don't, it will do the test on the first invocation
dnl below and so your pretty section titles won't work as well
dnl AC_CHECK_HEADERS(iostream)

dnl Small macro to print out a nice, pretty section title.
define(SECTION_TITLE,
[
	dnl echo >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
	echo '   $1' >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
])


SECTION_TITLE([Checking for netcdf4 libraries])


AX_LIB_NETCDF4()
if test "$with_netcdf4" = "yes"; then
	dnl This add to the NETCDF4_LIBS was needed to compile
	AC_SUBST(NETCDF4_VERSION)
	AC_SUBST(NETCDF4_CC)
	AC_SUBST(NETCDF4_CFLAGS)
	AC_SUBST(NETCDF4_CPPFLAGS)
	AC_SUBST(NETCDF4_LDFLAGS)
	AC_SUBST(NETCDF4_LIBS)
	AC_SUBST(NETCDF4_FC)
	AC_SUBST(NETCDF4_FFLAGS)
	AC_SUBST(NETCDF4_FLIBS)
	AC_DEFINE(HAVE_NETCDF4)
	options="$options netcdf4"
else
	AC_MSG_ERROR([Unable to find NetCDF4, we need NetCDF4 for Fortran (nf-config).])
fi



dnl Now process the options strings. Essentially, we want two lists
dnl one for the options present (which we have) and one for the options
dnl missing (which we don't)

SECTION_TITLE([Configuration results])

echo "Options:" >& AS_MESSAGE_FD
echo "$options" >& AS_MESSAGE_FD
echo >& AS_MESSAGE_FD

echo "Missing options:" >& AS_MESSAGE_FD
echo "$options" "$all_options" | tr ' ' '\n' | sort | uniq -u | tr '\n' ' ' >& AS_MESSAGE_FD



dnl This part will be used to establish R variables

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

FC=`"${R_HOME}/bin/R" CMD config FC`
FCFLAGS=`"${R_HOME}/bin/R" CMD config FCFLAGS`
SHLIB_FFLAGS=`"${R_HOME}/bin/R" CMD config SHLIB_FFLAGS`
FPICFLAGS=`"${R_HOME}/bin/R" CMD config FPICFLAGS`
SAFE_FFLAGS=`"${R_HOME}/bin/R" CMD config SAFE_FFLAGS`


FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
SAFE_FFLAGS=`"${R_HOME}/bin/R" CMD config SAFE_FFLAGS`



echo -e "\n" >& AS_MESSAGE_FD
SECTION_TITLE([Important variables])
if test "$with_netcdf4_fortran" = "yes"; then
	echo "NETCDF4_VERSION=$NETCDF4_VERSION" >& AS_MESSAGE_FD
	echo "NETCDF4_FC=$NETCDF4_FC" >& AS_MESSAGE_FD

	echo "NETCDF4_FLIBS=$NETCDF4_FLIBS" >& AS_MESSAGE_FD
	echo "NETCDF4_FFLAGS=$NETCDF4_FFLAGS" >& AS_MESSAGE_FD
	echo "" >& AS_MESSAGE_FD
	echo "From R session:" >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >& AS_MESSAGE_FD
	echo "FC=$FC" >& AS_MESSAGE_FD
	echo "FCFLAGS=$FCFLAGS" >& AS_MESSAGE_FD
	echo "SHLIB_FFLAGS=$SHLIB_FFLAGS" >& AS_MESSAGE_FD
	echo "FPICFLAGS=$FPICFLAGS" >& AS_MESSAGE_FD
	echo "SAFE_FFLAGS=$SAFE_FFLAGS" >& AS_MESSAGE_FD
	echo "FLIBS=$FLIBS" >& AS_MESSAGE_FD
	echo "BLAS_LIBS=$BLAS_LIBS" >& AS_MESSAGE_FD
	echo "LAPACK_LIBS=$LAPACK_LIBS" >& AS_MESSAGE_FD
	echo "" >& AS_MESSAGE_FD
fi

dnl Process Makevars.in to make Makevars
AC_OUTPUT(src/Makevars)
