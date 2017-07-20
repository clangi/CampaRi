AC_INIT([CampaRi],[version 0.8.8],[dgarolini@gmail.com])

dnl Small macro to print out a nice, pretty section title.
define(SECTION_TITLE,
[
	dnl echo >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
	echo '   $1' >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
])


dnl Writing the version to screen
SECTION_TITLE([CampaRi package, AC_PACKAGE_VERSION])
#AC_LANG(Fortran)
AC_PROG_FC
AC_PATH_PROG(WHICH, which)
AC_PATH_PROG(ECHO, echo)
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(BASH, bash)
AC_PATH_PROG(RM, rm)
AC_PATH_PROG(MV, mv)
AC_PATH_PROG(MAKE, make)
AC_PATH_PROG(CD, cd)
if test "$CD" = "no" -o "$CD" = ""; then
  CD="cd"
fi
AC_PATH_PROG(TAR, tar)
AC_PATH_PROG(MKDIR, mkdir)
AC_PATH_PROG(PWD, pwd)
AC_PATH_PROG(DIRNAME, dirname)
AC_PATH_PROG(BASENAME, basename)

m4_include([tools/m4/pkg.m4])
PKG_PROG_PKG_CONFIG([0.20])
# AC_PATH_PROG(PKG_CONFIG, pkg-config)

dnl list of all possible optional components
all_options="netcdf4 campari"

dnl When a component is found, append it here
options=

dnl Force the compiler to run once and do all the basic checks
dnl if you don't, it will do the test on the first invocation
dnl below and so your pretty section titles won't work as well
dnl AC_CHECK_HEADERS(iostream)

SECTION_TITLE([Checking for netcdf4 libraries])

dnl include m4 macro for netcdf4 check
m4_include([tools/m4/ax_lib_netcdf4.m4])

AX_LIB_NETCDF4()
if test "$with_netcdf4_fortran" = "yes"; then
	dnl This add to the NETCDF4_LIBS was needed to compile

  # test whether netcdf target compiler is the same of the selected one
  AC_MSG_CHECKING([whether netcdf target compiler is the same compiler that has been selected])
  if test `$BASENAME "${NETCDF4_FC}"` = `$BASENAME "${FC}"`; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
		AC_MSG_CHECKING([if netcdf has been built with cmake (not supported for nf-config yet)])
		testing_cmake_exception=`$ECHO "${NETCDF4_FC}" | $GREP cmake`
		if test -n testing_cmake_exception; then
			AC_MSG_RESULT([yes])
			AC_MSG_NOTICE([NetCDF has been installed using cmake which still lacks a proper nf-config linking])
			AC_MSG_NOTICE([The message on the NETCDF4 FC was ${NETCDF4_FC}])
      AC_MSG_WARN([If you compiled NETCDF4 with a different compiler than $FC the compilation will crash])
			AC_MSG_CHECKING([again if netcdf was built with Fortran bindings])
			if test -z "$NC_CONFIG"; then
				AC_MSG_RESULT([no])
				AC_MSG_WARN([We were not able to find nc-config])
				with_netcdf4_fortran=no
			else
				did_it_has_fortran=`$NC_CONFIG --has-fortran`
				if test -n "$did_it_has_fortran"; then
					AC_MSG_RESULT([yes])
					AC_MSG_NOTICE([Trying to bind netcdf with usual libraries (guessing -lnetcdff)])
					NETCDF4_FC="$FC"
					NETCDF4_FFLAGS="${NETCDF4_CFLAGS}"
					NETCDF4_FLIBS="-lnetcdff ${NETCDF4_LIBS}"
				else
					AC_MSG_RESULT([no])
					AC_MSG_WARN([No Fortran bindings found. Option disabled])
					with_netcdf4_fortran=no
				fi
			fi
		else
		 	AC_MSG_RESULT([no])
			AC_MSG_WARN([

  ATTENTION!! The target compiler used for netcdf4 is different from the
  selected one for campari. Please consider reinstalling netcdf4 using
  the correct --exec-prefix.
])
			with_netcdf4_fortran=no
		fi
  fi
else
	AC_MSG_WARN([

	Unable to find NetCDF4 and set Fortran bindings. Please be
	sure that command nc-config(nf-config for recent versions) exists
	and it is not competing with previous versions (nc-config should not
	have fortran bindings if nf-config is provided).
])
	with_netcdf4_fortran=no
fi

if test "$with_netcdf4_fortran" = "no"; then
  AC_MSG_WARN([

	ATTENTION!! It is not advisable to avoid the use netcdf memory handling.
])
AC_MSG_WARN([NETCDF Fortran bindings not find. Memory handling will be limited to R.])
	NETCDF4_CC=""
	NETCDF4_VERSION=""
	NETCDF4_CFLAGS=""
	NETCDF4_CPPFLAGS=""
	NETCDF4_LDFLAGS=""
	NETCDF4_LIBS=""
	NETCDF4_FC="$FC"
	NETCDF4_FFLAGS=""
	NETCDF4_FLIBS=""
	SOURCES="CampaRi_init.c m_variables_gen.f90 m_gutenberg.f90 m_clustering.f90 m_gen_nbls.f90 m_mst.f90 fprintf_wrapper.c m_hw_fprintf.f90 main_clu_adjl_mst.f90 utilities.f90 dist_clusters_utils.f90 gen_progind.f90 gen_manycuts.f90 contract_mst.f90"
	PREPROCFLAGS="-cpp"
else
	SOURCES="CampaRi_init.c m_variables_gen.f90 m_gutenberg.f90 m_clustering.f90 m_gen_nbls.f90 m_mst_netcdf.f90 m_mst.f90 fprintf_wrapper.c m_hw_fprintf.f90 main_clu_adjl_mst.f90 utilities_netcdf.f90 utilities.f90 dist_clusters_utils.f90 gen_progind.f90 gen_manycuts.f90 contract_mst.f90"
	options="$options netcdf4"
	PREPROCFLAGS="-cpp -DLINK_NETCDF"
fi

dnl creating the objects needed
OBJECTS=$($ECHO "$SOURCES" | $SED 's/.f90/.o/g' | $SED 's/\.c/.o/g' )

echo "$OBJECTS"
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

# R_PACKAGE_DIR=$PWD # TO DELETE
R_PACKAGE_DIR="." # TO DELETE
CAMPARI_HOME="${R_PACKAGE_DIR}/inst/campari_lib"
AC_MSG_CHECKING([if original version of CAMPARI must be installed])
AC_ARG_WITH(campari,
	[AS_HELP_STRING([--with-campari=yes/<dir>],
     [install original campari in <dir>. If yes <dir> is defaulted to ${R_PACKAGE_DIR}/inst/campari_lib])],,
  [with_campari=no])
AC_MSG_RESULT([${with_campari}])
if test "${with_campari}" != "no"; then
	if test "${with_campari}" = "yes"; then
		if test ! -d "${CAMPARI_HOME}" -o -z "${CAMPARI_HOME}"; then
			AC_MSG_ERROR([Impossible to set automatically CAMPARI_HOME. please consider inserting it manually (--with-campari).])
		fi
	else
	  case $with_campari in
	    */)
	    if test -d "$with_campari"; then
	      CAMPARI_HOME="$with_campari"
	    else
	      AC_MSG_ERROR([campari directory specified not present.])
	    fi
	    ;;
	    *) AC_MSG_ERROR([bad value $with_campari for --with_campari (remember to use / as last char).]);;
	  esac
	fi
fi

# ---------------------------
#    CAMPARI to install
# ---------------------------
if test "${with_campari}" != "no"; then
	AC_MSG_ERROR([


ATTENTION! The --with-campari option is still under development and it is not directly available.

])
	campari_installed="0"
	SECTION_TITLE([Original campari library installing...])
	AC_MSG_CHECKING([which folder will be used for CAMPARI_HOME])
	AC_MSG_RESULT([$CAMPARI_HOME])
	# test if files are in the right place
	if test -f "${CAMPARI_HOME}/source/chainsaw.f90" -a -d "${CAMPARI_HOME}/bin/" -a -d "${CAMPARI_HOME}/lib/" -a -d "${CAMPARI_HOME}/params/"; then
		AC_MSG_NOTICE([
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Installing Original Campari
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
])
		CONFIG_TO_RUN="${CAMPARI_HOME}/source/configure --with-campari-home=${CAMPARI_HOME}"
		AC_MSG_NOTICE([running the configuration file of campari using:
${CONFIG_TO_RUN}])
		# configuration script run
		if AC_TRY_EVAL(CONFIG_TO_RUN); then
			AC_MSG_NOTICE([
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Configuration completed succesfully
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
])
			AC_MSG_NOTICE([using make to install campari.])
			AC_MSG_NOTICE([or multi-threaded and mpi installation, please check campari configuration file.])
		  to_go_back_dir=$PWD
			AC_MSG_CHECKING([whether $to_go_back_dir is source])
		  if test "$to_go_back_dir" != "${CAMPARI_HOME}/source"; then
		    AC_MSG_RESULT([no])
		    doit="$CD ${CAMPARI_HOME}/source"
		    AC_MSG_NOTICE([Changing directory using:
		$doit])
		    if AC_TRY_EVAL(doit); then
		      AC_MSG_NOTICE([directory sucessfully changed in ${CAMPARI_HOME}/source/.])
		      AC_MSG_NOTICE([previous directory will be reset after configuration.])
		    else
		      AC_MSG_ERROR([Something went wrong during directory changing.])
		    fi
		  else
		    AC_MSG_RESULT([yes])
		  fi
			doit="$MAKE campari"
			AC_MSG_NOTICE([running make in this shape:
$doit
])
		if AC_TRY_EVAL($doit); then
			campari_installed=1
			options="$options campari"
			echo "FREAKING DONE"
		else
			campari_installed=0
		fi
		if test "$to_go_back_dir" != "${CAMPARI_HOME}/source/"; then
			AC_MSG_NOTICE([setting back original directory])
		  doit="$CD ${to_go_back_dir}"
		  if AC_TRY_EVAL(doit); then
		    AC_MSG_NOTICE([directory sucessfully changed to previous value: ${to_go_back_dir}.])
		  fi
		fi
		# unsuccessful configuration
		else
			campari_installed=0
		fi
		if test "${campari_installed}" = "0"; then
			AC_MSG_WARN([

--------------------------------------------------------
WARNING!! Installation of original campari did not end
well. If you want full functionality, consider to
install it separately and link it easily using the R
functions within the package CampaRi.
--------------------------------------------------------

])
			AC_MSG_NOTICE([the installation of the package will continue without this feature.
])
		fi
	# impossible to find the structure and source directory.
	else
		AC_MSG_ERROR([
Campari source files and structure not found in $CAMPARI_HOME.
Please copy them from inst/campari_lib.
])
	fi
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






echo "" >& AS_MESSAGE_FD
SECTION_TITLE([Important variables])
# netcdf4
# ---------------------------------------------------------------------------
echo "PREPROCFLAGS=$PREPROCFLAGS" >& AS_MESSAGE_FD
AC_SUBST(PREPROCFLAGS)
AC_SUBST(HAVE_NETCDF4)
AC_SUBST(NETCDF4_VERSION)
AC_SUBST(SOURCES)
AC_SUBST(OBJECTS)
# AC_SUBST(NETCDF4_CC)
# AC_SUBST(NETCDF4_CFLAGS)
# AC_SUBST(NETCDF4_CPPFLAGS)
AC_SUBST(NETCDF4_LDFLAGS)
AC_SUBST(NETCDF4_LIBS)
AC_SUBST(NETCDF4_FC)
AC_SUBST(NETCDF4_FFLAGS)
AC_SUBST(NETCDF4_FLIBS)
AC_SUBST(with_netcdf4_fortran)
if test "$with_netcdf4_fortran" = "yes"; then
  # echo "HAVE_NETCDF4=1" >& AS_MESSAGE_FD
	echo "NETCDF4_VERSION=$NETCDF4_VERSION" >& AS_MESSAGE_FD
	echo "NETCDF4_FC=$NETCDF4_FC" >& AS_MESSAGE_FD
	echo "NETCDF4_FFLAGS=$NETCDF4_FFLAGS" >& AS_MESSAGE_FD
	# echo "NETCDF4_FCLAGS=$NETCDF4_FCLAGS" >& AS_MESSAGE_FD # it is not set
  echo "NETCDF4_LDFLAGS(not used)=$NETCDF4_LDFLAGS" >& AS_MESSAGE_FD
	echo "NETCDF4_FLIBS=$NETCDF4_FLIBS" >& AS_MESSAGE_FD
	echo '--' >& AS_MESSAGE_FD
else
	echo "NO NETCDF4 ACTIVE" >& AS_MESSAGE_FD
	echo '--' >& AS_MESSAGE_FD
fi

echo ' ' >& AS_MESSAGE_FD
echo "From R session:" >& AS_MESSAGE_FD
echo ' ' >& AS_MESSAGE_FD
echo "FC=$FC" >& AS_MESSAGE_FD
echo "FCFLAGS=$FCFLAGS" >& AS_MESSAGE_FD
echo "SHLIB_FFLAGS=$SHLIB_FFLAGS" >& AS_MESSAGE_FD
echo "FPICFLAGS=$FPICFLAGS" >& AS_MESSAGE_FD
echo "SAFE_FFLAGS=$SAFE_FFLAGS" >& AS_MESSAGE_FD
echo "FLIBS=$FLIBS" >& AS_MESSAGE_FD
echo "BLAS_LIBS=$BLAS_LIBS" >& AS_MESSAGE_FD
echo "LAPACK_LIBS=$LAPACK_LIBS" >& AS_MESSAGE_FD
echo "R_PACKAGE_DIR=$R_PACKAGE_DIR" >& AS_MESSAGE_FD
echo "R_LIBRARY_DIR=$R_LIBRARY_DIR" >& AS_MESSAGE_FD
echo "R_PACKAGE_NAME=$R_PACKAGE_NAME" >& AS_MESSAGE_FD
echo "" >& AS_MESSAGE_FD



dnl Process Makevars.in to make Makevars
AC_OUTPUT(src/Makevars)
