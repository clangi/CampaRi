# ------------------------------------------------
#
#   		 CAMPARI v3 AUTOCONF
#
#
# by Davide Garolini
# ------------------------------------------------
#
# This file is generating the configure file using
# autoconf (>2.59).
#
#
# Basic needs:
# - A Fortran03 compiler
#
#
# Needed libraries:
# - netcdf4 (>=4.4.0)
# - XDR-library (with Fortran bindings)
# - FFTW3
# - Lapack library

# Usable compilers:
# - intel
# - sun # command not specified correctly
# - gnu
# - cray # command to be specified
#


# Running autoconf on this file will trigger a warning if
# autoconf is not at least the specified version.
AC_PREREQ([2.69])
AC_LANG(Fortran)
# Initialize with name, version, and support email address.
AC_INIT([CAMPARI],[version 3.0],[davide.garolini@uzh.ch])
# Create the VERSION file, which contains the package version from
# AC_INIT.

echo AC_PACKAGE_VERSION>VERSION
AC_SUBST(PACKAGE_VERSION)

dnl Small macro to print out a nice, pretty section title.
define(SECTION_TITLE,
[
        dnl echo >& AS_MESSAGE_FD
        echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
        echo '   $1' >& AS_MESSAGE_FD
        echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
])

dnl Writing the version to screen
SECTION_TITLE([AC_PACKAGE_NAME, AC_PACKAGE_VERSION])

# Keep libtool macros in an m4 directory.
# AC_CONFIG_MACRO_DIR([m4])
# This m4 seems to do actually nothing. It is not needed.

# Find out about the target we're building for.
#AC_CANONICAL_TARGET

# Check for the existence of this file before proceeding.
# AC_CONFIG_SRCDIR([parsekey.f90])


# -----------------------------------------
#   General language and compiler init
# -----------------------------------------

# NB on AM_PROG_CC_C_O: Test whether the Fortran compiler accepts the options -c and -o simultaneously, and define F77_NO_MINUS_C_MINUS_O or FC_NO_MINUS_C_MINUS_O, respectively, if it does not.

# Find fortran compiler automatically and setting many other useful commands
AC_MSG_NOTICE([finding standard Fortran compiler])
# AC_LANG(Fortran)
AC_PROG_FC
AC_PROG_F77($FC)
AC_PROG_INSTALL # failing because no install-sh (AC_command used in m4s)
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP
AC_PATH_PROG(WHICH, which)
AC_PATH_PROG(BASH, bash)
AC_PATH_PROG(RM, rm)
AC_PATH_PROG(MV, mv)
AC_PATH_PROG(CP, cp)
AC_PATH_PROG(MAKE, make)
AC_PATH_PROG(ECHO, echo)
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(EGREP, egrep)
AC_PATH_PROG(TR, tr)
AC_PATH_PROG(TOUCH, touch)
# why this is impossible to find???
AC_PATH_PROG(CD, cd)
if test "$CD" = "no" -o "$CD" = ""; then
  AC_MSG_NOTICE([CD variable is auto-assigned to cd because ac path prog could not find it.])
  CD="cd"
fi
AC_PATH_PROG(TAR, tar)
AC_PATH_PROG(MKDIR, mkdir)
AC_PATH_PROG(PWD, pwd)
AC_PATH_PROG(DIRNAME, dirname)
AC_PATH_PROG(BASENAME, basename)

# usually it is automatically linked when using pkg-config > 0.19
m4_include([m4/pkg.m4])
PKG_PROG_PKG_CONFIG([0.20]) # it sets PKG_CONFIG* stuff
# AC_PATH_PROG(PKG_CONFIG, pkg-config)


# Set output variable INSTALL to the name of a BSD-compatible install program,
# if one is found in the current PATH. Otherwise, set INSTALL to ‘dir/install-sh -c’,
# checking the directories specified to AC_CONFIG_AUX_DIR
# (or its default directories) to determine dir (see Output).Also set the variables
# INSTALL_PROGRAM and INSTALL_SCRIPT to ‘${INSTALL}’ and INSTALL_DATA to ‘${INSTALL} -m 644’.


# --------------------------------
#     Setting up CAMPARI_HOME
# --------------------------------

# Default
working_dir=$PWD
if test -d "${working_dir}/source"; then
  CAMPARI_HOME="${working_dir}"
else
  CAMPARI_HOME="$($DIRNAME "${working_dir}")"
fi

# Inserted value?
AC_MSG_NOTICE([determining CAMPARI home (can be set --with-campari-home))])
AC_ARG_WITH(campari-home,
	[AS_HELP_STRING([--with-campari-home=<dir>],
     [use campari home <dir>. Defaults to pwd])],
  [])
if test -n "$with_campari_home"; then
  case $with_campari_home in
    */*)
    if test -d "$with_campari_home"; then
      CAMPARI_HOME="$with_campari_home"
    else
      AC_MSG_ERROR([campari-home specified not present.])
    fi
    ;;
    *) AC_MSG_ERROR([bad value $with_campari_home for --with_campari_home.]);;
  esac
fi

# check if the structure of CAMPARI_HOME is coherent
if test -f "${CAMPARI_HOME}/source/chainsaw.f90" -a \
-d "${CAMPARI_HOME}/params/"; then
  if test ! -d "${CAMPARI_HOME}/bin/" -a ! -d "${CAMPARI_HOME}/lib/"; then
    AC_MSG_WARN([bin or lib sub-directories have not been found in the tree.
They will be created if possible.])
  fi
  AC_MSG_NOTICE([CAMPARI_HOME is set to $CAMPARI_HOME])
  AC_MSG_CHECKING([whether $working_dir is source])
  to_go_back_dir=$working_dir # directory to go back to if not in source
  if test "$to_go_back_dir" != "${CAMPARI_HOME}/source"; then
    AC_MSG_RESULT([no])
    doit="$CD ${CAMPARI_HOME}/source/"
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
else
	AC_MSG_ERROR([
  Campari source files and structure not found in $CAMPARI_HOME.
  Please consider if you have a complete version of campari and
  then cd in CAMPARI_HOME/source and run again the configuration
  file.
	])
fi

AC_MSG_CHECKING([if the main directory has bin sub-directory])
if test -d "${CAMPARI_HOME}/bin/"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  doit="$MKDIR ${CAMPARI_HOME}/bin/"
  AC_MSG_NOTICE([creating directory using:
$doit])
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([directory sucessfully created.])
  else
    AC_MSG_ERROR([Something went wrong during directory creation.])
  fi
fi

AC_MSG_CHECKING([if the main directory has lib sub-directory])
if test -d "${CAMPARI_HOME}/lib/"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  doit="$MKDIR ${CAMPARI_HOME}/lib/"
  AC_MSG_NOTICE([creating directory using:
$doit])
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([directory sucessfully created.])
  else
    AC_MSG_ERROR([Something went wrong during directory creation.])
  fi
fi

# checking existance of DEPENDENCIES and Makefile.local
AC_MSG_CHECKING([existance of DEPENDENCIES])
if test -e "${CAMPARI_HOME}/source/DEPENDENCIES"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  doit="$TOUCH ${CAMPARI_HOME}/source/DEPENDENCIES"
  AC_MSG_NOTICE([creating empty file using:
$doit])
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([Empty file sucessfully created.])
  else
    AC_MSG_ERROR([Something went wrong during empty file creation.])
  fi
fi
AC_MSG_CHECKING([existance of Makefile.local])
if test -e "${CAMPARI_HOME}/source/Makefile.local"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  doit="$TOUCH ${CAMPARI_HOME}/source/Makefile.local"
  AC_MSG_NOTICE([creating empty file using:
$doit])
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([Empty file sucessfully created.])
  else
    AC_MSG_ERROR([Something went wrong during empty file creation.])
  fi
fi

# -------------------------------------------
#       Checking the cput architecture
# -------------------------------------------

# (from ax_blas checks)
AC_FC_DUMMY_MAIN
AC_F77_DUMMY_MAIN

# (from blas)
AC_CANONICAL_HOST # to check other variables
AC_MSG_CHECKING([which host cpu is present])
AC_MSG_RESULT([$host_cpu])
if test -n "${host_cpu}"; then
  ARCH="$host_cpu"
else
  AC_MSG_ERROR([Could not find cpu architecture using ac_canonical_host])
fi

build_linux=no
build_windows=no
build_mac=no
AC_MSG_CHECKING([which host os is running])
AC_MSG_RESULT([$host_os])
# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=yes
        system_os="linux"
        # make_dep_script="make_dependencies_linux.sh"
        make_dep_script="make_dependencies.sh" # it got updated to be cross-platform
        ;;
    cygwin*|mingw*)
        build_windows=yes
        system_os="windows"
        ;;
    darwin*)
        build_mac=yes
        system_os="macOS"
        # make_dep_script="make_dependencies_mac.sh"
        make_dep_script="make_dependencies.sh"
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

if test "${host_os}" = "windows"; then
  AC_MSG_ERROR([Windows is not a configurable system for CAMPARI yet.])
fi

# counting number of cpus available:
m4_include([m4/ax_count_cpus.m4])
AX_COUNT_CPUS()

AC_MSG_CHECKING([if the bin directory has ${ARCH} sub-directory])
if test -d "${CAMPARI_HOME}/bin/${ARCH}"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  doit="$MKDIR ${CAMPARI_HOME}/bin/${ARCH}"
  AC_MSG_NOTICE([creating directory using:
$doit])
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([directory sucessfully created.])
  else
    AC_MSG_ERROR([Something went wrong during directory creation.])
  fi
fi

AC_MSG_CHECKING([if the lib directory has ${ARCH} sub-directory])
if test -d "${CAMPARI_HOME}/lib/${ARCH}"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  doit="$MKDIR ${CAMPARI_HOME}/lib/${ARCH}"
  AC_MSG_NOTICE([creating directory using:
$doit])
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([directory sucessfully created.])
  else
    AC_MSG_ERROR([Something went wrong during directory creation.])
  fi
fi


# list of all possible optional components
# ----------------------------------------
all_needed_options="fortran03 LAPACK"
all_flavour_options="MPI specific_compiler warning_flags debug_flags threads NetCDF4 xdr_lib hsl_lib FFTW3"
# When a component is found, append it here
needed_options=
flavour_options=
# ----------------------------------------


# If F03 is set, FC isn't, then set FC to whatever F03 was set
if test "${F03+set}" = set -a "${FC+set}" != set; then
   FC=$F03
fi

# If F90 is set, FC isn't, then set FC to whatever F90 was set to.
if test "${F90+set}" = set -a "${FC+set}" != set; then
   FC=$F90
fi

# If we couldn't find an F90 compiler, deal with it.
if test "${FC+set}" = set && test "x$FC" = x; then
   AC_MSG_FAILURE([Can't find F90 compiler.])
fi



# --------------------------------
#    MPI option
# --------------------------------

AC_LANG_PUSH(Fortran)
m4_include([m4/ax_prog_fc_mpi.m4])
m4_include([m4/ax_mpi.m4])

# If --with-mpi=auto is used, try to find MPI, but use standard FC compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard FC compiler instead.

SECTION_TITLE([MPI option compiler and flags])
AC_MSG_CHECKING(whether the mpi layer mode is active)
AC_ARG_ENABLE(mpi,
AS_HELP_STRING([--enable-mpi=yes/auto/no/<exe>],
               [enable mpi compilation, default: no]),
[case "${enableval}" in
        "" | yes | auto)
        mpi="auto";;
        no) mpi="no" ;;
        *) mpi="${enableval}" ;; esac],
[mpi="no"])
if test "${mpi}" != "no"; then
  AC_MSG_RESULT([yes])
  if test "${mpi}" = "auto"; then
    AX_MPI([AC_MSG_NOTICE([MPI support found with MPIFC=$MPIFC and MPILIBS=$MPILIBS])],
      [AC_MSG_WARN([initial search did not go well. MPI option disabled.])]
      mpi=no)
    # MPIFF="$MPIFC" #this will be done at last (after the fc checks)
    mpi="$MPIFC"
  else
    AC_MSG_WARN([mpi compiler inserted manually. It will be checked automatically.])
  fi
  AC_MSG_NOTICE([selected mpi compiler is: $mpi])

  # building the directory tree
  AC_MSG_CHECKING([if the lib directory has mpi sub-directory])
  if test -d "${CAMPARI_HOME}/lib/${ARCH}/mpi"; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    doit="$MKDIR ${CAMPARI_HOME}/lib/${ARCH}/mpi"
    AC_MSG_NOTICE([creating directory using:
  $doit])
    if AC_TRY_EVAL(doit); then
      AC_MSG_NOTICE([directory sucessfully created.])
    else
      AC_MSG_ERROR([Something went wrong during directory creation.])
    fi
  fi
else
  AC_MSG_RESULT([no])
fi
# AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
#   use_mpi=no
#   if test x"$with_mpi" = xyes; then
#     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
#   else
#     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
#   fi
# ])
AC_LANG_POP(Fortran)

# --------------------------------
#    specific compilers and flags
# --------------------------------
SECTION_TITLE([Compiler choice and specific flags settings])

# General FCFLAGS from campari (I don't know how to
# look for them. They are cross compilers)

# preprocessor directives
PREPROCFLAGS="-DLINK_LAPACK -DDISABLE_FLOAT "
# standard is gfortran (please consider that cpp is fpp in intel and so on)

# standard flags for different compilers
INTELDEFAULTS="-real-size 64 -qoverride-limits "
SUNDEFAULTS="-r8const -ftrap=common,no%overflow"
GNUDEFAULTS="-fall-intrinsics -ffpe-trap=invalid,zero  -fdefault-real-8 -fdefault-double-8"
CRAYDEFAULTS="-K trap=inv,divz -emfF -dp -sreal64 -J '.'"
PGIDEFAULTS="-Mallocatable=03 -Mfree -Mextend -Mbackslash -DPGI_FORTRAN"
ABSOFTDEFAULTS="-N113 -ffree -en -TENV:simd_imask -TENV:simd_zmask"
# autovectorization
# AUTOVEC="-msse4.2"
std_compiler="$FC"
avail_compilers=('absoft' 'gnu' 'cray' 'sun' 'intel' 'pgi')
STD_FCFLAGS='$FCFLAGS'

# enable specific compilation
# -------------------------------
AC_ARG_ENABLE(compiler,
AS_HELP_STRING([--enable-compiler=no/gnu/intel/(cray)/sun/pgi/<exe>],
               [enable specific compilation, default: no]),
[case "${enableval}" in
        no) compiler="no" ;;
        absoft|f95)  compiler="absoft" ;;
        gnu|gfortran) compiler="gnu" ;;
        cray) compiler="cray" ;;
        sun) compiler="sun" ;;
        intel|ifort) compiler="intel" ;;
        pgi|pgfortran) compiler="pgi" ;;
        *) compiler="${enableval}" ;; esac],
[compiler="no"])

AC_MSG_CHECKING([whether a specific compiler have been inserted])
AC_MSG_RESULT([$compiler])

# main flags settings
if test "$compiler" != "no"; then
AC_MSG_CHECKING([whether $FCFLAGS])
  if test "${mpi}" = "no"; then
    AC_MSG_NOTICE([ specific compilers will have default flags set.
To override specific flags please use the Makefile.local
file in the source directory or define the variable UFCFLAGS
(possible also using --with-user-fcflags)(it overrides FCFLAGS only).])

    flavour_options="${flavour_options} specific_compiler"
    # Main compiler switch
    case "${compiler}" in
      gnu) FC="gfortran"
      FCFLAGS="$GNUDEFAULTS"
      AUTOVEC="-march=native"
      OPENLOOPS="-funroll-loops"
      PREPROCFLAGS_specific="-cpp"
      OPTIMIZATION_LEVEL="-O3";;

      absoft) FC="f95"
      FCFLAGS="$ABSOFTDEFAULTS"
      PREPROCFLAGS_specific="-cpp"
      OPTIMIZATION_LEVEL="-O3"
      AUTOVEC="-march=host";;

      # todo
      cray) FC="cray"
      FCFLAGS="$CRAYDEFAULTS";;
      sun) FC="sun"
      FCFLAGS="$FCFLAGS $SUNDEFAULTS"
      OPTIMIZATION_LEVEL="-xO4";;
      # todo

      intel) FC="ifort"
      FCFLAGS="$INTELDEFAULTS"
      PREPROCFLAGS_specific="-fpp"
      OPTIMIZATION_LEVEL="-O3";;

      pgi) FC="pgfortran"
      FCFLAGS="$PGIDEFAULTS"
      PREPROCFLAGS_specific="-Mpreprocess"
      OPTIMIZATION_LEVEL="-O3";;

      *) FC="${compiler}"
      SP_COMPILER="$FC"
      FCFLAGS="$FCFLAGS";;
    esac
  else
    AC_MSG_ERROR([

      If you wish to use a specific compiler for the MPI compilation,
      please use the --enable-mpi=<exe> option and not the --enable-compiler
      option.
      ])
  fi
fi

# MPI
if test "${mpi}" != "no"; then
  AC_MSG_WARN([unwrapping the mpi layer to see original compiler.])
  AC_MSG_CHECKING([if -show can unveil it])
  wrapped_compiler=$(eval $mpi -show | $AWK '{print $1;}')
  save_FC="$FC"
  if test -n "${wrapped_compiler}"; then
    AC_MSG_RESULT([yes])
    AC_MSG_NOTICE([resulting compiler is $wrapped_compiler])
    compiler=$wrapped_compiler
    FC=$wrapped_compiler
    wrapper_flags=$(eval $mpi -show | $AWK '{$1= ""; print $0}')
    if test -n "${wrapper_flags}"; then
      mpi_wrapper_flags=$(eval $ECHO "${wrapper_flags}" | grep "mpi")
      if test -n "${mpi_wrapper_flags}"; then
        AC_MSG_NOTICE([found following wrapper flags for $wrapped_compiler mpi layer:
$wrapper_flags])
      else
        AC_MSG_WARN([mpi layer wrapper flags not correct (no mpi directive available).
        This is an aberration of the -show option or the mpi compiler is not properly set.
        The mpi mode will disabled. Consider manual insertion.])
        mpi="no"
        FC="${save_FC}"
      fi
    else
      AC_MSG_WARN([mpi layer wrapper flags not found. This is an aberration of
the -show option. The mpi mode will disabled. Consider manual insertion.])
      mpi="no"
      FC="${save_FC}"
    fi
  else
    AC_MSG_RESULT([no])
    AC_MSG_NOTICE([continuing with inserted/found mpi compiler $mpi])
    compiler=$mpi
    FC=$mpi
  fi
  AC_MSG_NOTICE([future compiler assigned to $mpi])
fi


# finding the compiler location
if test "${compiler}" != "no"; then
  # checking directly the existance of compiler name
	if test -z "${SP_COMPILER}"; then
		AC_PATH_PROG(SP_COMPILER, $FC)
    if test -z "${SP_COMPILER}"; then
      AC_MSG_WARN([
  Configuration PATH_PROG was not able to find $FC path. This is
  a possible aberration due to the full path insertion (e.g. with --enable-mpi).
  The conf will continue its checks also for this special checks.])
      SP_COMPILER=$FC
    fi
    compiler_basename=$(eval $BASENAME "${SP_COMPILER}")
  # case in which command of the compiler is inserted manually
  else
    AC_MSG_NOTICE([checking the existance of the inserted compiler])
    FC="${SP_COMPILER}"
    # here it would be possible also to use AC_PATH_PROG or which
    if test $(eval $DIRNAME "$FC") = "."; then
      AC_MSG_NOTICE([inserted command-compiler only. It will be searched in the system.])
      AC_PATH_PROG(SP_COMPILER, $FC)
      compiler_basename="$FC"
    else
      AC_MSG_NOTICE([checking if the inserted compiler is an existing (and x) file.])
      if test ! -x "${SP_COMPILER}"; then
        AC_MSG_ERROR([

  The specific compiler inserted ($SP_COMPILER) should be an
  existing specific exec file which was not found or was not x
  ])
      else
        AC_MSG_NOTICE([$SP_COMPILER exists and is x.])
        compiler_basename=$(eval $BASENAME "${SP_COMPILER}")
      fi
    fi
  fi
  AC_MSG_CHECKING([whether selected compiler is present in the system (ac_path_prog based)])
  if test -z "${SP_COMPILER}" ; then
    AC_MSG_RESULT([no])
    if test "${mpi}" = "no"; then
      AC_MSG_ERROR([Unable to locate the compiler ${FC} using ac_path_prog nor directly])
    else
      AC_MSG_WARN([
ATTENTION!! mpi compiler not found. The option will be turned off.])
      FC="$save_FC"
      mpi="no"
      mpifc_not_found="yes"
    fi
  else
    AC_MSG_RESULT([yes])
  fi

  # AC_PROG_FC([$FC]) # why is this asking itself if there is a GNU Fortran compiler?
  AC_MSG_CHECKING([which supported compiler was found])
  AC_MSG_RESULT([$FC])
  AC_MSG_CHECKING([which are its standard compiler flags])
  AC_MSG_RESULT([$FCFLAGS])
else
  compiler_basename="$FC"
fi

# specific redirection to the standard flags
if test -n "${compiler_basename}"; then
  AC_MSG_CHECKING([if inserted compiler corrisponds to default possibilities])
  if test "${compiler}" = "no"; then to_add_option=yes; fi
  case "${compiler_basename}" in
    ifort)
      FCFLAGS="$INTELDEFAULTS"
      PREPROCFLAGS_specific="-fpp"
      OPTIMIZATION_LEVEL="-O3"
      compiler="intel"
      AC_MSG_RESULT([$compiler]);;
    gfortran)
      FCFLAGS="$GNUDEFAULTS"
      AUTOVEC="-march=native"
      OPENLOOPS="-funroll-loops"
      PREPROCFLAGS_specific="-cpp"
      OPTIMIZATION_LEVEL="-O3"
      compiler="gnu"
      AC_MSG_RESULT([$compiler]);;
    pgf90|pgf95|pgfortran)
      FCFLAGS="$PGIDEFAULTS"
      compiler="pgi"
      PREPROCFLAGS_specific="-Mpreprocess"
      AC_MSG_RESULT([$compiler]);;
  esac
  if test "${compiler}" != "no" -a -n "${to_add_option}"; then
    flavour_options="${flavour_options} specific_compiler"
  fi
  if test "${compiler}" != "no"; then
    AC_MSG_NOTICE([Having selected a compiler (even if automatically
using the standard of the system) which was commonly
used and optimized by campari, it has been defined and
the specific fcflags have been used.])
    AC_MSG_CHECKING([which new compiler is present])
    AC_MSG_RESULT([$FC])
    AC_MSG_CHECKING([which new compiler flags are present])
    AC_MSG_RESULT([$FCFLAGS])
  fi
fi

if test -z "${PREPROCFLAGS_specific}"; then
  PREPROCFLAGS_specific="-cpp"
  AC_MSG_WARN([not having found any specific compiler already used,
the preprocessing directive is set automatically to
-cpp as gfortran standard.])
fi
PREPROCFLAGS="${PREPROCFLAGS_specific} ${PREPROCFLAGS}"

# adding the flags!
if test "${mpi}" != "no" -a -n "$wrapper_flags" -a -z "$mpifc_not_found"; then
  AC_MSG_NOTICE([The standard compiler have been overloaded with mpi layer of wrapper_flags:
$wrapper_flags])
  AC_MSG_CHECKING([if it is possible to set the mpi fcflags for $FC])
  wrapper_fcflags=$(eval $ECHO "${wrapper_flags}" | sed 's/-l.*//' | sed 's/-L.*//')
  if test -z "${wrapper_fcflags}"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([Then using it all as first layer (first call to flags and libs)])
    wrapper_fcflags="${wrapper_flags}"
    wrapper_libs=""
  else
    AC_MSG_RESULT([yes])
    AC_MSG_CHECKING([if it is possible to set the mpi libraries for $FC])
    wrapper_libs=$(eval $ECHO "${wrapper_flags}" | sed "s|${wrapper_fcflags}||")
    if test -z "${wrapper_libs}"; then
      AC_MSG_RESULT([no])
      AC_MSG_WARN([The mpi libraries will be held next to the compiler as it was not possible to split them from the other flags.])
      wrapper_libs=""
      wrapper_fcflags="${wrapper_flags}"
    else
      AC_MSG_RESULT([yes])
    fi
  fi
fi
echo $compiler
# warning on the loss of standard FCFLAGS
AC_MSG_NOTICE([if a specific compiler have been inserted, the standard FCFLAGS of the system will be lost])
case $compiler in
  absoft|gnu|cray|sun|intel|pgi)
    AC_MSG_WARN([the following standard FCFLAGS will be lost because you use a specific compiler (we defined best flags)])
    AC_MSG_NOTICE([lost flags: $STD_FCFLAGS]);;
esac

# enable fast computation!! i.e. strip of the extra flags and optimizaitons
AC_ARG_ENABLE(fast-compilation,
AS_HELP_STRING([--enable-fast-compilation=yes/no],
               [enable fast compilation, default: no]),
[case "${enableval}" in
        "" | yes)
        fast_comp="yes";;
        no) fast_comp="no" ;;
        *) fast_comp="no" ;; esac],
[fast_comp="no"])
if test "${fast_comp}" = "yes"; then
  FCFLAGS=''
  OPTIMIZATION_LEVEL=''
  AUTOVEC=''
  OPENLOOPS=''
fi

# mpi layer on top of the existing flags
if test "${mpi}" != "no" -a -z "$mpifc_not_found"; then
  FLIBS=
  FCLIBS=
  if test -n "${mpi_wrapper_flags}"; then
    AC_MSG_NOTICE([$FC original compiler (with mpi flags) will be used instead of $mpi.])
    MPIFF="$FC"
    FCFLAGS="${wrapper_fcflags} ${FCFLAGS}"
    MPILIBS="${MPILIBS} ${wrapper_libs}"
  else
    AC_MSG_NOTICE([As the wrapper flags does not contain mpi directives, we will try to use the inserted mpi compiler: $mpi.])
    MPIFF="$mpi"
  fi
  flavour_options="${flavour_options} MPI"
fi


# -------------------------------
# overriding flag definition
AC_ARG_WITH(user-fcflags,
	[AS_HELP_STRING([--with-user-fcflags=<whatever>],
     [override the fcflags variable completely. Consider using Makefile.local])],
  [UFCFLAGS=""])
if test -n "$with_user_fcflags"; then
UFCFLAGS="${with_user_fcflags}"
fi
AC_ARG_VAR([UFCFLAGS], [user defined compiler flags (it overrides FCFLAGS)])
if test "$UFCFLAGS" != ""; then
  FCFLAGS=$UFCFLAGS
fi

# -------------------------------
# enable warn flags
AC_ARG_ENABLE(warnflags,
AS_HELP_STRING([--enable-warnflags=yes/no],
               [enable warnflags, default: no. Use UFCFLAGS to insert specific flags.]),
[case "${enableval}" in
        no)  warnflags="no" ;;
        yes) warnflags="yes" ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-warnflags]) ;; esac],
[warnflags="no"])
if test "$warnflags" = "yes"; then
  WARNFLAGS="-Wintrinsics-std -Wsurprising -Warray-bounds -Waliasing -Wampersand\
   -Wunderflow -Walign-commons -Wintrinsic-shadow -Wtabs -Wline-truncation -Wcharacter-truncation"
  # -Warray-temporaries  -Wimplicit-interface  -Wconversion
  flavour_options="${flavour_options} warning_flags"
fi

# -------------------------------
# debugging flags
AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug=yes/no],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=yes ;;
             no)  debug=no ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=no])
AC_MSG_CHECKING([whether compiler flags for debugging must be used])
AC_MSG_RESULT([$debug])
# standard search of debug variables for compilationif test $debug != false; then
if test "$debug" = "yes"; then
  flavour_options="${flavour_options} debug_flags"
  AC_MSG_CHECKING([which debugging flags will be set])
  INTELDEBUG="-check bounds -check uninit -check arg_temp_created -check pointers -traceback" #-g -debug extended -debug-parameters all
  SUNDEBUG="-xcheck=init_local -g"
  GNUDEBUG="-fbacktrace -fbounds-check -fcheck-array-temporaries -g"
  CRAYDEBUG="-emfFD"
  PGIDEBUG="-g"

  if test "$compiler" != "no"; then
    case "${compiler}" in
            gnu) DBFLAGS=$GNUDEBUG ;;
            cray) DBFLAGS=$CRAYDEBUG;;
            sun) DBFLAGS=$SUNDEBUG ;;
            intel) DBFLAGS=$INTELDEBUG;;
            pgi) DBFLAGS=$PGIDEBUG;;
    esac
  else
    DBFLAGS="-g"
  fi
  AC_DEFINE([HAVE_DBFLAGS], [1], [Defined if you activated debug flags.])
  AC_MSG_RESULT([$DBFLAGS])
fi

# -------------------------------
# multi-threading flags
#
# re='^[0-9]+$'
AC_MSG_CHECKING([how many CPUs are available])
AC_MSG_RESULT([${CPU_COUNT}])
n_threads=${CPU_COUNT}
m4_include([m4/ax_openmp.m4])
AC_ARG_ENABLE(threads,
AS_HELP_STRING([--enable-threads=yes/no/-thread_flags/n_threads],
               [enable multi-threading, default: no]),
[case "${enableval}" in
             '' | yes) threads=yes
             THREADFLAGS=
             ;;
             no) threads=no ;;
             -*) threads=yes
                 already_inserted=yes
                 a=1
                 THREADFLAGS="${enableval}";;
             *)
              if test $enableval -gt 0 ; then
                AC_MSG_NOTICE([selected $enableval threads])
                n_threads=$enableval
                THREADFLAGS=
                threads=yes
              else
                AC_MSG_ERROR([bad value ${enableval} for --enable-threads])
              fi
             ;;
esac],
[threads=no])
AC_MSG_CHECKING([whether compiler flags for multi-threading must be used])
AC_MSG_RESULT([$threads])

# standard search of debug variables for compilationif test $debug != false; then
if test "$threads" = "yes"; then

  # checking sub-directory tree
  AC_MSG_CHECKING([if the lib directory has ${ARCH}/threads sub-directory])
  if test -d "${CAMPARI_HOME}/lib/${ARCH}/threads"; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    doit="$MKDIR ${CAMPARI_HOME}/lib/${ARCH}/threads"
    AC_MSG_NOTICE([creating directory using:
  $doit])
    if AC_TRY_EVAL(doit); then
      AC_MSG_NOTICE([directory sucessfully created.])
    else
      AC_MSG_ERROR([Something went wrong during directory creation.])
    fi
  fi

  if test "$mpi" != "no"; then
    # checking sub-directory tree
    AC_MSG_CHECKING([if the lib directory has ${ARCH}/threads sub-directory])
    if test -d "${CAMPARI_HOME}/lib/${ARCH}/mpi_threads"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      doit="$MKDIR ${CAMPARI_HOME}/lib/${ARCH}/mpi_threads"
      AC_MSG_NOTICE([creating directory using:
    $doit])
      if AC_TRY_EVAL(doit); then
        AC_MSG_NOTICE([directory sucessfully created.])
      else
        AC_MSG_ERROR([Something went wrong during directory creation.])
      fi
    fi
  fi


  if test "$compiler" != "no"; then
    if test "${already_inserted}" != "yes"; then
      case "${compiler}" in
                gnu) THREADFLAGS="-fopenmp"
                      found_here="yes";;
                f95) THREADFLAGS="-openmp"
                      found_here="yes";;
                # cray) DBFLAGS=$CRAYDEBUG;;
                # sun) DBFLAGS=$SUNDEBUG ;;
                intel) THREADFLAGS="-openmp"
                        found_here="yes";;
                # pgi) DBFLAGS=$PGIDEBUG;;
      esac
      if test "${found_here}" = "yes"; then
        AC_MSG_NOTICE([automatic threading flags set for $compiler compiler: $THREADFLAGS])
        already_inserted="yes"
        a=1
      else
        AC_MSG_NOTICE([compiler $compiler not recognized in the supported pool for automatic settings.])
      fi
    else
      AC_MSG_NOTICE([even if an inserted compiler is present, the compiler flag has already been inserted manually.])
    fi
  fi

  # checking the number of threads
  AC_MSG_CHECKING([how many threads will be used])
  AC_MSG_RESULT([$n_threads])
  flavour_options="${flavour_options} threads"
  AC_MSG_NOTICE([checking the OpenMP support with non-canonical ax openmp])
  AX_OPENMP(
    [AC_MSG_NOTICE([OpenMP support found.])],
    [AC_MSG_WARN([no OpenMP support has been found directly. If specific compiler has been used the flag will be set automatically])]
  )
  AC_MSG_NOTICE([(overloading variables) checking the OpenMP support with canonical ac openmp])
  AC_OPENMP()
  a=2
  if test "${already_inserted}" != "yes"; then
    AC_MSG_NOTICE([analysing the possible compiling languages for NON-posix support (only OpenMP)])
    if test -n "${OPENMP_FFLAGS}"; then # note that it is wrong, FCFLAGS should be used instead of FFLAGS(F77)
      AC_MSG_NOTICE([found correct Fortran OpenMP flag (OPENMP_FFLAGS).])
      THREADFLAGS="$OPENMP_FFLAGS"
      a=1
    else
      AC_MSG_WARN([OpenMP support not found for Fortran compiler (OPENMP_FFLAGS)])
      AC_MSG_NOTICE([checking other possible flags (OPENMP_CFLAGS)])
      if test -n "${OPENMP_CFLAGS}"; then
        AC_MSG_NOTICE([found OpenMP flag for C compiler (OPENMP_CFLAGS).])
        THREADFLAGS="${OPENMP_CFLAGS}"
      else
        AC_MSG_WARN([OpenMP support not found for C compiler (OPENMP_CFLAGS)])
        AC_MSG_NOTICE([checking other possible flags (OPENMP_CXXFLAGS)])
        if test -n "${OPENMP_CXXFLAGS}"; then
          AC_MSG_NOTICE([found OpenMP flag for CXX compiler (OPENMP_CXXFLAGS).])
          THREADFLAGS="${OPENMP_CXXFLAGS}"
        else
          AC_MSG_WARN([OpenMP support not found for CXX compiler (OPENMP_CXXFLAGS)])
          AC_MSG_WARN([POSIX standard (GNU/Linux) -lpthread not implemented. If specific compiler has been inserted it will be set automatically.])
          # THREADFLAGS="-lpthread"
          a=0
        fi
      fi
    fi
  else
    AC_MSG_NOTICE([fcflags for multi-threading inserted manually. It has been set to $THREADFLAGS])
  fi

  if test "$a" = "0"; then
    if test "$compiler" = "no"; then
      AC_MSG_ERROR([
        No openmp support found. [pay attention to not set --disable-openmp
        and --enable-multi-thread at the same time]
        Canonical UNIX links (posix: -lpthread) could be added but present
        implementation needs openmp support.
        ])

    fi
  elif test $a = 1; then
    # OpenMP compiled libraries (fftw3 libs)
    AC_MSG_NOTICE([found openmp support for Fortran compiler. Threaded fftw3 flags set to openmp option.])
  else
    AC_MSG_NOTICE([found openmp support but not for the Fortran compiler. Threaded fftw3 flags set to openmp option.])
  fi
  AC_MSG_NOTICE([threaded fftw3 will check using these links: $FFTW3_LIBS_TMP])
  AC_MSG_CHECKING([which multi-threading flags will be set ])
  AC_MSG_RESULT([$THREADFLAGS])

  if test "${n_threads}" -gt "${CPU_COUNT}"; then
    AC_MSG_WARN([Inserted more threads than physical cores. Trying to compile with inserted values])
    AC_MSG_NOTICE([trying compilation using:
${FC} ${FCFLAGS} ${THREADFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90 -o conftest])
    AC_MSG_CHECKING([whether multi-threading is compiling using ${n_threads} threads])
    cat <<EOF >conftest.f90
program blah
  implicit none
  integer tpi, OMP_GET_THREAD_NUM
  call OMP_SET_NUM_THREADS(${n_threads})
EOF
    cat<<'EOF' >>conftest.f90
!$OMP parallel private(tpi)
  tpi = OMP_GET_THREAD_NUM()
!$OMP CRITICAL(PRINT)
  write(*,*) tpi
!$OMP END CRITICAL(PRINT)
!$OMP END parallel
end
EOF
    doit='${FC} ${FCFLAGS} ${THREADFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90 -o conftest'
    if AC_TRY_EVAL(doit); then
      AC_MSG_RESULT([yes])
      # try to evaluate it
      doit='./conftest'
      `$doit` >& threads_assignments_conftest
      # output=`$doit`
      # AC_MSG_NOTICE([$output])
      AC_MSG_NOTICE([running script to detect non-functional behavior])
      if AC_TRY_EVAL(doit); then
        AC_MSG_NOTICE([evaluation successful for ${n_threads} threads])
        AC_MSG_NOTICE([please check manually the file threads_assignments_conftest to see thread assignments])
      else
        AC_MSG_WARN([number of threads probably not supported. The system will try the standard number of cores($CPU_COUNT)])
        n_threads=${CPU_COUNT}
      fi
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([
        Unable to compile conftest with multi-threads (compiler = $FC, flags = $THREADFLAGS).
        Please check your compiler support.
        ])
    fi
  fi
  if test "${n_threads}" -le "${CPU_COUNT}"; then
    AC_MSG_NOTICE([trying compilation using:
${FC} ${FCFLAGS} ${THREADFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90-o conftest])
    AC_MSG_CHECKING([whether multi-threading is compiling using ${n_threads} threads])
    cat <<EOF >conftest.f90
program blah
  implicit none
  integer tpi, OMP_GET_THREAD_NUM
  call OMP_SET_NUM_THREADS(${n_threads})
EOF
    cat<<'EOF' >>conftest.f90
!$OMP parallel private(tpi)
  tpi = OMP_GET_THREAD_NUM()
!$OMP CRITICAL(PRINT)
  write(*,*) tpi
!$OMP END CRITICAL(PRINT)
!$OMP END parallel
end
EOF
    doit='${FC} ${FCFLAGS} ${THREADFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90 -o conftest'
    if AC_TRY_EVAL(doit); then
      AC_MSG_RESULT([yes])
      # try to evaluate it
      doit='./conftest'
      `$doit` >& threads_assignments_conftest
      # output=`$doit`
      # AC_MSG_NOTICE([$output])
      AC_MSG_NOTICE([running script to detect non-functional behavior])
      if AC_TRY_EVAL(doit); then
        AC_MSG_NOTICE([evaluation successful for ${n_threads} threads])
        AC_MSG_NOTICE([please check manually the file threads_assignments_conftest to see thread assignments])
      else
        AC_MSG_WARN([number of threads probably not supported. The system will try the standard number of cores($CPU_COUNT)])
        n_threads=${CPU_COUNT}
      fi
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([
        Unable to compile conftest with multi-threads (compiler = $FC, flags = $THREADFLAGS).
        Please check your compiler support.
        ])
    fi
  fi
fi

# -------------------------------
# testing general compilation
if test "$compiler" != "no"; then
  AC_MSG_NOTICE([trying compilation using:
${FC} -c ${FCFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90])
  AC_MSG_CHECKING([whether selected compiler and compiler flags can compile])
  cat <<EOF >conftest.f90
module conftest_module
!USE ISO_C_BINDING
end module conftest_module
EOF
  doit='${FC} -c ${FCFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Unable to compile with the selected compiler: ${compiler}.])
  fi
  # this is not known why or not we should need these
  AC_MSG_CHECKING([if compiler is the std compiler in order to set FLIBS and FCLIBS])
  if test "$FC" != "${std_compiler}" || \
  test "${compiler_basename}" != "${std_compiler}"; then
    AC_MSG_RESULT([no])
    FLIBS=
    FCLIBS=
  else
    AC_MSG_RESULT([yes])
  fi
fi

# ---------------------------------
#  f03 f08 TS29113 support check
# ---------------------------------
# Set the default fortran builds; default is to build f03
nc_build_f90=no
nc_build_f03=yes

# Does the user want to test for Fortran 2003 compiler.
# Added to get around NAG evaluation compiler error turning off F03
AC_MSG_CHECKING([whether Fortran compiler is checked for ISO_C_BINDING support])
AC_ARG_ENABLE([f0308-support],
              [AS_HELP_STRING([--disable-f0308-support],
                              [disable check of ISO_C_BINDING support in Fortran compiler])],[],[enable_f0308_support=yes])
test "x$enable_f0308_support" = xno || enable_f0308_support=yes
AC_MSG_RESULT($enable_f0308_support)

# disable-f03
AC_MSG_CHECKING([whether F03 native code is desired])
AC_ARG_ENABLE([f03],
              [AS_HELP_STRING([--disable-f03],
                              [suppress Campari Fortran 2003 native code])],[],[enable_f03=yes])
test "x$enable_f03" = xno || enable_f03=yes
AC_MSG_RESULT([$enable_f03])


# next check to see if the Fortran compiler will support
# ISO_C_BINDING
nc_iso_c_binding=yes

# Some f90 compilers change the case of the mod file names. Others
# require special flags to be set to deal with .f90 files. Find out
# about this compiler.
m4_include([m4/ax_f90_module_flag.m4])
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
AC_LANG_POP(Fortran)
AX_F90_MODULE_FLAG
AC_SUBST(MOD_FLAG, ${ax_cv_f90_modflag})

AC_MSG_CHECKING([if f03compiler was disabled])
if test "${enable_f03}" != "yes"; then
  enable_f0308_support=no
  AC_MSG_RESULT([yes])
  AC_MSG_NOTICE([therefore disabling also f0308-support.])
else
  AC_MSG_RESULT([no])
fi

# f03 ISO_C_BINDING check
# -----------------------
if test "x$enable_f0308_support" = xyes; then
  case "${compiler}" in
    intel) FCFLAGS="$FCFLAGS -stand f03";;
    gnu) FCFLAGS="$FCFLAGS -std=f2008";;
  esac
  AC_MSG_CHECKING([if Fortran compiler supports Fortran 2003 ISO_C_BINDING])
  cat <<EOF >conftest.f90
module conftest_module
USE ISO_C_BINDING
end module conftest_module
EOF
  doit='$FC -c ${FCFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
     nc_iso_c_binding=yes
   	needed_options="$needed_options fortran03"
  else
     nc_iso_c_binding=no
  fi
  AC_MSG_RESULT([$nc_iso_c_binding])


  # f08 ISO_FORTRAN_ENV check
  # ------------------------
  # Compiler supports Fortran 2008 ISO_FORTRAN_ENV - to ask if it is needed
  AC_MSG_CHECKING([if Fortran compiler supports Fortran 2008 ISO_FORTRAN_ENV additions])
  cat <<EOF >conftest.f90
Program test
USE ISO_FORTRAN_ENV, ONLY: REAL32, REAL64, INT8, INT16, INT32, INT64
End Program
EOF
  doit='$FC -c ${FCFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
     nc_f2008_is0_env_additions=yes
     AC_DEFINE([HAVE_F2008], 1, [if true, use Fortran 2008 ISO_FORTRAN_ENV additions])
  else
     nc_f2008_is0_env_additions=no
  fi
  AC_MSG_RESULT([$nc_f2008_is0_env_additions])

  # Compiler supports TS29113 - unknown what it is, bound to F03
  # ---------------------------------------------------------------
  AC_MSG_CHECKING([if Fortran compiler supports TS29113 standard extension])
  cat <<EOF >conftest.f90
Program test
USE ISO_C_BINDING, ONLY: C_PTRDIFF_T
End Program
EOF
  doit='$FC -c ${FCFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
     nc_ts29113_support=yes
     AC_DEFINE([HAVE_TS29113_SUPPORT], 1, [if true, Fortran compiler supports TS29113 extension])
  else
     nc_ts29113_support=no
  fi
  AC_MSG_RESULT([$nc_ts29113_support])

  if test "x$nc_iso_c_binding" = xno; then
    if test "x$enable_compiler_recover" = xyes; then
      nc_build_f90=yes
      AC_MSG_WARN([ISO_C_BINDING not supported - defaulting to F90 interfaces])
    else
      AC_MSG_FAILURE([Can't use build F03 interfaces and cannot do compiler recovery])
    fi
  fi
fi

# f90 wrappers or f03 native code
# --------------------------------
## Do we want to use the fortran 90 wrappers instead of 2003 native code?

AC_MSG_CHECKING([if f0308-support was disabled])
if test "${enable_f0308_support}" != "yes"; then
  enable_f03=no
  AC_MSG_RESULT([yes])
  AC_MSG_NOTICE([therefore disabling also f03.])
else
  AC_MSG_RESULT([no])
fi

if test "x$enable_f03" = xyes ; then nc_build_f03=yes; fi

# Decide on f90 versus 2003
if test "x$nc_iso_c_binding" = xno; then
  #Cannot build f03 native code
  AC_MSG_NOTICE([Fortran 2003 interface not supported; Using fortran 90 interface])
  nc_build_f90=yes
fi

# Guarantee build state; f90 currently overrides f03
if test "x$nc_build_f90" = xyes ; then
  nc_build_f03=no
else
  nc_build_f03=yes
fi

# Turn off building old F77 interfaces if nc_build_f03 is still yes
test "x$nc_build_f03" = xyes && nc_build_f77=no

if test "x$nc_build_f90" = xno; then
   F90=
   # If F77 is set, and FC isn't, then set FC to whatever F77 was set to.
   AS_IF([test x${F77+set} = xset -a x${FC+set} = x],[FC=$F77])

fi

# See if the fortran 90 build is desired. If so, set some stuff
if test "x$nc_build_f90" = xyes; then
   AC_MSG_NOTICE([setting up Fortran 90])
   if test "${F90+set}" != set; then
      F90=$FC
   fi
   AC_MSG_CHECKING([what FC is set to])
   AC_MSG_RESULT([$FC])

   if test "x$nc_build_f03" = xno; then
     # Some f90 compilers change the case of the mod file names. Others
     # require special flags to be set to deal with .f90 files. Find out
     # about this compiler.
     AC_LANG_PUSH(Fortran)
     AC_FC_SRCEXT(f90)
     AC_LANG_POP(Fortran)
     AX_F90_MODULE_FLAG
     AC_SUBST(MOD_FLAG, ${ax_cv_f90_modflag})
   fi

   # Set the FCLIBS flag to help with shared libraries.
  AC_FC_LIBRARY_LDFLAGS
  AC_F77_LIBRARY_LDFLAGS
fi

# --------------------------------
# Fortran types check (from netcf4)
# --------------------------------
# Does the user want to check into fortran type information? TODO make it better
AC_MSG_CHECKING([whether fortran type sizes should be checked])
AC_ARG_ENABLE([fortran-type-check],
              [AS_HELP_STRING([--disable-fortran-type-check],
                              [cause the Fortran type sizes checks to be skipped])])
test "x$enable_fortran_type_check" = xno || enable_fortran_type_check=yes
AC_MSG_RESULT($enable_fortran_type_check)

m4_include([m4/ud_misc.m4])
# This will test the f77 compiler, as well as check some fortran types.
if test "$enable_fortran_type_check" = yes; then
   UD_FORTRAN_TYPES
   if test "x$ud_fc_failed" = xyes; then
      AC_MSG_ERROR([F77 compiler doesn't work, and compiler recovery disabled.])
   fi

   # This else is mysterious to me
else
   AC_DEFINE(NCBYTE_T, byte)
   AC_DEFINE(NCSHORT_T, integer*2)
   AC_DEFINE(NF_INT1_T, byte)
   AC_DEFINE(NF_INT2_T, integer*2)
   AC_DEFINE(NF_INT1_IS_C_SIGNED_CHAR, 1, [default])
   AC_DEFINE(NF_INT2_IS_C_SHORT, 1, [default])
   AC_DEFINE(NF_INT_IS_C_INT, 1, [default])
   AC_DEFINE(NF_REAL_IS_C_FLOAT, 1, [default])
   AC_DEFINE(NF_DOUBLEPRECISION_IS_C_DOUBLE, 1, [default])
fi

# ac_cv_prog_f90_uppercase_mod=no
#
# # Test whether F90 module names get capitalized.
# if test "x$nc_build_f90" = xyes; then
#    AC_PROG_FC_UPPERCASE_MOD #undefined macro
# fi
# if test "x$ac_cv_prog_f90_uppercase_mod" = xyes ; then
# AC_MSG_FAILURE([Upper Case Modules names not supported])
# fi


# this conditional is failing because it is a directive to automake
# AM_CONDITIONAL(BUILD_F03, [test "x$nc_build_f03" = xyes]) # Is this needed?



# ----------------------------------
#    user options and library dir
# ----------------------------------
SECTION_TITLE([Checking dependencies and libraries folder])
AC_MSG_NOTICE([loading dependencies])
AC_MSG_CHECKING([existance of file DEPENDENCIES in source])
if test -e "${CAMPARI_HOME}/source/DEPENDENCIES"; then
  AC_MSG_RESULT([yes])
  AC_MSG_CHECKING([existance of file ${make_dep_script} in source])
  if test -e "${CAMPARI_HOME}/source/${make_dep_script}"; then
    AC_MSG_RESULT([yes])
    AC_MSG_WARN([To manually update dependencies please run ${make_dep_script} in the source directory.])
    AC_MSG_NOTICE([specific location of ${make_dep_script} is: ${CAMPARI_HOME}/source/${make_dep_script}])
    AC_MSG_NOTICE([Running ${make_dep_script} using default bash in current directory (source).])
    $($BASH ${CAMPARI_HOME}/source/${make_dep_script} ${CAMPARI_HOME}/source/)
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([${make_dep_script} file missing from source directory])
  fi
else
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([DEPENDENCIES file missing from source directory])
fi

# look for empty dir -- NO MORE IN USE - DEPRECATED
# -------------
EXT_LIBS_DIR="${CAMPARI_HOME}/ext_libs/"
if test "${XDR_LIBS}" = "" ; then
AC_MSG_CHECKING([whether any library is present in the $CAMPARI_HOME/ext_libs/ folder])
  if test -d "$EXT_LIBS_DIR" ; then
    if test -n "$(ls -A $EXT_LIBS_DIR)" ; then
      AC_MSG_RESULT([libs found])
  echo "
  $EXT_LIBS_DIR is not Empty. If the system requirement are not available,
  we will look for special libraries in this directory (xdr libs).
  If you have these libraries installed somewhere specifically,
  please refer their locations using --with-<specific-lib>=<lib-location>
  (Check help for specific options).
  "
      ax_ext_lib=yes
    else
      AC_MSG_RESULT([empty dir])
      ax_ext_lib=no
    fi
  else
    AC_MSG_RESULT([dir not found])
    ax_ext_lib=no
  fi
fi

# -------------------------------
#       Checking NetCDF4 libs
# -------------------------------

SECTION_TITLE([Checking for netcdf4 libraries])

dnl include m4 macro for netcdf4 check
m4_include([m4/ax_lib_netcdf4.m4])

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
				if test -n "${did_it_has_fortran}"; then
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

ATTENTION!! NO NETCDF4 FUNCTIONALITY ENABLED. It is not advisable to avoid the use netcdf.

])
  NETCDF4_CC=""
  NETCDF4_VERSION=""
  NETCDF4_CFLAGS=""
  NETCDF4_CPPFLAGS=""
  NETCDF4_LDFLAGS=""
  NETCDF4_LIBS=""
  NETCDF4_FC="$FC"
  NETCDF4_FFLAGS=""
  NETCDF4_FLIBS=""
else
  flavour_options="$flavour_options NetCDF4"
  PREPROCFLAGS="${PREPROCFLAGS} -DLINK_NETCDF"
fi
# -------------------------------
#       Checking LAPACK libs
# -------------------------------
SECTION_TITLE([Checking for LAPACK libraries])
# To link with LAPACK, you should link with:
# $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS

m4_include([m4/ax_blas.m4])
m4_include([m4/ax_lapack.m4])
AX_BLAS([AC_MSG_NOTICE([BLAS libraries found.])],
[AC_MSG_ERROR([Unable to find blas libraries. Please consider linking them using --with-blas=<lib>.])])
AX_LAPACK([AC_MSG_NOTICE([LAPACK libraries found.])],
[AC_MSG_ERROR([Unable to find lapack libraries. Please consider linking them using --with-lapack=<lib>.])])
needed_options="${needed_options} LAPACK"
# -------------------------------
#       Checking xdr libs
# -------------------------------
SECTION_TITLE([Checking for xdr libraries])

XDR_PREFIX=""
XDR_LIBS=""

# with-xdr ?
# -------------
# AC_ARG_WITH(xdr,
# 	[AS_HELP_STRING([--with-xdr=<lib>], [use XDR library <lib>])],
#   [])
# case $with_xdr in
#   yes | "") with_xdr=yes ;;
#   *.a | *.so | *.so.* | *.o)
#   if test -e $with_xdr; then
#     XDR_PREFIX=$(dirname $with_xdr)
#   else
#     AC_MSG_WARN([library specified not present. It will be searched automatically])
#   fi
#   ;;
#   */*)
#   if test -d $with_xdr; then
#     XDR_PREFIX=$with_xdr
#   else
#     AC_MSG_WARN([library specified not present. It will be searched automatically])
#   fi
#   ;;
#   *) AC_MSG_ERROR([bad value $with_xdr for --with_xdr. It must be a location or a library file.]);;
# esac
# XDR_LIBS can be "", "*.a" or "-l*"



# this is the case in which we want to use the old library build (libxdrf.a)
AC_MSG_CHECKING([support for xdr library (old one)])
AC_ARG_WITH(xdr,
	[AS_HELP_STRING([--with-xdr=<lib>], [use XDR <lib> with the old format (only libxdrf.a is supported). This option is by default not set as the makefile will make and link by itself the library])],
  [])
case $with_xdr in
  yes)
    AC_MSG_RESULT([yes])
    AC_MSG_NOTICE([manual option selected. The lib will be searched in common locations.])
    with_xdr="yes"
  ;;
  *.a | *.so | *.so.* | *.o)
    AC_MSG_RESULT([yes])
    if test -e $with_xdr; then
      XDR_PREFIX=$(dirname $with_xdr)
      AC_MSG_NOTICE([inserted file/folder does exist. The lib will be checked and tested automatically.])
    else
      AC_MSG_WARN([library specified not present. It will be searched automatically])
    fi
    with_xdr="yes"
  ;;
  */*)
    AC_MSG_RESULT([yes])
    if test -d $with_xdr; then
      XDR_PREFIX=$with_xdr
      AC_MSG_NOTICE([inserted file/folder does exist. The lib will be checked and tested automatically.])
    else
      AC_MSG_WARN([library specified not present. It will be searched automatically])
    fi
    with_xdr="yes"
  ;;
  no)
    AC_MSG_RESULT([no])
    AC_MSG_WARN([xdr library support disabled. It will not be linked])
    with_xdr="no";;
  "")  # here there is also the standard
  AC_MSG_RESULT([no])
  xdr_default="yes"
  with_xdr="yes"
  AC_MSG_NOTICE([the xdr library will be automatically loaded if possible.])
  USE_INTERNAL_XDR="1"
  ;;
  *)
    AC_MSG_RESULT([not valid insertion])
    AC_MSG_ERROR([bad value $with-xdr for --with-xdr. It must be a location or a library file.])
  ;;
esac
# this is connecting the code relative to the old version (the only one supported)
with_xdr_old=$with_xdr
#
#
#
#
# ATTENTION NEW VERSION OF XDR NOT IMPLEMENTED
if test ! "$xdr_default" = "yes"; then
  # ATTENTION libxdrfile.a not implemented
  AC_MSG_WARN([
    !!WARNING!! - manual xdr library insertion (or general system search)
    Support for new versions of the xdr library is not available. It is
    generally advisable to use the built-in library. Indicatively the library
    needed is the old one which is based on libxdrf.a instead of libxdrfile.a.
    In the standard behaviour the Makefile will compile this library
    automatically and link it.])
else
  AC_MSG_WARN([
    xdr library will be compiled and organized in the make step.
    Please consider to compile them separately and test their
    functionalities using the --with-xdr option of the configure.])
fi
if test "$with_xdr_old" = "yes"; then
  flavour_options="$flavour_options xdr_lib"
  # AC_DEFINE([WANT_LIBXDR], [1], [Defined if you want to use xdr lib]) ! not working
  PREPROCFLAGS="${PREPROCFLAGS} -DLINK_XDR"
fi
#
#
#

# look for xdr libraries automatically
# -----------------------------
# AC_MSG_CHECKING([if xdr libraries have been found automatically])
# looking for specific files

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# NOTE: THIS IS NO MORE SUPPORTED AS A NEW VERSION OF THE LIBRARY XDR
# IT IS NOT NEEDED AT ALL (THIS OLD VERSION IS COMPILED ON THE FLY)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# if test "$with_xdr_old" = "no"; then
#
#   XDR_LIB_TOTEST="libxdrfile.a"
#   #if prefix is set by user then search for lib.a under that path
#   if test -n "${XDR_PREFIX}"; then
# 		for XDR_LIB_PATH in lib lib64
# 		do
# 		XDR_LIB_PATH_ABS=${XDR_PREFIX}/${XDR_LIB_PATH}
# 		AC_MSG_NOTICE([search for xdr library in: ${XDR_LIB_PATH_ABS}])
# 		if test -e ${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}; then
# 			XDR_LIBS="-L${XDR_LIB_PATH_ABS}"
# 			XDR_LIBS="${XDR_LIBS} -lxdrfile"
# 			XDR_CFLAGS="-I${XDR_PREFIX}/include/xdrfile/"
# 			AC_MSG_NOTICE([found XDR library in: ${XDR_LIB_PATH_ABS}])
# 			break;
# 		fi
# 		done
#   else
#   	# if prefix not provided, try pkg-config
#   	AC_MSG_NOTICE([No directory was specified for --with-xdr. Trying to find xdr using pkg-config.])
#     if test -n "${PKG_CONFIG}" ; then
#       XDR_CONFIG="${PKG_CONFIG} xdr"
#       XDR_CFLAGS="`${XDR_CONFIG} --cflags`"
#       XDR_LIBS="`${XDR_CONFIG} --libs`"
#     fi
#   	if test -z "$XDR_LIBS"; then
#     	AC_MSG_NOTICE([pkg-config was not able to find the xdr library.])
#     	AC_MSG_NOTICE([Trying some common locations.])
#   		#still not detected, then search some default standard location
#   		for XDR_PREFIX in /usr/local /usr #skip searching/usr since it is standard system path
#   		do
#   			for XDR_LIB_PATH in lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
#   			do
#   			XDR_LIB_PATH_ABS=${XDR_PREFIX}/${XDR_LIB_PATH}
#   			AC_MSG_NOTICE([search for xdr library in: ${XDR_LIB_PATH_ABS}])
#   			if test -e ${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}; then
#   				XDR_LIBS="-L${XDR_LIB_PATH_ABS}"
#   				XDR_LIBS="${XDR_LIBS} -lxdrfile"
#   				XDR_CFLAGS="-I${XDR_PREFIX}/include/xdrfile/"
#   				AC_MSG_NOTICE([found xdr library in: ${XDR_LIB_PATH_ABS}])
#   				break;
#   			fi
#   			done
#   		if test -n "${XDR_LIBS}"; then
#   			break
#   		fi
#   		done
#     fi
#
#     #fail to detect xdr_library, use the standard system path (which can be ommited)
#     if test "${XDR_LIBS}" = ""; then
#       AC_MSG_WARN([fail to detect xdr_library, trying with a macro in the
#         standard system files (compilation not granted).])
#       AC_SEARCH_LIBS(xdrfile_open, [xdrfile], [XDR_LIBS=$LIBS], [])
#     fi
#
#     # still anything found? installing..
#     if test "${XDR_LIBS}" = ""; then
#       AC_MSG_WARN([fail to detect xdr_library, looking for installation files])
#       AC_MSG_CHECKING([if xdr libraries can are internally installed  in campariv3/ext_libs/installed ])
#       a=0
#       if test "$ax_ext_lib" = "yes"; then
#         if test -d "${EXT_LIBS_DIR}/installed/xdr"; then
#           AC_MSG_RESULT([yes])
#           XDR_PREFIX="${EXT_LIBS_DIR}/installed/xdr"
#         else
#           AC_MSG_RESULT([no])
#           AC_MSG_CHECKING([if xdr libraries can be installed in ../ext_libs/installed/xdr ])
#           if test -e "${EXT_LIBS_DIR}/xdrfile*.tar.gz"; then #FAILING HERE!
#             # TODO CHECK THE AUTOMATIC INSTALLATION
#             # AC_MSG_RESULT([yes])
#             XDR_PREFIX="${EXT_LIBS_DIR}/installed/xdr"
#           else
#             a=1
#             AC_MSG_WARN([To install the xdrfile automatically the .tar.gz file is needed in campariv3/ext_libs])
#           fi
#         fi
#       else
#         a=1
#       fi
#       if test a = 1; then
#         AC_MSG_RESULT([no])
#         AC_MSG_ERROR([Impossible to find or install xdr libraries. Please check the
#           installation and eventually try to add --with-xdr=PATH manually.])
#         # THIS IS ONLY A WARNING BECAUSE WE NEED THE libxdrf.a ONLY FOR NOW.
#         # MUST BE SET AS ERROR AFTERWARDS
#       fi
#     fi
#   fi
# fi
# looking for the old library (NOW it is the default)

if test "$with_xdr_old" = "yes" -a "$xdr_default" != "yes"; then
  XDR_LIB_TOTEST="libxdrf.a"
  #if prefix is set by user then search for lib.a under that path
  if test -n "${XDR_PREFIX}"; then
    for XDR_LIB_PATH in x86_64 lib lib64
    do
    XDR_LIB_PATH_ABS=${XDR_PREFIX}/${XDR_LIB_PATH}
    AC_MSG_NOTICE([search for xdr library in: ${XDR_LIB_PATH_ABS}])
    if test -e ${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}; then
      XDR_LIBS="${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}"
      XDR_CFLAGS=""
      AC_MSG_NOTICE([found XDR library in: ${XDR_LIB_PATH_ABS}])
      break;
    fi
    done
  else
    possible_dir="${CAMPARI_HOME}/lib"
    possible_dir2="${CAMPARI_HOME}/source/xdr"
    possible_dir3="$($DIRNAME "${CAMPARI_HOME}")/xdr"
    # if prefix not provided, try pkg-config
    AC_MSG_NOTICE([No directory was specified for --with-xdr. Trying some common locations.])
    # #still not detected, then search some default standard location
    for XDR_PREFIX in $possible_dir $possible_dir2 $possible_dir3 /software/xdr /usr/local /usr #skip searching/usr since it is standard system path
    do
      for XDR_LIB_PATH in "" x86_64 lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
      do
      XDR_LIB_PATH_ABS=${XDR_PREFIX}/${XDR_LIB_PATH}
      AC_MSG_NOTICE([search for xdr library in: ${XDR_LIB_PATH_ABS}])
      if test -e ${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}; then
        XDR_LIBS="${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}"
        XDR_CFLAGS=""
        AC_MSG_NOTICE([found xdr library in: ${XDR_LIB_PATH_ABS}])
        break;
      fi
      done
    if test -n "${XDR_LIBS}"; then
      break
    fi
    done
  fi
  if test "${XDR_LIBS}" = ""; then
    AC_MSG_WARN([Unable to find old XDR library. Please consider properly --with-xdr=<lib>])
  fi
fi

xdr_check_default=$xdr_default
# NB: this check should be enforced enable-xdr-check=yes
AC_MSG_CHECKING([whether Fortran compiler is checked for XDR support])
AC_ARG_ENABLE([xdr-check],
  [AS_HELP_STRING([--disable-xdr-check],
    [(dev opt) disable check of xdr library support in Fortran compiler. ])],
    [],[enable_xdr_check=$xdr_check_default]) #ATTENTION TO modify TO YES
test "x$enable_xdr_check" = xno || enable_xdr_check=yes
AC_MSG_RESULT($enable_xdr_check)
if test "${enable_xdr_check}" = "yes" -a "${xdr_default}" = "yes"; then
  AC_MSG_WARN([
    The checks for xdr library are not available for the default option
    of compiling it on make run. The checks will not be considered.
    In case of malfunctioning, please compile separately the xdr library
    (e.g. using the make libxdrf.a in the source directory) and add the
    main lib (or source/xdr) directory to --with-xdr. This will run
    all the checks. ])
  enable_xdr_check="no"
fi
# now checking xdr_ftest
# ----------------------
# if test "x$enable_xdr_check" = xyes; then
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # THIS IS THE NEW XDR --------------------------no more needed-----
# -------------------------------------------------------------------
# -------------------------------------------------------------------
#   if test "$with_xdr_old" = "no"; then
#     # checking the header file
#     AC_MSG_CHECKING([whether it is possible to link xdrfile_open function])
#     AC_TRY_LINK_FUNC(xdrfile_open,[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
#     # checking the example tests (can be more than one, e.g. trr-support)
#     AC_MSG_CHECKING([if Fortran compiler supports xdr libraries ])
#
#     # NB: this examples are not working with more stringent compiler options
#     # doit='$FC -c ${FCFLAGS} lib_test/xdr.F90'
#     # AC_TRY_EVAL(doit)
#     # doit='$FC -c ${FCFLAGS} lib_test/example/read-and-write-xtc/read-and-write-xtc.f90'
#     # AC_TRY_EVAL(doit)
#     # doit='$FC ${FCFLAGS} *.o  ${XDR_LIBS} ${LIBS}'
#     LIBS_TMP="${LIBS}"
#     CFLAGS_TMP="${CFLAGS}"
#     AC_LANG_PUSH([C++])
#       LIBS="$LIBS $XDRLIBS" # it is not this
#       CFLAGS="$CFLAGS $XDR_CFLAGS"
#       AC_CHECK_HEADER(xdrfile.h, [a=0], [a=1])
#     AC_LANG_POP([C++])
#     LIBS="${LIBS_TMP}"
#     CFLAGS="${CFLAGS_TMP}"
#
#     if test $a = 1; then
#       AC_MSG_WARN([Could not find or include xdrfile.h])
#     fi
#     # checking compilation
#     AC_MSG_NOTICE([trying compilation using:
# ${FC} -c ${FCFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}])
#     AC_MSG_CHECKING([whether selected compiler and compiler flags can compile])
#     cat <<EOF >conftest.f90
# program testxtc
#
# parameter (maxatom=10000,maxx=3*maxatom)
# integer xd,xd2,natoms,step,ret,i
# real    time,box(9),x(maxx)
#
# call xdrfopen(xd,"lib_test/test.xtc","r",ret)
# print *,'opened lib_test/test.xtc, ret=',ret
# call xdrfopen(xd2,"lib_test/testout.xtc","w",ret)
# print *,'opened lib_test/testout.xtc, ret=',ret
#
# 10 call readxtc(xd,natoms,step,time,box,x,prec,ret)
#
# if ( ret .eq. 1 ) then
#   call writextc(xd2,natoms,step,time,box,x,prec,ret)
# else
#   print *,'Error reading xtc'
# endif
#
# stop
# end
# EOF
# # Strangely this system is compiling but the xdrfopen extention should be invisible.
# # when using libxdrfile.a. PLEASE CHANGE THE TEST ONCE CAMPARI HAS BEEN UPDATED
# # FOR THE NEW VERSION OF XDR
#     doit='${FC} -c ${FCFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}'
#     if AC_TRY_EVAL(doit); then
#       AC_MSG_RESULT([yes])
#       needed_options="$needed_options xdr"
#       AC_DEFINE([HAVE_LIBXDR], [1], [Defined if you have LIBXDR support])
#     else
#       AC_MSG_RESULT([no])
#       AC_MSG_WARN([Unable to compile conftest.f90 for xdr library support.
#         XDR checks unsuccesful. The library has not been linked properly.
#         Please check it and eventually try to add --with-xdr=PATH manually.])
#     fi

a=2
# now check for old version already installed
# --------------------------------------------
if test "x$enable_xdr_check" = xyes; then
  # checking the header file
  AC_MSG_CHECKING([whether it is possible to link xdrfopen function])
  AC_TRY_LINK_FUNC(xdrfopen,[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
  # checking the example tests (can be more than one, e.g. trr-support)
  AC_MSG_NOTICE([finding header file xdrf.h])
  LIBS_TMP="${LIBS}"
  CFLAGS_TMP="${CFLAGS}"
  AC_LANG_PUSH([C++])
    LIBS="${LIBS} ${XDR_LIBS}" # it is not this
    CFLAGS="${CFLAGS} ${XDR_CFLAGS}"
    AC_CHECK_HEADER(${XDR_LIB_PATH_ABS}/xdrf.h, [a=0], [a=1])
  AC_LANG_POP([C++])
  LIBS="${LIBS_TMP}"
  CFLAGS="${CFLAGS_TMP}"
  if test $a = 1; then
    AC_MSG_WARN([Could not find or include xdrf.h])
  fi

  # checking compilation
  AC_MSG_NOTICE([trying compilation using:
${FC} -c ${FCFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}])
  AC_MSG_CHECKING([whether selected compiler and compiler flags can compile])
  cat <<EOF >conftest.f90
program testxtc

   parameter (maxatom=10000,maxx=3*maxatom)
   integer xd,xd2,natoms,step,ret,i
   real    time,box(9),x(maxx)

   call xdrfopen(xd,"lib_test/test.xtc","r",ret)
   print *,'opened lib_test/test.xtc, ret=',ret
   call xdrfopen(xd2,"lib_test/testout.xtc","w",ret)
   print *,'opened lib_test/testout.xtc, ret=',ret

10 call readxtc(xd,natoms,step,time,box,x,prec,ret)

   if ( ret .eq. 1 ) then
      call writextc(xd2,natoms,step,time,box,x,prec,ret)
   else
      print *,'Error reading xtc'
   endif

   stop
end
EOF
  doit='${FC} -c ${FCFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}'
  if AC_TRY_EVAL(doit); then
    AC_MSG_RESULT([yes])
   	needed_options="$needed_options xdr_old"
  else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([Unable to compile conftest.f90 for xdr library support.
      XDR checks unsuccesful. The library has not been linked properly.
      Please check it and eventually try to add --with-xdr-old=PATH manually.
      By default the library will be compiled and linked during the makefile run])
  fi
else
 	needed_options="$needed_options xdr library"
  AC_MSG_WARN([The abolishment of xdr checks are not safe and usually are restricted to dev.])
fi

# -------------------------------
#       Checking hsl libs
# -------------------------------
SECTION_TITLE([Checking for hsl libraries])



# this is the case in which we want to use the old library build (libxdrf.a)
AC_MSG_CHECKING([support for hsl library])
AC_ARG_WITH(hsl,
	[AS_HELP_STRING([--with-hsl=<lib>], [use hsl <lib>. This option is by default not set as the makefile will make and link by itself the library])],
  [])
case $with_hsl in
  yes)
    AC_MSG_RESULT([yes])
    AC_MSG_NOTICE([manual option selected. The lib will be searched in common locations.])
    with_hsl="yes";;
  *.a | *.so | *.so.* | *.o)
    AC_MSG_RESULT([yes])
    if test -e $with_hsl; then
      HSL_PREFIX=$(dirname $with_hsl)
      AC_MSG_NOTICE([inserted file/folder does exist. The lib will be checked and tested automatically.])
    else
      AC_MSG_WARN([library specified not present. It will be searched automatically])
    fi
    with_hsl="yes";;
  */*)
    AC_MSG_RESULT([yes])
    if test -d $with_hsl; then
      HSL_PREFIX=$with_hsl
      AC_MSG_NOTICE([inserted file/folder does exist. The lib will be checked and tested automatically.])
    else
      AC_MSG_WARN([library specified not present. It will be searched automatically])
    fi
    with_hsl="yes";;
  no)
    AC_MSG_RESULT([no])
    AC_MSG_WARN([hsl library support disabled. It will not be linked])
    with_hsl="no";;
  "")  # here there is also the standard
    AC_MSG_RESULT([no])
    hsl_default="yes"
    with_hsl="yes"
    USE_INTERNAL_HSL=1
    AC_MSG_NOTICE([the hsl library will be automatically loaded if possible.]);;
  *)
    AC_MSG_RESULT([not valid insertion])
    AC_MSG_ERROR([bad value $with-hsl for --with-hsl. It must be a location or a library file.]);;
esac
if test ! "$hsl_default" = "yes"; then
  AC_MSG_WARN([

    !!WARNING!! - manual hsl library insertion (or general system search)
    Indicatively the library needed is searched using libhsl.a as
    reference. In the standard behaviour the Makefile will compile this
    library automatically and link it.
])
else
  AC_MSG_WARN([
    hsl library will be compiled and organized in the make step.
    Please consider to compile them separately and test their
    functionalities using the --with-hsl option of the configure.])
fi
if test "$with_hsl" = "yes"; then
  flavour_options="$flavour_options hsl_lib"
  # AC_DEFINE([WANT_LIBHSL], [1], [Defined if you want to use hsl lib])
  PREPROCFLAGS="${PREPROCFLAGS} -DLINK_HSL"
fi


if test "$with_hsl" = "yes" -a "$hsl_default" != "yes"; then
  HSL_LIB_TOTEST="libhsl.a
"
  #if prefix is set by user then search for lib.a under that path
  if test -n "${HSL_PREFIX}"; then
    for HSL_LIB_PATH in x86_64 lib lib64
    do
    HSL_LIB_PATH_ABS=${HSL_PREFIX}/${HSL_LIB_PATH}
    AC_MSG_NOTICE([search for hsl library in: ${HSL_LIB_PATH_ABS}])
    if test -e ${HSL_LIB_PATH_ABS}/${HSL_LIB_TOTEST}; then
      HSL_LIBS="${HSL_LIB_PATH_ABS}/${HSL_LIB_TOTEST}"
      HSL_CFLAGS=""
      AC_MSG_NOTICE([found HSL library in: ${HSL_LIB_PATH_ABS}])
      break;
    fi
    done
  else
    possible_dir="${CAMPARI_HOME}/lib"
    possible_dir2="${CAMPARI_HOME}/source/hsl"
    possible_dir3="$($DIRNAME "${CAMPARI_HOME}")/hsl"
    # if prefix not provided, try pkg-config
    AC_MSG_NOTICE([No directory was specified for --with-hsl. Trying some common locations.])
    # #still not detected, then search some default standard location
    for HSL_PREFIX in $possible_dir $possible_dir2 $possible_dir3 /software/hsl /usr/local /usr #skip searching/usr since it is standard system path
    do
      for HSL_LIB_PATH in "" x86_64 lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
      do
      HSL_LIB_PATH_ABS=${HSL_PREFIX}/${HSL_LIB_PATH}
      AC_MSG_NOTICE([search for hsl library in: ${HSL_LIB_PATH_ABS}])
      if test -e ${HSL_LIB_PATH_ABS}/${HSL_LIB_TOTEST}; then
        HSL_LIBS="${HSL_LIB_PATH_ABS}/${HSL_LIB_TOTEST}"
        HSL_CFLAGS=""
        AC_MSG_NOTICE([found hsl library in: ${HSL_LIB_PATH_ABS}])
        break;
      fi
      done
    if test -n "${HSL_LIBS}"; then
      break
    fi
    done
  fi
  if test "${HSL_LIBS}" = ""; then
    AC_MSG_WARN([Unable to find HSL library. Please consider properly --with-hsl=<lib>])
  fi
fi

AC_MSG_WARN([The checks for hsl library are not available.])


# -------------------------------
#       Checking FFTW3 libs
# -------------------------------
SECTION_TITLE([Checking for FFTW3 libraries])

FFTW3_PREFIX=""
FFTW3_LIBS=""

# with-fftw3 ?
# -------------
AC_ARG_WITH(fftw3,
	[AS_HELP_STRING([--with-fftw3=<lib>], [use FFTW3 library <lib>])],
  [])
case $with_fftw3 in
  yes | "") with_fftw3=yes ;;
  *.a | *.so | *.so.* | *.o)
  if test -e $with_fftw3; then
    FFTW3_PREFIX=$($DIRNAME "${with_fftw3}")
  else
    AC_MSG_WARN([the library specified was not present in the specified directory. It will be searched automatically])
  fi
  FFTW3_LIB_TOTEST_IT=$($BASENAME "${with_fftw3}") # I don't really need it
  with_fftw3=yes
  AC_MSG_WARN([
  The insertion of a library file (.a) will not use directly
  the inserted file name but it will use the directory to find
  the usual library files (libfftw3.a and libfftw3_omp/threads.a if
  threads are enabled). No system links will be added in this way.])
  ;;
  */*)
  if test -d $with_fftw3; then
    FFTW3_PREFIX=$with_fftw3
    with_fftw3=yes
  else
    AC_MSG_WARN([library specified not present. It will be searched automatically])
  fi
  ;;
  no) AC_MSG_WARN([Every fftw3 functionality will not be used. CAMPARI will have reduced functionalities.]);;
  *) AC_MSG_ERROR([bad value $with_fftw3 for --with-FFTW3. It must be a location or a library file.]);;
esac
# FFTW3_LIBS can be "", "*.a" or "-l

# look for FFTW3 libraries automatically
# -----------------------------
# AC_MSG_CHECKING([if FFTW3 libraries have been found automatically])
# looking for specific files
if test "${with_fftw3}" != "no"; then
  # if lib.a was not inserted
  a=2
  if test -z "$FFTW3_PREFIX" -o -n "$FFTW3_LIB_TOTEST_IT"; then
    FFTW3_LIBS_TMP="-lfftw3"
    if test "$threads" = "yes"; then
      FFTW3_LIBS_TMP="-lfftw3 -lm"
    fi
  else
    FFTW3_LIBS_TMP=
  fi

  FFTW3_LIB_TOTEST_IT="libfftw3.a"
  # we need to have two dependencies with multi-threads
  if test "$threads" = "yes"; then
    FFTW3_LIB_TOTEST_IT="${FFTW3_LIB_TOTEST_IT} libfftw3_omp.a libfftw3_threads.a" # array
    a=0
  fi

  # if prefix not provided, try pkg-config (if no direct lib.a has been inserted)
  if test -n "${FFTW3_LIBS_TMP}" -a -z "${FFTW3_PREFIX}"; then
    AC_MSG_NOTICE([No directory was specified for --with-fftw3. Trying to find fftw3 using pkg-config.])
    if ! test -z "${PKG_CONFIG}" ; then
      FFTW3_CONFIG="${PKG_CONFIG} fftw3"
      FFTW3_CFLAGS="`${FFTW3_CONFIG} --cflags`" # to change to fortran flags
      FFTW3_LIBS="`${FFTW3_CONFIG} --libs`"
      # if test "$threads" = "yes"; then
      AC_MSG_NOTICE([multi-threaded fftw3 pkg-config search not available])
      FFTW3_LIBS=""
      # fi
    else
      AC_MSG_NOTICE([pkg-config was not able to find the FFTW3 library.])
    fi
  fi


  if test "$FFTW3_LIBS" = ""; then
    # if you specified directly the lib.a it will be looked into common locations
    #if prefix is set by user then search for lib.a under that path
    if test -z "${FFTW3_PREFIX}"; then
      AC_MSG_NOTICE([No directory was specified for --with-fftw3. Looking for it in usual locations])
      FFTW3_PREFIX=""
    fi
    AC_MSG_NOTICE([Looking for the library directories and files])
    for FFTW3_LIB_TOTEST in $FFTW3_LIB_TOTEST_IT; do
  	  for FFTW3_PREFIX in "${FFTW3_PREFIX}" /usr /usr/local; do
    		for FFTW3_LIB_PATH in "" lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu; do
    		FFTW3_LIB_PATH_ABS=${FFTW3_PREFIX}/${FFTW3_LIB_PATH}
    		AC_MSG_NOTICE([search for ${FFTW3_LIB_TOTEST} file in: ${FFTW3_LIB_PATH_ABS}])
    		if test -e ${FFTW3_LIB_PATH_ABS}/${FFTW3_LIB_TOTEST}; then
          if test "$a" != "2"; then
            a=$((a+1)) # increment due to the fact that I need ALSo libfftw3
          fi
          if test -n "${FFTW3_LIBS_TMP}"; then
            if test "${FFTW3_LIB_TOTEST}" = "libfftw3_omp.a"; then
              FFTW3_LIBS_TMP="${FFTW3_LIBS_TMP} -lfftw3_omp"
            elif test "${FFTW3_LIB_TOTEST}" = "libfftw3_threads.a"; then
              FFTW3_LIBS_TMP="${FFTW3_LIBS_TMP} -lfftw3_threads"
            fi
            if test "${a}" = "2"; then
              FFTW3_LIBS="${FFTW3_LIBS} ${FFTW3_LIBS_TMP}"
              FFTW3_LIBS="${FFTW3_LIBS}  -L${FFTW3_LIB_PATH_ABS}"
              FFTW3_CFLAGS="${FFTW3_CFLAGS} -I`$ECHO "${FFTW3_PREFIX}" | $SED 's|lib.*||'`/include/"
            fi
          else
            FFTW3_LIBS="${FFTW3_LIBS} ${FFTW3_LIB_PATH_ABS}/${FFTW3_LIB_TOTEST}"
          fi
    			AC_MSG_NOTICE([found ${FFTW3_LIB_TOTEST} file in: ${FFTW3_LIB_PATH_ABS}])
    			break;
    		fi
    		done
        if test -n "${FFTW3_LIBS}" ; then
          if test "$a" = "1" -o "$a" = "2"; then
            break
          fi
        fi
      done
      if test -n "${FFTW3_LIBS}" -a "$a" = "2"; then
        break
      fi
    done
  fi

  # checking the include directory
  AC_ARG_VAR(fftw3_includedir, [include directory of fftw3 library])
  if test -z "${fftw3_includedir}"; then
    fftw3_includedir="`$ECHO "${FFTW3_LIB_PATH_ABS}" | $SED 's|lib.*||'`/include/"
  fi
  AC_MSG_CHECKING([whether the include directory is correctly set. Looking for fftw3.h])
  if test -e "${fftw3_includedir}/fftw3.h"; then
    AC_MSG_RESULT([found])
    FFTW3_CFLAGS="${FFTW3_CFLAGS} -I${fftw3_includedir}"
  else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([Please specify the include directory of the fftw3 library])
  fi

  #fail to detect FFTW3_library, use the standard system path (which can be ommited)
  if test -z "${FFTW3_LIBS}"; then
    if test "${multi-thread}" = "yes"; then
      AC_MSG_WARN([fail to detect fftw3 library with multi_threading, looking the standard system path using function: fftw_init_threads])
      AC_SEARCH_LIBS(fftw_init_threads, [FFTW3file], [FFTW3_LIBS=$LIBS], [])
    else
      AC_MSG_WARN([fail to detect FFTW3_library, looking the standard system path using function: fftw_alloc_complex])
      AC_SEARCH_LIBS(fftw_alloc_complex, [FFTW3file], [FFTW3_LIBS=$LIBS], [])
    fi
  fi

  # still anything found? installing..
  if test -z "${FFTW3_LIBS}"; then
    AC_MSG_CHECKING([if FFTW3 libraries can be installed automatically from ext_libs])
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Even if the fftw3 option was active, no library was found.])
  fi
fi

# NB: this check should be enforced enable-FFTW3-check=yes
# ------------------------------
AC_MSG_CHECKING([whether Fortran compiler is checked for FFTW3 support])
AC_ARG_ENABLE([fftw3-check],
  [AS_HELP_STRING([--disable-fftw3-check],
    [(dev opt) disable check of fftw3 library support in Fortran compiler. ])],
    [],[enable_fftw3_check=yes])
# if you disabled the library no check follows
if test "${with_fftw3}" = "no"; then
  enable_fftw3_check=no
  FFTW3_LIBS=
  FFTW3_CFLAGS=
fi
AC_MSG_RESULT($enable_fftw3_check)

# now checking FFTW3_ftest
# ----------------------
if test "x$enable_fftw3_check" = xyes; then
  # checking the header file

  if test -n "${FFTW3_LIBS}"; then
    # checking threads support for fftw3
    if test "${threads}" = "yes"; then
      AC_MSG_NOTICE([trying compilation for multi-threaded fftw3 using:
${FC} ${FCFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} ${THREADFLAGS} conftest.f90 ${FFTW3_LIBS}])
      AC_MSG_CHECKING([whether multi-thread support is available for fftw3])
        cat <<EOF >conftest.f90
program test
  implicit none
  include 'fftw3.f'

  integer :: n, i
  integer :: plan
  complex, allocatable :: in(:), out(:)
  integer iret
  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(${n_threads})

  write(*,*) 'Input data:'
  n = 10
  allocate(in(n))
  allocate(out(n))
  do i=1,n
    in(i) = cmplx(i,0.0)
    write(*,*) in(i)
  enddo

  ! Forward Fourier transform
  call dfftw_plan_dft_1d(plan,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

  write(*,*) 'Fourier transform of the input data:'
  do i=1,n
    write(*,*) out(i)
  enddo

  ! Inverse Fourier transform
  call dfftw_plan_dft_1d(plan,n,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

  write(*,*) 'Recovered data from inverse Fourier transform:'
  do i=1,n
    write(*,*) real(in(i)/n)
  enddo
end program test
EOF
      doit='${FC} ${FCFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} ${THREADFLAGS} conftest.f90 ${FFTW3_LIBS}'
      if AC_TRY_EVAL(doit); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Unable to compile fftw3 with threads])
      fi

    # standard NO multi-threaded fft3
    else
      AC_MSG_NOTICE([trying compilation for fftw3 using:
${FC} ${FCFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} ${THREADFLAGS} conftest.f90 ${FFTW3_LIBS}])
      AC_MSG_CHECKING([whether fftw3 is supported by specific compilers])
        cat <<EOF >conftest.f90
program test
  implicit none
  include 'fftw3.f'

  integer :: n, i
  integer :: plan
  complex, allocatable :: in(:), out(:)
  write(*,*) 'Input data:'
  n = 10
  allocate(in(n))
  allocate(out(n))
  do i=1,n
    in(i) = cmplx(i,0.0)
    write(*,*) in(i)
  enddo

  ! Forward Fourier transform
  call dfftw_plan_dft_1d(plan,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

  write(*,*) 'Fourier transform of the input data:'
  do i=1,n
    write(*,*) out(i)
  enddo

  ! Inverse Fourier transform
  call dfftw_plan_dft_1d(plan,n,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

  write(*,*) 'Recovered data from inverse Fourier transform:'
  do i=1,n
    write(*,*) real(in(i)/n)
  enddo
end program test
EOF
      doit='${FC} ${FCFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} conftest.f90 ${FFTW3_LIBS}'
      if AC_TRY_EVAL(doit); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Unable to compile fftw3.])
      fi
    fi
  # FFTW3_LIBS neither assigned - bug check
  else
    AC_MSG_ERROR([Variable FFTW3_LIBS must be assigned to compile FFTW3 dependent source code.])
  fi
else
  AC_MSG_WARN([The abolishment of FFTW3 checks are not safe and usually are restricted to dev.])
fi
if test "${with_fftw3}" != "no"; then
 	flavour_options="$flavour_options FFTW3"
  PREPROCFLAGS="${PREPROCFLAGS} -DLINK_FFTW"
fi
# ----------------------------------
# ----------------------------------
# final testing of general compilation (if libs have diff target compiler)
# if test "$compiler" != "no"; then
#   AC_MSG_NOTICE([trying compilation using:
# ${FC} -c ${FCFLAGS} ${NETCDF4_FFLAGS}  ${FCFLAGS_f90} conftest.f90 ${NETCDF4_FLIBS}])
#   AC_MSG_CHECKING([whether selected compiler and compiler flags can compile])
#   cat <<EOF >conftest.f90
# module conftest_module
# USE ISO_C_BINDING
# use netcdf
# end module conftest_module
# EOF
#   doit='${FC} -c ${FCFLAGS} ${NETCDF4_FFLAGS} ${FCFLAGS_f90} conftest.f90 ${NETCDF4_FLIBS}'
#   if AC_TRY_EVAL(doit); then
#     AC_MSG_RESULT([yes])
#   else
#     AC_MSG_RESULT([no])
#     AC_MSG_ERROR([Unable to compile with the selected compiler: ${compiler}.])
#   fi
#   FLIBS=
#   FCLIBS=
# fi
# ----------------------------------
# ----------------------------------
# ----------------------------------
# ----------------------------------


# ----------------------------------
#       some last cleanings
# ----------------------------------
if test -e "a.out"; then
  $RM "a.out"
fi
if test -e "threads_assignments_conftest"; then
  $RM "threads_assignments_conftest"
fi


# ----------------------------------
#         Printing results
# ----------------------------------

dnl Now process the options strings. Essentially, we want two lists
dnl one for the options present (which we have) and one for the options
dnl missing (which we don't)

SECTION_TITLE([Configuration results])

echo "Needed options:" >& AS_MESSAGE_FD
echo "$needed_options" >& AS_MESSAGE_FD
echo >& AS_MESSAGE_FD

echo "Flavour options:" >& AS_MESSAGE_FD
echo "${flavour_options}" >& AS_MESSAGE_FD
echo >& AS_MESSAGE_FD

echo "Missing options:" >& AS_MESSAGE_FD
echo "$flavour_options" "$all_flavour_options" | tr ' ' '\n' | sort | uniq -u | tr '\n' ' ' >& AS_MESSAGE_FD
echo >& AS_MESSAGE_FD
# echo >& AS_MESSAGE_FD

# ----------------------------------
#        Var substitutions
# ----------------------------------
echo " " >& AS_MESSAGE_FD
SECTION_TITLE([Important variables])

# general settings
# ---------------------------------------------------------------------------
AC_SUBST(CAMPARI_HOME)
AC_SUBST(ARCH)
echo "CAMPARI_HOME=$CAMPARI_HOME" >& AS_MESSAGE_FD
echo "ARCH=$ARCH" >& AS_MESSAGE_FD
echo "SYSTEM_OS=${system_os}" >& AS_MESSAGE_FD
echo '--' >& AS_MESSAGE_FD

# general compiler
# ---------------------------------------------------------------------------
AC_SUBST(FC)
AC_SUBST(CC)
echo "FC=$FC" >& AS_MESSAGE_FD
echo "CC=$CC" >& AS_MESSAGE_FD
echo '--' >& AS_MESSAGE_FD
# echo "F77=$F77" >& AS_MESSAGE_FD

# compiler flags
# ---------------------------------------------------------------------------
AC_SUBST(FCFLAGS)
AC_SUBST(FFLAGS)
AC_SUBST(PREPROCFLAGS)
AC_SUBST(AUTOVEC)
AC_SUBST(OPENLOOPS)
AC_SUBST(OPTIMIZATION_LEVEL)
echo "FCFLAGS=$FCFLAGS" >& AS_MESSAGE_FD
# echo "FFLAGS=$FFLAGS" >& AS_MESSAGE_FD # automatically set to -g -O2 but not used
# echo "FCLAGS=$FCLAGS" >& AS_MESSAGE_FD # this does NOT exist
# echo "FCFLAGS_f90=$FCFLAGS_f90" >& AS_MESSAGE_FD
# echo "CFLAGS=$CFLAGS" >& AS_MESSAGE_FD
# echo "CFFLAGS=$CFFLAGS" >& AS_MESSAGE_FD
# echo "SHLIB_FFLAGS=$SHLIB_FFLAGS" >& AS_MESSAGE_FD
# echo "FPICFLAGS=$FPICFLAGS" >& AS_MESSAGE_FD
# echo "SAFE_FFLAGS=$SAFE_FFLAGS" >& AS_MESSAGE_FD
echo "PREPROCFLAGS=$PREPROCFLAGS" >& AS_MESSAGE_FD
echo "AUTOVEC=$AUTOVEC" >& AS_MESSAGE_FD
echo "OPENLOOPS=$OPENLOOPS" >& AS_MESSAGE_FD
echo "OPTIMIZATION_LEVEL=$OPTIMIZATION_LEVEL" >& AS_MESSAGE_FD
echo '--' >& AS_MESSAGE_FD

# special Flags
# ---------------------------------------------------------------------------
AC_SUBST(WARNFLAGS)
AC_SUBST(DBFLAGS)
if test "${warnflags}" = "yes"; then
  echo "WARNFLAGS=$WARNFLAGS" >& AS_MESSAGE_FD
  g=1
fi
if test "${debug}" = "yes"; then
  echo "DBFLAGS=$DBFLAGS" >& AS_MESSAGE_FD
  g=1
fi
if test "$g" = "1"; then echo '--' >& AS_MESSAGE_FD; fi


# multi_threading
# ---------------------------------------------------------------------------
AC_SUBST(THREADFLAGS)
if test "$threads" = "yes"; then
  echo "OPENMP_CFLAGS=$OPENMP_CFLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_CXXFLAGS=$OPENMP_CXXFLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_F77FLAGS=$OPENMP_F77FLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_FFLAGS=$OPENMP_FFLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_FCLAGS=$OPENMP_FCLAGS" >& AS_MESSAGE_FD
  echo "THREADFLAGS=$THREADFLAGS" >& AS_MESSAGE_FD
  echo '--' >& AS_MESSAGE_FD
fi

# mpi-related
# ---------------------------------------------------------------------------
AC_SUBST(HAVE_MPI)
AC_SUBST(MPIFF)
AC_SUBST(MPIFFLAGS)
AC_SUBST(MPILIBS)

if test "$mpi" != "no"; then
  echo "MPICC=$MPICC" >& AS_MESSAGE_FD
  echo "MPICXX=$MPICXX" >& AS_MESSAGE_FD
  echo "MPIF77=$MPIF77" >& AS_MESSAGE_FD
  echo "MPIFC=$MPIFC" >& AS_MESSAGE_FD
  echo "MPIFF(this will be substituted)=$MPIFF" >& AS_MESSAGE_FD
  echo "wrapper_fcflags(indicate the mpi layer on the compiler wrapper - only compiler flags)=$wrapper_fcflags" >& AS_MESSAGE_FD
  echo "MPILIBS(it comes from the wrapper plus specific libs)=$MPILIBS" >& AS_MESSAGE_FD
  echo '--' >& AS_MESSAGE_FD
fi

# netcdf4
# ---------------------------------------------------------------------------
AC_SUBST(HAVE_NETCDF4)
AC_SUBST(NETCDF4_VERSION)
# AC_SUBST(NETCDF4_CC)
# AC_SUBST(NETCDF4_CFLAGS)
# AC_SUBST(NETCDF4_CPPFLAGS)
AC_SUBST(NETCDF4_LDFLAGS)
AC_SUBST(NETCDF4_LIBS)
AC_SUBST(NETCDF4_FC)
AC_SUBST(NETCDF4_FFLAGS)
AC_SUBST(NETCDF4_FLIBS)
if test "$with_netcdf4_fortran" = "yes"; then
  # echo "HAVE_NETCDF4=1" >& AS_MESSAGE_FD
	echo "NETCDF4_VERSION=$NETCDF4_VERSION" >& AS_MESSAGE_FD
	echo "NETCDF4_FC=$NETCDF4_FC" >& AS_MESSAGE_FD
	echo "NETCDF4_FFLAGS=$NETCDF4_FFLAGS" >& AS_MESSAGE_FD
	# echo "NETCDF4_FCLAGS=$NETCDF4_FCLAGS" >& AS_MESSAGE_FD # it is not set
  echo "NETCDF4_LDFLAGS(not used)=$NETCDF4_LDFLAGS" >& AS_MESSAGE_FD
  echo "NETCDF4_FLIBS=$NETCDF4_FLIBS" >& AS_MESSAGE_FD
  echo '--' >& AS_MESSAGE_FD
fi

# lapack and blas
# ---------------------------------------------------------------------------
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
echo "BLAS_LIBS=$BLAS_LIBS" >& AS_MESSAGE_FD
echo "LAPACK_LIBS=$LAPACK_LIBS" >& AS_MESSAGE_FD
echo '--' >& AS_MESSAGE_FD

# general libs
# ---------------------------------------------------------------------------
AC_SUBST(LIBS)
AC_SUBST(FLIBS)
AC_SUBST(FCLIBS)
echo "LIBS=$LIBS" >& AS_MESSAGE_FD
echo "FLIBS=$FLIBS" >& AS_MESSAGE_FD
echo "FCLIBS=$FCLIBS" >& AS_MESSAGE_FD
echo '--' >& AS_MESSAGE_FD

# fftw3 lib
# ---------------------------------------------------------------------------
AC_SUBST(HAVE_FFTW3)
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)
if test "${with_fftw3}" != "no"; then
  echo "FFTW3_CFLAGS=$FFTW3_CFLAGS" >& AS_MESSAGE_FD
  echo "FFTW3_LIBS=$FFTW3_LIBS" >& AS_MESSAGE_FD
  echo '--' >& AS_MESSAGE_FD
fi

# xdr lib
# ---------------------------------------------------------------------------
AC_SUBST(USE_INTERNAL_XDR)
AC_SUBST(XDR_CFLAGS)
AC_SUBST(XDR_LIBS)
if test "$with_xdr" = "yes"; then
  if test "$xdr_default" != "yes"; then
    echo "XDR_CFLAGS=$XDR_CFLAGS" >& AS_MESSAGE_FD
    echo "XDR_LIBS=$XDR_LIBS" >& AS_MESSAGE_FD
    echo '--' >& AS_MESSAGE_FD
  else
    echo 'The xdr library will be set and linked during the make step.' >& AS_MESSAGE_FD
    echo '--' >& AS_MESSAGE_FD
  fi
fi
# hsl lib
# ---------------------------------------------------------------------------
AC_SUBST(USE_INTERNAL_HSL)
AC_SUBST(HSL_CFLAGS) # blas
AC_SUBST(HSL_LIBS) # blas
if test "$with_hsl" = "yes"; then
  if test "$hsl_default" != "yes"; then
    echo "HSL_CFLAGS=$HSL_CFLAGS" >& AS_MESSAGE_FD
    echo "HSL_LIBS=$HSL_LIBS" >& AS_MESSAGE_FD
    echo '--' >& AS_MESSAGE_FD
  else
    echo 'The hsl library will be set and linked during the make step.' >& AS_MESSAGE_FD
    echo '--' >& AS_MESSAGE_FD
  fi
  echo "" >& AS_MESSAGE_FD
fi


AC_MSG_NOTICE([generating makefiles])

dnl Process Makefile.in to create Makefile
AC_CONFIG_FILES([Makefile])

if test "$to_go_back_dir" != "${CAMPARI_HOME}/source/"; then
  AC_MSG_NOTICE([setting back original directory])
  doit="$CD ${to_go_back_dir}"
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([directory sucessfully changed to previous value: ${to_go_back_dir}.])
  fi
fi

AC_OUTPUT
