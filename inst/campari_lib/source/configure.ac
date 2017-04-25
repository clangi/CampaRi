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
echo -n AC_PACKAGE_VERSION>VERSION
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

# Find out about the host we're building on.
#AC_CANONICAL_HOST

# Find out about the target we're building for.
#AC_CANONICAL_TARGET

# Check for the existence of this file before proceeding.
# AC_CONFIG_SRCDIR([parsekey.f90])






# ----------------------
#   Compilers checks
# ----------------------

# NB on AM_PROG_CC_C_O: Test whether the Fortran compiler accepts the options -c and -o simultaneously, and define F77_NO_MINUS_C_MINUS_O or FC_NO_MINUS_C_MINUS_O, respectively, if it does not.

# Find fortran compiler. If FC isn't set after AC_PROG_FC, we couldn't
# find a fortran compiler. If the user requested the fortran API, and
# has told us not to recover from a missing compiler, then bail out
# right here.
AC_MSG_NOTICE([finding standard Fortran compiler])
# AC_LANG(Fortran)
AC_PROG_FC
AC_PROG_F77($FC)
AC_PROG_INSTALL # failing because no install-sh
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP
AC_PATH_PROG(WHICH, which)
AC_PATH_PROG(BASH, bash)
AC_PATH_PROG(RM, rm)
AC_PATH_PROG(MV, mv)
AC_PATH_PROG(MAKE, make)
# why this is impossible to find???
AC_PATH_PROG(CD, cd)
if test "$CD" = "no" -o "$CD" = ""; then
  CD="cd"
fi
AC_PATH_PROG(TAR, tar)
AC_PATH_PROG(MKDIR, mkdir)
AC_PATH_PROG(PWD, pwd)
AC_PATH_PROG(DIRNAME, dirname)
AC_PATH_PROG(PKG_CONFIG, pkg-config)



# Set output variable INSTALL to the name of a BSD-compatible install program,
# if one is found in the current PATH. Otherwise, set INSTALL to ‘dir/install-sh -c’,
# checking the directories specified to AC_CONFIG_AUX_DIR
# (or its default directories) to determine dir (see Output).Also set the variables INSTALL_PROGRAM and INSTALL_SCRIPT to ‘${INSTALL}’ and INSTALL_DATA to ‘${INSTALL} -m 644’.



# Setting up CAMPARI_HOME
working_dir=$PWD
if test -d "${working_dir}/source"; then
  CAMPARI_HOME="${working_dir}"
else
  CAMPARI_HOME="$($DIRNAME "${working_dir}")"
fi
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
if test -f "${CAMPARI_HOME}/source/chainsaw.f90" -a -d "${CAMPARI_HOME}/bin/" -a -d "${CAMPARI_HOME}/lib/" -a -d "${CAMPARI_HOME}/params/"; then
  AC_MSG_NOTICE([CAMPARI_HOME is set to $CAMPARI_HOME])
  AC_MSG_CHECKING([whether $working_dir is source])
  to_go_back_dir=$working_dir
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

# This checks are used by AX_BLAS but need more TODO to comprehend them
AC_FC_DUMMY_MAIN
AC_F77_DUMMY_MAIN
# from blas
AC_CANONICAL_HOST # to check other variables
AC_MSG_CHECKING([which host cpu is present])
AC_MSG_RESULT([$host_cpu])
if test -n "${host_cpu}"; then
  ARCH="$host_cpu"
else
  AC_MSG_ERROR([Could not find cpu architecture using ac_canonical_host])
fi

# counting number of cpus available:
m4_include([m4/ax_count_cpus.m4])
AX_COUNT_CPUS()



dnl list of all possible optional components
all_needed_options="NetCDF4 LAPACK xdr FFTW3 fortran03"
all_flavour_options="specific_compiler debug_flags warning_flags old_xdr_lib multi-thread"


dnl When a component is found, append it here
needed_options=
flavour_options=

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
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: no])
],,[with_mpi=no])

if test ${with_mpi} != no; then
  SECTION_TITLE([MPI option compiler and flags])

  AX_MPI([AC_MSG_NOTICE([MPI support found. Following the possible variables and compilers
MPICC=$MPICC
MPICXX=$MPICXX
MPIF77=$MPIF77
MPIFC=$MPIFC
MPILIBS=$MPILIBS
FC=$FC])],
    [AC_MSG_WARN([initial search did not go well.])])
  MPILFLAGS=-O3
  MPIFF="$MPIFC"
  FC="$MPIFF"
  with_mpi=yes
  #
  # AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  #   use_mpi=no
  #   if test x"$with_mpi" = xyes; then
  #     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  #   else
  #     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  #   fi
  # ])
fi
AC_LANG_POP(Fortran)

# --------------------------------
#    specific compilers and flags
# --------------------------------
SECTION_TITLE([Compiler choice and specific flags settings])

# General FFLAGS from campari (I don't know how to
# look for them. They are cross compilers)

# preprocessor directives
PREPROCFLAGS="-cpp -DLINK_LAPACK -DLINK_NETCDF -DDISABLE_FLOAT -DLINK_FFTW -DLINK_XDR"

# autovectorization
AUTOVEC="-msse4.2"
# loop optimization
OPENLOOPS="-funroll-loops"
# misterious lflags
# LFLAGS=-O3


# enable specific compilation
AC_ARG_ENABLE(compiler,
AS_HELP_STRING([--enable-compiler=no/gnu/intel/cray/sun/pgi],
               [enable specific compilation, default: no]),
[case "${enableval}" in
        no)  compiler="no" ;;
        gnu) compiler="gnu" ;;
        cray) compiler="cray" ;;
        sun) compiler="sun" ;;
        intel) compiler="intel" ;;
        pgi) compiler="pgi" ;;
        *) compiler="${enableval}" ;; esac],
[compiler="no"])



AC_MSG_CHECKING([whether a specific compiler have been inserted])
AC_MSG_RESULT([$compiler])

if test "$compiler" != "no"; then
  echo "
  WARNING: if a specific compiler has been inserted manually, then the flags
  will be selected from standard values. To override specific flags please use
  the Makefile.local file in the source directory or the --use
  "
  flavour_options="${flavour_options} specific_compiler"

  # standard flags for different compilers
  INTELDEFAULTS="-O3 -real-size 64 -qoverride-limits -stand f03 -fpp"
  SUNDEFAULTS="-r8const -fpp -xO4 -ftrap=common,no%overflow"
  GNUDEFAULTS="-cpp -O3 -fall-intrinsics -ffpe-trap=invalid,zero -std=f2008 -fdefault-real-8 -fdefault-double-8"
  CRAYDEFAULTS="-K trap=inv,divz -emfF -dp -sreal64 -J '.'"
  PGIDEFAULTS="-Mallocatable=03 -Mfree -Mextend -Mpreprocess -Mbackslash -DPGI_FORTRAN"
  if test "$with_mpi" != "no"; then
    AC_MSG_WARN([
      Inserting a custom compiler with mpi support is not fully available
      (only the automatic flags will be inserted).
      Please consider adding manually flags and compiler
      in the Makefile.local])
    save_FC="$FC"
  fi
  case "${compiler}" in
          gnu) FC="gfortran"
          FCFLAGS=$GNUDEFAULTS ;;
          cray) FC="cray"
          FCFLAGS=$CRAYDEFAULTS;;
          sun) FC="sun"
          FCFLAGS=$SUNDEFAULTS ;;
          intel) FC="ifort"
          FCFLAGS=$INTELDEFAULTS;;
          pgi) FC="pgi"
          FCFLAGS=$PGIDEFAULTS;;
          *) FC="${compiler}"
          FCFLAGS="";;
  esac
  if test "$with_mpi" != "no"; then
    FC="${save_FC}"
  fi
  # checking directly the existance of compiler name
	if test -z "${SP_COMPILER}"; then
		AC_PATH_PROG(SP_COMPILER, $FC)
  fi
  AC_MSG_CHECKING([whether selected compiler is present in the system (ac_path_prog based)])
  if test -z "${SP_COMPILER}" ; then
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Unable to locate the compiler ${FC} using ac_path_prog])
  else
    AC_MSG_RESULT([yes])
  fi

  # AC_PROG_FC([$FC]) # why is this asking itself if there is a GNU Fortran compiler?
  AC_MSG_CHECKING([which new compiler is present])
  AC_MSG_RESULT([$FC])
  AC_MSG_CHECKING([which new compiler flags are present])
  AC_MSG_RESULT([$FCFLAGS])
fi

# overriding flag definition
AC_ARG_VAR([UFCFLAGS], [user defined compiler flags (it overrides FCFLAGS)])
if test "$UFCFLAGS" != ""; then
  FCFLAGS=$UFCFLAGS
fi

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
  WARNFLAGS="-Wintrinsics-std -Wsurprising -Warray-bounds -Waliasing -Wampersand -Wunderflow -Walign-commons -Wintrinsic-shadow -Wtabs -Wline-truncation -Wcharacter-truncation"
  # -Warray-temporaries  -Wimplicit-interface  -Wconversion
  flavour_options="${flavour_options} warning_flags"
else
  WARNFLAGS=
fi

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
DBFLAGS=
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
  AC_MSG_RESULT([$DBFLAGS])
fi

# threadflags
re='^[0-9]+$'
AC_MSG_CHECKING([how many CPUs are available])
AC_MSG_RESULT([${CPU_COUNT}])
n_threads=${CPU_COUNT}
m4_include([m4/ax_openmp.m4])
AC_ARG_ENABLE(multi-thread,
AS_HELP_STRING([--enable-multi-thread=yes/no],
               [enable multi-threading, default: no]),
[case "${enableval}" in
             '' | yes) multi_thread=yes ;;
             no)  multi_thread=no ;;
             *)
              if test $enableval -gt 0 ; then
                AC_MSG_NOTICE([selected $enableval threads])
                n_threads=$enableval
                multi_thread=yes
              else
                AC_MSG_ERROR([bad value ${enableval} for --enable-multi-thread])
              fi
             ;;
esac],
[multi_thread=no])
AC_MSG_CHECKING([whether compiler flags for multi-threading must be used])
AC_MSG_RESULT([$multi_thread])

THREADFLAGS=
# standard search of debug variables for compilationif test $debug != false; then
if test "$multi_thread" = "yes"; then
  AC_MSG_CHECKING([how many threads will be used])
  AC_MSG_RESULT([$n_threads])
  flavour_options="${flavour_options} multi-thread"
  AC_MSG_NOTICE([checking the OpenMP support with non-canonical ax openmp])
  AX_OPENMP(
    [AC_MSG_NOTICE([OpenMP support found.])],
    [AC_MSG_WARN([no OpenMP support has been found directly. If specific compiler has been used the flag will be set automatically])]
  )
  AC_MSG_NOTICE([(overloading variables) checking the OpenMP support with canonical ac openmp])
  AC_OPENMP()
  a=2
  AC_MSG_NOTICE([analysing the possible compiling languages for NON-posix support (only OpenMP)])
  if test -n "${OPENMP_FFLAGS}"; then
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
  if test "$compiler" != "no"; then
    case "${compiler}" in
              gnu) THREADFLAGS="-fopenmp" ;;
              # cray) DBFLAGS=$CRAYDEBUG;;
              # sun) DBFLAGS=$SUNDEBUG ;;
              intel) THREADFLAGS="-openmp";;
              # pgi) DBFLAGS=$PGIDEBUG;;
    esac
    AC_MSG_NOTICE([automatic threading flags set for $compiler compiler: $THREADFLAGS])
    a=1
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
    FFTW3_LIBS_TMP="-lfftw3_omp -lfftw3 -lm"
  else
    AC_MSG_NOTICE([found openmp support but not for the Fortran compiler. Threaded fftw3 flags set to openmp option.])
    FFTW3_LIBS_TMP="-lfftw3_omp -lfftw3 -lm"
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
      `$doit` >& multi_thread_assignments_conftest
      # output=`$doit`
      # AC_MSG_NOTICE([$output])
      AC_MSG_NOTICE([running script to detect non-functional behavior])
      if AC_TRY_EVAL(doit); then
        AC_MSG_NOTICE([evaluation successful for ${n_threads} threads])
        AC_MSG_NOTICE([please check manually the file multi_thread_assignments_conftest to see thread assignments])
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
      `$doit` >& multi_thread_assignments_conftest
      # output=`$doit`
      # AC_MSG_NOTICE([$output])
      AC_MSG_NOTICE([running script to detect non-functional behavior])
      if AC_TRY_EVAL(doit); then
        AC_MSG_NOTICE([evaluation successful for ${n_threads} threads])
        AC_MSG_NOTICE([please check manually the file multi_thread_assignments_conftest to see thread assignments])
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
  FLIBS=
  FCLIBS=
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
AC_ARG_ENABLE([f03-compiler-check],
              [AS_HELP_STRING([--disable-f03-compiler-check],
                              [disable check of ISO_C_BINDING support in Fortran compiler])],[],[enable_f03_compiler_check=yes])
test "x$enable_f03_compiler_check" = xno || enable_f03_compiler_check=yes
AC_MSG_RESULT($enable_f03_compiler_check)

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


# f03 ISO_C_BINDING check
# -----------------------
if test "x$enable_f03_compiler_check" = xyes; then
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
AC_MSG_CHECKING([whether F03 native code is desired])
AC_ARG_ENABLE([f03],
              [AS_HELP_STRING([--disable-f03],
                              [suppress Campari Fortran 2003 native code])],[],[enable_f03=yes])
test "x$enable_f03" = xno || enable_f03=yes
AC_MSG_RESULT([$enable_f03])

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
  AC_MSG_CHECKING([existance of file make_dependencies.sh in source])
  if test -e "${CAMPARI_HOME}/source/make_dependencies.sh"; then
    AC_MSG_RESULT([yes])
    AC_MSG_WARN([To manually update dependencies please run make_dependencies.sh . in the source directory.])
    AC_MSG_NOTICE([Running make_dependencies.sh . using default bash])
    $BASH make_dependencies.sh .
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([make_dependencies.sh file missing from source directory])
  fi
else
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([DEPENDENCIES file missing from source directory])
fi

# look for empty dir
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
	needed_options="$needed_options NetCDF4"
else
	AC_MSG_ERROR([Unable to find NetCDF4 (version >=4.4.0), we need NetCDF4 specifically for Fortran.
	Please be sure that command nf-config exists and it is not competing with previous versions (nc-config must have NO fortran bindings).])
fi

# test whether netcdf target compiler is the same of the selected one
AC_MSG_CHECKING([whether netcdf target compiler is the same compiler that has been selected])
if test "${NETCDF4_FC}" = "${FC}"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([
    ATTENTION!! The target compiler used for netcdf4 is different from the
    selected one for campari. Please consider reinstalling netcdf4 using
    the correct --exec-prefix. ])
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

# -------------------------------
#       Checking xdr libs
# -------------------------------
SECTION_TITLE([Checking for xdr libraries])

XDR_PREFIX=""
XDR_LIBS=""

# with-xdr ?
# -------------
AC_ARG_WITH(xdr,
	[AS_HELP_STRING([--with-xdr=<lib>], [use XDR library <lib>])],
  [])
case $with_xdr in
  yes | "") with_xdr=yes ;;
  *.a | *.so | *.so.* | *.o)
  if test -e $with_xdr; then
    XDR_PREFIX=$(dirname $with_xdr)
  else
    AC_MSG_WARN([library specified not present. It will be searched automatically])
  fi
  ;;
  */*)
  if test -d $with_xdr; then
    XDR_PREFIX=$with_xdr
  else
    AC_MSG_WARN([library specified not present. It will be searched automatically])
  fi
  ;;
  *) AC_MSG_ERROR([bad value $with_xdr for --with_xdr. It must be a location or a library file.]);;
esac
# XDR_LIBS can be "", "*.a" or "-l*"



# this is the case in which we want to use the old library build (libxdrf.a)
AC_MSG_CHECKING([search for old xdr library])
AC_ARG_WITH(xdr-old,
	[AS_HELP_STRING([--with-xdr-old=<lib>], [use XDR <lib> with the old format (libxdrf.a)])],
  [])
case $with_xdr_old in
  yes)
    AC_MSG_RESULT([yes])
    with_xdr_old=yes
  ;;
  *.a | *.so | *.so.* | *.o)
    AC_MSG_RESULT([yes])
    if test -e $with_xdr_old; then
      XDR_PREFIX=$(dirname $with_xdr_old)
    else
      AC_MSG_WARN([library specified not present. It will be searched automatically])
    fi
    with_xdr_old="yes"
  ;;
  */*)
    AC_MSG_RESULT([yes])
    if test -d $with_xdr_old; then
      XDR_PREFIX=$with_xdr_old
    else
      AC_MSG_WARN([library specified not present. It will be searched automatically])
    fi
    with_xdr_old="yes"
  ;;
  no | "")
    AC_MSG_RESULT([no])
    with_xdr_old="no"
  ;;
  *)
    AC_MSG_RESULT([not valid insertion])
    AC_MSG_ERROR([bad value $with-xdr-old for --with-xdr-old. It must be a location or a library file.])
  ;;
esac
a=2

#
#
#
#
#
#
#
# ATTENTION NEW VERSION OF XDR NOT IMPLEMENTED
if test "$with_xdr_old" = "no"; then
  # ATTENTION libxdrfile.a not implemented
  AC_MSG_WARN([

    !!WARNING!!
    Support for new versions of the xdr library is not available. With new
    versions the main functions have been rewritten with reference to c_pointers
    which implementation in fortran is not trivial. Indicatively the library
    needed is called libxdrf.a instead of libxdrfile.a Please send an email
    to the maintainers to have further information. We will resolve this
    incompatibility as soon as possible.
    ])
  with_xdr_old="yes"
fi
#
#
#
#
#
#
#

# look for xdr libraries automatically
# -----------------------------
# AC_MSG_CHECKING([if xdr libraries have been found automatically])
# looking for specific files
if test "$with_xdr_old" = "no"; then

  XDR_LIB_TOTEST="libxdrfile.a"
  #if prefix is set by user then search for lib.a under that path
  if test -n "${XDR_PREFIX}"; then
		for XDR_LIB_PATH in lib lib64
		do
		XDR_LIB_PATH_ABS=${XDR_PREFIX}/${XDR_LIB_PATH}
		AC_MSG_NOTICE([search for xdr library in: ${XDR_LIB_PATH_ABS}])
		if test -e ${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}; then
			XDR_LIBS="-L${XDR_LIB_PATH_ABS}"
			XDR_LIBS="${XDR_LIBS} -lxdrfile"
			XDR_CFLAGS="-I${XDR_PREFIX}/include/xdrfile/"
			AC_MSG_NOTICE([found XDR library in: ${XDR_LIB_PATH_ABS}])
			break;
		fi
		done
  else
  	# if prefix not provided, try pkg-config
  	AC_MSG_NOTICE([No directory was specified for --with-xdr. Trying to find xdr using pkg-config.])
    if ! test -z "${PKG_CONFIG}" ; then
      XDR_CONFIG="${PKG_CONFIG} xdr"
      XDR_CFLAGS="`${XDR_CONFIG} --cflags`"
      XDR_LIBS="`${XDR_CONFIG} --libs`"
    fi
  	if test -z "$XDR_LIBS"; then
    	AC_MSG_NOTICE([pkg-config was not able to find the xdr library.])
    	AC_MSG_NOTICE([Trying some common locations.])
  		#still not detected, then search some default standard location
  		for XDR_PREFIX in /usr/local /usr #skip searching/usr since it is standard system path
  		do
  			for XDR_LIB_PATH in lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
  			do
  			XDR_LIB_PATH_ABS=${XDR_PREFIX}/${XDR_LIB_PATH}
  			AC_MSG_NOTICE([search for xdr library in: ${XDR_LIB_PATH_ABS}])
  			if test -e ${XDR_LIB_PATH_ABS}/${XDR_LIB_TOTEST}; then
  				XDR_LIBS="-L${XDR_LIB_PATH_ABS}"
  				XDR_LIBS="${XDR_LIBS} -lxdrfile"
  				XDR_CFLAGS="-I${XDR_PREFIX}/include/xdrfile/"
  				AC_MSG_NOTICE([found xdr library in: ${XDR_LIB_PATH_ABS}])
  				break;
  			fi
  			done
  		if test -n "${XDR_LIBS}"; then
  			break
  		fi
  		done
    fi

    #fail to detect xdr_library, use the standard system path (which can be ommited)
    if test "${XDR_LIBS}" = ""; then
      AC_MSG_WARN([fail to detect xdr_library, trying with a macro in the
        standard system files (compilation not granted).])
      AC_SEARCH_LIBS(xdrfile_open, [xdrfile], [XDR_LIBS=$LIBS], [])
    fi

    # still anything found? installing..
    if test "${XDR_LIBS}" = ""; then
      AC_MSG_WARN([fail to detect xdr_library, looking for installation files])
      AC_MSG_CHECKING([if xdr libraries can are internally installed  in campariv3/ext_libs/installed ])
      a=0
      if test "$ax_ext_lib" = "yes"; then
        if test -d "${EXT_LIBS_DIR}/installed/xdr"; then
          AC_MSG_RESULT([yes])
          XDR_PREFIX="${EXT_LIBS_DIR}/installed/xdr"
        else
          AC_MSG_RESULT([no])
          AC_MSG_CHECKING([if xdr libraries can be installed in ../ext_libs/installed/xdr ])
          if test -e "${EXT_LIBS_DIR}/xdrfile*.tar.gz"; then #FAILING HERE!
            # TODO CHECK THE AUTOMATIC INSTALLATION
            # AC_MSG_RESULT([yes])
            XDR_PREFIX="${EXT_LIBS_DIR}/installed/xdr"
          else
            a=1
            AC_MSG_WARN([To install the xdrfile automatically the .tar.gz file is needed in campariv3/ext_libs])
          fi
        fi
      else
        a=1
      fi
      if test a = 1; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Impossible to find or install xdr libraries. Please check the
          installation and eventually try to add --with-xdr=PATH manually.])
        # THIS IS ONLY A WARNING BECAUSE WE NEED THE libxdrf.a ONLY FOR NOW.
        # MUST BE SET AS ERROR AFTERWARDS
      fi
    fi
  fi
fi
# looking for the old library
if test "$with_xdr_old" = "yes"; then
  flavour_options="$flavour_options old_xdr_lib"
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
    possible_dir="$($DIRNAME "${CAMPARI_HOME}")/xdr"
    # if prefix not provided, try pkg-config
    AC_MSG_NOTICE([No directory was specified for --with-xdr-old. Trying some common locations.])
    # #still not detected, then search some default standard location
    for XDR_PREFIX in $possible_dir /software/xdr /usr/local /usr #skip searching/usr since it is standard system path
    do
      for XDR_LIB_PATH in x86_64 lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
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
    AC_MSG_ERROR([Unable to find old XDR library. Please consider properly --with-xdr-old=<lib>])
  fi
fi


# NB: this check should be enforced enable-xdr-check=yes
AC_MSG_CHECKING([whether Fortran compiler is checked for XDR support])
AC_ARG_ENABLE([xdr-check],
  [AS_HELP_STRING([--disable-xdr-check],
    [(dev opt) disable check of xdr library support in Fortran compiler. ])],
    [],[enable_xdr_check=yes]) #ATTENTION TO modify TO YES
test "x$enable_xdr_check" = xno || enable_xdr_check=yes
AC_MSG_RESULT($enable_xdr_check)

# now checking xdr_ftest
# ----------------------
if test "x$enable_xdr_check" = xyes; then
  if test "$with_xdr_old" = "no"; then
    # checking the header file
    AC_MSG_CHECKING([whether it is possible to link xdrfile_open function])
    AC_TRY_LINK_FUNC(xdrfile_open,[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
    # checking the example tests (can be more than one, e.g. trr-support)
    AC_MSG_CHECKING([if Fortran compiler supports xdr libraries ])

    # NB: this examples are not working with more stringent compiler options
    # doit='$FC -c ${FCFLAGS} lib_test/xdr.F90'
    # AC_TRY_EVAL(doit)
    # doit='$FC -c ${FCFLAGS} lib_test/example/read-and-write-xtc/read-and-write-xtc.f90'
    # AC_TRY_EVAL(doit)
    # doit='$FC ${FCFLAGS} *.o  ${XDR_LIBS} ${LIBS}'
    LIBS_TMP="${LIBS}"
    CFLAGS_TMP="${CFLAGS}"
    AC_LANG_PUSH([C++])
      LIBS="$LIBS $XDRLIBS" # it is not this
      CFLAGS="$CFLAGS $XDR_CFLAGS"
      AC_CHECK_HEADER(xdrfile.h, [a=0], [a=1])
    AC_LANG_POP([C++])
    LIBS="${LIBS_TMP}"
    CFLAGS="${CFLAGS_TMP}"

    if test $a = 1; then
      AC_MSG_WARN([Could not find or include xdrfile.h])
    fi
    # checking compilation
    AC_MSG_NOTICE([trying compilation using:
${FC} -c ${FCFLAGS} ${DBFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}])
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
# Strangely this system is compiling but the xdrfopen extention should be invisible.
# when using libxdrfile.a. PLEASE CHANGE THE TEST ONCE CAMPARI HAS BEEN UPDATED
# FOR THE NEW VERSION OF XDR
    doit='${FC} -c ${FCFLAGS} ${DBFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}'
    if AC_TRY_EVAL(doit); then
      AC_MSG_RESULT([yes])
      needed_options="$needed_options xdr"
      AC_DEFINE([HAVE_LIBXDR], [1], [Defined if you have LIBXDR support])
    else
      AC_MSG_RESULT([no])
      AC_MSG_WARN([Unable to compile conftest.f90 for xdr library support.
        XDR checks unsuccesful. The library has not been linked properly.
        Please check it and eventually try to add --with-xdr=PATH manually.])
    fi

# now check for old version already installed
# --------------------------------------------
  elif test "$with_xdr_old" = "yes"; then
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
${FC} -c ${FCFLAGS} ${DBFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}])
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
    doit='${FC} -c ${FCFLAGS} ${DBFLAGS} ${FCFLAGS_f90} ${XDR_CFLAGS} conftest.f90 ${XDR_LIBS}'
    if AC_TRY_EVAL(doit); then
      AC_MSG_RESULT([yes])
     	needed_options="$needed_options xdr_old"
      AC_DEFINE([HAVE_LIBXDR], [1], [Defined if you have XDR support])
    else
      AC_MSG_RESULT([no])
      AC_MSG_WARN([Unable to compile conftest.f90 for xdr library support.
        XDR checks unsuccesful. The library has not been linked properly.
        Please check it and eventually try to add --with-xdr-old=PATH manually.])
    fi
  fi
else
 	needed_options="$needed_options xdr_old"
  AC_MSG_WARN([The abolishment of xdr checks are not safe and usually are restricted to dev.])
fi


# -------------------------------
#       Checking FFTW3 libs
# -------------------------------
SECTION_TITLE([Checking for FFTW3 libraries])

FFTW3_PREFIX=""
FFTW3_LIBS=""

# with-FFTW3 ?
# -------------
AC_ARG_WITH(FFTW3,
	[AS_HELP_STRING([--with-FFTW3=<lib>], [use FFTW3 library <lib>])],
  [])
case $with_FFTW3 in
  yes | "") with_FFTW3=yes ;;
  *.a | *.so | *.so.* | *.o)
  if test -e $with_FFTW3; then
    FFTW3_PREFIX=$(dirname $with_FFTW3)
  else
    AC_MSG_WARN([library specified not present. It will be searched automatically])
  fi
  ;;
  */*)
  if test -d $with_FFTW3; then
    FFTW3_PREFIX=$with_FFTW3
  else
    AC_MSG_WARN([library specified not present. It will be searched automatically])
  fi
  ;;
  *) AC_MSG_ERROR([bad value $with_FFTW3 for --with_FFTW3. It must be a location or a library file.]);;
esac
# FFTW3_LIBS can be "", "*.a" or "-l*"


a=2
# look for FFTW3 libraries automatically
# -----------------------------
# AC_MSG_CHECKING([if FFTW3 libraries have been found automatically])
# looking for specific files
FFTW3_LIB_TOTEST="libfftw3.a"
#if prefix is set by user then search for lib.a under that path
if test -n "${FFTW3_PREFIX}"; then
		for FFTW3_LIB_PATH in lib lib64
		do
		FFTW3_LIB_PATH_ABS=${FFTW3_PREFIX}/${FFTW3_LIB_PATH}
		AC_MSG_NOTICE([search for FFTW3 library in: ${FFTW3_LIB_PATH_ABS}])
		if test -e ${FFTW3_LIB_PATH_ABS}/${FFTW3_LIB_TOTEST}; then
			FFTW3_LIBS="-L${FFTW3_LIB_PATH_ABS}"
      if test "$multi_thread" = "yes"; then
        FFTW3_LIBS="${FFTW3_LIBS} ${FFTW3_LIBS_TMP}"
      else
  			FFTW3_LIBS="${FFTW3_LIBS} -lfft3"
      fi
			FFTW3_CFLAGS="-I${FFTW3_PREFIX}/include/"
			AC_MSG_NOTICE([found FFTW3 library in: ${FFTW3_LIB_PATH_ABS}])
			break;
		fi
		done
else
	# if prefix not provided, try pkg-config
	AC_MSG_NOTICE([No directory was specified for --with-FFTW3. Trying to find fftw3 using pkg-config.])
  if ! test -z "${PKG_CONFIG}" ; then
    FFTW3_CONFIG="${PKG_CONFIG} fftw3"
    FFTW3_CFLAGS="`${FFTW3_CONFIG} --cflags`" # to change to fortran flags
    FFTW3_LIBS="`${FFTW3_CONFIG} --libs`"
    # if test "$multi_thread" = "yes"; then
    AC_MSG_NOTICE([multi-threaded fftw3 pkg-config search not available])
    FFTW3_LIBS=""
    # fi
  fi
  if test "$FFTW3_LIBS" = ""; then
  	AC_MSG_NOTICE([pkg-config was not able to find the FFTW3 library.])
  	AC_MSG_NOTICE([Trying some common locations.])
		#still not detected, then search some default standard location
		for FFTW3_PREFIX in /usr /usr/local #skip searching /usr since it is standard system path
		do
			for FFTW3_LIB_PATH in lib lib64 lib/x86_64-linux-gnu lib64/x86_64-linux-gnu
			do
			FFTW3_LIB_PATH_ABS=${FFTW3_PREFIX}/${FFTW3_LIB_PATH}
			AC_MSG_NOTICE([search for FFTW3 library in: ${FFTW3_LIB_PATH_ABS}])
			if test -e ${FFTW3_LIB_PATH_ABS}/${FFTW3_LIB_TOTEST}; then
				FFTW3_LIBS="-L${FFTW3_LIB_PATH_ABS}"
        if test "$multi_thread" = "yes"; then
          FFTW3_LIBS="${FFTW3_LIBS} ${FFTW3_LIBS_TMP}"
        else
    			FFTW3_LIBS="${FFTW3_LIBS} -lfftw3"
        fi
				FFTW3_CFLAGS="-I${FFTW3_PREFIX}/include/"
				AC_MSG_NOTICE([found FFTW3 library in: ${FFTW3_LIB_PATH_ABS}])
				break;
			fi
			done
		if test -n "${FFTW3_LIBS}"; then
			break
		fi
		done
	fi
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
  AC_MSG_ERROR([
    A working fftw3 library must be provided.
    ])
fi


# NB: this check should be enforced enable-FFTW3-check=yes
# ------------------------------
AC_MSG_CHECKING([whether Fortran compiler is checked for FFTW3 support])
AC_ARG_ENABLE([FFTW3-check],
  [AS_HELP_STRING([--disable-FFTW3-check],
    [(dev opt) disable check of FFTW3 library support in Fortran compiler. ])],
    [],[enable_FFTW3_check=yes])
test "x$enable_FFTW3_check" = xno || enable_FFTW3_check=yes
AC_MSG_RESULT($enable_FFTW3_check)

# now checking FFTW3_ftest
# ----------------------
if test "x$enable_FFTW3_check" = xyes; then
  # checking the header file

  if test -n "${FFTW3_LIBS}"; then
    # checking multi_thread support for fftw3
    if test "${multi_thread}" = "yes"; then
      AC_MSG_NOTICE([trying compilation for multi-threaded fftw3 using:
${FC} ${FCFLAGS} ${DBFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} ${THREADFLAGS} conftest.f90 ${FFTW3_LIBS}])
      AC_MSG_CHECKING([whether NEEDED multi-thread support is available for fftw3])
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
      doit='${FC} ${FCFLAGS} ${DBFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} ${THREADFLAGS} conftest.f90 ${FFTW3_LIBS}'
      if AC_TRY_EVAL(doit); then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Unable to compile fftw3 with multi_thread])
      fi

    # standard NO multi-threaded fft3
    else
      AC_MSG_NOTICE([trying compilation for fftw3 using:
${FC} ${FCFLAGS} ${DBFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} ${THREADFLAGS} conftest.f90 ${FFTW3_LIBS}])
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
      doit='${FC} ${FCFLAGS} ${DBFLAGS} ${FFTW3_CFLAGS} ${FCFLAGS_f90} conftest.f90 ${FFTW3_LIBS}'
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
 	needed_options="$needed_options FFTW3"
  AC_MSG_WARN([The abolishment of FFTW3 checks are not safe and usually are restricted to dev.])
fi


# ----------------------------------
# ----------------------------------
# final testing of general compilation (if libs have diff target compiler)
# if test "$compiler" != "no"; then
#   AC_MSG_NOTICE([trying compilation using:
# ${FC} -c ${FCFLAGS} ${NETCDF4_FFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90 ${NETCDF4_FLIBS}])
#   AC_MSG_CHECKING([whether selected compiler and compiler flags can compile])
#   cat <<EOF >conftest.f90
# module conftest_module
# USE ISO_C_BINDING
# use netcdf
# end module conftest_module
# EOF
#   doit='${FC} -c ${FCFLAGS} ${NETCDF4_FFLAGS} ${DBFLAGS} ${FCFLAGS_f90} conftest.f90 ${NETCDF4_FLIBS}'
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

echo -e "\n" >& AS_MESSAGE_FD
SECTION_TITLE([Important variables])

# echo '        ' | sed -e's/./-/g' >& AS_MESSAGE_FD
echo "ARCH=$ARCH" >& AS_MESSAGE_FD
echo "FC=$FC" >& AS_MESSAGE_FD
echo "F77=$F77" >& AS_MESSAGE_FD
echo "FCFLAGS=$FCFLAGS" >& AS_MESSAGE_FD
if test "$with_mpi" != "no"; then
  echo "MPICC=$MPICC" >& AS_MESSAGE_FD
  echo "MPICXX=$MPICXX" >& AS_MESSAGE_FD
  echo "MPIF77=$MPIF77" >& AS_MESSAGE_FD
  echo "MPIFC=$MPIFC" >& AS_MESSAGE_FD
  echo "MPILIBS=$MPILIBS" >& AS_MESSAGE_FD
fi
# echo "CFFLAGS=$CFFLAGS" >& AS_MESSAGE_FD
# echo "CFLAGS=$CFLAGS" >& AS_MESSAGE_FD
echo "FFLAGS=$FFLAGS" >& AS_MESSAGE_FD
echo "PREPROCFLAGS=$PREPROCFLAGS" >& AS_MESSAGE_FD
echo "AUTOVEC=$AUTOVEC" >& AS_MESSAGE_FD
echo "OPENLOOPS=$OPENLOOPS" >& AS_MESSAGE_FD
if test "$multi_thread" = "yes"; then
  echo "OPENMP_CFLAGS=$OPENMP_CFLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_CXXFLAGS=$OPENMP_CXXFLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_F77FLAGS=$OPENMP_F77FLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_FFLAGS=$OPENMP_FFLAGS" >& AS_MESSAGE_FD
  echo "OPENMP_FCLAGS=$OPENMP_FCLAGS" >& AS_MESSAGE_FD
  echo "THREADFLAGS=$THREADFLAGS" >& AS_MESSAGE_FD
fi
echo "WARNFLAGS=$WARNFLAGS" >& AS_MESSAGE_FD
echo "LDFLAGS=$LDFLAGS" >& AS_MESSAGE_FD
echo "LFLAGS=$LFLAGS" >& AS_MESSAGE_FD
echo "FCFLAGS_f90=$FCFLAGS_f90" >& AS_MESSAGE_FD
echo "DBFLAGS=$DBFLAGS" >& AS_MESSAGE_FD
echo "XDR_CFLAGS=$XDR_CFLAGS" >& AS_MESSAGE_FD
echo "FFTW3_CFLAGS=$FFTW3_CFLAGS" >& AS_MESSAGE_FD
if test "$with_netcdf4_fortran" = "yes"; then
  echo '        ' | sed -e's/./-/g' >& AS_MESSAGE_FD
	echo "NETCDF4_VERSION=$NETCDF4_VERSION" >& AS_MESSAGE_FD
	echo "NETCDF4_FC=$NETCDF4_FC" >& AS_MESSAGE_FD
	echo "NETCDF4_FFLAGS=$NETCDF4_FFLAGS" >& AS_MESSAGE_FD
  echo "NETCDF4_FLIBS=$NETCDF4_FLIBS" >& AS_MESSAGE_FD
	# echo "" >& AS_MESSAGE_FD
  echo '        ' | sed -e's/./-/g' >& AS_MESSAGE_FD
fi
# echo "SHLIB_FFLAGS=$SHLIB_FFLAGS" >& AS_MESSAGE_FD
# echo "FPICFLAGS=$FPICFLAGS" >& AS_MESSAGE_FD
# echo "SAFE_FFLAGS=$SAFE_FFLAGS" >& AS_MESSAGE_FD
echo "LAPACK_LIBS=$LAPACK_LIBS" >& AS_MESSAGE_FD
echo "FFTW3_LIBS=$FFTW3_LIBS" >& AS_MESSAGE_FD
echo "BLAS_LIBS=$BLAS_LIBS" >& AS_MESSAGE_FD
echo "XDR_LIBS=$XDR_LIBS" >& AS_MESSAGE_FD
echo "LIBS=$LIBS" >& AS_MESSAGE_FD
echo "FLIBS=$FLIBS" >& AS_MESSAGE_FD
echo "FCLIBS=$FCLIBS" >& AS_MESSAGE_FD
echo "" >& AS_MESSAGE_FD


# ----------------------------------
#        Var substitutions
# ----------------------------------

# general settings
AC_SUBST(CAMPARI_HOME)
AC_SUBST(ARCH)

# general compiler flags
AC_SUBST(FC)
AC_SUBST(FCFLAGS)
AC_SUBST(DBFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LFLAGS)
AC_SUBST(FFFLAGS)
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(XDR_CFLAGS)
AC_SUBST(WARNFLAGS)
AC_SUBST(PREPROCFLAGS)
AC_SUBST(AUTOVEC)
AC_SUBST(OPENLOOPS)
AC_SUBST(THREADFLAGS)

# mpi-related
AC_SUBST(MPIFF)
AC_SUBST(MPILFLAGS)

# netcdf4
AC_SUBST(NETCDF4_VERSION)
AC_SUBST(NETCDF4_CC)
AC_SUBST(NETCDF4_CFLAGS)
AC_SUBST(NETCDF4_CPPFLAGS)
AC_SUBST(NETCDF4_LDFLAGS)
AC_SUBST(NETCDF4_LIBS)
AC_SUBST(NETCDF4_FC)
AC_SUBST(NETCDF4_FFLAGS)
AC_SUBST(NETCDF4_FLIBS)

# lapack and blas
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)

# general libs
AC_SUBST(LIBS)
AC_SUBST(FLIBS)
AC_SUBST(FCLIBS)
AC_SUBST(XDR_LIBS)
AC_SUBST(FFTW3_LIBS)


AC_MSG_NOTICE([generating makefiles])

dnl Process Makefile.in to create Makefile
AC_CONFIG_FILES([Makefile])

AC_MSG_NOTICE([setting back original directory])
if test "$to_go_back_dir" != "${CAMPARI_HOME}/source/"; then
  doit="$CD ${to_go_back_dir}"
  if AC_TRY_EVAL(doit); then
    AC_MSG_NOTICE([directory sucessfully changed to previous value: ${to_go_back_dir}.])
  fi
fi

AC_OUTPUT
