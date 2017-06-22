# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_lib_netcdf4.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_NETCDF4([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of the NetCDF v4 library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   macro will call nc-config to check the output of the '--has-pnetcdf'
#   option and error out if the requested parallel isn't supported.
#
#   If the optional argument is omitted, no check is made to see if NetCDF
#   has parallel support.
#
#   The macro adds a --with-netcdf4 option accepting one of three values:
#
#     no   - do not check for the NetCDF4 library.
#     yes  - do check for NetCDF4 library in standard locations.
#     path - installation prefix for NetCDF version 4.
#
#   If NetCDF4 is successfully found, this macro calls
#
#     AC_SUBST(NETCDF4_VERSION)
#     AC_SUBST(NETCDF4_CC)
#     AC_SUBST(NETCDF4_CFLAGS)
#     AC_SUBST(NETCDF4_CPPFLAGS)
#     AC_SUBST(NETCDF4_LDFLAGS)
#     AC_SUBST(NETCDF4_LIBS)
#     AC_SUBST(NETCDF4_FC)
#     AC_SUBST(NETCDF4_FFLAGS)
#     AC_SUBST(NETCDF4_FLIBS)
#     AC_DEFINE(HAVE_NETCDF4)
#
#   It also sets
#
#     with_netcdf4="yes"
#     with_netcdf4_fortran="yes"    (if NetCDF has Fortran support)
#     with_netcdf4_parallel="yes"   (if NetCDF has MPI support)
#
#   If NetCDF4 is disabled or not found, this macros sets
#
#     with_netcdf4="no"
#     with_netcdf4_fortran="no"
#
#   Note it does not set with_netcdf4_parallel in this case.
#
#   Your configuration script can test $with_netcdf4 to take any further
#   actions. NETCDF4_{C,CPP,LD}FLAGS may be used when building with C or
#   C++. NETCDF4_F{FLAGS,LIBS} and NETCDF4_LDFLAGS should be used when
#   building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for NetCDF4 support
#        AX_LIB_NETCDF4()
#
#     2) dnl Check for serial NetCDF4 support
#        AX_LIB_NETCDF4([serial])
#
#     3) dnl Check for parallel NetCDF4 support
#        AX_LIB_NETCDF4([parallel])
#
#   One could test $with_netcdf4 for the outcome or display it as follows
#
#     echo "NetCDF v4 support:  $with_netcdf4"
#
#   One could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that NetCDF v4 was built with:
#
#     AX_LIB_NETCDF4([parallel])
#     if test "$with_netcdf4" = "yes"; then
#             CC="$NETCDF4_CC"
#     else
#             AC_MSG_ERROR([Unable to find NetCDF4, we need parallel NetCDF4.])
#     fi
#
# LICENSE
#
#   Copyright (c) 2016 Timothy Brown <tbrown@freeshell.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 2

AC_DEFUN([AX_LIB_NETCDF4], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
    netcdf4_requested_mode="serial"
elif test "m4_normalize(m4_default([$1],[]))" = "serial"  ; then
    netcdf4_requested_mode="serial"
elif test "m4_normalize(m4_default([$1],[]))" = "parallel"; then
    netcdf4_requested_mode="parallel"
else
    AC_MSG_ERROR([
Unrecognized value for AX[]_LIB_NETCDF4 within configure.ac.
If supplied, argument 1 must be either 'serial' or 'parallel'.
])
fi

dnl Add a default --with-netcdf4 configuration option.
AC_ARG_WITH([netcdf4],
  AS_HELP_STRING(
    [--with-netcdf4=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [serial],   [base directory of serial NetCDF4 installation],
            [parallel], [base directory of parallel NetCDF4 installation],
            [base directory of NetCDF4 installation])
  ),
  [if test "$withval" = "no"; then
     with_netcdf4="no"
   elif test "$withval" = "yes"; then
     with_netcdf4="yes"
   else
     with_netcdf4="yes"
     NETCDF4_PREFIX="${withval}"
     NC_CONFIG_TO_FIND="nc-config"
     #if prefix is set by user then search for nc-config
     if test -d "${NETCDF4_PREFIX}"; then
       for NC_CONFIG_PATH in bin /
       do
       NETCDF4_PATH_ABS=${NETCDF4_PREFIX}/${NC_CONFIG_PATH}
       AC_MSG_NOTICE([search for nc-config in: ${NETCDF4_PATH_ABS}])
       if test ! -h ${NETCDF4_PATH_ABS}; then
         if test -f ${NETCDF4_PATH_ABS}/${NC_CONFIG_TO_FIND} || \
         test -x ${NETCDF4_PATH_ABS}/${NC_CONFIG_TO_FIND}; then
           NC_CONFIG="${NETCDF4_PATH_ABS}/${NC_CONFIG_TO_FIND}"
           AC_MSG_NOTICE([found nc-config in: ${NETCDF4_PATH_ABS}])
           break;
         fi
       fi
       done
     else
       AC_MSG_ERROR([Inserted netcdf4 folder does not contain any nc-config to set flags.])
     fi
   fi],
   [with_netcdf4="yes"]
)

dnl Set defaults to blank
NETCDF4_CC=""
NETCDF4_VERSION=""
NETCDF4_CFLAGS=""
NETCDF4_CPPFLAGS=""
NETCDF4_LDFLAGS=""
NETCDF4_LIBS=""
NETCDF4_FC=""
NETCDF4_FFLAGS=""
NETCDF4_FLIBS=""

dnl Try and find NetCDF4 tools and options.
if test "$with_netcdf4" = "yes"; then
  if test -z "$NC_CONFIG"; then
    dnl Check to see if NC_CONFIG is in the path.
    AC_PATH_PROGS([NC_CONFIG], [nc-config], [])
    NETCDF4_PREFIX=$(AS_DIRNAME([$(AS_DIRNAME(["$NC_CONFIG"]))]))
  else
    AC_MSG_CHECKING([Using provided NetCDF4 prefix])
    AC_MSG_RESULT([$NC_CONFIG])
  fi

  AC_MSG_CHECKING([for NetCDF4 libraries])

  if test ! -f "$NC_CONFIG" || test ! -x "$NC_CONFIG"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([

Unable to locate NetCDF4 compilation helper script 'nc-config'.
Please specify --with-netcdf4=<LOCATION> as the full path prefix
where NetCDF4 has been installed.
NetCDF4 support is being disabled (equivalent to --with-netcdf4=no).
])
  with_netcdf4="no"
  with_netcdf4_fortran="no"
  else
    dnl Get the actual compiler used
    NETCDF4_CC=$(eval $NC_CONFIG --cc | $AWK '{print $[]1}')
    if test "$NETCDF4_CC" = "ccache"; then
      NETCDF4_CC=$(eval $NC_CONFIG --cc | $AWK '{print $[]2}')
    fi

    dnl Look for version
    NETCDF4_VERSION=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')

    dnl Look for the CFLAGS
    NETCDF4_CFLAGS=$(eval $NC_CONFIG --cflags)

    dnl Look for the LIBS and LDFLAGS
    NETCDF4_tmp_clibs=$(eval $NC_CONFIG --libs)

    dnl Sort out the tmp libs based on their prefixes
    for arg in $NETCDF4_tmp_clibs ; do
      case "$arg" in
        -L*) echo $NETCDF4_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
               || NETCDF4_LDFLAGS="$arg $NETCDF4_LDFLAGS"
          ;;
        -l*) echo $NETCDF4_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
               || NETCDF4_LIBS="$arg $NETCDF4_LIBS"
          ;;
      esac
    done
    # echo "libs=$NETCDF4_LIBS"
    # echo "cflags=$NETCDF4_CFLAGS"
    AC_MSG_RESULT([yes (version $[NETCDF4_VERSION])])

    dnl See if we need (and have) parallel support
    if test "$netcdf4_requested_mode" = "parallel" ; then
      with_netcdf4_parallel=$(eval $NC_CONFIG --has-pnetcdf)
      if test "$with_netcdf4_parallel" = "no" ; then
        AC_MSG_ERROR([
parallel NetCDF4 is not supported (while it was requested)
])
      fi
    fi



    AC_MSG_NOTICE([looking for nf-config (fortran explicit bindings)])
    NF_CONFIG_TO_FIND="nf-config"
    if test -d "${NETCDF4_PREFIX}"; then
      for NF_CONFIG_PATH in bin /
      do
      NETCDF4_PATH_ABS=${NETCDF4_PREFIX}/${NF_CONFIG_PATH}
      AC_MSG_NOTICE([search for nf-config in: ${NETCDF4_PATH_ABS}])
      if test ! -h ${NETCDF4_PATH_ABS}; then
        if test -f ${NETCDF4_PATH_ABS}/${NF_CONFIG_TO_FIND} || \
        test -x ${NETCDF4_PATH_ABS}/${NF_CONFIG_TO_FIND}; then
          NF_CONFIG="${NETCDF4_PATH_ABS}/${NF_CONFIG_TO_FIND}"
          AC_MSG_NOTICE([found nf-config in: ${NETCDF4_PATH_ABS}])
          break;
        fi
      fi
      done
    else
      AC_MSG_WARN([Inserted netcdf4 folder does not contain any nf-config to set fortran flags.])
    fi

    # checking fortran bindings
    AC_MSG_CHECKING([for matching NetCDF4 Fortran libraries using nf-config])
    if test ! -f "$NF_CONFIG" || test ! -x "$NF_CONFIG"; then
      AC_MSG_RESULT([no])
      if test -f "$NC_CONFIG" || test -x "$NC_CONFIG"; then
        AC_MSG_CHECKING([whether nc-config has still the Fortran bindings])
        NETCDF4_FC=$(eval $NC_CONFIG --fc | $AWK '{print $[]1}')
        if test -n "${NETCDF4_FC}" || test "${NETCDF4_FC}" != ""; then
          AC_MSG_RESULT([yes])
          NETCDF_FVERSION=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')
          AC_MSG_RESULT([yes (version $[NETCDF_FVERSION])])
          NETCDF4_FC=$(eval $NC_CONFIG --fc | $AWK '{print $[]1}')
          if test "$NETCDF4_FC" = "ccache"; then
            NETCDF4_FC=$(eval $NC_CONFIG --fc | $AWK '{print $[]2}')
          fi
          dnl Look for the FFLAGS
          NETCDF4_FFLAGS=$(eval $NC_CONFIG --fflags)

          dnl Look for the FLIBS and LDFLAGS
          NETCDF4_tmp_flibs=$(eval $NC_CONFIG --flibs)

          dnl Sort out the tmp libs based on their prefixes
          for arg in $NETCDF4_tmp_flibs ; do
            case "$arg" in
              -L*) echo $NETCDF4_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                     || NETCDF4_LDFLAGS="$arg $NETCDF4_LDFLAGS"
                ;;
              -l*) echo $NETCDF4_FLIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                     || NETCDF4_FLIBS="$arg $NETCDF4_FLIBS"
                ;;
            esac
          done
          with_netcdf4_fortran="yes"
        else
          AC_MSG_RESULT([no])
          AC_MSG_WARN([No Fortran bindings found (not even in nc-config).])
          with_netcdf4_fortran="no"
        fi
      else
        AC_MSG_WARN([No Fortran bindings found (neither nc-config exe present).])
        with_netcdf4_fortran="no"
      fi
      # echo "libs=$NETCDF4_LIBS"
      # echo "cfalgs=$NETCDF4_CFLAGS"
      # echo "glibs=$NETCDF4_FLIBS"
      # echo "fflags=$NETCDF4_FFLAGS"
    else
      NETCDF_FVERSION=$(eval $NF_CONFIG --version | $AWK '{print $[]2}')
      AC_MSG_RESULT([yes (version $[NETCDF_FVERSION])])
      NETCDF4_FC=$(eval $NF_CONFIG --fc | $AWK '{print $[]1}')
      if test "$NETCDF4_FC" = "ccache"; then
        NETCDF4_FC=$(eval $NF_CONFIG --fc | $AWK '{print $[]2}')
      fi
      dnl Look for the FFLAGS
      NETCDF4_FFLAGS=$(eval $NF_CONFIG --fflags)

      dnl Look for the FLIBS and LDFLAGS
      NETCDF4_tmp_flibs=$(eval $NF_CONFIG --flibs)

      dnl Sort out the tmp libs based on their prefixes
      for arg in $NETCDF4_tmp_flibs ; do
        case "$arg" in
          -L*) echo $NETCDF4_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                 || NETCDF4_LDFLAGS="$arg $NETCDF4_LDFLAGS"
            ;;
          -l*) echo $NETCDF4_FLIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                 || NETCDF4_FLIBS="$arg $NETCDF4_FLIBS"
            ;;
        esac
      done
      with_netcdf4_fortran="yes"
    fi

    # check if include dir is WRONG in nc-config
    AC_MSG_NOTICE([checking correct position of include directory])
    if test -f "$NF_CONFIG" || test -x "$NF_CONFIG"; then
      NX_CONFIG="$NF_CONFIG"
    elif test -f "$NC_CONFIG" || test -x "$NC_CONFIG"; then
      NX_CONFIG="$NC_CONFIG"
    else
      AC_MSG_WARN([no executable command present (i.e. no nc-config nor nf-config).])
    fi
    if test -n "${NX_CONFIG}"; then
      NETCDF4_INCLUDEDIR=$(eval $NX_CONFIG --includedir)
      if test ! -d "${NETCDF4_INCLUDEDIR}" || \
      test -z "$(ls -A ${NETCDF4_INCLUDEDIR})" || \
      test ! -e "${NETCDF4_INCLUDEDIR}/netcdf.h"; then
        AC_MSG_WARN([

The include directory is empty or it does
not exists or it does not contain netcdf.h.
])
        AC_MSG_NOTICE([looking for correct options using pkg-config.])
        if test -n "${PKG_CONFIG}" ; then
          NC_PKG_CONFIG=no
          module_name="netcdf"
          # if test $(${PKG_CONFIG} --libs "netcdf >= 4.4"); then
          #   module_name="netcdf-fortran"
          #   AC_MSG_NOTICE([found version higher than 4.4, the netcdf-fortran will be checked.])
          # fi
          for module_name in netcdf-fortran netcdf; do
            AC_MSG_CHECKING([${module_name}.pc in system locations])
            PKG_CHECK_EXISTS(${module_name}, [
              AC_MSG_RESULT([yes])
              NC_PKG_CONFIG=yes],
              [AC_MSG_RESULT([no])])
            if test "${NC_PKG_CONFIG}" != "no"; then
              PKG_CHECK_MODULES(NETCDF4, ${module_name},[
                AC_MSG_NOTICE([new libs and flags assigned to NETCDF4])
                NC_PKG_CONFIG=yes
                break;
                ],[
                AC_MSG_WARN([still nothing found for libs and flags])
                NC_PKG_CONFIG=no])
            fi
          done
          # still not found... manual setting of libs and flags
          if test "${NC_PKG_CONFIG}" = "no"; then
            # looking for netcdf.pc
            AC_MSG_NOTICE([looking directly for netcdf.pc in common locations.])
            for PKG_CONFIG_TO_FIND in netcdf-fortran.pc netcdf.pc
            do
            if test -d "${NETCDF4_PREFIX}"; then
              for PKG_CONFIG_PATH in / lib lib/pkgconfig; do
                NETCDF4_PATH_ABS=${NETCDF4_PREFIX}/${PKG_CONFIG_PATH}
                AC_MSG_NOTICE([search for $PKG_CONFIG_TO_FIND in: ${NETCDF4_PATH_ABS}])
                if test -f ${NETCDF4_PATH_ABS}/${PKG_CONFIG_TO_FIND}; then
                  PKG_CONFIG_NC="${NETCDF4_PATH_ABS}/${PKG_CONFIG_TO_FIND}"
                  AC_MSG_NOTICE([found $PKG_CONFIG_TO_FIND in: ${NETCDF4_PATH_ABS}])
                  break;
                fi
              done
            fi
            if test -n "${PKG_CONFIG_NC}"; then break; fi
            done
            # now let's use this netcdf*.pc
            if test -n "${PKG_CONFIG_NC}"; then
              NETCDF4_tmp_flibs=$(eval $PKG_CONFIG --libs "${PKG_CONFIG_NC}")
              NETCDF4_CFLAGS=$(eval $PKG_CONFIG --cflags "${PKG_CONFIG_NC}")
              includedir=$(eval $PKG_CONFIG --cflags-only-I "${PKG_CONFIG_NC}" | $SED "s/-I//1")
              for arg in $NETCDF4_tmp_flibs ; do
                case "$arg" in
                  -L*) echo $NETCDF4_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                         || NETCDF4_LDFLAGS="$arg"
                    ;;
                  -l*) echo $NETCDF4_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                         || NETCDF4_LIBS="$arg"
                    ;;
                esac
              done
              # echo "libs=$NETCDF4_LIBS"
              # echo "cfalgs=$NETCDF4_CFLAGS"
              # echo "glibs=$NETCDF4_FFLAGS"
              # echo "fflags=$NETCDF4_FFLAGS"
              # echo "includedir=$includedir"
              if test -n "$includedir" || \
              test -d "${includedir}" || \
              test -n "$(ls -A ${includedir})" || \
              test -e "${includedir}/netcdf.h"; then
                AC_MSG_NOTICE([new include directory found in $includedir])
                with_netcdf4="yes"
                AC_MSG_WARN([
the fortran flags (NETCDF4_FLIBS and NETCDF4_FFLAGS)
will be taken fron the C flags because for now pkg-config
is not set to look for them automatically.
])
                NETCDF4_FFLAGS="$NETCDF4_CFLAGS"
                NETCDF4_FLIBS="$NETCDF4_LIBS"
                with_netcdf4_fortran="yes"
              else
                with_netcdf4="no"
                with_netcdf4_fortran="no"
                AC_MSG_WARN([
pkg-config could not set a new include directory
automatically. This feature ended in a dead point
because include dir (which was wrongly set in nc-config)
was not found neither in netcdf*.pc.
Netcdf4 functionality will be turned off.
])
              fi
            else
              with_netcdf4="no"
              with_netcdf4_fortran="no"
              AC_MSG_WARN([
It was impossible to find netcdf*.pc. Not having the correct include,
bin and lib folder, the netcdf4 support will be disabled.
])
            fi
          fi
        # no pkg-config command found.
        else
          with_netcdf4="no"
          with_netcdf4_fortran="no"
          AC_MSG_WARN([
pkg-config not found. Not having the correct include,
bin and lib folder, the netcdf4 support will be disabled.
])
        fi
      else
        AC_MSG_NOTICE([includedir found correctly])
      fi
    fi

    dnl See if we can compile
    ax_lib_netcdf4_save_CC=$CC
    ax_lib_netcdf4_save_CPPFLAGS=$CPPFLAGS
    ax_lib_netcdf4_save_LIBS=$LIBS
    ax_lib_netcdf4_save_LDFLAGS=$LDFLAGS
    CC=$NETCDF4_CC
    CFLAGS=$NETCDF4_CFLAGS
    LIBS=$NETCDF4_LIBS
    LDFLAGS=$NETCDF4_LDFLAGS
    AC_CHECK_HEADER([netcdf.h], [ac_cv_netcdf4_h=yes], [ac_cv_netcdf4_h=no])
    AC_CHECK_LIB([netcdf], [nc_create], [ac_cv_libnetcdf4=yes],
                 [ac_cv_libnetcdf4=no])
    if test "$ac_cv_netcdf4_h" = "no" && \
       test "$ac_cv_libnetcdf4" = "no" ; then
        AC_MSG_WARN([Unable to compile NetCDF4 test program])
    fi

    CC=$ax_lib_netcdf4_save_CC
    CFLAGS=$ax_lib_netcdf4_save_CFLAGS
    LIBS=$ax_lib_netcdf4_save_LIBS
    LDFLAGS=$ax_lib_hdf5_save_LDFLAGS

    AC_SUBST([NETCDF4_VERSION])
    AC_SUBST([NETCDF4_CC])
    AC_SUBST([NETCDF4_CFLAGS])
    AC_SUBST([NETCDF4_LDFLAGS])
    AC_SUBST([NETCDF4_LIBS])
    AC_SUBST([NETCDF4_FC])
    AC_SUBST([NETCDF4_FFLAGS])
    AC_SUBST([NETCDF4_FLIBS])
    AC_DEFINE([HAVE_NETCDF4], [1], [Defined if you have NETCDF4 support])
  fi
fi
])
