dnl UD macros for campariv3 configure (inspired by netcdf4)



dnl Convert a string to all uppercase.
dnl
define([uppercase],
[translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)])

dnl
dnl Check for an nm(1) utility.
dnl
AC_DEFUN([UD_PROG_NM],
[
    case "${NM-unset}" in
	unset) AC_CHECK_PROGS(NM, nm, nm) ;;
	*) AC_CHECK_PROGS(NM, $NM nm, nm) ;;
    esac
    AC_MSG_CHECKING(nm flags)
    case "${NMFLAGS-unset}" in
	unset) NMFLAGS= ;;
    esac
    AC_MSG_RESULT($NMFLAGS)
    AC_SUBST(NMFLAGS)
])


dnl Get information about Fortran data types.
dnl
AC_DEFUN([UD_FORTRAN_TYPES],
[
#    AC_REQUIRE([UD_PROG_FC])
    case "$FC" in
    '')
	;;
    *)
	AC_REQUIRE([UD_CHECK_FCALLSCSUB])
dnl	UD_CHECK_FORTRAN_TYPE(NF_INT1_T, byte integer*1 "integer(kind(1))")
dnl	UD_CHECK_FORTRAN_TYPE(NF_INT2_T, integer*2 "integer(kind(2))")
	UD_CHECK_FORTRAN_TYPE(NF_INT1_T, byte integer*1 "integer(kind=1)" "integer(selected_int_kind(2))")
	UD_CHECK_FORTRAN_TYPE(NF_INT2_T, integer*2 "integer(kind=2)" "integer(selected_int_kind(4))")

	case "${NF_INT1_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT1_T, "signed char", INT1)
		UD_CHECK_CTYPE_FORTRAN($NF_INT1_T, "short", INT1)
		UD_CHECK_CTYPE_FORTRAN($NF_INT1_T, "int", INT1)
		UD_CHECK_CTYPE_FORTRAN($NF_INT1_T, "long", INT1)
		;;
	esac
	case "${NF_INT2_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT2_T, short, INT2)
		UD_CHECK_CTYPE_FORTRAN($NF_INT2_T, int, INT2)
		UD_CHECK_CTYPE_FORTRAN($NF_INT2_T, long, INT2)
		;;
	esac
	UD_CHECK_CTYPE_FORTRAN(integer, int long, INT)
	UD_CHECK_CTYPE_FORTRAN(real, float double, REAL)
	UD_CHECK_CTYPE_FORTRAN(doubleprecision, double float, DOUBLEPRECISION)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCBYTE_T, byte integer*1 integer, byte)
	UD_CHECK_FORTRAN_NCTYPE(NCBYTE_T, byte integer*1 "integer(kind=1)" "integer(selected_int_kind(2))" integer, byte)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCSHORT_T, integer*2 integer, short)
	UD_CHECK_FORTRAN_NCTYPE(NCSHORT_T, integer*2 "integer(kind=2)" "integer(selected_int_kind(4))" integer, short)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_SHORT_T, $NCSHORT_T, short, SHRT_MIN, SHRT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCLONG_T, integer*4 integer, long)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_INT_T, integer, int, INT_MIN, INT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCFLOAT_T, real*4 real, float)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_FLOAT_T, $NCFLOAT_T, float, FLT_MIN, FLT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCDOUBLE_T, real*8 doubleprecision real, double)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_DOUBLE_T, $NCDOUBLE_T, double, DBL_MIN, DBL_MAX)
	;;
    esac
])

dnl Check for the name format of a Fortran-callable C routine.
dnl
dnl UD_CHECK_FCALLSCSUB
AC_DEFUN([UD_CHECK_FCALLSCSUB],
[
#    AC_REQUIRE([UD_PROG_FC])
    case "$FC" in
	'') ;;
	*)
	    AC_REQUIRE([UD_PROG_NM])
	    AC_BEFORE([UD_CHECK_FORTRAN_CTYPE])
	    AC_BEFORE([UD_CHECK_CTYPE_FORTRAN])
	    AC_MSG_CHECKING(for C-equivalent to Fortran routine \"SUB\")
	    cat >conftest.f <<\EOF
              call sub()
              end
EOF
	    doit='$FC -c ${FFLAGS} conftest.f'
	    if AC_TRY_EVAL(doit); then
		FCALLSCSUB=`$NM $NMFLAGS conftest.o | awk '
		    /SUB_/{print "SUB_";exit}
		    /SUB/ {print "SUB"; exit}
		    /sub_/{print "sub_";exit}
		    /sub/ {print "sub"; exit}'`
		case "$FCALLSCSUB" in
		    '') AC_MSG_ERROR(not found)
			;;
		    *)  AC_MSG_RESULT($FCALLSCSUB)
			;;
		esac
	    else
		AC_MSG_ERROR(Could not compile conftest.f)
	    fi
	    rm -f conftest*
	    ;;
    esac
])


dnl Check for a C type equivalent to a Fortran type.
dnl
dnl UD_CHECK_CTYPE_FORTRAN(ftype, ctypes, fmacro_root)
dnl
AC_DEFUN([UD_CHECK_CTYPE_FORTRAN],
[
    cat >conftestf.f <<EOF
           $1 values(4)
           data values /-1, -2, -3, -4/
           call sub(values)
           end
EOF
    for ctype in $2; do
	AC_MSG_CHECKING(if Fortran \"$1\" is C \"$ctype\")
	cat >conftest.c <<EOF
            #include <stdlib.h>
	    void $FCALLSCSUB(values)
		$ctype values[[4]];
	    {
		exit(values[[1]] != -2 || values[[2]] != -3);
	    }
EOF
	doit='$CC -c ${CPPFLAGS} ${CFLAGS} conftest.c'
	if AC_TRY_EVAL(doit); then
	    doit='$FC ${FFLAGS} -c conftestf.f'
	    if AC_TRY_EVAL(doit); then
	        doit='$FC -o conftest ${FFLAGS} ${LDFLAGS} conftestf.o conftest.o ${FLIBS} ${LIBS}'
	        if AC_TRY_EVAL(doit); then
		    doit=./conftest
		    if AC_TRY_EVAL(doit); then
		        AC_MSG_RESULT(yes)
		        cname=`echo $ctype | tr ' abcdefghijklmnopqrstuvwxyz' \
			    _ABCDEFGHIJKLMNOPQRSTUVWXYZ`
		        AC_DEFINE_UNQUOTED(NF_$3[]_IS_C_$cname, [1], [fortran to c conversion])
		        break
		    else
		        AC_MSG_RESULT(no)
		    fi
	        else
		    AC_MSG_ERROR(Could not link conftestf.o and conftest.o)
	        fi
	    else
		AC_MSG_ERROR(Could not compile conftestf.f)
	    fi
	else
	    AC_MSG_ERROR(Could not compile conftest.c)
	fi
    done
    rm -f conftest*
])


dnl Check for a C type equivalent to a Fortran type.
dnl
dnl UD_CHECK_CTYPE_FORTRAN(ftype, ctypes, fmacro_root)
dnl
AC_DEFUN([UD_CHECK_CTYPE_FORTRAN],
[
    cat >conftestf.f <<EOF
           $1 values(4)
           data values /-1, -2, -3, -4/
           call sub(values)
           end
EOF
    for ctype in $2; do
	AC_MSG_CHECKING(if Fortran \"$1\" is C \"$ctype\")
	cat >conftest.c <<EOF
            #include <stdlib.h>
	    void $FCALLSCSUB(values)
		$ctype values[[4]];
	    {
		exit(values[[1]] != -2 || values[[2]] != -3);
	    }
EOF
	doit='$CC -c ${CPPFLAGS} ${CFLAGS} conftest.c'
	if AC_TRY_EVAL(doit); then
	    doit='$FC ${FFLAGS} -c conftestf.f'
	    if AC_TRY_EVAL(doit); then
	        doit='$FC -o conftest ${FFLAGS} ${LDFLAGS} conftestf.o conftest.o ${FLIBS} ${LIBS}'
	        if AC_TRY_EVAL(doit); then
		    doit=./conftest
		    if AC_TRY_EVAL(doit); then
		        AC_MSG_RESULT(yes)
		        cname=`echo $ctype | tr ' abcdefghijklmnopqrstuvwxyz' \
			    _ABCDEFGHIJKLMNOPQRSTUVWXYZ`
		        AC_DEFINE_UNQUOTED(NF_$3[]_IS_C_$cname, [1], [fortran to c conversion])
		        break
		    else
		        AC_MSG_RESULT(no)
		    fi
	        else
		    AC_MSG_ERROR(Could not link conftestf.o and conftest.o)
	        fi
	    else
		AC_MSG_ERROR(Could not compile conftestf.f)
	    fi
	else
	    AC_MSG_ERROR(Could not compile conftest.c)
	fi
    done
    rm -f conftest*
])


dnl Check for a Fortran type equivalent to a netCDF type.
dnl
dnl UD_CHECK_FORTRAN_NCTYPE(forttype, possibs, nctype)
dnl
AC_DEFUN([UD_CHECK_FORTRAN_NCTYPE],
[
    AC_MSG_CHECKING(for Fortran-equivalent to netCDF \"$3\")
    for type in $2; do
	cat >conftest.f <<EOF
               $type foo
               end
EOF
	doit='$FC -c ${FFLAGS} conftest.f'
	if AC_TRY_EVAL(doit); then
	    break;
	fi
    done
    rm -f conftest.f conftest.o
    AC_DEFINE_UNQUOTED($1, $type, [type definition])
    AC_MSG_RESULT($type)
    $1=$type
])


dnl Check for a Fortran data type.
dnl
dnl UD_CHECK_FORTRAN_TYPE(varname, ftypes)
dnl
AC_DEFUN([UD_CHECK_FORTRAN_TYPE],
[
    for ftype in $2; do
	AC_MSG_CHECKING(for Fortran \"$ftype\")
	cat >conftest.f <<EOF
      subroutine sub(value)
      $ftype value
      end
EOF
	doit='$FC -c ${FFLAGS} conftest.f'
	if AC_TRY_EVAL(doit); then
	    AC_MSG_RESULT(yes)
	    $1=$ftype
	    AC_DEFINE_UNQUOTED($1, $ftype, [type thing])
	    break
	else
	    AC_MSG_RESULT(no)
	fi
    done
    rm -f conftest*
])
