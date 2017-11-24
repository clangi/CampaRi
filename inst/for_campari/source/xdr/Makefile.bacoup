# This make file is part of the xdrf package.
#
# (C) 1995 Frans van Hoesel, hoesel@chem.rug.nl
# Define ARCHXDR to be one the of architectures found in the conf directory.

# Set FORTRAN compiler and flags for ARCHXDR
FFLAGS  = $(XDRFFLAGS)
ifeq ($(system_os),mac)
FFLAGS  = "$(XDRFFLAGS) -I${SRC_DIR}/xdr/malloclib/"
endif

# Set C compiler and flags for ARCHXDR
CFLAGS 	= $(XDRCFLAGS)

# Set HASRANLIB to t if your system needs ranlib
HASRANLIB = f

# Set LIBS to libraries to link your c program with
# for SUN Solaris 2:
# LIBS = -lnsl
LIBS =

SHELL	= /bin/sh
# set RM to rm -f if you want to remove the intermediate C source file
# generated for the FORTRAN interface library
RMCMD 	= rm -f
#RMCMD	= /bin/true
M4 	= m4
M4FILE	= conf/$(ARCHXDR).m4

LFOBS 	= ${LIB_DIR}/${ARCH}/libxdrf.o
LOBS	= ${LIB_DIR}/${ARCH}/xdr3dfcoord.o

default: ${LIB_DIR}/${ARCH}/libxdrf.a ${LIB_DIR}/${ARCH}/ctest ${LIB_DIR}/${ARCH}/ftest

${LIB_DIR}/${ARCH}/libxdrf.a:  $(LFOBS) ${LIB_DIR}/${ARCH}/ftocstr.o
	ar cr ${LIB_DIR}/${ARCH}/libxdrf.a $?
	case x$(HASRANLIB) in xt ) echo ranlib; ranlib ${LIB_DIR}/${ARCH}/libxdrf.a ;; esac

${LIB_DIR}/${ARCH}/ctest:	ctest.c ${LIB_DIR}/${ARCH}/libxdrf.a
	$(CC) -o ${LIB_DIR}/${ARCH}/ctest $(CFLAGS) ctest.c ${LIB_DIR}/${ARCH}/libxdrf.a -lm $(LIBS)

${LIB_DIR}/${ARCH}/ftest:	ftest.f ${LIB_DIR}/${ARCH}/libxdrf.a
	$(FF) -o ${LIB_DIR}/${ARCH}/ftest $(FFLAGS) ftest.f ${LIB_DIR}/${ARCH}/libxdrf.a $(LIBS)
clean:
	rm -f $(LFOBS) $(LOBS) ${LIB_DIR}/${ARCH}/ftocstr.o ${LIB_DIR}/${ARCH}/libxdrf.a ${LIB_DIR}/${ARCH}/ftest ${LIB_DIR}/${ARCH}/ctest

tidy:
	rm -f $(LOBS) $(LFOBS)

${LIB_DIR}/${ARCH}/ftocstr.o: ftocstr.c
	$(CC) $(CFLAGS) -c ftocstr.c -o ${LIB_DIR}/${ARCH}/ftocstr.o

${LIB_DIR}/${ARCH}/libxdrf.o:	libxdrf.m4 $(M4FILE)
	$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c
	$(CC) $(CFLAGS) -c libxdrf.c -o ${LIB_DIR}/${ARCH}/libxdrf.o
	$(RMCMD) libxdrf.c

conf/.m4:
	@echo "ERROR: you didn't set ARCHXDR in the Makefile"
