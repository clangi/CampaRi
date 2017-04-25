#--------------------------------------------------------------------------#
# LICENSE INFO:                                                            #
#--------------------------------------------------------------------------#
#    This file is part of CAMPARI.                                         #
#                                                                          #
#    Version 2.0                                                           #
#                                                                          #
#    Copyright (C) 2014, The CAMPARI development team (current and former  #
#                        contributors)                                     #
#                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang #
#                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     #
#                        Nicholas Lyle, Nicolas Bloechliger                #
#                                                                          #
#    Website: http://sourceforge.net/projects/campari/                     #
#                                                                          #
#    CAMPARI is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation, either version 3 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    CAMPARI is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      #
#--------------------------------------------------------------------------#
# AUTHORSHIP INFO:                                                         #
#--------------------------------------------------------------------------#
#                                                                          #
# MAIN AUTHOR:   Andreas Vitalis                                           #
# CONTRIBUTIONS: Adam Steffen, Albert Mao                                  #
#                                                                          #
#--------------------------------------------------------------------------#

# standard
# CAMPARI_HOME=/software/campariv3
CAMPARI_HOME=/home/dgarolini/projects/CampaR/CampaRi/inst/campari_lib

# x86_64 locale (Default)
ARCH=x86_64
BIN_DIR=${CAMPARI_HOME}/bin
LIB_DIR=${CAMPARI_HOME}/lib
SRC_DIR=${CAMPARI_HOME}/source
FF= gfortran
# the following seems not used at all
# LLFLAGS=-O3
EXTRA_LIBS=  -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/opt/intel/compilers_and_libraries_2017.1.132/linux/ipp/lib/intel64 -L/opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2017.1.132/linux/tbb/lib/intel64/gcc4.7 -L/opt/intel/compilers_and_libraries_2017.1.132/linux/daal/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2017.1.132/linux/daal/../tbb/lib/intel64_lin/gcc4.4 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/opt/intel/compilers_and_libraries_2017.1.132/linux/ipp/lib/intel64 -L/opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2017.1.132/linux/tbb/lib/intel64/gcc4.7 -L/opt/intel/compilers_and_libraries_2017.1.132/linux/daal/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2017.1.132/linux/daal/../tbb/lib/intel64_lin/gcc4.4 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lquadmath -lcblas -lf77blas -latlas -llapack /software/xdr/x86_64/libxdrf.a -L/usr/lib/x86_64-linux-gnu -lfftw3 -lnetcdff 

MV=/bin/mv
RM=/bin/rm

COMPDEFAULTS=-g -O2 -g -O2    -I/usr/include/ 

# Gfortran GNU Fortran compiler
LFLAGS = 
FFLAGS = -cpp -DLINK_LAPACK -DLINK_NETCDF -DDISABLE_FLOAT -DLINK_FFTW -DLINK_XDR ${COMPDEFAULTS} -msse4.2 -funroll-loops
THREADFLAGS = 
# -msse4.2 is for autovectorization
# -funroll-loops Unroll loops whose number of iterations can be determined at compile time or upon entry to the loop.

MPIFF=
MPILFLAGS=
MPIEXTRA_LIBS=
MPIFFLAGS=
MPIEXTRA_LIBS=${EXTRA_LIBS}

# MPIFF=mpif90
# MPILFLAGS=-O3
# MPIEXTRA_LIBS=
# MPIFF = mpif90
# MPILFLAGS = -O3
# MPIFFLAGS = -DLINK_LAPACK -DLINK_NETCDF -DDISABLE_FLOAT -DLINK_XDR -DLINK_FFTW  ${WARNFLAGS} -msse4.2 -funroll-loops -I/usr/include
# MPIEXTRA_LIBS=/software/xdr/x86_64/libxdrf.a -lfftw3 -llapack -lnetcdf -lnetcdff




# Please define any local settings in the file Makefile.local and avoid certain strings in pathnames, as they
# are run through patsubst (see below)
include Makefile.local # for overwriting





# Include the script-generated module dependencies: note this uses LIB_DIR, SRC_DIR and ARCH
include DEPENDENCIES

# to add a module, add the module here, add the module in make_dependencies.sh, re-run make_dependencies, and compile
PREPREMODS=accept aminos atoms clusters commline contacts cutoffs diffrac dipolavg distrest dssps ems energies ewalds forces fos fyoc grandensembles grids inter interfaces ionize iounit keys math mcgrid mcsums mini molecule movesets mpistuff ncdm paircorr params pdb polyavg polypep sequen shakeetal system tabpot threads torsn ujglobals units wl zmatrix

PREMODS=$(addsuffix .f90,$(PREPREMODS))
MODULES=$(patsubst %.f90,$(SRC_DIR)/mod_%.f90,$(PREMODS))
MODOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/%.o,$(PREMODS))
MPIMODOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/mpi/%.o,$(PREMODS))
THREADMODOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/threads/%.o,$(PREMODS))
MPITHREADMODOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/mpi_threads/%.o,$(PREMODS))

print:
	${MODOBJS}

# to add a source, add the source here, add the source in make_dependencies.sh, re-run make_dependencies, and compile
PREPRESRCS=accsim allocate assignprm backup boundary cartld cartmd chainsaw clustering clustering_utils conrot constraint_solvers dssp emicroscopy energy energy_wrap ensemble ewald flow fmcscgrid fmsmcmpi force force_wrap fyczmat getkey graph_algorithms holes initial inner_loops inner_loops_en inner_loops_imp intmd intld makeio makepept math_utils mcmove mcstat minimize ncdm_proc nucconrot parsefiles parsekey particlefluctuation polar polymer proteus prtpdb readfyc readgrid readpdb readprm restart rigidmoves sanity_checks sav sav_auxil setconf sidechain string_utils structure summary thread_utils titrate topology torconrot torsion ujconrot ujsugar_pucker unbond utilities wanglandau

PRESRCS=$(addsuffix .f90,$(PREPRESRCS))
SOURCES=$(patsubst %.f90,${SRC_DIR}/%.f90,$(PRESRCS))
OBJECTS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/%.o,$(PRESRCS))
MPIOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/mpi/%.o,$(PRESRCS))
THREADOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/threads/%.o,$(PRESRCS))
MPITHREADOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/mpi_threads/%.o,$(PRESRCS))

# our strategy is to keep objects completely out of the SRC-directory (not even temporary) and to use a sub-directory for the MPI-version
# that means all the targets are only defined in that target directory, which is necessary in order to be able to address, maintain, and
# remove them independently.
# note that this includes the module interface definitions (.mod), which are unfortunately not compiler-independent. but this way it at least
# allows us to do preprocessing in modules without defining separate modules
# also note that dependencies are and should be exclusively handled through DEPENDENCIES
${OBJECTS} :
	${FF} ${FFLAGS} -c ${COMPDEFAULTS} -I${LIB_DIR}/${ARCH} $(patsubst ${LIB_DIR}/${ARCH}/%.o,${SRC_DIR}/%.f90,$@) -o $@

${MPIOBJS} :
	${MPIFF} ${MPIFFLAGS} -c ${COMPDEFAULTS} -DENABLE_MPI -I${LIB_DIR}/${ARCH}/mpi $(patsubst ${LIB_DIR}/${ARCH}/mpi/%.o,${SRC_DIR}/%.f90,$@) -o $@

${THREADOBJS} :
	${FF} ${FFLAGS} -c ${COMPDEFAULTS} ${THREADFLAGS} -DENABLE_THREADS -I${LIB_DIR}/${ARCH}/threads $(patsubst ${LIB_DIR}/${ARCH}/threads/%.o,${SRC_DIR}/%.f90,$@) -o $@

${MPITHREADOBJS} :
	${MPIFF} ${MPIFFLAGS} -c ${COMPDEFAULTS} ${THREADFLAGS} -DENABLE_MPI -DENABLE_THREADS -I${LIB_DIR}/${ARCH}/mpi_threads $(patsubst ${LIB_DIR}/${ARCH}/mpi_threads/%.o,${SRC_DIR}/%.f90,$@) -o $@


${MODOBJS} :
	${FF} ${FFLAGS} -c ${COMPDEFAULTS} -I${LIB_DIR}/${ARCH} $(patsubst ${LIB_DIR}/${ARCH}/%.o,${SRC_DIR}/mod_%.f90,$@) -o $@
	${MV} $(patsubst ${LIB_DIR}/${ARCH}/%.o,${SRC_DIR}/%.mod,$@) $(patsubst %.o,%.mod,$@)

${MPIMODOBJS} :
	${MPIFF} ${MPIFFLAGS} -c ${COMPDEFAULTS} -I${LIB_DIR}/${ARCH}/mpi -DENABLE_MPI $(patsubst ${LIB_DIR}/${ARCH}/mpi/%.o,${SRC_DIR}/mod_%.f90,$@) -o $@
	${MV} $(patsubst ${LIB_DIR}/${ARCH}/mpi/%.o,${SRC_DIR}/%.mod,$@) $(patsubst %.o,%.mod,$@)

${THREADMODOBJS} :
	${FF} ${FFLAGS} -c ${COMPDEFAULTS} ${THREADFLAGS} -I${LIB_DIR}/${ARCH}/threads -DENABLE_THREADS $(patsubst ${LIB_DIR}/${ARCH}/threads/%.o,${SRC_DIR}/mod_%.f90,$@) -o $@
	${MV} $(patsubst ${LIB_DIR}/${ARCH}/threads/%.o,${SRC_DIR}/%.mod,$@) $(patsubst %.o,%.mod,$@)

${MPITHREADMODOBJS} :
	${MPIFF} ${MPIFFLAGS} -c ${COMPDEFAULTS} ${THREADFLAGS} -I${LIB_DIR}/${ARCH}/mpi_threads -DENABLE_MPI -DENABLE_THREADS $(patsubst ${LIB_DIR}/${ARCH}/mpi_threads/%.o,${SRC_DIR}/mod_%.f90,$@) -o $@
	${MV} $(patsubst ${LIB_DIR}/${ARCH}/mpi_threads/%.o,${SRC_DIR}/%.mod,$@) $(patsubst %.o,%.mod,$@)



# library assembly is trivial (note that the name of the target doesn't match the target -> forced execution every time!)
library: $(OBJECTS)
	ar -rclvs ${LIB_DIR}/${ARCH}/lcampari.a ${MODOBJS} ${OBJECTS}

library_mpi: $(MPIOBJS)
	ar -rclvs ${LIB_DIR}/${ARCH}/mpi/lcampari_mpi.a ${MPIMODOBJS} ${MPIOBJS}

library_threads: $(THREADOBJS)
	ar -rclvs ${LIB_DIR}/${ARCH}/threads/lcampari_threads.a ${THREADMODOBJS} ${THREADOBJS}

library_mpi_threads: $(MPITHREADOBJS)
	ar -rclvs ${LIB_DIR}/${ARCH}/mpi_threads/lcampari_mpi_threads.a ${MPITHREADMODOBJS} ${MPITHREADOBJS}


# linking hopefully as well (see above note, true here due to target location!)
campari: library
	${FF} $(LFLAGS) -o ${BIN_DIR}/${ARCH}/campari ${LIB_DIR}/${ARCH}/chainsaw.o ${LIB_DIR}/${ARCH}/lcampari.a ${EXTRA_LIBS}

campari_mpi: library_mpi
	${MPIFF} $(MPILFLAGS) -DENABLE_MPI -o ${BIN_DIR}/${ARCH}/campari_mpi ${LIB_DIR}/${ARCH}/mpi/chainsaw.o ${LIB_DIR}/${ARCH}/mpi/lcampari_mpi.a ${MPIEXTRA_LIBS}

campari_threads: library_threads
	${FF} $(LFLAGS) ${THREADFLAGS} -DENABLE_THREADS -o ${BIN_DIR}/${ARCH}/campari_threads ${LIB_DIR}/${ARCH}/threads/chainsaw.o ${LIB_DIR}/${ARCH}/threads/lcampari_threads.a ${EXTRA_LIBS}

campari_mpi_threads: library_mpi_threads
	${MPIFF} $(MPILFLAGS) ${THREADFLAGS} -DENABLE_MPI -DENABLE_THREADS -o ${BIN_DIR}/${ARCH}/campari_mpi_threads ${LIB_DIR}/${ARCH}/mpi_threads/chainsaw.o ${LIB_DIR}/${ARCH}/mpi_threads/lcampari_mpi_threads.a ${EXTRA_LIBS}


# some fake targets, which are really clean-up commands
objclean:
	for i in ${LIB_DIR}/${ARCH}/*.o; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi/*.o; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/threads/*.o; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi_threads/*.o; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/*.mod; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi/*.mod; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/threads/*.mod; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi_threads/*.mod; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi_threads/lcampari_mpi_threads.a ${LIB_DIR}/${ARCH}/mpi/lcampari_mpi.a ${LIB_DIR}/${ARCH}/threads/lcampari_threads.a ${LIB_DIR}/${ARCH}/lcampari.a; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done

clean: objclean
	for i in ${BIN_DIR}/${ARCH}/campari*; do (if [ -e $$i ]; then ${RM} $$i; fi;) done
