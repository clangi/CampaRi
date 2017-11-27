#!/bin/bash
#--------------------------------------------------------------------------#
# LICENSE INFO:                                                            #
#--------------------------------------------------------------------------#
#    This file is part of CAMPARI.                                         #
#                                                                          #
#    Version 3.0                                                           #
#                                                                          #
#    Copyright (C) 2017, The CAMPARI development team (current and former  #
#                        contributors)                                     #
#                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang #
#                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     #
#                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  #
#                        Davide Garolini, Jiri Vymetal                     #
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
# CONTRIBUTIONS: Adam Steffen, Albert Mao, Davide Garolini                 #
#                                                                          #
#--------------------------------------------------------------------------#

SRC_DIR=$1


if test "${SRC_DIR}" = "" -o -z "${SRC_DIR}"; then
  SRC_DIR="."
fi


MODULES="accept aminos atoms clusters commline contacts cutoffs diffrac dipolavg distrest dssps ems energies ewalds forces fos fyoc grandensembles grids inter interfaces ionize iounit keys math mcgrid mcsums mini molecule movesets mpistuff ncdm paircorr params pdb polyavg polypep sequen shakeetal system tabpot threads torsn ujglobals units wl zmatrix"

PRESOURCES="accsim allocate assignprm backup boundary cartld cartmd chainsaw clustering clustering_utils conrot constraint_solvers datasaw dssp emicroscopy energy energy_wrap ensemble ewald flow fmcscgrid fmsmcmpi force force_wrap fyczmat getkey graph_algorithms hacking holes initial inner_loops inner_loops_en inner_loops_imp intld intmd makeio makepept math_utils mcmove mcstat minimize ncdm_proc nucconrot parsefiles parsekey particlefluctuation polar polymer proteus prtpdb readfyc readgrid readpdb readprm restart rigidmoves sanity_checks sav sav_auxil setconf sidechain string_utils structure summary thread_utils titrate topology torconrot torsion ujconrot ujsugar_pucker unbond utilities wanglandau"

SOURCES=""
for i in $PRESOURCES;
do
  SOURCES=${SOURCES}"${i}.f90 "
done


MODOS=""
for i in $MODULES;
do
  MODOS=${MODOS}"mod_${i}.f90 "
done

echo "DEPENDENCIES: Found the following source files in the specified directory

${MODOS}"
echo ' '
echo 'DEPENDENCIES: loading modules relationships...'
rm ${SRC_DIR}/DEPENDENCIES
for i in $MODULES;
do
  echo \${LIB_DIR}/\${ARCH}/${i}.o \${LIB_DIR}/\${ARCH}/mpi/${i}.o \${LIB_DIR}/\${ARCH}/threads/${i}.o \${LIB_DIR}/\${ARCH}/mpi_threads/${i}.o \${LIB_DIR}/\${ARCH}/${i}.mod \${LIB_DIR}/\${ARCH}/mpi/${i}.mod \${LIB_DIR}/\${ARCH}/threads/${i}.mod \${LIB_DIR}/\${ARCH}/mpi_threads/${i}.mod: \${SRC_DIR}/mod_${i}.f90 >> ${SRC_DIR}/DEPENDENCIES
  DEPS=`grep -H " use ${i}" ${SOURCES} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS2=`grep -H " use ${i}" ${SOURCES} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/mpi/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS5=`grep -H " use ${i}" ${SOURCES} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/threads/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS7=`grep -H " use ${i}" ${SOURCES} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/mpi_threads/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS3=`grep -H " use ${i}" ${MODOS} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS4=`grep -H " use ${i}" ${MODOS} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/mpi/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS6=`grep -H " use ${i}" ${MODOS} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/threads/" $1}}' | sed -e "s/\.f90/.o/"`
  DEPS8=`grep -H " use ${i}" ${MODOS} | sed -e 's/:/ /' -e 's/,ON.*//' -e 's/, ON.*//' -e 's/:/ /' -e 's/,on.*//' -e 's/, on.*//' | sort -u | awk -v ii=${i} '{if (($2 == "use") && ($3 == ii)) {print "${LIB_DIR}/${ARCH}/mpi_threads/" $1}}' | sed -e "s/\.f90/.o/"`
  echo ${DEPS3} ${DEPS}: \${LIB_DIR}/\${ARCH}/${i}.o \${LIB_DIR}/\${ARCH}/${i}.mod >> ${SRC_DIR}/DEPENDENCIES
  echo ${DEPS4} ${DEPS2}: \${LIB_DIR}/\${ARCH}/mpi/${i}.o \${LIB_DIR}/\${ARCH}/mpi/${i}.mod >> ${SRC_DIR}/DEPENDENCIES
  echo ${DEPS6} ${DEPS5}: \${LIB_DIR}/\${ARCH}/threads/${i}.o \${LIB_DIR}/\${ARCH}/threads/${i}.mod >> ${SRC_DIR}/DEPENDENCIES
  echo ${DEPS8} ${DEPS7}: \${LIB_DIR}/\${ARCH}/mpi_threads/${i}.o \${LIB_DIR}/\${ARCH}/mpi_threads/${i}.mod >> ${SRC_DIR}/DEPENDENCIES
done

echo "DEPENDENCIES: Found the following source files in the specified directory

${SOURCES}"
echo ' '
echo 'DEPENDENCIES: loading sources relationships...'
for i in $SOURCES;
do
  echo \${LIB_DIR}/\${ARCH}/${i%.f90}.o \${LIB_DIR}/\${ARCH}/mpi/${i%.f90}.o \${LIB_DIR}/\${ARCH}/threads/${i%.f90}.o \${LIB_DIR}/\${ARCH}/mpi_threads/${i%.f90}.o: \${SRC_DIR}/${i} >> ${SRC_DIR}/DEPENDENCIES
done

echo 'DEPENDENCIES: loading macros.i relationships...'
MDEPS=`grep -H -l1 "#include \"macros.i\"" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/" $1}' | sed -e "s/\.f90/.o/"`
MDEPS2=`grep -H -l1 "#include \"macros.i\"" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/mpi/" $1}' | sed -e "s/\.f90/.o/"`
MDEPS3=`grep -H -l1 "#include \"macros.i\"" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/threads/" $1}' | sed -e "s/\.f90/.o/"`
MDEPS4=`grep -H -l1 "#include \"macros.i\"" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/mpi_threads/" $1}' | sed -e "s/\.f90/.o/"`
echo ${MDEPS} ${MDEPS2} ${MDEPS3} ${MDEPS4}: \${SRC_DIR}/macros.i >> ${SRC_DIR}/DEPENDENCIES

MDEPS=`grep -H -l1 "#include \"macros.i\"" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/" $1}' | sed -e "s/\.f90/.o/"`
MDEPS2=`grep -H -l1 "#include \"macros.i\"" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/mpi/" $1}' | sed -e "s/\.f90/.o/"`
MDEPS3=`grep -H -l1 "#include \"macros.i\"" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/threads/" $1}' | sed -e "s/\.f90/.o/"`
MDEPS4=`grep -H -l1 "#include \"macros.i\"" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/mpi_threads/" $1}' | sed -e "s/\.f90/.o/"`
echo ${MDEPS} ${MDEPS2} ${MDEPS3} ${MDEPS4}: \${SRC_DIR}/macros.i >> ${SRC_DIR}/DEPENDENCIES
