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
# CONTRIBUTIONS: Davide Garolini                                           #
#                                                                          #
#--------------------------------------------------------------------------#

OBJEXT=o

libhsl_a_OBJECTS1 = eb13s.$(OBJEXT) eb13d.$(OBJEXT) ma48_sdeps.$(OBJEXT) eb13_sdeps.$(OBJEXT) eb13_ddeps.$(OBJEXT) ma48_ddeps.$(OBJEXT)
libhsl_a_OBJECTS2 = common90.$(OBJEXT) ma48_sdeps90.$(OBJEXT) ma48_ddeps90.$(OBJEXT) hsl_ma48s.$(OBJEXT) hsl_ma48d.$(OBJEXT) hsl_ma48s_ciface.$(OBJEXT) hsl_ma48d_ciface.$(OBJEXT)

nodist_modules_HEADERS = hsl_ma48_single.mod hsl_ma48_double.mod hsl_zd11_single.mod hsl_zd11_double.mod 

libhsl_ma48_a_SOURCES = \
	eb13_sdeps.f \
	eb13_ddeps.f \
	eb13s.f \
	eb13d.f \
	hsl_ma48s.f90 \
	hsl_ma48d.f90 \
	common90.f90 \
	ma48_sdeps.f \
	ma48_sdeps90.f90 \
	ma48_ddeps.f \
	ma48_ddeps90.f90 \
	hsl_ma48s_ciface.f90 \
	hsl_ma48d_ciface.f90


$(patsubst %,${LIB_DIR}/${ARCH}/%,$(libhsl_a_OBJECTS1)) : $(libhsl_ma48_a_SOURCES)
	$(FF) -c ${HSLFFLAGS} $(patsubst ${LIB_DIR}/${ARCH}/%.$(OBJEXT),${SRC_DIR}/hsl/%.f,$@) -o $@

$(patsubst %,${LIB_DIR}/${ARCH}/%,$(libhsl_a_OBJECTS2)) : $(libhsl_ma48_a_SOURCES)
	$(FF) -c ${HSLFFLAGS} $(patsubst ${LIB_DIR}/${ARCH}/%.$(OBJEXT),${SRC_DIR}/hsl/%.f90,$@) -o $@

${LIB_DIR}/${ARCH}/hsl_ma48_single.mod: hsl_ma48_single.mod
	$(MV) hsl_ma48_single.mod ${LIB_DIR}/${ARCH}/
hsl_ma48_single.mod:
	$(FF) -c ${HSLFFLAGS} hsl_ma48s.f90 -o ${LIB_DIR}/${ARCH}/hsl_ma48s.$(OBJEXT)
${LIB_DIR}/${ARCH}/hsl_ma48_double.mod: hsl_ma48_double.mod
	$(MV) hsl_ma48_double.mod ${LIB_DIR}/${ARCH}/
hsl_ma48_double.mod:
	$(FF) -c ${HSLFFLAGS} hsl_ma48d.f90 -o ${LIB_DIR}/${ARCH}/hsl_ma48d.$(OBJEXT)
${LIB_DIR}/${ARCH}/hsl_zd11_single.mod: hsl_zd11_single.mod
	$(MV) hsl_zd11_single.mod ${LIB_DIR}/${ARCH}/
hsl_zd11_single.mod:
	$(FF) -c ${HSLFFLAGS} ma48_sdeps90.f90 -o ${LIB_DIR}/${ARCH}/ma48_sdeps90.$(OBJEXT)
${LIB_DIR}/${ARCH}/hsl_zd11_double.mod: hsl_zd11_double.mod
	$(MV) hsl_zd11_double.mod ${LIB_DIR}/${ARCH}/
hsl_zd11_double.mod:
	$(FF) -c ${HSLFFLAGS} ma48_ddeps90.f90 -o ${LIB_DIR}/${ARCH}/ma48_ddeps90.$(OBJEXT)

${LIB_DIR}/${ARCH}/libhsl.a: $(patsubst %,${LIB_DIR}/${ARCH}/%,${libhsl_a_OBJECTS1}) $(patsubst %,${LIB_DIR}/${ARCH}/%,${libhsl_a_OBJECTS2}) $(patsubst %.mod,${LIB_DIR}/${ARCH}/%.mod,$(nodist_modules_HEADERS))
	ar -rvu ${LIB_DIR}/${ARCH}/libhsl.a  $(patsubst %,${LIB_DIR}/${ARCH}/%,${libhsl_a_OBJECTS1}) $(patsubst %,${LIB_DIR}/${ARCH}/%,${libhsl_a_OBJECTS2})

clean_hsl:
	for i in $(patsubst %,${LIB_DIR}/${ARCH}/%,$(libhsl_a_OBJECTS1)); do if [ -e $${i} ]; then ${RM} $${i}; fi; done
	for i in $(patsubst %,${LIB_DIR}/${ARCH}/%,$(libhsl_a_OBJECTS2)); do if [ -e $${i} ]; then ${RM} $${i}; fi; done
	for i in $(patsubst %,${LIB_DIR}/${ARCH}/%,$(nodist_modules_HEADERS)); do if [ -e $${i} ]; then ${RM} $${i}; fi; done
	for i in $(patsubst %,${LIB_DIR}/${ARCH}/mpi/%,$(nodist_modules_HEADERS)); do if [ -e $${i} ]; then ${RM} $${i}; fi; done
	for i in $(patsubst %,${LIB_DIR}/${ARCH}/threads/%,$(nodist_modules_HEADERS)); do if [ -e $${i} ]; then ${RM} $${i}; fi; done
	for i in $(patsubst %,${LIB_DIR}/${ARCH}/mpi_threads/%,$(nodist_modules_HEADERS)); do if [ -e $${i} ]; then ${RM} $${i}; fi; done
