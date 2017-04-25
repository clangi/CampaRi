#!/bin/bash
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
#                                                                          #
#--------------------------------------------------------------------------#

if [ "$1" == "--help" ] || [ "$1" == "-h" ];
then
  echo -e "\nThis tool renumbers residues in PDB-files while trying to preserve all other information.\n\nUsage:\n./rersnmb.sh PDB_IN PDB_OUT RS_START\n\nPDB_IN  : input PDB-file\nPDB_OUT : output PDB-file\nRS_START: value for first residue in renumbered list\n\n--------------------------------------------\n";
  exit 1
fi

if [[ "$1" == -* ]] || [[ "$2" == -* ]] || [[ "$3" == -* ]];
then
  echo -e "\nFatal. Unrecognized option. Use --help for usage information.\n";
  exit -1
fi

if [ -z "$1" ] || [ -z "$2" ];
then
  echo -e "\nFatal. Too few arguments. Use --help for usage information.\n";
  exit -1
fi

if [ -r "$1" ] && [ -f "$1" ];
then
  dummy="dummy";# everything ok
else
  echo -e "\nFatal. Cannot open/read input PDB-file. Use --help for usage information.\n";
  exit -1
fi

if [ -r "$2" ] && [ -f "$2" ];
then
  echo -e "\nWarning. Output PDB-file exists.\n";
  if [ -w "$2" ];
  then
    echo -e "And will be overwritten.\n";
  else
    echo -e "\nFatal. Cannot write to output PDB-file. Use --help for usage information.\n";
    exit -1
  fi
else
  touch $2;
  zero=0;
  if [ $? -eq $zero ];
  then
    rm $2;
  else
    echo -e "\nFatal. Cannot write output PDB-file. Use --help for usage information.\n";
    exit -1
  fi
fi

declare -i frs
frs=$3;

cut -b 1-6 $1 > tmp1;
cut -b 7-22 $1 > tmp5;
cut -b 23-26 $1 > tmp2;
cut -b 27-100 $1 > tmp6;
paste -d " " tmp1 tmp2 > tmp3;
rm tmp1 tmp2;

#awk -v k=$3 -v curk=0 '{if (($1 == "ATOM") || ($1 == "HETATM")) {print;} }' tmp3 > tmp4;
#cut -b 24-27 tmp4 > tmp5;
#cut -b 1-23 tmp4 > tmp6;

i=`expr $frs - 1`
awk -v k=$i -v l=$i -v curk=0 '{if ($1 == "MODEL") {k=l}; if (($1 == "ATOM") || ($1 == "HETATM") || ($1 == "TER")) {if ($2 != curk) {k = k + 1; curk = $2;}; printf("%6-s %4d\n",$1,k);} else {print;}}' tmp3 > tmp4;
cut -b 1-6 tmp4 > tmp1;
cut -b 8-11 tmp4 > tmp2;

paste -d "" tmp1 tmp5 tmp2 tmp6 > $2;
rm tmp1 tmp2 tmp3 tmp4 tmp5 tmp6;


#cut -b 28-100 tmp4 > tmp9;
#paste -d "" tmp8 tmp9;

