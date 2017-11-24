#!/bin/bash
############################################################################
# LICENSE INFO:                                                            #
############################################################################
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
############################################################################
# AUTHORSHIP INFO:                                                         #
############################################################################
#                                                                          #
# MAIN AUTHOR:   Andreas Vitalis                                           #
#                                                                          #
############################################################################


if [ "$1" == "--help" ] || [ "$1" == "-h" ];
then
  echo -e "\nThis tool tries to extract a CAMPARI-compatible sequence file from (mostly) the header records of an input PDB file.\nRuntime-wise, it is extremely inefficient for the simplicity of the task.\nIt relies on sed/awk/grep, and will produce warnings for cases where additional steps are required.\nUsage:\n./gen_seqfile.sh PDB_IN SEQ_OUT XYZ WAT\n\nPDB_IN  : input PDB-file\nSEQ_OUT : output sequence file\nXYZ     : default choice for histidine residues (HID/HIE/HIP)\nWAT     : choice of water model (T3P/T4P/T4E/T5P/SPC); this is optional (if not provided, crystal waters are ignored)\n--------------------------------------------\n";
  exit 1
fi

if [[ "$1" == -* ]] || [[ "$2" == -* ]] || [[ "$3" == -* ]];
then
  echo -e "\nFatal. Unrecognized option. Use --help for usage information.\n";
  exit -1
fi

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ];
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
  echo -e "\nWarning. Output sequence file exists.\n";
  if [ -w "$2" ];
  then
    echo -e "And will be overwritten.\n";
  else
    echo -e "\nFatal. Cannot write to output sequence file. Use --help for usage information.\n";
    exit -1
  fi
else
  touch $2;
  zero=0;
  if [ $? -eq $zero ];
  then
    rm $2;
  else
    echo -e "\nFatal. Cannot write output sequence file. Use --help for usage information.\n";
    exit -1
  fi
fi

tmpfn=$2.lfdkgjljtl409580
grep -F 'SEQRES ' $1 | awk -v cl="?" -v vhis="$3" -v npc=0 -v ncp2=0 '{if ($1 != "SEQRES") {next}; if (cl != $3) {npc = 0}; for (k=5;k<=NF;k++) {suf1=""; if (npc == 0) {suf1="_N"}; suf2=""; npc = npc + 1;  npc2 = npc2 + 1; if (npc == $4) {suf2="_C"}; if ($k == "HIS") {print vhis  " " suf1 suf2 " " npc " " $3 " " npc2} else {print $k " " suf1 suf2 " " npc " " $3 " " npc2}}; cl=$3}' > ${tmpfn}

tmpfn2=$2.l09456lkjdgljk243
tmpfn3=$2.l0980d98sdfkjh2
tmpfn4=$2.kjhkhrtiuy4kjh3450

grep -F 'SEQRES ' $1 | awk -v cl="?" -v vhis="$3" -v npc=0 -v ncp2=0 '{if ($1 != "SEQRES") {next}; if (cl != $3) {print $3 "_" $4 "_" $5; cl=$3;}}' > ${tmpfn2}
for i in `cat $tmpfn2`; 
do
  chn=`echo ${i} | sed -e 's/_/ /g' | awk '{print $1}'`
  grep -v -F "ANISOU" $1 | grep  "CA  [A-Z][A-Z][A-Z] ${chn}" | sed -e "s/ ${chn}\([0-9]\)/  \1/" -e "s/ ${chn} /   /" | awk '{printf("%3s%07d_",$4,$5)}' > ${tmpfn3}
  chl=`echo ${i} | sed -e 's/_/ /g' | awk '{print $2}'`
  chshf=`grep -v -F "ANISOU" $1 | grep "CA  [A-Z][A-Z][A-Z] ${chn}" | head -n 1 | sed -e "s/ ${chn}\([0-9]\)/  \1/" -e "s/ ${chn} /   /" | awk '{print $5-1}'`
  first=`echo ${i} | sed -e 's/_/ /g' | awk '{print $3}'`
#  echo $chn $chl $first $chshf
  bla=`grep -F "REMARK 465     $first $chn" $1 | head -n1`
  if [ -n "${bla}" ]; then chshf_try=`echo "$bla" | awk '{print $NF-1}'`; if [[ $chsf_try -lt $chshf ]]; then chshf=$chshf_try; fi; fi
#  echo $chshf
  lo=`echo "${chshf} - 20" | bc`
  hi=`echo "${chshf} - 20" | bc`
  
  if [[ ${chl} -gt 7 ]]; then
    kk=`echo "${chl}-6" | bc`;
    for k in `seq 2 ${kk}`;
    do
      for ll in `seq $chshf $hi` `seq $chshf -1 $lo`;
        do
        patt=`awk -v a1=${chn} -v a2=${k} -v a3=0 -v a4=${ll} '{if ($3 == a1) {a3=a3+1; if (a3 >= a2) {printf("%3s%07d_",$1,$2+a4)}; if (a3 == a2+5) {exit}}}' $tmpfn`
        answ=`grep -F "$patt" ${tmpfn3}`
        if [ -n "${answ}" ];
        then
          chshf=${ll};
          break 2;
        fi
      done
    done
  fi
#  for k in `seq 1 ${k}-
  echo "${chn}_${chl}_${chshf}" >> ${tmpfn4}
done

awk '{if ($1 == "A") {nm="RPA"} else if ($1 == "G") {nm="RPG"} else if ($1 == "C") {nm="RPC"} else if ($1 == "U") {nm="RPU"} else if ($1 == "DG") {nm="DPG"} else if ($1 == "DA") {nm="DPA"} else if ($1 == "DC") {nm="DPC"} else if ($1 == "DT") {nm="DPT"} else if ($1 == "G") {nm="RPG"} else {nm=$1;} if (substr($2,1,1) == "_") {print nm $2 " " $3 " " $4 " " $5;} else {print nm " " $2 " " $3 " " $4;}}' ${tmpfn} > $2
mv $2 ${tmpfn}

sed -i -e 's/ACE_N /ACE /' -e 's/FOR_N /FOR /' -e 's/NME_C /NME /' -e 's/NH2_C /NH2 /' ${tmpfn}

for i in `grep "RP[A,T,G,C,U][ ,_]" $tmpfn | awk '{print $3}' | sort -u`; 
do 
  nops=`awk -v lt=${i} -v c1=0 -v c2=0 '{if (($1 == "ATOM") && ($5 == lt)) {if ($3 == "P") {c1=c1+1;}} print c1}' $1 | tail -n1`
  noO2s=`awk -v lt=${i} -v c1=0 -v c2=0 '{if (($1 == "ATOM") && ($5 == lt)) {if ((substr($3,1,2) == "O2") && (length($3) > 2)) {if ((substr($3,3,1) !~ /[[:alnum:]]/) && (substr($3,3,1) != " ")) {c2=c2+1;}}} print c2}' $1 | tail -n1`
  expe=`grep "RP[A,T,G,C,U][ ,_]" $tmpfn | awk -v ch=${i} '{if ($3 == ch) {print $3}}' | wc | awk '{print $1}'`
  expe2=`grep "RP[A,T,G,C,U][ ,_]" $tmpfn | awk -v ch=${i} '{if ($3 == ch) {print $3}}' | wc | awk '{print $1-1}'`
  if [[ $nops -eq $expe2 ]]; then awk -v ch=$i '{if (($3 == ch) && (substr($1,1,2) == "RP") && (substr($1,4,2) == "_N")) {sub("RP","RI",$1); sub("_N","",$1); print} else {print;}}' $tmpfn > $2; mv $2 $tmpfn ; fi 
  if [[ $noO2s -eq 0 ]]; then awk -v ch=$i '{if (($3 == ch) && (substr($1,1,2) == "RP")) {sub("RP","DP",$1); print} else {print;}}' $tmpfn > $2; mv $2 $tmpfn; fi
done

for i in `grep "DP[A,T,G,C,U][ ,_]" $tmpfn | awk '{print $3}' | sort -u`; 
do 
  nops=`awk -v lt=${i} -v c1=0 -v c2=0 '{if (($1 == "ATOM") && ($5 == lt)) {if ($3 == "P") {c1=c1+1;}} print c1}' $1 | tail -n1`
  expe=`grep "DP[A,T,G,C,U][ ,_]" $tmpfn | awk -v ch=${i} '{if ($3 == ch) {print $3}}' | wc | awk '{print $1}'`
  expe2=`grep "DP[A,T,G,C,U][ ,_]" $tmpfn | awk -v ch=${i} '{if ($3 == ch) {print $3}}' | wc | awk '{print $1-1}'`
  if [[ $nops -eq $expe2 ]]; then awk -v ch=$i '{if (($3 == ch) && (substr($1,1,2) == "DP") && (substr($1,4,2) == "_N")) {sub("DP","DI",$1); sub("_N","",$1); print} else {print;}}' $tmpfn > $2; mv $2 $tmpfn ; fi 
  if [[ $noO2s -gt 0 ]]; then awk -v ch=$i '{if (($3 == ch) && (substr($1,1,2) == "DP")) {sub("DP","RP",$1); print} else {print;}}' $tmpfn > $2; mv $2 $tmpfn; fi
done

for i in `grep -F 'SSBOND ' $1 | awk '{print $5 "_" $4 "_" $8 "_" $7}'`;
do 
  ix1=`echo "${i}" | sed -e 's/_/ /g' | awk '{print $1}'`;
  ix2=`echo "${i}" | sed -e 's/_/ /g' | awk '{print $3}'`;
  chn1=`echo "${i}" | sed -e 's/_/ /g' | awk '{print $2}'`;
  chn2=`echo "${i}" | sed -e 's/_/ /g' | awk '{print $4}'`;
  CSHF1=`sed -e 's/_/ /g' $tmpfn4 | awk -v c=$chn1 '{if ($1 == c) {print $3}}'`
  CSHF2=`sed -e 's/_/ /g' $tmpfn4 | awk -v c=$chn2 '{if ($1 == c) {print $3}}'`
 
  for nm in CYS_N_C CYS_N CYS_C CYS_D CYS;
  do 
    if [ -z "`grep -F "${nm} " ${tmpfn} | grep -F " $chn1"`" ]; then continue; fi
    for nm2 in CYS_N_C CYS_N CYS_C CYS_D CYS;
    do 
      if [ -z "`grep -F "${nm2} " ${tmpfn}`" ]; then continue; fi
      st1=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm} -v s=$CSHF1 '{print tx " " $1-s " " $2}'`
      if [ "${chn1}" == "${chn2}" ];
      then
        for offset2 in 0 1 -1 2 -2 3 -3 4 -4 5 -5;
        do
          offset=`echo "${offset2}+${CSHF1}" | bc`
          st1=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm} -v s=${offset} '{print tx " " $1-s " " $2}'`
          if [ -z "`grep -F "${st1}" ${tmpfn}`" ];
          then
            continue; #if [ "${nm}" == "CYS" ] && [ "${nm2}" == "CYS" ]; then echo -e "Error. SSBOND entry giving indeces and chains $ix1, $ix2, $chn1, $chn2 respectively, cannot be matched. Skipped."; fi
          fi
#         let us hope for the best
#          if [ "${nm}" == "CYS" ] && [ "${nm2}" == "CYS" ]; then echo -e "Warning. SSBOND entry giving indeces and chains $ix1, $ix2, $chn1, $chn2 respectively, could not be matched with DBREF entry. Assuming numbering from 1."; fi
          rept=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm} -v s=${offset} '{print tx " " $3-s " " $1 " " $2 " " $4}'`
          st2=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm2} -v s=${offset} '{print tx " " $3-s " " $4}'`
          if [ -z "`grep -F "${st2}" ${tmpfn}`" ];
          then
            continue; # if [ "${nm}" == "CYS" ] && [ "${nm2}" == "CYS" ]; then echo -e "Error. This does not work either. Skipped."; fi
          else
            rept="${rept} "`grep -F "${st2}" ${tmpfn} | awk '{print $4 " " $5}'`
            sed -i -e "s/${st1}/${rept}/" $tmpfn; break 3;
          fi
        done
      else
#       try with different offsets (first case should trigger usually)
        for offset2 in 0 1 -1 2 -2 3 -3 4 -4 5 -5;
        do
          for offset4 in 0 1 -1 2 -2 3 -3 4 -4 5 -5;
          do
            offset=`echo "${offset2}+${CSHF1}" | bc`
            offset3=`echo "${offset4}+${CSHF2}" | bc`
            st1=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm} -v s=${offset} '{print tx " " $1-s " " $2}'`
            st2=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm2} -v s=${offset3} '{print tx " " $3-s " " $4}'`
            if [ -z "`grep -F "${st1}" ${tmpfn}`" ];
            then
              break 1; # if [ "${nm}" == "CYS" ] && [ "${nm2}" == "CYS" ]; then echo -e "Error. SSBOND entry giving indeces and chains $ix1, $ix2, $chn1, $chn2 respectively, cannot be matched. Skipped."; fi
            fi
            if [ -z "`grep -F "${st2}" ${tmpfn}`" ];
            then
              continue; # if [ "${nm}" == "CYS" ] && [ "${nm2}" == "CYS" ]; then echo -e "Error. SSBOND entry giving indeces and chains $ix1, $ix2, $chn1, $chn2 respectively, cannot be matched. Skipped."; fi
            fi
#           if we're still here, we have a match (hope it's correct)
            rept=`echo "${i}" | sed -e 's/_/ /g' | awk -v tx=${nm} -v s=${offset} '{print tx " " $3-s " " $1 " " $2 " " $4}'`
            rept="${rept} "`grep -F "${st2}" ${tmpfn} | awk '{print $4 " " $5}'`
            sed -i -e "s/${st1}/${rept}/" $tmpfn; break 4;
          done
        done
      fi
    done
  done
done

nthere=`grep -F 'SSBOND ' $1 | awk -v k=0 '{if ($1 == "SSBOND") {k=k+1; print k}}' | tail -n1`

if [[ ${nthere} -gt 0 ]];
then
  nproc=`awk -vk=0 '{if (NF > 4) {k=k+1;} print k}' ${tmpfn} | tail -n1`
  echo -e "It seems that ${nproc}/${nthere} SSBOND records could be processed."
fi

#mv $tmpfn $2


for i in `grep -F 'HET  ' $1 | awk '{if ($1 == "HET") {print}}' | awk '{if (NF == 5) {print $2 "_" $3 "_" $4} else {print $2 "_" substr($3,1,1) "_" substr($3,2)}}'`;
do
  chn=`echo ${i} | sed -e 's/_/ /g' | awk '{print $2}'`
  rsn=`echo ${i} | sed -e 's/_/ /g' | awk '{print $1;}'`
  resn=`echo ${i} | sed -e 's/_/ /g' | awk '{if (length($1) == 3) {print $1} else if (length($1) == 2) {print $1 "X"} else if (length($1) == 1) {print $1 "XX"} else {print substr($1,1,3)}}'`
  resnl=`echo ${i} | sed -e 's/_/ /g' | awk '{print length($1);}'`

   
  inseq=`awk -v c1=$chn -v r1=$resn '{if ((substr($1,1,3) == r1) && ($3 == c1)) {print;}}' $tmpfn`
  if [[ -z "$inseq" ]]; then
    echo -e "Residue $rsn, which is described as ..."
    grep -F 'HETNAM ' $1 | grep ${rsn} | awk -v ss=${rsn} '{if (($2 == ss) || ($3 == ss)) {print}}'
    echo -e "... is considered unsupported.\nIf it is a residue supported natively by CAMPARI, it is much better to change names in both sequence and PDB files."
    inlnk=`grep -F "LINK     " $1 | awk -v rr=$rsn '{if (((NF == 12) || ((NF == 11) && (length($4 <= 4)))) && (($7 == rr) || ($3 == rr))) {print $3 "_" $7} else if (((NF == 10) || ((NF == 11) && (length($4 > 4)))) && (($6 == rr) || ($3 == rr))) {print $3 "_" $6}}'`
    if [[ -n "$inlnk" ]];
    then
      echo -e "Warning. Residue name $rsn appears to be implicated in chemical links.\nThis happens for example with glycoproteins and will almost certainly lead to errors down the road. Glycosylated residues, oligo-sugars, or similar entities must be setup manually." 
    fi
    if [[ "$resnl" -ne 3 ]];
    then
      echo -e "Warning. Because all additional residues except water are assumed unsupported and require 3-letter codes, the residue name $rsn has been adjusted to $resn.\nThis has to be done also in the PDB file for the information to be used."
    fi
    echo ${resn}_N_C >> ${tmpfn};
  else
    echo -e "Residue $rsn, which is described as ..."
    grep -F 'HETNAM ' $1 | grep ${rsn} | awk -v ss=${rsn} '{if (($2 == ss) || ($3 == ss)) {print}}'
    echo -e "... is considered unsupported and embedded in a recognized polymer chain.\nIf it is a residue supported natively by CAMPARI and names differ, it is much better to adopt names in both sequence and PDB files."
  fi
done

#waters
nwat=`grep -F "FORMUL " $1 | awk '{if ($3 == "HOH") {if (substr($4,1,1) == "*") {for (k=2;k<=10;k++) {if (substr($4,k,1) == "(") {print substr($4,2,k-2)}}} else {print "1";}}}'`
if [[ -n "$4" ]]; then watnam=`printf "%3.3s" $4`;
  if [[ -n "${nwat}" ]];
  then
  
    if [[ "$watnam" != "T3P" ]] && [[ "$watnam" != "SPC" ]] && [[ "$watnam" != "T4P" ]] && [[ "$watnam" != "T4E" ]] && [[ "$watnam" != "T5P" ]];
    then
      echo "Warning. The selected water residue name is not a natively supported residue type in CAMPARI."
    fi
    if [[ "$watnam" != "T3P" ]] && [[ "$watnam" != "SPC" ]];
    then
      echo "Warning. Water residue names in PDB file will have to be manually adjusted from HOH to $watnam (automatic adjustment supported only for T3P/SPC in conjunction with FMCSC_PDB_R_CONV)."
    fi
    for i in `seq 1 $nwat`; do echo $watnam >> ${tmpfn}; done;
  fi  
fi

echo "END"  >> ${tmpfn}

grep -F "REMARK 465" $1 > $tmpfn3
rm $tmpfn2
for i in `awk '{if (NF > 4) {print $1 "*" $3 "*" $4 "*" $6} else {print $1 "*" $2 "*" $3;}}' $tmpfn`;
do
  chn=`echo ${i} | sed -e 's/\*/ /g' | awk '{print $3}'`
  rsi=`echo ${i} | sed -e 's/\*/ /g' | awk '{print $2}'`
  rsn=`echo "${i}" | sed -e 's/\*/ /g' | awk '{print substr($1,1,3)}'`
  rsnr=`echo "${i}" | sed -e 's/\*/ /g' | awk '{print $1}'`
  if [[ "${rsn}" == "$3" ]]; then rsn=HIS; fi
  poss=`echo "${i}" | sed -e 's/\*/ /g' | awk '{print $4}'`
 
  if [[ -n "${chn}" ]]; then
    err=1;
    shf=`grep -F "${chn}_" $tmpfn4 | sed -e 's/_/ /g' | awk '{print $3}'`
    if [[ -n "`awk -v c1=$rsn -v c2=$chn -v c3=$rsi -v c4=$shf '{if (($3 == c1) && ($4 == c2) && ($5 == c3+c4)) {print}}' $tmpfn3`" ]];
    then
      echo "$rsnr MISSING $poss" >> $tmpfn2;       
    else
#      echo $rsn $rsi $chn $shf
      awk -v c1=$rsn -v c2=$chn -v c3=$rsi -v c4=$shf '{if (($3 == c1) && ($4 == c2) && ($5 == c3+c4)) {print}}' $tmpfn3
      echo "$rsnr $poss" >> $tmpfn2;
    fi
  else
    echo "$rsnr $poss" >> $tmpfn2;
#    echo ${i} $rsn $rsi $chn $shf
  fi
done

rm $tmpfn4
awk -v nnr=1 -v nr=0 -v isv=0 '{if ($2 == "MISSING") {if ((nr == 0) && (substr($1,4,2) != "_N")) {isv=1;} nr=nr+1; if (substr($1,4,2) == "_C") {if (isv == 1) {print nnr,nr,-NR}; isv=0; nr=0; nnr=nnr+1}} else if (nr > 0) {if (isv == 1) {print nnr,nr,NR}; isv=0;nnr=nnr+1;nr=0;}}' $tmpfn2 > $tmpfn4

if [[ -s $tmpfn4 ]];
then
  awk -v nr=0 -v hlp=$tmpfn4 -v isv=0 '{if ($2 == "MISSING") {if ((nr==0) && (substr($1,4,2) != "_N")) {isv=1; getline hlpl<hlp; split(hlpl,hlpi," ");} nr=nr+1; if (substr($1,4,2) == "_C") {isv=0; nr=0}; if ((hlpi[3] > 0) && (isv == 1) && (((nr == hlpi[2]-1) && (hlpi[2] > 1)) || ((nr == hlpi[2]) && (hlpi[2] == 1)))) {print "NH2"} else if ((hlpi[3] > 0) && (isv == 1) && (nr == hlpi[2])) {print "FOR";} else {print $1 " " $3}} else if (nr > 0) {if ((nr == 1) && (isv == 1)) {print "FOR"} else {print $1 " " $2}; isv=0;nr=0} else {print $1 " " $2}}' $tmpfn2 > $2

  awk -v nr=0 -v hlp=$tmpfn4 -v isv=0 '{if ($2 == "MISSING") {if ((nr==0) && (substr($1,4,2) != "_N")) {isv=1; getline hlpl<hlp; split(hlpl,hlpi," ");} nr=nr+1; if (substr($1,4,2) == "_C") {isv=0; nr=0;}; if ((hlpi[3] > 0) && (isv == 1) && (((nr == hlpi[2]-1) && (hlpi[2] > 1)) || ((nr == hlpi[2]) && (hlpi[2] == 1)))) {print "Warning. Replacing missing residue " $1 "-" NR " with NH2 (artificial chain break due to missing loop of " nr+1 " residues)."} else if ((hlpi[3] > 0) && (isv == 1) && (nr == hlpi[2])) {print "Warning. Replacing missing residue " $1 "-" NR " with FOR.";}} else if (nr > 0) {if ((nr == 1) && (isv == 1)) {print "Warning. Replacing " $1 "-" NR " with FOR. This residue will need to be deleted manually from the PDB file."}; isv=0;nr=0;}}' $tmpfn2
else
  awk '{if ($2 == "MISSING") {print $1 " " $3} else {print $1 " " $2}}' $tmpfn2 > $2
fi
ncys=`grep -c -F "CYS" $2`;
nsg=`awk -v k=0 '{if (($1 == "ATOM") && ($3 == "SG") && ($4 == "CYS")) {k=k+1; print k}; if (($1 == "ENDMDL") || ($1 == "END")) {exit;}}' $1 | tail -n1`

if [[ ${ncys} -gt ${nsg} ]];
then
  echo -e "There are ${ncys} cysteine residues yet only ${nsg} sulfur atoms (SG). This implies that disulfide bonds could theoretically be missing."
elif [[ ${ncys} -lt ${nsg} ]];
then
  echo -e "There are ${ncys} cysteine residues yet ${nsg} sulfur atoms (SG). This implies a misinterpretation of the input file."
fi

rm $tmpfn4 $tmpfn2 $tmpfn3 $tmpfn

