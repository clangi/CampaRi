!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 2.0                                                           !
!                                                                          !
!    Copyright (C) 2014, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger                !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    CAMPARI is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Andreas Vitalis                                           !
! CONTRIBUTIONS: Adam Steffen                                              !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!
!-----------------------------------------------------------------------
!
! this routine corrects certain ambiguous names used for PDB atoms
! (but by no means all) to improve read-in success
! this list can always be appended
!
subroutine pdbcorrect(resname,nm,rs,rsshft,rsshftl)
!
  use pdb
  use system
  use sequen, ONLY: seqflag
!
  implicit none
!
  integer rs,rsshft,rsshftl(*)
  character(3) resname
  character(4) nm
  logical nucat
!
  nucat = .false.
!
! convert residue names to internal standard
!
  if (resname.eq.'NAC') resname = 'NME'
  if (resname.eq.'HSD') resname = 'HID'
  if (resname.eq.'HSE') resname = 'HIE'
  if (resname.eq.'HSP') resname = 'HIP'
  if (resname.eq.'SOD') resname = 'NA+'
  if (resname.eq.'CLA') resname = 'CL-'
  if (resname.eq.'POT') resname = 'K+ '
  if (resname.eq.'CES') resname = 'CS+'
  if (resname.eq.'NA ') resname = 'NA+'
  if (resname.eq.'CL ') resname = 'CL-'
  if (resname.eq.'K  ') resname = 'K+ '
  if (resname.eq.'CS ') resname = 'CS+'
  if (resname.eq.'ACT') resname = 'AC-'
  if (resname.eq.'HOH') then
    if (nm.eq." O  ") nm = " OW "
    if (pdb_convention(2).eq.2) then
      resname = 'SPC'
    else
      resname = 'T3P'
    end if
  end if
!
! nucleotide pdb-naming is a big mess
! use O2* to flag ribonucleotides (fatal if missing!) for straight nondescript pdb
  if (((nm.eq." O2'").OR.(nm.eq." O2*")).AND.&
 &    ((resname.eq.'  T').OR.(resname.eq.'  C').OR.(resname.eq.'  U').OR.(resname.eq.'  G').OR.(resname.eq.'  A'))) then
    resname(1:2) = 'R '
    nucat = .true.
  else if ((resname.eq.'DT5').OR.(resname.eq.'DG5').OR.(resname.eq.'DA5').OR.(resname.eq.'DU5').OR.(resname.eq.'DC5')) then
    resname(3:3) = resname(2:2)
    resname(1:2) = 'DI'
    nucat = .true.
  else if ((resname.eq.'RT5').OR.(resname.eq.'RG5').OR.(resname.eq.'RA5').OR.(resname.eq.'RU5').OR.(resname.eq.'RC5')) then
    resname(3:3) = resname(2:2)
    resname(1:2) = 'RI'
    nucat = .true.
! AMBER and GROMOS conversions  
  else if ((((resname.eq.' DT').OR.(resname.eq.'DT ').OR.(resname.eq.'DT3').OR.(resname.eq.'DTN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'DTH').AND.(pdb_convention(2).eq.2)).OR.(resname.eq.'  T')) then
    resname = 'DPT'
    nucat = .true.
  else if ((((resname.eq.' DG').OR.(resname.eq.'DG ').OR.(resname.eq.'DG3').OR.(resname.eq.'DGN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'DGU').AND.(pdb_convention(2).eq.2)).OR.(resname.eq.'  G')) then
    resname = 'DPG'
    nucat = .true.
  else if ((((resname.eq.' DA').OR.(resname.eq.'DA ').OR.(resname.eq.'DA3').OR.(resname.eq.'DAN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'DAD').AND.(pdb_convention(2).eq.2)).OR.(resname.eq.'  A')) then
    resname = 'DPA'
    nucat = .true.
  else if ((((resname.eq.' DU').OR.(resname.eq.'DU ').OR.(resname.eq.'DU3').OR.(resname.eq.'DUN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'DUR').AND.(pdb_convention(2).eq.2)).OR.(resname.eq.'  U')) then
    resname = 'DPU'
    nucat = .true.
  else if ((((resname.eq.' DC').OR.(resname.eq.'DC ').OR.(resname.eq.'DC3').OR.(resname.eq.'DCN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'DCY').AND.(pdb_convention(2).eq.2)).OR.(resname.eq.'  C')) then
    resname = 'DPC'
    nucat = .true.
  else if ((((resname.eq.' RT').OR.(resname.eq.'RT ').OR.(resname.eq.'RT3').OR.(resname.eq.'RTN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'THY').AND.(pdb_convention(2).eq.2))) then
    resname = 'RPT'
    nucat = .true.
  else if ((((resname.eq.' RG').OR.(resname.eq.'RG ').OR.(resname.eq.'RG3').OR.(resname.eq.'RGN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'GUA').AND.(pdb_convention(2).eq.2))) then
    resname = 'RPG'
    nucat = .true.
  else if ((((resname.eq.' RA').OR.(resname.eq.'RA ').OR.(resname.eq.'RA3').OR.(resname.eq.'RAN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'ADE').AND.(pdb_convention(2).eq.2))) then
    resname = 'RPA'
    nucat = .true.
  else if ((((resname.eq.' RU').OR.(resname.eq.'RU ').OR.(resname.eq.'RU3').OR.(resname.eq.'RUN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'URA').AND.(pdb_convention(2).eq.2))) then
    resname = 'RPU'
    nucat = .true.
  else if ((((resname.eq.' RC').OR.(resname.eq.'RC ').OR.(resname.eq.'RC3').OR.(resname.eq.'RCN'))&
 &.AND.(pdb_convention(2).eq.4)).OR.((resname.eq.'CYT').AND.(pdb_convention(2).eq.2))) then
    resname = 'RPC'
    nucat = .true.
! CHARMM conversion
  else if (pdb_convention(2).eq.3) then
    if ((resname.eq.'THY').AND.(nm.eq." O2'")) then
      resname = 'R T'
      nm = ' O2*'
      nucat = .true.
    else if ((resname.eq.'GUA').AND.(nm.eq." O2'")) then
      resname = 'R G'
      nm = ' O2*'
      nucat = .true.
    else if ((resname.eq.'ADE').AND.(nm.eq." O2'")) then
      resname = 'R A'
      nm = ' O2*'
      nucat = .true.
    else if ((resname.eq.'URA').AND.(nm.eq." O2'")) then
      resname = 'R U'
      nm = ' O2*'
      nucat = .true.
    else if ((resname.eq.'CYT').AND.(nm.eq." O2'")) then
      resname = 'R C'
      nm = ' O2*'
      nucat = .true.
    else if (resname.eq.'THY') then
      resname = 'DPT'
      nucat = .true.
    else if (resname.eq.'GUA') then
      resname = 'DPG'
      nucat = .true.
    else if (resname.eq.'ADE') then
      resname = 'DPA'
      nucat = .true.
    else if (resname.eq.'URA') then
      resname = 'DPU'
      nucat = .true.
    else if (resname.eq.'CYT') then
      resname = 'DPC'
      nucat = .true.
    end if
  else if ((resname.eq.'DPG').OR.(resname.eq.'DPA').OR.(resname.eq.'DPC').OR.(resname.eq.'DPT').OR.(resname.eq.'DPU').OR.&
 &         (resname.eq.'RPG').OR.(resname.eq.'RPA').OR.(resname.eq.'RPC').OR.(resname.eq.'RPT').OR.(resname.eq.'RPU').OR.&
 &         (resname.eq.'D5P').OR.(resname.eq.'R5P').OR.(resname.eq.'DIB').OR.(resname.eq.'RIB').OR.&
 &         (resname.eq.'DIG').OR.(resname.eq.'DIA').OR.(resname.eq.'DIC').OR.(resname.eq.'DIT').OR.(resname.eq.'DIU').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'RIA').OR.(resname.eq.'RIC').OR.(resname.eq.'RIT').OR.(resname.eq.'RIU')) then
    nucat = .true.
  end if
!
! polypeptides (note that HIS is a compatibility addition here, it is not an actual residue type)
  if ((resname.eq.'NME').OR.(resname.eq.'GLY').OR.&
 &    (resname.eq.'ALA').OR.(resname.eq.'AIB').OR.&
 &    (resname.eq.'ABA').OR.(resname.eq.'NVA').OR.&
 &    (resname.eq.'NLE').OR.(resname.eq.'DAB').OR.&
 &    (resname.eq.'ORN').OR.(resname.eq.'HYP').OR.&
 &    (resname.eq.'LEU').OR.(resname.eq.'PRO').OR.&
 &    (resname.eq.'VAL').OR.(resname.eq.'ILE').OR.&
 &    (resname.eq.'SER').OR.(resname.eq.'THR').OR.&
 &    (resname.eq.'MET').OR.(resname.eq.'CYS').OR.&
 &    (resname.eq.'TRP').OR.(resname.eq.'HIE').OR.&
 &    (resname.eq.'PHE').OR.(resname.eq.'HID').OR.&
 &    (resname.eq.'TYR').OR.(resname.eq.'HIP').OR.&
 &    (resname.eq.'ASP').OR.(resname.eq.'ASN').OR.&
 &    (resname.eq.'GLU').OR.(resname.eq.'GLN').OR.&
 &    (resname.eq.'ARG').OR.(resname.eq.'LYS').OR.&
 &    (resname.eq.'CYX').OR.(resname.eq.'LYD').OR.&
 &    (resname.eq.'SEP').OR.(resname.eq.'ASH').OR.&
 &    (resname.eq.'TPO').OR.(resname.eq.'GLH').OR.&
 &    (resname.eq.'PTR').OR.(resname.eq.'TYO').OR.&
 &    (resname.eq.'PCA').OR.(resname.eq.'KAC').OR.&
 &    (resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &    (resname.eq.'KM3').OR.(resname.eq.'HIS')) then
    if (nm.eq.' HN ')  nm = ' H  '
    if (nm.eq.' HN1')  nm = '1H  '
    if (nm.eq.' HN2')  nm = '2H  '
    if (nm.eq.' HN3')  nm = '3H  '
    if (nm.eq.' H1 ')  nm = '1H  '
    if (nm.eq.' H2 ')  nm = '2H  '
    if (nm.eq.' H3 ')  nm = '3H  '
    if (nm.eq.' HT1')  nm = '1H  '
    if (nm.eq.' HT2')  nm = '2H  '
    if (nm.eq.' HT3')  nm = '3H  '
    if (nm.eq.' OT1')  nm = '1OXT'
    if (nm.eq.'OCT1')  nm = '1OXT'
    if (nm.eq.' OC1')  nm = '1OXT'
    if ((resname.ne.'SEP').AND.(resname.ne.'PTR').AND.(resname.ne.'TYO').AND.(resname.ne.'TPO')) then
      if (nm.eq.' OT ')  nm = '1OXT'
    end if
    if (nm.eq.' O1 ')  nm = '1OXT'
    if (nm.eq.' O2 ')  nm = '2OXT'
    if (nm.eq.' OT2')  nm = '2OXT'
    if (nm.eq.'OCT2')  nm = '2OXT'
    if (nm.eq.' OC2')  nm = '2OXT'
    if (pdb_convention(2).eq.3) then
      if (nm.eq.' CAY') then
        nm = ' CH3'
        resname = 'ACE'
      end if
      if (nm.eq.' CY ') then
!       note we're using individual atoms for establishing that there is a cap patch (not perfect but oh well)
        nm = ' C  '
        resname = 'ACE'
        rsshft = rsshft + 1
        rsshftl(rsshft) = rs
      end if
      if (nm.eq.' OY ') then
        nm = ' O  '
        resname = 'ACE'
      end if
      if (nm.eq.' HY1') then
        nm = '1H  '
        resname = 'ACE'
      end if
      if (nm.eq.' HY2') then
        nm = '2H  '
        resname = 'ACE'
      end if
      if (nm.eq.' HY3') then
        nm = '3H  '
        resname = 'ACE'
      end if
      if (nm.eq.' NT ') then
        nm = ' N  '
        resname = 'NME'
        rsshft = rsshft + 1
        rsshftl(rsshft) = rs
      end if
      if (nm.eq.' CAT') then
        nm = ' CH3'
        resname = 'NME'
      end if
      if (nm.eq.' HNT') then
        nm = ' H  '
        resname = 'NME'
      end if
      if (nm.eq.'HAT1') then
        nm = '1H  '
        resname = 'NME'
      end if
      if (nm.eq.'HAT2') then
        nm = '2H  '
        resname = 'NME'
      end if
      if (nm.eq.'HAT3') then
        nm = '3H  '
        resname = 'NME'
      end if
      if (nm.eq.' NT2') then
        nm = ' N  '
        resname = 'NH2'
        rsshft = rsshft + 1
        rsshftl(rsshft) = rs
      end if
      if (nm.eq.'HT21') then
        nm = '1H  '
        resname = 'NH2'
      end if
      if (nm.eq.'HT22') then
        nm = '2H  '
        resname = 'NH2'
      end if
    end if
  end if
!
! now for specific residue type
! caps: ACE
  if (resname.eq.'ACE') then
    if (nm.eq.' CY ')  nm = ' C  '
    if (nm.eq.' CAY')  nm = ' CH3'
    if (nm.eq.' CA ')  nm = ' CH3'
    if (nm.eq.' OY ')  nm = ' O  '
    if (nm.eq.' HY1')  nm = '1H  '
    if (nm.eq.' HY2')  nm = '2H  '
    if (nm.eq.' HY3')  nm = '3H  '
    if (nm.eq.'HH31')  nm = '1H  '
    if (nm.eq.'HH32')  nm = '2H  '
    if (nm.eq.'HH33')  nm = '3H  '
    if (nm.eq.'1HH3')  nm = '1H  '
    if (nm.eq.'2HH3')  nm = '2H  '
    if (nm.eq.'3HH3')  nm = '3H  '
    if (nm.eq.' H1 ')  nm = '1H  '
    if (nm.eq.' H2 ')  nm = '2H  '
    if (nm.eq.' H3 ')  nm = '3H  '
!
! FOR
  else if (resname.eq.'FOR') then
    if (nm.eq.' CY ')  nm = ' C  '
    if (nm.eq.' OY ')  nm = ' O  '
    if (nm.eq.' HY ')  nm = ' H  '
!
! NME
  else if (resname.eq.'NME') then
    if (nm.eq.' NT ')  nm = ' N  '
    if (nm.eq.' CT ')  nm = ' CH3'
    if (nm.eq.' CAT')  nm = ' CH3'
    if (nm.eq.' CA ')  nm = ' CH3'
    if (nm.eq.' HNT')  nm = ' H  '
    if (nm.eq.' HT1')  nm = '1H  '
    if (nm.eq.' HT2')  nm = '2H  '
    if (nm.eq.' HT3')  nm = '3H  '
    if (nm.eq.'HAT1')  nm = '1H  '
    if (nm.eq.'HAT2')  nm = '2H  '
    if (nm.eq.'HAT3')  nm = '3H  '
    if (nm.eq.'HH31')  nm = '1H  '
    if (nm.eq.'HH32')  nm = '2H  '
    if (nm.eq.'HH33')  nm = '3H  '
    if (nm.eq.'1HH3')  nm = '1H  '
    if (nm.eq.'2HH3')  nm = '2H  '
    if (nm.eq.'3HH3')  nm = '3H  '
    if (nm.eq.' H1 ')  nm = '1H  '
    if (nm.eq.' H2 ')  nm = '2H  '
    if (nm.eq.' H3 ')  nm = '3H  '
!
! NH2
  else if (resname.eq.'NH2') then
    if (nm.eq.' NT ')  nm = ' N  '
    if (nm.eq.' HT1')  nm = '1H  '
    if (nm.eq.' HT2')  nm = '2H  '
    if (nm.eq.' H1 ')  nm = '1H  '
    if (nm.eq.' H2 ')  nm = '2H  '
!
! note that standard amino acids sometimes use 2 and 3 as indices for CH2 units
! fix by shifting "3" to "1"
! GLY
  else if (resname.eq.'GLY') then
    if (nm.eq.' HA1')  nm = '1HA '
    if (nm.eq.' HA2')  nm = '2HA '
    if (nm.eq.' HA3')  nm = '1HA '
    if (nm.eq.'3HA ')  nm = '1HA '
!
! ALA
  else if (resname.eq.'ALA') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '3HB '
!
! VAL
  else if (resname.eq.'VAL') then
    if (nm.eq.'HG11')  nm = '1HG1'
    if (nm.eq.'HG12')  nm = '2HG1'
    if (nm.eq.'HG13')  nm = '3HG1'
    if (nm.eq.'HG21')  nm = '1HG2'
    if (nm.eq.'HG22')  nm = '2HG2'
    if (nm.eq.'HG23')  nm = '3HG2'
!
! LEU
  else if (resname.eq.'LEU') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.'HD11')  nm = '1HD1'
    if (nm.eq.'HD12')  nm = '2HD1'
    if (nm.eq.'HD13')  nm = '3HD1'
    if (nm.eq.'HD21')  nm = '1HD2'
    if (nm.eq.'HD22')  nm = '2HD2'
    if (nm.eq.'HD23')  nm = '3HD2'
!
! ILE
  else if (resname.eq.'ILE') then
    if (nm.eq.' CD ')  nm = ' CD1'
    if (nm.eq.'HG11')  nm = '1HG1'
    if (nm.eq.'HG12')  nm = '2HG1'
    if (nm.eq.'HG13')  nm = '1HG1'
    if (nm.eq.'3HG1')  nm = '1HG1'
    if (nm.eq.'HG21')  nm = '1HG2'
    if (nm.eq.'HG22')  nm = '2HG2'
    if (nm.eq.'HG23')  nm = '3HG2'
    if (nm.eq.'HD11')  nm = '1HD1'
    if (nm.eq.' HD1')  nm = '1HD1'
    if (nm.eq.'HD12')  nm = '2HD1'
    if (nm.eq.' HD2')  nm = '2HD1'
    if (nm.eq.'HD13')  nm = '3HD1'
    if (nm.eq.' HD3')  nm = '3HD1'
!
! PHE
  else if (resname.eq.'PHE') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
!
! PRO
  else if (resname.eq.'PRO') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '1HD '
    if (nm.eq.'3HD ')  nm = '1HD '
!
! MET
  else if (resname.eq.'MET') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HE1')  nm = '1HE '
    if (nm.eq.' HE2')  nm = '2HE '
    if (nm.eq.' HE3')  nm = '3HE '
!
! AIB
  else if (resname.eq.'AIB') then
    if (nm.eq.'HB11')  nm = '1HB1'
    if (nm.eq.'HB12')  nm = '2HB1'
    if (nm.eq.'HB13')  nm = '3HB1'
    if (nm.eq.'HB21')  nm = '1HB2'
    if (nm.eq.'HB22')  nm = '2HB2'
    if (nm.eq.'HB23')  nm = '3HB2'
!
! SER
  else if (resname.eq.'SER') then
    if (nm.eq.' OG1')  nm = ' OG '
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = ' HG '
    if (nm.eq.' HOG')  nm = ' HG '
!
! THR
  else if (resname.eq.'THR') then
    if (nm.eq.' OG ')  nm = ' OG1'
    if (nm.eq.' CG ')  nm = ' CG2'
    if (nm.eq.' HOG')  nm = ' HG1'
    if (nm.eq.'HG21')  nm = '1HG2'
    if (nm.eq.'HG22')  nm = '2HG2'
    if (nm.eq.'HG23')  nm = '3HG2'
!
! CYS
  else if ((resname.eq.'CYS').OR.(resname.eq.'CYX')) then
    if (nm.eq.' SG1')  nm = ' SG '
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = ' HG '
    if (nm.eq.' HSG')  nm = ' HG '
!
! ASN
  else if (resname.eq.'ASN') then
    if (nm.eq.' OD ')  nm = ' OD1'
    if (nm.eq.' ND ')  nm = ' ND2'
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.'HD21')  nm = '1HD2'
    if (nm.eq.'HND1')  nm = '1HD2'
    if (nm.eq.'HD22')  nm = '2HD2'
    if (nm.eq.'HND2')  nm = '2HD2'
!
! GLN
  else if (resname.eq.'GLN') then
    if (nm.eq.' OE ')  nm = ' OE1'
    if (nm.eq.' NE ')  nm = ' NE2'
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.'HE21')  nm = '1HE2'
    if (nm.eq.'HNE1')  nm = '1HE2'
    if (nm.eq.'HE22')  nm = '2HE2'
    if (nm.eq.'HNE2')  nm = '2HE2'
!
! TYR/TYO
  else if ((resname.eq.'TYR').OR.(resname.eq.'TYO')) then
    if (nm.eq.' HOH')  nm = ' HH '
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if ((resname.eq.'TYO').AND.(nm.eq.' OH ')) nm = ' OZ '
!
! TRP
  else if (resname.eq.'TRP') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HNE')  nm = ' HE1'
!
! HID
  else if (resname.eq.'HID') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HD ')  nm = ' HD2'
    if (nm.eq.' HE ')  nm = ' HE1'
    if (nm.eq.' HND')  nm = ' HD1'
!
! HIE
  else if (resname.eq.'HIE') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HD ')  nm = ' HD2'
    if (nm.eq.' HE ')  nm = ' HE1'
    if (nm.eq.' HNE')  nm = ' HE2'
!
! ASP/ASH
  else if ((resname.eq.'ASP').OR.(resname.eq.'ASH')) then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (resname.eq.'ASH') then
      if (nm.eq.' HO ') nm = ' HD2'
      if (nm.eq.' OH ') nm = ' OD2'
    end if
!
! GLU/GLH
  else if ((resname.eq.'GLU').OR.(resname.eq.'GLH')) then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (resname.eq.'GLH') then
      if (nm.eq.' HO ') nm = ' HE2'
      if (nm.eq.' OH ') nm = ' OE2'
    end if
!
! LYS/LYD
  else if ((resname.eq.'LYS').OR.(resname.eq.'LYD').OR.(resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &         (resname.eq.'KM2').OR.(resname.eq.'KM3')) then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '1HD '
    if (nm.eq.'3HD ')  nm = '1HD '
    if (nm.eq.' HE1')  nm = '1HE '
    if (nm.eq.' HE2')  nm = '2HE '
    if (nm.eq.' HE3')  nm = '1HE '
    if (nm.eq.'3HE ')  nm = '1HE '
    if (nm(1:3).eq.'HNZ') nm(1:3) = ' HZ'
    if (resname.eq.'KAC') then
      if (nm.eq.' HI1')  nm = '1HI '
      if (nm.eq.' HI2')  nm = '2HI '
      if (nm.eq.' HI3')  nm = '3HI '
    else if (resname.eq.'KM1') then
      if (nm.eq.' HH1')  nm = '1HH '
      if (nm.eq.' HH2')  nm = '2HH '
      if (nm.eq.' HH3')  nm = '3HH '
      if (nm.eq.' HZ1')  nm = '1HZ '
      if (nm.eq.' HZ2')  nm = '2HZ '
      if (nm.eq.' HZ3')  nm = '1HZ '
      if (nm.eq.'3HZ ')  nm = '1HZ '
    else if (resname.eq.'KM2') then
      if (nm.eq.'HH11')  nm = '1HH1'
      if (nm.eq.'HH12')  nm = '2HH1'
      if (nm.eq.'HH13')  nm = '3HH1'
      if (nm.eq.'HH21')  nm = '1HH2'
      if (nm.eq.'HH22')  nm = '2HH2'
      if (nm.eq.'HH23')  nm = '3HH2'
      if (nm.eq.' HZ1')  nm = ' HZ '
      if (nm.eq.' HZ3')  nm = ' HZ '
    else if (resname.eq.'KM3') then
      if (nm.eq.'HH11')  nm = '1HH1'
      if (nm.eq.'HH12')  nm = '2HH1'
      if (nm.eq.'HH13')  nm = '3HH1'
      if (nm.eq.'HH21')  nm = '1HH2'
      if (nm.eq.'HH22')  nm = '2HH2'
      if (nm.eq.'HH23')  nm = '3HH2'
      if (nm.eq.'HH31')  nm = '1HH3'
      if (nm.eq.'HH32')  nm = '2HH3'
      if (nm.eq.'HH33')  nm = '3HH3'
    else
      if (nm.eq.' HZ1')  nm = '1HZ '
      if (nm.eq.' HZ2')  nm = '2HZ '
      if (resname.eq.'LYS') then
        if (nm.eq.' HZ3')  nm = '3HZ '
      else if (resname.eq.'LYD') then
        if (nm.eq.' HZ3')  nm = '1HZ '
        if (nm.eq.'3HZ ')  nm = '1HZ '
      end if
    end if
!
! ARG
  else if (resname.eq.'ARG') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '1HD '
    if (nm.eq.'3HD ')  nm = '1HD '
    if (nm.eq.'HH11')  nm = '1HH1'
    if (nm.eq.'HN11')  nm = '1HH1'
    if (nm.eq.'HH12')  nm = '2HH1'
    if (nm.eq.'HN12')  nm = '2HH1'
    if (nm.eq.'HH21')  nm = '1HH2'
    if (nm.eq.'HN21')  nm = '1HH2'
    if (nm.eq.'HH22')  nm = '2HH2'
    if (nm.eq.'HN22')  nm = '2HH2'
!
! HIP/HIS (for the fuzzy HIS, superset of HID/E/P)
  else if ((resname.eq.'HIP').OR.(resname.eq.'HIS')) then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HD ')  nm = ' HD2'
    if (nm.eq.' HE ')  nm = ' HE1'
    if (nm.eq.' HND')  nm = ' HD1'
    if (nm.eq.' HNE')  nm = ' HE2'
!
! SEP
  else if (resname.eq.'SEP') then
    if (nm.eq.' OG1')  nm = ' OG '
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' OP1')  nm = ' O1P'
    if (nm.eq.' OP2')  nm = ' O2P'
    if (nm.eq.' HO ')  nm = ' HOP'
    if (nm.eq.' HT ')  nm = ' HOP'
    if (nm.eq.' OT ')  nm = ' OH '
!
! TPO
  else if (resname.eq.'TPO') then
    if (nm.eq.' OG ')  nm = ' OG1'
    if (nm.eq.' CG ')  nm = ' CG2'
    if (nm.eq.'HG21')  nm = '1HG2'
    if (nm.eq.'HG22')  nm = '2HG2'
    if (nm.eq.'HG23')  nm = '3HG2'
    if (nm.eq.' OP1')  nm = ' O1P'
    if (nm.eq.' OP2')  nm = ' O2P'
    if (nm.eq.' HO ')  nm = ' HOP'
    if (nm.eq.' HT ')  nm = ' HOP'
    if (nm.eq.' OT ')  nm = ' OH '
!
! PTR
  else if (resname.eq.'PTR') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' OP1')  nm = ' O1P'
    if (nm.eq.' OP2')  nm = ' O2P'
    if (nm.eq.' HO ')  nm = ' HOP'
    if (nm.eq.' HT ')  nm = ' HOP'
    if (nm.eq.' OT ')  nm = ' OH '
!
! ORN
  else if (resname.eq.'ORN') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '1HD '
    if (nm.eq.'3HD ')  nm = '1HD '
    if (nm.eq.' HE1')  nm = '1HE '
    if (nm.eq.'HNE1')  nm = '1HE '
    if (nm.eq.' HE2')  nm = '2HE '
    if (nm.eq.'HNE2')  nm = '2HE '
    if (nm.eq.' HE3')  nm = '3HE '
    if (nm.eq.'HNE3')  nm = '3HE '
!
! ABA
  else if (resname.eq.'ABA') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '3HG '
!
! NVA
  else if (resname.eq.'NVA') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '3HD '
!
! NLE
  else if (resname.eq.'NVA') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '1HD '
    if (nm.eq.'3HD ')  nm = '1HD '
    if (nm.eq.' HE1')  nm = '1HE '
    if (nm.eq.' HE2')  nm = '2HE '
    if (nm.eq.' HE3')  nm = '3HE '
!
! DAB
  else if (resname.eq.'DAB') then
    if (nm.eq.' HB1')  nm = '1HB '
    if (nm.eq.' HB2')  nm = '2HB '
    if (nm.eq.' HB3')  nm = '1HB '
    if (nm.eq.'3HB ')  nm = '1HB '
    if (nm.eq.' HG1')  nm = '1HG '
    if (nm.eq.' HG2')  nm = '2HG '
    if (nm.eq.' HG3')  nm = '1HG '
    if (nm.eq.'3HG ')  nm = '1HG '
    if (nm.eq.' HD1')  nm = '1HD '
    if (nm.eq.'HND1')  nm = '1HD '
    if (nm.eq.' HD2')  nm = '2HD '
    if (nm.eq.'HND2')  nm = '2HD '
    if (nm.eq.' HD3')  nm = '3HD '
    if (nm.eq.'HND3')  nm = '3HD '
!
! small molecules
! URE
  else if (resname.EQ.'URE') then
    if (nm.eq.' N1U')  nm = ' N1 '
    if (nm.eq.' N2U')  nm = ' N2 '
    if (nm.eq.' CU ')  nm = ' C '
    if (nm.eq.' OU ')  nm = ' O  '
    if (nm.eq.' H11')  nm = '1HN1'
    if (nm.eq.' H12')  nm = '2HN1'
    if (nm.eq.' H21')  nm = '1HN2'
    if (nm.eq.' H22')  nm = '2HN2'
    if (nm.eq.'H11U')  nm = '1HN1'
    if (nm.eq.'H12U')  nm = '2HN1'
    if (nm.eq.'H21U')  nm = '1HN2'
    if (nm.eq.'H22U')  nm = '2HN2'
    if (nm.eq.'HN11')  nm = '1HN1'
    if (nm.eq.'HN12')  nm = '2HN1'
    if (nm.eq.'HN21')  nm = '1HN2'
    if (nm.eq.'HN22')  nm = '2HN2'
    if (nm.eq.'1H1 ')  nm = '1HN1'
    if (nm.eq.'2H1 ')  nm = '2HN1'
    if (nm.eq.'1H2 ')  nm = '1HN2'
    if (nm.eq.'2H2 ')  nm = '2HN2'
!
! GDN
  else if (resname.EQ.'GDN') then
    if (nm.eq.' H11')  nm = '1HN1'
    if (nm.eq.' H12')  nm = '2HN1'
    if (nm.eq.' H21')  nm = '1HN2'
    if (nm.eq.' H22')  nm = '2HN2'
    if (nm.eq.' H31')  nm = '1HN3'
    if (nm.eq.' H32')  nm = '2HN3'
    if (nm.eq.'HN11')  nm = '1HN1'
    if (nm.eq.'HN12')  nm = '2HN1'
    if (nm.eq.'HN21')  nm = '1HN2'
    if (nm.eq.'HN22')  nm = '2HN2'
    if (nm.eq.'HN31')  nm = '1HN3'
    if (nm.eq.'HN32')  nm = '2HN3'
    if (nm.eq.'1H1 ')  nm = '1HN1'
    if (nm.eq.'2H1 ')  nm = '2HN1'
    if (nm.eq.'1H2 ')  nm = '1HN2'
    if (nm.eq.'2H2 ')  nm = '2HN2'
    if (nm.eq.'1H3 ')  nm = '1HN3'
    if (nm.eq.'2H3 ')  nm = '2HN3'
!
! NH4
  else if (resname.eq.'NH4') then
    if (nm.eq.'1H  ')  nm = '1HN '
    if (nm.eq.'2H  ')  nm = '2HN '
    if (nm.eq.'3H  ')  nm = '3HN '
    if (nm.eq.'4H  ')  nm = '4HN '
    if (nm.eq.' H1 ')  nm = '1HN '
    if (nm.eq.' H2 ')  nm = '2HN '
    if (nm.eq.' H3 ')  nm = '3HN '
    if (nm.eq.' H4 ')  nm = '4HN '
    if (nm.eq.' HN1')  nm = '1HN '
    if (nm.eq.' HN2')  nm = '2HN '
    if (nm.eq.' HN3')  nm = '3HN '
    if (nm.eq.' HN4')  nm = '4HN '
!
! Acetate
  else if (resname.eq.'AC-') then
    if (nm.eq.' H1 ')  nm = '1H  '
    if (nm.eq.' H2 ')  nm = '2H  '
    if (nm.eq.' H3 ')  nm = '3H  '
    if (nm.eq.' HC1')  nm = '1H  '
    if (nm.eq.' HC2')  nm = '2H  '
    if (nm.eq.' HC3')  nm = '3H  '
!
! Methylammonium
  else if (resname.eq.'1MN') then
    if (nm.eq.' H1 ')  nm = '1HN '
    if (nm.eq.' H2 ')  nm = '2HN '
    if (nm.eq.' H3 ')  nm = '3HN '
    if (nm.eq.'1H  ')  nm = '1HN '
    if (nm.eq.'2H  ')  nm = '2HN '
    if (nm.eq.'3H  ')  nm = '3HN '
    if (nm.eq.' HC1')  nm = '1HT '
    if (nm.eq.' HC2')  nm = '2HT '
    if (nm.eq.' HC3')  nm = '3HT '
    if (nm.eq.' HT1')  nm = '1HT '
    if (nm.eq.' HT2')  nm = '2HT '
    if (nm.eq.' HT3')  nm = '3HT '
!
! Dimethylammonium
  else if (resname.eq.'2MN') then
    if (nm.eq.' H1 ')  nm = '1HN '
    if (nm.eq.' H2 ')  nm = '2HN '
    if (nm.eq.'1H  ')  nm = '1HN '
    if (nm.eq.'2H  ')  nm = '2HN '
    if (nm.eq.'HC11')  nm = '1HT1'
    if (nm.eq.'HC12')  nm = '2HT1'
    if (nm.eq.'HC13')  nm = '3HT1'
    if (nm.eq.'HT11')  nm = '1HT1'
    if (nm.eq.'HT12')  nm = '2HT1'
    if (nm.eq.'HT13')  nm = '3HT1'
    if (nm.eq.'HC21')  nm = '1HT2'
    if (nm.eq.'HC22')  nm = '2HT2'
    if (nm.eq.'HC23')  nm = '3HT2'
    if (nm.eq.'HT21')  nm = '1HT2'
    if (nm.eq.'HT22')  nm = '2HT2'
    if (nm.eq.'HT23')  nm = '3HT2'
!
! Inorganic ions
  else if ((resname.eq.'K+ ').OR.(resname.eq.'NA+').OR.(resname.eq.'CL-').OR.(resname.eq.'CS+')) then
    if (nm.eq.' SOD') nm = ' NA '
    if (nm.eq.' POT') nm = ' K  '
    if (nm.eq.' CES') nm = ' CS '
    if (nm.eq.' CLA') nm = ' CL '
!
! Waters
  else if ((resname.eq.'SPC').OR.(resname.eq.'T3P')) then
    if (nm.eq.' HW1')  nm = '1HW '
    if (nm.eq.' HW2')  nm = '2HW '
    if (nm.eq.'1H  ')  nm = '1HW '
    if (nm.eq.'2H  ')  nm = '2HW '
    if (nm.eq.' H1 ')  nm = '1HW '
    if (nm.eq.' H2 ')  nm = '2HW '
    if (nm.eq.' OH2')  nm = ' OW '
!
! Tip4/5p
  else if ((resname.eq.'T4P').OR.(resname.eq.'T4E').OR.(resname.eq.'T5P')) then
    if (nm.eq.' HW1')  nm = '1HW '
    if (nm.eq.' HW2')  nm = '2HW '
    if (nm.eq.'1H  ')  nm = '1HW '
    if (nm.eq.'2H  ')  nm = '2HW '
    if (nm.eq.' H1 ')  nm = '1HW '
    if (nm.eq.' H2 ')  nm = '2HW '
    if (nm.eq.' OH2')  nm = ' OW '
    if (nm.eq.' OM ')  nm = ' MW '
    if ((resname.ne.'T5P').AND.(nm.eq.' HW3'))  nm = ' MW '
    if (nm.eq.' LP1')  nm = '1LP '
    if (nm.eq.' LP2')  nm = '2LP '
    if ((resname.eq.'T5P').AND.(nm.eq.' EP1'))  nm = '1LP '
    if ((resname.ne.'T5P').AND.(nm.eq.' EP1'))  nm = ' MW '
    if (nm.eq.' EP2')  nm = '2LP '
    if ((resname.eq.'T5P').AND.(nm.eq.'1EP '))  nm = '1LP '
    if ((resname.ne.'T5P').AND.(nm.eq.'1EP '))  nm = ' MW '
    if (nm.eq.'2EP ')  nm = '2LP '
    if (nm.eq.' EPW')  nm = ' MW '
    if (nm.eq.' EP ')  nm = ' MW '
    if ((resname.eq.'T5P').AND.(nm.eq.' MW1'))  nm = '1LP '
    if ((resname.ne.'T5P').AND.(nm.eq.' MW1'))  nm = ' MW '
    if (nm.eq.' MW2')  nm = '2LP '
    if ((resname.eq.'T5P').AND.(nm.eq.'1MW '))  nm = '1LP '
    if ((resname.ne.'T5P').AND.(nm.eq.'1MW '))  nm = ' MW '
    if (nm.eq.'2MW ')  nm = '2LP '
    if ((resname.eq.'T5P').AND.(nm.eq.' HW3'))  nm = '1LP '
    if (nm.eq.' HW4')  nm = '2LP '
!
! methanol
  else if (resname.eq.'MOH') then
    if (nm.eq.'CMET')  nm = ' CT '
    if (nm.eq.'OMET')  nm = ' O  '
    if (nm.eq.'HMET')  nm = ' HO '
    if (nm.eq.' HC1')  nm = '1H  '
    if (nm.eq.' HC2')  nm = '2H  '
    if (nm.eq.' HC3')  nm = '3H  '
    if (nm.eq.' HT1')  nm = '1H  '
    if (nm.eq.' HT2')  nm = '2H  '
    if (nm.eq.' HT3')  nm = '3H  '
    if (nm.eq.'1HC ')  nm = '1H  '
    if (nm.eq.'2HC ')  nm = '2H  '
    if (nm.eq.'3HC ')  nm = '3H  '
!
  end if
!
! nucleotides
  if (nucat.EQV..true.) then
!   sugar
    if (nm.eq." H5'") nm = '1H5*'
    if (nm.eq."H5''") nm = '2H5*'
    if (nm.eq."H5'1") nm = '1H5*'
    if (nm.eq."H5'2") nm = '2H5*'
    if (nm.eq."1H5'") nm = '1H5*'
    if (nm.eq."2H5'") nm = '2H5*'
!   DNA version of H2 hydrogens
    if (nm.eq." H2'") nm = '1H2*'
    if (nm.eq."H2''") nm = '2H2*'
    if (nm.eq."H2'2") nm = '2H2*'
    if (nm.eq."2H2'") nm = '2H2*'
    if (resname(1:1).eq.'D') then
      if (nm.eq."H2'1") nm = '1H2*'
      if (nm.eq."1H2'") nm = '1H2*'
    else if (resname(1:1).eq.'R') then
      if (nm.eq."H2'1") nm = ' H2*'
      if (nm.eq."1H2'") nm = ' H2*'
    end if
!   rest
    if (nm.eq."2HO'") nm = '2HO*'
    if (ua_model.gt.0) then
      if (nm.eq.' H2*') nm = '2HO*'
      if (nm.eq." H2'") nm = '2HO*'
    end if
    if (ua_model.gt.0) then
      if (nm.eq.' H3*') nm = '3HO*'
      if (nm.eq." H3'") nm = '3HO*'
    end if
    if (nm.eq."HO2'") nm = '2HO*'
    if (nm.eq."HO'2") nm = '2HO*'
    if (nm.eq." H1'") nm = ' H1*'
    if (nm.eq." H3'") nm = ' H3*'
    if (nm.eq." H4'") nm = ' H4*'
    if (nm.eq." C1'") nm = ' C1*'
    if (nm.eq." C2'") nm = ' C2*'
    if (nm.eq." O2'") nm = ' O2*'
    if (nm.eq." C3'") nm = ' C3*'
    if (nm.eq." C4'") nm = ' C4*'
    if (nm.eq." C5'") nm = ' C5*'
    if (nm.eq." O4'") nm = ' O4*'
!   phosphate
    if (nm.eq." O5'") nm = ' O5*'
    if (nm.eq." OP1") nm = ' O1P'
    if (nm.eq." OP2") nm = ' O2P'
    if (nm.eq." O3'") nm = ' O3*'
    if (nm.eq." C3'") nm = ' C3*'
    if ((pdb_convention(2).eq.3).OR.(pdb_convention(2).eq.4)) then
!      if ((resname(2:2).eq."I").OR.((pdb_convention(2).eq.3).AND.(seqflag(rs).eq.24))) then
!        if (nm.eq.' H5T') nm = '5HO*'
!      else
        if (nm.eq.' H5T') nm = 'XHOP'
!      end if
      if (nm.eq.' O5T') nm = ' O3*'
      if (nm.eq.' H3T') nm = '3HO*'
      if (nm.eq.'HO5*') nm = '5HO*'
      if (nm.eq.'HO3*') nm = '3HO*'
      if (nm.eq.' HO5') nm = '5HO*'
      if (nm.eq.' HO3') nm = '3HO*'
    end if
!   bases
    if (nm.eq.' H21')  nm = '1H2 '
    if (nm.eq.' H22')  nm = '2H2 '
    if (nm.eq.' H41')  nm = '1H4 '
    if (nm.eq.' H42')  nm = '2H4 '
    if (nm.eq.' H61')  nm = '1H6 '
    if (nm.eq.' H62')  nm = '2H6 '
    if (nm.eq.' H51')  nm = '1H5M'
    if (nm.eq.' H52')  nm = '2H5M'
    if (nm.eq.' H53')  nm = '3H5M'
    if (resname(3:3).eq.'T') then
      if (nm.eq.' C7 ')  nm = ' C5M'
      if (nm.eq.' H71')  nm = '1H5M'
      if (nm.eq.' H72')  nm = '2H5M'
      if (nm.eq.' H73')  nm = '3H5M'
      if (nm.eq.'1H7 ')  nm = '1H5M'
      if (nm.eq.'2H7 ')  nm = '2H5M'
      if (nm.eq.'3H7 ')  nm = '3H5M'
    end if
  end if
!
end
!
!----------------------------------------------------------------------------------------------------
!
function seqt_from_name(iname,start)
!
  use aminos
  use pdb, ONLY: pdb_unknowns,n_pdbunk
!
  implicit none
!
  integer i,seqt_from_name,start
  character(3) iname
!
  do i=1,MAXAMINO
    if (i.eq.26) cycle
    if (iname.eq.amino(i)) then
      seqt_from_name = i
      return
    end if
  end do
  if (iname.eq.'HIS') then
    seqt_from_name = 99999
    return
  end if
!
  do i=start,n_pdbunk
    if (iname.eq.pdb_unknowns(i)(1:3)) then
      seqt_from_name = -i
      start = i+1
      return
    end if
  end do
  seqt_from_name = 0
!
end
!
!--------------------------------------------------------------------------------------
!
! this function reads in a pdb-file and attempts to extract values for CAMPARI's 
! native degrees of freedom from the xyz-information
! because of the rigidity assumption, this is a non-trivial task and will
! often lead to xyz-mismatches which become increasingly worse the longer the
! polymers are
! straight xyz-porting and re-building of the Z-matrix is done by FMCSC_readpdb2()
!
subroutine FMSMC_readpdb()
!
  use iounit
  use math
  use atoms
  use sequen
  use polypep
  use fyoc
  use aminos
  use molecule
  use zmatrix
  use pdb
  use system
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer ipdb,nl,i,j,k,sta,sto,ats,ate,t1,t2,freeunit,rs,rs2,flagrsshft,st1,st2
  integer stat,strs,jupp,firs,rsshft,ii,jj,currs,lastrs
  character(MAXSTRLEN), ALLOCATABLE:: pdbl(:)
  character(6), ALLOCATABLE:: kw(:)
  character(6) keyw
  character(4), ALLOCATABLE:: pdbnm(:)
  character(3), ALLOCATABLE:: rsnam(:)
  character(3) resname
  character(25) atnmb
  character(1) extrabit,refbit,extrabit2,refbit2
  character, ALLOCATABLE:: chnam(:)
  character(MAXSTRLEN) fn
  character(MAXSTRLEN+50) fpdbstr
  integer, ALLOCATABLE:: rsnmb(:),rsshftl(:),lix(:)
  integer ok(120),ok2(120)
  RTYPE, ALLOCATABLE:: xr(:),yr(:),zr(:)
  type(t_pdbxyzs) pdbxyzs
  logical exists,ribon,foundp
#ifdef ENABLE_MPI
  integer tslash
  character(re_aux(10)) xpont
#endif
!
! allocation first
  allocate(pdbl(3*sum(natres(1:nseq))+2000))
  allocate(pdbnm(3*sum(natres(1:nseq))+100))
  allocate(rsnam(3*sum(natres(1:nseq))+100))
  allocate(chnam(3*sum(natres(1:nseq))+100))
  allocate(kw(3*sum(natres(1:nseq))+100))
  allocate(rsnmb(3*sum(natres(1:nseq))+100))
  allocate(rsshftl(3*sum(natres(1:nseq))+100))
  allocate(xr(3*sum(natres(1:nseq))+100))
  allocate(yr(3*sum(natres(1:nseq))+100))
  allocate(zr(3*sum(natres(1:nseq))+100))
  allocate(lix(3*sum(natres(1:nseq))+100))
!
  do j=1,10
    ok2(j) = 0
  end do
!
  ipdb = freeunit()
  call strlims(pdbinfile,t1,t2)
#ifdef ENABLE_MPI
  if (pdb_mpimany.EQV..true.) then
    do i=t2,t1,-1
      if (pdbinfile(i:i).eq.SLASHCHAR) exit
    end do
    tslash = i
    call int2str(myrank,xpont,re_aux(10))
    if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
      fn =  pdbinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//pdbinfile((tslash+1):t2)
    else
      fn =  'N_'//xpont(1:re_aux(10))//'_'//pdbinfile(t1:t2)
    end if
    call strlims(fn,t1,t2)
  else
    fn(t1:t2) = pdbinfile(t1:t2)
  end if
#else
  fn(t1:t2) = pdbinfile(t1:t2)
#endif
  inquire (file=fn(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading from pdb (',fn(t1:t2),').'
    call fexit()
  end if
!
  open (unit=ipdb,file=fn(t1:t2),status='old')
!
  k = 0
  do while (1.eq.1)
    k = k + 1
    if (k.ge.(3*sum(natres(1:nseq))+2000)) then
      write(ilog,*) 'Fatal. Spec.d input file is too large for reading from pdb (',fn(t1:t2),').'
      write(ilog,*) 'If it is a trajectory file, cut out the snapshot to be used.'
      call fexit()
    end if
 1356 format(FORM_MAXSTRLEN)
    read(ipdb,1356,err=30,end=30) pdbl(k)        
  end do
 30   k = k - 1
  close(unit=ipdb)
!
  nl = k
  sta = -1
  sto = -1
  do i=1,nl
    read(pdbl(i),'(a6)') keyw
    if ((keyw.eq.'MODEL ').AND.(sta.eq.-1)) then
      write(ilog,*) 'WARNING. PDB fle seems to contain multiple models. Using first one only (',fn(t1:t2),').'
      sta = i
    else if ((keyw.eq.'ENDMDL').AND.(sta.ne.-1)) then
      sto = i
      exit
    else if (keyw.eq.'ENDMDL') then
      write(ilog,*) 'Fatal. Encountered ENDMDL before MODEL-keyword while reading from PDB (',fn(t1:t2),').'
      call fexit()
    end if
  end do
!
  ats = -1
  ate = -1
  if ((sta.ne.-1).AND.(sto.ne.-1)) then
    do i=sta,sto
      read(pdbl(i),'(a6)') keyw
      if (((keyw.eq.'ATOM  ').OR.&
 &         (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
        ats = i
      else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &         (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
        ate = i-1
        exit
      end if
    end do
  else if ((sta.eq.-1).AND.(sto.eq.-1)) then
    do i=1,nl
      read(pdbl(i),'(a6)') keyw
      if (((keyw.eq.'ATOM  ').OR.&
 &         (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
        ats = i
      else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &         (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
        ate = i-1
        exit
      end if
    end do
  else
    write(ilog,*) 'Fatal. Could not find ENDMDL keyword despite MODEL keyword while reading from PDB (',fn(t1:t2),').'
    call fexit()
  end if
  if ((ats.le.0).AND.(ate.le.0)) then
    write(ilog,*) 'Fatal. PDB file (',fn(t1:t2),') appears to contain no coordinate information whatsoever. Please check input.'
    call fexit()
  end if
!
! this gives us the frame, within which to operate
  if ((ate-ats+1).ge.(3*sum(natres(1:nseq))+100)) then
    write(ilog,*) 'Fatal. PDB-File has too many atoms. Check for unwanted water molecules or a similar problem (',fn(t1:t2),').'
    call fexit()
  end if
!
  
  call strlims(pdb_formstr(2),st1,st2)
  fpdbstr = '('//pdb_formstr(2)(st1:st2)//')'
    34 format (a6,a5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3) ! the default format
!
  k = 0
  rsshft = 0
  rsshftl(:) = 0
  currs = huge(currs)
  refbit = ' '
  refbit2 = ' '
  do i=ats,ate
    k = k + 1
    read(pdbl(i),'(a6)') keyw
    if ((keyw.eq.'TER   ').OR.(keyw.eq.'ANISOU')) then
      k = k - 1
      cycle
    end if
    read(pdbl(i),fpdbstr) kw(k),atnmb,pdbnm(k),extrabit,rsnam(k),chnam(k),rsnmb(k),extrabit2,xr(k),yr(k),zr(k)
    lix(k) = i
!   the currs/lastrs construct is (at least sporadically) able to read in pdb-files in which
!   the residue numbering is off as long as it is a number and the number for a new residue is always
!   different from the previous one
    if (k.eq.1) lastrs = rsnmb(k)-1
    if (currs.ne.huge(currs)) then
      if ((rsnmb(k).ne.currs).OR.(extrabit2.ne.refbit2)) then
        refbit = ' '
        if (currs.ne.lastrs+1) currs = lastrs+1
        lastrs = currs
        currs = rsnmb(k)
        if (rsnmb(k).ne.lastrs+1) then
          rsnmb(k) = lastrs + 1
        end if
      else
        if (rsnmb(k).ne.lastrs+1) rsnmb(k) = lastrs+1
      end if
    end if
    if (k.eq.1) currs = rsnmb(k)
    if ((refbit.eq.' ').AND.(extrabit.ne.' ')) refbit = extrabit
    refbit2 = extrabit2
    call toupper(kw(k))
!   drop the records for TER. ANISOU, and any alternate conformations not found first
    if ((kw(k).eq.'TER   ').OR.(kw(k).eq.'ANISOU').OR.((extrabit.ne.' ').AND.(extrabit.ne.refbit))) then
      k = k - 1
      cycle
    end if
    call toupper(pdbnm(k))
    call toupper(rsnam(k))
    call toupper(chnam(k))
    call pdbcorrect(rsnam(k),pdbnm(k),rsnmb(k),rsshft,rsshftl)
  end do
!
  ii = 1
  do j=1,rsshft
    flagrsshft = 0
    do i=ii,k
      if (rsnmb(i).eq.rsshftl(j)) then
        do jj=i,k
          if (rsnmb(jj).ne.rsshftl(j)) then
            ii = jj
            exit
          end if
          if ((rsnam(jj).eq.'NH2').OR.(rsnam(jj).eq.'NME')) then
            rsnmb(jj) = rsnmb(jj) + 1
            if (flagrsshft.eq.0) flagrsshft = 1
            if (flagrsshft.eq.-1) flagrsshft = 2
          end if
          if ((rsnam(jj).eq.'FOR').OR.(rsnam(jj).eq.'ACE')) then
            rsnmb(jj) = rsnmb(jj) - 1
            if (flagrsshft.eq.0) flagrsshft = -1
            if (flagrsshft.eq.1) flagrsshft = 2
          end if
          if (jj.eq.k) ii = k + 1
        end do
        if (flagrsshft.eq.1) then
          if (ii.le.k) then
            rsnmb(ii:k) = rsnmb(ii:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 1
          end if
        else if (flagrsshft.eq.-1) then
          if (i.le.k) then
            rsnmb(i:k) = rsnmb(i:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 1
          end if
        else if (flagrsshft.eq.2) then
          if (i.le.k) then
            rsnmb(i:k) = rsnmb(i:k) + 1
          end if
          if (ii.le.k) then
            rsnmb(ii:k) = rsnmb(ii:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 2
          end if
        end if
        exit ! skip out, go to next shft
      end if
    end do
  end do
!
  jupp = -1000
  do i=1,k
    if (rsnmb(i).ne.jupp) then
      if (i.gt.1) then
        ribon = .false.
        foundp = .false.
        do j=firs,i-1
          if (pdbnm(j).eq.' P  ') then ! crucial that P is always called P
            foundp = .true.
            exit
          end if
        end do
        do j=firs,i-1
          if (rsnam(j)(1:2).eq.'R ') then
            ribon = .true.
          end if
        end do
        if (ribon.EQV..true.) then
          do j=firs,i-1
!           hopeful
            if (foundp.EQV..true.) then
              rsnam(j)(1:2) = 'RP'
            else
              rsnam(j)(1:2) = 'RI'
            end if
            if (pdb_convention(2).eq.3) then
              if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = '2HO*'
              if (pdbnm(j)(1:4).eq.'2H2*') pdbnm(j)(1:4) = ' H2*'
            else
              if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = ' H2*'
            end if
          end do
        else if (rsnam(firs)(1:2).eq.'DP') then
          if (foundp.EQV..false.) then
            do j=firs,i-1
              rsnam(j)(2:2) = 'I'
            end do
          end if
        end if
        if ((ribon.EQV..true.).OR.(rsnam(firs)(1:2).eq.'DP')) then
          if (foundp.EQV..true.) then
            do j=firs,i-1          
              if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = '5HO*'
            end do
          else
            do j=firs,i-1  
              if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = ' HOP'
            end do
          end if
        end if
      end if
      firs = i
      jupp = rsnmb(i)
    end if
    if (i.eq.k) then
      ribon = .false.
      foundp = .false.
      do j=firs,k
        if (pdbnm(j).eq.' P  ') then ! crucial that P is always called P
          foundp = .true.
          exit
        end if
      end do
      do j=firs,k
        if (rsnam(j)(1:2).eq.'R ') then
          ribon = .true.
        end if
      end do
      if (ribon.EQV..true.) then
        do j=firs,k
          if (foundp.EQV..true.) then
            rsnam(j)(1:2) = 'RP'
          else
            rsnam(j)(1:2) = 'RI'
          end if
          if (pdb_convention(2).eq.3) then
            if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = '2HO*'
            if (pdbnm(j)(1:4).eq.'2H2*') pdbnm(j)(1:4) = ' H2*'
          else
            if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = ' H2*'
          end if
        end do
      else if (rsnam(firs)(1:2).eq.'DP') then
        if (foundp.EQV..false.) then
          do j=firs,k
            rsnam(j)(2:2) = 'I'
          end do
        end if
      end if
      if ((ribon.EQV..true.).OR.(rsnam(firs)(1:2).eq.'DP')) then
        if (foundp.EQV..true.) then
          do j=firs,i-1          
            if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = '5HO*'
          end do
        else
          do j=firs,i-1  
            if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = ' HOP'
          end do
        end if
      end if
    end if
  end do
!
  nl = k
!
  call adjust_PONuc(nl,rsnmb(1:nl),pdbnm(1:nl),rsnam(1:nl),xr(1:nl),yr(1:nl),zr(1:nl))
!
! this function reorders the input arrays (except pdbl) so that everything is contiguous per residue
  call reorder_PDB_byrsnmb(nl,kw(1:nl),rsnmb(1:nl),rsnam(1:nl),pdbnm(1:nl),xr(1:nl),yr(1:nl),zr(1:nl),&
 &lix(1:nl),chnam(1:nl))
!
! now let's check whether the sequence matches
  call strlims(seqfile,st1,st2)
  rs = 1
  resname = amino(seqtyp(rs))
  i = 1
  do while ((resname.ne.rsnam(i)).AND.&
 & (.NOT.((rsnam(i).eq.'HIS').AND.((resname.eq.'HID').OR.(resname.eq.'HIE').OR.(resname.eq.'HIP'))))) 
    i = i+1
    if (i.gt.k) then
      write(ilog,*) 'Fatal. Cannot find first residue (',resname,') &
 &in PDB-file (',fn(t1:t2),').'
      call fexit()
    end if
  end do
  if (rsnam(i).eq.'HIS') then
    do j=i+1,k
      if (j.gt.1) then
        if (rsnmb(j).ne.rsnmb(j-1)) exit
      end if
      rsnam(j) = resname
    end do
  end if
  strs = rsnmb(i)
  rs2 = strs
  stat = i
  do while (rs.lt.nseq)
    rs = rs + 1
    rs2 = rs2 + 1
    resname = amino(seqtyp(rs))
    do i=1,k
      if (rsnmb(i).eq.rs2) then
        if ((rsnam(i).eq.resname).OR.&
 &  ((rsnam(i).eq.'HIS').AND.((resname.eq.'HID').OR.(resname.eq.'HIE').OR.(resname.eq.'HIP')))) then
          if (rsnam(i).eq.'HIS') then
            do j=i+1,k
              if (j.gt.1) then
                if (rsnmb(j).ne.rsnmb(j-1)) exit
              end if
              rsnam(j) = resname
            end do
          end if
          exit
        else
          write(ilog,*) 'Fatal. Sequences in sequence file and PDB-file do not match.'
          write(ilog,*) 'Seq. file: ',seqfile(st1:st2),', PDB-file: ',fn(t1:t2)
          write(ilog,*) 'Position ',rs,' in seq. file (',resname,') and position ',&
 &                       rsnmb(i),' in PDB-file (',rsnam(i),')!'
          call fexit()
        end if
      end if
    end do
  end do
!
! if we're still here, we can finally proceed to extract the rigid body-coordinates and torsional angles
  rs = 1
  rs2 = strs
  do while (rs.le.nseq)
    ok(:) = 0
    
    resname = amino(seqtyp(rs))
    do i=stat,k
      if (rsnmb(i).ne.rs2) then
        stat = i
        exit
      end if
!
!     here we're simply collecting the xyz-coordinates of relevant atoms (as matched
!     up by pdb atom-name) into a specialized array (pdbxyzs -> pdb.i)
!
      call readpdb1_collect(stat,i,pdbnm,xr,yr,zr,pdbxyzs,ok)
!
    end do

!   shift for pdb-format nucleotides (3'-O is part of sugar-res.)
    if ((ok(66).eq.0).AND.(ok2(7).eq.1)) then
      ok(66) = 1
      ok2(7) = 0
      do j=1,3
        pdbxyzs%no3(j) = pdbxyzs%no32l(j)
      end do
    else if ((ok(66).eq.0).AND.(ok2(6).eq.1)) then
      if ((seqtyp(rs-1).ge.76).AND.(seqtyp(rs-1).le.87)) then
        ok(66) = 1
        ok2(6) = 0
        do j=1,3
          pdbxyzs%no3(j) = pdbxyzs%no3l(j)
        end do
      end if
    end if
!
!   now with the presence-flags (ok), and the xyz-information all present, we can
!   proceed to extract the internal coordinates from the latter
!
!   rigid-body coordinates first
    call readpdb1_rbs(rs,resname,pdbxyzs,ok)
    if (maxval(ok).le.-1) then
      lastrs = rs - 1
      exit
    else
      lastrs = rs
    end if
!
!
!   now chi-angles
    call readpdb1_chis(rs,resname,pdbxyzs,ok,ok2)
!
!
!   now omega-angles for caps
    call readpdb1_omegas(rs,resname,pdbxyzs,ok,ok2)
!
!   cycle out immed. for all single-residue molecules and NME/NH2
    if (((rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1).eq.0).AND.&
 &       (seqpolty(rs).ne.'N').AND.(seqpolty(rs).ne.'P')).OR.&
 &      (resname.eq.'NH2').OR.(resname.eq.'NME')) then
      rs = rs + 1
      rs2 = rs2 + 1
      cycle
    end if

!   and remaining backbone angles
    call readpdb1_bbs(rs,resname,pdbxyzs,ok,ok2)
!
!
!   finally, shift the backwards counters
    do j=1,3
      if (resname.eq.'FOR') then
        if (ok(20).eq.1) then
          pdbxyzs%cl(j) = pdbxyzs%c(j)
          ok2(3) = 1
        else
          ok2(3) = 0
        end if
        if (ok(27).eq.1) then
          pdbxyzs%ol(j) = pdbxyzs%o(j)
          ok2(9) = 1
        else
          ok2(9) = 0
        end if
!       in UA, formyl is only CO
        if (ua_model.eq.0) then
          if (ok(29).eq.1) then
            pdbxyzs%cal(j) = pdbxyzs%h(j)
            ok2(2) = 1
          else
            ok2(2) = 0
          end if
        end if
      else if (resname.eq.'ACE') then
        if (ok(20).eq.1) then
          pdbxyzs%cl(j) = pdbxyzs%c(j)
          ok2(3) = 1
        else
          ok2(3) = 0
        end if
        if (ok(30).eq.1) then
          pdbxyzs%cal(j) = pdbxyzs%ch3(j)
          ok2(2) = 1
        else
          ok2(2) = 0
        end if
        if (ok(27).eq.1) then
          pdbxyzs%ol(j) = pdbxyzs%o(j)
          ok2(9) = 1
        else
          ok2(9) = 0
        end if
      else if (seqpolty(rs).eq.'N') then
        if (ok(65).eq.1) then
          pdbxyzs%nc3l(j) = pdbxyzs%nc3(j)
          ok2(4) = 1
        else
          ok2(4) = 0
        end if
        if (ok(64).eq.1) then
          pdbxyzs%nc4l(j) = pdbxyzs%nc4(j)
          ok2(5) = 1
        else
          ok2(5) = 0
        end if
        if (ok(66).eq.1) then
          pdbxyzs%no3l(j) = pdbxyzs%no3(j)
          ok2(6) = 1
        else
          ok2(6) = 0
        end if
        if (ok(67).eq.1) then
          pdbxyzs%no32l(j) = pdbxyzs%no32(j)
          ok2(7) = 1
        else
          ok2(7) = 0
        end if
        if (ok(63).eq.1) then
          pdbxyzs%nc5l(j) = pdbxyzs%nc5(j)
          ok2(8) = 1
        else
          ok2(8) = 0
        end if
      else if (seqpolty(rs).eq.'P') then
        if (ok(18).eq.1) then
          pdbxyzs%nl(j) = pdbxyzs%n(j)
          ok2(1) = 1
        else
          ok2(1) = 0
        end if
        if (ok(19).eq.1) then
          pdbxyzs%cal(j) = pdbxyzs%ca(j)
          ok2(2) = 1
        else
          ok2(2) = 0
        end if
        if (ok(20).eq.1) then
          pdbxyzs%cl(j) = pdbxyzs%c(j)
          ok2(3) = 1
        else
          ok2(3) = 0
        end if
        if (ok(27).eq.1) then
          pdbxyzs%ol(j) = pdbxyzs%o(j)
          ok2(9) = 1
        else
          ok2(9) = 0
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported polymer type i&
 &n extracting backbone angles from pdb.'
        call fexit()
      end if
    end do
!
!   eliminate back references for end of molecule
    if (rs.eq.rsmol(molofrs(rs),2)) then
      ok2(:) = 0
    end if
!
    rs = rs + 1
    rs2 = rs2 + 1
!
  end do
!
! finally, we of course have to build the structure from the new values
  call fyczmat()
  if (lastrs.le.0) then
    write(ilog,*) 'Fatal. Failed to read in any structural information for system from pdb-file ',fn(t1:t2),'.'
    call fexit()
  else
    if (allocated(pdb_didread).EQV..true.) deallocate(pdb_didread)
    allocate(pdb_didread(molofrs(lastrs),3))
    do i=1,molofrs(lastrs)
      pdb_didread(i,1) = rsmol(i,1)
      pdb_didread(i,2) = min(lastrs,rsmol(i,2))
    end do
  end if
!
  do i=1,nmol
    call makexyz_formol(i)
    call update_rigid(i)
    call update_rigidm(i)
    call update_image(i)
  end do
!
! deallocation last
  deallocate(pdbl)
  deallocate(pdbnm)
  deallocate(rsnam)
  deallocate(chnam)
  deallocate(kw)
  deallocate(rsnmb)
  deallocate(xr)
  deallocate(yr)
  deallocate(zr)
  deallocate(lix)
  write(ilog,*)
!
end
!
!----------------------------------------------------------------------------------
!
subroutine readpdb1_collect(stat,i,pdbnm,xr,yr,zr,xyzs,ok)
!
  use iounit
  use pdb
  use atoms
  use sequen
  use molecule
!
  implicit none
!
  integer i,stat
  character(4) pdbnm(3*sum(natres(1:nseq))+100)
  integer ok(120)
  RTYPE xr(3*sum(natres(1:nseq))+100),yr(3*sum(natres(1:nseq))+100),zr(3*sum(natres(1:nseq))+100)
  type(t_pdbxyzs) xyzs
!
  if (pdbnm(i).eq.' N  ') then
    xyzs%n(1) = xr(i)
    xyzs%n(2) = yr(i)
    xyzs%n(3) = zr(i)
    ok(18) = 1
  else if (pdbnm(i).eq.' CA ') then
    xyzs%ca(1) = xr(i)
    xyzs%ca(2) = yr(i)
    xyzs%ca(3) = zr(i)
    ok(19) = 1
  else if (pdbnm(i).eq.' C  ') then
    xyzs%c(1) = xr(i)
    xyzs%c(2) = yr(i)
    xyzs%c(3) = zr(i) 
    ok(20) = 1
  else if (pdbnm(i).eq.' O  ') then
    xyzs%o(1) = xr(i)
    xyzs%o(2) = yr(i)
    xyzs%o(3) = zr(i) 
    ok(27) = 1
  else if (pdbnm(i).eq.' H  ') then
    xyzs%h(1) = xr(i)
    xyzs%h(2) = yr(i)
    xyzs%h(3) = zr(i) 
    ok(29) = 1
  else if (pdbnm(i).eq.' CH3') then
    xyzs%ch3(1) = xr(i)
    xyzs%ch3(2) = yr(i)
    xyzs%ch3(3) = zr(i)
    ok(30) = 1
  else if (pdbnm(i).eq.' CB ') then
    xyzs%cb(1) = xr(i)
    xyzs%cb(2) = yr(i)
    xyzs%cb(3) = zr(i)
    ok(1) = 1
  else if (pdbnm(i).eq.' CG ') then
    xyzs%cg(1) = xr(i)
    xyzs%cg(2) = yr(i)
    xyzs%cg(3) = zr(i)
    ok(2) = 1
  else if (pdbnm(i).eq.' CD ') then
    xyzs%cd(1) = xr(i)
    xyzs%cd(2) = yr(i)
    xyzs%cd(3) = zr(i)
    ok(3) = 1
  else if (pdbnm(i).eq.' CE ') then
    xyzs%ce(1) = xr(i)
    xyzs%ce(2) = yr(i)
    xyzs%ce(3) = zr(i)
    ok(4) = 1
  else if (pdbnm(i).eq.' CG ') then
    xyzs%cg(1) = xr(i)
    xyzs%cg(2) = yr(i)
    xyzs%cg(3) = zr(i)
    ok(5) = 1
  else if (pdbnm(i).eq.' CG1') then
    xyzs%cg1(1) = xr(i)
    xyzs%cg1(2) = yr(i)
    xyzs%cg1(3) = zr(i)
    ok(6) = 1
  else if (pdbnm(i).eq.' CD1') then
    xyzs%cd1(1) = xr(i)
    xyzs%cd1(2) = yr(i)
    xyzs%cd1(3) = zr(i)
    ok(7) = 1
  else if (pdbnm(i).eq.' SG ') then
    xyzs%sg(1) = xr(i)
    xyzs%sg(2) = yr(i)
    xyzs%sg(3) = zr(i)
    ok(8) = 1
  else if (pdbnm(i).eq.' OG ') then
    xyzs%og(1) = xr(i)
    xyzs%og(2) = yr(i)
    xyzs%og(3) = zr(i)
    ok(9) = 1
  else if (pdbnm(i).eq.' HG ') then
    xyzs%hg(1) = xr(i)
    xyzs%hg(2) = yr(i)
    xyzs%hg(3) = zr(i)
    ok(10) = 1
  else if (pdbnm(i).eq.' HG1') then
    xyzs%hg1(1) = xr(i)
    xyzs%hg1(2) = yr(i)
    xyzs%hg1(3) = zr(i)
    ok(11) = 1
  else if (pdbnm(i).eq.' SD ') then
    xyzs%sd(1) = xr(i)
    xyzs%sd(2) = yr(i)
    xyzs%sd(3) = zr(i)
    ok(12) = 1
  else if (pdbnm(i).eq.' OD1') then
    xyzs%od1(1) = xr(i)
    xyzs%od1(2) = yr(i)
    xyzs%od1(3) = zr(i)
    ok(13) = 1
  else if (pdbnm(i).eq.' ND1') then
    xyzs%nd1(1) = xr(i)
    xyzs%nd1(2) = yr(i)
    xyzs%nd1(3) = zr(i)
    ok(14) = 1
  else if (pdbnm(i).eq.' CZ ') then
    xyzs%cz(1) = xr(i)
    xyzs%cz(2) = yr(i)
    xyzs%cz(3) = zr(i)
    ok(15) = 1
  else if (pdbnm(i).eq.' NE ') then
    xyzs%ne(1) = xr(i)
    xyzs%ne(2) = yr(i)
    xyzs%ne(3) = zr(i)
    ok(16) = 1
  else if (pdbnm(i).eq.' NZ ') then
    xyzs%nz(1) = xr(i)
    xyzs%nz(2) = yr(i)
    xyzs%nz(3) = zr(i)
    ok(17) = 1
  else if (pdbnm(i).eq.' OE1') then
    xyzs%oe1(1) = xr(i)
    xyzs%oe1(2) = yr(i)
    xyzs%oe1(3) = zr(i)
    ok(21) = 1
  else if (pdbnm(i).eq.' OG1') then
    xyzs%og1(1) = xr(i)
    xyzs%og1(2) = yr(i)
    xyzs%og1(3) = zr(i)
    ok(22) = 1
  else if (pdbnm(i).eq.' OW ') then
    xyzs%ow(1) = xr(i)
    xyzs%ow(2) = yr(i)
    xyzs%ow(3) = zr(i)
    ok(23) = 1
  else if (pdbnm(i).eq.'1HW ') then
    xyzs%hw1(1) = xr(i)
    xyzs%hw1(2) = yr(i)
    xyzs%hw1(3) = zr(i)
    ok(24) = 1
  else if (pdbnm(i).eq.'2HW ') then
    xyzs%hw2(1) = xr(i)
    xyzs%hw2(2) = yr(i)
    xyzs%hw2(3) = zr(i)
    ok(25) = 1
  else if (pdbnm(i).eq.' N1 ') then
    xyzs%n1(1) = xr(i)
    xyzs%n1(2) = yr(i)
    xyzs%n1(3) = zr(i)
    ok(28) = 1
  else if (pdbnm(i).eq.' CT ') then
    xyzs%ct(1) = xr(i)
    xyzs%ct(2) = yr(i)
    xyzs%ct(3) = zr(i)
    ok(31) = 1
  else if (pdbnm(i).eq.'1HT ') then
    xyzs%ht1(1) = xr(i)
    xyzs%ht1(2) = yr(i)
    xyzs%ht1(3) = zr(i)
    ok(32) = 1
  else if (pdbnm(i).eq.'2HT ') then
    xyzs%ht2(1) = xr(i)
    xyzs%ht2(2) = yr(i)
    xyzs%ht2(3) = zr(i)
    ok(33) = 1
  else if (pdbnm(i).eq.' C1 ') then
    xyzs%c1(1) = xr(i)
    xyzs%c1(2) = yr(i)
    xyzs%c1(3) = zr(i)
    ok(34) = 1
  else if (pdbnm(i).eq.' C21') then
    xyzs%c21(1) = xr(i)
    xyzs%c21(2) = yr(i)
    xyzs%c21(3) = zr(i)
    ok(35) = 1
  else if (pdbnm(i).eq.'1HCT') then
    xyzs%hct1(1) = xr(i)
    xyzs%hct1(2) = yr(i)
    xyzs%hct1(3) = zr(i)
    ok(36) = 1
  else if (pdbnm(i).eq.'1HNT') then
    xyzs%hnt1(1) = xr(i)
    xyzs%hnt1(2) = yr(i)
    xyzs%hnt1(3) = zr(i)
    ok(37) = 1
  else if (pdbnm(i).eq.' HO ') then
    xyzs%ho(1) = xr(i)
    xyzs%ho(2) = yr(i)
    xyzs%ho(3) = zr(i)
    ok(38) = 1
  else if (pdbnm(i).eq.' CCT') then
    xyzs%cct(1) = xr(i)
    xyzs%cct(2) = yr(i)
    xyzs%cct(3) = zr(i)
    ok(39) = 1
  else if (pdbnm(i).eq.' CNT') then
    xyzs%cnt(1) = xr(i)
    xyzs%cnt(2) = yr(i)
    xyzs%cnt(3) = zr(i)
    ok(40) = 1
  else if (pdbnm(i).eq.' C31') then
    xyzs%c31(1) = xr(i)
    xyzs%c31(2) = yr(i)
    xyzs%c31(3) = zr(i)
    ok(41) = 1
  else if (pdbnm(i).eq.' C4 ') then
    xyzs%c4(1) = xr(i)
    xyzs%c4(2) = yr(i)
    xyzs%c4(3) = zr(i)
    ok(42) = 1
  else if (pdbnm(i).eq.' CBT') then
    xyzs%cbt(1) = xr(i)
    xyzs%cbt(2) = yr(i)
    xyzs%cbt(3) = zr(i)
    ok(43) = 1
  else if (pdbnm(i).eq.'1HBT') then
    xyzs%hbt1(1) = xr(i)
    xyzs%hbt1(2) = yr(i)
    xyzs%hbt1(3) = zr(i)
    ok(44) = 1
  else if (pdbnm(i).eq.' C2 ') then
    xyzs%c2(1) = xr(i)
    xyzs%c2(2) = yr(i)
    xyzs%c2(3) = zr(i)
    ok(45) = 1
  else if (pdbnm(i).eq.'1H1 ') then
    xyzs%h11(1) = xr(i)
    xyzs%h11(2) = yr(i)
    xyzs%h11(3) = zr(i)
    ok(46) = 1
  else if (pdbnm(i).eq.'1H2 ') then
    xyzs%h21(1) = xr(i)
    xyzs%h21(2) = yr(i)
    xyzs%h21(3) = zr(i)
    ok(47) = 1
  else if (pdbnm(i).eq.'2OXT') then
    xyzs%oxt(1) = xr(i)
    xyzs%oxt(2) = yr(i)
    xyzs%oxt(3) = zr(i)
    ok(48) = 1
  else if (pdbnm(i).eq.' OH ') then
    xyzs%oh(1) = xr(i)
    xyzs%oh(2) = yr(i)
    xyzs%oh(3) = zr(i)
    ok(49) = 1
  else if (pdbnm(i).eq.' HH ') then
    xyzs%hh(1) = xr(i)
    xyzs%hh(2) = yr(i)
    xyzs%hh(3) = zr(i)
    ok(50) = 1
  else if (pdbnm(i).eq.' CE2') then
    xyzs%ce2(1) = xr(i)
    xyzs%ce2(2) = yr(i)
    xyzs%ce2(3) = zr(i)
    ok(51) = 1
  else if (pdbnm(i).eq.'1H  ') then
    xyzs%h1(1) = xr(i)
    xyzs%h1(2) = yr(i)
    xyzs%h1(3) = zr(i)
    ok(52) = 1
  else if (pdbnm(i).eq.' ND ') then
    xyzs%nd(1) = xr(i)
    xyzs%nd(2) = yr(i)
    xyzs%nd(3) = zr(i)
    ok(53) = 1
  else if (pdbnm(i).eq.'1HN ') then
    xyzs%hn1(1) = xr(i)
    xyzs%hn1(2) = yr(i)
    xyzs%hn1(3) = zr(i)
    ok(54) = 1
  else if (pdbnm(i).eq.'2HN ') then
    xyzs%hn2(1) = xr(i)
    xyzs%hn2(2) = yr(i)
    xyzs%hn2(3) = zr(i)
    ok(55) = 1
  else if (pdbnm(i).eq.'1OXT') then
    xyzs%oxt1(1) = xr(i)
    xyzs%oxt1(2) = yr(i)
    xyzs%oxt1(3) = zr(i)
    ok(56) = 1
  else if (pdbnm(i).eq.' N2 ') then
    xyzs%n2(1) = xr(i)
    xyzs%n2(2) = yr(i)
    xyzs%n2(3) = zr(i)
    ok(57) = 1
  else if (pdbnm(i).eq.' C6 ') then
    xyzs%c6(1) = xr(i)
    xyzs%c6(2) = yr(i)
    xyzs%c6(3) = zr(i)
    ok(58) = 1
  else if (pdbnm(i).eq.' N9 ') then
    xyzs%n9(1) = xr(i)
    xyzs%n9(2) = yr(i)
    xyzs%n9(3) = zr(i)
    ok(59) = 1
  else if (pdbnm(i).eq.' C8 ') then
    xyzs%c8(1) = xr(i)
    xyzs%c8(2) = yr(i)
    xyzs%c8(3) = zr(i)
    ok(60) = 1
  else if (pdbnm(i)(2:4).eq.'O3*') then
    if (ok(61).eq.1) then
      if (((xr(i)-xyzs%np(1))**2+(yr(i)-xyzs%np(2))**2+(zr(i)-xyzs%np(3))**2).gt.4.0) then
        xyzs%no32(1) = xr(i)
        xyzs%no32(2) = yr(i)
        xyzs%no32(3) = zr(i)
        ok(67) = 1
      else
        xyzs%no3(1) = xr(i)
        xyzs%no3(2) = yr(i)
        xyzs%no3(3) = zr(i)
        ok(66) = 1
      end if
!   if P not yet defined we most likely have CAMPARI convention and phosphate O3*
!   or 5'-terminal nucleoside
    else
      xyzs%no3(1) = xr(i)
      xyzs%no3(2) = yr(i)
      xyzs%no3(3) = zr(i)
      ok(66) = 1
    end if
  else if (pdbnm(i).eq.' P  ') then
    xyzs%np(1) = xr(i)
    xyzs%np(2) = yr(i)
    xyzs%np(3) = zr(i)
    ok(61) = 1
  else if (pdbnm(i).eq.' O5*') then
    xyzs%no5(1) = xr(i)
    xyzs%no5(2) = yr(i)
    xyzs%no5(3) = zr(i)
    ok(62) = 1
  else if (pdbnm(i).eq.' C5*') then
    xyzs%nc5(1) = xr(i)
    xyzs%nc5(2) = yr(i)
    xyzs%nc5(3) = zr(i)
    ok(63) = 1
  else if (pdbnm(i).eq.' C4*') then
    xyzs%nc4(1) = xr(i)
    xyzs%nc4(2) = yr(i)
    xyzs%nc4(3) = zr(i)
    ok(64) = 1
  else if (pdbnm(i).eq.' C3*') then
    xyzs%nc3(1) = xr(i)
    xyzs%nc3(2) = yr(i)
    xyzs%nc3(3) = zr(i)
    ok(65) = 1
  else if (pdbnm(i).eq.' HOP') then
    xyzs%nhop(1) = xr(i)
    xyzs%nhop(2) = yr(i)
    xyzs%nhop(3) = zr(i)
    ok(68) = 1
  else if (pdbnm(i).eq.'3HO*') then
    xyzs%nho3(1) = xr(i)
    xyzs%nho3(2) = yr(i)
    xyzs%nho3(3) = zr(i)
    ok(69) = 1
  else if (pdbnm(i).eq.' C1*') then
    xyzs%nc1(1) = xr(i)
    xyzs%nc1(2) = yr(i)
    xyzs%nc1(3) = zr(i)
    ok(70) = 1
  else if (pdbnm(i).eq.' C2*') then
    xyzs%nc2(1) = xr(i)
    xyzs%nc2(2) = yr(i)
    xyzs%nc2(3) = zr(i)
    ok(71) = 1
  else if (pdbnm(i).eq.' O2*') then
    xyzs%no2(1) = xr(i)
    xyzs%no2(2) = yr(i)
    xyzs%no2(3) = zr(i)
    ok(72) = 1
  else if (pdbnm(i).eq.'2HO*') then
    xyzs%nho2(1) = xr(i)
    xyzs%nho2(2) = yr(i)
    xyzs%nho2(3) = zr(i)
    ok(73) = 1
  else if (pdbnm(i).eq.' O1*') then
    xyzs%no1s(1) = xr(i)
    xyzs%no1s(2) = yr(i)
    xyzs%no1s(3) = zr(i)
    ok(74) = 1
  else if (pdbnm(i).eq.'1HO*') then
    xyzs%nho1(1) = xr(i)
    xyzs%nho1(2) = yr(i)
    xyzs%nho1(3) = zr(i)
    ok(75) = 1
  else if (pdbnm(i).eq.' O4*') then
    xyzs%no4(1) = xr(i)
    xyzs%no4(2) = yr(i)
    xyzs%no4(3) = zr(i)
    ok(76) = 1
  else if (pdbnm(i).eq.'5HO*') then
    xyzs%nho5(1) = xr(i)
    xyzs%nho5(2) = yr(i)
    xyzs%nho5(3) = zr(i)
    ok(77) = 1
  else if (pdbnm(i).eq.'1H5M') then
    xyzs%thyh1(1) = xr(i)
    xyzs%thyh1(2) = yr(i)
    xyzs%thyh1(3) = zr(i)
    ok(78) = 1
  else if (pdbnm(i).eq.' C5M') then
    xyzs%c5m(1) = xr(i)
    xyzs%c5m(2) = yr(i)
    xyzs%c5m(3) = zr(i)
    ok(79) = 1
  else if (pdbnm(i).eq.' C5 ') then
    xyzs%c5(1) = xr(i)
    xyzs%c5(2) = yr(i)
    xyzs%c5(3) = zr(i)
    ok(80) = 1
  else if (pdbnm(i).eq.' S  ') then
    xyzs%s(1) = xr(i)
    xyzs%s(2) = yr(i)
    xyzs%s(3) = zr(i)
    ok(81) = 1
  else if (pdbnm(i).eq.' C3 ') then
    xyzs%c3(1) = xr(i)
    xyzs%c3(2) = yr(i)
    xyzs%c3(3) = zr(i)
    ok(82) = 1
  else if (pdbnm(i).eq.' CT1') then
    xyzs%ct1(1) = xr(i)
    xyzs%ct1(2) = yr(i)
    xyzs%ct1(3) = zr(i)
    ok(83) = 1
  else if (pdbnm(i).eq.' CT2') then
    xyzs%ct2(1) = xr(i)
    xyzs%ct2(2) = yr(i)
    xyzs%ct2(3) = zr(i)
    ok(84) = 1
  else if (pdbnm(i).eq.' CB1') then
    xyzs%cb1(1) = xr(i)
    xyzs%cb1(2) = yr(i)
    xyzs%cb1(3) = zr(i)
    ok(85) = 1
  else if (pdbnm(i).eq.' CB2') then
    xyzs%cb2(1) = xr(i)
    xyzs%cb2(2) = yr(i)
    xyzs%cb2(3) = zr(i)
    ok(86) = 1
  else if (pdbnm(i).eq.' HS ') then
    xyzs%hs(1) = xr(i)
    xyzs%hs(2) = yr(i)
    xyzs%hs(3) = zr(i)
    ok(87) = 1
  else if (pdbnm(i).eq.'1HT1') then
    xyzs%ht11(1) = xr(i)
    xyzs%ht11(2) = yr(i)
    xyzs%ht11(3) = zr(i)
    ok(88) = 1
  else if (pdbnm(i).eq.'1HT2') then
    xyzs%ht12(1) = xr(i)
    xyzs%ht12(2) = yr(i)
    xyzs%ht12(3) = zr(i)
    ok(89) = 1
  else if (pdbnm(i).eq.'1HT3') then
    xyzs%ht13(1) = xr(i)
    xyzs%ht13(2) = yr(i)
    xyzs%ht13(3) = zr(i)
    ok(90) = 1
  else if (pdbnm(i).eq.' CT3') then
    xyzs%ct3(1) = xr(i)
    xyzs%ct3(2) = yr(i)
    xyzs%ct3(3) = zr(i)
    ok(91) = 1
  else if (pdbnm(i).eq.' O1 ') then
    xyzs%o1(1) = xr(i)
    xyzs%o1(2) = yr(i)
    xyzs%o1(3) = zr(i)
    ok(92) = 1
  else if (pdbnm(i).eq.' O2 ') then
    xyzs%o2(1) = xr(i)
    xyzs%o2(2) = yr(i)
    xyzs%o2(3) = zr(i)
    ok(93) = 1
  else if (pdbnm(i).eq.' CL ') then
    xyzs%chl(1) = xr(i)
    xyzs%chl(2) = yr(i)
    xyzs%chl(3) = zr(i)
    ok(94) = 1
  else if (pdbnm(i).eq.' C11') then
    xyzs%c11(1) = xr(i)
    xyzs%c11(2) = yr(i)
    xyzs%c11(3) = zr(i)
    ok(95) = 1
  else if (pdbnm(i).eq.' C12') then
    xyzs%c12(1) = xr(i)
    xyzs%c12(2) = yr(i)
    xyzs%c12(3) = zr(i)
    ok(96) = 1
  else if (pdbnm(i).eq.' OE2') then
    xyzs%oe2(1) = xr(i)
    xyzs%oe2(2) = yr(i)
    xyzs%oe2(3) = zr(i)
    ok(97) = 1
  else if (pdbnm(i).eq.' OD2') then
    xyzs%od2(1) = xr(i)
    xyzs%od2(2) = yr(i)
    xyzs%od2(3) = zr(i)
    ok(98) = 1
  else if (pdbnm(i).eq.' HE2') then
    xyzs%he2(1) = xr(i)
    xyzs%he2(2) = yr(i)
    xyzs%he2(3) = zr(i)
    ok(99) = 1
  else if (pdbnm(i).eq.' HD2') then
    xyzs%hd2(1) = xr(i)
    xyzs%hd2(2) = yr(i)
    xyzs%hd2(3) = zr(i)
    ok(100) = 1
  else if (pdbnm(i).eq.' OZ ') then
    xyzs%oz(1) = xr(i)
    xyzs%oz(2) = yr(i)
    xyzs%oz(3) = zr(i)
    ok(101) = 1
  else if (pdbnm(i).eq.'1HZ ') then
    xyzs%hnz1(1) = xr(i)
    xyzs%hnz1(2) = yr(i)
    xyzs%hnz1(3) = zr(i)
    ok(102) = 1
  else if (pdbnm(i).eq.' CH1') then
    xyzs%ch1(1) = xr(i)
    xyzs%ch1(2) = yr(i)
    xyzs%ch1(3) = zr(i)
    ok(103) = 1
  else if (pdbnm(i).eq.' CH ') then
    xyzs%ch(1) = xr(i)
    xyzs%ch(2) = yr(i)
    xyzs%ch(3) = zr(i)
    ok(104) = 1
  else if (pdbnm(i).eq.' CI ') then
    xyzs%ci(1) = xr(i)
    xyzs%ci(2) = yr(i)
    xyzs%ci(3) = zr(i)
    ok(105) = 1
  else if (stat.eq.i) then
    xyzs%no1(1) = xr(i)
    xyzs%no1(2) = yr(i)
    xyzs%no1(3) = zr(i)
    ok(26) = 1
  end if
!
end
!
!----------------------------------------------------------------------------------
!
! this one transfers over rigid-body coords
! note we allow some exceptions for incompleteness for OH2, CH4, and NH4
!
subroutine readpdb1_rbs(rs,resname,xyzs,ok)
!
  use iounit
  use pdb
  use atoms
  use sequen
  use molecule
  use polypep
  use aminos
  use system
!
  implicit none
!
  character(3) resname
  integer rs
  integer ok(120)
  type(t_pdbxyzs) xyzs
!
  if ((resname.eq.'NA+').OR.&
 &      (resname.eq.'K+ ').OR.(resname.eq.'BR-').OR.&
 &      (resname.eq.'CS+').OR.(resname.eq.'I- ')) then
    if (ok(26).ne.1) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%no1(1)
    y(at(rs)%bb(1)) = xyzs%no1(2)
    z(at(rs)%bb(1)) = xyzs%no1(3)
  else if (resname.eq.'CL-') then
    if (ok(94).ne.1) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%chl(1)
    y(at(rs)%bb(1)) = xyzs%chl(2)
    z(at(rs)%bb(1)) = xyzs%chl(3)
  else if (resname.eq.'O2 ') then
    if ((ok(92).ne.1).OR.(ok(93).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%o1(1)
    y(at(rs)%bb(1)) = xyzs%o1(2)
    z(at(rs)%bb(1)) = xyzs%o1(3)
    x(at(rs)%bb(2)) = xyzs%o2(1)
    y(at(rs)%bb(2)) = xyzs%o2(2)
    z(at(rs)%bb(2)) = xyzs%o2(3)
  else if ((resname.eq.'SPC').OR.(resname.eq.'T3P').OR.&
 &         (resname.eq.'T4P').OR.(resname.eq.'T4E').OR.&
 &         (resname.eq.'T5P')) then
    if ((ok(23).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ow(1)
    y(at(rs)%bb(1)) = xyzs%ow(2)
    z(at(rs)%bb(1)) = xyzs%ow(3)
    if ((ok(24).ne.1).OR.(ok(25).ne.1)) then
      call makexyz_formol3(molofrs(rs)) ! create a proper reference frame
    else
      x(at(rs)%bb(2)) = xyzs%hw1(1)
      y(at(rs)%bb(2)) = xyzs%hw1(2)
      z(at(rs)%bb(2)) = xyzs%hw1(3)
      x(at(rs)%bb(3)) = xyzs%hw2(1)
      y(at(rs)%bb(3)) = xyzs%hw2(2)
      z(at(rs)%bb(3)) = xyzs%hw2(3)
    end if
  else if (resname.eq.'NH4') then
    if (ok(18).ne.1) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%n(1)
    y(at(rs)%bb(1)) = xyzs%n(2)
    z(at(rs)%bb(1)) = xyzs%n(3)
    if ((ok(54).ne.1).OR.(ok(55).ne.1)) then
      call makexyz_formol3(molofrs(rs)) ! create a proper reference frame
    else
      x(at(rs)%bb(2)) = xyzs%hn1(1)
      y(at(rs)%bb(2)) = xyzs%hn1(2)
      z(at(rs)%bb(2)) = xyzs%hn1(3)
      x(at(rs)%bb(3)) = xyzs%hn2(1)
      y(at(rs)%bb(3)) = xyzs%hn2(2)
      z(at(rs)%bb(3)) = xyzs%hn2(3)
    end if
  else if (resname.eq.'NO3') then
    if ((ok(18).ne.1).OR.(ok(92).ne.1).OR.(ok(93).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%n(1)
    y(at(rs)%bb(1)) = xyzs%n(2)
    z(at(rs)%bb(1)) = xyzs%n(3)
    x(at(rs)%bb(2)) = xyzs%o1(1)
    y(at(rs)%bb(2)) = xyzs%o1(2)
    z(at(rs)%bb(2)) = xyzs%o1(3)
    x(at(rs)%bb(3)) = xyzs%o2(1)
    y(at(rs)%bb(3)) = xyzs%o2(2)
    z(at(rs)%bb(3)) = xyzs%o2(3)
  else if (resname.eq.'LCP') then
    if ((ok(94).ne.1).OR.(ok(92).ne.1).OR.(ok(93).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%chl(1)
    y(at(rs)%bb(1)) = xyzs%chl(2)
    z(at(rs)%bb(1)) = xyzs%chl(3)
    x(at(rs)%bb(2)) = xyzs%o1(1)
    y(at(rs)%bb(2)) = xyzs%o1(2)
    z(at(rs)%bb(2)) = xyzs%o1(3)
    x(at(rs)%bb(3)) = xyzs%o2(1)
    y(at(rs)%bb(3)) = xyzs%o2(2)
    z(at(rs)%bb(3)) = xyzs%o2(3)
  else if (resname.eq.'AC-') then
    if ((ok(20).ne.1).OR.(ok(48).ne.1).OR.(ok(56).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%c(1)
    y(at(rs)%bb(1)) = xyzs%c(2)
    z(at(rs)%bb(1)) = xyzs%c(3)
    x(at(rs)%bb(2)) = xyzs%oxt1(1)
    y(at(rs)%bb(2)) = xyzs%oxt1(2)
    z(at(rs)%bb(2)) = xyzs%oxt1(3)
    x(at(rs)%bb(3)) = xyzs%oxt(1)
    y(at(rs)%bb(3)) = xyzs%oxt(2)
    z(at(rs)%bb(3)) = xyzs%oxt(3)

  else if (resname.eq.'AC-') then
    if ((ok(20).ne.1).OR.(ok(48).ne.1).OR.(ok(56).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%c(1)
    y(at(rs)%bb(1)) = xyzs%c(2)
    z(at(rs)%bb(1)) = xyzs%c(3)
    x(at(rs)%bb(2)) = xyzs%oxt1(1)
    y(at(rs)%bb(2)) = xyzs%oxt1(2)
    z(at(rs)%bb(2)) = xyzs%oxt1(3)
    x(at(rs)%bb(3)) = xyzs%oxt(1)
    y(at(rs)%bb(3)) = xyzs%oxt(2)
    z(at(rs)%bb(3)) = xyzs%oxt(3)
  else if (resname.eq.'GDN') then
    if ((ok(20).ne.1).OR.(ok(28).ne.1).OR.(ok(57).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%c(1)
    y(at(rs)%bb(1)) = xyzs%c(2)
    z(at(rs)%bb(1)) = xyzs%c(3)
    x(at(rs)%bb(2)) = xyzs%n1(1)
    y(at(rs)%bb(2)) = xyzs%n1(2)
    z(at(rs)%bb(2)) = xyzs%n1(3)
    x(at(rs)%bb(3)) = xyzs%n2(1)
    y(at(rs)%bb(3)) = xyzs%n2(2)
    z(at(rs)%bb(3)) = xyzs%n2(3)
  else if (resname.eq.'1MN') then
    if ((ok(18).ne.1).OR.(ok(31).ne.1).OR.(ok(54).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%n(1)
    y(at(rs)%bb(1)) = xyzs%n(2)
    z(at(rs)%bb(1)) = xyzs%n(3)
    x(at(rs)%bb(2)) = xyzs%ct(1)
    y(at(rs)%bb(2)) = xyzs%ct(2)
    z(at(rs)%bb(2)) = xyzs%ct(3)
    x(at(rs)%bb(3)) = xyzs%hn1(1)
    y(at(rs)%bb(3)) = xyzs%hn1(2)
    z(at(rs)%bb(3)) = xyzs%hn1(3)
  else if (resname.eq.'2MN') then
    if ((ok(18).ne.1).OR.(ok(83).ne.1).OR.(ok(84).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct1(1)
    y(at(rs)%bb(1)) = xyzs%ct1(2)
    z(at(rs)%bb(1)) = xyzs%ct1(3)
    x(at(rs)%bb(2)) = xyzs%n(1)
    y(at(rs)%bb(2)) = xyzs%n(2)
    z(at(rs)%bb(2)) = xyzs%n(3)
    x(at(rs)%bb(3)) = xyzs%ct2(1)
    y(at(rs)%bb(3)) = xyzs%ct2(2)
    z(at(rs)%bb(3)) = xyzs%ct2(3)
  else if (resname.eq.'MOH') then
    if (ua_model.eq.0) then
      if ((ok(31).ne.1).OR.(ok(27).ne.1).OR.(ok(32).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
    else
      if ((ok(31).ne.1).OR.(ok(27).ne.1).OR.(ok(38).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
    end if
    x(at(rs)%bb(1)) = xyzs%ct(1)
    y(at(rs)%bb(1)) = xyzs%ct(2)
    z(at(rs)%bb(1)) = xyzs%ct(3)
    x(at(rs)%sc(1)) = xyzs%o(1)
    y(at(rs)%sc(1)) = xyzs%o(2)
    z(at(rs)%sc(1)) = xyzs%o(3)
    if (ua_model.eq.0) then
      x(at(rs)%bb(2)) = xyzs%ht1(1)
      y(at(rs)%bb(2)) = xyzs%ht1(2)
      z(at(rs)%bb(2)) = xyzs%ht1(3)
    else
      x(at(rs)%sc(2)) = xyzs%ho(1)
      y(at(rs)%sc(2)) = xyzs%ho(2)
      z(at(rs)%sc(2)) = xyzs%ho(3)
    end if
  else if (resname.eq.'EOH') then
    if ((ok(31).ne.1).OR.(ok(1).ne.1).OR.(ok(27).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct(1)
    y(at(rs)%bb(1)) = xyzs%ct(2)
    z(at(rs)%bb(1)) = xyzs%ct(3)
    x(at(rs)%bb(2)) = xyzs%cb(1)
    y(at(rs)%bb(2)) = xyzs%cb(2)
    z(at(rs)%bb(2)) = xyzs%cb(3)
    x(at(rs)%bb(3)) = xyzs%o(1)
    y(at(rs)%bb(3)) = xyzs%o(2)
    z(at(rs)%bb(3)) = xyzs%o(3)
  else if (resname.eq.'MSH') then
    if (ua_model.eq.0) then
      if ((ok(31).ne.1).OR.(ok(81).ne.1).OR.(ok(32).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
    else if (ua_model.eq.1) then
      if ((ok(31).ne.1).OR.(ok(81).ne.1).OR.(ok(87).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
    else
      if ((ok(31).ne.1).OR.(ok(81).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
    end if
    x(at(rs)%bb(1)) = xyzs%ct(1)
    y(at(rs)%bb(1)) = xyzs%ct(2)
    z(at(rs)%bb(1)) = xyzs%ct(3)
    x(at(rs)%sc(1)) = xyzs%s(1)
    y(at(rs)%sc(1)) = xyzs%s(2)
    z(at(rs)%sc(1)) = xyzs%s(3)
    if (ua_model.eq.0) then
      x(at(rs)%bb(2)) = xyzs%ht1(1)
      y(at(rs)%bb(2)) = xyzs%ht1(2)
      z(at(rs)%bb(2)) = xyzs%ht1(3)
    else if (ua_model.eq.1) then
      x(at(rs)%sc(2)) = xyzs%hs(1)
      y(at(rs)%sc(2)) = xyzs%hs(2)
      z(at(rs)%sc(2)) = xyzs%hs(3)
    end if
  else if ((resname.eq.'PPA').OR.(resname.eq.'ACA').OR.&
 &         (resname.eq.'NMF').OR.(resname.eq.'NMA').OR.&
 &         (resname.eq.'FOA').OR.(resname.eq.'DMA')) then
    if ((ok(20).ne.1).OR.(ok(27).ne.1).OR.(ok(18).ne.1)) then
    end if
    x(at(rs)%bb(1)) = xyzs%c(1)
    y(at(rs)%bb(1)) = xyzs%c(2)
    z(at(rs)%bb(1)) = xyzs%c(3)
    x(at(rs)%bb(2)) = xyzs%o(1)
    y(at(rs)%bb(2)) = xyzs%o(2)
    z(at(rs)%bb(2)) = xyzs%o(3)
    x(at(rs)%bb(3)) = xyzs%n(1)
    y(at(rs)%bb(3)) = xyzs%n(2)
    z(at(rs)%bb(3)) = xyzs%n(3)
  else if (resname.eq.'URE') then
    if ((ok(20).ne.1).OR.(ok(27).ne.1).OR.(ok(28).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%c(1)
    y(at(rs)%bb(1)) = xyzs%c(2)
    z(at(rs)%bb(1)) = xyzs%c(3)
    x(at(rs)%bb(2)) = xyzs%o(1)
    y(at(rs)%bb(2)) = xyzs%o(2)
    z(at(rs)%bb(2)) = xyzs%o(3)
    x(at(rs)%bb(3)) = xyzs%n1(1)
    y(at(rs)%bb(3)) = xyzs%n1(2)
    z(at(rs)%bb(3)) = xyzs%n1(3)
   else if (resname.eq.'CH4') then
    if (ua_model.eq.0) then
      if (ok(31).ne.1) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(at(rs)%bb(1)) = xyzs%ct(1)
      y(at(rs)%bb(1)) = xyzs%ct(2)
      z(at(rs)%bb(1)) = xyzs%ct(3)
      if ((ok(32).ne.1).OR.(ok(33).ne.1)) then
        call makexyz_formol3(molofrs(rs)) ! create a proper reference frame
      else
        x(at(rs)%bb(2)) = xyzs%ht1(1)
        y(at(rs)%bb(2)) = xyzs%ht1(2)
        z(at(rs)%bb(2)) = xyzs%ht1(3)
        x(at(rs)%bb(3)) = xyzs%ht2(1)
        y(at(rs)%bb(3)) = xyzs%ht2(2)
        z(at(rs)%bb(3)) = xyzs%ht2(3)
      end if
    else
      if (ok(31).ne.1) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(at(rs)%bb(1)) = xyzs%ct(1)
      y(at(rs)%bb(1)) = xyzs%ct(2)
      z(at(rs)%bb(1)) = xyzs%ct(3)
    end if
  else if ((resname.eq.'PRP').OR.(resname.eq.'IBU')) then
    if ((ok(83).ne.1).OR.(ok(1).ne.1).OR.(ok(84).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct1(1)
    y(at(rs)%bb(1)) = xyzs%ct1(2)
    z(at(rs)%bb(1)) = xyzs%ct1(3)
    x(at(rs)%bb(2)) = xyzs%cb(1)
    y(at(rs)%bb(2)) = xyzs%cb(2)
    z(at(rs)%bb(2)) = xyzs%cb(3)
    x(at(rs)%bb(3)) = xyzs%ct2(1)
    y(at(rs)%bb(3)) = xyzs%ct2(2)
    z(at(rs)%bb(3)) = xyzs%ct2(3)
  else if (resname.eq.'NBU') then
    if ((ok(83).ne.1).OR.(ok(85).ne.1).OR.(ok(86).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct1(1)
    y(at(rs)%bb(1)) = xyzs%ct1(2)
    z(at(rs)%bb(1)) = xyzs%ct1(3)
    x(at(rs)%bb(2)) = xyzs%cb1(1)
    y(at(rs)%bb(2)) = xyzs%cb1(2)
    z(at(rs)%bb(2)) = xyzs%cb1(3)
    x(at(rs)%bb(3)) = xyzs%cb2(1)
    y(at(rs)%bb(3)) = xyzs%cb2(2)
    z(at(rs)%bb(3)) = xyzs%cb2(3)
  else if (resname.eq.'EMT') then
    if ((ok(40).ne.1).OR.(ok(81).ne.1).OR.(ok(1).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%cnt(1)
    y(at(rs)%bb(1)) = xyzs%cnt(2)
    z(at(rs)%bb(1)) = xyzs%cnt(3)
    x(at(rs)%bb(2)) = xyzs%s(1)
    y(at(rs)%bb(2)) = xyzs%s(2)
    z(at(rs)%bb(2)) = xyzs%s(3)
    x(at(rs)%bb(3)) = xyzs%cb(1)
    y(at(rs)%bb(3)) = xyzs%cb(2)
    z(at(rs)%bb(3)) = xyzs%cb(3)
  else if ((resname.eq.'PCR').OR.(resname.eq.'TOL')) then
    if ((ok(31).ne.1).OR.(ok(34).ne.1).OR.(ok(35).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct(1)
    y(at(rs)%bb(1)) = xyzs%ct(2)
    z(at(rs)%bb(1)) = xyzs%ct(3)
    x(at(rs)%bb(2)) = xyzs%c1(1)
    y(at(rs)%bb(2)) = xyzs%c1(2)
    z(at(rs)%bb(2)) = xyzs%c1(3)
    x(at(rs)%bb(3)) = xyzs%c21(1)
    y(at(rs)%bb(3)) = xyzs%c21(2)
    z(at(rs)%bb(3)) = xyzs%c21(3)
  else if (resname.eq.'BEN') then
    if ((ok(82).ne.1).OR.(ok(34).ne.1).OR.(ok(45).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%c1(1)
    y(at(rs)%bb(1)) = xyzs%c1(2)
    z(at(rs)%bb(1)) = xyzs%c1(3)
    x(at(rs)%bb(2)) = xyzs%c2(1)
    y(at(rs)%bb(2)) = xyzs%c2(2)
    z(at(rs)%bb(2)) = xyzs%c2(3)
    x(at(rs)%bb(3)) = xyzs%c3(1)
    y(at(rs)%bb(3)) = xyzs%c3(2)
    z(at(rs)%bb(3)) = xyzs%c3(3)
  else if (resname.eq.'NAP') then
    if ((ok(95).ne.1).OR.(ok(96).ne.1).OR.(ok(35).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%c11(1)
    y(at(rs)%bb(1)) = xyzs%c11(2)
    z(at(rs)%bb(1)) = xyzs%c11(3)
    x(at(rs)%bb(2)) = xyzs%c12(1)
    y(at(rs)%bb(2)) = xyzs%c12(2)
    z(at(rs)%bb(2)) = xyzs%c12(3)
    x(at(rs)%bb(3)) = xyzs%c21(1)
    y(at(rs)%bb(3)) = xyzs%c21(2)
    z(at(rs)%bb(3)) = xyzs%c21(3)
  else if ((resname.eq.'IME').OR.(resname.eq.'IMD')) then
    if ((ok(31).ne.1).OR.(ok(34).ne.1).OR.(ok(57).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct(1)
    y(at(rs)%bb(1)) = xyzs%ct(2)
    z(at(rs)%bb(1)) = xyzs%ct(3)
    x(at(rs)%bb(2)) = xyzs%c1(1)
    y(at(rs)%bb(2)) = xyzs%c1(2)
    z(at(rs)%bb(2)) = xyzs%c1(3)
    x(at(rs)%bb(3)) = xyzs%n2(1)
    y(at(rs)%bb(3)) = xyzs%n2(2)
    z(at(rs)%bb(3)) = xyzs%n2(3)
  else if (resname.eq.'MIN') then
    if ((ok(31).ne.1).OR.(ok(82).ne.1).OR.(ok(45).ne.1)) then
      call readpdb_bailout(rs,resname,ok)
      return
    end if
    x(at(rs)%bb(1)) = xyzs%ct(1)
    y(at(rs)%bb(1)) = xyzs%ct(2)
    z(at(rs)%bb(1)) = xyzs%ct(3)
    x(at(rs)%bb(2)) = xyzs%c3(1)
    y(at(rs)%bb(2)) = xyzs%c3(2)
    z(at(rs)%bb(2)) = xyzs%c3(3)
    x(at(rs)%bb(3)) = xyzs%c2(1)
    y(at(rs)%bb(3)) = xyzs%c2(2)
    z(at(rs)%bb(3)) = xyzs%c2(3)
  else if (rsmol(molofrs(rs),1).eq.rs) then
    if (resname.eq.'FOR') then
      if (ua_model.eq.0) then
        if ((ok(20).ne.1).OR.(ok(27).ne.1).OR.(ok(29).ne.1)) then
          call readpdb_bailout(rs,resname,ok)
          return
        end if
      else
        if ((ok(20).ne.1).OR.(ok(27).ne.1)) then
          call readpdb_bailout(rs,resname,ok)
          return
        end if
      end if
      x(at(rs)%bb(1)) = xyzs%c(1)
      y(at(rs)%bb(1)) = xyzs%c(2)
      z(at(rs)%bb(1)) = xyzs%c(3)
      x(at(rs)%bb(2)) = xyzs%o(1)
      y(at(rs)%bb(2)) = xyzs%o(2)
      z(at(rs)%bb(2)) = xyzs%o(3)
      if (ua_model.eq.0) then
        x(at(rs)%bb(3)) = xyzs%h(1)
        y(at(rs)%bb(3)) = xyzs%h(2)
        z(at(rs)%bb(3)) = xyzs%h(3)
      end if
    else if (resname.eq.'ACE') then
      if ((ok(20).ne.1).OR.(ok(27).ne.1).OR.(ok(30).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(at(rs)%bb(1)) = xyzs%ch3(1)
      y(at(rs)%bb(1)) = xyzs%ch3(2)
      z(at(rs)%bb(1)) = xyzs%ch3(3)
      x(at(rs)%bb(2)) = xyzs%c(1)
      y(at(rs)%bb(2)) = xyzs%c(2)
      z(at(rs)%bb(2)) = xyzs%c(3)
      x(at(rs)%sc(1)) = xyzs%o(1)
      y(at(rs)%sc(1)) = xyzs%o(2)
      z(at(rs)%sc(1)) = xyzs%o(3)
    else if ((resname.eq.'URA').OR.(resname.eq.'CYT').OR.&
 &           (resname.eq.'THY')) then
      if ((ok(28).ne.1).OR.(ok(45).ne.1).OR.(ok(58).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(at(rs)%bb(1)) = xyzs%n1(1)
      y(at(rs)%bb(1)) = xyzs%n1(2)
      z(at(rs)%bb(1)) = xyzs%n1(3)
      x(at(rs)%bb(2)) = xyzs%c2(1)
      y(at(rs)%bb(2)) = xyzs%c2(2)
      z(at(rs)%bb(2)) = xyzs%c2(3)
      x(at(rs)%bb(3)) = xyzs%c6(1)
      y(at(rs)%bb(3)) = xyzs%c6(2)
      z(at(rs)%bb(3)) = xyzs%c6(3)
    else if ((resname.eq.'ADE').OR.(resname.eq.'GUA').OR.(resname.eq.'PUR')) then
      if ((ok(59).ne.1).OR.(ok(42).ne.1).OR.(ok(60).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(at(rs)%bb(1)) = xyzs%n9(1)
      y(at(rs)%bb(1)) = xyzs%n9(2)
      z(at(rs)%bb(1)) = xyzs%n9(3)
      x(at(rs)%bb(2)) = xyzs%c4(1)
      y(at(rs)%bb(2)) = xyzs%c4(2)
      z(at(rs)%bb(2)) = xyzs%c4(3)
      x(at(rs)%bb(3)) = xyzs%c8(1)
      y(at(rs)%bb(3)) = xyzs%c8(2)
      z(at(rs)%bb(3)) = xyzs%c8(3)
    else if ((resname.eq.'ALA').OR.(resname.eq.'AIB').OR.&
 &           (resname.eq.'GLY').OR.(resname.eq.'ABA').OR.&
 &           (resname.eq.'NVA').OR.(resname.eq.'PRO').OR.&
 &           (resname.eq.'NLE').OR.(resname.eq.'LEU').OR.&
 &           (resname.eq.'ILE').OR.(resname.eq.'VAL').OR.&
 &           (resname.eq.'MET').OR.(resname.eq.'CYS').OR.&
 &           (resname.eq.'THR').OR.(resname.eq.'HIE').OR.&
 &           (resname.eq.'HID').OR.(resname.eq.'TRP').OR.&
 &           (resname.eq.'TYR').OR.(resname.eq.'PHE').OR.&
 &           (resname.eq.'SER').OR.(resname.eq.'ASN').OR.&
 &           (resname.eq.'ASP').OR.(resname.eq.'GLN').OR.&
 &           (resname.eq.'GLU').OR.(resname.eq.'ARG').OR.&
 &           (resname.eq.'LYS').OR.(resname.eq.'ORN').OR.&
 &           (resname.eq.'LYD').OR.(resname.eq.'ASH').OR.&
 &           (resname.eq.'SEP').OR.(resname.eq.'GLH').OR.&
 &           (resname.eq.'TPO').OR.(resname.eq.'CYX').OR.&
 &           (resname.eq.'PTR').OR.(resname.eq.'TYO').OR.&
 &           (resname.eq.'DAB').OR.(resname.eq.'HIP').OR.&
 &           (resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &           (resname.eq.'KM2').OR.(resname.eq.'KM3')) then
      if ((ok(18).ne.1).OR.(ok(19).ne.1).OR.(ok(20).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(at(rs)%bb(1)) = xyzs%n(1)
      y(at(rs)%bb(1)) = xyzs%n(2)
      z(at(rs)%bb(1)) = xyzs%n(3)
      x(at(rs)%bb(2)) = xyzs%ca(1)
      y(at(rs)%bb(2)) = xyzs%ca(2)
      z(at(rs)%bb(2)) = xyzs%ca(3)
      x(at(rs)%bb(3)) = xyzs%c(1)
      y(at(rs)%bb(3)) = xyzs%c(2)
      z(at(rs)%bb(3)) = xyzs%c(3)
    else if ((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)) then
      if ((ok(62).ne.1).OR.(ok(63).ne.1).OR.(ok(64).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(nuci(rs,1)) = xyzs%no5(1)
      y(nuci(rs,1)) = xyzs%no5(2)
      z(nuci(rs,1)) = xyzs%no5(3)
      x(nuci(rs,2)) = xyzs%nc5(1)
      y(nuci(rs,2)) = xyzs%nc5(2)
      z(nuci(rs,2)) = xyzs%nc5(3)
      x(nuci(rs,3)) = xyzs%nc4(1)
      y(nuci(rs,3)) = xyzs%nc4(2)
      z(nuci(rs,3)) = xyzs%nc4(3)
    else if (seqpolty(rs).eq.'N') then
      if ((ok(61).ne.1).OR.(ok(62).ne.1).OR.(ok(66).ne.1)) then
        call readpdb_bailout(rs,resname,ok)
        return
      end if
      x(nuci(rs,1)) = xyzs%no3(1)
      y(nuci(rs,1)) = xyzs%no3(2)
      z(nuci(rs,1)) = xyzs%no3(3)
      x(nuci(rs,2)) = xyzs%np(1)
      y(nuci(rs,2)) = xyzs%np(2)
      z(nuci(rs,2)) = xyzs%np(3)
      x(nuci(rs,3)) = xyzs%no5(1)
      y(nuci(rs,3)) = xyzs%no5(2)
      z(nuci(rs,3)) = xyzs%no5(3)
    else
      write(ilog,*) 'Fatal. Unknown small molecule or N-terminal residue ',rs,' ('&
 &,resname,') in readpdb1_rbs(...).'
      call fexit() 
    end if
  end if
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine readpdb_bailout(rs,resname,ok)
!
  use iounit
  use sequen
!
  implicit none
!
  integer ok(120),rs
  character(3) resname
!
  if (globrandomize.ne.2) then
    write(ilog,*) 'Warning. Cannot read-in coordinates for residue ',rs,' (',resname,'). If possible, &
 &coordinate information for ALL remaining residues will be generated via setting for FMCSC_RANDOMIZE.'
  end if
!
  ok(:) = -1
!
end
!
!----------------------------------------------------------------------------------
!
! this one extracts chi-angles
!
subroutine readpdb1_chis(rs,resname,xyzs,ok,ok2)
!
  use iounit
  use pdb
  use atoms
  use sequen
  use molecule
  use polypep
  use aminos
  use fyoc
  use math
  use system
  use zmatrix
!
  implicit none
!
  character(3) resname
  integer rs,shf
  integer ok(120),ok2(120)
  RTYPE tor
  type(t_pdbxyzs) xyzs
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  if (resname.eq.'PRO') then
!   we edit the Z-matrix directly (phi is set below)
    if (rs.gt.rsmol(molofrs(rs),1)) then
      if (ok(18)+ok(19)+ok(1)+ok(2)+ok(3)+ok2(3).eq.6) then
        call dihed(xyzs%cl,xyzs%n,xyzs%ca,xyzs%cb,tor)
        ztor(at(rs)%sc(2-shf)) = radian*tor
        call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
        ztor(at(rs)%sc(3-shf)) = radian*tor
        call dihed(xyzs%cb,xyzs%ca,xyzs%n,xyzs%cd,tor)
        ztor(at(rs)%sc(4-shf)) = radian*tor
        call bondang2(xyzs%n,xyzs%ca,xyzs%cb,tor)
        bang(at(rs)%sc(2-shf)) = radian*tor
        call bondang2(xyzs%ca,xyzs%cb,xyzs%cg,tor)
        bang(at(rs)%sc(3-shf)) = radian*tor
        call bondang2(xyzs%ca,xyzs%n,xyzs%cd,tor)
        bang(at(rs)%sc(4-shf)) = radian*tor
      end if
    else
      if (ok(18)+ok(19)+ok(1)+ok(2)+ok(3)+ok(52).eq.6) then
        call dihed(xyzs%h1,xyzs%n,xyzs%ca,xyzs%cb,tor)
        ztor(at(rs)%sc(2-shf)) = radian*tor
        call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
        ztor(at(rs)%sc(3-shf)) = radian*tor
        call dihed(xyzs%cb,xyzs%ca,xyzs%n,xyzs%cd,tor)
        ztor(at(rs)%sc(4-shf)) = radian*tor
        call bondang2(xyzs%n,xyzs%ca,xyzs%cb,tor)
        bang(at(rs)%sc(2-shf)) = radian*tor
        call bondang2(xyzs%ca,xyzs%cb,xyzs%cg,tor)
        bang(at(rs)%sc(3-shf)) = radian*tor
        call bondang2(xyzs%ca,xyzs%n,xyzs%cd,tor)
        bang(at(rs)%sc(4-shf)) = radian*tor
      end if
    end if
  else if (nchi(rs).eq.0) then
!   do nothing 
  else if (resname.eq.'PRP') then
    if (ua_model.eq.0) then
      if (ok(84)+ok(1)+ok(83)+ok(88).eq.4) then
        call dihed(xyzs%ct2,xyzs%cb,xyzs%ct1,xyzs%ht11,tor)
        chi(1,rs) = radian*tor
      end if
      if (ok(84)+ok(1)+ok(83)+ok(89).eq.4) then
        call dihed(xyzs%ct1,xyzs%cb,xyzs%ct2,xyzs%ht12,tor)
        chi(2,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'1MN') then
    if (ua_model.eq.0) then
      if (ok(18)+ok(32)+ok(31)+ok(54).eq.4) then
        call dihed(xyzs%hn1,xyzs%n,xyzs%ct,xyzs%ht1,tor)
        chi(1,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'2MN') then
    if (ua_model.eq.0) then
      if (ok(18)+ok(88)+ok(83)+ok(54).eq.4) then
        call dihed(xyzs%hn1,xyzs%n,xyzs%ct1,xyzs%ht11,tor)
        chi(1,rs) = radian*tor
      end if
      if (ok(18)+ok(89)+ok(84)+ok(54).eq.4) then
        call dihed(xyzs%hn1,xyzs%n,xyzs%ct2,xyzs%ht12,tor)
        chi(2,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'NBU') then
    if (ok(83)+ok(85)+ok(86)+ok(84).eq.4) then
      call dihed(xyzs%ct1,xyzs%cb1,xyzs%cb2,xyzs%ct2,tor)
      chi(1,rs) = radian*tor
    end if
    if (ua_model.eq.0) then
      if (ok(86)+ok(85)+ok(83)+ok(88).eq.4) then
        call dihed(xyzs%cb2,xyzs%cb1,xyzs%ct1,xyzs%ht11,tor)
        chi(2,rs) = radian*tor
      end if
      if (ok(85)+ok(86)+ok(84)+ok(89).eq.4) then
        call dihed(xyzs%cb1,xyzs%cb2,xyzs%ct2,xyzs%ht12,tor)
        chi(3,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'EMT') then
    if (ok(40)+ok(81)+ok(1)+ok(39).eq.4) then
      call dihed(xyzs%cnt,xyzs%s,xyzs%cb,xyzs%cct,tor)
      chi(1,rs) = radian*tor
    end if
    if (ua_model.eq.0) then
      if (ok(1)+ok(81)+ok(40)+ok(37).eq.4) then
        call dihed(xyzs%cb,xyzs%s,xyzs%cnt,xyzs%hnt1,tor)
        chi(2,rs) = radian*tor
      end if
      if (ok(81)+ok(1)+ok(39)+ok(36).eq.4) then
        call dihed(xyzs%s,xyzs%cb,xyzs%cct,xyzs%hct1,tor)
        chi(3,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'IBU') then
    if (ok(84)+ok(1)+ok(83)+ok(88).eq.4) then
      call dihed(xyzs%ct2,xyzs%cb,xyzs%ct1,xyzs%ht11,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(91)+ok(1)+ok(84)+ok(89).eq.4) then
      call dihed(xyzs%ct3,xyzs%cb,xyzs%ct2,xyzs%ht12,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(84)+ok(1)+ok(91)+ok(90).eq.4) then
      call dihed(xyzs%ct2,xyzs%cb,xyzs%ct3,xyzs%ht13,tor)
      chi(3,rs) = radian*tor
    end if
  else if (resname.eq.'AC-') then
    if (ok(52)+ok(30)+ok(20)+ok(56).eq.4) then
      call dihed(xyzs%h1,xyzs%ch3,xyzs%c,xyzs%oxt1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'MOH') then
    if (ok(32)+ok(31)+ok(27)+ok(38).eq.4) then
      call dihed(xyzs%ht1,xyzs%ct,xyzs%o,xyzs%ho,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'MSH') then
    if (ok(32)+ok(31)+ok(81)+ok(87).eq.4) then
      call dihed(xyzs%ht1,xyzs%ct,xyzs%s,xyzs%hs,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'EOH') then
    if (ua_model.eq.0) then
      if (ok(27)+ok(1)+ok(31)+ok(32).eq.4) then
        call dihed(xyzs%o,xyzs%cb,xyzs%ct,xyzs%ht1,tor)
        chi(1,rs) = radian*tor
      end if
      if (ok(31)+ok(1)+ok(27)+ok(38).eq.4) then
        call dihed(xyzs%ct,xyzs%cb,xyzs%o,xyzs%ho,tor)
        chi(2,rs) = radian*tor
      end if
    else
      if (ok(31)+ok(1)+ok(27)+ok(38).eq.4) then
        call dihed(xyzs%ct,xyzs%cb,xyzs%o,xyzs%ho,tor)
        chi(1,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'NMF') then
    if (ok(20)+ok(18)+ok(40)+ok(37).eq.4) then
      call dihed(xyzs%c,xyzs%n,xyzs%cnt,xyzs%hnt1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'NMA') then
    if (ok(20)+ok(18)+ok(40)+ok(37).eq.4) then
      call dihed(xyzs%c,xyzs%n,xyzs%cnt,xyzs%hnt1,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(27)+ok(20)+ok(39)+ok(36).eq.4) then
      call dihed(xyzs%o,xyzs%c,xyzs%cct,xyzs%hct1,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'ACA') then
    if (ok(27)+ok(20)+ok(39)+ok(36).eq.4) then
      call dihed(xyzs%o,xyzs%c,xyzs%cct,xyzs%hct1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'PPA') then
    if (ok(27)+ok(20)+ok(39)+ok(43).eq.4) then
      call dihed(xyzs%o,xyzs%c,xyzs%cct,xyzs%cbt,tor)
      chi(1,rs) = radian*tor
    end if
    if (ua_model.eq.0) then
      if (ok(20)+ok(39)+ok(43)+ok(44).eq.4) then
        call dihed(xyzs%c,xyzs%cct,xyzs%cbt,xyzs%hbt1,tor)
        chi(2,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'DMA') then
    if (ok(20)+ok(18)+ok(34)+ok(46).eq.4) then
      call dihed(xyzs%c,xyzs%n,xyzs%c1,xyzs%h11,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(20)+ok(18)+ok(45)+ok(47).eq.4) then
      call dihed(xyzs%c,xyzs%n,xyzs%c2,xyzs%h21,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(18)+ok(20)+ok(39)+ok(36).eq.4) then
      call dihed(xyzs%n,xyzs%c,xyzs%cct,xyzs%hct1,tor)
      chi(3,rs) = radian*tor
    end if
  else if (resname.eq.'TOL') then
    if (ok(35)+ok(34)+ok(31)+ok(32).eq.4) then
      call dihed(xyzs%c21,xyzs%c1,xyzs%ct,xyzs%ht1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'PCR') then
    if (ua_model.eq.0) then
      if (ok(35)+ok(34)+ok(31)+ok(36).eq.4) then
        call dihed(xyzs%c21,xyzs%c1,xyzs%ct,xyzs%hct1,tor)
        chi(1,rs) = radian*tor
      end if
      if (ok(41)+ok(42)+ok(27)+ok(38).eq.4) then
        call dihed(xyzs%c31,xyzs%c4,xyzs%o,xyzs%ho,tor)
        chi(2,rs) = radian*tor
      end if
    else
      if (ok(41)+ok(42)+ok(27)+ok(38).eq.4) then
        call dihed(xyzs%c31,xyzs%c4,xyzs%o,xyzs%ho,tor)
        chi(1,rs) = radian*tor
      end if
    end if
  else if ((resname.eq.'IME').OR.(resname.eq.'IMD')) then
    if (ok(57)+ok(34)+ok(31)+ok(32).eq.4) then
      call dihed(xyzs%n2,xyzs%c1,xyzs%ct,xyzs%ht1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'MIN') then
    if (ok(45)+ok(82)+ok(31)+ok(32).eq.4) then
      call dihed(xyzs%c2,xyzs%c3,xyzs%ct,xyzs%ht1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'ABA') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'NVA') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(3).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'NLE') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(3).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(3)+ok(4).eq.4) then
      call dihed(xyzs%cb,xyzs%cg,xyzs%cd,xyzs%ce,tor)
      chi(3,rs) = radian*tor
    end if
  else if (resname.eq.'VAL') then
    if (ok(18)+ok(19)+ok(1)+ok(6).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg1,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'LEU') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(7).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd1,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'ILE') then
    if (ok(18)+ok(19)+ok(1)+ok(6).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg1,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(6)+ok(7).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg1,xyzs%cd1,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'MET') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(12).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%sd,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(12)+ok(4).eq.4) then
      call dihed(xyzs%cb,xyzs%cg,xyzs%sd,xyzs%ce,tor)
      chi(3,rs) = radian*tor
    end if
  else if (resname.eq.'SER') then
    if (ok(18)+ok(19)+ok(1)+ok(9).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%og,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(9)+ok(10).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%og,xyzs%hg,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'SEP') then
    if (ok(18)+ok(19)+ok(1)+ok(9).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%og,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(9)+ok(61).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%og,xyzs%np,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(9)+ok(61)+ok(49).eq.4) then
      call dihed(xyzs%cb,xyzs%og,xyzs%np,xyzs%oh,tor)
      chi(3,rs) = radian*tor
    end if
    if (ok(9)+ok(61)+ok(49)+ok(68).eq.4) then
      call dihed(xyzs%og,xyzs%np,xyzs%oh,xyzs%nhop,tor)
      chi(4,rs) = radian*tor
    end if
  else if (resname.eq.'CYS') then
    if (ok(18)+ok(19)+ok(1)+ok(8).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%sg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(8)+ok(10).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%sg,xyzs%hg,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'CYX') then
    if (ok(18)+ok(19)+ok(1)+ok(8).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%sg,tor)
      chi(1,rs) = radian*tor
    end if
  else if (resname.eq.'THR') then
    if (ok(18)+ok(19)+ok(1)+ok(22).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%og1,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(22)+ok(11).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%og1,xyzs%hg1,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'TPO') then
    if (ok(18)+ok(19)+ok(1)+ok(22).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%og1,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(22)+ok(61).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%og1,xyzs%np,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(22)+ok(61)+ok(49).eq.4) then
      call dihed(xyzs%cb,xyzs%og1,xyzs%np,xyzs%oh,tor)
      chi(3,rs) = radian*tor
    end if
    if (ok(22)+ok(61)+ok(49)+ok(68).eq.4) then
      call dihed(xyzs%og1,xyzs%np,xyzs%oh,xyzs%nhop,tor)
      chi(4,rs) = radian*tor
    end if
  else if ((resname.eq.'ASN').OR.(resname.eq.'ASP').OR.(resname.eq.'ASH')) then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(13).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%od1,tor)
      chi(2,rs) = radian*tor
    end if
    if (resname.eq.'ASH') then
      if (ok(1)+ok(2)+ok(98)+ok(100).eq.4) then
        call dihed(xyzs%cb,xyzs%cg,xyzs%od2,xyzs%hd2,tor)
        chi(3,rs) = radian*tor
      end if
    end if
  else if ((resname.eq.'GLN').OR.(resname.eq.'GLU').OR.(resname.eq.'GLH')) then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(3).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(3)+ok(21).eq.4) then
      call dihed(xyzs%cb,xyzs%cg,xyzs%cd,xyzs%oe1,tor)
      chi(3,rs) = radian*tor
    end if
    if (resname.eq.'GLH') then
      if (ok(2)+ok(3)+ok(97)+ok(99).eq.4) then
        call dihed(xyzs%cg,xyzs%cd,xyzs%oe2,xyzs%he2,tor)
        chi(4,rs) = radian*tor
      end if
    end if
  else if ((resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &         (resname.eq.'HIP')) then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(14).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%nd1,tor)
      chi(2,rs) = radian*tor
    end if
  else if ((resname.eq.'PHE').OR.(resname.eq.'TYR').OR.&
 &         (resname.eq.'TRP').OR.(resname.eq.'TYO')) then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(7).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd1,tor)
      chi(2,rs) = radian*tor
    end if
    if (resname.eq.'TYR') then
      if (ok(49)+ok(15)+ok(50)+ok(51).eq.4) then
        call dihed(xyzs%ce2,xyzs%cz,xyzs%oh,xyzs%hh,tor)
        chi(3,rs) = radian*tor
      end if
    end if
  else if (resname.eq.'PTR') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(7).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd1,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(49)+ok(15)+ok(101)+ok(61).eq.4) then
      call dihed(xyzs%ce2,xyzs%cz,xyzs%oz,xyzs%np,tor)
      chi(3,rs) = radian*tor
    end if
    if (ok(15)+ok(101)+ok(61)+ok(49).eq.4) then
      call dihed(xyzs%cz,xyzs%oz,xyzs%np,xyzs%oh,tor)
      chi(4,rs) = radian*tor
    end if
    if (ok(101)+ok(61)+ok(49)+ok(68).eq.4) then
      call dihed(xyzs%oz,xyzs%np,xyzs%oh,xyzs%nhop,tor)
      chi(5,rs) = radian*tor
    end if
  else if ((resname.eq.'LYS').OR.(resname.eq.'LYD').OR.(resname.eq.'KAC').OR.(resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &         (resname.eq.'KM3')) then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(3).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(3)+ok(4).eq.4) then
      call dihed(xyzs%cb,xyzs%cg,xyzs%cd,xyzs%ce,tor)
      chi(3,rs) = radian*tor 
    end if
    if (ok(1)+ok(2)+ok(3)+ok(17).eq.4) then
      call dihed(xyzs%cg,xyzs%cd,xyzs%ce,xyzs%nz,tor)
      chi(4,rs) = radian*tor
    end if
    if (resname.eq.'LYD') then
      if (ok(2)+ok(3)+ok(17)+ok(102).eq.4) then
        call dihed(xyzs%cd,xyzs%ce,xyzs%nz,xyzs%hnz1,tor)
        chi(5,rs) = radian*tor
      end if
    else if ((resname.eq.'KM2').OR.(resname.eq.'KM3')) then
      if (ok(2)+ok(3)+ok(17)+ok(103).eq.4) then
        call dihed(xyzs%cd,xyzs%ce,xyzs%nz,xyzs%ch1,tor)
        chi(5,rs) = radian*tor
      end if
    else if ((resname.eq.'KAC').OR.(resname.eq.'KM1')) then
      if (ok(2)+ok(3)+ok(17)+ok(104).eq.4) then
        call dihed(xyzs%cd,xyzs%ce,xyzs%nz,xyzs%ch,tor)
        chi(5,rs) = radian*tor
      end if
      if (resname.eq.'KAC') then
        if (ok(3)+ok(17)+ok(104)+ok(105).eq.4) then
          call dihed(xyzs%ce,xyzs%nz,xyzs%ch,xyzs%ci,tor)
          chi(6,rs) = radian*tor
        end if
      end if
    end if
  else if (resname.eq.'ORN') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(3).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(3)+ok(16).eq.4) then
      call dihed(xyzs%cb,xyzs%cg,xyzs%cd,xyzs%ne,tor)
      chi(3,rs) = radian*tor 
    end if
  else if (resname.eq.'DAB') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(53).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%nd,tor)
      chi(2,rs) = radian*tor
    end if
  else if (resname.eq.'ARG') then
    if (ok(18)+ok(19)+ok(1)+ok(2).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%cb,xyzs%cg,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(19)+ok(1)+ok(2)+ok(3).eq.4) then
      call dihed(xyzs%ca,xyzs%cb,xyzs%cg,xyzs%cd,tor)
      chi(2,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(3)+ok(16).eq.4) then
      call dihed(xyzs%cb,xyzs%cg,xyzs%cd,xyzs%ne,tor)
      chi(3,rs) = radian*tor
    end if
    if (ok(1)+ok(2)+ok(16)+ok(15).eq.4) then
      call dihed(xyzs%cg,xyzs%cd,xyzs%ne,xyzs%cz,tor)
      chi(4,rs) = radian*tor
    end if
  else if ((resname.eq.'R5P').OR.(resname.eq.'RIB')) then
    if (ok(65)+ok(71)+ok(72)+ok(73).eq.4) then
      call dihed(xyzs%nc3,xyzs%nc2,xyzs%no2,xyzs%nho2,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(74)+ok(75)+ok(70)+ok(71).eq.4) then
      call dihed(xyzs%nc2,xyzs%nc1,xyzs%no1s,xyzs%nho1,tor)
      chi(2,rs) = radian*tor
    end if
  else if ((resname.eq.'D5P').OR.(resname.eq.'DIB')) then
    if (ok(74)+ok(75)+ok(70)+ok(71).eq.4) then
      call dihed(xyzs%nc2,xyzs%nc1,xyzs%no1s,xyzs%nho1,tor)
      chi(1,rs) = radian*tor
    end if
  else if ((resname.eq.'RPG').OR.(resname.eq.'RPA').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'RIA')) then
    if (ok(65)+ok(71)+ok(72)+ok(73).eq.4) then
      call dihed(xyzs%nc3,xyzs%nc2,xyzs%no2,xyzs%nho2,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(42)+ok(59)+ok(70)+ok(71).eq.4) then
      call dihed(xyzs%nc2,xyzs%nc1,xyzs%n9,xyzs%c4,tor)
      chi(2,rs) = radian*tor
    end if
  else if ((resname.eq.'DPG').OR.(resname.eq.'DPA').OR.&
 &         (resname.eq.'DIG').OR.(resname.eq.'DIA')) then
    if (ok(42)+ok(59)+ok(70)+ok(71).eq.4) then
      call dihed(xyzs%nc2,xyzs%nc1,xyzs%n9,xyzs%c4,tor)
      chi(1,rs) = radian*tor
    end if
  else if ((resname.eq.'RPU').OR.(resname.eq.'RPC').OR.&
 &         (resname.eq.'RPT').OR.(resname.eq.'RIU').OR.&
 &         (resname.eq.'RIC').OR.(resname.eq.'RIT')) then
    if (ok(65)+ok(71)+ok(72)+ok(73).eq.4) then
      call dihed(xyzs%nc3,xyzs%nc2,xyzs%no2,xyzs%nho2,tor)
      chi(1,rs) = radian*tor
    end if
    if (ok(45)+ok(28)+ok(70)+ok(71).eq.4) then
      call dihed(xyzs%nc2,xyzs%nc1,xyzs%n1,xyzs%c2,tor)
      chi(2,rs) = radian*tor
    end if
    if ((resname.eq.'RPT').OR.(resname.eq.'RIT')) then
      if (ua_model.eq.0) then
        if (ok(42)+ok(80)+ok(79)+ok(78).eq.4) then
          call dihed(xyzs%c4,xyzs%c5,xyzs%c5m,xyzs%thyh1,tor)
          chi(3,rs) = radian*tor
        end if
      end if
    end if
  else if ((resname.eq.'DPU').OR.(resname.eq.'DPC').OR.&
 &         (resname.eq.'DPT').OR.(resname.eq.'DIU').OR.&
 &         (resname.eq.'DIC').OR.(resname.eq.'DIT')) then
    if (ok(45)+ok(28)+ok(70)+ok(71).eq.4) then
      call dihed(xyzs%nc2,xyzs%nc1,xyzs%n1,xyzs%c2,tor)
      chi(1,rs) = radian*tor
    end if
    if ((resname.eq.'DPT').OR.(resname.eq.'DIT')) then
      if (ua_model.eq.0) then
        if (ok(42)+ok(80)+ok(79)+ok(78).eq.4) then
          call dihed(xyzs%c4,xyzs%c5,xyzs%c5m,xyzs%thyh1,tor)
          chi(2,rs) = radian*tor
        end if
      end if
    end if
  end if
!
end
!
!----------------------------------------------------------------------------------
!
! this one just extracts omega angles for caps/amides
!
subroutine readpdb1_omegas(rs,resname,xyzs,ok,ok2)
!
  use iounit
  use pdb
  use atoms
  use sequen
  use molecule
  use polypep
  use aminos
  use fyoc
  use math
!
  implicit none
!
  character(3) resname
  integer rs
  integer ok(120),ok2(120)
  RTYPE tor
  type(t_pdbxyzs) xyzs

! take care of omega angles for NME / NMA / NMF
  if (resname.eq.'NME') then
    if (ok2(2)+ok2(3)+ok(18)+ok(30).eq.4) then
      call dihed(xyzs%cal,xyzs%cl,xyzs%n,xyzs%ch3,tor)
      omega(rs) = radian*tor
!      write(ilog,*) 'OMEGA ',rs,': ',omega(rs)
    end if
  else if (resname.eq.'NMF') then
    if (ok(27)+ok(20)+ok(18)+ok(40).eq.4) then
      call dihed(xyzs%o,xyzs%c,xyzs%n,xyzs%cnt,tor)
      omega(rs) = radian*tor
    end if
  else if (resname.eq.'NMA') then
    if (ok(27)+ok(20)+ok(18)+ok(40).eq.4) then
      call dihed(xyzs%o,xyzs%c,xyzs%n,xyzs%cnt,tor)
      omega(rs) = radian*tor
    end if
  end if
!
end
!
!----------------------------------------------------------------------------------
!
! this one assigns general backbone angles
!
subroutine readpdb1_bbs(rs,resname,xyzs,ok,ok2)
!
  use iounit
  use pdb
  use atoms
  use sequen
  use molecule
  use polypep
  use aminos
  use fyoc
  use math
  use zmatrix
  use system
!
  implicit none
!
  character(3) resname
  integer rs
  integer ok(120),ok2(120)
  RTYPE tor
  type(t_pdbxyzs) xyzs
!
!
! and now remaining backbone angles
!
  if ((rs.eq.rsmol(molofrs(rs),2)).AND.(yline(rs).gt.0)) then
!   C-terminal psi-angles are different
    if (ok(18)+ok(19)+ok(20)+ok(48).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%c,xyzs%oxt1,tor)
      if (tor.lt.0.0) then
        psi(rs) = radian*tor + 180.0
      else
        psi(rs) = radian*tor - 180.0
      end if
    end if
!   standard psi
  else if ((yline(rs).gt.0).AND.(resname.ne.'ACE')) then
    if (ok(18)+ok(19)+ok(20)+ok(27).eq.4) then
      call dihed(xyzs%n,xyzs%ca,xyzs%c,xyzs%o,tor)
      if (tor.lt.0.0) then
        psi(rs) = radian*tor + 180.0
      else
        psi(rs) = radian*tor - 180.0
      end if
!    write(ilog,*) 'PSI ',rs-1,': ',psi(rs-1)
    end if
  end if
  if ((fline(rs).gt.0).AND.(rs.eq.rsmol(molofrs(rs),1))) then
!   N-terminal phi is different
    if (ok(20)+ok(19)+ok(18)+ok(52).eq.4) then
      call dihed(xyzs%c,xyzs%ca,xyzs%n,xyzs%h1,tor)
      phi(rs) = radian*tor
!     for proline, we have to set Z-matrix directly
      if (resname.eq.'PRO') then 
        ztor(fline(rs)) = radian*tor
      end if
    end if
  else if (fline(rs).gt.0) then
!   standard phi
    if (ok2(3)+ok(18)+ok(19)+ok(20).eq.4) then
      call dihed(xyzs%cl,xyzs%n,xyzs%ca,xyzs%c,tor)
      phi(rs) = radian*tor
!      write(ilog,*) 'PHI ',rs,': ',phi(rs)
    end if
  end if
  if (wline(rs).gt.0) then
    if (seqtyp(rs-1).eq.28) then
      if (ua_model.gt.0) then
        if (ok(18).ne.1) then
           write(ilog,*) 'Fatal. Cannot read-in coordinates for nitrogen in resid&
 &ue ',rs,' (',resname,') which is required for molecule ',molofrs(rs),'.'
          call fexit()
        end if
        if (ok2(9)+ok2(3)+ok(18)+ok(19).eq.4) then
          call dihed(xyzs%ol,xyzs%cl,xyzs%n,xyzs%ca,tor)
          omega(rs) = radian*tor
        end if
!       in addition, we need to catch up with the fact that FOR only read in
!       2 RB reference atoms
        x(at(rs)%bb(1)) = xyzs%n(1)
        y(at(rs)%bb(1)) = xyzs%n(2)
        z(at(rs)%bb(1)) = xyzs%n(3)
      else
        if (ok2(2)+ok2(3)+ok(18)+ok(19).eq.4) then
          call dihed(xyzs%cal,xyzs%cl,xyzs%n,xyzs%ca,tor)
          omega(rs) = radian*tor
        end if
      end if
    else
      if (ok2(2)+ok2(3)+ok(18)+ok(19).eq.4) then
        call dihed(xyzs%cal,xyzs%cl,xyzs%n,xyzs%ca,tor)
        omega(rs) = radian*tor
!      write(ilog,*) 'OMEGA ',rs,': ',omega(rs)
      end if
    end if
  end if
!   5'-terminal nucleosides
  if ((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)) then ! DIB, RIB, RIX, or DIX
    if (nucsline(1,rs).gt.0) then
      if (ok(63)+ok(62)+ok(65)+ok(77).eq.4) then
        call dihed(xyzs%nc4,xyzs%nc5,xyzs%no5,xyzs%nho5,tor)
        nucs(1,rs) = radian*tor
      end if
    end if
    if (nucsline(2,rs).gt.0) then
      if (ok(65)+ok(62)+ok(63)+ok(64).eq.4) then
        call dihed(xyzs%no5,xyzs%nc5,xyzs%nc4,xyzs%nc3,tor)
        nucs(2,rs) = radian*tor
      end if
    end if
!   now the sugar, note that we edit the Z-matrix directly
!   and that C5-C4-C3-O3(+1) is dealt with below
    if (ok(63)+ok(64)+ok(65)+ok(70)+ok(71)+ok(76).eq.6) then
      call dihed(xyzs%nc5,xyzs%nc4,xyzs%nc3,xyzs%nc2,tor)
      ztor(at(rs)%sc(1)) = radian*tor
      call dihed(xyzs%nc4,xyzs%nc3,xyzs%nc2,xyzs%nc1,tor)
      ztor(at(rs)%sc(2)) = radian*tor
      call dihed(xyzs%nc2,xyzs%nc3,xyzs%nc4,xyzs%no4,tor)
      ztor(at(rs)%sc(3)) = radian*tor
      call bondang2(xyzs%nc4,xyzs%nc3,xyzs%nc2,tor)
      bang(at(rs)%sc(1)) = radian*tor
      call bondang2(xyzs%nc3,xyzs%nc2,xyzs%nc1,tor)
      bang(at(rs)%sc(2)) = radian*tor
      call bondang2(xyzs%nc3,xyzs%nc4,xyzs%no4,tor)
      bang(at(rs)%sc(3)) = radian*tor
    end if
!   note we can't assign further nucs yet
  else if (seqpolty(rs).eq.'N') then
    if (nucsline(1,rs).gt.0) then
!     regular 5'-P-nucleotides
      if (rs.eq.rsmol(molofrs(rs),1)) then
        if (ok(68)+ok(66)+ok(61)+ok(62).eq.4) then
          call dihed(xyzs%no5,xyzs%np,xyzs%no3,xyzs%nhop,tor)
          nucs(1,rs) = radian*tor
        end if
      else
        if (ok2(4)+ok(66)+ok(61)+ok(62).eq.4) then
          call dihed(xyzs%nc3l,xyzs%no3,xyzs%np,xyzs%no5,tor)
          nucs(1,rs) = radian*tor
        end if
      end if
    end if
    if (nucsline(2,rs).gt.0) then
      if (ok(66)+ok(61)+ok(62)+ok(63).eq.4) then
        call dihed(xyzs%no3,xyzs%np,xyzs%no5,xyzs%nc5,tor)
        nucs(2,rs) = radian*tor
      end if
    end if
    if (nucsline(3,rs).gt.0) then
      if (ok(61)+ok(62)+ok(63)+ok(64).eq.4) then
        call dihed(xyzs%np,xyzs%no5,xyzs%nc5,xyzs%nc4,tor)
        nucs(3,rs) = radian*tor
      end if
    end if
    if (nucsline(4,rs).gt.0) then
      if (ok(65)+ok(62)+ok(63)+ok(64).eq.4) then
        call dihed(xyzs%no5,xyzs%nc5,xyzs%nc4,xyzs%nc3,tor)
        nucs(4,rs) = radian*tor
      end if
    end if
!   now the sugar, note that we edit the Z-matrix directly
!   and that C5-C4-C3-O3(+1) is dealt with below
    if (ok(63)+ok(64)+ok(65)+ok(70)+ok(71)+ok(76).eq.6) then
      call dihed(xyzs%nc5,xyzs%nc4,xyzs%nc3,xyzs%nc2,tor)
      ztor(at(rs)%sc(1)) = radian*tor
      call dihed(xyzs%nc4,xyzs%nc3,xyzs%nc2,xyzs%nc1,tor)
      ztor(at(rs)%sc(2)) = radian*tor
      call dihed(xyzs%nc2,xyzs%nc3,xyzs%nc4,xyzs%no4,tor)
      ztor(at(rs)%sc(3)) = radian*tor
      call bondang2(xyzs%nc4,xyzs%nc3,xyzs%nc2,tor)
      bang(at(rs)%sc(1)) = radian*tor
      call bondang2(xyzs%nc3,xyzs%nc2,xyzs%nc1,tor)
      bang(at(rs)%sc(2)) = radian*tor
      call bondang2(xyzs%nc3,xyzs%nc4,xyzs%no4,tor)
      bang(at(rs)%sc(3)) = radian*tor
    end if
!   note we can't assign further nucs yet
!
    if (rs.gt.rsmol(molofrs(rs),1)) then
!    the following is the bond around the sugar ring
!    we edit the Z-matrix directly
      if (seqpolty(rs-1).eq.'N') then
        if (ok2(5)+ok2(4)+ok2(8)+ok(66).eq.4) then
          call dihed(xyzs%nc5l,xyzs%nc4l,xyzs%nc3l,xyzs%no3,tor)
          ztor(nuci(rs,1)) = radian*tor
        end if
      end if
      if ((seqtyp(rs-1).ge.76).AND.(seqtyp(rs-1).le.87)) then
        if (nucsline(3,rs-1).gt.0) then
          if (ok2(5)+ok2(4)+ok(66)+ok(61).eq.4) then
            call dihed(xyzs%nc4l,xyzs%nc3l,xyzs%no3,xyzs%np,tor)
            nucs(3,rs-1) = radian*tor
          end if
        end if
      else if (seqpolty(rs-1).eq.'N') then
        if (nucsline(5,rs-1).gt.0) then
          if (ok2(5)+ok2(4)+ok(66)+ok(61).eq.4) then
            call dihed(xyzs%nc4l,xyzs%nc3l,xyzs%no3,xyzs%np,tor)
            nucs(5,rs-1) = radian*tor
          end if
        end if
      end if
    end if
!   the following is the bond around the sugar ring
!   we edit the Z-matrix directly
    if (rs.eq.rsmol(molofrs(rs),2)) then
      if (ok(64)+ok(65)+ok(67)+ok(69).eq.4) then
        call dihed(xyzs%nc5,xyzs%nc4,xyzs%nc3,xyzs%no32,tor)
        ztor(at(rs)%bb(9)) = radian*tor
      end if
      if (nucsline(5,rs).gt.0) then
        if (ok(64)+ok(65)+ok(67)+ok(69).eq.4) then
          call dihed(xyzs%nc4,xyzs%nc3,xyzs%no32,xyzs%nho3,tor)
          nucs(5,rs) = radian*tor
        end if
      end if
    end if
!
  end if
!
end
!  
!----------------------------------------------------------------------------------
!
! if for analysis purposes we "simulate" through a pdb-trajectory, the data
! can either be available in single-file format or in separate snapshot files
! the latter is set up here (mostly extraction of format and sanity checks)
!
subroutine setup_pdbsnaps()
!
  use iounit
  use pdb
  use atoms
  use sequen
  use system
  use params
  use clusters
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer i,n1,dummy,j,k,dummy2,ipdb,t1,t2,t3,t4,t5,t6,freeunit,dummy3,azero
  logical exists
  character(80) string
  character(MAXSTRLEN) fn,pdbsnap
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
  integer tslash
#endif
!
  if ((use_pdb_template.EQV..true.).AND.(align%yes.EQV..false.).AND.(n_pdbunk.le.0)) then
    write(ilog,*) 'Warning. Disabling use of pdb-template due to analyzing data that are already in &
 &annotated (pdb) format.'
    use_pdb_template = .false.
  else if (use_pdb_template.EQV..true.) then
    use_pdb_template = .false.
  end if
!
  azero = 0
  ipdb = freeunit()
  call strlims(pdbinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (pdbinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  pdbinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//pdbinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//pdbinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = pdbinfile(t1:t2)
#endif
  inquire (file=fn(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &first pdb-snapshot (',fn(t1:t2),') in setup_pdbsnaps().'
    call fexit()
  end if
!
! search for numbers in string (after last forward slash)
  do i=t2,t1,-1
    if (pdbinfile(i:i).eq.'/') exit
  end do
  if (i.lt.t1) i = t1
  call extract_int2(pdbinfile,n1,i,k,pdb_format(3))
  if (i.eq.t1) then
    write(ilog,*) 'Fatal. Base filename for reading pdb-snapshots ('&
 &,pdbinfile(t1:t2),') contains no number. Uninterpretable!'
    call fexit()
  end if
  j = i
  call extract_int2(pdbinfile,dummy,j,dummy2,dummy3)
  if (j.ne.i) then
    write(ilog,*) 'Fatal. Base filename for reading pdb-snapshots ('&
 &,pdbinfile(t1:t2),') contains 2 or more numbers. Uninterpretable!'
    call fexit()
  end if
!
! store the ending
  if (i.le.t2) then
    pdb_suff(1:t2-i+1) = pdbinfile(i:t2)
    pdb_format(2) = t2-i+1
  end if
!
! store the beginning
  if (k.gt.t1) then
#ifdef ENABLE_MPI
    do i=k-1,t1,-1
      if (pdbinfile(i:i).eq.SLASHCHAR) exit
    end do
    tslash = i
    call int2str(myrank,xpont,re_aux(10))
    if ((tslash.ge.t1).AND.(tslash.lt.(k-1))) then
      pdb_pref = pdbinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//pdbinfile((tslash+1):(k-1))
    else
      pdb_pref =  'N_'//xpont(1:re_aux(10))//'_'//pdbinfile(t1:t2)
    end if
    call strlims(pdb_pref,t1,t2)
    if (t1.ne.1) then
      write(ilog,*) 'Fatal. Base filename for reading pdb-snapshots (',pdbinfile(t1:t2),') is uninterpretable for MPI processing!'
      call fexit()
    end if
    pdb_format(1) = t2-t1+1
#else
    pdb_pref(1:(k-t1)) = pdbinfile(t1:(k-1))
    pdb_format(1) = k-t1
#endif
  end if
!
! store the offset
  pdb_format(4) = n1 - 1
!
  if (be_unsafe.EQV..false.) then
!   probe existence of all the required files
    do i=n1,n1+nsim-1
      if (pdb_format(3).gt.0) then
        call int2str(i,string,pdb_format(3))
        call strlims(string,t3,t4)
        if ((pdb_format(1).gt.0).AND.(pdb_format(2).gt.0)) then
          pdbsnap = pdb_pref(1:pdb_format(1))//string(1:pdb_format(3))&
   &//pdb_suff(1:pdb_format(2))
        else if (pdb_format(1).gt.0) then
          pdbsnap = pdb_pref(1:pdb_format(1))//string(1:pdb_format(3))
        else if (pdb_format(2).gt.0) then
          pdbsnap = string(1:pdb_format(3))//pdb_suff(1:pdb_format(2))
        else
          pdbsnap = string(1:pdb_format(3))
        end if
      else
        call int2str(i,string,azero)
        call strlims(string,t3,t4)
        if ((pdb_format(1).gt.0).AND.(pdb_format(2).gt.0)) then
          pdbsnap = pdb_pref(1:pdb_format(1))//string(t3:t4)&
   &//pdb_suff(1:pdb_format(2))
        else if (pdb_format(1).gt.0) then
          pdbsnap = pdb_pref(1:pdb_format(1))//string(t3:t4)
        else if (pdb_format(2).gt.0) then
          pdbsnap = string(t3:t4)//pdb_suff(1:pdb_format(2))
        else
          pdbsnap = string(t3:t4)
        end if
      end if
      ipdb = freeunit()
      call strlims(pdbsnap,t5,t6)
      inquire (file=pdbsnap(t5:t6),exist=exists)
      if (exists.EQV..false.) then
        write(ilog,*) 'Fatal. Cannot open snapshot with number ',i,' (&
   &',pdbsnap(t5:t6),') requested through number of steps.'
        call fexit()
      end if
    end do
  end if ! be_unsafe is false
!
! and allocate memory for global arrays
  if (allocated(xpdb).EQV..true.) deallocate(xpdb)
  if (allocated(ypdb).EQV..true.) deallocate(ypdb)
  if (allocated(zpdb).EQV..true.) deallocate(zpdb)
!
end
!  
!----------------------------------------------------------------------------------
!
#ifdef LINK_XDR
!
! if for analysis purposes we "simulate" through an xtc-trajectory, we'll
! make sure that the file exists
! we will, however, not ensure that the file is sane (with respect to nsim)
! we do check whether the number of atoms is correct.
!
subroutine setup_xtctraj()
!
  use iounit
  use pdb
  use mcsums
  use atoms
  use sequen
  use clusters
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer t2,t1,i,azero,freeunit,magic,natoms,kk,kkk,ret
  logical exists,dofatal
  character rmode
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
  integer tslash
#endif
!
  dofatal = .true.
  azero = 0
!
! note this will temporarily abuse ipdbtraj which is why we do it first
  if (use_pdb_template.EQV..true.) then
!   
    call setup_pdbtemplate()
!   init map
    do i=1,n
      pdbmap(i) = -1
    end do
!   fake-read template file to get map
    call FMSMC_readpdb2(azero)
    close(unit=ipdbtraj) ! opened in setup_pdbtemplate
    call strlims(pdbtmplfile,t1,t2)
    do i=1,n
      if ((pdbmap(i).le.0).OR.(pdbmap(i).gt.n)) then
        write(ilog,*) 'Fatal. Incomplete map extracted from pdb-template (',pdbtmplfile(t1:t2),').&
  & Check for unusual atom naming in pdb-file and/or setting for FMCSC_PDB_R_CONV.'
        call fexit()
      end if
    end do
!   the template may serve a double fxn for extracting reference coordinates
    if (align%yes.EQV..true.) then
      call init_structalign()
    end if
!
  end if
!
  ipdbtraj = freeunit()
  call strlims(xtcinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (xtcinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  xtcinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//xtcinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//xtcinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = xtcinfile(t1:t2)
#endif
  inquire (file=fn(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from xtc (',fn(t1:t2),') in setup_xtctraj().'
    call fexit()
  end if
!
! open
  rmode = "r"
  call xdrfopen(ipdbtraj,fn(t1:t2),rmode,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
!
! get the header-line
  call xdrfint(ipdbtraj,magic,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  call xdrfint(ipdbtraj,natoms,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  if (natoms.ne.n) then
    write(ilog,*) 'Fatal. XTC-file (',fn(t1:t2),') does not m&
 &atch the number of atoms defined by the sequence-file.'
    call fexit()
  end if
  call xdrfint(ipdbtraj,kk,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  if (kk.ne.1) then
    write(ilog,*) 'Warning. XTC-file does not seem to begin with a s&
 &napshot labelled as the first one (#1). The file might be corrupte&
 &d or might have been altered.'
  end if
  call xdrfint(ipdbtraj,kkk,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
!
! close and re-open
  call xdrfclose(ipdbtraj,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  call xdrfopen(ipdbtraj,fn(t1:t2),rmode,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
!
end
!
!-----------------------------------------------------------------------
!
subroutine close_xtctraj()
!
  use mcsums
  use pdb
!
  implicit none
!
  integer ret
  logical dofatal
!
  dofatal = .true.
!
  if (allocated(xpdb).EQV..true.) deallocate(xpdb)
  if (allocated(ypdb).EQV..true.) deallocate(ypdb)
  if (allocated(zpdb).EQV..true.) deallocate(zpdb)
  if (allocated(pdbmap).EQV..true.) deallocate(pdbmap)
!
  call xdrfclose(ipdbtraj,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
!
end
!  
!-----------------------------------------------------------------------
!
subroutine fail_readxtc(ret,dofatal)
!
  use iounit
  use sequen
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  logical dofatal
  integer ret,t1,t2
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  integer i,tslash
  character(re_aux(10)) xpont
#endif
!
  call strlims(xtcinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (xtcinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  xtcinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//xtcinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//xtcinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = xtcinfile(t1:t2)
#endif
  if (dofatal.EQV..true.) then
    write(ilog,*) 'Fatal. Reading from xtc-file (',fn(t1:t2),')&
 &exited with an error (got ',ret,'). Sequence mismatch? Number of s&
 &napshots incorrect?'
    call fexit()
  else
    write(ilog,*) 'Fatal. Reading from xtc-file (',fn(t1:t2),')&
 &exited with an error (got ',ret,'). Number of snapshots incorrect?'
  end if
!
end
!
#endif
!  
!----------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
! if for analysis purposes we "simulate" through a netcdf-trajectory, we'll
! make sure that the file exists
! we will, however, not ensure that the file is sane (with respect to nsim)
! we do check whether the number of atoms is correct.
!
subroutine setup_netcdftraj()
!
  use iounit
  use pdb
  use mcsums
  use atoms
  use sequen
  use netcdf
  use system
  use clusters
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer t2,t1,i,athree,azero,nframes,natoms,nspats(3),ret,fndds(5),dimlen
  logical exists
  character(MAXSTRLEN) trystr,ucstr,fn
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
  integer tslash
#endif 
!
  azero = 0
  athree = 3
!
! note this will temporarily abuse ipdbtraj which is why we do it first
  if (use_pdb_template.EQV..true.) then
!   
    call setup_pdbtemplate()
!   init map
    do i=1,n
      pdbmap(i) = -1
    end do
!   fake-read template file to get map
    call FMSMC_readpdb2(azero)
    close(unit=ipdbtraj) ! opened in setup_pdbtemplate
    call strlims(pdbtmplfile,t1,t2)
    do i=1,n
      if ((pdbmap(i).le.0).OR.(pdbmap(i).gt.n)) then
        write(ilog,*) 'Fatal. Incomplete map extracted from pdb-template (',pdbtmplfile(t1:t2),').&
  & Check for unusual atom naming in pdb-file and/or setting for FMCSC_PDB_R_CONV.'
        call fexit()
      end if
    end do
!   the template may serve a double fxn for extracting reference coordinates
    if (align%yes.EQV..true.) then
      call init_structalign()
    end if
!
  end if
!
  call strlims(netcdfinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (netcdfinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  netcdfinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//netcdfinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//netcdfinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = netcdfinfile(t1:t2)
#endif
  inquire (file=fn(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from NetCDF (',fn(t1:t2),') in setup_netcdftraj().'
    call fexit()
  end if
!
! open
 44 format('Warning. Ambigous dimensions in NetCDF file. Encountered ',a,' twice but keeping&
 &only first.')
  call check_fileionetcdf( nf90_open(fn(t1:t2), NF90_NOWRITE, ipdbtraj) )
!
! find the necessary dimensions: three are required
  fndds(:) = 0
  nspats(3) = 0
  nframes = 0
  do i=1,NF90_MAX_DIMS
    ret = nf90_inquire_dimension(ipdbtraj,i,trystr,dimlen)
    if (ret.eq.NF90_NOERR) then
      ucstr(1:15) = trystr(1:15)
      call toupper(ucstr(1:15))
      if (ucstr(1:5).eq.'FRAME') then
        if (fndds(1).eq.1) then 
          write(ilog,44) ucstr(1:5)
        else
          nframes = dimlen
          fndds(1) = 1
          netcdf_ids(21) = i
        end if
      else if (ucstr(1:4).eq.'ATOM') then
        if (fndds(2).eq.1) then 
          write(ilog,44) ucstr(1:4)
        else
          natoms = dimlen
          fndds(2) = 1
          netcdf_ids(22) = i
        end if
      else if (ucstr(1:7).eq.'SPATIAL') then
        if (fndds(3).eq.1) then 
          write(ilog,44) ucstr(1:7)
        else
          nspats(1) = dimlen
          fndds(3) = 1
          netcdf_ids(23) = i
        end if
      else if (ucstr(1:12).eq.'CELL_SPATIAL') then
        if (fndds(4).eq.1) then 
          write(ilog,44) ucstr(1:12)
        else
          nspats(2) = dimlen
          fndds(4) = 1
          netcdf_ids(24) = i
        end if
      else if (ucstr(1:12).eq.'CELL_ANGULAR') then
        if (fndds(5).eq.1) then 
          write(ilog,44) ucstr(1:12)
        else
          nspats(3) = dimlen
          fndds(5) = 1
          netcdf_ids(25) = i
        end if
      end if
    else if (ret.eq.NF90_EBADDIM) then
!     do nothing
    else ! get us out of here
     call check_fileionetcdf( nf90_inquire_dimension(ipdbtraj,i,trystr,dimlen) )
    end if
  end do
!
  if (natoms.ne.n) then
    write(ilog,*) 'Fatal. NetCDF-file (',fn(t1:t2),') does not m&
 &atch the number of atoms defined by the sequence-file.'
    call fexit()
  end if
  if (nspats(1).ne.3) then
    write(ilog,*) 'Fatal. NetCDF-file (',fn(t1:t2),') seems to encode &
 &non-3-dimensional trajectory data which is currently not supported by CAMPARI.'
    call fexit()
  end if
  if ((nspats(2).ne.3).OR.(nspats(3).ne.3)) then
    write(ilog,*) 'Fatal. NetCDF-file (',fn(t1:t2),') does not comply with &
 &AMBER format for unit cell declarations (dimensioanlity).'
    call fexit()
  end if
  if (nframes.lt.nsim) then
    if (nframes.lt.1) then
      write(ilog,*) 'Fatal. NetCDF-file (',fn(t1:t2),') has no frame &
 &data (empty containers).'
      call fexit()
    else
      write(ilog,*) 'Warning. NetCDF-file (',fn(t1:t2),') does not m&
 &atch the number of snapshots requested for analysis (',nframes,' vs. ',nsim,').'
#ifdef ENABLE_MPI
      i = nframes
      call MPI_SyncInt(i,athree)
      nsim = i
      write(ilog,*) 'Across all MPI processes, the maximum possible number of frames was found to be ',nsim,'.'
      if (i.lt.nframes) write(ilog,*) 'This is unlikely to be a meaningful calculation.'
#else
      write(ilog,*) 'Reducing analysis frames accordingly.'
      nsim = nframes
#endif
    end if
  end if
!
! now find the necessary variables, only coordinates are required
  fndds(:) = 0
  ret = nf90_inq_varid(ipdbtraj,"coordinates",netcdf_ids(32))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "coordinates" not found in NetCDF-file (',&
 &fn(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ipdbtraj,"coordinates",netcdf_ids(32)) )
  end if
!
  ret = nf90_inq_varid(ipdbtraj,"cell_lengths",fndds(1))
  if (ret.eq.NF90_NOERR) then
     netcdf_ids(33) = fndds(1)
  else if (ret.eq.NF90_ENOTVAR) then
!    write(ilog,*) 'Warning. Variable "cell_lengths" not found in NetCDF-file (',&
! &fn(t1:t2),'). Assuming fixed volume conditions.'
!    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ipdbtraj,"cell_lengths",fndds(1)) )
  end if
!
  ret = nf90_inq_varid(ipdbtraj,"cell_angles",fndds(1))
  if (ret.eq.NF90_NOERR) then
     netcdf_ids(34) = fndds(1)
  else if (ret.eq.NF90_ENOTVAR) then
!    write(ilog,*) 'Warning. Variable "cell_angles" not found in NetCDF-file (',&
! &fn(t1:t2),'). Assuming fixed volume conditions.'
!    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ipdbtraj,"cell_lengths",fndds(1)) )
  end if
!
! if we're still here, we should have a sane file
!
end
!
!-----------------------------------------------------------------------
!
subroutine close_netcdftraj()
!
  use mcsums
  use pdb
  use netcdf
!
  implicit none
!
  if (allocated(xpdb).EQV..true.) deallocate(xpdb)
  if (allocated(ypdb).EQV..true.) deallocate(ypdb)
  if (allocated(zpdb).EQV..true.) deallocate(zpdb)
  if (allocated(pdbmap).EQV..true.) deallocate(pdbmap)
!
  call check_fileionetcdf( nf90_close(ipdbtraj) )
!
end
!
#endif
!  
!-----------------------------------------------------------------------
!
! if for analysis purposes we "simulate" through a dcd-trajectory, we'll
! make sure that the file exists
! we will, however, not ensure that the file is sane 
! we do check whether the number of atoms is correct.
!
subroutine setup_dcdtraj()
!
  use iounit
  use pdb
  use mcsums
  use atoms
  use sequen
  use system
  use clusters
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer i,t2,t1,freeunit,natoms,iis(100),iomess,azero
  logical exists,dofatal
  real(KIND=4) fakedt
  character(4) shortti
  character(80) longti,timest
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  integer tslash
  character(re_aux(10)) xpont
#endif
!
  dofatal = .true.
  azero = 0
  iis(:) = -1
!
! note this will temporarily abuse ipdbtraj which is why we do it first
  if (use_pdb_template.EQV..true.) then
!   
    call setup_pdbtemplate()
!   init map
    do i=1,n
      pdbmap(i) = -1
    end do
!   fake-read template file to get map
    call FMSMC_readpdb2(azero)
    close(unit=ipdbtraj) ! opened in setup_pdbtemplate
    call strlims(pdbtmplfile,t1,t2)
    do i=1,n
      if ((pdbmap(i).le.0).OR.(pdbmap(i).gt.n)) then
        write(ilog,*) 'Fatal. Incomplete map extracted from pdb-template (',pdbtmplfile(t1:t2),').&
  & Check for unusual atom naming in pdb-file and/or setting for FMCSC_PDB_R_CONV.'
        call fexit()
      end if
    end do
!   the template may serve a double fxn for extracting reference coordinates
    if (align%yes.EQV..true.) then
      call init_structalign()
    end if
!
  end if
!
  ipdbtraj = freeunit()
  call strlims(dcdinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (dcdinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  dcdinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//dcdinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//dcdinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = dcdinfile(t1:t2)
#endif
  inquire (file=fn(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from dcd (',fn(t1:t2),') in setup_dcdtraj().'
    call fexit()
  end if
!
! open
  open(unit=ipdbtraj,status='old',form='unformatted',file=fn(t1:t2),iostat=iomess)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
!
! get the header
  read(ipdbtraj,iostat=iomess) shortti,iis(1),iis(2),iis(3),iis(4),iis(5),iis(6),iis(7),iis(8),iis(9),&
 &        fakedt,iis(11),iis(12),iis(13),iis(14),iis(15),iis(16),iis(17),iis(18),iis(19),iis(20)
  if (iis(9).ne.0) then
    write(ilog,*) 'Fatal. Header of DCD-trajectory file (',fn(t1:t2),') indicates usage&
 & of the "fixed atoms" feature. This is currently not supported.'
    call fexit()
  end if
  if (iis(11).eq.1) then
!   do nothing
  else if (iis(11).eq.0) then
    dcd_withbox = .false.
  else
    write(ilog,*) 'Warning. Header of DCD-trajectory file (',fn(t1:t2),') is unclear on&
 & the inclusion of box coordinates with each frame. This run may cause unexpected crashes since&
 & it assumes box coordinates to be contained.'
  end if
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  read(ipdbtraj,iostat=iomess) iis(19),longti,timest
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  read(ipdbtraj,iostat=iomess) natoms
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  if (natoms.ne.n) then
    write(ilog,*) 'Fatal. DCD-file (',fn(t1:t2),') does not m&
 &atch the number of atoms defined by the sequence-file.'
    call fexit()
  end if
  if (pdb_analyze.EQV..true.) then
    if (iis(1).lt.nsim) then
      write(ilog,*) 'Warning. DCD-file (',fn(t1:t2),') reports having fewer &
 &snapshots (',iis(1),') then requested for analysis. This run may cause unexpected crashes.'
!      nsim = iis(1)
    end if
  end if
!
! close and re-open
  close(unit=ipdbtraj,iostat=iomess)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  open(unit=ipdbtraj,status='old',form='unformatted',file=fn(t1:t2),iostat=iomess)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
!
!
end
!
!-----------------------------------------------------------------------
!
subroutine close_dcdtraj()
!
  use mcsums
  use pdb
!
  implicit none
!
  integer iomess
  logical dofatal
!
  dofatal = .true.
!
  if (allocated(xpdb).EQV..true.) deallocate(xpdb)
  if (allocated(ypdb).EQV..true.) deallocate(ypdb)
  if (allocated(zpdb).EQV..true.) deallocate(zpdb)
  if (allocated(pdbmap).EQV..true.) deallocate(pdbmap)
!
  close(unit=ipdbtraj,iostat=iomess)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
!
end
!  
!-----------------------------------------------------------------------
!
subroutine fail_readdcd(iomess,dofatal)
!
  use iounit
  use sequen
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  logical dofatal
  integer t1,t2,iomess
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  integer i,tslash
  character(re_aux(10)) xpont
#endif
!
  call strlims(dcdinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (dcdinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  dcdinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//dcdinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//dcdinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = dcdinfile(t1:t2)
#endif
!
  if (dofatal.EQV..true.) then
    write(ilog,*) 'Fatal. Reading from dcd-file (',fn(t1:t2),') &
 &exited with an error (got ',iomess,'). Sequence mismatch? Number of s&
 &napshots incorrect?'
    call fexit()
  else
    write(ilog,*) 'Warning. Reading from dcd-file (',fn(t1:t2),') &
 &exited with an error (got ',iomess,'). Number of snapshots incorrect?'
  end if
!
end
!
!----------------------------------------------------------------------------------
!
! if for analysis purposes we "simulate" through a pdb-trajectory, the data
! can either be available in single-file format or in separate snapshot files
! the former is set up here (mostly sanity checks to MODEL/ENDMDL structure)
!
subroutine setup_pdbtraj()
!
  use iounit
  use atoms
  use pdb
  use mcsums
  use sequen
  use system
  use params
  use clusters
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer nl,i,k,t2,t1,c1,c2,freeunit
  integer, ALLOCATABLE:: sta(:),sto(:) 
  character(MAXSTRLEN) pdbline
  character(6) keyw
  logical exists
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
  integer tslash
#endif
!
  if ((use_pdb_template.EQV..true.).AND.(align%yes.EQV..false.).AND.(n_pdbunk.le.0)) then
    write(ilog,*) 'Warning. Disabling use of pdb-template due to analyzing data that are already in &
 &annotated (pdb) format.'
    use_pdb_template = .false.
  else if (use_pdb_template.EQV..true.) then
    use_pdb_template = .false.
  end if
!
  ipdbtraj = freeunit()
  call strlims(pdbinfile,t1,t2)
#ifdef ENABLE_MPI
  if ((pdb_mpimany.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
    do i=t2,t1,-1
      if (pdbinfile(i:i).eq.SLASHCHAR) exit
    end do
    tslash = i
    call int2str(myrank,xpont,re_aux(10))
    if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
      fn =  pdbinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//pdbinfile((tslash+1):t2)
    else
      fn =  'N_'//xpont(1:re_aux(10))//'_'//pdbinfile(t1:t2)
    end if
    call strlims(fn,t1,t2)
    if ((pdb_analyze.EQV..false.).AND.(pdb_readmode.ne.1).AND.(fycxyz.eq.1)) then
      write(ilog,*) 'Warning. Reading multiple pdb-files for parallel simulations can be problematic if the coordinates are &
 &of low precision. This is because covalent geometries will differ slightly between replicas. It is recommended to &
 &use FMCSC_PDB_INPUTSTRING (and FMCSC_PDB_OUTPUTSTRING) to circumvent this issue whenever possible.'
    end if
  else
    fn(t1:t2) = pdbinfile(t1:t2)
  end if
#else
  fn(t1:t2) = pdbinfile(t1:t2)
#endif
  inquire (file=fn(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from pdb (',fn(t1:t2),') in setup_pdbtraj().'
    call fexit()
  end if
!
  if (be_unsafe.EQV..false.) then
    open (unit=ipdbtraj,file=fn(t1:t2),status='old')
    allocate(sta(nsim+10))
    allocate(sto(nsim+10))
    k = 0
    sta(:) = -1
    sto(:) = -1
    c1 = 0
    c2 = 0
    do while (1.eq.1)
      k = k + 1
 1356 format(FORM_MAXSTRLEN)
      read(ipdbtraj,1356,err=30,end=30) pdbline
      read(pdbline,'(a6)') keyw
      if (keyw.eq.'MODEL ') then
        c1 = c1 + 1
        sta(c1) = k
      else if (keyw.eq.'ENDMDL') then
        c2 = c2  + 1
        sto(c2) = k
      end if
!     if necessary snapshots found exit early
      if ((c2.ge.nsim).AND.(c1.ge.nsim)) exit
!     if inconsistent file avoid seg-fault
      if ((c2.ge.(nsim+10)).OR.(c1.ge.(nsim+10))) exit
    end do
   30   k = k - 1
    nl = k
    close(unit=ipdbtraj)
!
!   a few sanity checks
    if ((c1.gt.0).OR.(c2.gt.0)) then
      if ((c1.ne.c2).OR.(sto(1).le.sta(1))) then
        write(ilog,*) 'Fatal. Bad format (check MODEL, ENDMDL) in pdb-&
   &file (',fn(t1:t2),') while executing setup_pdbtraj().'
        deallocate(sta)
        deallocate(sto)
        call fexit()
      end if
      do i=2,c1
        if ((sto(i).le.sta(i)).OR.(sta(i).le.sto(i-1))) then
          write(ilog,*) 'Fatal. Bad format (check MODEL, ENDMDL) in pd&
   &b-file (',fn(t1:t2),') while executing setup_pdbtraj().'
          deallocate(sta)
          deallocate(sto)
          call fexit()
        end if
      end do
      if (pdb_analyze.EQV..true.) then
        if (nsim.gt.c1) then
        write(ilog,*) 'Fatal. Not enough snapshots (requested ',nsim,'&
   &) in file (',fn(t1:t2),') while executing setup_pdbtraj().'
          deallocate(sta)
          deallocate(sto)
          call fexit()
        end if
      end if
    else
      if (nl.gt.(sum(natres(1:nseq))+2000)) then
        write(ilog,*) 'Fatal. Single model pdb-file too large (',&
   &fn(t1:t2),') while executing setup_pdbtraj().'
        deallocate(sta)
        deallocate(sto)
        call fexit()
      else
        sta(1) = 1
        sto(1) = nl
      end if
    end if
!   deallocate temporary arrays
    deallocate(sta)
    deallocate(sto)
  end if ! be_unsafe is false
!
! if all is well open the file for subsequent reading
  ipdbtraj = freeunit()
  open (unit=ipdbtraj,file=fn(t1:t2),status='old')
!
! and allocate memory for global arrays
  if (allocated(xpdb).EQV..true.) deallocate(xpdb)
  if (allocated(ypdb).EQV..true.) deallocate(ypdb)
  if (allocated(zpdb).EQV..true.) deallocate(zpdb)
!
end
!  
!----------------------------------------------------------------------------------
!
! closing of single-file trajectory and memory deallocation
!
subroutine close_pdbtraj()
!
  use iounit
  use pdb
  use mcsums
!
  implicit none
!
  if ((allocated(xpdb).EQV..false.).OR.&
 &    (allocated(ypdb).EQV..false.).OR.&
 &    (allocated(zpdb).EQV..false.)) then
    write(ilog,*) 'Fatal. This is a bug. PDB-coordinate arrays were &
 &not properly allocated, or prematurely freed. Please fix or repor&
 &t this problem.'
    call fexit()
  else
    deallocate(xpdb)
    deallocate(ypdb)
    deallocate(zpdb)
  end if
  close(unit=ipdbtraj)
!
end
!  
!----------------------------------------------------------------------------------
!
! if for analysis purposes we "simulate" through a xtc/dcd-trajectory, we offer
! the chance to read-in a complementary pdb-template which is set up here
!
subroutine setup_pdbtemplate()
!
  use iounit
  use atoms
  use pdb
  use mcsums
  use sequen
  use system
!
  implicit none
!
  integer t2,t1,freeunit
  logical exists
!
  ipdbtraj = freeunit()
  call strlims(pdbtmplfile,t1,t2)
  inquire (file=pdbtmplfile(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from pdb (',pdbtmplfile(t1:t2),') in setup_pdbtemplate().'
    call fexit()
  end if
!
  ipdbtraj = freeunit()
  open (unit=ipdbtraj,file=pdbtmplfile(t1:t2),status='old')
!
! and allocate memory for global map array
  allocate(pdbmap(n))
!
end
!
!-----------------------------------------------------------------------------------------------
!
! this is an auxiliary routine only used if an analysis of a trajectory in pdb-format
! uses an additional pdb-file (template for this)
!
subroutine read_pdb_template()
!
  use iounit
  use atoms
  use pdb
  use clusters
  use mcsums
!
  implicit none
!
  integer i,azero,t1,t2
!
  if (use_pdb_template.EQV..false.) return
!
  azero = 0
!
  call setup_pdbtemplate()
! init map
  do i=1,n
    pdbmap(i) = -1
  end do
! fake-read template file to get map
  call FMSMC_readpdb2(azero)
  close(unit=ipdbtraj) ! opened in setup_pdbtemplate
  call strlims(pdbtmplfile,t1,t2)
  do i=1,n
    if ((pdbmap(i).le.0).OR.(pdbmap(i).gt.n)) then
      write(ilog,*) 'Fatal. Incomplete map extracted from pdb-template (',pdbtmplfile(t1:t2),').&
& Check for unusual atom naming in pdb-file and/or setting for FMCSC_PDB_R_CONV.'
      call fexit()
    end if
  end do
! the template may serve a double fxn for extracting reference coordinates
  if (align%yes.EQV..true.) then
    call init_structalign()
  end if
!
end
!
!------------------------------------------------------------------------------------------
!
! this is the main function to process PDB files (FMCSC_PDB_READMODE is 2/3)
! it uses coordinate information to derive the complete Z matrix 
! the argument specifies the model to be used, and a value of <=0 is a special indicator
! that delineates the read-in at the beginning of a run from an analysis task
!
subroutine FMSMC_readpdb2(mdl)
!
  use iounit
  use sequen
  use polypep
  use aminos
  use atoms
  use pdb
  use mcsums
  use molecule
  use zmatrix
  use interfaces
  use system
  use clusters, ONLY: align
#ifdef ENABLE_MPI
  use mpistuff
#endif
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer ipdb,nl,i,j,k,ats,ate,freeunit,rs,rs2,t1,t3,jj,ii,kk,currs,lastrs,st1,st2,lastmol,pdb_last,iomess
  integer stat,strs,mdl,sta,sto,imol,jupp,firs,azero,rsshft,flagrsshft,unker,seqt_from_name,match(10)
  character(80) string
  character(MAXSTRLEN) pdbsnap,fn,testl
  character(MAXSTRLEN+50) fpdbstr
  character(6) keyw
  character(3) resname
  character(MAXSTRLEN), ALLOCATABLE:: pdbl(:)
  character(6), ALLOCATABLE:: kw(:)
  character(4), ALLOCATABLE:: pdbnm(:)
  character(3), ALLOCATABLE:: rsnam(:)
  character(25) atnmb
  character(1) extrabit,refbit,extrabit2,refbit2
  character, ALLOCATABLE:: chnam(:)
  integer, ALLOCATABLE:: rsnmb(:),lix(:),rsshftl(:),seqinfile(:,:)
  RTYPE refco(9),newco(9),algo(13)
  logical exists,foundit,ribon,foundp
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
  integer tslash
#endif
!
  azero = 0
 1356 format(FORM_MAXSTRLEN)
!
! if the file format is that of many independent files
  if (pdb_fileformat.eq.2) then
    if (pdb_format(3).gt.0) then
      call int2str(mdl+pdb_format(4),string,pdb_format(3))
      if ((pdb_format(1).gt.0).AND.(pdb_format(2).gt.0)) then
        pdbsnap = pdb_pref(1:pdb_format(1))//string(1:pdb_format(3))&
 &//pdb_suff(1:pdb_format(2))
      else if (pdb_format(1).gt.0) then
        pdbsnap = pdb_pref(1:pdb_format(1))//string(1:pdb_format(3))
      else if (pdb_format(2).gt.0) then
        pdbsnap = string(1:pdb_format(3))//pdb_suff(1:pdb_format(2))
      else
        pdbsnap = string(1:pdb_format(3))
      end if
    else
      call int2str(mdl+pdb_format(4),string,azero)
      call strlims(string,t1,t3)
      if ((pdb_format(1).gt.0).AND.(pdb_format(2).gt.0)) then
        pdbsnap = pdb_pref(1:pdb_format(1))//string(t1:t3)&
 &//pdb_suff(1:pdb_format(2))
      else if (pdb_format(1).gt.0) then
        pdbsnap = pdb_pref(1:pdb_format(1))//string(t1:t3)
      else if (pdb_format(2).gt.0) then
        pdbsnap = string(t1:t3)//pdb_suff(1:pdb_format(2))
      else
        pdbsnap = string(t1:t3)
      end if
    end if
    ipdb = freeunit()
    call strlims(pdbsnap,t1,t3)
    inquire (file=pdbsnap(t1:t3),exist=exists)
    if (exists.EQV..false.) then
      if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
        write(ilog,*) 'Warning. Cannot open PDB snapshot (',pdbsnap(t1:t3),').'
        mdl = -mdl
        return
      else
        write(ilog,*) 'Fatal. Cannot open PDB snapshot (',pdbsnap(t1:t3),').'
        call fexit()
      end if
    end if
    open (unit=ipdb,file=pdbsnap(t1:t3),status='old')
!
    k = 0
    do while (1.eq.1)
      k = k + 1
      read(ipdb,1356,iostat=iomess) testl
      if (iomess.eq.IOSTAT_END) then
        k = k - 1
        exit
      end if
      if ((testl(1:4).eq.'END ').OR.(testl(1:6).eq.'ENDMDL')) exit
    end do
    if (allocated(xpdb).EQV..false.) allocate(xpdb(k+10))
    if (allocated(ypdb).EQV..false.) allocate(ypdb(k+10))
    if (allocated(zpdb).EQV..false.) allocate(zpdb(k+10))
    if (size(xpdb).lt.(k+10)) then
      deallocate(xpdb)
      deallocate(ypdb)
      deallocate(zpdb)
      allocate(xpdb(k+10))
      allocate(ypdb(k+10))
      allocate(zpdb(k+10))
    end if
    allocate(pdbl(k+10))
    allocate(pdbnm(k+10))
    allocate(chnam(k+10))
    allocate(rsnam(k+10))
    allocate(kw(k+10))
    allocate(rsnmb(k+10))
    allocate(lix(k+10))
    allocate(rsshftl(k+10))
    rewind(unit=ipdb)
    do i=1,k
      read(ipdb,1356) pdbl(i)
    end do
    close(unit=ipdb)
!
    nl = k
    sta = -1
    sto = -1
    do i=1,nl
      read(pdbl(i),'(a6)') keyw
      if ((keyw.eq.'MODEL ').AND.(sta.eq.-1)) then
!        write(ilog,*) 'WARNING. PDB-File seems to contain multiple m
! &odels. Using first one only (',pdbsnap(1:t3),').'
        sta = i
      else if ((keyw.eq.'ENDMDL').AND.(sta.ne.-1)) then
        sto = i
        exit
      else if (keyw.eq.'ENDMDL') then
        if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
          write(ilog,*) 'Warning. Encountered ENDMDL before MODEL-keyword while reading from PDB (',pdbsnap(t1:t3),').'
          mdl = -mdl
          return
        else
          write(ilog,*) 'Fatal. Encountered ENDMDL before MODEL-keyword while reading from PDB (',pdbsnap(t1:t3),').'
          call fexit()
        end if
      end if
    end do
!
    ats = -1
    ate = -1
    if ((sta.ne.-1).AND.(sto.ne.-1)) then
      do i=sta,sto
        read(pdbl(i),'(a6)') keyw
        if (((keyw.eq.'ATOM  ').OR.&
 &           (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
          ats = i
        else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &           (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
          ate = i-1
          exit
        end if
      end do
    else if ((sta.eq.-1).AND.(sto.eq.-1)) then
      do i=1,nl
        read(pdbl(i),'(a6)') keyw
        if (((keyw.eq.'ATOM  ').OR.&
 &           (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
          ats = i
        else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &           (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
          ate = i-1
          exit
        end if
      end do
      if ((ats.gt.0).AND.(ate.eq.-1)) then
        ate = nl
      end if
    else
      if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
        write(ilog,*) 'Warning. Could not find ENDMDL keyword despite MODEL keyword while reading from PDB (',pdbsnap(t1:t3),').'
        mdl = -mdl
        return
      else
        write(ilog,*) 'Fatal. Could not find ENDMDL keyword despite MODEL keyword while reading from PDB (',pdbsnap(t1:t3),').'
        call fexit()
      end if
    end if
    if ((ats.le.0).AND.(ate.le.0)) then
      if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
        write(ilog,*) 'Warning. PDB file (',pdbsnap(t1:t3),') appears to contain no coordinate information whatsoever. &
 &Please check input.'
        mdl = -mdl
        return
      else
        write(ilog,*) 'Fatal. PDB file (',pdbsnap(t1:t3),') appears to contain no coordinate information whatsoever. &
 &Please check input.'
        call fexit()
      end if
    end if

    fn(t1:t3) = pdbsnap(t1:t3)
!
! else the pdb-file is a large trajectory file and we'll just keep reading it
  else
!
    call strlims(pdbinfile,t1,t3)
#ifdef ENABLE_MPI
    if (((pdb_mpimany.EQV..true.).OR.(pdb_analyze.EQV..true.)).AND.(mdl.ge.1)) then
      if ((pdb_analyze.EQV..false.).AND.(fycxyz.eq.1)) then
        write(ilog,*) 'Warning. Reading multiple pdb-files for parallel simulations can be problematic if the coordinates are &
 &of low precision. This is because covalent geometries will differ slightly between replicas. It is recommended to &
 &use FMCSC_PDB_INPUTSTRING (and FMCSC_PDB_OUTPUTSTRING) to circumvent this issue whenever possible.'
      end if
      do i=t3,t1,-1
        if (pdbinfile(i:i).eq.SLASHCHAR) exit
      end do
      tslash = i
      call int2str(myrank,xpont,re_aux(10))
      if ((tslash.ge.t1).AND.(tslash.lt.t3)) then
        fn =  pdbinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//pdbinfile((tslash+1):t3)
      else
        fn =  'N_'//xpont(1:re_aux(10))//'_'//pdbinfile(t1:t3)
      end if
      call strlims(fn,t1,t3)
    else
      fn(t1:t3) = pdbinfile(t1:t3)
    end if
#else
    fn(t1:t3) = pdbinfile(t1:t3)
#endif
    if (mdl.le.0) then ! flag for template
      call strlims(pdbtmplfile,t1,t3)
      fn(t1:t3) = pdbtmplfile(t1:t3)
    end if
    sta = -1
    sto = -1    
    k = 0
    do while (1.eq.1)
      k = k + 1
      read(ipdbtraj,1356,iostat=iomess) testl
      if (iomess.eq.IOSTAT_END) then
        k = k
        exit
      end if
      if ((testl(1:4).eq.'END ').OR.(testl(1:6).eq.'ENDMDL')) exit
    end do
    do i=1,k
      backspace(unit=ipdbtraj)
    end do
    if (iomess.eq.IOSTAT_END) then
      k = k - 1
    end if
    allocate(pdbl(k+10))
    allocate(pdbnm(k+10))
    allocate(chnam(k+10))
    allocate(rsnam(k+10))
    allocate(kw(k+10))
    allocate(rsnmb(k+10))
    allocate(lix(k+10))
    allocate(rsshftl(k+10))
    if (allocated(xpdb).EQV..false.) allocate(xpdb(k+10))
    if (allocated(ypdb).EQV..false.) allocate(ypdb(k+10))
    if (allocated(zpdb).EQV..false.) allocate(zpdb(k+10))
    if (size(xpdb).lt.(k+10)) then
      deallocate(xpdb)
      deallocate(ypdb)
      deallocate(zpdb)
      allocate(xpdb(k+10))
      allocate(ypdb(k+10))
      allocate(zpdb(k+10))
    end if
    do i=1,k
      read(ipdbtraj,1356) pdbl(i)
      read(pdbl(i),'(a6)') keyw
      if (keyw.eq.'MODEL ') then
        sta = i
      else if (keyw.eq.'ENDMDL') then
        sto = i
        exit
      else if (keyw(1:4).eq.'END ') then
        if (pdb_analyze.EQV..true.) then
          write(ilog,*) 'Warning. Discovered END keyword in PDB file ',fn(t1:t3),' when expecting only MODEL/ENDMDL syntax.'
        end if
        sto = i
        exit
      end if 
    end do
!
    nl = k
    if (sto.eq.-1) then
      if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then ! should always trigger the first
        write(ilog,*) 'Warning.  Not enough snapshots (requested ',mdl,') in file (',fn(t1:t3),') while executing &
 &FMSMC_readpdb2().'
        mdl = -mdl
        return
      else if ((mdl.eq.1).AND.(pdb_analyze.EQV..true.)) then
        write(ilog,*) 'Fatal. Trajectory file (',fn(t1:t3),') contains no usable models or is unterminated (END/ENDMDL).'
        call fexit()
      else
        write(ilog,*) 'Fatal. PDB input from file (',fn(t1:t3),') cannot be used (is empty or unterminated (END/ENDMDL)).'
      end if
    end if
!   warn if MODEL is missing
    if ((sta.eq.-1).AND.(pdb_analyze.EQV..true.)) then
      write(ilog,*) 'Warning. When trying to read the ',mdl,'th snapshot from ',fn(t1:t3),', no MODEL keyword was found yet &
 &an END or ENDMDL was present.'
    end if
!
    ats = -1
    ate = -1
    do i=1,nl
      read(pdbl(i),'(a6)') keyw
      if (((keyw.eq.'ATOM  ').OR.&
 &         (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
        ats = i
      else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &         (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
        ate = i-1
        exit
      end if
    end do
    if ((ats.gt.0).AND.(ate.eq.-1)) then
      ate = n
    end if
    if ((ats.le.0).AND.(ate.le.0)) then
      if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
        write(ilog,*) 'Warning. Snapshot read from PDB file (',fn(t1:t3),') is empty or absent. &
 &Please check input and setting for FMCSC_NRSTEPS.'
        mdl = -mdl
        return
      else
        write(ilog,*) 'Fatal. PDB file (',fn(t1:t3),') appears to contain no coordinate information whatsoever. &
 &Please check input.'
        call fexit()
      end if
    end if
!
  end if
!
!
  call strlims(pdb_formstr(2),st1,st2)
  fpdbstr = '('//pdb_formstr(2)(st1:st2)//')'
!
  k = 0
  rsshft = 0
  rsshftl(:) = 0
  currs = huge(currs)
  refbit = ' '
  refbit2 = ' '
  do i=ats,ate
    k = k + 1
    read(pdbl(i),'(a6)') keyw
    if ((keyw.eq.'TER   ').OR.(keyw.eq.'ANISOU')) then
      k = k - 1
      cycle
    end if
    read(pdbl(i),fpdbstr) kw(k),atnmb,pdbnm(k),extrabit,rsnam(k),chnam(k),rsnmb(k),extrabit2&
 &,xpdb(k),ypdb(k),zpdb(k)
    lix(k) = i
!   the currs/lastrs construct is (at least sporadically) able to read in pdb-files in which
!   the residue numbering is off as long as it is a number and either the number for a new residue is always
!   different from the previous one or the extrabit2 (insertion code) differs
    if (k.eq.1) lastrs = rsnmb(k)-1
    if (currs.ne.huge(currs)) then
      if ((rsnmb(k).ne.currs).OR.(extrabit2.ne.refbit2)) then
        refbit = ' '
        if (currs.ne.lastrs+1) currs = lastrs+1
        lastrs = currs
        currs = rsnmb(k)
        if (rsnmb(k).ne.lastrs+1) then
          rsnmb(k) = lastrs + 1
        end if
      else
        if (rsnmb(k).ne.lastrs+1) rsnmb(k) = lastrs+1
      end if
    end if
    if (k.eq.1) currs = rsnmb(k)
    if ((refbit.eq.' ').AND.(extrabit.ne.' ')) refbit = extrabit
    refbit2 = extrabit2
    call toupper(kw(k))
!   drop the records for TER. ANISOU, and any alternate conformations not found first
    if ((kw(k).eq.'TER   ').OR.(kw(k).eq.'ANISOU').OR.((extrabit.ne.' ').AND.(extrabit.ne.refbit))) then
      k = k - 1
      cycle
    end if
    call toupper(pdbnm(k))
    call toupper(rsnam(k))
    call toupper(chnam(k))
    call pdbcorrect(rsnam(k),pdbnm(k),rsnmb(k),rsshft,rsshftl)
  end do
!
  ii = 1
  do j=1,rsshft
    flagrsshft = 0
    do i=ii,k
      if (rsnmb(i).eq.rsshftl(j)) then
        do jj=i,k
          if (rsnmb(jj).ne.rsshftl(j)) then
            ii = jj
            exit
          end if
          if ((rsnam(jj).eq.'NH2').OR.(rsnam(jj).eq.'NME')) then
            rsnmb(jj) = rsnmb(jj) + 1
            if (flagrsshft.eq.0) flagrsshft = 1
            if (flagrsshft.eq.-1) flagrsshft = 2
          end if
          if ((rsnam(jj).eq.'FOR').OR.(rsnam(jj).eq.'ACE')) then
            rsnmb(jj) = rsnmb(jj) - 1
            if (flagrsshft.eq.0) flagrsshft = -1
            if (flagrsshft.eq.1) flagrsshft = 2
          end if
          if (jj.eq.k) ii = k + 1
        end do
        if (flagrsshft.eq.1) then
          if (ii.le.k) then
            rsnmb(ii:k) = rsnmb(ii:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 1
          end if
        else if (flagrsshft.eq.-1) then
          if (i.le.k) then
            rsnmb(i:k) = rsnmb(i:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 1
          end if
        else if (flagrsshft.eq.2) then
          if (i.le.k) then
            rsnmb(i:k) = rsnmb(i:k) + 1
          end if
          if (ii.le.k) then
            rsnmb(ii:k) = rsnmb(ii:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 2
          end if
        end if
        exit ! skip out, go to next shft
      end if
    end do
  end do
!
  jupp = -1000
  do i=1,k
    if (rsnmb(i).ne.jupp) then
      if (i.gt.1) then
        ribon = .false.
        foundp = .false.
        do j=firs,i-1
          if (pdbnm(j).eq.' P  ') then ! crucial that P is always called P
            foundp = .true.
            exit
          end if
        end do
        do j=firs,i-1
          if (rsnam(j)(1:2).eq.'R ') then
            ribon = .true.
          end if
        end do
        if (ribon.EQV..true.) then
          do j=firs,i-1
!           hopeful
            if (foundp.EQV..true.) then
              rsnam(j)(1:2) = 'RP'
            else
              rsnam(j)(1:2) = 'RI'
            end if
            if (pdb_convention(2).eq.3) then
              if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = '2HO*'
              if (pdbnm(j)(1:4).eq.'2H2*') pdbnm(j)(1:4) = ' H2*'
            else
              if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = ' H2*'
            end if
          end do
        else if (rsnam(firs)(1:2).eq.'DP') then
          if (foundp.EQV..false.) then
            do j=firs,i-1
              rsnam(j)(2:2) = 'I'
            end do
          end if
        end if
        if ((ribon.EQV..true.).OR.(rsnam(firs)(1:2).eq.'DP')) then
          if (foundp.EQV..true.) then
            do j=firs,i-1          
              if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = '5HO*'
            end do
          else
            do j=firs,i-1  
              if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = ' HOP'
            end do
          end if
        end if
      end if
      firs = i
      jupp = rsnmb(i)
    end if
    if (i.eq.k) then
      ribon = .false.
      foundp = .false.
      do j=firs,k
        if (pdbnm(j).eq.' P  ') then ! crucial that P is always called P
          foundp = .true.
          exit
        end if
      end do
      do j=firs,k
        if (rsnam(j)(1:2).eq.'R ') then
          ribon = .true.
        end if
      end do
      if (ribon.EQV..true.) then
        do j=firs,k
          if (foundp.EQV..true.) then
            rsnam(j)(1:2) = 'RP'
          else
            rsnam(j)(1:2) = 'RI'
          end if
          if (pdb_convention(2).eq.3) then
            if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = '2HO*'
            if (pdbnm(j)(1:4).eq.'2H2*') pdbnm(j)(1:4) = ' H2*'
          else
            if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = ' H2*'
          end if
        end do
      else if (rsnam(firs)(1:2).eq.'DP') then
        if (foundp.EQV..false.) then
          do j=firs,k
            rsnam(j)(2:2) = 'I'
          end do
        end if
      end if
      if ((ribon.EQV..true.).OR.(rsnam(firs)(1:2).eq.'DP')) then
        if (foundp.EQV..true.) then
          do j=firs,i-1          
            if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = '5HO*'
          end do
        else
          do j=firs,i-1  
            if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = ' HOP'
          end do
        end if
      end if
    end if
  end do
!
  nl = k
!
  call adjust_PONuc(nl,rsnmb(1:nl),pdbnm(1:nl),rsnam(1:nl),xpdb(1:nl),ypdb(1:nl),zpdb(1:nl))
!
! this function reorders the input arrays (except pdbl) so that everything is contiguous per residue
  call reorder_PDB_byrsnmb(nl,kw(1:nl),rsnmb(1:nl),rsnam(1:nl),pdbnm(1:nl),xpdb(1:nl),ypdb(1:nl),zpdb(1:nl),&
 &lix(1:nl),chnam(1:nl))
!
  if (allocated(pdbmap2).EQV..false.) allocate(pdbmap2(n))
  pdbmap2(:) = 0
!
  if (pdb_readmode.eq.3) then ! structural input is sparse (parts are potentially missing)
!
!   now let's locate the sequence
    kk = rsnmb(nl)-rsnmb(1)+1
    allocate(seqinfile(kk,6))
    seqinfile(:,:) = 0
    seqinfile(1,2) = 1
    j = 0
    jj = 1
    do i=2,nl
      if (rsnmb(i).ne.rsnmb(i-1)) then
        j = j + 1
        seqinfile(j,1) = seqt_from_name(rsnam(i-1),jj)
        seqinfile(j+1,2) = i
      end if
    end do
    j = j + 1
    seqinfile(j,1) = seqt_from_name(rsnam(nl),jj)
    do unker=1,n_pdbunk
      rs = pdb_unkbnd(unker,3)
      seqtyp(rs) = -unker ! so specific comparison works
    end do
!
    stat = 1
    do imol=1,nmol
      lastmol = imol
      match(:) = 0
      rs= rsmol(imol,1)
      do while (rs.le.rsmol(imol,2))
        if ((seqtyp(rs).eq.seqinfile(stat,1)).OR.&
 &      ((seqinfile(stat,1).eq.99999).AND.(seqtyp(rs).ge.13).AND.(seqtyp(rs).le.15))) then
          if (rsmol(imol,1).eq.rsmol(imol,2)) then ! single-res. molecule
            match(1) = rs
            match(2) = 1
            match(3) = stat
            exit
          end if
          j = rs+1
          jj = stat+1
          if ((j.gt.rsmol(imol,2)).OR.(jj.gt.kk)) then
!           we have no chance to extend the molecule
            match(1) = rs
            match(2) = 1
            match(3) = stat
            exit
          end if
          do while ((j.le.rsmol(imol,2)).AND.(jj.le.kk))
            if (((seqtyp(j).ne.seqinfile(jj,1)).AND.&
 &               (.NOT.((seqinfile(jj,1).eq.99999).AND.(seqtyp(j).ge.13).AND.(seqtyp(j).le.15)))).OR.&
                (j.eq.rsmol(imol,2)).OR.(jj.eq.kk)) then
              if ((j-rs).gt.match(2)) then
                match(1) = rs
                if ((seqtyp(j).ne.seqinfile(jj,1)).AND.&
 &               (.NOT.((seqinfile(jj,1).eq.99999).AND.(seqtyp(j).ge.13).AND.(seqtyp(j).le.15)))) then
                  match(2) = j-rs
                else
                  match(2) = j-rs+1
                end if
                match(3) = stat
              end if
              exit
            end if
            j = j + 1
            jj = jj + 1
          end do
        end if
        if (match(2).ge.(rsmol(imol,2)-rs)) exit ! if the matched segment is n-1 or n, it cannot be improved
        rs = rs + 1
      end do
      if ((match(1).le.0).AND.(imol.gt.1)) then
        write(ilog,*) 'Warning. Cannot find any candidate structural information for molecule ',&
 &imol,'. All structural input is stopped at ',seqinfile(imol-1,4),'th ATOM/HETATM record in input PDB-file (',fn(t1:t3),'). &
 &Molecules with no structural information present should always appear at the end of sequence input.'
        lastmol = imol-1
        exit
      else if (match(1).le.0) then
        write(ilog,*) 'Fatal. Cannot find any structural information for first molecule. This can indicate superfluous atom &
 &records at the beginning of the PDB file (in FMCSC_PDB_READMODE 3, all superfluous information should be discarded), a &
 &complete mismatch in sequences (e.g., by the incorrect placement of a molecule to be built randomly at the beginning of &
 & sequence input), or an empty or corrupt input file.'
        call fexit()
      else
        seqinfile(imol,3) = seqinfile(match(3),2)
        if ((match(3)+match(2)).gt.size(seqinfile,dim=1)) then
          seqinfile(imol,4) = nl
        else
          seqinfile(imol,4) = seqinfile(match(3)+match(2),2)-1
        end if
        seqinfile(imol,5) = match(1)
        seqinfile(imol,6) = match(1)+match(2)-1
        if (match(2).ne.(rsmol(imol,2)-rsmol(imol,1)+1)) then
 456 format(1x,'Molecule ',i7,' (res. ',i7,'-',i7,') was partially matched &
 &(res. ',i7,'-',i7,'), and the following PDB records enclose the matched &
 &region:'/'====================================================================='/(3x,a)/(3x,a)/&
 &         '=====================================================================')
          write(ilog,456) imol,rsmol(imol,1),rsmol(imol,2),match(1),match(1)+match(2)-1,&
 &trim(pdbl(lix(seqinfile(imol,3)))),trim(pdbl(lix(seqinfile(imol,4))))
        end if
        if ((seqinfile(imol,4).eq.nl).AND.(imol.lt.nmol)) then
          write(ilog,*) 'All structural input is exhausted in reading from input PDB-file (',fn(t1:t3),'). Molecule(s) &
 &',imol+1,' to ',nmol,' will be constructed entirely via fallback option for FMCSC_RANDOMIZE.'
        end if
        pdb_last = match(1)+match(2)-1
        stat = stat + match(2)
      end if
!
      rs = match(1)
      rs2 = match(3)
      if (seqtyp(rs).gt.0) then
        resname = amino(seqtyp(rs))
      else
        resname = 'UNK'
      end if
      if (match(1).ne.rsmol(imol,1)) then
        if (seqtyp(match(1)).lt.0) then
          do ii=1,3
            refco(3*ii-2) = x(at(rs)%bb(1)+ii-1)
            refco(3*ii-1) = y(at(rs)%bb(1)+ii-1)
            refco(3*ii) = z(at(rs)%bb(1)+ii-1)
          end do
        else
          do ii=1,3
            refco(3*ii-2) = x(at(rs)%bb(ii))
            refco(3*ii-1) = y(at(rs)%bb(ii))
            refco(3*ii) = z(at(rs)%bb(ii))
          end do
        end if
      end if
      j = -1
      do i=seqinfile(imol,3),seqinfile(imol,4)
        if (i.gt.seqinfile(imol,3)) then
          if (rsnmb(i).ne.rsnmb(i-1)) then
            rs = rs + 1
            if (seqtyp(rs).gt.0) then
              resname = amino(seqtyp(rs))
            else
              resname = 'UNK'
            end if
            if (j.gt.0) then
              if ((j+1).ne.at(rs-1)%nbb+at(rs-1)%nsc) then
                write(ilog,*) 'Fatal. Structural information for unsupported residue with name ',rsnam(i-1),' read from &
 &file (FMCSC_PDBFILE: ',fn(t1:t3),') is inconsistent with that inferred from template (FMCSC_PDB_TEMPLATE).'
                call fexit()
              end if
            end if
            j = -1
          end if
        end if
!
        if (seqtyp(rs).lt.0) then
          j = j + 1
          call trfxyz(i,at(rs)%bb(1)+j)
          cycle
        end if
!
!       the transfer wrapper functions are pretty simple (and tedious) scanning lists
!       atoms have previously been reordered to allow sequential processing (in residue order)
!       even if the original pdb-input list is out of order (CHARMM caps for instance)
!       this will also set pdbmap correctly 
!
!       single residues molecules first
        if (rsmol(molofrs(rs),1).eq.rsmol(molofrs(rs),2)) then
!
          call trf_sms(i,rs,resname,pdbnm(i))
!
!       now the backbone of N-terminal residues
        else if (rs.eq.rsmol(molofrs(rs),1)) then
!
          call trf_nts(i,rs,resname,pdbnm(i))
!
!       now the backbone of C-terminal residues
        else if (rs.eq.rsmol(molofrs(rs),2)) then
!
          call trf_cts(i,rs,resname,pdbnm(i))
!
!       finally, regular in-the-middle-of-the-chain residues
        else
!
          call trf_bbs(i,rs,resname,pdbnm(i))
!
        end if
!
!       now the sidechains
!
        call trf_scs(i,rs,resname,pdbnm(i))

      end do
      if (match(1).ne.rsmol(imol,1)) then
        rs = match(1)
!       now we try to move the pre-built piece over by alignment
!       for this, one should ideally search for consecutively bonded atoms that were successfully read-in
!       since the building will create unwanted results in any case if backbone atoms are missing, we just
!       use the first 3 of those (there always must be at least 3 in a polymer residue!)
        if (seqtyp(match(1)).lt.0) then
          do ii=1,3
            newco(3*ii-2) = x(at(rs)%bb(1)+ii-1)
            newco(3*ii-1) = y(at(rs)%bb(1)+ii-1)
            newco(3*ii) = z(at(rs)%bb(1)+ii-1)
          end do
        else
          do ii=1,3
            newco(3*ii-2) = x(at(rs)%bb(ii))
            newco(3*ii-1) = y(at(rs)%bb(ii))
            newco(3*ii) = z(at(rs)%bb(ii))
          end do
        end if
        ii = 3
        call align_3D(ii,refco,newco,algo(1:3),algo(4:7),algo(8:10),algo(11:13))
        do j=atmol(imol,1),at(rs)%bb(1)-1
          call quat_conjugate4(algo(4:7),j,algo(8:10))
        end do
        x(atmol(imol,1):(at(rs)%bb(1)-1)) = x(atmol(imol,1):(at(rs)%bb(1)-1)) + algo(1)
        y(atmol(imol,1):(at(rs)%bb(1)-1)) = y(atmol(imol,1):(at(rs)%bb(1)-1)) + algo(2)
        z(atmol(imol,1):(at(rs)%bb(1)-1)) = z(atmol(imol,1):(at(rs)%bb(1)-1)) + algo(3)
      end if
      if (seqinfile(imol,4).eq.nl) exit
    end do
    if (allocated(pdb_didread).EQV..true.) deallocate(pdb_didread)
    allocate(pdb_didread(lastmol,3))
    pdb_didread(1:lastmol,1) = seqinfile(1:lastmol,5)
    pdb_didread(1:lastmol,2) = seqinfile(1:lastmol,6)
    deallocate(seqinfile)
    do unker=1,n_pdbunk
      rs = pdb_unkbnd(unker,3)
      seqtyp(rs) = 26 ! restore
    end do
!
  else if (pdb_readmode.eq.2) then
!
!   now let's check whether the sequence matches
    rs = 1
    unker = 0
    if ((seqtyp(rs).eq.26).OR.(seqtyp(rs).le.0)) then
      unker = unker + 1
      resname = pdb_unknowns(unker)
    else
      resname = amino(seqtyp(rs))
    end if
    i = 1
    do while ((resname.ne.rsnam(i)).AND.&
 & (.NOT.((rsnam(i).eq.'HIS').AND.((resname.eq.'HID').OR.(resname.eq.'HIE').OR.(resname.eq.'HIP')))))
      i = i+1
      if (i.gt.nl) then
        if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
          write(ilog,*) 'Warning. Cannot find first residue (',resname,') in PDB-file (',fn(t1:t3),').'
          mdl = -mdl
          return
        else
          write(ilog,*) 'Fatal. Cannot find first residue (',resname,') in PDB-file (',fn(t1:t3),').'
          call fexit()
        end if
      end if
    end do
    if (rsnam(i).eq.'HIS') then
      do j=i+1,nl
        if (j.gt.1) then
          if (rsnmb(j).ne.rsnmb(j-1)) exit
        end if
        rsnam(j) = resname
      end do
    end if
    strs = rsnmb(i)
    rs2 = strs
    stat = i
    call strlims(seqfile,st1,st2)
    pdb_last = nseq
    do while (rs.lt.nseq)
      rs = rs + 1
      rs2 = rs2 + 1
      if ((seqtyp(rs).eq.26).OR.(seqtyp(rs).le.0)) then
        unker = unker + 1
        resname = pdb_unknowns(unker)
      else
        resname = amino(seqtyp(rs))
      end if
      foundit = .false.
      do i=stat,nl
        if (rsnmb(i).eq.rs2) then
          if ((rsnam(i).eq.resname).OR.&
 &  ((rsnam(i).eq.'HIS').AND.((resname.eq.'HID').OR.(resname.eq.'HIE').OR.(resname.eq.'HIP')))) then
            foundit = .true.
            if (rsnam(i).eq.'HIS') then
              do j=i+1,nl
                if (j.gt.1) then
                  if (rsnmb(j).ne.rsnmb(j-1)) exit
                end if
                rsnam(j) = resname
              end do
            end if
            exit
          else
            if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
              write(ilog,*) 'Warning. Sequences in sequence file and PDB-file do not match.'
              write(ilog,*) 'Seq. file: ',seqfile(st1:st2),', PDB-file: ',fn(t1:t3)
              write(ilog,*) 'Position ',rs,' in seq. file (',resname,') and position ',&
   &                       rsnmb(i),' in PDB-file (',rsnam(i),')!'
              mdl = -mdl
              return
            else
              write(ilog,*) 'Fatal. Sequences in sequence file and PDB-file do not match.'
              write(ilog,*) 'Seq. file: ',seqfile(st1:st2),', PDB-file: ',fn(t1:t3)
              write(ilog,*) 'Position ',rs,' in seq. file (',resname,') and position ',&
   &                       rsnmb(i),' in PDB-file (',rsnam(i),')!'
              call fexit()
            end if
          end if
        end if
      end do
      if (foundit.EQV..false.) then
        if ((pdb_analyze.EQV..true.).OR.(fn(t1:t3).eq.pdbtmplfile(t1:t3)).OR.(rs.le.1)) then
          if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
            write(ilog,*) 'Warning. Residue numbering in input PDB-file ',fn(t1:t3),' is d&
   &iscontinuous or residues are missing. Please fix this first.'
            mdl = -mdl
            return
          else
            write(ilog,*) 'Fatal. Residue numbering in input PDB-file ',fn(t1:t3),' is d&
   &iscontinuous or residues are missing. Please fix this first.'
            call fexit()
          end if
        end if
        write(ilog,*) 'Warning. Residue numbering in input PDB-file is d&
   &iscontinuous or residues are missing. Stopping all read-in of structural information at &
   &residue ',rs,'. If possible, missing information will be generated via settings for FMCSC_RANDOMIZE.'
        pdb_last = rs - 1
        exit
      end if
    end do
!
!   if we're still here, we can finally proceed to extract the xyz-coordinates
!
    rs = 1
    rs2 = strs
    i = at(pdb_last)%bb(1)+at(pdb_last)%nbb+at(pdb_last)%nsc-1
    x(1:i) = 0.0
    y(1:i) = 0.0
    z(1:i) = 0.0
    do while (rs.le.pdb_last)
      j = -1
      resname = amino(seqtyp(rs))
      do i=stat,nl
!
        if (rsnmb(i).ne.rs2) then
          stat = i
          exit
        end if
!
        if (i.gt.stat) then
          if (rsnmb(i).ne.rsnmb(i-1)) then
            if (j.gt.0) then
              if (j.ne.at(rs-1)%nbb+at(rs-1)%nsc) then
                if ((mdl.gt.1).AND.(pdb_analyze.EQV..true.)) then
                  write(ilog,*) 'Warning. Structural information for unsupported residue with name ',rsnam(i-1),' read from &
 &file (FMCSC_PDBFILE: ',fn(t1:t3),') is inconsistent with that inferred from template (FMCSC_PDB_TEMPLATE).'
                  mdl = -mdl
                  return
                else
                  write(ilog,*) 'Fatal. Structural information for unsupported residue with name ',rsnam(i-1),' read from &
 &file (FMCSC_PDBFILE: ',fn(t1:t3),') is inconsistent with that inferred from template (FMCSC_PDB_TEMPLATE).'
                  call fexit()
                end if
              end if
            end if
          end if
        end if
!
        if (seqtyp(rs).eq.26) then
          j = j + 1
          call trfxyz(i,at(rs)%bb(1)+j)
          cycle
        end if
!
!       the transfer wrapper functions are pretty simple (and tedious) scanning lists
!       atoms have previously been reordered to allow sequential processing (in residue order)
!       even if the original pdb-input list is out of order (CHARMM caps for instance)
!       this will also set pdbmap correctly 
!
!       single residues molecules first
        if (rsmol(molofrs(rs),1).eq.rsmol(molofrs(rs),2)) then
!
          call trf_sms(i,rs,resname,pdbnm(i))
!
!       now the backbone of N-terminal residues
        else if (rs.eq.rsmol(molofrs(rs),1)) then
!
          call trf_nts(i,rs,resname,pdbnm(i))
!
!       now the backbone of C-terminal residues
        else if (rs.eq.rsmol(molofrs(rs),2)) then
!
          call trf_cts(i,rs,resname,pdbnm(i))
!
!       finally, regular in-the-middle-of-the-chain residues
        else
!
          call trf_bbs(i,rs,resname,pdbnm(i))
!
        end if
!
!       now the sidechains
!
        call trf_scs(i,rs,resname,pdbnm(i))
!
!
      end do
      rs = rs + 1
      rs2 = rs2 + 1
    end do
    imol = molofrs(rs-1)
    if (allocated(pdb_didread).EQV..true.) deallocate(pdb_didread)
    allocate(pdb_didread(imol,3))
    do i=1,imol-1
      pdb_didread(i,1:2) = rsmol(i,1:2)
    end do
    pdb_didread(imol,1) = rsmol(i,1)
    pdb_didread(imol,2) = pdb_last
!
  end if ! else if (pdb_readmode.eq.2)
!
! generate appropriate internal coordinates from the read-in Cartesian values (unless we are
! just reading the template)
  if (mdl.gt.0) then
    k = 0
    kk = 0
    jj = size(pdb_didread,dim=1)
    do imol=1,nmol
      if (((pdbmap2(atmol(imol,1)).le.0).OR.(pdbmap2(min(atmol(imol,2),atmol(imol,1)+1)).le.0).OR.&
 &         (pdbmap2(min(atmol(imol,2),atmol(imol,1)+2)).le.0)).AND.&
 &        ((pdbmap2(atmol(imol,1)).gt.0).OR.(pdbmap2(min(atmol(imol,2),atmol(imol,1)+1)).gt.0).OR.&
 &         (pdbmap2(min(atmol(imol,2),atmol(imol,1)+2)).gt.0))) then
!       atoms from reference frame are missing
        call makexyz_formol4(imol,(pdbmap2(atmol(imol,1)).le.0),&
 &                                (pdbmap2(min(atmol(imol,2),atmol(imol,1)+1)).le.0),&
 &                                (pdbmap2(min(atmol(imol,2),atmol(imol,1)+2)).le.0))
        if ((rsmol(imol,1).ne.rsmol(imol,2)).AND.(mdl.eq.1)) then
          if (pdb_analyze.EQV..false.) then
            write(ilog,'(1x,a,i8,a)') 'Warning. Only 1 or 2 atoms of the first 3 could be read for molecule ',imol,', which &
 &is unlikely to create a satisfactory starting structure. Check for naming problems or missing atoms. In case of the latter, &
 &it may be better to delete all information pertaining to the incomplete residue.'
          else
            write(ilog,*) 'Warning. Only 1 or 2 atoms of the first 3 could be read for molecule ',imol,'. Check &
 &for naming problems or missing atoms. Any further warnings of this type are omitted.'
          end if
        else if (atmol(imol,2).gt.atmol(imol,1)+1) then
          k = k + 1
        end if
      else if ((pdbmap2(atmol(imol,1)).le.0).AND.(pdbmap2(min(atmol(imol,2),atmol(imol,1)+1)).le.0).AND.&
 &         (pdbmap2(min(atmol(imol,2),atmol(imol,1)+2)).le.0)) then
        if ((mdl.eq.1).AND.(imol.le.jj).AND.(kk.eq.0)) then
          if (imol.gt.size(pdb_didread,dim=1)) then
            write(ilog,'(1x,a,i8,a)') 'Warning. Reference frame could not be read for molecule ',imol,', which &
 &will create arbitrary absolute coordinates.'
            kk = kk + 1
          else if (pdb_didread(imol,2).eq.rsmol(imol,1)) then
            write(ilog,'(1x,a,i8,a)') 'Warning. Reference frame could not be read for molecule ',imol,', which &
 &will likely create arbitrary absolute coordinates.'
            kk = kk + 1
          end if
          if (kk.eq.1) then
            if (pdb_analyze.EQV..false.) then
              write(ilog,*) 'If this is unexpected, check for naming problems or missing atoms. &
 &Any further warnings of this type are omitted.'
            else
              write(ilog,*) 'Check sequence input. Any further warnings of this type are omitted.'
            end if
          end if
        end if
      end if
      call genzmat_f(imol,mdl)
    end do
    if ((k.gt.0).AND.(mdl.eq.1)) then
      write(ilog,'(1x,a)') 'Warning. At least one single-residue molecule was processed from structural input with &
 &partial coordinates for its reference frame (first 3 atoms). This can shift input coordinates in these molecules.'
      if (pdb_analyze.EQV..true.) write(ilog,*) 'Any further warnings of this type are omitted.'
    end if
!     
!   and transfer them into FMSMC-arrays (which are the ones used in analysis functions)
    call zmatfyc()
!
!   and re-build xyz so that hydrogen corrections are accounted for
    do imol=1,nmol
      call makexyz_formol(imol)
    end do
!
    do imol=1,nmol
      call update_rigid(imol)
      call update_rigidm(imol)
    end do
    if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
      if (align%yes.EQV..true.) then
        do imol=1,nmol
          call update_image2(imol)
        end do
      else
        do imol=1,nmol
          call update_image(imol)
        end do
      end if
    end if
  end if
!
  deallocate(pdbl)
  deallocate(pdbnm)
  deallocate(chnam)
  deallocate(rsnam)
  deallocate(kw)
  deallocate(rsnmb)
  deallocate(lix)
  if (mdl.le.0) write(ilog,*) 
!
end
!
!-----------------------------------------------------------------------
!
subroutine trf_sms(i,rs,resname,pdbnm)
!
  use pdb
  use iounit
  use polypep
  use sequen
  use aminos
  use system
  use molecule
!
  implicit none
!
  integer i,rs,shf,shf2,shf3,shf4,shfx1,shfx2
  character(3) resname
  character(4) pdbnm
  character(1) poltyc
!
  shf = 0
  shf2 = 0
  shf3 = 0
  shfx1 = 0
  shfx2 = 0
  shf4 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
    shf4 = 4
    if (ua_model.eq.2) then
      shfx1 = 1
      shfx2 = 2
    end if
  end if
  if (seqtyp(rs).lt.0) then
    poltyc = aminopolty(seqtyp(26))
  else
    poltyc = aminopolty(seqtyp(rs))
  end if 
!
  if (poltyc.eq.'P') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.'1OXT') call trfxyz(i,oi(rs))
    if (pdbnm.eq.'2OXT') call trfxyz(i,at(rs)%bb(5))
!   PRO, HYP 
    if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.32)) then
      if (moltermid(molofrs(rs),1).eq.1) then
        if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(6))
        if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(7))
      else if (moltermid(molofrs(rs),1).eq.2) then
        if (pdbnm.eq.' H  ') call trfxyz(i,at(rs)%bb(6))
      else
        write(ilog,*) 'Unsupported N-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
      end if
      if (ua_model.eq.0) then
        if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
      end if
    else
      if (moltermid(molofrs(rs),1).eq.1) then
        if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(6))
        if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(7))
        if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%bb(8))
      else if (moltermid(molofrs(rs),1).eq.2) then
        if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(6))
        if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(7))
       else
        write(ilog,*) 'Unsupported N-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
      end if
      if (ua_model.eq.0) then
        if (resname.eq.'GLY') then
          if (pdbnm.eq.'1HA ') call trfxyz(i,at(rs)%sc(1))
          if (pdbnm.eq.'2HA ') call trfxyz(i,at(rs)%sc(2))
        else if (resname.ne.'AIB') then
          if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
        end if
      end if
    end if
  else if (poltyc.eq.'N') then
    if ((resname.eq.'D5P').OR.(resname.eq.'DPC').OR.&
   &       (resname.eq.'DPU').OR.(resname.eq.'DPT').OR.&
   &       (resname.eq.'DPA').OR.(resname.eq.'DPG')) then
      if (pdbnm.eq.' O3*') call trfxyz(i,nuci(rs,1))
      if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
      if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
      if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
      if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
      if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
      if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2O3*') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(10))
        if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(11))
        if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(12))
        if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(13))
        if (pdbnm.eq.'1H2*') call trfxyz(i,at(rs)%sc(4))
        if (pdbnm.eq.'2H2*') call trfxyz(i,at(rs)%sc(5))
        if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
      end if
      if (moltermid(molofrs(rs),1).eq.1) then
        if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%bb(14-shf4))
      else
        write(ilog,*) 'Unsupported 5"-terminus type for residue &
   &',rs,' in PDB-file reader. Offending mode is ',&
   &moltermid(molofrs(rs),1),'.'
      end if
      if (moltermid(molofrs(rs),2).eq.1) then
        if (pdbnm.eq.'3HO*') call trfxyz(i,at(rs)%bb(15-shf4))
      else
        write(ilog,*) 'Unsupported 3"-terminus type for residue &
   &',rs,' in PDB-file reader. Offending mode is ',&
   &moltermid(molofrs(rs),2),'.'
      end if
    else if ((resname.eq.'R5P').OR.(resname.eq.'RPC').OR.&
   &    (resname.eq.'RPU').OR.(resname.eq.'RPT').OR.&
   &    (resname.eq.'RPA').OR.(resname.eq.'RPG')) then
      if (pdbnm.eq.' O3*') call trfxyz(i,nuci(rs,1))
      if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
      if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
      if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
      if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
      if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
      if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2O3*') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.' O2*') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HO*') call trfxyz(i,at(rs)%sc(7-shf2))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(10))
        if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(11))
        if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(12))
        if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(13))
        if (pdbnm.eq.' H2*') call trfxyz(i,at(rs)%sc(5))
        if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
      end if
      if (moltermid(molofrs(rs),1).eq.1) then
        if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%bb(14-shf4))
      else
        write(ilog,*) 'Unsupported 5"-terminus type for residue &
   &',rs,' in PDB-file reader. Offending mode is ',&
   &moltermid(molofrs(rs),1),'.'
      end if
      if (moltermid(molofrs(rs),2).eq.1) then
        if (pdbnm.eq.'3HO*') call trfxyz(i,at(rs)%bb(15-shf4))
      else
        write(ilog,*) 'Unsupported 3"-terminus type for residue &
   &',rs,' in PDB-file reader. Offending mode is ',&
   &moltermid(molofrs(rs),2),'.'
      end if
    else
      write(ilog,*) 'Fatal. The PDB-reader in mode 2 (incl. PDB-analys&
 &is mode) does not support free nucleotides of the selected residue type yet. Check back later.'
      call fexit()
    end if
  end if
  if (resname.eq.'NA+') then
    if (pdbnm.eq.' NA ') call trfxyz(i,at(rs)%bb(1))
  else if (resname.eq.'CL-') then
    if (pdbnm.eq.' CL ') call trfxyz(i,at(rs)%bb(1))
  else if (resname.eq.'K+ ') then
    if (pdbnm.eq.' K  ') call trfxyz(i,at(rs)%bb(1))
  else if (resname.eq.'BR-') then
    if (pdbnm.eq.' BR ') call trfxyz(i,at(rs)%bb(1))
  else if (resname.eq.'CS+') then
    if (pdbnm.eq.' CS ') call trfxyz(i,at(rs)%bb(1))
  else if (resname.eq.'I- ') then
    if (pdbnm.eq.' I  ') call trfxyz(i,at(rs)%bb(1))
  else if (resname.eq.'O2 ') then
    if (pdbnm.eq.' O1 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%bb(2))
  else if ((resname.eq.'SPC').OR.(resname.eq.'T3P')) then
    if (pdbnm.eq.' OW ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.'1HW ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'2HW ') call trfxyz(i,at(rs)%bb(3))
  else if ((resname.eq.'T4P').OR.(resname.eq.'T4E')) then
    if (pdbnm.eq.' OW ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.'1HW ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'2HW ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' MW ') call trfxyz(i,at(rs)%bb(4))
  else if (resname.eq.'T5P') then
    if (pdbnm.eq.' OW ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.'1HW ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'2HW ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'1LP ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'2LP ') call trfxyz(i,at(rs)%bb(5))
  else if (resname.eq.'NO3') then
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' O3 ') call trfxyz(i,at(rs)%bb(4))
  else if (resname.eq.'LCP') then
    if (pdbnm.eq.' CL ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' O3 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O4 ') call trfxyz(i,at(rs)%bb(5))
  else if (resname.eq.'AC-') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.'1OXT') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'2OXT') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' CH3') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%sc(4))
    end if
  else if (resname.eq.'GDN') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N2 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'1HN1') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.'2HN1') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.'1HN2') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.'2HN2') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.'1HN3') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.'2HN3') call trfxyz(i,at(rs)%bb(10))
  else if (resname.eq.'1MN') then
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'1HN ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'2HN ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'3HN ') call trfxyz(i,at(rs)%bb(5))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%sc(3))
    end if
  else if (resname.eq.'2MN') then
    if (pdbnm.eq.' CT1') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' CT2') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'1HN ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'2HN ') call trfxyz(i,at(rs)%bb(5))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT1') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT1') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT1') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'1HT2') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HT2') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'3HT2') call trfxyz(i,at(rs)%sc(6))
    end if
  else if (resname.eq.'URE') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'1HN1') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.'2HN1') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' N2 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'1HN2') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.'2HN2') call trfxyz(i,at(rs)%bb(8))
  else if (resname.eq.'NMF') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(3))
    if ((pdbnm.eq.' HN ').OR.(pdbnm.eq.' H  ')) then
      call trfxyz(i,at(rs)%bb(5-shf))
    end if
    if (pdbnm.eq.' CNT') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HC ') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.'1HNT') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'2HNT') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'3HNT') call trfxyz(i,at(rs)%sc(4))
    end if
  else if (resname.eq.'NMA') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(3))
    if ((pdbnm.eq.' HN ').OR.(pdbnm.eq.' H  ')) then
      call trfxyz(i,at(rs)%bb(4))
    end if
    if (pdbnm.eq.' CNT') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' CCT') call trfxyz(i,at(rs)%sc(2))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HNT') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'2HNT') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'3HNT') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'1HCT') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HCT') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'3HCT') call trfxyz(i,at(rs)%sc(8))
    end if
  else if (resname.eq.'ACA') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'1HN ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'2HN ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' CCT') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HCT') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'2HCT') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'3HCT') call trfxyz(i,at(rs)%sc(4))
    end if
  else if (resname.eq.'PPA') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'1HN ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'2HN ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' CCT') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' CBT') call trfxyz(i,at(rs)%sc(2))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HCT') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'2HCT') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'1HBT') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'2HBT') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'3HBT') call trfxyz(i,at(rs)%sc(7))
    end if
  else if (resname.eq.'CH4') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%bb(2))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%bb(3))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.'4HT ') call trfxyz(i,at(rs)%bb(5))
    end if
  else if (resname.eq.'NH4') then
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.'1HN ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'2HN ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'3HN ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.'4HN ') call trfxyz(i,at(rs)%bb(5))
  else if (resname.eq.'FOA') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HC ') call trfxyz(i,at(rs)%bb(4))
    end if
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.'1HN ') call trfxyz(i,at(rs)%bb(5-shf))
    if (pdbnm.eq.'2HN ') call trfxyz(i,at(rs)%bb(6-shf))
  else if (resname.eq.'MOH') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%bb(2))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%bb(3))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%bb(4))
    end if
    if (pdbnm.eq.' HO ') call trfxyz(i,at(rs)%sc(2))
  else if (resname.eq.'DMA') then
    if (pdbnm.eq.' C  ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N  ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C1 ') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' CCT') call trfxyz(i,at(rs)%sc(3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H1 ') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2H1 ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'3H1 ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1H2 ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2H2 ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'3H2 ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'1HCT') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'2HCT') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'3HCT') call trfxyz(i,at(rs)%sc(12))
    end if
  else if (resname.eq.'PCR') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C21') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C22') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' C31') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C32') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(8))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HCT') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HCT') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HCT') call trfxyz(i,at(rs)%sc(3))
    end if
    if (pdbnm.eq.' HO ') call trfxyz(i,at(rs)%sc(4-shf3))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H21') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.' H22') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H31') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H32') call trfxyz(i,at(rs)%bb(12))
    end if
  else if (resname.eq.'PRP') then
    if (pdbnm.eq.' CT1') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' CT2') call trfxyz(i,at(rs)%bb(3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT1') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT1') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT1') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'1HT2') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HT2') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'3HT2') call trfxyz(i,at(rs)%sc(6))
    end if
  else if (resname.eq.'NBU') then
    if (pdbnm.eq.' CT1') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' CB1') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' CB2') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' CT2') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT1') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'2HT1') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'3HT1') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'1HB1') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.'2HB1') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'1HB2') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'2HB2') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HT2') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HT2') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'3HT2') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'IBU') then
    if (pdbnm.eq.' CT1') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' CT2') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' CT3') call trfxyz(i,at(rs)%bb(4))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT1') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT1') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT1') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.' HB ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'1HT2') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HT2') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'3HT2') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HT3') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HT3') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'3HT3') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'BEN') then
    if (pdbnm.eq.' C1 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C3 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(6))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%bb(7))
      if (pdbnm.eq.' H2 ') call trfxyz(i,at(rs)%bb(8))
      if (pdbnm.eq.' H3 ') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.' H4 ') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H5 ') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%bb(12))
    end if
  else if (resname.eq.'TOL') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C21') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C22') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' C31') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C32') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(7))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%sc(3))
    end if
    if (ua_model.lt.2) then
      if (pdbnm.eq.'1H2 ') call trfxyz(i,at(rs)%bb(8))
      if (pdbnm.eq.'2H2 ') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.'1H3 ') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.'2H3 ') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H4 ') call trfxyz(i,at(rs)%bb(12))
    end if
  else if (resname.eq.'MSH') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' S  ') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%bb(2))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%bb(3))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%bb(4))
    end if
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HS ') call trfxyz(i,at(rs)%sc(2))
    end if
  else if (resname.eq.'EOH') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' O  ') call trfxyz(i,at(rs)%bb(3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%bb(5))
    end if
    if (pdbnm.eq.' HO ') call trfxyz(i,at(rs)%sc(4-shf3))
  else if (resname.eq.'EMT') then
    if (pdbnm.eq.' CNT') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' S  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' CCT') call trfxyz(i,at(rs)%sc(1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HNT') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'2HNT') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'3HNT') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HCT') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HCT') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'3HCT') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'IMD') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N2 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' C3 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(6))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%sc(3))
    end if
    if (pdbnm.eq.' HN2') call trfxyz(i,at(rs)%bb(7))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HC2') call trfxyz(i,at(rs)%bb(8))
      if (pdbnm.eq.' HC3') call trfxyz(i,at(rs)%bb(9))
    end if
  else if (resname.eq.'IME') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C1 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' N2 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' C3 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(6))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%sc(3))
    end if
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HC2') call trfxyz(i,at(rs)%bb(7))
      if (pdbnm.eq.' HC3') call trfxyz(i,at(rs)%bb(8))
    end if
    if (pdbnm.eq.' HN3') call trfxyz(i,at(rs)%bb(9-shfx2))
  else if (resname.eq.'NAP') then
    if (pdbnm.eq.' C11') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C12') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C21') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C22') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' C31') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C32') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C23') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' C24') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' C33') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.' C34') call trfxyz(i,at(rs)%bb(10))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H21') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H22') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H31') call trfxyz(i,at(rs)%bb(13))
      if (pdbnm.eq.' H32') call trfxyz(i,at(rs)%bb(14))
      if (pdbnm.eq.' H23') call trfxyz(i,at(rs)%bb(15))
      if (pdbnm.eq.' H24') call trfxyz(i,at(rs)%bb(16))
      if (pdbnm.eq.' H33') call trfxyz(i,at(rs)%bb(17))
      if (pdbnm.eq.' H34') call trfxyz(i,at(rs)%bb(18))
    end if
  else if (resname.eq.'MIN') then
    if (pdbnm.eq.' CT ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C3 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C3A') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C7A') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' C7 ') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(10))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HT ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HT ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3HT ') call trfxyz(i,at(rs)%sc(3))
    end if
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H2 ') call trfxyz(i,at(rs)%bb(11))
    end if
    if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%bb(12-shfx1))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H4 ') call trfxyz(i,at(rs)%bb(13))
      if (pdbnm.eq.' H7 ') call trfxyz(i,at(rs)%bb(14))
      if (pdbnm.eq.' H5 ') call trfxyz(i,at(rs)%bb(15))
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%bb(16))
    end if
  else if (resname.eq.'CYT') then
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' N4 ') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%bb(9))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H5 ') call trfxyz(i,at(rs)%bb(11))
    end if
    if (pdbnm.eq.'1H4 ') call trfxyz(i,at(rs)%bb(12-shfx2))
    if (pdbnm.eq.'2H4 ') call trfxyz(i,at(rs)%bb(13-shfx2))
  else if (resname.eq.'URA') then
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' O4 ') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%bb(9))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%bb(10))
    end if
    if (pdbnm.eq.' H3 ') call trfxyz(i,at(rs)%bb(11-shfx1))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H5 ') call trfxyz(i,at(rs)%bb(12))
    end if
  else if (resname.eq.'THY') then
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' C5M') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' O4 ') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%bb(10))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%bb(11))
    end if
    if (pdbnm.eq.' H3 ') call trfxyz(i,at(rs)%bb(12-shfx1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5M') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2H5M') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'3H5M') call trfxyz(i,at(rs)%sc(3))
    end if
  else if (resname.eq.'PUR') then
    if (pdbnm.eq.' N9 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C8 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' N7 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.' H9 ') call trfxyz(i,at(rs)%bb(10))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H8 ') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H2 ') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%bb(13))
    end if
  else if (resname.eq.'ADE') then
    if (pdbnm.eq.' N9 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C8 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' N7 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.' N6 ') call trfxyz(i,at(rs)%bb(10))
    if (pdbnm.eq.' H9 ') call trfxyz(i,at(rs)%bb(11))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H8 ') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H2 ') call trfxyz(i,at(rs)%bb(13))
    end if
    if (pdbnm.eq.'1H6 ') call trfxyz(i,at(rs)%bb(14-shfx2))
    if (pdbnm.eq.'2H6 ') call trfxyz(i,at(rs)%bb(15-shfx2))
  else if (resname.eq.'GUA') then
    if (pdbnm.eq.' N9 ') call trfxyz(i,at(rs)%bb(1))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.' C8 ') call trfxyz(i,at(rs)%bb(3))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' N7 ') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%bb(6))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%bb(7))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%bb(8))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.' N2 ') call trfxyz(i,at(rs)%bb(10))
    if (pdbnm.eq.' O6 ') call trfxyz(i,at(rs)%bb(11))
    if (pdbnm.eq.' H9 ') call trfxyz(i,at(rs)%bb(12))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H8 ') call trfxyz(i,at(rs)%bb(13))
    end if
    if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%bb(14-shfx1))
    if (pdbnm.eq.'1H2 ') call trfxyz(i,at(rs)%bb(15-shfx1))
    if (pdbnm.eq.'2H2 ') call trfxyz(i,at(rs)%bb(16-shfx1))
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine trf_nts(i,rs,resname,pdbnm)
!
  use pdb
  use iounit
  use polypep
  use sequen
  use aminos
  use molecule
  use system
  use atoms
!
  implicit none
!
  integer i,rs,shf2,shf4
  character(3) resname
  character(4) pdbnm
!
  shf2 = 0
  shf4 = 0
  if (ua_model.gt.0) then
    shf2 = 2
    shf4 = 4
  end if
!
  if (resname.eq.'ACE') then
    if (pdbnm.eq.' CH3') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%sc(2))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%sc(4))
    end if
  else if (resname.eq.'FOR') then
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' H  ') call trfxyz(i,at(rs)%bb(3))
    end if
  else if (resname.eq.'GLY') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (moltermid(molofrs(rs),1).eq.1) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
      if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%bb(7))
    else if (moltermid(molofrs(rs),1).eq.2) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
    else
      write(ilog,*) 'Unsupported N-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HA ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HA ') call trfxyz(i,at(rs)%sc(2))
    end if
  else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.32)) then ! PRO,HYP
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (moltermid(molofrs(rs),1).eq.1) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
    else if (moltermid(molofrs(rs),1).eq.2) then
      if (pdbnm.eq.' H  ') call trfxyz(i,at(rs)%bb(5))
    else
      write(ilog,*) 'Unsupported N-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
    end if
  else if (resname.eq.'AIB') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (moltermid(molofrs(rs),1).eq.1) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
      if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%bb(7))
    else if (moltermid(molofrs(rs),1).eq.2) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
     else
      write(ilog,*) 'Unsupported N-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
    end if

  else if ((resname.eq.'D5P').OR.(resname.eq.'DPC').OR.&
 &       (resname.eq.'DPU').OR.(resname.eq.'DPT').OR.&
 &       (resname.eq.'DPA').OR.(resname.eq.'DPG')) then
    if (nuci(rs+1,1).gt.0) then
      if (pdbnm(2:4).eq.'O3*') then
!       possible shift to rs+1 is handled previously in adjust_PONuc(...) and reorder_PDB_byrsnmb(...)
        call trfxyz(i,nuci(rs,1))
      end if
    else ! we must assume self-residue
      if (pdbnm(2:4).eq.'O3*') then
        call trfxyz(i,nuci(rs,1))
      end if
    end if
    if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.'1H2*') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
    if (moltermid(molofrs(rs),1).eq.1) then
      if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%bb(13-shf4))
    else
      write(ilog,*) 'Unsupported 5"-terminus type for residue &
 &',rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
    end if
  else if ((resname.eq.'DIB').OR.(resname.eq.'DIC').OR.&
 &       (resname.eq.'DIU').OR.(resname.eq.'DIT').OR.&
 &       (resname.eq.'DIA').OR.(resname.eq.'DIG')) then
    if ((nuci(rs+1,1).gt.0)) then ! .AND.((i-stat).gt.1)) then
!     note that in normal pdb-files the O3* is part of the residue
!     which makes up the sugar, not the one making up the following Pi
      if (pdbnm(2:4).eq.'O3*') call trfxyz(i,nuci(rs+1,1))
    else
      if (pdbnm(2:4).eq.'O3*') then
        write(ilog,*) "Fatal. O3* partitioning in nucleotide pdb-files must conform &
 &to CAMPARI standards if 5'-nucleosides are used in conjunction with unknown residues."
        call fexit()
      end if
    end if
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,1))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (pdbnm.eq.'5HO*') call trfxyz(i,at(rs)%bb(9-shf4))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(6))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(7))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(8))
      if (pdbnm.eq.'1H2*') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
  else if ((resname.eq.'R5P').OR.(resname.eq.'RPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'RPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'RPG')) then
    if (nuci(rs+1,1).gt.0) then
      if (pdbnm(2:4).eq.'O3*') then
!       possible shift to rs+1 is handled previously in adjust_PONuc(...) and reorder_PDB_byrsnmb(...)
        call trfxyz(i,nuci(rs,1))
      end if
    else ! we must assume self-residue
      if (pdbnm(2:4).eq.'O3*') then
        call trfxyz(i,nuci(rs,1))
      end if
    end if
    if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (pdbnm.eq.' O2*') call trfxyz(i,at(rs)%sc(4))
    if (pdbnm.eq.'2HO*') call trfxyz(i,at(rs)%sc(7-shf2))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
    if (moltermid(molofrs(rs),1).eq.1) then
      if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%bb(13-shf4))
    else
      write(ilog,*) 'Unsupported 5"-terminus type for residue &
 &',rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
    end if
  else if ((resname.eq.'RIB').OR.(resname.eq.'RIC').OR.&
 &    (resname.eq.'RIU').OR.(resname.eq.'RIT').OR.&
 &    (resname.eq.'RIA').OR.(resname.eq.'RIG')) then
    if ((nuci(rs+1,1).gt.0)) then !!.AND.((i-stat).gt.1)) then
!     see above
      if (pdbnm(2:4).eq.'O3*') call trfxyz(i,nuci(rs+1,1))
    else
      if (pdbnm(2:4).eq.'O3*') then
        write(ilog,*) "Fatal. O3* partitioning in nucleotide pdb-files must conform &
 &to CAMPARI standards if 5'-nucleosides are used in conjunction with unknown residues."
        call fexit()
      end if
    end if
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,1))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (pdbnm.eq.' O2*') call trfxyz(i,at(rs)%sc(4))
    if (pdbnm.eq.'2HO*') call trfxyz(i,at(rs)%sc(7-shf2))
    if (pdbnm.eq.'5HO*') call trfxyz(i,at(rs)%bb(9-shf4))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(6))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(7))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(8))
      if (pdbnm.eq.' H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
  else
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (moltermid(molofrs(rs),1).eq.1) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
      if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%bb(7))
    else if (moltermid(molofrs(rs),1).eq.2) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(6))
    else
      write(ilog,*) 'Unsupported N-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),1),'.'
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine trf_cts(i,rs,resname,pdbnm)
!
  use pdb
  use iounit
  use polypep
  use sequen
  use aminos
  use molecule
  use system
  use atoms
!
  implicit none
!
  integer i,rs,shf2,shf4
  character(3) resname
  character(4) pdbnm
!
  shf2 = 0
  shf4 = 0
  if (ua_model.gt.0) then
    shf2 = 2
    shf4 = 4
  end if
!
  if (resname.eq.'NME') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CH3') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(4))
      if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(5))
      if (pdbnm.eq.'3H  ') call trfxyz(i,at(rs)%bb(6))
    end if
  else if (resname.eq.'NH2') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.'1H  ') call trfxyz(i,at(rs)%bb(2))
    if (pdbnm.eq.'2H  ') call trfxyz(i,at(rs)%bb(3))
  else if ((resname.eq.'D5P').OR.(resname.eq.'DPC').OR.&
 &       (resname.eq.'DPU').OR.(resname.eq.'DPT').OR.&
 &       (resname.eq.'DPA').OR.(resname.eq.'DPG')) then
    if (pdbnm(1:4).eq.' O3*') then
      call trfxyz(i,nuci(rs,1))
    end if
    if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (pdbnm.eq.'2O3*') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.'3HO*') call trfxyz(i,at(rs)%bb(14-shf4))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(13))
      if (pdbnm.eq.'1H2*') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
  else if ((resname.eq.'R5P').OR.(resname.eq.'RPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'RPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'RPG')) then
    if (pdbnm(1:4).eq.' O3*') then
      call trfxyz(i,nuci(rs,1))
    end if
    if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (pdbnm.eq.' O2*') call trfxyz(i,at(rs)%sc(4))
    if (pdbnm.eq.'2O3*') call trfxyz(i,at(rs)%bb(9))
    if (pdbnm.eq.'2HO*') call trfxyz(i,at(rs)%sc(7-shf2))
    if (pdbnm.eq.'3HO*') call trfxyz(i,at(rs)%bb(14-shf4))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(13))
      if (pdbnm.eq.' H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
  else
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.'1OXT') call trfxyz(i,oi(rs))
    if (pdbnm.eq.'2OXT') call trfxyz(i,at(rs)%bb(5))
    if (moltermid(molofrs(rs),2).ne.1) then
      write(ilog,*) 'Unsupported C-terminus type for residue '&
 &,rs,' in PDB-file reader. Offending mode is ',&
 &moltermid(molofrs(rs),2),'.'
       call fexit()
    end if
    if (resname.eq.'GLY') then
      if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1HA ') call trfxyz(i,at(rs)%sc(1))
        if (pdbnm.eq.'2HA ') call trfxyz(i,at(rs)%sc(2))
      end if
    else if (resname.eq.'AIB') then
      if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs)) 
    else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.32)) then ! PRO, HYP
      if (ua_model.eq.0) then
        if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
      end if
    else
      if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs))
      if (ua_model.eq.0) then
        if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
      end if
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine trf_bbs(i,rs,resname,pdbnm)
!
  use pdb
  use iounit
  use polypep
  use sequen
  use aminos
  use molecule
  use system
  use atoms
!
  implicit none
!
  integer i,rs,shf2
  character(3) resname
  character(4) pdbnm
!
  shf2 = 0
  if (ua_model.gt.0) then
    shf2 =2
  end if
!
  if ((resname.eq.'D5P').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'DPU').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'DPA').OR.(resname.eq.'DPG')) then
    if (nuci(rs+1,1).gt.0) then
      if (pdbnm(2:4).eq.'O3*') then
!       possible shift to rs+1 is handled previously in adjust_PONuc(...) and reorder_PDB_byrsnmb(...)
        call trfxyz(i,nuci(rs,1))
      end if
    else ! we must assume self-residue
      if (pdbnm(2:4).eq.'O3*') then
        call trfxyz(i,nuci(rs,1))
      end if
    end if
    if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.'1H2*') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
  else if ((resname.eq.'R5P').OR.(resname.eq.'RPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'RPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'RPG')) then
    if (nuci(rs+1,1).gt.0) then
      if (pdbnm(2:4).eq.'O3*') then
!       possible shift to rs+1 is handled previously in adjust_PONuc(...) and reorder_PDB_byrsnmb(...)
        call trfxyz(i,nuci(rs,1))
      end if
    else ! we must assume self-residue
      if (pdbnm(2:4).eq.'O3*') then
        call trfxyz(i,nuci(rs,1))
      end if
    end if
    if (pdbnm.eq.' P  ') call trfxyz(i,nuci(rs,2))
    if (pdbnm.eq.' O5*') call trfxyz(i,nuci(rs,3))
    if (pdbnm.eq.' C5*') call trfxyz(i,nuci(rs,4))
    if (pdbnm.eq.' C4*') call trfxyz(i,nuci(rs,5))
    if (pdbnm.eq.' C3*') call trfxyz(i,nuci(rs,6))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%bb(4))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%bb(5))
    if (pdbnm.eq.' C2*') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' C1*') call trfxyz(i,at(rs)%sc(2))
    if (pdbnm.eq.' O4*') call trfxyz(i,at(rs)%sc(3))
    if (pdbnm.eq.' O2*') call trfxyz(i,at(rs)%sc(4))
    if (pdbnm.eq.'2HO*') call trfxyz(i,at(rs)%sc(7-shf2))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5*') call trfxyz(i,at(rs)%bb(9))
      if (pdbnm.eq.'2H5*') call trfxyz(i,at(rs)%bb(10))
      if (pdbnm.eq.' H4*') call trfxyz(i,at(rs)%bb(11))
      if (pdbnm.eq.' H3*') call trfxyz(i,at(rs)%bb(12))
      if (pdbnm.eq.' H2*') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.' H1*') call trfxyz(i,at(rs)%sc(6))
    end if
  else if (resname.eq.'GLY') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HA ') call trfxyz(i,at(rs)%sc(1))
      if (pdbnm.eq.'2HA ') call trfxyz(i,at(rs)%sc(2))
    end if
  else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.32)) then ! PRO, HYP
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
    end if
  else if (resname.eq.'AIB') then
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs))
  else
    if (pdbnm.eq.' N  ') call trfxyz(i,ni(rs))
    if (pdbnm.eq.' CA ') call trfxyz(i,cai(rs))
    if (pdbnm.eq.' C  ') call trfxyz(i,ci(rs))
    if (pdbnm.eq.' O  ') call trfxyz(i,oi(rs))
    if (pdbnm.eq.' H  ') call trfxyz(i,hni(rs))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HA ') call trfxyz(i,at(rs)%sc(1))
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine trf_scs(i,rs,resname,pdbnm)
!
  use pdb
  use iounit
  use polypep
  use sequen
  use aminos
  use molecule
  use system
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer i,rs,kk,shf,shf2,shf3,shf5,shf7,shf9,shfx1,shfx2,shfx4,vshf
  character(3) resname
  character(4) pdbnm
!
  shf = 0
  shf2 = 0
  shf3 = 0
  shf5 = 0
  shf7 = 0
  shf9 = 0
  shfx1 = 0
  shfx2 = 0
  shfx4 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
    shf5 = 5
    shf7 = 7
    shf9 = 9
    if (ua_model.eq.2) then
      shfx1 = 1
      shfx2 = 2
      shfx4 = 4
    end if
  end if
!
  if (resname.eq.'ALA') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'3HB ') call trfxyz(i,at(rs)%sc(5))
    end if
  else if (resname.eq.'AIB') then
    if (pdbnm.eq.' CB1') call trfxyz(i,at(rs)%sc(1))
    if (pdbnm.eq.' CB2') call trfxyz(i,at(rs)%sc(2))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB1') call trfxyz(i,at(rs)%sc(3))
      if (pdbnm.eq.'2HB1') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'3HB1') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'1HB2') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB2') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'3HB2') call trfxyz(i,at(rs)%sc(8))
    end if
  else if (resname.eq.'ABA') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'3HG ') call trfxyz(i,at(rs)%sc(8))
    end if
  else if (resname.eq.'VAL') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG1') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CG2') call trfxyz(i,at(rs)%sc(4-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'1HG1') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HG1') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'3HG1') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HG2') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HG2') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'3HG2') call trfxyz(i,at(rs)%sc(11))
    end if
  else if (resname.eq.'NVA') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'3HD ') call trfxyz(i,at(rs)%sc(11))
    end if
  else if (resname.eq.'LEU') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.' HG ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HD1') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HD1') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'3HD1') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'1HD2') call trfxyz(i,at(rs)%sc(12))
      if (pdbnm.eq.'2HD2') call trfxyz(i,at(rs)%sc(13))
      if (pdbnm.eq.'3HD2') call trfxyz(i,at(rs)%sc(14))
    end if
  else if (resname.eq.'ILE') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG1') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CG2') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD1') call trfxyz(i,at(rs)%sc(5-shf))
    if (ua_model.eq.0) then
    if (pdbnm.eq.' HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HG1') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HG1') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HG2') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HG2') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'3HG2') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'1HD1') call trfxyz(i,at(rs)%sc(12))
      if (pdbnm.eq.'2HD1') call trfxyz(i,at(rs)%sc(13))
      if (pdbnm.eq.'3HD1') call trfxyz(i,at(rs)%sc(14))
    end if
  else if (resname.eq.'NLE') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CE ') call trfxyz(i,at(rs)%sc(5-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'1HE ') call trfxyz(i,at(rs)%sc(12))
      if (pdbnm.eq.'2HE ') call trfxyz(i,at(rs)%sc(13))
      if (pdbnm.eq.'3HE ') call trfxyz(i,at(rs)%sc(14))
    end if
  else if (resname.eq.'PRO') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(10))
    end if
  else if (resname.eq.'MET') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' SD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CE ') call trfxyz(i,at(rs)%sc(5-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'1HE ') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'2HE ') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'3HE ') call trfxyz(i,at(rs)%sc(12))
    end if
  else if (resname.eq.'PHE') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' CE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' CE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' CZ ') call trfxyz(i,at(rs)%sc(8-shf))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD1') call trfxyz(i,at(rs)%sc(11-shf3))
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(12-shf3))
      if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(13-shf3))
      if (pdbnm.eq.' HE2') call trfxyz(i,at(rs)%sc(14-shf3))
      if (pdbnm.eq.' HZ ') call trfxyz(i,at(rs)%sc(15-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(10))
    end if
  else if ((resname.eq.'TYR').OR.(resname.eq.'TYO')) then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' CE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' CE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' CZ ') call trfxyz(i,at(rs)%sc(8-shf))
    if (resname.eq.'TYR') then
      if (pdbnm.eq.' OH ') call trfxyz(i,at(rs)%sc(9-shf))
      if (pdbnm.eq.' HH ') call trfxyz(i,at(rs)%sc(16-shf3-shfx4))
    else if (resname.eq.'TYO') then
      if (pdbnm.eq.' OZ ') call trfxyz(i,at(rs)%sc(9-shf))
    end if
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD1') call trfxyz(i,at(rs)%sc(12-shf3))
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(13-shf3))
      if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(14-shf3))
      if (pdbnm.eq.' HE2') call trfxyz(i,at(rs)%sc(15-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(11))
    end if
  else if (resname.eq.'PTR') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' CE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' CE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' CZ ') call trfxyz(i,at(rs)%sc(8-shf))
    if (pdbnm.eq.' OZ ') call trfxyz(i,at(rs)%sc(9-shf))
    if (pdbnm.eq.' P  ') call trfxyz(i,at(rs)%sc(10-shf))
    if (pdbnm.eq.' OH ') call trfxyz(i,at(rs)%sc(11-shf))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%sc(12-shf))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%sc(13-shf))
    if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%sc(20-shf3-shfx4))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD1') call trfxyz(i,at(rs)%sc(16-shf3))
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(17-shf3))
      if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(18-shf3))
      if (pdbnm.eq.' HE2') call trfxyz(i,at(rs)%sc(19-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(14))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(15))
    end if
  else if (resname.eq.'HIE') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' ND1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' CE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' NE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' HE2') call trfxyz(i,at(rs)%sc(12-shf3-shfx2))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(10-shf3))
      if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(11-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'HID') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' ND1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' HD1') call trfxyz(i,at(rs)%sc(10-shf3))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' CE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' NE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(11-shf3))
      if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(12-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'HIP') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' ND1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' HD1') call trfxyz(i,at(rs)%sc(10-shf3))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' CE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' NE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' HE2') call trfxyz(i,at(rs)%sc(13-shf3-shfx2))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(11-shf3))
      if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(12-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'TRP') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' NE1') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' HE1') call trfxyz(i,at(rs)%sc(15-shf3-shfx1))
    if (pdbnm.eq.' CE2') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' CE3') call trfxyz(i,at(rs)%sc(8-shf))
    if (pdbnm.eq.' CZ2') call trfxyz(i,at(rs)%sc(9-shf))
    if (pdbnm.eq.' CZ3') call trfxyz(i,at(rs)%sc(10-shf))
    if (pdbnm.eq.' CH2') call trfxyz(i,at(rs)%sc(11-shf))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' HD1') call trfxyz(i,at(rs)%sc(14-shf3))
      if (pdbnm.eq.' HE3') call trfxyz(i,at(rs)%sc(16-shf3))
      if (pdbnm.eq.' HZ2') call trfxyz(i,at(rs)%sc(17-shf3))
      if (pdbnm.eq.' HZ3') call trfxyz(i,at(rs)%sc(18-shf3))
      if (pdbnm.eq.' HH2') call trfxyz(i,at(rs)%sc(19-shf3))
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(12))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(13))
    end if
  else if (resname.eq.'SER') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' OG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' HG ') call trfxyz(i,at(rs)%sc(6-shf3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(5))
    end if
  else if (resname.eq.'SEP') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' OG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' P  ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' OH ') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%sc(10-shf3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(9))
    end if
  else if ((resname.eq.'CYS').OR.(resname.eq.'CYX')) then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' SG ') call trfxyz(i,at(rs)%sc(3-shf))
    if ((ua_model.lt.2).AND.(resname.eq.'CYS')) then
      if (pdbnm.eq.' HG ') then
        if (disulf(rs).gt.0) then
          if (crlk_idx(rs).gt.0) then
            if ((crosslink(crlk_idx(rs))%itstype.gt.2).OR.(crosslink(crlk_idx(rs))%itstype.le.0)) then
              if (at(rs)%nsc.ge.(6-shf3)) call trfxyz(i,at(rs)%sc(6-shf3))
            end if
          else
            if (at(rs)%nsc.ge.(6-shf3)) call trfxyz(i,at(rs)%sc(6-shf3))
          end if
!        else if (disulf(rs).gt.0) then
!          do kk=1,n_crosslinks
!            if (((crosslink(kk)%rsnrs(1).eq.rs).OR.(crosslink(kk)%rsnrs(2).eq.rs)).AND.&
! &              ((crosslink(kk)%itstype.gt.2).OR.(crosslink(kk)%itstype.le.0))) then
!                call trfxyz(i,at(rs)%sc(6-shf3))
!              exit
!            end if
!          end do
        else if (at(rs)%nsc.ge.(6-shf3)) then ! crude but needed if this function is used before crlk data str are set up
          call trfxyz(i,at(rs)%sc(6-shf3))
        end if
      end if
    end if
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(4))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(5))
    end if
  else if (resname.eq.'THR') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' OG1') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' HG1') call trfxyz(i,at(rs)%sc(6-shf2))
    if (pdbnm.eq.' CG2') call trfxyz(i,at(rs)%sc(4-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'1HG2') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HG2') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'3HG2') call trfxyz(i,at(rs)%sc(9))
    end if
  else if (resname.eq.'TPO') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' OG1') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CG2') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' P  ') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' OH ') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' O1P') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.' O2P') call trfxyz(i,at(rs)%sc(8-shf))
    if (pdbnm.eq.' HOP') call trfxyz(i,at(rs)%sc(13-shf3-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.' HB ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'1HG2') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'2HG2') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'3HG2') call trfxyz(i,at(rs)%sc(12))
    end if
  else if (resname.eq.'ASN') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' OD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' ND2') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.'1HD2') call trfxyz(i,at(rs)%sc(8-shf3))
    if (pdbnm.eq.'2HD2') call trfxyz(i,at(rs)%sc(9-shf3))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(7))
    end if
  else if (resname.eq.'GLN') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' OE1') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' NE2') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.'1HE2') call trfxyz(i,at(rs)%sc(11-shf5))
    if (pdbnm.eq.'2HE2') call trfxyz(i,at(rs)%sc(12-shf5))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(10))
    end if
  else if ((resname.eq.'ASP').OR.(resname.eq.'ASH')) then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' OD1') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' OD2') call trfxyz(i,at(rs)%sc(5-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(7))
    end if
    if (resname.eq.'ASH') then
      if (pdbnm.eq.' HD2') call trfxyz(i,at(rs)%sc(8-shf3))
    end if
  else if ((resname.eq.'GLU').OR.(resname.eq.'GLH')) then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' OE1') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' OE2') call trfxyz(i,at(rs)%sc(6-shf))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(10))
    end if
    if (resname.eq.'GLH') then
      if (pdbnm.eq.' HE2') call trfxyz(i,at(rs)%sc(11-shf5))
    end if
  else if ((resname.eq.'LYS').OR.(resname.eq.'LYD').OR.(resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &         (resname.eq.'KM2').OR.(resname.eq.'KM3')) then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' CE ') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' NZ ') call trfxyz(i,at(rs)%sc(6-shf))
    if (ua_model.eq.0) then
      vshf = 0
      if ((resname.eq.'KAC').OR.(resname.eq.'KM3')) vshf = 3
      if (resname.eq.'KM2') vshf = 2
      if (resname.eq.'KM1') vshf = 1
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(7+vshf))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(8+vshf))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(9+vshf))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(10+vshf))
      if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(11+vshf))
      if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(12+vshf))
      if (pdbnm.eq.'1HE ') call trfxyz(i,at(rs)%sc(13+vshf))
      if (pdbnm.eq.'2HE ') call trfxyz(i,at(rs)%sc(14+vshf))
    end if
    if ((resname.eq.'LYS').OR.(resname.eq.'LYD')) then
      if (pdbnm.eq.'1HZ ') call trfxyz(i,at(rs)%sc(15-shf9))
      if (pdbnm.eq.'2HZ ') call trfxyz(i,at(rs)%sc(16-shf9))
      if ((resname.eq.'LYS').AND.(pdbnm.eq.'3HZ ')) call trfxyz(i,at(rs)%sc(17-shf9))
    else if (resname.eq.'KAC') then
      if (pdbnm.eq.' CH ') call trfxyz(i,at(rs)%sc(7-shf))
      if (pdbnm.eq.' CI ') call trfxyz(i,at(rs)%sc(8-shf))
      if (pdbnm.eq.' OH ') call trfxyz(i,at(rs)%sc(9-shf))
      if (pdbnm.eq.' HZ ') call trfxyz(i,at(rs)%sc(18-shf9))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1HI ') call trfxyz(i,at(rs)%sc(19-shf9))
        if (pdbnm.eq.'2HI ') call trfxyz(i,at(rs)%sc(20-shf9))
        if (pdbnm.eq.'3HI ') call trfxyz(i,at(rs)%sc(21-shf9))
      end if
    else if (resname.eq.'KM1') then
      if (pdbnm.eq.' CH ') call trfxyz(i,at(rs)%sc(7-shf))
      if (pdbnm.eq.'1HZ ') call trfxyz(i,at(rs)%sc(16-shf9))
      if (pdbnm.eq.'2HZ ') call trfxyz(i,at(rs)%sc(17-shf9))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1HH ') call trfxyz(i,at(rs)%sc(18-shf9))
        if (pdbnm.eq.'2HH ') call trfxyz(i,at(rs)%sc(19-shf9))
        if (pdbnm.eq.'3HH ') call trfxyz(i,at(rs)%sc(20-shf9))
      end if
    else if (resname.eq.'KM2') then
      if (pdbnm.eq.' CH1') call trfxyz(i,at(rs)%sc(7-shf))
      if (pdbnm.eq.' CH2') call trfxyz(i,at(rs)%sc(8-shf))
      if (pdbnm.eq.' HZ ') call trfxyz(i,at(rs)%sc(17-shf9))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1HH1') call trfxyz(i,at(rs)%sc(18-shf9))
        if (pdbnm.eq.'2HH1') call trfxyz(i,at(rs)%sc(19-shf9))
        if (pdbnm.eq.'3HH1') call trfxyz(i,at(rs)%sc(20-shf9))
        if (pdbnm.eq.'1HH2') call trfxyz(i,at(rs)%sc(21-shf9))
        if (pdbnm.eq.'2HH2') call trfxyz(i,at(rs)%sc(22-shf9))
        if (pdbnm.eq.'3HH2') call trfxyz(i,at(rs)%sc(23-shf9))
      end if
    else if (resname.eq.'KM3') then
      if (pdbnm.eq.' CH1') call trfxyz(i,at(rs)%sc(7-shf))
      if (pdbnm.eq.' CH2') call trfxyz(i,at(rs)%sc(8-shf))
      if (pdbnm.eq.' CH3') call trfxyz(i,at(rs)%sc(9-shf))
      if (ua_model.eq.0) then
        if (pdbnm.eq.'1HH1') call trfxyz(i,at(rs)%sc(18-shf9))
        if (pdbnm.eq.'2HH1') call trfxyz(i,at(rs)%sc(19-shf9))
        if (pdbnm.eq.'3HH1') call trfxyz(i,at(rs)%sc(20-shf9))
        if (pdbnm.eq.'1HH2') call trfxyz(i,at(rs)%sc(21-shf9))
        if (pdbnm.eq.'2HH2') call trfxyz(i,at(rs)%sc(22-shf9))
        if (pdbnm.eq.'3HH2') call trfxyz(i,at(rs)%sc(23-shf9))
        if (pdbnm.eq.'1HH3') call trfxyz(i,at(rs)%sc(24-shf9))
        if (pdbnm.eq.'2HH3') call trfxyz(i,at(rs)%sc(25-shf9))
        if (pdbnm.eq.'3HH3') call trfxyz(i,at(rs)%sc(26-shf9))
      end if
    end if
  else if (resname.eq.'ORN') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' NE ') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.'1HE ') call trfxyz(i,at(rs)%sc(12-shf7))
    if (pdbnm.eq.'2HE ') call trfxyz(i,at(rs)%sc(13-shf7))
    if (pdbnm.eq.'3HE ') call trfxyz(i,at(rs)%sc(14-shf7))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(8))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(11))
    end if
  else if (resname.eq.'DAB') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' ND ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(9-shf5))
    if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(10-shf5))
    if (pdbnm.eq.'3HD ') call trfxyz(i,at(rs)%sc(11-shf5))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(5))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(6))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(7))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(8))
    end if
  else if (resname.eq.'ARG') then
    if (pdbnm.eq.' CB ') call trfxyz(i,at(rs)%sc(2-shf))
    if (pdbnm.eq.' CG ') call trfxyz(i,at(rs)%sc(3-shf))
    if (pdbnm.eq.' CD ') call trfxyz(i,at(rs)%sc(4-shf))
    if (pdbnm.eq.' NE ') call trfxyz(i,at(rs)%sc(5-shf))
    if (pdbnm.eq.' HE ') call trfxyz(i,at(rs)%sc(15-shf7))
    if (pdbnm.eq.' CZ ') call trfxyz(i,at(rs)%sc(6-shf))
    if (pdbnm.eq.' NH1') call trfxyz(i,at(rs)%sc(7-shf))
    if (pdbnm.eq.'1HH1') call trfxyz(i,at(rs)%sc(16-shf7))
    if (pdbnm.eq.'2HH1') call trfxyz(i,at(rs)%sc(17-shf7))
    if (pdbnm.eq.' NH2') call trfxyz(i,at(rs)%sc(8-shf))
    if (pdbnm.eq.'1HH2') call trfxyz(i,at(rs)%sc(18-shf7))
    if (pdbnm.eq.'2HH2') call trfxyz(i,at(rs)%sc(19-shf7))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1HB ') call trfxyz(i,at(rs)%sc(9))
      if (pdbnm.eq.'2HB ') call trfxyz(i,at(rs)%sc(10))
      if (pdbnm.eq.'1HG ') call trfxyz(i,at(rs)%sc(11))
      if (pdbnm.eq.'2HG ') call trfxyz(i,at(rs)%sc(12))
      if (pdbnm.eq.'1HD ') call trfxyz(i,at(rs)%sc(13))
      if (pdbnm.eq.'2HD ') call trfxyz(i,at(rs)%sc(14))
    end if
  else if ((resname.eq.'R5P').OR.(resname.eq.'D5P').OR.&
 &         (resname.eq.'RIB').OR.(resname.eq.'DIB')) then
    if ((resname.eq.'R5P').OR.(resname.eq.'RIB')) kk = 7 - shf2
    if ((resname.eq.'D5P').OR.(resname.eq.'DIB')) kk = 6 - shf3
    if (pdbnm.eq.' O1*') call trfxyz(i,at(rs)%sc(kk+1))
    if (pdbnm.eq.'1HO*') call trfxyz(i,at(rs)%sc(kk+2))
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &         (resname.eq.'RIC').OR.(resname.eq.'DIC')) then
    if ((resname.eq.'RPC').OR.(resname.eq.'RIC')) kk = 7 - shf2
    if ((resname.eq.'DPC').OR.(resname.eq.'DIC')) kk = 6 - shf3
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%sc(kk+1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%sc(kk+2))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%sc(kk+3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%sc(kk+4))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%sc(kk+5))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%sc(kk+6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%sc(kk+7))
    if (pdbnm.eq.' N4 ') call trfxyz(i,at(rs)%sc(kk+8))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%sc(kk+9))
      if (pdbnm.eq.' H5 ') call trfxyz(i,at(rs)%sc(kk+10))
    end if
    if (pdbnm.eq.'1H4 ') call trfxyz(i,at(rs)%sc(kk+11-shfx2))
    if (pdbnm.eq.'2H4 ') call trfxyz(i,at(rs)%sc(kk+12-shfx2))
  else if ((resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &         (resname.eq.'RIU').OR.(resname.eq.'DIU')) then
    if ((resname.eq.'RPU').OR.(resname.eq.'RIU')) kk = 7 - shf2
    if ((resname.eq.'DPU').OR.(resname.eq.'DIU')) kk = 6 - shf3
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%sc(kk+1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%sc(kk+2))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%sc(kk+3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%sc(kk+4))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%sc(kk+5))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%sc(kk+6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%sc(kk+7))
    if (pdbnm.eq.' O4 ') call trfxyz(i,at(rs)%sc(kk+8))
    if (pdbnm.eq.' H3 ') call trfxyz(i,at(rs)%sc(kk+10-shfx1))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%sc(kk+9))
      if (pdbnm.eq.' H5 ') call trfxyz(i,at(rs)%sc(kk+11))
    end if
  else if ((resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &         (resname.eq.'RIT').OR.(resname.eq.'DIT')) then
    if ((resname.eq.'RPT').OR.(resname.eq.'RIT')) kk = 7 - shf2
    if ((resname.eq.'DPT').OR.(resname.eq.'DIT')) kk = 6 - shf3
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%sc(kk+1))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%sc(kk+2))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%sc(kk+3))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%sc(kk+4))
    if (pdbnm.eq.' O2 ') call trfxyz(i,at(rs)%sc(kk+5))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%sc(kk+6))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%sc(kk+7))
    if (pdbnm.eq.' C5M') call trfxyz(i,at(rs)%sc(kk+8))
    if (pdbnm.eq.' O4 ') call trfxyz(i,at(rs)%sc(kk+9))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H6 ') call trfxyz(i,at(rs)%sc(kk+10))
    end if
    if (pdbnm.eq.' H3 ') call trfxyz(i,at(rs)%sc(kk+11-shfx1))
    if (ua_model.eq.0) then
      if (pdbnm.eq.'1H5M') call trfxyz(i,at(rs)%sc(kk+12))
      if (pdbnm.eq.'2H5M') call trfxyz(i,at(rs)%sc(kk+13))
      if (pdbnm.eq.'3H5M') call trfxyz(i,at(rs)%sc(kk+14))
    end if
  else if ((resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &         (resname.eq.'RIA').OR.(resname.eq.'DIA')) then
    if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) kk = 7 - shf2
    if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) kk = 6 - shf3
    if (pdbnm.eq.' N9 ') call trfxyz(i,at(rs)%sc(kk+1))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%sc(kk+2))
    if (pdbnm.eq.' C8 ') call trfxyz(i,at(rs)%sc(kk+3))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%sc(kk+4))
    if (pdbnm.eq.' N7 ') call trfxyz(i,at(rs)%sc(kk+5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%sc(kk+6))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%sc(kk+7))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%sc(kk+8))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%sc(kk+9))
    if (pdbnm.eq.' N6 ') call trfxyz(i,at(rs)%sc(kk+10))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H8 ') call trfxyz(i,at(rs)%sc(kk+11))
      if (pdbnm.eq.' H2 ') call trfxyz(i,at(rs)%sc(kk+12))
    end if
    if (pdbnm.eq.'1H6 ') call trfxyz(i,at(rs)%sc(kk+13-shfx2))
    if (pdbnm.eq.'2H6 ') call trfxyz(i,at(rs)%sc(kk+14-shfx2))
  else if ((resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'DIG')) then
    if ((resname.eq.'RPG').OR.(resname.eq.'RIG')) kk = 7 - shf2
    if ((resname.eq.'DPG').OR.(resname.eq.'DIG')) kk = 6 - shf3
    if (pdbnm.eq.' N9 ') call trfxyz(i,at(rs)%sc(kk+1))
    if (pdbnm.eq.' C4 ') call trfxyz(i,at(rs)%sc(kk+2))
    if (pdbnm.eq.' C8 ') call trfxyz(i,at(rs)%sc(kk+3))
    if (pdbnm.eq.' C5 ') call trfxyz(i,at(rs)%sc(kk+4))
    if (pdbnm.eq.' N7 ') call trfxyz(i,at(rs)%sc(kk+5))
    if (pdbnm.eq.' N3 ') call trfxyz(i,at(rs)%sc(kk+6))
    if (pdbnm.eq.' C2 ') call trfxyz(i,at(rs)%sc(kk+7))
    if (pdbnm.eq.' C6 ') call trfxyz(i,at(rs)%sc(kk+8))
    if (pdbnm.eq.' N1 ') call trfxyz(i,at(rs)%sc(kk+9))
    if (pdbnm.eq.' N2 ') call trfxyz(i,at(rs)%sc(kk+10))
    if (pdbnm.eq.' O6 ') call trfxyz(i,at(rs)%sc(kk+11))
    if (ua_model.lt.2) then
      if (pdbnm.eq.' H8 ') call trfxyz(i,at(rs)%sc(kk+12))
    end if
    if (pdbnm.eq.' H1 ') call trfxyz(i,at(rs)%sc(kk+13-shfx1))
    if (pdbnm.eq.'1H2 ') call trfxyz(i,at(rs)%sc(kk+14-shfx1))
    if (pdbnm.eq.'2H2 ') call trfxyz(i,at(rs)%sc(kk+15-shfx1))
  end if
!
end
!  
!----------------------------------------------------------------------------------
!     
subroutine trfxyz(i,j)
!
  use pdb
  use atoms
!
  implicit none
!
  integer i,j
!
  x(j) = xpdb(i)
  y(j) = ypdb(i)
  z(j) = zpdb(i)
  if (use_pdb_template.EQV..true.) pdbmap(j) = i
  if (pdbmap2(j).le.0) then
    pdbmap2(j) = i
  else if (pdbmap2(j).ne.i) then
!   this is weird -> possibly warn/exit
    pdbmap2(j) = i
  end if
!
end
!
!----------------------------------------------------------------------------------
!
! reorder list of input PDB lines such that nucleotides become unproblematic
!
subroutine adjust_PONuc(pdbsz,rsnmbin,pdbnm,rsnam,xpin,ypin,zpin)
!
  use pdb
  use sequen
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: pdbsz
  character(4), INTENT(INOUT):: pdbnm(pdbsz)
  character(3), INTENT(INOUT):: rsnam(pdbsz)
  RTYPE, INTENT(IN):: xpin(pdbsz),ypin(pdbsz),zpin(pdbsz)
  integer, INTENT(INOUT):: rsnmbin(pdbsz)
!
  integer i,loix,oix,pix
  RTYPE thresh,pv(3),lpv(3),ov(3),lov(3),d2
  logical havep,haveo3star,hadp,hado3star,switched
  integer, ALLOCATABLE:: rsnmb(:)
!
  if (pdbsz.le.0) return
!
  havep = .false.
  do i=1,nseq
    if ((seqpolty(i).eq.'N').OR.(seqtyp(i).lt.0)) then
      havep = .true.
      exit
    end if
  end do
  if (havep.EQV..false.) return
!
  allocate(rsnmb(pdbsz))
  rsnmb(:) = rsnmbin(1:pdbsz)
  thresh = 1.9*1.9
!
  hadp = .false.
  havep = .false.
  hado3star = .false.
  haveo3star = .false.
  oix = 0
  loix = 0
  pix = 0
  do i=1,pdbsz+1
    switched = .false.
    if (i.gt.pdbsz) then
      switched = .true.
    else if (i.gt.1) then
      if (rsnmb(i).ne.rsnmb(i-1)) then
        switched = .true.
      end if
    end if
    if (switched.EQV..true.) then
      if ((havep.EQV..true.).AND.(haveo3star.EQV..true.)) then
        d2 = sum((pv(:)-ov(:))*(pv(:)-ov(:)))
        if (d2.le.thresh) then
!         this means we are most likely in CAMPARI convention
          havep = .false.
          haveo3star = .false.
        end if
      end if
      if ((havep.EQV..true.).AND.(hado3star.EQV..true.)) then
        d2 = sum((pv(:)-lov(:))*(pv(:)-lov(:)))
        if (d2.le.thresh) then
          rsnmbin(loix) = rsnmbin(loix) + 1
          rsnam(loix) = rsnam(pix)
          if (pdbnm(loix).eq." O3'") pdbnm(loix) = " O3*"
          havep = .false.
        end if
      end if
      hadp = havep
      havep = .false.
      hado3star = haveo3star
      haveo3star = .false.
    end if
    if (i.le.pdbsz) then
      if (pdbnm(i).eq.' P  ') then
        if (havep.EQV..true.) then
          write(ilog,*) 'In PDB input, a residue should never contain 2 atoms with exactly the same name. This is fatal &
   &for certain reference atoms (here " P  "). Note that some names are converted automatically, which may give rise to &
   &this error.'
          call fexit()
        end if
        havep = .true.
        pix = i
        lpv(:) = pv(:)
        pv(1) = xpin(i)
        pv(2) = ypin(i)
        pv(3) = zpin(i)
      else if ((pdbnm(i).eq.' O3*').OR.(pdbnm(i).eq." O3'")) then
        if (haveo3star.EQV..true.) then
          write(ilog,*) 'In PDB input, a residue should never contain 2 atoms with exactly the same name. This is fatal &
   &for certain reference atoms (here " O3*"). Note that some names are converted automatically, which may give rise to &
   &this error.'
          call fexit()
        end if
        haveo3star = .true.
        if (oix.gt.0) then
          loix = oix
          lov(:) = ov(:)
        end if
        oix = i
        ov(1) = xpin(i)
        ov(2) = ypin(i)
        ov(3) = zpin(i)
      end if
    end if
  end do
!
end
!
!-----------------------------------------------------------------------------------
!
subroutine reorder_PDB_byrsnmb(nl,kw,rsnmb,rsnam,pdbnm,xpin,ypin,zpin,lix,chnam)
!
  use interfaces
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: nl
  integer, INTENT(INOUT):: rsnmb(nl),lix(nl)
  character(3), INTENT(INOUT):: rsnam(nl)
  character(4), INTENT(INOUT):: pdbnm(nl)
  character(6), INTENT(INOUT):: kw(nl)
  character, INTENT(INOUT):: chnam(nl)
  RTYPE, INTENT(INOUT):: xpin(nl),ypin(nl),zpin(nl)
!
  character(6), ALLOCATABLE:: hlp(:)
  integer, ALLOCATABLE:: iv1(:),iv2(:),iv3(:),iv4(:)
  RTYPE, ALLOCATABLE:: rv(:)
  character, ALLOCATABLE:: cv(:)
  integer i,ii,jj,k,l,pix
  logical upornot
  RTYPE thresh,d21,d22
!
  allocate(iv4(nl))
  allocate(iv1(nl))
  allocate(iv2(nl))
  allocate(iv3(nl))
  iv1(:) = rsnmb(1:nl)
  do i=1,nl
    iv3(i) = i
  end do
  upornot = .true.
  i = nl
  ii = 1
  jj = nl
  call merge_sort(i,upornot,iv1,iv2,ii,jj,iv3,iv4)
  deallocate(iv1)
  allocate(rv(nl))
  rv(1:nl) = xpin(1:nl)
  do i=1,nl
    xpin(i) = rv(iv4(i))
  end do
  rv(1:nl) = ypin(1:nl)
  do i=1,nl
    ypin(i) = rv(iv4(i))
  end do
  rv(1:nl) = zpin(1:nl)
  do i=1,nl
    zpin(i) = rv(iv4(i))
  end do
  deallocate(rv)
  allocate(hlp(nl))
  do i=1,nl
    hlp(i)(1:3) = rsnam(i)
  end do
  do i=1,nl
    rsnam(i) = hlp(iv4(i))(1:3)
  end do
  hlp(1:nl)(1:4) = pdbnm(1:nl)
  do i=1,nl
    pdbnm(i) = hlp(iv4(i))(1:4)
  end do
  hlp(1:nl)(1:6) = kw(1:nl)
  do i=1,nl
    kw(i) = hlp(iv4(i))(1:6)
  end do
  deallocate(hlp)
  allocate(cv(nl))
  cv(1:nl) = chnam(1:nl)
  do i=1,nl
    chnam(i) = cv(iv4(i))
  end do
  deallocate(cv)
  ii = 0
  do i=1,nl
    if (i.gt.1) then
      if (iv2(i)-iv2(i-1).gt.1) then
        ii = ii + 1 - (iv2(i)-iv2(i-1))
      end if
    end if
    rsnmb(i) = iv2(i) + ii
  end do
  iv3(1:nl) = lix(1:nl)
  do i=1,nl
    lix(i) = iv3(iv4(i))
  end do
  deallocate(iv2)
  deallocate(iv3)
  deallocate(iv4)
!
! check for duplicated names
  thresh = 1.9*1.9
  ii = 1
  jj = 0
  do i=1,nl+1
    upornot = .false.
    if (i.gt.nl) then
      upornot = .true.
    else if (i.gt.1) then
      if (rsnmb(i).ne.rsnmb(i-1)) then
        upornot = .true.
      end if
    end if
    if (upornot.EQV..true.) then
      jj = i-1
      pix = 0
      do k=ii,jj
        if (pdbnm(k).eq.' P  ') then
          pix = k
          exit
        end if
      end do
      do k=ii,jj
        do l=k+1,jj
          if (pix.gt.0) then
            if (((pdbnm(k).eq.' O3*').OR.(pdbnm(k).eq." O3'").OR.(pdbnm(k).eq."2O3*").OR.(pdbnm(k).eq."2O3'")).AND.&
 &              ((pdbnm(l).eq.' O3*').OR.(pdbnm(l).eq." O3'").OR.(pdbnm(l).eq."2O3*").OR.(pdbnm(l).eq."2O3'"))) then
              d21 = (xpin(k)-xpin(pix))**2+(ypin(k)-ypin(pix))**2+(zpin(k)-zpin(pix))**2
              d22 = (xpin(l)-xpin(pix))**2+(ypin(l)-ypin(pix))**2+(zpin(l)-zpin(pix))**2
              if ((d21.le.thresh).AND.(d22.gt.thresh)) then
                pdbnm(k) = ' O3*'
                pdbnm(l) = '2O3*'
              else if ((d21.gt.thresh).AND.(d22.le.thresh)) then
                pdbnm(k) = '2O3*'
                pdbnm(l) = ' O3*'
              end if
            end if
          end if
          if (pdbnm(k).eq.pdbnm(l)) then
            write(ilog,*) 'Fatal. PDB input should never contain identical atom names within a single residue.'
            write(ilog,*) 'Encountered for name ',pdbnm(k),' (reordered records ',k,' and ',l,').'
            call fexit()
          end if
        end do
      end do
      ii = i
    end if
  end do
!
end
!
!----------------------------------------------------------------------------------
!
! here we're pursuing a strategy to use the geometry fixated in the pdb-file
! and to simply simulate from there, only adjusting the actual degrees of freedom,
! i.e., rigid body coordinates and dihedrals
!
subroutine FMSMC_readpdb3()
!
  use pdb
!
  implicit none
!
  integer k,aone
!
  k = pdb_fileformat !not that it matters
  aone = 1
  pdb_fileformat = 1
!
  call setup_pdbtraj()
  call FMSMC_readpdb2(aone)
!
  pdb_fileformat = k
!
end
!
!---------------------------------------------------------------------------------------
!
#ifdef LINK_XDR
!
! in reading xtc-files we make the (very strong) assumption that the order and number
! of atoms in the xtc-file matches that obtained by the builder EXACTLY.
! this functionality is useful only for CAMPARI-generated files unless a PDB template
! is used
! note we (still) assume NVT everywhere
!
subroutine FMSMC_readxtc(mdl)
! 
  use iounit
  use pdb
  use atoms
  use mcsums
  use molecule
  use zmatrix
  use sequen
  use clusters, ONLY: align
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  logical dofatal
  integer i,j,ret,natoms,magic,kk,kkk,t1,t2,mdl,imol
  real(KIND=4), ALLOCATABLE:: coords(:)
  real(KIND=4) boxer(3,3),xtcprec
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  integer tslash
  character(re_aux(10)) xpont
#endif
!
  if (mdl.eq.1) then
    dofatal = .true.
  else
    dofatal = .false.
  end if
!
  call strlims(xtcinfile,t1,t2)
#ifdef ENABLE_MPI
  do i=t2,t1,-1
    if (xtcinfile(i:i).eq.SLASHCHAR) exit
  end do
  tslash = i
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  xtcinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//xtcinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//xtcinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = xtcinfile(t1:t2)
#endif
!
! get the header-line
  call xdrfint(ipdbtraj,magic,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  if ((dofatal.EQV..false.).AND.(ret.ne.1)) then
    mdl = -mdl
    return
  end if
  call xdrfint(ipdbtraj,natoms,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  if ((dofatal.EQV..false.).AND.(ret.ne.1)) then
    mdl = -mdl
    return
  end if
  if (natoms.ne.n) then
    write(ilog,*) 'Fatal. XTC-file (',fn(t1:t2),') does not m&
 &atch the number of atoms defined by the sequence-file.'
    call fexit()
  end if
  call xdrfint(ipdbtraj,kk,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  if ((dofatal.EQV..false.).AND.(ret.ne.1)) then
    mdl = -mdl
    return
  end if
  if (kk.ne.mdl) then
!    write(ilog,*) 'Warning. XTC-file might have inconsistently numbe
! &red snapshots. The file might be corrupt.'
  end if
  call xdrfint(ipdbtraj,kkk,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal) 
  if ((dofatal.EQV..false.).AND.(ret.ne.1)) then
    mdl = -mdl
    return
  end if
!  write(*,*) kk,kkk,magic,natoms
!
! get box - note that since there's no NPT yet, this information is just flushed (WARNING)
  do i=1,3
    do j=1,3
      call xdrffloat(ipdbtraj,boxer(i,j),ret)
    end do
  end do
!
! get and transfer the coords
  allocate(coords(3*n))
  xtcprec = xtc_prec
  call xdrf3dfcoord(ipdbtraj,coords,natoms,xtcprec,ret)
  if (ret.ne.1) call fail_readxtc(ret,dofatal)
  if ((dofatal.EQV..false.).AND.(ret.ne.1)) then
    mdl = -mdl
    return
  end if
  if (use_pdb_template.EQV..true.) then
   do i=1,n
      x(i) = 10.0*coords(3*pdbmap(i)-2)
      y(i) = 10.0*coords(3*pdbmap(i)-1)
      z(i) = 10.0*coords(3*pdbmap(i))
    end do
  else
    do i=1,n
      x(i) = 10.0*coords(3*i-2)
      y(i) = 10.0*coords(3*i-1)
      z(i) = 10.0*coords(3*i)
    end do
  end if
  deallocate(coords)
!
! now generate appropriate internal coordinates from the read-in Cartesian values
  do imol=1,nmol
    call genzmat(imol)
  end do
!     
! and transfer them into FMSMC-arrays (which are the ones used in analysis functions)
  call zmatfyc2()
  if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
    do imol=1,nmol
      call makexyz_formol(imol)
    end do
  end if
!
! and finally update rigid-body coordinates
  do imol=1,nmol
    call update_rigid(imol)
    call update_rigidm(imol)
  end do
  if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
    if (align%yes.EQV..true.) then
      do imol=1,nmol
        call update_image2(imol)
      end do
    else
      do imol=1,nmol
        call update_image(imol)
      end do
    end if
  end if
!
end
!
#endif
!
!-----------------------------------------------------------------------
!
! in reading dcd-files we make analogous assumptions
!
subroutine FMSMC_readdcd(mdl)
! 
  use iounit
  use pdb
  use atoms
  use mcsums
  use molecule
  use zmatrix
  use sequen
  use clusters, ONLY: align
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer j,natoms,t1,t2,mdl,imol,iis(100),iomess
  real(KIND=4), ALLOCATABLE:: coords(:)
  real(KIND=8) boxer(6),fakedt
  character(4) shortti
  character(80) timest,longti
  character(MAXSTRLEN) fn
  logical dofatal
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
  integer tslash
#endif
!
  dofatal = .true.
  call strlims(dcdinfile,t1,t2)
#ifdef ENABLE_MPI
  do j=t2,t1,-1
    if (dcdinfile(j:j).eq.SLASHCHAR) exit
  end do
  tslash = j
  call int2str(myrank,xpont,re_aux(10))
  if ((tslash.ge.t1).AND.(tslash.lt.t2)) then
    fn =  dcdinfile(t1:tslash)//'N_'//xpont(1:re_aux(10))//'_'//dcdinfile((tslash+1):t2)
  else
    fn =  'N_'//xpont(1:re_aux(10))//'_'//dcdinfile(t1:t2)
  end if
  call strlims(fn,t1,t2)
#else
  fn(t1:t2) = dcdinfile(t1:t2)
#endif
!
  if (mdl.eq.1) then
!   get the header
    read(ipdbtraj,iostat=iomess) shortti,iis(1),iis(2),iis(3),iis(4),iis(5),iis(6),iis(7),iis(8),iis(9),&
 &        fakedt,iis(10),iis(11),iis(12),iis(13),iis(14),iis(15),iis(16),iis(17),iis(18)
    if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
    read(ipdbtraj,iostat=iomess) iis(19),longti,timest
    if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
    read(ipdbtraj,iostat=iomess) natoms
    if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
    if (natoms.ne.n) then
      write(ilog,*) 'Fatal. DCD-file (',fn(t1:t2),') does not m&
 &atch the number of atoms defined by the sequence-file.'
      call fexit()
    end if
  end if
!
  if (mdl.gt.1) dofatal = .false.
!
! get box - note that since there's no NPT yet, this information is just flushed (WARNING)
  if (dcd_withbox.EQV..true.) then
    read(ipdbtraj,iostat=iomess) (boxer(j),j=1,6)
    if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
    if ((dofatal.EQV..false.).AND.(iomess.ne.0)) then
      mdl = -mdl
      return
    end if
  end if
!
! get and transfer the coords
  allocate(coords(n))
  read(ipdbtraj,iostat=iomess) (coords(j),j=1,n)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  if ((dofatal.EQV..false.).AND.(iomess.ne.0)) then
    mdl = -mdl
    return
  end if
  if (use_pdb_template.EQV..true.) then
    x(1:n) = 1.0*coords(pdbmap(1:n))
  else
    x(1:n) = 1.0*coords(1:n)
  end if
  read(ipdbtraj,iostat=iomess) (coords(j),j=1,n)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  if ((dofatal.EQV..false.).AND.(iomess.ne.0)) then
    mdl = -mdl
    return
  end if
  if (use_pdb_template.EQV..true.) then
    y(1:n) = 1.0*coords(pdbmap(1:n))
  else
    y(1:n) = 1.0*coords(1:n)
  end if
  read(ipdbtraj,iostat=iomess) (coords(j),j=1,n)
  if (iomess.ne.0) call fail_readdcd(iomess,dofatal)
  if ((dofatal.EQV..false.).AND.(iomess.ne.0)) then
    mdl = -mdl
    return
  end if
  if (use_pdb_template.EQV..true.) then
    z(1:n) = 1.0*coords(pdbmap(1:n))
  else
    z(1:n) = 1.0*coords(1:n)
  end if
  deallocate(coords)
!
! now generate appropriate internal coordinates from the read-in Cartesian values
  do imol=1,nmol
    call genzmat(imol)
  end do
!     
! and transfer them into FMSMC-arrays (which are the ones used in analysis functions)
  call zmatfyc2()
  if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
    do imol=1,nmol
      call makexyz_formol(imol)
    end do
  end if
!
! and finally update rigid-body coordinates
  do imol=1,nmol
    call update_rigid(imol)
    call update_rigidm(imol)
  end do
  if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
    if (align%yes.EQV..true.) then
      do imol=1,nmol
        call update_image2(imol)
      end do
    else
      do imol=1,nmol
        call update_image(imol)
      end do
    end if
  end if
!
end
!
!-------------------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
! see above
!
subroutine FMSMC_readnetcdf(mdl)
! 
  use iounit
  use pdb
  use atoms
  use mcsums
  use molecule
  use zmatrix
  use sequen
  use netcdf
  use clusters, ONLY: align
!
  implicit none
!
  integer i,mdl,imol
  real(KIND=4), ALLOCATABLE:: coords(:)
!  RTYPE box1(3), box2(3)
!
! for now we will not even read out the box coordinates since they would be flushed anyway
! WARNING: this needs to change for NPT support
!
  allocate(coords(3*n))
  call check_fileionetcdf( nf90_get_var(ipdbtraj, netcdf_ids(32), coords(1:n), &
 &                                       start = (/ 1,1,mdl /) , count = (/ 1,n,1 /)) )
  call check_fileionetcdf( nf90_get_var(ipdbtraj, netcdf_ids(32), coords(n+1:2*n),&
 &                                        start = (/ 2,1,mdl /) , count = (/ 1,n,1 /)) )
  call check_fileionetcdf( nf90_get_var(ipdbtraj, netcdf_ids(32), coords(2*n+1:3*n),&
 &                                        start = (/ 3,1,mdl /) , count = (/ 1,n,1 /)) )
!
  if (use_pdb_template.EQV..true.) then
   do i=1,n
      x(i) = 1.0*coords(pdbmap(i))
      y(i) = 1.0*coords(n+pdbmap(i))
      z(i) = 1.0*coords(2*n+pdbmap(i))
    end do
  else
    x(1:n) = 1.0*coords(1:n)
    y(1:n) = 1.0*coords(n+1:2*n)
    z(1:n) = 1.0*coords(2*n+1:3*n)
  end if
  deallocate(coords)
!
! now generate appropriate internal coordinates from the read-in Cartesian values
  do imol=1,nmol
    call genzmat(imol)
  end do
!     
! and transfer them into FMSMC-arrays (which are the ones used in analysis functions)
  call zmatfyc2()
  if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
    do imol=1,nmol
      call makexyz_formol(imol)
    end do
  end if
!
! and finally update rigid-body coordinates
  do imol=1,nmol
    call update_rigid(imol)
    call update_rigidm(imol)
  end do
  if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
    if (align%yes.EQV..true.) then
      do imol=1,nmol
        call update_image2(imol)
      end do
    else
      do imol=1,nmol
        call update_image(imol)
      end do
    end if
  end if
!
end
!
#endif
!
!------------------------------------------------------------------------------------
!
subroutine infer_from_pdb(mode)
!
  use pdb
  use iounit
  use sequen
  use polypep
  use atoms
  use aminos
  use params
  use zmatrix
  use molecule
  use math
  use interfaces
  use system
!
  implicit none
!
  integer ipdb,imol,nl,i,j,k,ats,ate,freeunit,rs,t1,t3,ii,currs,lastrs,mode,st1,st2,rsl,rsh,jjj,rssz,jmol
  integer stat,sta,sto,azero,rsshft,unkcnt,jupp,jj,kk,flagrsshft,firs,btext,seqt_from_name,match(10)
  character(MAXSTRLEN) string
  character(MAXSTRLEN) pdbsnap
  character(MAXSTRLEN+50) fpdbstr
  character(6) keyw
  character(3) resname
  character(MAXSTRLEN), ALLOCATABLE:: pdbl(:)
  character(6), ALLOCATABLE:: kw(:)
  character(4), ALLOCATABLE:: pdbnm(:)
  character(3), ALLOCATABLE:: rsnam(:),bexc(:)
  character(25) atnmb
  character(1) extrabit,refbit,extrabit2,refbit2
  character, ALLOCATABLE:: chnam(:)
  integer, ALLOCATABLE:: rsnmb(:),rsshftl(:),iv1(:),iv2(:),iv3(:),lix(:),seqinfile(:,:)
  logical exists,ribon,foundp
  logical, ALLOCATABLE:: lwarn(:)
  RTYPE, ALLOCATABLE:: xyzdum(:,:),d2max(:)
  RTYPE d2
!
  azero = 0
!
  if (pdb_readmode.eq.1) then
    write(ilog,*) 'Warning. Runs with unsupported residues require the direct read-in of Cartesian coordinates &
 &from pdb-file(s). Changing assumed setting for FMCSC_PDB_READMODE.'
    if (pdb_analyze.EQV..true.) then
      pdb_readmode = 2
    else
      pdb_readmode = 3
    end if
  end if
  ipdb = freeunit()
  if ((pdb_analyze.EQV..true.).AND.(pdb_fileformat.ge.3).AND.(pdb_fileformat.le.5)) then
    if (use_pdb_template.EQV..false.) then
      write(ilog,*) 'Fatal. Analysis runs on binary trajectory files with unsupported residues require the use of an auxiliary &
 &pdb-file (via FMCSC_PDB_TEMPLATE) for inference.'
      call fexit()
    else
      pdbsnap = pdbtmplfile
    end if
  else if (pdb_analyze.EQV..false.) then
    if ((use_pdb_template.EQV..true.).AND.(pdbinput.EQV..true.)) then
      if (mode.eq.1) then
        write(ilog,*) 'Warning. Specified both a template and an input structure for inferring topology &
 &in a run featuring unsupported residues. Using template only to infer topology.'
      end if
      pdbsnap = pdbtmplfile
    else if (use_pdb_template.EQV..true.) then
      pdbsnap = pdbtmplfile
    else if (pdbinput.EQV..true.) then
      pdbsnap = pdbinfile
    else
      write(ilog,*) 'Fatal. Runs with unsupported residues require the use of an auxiliary &
 &pdb-file (either via FMCSC_PDB_TEMPLATE or via FMCSC_PDBFILE) for inference.'
      call fexit()
    end if
  else
    if (use_pdb_template.EQV..true.) then
      if (mode.eq.1) then
        write(ilog,*) 'Warning. Specified both a template and a pdb trajectory for inferring topology &
 &in an analysis run featuring unsupported residues. Using template only to infer topology.'
      end if
      pdbsnap = pdbtmplfile
    else
      pdbsnap = pdbinfile
    end if
  end if
  call strlims(pdbsnap,t1,t3)
  inquire (file=pdbsnap(t1:t3),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open pdb-file required to support unknown residues in seq. file (&
 &provided: ',pdbsnap(t1:t3),').'
    call fexit()
  end if
  open (unit=ipdb,file=pdbsnap(t1:t3),status='old')
!
  k = 0
  do while (1.eq.1)
    k = k + 1
 1356 format(FORM_MAXSTRLEN)
    read(ipdb,1356,err=30,end=30) string   
  end do
 30   k = k - 1
  close(unit=ipdb)
  allocate(pdbl(k+10))
!
  open (unit=ipdb,file=pdbsnap(t1:t3),status='old')
!
  k = 0
  do while (1.eq.1)
    k = k + 1
    read(ipdb,1356,err=31,end=31) pdbl(k)        
  end do
 31   k = k - 1
  close(unit=ipdb)
!
  allocate(pdbnm(k+10))
  allocate(chnam(k+10))
  allocate(rsnam(k+10))
  allocate(kw(k+10))
  allocate(rsnmb(k+10))
  allocate(lix(k+10))
  allocate(rsshftl(k+10))
  allocate(xyzdum(3,k+10))
!
  nl = k
  sta = -1
  sto = -1
  do i=1,nl
    read(pdbl(i),'(a6)') keyw
    if ((keyw.eq.'MODEL ').AND.(sta.eq.-1)) then
      sta = i
    else if ((keyw.eq.'ENDMDL').AND.(sta.ne.-1)) then
      sto = i
      exit
    else if (keyw.eq.'ENDMDL') then
      write(ilog,*) 'Fatal. Encountered ENDMDL before MODEL-keywor&
 &d while reading from PDB (',pdbsnap(t1:t3),').'
      call fexit()
    end if
  end do
!
  ats = -1
  ate = -1
  if ((sta.ne.-1).AND.(sto.ne.-1)) then
    do i=sta,sto
      read(pdbl(i),'(a6)') keyw
      if (((keyw.eq.'ATOM  ').OR.&
 &         (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
        ats = i
      else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &         (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
        ate = i-1
        exit
      end if
    end do
  else if ((sta.eq.-1).AND.(sto.eq.-1)) then
    do i=1,nl
      read(pdbl(i),'(a6)') keyw
      if (((keyw.eq.'ATOM  ').OR.&
 &         (keyw.eq.'HETATM')).AND.(ats.eq.-1)) then
        ats = i
      else if (((keyw.ne.'ATOM  ').AND.(keyw.ne.'TER   ').AND.&
 &         (keyw.ne.'HETATM').AND.(keyw.ne.'ANISOU')).AND.(ats.ne.-1)) then
        ate = i-1
        exit
      end if
    end do
  else
    write(ilog,*) 'Fatal. Could not find ENDMDL keyword despite MO&
 &DEL keyword while reading from PDB (',pdbsnap(t1:t3),').'
    call fexit()
  end if
  if ((ats.le.0).AND.(ate.le.0)) then
    write(ilog,*) 'Fatal. PDB file (',pdbsnap(t1:t3),') appears to contain no coordinate information whatsoever. &
 &Please check input.'
    call fexit()
  end if
!
  call strlims(pdb_formstr(2),st1,st2)
  fpdbstr = '('//pdb_formstr(2)(st1:st2)//')'
!
! this gives us the frame, within which to operate
  k = 0
  rsshft = 0
  rsshftl(:) = 0
  currs = huge(currs)
  refbit = ' '
  refbit2 = ' '
  do i=ats,ate
    k = k + 1
    read(pdbl(i),'(a6)') keyw
    if ((keyw.eq.'TER   ').OR.(keyw.eq.'ANISOU')) then
      k = k - 1
      cycle
    end if
    read(pdbl(i),fpdbstr) kw(k),atnmb,pdbnm(k),extrabit,rsnam(k),chnam(k),rsnmb(k),extrabit2&
 &,xyzdum(1,k),xyzdum(2,k),xyzdum(3,k)
    lix(k) = i
!    
!   the currs/lastrs construct is (at least sporadically) able to read in pdb-files in which
!   the residue numbering is off as long as it is a number and the number for a new residue is always
!   different from the previous one
    if (k.eq.1) lastrs = rsnmb(k)-1
    if (currs.ne.huge(currs)) then
      if ((rsnmb(k).ne.currs).OR.(extrabit2.ne.refbit2)) then
        refbit = ' '
        if (currs.ne.lastrs+1) currs = lastrs+1
        lastrs = currs
        currs = rsnmb(k)
        if (rsnmb(k).ne.lastrs+1) then
          rsnmb(k) = lastrs + 1
        end if
      else
        if (rsnmb(k).ne.lastrs+1) rsnmb(k) = lastrs+1
      end if
    end if
    if (k.eq.1) currs = rsnmb(k)
    if ((refbit.eq.' ').AND.(extrabit.ne.' ')) refbit = extrabit
    refbit2 = extrabit2
    call toupper(kw(k))
!   drop the records for TER
    if ((kw(k).eq.'TER   ').OR.(kw(k).eq.'ANISOU').OR.((extrabit.ne.' ').AND.(extrabit.ne.refbit))) then
      k = k - 1
      cycle
    end if
    call toupper(pdbnm(k))
    call toupper(rsnam(k))
    call toupper(chnam(k))
    call pdbcorrect(rsnam(k),pdbnm(k),rsnmb(k),rsshft,rsshftl)
  end do
!
  nl = k
!
  ii = 1
  do j=1,rsshft
    flagrsshft = 0
    do i=ii,k
      if (rsnmb(i).eq.rsshftl(j)) then
        do jj=i,k
          if (rsnmb(jj).ne.rsshftl(j)) then
            ii = jj
            exit
          end if
          if ((rsnam(jj).eq.'NH2').OR.(rsnam(jj).eq.'NME')) then
            rsnmb(jj) = rsnmb(jj) + 1
            if (flagrsshft.eq.0) flagrsshft = 1
            if (flagrsshft.eq.-1) flagrsshft = 2
          end if
          if ((rsnam(jj).eq.'FOR').OR.(rsnam(jj).eq.'ACE')) then
            rsnmb(jj) = rsnmb(jj) - 1
            if (flagrsshft.eq.0) flagrsshft = -1
            if (flagrsshft.eq.1) flagrsshft = 2
          end if
          if (jj.eq.k) ii = k + 1
        end do
        if (flagrsshft.eq.1) then
          if (ii.le.k) then
            rsnmb(ii:k) = rsnmb(ii:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 1
          end if
        else if (flagrsshft.eq.-1) then
          if (i.le.k) then
            rsnmb(i:k) = rsnmb(i:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 1
          end if
        else if (flagrsshft.eq.2) then
          if (i.le.k) then
            rsnmb(i:k) = rsnmb(i:k) + 1
          end if
          if (ii.le.k) then
            rsnmb(ii:k) = rsnmb(ii:k) + 1
          end if
          if (j.lt.rsshft) then
            rsshftl(j+1:rsshft) = rsshftl(j+1:rsshft) + 2
          end if
        end if
        exit ! skip out, go to next shft
      end if
    end do
  end do
!
  jupp = -1000
  do i=1,k
    if (rsnmb(i).ne.jupp) then
      if (i.gt.1) then
        ribon = .false.
        foundp = .false.
        do j=firs,i-1
          if (pdbnm(j).eq.' P  ') then ! crucial that P is always called P
            foundp = .true.
            exit
          end if
        end do
        do j=firs,i-1
          if (rsnam(j)(1:2).eq.'R ') then
            ribon = .true.
          end if
        end do
        if (ribon.EQV..true.) then
          do j=firs,i-1
!           hopeful
            if (foundp.EQV..true.) then
              rsnam(j)(1:2) = 'RP'
            else
              rsnam(j)(1:2) = 'RI'
            end if
            if (pdb_convention(2).eq.3) then
              if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = '2HO*'
              if (pdbnm(j)(1:4).eq.'2H2*') pdbnm(j)(1:4) = ' H2*'
            else
              if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = ' H2*'
            end if
          end do
        else if (rsnam(firs)(1:2).eq.'DP') then
          if (foundp.EQV..false.) then
            do j=firs,i-1
              rsnam(j)(2:2) = 'I'
            end do
          end if
        end if
        if ((ribon.EQV..true.).OR.(rsnam(firs)(1:2).eq.'DP')) then
          if (foundp.EQV..true.) then
            do j=firs,i-1          
              if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = '5HO*'
            end do
          else
            do j=firs,i-1  
              if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = ' HOP'
            end do
          end if
        end if
      end if
      firs = i
      jupp = rsnmb(i)
    end if
    if (i.eq.k) then
      ribon = .false.
      foundp = .false.
      do j=firs,i-1
        if (pdbnm(j).eq.' P  ') then ! crucial that P is always called P
          foundp = .true.
          exit
        end if
      end do
      do j=firs,k
        if (rsnam(j)(1:2).eq.'R ') then
          ribon = .true.
        end if
      end do
      if (ribon.EQV..true.) then
        do j=firs,k
          if (foundp.EQV..true.) then
            rsnam(j)(1:2) = 'RP'
          else
            rsnam(j)(1:2) = 'RI'
          end if
          if (pdb_convention(2).eq.3) then
            if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = '2HO*'
            if (pdbnm(j)(1:4).eq.'2H2*') pdbnm(j)(1:4) = ' H2*'
          else
            if (pdbnm(j)(1:4).eq.'1H2*') pdbnm(j)(1:4) = ' H2*'
          end if
        end do
      else if (rsnam(firs)(1:2).eq.'DP') then
        if (foundp.EQV..false.) then
          do j=firs,k
            rsnam(j)(2:2) = 'I'
          end do
        end if
      end if
      if ((ribon.EQV..true.).OR.(rsnam(firs)(1:2).eq.'DP')) then
        if (foundp.EQV..true.) then
          do j=firs,i-1          
            if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = '5HO*'
          end do
        else
          do j=firs,i-1  
            if (pdbnm(j)(1:4).eq.'XHOP') pdbnm(j)(1:4) = ' HOP'
          end do
        end if
      end if
    end if
  end do
  call adjust_PONuc(nl,rsnmb(1:nl),pdbnm(1:nl),rsnam(1:nl),xyzdum(1,1:nl),xyzdum(2,1:nl),xyzdum(3,1:nl))
!
! this function reorders the input arrays (except pdbl) so that everything is contiguous per residue
  call reorder_PDB_byrsnmb(nl,kw(1:nl),rsnmb(1:nl),rsnam(1:nl),pdbnm(1:nl),xyzdum(1,1:nl),xyzdum(2,1:nl),xyzdum(3,1:nl),&
 &lix(1:nl),chnam(1:nl))
!
  if (mode.eq.1) then
!
    if (allocated(pdb_unkbnd).EQV..true.) then
      deallocate(pdb_unkbnd)
    end if
    allocate(pdb_unkbnd(n_pdbunk,4))
    pdb_unkbnd(:,:) = 0
    kk = 0
    do j=1,nseq
      if (seqtyp(j).le.0) then
        kk = kk + 1
        pdb_unkbnd(kk,3) = j
      end if
    end do
    kk = 1
    allocate(lwarn(n_pdbunk))
    lwarn(:) = .false.
    jj = 1
    do i=1,nl
      foundp = .false.
      if (i.eq.1) foundp = .true.
      if (i.gt.1) then
        if (rsnmb(i).ne.rsnmb(i-1)) then
          foundp = .true.
        end if
      end if
      if (foundp.EQV..true.) then
        if (pdb_unkbnd(kk,1).gt.0) then
          lwarn(kk) = .true.
          pdb_unkbnd(kk,2) = i-1
          kk = kk + 1
          if (kk.gt.n_pdbunk) exit
        end if
        if ((rsnam(i).ne.pdb_unknowns(kk)).AND.(rsmol(molofrs(pdb_unkbnd(kk,3)),1).eq.rsmol(molofrs(pdb_unkbnd(kk,3)),2)).AND.&
 &          (seqt_from_name(rsnam(i),jj).le.0)) then
!         we have a mismatch, which can mean context or superfluous information in the pdb but also
!         duplicated unsupported small molecules in the sequence -> check
          exists = .true.
          do while (exists.EQV..true.)
            exists = .false.
            do j=1,kk-1
              if ((pdb_unknowns(kk).eq.pdb_unknowns(j)).AND.(lwarn(j).EQV..true.)) then
                exists = .true.
                exit
              end if
            end do
            if (exists.EQV..false.) exit
            kk = kk + 1
            if (kk.gt.n_pdbunk) exit
            if (rsmol(molofrs(pdb_unkbnd(kk,3)),1).ne.rsmol(molofrs(pdb_unkbnd(kk,3)),2)) exit
            if (rsnam(i).eq.pdb_unknowns(kk)) exit
          end do
          if (kk.gt.n_pdbunk) exit
        end if
      end if
      if ((i.eq.nl).AND.(pdb_unkbnd(kk,1).gt.0).AND.(pdb_unkbnd(kk,2).le.0)) then ! last residue must be able to close
        lwarn(kk) = .true.
        pdb_unkbnd(kk,2) = i
      else if (foundp.EQV..true.) then
        if (rsnam(i).eq.pdb_unknowns(kk)) then
          pdb_unkbnd(kk,1) = i
        end if
      end if
    end do
    lwarn(:) = .false.
    do i=n_pdbunk,1,-1
      if ((pdb_unkbnd(i,1).le.0).OR.(pdb_unkbnd(i,2).le.0)) then
        foundp = .false.
        jj = molofrs(pdb_unkbnd(i,3))
        if (rsmol(jj,1).eq.rsmol(jj,2)) then
          do j=1,i-1
            if (pdb_unknowns(i).eq.pdb_unknowns(j)) then
              if (minval(pdb_unkbnd(j,1:2)).gt.0) then
                if (lwarn(j).EQV..false.) then
                  write(ilog,'(1x,a,a3,a,a,a)') 'Remark. Template for unsupported small molecule ',pdb_unknowns(i),' is used &
 &repeatedly from input file ',pdbsnap(t1:t3),'.'
                end if
                lwarn(j) = .true.
                pdb_unkbnd(i,1) = -1
                pdb_unkbnd(i,2) = j
                foundp = .true.
                exit
              end if
            end if       
          end do
        end if
      end if
      if ((pdb_unkbnd(i,1).eq.0).OR.(pdb_unkbnd(i,2).eq.0)) then
        write(ilog,'(1x,3(a),i7,3(a))') 'Fatal. Cannot find coordinates for unsupported residue with name ',pdb_unknowns(i),' &
 &(sequence position ',pdb_unkbnd(i,3),') in &
 &input PDB file (',pdbsnap(t1:t3),'). This may be related to the order in which residues appear in the template. Only &
 &single-residue molecules can be duplicated automatically from the template.'
        do j=1,i-1
          if ((pdb_unkbnd(j,1).eq.0).OR.(pdb_unkbnd(j,2).eq.0)) exit
        end do
        if (i.ne.j) then
          write(ilog,'(1x,3(a),i7,a)') 'The first unsupported residue to miss information has name ',pdb_unknowns(j),' &
 &and sequence position ',pdb_unkbnd(j,3),'.'
        end if
        call fexit()
      end if
    end do
    deallocate(lwarn)
    allocate(d2max(nl))
    do i=1,n_pdbunk
      rs = pdb_unkbnd(i,3)
      d2max(:) = 1.5 ! not all values are set later
      if (pdb_unkbnd(i,1).eq.-1) then
        jj = pdb_unkbnd(i,2)
        natres(rs) = pdb_unkbnd(jj,2)-pdb_unkbnd(jj,1)+1
      else
        natres(rs) = pdb_unkbnd(i,2)-pdb_unkbnd(i,1)+1
        jj = i
      end if 
      do j=pdb_unkbnd(jj,1),pdb_unkbnd(jj,2)
        d2max(j) = 0.0
        do k=pdb_unkbnd(jj,1),pdb_unkbnd(jj,2)
          call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2)
          if (d2.gt.d2max(j)) d2max(j) = d2
        end do
      end do
      d2 = HUGE(d2)
      do j=pdb_unkbnd(jj,1),pdb_unkbnd(jj,2)
        if (d2max(j).lt.d2) then
          d2 = d2max(j)
          refat(rs) = j ! temporary assignment -> may need to be shifted in call with mode 2!
          resrad(rs) = 1.5*sqrt(d2max(j))
        end if
      end do
    end do
!
  else if (mode.eq.2) then
!   now let's locate the sequence
    rssz = rsnmb(nl)-rsnmb(1)+1
    kk = max(rssz,maxval(molofrs(pdb_unkbnd(:,3))))
    allocate(seqinfile(kk,6))
    seqinfile(:,:) = 0
    seqinfile(1,2) = 1
    j = 0 
    jj = 1
    do i=2,nl
      if (rsnmb(i).ne.rsnmb(i-1)) then
        j = j + 1
        seqinfile(j,1) = seqt_from_name(rsnam(i-1),jj)
        seqinfile(j+1,2) = i
      end if
    end do
    j = j + 1
    seqinfile(j,1) = seqt_from_name(rsnam(k),jj)
!
    do unkcnt=1,n_pdbunk
      rs = pdb_unkbnd(unkcnt,3)
      seqtyp(rs) = -unkcnt ! so comparison works
    end do
!
!   initialization and prep
    ribon = use_pdb_template
    use_pdb_template = .true.
    allocate(pdbmap(max(nl,n)))
    pdbmap(:) = 0
    allocate(xpdb(max(nl,n)))
    allocate(ypdb(max(nl,n)))
    allocate(zpdb(max(nl,n)))
    xpdb(1:n) = x(1:n)
    ypdb(1:n) = y(1:n)
    zpdb(1:n) = z(1:n)
    allocate(pdbmap2(n))
    stat = 1
!
    do unkcnt=1,n_pdbunk
      match(:) = 0
      rs = pdb_unkbnd(unkcnt,3)
!     shift refat
      if (pdb_unkbnd(unkcnt,1).eq.-1) then
        refat(rs) = refat(rs) + at(rs)%bb(1) - pdb_unkbnd(pdb_unkbnd(unkcnt,2),1)
      else
        refat(rs) = refat(rs) + at(rs)%bb(1) - pdb_unkbnd(unkcnt,1)
      end if
      imol = molofrs(rs)
      rsl = max(rsmol(imol,1),rs-1)
      rsh = min(rsmol(imol,2),rs+1)
      jj = stat
      j = rsl
      if (pdb_unkbnd(unkcnt,1).eq.-1) then
!       small molecule repeated -> take information from first instance
        jmol = molofrs(pdb_unkbnd(pdb_unkbnd(unkcnt,2),3))
        j = -1
        do i=seqinfile(jmol,3),seqinfile(jmol,4)
          j = j + 1
          call trfxyz(i,at(rs)%bb(1)+j)
        end do
        cycle
      end if
      do while (jj.le.rssz)
        if (seqtyp(rsl).eq.seqinfile(jj,1)) then
          jjj = jj
          do while ((j.le.rsh).AND.(jjj.le.rssz))
            if ((seqtyp(j).ne.seqinfile(jjj,1)).OR.(j.eq.rsh).OR.(jjj.eq.rssz)) then
              if (seqtyp(j).ne.seqinfile(jjj,1)) then
                j = rsl
                exit
              end if
              if (j.lt.rsh) then
                write(ilog,*) 'Fatal. Structural input from PDB file ',pdbsnap(t1:t3),' must not end in an unsupported &
 &residue that is not also the last residue in its molecule.'
                call fexit()
              end if
              match(1) = rsl
              match(2) = rsh-rsl+1
              match(3) = jj
              exit
            end if
            jjj = jjj + 1
            j = j + 1
          end do
        end if
        if (match(1).gt.0) exit
        jj = jj + 1
      end do
      if (match(1).le.0) then
        write(ilog,'(1x,a,a,a,a,i7,a)') 'Fatal. Input PDB file (',pdbsnap(t1:t3),') does not appear to contain exact sequence &
 &context required to infer topology of all unsupported residues (failed for ',pdb_unknowns(unkcnt),pdb_unkbnd(unkcnt,3),').'
        call fexit()
      end if
      seqinfile(imol,3) = seqinfile(match(3),2)
      if ((match(3)+match(2)).gt.rssz) then
        seqinfile(imol,4) = k
      else
        seqinfile(imol,4) = seqinfile(match(3)+match(2),2)-1
      end if
      seqinfile(imol,5) = match(1)
      seqinfile(imol,6) = match(1)+match(2)-1
      stat = match(3) + 1
      if ((rs.eq.rsl).AND.(rsh.gt.rsl)) stat = match(3)
 456 format(1x,'Unsupported residue ',i7,' was found and matched in a stretch &
 &enclosed by the following PDB records:'/'====================================================================='/(3x,a)/(3x,a)/&
 &         '=====================================================================')
      write(ilog,456) pdb_unkbnd(unkcnt,3),trim(pdbl(lix(seqinfile(imol,3)))),trim(pdbl(lix(seqinfile(imol,4))))
!
      rs = rsl
      if (seqtyp(rs).gt.0) then
        resname = amino(seqtyp(rs))
      else
        resname = 'UNK'
      end if
      j = -1
      do i=seqinfile(imol,3),seqinfile(imol,4)
        if (i.gt.seqinfile(imol,3)) then
          if (rsnmb(i).ne.rsnmb(i-1)) then
            rs = rs + 1
            if (seqtyp(rs).gt.0) then
              resname = amino(seqtyp(rs))
            else
              resname = 'UNK'
            end if
            j = -1
          end if
        end if
!
        if (seqtyp(rs).lt.0) then
          j = j + 1
          call trfxyz(i,at(rs)%bb(1)+j)
          cycle
        end if
!
!       single residues molecules first
        if (rsmol(molofrs(rs),1).eq.rsmol(molofrs(rs),2)) then
          call trf_sms(i,rs,resname,pdbnm(i))
!       now the backbone of N-terminal residues
        else if (rs.eq.rsmol(molofrs(rs),1)) then
          call trf_nts(i,rs,resname,pdbnm(i))
!       now the backbone of C-terminal residues
        else if (rs.eq.rsmol(molofrs(rs),2)) then
          call trf_cts(i,rs,resname,pdbnm(i))
!       finally, regular in-the-middle-of-the-chain residues
        else
          call trf_bbs(i,rs,resname,pdbnm(i))
        end if
!       now the sidechains
        call trf_scs(i,rs,resname,pdbnm(i))
      end do
    end do
    deallocate(seqinfile)
!
    allocate(d2max(nl))
!
    btext = 0
    do i=1,n_pdbunk
      rs = pdb_unkbnd(i,3)
      imol = molofrs(rs)
      rsl = max(rsmol(imol,1),rs-1)
      rsh = min(rsmol(imol,2),rs+1)
      foundp = .false.
      do j=at(rsl)%bb(1),at(rsh)%bb(1)+at(rsh)%nbb+at(rsh)%nsc-1
        if ((pdbmap(j).le.0).AND.(foundp.EQV..false.)) then
          write(ilog,*) 'Warning. When using unsupported residues, it is strongly recommended that &
 &atoms in those supported residues in the input file (',pdbsnap(t1:t3),'), which are directly adjacent and in the &
 &same molecule, can be matched exactly to internal CAMPARI representation.'
          write(ilog,'(1x,a,i7,3(a))') 'A fatal error may occur below. Failed reads for position ',rs,' (',pdb_unknowns(i),') are:'
          foundp = .true.
        end if
        if (pdbmap(j).le.0) then
          write(ilog,'(1x,a,i8,3(a),i7)') 'Atom ',j,' with biotype name ',bio_code(b_type(j)),' in residue ',atmres(j)
        end if
      end do
!      if (foundp.EQV..true.) call fexit()
    end do
    do i=1,n_pdbunk
      rs = pdb_unkbnd(i,3)
      foundp = .false.
      if (pdb_unkbnd(i,1).gt.0) then
        do j=1,i-1
          if (pdb_unkbnd(j,1).le.0) cycle
          if (rsnam(pdb_unkbnd(j,1))(1:3).eq.rsnam(pdb_unkbnd(i,1))(1:3)) then
            foundp = .true.
            exit
          end if
        end do
        if (foundp.EQV..false.) btext = btext + pdb_unkbnd(i,2) - pdb_unkbnd(i,1) + 1
      end if
!
!     we want to hijack atnam(j)(3:3), so keep clean
      if (pdb_unkbnd(i,1).eq.-1) then
        jj = pdb_unkbnd(i,2)
        pdb_unkbnd(i,1:2) = pdb_unkbnd(jj,1:2)
      end if
      do jj=pdb_unkbnd(i,1),pdb_unkbnd(i,2)
        j = at(rs)%bb(1) + jj - pdb_unkbnd(i,1)
        if ((pdbnm(jj)(1:3).eq.'BR ').OR.(pdbnm(jj)(1:3).eq.' BR')) then
          mass(j) = 79.904
          d2max(jj) = 2.6
          atnam(j) = 'BR '
        else if ((pdbnm(jj)(1:3).eq.'AL ').OR.(pdbnm(jj)(1:3).eq.' AL')) then
          mass(j) = 26.98
          d2max(jj) = 2.6
          atnam(j) = 'AL '
        else if (pdbnm(jj)(1:3).eq.'CU ') then ! avoid clashes with C*
          mass(j) = 63.55
          d2max(jj) = 2.8
          atnam(j) = 'CU '
        else if ((pdbnm(jj)(1:3).eq.'MG ').OR.(pdbnm(jj)(1:3).eq.' MG')) then
          mass(j) = 24.31
          d2max(jj) = 3.0
          atnam(j) = 'MG '
        else if (pdbnm(jj)(1:3).eq.'CA ') then ! avoid clashes with C*
          mass(j) = 40.08
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'CA '
        else if ((pdbnm(jj)(1:3).eq.'FE ').OR.(pdbnm(jj)(1:3).eq.' FE')) then
          mass(j) =  55.85
          d2max(jj) = 2.8
          atnam(j) = 'FE '
        else if ((pdbnm(jj)(1:3).eq.'F  ').OR.(pdbnm(jj)(1:3).eq.' F ')) then
          mass(j) = 18.998
          d2max(jj) = 1.55
          atnam(j) = 'F  '
        else if ((pdbnm(jj)(1:3).eq.'I  ').OR.(pdbnm(jj)(1:3).eq.' I ')) then
          mass(j) = 126.905
          d2max(jj) = 3.0
          atnam(j) = 'I  '
        else if ((pdbnm(jj)(1:3).eq.'MN ').OR.(pdbnm(jj)(1:3).eq.' MN')) then
          mass(j) = 54.94
          d2max(jj) = 3.0
          atnam(j) = 'MN '
        else if (pdbnm(jj)(1:3).eq.'CO ') then ! avoid clashes with C*
          mass(j) = 58.93
          d2max(jj) = 2.7 
          atnam(j) = 'CO '
        else if (pdbnm(jj)(1:3).eq.'CR ') then ! avoid clashes with C*
          mass(j) = 52.00
          d2max(jj) = 3.0
          atnam(j) = 'CR '
        else if ((pdbnm(jj)(1:3).eq.'ZN ').OR.(pdbnm(jj)(1:3).eq.' ZN')) then
          mass(j) = 65.41
          d2max(jj) = 2.6
          atnam(j) = 'ZN '
        else if (pdbnm(jj)(1:3).eq.'SI ') then ! avoid clashes with S*
          mass(j) = 28.09
          d2max(jj) = 2.5
          atnam(j) = 'SI '
        else if ((pdbnm(jj)(1:3).eq.'AS ').OR.(pdbnm(jj)(1:3).eq.' AS')) then
          mass(j) = 74.92
          d2max(jj) = 2.6
          atnam(j) = 'AS '
        else if (pdbnm(jj)(1:3).eq.'SE ') then ! avoid clashes with S*
          mass(j) = 78.96
          d2max(jj) = 2.6
          atnam(j) = 'SE '
        else if ((pdbnm(jj)(1:3).eq.'LI ').OR.(pdbnm(jj)(1:3).eq.' LI')) then
          mass(j) = 6.94
          d2max(jj) = 2.8 ! probably too small
          atnam(j) = 'LI '
        else if ((pdbnm(jj)(1:3).eq.'K  ').OR.(pdbnm(jj)(1:3).eq.' K ')) then
          mass(j) = 39.098
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'K  '
        else if (pdbnm(jj)(1:3).eq.'NA ') then ! avoid clashes with N*
          mass(j) = 22.990
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'NA '
        else if ((pdbnm(jj)(1:3).eq.'CL ').OR.(pdbnm(jj)(1:4).eq.' CL ')) then
          mass(j) = 35.453
          d2max(jj) = 2.2
          atnam(j) = 'CL '
        else if (pdbnm(jj)(1:2).eq.' P') then
          mass(j) = 30.974
          atnam(j) = 'P  '
          d2max(jj) = 2.0
        else if (pdbnm(jj)(1:2).eq.' S') then
          mass(j) = 32.066
          atnam(j) = 'S  '
          d2max(jj) = 2.1
        else if (pdbnm(jj)(1:2).eq.' O') then
          mass(j) = 15.999
          atnam(j) = 'O  '
          d2max(jj) = 1.6
        else if (pdbnm(jj)(1:2).eq.' N') then
          mass(j) = 14.007
          atnam(j) = 'N  '
          d2max(jj) = 1.6
        else if (pdbnm(jj)(1:2).eq.' C') then
          mass(j) = 12.011
          atnam(j) = 'C  '
          d2max(jj) = 1.7
        else if (pdbnm(jj)(1:2).eq.' H') then
          mass(j) = 1.008
          atnam(j) = 'H  '
          d2max(jj) = 1.3
        else if (pdbnm(jj)(1:2).eq.' F') then
          mass(j) = 18.998
          d2max(jj) = 1.55
          atnam(j) = 'F  '
        else if (pdbnm(jj)(1:2).eq.' I') then
          mass(j) = 126.905
          d2max(jj) = 3.0
          atnam(j) = 'I  '
        else if (pdbnm(jj)(1:2).eq.'CA') then ! avoid clashes with C*
          mass(j) = 40.08
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'CA '
        else if (pdbnm(jj)(1:2).eq.'CL') then
          mass(j) = 35.453
          d2max(jj) = 2.2
          atnam(j) = 'CL '
        else if (pdbnm(jj)(2:2).eq.'P') then
          mass(j) = 30.974
          atnam(j) = 'P  '
          d2max(jj) = 2.0
        else if (pdbnm(jj)(2:2).eq.'S') then
          mass(j) = 32.066
          atnam(j) = 'S  '
          d2max(jj) = 2.1
        else if (pdbnm(jj)(2:2).eq.'O') then
          mass(j) = 15.999
          atnam(j) = 'O  '
          d2max(jj) = 1.6
        else if (pdbnm(jj)(2:2).eq.'N') then
          mass(j) = 14.007
          atnam(j) = 'N  '
          d2max(jj) = 1.6
        else if (pdbnm(jj)(2:2).eq.'C') then
          mass(j) = 12.011
          atnam(j) = 'C  '
          d2max(jj) = 1.7
        else if (pdbnm(jj)(2:2).eq.'H') then
          mass(j) = 1.008
          atnam(j) = 'H  '
          d2max(jj) = 1.3
        else if (pdbnm(jj)(1:2).eq.'BR') then
          mass(j) = 79.904
          d2max(jj) = 2.6
          atnam(j) = 'BR '
        else if (pdbnm(jj)(1:2).eq.'AL') then
          mass(j) = 26.98
          d2max(jj) = 2.6
          atnam(j) = 'AL '
        else if (pdbnm(jj)(1:2).eq.'CU') then
          mass(j) = 63.55
          d2max(jj) = 2.8
          atnam(j) = 'CU '
        else if (pdbnm(jj)(1:2).eq.'MG') then
          mass(j) = 24.31
          d2max(jj) = 3.0
          atnam(j) = 'MG '
        else if (pdbnm(jj)(1:2).eq.'CA') then
          mass(j) = 40.08
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'CA '
        else if (pdbnm(jj)(1:2).eq.'FE') then
          mass(j) =  55.85
          d2max(jj) = 2.8
          atnam(j) = 'FE '
        else if (pdbnm(jj)(1:2).eq.'MN') then
          mass(j) = 54.94
          d2max(jj) = 3.0
          atnam(j) = 'MN '
        else if (pdbnm(jj)(1:2).eq.'CO') then
          mass(j) = 58.93
          d2max(jj) = 2.7 
          atnam(j) = 'CO '
        else if (pdbnm(jj)(1:2).eq.'CR') then
          mass(j) = 52.00
          d2max(jj) = 3.0
          atnam(j) = 'CR '
        else if (pdbnm(jj)(1:2).eq.'ZN') then
          mass(j) = 65.41
          d2max(jj) = 2.6
          atnam(j) = 'ZN '
        else if (pdbnm(jj)(1:2).eq.'SI') then
          mass(j) = 28.09
          d2max(jj) = 2.5
          atnam(j) = 'SI '
        else if (pdbnm(jj)(1:2).eq.'AS') then
          mass(j) = 74.92
          d2max(jj) = 2.6
          atnam(j) = 'AS '
        else if (pdbnm(jj)(1:2).eq.'SE') then
          mass(j) = 78.96
          d2max(jj) = 2.6
          atnam(j) = 'SE '
        else if (pdbnm(jj)(1:2).eq.'LI') then
          mass(j) = 6.94
          d2max(jj) = 2.8 ! probably too small
          atnam(j) = 'LI '
        else if (pdbnm(jj)(1:2).eq.'NA') then
          mass(j) = 22.990
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'NA '
        else if (pdbnm(jj)(1:1).eq.'F') then
          mass(j) = 18.998
          d2max(jj) = 1.55
          atnam(j) = 'F  '
        else if (pdbnm(jj)(1:1).eq.'I') then
          mass(j) = 126.905
          d2max(jj) = 3.0
          atnam(j) = 'I  '
        else if ((pdbnm(jj)(1:1).eq.'K').OR.(pdbnm(jj)(1:2).eq.' K')) then
          mass(j) = 39.098
          d2max(jj) = 3.0 ! max - should be higher
          atnam(j) = 'K  '
        else
          mass(j) = 10.0
          atnam(j) = '?  '
          d2max(jj) = 1.6
        end if
      end do
    end do
!
!   extend bio type list (note we are not extending bio_res; change if needed)
    allocate(bexc(n_biotyp))
    allocate(iv1(n_biotyp))
    allocate(iv2(n_biotyp))
    allocate(iv3(n_biotyp))
    iv1(:) = bio_ljtyp(:)
    iv2(:) = bio_ctyp(:)
    iv3(:) = bio_botyp(:)
    bexc(:) = bio_code(:)
    deallocate(bio_ljtyp)
    allocate(bio_ljtyp(n_biotyp+btext))
    bio_ljtyp(1:n_biotyp) = iv1(:)
    deallocate(bio_ctyp)
    allocate(bio_ctyp(n_biotyp+btext))
    bio_ctyp(1:n_biotyp) = iv2(:)
    deallocate(bio_botyp)
    allocate(bio_botyp(n_biotyp+btext))
    bio_botyp(1:n_biotyp) = iv3(:)
    deallocate(bio_code)
    allocate(bio_code(n_biotyp+btext))
    bio_code(1:n_biotyp) = bexc(:)
!
 666 format(' Reference atom and radius for unsupported residue ',i8,' (',a3,'/UNK) are ',i10,' and ',g8.3,'A respectively.') 
    do i=1,n_pdbunk
!
      rs = pdb_unkbnd(i,3)
!
      call infer_topology(i,nl,pdb_unkbnd,xyzdum(1:3,1:nl),d2max(1:nl),pdbnm(1:nl),rsnam(1:nl))
!
      call infer_polymer(i,nl,pdb_unkbnd,pdbnm(1:nl))
!
      seqtyp(rs) = 26
      foundp = .false.
      do j=1,i-1
        if ((pdb_unkbnd(i,1).eq.pdb_unkbnd(j,1)).AND.(pdb_unkbnd(i,2).eq.pdb_unkbnd(j,2))) then
          foundp = .true.
          exit
        end if
      end do
      if (foundp.EQV..false.) then
        write(ilog,666) rs,pdb_unknowns(i),refat(rs),resrad(rs)
        write(ilog,*)
      end if
    end do
!
    do i=1,n_pdbunk
      rs = pdb_unkbnd(i,3)
      call correct_topology(rs,nl,xyzdum(1:3,1:nl))
    end do
!
    deallocate(d2max)
    do i=1,n
      if (pdbmap(i).gt.0) then
        xref(i) = xyzdum(1,pdbmap(i))
        yref(i) = xyzdum(2,pdbmap(i))
        zref(i) = xyzdum(3,pdbmap(i))
      else
        xref(i) = x(i)
        yref(i) = y(i)
        zref(i) = z(i)
      end if
    end do
    deallocate(zpdb)
    deallocate(ypdb)
    deallocate(xpdb)
    deallocate(pdbmap)
    deallocate(pdbmap2)
    use_pdb_template = ribon
!
  end if
!
  deallocate(pdbl)
  deallocate(pdbnm)
  deallocate(chnam)
  deallocate(rsnam)
  deallocate(kw)
  deallocate(rsnmb)
  deallocate(rsshftl)
  deallocate(xyzdum)
  deallocate(lix)
!
end
!
!--------------------------------------------------------------------------------------------------
!
! this subroutine generally uses true numbering: the PDB-numbering is only recovered
! when accessing pdbnm or comparing to unkbnd
!
subroutine infer_polymer(i,nl,unkbnd,pdbnm)
!
  use pdb
  use zmatrix
  use sequen
  use molecule
  use math
  use params
  use atoms
  use polypep
  use iounit
  use system
  use fyoc
!
  implicit none
!
  integer, INTENT(IN) :: i,nl,unkbnd(n_pdbunk,4)
  character(4), INTENT(IN) :: pdbnm(nl)
  integer j,ji,k,ki,rs,imol,nnuci,shf,lastnuc,hnuc,vlnc
  logical bbone,tobbone
!
  rs = unkbnd(i,3)
  imol = molofrs(rs)
!
! single-molecule residue
  if ((rs.eq.rsmol(imol,1)).AND.(rs.eq.rsmol(imol,2))) then
!   for now, we will set no pointers for these (essentially rigid small molecules)
    do j=1,i-1
      if ((pdb_unkbnd(i,1).eq.pdb_unkbnd(j,1)).AND.(pdb_unkbnd(i,2).eq.pdb_unkbnd(j,2))) return
    end do
!
! N-terminal
  else if (rs.eq.rsmol(imol,1)) then
    hnuc = 0
    lastnuc = 0
    do j=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
      ji = j-at(rs)%bb(1)+unkbnd(i,1)
      if ((pdbnm(ji).eq.' O3*').OR.(pdbnm(ji).eq." O3'")) then
        nuci(rs,1) = j
      end if
      if ((nuci(rs,1).gt.0).AND.(hnuc.le.0)) then
        if (pdbnm(ji).eq.' HOP') then
          if (iz(1,j).eq.nuci(rs,1)) then
            hnuc = j
          end if
        end if
      end if
      if (nuci(rs,1).gt.0) then
        if (pdbnm(ji).eq.' P  ') then
          if (iz(1,j).eq.nuci(rs,1)) then
            nuci(rs,2) = j
          end if
        end if
      end if
      if ((pdbnm(ji).eq.' O5*').OR.(pdbnm(ji).eq." O5'")) then
        if (nuci(rs,2).gt.0) then
          if (iz(1,j).eq.nuci(rs,2)) then
            nuci(rs,3) = j
            lastnuc = 3
          end if
        else if (nuci(rs,1).le.0) then
          nuci(rs,1) = j
          lastnuc = 1
        end if
      end if
      if ((nuci(rs,1).gt.0).AND.(hnuc.le.0)) then
        if (pdbnm(ji).eq.'5HO*') then
          if (iz(1,j).eq.nuci(rs,1)) then
            hnuc = j
          end if
        end if
      end if
      if ((lastnuc.gt.0).AND.(lastnuc.le.5)) then
        if (nuci(rs,lastnuc).gt.0) then
          if ((pdbnm(ji).eq.' C5*').OR.(pdbnm(ji).eq." C5'")) then
            if (iz(1,j).eq.nuci(rs,lastnuc)) then
              lastnuc = lastnuc + 1
              nuci(rs,lastnuc) = j
            end if
          end if
        end if
        if ((nuci(rs,lastnuc).gt.0).AND.(lastnuc.lt.6)) then
         if ((pdbnm(ji).eq.' C4*').OR.(pdbnm(ji).eq." C4'")) then
            if (iz(1,j).eq.nuci(rs,lastnuc)) then
              lastnuc = lastnuc + 1
              nuci(rs,lastnuc) = j
            end if
          end if
        end if
        if  ((nuci(rs,lastnuc).gt.0).AND.(lastnuc.lt.6)) then
         if ((pdbnm(ji).eq.' C3*').OR.(pdbnm(ji).eq." C3'")) then
            if (iz(1,j).eq.nuci(rs,lastnuc)) then
              lastnuc = lastnuc + 1
              nuci(rs,lastnuc) = j
            end if
          end if
        end if
      end if
    end do
!   note that the last nuc-angle (C4-C3-O3(+1)-P(+1)) is set elsewhere (correct_topology(...))
    if ((lastnuc.eq.6).OR.(lastnuc.eq.4)) then
      seqpolty(rs) = 'N'
      if ((hnuc.gt.0).AND.(nuci(rs,1).gt.0).AND.(nuci(rs,2).gt.0).AND.(nuci(rs,3).gt.0)) then
        j = hnuc
        if ((iz(1,j).eq.nuci(rs,1)).AND.(iz(2,j).eq.nuci(rs,2)).AND.(iz(3,j).eq.nuci(rs,3))) then
          nnucs(rs) = nnucs(rs)+1
          nucsline(nnucs(rs),rs) = j
        end if
      end if
      if (lastnuc.eq.6) then
        if (nuci(rs,4).gt.0) then
          j = nuci(rs,4)
          if ((iz(1,j).eq.nuci(rs,3)).AND.(iz(2,j).eq.nuci(rs,2)).AND.(iz(3,j).eq.nuci(rs,1))) then
            nnucs(rs) = nnucs(rs)+1
            nucsline(nnucs(rs),rs) = j
          end if
        end if
        if (nuci(rs,5).gt.0) then
          j = nuci(rs,5)
          if ((iz(1,j).eq.nuci(rs,4)).AND.(iz(2,j).eq.nuci(rs,3)).AND.(iz(3,j).eq.nuci(rs,2))) then
            nnucs(rs) = nnucs(rs)+1
            nucsline(nnucs(rs),rs) = j
          end if
        end if
      end if
      if (nuci(rs,lastnuc).gt.0) then
        j = nuci(rs,lastnuc)
        if ((iz(1,j).eq.nuci(rs,lastnuc-1)).AND.(iz(2,j).eq.nuci(rs,lastnuc-2)).AND.(iz(3,j).eq.nuci(rs,lastnuc-3))) then
          nnucs(rs) = nnucs(rs)+1
          nucsline(nnucs(rs),rs) = j
        end if
      end if
    end if
    hnuc = 0
    do j=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
      ji = j-at(rs)%bb(1)+unkbnd(i,1)
      if (pdbnm(ji).eq.' N  ') then
        ni(rs) = j
      end if
      if (ni(rs).gt.0) then    
        if (pdbnm(ji).eq.'1H  ') then
          if (iz(1,j).eq.ni(rs)) then
            hnuc = j
          end if
        end if
        if (pdbnm(ji).eq.' CA ') then
          if (iz(1,j).eq.ni(rs)) then
            cai(rs) = j
          end if
        end if
      else
        if ((pdbnm(ji).eq.' CA ').OR.(pdbnm(ji).eq.' CH3')) then
          cai(rs) = j
        end if
      end if
      if (cai(rs).gt.0) then
        if (pdbnm(ji).eq.' C  ') then
          if (iz(1,j).eq.cai(rs)) then
            ci(rs) = j
          end if
        end if
        if ((pdbnm(ji)(2:4).eq.'HA ').AND.(ua_model.le.0)) then
          if (iz(1,j).eq.cai(rs)) then
            if (at(rs)%nsc.gt.0) then
              at(rs)%sc(2:(at(rs)%nsc+1)) = at(rs)%sc(1:at(rs)%nsc)
            end if
            at(rs)%sc(1) = j
            at(rs)%nsc = at(rs)%nsc + 1
            if (j-at(rs)%bb(1)+1.lt.at(rs)%nbb) then
              at(rs)%bb((j-at(rs)%bb(1)+1):(at(rs)%nbb-1)) = at(rs)%bb((j-at(rs)%bb(1)+2):at(rs)%nbb)
            end if
            at(rs)%nbb = at(rs)%nbb - 1
          end if
        end if
      else if (rs.lt.rsmol(imol,2)) then
        if (ni(rs+1).gt.0) then
          if (pdbnm(ji).eq.' C  ') then
            if (iz(1,ni(rs+1)).eq.j) then
              ci(rs) = j
            end if
          end if
        else if ((seqtyp(rs+1).eq.-1).AND.(ci(rs).le.0)) then
          if (pdbnm(ji).eq.' C  ') then
            ci(rs) = j
          end if
        end if
      end if
      if (ci(rs).gt.0) then
        if ((pdbnm(ji).eq.' O  ').OR.(pdbnm(ji).eq.'1OXT')) then
          if (iz(1,j).eq.ci(rs)) then
            oi(rs) = j
          end if
        end if
      else if (rs.lt.rsmol(imol,2)) then
        if (ni(rs+1).gt.0) then
          if (pdbnm(ji).eq.' O  ') then
            if (iz(1,j).eq.iz(1,ni(rs+1))) then
              oi(rs) = j
            end if
          end if
        end if
      end if
    end do
    if ((hnuc.gt.0).AND.(hnuc.gt.(at(rs)%bb(1)+2))) then
      j = hnuc
      if ((iz(1,j).eq.ni(rs)).AND.(iz(2,j).eq.cai(rs)).AND.(iz(3,j).eq.ci(rs))) then
        fline(rs) = j
      end if
    end if
    if ((oi(rs).gt.0).AND.(oi(rs).gt.unkbnd(i,1)+2)) then
      j = oi(rs)
      if ((iz(1,j).eq.ci(rs)).AND.(iz(2,j).eq.cai(rs)).AND.(iz(3,j).eq.ni(rs))) then
        yline(rs) = j
      end if
    end if
    if ((ni(rs).gt.0).AND.(cai(rs).gt.0).AND.(ci(rs).gt.0).AND.(oi(rs).gt.0)) seqpolty(rs) = 'P'
!
! middle chain or C-terminal
  else 
    do j=1,6
      if (nuci(rs-1,j).le.0) exit
    end do
    nnuci = j-1
    if (nnuci.eq.6) then
      shf = 0
    else if (nnuci.eq.4) then
      shf = 2
    end if
    if (nnuci.gt.0) then
      do j=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
        ji = j-at(rs)%bb(1)+unkbnd(i,1)
        if ((pdbnm(ji).eq.' O3*').OR.(pdbnm(ji).eq." O3'")) then
          if (iz(1,j).eq.nuci(rs-1,nnuci)) then
            nuci(rs,1) = j
          end if
        end if
        if (nuci(rs,1).gt.0) then
          if (pdbnm(ji).eq.' P  ') then
            if (iz(1,j).eq.nuci(rs,1)) then
              nuci(rs,2) = j
            end if
          end if
        end if
        if (nuci(rs,2).gt.0) then
          if ((pdbnm(ji).eq.' O5*').OR.(pdbnm(ji).eq." O5'")) then
            if (iz(1,j).eq.nuci(rs,2)) then
              nuci(rs,3) = j
            end if
          end if
        end if
        if (nuci(rs,3).gt.0) then
          if ((pdbnm(ji).eq.' C5*').OR.(pdbnm(ji).eq." C5'")) then
            if (iz(1,j).eq.nuci(rs,3)) then
              nuci(rs,4) = j
            end if
          end if
        end if
        if (nuci(rs,4).gt.0) then
          if ((pdbnm(ji).eq.' C4*').OR.(pdbnm(ji).eq." C4'")) then
            if (iz(1,j).eq.nuci(rs,4)) then
              nuci(rs,5) = j
            end if
          end if
        end if
        if (nuci(rs,5).gt.0) then
          if ((pdbnm(ji).eq.' C3*').OR.(pdbnm(ji).eq." C3'")) then
            if (iz(1,j).eq.nuci(rs,5)) then
              nuci(rs,6) = j
            end if
          end if
        end if
!       note that the last nuc-angle (C4-C3-O3(+1)-P(+1)) is set elsewhere (correct_topology(...))
        if ((nnuci.eq.6).OR.(nnuci.eq.4)) then
          if ((nuci(rs,3).gt.0).AND.(j.eq.nuci(rs,3))) then
            if ((iz(1,j).eq.nuci(rs,2)).AND.(iz(2,j).eq.nuci(rs,1)).AND.(iz(3,j).eq.nuci(rs-1,6-shf))) then
              nnucs(rs) = nnucs(rs)+1
              nucsline(nnucs(rs),rs) = j
            end if
          end if
        end if
        if ((nuci(rs,4).gt.0).AND.(j.eq.nuci(rs,4))) then
          if ((iz(1,j).eq.nuci(rs,3)).AND.(iz(2,j).eq.nuci(rs,2)).AND.(iz(3,j).eq.nuci(rs,1))) then
            nnucs(rs) = nnucs(rs)+1
            nucsline(nnucs(rs),rs) = j
          end if
        end if
        if ((nuci(rs,5).gt.0).AND.(j.eq.nuci(rs,5))) then
          if ((iz(1,j).eq.nuci(rs,4)).AND.(iz(2,j).eq.nuci(rs,3)).AND.(iz(3,j).eq.nuci(rs,2))) then
            nnucs(rs) = nnucs(rs)+1
            nucsline(nnucs(rs),rs) = j
          end if
        end if
        if ((nuci(rs,6).gt.0).AND.(j.eq.nuci(rs,6))) then
          if ((iz(1,j).eq.nuci(rs,5)).AND.(iz(2,j).eq.nuci(rs,4)).AND.(iz(3,j).eq.nuci(rs,3))) then
            nnucs(rs) = nnucs(rs)+1
            nucsline(nnucs(rs),rs) = j
          end if
        end if
      end do
    end if
!
    do j=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
      ji = j-at(rs)%bb(1)+unkbnd(i,1)
      if (ci(rs-1).gt.0) then
        if (pdbnm(ji).eq.' N  ') then
          if (iz(1,j).eq.ci(rs-1)) then
            ni(rs) = j
          end if
        end if
      end if
      if (ni(rs).gt.0) then    
        if (pdbnm(ji).eq.' HN ') then
          if (iz(1,j).eq.ni(rs)) then
            hni(rs) = j
          end if
        end if
        if (pdbnm(ji).eq.' CA ') then
          if (iz(1,j).eq.ni(rs)) then
            cai(rs) = j
          end if
        end if
      end if
      if (cai(rs).gt.0) then
        if (pdbnm(ji).eq.' C  ') then
          if (iz(1,j).eq.cai(rs)) then
            ci(rs) = j
          end if
        end if
        if ((pdbnm(ji)(2:4).eq.'HA ').AND.(ua_model.le.0)) then
          if (iz(1,j).eq.cai(rs)) then
            if (at(rs)%nsc.gt.0) then
              at(rs)%sc(2:(at(rs)%nsc+1)) = at(rs)%sc(1:at(rs)%nsc)
            end if
            at(rs)%sc(1) = j
            at(rs)%nsc = at(rs)%nsc + 1
            if ((j-at(rs)%bb(1)+1).lt.at(rs)%nbb) then
              at(rs)%bb((j-at(rs)%bb(1)+1):(at(rs)%nbb-1)) = at(rs)%bb((j-at(rs)%bb(1)+2):at(rs)%nbb)
            end if
            at(rs)%nbb = at(rs)%nbb - 1
          end if
        end if
      else if (rs.lt.rsmol(imol,2)) then
        if (ni(rs+1).gt.0) then
          if (pdbnm(ji).eq.' C  ') then
            if (iz(1,ni(rs+1)).eq.j) then
              ci(rs) = j
            end if
          end if
        end if
      end if
      if (ci(rs).gt.0) then
        if ((pdbnm(ji).eq.' O  ').OR.(pdbnm(ji).eq.'1OXT')) then
          if (iz(1,j).eq.ci(rs)) then
            oi(rs) = j
          end if
        end if
      else if (rs.lt.rsmol(imol,2)) then
        if (ni(rs+1).gt.0) then
          if (pdbnm(ji).eq.' O  ') then
            if (iz(1,j).eq.iz(1,ni(rs+1))) then
              oi(rs) = j
            end if
          end if
        end if
      end if
      if ((ci(rs).gt.0).AND.(yline(rs).gt.0).AND.(rs.eq.rsmol(imol,2))) then
        if ((pdbnm(ji).eq.'2OXT')) then
          if ((iz(1,j).eq.ci(rs)).AND.(iz(2,j).eq.cai(rs)).AND.(iz(3,j).eq.ni(rs))) then
            yline2(rs) = j
          end if
        end if
      end if
      if ((cai(rs).gt.0).AND.(j.eq.cai(rs))) then
        if ((iz(1,j).eq.ni(rs)).AND.(iz(2,j).eq.ci(rs-1))) then
          wline(rs) = j
        end if
      end if
      if ((ci(rs).gt.0).AND.(j.eq.ci(rs))) then
        if ((iz(1,j).eq.cai(rs)).AND.(iz(2,j).eq.ni(rs)).AND.(iz(3,j).eq.ci(rs-1))) then
          fline(rs) = j
        end if
      end if
      if ((hni(rs).gt.0).AND.(j.eq.hni(rs))) then
        if ((iz(1,j).eq.ni(rs)).AND.(iz(2,j).eq.cai(rs)).AND.(iz(3,j).eq.ci(rs))) then
          fline2(rs) = j
        end if
      end if
      if ((oi(rs).gt.0).AND.(j.eq.oi(rs))) then
        if ((iz(1,j).eq.ci(rs)).AND.(iz(2,j).eq.cai(rs)).AND.(iz(3,j).eq.ni(rs))) then
          yline(rs) = j
        end if
      end if
      if ((ni(rs).gt.0).AND.(j.eq.ni(rs)).AND.(yline(rs-1).gt.0)) then
        if (seqtyp(rs-1).eq.27) then
          if ((iz(1,j).eq.ci(rs-1)).AND.(iz(2,j).eq.cai(rs-1)).AND.(iz(3,j).eq.oi(rs-1))) then
            yline2(rs-1) = j
          end if
        else if ((seqtyp(rs-1).eq.28).AND.(ua_model.le.0)) then
          if ((iz(1,j).eq.ci(rs-1)).AND.(iz(2,j).eq.at(rs-1)%bb(3)).AND.(iz(3,j).eq.oi(rs-1))) then
            yline2(rs-1) = j
          end if
        else if ((seqtyp(rs-1).eq.28).AND.(ua_model.gt.0)) then
          if ((iz(1,j).eq.ci(rs-1)).AND.(iz(2,j).eq.oi(rs-1))) then
            yline2(rs-1) = j
          end if
        else
          if ((iz(1,j).eq.ci(rs-1)).AND.(iz(2,j).eq.cai(rs-1)).AND.(iz(3,j).eq.ni(rs-1))) then
            yline2(rs-1) = j
          end if
        end if
      end if
    end do
    if ((ni(rs).gt.0).AND.(cai(rs).gt.0).AND.(ci(rs).gt.0).AND.(oi(rs).gt.0)) seqpolty(rs) = 'P'
    if ((nuci(rs,6).gt.0).OR.((nuci(rs,4).gt.0).AND.(nuci(rs,5).le.0).AND.(rs.eq.rsmol(molofrs(rs),1)))) seqpolty(rs) = 'N'
    if ((yline2(rs).gt.0).AND.(yline(rs).le.0)) yline2(rs) = -1000
    if ((seqpolty(rs).ne.'P').AND.(yline(rs).gt.0)) then
      yline(rs) = -1000
      yline2(rs) = -1000
    end if
    if ((seqpolty(rs).ne.'P').AND.(fline(rs).gt.0)) then
      fline(rs) = -1000
      fline2(rs) = -1000
    end if
  end if
!
  if (seqpolty(rs).eq.'N') then
    do j=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
      ji = j-at(rs)%bb(1)+unkbnd(i,1)
      bbone = .false.
      tobbone = .false.
      if (pdbnm(ji).eq.'2O3*') cycle
      if (iz(1,j).gt.0) then
        if ((iz(1,j)-at(rs)%bb(1)+unkbnd(i,1).gt.0).AND.(iz(1,j)-at(rs)%bb(1)+unkbnd(i,1).le.nl)) then
          if (pdbnm(iz(1,j)-at(rs)%bb(1)+unkbnd(i,1)).eq.'2O3*') tobbone = .true.
        end if
      end if
      do k=1,6
        if (j.eq.nuci(rs,k)) then
          bbone = .true.
          exit
        end if
        if (iz(1,j).eq.nuci(rs,k)) then
          tobbone = .true.
        end if
      end do  
      if (bbone.EQV..true.) cycle
      vlnc = 0
      do k=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
        ki = k-at(rs)%bb(1)+unkbnd(i,1)
        if (j.eq.k) cycle
        if ((iz(1,j).eq.k).OR.(iz(1,k).eq.j)) then
          vlnc = vlnc + 1
        end if
      end do
      do k=1,nadd
        if ((j.eq.iadd(1,k)).OR.(j.eq.iadd(2,k))) then
          vlnc = vlnc + 1
        end if
      end do
      if ((vlnc.le.1).AND.(tobbone.EQV..true.)) cycle
      at(rs)%sc(at(rs)%nsc+1) = j
      at(rs)%nsc = at(rs)%nsc + 1
      if (at(rs)%nbb.gt.1) then
        do k=1,at(rs)%nbb
          if (at(rs)%bb(k).eq.j) then
            at(rs)%bb(k:(at(rs)%nbb-1)) = at(rs)%bb((k+1):at(rs)%nbb)
            exit
          end if
        end do
      end if
      at(rs)%nbb = at(rs)%nbb - 1
    end do
  else if (seqpolty(rs).eq.'P') then
    do j=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
      ji = j-at(rs)%bb(1)+unkbnd(i,1)
      tobbone = .false.
      if (j.eq.ci(rs)) cycle
      if (j.eq.ni(rs)) cycle
      if (j.eq.cai(rs)) cycle
      if (j.eq.oi(rs)) cycle
      if (j.eq.hni(rs)) cycle
      if (iz(1,j).eq.ci(rs)) tobbone = .true.
      if (iz(1,j).eq.cai(rs)) tobbone = .true.
      if (iz(1,j).eq.ni(rs)) tobbone = .true.
      vlnc = 0
      do k=at(rs)%bb(1),at(rs)%bb(1)+unkbnd(i,2)-unkbnd(i,1)
        ki = k-at(rs)%bb(1)+unkbnd(i,1)
        if (j.eq.k) cycle
        if ((iz(1,j).eq.k).OR.(iz(1,k).eq.j)) then
          vlnc = vlnc + 1
        end if
      end do
      do k=1,nadd
        if ((j.eq.iadd(1,k)).OR.(j.eq.iadd(2,k))) then
          vlnc = vlnc + 1
        end if
      end do
      if ((vlnc.le.1).AND.(tobbone.EQV..true.)) cycle
!     there is a potential problem here if the geometry is not good and HA gets linked up
!     with another atom -> have to skip eliminating it from bb-list as it is already in sc-list
      if (j.eq.at(rs)%sc(1)) cycle
!     for the remainder -> eliminate from backbone list (must be there)
      at(rs)%sc(at(rs)%nsc+1) = j
      at(rs)%nsc = at(rs)%nsc + 1
      if (at(rs)%nbb.gt.1) then
        do k=1,at(rs)%nbb-1
          if (at(rs)%bb(k).eq.j) then
            at(rs)%bb(k:(at(rs)%nbb-1)) = at(rs)%bb((k+1):at(rs)%nbb)
            exit
          end if
        end do
      end if
      at(rs)%nbb = at(rs)%nbb - 1
    end do
  end if
!  write(*,*) 'set for ',rs,': ',at(rs)%nbb,at(rs)%bb(1:at(rs)%nbb),at(rs)%nsc,at(rs)%sc(1:at(rs)%nsc)
!  write(*,*) fline(rs),fline2(rs),yline(rs),yline2(max(rs-1,1)),yline2(rs),wline(rs),ni(rs),cai(rs),&
! &ci(rs),hni(rs),oi(rs),at(rs)%sc(1)
!  write(*,*) nuci(rs,1:6),nucsline(1:nnucs(rs),rs),nucsline(1:nnucs(max(rs-1,1)),max(rs-1,1)),seqpolty(rs)
  if (seqpolty(rs).eq.'P') then
    write(ilog,*) 'Identified unsupported residue ',pdb_unknowns(i),' (#',rs,') as a polypeptide residue.'
  else if (seqpolty(rs).eq.'N') then
    write(ilog,*) 'Identified unsupported residue ',pdb_unknowns(i),' (#',rs,') as a polynucleotide residue.'
  else
    write(ilog,*) 'Unsupported residue ',pdb_unknowns(i),' (#',rs,') does not conform to a supported polymer type.'
  end if
!
end
!
!------------------------------------------------------------------------------------------------
!

