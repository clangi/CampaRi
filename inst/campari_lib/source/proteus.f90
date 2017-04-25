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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
! the proteus-routines do the bulk of setting up the molecular geometries
! the pointers for the active degrees of freedom are set as well
! sidechain is the corresponding routine for sidechains (see below)
!
!--------------------------------------------------------------------
!
subroutine proteus_init(mode)
!
  use fyoc
  use atoms
  use movesets
  use sequen
!
  implicit none
!
  integer mode
!
  if (mode.eq.1) then
!   the currently active (i.e., next) atom
    n = 1
!   initialize move set lists
    fylst%nr = 0
    wlst%nr = 0
    chilst%nr = 0
    nuclst%nr = 0
    puclst%nr = 0
    nucpuclst%nr = 0
  else if (mode.eq.2) then
!   finalize atom array
    n = n - 1
!   transfer move set lists, initialize weights as flat
    if (fylst%nr.gt.0) then
      allocate(fylst%idx(fylst%nr))
      allocate(fylst%wt(fylst%nr))
      fylst%idx(1:fylst%nr) = resfy(1:fylst%nr)
      fylst%wt(1:fylst%nr) = 1.0*resfy(1:fylst%nr)
    end if
    deallocate(resfy)
    if (wlst%nr.gt.0) then
      allocate(wlst%idx(wlst%nr))
      allocate(wlst%wt(wlst%nr))
      wlst%idx(1:wlst%nr) = resw(1:wlst%nr)
      wlst%wt(1:wlst%nr) = 1.0*resw(1:wlst%nr)
    end if
    deallocate(resw)
    if (chilst%nr.gt.0) then
      allocate(chilst%idx(chilst%nr))
      allocate(chilst%wt(chilst%nr))
      chilst%idx(1:chilst%nr) = reschi(1:chilst%nr)
      chilst%wt(1:chilst%nr) = 1.0*reschi(1:chilst%nr)
    end if
    deallocate(reschi)
    if (nuclst%nr.gt.0) then
      allocate(nuclst%idx(nuclst%nr))
      allocate(nuclst%wt(nuclst%nr))
      nuclst%idx(1:nuclst%nr) = resnuc(1:nuclst%nr)
      nuclst%wt(1:nuclst%nr) = 1.0*resnuc(1:nuclst%nr)
    end if
    deallocate(resnuc)
    if (puclst%nr.gt.0) then
      allocate(puclst%idx(puclst%nr))
      allocate(puclst%wt(puclst%nr))
      puclst%idx(1:puclst%nr) = respuc(1:puclst%nr)
      puclst%wt(1:puclst%nr) = 1.0*respuc(1:puclst%nr)
    end if
    deallocate(respuc)
    if (nucpuclst%nr.gt.0) then
      allocate(nucpuclst%idx(nucpuclst%nr))
      allocate(nucpuclst%wt(nucpuclst%nr))
      nucpuclst%idx(1:nucpuclst%nr) = resnucpuc(1:nucpuclst%nr)
      nucpuclst%wt(1:nucpuclst%nr) = 1.0*resnucpuc(1:nucpuclst%nr)
    end if
    deallocate(resnucpuc)
  end if 
!     
end
!
!--------------------------------------------------------------------
!
! torsional setup
! for polypeptides the different angles are defined as follows
!
! all for residue i:
! omega (w): CA(i-1)->C(i-1)->N(i)->CA(i) (wline)
! phi (f): C(i-1)->N(i)->CA(i)->C(i)      (fline)
! dependent: C(i)->CA(i)->N(i)->HN(i)     (fline2)
! psi (y): N(i)->CA(i)->C(i)->O(i)        (yline)
! dependent: N(i)->CA(i)->C(i)->N(i+1)    (yline2)
!     
! note that this means that yline2 is only encountered
! in residue i+1
! all of these are potentially ill-defined if reference atoms
! are missing, although an attempt is made to overcome this (whenever possible)
! by pointing N(i),CA(i),... to appropriate OTHER atoms
! example: the methyl carbon in ACE is treated as CA
!
subroutine proteus(imol)
!
  use aminos
  use atoms
  use fyoc
  use iounit
  use polypep
  use sequen
  use molecule
  use system
  use movesets
!
  implicit none
!
  integer i,k,imol,nn1,nn2,nn3,nc1,nc21,nc22,c1
  integer ntyp(MAXAMINO),catyp(MAXAMINO),ctyp(MAXAMINO),nc2,v,shf
  integer hntyp(MAXAMINO),otyp(MAXAMINO),hatyp(MAXAMINO),j
  integer ptyp(12),o5typ(12),o3typ(12),optyp(12),h3typ(12)
  integer c5typ(12),c4typ(12),c3typ(12),h5typ(12),h4typ(12)
  integer th5typ(12),th3typ(12),to3typ(12),tc3typ(12)
  integer otze,c3m1,c4m1,c5m1,km1,km2,nnd1,ncd2,ncd1
  character(3) resname
  logical cyclic
  RTYPE dh(200),an(200),bo(200),dihi
!
! biopolymer atom types for amino acid backbone atoms
!
  data ntyp   /   1,   7,  15,  25,  37,  51,  61,  73,  92,&
 &              103, 118, 134, 157, 174, 190, 206, 216, 228,&
 &              240, 254, 267, 283, 300, 314, 321,   1,   0,&
 &                0,   0,   0, 385, 395, 408, 418, 430,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0, 519,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,1074,&
 &             1088,1116,1125,1100,1140,1159,1172,1187,1207,&
 &             1225,1243/
  data catyp  /   2,   8,  16,  26,  38,  52,  62,  74,  93,&
 &              104, 119, 135, 158, 175, 191, 207, 217, 229,&
 &              241, 255, 268, 284, 301, 315, 322,   2,   0,&
 &                0,   0,   0, 386, 396, 409, 419, 431,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0, 520,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,1075,&
 &             1089,1117,1126,1101,1141,1160,1173,1188,1208,&
 &             1226,1244/
  data ctyp   /   3,   9,  17,  27,  39,  53,  63,  75,  94,&
 &              105, 120, 136, 159, 176, 192, 208, 218, 230,&
 &              242, 256, 269, 285, 302, 316, 323,   3,   0,&
 &                0,   0,   0, 387, 397, 410, 420, 432,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0, 521,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,1076,&
 &             1090,1118,1127,1102,1142,1161,1174,1189,1209,&
 &             1227,1245/
  data hntyp  /   4,  10,  18,  28,  40,  54,  64,  76,   0,&
 &              106, 121, 137, 160, 177, 193, 209, 219, 231,&
 &              243, 257, 270, 286, 303, 317, 324,   4,   0,&
 &                0,   0,   0, 388,   0, 411, 421, 433,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0, 522,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,1077,&
 &             1091,1119,1128,1103,1143,1162,1175,1190,1210,&
 &             1228,1246/
  data otyp   /   5,  11,  19,  29,  41,  55,  65,  77,  95,&
 &              107, 122, 138, 161, 178, 194, 210, 220, 232,&
 &              244, 258, 271, 287, 304, 318, 325,   5,   0,&
 &                0,   0,   0, 389, 398, 412, 422, 434,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0, 523,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,1078,&
 &             1092,1120,1129,1104,1144,1163,1176,1191,1211,&
 &             1229,1247/
  data hatyp  /   6,  12,  20,  30,  42,  56,  66,  78,  96,&
 &              108, 123, 139, 162, 179, 195, 211, 221, 233,&
 &              245, 259, 272, 288, 305,   0, 326,   6,   0,&
 &                0,   0,   0, 390, 399, 413, 423, 435,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0, 524,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,   0,&
 &                0,   0,   0,   0,   0,   0,   0,   0,1079,&
 &             1093,1121,1130,1105,1145,1164,1177,1192,1212,&
 &             1230,1248/
!
! nucleic acid type arrays
  data ptyp   /724,671,778,834,892,698,645,751,806,863,624,936/
  data o5typ  /725,672,779,835,893,699,646,752,807,864,625,937/
  data o3typ  /726,673,780,836,894,700,647,753,808,865,626,938/
  data optyp  /727,674,781,837,895,701,648,754,809,866,627,939/
  data c5typ  /728,675,782,838,896,702,649,755,810,867,628,940/
  data c4typ  /730,677,784,840,898,704,651,757,812,869,630,942/
  data c3typ  /733,680,787,843,901,707,654,760,815,872,633,945/
  data h5typ  /729,676,783,839,897,703,650,756,811,868,629,941/
  data h4typ  /732,679,786,842,900,706,653,759,814,871,632,944/
  data h3typ  /734,681,788,844,902,708,655,761,816,873,635,947/
  data tc3typ /925,925,925,925,925,928,928,928,928,928,925,928/
  data to3typ /926,926,926,926,926,929,929,929,929,929,634,946/
  data th3typ /927,927,927,927,927,930,930,930,930,930,636,948/
  data th5typ /924,924,924,924,924,924,924,924,924,924,956,955/
!
!
! set the more common values - i realize this is terrible, but it is simply bad coding
! to put the values directly in the argument list of zatom()
! -> hence we need to pin them to a local variable somewhere, but that would be terribly
! long if done for each possible atom individually
! so we use these non-descript arrays to fill up with common values -> good for find+replace strategies
!
! bond lengths
  bo(1) = 0.0d0
  bo(2) = 1.09d0
  bo(3) = 1.01d0
  bo(4) = 1.522d0
  bo(5) = 1.529d0
  bo(6) = 0.9572d0   ! TIP3/4/5P O-H
  bo(7) = 1.329d0
  bo(8) = 1.01d0
  bo(9) = 1.525d0
  bo(10) = 1.25d0
  bo(11) = 1.02d0
  bo(12) = 1.08d0
  bo(13) = 1.229d0
  bo(14) = 1.335d0
  bo(15) = 0.7d0     ! TIP5P O-LP
  bo(16) = 1.458d0
  bo(17) = 1.34d0
  bo(18) = 1.4660d0
  bo(19) = 1.231d0   ! Engh/Huber C-O  
  bo(20) = 1.449d0
  bo(21) = 0.15d0    ! TIP4P O-M
  bo(22) = 1.0d0     ! SPC O-H
  bo(23) = 1.092d0   ! N-terminal proline N-H
  bo(24) = 1.4d0     ! OPLS-AA PCR arC-arC
  bo(25) = 1.41d0    ! OPLS-AA MOH C-O
  bo(26) = 0.945d0   ! OPLS-AA MOH/PCR O-H
  bo(27) = 1.51d0    ! CA-C in N/C-terminal standard residues + C-terminal PRO, also CT-C1 in OPLSAA-PCR
  bo(28) = 1.11d0    ! CA-HA in C-terminal standard residues + C-terminal PRO
  bo(29) = 1.50d0    ! N-CA in N-terminal AIB
  bo(30) = 1.364d0   ! C6-O in OPLSAA-PCR
! now the nucleic acid parameters (Parkinson et al. in Acta Cryst. D52 (1996) p57-64)
  bo(31) = 1.397d0   ! N1-C2 in CYT
  bo(32) = 1.240d0   ! C2-O2 in CYT
  bo(33) = 1.367d0   ! N1-C6 in CYT
  bo(34) = 1.353d0   ! C2-N3 in CYT
  bo(35) = 1.339d0   ! C5-C6 in CYT, C5-C6 in THY, C2-N1 in ADE
  bo(36) = 1.335d0   ! N3-C4 in CYT, C4-N4 in CYT, C6-N6 in ADE
  bo(37) = 1.381d0   ! N1-C2 in URA
  bo(38) = 1.219d0   ! C2-O2 in URA
  bo(39) = 1.375d0   ! N1-C6 in URA, N9-C4 in GUA
  bo(40) = 1.373d0   ! C2-N3 in URA, C2-N3 in THY, N9-C8 in ADE, C2-N1 in GUA
  bo(41) = 1.337d0   ! C5-C6 in URA
  bo(42) = 1.380d0   ! N3-C4 in URA and some TRP
  bo(43) = 1.232d0   ! C4-O4 in URA
  bo(44) = 1.376d0   ! N1-C2 in THY
  bo(45) = 1.220d0   ! C2-O2 in THY
  bo(46) = 1.378d0   ! N1-C6 in THY
  bo(47) = 1.382d0   ! N3-C4 in THY
  bo(48) = 1.228d0   ! C4-O4 in THY
  bo(49) = 1.496d0   ! C5-CT in THY
  bo(50) = 1.374d0   ! N9-C4 in ADE, N9-C8 in GUA
  bo(51) = 1.383d0   ! C4-C5 in ADE
  bo(52) = 1.311d0   ! C8-N7 in ADE
  bo(53) = 1.344d0   ! C4-N3 in ADE
  bo(54) = 1.331d0   ! N3-C2 in ADE
  bo(55) = 1.406d0   ! C5-C6 in ADE
  bo(56) = 1.379d0   ! C4-C5 in GUA
  bo(57) = 1.305d0   ! C8-N7 in GUA
  bo(58) = 1.350d0   ! C4-N3 in GUA and Duffy/van Gunsteren/Smith UREA C-N
  bo(59) = 1.323d0   ! N3-C2 in GUA
  bo(60) = 1.341d0   ! C2-N2 in GUA
  bo(61) = 1.237d0   ! C6-O6 in GUA
  bo(62) = 1.419d0   ! C5-C6 in GUA
  bo(63) = 1.607d0   ! O3P-P
  bo(64) = 1.593d0   ! O5P-P
  bo(65) = 1.485d0   ! OP-P
  bo(66) = 1.440d0   ! O5P-C5* (R/D)
  bo(67) = 1.510d0   ! C5*-C4* (R)
  bo(68) = 1.511d0   ! C5*-C4* (D)
  bo(69) = 1.524d0   ! C4*-C3* (R)
  bo(70) = 1.528d0   ! C4*-C3* (D)
  bo(71) = 1.423d0   ! C3*-O3* (R)
  bo(72) = 1.431d0   ! C3*-O3* (D)
! some more small molecule related stuff (from OPLS typically)
  bo(73) = 1.504d0  ! methylimidazole CT-C1
  bo(74) = 1.371d0  ! methylimidazole C1-C2
  bo(76) = 1.335d0  ! methylimidazole C3-N(CC)
  bo(77) = 1.343d0  ! methylimidazole C3-N(HC)
  bo(78) = 1.08d0   ! methylimidazole C-H (aromatic)
  bo(79) = 1.394d0  ! methylimidazole HC(C-N)CC
  bo(80) = 1.381d0  ! methylimidazole XC(C-N)CH and some TRP
  bo(81) = 1.810d0  ! thiol or thioether CS
  bo(82) = 1.336d0  ! thiol SH
  bo(83) = 1.495d0  ! MIN-OPLS
  bo(84) = 1.352d0  ! MIN-OPLS
  bo(85) = 1.459d0  ! MIN-OPLS
  bo(86) = 1.3796d0 ! IMD-fit from OPLS-minimzation
  bo(87) = 1.3708d0 ! IMD-fit from OPLS-minimzation
  bo(88) = 1.3425d0 ! IMD-fit from OPLS-minimzation
  bo(89) = 1.3952d0 ! IMD-fit from OPLS-minimzation
  bo(90) = 1.3969d0 ! IME-fit from OPLS-minimzation
  bo(91) = 1.3723d0 ! IME-fit from OPLS-minimzation
  bo(92) = 1.335d0  ! IME-fit from OPLS-minimzation
  bo(93) = 1.3791d0 ! IME-fit from OPLS-minimzation
  bo(94) = 1.265d0  ! C-O in Duffy/Jorgensen/van Gunsteren URE
  bo(95) = 1.320d0  ! C-N in Weerasinghe/Smith GDN
  bo(96) = 1.414d0  ! Cl-O in perchlorate anion
  bo(97) = 1.24d0   ! N-O in nitrate
  bo(98) = 1.494d0  ! NH2C2 in sec. ammonium ions  
  bo(99) = 1.489d0  ! NH3C in prim. ammonium ions
  bo(100) = 0.125d0 ! TIP4P-Ew O-M
  bo(101) = 1.36d0  ! shortest BL in napthalene
  bo(102) = 1.42d0  ! longest BL in naphthalene
  bo(103) = 1.21d0 ! molecular oxygen
!
! bond angles 
  an(1) = 0.0d0
  an(2) = 120.0d0    ! CCC, CCH in PCR, CNH in NH2/URE, NCH in FOA/FOR, O(CA)CN in NH2/NME, CCO in ACE, O(CA)CN
!                    ! in C-terminal residues, generally around idealized, planar sp2-centers
  an(3) = 109.47d0   ! generic tetrahedral (usually alkyl) CCH, HCH
  an(4) = 122.9d0    ! OCN in all OPLS-AA amides
  an(5) = 119.8d0    ! CNH in all OPLS-AA amides
  an(6) = 52.26d0    ! HOM in T4P
  an(7) = 104.52d0   ! HOH in T3/4/5P
  an(8) = 110.69477d0  ! HOM in T5P (determined numerically)
  an(9) = 110.0d0    ! NCAHA in N-terminal residues except PRO and uncharged glycine
  an(10) = 117.0d0   ! CACOXT in C-terminal residues
  an(11) = 121.7d0   ! CNCA in all in-chain and all standard C-terminal residues (+GLY)
  an(12) = 116.2d0   ! CACN in all in-chain and all (except PRO) C-terminal residues incl. NME,NH2 and some MIN
  an(13) = 121.0d0   ! CNCA in C-terminal PRO/AIB, CANHN in C-terminal GLY/AIB
  an(14) = 119.15d0  ! CNHN in all in-chain and all standard C-terminal residues and in NME
  an(15) = 111.6d0   ! NCAC in C-terminal AIB/PRO and HID ND1-CG-CD2 and some MIN
  an(16) = 120.6d0   ! CNCH3 in NME
  an(17) = 110.5d0   ! NCAC in all N-terminal residues, all in-chain residues, and all C-terminal residues except PRO/AIB
  an(18) = 110.8d0   ! NCAHA in N-terminal PRO
  an(19) = 112.7d0   ! CACN in C-terminal PRO/AIB
  an(20) = 122.5d0   ! CACO in N-terminal AIB
  an(21) = 123.0d0   ! CACO in N-terminal PRO
  an(22) = 120.8d0   ! CACO in remaining N-terminal residues except ACE
  an(23) = 110.7d0   ! NCAC in N-terminal
  an(24) = 120.5d0   ! CACO in in-chain residues
  an(25) = 120.4d0   ! OCCCT in OPLS-AA amides
  an(26) = 120.45d0  ! OCCCT in OPLS-AA NMA (modified for better planarity)
  an(27) = 116.65d0  ! NCCCT in OPLS-AA NMA (see above)
  an(28) = 116.05d0  ! NCH in OPLS-AA NMF (see above)
  an(29) = 121.9d0   ! CNCNT in OPLS-AA amides
  an(30) = 111.10d0  ! CCCTCBT in OPLS-AA PPA
  an(31) = 108.5d0   ! COH in OPLS-AA MOH
  an(32) = 113.0d0   ! COH in OPLS-AA PCR, presumed C-N-C in dimethylammonium
  an(33) = 118.4d0   ! CNHN in OPLS-AA sec. amides
  an(34) = 117.2d0   ! NCN in van Gunsteren URE (should be equal to Duffy/Jorgensen)
  an(35) = 121.05d0  ! OCH in OPLS-AA NMF (see above)
  an(36) = 107.9d0   ! CCAHA in all residues except caps and N-terminal PRO
  an(37) = 108.0d0   ! HNH in N-terminal GLY,AIB
  an(38) = 111.0d0!106.6d0   ! CCAHA in N-terminal PRO
  an(39) = 120.25d0  ! CANHN in NME 
  an(40) = 126.0d0   ! OCO in OPLS-AA carobxylate
  an(41) = 121.8d0   ! corrected CNCNT in OPLS-AA amides
! now the nucleic acid parameters (Parkinson et al. in Acta Cryst. D52 (1996) p57-64)
  an(42) = 120.3d0   ! C6N1C2 in CYT
  an(43) = 118.9d0   ! N1C2O2 in CYT
  an(44) = 119.2d0   ! N1C2N3 in CYT
  an(45) = 121.0d0   ! N1C6C5 in CYT
  an(46) = 118.0d0   ! N3C4N4 in CYT
  an(47) = 119.85d0  ! for H1 in CYT
  an(48) = 0.0
  an(49) = 0.0
  an(50) = 0.0
  an(51) = 121.0d0   ! C6N1C2 in URA, N1C6C5 in CYT
  an(52) = 114.9d0   ! N1C2N3 in URA
  an(53) = 127.0d0   ! C2N3C4 in URA
  an(54) = 122.8d0   ! N1C2O2 in URA
  an(55) = 122.7d0   ! N1C6C5 in URA and some MIN
  an(56) = 119.4d0   ! N3C4O4 in URA
  an(57) = 119.5d0   ! for H1 in URA, for H6 in CYT
  an(58) = 118.65d0  ! for H6 in URA
  an(59) = 116.5d0   ! for H3 in URA
  an(60) = 120.15d0  ! for H5 in URA
  an(61) = 121.3d0   ! C6N1C2 in THY, for H5 in CYT
  an(62) = 114.6d0   ! N1C2N3 in THY
  an(63) = 127.2d0   ! C2N3C4 in THY
  an(64) = 123.1d0   ! N1C2O2 in THY, for H8 in ADE and some MIN
  an(65) = 123.7d0   ! N1C6C5 in THY
  an(66) = 119.9d0   ! N3C4O4 in THY, N2C2N3 in GUA, C2N3C4 in CYT
  an(67) = 119.0d0   ! C4C5CT in THY
  an(68) = 119.35d0  ! for H1 in THY
  an(69) = 118.15d0  ! for H6 in THY
  an(70) = 116.4d0   ! for H3 in THY
  an(71) = 118.6d0   ! N1C6N6 in ADE
  an(72) = 105.8d0   ! N9C4C5 in ADE, C8N9C4 in ADE
  an(73) = 113.8d0   ! N9C8N7 in ADE
  an(74) = 127.4d0   ! N9C4N3 in ADE
  an(75) = 110.6d0   ! C4N3C2 in ADE
  an(76) = 117.0d0   ! C4C5C6 in ADE
  an(77) = 129.3d0   ! N3C2N1 in ADE
  an(78) = 115.35d0  ! for H2 in ADE
  an(79) = 127.1d0   ! for H9 in ADE
  an(80) = 106.4d0   ! C8N9C4 in GUA
  an(81) = 128.6d0   ! C5C6O6 in GUA and some MIN
  an(82) = 105.4d0   ! N9C4C5 in GUA
  an(83) = 113.1d0   ! N9C8N7 in GUA
  an(84) = 126.0d0   ! N9C4N3 in GUA
  an(85) = 111.9d0   ! C4N3C2 in GUA
  an(86) = 118.8d0   ! C4C5C6 in GUA
  an(87) = 123.9d0   ! N3C2N1 in GUA
  an(88) = 126.8d0   ! for H9 in GUA
  an(89) = 123.45d0  ! for H8 in GUA
  an(90) = 117.45d0  ! for H1 in GUA
  an(91) = 104.0d0   ! O3P-P-O5P
  an(92) = 108.1d0   ! O5P-P-OP
  an(93) = 107.4d0   ! O3P-P-OP
  an(94) = 120.9d0   ! P-O5P-C5* (R/D)
  an(95) = 110.2d0   ! O5P-C5*-C4* (R/D)
  an(96) = 115.5d0   ! C5*-C4*-C3* (R)
  an(97) = 114.7d0   ! C5*-C4*-C3* (D)
  an(98) = 110.6d0   ! C4*-C3*-O3*(P) (R)
  an(99) = 110.3d0   ! C4*-C3*-O3*(P) (D)
  an(100) = 119.7d0  ! C3*-O3*(P)-P (R/D)
! some more small molecule stuff (mostly OPLS)
  an(101) = 121.4d0  ! CNO in van Gunsteren URE (should be equal to Duffy/Jorgensen)
  an(102) = 121.6d0  ! CTCN(HC) in methylimidazole and some MIN
  an(103) = 124.5d0  ! CTCN(CC) in methylimidazole
  an(104) = 112.7d0  ! CCC standard in OPLS
  an(105) = 110.7d0  ! HCC standard in OPLS
  an(106) = 107.8d0  ! HCH standard in OPLS
  an(107) = 96.0d0   ! CSH in thiols
  an(108) = 114.7d0  ! CCS in sulfides
  an(109) = 98.9d0   ! CSC in sulfides
  an(110) = 107.3d0  !  in methylimidazole
  an(111) = 108.3d0  !  in methylimidazole
  an(112) = 126.35d0 !  in methylimidazole
  an(113) = 130.7d0  !  in methylimidazole
  an(114) = 108.7d0  !  in methylimidazole and some MIN
  an(115) = 128.2d0  !  in methylimidazole
  an(116) = 124.2d0  !  in methylimidazole (adjusted)
  an(117) = 105.3d0  !  in methylimidazole
  an(118) = 125.0d0  ! some MIN
  an(119) = 121.69d0 ! IMD-fit from OPLS-minimzation
  an(120) = 130.79d0 ! IMD-fit from OPLS-minimzation
  an(121) = 106.23d0 ! IMD-fit from OPLS-minimzation
  an(122) = 108.52d0 ! IMD-fit from OPLS-minimzation
  an(123) = 123.09d0 ! IME-fit from OPLS-minimzation
  an(124) = 129.29d0 ! IME-fit from OPLS-minimzation
  an(125) = 105.68d0 ! IME-fit from OPLS-minimzation
  an(126) = 108.14d0 ! IME-fit from OPLS-minimzation
  an(127) = 136.00d0 ! HNH in Smith GDN
  an(128) = 112.00d0 ! CNH in Smith GDN
  an(129) = 118.45d0 ! OCH in FOA
!
!
  dh(1) = 0.0d0
  dh(2) = 109.47d0
  dh(3) = 180.0d0
  dh(4) = 147.3d0    ! C5*-C4*-C3*-O3*(P) (C2*-endo R)
  dh(5) = 145.2d0    ! C5*-C4*-C3*-O3*(P) (C2*-endo D)
  dh(6) = 26.62d0    ! C5*-C4*-C3*-H3*(P) (C2*-endo R)
  dh(7) = 26.62d0    ! C5*-C4*-C3*-H3*(P) (C2*-endo D)
!
  cyclic = .false.
  ! WARNING: make this (imol): !!!if (globcyclic.EQV..true.) cyclic = .true.
!
! actually, this is currently not supported
!
! number of backbone and side-chain torsion angles for peptides
!
! initialize arrays for current chain
  do i = rsmol(imol,1),rsmol(imol,2)
    nchi(i) = 0
    at(i)%nbb = 0
    at(i)%nsc = 0
!    tauline(i) = 0
    ni(i) = 0
    hni(i) = 0
    cai(i) = 0
    ci(i) = 0
    oi(i) = 0
    nnucs(i) = 0
    do j=1,6
      nuci(i,j) = 0
    end do
  end do
!
! set atom counter and the first residue number and type
  atmol(imol,1) = n
  i = rsmol(imol,1)
  if ((seqtyp(i).gt.0).AND.(seqtyp(i).le.MAXAMINO)) then
    resname = amino(seqtyp(i))
  else if (seqtyp(i).eq.-1) then
    resname = '???'
  else
    write(ilog,*) 'Fatal. Encountered unsupported residue type in proteus(...). This is an &
 &omission bug.'
    call fexit()
  end if
    
!
!
! take care of free amino acids and other one-residue molecules
  if ((rsmol(imol,2)-rsmol(imol,1)).eq.0) then
!
    if (resname.eq.'NA+') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(444,bo(1),an(1),dh(1),0,0,0,0)
!
    else if (resname.eq.'CL-') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(445,bo(1),an(1),dh(1),0,0,0,0)
!
    else if (resname.eq.'K+ ') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(531,bo(1),an(1),dh(1),0,0,0,0)
!
    else if (resname.eq.'BR-') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(532,bo(1),an(1),dh(1),0,0,0,0)
!
    else if (resname.eq.'CS+') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(533,bo(1),an(1),dh(1),0,0,0,0)
!
    else if (resname.eq.'I- ') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(534,bo(1),an(1),dh(1),0,0,0,0)
!
    else if (resname.eq.'O2 ') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1073,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1073,bo(103),an(1),dh(1),n-1,0,0,0)
!
    else if (resname.eq.'LCP') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1046,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1047,bo(96),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1047,bo(96),an(3),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1047,bo(96),an(3),an(3),n-3,n-2,n-1,-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1047,bo(96),an(3),an(3),n-4,n-3,n-2,1)
!
    else if (resname.eq.'NO3') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1048,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1049,bo(97),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1049,bo(97),an(2),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1049,bo(97),an(2),an(2),n-3,n-2,n-1,-1)
!
    else if (resname.eq.'GDN') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(541,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
!      call z_add(542,bo(17),an(1),dh(1),n-1,0,0,0)
      call z_add(542,bo(95),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn2 = n
!      call z_add(542,bo(17),an(2),dh(1),n-2,n-1,0,0)
      call z_add(542,bo(95),an(2),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn3 = n
!      call z_add(542,bo(17),an(2),dh(3),n-3,n-2,n-1,0)
      call z_add(542,bo(95),an(2),dh(3),n-3,n-2,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(543,bo(8),an(2),dh(1),n-3,n-4,n-2,0)
      call z_add(543,bo(22),an(128),dh(1),n-3,n-4,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(543,bo(8),an(2),dh(3),n-4,n-5,n-1,0)
      call z_add(543,bo(22),an(128),dh(3),n-4,n-5,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(543,bo(8),an(2),dh(1),n-4,n-6,n-3,0)
      call z_add(543,bo(22),an(128),dh(1),n-4,n-6,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(543,bo(8),an(2),dh(3),n-5,n-7,n-1,0)
      call z_add(543,bo(22),an(128),dh(3),n-5,n-7,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(543,bo(8),an(2),dh(1),n-5,n-8,n-6,0)
      call z_add(543,bo(22),an(128),dh(1),n-5,n-8,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(543,bo(8),an(2),dh(3),n-6,n-9,n-1,0)
      call z_add(543,bo(22),an(128),dh(3),n-6,n-9,n-1,0)
!
    else if (resname.eq.'URE') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(446,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(447,bo(94),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(448,bo(58),an(101),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn2 = n
      call z_add(448,bo(58),an(34),dh(3),n-3,n-1,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(449,bo(8),an(2),dh(1),n-2,n-4,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(449,bo(8),an(2),dh(3),n-3,n-5,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(449,bo(8),an(2),dh(1),n-3,n-6,n-5,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(449,bo(8),an(2),dh(3),n-4,n-7,n-1,0)
!
    else if (resname.eq.'FOA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(491,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(492,bo(13),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(493,bo(14),an(4),dh(1),n-2,n-1,0,0)
      if (ua_model.eq.0) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(495,bo(2),an(129),an(129),n-3,n-2,n-1,1)
      else
        shf = 1
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(494,bo(8),an(5),dh(3),n-2+shf,n-4+shf,n-3+shf,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(494,bo(8),an(5),dh(3),n-3+shf,n-5+shf,n-1,0)
!
    else if (resname.eq.'NMF') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      c1 = n
      call z_add(454,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(455,bo(13),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(456,bo(14),an(4),dh(1),n-2,n-1,0,0)
      wline(i) = n
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(458,bo(20),an(41),omega(i),n-1,n-3,n-2,0)
      if (ua_model.eq.0) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(460,bo(2),an(35),dh(3),n-4,n-3,n-2,0)
      else
        shf = 1
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(457,bo(8),an(5),an(33),n-3+shf,n-5+shf,n-2+shf,-1)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(459,bo(2),an(3),chi(1,i),n-3+shf,n-4+shf,n-6+shf,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(459,bo(2),an(3),an(3),n-4+shf,n-5+shf,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(459,bo(2),an(3),an(3),n-5+shf,n-6+shf,n-2,-1)
      end if
!
    else if (resname.eq.'NMA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(461,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(462,bo(13),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(463,bo(14),an(4),dh(1),n-2,n-1,0,0)
      wline(i) = n
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(467,bo(20),an(41),omega(i),n-1,n-3,n-2,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(465,bo(4),an(26),dh(3),n-4,n-3,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(464,bo(8),an(5),an(33),n-3,n-5,n-2,-1)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(468,bo(2),an(3),chi(1,i),n-3,n-4,n-6,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(468,bo(2),an(3),an(3),n-4,n-5,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(468,bo(2),an(3),an(3),n-5,n-6,n-2,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(466,bo(2),an(3),chi(2,i),n-5,n-9,n-8,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(466,bo(2),an(3),an(3),n-6,n-10,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(466,bo(2),an(3),an(3),n-7,n-11,n-2,-1)
      end if
!
    else if (resname.eq.'DMA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(496,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(497,bo(13),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(498,bo(14),an(4),dh(1),n-2,n-1,0,0)
!     the Z-CNT
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(501,bo(20),an(29),dh(1),n-1,n-3,n-2,0)
!     and the E-CNT
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(503,bo(20),an(29),dh(3),n-2,n-4,n-1,0)
!     and the CCT
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(499,bo(4),an(25),dh(3),n-5,n-4,n-3,0)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(502,bo(2),an(3),chi(1,i),n-3,n-4,n-6,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(502,bo(2),an(3),an(3),n-4,n-5,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(502,bo(2),an(3),an(3),n-5,n-6,n-2,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(504,bo(2),an(3),chi(2,i),n-5,n-7,n-9,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(504,bo(2),an(3),an(3),n-6,n-8,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(504,bo(2),an(3),an(3),n-7,n-9,n-2,-1)
        chiline(3,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(500,bo(2),an(3),chi(3,i),n-7,n-12,n-10,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(500,bo(2),an(3),an(3),n-8,n-13,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(500,bo(2),an(3),an(3),n-9,n-14,n-2,-1)
      end if
!
    else if (resname.eq.'ACA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(469,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(470,bo(13),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(471,bo(14),an(4),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(472,bo(8),an(5),dh(3),n-1,n-3,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(472,bo(8),an(5),dh(3),n-2,n-4,n-1,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(473,bo(4),an(25),dh(3),n-5,n-4,n-3,0)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(474,bo(2),an(3),chi(1,i),n-1,n-6,n-5,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(474,bo(2),an(3),an(3),n-2,n-7,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(474,bo(2),an(3),an(3),n-3,n-8,n-2,-1)
      end if
!
    else if (resname.eq.'PPA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(475,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(476,bo(13),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(477,bo(14),an(4),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(478,bo(8),an(5),dh(3),n-1,n-3,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(478,bo(8),an(5),dh(3),n-2,n-4,n-1,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(479,bo(4),an(25),dh(3),n-5,n-4,n-3,0)
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(481,bo(5),an(30),chi(1,i),n-1,n-6,n-5,0)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(480,bo(2),an(3),an(3),n-2,n-7,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(480,bo(2),an(3),an(3),n-3,n-8,n-2,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(482,bo(2),an(3),chi(2,i),n-3,n-4,n-9,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(482,bo(2),an(3),an(3),n-4,n-5,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(482,bo(2),an(3),an(3),n-5,n-6,n-2,-1)
      end if
!
    else if (resname.eq.'SPC') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(450,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(451,bo(22),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(451,bo(22),an(3),dh(1),n-2,n-1,0,0)
!
    else if (resname.eq.'T3P') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(452,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(453,bo(6),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(453,bo(6),an(7),dh(1),n-2,n-1,0,0)
!
    else if (resname.eq.'T4P') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(483,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(484,bo(6),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(484,bo(6),an(7),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(485,bo(21),an(6),dh(1),n-3,n-2,n-1,0)
!
    else if (resname.eq.'T4E') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1050,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1051,bo(6),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1051,bo(6),an(7),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1052,bo(100),an(6),dh(1),n-3,n-2,n-1,0)
!
    else if (resname.eq.'T5P') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(486,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(487,bo(6),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(487,bo(6),an(7),dh(1),n-2,n-1,0,0)
!     WARNING!!! what is bond angle charge-site to H? has to make cs to cs tetrahedral
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(488,bo(15),an(8),an(8),n-3,n-2,n-1,1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(488,bo(15),an(8),an(8),n-4,n-3,n-2,-1)
!
    else if (resname.eq.'CH4') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(489,bo(1),an(1),dh(1),0,0,0,0)
      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(490,bo(2),an(1),dh(1),n-1,0,0,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(490,bo(2),an(3),dh(1),n-2,n-1,0,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(490,bo(2),an(3),an(3),n-3,n-2,n-1,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(490,bo(2),an(3),an(3),n-4,n-3,n-2,1)
      end if
!
    else if (resname.eq.'PRP') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(963,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(965,bo(5),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(963,bo(5),an(104),dh(1),n-1,n-2,0,0)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(964,bo(2),an(3),chi(1,i),n-3,n-2,n-1,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(964,bo(2),an(3),an(3),n-4,n-3,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(964,bo(2),an(3),an(3),n-5,n-4,n-2,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(966,bo(2),an(3),an(3),n-5,n-6,n-4,1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(966,bo(2),an(3),an(3),n-6,n-7,n-5,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(964,bo(2),an(3),chi(2,i),n-6,n-7,n-8,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(964,bo(2),an(3),an(3),n-7,n-8,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(964,bo(2),an(3),an(3),n-8,n-9,n-2,-1)
      end if
!
    else if (resname.eq.'2MN') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1044,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1042,bo(98),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1044,bo(98),an(32),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1043,bo(11),an(3),an(3),n-2,n-3,n-1,1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1043,bo(11),an(3),an(3),n-3,n-4,n-2,-1)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1045,bo(2),an(3),chi(1,i),n-5,n-4,n-3,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1045,bo(2),an(3),an(3),n-6,n-5,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1045,bo(2),an(3),an(3),n-7,n-6,n-2,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1045,bo(2),an(3),chi(2,i),n-6,n-7,n-8,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1045,bo(2),an(3),an(3),n-7,n-8,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1045,bo(2),an(3),an(3),n-8,n-9,n-2,-1)
      end if
!
    else if (resname.eq.'NBU') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(967,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(969,bo(5),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(969,bo(5),an(104),dh(1),n-1,n-2,0,0)
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(967,bo(5),an(104),chi(1,i),n-1,n-2,n-3,0)
      if (ua_model.eq.0) then
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(968,bo(2),an(3),chi(2,i),n-4,n-3,n-2,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(968,bo(2),an(3),an(3),n-5,n-4,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(968,bo(2),an(3),an(3),n-6,n-5,n-2,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(970,bo(2),an(3),an(3),n-6,n-7,n-5,1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(970,bo(2),an(3),an(3),n-7,n-8,n-6,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(970,bo(2),an(3),an(3),n-7,n-6,n-8,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(970,bo(2),an(3),an(3),n-8,n-7,n-9,-1)
        chiline(3,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(968,bo(2),an(3),chi(3,i),n-8,n-9,n-10,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(968,bo(2),an(3),an(3),n-9,n-10,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(968,bo(2),an(3),an(3),n-10,n-11,n-2,-1)
      end if
!
    else if (resname.eq.'IBU') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(971,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(973,bo(5),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(971,bo(5),an(104),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(971,bo(5),an(104),an(104),n-2,n-3,n-1,1)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),chi(1,i),n-4,n-3,n-2,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),an(3),n-5,n-4,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),an(3),n-6,n-5,n-2,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(974,bo(2),an(3),an(3),n-6,n-7,n-5,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),chi(2,i),n-6,n-7,n-5,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),an(3),n-7,n-8,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),an(3),n-8,n-9,n-2,-1)
        chiline(3,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),chi(3,i),n-8,n-10,n-9,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),an(3),n-9,n-11,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(972,bo(2),an(3),an(3),n-10,n-12,n-2,-1)
      end if
!
    else if (resname.eq.'EMT') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(994,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(996,bo(81),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(997,bo(81),an(109),dh(1),n-1,n-2,0,0)
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(999,bo(5),an(108),chi(1,i),n-1,n-2,n-3,0)
      if (ua_model.eq.0) then
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(995,bo(2),an(3),chi(2,i),n-4,n-3,n-2,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(995,bo(2),an(3),an(3),n-5,n-4,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(995,bo(2),an(3),an(3),n-6,n-5,n-2,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(998,bo(2),an(3),an(3),n-5,n-4,n-6,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(998,bo(2),an(3),an(3),n-6,n-5,n-7,-1)
        chiline(3,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1000,bo(2),an(3),chi(3,i),n-6,n-7,n-8,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1000,bo(2),an(3),an(3),n-7,n-8,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1000,bo(2),an(3),an(3),n-8,n-9,n-2,-1)
      end if
!
    else if (resname.eq.'NH4') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(535,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(536,bo(8),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(536,bo(8),an(3),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(536,bo(8),an(3),an(3),n-3,n-2,n-1,-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(536,bo(8),an(3),an(3),n-4,n-3,n-2,1)
!
    else if (resname.eq.'MOH') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(505,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(506,bo(25),an(1),dh(1),n-1,0,0,0)
      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(507,bo(2),an(3),dh(1),n-2,n-1,0,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(507,bo(2),an(3),an(3),n-3,n-2,n-1,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(507,bo(2),an(3),an(3),n-4,n-3,n-2,1)
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(508,bo(26),an(31),chi(1,i),n-4,n-5,n-3,0)
      else
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(508,bo(26),an(31),dh(1),n-1,n-2,0,0)
      end if
!
    else if (resname.eq.'EOH') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(984,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(986,bo(5),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(988,bo(25),an(3),dh(1),n-1,n-2,0,0)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(985,bo(2),an(3),chi(1,i),n-3,n-2,n-1,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(985,bo(2),an(3),an(3),n-4,n-3,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(985,bo(2),an(3),an(3),n-5,n-4,n-2,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(987,bo(2),an(3),an(3),n-5,n-6,n-4,1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(987,bo(2),an(3),an(3),n-6,n-7,n-5,-1)
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(989,bo(26),an(31),chi(2,i),n-6,n-7,n-8,0)
      else
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(989,bo(26),an(31),chi(2,i),n-1,n-2,n-3,0)
      end if
!
    else if (resname.eq.'1MN') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1038,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1040,bo(99),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1039,bo(11),an(3),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1039,bo(11),an(3),an(3),n-3,n-2,n-1,1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1039,bo(11),an(3),an(3),n-4,n-3,n-2,-1)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1041,bo(2),an(3),chi(1,i),n-4,n-5,n-3,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1041,bo(2),an(3),an(3),n-5,n-6,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(1041,bo(2),an(3),an(3),n-6,n-7,n-2,-1)
      end if
!
    else if (resname.eq.'MSH') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(990,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(991,bo(81),an(1),dh(1),n-1,0,0,0)
      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(992,bo(2),an(3),dh(1),n-2,n-1,0,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(992,bo(2),an(3),an(3),n-3,n-2,n-1,-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(992,bo(2),an(3),an(3),n-4,n-3,n-2,1)
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(993,bo(82),an(107),chi(1,i),n-4,n-5,n-3,0)
      else if (ua_model.eq.1) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(993,bo(82),an(107),dh(1),n-1,n-2,0,0)
      end if
!
    else if (resname.eq.'AC-') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(537,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(538,bo(10),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(538,bo(10),an(40),dh(1),n-2,n-1,0,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(539,bo(4),an(10),dh(3),n-3,n-2,n-1,0)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(540,bo(2),an(3),chi(1,i),n-1,n-4,n-3,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(540,bo(2),an(3),an(3),n-2,n-5,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(540,bo(2),an(3),an(3),n-3,n-6,n-2,1)
      end if
!
    else if (resname.eq.'BEN') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(1053,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc21 = n
      call z_add(1053,bo(24),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc22 = n
      call z_add(1053,bo(24),an(2),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(1053,bo(24),an(2),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(1053,bo(24),an(2),dh(1),n-2,n-4,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(1053,bo(24),an(2),dh(1),n-2,n-4,n-5,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1054,bo(12),an(2),dh(3),n-6,n-5,n-4,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1054,bo(12),an(2),dh(3),n-6,n-7,n-4,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(1054,bo(12),an(2),dh(3),n-6,n-8,n-4,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1054,bo(12),an(2),dh(3),n-6,n-8,n-4,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(1054,bo(12),an(2),dh(3),n-6,n-8,n-5,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1054,bo(12),an(2),dh(3),n-6,n-8,n-7,0)
      end if
!
    else if (resname.eq.'TOL') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(975,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(977,bo(27),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc21 = n
      call z_add(978,bo(24),an(2),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc22 = n
      call z_add(978,bo(24),an(2),dh(3),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(980,bo(24),an(2),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(980,bo(24),an(2),dh(1),n-2,n-4,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(982,bo(24),an(2),dh(1),n-2,n-4,n-5,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      if (ua_model.eq.0) then
        shf = 0
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(976,bo(2),an(3),chi(1,i),n-7,n-6,n-5,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(976,bo(2),an(3),an(3),n-8,n-7,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(976,bo(2),an(3),an(3),n-9,n-8,n-2,1)
      else
        shf = 3
      end if
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(979,bo(12),an(2),dh(3),n-8+shf,n-9+shf,n-7+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(979,bo(12),an(2),dh(3),n-8+shf,n-10+shf,n-9+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(981,bo(12),an(2),dh(3),n-8+shf,n-10+shf,n-11+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(981,bo(12),an(2),dh(3),n-8+shf,n-10+shf,n-12+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(983,bo(12),an(2),dh(3),n-8+shf,n-10+shf,n-12+shf,0)
      end if
!
    else if (resname.eq.'PCR') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(509,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(511,bo(27),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc21 = n
      call z_add(512,bo(24),an(2),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc22 = n
      call z_add(512,bo(24),an(2),an(2),n-2,n-3,n-1,-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(514,bo(24),an(2),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(514,bo(24),an(2),dh(1),n-2,n-4,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(516,bo(24),an(2),dh(1),n-2,n-4,n-5,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(517,bo(30),an(2),dh(3),n-1,n-3,n-5,0)
      if (ua_model.eq.0) then
        shf = 0
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(510,bo(2),an(3),chi(1,i),n-8,n-7,n-6,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(510,bo(2),an(3),an(3),n-9,n-8,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(510,bo(2),an(3),an(3),n-10,n-9,n-2,1)
      else
        shf = 3
      end if
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(513,bo(12),an(2),dh(3),n-9+shf,n-10+shf,n-8+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(513,bo(12),an(2),dh(3),n-9+shf,n-11+shf,n-10+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(515,bo(12),an(2),dh(3),n-9+shf,n-11+shf,n-12+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(515,bo(12),an(2),dh(3),n-9+shf,n-11+shf,n-13+shf,0)
      else
        shf = shf + 4
      end if
      if (ua_model.eq.0) then
        chiline(2,i) = n
      else
        chiline(1,i) = n
      end if
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(518,bo(26),an(32),chi(2,i),n-8+shf,n-9+shf,n-11+shf,0)
!
    else if (resname .eq. 'IMD') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(1001,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(1003,bo(73),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      nnd1 = n
      call z_add(1004,bo(86),an(119),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      ncd2 = n
      call z_add(1006,bo(87),an(120),dh(3),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = nnd1
      call z_add(1008,bo(88),an(121),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1010,bo(89),an(122),dh(1),n-2,n-4,n-3,0)
!     just a ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      if (ua_model.eq.0) then
        shf = 0
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1002,bo(2),an(3),chi(1,i),n-6,n-5,n-4,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1002,bo(2),an(3),an(3),n-7,n-6,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1002,bo(2),an(3),an(3),n-8,n-7,n-2,1)
      else
        shf = 3
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = nnd1
      call z_add(1005,bo(3),an(112),dh(3),n-7+shf,n-8+shf,n-6+shf,0)
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd2
        call z_add(1007,bo(78),an(115),dh(3),n-7+shf,n-9+shf,n-8+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = nnd1
        call z_add(1009,bo(78),an(116),dh(3),n-7+shf,n-9+shf,n-10+shf,0)
      end if
!
    else if (resname .eq. 'IME') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(1011,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(1013,bo(73),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      nnd1 = n
      call z_add(1014,bo(90),an(123),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      ncd2 = n
      call z_add(1015,bo(91),an(124),dh(3),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = nnd1
      call z_add(1017,bo(92),an(125),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1019,bo(93),an(126),dh(1),n-2,n-4,n-3,0)
!     just a ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      if (ua_model.eq.0) then
        shf = 0
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1012,bo(2),an(3),chi(1,i),n-6,n-5,n-4,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1012,bo(2),an(3),an(3),n-7,n-6,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1012,bo(2),an(3),an(3),n-8,n-7,n-2,1)
      else
        shf = 3
      end if
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd2
        call z_add(1016,bo(78),an(113),dh(3),n-6+shf,n-8+shf,n-7+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = nnd1
        call z_add(1018,bo(78),an(116),dh(3),n-6+shf,n-8+shf,n-7+shf,0)
      else
        shf = shf + 2
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1020,bo(3),an(112),dh(3),n-6+shf,n-8+shf,n-9+shf,0)
!
    else if (resname.eq.'NAP') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc21 = n
      call z_add(1055,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc22 = n
      call z_add(1055,bo(89),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(1056,bo(102),an(44),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(1056,bo(102),an(44),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(1057,bo(101),an(2),dh(1),n-2,n-4,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(1057,bo(101),an(2),dh(1),n-2,n-4,n-5,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(1056,bo(102),an(102),dh(3),n-6,n-4,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(1056,bo(102),an(102),dh(3),n-6,n-4,n-7,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc21
      call z_add(1057,bo(101),an(2),dh(3),n-2,n-8,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc22
      call z_add(1057,bo(101),an(2),dh(3),n-2,n-8,n-6,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1058,bo(12),an(2),dh(3),n-8,n-10,n-6,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(1058,bo(12),an(2),dh(3),n-8,n-10,n-6,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1059,bo(12),an(2),dh(3),n-8,n-10,n-7,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(1059,bo(12),an(2),dh(3),n-8,n-10,n-9,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1058,bo(12),an(2),dh(3),n-8,n-14,n-6,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(1058,bo(12),an(2),dh(3),n-8,n-14,n-6,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc21
        call z_add(1059,bo(12),an(2),dh(3),n-8,n-10,n-7,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc22
        call z_add(1059,bo(12),an(2),dh(3),n-8,n-10,n-9,0)
      end if
!
    else if (resname .eq. 'MIN') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(1021,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(1023,bo(83),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      ncd1 = n
      call z_add(1024,bo(84),an(118),dh(1),n-1,n-2,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      ncd2 = n
      call z_add(1026,bo(85),an(81),dh(3),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(1027,bo(80),an(114),dh(1),n-2,n-3,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1029,bo(42),an(15),dh(1),n-1,n-3,n-4,0)
!     just a ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-3,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1030,bo(24),an(12),dh(3),n-3,n-1,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1032,bo(24),an(55),dh(1),n-2,n-4,n-1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1034,bo(24),an(2),dh(1),n-2,n-5,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1036,bo(24),an(2),dh(1),n-2,n-4,n-6,0)
!     just a ring closure
      call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
      if (ua_model.eq.0) then
        shf = 0
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1022,bo(2),an(3),chi(1,i),n-10,n-9,n-8,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1022,bo(2),an(3),an(3),n-11,n-10,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(1022,bo(2),an(3),an(3),n-12,n-11,n-2,1)
      else
        shf = 3
      end if
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd1
        call z_add(1025,bo(12),125.35,dh(3),n-11+shf,n-12+shf,n-9+shf,0)
      else
        shf = shf + 1
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(1028,bo(3),122.65,dh(3),n-10+shf,n-12+shf,n-9+shf,0)
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd2
        call z_add(1031,bo(12),an(2),an(2),n-9+shf,n-12+shf,n-7+shf,1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd2
        call z_add(1033,bo(12),an(2),an(2),n-9+shf,n-11+shf,n-7+shf,1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd2
        call z_add(1035,bo(12),118.05,dh(3),n-9+shf,n-11+shf,n-8+shf,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        eqatm(n) = ncd2
        call z_add(1037,bo(12),121.4,dh(3),n-9+shf,n-11+shf,n-10+shf,0)
      end if
!
    else if (resname.eq.'CYT') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(558,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(560,bo(31),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc2 = n
      call z_add(568,bo(33),an(42),dh(1),n-2,n-1,0,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(562,bo(34),an(44),dh(1),n-2,n-3,n-1,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(561,bo(32),an(43),dh(3),n-3,n-4,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(566,bo(35),an(45),dh(1),n-3,n-5,n-4,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nn1
      call z_add(563,bo(36),an(66),dh(1),n-3,n-5,n-6,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nn1
      call z_add(564,bo(36),an(46),dh(3),n-1,n-4,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(559,bo(3),an(47),dh(3),n-8,n-7,n-5,0)
      if (ua_model.lt.2) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(569,bo(2),an(57),dh(3),n-7,n-9,n-8,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(567,bo(2),an(61),dh(3),n-5,n-8,n-10,0)
      else
        shf = 2
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(565,bo(2),an(2),dh(1),n-4+shf,n-5+shf,n-7+shf,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(565,bo(2),an(2),dh(3),n-5+shf,n-6+shf,n-1,0)
!
    else if (resname.eq.'URA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(570,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(572,bo(37),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc2 = n
      call z_add(580,bo(39),an(51),dh(1),n-2,n-1,0,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(574,bo(40),an(52),dh(1),n-2,n-3,n-1,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(573,bo(38),an(54),dh(3),n-3,n-4,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(578,bo(41),an(55),dh(1),n-3,n-5,n-4,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nn1
      call z_add(576,bo(42),an(53),dh(1),n-3,n-5,n-6,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nn1
      call z_add(577,bo(43),an(56),dh(3),n-1,n-4,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(571,bo(3),an(57),dh(3),n-8,n-7,n-5,0)
      if (ua_model.lt.2) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(581,bo(2),an(58),dh(3),n-7,n-9,n-8,0)
      else
        shf = 1
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(575,bo(3),an(59),dh(3),n-7+shf,n-9+shf,n-10+shf,0)
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(579,bo(2),an(60),dh(3),n-6,n-5,n-9,0)
      end if
!
    else if (resname.eq.'THY') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(582,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(584,bo(44),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc2 = n
      call z_add(593,bo(46),an(61),dh(1),n-2,n-1,0,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(586,bo(40),an(62),dh(1),n-2,n-3,n-1,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(585,bo(45),an(64),dh(3),n-3,n-4,n-2,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(590,bo(35),an(65),dh(1),n-3,n-5,n-4,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nn1
      call z_add(588,bo(47),an(63),dh(1),n-3,n-5,n-6,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(591,bo(49),an(67),dh(3),n-2,n-1,n-4,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nn1
      call z_add(589,bo(48),an(66),dh(3),n-2,n-5,n-7,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(583,bo(3),an(68),dh(3),n-9,n-8,n-6,0)
      if (ua_model.lt.2) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(594,bo(2),an(69),dh(3),n-8,n-10,n-9,0)
      else
        shf = 1
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(587,bo(2),an(70),dh(3),n-8+shf,n-10+shf,n-11+shf,0)
      if (ua_model.eq.0) then
        chiline(1,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(592,bo(2),an(3),chi(1,i),n-5,n-7,n-6,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(592,bo(2),an(3),an(3),n-6,n-8,n-1,-1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(592,bo(2),an(3),an(3),n-7,n-9,n-2,1)
      end if
!
    else if (resname.eq.'PUR') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(1071,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(1064,bo(50),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc2 = n
      call z_add(1069,bo(40),an(72),dh(1),n-2,n-1,0,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(1065,bo(51),an(72),dh(1),n-2,n-3,n-1,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(1068,bo(52),an(73),dh(1),n-2,n-4,n-3,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(1063,bo(53),an(74),dh(3),n-4,n-5,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(1061,bo(54),an(75),dh(3),n-1,n-5,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(1066,bo(55),an(76),dh(3),n-4,n-6,n-7,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(1060,bo(35),an(77),dh(1),n-2,n-3,n-7,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(1072,bo(3),an(79),dh(1),n-9,n-8,n-4,0)
      if (ua_model.lt.2) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(1070,bo(2),an(64),dh(3),n-8,n-10,n-9,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc1
        call z_add(1062,bo(2),an(78),dh(3),n-5,n-6,n-10,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc1
        call z_add(1067,bo(2),an(71),dh(3),n-5,n-4,n-9,0)
      end if
!
    else if (resname.eq.'ADE') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(607,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(599,bo(50),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc2 = n
      call z_add(605,bo(40),an(72),dh(1),n-2,n-1,0,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(600,bo(51),an(72),dh(1),n-2,n-3,n-1,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(604,bo(52),an(73),dh(1),n-2,n-4,n-3,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(598,bo(53),an(74),dh(3),n-4,n-5,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(596,bo(54),an(75),dh(3),n-1,n-5,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(601,bo(55),an(76),dh(3),n-4,n-6,n-7,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(595,bo(35),an(77),dh(1),n-2,n-3,n-7,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(602,bo(36),an(71),dh(3),n-2,n-1,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(608,bo(3),an(79),dh(1),n-10,n-9,n-5,0)
      if (ua_model.lt.2) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(606,bo(2),an(64),dh(3),n-9,n-11,n-10,0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc1
        call z_add(597,bo(2),an(78),dh(3),n-6,n-7,n-11,0)
      else
        shf = 2
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(603,bo(3),an(2),dh(3),n-4+shf,n-6+shf,n-5+shf,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(603,bo(3),an(2),dh(3),n-5+shf,n-7+shf,n-1,0)
!
    else if (resname.eq.'GUA') then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nn1 = n
      call z_add(622,bo(1),an(1),dh(1),0,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc1 = n
      call z_add(615,bo(39),an(1),dh(1),n-1,0,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      nc2 = n
      call z_add(620,bo(50),an(80),dh(1),n-2,n-1,0,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(616,bo(56),an(82),dh(1),n-2,n-3,n-1,0) 
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc2
      call z_add(619,bo(57),an(83),dh(1),n-2,n-4,n-3,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(614,bo(58),an(84),dh(3),n-4,n-5,n-3,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(611,bo(59),an(85),dh(3),n-1,n-5,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(617,bo(62),an(86),dh(3),n-4,n-6,n-7,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(609,bo(40),an(87),dh(1),n-2,n-3,n-7,0)
!     ring closure
      call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(612,bo(60),an(66),dh(3),n-3,n-4,n-8,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(618,bo(61),an(81),dh(1),n-3,n-7,n-6,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(623,bo(3),an(88),dh(3),n-11,n-10,n-8,0)
      if (ua_model.lt.2) then
        shf = 0
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        eqatm(n) = nc2
        call z_add(621,bo(2),an(89),dh(3),n-10,n-8,n-12,0)
      else
        shf = 1
      end if
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      eqatm(n) = nc1
      call z_add(610,bo(3),an(90),dh(3),n-5+shf,n-7+shf,n-8+shf,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(613,bo(3),an(2),dh(3),n-5+shf,n-8+shf,n-9+shf,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(613,bo(3),an(2),dh(3),n-6+shf,n-9+shf,n-1,0)
!
!   geometry-wise, ribonucleotides and deoxyribonucleotides are substantially
!   different, so we separate them even though the sequence of atoms is exactly
!   the same in the backbone
    else if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG').OR.(resname.eq.'R5P')) then
      v = seqtyp(i) - 65
      if (resname.eq.'R5P') v = 11
      nnucs(i) = 5
!
      nuci(i,1) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(923,bo(1),an(1),dh(1),0,0,0,0)
      nuci(i,2) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(ptyp(v),bo(63),an(1),dh(1),nuci(i,1),0,0,0)
      nuci(i,3) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(o5typ(v),bo(64),an(91),dh(1),&
 &                             nuci(i,2),nuci(i,1),0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                             nuci(i,2),nuci(i,3),nuci(i,1),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                             nuci(i,2),nuci(i,3),nuci(i,1),1)
      nuci(i,4) = n
      nucsline(2,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                             nuci(i,3),nuci(i,2),nuci(i,1),0)
      nuci(i,5) = n
      nucsline(3,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(c4typ(v),bo(67),an(95),nucs(3,i),&
 &                             nuci(i,4),nuci(i,3),nuci(i,2),0)
      nuci(i,6) = n
      nucsline(4,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(tc3typ(v),bo(69),an(96),nucs(4,i),&
 &                             nuci(i,5),nuci(i,4),nuci(i,3),0)
      if (moltermid(imol,2).eq.1) then
        otze = n
        nucsline(6,i) = n ! hijack
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(to3typ(v),bo(71),an(98),dh(4),&
 &                             nuci(i,6),nuci(i,5),nuci(i,4),0)
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                  nuci(i,5),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                  nuci(i,5),1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                  nuci(i,6),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),&
 &                                                  otze,1)
      end if
!     add the capping hydrogens
      if (moltermid(imol,1).eq.1) then
        nucsline(1,i) = n
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
!       WARNING: temporary values for bo/an
        call z_add(th5typ(v),bo(26),an(31),nucs(1,i),&
 &                                 nuci(i,1),nuci(i,2),nuci(i,3),0)
      else
        write(ilog,*) 'Fatal. Requested unknown N-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
        call fexit()
      end if
      if (moltermid(imol,2).eq.1) then
        nucsline(5,i) = n
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
!       WARNING: temporary values for bo/an
        call z_add(th3typ(v),bo(26),an(31),nucs(5,i),&
 &                                 otze,nuci(i,6),nuci(i,5),0)
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),&
 &                                                   nuci(i,6))
!
    else if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &           (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &           (resname.eq.'DPG').OR.(resname.eq.'D5P')) then
      v = seqtyp(i) - 65
      if (resname.eq.'D5P') v = 12
      nnucs(i) = 5
!
      nuci(i,1) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(923,bo(1),an(1),dh(1),0,0,0,0)
      nuci(i,2) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(ptyp(v),bo(63),an(1),dh(1),nuci(i,1),0,0,0)
      nuci(i,3) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(o5typ(v),bo(64),an(91),dh(1),&
 &                             nuci(i,2),nuci(i,1),0,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                             nuci(i,2),nuci(i,3),nuci(i,1),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                             nuci(i,2),nuci(i,3),nuci(i,1),1)
      nuci(i,4) = n
      nucsline(2,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                             nuci(i,3),nuci(i,2),nuci(i,1),0)
      nuci(i,5) = n
      nucsline(3,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(c4typ(v),bo(68),an(95),nucs(3,i),&
 &                             nuci(i,4),nuci(i,3),nuci(i,2),0)
      nuci(i,6) = n
      nucsline(4,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(tc3typ(v),bo(70),an(97),nucs(4,i),&
 &                             nuci(i,5),nuci(i,4),nuci(i,3),0)
      if (moltermid(imol,2).eq.1) then
        otze = n
        nucsline(6,i) = n ! hijack
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(to3typ(v),bo(72),an(99),dh(5),&
 &                             nuci(i,6),nuci(i,5),nuci(i,4),0)
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                  nuci(i,5),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                  nuci(i,5),1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                  nuci(i,6),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),&
 &                                                  otze,1)
      end if
!     add the capping hydrogens
      if (moltermid(imol,1).eq.1) then
        nucsline(1,i) = n
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
!       WARNING: temporary values for bo/an
        call z_add(th5typ(v),bo(26),an(31),nucs(1,i),&
 &                                 nuci(i,1),nuci(i,2),nuci(i,3),0)
      else
        write(ilog,*) 'Fatal. Requested unknown N-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
        call fexit()
      end if
      if (moltermid(imol,2).eq.1) then
        nucsline(5,i) = n
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
!       WARNING: temporary values for bo/an
        call z_add(th3typ(v),bo(26),an(31),nucs(5,i),&
 &                                 otze,nuci(i,6),nuci(i,5),0)
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),&
 &                                                   nuci(i,6))
!
    else if (resname.eq.'GLY') then
      ni(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(339,bo(1),an(1),dh(1),0,0,0,0)
      cai(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(555,bo(16),an(1),dh(1),ni(i),0,0,0)
      ci(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
!      tauline(i) = n
      atmres(n) = rsmol(imol,1)
      call z_add(361,bo(9),an(17),dh(1),cai(i),ni(i),0,0)
      oi(i) = n
      yline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      dihi = psi(i)-180.0d0
      call z_add(363,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      yline2(i) = n
      call z_add(363,bo(10),an(10),psi(i),&
 &              ci(i),cai(i),ni(i),0)
      if (moltermid(imol,2).eq.1) then
!       do nothing
      else if (moltermid(imol,2).eq.2) then
        write(ilog,*) 'Fatal. Uncharged C-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      if (moltermid(imol,1).eq.1) then
        fline(i) = n
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(342,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(342,bo(11),an(3),an(3),ni(i),cai(i),n-1,1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(342,bo(11),an(3),an(3),ni(i),cai(i),n-2,-1)
        if (ua_model.eq.0) then
          at(i)%nsc = at(i)%nsc + 1
          at(i)%sc(at(i)%nsc) = n
          atmres(n) = rsmol(imol,1)
          call z_add(344,bo(2),an(9),an(36),cai(i),ni(i),ci(i),-chiral(i))
          at(i)%nsc = at(i)%nsc + 1
          at(i)%sc(at(i)%nsc) = n
          atmres(n) = rsmol(imol,1)
          call z_add(344,bo(2),an(9),an(36),cai(i),ni(i),ci(i),chiral(i))
        end if
      else if (moltermid(imol,1).eq.2) then
        write(ilog,*) 'Fatal. Uncharged N-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown N-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
        call fexit()
      end if
!
    else if (resname.eq.'PRO') then
      ni(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(345,bo(1),an(1),dh(1),0,0,0,0)
      cai(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(556,bo(18),an(1),dh(1),ni(i),0,0,0)
      ci(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
!      tauline(i) = n
      atmres(n) = rsmol(imol,1)
      call z_add(367,bo(9),an(15),dh(1),cai(i),ni(i),0,0)
      oi(i) = n
      yline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      dihi = psi(i)-180.0d0
      call z_add(368,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
      yline2(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(368,bo(10),an(10),psi(i),&
 &                ci(i),cai(i),ni(i),0)
      if (moltermid(imol,2).eq.1) then
!       do nothing
      else if (moltermid(imol,2).eq.2) then
        write(ilog,*) 'Fatal. Uncharged C-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      fline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(348,bo(11),an(3),phi(i),ni(i),cai(i),ci(i)&
 &,0)
      if (moltermid(imol,1).eq.1) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(348,bo(11),an(3),an(3),ni(i),cai(i),n-1,1)
        if (ua_model.eq.0) then
          at(i)%nsc = at(i)%nsc + 1
          at(i)%sc(at(i)%nsc) = n
          atmres(n) = rsmol(imol,1)
          call z_add(350,bo(23),an(18),an(38),cai(i),ni(i),ci(i),-chiral(i))
        end if
      else if (moltermid(imol,1).eq.2) then
        write(ilog,*) 'Fatal. Uncharged N-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown N-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
        call fexit()
      end if
      call sidechain(resname,i,cai(i),ni(i),ci(i))
!
    else if (resname.eq.'AIB') then
      ni(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(544,bo(1),an(1),dh(1),0,0,0,0)
      cai(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(557,bo(29),an(1),dh(1),ni(i),0,0,0)
      ci(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(551,bo(27),an(23),dh(1),cai(i),ni(i),0,0)
      oi(i) = n
      yline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      dihi = psi(i)-180.0d0
      call z_add(553,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
      yline2(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(553,bo(10),an(10),psi(i),&
 &                ci(i),cai(i),ni(i),0)
      if (moltermid(imol,2).eq.1) then
!       do nothing
      else if (moltermid(imol,2).eq.2) then
        write(ilog,*) 'Fatal. Uncharged C-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      fline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(547,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(547,bo(11),an(3),an(37),ni(i),cai(i),n-1,1)
      if (moltermid(imol,1).eq.1) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(547,bo(11),an(3),an(37),ni(i),cai(i),n-2,-1)
      else if (moltermid(imol,1).eq.2) then
        write(ilog,*) 'Fatal. Uncharged N-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown N-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
        call fexit()
      end if
      call sidechain (resname,i,cai(i),ni(i),ci(i))
!
    else if ((resname.eq.'ALA').OR.(resname.eq.'MET').OR.&
 &           (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &           (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &           (resname.eq.'ABA').OR.&
 &           (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &           (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &           (resname.eq.'TRP').OR.(resname.eq.'PHE').OR.&
 &           (resname.eq.'TYR').OR.(resname.eq.'HIP').OR.&
 &           (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &           (resname.eq.'GLU').OR.(resname.eq.'ASP').OR.&
 &           (resname.eq.'HID').OR.(resname.eq.'HIE').OR.&
 &           (resname.eq.'ARG').OR.(resname.eq.'LYS').OR.&
 &           (resname.eq.'DAB').OR.(resname.eq.'ORN').OR.&
 &           (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &           (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &           (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &           (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &           (resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &           (resname.eq.'KM2').OR.(resname.eq.'KM3')) then
      ni(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(333,bo(1),an(1),dh(1),0,0,0,0)
      cai(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(554,bo(16),an(1),dh(1),ni(i),0,0,0)
      ci(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
!      tauline(i) = n
      atmres(n) = rsmol(imol,1)
      call z_add(355,bo(9),an(17),dh(1),cai(i),ni(i),0,0)
      oi(i) = n
      yline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      dihi = psi(i)-180.0d0
      call z_add(357,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
      yline2(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(357,bo(10),an(10),psi(i),&
 &                ci(i),cai(i),ni(i),0)
      if (moltermid(imol,2).eq.1) then
!       do nothing
      else if (moltermid(imol,2).eq.2) then
        write(ilog,*) 'Fatal. Uncharged C-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown C-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
        call fexit()
      end if
      fline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(336,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(336,bo(11),an(3),an(3),ni(i),cai(i),n-1,1)
      if (moltermid(imol,1).eq.1) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = rsmol(imol,1)
        call z_add(336,bo(11),an(3),an(3),ni(i),cai(i),n-2,-1)
        if (ua_model.eq.0) then
          at(i)%nsc = at(i)%nsc + 1
          at(i)%sc(at(i)%nsc) = n
          atmres(n) = rsmol(imol,1)
          call z_add(338,bo(2),an(9),an(36),cai(i),ni(i),ci(i),-chiral(i))
        end if
      else if (moltermid(imol,1).eq.2) then
        write(ilog,*) 'Fatal. Uncharged N-terminus is currently not &
 &supported for free amino acids (molecule ',imol,').'
        call fexit()
      else
        write(ilog,*) 'Fatal. Requested unknown N-terminus type for &
 &molecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
        call fexit()
      end if
      call sidechain(resname,i,cai(i),ni(i),ci(i))
!
    else if (seqtyp(i).eq.-1) then
      do j=n,n+natres(i)-1
        at(i)%bb(j-n+1) = j
      end do
      at(i)%nbb = natres(i)
      atmres(n:(n+natres(i)-1)) = i
      n = n + natres(i)
    end if
!
    atmol(imol,2) = n - 1
!
    if (nchi(i).gt.0) then
      chilst%nr = chilst%nr + 1
      reschi(chilst%nr) = i
    end if
!
    if ((fline(i).gt.0).OR.(yline(i).gt.0)) then
      fylst%nr = fylst%nr + 1
      resfy(fylst%nr) = i
    end if
!
    if (wline(i).gt.0) then
      wlst%nr = wlst%nr + 1
      resw(wlst%nr) = i
    end if
!
    if (nnucs(i).gt.0) then
      nuclst%nr = nuclst%nr + 1
      resnuc(nuclst%nr) = i
      if (aminopolty(seqtyp(i)).eq.'N') then
        nucpuclst%nr = nucpuclst%nr + 1
        resnucpuc(nucpuclst%nr) = i
      end if
    end if
!
    if (pucline(i).gt.0) then
      puclst%nr = puclst%nr + 1
      respuc(puclst%nr) = i
    end if
!
    return
!
  end if
!
! build the first residue as a glycine
!
  if (resname .eq. 'GLY') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(339,bo(1),an(1),dh(1),0,0,0,0)
    cai(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(340,bo(16),an(1),dh(1),ni(i),0,0,0)
    ci(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
!    tauline(i) = n
    atmres(n) = rsmol(imol,1)
    call z_add(341,bo(9),an(17),dh(1),cai(i),ni(i),0,0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    dihi = psi(i)-180.0d0
    call z_add(343,bo(19),an(22),dihi,ci(i),cai(i),ni(i),0)
    if (moltermid(imol,1).eq.1) then
      fline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(342,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(342,bo(11),an(3),an(3),ni(i),cai(i),n-1,1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(342,bo(11),an(3),an(3),ni(i),cai(i),n-2,-1)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(344,bo(2),an(9),an(36),cai(i),ni(i),ci(i),&
 &                                                       -chiral(i))
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(344,bo(2),an(9),an(36),cai(i),ni(i),ci(i),&
 &                                                        chiral(i))
      end if
    else if (moltermid(imol,1).eq.2) then
      fline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(342,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(342,bo(11),an(3),an(37),ni(i),cai(i),n-1,1)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(344,bo(2),an(3),an(3),&
 &                cai(i),ni(i),ci(i),-chiral(i))
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(344,bo(2),an(3),dh(2),&
 &              cai(i),ni(i),n-1,0)
      end if
    else
      write(ilog,*) 'Fatal. Requested unknown N-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
      call fexit()
    end if
!
!   build the first residue as a proline
!
  else if (resname .eq. 'PRO') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(345,bo(1),an(1),dh(1),0,0,0,0)
    cai(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(346,bo(18),an(1),dh(1),ni(i),0,0,0)
    ci(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
!    tauline(i) = n
    atmres(n) = rsmol(imol,1)
    call z_add(347,bo(9),an(17),dh(1),cai(i),ni(i),0,0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    dihi = psi(i)-180.0d0
    call z_add(349,bo(19),an(21),dihi,ci(i),cai(i),ni(i),0)
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(348,bo(11),an(3),phi(i),ni(i),cai(i),ci(i)&
 &,0)
    if (moltermid(imol,1).eq.1) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(348,bo(11),an(3),an(3),ni(i),cai(i),n-1,1)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(350,bo(23),an(18),an(38),&
 &              cai(i),ni(i),ci(i),-chiral(i))
      end if
    else if (moltermid(imol,1).eq.2) then
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(350,bo(23),an(18),an(38),&
 &              cai(i),ni(i),ci(i),-chiral(i))
      end if
    else
      write(ilog,*) 'Fatal. Requested unknown N-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
      call fexit()
    end if
    call sidechain(resname,i,cai(i),ni(i),ci(i))
!
!   build the first residue as a methylalanine
!
  else if (resname .eq. 'AIB') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(544,bo(1),an(1),dh(1),0,0,0,0)
    cai(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(545,bo(29),an(1),dh(1),ni(i),0,0,0)
    ci(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(546,bo(27),an(23),dh(1),cai(i),ni(i),0,0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    dihi = psi(i)-180.0d0
    call z_add(548,bo(19),an(20),dihi,ci(i),cai(i),ni(i),0)
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(547,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(547,bo(11),an(3),an(37),ni(i),cai(i),n-1,1)
    if (moltermid(imol,1).eq.1) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(547,bo(11),an(3),an(37),ni(i),cai(i),n-2,-1)
    else
      write(ilog,*) 'Fatal. Requested unknown N-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
      call fexit()
    end if
    call sidechain (resname,i,cai(i),ni(i),ci(i))
!
!   build the first residue as an N-terminal acetyl group
!   
  else if (resname .eq. 'ACE') then
    cai(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(372,bo(1),an(1),dh(1),0,0,0,0)
    ci(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(374,bo(9),an(1),dh(1),n-1,0,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = rsmol(imol,1)
    oi(i) = n
    call z_add(375,bo(19),an(2),dh(1),n-1,n-2,0,0)
    if (ua_model.eq.0) then
!      ni(i) = n
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
!    yline(i) = n
      psi(i) = 180.0d0
      call z_add(373,bo(2),an(3),psi(i)-180.0,n-3,n-2,n-1,0)
!    at(i)%nbb = at(i)%nbb + 1
!    at(i)%bb(at(i)%nbb) = n
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(373,bo(2),an(3),an(3),n-4,n-3,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,1)
      call z_add(373,bo(2),an(3),an(3),n-5,n-4,n-2,-1)
    end if
!
!   build the first residue as an N-terminal formyl group
!
  else if (resname .eq. 'FOR') then
    ci(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(376,bo(1),an(1),dh(1),0,0,0,0)
    oi(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(378,bo(19),an(1),dh(1),n-1,0,0,0)
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(377,bo(2),an(2),dh(1),n-2,n-1,0,0)
    end if
!
!   build the first residue for 5'-P-ribonucleotides
!
  else if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &         (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &         (resname.eq.'RPG').OR.(resname.eq.'R5P')) then
    v = seqtyp(i) - 65
    if (resname.eq.'R5P') v = 11
    nnucs(i) = 5
!
    nuci(i,1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(923,bo(1),an(1),dh(1),0,0,0,0)
    nuci(i,2) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(ptyp(v),bo(63),an(1),dh(1),nuci(i,1),0,0,0)
    nuci(i,3) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(o5typ(v),bo(64),an(91),dh(1),&
 &                       nuci(i,2),nuci(i,1),0,0) 
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                       nuci(i,2),nuci(i,3),nuci(i,1),-1)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                       nuci(i,2),nuci(i,3),nuci(i,1),1)
    nuci(i,4) = n
    nucsline(2,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                       nuci(i,3),nuci(i,2),nuci(i,1),0)
    nuci(i,5) = n
    nucsline(3,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c4typ(v),bo(67),an(95),nucs(3,i),&
 &                       nuci(i,4),nuci(i,3),nuci(i,2),0)
    nuci(i,6) = n
    nucsline(4,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c3typ(v),bo(69),an(96),nucs(4,i),&
 &                       nuci(i,5),nuci(i,4),nuci(i,3),0)
!
    call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))
!
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                nuci(i,6),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(h3typ(v),bo(2),an(3),dh(6),nuci(i,6),nuci(i,5),&
! &                                                nuci(i,4),0)
      call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),at(i)%sc(1),-1)
    end if
!   add the capping hydrogens
    if (moltermid(imol,1).eq.1) then
      nucsline(1,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!     WARNING: temporary values for bo/an
      call z_add(th5typ(v),bo(26),an(31),nucs(1,i),&
 &                                 nuci(i,1),nuci(i,2),nuci(i,3),0)
    else
      write(ilog,*) 'Fatal. Requested unknown N-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
      call fexit()
    end if
!    call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))
!
!   build the first residue for ribonucleosides (5'-terminal caps for nucleic acids)
!
  else if ((resname.eq.'RIC').OR.(resname.eq.'RIU').OR.&
 &         (resname.eq.'RIT').OR.(resname.eq.'RIA').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'RIB')) then
    v = seqtyp(i) - 77
    if (resname.eq.'RIB') v = 11
!
    nnucs(i) = 3
    nuci(i,1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(957,bo(1),an(1),dh(1),0,0,0,0)
    nuci(i,2) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(958,bo(66),an(1),dh(1),nuci(i,1),0,0,0)
    nuci(i,3) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c4typ(v),bo(67),an(95),dh(1),&
 &                       nuci(i,2),nuci(i,1),0,0)
    nuci(i,4) = n
    nucsline(2,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c3typ(v),bo(69),an(96),nucs(2,i),&
 &                       nuci(i,3),nuci(i,2),nuci(i,1),0)
!
    call nuc_sidechain (resname,i,nuci(i,1),nuci(i,2),nuci(i,3),nuci(i,4))
!
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,2),nuci(i,1),&
 &                                                nuci(i,3),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,2),nuci(i,1),&
 &                                                nuci(i,3),1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,3),nuci(i,2),&
 &                                                nuci(i,4),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(h3typ(v),bo(2),an(3),dh(6),nuci(i,4),nuci(i,3),&
! &                                                nuci(i,2),0)
      call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),at(i)%sc(1),-1)
    end if
!   add the capping hydrogens
    nucsline(1,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
!   WARNING: temporary values for bo/an
    call z_add(959,bo(26),an(31),nucs(1,i),&
 &                               nuci(i,1),nuci(i,2),nuci(i,3),0)
!    call nuc_sidechain (resname,i,nuci(i,1),nuci(i,2),nuci(i,3),nuci(i,4))
!
!  build the first residue for deoxyribonucleotides
!
  else if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &         (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &         (resname.eq.'DPG').OR.(resname.eq.'D5P')) then
    v = seqtyp(i) - 65
    if (resname.eq.'D5P') v = 12
    nnucs(i) = 5
!
    nuci(i,1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(923,bo(1),an(1),dh(1),0,0,0,0)
    nuci(i,2) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(ptyp(v),bo(63),an(1),dh(1),nuci(i,1),0,0,0)
    nuci(i,3) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(o5typ(v),bo(64),an(91),dh(1),&
 &                             nuci(i,2),nuci(i,1),0,0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                             nuci(i,2),nuci(i,3),nuci(i,1),-1)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                             nuci(i,2),nuci(i,3),nuci(i,1),1)
    nuci(i,4) = n
    nucsline(2,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                             nuci(i,3),nuci(i,2),nuci(i,1),0)
    nuci(i,5) = n
    nucsline(3,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c4typ(v),bo(68),an(95),nucs(3,i),&
 &                             nuci(i,4),nuci(i,3),nuci(i,2),0)
    nuci(i,6) = n
    nucsline(4,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c3typ(v),bo(70),an(97),nucs(4,i),&
 &                             nuci(i,5),nuci(i,4),nuci(i,3),0)
!
    call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))
!
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                nuci(i,6),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(h3typ(v),bo(2),an(3),dh(7),nuci(i,6),nuci(i,5),&
! &                                             nuci(i,4),0)
      call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),at(i)%sc(1),-1)
    end if
!   add the capping hydrogens
    if (moltermid(imol,1).eq.1) then
      nucsline(1,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!     WARNING: temporary values for bo/an
      call z_add(th5typ(v),bo(26),an(31),nucs(1,i),&
 &                                 nuci(i,1),nuci(i,2),nuci(i,3),0)
    else
      write(ilog,*) 'Fatal. Requested unknown N-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
      call fexit()
    end if
!
!    call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))
!
!   build the first residue for deoxyribonucleosides (5'-terminal caps for nucleic acids)
!
  else if ((resname.eq.'DIC').OR.(resname.eq.'DIU').OR.&
 &         (resname.eq.'DIT').OR.(resname.eq.'DIA').OR.&
 &         (resname.eq.'DIG').OR.(resname.eq.'DIB')) then
    v = seqtyp(i) - 77
    if (resname.eq.'DIB') v = 12
!
    nnucs(i) = 3
    nuci(i,1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(960,bo(1),an(1),dh(1),0,0,0,0)
    nuci(i,2) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(961,bo(66),an(1),dh(1),nuci(i,1),0,0,0)
    nuci(i,3) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c4typ(v),bo(68),an(95),dh(1),&
 &                       nuci(i,2),nuci(i,1),0,0)
    nuci(i,4) = n
    nucsline(2,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(c3typ(v),bo(70),an(97),nucs(2,i),&
 &                       nuci(i,3),nuci(i,2),nuci(i,1),0)
!
    call nuc_sidechain (resname,i,nuci(i,1),nuci(i,2),nuci(i,3),nuci(i,4))
!
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,2),nuci(i,1),&
 &                                                nuci(i,3),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,2),nuci(i,1),&
 &                                                nuci(i,3),1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,3),nuci(i,2),&
 &                                                nuci(i,4),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
!      call z_add(h3typ(v),bo(2),an(3),dh(7),nuci(i,4),nuci(i,3),&
! &                                                nuci(i,2),0)
      call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),at(i)%sc(1),-1)
    end if
!   add the capping hydrogens
    nucsline(1,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
!   WARNING: temporary values for bo/an
    call z_add(962,bo(26),an(31),nucs(1,i),&
 &                               nuci(i,1),nuci(i,2),nuci(i,3),0)
!    call nuc_sidechain (resname,i,nuci(i,1),nuci(i,2),nuci(i,3),nuci(i,4))
!
  else if (seqtyp(i).eq.-1) then
!
    do j=n,n+natres(i)-1
      at(i)%bb(j-n+1) = j
    end do
    at(i)%nbb = natres(i)
    atmres(n:(n+natres(i)-1)) = i
    n = n + natres(i)
!
!   build the first residue for all standard amino acids
!   
  else
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(333,bo(1),an(1),dh(1),0,0,0,0)
    cai(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(334,bo(16),an(1),dh(1),ni(i),0,0,0)
    ci(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
!    tauline(i) = n
    atmres(n) = rsmol(imol,1)
    call z_add(335,bo(9),an(17),dh(1),cai(i),ni(i),0,0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    dihi = psi(i)-180.0d0
    call z_add(337,bo(19),an(22),dihi,ci(i),cai(i),ni(i),0)
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(336,bo(11),an(3),phi(i),ni(i),cai(i),ci(i),0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,1)
    call z_add(336,bo(11),an(3),an(3),ni(i),cai(i),n-1,1)
    if (moltermid(imol,1).eq.1) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,1)
      call z_add(336,bo(11),an(3),an(3),ni(i),cai(i),n-2,-1)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(338,bo(2),an(9),an(36),cai(i),ni(i),ci(i),-chiral(i))
      end if
    else if (moltermid(imol,1).eq.2) then
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = rsmol(imol,1)
        call z_add(338,bo(2),an(9),an(36),cai(i),ni(i),ci(i),-chiral(i))
      end if
    else
      write(ilog,*) 'Fatal. Requested unknown N-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,1),'.'
      call fexit()
    end if
    call sidechain (resname,i,cai(i),ni(i),ci(i))
  end if
!
!
  if ((resname.ne.'ACE').AND.(resname.ne.'FOR').AND.&
 &    (resname.ne.'NME').AND.(resname.ne.'NH2').AND.&
 &    (resname.ne.'RPC').AND.(resname.ne.'RPU').AND.&
 &    (resname.ne.'RPT').AND.(resname.ne.'RPA').AND.&
 &    (resname.ne.'RPG').AND.(resname.ne.'DPC').AND.&
 &    (resname.ne.'DPU').AND.(resname.ne.'DPT').AND.&
 &    (resname.ne.'DPA').AND.(resname.ne.'DPG').AND.&
 &    (resname.ne.'R5P').AND.(resname.ne.'D5P').AND.&
 &    (resname.ne.'RIC').AND.(resname.ne.'RIU').AND.&
 &    (resname.ne.'RIT').AND.(resname.ne.'RIA').AND.&
 &    (resname.ne.'RIG').AND.(resname.ne.'DIC').AND.&
 &    (resname.ne.'DIU').AND.(resname.ne.'DIT').AND.&
 &    (resname.ne.'DIA').AND.(resname.ne.'DIG').AND.&
 &    (resname.ne.'RIB').AND.(resname.ne.'DIB')) then
    if ((fline(i).gt.0).OR.(yline(i).gt.0)) then
      fylst%nr = fylst%nr + 1
      resfy(fylst%nr) = i
    end if
  end if
!
! add all residues with sample-able sidechains to an array
! remember that the definition of sidechain requires this to
! be a branched part of the chain that is not connected to
! anything outside of its own residue
  if (nchi(i).gt.0) then
    chilst%nr = chilst%nr + 1
    reschi(chilst%nr) = i
  end if
!
! add all residues with non-peptide backbones (i.e., nucleic
! acids at the moment) to the corresponding array
  if (nnucs(i).gt.0) then
    nuclst%nr = nuclst%nr + 1
    resnuc(nuclst%nr) = i
    if (aminopolty(seqtyp(i)).eq.'N') then
      nucpuclst%nr = nucpuclst%nr + 1
      resnucpuc(nucpuclst%nr) = i
    end if
  end if
!
  if (pucline(i).gt.0) then
    puclst%nr = puclst%nr + 1
    respuc(puclst%nr) = i
  end if
!
!     
! build atoms for residues in the middle of the chain
!     
  do i = rsmol(imol,1)+1,rsmol(imol,2)-1
!   set residue type and name again
    k = seqtyp(i)
    if ((k.gt.0).AND.(k.le.MAXAMINO)) then
      resname = amino(k)
    else if (k.eq.-1) then
      resname = '???'
    else
      write(ilog,*) 'Fatal. Encountered unsupported residue type in proteus(...). This is an &
 &omission bug.'
      call fexit()
    end if
!
    if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG').OR.(resname.eq.'R5P')) then
      v = seqtyp(i) - 65
      if (resname.eq.'R5P') v = 11
      nnucs(i) = 5
!
      km1 = seqtyp(i-1)
      if ((km1.ge.76).AND.(km1.le.87)) then
        c5m1 =  nuci(i-1,2)
        c4m1 =  nuci(i-1,3)
        c3m1 =  nuci(i-1,4)
        km2 = 3
      else if ((km1.ge.64).AND.(km1.le.75)) then
        c5m1 =  nuci(i-1,4)
        c4m1 =  nuci(i-1,5)
        c3m1 =  nuci(i-1,6)
        km2 = 5
      else if (km1.eq.-1) then
        c5m1 = 0
        c4m1 = 0 
        c3m1 = 0
        km2 = 1
      else
        write(ilog,*) 'Fatal. Unsupported linkage in proteus(...).'
        call fexit()
      end if
      if ((km1.gt.0).AND.((c5m1.le.0).OR.(c4m1.le.0).OR.(c3m1.le.0))) then
        write(ilog,*) 'Fatal. Missing reference(s) in proteus(...).'
        call fexit()
      end if
!
      nuci(i,1) = n
      nucsline(6,i-1) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(o3typ(v),bo(71),an(98),dh(4),&
 &                  c3m1,c4m1,c5m1,0)
      nuci(i,2) = n
      nucsline(km2,i-1) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(ptyp(v),bo(63),an(100),nucs(km2,i-1),&
 &                  nuci(i,1),c3m1,c4m1,0)
      nuci(i,3) = n
      nucsline(1,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(o5typ(v),bo(64),an(91),nucs(1,i),&
 &                  nuci(i,2),nuci(i,1),c3m1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),1)
      nuci(i,4) = n
      nucsline(2,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                  nuci(i,3),nuci(i,2),nuci(i,1),0)
      nuci(i,5) = n
      nucsline(3,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(c4typ(v),bo(67),an(95),nucs(3,i),&
 &                  nuci(i,4),nuci(i,3),nuci(i,2),0)
      nuci(i,6) = n
      nucsline(4,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(c3typ(v),bo(69),an(96),nucs(4,i),&
 &                  nuci(i,5),nuci(i,4),nuci(i,3),0)
!
      call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))

      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                nuci(i,6),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
!        call z_add(h3typ(v),bo(2),an(3),dh(6),nuci(i,6),nuci(i,5),nuci(i,4),0)
        call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),at(i)%sc(1),-1)
      end if
! 
!      call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),&
! &                                                   nuci(i,6))
!
      if (nchi(i).gt.0) then
        chilst%nr = chilst%nr + 1
        reschi(chilst%nr) = i
      end if
!
      if (nnucs(i).gt.0) then
        nuclst%nr = nuclst%nr + 1
        resnuc(nuclst%nr) = i
        if (aminopolty(seqtyp(i)).eq.'N') then
          nucpuclst%nr = nucpuclst%nr + 1
          resnucpuc(nucpuclst%nr) = i
        end if
      end if
!
!    deoxyribonucleotides
!
    else if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &           (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &           (resname.eq.'DPG').OR.(resname.eq.'D5P')) then
      v = seqtyp(i) - 65
      if (resname.eq.'D5P') v = 12
      nnucs(i) = 5
!
      km1 = seqtyp(i-1)
      if ((km1.ge.76).AND.(km1.le.87)) then
        c5m1 =  nuci(i-1,2)
        c4m1 =  nuci(i-1,3)
        c3m1 =  nuci(i-1,4)
        km2 = 3       
      else if ((km1.ge.64).AND.(km1.le.75)) then
        c5m1 =  nuci(i-1,4)
        c4m1 =  nuci(i-1,5)
        c3m1 =  nuci(i-1,6)
        km2 = 5
      else if (km1.eq.-1) then
        c5m1 = 0
        c4m1 = 0 
        c3m1 = 0
        km2 = 1
      else
        write(ilog,*) 'Fatal. Unsupported linkage in proteus(...).'
        call fexit()
      end if
      if ((km1.gt.0).AND.((c5m1.le.0).OR.(c4m1.le.0).OR.(c3m1.le.0))) then
        write(ilog,*) 'Fatal. Missing reference(s) in proteus(...).'
        call fexit()
      end if
!
      nuci(i,1) = n
      nucsline(6,i-1) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(o3typ(v),bo(72),an(99),dh(5),&
 &                  c3m1,c4m1,c5m1,0)
      nuci(i,2) = n
      nucsline(km2,i-1) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(ptyp(v),bo(63),an(100),nucs(km2,i-1),&
 &                  nuci(i,1),c3m1,c4m1,0)
      nuci(i,3) = n
      nucsline(1,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(o5typ(v),bo(64),an(91),nucs(1,i),&
 &                  nuci(i,2),nuci(i,1),c3m1,0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),1)
      nuci(i,4) = n
      nucsline(2,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                  nuci(i,3),nuci(i,2),nuci(i,1),0)
      nuci(i,5) = n
      nucsline(3,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(c4typ(v),bo(68),an(95),nucs(3,i),&
 &                  nuci(i,4),nuci(i,3),nuci(i,2),0)
      nuci(i,6) = n
      nucsline(4,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
      call z_add(c3typ(v),bo(70),an(97),nucs(4,i),&
 &                  nuci(i,5),nuci(i,4),nuci(i,3),0)
!
      call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))
!
      if (ua_model.eq.0) then
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                nuci(i,6),-1)
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
!        call z_add(h3typ(v),bo(2),an(3),dh(7),nuci(i,6),nuci(i,5),&
! &                                                nuci(i,4),0)
        call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),at(i)%sc(1),-1)
      end if
! 
!      call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),nuci(i,6))
!
      if (nchi(i).gt.0) then
        chilst%nr = chilst%nr + 1
        reschi(chilst%nr) = i
      end if
!
      if (nnucs(i).gt.0) then
        nuclst%nr = nuclst%nr + 1
        resnuc(nuclst%nr) = i
        if (aminopolty(seqtyp(i)).eq.'N') then
          nucpuclst%nr = nucpuclst%nr + 1
          resnucpuc(nucpuclst%nr) = i
        end if
      end if
!
!   generic amino acid residues
    else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &          ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &          (seqtyp(i).eq.51).OR.&
 &          ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
!
!     backbone-N
      shf = 0
      if (seqtyp(i).eq.8) then
        if (disulf(i).ne.0) shf = 10
      end if
      ni(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
!     not all residues (in particular caps) might have the right ref. atoms
      if (seqtyp(i-1).eq.28) then
        if (ua_model.eq.0) then
          call z_add(ntyp(k)+shf,bo(7),an(2),dh(3),ci(i-1),at(i-1)%bb(3),oi(i-1),0)
        else
          call z_add(ntyp(k)+shf,bo(7),an(2),dh(1),ci(i-1),oi(i-1),0,0)
        end if
      else if (seqtyp(i-1).eq.27) then
!        call z_add(ntyp(k)+shf,bo(7),an(2),an(2),
! &              ci(i-1),cai(i-1),oi(i-1),-1)
        call z_add(ntyp(k)+shf,bo(7),an(2),dh(3),&
 &              ci(i-1),cai(i-1),oi(i-1),0)
      else
        yline2(i-1) = n
        call z_add(ntyp(k)+shf,bo(7),an(12),psi(i-1),&
 &              ci(i-1),cai(i-1),ni(i-1),0)
      end if
!
!     C-alpha
      cai(i) = n
      wline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i
!     not all residues (in particular caps) might have the right ref. atoms
      if (seqtyp(i-1).eq.28) then
        if (ua_model.eq.0) then
          call z_add(catyp(k)+shf,bo(16),an(11),omega(i),ni(i),ci(i-1),at(i-1)%bb(3),0)
        else
          call z_add(catyp(k)+shf,bo(16),an(11),omega(i),ni(i),ci(i-1),oi(i-1),0)
        end if
      else
        call z_add(catyp(k)+shf,bo(16),an(11),omega(i),&
 &              ni(i),ci(i-1),cai(i-1),0)
      end if
!
!     carbonyl carbon
      ci(i) = n
      fline(i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
!      tauline(i) = n
      atmres(n) = i
      call z_add(ctyp(k)+shf,bo(9),an(17),phi(i),&
 &              cai(i),ni(i),ci(i-1),0)
!
!     carbonyl oxygen
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = i    
      yline(i) = n
      oi(i) = n
      dihi = psi(i)-180.0d0
      call z_add(otyp(k)+shf,bo(19),an(24),dihi,ci(i),cai(i),ni(i),0)
!
!     evtl. amide hydrogen 
      if ((resname.ne.'PRO').AND.(resname.ne.'HYP').AND.&
 &          (resname.ne.'PCA')) then 
        hni(i) = n
        at(i)%nbb = at(i)%nbb + 1
        at(i)%bb(at(i)%nbb) = n
        atmres(n) = i
        if (hntyp(k) .ne. 0) then
          fline2(i) = n
        else
          fline2(i) = 0
        end if
        dihi = phi(i)-180.0d0
        call z_add(hntyp(k)+shf,bo(11),an(14),dihi,&
 &                             ni(i),cai(i),ci(i),0)
      end if
!
!     evtl. alpha-hydrogen
      if (ua_model.eq.0) then
        if (resname.ne.'AIB') then
          at(i)%nsc = at(i)%nsc + 1
          at(i)%sc(at(i)%nsc) = n
          atmres(n) = i
          call z_add(hatyp(k)+shf,bo(2),an(9),an(36),cai(i),ni(i),ci(i),-chiral(i))
        end if
      end if
!
!     in any case, this has to be a valid psi/phi/omega-residue
      fylst%nr = fylst%nr + 1
      resfy(fylst%nr) = i
      wlst%nr = wlst%nr + 1
      resw(wlst%nr) = i
!
!     finally build the sidechain
      call sidechain (resname,i,cai(i),ni(i),ci(i))
!
!     add all residues with sample-able sidechains to an array
!     remember that the definition of sidechain requires this to
!     be a branched part of the chain that is not connected to
!     anything outside of its own residue
      if (nchi(i).gt.0) then
        chilst%nr = chilst%nr + 1
        reschi(chilst%nr) = i
      end if
!
      if (pucline(i).gt.0) then
        puclst%nr = puclst%nr + 1
        respuc(puclst%nr) = i
      end if
!
    else if (seqtyp(i).eq.-1) then
!
      do j=n,n+natres(i)-1
        at(i)%bb(j-n+1) = j
      end do
      at(i)%nbb = natres(i)
      atmres(n:(n+natres(i)-1)) = i
      n = n + natres(i)
!
    else
!
      write(ilog,*) 'Fatal. Encountered residue within chain in prot&
 &eus(...), which is unrecognized. This is either a bug or a feature&
 & which is not yet supported.'
      call fexit()
!
    end if
!
  end do
!     
! set the number and type of the last residue 
! 
  i = rsmol(imol,2)
  k = seqtyp(rsmol(imol,2))
  if ((k.gt.0).AND.(k.le.MAXAMINO)) then
    resname = amino(k)
  else if (k.eq.-1) then
    resname = '???'
  else
    write(ilog,*) 'Fatal. Encountered unsupported residue type in proteus(...). This is an &
 &omission bug.'
    call fexit()
  end if
!
! build the last residue as a glycine
!   
  if (resname .eq. 'GLY') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
!   not all residues (in particular caps) might have the right ref. atoms
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(359,bo(7),an(2),an(2),ci(i-1),at(i-1)%bb(3),oi(i-1),-1)
      else
        call z_add(359,bo(7),an(2),dh(1),ci(i-1),oi(i-1),0,0)
      end if
    else if (seqtyp(i-1).eq.27) then
      call z_add(359,bo(7),an(2),dh(3),&
 &             ci(i-1),cai(i-1),oi(i-1),0)
    else
      yline2(i-1) = n
      call z_add(359,bo(7),an(12),psi(i-1),&
 &             ci(i-1),cai(i-1),ni(i-1),0)
    end if
    cai(i) = n
    wline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(360,bo(16),an(11),omega(i),ni(i),ci(i-1),at(i-1)%bb(3),0)
      else
        call z_add(360,bo(16),an(11),omega(i),ni(i),ci(i-1),oi(i-1),0)
      end if
    else
      call z_add(360,bo(16),an(11),omega(i),&
 &              ni(i),ci(i-1),cai(i-1),0)
    end if
    ci(i) = n
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(361,bo(9),an(17),phi(i),&
 &              cai(i),ni(i),ci(i-1),0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    dihi = psi(i)-180.0d0
    call z_add(363,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    yline2(i) = n
    call z_add(363,bo(10),an(10),psi(i),&
 &              ci(i),cai(i),ni(i),0)
    hni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    fline2(i) = n
    dihi = phi(i)-180.0d0
    call z_add(362,bo(11),an(13),dihi,ni(i),cai(i),ci(i),0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,2)
      call z_add(364,bo(28),an(3),an(36),cai(i),ni(i),ci(i),-chiral(i))
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,2)
      call z_add(364,bo(28),an(3),an(36),cai(i),ni(i),ci(i),chiral(i))
    end if
    if (moltermid(imol,2).eq.1) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
!   
!   build the last residue as a proline
!   
  else if (resname .eq. 'PRO') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
!   not all residues (in particular caps) might have the right ref. atoms
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(365,bo(17),an(2),an(2),ci(i-1),at(i-1)%bb(3),oi(i-1),-1)
      else
        call z_add(365,bo(17),an(2),dh(1),ci(i-1),oi(i-1),0,0)
      end if
    else if (seqtyp(i-1).eq.27) then
      call z_add(365,bo(17),an(2),dh(3),&
 &             ci(i-1),cai(i-1),oi(i-1),0)
    else
      yline2(i-1) = n
      call z_add(365,bo(17),an(19),psi(i-1),&
 &             ci(i-1),cai(i-1),ni(i-1),0)
    end if
    cai(i) = n
    wline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(366,bo(16),an(13),omega(i),ni(i),ci(i-1),at(i-1)%bb(3),0)
      else
        call z_add(366,bo(16),an(13),omega(i),ni(i),ci(i-1),oi(i-1),0)
      end if
    else
      call z_add(366,bo(16),an(13),omega(i),&
 &                ni(i),ci(i-1),cai(i-1),0)
    end if
    ci(i) = n
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(367,bo(27),an(15),phi(i),&
 &                cai(i),ni(i),ci(i-1),0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    dihi = psi(i)-180.0d0
    call z_add(368,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
    yline2(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(368,bo(10),an(10),psi(i),&
 &                ci(i),cai(i),ni(i),0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,2)
      call z_add(369,bo(28),an(3),an(36),cai(i),ni(i),ci(i),-chiral(i))
    end if
    if (moltermid(imol,2).eq.1) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
    call sidechain (resname,i,cai(i),ni(i),ci(i))
!   
!   build the last residue as a methylalanine
!   
  else if (resname .eq. 'AIB') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(549,bo(17),an(2),an(2),ci(i-1),at(i-1)%bb(3),oi(i-1),-1)
      else
        call z_add(549,bo(17),an(2),dh(1),ci(i-1),oi(i-1),0,0)
      end if
    else if (seqtyp(i-1).eq.27) then
      call z_add(549,bo(17),an(2),dh(3),&
 &             ci(i-1),cai(i-1),oi(i-1),0)
    else
      yline2(i-1) = n
      call z_add(549,bo(17),an(19),psi(i-1),&
 &                ci(i-1),cai(i-1),ni(i-1),0)
    end if
    cai(i) = n
    wline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(550,bo(16),an(13),omega(i),ni(i),ci(i-1),at(i-1)%bb(3),0)
      else
        call z_add(550,bo(16),an(13),omega(i),ni(i),ci(i-1),oi(i-1),0)
      end if
    else
      call z_add(550,bo(16),an(13),omega(i),&
 &                ni(i),ci(i-1),cai(i-1),0)
    end if
    ci(i) = n
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(551,bo(27),an(15),phi(i),&
 &                cai(i),ni(i),ci(i-1),0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    dihi = psi(i)-180.0d0
    call z_add(553,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
    yline2(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(553,bo(10),an(10),psi(i),&
 &                ci(i),cai(i),ni(i),0)
    hni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    fline2(i) = n
    dihi = phi(i)-180.0d0
    call z_add(552,bo(11),an(13),dihi,ni(i),cai(i),ci(i),0)
    if (moltermid(imol,2).eq.1) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
    call sidechain (resname,i,cai(i),ni(i),ci(i))
!   
!   build the last residue as a C-terminal N-methylamide
!   
  else if (resname .eq. 'NME') then
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    ni(i) = n
!   not all residues (inparticular caps) might have the right ref. atoms
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(379,bo(7),an(2),an(2),ci(i-1),at(i-1)%bb(3),oi(i-1),-1)
      else
        call z_add(379,bo(7),an(2),dh(1),ci(i-1),oi(i-1),0,0)
      end if
    else if (seqtyp(i-1).eq.27) then
       call z_add(379,bo(7),an(2),dh(3),ci(i-1),cai(i-1),oi(i-1),0)
    else
      yline2(i-1) = n
      call z_add(379,bo(7),an(12),psi(i-1),&
 &            ci(i-1),cai(i-1),ni(i-1),0)
    end if
    wline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    cai(i) = n
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(381,bo(16),an(16),omega(i),n-1,ci(i-1),at(i-1)%bb(3),0)
      else
        call z_add(381,bo(16),an(16),omega(i),n-1,ci(i-1),oi(i-1),0)
      end if
    else
      call z_add(381,bo(16),an(16),omega(i),n-1,ci(i-1),&
 &                                               cai(i-1),0)
    end if
    hni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
!   we use a weird dihedral to get rid of a requirement for using something like wline2
    call z_add(380,bo(11),an(14),dh(3),n-2,ci(i-1),n-1,0)
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(382,bo(2),an(3),dh(3),n-2,n-3,ci(i-1),0)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(382,bo(2),an(3),an(3),n-3,n-4,n-1,1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(382,bo(2),an(3),an(3),n-4,n-5,n-2,-1)
    end if
!   
!   build the last residue as a C-terminal amide
!   
  else if (resname .eq. 'NH2') then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
!   not all residues (in particular caps) might have the right ref. atoms
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(383,bo(7),an(2),an(2),ci(i-1),at(i-1)%bb(3),oi(i-1),-1)
      else
        call z_add(383,bo(7),an(2),dh(1),ci(i-1),oi(i-1),0,0)
      end if
    else if (seqtyp(i-1).eq.27) then
      call z_add(383,bo(7),an(2),an(2),&
 &            ci(i-1),cai(i-1),oi(i-1),-1)
    else
      yline2(i-1) = n
      call z_add(383,bo(7),an(12),psi(i-1),&
 &            ci(i-1),cai(i-1),ni(i-1),0)
    end if
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(384,bo(11),an(2),dh(3),&
 &            n-1,ci(i-1),oi(i-1),0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(384,bo(11),an(2),dh(3),&
 &            n-2,ci(i-1),n-1,0)
!   
!   build the last residue for ribonucleotides
!
  else if ((resname.eq.'RPC').OR.(resname.eq.'RPU').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG').OR.(resname.eq.'R5P')) then
    v = seqtyp(i) - 65
    if (resname.eq.'R5P') v = 11
    nnucs(i) = 5
!
    km1 = seqtyp(i-1)
    if ((km1.ge.76).AND.(km1.le.87)) then
      c5m1 =  nuci(i-1,2)
      c4m1 =  nuci(i-1,3)
      c3m1 =  nuci(i-1,4)
      km2 = 3
    else if ((km1.ge.64).AND.(km1.le.75)) then
      c5m1 =  nuci(i-1,4)
      c4m1 =  nuci(i-1,5)
      c3m1 =  nuci(i-1,6)
      km2 = 5
    else if (km1.eq.-1) then
      c5m1 = 0
      c4m1 = 0 
      c3m1 = 0
      km2 = 1
    else
      write(ilog,*) 'Fatal. Unsupported linkage in proteus(...).'
      call fexit()
    end if
    if ((km1.gt.0).AND.((c5m1.le.0).OR.(c4m1.le.0).OR.(c3m1.le.0))) then
      write(ilog,*) 'Fatal. Missing reference(s) in proteus(...).'
      call fexit()
    end if
!
    nuci(i,1) = n
    nucsline(6,i-1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(o3typ(v),bo(71),an(98),dh(4),&
 &                  c3m1,c4m1,c5m1,0)
    nuci(i,2) = n
    nucsline(km2,i-1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(ptyp(v),bo(63),an(100),nucs(km2,i-1),&
 &                  nuci(i,1),c3m1,c4m1,0)
    nuci(i,3) = n
    nucsline(1,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(o5typ(v),bo(64),an(91),nucs(1,i),&
 &                  nuci(i,2),nuci(i,1),c3m1,0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),-1)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),1)
    nuci(i,4) = n
    nucsline(2,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                  nuci(i,3),nuci(i,2),nuci(i,1),0)
    nuci(i,5) = n
    nucsline(3,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(c4typ(v),bo(67),an(95),nucs(3,i),&
 &                  nuci(i,4),nuci(i,3),nuci(i,2),0)
    nuci(i,6) = n
    nucsline(4,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(tc3typ(v),bo(69),an(96),nucs(4,i),&
 &                  nuci(i,5),nuci(i,4),nuci(i,3),0)
    if (moltermid(imol,2).eq.1) then
      otze = n
      nucsline(6,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(to3typ(v),bo(71),an(98),dh(4),&
 &                  nuci(i,6),nuci(i,5),nuci(i,4),0)
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                nuci(i,6),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),&
 &                                                otze,1)
    end if
!   add the capping hydrogens
    if (moltermid(imol,2).eq.1) then
      nucsline(5,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
!     WARNING: temporary values for bo/an
      call z_add(th3typ(v),bo(26),an(31),nucs(5,i),&
 &                                 otze,nuci(i,6),nuci(i,5),0)
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
!
    call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),&
 &                                                   nuci(i,6))
!   
!   build the last residue for deoxyribonucleotides
!
  else if ((resname.eq.'DPC').OR.(resname.eq.'DPU').OR.&
 &           (resname.eq.'DPT').OR.(resname.eq.'DPA').OR.&
 &           (resname.eq.'DPG').OR.(resname.eq.'D5P')) then
    v = seqtyp(i) - 65
    if (resname.eq.'D5P') v = 12
    nnucs(i) = 5
!
    km1 = seqtyp(i-1)
    if ((km1.ge.76).AND.(km1.le.87)) then
      c5m1 =  nuci(i-1,2)
      c4m1 =  nuci(i-1,3)
      c3m1 =  nuci(i-1,4)
      km2 = 3
    else if ((km1.ge.64).AND.(km1.le.75)) then
      c5m1 =  nuci(i-1,4)
      c4m1 =  nuci(i-1,5)
      c3m1 =  nuci(i-1,6)
      km2 = 5
    else if (km1.eq.-1) then
      c5m1 = 0
      c4m1 = 0 
      c3m1 = 0
      km2 = 1
    else
      write(ilog,*) 'Fatal. Unsupported linkage in proteus(...).'
      call fexit()
    end if
    if ((km1.gt.0).AND.((c5m1.le.0).OR.(c4m1.le.0).OR.(c3m1.le.0))) then
      write(ilog,*) 'Fatal. Missing reference(s) in proteus(...).'
      call fexit()
    end if
!
    nuci(i,1) = n
    nucsline(6,i-1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(o3typ(v),bo(72),an(99),dh(5),&
 &                  c3m1,c4m1,c5m1,0)
    nuci(i,2) = n
    nucsline(km2,i-1) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(ptyp(v),bo(63),an(100),nucs(km2,i-1),&
 &                  nuci(i,1),c3m1,c4m1,0)
    nuci(i,3) = n
    nucsline(1,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(o5typ(v),bo(64),an(91),nucs(1,i),&
 &                  nuci(i,2),nuci(i,1),c3m1,0)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),-1)
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(optyp(v),bo(65),an(92),an(93),&
 &                  nuci(i,2),nuci(i,3),nuci(i,1),1)
    nuci(i,4) = n
    nucsline(2,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(c5typ(v),bo(66),an(94),nucs(2,i),&
 &                  nuci(i,3),nuci(i,2),nuci(i,1),0)
    nuci(i,5) = n
    nucsline(3,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(c4typ(v),bo(68),an(95),nucs(3,i),&
 &                  nuci(i,4),nuci(i,3),nuci(i,2),0)
    nuci(i,6) = n
    nucsline(4,i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(tc3typ(v),bo(70),an(97),nucs(4,i),&
 &                  nuci(i,5),nuci(i,4),nuci(i,3),0)
    if (moltermid(imol,2).eq.1) then
      otze = n
      nucsline(6,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(to3typ(v),bo(72),an(99),dh(5),&
 &                  nuci(i,6),nuci(i,5),nuci(i,4),0)
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
    if (ua_model.eq.0) then
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h5typ(v),bo(2),an(3),an(3),nuci(i,4),nuci(i,3),&
 &                                                nuci(i,5),1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h4typ(v),bo(2),an(3),an(3),nuci(i,5),nuci(i,4),&
 &                                                nuci(i,6),-1)
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
      call z_add(h3typ(v),bo(2),an(3),an(3),nuci(i,6),nuci(i,5),&
 &                                                otze,1)
    end if
!   add the capping hydrogens
    if (moltermid(imol,2).eq.1) then
      nucsline(5,i) = n
      at(i)%nbb = at(i)%nbb + 1
      at(i)%bb(at(i)%nbb) = n
      atmres(n) = rsmol(imol,2)
!     WARNING: temporary values for bo/an
      call z_add(th3typ(v),bo(26),an(31),nucs(5,i),&
 &                                 otze,nuci(i,6),nuci(i,5),0)
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
!
    call nuc_sidechain (resname,i,nuci(i,3),nuci(i,4),nuci(i,5),&
 &                                                   nuci(i,6))
!
!  build the last residue for standard amino acids
!
  else if (seqtyp(i).gt.0) then
    ni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(353,bo(7),an(2),an(2),ci(i-1),at(i-1)%bb(3),oi(i-1),-1)
      else
        call z_add(353,bo(7),an(2),dh(1),ci(i-1),oi(i-1),0,0)
      end if
    else if (seqtyp(i-1).eq.27) then
      call z_add(353,bo(7),an(2),dh(3),&
 &             ci(i-1),cai(i-1),oi(i-1),0)
    else
      yline2(i-1) = n
      call z_add(353,bo(7),an(12),psi(i-1),&
 &              ci(i-1),cai(i-1),ni(i-1),0)
    end if
    cai(i) = n
    wline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    if (seqtyp(i-1).eq.28) then
      if (ua_model.eq.0) then
        call z_add(354,bo(16),an(11),omega(i),ni(i),ci(i-1),at(i-1)%bb(3),0)
      else
        call z_add(354,bo(16),an(11),omega(i),ni(i),ci(i-1),oi(i-1),0)
      end if
    else
      call z_add(354,bo(16),an(11),omega(i),&
 &              ni(i),ci(i-1),cai(i-1),0)
    end if
    ci(i) = n
    fline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
!    tauline(i) = n
    atmres(n) = rsmol(imol,2)
    call z_add(355,bo(9),an(17),phi(i),&
 &              cai(i),ni(i),ci(i-1),0)
    oi(i) = n
    yline(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    dihi = psi(i)-180.0d0
    call z_add(357,bo(10),an(10),dihi,ci(i),cai(i),ni(i),0)
    yline2(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    call z_add(357,bo(10),an(10),psi(i),&
 &                ci(i),cai(i),ni(i),0)
    hni(i) = n
    at(i)%nbb = at(i)%nbb + 1
    at(i)%bb(at(i)%nbb) = n
    atmres(n) = rsmol(imol,2)
    fline2(i) = n
    dihi = phi(i)-180.0d0
    call z_add(356,bo(11),an(14),dihi,ni(i),cai(i),ci(i),0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = rsmol(imol,2)
      call z_add(358,bo(12),an(3),an(36),cai(i),ni(i),ci(i),-chiral(i))
    end if
    if (moltermid(imol,2).eq.1) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Requested unknown C-terminus type for mo&
 &lecule ',imol,'. Offending mode is ',moltermid(imol,2),'.'
      call fexit()
    end if
    call sidechain (resname,i,cai(i),ni(i),ci(i))
!
  else if (seqtyp(i).eq.-1) then
!
    do j=n,n+natres(i)-1
      at(i)%bb(j-n+1) = j
    end do
    at(i)%nbb = natres(i)
    atmres(n:(n+natres(i)-1)) = i
    n = n + natres(i)
!
  end if

!
  if ((fline(i).gt.0).OR.(yline(i).gt.0)) then
    fylst%nr = fylst%nr + 1
    resfy(fylst%nr) = i
  end if
!
  if (wline(i).gt.0) then
    wlst%nr = wlst%nr + 1
    resw(wlst%nr) = i
  end if
!
  if (nchi(i).gt.0) then
    chilst%nr = chilst%nr + 1
    reschi(chilst%nr) = i
  end if
!
  if (nnucs(i).gt.0) then
    nuclst%nr = nuclst%nr + 1
    resnuc(nuclst%nr) = i
    if (aminopolty(seqtyp(i)).eq.'N') then
      nucpuclst%nr = nucpuclst%nr + 1
      resnucpuc(nucpuclst%nr) = i
    end if
  end if
!
  if (pucline(i).gt.0) then
    puclst%nr = puclst%nr + 1
    respuc(puclst%nr) = i
  end if
!
!
! set atom limits for molecule
  atmol(imol,2) = n - 1
!  
end
!
!--------------------------------------------------------------------
