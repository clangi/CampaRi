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
! CONTRIBUTIONS: Hoang Tran, Rohit Pappu, Marco Bacci                      !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
!
! resname  3-letter name of current amino accd residue
! i        number of the current amino acid residue
! ica      atom number of alpha carbon in residue i
! nn       atom number of amide nitrogen in residue i
! cc       atom number of carbonyl carbon in residue i
!
!
subroutine sidechain (resname,i,ica,nn,cc)
!
  use atoms
  use fyoc
  use polypep
  use sequen
  use molecule
  use iounit
  use system
!
  implicit none
!
  integer i,ica,nn,cc,ncz,ncg,ncd1,ncd2,nnd1,nne2,nnd2,ccm1,shf
  character(3) resname
! some lengths, angles and dihedrals
  RTYPE dh(100),an(150),bo(100),dih1
!
! see proteus.f for our strategy with regards to eliminating calls of zatom
! with the numbers provided directly as arguments
!
  bo(1) = 1.091d0   ! standard C-H
  bo(2) = 1.53d0    ! standard CH-CH2
  bo(3) = 1.521d0   ! standard CH-CH3
  bo(4) = 1.52d0    ! standard CH2-CH2 + CH2-C in ASN/GLN
  bo(5) = 1.54d0    ! standard CH-CH
  bo(6) = 1.513d0   ! CH2-CH3 in ABA,NVA,NLE,ILE
  bo(7) = 1.503d0   ! CD-CG in PRO
  bo(8) = 1.502d0   ! CB-CG in PHE
  bo(9) = 1.489d0   ! C-N in LYS,ORN,DAB
  bo(10) = 1.02d0   ! N-H in LYS,ORN,DAB,GLN,ASN,HIE,HID,HIP,ARG(all)
  bo(11) = 1.498d0  ! TRP CB-CG
  bo(12) = 1.365d0  ! TRP CB-CD1
  bo(13) = 1.433d0  ! TRP CB-CD2
  bo(14) = 1.374d0  ! TRP CD1-NE1, HIE CD2-NE2
  bo(15) = 1.370d0  ! TRP CE2-NE1 (ring closure is CE2-CD2)
  bo(16) = 1.398d0  ! TRP CD2-CE3
  bo(17) = 1.394d0  ! TRP CE2-CZ2
  bo(18) = 1.382d0  ! TRP CE3-CZ3, PHE/TYR CD-CE, PHE CE-CZ, HID CD2-NE2 
  bo(19) = 1.368d0  ! TRP CZ2-CH2 (ring closure is CH2-CZ3)
  bo(20) = 1.389d0  ! TYR CG-CD
  bo(21) = 1.384d0  ! PHE CG-CD
  bo(22) = 1.512d0  ! TYR CB-CG
  bo(23) = 1.08d0   ! ArC-H
  bo(24) = 1.05d0   ! TRP N-H
  bo(25) = 1.497d0  ! HIE/D/P CB-CG
  bo(26) = 1.356d0  ! HID/E CG-CD2
  bo(27) = 1.371d0  ! HIP CD2-NE2, HIE CG-ND1
  bo(28) = 1.319d0  ! HIE ND1-CE1
  bo(29) = 1.378d0  ! TYR CG-CZ, HIP/D CG-ND1
  bo(30) = 1.354d0  ! HIP CG-CD2
  bo(31) = 1.321d0  ! HIP ND1-CE1
  bo(32) = 1.345d0  ! HID ND1-CE1
  bo(33) = 1.376d0  ! TYR CZ-OH
  bo(34) = 0.97d0   ! TYR O-H
  bo(35) = 0.94d0   ! SER/THR O-H
  bo(36) = 1.417d0  ! SER CB-OG
  bo(37) = 1.433d0  ! THR CB-OG
  bo(38) = 1.808d0  ! CYS CB-SG
  bo(39) = 1.334d0  ! CYS S-H
  bo(40) = 1.822d0  ! CYS C-S
  bo(41) = 1.326d0  ! ARG CZ-NH
  bo(42) = 1.329d0  ! ARG NE-CZ
  bo(43) = 1.46d0   ! ARG CD-NE
  bo(44) = 1.791d0  ! MET S-CH3
  bo(45) = 1.803d0  ! MET CH2-S
  bo(46) = 1.516d0  ! ASN,GLN,GLU CH2-C(O)
  bo(47) = 1.231d0  ! ASN,GLN C-O
  bo(48) = 1.249d0  ! ASP,GLU C-O
  bo(49) = 1.328d0  ! ASN,GLN C-N
  bo(50) = 1.51d0   ! ASP CH2-C(O)
  bo(51) = 1.492d0  ! PRO CB-CG
  bo(52) = 0.0d0
  bo(53) = 1.437d0  ! PRO CD-N
  bo(54) = 1.47d0   ! C-N in neutral, aliphatic amines
  bo(55) = 1.44d0   ! phosphoester C-O
  bo(56) = 1.59d0   ! phosphoester (C-)O-P
  bo(57) = 1.57d0   ! hydrogen phosphate (H-)O-P
  bo(58) = 1.485d0  ! O=-P
  bo(59) = 1.34d0   ! arC-O(-P) ancedotal from high-res PDB
  bo(60) = 1.208d0  ! GLH,ASH CD=OE1 engh/huber. 1.206 in Acta Cryst. (1971). B27, 893-898 
  bo(61) = 1.304d0  ! GLH,ASH CD-OE2 engh/huber. 1.321 in Acta Cryst. (1971). B27, 893-898 
  bo(62) = 0.98d0   ! carboxylic acid O-H (1.011d0 in Acta Cryst. (1971). B27, 893-898)
  bo(63) = 1.30d0   ! estimate for phenolate anion C-O
  bo(64) = 1.75d0   ! estimate for thiolate anion C-S
  bo(65) = 1.494d0  ! NH2C2 in sec. ammonium ions
  bo(66) = 1.498d0  ! NHC3 in tert. ammonium ions  
  bo(67) = 1.500d0  ! NC4 in quat. ammonium ions 
  bo(68) = 1.458d0  ! Sec. amide N-C(al)
  bo(69) = 1.329d0  ! Sec. amide N-C(O)
  bo(70) = 1.525d0  ! Sec. amide C(O)-C(al)
  an(1) = 0.0d0
  an(2) = 109.47d0  ! standard tetrahedral a la HCH
  an(3) = 120.0d0   ! standard sp2-center (idealized)
  an(4) = 118.55d0  ! ASP,GLU CH2CO
  an(5) = 122.9d0   ! ASP,GLU OCO
  an(6) = 120.7d0   ! PHE CBCGCD and CGCDCE, TRP X-CE3-HE3
  an(7) = 110.4d0   ! ALA,ILE,AIB N-CH-CH3
  an(8) = 110.5d0   ! standard N-CA-CH2, ALA C-CA-CB, VAL CA-CB-CG, C-O-H in COOH
  an(9) = 110.1d0   ! standard C-CA-CB
  an(10) = 102.8d0  ! PRO,HYP,PCA ring angles except at N
  an(11) = 113.8d0  ! standard CH-CH2-CH3, PHE,HIE/D/P CACBCG
  an(12) = 113.9d0  ! TYR CACBCG
  an(13) = 120.95d0 ! TYR CBCGCD
  an(14) = 121.2d0  ! TYR CGCDCE
  an(15) = 119.6d0  ! TYR CDCECZ
  an(16) = 119.825d0! TYR CECZOH - slightly adjusted from 119.9 to preserve symmetry (is probably intended to be slightly oop)
  an(17) = 108.5d0  ! TYR CZOHHH
  an(18) = 113.6d0  ! TRP CACBCG
  an(19) = 126.9d0  ! TRP CBCGCD1
  an(20) = 126.8d0  ! TRP CACBCD2
  an(21) = 106.3d0  ! TRP CD1CGCD2
  an(22) = 110.2d0  ! TRP CGCD1NE1
  an(23) = 108.9d0  ! TRP CD1NE1CE2
  an(24) = 118.8d0  ! TRP CE2CD2CE3
  an(25) = 122.4d0  ! TRP CD2CE2CZ2
  an(26) = 118.6d0  ! TRP CD2CE3CZ3, PHE CDCGCD
  an(27) = 117.5d0  ! TRP CE2CZ2CH2
  an(28) = 124.9d0  ! TRP XCD1HD1 and C-C=O in COOH (Acta Cryst. (1971). B27, 893-898)
  an(29) = 125.55d0 ! TRP XNE1HE1
  an(30) = 110.0d0  ! GLY NCAHA2
  an(31) = 107.9d0  ! GLY CCAHA2
  an(32) = 111.5d0  ! standard N-CH-CH
  an(33) = 109.1d0  ! standard C-CH-CH
  an(34) = 106.9d0  ! SER/THR/HYP COH
  an(35) = 116.3d0  ! standard CH-CH2-CH
  an(36) = 110.7d0  ! standard CH2-CH-CH3
  an(37) = 110.8d0  ! standard CH3-CH-CH3
  an(38) = 111.1d0  ! SER CACBO
  an(39) = 109.6d0  ! THR CACBO
  an(40) = 114.4d0  ! CYS CACBS
  an(41) = 96.d0    ! CYS CBCSH
  an(42) = 109.3d0  ! THR OCBCG, HIP CGND1CE1, HIE ND1CGCD2
  an(43) = 122.7d0  ! HIP/D CBCGND1, ASN/GLN OCN 
  an(44) = 131.2d0  ! HIP CBCGCD2
  an(45) = 106.1d0  ! HIP ND1CGCD2
  an(46) = 111.0d0  ! PRO/HYP/PCA CCACB
  an(47) = 107.2d0  ! HIP CGCD2NE2
  an(48) = 125.35d0 ! HIP CGND1HD1
  an(49) = 126.4d0  ! HIP XCD2HD2
  an(50) = 125.8d0  ! HIP XCE1HE1
  an(51) = 125.5d0  ! HIP XNE2HE2 HID CGND1HD1
  an(52) = 129.1d0  ! HID/E CBCGCD2
  an(53) = 105.2d0  ! HID ND1CGCD2
  an(54) = 109.0d0  ! HID CGND1CE1
  an(55) = 109.5d0  ! HID CGCD2NE2
  an(56) = 125.25d0 ! HID CGCD2HD2
  an(57) = 125.3d0  ! HID ND1CE1HE1
  an(58) = 121.6d0  ! HIE CBCGND1
  an(59) = 105.6d0  ! HIE CGND1CE1
  an(60) = 106.5d0  ! HIE CGCD2NE2
  an(61) = 127.2d0  ! HIE CGCD2HD2
  an(62) = 124.1d0  ! HIE ND1CE1HE1
  an(63) = 126.65d0 ! HIE CD2NE2HE2
  an(64) = 114.1    ! general CH-CH2-CH2
  an(65) = 111.3d0  ! general CH2-CH2-CH2
  an(66) = 111.9d0  ! LYS/ORN/DAB CCN
  an(67) = 107.8d0  ! AIB C-C-CH3, ASP CACBCG
  an(68) = 112.6d0  ! GLU,GLN,ASN CACBCG
  an(69) = 120.9d0  ! GLN,ASN CCO
  an(70) = 116.4d0  ! GLN,ASN CCN
  an(71) = 118.1d0  ! TYR CDCGCD
  an(72) = 114.2  ! standard CH2-CH2-CH3
  an(73) = 112.0d0  ! ARG CGCDNE
  an(74) = 124.2d0  ! ARG CDNECZ
  an(75) = 120.15d0 ! ARG NECZNH
  an(76) = 119.7d0  ! ARG NHCZNH
  an(77) = 112.7d0  ! MET CBCGSD HYP CBCGO
  an(78) = 110.9d0  ! MET CGSDCE
  an(79) = 126.0d0  ! PCA CCO (sc)
  an(80) = 121.25d0 ! TRP XCZ2HZ2
  an(81) = 119.2889495554d0 ! TRP XCZ3HZ3
  an(82) = 119.3610504443d0 ! TRP XCH2HH2
  an(83) = 111.85   ! PRO CANCD
  an(84) = 110.3d0  ! CNH in neutral amines
  an(85) = 107.1d0  ! HNH in neutral amines
  an(86) = 106.7d0  ! (C)-O-P-O(H) see Acta Cryst. (1979). B35, 2749-2751 for rough guide on alkyl phosphate
  an(87) = 113.7d0  ! P-O-H
  an(88) = 110.0d0  ! (H-)O-P-=O 
  an(89) = 106.2d0  ! (C-)O-P=-O
  an(90) = 120.0d0  ! alC-O-P
  an(91) = 128.0d0  ! arC-O-P ancedotal from high-res PDB
  an(92) = 121.9d0  ! carboxylic acid O=C-O, Acta Cryst. (1971). B27, 893-898 (from ND)
  an(93) = 113.2d0  ! C-C-O(-H) in COOH, Acta Cryst. (1971). B27, 893-898 (from ND)
  an(94) = 113.0d0  ! C-N-C in secondary ammonium ions
  an(95) = 109.0d0  ! C-N-H in secondary ammonium ions
  an(96) = 108.0d0  ! unused
  an(97) = 111.0d0  ! C-N-C in tertiary ammonium ions
  an(98) = 108.0d0  ! C-N-H in tertiary ammonium ions
  an(99) = 121.7d0  ! C(al)-N-C(O) in sec. amides (from peptide) 
  an(100) = 119.15d0! C-N-H in sec. amides (from peptide)
  dh(1) = 0.0d0
  dh(2) = 180.0d0
  dh(3) = -79.7d0  ! HYP CACBCGO
  dh(4) = -179.30865917086  ! C(-1)NCACB for L-PRO puck chi1 = -28.0 
  dh(5) = 170.49134082914   ! C(-1)NCACB for L-PRO puck chi1 = 31.0

! initialize the number of atoms for each sidechain
  if (at(i)%nsc.ne.1) at(i)%nsc = 0
  nchi(i) = 0
!
! glycine residue  (GLY)
!
  if (resname .eq. 'GLY') then
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(6,bo(1),an(30),an(31),ica,nn,cc,chiral(i))
    end if
!
! alanine residue  (ALA)
!
  else if (resname .eq. 'ALA') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(13,bo(3),an(7),an(8),ica,nn,cc,chiral(i))
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(14,bo(1),an(2),chi(1,i),n-1,ica,nn,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(14,bo(1),an(2),an(2),n-2,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(14,bo(1),an(2),an(2),n-3,ica,n-1,1)
    end if
!
! alpha amino-n-butyric acid residue  (ABA)
!
  else if (resname .eq. 'ABA') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(414,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(416,bo(6),an(11),chi(1,i),n-1,ica,nn,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(415,bo(1),an(2),an(2),n-2,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(415,bo(1),an(2),an(2),n-3,ica,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(417,bo(1),an(2),dh(2),n-3,n-4,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(417,bo(1),an(2),an(2),n-4,n-5,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(417,bo(1),an(2),an(2),n-5,n-6,n-2,-1)
    end if
!
! norvaline residue  (NVA)
!
  else if (resname .eq. 'NVA') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(424,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(426,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(428,bo(6),an(72),chi(2,i),n-1,n-2,ica,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(425,bo(1),an(2),an(2),n-3,ica,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(425,bo(1),an(2),an(2),n-4,ica,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(427,bo(1),an(2),an(2),n-4,n-5,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(427,bo(1),an(2),an(2),n-5,n-6,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(429,bo(1),an(2),dh(2),n-5,n-6,n-7,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(429,bo(1),an(2),an(2),n-6,n-7,n-1,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(429,bo(1),an(2),an(2),n-7,n-8,n-2,1)
    end if
!
! norleucine residue  (NLE)
!
  else if (resname .eq. 'NLE') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(436,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(438,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(440,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(442,bo(6),an(72),chi(3,i),n-1,n-2,n-3,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(437,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(437,bo(1),an(2),an(2),n-5,ica,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(439,bo(1),an(2),an(2),n-5,n-6,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(439,bo(1),an(2),an(2),n-6,n-7,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(441,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(441,bo(1),an(2),an(2),n-7,n-8,n-6,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(443,bo(1),an(2),dh(2),n-7,n-8,n-9,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(443,bo(1),an(2),an(2),n-8,n-9,n-1,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(443,bo(1),an(2),an(2),n-9,n-10,n-2,1)
    end if
! 
! valine residue  (VAL) 
!
  else if (resname .eq. 'VAL') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(21,bo(5),an(32),an(33),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(23,bo(3),an(8),chi(1,i),n-1,ica,nn,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(23,bo(3),an(8),an(8),n-2,ica,n-1,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(22,bo(1),an(2),an(2),n-3,ica,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(24,bo(1),an(2),dh(2),n-3,n-4,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(24,bo(1),an(2),an(2),n-4,n-5,n-1,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(24,bo(1),an(2),an(2),n-5,n-6,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(24,bo(1),an(2),dh(2),n-5,n-7,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(24,bo(1),an(2),an(2),n-6,n-8,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(24,bo(1),an(2),an(2),n-7,n-9,n-2,-1)
    end if
!
! leucine residue  (LEU)
!
  else if (resname .eq. 'LEU') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(31,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(33,bo(2),an(35),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(35,bo(3),an(36),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(35,bo(3),an(36),an(37),n-2,n-3,n-1,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(32,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(32,bo(1),an(2),an(2),n-5,ica,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(34,bo(1),an(2),an(2),n-5,n-6,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(36,bo(1),an(2),dh(2),n-5,n-6,n-7,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(36,bo(1),an(2),an(2),n-6,n-7,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(36,bo(1),an(2),an(2),n-7,n-8,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(36,bo(1),an(2),dh(2),n-7,n-9,n-10,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(36,bo(1),an(2),an(2),n-8,n-10,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(36,bo(1),an(2),an(2),n-9,n-11,n-2,-1)
    end if
!
! isoleucine residue  (ILE) 
!
  else if (resname .eq. 'ILE') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(43,bo(5),an(32),an(33),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(45,bo(2),an(8),chi(1,i),n-1,ica,nn,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(47,bo(3),an(7),an(36),n-2,ica,n-1,1)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(49,bo(6),an(11),chi(2,i),n-2,n-3,ica,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(44,bo(1),an(2),an(2),n-4,ica,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(46,bo(1),an(2),an(2),n-4,n-5,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(46,bo(1),an(2),an(2),n-5,n-6,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(48,bo(1),an(2),dh(2),n-5,n-7,n-6,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(48,bo(1),an(2),an(2),n-6,n-8,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(48,bo(1),an(2),an(2),n-7,n-9,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(50,bo(1),an(2),dh(2),n-7,n-9,n-10,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(50,bo(1),an(2),an(2),n-8,n-10,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(50,bo(1),an(2),an(2),n-9,n-11,n-2,-1)
    end if
!
! serine residue  (SER) 
!
  else if (resname .eq. 'SER') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(57,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(59,bo(36),an(38),chi(1,i),n-1,ica,nn,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(58,bo(1),an(2),an(2),n-2,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(58,bo(1),an(2),an(2),n-3,ica,n-2,-1)
      chiline(2,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(60,bo(35),an(34),chi(2,i),n-3,n-4,ica,0)
    else
      chiline(2,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(60,bo(35),an(34),chi(2,i),n-1,n-2,ica,0)
    end if
!
! phosphoserine residue  (SEP) 
!
  else if (resname .eq. 'SEP') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1165,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1167,bo(55),an(38),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1168,bo(56),an(90),chi(2,i),n-1,n-2,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1170,bo(57),an(86),chi(3,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1169,bo(58),an(88),an(89),n-2,n-1,n-3,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1169,bo(58),an(88),an(89),n-3,n-2,n-4,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1166,bo(1),an(2),an(2),n-6,ica,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1166,bo(1),an(2),an(2),n-7,ica,n-6,-1)
      chiline(4,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1171,bo(35),an(87),chi(4,i),n-5,n-6,n-7,0)
    else
      chiline(4,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1171,bo(35),an(87),chi(4,i),n-3,n-4,n-5,0)
    end if
!
! threonine residue  (THR) 
!
  else if (resname .eq. 'THR') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(67,bo(5),an(32),an(33),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(69,bo(37),an(39),chi(1,i),n-1,ica,nn,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(71,bo(3),an(8),an(42),n-2,ica,n-1,1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(68,bo(1),an(2),an(2),n-3,ica,n-2,-1)
      chiline(2,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(70,bo(35),an(34),chi(2,i),n-3,n-4,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(72,bo(1),an(2),dh(2),n-3,n-5,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(72,bo(1),an(2),an(2),n-4,n-6,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(72,bo(1),an(2),an(2),n-5,n-7,n-2,-1)
    else
      chiline(2,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(70,bo(35),an(34),chi(2,i),n-2,n-3,ica,0)
    end if
!
!  phosphothreonine residue (TPO) 
!
  else if (resname .eq. 'TPO') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1178,bo(5),an(32),an(33),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1180,bo(55),an(39),chi(1,i),n-1,ica,nn,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1181,bo(3),an(8),an(42),n-2,ica,n-1,1)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1183,bo(56),an(90),chi(2,i),n-2,n-3,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1185,bo(57),an(86),chi(3,i),n-1,n-3,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1184,bo(58),an(88),an(89),n-2,n-1,n-4,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1184,bo(58),an(88),an(89),n-3,n-2,n-5,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1179,bo(1),an(2),an(2),n-7,ica,n-6,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1182,bo(1),an(2),dh(2),n-6,n-8,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1182,bo(1),an(2),an(2),n-7,n-9,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1182,bo(1),an(2),an(2),n-8,n-10,n-2,-1)
      chiline(4,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1186,bo(35),an(87),chi(4,i),n-7,n-8,n-10,0)
    else
      chiline(4,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1186,bo(35),an(87),chi(4,i),n-3,n-4,n-6,0)
    end if
!
! cysteine residue  (CYS) 
!
  else if (resname .eq. 'CYS') then
    if (disulf(i) .eq. 0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(79,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(81,bo(38),an(40),chi(1,i),n-1,ica,nn,0)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(80,bo(1),an(2),an(2),n-2,ica,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(80,bo(1),an(2),an(2),n-3,ica,n-2,-1) 
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(82,bo(39),an(41),chi(2,i),n-3,n-4,ica,0)
      else if (ua_model.eq.1) then
        chiline(2,i) = n
        nchi(i) = nchi(i) + 1
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(82,bo(39),an(41),chi(2,i),n-1,n-2,ica,0)
      end if
    else if (disulf(i) .gt. i) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(89,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(91,bo(40),an(40),chi(1,i),n-1,ica,nn,0)
      if (ua_model.eq.0) then
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(90,bo(1),an(2),an(2),n-2,ica,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(90,bo(1),an(2),an(2),n-3,ica,n-2,-1)
      end if
    else if (disulf(i) .lt. i) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(89,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(91,bo(40),an(40),chi(1,i),n-1,ica,nn,0)
      if (ua_model.eq.0) then
!       add the extra bond
        call z_add(-1,bo(52),an(1),dh(1),at(disulf(i))%sc(3),n-1,0,0)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(90,bo(1),an(2),an(2),n-2,ica,n-1,1)
        at(i)%nsc = at(i)%nsc + 1
        at(i)%sc(at(i)%nsc) = n
        atmres(n) = i
        call z_add(90,bo(1),an(2),an(2),n-3,ica,n-2,-1)
      else
!       add the extra bond
        call z_add(-1,bo(52),an(1),dh(1),at(disulf(i))%sc(2),n-1,0,0)
      end if
    end if
!
! deprotonated cysteine residue (CYX, thiolate anion) 
!
  else if (resname .eq. 'CYX') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1122,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1124,bo(64),an(40),chi(1,i),n-1,ica,nn,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1123,bo(1),an(2),an(2),n-2,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1123,bo(1),an(2),an(2),n-3,ica,n-2,-1) 
    end if
!
! proline residue  (PRO)
!
! with the original not Engh/Huber consistent params, the bond length for the old ring closure
! bond (N-C_gamma) is 1.7A. this should be 1.473A instead. with the correct Engh/Huber it's
! 1.59A. i read somewhere about inconsistency of Engh/Huber with pucker, so we'll go with
! the high-resolution X-ray-set used to derive the pucker instead.
! there the angles seem to be the following (bond lengths are the same):
! N-CH1E-CH2E, CH1E-CH2E-CH2P, CH2E-CH2P-CH2P, CH2P-CH2P-N = 102.0-103.0
! CH1E-N-CH2P, C-CH1E-CH2E = 110.0-111.0
! this is based on a very rough analysis but seems to work OK
!
  else if (resname .eq. 'PRO') then
    if (fline(i).gt.0) pucline(i) = fline(i)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    if (i.eq.rsmol(molofrs(i),1)) then
      if (i.eq.rsmol(molofrs(i),2)) then
        if ((moltermid(molofrs(i),1).eq.1).AND.(moltermid(molofrs(i),2).eq.1)) then
          ccm1 = at(i)%bb(6)
        else
          write(ilog,*) 'Fatal. Encountered unsupported C- and/or N-terminus (RES: ',i,') in &
&sidechain().'
          call fexit()
        end if
      else
        if (moltermid(molofrs(i),1).eq.1) then
          ccm1 = at(i)%bb(5)
        else
          write(ilog,*) 'Fatal. Encountered unsupported N-terminus (RES: ',i,') in &
&sidechain().'
          call fexit()
        end if
      end if
!     these are crude: note that with pucker sampling most d.o.f.s will relax
      if ((chiral(i).eq.1).AND.(chi(1,i).lt.0.0)) then
        dih1 = dh(4) + 64.0
      else if ((chiral(i).eq.1).AND.(chi(1,i).gt.0.0)) then
        dih1 = dh(5) + 57.0
      else if ((chiral(i).eq.-1).AND.(chi(1,i).lt.0.0)) then
        dih1 = -dh(5) + 58.0
      else if ((chiral(i).eq.-1).AND.(chi(1,i).gt.0.0)) then
        dih1 = -dh(4) + 53.0
      end if
    else
!     ci always set even for caps
      ccm1 = ci(i-1)
      if ((chiral(i).eq.1).AND.(chi(1,i).lt.0.0)) then
        dih1 = dh(4)
      else if ((chiral(i).eq.1).AND.(chi(1,i).gt.0.0)) then
        dih1 = dh(5)
      else if ((chiral(i).eq.-1).AND.(chi(1,i).lt.0.0)) then
        dih1 = -dh(5)
      else if ((chiral(i).eq.-1).AND.(chi(1,i).gt.0.0)) then
        dih1 = -dh(4)
      end if
    end if
    call z_add(97,bo(2),an(10),dih1,ica,nn,ccm1,0)
!    call z_add(97,bo(2),an(10),an(46),ica,nn,cc,chiral(i))
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    chiline(1,i) = n
!    nchi(i) = nchi(i) + 1
    call z_add(99,bo(51),an(10),chi(1,i),n-1,ica,nn,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    chiline(2,i) = n
!    nchi(i) = nchi(i) + 1
    if ((i.eq.rsmol(molofrs(i),1)).AND.&
&                  (i.eq.rsmol(molofrs(i),2))) then
      call z_add(351,bo(53),an(83),chi(2,i),nn,ica,n-2,0)
    else if (i.eq.rsmol(molofrs(i),1)) then
      call z_add(351,bo(53),an(83),chi(2,i),nn,ica,n-2,0)
    else if (i.eq.rsmol(molofrs(i),2)) then
      call z_add(370,bo(53),an(83),chi(2,i),nn,ica,n-2,0)
    else
      call z_add(101,bo(53),an(83),chi(2,i),nn,ica,n-2,0)
    end if
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(98,bo(1),an(2),an(2),n-3,ica,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(98,bo(1),an(2),an(2),n-4,ica,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(100,bo(1),an(2),an(2),n-4,n-5,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(100,bo(1),an(2),an(2),n-5,n-6,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      if ((i.eq.rsmol(molofrs(i),1)).AND.&
&                  (i.eq.rsmol(molofrs(i),2))) then
        call z_add(352,bo(1),an(2),an(2),n-5,n-6,nn,1)
      else if (i.eq.rsmol(molofrs(i),1)) then
        call z_add(352,bo(1),an(2),an(2),n-5,n-6,nn,1)
      else if (i.eq.rsmol(molofrs(i),2)) then
        call z_add(371,bo(1),an(2),an(2),n-5,n-6,nn,1)
      else
        call z_add(102,bo(1),an(2),an(2),n-5,n-6,nn,1)
      end if
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      if ((i.eq.rsmol(molofrs(i),1)).AND.&
&                  (i.eq.rsmol(molofrs(i),2))) then
        call z_add(352,bo(1),an(2),an(2),n-6,n-7,nn,-1)
      else if (i.eq.rsmol(molofrs(i),1)) then
        call z_add(352,bo(1),an(2),an(2),n-6,n-7,nn,-1)
      else if (i.eq.rsmol(molofrs(i),2)) then
        call z_add(371,bo(1),an(2),an(2),n-6,n-7,nn,-1)
      else
        call z_add(102,bo(1),an(2),an(2),n-6,n-7,nn,-1)
      end if
    end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! hydroxyproline residue  (HYP)
! Note:  this is the B-form poly-Hyp structure from Bansal '79.
!        currently for structure building purposes only.
!        If want to use for simulation, must take into account
!        puckering.  Also, the direction of the O-H hydroxyl bond
!        is optimized for intrachain hydrogen bonding.
!        HT - 16 Jun '04
!
  else if (resname .eq. 'HYP') then
    if (fline(i).gt.0) pucline(i) = fline(i)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(400,bo(2),an(10),an(46),ica,nn,cc,chiral(i))
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    chiline(1,i) = n
!    nchi(i) = nchi(i) + 1
    call z_add(402,bo(51),an(10),chi(1,i),n-1,ica,nn,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    chiline(2,i) = n
!    nchi(i) = nchi(i) + 1
    call z_add(404,bo(7),an(10),chi(2,i),n-1,n-2,ica,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),nn,n-1,0,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(401,bo(1),an(2),an(2),n-3,ica,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(401,bo(1),an(2),an(2),n-4,ica,n-3,-1)
!ccccccccccccccccgamma   carbon attachments      gamma, beta, delta
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(406,bo(37),an(77),dh(3),n-4,n-5,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(403,bo(1),an(2),an(2),n-5,n-6,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
!     careful here: this is the first chi, but because of hijacking we have the initial value in chi(3,i)
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      chi(1,i) = chi(3,i)
      call z_add(407,bo(35),an(34),chi(1,i),n-2,n-6,n-7,0)
!ccccccccccccccccend   gamma carbon attachments
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(405,bo(1),an(2),an(2),n-6,n-7,nn,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(405,bo(1),an(2),an(2),n-7,n-8,nn,-1)
    else
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(406,bo(37),an(77),dh(3),n-2,n-3,ica,0)
      chiline(1,i) = n
      nchi(i) = nchi(i) + 1
      chi(1,i) = chi(3,i)
      call z_add(407,bo(35),an(34),chi(1,i),n-1,n-3,n-4,0)
    end if
!ccccccccc end hydroxyproline cccccccccccccccccccccccccccccccccccccccccc

!
! phenylalanine residue  (PHE) 
!
  else if (resname .eq. 'PHE') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(109,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncg = n
    call z_add(111,bo(8),an(11),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd1 = n
    call z_add(112,bo(21),an(6),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
!    call z_add(112,bo(21),an(6),an(26),n-2,n-3,n-1,1)
    call z_add(112,bo(21),an(6),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(114,bo(18),an(6),dh(2),n-2,n-3,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(114,bo(18),an(6),dh(2),n-2,n-4,n-5,0)
!   note that this doesn't quite work: the chain closure bond and angle will be off
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(116,bo(18),an(3),dh(1),n-2,n-4,n-5,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(110,bo(1),an(2),an(2),n-7,ica,n-6,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(110,bo(1),an(2),an(2),n-8,ica,n-7,-1)
    else
      shf = 2
    end if
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(113,bo(23),an(3),dh(2),n-7+shf,n-8+shf,n-6+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(113,bo(23),an(3),dh(2),n-7+shf,n-9+shf,n-8+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(115,bo(23),an(3),dh(2),n-7+shf,n-9+shf,n-10+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(115,bo(23),an(3),dh(2),n-7+shf,n-9+shf,n-11+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(117,bo(23),an(3),dh(2),n-7+shf,n-9+shf,n-11+shf,0)
    end if
!
! tyrosine residue (TYR)
! 
  else if (resname .eq. 'TYR') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(124,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncg = n
    call z_add(126,bo(22),an(12),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd1 = n
    call z_add(127,bo(20),an(13),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
!    call z_add(127,bo(20),an(13),an(71),n-2,n-3,n-1,1)
    call z_add(127,bo(20),an(13),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(129,bo(18),an(14),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(129,bo(18),an(14),dh(1),n-2,n-4,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(131,bo(29),an(15),dh(1),n-2,n-4,n-5,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(132,bo(33),an(16),dh(2),n-1,n-2,n-4,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(125,bo(1),an(2),an(2),n-8,ica,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(125,bo(1),an(2),an(2),n-9,ica,n-8,-1)
    else
      shf = 2
    end if
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(128,bo(23),an(3),dh(2),n-8+shf,n-6+shf,n-4+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(128,bo(23),an(3),dh(2),n-8+shf,n-6+shf,n-5+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(130,bo(23),an(3),dh(2),n-8+shf,n-10+shf,n-11+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(130,bo(23),an(3),dh(2),n-8+shf,n-10+shf,n-12+shf,0)
    else
      shf = shf + 4
    end if
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(133,bo(34),an(17),chi(3,i),n-7+shf,n-8+shf,n-9+shf,0)
!
! deprotonated tyrosine residue (TYO, phenolate anion)
! 
  else if (resname .eq. 'TYO') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1131,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncg = n
    call z_add(1133,bo(22),an(12),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd1 = n
    call z_add(1134,bo(20),an(13),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
!    call z_add(127,bo(20),an(13),an(71),n-2,n-3,n-1,1)
    call z_add(1134,bo(20),an(13),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(1136,bo(18),an(14),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(1136,bo(18),an(14),dh(1),n-2,n-4,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(1138,bo(29),an(15),dh(1),n-2,n-4,n-5,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(1139,bo(63),an(16),dh(2),n-1,n-2,n-4,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1132,bo(1),an(2),an(2),n-8,ica,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1132,bo(1),an(2),an(2),n-9,ica,n-8,-1)
    else
      shf = 2
    end if
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(1135,bo(23),an(3),dh(2),n-8+shf,n-6+shf,n-4+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1135,bo(23),an(3),dh(2),n-8+shf,n-6+shf,n-5+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(1137,bo(23),an(3),dh(2),n-8+shf,n-10+shf,n-11+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1137,bo(23),an(3),dh(2),n-8+shf,n-10+shf,n-12+shf,0)
    end if
!
! phosphotyrosine residue (PTR) 
!
  else if (resname .eq. 'PTR') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1146,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncg = n
    call z_add(1148,bo(22),an(12),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd1 = n
    call z_add(1149,bo(20),an(13),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
    call z_add(1149,bo(20),an(13),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(1151,bo(18),an(14),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(1151,bo(18),an(14),dh(1),n-2,n-4,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(1153,bo(29),an(15),dh(1),n-2,n-4,n-5,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(1154,bo(59),an(16),dh(2),n-1,n-2,n-4,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1155,bo(56),an(91),chi(3,i),n-1,n-2,n-3,0)
    chiline(4,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1157,bo(57),an(86),chi(4,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1156,bo(58),an(88),an(89),n-2,n-1,n-3,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1156,bo(58),an(88),an(89),n-3,n-2,n-4,-1)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1147,bo(1),an(2),an(2),n-12,ica,n-11,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1147,bo(1),an(2),an(2),n-13,ica,n-12,-1)
    else
      shf = 2
    end if
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(1150,bo(23),an(3),dh(2),n-12+shf,n-10+shf,n-8+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1150,bo(23),an(3),dh(2),n-12+shf,n-10+shf,n-9+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(1152,bo(23),an(3),dh(2),n-12+shf,n-14+shf,n-15+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(1152,bo(23),an(3),dh(2),n-12+shf,n-14+shf,n-16+shf,0)
    else
      shf = shf + 4
    end if
    chiline(5,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1158,bo(35),an(78),chi(5,i),n-9+shf,n-10+shf,n-11+shf,0)
!
! tryptophan residue  (TRP) 
!
  else if (resname .eq. 'TRP') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(140,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(142,bo(11),an(18),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd1 = n
    call z_add(143,bo(12),an(19),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
!    call z_add(145,bo(13),an(20),an(21),n-2,n-3,n-1,1)
    call z_add(145,bo(13),an(20),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(146,bo(14),an(22),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(148,bo(15),an(23),dh(1),n-1,n-3,n-4,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-3,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(149,bo(16),an(24),dh(2),n-3,n-1,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(151,bo(17),an(25),dh(1),n-2,n-4,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(153,bo(18),an(26),dh(1),n-2,n-5,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(155,bo(19),an(27),dh(1),n-2,n-4,n-6,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(141,bo(1),an(2),an(2),n-10,ica,n-9,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(141,bo(1),an(2),an(2),n-11,ica,n-10,-1)
    else
      shf = 2
    end if
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd1
      call z_add(144,bo(23),an(28),dh(2),n-10+shf,n-11+shf,n-9+shf,0)
    else
      shf = shf + 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd1
    call z_add(147,bo(24),an(29),dh(2),n-9+shf,n-11+shf,n-12+shf,0)
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(150,bo(23),an(6),dh(2),n-8+shf,n-11+shf,n-9+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(152,bo(23),an(80),dh(2),n-8+shf,n-10+shf,n-12+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(154,bo(23),an(81),dh(2),n-8+shf,n-10+shf,n-13+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(156,bo(23),an(82),dh(2),n-8+shf,n-10+shf,n-12+shf,0)
    end if
!
! protonated histidine residue  (HIP)
!
! somehow the ring geometry isn't completely accurate in engh/huber
! adjusted bond length for CR1H-NH1 to 1.371 (from 1.374)
!
  else if (resname .eq. 'HIP') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(163,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(165,bo(25),an(11),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nnd1 = n
    call z_add(166,bo(29),an(43),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
!    call z_add(168,bo(30),an(44),an(45),n-2,n-3,n-1,1)
    call z_add(168,bo(30),an(44),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nnd1
    call z_add(170,bo(31),an(42),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(172,bo(27),an(47),dh(1),n-2,n-4,n-3,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(164,bo(1),an(2),an(2),n-6,ica,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(164,bo(1),an(2),an(2),n-7,ica,n-6,-1)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nnd1
    call z_add(167,bo(10),an(48),dh(1),n-6+shf,n-7+shf,n-8+shf,0)
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(169,bo(23),an(49),dh(2),n-6+shf,n-8+shf,n-7+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nnd1
      call z_add(171,bo(23),an(50),dh(2),n-6+shf,n-8+shf,n-9+shf,0)
    else
      shf = shf + 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(173,bo(10),an(51),dh(2),n-6+shf,n-8+shf,n-10+shf,0)
!
! histidine (HD only) residue  (HID)
!
! managed to build this with 'latest' parameters from X-PLOR, however, still incosistent
! with planarity (poor statistics for HID) -> only difference NH1-CRH-NR 109.5 instead of 111.7
! also note that according to Engh/Huber C-beta is out of the plane (i.e. C-gamma is not a perfect sp2)??
!
  else if (resname .eq. 'HID') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(180,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(182,bo(25),an(11),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nnd1 = n
    call z_add(183,bo(29),an(43),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
    call z_add(185,bo(26),an(52),an(53),n-2,n-3,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nnd1
    call z_add(187,bo(32),an(54),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(189,bo(18),an(55),dh(1),n-2,n-4,n-3,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(181,bo(1),an(2),an(2),n-6,ica,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(181,bo(1),an(2),an(2),n-7,ica,n-6,-1)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nnd1
    call z_add(184,bo(10),an(51),dh(2),n-6+shf,n-7+shf,n-5+shf,0)
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(186,bo(23),an(56),dh(2),n-6+shf,n-8+shf,n-7+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nnd1
      call z_add(188,bo(23),an(57),dh(2),n-6+shf,n-8+shf,n-9+shf,0)
    end if
!
! histidine (HE only) residue  (HIE)
!
! HIE is relatively well-defined in the complete Engh/Huber set (from X-PLOR)
! only minor problems, everything planar ...
!
  else if (resname .eq. 'HIE') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(196,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(198,bo(25),an(11),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nnd1 = n
    call z_add(199,bo(27),an(58),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncd2 = n
!    call z_add(200,bo(26),an(52),an(42),n-2,n-3,n-1,1)
    call z_add(200,bo(26),an(52),dh(2),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nnd1
    call z_add(202,bo(28),an(59),dh(1),n-2,n-3,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(204,bo(14),an(60),dh(1),n-2,n-4,n-3,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(197,bo(1),an(2),an(2),n-6,ica,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(197,bo(1),an(2),an(2),n-7,ica,n-6,-1)
    else
      shf = 2
    end if
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = ncd2
      call z_add(201,bo(23),an(61),dh(2),n-5+shf,n-7+shf,n-6+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nnd1
      call z_add(203,bo(23),an(62),dh(2),n-5+shf,n-7+shf,n-8+shf,0)
    else
      shf = shf + 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = ncd2
    call z_add(205,bo(10),an(63),dh(2),n-5+shf,n-7+shf,n-9+shf,0)
!
! aspartic acid residue  (ASP)
!
! not sure i found all necessary parameters in engh/huber
! upped CH2E-C-OC from 118.4 to 118.55 to make unit planar
!
  else if (resname .eq. 'ASP') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(212,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(214,bo(50),an(67),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(215,bo(48),an(4),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    call z_add(215,bo(48),an(4),an(5),n-2,n-3,n-1,1)
    call z_add(215,bo(48),an(4),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(213,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(213,bo(1),an(2),an(2),n-5,ica,n-4,-1)
    end if
!
! protonated (neutral) aspartic acid residue

  else if (resname .eq. 'ASH') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1094,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1096,bo(50),an(67),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1097,bo(60),an(28),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nnd2 = n
    call z_add(1098,bo(61),an(93),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1095,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1095,bo(1),an(2),an(2),n-5,ica,n-4,-1)
    else
      shf = 2
    end if
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1099,bo(62),an(8),chi(3,i),n-3+shf,n-5+shf,n-6+shf,0)
!
! asparagine residue  (ASN)
!
! the amide isn't planar when built according to engh/huber
! i increased the angles CH2E-C-O and NH2-C-O by 0.1 each to fix that
!
  else if (resname .eq. 'ASN') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(222,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(224,bo(46),an(68),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(225,bo(47),an(69),chi(2,i),n-1,n-2,ica,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nnd2 = n
!    call z_add(226,bo(49),an(70),an(43),n-2,n-3,n-1,1)
    call z_add(226,bo(49),an(70),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(223,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(223,bo(1),an(2),an(2),n-5,ica,n-4,-1)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(227,bo(10),an(3),dh(1),n-3+shf,n-5+shf,n-6+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(227,bo(10),an(3),dh(2),n-4+shf,n-6+shf,n-1,0)
!
! glutamic acid residue  (GLU)
!
! not sure i found all necessary parameters in engh/huber
! upped CH2E-C-OC from 118.4 to 118.55 to make unit planar
!
  else if (resname .eq. 'GLU') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(234,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(236,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(238,bo(46),an(68),chi(2,i),n-1,n-2,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(239,bo(48),an(4),chi(3,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    call z_add(239,bo(48),an(4),an(5),n-2,n-3,n-1,1)
    call z_add(239,bo(48),an(4),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(235,bo(1),an(2),an(2),n-5,ica,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(235,bo(1),an(2),an(2),n-6,ica,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(237,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(237,bo(1),an(2),an(2),n-7,n-8,n-6,-1)
    end if
!
! protonated glutamic acid residue (GLH)
!
  else if (resname .eq. 'GLH') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1080,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1082,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1084,bo(46),an(68),chi(2,i),n-1,n-2,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1085,bo(60),an(28),chi(3,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1086,bo(61),an(93),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1081,bo(1),an(2),an(2),n-5,ica,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1081,bo(1),an(2),an(2),n-6,ica,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1083,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1083,bo(1),an(2),an(2),n-7,n-8,n-6,-1) 
      shf = 0
    else
      shf = 4
    end if
    chiline(4,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1087,bo(62),an(86),chi(4,i),n-5+shf,n-7+shf,n-8+shf,0)
!
! glutamine residue  (GLN)
!
! the amide isn't planar when built according to engh/huber
! i increased the angles CH2E-C-O and NH2-C-O by 0.1 each to fix that
!
  else if (resname .eq. 'GLN') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(246,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(248,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(250,bo(46),an(68),chi(2,i),n-1,n-2,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(251,bo(47),an(69),chi(3,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nne2 = n
!    call z_add(252,bo(49),an(70),an(43),n-2,n-3,n-1,1)
    call z_add(252,bo(49),an(70),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(247,bo(1),an(2),an(2),n-5,ica,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(247,bo(1),an(2),an(2),n-6,ica,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(249,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(249,bo(1),an(2),an(2),n-7,n-8,n-6,-1)
    else
      shf = 4
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(253,bo(10),an(3),dh(1),n-5+shf,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(253,bo(10),an(3),dh(2),n-6+shf,n-8+shf,n-1,0)
!
! methionnine residue  (MET)
!
  else if (resname .eq. 'MET') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(260,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(262,bo(2),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(264,bo(45),an(77),chi(2,i),n-1,n-2,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(265,bo(44),an(78),chi(3,i),n-1,n-2,n-3,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(261,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(261,bo(1),an(2),an(2),n-5,ica,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(263,bo(1),an(2),an(2),n-5,n-6,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(263,bo(1),an(2),an(2),n-6,n-7,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(266,bo(1),an(2),dh(2),n-5,n-6,n-7,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(266,bo(1),an(2),an(2),n-6,n-7,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(266,bo(1),an(2),an(2),n-7,n-8,n-2,-1)
    end if
!
! lysine residue  (LYS)
!
  else if (resname .eq. 'LYS') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(273,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(275,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(277,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(279,bo(4),an(65),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(281,bo(9),an(66),chi(4,i),n-1,n-2,n-3,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(274,bo(1),an(2),an(2),n-5,ica,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(274,bo(1),an(2),an(2),n-6,ica,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(276,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(276,bo(1),an(2),an(2),n-7,n-8,n-6,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(278,bo(1),an(2),an(2),n-7,n-8,n-6,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(278,bo(1),an(2),an(2),n-8,n-9,n-7,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(280,bo(1),an(2),an(2),n-8,n-9,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(280,bo(1),an(2),an(2),n-9,n-10,n-8,-1)
    else
      shf = 8
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(282,bo(10),an(2),dh(2),n-9+shf,n-10+shf,n-11+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(282,bo(10),an(2),an(2),n-10+shf,n-11+shf,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(282,bo(10),an(2),an(2),n-11+shf,n-12+shf,n-2,-1)
!
!  deprotonated lysine residue  (LYD)
!
  else if (resname .eq. 'LYD') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1106,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1108,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1110,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1112,bo(4),an(65),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1114,bo(9),an(66),chi(4,i),n-1,n-2,n-3,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1107,bo(1),an(2),an(2),n-5,ica,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1107,bo(1),an(2),an(2),n-6,ica,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1109,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1109,bo(1),an(2),an(2),n-7,n-8,n-6,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1111,bo(1),an(2),an(2),n-7,n-8,n-6,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1111,bo(1),an(2),an(2),n-8,n-9,n-7,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1113,bo(1),an(2),an(2),n-8,n-9,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1113,bo(1),an(2),an(2),n-9,n-10,n-8,-1)
    else
      shf = 8
    end if
    nchi(i) = nchi(i) + 1
    chiline(5,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1115,bo(10),an(84),chi(5,i),n-9+shf,n-10+shf,n-11+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1115,bo(10),an(84),an(85),n-10+shf,n-11+shf,n-1,1)
!
!  acetylated lysine residue  (KAC)
!
  else if (resname .eq. 'KAC') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1193,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1195,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1197,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1199,bo(4),an(65),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1201,bo(68),an(66),chi(4,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(5,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1203,bo(69),an(99),chi(5,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(6,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1205,bo(70),an(3),chi(6,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1204,bo(47),an(3),dh(2),n-2,n-1,n-3,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1194,bo(1),an(2),an(2),n-8,ica,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1194,bo(1),an(2),an(2),n-9,ica,n-8,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1196,bo(1),an(2),an(2),n-9,n-10,n-8,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1196,bo(1),an(2),an(2),n-10,n-11,n-9,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1198,bo(1),an(2),an(2),n-10,n-11,n-9,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1198,bo(1),an(2),an(2),n-11,n-12,n-10,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1200,bo(1),an(2),an(2),n-11,n-12,n-10,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1200,bo(1),an(2),an(2),n-12,n-13,n-11,-1)
    else
      shf = 8
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1202,bo(10),an(100),dh(2),n-12+shf,n-11+shf,n-13+shf,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1206,bo(1),an(2),dh(2),n-11+shf,n-12+shf,n-13+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1206,bo(1),an(2),an(2),n-12+shf,n-13+shf,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1206,bo(1),an(2),an(2),n-13+shf,n-14+shf,n-2,-1)
    end if
!
!  monomethylated lysine residue  (KM1)
!
  else if (resname .eq. 'KM1') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1213,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1215,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1217,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1219,bo(4),an(65),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1221,bo(65),an(66),chi(4,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(5,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1223,bo(65),an(94),chi(5,i),n-1,n-2,n-3,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1214,bo(1),an(2),an(2),n-6,ica,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1214,bo(1),an(2),an(2),n-7,ica,n-6,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1216,bo(1),an(2),an(2),n-7,n-8,n-6,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1216,bo(1),an(2),an(2),n-8,n-9,n-7,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1218,bo(1),an(2),an(2),n-8,n-9,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1218,bo(1),an(2),an(2),n-9,n-10,n-8,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1220,bo(1),an(2),an(2),n-9,n-10,n-8,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1220,bo(1),an(2),an(2),n-10,n-11,n-9,-1)
    else
      shf = 8
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1222,bo(10),an(95),an(95),n-10+shf,n-11+shf,n-9+shf,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1222,bo(10),an(95),an(95),n-11+shf,n-12+shf,n-10+shf,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1224,bo(1),an(2),dh(1),n-11+shf,n-12+shf,n-13+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1224,bo(1),an(2),an(2),n-12+shf,n-13+shf,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1224,bo(1),an(2),an(2),n-13+shf,n-14+shf,n-2,-1)
    end if
!
!  dimethylated lysine residue  (KM2)
!
  else if (resname .eq. 'KM2') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1231,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1233,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1235,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1237,bo(4),an(65),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1239,bo(66),an(66),chi(4,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(5,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1241,bo(66),an(97),chi(5,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1241,bo(66),an(97),an(97),n-2,n-3,n-1,1)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1232,bo(1),an(2),an(2),n-7,ica,n-6,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1232,bo(1),an(2),an(2),n-8,ica,n-7,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1234,bo(1),an(2),an(2),n-8,n-9,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1234,bo(1),an(2),an(2),n-9,n-10,n-8,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1236,bo(1),an(2),an(2),n-9,n-10,n-8,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1236,bo(1),an(2),an(2),n-10,n-11,n-9,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1238,bo(1),an(2),an(2),n-10,n-11,n-9,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1238,bo(1),an(2),an(2),n-11,n-12,n-10,-1)
    else
      shf = 8
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1240,bo(10),an(98),an(98),n-11+shf,n-12+shf,n-10+shf,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1242,bo(1),an(2),dh(1),n-11+shf,n-12+shf,n-13+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1242,bo(1),an(2),an(2),n-12+shf,n-13+shf,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1242,bo(1),an(2),an(2),n-13+shf,n-14+shf,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1242,bo(1),an(2),dh(1),n-13+shf,n-15+shf,n-16+shf,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1242,bo(1),an(2),an(2),n-14+shf,n-16+shf,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1242,bo(1),an(2),an(2),n-15+shf,n-17+shf,n-2,-1)
    end if
!
!  trimethylated lysine residue  (KM3)
!
  else if (resname .eq. 'KM3') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1249,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1251,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1253,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1255,bo(4),an(65),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1257,bo(67),an(66),chi(4,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(5,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1258,bo(67),an(2),chi(5,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1258,bo(67),an(2),an(2),n-2,n-3,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(1258,bo(67),an(2),an(2),n-3,n-4,n-2,-1)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1250,bo(1),an(2),an(2),n-8,ica,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1250,bo(1),an(2),an(2),n-9,ica,n-8,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1252,bo(1),an(2),an(2),n-9,n-10,n-8,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1252,bo(1),an(2),an(2),n-10,n-11,n-9,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1254,bo(1),an(2),an(2),n-10,n-11,n-9,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1254,bo(1),an(2),an(2),n-11,n-12,n-10,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1256,bo(1),an(2),an(2),n-11,n-12,n-10,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1256,bo(1),an(2),an(2),n-12,n-13,n-11,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),dh(1),n-11,n-12,n-13,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),an(2),n-12,n-13,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),an(2),n-13,n-14,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),dh(1),n-13,n-15,n-16,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),an(2),n-14,n-16,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),an(2),n-15,n-17,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),dh(1),n-15,n-18,n-19,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),an(2),n-16,n-19,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(1259,bo(1),an(2),an(2),n-17,n-20,n-2,-1)
    end if
!
! arginine residue  (ARG)
!
! upped the NH1-C-NC2 by 0.15 to accomodate the NC2-C-NC2 of 119.7
!
  else if (resname .eq. 'ARG') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(289,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(291,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(293,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    nchi(i) = nchi(i) + 1
    chiline(3,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(295,bo(43),an(73),chi(3,i),n-1,n-2,n-3,0)
    nchi(i) = nchi(i) + 1
    chiline(4,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    ncz = n
    call z_add(297,bo(42),an(74),chi(4,i),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(298,bo(41),an(75),dh(2),n-1,n-2,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    call z_add(298,bo(41),an(75),an(76),n-2,n-3,n-1,1)
    call z_add(298,bo(41),an(75),dh(2),n-2,n-3,n-1,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(290,bo(1),an(2),an(2),n-7,ica,n-6,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(290,bo(1),an(2),an(2),n-8,ica,n-7,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(292,bo(1),an(2),an(2),n-8,n-9,n-7,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(292,bo(1),an(2),an(2),n-9,n-10,n-8,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(294,bo(1),an(2),an(2),n-9,n-10,n-8,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(294,bo(1),an(2),an(2),n-10,n-11,n-9,-1)
    else
      shf = 6
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!    call z_add(296,bo(10),an(3),an(3),n-10,n-11,n-9,1)
    call z_add(296,bo(10),an(3),dh(2),n-10+shf,n-11+shf,n-9+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(299,bo(10),an(3),dh(2),n-9+shf,n-10+shf,n-11+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(299,bo(10),an(3),dh(2),n-10+shf,n-11+shf,n-1,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(299,bo(10),an(3),dh(2),n-10+shf,n-12+shf,n-13+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(299,bo(10),an(3),dh(2),n-11+shf,n-13+shf,n-1,0)
!
! gamma residue  (GAM)
!
  else if (resname .eq. 'GAM') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(391,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(393,bo(6),an(11),chi(1,i),n-1,ica,nn,0)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(392,bo(1),an(2),an(2),n-2,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(392,bo(1),an(2),an(2),n-3,ica,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(394,bo(1),an(2),dh(2),n-3,n-4,ica,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(394,bo(1),an(2),an(2),n-4,n-5,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(394,bo(1),an(2),an(2),n-5,n-6,n-2,-1)
    end if
!
! ornithine residue  (ORN)
!
  else if (resname .eq. 'ORN') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(306,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(308,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(310,bo(4),an(65),chi(2,i),n-1,n-2,ica,0)
    chiline(3,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(312,bo(9),an(66),chi(3,i),n-1,n-2,n-3,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(307,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(307,bo(1),an(2),an(2),n-5,ica,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(309,bo(1),an(2),an(2),n-5,n-6,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(309,bo(1),an(2),an(2),n-6,n-7,n-5,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(311,bo(1),an(2),an(2),n-6,n-7,n-5,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(311,bo(1),an(2),an(2),n-7,n-8,n-6,-1)
    else
      shf = 6
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(313,bo(10),an(2),dh(2),n-7+shf,n-8+shf,n-9+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(313,bo(10),an(2),an(2),n-8+shf,n-9+shf,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(313,bo(10),an(2),an(2),n-9+shf,n-10+shf,n-2,-1)
!
! alpha,gamma-diamino butyric acid residue
!
  else if (resname .eq. 'DAB') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(525,bo(2),an(8),an(9),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(527,bo(4),an(64),chi(1,i),n-1,ica,nn,0)
    chiline(2,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(529,bo(9),an(66),chi(2,i),n-1,n-2,ica,0)
    if (ua_model.eq.0) then
      shf = 0    
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(526,bo(1),an(2),an(2),n-3,ica,n-2,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(526,bo(1),an(2),an(2),n-4,ica,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(528,bo(1),an(2),an(2),n-4,n-5,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(528,bo(1),an(2),an(2),n-5,n-6,n-4,-1)
    else
      shf = 4
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(530,bo(10),an(2),dh(2),n-5+shf,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(530,bo(10),an(2),an(2),n-6+shf,n-7+shf,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
   call z_add(530,bo(10),an(2),an(2),n-7+shf,n-8+shf,n-2,-1)
!
! methylalannne residue (AIB)
!
 else if (resname .eq. 'AIB') then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(319,bo(3),an(7),an(67),ica,nn,cc,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(319,bo(3),an(7),an(67),ica,nn,cc,-1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(320,bo(1),an(2),chi(1,i),n-2,ica,nn,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(320,bo(1),an(2),an(2),n-3,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(320,bo(1),an(2),an(2),n-4,ica,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(320,bo(1),an(2),chi(1,i),n-4,ica,nn,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(320,bo(1),an(2),an(2),n-5,ica,n-1,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(320,bo(1),an(2),an(2),n-6,ica,n-2,-1)
    end if
!
! pyroglutamic accd residue  (PCA)
!
  else if (resname .eq. 'PCA') then
    if (fline(i).gt.0) pucline(i) = fline(i)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(327,bo(2),an(10),an(46),ica,nn,cc,chiral(i))
    chiline(1,i) = n
    nchi(i) = nchi(i) + 1
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(329,bo(51),an(10),chi(1,i),n-1,ica,nn,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(331,bo(7),an(10),chi(2,i),n-1,n-2,ica,0)
!    at(i)%nsc = at(i)%nsc + 1
!    at(i)%sc(at(i)%nsc) = n
!    atmres(n) = i
!   just a ring closure
    call z_add(-1,bo(52),an(1),dh(1),nn,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(332,bo(47),an(79),an(79),n-1,nn,n-2,1)
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(328,bo(1),an(2),an(2),n-4,ica,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(328,bo(1),an(2),an(2),n-5,ica,n-4,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(330,bo(1),an(2),an(2),n-5,n-6,n-4,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(330,bo(1),an(2),an(2),n-6,n-7,n-5,-1)
    end if
!
! unknown residue  (UNK)
!
  else if (resname .eq. 'UNK') then
    if (ua_model.eq.0) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(6,bo(1),an(30),an(31),ica,nn,cc,chiral(i))
    end if
!
  else
!
    write(ilog,*) 'Fatal. Called sidechain(...) with unsupported res&
 &idue type. This is either a bug or a feature which is not yet supp&
 &orted (offending residue: ',resname,').'
    call fexit()
!
  end if
!
end
!
!------------------------------------------------------------------------
!
subroutine nuc_sidechain(resname,i,o5pi,c5i,c4i,c3i)
!
  use atoms
  use fyoc
  use polypep
  use sequen
  use molecule
  use iounit
  use system
!
  implicit none
!
  integer i,o5pi,c5i,c4i,c3i,shf
  integer nn1,nc1,nc2
  character(3) resname
! some lengths, angles and dihedrals
  RTYPE dh(100),an(100),bo(100)
!
! see proteus.f for our strategy with regards to eliminating calls of zatom
! with the numbers provided directly as arguments
!
! now the nucleic acid parameters (Parkinson et al. in Acta Cryst. D52 (1996) p57-64)
! nevermind the bad numbering ...
!
  bo(1) = 0.0d0
  bo(2) = 1.091d0
  bo(3) = 1.01d0
  bo(4) = 0.945d0
  bo(11) = 1.413d0   ! C2-O2 in R
  bo(12) = 1.525d0   ! C3-C2 in R
  bo(13) = 1.528d0   ! C2-C1 in R
  bo(14) = 1.453d0   ! C4-O4 in R
  bo(15) = 1.446d0   ! C4-O4 in D
  bo(16) = 1.518d0   ! C3-C2 in D
  bo(17) = 1.521d0   ! C2-C1 in D
  bo(18) = 1.420d0   ! C1-O1 in D/R (made up)
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
  bo(42) = 1.380d0   ! N3-C4 in URA
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
  bo(58) = 1.350d0   ! C4-N3 in GUA
  bo(59) = 1.323d0   ! N3-C2 in GUA
  bo(60) = 1.341d0   ! C2-N2 in GUA
  bo(61) = 1.237d0   ! C6-O6 in GUA
  bo(62) = 1.419d0   ! C5-C6 in GUA
  bo(63) = 1.459d0   ! C1*-N9 in R/D-GUA
  bo(64) = 1.462d0   ! C1*-N9 in R/D-ADE
  bo(65) = 1.473d0   ! C1*-N9 in R/D-THY
  bo(66) = 1.469d0   ! C1*-N9 in R/D-URA
  bo(67) = 1.470d0   ! C1*-N9 in R/D-CYT
!
!
  an(1) = 0.0d0
  an(2) = 120.0d0
  an(3) = 109.47d0
  an(4) = 108.5d0
  an(11) = 102.7d0   ! C4C3C2 in R and C3C2C1 in D
  an(12) = 101.5d0   ! C3C2C1 in R
  an(13) = 105.5d0   ! C3C4O4 in R
  an(14) = 113.3d0   ! C3C2O2 in R
  an(15) = 108.2d0   ! O4C1N1/9 in R
  an(16) = 103.2d0   ! C4C3C2 in D
  an(17) = 105.6d0   ! C3C4O4 in D
  an(18) = 107.8d0   ! O4C1N1/9 in D
  an(19) = 113.1d0   ! C2C1N1/9 in R
  an(20) = 115.0d0   ! C2C1N1/9 in D
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
  an(53) = 127.0d0   ! C2N3C4 in URA and C1*N9C8 in ADE
  an(54) = 122.8d0   ! N1C2O2 in URA
  an(55) = 122.7d0   ! N1C6C5 in URA
  an(56) = 119.4d0   ! N3C4O4 in URA
  an(57) = 119.5d0   ! for H1 in URA, for H6 in CYT
  an(58) = 118.65d0  ! for H6 in URA
  an(59) = 116.5d0   ! for H3 in URA
  an(60) = 120.15d0  ! for H5 in URA
  an(61) = 121.3d0   ! C6N1C2 in THY, for H5 in CYT
  an(62) = 114.6d0   ! N1C2N3 in THY
  an(63) = 127.2d0   ! C2N3C4 in THY
  an(64) = 123.1d0   ! N1C2O2 in THY, for H8 in ADE
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
  an(81) = 128.6d0   ! C5C6O6 in GUA
  an(82) = 105.4d0   ! N9C4C5 in GUA
  an(83) = 113.1d0   ! N9C8N7 in GUA
  an(84) = 126.0d0   ! N9C4N3 in GUA
  an(85) = 111.9d0   ! C4N3C2 in GUA
  an(86) = 118.8d0   ! C4C5C6 in GUA and C1*N1C2 in CYT
  an(87) = 123.9d0   ! N3C2N1 in GUA
  an(88) = 126.8d0   ! for H9 in GUA
  an(89) = 123.45d0  ! for H8 in GUA
  an(90) = 117.45d0  ! for H1 in GUA
  an(91) = 120.8d0   ! C1*N1C6 in CYT
  an(92) = 120.4d0   ! C1*N1C6 in THY
  an(93) = 118.2d0   ! C1*N1C2 in THY
  an(94) = 121.2d0   ! C1*N1C6 in URA
  an(95) = 117.7d0   ! C1*N1C2 in URA
  an(96) = 127.7d0   ! C1*N9C8 in ADE
  an(97) = 126.3d0   ! C1*N9C4 in ADE
  an(98) = 126.5d0   ! C1*N9C4 in GUA
!
!
  dh(1) = 0.0d0
  dh(3) = 180.0d0
  dh(11) = -96.6d0   ! C2C3C4C5 in C2'-endo-R
  dh(12) = -35.4d0   ! C1C2C3C4 in C2'-endo-R
  dh(13) = 24.2d0    ! C2C3C4O4 in C2'-endo-R
  dh(14) = -157.6d0  ! O2C2C3C4 in C2'-endo-R
  dh(15) = -142.3d0  ! C3C2C1N1/9 in C2'-endo-R/D
  dh(16) = -98.0d0   ! C2C3C4C5 in C2'-endo-D
  dh(17) = -33.1d0   ! C1C2C3C4 in C2'-endo-D
  dh(18) = 22.6d0    ! C2C3C4O4 in C2'-endo-D
!
! empty ribonucleotide
!
  if ((resname.eq.'R5P').OR.(resname.eq.'RIB')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(637,bo(12),an(11),dh(11),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(641,bo(13),an(12),dh(12),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(631,bo(14),an(13),dh(13),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(638,bo(11),an(14),an(14),n-3,c3i,n-2,1)
!    call z_add(638,bo(11),an(14),dh(14),n-3,c3i,c4i,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(639,bo(2),an(3),an(3),n-4,c3i,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(643,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 2
    end if
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(640,bo(4),an(4),chi(1,i),n-3+shf,n-6+shf,c3i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!   note the bond length is made up
    call z_add(642,bo(18),an(15),an(15),n-6+shf,n-7+shf,n-5+shf,-1)
!    call z_add(642,bo(18),an(15),dh(15),n-6+shf,n-5+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(644,bo(4),an(4),chi(2,i),n-1,n-7+shf,n-8+shf,0)
!
! empty deoxyribonucleotide
!
  else if ((resname .eq. 'D5P').OR.(resname.eq.'DIB')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(949,bo(16),an(16),dh(16),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(951,bo(17),an(11),dh(17),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(943,bo(15),an(17),dh(18),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(950,bo(2),an(3),an(3),n-3,c3i,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(950,bo(2),an(3),an(3),n-4,c3i,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(953,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 3
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
!   note the bond length is made up
    call z_add(952,bo(18),an(18),an(18),n-5+shf,n-6+shf,n-4+shf,-1)
!    call z_add(952,bo(18),an(18),dh(15),n-5+shf,n-4+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(954,bo(4),an(4),chi(1,i),n-1,n-6+shf,n-7+shf,0)
!
! cytidine (RPC)
!
  else if ((resname .eq. 'RPC').OR.(resname.eq.'RIC')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(682,bo(12),an(11),dh(11),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(685,bo(13),an(12),dh(12),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(678,bo(14),an(13),dh(13),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(683,bo(11),an(14),an(14),n-3,c3i,n-2,1)
!    call z_add(683,bo(11),an(14),dh(14),n-3,c3i,c4i,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(684,bo(2),an(3),an(3),n-4,c3i,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(686,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 2
    end if
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(931,bo(4),an(4),chi(1,i),n-3+shf,n-6+shf,c3i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(687,bo(67),an(19),an(15),n-6+shf,n-7+shf,n-5+shf,-1)
!    call z_add(687,bo(67),an(15),dh(15),n-6+shf,n-5+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(688,bo(31),an(86),chi(2,i),n-1,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(696,bo(33),an(91),an(42),n-2,n-8+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(690,bo(34),an(44),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(689,bo(32),an(43),dh(3),n-3,n-4,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(694,bo(35),an(45),dh(1),n-3,n-5,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(691,bo(36),an(66),dh(1),n-3,n-5,n-6,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(692,bo(36),an(46),dh(3),n-1,n-4,n-6,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(697,bo(2),an(57),dh(3),n-6,n-8,n-7,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(695,bo(2),an(61),dh(3),n-4,n-7,n-9,0)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(693,bo(2),an(2),dh(1),n-3+shf,n-4+shf,n-6+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(693,bo(2),an(2),dh(3),n-4+shf,n-5+shf,n-1,0)
!
! deoxy-cytidine (DPC)
!
  else if ((resname .eq. 'DPC').OR.(resname.eq.'DIC')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(656,bo(16),an(16),dh(16),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(658,bo(17),an(11),dh(17),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(652,bo(15),an(17),dh(18),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(657,bo(2),an(3),an(3),n-3,c3i,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(657,bo(2),an(3),an(3),n-4,c3i,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(659,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 3
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(660,bo(67),an(20),an(18),n-5+shf,n-6+shf,n-4+shf,-1)
!    call z_add(660,bo(67),an(18),dh(15),n-5+shf,n-4+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(661,bo(31),an(86),chi(1,i),n-1,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(669,bo(33),an(91),an(42),n-2,n-7+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(663,bo(34),an(44),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(662,bo(32),an(43),dh(3),n-3,n-4,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(667,bo(35),an(45),dh(1),n-3,n-5,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(664,bo(36),an(66),dh(1),n-3,n-5,n-6,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(665,bo(36),an(46),dh(3),n-1,n-4,n-6,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(670,bo(2),an(57),dh(3),n-6,n-8,n-7,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(668,bo(2),an(61),dh(3),n-4,n-7,n-9,0)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(666,bo(2),an(2),dh(1),n-3+shf,n-4+shf,n-6+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(666,bo(2),an(2),dh(3),n-4+shf,n-5+shf,n-1,0)
!
! uridine
!
  else if ((resname .eq. 'RPU').OR.(resname.eq.'RIU')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(735,bo(12),an(11),dh(11),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(738,bo(13),an(12),dh(12),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(731,bo(14),an(13),dh(13),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(736,bo(11),an(14),an(14),n-3,c3i,n-2,1)
!    call z_add(736,bo(11),an(14),dh(14),n-3,c3i,c4i,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(737,bo(2),an(3),an(3),n-4,c3i,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(739,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 2
    end if
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(932,bo(4),an(4),chi(1,i),n-3+shf,n-6+shf,c3i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(740,bo(66),an(19),an(15),n-6+shf,n-7+shf,n-5+shf,-1)
!    call z_add(740,bo(66),an(15),dh(15),n-6+shf,n-5+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(741,bo(37),an(95),chi(2,i),n-1,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(749,bo(39),an(94),an(51),n-2,n-8+shf,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(743,bo(40),an(52),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(742,bo(38),an(54),dh(3),n-3,n-4,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(747,bo(41),an(55),dh(1),n-3,n-5,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(745,bo(42),an(53),dh(1),n-3,n-5,n-6,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(746,bo(43),an(56),dh(3),n-1,n-4,n-6,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(750,bo(2),an(58),dh(3),n-6,n-8,n-7,0)
    else
      shf = 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(744,bo(3),an(59),dh(3),n-6+shf,n-8+shf,n-9+shf,0)
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(748,bo(2),an(60),dh(3),n-5,n-4,n-8,0)
    end if
!
! deoxy-uridine (DPU)
!
  else if ((resname .eq. 'DPU').OR.(resname.eq.'DIU')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(709,bo(16),an(16),dh(16),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(711,bo(17),an(11),dh(17),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(705,bo(15),an(17),dh(18),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(710,bo(2),an(3),an(3),n-3,c3i,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(710,bo(2),an(3),an(3),n-4,c3i,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(712,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 3
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(713,bo(66),an(20),an(18),n-5+shf,n-6+shf,n-4+shf,-1)
!    call z_add(713,bo(66),an(18),dh(15),n-5+shf,n-4+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(714,bo(37),an(95),chi(1,i),n-1,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(722,bo(39),an(94),an(51),n-2,n-7+shf,n-1,1)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(716,bo(40),an(52),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(715,bo(38),an(54),dh(3),n-3,n-4,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(720,bo(41),an(55),dh(1),n-3,n-5,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(718,bo(42),an(53),dh(1),n-3,n-5,n-6,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(719,bo(43),an(56),dh(3),n-1,n-4,n-6,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(723,bo(2),an(58),dh(3),n-6,n-8,n-7,0)
    else
      shf = 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(717,bo(3),an(59),dh(3),n-6+shf,n-8+shf,n-9+shf,0)
    if (ua_model.lt.2) then
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(721,bo(2),an(60),dh(3),n-5,n-4,n-8,0)
    end if
!
! thymidine
!
  else if ((resname .eq. 'RPT').OR.(resname.eq.'RIT')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(789,bo(12),an(11),dh(11),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(792,bo(13),an(12),dh(12),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(785,bo(14),an(13),dh(13),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(790,bo(11),an(14),an(14),n-3,c3i,n-2,1)
!    call z_add(790,bo(11),an(14),dh(14),n-3,c3i,c4i,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(791,bo(2),an(3),an(3),n-4,c3i,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(793,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 2
    end if
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(933,bo(4),an(4),chi(1,i),n-3+shf,n-6+shf,c3i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(794,bo(65),an(19),an(15),n-6+shf,n-7+shf,n-5+shf,-1)
!    call z_add(794,bo(65),an(15),dh(15),n-6+shf,n-5+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(795,bo(44),an(93),chi(2,i),n-1,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(804,bo(46),an(92),an(61),n-2,n-8+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(797,bo(40),an(62),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(796,bo(45),an(64),dh(3),n-3,n-4,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(801,bo(35),an(65),dh(1),n-3,n-5,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(799,bo(47),an(63),dh(1),n-3,n-5,n-6,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(802,bo(49),an(67),dh(3),n-2,n-1,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(800,bo(48),an(66),dh(3),n-2,n-5,n-7,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(805,bo(2),an(69),dh(3),n-7,n-9,n-8,0)
    else
      shf = 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(798,bo(2),an(70),dh(3),n-7+shf,n-9+shf,n-10+shf,0)
    if (ua_model.eq.0) then
      chiline(3,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(803,bo(2),an(3),chi(3,i),n-4,n-6,n-5,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(803,bo(2),an(3),an(3),n-5,n-7,n-1,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(803,bo(2),an(3),an(3),n-6,n-8,n-2,1)
    end if
!
! deoxy-thymidine (DPT)
!
  else if ((resname .eq. 'DPT').OR.(resname.eq.'DIT')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(762,bo(16),an(16),dh(16),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(764,bo(17),an(11),dh(17),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(758,bo(15),an(17),dh(18),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(763,bo(2),an(3),an(3),n-3,c3i,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(763,bo(2),an(3),an(3),n-4,c3i,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(765,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 3
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(766,bo(65),an(20),an(18),n-5+shf,n-6+shf,n-4+shf,-1)
!    call z_add(766,bo(65),an(18),dh(15),n-5+shf,n-4+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(767,bo(44),an(93),chi(1,i),n-1,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(776,bo(46),an(92),an(61),n-2,n-7+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(769,bo(40),an(62),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(768,bo(45),an(64),dh(3),n-3,n-4,n-2,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(773,bo(35),an(65),dh(1),n-3,n-5,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(771,bo(47),an(63),dh(1),n-3,n-5,n-6,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(774,bo(49),an(67),dh(3),n-2,n-1,n-4,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nn1
    call z_add(772,bo(48),an(66),dh(3),n-2,n-5,n-7,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(777,bo(2),an(69),dh(3),n-7,n-9,n-8,0)
    else
      shf = 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(770,bo(2),an(70),dh(3),n-7+shf,n-9+shf,n-10+shf,0)
    if (ua_model.eq.0) then
      chiline(2,i) = n
      nchi(i) = nchi(i) + 1
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(775,bo(2),an(3),chi(2,i),n-4,n-6,n-5,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(775,bo(2),an(3),an(3),n-5,n-7,n-1,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(775,bo(2),an(3),an(3),n-6,n-8,n-2,1)
    end if
!
! adenosine
!
  else if ((resname .eq. 'RPA').OR.(resname.eq.'RIA')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(845,bo(12),an(11),dh(11),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(848,bo(13),an(12),dh(12),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(841,bo(14),an(13),dh(13),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(846,bo(11),an(14),an(14),n-3,c3i,n-2,1)
!    call z_add(846,bo(11),an(14),dh(14),n-3,c3i,c4i,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(847,bo(2),an(3),an(3),n-4,c3i,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(849,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 2
    end if
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(934,bo(4),an(4),chi(1,i),n-3+shf,n-6+shf,c3i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(862,bo(64),an(19),an(15),n-6+shf,n-7+shf,n-5+shf,-1)
!    call z_add(862,bo(64),an(15),dh(15),n-6+shf,n-5+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(854,bo(50),an(97),chi(2,i),n-1,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(860,bo(40),an(96),an(72),n-2,n-8+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(855,bo(51),an(72),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(859,bo(52),an(73),dh(1),n-2,n-4,n-3,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(853,bo(53),an(74),dh(3),n-4,n-5,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(851,bo(54),an(75),dh(3),n-1,n-5,n-6,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(856,bo(55),an(76),dh(3),n-4,n-6,n-7,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(850,bo(35),an(77),dh(1),n-2,n-3,n-7,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(857,bo(36),an(71),dh(3),n-2,n-1,n-3,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(861,bo(2),an(64),dh(3),n-8,n-10,n-9,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc1
      call z_add(852,bo(2),an(78),dh(3),n-5,n-6,n-10,0)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(858,bo(3),an(2),dh(3),n-3+shf,n-5+shf,n-4+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(858,bo(3),an(2),dh(3),n-4+shf,n-6+shf,n-1,0)
!
! deoxy-adenosine (DPA)
!
  else if ((resname .eq. 'DPA').OR.(resname.eq.'DIA')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(817,bo(16),an(16),dh(16),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(819,bo(17),an(11),dh(17),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(813,bo(15),an(17),dh(18),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(818,bo(2),an(3),an(3),n-3,c3i,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(818,bo(2),an(3),an(3),n-4,c3i,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(820,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 3
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(833,bo(64),an(20),an(18),n-5+shf,n-6+shf,n-4+shf,-1)
!    call z_add(833,bo(64),an(18),dh(15),n-5+shf,n-4+shf,c4i,0)
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    call z_add(825,bo(50),an(97),chi(1,i),n-1,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(831,bo(40),an(96),an(72),n-2,n-7+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(826,bo(51),an(72),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(830,bo(52),an(73),dh(1),n-2,n-4,n-3,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(824,bo(53),an(74),dh(3),n-4,n-5,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(822,bo(54),an(75),dh(3),n-1,n-5,n-6,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(827,bo(55),an(76),dh(3),n-4,n-6,n-7,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(821,bo(35),an(77),dh(1),n-2,n-3,n-7,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(828,bo(36),an(71),dh(3),n-2,n-1,n-3,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc2
      call z_add(832,bo(2),an(64),dh(3),n-8,n-10,n-9,0)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc1
      call z_add(823,bo(2),an(78),dh(3),n-5,n-6,n-10,0)
    else
      shf = 2
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(829,bo(3),an(2),dh(3),n-3+shf,n-5+shf,n-4+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(829,bo(3),an(2),dh(3),n-4+shf,n-6+shf,n-1,0)
!
! guanosine (RPG)
!
  else if ((resname .eq. 'RPG').OR.(resname.eq.'RIG')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(903,bo(12),an(11),dh(11),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(906,bo(13),an(12),dh(12),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(899,bo(14),an(13),dh(13),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(904,bo(11),an(14),an(14),n-3,c3i,n-2,1)
!    call z_add(904,bo(11),an(14),dh(14),n-3,c3i,c4i,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(905,bo(2),an(3),an(3),n-4,c3i,n-3,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(907,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 2
    end if
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(935,bo(4),an(4),chi(1,i),n-3+shf,n-6+shf,c3i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(921,bo(63),an(19),an(15),n-6+shf,n-7+shf,n-5+shf,-1)
!    call z_add(921,bo(63),an(15),dh(15),n-6+shf,n-5+shf,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    nchi(i) = nchi(i) + 1
    chiline(2,i) = n
    call z_add(914,bo(39),an(98),chi(2,i),n-1,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(919,bo(50),an(53),an(80),n-2,n-8+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(915,bo(56),an(82),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(918,bo(57),an(83),dh(1),n-2,n-4,n-3,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(913,bo(58),an(84),dh(3),n-4,n-5,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(910,bo(59),an(85),dh(3),n-1,n-5,n-6,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(916,bo(62),an(86),dh(3),n-4,n-6,n-7,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(908,bo(40),an(87),dh(1),n-2,n-3,n-7,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(911,bo(60),an(66),dh(3),n-3,n-4,n-8,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(917,bo(61),an(81),dh(1),n-3,n-7,n-6,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc1
      call z_add(920,bo(2),an(89),dh(3),n-9,n-7,n-11,0)
    else
      shf = 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(909,bo(3),an(90),dh(3),n-4+shf,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(912,bo(3),an(2),dh(3),n-4+shf,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(912,bo(3),an(2),dh(3),n-5+shf,n-8+shf,n-1,0)
!
! deoxy-guanosine (DPG)
!
  else if ((resname .eq. 'DPG').OR.(resname.eq.'DIG')) then
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(874,bo(16),an(16),dh(16),c3i,c4i,c5i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(876,bo(17),an(11),dh(17),n-1,c3i,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(870,bo(15),an(17),dh(18),c4i,c3i,n-2,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    if (ua_model.eq.0) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(875,bo(2),an(3),an(3),n-3,c3i,n-2,-1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(875,bo(2),an(3),an(3),n-4,c3i,n-3,1)
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      call z_add(877,bo(2),an(3),an(3),n-4,n-5,n-3,1)
    else
      shf = 3
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nn1 = n
    call z_add(891,bo(63),an(20),an(18),n-5+shf,n-6+shf,n-4+shf,-1)
!    call z_add(891,bo(63),an(18),dh(15),n-5+shf,n-4+shf,c4i,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc1 = n
    nchi(i) = nchi(i) + 1
    chiline(1,i) = n
    call z_add(884,bo(39),an(98),chi(1,i),n-1,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    nc2 = n
    call z_add(889,bo(50),an(53),an(80),n-2,n-7+shf,n-1,1) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(885,bo(56),an(82),dh(1),n-2,n-3,n-1,0) 
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc2
    call z_add(888,bo(57),an(83),dh(1),n-2,n-4,n-3,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(883,bo(58),an(84),dh(3),n-4,n-5,n-3,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(880,bo(59),an(85),dh(3),n-1,n-5,n-6,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(886,bo(62),an(86),dh(3),n-4,n-6,n-7,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(878,bo(40),an(87),dh(1),n-2,n-3,n-7,0)
!   ring closure
    call z_add(-1,bo(1),an(1),dh(1),n-2,n-1,0,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(881,bo(60),an(66),dh(3),n-3,n-4,n-8,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(887,bo(61),an(81),dh(1),n-3,n-7,n-6,0)
    if (ua_model.lt.2) then
      shf = 0
      at(i)%nsc = at(i)%nsc + 1
      at(i)%sc(at(i)%nsc) = n
      atmres(n) = i
      eqatm(n) = nc1
      call z_add(890,bo(2),an(89),dh(3),n-9,n-7,n-11,0)
    else
      shf = 1
    end if
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    eqatm(n) = nc1
    call z_add(879,bo(3),an(90),dh(3),n-4+shf,n-6+shf,n-7+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(882,bo(3),an(2),dh(3),n-4+shf,n-7+shf,n-8+shf,0)
    at(i)%nsc = at(i)%nsc + 1
    at(i)%sc(at(i)%nsc) = n
    atmres(n) = i
    call z_add(882,bo(3),an(2),dh(3),n-5+shf,n-8+shf,n-1,0)
!
  else
!
    write(ilog,*) 'Fatal. Called nuc_sidechain(...) with unsupported&
 & residue type. This is either a bug or a feature which is not yet &
 &supported (offending residue: ',resname,').'
    call fexit()
!
  end if
!
end
!
