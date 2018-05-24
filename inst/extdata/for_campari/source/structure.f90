!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 3.0                                                           !
!                                                                          !
!    Copyright (C) 2017, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  !
!                        Davide Garolini, Jiri Vymetal                     !
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
! CONTRIBUTIONS: Rohit Pappu, Marco Bacci                                  !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ##########################################
! ##                                      ##
! ##    STRUCTURAL ANALYSIS AND FEATURE   ##
! ##         EXTRACTION ROUTINES          ##
! ##                                      ##
! ##########################################
!
!----------------------------------------------------------------------------
!
! set up the residue radii for topology-assisted cutoffs
!
subroutine setup_resrad()
!
  use iounit
  use sequen
  use aminos
!
  implicit none
!
  integer rs
  character(3) resname
!
! these values were recently revised to use a buffer of 0.5A for polymer residues
! small molecule values are only marginally buffered
!
  mrrd = 0.0
  do rs=1,nseq
    if (seqtyp(rs).le.0) cycle
    resname = amino(seqtyp(rs))
!   termini
    if (resname.EQ.'FOR') resrad(rs) = 1.7
    if (resname.EQ.'ACE') resrad(rs) = 2.9
    if (resname.EQ.'NH2') resrad(rs) = 1.5
    if (resname.EQ.'NME') resrad(rs) = 2.7
!   the odd ones
    if (resname.EQ.'GLY') resrad(rs) = 2.9
    if (resname.EQ.'PRO') resrad(rs) = 3.9
!   hydrophobics
    if (resname.EQ.'ALA') resrad(rs) = 2.9
    if (resname.EQ.'VAL') resrad(rs) = 4.0
    if (resname.EQ.'LEU') resrad(rs) = 5.2
    if (resname.EQ.'ILE') resrad(rs) = 5.2
    if (resname.EQ.'MET') resrad(rs) = 6.9
!   polars
    if (resname.EQ.'SER') resrad(rs) = 3.7
    if (resname.EQ.'THR') resrad(rs) = 4.0
    if (resname.EQ.'CYS') resrad(rs) = 4.3
    if (resname.EQ.'ASN') resrad(rs) = 5.1
    if (resname.EQ.'GLN') resrad(rs) = 6.4
    if (resname.EQ.'ASH') resrad(rs) = 4.8
    if (resname.EQ.'GLH') resrad(rs) = 6.2
    if (resname.EQ.'LYD') resrad(rs) = 7.6
    if (resname.EQ.'KAC') resrad(rs) = 10.1
!   chargeds
    if (resname.EQ.'ASP') resrad(rs) = 4.1
    if (resname.EQ.'GLU') resrad(rs) = 5.4
    if (resname.EQ.'LYS') resrad(rs) = 7.6
    if (resname.EQ.'KM1') resrad(rs) = 8.9
    if (resname.EQ.'KM2') resrad(rs) = 8.9
    if (resname.EQ.'KM3') resrad(rs) = 8.9
    if (resname.EQ.'ARG') resrad(rs) = 8.8
    if (resname.EQ.'TYO') resrad(rs) = 6.9
    if (resname.EQ.'CYX') resrad(rs) = 3.3
    if (resname.EQ.'SEP') resrad(rs) = 6.3
    if (resname.EQ.'TPO') resrad(rs) = 6.3
    if (resname.EQ.'PTR') resrad(rs) = 10.2
!   aromatics
    if (resname.EQ.'HIE') resrad(rs) = 6.2
    if (resname.EQ.'HID') resrad(rs) = 6.3
    if (resname.EQ.'HIP') resrad(rs) = 6.3
    if (resname.EQ.'PHE') resrad(rs) = 6.6
    if (resname.EQ.'TYR') resrad(rs) = 7.5
    if (resname.EQ.'TRP') resrad(rs) = 8.3
!   non-standards
    if (resname.EQ.'ABA') resrad(rs) = 4.0
    if (resname.EQ.'GAM') resrad(rs) = 4.0
    if (resname.EQ.'HYP') resrad(rs) = 5.4
    if (resname.EQ.'PCA') resrad(rs) = 5.0
    if (resname.EQ.'AIB') resrad(rs) = 2.9
    if (resname.EQ.'NVA') resrad(rs) = 5.2
    if (resname.EQ.'NLE') resrad(rs) = 6.5
    if (resname.EQ.'DAB') resrad(rs) = 5.2
    if (resname.EQ.'ORN') resrad(rs) = 6.5
!   nucleic acid residues (5' nucleotides require extra space -> parse_moltyp)
    if (resname.EQ.'RPU') resrad(rs) = 7.7
    if (resname.EQ.'DPU') resrad(rs) = 7.7
    if (resname.EQ.'RPC') resrad(rs) = 7.7
    if (resname.EQ.'DPC') resrad(rs) = 7.7
    if (resname.EQ.'RPT') resrad(rs) = 7.7
    if (resname.EQ.'DPT') resrad(rs) = 7.7
    if (resname.EQ.'RPA') resrad(rs) = 7.7
    if (resname.EQ.'DPA') resrad(rs) = 7.7
    if (resname.EQ.'RPG') resrad(rs) = 7.7
    if (resname.EQ.'DPG') resrad(rs) = 7.7
    if (resname.EQ.'R5P') resrad(rs) = 5.7
    if (resname.EQ.'D5P') resrad(rs) = 5.5
!   5'-nucleosides 
    if (resname.EQ.'RIU') resrad(rs) = 6.0
    if (resname.EQ.'DIU') resrad(rs) = 6.0
    if (resname.EQ.'RIC') resrad(rs) = 6.6
    if (resname.EQ.'DIC') resrad(rs) = 6.6
    if (resname.EQ.'RIT') resrad(rs) = 6.3
    if (resname.EQ.'DIT') resrad(rs) = 6.3
    if (resname.EQ.'RIA') resrad(rs) = 7.4
    if (resname.EQ.'DIA') resrad(rs) = 7.4
    if (resname.EQ.'RIG') resrad(rs) = 6.8
    if (resname.EQ.'DIG') resrad(rs) = 6.8
    if (resname.EQ.'RIB') resrad(rs) = 4.8
    if (resname.EQ.'DIB') resrad(rs) = 4.7
!   small molecules (only marginally buffered)
    if (resname.EQ.'URE') resrad(rs) = 2.1
    if (resname.EQ.'GDN') resrad(rs) = 2.1
    if (resname.EQ.'FOA') resrad(rs) = 2.1
    if (resname.EQ.'CH4') resrad(rs) = 1.1
    if (resname.EQ.'NH4') resrad(rs) = 1.1
    if (resname.EQ.'SPC') resrad(rs) = 1.1
    if (resname.EQ.'T3P') resrad(rs) = 1.1
    if (resname.EQ.'T4P') resrad(rs) = 1.1
    if (resname.EQ.'T4E') resrad(rs) = 1.1
    if (resname.EQ.'T5P') resrad(rs) = 1.1
    if (resname.EQ.'NMF') resrad(rs) = 2.3
    if (resname.EQ.'AC-') resrad(rs) = 2.3 
    if (resname.EQ.'NMA') resrad(rs) = 3.4
    if (resname.EQ.'ACA') resrad(rs) = 2.2
    if (resname.EQ.'PPA') resrad(rs) = 3.5
    if (resname.EQ.'NA+') resrad(rs) = 0.0
    if (resname.EQ.'CL-') resrad(rs) = 0.0
    if (resname.EQ.'K+ ') resrad(rs) = 0.0
    if (resname.EQ.'BR-') resrad(rs) = 0.0
    if (resname.EQ.'CS+') resrad(rs) = 0.0
    if (resname.EQ.'I- ') resrad(rs) = 0.0
    if (resname.EQ.'MOH') resrad(rs) = 2.0
    if (resname.EQ.'DMA') resrad(rs) = 3.5
    if (resname.EQ.'PCR') resrad(rs) = 4.7
    if (resname.EQ.'CYT') resrad(rs) = 3.8
    if (resname.EQ.'URA') resrad(rs) = 3.8
    if (resname.EQ.'THY') resrad(rs) = 4.1
    if (resname.EQ.'ADE') resrad(rs) = 3.8
    if (resname.EQ.'GUA') resrad(rs) = 4.3
    if (resname.EQ.'PRP') resrad(rs) = 2.2
    if (resname.EQ.'NBU') resrad(rs) = 3.5
    if (resname.EQ.'IBU') resrad(rs) = 2.2
    if (resname.EQ.'TOL') resrad(rs) = 4.0
    if (resname.EQ.'EOH') resrad(rs) = 3.3
    if (resname.EQ.'MSH') resrad(rs) = 2.4
    if (resname.EQ.'EMT') resrad(rs) = 3.8
    if (resname.EQ.'IMD') resrad(rs) = 3.3
    if (resname.EQ.'IME') resrad(rs) = 3.3
    if (resname.EQ.'MIN') resrad(rs) = 4.0
    if (resname.EQ.'LCP') resrad(rs) = 1.5
    if (resname.EQ.'NO3') resrad(rs) = 1.4
    if (resname.EQ.'1MN') resrad(rs) = 2.2
    if (resname.EQ.'2MN') resrad(rs) = 2.2
    if (resname.EQ.'BEN') resrad(rs) = 4.0
    if (resname.EQ.'NAP') resrad(rs) = 4.0
    if (resname.EQ.'PUR') resrad(rs) = 3.5
    if (resname.EQ.'O2 ') resrad(rs) = 1.25
!
    if (resrad(rs).gt.mrrd) mrrd = resrad(rs)
  end do
!
! for the atom numbers, a maximum titration state for normal titratable groups can be acommodated
  do rs=1,nseq
    if (seqtyp(rs).le.0) cycle
    resname = amino(seqtyp(rs))
!   peptide termini
    if (resname.EQ.'FOR') natres(rs) = 3
    if (resname.EQ.'ACE') natres(rs) = 6
    if (resname.EQ.'NH2') natres(rs) = 3
    if (resname.EQ.'NME') natres(rs) = 6
!   the odd ones
    if (resname.EQ.'GLY') natres(rs) = 10 + 1
    if (resname.EQ.'PRO') natres(rs) = 9 + 9
!   hydrophobics
    if (resname.EQ.'ALA') natres(rs) = 10 + 4
    if (resname.EQ.'VAL') natres(rs) = 10 + 10
    if (resname.EQ.'LEU') natres(rs) = 10 + 13
    if (resname.EQ.'ILE') natres(rs) = 10 + 13
    if (resname.EQ.'MET') natres(rs) = 10 + 11
!   polars
    if (resname.EQ.'SER') natres(rs) = 10 + 5
    if (resname.EQ.'THR') natres(rs) = 10 + 8
    if (resname.EQ.'CYS') natres(rs) = 10 + 5
    if (resname.EQ.'ASN') natres(rs) = 10 + 8
    if (resname.EQ.'GLN') natres(rs) = 10 + 11
    if (resname.EQ.'ASH') natres(rs) = 10 + 7
    if (resname.EQ.'GLH') natres(rs) = 10 + 10
    if (resname.EQ.'LYD') natres(rs) = 10 + 16
    if (resname.eq.'KAC') natres(rs) = 10 + 20
!   chargeds
    if (resname.EQ.'ASP') natres(rs) = 10 + 7
    if (resname.EQ.'GLU') natres(rs) = 10 + 10
    if (resname.EQ.'LYS') natres(rs) = 10 + 16
    if (resname.EQ.'KM1') natres(rs) = 10 + 19
    if (resname.EQ.'KM2') natres(rs) = 10 + 22
    if (resname.EQ.'KM3') natres(rs) = 10 + 25
    if (resname.EQ.'ARG') natres(rs) = 10 + 18
    if (resname.EQ.'TYO') natres(rs) = 10 + 15
    if (resname.EQ.'CYX') natres(rs) = 10 + 5
    if (resname.EQ.'SEP') natres(rs) = 10 + 10
    if (resname.EQ.'TPO') natres(rs) = 10 + 13
    if (resname.EQ.'PTR') natres(rs) = 10 + 20
!   aromatics
    if (resname.EQ.'HIE') natres(rs) = 10 + 11
    if (resname.EQ.'HID') natres(rs) = 10 + 11
    if (resname.EQ.'HIP') natres(rs) = 10 + 12
    if (resname.EQ.'PHE') natres(rs) = 10 + 14
    if (resname.EQ.'TYR') natres(rs) = 10 + 15
    if (resname.EQ.'TRP') natres(rs) = 10 + 18
!   non-standards
    if (resname.EQ.'ABA') natres(rs) = 10 + 7
    if (resname.EQ.'GAM') natres(rs) = 10 + 7
    if (resname.EQ.'HYP') natres(rs) = 9 + 10
    if (resname.EQ.'PCA') natres(rs) = 9 + 8
    if (resname.EQ.'AIB') natres(rs) = 9 + 8
    if (resname.EQ.'NVA') natres(rs) = 10 + 10
    if (resname.EQ.'NLE') natres(rs) = 10 + 13
    if (resname.EQ.'DAB') natres(rs) = 10 + 10
    if (resname.EQ.'ORN') natres(rs) = 10 + 13
!   nucleic acid residues
    if (resname.EQ.'RPU') natres(rs) = 13 + 20
    if (resname.EQ.'DPU') natres(rs) = 12 + 20
    if (resname.EQ.'RPC') natres(rs) = 13 + 21
    if (resname.EQ.'DPC') natres(rs) = 12 + 21
    if (resname.EQ.'RPT') natres(rs) = 13 + 23
    if (resname.EQ.'DPT') natres(rs) = 12 + 23
    if (resname.EQ.'RPA') natres(rs) = 13 + 23
    if (resname.EQ.'DPA') natres(rs) = 12 + 23
    if (resname.EQ.'RPG') natres(rs) = 13 + 24
    if (resname.EQ.'DPG') natres(rs) = 12 + 24
    if (resname.EQ.'R5P') natres(rs) = 13 + 12
    if (resname.EQ.'D5P') natres(rs) = 12 + 12
!   5-'nucleosides
    if (resname.EQ.'RIU') natres(rs) = 9 + 20
    if (resname.EQ.'DIU') natres(rs) = 8 + 20
    if (resname.EQ.'RIC') natres(rs) = 9 + 21
    if (resname.EQ.'DIC') natres(rs) = 8 + 21
    if (resname.EQ.'RIT') natres(rs) = 9 + 23
    if (resname.EQ.'DIT') natres(rs) = 8 + 23
    if (resname.EQ.'RIA') natres(rs) = 9 + 23
    if (resname.EQ.'DIA') natres(rs) = 8 + 23
    if (resname.EQ.'RIG') natres(rs) = 9 + 24
    if (resname.EQ.'DIG') natres(rs) = 8 + 24
    if (resname.EQ.'RIB') natres(rs) = 9 + 12
    if (resname.EQ.'DIB') natres(rs) = 8 + 12
!   small molecules
    if (resname.EQ.'URE') natres(rs) = 8
    if (resname.EQ.'GDN') natres(rs) = 10
    if (resname.EQ.'FOA') natres(rs) = 6
    if (resname.EQ.'CH4') natres(rs) = 5
    if (resname.EQ.'NH4') natres(rs) = 5
    if (resname.EQ.'SPC') natres(rs) = 3
    if (resname.EQ.'T3P') natres(rs) = 3
    if (resname.EQ.'T4P') natres(rs) = 4
    if (resname.EQ.'T4E') natres(rs) = 4
    if (resname.EQ.'T5P') natres(rs) = 5
    if (resname.EQ.'NMF') natres(rs) = 9
    if (resname.EQ.'AC-') natres(rs) = 8 
    if (resname.EQ.'NMA') natres(rs) = 12
    if (resname.EQ.'ACA') natres(rs) = 9
    if (resname.EQ.'PPA') natres(rs) = 12
    if (resname.EQ.'NA+') natres(rs) = 1
    if (resname.EQ.'CL-') natres(rs) = 1
    if (resname.EQ.'K+ ') natres(rs) = 1
    if (resname.EQ.'BR-') natres(rs) = 1
    if (resname.EQ.'CS+') natres(rs) = 1
    if (resname.EQ.'I- ') natres(rs) = 1
    if (resname.EQ.'MOH') natres(rs) = 6
    if (resname.EQ.'DMA') natres(rs) = 15
    if (resname.EQ.'PCR') natres(rs) = 16
    if (resname.EQ.'CYT') natres(rs) = 14
    if (resname.EQ.'URA') natres(rs) = 13
    if (resname.EQ.'THY') natres(rs) = 16
    if (resname.EQ.'ADE') natres(rs) = 16
    if (resname.EQ.'GUA') natres(rs) = 17
    if (resname.EQ.'PRP') natres(rs) = 11
    if (resname.EQ.'NBU') natres(rs) = 14
    if (resname.EQ.'IBU') natres(rs) = 14
    if (resname.EQ.'TOL') natres(rs) = 15
    if (resname.EQ.'EOH') natres(rs) = 9
    if (resname.EQ.'MSH') natres(rs) = 6
    if (resname.EQ.'EMT') natres(rs) = 12
    if (resname.EQ.'IME') natres(rs) = 13
    if (resname.EQ.'IMD') natres(rs) = 13
    if (resname.EQ.'MIN') natres(rs) = 21
    if (resname.EQ.'LCP') natres(rs) = 5
    if (resname.EQ.'NO3') natres(rs) = 4
    if (resname.EQ.'1MN') natres(rs) = 8
    if (resname.EQ.'2MN') natres(rs) = 11
    if (resname.EQ.'BEN') natres(rs) = 12
    if (resname.EQ.'NAP') natres(rs) = 18
    if (resname.EQ.'PUR') natres(rs) = 13
    if (resname.EQ.'O2 ') natres(rs) = 2
!
    if (natres(rs).le.0) then
      write(ilog,*) 'Fatal. Number of atoms in residue is zero for r&
 &esidue ',rs,' (',resname,'). This is an omission bug.'
      call fexit()
    end if
!
  end do
!
end

!
!-----------------------------------------------------------------------
!
! a few utility functions used by the amide PC calculator
! note that these functions crucially rely on the arrays in polypep.i
! as well as the order of atoms in the sidechains as in sidechain.f
!
subroutine get_COHN_bbbb(rs1,rs2,v1,v2,v3,v4)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
!
  implicit none
!
  RTYPE v1(3),v2(3),v3(3),v4(3)
  integer rs1,rs2
!
  if (((seqflag(rs1).eq.2).OR.(seqflag(rs1).eq.5).OR.(seqflag(rs1).eq.8).OR.(seqflag(rs1).eq.10)).AND.&
 &    (oi(rs1).gt.0).AND.(ci(rs1).gt.0)) then
    v1(1) = x(ci(rs1))
    v1(2) = y(ci(rs1))
    v1(3) = z(ci(rs1))
    v2(1) = x(oi(rs1))
    v2(2) = y(oi(rs1))
    v2(3) = z(oi(rs1))
  else
    write(ilog,*) 'Fatal. Called get_COHN_bbbb(...) with improper&
 & residue (Nr. ',rs1,': ',amino(seqtyp(rs1)),').'
    call fexit()
  end if
  if (((seqflag(rs2).eq.2).OR.(seqflag(rs2).eq.8).OR.(seqflag(rs2).eq.11)).AND.&
 &    (ni(rs2).gt.0).AND.(hni(rs2).gt.0)) then
    v3(1) = x(hni(rs2))
    v3(2) = y(hni(rs2))
    v3(3) = z(hni(rs2))
    v4(1) = x(ni(rs2))
    v4(2) = y(ni(rs2))
    v4(3) = z(ni(rs2))
  else
    write(ilog,*) 'Fatal. Called get_COHN_bbbb(...) with improper&
 & residue (Nr. ',rs2,': ',amino(seqtyp(rs2)),').'
    call fexit()
  end if
!
end
!
!----------------------------------
!
subroutine get_COHN_scbb(rs1,rs2,v1,v2,v3,v4)
!
  use iounit
  use sequen
  use polypep
  use atoms
  use aminos
  use system
!
  implicit none
!
  RTYPE v1(3),v2(3),v3(3),v4(3)
  integer rs1,rs2,shf
  character(3) resname
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  resname = amino(seqtyp(rs1))
  if (resname.eq.'GLN') then
    v1(1) = x(at(rs1)%sc(4-shf))
    v1(2) = y(at(rs1)%sc(4-shf))
    v1(3) = z(at(rs1)%sc(4-shf))
    v2(1) = x(at(rs1)%sc(5-shf))
    v2(2) = y(at(rs1)%sc(5-shf))
    v2(3) = z(at(rs1)%sc(5-shf))
  else if (resname.eq.'ASN') then
    v1(1) = x(at(rs1)%sc(3-shf))
    v1(2) = y(at(rs1)%sc(3-shf))
    v1(3) = z(at(rs1)%sc(3-shf))
    v2(1) = x(at(rs1)%sc(4-shf))
    v2(2) = y(at(rs1)%sc(4-shf))
    v2(3) = z(at(rs1)%sc(4-shf))
  else
    write(ilog,*) 'Fatal. Called get_COHN_scbb(...) with improper re&
 &sidue (Nr. ',rs1,': ',resname,').'
    call fexit()
  end if
  if (((seqflag(rs2).eq.2).OR.(seqflag(rs2).eq.8).OR.(seqflag(rs2).eq.11)).AND.&
 &    (ni(rs2).gt.0).AND.(hni(rs2).gt.0)) then
    v3(1) = x(hni(rs2))
    v3(2) = y(hni(rs2))
    v3(3) = z(hni(rs2))
    v4(1) = x(ni(rs2))
    v4(2) = y(ni(rs2))
    v4(3) = z(ni(rs2))
  else
    write(ilog,*) 'Fatal. Called get_COHN_scbb(...) with improper&
 & residue (Nr. ',rs2,': ',amino(seqtyp(rs2)),').'
    call fexit()
  end if
!
end
!
!----------------------------------
!
subroutine get_COHHN_bbsc(rs1,rs2,v1,v2,v3,v4,v5)
!
  use iounit
  use sequen
  use polypep
  use atoms
  use aminos
  use system
!
  implicit none
!
  RTYPE v1(3),v2(3),v3(3),v4(3),v5(3)
  integer rs1,rs2,shf,shf3,shf5
!
  shf = 0
  shf3 = 0
  shf5 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf3 = 3
    shf5 = 5
  end if
!
  if (((seqflag(rs1).eq.2).OR.(seqflag(rs1).eq.5).OR.(seqflag(rs1).eq.8).OR.(seqflag(rs1).eq.10)).AND.&
 &    (oi(rs1).gt.0).AND.(ci(rs1).gt.0)) then
    v1(1) = x(ci(rs1))
    v1(2) = y(ci(rs1))
    v1(3) = z(ci(rs1))
    v2(1) = x(oi(rs1))
    v2(2) = y(oi(rs1))
    v2(3) = z(oi(rs1))
  else
    write(ilog,*) 'Fatal. Called get_COHHN_bbsc(...) with improper r&
 &esidue (Nr. ',rs1,': ',amino(seqtyp(rs1)),')'
    call fexit()
  end if
  if (seqtyp(rs2).eq.19) then
!    write(ilog,*) 'Got ',attyp(at(rs2)%sc(11)),' as GLN-H.'
    v3(1) = x(at(rs2)%sc(11-shf5))
    v3(2) = y(at(rs2)%sc(11-shf5))
    v3(3) = z(at(rs2)%sc(11-shf5))
!    write(ilog,*) 'Got ',attyp(at(rs2)%sc(12)),' as GLN-H.'
    v4(1) = x(at(rs2)%sc(12-shf5))
    v4(2) = y(at(rs2)%sc(12-shf5))
    v4(3) = z(at(rs2)%sc(12-shf5))
!    write(ilog,*) 'Got ',attyp(at(rs2)%sc(6)),' as GLN-N.'
    v5(1) = x(at(rs2)%sc(6-shf))
    v5(2) = y(at(rs2)%sc(6-shf))
    v5(3) = z(at(rs2)%sc(6-shf))
  else if (seqtyp(rs2).eq.17) then
!    write(ilog,*) 'Got ',attyp(at(rs2)%sc(8)),' as ASN-H.'
    v3(1) = x(at(rs2)%sc(8-shf3))
    v3(2) = y(at(rs2)%sc(8-shf3))
    v3(3) = z(at(rs2)%sc(8-shf3))
!    write(ilog,*) 'Got ',attyp(at(rs2)%sc(9)),' as ASN-H.'
    v4(1) = x(at(rs2)%sc(9-shf3))
    v4(2) = y(at(rs2)%sc(9-shf3))
    v4(3) = z(at(rs2)%sc(9-shf3))
!    write(ilog,*) 'Got ',attyp(at(rs2)%sc(5)),' as ASN-N.'
    v5(1) = x(at(rs2)%sc(5-shf))
    v5(2) = y(at(rs2)%sc(5-shf))
    v5(3) = z(at(rs2)%sc(5-shf))
  else
    write(ilog,*) 'Fatal. Called get_OHHNC_bbsc(...) with improper r&
 &esidue (Nr. ',rs2,': ',amino(seqtyp(rs2)),')'
    call fexit()
  end if
!
end
!
!----------------------------------
!    
subroutine get_COHHN_scsc(rs1,rs2,v1,v2,v3,v4,v5)
!
  use iounit
  use polypep
  use sequen
  use atoms
  use aminos
  use system
!
  implicit none
!
  RTYPE v1(3),v2(3),v3(3),v4(3),v5(3)
  integer rs1,rs2,shf,shf3,shf5
!
  shf = 0
  shf3 = 0
  shf5 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf3 = 3
    shf5 = 5
  end if
!
  if (seqtyp(rs1).eq.19) then
    v1(1) = x(at(rs1)%sc(4-shf))
    v1(2) = y(at(rs1)%sc(4-shf))
    v1(3) = z(at(rs1)%sc(4-shf))
    v2(1) = x(at(rs1)%sc(5-shf))
    v2(2) = y(at(rs1)%sc(5-shf))
    v2(3) = z(at(rs1)%sc(5-shf))
  else if (seqtyp(rs1).eq.17) then
    v1(1) = x(at(rs1)%sc(3-shf))
    v1(2) = y(at(rs1)%sc(3-shf))
    v1(3) = z(at(rs1)%sc(3-shf))
    v2(1) = x(at(rs1)%sc(4-shf))
    v2(2) = y(at(rs1)%sc(4-shf))
    v2(3) = z(at(rs1)%sc(4-shf))
  else
    write(ilog,*) 'Fatal. Called get_COHHN_scsc(...) with improper r&
 &esidue (Nr. ',rs1,': ',amino(seqtyp(rs1)),').'
    call fexit()
  end if
  if (seqtyp(rs2).eq.19) then
    v3(1) = x(at(rs2)%sc(11-shf5))
    v3(2) = y(at(rs2)%sc(11-shf5))
    v3(3) = z(at(rs2)%sc(11-shf5))
    v4(1) = x(at(rs2)%sc(12-shf5))
    v4(2) = y(at(rs2)%sc(12-shf5))
    v4(3) = z(at(rs2)%sc(12-shf5))
    v5(1) = x(at(rs2)%sc(6-shf))
    v5(2) = y(at(rs2)%sc(6-shf))
    v5(3) = z(at(rs2)%sc(6-shf))
  else if (seqtyp(rs2).eq.17) then
    v3(1) = x(at(rs2)%sc(8-shf3))
    v3(2) = y(at(rs2)%sc(8-shf3))
    v3(3) = z(at(rs2)%sc(8-shf3))
    v4(1) = x(at(rs2)%sc(9-shf3))
    v4(2) = y(at(rs2)%sc(9-shf3))
    v4(3) = z(at(rs2)%sc(9-shf3))
    v5(1) = x(at(rs2)%sc(5-shf))
    v5(2) = y(at(rs2)%sc(5-shf))
    v5(3) = z(at(rs2)%sc(5-shf))
  else
    write(ilog,*) 'Fatal. Called get_COHHN_scsc(...) with improper r&
 &esidue (Nr. ',rs2,': ',amino(seqtyp(rs2)),').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------
!
subroutine get_CONHHCONHH_scsc(rs1,rs2,v1,v2,v3,v4,v5,v6,v7,v8,v9,&
 &v10)
!
  use iounit
  use polypep
  use sequen
  use atoms
  use aminos
  use system
!
  implicit none
!
  RTYPE v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3),v9(3)
  RTYPE v10(3)
  integer rs1,rs2,shf,shf3,shf5
!
  shf = 0
  shf3 = 0
  shf5 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf3 = 3
    shf5 = 5
  end if
!
  if (seqtyp(rs1).eq.19) then
    v1(1) = x(at(rs1)%sc(4-shf))
    v1(2) = y(at(rs1)%sc(4-shf))
    v1(3) = z(at(rs1)%sc(4-shf))
    v2(1) = x(at(rs1)%sc(5-shf))
    v2(2) = y(at(rs1)%sc(5-shf))
    v2(3) = z(at(rs1)%sc(5-shf))
    v3(1) = x(at(rs1)%sc(6-shf))
    v3(2) = y(at(rs1)%sc(6-shf))
    v3(3) = z(at(rs1)%sc(6-shf))
    v4(1) = x(at(rs1)%sc(11-shf5))
    v4(2) = y(at(rs1)%sc(11-shf5))
    v4(3) = z(at(rs1)%sc(11-shf5))
    v5(1) = x(at(rs1)%sc(12-shf5))
    v5(2) = y(at(rs1)%sc(12-shf5))
    v5(3) = z(at(rs1)%sc(12-shf5))
  else if (seqtyp(rs1).eq.17) then
    v1(1) = x(at(rs1)%sc(3-shf))
    v1(2) = y(at(rs1)%sc(3-shf))
    v1(3) = z(at(rs1)%sc(3-shf))
    v2(1) = x(at(rs1)%sc(4-shf))
    v2(2) = y(at(rs1)%sc(4-shf))
    v2(3) = z(at(rs1)%sc(4-shf))
    v3(1) = x(at(rs1)%sc(5-shf))
    v3(2) = y(at(rs1)%sc(5-shf))
    v3(3) = z(at(rs1)%sc(5-shf))
    v4(1) = x(at(rs1)%sc(8-shf3))
    v4(2) = y(at(rs1)%sc(8-shf3))
    v4(3) = z(at(rs1)%sc(8-shf3))
    v5(1) = x(at(rs1)%sc(9-shf3))
    v5(2) = y(at(rs1)%sc(9-shf3))
    v5(3) = z(at(rs1)%sc(9-shf3))
  else
    write(ilog,*) 'Fatal. Called get_CONHHCONHH_scsc(...) with impro&
 &per residue (Nr. ',rs1,': ',amino(seqtyp(rs1)),').'
    call fexit()
  end if
  if (seqtyp(rs2).eq.19) then
    v6(1) = x(at(rs2)%sc(4-shf))
    v6(2) = y(at(rs2)%sc(4-shf))
    v6(3) = z(at(rs2)%sc(4-shf))
    v7(1) = x(at(rs2)%sc(5-shf))
    v7(2) = y(at(rs2)%sc(5-shf))
    v7(3) = z(at(rs2)%sc(5-shf))
    v8(1) = x(at(rs2)%sc(6-shf))
    v8(2) = y(at(rs2)%sc(6-shf))
    v8(3) = z(at(rs2)%sc(6-shf))
    v9(1) = x(at(rs2)%sc(11-shf5))
    v9(2) = y(at(rs2)%sc(11-shf5))
    v9(3) = z(at(rs2)%sc(11-shf5))
    v10(1) = x(at(rs2)%sc(12-shf5))
    v10(2) = y(at(rs2)%sc(12-shf5))
    v10(3) = z(at(rs2)%sc(12-shf5))
  else if (seqtyp(rs2).eq.17) then
    v6(1) = x(at(rs2)%sc(3-shf))
    v6(2) = y(at(rs2)%sc(3-shf))
    v6(3) = z(at(rs2)%sc(3-shf))
    v7(1) = x(at(rs2)%sc(4-shf))
    v7(2) = y(at(rs2)%sc(4-shf))
    v7(3) = z(at(rs2)%sc(4-shf))
    v8(1) = x(at(rs2)%sc(5-shf))
    v8(2) = y(at(rs2)%sc(5-shf))
    v8(3) = z(at(rs2)%sc(5-shf))
    v9(1) = x(at(rs2)%sc(8-shf3))
    v9(2) = y(at(rs2)%sc(8-shf3))
    v9(3) = z(at(rs2)%sc(8-shf3))
    v10(1) = x(at(rs2)%sc(9-shf3))
    v10(2) = y(at(rs2)%sc(9-shf3))
    v10(3) = z(at(rs2)%sc(9-shf3))
  else
    write(ilog,*) 'Fatal. Called get_CONHHCONHH_scsc(...) with impro&
 &per residue (Nr. ',rs2,': ',amino(seqtyp(rs2)),').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------
!
subroutine get_CONHH_sc(rs1,v1,v2,v3,v4,v5)
!
  use iounit
  use polypep
  use sequen
  use atoms
  use aminos
  use system
!
  implicit none
!
  RTYPE v1(3),v2(3),v3(3),v4(3),v5(3)
  integer rs1,shf,shf3,shf5
!
  shf = 0
  shf3 = 0
  shf5 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf3 = 3
    shf5 = 5
  end if
!
  if (seqtyp(rs1).eq.19) then
    v1(1) = x(at(rs1)%sc(4-shf))
    v1(2) = y(at(rs1)%sc(4-shf))
    v1(3) = z(at(rs1)%sc(4-shf))
    v2(1) = x(at(rs1)%sc(5-shf))
    v2(2) = y(at(rs1)%sc(5-shf))
    v2(3) = z(at(rs1)%sc(5-shf))
    v3(1) = x(at(rs1)%sc(6-shf))
    v3(2) = y(at(rs1)%sc(6-shf))
    v3(3) = z(at(rs1)%sc(6-shf))
    v4(1) = x(at(rs1)%sc(11-shf5))
    v4(2) = y(at(rs1)%sc(11-shf5))
    v4(3) = z(at(rs1)%sc(11-shf5))
    v5(1) = x(at(rs1)%sc(12-shf5))
    v5(2) = y(at(rs1)%sc(12-shf5))
    v5(3) = z(at(rs1)%sc(12-shf5))
  else if (seqtyp(rs1).eq.17) then
    v1(1) = x(at(rs1)%sc(3-shf))
    v1(2) = y(at(rs1)%sc(3-shf))
    v1(3) = z(at(rs1)%sc(3-shf))
    v2(1) = x(at(rs1)%sc(4-shf))
    v2(2) = y(at(rs1)%sc(4-shf))
    v2(3) = z(at(rs1)%sc(4-shf))
    v3(1) = x(at(rs1)%sc(5-shf))
    v3(2) = y(at(rs1)%sc(5-shf))
    v3(3) = z(at(rs1)%sc(5-shf))
    v4(1) = x(at(rs1)%sc(8-shf3))
    v4(2) = y(at(rs1)%sc(8-shf3))
    v4(3) = z(at(rs1)%sc(8-shf3))
    v5(1) = x(at(rs1)%sc(9-shf3))
    v5(2) = y(at(rs1)%sc(9-shf3))
    v5(3) = z(at(rs1)%sc(9-shf3))
  else
    write(ilog,*) 'Fatal. Called get_CONHH_sc(...) with improper res&
 &idue (Nr. ',rs1,': ',amino(seqtyp(rs1)),').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------
!
!  ROUTINES NEEDED FOR ANALYSIS OF RESIDUE CONTACTS AND MOLECULAR CLUSTERS
!
!----------------------------------------------------------------------------
!
! this routine takes care of averaging contact information
! instantaneous values for total numbers of contacts are available
! through c1,c2
! Commented HACK: to give info on the closest BOUND_PARTNER
!
subroutine rescontacts()
!
  use iounit
  use atoms
  use polypep
  use sequen
  use contacts
  use molecule
  use cutoffs
  use mcsums
  use system
  use grandensembles
  use pdb
! use paircorr
!
  implicit none
!
  integer rs,rs2,c1,c2,i,j,k,ncoc,imol,jmol,ii,kk,cnt1,cnt2,cnt3,cnt4
  integer curclu,nhits,cobin!,bpcounter
  RTYPE mindis2,comdis2,com1(3),com2(3),cutmax,cdis2,dis2,coo,scal
  integer, ALLOCATABLE:: inclu(:),cluwt(:,:)
  logical doclu,inrange,ismember!,bpfound,doinst,finished
!
!  integer nbpres                  !number of residues to which binding partnes will be found. Should be read in input file.
!  integer, ALLOCATABLE:: bpres(:) !residue id of the res.s of which we want the binding partners. Should be read from input file.
!  RTYPE, ALLOCATABLE:: bpdata(:,:) 
!!
!  doinst = .false.       !for now everything is controlled as GENERAL_DIST.dat in PCCALC. I borrow from there.
!  finished = .false.
!  if (inst_gpc.gt.0) then
!    if (mod(n_pccalls,inst_gpc).eq.0) doinst = .true.
!  end if
!  if (doinst.EQV..true.) then
!    bpfound = .false.          !important
!    bpcounter = 1
!    nbpres = 42
!    allocate(bpres(nbpres))    !This should contain the list of residues we want to find the binding parter of
!    allocate(bpdata(nbpres,4)) !1: binding partner according to com-com distances. 2: binding partner acc. to min atom-atom dis
!                               !3: com-com distance to binding partner. 4: minimum atom-atom distance to binding partner
!  end if
!
  ncoc = nstep/contactcalc
  if (mod(ncoc,clucalc).eq.0) then
    doclu = .true.
  else
    doclu = .false.
  end if
!
  if (doclu.EQV..true.) then
    allocate(inclu(nsolutes))
    allocate(cluwt(nsolutes,2))
    cluwt(:,:) = 0  ! # of links and size per cluster
    inclu(:) = 0
  end if
!
  if (use_frame_weights.EQV..true.) then
    scal = 1.0*framewts(curframe)
  else
    scal = 1.0
  end if
!
! some initialization
  c1 = 0
  c2 = 0
  if (contact_cuts(1).gt.contact_cuts(2)) then
    cutmax = sqrt(contact_cuts(1))
  else
    cutmax = sqrt(contact_cuts(2))
  end if
!
!write(ilog,*) '----------------- initialization --------------------'
!  if (doinst.EQV..true.) then
!    do i=1,nbpres
!      bpres(i) = i             !this should be read from an input file. IT HAS TO BE SORTED!
!      bpdata(i,1) = 0.
!      bpdata(i,2) = 0.
!      bpdata(i,3) = cutmax      !bad coding. should be initialize elsewhere.
!      bpdata(i,4) = cutmax
!      !write(ilog,*) 'i,bpres(i),bpdata(1:4)',i,bpres(i),bpdata(i,:)
!    end do
!  end if
!  !write(ilog,*) '----------------------------------------------------'
!
  cnt1 = 0
  do i=1,nsolutes
    cnt3 = cnt1
    imol = solutes(i)
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) then
          cnt1 = cnt3 + rsmol(imol,2) - rsmol(imol,1) + 1
          cycle
        end if
      end if
    end if
    ncontactavg(i,i) = ncontactavg(i,i) + 1
    if (doclu.EQV..true.) then
      if (inclu(i).le.0) then
        inclu(i) = i
        cluwt(inclu(i),2) = 1
      end if
    end if
    do rs=rsmol(imol,1),rsmol(imol,2)
      cnt1 = cnt1 + 1
      ii = refat(rs)
      call rescom(rs,com1)
      do rs2=rs+contact_off,rsmol(imol,2)
        if (rs.eq.rs2) then
          contact_maps(cnt1,cnt1) = contact_maps(cnt1,cnt1) + scal
          c1 = c1 + 1
          c2 = c2 + 1
          cycle
        end if
        cdis2 = (cutmax+resrad(rs)+resrad(rs2))**2
        kk = refat(rs2)
        call dis_bound(ii,kk,dis2)
        if (dis2.lt.cdis2) then
!         now get the correct com-distance
          call rescom(rs2,com2)
          call dis_bound2(com1,com2,comdis2)
          call respmindis(rs,rs2,mindis2)
          if (mindis2.lt.contact_cuts(1)) then
            contact_maps(cnt1+rs2-rs,cnt1) = contact_maps(cnt1+rs2-rs,cnt1) + scal
            c1 = c1 + 1
          end if
          if (comdis2.lt.contact_cuts(2)) then
            contact_maps(cnt1,cnt1+rs2-rs) = contact_maps(cnt1,cnt1+rs2-rs) + scal
            c2 = c2 + 1
          end if
        end if
      end do
    end do
    cnt4 = cnt3 + rsmol(imol,2) - rsmol(imol,1) + 1
    do j=i+1,nsolutes
      cnt1 = cnt3
      inrange = .false.
      jmol = solutes(j)
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(jmol)).EQV..true.) then
          if (ismember(ispresent,jmol).EQV..false.) then
            cnt4 = cnt4 + rsmol(jmol,2) - rsmol(jmol,1) + 1
            cnt1 = cnt1 + rsmol(imol,2) - rsmol(imol,1) + 1
            cycle
          end if
        end if
      end if
      ncontactavg(i,j) = ncontactavg(i,j) + 1
      do rs=rsmol(imol,1),rsmol(imol,2)
!        if (doinst.EQV..true.) then
!          if (bpfound.EQV..true.) then
!            if (bpcounter.lt.nbpres) then
!              bpcounter = bpcounter + 1
!            else
!              finished = .true.
!            end if
!          end if
!          bpfound = .false.
!          if ((bpfound.EQV..true.).AND.(finished.EQV..true.)) then
!            write(ilog,*) 'Fatal. Found new bp reference residue once the search is complete. This is a bug.'
!            call fexit()
!          end if
!          if (rs.eq.bpres(bpcounter)) then
!            bpfound = .true.
!            !write(ilog,*) 'ref. residue found, bpcounter',bpcounter,bpres(bpcounter),rs
!          end if
!        end if
        cnt1 = cnt1 + 1
        cnt2 = cnt4
        ii = refat(rs)
        call rescom(rs,com1)
        do rs2=rsmol(jmol,1),rsmol(jmol,2)
          cnt2 = cnt2 + 1
          cdis2 = (cutmax+resrad(rs)+resrad(rs2))**2
          kk = refat(rs2)
          call dis_bound(ii,kk,dis2)
          if (dis2.lt.cdis2) then
            call rescom(rs2,com2)
            call dis_bound2(com1,com2,comdis2)
            call respmindis(rs,rs2,mindis2)
!            if (doinst.EQV..true.) then
!              if (bpfound.EQV..true.) then
!                if (bpdata(bpcounter,3).gt.dsqrt(comdis2)) then
!                  !write(ilog,*) 'Found closer partner com-com:',rs2,'from',bpdata(bpcounter,3),'to',dsqrt(comdis2),''
!                  bpdata(bpcounter,1) = 1.0*rs2
!                  bpdata(bpcounter,3) = dsqrt(comdis2)
!                end if
!                if (bpdata(bpcounter,4).gt.dsqrt(mindis2)) then
!                  !write(ilog,*) 'Found closer partner atom-atom:',rs2,'from',bpdata(bpcounter,4),'to',dsqrt(mindis2),''
!                  bpdata(bpcounter,2) = 1.0*rs2
!                  bpdata(bpcounter,4) = dsqrt(mindis2)
!                end if
!              end if
!            end if
            if (mindis2.lt.contact_cuts(1)) then
              contact_maps(cnt2,cnt1) = contact_maps(cnt2,cnt1) + scal
              c1 = c1 + 1
              inrange = .true.
            end if
            if (comdis2.lt.contact_cuts(2)) then
              contact_maps(cnt1,cnt2) = contact_maps(cnt1,cnt2) + scal
              c2 = c2 + 1
            end if
          end if
        end do
      end do
      cnt4 = cnt4 + rsmol(jmol,2) - rsmol(jmol,1) + 1
!     the cluster detection algorithm works generally without additional cost
!     for each molecule-molecule contact (<= res.-res.-contact) it does a number of integer assignment operations
!     the only potentially expensive operation are the searches when joining clusters of size > 1. these could be
!     implemented more efficiently using a cluster-member-list data structure with some overhead cost
!     this is not done here (yet)
      if (doclu.EQV..true.) then
        if (inrange.EQV..true.) then
!         isolated new mol. (add single link, increment size)
          if (inclu(j).le.0) then 
            inclu(j) = inclu(i)
            cluwt(inclu(i),2) = cluwt(inclu(i),2) + 1
            cluwt(inclu(i),1) = cluwt(inclu(i),1) + 1
!         already part of host cluster (just add the link)
          else if (inclu(i).eq.inclu(j)) then
            cluwt(inclu(i),1) = cluwt(inclu(i),1) + 1
!         part of a different cluster -> conjoin smartly
          else
!           host cluster is tiny -> quick-add (single link)
            if (cluwt(inclu(i),2).eq.1) then
              cluwt(inclu(i),:) = 0
              inclu(i) = inclu(j)
              cluwt(inclu(j),2) = cluwt(inclu(j),2) + 1
              cluwt(inclu(j),1) = cluwt(inclu(j),1) + 1
!           otherwise add smaller cluster
            else if (cluwt(inclu(i),2).ge.cluwt(inclu(j),2)) then
              curclu = inclu(j)
              nhits = 0
              do k=1,nsolutes
                if (inclu(k).eq.curclu) then
                  inclu(k) = inclu(i) 
                  nhits = nhits + 1
                end if
                if (nhits.eq.cluwt(inclu(j),2)) exit
              end do
              cluwt(inclu(i),2) = cluwt(inclu(i),2) + cluwt(curclu,2)
              cluwt(inclu(i),1) = cluwt(inclu(i),1) + cluwt(curclu,1) + 1
              cluwt(curclu,:) = 0
            else
              curclu = inclu(i)
              nhits = 0
              do k=1,nsolutes
                if (inclu(k).eq.curclu) then
                  inclu(k) = inclu(j)
                  nhits = nhits + 1
                end if
                if (nhits.eq.cluwt(inclu(i),2)) exit
              end do
              cluwt(inclu(j),2) = cluwt(inclu(j),2) + cluwt(curclu,2)
              cluwt(inclu(j),1) = cluwt(inclu(j),1) + cluwt(curclu,1) + 1
              cluwt(curclu,:) = 0
            end if
          end if
!          for i
!          if (clu_nbl(i)%nnbs.eq.0) then
!            allocate(clu_nbl(i)%nb(clu_nbl(i)%alsz))
!          else
!            if (clu_nbl(i)%nnbs.eq.clu_nbl(i)%alsz) call resz_clunbl(i)
!          end if
!          clu_nbl(i)%nnbs = clu_nbl(i)%nnbs + 1
!          clu_nbl(i)%nb(clu_nbl(i)%nnbs) = j
!!         for j
!          if (clu_nbl(j)%nnbs.eq.0) then
!            allocate(clu_nbl(j)%nb(clu_nbl(j)%alsz))
!            clu_nbl(j)%nnbs = clu_nbl(j)%nnbs + 1
!            clu_nbl(j)%nb(clu_nbl(j)%nnbs) = i
!          else
!            if (clu_nbl(j)%nnbs.eq.clu_nbl(j)%alsz) call resz_clunbl(j)
!          end if
!          clu_nbl(j)%nnbs = clu_nbl(j)%nnbs + 1
!          clu_nbl(j)%nb(clu_nbl(j)%nnbs) = i
        end if
      end if
    end do
  end do
!  if (doinst.EQV..true.) then
! !784 format(i9,1x,1000000(i9,1x,i9,1x,f6.3,1x,f6.3,1x))
!    !write(ibprtnr,784) nstep,(nint(bpdata(j,1)),nint(bpdata(j,2)),bpdata(j,3),bpdata(j,4),j=1,nbpres)
!  784 format(I9,1x,84(I4,1x),84(f6.2,1x))
!    write(ibprtnr,784) nstep,int(bpdata(:,1)),int(bpdata(:,2)),bpdata(:,3),bpdata(:,4)
!    deallocate(bpres)
!    deallocate(bpdata)
!  end if
!
! after collecting raw contact numbers, we can also update the histograms
  if (nsolutes.gt.0) then
    if (c1.ge.(nressolute*20)) then
      write(ilog,*) 'Warning: Too many contacts. Data point &
 & omitted. Increase size of contact_hists.'
    else
      contact_hists(1,c1+1) = contact_hists(1,c1+1) + scal
    end if
    if (c2.ge.(nressolute*20)) then
      write(ilog,*) 'Warning: Too many contacts. Data point&
 & omitted. Increase size of contact_hists.'
    else
      contact_hists(2,c2+1) = contact_hists(2,c2+1) + scal
    end if
  end if
!
! and - if need be - we can parse the arrays for molecular cluster analysis
  if (doclu.EQV..true.) then
    do i=1,nsolutes
      if (cluwt(i,2).gt.0) then
!       normalize by minimum connectivity
        coo = 1.0
        if (cluwt(i,2).gt.1) coo = (1.0*cluwt(i,1))/(1.0*(cluwt(i,2)-1))
!       add the counts to the respective histograms
        do j=1,nsolutes
          if (i.eq.inclu(j)) then
            clumol(j,cluwt(i,2)) = clumol(j,cluwt(i,2)) + scal
          end if
        end do
        clusys(cluwt(i,2)) = clusys(cluwt(i,2)) + scal*cluwt(i,2)
        cobin = max(1,min(floor((coo-1.0)/0.05) + 1,100))
        clucoo(cluwt(i,2),cobin) = clucoo(cluwt(i,2),cobin) + scal
      end if
    end do
    deallocate(cluwt)
    deallocate(inclu)
  end if
!
end
!
!-----------------------------------------------------------------
!
subroutine rescom(rs,com1)
!
  use atoms
  use polypep
!
  implicit none
!
  integer rs,k,kk
  RTYPE com1(3),tmass
! 
! COM 
  com1(:) = 0.0
  tmass = 0.0
  do k=1,at(rs)%nsc+at(rs)%nbb
    if (k.le.at(rs)%nbb) then
      kk = at(rs)%bb(k)
    else
      kk = at(rs)%sc(k-at(rs)%nbb)
    end if
    tmass = tmass + mass(kk)
    com1(1) = com1(1) + mass(kk)*x(kk)
    com1(2) = com1(2) + mass(kk)*y(kk)
    com1(3) = com1(3) + mass(kk)*z(kk)
  end do
  com1(1) = com1(1)/tmass
  com1(2) = com1(2)/tmass
  com1(3) = com1(3)/tmass
!
end
!
!--------------------------------------------------------------------------------------
!
subroutine respmindis(rs,rs2,mindis2)
!
  use polypep
!
  implicit none
!
  integer rs,rs2,l,ll,k,kk
  RTYPE mindis2,svec(3),dis2
!
  mindis2 = huge(mindis2)
! we will find the closest atom-atom distance
  call dis_bound_rs(rs,rs2,svec)
  do k=1,at(rs)%nsc+at(rs)%nbb
    if (k.le.at(rs)%nbb) then
      kk = at(rs)%bb(k)
    else
      kk = at(rs)%sc(k-at(rs)%nbb)
    end if
    do l=1,at(rs2)%nsc+at(rs2)%nbb
      if (l.le.at(rs2)%nbb) then
        ll = at(rs2)%bb(l)
      else
        ll = at(rs2)%sc(l-at(rs2)%nbb)
      end if
      call dis_bound4(kk,ll,svec,dis2)
      if (dis2.lt.mindis2) mindis2 = dis2
    end do
  end do
!
end
!
!--------------------------------------------------------------------------------------
!
subroutine prt_rescontacts()
!
  use iounit
  use sequen
  use contacts
  use mpistuff
  use system
  use molecule
  use aminos
  use pdb
!
  implicit none
!
  integer lasti,i,j,freeunit,iu,ii,jj,imol,cnt1,cnt2,cnt3,rs,rs2
  RTYPE tot1,tot2
  character(60) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  logical exists
  character(3), ALLOCATABLE:: rnams(:)
  integer, ALLOCATABLE:: rnmbs(:)
  RTYPE, ALLOCATABLE:: normer(:,:)
!
! normalize
  allocate(normer(nsolutes,nsolutes))
  if (use_frame_weights.EQV..false.) then
    normer(:,:) = 1.0*ncontactavg(:,:)
  else
    normer(:,:) = 0.0
    j = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),contactcalc).eq.0)) then
        j = j + 1
        normer(:,:) = normer(:,:) + framewts2(i)
      end if
    end do
    if (j.ne.ncontactavg(1,1)) then
      write(ilog,*) 'Warning. Contact maps have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in CONTACTMAP.dat or related files.'
    end if
  end if
  cnt1 = 0
  cnt3 = 0
  cnt2 = 0
  do i=1,nsolutes 
    cnt3 = cnt1
    if (ncontactavg(i,i).gt.0) then
      do rs=rsmol(solutes(i),1),rsmol(solutes(i),2)
        cnt1 = cnt1 + 1
        do rs2=rs+contact_off,rsmol(solutes(i),2)
          if (rs.eq.rs2) then
            contact_maps(cnt1,cnt1) = contact_maps(cnt1,cnt1)/normer(i,i)
          else
            contact_maps(cnt1+rs2-rs,cnt1) = contact_maps(cnt1+rs2-rs,cnt1)/normer(i,i)
            contact_maps(cnt1,cnt1+rs2-rs) = contact_maps(cnt1,cnt1+rs2-rs)/normer(i,i)
          end if
        end do
      end do
    end if
    cnt1 = cnt3
    do rs=rsmol(solutes(i),1),rsmol(solutes(i),2)
      cnt2 = cnt3 + rsmol(solutes(i),2) - rsmol(solutes(i),1) + 1
      cnt1 = cnt1 + 1
      do j=i+1,nsolutes
        if (ncontactavg(i,j).gt.0) then
          do rs2=rsmol(solutes(j),1),rsmol(solutes(j),2)
            cnt2 = cnt2 + 1
            contact_maps(cnt2,cnt1) = contact_maps(cnt2,cnt1)/normer(i,j)
            contact_maps(cnt1,cnt2) = contact_maps(cnt1,cnt2)/normer(i,j)
          end do
        end if
      end do
    end do
  end do
!
! allocate, get annotation arrays
  allocate(rnams(nressolute))
  allocate(rnmbs(nressolute))
  j = 0
  do imol=1,nmol
    if (is_solvent(imol).EQV..false.) then
      do i=rsmol(imol,1),rsmol(imol,2)
        j = j + 1
        rnams(j) = amino(seqtyp(i))
        rnmbs(j) = i
      end do
    end if
  end do
!
 45   format(2000000(i8,': ',a3,' |'))
 46   format(2000000(g14.7,1x))
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_CONTACTMAP.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'CONTACTMAP.dat'
  end if
#else
  fn = 'CONTACTMAP.dat'
#endif 
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  write(iu,45) (rnmbs(j),rnams(j),j=1,nressolute)
  do i=1,nressolute
    write(iu,46) (contact_maps(i,j),j=1,nressolute)
  end do
  close(unit=iu)
!     
 47   format(i10,1x,g14.7,1x,g14.7)
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_CONTACT_HISTOGRAMS.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'CONTACT_HISTOGRAMS.dat'
  end if
#else
  fn = 'CONTACT_HISTOGRAMS.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  lasti = 0
  do i=nressolute*20,1,-1
    if ((contact_hists(1,i).gt.0.0).OR.(contact_hists(2,i).gt.0.0))&
 &then
      lasti = i
      exit
    end if
  end do
  tot1 = sum(contact_hists(1,1:lasti))
  tot2 = sum(contact_hists(2,1:lasti))
  do i=1,lasti
    write(iu,47) i-1,contact_hists(1,i)/tot1,contact_hists(2,i)/tot2
  end do
  close(unit=iu)
!
! finally, allow print-out of cluster distributions as well
  if (contactcalc*clucalc.le.nsim) then
    call prt_molclusters()
  end if
!
  deallocate(rnams)
  deallocate(rnmbs)
  deallocate(normer)  
!
end
!
!----------------------------------------------------------------------------
!
subroutine resz_clunbl(isol)
!
  use contacts
!
  implicit none
!
  integer isol,bua(clu_nbl(isol)%nnbs),i
!
  do i=1,clu_nbl(isol)%nnbs
    bua(i) = clu_nbl(isol)%nb(i)
  end do
!
  clu_nbl(isol)%alsz = clu_nbl(isol)%alsz*4
  deallocate(clu_nbl(isol)%nb)
  allocate(clu_nbl(isol)%nb(clu_nbl(isol)%alsz))
!
  do i=1,clu_nbl(isol)%nnbs
    clu_nbl(isol)%nb(i) = bua(i)
  end do
!
end
!
!-----------------------------------------------------------------
!
subroutine prt_molclusters()
!
  use iounit
  use molecule
  use contacts
  use mpistuff
!
  implicit none
!
  integer i,j,freeunit,iu,maxclsz,ii,jj
  RTYPE tot
  character(60) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  logical exists
!
! first output for the overall cluster distribution
  tot = 0.0
  maxclsz = 0
  do i=1,nsolutes
    if ((clusys(i).gt.0).AND.(i.gt.maxclsz)) maxclsz = i
    tot = tot + clusys(i)
  end do
  if (tot.gt.0.0) clusys(1:maxclsz) = clusys(1:maxclsz)/tot
!
 56   format(i6,1x,g14.7)
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_CLUSTERS.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'CLUSTERS.dat'
  end if
#else
  fn = 'CLUSTERS.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  do i=1,nsolutes
    write(iu,56) i,clusys(i)
  end do
  close(unit=iu)
!
  do i=1,nsolutes
    tot = sum(clumol(i,1:maxclsz))
    if (tot.gt.0.0) clumol(i,1:maxclsz) = clumol(i,1:maxclsz)/tot
  end do
!
 57   format(i10,2000000(1x,g14.7))
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_MOLCLUSTERS.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'MOLCLUSTERS.dat'
  end if
#else
  fn = 'MOLCLUSTERS.dat'
#endif 
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  do i=1,nsolutes
    write(iu,57) solutes(i),(clumol(i,j),j=1,maxclsz)
  end do
  close(unit=iu)
!
  do i=1,maxclsz
    tot = sum(clucoo(i,:))
    if (tot.gt.0.0) clucoo(i,:) = clucoo(i,:)/tot
  end do
!
 58   format(g12.5,2000000(1x,g14.7))
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_COOCLUSTERS.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'COOCLUSTERS.dat'
  end if
#else
  fn = 'COOCLUSTERS.dat'
#endif
  call strlims(fn,ii,jj) 
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  do j=1,100
    write(iu,58) 1.0+(j-0.5)*0.05,(clucoo(i,j),i=1,maxclsz)
  end do
  close(unit=iu)
!
end
!
!
!----------------------------------------------------------------------------------
!
! ROUTINES NEEDED TO PERFORM PAIR CORRELATION ANALYSIS (DISTANCE DISTRIBUTIONS)
!
!----------------------------------------------------------------------------------
!
subroutine do_general_pc()
!
  use atoms
  use sequen
  use molecule
  use paircorr
  use system
  use pdb
  use grandensembles
  use mcsums
!
  implicit none
!
  integer i,pcbin,rs1,rs2,ati,atj,kk
  integer imol,jmol,kmol,lmol,imt,jmt,cnt,cnt2
  RTYPE dis,dis2,svec(3),incr
  logical ismember,doinst
  RTYPE, ALLOCATABLE:: disvec(:)
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  n_pccalls = n_pccalls + 1
  doinst = .false.
  kk = 0
  if (inst_gpc.gt.0) then
    if (mod(n_pccalls,inst_gpc).eq.0) doinst = .true.
  end if
  if (doinst.EQV..true.) then
    if (n_inst_gpc.le.0) call inst_gpc_header()
    allocate(disvec(n_inst_gpc))
    disvec(:) = -1.0
  end if
!
! note we're assuming do_rbc_pc was called prior to this and we have an up-to-date voli_pc
!
  do i=1,gpc%nlst
!   for specific pairs it's trivial
    if (gpc%lty(i).EQV..false.) then
      rs1 = atmres(gpc%lst(i,1))
      rs2 = atmres(gpc%lst(i,2))
      imol = molofrs(atmres(gpc%lst(i,1)))
      jmol = molofrs(atmres(gpc%lst(i,2)))
      kk = kk + 1
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
          if (ismember(ispresent,imol).EQV..false.) cycle
        end if
        if (ismember(fluctypes,moltypid(jmol)).EQV..true.) then
          if (ismember(ispresent,jmol).EQV..false.) cycle
        end if
      end if
      gpc%navg(i) = gpc%navg(i) + 1
      if (imol.eq.jmol) then
        dis2 = (x(gpc%lst(i,1)) - x(gpc%lst(i,2)))**2 + &
 &             (y(gpc%lst(i,1)) - y(gpc%lst(i,2)))**2 + &
 &             (z(gpc%lst(i,1)) - z(gpc%lst(i,2)))**2
        if (dis2.le.0.0) then
          dis = 0.0
        else
          dis = sqrt(dis2)
        end if
        if (doinst.EQV..true.) disvec(kk) = dis
        pcbin = max(1,floor(dis/pc_binsz) + 1)
!       we will mercilessly overstock the last bin in case of range exceptions
        if (pcbin.gt.MAXPCBINS) then
          gpc%pc(gpc%lst(i,3),MAXPCBINS) = gpc%pc(gpc%lst(i,3),MAXPCBINS) + incr
        else
          gpc%pc(gpc%lst(i,3),pcbin) = gpc%pc(gpc%lst(i,3),pcbin) + incr
        end if
      else
        call dis_bound_rs(rs1,rs2,svec)
        call dis_bound4(gpc%lst(i,1),gpc%lst(i,2),svec,dis2)
        if (doinst.EQV..true.) disvec(kk) = sqrt(dis2)
        pcbin = max(1,floor(sqrt(dis2)/pc_binsz) + 1)
!       we will mercilessly overstock the last bin in case of range exceptions
        if (pcbin.gt.MAXPCBINS) then
          gpc%pc(gpc%lst(i,3),MAXPCBINS) = gpc%pc(gpc%lst(i,3),MAXPCBINS) + incr/voli_pc(pcbin)
        else
          gpc%pc(gpc%lst(i,3),pcbin) = gpc%pc(gpc%lst(i,3),pcbin) + incr/voli_pc(pcbin)
        end if
      end if
    else
!   for type-generalized pairs it's a little more work
      imol = molofrs(atmres(gpc%lst(i,1)))
      jmol = molofrs(atmres(gpc%lst(i,2)))
      imt = an_grp_mol(imol)
      jmt = an_grp_mol(jmol)
!     intramolecular terms are relatively simple
      if (imol.eq.jmol) then
        cnt = 0
        do kmol=molangr(imt,1),nmol
          if (an_grp_mol(kmol).eq.imt) then
            cnt = cnt + 1
            kk = kk + 1
            if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
              if (ismember(fluctypes,imt).EQV..true.) then
                if (ismember(ispresent,kmol).EQV..false.) then
                  if (cnt.eq.molangr(imt,2)) exit
                  cycle
                end if
              end if
            end if
            gpc%navg(i) = gpc%navg(i) + 1
            ati = atmol(kmol,1)+gpc%lst(i,1)-atmol(imol,1)
            atj = atmol(kmol,1)+gpc%lst(i,2)-atmol(imol,1)
            dis2 = (x(ati) - x(atj))**2 + (y(ati) - y(atj))**2 + (z(ati) - z(atj))**2
            if (dis2.le.0.0) then
              dis = 0.0
            else
              dis = sqrt(dis2)
            end if
            if (doinst.EQV..true.) disvec(kk) = dis
            pcbin = max(1,floor(dis/pc_binsz) + 1)
            if (pcbin.gt.MAXPCBINS) then
              gpc%pc(gpc%lst(i,3),MAXPCBINS) = gpc%pc(gpc%lst(i,3),MAXPCBINS) + incr
            else
              gpc%pc(gpc%lst(i,3),pcbin) = gpc%pc(gpc%lst(i,3),pcbin) + incr
            end if
            if (cnt.eq.molangr(imt,2)) exit
          end if
        end do
!     now intermolecular terms for molecules of the same type
      else if (imt.eq.jmt) then
        cnt = 0
        cnt2 = 0
        do kmol=molangr(imt,1),nmol
          if (an_grp_mol(kmol).ne.imt) cycle
          cnt2 = cnt2 + 1
          if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
            if (ismember(fluctypes,imt).EQV..true.) then
              if (ismember(ispresent,kmol).EQV..false.) then
                kk = kk + molangr(jmt,2) - cnt2
                if ((gpc%lst(i,1)-atmol(imol,1)).ne.(gpc%lst(i,2)-atmol(jmol,1))) kk = kk + molangr(jmt,2) - cnt2
                if (cnt.eq.molangr(imt,2)) exit
                cycle
              end if
            end if
          end if
          do lmol=kmol+1,nmol
            if (an_grp_mol(lmol).eq.imt) then
              kk = kk + 1
              if ((gpc%lst(i,1)-atmol(imol,1)).ne.(gpc%lst(i,2)-atmol(jmol,1))) kk = kk + 1
              if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
                if (ismember(fluctypes,imt).EQV..true.) then
                  if (ismember(ispresent,lmol).EQV..false.) cycle
                end if
              end if
              gpc%navg(i) = gpc%navg(i) + 1
              ati = atmol(kmol,1)+gpc%lst(i,1)-atmol(imol,1)
              atj = atmol(lmol,1)+gpc%lst(i,2)-atmol(jmol,1)
              rs1 = atmres(ati)
              rs2 = atmres(atj)
              call dis_bound_rs(rs1,rs2,svec)
              call dis_bound4(ati,atj,svec,dis2)
              if (doinst.EQV..true.) disvec(kk) = sqrt(dis2)
              pcbin = max(1,floor(sqrt(dis2)/pc_binsz) + 1)
              if (pcbin.gt.MAXPCBINS) then
                gpc%pc(gpc%lst(i,3),MAXPCBINS) = gpc%pc(gpc%lst(i,3),MAXPCBINS) + incr/voli_pc(pcbin)
              else
                gpc%pc(gpc%lst(i,3),pcbin) = gpc%pc(gpc%lst(i,3),pcbin) + incr/voli_pc(pcbin)
              end if
              if ((gpc%lst(i,1)-atmol(imol,1)).ne.(gpc%lst(i,2)-atmol(jmol,1))) then
!               off-diagonal elements have a two-fold degeneracy due to
!               the symmetry of the matrix
                gpc%navg(i) = gpc%navg(i) + 1
                ati = atmol(kmol,1)+gpc%lst(i,2)-atmol(jmol,1)
                atj = atmol(lmol,1)+gpc%lst(i,1)-atmol(imol,1)
                rs1 = atmres(ati)
                rs2 = atmres(atj)
                call dis_bound_rs(rs1,rs2,svec)
                call dis_bound4(ati,atj,svec,dis2)
                if (doinst.EQV..true.) disvec(kk) = sqrt(dis2)
                pcbin = max(1,floor(sqrt(dis2)/pc_binsz) + 1)
                if (pcbin.gt.MAXPCBINS) then
                  gpc%pc(gpc%lst(i,3),MAXPCBINS) = gpc%pc(gpc%lst(i,3),MAXPCBINS) + incr/voli_pc(pcbin)
                else
                  gpc%pc(gpc%lst(i,3),pcbin) = gpc%pc(gpc%lst(i,3),pcbin) + incr/voli_pc(pcbin)
                end if
              end if
            end if
          end do
          cnt = cnt + 1
          if (cnt.eq.molangr(imt,2)) exit
        end do
!     and finally intermolecular terms for molecules of different types
      else
        cnt = 0
        do kmol=molangr(imt,1),nmol
          if (an_grp_mol(kmol).ne.imt) cycle
          if (ismember(fluctypes,imt).EQV..true.) then
            if (ismember(ispresent,kmol).EQV..false.) then
              kk = kk + molangr(jmt,2)
              if (cnt.eq.molangr(imt,2)) exit
              cycle
            end if
          end if
          cnt2 = 0
          do lmol=molangr(jmt,1),nmol
            if (an_grp_mol(lmol).eq.jmt) then
              kk = kk + 1
              if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
                if (ismember(fluctypes,jmt).EQV..true.) then
                  if (ismember(ispresent,lmol).EQV..false.) then
                    if (cnt2.eq.molangr(jmt,2)) exit
                    cycle
                  end if
                end if
              end if
              gpc%navg(i) = gpc%navg(i) + 1
              ati = atmol(kmol,1)+gpc%lst(i,1)-atmol(imol,1)
              atj = atmol(lmol,1)+gpc%lst(i,2)-atmol(jmol,1)
              rs1 = atmres(ati)
              rs2 = atmres(atj)
              call dis_bound_rs(rs1,rs2,svec)
              call dis_bound4(ati,atj,svec,dis2)
              if (doinst.EQV..true.) disvec(kk) = sqrt(dis2)
              pcbin = max(1,floor(sqrt(dis2)/pc_binsz) + 1)
              if (pcbin.gt.MAXPCBINS) then
                gpc%pc(gpc%lst(i,3),MAXPCBINS) = gpc%pc(gpc%lst(i,3),MAXPCBINS) + incr/voli_pc(pcbin)
              else
                gpc%pc(gpc%lst(i,3),pcbin) = gpc%pc(gpc%lst(i,3),pcbin) + incr/voli_pc(pcbin)
              end if
              cnt2 = cnt2 + 1
            end if
            if (cnt2.eq.molangr(jmt,2)) exit
          end do
          if (cnt.eq.molangr(imt,2)) exit
        end do
      end if ! what type of type-generalization
    end if !whether or not type-generalized
  end do
!
  if (doinst.EQV..true.) then
 666 format(i11,1x,100000000(g11.5,1x))
    write(igpc,666) nstep,disvec(:)
    deallocate(disvec)
  end if
!
end
!
!---------------------------------------------------------------------
!
subroutine inst_gpc_header()
!
  use paircorr
  use mcsums
  use system
  use molecule
  use sequen
  use atoms
!
  implicit none
!
  integer i,imol,jmol,kmol,lmol,imt,jmt,cnt,cnt2,ati(3),kk,strl
  character(12), ALLOCATABLE:: headstr(:)
  character(8) sa1
!
  if ((inst_gpc.le.0).OR.(inst_gpc.gt.(nsim/pccalc))) return
!
  strl = 1
  sa1(:) = ' '
 555 format(100000001(a12))
!
  do kk=1,3
    n_inst_gpc = 0
!
    do i=1,gpc%nlst
!     for specific pairs it's trivial
      if (gpc%lty(i).EQV..false.) then
        ati(2) = gpc%lst(i,1)
        ati(3) = gpc%lst(i,2)
        n_inst_gpc = n_inst_gpc + 1
        if (kk.gt.1) then
          call int2str(ati(kk),sa1,strl)
          headstr(n_inst_gpc) = '| '//sa1//' |'
          strl = 1
        end if
      else
!     for type-generalized pairs it's a little more work
        imol = molofrs(atmres(gpc%lst(i,1)))
        jmol = molofrs(atmres(gpc%lst(i,2)))
        imt = moltypid(imol)
        jmt = moltypid(jmol)
!       intramolecular terms are relatively simple
        if (imol.eq.jmol) then
          cnt = 0
          do kmol=moltyp(imt,1),nmol
            if (moltypid(kmol).eq.imt) then
              cnt = cnt + 1
              ati(2) = atmol(kmol,1)+gpc%lst(i,1)-atmol(imol,1)
              ati(3) = atmol(kmol,1)+gpc%lst(i,2)-atmol(imol,1)
              n_inst_gpc = n_inst_gpc + 1
              if (kk.gt.1) then
                call int2str(ati(kk),sa1,strl)
                headstr(n_inst_gpc) = '| '//sa1//' |'
                strl = 1
              end if
              if (cnt.eq.moltyp(imt,2)) exit
            end if
          end do
!       now intermolecular terms for molecules of the same type
        else if (imt.eq.jmt) then
          cnt = 0
          do kmol=moltyp(imt,1),nmol
            if (moltypid(kmol).ne.imt) cycle
            do lmol=kmol+1,nmol
              if (moltypid(lmol).eq.imt) then
                ati(2) = atmol(kmol,1)+gpc%lst(i,1)-atmol(imol,1)
                ati(3) = atmol(lmol,1)+gpc%lst(i,2)-atmol(jmol,1)
                n_inst_gpc = n_inst_gpc + 1
                if (kk.gt.1) then
                  call int2str(ati(kk),sa1,strl)
                  headstr(n_inst_gpc) = '| '//sa1//' |'
                  strl = 1
                end if
                if ((gpc%lst(i,1)-atmol(imol,1)).ne.(gpc%lst(i,2)-atmol(jmol,1))) then
                  ati(2) = atmol(kmol,1)+gpc%lst(i,2)-atmol(jmol,1)
                  ati(3) = atmol(lmol,1)+gpc%lst(i,1)-atmol(imol,1)
                  n_inst_gpc = n_inst_gpc + 1
                  if (kk.gt.1) then
                    call int2str(ati(kk),sa1,strl)
                    headstr(n_inst_gpc) = '| '//sa1//' |'
                    strl = 1
                  end if
                end if
              end if
            end do
            cnt = cnt + 1
            if (cnt.eq.moltyp(imt,2)) exit
          end do
!       and finally intermolecular terms for molecules of different types
        else
          cnt = 0
          do kmol=moltyp(imt,1),nmol
            if (moltypid(kmol).ne.imt) cycle
            cnt2 = 0
            do lmol=moltyp(jmt,1),nmol
              if (moltypid(lmol).eq.jmt) then
                ati(2) = atmol(kmol,1)+gpc%lst(i,1)-atmol(imol,1)
                ati(3) = atmol(lmol,1)+gpc%lst(i,2)-atmol(jmol,1)
                n_inst_gpc = n_inst_gpc + 1
                if (kk.gt.1) then
                  call int2str(ati(kk),sa1,strl)
                  headstr(n_inst_gpc) = '| '//sa1//' |'
                  strl = 1
                end if
                cnt2 = cnt2 + 1
              end if
              if (cnt2.eq.moltyp(jmt,2)) exit
            end do
            if (cnt.eq.moltyp(imt,2)) exit
          end do
        end if ! what type of type-generalization
      end if !whether or not type-generalized
    end do  
    if (kk.eq.1) then
      allocate(headstr(n_inst_gpc))
    else if (kk.ge.2) then
      write(igpc,555) '#           ',headstr
    end if
    if (kk.eq.3) deallocate(headstr)
!
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine prt_general_pc()
!
  use iounit
  use paircorr
  use mpistuff
  use sequen
  use atoms
  use pdb
  use system
!
  implicit none
!
  integer i,j,maxb,iu,freeunit,ii,jj,imol,jmol
  integer(8) k
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  character(60) fn
  logical exists,prtwarn
  RTYPE tt,netcnt
  RTYPE, ALLOCATABLE:: normer(:)
  integer(8), ALLOCATABLE:: aider(:)
!
  maxb = -1
  do i=MAXPCBINS,1,-1
    do j=1,gpc%nos
      if (gpc%pc(j,i).gt.0.0) then
        maxb = i 
        exit
      end if
    end do
    if (maxb.ne.-1) exit
  end do
!
  if (maxb.eq.-1) then
    write(ilog,*) 'Warning: General PC analysis yielded all zeros. O&
 &mitting write-out of data (frequency setting?).'
    return
  end if
!
  allocate(normer(gpc%nlst))
  if (use_frame_weights.EQV..true.) then
    prtwarn = .false.
    allocate(aider(gpc%nlst))
    aider(:) = 1
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),pccalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)
      end if
    end do
    do i=1,gpc%nlst
      aider(i) = gpc%navg(i)/k
      if (mod(gpc%navg(i),k).ne.0) prtwarn = .true.
    end do
    normer(:)= normer(:)*aider(:)
    if (prtwarn.EQV..true.) then
      write(ilog,*) 'Warning. General pair correlation functions have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in GENERAL_PC.dat or analogous files.'
    end if
    deallocate(aider)
  else
    normer(:) = 1.0*gpc%navg(1:gpc%nlst)
  end if
!
  do j=1,gpc%nos
    netcnt = 0.0
    tt = sum(gpc%pc(j,1:maxb))
    if (tt.le.0.0) cycle
    do i=1,gpc%nlst
      if (gpc%lst(i,3).eq.j) then
        imol = molofrs(atmres(gpc%lst(i,1)))
        jmol = molofrs(atmres(gpc%lst(i,2)))
        if (jmol.eq.imol) then
          exit
        else
          netcnt = netcnt + normer(i)
        end if
      end if
    end do
    if (imol.eq.jmol) then
!     intramolecular 
      do i=1,maxb
        gpc%pc(j,i) = gpc%pc(j,i)/tt
      end do
    else
!     intermolecular
      do i=1,maxb
        gpc%pc(j,i) = gpc%pc(j,i)/netcnt
      end do
    end if
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_GENERAL_PC.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'GENERAL_PC.dat'
  end if
#else
  fn = 'GENERAL_PC.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
 48   format(f11.4,1x,100000000(g11.5,1x))
  do i=1,maxb
    write(iu,48) (i-0.5)*pc_binsz,(gpc%pc(j,i),j=1,gpc%nos)
  end do
  close(unit=iu)
!
  deallocate(normer)
!
end
!
!---------------------------------------------------------------------
!
! this generic PC calculator for molecular centroids can be very useful in simulations
! of small molecules to get a generic idea of preferential interactions in the system.
! since it computes only intermolecular distances, it does observe periodic boundary
! conditions
!
subroutine do_rbc_pc()
!
  use molecule
  use paircorr
  use system
  use math
  use grandensembles
  use movesets
  use pdb
!
  implicit none
!
  integer i,j,k,putit,imol,jmol,imt,jmt
  RTYPE comi(3),comj(3),dis2,incr
  logical ismember
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! re-compute volume element (note general_pc may rely on this)
    call get_voli_pc()
  end if
!
  if (nsolutes.le.1) return
!
   if (use_dyn.EQV..true.) then
    if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      if (in_dyncyc.EQV..false.) then
        do i=1,nsolutes
          imol = solutes(i)
          call update_rigidm(imol)
        end do
      end if
    end if
  end if
!
  do i=1,nsolutes
    imol = solutes(i)
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    imt = an_grp_mol(imol)
    do j=i+1,nsolutes
      jmol = solutes(j)
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(jmol)).EQV..true.) then
          if (ismember(ispresent,jmol).EQV..false.) cycle
        end if
      end if
      jmt = an_grp_mol(jmol)
      if (imt.lt.jmt) then
        npcavg(imt,jmt) = npcavg(imt,jmt) + 1
      else
        npcavg(jmt,imt) = npcavg(jmt,imt) + 1
      end if
      do k=1,3
        if (use_dyn.EQV..true.) then
          comi(k) = comm(imol,k)
          comj(k) = comm(jmol,k)
        else
          comi(k) = com(imol,k)
          comj(k) = com(jmol,k)
        end if
      end do
      call dis_bound2(comi,comj,dis2)
      putit = max(1,floor(sqrt(dis2)/pc_binsz) + 1)
      if (putit.gt.MAXPCBINS) putit = MAXPCBINS
      if (jmt.gt.imt) then
        rbc_pc(imt,jmt,putit) = rbc_pc(imt,jmt,putit) + incr/voli_pc(putit)
      else
        rbc_pc(jmt,imt,putit) = rbc_pc(jmt,imt,putit) + incr/voli_pc(putit)
      end if
    end do
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine get_voli_pc()
!
  use molecule
  use paircorr
  use system
  use math
  use iounit
!
  implicit none
!
  integer i,j,k,l
  RTYPE bndr3,iva,iva2,iva3,iva4,iva6,ivb,ivb2,ivb3,ivb4,ivb6,bndr6t,bndr3t,bndr2t,bndr0t
  RTYPE ivav(3),ivbv(3),ivav2(3),ivbv2(3),ivav3(3),ivbv3(3),volt,minb,maxb,pi43
  RTYPE termx1,termx2,termk1,termk2,terml1,terml2,termkl1,termkl2
!
  pi43 = (4./3.)*pi
!
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
!     PBC with rectangular boxes (messy)
      minb = bnd_params(1)
      maxb = bnd_params(1)
      do i=2,3
        if (bnd_params(i).lt.minb) minb = bnd_params(i)
        if (bnd_params(i).gt.maxb) maxb = bnd_params(i)
      end do
      minb = 0.5*minb
      maxb = 0.5*sqrt(sum(bnd_params(1:3)*bnd_params(1:3)))
      volt = bnd_params(1)*bnd_params(2)*bnd_params(3)
      do i=1,MAXPCBINS
        iva = (i-1)*pc_binsz
        ivb = (i)*pc_binsz
        if (ivb.le.minb) then
          iva3 = iva*iva*iva
          ivb3 = ivb*ivb*ivb
          voli_pc(i) = pi43*(ivb3-iva3)
        else
          if ((ivb.gt.maxb).OR.(i.eq.MAXPCBINS)) ivb = maxb
          iva3 = iva*iva*iva
          ivb3 = ivb*ivb*ivb
          voli_pc(i) = pi43*(ivb3-iva3)
          ivbv = 0.5*bnd_params(1:3)/ivb
          ivav = 0.5*bnd_params(1:3)/iva
          ivav2(1:3) = ivav(1:3)*ivav(1:3)
          ivbv2(1:3) = ivbv(1:3)*ivbv(1:3)
          ivav3(1:3) = ivav(1:3)*ivav2(1:3)
          ivbv3(1:3) = ivbv(1:3)*ivbv2(1:3)
!         subtract spherical caps in all three directions (2 each)
          do k=1,3
            if ((ivbv(k).lt.1.0).AND.(ivav(k).lt.1.0)) then
!              voli_pc(i) = voli_pc(i) - ((0.5*(iva+ivb))**3)*0.5*pi43*(2.0 + (ivbv3(k) - ivav3(k)) + 3.0*(ivbv(k)-ivav(k)))
              voli_pc(i) = voli_pc(i) - (ivb**3)*0.5*pi43*(2.0 + ivbv3(k) - 3.0*ivbv(k)) + &
 &                                      (iva**3)*0.5*pi43*(2.0 + ivav3(k) - 3.0*ivav(k))
            else if (ivbv(k).lt.1.0) then
               voli_pc(i) = voli_pc(i) - (ivb**3)*0.5*pi43*(2.0 + ivbv3(k) - 3.0*ivbv(k))
            end if
          end do
!         add back spherical wedges in all three pairs of directions
          do k=1,3
            do l=k+1,3     
              if (((ivbv2(k)+ivbv2(l)).lt.1.0).AND.((ivav2(k)+ivav2(l)).lt.1.0)) then
                termk1 = -0.5*(3.0*ivav(k) - ivav3(k))*acos(ivav(l)/sqrt(1.0-ivav2(k)))
                termk2 = -0.5*(3.0*ivbv(k) - ivbv3(k))*acos(ivbv(l)/sqrt(1.0-ivbv2(k)))
                terml1 = -0.5*(3.0*ivav(l) - ivav3(l))*acos(ivav(k)/sqrt(1.0-ivav2(l)))
                terml2 = -0.5*(3.0*ivbv(l) - ivbv3(l))*acos(ivbv(k)/sqrt(1.0-ivbv2(l)))
                termkl1 = acos(ivav(l)*ivav(k)/((sqrt(1.0-ivav2(k)))*(sqrt(1.0-ivav2(l)))))
                termkl2 = acos(ivbv(l)*ivbv(k)/((sqrt(1.0-ivbv2(k)))*(sqrt(1.0-ivbv2(l)))))
                termx1 = ivav(l)*ivav(k)*sqrt(1.0 - ivav2(k) - ivav2(l)) 
                termx2 = ivbv(l)*ivbv(k)*sqrt(1.0 - ivbv2(k) - ivbv2(l))
                voli_pc(i) = voli_pc(i) + (ivb**3)*(8./3.)*(termkl2 + termk2 + terml2 + termx2) - &
     &                                    (iva**3)*(8./3.)*(termkl1 + termk1 + terml1 + termx1)
              else if ((ivbv2(k)+ivbv2(l)).lt.1.0) then
                termk2 = -0.5*(3.0*ivbv(k) - ivbv3(k))*acos(ivbv(l)/sqrt(1.0-ivbv2(k)))
                terml2 = -0.5*(3.0*ivbv(l) - ivbv3(l))*acos(ivbv(k)/sqrt(1.0-ivbv2(l)))
                termkl2 = acos(ivbv(l)*ivbv(k)/((sqrt(1.0-ivbv2(k)))*(sqrt(1.0-ivbv2(l)))))
                termx2 = ivbv(l)*ivbv(k)*sqrt(1.0 - ivbv2(k) - ivbv2(l))
                voli_pc(i) = voli_pc(i) + (ivb**3)*(8./3.)*(termkl2 + termk2 + terml2 + termx2)
              end if
            end do
          end do
        end if
        voli_pc(i) = voli_pc(i)/volt
        if (i*pc_binsz.ge.maxb) exit
      end do
      do j=1,MAXPCBINS
        if (voli_pc(j).le.0.0) exit
      end do
      if (j.lt.MAXPCBINS) then
        do i=j,MAXPCBINS
          voli_pc(i) = voli_pc(max(1,j-1))
        end do
      end if
    else
      write(ilog,*) 'Warning. Analytical volume element for pair correlation data not yet available for chosen system container &
 &(in get_voli_pc()). Using flat element (will print out raw distance distributions for all cases).'
      voli_pc(:) = 1.0
    end if
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   wall-type BC with spherical walls
    if (bnd_shape.eq.2) then
      bndr3 = bnd_params(5)*bnd_params(4)
      bndr6t = 0.1875/(bndr3*bndr3)
      bndr3t = (16.0/3.0)*bndr3
      bndr2t = 3.0*bnd_params(5)
      bndr0t = (1./6.)
      do i=1,MAXPCBINS
        iva = (i-1)*pc_binsz
        iva2 = iva*iva
        iva4 = iva2*iva2
        iva6 = iva4*iva2
        ivb = (i)*pc_binsz
        if ((ivb.gt.2.0*bnd_params(4)).OR.(i.eq.MAXPCBINS)) ivb = 2.0*bnd_params(4)
        ivb2 = ivb*ivb
        ivb4 = ivb2*ivb2
        ivb6 = ivb2*ivb4
        voli_pc(i) = bndr6t*(bndr0t*(ivb6-iva6) - bndr2t*(ivb4-iva4) + bndr3t*(ivb2*ivb-iva2*iva))
        if (i*pc_binsz.gt.2.0*bnd_params(4)) then
          exit
        end if
      end do
      do j=1,MAXPCBINS
        if (voli_pc(j).le.0.0) exit
      end do
      if (j.lt.MAXPCBINS) then
        do i=j,MAXPCBINS
          voli_pc(i) = voli_pc(max(1,j-1))
        end do
      end if
    else
      write(ilog,*) 'Warning. Analytical volume element for pair correlation data not yet available for chosen system container &
 &(in get_voli_pc()). Using flat element (will print out raw distance distributions for all cases).'
      voli_pc(:) = 1.0
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in get_voli_pc() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------
!
subroutine prt_rbc_pc()
!
  use molecule
  use paircorr
  use mpistuff
  use system
  use pdb
  use iounit
!
  implicit none
!
  integer i,k,freeunit,iu,howmany,maxk,imt,jmt,ii,jj
  RTYPE, ALLOCATABLE:: vals(:)
  integer, ALLOCATABLE:: tys(:,:),aider(:,:)
  logical, ALLOCATABLE:: prtmat(:,:)
  RTYPE pcb,bndr0t
  character(60) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  logical exists,supvol
  RTYPE, ALLOCATABLE:: normer(:,:)
!
  if (nsolutes.le.1) return
!
  allocate(vals(nangrps*nangrps))
  allocate(tys(nangrps*nangrps,2))
  allocate(prtmat(nangrps,nangrps))
!
! init
  do imt=1,nangrps
    do jmt=imt,nangrps
      prtmat(imt,jmt) = .false.
    end do
  end do
  maxk = 0
!
  allocate(normer(nangrps,nangrps))
  if (use_frame_weights.EQV..true.) then
    allocate(aider(nangrps,nangrps))
    normer(:,:) = 0.0
    k = 0
    do imt=1,nangrps
      do jmt=1,nangrps
        if (imt.eq.jmt) then
          aider(imt,jmt) = molangr(imt,2)*(molangr(imt,2)-1)/2
        else
          aider(imt,jmt) =  molangr(imt,2)*molangr(jmt,2)
        end if
      end do
    end do
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),pccalc).eq.0)) then
        k = k + 1
        normer(:,:) = normer(:,:) + framewts2(i)*aider(:,:)
      end if
    end do
    if ((molangr(1,2)*(molangr(1,2)-1)/2).ne.(npcavg(1,1)/k)) then
      write(ilog,*) 'Warning. Rigid-body pair correlation functions have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in RBC_PC.dat or analogous files.'
    end if
    deallocate(aider)
  else
    normer(:,:) = 1.0*npcavg(:,:)
  end if
!
! normalize the raw counts
  do imt=1,nangrps
    do jmt=imt,nangrps
      jj = 0
      do k=1,MAXPCBINS
        if (rbc_pc(imt,jmt,k).gt.0.0) jj = k
        if ((rbc_pc(imt,jmt,k).gt.0.0).AND.(k.gt.maxk)) maxk = k
      end do
      if (jj.gt.0) then
        prtmat(imt,jmt) = .true.
        rbc_pc(imt,jmt,:) = rbc_pc(imt,jmt,:)/normer(imt,jmt)
      end if
    end do
  end do
!
  howmany = 0
  do imt=1,nangrps
    do jmt=imt,nangrps
      if (prtmat(imt,jmt).EQV..true.) then
        howmany = howmany + 1
        tys(howmany,1) = imt
        tys(howmany,2) = jmt
      end if
    end do
  end do
!
  if (howmany.eq.0) then
    deallocate(vals)
    deallocate(tys)
    deallocate(prtmat)
    deallocate(normer)
    return
  end if
!
 67   format('#    Distance |',' Vol. Element |',4000000('|',i3,' vs',i3,'    |'))
 66   format(f14.7,1x,4000001(g14.7,1x))
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_RBC_PC.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'RBC_PC.dat'
  end if
#else
  fn = 'RBC_PC.dat'
#endif
  call strlims(fn,ii,jj) 
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
!
  howmany = 0
  do imt=1,nangrps
    do jmt=imt,nangrps
      if (prtmat(imt,jmt).EQV..true.) then
        howmany = howmany + 1
        tys(howmany,1) = imt
        tys(howmany,2) = jmt
      end if
    end do
  end do
  write(iu,67) (tys(i,1),tys(i,2),i=1,howmany)
  supvol = .false.
  do k=1,MAXPCBINS!maxk 
    howmany = 0
    pcb = (k-0.5)*pc_binsz
    if (bnd_type.eq.1) then
      if (bnd_shape.eq.1) then
        bndr0t = 0.5*sqrt(sum(bnd_params(1:3)*bnd_params(1:3)))
        if ((k.gt.maxk).AND.((k-1)*pc_binsz.ge.bndr0t)) exit
        if ((k-1)*pc_binsz.ge.bndr0t) supvol = .true.
      else if (bnd_shape.eq.3) then
        bndr0t = sqrt(bnd_params(5)+0.25*bnd_params(6)**2)
        if ((k.gt.maxk).AND.((k-1)*pc_binsz.ge.bndr0t)) exit
        if ((k-1)*pc_binsz.ge.bndr0t) supvol = .true.
      end if
    else if ((bnd_type.eq.2).OR.(bnd_type.eq.3).OR.(bnd_type.eq.4)) then
      if (bnd_shape.eq.1) then
        bndr0t = sqrt(sum(bnd_params(1:3)*bnd_params(1:3)))
        if ((k.gt.maxk).AND.((k-1)*pc_binsz.ge.bndr0t)) exit
        if ((k-1)*pc_binsz.ge.bndr0t) supvol = .true.
      else if (bnd_shape.eq.2) then
        if ((k.gt.maxk).AND.((k-1)*pc_binsz.ge.2.0*bnd_params(4))) exit
        if ((k-1)*pc_binsz.ge.2.0*bnd_params(4)) supvol = .true.
      else if (bnd_shape.eq.3) then
        bndr0t = sqrt(bnd_params(5)+bnd_params(6)**2)
        if ((k.gt.maxk).AND.((k-1)*pc_binsz.ge.bndr0t)) exit
        if ((k-1)*pc_binsz.ge.bndr0t) supvol = .true.
      end if
    end if
    do imt=1,nangrps
      do jmt=imt,nangrps
        if (prtmat(imt,jmt).EQV..true.) then
          howmany = howmany + 1
          vals(howmany) = rbc_pc(imt,jmt,k)
        end if
      end do
    end do
    if (supvol.EQV..true.) then
      write(iu,66) pcb,0.0,(vals(i),i=1,howmany)
    else
      write(iu,66) pcb,voli_pc(k),(vals(i),i=1,howmany)
    end if
  end do
!
  close(unit=iu)
!
  deallocate(normer)
  deallocate(vals)
  deallocate(tys)
  deallocate(prtmat)
!
end
!
!---------------------------------------------------------------------
!
! the generic amide PC calculator is just moderately useful in most settings
! except clean polyamide systems.
! it is working purely on intramolecular distances and henceforth observes no periodic
! boundary conditions
!
subroutine do_amid_pc()
!
  use atoms
  use polypep
  use sequen
  use molecule
  use paircorr
  use grandensembles
  use system
  use pdb
!
  implicit none
!
  integer i,j,puta,putb,putc,imol
  RTYPE d2a,d2b,d2c,v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3),v9(3)
  RTYPE v1(3),v10(3),incr
  logical ismember
!
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  if ((ampc%do_bs.EQV..false.).AND.(ampc%do_ss.EQV..false.)&
 &               .AND.(ampc%do_bb.EQV..false.)) return
!
  do imol=1,nmol
!
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
    do i=rsmol(imol,1),rsmol(imol,2)
!
!     the intra-residue is by definition bb-sc
      if ((ampc%do_bs.EQV..true.).AND.(i.gt.rsmol(imol,1)).AND.&
 &        ((seqtyp(i).eq.17).OR.(seqtyp(i).eq.19))) then 
!     sidechain to "own" backbone (peptide bond connecting residues i-1,i)
        if ((ci(i-1).gt.0).AND.(oi(i-1).gt.0)) then
          call get_CONHH_sc(i,v1,v2,v3,v4,v5)
!         middle of C-N bond vs. middle of C-N bond
          d2a = (0.5*(v1(1) + v3(1) - x(ci(i-1))-x(ni(i))))**2&
 &            + (0.5*(v1(2) + v3(2) - y(ci(i-1))-y(ni(i))))**2&
 &            + (0.5*(v1(3) + v3(3) - z(ci(i-1))-z(ni(i))))**2
          puta = max(1,floor(sqrt(d2a)/pc_binsz) + 1)
          if (puta.gt.MAXPCBINS) puta = MAXPCBINS
          ampc%bs(1,puta) = ampc%bs(1,puta) + incr
!         bb oxygen vs. sc hydrogen on opposite side of sc oxygen ("trans"-hydrogen)
          d2b = (x(oi(i-1)) - v4(1))**2&
 &            + (y(oi(i-1)) - v4(2))**2&
 &            + (z(oi(i-1)) - v4(3))**2
          putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
          if (putb.gt.MAXPCBINS) putb = MAXPCBINS
          ampc%bs(2,putb) = ampc%bs(2,putb) + incr
!         bb oxygen vs. sc hydrogen on same side as sc oxygen ("cis"-hydrogen)
          d2b = (v5(1) - x(oi(i-1)))**2&
 &            + (v5(2) - y(oi(i-1)))**2&
 &            + (v5(3) - z(oi(i-1)))**2
          putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
          if (putb.gt.MAXPCBINS) putb = MAXPCBINS
          ampc%bs(3,putb) = ampc%bs(3,putb) + incr
!         sc oxygen vs. bb hydrogen
          d2b = (v2(1) - x(hni(i)))**2&
 &            + (v2(2) - y(hni(i)))**2&
 &            + (v2(3) - z(hni(i)))**2
          putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
          if (putb.gt.MAXPCBINS) putb = MAXPCBINS
          ampc%bs(4,putb) = ampc%bs(4,putb) + incr
        end if
      end if
!
      do j=i+1,rsmol(imol,2)
!
!       backbone-backbone
!
        if (ampc%do_bb.EQV..true.) then
!
          if ((i.gt.rsmol(imol,1)).AND.((seqpolty(i).eq.'P').OR.(seqpolty(j).eq.'P'))) then
!
!           middle of C-N bond vs. middle of C-N bond
            if ((ci(i-1).gt.0).AND.(ci(j-1).gt.0).AND.(ni(i).gt.0)&
 &              .AND.(ni(j).gt.0)) then
              d2a=(0.5*(x(ci(i-1))+x(ni(i))-x(ci(j-1))-x(ni(j))))**2&
 &              + (0.5*(y(ci(i-1))+y(ni(i))-y(ci(j-1))-y(ni(j))))**2&
 &              + (0.5*(z(ci(i-1))+z(ni(i))-z(ci(j-1))-z(ni(j))))**2
              puta = max(1,floor(sqrt(d2a)/pc_binsz) + 1)
              if (puta.gt.MAXPCBINS) puta = MAXPCBINS
              ampc%bb(1,puta) = ampc%bb(1,puta) + incr
            end if
 
!         oxygen vs. hydrogen
            if ((hni(i).gt.0).AND.(oi(j-1).gt.0)) then
              d2b = (x(hni(i)) - x(oi(j-1)))**2&
 &                + (y(hni(i)) - y(oi(j-1)))**2&
 &                + (z(hni(i)) - z(oi(j-1)))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bb(2,putb) = ampc%bb(2,putb) + incr
            end if
            if ((hni(j).gt.0).AND.(oi(i-1).gt.0)) then   
              d2c = (x(hni(j)) - x(oi(i-1)))**2&
 &                + (y(hni(j)) - y(oi(i-1)))**2&
 &                + (z(hni(j)) - z(oi(i-1)))**2
              putc = max(1,floor(sqrt(d2c)/pc_binsz) + 1)
              if (putc.gt.MAXPCBINS) putc = MAXPCBINS
              ampc%bb(2,putc) = ampc%bb(2,putc) + incr
            end if
!
          end if
        end if
!
!       sidechain-sidechain
!
        if (ampc%do_ss.EQV..true.) then
!
          if (((seqtyp(i).eq.17).OR.(seqtyp(i).eq.19)).AND.&
 &          ((seqtyp(j).eq.17).OR.(seqtyp(j).eq.19))) then
            call get_CONHHCONHH_scsc&
 &                  (i,j,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
!           middle of C-N bond vs. middle of C-N bond
            d2a = (0.5*(v1(1) + v3(1) - v6(1) - v8(1)))**2&
 &              + (0.5*(v1(2) + v3(2) - v6(2) - v8(2)))**2&
 &              + (0.5*(v1(3) + v3(3) - v6(3) - v8(3)))**2
            puta = max(1,floor(sqrt(d2a)/pc_binsz) + 1)
            if (puta.gt.MAXPCBINS) puta = MAXPCBINS
            ampc%ss(1,puta) = ampc%ss(1,puta) + incr
!           oxygen vs. hydrogen on opposite side of oxygen ("trans"-hydrogen)
            d2b = (v2(1) - v9(1))**2&
 &              + (v2(2) - v9(2))**2&
 &              + (v2(3) - v9(3))**2
            d2c = (v4(1) - v7(1))**2&
 &              + (v4(2) - v7(2))**2&
 &              + (v4(3) - v7(3))**2
            putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
            if (putb.gt.MAXPCBINS) putb = MAXPCBINS
            putc = max(1,floor(sqrt(d2c)/pc_binsz) + 1)
            if (putc.gt.MAXPCBINS) putc = MAXPCBINS
            ampc%ss(2,putb) = ampc%ss(2,putb) + incr
            ampc%ss(2,putc) = ampc%ss(2,putc) + incr
!           oxygen vs. hydrogen on same side as oxygen ("cis"-hydrogen)
            d2b = (v2(1) - v10(1))**2&
 &              + (v2(2) - v10(2))**2&
 &              + (v2(3) - v10(3))**2
            d2c = (v5(1) - v7(1))**2&
 &              + (v5(2) - v7(2))**2&
 &              + (v5(3) - v7(3))**2
            putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
            if (putb.gt.MAXPCBINS) putb = MAXPCBINS
            putc = max(1,floor(sqrt(d2c)/pc_binsz) + 1)
            if (putc.gt.MAXPCBINS) putc = MAXPCBINS
            ampc%ss(3,putb) = ampc%ss(3,putb) + incr
            ampc%ss(3,putc) = ampc%ss(3,putc) + incr
          end if
        end if
!
!     i-sidechain-j-backbone
!
        if (ampc%do_bs.EQV..true.) then
!
          if ((seqtyp(i).eq.17).OR.(seqtyp(i).eq.19)) then
            call get_CONHH_sc(i,v1,v2,v3,v4,v5)
!           middle of C-N bond vs. middle of C-N bond
            if ((ci(j-1).gt.0).AND.(ni(j).gt.0)) then
              d2a = (0.5*(v1(1) + v3(1) - x(ci(j-1))-x(ni(j))))**2&
 &                + (0.5*(v1(2) + v3(2) - y(ci(j-1))-y(ni(j))))**2&
 &                + (0.5*(v1(3) + v3(3) - z(ci(j-1))-z(ni(j))))**2
              puta = max(1,floor(sqrt(d2a)/pc_binsz) + 1)
              if (puta.gt.MAXPCBINS) puta = MAXPCBINS
              ampc%bs(1,puta) = ampc%bs(1,puta) + incr
            end if
!           bb oxygen vs. sc hydrogen on opposite side of sc oxygen ("trans"-hydrogen)
            if (oi(j-1).gt.0) then
              d2b = (x(oi(j-1)) - v4(1))**2&
 &                + (y(oi(j-1)) - v4(2))**2&
 &                + (z(oi(j-1)) - v4(3))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bs(2,putb) = ampc%bs(2,putb) + incr
!             bb oxygen vs. sc hydrogen on same side as sc oxygen ("cis"-hydrogen)
              d2b = (v5(1) - x(oi(j-1)))**2&
 &                + (v5(2) - y(oi(j-1)))**2&
 &                + (v5(3) - z(oi(j-1)))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bs(3,putb) = ampc%bs(3,putb) + incr
            end if
!           sc oxygen vs. bb hydrogen
            if (hni(j).gt.0) then
              d2b = (v2(1) - x(hni(j)))**2&
 &                + (v2(2) - y(hni(j)))**2&
 &                + (v2(3) - z(hni(j)))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bs(4,putb) = ampc%bs(4,putb) + incr
            end if
          end if
!
!         i-backbone-j-sidechain
!
          if (((seqtyp(j).eq.17).OR.(seqtyp(j).eq.19)).AND.&
 &            (i.gt.rsmol(imol,1))) then
            call get_CONHH_sc(j,v1,v2,v3,v4,v5)
!           middle of C-N bond vs. middle of C-N bond
            if ((ci(i-1).gt.0).AND.(ni(i).gt.0)) then
              d2a = (0.5*(v1(1) + v3(1) - x(ci(i-1))-x(ni(i))))**2&
 &                + (0.5*(v1(2) + v3(2) - y(ci(i-1))-y(ni(i))))**2&
 &                + (0.5*(v1(3) + v3(3) - z(ci(i-1))-z(ni(i))))**2
              puta = max(1,floor(sqrt(d2a)/pc_binsz) + 1)
              if (puta.gt.MAXPCBINS) puta = MAXPCBINS
              ampc%bs(1,puta) = ampc%bs(1,puta) + incr
            end if
!           bb oxygen vs. sc hydrogen on opposite side of sc oxygen ("trans"-hydrogen)
            if (oi(i-1).gt.0) then
              d2b = (x(oi(i-1)) - v4(1))**2&
 &                + (y(oi(i-1)) - v4(2))**2&
 &                + (z(oi(i-1)) - v4(3))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bs(2,putb) = ampc%bs(2,putb) + incr
!             bb oxygen vs. sc hydrogen on same side as sc oxygen ("cis"-hydrogen)
              d2b = (v5(1) - x(oi(i-1)))**2&
 &                + (v5(2) - y(oi(i-1)))**2&
 &                + (v5(3) - z(oi(i-1)))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bs(3,putb) = ampc%bs(3,putb) + incr
            end if
!           sc oxygen vs. bb hydrogen
            if (hni(i).gt.0) then
              d2b = (v2(1) - x(hni(i)))**2&
 &                + (v2(2) - y(hni(i)))**2&
 &                + (v2(3) - z(hni(i)))**2
              putb = max(1,floor(sqrt(d2b)/pc_binsz) + 1)
              if (putb.gt.MAXPCBINS) putb = MAXPCBINS
              ampc%bs(4,putb) = ampc%bs(4,putb) + incr
            end if
!
          end if
        end if
!
      end do
    end do
!
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine prt_amid_pc()
!
  use paircorr
  use mpistuff
!
  implicit none
!
  integer lastbin,i,j,freeunit,iu,ii,jj
  RTYPE t1,t2,t3,t4,t5,t6,t7,t8,t9,togoout(9)
  character(60) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  logical exists
!
! sanity check
  if ((ampc%do_bs.EQV..false.).AND.(ampc%do_ss.EQV..false.)&
 &               .AND.(ampc%do_bb.EQV..false.)) return
!
! normalize the raw counts
  lastbin = 0
  if (ampc%do_bb.EQV..true.) then
    t1 = sum(ampc%bb(1,:))
    t2 = sum(ampc%bb(2,:))
  end if
  if (ampc%do_ss.EQV..true.) then
    t3 = sum(ampc%ss(1,:))
    t4 = sum(ampc%ss(2,:))
    t5 = sum(ampc%ss(3,:))
  end if
  if (ampc%do_bs.EQV..true.) then
    t6 = sum(ampc%bs(1,:))
    t7 = sum(ampc%bs(2,:))
    t8 = sum(ampc%bs(3,:))
    t9 = sum(ampc%bs(4,:))
  end if
  do i=1,MAXPCBINS
    if (ampc%do_bb.EQV..true.) then
      if ((ampc%bb(1,i).gt.0.0).OR.(ampc%bb(2,i).gt.0.0)) lastbin = i
    end if
    if (ampc%do_ss.EQV..true.) then
      if ((ampc%ss(1,i).gt.0.0).OR.(ampc%ss(2,i).gt.0.0).OR.(ampc%ss(3,i).gt.0.0)) lastbin = i
    end if
    if (ampc%do_bs.EQV..true.) then
      if ((ampc%bs(2,i).gt.0.0).OR.(ampc%bs(3,i).gt.0.0).OR.&
 &        (ampc%bs(4,i).gt.0.0).OR.(ampc%bs(1,i).gt.0.0)) lastbin = i
    end if
  end do
!
! another sanity check (should be redundant)
  if (lastbin.eq.0) return
!
  if (ampc%do_bb.EQV..true.) then
    if (t1.gt.0.0) ampc%bb(1,:) = ampc%bb(1,:)/t1
    if (t2.gt.0.0) ampc%bb(2,:) = ampc%bb(2,:)/t2
  end if
  if (ampc%do_ss.EQV..true.) then
    if (t3.gt.0.0) ampc%ss(1,:) = ampc%ss(1,:)/t3
    if (t4.gt.0.0) ampc%ss(2,:) = ampc%ss(2,:)/t4
    if (t5.gt.0.0) ampc%ss(3,:) = ampc%ss(3,:)/t5
  end if
  if (ampc%do_bs.EQV..true.) then
    if (t6.gt.0.0) ampc%bs(1,:) = ampc%bs(1,:)/t6
    if (t7.gt.0.0) ampc%bs(2,:) = ampc%bs(2,:)/t7
    if (t8.gt.0.0) ampc%bs(3,:) = ampc%bs(3,:)/t8
    if (t9.gt.0.0) ampc%bs(4,:) = ampc%bs(4,:)/t9
  end if
!
! now let's go ahead and print, at the moment the conversion to a PC fxn has to happen outside
! because of the most likely numerically obtained prior/reference
! 
 46   format(f14.7,1x,9(g14.7,1x))
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod//'_AMIDES_PC.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'AMIDES_PC.dat'
  end if
#else
  fn = 'AMIDES_PC.dat'
#endif 
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  do i=1,lastbin
    if (ampc%do_bb.EQV..true.) then
      togoout(1:2) = ampc%bb(1:2,i)
    else
      togoout(1:2) = 0.0d0
    end if
    if (ampc%do_ss.EQV..true.) then
      togoout(3:5) = ampc%ss(1:3,i)
    else
      togoout(3:5) = 0.0d0
    end if
    if (ampc%do_bs.EQV..true.) then
      togoout(6:9) = ampc%bs(1:4,i)
    else
      togoout(6:9) = 0.0d0
    end if
    write(iu,46) (i-0.5)*pc_binsz,(togoout(j),j=1,9)
  end do
  close(unit=iu)
!
end
!
!---------------------------------------------------------------------
!
! ROUTINES FOR COVARIANCE ANALYSIS
!
!---------------------------------------------------------------------
!
subroutine setup_covar_int()
!
  use iounit
  use fyoc
  use torsn
  use molecule
  use mpistuff
  use sequen
  use system
  use mcsums
  use aminos
!
  implicit none
!
  integer dc(nangrps),rs,mt
!
 55   format('# Analysis group ',i3,' (Ref.: Mol. ',i5,' / Res. ',i5,'-',&
 &i5,') has ',i5,' torsions / angles.')
!
  do mt=1,nangrps
    dc(mt) = 0
    if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
    do rs=rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
!
      if (notors(rs).EQV..true.) cycle
!
      if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) dc(mt) = dc(mt) + 1
      if (yline(rs).gt.0) dc(mt) = dc(mt) + 1
      if (wline(rs).gt.0) dc(mt) = dc(mt) + 1
      dc(mt) = dc(mt) + nnucs(rs)
      dc(mt) = dc(mt) + nchi(rs)
      if ((pucline(rs).gt.0).OR.((seqpolty(rs).eq.'N').AND.(nucsline(6,rs).gt.0))) then
        dc(mt) = dc(mt) + 7
      end if
    end do
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(myrank.eq.0)) then
      write(itrcv(mt),55) mt,molangr(mt,1),rsmol(molangr(mt,1),1),&
 &rsmol(molangr(mt,1),2),dc(mt)
    else if (use_REMC.EQV..true.) then
      write(itrcv(mt),55) mt,molangr(mt,1),rsmol(molangr(mt,1),1),&
 &rsmol(molangr(mt,1),2),dc(mt)
    end if
#else
    write(itrcv(mt),55) mt,molangr(mt,1),rsmol(molangr(mt,1),1),&
 &rsmol(molangr(mt,1),2),dc(mt)
#endif
  end do
!
  if (covmode.eq.3) then
    if (ntorlcs.le.0) then
      write(ilog,*) 'Requested LCT covariance analysis (mode ',&
 &covmode,'), but no LCs are provided. Turning off.'
      covcalc = nsim+1
!     manually shutdown
#ifdef ENABLE_MPI
      if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
          do mt=1,nangrps
            close(unit=itrcv(mt),status='delete')
          end do
        end if
      else if (use_REMC.EQV..true.) then
        do mt=1,nangrps
          close(unit=itrcv(mt),status='delete')
        end do
      end if
#else
      do mt=1,nangrps
        close(unit=itrcv(mt),status='delete')
      end do
#endif
      return
    end if
  end if
!
end
!
!---------------------------------------------------------------------
!
subroutine do_covar_int()
!
  use iounit
  use sequen
  use fyoc
  use math
  use molecule
  use torsn
  use movesets
  use mpistuff
  use energies
  use grandensembles
  use system
  use aminos
  use zmatrix
  use polypep
!
  implicit none
!
  integer i,j,dc,rs,mt,imol,shf,shfbu
  RTYPE dum
  logical ismember
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
 64   format(2000000(g11.4,1x))
!
! increase counter
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    do imol=1,nmol
      if (do_pol(moltypid(imol)).EQV..false.) cycle
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
      mt = an_grp_mol(imol)
      ncovaravg(mt) = ncovaravg(mt) + 1
    end do
  else
    do mt=1,nangrps
      if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
      ncovaravg(mt) = ncovaravg(mt) + molangr(mt,2)
    end do
  end if
!
  do imol=1,nmol
    dc = 0
    mt = moltypid(imol)
    if (do_pol(mt).EQV..false.) cycle
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!    
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      if (notors(rs).EQV..true.) cycle
!
      if (wline(rs).gt.0) then
        dc = dc + 1
        curtvec(dc) = omega(rs)
      end if
      if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
        dc = dc + 1
        curtvec(dc) = phi(rs)
      end if
      if (yline(rs).gt.0) then
        dc = dc + 1
        curtvec(dc) = psi(rs)
      end if
      do i=1,nnucs(rs)
        dc = dc + 1
        curtvec(dc) = nucs(i,rs)
      end do
      do i=1,nchi(rs)
        dc = dc + 1
        curtvec(dc) = chi(i,rs)
      end do
      if ((pucline(rs).gt.0).OR.((seqpolty(rs).eq.'N').AND.(nucsline(6,rs).gt.0))) then
        if (seqpolty(rs).eq.'N') then
          shfbu = shf
          shf = 1
        end if
        if (seqpolty(rs).eq.'N') then
          dc = dc + 1
          curtvec(dc) = ztor(nucsline(6,rs))
        else
          dc = dc + 1
          curtvec(dc) = ztor(fline(rs))
        end if
        dc = dc + 1
        curtvec(dc) = ztor(at(rs)%sc(2-shf))
        dc = dc + 1
        curtvec(dc) = ztor(at(rs)%sc(3-shf))
        dc = dc + 1
        curtvec(dc) = ztor(at(rs)%sc(4-shf))
        dc = dc + 1
        curtvec(dc) = bang(at(rs)%sc(2-shf))
        dc = dc + 1
        curtvec(dc) = bang(at(rs)%sc(3-shf))
        dc = dc + 1
        curtvec(dc) = bang(at(rs)%sc(4-shf))
        if (seqpolty(rs).eq.'N') shf = shfbu
      end if
    end do
    if (covmode.eq.1) then
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        write(itrcv(an_grp_mol(imol)),64) (curtvec(i)/RADIAN,i=1,dc)
      else if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
          call MPI_AVGWriteTRCV(an_grp_mol(imol),dc,curtvec)
        else 
          call MPI_AVGSendTRCV(an_grp_mol(imol),dc,curtvec)
        end if
      end if
#else
      write(itrcv(an_grp_mol(imol)),64) (curtvec(i)/RADIAN,i=1,dc)
#endif
    else if (covmode.eq.2) then
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        write(itrcv(an_grp_mol(imol)),64) (cos(curtvec(i)/RADIAN),sin(curtvec(i)/RADIAN)&
 &                                                      ,i=1,dc)
      else if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
          call MPI_AVGWriteTRCV(an_grp_mol(imol),dc,curtvec)
        else 
          call MPI_AVGSendTRCV(an_grp_mol(imol),dc,curtvec)
        end if
      end if
#else
      write(itrcv(an_grp_mol(imol)),64) (cos(curtvec(i)/RADIAN),sin(curtvec(i)/RADIAN)&
 &                                                      ,i=1,dc)
#endif
    else if (covmode.eq.3) then
      do i=1,ntorlcs
        dum = 0.0
        if (use_lctmoves.EQV..false.) then
          do j=1,dc
            if (torlcmode.eq.1) then
              dum = dum + torlc_coeff(i,2*j-1)*cos(curtvec(j)/RADIAN)&
 &                      + torlc_coeff(i,2*j)*sin(curtvec(j)/RADIAN)
            else if (torlcmode.eq.2) then
              dum = dum + torlc_coeff(i,j)*(curtvec(j)/RADIAN)
            end if
          end do
          curtvec(i) = dum
        else
          curtvec(i) = lct(i)
        end if
      end do
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        write(itrcv(an_grp_mol(imol)),64) (curtvec(i),i=1,ntorlcs)
      else if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
          call MPI_AVGWriteTRCV(an_grp_mol(imol),ntorlcs,curtvec)
        else 
          call MPI_AVGSendTRCV(an_grp_mol(imol),ntorlcs,curtvec)
        end if
      end if
#else
      write(itrcv(an_grp_mol(imol)),64) (curtvec(i),i=1,ntorlcs)
#endif
    end if
!
  end do
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (mpi_cnt_trcv.eq.(mpi_nodes-1)) then
      mpi_cnt_trcv = 0
    else
      mpi_cnt_trcv = mpi_cnt_trcv + 1
    end if
  end if
#endif
!  
end
!
!------------------------------------------------------------------------
!
! this routine is a simple wrapper to align the current coordinates to
! a reference set and uses align_3D(...)
!
subroutine struct_align()
!
  use atoms
  use clusters
  use pdb
  use iounit
  use molecule
  use sequen
  use polyavg
  use mcsums
  use system, ONLY: bnd_type
!
  implicit none
!
  integer i,imol,lmol
  RTYPE cto(3),ctn(3),tvec(3),qrot(4),rmsdv,rmsdv2,shcom(3)
!
  if (align%haveset.EQV..false.) then
    if (use_pdb_template.EQV..true.) then
      write(ilog,*) 'Fatal. Alignment set has no coordinates assigned even though &
 & a template is being used. This is a bug.'
      call fexit()
    end if
!   the case for sequential alignment needs to be seeded
    lmol = 0
    tvec(:) = 0.0
    do i=1,align%nr
      imol = molofrs(atmres(align%set(i))) ! list is sorted
      if (align%mmol.gt.0) then
        if (lmol.ne.imol) then
          call shift_bound3(align%mmol,imol,tvec)
          shcom(:) = tvec(:) ! conversion
        else
          shcom(:) = tvec(:)
        end if
        lmol = imol
      else
        shcom(:) = 0.0
      end if
      align%refxyz(3*i-2) = x(align%set(i)) + shcom(1)
      align%refxyz(3*i-1) = y(align%set(i)) + shcom(2)
      align%refxyz(3*i) = z(align%set(i)) + shcom(3)
    end do
    if (align%diffnr.gt.0) then
!     if needed, also seed difference coordinates into hlpxyz
      lmol = 0
      tvec(:) = 0.0
      do i=1,align%diffnr
        imol = molofrs(atmres(align%diffset(i))) ! list is sorted
        if (align%mmol.gt.0) then
          if (lmol.ne.imol) then
            call shift_bound3(align%mmol,imol,tvec)
            shcom(:) = tvec(:) ! conversion
          else
            shcom(:) = tvec(:)
          end if
          lmol = imol
        else
          shcom(:) = 0.0
        end if
        align%hlpxyz(3*i-2) = x(align%diffset(i)) + shcom(1)
        align%hlpxyz(3*i-1) = y(align%diffset(i)) + shcom(2)
        align%hlpxyz(3*i) = z(align%diffset(i)) + shcom(3)
      end do
    end if
    align%haveset = .true.
    return
  end if
!
  if ((align%mmol.gt.0).AND.(nmol.gt.1).AND.(bnd_type.eq.1)) then
!   bring all coordinates into reference frame prior to alignment if so desired
    lmol = 0
    tvec(:) = 0.0
    do i=1,n
      imol = molofrs(atmres(i)) ! list is sorted
      if (lmol.ne.imol) then
        call shift_bound3(align%mmol,imol,tvec)
        shcom(:) = tvec(:) ! conversion
      else
        shcom(:) = tvec(:)
      end if
      lmol = imol
      x(i) = x(i) + shcom(1)
      y(i) = y(i) + shcom(2)
      z(i) = z(i) + shcom(3)
    end do
  end if
! collect new PBC-corrected coordinates for alignment set
  do i=1,align%nr
    align%curxyz(3*i-2) = x(align%set(i))
    align%curxyz(3*i-1) = y(align%set(i))
    align%curxyz(3*i) = z(align%set(i))
  end do
!
  call align_3D(align%nr,align%curxyz,align%refxyz,tvec,qrot,cto,ctn)
!
! rotate in place (all!)
  do i=1,n
    call quat_conjugate4(qrot,i,cto)
  end do
! translate
  x(1:n) = x(1:n) + tvec(1)
  y(1:n) = y(1:n) + tvec(2)
  z(1:n) = z(1:n) + tvec(3)
  do i=1,nmol
    call update_rigid(i)
  end do
!
! compute and print RMSD - note that this destroys the coordinates prior to alignment that are still in align%curxyz
  if (align%instrmsd.EQV..true.) then
    if (align%diffnr.gt.0) then
      rmsdv2 = 0.0
      do i=1,align%diffnr
        rmsdv2 = rmsdv2 + (align%hlpxyz(3*i-2) - (x(align%diffset(i))))**2 + &
 &                        (align%hlpxyz(3*i-1) - (y(align%diffset(i))))**2 + &
 &                        (  align%hlpxyz(3*i) - (z(align%diffset(i))))**2
      end do 
      rmsdv2 = sqrt(rmsdv2/(1.0*align%diffnr))
    end if
    do i=1,align%nr
      align%curxyz(3*i-2) = x(align%set(i))
      align%curxyz(3*i-1) = y(align%set(i))
      align%curxyz(3*i) = z(align%set(i))
    end do
    rmsdv = sqrt(sum((align%curxyz(1:3*align%nr)-align%refxyz(1:3*align%nr))**2)/(1.0*align%nr))
 546 format(f13.5,1x,f13.5)
    if (align%diffnr.eq.0) then
      write(irmsd,546) rmsdv
    else
      write(irmsd,546) rmsdv,rmsdv2
    end if
  end if
!
  if (align%refset.EQV..false.) then
    do i=1,align%nr
      align%refxyz(3*i-2) = x(align%set(i))
      align%refxyz(3*i-1) = y(align%set(i))
      align%refxyz(3*i) = z(align%set(i))
    end do
    if (align%diffnr.gt.0) then
      do i=1,align%diffnr
        align%hlpxyz(3*i-2) = x(align%diffset(i))
        align%hlpxyz(3*i-1) = y(align%diffset(i))
        align%hlpxyz(3*i) = z(align%diffset(i))
      end do
    end if
  end if
!
end
!
!
!---------------------------------------------------------------------------------------
!
subroutine init_structalign()
!
  use molecule
  use clusters, ONLY: align,cfilen
  use atoms, ONLY: x,y,z,atmres,n
  use sequen, ONLY: molofrs
  use iounit
!
  implicit none
!
  integer imol,lmol,i,iu,t1,t2,freeunit,nlst
  integer, ALLOCATABLE:: tmplst(:)
  RTYPE tvec(3),shcom(3)
  logical exists,dum(2)
!
  do imol=1,nmol
    com(imol,1) = sum(x(atmol(imol,1):atmol(imol,2)))
    com(imol,2) = sum(y(atmol(imol,1):atmol(imol,2)))
    com(imol,3) = sum(z(atmol(imol,1):atmol(imol,2)))
    com(imol,:) = com(imol,:)/dble(atmol(imol,2)-atmol(imol,1)+1)
  end do
  if (align%yes.EQV..true.) then
    lmol = 0
    tvec(:) = 0.0
    do i=1,align%nr
      imol = molofrs(atmres(align%set(i))) ! list is sorted
      if (align%mmol.gt.0) then
        if (lmol.ne.imol) then
          call shift_bound3(align%mmol,imol,tvec)
          shcom(:) = tvec(:) ! conversion
        else
          shcom(:) = tvec(:)
        end if
        lmol = imol
      else
        shcom(:) = 0.0
      end if
      align%refxyz(3*i-2) = x(align%set(i)) + shcom(1)
      align%refxyz(3*i-1) = y(align%set(i)) + shcom(2)
      align%refxyz(3*i) = z(align%set(i)) + shcom(3)
    end do
    align%haveset = .true.
    align%refset = .true.
!   do we have a separate distance set
    if (align%instrmsd.EQV..true.) then
      call strlims(cfilen,t1,t2)
      inquire(file=cfilen(t1:t2),exist=exists)
      if (exists.EQV..false.) then
        write(ilog,*) 'Using the same distance and alignments sets for RMSD calculation.'
      else
        write(ilog,*) 'CAMPARI will try to use the index information in ',cfilen(t1:t2),' as distance set.'
        iu = freeunit()
        allocate(tmplst(n))
        open(unit=iu,file=cfilen(t1:t2),status='old')
        lmol = n
        nlst = 0
        dum(:) = .true.
        call read_atmidx(iu,nlst,tmplst,lmol,dum(1),dum(2))
        close(unit=iu)
!       sorted and nonredundant list -> just copy
        if (nlst.gt.0) then
          align%diffnr = nlst
          allocate(align%diffset(nlst))
          align%diffset(1:nlst) = tmplst(1:nlst)
          allocate(align%hlpxyz(3*nlst))
          deallocate(tmplst)
!         if needed, store original coordinates in hlpxyz
          lmol = 0
          tvec(:) = 0.0
          do i=1,align%diffnr
            imol = molofrs(atmres(align%diffset(i))) ! list is sorted
            if (align%mmol.gt.0) then
              if (lmol.ne.imol) then
                call shift_bound3(align%mmol,imol,tvec)
                shcom(:) = tvec(:) ! conversion
              else
                shcom(:) = tvec(:)
              end if
              lmol = imol
            else
              shcom(:) = 0.0
            end if
            align%hlpxyz(3*i-2) = x(align%diffset(i)) + shcom(1)
            align%hlpxyz(3*i-1) = y(align%diffset(i)) + shcom(2)
            align%hlpxyz(3*i) = z(align%diffset(i)) + shcom(3)
          end do
        else
          write(ilog,*) 'Warning. Index extraction for RMSD difference set from ',cfilen(t1:t2),' failed.'
        end if
      end if
    end if
  end if
!
end
!
!-------------------------------------------------------------------------
!
subroutine init_aligndisset()
!
  use molecule
  use clusters, ONLY: align,cfilen
  use atoms, ONLY: x,y,z,atmres,n
  use sequen, ONLY: molofrs
  use iounit
!
  implicit none
!
  integer imol,lmol,i,iu,t1,t2,freeunit,nlst
  integer, ALLOCATABLE:: tmplst(:)
  RTYPE tvec(3),shcom(3)
  logical exists,dum(2)
!
  do imol=1,nmol
    com(imol,1) = sum(x(atmol(imol,1):atmol(imol,2)))
    com(imol,2) = sum(y(atmol(imol,1):atmol(imol,2)))
    com(imol,3) = sum(z(atmol(imol,1):atmol(imol,2)))
    com(imol,:) = com(imol,:)/dble(atmol(imol,2)-atmol(imol,1)+1)
  end do
!
  call strlims(cfilen,t1,t2)
  inquire(file=cfilen(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Using the same distance and alignments sets for RMSD calculation.'
  else
    write(ilog,*) 'CAMPARI will try to use the index information in ',cfilen(t1:t2),' as distance set.'
    iu = freeunit()
    allocate(tmplst(n))
    open(unit=iu,file=cfilen(t1:t2),status='old')
    lmol = n
    nlst = 0
    dum(:) = .true.
    call read_atmidx(iu,nlst,tmplst,lmol,dum(1),dum(2))
    close(unit=iu)
!   sorted and nonredundant list -> just copy
    if (nlst.gt.0) then
      align%diffnr = nlst
      allocate(align%diffset(nlst))
      align%diffset(1:nlst) = tmplst(1:nlst)
      allocate(align%hlpxyz(3*nlst))
      deallocate(tmplst)
!     if needed, store original coordinates in hlpxyz
      lmol = 0
      tvec(:) = 0.0
      do i=1,align%diffnr
        imol = molofrs(atmres(align%diffset(i))) ! list is sorted
        if (align%mmol.gt.0) then
          if (lmol.ne.imol) then
            call shift_bound3(align%mmol,imol,tvec)
            shcom(:) = tvec(:) ! conversion
          else
            shcom(:) = tvec(:)
          end if
          lmol = imol
        else
          shcom(:) = 0.0
        end if
        align%hlpxyz(3*i-2) = x(align%diffset(i)) + shcom(1)
        align%hlpxyz(3*i-1) = y(align%diffset(i)) + shcom(2)
        align%hlpxyz(3*i) = z(align%diffset(i)) + shcom(3)
      end do
    else
      write(ilog,*) 'Warning. Index extraction for RMSD difference set from ',cfilen(t1:t2),' failed.'
    end if
  end if

end
!
!---------------------------------------------------------------------------
!
! #################################################
! ##                                             ##
! ## subroutine hbond -- backbone hydrogen bonds ##
! ##                                             ##
! #################################################
!
!
! "hbond" scans a single chain pdb file for all backbone
! hydrogen bonds based on  distance criteria
!
!
!  subroutine hbond(wherehb)
!c
!  use polypep
!  use saw
!  use sequen
!c
!  implicit none
!c
!  integer i,j,k,diff,diff2
!  integer bid,cid,dd,ddp,hid,nid,oid,xx
!  integer wherehb(100)
!  RTYPE ang1,ang2,getbang,getblen,getztor,ptorsn
!  RTYPE ondist,tau1,tau2
!  logical accpt(nseq),donor(nseq)
!  logical plane1,plane2
!c
!  xmax = nseq -1
!  do i = 1,nseq
!    accpt(i) = .false.
!    donor(i) = .false.
!  end do
!c     
!  nhbond = 0
!  plane1 = .false.
!  plane2 = .false.
!  do xx = 3,xmax
!    do i = 2,nseq-xx
!      if(accpt(i)) go to 20
!      do j = i+xx, nseq-1
!        if((j-i) .ne. xx) go to 10
!        if(donor(j)) go to 10
!        oid = oi(i)
!        cid = ci(i)
!        bid = ni(i+1)
!        dd  = cai(j)
!        hid = hni(j)
!        nid = ni(j)
!        ddp = ci(j-1)
!c           
!        ondist = getblen(nid,oid)
!        ang1 = getbang(nid,oid,cid)
!        ang2 = getbang(oid,nid,dd)                    
!        tau1 = getztor(nid,oid,cid,bid)
!        plane1 = .false.
!        if(tau1 .lt. 0.0d0 .and. tau1 .gt. -90.0d0) then 
!          plane1 = .true. 
!        else if(tau1 .gt. 0.0d0 .and. tau1 .lt. 90.0d0) then 
!          plane1 = .true.
!        end if 
!        plane2 = .false.
!        tau2 = ptorsn(ddp,dd,nid,oid)
!        plane2 = .false.
!        if(tau2 .ge. -60.0d0 .and. tau2 .lt. 0.0d0) then
!          plane2 = .true. 
!        else if(tau2 .ge. 0.0d0 .and. tau2 .le. 60.0d0) then 
!          plane2 = .true.
!        end if
!c           
!        if(ondist .le. 5.0d0) then 
!
!c            if(ondist .le. 3.2d0 .and. ang1 .gt. 90.0d0 
!c     &        .and. ang2 .gt. 90.0d0 .and. plane1 
!c     &        .and. plane2) then 
!          nhbond = nhbond + 1
!          accpt(i) = .true.
!          donor(j) = .true.
!          wherehb((nhbond*2)-1)=i
!          wherehb(nhbond*2)=j
!        end if
! 10       continue
!      end do
! 20     continue
!    end do
!  end do
!  return
!  end
