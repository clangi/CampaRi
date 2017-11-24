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
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ##########################################
! ##                                      ##
! ##  SECONDARY STRUCTURAL KABSCH-SANDER  ##
! ##         ASSIGNMENT ROUTINES          ##
! ##                                      ##
! ##########################################
!
!
! since string operations are stupid we use an integer code internally:
! first column (global assignment):
!    0: unassigned ( )
!    1: single-paired parallel strand (E)
!    2: single-paired antiparallel strand (E)
!    3: double-paired all-parallel strand in sheet (E)
!    4: double-paired mixed parallel-antiparallel strand in sheet (E)
!    5: double-paired all-antiparallel strand in sheet (E)
!    6: canonical alpha-helix (H)
!    7: 3/10-helix (G)
!    8: PI-helix (I)
!    9: Generic turn (3,4,5-registry and NOT anything else above) (T)
!   10: Generic bend (bend and NOT anything else above) (S)
!
!  second column (i to i+3 H-bond registry)
!    0: unassigned ( )
!    1: neither donor nor acceptor but surrounded by i-i+3 H-bond fork (3)
!    2: just donor of H-bond to i-3 (<)
!    3: just acceptor for H-bond from i+3 (>)
!    4: donor of H-bond to i-3 and acceptor for H-bond from i+3 (X)
!
!  third column (i to i+4 H-bond registry)
!    0: unassigned ( )
!    1: neither donor nor acceptor but surrounded by i-i+4 H-bond fork (4)
!    2: just donor of H-bond to i-4 (<)
!    3: just acceptor for H-bond from i+4 (>)
!    4: donor of H-bond to i-4 and acceptor for H-bond from i+4 (X)
!
!  fourth column (i to i+5 H-bond registry)
!    0: unassigned ( )
!    1: neither donor nor acceptor but surrounded by i-i+5 H-bond fork (5)
!    2: just donor of H-bond to i-5 (<)
!    3: just acceptor for H-bond from i+5 (>)
!    4: donor of H-bond to i-5 and acceptor for H-bond from i+5 (X)
!
!  fifth column (chain bending)
!    0: unassigned ( )
!   10: chain bends by more than 70 degrees over residues i-2,i,i+2 (S)
!
!  sixth column (chain chirality)
!    0: unassigned ( )
!    1: negative chirality along residues i-1,i,i+1,i+2 as is typical for twisted sheets (-)
!    2: positive chirality along residues i-1,i,i+1,i+2 as is typical for right-handed helices (+)
!
!  the seventh and eigth columns are currently not in use
!
!-----------------------------------------------------------------------
!
subroutine setup_dssp()
!
  use iounit
  use polypep
  use aminos
  use sequen
  use dssps
  use system
  use mpistuff
  use energies
  use molecule
!
  implicit none
!
  integer i,aone,k,imol,j
!
  aone = 1
  pep_sz = 0
  imol = 0
  npepmol = 0
!
  do i=1,nseq
!    if (seqpolty(i).eq.'P') then
    if ((seqpolty(i).eq.'P').OR.((seqtyp(i).ge.27).AND.(seqtyp(i).le.30)).OR.&
 &      ((seqtyp(i).eq.26).AND.(i.eq.rsmol(molofrs(i),1)).AND.(oi(i).gt.0).AND.(ci(i).gt.0)).OR.&
 &      ((seqtyp(i).eq.26).AND.(i.eq.rsmol(molofrs(i),2)).AND.(ni(i).gt.0).AND.(hni(i).gt.0)) ) then
      pep_sz = pep_sz + 1
      if (molofrs(i).ne.imol) then
        npepmol = npepmol + 1
        imol = molofrs(i)
      end if
    end if
  end do
!
  if (pep_sz.le.0) then
    write(ilog,*) 'Warning. No polypeptide residues found. Turning D&
 &SSP analysis and/or restraint potential off.'
    dsspcalc = nsim + 1
    use_DSSP = .false.
    scale_DSSP = 0.0
    return
  end if
!
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (inst_dssp.EQV..true.) then
      write(ilog,*) 'Warning. Instantaneous output for DSSP-analysis&
 & is currently not available in MPI averaging runs. Disabled.'
      inst_dssp = .false.
    end if
  end if
#endif
!
  call allocate_dssps(aone)
!
  k = 0
  imol = 0
  j = 0
  do i=1,nseq
!    if (seqpolty(i).eq.'P') then
    if ((seqpolty(i).eq.'P').OR.((seqtyp(i).ge.27).AND.(seqtyp(i).le.30)).OR.&
 &      ((seqtyp(i).eq.26).AND.(i.eq.rsmol(molofrs(i),1)).AND.(oi(i).gt.0).AND.(ci(i).gt.0)).OR.&
 &      ((seqtyp(i).eq.26).AND.(i.eq.rsmol(molofrs(i),2)).AND.(ni(i).gt.0).AND.(hni(i).gt.0)) ) then
      k = k + 1
      peprs(i) = k
      rspep(k) = i
      if (molofrs(i).ne.imol) then
        j = j + 1
        pepmolsz(j) = 1
        pepmol(molofrs(i)) = j
        molpep(j,1) = k
        molpep(j,2) = k
        imol = molofrs(i)
      else
        pepmolsz(j) = pepmolsz(j) + 1
        molpep(j,2) = k
      end if
    end if
  end do
!
! setup natural bounds for HB-score histograms
  if ((par_DSSP(2).le.par_DSSP(1)).OR.&
 &    (par_DSSP(1).lt.par_DSSP(3)).OR.&
 &    (par_DSSP(2).le.par_DSSP(3))) then
    write(ilog,*) 'Warning. Hydrogen-bond energy cutoffs are ill-def&
 &ined. Adjusting to -0.5, -2.5, and -4.0 kcal/mol.'
    par_DSSP(2) = -0.5
    par_DSSP(1) = -2.5
    par_DSSP(3) = -4.0
  end if
  if (par_DSSP2(2).eq.3) then
    par_DSSP(5) = par_DSSP(3)/par_DSSP(1)
    par_DSSP(6) = par_DSSP(5)/100.0
  else
    par_DSSP(5) = 1.0
    par_DSSP(6) = 0.01
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this routine assembles hydrogen bond lists (in hbl), but also
! a network matrix (hb_mat). for the latter the codes are (always rs1 < rs2):
! 1: rs1 is acceptor only, rs2 donor only
! 2: rs2 is donor only, rs1 acceptor only
! 3: rs1 and rs2 are both both acceptor and donor
!
subroutine get_bbhb(ps1,ps2)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
!
  implicit none
!
  integer rs1,rs2,ps1,ps2,rsbu
  RTYPE don1,doh1,dcn1,dch1,don2,doh2,dcn2,dch2,qfac,d2
  RTYPE hben1,hben2,svec(3),dvec1(3),dvec2(3),dvec3(3),dvec4(3)
  RTYPE dvec5(3),dvec6(3),dvec7(3),dvec8(3),dis
!
  rs1 = rspep(ps1)
  rs2 = rspep(ps2)
!  if (((seqpolty(rs1).ne.'P').AND.((seqtyp(rs1).gt.30).OR.(seqtyp(rs1).lt.27))).OR.&
! &    ((seqpolty(rs2).ne.'P').AND.((seqtyp(rs2).gt.30).OR.(seqtyp(rs2).lt.27))).OR.&
  if (abs(rs1-rs2).le.1) then
!   HB-assignments only work for pairs of non-neighbor polypeptide residues
    return
  end if
!
  qfac = -27.888 ! A*kcal/mol
!
! to be on the safe side ...
  if (rs1.gt.rs2) then
    rsbu = rs1
    rs1 = rs2
    rs2 = rsbu
  end if 
!
  svec(:) = 0.0
  call dis_bound_rs3(rs1,rs2,svec,dis)
! use distance screening to save CPU time
! this of course makes assumptions about the energy cutoff values for the
! H-bond formula below ...
  if (dis.gt.par_DSSP(4)) return
!
! we have to screen for PBC as well as for proline which cannot be a donor
!  if (hni(rs2).gt.0) then
   if ((hni(rs2).gt.0).AND.(oi(rs1).gt.0)) then
    call dis_bound5(oi(rs1),ni(rs2),svec,d2,dvec1) 
    don1 = sqrt(1./d2)
    call dis_bound5(oi(rs1),hni(rs2),svec,d2,dvec2)
    doh1 = sqrt(1./d2)
    call dis_bound5(ci(rs1),hni(rs2),svec,d2,dvec3)
    dch1 = sqrt(1./d2)
    call dis_bound5(ci(rs1),ni(rs2),svec,d2,dvec4)
    dcn1 = sqrt(1./d2)
  end if
!  if (hni(rs1).gt.0) then
  if ((hni(rs1).gt.0).AND.(oi(rs2).gt.0)) then
    call dis_bound5(ni(rs1),oi(rs2),svec,d2,dvec5)
    don2 = sqrt(1./d2)
    call dis_bound5(hni(rs1),oi(rs2),svec,d2,dvec6)
    doh2 = sqrt(1./d2)
    call dis_bound5(hni(rs1),ci(rs2),svec,d2,dvec7)
    dch2 = sqrt(1./d2)
    call dis_bound5(ni(rs1),ci(rs2),svec,d2,dvec8)
    dcn2 = sqrt(1./d2)
  end if
!
  if ((hni(rs2).gt.0).AND.(oi(rs1).gt.0)) then
    hben1 = qfac*(doh1 + dcn1 - don1 - dch1)
  else
    hben1 = 0.0 ! remember par_DSSP(2) is always negative
  end if
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = rs2
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = rs1
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   update res-res-flag
    if (hb_mat(peprs(rs1),peprs(rs2)).gt.0) then
      hb_mat(peprs(rs1),peprs(rs2)) = 3
    else
      hb_mat(peprs(rs1),peprs(rs2)) = 1
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
  if ((hni(rs1).gt.0).AND.(oi(rs2).gt.0)) then
    hben2 = qfac*(doh2 + dcn2 - don2 - dch2)
  else
    hben2 = 0.0
  end if
  if (hben2.lt.par_DSSP(2)) then
!   add donor hydrogen bond for rs1
    if (hbl(peprs(rs1))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%ndhbs = hbl(peprs(rs1))%ndhbs + 1
      hbl(peprs(rs1))%dhb(hbl(peprs(rs1))%ndhbs) = rs2
      hbl(peprs(rs1))%endhb(hbl(peprs(rs1))%ndhbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
!   add acceptor hydrogen bond for rs2
    if (hbl(peprs(rs2))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%nahbs = hbl(peprs(rs2))%nahbs + 1
      hbl(peprs(rs2))%ahb(hbl(peprs(rs2))%nahbs) = rs1
      hbl(peprs(rs2))%enahb(hbl(peprs(rs2))%nahbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
!   update res-res-flag
    if (hb_mat(peprs(rs1),peprs(rs2)).gt.0) then
      hb_mat(peprs(rs1),peprs(rs2)) = 3
    else
      hb_mat(peprs(rs1),peprs(rs2)) = 2
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same routine which also computes the derivatives of all H-bond energies
! with respect to the reference atoms
!
subroutine get_bbhb_der(ps1,ps2)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
!
  implicit none
!
  integer rs1,rs2,ps1,ps2,j,rsbu
  RTYPE id1(4),id2(4),d2(4),qfac,dis
  RTYPE hben1,hben2,svec(3),dvec(3,4),dhben1(3,4),dhben2(3,4)
!
  rs1 = rspep(ps1)
  rs2 = rspep(ps2)
!  if (((seqpolty(rs1).ne.'P').AND.((seqtyp(rs1).gt.30).OR.(seqtyp(rs1).lt.27))).OR.&
! &    ((seqpolty(rs2).ne.'P').AND.((seqtyp(rs2).gt.30).OR.(seqtyp(rs2).lt.27))).OR.&
  if (abs(rs1-rs2).le.1) then
!   HB-assignments only work for pairs of polypeptides residues
    return
  end if
!
  qfac = -27.888 ! A*kcal/mol
!
! to be on the safe side ...
  if (rs1.gt.rs2) then
    rsbu = rs1
    rs1 = rs2
    rs2 = rsbu
  end if 
!
  svec(:) = 0.0
  call dis_bound_rs3(rs1,rs2,svec,dis)
! use distance screening to save CPU time
! this of course makes assumptions about the energy cutoff values for the
! H-bond formula below ...
  if (dis.gt.par_DSSP(4)) return
!
! we have to screen for PBC as well as for proline which cannot be a donor
!  if (hni(rs2).gt.0) then
  if ((hni(rs2).gt.0).AND.(oi(rs1).gt.0)) then
    call dis_bound5(oi(rs1),ni(rs2),svec,d2(1),dvec(:,1)) 
    call dis_bound5(oi(rs1),hni(rs2),svec,d2(2),dvec(:,2))
    call dis_bound5(ci(rs1),hni(rs2),svec,d2(3),dvec(:,3))
    call dis_bound5(ci(rs1),ni(rs2),svec,d2(4),dvec(:,4))
    id2(:) = 1./d2(:)
    id1(:) = sqrt(id2(:))
    hben1 = qfac*(id1(2) + id1(4) - id1(1) - id1(3))
    do j=1,4
      dhben1(:,j) = -qfac*dvec(:,j)*id2(j)*id1(j)
    end do
    dhben1(:,1) = -dhben1(:,1)
    dhben1(:,3) = -dhben1(:,3)
  else
    hben1 = 0.0 ! remember par_DSSP(2) is always negative
  end if
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = rs2
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
      if (hben1.lt.par_DSSP(3)) then
        hbl(peprs(rs1))%denahb(:,hbl(peprs(rs1))%nahbs) = 0.0
      else
!       force on the acceptor oxygen (donor atoms are stored in the donor HB)
        hbl(peprs(rs1))%denahb(1:3,hbl(peprs(rs1))%nahbs) = &
 &           - dhben1(:,1) - dhben1(:,2)
!       force on the acceptor carbon
        hbl(peprs(rs1))%denahb(4:6,hbl(peprs(rs1))%nahbs) = &
 &           - dhben1(:,3) - dhben1(:,4)
      end if
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = rs1
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
      if (hben1.lt.par_DSSP(3)) then
        hbl(peprs(rs2))%dendhb(:,hbl(peprs(rs2))%ndhbs) = 0.0
      else
!       force on the donor nitrogen (acc. atoms are stored in the acc. HB)
        hbl(peprs(rs2))%dendhb(1:3,hbl(peprs(rs2))%ndhbs) = &
 &            dhben1(:,1) + dhben1(:,4)
!       force on the donor hydrogen
        hbl(peprs(rs2))%dendhb(4:6,hbl(peprs(rs2))%ndhbs) = &
 &            dhben1(:,2) + dhben1(:,3)
      end if
    end if
!   update res-res-flag
    if (hb_mat(peprs(rs1),peprs(rs2)).gt.0) then
      hb_mat(peprs(rs1),peprs(rs2)) = 3
    else
      hb_mat(peprs(rs1),peprs(rs2)) = 1
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
!  if (hni(rs1).gt.0) then
  if ((hni(rs1).gt.0).AND.(oi(rs2).gt.0)) then
    call dis_bound5(ni(rs1),oi(rs2),svec,d2(1),dvec(:,1))
    call dis_bound5(hni(rs1),oi(rs2),svec,d2(2),dvec(:,2))
    call dis_bound5(hni(rs1),ci(rs2),svec,d2(3),dvec(:,3))
    call dis_bound5(ni(rs1),ci(rs2),svec,d2(4),dvec(:,4))
    id2(:) = 1./d2(:)
    id1(:) = sqrt(id2(:))
    hben2 = qfac*(id1(2) + id1(4) - id1(1) - id1(3))
    do j=1,4
      dhben2(:,j) = -qfac*dvec(:,j)*id2(j)*id1(j)
    end do
    dhben2(:,1) = -dhben2(:,1)
    dhben2(:,3) = -dhben2(:,3)
  else
    hben2 = 0.0
  end if
  if (hben2.lt.par_DSSP(2)) then
!   add donor hydrogen bond for rs1
    if (hbl(peprs(rs1))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%ndhbs = hbl(peprs(rs1))%ndhbs + 1
      hbl(peprs(rs1))%dhb(hbl(peprs(rs1))%ndhbs) = rs2
      hbl(peprs(rs1))%endhb(hbl(peprs(rs1))%ndhbs) = &
 &                                    max(hben2,par_DSSP(3))
      if (hben2.lt.par_DSSP(3)) then
        hbl(peprs(rs1))%dendhb(:,hbl(peprs(rs1))%ndhbs) = 0.0
      else
!       force on the donor nitrogen (acceptor atoms are stored in the acc. HB)
        hbl(peprs(rs1))%dendhb(1:3,hbl(peprs(rs1))%ndhbs) = &
 &            -dhben2(:,1) - dhben2(:,4)
!       force on the donor hydrogen
        hbl(peprs(rs1))%dendhb(4:6,hbl(peprs(rs1))%ndhbs) = &
 &            -dhben2(:,2) - dhben2(:,3)
      end if
    end if
!   add acceptor hydrogen bond for rs2
    if (hbl(peprs(rs2))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%nahbs = hbl(peprs(rs2))%nahbs + 1
      hbl(peprs(rs2))%ahb(hbl(peprs(rs2))%nahbs) = rs1
      hbl(peprs(rs2))%enahb(hbl(peprs(rs2))%nahbs) = &
 &                                    max(hben2,par_DSSP(3))
      if (hben2.lt.par_DSSP(3)) then
        hbl(peprs(rs2))%denahb(:,hbl(peprs(rs2))%nahbs) = 0.0
      else
!       force on the acceptor oxygen (donor atoms are stored in the donor HB)
        hbl(peprs(rs2))%denahb(1:3,hbl(peprs(rs2))%nahbs) = &
 &             dhben2(:,1) + dhben2(:,2)
!       force on the acceptor carbon
        hbl(peprs(rs2))%denahb(4:6,hbl(peprs(rs2))%nahbs) = &
 &             dhben2(:,3) + dhben2(:,4)
      end if
    end if
!   update res-res-flag
    if (hb_mat(peprs(rs1),peprs(rs2)).gt.0) then
      hb_mat(peprs(rs1),peprs(rs2)) = 3
    else
      hb_mat(peprs(rs1),peprs(rs2)) = 2
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine assign_bridge(pps1,pps2)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: pps1,pps2
!
  integer ps1,ps2,i,hbrpty
  RTYPE eant,epar
  logical foundit
!
  ps1 = min(pps1,pps2)
  ps2 = max(pps1,pps2)
!
  if ((ps1.eq.1).OR.(ps1.eq.pep_sz).OR.&
 &    (ps2.eq.1).OR.(ps2.eq.pep_sz)) return
!
  if ((molofrs(rspep(ps1)).ne.molofrs(rspep(ps1+1))).OR.&
 &    (molofrs(rspep(ps1)).ne.molofrs(rspep(ps1-1))).OR.&
 &    (molofrs(rspep(ps2)).ne.molofrs(rspep(ps2+1))).OR.&
 &    (molofrs(rspep(ps2)).ne.molofrs(rspep(ps2-1)))) return
!
  if (abs(ps1-ps2).lt.2.9) then
    write(ilog,*) 'Fatal. Residues have to be at least separated by &
 &two more peptide(!) residues in assign_bridge(...).'
    call fexit()
  end if
!
  hbrpty = 0
  eant = 0.0
  epar = 0.0
!
! the algorithm below will only work if polypeptide residues are
! continuous, i.e., there's no backbone mutations within a given
! chain
  if ((rspep(ps1).gt.rsmol(molofrs(rspep(ps1)),1)).AND.&
 &    (rspep(ps2).lt.rsmol(molofrs(rspep(ps2)),2))) then
!   and antiparallel sheet is described by either the residue being an
!   H-bonded one (reciprocal) or by the surrounding residues having to make
!   the proper H-bonds
    if (((hb_mat(ps1+1,ps2-1).ge.2).AND.&
 &   ((hb_mat(ps1-1,ps2+1).eq.1).OR.(hb_mat(ps1-1,ps2+1).eq.3))).OR.&
 &   (hb_mat(ps1,ps2).eq.3)) then
!      write(*,*) 'found ',rspep(ps1),rspep(ps2),' anti'
      hbrpty = 2
!   a parallel sheet is described by alternate registry involving the next
!   neighbors on either strand
    else if (((hb_mat(ps1+1,ps2).ge.2).AND.&
 &  ((hb_mat(ps1-1,ps2).eq.1).OR.(hb_mat(ps1-1,ps2).eq.3))).OR.&
 &           ((hb_mat(ps1,ps2-1).ge.2).AND.&
 &  ((hb_mat(ps1,ps2+1).eq.1).OR.(hb_mat(ps1,ps2+1).eq.3)))) then
      hbrpty = 1
!      write(*,*) 'found ',rspep(ps1),rspep(ps2),' para'
    else
!     no sheet residue obviously
      return
    end if
  else
    return
  end if
!
! now if possible collapse the new-found sheet-element into an existing bridge structure
  foundit = .false.
  do i=1,bdgtb%nrs
!   i is the right continuation residue and strand-pair is correct orientation?
    if ((bdgtb%lims(i,2)+1.eq.ps1).AND.&
 &      (hbrpty.eq.bdgtb%ty(i))) then
!     j is the right continuation residue for the two cases?
      if (hbrpty.eq.2) then
        if (ps2.eq.bdgtb%lims(i,3)-1) then
          foundit = .true.
          bdgtb%lims(i,2) = bdgtb%lims(i,2) + 1
          bdgtb%lims(i,3) = bdgtb%lims(i,3) - 1
          if (dssp_ass(ps1,1).eq.1) then
            dssp_ass(ps1,1) = 4
          else if (dssp_ass(ps1,1).eq.2) then
            dssp_ass(ps1,1) = 5
          else
            dssp_ass(ps1,1) = 2
          end if
          if (dssp_ass(ps2,1).eq.1) then
            dssp_ass(ps2,1) = 4
          else if (dssp_ass(ps2,1).eq.2) then
            dssp_ass(ps2,1) = 5
          else
            dssp_ass(ps2,1) = 2              
          end if
        end if
      else if (hbrpty.eq.1) then
        if (ps2.eq.bdgtb%lims(i,4)+1) then
          foundit = .true.
          bdgtb%lims(i,2) = bdgtb%lims(i,2) + 1
          bdgtb%lims(i,4) = bdgtb%lims(i,4) + 1
          if (dssp_ass(ps1,1).eq.1) then
            dssp_ass(ps1,1) = 3
          else if (dssp_ass(ps1,1).eq.2) then
            dssp_ass(ps1,1) = 4
          else
            dssp_ass(ps1,1) = 1
          end if
          if (dssp_ass(ps2,1).eq.1) then
            dssp_ass(ps2,1) = 3
          else if (dssp_ass(ps2,1).eq.2) then
            dssp_ass(ps2,1) = 4
          else
            dssp_ass(ps2,1) = 1              
          end if
        end if
      end if
    end if
  end do
!
! if not possible, make a new one
  if (foundit.EQV..false.) then
    bdgtb%nrs = bdgtb%nrs + 1
    bdgtb%lims(bdgtb%nrs,1) = ps1
    bdgtb%lims(bdgtb%nrs,2) = ps1
    bdgtb%lims(bdgtb%nrs,3) = ps2
    bdgtb%lims(bdgtb%nrs,4) = ps2
    bdgtb%ty(bdgtb%nrs) = hbrpty
    if (hbrpty.eq.1) then
      dssp_ass(ps1,1) = 1
      dssp_ass(ps2,1) = 1
    else if (hbrpty.eq.2) then
      dssp_ass(ps1,1) = 2
      dssp_ass(ps2,1) = 2
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this subroutine accumulates the total hydrogen bond energy for each bridge
! it fundamentally relies on the assumed structure of the bridgetabel for
! things like the strand in lims(1:2) always being N-terminal to the one
! in lims(3:4)
!
subroutine weigh_bridges(esc,wesc,escm,wescm,normby)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer ps1,ps2,i,j,allhb,imol,pmol1,pmol2,normby
  RTYPE en1,esc,wesc,escm(npepmol),wescm(npepmol),allhbm(npepmol)
  logical fi(3)
  esc = 0.0
  escm(:) = 0.0
  wesc = 0.0
  wescm(:) = 0.0
  allhb = 0
  allhbm(:) = 0.0
!
  do i=1,bdgtb%nrs
!    write(*,*) i,':'
!    write(*,*) bdgtb%ty(i)
!    write(*,*) bdgtb%lims(i,1:2)
!    write(*,*) bdgtb%lims(i,3:4)
    en1 = 0.0
    if (bdgtb%ty(i).eq.2) then
      do ps1=bdgtb%lims(i,1),bdgtb%lims(i,2)
        ps2 = bdgtb%lims(i,4)-(ps1-bdgtb%lims(i,1))
        if (ps1.gt.ps2) then
          write(ilog,*) 'Fatal. Bridge structure was created with wrong&
 & strand order. This is most likely a bug in assign_bridge().'
          call fexit()
        end if
        if (hb_mat(ps1,ps2).eq.3) then
!         directly bonded residue pair in antiparallel sheet
          fi(:) = .false.
          do j=1,hbl(ps1)%nahbs
            if (hbl(ps1)%ahb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%enahb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1)%enahb(j)
              end if
              fi(1) = .true.
            end if
          end do
          do j=1,hbl(ps1)%ndhbs
            if (hbl(ps1)%dhb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%endhb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1)%endhb(j)
              end if
              fi(2) = .true.
            end if
          end do
          if ((fi(1).EQV..false.).OR.(fi(2).EQV..false.)) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
!       non-bonded residue pair: only find half-contribution if terminal
        else if (ps1.eq.bdgtb%lims(i,1)) then
          fi(3) = .false.
          do j=1,hbl(ps1-1)%nahbs
            if (hbl(ps1-1)%ahb(j).eq.rspep(ps2+1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1-1)%enahb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1-1)%enahb(j)
              end if
              fi(3) = .true.
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        else if (ps1.eq.bdgtb%lims(i,2)) then
          fi(3) = .false.
          do j=1,hbl(ps1+1)%ndhbs
            if (hbl(ps1+1)%dhb(j).eq.rspep(ps2-1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1+1)%endhb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1+1)%endhb(j)
              end if
              fi(3) = .true.
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        end if
      end do
    else if (bdgtb%ty(i).eq.1) then
      do ps1=bdgtb%lims(i,1),bdgtb%lims(i,2)
        ps2 = bdgtb%lims(i,3)+(ps1-bdgtb%lims(i,1))
        if (ps1.gt.ps2) then
          write(ilog,*) 'Fatal. Bridge structure was created with wrong&
 & strand order. This is most likely a bug in assign_bridge().'
          call fexit()
        end if
        if (((hb_mat(ps1,ps2+1).eq.1).OR.(hb_mat(ps1,ps2+1).eq.3))&
 &     .AND.(hb_mat(ps1,ps2-1).ge.2)) then
!         residue in parallel sheet for which polar atoms of ps1 point toward other strand
          fi(:) = .false.
          do j=1,hbl(ps1)%nahbs
            if (hbl(ps1)%ahb(j).eq.rspep(ps2+1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%enahb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1)%enahb(j)
              end if
              fi(1) = .true.
            end if
          end do
          do j=1,hbl(ps1)%ndhbs
            if (hbl(ps1)%dhb(j).eq.rspep(ps2-1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%endhb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1)%endhb(j)
              end if
              fi(2) = .true.
            end if
          end do
          if ((fi(1).EQV..false.).OR.(fi(2).EQV..false.)) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
!        else if (((hb_mat(ps1,ps2+1).eq.1)
! &            .OR.(hb_mat(ps1,ps2+1).eq.3))
! &          .AND.(hb_mat(ps1,ps2-1).ge.2)) then
!         residue in parallel sheet for which polar atoms of ps1 point away from the other strand

!       non-bonded residue pair: only find half-contribution if terminal
        else if (ps1.eq.bdgtb%lims(i,1)) then
          fi(3) = .false.
          do j=1,hbl(ps1-1)%nahbs
            if (hbl(ps1-1)%ahb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1-1)%enahb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1-1)%enahb(j)
              end if
              fi(3) = .true.
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        else if (ps1.eq.bdgtb%lims(i,2)) then
          fi(3) = .false.
          do j=1,hbl(ps1+1)%ndhbs
            if (hbl(ps1+1)%dhb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1+1)%endhb(j),par_DSSP(1))
              else
                en1 = en1 + hbl(ps1+1)%endhb(j)
              end if
              fi(3) = .true.
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        end if
      end do
    end if
    wesc = wesc + en1
    allhb = allhb + (bdgtb%lims(i,2)-bdgtb%lims(i,1)+2)
    pmol1 = pepmol(molofrs(rspep(bdgtb%lims(i,1))))
    pmol2 = pepmol(molofrs(rspep(bdgtb%lims(i,3))))
    if (pmol1.eq.pmol2) then
      wescm(pmol1) = wescm(pmol1) + en1
      allhbm(pmol1)= allhbm(pmol1) + 1.0*&
 &            (bdgtb%lims(i,2) - bdgtb%lims(i,1) + 2)
    else
!     for cross-chain bridges, split symmetrically 
      wescm(pmol1) = wescm(pmol1) + 0.5*en1
      allhbm(pmol1)= allhbm(pmol1) + 0.5*&
 &            (bdgtb%lims(i,2) - bdgtb%lims(i,1) + 2)
      wescm(pmol2) = wescm(pmol2) + 0.5*en1
      allhbm(pmol2)= allhbm(pmol2) + 0.5*&
 &            (bdgtb%lims(i,2) - bdgtb%lims(i,1) + 2)
    end if
    bdgtb%phbe(i) = en1/((bdgtb%lims(i,2)-bdgtb%lims(i,1)+2))
!    write(*,*) en1/(1.0*(bdgtb%lims(i,2)-bdgtb%lims(i,1)+2))
  end do
  do i=1,pep_sz
    if ((dssp_ass(i,1).ge.1).AND.(dssp_ass(i,1).le.5)) then
      esc = esc + 1.0
      pmol1 = pepmol(molofrs(rspep(i)))
      escm(pmol1) = escm(pmol1) + 1.0
    end if
  end do
  if (normby.gt.0) esc = esc/(1.0*normby)
  if (allhb.gt.0) then
    if (par_DSSP2(2).eq.2) then
      wesc = min(1.0d0,wesc/(1.0*par_DSSP(1)*allhb))
    else
      wesc = wesc/(1.0*par_DSSP(1)*allhb)
    end if
  else
    wesc = 0.0
  end if
!
  do imol=1,npepmol
    if (allhbm(imol).gt.0.0) then 
      if (par_DSSP2(2).eq.2) then
        wescm(imol) = min(1.0d0,wescm(imol)/(par_DSSP(1)*allhbm(imol)))
      else
        wescm(imol) = wescm(imol)/(par_DSSP(1)*allhbm(imol))
      end if
    else
      wescm(imol) = 0.0
    end if
    escm(imol) = escm(imol)/(1.0*pepmolsz(imol))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! same thing only for derivatives as well (in wesc_der)
!
subroutine weigh_bridges_der(esc,wesc,escm,wescm,normby)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer ps1,ps2,i,j,l,allhb,imol,pmol1,pmol2,normby
  RTYPE en1,esc,wesc,escm(npepmol),wescm(npepmol),allhbm(npepmol),dum
  logical fi(3)
  esc = 0.0
  escm(:) = 0.0
  wesc = 0.0
  wescm(:) = 0.0
  allhb = 0
  allhbm(:) = 0.0
  wesc_der(:,:,:) = 0.0
!
  do i=1,bdgtb%nrs
!    write(*,*) i,':'
!    write(*,*) bdgtb%ty(i)
!    write(*,*) bdgtb%lims(i,1:2)
!    write(*,*) bdgtb%lims(i,3:4)
    en1 = 0.0
    if (bdgtb%ty(i).eq.2) then
      do ps1=bdgtb%lims(i,1),bdgtb%lims(i,2)
        ps2 = bdgtb%lims(i,4)-(ps1-bdgtb%lims(i,1))
        if (ps1.gt.ps2) then
          write(ilog,*) 'Fatal. Bridge structure was created with wrong&
 & strand order. This is most likely a bug in assign_bridge().'
          call fexit()
        end if
        if (hb_mat(ps1,ps2).eq.3) then
!         directly bonded residue pair in antiparallel sheet
          fi(:) = .false.
          do j=1,hbl(ps1)%nahbs
            if (hbl(ps1)%ahb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%enahb(j),par_DSSP(1))
                if (hbl(ps1)%enahb(j).ge.par_DSSP(1)) then
                  wesc_der(:,1,ps1) = wesc_der(:,1,ps1) +&
 &                   hbl(ps1)%denahb(1:3,j)
                  wesc_der(:,2,ps1) = wesc_der(:,2,ps1) +&
 &                   hbl(ps1)%denahb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1)%enahb(j)
                wesc_der(:,1,ps1) = wesc_der(:,1,ps1) +&
 &                   hbl(ps1)%denahb(1:3,j)
                wesc_der(:,2,ps1) = wesc_der(:,2,ps1) +&
 &                   hbl(ps1)%denahb(4:6,j)
              end if
              fi(1) = .true.
              do l=1,hbl(ps2)%ndhbs
                if (hbl(ps2)%dhb(l).eq.rspep(ps1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2)%endhb(l).ge.par_DSSP(1)) then
                      wesc_der(:,3,ps2) = wesc_der(:,3,ps2) +&
 &                   hbl(ps2)%dendhb(1:3,l)
                      wesc_der(:,4,ps2) = wesc_der(:,4,ps2) +&
 &                   hbl(ps2)%dendhb(4:6,l)
                    end if
                  else
                    wesc_der(:,3,ps2) = wesc_der(:,3,ps2) +&
 &                   hbl(ps2)%dendhb(1:3,l)
                    wesc_der(:,4,ps2) = wesc_der(:,4,ps2) +&
 &                   hbl(ps2)%dendhb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          do j=1,hbl(ps1)%ndhbs
            if (hbl(ps1)%dhb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%endhb(j),par_DSSP(1))
                if (hbl(ps1)%endhb(j).ge.par_DSSP(1)) then
                  wesc_der(:,3,ps1) = wesc_der(:,3,ps1) +&
 &                   hbl(ps1)%dendhb(1:3,j)
                  wesc_der(:,4,ps1) = wesc_der(:,4,ps1) +&
 &                   hbl(ps1)%dendhb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1)%endhb(j)
                wesc_der(:,3,ps1) = wesc_der(:,3,ps1) +&
 &                   hbl(ps1)%dendhb(1:3,j)
                wesc_der(:,4,ps1) = wesc_der(:,4,ps1) +&
 &                   hbl(ps1)%dendhb(4:6,j)
              end if
              fi(2) = .true.
              do l=1,hbl(ps2)%nahbs
                if (hbl(ps2)%ahb(l).eq.rspep(ps1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2)%enahb(l).ge.par_DSSP(1)) then
                      wesc_der(:,1,ps2) = wesc_der(:,1,ps2) +&
 &                   hbl(ps2)%denahb(1:3,l)
                      wesc_der(:,2,ps2) = wesc_der(:,2,ps2) +&
 &                   hbl(ps2)%denahb(4:6,l)
                    end if
                  else
                    wesc_der(:,1,ps2) = wesc_der(:,1,ps2) +&
 &                   hbl(ps2)%denahb(1:3,l)
                    wesc_der(:,2,ps2) = wesc_der(:,2,ps2) +&
 &                   hbl(ps2)%denahb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          if ((fi(1).EQV..false.).OR.(fi(2).EQV..false.)) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
!       non-bonded residue pair: only find half-contribution if terminal
        else if (ps1.eq.bdgtb%lims(i,1)) then
          fi(3) = .false.
          do j=1,hbl(ps1-1)%nahbs
            if (hbl(ps1-1)%ahb(j).eq.rspep(ps2+1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1-1)%enahb(j),par_DSSP(1))
                if (hbl(ps1-1)%enahb(j).ge.par_DSSP(1)) then
                  wesc_der(:,1,ps1-1) = wesc_der(:,1,ps1-1) +&
 &                   hbl(ps1-1)%denahb(1:3,j)
                  wesc_der(:,2,ps1-1) = wesc_der(:,2,ps1-1) +&
 &                   hbl(ps1-1)%denahb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1-1)%enahb(j)
                wesc_der(:,1,ps1-1) = wesc_der(:,1,ps1-1) +&
 &                   hbl(ps1-1)%denahb(1:3,j)
                wesc_der(:,2,ps1-1) = wesc_der(:,2,ps1-1) +&
 &                   hbl(ps1-1)%denahb(4:6,j)
              end if
              fi(3) = .true.
              do l=1,hbl(ps2+1)%ndhbs
                if (hbl(ps2+1)%dhb(l).eq.rspep(ps1-1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2+1)%endhb(l).ge.par_DSSP(1)) then
                      wesc_der(:,3,ps2+1) = wesc_der(:,3,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(1:3,l)
                      wesc_der(:,4,ps2+1) = wesc_der(:,4,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(4:6,l)
                    end if
                  else
                    wesc_der(:,3,ps2+1) = wesc_der(:,3,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(1:3,l)
                    wesc_der(:,4,ps2+1) = wesc_der(:,4,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        else if (ps1.eq.bdgtb%lims(i,2)) then
          fi(3) = .false.
          do j=1,hbl(ps1+1)%ndhbs
            if (hbl(ps1+1)%dhb(j).eq.rspep(ps2-1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1+1)%endhb(j),par_DSSP(1))
                if (hbl(ps1+1)%endhb(j).ge.par_DSSP(1)) then
                  wesc_der(:,3,ps1+1) = wesc_der(:,3,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(1:3,j)
                  wesc_der(:,4,ps1+1) = wesc_der(:,4,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1+1)%endhb(j)
                wesc_der(:,3,ps1+1) = wesc_der(:,3,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(1:3,j)
                wesc_der(:,4,ps1+1) = wesc_der(:,4,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(4:6,j)
              end if
              fi(3) = .true.
              do l=1,hbl(ps2-1)%nahbs
                if (hbl(ps2-1)%ahb(l).eq.rspep(ps1+1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2-1)%enahb(l).ge.par_DSSP(1)) then
                      wesc_der(:,1,ps2-1) = wesc_der(:,1,ps2-1) +&
 &                   hbl(ps2-1)%denahb(1:3,l)
                      wesc_der(:,2,ps2-1) = wesc_der(:,2,ps2-1) +&
 &                   hbl(ps2-1)%denahb(4:6,l)
                    end if
                  else
                    wesc_der(:,1,ps2-1) = wesc_der(:,1,ps2-1) +&
 &                   hbl(ps2-1)%denahb(1:3,l)
                    wesc_der(:,2,ps2-1) = wesc_der(:,2,ps2-1) +&
 &                   hbl(ps2-1)%denahb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        end if
      end do
    else if (bdgtb%ty(i).eq.1) then
      do ps1=bdgtb%lims(i,1),bdgtb%lims(i,2)
        ps2 = bdgtb%lims(i,3)+(ps1-bdgtb%lims(i,1))
        if (ps1.gt.ps2) then
          write(ilog,*) 'Fatal. Bridge structure was created with wrong&
 & strand order. This is most likely a bug in assign_bridge().'
          call fexit()
        end if
        if (((hb_mat(ps1,ps2+1).eq.1).OR.(hb_mat(ps1,ps2+1).eq.3))&
 &     .AND.(hb_mat(ps1,ps2-1).ge.2)) then
!         residue in parallel sheet for which polar atoms of ps1 point toward other strand
          fi(:) = .false.
          do j=1,hbl(ps1)%nahbs
            if (hbl(ps1)%ahb(j).eq.rspep(ps2+1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%enahb(j),par_DSSP(1))
                if (hbl(ps1)%enahb(j).ge.par_DSSP(1)) then
                  wesc_der(:,1,ps1) = wesc_der(:,1,ps1) +&
 &                   hbl(ps1)%denahb(1:3,j)
                  wesc_der(:,2,ps1) = wesc_der(:,2,ps1) +&
 &                   hbl(ps1)%denahb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1)%enahb(j)
                wesc_der(:,1,ps1) = wesc_der(:,1,ps1) +&
 &                   hbl(ps1)%denahb(1:3,j)
                wesc_der(:,2,ps1) = wesc_der(:,2,ps1) +&
 &                   hbl(ps1)%denahb(4:6,j)
              end if
              fi(1) = .true.
              do l=1,hbl(ps2+1)%ndhbs
                if (hbl(ps2+1)%dhb(l).eq.rspep(ps1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2+1)%endhb(l).ge.par_DSSP(1)) then
                      wesc_der(:,3,ps2+1) = wesc_der(:,3,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(1:3,l)
                      wesc_der(:,4,ps2+1) = wesc_der(:,4,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(4:6,l)
                    end if
                  else
                    wesc_der(:,3,ps2+1) = wesc_der(:,3,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(1:3,l)
                    wesc_der(:,4,ps2+1) = wesc_der(:,4,ps2+1) +&
 &                   hbl(ps2+1)%dendhb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          do j=1,hbl(ps1)%ndhbs
            if (hbl(ps1)%dhb(j).eq.rspep(ps2-1)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1)%endhb(j),par_DSSP(1))
                if (hbl(ps1)%endhb(j).ge.par_DSSP(1)) then
                  wesc_der(:,3,ps1) = wesc_der(:,3,ps1) +&
 &                   hbl(ps1)%dendhb(1:3,j)
                  wesc_der(:,4,ps1) = wesc_der(:,4,ps1) +&
 &                   hbl(ps1)%dendhb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1)%endhb(j)
                wesc_der(:,3,ps1) = wesc_der(:,3,ps1) +&
 &                   hbl(ps1)%dendhb(1:3,j)
                wesc_der(:,4,ps1) = wesc_der(:,4,ps1) +&
 &                   hbl(ps1)%dendhb(4:6,j)
              end if
              fi(2) = .true.
              do l=1,hbl(ps2-1)%nahbs
                if (hbl(ps2-1)%ahb(l).eq.rspep(ps1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2-1)%enahb(l).ge.par_DSSP(1)) then
                      wesc_der(:,1,ps2-1) = wesc_der(:,1,ps2-1) +&
 &                   hbl(ps2-1)%denahb(1:3,l)
                      wesc_der(:,2,ps2-1) = wesc_der(:,2,ps2-1) +&
 &                   hbl(ps2-1)%denahb(4:6,l)
                    end if
                  else
                    wesc_der(:,1,ps2-1) = wesc_der(:,1,ps2-1) +&
 &                   hbl(ps2-1)%denahb(1:3,l)
                    wesc_der(:,2,ps2-1) = wesc_der(:,2,ps2-1) +&
 &                   hbl(ps2-1)%denahb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          if ((fi(1).EQV..false.).OR.(fi(2).EQV..false.)) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
!        else if (((hb_mat(ps1,ps2+1).eq.1)
! &            .OR.(hb_mat(ps1,ps2+1).eq.3))
! &          .AND.(hb_mat(ps1,ps2-1).ge.2)) then
!         residue in parallel sheet for which polar atoms of ps1 point away from the other strand

!       non-bonded residue pair: only find half-contribution if terminal
        else if (ps1.eq.bdgtb%lims(i,1)) then
          fi(3) = .false.
          do j=1,hbl(ps1-1)%nahbs
            if (hbl(ps1-1)%ahb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1-1)%enahb(j),par_DSSP(1))
                if (hbl(ps1-1)%enahb(j).ge.par_DSSP(1)) then
                  wesc_der(:,1,ps1-1) = wesc_der(:,1,ps1-1) +&
 &                   hbl(ps1-1)%denahb(1:3,j)
                  wesc_der(:,2,ps1-1) = wesc_der(:,2,ps1-1) +&
 &                   hbl(ps1-1)%denahb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1-1)%enahb(j)
                wesc_der(:,1,ps1-1) = wesc_der(:,1,ps1-1) +&
 &                   hbl(ps1-1)%denahb(1:3,j)
                wesc_der(:,2,ps1-1) = wesc_der(:,2,ps1-1) +&
 &                   hbl(ps1-1)%denahb(4:6,j)
              end if
              fi(3) = .true.
              do l=1,hbl(ps2)%ndhbs
                if (hbl(ps2)%dhb(l).eq.rspep(ps1-1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2)%endhb(l).ge.par_DSSP(1)) then
                      wesc_der(:,3,ps2) = wesc_der(:,3,ps2) +&
 &                   hbl(ps2)%dendhb(1:3,l)
                      wesc_der(:,4,ps2) = wesc_der(:,4,ps2) +&
 &                   hbl(ps2)%dendhb(4:6,l)
                    end if
                  else
                    wesc_der(:,3,ps2) = wesc_der(:,3,ps2) +&
 &                   hbl(ps2)%dendhb(1:3,l)
                    wesc_der(:,4,ps2) = wesc_der(:,4,ps2) +&
 &                   hbl(ps2)%dendhb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        else if (ps1.eq.bdgtb%lims(i,2)) then
          fi(3) = .false.
          do j=1,hbl(ps1+1)%ndhbs
            if (hbl(ps1+1)%dhb(j).eq.rspep(ps2)) then
              if (par_DSSP2(2).eq.1) then
                en1 = en1 + max(hbl(ps1+1)%endhb(j),par_DSSP(1))
                if (hbl(ps1+1)%endhb(j).ge.par_DSSP(1)) then
                  wesc_der(:,3,ps1+1) = wesc_der(:,3,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(1:3,j)
                  wesc_der(:,4,ps1+1) = wesc_der(:,4,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(4:6,j)
                end if
              else
                en1 = en1 + hbl(ps1+1)%endhb(j)
                wesc_der(:,3,ps1+1) = wesc_der(:,3,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(1:3,j)
                wesc_der(:,4,ps1+1) = wesc_der(:,4,ps1+1) +&
 &                   hbl(ps1+1)%dendhb(4:6,j)
              end if
              fi(3) = .true.
              do l=1,hbl(ps2)%nahbs
                if (hbl(ps2)%ahb(l).eq.rspep(ps1+1)) then
                  if (par_DSSP2(2).eq.1) then
                    if (hbl(ps2)%enahb(l).ge.par_DSSP(1)) then
                      wesc_der(:,1,ps2) = wesc_der(:,1,ps2) +&
 &                   hbl(ps2)%denahb(1:3,l)
                      wesc_der(:,2,ps2) = wesc_der(:,2,ps2) +&
 &                   hbl(ps2)%denahb(4:6,l)
                    end if
                  else
                    wesc_der(:,1,ps2) = wesc_der(:,1,ps2) +&
 &                   hbl(ps2)%denahb(1:3,l)
                    wesc_der(:,2,ps2) = wesc_der(:,2,ps2) +&
 &                   hbl(ps2)%denahb(4:6,l)
                  end if
                end if
              end do
            end if
          end do
          if (fi(3).EQV..false.) then
            write(ilog,*) 'Fatal. Bridge-structure assembled is inco&
 &nsistent with H-bond pattern in weigh_bridges(). This is a bug.'
            call fexit()
          end if
        end if
      end do
    end if
    wesc = wesc + en1
    allhb = allhb + (bdgtb%lims(i,2)-bdgtb%lims(i,1)+2)
    pmol1 = pepmol(molofrs(rspep(bdgtb%lims(i,1))))
    pmol2 = pepmol(molofrs(rspep(bdgtb%lims(i,3))))
    if (pmol1.eq.pmol2) then
      wescm(pmol1) = wescm(pmol1) + en1
      allhbm(pmol1)= allhbm(pmol1) + 1.0*&
 &            (bdgtb%lims(i,2) - bdgtb%lims(i,1) + 2)
    else
!     for cross-chain bridges, split symmetrically 
      wescm(pmol1) = wescm(pmol1) + 0.5*en1
      allhbm(pmol1)= allhbm(pmol1) + 0.5*&
 &            (bdgtb%lims(i,2) - bdgtb%lims(i,1) + 2)
      wescm(pmol2) = wescm(pmol2) + 0.5*en1
      allhbm(pmol2)= allhbm(pmol2) + 0.5*&
 &            (bdgtb%lims(i,2) - bdgtb%lims(i,1) + 2)
    end if
    bdgtb%phbe(i) = en1/((bdgtb%lims(i,2)-bdgtb%lims(i,1)+2))
!    write(*,*) en1/(1.0*(bdgtb%lims(i,2)-bdgtb%lims(i,1)+2))
  end do
  do i=1,pep_sz
    if ((dssp_ass(i,1).ge.1).AND.(dssp_ass(i,1).le.5)) then
      esc = esc + 1.0
      pmol1 = pepmol(molofrs(rspep(i)))
      escm(pmol1) = escm(pmol1) + 1.0
    end if
  end do
  if (normby.gt.0) esc = esc/(1.0*normby)
    
  if (allhb.gt.0) then
    if (par_DSSP2(2).eq.2) then
      wesc = wesc/(1.0*par_DSSP(1)*allhb)
      if (wesc.le.1.0) then
        wesc_der(:,:,:) = wesc_der(:,:,:)/(1.0*par_DSSP(1)*allhb)
      else
        wesc = 1.0
        wesc_der(:,:,:) = 0.0
      end if
    else
      wesc = wesc/(1.0*par_DSSP(1)*allhb)
      wesc_der(:,:,:) = wesc_der(:,:,:)/(1.0*par_DSSP(1)*allhb)
    end if
  else
    wesc = 0.0
  end if
!
  do imol=1,npepmol
    if (allhbm(imol).gt.0.0) then 
      if (par_DSSP2(2).eq.2) then
        wescm(imol) = min(1.0d0,wescm(imol)/(par_DSSP(1)*allhbm(imol)))
      else
        wescm(imol) = wescm(imol)/(par_DSSP(1)*allhbm(imol))
      end if
    else
      wescm(imol) = 0.0
    end if
    escm(imol) = escm(imol)/(1.0*pepmolsz(imol))
  end do
!
! now generate the (more or less) smooth E-score from the fractional
! HB energies (depending on par_DSSP2(2)) and the actual E-fraction
  if (par_DSSP2(2).eq.3) then
    dum = min(1.0d0,(wesc**(1.0/par_DSSP2(1)))*esc)
    if (dum.gt.1.0) then
      wesc_der(:,:,:) = 0.0
    else if (dum.gt.0.0) then
      wesc_der(:,:,:) = (1.0/par_DSSP2(1))*(dum/wesc)*wesc_der(:,:,:)
    end if
    wesc = dum
  else
    if (wesc.gt.0.0) then
      wesc_der(:,:,:) = esc*(1.0/par_DSSP2(1))*(wesc**((1.0/par_DSSP2(1)) - 1.0))*wesc_der(:,:,:)
      wesc = (wesc**(1.0/par_DSSP2(1)))*esc
    end if
  end if
!
  do i=1,npepmol
    if (par_DSSP2(2).eq.3) then
      wescm(i) = min(1.0d0,(wescm(i)**(1.0/par_DSSP2(1)))*escm(i))
    else
      wescm(i) = (wescm(i)**(1.0/par_DSSP2(1)))*escm(i)
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine assign_helixturns()
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer i,j,k,kk,klst(3),hmk(3)
  logical dropout
!
  klst(1) = 4
  klst(2) = 3
  klst(3) = 5
  hmk(1) = 6
  hmk(2) = 7
  hmk(3) = 8
!
! assign turns
  do i=1,pep_sz
    do k=3,5
      if (i.le.pep_sz-k) then
        if (molofrs(rspep(i)).eq.molofrs(rspep(i+k))) then
!         test for proper hydrogen bond
          if ((hb_mat(i,i+k).eq.1).OR.(hb_mat(i,i+k).eq.3)) then
            dssp_ass(i+k,k-1) = 2
            do j=i+1,i+k-1
              if (dssp_ass(j,k-1).eq.0) then
                dssp_ass(j,k-1) = 1
              end if
            end do
            if (dssp_ass(i,k-1).eq.2) then
              dssp_ass(i,k-1) = 4
            else
              dssp_ass(i,k-1) = 3
            end if
          end if
        end if
      end if
    end do
  end do
!
! now assign helices
  do i=1,pep_sz
    do kk=1,3
      k = klst(kk)
      if ((i.le.pep_sz-k+1).AND.(i.ge.2)) then
        if (molofrs(rspep(i-1)).eq.molofrs(rspep(i+k-1))) then
          if ((dssp_ass(i,k-1).ge.3).AND.(dssp_ass(i-1,k-1).ge.3)) then
            dropout = .false.
            do j=i,i+k-1
              if ((dssp_ass(j,1).ne.hmk(kk)).AND.&
 &                              (dssp_ass(j,1).ne.0)) then
                dropout = .true.
                exit
              end if
            end do
            if (dropout.EQV..true.) exit
            do j=i,i+k-1
              dssp_ass(j,1) = hmk(kk)
            end do
!           now exit the kk-loop since other helix assignments have
!           lower priority
            exit
          end if
        end if
      end if
    end do
  end do
!
! and finally classify remaining residues as either turn or bend
  do i=1,pep_sz
    do k=3,5
      if ((i.le.pep_sz-k).AND.(dssp_ass(i,1).eq.0)) then
        if (molofrs(rspep(i)).eq.molofrs(rspep(i+k))) then
          do j=1,k-1
            if (i.gt.j) then
              if (dssp_ass(i-j,k-1).ge.3) then
                dssp_ass(i,1) = 9
              end if
            end if
          end do
        end if
      end if
    end do
    if (dssp_ass(i,1).eq.0) then
      dssp_ass(i,1) = dssp_ass(i,5) ! i.e., 10 (S)
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine weigh_helices(hsc,whsc,hscm,whscm,normby)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer i,j,k,hstart,hend,hnb,allhb,hseg,inhs,imol,pmol,normby
  integer inhsm(npepmol),allhbm(npepmol)
  RTYPE hhbe,hsc,whsc,hscm(npepmol),whscm(npepmol)
!
! now assign helices
  hseg = 0
  inhs = 0
  allhb = 0
  allhbm(:) = 0
  whscm(:) = 0.0
  whsc = 0.0
  inhs = 0
  inhsm(:) = 0
  do i=1,pep_sz
    pmol = pepmol(molofrs(rspep(i)))
!   helix segments never extend into terminal residues
    if ((hseg.eq.0).AND.(dssp_ass(i,1).eq.6)) then
      hstart = i
      hhbe = 0.0
      hnb = 0
      hseg = 1
    else if ((hseg.eq.1).AND.(dssp_ass(i,1).ne.6)) then
      hend = i-1
      do k=hstart-1,hend-3
!        foundit = .false.
        do j=1,hbl(k)%nahbs
          if (hbl(k)%ahb(j).eq.rspep(k+4)) then
            if (par_DSSP2(2).eq.1) then
              hhbe = hhbe + max(hbl(k)%enahb(j),par_DSSP(1))
            else
              hhbe = hhbe + hbl(k)%enahb(j)
            end if
!            foundit = .true.
            exit
          end if
        end do
!       we could use this to diagnose integrity of helices but don't currently
!        if (foundit.EQV..false.) then
!          write(ilog,*) 'Warning: Missing i->i+4 HB at ',rspep(k)
!        end if
      end do
      whsc = whsc + hhbe
      whscm(pmol) = whscm(pmol) + hhbe
      allhb = allhb + hend - 1 - hstart
      allhbm(pmol) = allhbm(pmol) + hend - 1 - hstart
      inhs = inhs + hend - hstart + 1
      inhsm(pmol) = inhsm(pmol) + hend - hstart + 1
      hhbe = 0.0
      hnb = 0
      hseg = 0
    end if
  end do
!
  if (allhb.gt.0) then 
    if (par_DSSP2(2).eq.2) then
      whsc = min(1.0d0,whsc/(1.0*par_DSSP(1)*allhb))
    else
      whsc = whsc/(1.0*par_DSSP(1)*allhb)
    end if
  else
    whsc = 0.0
  end if
  if (normby.gt.0) hsc = 1.0*inhs/(1.0*normby)
!
  do imol=1,npepmol
    if (allhbm(imol).gt.0) then 
      if (par_DSSP2(2).eq.2) then
        whscm(imol) = min(1.0d0,whscm(imol)/(1.0*par_DSSP(1)*allhbm(imol)))
      else
        whscm(imol) = whscm(imol)/(1.0*par_DSSP(1)*allhbm(imol))
      end if
    else
      whscm(imol) = 0.0
    end if
    hscm(imol) = 1.0*inhsm(imol)/(1.0*pepmolsz(imol))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! note that this fxn initializes the global whsc_der to zero
! it provides the derivatives of whsc with respect to the positions
! of the relevant donor and acceptor atoms on all the peptide residues
!
subroutine weigh_helices_der(hsc,whsc,hscm,whscm,normby)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer i,j,k,l,hstart,hend,hnb,allhb,hseg,inhs,imol,pmol,normby
  integer inhsm(npepmol),allhbm(npepmol)
  RTYPE hhbe,hsc,whsc,hscm(npepmol),whscm(npepmol),dum
!
! now assign helices
  hseg = 0
  inhs = 0
  allhb = 0
  allhbm(:) = 0
  whscm(:) = 0.0
  whsc = 0.0
  inhs = 0
  inhsm(:) = 0
  whsc_der(:,:,:) = 0.0
  do i=1,pep_sz
    pmol = pepmol(molofrs(rspep(i)))
!   helix segments never extend into terminal residues
    if ((hseg.eq.0).AND.(dssp_ass(i,1).eq.6)) then
      hstart = i
      hhbe = 0.0
      hnb = 0
      hseg = 1
    else if ((hseg.eq.1).AND.(dssp_ass(i,1).ne.6)) then
      hend = i-1
      do k=hstart-1,hend-3
!        foundit = .false.
        do j=1,hbl(k)%nahbs
          if (hbl(k)%ahb(j).eq.rspep(k+4)) then
            if (par_DSSP2(2).eq.1) then
              hhbe = hhbe + max(hbl(k)%enahb(j),par_DSSP(1))
              if (hbl(k)%enahb(j).ge.par_DSSP(1)) then
                whsc_der(:,1,k) = whsc_der(:,1,k) + &
 &                    hbl(k)%denahb(1:3,j)
                whsc_der(:,2,k) = whsc_der(:,2,k) + &
 &                    hbl(k)%denahb(4:6,j)
              end if
            else
              hhbe = hhbe + hbl(k)%enahb(j)
              whsc_der(:,1,k) = whsc_der(:,1,k) + &
 &                    hbl(k)%denahb(1:3,j)
              whsc_der(:,2,k) = whsc_der(:,2,k) + &
 &                    hbl(k)%denahb(4:6,j)
            end if
            do l=1,hbl(k+4)%ndhbs
              if (hbl(k+4)%dhb(l).eq.rspep(k)) then
                if (par_DSSP2(2).eq.1) then
                  if (hbl(k+4)%endhb(l).ge.par_DSSP(1)) then
                    whsc_der(:,3,k+4) = whsc_der(:,3,k+4) + &
 &                    hbl(k+4)%dendhb(1:3,l)
                    whsc_der(:,4,k+4) = whsc_der(:,4,k+4) + &
 &                    hbl(k+4)%dendhb(4:6,l)
                  end if
                else
                  whsc_der(:,3,k+4) = whsc_der(:,3,k+4) + &
 &                    hbl(k+4)%dendhb(1:3,l)
                  whsc_der(:,4,k+4) = whsc_der(:,4,k+4) + &
 &                    hbl(k+4)%dendhb(4:6,l)
                end if
              end if
            end do
!            foundit = .true.
            exit
          end if
        end do
!       we could use this to diagnose integrity of helices but don't currently
!        if (foundit.EQV..false.) then
!          write(ilog,*) 'Warning: Missing i->i+4 HB at ',rspep(k)
!        end if
      end do
      whsc = whsc + hhbe
      whscm(pmol) = whscm(pmol) + hhbe
      allhb = allhb + hend - 1 - hstart
      allhbm(pmol) = allhbm(pmol) + hend - 1 - hstart
      inhs = inhs + hend - hstart + 1
      inhsm(pmol) = inhsm(pmol) + hend - hstart + 1
      hhbe = 0.0
      hnb = 0
      hseg = 0
    end if
  end do
!
  if (allhb.gt.0) then 
    if (par_DSSP2(2).eq.2) then
      whsc = whsc/(1.0*par_DSSP(1)*allhb)
      if (whsc.gt.1.0) then
        whsc_der(:,:,:) = 0.0
        whsc = 1.0
      else
        whsc_der(:,:,:) = whsc_der(:,:,:)/(1.0*par_DSSP(1)*allhb)
      end if
    else
      whsc = whsc/(1.0*par_DSSP(1)*allhb)
      whsc_der(:,:,:) = whsc_der(:,:,:)/(1.0*par_DSSP(1)*allhb)
    end if
  else
    whsc = 0.0
  end if
  if (normby.gt.0) hsc = 1.0*inhs/(1.0*normby)
!
  do imol=1,npepmol
    if (allhbm(imol).gt.0) then 
      if (par_DSSP2(2).eq.2) then
        whscm(imol) = min(1.0d0,whscm(imol)/(1.0*par_DSSP(1)*allhbm(imol)))
      else
        whscm(imol) = whscm(imol)/(1.0*par_DSSP(1)*allhbm(imol))
      end if
    else
      whscm(imol) = 0.0
    end if
    hscm(imol) = 1.0*inhsm(imol)/(1.0*pepmolsz(imol))
  end do
!
! now generate the (more or less) smooth H-score from the fractional
! HB energies (depending on par_DSSP2(2)) and the actual H-fraction
  if (par_DSSP2(2).eq.3) then
    dum = min(1.0d0,(whsc**(1.0/par_DSSP2(1)))*hsc)
    if (dum.gt.1.0) then
      whsc_der(:,:,:) = 0.0
    else if (dum.gt.0.0) then
      whsc_der(:,:,:) = (1.0/par_DSSP2(1))*(dum/whsc)*whsc_der(:,:,:)
    end if
    whsc = dum
  else
    if (whsc.gt.0.0) then
      whsc_der(:,:,:) = hsc*(1.0/par_DSSP2(1))*(whsc**((1.0/par_DSSP2(1)) - 1.0))*whsc_der(:,:,:)
      whsc = (whsc**(1.0/par_DSSP2(1)))*hsc
    end if
  end if
!
  do i=1,npepmol
    if (par_DSSP2(2).eq.3) then
      whscm(i) = min(1.0d0,(whscm(i)**(1.0/par_DSSP2(1)))*hscm(i))
    else
      whscm(i) = (whscm(i)**(1.0/par_DSSP2(1)))*hscm(i)
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine assign_curves()
!
  use iounit
  use atoms
  use polypep
  use math
  use sequen
  use dssps
  use molecule
!
  implicit none
!
  integer i
  RTYPE tor,getztor,ang,getbang
!
  do i=2,pep_sz-2
!   bend assignment needs to be cycled for termini incl. special exclusion for FOR/NH2 caps 
    if (molofrs(rspep(i-1)).eq.molofrs(rspep(i+2)).AND.(cai(rspep(i+2)).gt.0).AND.(cai(rspep(i-1)).gt.0)) then
      tor = getztor(cai(rspep(i-1)),cai(rspep(i)),&
 &                  cai(rspep(i+1)),cai(rspep(i+2)))
      if (tor.lt.0.0) then
        dssp_ass(i,6) = 1
      else
        dssp_ass(i,6) = 2
      end if
      if (i.gt.2) then
        if (cai(rspep(i-2)).gt.0) then
          ang = getbang(cai(rspep(i-2)),cai(rspep(i)),&
 &                    cai(rspep(i+2)))
          if (ang.lt.110.0) then
            dssp_ass(i,5) = 10
          end if
        end if
      end if
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine update_dssp(esc,wesc,esch,hsc,whsc,hsch,escm,wescm,hscm,whscm,intraflag)
!
  use iounit
  use sequen
  use dssps
  use atoms
  use polypep
  use mcsums
  use system
  use grandensembles
  use molecule
!
  implicit none
!
  integer i,j,ps1,ps2,normby,lastp,pmol
  RTYPE esch,hsch
  RTYPE esc,escm(npepmol),wesc,wescm(npepmol)
  RTYPE hsc,hscm(npepmol),whsc,whscm(npepmol)
  logical intraflag,ismember
!
  hb_mat(:,:) = 0
  hbl(:)%nahbs = 0
  hbl(:)%ndhbs = 0
  dssp_ass(:,:) = 0
!
! now go through the various sub-analyses necessary for DSSP ...
!
  call assign_curves()
!
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    if (intraflag.EQV..true.) then
      normby = pep_sz
      do ps1=1,pep_sz
        pmol = pepmol(molofrs(rspep(ps1)))
        do ps2=ps1+1,pep_sz
          if (molofrs(rspep(ps2)).ne.pmol) exit
          call get_bbhb(ps1,ps2)
        end do
        lastp = ps2
        if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
          if (ismember(ispresent,pmol).EQV..false.) cycle
        end if
        do ps2=lastp,pep_sz
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps2)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps2))).EQV..false.) cycle
          end if
          call get_bbhb(ps1,ps2)
        end do
      end do
    else
      normby = 0
      do ps1=1,pep_sz
        pmol = pepmol(molofrs(rspep(ps1)))
        if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
          if (ismember(ispresent,pmol).EQV..false.) cycle
        end if
        normby = normby + 1
        do ps2=ps1+1,pep_sz
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps2)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps2))).EQV..false.) cycle
          end if
          call get_bbhb(ps1,ps2)
        end do
      end do
    end if
  else
    normby = pep_sz
    do ps1=1,pep_sz
      do ps2=ps1+1,pep_sz
        call get_bbhb(ps1,ps2)
      end do
    end do
  end if
!
! now assemble the bridge structure
  bdgtb%nrs = 0
  do i=1,pep_sz
    do j=i+3,pep_sz
      call assign_bridge(i,j)
    end do
  end do
!
  call weigh_bridges(esc,esch,escm,wescm,normby)
!
  call assign_helixturns()
!
  call weigh_helices(hsc,hsch,hscm,whscm,normby)
!
! now generate the (more or less) smooth E/H-scores from the fractional
! HB energies (depending on par_DSSP2(2)) and the actual E/H-fractions
  if (par_DSSP2(2).eq.3) then
    wesc = min(1.0d0,(esch**(1.0/par_DSSP2(1)))*esc)
    whsc = min(1.0d0,(hsch**(1.0/par_DSSP2(1)))*hsc)
  else
    wesc = (esch**(1.0/par_DSSP2(1)))*esc
    whsc = (hsch**(1.0/par_DSSP2(1)))*hsc
  end if
!
  do i=1,npepmol
    if (par_DSSP2(2).eq.3) then
      wescm(i) = min(1.0d0,(wescm(i)**(1.0/par_DSSP2(1)))*escm(i))
      whscm(i) = min(1.0d0,(whscm(i)**(1.0/par_DSSP2(1)))*hscm(i))
    else
      wescm(i) = (wescm(i)**(1.0/par_DSSP2(1)))*escm(i)
      whscm(i) = (whscm(i)**(1.0/par_DSSP2(1)))*hscm(i)
    end if
  end do
!
end

!-----------------------------------------------------------------------
!
subroutine update_dssp_der(esc,wesc,hsc,whsc,&
 &                   escm,wescm,hscm,whscm,intraflag)
!
  use iounit
  use sequen
  use dssps
  use atoms
  use polypep
  use mcsums
  use system
  use grandensembles
  use molecule
!
  implicit none
!
  integer i,j,ps1,ps2,normby,pmol,lastp
  RTYPE esc,escm(npepmol),wesc,wescm(npepmol)
  RTYPE hsc,hscm(npepmol),whsc,whscm(npepmol)
  logical ismember,intraflag
!
  hb_mat(:,:) = 0
  hbl(:)%nahbs = 0
  hbl(:)%ndhbs = 0
  dssp_ass(:,:) = 0
!
! now go through the various sub-analyses necessary for DSSP ...
!
  call assign_curves()
!
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    if (intraflag.EQV..true.) then
      normby = pep_sz
      do ps1=1,pep_sz
        pmol = pepmol(molofrs(rspep(ps1)))
        do ps2=ps1+1,pep_sz
          if (molofrs(rspep(ps2)).ne.pmol) exit
          call get_bbhb_der(ps1,ps2)
        end do
        lastp = ps2
        if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
          if (ismember(ispresent,pmol).EQV..false.) cycle
        end if
        do ps2=lastp,pep_sz
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps2)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps2))).EQV..false.) cycle
          end if
          call get_bbhb_der(ps1,ps2)
        end do
      end do
    else
      normby = 0
      do ps1=1,pep_sz
        pmol = pepmol(molofrs(rspep(ps1)))
        if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
          if (ismember(ispresent,pmol).EQV..false.) cycle
        end if
        normby = normby + 1
        do ps2=ps1+1,pep_sz
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps2)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps2))).EQV..false.) cycle
          end if
          call get_bbhb_der(ps1,ps2)
        end do
      end do
    end if
  else
    normby = pep_sz
    do ps1=1,pep_sz
      do ps2=ps1+1,pep_sz
        call get_bbhb_der(ps1,ps2)
      end do
    end do
  end if
!
! now assemble the bridge structure
  bdgtb%nrs = 0
  do i=1,pep_sz
    do j=i+3,pep_sz
      call assign_bridge(i,j)
    end do
  end do
!
  call weigh_bridges_der(esc,wesc,escm,wescm,normby)
!
  call assign_helixturns()
!
  call weigh_helices_der(hsc,whsc,hscm,whscm,normby)
!
end
!
!-----------------------------------------------------------------------
!
subroutine do_dssp()
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use mcsums
  use system
  use molecule
  use grandensembles
  use pdb
!
  implicit none
!
  integer i,j,k,curseg,segl,jj1,jj2,kk1,kk2,pmol
  character dsspstring(8)
  character(3) rnam
  RTYPE en1,en2,esc,wesc,hsc,whsc,hsch,esch,incr
  RTYPE escm(npepmol),wescm(npepmol),hscm(npepmol),whscm(npepmol)
  logical afalse,ismember,endseg,fndit(npepmol)
!
  fndit(:) = .false.
  afalse = .false.
  call update_dssp(esc,wesc,esch,hsc,whsc,hsch,escm,wescm,hscm,whscm,afalse)
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  19  format('##################################  STEP ',i13,'  ',&
 &       '##################################')
  20  format('# NR.     AA   DSSP-ASSIGNMENT   TOP N-H--O-C   TOP C-O--H&
 &-N   SEC N-H--O-C   SEC C-O--H-N')
  21  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &4(i5,f9.2,1x))
  22  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &3(i5,f9.2,1x))
  23  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &2(i5,f9.2,1x),15x,i5,f9.2)
  24  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &2(i5,f9.2,1x))
  26  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &15x,i5,f9.2)
  25  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &i5,f9.2)
  27  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')')
  28  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &i5,f9.2,1x,15x,i5,f9.2)
  29  format(i8,1x,'(',a3,'):',2x,a1,1x,5(a1),1x,'(',i2,')',1x,&
 &15x,i5,f9.2,1x,15x,i5,f9.2)
  30  format('Summary: E-Score:',f8.4,' H-Score: ',f8.4,' E-Fraction: ',&
 &f8.4,' H-Fraction: ',f8.4)
!
  if (inst_dssp.EQV..true.) then
    write(idssp,19) nstep
    write(idssp,20)
    do i=1,pep_sz
      pmol = pepmol(molofrs(rspep(i)))
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
          if (ismember(ispresent,pmol).EQV..false.) cycle
        end if
      end if
      dsspstring(1:8) = ' '
      do j=1,8
        dsspstring(j:j) = dssp_map(j,dssp_ass(i,j)+1)
      end do
!     find lowest-energy H-bonds
      en1 = 0.0
      en2 = 0.0
      jj1 = 0
      jj2 = 0
      do j=1,hbl(i)%nahbs
        if (hbl(i)%enahb(j).lt.en1) then
          jj2 = jj1
          en2 = en1
          jj1 = j
          en1 = hbl(i)%enahb(j)
        else if (hbl(i)%enahb(j).lt.en2) then
          jj2 = j
          en2 = hbl(i)%enahb(j)
        end if
      end do
      en1 = 0.0
      en2 = 0.0
      kk1 = 0
      kk2 = 0
      do k=1,hbl(i)%ndhbs
        if (hbl(i)%endhb(k).lt.en1) then
          kk2 = kk1
          en2 = en1
          kk1 = k
          en1 = hbl(i)%endhb(k)
        else if (hbl(i)%endhb(k).lt.en2) then
          kk2 = k
          en2 = hbl(i)%endhb(k)
        end if
      end do
      rnam = amino(seqtyp(rspep(i)))
      if ((hbl(i)%nahbs.ge.2).AND.(hbl(i)%ndhbs.ge.2)) then
        write(idssp,21) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%dhb(kk1)-rspep(i),&
 &hbl(i)%endhb(kk1),hbl(i)%ahb(jj1)-rspep(i),hbl(i)%enahb(jj1),&
 &hbl(i)%dhb(kk2)-rspep(i),hbl(i)%endhb(kk2),&
 &hbl(i)%ahb(jj2)-rspep(i),hbl(i)%enahb(jj2)
      else if ((hbl(i)%nahbs.ge.1).AND.(hbl(i)%ndhbs.ge.2)) then
        write(idssp,22) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%dhb(kk1)-rspep(i),&
 &hbl(i)%endhb(kk1),hbl(i)%ahb(jj1)-rspep(i),hbl(i)%enahb(jj1),&
 &hbl(i)%dhb(kk2)-rspep(i),hbl(i)%endhb(kk2)
      else if ((hbl(i)%nahbs.ge.2).AND.(hbl(i)%ndhbs.ge.1)) then
        write(idssp,23) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%dhb(kk1)-rspep(i),&
 &hbl(i)%endhb(kk1),hbl(i)%ahb(jj1)-rspep(i),hbl(i)%enahb(jj1),&
 &hbl(i)%ahb(jj2)-rspep(i),hbl(i)%enahb(jj2)
      else if (hbl(i)%ndhbs.ge.2) then
        write(idssp,28) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%dhb(kk1)-rspep(i),&
 &hbl(i)%endhb(kk1),hbl(i)%dhb(kk2)-rspep(i),hbl(i)%endhb(kk2)
      else if (hbl(i)%nahbs.ge.2) then
        write(idssp,29) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%ahb(jj1)-rspep(i),&
 &hbl(i)%enahb(jj1),hbl(i)%ahb(jj2)-rspep(i),hbl(i)%enahb(jj2)
      else if ((hbl(i)%nahbs.ge.1).AND.(hbl(i)%ndhbs.ge.1)) then
        write(idssp,24) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%dhb(kk1)-rspep(i),&
 &hbl(i)%endhb(kk1),hbl(i)%ahb(jj1)-rspep(i),hbl(i)%enahb(jj1)
      else if (hbl(i)%ndhbs.ge.1) then
        write(idssp,25) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%dhb(kk1)-rspep(i),&
 &hbl(i)%endhb(kk1)
      else if (hbl(i)%nahbs.ge.1) then
        write(idssp,26) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1),hbl(i)%ahb(jj1)-rspep(i),&
 &hbl(i)%enahb(jj1)
      else
        write(idssp,27) rspep(i),rnam,dsspstring(1:1),&
 &dsspstring(2:2),dsspstring(3:3),dsspstring(4:4),dsspstring(5:5),&
 &dsspstring(6:6),dssp_ass(i,1)
      end if 
    end do 
  end if
!
! accumulate statistics
  curseg = 0
  segl = 0
  do i=1,pep_sz
!   skip non-presents
!  (note this requires segments to be ended at the last residue not the first one of the next)
    pmol = pepmol(molofrs(rspep(i)))
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
        if (ismember(ispresent,pmol).EQV..false.) cycle
      end if
    end if
!   count
    if (fndit(pmol).EQV..false.) then
      dssp_cnt(pmol) = dssp_cnt(pmol) + 1 
      fndit(pmol) = .true.
    end if
!   last peptide residue in a molecule always ends a segment
    endseg = .false.
    if (i.eq.pep_sz) endseg = .true.
    if (endseg.EQV..false.) then
      if (molofrs(rspep(i)).ne.molofrs(rspep(i+1))) endseg = .true.
    end if
    if (endseg.EQV..true.) then
      if (curseg.ne.0) then
        if (dssp_ass(i,1).eq.curseg) then
          segl = segl + 1
          do j=i-segl,i
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
        else
          do j=i-segl,i-1
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
          if (dssp_ass(i,1).gt.0) then
            curseg = dssp_ass(i,1)
            segl = 1
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end if
        end if
        curseg = 0
        segl = 0
      else
        if (dssp_ass(i,1).gt.0) then
          curseg = dssp_ass(i,1)
          segl = 1
          dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
        end if
        curseg = 0
        segl = 0
      end if
    else
      if (curseg.ne.0) then
        if (dssp_ass(i,1).eq.curseg) then
          segl = segl + 1
        else if (dssp_ass(i,1).ne.0) then
          do j=i-segl,i-1
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
          curseg = dssp_ass(i,1)
          segl = 1
        else
          do j=i-segl,i-1
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
          curseg = 0
          segl = 0
        end if
      else
        if (dssp_ass(i,1).ne.0) then
!         start a new segment
          curseg = dssp_ass(i,1)
          segl = 1
        else
!         do nothing
        end if
      end if          
    end if
  end do
!
! now re-parse for the low resolution H-bonded strand-segment (11)
  curseg = 0
  segl = 0
  do i=1,pep_sz
!   skip non-presents
!  (note this requires segments to be ended at the last residue not the first one of the next)
    pmol = pepmol(molofrs(rspep(i)))
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(pmol)).EQV..true.) then
        if (ismember(ispresent,pmol).EQV..false.) cycle
      end if
    end if
!   do not increment counter in dssp_cnt again (already done above)
!   last peptide residue in a molecule always ends a segment
    endseg = .false.
    if (i.eq.pep_sz) endseg = .true.
    if (endseg.EQV..false.) then
      if (molofrs(rspep(i)).ne.molofrs(rspep(i+1))) endseg = .true.
    end if
    if (endseg.EQV..true.) then
      if (curseg.ne.0) then
        if ((dssp_ass(i,1).ge.1).AND.(dssp_ass(i,1).le.5)) then
          segl = segl + 1
          do j=i-segl,i
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
        else
          do j=i-segl,i-1
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
        end if
        curseg = 0
        segl = 0
      else
        if ((dssp_ass(i,1).ge.1).AND.(dssp_ass(i,1).le.5)) then
          curseg = 11
          segl = 1
          dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
        end if
        curseg = 0
        segl = 0
      end if
    else
      if (curseg.ne.0) then
        if ((dssp_ass(i,1).ge.1).AND.(dssp_ass(i,1).le.5)) then
          segl = segl + 1
        else 
          do j=i-segl,i-1
            dssp_perrs(j,curseg,segl) = dssp_perrs(j,curseg,segl) + incr
          end do
          curseg = 0
          segl = 0
        end if
      else
        if ((dssp_ass(i,1).ge.1).AND.(dssp_ass(i,1).le.5)) then
!         start a new segment
          curseg = 11
          segl = 1
        else
!         do nothing
        end if
      end if          
    end if
  end do
!
! increment 1D-histograms
  jj1 = min(floor(hsc/0.01) + 1,100)
  jj2 = min(floor(esc/0.01) + 1,100)
  dssp_hists(1,jj1) = dssp_hists(1,jj1) + incr
  dssp_hists(2,jj2) = dssp_hists(2,jj2) + incr
  jj1 = min(floor(hsch/par_DSSP(6)) + 1,100)
  jj2 = min(floor(esch/par_DSSP(6)) + 1,100)
  dssp_hists(3,jj1) = dssp_hists(3,jj1) + incr
  dssp_hists(4,jj2) = dssp_hists(4,jj2) + incr
!
! finally increment 2D-histogram
  jj1 = min(floor(whsc/0.01) + 1,100)
  jj2 = min(floor(wesc/0.01) + 1,100)
  dssp_2dhist(jj1,jj2) = dssp_2dhist(jj1,jj2) + incr
!
  if (inst_dssp.EQV..true.) then
    write(idssp,30) wesc,whsc,esc,hsc
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine prt_dssp()
!
  use iounit
  use sequen
  use molecule
  use dssps
  use mpistuff
  use pdb
  use system
!
  implicit none
!
  integer i,j,k,maxi,freeunit,it,it2,ii,jj
  logical exists,fly
  character(60) fn
  RTYPE tot(5)
  RTYPE, ALLOCATABLE:: normer(:),segall(:,:)
#ifdef ENABLE_MPI
  character(re_aux(10)) xpont
#endif
!
 27 format(11i10)
 28 format(11(g14.6,1x))
!
  if (use_frame_weights.EQV..false.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,re_aux(10))
      fn =  'N_'//xpont(1:re_aux(10))//'_DSSP_RES.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'DSSP_RES.dat'
    end if
#else
    fn = 'DSSP_RES.dat'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      it = freeunit()
      open(unit=it,file=fn(ii:jj),status='old')
      close(unit=it,status='delete')
    end if
    it=freeunit()
    open(unit=it,file=fn(ii:jj),status='new')
    maxi = 1
    do i=2,pep_sz
      fly = .false.
      do j=1,pep_sz
        do k=1,11
          if (dssp_perrs(j,k,i).gt.0.0) then
            maxi = i
            fly = .true.
            exit
          end if
        end do
        if (fly.EQV..true.) exit
      end do
    end do
   902  format(i8,100000000(1x,i10))
    do j=1,pep_sz
      write(it,902) rspep(j),((nint(dssp_perrs(j,k,i)), k=1,11), i=1,maxi)
    end do
    close(unit=it)
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,re_aux(10))
    fn =  'N_'//xpont(1:re_aux(10))//'_DSSP_NORM_RES.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'DSSP_NORM_RES.dat'
  end if
#else
  fn = 'DSSP_NORM_RES.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
  maxi = 1
  do i=2,pep_sz
    fly = .false.
    do j=1,pep_sz
      do k=1,11
        if (dssp_perrs(j,k,i).gt.0.0) then
          maxi = i
          fly = .true.
          exit
        end if
      end do
      if (fly.EQV..true.) exit
    end do
  end do
!
 903  format(i8,100000000(1x,g12.5))
  allocate(normer(npepmol))
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),dsspcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)
      end if
    end do
    if (k.ne.dssp_cnt(1)) then
      write(ilog,*) 'Warning. DSSP averages and distributions have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all DSSP-related output (DSSP_RES.dat, DSSP_HIST.dat, etc. ...).'
    end if
  else
    normer(:) = 1.0*dssp_cnt(:)
  end if
  do j=1,pep_sz
    if (normer(pepmol(molofrs(rspep(j)))).le.0.0) cycle
    write(it,903) rspep(j),((dssp_perrs(j,k,i)/normer(pepmol(molofrs(rspep(j)))), k=1,11), i=1,maxi)
  end do
  close(unit=it)
!
! integrate res.-specific information to get per molecule type information
!
  if (use_frame_weights.EQV..false.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,re_aux(10))
      fn =  'N_'//xpont(1:re_aux(10))//'_DSSP.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'DSSP.dat'
    end if
#else
    fn = 'DSSP.dat'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      it = freeunit()
      open(unit=it,file=fn(ii:jj),status='old')
      close(unit=it,status='delete')
    end if
    it=freeunit()
    open(unit=it,file=fn(ii:jj),status='new')
  end if
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,re_aux(10))
    fn =  'N_'//xpont(1:re_aux(10))//'_DSSP_NORM.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'DSSP_NORM.dat'
  end if
#else
  fn = 'DSSP_NORM.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it2 = freeunit()
    open(unit=it2,file=fn(ii:jj),status='old')
    close(unit=it2,status='delete')
  end if
  it2=freeunit()
  open(unit=it2,file=fn(ii:jj),status='new')
!
  allocate(segall(11,maxi))
  segall(:,:) = 0.0 
!
  do k=1,11
    do j=1,maxi
      do i=1,pep_sz
        if (dssp_perrs(i,k,j).gt.0.0) then
          segall(k,j) = segall(k,j) + dssp_perrs(i,k,j)/normer(pepmol(molofrs(rspep(i))))
        end if
      end do
    end do
  end do
  do j=1,maxi
    do k=1,11
      segall(k,j) = segall(k,j)/(1.0*j)
    end do
  end do
  if (use_frame_weights.EQV..false.) then
    do j=1,maxi
      write(it,27) (nint(sum(dssp_perrs(1:pep_sz,k,j))/(1.0*j)),k=1,11)
    end do
  end if
  do j=1,maxi
    write(it2,28) (segall(k,j),k=1,11)
  end do
!
  deallocate(segall)
  if (use_frame_weights.EQV..false.) close(unit=it)
  close(unit=it2)
!
! print out E/H-score histograms
!
 47   format(6(g15.7,1x))
 48   format(100(g14.6,1x))
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,re_aux(10))
    fn =  'N_'//xpont(1:re_aux(10))//'_DSSP_HIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'DSSP_HIST.dat'
  end if
#else
  fn = 'DSSP_HIST.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
!
  do k=1,4
    tot(k) = sum(dssp_hists(k,:))
    if (tot(k).gt.0.0) dssp_hists(k,:) = dssp_hists(k,:)/tot(k)
  end do
!
  if (maxval(tot(1:4)).gt.0.0) then
    do i=1,100
      write(it,47) (i-0.5)*0.01,(i-0.5)*par_DSSP(6),(dssp_hists(k,i),k=1,4)
    end do
  end if
  close(unit=it)
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,re_aux(10))
    fn =  'N_'//xpont(1:re_aux(10))//'_DSSP_EH_HIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'DSSP_EH_HIST.dat'
  end if
#else
  fn = 'DSSP_EH_HIST.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
!
  tot(5) = sum(dssp_2dhist(:,:))
  if (tot(5).gt.0.0) then
    dssp_2dhist(:,:) = dssp_2dhist(:,:)/tot(5)
    do i=1,100
      write(it,48) (dssp_2dhist(k,i),k=1,100)
    end do
  end if
  close(unit=it)
!
  deallocate(normer)
! 
end
!
!-----------------------------------------------------------------------
!
! this routine assembles hydrogen bond lists (in hbl), but also
! a network matrix (hb_mat). for the latter the codes are (always rs1 < rs2):
! 1: rs1 is acceptor only, rs2 donor only
! 2: rs2 is donor only, rs1 acceptor only
! 3: rs1 and rs2 are both both acceptor and donor
!
subroutine get_bbhb_gln(ps1,ps2)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use system
  use sequen
  use dssps
!
  implicit none
!
  integer rs1,rs2,ps1,ps2,rsbu,shf,shf2
  integer nsc1,nsc2,osc1,osc2,hsc1,hsc2,csc1,csc2
  RTYPE don1,doh1,dcn1,dch1,don2,doh2,dcn2,dch2,qfac,d2
  RTYPE hben1,hben2,svec(3),dvec1(3),dvec2(3),dvec3(3),dvec4(3)
  RTYPE dvec5(3),dvec6(3),dvec7(3),dvec8(3),dis,doh21,dch21
  RTYPE doh11,dch11,doh12,dch12,dch22,doh22
!
  shf = 0
  shf2 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 5
  end if
!
  rs1 = rspep(ps1)
  rs2 = rspep(ps2)
!
  nsc1 = at(rs1)%sc(6-shf)
  nsc2 = at(rs2)%sc(6-shf)
  osc1 = at(rs1)%sc(5-shf)
  osc2 = at(rs2)%sc(5-shf)
  csc1 = at(rs1)%sc(4-shf)
  csc2 = at(rs2)%sc(4-shf)
  hsc1 = at(rs1)%sc(11-shf2)
  hsc2 = at(rs2)%sc(11-shf2)
!
  if ((seqpolty(rs1).ne.'P').OR.&
 &    (seqpolty(rs2).ne.'P')) then
!   HB-assignments only work for pairs of polypeptides residues
    return
  end if
!
  qfac = -27.888 ! A*kcal/mol
!
! to be on the safe side ...
  if (rs1.gt.rs2) then
    rsbu = rs1
    rs1 = rs2
    rs2 = rsbu
  end if 
!
  svec(:) = 0.0
  call dis_bound_rs3(rs1,rs2,svec,dis)
! use distance screening to save CPU time
! this of course makes assumptions about the energy cutoff values for the
! H-bond formula below ...
  if (dis.gt.15.0) return
!
  if (abs(rs2-rs1).ge.2) then
! we have to screen for PBC as well as for proline which cannot be a donor
  if (hni(rs2).gt.0) then
    call dis_bound5(oi(rs1),ni(rs2),svec,d2,dvec1) 
    don1 = sqrt(1./d2)
    call dis_bound5(oi(rs1),hni(rs2),svec,d2,dvec2)
    doh1 = sqrt(1./d2)
    call dis_bound5(ci(rs1),hni(rs2),svec,d2,dvec3)
    dch1 = sqrt(1./d2)
    call dis_bound5(ci(rs1),ni(rs2),svec,d2,dvec4)
    dcn1 = sqrt(1./d2)
    hben1 = qfac*(doh1 + dcn1 - don1 - dch1)
  else
    hben1 = 0.0 ! remember par_DSSP(2) is always negative
  end if
!
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = 1
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = 1
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
  if (hni(rs1).gt.0) then
    call dis_bound5(ni(rs1),oi(rs2),svec,d2,dvec5)
    don2 = sqrt(1./d2)
    call dis_bound5(hni(rs1),oi(rs2),svec,d2,dvec6)
    doh2 = sqrt(1./d2)
    call dis_bound5(hni(rs1),ci(rs2),svec,d2,dvec7)
    dch2 = sqrt(1./d2)
    call dis_bound5(ni(rs1),ci(rs2),svec,d2,dvec8)
    dcn2 = sqrt(1./d2)
    hben2 = qfac*(doh2 + dcn2 - don2 - dch2)
  else
    hben2 = 0.0
  end if
!
  if (hben2.lt.par_DSSP(2)) then
!   add donor hydrogen bond for rs1
    if (hbl(peprs(rs1))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%ndhbs = hbl(peprs(rs1))%ndhbs + 1
      hbl(peprs(rs1))%dhb(hbl(peprs(rs1))%ndhbs) = 1
      hbl(peprs(rs1))%endhb(hbl(peprs(rs1))%ndhbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
!   add acceptor hydrogen bond for rs2
    if (hbl(peprs(rs2))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%nahbs = hbl(peprs(rs2))%nahbs + 1
      hbl(peprs(rs2))%ahb(hbl(peprs(rs2))%nahbs) = 1
      hbl(peprs(rs2))%enahb(hbl(peprs(rs2))%nahbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
  end if
!
  if (abs(rs2-rs1).ge.1) then
! sidechain(acceptor)-sidechain(donor: 2x)
  call dis_bound5(osc1,nsc2,svec,d2,dvec1) 
  don1 = sqrt(1./d2)
  call dis_bound5(osc1,hsc2,svec,d2,dvec2)
  doh11 = sqrt(1./d2)
  call dis_bound5(csc1,hsc2,svec,d2,dvec3)
  dch11 = sqrt(1./d2)
  call dis_bound5(osc1,hsc2+1,svec,d2,dvec2)
  doh12 = sqrt(1./d2)
  call dis_bound5(csc1,hsc2+1,svec,d2,dvec3)
  dch12 = sqrt(1./d2)
  call dis_bound5(csc1,nsc2,svec,d2,dvec4)
  dcn1 = sqrt(1./d2)
  hben1 = qfac*(doh11 + dcn1 - don1 - dch11)
!
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = 2
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = 2
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
  hben1 = qfac*(doh12 + dcn1 - don1 - dch12)
!
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = 2
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = 2
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
! sidechain(donor: 2x)-sidechain(acceptor)
  call dis_bound5(nsc1,osc2,svec,d2,dvec5)
  don2 = sqrt(1./d2)
  call dis_bound5(hsc1,osc2,svec,d2,dvec6)
  doh21 = sqrt(1./d2)
  call dis_bound5(hsc1,csc2,svec,d2,dvec7)
  dch21 = sqrt(1./d2)
  call dis_bound5(hsc1+1,osc2,svec,d2,dvec6)
  doh22 = sqrt(1./d2)
  call dis_bound5(hsc1+1,csc2,svec,d2,dvec7)
  dch22 = sqrt(1./d2)
  call dis_bound5(nsc1,csc2,svec,d2,dvec8)
  dcn2 = sqrt(1./d2)
  hben2 = qfac*(doh21 + dcn2 - don2 - dch21)
!
  if (hben2.lt.par_DSSP(2)) then
!   add donor hydrogen bond for rs1
    if (hbl(peprs(rs1))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%ndhbs = hbl(peprs(rs1))%ndhbs + 1
      hbl(peprs(rs1))%dhb(hbl(peprs(rs1))%ndhbs) = 2
      hbl(peprs(rs1))%endhb(hbl(peprs(rs1))%ndhbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
!   add acceptor hydrogen bond for rs2
    if (hbl(peprs(rs2))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%nahbs = hbl(peprs(rs2))%nahbs + 1
      hbl(peprs(rs2))%ahb(hbl(peprs(rs2))%nahbs) = 2
      hbl(peprs(rs2))%enahb(hbl(peprs(rs2))%nahbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
  hben2 = qfac*(doh22 + dcn2 - don2 - dch22)
!
  if (hben2.lt.par_DSSP(2)) then
!   add donor hydrogen bond for rs1
    if (hbl(peprs(rs1))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%ndhbs = hbl(peprs(rs1))%ndhbs + 1
      hbl(peprs(rs1))%dhb(hbl(peprs(rs1))%ndhbs) = 2
      hbl(peprs(rs1))%endhb(hbl(peprs(rs1))%ndhbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
!   add acceptor hydrogen bond for rs2
    if (hbl(peprs(rs2))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%nahbs = hbl(peprs(rs2))%nahbs + 1
      hbl(peprs(rs2))%ahb(hbl(peprs(rs2))%nahbs) = 2
      hbl(peprs(rs2))%enahb(hbl(peprs(rs2))%nahbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
! sidechain(donor: 2x)-backbone(acceptor)
  call dis_bound5(oi(rs1),nsc2,svec,d2,dvec1) 
  don1 = sqrt(1./d2)
  call dis_bound5(oi(rs1),hsc2,svec,d2,dvec2)
  doh11 = sqrt(1./d2)
  call dis_bound5(ci(rs1),hsc2,svec,d2,dvec3)
  dch11 = sqrt(1./d2)
  call dis_bound5(oi(rs1),hsc2+1,svec,d2,dvec2)
  doh12 = sqrt(1./d2)
  call dis_bound5(ci(rs1),hsc2+1,svec,d2,dvec3)
  dch12 = sqrt(1./d2)
  call dis_bound5(ci(rs1),nsc2,svec,d2,dvec4)
  dcn1 = sqrt(1./d2)
  hben1 = qfac*(doh11 + dcn1 - don1 - dch11)
!
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = 3
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = 3
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
  hben1 = qfac*(doh12 + dcn1 - don1 - dch12)
!
  if (hben1.lt.par_DSSP(2)) then
!   add acceptor hydrogen bond for rs1
    if (hbl(peprs(rs1))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%nahbs = hbl(peprs(rs1))%nahbs + 1
      hbl(peprs(rs1))%ahb(hbl(peprs(rs1))%nahbs) = 3
      hbl(peprs(rs1))%enahb(hbl(peprs(rs1))%nahbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
!   add donor hydrogen bond for rs2
    if (hbl(peprs(rs2))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%ndhbs = hbl(peprs(rs2))%ndhbs + 1
      hbl(peprs(rs2))%dhb(hbl(peprs(rs2))%ndhbs) = 3
      hbl(peprs(rs2))%endhb(hbl(peprs(rs2))%ndhbs) = &
 &                                    max(hben1,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
  end if
!
! backbone(donor)-sidechain(acceptor)
  if (hni(rs1).gt.0) then
    call dis_bound5(ni(rs1),osc2,svec,d2,dvec5)
    don2 = sqrt(1./d2)
    call dis_bound5(hni(rs1),osc2,svec,d2,dvec6)
    doh2 = sqrt(1./d2)
    call dis_bound5(hni(rs1),csc2,svec,d2,dvec7)
    dch2 = sqrt(1./d2)
    call dis_bound5(ni(rs1),csc2,svec,d2,dvec8)
    dcn2 = sqrt(1./d2)
    hben2 = qfac*(doh2 + dcn2 - don2 - dch2)
  else
    hben2 = 0.0
  end if
!
  if (hben2.lt.par_DSSP(2)) then
!   add donor hydrogen bond for rs1
    if (hbl(peprs(rs1))%ndhbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around donor &
 &of residue ',rs1,' exceeds ',MAXHBS,'. This is suggestive of bad p&
 &arameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs1))%ndhbs = hbl(peprs(rs1))%ndhbs + 1
      hbl(peprs(rs1))%dhb(hbl(peprs(rs1))%ndhbs) = 4
      hbl(peprs(rs1))%endhb(hbl(peprs(rs1))%ndhbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
!   add acceptor hydrogen bond for rs2
    if (hbl(peprs(rs2))%nahbs.eq.MAXHBS) then
      write(ilog,*) 'Warning. Number of hydrogen bonds around accept&
 &or of residue ',rs2,' exceeds ',MAXHBS,'. This is suggestive of ba&
 &d parameters or a bug (in get_bbhb(...)).'
    else
      hbl(peprs(rs2))%nahbs = hbl(peprs(rs2))%nahbs + 1
      hbl(peprs(rs2))%ahb(hbl(peprs(rs2))%nahbs) = 4
      hbl(peprs(rs2))%enahb(hbl(peprs(rs2))%nahbs) = &
 &                                    max(hben2,par_DSSP(3))
    end if
  else
!   do nothing: hydrogen bond criteria not met
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine hbs_intra(counts)
!
  use iounit
  use polypep
  use atoms
  use aminos
  use sequen
  use dssps
  use mcsums
  use system
  use grandensembles
  use molecule
!
  implicit none
!
  integer i,k
  integer ps1,ps2
  logical ismember
  RTYPE counts(16)
!
  hbl(:)%nahbs = 0
  hbl(:)%ndhbs = 0
!
  do ps1=1,pep_sz
    do ps2=ps1,pep_sz
      if (molofrs(rspep(ps1)).eq.molofrs(rspep(ps2))) then
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps1)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps1))).EQV..false.) cycle
          end if
        end if
        call get_bbhb_gln(ps1,ps2)
      end if
    end do
  end do
!
  do i=1,pep_sz
    do k=1,hbl(i)%nahbs
      if (hbl(i)%ahb(k).eq.1) counts(1) = counts(1) + 1.0/pep_sz
      if (hbl(i)%ahb(k).eq.2) counts(3) = counts(3) + 1.0/pep_sz
      if (hbl(i)%ahb(k).eq.3) counts(5) = counts(5) + 1.0/pep_sz
      if (hbl(i)%ahb(k).eq.4) counts(7) = counts(7) + 1.0/pep_sz
    end do
    do k=1,hbl(i)%ndhbs
      if (hbl(i)%dhb(k).eq.1) counts(2) = counts(2) + 1.0/pep_sz
      if (hbl(i)%dhb(k).eq.2) counts(4) = counts(4) + 1.0/pep_sz
      if (hbl(i)%dhb(k).eq.3) counts(6) = counts(6) + 1.0/pep_sz
      if (hbl(i)%dhb(k).eq.4) counts(8) = counts(8) + 1.0/pep_sz
    end do
  end do
!  counts(9) = counts(9) + sum(hbl(1:pep_sz)%nahbs)/(2.0*pep_sz)
!  counts(10) = counts(10) + sum(hbl(1:pep_sz)%ndhbs)/(3.0*pep_sz)
!
  hbl(:)%nahbs = 0
  hbl(:)%ndhbs = 0
!
  do ps1=1,pep_sz
    do ps2=ps1,pep_sz
      if (molofrs(rspep(ps1)).ne.molofrs(rspep(ps2))) then
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps1)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps1))).EQV..false.) cycle
          end if
          if (ismember(fluctypes,moltypid(molofrs(rspep(ps2)))).EQV..true.) then
            if (ismember(ispresent,molofrs(rspep(ps2))).EQV..false.) cycle
          end if
        end if
        call get_bbhb_gln(ps1,ps2)
      end if
    end do
  end do
!
  do i=1,pep_sz
    do k=1,hbl(i)%nahbs
      if (hbl(i)%ahb(k).eq.1) counts(9) = counts(9) + 1.0/pep_sz
      if (hbl(i)%ahb(k).eq.2) counts(11) = counts(11) + 1.0/pep_sz
      if (hbl(i)%ahb(k).eq.3) counts(13) = counts(13) + 1.0/pep_sz
      if (hbl(i)%ahb(k).eq.4) counts(15) = counts(15) + 1.0/pep_sz
    end do
    do k=1,hbl(i)%ndhbs
      if (hbl(i)%dhb(k).eq.1) counts(10) = counts(10) + 1.0/pep_sz
      if (hbl(i)%dhb(k).eq.2) counts(12) = counts(12) + 1.0/pep_sz
      if (hbl(i)%dhb(k).eq.3) counts(14) = counts(14) + 1.0/pep_sz
      if (hbl(i)%dhb(k).eq.4) counts(16) = counts(16) + 1.0/pep_sz
    end do
  end do
!
end
