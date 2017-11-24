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
#define VCI_SIZE 92
#define VCI_SIZEM1 91
#define VSN_RESIZE 20
!
!
!------------------------------------------------------------------------------
!
! the routine to compute energies and forces for a pair of residues covering terms
! IPP, ATTLJ, WCA, and the SAV increments and NB list generation
! this should not be called with pairs of crosslinked residues
!
! this routine is not thread-safe per se as it operates on globals, specifically
! svte
!
! the analogous routine, Vforce_rsp2 below, uses only passed arrays
!
subroutine Vforce_rsp(evec,rs1,rs2,cut,ca_f,si)
!
  use iounit
  use inter
  use atoms
  use params
  use energies
  use sequen
  use polypep
  use math
  use cutoffs
  use grandensembles
  use system
  use forces
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  integer rs,i,j,k,ii,kk,imol1,imol2,alcsz
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE idvn1(at(rs1)%na*at(rs2)%na,3)
  RTYPE tmaxd(at(rs1)%na*at(rs2)%na)
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE evec(MAXENERGYTERMS),incr
  RTYPE svec(3),efvoli,efvolk,datri,datrk,pfac,rt23,pwc2,pwc3,derinc,olaps(3),derincs(3)
  logical ismember,doid2
  RTYPE for_tmp(at(rs1)%na*at(rs2)%na,3)
  RTYPE ca_f(n,3)
  type(t_sisa) si(n)
!
! no intialization needed: note that on the calling side it requires a lot of
! care to handle fxns which "silently" increment arguments
!
  if (ideal_run.EQV..true.) then
    return
  end if
!  
! Prevent non-present molecules from interacting with each other
! or with present molecules.  Note that intra-nonpresent-molecule
! interactions are still allowed, in order to correctly simulate
! an infinite-dilution implicit solvent reference state particle
! bath -> partially redundant overlap with use_FEG
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    imol1 = molofrs(rs1)
    imol2 = molofrs(rs2)
    if ((imol1.ne.imol2).AND.((.not.ismember(ispresent,imol1)).OR. &
 &   (.not.ismember(ispresent,imol2)))) return
  end if
!
! in FEG-1 we branch out if precisely one of the residues is ghosting
! always remember that the interaction of two ghost-residues (incl self!)
! is the full Hamiltonian, not the ghosted one
! in FEG-2 we branch out if either residue is ghosted
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
   if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Vforce_rsp_feg(evec,rs1,rs2,cut,ca_f)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Vforce_rsp_feg(evec,rs1,rs2,cut,ca_f)
        return
      end if
    end if
  end if
!
  doid2 = .false.
  if  ((use_attLJ.EQV..true.).OR.&
 &     ((use_IPP.EQV..true.).AND.(use_hardsphere.EQV..false.)).OR.&
 &      (use_WCA.EQV..true.)) then
    doid2 = .true.
  end if
  for_tmp(:,:) = 0.0
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
!   intra-residue interactions are never subject to PBC (residue size should be within ~10A)
    rs = rs1
!   first set the necessary atom-specific parameters (mimimal)
    k = 1
    do i=1,nrsintra(rs)
      ii = iaa(rs)%atin(i,1)
      kk = iaa(rs)%atin(i,2)
      dvec(k,1) = x(kk) - x(ii)
      dvec(k,2) = y(kk) - y(ii)
      dvec(k,3) = z(kk) - z(ii)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsin(i)*fudge(rs)%rsin_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
          terms(k) = fudge(rs)%rsin_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
          idx2(k,1) = ii
          idx2(k,2) = kk
          k = k + 1
        end if
      end if
    end do
  else if (abs(rs1-rs2).eq.1) then
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC 
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
    k = 1
    do i=1,nrsnb(rs)
      ii = iaa(rs)%atnb(i,1)
      kk = iaa(rs)%atnb(i,2)
      dvec(k,1) = x(kk) - x(ii) + svec(1)
      dvec(k,2) = y(kk) - y(ii) + svec(2)
      dvec(k,3) = z(kk) - z(ii) + svec(3)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsnb(i)*fudge(rs)%rsnb_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
          terms(k) = fudge(rs)%rsnb_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
          idx2(k,1) = ii
          idx2(k,2) = kk
          k = k + 1
        end if
      end if
    end do
  else
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!   loop over complete set of atom-atom interactions in distant residues
    k = 1
    do i=1,at(rs1)%na
      if (i.le.at(rs1)%nbb) then
        ii = at(rs1)%bb(i)
      else
        ii = at(rs1)%sc(i-at(rs1)%nbb)
      end if
      do j=1,at(rs2)%na
        if (j.le.at(rs2)%nbb) then
          kk = at(rs2)%bb(j)
        else
          kk = at(rs2)%sc(j-at(rs2)%nbb)
        end if
        dvec(k,1) = x(kk) - x(ii) + svec(1)
        dvec(k,2) = y(kk) - y(ii) + svec(2)
        dvec(k,3) = z(kk) - z(ii) + svec(3)
        d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
          term0(k) = lj_eps(attyp(ii),attyp(kk))
          if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
            terms(k) = lj_sig(attyp(ii),attyp(kk))
            if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
    end do
  end if
  k = k - 1
  alcsz = k
  if (alcsz.gt.0) then
    if (doid2.EQV..true.) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
    end if
    if (use_IPP.EQV..true.) then
      term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf
      term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
      evec(1) = evec(1) + scale_IPP*4.0*&
 &                sum(term2(1:alcsz))
      pfac = scale_IPP*4.0*nindx
      term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
      for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term3(1:alcsz)*dvec(1:alcsz,1)
      for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term3(1:alcsz)*dvec(1:alcsz,2)
      for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term3(1:alcsz)*dvec(1:alcsz,3)
    end if
    if (use_attLJ.EQV..true.) then
      term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3
      term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
      evec(3) = evec(3) - scale_attLJ*4.0*&
 &                sum(term2(1:alcsz))
      pfac = scale_attLJ*24.0
      term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
      for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term3(1:alcsz)*dvec(1:alcsz,1)
      for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term3(1:alcsz)*dvec(1:alcsz,2)
      for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term3(1:alcsz)*dvec(1:alcsz,3)
    end if
    if (use_IMPSOLV.EQV..true.) then
      d1(1:alcsz) = sqrt(d2(1:alcsz))
      idvn1(1:alcsz,1) = dvec(1:alcsz,1)/d1(1:alcsz)
      idvn1(1:alcsz,2) = dvec(1:alcsz,2)/d1(1:alcsz)
      idvn1(1:alcsz,3) = dvec(1:alcsz,3)/d1(1:alcsz)
      tmaxd(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1)
      if (par_IMPSOLV2(1).eq.1) then
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          if (d1(k).lt.tmaxd(k)) then
            efvoli = atsavred(ii)*atvol(ii)
            efvolk = atsavred(kk)*atvol(kk)
            datri = 2.0*atr(ii)
            datrk = 2.0*atr(kk)
            if (d1(k).gt.(tmaxd(k)-datrk)) then
              incr = -(tmaxd(k)-d1(k))/datrk
              derinc = -efvolk/datrk
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = efvolk/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
            svte(ii) = svte(ii) + incr*efvolk
            if (derinc.ne.0.0) then
              if ((si(ii)%nix+1).gt.si(ii)%alsz) call Vsisa_resize(ii,si,VSN_RESIZE)
              si(ii)%nix = si(ii)%nix + 1
              si(ii)%ix(si(ii)%nix) = kk
              si(ii)%dr(1,si(ii)%nix) = idvn1(k,1)*derinc
              si(ii)%dr(2,si(ii)%nix) = idvn1(k,2)*derinc
              si(ii)%dr(3,si(ii)%nix) = idvn1(k,3)*derinc
            end if
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
              derinc = efvoli/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
              derinc = -efvoli/termr(k)
            else
              incr = -1.0
              derinc = 0.0
            end if
            svte(kk) = svte(kk) + incr*efvoli
            if (derinc.ne.0.0) then
              if ((si(kk)%nix+1).gt.si(kk)%alsz) call Vsisa_resize(kk,si,VSN_RESIZE)
              si(kk)%nix = si(kk)%nix + 1
              si(kk)%ix(si(kk)%nix) = ii
              si(kk)%dr(1,si(kk)%nix) = idvn1(k,1)*derinc
              si(kk)%dr(2,si(kk)%nix) = idvn1(k,2)*derinc
              si(kk)%dr(3,si(kk)%nix) = idvn1(k,3)*derinc
            end if
          end if
        end do
      else if (par_IMPSOLV2(1).eq.2) then
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          if (d1(k).lt.tmaxd(k)) then
            call sphere_overlaps_wder(atr(ii),atr(kk),d1(k),par_IMPSOLV(1),atvol(ii),atvol(kk),olaps,derincs)
            if (olaps(2).gt.0.0) then
              svte(ii) = svte(ii) - (olaps(2)-olaps(1))*atsavred(kk)
              if ((derincs(2).ne.0.0).OR.(derincs(1).ne.0.0)) then
                if ((si(ii)%nix+1).gt.si(ii)%alsz) call Vsisa_resize(ii,si,VSN_RESIZE)
                si(ii)%nix = si(ii)%nix + 1
                si(ii)%ix(si(ii)%nix) = kk
                derinc = (derincs(2)-derincs(1))*atsavred(kk)
                si(ii)%dr(1,si(ii)%nix) = idvn1(k,1)*derinc
                si(ii)%dr(2,si(ii)%nix) = idvn1(k,2)*derinc
                si(ii)%dr(3,si(ii)%nix) = idvn1(k,3)*derinc
              end if
            end if
            if (olaps(3).gt.0.0) then
              svte(kk) = svte(kk) - (olaps(3)-olaps(1))*atsavred(ii)
              if ((derincs(3).ne.0.0).OR.(derincs(1).ne.0.0)) then
                if ((si(kk)%nix+1).gt.si(kk)%alsz) call Vsisa_resize(kk,si,VSN_RESIZE)
                si(kk)%nix = si(kk)%nix + 1
                si(kk)%ix(si(kk)%nix) = ii
                derinc = (derincs(1)-derincs(3))*atsavred(ii)
                si(kk)%dr(1,si(kk)%nix) = idvn1(k,1)*derinc
                si(kk)%dr(2,si(kk)%nix) = idvn1(k,2)*derinc
                si(kk)%dr(3,si(kk)%nix) = idvn1(k,3)*derinc
              end if
            end if
          end if
        end do
      end if
    end if
    if (use_WCA.EQV..true.) then
      pfac = scale_WCA*par_WCA(2)
      rt23 = ROOT26*ROOT26
      pwc2 = par_WCA(1)*par_WCA(1)
      pwc3 = 24.0*scale_WCA
      do k=1,alcsz
        if (term0(k).gt.0.0) then
          if (d2(k).lt.rt23*terms(k)) then
            term1(k) = (terms(k)*id2(k))**3
            term2(k) = term1(k)**2
            term3(k) = 4.0*(term2(k) - term1(k)) + 1.0
            evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                             (term3(k) - par_WCA(2))
            term3(k) = -pwc3*(2.0*term2(k) - term1(k))*term0(k)*id2(k)
            for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
            for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
            for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
          else if (d2(k).lt.pwc2*terms(k)) then
            term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
            term2(k) = 0.5*cos(term1(k)) - 0.5
            evec(5) = evec(5) + pfac*term0(k)*term2(k)
            term3(k) = -sin(term1(k))*pfac*term0(k)*par_WCA(3)/terms(k)
            for_tmp(k,1) = for_tmp(k,1) + term3(k)*dvec(k,1)
            for_tmp(k,2) = for_tmp(k,2) + term3(k)*dvec(k,2)
            for_tmp(k,3) = for_tmp(k,3) + term3(k)*dvec(k,3)
          end if
        end if
      end do
    end if
    do j=1,3
      do k=1,alcsz
        ii = idx2(k,1)
        kk = idx2(k,2)
        ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
        ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
      end do
    end do
  end if
!
end
!
!--------------------------------------------------------------------------
!
! same as Vforce_rsp when ghosting is enabled (called from within Vforce_rsp)
!
subroutine Vforce_rsp_feg(evec,rs1,rs2,cut,ca_f)
!
  use iounit
  use inter
  use atoms
  use params
  use energies
  use sequen
  use polypep
  use math
  use cutoffs
  use grandensembles
  use system
  use forces
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  integer rs,i,j,k,ii,kk,alcsz
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE term4(at(rs1)%na*at(rs2)%na)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac
  logical doid2
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE for_tmp(at(rs1)%na*at(rs2)%na,3)
  RTYPE ca_f(n,3)

! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
!
  doid2 = .false.
  if  ((use_FEGS(3).EQV..true.).OR.(use_FEGS(1).EQV..true.)) then
    doid2 = .true.
  else
!   currently (this might change), only these two terms supported
    return
  end if
  for_tmp(:,:) = 0.0
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
  if (rs1.eq.rs2) then
!   intra-residue interactions are never subject to PBC (residue size should be within ~10A)
    rs = rs1
!   first set the necessary atom-specific parameters (mimimal)
    k = 1
    do i=1,nrsintra(rs)
      ii = iaa(rs)%atin(i,1)
      kk = iaa(rs)%atin(i,2)
      dvec(k,1) = x(kk) - x(ii)
      dvec(k,2) = y(kk) - y(ii)
      dvec(k,3) = z(kk) - z(ii)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsin(i)*fudge(rs)%rsin_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if (term0(k).gt.0.0) then
          terms(k) = fudge(rs)%rsin_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          idx2(k,1) = ii
          idx2(k,2) = kk
          k = k + 1
        end if
      end if
    end do
!   correction terms must not be encountered here
    if (use_CORR.EQV..true.) then
      write(ilog,*) 'Fatal. This is a bug. In the FEG-mode 1, intram&
 &olecular contributions are to be omitted, and in FEG-mode 2 they o&
 &ught to be explicitly excluded. Still encountered computation of t&
 &orsional correction terms. Please report this problem.'
      call fexit()
    end if
  else if (abs(rs1-rs2).eq.1) then
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC 
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
    k = 1
    do i=1,nrsnb(rs)
      ii = iaa(rs)%atnb(i,1)
      kk = iaa(rs)%atnb(i,2)
      dvec(k,1) = x(kk) - x(ii) + svec(1)
      dvec(k,2) = y(kk) - y(ii) + svec(2)
      dvec(k,3) = z(kk) - z(ii) + svec(3)
      d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
      if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
        term0(k) = fudge(rs)%rsnb(i)*fudge(rs)%rsnb_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if (term0(k).gt.0.0) then
          terms(k) = fudge(rs)%rsnb_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          idx2(k,1) = ii
          idx2(k,2) = kk
          k = k + 1
        end if
      end if
    end do
  else
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!   loop over complete set of atom-atom interactions in distant residues
    k = 1
    do i=1,at(rs1)%na
      if (i.le.at(rs1)%nbb) then
        ii = at(rs1)%bb(i)
      else
        ii = at(rs1)%sc(i-at(rs1)%nbb)
      end if
      do j=1,at(rs2)%na
        if (j.le.at(rs2)%nbb) then
          kk = at(rs2)%bb(j)
        else
          kk = at(rs2)%sc(j-at(rs2)%nbb)
        end if
        dvec(k,1) = x(kk) - x(ii) + svec(1)
        dvec(k,2) = y(kk) - y(ii) + svec(2)
        dvec(k,3) = z(kk) - z(ii) + svec(3)
        d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.&
 &      (cut.EQV..false.)) then
          term0(k) = lj_eps(attyp(ii),attyp(kk))
          if (term0(k).gt.0.0) then
            terms(k) = lj_sig(attyp(ii),attyp(kk))
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
    end do
  end if
!
  k = k - 1
  alcsz = k
  if (alcsz.gt.0) then
    if (doid2.EQV..true.) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
    end if
    if (doid2.EQV..true.) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3
    end if
    if (use_FEGS(1).EQV..true.) then
      term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
      term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz)
      evec(1) = evec(1) + par_FEG2(1)*4.0*sum(term3(1:alcsz))
      pfac = par_FEG2(1)*48.0
      term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
      for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) - term4(1:alcsz)*dvec(1:alcsz,1)
      for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) - term4(1:alcsz)*dvec(1:alcsz,2)
      for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) - term4(1:alcsz)*dvec(1:alcsz,3)
    end if
    if (use_FEGS(3).EQV..true.) then
      term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
      term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)
      evec(3) = evec(3) - par_FEG2(5)*4.0*sum(term3(1:alcsz))
      pfac = par_FEG2(5)*24.0
      term4(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
      for_tmp(1:alcsz,1) = for_tmp(1:alcsz,1) + term4(1:alcsz)*dvec(1:alcsz,1)
      for_tmp(1:alcsz,2) = for_tmp(1:alcsz,2) + term4(1:alcsz)*dvec(1:alcsz,2)
      for_tmp(1:alcsz,3) = for_tmp(1:alcsz,3) + term4(1:alcsz)*dvec(1:alcsz,3)
    end if
    do j=1,3
      do k=1,alcsz
        ii = idx2(k,1)
        kk = idx2(k,2)
        ca_f(ii,j) = ca_f(ii,j) + for_tmp(k,j)
        ca_f(kk,j) = ca_f(kk,j) - for_tmp(k,j)
      end do
    end do
  end if
!
end
!
!-----------------------------------------------------------------------------------
!
! the generally superior routine to operate on the transpose array covering the cases of crosslinks, GC presence, FEG,
! and neighbor exclusions directly 
!
subroutine Vforce_rsp2(evec,rs1,rs2,cut,ca_f,s_sv,si)
!
  use molecule
  use iounit
  use energies
  use atoms
  use grandensembles
  use system
  use sequen
  use polypep
  use cutoffs
  use params
  use forces, ONLY: t_sisa
  use math
  use fyoc, ONLY: disulf
  use inter
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  integer imol1,imol2,i,j,k,ii,kk,alcsz,lk,ll,rlo
  logical isfeg,ismember,doid2
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n),s_sv(n),efvolk,efvoli,datri,datrk,pfac,rt23,pwc2,pwc3
  type(t_sisa) si(n)
!
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE svec(3),tvec(3),incr,derinc,olaps(3),derincs(3)
  RTYPE dvec(3,at(rs1)%na*at(rs2)%na),for_tmp(3,at(rs1)%na*at(rs2)%na)
  RTYPE d1(at(rs1)%na*at(rs2)%na),d2(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na),terms(at(rs1)%na*at(rs2)%na),termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na),term2(at(rs1)%na*at(rs2)%na),term3(at(rs1)%na*at(rs2)%na)

!
  imol1 = molofrs(rs1)
  imol2 = molofrs(rs2)
!
  if (ideal_run.EQV..true.) return
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    if ((imol1.ne.imol2).AND.((.not.ismember(ispresent,imol1)).OR.(.not.ismember(ispresent,imol2)))) return
  end if
!
  isfeg = .false.
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
      if (((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &        ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) isfeg = .true.
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.)) isfeg = .true.
    end if
    if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
  end if
!
  call dis_bound_rs(rs1,rs2,svec)
!
  k = 1
!
  if (isfeg.EQV..false.) then
!
    doid2 = .false.
    if ((use_attLJ.EQV..true.).OR.((use_IPP.EQV..true.).AND.(use_hardsphere.EQV..false.)).OR.(use_WCA.EQV..true.)) then
      doid2 = .true.
    end if
!
    if (max(rs2,rs1).eq.disulf(min(rs1,rs2))) then
      lk = crlk_idx(rs1)
      ll = 1 ! we can step through the crosslink excludes because they are sorted exactly the same
      do i=1,at(rs1)%na
        if (i.le.at(rs1)%nbb) then
          ii = at(rs1)%bb(i)
        else
          ii = at(rs1)%sc(i-at(rs1)%nbb)
        end if
        do j=1,at(rs2)%na
          if (j.le.at(rs2)%nbb) then
            kk = at(rs2)%bb(j)
          else
            kk = at(rs2)%sc(j-at(rs2)%nbb)
          end if
          if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) then
            if (crosslink(lk)%is14in(ll).EQV..false.) then
              ll = ll + 1
              cycle
            end if
          end if
          dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii) + svec(:)
          d2(k) = sum(dvec(:,k)**2) 
          if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
!           we must not have any special exclusion rules for crosslinks
            if (ll.le.crosslink(lk)%nrsin) then
              if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) then ! must be a 1-4
                term0(k) = fudge_st_14*lj_eps_14(attyp(ii),attyp(kk))
                terms(k) = lj_sig_14(attyp(ii),attyp(kk))
                ll = ll + 1
              else
                term0(k) = lj_eps(attyp(ii),attyp(kk))
                terms(k) = lj_sig(attyp(ii),attyp(kk))
              end if
            else
              term0(k) = lj_eps(attyp(ii),attyp(kk))
              terms(k) = lj_sig(attyp(ii),attyp(kk))
            end if
            if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
              if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
              idx2(k,1) = ii
              idx2(k,2) = kk
              k = k + 1
            end if
          else if (ll.le.crosslink(lk)%nrsin) then
            if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) ll = ll + 1
          end if
        end do
      end do
!
    else if (rs1.eq.rs2) then
!
      do i=1,nrsintra(rs1)
        ii = iaa(rs1)%atin(i,1)
        kk = iaa(rs1)%atin(i,2)
        dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii)
        d2(k) = sum(dvec(:,k)**2) 
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
          term0(k) = fudge(rs1)%rsin(i)*fudge(rs1)%rsin_lje(i)
          if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
            terms(k) = fudge(rs1)%rsin_ljs(i)
            if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
!
    else if (abs(rs1-rs2).eq.1) then
!
      rlo = min(rs1,rs2)
      if (rlo.eq.rs2) svec(:) = -svec(:)
      do i=1,nrsnb(rlo)
        ii = iaa(rlo)%atnb(i,1)
        kk = iaa(rlo)%atnb(i,2)
        dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii) + svec(:)
        d2(k) = sum(dvec(:,k)**2) 
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
          term0(k) = fudge(rlo)%rsnb(i)*fudge(rlo)%rsnb_lje(i)
          if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
            terms(k) = fudge(rlo)%rsnb_ljs(i)
            if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
!
    else
!
      do i=1,at(rs1)%na
        if (i.le.at(rs1)%nbb) then
          ii = at(rs1)%bb(i)
        else
          ii = at(rs1)%sc(i-at(rs1)%nbb)
        end if
        do j=1,at(rs2)%na
          if (j.le.at(rs2)%nbb) then
            kk = at(rs2)%bb(j)
          else
            kk = at(rs2)%sc(j-at(rs2)%nbb)
          end if
          dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii) + svec(:)
          d2(k) = sum(dvec(:,k)**2) 
          if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
            term0(k) = lj_eps(attyp(ii),attyp(kk))
            if ((term0(k).gt.0.0).OR.(use_IMPSOLV.EQV..true.)) then
              terms(k) = lj_sig(attyp(ii),attyp(kk))
              if (use_IMPSOLV.EQV..true.) termr(k) = atr(ii)+atr(kk)
              idx2(k,1) = ii
              idx2(k,2) = kk
              k = k + 1
            end if
          end if
        end do
      end do
    end if
!
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      for_tmp(:,1:alcsz) = 0.0
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
      end if
      if (use_IPP.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(1) = evec(1) + scale_IPP*4.0*sum(term2(1:alcsz))
        pfac = scale_IPP*4.0*nindx
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        do k=1,alcsz
          for_tmp(:,k) = for_tmp(:,k) - term3(k)*dvec(:,k)
        end do
      end if
      if (use_attLJ.EQV..true.) then
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3
        term2(1:alcsz) = term0(1:alcsz)*term1(1:alcsz)
        evec(3) = evec(3) - scale_attLJ*4.0*sum(term2(1:alcsz))
        pfac = scale_attLJ*24.0
        term3(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*pfac
        do k=1,alcsz
          for_tmp(:,k) = for_tmp(:,k) + term3(k)*dvec(:,k)
        end do
      end if
      if (use_IMPSOLV.EQV..true.) then
        term2(1:alcsz) = termr(1:alcsz) + par_IMPSOLV(1)
        term3(1:alcsz) = term2(1:alcsz)*term2(1:alcsz)
        if (par_IMPSOLV2(1).eq.1) then
          do k=1,alcsz
            ii = idx2(k,1)
            kk = idx2(k,2)
            if (d2(k).lt.term3(k)) then
              d1(k) = sqrt(d2(k))
              efvoli = atsavred(ii)*atvol(ii)
              efvolk = atsavred(kk)*atvol(kk)
              datri = 2.0*atr(ii)
              datrk = 2.0*atr(kk)
              tvec(:) = dvec(:,k)/d1(k)
              if (d1(k).gt.(term2(k)-datrk)) then
                incr = -(term2(k)-d1(k))/datrk
                derinc = -efvolk/datrk
              else if (d1(k).lt.termr(k)) then
                incr = -d1(k)/termr(k)
                derinc = efvolk/termr(k)
              else
                incr = -1.0
                derinc = 0.0
              end if
              s_sv(ii) = s_sv(ii) + incr*efvolk
              if (derinc.ne.0.0) then
                if ((si(ii)%nix+1).gt.si(ii)%alsz) call Vsisa_resize(ii,si,VSN_RESIZE)
                si(ii)%nix = si(ii)%nix + 1
                si(ii)%ix(si(ii)%nix) = kk
                si(ii)%dr(:,si(ii)%nix) = derinc*tvec(:)
              end if
              if (d1(k).gt.(term2(k)-datri)) then
                incr = -(term2(k)-d1(k))/datri
                derinc = efvoli/datri
              else if (d1(k).lt.termr(k)) then
                incr = -d1(k)/termr(k)
                derinc = -efvoli/termr(k)
              else
                incr = -1.0
                derinc = 0.0
              end if
              s_sv(kk) = s_sv(kk) + incr*efvoli
              if (derinc.ne.0.0) then
                if ((si(kk)%nix+1).gt.si(kk)%alsz) call Vsisa_resize(kk,si,VSN_RESIZE)
                si(kk)%nix = si(kk)%nix + 1
                si(kk)%ix(si(kk)%nix) = ii
                si(kk)%dr(:,si(kk)%nix) = derinc*tvec(:)
              end if
            end if
          end do
        else if (par_IMPSOLV2(1).eq.2) then
          do k=1,alcsz
            ii = idx2(k,1)
            kk = idx2(k,2)
            if (d2(k).lt.term3(k)) then
              d1(k) = sqrt(d2(k))
              tvec(:) = dvec(:,k)/d1(k)
              call sphere_overlaps_wder(atr(ii),atr(kk),d1(k),par_IMPSOLV(1),atvol(ii),atvol(kk),olaps,derincs)
              if (olaps(2).gt.0.0) then
                s_sv(ii) = s_sv(ii) - (olaps(2)-olaps(1))*atsavred(kk)
                if ((derincs(2).ne.0.0).OR.(derincs(1).ne.0.0)) then
                  if ((si(ii)%nix+1).gt.si(ii)%alsz) call Vsisa_resize(ii,si,VSN_RESIZE)
                  si(ii)%nix = si(ii)%nix + 1
                  si(ii)%ix(si(ii)%nix) = kk
                  derinc = (derincs(2)-derincs(1))*atsavred(kk)
                  si(ii)%dr(:,si(ii)%nix) = derinc*tvec(:)
                end if
              end if
              if (olaps(3).gt.0.0) then
                s_sv(kk) = s_sv(kk) - (olaps(3)-olaps(1))*atsavred(ii)
                if ((derincs(3).ne.0.0).OR.(derincs(1).ne.0.0)) then
                  if ((si(kk)%nix+1).gt.si(kk)%alsz) call Vsisa_resize(kk,si,VSN_RESIZE)
                  si(kk)%nix = si(kk)%nix + 1
                  si(kk)%ix(si(kk)%nix) = ii
                  derinc = (derincs(1)-derincs(3))*atsavred(ii)
                  si(kk)%dr(:,si(kk)%nix) = derinc*tvec(:)
                end if
              end if
            end if
          end do
        end if
      end if
      if (use_WCA.EQV..true.) then
        pfac = scale_WCA*par_WCA(2)
        rt23 = ROOT26*ROOT26
        pwc2 = par_WCA(1)*par_WCA(1)
        pwc3 = 24.0*scale_WCA
        do k=1,alcsz
          ii = idx2(k,1)
          kk = idx2(k,2)
          if (term0(k).gt.0.0) then
            if (d2(k).lt.rt23*terms(k)) then
              term1(k) = (terms(k)*id2(k))**3
              term2(k) = term1(k)**2
              term3(k) = 4.0*(term2(k) - term1(k)) + 1.0
              evec(5) = evec(5) + scale_WCA*term0(k)*(term3(k) - par_WCA(2))
              term3(k) = -pwc3*(2.0*term2(k) - term1(k))*term0(k)*id2(k)
              for_tmp(:,k) = for_tmp(:,k) + term3(k)*dvec(:,k)
            else if (d2(k).lt.pwc2*terms(k)) then
              term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
              term2(k) = 0.5*cos(term1(k)) - 0.5
              evec(5) = evec(5) + pfac*term0(k)*term2(k)
              term3(k) = -sin(term1(k))*pfac*term0(k)*par_WCA(3)/terms(k)
              for_tmp(:,k) = for_tmp(:,k) + term3(k)*dvec(:,k)
            end if
          end if
        end do
      end if
    end if
!
  else
!
    doid2 = .false.
    if  ((use_FEGS(3).EQV..true.).OR.(use_FEGS(1).EQV..true.)) then
      doid2 = .true.
    else
!     currently (this might change), only these two terms supported
      return
    end if
!
    if (max(rs2,rs1).eq.disulf(min(rs1,rs2))) then
      lk = crlk_idx(rs1)
      ll = 1 ! we can step through the crosslink excludes because they are sorted exactly the same
      do i=1,at(rs1)%na
        if (i.le.at(rs1)%nbb) then
          ii = at(rs1)%bb(i)
        else
          ii = at(rs1)%sc(i-at(rs1)%nbb)
        end if
        do j=1,at(rs2)%na
          if (j.le.at(rs2)%nbb) then
            kk = at(rs2)%bb(j)
          else
            kk = at(rs2)%sc(j-at(rs2)%nbb)
          end if
          if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) then
            if (crosslink(lk)%is14in(ll).EQV..false.) then
              ll = ll + 1
              cycle
            end if
          end if
          dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii) + svec(:)
          d2(k) = sum(dvec(:,k)**2) 
          if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
!           we must not have any special exclusion rules for crosslinks
            if (ll.le.crosslink(lk)%nrsin) then
              if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) then ! must be a 1-4
                term0(k) = fudge_st_14*lj_eps_14(attyp(ii),attyp(kk))
                terms(k) = lj_sig_14(attyp(ii),attyp(kk))
                ll = ll + 1
              else
                term0(k) = lj_eps(attyp(ii),attyp(kk))
                terms(k) = lj_sig(attyp(ii),attyp(kk))
              end if
            else
              term0(k) = lj_eps(attyp(ii),attyp(kk))
              terms(k) = lj_sig(attyp(ii),attyp(kk))
            end if
            if (term0(k).gt.0.0) then
              idx2(k,1) = ii
              idx2(k,2) = kk
              k = k + 1
            end if
          else if (ll.le.crosslink(lk)%nrsin) then
            if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) ll = ll + 1
          end if
        end do
      end do
!
    else if (rs1.eq.rs2) then
!
      do i=1,nrsintra(rs1)
        ii = iaa(rs1)%atin(i,1)
        kk = iaa(rs1)%atin(i,2)
        dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii)
        d2(k) = sum(dvec(:,k)**2) 
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
          term0(k) = fudge(rs1)%rsin(i)*fudge(rs1)%rsin_lje(i)
          if (term0(k).gt.0.0) then
            terms(k) = fudge(rs1)%rsin_ljs(i)
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
!
    else if (abs(rs1-rs2).eq.1) then
!
      rlo = min(rs1,rs2)
      if (rlo.eq.rs2) svec(:) = -svec(:)
      do i=1,nrsnb(rlo)
        ii = iaa(rlo)%atnb(i,1)
        kk = iaa(rlo)%atnb(i,2)
        dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii) + svec(:)
        d2(k) = sum(dvec(:,k)**2) 
        if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
          term0(k) = fudge(rlo)%rsnb(i)*fudge(rlo)%rsnb_lje(i)
          if (term0(k).gt.0.0) then
            terms(k) = fudge(rlo)%rsnb_ljs(i)
            idx2(k,1) = ii
            idx2(k,2) = kk
            k = k + 1
          end if
        end if
      end do
!
    else
!
      do i=1,at(rs1)%na
        if (i.le.at(rs1)%nbb) then
          ii = at(rs1)%bb(i)
        else
          ii = at(rs1)%sc(i-at(rs1)%nbb)
        end if
        do j=1,at(rs2)%na
          if (j.le.at(rs2)%nbb) then
            kk = at(rs2)%bb(j)
          else
            kk = at(rs2)%sc(j-at(rs2)%nbb)
          end if
          dvec(:,k) = atinfo(1:3,kk) - atinfo(1:3,ii) + svec(:)
          d2(k) = sum(dvec(:,k)**2) 
          if (((cut.EQV..true.).AND.(d2(k).lt.mcnb_cutoff2)).OR.(cut.EQV..false.)) then
            term0(k) = lj_eps(attyp(ii),attyp(kk))
            if (term0(k).gt.0.0) then
              terms(k) = lj_sig(attyp(ii),attyp(kk))
              idx2(k,1) = ii
              idx2(k,2) = kk
              k = k + 1
            end if
          end if
        end do
      end do
    end if
!
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      for_tmp(:,1:alcsz) = 0.0
      if (doid2.EQV..true.) then
        id2(1:alcsz) = 1.0/d2(1:alcsz)
        term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3
      end if
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz)
        evec(1) = evec(1) + par_FEG2(1)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(1)*48.0
        termr(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        do k=1,alcsz
          for_tmp(:,k) = for_tmp(:,k) - termr(k)*dvec(:,k)
        end do
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        term3(1:alcsz) = term0(1:alcsz)*term2(1:alcsz)
        evec(3) = evec(3) - par_FEG2(5)*4.0*sum(term3(1:alcsz))
        pfac = par_FEG2(5)*24.0
        termr(1:alcsz) = pfac*term3(1:alcsz)*term2(1:alcsz)*term1(1:alcsz)*id2(1:alcsz)
        do k=1,alcsz
          for_tmp(:,k) = for_tmp(:,k) + termr(k)*dvec(:,k)
        end do
      end if
    end if
!
  end if
!
  do k=1,alcsz
    ca_f(:,idx2(k,1)) = ca_f(:,idx2(k,1)) + for_tmp(:,k)
    ca_f(:,idx2(k,2)) = ca_f(:,idx2(k,2)) - for_tmp(:,k)
  end do
!
end
!
!---------------------------------------------------------------------------------
!
! just the LJ potential with support for implicit solvent models
!
subroutine Vforce_LJIMP_C(evec,rs1,ca_f,s_sv,si,lora,hira3,hira2,which)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,lo,k,i,shii,jj,ki,kj,imol,chnksz,chnk,savnbs,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n),s_sv(n)
  type(t_sisa) si(n)
!
  RTYPE term6(VCI_SIZE),xt0(VCI_SIZE),dvec(VCI_SIZE,3)
  RTYPE foin(VCI_SIZE,3),d1(VCI_SIZE),xt1(VCI_SIZE),xt2(VCI_SIZE),xt3(VCI_SIZE),id2(VCI_SIZE),xt4(VCI_SIZE)
  integer savidx(VCI_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(7,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  integer ka(hira2),kb(hira2),kc(hira2)
  RTYPE term0(hira2),termv(hira2),termx(hira2)
  integer sh(hira3-lora+1),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8
  RTYPE incs(3),derincs(3)
!  
! probably seg-faulted already if hira is 0
  hira = hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atsavinfo(1,ii) + par_IMPSOLV(1)
    ixis(5,j) = atsavinfo(2,ii) ! scaled volume or volume reduction factor (depends on par_IMPSOLV2(1))
    ixis(6,j) = atsavinfo(1,ii) ! radius
    if (par_IMPSOLV2(1).eq.1) then
      ixis(7,j) = 2.0*atsavinfo(1,ii)
    else if (par_IMPSOLV2(1).eq.2) then
      ixis(7,j) = atsavinfo(3,ii) ! volume
    end if
  end do
  ii = refat(rs1)
!
  if ((which.eq.1).OR.(which.eq.3)) then ! standard range, recompute shift vectors
    if (which.eq.1) then
      k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    else
      k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira3)
    end if
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else if (which.eq.4) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
! reuse parameters
  lop1 = 4.0*scale_IPP
  lop2 = 4.0*scale_attLJ
  lop3 = electric*scale_POLAR
!
! this loop has data dependency and a conditional -> as few operations as possible in here
  lo = 1
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
  do k=1,hira
    do j=k1(k),k1(k)+sh(k)
      jj = ka(k) + j - k1(k)
      kb(jj) = attyp(j)
      kc(jj) = j
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
      term0(jj) = atsavinfo(1,j) ! radius
      termv(jj) = atsavinfo(2,j) ! scaled volume or volume red. factor
      if (par_IMPSOLV2(1).eq.2) termx(jj) = atsavinfo(3,j) ! volume
    end do
  end do
!
  do chnk=1,hira2,VCI_SIZE
    hi = min(chnk + VCI_SIZEM1,hira2)
    chnksz = hi-chnk+1
    dvec(1:chnksz,1) = dveci(chnk:hi,1)
    dvec(1:chnksz,2) = dveci(chnk:hi,2)
    dvec(1:chnksz,3) = dveci(chnk:hi,3)
    do i=1,at(rs1)%na
!     inverse distance squared and LJ potential
      xt1(1:chnksz) = ((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                   (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                   (dvec(1:chnksz,3)-ixis(3,i))**2)
!      d1(1:chnksz) = sqrt(xt1(1:chnksz))
      id2(1:chnksz) = 1.0/xt1(1:chnksz)
      term6(1:chnksz) = (lj_sig(kb(chnk:hi),k0(i,1))*id2(1:chnksz))**3
      foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
      foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
      term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
      evec(1) = evec(1) + sum(foin(1:chnksz,3))
      evec(3) = evec(3) - sum(term6(1:chnksz))
      term6(1:chnksz) = id2(1:chnksz)*(12.0*foin(1:chnksz,3)-6.0*term6(1:chnksz))
!     force 
      foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
      foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
      foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
      ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
      ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
      ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
      for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
      for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
      for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
!
      term6(1:chnksz) = (term0(chnk:hi) + ixis(4,i))
      xt0(1:chnksz) = term6(1:chnksz)**2
      savnbs = 0
      do k=1,chnksz
        if (xt1(k).lt.xt0(k)) then
          ki = chnk+k-1
          savnbs = savnbs + 1
          savidx(savnbs) = kc(ki)
!         distance vector
          foin(savnbs,1) = (dvec(k,1)-ixis(1,i))
          foin(savnbs,2) = (dvec(k,2)-ixis(2,i))
          foin(savnbs,3) = (dvec(k,3)-ixis(3,i))
!         inverse distance squared
          id2(savnbs) = id2(k)
          d1(savnbs) = xt1(k)
!         maximum permissible distance (effective cutoff)
          term6(savnbs) = term6(k)
!         size/volume parameters for atoms k
          if (par_IMPSOLV2(1).eq.1) then
            xt0(savnbs) = termv(ki) ! scaled volume
            xt1(savnbs) = term0(ki) ! radius
          else if (par_IMPSOLV2(1).eq.2) then
            xt0(savnbs) = termv(ki) ! volume red factor
            xt1(savnbs) = term0(ki) ! radius
            xt4(savnbs) = termx(ki) ! volume
          end if
        end if
      end do 
      if (savnbs.eq.0) cycle
!
      if ((savnbs+si(k0(i,2))%nix).gt.si(k0(i,2))%alsz) then
        k = max(VSN_RESIZE,savnbs+si(k0(i,2))%nix-si(k0(i,2))%alsz) + 1
        call Vsisa_resize(k0(i,2),si,k)
      end if
      if (par_IMPSOLV2(1).eq.1) then
        d1(1:savnbs) = sqrt(d1(1:savnbs))           ! distance
        xt2(1:savnbs) = (ixis(6,i) + xt1(1:savnbs)) ! sum of radii
        xt3(1:savnbs) = id2(1:savnbs)*d1(1:savnbs)  ! inverse distance 
        do k=1,savnbs
          derincs(2) = 0.0
          incs(2) = -1.0
          lop4 = 2.0*xt1(k)
          if (d1(k).gt.(term6(k)-lop4)) then
            lop5 = 1.0/lop4
            incs(1) = -(term6(k)-d1(k))*lop5
            derincs(1) = -xt0(k)*lop5
          else if (d1(k).lt.xt2(k)) then
            lop6 = 1.0/xt2(k) ! rare
            incs(1) = -d1(k)*lop6
            derincs(1) = xt0(k)*lop6
            incs(2) = incs(1)
            derincs(2) = -ixis(5,i)*lop6
          else
            incs(1) = -1.0
            derincs(1) = 0.0
          end if
          if (d1(k).gt.(term6(k)-ixis(7,i))) then
            lop7 = 1.0/ixis(7,i)
            incs(2) = -(term6(k)-d1(k))*lop7
            derincs(2) = ixis(5,i)*lop7
          end if
          s_sv(k0(i,2)) = s_sv(k0(i,2)) + incs(1)*xt0(k)
          s_sv(savidx(k)) = s_sv(savidx(k)) + incs(2)*ixis(5,i)
          if (derincs(1).ne.0.0) then
            lop8 = derincs(1)*xt3(k)
            si(k0(i,2))%nix = si(k0(i,2))%nix + 1
            si(k0(i,2))%ix(si(k0(i,2))%nix) = savidx(k)
            si(k0(i,2))%dr(1,si(k0(i,2))%nix) = foin(k,1)*lop8
            si(k0(i,2))%dr(2,si(k0(i,2))%nix) = foin(k,2)*lop8
            si(k0(i,2))%dr(3,si(k0(i,2))%nix) = foin(k,3)*lop8
          end if
          if (derincs(2).ne.0.0) then
            lop8 = derincs(2)*xt3(k)
            if ((si(savidx(k))%nix+1).gt.si(savidx(k))%alsz) call Vsisa_resize(savidx(k),si,VSN_RESIZE)
            si(savidx(k))%nix = si(savidx(k))%nix + 1
            si(savidx(k))%ix(si(savidx(k))%nix) = k0(i,2)
            si(savidx(k))%dr(1,si(savidx(k))%nix) = foin(k,1)*lop8
            si(savidx(k))%dr(2,si(savidx(k))%nix) = foin(k,2)*lop8
            si(savidx(k))%dr(3,si(savidx(k))%nix) = foin(k,3)*lop8
          end if
        end do
      else if (par_IMPSOLV2(1).eq.2) then
        d1(1:savnbs) = sqrt(d1(1:savnbs))           ! distance
        xt3(1:savnbs) = id2(1:savnbs)*d1(1:savnbs)  ! inverse distance 
        do k=1,savnbs
          call sphere_overlaps_wder(ixis(6,i),xt1(k),d1(k),par_IMPSOLV(1),ixis(7,i),xt4(k),incs,derincs)
          if (incs(2).gt.0.0) then
            s_sv(k0(i,2)) = s_sv(k0(i,2)) - (incs(2)-incs(1))*xt0(k)
            if ((derincs(1).ne.0.0).OR.(derincs(2).ne.0.0)) then
              lop8 = (derincs(2)-derincs(1))*xt0(k)*xt3(k)
              si(k0(i,2))%nix = si(k0(i,2))%nix + 1
              si(k0(i,2))%ix(si(k0(i,2))%nix) = savidx(k)
              si(k0(i,2))%dr(1,si(k0(i,2))%nix) = foin(k,1)*lop8
              si(k0(i,2))%dr(2,si(k0(i,2))%nix) = foin(k,2)*lop8
              si(k0(i,2))%dr(3,si(k0(i,2))%nix) = foin(k,3)*lop8
            end if
          end if
          if (incs(3).gt.0.0) then
            s_sv(savidx(k)) = s_sv(savidx(k)) - (incs(3)-incs(1))*ixis(5,i)
            if ((derincs(1).ne.0.0).OR.(derincs(3).ne.0.0)) then
              lop8 = (derincs(1)-derincs(3))*ixis(5,i)*xt3(k)
              if ((si(savidx(k))%nix+1).gt.si(savidx(k))%alsz) call Vsisa_resize(savidx(k),si,VSN_RESIZE)
              si(savidx(k))%nix = si(savidx(k))%nix + 1
              si(savidx(k))%ix(si(savidx(k))%nix) = k0(i,2)
              si(savidx(k))%dr(1,si(savidx(k))%nix) = foin(k,1)*lop8
              si(savidx(k))%dr(2,si(savidx(k))%nix) = foin(k,2)*lop8
              si(savidx(k))%dr(3,si(savidx(k))%nix) = foin(k,3)*lop8
            end if
          end if
        end do
      end if
    end do
  end do
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the vectorized version of force_rsp_long
! note that unlike in en_rsp_long, long-range electrostatics is handled externally(!)
!
subroutine Vforce_rsp_long(evec,rs1,rs2,ca_f,sum_s)
!
  use iounit
  use energies
  use atoms
  use polypep
  use inter
  use units
  use tabpot
  use molecule
  use sequen
  use cutoffs
  use system
  use grandensembles
  use forces
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer rs,i,j,k,ii,kk,rhi,rlo,imol1,imol2,alcsz
  RTYPE sum_s(n)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac,pfac2,pfac3,pfac4,dum1,dd1,dd2,hlp1,hlp2,hlp3,hlp4,dhlp1,dhlp2,dhlp3
  RTYPE d2(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE term4(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na),termqi(at(rs1)%na*at(rs2)%na),termqk(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na),terms1(at(rs1)%na*at(rs2)%na),terms2(at(rs1)%na*at(rs2)%na)
  RTYPE for_tmpi(at(rs1)%na*at(rs2)%na,3)
  RTYPE for_tmpk(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdi(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdk(at(rs1)%na*at(rs2)%na,3)
  RTYPE termf(at(rs1)%na*at(rs2)%na,3)
  integer tbin(at(rs1)%na*at(rs2)%na)
  logical ismember
  RTYPE ca_f(n,3)
!
  if ((ideal_run.EQV..true.).OR.&
 &    ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.))) then
    return
  end if
!  
! Prevent non-present molecules from interacting with each other
! or with present molecules.  Note that intra-nonpresent-molecule
! interactions are still allowed, in order to correctly simulate
! an infinite-dilution implicit solvent reference state particle
! bath -> partially redundant overlap with use_FEG
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    imol1 = molofrs(rs1)
    imol2 = molofrs(rs2)
    if ((imol1.ne.imol2).AND.((.not.ismember(ispresent,imol1)).OR. &
 &   (.not.ismember(ispresent,imol2)))) return
  end if
!
! in FEG we branch out if precisely one of the residues is ghosting
! always remember that the interaction of two ghost-residues (incl self!)
! is the full Hamiltonian, not the ghosted one
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
   if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Vforce_rsp_long_feg(evec,rs1,rs2,ca_f,sum_s)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Vforce_rsp_long_feg(evec,rs1,rs2,ca_f,sum_s)
        return
      end if
    end if
  end if
!
  svec(:) = 0.0
!
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
    rs = rs1
    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*scale_POLAR
      alcsz = nrpolintra(rs)
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2) 
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(:,ii)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(:,ii)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)
              tsdi(i,1:3) = scrq_dr(:,ii)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)
              tsdi(i,1:3) = scrq_dr(:,ii)
              termr(i) = atr(ii)+atr(kk)
            end do
          end if
        else
          do i=1,nrpolintra(rs)
            ii = iaa(rs)%polin(i,1)
            kk = iaa(rs)%polin(i,2)
            dvec(i,1) = x(kk) - x(ii)
            dvec(i,2) = y(kk) - y(ii)
            dvec(i,3) = z(kk) - z(ii)
            termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
!         with a minor if-statement, the normal models (1,2,5,6) can be handled together
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              sum_s(ii) = sum_s(ii) + termqi(i)
              sum_s(kk) = sum_s(kk) + termqk(i)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3) + for_tmpi(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3) + for_tmpk(i,1:3)
            end do
!         the simplest one is always the straight distance-dependence (only complexity
!         comes from contact dielectric)
          else if (scrq_model.eq.4) then
            do i=1,nrpolintra(rs)
              if (d1(i).ge.terms(i)) then
                term1(i) = id1(i)
                term0(i) = 2.0*id1(i)
              else
                term1(i) = 1.0/terms(i)
                term0(i) = term1(i)
              end if
            end do
            term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(term2(1:alcsz))
            k = 0
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3)
            end do
!         the remaining models (3,7,8,9) have spliced-in distance-dependencies
!         this leads to bulkier-looking code but does introduce that much extra complexity compared
!         to the respective reference models (1,5,6,2 in that order)
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              if (abs(term1(i)).gt.abs(term2(i))) then
                term4(i) = term1(i)
                term3(i) = 1.0
              else
                dd1 = (d1(i) - termr(i))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(i) = (1.0 - par_IMPSOLV(9))
                  term4(i) = term3(i)*term1(i) + par_IMPSOLV(9)*term2(i)
                  term0(i) = term3(i)*term0(i)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(i) = (1.0 - dum1)
                  term4(i) = term3(i)*term1(i) + dum1*term2(i)
                  term0(i) = term3(i)*term0(i)
                  dum1 = -scale_POLAR*id2(i)*(term2(i)-term1(i))*pfac4
                  for_tmpi(i,1) = dum1*dvec(i,1)
                  for_tmpi(i,2) = dum1*dvec(i,2)
                  for_tmpi(i,3) = dum1*dvec(i,3)
                else
                  term4(i) = term1(i)
                  term3(i) = 1.0
                end if
              end if
              sum_s(ii) = sum_s(ii) + term0(i)*termqi(i)
              sum_s(kk) = sum_s(kk) + term0(i)*termqk(i)
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            k = 0
            do i=1,nrpolintra(rs)
              ii = iaa(rs)%polin(i,1)
              kk = iaa(rs)%polin(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3) + for_tmpi(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3) + for_tmpk(i,1:3)
            end do
          end if
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          do i=1,nrpolintra(rs)
            ii = iaa(rs)%polin(i,1)
            kk = iaa(rs)%polin(i,2)
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
          end do
        end if
      end if
!
    end if
!
!   note that this might be terribly inefficient if POLAR and TABUL are used concurrently
    if (use_TABUL.EQV..true.) then
!
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        k = 0
        alcsz = tbp%rsvec(rlo)-tbp%rsmat(rlo,rhi)+1
        do i=tbp%rsmat(rlo,rhi),tbp%rsvec(rlo)
          ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
        term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsvec(rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
            hlp1 = 2.0*term3(k) - 3.0*term2(k)
            hlp2 = term3(k) - term2(k)
            hlp3 = hlp2 - term2(k) + term1(k)
            dhlp1 = 6.0*term2(k) - 6.0*term1(k)
            dhlp2 = 3.0*term2(k) - 2.0*term1(k)
            dhlp3 = dhlp2 - 2.0*term1(k) + 1.0
            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
            ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k,1:3)
          else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!
!
!
! neighboring residues
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC evtl.y
!
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*scale_POLAR
      alcsz = nrpolnb(rs)
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2) 
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(:,ii)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms(i) = scrq(ii)*scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)*scrq(ii)
              tsdi(i,1:3) = scrq_dr(:,ii)*scrq(kk)
              termqk(i) = termq(i)*scrq(ii)
              termqi(i) = termq(i)*scrq(kk)
              termr(i) = atr(ii)+atr(kk)
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)
              tsdi(i,1:3) = scrq_dr(:,ii)
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              dvec(i,1) = x(kk) - x(ii)
              dvec(i,2) = y(kk) - y(ii)
              dvec(i,3) = z(kk) - z(ii)
              termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
              terms1(i) = scrq(ii)
              terms2(i) = scrq(kk)
              tsdk(i,1:3) = scrq_dr(:,kk)
              tsdi(i,1:3) = scrq_dr(:,ii)
              termr(i) = atr(ii)+atr(kk)
            end do
          end if
        else
          do i=1,nrpolnb(rs)
            ii = iaa(rs)%polnb(i,1)
            kk = iaa(rs)%polnb(i,2)
            dvec(i,1) = x(kk) - x(ii)
            dvec(i,2) = y(kk) - y(ii)
            dvec(i,3) = z(kk) - z(ii)
            termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              sum_s(ii) = sum_s(ii) + termqi(i)
              sum_s(kk) = sum_s(kk) + termqk(i)
              ca_f(ii,1:3) = ca_f(ii,1:3) + for_tmpi(i,1:3) - termf(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + for_tmpk(i,1:3) + termf(i,1:3)
            end do
          else if (scrq_model.eq.4) then
            do i=1,nrpolnb(rs)
              if (d1(i).ge.terms(i)) then
                term1(i) = id1(i)
                term0(i) = 2.0*id1(i)
              else
                term1(i) = 1.0/terms(i)
                term0(i) = term1(i)
              end if
            end do
            term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(term2(1:alcsz))
            k = 0
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3)
            end do
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              if (abs(term1(i)).gt.abs(term2(i))) then
                term4(i) = term1(i)
                term3(i) = 1.0
              else
                dd1 = (d1(i) - termr(i))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(i) = (1.0 - par_IMPSOLV(9))
                  term4(i) = term3(i)*term1(i) + par_IMPSOLV(9)*term2(i)
                  term0(i) = term3(i)*term0(i)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(i) = (1.0 - dum1)
                  term4(i) = term3(i)*term1(i) + dum1*term2(i)
                  term0(i) = term3(i)*term0(i)
                  dum1 = -scale_POLAR*id2(i)*(term2(i)-term1(i))*pfac4
                  for_tmpi(i,1) = dum1*dvec(i,1)
                  for_tmpi(i,2) = dum1*dvec(i,2)
                  for_tmpi(i,3) = dum1*dvec(i,3)
                else
                  term4(i) = term1(i)
                  term3(i) = 1.0
                end if
              end if
              sum_s(ii) = sum_s(ii) + term0(i)*termqi(i)
              sum_s(kk) = sum_s(kk) + term0(i)*termqk(i)
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            k = 0
            do i=1,nrpolnb(rs)
              ii = iaa(rs)%polnb(i,1)
              kk = iaa(rs)%polnb(i,2)
              ca_f(ii,1:3) = ca_f(ii,1:3) - termf(i,1:3) + for_tmpi(i,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + termf(i,1:3) + for_tmpk(i,1:3)
            end do
          end if
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          do i=1,nrpolnb(rs)
            ii = iaa(rs)%polnb(i,1)
            kk = iaa(rs)%polnb(i,2)
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
          end do
        end if
      end if
!
    end if
!
    if (use_TABUL.EQV..true.) then
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        k = 0
        alcsz = tbp%rsmat(rhi,rlo) - tbp%rsmat(rlo,rhi) + 1
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
!       note we need to correct: the svec is adjusted such that the
!       sign is correct when using the polnb(rs)-array, but that TAB
!       assumes rs1,rs2 order
        if (rs.eq.rs2) then
          dvec(1:alcsz,1) = dvec(1:alcsz,1) - svec(1)
          dvec(1:alcsz,2) = dvec(1:alcsz,2) - svec(2)
          dvec(1:alcsz,3) = dvec(1:alcsz,3) - svec(3)
        else
          dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
          dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
          dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        end if
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
        term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)

            hlp1 = 2.0*term3(k) - 3.0*term2(k)
            hlp2 = term3(k) - term2(k)
            hlp3 = hlp2 - term2(k) + term1(k)
            dhlp1 = 6.0*term2(k) - 6.0*term1(k)
            dhlp2 = 3.0*term2(k) - 2.0*term1(k)
            dhlp3 = dhlp2 - 2.0*term1(k) + 1.0
            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
            ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k,1:3)
          else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!
!
!
! all other residue pairs
!
 67   format(4g14.6)
  else
!
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!
    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*scale_POLAR
      alcsz = at(rs1)%npol*at(rs2)%npol
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms(k) = scrq(ii)*scrq(kk)
                tsdk(k,1:3) = scrq_dr(:,kk)*scrq(ii)
                tsdi(k,1:3) = scrq_dr(:,ii)*scrq(kk)
                termqk(k) = termq(k)*scrq(ii)
                termqi(k) = termq(k)*scrq(kk)
              end do
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms(k) = atr(ii)+atr(kk)
              end do
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms(k) = scrq(ii)*scrq(kk)
                tsdk(k,1:3) = scrq_dr(:,kk)*scrq(ii)
                tsdi(k,1:3) = scrq_dr(:,ii)*scrq(kk)
                termqk(k) = termq(k)*scrq(ii)
                termqi(k) = termq(k)*scrq(kk)
                termr(k) = atr(ii)+atr(kk)
              end do
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms1(k) = scrq(ii)
                terms2(k) = scrq(kk)
                tsdk(k,1:3) = scrq_dr(:,kk)
                tsdi(k,1:3) = scrq_dr(:,ii)
              end do
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                dvec(k,1) = x(kk) - x(ii)
                dvec(k,2) = y(kk) - y(ii)
                dvec(k,3) = z(kk) - z(ii)
                termq(k) = atq(ii)*atq(kk)
                terms1(k) = scrq(ii)
                terms2(k) = scrq(kk)
                tsdk(k,1:3) = scrq_dr(:,kk)
                tsdi(k,1:3) = scrq_dr(:,ii)
                termr(k) = atr(ii)+atr(kk)
              end do
            end do
          end if
        else
          do i=1,at(rs1)%npol
            ii = at(rs1)%pol(i)
            do j=1,at(rs2)%npol
              kk = at(rs2)%pol(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
            end do
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                sum_s(ii) = sum_s(ii) + termqi(k)
                sum_s(kk) = sum_s(kk) + termqk(k)
                ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3) + for_tmpi(k,1:3)
                ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3) + for_tmpk(k,1:3)
              end do
            end do
          else if (scrq_model.eq.4) then
            do i=1,at(rs1)%npol
              do j=1,at(rs2)%npol
                k = k + 1
                if (d1(k).ge.terms(k)) then
                  term1(k) = id1(k)
                  term0(k) = 2.0*id1(k)
                else
                  term1(k) = 1.0/terms(k)
                  term0(k) = term1(k)
                end if
              end do
            end do
            term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(term2(1:alcsz))
            k = 0
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3)
                ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3)
              end do
            end do
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                if (abs(term1(k)).gt.abs(term2(k))) then
                  term4(k) = term1(k)
                  term3(k) = 1.0
                else
                  dd1 = (d1(k) - termr(k))*pfac3
                  dd2 = 1.0 - dd1
                  if (dd1.lt.0.0) then
                    term3(k) = (1.0 - par_IMPSOLV(9))
                    term4(k) = term3(k)*term1(k) + par_IMPSOLV(9)*term2(k)
                    term0(k) = term3(k)*term0(k)
                  else if (dd2.gt.0.0) then
                    dum1 = par_IMPSOLV(9)*dd2
                    term3(k) = (1.0 - dum1)
                    term4(k) = term3(k)*term1(k) + dum1*term2(k)
                    term0(k) = term3(k)*term0(k)
                    dum1 = -scale_POLAR*id2(k)*(term2(k)-term1(k))*pfac4
                    for_tmpi(k,1) = dum1*dvec(k,1)
                    for_tmpi(k,2) = dum1*dvec(k,2)
                    for_tmpi(k,3) = dum1*dvec(k,3)
                  else
                    term4(k) = term1(k)
                    term3(k) = 1.0
                  end if
                end if
                sum_s(ii) = sum_s(ii) + term0(k)*termqi(k)
                sum_s(kk) = sum_s(kk) + term0(k)*termqk(k)
              end do
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(term1(1:alcsz))
            term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
            k = 0
            do i=1,at(rs1)%npol
              ii = at(rs1)%pol(i)
              do j=1,at(rs2)%npol
                kk = at(rs2)%pol(j)
                k = k + 1
                ca_f(ii,1:3) = ca_f(ii,1:3) - termf(k,1:3) + for_tmpi(k,1:3)
                ca_f(kk,1:3) = ca_f(kk,1:3) + termf(k,1:3) + for_tmpk(k,1:3)
              end do
            end do
          end if
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          k = 0
          do i=1,at(rs1)%npol
            ii = at(rs1)%pol(i)
            do j=1,at(rs2)%npol
              kk = at(rs2)%pol(j)
              k = k + 1
              ca_f(ii,1:3) = ca_f(ii,1:3) - term2(k)*dvec(k,1:3)
              ca_f(kk,1:3) = ca_f(kk,1:3) + term2(k)*dvec(k,1:3)
            end do
          end do
        end if
      end if
!
    end if
!
    if (use_TABUL.EQV..true.) then
!
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        alcsz = tbp%rsmat(rhi,rlo) - tbp%rsmat(rlo,rhi) + 1
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
        term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
!        id1(1:alcsz) = 1.0/d1(1:alcsz)
!        tbin(1:alcsz) = floor((d1(1:alcsz)-tbp%dis(1))/tbp%res) + 1
!        pfac = scale_TABUL*(1.0/tbp%res)
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)

            hlp1 = 2.0*term3(k) - 3.0*term2(k)
            hlp2 = term3(k) - term2(k)
            hlp3 = hlp2 - term2(k) + term1(k)
            dhlp1 = 6.0*term2(k) - 6.0*term1(k)
            dhlp2 = 3.0*term2(k) - 2.0*term1(k)
            dhlp3 = dhlp2 - 2.0*term1(k) + 1.0

            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))

!            term1(k) = (tbp%dis(tbin(k)+1) - d1(k))*tbp%pot(tbp%lst(i,3),tbin(k))
!            term2(k) = (d1(k) - tbp%dis(tbin(k)))*tbp%pot(tbp%lst(i,3),tbin(k)+1)
!            term3(k) = pfac*(term1(k)+term2(k))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
!            term4(k) = id1(k)*pfac*(tbp%pot(tbp%lst(i,3),tbin(k)+1)-tbp%pot(tbp%lst(i,3),tbin(k)))

            ca_f(ii,1:3) = ca_f(ii,1:3) + hlp4*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) - hlp4*dvec(k,1:3)
         else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!       
  end if
!
!
end
!
!--------------------------------------------------------------------------
!
! this subroutine computes all monopole-monopole and monopole-dipole, but not dipole-dipole terms
! between rs1 and rs2
!
subroutine Vforce_rsp_long_lrel(evec,rs1,rs2,svec,ca_f,sum_s)
!
  use atoms
  use polypep
  use energies
  use molecule
  use sequen
  use cutoffs
  use units
  use forces
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  RTYPE, INTENT(IN):: svec(3)
!
  integer i,j,k,ii,kk,alcsz,g1,g2
  RTYPE sum_s(n)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE pfac,pfac2,pfac3,pfac4,dum1,dd1,dd2
  RTYPE d2(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE term3(at(rs1)%na*at(rs2)%na)
  RTYPE term4(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na),termqi(at(rs1)%na*at(rs2)%na),termqk(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na),terms1(at(rs1)%na*at(rs2)%na),terms2(at(rs1)%na*at(rs2)%na)
  RTYPE for_tmpi(at(rs1)%na*at(rs2)%na,3)
  RTYPE for_tmpk(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdi(at(rs1)%na*at(rs2)%na,3)
  RTYPE tsdk(at(rs1)%na*at(rs2)%na,3)
  RTYPE termf(at(rs1)%na*at(rs2)%na,3)
  RTYPE ca_f(3,n)
  logical dogh
!
  if ((abs(rs2-rs1).eq.1).AND.(molofrs(rs1).eq.molofrs(rs2))) then
    write(ilog,*) 'Warning. Neighboring residues in the same molecule are beyond mid-range &
 &cutoff. This indicates an unstable simulation or a bug. Expect run to crash soon.'
  end if
!
  if (use_FEG.EQV..true.) then
    if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &    ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
      dogh = .false.
    else
      dogh = .true.
    end if
    if ((dogh.EQV..true.).AND.(use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
      write(ilog,*) 'Fatal. Missing support for combination of FEG and the ABSINTH&
 & implicit solvation model in Vforce_rsp_long_lrel(...). This is an omission bug.'
      call fexit()
    end if
  else
    dogh = .false.
  end if
!
  if ((use_POLAR.EQV..true.).AND.(lrel_md.eq.5)) then
!
!   first set the necessary atom-specific parameters (mimimal)
    k = 0
    pfac = electric*scale_POLAR
    alcsz = at(rs1)%npol*at(rs2)%npol
    if (alcsz.gt.0) then
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle ! dipole-dipole
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms(k) = scrq(ii)*scrq(kk)
              tsdk(k,1:3) = scrq_dr(:,kk)*scrq(ii)
              tsdi(k,1:3) = scrq_dr(:,ii)*scrq(kk)
              termqk(k) = termq(k)*scrq(ii)
              termqi(k) = termq(k)*scrq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
          pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms(k) = atr(ii)+atr(kk)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          pfac3 = 1./par_IMPSOLV(1)
          pfac4 = par_IMPSOLV(9)*pfac3
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms(k) = scrq(ii)*scrq(kk)
              tsdk(k,1:3) = scrq_dr(:,kk)*scrq(ii)
              tsdi(k,1:3) = scrq_dr(:,ii)*scrq(kk)
              termqk(k) = termq(k)*scrq(ii)
              termqi(k) = termq(k)*scrq(kk)
              termr(k) = atr(ii)+atr(kk)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms1(k) = scrq(ii)
              terms2(k) = scrq(kk)
              tsdk(k,1:3) = scrq_dr(:,kk)
              tsdi(k,1:3) = scrq_dr(:,ii)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          pfac3 = 1./par_IMPSOLV(1)
          pfac4 = par_IMPSOLV(9)*pfac3
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              termq(k) = atq(ii)*atq(kk)
              terms1(k) = scrq(ii)
              terms2(k) = scrq(kk)
              tsdk(k,1:3) = scrq_dr(:,kk)
              tsdi(k,1:3) = scrq_dr(:,ii)
              termr(k) = atr(ii)+atr(kk)
            end do
            end do
          end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
        do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
          if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
          do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            termq(k) = atq(ii)*atq(kk)
          end do
          end do
        end do
        end do
      end if
!     redefine alcsz
      alcsz = k
!     now perform the vectorizable bulk operations (maximal)
      dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
      dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
      dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
      d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                  dvec(1:alcsz,3)**2
      d1(1:alcsz) = sqrt(d2(1:alcsz))
      id1(1:alcsz) = 1.0/d1(1:alcsz)
      id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!     finish up -> screening model specificity is high
      if (use_IMPSOLV.EQV..true.) then
        k = 0
        if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
            call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
            tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
            tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
            tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
            tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
            tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
            tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
            termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
            termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
          end if
          term0(1:alcsz) = pfac*id1(1:alcsz)
          termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
          termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
          term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
          for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
          for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
          for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
          for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
          for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
          for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
          term1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
          termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
          termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              sum_s(ii) = sum_s(ii) + termqi(k)
              sum_s(kk) = sum_s(kk) + termqk(k)
              ca_f(:,ii) = ca_f(:,ii) - termf(k,1:3) + for_tmpi(k,1:3)
              ca_f(:,kk) = ca_f(:,kk) + termf(k,1:3) + for_tmpk(k,1:3)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
          do k=1,alcsz
            if (d1(k).ge.terms(k)) then
              term1(k) = id1(k)
              term0(k) = 2.0*id1(k)
            else
              term1(k) = 1.0/terms(k)
              term0(k) = term1(k)
            end if
          end do
          term2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
          termqi(1:alcsz) = term2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
          term2(1:alcsz) = term2(1:alcsz)*term1(1:alcsz)
          termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
          termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
          termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
          evec(6) = evec(6) + sum(term2(1:alcsz))
          k = 0
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              ca_f(:,ii) = ca_f(:,ii) - termf(k,1:3)
              ca_f(:,kk) = ca_f(:,kk) + termf(k,1:3)
            end do
            end do
          end do
          end do
        else
          if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
            call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
            call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
            tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
            tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
            tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
            tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
            tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
            tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
            termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
            termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
          end if
          term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
          term1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
          term2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
          for_tmpi(1:alcsz,1:3) = 0.0
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              if (abs(term1(k)).gt.abs(term2(k))) then
                term4(k) = term1(k)
                term3(k) = 1.0
              else
                dd1 = (d1(k) - termr(k))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(k) = (1.0 - par_IMPSOLV(9))
                  term4(k) = term3(k)*term1(k) + par_IMPSOLV(9)*term2(k)
                  term0(k) = term3(k)*term0(k)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(k) = (1.0 - dum1)
                  term4(k) = term3(k)*term1(k) + dum1*term2(k)
                  term0(k) = term3(k)*term0(k)
                  dum1 = -scale_POLAR*id2(k)*(term2(k)-term1(k))*pfac4
                  for_tmpi(k,1) = dum1*dvec(k,1)
                  for_tmpi(k,2) = dum1*dvec(k,2)
                  for_tmpi(k,3) = dum1*dvec(k,3)
                else
                  term4(k) = term1(k)
                  term3(k) = 1.0
                end if
              end if
              sum_s(ii) = sum_s(ii) + term0(k)*termqi(k)
              sum_s(kk) = sum_s(kk) + term0(k)*termqk(k)
            end do
            end do
          end do
          end do
          term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
          for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
          for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
          for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
          for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
          for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
          for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
          term1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
          termf(1:alcsz,1) = term2(1:alcsz)*dvec(1:alcsz,1)
          termf(1:alcsz,2) = term2(1:alcsz)*dvec(1:alcsz,2)
          termf(1:alcsz,3) = term2(1:alcsz)*dvec(1:alcsz,3)
          k = 0
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              ca_f(:,ii) = ca_f(:,ii) - termf(k,1:3) + for_tmpi(k,1:3)
              ca_f(:,kk) = ca_f(:,kk) + termf(k,1:3) + for_tmpk(k,1:3)
            end do
            end do
          end do
          end do
        end if
      else
        if (dogh.EQV..true.) then
          term1(1:alcsz) = electric*par_FEG2(9)*termq(1:alcsz)/(par_FEG2(10) + d1(1:alcsz))
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)/(par_FEG2(10) + d1(1:alcsz))
        else
          term1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(term1(1:alcsz))
          term2(1:alcsz) = term1(1:alcsz)*id2(1:alcsz)
        end if
        k = 0
        do g1=1,at(rs1)%ndpgrps
        do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
          if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle
          do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            ca_f(:,ii) = ca_f(:,ii) - term2(k)*dvec(k,1:3)
            ca_f(:,kk) = ca_f(:,kk) + term2(k)*dvec(k,1:3)
          end do
          end do
        end do
        end do
      end if
    end if
!
  else if ((use_POLAR.EQV..true.).AND.(lrel_md.eq.4)) then
!
!   this is handled through cgrp_ia instead
!
  end if
!
end
!
!---------------------------------------------------------------------------------------------
!
subroutine Vgenmu(r1,r2,i_genmu,ro,alcsz)
!
  implicit none
!
  integer i_genmu,alcsz
  RTYPE r1(alcsz),r2(alcsz),ro(alcsz),dum1
!
  if (i_genmu.eq.0) then
    ro(1:alcsz) = sqrt(r1(1:alcsz)*r2(1:alcsz))
  else if (i_genmu.eq.1) then
    ro(1:alcsz) = 0.5*(r1(1:alcsz) + r2(1:alcsz))
  else if (i_genmu.eq.-1) then
    ro(1:alcsz) = 1.0/(0.5*(1./r1(1:alcsz) + 1./r2(1:alcsz)))
  else if (i_genmu.eq.2) then
    ro(1:alcsz) = sqrt(0.5*(r1(1:alcsz)*r1(1:alcsz) + r2(1:alcsz)*r2(1:alcsz)))
  else
    dum1 = 1./(1.*i_genmu)
    ro(1:alcsz) = (0.5*(r1(1:alcsz)**(i_genmu) + r2(1:alcsz)**(i_genmu)))**(dum1)
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine Vgenmu_dr(r1,r2,i_genmu,dr1,dr2,alcsz)
!
  implicit none
!
  integer i_genmu,alcsz
  RTYPE dr1(alcsz),dr2(alcsz),r1(alcsz),r2(alcsz),help(alcsz),dum1,dum3
!
  if (i_genmu.eq.0) then
    help(1:alcsz) = 1.0/(2.0*sqrt(r1(1:alcsz)*r2(1:alcsz)))
    dr1(1:alcsz) = help(1:alcsz)*r2(1:alcsz)
    dr2(1:alcsz) = dr1(1:alcsz)*r1(1:alcsz)/r2(1:alcsz)
  else if (i_genmu.eq.1) then
    dr1(1:alcsz) = 0.5
    dr2(1:alcsz) = dr1(1:alcsz)
  else if (i_genmu.eq.-1) then
    help(1:alcsz) = 1./r1(1:alcsz) + 1./r2(1:alcsz)
    dr1(1:alcsz) = 1./(((r1(1:alcsz)*r1(1:alcsz))*help(1:alcsz))*(0.5*help(1:alcsz)))
    dr2(1:alcsz) = 1./(((r2(1:alcsz)*r2(1:alcsz))*help(1:alcsz))*(0.5*help(1:alcsz)))
  else if (i_genmu.eq.2) then
    help(1:alcsz) = 1.0/(2.0*sqrt(0.5*(r1(1:alcsz)*r1(1:alcsz) + r2(1:alcsz)*r2(1:alcsz))))
    dr1(1:alcsz) = r1(1:alcsz)*help(1:alcsz)
    dr2(1:alcsz) = r2(1:alcsz)*help(1:alcsz) ! dr1(1:alcsz)/r1(1:alcsz)
  else
    help(1:alcsz) = 0.5*(r1(1:alcsz)**i_genmu + r2(1:alcsz)**i_genmu)
    dum1 = 1./(1.*i_genmu) - 1.0
    dum3 = 1.*i_genmu - 1.0
    dr1(1:alcsz) = 0.5*((help(1:alcsz))**dum1)*((r1(1:alcsz)**dum3))
    dr2(1:alcsz) = 0.5*((help(1:alcsz))**dum1)*((r2(1:alcsz)**dum3))
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the vectorized version of force_rsp_long_feg
! this does not currently support IMPSOLV (note Ven_rsp_long does)
! support is very easy to add for screening models 1,2,5,6, but more annoying for the
! other ones
!
subroutine Vforce_rsp_long_feg(evec,rs1,rs2,ca_f,sum_s)
!
  use iounit
  use energies
  use atoms
  use polypep
  use inter
  use units
  use tabpot
  use molecule
  use sequen
  use cutoffs
  use system
  use grandensembles
  use forces
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer rs,i,j,k,ii,kk,alcsz
  RTYPE sum_s(n)
  RTYPE evec(MAXENERGYTERMS)
  RTYPE svec(3),pfac
  RTYPE d2(at(rs1)%na*at(rs2)%na),id1s(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na)
  RTYPE ca_f(n,3)
!
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
!
  svec(:) = 0.0
!
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
    rs = rs1
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = nrpolintra(rs)
      if (alcsz.gt.0) then
        do i=1,nrpolintra(rs)
          ii = iaa(rs)%polin(i,1)
          kk = iaa(rs)%polin(i,2)
          dvec(i,1) = x(kk) - x(ii)
          dvec(i,2) = y(kk) - y(ii)
          dvec(i,3) = z(kk) - z(ii)
          termq(i) = fudge(rs)%elin(i)*atq(ii)*atq(kk)
        end do
!       now perform the vectorizable bulk operations (maximal)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        term1(1:alcsz) = pfac*termq(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + sum(term1(1:alcsz))
        term2(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)*id1(1:alcsz)
        do i=1,nrpolintra(rs)
          ii = iaa(rs)%polin(i,1)
          kk = iaa(rs)%polin(i,2)
          ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
          ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
        end do
      end if
!
    end if
!
!
!
! neighboring residues
!
  else if (abs(rs1-rs2).eq.1) then
!
!   there is no topological requirement for neighboring residues to be close (pure sequence)
!   so we have to check for BC evtl.y
!
    if (rs1.gt.rs2) then
      rs = rs2
      call dis_bound_rs(rs2,rs1,svec)
    else
      rs = rs1
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = nrpolnb(rs)
      if (alcsz.gt.0) then
        do i=1,nrpolnb(rs)
          ii = iaa(rs)%polnb(i,1)
          kk = iaa(rs)%polnb(i,2)
          dvec(i,1) = x(kk) - x(ii)
          dvec(i,2) = y(kk) - y(ii)
          dvec(i,3) = z(kk) - z(ii)
          termq(i) = fudge(rs)%elnb(i)*atq(ii)*atq(kk)
        end do
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        term1(1:alcsz) = pfac*termq(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + sum(term1(1:alcsz))
        term2(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)*id1(1:alcsz)
        do i=1,nrpolnb(rs)
          ii = iaa(rs)%polnb(i,1)
          kk = iaa(rs)%polnb(i,2)
          ca_f(ii,1:3) = ca_f(ii,1:3) - term2(i)*dvec(i,1:3)
          ca_f(kk,1:3) = ca_f(kk,1:3) + term2(i)*dvec(i,1:3)
        end do
      end if
!
    end if
!
!
!
! all other residue pairs
!
 67   format(4g14.6)
  else
!
!   there is no topological relationship for remaining residues -> always check BC
    call dis_bound_rs(rs1,rs2,svec)
!
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
      alcsz = at(rs1)%npol*at(rs2)%npol
      if (alcsz.gt.0) then
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            kk = at(rs2)%pol(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            termq(k) = atq(ii)*atq(kk)
          end do
        end do
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id1s(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        term1(1:alcsz) = pfac*termq(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + sum(term1(1:alcsz))
        term2(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)*id1(1:alcsz)
        k = 0
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            kk = at(rs2)%pol(j)
            k = k + 1
            ca_f(ii,1:3) = ca_f(ii,1:3) - term2(k)*dvec(k,1:3)
            ca_f(kk,1:3) = ca_f(kk,1:3) + term2(k)*dvec(k,1:3)
          end do
        end do
      end if
!
    end if
!
!       
  end if
!
!
end
!
!-----------------------------------------------------------------------------------------------------
!
! the generally more useful replacement of Vforce_rsp_long (handles FEG, crosslinks, GC presence, and SR exclusions directly)
! call with rs1 <= rs2
!
subroutine Vforce_rsp_long2(evec,rs1,rs2,ca_f,sum_s)
!
  use molecule
  use iounit
  use energies
  use atoms
  use grandensembles
  use system
  use sequen
  use polypep
  use cutoffs
  use params
  use math
  use forces
  use units
  use tabpot
  use fyoc, ONLY: disulf
  use inter
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer imol1,imol2,i,j,k,ii,kk,alcsz,rhi,rlo,ll,lk
  logical isfeg,ismember
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n),sum_s(n)
!
  integer tbin(at(rs1)%na*at(rs2)%na),idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE svec(3),pfac,pfac2,pfac3,pfac4,dd1,dd2,dum1,hlp1,hlp2,hlp3,dhlp1,dhlp2,dhlp3,hlp4
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3),tsdi(at(rs1)%na*at(rs2)%na,3),tsdk(at(rs1)%na*at(rs2)%na,3)
  RTYPE for_tmpi(at(rs1)%na*at(rs2)%na,3),for_tmpk(at(rs1)%na*at(rs2)%na,3),termpre(at(rs1)%na*at(rs2)%na)
  RTYPE d1(at(rs1)%na*at(rs2)%na),d2(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE termq(at(rs1)%na*at(rs2)%na),terms(at(rs1)%na*at(rs2)%na),termqk(at(rs1)%na*at(rs2)%na),termqi(at(rs1)%na*at(rs2)%na)
  RTYPE terms1(at(rs1)%na*at(rs2)%na),terms2(at(rs1)%na*at(rs2)%na),termr(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na),termf(at(rs1)%na*at(rs2)%na,3),term3(at(rs1)%na*at(rs2)%na),term4(at(rs1)%na*at(rs2)%na)
!
  imol1 = molofrs(rs1)
  imol2 = molofrs(rs2)
!
  if (ideal_run.EQV..true.) return
!
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    if ((imol1.ne.imol2).AND.((.not.ismember(ispresent,imol1)).OR.(.not.ismember(ispresent,imol2)))) return
  end if
!
  isfeg = .false.
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
      if (((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &        ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) isfeg = .true.
    else if (fegmode.eq.2) then
      if ((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.)) isfeg = .true.
    end if
    if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) return
  end if
!
  if (rs1.ne.rs2) then
    call dis_bound_rs(rs1,rs2,svec)
  else
    svec(:) = 0.0
  end if
!
  if (isfeg.EQV..false.) then
!
    if ((use_TABUL.EQV..false.).AND.(use_POLAR.EQV..false.)) return

    if (use_POLAR.EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      if (max(rs1,rs2).eq.disulf(min(rs1,rs2))) then
        ll = 1
        lk = crlk_idx(rs1)
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            kk = at(rs2)%pol(j)
            if (ll.le.crosslink(lk)%nrspol) then 
              if ((ii.eq.crosslink(lk)%exclpol(ll,1)).AND.(kk.eq.crosslink(lk)%exclpol(ll,2))) then
                if (crosslink(lk)%is14pol(ll).EQV..false.) then
                  ll = ll + 1
                  cycle
                end if
              end if
            end if
            k = k + 1
            idx2(k,1) = ii
            idx2(k,2) = kk
            termpre(k) = 1.0
            if (ll.le.crosslink(lk)%nrspol) then 
              if ((ii.eq.crosslink(lk)%exclpol(ll,1)).AND.(kk.eq.crosslink(lk)%exclpol(ll,2))) then
                termpre(k) = fudge_el_14
                ll = ll + 1
              end if
            end if
          end do
        end do
      else if (rs1.eq.rs2) then
        do i=1,nrpolintra(rs1)
          k = k + 1
          idx2(k,1:2) = iaa(rs1)%polin(i,1:2)
          termpre(k) = fudge(rs1)%elin(i)
        end do
      else if (abs(rs1-rs2).le.1) then
        rlo = min(rs1,rs2)
        if (rlo.eq.rs2) svec(:) = -svec(:)
        do i=1,nrpolnb(rlo)
          k = k + 1
          idx2(k,1:2) = iaa(rlo)%polnb(i,1:2)
          termpre(k) = fudge(rlo)%elnb(i)
        end do
      else
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            k = k + 1
            idx2(k,1) = ii
            idx2(k,2) = at(rs2)%pol(j)
            termpre(k) = 1.0
          end do
        end do
      end if
      alcsz = k
      pfac = electric*scale_POLAR
      if (alcsz.gt.0) then
        if (use_IMPSOLV.EQV..true.) then
          if (scrq_model.le.2) then
            do k=1,alcsz
              dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
              dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
              dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
              termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
              terms(k) = scrq(idx2(k,1))*scrq(idx2(k,2))
              tsdk(k,1:3) = scrq_dr(:,idx2(k,2))*scrq(idx2(k,1))
              tsdi(k,1:3) = scrq_dr(:,idx2(k,1))*scrq(idx2(k,2))
              termqk(k) = termq(k)*scrq(idx2(k,1))
              termqi(k) = termq(k)*scrq(idx2(k,2))
            end do
          else if (scrq_model.eq.4) then
            pfac2 = par_IMPSOLV(8)*electric*scale_POLAR
            do k=1,alcsz
              dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
              dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
              dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
              termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
              terms(k) = atr(idx2(k,1))+atr(idx2(k,2))
            end do
          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do k=1,alcsz
              dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
              dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
              dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
              termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
              terms(k) = scrq(idx2(k,1))*scrq(idx2(k,2))
              tsdk(k,1:3) = scrq_dr(:,idx2(k,2))*scrq(idx2(k,1))
              tsdi(k,1:3) = scrq_dr(:,idx2(k,1))*scrq(idx2(k,2))
              termqk(k) = termq(k)*scrq(idx2(k,1))
              termqi(k) = termq(k)*scrq(idx2(k,2))
              termr(k) = atr(idx2(k,1))+atr(idx2(k,2))
            end do
          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            do k=1,alcsz
              dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
              dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
              dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
              termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
              terms1(k) = scrq(idx2(k,1))
              terms2(k) = scrq(idx2(k,2))
              tsdk(k,1:3) = scrq_dr(:,idx2(k,2))
              tsdi(k,1:3) = scrq_dr(:,idx2(k,1))
            end do
          else ! 7,8
            pfac2 = par_IMPSOLV(8)*electric
            pfac3 = 1./par_IMPSOLV(1)
            pfac4 = par_IMPSOLV(9)*pfac3
            do k=1,alcsz
              dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
              dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
              dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
              termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
              terms1(k) = scrq(idx2(k,1))
              terms2(k) = scrq(idx2(k,2))
              tsdk(k,1:3) = scrq_dr(:,idx2(k,2))
              tsdi(k,1:3) = scrq_dr(:,idx2(k,1))
              termr(k) = atr(idx2(k,1))+atr(idx2(k,2))
            end do
          end if
        else
          do k=1,alcsz
            dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
            dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
            dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
            termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
          end do
        end if
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        id2(1:alcsz) = id1(1:alcsz)*id1(1:alcsz)
!       finish up -> screening model specificity is high
        if (use_IMPSOLV.EQV..true.) then
          k = 0
          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
            if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz)
            termqk(1:alcsz) = termqk(1:alcsz)*term0(1:alcsz)
            termqi(1:alcsz) = termqi(1:alcsz)*term0(1:alcsz)
            term0(1:alcsz) = term0(1:alcsz)*termq(1:alcsz)
            for_tmpi(1:alcsz,1) = term0(1:alcsz)*tsdi(1:alcsz,1)
            for_tmpi(1:alcsz,2) = term0(1:alcsz)*tsdi(1:alcsz,2)
            for_tmpi(1:alcsz,3) = term0(1:alcsz)*tsdi(1:alcsz,3)
            for_tmpk(1:alcsz,1) = term0(1:alcsz)*tsdk(1:alcsz,1)
            for_tmpk(1:alcsz,2) = term0(1:alcsz)*tsdk(1:alcsz,2)
            for_tmpk(1:alcsz,3) = term0(1:alcsz)*tsdk(1:alcsz,3)
            terms1(1:alcsz) = term0(1:alcsz)*terms(1:alcsz)
            evec(6) = evec(6) + sum(terms1(1:alcsz))
            terms2(1:alcsz) = terms1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = terms2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = terms2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = terms2(1:alcsz)*dvec(1:alcsz,3)
            do k=1,alcsz
              sum_s(idx2(k,1)) = sum_s(idx2(k,1)) + termqi(k)
              sum_s(idx2(k,2)) = sum_s(idx2(k,2)) + termqk(k)
              ca_f(:,idx2(k,1)) = ca_f(:,idx2(k,1)) - termf(k,1:3) + for_tmpi(k,1:3)
              ca_f(:,idx2(k,2)) = ca_f(:,idx2(k,2)) + termf(k,1:3) + for_tmpk(k,1:3)
            end do
          else if (scrq_model.eq.4) then
            do k=1,alcsz
              if (d1(k).ge.terms(k)) then
                terms1(k) = id1(k)
                term0(k) = 2.0*id1(k)
              else
                terms1(k) = 1.0/terms(k)
                term0(k) = terms1(k)
              end if
            end do
            terms2(1:alcsz) = pfac2*termq(1:alcsz)*id1(1:alcsz)
            termqi(1:alcsz) = terms2(1:alcsz)*id2(1:alcsz)*term0(1:alcsz)
            terms2(1:alcsz) = terms2(1:alcsz)*terms1(1:alcsz)
            termf(1:alcsz,1) = termqi(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = termqi(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = termqi(1:alcsz)*dvec(1:alcsz,3)
            evec(6) = evec(6) + sum(terms2(1:alcsz))
            do k=1,alcsz
              ca_f(:,idx2(k,1)) = ca_f(:,idx2(k,1)) - termf(k,1:3)
              ca_f(:,idx2(k,2)) = ca_f(:,idx2(k,2)) + termf(k,1:3)
            end do
          else
            if ((scrq_model.eq.7).OR.(scrq_model.eq.8)) then
              call Vgenmu(terms1(1:alcsz),terms2(1:alcsz),i_sqm,terms(1:alcsz),alcsz)
              call Vgenmu_dr(terms1(1:alcsz),terms2(1:alcsz),i_sqm,termqi(1:alcsz),termqk(1:alcsz),alcsz)
              tsdk(1:alcsz,1) = tsdk(1:alcsz,1)*termqk(1:alcsz)
              tsdi(1:alcsz,1) = tsdi(1:alcsz,1)*termqi(1:alcsz)
              tsdk(1:alcsz,2) = tsdk(1:alcsz,2)*termqk(1:alcsz)
              tsdi(1:alcsz,2) = tsdi(1:alcsz,2)*termqi(1:alcsz)
              tsdk(1:alcsz,3) = tsdk(1:alcsz,3)*termqk(1:alcsz)
              tsdi(1:alcsz,3) = tsdi(1:alcsz,3)*termqi(1:alcsz)
              termqk(1:alcsz) = termq(1:alcsz)*termqk(1:alcsz)
              termqi(1:alcsz) = termq(1:alcsz)*termqi(1:alcsz)
            end if
            term0(1:alcsz) = pfac*id1(1:alcsz) ! cc4 without atq-terms
            terms1(1:alcsz) = electric*termq(1:alcsz)*terms(1:alcsz) ! ff1
            terms2(1:alcsz) = pfac2*termq(1:alcsz)/termr(1:alcsz) ! ff2
            for_tmpi(1:alcsz,1:3) = 0.0
            do k=1,alcsz
              if (abs(terms1(k)).gt.abs(terms2(k))) then
                term4(k) = terms1(k)
                term3(k) = 1.0
              else
                dd1 = (d1(k) - termr(k))*pfac3
                dd2 = 1.0 - dd1
                if (dd1.lt.0.0) then
                  term3(k) = (1.0 - par_IMPSOLV(9))
                  term4(k) = term3(k)*terms1(k) + par_IMPSOLV(9)*terms2(k)
                  term0(k) = term3(k)*term0(k)
                else if (dd2.gt.0.0) then
                  dum1 = par_IMPSOLV(9)*dd2
                  term3(k) = (1.0 - dum1)
                  term4(k) = term3(k)*terms1(k) + dum1*terms2(k)
                  term0(k) = term3(k)*term0(k)
                  dum1 = -scale_POLAR*id2(k)*(terms2(k)-terms1(k))*pfac4
                  for_tmpi(k,1) = dum1*dvec(k,1)
                  for_tmpi(k,2) = dum1*dvec(k,2)
                  for_tmpi(k,3) = dum1*dvec(k,3)
                else
                  term4(k) = terms1(k)
                  term3(k) = 1.0
                end if
              end if
              sum_s(idx2(k,1)) = sum_s(idx2(k,1)) + term0(k)*termqi(k)
              sum_s(idx2(k,2)) = sum_s(idx2(k,2)) + term0(k)*termqk(k)
            end do
            term0(1:alcsz) = pfac*term3(1:alcsz)*termq(1:alcsz) ! overwrite of cc4
            for_tmpk(1:alcsz,1) = -for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdk(1:alcsz,1)*id1(1:alcsz)
            for_tmpk(1:alcsz,2) = -for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdk(1:alcsz,2)*id1(1:alcsz)
            for_tmpk(1:alcsz,3) = -for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdk(1:alcsz,3)*id1(1:alcsz)
            for_tmpi(1:alcsz,1) = for_tmpi(1:alcsz,1) + term0(1:alcsz)*tsdi(1:alcsz,1)*id1(1:alcsz)
            for_tmpi(1:alcsz,2) = for_tmpi(1:alcsz,2) + term0(1:alcsz)*tsdi(1:alcsz,2)*id1(1:alcsz)
            for_tmpi(1:alcsz,3) = for_tmpi(1:alcsz,3) + term0(1:alcsz)*tsdi(1:alcsz,3)*id1(1:alcsz)
            terms1(1:alcsz) = scale_POLAR*term4(1:alcsz)*id1(1:alcsz)
            evec(6) = evec(6) + sum(terms1(1:alcsz))
            terms2(1:alcsz) = terms1(1:alcsz)*id2(1:alcsz)
            termf(1:alcsz,1) = terms2(1:alcsz)*dvec(1:alcsz,1)
            termf(1:alcsz,2) = terms2(1:alcsz)*dvec(1:alcsz,2)
            termf(1:alcsz,3) = terms2(1:alcsz)*dvec(1:alcsz,3)
            do k=1,alcsz
              ca_f(:,idx2(k,1)) = ca_f(:,idx2(k,1)) - termf(k,1:3) + for_tmpi(k,1:3)
              ca_f(:,idx2(k,2)) = ca_f(:,idx2(k,2)) + termf(k,1:3) + for_tmpk(k,1:3)
            end do
          end if
        else
          terms1(1:alcsz) = pfac*termq(1:alcsz)*id1(1:alcsz)
          evec(6) = evec(6) + sum(terms1(1:alcsz))
          terms2(1:alcsz) = terms1(1:alcsz)*id2(1:alcsz)
          do k=1,alcsz
            ca_f(:,idx2(k,1)) = ca_f(:,idx2(k,1)) - terms2(k)*dvec(k,1:3)
            ca_f(:,idx2(k,2)) = ca_f(:,idx2(k,2)) + terms2(k)*dvec(k,1:3)
          end do
        end if
      end if
!
    end if
!
    if (use_TABUL.EQV..true.) then
!
      if (tbp%rsmat(rs1,rs2).gt.0) then
        rhi = max(rs1,rs2)
        rlo = min(rs1,rs2)
        alcsz = tbp%rsmat(rhi,rlo) - tbp%rsmat(rlo,rhi) + 1
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
          kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
          k = k + 1
          dvec(k,1) = x(kk) - x(ii)
          dvec(k,2) = y(kk) - y(ii)
          dvec(k,3) = z(kk) - z(ii)
        end do
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        pfac = 1.0/tbp%res
        term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
        terms1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
        terms2(1:alcsz) = terms1(1:alcsz)*terms1(1:alcsz) ! quadratic
        termq(1:alcsz) = terms2(1:alcsz)*terms1(1:alcsz) ! cubic
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        pfac2 = pfac*scale_TABUL
!        id1(1:alcsz) = 1.0/d1(1:alcsz)
!        tbin(1:alcsz) = floor((d1(1:alcsz)-tbp%dis(1))/tbp%res) + 1
!        pfac = scale_TABUL*(1.0/tbp%res)
        k = 0
        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
          k = k + 1
          if (tbin(k).ge.tbp%bins) then
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
          else if (tbin(k).ge.1) then
            ii = atmol(molofrs(rlo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
            kk = atmol(molofrs(rhi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)

            hlp1 = 2.0*termq(k) - 3.0*terms2(k)
            hlp2 = termq(k) - terms2(k)
            hlp3 = hlp2 - terms2(k) + terms1(k)
            dhlp1 = 6.0*terms2(k) - 6.0*terms1(k)
            dhlp2 = 3.0*terms2(k) - 2.0*terms1(k)
            dhlp3 = dhlp2 - 2.0*terms1(k) + 1.0

            hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))

!            terms1(k) = (tbp%dis(tbin(k)+1) - d1(k))*tbp%pot(tbp%lst(i,3),tbin(k))
!            terms2(k) = (d1(k) - tbp%dis(tbin(k)))*tbp%pot(tbp%lst(i,3),tbin(k)+1)
!            termq(k) = pfac*(terms1(k)+terms2(k))
            evec(9) = evec(9) + scale_TABUL*hlp4
            hlp4 = pfac2*id1(k)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k)) - tbp%pot(tbp%lst(i,3),tbin(k)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1)))
!            term4(k) = id1(k)*pfac*(tbp%pot(tbp%lst(i,3),tbin(k)+1)-tbp%pot(tbp%lst(i,3),tbin(k)))

            ca_f(:,ii) = ca_f(:,ii) + hlp4*dvec(k,1:3)
            ca_f(:,kk) = ca_f(:,kk) - hlp4*dvec(k,1:3)
         else
            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
          end if
        end do
      end if
    end if
!
  else
!
    if (use_FEGS(6).EQV..true.) then
!
!     first set the necessary atom-specific parameters (mimimal)
      k = 0
      pfac = electric*par_FEG2(9)
!
      if (max(rs1,rs2).eq.disulf(min(rs1,rs2))) then
        ll = 1
        lk = crlk_idx(rs1)
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            kk = at(rs2)%pol(j)
            if (ll.le.crosslink(lk)%nrspol) then 
              if ((ii.eq.crosslink(lk)%exclpol(ll,1)).AND.(kk.eq.crosslink(lk)%exclpol(ll,2))) then
                if (crosslink(lk)%is14pol(ll).EQV..false.) then
                  ll = ll + 1
                  cycle
                end if
              end if
            end if
            k = k + 1
            idx2(k,1) = ii
            idx2(k,2) = kk
            termpre(k) = 1.0
            if (ll.le.crosslink(lk)%nrspol) then 
              if ((ii.eq.crosslink(lk)%exclpol(ll,1)).AND.(kk.eq.crosslink(lk)%exclpol(ll,2))) then
                termpre(k) = fudge_el_14
                ll = ll + 1
              end if
            end if
          end do
        end do
      else if (rs1.eq.rs2) then
        do i=1,nrpolintra(rs1)
          k = k + 1
          idx2(k,1:2) = iaa(rs1)%polin(i,1:2)
          termpre(k) = fudge(rs1)%elin(i)
        end do
      else if (abs(rs1-rs2).le.1) then
        rlo = min(rs1,rs2)
        if (rlo.eq.rs2) svec(:) = -svec(:)
        do i=1,nrpolnb(rlo)
          k = k + 1
          idx2(k,1:2) = iaa(rlo)%polnb(i,1:2)
          termpre(k) = fudge(rlo)%elnb(i)
        end do
      else
        do i=1,at(rs1)%npol
          ii = at(rs1)%pol(i)
          do j=1,at(rs2)%npol
            k = k + 1
            idx2(k,1) = ii
            idx2(k,2) = at(rs2)%pol(j)
            termpre(k) = 1.0
          end do
        end do
      end if
      alcsz = k
      if (alcsz.gt.0) then
        do k=1,alcsz
          dvec(k,1) = x(idx2(k,2)) - x(idx2(k,1))
          dvec(k,2) = y(idx2(k,2)) - y(idx2(k,1))
          dvec(k,3) = z(idx2(k,2)) - z(idx2(k,1))
          termq(k) = termpre(k)*atq(idx2(k,1))*atq(idx2(k,2))
        end do
!       now perform the vectorizable bulk operations (maximal)
        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + &
 &                    dvec(1:alcsz,3)**2
        d1(1:alcsz) = sqrt(d2(1:alcsz))
        id1(1:alcsz) = 1.0/d1(1:alcsz)
        termqi(1:alcsz) = 1.0/(d1(1:alcsz) + par_FEG2(10))
        terms1(1:alcsz) = pfac*termq(1:alcsz)*termqi(1:alcsz)
        evec(6) = evec(6) + sum(terms1(1:alcsz))
        termqk(1:alcsz) = terms1(1:alcsz)*termqi(1:alcsz)*id1(1:alcsz)
        do k=1,alcsz
          ca_f(:,idx2(k,1)) = ca_f(:,idx2(k,1)) - termqk(k)*dvec(k,1:3)
          ca_f(:,idx2(k,2)) = ca_f(:,idx2(k,2)) + termqk(k)*dvec(k,1:3)
        end do
      end if
    end if
!       
  end if
!
!
end
!
!-----------------------------------------------------------------------
!
! this routine exploits that some data can be precalculated on a per-atom (rather than per atom pair) basis
! r1: first charge, transformed to take mean
! r1b: first charge, inner derivative
! r2/r2b: the same as vectors for the other charges
!
! this fxn should not be called with i_genmu being 1 (redundant)
!
subroutine Vgenmu_all3(r1,r1b,r2,r2b,ro,dr1,dr2,alcsz)
!
  use iounit, ONLY: ilog
  use energies, ONLY: i_sqm
!
  implicit none
!
  integer,INTENT(IN):: alcsz
  RTYPE, INTENT(IN):: r1,r1b,r2(alcsz),r2b(alcsz)
  RTYPE, INTENT(OUT):: dr1(alcsz),dr2(alcsz),ro(alcsz)
!
  RTYPE help(alcsz),dum1
!
  if ((r1.le.0.0).AND.(i_sqm.ne.1)) then
    write(ilog,*) 'Fatal. Called Vgenmu_all3(...) with corrupted arguments.'
    call fexit()
  end if
!
  if (i_sqm.eq.0) then ! cheap
    ro(1:alcsz) = r1*r2(1:alcsz)         ! geometric mean
    help(1:alcsz) = 0.5*r1b*r2b(1:alcsz) ! half of inverse of ro
    dr1(1:alcsz) = help(1:alcsz)*r2(1:alcsz)*r2(1:alcsz)
    dr2(1:alcsz) = help(1:alcsz)*r1*r1
  else if (i_sqm.eq.-1) then ! one division
    ro(1:alcsz) = 1.0/(r1+r2(1:alcsz))   ! denominator always > 0.
    help(1:alcsz) = -ro(1:alcsz)*ro(1:alcsz)
    dr1(1:alcsz) = help(1:alcsz)*r1b
    dr2(1:alcsz) = help(1:alcsz)*r2b(1:alcsz)
  else if (i_sqm.eq.-2) then ! one inverse square root
    ro(1:alcsz) = sqrt(1.0/(r1+r2(1:alcsz))) ! denominator always > 0.
    help(1:alcsz) = -0.5*ro(1:alcsz)*ro(1:alcsz)*ro(1:alcsz)
    dr1(1:alcsz) = help(1:alcsz)*r1b
    dr2(1:alcsz) = help(1:alcsz)*r2b(1:alcsz)
  else if (i_sqm.eq.2) then  ! one square root and one division
    ro(1:alcsz) = sqrt(r1+r2(1:alcsz))
    help(1:alcsz) = 0.5/ro(1:alcsz)
    dr1(1:alcsz) = help(1:alcsz)*r1b
    dr2(1:alcsz) = help(1:alcsz)*r2b(1:alcsz)
  else if (i_sqm.eq.1) then  ! nothing to do here really
    ro(1:alcsz) = r1 + r2(1:alcsz)
    dr1(1:alcsz) = r1b
    dr2(1:alcsz) = r2b(1:alcsz)
  else                       ! calls to pow with FP exponents -> very expensive
    dum1 = (1.0/(1.0*i_sqm))
    ro(1:alcsz) = (r1+r2(1:alcsz))**dum1
    help(1:alcsz) = dum1*ro(1:alcsz)/(r1+r2(1:alcsz))
    dr1(1:alcsz) = help(1:alcsz)*r1b
    dr2(1:alcsz) = help(1:alcsz)*r2b(1:alcsz)
  end if
!
end
!
!---------------------------------------------------------------------------------------------
!
! these are the models with purely environmental screening
! note that group consistency is handled entirely in force_setup_scrqs2(...) in force.f90
! the same is true for preprocessing for models 5/6
!
subroutine Vforce_PSCRM1256_C(evec,rs1,ca_f,sum_s,lora,hira3,hira2,which)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,lo,k,i,shii,jj,ki,kj,imol,chnksz,chnk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n),sum_s(n)
!
  RTYPE term6(VCI_SIZE),xt0(VCI_SIZE),dvec(VCI_SIZE,3),termsqd(VCI_SIZE,3)
  RTYPE foin(VCI_SIZE,3),xt1(VCI_SIZE),id1(VCI_SIZE),xt2(VCI_SIZE),xt3(VCI_SIZE),xt5(VCI_SIZE),xt4(VCI_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3),sum_k(hira2)
  RTYPE ixis(9,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  integer ka(hira2)
  RTYPE term0(hira2),termx(hira2),termx2(hira2),termq(3,hira2)
  integer sh(hira3-lora+1),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE lop1,lop3,lop4
  RTYPE incri
!  
! probably seg-faulted already if hira is 0
  hira=hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  sum_k(:) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
    if (scrq_model.le.2) then
      ixis(5,j) = scrq(ii)
      ixis(6:8,j) = scrq_dr(:,ii)
    else
      ixis(5:6,j) = atsavinfo(1:2,ii)
      ixis(7:9,j) = scrq_dr(:,ii)
    end if
  end do
  ii = refat(rs1)
!
  if ((which.eq.1).OR.(which.eq.3)) then ! standard range, recompute shift vectors
    if (which.eq.1) then
      k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    else
      k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira3)
    end if
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else if (which.eq.4) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
! reuse parameters
  lop3 = electric*scale_POLAR
!
! this loop has data dependency and a conditional -> as few operations as possible in here
  lo = 1
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
!
  if (scrq_model.le.2) then
    do k=1,hira
      do j=k1(k),k1(k)+sh(k)
        jj = ka(k) + j - k1(k)
        dveci(jj,1) = srav(k,1) + atinfo(1,j)
        dveci(jj,2) = srav(k,2) + atinfo(2,j)
        dveci(jj,3) = srav(k,3) + atinfo(3,j)
        term0(jj) = atinfo(4,j)
        termx(jj) = scrq(j)
        termq(:,jj) = scrq_dr(:,j)
      end do
    end do
    do chnk=1,hira2,VCI_SIZE
      hi = min(chnk + VCI_SIZEM1,hira2)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      termsqd(1:chnksz,:) = transpose(termq(:,chnk:hi))
      termsqd(1:chnksz,1) = termsqd(1:chnksz,1)*term0(chnk:hi)
      termsqd(1:chnksz,2) = termsqd(1:chnksz,2)*term0(chnk:hi)
      termsqd(1:chnksz,3) = termsqd(1:chnksz,3)*term0(chnk:hi)
      xt0(1:chnksz) = term0(chnk:hi)
      xt2(1:chnksz) = termx(chnk:hi)
      do i=1,at(rs1)%na
        if (ixis(4,i).eq.0.0) cycle
        lop1 = lop3*ixis(4,i)
        lop4 = lop1*ixis(5,i)
!       inverse distance squared
        xt1(1:chnksz) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
     &                       (dvec(1:chnksz,2)-ixis(2,i))**2 + &
     &                       (dvec(1:chnksz,3)-ixis(3,i))**2)
        id1(1:chnksz) = sqrt(xt1(1:chnksz))
        term6(1:chnksz) = id1(1:chnksz)*lop1*xt0(1:chnksz)
        sum_k(chnk:hi) = sum_k(chnk:hi) + ixis(5,i)*term6(1:chnksz)
        term6(1:chnksz) = term6(1:chnksz)*xt2(1:chnksz)
        incri = sum(term6(1:chnksz))
        sum_s(k0(i,2)) = sum_s(k0(i,2)) + incri
        evec(6) = evec(6) + ixis(5,i)*sum(term6(1:chnksz))
!       force 
        xt1(1:chnksz) = term6(1:chnksz)*xt1(1:chnksz)*ixis(5,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*xt1(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*xt1(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*xt1(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1)) + ixis(6,i)*incri
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2)) + ixis(7,i)*incri
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3)) + ixis(8,i)*incri
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1) + lop4*termsqd(1:chnksz,1)*id1(1:chnksz)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2) + lop4*termsqd(1:chnksz,2)*id1(1:chnksz)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3) + lop4*termsqd(1:chnksz,3)*id1(1:chnksz)
      end do
    end do
  else
    do k=1,hira
      do j=k1(k),k1(k)+sh(k)
        jj = ka(k) + j - k1(k)
        dveci(jj,1) = srav(k,1) + atinfo(1,j)
        dveci(jj,2) = srav(k,2) + atinfo(2,j)
        dveci(jj,3) = srav(k,3) + atinfo(3,j)
        term0(jj) = atinfo(4,j)
        termx(jj) = atsavinfo(1,j)
        termx2(jj) = atsavinfo(2,j)
        termq(:,jj) = scrq_dr(:,j)
      end do
    end do
    do chnk=1,hira2,VCI_SIZE
      hi = min(chnk + VCI_SIZEM1,hira2)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      termsqd(1:chnksz,:) = transpose(termq(:,chnk:hi))
      xt0(1:chnksz) = term0(chnk:hi)
      xt2(1:chnksz) = termx(chnk:hi)
      xt3(1:chnksz) = termx2(chnk:hi)
      do i=1,at(rs1)%na
        if (ixis(4,i).eq.0.0) cycle
        lop1 = lop3*ixis(4,i)
        xt1(1:chnksz) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
     &                       (dvec(1:chnksz,2)-ixis(2,i))**2 + &
     &                       (dvec(1:chnksz,3)-ixis(3,i))**2)
        id1(1:chnksz) = sqrt(xt1(1:chnksz))
!       screened Coulomb potential
        if (i_sqm.eq.1) then
          term6(1:chnksz) = xt2(1:chnksz) + ixis(5,i)
          xt4(1:chnksz) = 0.5
          xt5(1:chnksz) = 0.5
        else
          call Vgenmu_all3(ixis(5,i),ixis(6,i),xt2(1:chnksz),xt3(1:chnksz),term6(1:chnksz),xt5(1:chnksz),xt4(1:chnksz),chnksz)
        end if
        id1(1:chnksz) = xt0(1:chnksz)*id1(1:chnksz)*lop1 ! term6(1:chnksz)
        incri = sum(id1(1:chnksz)*xt5(1:chnksz))
        sum_s(k0(i,2)) = sum_s(k0(i,2)) + incri
        xt4(1:chnksz) = id1(1:chnksz)*xt4(1:chnksz)
        sum_k(chnk:hi) = sum_k(chnk:hi) + xt4(1:chnksz)
        term6(1:chnksz) = term6(1:chnksz)*id1(1:chnksz)
        evec(6) = evec(6) + sum(term6(1:chnksz))
        xt5(1:chnksz) = xt1(1:chnksz)*term6(1:chnksz)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*xt5(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*xt5(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*xt5(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1)) + ixis(7,i)*incri
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2)) + ixis(8,i)*incri
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3)) + ixis(9,i)*incri
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1) + termsqd(1:chnksz,1)*xt4(1:chnksz)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2) + termsqd(1:chnksz,2)*xt4(1:chnksz)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3) + termsqd(1:chnksz,3)*xt4(1:chnksz)
      end do
    end do
  end if
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
    sum_s(k1(k):ki) = sum_s(k1(k):ki) + sum_k(ka(k):kj)
  end do
!
end
!
!-------------------------------------------------------------------
!
! the fxn for purely distance-dependent screening is notably simpler
!
subroutine Vforce_PSCRM4_C(evec,rs1,ca_f,lora,hira3,hira2,which)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,k,i,shii,jj,ki,kj,imol,chnksz,chnk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VCI_SIZE),xt0(VCI_SIZE),dvec(VCI_SIZE,3),id2(VCI_SIZE)
  RTYPE foin(VCI_SIZE,3),xt1(VCI_SIZE),xt2(VCI_SIZE),xt3(VCI_SIZE)
!
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(5,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  integer ka(hira2)
  RTYPE term0(hira2),termx(hira2)
  integer sh(hira3-lora+1),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE lop1,lop3
!  
! probably seg-faulted already if hira is 0
  hira = hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
    ixis(5,j) = atr(ii)
  end do
  ii = refat(rs1)
!
  if ((which.eq.1).OR.(which.eq.3)) then ! standard range, recompute shift vectors
    if (which.eq.1) then
      k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    else
      k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira3)
    end if
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else if (which.eq.4) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
!
  lop3 = par_IMPSOLV(8)*electric*scale_POLAR
  do k=1,hira
    do j=k1(k),k1(k)+sh(k)
      jj = ka(k) + j - k1(k)
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
      term0(jj) = atinfo(4,j)
      termx(jj) = atr(j)
    end do
  end do
  do chnk=1,hira2,VCI_SIZE
    hi = min(chnk + VCI_SIZEM1,hira2)
    chnksz = hi-chnk+1
    dvec(1:chnksz,1) = dveci(chnk:hi,1)
    dvec(1:chnksz,2) = dveci(chnk:hi,2)
    dvec(1:chnksz,3) = dveci(chnk:hi,3)
    xt0(1:chnksz) = term0(chnk:hi)
    xt2(1:chnksz) = termx(chnk:hi)
    do i=1,at(rs1)%na
      if (ixis(4,i).eq.0.0) cycle
      lop1 = lop3*ixis(4,i)
!     distance squared
      xt1(1:chnksz) = ((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                   (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                   (dvec(1:chnksz,3)-ixis(3,i))**2)
      id2(1:chnksz) = 1.0/xt1(1:chnksz)
      term6(1:chnksz) = (ixis(5,i)+xt2(1:chnksz))**2
      xt3(1:chnksz) = lop1*xt0(1:chnksz)*id2(1:chnksz)
      do k=1,chnksz
        if (xt1(k).lt.term6(k)) then
          xt3(k) = xt3(k)*sqrt(xt1(k)/term6(k))
          id2(k) = 0.5*id2(k)
        end if
      end do
      evec(6) = evec(6) + sum(xt3(1:chnksz))
!     force 
      xt1(1:chnksz) = 2.0*xt3(1:chnksz)*id2(1:chnksz)
      foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*xt1(1:chnksz)
      foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*xt1(1:chnksz)
      foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*xt1(1:chnksz)
      ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
      ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
      ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
      for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
      for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
      for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
    end do
  end do
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
  end do
!
end
!
!---------------------------------------------------------------------------------------------
!
subroutine Vforce_PSCRM3789_C(evec,rs1,ca_f,sum_s,lora,hira3,hira2,which)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,lo,k,i,shii,jj,ki,kj,imol,chnksz,chnk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n),sum_s(n)
!
  RTYPE termg(VCI_SIZE),terme(VCI_SIZE),xt0(VCI_SIZE),dvec(VCI_SIZE,3),termsqd(VCI_SIZE,3)
  RTYPE foin(VCI_SIZE,3),xt1(VCI_SIZE),id1(VCI_SIZE),d1(VCI_SIZE),xt2(VCI_SIZE),xt3(VCI_SIZE),xt5(VCI_SIZE),xt4(VCI_SIZE)
  RTYPE xt6(VCI_SIZE),xt7(VCI_SIZE),xt8(VCI_SIZE),xt9(VCI_SIZE),xt10(VCI_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3),sum_k(hira2)
  RTYPE ixis(10,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  integer ka(hira2)
  RTYPE term0(hira2),termx(hira2),termx2(hira2),termq(3,hira2),termx3(hira2)
  integer sh(hira3-lora+1),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8,lop10
  RTYPE incri
!  
! probably seg-faulted already if hira is 0
  hira = hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
  sum_k(:) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
    if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
      ixis(5,j) = scrq(ii)
      ixis(6:8,j) = scrq_dr(:,ii)
      ixis(9,j) = atr(ii)
    else
      ixis(5:6,j) = atsavinfo(1:2,ii)
      ixis(7:9,j) = scrq_dr(:,ii)
      ixis(10,j) = atr(ii)
    end if
  end do
  ii = refat(rs1)
!
  if ((which.eq.1).OR.(which.eq.3)) then ! standard range, recompute shift vectors
    if (which.eq.1) then
      k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    else
      k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira3)
    end if
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else if (which.eq.4) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
! reuse parameters
  lop2 = par_IMPSOLV(8)
  lop3 = electric*scale_POLAR
  lop5 = 1.0/par_IMPSOLV(1)
  lop6 = par_IMPSOLV(9)*lop5
  lop7 = par_IMPSOLV(9)
  lop10 = 1.0 - par_IMPSOLV(9)
!
! this loop has data dependency and a conditional -> as few operations as possible in here
  lo = 1
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
!
  if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
    do k=1,hira
      do j=k1(k),k1(k)+sh(k)
        jj = ka(k) + j - k1(k)
        dveci(jj,1) = srav(k,1) + atinfo(1,j)
        dveci(jj,2) = srav(k,2) + atinfo(2,j)
        dveci(jj,3) = srav(k,3) + atinfo(3,j)
        term0(jj) = atinfo(4,j)
        termx(jj) = scrq(j)
        termq(:,jj) = scrq_dr(:,j)
        termx2(jj) = atr(j)
      end do
    end do
    do chnk=1,hira2,VCI_SIZE
      hi = min(chnk + VCI_SIZEM1,hira2)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      termsqd(1:chnksz,:) = transpose(termq(:,chnk:hi))
!      termsqd(1:chnksz,1) = termsqd(1:chnksz,1)*term0(chnk:hi)
!      termsqd(1:chnksz,2) = termsqd(1:chnksz,2)*term0(chnk:hi)
!      termsqd(1:chnksz,3) = termsqd(1:chnksz,3)*term0(chnk:hi)
      xt0(1:chnksz) = term0(chnk:hi)
      xt2(1:chnksz) = termx(chnk:hi)
      xt4(1:chnksz) = termx2(chnk:hi)
      do i=1,at(rs1)%na
        if (ixis(4,i).eq.0.0) cycle
        lop1 = lop3*ixis(4,i)
        lop4 = lop1*ixis(5,i)
        xt3(1:chnksz) = ((dvec(1:chnksz,1)-ixis(1,i))**2 + &
     &                   (dvec(1:chnksz,2)-ixis(2,i))**2 + &
     &                   (dvec(1:chnksz,3)-ixis(3,i))**2)
        d1(1:chnksz) = sqrt(xt3(1:chnksz))
        id1(1:chnksz) = 1.0/d1(1:chnksz)
        xt1(1:chnksz) = id1(1:chnksz)*id1(1:chnksz)
        termg(1:chnksz) = id1(1:chnksz)*lop1*xt0(1:chnksz) ! gas-phase Coulomb
        xt3(1:chnksz) = lop5*(d1(1:chnksz) - ixis(9,i) - xt4(1:chnksz))
        xt5(1:chnksz) = 1.0 - xt3(1:chnksz)
        xt6(1:chnksz) = 1.0
        terme(1:chnksz) = ixis(5,i)*xt2(1:chnksz) ! environmental screening factor 
        xt8(1:chnksz) = 0.0
        do k=1,chnksz
          if (xt5(k).gt.0.0) then ! within first shell or closer
            lop8 = lop2/(ixis(9,i)+xt4(k))
            if (terme(k).gt.lop8) then
!             do nothing
            else    
              if (xt3(k).lt.0.0) then ! within contact radius
                xt6(k) = lop10
                terme(k) = terme(k)*lop10 + lop7*lop8
              else if (xt5(k).gt.0.0) then ! within first shell
                xt8(k) = lop6*id1(k)*(lop8-terme(k))
                xt6(k) = 1.0 - lop7*xt5(k)
                terme(k) = xt6(k)*terme(k) + lop7*xt5(k)*lop8
              end if
            end if
          end if
        end do
        xt5(1:chnksz) = termg(1:chnksz)*xt6(1:chnksz)
        sum_k(chnk:hi) = sum_k(chnk:hi) + ixis(5,i)*xt5(1:chnksz)
        incri = sum(xt5(1:chnksz)*xt2(1:chnksz))
        sum_s(k0(i,2)) = sum_s(k0(i,2)) + incri
        xt8(1:chnksz) = termg(1:chnksz)*xt8(1:chnksz)
        termg(1:chnksz) = termg(1:chnksz)*terme(1:chnksz)
        evec(6) = evec(6) + sum(termg(1:chnksz))
!       force 
        xt1(1:chnksz) = termg(1:chnksz)*xt1(1:chnksz) + xt8(1:chnksz)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*xt1(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*xt1(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*xt1(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1)) + ixis(6,i)*incri
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2)) + ixis(7,i)*incri
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3)) + ixis(8,i)*incri
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1) + xt5(1:chnksz)*termsqd(1:chnksz,1)*ixis(5,i)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2) + xt5(1:chnksz)*termsqd(1:chnksz,2)*ixis(5,i)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3) + xt5(1:chnksz)*termsqd(1:chnksz,3)*ixis(5,i)
      end do
    end do
  else
    do k=1,hira
      do j=k1(k),k1(k)+sh(k)
        jj = ka(k) + j - k1(k)
        dveci(jj,1) = srav(k,1) + atinfo(1,j)
        dveci(jj,2) = srav(k,2) + atinfo(2,j)
        dveci(jj,3) = srav(k,3) + atinfo(3,j)
        term0(jj) = atinfo(4,j)
        termx(jj) = atsavinfo(1,j)
        termx2(jj) = atsavinfo(2,j)
        termq(:,jj) = scrq_dr(:,j)
        termx3(jj) = atr(j)
      end do
    end do
    do chnk=1,hira2,VCI_SIZE
      hi = min(chnk + VCI_SIZEM1,hira2)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      termsqd(1:chnksz,:) = transpose(termq(:,chnk:hi))
      xt0(1:chnksz) = term0(chnk:hi)
      xt2(1:chnksz) = termx(chnk:hi)
      xt3(1:chnksz) = termx2(chnk:hi)
      xt7(1:chnksz) = termx3(chnk:hi)
      do i=1,at(rs1)%na
        if (ixis(4,i).eq.0.0) cycle
        lop1 = lop3*ixis(4,i)
        xt4(1:chnksz) = ((dvec(1:chnksz,1)-ixis(1,i))**2 + &
     &                   (dvec(1:chnksz,2)-ixis(2,i))**2 + &
     &                   (dvec(1:chnksz,3)-ixis(3,i))**2)
        d1(1:chnksz) = sqrt(xt4(1:chnksz))
        id1(1:chnksz) = 1.0/d1(1:chnksz)
        xt1(1:chnksz) = id1(1:chnksz)*id1(1:chnksz)
!       screened Coulomb potential
        if (i_sqm.eq.1) then
          terme(1:chnksz) = xt2(1:chnksz) + ixis(5,i)
          xt5(1:chnksz) = 0.5
          xt4(1:chnksz) = 0.5
        else
          call Vgenmu_all3(ixis(5,i),ixis(6,i),xt2(1:chnksz),xt3(1:chnksz),terme(1:chnksz),xt5(1:chnksz),xt4(1:chnksz),chnksz)
        end if
        termg(1:chnksz) = id1(1:chnksz)*lop1*xt0(1:chnksz) ! gas-phase Coulomb
        xt9(1:chnksz) = lop5*(d1(1:chnksz) - ixis(10,i) - xt7(1:chnksz))
        xt10(1:chnksz) = 1.0 - xt9(1:chnksz)
        xt6(1:chnksz) = 1.0
        xt8(1:chnksz) = 0.0
        do k=1,chnksz
          if (xt10(k).gt.0.0) then ! within first shell or closer
            lop8 = lop2/(ixis(10,i)+xt7(k))
            if (terme(k).gt.lop8) then
!             do nothing
            else    
              if (xt9(k).lt.0.0) then ! within contact radius
                xt6(k) = lop10
                terme(k) = terme(k)*lop10 + lop7*lop8
              else if (xt10(k).gt.0.0) then ! within first shell
                xt8(k) = lop6*id1(k)*(lop8-terme(k))
                xt6(k) = 1.0 - lop7*xt10(k)
                terme(k) = xt6(k)*terme(k) + lop7*xt10(k)*lop8
              end if
            end if
          end if
        end do
        xt9(1:chnksz) = termg(1:chnksz)*xt6(1:chnksz)
        sum_k(chnk:hi) = sum_k(chnk:hi) + xt9(1:chnksz)*xt4(1:chnksz)
        incri = sum(xt9(1:chnksz)*xt5(1:chnksz))
        sum_s(k0(i,2)) = sum_s(k0(i,2)) + incri
        xt8(1:chnksz) = termg(1:chnksz)*xt8(1:chnksz)
        termg(1:chnksz) = termg(1:chnksz)*terme(1:chnksz)
        evec(6) = evec(6) + sum(termg(1:chnksz))
        xt1(1:chnksz) = termg(1:chnksz)*xt1(1:chnksz) + xt8(1:chnksz)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*xt1(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*xt1(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*xt1(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1)) + ixis(7,i)*incri
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2)) + ixis(8,i)*incri
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3)) + ixis(9,i)*incri
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1) + termsqd(1:chnksz,1)*xt9(1:chnksz)*xt4(1:chnksz)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2) + termsqd(1:chnksz,2)*xt9(1:chnksz)*xt4(1:chnksz)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3) + termsqd(1:chnksz,3)*xt9(1:chnksz)*xt4(1:chnksz)
      end do
    end do
  end if
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
    sum_s(k1(k):ki) = sum_s(k1(k):ki) + sum_k(ka(k):kj)
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_dsavdr(sdr,tpi)
!
  use cutoffs
  use forces
  use atoms, ONLY: n
#ifdef ENABLE_THREADS
  use threads
  use iounit
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer i,l,sta,sto
  RTYPE sdr(3,n)
!
#ifdef ENABLE_THREADS
  integer j,tpx,ii
  integer(KIND=8) ttimer
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, Vforce_dsavdr(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    if (thr_dlb(15,1).gt.0) then
      if (tpi.eq.1) thr_dlb(15,2) = thr_dlb(15,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(33,tpi) = thr_timings(33,tpi) + ttimer
    end if
    sta = thr_limits(85,tpi)
    sto = thr_limits(86,tpi)
    do i=sta,sto
      ii = sisa(i)%nix + sum(thr_sisa(i,1:thrdat%maxn)%nix)
      if (ii.gt.sisa(i)%alsz) then
        l = ii-sisa(i)%alsz
        call Vsisa_resize(i,sisa,l)
      end if
      do tpx=1,thrdat%maxn
        if (thr_sisa(i,tpx)%nix.gt.0) then
          j = sisa(i)%nix + 1
          sisa(i)%nix = sisa(i)%nix + thr_sisa(i,tpx)%nix
          sisa(i)%dr(:,j:sisa(i)%nix) = thr_sisa(i,tpx)%dr(:,1:thr_sisa(i,tpx)%nix)
          sisa(i)%ix(j:sisa(i)%nix) = thr_sisa(i,tpx)%ix(1:thr_sisa(i,tpx)%nix)
        end if
      end do
    end do
  else
    sta = 1
    sto = n
  end if
#else
  sta = 1
  sto = n
#endif
  do i=sta,sto
    do l=1,sisa(i)%nix
      sdr(:,i) = sdr(:,i) - sisa(i)%dr(:,l)
    end do
  end do
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (thr_dlb(15,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(34,tpi) = thr_timings(34,tpi) + ttimer
    end if
  end if
#endif
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_dcbdscrq(rs,ca_f)
!
  use polypep
  use energies
  use forces
  use atoms
!
  implicit none
!
  integer, INTENT(IN):: rs
!
  integer i,ii,ll,l
  RTYPE ca_f(3,n),totss
!
  do i=1,at(rs)%npol
    ii = at(rs)%pol(i)
    totss = sum_scrcbs(ii)+sum_scrcbs_tr(ii)+sum_scrcbs_lr(ii) 
    do l=1,sisa(ii)%nix
      ll = sisa(ii)%ix(l)
      ca_f(:,ll) = ca_f(:,ll) + totss*sisa(ii)%qr(:,l)
    end do
  end do
!
end
!
!
!-------------------------------------------------------------------
!
subroutine Vsisa_resize(ii,si,inc)
!
  use atoms, ONLY: n
  use forces
  use energies, ONLY: use_POLAR
!
  implicit none
!
  integer, INTENT(IN):: ii,inc
!
  type(t_sisa) si(n)
!
  integer sibi(max(si(ii)%nix,1)),oldn
  RTYPE sibu(3,max(si(ii)%nix,1))
!
  oldn = si(ii)%alsz
  if (si(ii)%nix.gt.0) then ! negative values not recommended
    sibi(:) = si(ii)%ix(1:si(ii)%nix)
    sibu(:,:) = si(ii)%dr(:,1:si(ii)%nix)
  end if
  if (inc.ne.0) then !
    si(ii)%alsz = max(si(ii)%nix,si(ii)%alsz + inc)
  else
    si(ii)%alsz = si(ii)%nix ! reduce to minimum size
  end if
  if (oldn.gt.0) deallocate(si(ii)%dr)
  if (si(ii)%alsz.gt.0) then
    allocate(si(ii)%dr(3,si(ii)%alsz))
    if (si(ii)%nix.gt.0) si(ii)%dr(:,1:si(ii)%nix) = sibu(:,:)
  end if
  if (oldn.gt.0) deallocate(si(ii)%ix)
  if (si(ii)%alsz.gt.0) then
    allocate(si(ii)%ix(si(ii)%alsz))
    if (si(ii)%nix.gt.0) si(ii)%ix(1:si(ii)%nix) = sibi(:)
  end if
  if (use_POLAR.EQV..true.) then
    if (si(ii)%nix.gt.0) sibu(:,:) = si(ii)%qr(:,1:si(ii)%nix)
    if (oldn.gt.0) deallocate(si(ii)%qr)
    if (si(ii)%alsz.gt.0) then
      allocate(si(ii)%qr(3,si(ii)%alsz))
      if (si(ii)%nix.gt.0) si(ii)%qr(:,1:si(ii)%nix) = sibu(:,:)
    end if
  end if
!
end
!
!--------------------------------------------------------------------
!
