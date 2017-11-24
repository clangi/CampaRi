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
!------------------------------------------------------------------------------
!
! the vectorized version of en_rsp:
! note that for thread-safe execution we assume that multiple instances
! of sv_in(:) and evec(:) are passed
! note that the speed-up hinges crucially on the arrays being statically 
! allocated at call-time (for some reason dynamic allocation is terribly slow),
! and on defining the range -> 1:alcsz
! it is expected that this function is called based on a neighbor list or similar
! (it is not necessarily appropriate to obtain the intended cutoffs simply calling
! this function for every residue pair)
!
subroutine Ven_rsp(evec,rs1,rs2,sv_in,cut)
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
  use fyoc
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  RTYPE evec(MAXENERGYTERMS),sv_in(n)
!
  integer rs,i,j,k,ii,kk,imol1,imol2,alcsz,lk,ll
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE tmaxd(at(rs1)%na*at(rs2)%na)
  integer idx2(at(rs1)%na*at(rs2)%na,2)
  RTYPE incr,svec(3),efvoli,efvolk,datri,datrk,pfac,rt23,pwc2,olaps(3)
  logical ismember,doid2,docut
  integer olap(at(rs1)%na*at(rs2)%na)
  RTYPE d2diff(at(rs1)%na*at(rs2)%na)
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
        call Ven_rsp_feg(evec,rs1,rs2,cut)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Ven_rsp_feg(evec,rs1,rs2,cut)
        return
      end if
    end if
  end if
!
  doid2 = .false.
  if  ((use_attLJ.EQV..true.).OR.&
 &     ((use_IPP.EQV..true.).AND.(use_hardsphere.EQV..false.))) then
    doid2 = .true.
  end if
  docut = cut
  if (mccutmode.eq.2) docut = .false.
!
! to set up, consider four different cases: intra-residue (rs1 == rs2), neighbors in seq., crosslinked neighbors, or others
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
      if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
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
  else if ((abs(rs1-rs2).eq.1).AND.(max(rs2,rs1).ne.disulf(min(rs1,rs2)))) then
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
      if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
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
    if (max(rs2,rs1).eq.disulf(min(rs1,rs2))) then
!     crosslink requires adjustment of LJ parameters and exclusions
      call dis_bound_rs(rs1,rs2,svec) ! shift vector may be nonzero for intermolecular crosslinks in PBC
      k = 1
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
          dvec(k,1) = x(kk) - x(ii) + svec(1) 
          dvec(k,2) = y(kk) - y(ii) + svec(2)
          dvec(k,3) = z(kk) - z(ii) + svec(3)
          d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
          if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
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
    else
!     there is no topological relationship for remaining residues -> always check BC
      call dis_bound_rs(rs1,rs2,svec)
!     loop over complete set of atom-atom interactions in distant residues
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
          if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
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
  end if
!
! now do the actual potentials
  k = k - 1
  alcsz = k
  if (alcsz.gt.0) then
    if (doid2.EQV..true.) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
    end if
    if (use_IPP.EQV..true.) then
      if (use_hardsphere.EQV..true.) then
        d2diff(1:alcsz) = max((terms(1:alcsz)-d2(1:alcsz))/terms(1:alcsz),0.0d0)
        olap(1:alcsz) = ceiling(d2diff(1:alcsz))
        term1(1:alcsz) = dble(olap(1:alcsz))
        evec(1) = evec(1) + screenbarrier*sum(term1(1:alcsz))
      else
        term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**nhalf
        evec(1) = evec(1) + scale_IPP*4.0*&
 &                sum(term0(1:alcsz)*term1(1:alcsz))
      end if
    end if
    if (use_attLJ.EQV..true.) then
      term1(1:alcsz) = (terms(1:alcsz)*id2(1:alcsz))**3
      evec(3) = evec(3) - scale_attLJ*4.0*&
 &                sum(term0(1:alcsz)*term1(1:alcsz))
    end if
    if (use_IMPSOLV.EQV..true.) then
      d1(1:alcsz) = sqrt(d2(1:alcsz))
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
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
!            if ((ii.eq.714).AND.(kk.eq.736)) write(0,*) d1(k),-incr*efvolk
            sv_in(ii) = sv_in(ii) + incr*efvolk
            if (d1(k).gt.(tmaxd(k)-datri)) then
              incr = -(tmaxd(k)-d1(k))/datri
            else if (d1(k).lt.termr(k)) then
              incr = -d1(k)/termr(k)
            else
              incr = -1.0
            end if
            sv_in(kk) = sv_in(kk) + incr*efvoli
!            if ((ii.eq.714).AND.(kk.eq.736)) write(0,*) d1(k),-incr*efvoli
          end if
        end do
      else if (par_IMPSOLV2(1).eq.2) then
        do k=1,alcsz
          if (d1(k).lt.tmaxd(k)) then
            ii = idx2(k,1)
            kk = idx2(k,2)
            call sphere_overlaps(atr(ii),atr(kk),d1(k),par_IMPSOLV(1),atvol(ii),atvol(kk),olaps(1:3))
!              if ((ii.eq.714).AND.(kk.eq.736)) write(0,*) d1(k),(olaps(2)-olaps(1))*atsavred(kk),(olaps(3)-olaps(1))*atsavred(ii)
            if (olaps(2).gt.0.0) then
              sv_in(ii) = sv_in(ii) - (olaps(2)-olaps(1))*atsavred(kk)
            end if
            if (olaps(3).gt.0.0) then
              sv_in(kk) = sv_in(kk) - (olaps(3)-olaps(1))*atsavred(ii)
            end if
          end if
        end do
      end if
    end if
    if (use_WCA.EQV..true.) then
      pfac = scale_WCA*par_WCA(2)
      rt23 = ROOT26*ROOT26
      pwc2 = par_WCA(1)*par_WCA(1)
      do k=1,alcsz
        if (term0(k).gt.0.0) then
          if (d2(k).lt.rt23*terms(k)) then
            term1(k) = (terms(k)/d2(k))**3
            term1(k) = 4.0*(term1(k)**2 - term1(k)) + 1.0
            evec(5) = evec(5) + scale_WCA*term0(k)*&
 &                             (term1(k) - par_WCA(2))
          else if (d2(k).lt.pwc2*terms(k)) then
            term1(k) = par_WCA(3)*d2(k)/terms(k) + par_WCA(4)
            term1(k) = 0.5*cos(term1(k)) - 0.5
            evec(5) = evec(5) + pfac*term0(k)*term1(k)
          end if
        end if
      end do
    end if
  end if
!
  if (rs1.eq.rs2) then
!   add correction terms
    if (use_CORR.EQV..true.) then
      call e_corrector(rs,evec)
    end if
  end if
!
!
!
end

!
!-----------------------------------------------------------------------
!
! the same as Ven_rsp with support for ghosted interactions
!
subroutine Ven_rsp_feg(evec,rs1,rs2,cut)
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
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  integer rs,i,j,k,ii,kk,alcsz,lk,ll
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE id2(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE terms(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  RTYPE evec(MAXENERGYTERMS),svec(3)
  logical docut
! 
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
! note that if fegmode is 1, only one of the two is ghosted, otherwise background Hamiltonian
! if fegmode is 2 and both are ghosted, we ensure the fully de-coupled ghost always remains
! fully ghosted so as to not pick up spurious interactions with increased coupling
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) &
 &  return
!
  docut = cut
  if (mccutmode.eq.2) docut = .false.
!
! three different cases: intra-residue (rs1 == rs2), neighbors in seq., or others
!
  if (rs1.eq.rs2) then
!
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
      if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
        term0(k) = fudge(rs)%rsin(i)*fudge(rs)%rsin_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if (term0(k).gt.0.0) then
          terms(k) = fudge(rs)%rsin_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        evec(1) = evec(1) + par_FEG2(1)*4.0*&
 &      sum(term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz))
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        evec(3) = evec(3) - par_FEG2(5)*4.0*&
 &                  sum(term0(1:alcsz)*term2(1:alcsz))
      end if
!
    end if
!
!   add correction terms
    if (use_CORR.EQV..true.) then
      call e_corrector(rs,evec)
    end if
!
!
!
  else if ((abs(rs1-rs2).eq.1).AND.(max(rs2,rs1).ne.disulf(min(rs1,rs2)))) then
!
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
      if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
        term0(k) = fudge(rs)%rsnb(i)*fudge(rs)%rsnb_lje(i) !lj_eps(attyp(ii),attyp(kk))
        if (term0(k).gt.0.0) then
          terms(k) = fudge(rs)%rsnb_ljs(i) !lj_sig(attyp(ii),attyp(kk))
          k = k + 1
        end if
      end if
    end do
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        evec(1) = evec(1) + par_FEG2(1)*4.0*&
 &      sum(term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz))
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        evec(3) = evec(3) - par_FEG2(5)*4.0*&
 &                  sum(term0(1:alcsz)*term2(1:alcsz))
      end if
!
    end if
!
!
!
  else
!
    if (max(rs2,rs1).eq.disulf(min(rs1,rs2))) then
!
      call dis_bound_rs(rs1,rs2,svec) ! shift vector may be nonzero for intermolecular crosslinks in PBC
      k = 1
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
          dvec(k,1) = x(kk) - x(ii) + svec(1)
          dvec(k,2) = y(kk) - y(ii) + svec(2)
          dvec(k,3) = z(kk) - z(ii) + svec(3)
          d2(k) = dvec(k,1)**2 + dvec(k,2)**2 + dvec(k,3)**2
          if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
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
              k = k + 1
            end if
          else if (ll.le.crosslink(lk)%nrsin) then
            if ((ii.eq.crosslink(lk)%exclin(ll,1)).AND.(kk.eq.crosslink(lk)%exclin(ll,2))) ll = ll + 1
          end if
        end do
      end do
!
    else
  !   there is no topological relationship for remaining residues -> always check BC
      call dis_bound_rs(rs1,rs2,svec)
!     loop over complete set of atom-atom interactions in distant residues
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
          if ((d2(k).lt.mcnb_cutoff2).OR.(docut.EQV..false.)) then
            term0(k) = lj_eps(attyp(ii),attyp(kk))
            if (term0(k).gt.0.0) then
              terms(k) = lj_sig(attyp(ii),attyp(kk))
              k = k + 1
            end if
          end if
        end do
      end do
    end if
    k = k - 1
    alcsz = k
    if (alcsz.gt.0) then
      id2(1:alcsz) = 1.0/d2(1:alcsz)
      term1(1:alcsz) = (d2(1:alcsz)/terms(1:alcsz))**3 
      if (use_FEGS(1).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(2))
        evec(1) = evec(1) + par_FEG2(1)*4.0*&
 &      sum(term0(1:alcsz)*term2(1:alcsz)*term2(1:alcsz))
      end if
      if (use_FEGS(3).EQV..true.) then
        term2(1:alcsz) = 1.0/(term1(1:alcsz) + par_FEG2(6))
        evec(3) = evec(3) - par_FEG2(5)*4.0*&
 &                  sum(term0(1:alcsz)*term2(1:alcsz))
      end if
!
    end if
!
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the vectorized version of en_rsp_long
! note that the speed-up hinges crucially on the arrays being statically 
! allocated at call-time (for some reason dynamic allocation is terribly slow),
! and on defining the range -> 1:alcsz
!
subroutine Ven_rsp_long(evec,rs1,rs2,cut)
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
  use fyoc
  use grandensembles
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  integer rs,i,j,k,ii,kk,ll,lk,rhi,rlo,imol1,imol2,alcsz,rsl,rsh
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,hlp1,hlp2,hlp3
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  integer tbin(at(rs1)%na*at(rs2)%na),idx2(at(rs1)%na*at(rs2)%na,2)
  logical ismember
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
        call Ven_rsp_long_feg(evec,rs1,rs2,cut)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Ven_rsp_long_feg(evec,rs1,rs2,cut)
        return
      end if
    end if
  end if
!
! first get image shift vector
  svec(:) = 0.0
  rs = 0
  if (rs1.eq.rs2) then
!   do nothing
  else if ((abs(rs1-rs2).eq.1).AND.(max(rs1,rs2).ne.disulf(min(rs1,rs2)))) then
    rs = min(rs1,rs2)
    if (rs1.gt.rs2) then
      call dis_bound_rs(rs2,rs1,svec)
    else
      call dis_bound_rs(rs1,rs2,svec)
    end if
  else if (max(rs1,rs2).eq.disulf(min(rs1,rs2))) then
    call dis_bound_rs(min(rs1,rs2),max(rs1,rs2),svec)
  else
    call dis_bound_rs(rs1,rs2,svec)
  end if
!
!
  if (use_POLAR.EQV..true.) then
    k = 0
    pfac = electric*scale_POLAR
    if (rs1.eq.rs2) then
      k = nrpolintra(rs1)
      do i=1,nrpolintra(rs1)
        idx2(i,1) = iaa(rs1)%polin(i,1)
        idx2(i,2) = iaa(rs1)%polin(i,2)
        term0(i) = fudge(rs1)%elin(i)
      end do
    else if ((abs(rs1-rs2).eq.1).AND.(max(rs1,rs2).ne.disulf(min(rs1,rs2)))) then
      k = nrpolnb(rs)
      do i=1,nrpolnb(rs)
        idx2(i,1) = iaa(rs)%polnb(i,1)
        idx2(i,2) = iaa(rs)%polnb(i,2)
        term0(i) = fudge(rs)%elnb(i)
      end do
    else if (max(rs1,rs2).eq.disulf(min(rs1,rs2))) then
      rsl = min(rs1,rs2)
      rsh = max(rs1,rs2)
      ll = 1
      lk = crlk_idx(rsl)
      do i=1,at(rsl)%npol
        ii = at(rsl)%pol(i)
        do j=1,at(rsh)%npol
          kk = at(rsh)%pol(j)
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
          term0(k) = 1.0
          if (ll.le.crosslink(lk)%nrspol) then 
            if ((ii.eq.crosslink(lk)%exclpol(ll,1)).AND.(kk.eq.crosslink(lk)%exclpol(ll,2))) then
              term0(k) = fudge_el_14
              ll = ll + 1
            end if
          end if
        end do
      end do
    else
      do i=1,at(rs1)%npol
        ii = at(rs1)%pol(i)
        do j=1,at(rs2)%npol
          k = k + 1
          idx2(k,1) = ii
          idx2(k,2) = at(rs2)%pol(j)
          term0(k) = 1.0
        end do
      end do
    end if
!
    dvec(1:k,1) = x(idx2(1:k,2)) - x(idx2(1:k,1)) + svec(1)
    dvec(1:k,2) = y(idx2(1:k,2)) - y(idx2(1:k,1)) + svec(2)
    dvec(1:k,3) = z(idx2(1:k,2)) - z(idx2(1:k,1)) + svec(3)

    if (use_IMPSOLV.EQV..true.) then
      if (((scrq_model.ge.7).AND.(scrq_model.le.9)).OR.((scrq_model.ge.3).AND.(scrq_model.le.4))) then
        termr(1:k) = atr(idx2(1:k,1)) + atr(idx2(1:k,2))
        term2(1:k) = atq(idx2(1:k,1))*atq(idx2(1:k,2))
      end if
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        do i=1,k
          term1(i) = atq(idx2(i,1))*atq(idx2(i,2))*genmu(scrq(idx2(i,1)),scrq(idx2(i,2)),i_sqm)
        end do 
      end if
      if (((scrq_model.ge.1).AND.(scrq_model.le.3)).OR.(scrq_model.eq.9)) then
        term1(1:k) = atq(idx2(1:k,1))*atq(idx2(1:k,2))*scrq(idx2(1:k,1))*scrq(idx2(1:k,2))
      end if
    else
      term0(1:k) = term0(1:k)*atq(idx2(1:k,1))*atq(idx2(1:k,2))
    end if
    d2(1:k) = dvec(1:k,1)**2 + dvec(1:k,2)**2 + dvec(1:k,3)**2
    d1(1:k) = sqrt(d2(1:k))
    id1(1:k) = 1.0/d1(1:k)
    if (use_IMPSOLV.EQV..true.) then
      if (((scrq_model.ge.7).AND.(scrq_model.le.9)).OR.(scrq_model.eq.3)) then
        term2(1:k) = electric*par_IMPSOLV(8)*term2(1:k)/termr(1:k)
        term1(1:k) = electric*term1(1:k)
        hlp2 = 1.0/par_IMPSOLV(1)
        do i=1,k
          if (abs(term1(i)).gt.abs(term2(i))) then
            term1(i) = term0(i)*term1(i)
          else
!           dissociated
            if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
              term1(i) = term0(i)*term1(i)
!           close contact
            else if (d1(i).lt.termr(i)) then
              term1(i) = term0(i)*(par_IMPSOLV(9)*(term2(i)-term1(i)) + term1(i))
!           first-shell regime
            else
              hlp1 = (termr(i)+par_IMPSOLV(1)-d1(i))*hlp2
              term1(i) = term0(i)*((1.0-par_IMPSOLV(9)*hlp1)*term1(i) + par_IMPSOLV(9)*hlp1*term2(i))
            end if
          end if
        end do
        evec(6) = evec(6) + scale_POLAR*sum(term1(1:k)*id1(1:k))
      else if (scrq_model.eq.4) then
        do i=1,k
          if (d1(i).ge.termr(i)) then
            term0(i) = term0(i)*term2(i)*id1(i)
          else
            term0(i) = term0(i)*term2(i)/termr(i)
          end if
        end do
        evec(6) = evec(6) + pfac*par_IMPSOLV(8)*sum(term0(1:k)*id1(1:k))
      else if (((scrq_model.ge.1).AND.(scrq_model.le.2)).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
        evec(6) = evec(6) + pfac*sum(term0(1:k)*term1(1:k)*id1(1:k))
      end if
    else
      evec(6) = evec(6) + pfac*sum(term0(1:k)*id1(1:k))
    end if
  end if
!
! note that this might be terribly inefficient if POLAR and TABUL are used concurrently
  if (use_TABUL.EQV..true.) then
!
    if (tbp%rsmat(rs1,rs2).gt.0) then
      if ((rs1.ne.rs2).AND.(rs.eq.rs2)) svec(:) = -1.0*svec(:)
      rhi = max(rs1,rs2)
      rlo = min(rs1,rs2)
      k = 0
      if (rs1.eq.rs2) then
        alcsz = tbp%rsvec(rlo)-tbp%rsmat(rlo,rhi)+1
      else
        alcsz = tbp%rsmat(rhi,rlo)-tbp%rsmat(rlo,rhi)+1
      end if
      do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rlo,rhi)+alcsz-1
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
      termr(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
      tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(termr(1:alcsz))))
      term0(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
      term1(1:alcsz) = term0(1:alcsz)*term0(1:alcsz) ! quadratic
      term2(1:alcsz) = term1(1:alcsz)*term0(1:alcsz) ! cubic
      k = 0
      do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rlo,rhi)+alcsz-1
        k = k + 1
        if (tbin(k).ge.tbp%bins) then
          evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
        else if (tbin(k).ge.1) then
          hlp1 = 2.0*term2(k) - 3.0*term1(k)
          hlp2 = term2(k) - term1(k)
          hlp3 = hlp2 - term1(k) + term0(k)
          termr(k) = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
 &                 tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
          evec(9) = evec(9) + scale_TABUL*termr(k)
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
!  else
!
!!   there is no topological relationship for remaining residues -> always check BC
!    call dis_bound_rs(rs1,rs2,svec)
!
!    if (use_POLAR.EQV..true.) then
!
!!     first set the necessary atom-specific parameters (mimimal)
!      k = 0
!      pfac = electric*scale_POLAR
!      alcsz = at(rs1)%npol*at(rs2)%npol
!      if (alcsz.gt.0) then
!        if (use_IMPSOLV.EQV..true.) then
!          if (scrq_model.le.2) then
!            do g1=1,at(rs1)%ndpgrps
!              do g2=1,at(rs2)%ndpgrps
!              if ((at(rs1)%dpgrp(g1)%nc.ne.0).AND.(at(rs2)%dpgrp(g2)%nc.ne.0)) then
!                hlp3 = 0.0
!!              else if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) then
!!                hlp3 = 0.0
!              else
!                hlp3 = at(rs1)%dpgrp(g1)%tsavq*at(rs2)%dpgrp(g2)%tsavq
!              end if
!              do i=1,at(rs1)%dpgrp(g1)%nats
!              ii = at(rs1)%dpgrp(g1)%ats(i)
!              do j=1,at(rs2)%dpgrp(g2)%nats
!                kk = at(rs2)%dpgrp(g2)%ats(j)
!
!
!!            do i=1,at(rs1)%npol
!!              ii = at(rs1)%pol(i)
!!              do j=1,at(rs2)%npol
!!                kk = at(rs2)%pol(j)
!                k = k + 1
!                dvec(k,1) = x(kk) - x(ii)! + svec(1)
!                dvec(k,2) = y(kk) - y(ii)! + svec(2)
!                dvec(k,3) = z(kk) - z(ii)! + svec(3)
!                term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk) - hlp3
!              end do
!              end do
!            end do
!            end do
!          else if (scrq_model.eq.4) then
!            do i=1,at(rs1)%npol
!              ii = at(rs1)%pol(i)
!              do j=1,at(rs2)%npol
!                kk = at(rs2)%pol(j)
!                k = k + 1
!                dvec(k,1) = x(kk) - x(ii)! + svec(1)
!                dvec(k,2) = y(kk) - y(ii)! + svec(2)
!                dvec(k,3) = z(kk) - z(ii)! + svec(3)
!                term0(k) = atq(ii)*atq(kk)
!                term2(k) = 1.0/(atr(ii)+atr(kk))
!              end do
!            end do
!          else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
!            pfac2 = par_IMPSOLV(8)*electric
!            do i=1,at(rs1)%npol
!              ii = at(rs1)%pol(i)
!              do j=1,at(rs2)%npol
!                kk = at(rs2)%pol(j)
!                k = k + 1
!                dvec(k,1) = x(kk) - x(ii)! + svec(1)
!                dvec(k,2) = y(kk) - y(ii)! + svec(2)
!                dvec(k,3) = z(kk) - z(ii)! + svec(3)
!                term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
!                termr(k) = atr(ii)+atr(kk)
!                term2(k) = atq(ii)*atq(kk)/termr(k)
!              end do
!            end do
!          else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
!            pfac2 = par_IMPSOLV(8)*electric
!            do i=1,at(rs1)%npol
!              ii = at(rs1)%pol(i)
!              do j=1,at(rs2)%npol
!                kk = at(rs2)%pol(j)
!                k = k + 1
!                dvec(k,1) = x(kk) - x(ii)! + svec(1)
!                dvec(k,2) = y(kk) - y(ii)! + svec(2)
!                dvec(k,3) = z(kk) - z(ii)! + svec(3)
!                term0(k) = atq(ii)*atq(kk)*&
! &                  genmu(scrq(ii),scrq(kk),i_sqm)
!              end do
!            end do
!          else ! 7,8
!            pfac2 = par_IMPSOLV(8)*electric
!            do i=1,at(rs1)%npol
!              ii = at(rs1)%pol(i)
!              do j=1,at(rs2)%npol
!                kk = at(rs2)%pol(j)
!                k = k + 1
!                dvec(k,1) = x(kk) - x(ii)! + svec(1)
!                dvec(k,2) = y(kk) - y(ii)! + svec(2)
!                dvec(k,3) = z(kk) - z(ii)! + svec(3)
!                term0(k) = atq(ii)*atq(kk)*&
! &                    genmu(scrq(ii),scrq(kk),i_sqm)
!                termr(k) = atr(ii)+atr(kk)
!                term2(k) = atq(ii)*atq(kk)/termr(k)
!              end do
!            end do
!          end if
!        else
!          do i=1,at(rs1)%npol
!            ii = at(rs1)%pol(i)
!            do j=1,at(rs2)%npol
!              kk = at(rs2)%pol(j)
!              k = k + 1
!              dvec(k,1) = x(kk) - x(ii)! + svec(1)
!              dvec(k,2) = y(kk) - y(ii)! + svec(2)
!              dvec(k,3) = z(kk) - z(ii)! + svec(3)
!              term0(k) = atq(ii)*atq(kk)
!            end do
!          end do
!        end if
!!       now perform the vectorizable bulk operations (maximal)
!        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
!        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
!        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
!        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
!        d1(1:alcsz) = sqrt(d2(1:alcsz))
!        id1(1:alcsz) = 1.0/d1(1:alcsz)
!!       finish up -> screening model specificity is high
!        if (use_IMPSOLV.EQV..true.) then
!          k = 0
!          if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.&
! &            (scrq_model.eq.6)) then
!            term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
!            evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
!          else if (scrq_model.eq.4) then
!            do i=1,at(rs1)%npol
!              do j=1,at(rs2)%npol
!                k = k + 1
!                if (id1(k).le.term2(k)) then
!                  term1(k) = term0(k)*id1(k)
!                else
!                  term1(k) = term0(k)*term2(k) 
!                end if
!              end do
!            end do
!            evec(6) = evec(6) + par_IMPSOLV(8)*pfac*&
! &                       sum(id1(1:alcsz)*term1(1:alcsz))
!          else
!            term2(1:alcsz) = term2(1:alcsz)*pfac2
!            term0(1:alcsz) = term0(1:alcsz)*electric
!            do i=1,at(rs1)%npol
!              do j=1,at(rs2)%npol
!                k = k + 1
!                if (abs(term0(k)).gt.abs(term2(k))) then
!                  term1(k) = term0(k)
!                else
!!                 dissociated
!                  if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
!                    term1(k) = term0(k)
!!                 close contact
!                  else if (d1(k).lt.termr(k)) then
!                    term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + &
! &                                        term0(k)
!!                 first-shell regime
!                  else
!                    term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/&
! &                             par_IMPSOLV(1)
!                    term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k)&
! &                   + par_IMPSOLV(9)*term1(k)*term2(k) 
!                  end if
!                end if
!              end do
!            end do
!            term1(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)
!            evec(6) = evec(6) + scale_POLAR*sum(term1(1:alcsz))
!          end if
!        else
!          term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
!          evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
!        end if
!      end if
!
!    end if
!
!    if (use_TABUL.EQV..true.) then
!
!      if (tbp%rsmat(rs1,rs2).gt.0) then
!        rhi = max(rs1,rs2)
!        rlo = min(rs1,rs2)
!        alcsz = tbp%rsmat(rhi,rlo) - tbp%rsmat(rlo,rhi) + 1
!        k = 0
!        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
!          ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
!          kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
!          k = k + 1
!          dvec(k,1) = x(kk) - x(ii)! + svec(1)
!          dvec(k,2) = y(kk) - y(ii)! + svec(2)
!          dvec(k,3) = z(kk) - z(ii)! + svec(3)
!        end do
!        dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
!        dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
!        dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
!        d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
!        d1(1:alcsz) = sqrt(d2(1:alcsz))
!        pfac = 1.0/tbp%res
!        termr(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
!        tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(termr(1:alcsz))))
!        term0(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
!        term1(1:alcsz) = term0(1:alcsz)*term0(1:alcsz) ! quadratic
!        term2(1:alcsz) = term1(1:alcsz)*term0(1:alcsz) ! cubic
!        k = 0
!        do i=tbp%rsmat(rlo,rhi),tbp%rsmat(rhi,rlo)
!          k = k + 1
!          if (tbin(k).ge.tbp%bins) then
!            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
!          else if (tbin(k).ge.1) then
!            hlp1 = 2.0*term2(k) - 3.0*term1(k)
!            hlp2 = term2(k) - term1(k)
!            hlp3 = hlp2 - term1(k) + term0(k)
!            termr(k) = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k)+1) + &
! &                   tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k)+1))
!            evec(9) = evec(9) + scale_TABUL*termr(k)
!          else
!            evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
!          end if
!        end do
!      end if
!    end if
!!       
!  end if
!
!
end
!
!-----------------------------------------------------------------------
!
subroutine Ven_rsp_long_feg(evec,rs1,rs2,cut)
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
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
  logical, INTENT(IN):: cut
!
  integer rs,i,j,k,ii,kk,ll,lk
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,hlp1,hlp2
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE id1s(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  integer idx2(at(rs1)%na*at(rs2)%na,2)
!
! we have a potential override to cover in par_FEG3, which allows residues to be
! fully de-coupled at all times
! note that if fegmode is 1, only one of the two is ghosted, otherwise background Hamiltonian
! if fegmode is 2 and both are ghosted, we ensure the fully de-coupled ghost always remains
! fully ghosted so as to not pick up spurious interactions with increased coupling
  if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.))&
 &  return
!
! first get image shift vector
  svec(:) = 0.0
  rs = 0
  if (rs1.eq.rs2) then
!   do nothing
  else if ((abs(rs1-rs2).eq.1).AND.(max(rs1,rs2).ne.disulf(min(rs1,rs2)))) then
    rs = min(rs1,rs2)
    if (rs1.gt.rs2) then
      call dis_bound_rs(rs2,rs1,svec)
    else
      call dis_bound_rs(rs1,rs2,svec)
    end if
  else if (max(rs1,rs2).eq.disulf(min(rs1,rs2))) then
    call dis_bound_rs(rs1,rs2,svec)
  else
    call dis_bound_rs(rs1,rs2,svec)
  end if
!
!
  if (use_FEGS(6).EQV..true.) then
    k = 0
    pfac = electric*par_FEG2(9)
    if (rs1.eq.rs2) then
      k = nrpolintra(rs1)
      do i=1,nrpolintra(rs1)
        idx2(i,1) = iaa(rs1)%polin(i,1)
        idx2(i,2) = iaa(rs1)%polin(i,2)
        term0(i) = fudge(rs1)%elin(i)
      end do
    else if ((abs(rs1-rs2).eq.1).AND.(max(rs1,rs2).ne.disulf(min(rs1,rs2)))) then
      k = nrpolnb(rs)
      do i=1,nrpolnb(rs)
        idx2(i,1) = iaa(rs)%polnb(i,1)
        idx2(i,2) = iaa(rs)%polnb(i,2)
        term0(i) = fudge(rs)%elnb(i)
      end do
    else if (max(rs1,rs2).eq.disulf(min(rs1,rs2))) then
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
          term0(k) = 1.0
          if (ll.le.crosslink(lk)%nrspol) then 
            if ((ii.eq.crosslink(lk)%exclpol(ll,1)).AND.(kk.eq.crosslink(lk)%exclpol(ll,2))) then
              term0(k) = fudge_el_14
              ll = ll + 1
            end if
          end if
        end do
      end do
    else
      do i=1,at(rs1)%npol
        ii = at(rs1)%pol(i)
        do j=1,at(rs2)%npol
          k = k + 1
          idx2(k,1) = ii
          idx2(k,2) = at(rs2)%pol(j)
          term0(k) = 1.0
        end do
      end do
    end if
!
    dvec(1:k,1) = x(idx2(1:k,2)) - x(idx2(1:k,1)) + svec(1)
    dvec(1:k,2) = y(idx2(1:k,2)) - y(idx2(1:k,1)) + svec(2)
    dvec(1:k,3) = z(idx2(1:k,2)) - z(idx2(1:k,1)) + svec(3)

    if (use_FEGS(4).EQV..true.) then
      if (((scrq_model.ge.7).AND.(scrq_model.le.9)).OR.((scrq_model.ge.3).AND.(scrq_model.le.4))) then
        termr(1:k) = atr(idx2(1:k,1)) + atr(idx2(1:k,2))
        term2(1:k) = atq(idx2(1:k,1))*atq(idx2(1:k,2))
      end if
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        do i=1,k
          term1(i) = atq(idx2(i,1))*atq(idx2(i,2))*genmu(scrq(idx2(i,1)),scrq(idx2(i,2)),i_sqm)
        end do 
      end if
      if (((scrq_model.ge.1).AND.(scrq_model.le.3)).OR.(scrq_model.eq.9)) then
        term1(1:k) = atq(idx2(1:k,1))*atq(idx2(1:k,2))*scrq(idx2(1:k,1))*scrq(idx2(1:k,2))
      end if
    else
      term0(1:k) = term0(1:k)*atq(idx2(1:k,1))*atq(idx2(1:k,2))
    end if
    d2(1:k) = dvec(1:k,1)**2 + dvec(1:k,2)**2 + dvec(1:k,3)**2
    d1(1:k) = sqrt(d2(1:k))
    id1s(1:k) = 1.0/(d1(1:k)+par_FEG2(10))
    if (use_FEGS(4).EQV..true.) then
      if (((scrq_model.ge.7).AND.(scrq_model.le.9)).OR.(scrq_model.eq.3)) then
        term2(1:k) = electric*par_IMPSOLV(8)*term2(1:k)/termr(1:k)
        term1(1:k) = electric*term1(1:k)
        hlp2 = 1.0/par_IMPSOLV(1)
        do i=1,k
          if (abs(term1(i)).gt.abs(term2(i))) then
            term1(i) = term0(i)*term1(i)
          else
!           dissociated
            if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
              term1(i) = term0(i)*term1(i)
!           close contact
            else if (d1(i).lt.termr(i)) then
              term1(i) = term0(i)*(par_IMPSOLV(9)*(term2(i)-term1(i)) + term1(i))
!           first-shell regime
            else
              hlp1 = (termr(i)+par_IMPSOLV(1)-d1(i))*hlp2
              term1(i) = term0(i)*((1.0-par_IMPSOLV(9)*hlp1)*term1(i) + par_IMPSOLV(9)*hlp1*term2(i))
            end if
          end if
        end do
        evec(6) = evec(6) + par_FEG2(9)*sum(term1(1:k)*id1s(1:k))
      else if (scrq_model.eq.4) then
        id1(1:k) = 1.0/d1(1:k)
        do i=1,k
          if (d1(i).ge.termr(i)) then
            term0(i) = term0(i)*term2(i)*id1(i)
          else
            term0(i) = term0(i)*term2(i)/termr(i)
          end if
        end do
        evec(6) = evec(6) + pfac*par_IMPSOLV(8)*sum(term0(1:k)*id1s(1:k))
      else if (((scrq_model.ge.1).AND.(scrq_model.le.2)).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
        evec(6) = evec(6) + pfac*sum(term0(1:k)*term1(1:k)*id1s(1:k))
      end if
    else
      evec(6) = evec(6) + pfac*sum(term0(1:k)*id1s(1:k))
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! note that this routine does NOT support crosslink corrections (we assume two crosslinked
! residues are never so far apart as to trigger the long-range corrections)
!
subroutine Ven_rsp_lrel(evec,rs1,rs2)
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
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer i,j,k,ii,kk,alcsz
  integer g1,g2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,pfac2,hlp3
  RTYPE d2(at(rs1)%na*at(rs2)%na)
  RTYPE dvec(at(rs1)%na*at(rs2)%na,3)
  RTYPE d1(at(rs1)%na*at(rs2)%na),id1(at(rs1)%na*at(rs2)%na)
  RTYPE term0(at(rs1)%na*at(rs2)%na)
  RTYPE termr(at(rs1)%na*at(rs2)%na)
  RTYPE term1(at(rs1)%na*at(rs2)%na)
  RTYPE term2(at(rs1)%na*at(rs2)%na)
  logical ismember
!
  if ((molofrs(rs1).eq.molofrs(rs2)).AND.(abs(rs1-rs2).le.1)) then
    write(ilog,*) 'Fatal. Called Ven_rsp_lrel(...) with residues that are topologically too close. This is either a bug or &
 &can possibly arise due to an abnormally long linkage.'
  end if
!
! Prevent non-present molecules from interacting with each other
! or with present molecules.  Note that intra-nonpresent-molecule
! interactions are still allowed, in order to correctly simulate
! an infinite-dilution implicit solvent reference state particle
! bath -> partially redundant overlap with use_FEG
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    ii = molofrs(rs1)
    kk = molofrs(rs2)
    if ((ii.ne.kk).AND.((.not.ismember(ispresent,ii)).OR.(.not.ismember(ispresent,kk)))) return
  end if
!
! branch out for ghosting
  if (use_FEG.EQV..true.) then
    if (fegmode.eq.1) then
   if(((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 & ((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..true.))) then
        call Ven_rsp_lrel_feg(evec,rs1,rs2)
        return
      end if
    else if (fegmode.eq.2) then
      if((par_FEG(rs1).EQV..true.).OR.(par_FEG(rs2).EQV..true.))then
        call Ven_rsp_lrel_feg(evec,rs1,rs2)
        return
      end if
    end if
  end if
!
  if (use_POLAR.EQV..true.) then
!
    if (rs1.gt.rs2) then
      call dis_bound_rs(rs2,rs1,svec)
    else
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
    pfac = electric*scale_POLAR
    if (lrel_mc.eq.1) then
      k = 0
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).AND.(at(rs2)%dpgrp(g2)%nc.eq.0)) cycle ! dipole-dipole
            if ((at(rs1)%dpgrp(g1)%nc.eq.0).OR.(at(rs2)%dpgrp(g2)%nc.eq.0)) then
              hlp3 = at(rs1)%dpgrp(g1)%tsavq*at(rs2)%dpgrp(g2)%tsavq
            else
              hlp3 = 0.0
            end if
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk) - hlp3
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
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
              term0(k) = atq(ii)*atq(kk)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
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
            term0(k) = atq(ii)*atq(kk)
          end do
          end do
        end do
        end do
      end if
      alcsz = k
    else if (lrel_mc.eq.2) then
      k = 0
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = atq(ii)*atq(kk) 
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = atq(ii)*atq(kk)
            end do
          end do
          end do
        end do
      end if
      alcsz = k
!   reduced monopole-monopole
    else if (lrel_mc.eq.3) then
      k = 0
      if (use_IMPSOLV.EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle ! dipole-dipole or dipole-monopole
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
          end do
        end do
      end if
      alcsz = k
    end if ! which lrel_mc
!
!   now perform the vectorizable bulk operations (maximal)
    dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
    dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
    dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
    d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
    d1(1:alcsz) = sqrt(d2(1:alcsz))
    id1(1:alcsz) = 1.0/d1(1:alcsz)
!
!   finish up -> screening model specificity is high
    if (use_IMPSOLV.EQV..true.) then
      k = 0
      if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
        term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
        evec(6) = evec(6) + pfac*sum(term1(1:alcsz))

      else if (scrq_model.eq.4) then
        do k=1,alcsz
          if (id1(k).le.term2(k)) then
            term1(k) = term0(k)*id1(k)
          else
            term1(k) = term0(k)*term2(k) 
          end if
        end do           
        evec(6) = evec(6) + par_IMPSOLV(8)*pfac*sum(id1(1:alcsz)*term1(1:alcsz))
      else   ! 3,7,8,9
        term2(1:alcsz) = term2(1:alcsz)*pfac2
        term0(1:alcsz) = term0(1:alcsz)*electric
        do k=1,alcsz
          if (abs(term0(k)).gt.abs(term2(k))) then
            term1(k) = term0(k)
          else
!        dissociated
            if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
              term1(k) = term0(k)
!           close contact
            else if (d1(k).lt.termr(k)) then
              term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + term0(k)
!           first-shell regime
            else
              term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/par_IMPSOLV(1)
              term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k) + par_IMPSOLV(9)*term1(k)*term2(k) 
            end if
          end if
        end do
        term1(1:alcsz) = term1(1:alcsz)*id1(1:alcsz)
        evec(6) = evec(6) + scale_POLAR*sum(term1(1:alcsz))
      end if
    else
      term1(1:alcsz) = term0(1:alcsz)*id1(1:alcsz)
      evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
    end if
!
  end if
!       
end
!
!-----------------------------------------------------------------------
!
! see above
!
subroutine Ven_rsp_lrel_feg(evec,rs1,rs2)
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
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer i,j,k,ii,kk,alcsz
  integer g1,g2
  RTYPE evec(MAXENERGYTERMS)
  RTYPE genmu,svec(3),pfac,pfac2
  RTYPE d2(at(rs1)%npol*at(rs2)%npol)
  RTYPE dvec(at(rs1)%npol*at(rs2)%npol,3)
  RTYPE d1(at(rs1)%npol*at(rs2)%npol),id1(at(rs1)%npol*at(rs2)%npol)
  RTYPE id1s(at(rs1)%npol*at(rs2)%npol)
  RTYPE term0(at(rs1)%npol*at(rs2)%npol)
  RTYPE termr(at(rs1)%npol*at(rs2)%npol)
  RTYPE term1(at(rs1)%npol*at(rs2)%npol)
  RTYPE term2(at(rs1)%npol*at(rs2)%npol)
!
  if (use_FEGS(6).EQV..true.) then
!
    if (rs1.gt.rs2) then
      call dis_bound_rs(rs2,rs1,svec)
    else
      call dis_bound_rs(rs1,rs2,svec)
    end if
!
    pfac = electric*par_FEG2(9)
    if (lrel_mc.eq.1) then
      k = 0
      if (use_FEGS(4).EQV..true.) then
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
            end do
            end do
          end do
          end do
        else if (scrq_model.eq.4) then
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
              term0(k) = atq(ii)*atq(kk)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*scrq(ii)*scrq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
            end do
            end do
          end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
            end do
            end do
          end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
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
              term0(k) = atq(ii)*atq(kk)*genmu(scrq(ii),scrq(kk),i_sqm)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
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
            term0(k) = atq(ii)*atq(kk)
          end do
          end do
        end do
        end do
      end if
      alcsz = k
    else if (lrel_mc.eq.2) then
      k = 0
      if (use_FEGS(4).EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
              end do
            end do
            end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = atq(ii)*atq(kk) 
              term2(k) = 1.0/(atr(ii)+atr(kk))
              end do
            end do
            end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
            end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
              end do
            end do
            end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            do i=1,at(rs1)%dpgrp(g1)%nats
            ii = at(rs1)%dpgrp(g1)%ats(i)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              do j=1,at(rs2)%dpgrp(g2)%nats
              kk = at(rs2)%dpgrp(g2)%ats(j)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*atq(ii)*atq(kk)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = atq(ii)*atq(kk)/termr(k)
              end do
            end do
            end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          do i=1,at(rs1)%dpgrp(g1)%nats
          ii = at(rs1)%dpgrp(g1)%ats(i)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            do j=1,at(rs2)%dpgrp(g2)%nats
            kk = at(rs2)%dpgrp(g2)%ats(j)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = atq(ii)*atq(kk)
            end do
          end do
          end do
        end do
      end if
      alcsz = k
    else if (lrel_mc.eq.3) then
      k = 0
      if (use_FEGS(4).EQV..true.) then
        if (scrq_model.le.2) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if (scrq_model.eq.4) then
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              term2(k) = 1.0/(atr(ii)+atr(kk))
            end do
          end do
        else if ((scrq_model.eq.3).OR.(scrq_model.eq.9)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = scrq(ii)*scrq(kk)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else if ((scrq_model.eq.5).OR.(scrq_model.eq.6)) then
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        else ! 7,8
          pfac2 = par_IMPSOLV(8)*electric
          do g1=1,at(rs1)%ndpgrps
            if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
            ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
            do g2=1,at(rs2)%ndpgrps
              if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
              kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
              k = k + 1
              dvec(k,1) = x(kk) - x(ii)
              dvec(k,2) = y(kk) - y(ii)
              dvec(k,3) = z(kk) - z(ii)
              term0(k) = genmu(scrq(ii),scrq(kk),i_sqm)*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
              termr(k) = atr(ii)+atr(kk)
              term2(k) = (1.0/termr(k))*cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
            end do
          end do
        end if
      else
        do g1=1,at(rs1)%ndpgrps
          if (at(rs1)%dpgrp(g1)%nc.eq.0) cycle
          ii = cglst%it(at(rs1)%dpgrp(g1)%cgn)
          do g2=1,at(rs2)%ndpgrps
            if (at(rs2)%dpgrp(g2)%nc.eq.0) cycle
            kk = cglst%it(at(rs2)%dpgrp(g2)%cgn)
            k = k + 1
            dvec(k,1) = x(kk) - x(ii)
            dvec(k,2) = y(kk) - y(ii)
            dvec(k,3) = z(kk) - z(ii)
            term0(k) = cglst%tc(at(rs1)%dpgrp(g1)%cgn)*cglst%tc(at(rs2)%dpgrp(g2)%cgn)
          end do
        end do
      end if
      alcsz = k
    end if ! which lrel_mc
!
!   now perform the vectorizable bulk operations (maximal)
    dvec(1:alcsz,1) = dvec(1:alcsz,1) + svec(1)
    dvec(1:alcsz,2) = dvec(1:alcsz,2) + svec(2)
    dvec(1:alcsz,3) = dvec(1:alcsz,3) + svec(3)
    d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
    d1(1:alcsz) = sqrt(d2(1:alcsz))
    id1(1:alcsz) = 1.0/d1(1:alcsz)
    id1s(1:alcsz) = 1.0/(d1(1:alcsz)+par_FEG2(10))
!   finish up -> screening model specificity is high
    if (use_FEGS(4).EQV..true.) then
      k = 0
      if ((scrq_model.le.2).OR.(scrq_model.eq.5).OR.(scrq_model.eq.6)) then
        term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
      else if (scrq_model.eq.4) then
        do k=1,alcsz
          if (id1(k).le.term2(k)) then
            term1(k) = term0(k)*id1(k)
          else
           term1(k) = term0(k)*term2(k) 
          end if
        end do           
        evec(6) = evec(6) + par_IMPSOLV(8)*pfac*sum(id1s(1:alcsz)*term1(1:alcsz))
      else    !  #3,7,8,9
        term2(1:alcsz) = term2(1:alcsz)*pfac2
        term0(1:alcsz) = term0(1:alcsz)*electric
        do k=1,alcsz
          if (abs(term0(k)).gt.abs(term2(k))) then
            term1(k) = term0(k)
          else
!        dissociated
            if ((termr(k)+par_IMPSOLV(1)).le.d1(k)) then
              term1(k) = term0(k)
!           close contact
            else if (d1(k).lt.termr(k)) then
              term1(k) = par_IMPSOLV(9)*(term2(k)-term0(k)) + term0(k)
!           first-shell regime
            else
              term1(k) = (termr(k)+par_IMPSOLV(1)-d1(k))/par_IMPSOLV(1)
              term1(k)= (1.0-par_IMPSOLV(9)*term1(k))*term0(k) + par_IMPSOLV(9)*term1(k)*term2(k) 
            end if
          end if
        end do
        term1(1:alcsz) = term1(1:alcsz)*id1s(1:alcsz)
        evec(6) = evec(6) + par_FEG2(9)*sum(term1(1:alcsz))
      end if
    else
      term1(1:alcsz) = term0(1:alcsz)*id1s(1:alcsz)
      evec(6) = evec(6) + pfac*sum(term1(1:alcsz))
    end if
!
  end if
!       
end
!
!--------------------------------------------------------------------------------------------
!
! for residues with just a single atom, many aspects of the residue pair functions become ineffective;
! however, large numbers of polar interactions with single atoms need to be computed with explicit ion baths
! and LREL_MC not 4.
! this is what this function is for: it must only be called with residues rsl:(rsl+nrs-1) forming a contiguous stretch
! of monoatomic ions with POLAR turned on
!
! WARNING: function is untested
!
subroutine Ven_ipl_long(evec,rs1,rsl,nrs)
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
  integer, INTENT(IN):: rs1,rsl,nrs
!
  integer rs,rsh,ii,iii,j,k,i,imol,imol2,lons,lolst(at(rs1)%na),g1
  logical ismember,ispresent,isinit
!
  RTYPE evec(MAXENERGYTERMS)
!
  RTYPE ixis(6),genmu
  RTYPE term0(nrs),termr(nrs),term1(nrs),term2(nrs)
  RTYPE drav(nrs,3),srav(nrs,3),dvec(nrs,3),d1(nrs),d2(nrs),id1(nrs)
  RTYPE lop1,lop3,pfac,hlp1,hlp2
!
  if ((ideal_run.EQV..true.).OR.((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.))) return
  if (nrs.le.0) return
!
  imol = molofrs(rs1)
!
! in GC, we can safely exit if rs1 is not physically present (no intra-residue interactions possible)
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    if (ismember(ispresent,imol).EQV..false.) return
  end if
!
  rsh = rsl + nrs - 1
!
  drav(1:nrs,1) = x(rsinfo(rsl,1):rsinfo(rsh,1)) - x(refat(rs1))
  drav(1:nrs,2) = y(rsinfo(rsl,1):rsinfo(rsh,1)) - y(refat(rs1))
  drav(1:nrs,3) = z(rsinfo(rsl,1):rsinfo(rsh,1)) - z(refat(rs1))
! PBC
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        lop1 = 1.0/bnd_params(j)
        lop3 = -bnd_params(j)
        srav(1:nrs,j) = lop3*anint(drav(1:nrs,j)*lop1)
      end do
    else if (bnd_shape.eq.3) then
      srav(1:nrs,1:2) = 0.0
      lop1 = 1.0/bnd_params(6)
      lop3 = -bnd_params(6)
      do j=3,3
        srav(1:nrs,j) = lop3*anint(drav(1:nrs,j)*lop1)
      end do
    end if
  else
    srav(:,:) = 0.0
  end if
!
  term0(1:nrs) = atq(rsinfo(rsl,1):rsinfo(rsh,1))
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    do rs=rsl,rsh
      k = rs - rsl + 1
      imol2 = molofrs(rs)
      if (ismember(ispresent,imol2).EQV..false.) term0(k) = 0.0
    end do
  end if
  isinit = .false.
  if ((rsl.le.rs1).AND.(rsh.ge.rs1)) then
    isinit = .true.
    k = rs1 - rsl + 1
    drav(k,1) = 1.0e9 ! make sure no FPE
    term0(k) = 0.0 ! make sure term is 0.0
  end if
!
! all interaction partners for rs1 are charged, hence the selection for atoms in rs1 is always the same
! lrel_mc 1: all polar; lrel_mc 2: all in charge groups with a net charge; lrel_mc 3: centers of charge groups w/ net charge
  lons = 0
  if (lrel_mc.eq.3) then
    do g1=1,at(rs1)%ndpgrps
      if (at(rs1)%dpgrp(g1)%nc.ne.0) then
        lons = lons + 1
        lolst = cglst%it(at(rs1)%dpgrp(g1)%cgn)
      end if
    end do
  else if (lrel_mc.eq.2) then
    do g1=1,at(rs1)%ndpgrps
      if (at(rs1)%dpgrp(g1)%nc.ne.0) then
        do j=1,at(rs1)%dpgrp(g1)%nats
          lons = lons + 1
          lolst(lons) = at(rs1)%dpgrp(g1)%ats(j)
        end do
      end if
    end do
  else if (lrel_mc.eq.1) then
    do j=1,at(rs1)%npol
      lons = lons + 1
      lolst(lons) = at(rs1)%pol(j)
    end do   
  end if
  do iii=1,lons
    ii = lolst(iii) ! rsinfo(rs1,1),rsinfo(rs1,1)+at(rs1)%na-1
!
    ixis(1) = x(ii)
    ixis(2) = y(ii)
    ixis(3) = z(ii)
    if (isinit.EQV..true.) then
      dvec(1:nrs,1) = drav(1:nrs,1) + srav(1:nrs,1)
      dvec(1:nrs,2) = drav(1:nrs,2) + srav(1:nrs,2)
      dvec(1:nrs,3) = drav(1:nrs,3) + srav(1:nrs,3)
    else
      dvec(1:nrs,1) = x(rsinfo(rsl,1):rsinfo(rsh,1)) - ixis(1) + srav(1:nrs,1)
      dvec(1:nrs,2) = y(rsinfo(rsl,1):rsinfo(rsh,1)) - ixis(2) + srav(1:nrs,2)
      dvec(1:nrs,3) = z(rsinfo(rsl,1):rsinfo(rsh,1)) - ixis(3) + srav(1:nrs,3)
    end if
    ixis(4) = atq(ii)
    ixis(5) = atr(ii)
    ixis(6) = scrq(ii)
    pfac = ixis(4)*electric*scale_POLAR
    if (use_IMPSOLV.EQV..true.) then
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        do i=1,k
          term1(i) = term0(i)*genmu(ixis(6),scrq(rsinfo(rsl,1)+i-1),i_sqm)
        end do
      end if
      if (((scrq_model.ge.1).AND.(scrq_model.le.3)).OR.(scrq_model.eq.9)) then
        term1(1:nrs) = term0(1:nrs)*ixis(6)*scrq(rsinfo(rsl,1):rsinfo(rsh,1))
      end if
      if (((scrq_model.ge.7).AND.(scrq_model.le.9)).OR.(scrq_model.eq.3)) then
        termr(1:nrs) = ixis(5) + atr(rsinfo(rsl,1):rsinfo(rsh,1))
        term2(1:nrs) = term0(1:nrs)
      else if (scrq_model.eq.4) then
        termr(1:nrs) = ixis(5) + atr(rsinfo(rsl,1):rsinfo(rsh,1))
      end if
    end if
    d2(1:nrs) = dvec(1:nrs,1)**2 + dvec(1:nrs,2)**2 + dvec(1:nrs,3)**2
    d1(1:nrs) = sqrt(d2(1:nrs))
    id1(1:nrs) = 1.0/d1(1:nrs)
    if (use_IMPSOLV.EQV..true.) then
      if (((scrq_model.ge.7).AND.(scrq_model.le.9)).OR.(scrq_model.eq.3)) then
        term2(1:nrs) = par_IMPSOLV(8)*term2(1:nrs)/termr(1:nrs)
        hlp2 = 1.0/par_IMPSOLV(1)
        do i=1,k
          if (abs(term1(i)).le.abs(term2(i))) then
!           dissociated
            if ((termr(i)+par_IMPSOLV(1)).le.d1(i)) then
  !           do nothing
!           close contact
            else if (d1(i).lt.termr(i)) then
              term1(i) = par_IMPSOLV(9)*(term2(i)-term1(i)) + term1(i)
!           first-shell regime
            else
              hlp1 = (termr(i)+par_IMPSOLV(1)-d1(i))*hlp2
              term1(i) = (1.0-par_IMPSOLV(9)*hlp1)*term1(i) + par_IMPSOLV(9)*hlp1*term2(i)
            end if
          end if
        end do
        evec(6) = evec(6) + pfac*sum(term1(1:nrs)*id1(1:nrs))
      else if (scrq_model.eq.4) then
        do i=1,k
          if (d1(i).ge.termr(i)) then
            term2(i) = term0(i)*id1(i)
          else
            term2(i) = term0(i)/termr(i)
          end if
        end do
        evec(6) = evec(6) + pfac*par_IMPSOLV(8)*sum(term2(1:nrs)*id1(1:nrs))
      else if (((scrq_model.ge.1).AND.(scrq_model.le.2)).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
        evec(6) = evec(6) + pfac*sum(term1(1:nrs)*id1(1:nrs))
      end if
    else
      evec(6) = evec(6) + pfac*sum(term0(1:nrs)*id1(1:nrs))
    end if
  end do
!
end
!
!-----------------------------------------------------------------------------------------------------------------
!
