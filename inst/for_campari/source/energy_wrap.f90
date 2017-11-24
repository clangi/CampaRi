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
!------------------------------------------------------------------------
!
!
!             #################################
!             #                               #
!             # FIRST: WRAPPER ROUTINES       #
!             #        FOR GLOBAL ENERGIES    #
!             #                               #
!             #################################
!
!
!------------------------------------------------------------------------
!
! Compute total energy as used in MC runs wo/ cutoffs
!
subroutine energy(evec,eout,tpi)
!
  use sequen
  use energies
  use molecule
  use atoms
  use ems, ONLY: curmassm
#ifdef ENABLE_THREADS
  use threads
#endif
  use wl
  use tabpot
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer i,j,imol,rs,sta,sto,azero,aone,athree,asix
#ifdef ENABLE_THREADS
  integer(KIND=8) ttimer
#endif
  RTYPE evec(MAXENERGYTERMS),eout
  logical sayno
  RTYPE evec_thr(MAXENERGYTERMS)
!
  evec_thr(:) = 0.0
  if (tpi.le.1) then
    evec(:) = 0.0
    if (use_EMICRO.EQV..true.) curmassm = .false.
  end if
#ifdef ENABLE_THREADS
  if ((tpi.gt.0).AND.(use_IMPSOLV.EQV..true.)) thr_svte(:,tpi) = 0.0
#endif
  sayno = .false.
  azero = 0
  aone = 1
  athree = 3
  asix = 6
!
  if (ideal_run.EQV..false.) then
#ifdef ENABLE_THREADS
    if (use_IMPSOLV.EQV..true.) call init_svte_threads(athree,tpi)
#else
    if (use_IMPSOLV.EQV..true.) call init_svte(athree)
#endif
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    do i=thr_limits(89,tpi),thr_limits(90,tpi)
      do j=i,nseq
        call Ven_rsp(evec_thr(:),i,j,thr_svte(:,tpi),sayno)
      end do
    end do
  else
    do i=1,nseq
      do j=i,nseq
        call Ven_rsp(evec_thr(:),i,j,svte,sayno)
      end do
    end do
  end if
#else
  do i=1,nseq
    do j=i,nseq
      call Ven_rsp(evec_thr(:),i,j,svte,sayno)
    end do
  end do
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
  if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call init_svte_threads(asix,tpi)
!$OMP BARRIER
#else
    call init_svte(asix)
#endif
    if (use_POLAR.EQV..true.) then
      call setup_scrqs2(tpi)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
  if (tpi.gt.0) then
    if ((use_TABUL.EQV..true.).AND.(use_POLAR.EQV..false.)) then
      do i=tpi,nseq,thrdat%maxn
        do j=i,nseq
          if (tbp%rsmat(i,j).gt.0) call Ven_rsp_long(evec_thr(:),i,j,sayno)
        end do
      end do
    else
      do i=thr_limits(91,tpi),thr_limits(92,tpi)
        do j=i,nseq
          call Ven_rsp_long(evec_thr(:),i,j,sayno)
        end do
      end do
    end if
  else
    do i=1,nseq
      do j=i,nseq
        call Ven_rsp_long(evec_thr(:),i,j,sayno)
      end do
    end do
  end if
#else
  do i=1,nseq
    do j=i,nseq
      call Ven_rsp_long(evec_thr(:),i,j,sayno)
    end do
  end do
#endif
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (thr_dlb(6,1).gt.0) then
      if (tpi.eq.1) thr_dlb(6,2) = thr_dlb(6,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(13,tpi) = thr_timings(13,tpi) + ttimer
    end if
    sta = thr_limits(37,tpi)
    sto = thr_limits(38,tpi)
  else
    sta = 1
    sto = nseq
  end if
#else
  sta = 1
  sto = nseq
#endif
  do rs=sta,sto
    if (use_BOND(1).EQV..true.) call en_bonds(rs,evec_thr(:))
    if (use_BOND(2).EQV..true.) call en_angles(rs,evec_thr(:))
    if (use_BOND(3).EQV..true.) call en_impropers(rs,evec_thr(:))
    if (use_BOND(4).EQV..true.) call en_torsions(rs,evec_thr(:))
    if (use_BOND(5).EQV..true.) call en_cmap(rs,evec_thr(:))
    if (use_IMPSOLV.EQV..true.) call en_freesolv(evec_thr(:),rs)
    call e_boundary_rs(rs,evec_thr(:),aone)
    if (use_TOR.EQV..true.) then
      if (par_TOR2(rs).gt.0) call en_torrs(evec_thr(:),rs)
    end if
  end do
#ifdef ENABLE_THREADS
  if ((thr_dlb(6,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(14,tpi) = thr_timings(14,tpi) + ttimer
  end if
#endif
!
#ifdef ENABLE_THREADS
!$OMP SECTIONS
!$OMP SECTION
#endif
  if (use_ZSEC.EQV..true.) then
    do imol=1,nmol
      call en_zsec_gl(imol,evec_thr(:))
    end do
  end if
!
#ifdef ENABLE_THREADS
!$OMP SECTION
#endif
  if (use_POLY.EQV..true.) then
    do imol=1,nmol
      call en_poly_gl(imol,evec_thr(:),azero)
    end do
  end if
!
#ifdef ENABLE_THREADS
!$OMP SECTION
#endif
  if (use_DSSP.EQV..true.) then
    call en_dssp_gl(evec_thr(:))
  end if
#ifdef ENABLE_THREADS
!$OMP END SECTIONS
#endif
!
  if (use_EMICRO.EQV..true.) then
    call en_emicro_gl(evec_thr(:),azero,tpi)
  end if
!
  if (use_DREST.EQV..true.) then
    call edrest(evec_thr(:),tpi)
  end if
!
#ifdef ENABLE_THREADS
!$OMP CRITICAL(MC_ENGLOB_COLLECT)
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(MC_ENGLOB_COLLECT)
#endif
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
  if (do_accelsim.EQV..true.) then
    call els_manage_justE(evec,1)
    evec(:) = evec(:) + hmjam%boosts(:,1)
  end if
  eout = sum(evec(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!-----------------------------------------------------------------------
!
! Compute total energy as used in MC runs w/ cutoffs
!
subroutine energy3(evec,nblup,eout,tpi)
!
  use sequen
  use energies
  use molecule
  use atoms
  use cutoffs
  use ems, ONLY: curmassm
#ifdef ENABLE_THREADS
  use threads
#endif
  use wl
  use ems, ONLY: curmassm
!
  implicit none
!
  integer, INTENT(IN):: tpi
  logical, INTENT(IN):: nblup
!
  integer i,j,imol,rs,sta,sto,azero,aone,athree,asix,rsl,rsh,irs,frs,ix,jj
#ifdef ENABLE_THREADS
  integer(KIND=8) ttimer
  integer stas(2),stos(2),atwo
#endif
  RTYPE evec(MAXENERGYTERMS),eout
  logical sayyes
  RTYPE evec_thr(MAXENERGYTERMS)
!
  evec_thr(:) = 0.0
  if (tpi.le.1) then
    evec(:) = 0.0
    if (use_EMICRO.EQV..true.) curmassm = .false.
  end if
#ifdef ENABLE_THREADS
  atwo = 2
  if ((tpi.gt.0).AND.(use_IMPSOLV.EQV..true.)) thr_svte(:,tpi) = 0.0
#endif
  sayyes = .true.
  azero = 0
  aone = 1
  athree = 3
  asix = 6
  ix = 1
!
#ifdef ENABLE_THREADS
  if (use_IMPSOLV.EQV..true.) call init_svte_threads(athree,tpi)
#else
  if (use_IMPSOLV.EQV..true.) call init_svte(athree)
#endif
!
  if ((ideal_run.EQV..false.).AND.(nblup.EQV..true.)) then
    if (use_cutoffs.EQV..true.) then
      rsl = 1
      rsh = nseq
      frs = nseq
      irs = 1
      jj = 0
      call respairs_new(rsl,rsh,ix,irs,frs,jj,sayyes,azero,tpi)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  end if
!
#ifdef ENABLE_THREADS
  if (use_cutoffs.EQV..true.) then
    if (tpi.gt.0) then
      call get_thread_loop_bounds(ix,aone,stas,stos,tpi)
      do jj=stas(1),stos(1)
        sta = 1
        sto = thr_rsp_nbl(jj)%srnrs
        if (jj.eq.stas(1)) sta = stas(2)
        if (jj.eq.stos(1)) sto = stos(2)
        do j=sta,sto
          call Ven_rsp(evec_thr(:),thr_rsp_nbl(jj)%sr(j,1),thr_rsp_nbl(jj)%sr(j,2),thr_svte(:,tpi),use_cutoffs)
        end do
      end do
    else
      do j=1,rsp_nbl(ix)%srnrs
        call Ven_rsp(evec_thr(:),rsp_nbl(ix)%sr(j,1),rsp_nbl(ix)%sr(j,2),svte,use_cutoffs)
      end do
    end if
  else if (ideal_run.EQV..false.) then
    if (tpi.gt.0) then
      if (thr_dlb(3,1).gt.0) then
        if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
        call System_Clock(count=ttimer)
        thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
      end if
      do i=thr_limits(9,tpi),thr_limits(10,tpi)
        do j=i,nseq
          call Ven_rsp(evec_thr(:),i,j,thr_svte(:,tpi),use_cutoffs)
        end do
      end do
      if (thr_dlb(3,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer
      end if
    else
      do i=1,nseq
        do j=i,nseq
          call Ven_rsp(evec_thr(:),i,j,svte,use_cutoffs)
        end do
      end do
    end if
  end if
#else
  if (use_cutoffs.EQV..true.) then
    do j=1,rsp_nbl(ix)%srnrs
      call Ven_rsp(evec_thr(:),rsp_nbl(ix)%sr(j,1),rsp_nbl(ix)%sr(j,2),svte,use_cutoffs)
    end do
  else if (ideal_run.EQV..false.) then
    do i=1,nseq
      do j=i,nseq
        call Ven_rsp(evec_thr(:),i,j,svte,use_cutoffs)
      end do
    end do
  end if
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
  if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    call init_svte_threads(asix,tpi)
!$OMP BARRIER
#else
    call init_svte(asix)
#endif
    if (use_POLAR.EQV..true.) then
      call setup_scrqs2(tpi)
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
  if (use_cutoffs.EQV..true.) then
    if (tpi.gt.0) then
      call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
      do jj=stas(1),stos(1)
        sta = 1
        sto = thr_rsp_nbl(jj)%trnrs
        if (jj.eq.stas(1)) sta = stas(2)
        if (jj.eq.stos(1)) sto = stos(2)
        do j=sta,sto
          call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
        end do
      end do
      call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
      do jj=stas(1),stos(1)
        sta = 1
        sto = thr_rsp_nbl(jj)%lrnrs
        if (jj.eq.stas(1)) sta = stas(2)
        if (jj.eq.stos(1)) sto = stos(2)
        do j=sta,sto
          call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
        end do
      end do
    else
      do j=1,rsp_nbl(ix)%trnrs
        call Ven_rsp_long(evec_thr(:),rsp_nbl(ix)%tr(j,1),rsp_nbl(ix)%tr(j,2),use_cutoffs)
      end do
      do j=1,rsp_nbl(ix)%lrnrs
        call Ven_rsp_lrel(evec_thr(:),rsp_nbl(ix)%lr(j,1),rsp_nbl(ix)%lr(j,2))
      end do
    end if
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
    if (tpi.gt.0) then
      if (thr_dlb(7,1).gt.0) then
        if (tpi.eq.1) thr_dlb(7,2) = thr_dlb(7,2) + 1
        call System_Clock(count=ttimer)
        thr_timings(17,tpi) = thr_timings(17,tpi) + ttimer
      end if
      do i=thr_limits(41,tpi),thr_limits(42,tpi)
        do j=i,nseq
          call Ven_rsp_long(evec_thr(:),i,j,use_cutoffs)
        end do
      end do
      if (thr_dlb(7,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(18,tpi) = thr_timings(18,tpi) + ttimer
      end if
    else
      do i=1,nseq
        do j=i,nseq
          call Ven_rsp_long(evec_thr(:),i,j,use_cutoffs)
        end do
      end do
    end if
  end if
#else
  if (use_cutoffs.EQV..true.) then
    do j=1,rsp_nbl(ix)%trnrs
      call Ven_rsp_long(evec_thr(:),rsp_nbl(ix)%tr(j,1),rsp_nbl(ix)%tr(j,2),use_cutoffs)
    end do
    do j=1,rsp_nbl(ix)%lrnrs
      call Ven_rsp_lrel(evec_thr(:),rsp_nbl(ix)%lr(j,1),rsp_nbl(ix)%lr(j,2))
    end do
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
    do i=1,nseq
      do j=i,nseq
        call Ven_rsp_long(evec_thr(:),i,j,use_cutoffs)
      end do
    end do
  end if
#endif
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (thr_dlb(6,1).gt.0) then
      if (tpi.eq.1) thr_dlb(6,2) = thr_dlb(6,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(13,tpi) = thr_timings(13,tpi) + ttimer
    end if
    sta = thr_limits(37,tpi)
    sto = thr_limits(38,tpi)
  else
    sta = 1
    sto = nseq
  end if
#else
  sta = 1
  sto = nseq
#endif
  do rs=sta,sto
    if (use_BOND(1).EQV..true.) call en_bonds(rs,evec_thr(:))
    if (use_BOND(2).EQV..true.) call en_angles(rs,evec_thr(:))
    if (use_BOND(3).EQV..true.) call en_impropers(rs,evec_thr(:))
    if (use_BOND(4).EQV..true.) call en_torsions(rs,evec_thr(:))
    if (use_BOND(5).EQV..true.) call en_cmap(rs,evec_thr(:))
    if (use_IMPSOLV.EQV..true.) call en_freesolv(evec_thr(:),rs)
    call e_boundary_rs(rs,evec_thr(:),aone)
    if (use_TOR.EQV..true.) then
      if (par_TOR2(rs).gt.0) call en_torrs(evec_thr(:),rs)
    end if
  end do
#ifdef ENABLE_THREADS
  if ((thr_dlb(6,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(14,tpi) = thr_timings(14,tpi) + ttimer
  end if
#endif
!
#ifdef ENABLE_THREADS
!$OMP SECTIONS
!$OMP SECTION
#endif
  if (use_ZSEC.EQV..true.) then
    do imol=1,nmol
      call en_zsec_gl(imol,evec_thr(:))
    end do
  end if
!
#ifdef ENABLE_THREADS
!$OMP SECTION
#endif
  if (use_POLY.EQV..true.) then
    do imol=1,nmol
      call en_poly_gl(imol,evec_thr(:),azero)
    end do
  end if
!
#ifdef ENABLE_THREADS
!$OMP SECTION
#endif
  if (use_DSSP.EQV..true.) then
    call en_dssp_gl(evec_thr(:))
  end if
#ifdef ENABLE_THREADS
!$OMP END SECTIONS
#endif
!
  if (use_EMICRO.EQV..true.) then
    call en_emicro_gl(evec_thr(:),azero,tpi)
  end if
!
  if (use_DREST.EQV..true.) then
    call edrest(evec_thr(:),tpi)
  end if
!
#ifdef ENABLE_THREADS
!$OMP CRITICAL(MC_ENGLOB_COLLECT)
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(MC_ENGLOB_COLLECT)
#endif
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
  if (do_accelsim.EQV..true.) then
    call els_manage_justE(evec,1)
    evec(:) = evec(:) + hmjam%boosts(:,1)
  end if
  eout = sum(evec(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!------------------------------------------------------------------------
!
!
!             #################################
!             #                               #
!             # SECOND: WRAPPER ROUTINES      #
!             #   FOR INCREMENTAL ENERGIES    #
!             #                               #
!             #################################
!
!
! a simple wrapper routine that handles short-range energy calculations for sidechain
! moves (only relevant nonbonded terms are computed ...)
!
subroutine chi_energy_short(rs,evec,cut,mode,tpi2)
!
  use sequen
  use iounit
  use energies
  use cutoffs
  use molecule
  use fyoc
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi2
  logical, INTENT(IN):: cut
!
  RTYPE evec(MAXENERGYTERMS)
!
  integer aone,atwo,azero,tpi,emmode,j,imol,sta,sto,ix,frs,irs,rsl,rsh,jj
  logical iflg
!
#ifdef ENABLE_THREADS
  integer stas(2),stos(2)
  integer tpx
  logical OMP_IN_PARALLEL
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
#ifdef ENABLE_THREADS
  if ((tpi2.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, chi_energy_short(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  evec_thr(:) = 0.0
!
  iflg = .true.
  azero = 0 
  aone = 1
  atwo = 2
  if (mode.eq.0) then
    ix = 2
#ifdef ENABLE_THREADS
    tpx = thrdat%maxn
#endif
  else
    ix = 1
#ifdef ENABLE_THREADS
    tpx = 0
#endif
  end if
  emmode = mode
  if (emmode.eq.0) emmode = 2
  imol = molofrs(rs)
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    tpi = tpi2
  else
    tpi = 1
  end if
#else
  tpi = 1
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
#ifdef ENABLE_THREADS
      call init_svte_threads(azero,tpi2)
#else
      call init_svte(azero)
#endif
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
#ifdef ENABLE_THREADS
      call init_svte_threads(aone,tpi2)
#else
      call init_svte(aone)
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if ((ideal_run.EQV..false.).AND.(cut.EQV..true.)) then
    rsl = rs
    rsh = rs
    frs = rs
    irs = rs
    jj = rs
    call respairs_new(rsl,rsh,ix,irs,frs,jj,iflg,azero,tpi2)
  end if

  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi2.gt.0) then
        call get_thread_loop_bounds(ix,aone,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%srnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp(evec_thr(:),thr_rsp_nbl(jj)%sr(j,1),thr_rsp_nbl(jj)%sr(j,2),thr_svte(:,tpi),use_cutoffs)
          end do
        end do
      else
#endif
      do j=1,rsp_nbl(ix)%srnrs
        call Ven_rsp(evec_thr(:),rsp_nbl(ix)%sr(j,1),rsp_nbl(ix)%sr(j,2),svte,use_cutoffs)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#else
      do j=1,nseq
        call Ven_rsp(evec_thr(:),min(rs,j),max(rs,j),svte,cut)
      end do
#endif
    end if
!
  else
!   do nothing
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
! the boundary terms depend only on coordinates
  call e_boundary_rs(rs,evec_thr(:),mode)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  sta = max(rsmol(imol,1),rs-1)
  sto = min(rsmol(imol,2),rs+1)
  if (use_BOND(1).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_bonds(j,evec_thr(:))
      if ((disulf(j).gt.0).AND.((disulf(j).lt.sta).OR.(disulf(j).gt.sto))) then
        call en_bonds(disulf(j),evec_thr(:))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(2).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_angles(j,evec_thr(:))
      if ((disulf(j).gt.0).AND.((disulf(j).lt.sta).OR.(disulf(j).gt.sto))) then
        call en_angles(disulf(j),evec_thr(:))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(3).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_impropers(j,evec_thr(:))
      if ((disulf(j).gt.0).AND.((disulf(j).lt.sta).OR.(disulf(j).gt.sto))) then
        call en_impropers(disulf(j),evec_thr(:))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(4).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_torsions(j,evec_thr(:))
      if ((disulf(j).gt.0).AND.((disulf(j).lt.sta).OR.(disulf(j).gt.sto))) then
        call en_torsions(disulf(j),evec_thr(:))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(5).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_cmap(j,evec_thr(:))
      if ((disulf(j).gt.0).AND.((disulf(j).lt.sta).OR.(disulf(j).gt.sto))) then
        call en_cmap(disulf(j),evec_thr(:))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
! due to likely load imbalance, we can/should pack some more asynchronously calculated terms
  if (use_POLY.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call en_poly_gl(imol,evec_thr(:),mode)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
! a few bias terms are calculated prior to move evaluation
  if (use_DREST.EQV..true.) call edrest(evec_thr(:),tpi2) ! threaded and safe
!
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evec_thr(:),emmode,tpi2) ! threaded and safe
!
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV and populate rs_vec and rs_vec2
#ifdef ENABLE_THREADS
      call init_svte_threads(atwo,tpi2)
#else
      call init_svte(atwo)
#endif
    end if
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (use_IMPSOLV.EQV..true.) rs_vec(rs) = 1
!$OMP END SINGLE ! implied barrier
#endif
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for sidechain
! moves (only relevant terms are computed ...)
!
subroutine chi_energy_long(rs,evec,cut,mode,tpi)
!
  use iounit
  use sequen
  use energies
  use cutoffs
  use mcsums, ONLY: nstep
  use fyoc, ONLY: disulf
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,tpi,mode
  logical, INTENT(IN):: cut
!
  integer j,i,k,irs,frs,ix,azero,athree,afive
  logical iflg
  RTYPE evec(MAXENERGYTERMS)
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),sta,sto,atwo,jj
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
  evec_thr(:) = 0.0
!
  iflg = .false.
  afive = 5
  azero = 0
#ifdef ENABLE_THREADS
  atwo = 2
#endif
  athree = 3
  if (mode.eq.1) then
    ix = 2
  else
    ix = 1
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(tpi)
  end if
!     
  if (use_TOR.EQV..true.) then
    if (par_TOR2(rs).gt.0) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call en_torrs(evec_thr(:),rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
    if (disulf(rs).gt.0) then
      if (par_TOR2(disulf(rs)).gt.0) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call en_torrs(evec_thr(:),disulf(rs))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
      end if
    end if
  end if
! 
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    irs = rs
    frs = rs
    if (cut.EQV..true.) then
      if ((use_POLAR.EQV..true.).AND.(scrq_model.ne.4)) then
        j = 1
        call respairs_append(ix,irs,frs,rs,j,tpi)
      end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
        do j=thr_limits(3,tpi),thr_limits(4,tpi)
          if ((rs_vec(j).eq.1).OR.(rs_vec2(j).eq.1)) call en_freesolv(evec_thr(:),j)
        end do
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
!
      do i=1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:),i)
        end if
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=1,nseq
        do j=i,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
        call en_freesolv(evec_thr(:),i)
      end do
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      sta = thr_limits(3,tpi)
      sto = thr_limits(4,tpi)
    else
      sta = 1
      sto = nseq
    end if
    do i=sta,sto
#else
    do i=1,nseq
#endif
      if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
        call en_freesolv(evec_thr(:),i)
      end if
    end do
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do j=1,nseq
        call Ven_rsp_long(evec_thr(:),rs,j,cut)
      end do
    end if
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
#ifdef ENABLE_THREADS
      call init_svte_threads(afive,tpi)
#else
      call init_svte(afive)
#endif
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        rs_vec(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
        rs_vec2(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
      else
        rs_vec(:) = 0
        rs_vec2(:) = 0
      end if
#else
      rs_vec(:) = 0
      rs_vec2(:) = 0
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!----------------------------------------------------------------------------------
!
! this is the core routine for incremental SR nonbonded energies (ATTLJ, WCA, IPP) and
! direct interaction nb-list population for a swiveling segment comprised
! of residues irs-frs with rs being pertubed internally
!
! it handles all bonded and several bias potentials (ZSEC, DSSP, EMICRO, DREST, POLY) and
! the (soft-wall) boundary term
!
subroutine swivel_energy_short(rs,irs,frs,evec,cut,mode,tpi2)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
  use fyoc
  use atoms, ONLY: svte
  use dssps, ONLY: pepmol
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,irs,frs,mode,tpi2
  logical, INTENT(IN):: cut
!
  RTYPE evec(MAXENERGYTERMS)
!
  integer aone,atwo,azero,tpi,emmode,i,j,imol,sta,sto,ix,rsl,rsh
  logical iflg
!
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),incr
  integer tpx,jj
  logical OMP_IN_PARALLEL
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
#ifdef ENABLE_THREADS
  if ((tpi2.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, swivel_energy_short(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  evec_thr(:) = 0.0
!
  iflg = .true.
  azero = 0 
  aone = 1
  atwo = 2
  if (mode.eq.0) then
    ix = 2
#ifdef ENABLE_THREADS
    tpx = thrdat%maxn
#endif
  else
    ix = 1
#ifdef ENABLE_THREADS
    tpx = 0
#endif
  end if
  emmode = mode
  if (emmode.eq.0) emmode = 2
  imol = molofrs(rs)
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    tpi = tpi2
  else
    tpi = 1
  end if
#else
  tpi = 1
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
#ifdef ENABLE_THREADS
      call init_svte_threads(azero,tpi2)
#else
      call init_svte(azero)
#endif
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
#ifdef ENABLE_THREADS
      call init_svte_threads(aone,tpi2)
#else
      call init_svte(aone)
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if ((ideal_run.EQV..false.).AND.(cut.EQV..true.)) then
    rsl = irs
    rsh = frs
    call respairs_new(rsl,rsh,ix,irs,frs,rs,iflg,aone,tpi2)
  end if
!
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then!
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi2.gt.0) then
        call get_thread_loop_bounds(ix,aone,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%srnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp(evec_thr(:),thr_rsp_nbl(jj)%sr(j,1),thr_rsp_nbl(jj)%sr(j,2),thr_svte(:,tpi),use_cutoffs)
          end do
        end do
      else
#endif
      do j=1,rsp_nbl(ix)%srnrs
        call Ven_rsp(evec_thr(:),rsp_nbl(ix)%sr(j,1),rsp_nbl(ix)%sr(j,2),svte,use_cutoffs)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#else
      do j=irs,frs
        do sta=1,irs-1
          call Ven_rsp(evec_thr(:),sta,j,svte,cut)
        end do
        do sta=frs+1,nseq
          call Ven_rsp(evec_thr(:),j,sta,svte,cut)
        end do
        call Ven_rsp(evec_thr(:),min(rs,j),max(rs,j),svte,cut)
      end do
#endif
    end if
!
  else
!   do nothing
  end if
!
  sta = max(rsmol(imol,1),rs-1)
  sto = min(rsmol(imol,2),rs+1)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  do i=1,n_crosslinks
    if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &      ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
    if (((crosslink(i)%rsnrs(1).gt.irs).AND.(crosslink(i)%rsnrs(1).lt.frs)).AND.&
 &      ((crosslink(i)%rsnrs(2).gt.irs).AND.(crosslink(i)%rsnrs(2).lt.frs))) cycle
    if ((crosslink(i)%rsnrs(1).gt.sto).OR.(crosslink(i)%rsnrs(1).lt.sta)) then
      if (use_BOND(1).EQV..true.) call en_bonds(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(2).EQV..true.) call en_angles(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(3).EQV..true.) call en_impropers(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(4).EQV..true.) call en_torsions(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(5).EQV..true.) call en_cmap(crosslink(i)%rsnrs(1),evec_thr(:))
    end if
    if ((crosslink(i)%rsnrs(2).gt.sto).OR.(crosslink(i)%rsnrs(2).lt.sta)) then
      if (use_BOND(1).EQV..true.) call en_bonds(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(2).EQV..true.) call en_angles(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(3).EQV..true.) call en_impropers(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(4).EQV..true.) call en_torsions(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(5).EQV..true.) call en_cmap(crosslink(i)%rsnrs(2),evec_thr(:))
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
!
  if (use_BOND(1).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_bonds(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(2).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_angles(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(3).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_impropers(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(4).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_torsions(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(5).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_cmap(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
! due to likely load imbalance, we can/should pack some more asynchronously calculated terms
  if (use_POLY.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call en_poly_gl(imol,evec_thr(:),mode)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if ((use_ZSEC.EQV..true.).AND.(seqpolty(rs).eq.'P')) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call en_zsec_gl(imol,evec_thr(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_DSSP.EQV..true.) then
    if (pepmol(imol).gt.0) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call en_dssp_gl(evec_thr(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
! boundary term can involve a fair number of terms
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    sta = irs+tpi2-1
    incr = thrdat%maxn
  else
    sta = irs
    incr = 1
  end if
  do j=sta,frs,incr
#else
  do j=irs,frs
#endif
    call e_boundary_rs(j,evec_thr(:),mode)
  end do

! a few bias terms are calculated prior to move evaluation
  if (use_DREST.EQV..true.) call edrest(evec_thr(:),tpi2) ! threaded and safe
!
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evec_thr(:),emmode,tpi2) ! threaded and safe
!
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV and populate rs_vec and rs_vec2
#ifdef ENABLE_THREADS
      call init_svte_threads(atwo,tpi2)
#else
      call init_svte(atwo)
#endif
    end if
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (use_IMPSOLV.EQV..true.) rs_vec(rs) = 1
!$OMP END SINGLE ! implied barrier
#endif
!
end
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for dedicated
! pivot moves of various types - calls swivel_energy_short
!
subroutine pivot_energy_short(rs,evec,cut,mode,ct,tpi2)
!
  use sequen, ONLY: molofrs
  use energies, ONLY: MAXENERGYTERMS
  use molecule, ONLY: rsmol
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi2
  logical, INTENT(IN):: cut,ct
!
  RTYPE evec(MAXENERGYTERMS)
!
  integer frs,irs,jj,imol
!
  imol = molofrs(rs)
  if (ct.EQV..true.) then ! the stretch (N)~irs->frs~(C) is the moving arm
    frs = rs
    irs = rsmol(imol,1)
  else
    frs = rsmol(imol,2)
    irs = rs
  end if
  jj = rs
  call swivel_energy_short(jj,irs,frs,evec,cut,mode,tpi2)
!
end
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for single torsion
! dihedral moves using rotation lists - calls swivel_energy_short
!
subroutine other_energy_short(rs,irs,frs,evec,cut,mode,tpi2)
!
  use energies, ONLY: MAXENERGYTERMS
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi2,irs,frs
  logical, INTENT(IN):: cut
!
  RTYPE evec(MAXENERGYTERMS)
!
  call swivel_energy_short(rs,irs,frs,evec,cut,mode,tpi2)
!
end
!
!----------------------------------------------------------------------------------
!
! a copy of swivel_energy_sort / pivot_energy_short for special purpose application in randomization
!
subroutine randomize_energy_short(rs,evec,cut,ct,molonly,nmr,mrdat,iacnt,mode)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
  use fyoc
  use atoms, ONLY: svte
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,nmr,mrdat(nmr,2),mode
  logical, INTENT(IN):: cut,ct,molonly
!
  RTYPE evec(MAXENERGYTERMS)
!
  integer iacnt,aone,atwo,azero,tpi,i,j,imol,sta,sto,irs,frs,jmol,rsl,rsh
  logical ldum
!
  RTYPE evec_thr(MAXENERGYTERMS),dis,dis2
!
  imol = molofrs(rs)
  if (mode.eq.1) then ! stretch
    if (ct.EQV..true.) then ! the stretch (N)~irs->frs~(C) is the moving arm
      frs = rs
      irs = rsmol(imol,1)
    else
      frs = rsmol(imol,2)
      irs = rs
    end if
  else ! single residue
    irs = rs
    frs = rs
  end if
!
  iacnt = 0
  evec_thr(:) = 0.0
!
  ldum = .false.
  azero = 0 
  aone = 1
  atwo = 2
  imol = molofrs(rs)
  tpi = 1
!
  if ((use_IPP.EQV..true.).OR.(use_FEGS(1).EQV..true.)) then
    if (molonly.EQV..true.) then
      rsl = rsmol(imol,1)
      rsh = min(mrdat(nmr,2),rsmol(imol,2))
    else
      rsl = 1
      rsh = mrdat(nmr,2)
    end if
    do j=irs,frs
      do sta=rsl,irs-1
        call dis_bound(refat(j),refat(sta),dis2)
        dis = sqrt(dis2)
        if (dis.lt.(mcnb_cutoff+resrad(j)+resrad(sta))) then
          iacnt = iacnt + 1
          call Ven_rsp(evec_thr(:),sta,j,svte,cut)
        end if
      end do
      do sta=frs+1,rsh
        jmol = molofrs(sta)
!       cycle with not yet built tails of additional molecules or with unbuilt C-terminal tail of same mol.
        if (((imol.ne.jmol).AND.(sta.gt.mrdat(jmol,2).OR.(sta.lt.mrdat(jmol,1)))).OR.&
 &          ((imol.eq.jmol).AND.((rs.lt.mrdat(imol,1)).AND.(sta.gt.mrdat(imol,2))))) then
          cycle
        end if
        call dis_bound(refat(j),refat(sta),dis2)
        dis = sqrt(dis2)
        if (dis.lt.(mcnb_cutoff+resrad(j)+resrad(sta))) then
          iacnt = iacnt + 1
          call Ven_rsp(evec_thr(:),j,sta,svte,cut)
        end if
      end do
      call dis_bound(refat(j),refat(rs),dis2)
      dis = sqrt(dis2)
      if (dis.lt.(mcnb_cutoff+resrad(j)+resrad(rs))) then
        iacnt = iacnt + 1
        call Ven_rsp(evec_thr(:),min(rs,j),max(rs,j),svte,cut)
      end if
    end do
  else
    iacnt = frs-irs+1
  end if
!
  sta = max(rsmol(imol,1),rs-1)
  sto = min(rsmol(imol,2),rs+1)
  do i=1,n_crosslinks
    if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &      ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
    if (((crosslink(i)%rsnrs(1).gt.irs).AND.(crosslink(i)%rsnrs(1).lt.frs)).AND.&
 &      ((crosslink(i)%rsnrs(2).gt.irs).AND.(crosslink(i)%rsnrs(2).lt.frs))) cycle
!   if a crosslink references an as-of-yet unbuilt stretch, it also needs to be discarded
    if ((crosslink(i)%rsnrs(1).ge.irs).AND.(crosslink(i)%rsnrs(1).le.frs)) then
      jmol = molofrs(crosslink(i)%rsnrs(2))
      if ((jmol.gt.imol).AND.(jmol.le.nmr)) then
        if ((crosslink(i)%rsnrs(2).lt.mrdat(jmol,1)).OR.(crosslink(i)%rsnrs(2).gt.mrdat(jmol,2))) then
          cycle
        end if
      else if (jmol.gt.imol) then
        cycle
      else if ((imol.eq.jmol).AND.(imol.le.nmr)) then
        if ((crosslink(i)%rsnrs(1).lt.mrdat(imol,1)).AND.(crosslink(i)%rsnrs(2).gt.mrdat(imol,2))) then
          cycle
        end if
      else if (imol.eq.jmol) then
!        should not happen
      else
!        do nothing
      end if
    else if ((crosslink(i)%rsnrs(2).ge.irs).AND.(crosslink(i)%rsnrs(2).le.frs)) then
      jmol = molofrs(crosslink(i)%rsnrs(1))
      if ((jmol.gt.imol).AND.(jmol.le.nmr)) then
        if ((crosslink(i)%rsnrs(1).lt.mrdat(jmol,1)).OR.(crosslink(i)%rsnrs(1).gt.mrdat(jmol,2))) then
          cycle
        end if
      else if (jmol.gt.imol) then
        cycle
      else if ((imol.eq.jmol).AND.(imol.le.nmr)) then
        if ((crosslink(i)%rsnrs(2).lt.mrdat(imol,1)).AND.(crosslink(i)%rsnrs(1).gt.mrdat(imol,2))) then
          cycle
        end if
      else if (imol.eq.jmol) then
!        should not happen
      else
!        do nothing
      end if
    end if
    if ((crosslink(i)%rsnrs(1).gt.sto).OR.(crosslink(i)%rsnrs(1).lt.sta)) then
      if (use_BOND(1).EQV..true.) call en_bonds(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(2).EQV..true.) call en_angles(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(3).EQV..true.) call en_impropers(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(4).EQV..true.) call en_torsions(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(5).EQV..true.) call en_cmap(crosslink(i)%rsnrs(1),evec_thr(:))
    end if
    if ((crosslink(i)%rsnrs(2).gt.sto).OR.(crosslink(i)%rsnrs(2).lt.sta)) then
      if (use_BOND(1).EQV..true.) call en_bonds(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(2).EQV..true.) call en_angles(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(3).EQV..true.) call en_impropers(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(4).EQV..true.) call en_torsions(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(5).EQV..true.) call en_cmap(crosslink(i)%rsnrs(2),evec_thr(:))
    end if
  end do
!
  if (use_BOND(1).EQV..true.) then
    do j=sta,sto
      if (disulf(j).gt.j) then
        jmol = molofrs(disulf(j))
        if (jmol.gt.nmr) then
          ldum = .true.
        else if (imol.eq.jmol) then
          if ((j.lt.mrdat(imol,1)).AND.(disulf(j).gt.mrdat(imol,2))) ldum = .true.
        else if ((disulf(j).lt.mrdat(jmol,1)).OR.(disulf(j).gt.mrdat(jmol,2))) then
          ldum = .true.
        end if
      end if
      call en_bonds_alt(j,evec_thr(:),ldum)
      ldum = .false.
    end do
  end if
  if (use_BOND(2).EQV..true.) then
    do j=sta,sto
      if (disulf(j).gt.j) then
        jmol = molofrs(disulf(j))
        if (jmol.gt.nmr) then
          ldum = .true.
        else if (imol.eq.jmol) then
          if ((j.lt.mrdat(imol,1)).AND.(disulf(j).gt.mrdat(imol,2))) ldum = .true.
        else if ((disulf(j).lt.mrdat(jmol,1)).OR.(disulf(j).gt.mrdat(jmol,2))) then
          ldum = .true.
        end if
      end if
      call en_angles_alt(j,evec_thr(:),ldum)
      ldum = .false.
    end do
  end if
  if (use_BOND(3).EQV..true.) then
    do j=sta,sto
      call en_impropers(j,evec_thr(:))
    end do
  end if
  if (use_BOND(4).EQV..true.) then
    do j=sta,sto
      call en_torsions(j,evec_thr(:))
    end do
  end if
  if (use_BOND(5).EQV..true.) then
    do j=sta,sto
      call en_cmap(j,evec_thr(:))
    end do
  end if
!
! boundary term can involve a fair number of terms
  do j=irs,frs
    call e_boundary_rs(j,evec_thr(:),aone)
  end do

! a few bias terms are calculated prior to move evaluation
  if (use_DREST.EQV..true.) call edrest(evec_thr(:),azero)
!
  evec(:) = evec(:) + evec_thr(:)
!
end

!
!--------------------------------------------------------------------------
!
! this is the core routine for incremental nonbonded LR energies (POLAR, TABUL) and
! potential appending of interaction nb-lists by indirect terms for a swiveling segment comprised
! of residues irs-frs with rs being pertubed internally
!
! it also handles TOR and the actual calculation of IMPSOLV
!
subroutine swivel_energy_long(rs,irs,frs,evec,cut,mode,tpi)
!
  use iounit
  use sequen
  use energies
  use cutoffs
  use molecule, ONLY: rsmol
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,tpi,mode,irs,frs
  logical, INTENT(IN):: cut
!
  integer j,i,k,ix,azero,afive,imol,sta,sto
  logical iflg
  RTYPE evec(MAXENERGYTERMS)
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),atwo,athree,jj
  logical OMP_IN_PARALLEL
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
#ifdef ENABLE_THREADS
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, swivel_energy_long(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  evec_thr(:) = 0.0
  imol = molofrs(rs)
!
  iflg = .false.
  azero = 0
  afive = 5
#ifdef ENABLE_THREADS
  atwo = 2
  athree = 3
#endif
  if (mode.eq.1) then
    ix = 2
  else
    ix = 1
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(tpi)
  end if
!     
  if (use_TOR.EQV..true.) then
    sta = max(rsmol(imol,1),rs-1)
    sto = min(rsmol(imol,2),rs+1)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do i=1,n_crosslinks
      if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &        ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
      if (((crosslink(i)%rsnrs(1).gt.sto).OR.(crosslink(i)%rsnrs(1).lt.sta)).AND.(par_TOR2(crosslink(i)%rsnrs(1)).gt.0)) then
        call en_torrs(evec_thr(:),crosslink(i)%rsnrs(1))
      end if
      if (((crosslink(i)%rsnrs(2).gt.sto).OR.(crosslink(i)%rsnrs(2).lt.sta)).AND.(par_TOR2(crosslink(i)%rsnrs(2)).gt.0)) then
        call en_torrs(evec_thr(:),crosslink(i)%rsnrs(2))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      if (par_TOR2(j).gt.0) call en_torrs(evec_thr(:),j)
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
! 
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.))) then
!
    if (cut.EQV..true.) then
      if ((use_POLAR.EQV..true.).AND.(scrq_model.ne.4)) then
        j = 0
        call respairs_append(ix,irs,frs,rs,j,tpi)
      end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
        do j=thr_limits(3,tpi),thr_limits(4,tpi)
          if ((rs_vec(j).eq.1).OR.(rs_vec2(j).eq.1)) call en_freesolv(evec_thr(:),j)
        end do
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
!
      do i=1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:),i)
        end if
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=irs,frs
        do j=1,irs-1
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
        do j=i,frs
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=1,irs-1
        do j=i,irs-1
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=frs+1,nseq
        do j=i,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=1,nseq
        call en_freesolv(evec_thr(:),i)
      end do
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      sta = thr_limits(3,tpi)
      sto = thr_limits(4,tpi)
    else
      sta = 1
      sto = nseq
    end if
    do i=sta,sto
#else
    do i=1,nseq
#endif
      if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
        call en_freesolv(evec_thr(:),i)
      end if
    end do
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=irs,frs
        do j=1,irs-1
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
        call Ven_rsp_long(evec_thr(:),min(rs,i),max(rs,i),cut)
      end do
    end if
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
#ifdef ENABLE_THREADS
      call init_svte_threads(afive,tpi)
#else
      call init_svte(afive)
#endif
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        rs_vec(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
        rs_vec2(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
      else
        rs_vec(:) = 0
        rs_vec2(:) = 0
      end if
#else
      rs_vec(:) = 0
      rs_vec2(:) = 0
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for dedicated
! pivot moves of various types - calls swivel_energy_long
!
subroutine pivot_energy_long(rs,evec,cut,mode,ct,tpi)
!
  use sequen, ONLY: molofrs
  use energies, ONLY: MAXENERGYTERMS
  use molecule, ONLY: rsmol
!
  implicit none
!
  integer, INTENT(IN):: rs,tpi,mode
  logical, INTENT(IN):: cut,ct
!
  RTYPE evec(MAXENERGYTERMS)
  integer imol,frs,irs,jj
!
  imol = molofrs(rs)
!
  jj = rs
  if (ct.EQV..true.) then ! the stretch (N)~irs->frs~(C) is the moving arm
    frs = rs
    irs = rsmol(imol,1)
  else
    frs = rsmol(imol,2)
    irs = rs
  end if
!
  call swivel_energy_long(jj,irs,frs,evec,cut,mode,tpi)
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for single torsion
! dihedral moves using rotation lists - calls swivel_energy_long
!
subroutine other_energy_long(rs,irs,frs,evec,cut,mode,tpi)
!
  use energies, ONLY: MAXENERGYTERMS
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi,irs,frs
  logical, INTENT(IN):: cut
!
  RTYPE evec(MAXENERGYTERMS)
!
  call swivel_energy_long(rs,irs,frs,evec,cut,mode,tpi)
!
end

!
!----------------------------------------------------------------------------------
!
! this is the core routine for incremental SR nonbonded energies (ATTLJ, WCA, IPP) and
! direct interaction nb-list population for a flexible segment comprised
! of residues irs-frs 
!
! it handles all bonded and several bias potentials (ZSEC, DSSP, EMICRO, DREST, POLY) and
! the (soft-wall) boundary term
!
! the only difference to swivel_energy_short is the inclusion of all intra-stretch interactions
!
subroutine stretch_energy_short(irs,frs,evec,cut,mode,tpi2)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
  use fyoc
  use atoms, ONLY: svte
  use dssps, ONLY: pepmol
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: irs,frs,mode,tpi2
  logical, INTENT(IN):: cut
!
  RTYPE evec(MAXENERGYTERMS)
!
  integer aone,atwo,azero,tpi,emmode,i,j,imol,sta,sto,ix,rsl,rsh
  logical iflg
!
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),incr
  integer tpx,jj
  logical OMP_IN_PARALLEL
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
#ifdef ENABLE_THREADS
  if ((tpi2.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, stretch_energy_short(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  evec_thr(:) = 0.0
!
  iflg = .true.
  azero = 0 
  aone = 1
  atwo = 2
  if (mode.eq.0) then
    ix = 2
#ifdef ENABLE_THREADS
    tpx = thrdat%maxn
#endif
  else
    ix = 1
#ifdef ENABLE_THREADS
    tpx = 0
#endif
  end if
  emmode = mode
  if (emmode.eq.0) emmode = 2
  imol = molofrs(irs)
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    tpi = tpi2
  else
    tpi = 1
  end if
#else
  tpi = 1
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
#ifdef ENABLE_THREADS
      call init_svte_threads(azero,tpi2)
#else
      call init_svte(azero)
#endif
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
#ifdef ENABLE_THREADS
      call init_svte_threads(aone,tpi2)
#else
      call init_svte(aone)
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
  if ((ideal_run.EQV..false.).AND.(cut.EQV..true.)) then
    rsl = irs
    rsh = frs
    j = irs
    call respairs_new(rsl,rsh,ix,irs,frs,j,iflg,azero,tpi2)
  end if
!
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi2.gt.0) then
        call get_thread_loop_bounds(ix,aone,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%srnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp(evec_thr(:),thr_rsp_nbl(jj)%sr(j,1),thr_rsp_nbl(jj)%sr(j,2),thr_svte(:,tpi),use_cutoffs)
          end do
        end do
      else
#endif
      do j=1,rsp_nbl(ix)%srnrs
        call Ven_rsp(evec_thr(:),rsp_nbl(ix)%sr(j,1),rsp_nbl(ix)%sr(j,2),svte,use_cutoffs)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#else
      do j=irs,frs
        do sta=1,irs-1
          call Ven_rsp(evec_thr(:),sta,j,svte,cut)
        end do
        do sta=frs+1,nseq
          call Ven_rsp(evec_thr(:),j,sta,svte,cut)
        end do
        do sta=j,frs
          call Ven_rsp(evec_thr(:),j,sta,svte,cut)
        end do
      end do
#endif
    end if
!
  else
!   do nothing
  end if
!
  sta = max(rsmol(imol,1),irs-1)
  sto = min(rsmol(imol,2),frs+1)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  do i=1,n_crosslinks
    if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &      ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
    if ((crosslink(i)%rsnrs(1).gt.sto).OR.(crosslink(i)%rsnrs(1).lt.sta)) then
      if (use_BOND(1).EQV..true.) call en_bonds(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(2).EQV..true.) call en_angles(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(3).EQV..true.) call en_impropers(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(4).EQV..true.) call en_torsions(crosslink(i)%rsnrs(1),evec_thr(:))
      if (use_BOND(5).EQV..true.) call en_cmap(crosslink(i)%rsnrs(1),evec_thr(:))
    end if
    if ((crosslink(i)%rsnrs(2).gt.sto).OR.(crosslink(i)%rsnrs(2).lt.sta)) then
      if (use_BOND(1).EQV..true.) call en_bonds(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(2).EQV..true.) call en_angles(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(3).EQV..true.) call en_impropers(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(4).EQV..true.) call en_torsions(crosslink(i)%rsnrs(2),evec_thr(:))
      if (use_BOND(5).EQV..true.) call en_cmap(crosslink(i)%rsnrs(2),evec_thr(:))
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
!
  if (use_BOND(1).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_bonds(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(2).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_angles(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(3).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_impropers(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(4).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_torsions(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_BOND(5).EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=sta,sto
      call en_cmap(j,evec_thr(:))
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
! due to likely load imbalance, we can/should pack some more asynchronously calculated terms
  if (use_POLY.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call en_poly_gl(imol,evec_thr(:),mode)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  rsl = (irs+frs)/2 ! to avoid errors from caps
  if ((use_ZSEC.EQV..true.).AND.(seqpolty(rsl).eq.'P')) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call en_zsec_gl(imol,evec_thr(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if (use_DSSP.EQV..true.) then
    if (pepmol(imol).gt.0) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call en_dssp_gl(evec_thr(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
! boundary term can involve a fair number of terms
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    sta = irs+tpi2-1
    incr = thrdat%maxn
  else
    sta = irs
    incr = 1
  end if
  do j=sta,frs,incr
#else
  do j=irs,frs
#endif
    call e_boundary_rs(j,evec_thr(:),mode)
  end do

! a few bias terms are calculated prior to move evaluation
  if (use_DREST.EQV..true.) call edrest(evec_thr(:),tpi2) ! threaded and safe
!
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evec_thr(:),emmode,tpi2) ! threaded and safe
!
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV and populate rs_vec and rs_vec2
#ifdef ENABLE_THREADS
      call init_svte_threads(atwo,tpi2)
#else
      call init_svte(atwo)
#endif
    end if
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
  if (use_IMPSOLV.EQV..true.) rs_vec(irs:frs) = 1
!$OMP END SINGLE ! implied barrier
#endif
!
end

!
!--------------------------------------------------------------------------
!
! this is the core routine for incremental nonbonded LR energies (POLAR, TABUL) and
! potential appending of interaction nb-lists by indirect terms for a flexible stretch comprised
! of residues irs-frs moving
!
! it also handles TOR and the actual calculation of IMPSOLV
!
! the only difference to swivel_energy_long is the treatment of intra-stretch interactions
!
subroutine stretch_energy_long(irs,frs,evec,cut,mode,tpi)
!
  use iounit
  use sequen
  use energies
  use cutoffs
  use molecule, ONLY: rsmol
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,mode,irs,frs
  logical, INTENT(IN):: cut
!
  integer j,i,k,ix,azero,afive,imol,sta,sto
  logical iflg
  RTYPE evec(MAXENERGYTERMS)
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),atwo,athree,jj
  logical OMP_IN_PARALLEL
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
#ifdef ENABLE_THREADS
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, stretch_energy_long(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  evec_thr(:) = 0.0
  imol = molofrs(irs)
!
  iflg = .false.
  azero = 0
  afive = 5
#ifdef ENABLE_THREADS
  atwo = 2
  athree = 3
#endif
  if (mode.eq.1) then
    ix = 2
  else
    ix = 1
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(tpi)
  end if
!     
  if (use_TOR.EQV..true.) then
    sta = max(rsmol(imol,1),irs-1)
    sto = min(rsmol(imol,2),frs+1)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do i=1,n_crosslinks
      if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &        ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
      if (((crosslink(i)%rsnrs(1).gt.sto).OR.(crosslink(i)%rsnrs(1).lt.sta)).AND.(par_TOR2(crosslink(i)%rsnrs(1)).gt.0)) then
        call en_torrs(evec_thr(:),crosslink(i)%rsnrs(1))
      end if
      if (((crosslink(i)%rsnrs(2).gt.sto).OR.(crosslink(i)%rsnrs(2).lt.sta)).AND.(par_TOR2(crosslink(i)%rsnrs(2)).gt.0)) then
        call en_torrs(evec_thr(:),crosslink(i)%rsnrs(2))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
#ifdef ENABLE_THREADS
    do j=sta+tpi-1,sto,thrdat%maxn
#else
    do j=sta,sto
#endif
      if (par_TOR2(j).gt.0) call en_torrs(evec_thr(:),j)
    end do
  end if
! 
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.))) then
!
    if (cut.EQV..true.) then
      if ((use_POLAR.EQV..true.).AND.(scrq_model.ne.4)) then
        j = 1
        call respairs_append(ix,irs,frs,azero,j,tpi)
      end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
        do j=thr_limits(3,tpi),thr_limits(4,tpi)
          if ((rs_vec(j).eq.1).OR.(rs_vec2(j).eq.1)) call en_freesolv(evec_thr(:),j)
        end do
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
!
      do i=1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:),i)
        end if
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=irs,frs
        do j=1,irs-1
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
        do j=i,frs
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=1,irs-1
        do j=i,irs-1
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=frs+1,nseq
        do j=i,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=1,nseq
        call en_freesolv(evec_thr(:),i)
      end do
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      sta = thr_limits(3,tpi)
      sto = thr_limits(4,tpi)
    else
      sta = 1
      sto = nseq
    end if
    do i=sta,sto
#else
    do i=1,nseq
#endif
      if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
        call en_freesolv(evec_thr(:),i)
      end if
    end do
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=irs,frs
        do j=1,irs-1
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
        do j=i,frs
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
    end if
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
#ifdef ENABLE_THREADS
      call init_svte_threads(afive,tpi)
#else
      call init_svte(afive)
#endif
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        rs_vec(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
        rs_vec2(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
      else
        rs_vec(:) = 0
        rs_vec2(:) = 0
      end if
#else
      rs_vec(:) = 0
      rs_vec2(:) = 0
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end

!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine cr_energy_short(rsi,rsf,evec,cut,mode,ct)
!
  use sequen
  use iounit
  use energies
  use cutoffs
  use molecule
  use fyoc
  use atoms,ONLY: svte
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,k,mode,irs,frs,imol,sta,sto,sta2,sto2,incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut,ct
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
!
!
! set residue range (lever arm) based on whether N- or C-aligned (the moving parts change!)
  if (ct.EQV..true.) then
    irs = rsmol(molofrs(rsi),1)
    frs = rsf
  else
    irs = rsi
    frs = rsmol(molofrs(rsi),2)
  end if
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
  sta = irs
  sto = frs
  sta2 = max(rsmol(molofrs(rsi),1),irs-1)
  sto2 = min(rsmol(molofrs(rsi),2),frs+1)
  incr = 1
  tpi = 1
  do j=sta,sto,incr
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
  end do
  if (n_crosslinks.gt.0) then
    do j=sta2,sto2,incr
      if (disulf(j).gt.0) then
        if ((j.le.(rsf+1)).AND.(j.ge.(rsi-1)).AND.&
 &          ((disulf(j).gt.(rsf+1)).OR.(disulf(j).lt.(rsi-1)))) then
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        else if ((disulf(j).gt.sto2).OR.(disulf(j).lt.sta2)) then
          if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
        end if
      end if
    end do
  end if
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
  imol = molofrs(rsi)
  sta = max(rsmol(imol,1),rsi-1)
  sto = min(rsmol(imol,2),rsf+1)
  tpi = 1
  incr = 1
  do j=sta,sto,incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
  end do
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do j=irs,frs
          call grd_respairs(j)
        end do
      else if (use_rescrit.EQV..true.) then
        do j=irs,frs
          call respairs(j)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in cr_energy_short&
 &(...).'
        call fexit()
      end if
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
      do k=sta,sto,incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
            end if
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
            end if
          end do
        end if
      end do
      do j=irs,frs
        call clear_rsp2(j)
      end do
!
    else
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
      do k=sta,sto,incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
          end do
        end if
      end do
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine cr_energy_long(rsi,rsf,evec,cut,mode,ct)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,mode,i,k,l,rs,nrsl,rsl(nseq),frs,irs,sta(2),sto(2),incr,tpi,azero,imol
  RTYPE evec(MAXENERGYTERMS)
  logical cut,ct,inrsl(nseq)
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
  imol = molofrs(rsi)
!
  if (ct.EQV..true.) then
    irs = rsmol(imol,1)
    frs = rsf
  else
    irs = rsi
    frs = rsmol(imol,2)
  end if
!
  azero = 0
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(azero)
  end if
!
  nrsl = 0
!
  if (use_TOR.EQV..true.) then
    sta(1) = max(rsmol(imol,1),rsi-1)
    sto(1) = min(rsmol(imol,2),rsf+1)
    do i=1,n_crosslinks
      if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &        ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
      if (((crosslink(i)%rsnrs(1).gt.sto(1)).OR.(crosslink(i)%rsnrs(1).lt.sta(1))).AND.(par_TOR2(crosslink(i)%rsnrs(1)).gt.0)) then
        call en_torrs(evec_thr(:,1),crosslink(i)%rsnrs(1))
      end if
      if (((crosslink(i)%rsnrs(2).gt.sto(1)).OR.(crosslink(i)%rsnrs(2).lt.sta(1))).AND.(par_TOR2(crosslink(i)%rsnrs(2)).gt.0)) then
        call en_torrs(evec_thr(:,1),crosslink(i)%rsnrs(2))
      end if
    end do
    do j=sta(1),sto(1)
      if (par_TOR2(j).gt.0) call en_torrs(evec_thr(:,1),j)
    end do
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do rs=1,irs-1
       if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      do rs=irs,frs
        nrsl = nrsl + 1
        rsl(nrsl) = rs
        inrsl(rs) = .true.
      end do
      do rs=frs+1,nseq
        if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in cr_energy_lon&
 &g(...).'
        call fexit()
      end if
!
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr ! 1,nseq
!         we use the inverse map to circumvent problems with double counting
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
!       and do the rest here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr !i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
!
      do i=sta(1),sto(1),incr !1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
!
    else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
    tpi = 1
    sta(1) = 1
    sto(1) = nseq
    incr = 1
    do j=sta(1),sto(1),incr !1,nseq
      call en_freesolv(evec_thr(:,tpi),j)
    end do
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in pivot_energy_&
 &long(...).'
        call fexit()
      end if
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      do k=sta(1),sto(1),incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        end if
      end do
      do j=irs,frs
        call clear_rsp2(j)
      end do
!
    else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      do k=sta(1),sto(1),incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else 
!         CR-stretch to rest of lever arm
          do j=rsi,rsf
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        end if
      end do
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for truly
! local CR moves (only relevant terms are computed ...)
! for long peptides, the complexity is significantly reduced when compared to pivot
! moves
!
subroutine ujcr_energy_short(rsi,rsf,evec,cut,mode)
!
  use sequen
  use iounit
  use energies
  use cutoffs
  use molecule
  use fyoc
  use atoms,ONLY: svte
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,k,mode,irs,frs,imol,sta,sto,incr,tpi
  RTYPE evec(MAXENERGYTERMS)
  logical cut
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
!
  irs = rsi
  frs = rsf
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
  sta = irs
  sto = frs
  tpi = 1
  incr = 1
  do j=sta,sto,incr
    call e_boundary_rs(j,evec_thr(:,tpi),mode)
  end do
!
! bonded terms may change but are restricted to a truly local effect (note -1/+1-padding, though)
  imol = molofrs(rsi)
  sta = max(rsmol(imol,1),rsi-1)
  sto = min(rsmol(imol,2),rsf+1)
  tpi = 1
  incr = 1
  do j=sta,sto,incr
    if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:,tpi))
    if (disulf(j).gt.0) then
      if ((disulf(j).lt.sta).OR.(disulf(j).gt.sto)) then
        if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:,tpi))
        if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:,tpi))
        if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:,tpi))
        if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:,tpi))
        if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:,tpi))
      end if
    end if
  end do
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do j=irs,frs
          call grd_respairs(j)
        end do
      else if (use_rescrit.EQV..true.) then
        do j=irs,frs
          call respairs(j)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in ujcr_energy_sho&
 &rt(...).'
        call fexit()
      end if
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
      do k=sta,sto,incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(max(j,k),min(j,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
            end if
          end do
        end if
      end do
      do j=irs,frs
        call clear_rsp2(j)
      end do
    else
      tpi = 1
      sta = 1
      sto = nseq
      incr = 1
      do k=sta,sto,incr ! 1nseq
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp(evec_thr(:,tpi),min(j,k),max(j,k),svte,cut)
          end do
        end if
      end do
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for pivot
! moves (only relevant terms are computed ...)
!
subroutine ujcr_energy_long(rsi,rsf,evec,cut,mode)
!
  use iounit
  use sequen
  use cutoffs
  use energies
  use molecule
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rsi,rsf,j,mode,i,k,l,rs,nrsl,rsl(nseq),frs,irs,sta(2),sto(2),incr,tpi,azero,imol
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
  imol = molofrs(rsi)
!
  frs = rsf
  irs = rsi
!
  azero = 0
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(azero)
  end if
!
  nrsl = 0
!
  if (use_TOR.EQV..true.) then
    sta(1) = max(rsmol(imol,1),rsi-1)
    sto(1) = min(rsmol(imol,2),rsf+1)
    do i=1,n_crosslinks
      if (((crosslink(i)%rsnrs(1).lt.rsi).OR.(crosslink(i)%rsnrs(1).gt.rsf)).AND.&
 &        ((crosslink(i)%rsnrs(2).lt.rsi).OR.(crosslink(i)%rsnrs(2).gt.rsf))) cycle
      if (((crosslink(i)%rsnrs(1).gt.sto(1)).OR.(crosslink(i)%rsnrs(1).lt.sta(1))).AND.(par_TOR2(crosslink(i)%rsnrs(1)).gt.0)) then
        call en_torrs(evec_thr(:,1),crosslink(i)%rsnrs(1))
      end if
      if (((crosslink(i)%rsnrs(2).gt.sto(1)).OR.(crosslink(i)%rsnrs(2).lt.sta(1))).AND.(par_TOR2(crosslink(i)%rsnrs(2)).gt.0)) then
        call en_torrs(evec_thr(:,1),crosslink(i)%rsnrs(2))
      end if
    end do
    do j=sta(1),sto(1)
      if (par_TOR2(j).gt.0) call en_torrs(evec_thr(:,1),j)
    end do
  end if
!
  if (use_IMPSOLV.EQV..true.) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if (cut.EQV..true.) then
      inrsl(:) = .false.
      do rs=1,irs-1
       if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      do rs=irs,frs
        nrsl = nrsl + 1
        rsl(nrsl) = rs
        inrsl(rs) = .true.
      end do
      do rs=frs+1,nseq
        if (rs_vec(rs).eq.1) then
          nrsl = nrsl + 1
          rsl(nrsl) = rs
          inrsl(rs) = .true.
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in ujcr_energy_l&
 &ong(...).'
        call fexit()
      end if
!
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr ! 1,nseq
!         we solve the issue of double counting by using the inverse map
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
!       and deal with the rest here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr ! i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      do i=sta(1),sto(1),incr ! 1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
!
    else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      do j=sta(1),sto(1),incr ! 1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
    end if
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if (cut.EQV..true.) then
      if (use_mcgrid.EQV..true.) then
        do i=irs,frs
          call grd_respairs(i)
        end do
      else if (use_rescrit.EQV..true.) then
        do i=irs,frs
          call respairs(i)
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in ujcr_energy_l&
 &ong(...).'
        call fexit()
      end if
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      do k=sta(1),sto(1),incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            if (rsp_mat(min(j,k),max(j,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
            end if
          end do
        end if
      end do
      do j=irs,frs
        call clear_rsp2(j)
      end do
    else
      tpi = 1
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      do k=sta(1),sto(1),incr
!       non-moving parts to moving parts
        if ((k.lt.irs).OR.(k.gt.frs)) then
          do j=irs,frs
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        else if ((k.ge.rsi).AND.(k.le.rsf)) then
!         intra CR-stretch
          do j=k,rsf
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),cut)
          end do
        end if
      end do
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a simple wrapper routine that handles short-range energy calculations for rigid body
! moves (only relevant terms are computed ...)
!
subroutine rigid_energy_short(imol,evec,cut,mode,tpi2)
!
  use sequen
  use iounit
  use energies
  use molecule
  use cutoffs
  use fyoc
  use atoms,ONLY: svte
  use dssps, ONLY: pepmol
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: imol,mode,tpi2
  logical, INTENT(IN):: cut
!
  integer rs,j,aone,atwo,azero,tpi,ix,emmode,irs,frs,rsl,rsh
  RTYPE evec(MAXENERGYTERMS)
  logical iflg
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),tpx,sta,sto,incr,jj
  logical OMP_IN_PARALLEL
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
#ifdef ENABLE_THREADS
  if ((tpi2.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, rigid_energy_short(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  evec_thr(:) = 0.0
!
  iflg = .true.
  azero = 0 
  aone = 1
  atwo = 2
  if (mode.eq.0) then
    ix = 2
#ifdef ENABLE_THREADS
    tpx = thrdat%maxn
#endif
  else
    ix = 1
#ifdef ENABLE_THREADS
    tpx = 0
#endif
  end if
  emmode = mode
  if (emmode.eq.0) emmode = 2
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    tpi = tpi2
  else
    tpi = 1
  end if
#else
  tpi = 1
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
#ifdef ENABLE_THREADS
      call init_svte_threads(azero,tpi2)
#else
      call init_svte(azero)
#endif
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
#ifdef ENABLE_THREADS
      call init_svte_threads(aone,tpi2)
#else
      call init_svte(aone)
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  frs = rsmol(imol,2)
  irs = rsmol(imol,1)
  if ((ideal_run.EQV..false.).AND.(cut.EQV..true.)) then
    rsl = irs
    rsh = frs
    call respairs_new(rsl,rsh,ix,irs,frs,azero,iflg,aone,tpi2)
  end if
!
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi2.gt.0) then
        call get_thread_loop_bounds(ix,aone,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%srnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp(evec_thr(:),thr_rsp_nbl(jj)%sr(j,1),thr_rsp_nbl(jj)%sr(j,2),thr_svte(:,tpi),use_cutoffs)
          end do
        end do
      else
#endif
      do j=1,rsp_nbl(ix)%srnrs
        call Ven_rsp(evec_thr(:),rsp_nbl(ix)%sr(j,1),rsp_nbl(ix)%sr(j,2),svte,use_cutoffs)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#else
      do j=irs,frs
        do rs=1,irs-1
          call Ven_rsp(evec_thr(:),rs,j,svte,cut)
        end do
        do rs=frs+1,nseq
          call Ven_rsp(evec_thr(:),j,rs,svte,cut)
        end do
      end do
#endif
    end if
!
  else
!   do nothing
  end if
!
  if (n_crosslinks.gt.0) then ! intermolecular only
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=irs,frs
      if (disulf(j).gt.0) then
        if ((disulf(j).gt.frs).OR.(disulf(j).lt.irs)) then
          if (use_BOND(1).EQV..true.) call en_bonds(j,evec_thr(:))
          if (use_BOND(2).EQV..true.) call en_angles(j,evec_thr(:))
          if (use_BOND(3).EQV..true.) call en_impropers(j,evec_thr(:))
          if (use_BOND(4).EQV..true.) call en_torsions(j,evec_thr(:))
          if (use_BOND(5).EQV..true.) call en_cmap(j,evec_thr(:))
          if (use_BOND(1).EQV..true.) call en_bonds(disulf(j),evec_thr(:))
          if (use_BOND(2).EQV..true.) call en_angles(disulf(j),evec_thr(:))
          if (use_BOND(3).EQV..true.) call en_impropers(disulf(j),evec_thr(:))
          if (use_BOND(4).EQV..true.) call en_torsions(disulf(j),evec_thr(:))
          if (use_BOND(5).EQV..true.) call en_cmap(disulf(j),evec_thr(:))
        end if
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
!
  if (use_DSSP.EQV..true.) then
    if (pepmol(imol).gt.0) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call en_dssp_gl(evec_thr(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
  end if

! boundary term can involve a fair number of terms
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    sta = irs+tpi2-1
    incr = thrdat%maxn
  else
    sta = irs
    incr = 1
  end if
  do j=sta,frs,incr
#else
  do j=irs,frs
#endif
    call e_boundary_rs(j,evec_thr(:),mode)
  end do
!
! a few bias terms are calculated prior to move evaluation
  if (use_DREST.EQV..true.) call edrest(evec_thr(:),tpi2) ! threaded and safe
!
  if (use_EMICRO.EQV..true.) call en_emicro_gl(evec_thr(:),emmode,tpi2) ! threaded and safe
!
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV and populate rs_vec and rs_vec2
#ifdef ENABLE_THREADS
      call init_svte_threads(atwo,tpi2)
#else
      call init_svte(atwo)
#endif
    end if
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!--------------------------------------------------------------------------
!
! the same for clusters of molecules
!
subroutine clurb_energy_short(imols,mmol,fmols,evec,cut,mode)
!
  use sequen
  use iounit
  use energies
  use molecule
  use cutoffs
  use system
  use fyoc
  use atoms,ONLY: svte
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer rs,imol,j,mode,i,mmol,imols(nmol),fmols(nmol),k,sta,sto,incr,tpi,rsl(nseq),nrsl
  RTYPE evec(MAXENERGYTERMS)
  logical cut
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
!
! the boundary terms depend only on coordinates -> can be taken care of immediately
  do i=1,mmol 
    imol = imols(i)
    do j=rsmol(imol,1),rsmol(imol,2)
      call e_boundary_rs(j,evec,mode)
    end do
  end do
! be safe
  tpi = 1
  do j=1,n_crosslinks
    if (use_BOND(1).EQV..true.) call en_bonds(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(crosslink(j)%rsnrs(1),evec_thr(:,tpi))
    if (use_BOND(1).EQV..true.) call en_bonds(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(2).EQV..true.) call en_angles(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(3).EQV..true.) call en_impropers(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(4).EQV..true.) call en_torsions(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
    if (use_BOND(5).EQV..true.) call en_cmap(crosslink(j)%rsnrs(2),evec_thr(:,tpi))
  end do
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     zero out svte and svbu (temporary arrays for changes in SAV)
      call init_svte(0)
    else if (mode.eq.1) then
!     store svte (the changes induced by the sampling) in svbu, zero out svte again
      call init_svte(1)
    end if
  end if
!
! let's only do all this setup work if necessary
  if ((use_IMPSOLV.EQV..true.).OR.(use_WCA.EQV..true.).OR.(use_IPP.EQV..true.).OR.&
 &    (use_attLJ.EQV..true.).OR.(use_CORR.EQV..true.).OR.(use_FEGS(1).EQV..true.).OR.&
 &    (use_FEGS(3).EQV..true.)) then
!
!   assemble list (makes code below easier)
    nrsl = 0
    do imol=1,nmol
      if (fmols(imol).eq.1) then
        if (bnd_type.eq.1) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
          end do
        end if
      else
        do rs=rsmol(imol,1),rsmol(imol,2)
          nrsl = nrsl + 1
          rsl(nrsl) = rs
        end do
      end if
    end do
!
    if ((cut.EQV..true.).AND.((use_rescrit.EQV..true.).OR.&
 &         (use_mcgrid.EQV..true.))) then
      if (use_mcgrid.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call grd_respairs(rs)
          end do
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call respairs(rs)
          end do
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in clurb_energy_sh&
 &ort(...).'
        call fexit()
      end if
!
      sta = 1
      sto = nrsl
      incr = 1
      tpi = 1
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta,sto,incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            if (rsp_mat(max(rs,k),min(rs,k)).ne.0) then
              call Ven_rsp(evec_thr(:,tpi),min(rs,k),max(rs,k),svte,use_cutoffs)
            end if
          end do
        end do
      end do
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          call clear_rsp2(rs)
        end do
      end do
    else
!
      sta = 1
      sto = nrsl
      incr = 1
      tpi = 1
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta,sto,incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            call Ven_rsp(evec_thr(:,tpi),min(rs,k),max(rs,k),svte,use_cutoffs)
          end do
        end do
      end do
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.1) then
!     get the difference in SAV
      call init_svte(2)
!     now use that difference to determine which dependent interactions we have to re-compute
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! a simple wrapper routine that handles long-range energy calculations for rigid body
! moves (only relevant terms are computed ...)
!
subroutine rigid_energy_long(imol,evec,cut,mode,tpi)
!
  use iounit
  use sequen
  use energies
  use cutoffs
  use molecule, ONLY: rsmol
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,mode,imol
  logical, INTENT(IN):: cut
!
  integer j,i,k,irs,frs,ix,azero,afive
  logical iflg
  RTYPE evec(MAXENERGYTERMS)
#ifdef ENABLE_THREADS
  integer stas(2),stos(2),sta,sto,atwo,athree,jj
#endif
  RTYPE evec_thr(MAXENERGYTERMS)
!
  evec_thr(:) = 0.0
!
  iflg = .false.
  azero = 0
  afive = 5
#ifdef ENABLE_THREADS
  atwo = 2
  athree = 3
#endif
  if (mode.eq.1) then
    ix = 2
  else
    ix = 1
  end if
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(tpi)
  end if
!     
  frs = rsmol(imol,2)
  irs = rsmol(imol,1)
!
  if (use_TOR.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do i=1,n_crosslinks
      if (((crosslink(i)%rsnrs(1).lt.irs).OR.(crosslink(i)%rsnrs(1).gt.frs)).AND.&
 &        ((crosslink(i)%rsnrs(2).lt.irs).OR.(crosslink(i)%rsnrs(2).gt.frs))) cycle
      if (par_TOR2(crosslink(i)%rsnrs(1)).gt.0) then
        call en_torrs(evec_thr(:),crosslink(i)%rsnrs(1))
      end if
      if (par_TOR2(crosslink(i)%rsnrs(2)).gt.0) then
        call en_torrs(evec_thr(:),crosslink(i)%rsnrs(2))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
! 
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.))) then
!
    if (cut.EQV..true.) then
      if ((use_POLAR.EQV..true.).AND.(scrq_model.ne.4)) then
        j = 0
        call respairs_append(ix,irs,frs,azero,j,tpi)
      end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
        do j=thr_limits(3,tpi),thr_limits(4,tpi)
          if ((rs_vec(j).eq.1).OR.(rs_vec2(j).eq.1)) call en_freesolv(evec_thr(:),j)
        end do
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
!
      do i=1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:),i)
        end if
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=irs,frs
        do j=1,irs-1
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
        do j=i,frs
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=1,irs-1
        do j=i,irs-1
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=frs+1,nseq
        do j=i,nseq
          if ((rs_vec(i).le.0).AND.(rs_vec(j).le.0)) cycle
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
      do i=1,nseq
        call en_freesolv(evec_thr(:),i)
      end do
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      sta = thr_limits(3,tpi)
      sto = thr_limits(4,tpi)
    else
      sta = 1
      sto = nseq
    end if
    do i=sta,sto
#else
    do i=1,nseq
#endif
      if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
        call en_freesolv(evec_thr(:),i)
      end if
    end do
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
    if (cut.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call get_thread_loop_bounds(ix,atwo,stas,stos,tpi)
        do jj=stas(1),stos(1)
          sta = 1
          sto = thr_rsp_nbl(jj)%trnrs
          if (jj.eq.stas(1)) sta = stas(2)
          if (jj.eq.stos(1)) sto = stos(2)
          do j=sta,sto
            call Ven_rsp_long(evec_thr(:),thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),use_cutoffs)
          end do
        end do
        if (use_POLAR.EQV..true.) then
          call get_thread_loop_bounds(ix,athree,stas,stos,tpi)
          do jj=stas(1),stos(1)
            sta = 1
            sto = thr_rsp_nbl(jj)%lrnrs
            if (jj.eq.stas(1)) sta = stas(2)
            if (jj.eq.stos(1)) sto = stos(2)
            do j=sta,sto
              call Ven_rsp_lrel(evec_thr(:),thr_rsp_nbl(jj)%lr(j,1),thr_rsp_nbl(jj)%lr(j,2))
            end do
          end do
        end if
      else
#endif
      do i=1,rsp_nbl(ix)%trnrs
        j = rsp_nbl(ix)%tr(i,1)
        k = rsp_nbl(ix)%tr(i,2)
        call Ven_rsp_long(evec_thr(:),j,k,use_cutoffs)
      end do
      do i=1,rsp_nbl(ix)%lrnrs
        j = rsp_nbl(ix)%lr(i,1)
        k = rsp_nbl(ix)%lr(i,2)
        call Ven_rsp_lrel(evec_thr(:),j,k)
      end do
#ifdef ENABLE_THREADS
      end if
#endif
    else
#ifdef ENABLE_THREADS
      call fexit()
#endif
      do i=irs,frs
        do j=1,irs-1
          call Ven_rsp_long(evec_thr(:),j,i,cut)
        end do
        do j=frs+1,nseq
          call Ven_rsp_long(evec_thr(:),i,j,cut)
        end do
      end do
    end if
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:)
!$OMP BARRIER
!$OMP SINGLE
  evec(:) = evec(:) + sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2)
!$OMP END SINGLE NOWAIT
  else
#endif
  evec(:) = evec(:) + evec_thr(:)
#ifdef ENABLE_THREADS
  end if
#endif
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
#ifdef ENABLE_THREADS
      call init_svte_threads(afive,tpi)
#else
      call init_svte(afive)
#endif
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        rs_vec(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
        rs_vec2(thr_limits(3,tpi):thr_limits(4,tpi)) = 0
      else
        rs_vec(:) = 0
        rs_vec2(:) = 0
      end if
#else
      rs_vec(:) = 0
      rs_vec2(:) = 0
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!------------------------------------------------------------------------
!
! the same for clusters of molecules
!
subroutine clurb_energy_long(imols,mmol,fmols,evec,cut,mode)
!
  use iounit
  use sequen
  use energies
  use molecule
  use cutoffs
  use system
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer imol,rs,j,mode,i,k,l,rsl(nseq),nrsl,sta(2),sto(2),incr,tpi,azero
  integer mmol,imols(nmol),fmols(nmol)
  RTYPE evec(MAXENERGYTERMS)
  logical cut,inrsl(nseq)
  RTYPE evec_thr(MAXENERGYTERMS,2)
!
  evec_thr(:,1) = 0.0
!
  azero = 0
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    call setup_scrqs2(azero)
  end if
!
  if (use_TOR.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do i=1,n_crosslinks
      if (par_TOR2(crosslink(i)%rsnrs(1)).gt.0) then
        call en_torrs(evec_thr(:,1),crosslink(i)%rsnrs(1))
      end if
      if (par_TOR2(crosslink(i)%rsnrs(2)).gt.0) then
        call en_torrs(evec_thr(:,1),crosslink(i)%rsnrs(2))
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
!
  nrsl = 0
!
  if ((use_IMPSOLV.EQV..true.).AND.((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.))) then
!
!   the problematic part is that the interaction of two residues can
!   indeed be changed by a third residue, which is changing conformation.
!   it is henceforth not safe to compute just the terms that change
!   by virtue of the conformational change.
!   to detect the implicitly affected terms we're using the temporary SAV
!   array (i.e., any atom whose SAV was affected by the conformational
!   change potentially interacts differently with all other atoms now)
!   through the rs_vec-array which is created during the last init_svte-call
!   in the short-range energy computation
!   note, however, that just like for every other interaction, the residue
!   based cutoff criterion still has to be fulfilled
!
    if ((cut.EQV..true.).AND.((use_rescrit.EQV..true.).OR.&
 &       (use_mcgrid.EQV..true.))) then
      inrsl(:) = .false.
      do imol=1,nmol
        if (fmols(imol).eq.1) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
            inrsl(rs) = .true.
          end do
        else
          do rs=rsmol(imol,1),rsmol(imol,2)
            if (rs_vec(rs).eq.1) then
              nrsl = nrsl + 1
              rsl(nrsl) = rs
              inrsl(rs) = .true.
            end if
          end do
        end if
      end do
      if (use_mcgrid.EQV..true.) then
        do i=1,nrsl
          call grd_respairs(rsl(i))
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,nrsl
          call respairs(rsl(i))
        end do
      end if
!
      sta(1) = 1
      sto(1) = nseq
      sto(2) = nrsl
      incr = 1
      tpi = 1
      do i=1,nrsl
        j = rsl(i)
        do k=sta(1),sto(1),incr !1,nseq
!         in order to avoid double counting, we use the inverse map to skip out if both residues are in rsl
          if (inrsl(k).EQV..true.) cycle
!
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
!       and handle the interactions within rsl here
        sta(2) = i + tpi - 1
        do l=sta(2),sto(2),incr !i,nrsl
          k = rsl(l)
          if (rsp_mat(min(j,k),max(j,k)).ne.0) then
            call Ven_rsp_long(evec_thr(:,tpi),min(j,k),max(j,k),use_cutoffs)
          end if
        end do
      end do
      do i=sta(1),sto(1),incr !1,nseq
        if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
          call en_freesolv(evec_thr(:,tpi),i)
        end if
      end do
!
      do i=1,nrsl
        call clear_rsp2(rsl(i))
      end do
!
    else ! meaning no cutoffs
!     this a disaster of course, use of cutoffs essential when doing IMPSOLV
      sta(1) = 1
      sto(1) = nseq
      incr = 1
      tpi = 1
      do j=sta(1),sto(1),incr !1,nseq
        do k=j,nseq
          call Ven_rsp_long(evec_thr(:,tpi),j,k,cut)
        end do
        call en_freesolv(evec_thr(:,tpi),j)
      end do
    end if
!
  else if (use_IMPSOLV.EQV..true.) then
    sta(1) = 1
    sto(1) = nseq
    incr = 1
    tpi = 1
    do i=sta(1),sto(1),incr !1,nseq
      if ((rs_vec(i).eq.1).OR.(rs_vec2(i).eq.1)) then
        call en_freesolv(evec_thr(:,tpi),i)
      end if
    end do
!
  else if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.).OR.(use_FEGS(6).EQV..true.)) then
    if ((cut.EQV..true.).AND.(use_rescrit.EQV..true.).OR.&
 &                      (use_mcgrid.EQV..true.)) then
      if (use_mcgrid.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call grd_respairs(rs)
          end do
        end do
      else if (use_rescrit.EQV..true.) then
        do i=1,mmol 
          imol = imols(i)
          do rs=rsmol(imol,1),rsmol(imol,2)
            call respairs(rs)
          end do
        end do
      else
        write(ilog,*) 'Fatal. Undefined cutoff mode in clurb_energy_&
 &long(...).'
        call fexit()
      end if
!     we'll re-interpret and re-utilize rsl here
!     this includes the potential pitfall of having to recompute internal cluster energies
!     if PBC are used
      nrsl = 0
      do imol=1,nmol
        if (fmols(imol).eq.1) then
          if (bnd_type.eq.1) then
            do rs=rsmol(imol,1),rsmol(imol,2)
              nrsl = nrsl + 1
              rsl(nrsl) = rs
            end do
          end if
        else
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
          end do
        end if
      end do
      sta(1) = 1
      sto(1) = nrsl
      incr = 1
      tpi = 1
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta(1),sto(1),incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            if (rsp_mat(min(rs,k),max(rs,k)).ne.0) then
              call Ven_rsp_long(evec_thr(:,tpi),min(rs,k),max(rs,k),use_cutoffs)
            end if
          end do
        end do
      end do
!
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          call clear_rsp2(rs)
        end do
      end do
    else
!     ditto
      nrsl = 0
      do imol=1,nmol
        if (fmols(imol).eq.1) then
          if (bnd_type.eq.1) then
            do rs=rsmol(imol,1),rsmol(imol,2)
              nrsl = nrsl + 1
              rsl(nrsl) = rs
            end do
          end if
        else
          do rs=rsmol(imol,1),rsmol(imol,2)
            nrsl = nrsl + 1
            rsl(nrsl) = rs
          end do
        end if
      end do
      sta(1) = 1
      sto(1) = nrsl
      incr = 1
      tpi = 1
      do i=1,mmol 
        imol = imols(i)
        do rs=rsmol(imol,1),rsmol(imol,2)
          do j=sta(1),sto(1),incr
            k = rsl(j)
            if (fmols(molofrs(k)).eq.1) then
              if (imol.ge.molofrs(k)) cycle
            end if
            call Ven_rsp_long(evec_thr(:,tpi),min(rs,k),max(rs,k),use_cutoffs)
          end do
        end do
      end do
    end if
!
  else
!   do nothing
  end if
!
  evec(:) = evec(:) + evec_thr(:,1)
!
  if (use_IMPSOLV.EQV..true.) then
    if (mode.eq.0) then
!     now update the atsav-array using svte (to be ready for posterior long-range calculation)
      call init_svte(5)
    else if (mode.eq.1) then
!     clear out the residue-svte-monitoring-vector
      do i=1,nseq
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end do
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
!             #################################
!             #                               #
!             # THIRD: OTHER ROUTINES         #
!             #                               #
!             #################################
!
!-------------------------------------------------------------------------------
!
! this simple wrapper routine sets the system's Hamiltonian to different
! conditions allow swap moves in REMC ...
! note that lamenergy is implictily threaded (if so compiled) through the call to the
! main time-consuming function energy3
!
subroutine lamenergy(rve,fve,tpi)
!
  use iounit
  use system
  use energies
  use units
  use mpistuff
  use molecule
  use dssps
  use ems
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  RTYPE rve(mpi_nodes),fve(mpi_nodes)
#ifdef ENABLE_MPI
  integer i,j,imol,which,i_start,i_end
  RTYPE evec(MAXENERGYTERMS)
  RTYPE vbu(MAXREDIMS),eva,evaref,dum,dum2
  logical needsavup,needemup,badflg(MAXREDIMS),atrue,afalse
  integer vbui(MAXREDIMS)
#endif
!
#ifdef ENABLE_MPI
!
  needsavup = .false.
  needemup = .false.
  atrue = .true.
  afalse = .false.
  dum2 = 0.0
  evaref = esave
  evec(:) = esterms(:)
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  rve(:) = 0.0
  fve(:) = 0.0
!
  do i=1,re_conddim
    if (re_types(i).eq.1) vbu(i)=kelvin
    if (re_types(i).eq.2) vbu(i)=scale_IPP
    if (re_types(i).eq.3) vbu(i)=scale_attLJ
    if (re_types(i).eq.4) vbu(i)=scale_WCA
    if (re_types(i).eq.5) vbu(i)=scale_POLAR
    if (re_types(i).eq.6) vbu(i)=scale_IMPSOLV
    if (re_types(i).eq.7) vbu(i)=par_IMPSOLV(2)
    if (re_types(i).eq.8) vbu(i)=scale_TOR
    if (re_types(i).eq.9) vbu(i)=scale_ZSEC
    if (re_types(i).eq.10) vbu(i)=par_ZSEC(1)
    if (re_types(i).eq.11) vbu(i)=par_ZSEC(3)
    if (re_types(i).eq.12) vbui(i)=scrq_model
    if (re_types(i).eq.13) vbu(i)=par_IMPSOLV(3)
    if (re_types(i).eq.14) vbu(i)=par_IMPSOLV(4)
    if (re_types(i).eq.15) vbu(i)=par_IMPSOLV(6)
    if (re_types(i).eq.16) vbu(i)=par_IMPSOLV(7)
    if (re_types(i).eq.17) vbu(i)=par_IMPSOLV(8)
    if (re_types(i).eq.18) vbui(i)=i_sqm
    if (re_types(i).eq.19) vbu(i)=par_IMPSOLV(9)
    if (re_types(i).eq.20) vbu(i)=scale_FEGS(1)
    if (re_types(i).eq.21) vbu(i)=scale_FEGS(3)
    if (re_types(i).eq.22) vbu(i)=scale_FEGS(6)
    if (re_types(i).eq.23) vbu(i)=scale_TABUL
    if (re_types(i).eq.24) vbu(i)=scale_POLY
    if (re_types(i).eq.25) vbu(i)=scale_DREST
    if (re_types(i).eq.26) vbu(i)=scale_FEGS(15)
    if (re_types(i).eq.27) vbu(i)=scale_FEGS(16)
    if (re_types(i).eq.28) vbu(i)=scale_FEGS(17)
    if (re_types(i).eq.29) vbu(i)=scale_FEGS(18)
    if (re_types(i).eq.30) vbu(i)=par_DSSP(9)
    if (re_types(i).eq.31) vbu(i)=par_DSSP(7)
!   do nothing for 32
    if (re_types(i).eq.33) vbu(i)=scale_EMICRO
    if (re_types(i).eq.34) vbu(i)=emthreshdensity
  end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
! note that MPI_REMaster(...) temporarily changes reol_all for the
! the actual exchange overlap calculations if need be (i.e., if re_nbmode is 1)
  if (reol_all.EQV..true.) then
    i_start = 1
    i_end = re_conditions
  else
    i_start = max((myrank+1)-1,1)
    i_end = min((myrank+1)+1,re_conditions)
  end if
!
  do i=i_start,i_end
!   set the different conditions
!   remember that we don't change the use_XX-flags
!   this implies, however, that for scale_XX = 0.0, we compute
!   a bunch of terms all multiplied by 0.0. in order to preserve
!   this information (to compute derivatives with respect to scale_XX),
!   we therefore use a little detour (set to 1.0, subtract out)
    eva = 0.0
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    do j=1,re_conddim
      badflg(j) = .false.
      if (re_types(j).eq.1) then
        kelvin = re_mat(i,j)
        invtemp = 1.0/(gasconst*kelvin)
      else if (re_types(j).eq.2) then
        scale_IPP   = re_mat(i,j)
        if (scale_IPP.le.0.0) then
          badflg(j) = .true.
          scale_IPP = 1.0
        end if
      else if (re_types(j).eq.3) then
        scale_attLJ = re_mat(i,j)
        if (scale_attLJ.le.0.0) then
          badflg(j) = .true.
          scale_attLJ = 1.0
        end if
      else if (re_types(j).eq.4) then
        scale_WCA   = re_mat(i,j)
        if (scale_WCA.le.0.0) then
          badflg(j) = .true.
          scale_WCA = 1.0
        end if
      else if (re_types(j).eq.5) then
        scale_POLAR = re_mat(i,j)
        if (scale_POLAR.le.0.0) then
          badflg(j) = .true.
          scale_POLAR = 1.0
        end if
      else if (re_types(j).eq.6) then
        scale_IMPSOLV = re_mat(i,j)
        if (scale_IMPSOLV.le.0.0) then
          badflg(j) = .true.
          scale_IMPSOLV = 1.0
        end if
      else if (re_types(j).eq.7) then
!       note that the relevance of this relies entirely on whether use_IMPSOLV is true
        par_IMPSOLV(2) = re_mat(i,j)
        if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
          coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
        else
          coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
        end if
      else if (re_types(j).eq.8) then
        scale_TOR = re_mat(i,j)
        if (scale_TOR.le.0.0) then
          badflg(j) = .true.
          scale_TOR = 1.0
        end if
      else if (re_types(j).eq.9) then
        scale_ZSEC = re_mat(i,j)
        if (scale_ZSEC.le.0.0) then
          badflg(j) = .true.
          scale_ZSEC = 1.0
        end if
      else if (re_types(j).eq.10) then
!       note that the relevance of this hinges on scale_ZSEC
        par_ZSEC(1) = re_mat(i,j)
      else if (re_types(j).eq.11) then
!       note that the relevance of this hinges on scale_ZSEC
        par_ZSEC(3) = re_mat(i,j)
      else if (re_types(j).eq.12) then
!       note that the relevance of this hinges on scale_IMPSOLV
        scrq_model = nint(re_mat(i,j))
        if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
          coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
        else
          coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
        end if
      else if (re_types(j).eq.13) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(3) = re_mat(i,j)
      else if (re_types(j).eq.14) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(4) = re_mat(i,j)
      else if (re_types(j).eq.15) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(6) = re_mat(i,j)
      else if (re_types(j).eq.16) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(7) = re_mat(i,j)
      else if (re_types(j).eq.17) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        par_IMPSOLV(8) = 1./re_mat(i,j)
      else if (re_types(j).eq.18) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        i_sqm = nint(re_mat(i,j))
      else if (re_types(j).eq.19) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        par_IMPSOLV(9) = re_mat(i,j)
      else if (re_types(j).eq.20) then
        scale_FEGS(1) = re_mat(i,j)
        call setup_parFEG(1)
!        if (scale_FEGS(1).le.0.0) then
!          badflg(j) = .true.
!          scale_FEGS(1) = 1.0
!        end if
      else if (re_types(j).eq.21) then
        scale_FEGS(3) = re_mat(i,j)
        call setup_parFEG(3)
!        if (scale_FEGS(3).le.0.0) then
!          badflg(j) = .true.
!          scale_FEGS(3) = 1.0
!        end if
      else if (re_types(j).eq.22) then
        scale_FEGS(6) = re_mat(i,j)
        call setup_parFEG(6)
!        if (scale_FEGS(6).le.0.0) then
!          badflg(j) = .true.
!          scale_FEGS(6) = 1.0
!        end if
      else if (re_types(j).eq.23) then
        scale_TABUL = re_mat(i,j)
        if (scale_TABUL.le.0.0) then
          badflg(j) = .true.
          scale_TABUL = 1.0
        end if
      else if (re_types(j).eq.24) then
        scale_POLY = re_mat(i,j)
        if (scale_POLY.le.0.0) then
          badflg(j) = .true.
          scale_POLY = 1.0
        end if
      else if (re_types(j).eq.25) then
        scale_DREST = re_mat(i,j)
        if (scale_DREST.le.0.0) then
          badflg(j) = .true.
          scale_DREST = 1.0
        end if
      else if (re_types(j).eq.26) then
        scale_FEGS(15) = re_mat(i,j)
      else if (re_types(j).eq.27) then
        scale_FEGS(16) = re_mat(i,j)
      else if (re_types(j).eq.28) then
        scale_FEGS(17) = re_mat(i,j)
      else if (re_types(j).eq.29) then
        scale_FEGS(18) = re_mat(i,j)
      else if (re_types(j).eq.30) then
!       note that the relevance of this hinges on scale_DSSP
        par_DSSP(9) = re_mat(i,j)
      else if (re_types(j).eq.31) then
!       note that the relevance of this hinges on scale_DSSP
        par_DSSP(7) = re_mat(i,j)
!     do nothing for 32
      else if (re_types(j).eq.33) then
        scale_EMICRO = re_mat(i,j)
        if (scale_EMICRO.le.0.0) then
          badflg(j) = .true.
          scale_EMICRO = 1.0
        end if
      else if (re_types(j).eq.34) then
        emthreshdensity = re_mat(i,j)
        needemup = .true.
      end if
!
    end do
    if (needsavup.EQV..true.) call absinth_savprm()
    if (needemup.EQV..true.) then
      call scale_emmap(emthreshdensity,dum,dum2)
      call precompute_diff_emmap()
    end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!   now compute the current energy at this condition with
!   our structure
    if (Tonly.EQV..true.) then
      eva = evaref
    else
      if (i.eq.i_start) then
        call energy3(esterms,atrue,esave,tpi)
        eva = esave
      else
        call energy3(esterms,afalse,esave,tpi)
        eva = esave
      end if
    end if
!   the computation of derivatives is really only meaningful when
!   noTI is false
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    fve(i) = 0.0
    do j=1,re_conddim
      if (badflg(j).EQV..true.) then
        if (re_types(j).eq.2) then
          fve(i) = fve(i) + esterms(1)
          esterms(1) = 0.0
        else if (re_types(j).eq.3) then
          fve(i) = fve(i) + esterms(3)
          esterms(3) = 0.0
        else if (re_types(j).eq.4) then
          fve(i) = fve(i) + esterms(5)
          esterms(5) = 0.0
        else if (re_types(j).eq.5) then
          fve(i) = fve(i) + esterms(6)
          esterms(6) = 0.0
        else if (re_types(j).eq.6) then
          fve(i) = fve(i) + esterms(4)
          esterms(4) = 0.0
        else if (re_types(j).eq.8) then
          fve(i) = fve(i) + esterms(7)
          esterms(7) = 0.0
        else if (re_types(j).eq.9) then
          fve(i) = fve(i) + esterms(8)
          esterms(8) = 0.0
        else if (re_types(j).eq.23) then
          fve(i) = fve(i) + esterms(9)
          esterms(9) = 0.0
        else if (re_types(j).eq.24) then
          fve(i) = fve(i) + esterms(14)
          esterms(14) = 0.0
        else if (re_types(j).eq.25) then
          fve(i) = fve(i) + esterms(10)
          esterms(10) = 0.0
        else if (re_types(j).eq.33) then
          fve(i) = fve(i) + esterms(21)
          esterms(21) = 0.0
        end if
      else
        if (re_types(j).eq.2) then
          fve(i) = fve(i) + esterms(1)/scale_IPP
        else if (re_types(j).eq.3) then
          fve(i) = fve(i) + esterms(3)/scale_attLJ
        else if (re_types(j).eq.4) then
          fve(i) = fve(i) + esterms(5)/scale_WCA
        else if (re_types(j).eq.5) then
          fve(i) = fve(i) + esterms(6)/scale_POLAR
        else if (re_types(j).eq.6) then
          fve(i) = fve(i) + esterms(4)/scale_IMPSOLV
        else if (re_types(j).eq.8) then
          fve(i) = fve(i) + esterms(7)/scale_TOR
        else if (re_types(j).eq.9) then
          fve(i) = fve(i) + esterms(8)/scale_ZSEC
        else if (re_types(j).eq.10) then
          which = 1
          do imol=1,nmol
            call der_zsec_gl(imol,fve(i),which)
          end do
        else if (re_types(j).eq.11) then
          which = 2
          do imol=1,nmol
            call der_zsec_gl(imol,fve(i),which)
          end do
        else if (re_types(j).eq.23) then
          fve(i) = fve(i) + esterms(9)/scale_TABUL
        else if (re_types(j).eq.24) then
          fve(i) = fve(i) + esterms(14)/scale_POLY
        else if (re_types(j).eq.25) then
          fve(i) = fve(i) + esterms(10)/scale_DREST
!        else if (re_types(j).eq.30) then
!          which = 2
!          do imol=1,nmol
!            call der_dssp_gl(imol,fve(i),which)
!          end do
!        else if (re_types(j).eq.31) then
!          which = 1
!          do imol=1,nmol
!            call der_dssp_gl(imol,fve(i),which)
!          end do
        else if (re_types(j).eq.33) then
          fve(i) = fve(i) + esterms(21)/scale_EMICRO
        end if
      end if
    end do
!   now recover the actual energy at this condition (see above)
    rve(i) = eva*invtemp
    needsavup = .false.
    needemup = .false.
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  do i=1,re_conddim
    if (re_types(i).eq.1) then 
      kelvin = vbu(i)
      invtemp = 1.0/(gasconst*kelvin)
    else if (re_types(i).eq.2) then
      scale_IPP   = vbu(i)
    else if (re_types(i).eq.3) then
      scale_attLJ = vbu(i)
    else if (re_types(i).eq.4) then
      scale_WCA   = vbu(i)
    else if (re_types(i).eq.5) then
      scale_POLAR = vbu(i)
    else if (re_types(i).eq.6) then
      scale_IMPSOLV = vbu(i)
    else if (re_types(i).eq.7) then
      par_IMPSOLV(2) = vbu(i)
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
      else
        coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
      end if
    else if (re_types(i).eq.8) then
      scale_TOR = vbu(i)
    else if (re_types(i).eq.9) then
      scale_ZSEC = vbu(i)
    else if (re_types(i).eq.10) then
      par_ZSEC(1) = vbu(i)
    else if (re_types(i).eq.11) then
      par_ZSEC(3) = vbu(i)
    else if (re_types(i).eq.12) then
      scrq_model = vbui(i)
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
      else
        coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
      end if
    else if (re_types(i).eq.13) then
      needsavup = .true.
      par_IMPSOLV(3) = vbu(i)
    else if (re_types(i).eq.14) then
      needsavup = .true.
      par_IMPSOLV(4) = vbu(i)
    else if (re_types(i).eq.15) then
      needsavup = .true.
      par_IMPSOLV(6) = vbu(i)
    else if (re_types(i).eq.16) then
      needsavup = .true.
      par_IMPSOLV(7) = vbu(i)
    else if (re_types(i).eq.17) then
      par_IMPSOLV(8) = vbu(i)
    else if (re_types(i).eq.18) then
      i_sqm = vbui(i)
    else if (re_types(i).eq.19) then
      par_IMPSOLV(9) = vbu(i)
    else if (re_types(i).eq.20) then
      scale_FEGS(1) = vbu(i)
      call setup_parFEG(1) 
    else if (re_types(i).eq.21) then
      scale_FEGS(3) = vbu(i)
      call setup_parFEG(3)
    else if (re_types(i).eq.22) then
      scale_FEGS(6) = vbu(i)
      call setup_parFEG(6)
    else if (re_types(i).eq.23) then
      scale_TABUL = vbu(i)
    else if (re_types(i).eq.24) then
      scale_POLY  = vbu(i)
    else if (re_types(i).eq.25) then
      scale_DREST = vbu(i)
    else if (re_types(i).eq.26) then
      scale_FEGS(15) = vbu(i)
    else if (re_types(i).eq.27) then
      scale_FEGS(16) = vbu(i)
    else if (re_types(i).eq.28) then
      scale_FEGS(17) = vbu(i)
    else if (re_types(i).eq.29) then
      scale_FEGS(18) = vbu(i)
    else if (re_types(i).eq.30) then
      par_DSSP(9) = vbu(i)
    else if (re_types(i).eq.31) then
      par_DSSP(7) = vbu(i)
!   do nothing for 32
    else if (re_types(i).eq.33) then
      scale_EMICRO = vbu(i)
    else if (re_types(i).eq.34) then
      emthreshdensity = vbu(i)
      needemup = .true.
    end if
  end do
  if (needsavup.EQV..true.) call absinth_savprm()
  if (needemup.EQV..true.) then
    call scale_emmap(emthreshdensity,dum,dum2)
    call precompute_diff_emmap()
  end if
  esave = evaref
  esterms(:) = evec(:)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
#else
!
  write(ilog,*) 'Fatal. Called lamenergy(...) in non MPI-calculation&
 &. This is most definitely a bug.'
  call fexit()
#endif
end
!
!-----------------------------------------------------------------------------
!
subroutine setup_parFEG(which)
!
  use energies
  use iounit 
!
  implicit none
!
  integer which
!
  if (which.eq.1) then
    if (scale_FEGS(1).ge.1.0) then
      par_FEG2(1) = scale_FEGS(1)
      par_FEG2(2) = 0.0
    else
      if (fegljmode.eq.1) then
        par_FEG2(1) = scale_FEGS(1)
        par_FEG2(2) = 0.0
      else if (fegljmode.eq.2) then
        par_FEG2(1) = scale_FEGS(1)**par_FEG2(8)
        par_FEG2(2) = par_FEG2(3)*(1.0-scale_FEGS(1)**par_FEG2(7))
      else if (fegljmode.eq.3) then
        par_FEG2(1) = (1.0-exp(-par_FEG2(8)*scale_FEGS(1)))/&
 &                                 (1.0-exp(-par_FEG2(8)))
        par_FEG2(2) = par_FEG2(3)*(1.0-scale_FEGS(1))**par_FEG2(7)
      else
        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-LJ (',&
 &fegljmode,').'
        call fexit()
      end if
    end if
  else if (which.eq.3) then
    if (scale_FEGS(3).ge.1.0) then
      par_FEG2(5) = scale_FEGS(3)
      par_FEG2(6) = 0.0
    else
      if (fegljmode.eq.1) then
        par_FEG2(5) = scale_FEGS(3)
        par_FEG2(6) = 0.0
      else if (fegljmode.eq.2) then
        par_FEG2(5) = scale_FEGS(3)**par_FEG2(8)
        par_FEG2(6) = par_FEG2(3)*(1.0-scale_FEGS(3)**par_FEG2(7))
      else if (fegljmode.eq.3) then
        par_FEG2(5) = (1.0-exp(-par_FEG2(8)*scale_FEGS(3)))/&
 &                                 (1.0-exp(-par_FEG2(8)))
        par_FEG2(6) = par_FEG2(3)*(1.0-scale_FEGS(3))**par_FEG2(7)
      else
        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-LJ (',&
 &fegljmode,').'
        call fexit()
      end if
    end if
  else if (which.eq.6) then
    if (scale_FEGS(6).ge.1.0) then
      par_FEG2(9) = scale_FEGS(6)
      par_FEG2(10) = 0.0
    else
      if (fegcbmode.eq.1) then
        par_FEG2(9) = scale_FEGS(6)
        par_FEG2(10) = 0.0
      else if (fegcbmode.eq.2) then
        par_FEG2(9) = scale_FEGS(6)**par_FEG2(12)
        par_FEG2(10) = par_FEG2(4)*(1.0-scale_FEGS(3)**par_FEG2(11))
      else
        write(ilog,*) 'Fatal. Got bad scaling mode for FEG-Cb (',&
 &fegcbmode,').'
        call fexit()
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called setup_parFEG(...) for unsupported p&
 &tential term (got ',which,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a simple wrapper to summarize torsional potentials applied to individual
! backbone torsions (handled through function en_torrs)
! currently NOT IN USE
!
function etorgl()
!
  use sequen
  use energies
!
  implicit none
!
  integer i
  RTYPE etorgl,ev(MAXENERGYTERMS)
!
  etorgl = 0.0d0
!
  do i = 1,nseq
    if (par_TOR2(i).gt.0) then
      call en_torrs(ev,i)
      etorgl = etorgl + ev(7)
    end if
  end do
!
end
!
!------------------------------------------------------------------------------------
!

