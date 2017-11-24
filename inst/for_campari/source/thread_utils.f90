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
!---------------------------------------------------------------------------
!
! a subroutine to set some initial limits for thread-based computations
!
#ifdef ENABLE_THREADS
!
subroutine threads_init()
!
  use sequen
  use atoms
  use threads
  use mcgrid
  use shakeetal
  use iounit
  use cutoffs
  use ewalds
  use molecule
  use system
  use energies
  use distrest
  use ems
  use movesets
  use forces
  use polypep
  use, INTRINSIC:: ISO_FORTRAN_ENV
!
  implicit none
!
  integer i,k,lastval(100),tpi,tpn,tpn2,shf2,stx(4),aone,afive
#ifdef ORACLE_FORTRAN
  integer(KIND=8) iacnt,lval2,tmps(2),kksum
#else
  integer(KIND=INT64) iacnt,lval2,tmps(2),kksum
#endif
  integer, ALLOCATABLE:: testmap(:)
!
  aone = 1
  afive = 5
  tpn = thrdat%maxn
  tpn2 = tpn
  shf2 = 0
  lastval(:) = 0
  lval2 = 0
#ifdef ORACLE_FORTRAN
  iacnt = int(n,KIND=8)*(int(n,KIND=8)-1_8)/2_8
#else
  iacnt = int(n,KIND=INT64)*(int(n,KIND=INT64)-1_INT64)/2_INT64
#endif
  do i=1,tpn
!   1-4 are static limits in the two most important hierarchy levels: atoms and residues
    thr_limits(1,i) = min(n,lastval(1)) + 1
    lastval(1) = nint(i*(1.0*n/(1.0*tpn)))
    thr_limits(2,i) = max(1,lastval(1))
    thr_limits(3,i) = min(nseq,lastval(2)) + 1
    lastval(2) = nint(i*(1.0*nseq/(1.0*tpn)))
    thr_limits(4,i) = max(1,lastval(2))
!   5,6 are dynamic bounds for grid-based NBL search
    thr_limits(5,i) = min(grid%pnts,lastval(3)) + 1
    lastval(3) = nint(i*(1.0*grid%pnts/(1.0*tpn)))
    thr_limits(6,i) = max(1,lastval(3))
    if (grid%pnts.le.0) thr_limits(6,i) = 0
!   7,8 are dynamic limits for NBL generation
    thr_limits(7,i) = min(nseq,lastval(4)) + 1
    lastval(4) = nint(i*(1.0*nseq/(1.0*tpn)))
    thr_limits(8,i) = max(1,lastval(4))
!   9,10 are dynamic limits for interactions (global or just IPP/ATTLJ/WCA/SAV) when no twin-range is considered
    tpn2 = tpn
!    if (thrdat%subnrs(1).gt.0) then 
!      tpn2 = thrdat%maxn - thrdat%subnrs(1)
!    end if
!    if (i.le.tpn2) then
!      thr_limits(9,i) = min(nseq,lastval(5)) + 1
!      lastval(5) = nint(i*(1.0*nseq/(1.0*tpn2)))
!      thr_limits(10,i) = max(1,lastval(5))
!    else
!      thr_limits(9,i) = 1
!      thr_limits(10,i) = 0
!    end if
    thr_limits(9,i) = min(nseq,lastval(5)) + 1
    lastval(5) = nint(i*(1.0*nseq/(1.0*tpn))) 
    thr_limits(10,i) = max(1,lastval(5))
!   11-26 are various limits for holonomic constraint solvers
    if ((settle_tip3ps.gt.0).AND.(i.gt.shf2)) then
      thr_limits(11,i) = min(settle_tip3ps,lastval(6)) + 1
      lastval(6) = nint((i-shf2)*(1.0*settle_tip3ps/(1.0*tpn2)))
      thr_limits(12,i) = max(1,lastval(6))
    else
      thr_limits(11,i) = 1
      thr_limits(12,i) = 0
    end if
    if ((settle_spcs.gt.0).AND.(i.gt.shf2)) then
      thr_limits(13,i) = min(settle_spcs,lastval(7)) + 1
      lastval(7) = nint((i-shf2)*(1.0*settle_spcs/(1.0*tpn2)))
      thr_limits(14,i) = max(1,lastval(7))
    else
      thr_limits(13,i) = 1
      thr_limits(14,i) = 0
    end if
    if ((settle_tip4ps.gt.0).AND.(i.gt.shf2)) then
      thr_limits(15,i) = min(settle_tip4ps,lastval(8)) + 1
      lastval(8) = nint((i-shf2)*(1.0*settle_tip4ps/(1.0*tpn2)))
      thr_limits(16,i) = max(1,lastval(8))
    else
      thr_limits(15,i) = 1
      thr_limits(16,i) = 0
    end if
    if ((settle_tip5ps.gt.0).AND.(i.gt.shf2)) then
      thr_limits(17,i) = min(settle_tip5ps,lastval(9)) + 1
      lastval(9) = nint((i-shf2)*(1.0*settle_tip5ps/(1.0*tpn2)))
      thr_limits(18,i) = max(1,lastval(9))
    else
      thr_limits(17,i) = 1
      thr_limits(18,i) = 0
    end if
    if ((settle_tip4pes.gt.0).AND.(i.gt.shf2)) then
      thr_limits(19,i) = min(settle_tip4pes,lastval(10)) + 1
      lastval(10) = nint((i-shf2)*(1.0*settle_tip4pes/(1.0*tpn2)))
      thr_limits(20,i) = max(1,lastval(10))
    else
      thr_limits(19,i) = 1
      thr_limits(20,i) = 0
    end if
    if ((settle_rest.gt.0).AND.(i.gt.shf2)) then
      thr_limits(21,i) = min(settle_rest,lastval(11)) + 1
      lastval(11) = nint((i-shf2)*(1.0*settle_rest/(1.0*tpn2)))
      thr_limits(22,i) = max(1,lastval(11))
    else
      thr_limits(21,i) = 1
      thr_limits(22,i) = 0
    end if
    thr_limits(23,i) = min(shake_cnt,lastval(12)) + 1
    lastval(12) = nint(i*(1.0*shake_cnt/(1.0*tpn)))
    thr_limits(24,i) = max(1,lastval(12))
    if (shake_cnt.le.0) thr_limits(24,i) = 0
    thr_limits(25,i) = min(cart_cons_grps,lastval(13)) + 1
    lastval(13) = nint(i*(1.0*cart_cons_grps/(1.0*tpn)))
    thr_limits(26,i) = max(1,lastval(13))
    if (cart_cons_grps.le.0) thr_limits(26,i) = 0
!   27-30 are various bounds for LR electrostatics treatments
    if ((lrel_md.eq.1).OR.(lrel_md.eq.3)) then
      thr_limits(27,i) = 1
      thr_limits(28,i) = 0
    else if (lrel_md.eq.2) then
      tpn2 = tpn
      if (thrdat%subnrs(1).gt.0) then
        tpn2 = thrdat%subnrs(1)
        if (i.le.tpn2) then
          thr_limits(27,i+tpn-thrdat%subnrs(1)) = min(nseq,lastval(14)) + 1
          lastval(14) = nint(i*(1.0*nseq/(1.0*tpn2)))
          thr_limits(28,i+tpn-thrdat%subnrs(1)) = max(1,lastval(14))
        else
          thr_limits(27,i-thrdat%subnrs(1)) = 1
          thr_limits(28,i-thrdat%subnrs(1)) = 0
        end if
      else 
        thr_limits(27,i) = min(nseq,lastval(14)) + 1
        lastval(14) = nint(i*(1.0*nseq/(1.0*tpn)))
        thr_limits(28,i) = max(1,lastval(14))
      end if
    else if ((lrel_md.eq.4).OR.(lrel_md.eq.5)) then
      tpn2 = tpn
      if (thrdat%subnrs(1).gt.0) then
        tpn2 = thrdat%subnrs(1)
        if ((i.le.tpn2).AND.(cglst%ncrs.gt.0)) then
          thr_limits(27,i+tpn-thrdat%subnrs(1)) = min(cglst%ncrs,lastval(14)) + 1
          lastval(14) = nint(i*(1.0*cglst%ncrs/(1.0*tpn2)))
          thr_limits(28,i+tpn-thrdat%subnrs(1)) = max(1,lastval(14))
        else
          thr_limits(27,i-thrdat%subnrs(1)) = 1
          thr_limits(28,i-thrdat%subnrs(1)) = 0
        end if
      else if (cglst%ncrs.gt.0) then
        thr_limits(27,i) = min(cglst%ncrs,lastval(14)) + 1
        lastval(14) = nint(i*(1.0*cglst%ncrs/(1.0*tpn)))
        thr_limits(28,i) = max(1,lastval(14))
      else
        thr_limits(27,i) = 1
        thr_limits(28,i) = 0
      end if
    end if
    if ((lrel_md.eq.2).AND.(ewald_mode.eq.1)) then
      if (ewald_mode.eq.2) then
        k = kdims(3,2)-kdims(3,1)+1 
      else
        k = kdims(3,2)
      end if
      if (use_dyn.EQV..true.) then
        if ((mod(k,tpn).ne.0).AND.(i.eq.1)) write(ilog,*) 'Warning. For better thread load balancing when using PME, consider &
 &having a number of mesh points in the z-dimension that is an integer multiple of the number of threads.'
      end if
      tpn2 = tpn
      if (thrdat%subnrs(1).gt.0) then 
        tpn2 = thrdat%subnrs(1)
        if (i.le.tpn2) then
          thr_limits(29,i+tpn-thrdat%subnrs(1)) = min(k,lastval(15)) + 1
          lastval(15) = nint(i*(1.0*k/(1.0*tpn2)))
          thr_limits(30,i+tpn-thrdat%subnrs(1)) = max(1,lastval(15))
        else
          thr_limits(29,i-thrdat%subnrs(1)) = 1
          thr_limits(30,i-thrdat%subnrs(1)) = 0
        end if
      else
        thr_limits(29,i) = min(k,lastval(15)) + 1
        lastval(15) = nint(i*(1.0*k/(1.0*tpn)))
        thr_limits(30,i) = max(1,lastval(15))
      end if
    else
      thr_limits(29,i) = 1
      thr_limits(30,i) = 0
    end if
!   31,32 are dynamic limits for interactions (global or just IPP/ATTLJ/WCA/SAV) when twin-range IS considered
    thr_limits(31,i) = min(nseq,lastval(16)) + 1
    lastval(16) = nint(i*(1.0*nseq/(1.0*tpn)))
    thr_limits(32,i) = max(1,lastval(16))
!   33,34 are dynamic limits for interactions (just POLAR/TABUL) when twin-range IS considered
    thr_limits(33,i) = min(nseq,lastval(17)) + 1
    lastval(17) = nint(i*(1.0*nseq/(1.0*tpn)))
    thr_limits(34,i) = max(1,lastval(17))
!   35,36 are dynamic limits for molecule-based nontrivial operations
    thr_limits(35,i) = min(nmol,lastval(18)) + 1
    lastval(18) = nint(i*(1.0*nmol/(1.0*tpn)))
    thr_limits(36,i) = max(1,lastval(18))
!   37-40 are dynamic limits for residue-based O(N) interactions (bonded and so on)
    call get_thread_loop_bounds(aone,afive,stx(1:2),stx(3:4),i)
    thr_limits(37,i) = stx(1)
    thr_limits(38,i) = stx(3)
    thr_limits(39,i) = min(n,lastval(20)) + 1
    lastval(20) = nint(i*(1.0*n/(1.0*tpn)))
    thr_limits(40,i) = max(1,lastval(20))
!   41,42 are dynamic limits for interactions (just POLAR/TABUL) when twin-range is not considered
    thr_limits(41,i) = min(nseq,lastval(21)) + 1
    lastval(21) = nint(i*(1.0*nseq/(1.0*tpn)))
    thr_limits(42,i) = max(1,lastval(21))
!   43,44 are dynamic limits for screened charge propagation
    call get_thread_loop_bounds(aone,afive,stx(1:2),stx(3:4),i)
    thr_limits(43,i) = stx(1)
    thr_limits(44,i) = stx(3)
!   45,46 are dynamic limits for residue operations scaling with the number of internal d.o.f.
    thr_limits(45,i) = min(nseq,lastval(23)) + 1
    lastval(23) = nint(i*(1.0*nseq/(1.0*tpn)))
    thr_limits(46,i) = max(1,lastval(23))
!   47-58 are various bounds for SC_DREST and SC_EMICRO
    if (ndrest.gt.0) then
      thr_limits(47,i) = min(ndrest,lastval(24)) + 1
      lastval(24) = nint(i*(1.0*ndrest/(1.0*tpn)))
      thr_limits(48,i) = max(1,lastval(24))
    else
      thr_limits(47,i) = 1
      thr_limits(48,i) = 0
    end if
    thr_limits(49,i) = min(emgrid%dim(3),lastval(25)) + 1
    lastval(25) = nint(i*(1.0*emgrid%dim(3)/(1.0*tpn)))
    thr_limits(50,i) = max(1,lastval(25))
    if ((use_EMICRO.EQV..true.).OR.(emcalc.le.nsim)) then
      k = emcoarsegrid%dim(2)
      if (allocated(emcoarsegrid%lblk).EQV..true.) k = emcoarsegrid%dimc(3)
      thr_limits(51,i) = min(k,lastval(26)) + 1
      lastval(26) = nint(i*(1.0*k/(1.0*tpn)))
      thr_limits(52,i) = max(1,lastval(26))
    else
      thr_limits(51,i) = 1
      thr_limits(52,i) = 0
    end if
    thr_limits(53,i) = min(n,lastval(27)) + 1
    lastval(27) = nint(i*(1.0*n/(1.0*tpn)))
    thr_limits(54,i) = max(1,lastval(27))
    if ((use_EMICRO.EQV..true.).OR.(emcalc.le.nsim)) then
      thr_limits(55,i) = min(emgrid%dim(3),lastval(28)) + 1
      lastval(28) = nint(i*(1.0*emgrid%dim(3)/(1.0*tpn)))
      thr_limits(56,i) = max(1,lastval(28))
    else
      thr_limits(55,i) = 1
      thr_limits(56,i) = 0
    end if
    if (use_EMICRO.EQV..true.) then
      thr_limits(57,i) = min(emcoarsegrid%dim(3),lastval(29)) + 1
      lastval(29) = nint(i*(1.0*emcoarsegrid%dim(3)/(1.0*tpn)))
      thr_limits(58,i) = max(1,lastval(29))
    else
      thr_limits(57,i) = 1
      thr_limits(58,i) = 0
    end if
    thr_limits(59,i) = min(nmol,lastval(30)) + 1 ! for fixed cost per mol
    lastval(30) = nint(i*(1.0*nmol/(1.0*tpn)))
    thr_limits(60,i) = max(1,lastval(30))
    thr_limits(61,i) = min(nmol,lastval(31)) + 1 ! for cost ~dof per mol with DLB-12
    lastval(31) = nint(i*(1.0*nmol/(1.0*tpn)))
    thr_limits(62,i) = max(1,lastval(31))
    if (unslst%nr.gt.0) then
      thr_limits(63,i) = min(unslst%nr,lastval(32)) + 1 ! fixed schedule in unslst%nr
      lastval(32) = nint(i*(1.0*unslst%nr/(1.0*tpn)))
      thr_limits(64,i) = max(1,lastval(32))
    else
      thr_limits(63,i) = 1
      thr_limits(64,i) = 0
    end if
    if (unklst%nr.gt.0) then
      thr_limits(65,i) = min(unklst%nr,lastval(33)) + 1 ! fixed schedule in unklst%nr
      lastval(33) = nint(i*(1.0*unklst%nr/(1.0*tpn)))
      thr_limits(66,i) = max(1,lastval(33))
    else
      thr_limits(65,i) = 1
      thr_limits(66,i) = 0
    end if
    do k=67,83,2
      thr_limits(k,i) = 1 ! reserved for temporary bounds
      thr_limits(k+1,i) = 0
    end do
!   85,86 are dynamic limits for the fxn Vforce_dsavdr
    thr_limits(85,i) = min(n,lastval(43)) + 1
    lastval(43) = nint(i*(1.0*n/(1.0*tpn)))
    thr_limits(86,i) = max(1,lastval(43))
!   87,88 are dynamic sampling weights for LREL_MD 4/5 forces
    thr_limits(87,i) = (i-1)*1000+1
    thr_limits(88,i) = i*1000
!   89,90 are static weights for complete N^2 interactions
#ifdef ORACLE_FORTRAN
    tmps(1) =  -min(iacnt,lval2) + 1 !thr_limits(89,i) = -min(iacnt,nint((i-1)*((1.0*iacnt)/(1.0*tpn)))) + 1
    lval2 = nint(int(i,KIND=8)*((1.0*iacnt)/(1.0*tpn)),KIND=8) ! nint(i*((1.0*iacnt)/(1.0*tpn))
    tmps(2) = -max(1_8,lval2) !thr_limits(90,i) = -max(1,nint(i_INT64*((1.0*iacnt)/(1.0*tpn)),KIND=8))
#else
    tmps(1) =  -min(iacnt,lval2) + 1 !thr_limits(89,i) = -min(iacnt,nint((i-1)*((1.0*iacnt)/(1.0*tpn)))) + 1
    lval2 = nint(int(i,KIND=INT64)*((1.0*iacnt)/(1.0*tpn)),KIND=INT64) ! nint(i*((1.0*iacnt)/(1.0*tpn))
    tmps(2) = -max(1_INT64,lval2) !thr_limits(90,i) = -max(1,nint(i_INT64*((1.0*iacnt)/(1.0*tpn)),KIND=8))
#endif
    kksum = 0
    
    do k=1,nseq
      kksum = kksum + at(k)%na*(n-(at(k)%bb(1)+at(k)%na-1)) + at(k)%na*(at(k)%na-1)/2
      if (tmps(1).le.1) then
        if (kksum.ge.-tmps(1)) then
          thr_limits(89,i) = min(nseq,k + 1)
          tmps(1) = -tmps(1) + 2
        end if
      end if
      if (tmps(2).lt.0) then
        if (kksum.gt.-tmps(2)) then
          thr_limits(90,i) = k
          tmps(2) = -tmps(2)
        end if
      end if
    end do
    if (i.eq.1) thr_limits(89,i) = 1
    if (tmps(2).lt.0) then
      thr_limits(90,i) = thr_limits(89,i)
    end if
    if (i.eq.tpn) thr_limits(90,i) = nseq
    if ((use_TABUL.EQV..true.).AND.(use_POLAR.EQV..false.)) then
      thr_limits(91,i) = 1
      thr_limits(92,i) = 0
    else
      thr_limits(91,i) = thr_limits(89,i)
      thr_limits(92,i) = thr_limits(90,i)
    end if
  end do
  thr_dlb(:,:) = -1
  thr_dlb(1:16,1) = 1
!
! use -1 to indicate permanent deactivation of load measuring or balancing for a specific block
  if (use_mcgrid.EQV..false.) thr_dlb(1,1) = -1
  if (use_cutoffs.EQV..false.) thr_dlb(2,1) = -1
  if (ideal_run.EQV..true.) then
    thr_dlb(3,1) = -1
    thr_dlb(13,1) = -1
    thr_dlb(16,1) = -1
  end if
  if ((no_shake.EQV..true.).OR.(fycxyz.eq.1)) thr_dlb(4,1) = -1
  if ((lrel_md.eq.1).OR.(lrel_md.eq.3)) thr_dlb(5,1) = -1
  if ((ideal_run.EQV..true.).OR.(is_fegplj.EQV..true.).OR.(is_fegprflj.EQV..true.).OR.(is_pewlj.EQV..true.).OR.&
 &    (is_plj.EQV..true.).OR.(is_prflj.EQV..true.).OR.(is_tab.EQV..true.).OR.(is_ev.EQV..true.).OR.&
 &    (is_lj.EQV..true.)) thr_dlb(7,1) = -1
  if ((use_IMPSOLV.EQV..false.).OR.(use_POLAR.EQV..false.)) thr_dlb(8,1) = -1
  if ((use_EMICRO.EQV..false.).AND.(emcalc.gt.nsim)) thr_dlb(11,1) = -1
  if (fycxyz.eq.2) thr_dlb(12,1) = -1
  if (use_cutoffs.EQV..false.) thr_dlb(13:14,1) = -1
  if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) thr_dlb(14,1) = -1
  if ((use_BOND(1).EQV..false.).AND.(use_BOND(2).EQV..false.).AND.(use_BOND(3).EQV..false.).AND.(use_BOND(4).EQV..false.).AND.&
 &    (use_BOND(5).EQV..false.).AND.(use_IMPSOLV.EQV..false.).AND.(.NOT.((bnd_type.eq.1).AND.(bnd_shape.eq.1)))) then
    if (bnd_type.eq.3) then
      thr_dlb(6,1) = -1
      thr_limits(37:38,1:tpn) = thr_limits(3:4,1:tpn)
    end if
  end if
  if (.NOT.((use_cutoffs.EQV..true.).AND.(use_dyn.EQV..true.).AND.((lrel_md.eq.4).OR.(lrel_md.eq.5)).AND.&
 &   (use_POLAR.EQV..true.))) then
    thr_dlb(16,1) = -1
  end if
  if (use_IMPSOLV.EQV..false.) thr_dlb(15,1) = -1
!
  thr_timings(:,:) = 0
  lastval(:) = 1
  do i=2,tpn
    do k=1,46
      if ((k.gt.33).AND.(k.lt.43)) cycle
      if (thr_limits(2*k,i).ge.thr_limits(2*k-1,i)) then
        if (thr_limits(2*k-1,i).le.thr_limits(2*k,lastval(k))) then
!          write(*,*) k,lastval(k),thr_limits(2*k-1,i),thr_limits(2*k,i),thr_limits(2*k-1,i),thr_limits(2*k,lastval(k))
          thr_limits(2*k-1,i) = thr_limits(2*k-1,i) + 1
        else
          lastval(k) = i
        end if
      end if
    end do
  end do
  do k=1,46
    if ((k.gt.33).AND.(k.lt.43)) cycle
    if (thr_limits(2*k,tpn).le.0) cycle
    allocate(testmap(thr_limits(2*k,tpn)))
    testmap(:) = 0
    do i=1,tpn
      if (thr_limits(2*k,i).ge.thr_limits(2*k-1,i)) then
        testmap(thr_limits(2*k-1,i):thr_limits(2*k,i)) = i
      end if
    end do
    if (minval(testmap).le.0) then
      write(ilog,'(1x,a,i3,a,i8,a)') 'Fatal. In setting up initial bounds for threads, category ',k,&
 &' with an apparent maximum of ',thr_limits(2*k,tpn),' contains one or more elements not accounted for. This is a bug.'
      call fexit()
    end if
    deallocate(testmap)
  end do

 567 format(a,' bounds ',i4,': ',10000(i7,1x))
  if (thrdat%verbosity.gt.2) then
    write(ithread,*) 'General initial bounds for individual threads:'
    do k=1,46
      if ((k.gt.33).AND.(k.lt.43)) cycle
      write(ithread,567) 'Lower',k,(thr_limits(2*k-1,i),i=1,tpn)
      write(ithread,567) 'Upper',k,(thr_limits(2*k,i),i=1,tpn)
    end do
    write(ithread,*) 
  end if
!
! the rest are pointers not bounds
  if (fycxyz.eq.2) then ! constraints and masses must be set in stone
    do tpi=1,tpn
      thr_limits(100,tpi) = 0
      do i=1,n
        if (mass(i).le.0.0) cycle
        do k=1,3
          if (cart_frz(i,k).EQV..true.) cycle
          if (i.ge.thr_limits(1,tpi)) exit
          thr_limits(100,tpi) = thr_limits(100,tpi) + 1
        end do
      end do
      thr_limits(101,tpi) = 0 ! the same thing without mass cycle
      do i=1,n
        do k=1,3
          if (cart_frz(i,k).EQV..true.) cycle
          if (i.ge.thr_limits(1,tpi)) exit
          thr_limits(101,tpi) = thr_limits(101,tpi) + 1
        end do
      end do
    end do
  end if
!  write(*,*) thr_dlb(1:11,1)
!
  if ((use_dyn.EQV..true.).AND.(use_mcgrid.EQV..true.).AND.(use_cutoffs.EQV..true.).AND.(use_POLAR.EQV..true.).AND.&
 &    (cglst%ncrs.gt.0)) then
    if (lrel_md.eq.5) then
      allocate(thr_lhlper(cglst%ncrs,nseq))
      thr_lhlper(:,:) = .false.
    else if (lrel_md.eq.4) then
      allocate(thr_lhlper(cglst%ncrs,cglst%ncrs))
      thr_lhlper(:,:) = .false.
    end if
  end if
  if ((use_mcgrid.EQV..true.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
    do i=1,thrdat%maxn
      allocate(thr_rsp_nbl(i)%tmpl(nseq))
    end do
  end if
!
end
!
!---------------------------------------------------------------------------------------------------------
!
subroutine threads_dlb(istep)
!
  use sequen
  use atoms
  use threads
  use mcgrid
  use iounit
  use shakeetal
  use cutoffs
  use molecule
  use ems
  use ewalds
  use system
  use energies, ONLY: ideal_run,is_tab
  use polypep, ONLY: at
!
  implicit none
!
  integer, INTENT(IN):: istep
!
  integer i,k,inc,totcs,shksz(7),allinc,tpn,tpi,ixes(6)
  RTYPE fracl(thrdat%maxn),libs(20),ttts(20),rinc,thrrat!,allms
  logical atrue
!
  atrue = .true.
!
!  tpn = thrdat%maxn
!  allms = 0.
!  if (thr_dlb(1,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(2,1:tpn)) - sum(thr_timings(1,1:tpn)))/(1.0*thr_dlb(1,2))
!  if (thr_dlb(2,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(4,1:tpn)) - sum(thr_timings(3,1:tpn)))/(1.0*thr_dlb(2,2))
!  if (thr_dlb(3,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(6,1:tpn)) - sum(thr_timings(5,1:tpn)))/(1.0*thr_dlb(3,2))
!  if (thr_dlb(4,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(10,1:tpn)) - sum(thr_timings(9,1:tpn)))/(1.0*thr_dlb(4,2))
!  if (thr_dlb(5,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(12,1:tpn)) - sum(thr_timings(11,1:tpn)))/(1.0*thr_dlb(5,2))
!  if (thr_dlb(6,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(14,1:tpn)) - sum(thr_timings(13,1:tpn)))/(1.0*thr_dlb(6,2))
!  if (thr_dlb(7,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(18,1:tpn)) - sum(thr_timings(17,1:tpn)))/(1.0*thr_dlb(7,2))
!  if (thr_dlb(8,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(20,1:tpn)) - sum(thr_timings(19,1:tpn)))/(1.0*thr_dlb(8,2))
!  if (thr_dlb(9,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(22,1:tpn)) - sum(thr_timings(21,1:tpn)))/(1.0*thr_dlb(9,2))
!  if (thr_dlb(10,2).gt.0) allms = allms + 1.0e3*(sum(thr_timings(24,1:tpn)) - sum(thr_timings(23,1:tpn)))/(1.0*thr_dlb(10,2))
!  write(ithread,*) 'Net cost per step of monitored loads is ~',allms/(1.0*tpn)/thrdat%rate,'ms.'
!
  thrrat = 0.05/(1.0*thrdat%maxn)
!
  if ((thr_dlb(1,1).gt.0).AND.(thr_dlb(1,2).eq.1*thrdat%dlbscale)) then ! DLB for grid scan
    allinc = 0
    thr_timings(101,1) = (sum(thr_timings(2,1:thrdat%maxn)) - sum(thr_timings(1,1:thrdat%maxn)))/thr_dlb(1,2)
    ttts(1) = 1.0e3*thr_timings(101,1)
    if (thr_timings(101,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(2,1:thrdat%maxn)-thr_timings(1,1:thrdat%maxn))/(1.0*thr_timings(101,1)*thr_dlb(1,2))
      libs(1) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*grid%pnts*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(1,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(6,i)-thr_limits(5,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(5,i-1)-thr_limits(6,i-1)-1,min(-1,int(max(-thrrat*thr_limits(6,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(6,i)-thr_limits(5,i)+1,max(1,int(min(thrrat*thr_limits(6,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
          if ((thr_limits(6,i).gt.thr_limits(5,i)).AND.((inc.gt.0).OR.(thr_limits(6,i-1).gt.thr_limits(5,i-1)))) then
            thr_limits(5,i) = thr_limits(5,i) + inc
            thr_limits(6,i-1) = thr_limits(6,i-1) + inc
          end if
        end do
      end if
      thr_timings(1:2,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(1,1) = 0
!        write(*,756) 'GBL',fracl(:)*thr_timings(101,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final bounds for grid-based neighbor scan are:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Grid points ',thr_limits(5,k),thr_limits(6,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for grid neighbor search is currently at ',libs(1),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',ttts(1)*0.01*libs(1)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(1,1) = -1
    end if
    thr_dlb(1,2) = 0
  end if
!
756 format(100(g10.3,1x))
  if ((thr_dlb(2,1).gt.0).AND.(thr_dlb(2,2).eq.1*thrdat%dlbscale)) then ! DLB for residues for NB generation
    allinc = 0
    thr_timings(102,1) = (sum(thr_timings(4,1:thrdat%maxn)) - sum(thr_timings(3,1:thrdat%maxn)))/thr_dlb(2,2)
    ttts(2) = 1.0e3*thr_timings(102,1)
    if (thr_timings(102,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(4,1:thrdat%maxn)-thr_timings(3,1:thrdat%maxn))/(1.0*thr_timings(102,1)*thr_dlb(2,2))
      libs(2) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*nseq*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(2,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(8,i)-thr_limits(7,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(7,i-1)-thr_limits(8,i-1)-1,min(-1,int(max(-thrrat*thr_limits(8,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(8,i)-thr_limits(7,i)+1,max(1,int(min(thrrat*thr_limits(8,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
          if ((thr_limits(8,i).gt.thr_limits(7,i)).AND.((inc.gt.0).OR.(thr_limits(8,i-1).gt.thr_limits(7,i-1)))) then
            thr_limits(7,i) = thr_limits(7,i) + inc
            thr_limits(8,i-1) = thr_limits(8,i-1) + inc
          end if
        end do
      end if
      thr_timings(3:4,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(2,1) = 0
!        write(*,756) 'NBL',fracl(:)*thr_timings(102,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final bounds for neighbor list construction are:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Residues ',thr_limits(7,k),thr_limits(8,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for neighbor list creation is currently at ',libs(2),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(2)*0.01*ttts(2)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(2,1) = -1
    end if
    thr_dlb(2,2) = 0
  end if
!
  if ((thr_dlb(3,1).gt.0).AND.(thr_dlb(3,2).eq.1*thrdat%dlbscale)) then ! DLB for residues for NBF
    allinc = 0
    tpn = thrdat%maxn - thrdat%subnrs(1)
    thr_timings(103,1) = (sum(thr_timings(6,1:tpn)) - sum(thr_timings(5,1:tpn)))/thr_dlb(3,2)
    ttts(3) = 1.0e3*thr_timings(103,1)
    if (thr_timings(103,1).gt.0) then
      fracl(1:tpn) = 1.0*(thr_timings(6,1:tpn)-thr_timings(5,1:tpn))/(1.0*thr_timings(103,1)*thr_dlb(3,2))
!      write(*,*) fracl(1:tpn),1.0e3*(thr_timings(6,1:tpn)-thr_timings(5,1:tpn))/thr_dlb(3,2)/thrdat%rate
      libs(3) = 100.0*(maxval(fracl)-(1.0/tpn))
      if ((is_tab.EQV..true.).OR.(use_cutoffs.EQV..false.).OR.(ideal_run.EQV..true.)) then
        if (int(0.5*nseq*(maxval(fracl)-minval(fracl))).le.0) then
          thr_dlb(3,1) = 0
        else
          do i=2,tpn
            rinc = 0.5*(thr_limits(10,i)-thr_limits(9,i)+1)*(fracl(i)-fracl(i-1))
            if (rinc.le.0.0) then
              inc = max(thr_limits(9,i-1)-thr_limits(10,i-1)-1,min(-1,int(max(-thrrat*thr_limits(10,thrdat%maxn),rinc))))
            else
              inc = min(thr_limits(10,i)-thr_limits(9,i)+1,max(1,int(min(thrrat*thr_limits(10,thrdat%maxn),rinc))))
            end if
            allinc = allinc + abs(inc)
            if ((thr_limits(10,i).gt.thr_limits(9,i)).AND.((inc.gt.0).OR.(thr_limits(10,i-1).gt.thr_limits(9,i-1)))) then
              thr_limits(9,i) = thr_limits(9,i) + inc
              thr_limits(10,i-1) = thr_limits(10,i-1) + inc
            end if
          end do
        end if
        thr_timings(5:6,1:tpn) = 0
        if (allinc.le.0) then
          thr_dlb(3,1) = 0
!          write(*,756) 'NBF',fracl(1:tpn)*thr_timings(103,1)*1.0e-3
          if (thrdat%verbosity.gt.1) then
            write(ithread,*) 'Final dynamic bounds for nonbonded force are:'
            do k=1,tpn
              write(ithread,751) k,'Residues ',thr_limits(9,k),thr_limits(10,k)
            end do
            write(ithread,*) 
          end if
        else if (thrdat%verbosity.gt.3) then
          write(ithread,752) 'Load imbalance for nonbonded force/energy is currently at ',libs(3),'%.'
          write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(3)*0.01*ttts(3)/thrdat%rate,'ms.'
        end if
      else
        ixes(1) = 9
        ixes(2) = 10
        ixes(3) = 69
        ixes(4) = 70
        ixes(5) = 71
        ixes(6) = 72
        call shift_thread_loop_bounds2(thrrat,fracl,atrue,ixes)
        thr_timings(5:6,1:thrdat%maxn) = 0
        if (thrdat%verbosity.gt.3) then
          write(ithread,752) 'Load imbalance for NBF without TR-update is currently at ',libs(3),'%.'
          write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(3)*0.01*ttts(3)/thrdat%rate,'ms.'
        end if
      end if
    else
      thr_dlb(3,1) = -1
    end if
    thr_dlb(3,2) = 0
  end if
!
  if ((thr_dlb(4,1).gt.0).AND.(thr_dlb(4,2).eq.1*thrdat%dlbscale)) then ! DLB for holonomic constraints
    shksz(1) = settle_tip3ps
    shksz(2) = settle_spcs
    shksz(3) = settle_tip4ps
    shksz(4) = settle_tip5ps
    shksz(5) = settle_tip4pes
    shksz(6) = settle_rest
    shksz(7) = shake_cnt
!
    totcs = 3*sum(shksz(1:6))
!
    thr_timings(104,1) = (sum(thr_timings(10,1:thrdat%maxn)) - sum(thr_timings(9,1:thrdat%maxn)))/thr_dlb(4,2)
    if ((totcs.le.0).AND.(shake_cnt.le.(nccgs+1))) then ! there is nothing to balance
      thr_dlb(4,1) = 0
      thr_timings(104,1) = 0.0
    end if
    if (shake_cnt.gt.0) then
      totcs = totcs + sum(constraints(1:shake_cnt)%nr+constraints(1:shake_cnt)%nr3+constraints(1:shake_cnt)%nr4)
    end if
!
    ttts(4) = 1.0e3*thr_timings(104,1)
    allinc = 0
    if (thr_timings(104,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(10,1:thrdat%maxn)-thr_timings(9,1:thrdat%maxn))/(1.0*thr_timings(104,1)*thr_dlb(4,2))
      libs(4) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(4,1) = 0
      else
        do k=1,7
          do i=2,thrdat%maxn
            if (shksz(k).gt.0) then
              rinc = 0.5*(thr_limits(10+2*k,i)-thr_limits(9+2*k,i)+1)*(fracl(i)-fracl(i-1))
              if (rinc.le.0.0) then
                inc = max(thr_limits(9+2*k,i-1)-thr_limits(10+2*k,i-1)-1,&
 &                        min(-1,int(max(-thrrat*thr_limits(10+2*k,thrdat%maxn),rinc))))
              else
                inc = min(thr_limits(10+2*k,i)-thr_limits(9+2*k,i)+1,max(1,int(min(thrrat*thr_limits(10+2*k,thrdat%maxn),rinc))))
              end if
              allinc = allinc + abs(inc)
!             both legitimate ranges of size at least 2
              if ((thr_limits(10+2*k,i).gt.thr_limits(9+2*k,i)).AND.&
 &                ((inc.gt.0).OR.(thr_limits(10+2*k,i-1).gt.thr_limits(9+2*k,i-1)))) then
                thr_limits(9+2*k,i) = thr_limits(9+2*k,i) + inc
                thr_limits(10+2*k,i-1) = thr_limits(10+2*k,i-1) + inc
              end if
            end if
          end do
        end do
      end if
      thr_timings(9:10,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(4,1) = 0
!        write(*,756) 'HOL',fracl(:)*thr_timings(104,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for holonomic constraint groups are:'
          do k=1,7
            if (shksz(k).le.0) cycle
            if (k.eq.1) write(ithread,*) 'SETTLE for Tip3P Waters'
            if (k.eq.2) write(ithread,*) 'SETTLE for SPC Waters'
            if (k.eq.3) write(ithread,*) 'SETTLE for Tip4P Waters'
            if (k.eq.4) write(ithread,*) 'SETTLE for Tip5P Waters'
            if (k.eq.5) write(ithread,*) 'SETTLE for Tip4P-Ew Waters'
            if (k.eq.6) write(ithread,*) 'SETTLE for other groups'
            if (k.eq.7) write(ithread,*) 'Constraint groups with iterative or approximate solver'
            do i=1,thrdat%maxn
              write(ithread,751) i,'Constraint groups ',thr_limits(9+2*k,i),thr_limits(10+2*k,i)
            end do
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for holonomic constraints is currently at ',libs(4),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(4)*0.01*ttts(4)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(4,1) = -1
    end if
    thr_dlb(4,2) = 0
  end if

  if ((thr_dlb(5,1).gt.0).AND.(thr_dlb(5,2).eq.1*thrdat%dlbscale)) then ! DLB for LR electrostatics
!
    if (thrdat%subnrs(1).le.0) then
      tpi = 1
      tpn = thrdat%maxn
    else
      tpi = thrdat%maxn - thrdat%subnrs(1) + 1
      tpn = thrdat%maxn
    end if
    allinc = 0
    totcs = 0
    if (lrel_md.eq.2) totcs = nseq
    if ((lrel_md.eq.4).OR.(lrel_md.eq.5))  totcs = cglst%ncrs
    thr_timings(105,1) = (sum(thr_timings(12,tpi:tpn)) - sum(thr_timings(11,tpi:tpn)))/thr_dlb(5,2)
    ttts(5) = 1.0e3*thr_timings(105,1)
    if (thr_timings(105,1).gt.0) then
      fracl(tpi:tpn) = 1.0*(thr_timings(12,tpi:tpn)-thr_timings(11,tpi:tpn))/(1.0*thr_timings(105,1)*thr_dlb(5,2))
      libs(5) = 100.0*(maxval(fracl)-(1.0/(tpn-tpi+1)))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(5,1) = 0
      else
        do i=tpi+1,tpn
          rinc = 0.5*(thr_limits(28,i)-thr_limits(27,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(27,i-1)-thr_limits(28,i-1)-1,min(-1,int(max(-thrrat*thr_limits(28,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(28,i)-thr_limits(27,i)+1,max(1,int(min(thrrat*thr_limits(28,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(28,i).gt.thr_limits(27,i)).AND.((inc.gt.0).OR.(thr_limits(28,i-1).gt.thr_limits(27,i-1)))) then
            thr_limits(27,i) = thr_limits(27,i) + inc
            thr_limits(28,i-1) = thr_limits(28,i-1) + inc
          end if
        end do
      end if
      thr_timings(11:12,tpi:tpn) = 0
      if (allinc.le.0) then
        thr_dlb(5,1) = 0
!         write(*,756) 'LRE',fracl(tpi:tpn)*thr_timings(105,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for long-range electrostatics are:'
          do k=tpi,tpn
            if (lrel_md.eq.4) then
              write(ithread,751) k+thrdat%subnrs(1),'Monopoles ',thr_limits(27,k),thr_limits(28,k)
            else
              write(ithread,751) k+thrdat%subnrs(1),'Residues ',thr_limits(27,k),thr_limits(28,k)
            end if
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        if (lrel_md.eq.2) then
          write(ithread,752) 'Load imbalance for Ewald/PME atomic operations is currently at ',libs(5),'%.'
          write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(5)*0.01*ttts(5)/thrdat%rate,'ms.'
        else
          write(ithread,752) 'Load imbalance for LREL_MD 4/5 neighbor list appending is currently at ',libs(5),'%.'
          write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(5)*0.01*ttts(5)/thrdat%rate,'ms.'
        end if
      end if
    else
     thr_dlb(5,1) = -1
    end if
    thr_dlb(5,2) = 0
  end if
!
  if ((thr_dlb(6,1).gt.0).AND.(thr_dlb(6,2).eq.1*thrdat%dlbscale)) then ! DLB for other interactions parsed by residue
!
    allinc = 0
    totcs = n
    thr_timings(106,1) = (sum(thr_timings(14,1:thrdat%maxn)) - sum(thr_timings(13,1:thrdat%maxn)))/thr_dlb(6,2)
    ttts(6) = 1.0e3*thr_timings(106,1)
    if (thr_timings(106,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(14,1:thrdat%maxn)-thr_timings(13,1:thrdat%maxn))/(1.0*thr_timings(106,1)*thr_dlb(6,2))
!     write(*,*) fracl(1:thrdat%maxn),1.0e3*(thr_timings(14,1:thrdat%maxn)-thr_timings(13,1:thrdat%maxn))/thr_dlb(6,2)/thrdat%rate
      libs(6) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(6,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(40,i)-thr_limits(39,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(39,i-1)-thr_limits(40,i-1)-1,min(-1,int(max(-thrrat*thr_limits(40,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(40,i)-thr_limits(39,i)+1,max(1,int(min(thrrat*thr_limits(40,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(40,i).gt.thr_limits(39,i)).AND.((inc.gt.0).OR.(thr_limits(40,i-1).gt.thr_limits(39,i-1)))) then
            thr_limits(39,i) = thr_limits(39,i) + inc
            thr_limits(40,i-1) = thr_limits(40,i-1) + inc
            k = 0
            if ((inc.gt.0).AND.(thr_limits(38,i).ge.thr_limits(37,i))) then
              do while (at(thr_limits(37,i)+k)%bb(1).le.thr_limits(39,i))
                if ((thr_limits(37,i)+k).eq.thr_limits(38,i)) exit
                k = k + 1
              end do
            else if (thr_limits(37,i).gt.0) then
              do while (at(thr_limits(37,i)+k)%bb(1).gt.thr_limits(39,i))
                if ((thr_limits(37,i)+k).eq.(thr_limits(37,i-1)+1)) exit
                k = k - 1
              end do
            end if
            thr_limits(37,i) = thr_limits(37,i) + k
            thr_limits(38,i-1) = thr_limits(38,i-1) + k
          end if
        end do
      end if
      thr_timings(13:14,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(6,1) = 0
!         write(*,756) 'BON',fracl(:)*thr_timings(106,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for bonded and similar interactions are:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Residues ',thr_limits(37,k),thr_limits(38,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for other residue-based interactions is currently at ',libs(6),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(6)*0.01*ttts(6)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(6,1) = -1
    end if
    thr_dlb(6,2) = 0
  end if
!
  if ((thr_dlb(7,1).gt.0).AND.(thr_dlb(7,2).eq.1*thrdat%dlbscale)) then ! DLB for possible second round of NB forces without TR
!
    allinc = 0
    totcs = nseq
    thr_timings(107,1) = (sum(thr_timings(18,1:thrdat%maxn)) - sum(thr_timings(17,1:thrdat%maxn)))/thr_dlb(7,2)
    ttts(7) = 1.0e3*thr_timings(107,1)
    if (thr_timings(107,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(18,1:thrdat%maxn)-thr_timings(17,1:thrdat%maxn))/(1.0*thr_timings(107,1)*thr_dlb(7,2))
      libs(7) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if ((is_tab.EQV..true.).OR.(use_cutoffs.EQV..false.).OR.(ideal_run.EQV..true.)) then
        if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
          thr_dlb(7,1) = 0
        else
          do i=2,thrdat%maxn
            rinc = 0.5*(thr_limits(42,i)-thr_limits(41,i)+1)*(fracl(i)-fracl(i-1))
            if (rinc.le.0.0) then
              inc = max(thr_limits(41,i-1)-thr_limits(42,i-1)-1,min(-1,int(max(-thrrat*thr_limits(42,thrdat%maxn),rinc))))
            else
              inc = min(thr_limits(42,i)-thr_limits(41,i)+1,max(1,int(min(thrrat*thr_limits(42,thrdat%maxn),rinc))))
            end if
            allinc = allinc + abs(inc)
!            both legitimate ranges of size at least 2
            if ((thr_limits(42,i).gt.thr_limits(41,i)).AND.((inc.gt.0).OR.(thr_limits(42,i-1).gt.thr_limits(41,i-1)))) then
              thr_limits(41,i) = thr_limits(41,i) + inc
              thr_limits(42,i-1) = thr_limits(42,i-1) + inc
            end if
          end do
        end if
        thr_timings(17:18,1:thrdat%maxn) = 0
        if (allinc.le.0) then
          thr_dlb(7,1) = 0
!           write(*,756) 'LFO',fracl(:)*thr_timings(107,1)*1.0e-3
          if (thrdat%verbosity.gt.1) then
            write(ithread,*) 'Final dynamic bounds for long-range nonbonded interactions are:'
            do k=1,thrdat%maxn
              write(ithread,751) k,'Residues ',thr_limits(41,k),thr_limits(42,k)
            end do
            write(ithread,*) 
          end if
        else if (thrdat%verbosity.gt.3) then
          write(ithread,752) 'Load imbalance for long-range nonbonded interactions is currently at ',libs(7),'%.'
          write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(7)*0.01*ttts(7)/thrdat%rate,'ms.'
        end if
      else
        ixes(1) = 41
        ixes(2) = 42
        ixes(3) = 81
        ixes(4) = 82
        ixes(5) = 83
        ixes(6) = 84
        call shift_thread_loop_bounds2(thrrat,fracl,atrue,ixes)
        thr_timings(17:18,1:thrdat%maxn) = 0
        if (thrdat%verbosity.gt.3) then
          write(ithread,752) 'Load imbalance for LR-NBF without TR update is currently at ',libs(3),'%.'
          write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(3)*0.01*ttts(3)/thrdat%rate,'ms.'
        end if
      end if
    else
      thr_dlb(7,1) = -1
    end if
    thr_dlb(7,2) = 0
  end if
!
  if ((thr_dlb(8,1).gt.0).AND.(thr_dlb(8,2).eq.1*thrdat%dlbscale)) then ! DLB for screened charge derivative propagation
!
    allinc = 0
    totcs = nseq
    thr_timings(108,1) = (sum(thr_timings(20,1:thrdat%maxn)) - sum(thr_timings(19,1:thrdat%maxn)))/thr_dlb(8,2)
    ttts(8) = 1.0e3*thr_timings(108,1)
    if (thr_timings(108,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(20,1:thrdat%maxn)-thr_timings(19,1:thrdat%maxn))/(1.0*thr_timings(108,1)*thr_dlb(8,2))
      libs(8) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(8,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(44,i)-thr_limits(43,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(43,i-1)-thr_limits(44,i-1)-1,min(-1,int(max(-thrrat*thr_limits(44,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(44,i)-thr_limits(43,i)+1,max(1,int(min(thrrat*thr_limits(44,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(44,i).gt.thr_limits(43,i)).AND.((inc.gt.0).OR.(thr_limits(44,i-1).gt.thr_limits(43,i-1)))) then
            thr_limits(43,i) = thr_limits(43,i) + inc
            thr_limits(44,i-1) = thr_limits(44,i-1) + inc
          end if
        end do
      end if
      thr_timings(19:20,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(8,1) = 0
!         write(*,756) 'SCQ',fracl(:)*thr_timings(108,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for propagation of screened charge derivatives are:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Residues ',thr_limits(43,k),thr_limits(44,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for propagation of screened charge derivatives is currently at ',libs(8),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(8)*0.01*ttts(8)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(8,1) = -1
    end if
    thr_dlb(8,2) = 0
  end if
!
  if ((thr_dlb(9,1).gt.0).AND.(thr_dlb(9,2).eq.1*thrdat%dlbscale)) then ! DLB for residue-based operations 
!                                                                         scaling with # of internal d.o.f.s
    allinc = 0
    totcs = nseq
    thr_timings(109,1) = (sum(thr_timings(22,1:thrdat%maxn)) - sum(thr_timings(21,1:thrdat%maxn)))/thr_dlb(9,2)
    ttts(9) = 1.0e3*thr_timings(109,1)
    if (thr_timings(109,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(22,1:thrdat%maxn)-thr_timings(21,1:thrdat%maxn))/(1.0*thr_timings(109,1)*thr_dlb(9,2))
      libs(9) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
!         write(*,756) 'prior ',thr_limits(45:46,1:thrdat%maxn),sum(thr_limits(46,1:thrdat%maxn)-thr_limits(45,1:thrdat%maxn))
!         write(*,756) 'DOF',fracl(:)*thr_timings(109,1)*1.0e-3
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(9,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(46,i)-thr_limits(45,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(45,i-1)-thr_limits(46,i-1)-1,min(-1,int(max(-thrrat*thr_limits(46,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(46,i)-thr_limits(45,i)+1,max(1,int(min(thrrat*thr_limits(46,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(46,i).gt.thr_limits(45,i)).AND.((inc.gt.0).OR.(thr_limits(46,i-1).gt.thr_limits(45,i-1)))) then
            thr_limits(45,i) = thr_limits(45,i) + inc
            thr_limits(46,i-1) = thr_limits(46,i-1) + inc
          end if
        end do
!        inc = 0
!        do i=1,thrdat%maxn
!          do k=thr_limits(45,i),thr_limits(46,i)
!            if ((k.le.0).OR.(k.gt.nseq)) inc = inc + 1
!          end do
!1        end do
!        if ((sum(thr_limits(46,1:thrdat%maxn)-thr_limits(45,1:thrdat%maxn)).ne.(nseq-thrdat%maxn)).OR.(inc.gt.0)) then
!         write(*,756) 'after ',thr_limits(45:46,1:thrdat%maxn),sum(thr_limits(46,1:thrdat%maxn)-thr_limits(45,1:thrdat%maxn)),inc
!        end if
      end if
      thr_timings(21:22,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(9,1) = 0
!         write(*,756) 'DOF',fracl(:)*thr_timings(109,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for operations scaling with the number of native, internal degrees of freedom:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Residues ',thr_limits(45,k),thr_limits(46,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for operations scaling with the number of native, internal &
 &degrees of freedom is currently at ',libs(9),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(9)*0.01*ttts(9)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(9,1) = -1
    end if
    thr_dlb(9,2) = 0
  end if
!
  if ((thr_dlb(10,1).gt.0).AND.(thr_dlb(10,2).eq.1*thrdat%dlbscale)) then ! DLB for molecules
    allinc = 0
    totcs = nmol
    thr_timings(110,1) = (sum(thr_timings(24,1:thrdat%maxn)) - sum(thr_timings(23,1:thrdat%maxn)))/thr_dlb(10,2)
    ttts(10) = 1.0e3*thr_timings(110,1)
    if (thr_timings(110,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(24,1:thrdat%maxn)-thr_timings(23,1:thrdat%maxn))/&
 &                                (1.0*thr_timings(110,1)*thr_dlb(10,2))
      libs(10) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(10,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(36,i)-thr_limits(35,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(35,i-1)-thr_limits(36,i-1)-1,min(-1,int(max(-thrrat*thr_limits(36,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(36,i)-thr_limits(35,i)+1,max(1,int(min(thrrat*thr_limits(36,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(36,i).gt.thr_limits(35,i)).AND.((inc.gt.0).OR.(thr_limits(36,i-1).gt.thr_limits(35,i-1)))) then
            thr_limits(35,i) = thr_limits(35,i) + inc
            thr_limits(36,i-1) = thr_limits(36,i-1) + inc
          end if
        end do
      end if
      thr_timings(23:24,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(10,1) = -1  ! never re-enable
!         write(*,756) 'MOL',fracl(:)*thr_timings(110,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for molecule operations scaling with the number of atoms:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Molecules ',thr_limits(35,k),thr_limits(36,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for molecule operations scaling with the number of atoms is currently at ',libs(10),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(10)*0.01*ttts(10)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(10,1) = -1
    end if
    thr_dlb(10,2) = 0
  end if
!
  if ((thr_dlb(11,1).gt.0).AND.(thr_dlb(11,2).eq.1*thrdat%dlbscale)) then ! DLB for spatial density #1
    allinc = 0
    totcs = emgrid%dim(3)
    thr_timings(111,1) = (sum(thr_timings(26,1:thrdat%maxn)) - sum(thr_timings(25,1:thrdat%maxn)))/thr_dlb(11,2)
    ttts(11) = 1.0e3*thr_timings(111,1)
    if (thr_timings(111,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(26,1:thrdat%maxn)-thr_timings(25,1:thrdat%maxn))/&
 &                                (1.0*thr_timings(111,1)*thr_dlb(11,2))
      libs(11) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(11,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(56,i)-thr_limits(55,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(55,i-1)-thr_limits(56,i-1)-1,min(-1,int(max(-thrrat*thr_limits(56,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(56,i)-thr_limits(55,i)+1,max(1,int(min(thrrat*thr_limits(56,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(56,i).gt.thr_limits(55,i)).AND.((inc.gt.0).OR.(thr_limits(56,i-1).gt.thr_limits(55,i-1)))) then
            thr_limits(55,i) = thr_limits(55,i) + inc
            thr_limits(56,i-1) = thr_limits(56,i-1) + inc
          end if
        end do
      end if
      thr_timings(25:26,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(11,1) = 0
!         write(*,756) 'EMR',fracl(:)*thr_timings(111,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for z-dimension of spatial density analysis lattice:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Slices ',thr_limits(55,k),thr_limits(56,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for z-dimension of spatial density analysis lattice is currently at ',libs(11),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(11)*0.01*ttts(11)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(11,1) = -1
    end if
    thr_dlb(11,2) = 0
  end if
!
  if ((thr_dlb(12,1).gt.0).AND.(thr_dlb(12,2).eq.1*thrdat%dlbscale)) then ! DLB for molecules
    allinc = 0
    totcs = nmol
    thr_timings(112,1) = (sum(thr_timings(28,1:thrdat%maxn)) - sum(thr_timings(27,1:thrdat%maxn)))/thr_dlb(12,2)
    ttts(12) = 1.0e3*thr_timings(112,1)
    if (thr_timings(112,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(28,1:thrdat%maxn)-thr_timings(27,1:thrdat%maxn))/&
 &                               (1.0*thr_timings(112,1)*thr_dlb(12,2))
      libs(12) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(12,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(62,i)-thr_limits(61,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(61,i-1)-thr_limits(62,i-1)-1,min(-1,int(max(-thrrat*thr_limits(62,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(62,i)-thr_limits(61,i)+1,max(1,int(min(thrrat*thr_limits(62,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
!          both legitimate ranges of size at least 2
          if ((thr_limits(62,i).gt.thr_limits(61,i)).AND.((inc.gt.0).OR.(thr_limits(62,i-1).gt.thr_limits(61,i-1)))) then
            thr_limits(61,i) = thr_limits(61,i) + inc
            thr_limits(62,i-1) = thr_limits(62,i-1) + inc
          end if
        end do
      end if
      thr_timings(27:28,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(12,1) = -1 ! never re-enable
!         write(*,756) 'MDF',fracl(:)*thr_timings(112,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for molecule operations scaling with the number of dofs:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Molecules ',thr_limits(61,k),thr_limits(62,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for molecule operations scaling with the number of dofs is currently at ',libs(12),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(12)*0.01*ttts(12)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(12,1) = -1
    end if
    thr_dlb(12,2) = 0
  end if
!
  if ((thr_dlb(13,1).gt.0).AND.(thr_dlb(13,2).eq.1*thrdat%dlbscale)) then ! DLB for NBF with TRU
    allinc = 0
    thr_timings(113,1) = (sum(thr_timings(30,1:thrdat%maxn)) - sum(thr_timings(29,1:thrdat%maxn)))/thr_dlb(13,2)
    ttts(13) = 1.0e3*thr_timings(113,1)
    if (thr_timings(113,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(30,1:thrdat%maxn)-thr_timings(29,1:thrdat%maxn))/&
 &                               (1.0*thr_timings(113,1)*thr_dlb(13,2))
!      write(*,*) fracl(1:thrdat%maxn),1.0e3*(thr_timings(30,1:thrdat%maxn)-thr_timings(29,1:thrdat%maxn))/thr_dlb(13,2)/thrdat%rate
      libs(13) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      ixes(1) = 31
      ixes(2) = 32
      ixes(3) = 73
      ixes(4) = 74
      ixes(5) = 75
      ixes(6) = 76
      call shift_thread_loop_bounds2(thrrat,fracl,atrue,ixes)
      thr_timings(29:30,1:thrdat%maxn) = 0
      if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for NBF with TR update is currently at ',libs(13),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(13)*0.01*ttts(13)/thrdat%rate,'ms.'
      end if
!      thr_dlb(13,1) = -1
    end if
    thr_dlb(13,2) = 0
  end if
!
  if ((thr_dlb(14,1).gt.0).AND.(thr_dlb(14,2).eq.1*thrdat%dlbscale)) then ! DLB for LR-NBF only with TR
    allinc = 0
    thr_timings(114,1) = (sum(thr_timings(32,1:thrdat%maxn)) - sum(thr_timings(31,1:thrdat%maxn)))/thr_dlb(14,2)
    ttts(14) = 1.0e3*thr_timings(114,1)
    if (thr_timings(114,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(32,1:thrdat%maxn)-thr_timings(31,1:thrdat%maxn))/&
 &                               (1.0*thr_timings(114,1)*thr_dlb(14,2))
!      write(*,*) fracl(1:thrdat%maxn),1.0e3*(thr_timings(32,1:thrdat%maxn)-thr_timings(31,1:thrdat%maxn))/thr_dlb(14,2)/thrdat%rate
      libs(14) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      ixes(1) = 33
      ixes(2) = 34
      ixes(3) = 77
      ixes(4) = 78
      ixes(5) = 79
      ixes(6) = 80
      call shift_thread_loop_bounds2(thrrat,fracl,atrue,ixes)
      thr_timings(31:32,1:thrdat%maxn) = 0
      if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for nonbonded LR-only w/ TR force/energy is currently at ',libs(14),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(14)*0.01*ttts(14)/thrdat%rate,'ms.'
      end if
!      thr_dlb(14,1) = -1
    end if
    thr_dlb(14,2) = 0
  end if
!
  if ((thr_dlb(15,1).gt.0).AND.(thr_dlb(15,2).eq.1*thrdat%dlbscale)) then ! DLB for Vforce_dsavdr
!
    allinc = 0
    totcs = n
    thr_timings(115,1) = (sum(thr_timings(34,1:thrdat%maxn)) - sum(thr_timings(33,1:thrdat%maxn)))/thr_dlb(15,2)
    ttts(15) = 1.0e3*thr_timings(115,1)
    if (thr_timings(115,1).gt.0) then
      fracl(1:thrdat%maxn) =1.0*(thr_timings(34,1:thrdat%maxn)-thr_timings(33,1:thrdat%maxn))/(1.0*thr_timings(115,1)*thr_dlb(15,2))
!     write(*,*) fracl(1:thrdat%maxn),1.0e3*(thr_timings(34,1:thrdat%maxn)-thr_timings(33,1:thrdat%maxn))/thr_dlb(15,2)/thrdat%rate
      libs(15) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(15,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(86,i)-thr_limits(85,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(85,i-1)-thr_limits(86,i-1)-1,min(-1,int(max(-thrrat*thr_limits(86,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(86,i)-thr_limits(85,i)+1,max(1,int(min(thrrat*thr_limits(86,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
          if ((thr_limits(86,i).gt.thr_limits(85,i)).AND.((inc.gt.0).OR.(thr_limits(86,i-1).gt.thr_limits(85,i-1)))) then
            thr_limits(85,i) = thr_limits(85,i) + inc
            thr_limits(86,i-1) = thr_limits(86,i-1) + inc
          end if
        end do
      end if
      thr_timings(33:34,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(15,1) = 0
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for solvation derivative propagation are:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Atoms ',thr_limits(85,k),thr_limits(86,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for Vforce_dsavdr is currently at ',libs(15),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(15)*0.01*ttts(15)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(15,1) = -1
    end if
    thr_dlb(15,2) = 0
  end if
!
  if ((thr_dlb(16,1).gt.0).AND.(thr_dlb(16,2).eq.1*thrdat%dlbscale)) then ! DLB for LREL_MD 4/5
!
    allinc = 0
    totcs = thrdat%maxn*1000 
    thr_timings(116,1) = (sum(thr_timings(36,1:thrdat%maxn)) - sum(thr_timings(35,1:thrdat%maxn)))/thr_dlb(16,2)
    ttts(16) = 1.0e3*thr_timings(116,1)
    if (thr_timings(116,1).gt.0) then
      fracl(1:thrdat%maxn) = 1.0*(thr_timings(36,1:thrdat%maxn)-thr_timings(35,1:thrdat%maxn))/&
 &(1.0*thr_timings(116,1)*thr_dlb(16,2))
!     write(*,*) fracl(1:thrdat%maxn),1.0e3*(thr_timings(36,1:thrdat%maxn)-thr_timings(35,1:thrdat%maxn))/thr_dlb(16,2)/thrdat%rate
      libs(16) = 100.0*(maxval(fracl)-(1.0/thrdat%maxn))
      if (int(0.5*totcs*(maxval(fracl)-minval(fracl))).le.0) then
        thr_dlb(16,1) = 0
      else
        do i=2,thrdat%maxn
          rinc = 0.5*(thr_limits(88,i)-thr_limits(87,i)+1)*(fracl(i)-fracl(i-1))
          if (rinc.le.0.0) then
            inc = max(thr_limits(87,i-1)-thr_limits(88,i-1)-1,min(-1,int(max(-thrrat*thr_limits(88,thrdat%maxn),rinc))))
          else
            inc = min(thr_limits(88,i)-thr_limits(87,i)+1,max(1,int(min(thrrat*thr_limits(88,thrdat%maxn),rinc))))
          end if
          allinc = allinc + abs(inc)
          if ((thr_limits(88,i).gt.thr_limits(87,i)).AND.((inc.gt.0).OR.(thr_limits(88,i-1).gt.thr_limits(87,i-1)))) then
            thr_limits(87,i) = thr_limits(87,i) + inc
            thr_limits(88,i-1) = thr_limits(88,i-1) + inc
          end if
        end do
      end if
      thr_timings(35:36,1:thrdat%maxn) = 0
      if (allinc.le.0) then
        thr_dlb(16,1) = 0
!         write(*,756) 'BON',fracl(:)*thr_timings(116,1)*1.0e-3
        if (thrdat%verbosity.gt.1) then
          write(ithread,*) 'Final dynamic bounds for LREL_MD 4/5 force evaluations are:'
          do k=1,thrdat%maxn
            write(ithread,751) k,'Residues ',thr_limits(37,k),thr_limits(38,k)
          end do
          write(ithread,*) 
        end if
      else if (thrdat%verbosity.gt.3) then
        write(ithread,752) 'Load imbalance for long-range electrostatic force (LREL_MD 4/5) is currently at ',libs(16),'%.'
        write(ithread,752) 'Possible savings per time step with optimal balance would be ',libs(16)*0.01*ttts(16)/thrdat%rate,'ms.'
      end if
    else
      thr_dlb(16,1) = -1
    end if
    thr_dlb(16,2) = 0
  end if
!
!
 751 format(2x,'Thread ',i4,': ',a,i8,' to ',i8)
 753 format(2x,'Thread ',i4,': ',a,i8,' to ',i8,' and ',a,i8,' to ',i8)
 752 format(a,g10.3,a)
!
  do i=1,20
    if (thr_dlb(i,1).gt.0) then
      if (thr_dlb(i,2).eq.0) then
        thr_dlb(i,3) = thr_dlb(i,3) + 1*thrdat%dlbscale
        if (thr_dlb(i,3).ge.thrdat%dlblen) then ! disable
          thr_dlb(i,1:3) = 0
        end if
      end if
    else if (thr_dlb(i,1).eq.0) then
      thr_dlb(i,3) = thr_dlb(i,3) + 1
      if (thr_dlb(i,3).ge.thrdat%dlbfreq) then ! re-enable
        thr_dlb(i,2:3) = 0
        thr_dlb(i,1) = 1
      end if
    end if
  end do
!
  if ((lrel_md.eq.2).AND.(use_dyn.EQV..true.)) then
    i = mod(istep,thrdat%dlbfreq)
    if ((i.eq.0).AND.(thr_timings(105,1).gt.0).AND.(thr_timings(103,1).gt.0)) then
      if (thrdat%maxn.gt.1) then
        if (thrdat%subnrs(1).eq.0) then
          inc = 0 ! nint(thrdat%maxn*1.0*thr_timings(105,1)/(1.0*(thr_timings(103,1)+thr_timings(105,1)))) ! guess
        else
          inc = thrdat%subnrs(1) ! nint(thrdat%maxn*1.0*thr_timings(105,1)/(1.0*(thr_timings(103,1)+thr_timings(105,1))))
        end if
        if (inc.ne.thrdat%subnrs(1)) then
          totcs = thrdat%maxn - inc
          do i=1,thrdat%maxn
            if (use_waterloops.EQV..false.) then
              if (i.le.totcs) then  
                thr_limits(9,i) = min(nseq,nint((i-1)*(1.0*nseq/(1.0*totcs)))) + 1
                thr_limits(10,i) = max(1,nint(i*(1.0*nseq/(1.0*totcs))))
              else
                thr_limits(9,i) = 1
                thr_limits(10,i) = 0
              end if
            else
              if (i.le.totcs) then  
                k = rsw1 - 1
                thr_limits(31,i) = min(k,nint((i-1)*(1.0*k/(1.0*totcs)))) + 1
                thr_limits(32,i) = max(1,nint(i*(1.0*k/(1.0*totcs))))
                if (k.le.0) thr_limits(32,i) = 0
                k = nseq - rsw1 + 1
                thr_limits(33,i) = min(k,nint((i-1)*(1.0*k/(1.0*totcs)))) + rsw1
                thr_limits(34,i) = max(1,nint(i*(1.0*k/(1.0*totcs)))) + rsw1 - 1
              else
                thr_limits(31,i) = 1
                thr_limits(32,i) = 0
                thr_limits(33,i) = 1
                thr_limits(34,i) = 0
              end if
              if (i.le.inc) then
                thr_limits(27,i+totcs) = min(nseq,nint((i-1)*(1.0*nseq/(1.0*inc)))) + 1
                thr_limits(28,i+totcs) = max(1,nint(i*(1.0*nseq/(1.0*inc))))
                if (ewald_mode.eq.2) then
                  k = kdims(3,2)-kdims(3,1)+1 
                else
                  k = kdims(3,2)
                end if
                thr_limits(29,i+totcs) = min(k,nint((i-1)*(1.0*k/(1.0*inc)))) + 1
                thr_limits(30,i+totcs) = max(1,nint(i*(1.0*k/(1.0*inc))))
              else
                thr_limits(27,i-inc) = 1
                thr_limits(28,i-inc) = 0
                thr_limits(29,i-inc) = 1
                thr_limits(30,i-inc) = 0
              end if
            end if
          end do
          write(ithread,*) 'Switching from ',thrdat%subnrs(1),' to ',inc,' dedicated PME nodes.'
          if (ewald_mode.eq.1) call pme_plan_threads(inc)
          thrdat%subnrs(1) = inc
        else
!         do nothing
        end if
      end if
    end if
  end if     
!
end
!
!-------------------------------------------------------------------------------------------------------------
!
! a subroutine to deal with incremental energy calculation loop bounds in MC
! this is difficult/impossible to DLB as the "dimensions" of the
! loop change unpredictably at runtime; therefore, we use an explicit heuristic based on interaction numbers
! a fourth, unrelated mode is provided -> clusters and numbers of outgoing links
! a fifth mode as well -> residues weighted by numbers of atoms
!
subroutine get_thread_loop_bounds(mode,lrsr,stas,stos,tpi)
!
  use threads
  use sequen
  use cutoffs
  use clusters
  use polypep, ONLY: at
  use atoms, ONLY: n,atmres
!
  implicit none
!
  integer, INTENT(IN):: lrsr,mode,tpi
!
  integer(KIND=8) allia,allrs,staia,stoia,imyia,ipria
  integer stas(2),stos(2),shf,i,ii,j,laf,fif
  RTYPE myia
!
  if ((mode.ne.1).AND.(mode.ne.2)) then
    call fexit()
  end if
!
  shf = 0
  if (mode.eq.2) shf = thrdat%maxn
  stas(:) = 0
  stos(:) = -1
!
! short-range list
  if (lrsr.eq.1) then
    allia = 0
    allrs = 0
    fif = 0
    laf = 0
    do i=1+shf,thrdat%maxn+shf
      if (thr_rsp_nbl(i)%srnrs.le.0) cycle
      allia = allia + thr_rsp_nbl(i)%sr(thr_rsp_nbl(i)%srnrs,3)
      allrs = allrs + thr_rsp_nbl(i)%srnrs
!      if (tpi.le.1) write(*,*) i,thr_rsp_nbl(i)%srnrs,thr_rsp_nbl(i)%sr(thr_rsp_nbl(i)%srnrs,3)
    end do
    if (allrs.le.thrdat%maxn) then
      do i=1+shf,thrdat%maxn+shf
        do j=1,thr_rsp_nbl(i)%srnrs
          laf = laf + 1
          if (laf.eq.tpi) then
            stas(1) = i
            stos(1) = i
            stas(2) = j
            stos(2) = j
            exit
          end if
        end do
        if (laf.eq.tpi) exit
      end do
      return
    end if
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    allia = 0
    do i=1+shf,thrdat%maxn+shf
      if (thr_rsp_nbl(i)%srnrs.le.0) cycle
      laf = i
      if (fif.le.0) fif = i
      allia = allia + thr_rsp_nbl(i)%sr(thr_rsp_nbl(i)%srnrs,3)
      if ((allia.gt.nint((tpi-1)*myia,KIND=8)).AND.(stas(1).le.0)) then
        stas(1) = i
      end if
      if ((allia.ge.nint(tpi*myia,KIND=8)).AND.(stos(1).le.-1)) then
        stos(1) = i
      end if
    end do
    if (tpi.eq.1) stas(1) = fif
    if (tpi.eq.thrdat%maxn) stos(1) = laf
    allia = 0
    do i=1+shf,stos(1)
      if (thr_rsp_nbl(i)%srnrs.le.0) cycle
      if ((i.ge.stas(1)).AND.(stas(2).le.0)) then
        j = 1
        do while (j.le.thr_rsp_nbl(i)%srnrs)
          if (((allia+thr_rsp_nbl(i)%sr(j,3)).gt.nint((tpi-1)*myia,KIND=8)).AND.&
 &            ((allia+thr_rsp_nbl(i)%sr(j,3)).le.nint(tpi*myia,KIND=8))) then
            stas(2) = j
            exit
          end if
          j = j + 1
        end do
      end if
      if ((stas(2).gt.0).AND.(i.eq.stos(1))) then
        j = 1
        do while (j.le.thr_rsp_nbl(i)%srnrs)
          if ((allia+thr_rsp_nbl(i)%sr(j,3)).gt.nint(tpi*myia,KIND=8)) then
            stos(2) = j-1
            exit
          else if ((allia+thr_rsp_nbl(i)%sr(j,3)).eq.nint(tpi*myia,KIND=8)) then
            stos(2) = j
            exit
          end if
          j = j + 1
        end do
      end if
      allia = allia + thr_rsp_nbl(i)%sr(thr_rsp_nbl(i)%srnrs,3)
    end do
    if (tpi.eq.1) stas(2) = 1
    if (tpi.eq.thrdat%maxn) stos(2) = thr_rsp_nbl(laf)%srnrs
    if (stos(1).gt.stas(1)) then
      if ((stos(2).le.0).AND.(stas(2).ge.1)) then
        stos(1) = stos(1) - 1
        if (stos(1).gt.shf) stos(2) = thr_rsp_nbl(stos(1))%srnrs
      end if
    end if
! mid-range list
  else if (lrsr.eq.2) then
    allia = 0
    allrs = 0
    fif = 0
    laf = 0
    do i=1+shf,thrdat%maxn+shf
      if (thr_rsp_nbl(i)%trnrs.le.0) cycle
      allia = allia + thr_rsp_nbl(i)%tr(thr_rsp_nbl(i)%trnrs,3)
      allrs = allrs + thr_rsp_nbl(i)%trnrs
    end do
    if (allrs.le.thrdat%maxn) then
      do i=1+shf,thrdat%maxn+shf
        do j=1,thr_rsp_nbl(i)%trnrs
          laf = laf + 1
          if (laf.eq.tpi) then
            stas(1) = i
            stos(1) = i
            stas(2) = j
            stos(2) = j
            exit
          end if
        end do
        if (laf.eq.tpi) exit
      end do
      return
    end if
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    allia = 0
    do i=1+shf,thrdat%maxn+shf
      if (thr_rsp_nbl(i)%trnrs.le.0) cycle
      laf = i
      if (fif.le.0) fif = i
      allia = allia + thr_rsp_nbl(i)%tr(thr_rsp_nbl(i)%trnrs,3)
      if ((allia.gt.nint((tpi-1)*myia,KIND=8)).AND.(stas(1).le.0)) then
        stas(1) = i
      end if
      if ((allia.ge.nint(tpi*myia,KIND=8)).AND.(stos(1).le.-1)) then
        stos(1) = i
      end if
    end do
    if (tpi.eq.1) stas(1) = fif
    if (tpi.eq.thrdat%maxn) stos(1) = laf
    allia = 0
    do i=1+shf,stos(1)
      if (thr_rsp_nbl(i)%trnrs.le.0) cycle
      if ((i.ge.stas(1)).AND.(stas(2).le.0)) then
        j = 1
        do while (j.le.thr_rsp_nbl(i)%trnrs)
          if (((allia+thr_rsp_nbl(i)%tr(j,3)).gt.nint((tpi-1)*myia,KIND=8)).AND.&
 &            ((allia+thr_rsp_nbl(i)%tr(j,3)).le.nint(tpi*myia,KIND=8))) then
            stas(2) = j
            exit
          end if
          j = j + 1
        end do
      end if
      if ((stas(2).gt.0).AND.(i.eq.stos(1))) then
        j = 1
        do while (j.le.thr_rsp_nbl(i)%trnrs)
          if ((allia+thr_rsp_nbl(i)%tr(j,3)).gt.nint(tpi*myia,KIND=8)) then
            stos(2) = j-1
            exit
          else if ((allia+thr_rsp_nbl(i)%tr(j,3)).eq.nint(tpi*myia,KIND=8)) then
            stos(2) = j
            exit
          end if
          j = j + 1
        end do
      end if
      allia = allia + thr_rsp_nbl(i)%tr(thr_rsp_nbl(i)%trnrs,3)
    end do
    if (tpi.eq.1) stas(2) = 1
    if (tpi.eq.thrdat%maxn) stos(2) = thr_rsp_nbl(laf)%trnrs
    if (stos(1).gt.stas(1)) then
      if ((stos(2).le.0).AND.(stas(2).ge.1)) then
        stos(1) = stos(1) - 1
        if (stos(1).gt.shf) stos(2) = thr_rsp_nbl(stos(1))%trnrs
      end if
    end if
! long-range list
  else if (lrsr.eq.3) then
    allia = 0
    allrs = 0
    fif = 0
    laf = 0
    do i=1+shf,thrdat%maxn+shf
      if (thr_rsp_nbl(i)%lrnrs.le.0) cycle
      allia = allia + thr_rsp_nbl(i)%lr(thr_rsp_nbl(i)%lrnrs,3)
      allrs = allrs + thr_rsp_nbl(i)%lrnrs
!      if (tpi.le.1) write(*,*) i,thr_rsp_nbl(i)%lrnrs,thr_rsp_nbl(i)%lr(thr_rsp_nbl(i)%lrnrs,3)
    end do
    if (allrs.le.thrdat%maxn) then
      do i=1+shf,thrdat%maxn+shf
        do j=1,thr_rsp_nbl(i)%lrnrs
          laf = laf + 1
          if (laf.eq.tpi) then
            stas(1) = i
            stos(1) = i
            stas(2) = j
            stos(2) = j
            exit
          end if
        end do
        if (laf.eq.tpi) exit
      end do
      return
    end if
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    allia = 0
    do i=1+shf,thrdat%maxn+shf
      if (thr_rsp_nbl(i)%lrnrs.le.0) cycle
      laf = i
      if (fif.le.0) fif = i
      allia = allia + thr_rsp_nbl(i)%lr(thr_rsp_nbl(i)%lrnrs,3)
      if ((allia.gt.nint((tpi-1)*myia,KIND=8)).AND.(stas(1).le.0)) then
        stas(1) = i
      end if
      if ((allia.ge.nint(tpi*myia,KIND=8)).AND.(stos(1).le.-1)) then
        stos(1) = i
      end if
    end do
    if (tpi.eq.1) stas(1) = fif
    if (tpi.eq.thrdat%maxn) stos(1) = laf
    allia = 0
    do i=1+shf,stos(1)
      if (thr_rsp_nbl(i)%lrnrs.le.0) cycle
      if ((i.ge.stas(1)).AND.(stas(2).le.0)) then
        j = 1
        do while (j.le.thr_rsp_nbl(i)%lrnrs)
          if (((allia+thr_rsp_nbl(i)%lr(j,3)).gt.nint((tpi-1)*myia,KIND=8)).AND.&
 &            ((allia+thr_rsp_nbl(i)%lr(j,3)).le.nint(tpi*myia,KIND=8))) then
            stas(2) = j
            exit
          end if
          j = j + 1
        end do
      end if
      if ((stas(2).gt.0).AND.(i.eq.stos(1))) then
        j = 1
        do while (j.le.thr_rsp_nbl(i)%lrnrs)
          if ((allia+thr_rsp_nbl(i)%lr(j,3)).gt.nint(tpi*myia,KIND=8)) then
            stos(2) = j-1
            exit
          else if ((allia+thr_rsp_nbl(i)%lr(j,3)).eq.nint(tpi*myia,KIND=8)) then
            stos(2) = j
            exit
          end if
          j = j + 1
        end do
      end if
      allia = allia + thr_rsp_nbl(i)%lr(thr_rsp_nbl(i)%lrnrs,3)
    end do
    if (tpi.eq.1) stas(2) = 1
    if (tpi.eq.thrdat%maxn) stos(2) = thr_rsp_nbl(laf)%lrnrs
    if (stos(1).gt.stas(1)) then
      if ((stos(2).le.0).AND.(stas(2).ge.1)) then
        stos(1) = stos(1) - 1
        if (stos(1).gt.shf) stos(2) = thr_rsp_nbl(stos(1))%lrnrs
      end if
    end if
! cluster list
  else if (lrsr.eq.4) then
    allia = 0
    allrs = 0
    fif = 0
    laf = 0
    if (mode.eq.1) then
      do i=1,nstruccls
        allia = allia + max(1,scluster(i)%nb) ! sum(scluster(1:nstruccls)%nb) -> cannot respect KIND=8 for allia as written
      end do
    else if (mode.eq.2) then
      do i=1,nstruccls
        allia = allia + scluster(i)%nmbrs ! sum(scluster(1:nstruccls)%nmbrs)
      end do
    end if
    if (nstruccls.le.thrdat%maxn) then
      do i=1,nstruccls
        laf = laf + 1
        if (laf.eq.tpi) then
          stas(1) = i
          stos(1) = i
          exit
        end if
        if (laf.eq.tpi) exit
      end do
      return
    end if
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    allia = 0
    do i=1,nstruccls
      if (mode.eq.1) then
        laf = i
        if (fif.le.0) fif = i
        allia = allia + max(1,scluster(i)%nb)
      else if (mode.eq.2) then
        if (scluster(i)%nmbrs.le.0) cycle
        laf = i
        if (fif.le.0) fif = i
        allia = allia + scluster(i)%nmbrs
      end if
      if ((allia.gt.nint((tpi-1)*myia,KIND=8)).AND.(stas(1).le.0)) then
        stas(1) = i
      end if
      if ((allia.ge.nint(tpi*myia,KIND=8)).AND.(stos(1).le.-1)) then
        stos(1) = i
        if (allia.gt.nint(tpi*myia,KIND=8)) stos(1) = stos(1) - 1
      end if
    end do
    if (tpi.eq.1) stas(1) = fif
    if (tpi.eq.thrdat%maxn) stos(1) = laf
! residues by atoms
  else if (lrsr.eq.5) then
    allia = n
    if (nseq.le.thrdat%maxn) then
      do i=1,nseq
        if (i.eq.tpi) then
          stas(1) = i
          stos(1) = i
          exit
        end if
        if (i.eq.tpi) exit
      end do
      return
    end if
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    allia = 0
    do i=1,nseq
      allia = allia + at(i)%na
      if ((allia.gt.nint((tpi-1)*myia,KIND=8)).AND.(stas(1).le.0)) then
        stas(1) = i
      end if
      if ((allia.ge.nint(tpi*myia,KIND=8)).AND.(stos(1).le.-1)) then
        stos(1) = i
        if (allia.gt.nint(tpi*myia,KIND=8)) stos(1) = stos(1) - 1
      end if
    end do
    if (tpi.eq.1) stas(1) = 1
    if (tpi.eq.thrdat%maxn) stos(1) = nseq
! long-range list for dynamics (rs_nbl)
  else if (lrsr.eq.6) then
    fif = 0
    if (lrel_md.eq.5) then
      allrs = 0 ! sum(rs_nbl(1:nseq)%nnblrs)
      allia = 0
      do i=1,nseq
        allrs = allrs + rs_nbl(i)%nnblrs
        allia = allia + at(i)%npol*rs_nbl(i)%nnblrats ! sum(at(1:nseq)%npol*rs_nbl(1:nseq)%nnblrats)
      end do
      if (allrs.le.thrdat%maxn) then
        laf = 0
        do i=1,nseq
          do j=1,rs_nbl(i)%nnblrs
            laf = laf + 1
            if (laf.eq.tpi) then
              stas(1) = i
              stos(1) = i
              stas(2) = j
              stos(2) = j
              exit
            end if
          end do
          if (laf.eq.tpi) exit
        end do
        return
      end if
      imyia = nint(thr_limits(88,tpi)*allia/(thrdat%maxn*1000.),KIND=8)
      if (tpi.gt.1) then
        ipria = nint(thr_limits(88,tpi-1)*allia/(thrdat%maxn*1000.),KIND=8)
      else
        ipria = 0
      end if
      allia = 0
      do i=1,nseq
        if (rs_nbl(i)%nnblrs.le.0) cycle
        allia = allia + at(i)%npol*rs_nbl(i)%nnblrats
        if ((allia.gt.ipria).AND.(stas(1).le.0)) then
          staia = allia - at(i)%npol*rs_nbl(i)%nnblrats
          stas(1) = i
        end if
        if ((allia.ge.imyia).AND.(stos(1).le.-1)) then
          stoia = allia - at(i)%npol*rs_nbl(i)%nnblrats
          stos(1) = i
        end if
      end do
      if (tpi.eq.1) stas(1) = 1
      if (tpi.eq.thrdat%maxn) stos(1) = nseq
      if (stas(1).gt.0) then
        allia = staia
        do j=1,rs_nbl(stas(1))%nnblrs
          allia = allia + at(stas(1))%npol*at(rs_nbl(stas(1))%nblr(j))%npol
          if ((allia.gt.ipria).AND.(allia.le.imyia)) then
            stas(2) = j
            exit
          end if
        end do
      end if
      if (tpi.eq.1) stas(2) = 1
      if ((stas(2).gt.0).AND.(stos(1).gt.0)) then
        allia = stoia
        do j=1,rs_nbl(stos(1))%nnblrs
          allia = allia + at(stos(1))%npol*at(rs_nbl(stos(1))%nblr(j))%npol
          if (allia.gt.imyia) then
            stos(2) = j-1
            exit
          else if (allia.eq.imyia) then
            stos(2) = j
            exit
          end if
        end do
      end if
      if (tpi.eq.1) stas(2) = 1
      if (tpi.eq.thrdat%maxn) stos(2) = rs_nbl(nseq)%nnblrs
    else if (lrel_md.eq.4) then
      allrs = 0
      do ii=1,cglst%ncs
        allrs = allrs + rs_nbl(atmres(cglst%it(ii)))%nnblrs
      end do
      allia = 0 ! sum(at(1:nseq)%nncgrps*rs_nbl(1:nseq)%nnblrats)
      do i=1,nseq
        allia = allia + at(i)%nncgrps*rs_nbl(i)%nnblrats
      end do 
      if (allrs.le.thrdat%maxn) then
        laf = 0
        do ii=1,cglst%ncs
          i = atmres(cglst%it(ii))
          do j=1,rs_nbl(i)%nnblrs
            laf = laf + 1
            if (laf.eq.tpi) then
              stas(1) = ii
              stos(1) = ii
              stas(2) = j
              stos(2) = j
              exit
            end if
          end do
          if (laf.eq.tpi) exit
        end do
        return
      end if
      imyia = nint(thr_limits(88,tpi)*allia/(thrdat%maxn*1000.),KIND=8)
      if (tpi.gt.1) then
        ipria = nint(thr_limits(88,tpi-1)*allia/(thrdat%maxn*1000.),KIND=8)
      else
        ipria = 0
      end if
      allia = 0
      do ii=1,cglst%ncs
        i = atmres(cglst%it(ii))
        if (rs_nbl(i)%nnblrs.le.0) cycle
        allia = allia + rs_nbl(i)%nnblrats
        if ((allia.gt.ipria).AND.(stas(1).le.0)) then
          staia = allia - rs_nbl(i)%nnblrats
          stas(1) = ii
        end if
        if ((allia.ge.imyia).AND.(stos(1).le.-1)) then
          stoia = allia - rs_nbl(i)%nnblrats
          stos(1) = ii
        end if
      end do
      if (tpi.eq.1) stas(1) = 1
      if (tpi.eq.thrdat%maxn) stos(1) = cglst%ncs
      if (stas(1).gt.0) then
        allia = staia
        ii = atmres(cglst%it(stas(1)))
        do j=1,rs_nbl(ii)%nnblrs
          allia = allia + at(rs_nbl(ii)%nblr(j))%nncgrps
          if ((allia.gt.ipria).AND.(allia.le.imyia)) then
            stas(2) = j
            exit
          end if
        end do
      end if
      if (tpi.eq.1) stas(2) = 1
      if ((stas(2).gt.0).AND.(stos(1).gt.0)) then
        allia = stoia
        ii = atmres(cglst%it(stos(1)))
        do j=1,rs_nbl(ii)%nnblrs
          allia = allia + at(rs_nbl(ii)%nblr(j))%nncgrps
          if (allia.gt.imyia) then
            stos(2) = j-1
            exit
          else if (allia.eq.imyia) then
            stos(2) = j
            exit
          end if
        end do
      end if
      if ((stas(2).gt.stos(2)).AND.(stas(1).eq.stos(1))) then ! no interactions
        stos(1) = 0
        stos(2) = -1
        return
      end if
      if (tpi.eq.1) stas(2) = 1
      if (tpi.eq.thrdat%maxn) stos(2) = rs_nbl(atmres(cglst%it(cglst%ncs)))%nnblrs
    end if
  end if
!
!  if (tpi.le.1) write(*,*) mode,lrsr,allrs,allia
! 55 format(i3,a2,10(i6,1x),g10.4)
!
end
!
!
!-------------------------------------------------------------------------------------------------------------
!
! a subroutine to deal with deterministic loop bounds for global forces/energies in dynamics
! this requires three levels: i) residue (for rs_nbl), ii) which type of rs_nbl subset, iii) index into latter
!
subroutine get_thread_loop_bounds2(fnbl,ixes,tpi)
!
  use threads
  use cutoffs
  use sequen, ONLY: nseq
  use polypep, ONLY: at
  use inter, ONLY: nrsintra,nrsnb,nrexpolin,nrexpolnb
  use energies, ONLY: use_POLAR
  use iounit, ONLY: ilog
!
  implicit none
!
  integer, INTENT(IN):: tpi,ixes(6)
  logical, INTENT(IN):: fnbl
!
  integer(KIND=8) allia,tmpia,lowia,fif,laf,hiwia,stbu3
  integer i,j,k,fiia,laia,stbu(2),topos,ii
  RTYPE myia
!
  laia = 0
  fiia = 0
  thr_limits(ixes(1),tpi) = 0
  thr_limits(ixes(3),tpi) = 0
  thr_limits(ixes(5),tpi) = 0
  thr_limits(ixes(2),tpi) = -1
  thr_limits(ixes(4),tpi) = -1
  thr_limits(ixes(6),tpi) = -1
!
! multipliers to adjust for:
! 1) NB interactions (total vs. LR vs. SR)
! 2) SR-topology interactions (rs_vec)
! 3) water loops

  if (fnbl.EQV..true.) then
!
    allia = 0
    do i=1,nseq
      allia = allia + (rs_nbl(i)%nnbtrats+rs_nbl(i)%ngnbtrats+rs_nbl(i)%nwnbtrats+&
 &             rs_nbl(i)%nnbats+rs_nbl(i)%ngnbats+rs_nbl(i)%nwnbats)*at(i)%na +  &
 &             nrsintra(i) + nrsnb(i)
    end do
!    allia = sum((rs_nbl(1:nseq)%nnbtrats+rs_nbl(1:nseq)%ngnbtrats+rs_nbl(1:nseq)%nwnbtrats+&
! &               rs_nbl(1:nseq)%nnbats+rs_nbl(1:nseq)%ngnbats+rs_nbl(1:nseq)%nwnbats)*at(1:nseq)%na) +  &
! &          sum(nrsintra(1:nseq)) + sum(nrsnb(1:nseq))
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
      allia = allia + sum(nrexpolin(1:nseq)) + sum(nrexpolnb(1:nseq))
    end if
!
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    lowia = nint((tpi-1)*myia,KIND=8)
    tmpia = 0
!   first loop: set residue bounds and loop index
    do i=1,nseq
      stbu(1) = thr_limits(ixes(1),tpi)
      stbu(2) = thr_limits(ixes(2),tpi)
      stbu3 = tmpia
      topos = nrsintra(i) + nrsnb(i)
      if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) topos = topos + nrexpolin(i) + nrexpolnb(i)
!     determine which residue and which list we fall into
      if (topos.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + topos
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 1
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 1
        end if
      end if
      if (rs_nbl(i)%nnbtrs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%nnbtrats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 2
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 2
        end if
      end if
      if (rs_nbl(i)%ngnbtrs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%ngnbtrats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 3
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 3
        end if
      end if
      if (rs_nbl(i)%nwnbtrs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%nwnbtrats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 4
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 4
        end if
      end if
      if (rs_nbl(i)%nnbs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%nnbats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 5
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 5
        end if
      end if
      if (rs_nbl(i)%ngnbs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%ngnbats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 6
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 6
        end if
      end if
      if (rs_nbl(i)%nwnbs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%nwnbats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 7
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 7
        end if
      end if
!     store interaction cumsum for first and last residues
      if (thr_limits(ixes(1),tpi).ne.stbu(1)) fif = stbu3
      if (thr_limits(ixes(2),tpi).ne.stbu(2)) laf = stbu3
    end do
!   if no interactions whatsoever, invalidate residue bounds and exit
    if (fiia.le.0) then
      thr_limits(ixes(1),tpi) = 1
      thr_limits(ixes(2),tpi) = 0
      return
    end if
    if (tpi.eq.1) then
      thr_limits(ixes(1),tpi) = 1 ! fiia
    end if
    if (tpi.eq.thrdat%maxn) then
      thr_limits(ixes(2),tpi) = nseq ! max(nseq-1,laia) ! laia
    end if
!
!   the starting bound is adjusted (+1) in the innermost bound, thus requiring clean-up
!   if outer bound(s) shift 
    j = thr_limits(ixes(1),tpi)
    topos = nrsintra(j) + nrsnb(j)
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) topos = topos + nrexpolin(j) + nrexpolnb(j)
    if (thr_limits(ixes(3),tpi).eq.1) then
      if (tpi.eq.1) then
        tmpia = fif
      else
        tmpia = fif + topos
        if (fif.ne.lowia) then
!         clean-up
          if (rs_nbl(j)%nnbtrs.gt.0) then
            thr_limits(ixes(3),tpi) = 2
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%ngnbtrs.gt.0) then
            thr_limits(ixes(3),tpi) = 3
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbtrs.gt.0) then
            thr_limits(ixes(3),tpi) = 4
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 5
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.2) then
      tmpia = fif + topos
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%nnbtrs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%nbtr(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%nnbtrs) then
          if (rs_nbl(j)%ngnbtrs.gt.0) then
            thr_limits(ixes(3),tpi) = 3
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbtrs.gt.0) then
            thr_limits(ixes(3),tpi) = 4
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 5
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.3) then
      tmpia = fif + topos + at(j)%na*rs_nbl(j)%nnbtrats
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%ngnbtrs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%gnbtr(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%ngnbtrs) then
          if (rs_nbl(j)%nwnbtrs.gt.0) then
            thr_limits(ixes(3),tpi) = 4
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 5
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
     end if
   else if (thr_limits(ixes(3),tpi).eq.4) then
      tmpia = fif + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats)
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%nwnbtrs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%wnbtr(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%nwnbtrs) then
          if (rs_nbl(j)%nnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 5
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.5) then
      tmpia = fif + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats)
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%nnbs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%nb(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%nnbs) then
          if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.6) then
      tmpia = fif + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats+rs_nbl(j)%nnbats)
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%ngnbs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%gnb(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%ngnbs) then
          if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.7) then
      tmpia = fif + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats+&
 &                                                     rs_nbl(j)%nnbats+rs_nbl(j)%ngnbats)
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%nwnbs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%wnb(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%nwnbs) then
          if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    end if
    if (thr_limits(ixes(3),tpi).eq.1) thr_limits(ixes(5),tpi) = 1
!
!   the ending bound is kept as is
    j = thr_limits(ixes(2),tpi)
    topos = nrsintra(j) + nrsnb(j)
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) topos = topos + nrexpolin(j) + nrexpolnb(j)
    if (thr_limits(ixes(4),tpi).eq.1) then
      tmpia = laf + topos
    else if (thr_limits(ixes(4),tpi).eq.2) then
      tmpia = laf + topos
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%nnbtrs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%nbtr(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.3) then
      tmpia = laf + topos + at(j)%na*rs_nbl(j)%nnbtrats
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%ngnbtrs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%gnbtr(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.4) then
      tmpia = laf + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats)
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%nwnbtrs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%wnbtr(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.5) then
      tmpia = laf + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats)
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%nnbs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%nb(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.6) then
      tmpia = laf + topos + at(j)%na*(rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats+rs_nbl(j)%nnbats)
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%ngnbs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%gnb(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.7) then
      tmpia = laf + topos + at(j)%na*(rs_nbl(j)%nnbats+rs_nbl(j)%ngnbats+rs_nbl(j)%nnbtrats+&
 &rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats)
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%nwnbs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%wnb(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    end if
!
! check only SR lists if neighbor list update was not performed (5: standard, 6: FEG, 7: water)
!
  else
!
    allia = 0
    do i=1,nseq
      allia = allia + (rs_nbl(i)%nnbats+rs_nbl(i)%ngnbats+rs_nbl(i)%nwnbats)*at(i)%na + nrsintra(i) + nrsnb(i)
    end do
!    allia = sum((rs_nbl(1:nseq)%nnbats+rs_nbl(1:nseq)%ngnbats+rs_nbl(1:nseq)%nwnbats)*at(1:nseq)%na) +  &
! &          sum(nrsintra(1:nseq)) + sum(nrsnb(1:nseq))
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
      allia = allia + sum(nrexpolin(1:nseq)) + sum(nrexpolnb(1:nseq))
    end if
!
    myia = (1.0*allia)/(1.0*thrdat%maxn)
    lowia = nint((tpi-1)*myia,KIND=8)
    hiwia = nint((tpi)*myia,KIND=8)
    tmpia = 0
!   first loop: set residue bounds and loop index
    do i=1,nseq
      stbu(1) = thr_limits(ixes(1),tpi)
      stbu(2) = thr_limits(ixes(2),tpi)
      stbu3 = tmpia
      topos = nrsintra(i) + nrsnb(i)
      if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) topos = topos + nrexpolin(i) + nrexpolnb(i)
!     determine which residue and which list we fall into
      if (topos.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + topos
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 1
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 1
        end if
      end if
      if (rs_nbl(i)%nnbs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%nnbats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 5
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 5
        end if
      end if
      if (rs_nbl(i)%ngnbs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%ngnbats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 6
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 6
        end if
      end if
      if (rs_nbl(i)%nwnbs.gt.0) then
        if (fiia.le.0) fiia = i
        laia = i
        tmpia = tmpia + rs_nbl(i)%nwnbats*at(i)%na
        if ((tmpia.gt.lowia).AND.(thr_limits(ixes(1),tpi).le.0)) then
          thr_limits(ixes(1),tpi) = i
          thr_limits(ixes(3),tpi) = 7
        end if
        if ((tmpia.ge.nint(tpi*myia,KIND=8)).AND.(thr_limits(ixes(2),tpi).le.-1)) then
          thr_limits(ixes(2),tpi) = i
          thr_limits(ixes(4),tpi) = 7
        end if
      end if
!     store interaction cumsum for first and last residues
      if (thr_limits(ixes(1),tpi).ne.stbu(1)) fif = stbu3
      if (thr_limits(ixes(2),tpi).ne.stbu(2)) laf = stbu3
    end do
!   if no interactions whatsoever, invalidate residue bounds and exit
    if (fiia.le.0) then
      thr_limits(ixes(1),tpi) = 1
      thr_limits(ixes(2),tpi) = 0
      return
    end if
    if (tpi.eq.1) then
      thr_limits(ixes(1),tpi) = 1 ! fiia
    end if
    if (tpi.eq.thrdat%maxn) then
      thr_limits(ixes(2),tpi) = nseq ! max(nseq-1,laia) ! laia
    end if
!
!   the starting bound is adjusted (+1) in the innermost bound, thus requiring clean-up
!   if outer bound(s) shift
    j = thr_limits(ixes(1),tpi)
    topos = nrsintra(j) + nrsnb(j)
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) topos = topos + nrexpolin(j) + nrexpolnb(j)
    if (thr_limits(ixes(3),tpi).eq.1) then
      if (tpi.eq.1) then
        tmpia = fif
      else
        tmpia = fif + topos
        if (fif.ne.lowia) then
!         clean-up
          if (rs_nbl(j)%nnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 5
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.5) then
      tmpia = fif + topos
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%nnbs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%nb(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%nnbs) then
          if (rs_nbl(j)%ngnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 6
            thr_limits(ixes(5),tpi) = 1
          else if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.6) then
      tmpia = fif + topos + at(j)%na*rs_nbl(j)%nnbats
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%ngnbs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%gnb(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%ngnbs) then
          if (rs_nbl(j)%nwnbs.gt.0) then
            thr_limits(ixes(3),tpi) = 7
            thr_limits(ixes(5),tpi) = 1
          else if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    else if (thr_limits(ixes(3),tpi).eq.7) then
      tmpia = fif + topos + at(j)%na*(rs_nbl(j)%nnbats+rs_nbl(j)%ngnbats)
      if (tmpia.eq.lowia) then
        thr_limits(ixes(5),tpi) = 1
      else
        thr_limits(ixes(5),tpi) = -1
        do i=1,rs_nbl(j)%nwnbs
          tmpia = tmpia + at(j)%na*at(rs_nbl(j)%wnb(i))%na
          if (tmpia.ge.lowia) then
            thr_limits(ixes(5),tpi) = i+1
            exit
          end if
        end do
!       clean-up
        if (thr_limits(ixes(5),tpi).gt.rs_nbl(j)%nwnbs) then
          if (j.lt.nseq) then
            thr_limits(ixes(1),tpi) = j+1
            thr_limits(ixes(3),tpi) = 1
          end if
        end if
      end if
    end if
    if (thr_limits(ixes(3),tpi).eq.1) thr_limits(ixes(5),tpi) = 1
!
!   the ending bound is kept as is
    j = thr_limits(ixes(2),tpi)
    topos = nrsintra(j) + nrsnb(j)
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) topos = topos + nrexpolin(j) + nrexpolnb(j)
    if (thr_limits(ixes(4),tpi).eq.1) then
      tmpia = laf + topos
    else if (thr_limits(ixes(4),tpi).eq.5) then
      tmpia = laf + topos
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%nnbs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%nb(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.6) then
      tmpia = laf + topos + at(j)%na*rs_nbl(j)%nnbats
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%ngnbs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%gnb(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    else if (thr_limits(ixes(4),tpi).eq.7) then
      tmpia = laf + topos + at(j)%na*(rs_nbl(j)%nnbats+rs_nbl(j)%ngnbats)
      thr_limits(ixes(6),tpi) = -1
      do i=1,rs_nbl(j)%nwnbs
        tmpia = tmpia + at(j)%na*at(rs_nbl(j)%wnb(i))%na
        if (tmpia.ge.nint(tpi*myia,KIND=8)) then
          thr_limits(ixes(6),tpi) = i
          exit
        end if
      end do
    end if
  end if
!$OMP BARRIER
!$OMP SINGLE
!  do i=1,thrdat%maxn
!    write(*,55) i,'X ',thr_limits(ixes(:),i)
!  end do
  k = 1
! recover that first thread always does first interaction
  do k=1,thrdat%maxn
    if (k.gt.1) then
      if ((thr_limits(ixes(1),k).eq.1).AND.(thr_limits(ixes(3),k).eq.1)) then
        thr_limits(ixes(3),k) = thr_limits(ixes(3),k) + 1
      end if
    end if
  end do
! repair empty sets
  do k=1,thrdat%maxn
    if (thr_limits(ixes(1),k).gt.thr_limits(ixes(2),k)) then
      if (k.eq.1) then
        thr_limits(ixes(1:2),k) = 1
        thr_limits(ixes(3:4),k) = 1
        thr_limits(ixes(5),k) = 1
        thr_limits(ixes(6),k) = -1
      else
        thr_limits(ixes(1:2),k) = thr_limits(ixes(2),k-1)
        if (thr_limits(ixes(4),k-1).gt.1) then
          thr_limits(ixes(3:4),k) = thr_limits(ixes(4),k-1)
        else if (fnbl.EQV..true.) then
          thr_limits(ixes(3:4),k) = 2
        else
          thr_limits(ixes(3:4),k) = 5
        end if
        thr_limits(ixes(5),k) = 1
        thr_limits(ixes(6),k) = 0
      end if
    end if
    if ((thr_limits(ixes(1),k).eq.thr_limits(ixes(2),k)).AND.(thr_limits(ixes(3),k).gt.thr_limits(ixes(4),k))) then
      if (k.eq.1) then
        thr_limits(ixes(3:4),k) = 1
        thr_limits(ixes(5),k) = 1
        thr_limits(ixes(6),k) = -1
      else
        if (thr_limits(ixes(4),k-1).gt.1) then
          thr_limits(ixes(3:4),k) = thr_limits(ixes(4),k-1)
        else if (fnbl.EQV..true.) then
          thr_limits(ixes(3:4),k) = 2
        else
          thr_limits(ixes(3:4),k) = 5
        end if
        thr_limits(ixes(5),k) = 1
        thr_limits(ixes(6),k) = 0
      end if
    end if
  end do
!  do i=1,thrdat%maxn
!    write(*,55) i,'Y ',thr_limits(ixes(:),i)
!  end do
  k = 1
! check rest
  do while (k.le.thrdat%maxn)
    i = k
    if (k.eq.1) then
      k = k + 1
      cycle
    end if
    if (thr_limits(ixes(1),k).lt.thr_limits(ixes(2),k-1)) then
      write(ilog,*) 'Fatal. Incoherent residue limits in get_thread_loop_bounds2(...). This is a bug.'
      call fexit()
    end if
    if ((thr_limits(ixes(3),k).eq.thr_limits(ixes(4),i)).AND.(thr_limits(ixes(1),k).eq.thr_limits(ixes(2),i)).AND.&
 &      (thr_limits(ixes(3),k).gt.1).AND.(thr_limits(ixes(5),k).gt.thr_limits(ixes(6),k))) then
!     this is a legitimate empty specification -> don't fix further
      thr_limits(ixes(5),k) = 1
      thr_limits(ixes(6),k) = 0
      k = k + 1
      cycle
    end if
    ii = i
    do while ((thr_limits(ixes(1),ii).eq.thr_limits(ixes(1),ii-1)).AND.(thr_limits(ixes(2),ii).eq.thr_limits(ixes(2),ii-1)).AND.&
 &      (thr_limits(ixes(3),ii).eq.thr_limits(ixes(3),ii-1)).AND.(thr_limits(ixes(4),ii).eq.thr_limits(ixes(4),ii-1)).AND.&
 &      (thr_limits(ixes(5),ii).eq.thr_limits(ixes(5),ii-1)).AND.(thr_limits(ixes(6),ii).eq.thr_limits(ixes(6),ii-1)))
       ii = ii + 1
       if (ii.gt.thrdat%maxn) exit
    end do
    if (ii-1.ge.i) then
      do j=i,ii-1
        if ((thr_limits(ixes(3),k-1).eq.1).AND.(thr_limits(ixes(4),k-1).eq.1)) then
          if (fnbl.EQV..true.) then
            thr_limits(ixes(3:4),j) = 2
          else
            thr_limits(ixes(3:4),j) = 5
          end if
        end if
        thr_limits(ixes(5),j) = 1
        thr_limits(ixes(6),j) = 0
      end do
    end if
    k = k + 1
  end do
!  do i=1,thrdat%maxn
!    write(*,55) i,'X ',thr_limits(ixes(:),i)
!  end do
!$OMP END SINGLE
!
 55 format(i3,a2,10(i7,1x),g10.4)
!
end
!
!
!---------------------------------------------------------------------------------------------------
!
subroutine check_thread_loop_bounds2(fnbl,ixes)
!
  use threads
  use cutoffs
  use sequen, ONLY: nseq
  use polypep, ONLY: at
  use inter, ONLY: nrsintra,nrsnb,nrexpolin,nrexpolnb
  use energies, ONLY: use_POLAR
!
  implicit none
!
  logical, INTENT(IN):: fnbl
  integer, INTENT(IN):: ixes(6)
!
  integer(KIND=8) tmpia,allia
  integer k,kk,j,i,bdums(19)
!
 55 format(i3,a2,10(i7,1x),g10.4)
!
! check bounds for debugging
  if (fnbl.EQV..true.) then
    tmpia = 0
    do i=1,nseq
      tmpia = tmpia + (rs_nbl(i)%nnbtrats+rs_nbl(i)%ngnbtrats+rs_nbl(i)%nwnbtrats+&
 &             rs_nbl(i)%nnbats+rs_nbl(i)%ngnbats+rs_nbl(i)%nwnbats)*at(i)%na +  &
 &             nrsintra(i) + nrsnb(i)
    end do
!    tmpia = sum((rs_nbl(1:nseq)%nnbtrats+rs_nbl(1:nseq)%ngnbtrats+rs_nbl(1:nseq)%nwnbtrats+&
! &             rs_nbl(1:nseq)%nnbats+rs_nbl(1:nseq)%ngnbats+rs_nbl(1:nseq)%nwnbats)*at(1:nseq)%na) +  &
! &        sum(nrsintra(1:nseq)) + sum(nrsnb(1:nseq))
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
      tmpia = tmpia + sum(nrexpolin(1:nseq)) + sum(nrexpolnb(1:nseq))
    end if
  else
    tmpia = 0
    do i=1,nseq
      tmpia = tmpia + (rs_nbl(i)%nnbats+rs_nbl(i)%ngnbats+rs_nbl(i)%nwnbats)*at(i)%na + nrsintra(i) + nrsnb(i)
    end do
!    tmpia = sum((rs_nbl(1:nseq)%nnbats+rs_nbl(1:nseq)%ngnbats+rs_nbl(1:nseq)%nwnbats)*at(1:nseq)%na) +  &
! &        sum(nrsintra(1:nseq)) + sum(nrsnb(1:nseq))
    if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
      tmpia = tmpia + sum(nrexpolin(1:nseq)) + sum(nrexpolnb(1:nseq))
    end if
  end if
!
  allia = 0
  do i=1,thrdat%maxn
    if (fnbl.EQV..true.) then
      do j=thr_limits(ixes(1),i)+1,thr_limits(ixes(2),i)-1
        allia = allia +  nrsintra(j) + nrsnb(j) + at(j)%na*(rs_nbl(j)%nnbats+rs_nbl(j)%ngnbats+rs_nbl(j)%nwnbats+&
                                                          rs_nbl(j)%nnbtrats+rs_nbl(j)%ngnbtrats+rs_nbl(j)%nwnbtrats)
        if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
          allia = allia + nrexpolin(j) + nrexpolnb(j)
        end if
      end do
    else
      do j=thr_limits(ixes(1),i)+1,thr_limits(ixes(2),i)-1
        allia = allia +  nrsintra(j) + nrsnb(j) + at(j)%na*(rs_nbl(j)%nnbats+rs_nbl(j)%ngnbats+rs_nbl(j)%nwnbats)
        if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
          allia = allia + nrexpolin(j) + nrexpolnb(j)
        end if
      end do
    end if
    do kk=1,2
      if ((kk.eq.2).AND.(thr_limits(ixes(1),i).eq.thr_limits(ixes(2),i))) exit
      k = thr_limits(ixes(kk),i)
      call get_thread_loop_inst(i,fnbl,k,ixes,bdums)
      do j=bdums(1),bdums(2)
        allia = allia + at(k)%na*at(rs_nbl(k)%nbtr(j))%na
      end do
      do j=bdums(3),bdums(4)
        allia = allia + at(k)%na*at(rs_nbl(k)%gnbtr(j))%na
      end do
      do j=bdums(5),bdums(6)
        allia = allia + at(k)%na*at(rs_nbl(k)%wnbtr(j))%na
      end do
      do j=bdums(7),bdums(8)
        allia = allia + at(k)%na*at(rs_nbl(k)%nb(j))%na
      end do
      do j=bdums(9),bdums(10)
        allia = allia + at(k)%na*at(rs_nbl(k)%gnb(j))%na
      end do
      do j=bdums(11),bdums(12)
        allia = allia + at(k)%na*at(rs_nbl(k)%wnb(j))%na
      end do
      if (bdums(19).eq.1) then
        allia = allia + nrsintra(k) + nrsnb(k)
        if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) then
          allia = allia + nrexpolin(k) + nrexpolnb(k)
        end if
      end if
    end do
  end do
  if (allia.ne.tmpia) then
    if (fnbl.EQV..true.) then
      do i=1,nseq
        write(*,55) i,': ',rs_nbl(i)%nnbtrs,rs_nbl(i)%ngnbtrs,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nnbs,rs_nbl(i)%ngnbs,rs_nbl(i)%nwnbs
      end do
    else
      do i=1,nseq
        write(*,55) i,': ',rs_nbl(i)%nnbs,rs_nbl(i)%ngnbs,rs_nbl(i)%nwnbs
      end do
    end if
    do i=1,thrdat%maxn
      write(*,55) i,': ',thr_limits(ixes(1),i),thr_limits(ixes(3),i),thr_limits(ixes(5),i),&
 &             thr_limits(ixes(2),i),thr_limits(ixes(4),i),thr_limits(ixes(6),i),rs_nbl(thr_limits(ixes(2),i))%nnbs,tmpia,allia
    end do
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------------------------------------------------
!
! for the multi-level NB-list setup, do DLB
!
subroutine shift_thread_loop_bounds2(thrrat,fracl,fnbl,ixes)
!
  use threads
  use cutoffs
  use iounit, ONLY: ilog
  use sequen, ONLY: nseq
!
  implicit none
!
  logical, INTENT(IN):: fnbl
  integer, INTENT(IN):: ixes(6)
  RTYPE, INTENT(IN):: fracl(thrdat%maxn),thrrat
!
  integer fixed_inc,tpi,inc,k,klo,khi,maxv(7)
  RTYPE rinc
  logical emptyset
!
! lower bound 
!
 45 format(100(i3,1x))
  fixed_inc = 500
!
  do tpi=2,thrdat%maxn
    rinc = 0.5*nseq*min(nseq,50)*(fracl(tpi)-fracl(tpi-1))/thrdat%maxn
    if (rinc.le.(-1*fixed_inc)) then
      inc = -1*fixed_inc
    else if (rinc.gt.fixed_inc) then
      inc = fixed_inc
    else
      inc = nint(rinc)
    end if
!
    if (inc.lt.0) then ! tpi-1 has done more work
!
      maxv(1) = 0
      if (fnbl.EQV..true.) then
        maxv(2) = rs_nbl(thr_limits(ixes(2),tpi-1))%nnbtrs
        maxv(3) = rs_nbl(thr_limits(ixes(2),tpi-1))%ngnbtrs
        maxv(4) = rs_nbl(thr_limits(ixes(2),tpi-1))%nwnbtrs
      else
        maxv(2:4) = 0 
      end if
      maxv(5) = rs_nbl(thr_limits(ixes(2),tpi-1))%nnbs
      maxv(6) = rs_nbl(thr_limits(ixes(2),tpi-1))%ngnbs
      maxv(7) = rs_nbl(thr_limits(ixes(2),tpi-1))%nwnbs
      if (tpi.gt.2) then
        if (thr_limits(ixes(2),tpi-1).eq.thr_limits(ixes(1),tpi-1)) then
          if (thr_limits(ixes(2),tpi-2).eq.thr_limits(ixes(1),tpi-1)) then
            klo = thr_limits(ixes(4),tpi-2)+1
          else
            klo = 2
          end if
        else
          klo = 2
        end if
      else
        klo = 2
      end if
      emptyset = .false.
      if ((thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi)).AND.(thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi)).AND.&
 &        (thr_limits(ixes(3),tpi).gt.1)) then
        emptyset = .true.
      end if
!
!     zeroth case: single topo in single residue
      if ((thr_limits(ixes(3),tpi-1).eq.thr_limits(ixes(4),tpi-1)).AND.&
 &        (thr_limits(ixes(2),tpi-1).eq.thr_limits(ixes(1),tpi-1)).AND.(thr_limits(ixes(3),tpi-1).eq.1)) then
!        write(*,*) tpi,' n0'
        if (thr_limits(ixes(2),tpi-1).gt.1) then
          thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1)
          thr_limits(ixes(3),tpi) = 1 
          thr_limits(ixes(5),tpi) = 1
          if (emptyset.EQV..true.) then ! if the prior set was empty, readjust also terminal bounds
            thr_limits(ixes(2),tpi) = thr_limits(ixes(1),tpi)
            thr_limits(ixes(4),tpi) = 1 
            thr_limits(ixes(6),tpi) = -1
          end if
!         if possible, create empty set for tpi-1 on last NBL of residue prior
          thr_limits(ixes(3:4),tpi-1) = 7
          thr_limits(ixes(6),tpi-1) = 0
          thr_limits(ixes(5),tpi-1) = 1
          thr_limits(ixes(1:2),tpi-1) =  thr_limits(ixes(1),tpi)-1
        end if
!     first case: single NBL in single residue -> must observe both limits
      else if ((thr_limits(ixes(3),tpi-1).eq.thr_limits(ixes(4),tpi-1)).AND.&
 &             (thr_limits(ixes(2),tpi-1).eq.thr_limits(ixes(1),tpi-1))) then
!        write(*,*) tpi,' n1'
        if ((thr_limits(ixes(6),tpi-1)+inc).ge.thr_limits(ixes(5),tpi-1)) then
          thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1) ! may be redundant
          thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi-1) ! may be redundant
          thr_limits(ixes(6),tpi-1) = thr_limits(ixes(6),tpi-1) + inc
          thr_limits(ixes(5),tpi) = thr_limits(ixes(6),tpi-1) + 1
        else
          if (maxv(thr_limits(ixes(4),tpi-1)).gt.0) then
            thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1) ! may be redundant
            thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi-1) ! may be redundant
            thr_limits(ixes(5),tpi) = 1
          end if
          k = thr_limits(ixes(4),tpi-1)-1
          do while (k.ge.klo) 
            if (maxv(k).gt.0) then
              thr_limits(ixes(3:4),tpi-1) = k
              thr_limits(ixes(5),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = maxv(k)
              exit
            end if
            k = k - 1
          end do
          if (k.lt.klo) then
            if (thr_limits(ixes(3),tpi-1).eq.1) then
              thr_limits(ixes(4),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = -1
            else ! make an empty list (tpi-1 is being disabled)
              thr_limits(ixes(3),tpi-1) = thr_limits(ixes(4),tpi-1)
              thr_limits(ixes(5),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = 0
            end if
          end if
        end if
!     second case: cut within the same list (no topos!) 
      else if ((thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi-1)).AND.&
 &             (thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi-1))) then
!        write(*,*) tpi,' n2'
!       shift up NBL or whole residue 
        if ((thr_limits(ixes(6),tpi-1)+inc).lt.1) then
          thr_limits(ixes(5),tpi) = 1
          k = thr_limits(ixes(4),tpi-1)-1
          do while (k.ge.klo) 
            if (maxv(k).gt.0) then
              thr_limits(ixes(4),tpi-1) = k
              thr_limits(ixes(6),tpi-1) = maxv(k)
              exit
            end if
            k = k - 1
          end do
          if (k.lt.klo) then
            if ((thr_limits(ixes(1),tpi-1).lt.thr_limits(ixes(2),tpi-1)).OR.(thr_limits(ixes(3),tpi-1).eq.1)) then
              thr_limits(ixes(4),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = -1
            else ! make an empty list (tpi-1 is being disabled)
              thr_limits(ixes(3),tpi-1) = thr_limits(ixes(4),tpi-1)
              thr_limits(ixes(5),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = 0
            end if
          end if
!       simply adjust bounds for both
        else
          thr_limits(ixes(6),tpi-1) = thr_limits(ixes(6),tpi-1) + inc
          thr_limits(ixes(5),tpi) = thr_limits(ixes(6),tpi-1) + 1
        end if
!     third case: cut between lists but on same residue -> donate 
      else if (thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi-1)) then
!        write(*,*) tpi,' n3'
        if (thr_limits(ixes(4),tpi-1).eq.1) then
          if (thr_limits(ixes(1),tpi-1).eq.thr_limits(ixes(2),tpi-1)) then ! tpi-1 only had topos
            if (thr_limits(ixes(1),tpi-1).gt.1) then
              thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1)
              thr_limits(ixes(3),tpi) = 1 
              thr_limits(ixes(5),tpi) = 1
              if (emptyset.EQV..true.) then ! if the prior tpi-set was empty, readjust also terminal bounds
                thr_limits(ixes(2),tpi) = thr_limits(ixes(1),tpi)
                thr_limits(ixes(4),tpi) = 1 
                thr_limits(ixes(6),tpi) = -1
              end if
!             if possible, create empty set for tpi-1 on last NBL of residue prior
              thr_limits(ixes(3:4),tpi-1) = 7
              thr_limits(ixes(6),tpi-1) = 0
              thr_limits(ixes(5),tpi-1) = 1
              thr_limits(ixes(1:2),tpi-1) = thr_limits(ixes(1),tpi)-1
            end if
          else
            thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1)
            thr_limits(ixes(3),tpi) = 1 
            thr_limits(ixes(5),tpi) = 1
            if (emptyset.EQV..true.) then ! if the prior set was empty, readjust also terminal bounds
              thr_limits(ixes(2),tpi) = thr_limits(ixes(1),tpi)
              thr_limits(ixes(4),tpi) = 1 
              thr_limits(ixes(6),tpi) = -1
            end if
            thr_limits(ixes(4),tpi-1) = 7
            thr_limits(ixes(6),tpi-1) = rs_nbl(thr_limits(ixes(1),tpi)-1)%nwnbs
            thr_limits(ixes(2),tpi-1) = thr_limits(ixes(1),tpi)-1
          end if
        else
          thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi-1)
          if (maxv(thr_limits(ixes(4),tpi-1))+inc.le.0) then
            k = thr_limits(ixes(4),tpi-1)-1
            do while (k.ge.klo) 
              if (maxv(k).gt.0) then
                thr_limits(ixes(4),tpi-1) = k
                thr_limits(ixes(6),tpi-1) = maxv(k)
                exit
              end if
              k = k - 1
            end do
            if (k.lt.klo) then
              if ((thr_limits(ixes(1),tpi-1).lt.thr_limits(ixes(2),tpi-1)).OR.(thr_limits(ixes(3),tpi-1).eq.1)) then
                thr_limits(ixes(4),tpi-1) = 1
                thr_limits(ixes(6),tpi-1) = -1
              else ! make an empty list (tpi-1 is being disabled)
                thr_limits(ixes(3),tpi-1) = thr_limits(ixes(4),tpi-1)
                thr_limits(ixes(5),tpi-1) = 1
                thr_limits(ixes(6),tpi-1) = 0
              end if
            end if
          else
            thr_limits(ixes(6),tpi-1) = maxv(thr_limits(ixes(4),tpi-1)) + inc
            thr_limits(ixes(5),tpi) = thr_limits(ixes(6),tpi-1) + 1
          end if
        end if
!     fourth case: cut was between residues already -> donate some to tpi
      else
!        write(*,*) tpi,' n4'
        if (maxv(thr_limits(ixes(4),tpi-1))+inc.gt.0) then
          thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1)
          thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi-1)
          thr_limits(ixes(6),tpi-1) = maxv(thr_limits(ixes(4),tpi-1))+inc
          thr_limits(ixes(5),tpi) = thr_limits(ixes(6),tpi-1) + 1
        else if (thr_limits(ixes(4),tpi-1).eq.1) then
          if (thr_limits(ixes(2),tpi-1).gt.1) then
            thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1)
            thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi-1)
            if (thr_limits(ixes(1),tpi-1).eq.thr_limits(ixes(2),tpi-1)) then
              thr_limits(ixes(1:2),tpi-1) = thr_limits(ixes(1:2),tpi-1) - 1
              thr_limits(ixes(4),tpi-1) = 7
              thr_limits(ixes(6),tpi-1) = 0
            else
              thr_limits(ixes(2),tpi-1) = thr_limits(ixes(2),tpi-1) - 1
              thr_limits(ixes(4),tpi-1) = 7
              thr_limits(ixes(6),tpi-1) = rs_nbl(thr_limits(ixes(2),tpi-1))%nwnbs
            end if
          end if
        else
          thr_limits(ixes(1),tpi) = thr_limits(ixes(2),tpi-1)
          thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi-1)
          thr_limits(ixes(5),tpi) = 1
          k = thr_limits(ixes(4),tpi-1)-1
          do while (k.ge.klo) 
            if (maxv(k).gt.0) then
              thr_limits(ixes(4),tpi-1) = k
              thr_limits(ixes(6),tpi-1) = maxv(k)
              exit
            end if
            k = k - 1
          end do
          if (k.lt.klo) then
            if ((thr_limits(ixes(1),tpi-1).lt.thr_limits(ixes(2),tpi-1)).OR.(thr_limits(ixes(3),tpi-1).eq.1)) then
              thr_limits(ixes(4),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = -1
            else ! make an empty list (tpi-1 is being disabled)
              thr_limits(ixes(3),tpi-1) = thr_limits(ixes(4),tpi-1)
              thr_limits(ixes(5),tpi-1) = 1
              thr_limits(ixes(6),tpi-1) = 0
            end if
          end if
        end if
      end if 
!
    else if (inc.gt.0) then ! inc is positive, i.e., tpi has done more work
!
      maxv(1) = 0
      if (fnbl.EQV..true.) then
        maxv(2) = rs_nbl(thr_limits(ixes(1),tpi))%nnbtrs
        maxv(3) = rs_nbl(thr_limits(ixes(1),tpi))%ngnbtrs
        maxv(4) = rs_nbl(thr_limits(ixes(1),tpi))%nwnbtrs
      else
        maxv(2:4) = 0 
      end if
      maxv(5) = rs_nbl(thr_limits(ixes(1),tpi))%nnbs
      maxv(6) = rs_nbl(thr_limits(ixes(1),tpi))%ngnbs
      maxv(7) = rs_nbl(thr_limits(ixes(1),tpi))%nwnbs
      khi = 7
      if (tpi.lt.thrdat%maxn) then
        if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) then
          if (thr_limits(ixes(1),tpi+1).eq.thr_limits(ixes(2),tpi)) then
            khi = thr_limits(ixes(3),tpi+1)-1
          else
            khi = 7
          end if
        else
          khi = 7
        end if
      end if
!     zeroth case: single topos
      if ((thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi)).AND.(thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)).AND.&
 &        (thr_limits(ixes(3),tpi).eq.1)) then
!        write(*,*) tpi,' p0'
        thr_limits(ixes(2),tpi-1) = thr_limits(ixes(1),tpi)
        thr_limits(ixes(4),tpi-1) = 1 
        thr_limits(ixes(6),tpi-1) = -1
        k = thr_limits(ixes(3),tpi)+1
        do while (k.le.khi) 
          if (maxv(k).gt.0) then
            thr_limits(ixes(3:4),tpi) = k
            thr_limits(ixes(5),tpi) = 1
            thr_limits(ixes(6),tpi) = 0
            exit
          end if
          k = k + 1
        end do
        if (k.gt.khi) then
          if (fnbl.EQV..true.) then
            thr_limits(ixes(3:4),tpi) = 2
          else
            thr_limits(ixes(3:4),tpi) = 5
          end if
          thr_limits(ixes(5),tpi) = 1
          thr_limits(ixes(6),tpi) = 0
        end if
!     first case: single set in single residue -> must observe both limits
      else if ((thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi)).AND.(thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi))) then
!        write(*,*) tpi,' p1'
        if (((thr_limits(ixes(5),tpi)+inc).le.thr_limits(ixes(6),tpi))) then
          thr_limits(ixes(2),tpi-1) = thr_limits(ixes(1),tpi) ! may be redundant
          thr_limits(ixes(4),tpi-1) = thr_limits(ixes(3),tpi) ! may be redundant
          thr_limits(ixes(5),tpi) = thr_limits(ixes(5),tpi)+inc
          thr_limits(ixes(6),tpi-1) = thr_limits(ixes(5),tpi) - 1
        else
          if (maxv(thr_limits(ixes(3),tpi)).gt.0) then
            thr_limits(ixes(2),tpi-1) = thr_limits(ixes(1),tpi) ! may be redundant
            thr_limits(ixes(4),tpi-1) = thr_limits(ixes(3),tpi) ! may be redundant
            thr_limits(ixes(6),tpi-1) = maxv(thr_limits(ixes(4),tpi-1))
          end if
!         must create empty set for tpi (in-place)
          thr_limits(ixes(5),tpi) = 1
          thr_limits(ixes(6),tpi) = 0
        end if
!     second case: cut within the same list (no topos!) 
      else if ((thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi-1)).AND.&
 &             (thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi-1))) then
!        write(*,*) tpi,' p2'
!       shift up NBL or whole residue 
        if ((thr_limits(ixes(5),tpi)+inc).gt.maxv(thr_limits(ixes(3),tpi))) then
          thr_limits(ixes(4),tpi-1) = thr_limits(ixes(3),tpi)
          thr_limits(ixes(6),tpi-1) = maxv(thr_limits(ixes(4),tpi-1))
          k = thr_limits(ixes(3),tpi)+1
          do while (k.le.khi) 
            if (maxv(k).gt.0) then
              thr_limits(ixes(3),tpi) = k
              thr_limits(ixes(5),tpi) = 1
              exit
            end if
            k = k + 1
          end do
          if (k.gt.khi) then
            if (thr_limits(ixes(2),tpi).gt.thr_limits(ixes(1),tpi)) then
              thr_limits(ixes(1),tpi) = thr_limits(ixes(1),tpi) + 1
              thr_limits(ixes(3),tpi) = 1
              thr_limits(ixes(5),tpi) = 1
              if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) thr_limits(ixes(6),tpi) = -1
            else ! create an empty list in-place
              thr_limits(ixes(4),tpi) = thr_limits(ixes(3),tpi)
              thr_limits(ixes(5),tpi) = 1
              thr_limits(ixes(6),tpi) = 0
            end if
          end if
!       simply adjust bounds for both
        else
          thr_limits(ixes(5),tpi) = thr_limits(ixes(5),tpi)+inc
          thr_limits(ixes(6),tpi-1) = thr_limits(ixes(5),tpi) - 1
        end if
!     third case: cut between lists but on same residue -> donate 
      else if (thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi-1)) then
!        write(*,*) tpi,' p3'
        thr_limits(ixes(4),tpi-1) = thr_limits(ixes(3),tpi)
        thr_limits(ixes(6),tpi-1) = min(inc,maxv(thr_limits(ixes(3),tpi)))
        if (inc.ge.maxv(thr_limits(ixes(3),tpi))) then
          k = thr_limits(ixes(3),tpi) + 1
          do while (k.le.khi) 
            if (maxv(k).gt.0) then
              thr_limits(ixes(3),tpi) = k
              thr_limits(ixes(5),tpi) = 1
              exit
            end if
            k = k + 1
          end do
          if (k.gt.khi) then
            if (thr_limits(ixes(2),tpi).gt.thr_limits(ixes(1),tpi)) then
              thr_limits(ixes(1),tpi) = thr_limits(ixes(1),tpi) + 1
              thr_limits(ixes(3),tpi) = 1
              thr_limits(ixes(5),tpi) = 1
              if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) thr_limits(ixes(6),tpi) = -1
            else
              thr_limits(ixes(4),tpi) = thr_limits(ixes(3),tpi) ! create an empty list on same place
              thr_limits(ixes(5),tpi) = 1
              thr_limits(ixes(6),tpi) = 0
            end if
          end if
        else
          thr_limits(ixes(5),tpi) = inc+1
        end if
!     fourth case: cut was between residues already -> donate topos to tpi-1
      else
!        write(*,*) tpi,' p4'
        thr_limits(ixes(2),tpi-1) = thr_limits(ixes(1),tpi)
        thr_limits(ixes(4),tpi-1) = 1 
        thr_limits(ixes(6),tpi-1) = -1
        k = 2
        do while (k.le.khi) 
          if (maxv(k).gt.0) then
            thr_limits(ixes(3),tpi) = k
            thr_limits(ixes(5),tpi) = 1
            exit
          end if
          k = k + 1
        end do
        if (k.gt.khi) then
          if (thr_limits(ixes(2),tpi).gt.thr_limits(ixes(1),tpi)) then
            thr_limits(ixes(1),tpi) = thr_limits(ixes(1),tpi) + 1
            thr_limits(ixes(3),tpi) = 1
            thr_limits(ixes(5),tpi) = 1
            if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) thr_limits(ixes(6),tpi) = -1
          else
            if (fnbl.EQV..true.) then
              thr_limits(ixes(3),tpi) = 2
            else
              thr_limits(ixes(3),tpi) = 5
            end if
            thr_limits(ixes(5),tpi) = 1
            if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) thr_limits(ixes(6),tpi) = 0
          end if
        end if
      end if 
    end if
!
    if ((thr_limits(ixes(1),tpi-1).eq.thr_limits(ixes(2),tpi-1)).AND.(thr_limits(ixes(3),tpi-1).gt.thr_limits(ixes(4),tpi-1))) then
      write(ilog,*) 'Fatal. Subroutine shift_thread_loop_bounds2(...) created an incorrect empty set for thread ',tpi-1,'. This is &
 &a bug.'
      call fexit()
    end if
    if (tpi.ge.3) then
      if ((thr_limits(ixes(1),tpi-1).eq.thr_limits(ixes(2),tpi-2)).AND.&
 &(thr_limits(ixes(3),tpi-1).lt.thr_limits(ixes(4),tpi-2))) then
      write(ilog,*) 'Fatal. Subroutine shift_thread_loop_bounds2(...) created overlapping sets for threads ',tpi-1,'/',tpi-2,&
 &'. This is a bug.'
      call fexit()
      end if
    end if
!
!   correct lower list bound where needed
    do k=tpi-2,1,-1
      if ((thr_limits(ixes(4),k).eq.thr_limits(ixes(3),tpi-1)).AND.&
 & (thr_limits(ixes(2),k).eq.thr_limits(ixes(1),tpi-1)).AND.(thr_limits(ixes(6),k).ge.thr_limits(ixes(5),tpi-1))) then
        thr_limits(ixes(5),tpi-1) = thr_limits(ixes(6),k) + 1
        exit
      else if (thr_limits(ixes(2),k).lt.thr_limits(ixes(1),tpi-1)) then
        exit
      end if
    end do
!
!   check single list bounds
    if ((thr_limits(ixes(3),tpi-1).eq.thr_limits(ixes(4),tpi-1)).AND.&
 &      (thr_limits(ixes(2),tpi-1).eq.thr_limits(ixes(1),tpi-1))) then
      if (thr_limits(ixes(6),tpi-1).lt.thr_limits(ixes(5),tpi-1)) then
        if (thr_limits(ixes(3),tpi-1).eq.1) then
          thr_limits(ixes(6),tpi-1) = -1
          thr_limits(ixes(5),tpi-1) = 1
        else
          thr_limits(ixes(6),tpi-1) = 0
          thr_limits(ixes(5),tpi-1) = 1
        end if
        cycle
      end if
    end if
!
    do k=tpi-2,1,-1
      if ((thr_limits(ixes(3),k).eq.thr_limits(ixes(4),k)).AND.(thr_limits(ixes(1),k).eq.thr_limits(ixes(2),k)).AND.&
 &        (thr_limits(ixes(6),k).lt.thr_limits(ixes(5),k))) cycle
      if (thr_limits(ixes(2),k).lt.thr_limits(ixes(1),tpi-1)) exit
      if ((thr_limits(ixes(3),tpi-1).eq.thr_limits(ixes(4),k)).AND.(thr_limits(ixes(1),tpi-1).eq.thr_limits(ixes(2),k)).AND.&
 &        (thr_limits(ixes(5),tpi-1).lt.thr_limits(ixes(6),k))) then
        do klo=1,tpi-1
          write(*,45) thr_limits(ixes(1:6),klo)
        end do
        write(ilog,*) 'Fatal. Subroutine shift_thread_loop_bounds2(...) created overlapping sets for threads ',tpi-1,'/',k,&
 &'. This is a bug.'
        call fexit()
      end if
    end do

  end do
!
end
!
!-----------------------------------------------------------------------------------------------------------------------
!
! because NBL-sizes change over time, occasionally this may produce "empty" loop limits
! in unexpected cases
!
subroutine get_thread_loop_inst(tpi,fnbl,ix,ixes,outix)
!
  use threads
  use cutoffs
  use polypep, ONLY: at
!
  implicit none
!
  integer, INTENT(IN):: tpi,ix,ixes(6)
  logical, INTENT(IN):: fnbl
!
  integer i,outix(19),lo,hi
!
  lo = thr_limits(ixes(1),tpi)
  hi = thr_limits(ixes(2),tpi)
  if (fnbl.EQV..true.) then
    outix(:) = 1
    outix(2) = rs_nbl(ix)%nnbtrs
    outix(4) = rs_nbl(ix)%ngnbtrs
    outix(6) = rs_nbl(ix)%nwnbtrs
    outix(8) = rs_nbl(ix)%nnbs
    outix(10) = rs_nbl(ix)%ngnbs
    outix(12) = rs_nbl(ix)%nwnbs
    outix(13) = rs_nbl(ix)%nnbtrats
    outix(14) = rs_nbl(ix)%ngnbtrats
    outix(15) = rs_nbl(ix)%nwnbtrats
    outix(16) = rs_nbl(ix)%nnbats
    outix(17) = rs_nbl(ix)%ngnbats
    outix(18) = rs_nbl(ix)%nwnbats
    if ((ix.eq.lo).AND.(ix.eq.hi)) then
      if (thr_limits(ixes(3),tpi).eq.2) outix(1) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.3) outix(3) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.4) outix(5) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.5) outix(7) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.6) outix(9) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.7) outix(11) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(4),tpi).eq.2) outix(2) = min(outix(2),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.3) outix(4) = min(outix(4),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.4) outix(6) = min(outix(6),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.5) outix(8) = min(outix(8),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.6) outix(10) = min(outix(10),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.7) outix(12) = min(outix(12),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(3),tpi).gt.2) outix(2) = 0
      if (thr_limits(ixes(3),tpi).gt.3) outix(4) = 0
      if (thr_limits(ixes(3),tpi).gt.4) outix(6) = 0
      if (thr_limits(ixes(3),tpi).gt.5) outix(8) = 0
      if (thr_limits(ixes(3),tpi).gt.6) outix(10) = 0
      if (thr_limits(ixes(4),tpi).lt.2) outix(2) = 0
      if (thr_limits(ixes(4),tpi).lt.3) outix(4) = 0
      if (thr_limits(ixes(4),tpi).lt.4) outix(6) = 0
      if (thr_limits(ixes(4),tpi).lt.5) outix(8) = 0
      if (thr_limits(ixes(4),tpi).lt.6) outix(10) = 0
      if (thr_limits(ixes(4),tpi).lt.7) outix(12) = 0
      if (thr_limits(ixes(3),tpi).ne.1) outix(19) = 0
    else if (ix.eq.lo) then
      if (thr_limits(ixes(3),tpi).eq.2) outix(1) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.3) outix(3) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.4) outix(5) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.5) outix(7) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.6) outix(9) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.7) outix(11) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).gt.2) outix(2) = 0
      if (thr_limits(ixes(3),tpi).gt.3) outix(4) = 0
      if (thr_limits(ixes(3),tpi).gt.4) outix(6) = 0
      if (thr_limits(ixes(3),tpi).gt.5) outix(8) = 0
      if (thr_limits(ixes(3),tpi).gt.6) outix(10) = 0
      if (thr_limits(ixes(3),tpi).ne.1) outix(19) = 0
    else if (ix.eq.hi) then
      if (thr_limits(ixes(4),tpi).eq.2) outix(2) = min(outix(2),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.3) outix(4) = min(outix(4),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.4) outix(6) = min(outix(6),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.5) outix(8) = min(outix(8),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.6) outix(10) = min(outix(10),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.7) outix(12) = min(outix(12),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).lt.2) outix(2) = 0
      if (thr_limits(ixes(4),tpi).lt.3) outix(4) = 0
      if (thr_limits(ixes(4),tpi).lt.4) outix(6) = 0
      if (thr_limits(ixes(4),tpi).lt.5) outix(8) = 0
      if (thr_limits(ixes(4),tpi).lt.6) outix(10) = 0
      if (thr_limits(ixes(4),tpi).lt.7) outix(12) = 0
    end if
  else
    outix(1:12) = 1
    outix(2) = 0
    outix(4) = 0
    outix(6) = 0
    outix(19) = 1
    outix(8) = rs_nbl(ix)%nnbs
    outix(10) = rs_nbl(ix)%ngnbs
    outix(12) = rs_nbl(ix)%nwnbs
    outix(13:15) = 0
    outix(16) = rs_nbl(ix)%nnbats
    outix(17) = rs_nbl(ix)%ngnbats
    outix(18) = rs_nbl(ix)%nwnbats
    if ((ix.eq.lo).AND.(ix.eq.hi)) then
      if (thr_limits(ixes(3),tpi).eq.5) outix(7) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.6) outix(9) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.7) outix(11) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(4),tpi).eq.5) outix(8) = min(outix(8),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.6) outix(10) = min(outix(10),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.7) outix(12) = min(outix(12),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(3),tpi).gt.5) outix(8) = 0
      if (thr_limits(ixes(3),tpi).gt.6) outix(10) = 0
      if (thr_limits(ixes(4),tpi).lt.5) outix(8) = 0
      if (thr_limits(ixes(4),tpi).lt.6) outix(10) = 0
      if (thr_limits(ixes(4),tpi).lt.7) outix(12) = 0
      if (thr_limits(ixes(3),tpi).ne.1) outix(19) = 0
    else if (ix.eq.lo) then
      if (thr_limits(ixes(3),tpi).eq.5) outix(7) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.6) outix(9) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).eq.7) outix(11) = thr_limits(ixes(5),tpi)
      if (thr_limits(ixes(3),tpi).gt.5) outix(8) = 0
      if (thr_limits(ixes(3),tpi).gt.6) outix(10) = 0
      if (thr_limits(ixes(3),tpi).ne.1) outix(19) = 0
    else if (ix.eq.hi) then
      if (thr_limits(ixes(4),tpi).eq.5) outix(8) = min(outix(8),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.6) outix(10) = min(outix(10),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).eq.7) outix(12) = min(outix(12),thr_limits(ixes(6),tpi))
      if (thr_limits(ixes(4),tpi).lt.5) outix(8) = 0
      if (thr_limits(ixes(4),tpi).lt.6) outix(10) = 0
      if (thr_limits(ixes(4),tpi).lt.7) outix(12) = 0
    end if
  end if  
!
  if (fnbl.EQV..true.) then
    if ((outix(2).gt.0).AND.((outix(1).ne.1).OR.(outix(2).ne.rs_nbl(ix)%nnbtrs))) then
      outix(13) = 0
      do i=outix(1),outix(2)
        outix(13) = outix(13) + at(rs_nbl(ix)%nbtr(i))%na
      end do
    else
      outix(13) = rs_nbl(ix)%nnbtrats
    end if
    if ((outix(4).gt.0).AND.((outix(3).ne.1).OR.(outix(4).ne.rs_nbl(ix)%ngnbtrs))) then
      outix(14) = 0
      do i=outix(3),outix(4)
        outix(14) = outix(14) + at(rs_nbl(ix)%gnbtr(i))%na
      end do
    else
      outix(14) = rs_nbl(ix)%ngnbtrats
    end if
    if ((outix(6).gt.0).AND.((outix(5).ne.1).OR.(outix(6).ne.rs_nbl(ix)%nwnbtrs))) then
      outix(15) = 0
      do i=outix(5),outix(6)
        outix(15) = outix(15) + at(rs_nbl(ix)%wnbtr(i))%na
      end do
    else
      outix(15) = rs_nbl(ix)%nwnbtrats
    end if
  end if
  if ((outix(8).gt.0).AND.((outix(7).ne.1).OR.(outix(8).ne.rs_nbl(ix)%nnbs))) then
    outix(16) = 0
    do i=outix(7),outix(8)
      outix(16) = outix(16) + at(rs_nbl(ix)%nb(i))%na
    end do
  else
    outix(16) = rs_nbl(ix)%nnbats
  end if
  if ((outix(10).gt.0).AND.((outix(9).ne.1).OR.(outix(10).ne.rs_nbl(ix)%ngnbs))) then
    outix(17) = 0
    do i=outix(9),outix(10)
      outix(17) = outix(17) + at(rs_nbl(ix)%gnb(i))%na
    end do
  else
    outix(17) = rs_nbl(ix)%ngnbats
  end if
  if ((outix(12).gt.0).AND.((outix(11).ne.1).OR.(outix(12).ne.rs_nbl(ix)%nwnbs))) then
    outix(18) = 0
    do i=outix(11),outix(12)
      outix(18) = outix(18) + at(rs_nbl(ix)%wnb(i))%na
    end do
  else
    outix(18) = rs_nbl(ix)%nwnbats
  end if
!
end
! 
!-----------------------------------------------------------------------------------------------------------
!
! the problem with NB-loop bounds is that the possible extent for the next step is not
! known until right before the loop is encountered as NBLs can change size and, in particular
! switch from empty to non-empty and vice versa
!
! it is thus necessary to constantly adjust the boudns in ixes to avoid errors and confusion
! generally, we adjust upper limits to maximally allowed range given fixed starting points for next
! thread (this can be a reduction!) first
! then, there is a barrier
! finally, lower bounds are pulled in as needed
!
subroutine fix_thread_loop_bounds2(fnbl,ixes,tpi)
!
  use threads
  use cutoffs
  use sequen, ONLY: nseq
  use iounit, ONLY: ilog
!
  implicit none
!
  integer, INTENT(IN):: tpi,ixes(6)
  logical, INTENT(IN):: fnbl
!
  integer k,k2,k3,maxv(7)
  logical needend,emptyset
!
  emptyset = .false.
  if ((thr_limits(ixes(2),tpi).gt.0).AND.(thr_limits(ixes(2),tpi).le.nseq)) then
    maxv(1) = 0
    if (fnbl.EQV..true.) then
      maxv(2) = rs_nbl(thr_limits(ixes(2),tpi))%nnbtrs
      maxv(3) = rs_nbl(thr_limits(ixes(2),tpi))%ngnbtrs
      maxv(4) = rs_nbl(thr_limits(ixes(2),tpi))%nwnbtrs
    else
      maxv(2:4) = 0 
    end if
    maxv(5) = rs_nbl(thr_limits(ixes(2),tpi))%nnbs
    maxv(6) = rs_nbl(thr_limits(ixes(2),tpi))%ngnbs
    maxv(7) = rs_nbl(thr_limits(ixes(2),tpi))%nwnbs
  else
    write(ilog,*) 'Fatal. Encountered an illegal last residue # for thread ',tpi,' (',thr_limits(ixes(2),tpi),&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
    call fexit()
  end if
!
! shrink bound
  if (thr_limits(ixes(4),tpi).gt.1) then
!    thr_limits(ixes(6),tpi) = min(thr_limits(ixes(6),tpi),maxv(thr_limits(ixes(4),tpi)))
  end if
!
  emptyset = .false.
  if ((thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)).AND.(thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi))) then
    if ((thr_limits(ixes(3),tpi).gt.1).AND.(thr_limits(ixes(6),tpi).lt.thr_limits(ixes(5),tpi))) then
      if (tpi.eq.thrdat%maxn) then
        emptyset = .true.
      else
        if (thr_limits(ixes(1),tpi+1).eq.(thr_limits(ixes(2),tpi)+1)) then
          emptyset = .true.
        else if (thr_limits(ixes(3),tpi+1).eq.(thr_limits(ixes(4),tpi)+1)) then
          emptyset = .true.
        end if
      end if
    end if
  end if
!  if (emptyset.EQV..true.) write(*,*) tpi,' is eligible and empty.'
  needend = .true.
  if (tpi.lt.thrdat%maxn) then
    if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi+1)) then
      if (thr_limits(ixes(4),tpi).eq.thr_limits(ixes(3),tpi+1)) then
        if (thr_limits(ixes(4),tpi).eq.1) then
          write(ilog,*) 'Fatal. Encountered overlapping sets (3) for threads ',tpi,'/',tpi+1,&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
          call fexit()
        end if
        if ((maxv(thr_limits(ixes(4),tpi)).gt.0).AND.(thr_limits(ixes(6),tpi).gt.0)) then
          thr_limits(ixes(6),tpi) = min(thr_limits(ixes(6),tpi),maxv(thr_limits(ixes(4),tpi)))
          needend = .false.
        end if
      else if (thr_limits(ixes(4),tpi).eq.(thr_limits(ixes(3),tpi+1)-1)) then
        if (maxv(thr_limits(ixes(4),tpi)).gt.0) then
          thr_limits(ixes(6),tpi) = maxv(thr_limits(ixes(4),tpi))
          needend = .false.
        end if
      else if (thr_limits(ixes(4),tpi).lt.thr_limits(ixes(3),tpi+1)) then
        if ((maxv(thr_limits(ixes(4),tpi)).gt.0).AND.(emptyset.EQV..true.)) then
!         reenable an eligible empty set (i.e., those not followed by another set on the same list) if it has become nonempty
!          write(*,*) 'reenabling for ',tpi,thr_limits(ixes(4),tpi),thr_limits(ixes(6),tpi)
          thr_limits(ixes(6),tpi) = maxv(thr_limits(ixes(4),tpi))
          needend = .false.
        end if
!        otherwise do nothing: we have to redetermine end bounds 4/6
      else
        write(ilog,*) 'Fatal. Encountered overlapping sets (2) for threads ',tpi,'/',tpi+1,&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
        call fexit()
      end if
    else if (thr_limits(ixes(2),tpi).lt.thr_limits(ixes(1),tpi+1)) then
      if ((maxv(thr_limits(ixes(4),tpi)).gt.0).AND.(emptyset.EQV..true.)) then
!       reenable an eligible empty set (i.e., those not followed by another set on the same list) if it has become nonempty
!        write(*,*) 'reenabling2 for ',tpi,thr_limits(ixes(4),tpi),thr_limits(ixes(6),tpi)
        thr_limits(ixes(6),tpi) = maxv(thr_limits(ixes(4),tpi))
        needend = .false.
      end if
!     do nothing: we have to redetermine end bounds 4/6
    else
      write(ilog,*) 'Fatal. Encountered overlapping sets (1) for threads ',tpi,'/',tpi+1,&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
      call fexit()
    end if
  end if
  if (needend.EQV..true.) then
    if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) then
      k2 = max(2,thr_limits(ixes(3),tpi))
    else
      k2 = 2
    end if
    if (tpi.eq.thrdat%maxn) then
      k3 = 7
    else
      if (thr_limits(ixes(2),tpi).lt.thr_limits(ixes(1),tpi+1)) then
        k3 = 7
      else ! implies equal rs-#
        k3 = thr_limits(ixes(3),tpi+1)-1
      end if
    end if
    do k=k3,k2,-1
      if (maxv(k).gt.0) then
        thr_limits(ixes(4),tpi) = k
        thr_limits(ixes(6),tpi) = maxv(thr_limits(ixes(4),tpi))
        needend = .false.
        exit
      end if
    end do
  end if
! still only topos in single res.
  if ((thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)).AND.(thr_limits(ixes(4),tpi).eq.1)) then
    needend = .false.
  end if
  if ((needend.EQV..true.).AND.(thr_limits(ixes(2),tpi).gt.thr_limits(ixes(1),tpi))) then 
    thr_limits(ixes(4),tpi) = 1
    thr_limits(ixes(6),tpi) = -1
    needend = .false.
  end if
  if (needend.EQV..true.) then ! this means no possible or only unexpected interactions for this thread
    thr_limits(ixes(4),tpi) = thr_limits(ixes(3),tpi) ! 1:2 the same is implied
    if (tpi.lt.thrdat%maxn) then
      if (thr_limits(ixes(3),tpi+1).eq.thr_limits(ixes(4),tpi)) then
        if ((thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi)).AND.(thr_limits(ixes(5),tpi).le.thr_limits(ixes(6),tpi))) then
           thr_limits(ixes(6),tpi) = thr_limits(ixes(5),tpi+1)-1
!          thr_limits(ixes(6),tpi) = maxv(thr_limits(ixes(4),tpi))
!          thr_limits(ixes(5),tpi+1) = maxv(thr_limits(ixes(4),tpi)) + 1
        else
          thr_limits(ixes(6),tpi) = 0
        end if
      else
        thr_limits(ixes(6),tpi) = 0
      end if
    else
      thr_limits(ixes(6),tpi) = maxv(thr_limits(ixes(4),tpi))
    end if
  end if
!
!$OMP BARRIER
!
  if ((thr_limits(ixes(1),tpi).gt.0).AND.(thr_limits(ixes(1),tpi).le.nseq)) then
    maxv(1) = 0
    if (fnbl.EQV..true.) then
      maxv(2) = rs_nbl(thr_limits(ixes(1),tpi))%nnbtrs
      maxv(3) = rs_nbl(thr_limits(ixes(1),tpi))%ngnbtrs
      maxv(4) = rs_nbl(thr_limits(ixes(1),tpi))%nwnbtrs
    else
      maxv(2:4) = 0 
    end if
    maxv(5) = rs_nbl(thr_limits(ixes(1),tpi))%nnbs
    maxv(6) = rs_nbl(thr_limits(ixes(1),tpi))%ngnbs
    maxv(7) = rs_nbl(thr_limits(ixes(1),tpi))%nwnbs
  else
    write(ilog,*) 'Fatal. Encountered an illegal first residue # for thread ',tpi,' (',thr_limits(ixes(1),tpi),&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
    call fexit()
  end if

  needend = .true.
  if (tpi.gt.1) then
    if (thr_limits(ixes(1),tpi).gt.thr_limits(ixes(2),tpi-1)) then
      thr_limits(ixes(3),tpi) = 1
      thr_limits(ixes(5),tpi) = 1
      needend = .false.
    else if (thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi-1)) then
      if (thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi-1)) then
        if (thr_limits(ixes(4),tpi-1).eq.1) then
          write(ilog,*) 'Fatal. Encountered overlapping sets (4) for threads ',tpi-1,'/',tpi,&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
          call fexit()
        end if
        if (thr_limits(ixes(6),tpi-1).gt.0) then
          if ((maxv(thr_limits(ixes(3),tpi)).gt.thr_limits(ixes(6),tpi-1))) then
            thr_limits(ixes(5),tpi) = thr_limits(ixes(6),tpi-1) + 1
            needend = .false.
          end if
        else
          k = tpi-1
          do while (k.ge.1)
            if ((thr_limits(ixes(4),k).ne.thr_limits(ixes(3),tpi)).OR.&
 &              (thr_limits(ixes(2),k).ne.thr_limits(ixes(1),tpi))) then
              if (maxv(thr_limits(ixes(3),tpi)).gt.0) then
!                write(*,*) 'hitit',tpi
                thr_limits(ixes(5),tpi) = 1
                needend = .false.
              end if
              exit
            end if
            if (thr_limits(ixes(6),k).gt.0) then
              if (maxv(thr_limits(ixes(3),tpi)).gt.thr_limits(ixes(6),k)) then
                thr_limits(ixes(5),tpi) = thr_limits(ixes(6),k) + 1
                needend = .false.
              end if
              exit
            end if
            k = k - 1
          end do
        end if
      else if (thr_limits(ixes(4),tpi-1).lt.thr_limits(ixes(3),tpi)) then
        k = thr_limits(ixes(3),tpi)
        do while (k.gt.thr_limits(ixes(4),tpi-1))
          if (maxv(k).gt.0) then
            thr_limits(ixes(5),tpi) = 1
            thr_limits(ixes(3),tpi) = k
            needend = .false.
          end if
          k = k - 1
        end do
      else
        write(ilog,*) 'Fatal. Encountered overlapping sets (5) for threads ',tpi-1,'/',tpi,&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
        call fexit()
      end if
    else
      write(ilog,*) 'Fatal. Encountered overlapping sets (6) for threads ',tpi-1,'/',tpi,&
 &') in fix_thread_loop_bounds2(...). This is a bug.'
      call fexit()
    end if
  end if
! only topos in single res.
  if ((thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)).AND.(thr_limits(ixes(4),tpi).eq.1)) then
    needend = .false.
  end if
  if (tpi.eq.1) needend = .false.
  if (needend.EQV..true.) then
    if (thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)) then
      k3 = thr_limits(ixes(4),tpi)
    else
      k3 = 7
    end if
    if (tpi.eq.1) then
      k2 = 2
    else
      if (thr_limits(ixes(1),tpi).eq.thr_limits(ixes(2),tpi-1)) then
        k2 = thr_limits(ixes(4),tpi-1)+1
      end if
    end if
    do k=k2,k3
      if (maxv(k).gt.0) then
        thr_limits(ixes(3),tpi) = k
        thr_limits(ixes(5),tpi) = 1
        needend = .false.
        exit
      end if
    end do
  end if
  if (needend.EQV..true.) then
    if (thr_limits(ixes(1),tpi).gt.thr_limits(ixes(2),tpi-1)) then
      thr_limits(ixes(3),tpi) = 1
      thr_limits(ixes(5),tpi) = 1
      needend = .false.
    end if
  end if
  if ((needend.EQV..true.).AND.(thr_limits(ixes(2),tpi).gt.thr_limits(ixes(1),tpi))) then
    if (thr_limits(ixes(4),tpi-1).eq.7) then ! have to shift first residue bound
      thr_limits(ixes(1),tpi) = thr_limits(ixes(1),tpi) + 1
      thr_limits(ixes(3),tpi) = 1
      thr_limits(ixes(5),tpi) = 1
    else
      thr_limits(ixes(3),tpi) = max(2,thr_limits(ixes(4),tpi-1)+1) ! lowest empty list
      thr_limits(ixes(5),tpi) = 1
    end if
    needend = .false.
  end if
  emptyset = .false.
  if (needend.EQV..true.) then ! this means no possible interactions for this thread
    thr_limits(ixes(3),tpi) = thr_limits(ixes(4),tpi)
    thr_limits(ixes(5),tpi) = 1
    if (tpi.gt.1) then
      if (thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi-1)) emptyset = .true.
    end if
    if (maxv(thr_limits(ixes(3),tpi)).le.0) emptyset = .true.
  end if
!
!$OMP BARRIER
  if (emptyset.EQV..true.) then
!    write(*,*) 'resetting ',tpi
    thr_limits(ixes(6),tpi) = 0
  end if
! ensure consistency in representing topos and empty lists
  if ((thr_limits(ixes(2),tpi).eq.thr_limits(ixes(1),tpi)).AND.(thr_limits(ixes(3),tpi).eq.thr_limits(ixes(4),tpi))) then
    if (thr_limits(ixes(3),tpi).eq.1) then
      thr_limits(ixes(5),tpi) = 1
      thr_limits(ixes(6),tpi) = -1
    else if (thr_limits(ixes(5),tpi).gt.thr_limits(ixes(6),tpi)) then
      thr_limits(ixes(5),tpi) = 1
      thr_limits(ixes(6),tpi) = 0
      if (tpi.eq.thrdat%maxn) then
        do k=7,thr_limits(ixes(4),tpi)+1,-1
          if (maxv(k).gt.0) then
            thr_limits(ixes(4),tpi) = k
            thr_limits(ixes(6),tpi) = maxv(k)
            thr_limits(ixes(5),tpi) = maxv(thr_limits(ixes(3),tpi))+1 ! remain disabled
            exit
          end if
        end do
      else
        if (thr_limits(ixes(1),tpi+1).gt.thr_limits(ixes(2),tpi)) then
          do k=7,thr_limits(ixes(4),tpi)+1,-1
            if (maxv(k).gt.0) then
              thr_limits(ixes(4),tpi) = k
              thr_limits(ixes(6),tpi) = maxv(k)
              thr_limits(ixes(5),tpi) = maxv(thr_limits(ixes(3),tpi))+1 ! remain disabled
              exit
            end if
          end do
        else if (thr_limits(ixes(3),tpi+1).gt.(thr_limits(ixes(4),tpi)+1)) then
          do k=thr_limits(ixes(3),tpi+1)-1,thr_limits(ixes(4),tpi)+1,-1
            if (maxv(k).gt.0) then
              thr_limits(ixes(4),tpi) = k
              thr_limits(ixes(6),tpi) = maxv(k)
              thr_limits(ixes(5),tpi) = maxv(thr_limits(ixes(3),tpi))+1 ! remain disabled
              exit
            end if
          end do
        end if
      end if
    end if
  end if
!
end
!
!-----------------------------------------------------------------------------------------------------------
!
! a routine to be called by all threads to get safe bounds for a specified interval
!
subroutine threads_bounds(lon,hin,tpi,tpn,bnds)
!
  implicit none
!
  integer, INTENT(IN):: lon,hin,tpi,tpn
!
  integer ton,bnds(2),hlp
  RTYPE rhlp
!
  ton = hin - lon + 1
!
  if (tpi.le.0) then
    bnds(1) = lon
    bnds(2) = hin
  else
    if (mod(ton,tpn).eq.0) then
      hlp = ton/tpn
      bnds(1) = (tpi-1)*hlp + lon
      bnds(2) = tpi*hlp + lon - 1
    else if (ton.lt.tpn) then
      if (tpi.gt.ton) then
        bnds(2) = 0
        bnds(1) = 1
      else
        bnds(1:2) = lon + tpi - 1
      end if
    else
      rhlp = (1.0*ton)/(1.0*tpn)
      bnds(1) = nint((tpi-1)*rhlp) + lon
      bnds(2) = nint(tpi*rhlp) + lon - 1
    end if
  end if
!
end
!
!-------------------------------------------------------------------------------------------------------------
!
! a subroutine to deal with holonomic constraint groups that need internal parallelization
!
subroutine threads_setup_shake()
!
  use threads
  use shakeetal
  use iounit
  use sequen
  use atoms
!
  implicit none
!
  integer tpn,i,ii,iii,k,tpx,kk,totc1,totc2,lastval(10),j,jj
  integer, ALLOCATABLE:: maxsn(:),minsn(:),brks(:)
  logical have_warned 
!
  have_warned = .false.
  nccgs = 0
  if (thrdat%suppflag(2).EQV..true.) return
!
  tpn = thrdat%maxn ! -1
  totc1 = 3*(settle_tip3ps+settle_tip4ps+settle_tip4pes+settle_tip5ps+settle_spcs+settle_rest)
  if (shake_cnt.gt.0) then
    totc2 = sum(constraints(1:shake_cnt)%nr+constraints(1:shake_cnt)%nr3+constraints(1:shake_cnt)%nr4)
    do iii=1,2
      do ii=1,shake_cnt
        i = nonsettlelst(ii)
        if (1.0*(constraints(i)%nr+constraints(i)%nr3+constraints(i)%nr4)/(1.0*(totc1+totc2)).gt.(0.1/tpn)) then
          if ((iii.eq.1).AND.(have_warned.EQV..false.)) then
            write(ilog,*) 'Warning. Structure of holonomic constraints is such that a single constraint group &
   &is too large to make parallelization by distributing constraint groups across threads efficient.'
            have_warned = .true.
          end if
          if ((constraints(i)%nr3.le.0).AND.(constraints(i)%nr4.le.0).AND.(constraints(i)%nr.gt.10*tpn)) then
            allocate(maxsn(constraints(i)%nr))
            allocate(minsn(constraints(i)%nr))
            allocate(brks(constraints(i)%nr))
            maxsn(1) = maxval(constraints(i)%idx(1,:))
            do k=2,constraints(i)%nr
              maxsn(k) = max(maxsn(k-1),maxval(constraints(i)%idx(k,:)))
            end do
            minsn(constraints(i)%nr) = minval(constraints(i)%idx(constraints(i)%nr,:))
            do k=constraints(i)%nr-1,1,-1
              minsn(k) = min(minsn(k+1),minval(constraints(i)%idx(k,:)))
            end do
            kk = 0
            do k=2,constraints(i)%nr-1
              if ((minval(constraints(i)%idx(k,:)).le.maxsn(k-1)).AND.(maxval(constraints(i)%idx(k,:)).gt.maxsn(k-1)).AND.&
   &              (minval(constraints(i)%idx(k,:)).lt.minsn(k+1)).AND.(maxval(constraints(i)%idx(k,:)).ge.minsn(k+1)).AND.&
   &              (maxsn(k-1).lt.minsn(k+1)).AND.(constraints(i)%massless(k).EQV..false.))  then
                if (molofrs(atmres(constraints(i)%idx(k,1))).ne.molofrs(atmres(constraints(i)%idx(k,2)))) then
                  if (iii.eq.1) then
                    write(ilog,*) 'Warning. A constraint group is ineligible for internal parallelization due to &
 &the presence of intermolecular constraints. This can cause performance loss and may be fixed in the future.'
                  end if
                  kk = 0
                  exit
                end if
                if (kk.le.0) then
                  kk = kk + 1
                  brks(kk) = k
!               we need to avoid creating small  constraint units that span both directions in their protected set
                else ! if (brks(kk).lt.(k-1)) then ! single constraint must fail
                  j = 1
                  do jj=1,kk
                    j = brks(jj)
                    if ((constraints(i)%idx(k,1).eq.constraints(i)%idx(j,1)).OR.&
 &(constraints(i)%idx(k,1).eq.constraints(i)%idx(j,2)).OR.(constraints(i)%idx(k,2).eq.constraints(i)%idx(j,1)).OR.&
 &(constraints(i)%idx(k,2).eq.constraints(i)%idx(j,2))) then
                      j = -1
                      exit
                    end if
                  end do
                  if (j.gt.0) then
                    kk = kk + 1
                    brks(kk) = k
                  end if
                end if
              end if
            end do
            if ((kk.gt.0).AND.(constraints(i)%nats.gt.10*tpn)) then
              nccgs = nccgs + 1
            end if
            if ((kk.gt.0).AND.(constraints(i)%nats.gt.10*tpn).AND.(iii.eq.2)) then
              tpx = 1
              ccg_limits(nccgs,1,1) = 1
              ccg_limits(nccgs,2,min(tpn,kk+1)) = constraints(i)%nr
              do k=1,kk
                if (tpx.eq.(kk+1)) exit
                if (brks(k).gt.(tpx*1.0*constraints(i)%nr/(1.0*min(tpn,kk+1)))) then
                  ccg_limits(nccgs,2,tpx) = brks(k)
                  ccg_limits(nccgs,1,tpx+1) = brks(k) + 1
                  tpx = tpx + 1
                else if ((min(tpn,kk+1)-tpx+1).gt.(kk-k+1)) then
                  ccg_limits(nccgs,2,tpx) = brks(k)
                  ccg_limits(nccgs,1,tpx+1) = brks(k) + 1
                  tpx = tpx + 1
                end if
              end do
              do k=1,tpn
                if (k.gt.(kk+1)) then
                  ccg_limits(nccgs,1,k) = 1
                  ccg_limits(nccgs,2,k) = 0
                end if
                ccg_limits(nccgs,3,k) = min(constraints(i)%nats,nint((k-1)*(1.0*constraints(i)%nats/(1.0*tpn)))) + 1
                ccg_limits(nccgs,4,k) = max(1,nint(k*(1.0*constraints(i)%nats/(1.0*tpn))))
                ccg_limits(nccgs,5,k) = i
                ccg_limits(nccgs,6,k) = ii
              end do
            end if
            deallocate(brks)
            deallocate(maxsn)
            deallocate(minsn)
          end if
        end if
      end do
      if ((iii.eq.1).AND.(nccgs.gt.0)) then
        allocate(ccg_limits(nccgs,6,tpn))
        ccg_limits(:,2,:) = 0
        ccg_limits(:,1,:) = 1
        ccg_limits(:,4,:) = 0
        ccg_limits(:,3,:) = 1
        nccgs = 0
      end if
    end do
  end if
  if (have_warned.EQV..true.) then
    if (nccgs.gt.0) then
      write(ilog,*) 'As a result, ',nccgs,' constraint groups were scheduled for internal parallelization.'
    else
      write(ilog,*) 'Despite unfavorable distribution of holonomic constraints across groups, group-internal &
 &parallelization was not deemed favorable/possible. This may indicate a bad ratio of the number of the threads &
 &and the size of the system.'
    end if
  end if 
!
  do iii=1,nccgs
    lastval(:) = 1
    do k=2,tpn
      do ii=1,2
        if (ccg_limits(iii,2*ii,k).ge.ccg_limits(iii,2*ii-1,k)) then
          if (ccg_limits(iii,2*ii-1,k).le.ccg_limits(iii,2*ii,lastval(ii))) then
            ccg_limits(iii,2*ii-1,k) = ccg_limits(iii,2*ii-1,k) + 1
          else
            lastval(ii) = k
          end if
        end if
      end do
    end do
  end do
!
!  do ii=1,nccgs
!    do k=1,tpn
!      write(*,*)' Thread ',k,'Mol ',ccg_limits(ii,5,k)
!      write(*,*) ccg_limits(ii,1:6,k)
!    end do
!  end do
!
end

!
!-------------------------------------------------------------------------------------------------------------
!
! a subroutine to deal with molecules that need internal parallelization for various operations
!
subroutine threads_setup_molecules()
!
  use threads
  use molecule
  use iounit
  use sequen
  use atoms
  use forces
  use system
  use zmatrix
  use polypep
  use shakeetal
  use clusters, ONLY: cdis_crit,cstorecalc
!
  implicit none
!
  integer tpn,k,kk,i,ii,iii,iiii,totc1,totc2,totc3,totc4,totc5,lastval(10),wkdon,blocksize,imol,curlevel,j,jj,ati
  integer maxblocksize,ixtl
  integer, ALLOCATABLE:: inlevel(:),inthread(:),threadcnts(:,:),ixtrace(:)
  logical hetero,placedit,have_warned
!
  tpn = thrdat%maxn ! -1
!
! this report is placed here so we can write it to ithread
 567 format(a,' bounds ',i4,': ',10000(i7,1x))
  if (nccgs.gt.0) then
    if (thrdat%verbosity.gt.2) then
      write(ithread,*) 'Bounds for individual constraint groups (SHAKE) by threads:'
      write(ithread,*) 
      do kk=1,nccgs
        write(ithread,*) 'Constraint group with ',constraints(ccg_limits(kk,5,1))%nats,' atoms:'
        do k=1,2
          write(ithread,567) 'Lower',k,(ccg_limits(kk,2*k-1,i),i=1,tpn)
          write(ithread,567) 'Upper',k,(ccg_limits(kk,2*k,i),i=1,tpn)
        end do
        write(ithread,*) 
      end do
    end if
  end if
!
  nmlgs = 0
!  thrdat%suppflag(1) = .true.
  if (thrdat%suppflag(1).EQV..true.) then
    allocate(thr_mlgix(nmol))
    thr_mlgix(:) = 0
    return
  end if
!
  totc1 = n
  have_warned = .false.
  do iii=1,2
    do ii=1,nmol
      totc2 = atmol(ii,2)-atmol(ii,1) + 1
      totc3 = rsmol(ii,2)-rsmol(ii,1) + 1
      if (allocated(dc_di).EQV..true.) then
        if (allocated(dc_di(ii)%valrecurs).EQV..true.) then
          totc4 = size(dc_di(ii)%valrecurs(1,:))
        else
          totc4 = 0
        end if
        totc5 = dc_di(ii)%maxntor
      else
        totc4 = 0
        totc5 = 0
      end if
      if ((((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))).OR.&
 &        ((1.0*totc2)/(1.0*(n)).gt.(0.1/tpn))) then
        if ((iii.eq.1).AND.(have_warned.EQV..false.).AND.((1.0*totc2)/(1.0*(n)).gt.(0.1/tpn))) then
          write(ilog,*) 'Warning. Distribution of atoms across molecules is such that a single molecule &
   &is too large to make simple parallelization by distributing molecules across threads efficient.'
          have_warned = .true.
        end if
        if (totc2.gt.10*tpn) then ! at leat 10 atoms per thread
          if (iii.eq.1) then
            nmlgs = nmlgs + 1
          else if (iii.eq.2) then
            nmlgs = nmlgs + 1
            do k=1,tpn
              mlg_limits(nmlgs,1,k) = min(totc2,nint((k-1)*(1.0*totc2/(1.0*tpn)))) + atmol(ii,1)
              mlg_limits(nmlgs,2,k) = max(1,nint(k*(1.0*totc2/(1.0*tpn)))) + atmol(ii,1) - 1
              mlg_limits(nmlgs,3,k) = min(totc3,nint((k-1)*(1.0*totc3/(1.0*tpn)))) + rsmol(ii,1)
              mlg_limits(nmlgs,4,k) = max(1,nint(k*(1.0*totc3/(1.0*tpn)))) + rsmol(ii,1) - 1
              mlg_limits(nmlgs,5:6,k) = ii
              if (totc4.gt.0) then
                mlg_limits(nmlgs,7,k) = min(totc4,nint((k-1)*(1.0*totc4/(1.0*tpn)))) + 1
                mlg_limits(nmlgs,8,k) = max(1,nint(k*(1.0*totc4/(1.0*tpn))))
              end if
              if (totc5.gt.0) then
                mlg_limits(nmlgs,9,k) = min(totc5,nint((k-1)*(1.0*totc5/(1.0*tpn)))) + 1
                mlg_limits(nmlgs,10,k) = max(1,nint(k*(1.0*totc5/(1.0*tpn))))
              end if
              mlg_limits(nmlgs,13,k) = min(totc5+3,nint((k-1)*(1.0*(totc5+3)/(1.0*tpn)))) - 2
              mlg_limits(nmlgs,14,k) = max(1,nint(k*(1.0*(totc5+3)/(1.0*tpn)))) - 3
            end do
          end if
        end if
      end if
    end do
    if ((iii.eq.1).AND.(nmlgs.gt.0)) then
      allocate(mlg_limits(nmlgs,14,tpn))
      mlg_limits(:,2,:) = 0
      mlg_limits(:,1,:) = 1
      mlg_limits(:,4,:) = 0
      mlg_limits(:,3,:) = 1
      mlg_limits(:,8,:) = 0
      mlg_limits(:,7,:) = 1
      mlg_limits(:,10,:) = 0
      mlg_limits(:,9,:) = 1
      mlg_limits(:,11,:) = 0
      mlg_limits(:,12,:) = 0
      mlg_limits(:,13,:) = 1
      mlg_limits(:,14,:) = 0
      nmlgs = 0
    end if
  end do
  if (have_warned.EQV..true.) then
    if (nmlgs.gt.0) then
      write(ilog,*) 'As a result, ',nmlgs,' molecules were scheduled for internal parallelization.'
    else
      write(ilog,*) 'Despite unfavorable distribution of atoms across molecules, molecule-internal &
 &parallelization was not deemed favorable. This indicates a bad ratio of the number of the threads &
 &and the size of the system.'
    end if
  end if
!
! try to reorder tasks for big molecules that allows internal parallelization of makexyz_formol
  maxblocksize = 1
  do iiii=1,2
    do iii=1,nmlgs
      imol = mlg_limits(iii,5,1)
      allocate(inlevel(n))
      allocate(inthread(n))
      allocate(threadcnts(n,tpn))
      allocate(ixtrace(atmol(imol,2) -atmol(imol,1) + 1))
!
!     first, generate a minimal trace from the last nonterminal atom in the last residue to the very first atom
      ixtl = 0
      ii = atmol(imol,2)
      do while (ii.gt.atmol(imol,1))
        if (ixtl.eq.0) then
          if (n12(ii).gt.1) then
            ixtl = ixtl + 1
            ixtrace(ixtl) = ii
          else
            ii = ii - 1
            cycle
          end if
        else
          ixtl = ixtl + 1
          ixtrace(ixtl) = ii
        end if
        ii = minval(i12(1:n12(ii),ii))
      end do
      ii = 1
      do while (ii.lt.ixtl)
        if (ixtrace(ii).le.(atmol(imol,1)+2)) exit
        do j=1,3
          placedit = .false.
          do k=ii+1,ixtl
            if ((iz(j,ixtrace(ii)).eq.ixtrace(k)).OR.(iz(j,ixtrace(ii)).le.0)) then
              placedit = .true.
              exit
            else if (abs(atmres(ixtrace(k))-atmres(ixtrace(ii))).gt.2) then
              exit
            end if
          end do
          if (placedit.EQV..false.) then
            do k=ixtl,ii,-1
              if (iz(j,ixtrace(ii)).lt.ixtrace(k)) then
                ixtrace(k+1) = iz(j,ixtrace(ii))
                exit
              end if
              ixtrace(k+1) = ixtrace(k)
            end do
            ixtl = ixtl + 1
          end if
        end do
        ii = ii + 1
      end do
!
!     smaller block sizes create better load balancing but require more barriers
!     set a reasonable blocksize based on molecule size (the block size does not change data dependency)
!     caveat: a strictly linear polymer (e.g., UA PEG) is never parallelizable according to the strategy below
      blocksize =  2*max(min(10,nint((1.0*(atmol(imol,2)-atmol(imol,1)+1))/(1.0*ixtl))),3)
      maxblocksize = max(maxblocksize,blocksize)
!
!      allocate memory
      inthread(:) = 0
      threadcnts(:,:) = 0
      inlevel(:) = atmol(imol,2) -atmol(imol,1) + 2
      inlevel(atmol(imol,1):(atmol(imol,1)+2)) = 1
      inthread(atmol(imol,1):(atmol(imol,1)+2)) = 1
      curlevel = 1
!     in the first pass, follow just the backbone, and just use the master thread
      tpn = 1
      do while (1.eq.1)
        wkdon = 0
        do ati=ixtl,1,-1
          ii = ixtrace(ati)
          if (ii.le.(atmol(imol,1)+2)) cycle
          if (inthread(ii).gt.0) cycle
          wkdon = wkdon + 1
          k = maxval(inlevel(iz(1:3,ii)))
!           write(*,*) 'wrkin ',ii,blocksize,k,iz(1:3,ii),inlevel(iz(1:3,ii))
          if (k.lt.curlevel) then
            jj = blocksize + 1
            placedit = .false.
            hetero = .false.
            jj = -1
            do j=1,3
              if (inlevel(iz(j,ii)).eq.k) then
                if (jj.eq.-1) then 
                  jj = inthread(iz(j,ii))
                else
                  if (jj.ne.inthread(iz(j,ii))) then
                    hetero = .true.
                    exit
                  end if
                end if
              end if
            end do
            if ((hetero.EQV..false.).AND.(threadcnts(k,jj).lt.blocksize)) then
              inthread(ii) = jj
              inlevel(ii) =  k
              threadcnts(k,jj) = threadcnts(k,jj) + 1
              if (iiii.eq.2) then
                mlg_ixs(1,k,iii,jj) = mlg_ixs(1,k,iii,jj) + 1
                mlg_ixs(mlg_ixs(1,k,iii,jj)+1,k,iii,jj) = ii
              end if
!              write(*,*) 'putting atom0 ',ii,' in ',jj,k
            else
              do kk=k+1,curlevel
                do j=1,tpn
                  if (threadcnts(kk,j).lt.blocksize) then
                    inthread(ii) = j
                    inlevel(ii) =  kk
                    threadcnts(kk,j) = threadcnts(kk,j) + 1
                    if (iiii.eq.2) then
                      mlg_ixs(1,kk,iii,j) = mlg_ixs(1,kk,iii,j) + 1
                      mlg_ixs(mlg_ixs(1,kk,iii,j)+1,kk,iii,j) = ii
                    end if
!                    write(*,*) 'putting atom1 ',ii,' in ',j,kk
                    placedit = .true.
                    exit
                  end if
                end do 
                if (placedit.EQV..true.) exit
              end do
            end if
          else if (k.eq.curlevel) then 
            hetero = .false.
            jj = -1
            do j=1,3
              if (inlevel(iz(j,ii)).eq.curlevel) then
                if (jj.eq.-1) then 
                  jj = inthread(iz(j,ii))
                else
                  if (jj.ne.inthread(iz(j,ii))) then
                    hetero = .true.
                    exit
                  end if
                end if
              end if
            end do
            if (hetero.EQV..true.) then
!             do nothing
            else if (threadcnts(curlevel,jj).lt.blocksize) then
              inthread(ii) = jj
              inlevel(ii) =  curlevel
              threadcnts(curlevel,jj) = threadcnts(curlevel,jj) + 1
              if (iiii.eq.2) then
                mlg_ixs(1,curlevel,iii,jj) = mlg_ixs(1,curlevel,iii,jj) + 1
                mlg_ixs(mlg_ixs(1,curlevel,iii,jj)+1,curlevel,iii,jj) = ii
              end if
!                 write(*,*) 'putting atom2 ',ii,' in ',jj,curlevel
            end if
          else
          end if
        end do
        curlevel = curlevel + 1
!        write(*,*) wkdon,curlevel
!       if the second condition triggers, we have a problem as it means that the backbone Z-matrix
!       relied on a terminal atom somewhere
        if ((wkdon.eq.0).OR.(curlevel.gt.(atmol(imol,2)-atmol(imol,1)+1))) exit
      end do
      do ii=atmol(imol,1)+3,atmol(imol,2)
        if (inlevel(ii).lt.(atmol(imol,2) -atmol(imol,1) + 2)) then
          mlg_limits(iii,11,1) = max(mlg_limits(iii,11,1),inlevel(ii))
        end if
      end do
      tpn = thrdat%maxn
!     in the second pass, consider all atoms with at least two bonds (this will mostly delay hydrogens)
      do while (1.eq.1)
        wkdon = 0
        do ii=atmol(imol,1)+3,atmol(imol,2)
          if (inthread(ii).gt.0) cycle
          if (n12(ii).le.1) cycle
          wkdon = wkdon + 1
          k = maxval(inlevel(iz(1:3,ii)))
          if (k.lt.curlevel) then
            jj = blocksize + 1
            placedit = .false.
            hetero = .false.
            jj = -1
            do j=1,3
              if (inlevel(iz(j,ii)).eq.k) then
                if (jj.eq.-1) then 
                  jj = inthread(iz(j,ii))
                else
                  if (jj.ne.inthread(iz(j,ii))) then
                    hetero = .true.
                    exit
                  end if
                end if
              end if
            end do
            if ((hetero.EQV..false.).AND.(threadcnts(k,jj).lt.blocksize)) then
              inthread(ii) = jj
              inlevel(ii) =  k
              threadcnts(k,jj) = threadcnts(k,jj) + 1
              if (iiii.eq.2) then
                mlg_ixs(1,k,iii,jj) = mlg_ixs(1,k,iii,jj) + 1
                mlg_ixs(mlg_ixs(1,k,iii,jj)+1,k,iii,jj) = ii
              end if
!              write(*,*) 'putting atom0 ',ii,' in ',jj,k
            else
              do kk=k+1,curlevel
                do j=1,tpn
                  if (threadcnts(kk,j).lt.blocksize) then
                    inthread(ii) = j
                    inlevel(ii) =  kk
                    threadcnts(kk,j) = threadcnts(kk,j) + 1
                    if (iiii.eq.2) then
                      mlg_ixs(1,kk,iii,j) = mlg_ixs(1,kk,iii,j) + 1
                      mlg_ixs(mlg_ixs(1,kk,iii,j)+1,kk,iii,j) = ii
                    end if
!                    write(*,*) 'putting atom1 ',ii,' in ',j,kk
                    placedit = .true.
                    exit
                  end if
                end do 
                if (placedit.EQV..true.) exit
              end do
            end if
          else if (k.eq.curlevel) then 
            hetero = .false.
            jj = -1
            do j=1,3
              if (inlevel(iz(j,ii)).eq.curlevel) then
                if (jj.eq.-1) then 
                  jj = inthread(iz(j,ii))
                else
                  if (jj.ne.inthread(iz(j,ii))) then
                    hetero = .true.
                    exit
                  end if
                end if
              end if
            end do
            if (hetero.EQV..true.) then
            else if (threadcnts(curlevel,jj).lt.blocksize) then
              inthread(ii) = jj
              inlevel(ii) =  curlevel
              threadcnts(curlevel,jj) = threadcnts(curlevel,jj) + 1
              if (iiii.eq.2) then
                mlg_ixs(1,curlevel,iii,jj) = mlg_ixs(1,curlevel,iii,jj) + 1
                mlg_ixs(mlg_ixs(1,curlevel,iii,jj)+1,curlevel,iii,jj) = ii
              end if
!                 write(*,*) 'putting atom2 ',ii,' in ',jj,curlevel
            end if
          else
          end if
        end do
        curlevel = curlevel + 1
        if ((wkdon.eq.0).OR.(curlevel.gt.(atmol(imol,2)-atmol(imol,1)+1))) exit
      end do
      do ii=atmol(imol,1)+3,atmol(imol,2)
        if (inlevel(ii).lt.(atmol(imol,2) -atmol(imol,1) + 2)) then
          mlg_limits(iii,11,1) = max(mlg_limits(iii,11,1),inlevel(ii))
        end if
      end do
      do while (1.eq.1)
        wkdon = 0
        do ii=atmol(imol,1)+3,atmol(imol,2)
          if (inthread(ii).gt.0) cycle
          wkdon = wkdon + 1
          k = maxval(inlevel(iz(1:3,ii)))
          if (k.lt.curlevel) then
            jj = blocksize + 1
            placedit = .false.
            hetero = .false.
            jj = -1
            do j=1,3
              if (inlevel(iz(j,ii)).eq.k) then
                if (jj.eq.-1) then 
                  jj = inthread(iz(j,ii))
                else
                  if (jj.ne.inthread(iz(j,ii))) then
                    hetero = .true.
                    exit
                  end if
                end if
              end if
            end do
            if ((hetero.EQV..false.).AND.(threadcnts(k,jj).lt.blocksize)) then
              inthread(ii) = jj
              inlevel(ii) =  k
              threadcnts(k,jj) = threadcnts(k,jj) + 1
              if (iiii.eq.2) then
                mlg_ixs(1,k,iii,jj) = mlg_ixs(1,k,iii,jj) + 1
                mlg_ixs(mlg_ixs(1,k,iii,jj)+1,k,iii,jj) = ii
              end if
!              write(*,*) 'putting atom0 ',ii,' in ',jj,k
            else
              do kk=k+1,curlevel
                do j=1,tpn
                  if (threadcnts(kk,j).lt.blocksize) then
                    inthread(ii) = j
                    inlevel(ii) =  kk
                    threadcnts(kk,j) = threadcnts(kk,j) + 1
                      if (iiii.eq.2) then
                        mlg_ixs(1,kk,iii,j) = mlg_ixs(1,kk,iii,j) + 1
                        mlg_ixs(mlg_ixs(1,kk,iii,j)+1,kk,iii,j) = ii
                      end if
!                    write(*,*) 'putting atom1 ',ii,' in ',j,kk
                    placedit = .true.
                    exit
                  end if
                end do 
                if (placedit.EQV..true.) exit
              end do
            end if
          else if (k.eq.curlevel) then 
            hetero = .false.
            jj = -1
            do j=1,3
              if (inlevel(iz(j,ii)).eq.curlevel) then
                if (jj.eq.-1) then 
                  jj = inthread(iz(j,ii))
                else
                  if (jj.ne.inthread(iz(j,ii))) then
                    hetero = .true.
                    exit
                  end if
                end if
              end if
            end do
            if (hetero.EQV..true.) then
            else if (threadcnts(curlevel,jj).lt.blocksize) then
              inthread(ii) = jj
              inlevel(ii) =  curlevel
              threadcnts(curlevel,jj) = threadcnts(curlevel,jj) + 1
              if (iiii.eq.2) then
                mlg_ixs(1,curlevel,iii,jj) = mlg_ixs(1,curlevel,iii,jj) + 1
                mlg_ixs(mlg_ixs(1,curlevel,iii,jj)+1,curlevel,iii,jj) = ii
              end if
!                 write(*,*) 'putting atom2 ',ii,' in ',jj,curlevel
            end if
          else
          end if
        end do
        curlevel = curlevel + 1
        if (wkdon.eq.0) exit
      end do
      do ii=atmol(imol,1)+3,atmol(imol,2)
        if (inlevel(ii).lt.(atmol(imol,2) -atmol(imol,1) + 2)) then
          mlg_limits(iii,11,1) = max(mlg_limits(iii,11,1),inlevel(ii))
        end if
      end do
      deallocate(ixtrace)
      deallocate(inlevel)
      deallocate(inthread)
      deallocate(threadcnts)
      if (iiii.eq.2) then
        mlg_limits(iii,11,2:thrdat%maxn) = mlg_limits(iii,11,1)
      end if 
    end do
    if ((iiii.eq.1).AND.(nmlgs.gt.0)) then
      allocate(mlg_ixs(maxblocksize+1,maxval(mlg_limits(1:nmlgs,11,1)),nmlgs,thrdat%maxn))
      mlg_ixs(1,:,:,:) = 0 
    end if
  end do
!
  allocate(thr_mlgix(nmol))
  thr_mlgix(:) = 0
  do iii=1,nmlgs
    thr_mlgix(mlg_limits(iii,5,1)) = iii
    lastval(:) = 1
    do k=2,tpn
      do ii=1,5
        if (ii.eq.3) cycle
        if (mlg_limits(iii,2*ii,k).ge.mlg_limits(iii,2*ii-1,k)) then
          if (mlg_limits(iii,2*ii-1,k).le.mlg_limits(iii,2*ii,lastval(ii))) then
            mlg_limits(iii,2*ii-1,k) = mlg_limits(iii,2*ii-1,k) + 1
          else
            lastval(ii) = k
          end if
        end if
      end do
    end do
  end do
  do iii=1,nmlgs
    mlg_limits(iii,12,1) = thrdat%maxn
    do i=1,thrdat%maxn
      if (maxval(mlg_ixs(1,1:mlg_limits(iii,11,1),iii,i)).le.0) then
        mlg_limits(iii,12,1) = i-1
        exit
      end if
    end do
  end do
!
  if ((nmlgs.gt.0).AND.(((use_dyn.EQV..true.).AND.(fycxyz.eq.1)).OR.&
 &                      ((cstorecalc.le.nsim).AND.((cdis_crit.eq.2).OR.(cdis_crit.eq.4))))) then
    allocate(thr_cart2f_hlp(22,maxval(mlg_limits(:,8,:)),tpn))
  end if
!
  if (nmlgs.gt.0) then
    if (thrdat%verbosity.gt.2) then
      write(ithread,*) 'Bounds for individual molecules by threads:'
      write(ithread,*) 
      do kk=1,nmlgs
        write(ithread,*) 'Molecule ',mlg_limits(kk,5,1),' with ',&
 &             atmol(mlg_limits(kk,5,1),2)-atmol(mlg_limits(kk,5,1),1)+1,' atoms:'
        k = 0
        do j=1,7
!          if ((j.eq.3).OR.(j.eq.6)) cycle
          k = k + 1
          write(ithread,567) 'Lower',k,(mlg_limits(kk,2*j-1,i),i=1,tpn)
          write(ithread,567) 'Upper',k,(mlg_limits(kk,2*j,i),i=1,tpn)
        end do
        write(ithread,*) 
      end do
    end if
  end if
!
end

!
!--------------------------------------------------------------------------------------------
!
subroutine threads_sanity_checks()
!
  use threads
  use energies
  use system
  use cutoffs
  use iounit
  use movesets
  use sequen
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
#ifdef ENABLE_MPI
!
  integer j
#endif
!
! some compatibility sanity checks
  if (use_CORR.EQV..true.) then
    write(ilog,*) 'Fatal. The use of the (quasi-obsolete) correction potential (FMCSC_SC_EXTRA) is not &
 &supported when using multi-threaded code.'
    call fexit()
  end if
!
  if ((use_DSSP.EQV..true.).OR.(use_ZSEC.EQV..true.)) then
    write(ilog,*) 'Warning. When using multi-threaded code, the computation of secondary structure restraint &
 &potentials (FMCSC_SC_ZSEC and/or FMCSC_SC_DSSP) is not (yet) parallelized in any way.'
  end if
!
  if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
    if ((ens%flag.eq.6).AND.(particleflucfreq.gt.0.0)) then
      write(ilog,*) 'Fatal. When using multi-threaded code, the semigrand ensemble is currently not supported. &
 &Please check back later.'
      call fexit()
    end if
    if (use_cutoffs.EQV..false.) then
      write(ilog,*) 'Fatal. The combination of explicit disabling of cutoffs (FMCSC_CUTOFFMODE 1) and MC-based sampling methods &
 &is not yet supported when using multi-threaded code. Please manually increase cutoffs instead so that they exceed &
 &all possible system dimensions, which achieves the same effect (keywords FMCSC_NBCUTOFF and FMCSC_ELCUTOFF).'
      call fexit()
    end if
    if ((have_djcr.EQV..true.).OR.(have_sjcr.EQV..true.).OR.(have_docr.EQV..true.).OR.(have_nuccr.EQV..true.)) then
      write(ilog,*) 'Fatal. Most types of concerted rotation moves are not yet supported when using multi-threaded code. &
 &Please check back later.'
      call fexit()
    end if
    if (have_clurb.EQV..true.) then
      write(ilog,*) 'Fatal. Cluster rigid body moves are not yet supported when using multi-threaded code. &
 &Please check back later.'
      call fexit()
    end if
  end if  
!
  if (((dyn_mode.eq.3).OR.(dyn_mode.eq.7)).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Fatal. The Langevin integrator in internal coordinate space is currently not supported when &
 &using multi-threaded code.'
    call fexit()
  end if 
!
#ifdef ENABLE_MPI
!
  if ((use_REMC.EQV..true.).AND.(re_freq.le.nsim)) then
    do j=1,re_conddim
      if ((re_types(j).eq.34).OR.((re_types(j).ge.13).AND.(re_types(j).le.16))) then
        write(ilog,*) 'Warning. One or more of the dimensions selected for replica exchange require setup computation that &
 &are not thread-parallelized in any way when using hybrid OpenMP/MPI code.'
        exit
      end if
    end do
  end if
#endif
!
  if ((thrdat%test.EQV..true.).AND.(do_restart.EQV..true.)) then
    write(ilog,*) 'Warning. Testing of multi-threaded routines (FMCSC_THREADS_TEST) is not possible in restarted runs. &
 &Request ignored.'
    thrdat%test = .false.
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
! a routine to emulate a REDUCE clause (sum) on a thread-private vector vec written to thr_rhlper(ix1:(ix1+vecl-1),ix2)
!
subroutine thr_combinevec(vec,vecl,ix1,ix2,k,tpi)
!
  use threads, ONLY: thr_limits,thrdat,thr_rhlper
!
  implicit none
!
  integer, INTENT(IN):: vecl,tpi,k,ix1,ix2
  RTYPE, INTENT(IN):: vec(vecl)
!
  integer i,j,OMP_GET_NUM_THREADS,tpn,ixl,ixn
!
  tpn = OMP_GET_NUM_THREADS()
  ixn = ix1 + vecl - 1
!$OMP BARRIER
  if ((1.0*vecl)/(1.0*tpn).le.1000.) then
!$OMP CRITICAL(VEC_UP)
    thr_rhlper(ix1:ixn,ix2) = thr_rhlper(ix1:ixn,ix2) + vec(1:vecl)
!$OMP END CRITICAL(VEC_UP)
!$OMP BARRIER
  else
    if (thr_limits(k+1,tpi).le.0) then ! all intervals must be valid due to condition above -> if so, don't recheck
      i = 1
      j = vecl
      call threads_bounds(i,j,tpi,tpn,thr_limits(k:(k+1),tpi))
    end if
!$OMP BARRIER
    do i=0,thrdat%maxn-1
      j = mod(i+tpi,thrdat%maxn) + 1
      ixl = thr_limits(k,j) + ix1 - 1
      ixn = thr_limits(k+1,j) + ix1 - 1
      thr_rhlper(ixl:ixn,ix2) = thr_rhlper(ixl:ixn,ix2) + vec(thr_limits(k,j):thr_limits(k+1,j))
!$OMP BARRIER
    end do
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
! the same as thr_combinevec for an integer vector (written to thr_hlper(ix1:(ix1+vecl-1),ix2))
!
subroutine thr_combineivec(vec,vecl,ix1,ix2,k,tpi)
!
  use threads, ONLY: thr_limits,thrdat,thr_hlper
!
  implicit none
!
  integer, INTENT(IN):: vecl,tpi,k,ix1,ix2,vec(vecl)
!
  integer i,j,OMP_GET_NUM_THREADS,tpn,ixl,ixn
!
  tpn = OMP_GET_NUM_THREADS()
  ixn = ix1 + vecl - 1
!$OMP BARRIER
  if ((1.0*vecl)/(1.0*tpn).le.1000.) then
!$OMP CRITICAL(VEC_UP)
    thr_hlper(ix1:ixn,ix2) = thr_hlper(ix1:ixn,ix2) + vec(1:vecl)
!$OMP END CRITICAL(VEC_UP)
!$OMP BARRIER
  else
    if (thr_limits(k+1,tpi).le.0) then ! all intervals must be valid due to condition above -> if so, don't recheck
      i = 1
      j = vecl
      call threads_bounds(i,j,tpi,tpn,thr_limits(k:(k+1),tpi))
    end if
!$OMP BARRIER
    do i=0,thrdat%maxn-1
      j = mod(i+tpi,thrdat%maxn) + 1
      ixl = thr_limits(k,j) + ix1 - 1
      ixn = thr_limits(k+1,j) + ix1 - 1
      thr_hlper(ixl:ixn,ix2) = thr_hlper(ixl:ixn,ix2) + vec(thr_limits(k,j):thr_limits(k+1,j))
!$OMP BARRIER
    end do
  end if
!
end
!
!------------------------------------------------------------------------------------------------------
!
! a routine to perform a reduce operation on a thread-private FP variable written to thr_rhlper(ix1,ix2)
!
subroutine thr_reduction(mode,val,ix1,ix2)
!
  use threads, ONLY: thr_rhlper
!
  implicit none
!
  integer, INTENT(IN):: mode,ix1,ix2
  RTYPE, INTENT(IN):: val
!
!$OMP BARRIER
  if (mode.eq.-1) then ! minimum
!$OMP CRITICAL(SCALAR_UP)
    thr_rhlper(ix1,ix2) = min(thr_rhlper(ix1,ix2),val)
!$OMP END CRITICAL(SCALAR_UP)
  else if (mode.eq.0) then ! sum
!$OMP CRITICAL(SCALAR_UP)
    thr_rhlper(ix1,ix2) = thr_rhlper(ix1,ix2) + val
!$OMP END CRITICAL(SCALAR_UP)
  else if (mode.eq.1) then ! maximum
!$OMP CRITICAL(SCALAR_UP)
    thr_rhlper(ix1,ix2) = max(thr_rhlper(ix1,ix2),val)
!$OMP END CRITICAL(SCALAR_UP)
  end if
!$OMP BARRIER
!
end
!
!------------------------------------------------------------------------------------------------------
!
! a routine to perform a reduce operation on a thread-private FP vector buffered on thr_rutil or thr_rutil2
! and copied back into possibly thread-private FP vector
!
subroutine thr_freduction(mode,val,valsz,ix1,tpi)
!
  use threads, ONLY: thr_rutil,thr_rutil2,thrdat,thr_limits
!
  implicit none
!
  integer, INTENT(IN):: mode,ix1,valsz,tpi
  RTYPE, INTENT(INOUT):: val(valsz)
!
  integer i,j,ixn,tpn,OMP_GET_NUM_THREADS,which
!
  tpn = OMP_GET_NUM_THREADS()
  ixn = ix1 + valsz - 1
  if (size(thr_rutil,dim=1).ge.(ix1+valsz-1)) then
    thr_rutil(ix1:ixn,tpi) = val(1:valsz)
    which = 1
  else if (allocated(thr_rutil2).EQV..true.) then
    if ((size(thr_rutil2,dim=1).ge.ixn).AND.(size(thr_rutil2,dim=2).ge.tpn)) then
      thr_rutil2(ix1:ixn,tpi) = val(1:valsz)
      which = 2
    else
!     be radical - this should not happen often
!$OMP BARRIER
      if (tpi.le.1) then
        deallocate(thr_rutil2)
        allocate(thr_rutil2(ixn,tpn))
      end if
!$OMP BARRIER
      thr_rutil2(ix1:ixn,tpi) = val(1:valsz)
      which = 2
    end if
  else
!$OMP BARRIER
    if (tpi.le.1) then
      allocate(thr_rutil2(ixn,tpn))
    end if
!$OMP BARRIER
    thr_rutil2(ix1:ixn,tpi) = val(1:valsz)
    which = 2
  end if
  i = 1
  j = valsz
  call threads_bounds(i,j,tpi,tpn,thr_limits(67:68,tpi))
  thr_limits(67:68,tpi) = thr_limits(67:68,tpi) + ix1 - 1
!$OMP BARRIER
  if (which.eq.1) then
    if (mode.eq.-1) then ! minimum
      do i=thr_limits(67,tpi),thr_limits(68,tpi)
        thr_rutil(i,1) = minval(thr_rutil(i,1:tpn))
      end do
    else if (mode.eq.0) then ! sum
      do i=thr_limits(67,tpi),thr_limits(68,tpi)
        thr_rutil(i,1) = sum(thr_rutil(i,1:tpn))
      end do
    else if (mode.eq.1) then ! maximum
      do i=thr_limits(67,tpi),thr_limits(68,tpi)
        thr_rutil(i,1) = maxval(thr_rutil(i,1:tpn))
      end do
    end if
  else if (which.eq.2) then
    if (mode.eq.-1) then ! minimum
      do i=thr_limits(67,tpi),thr_limits(68,tpi)
        thr_rutil2(i,1) = minval(thr_rutil2(i,1:tpn))
      end do
    else if (mode.eq.0) then ! sum
      do i=thr_limits(67,tpi),thr_limits(68,tpi)
        thr_rutil2(i,1) = sum(thr_rutil2(i,1:tpn))
      end do
    else if (mode.eq.1) then ! maximum
      do i=thr_limits(67,tpi),thr_limits(68,tpi)
        thr_rutil2(i,1) = maxval(thr_rutil2(i,1:tpn))
      end do
    end if
  end if
  thr_limits(67,tpi) = 1
  thr_limits(68,tpi) = 0
!$OMP BARRIER
  if (which.eq.1) then
    val(1:valsz) = thr_rutil(ix1:ixn,1)
  else if (which.eq.2) then
    val(1:valsz) = thr_rutil2(ix1:ixn,1)
  end if
!
end

!------------------------------------------------------------------------------------------------------
!
! a routine to perform a reduce operation on a thread-private FP vector buffered on thr_iutil2
! and copied back into possibly thread-private FP vector
!
subroutine thr_ireduction(mode,val,valsz,ix1,tpi)
!
  use threads, ONLY: thr_iutil2,thrdat,thr_limits
!
  implicit none
!
  integer, INTENT(IN):: mode,ix1,valsz,tpi
  integer, INTENT(INOUT):: val(valsz)
!
  integer i,j,ixn,tpn,OMP_GET_NUM_THREADS
!
  tpn = OMP_GET_NUM_THREADS()
  ixn = ix1 + valsz - 1
  if (allocated(thr_iutil2).EQV..true.) then
    if ((size(thr_iutil2,dim=1).ge.ixn).AND.(size(thr_iutil2,dim=2).ge.tpn)) then
      thr_iutil2(ix1:ixn,tpi) = val(1:valsz)
    else
!     be radical - this should not happen often
!$OMP BARRIER
      if (tpi.le.1) then
        deallocate(thr_iutil2)
        allocate(thr_iutil2(ixn,tpn))
      end if
!$OMP BARRIER
      thr_iutil2(ix1:ixn,tpi) = val(1:valsz)
    end if
  else
!$OMP BARRIER
    if (tpi.le.1) then
      allocate(thr_iutil2(ixn,tpn))
    end if
!$OMP BARRIER
    thr_iutil2(ix1:ixn,tpi) = val(1:valsz)
  end if

  i = 1
  j = valsz
  call threads_bounds(i,j,tpi,tpn,thr_limits(67:68,tpi))
  thr_limits(67:68,tpi) = thr_limits(67:68,tpi) + ix1 - 1
!$OMP BARRIER
  if (mode.eq.-1) then ! minimum
    do i=thr_limits(67,tpi),thr_limits(68,tpi)
      thr_iutil2(i,1) = minval(thr_iutil2(i,1:tpn))
    end do
  else if (mode.eq.0) then ! sum
    do i=thr_limits(67,tpi),thr_limits(68,tpi)
      thr_iutil2(i,1) = sum(thr_iutil2(i,1:tpn))
    end do
  else if (mode.eq.1) then ! maximum
    do i=thr_limits(67,tpi),thr_limits(68,tpi)
      thr_iutil2(i,1) = maxval(thr_iutil2(i,1:tpn))
    end do
  end if
  thr_limits(67,tpi) = 1
  thr_limits(68,tpi) = 0
!$OMP BARRIER
  val(1:valsz) = thr_iutil2(ix1:ixn,1)
!
end
!
!----------------------------------------------------------------------------------------------
!
! job_flags indicates the tasks to be performed:
!          (1): compute comm and rgv
!          (2): populate Z-matrix
!          (3): update center of mass velocity
!          (4): shift periodic image
!
subroutine molops_threads(idx,tpi,job_flags)
!
  use iounit
  use molecule
  use atoms
  use forces
  use threads
  use zmatrix
  use system
!
  implicit none
!
  integer, INTENT(IN):: tpi,idx
  logical, INTENT(IN):: job_flags(4)
!
  integer i,imol,i1,i2,i3,i4,j
  RTYPE tvec(3),getbang,getztor
  RTYPE dum3(3),hlp1,hlp2,hlp3,hlp4,hlp5,hlp6,immol
!
  imol = mlg_limits(idx,5,tpi)
!
  immol = 1.0/molmass(moltypid(imol))
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
!$OMP SINGLE
    if (job_flags(1).EQV..true.) then
      comm(imol,1) = x(atmol(imol,1))
      comm(imol,2) = y(atmol(imol,1))
      comm(imol,3) = z(atmol(imol,1))
    end if
    if (job_flags(2).EQV..true.) then
      blen(atmol(imol,1)) = 0.0
      bang(atmol(imol,1)) = 0.0
      ztor(atmol(imol,1)) = 0.0
    end if
    if (job_flags(3).EQV..true.) then
      dc_di(imol)%v(1) = cart_v(atmol(imol,1),1)
      dc_di(imol)%v(2) = cart_v(atmol(imol,1),2)
      dc_di(imol)%v(3) = cart_v(atmol(imol,1),3)
    end if
    if ((job_flags(4).EQV..true.).AND.(bnd_type.eq.1)) then
      if (bnd_shape.eq.1) then
        tvec(:) = -bnd_params(1:3)*anint(comm(imol,:)/bnd_params(1:3))
        x(atmol(imol,1)) = x(atmol(imol,1)) + tvec(1)
        y(atmol(imol,1)) = y(atmol(imol,1)) + tvec(2)
        z(atmol(imol,1)) = z(atmol(imol,1)) + tvec(3)
      else if (bnd_shape.eq.3) then
        tvec(1:2) = 0.0
        tvec(3) = -bnd_params(6)*anint(comm(imol,3)/bnd_params(6))
        z(atmol(imol,1)) = z(atmol(imol,1)) + tvec(3)
      end if
      comm(imol,:) = comm(imol,:) + tvec(:)
    end if
!$OMP END SINGLE
    return
  end if
!
  if (job_flags(2).EQV..true.) then
    do j=mlg_limits(idx,1,tpi),mlg_limits(idx,2,tpi)
      i1 = iz(1,j)
      i2 = iz(2,j)
      i3 = iz(3,j)
      i4 = iz(4,j)
      if (i1.gt.0) then
        blen(j) = sqrt((x(i1)-x(j))**2 + (y(i1)-y(j))**2 + (z(i1)-z(j))**2)
        if (i2.gt.0) then
          bang(j) = getbang(j,i1,i2)
          if (i3.gt.0) then
            if (i4.eq.0) then
              ztor(j) = getztor(j,i1,i2,i3)
            else
!             this case should be avoided
              ztor(j) = getbang(j,i1,i3)
              hlp1 = getztor(j,i1,i2,i3)
              if (hlp1.gt.0.0) then
                iz(4,j) = -1
              else
                iz(4,j) = 1
              end if
            end if
          else
            ztor(j) = 0.0
          end if
        else
          bang(j) = 0.0
          ztor(j) = 0.0
        end if
      else
        blen(j) = 0.0
        bang(j) = 0.0
        ztor(j) = 0.0
      end if
    end do
  end if
!
  if ((job_flags(1).EQV..true.).OR.(job_flags(3).EQV..true.)) then
!
!$OMP SINGLE
    if (job_flags(1).EQV..true.) then
      rgpcsm(imol,:,:) = 0.0
      rgpcsm(imol,1,1) = 1.0
      rgpcsm(imol,2,2) = 1.0
      rgpcsm(imol,3,3) = 1.0
      comm(imol,:) = 0.0
      rgvm(imol) = 0.0
    end if
    if (job_flags(3).EQV..true.) dc_di(imol)%v(1:3) = 0.0
!$OMP END SINGLE NOWAIT
!
    if (job_flags(1).EQV..true.) then
      thr_rutil(1,tpi)= sum(mass(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi))*x(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)))
      thr_rutil(2,tpi)= sum(mass(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi))*y(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)))
      thr_rutil(3,tpi)= sum(mass(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi))*z(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)))
    end if
    if (job_flags(3).EQV..true.) then
      thr_rutil(4,tpi)=sum(mass(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi))*cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),1))
      thr_rutil(5,tpi)=sum(mass(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi))*cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),2))
      thr_rutil(6,tpi)=sum(mass(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi))*cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),3))
    end if
!
!$OMP BARRIER
!$OMP SINGLE
    if (job_flags(1).EQV..true.) comm(imol,1:3) = comm(imol,1:3) + immol*sum(thr_rutil(1:3,1:thrdat%maxn),dim=2)
    if (job_flags(3).EQV..true.) dc_di(imol)%v(1:3) = dc_di(imol)%v(1:3) + immol*sum(thr_rutil(4:6,1:thrdat%maxn),dim=2)
!$OMP END SINGLE
!
    if (job_flags(1).EQV..true.) then
      dum3(:) = comm(imol,:)
      thr_rutil(7:9,tpi) = 0.0
!
!   update inertia tensor
      do i=mlg_limits(idx,1,tpi),mlg_limits(idx,2,tpi)
        hlp1 = x(i)-dum3(1)
        hlp2 = y(i)-dum3(2)
        hlp3 = z(i)-dum3(3)
        hlp4 = hlp1*hlp1
        hlp5 = hlp2*hlp2
        hlp6 = hlp3*hlp3
        thr_rutil(7,tpi) = thr_rutil(7,tpi) + mass(i)*(hlp5+hlp6)
        thr_rutil(8,tpi) = thr_rutil(8,tpi) + mass(i)*(hlp4+hlp6)
        thr_rutil(9,tpi) = thr_rutil(9,tpi) + mass(i)*(hlp4+hlp5)
      end do
    end if
  end if
!
  if (job_flags(4).EQV..true.) then
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        tvec(:) = -bnd_params(1:3)*anint((comm(imol,:)-(bnd_params(4:6)+0.5*bnd_params(1:3)))/bnd_params(1:3))
        x(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) = x(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) + tvec(1)
        y(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) = y(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) + tvec(2)
        z(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) = z(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) + tvec(3)
      else if (bnd_shape.eq.3) then
        tvec(1:2) = 0.0
        tvec(3) = -bnd_params(6)*anint((comm(imol,3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6))
        z(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) = z(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi)) + tvec(3)
      end if
    end if
  end if
!
  if ((job_flags(1).EQV..true.).OR.((job_flags(4).EQV..true.).AND.(bnd_type.eq.1))) then
!$OMP BARRIER
!$OMP SINGLE
    rgvm(imol) = rgvm(imol) + sum(thr_rutil(7,1:thrdat%maxn)+thr_rutil(8,1:thrdat%maxn)+thr_rutil(9,1:thrdat%maxn))
    if (job_flags(1).EQV..true.) rgvm(imol) = sqrt(rgvm(imol))
    if ((job_flags(4).EQV..true.).AND.(bnd_type.eq.1)) then
      comm(imol,:) = comm(imol,:) + tvec(:)
      com(imol,:) = com(imol,:) + tvec(:)
    end if
!$OMP END SINGLE NOWAIT
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
! a similar routine operating using no mass-weighting
! job_flags indicates the tasks to be performed:
!          (1): compute com and gyration tensor
!          (2): populate Z-matrix
!          (3): diagonalize gyration tensor and populate rgevs (SINGLE only)
!          (4): shift periodic image
!
subroutine molops_threads_geo(idx,tpi2,job_flags)
!
  use iounit
  use molecule
  use atoms
  use forces
  use threads
  use zmatrix
  use system
!
  implicit none
!
  integer, INTENT(IN):: tpi2,idx
  logical, INTENT(IN):: job_flags(4)
!
  integer i,imol,i1,i2,i3,i4,j
  RTYPE tvec(3),getbang,getztor,rgten(3,3),evmat(3,3)
  RTYPE dum1(3),dum3(3),hlp1,hlp2,hlp3,hlp4,hlp5,hlp6,immol
!
  imol = mlg_limits(idx,5,tpi2)
!
  if (job_flags(1).EQV..true.) immol = 1.0/(dble(atmol(imol,2)-atmol(imol,1)+1))
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
!$OMP SINGLE
    if (job_flags(1).EQV..true.) then
      com(imol,1) = x(atmol(imol,1))
      com(imol,2) = y(atmol(imol,1))
      com(imol,3) = z(atmol(imol,1))
    end if
    if (job_flags(2).EQV..true.) then
      blen(atmol(imol,1)) = 0.0
      bang(atmol(imol,1)) = 0.0
      ztor(atmol(imol,1)) = 0.0
    end if
    if ((job_flags(4).EQV..true.).AND.(bnd_type.eq.1)) then
      if (bnd_shape.eq.1) then
        tvec(:) = -bnd_params(1:3)*anint(com(imol,:)/bnd_params(1:3))
        x(atmol(imol,1)) = x(atmol(imol,1)) + tvec(1)
        y(atmol(imol,1)) = y(atmol(imol,1)) + tvec(2)
        z(atmol(imol,1)) = z(atmol(imol,1)) + tvec(3)
      else if (bnd_shape.eq.3) then
        tvec(1:2) = 0.0
        tvec(3) = -bnd_params(6)*anint(com(imol,3)/bnd_params(6))
        z(atmol(imol,1)) = z(atmol(imol,1)) + tvec(3)
      end if
      com(imol,:) = com(imol,:) + tvec(:)
    end if
!$OMP END SINGLE
    return
  end if
!
  if (job_flags(2).EQV..true.) then
    do j=mlg_limits(idx,1,tpi2),mlg_limits(idx,2,tpi2)
      i1 = iz(1,j)
      i2 = iz(2,j)
      i3 = iz(3,j)
      i4 = iz(4,j)
      if (i1.gt.0) then
        blen(j) = sqrt((x(i1)-x(j))**2 + (y(i1)-y(j))**2 + (z(i1)-z(j))**2)
        if (i2.gt.0) then
          bang(j) = getbang(j,i1,i2)
          if (i3.gt.0) then
            if (i4.eq.0) then
              ztor(j) = getztor(j,i1,i2,i3)
            else
!             this case should be avoided
              ztor(j) = getbang(j,i1,i3)
              hlp1 = getztor(j,i1,i2,i3)
              if (hlp1.gt.0.0) then
                iz(4,j) = -1
              else
                iz(4,j) = 1
              end if
            end if
          else
            ztor(j) = 0.0
          end if
        else
          bang(j) = 0.0
          ztor(j) = 0.0
        end if
      else
        blen(j) = 0.0
        bang(j) = 0.0
        ztor(j) = 0.0
      end if
    end do
  end if
!
  if (job_flags(1).EQV..true.) then
!
    thr_rutil(1,tpi2) = sum(x(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)))
    thr_rutil(2,tpi2) = sum(y(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)))
    thr_rutil(3,tpi2) = sum(z(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)))
!
!$OMP BARRIER
!$OMP SINGLE
    com(imol,1:3) = immol*sum(thr_rutil(1:3,1:thrdat%maxn),dim=2)
!$OMP END SINGLE
!
    dum3(:) = com(imol,:)
    thr_rutil(4:9,tpi2) = 0.0
!
!   update gyration tensor
    do i=mlg_limits(idx,1,tpi2),mlg_limits(idx,2,tpi2)
      hlp1 = x(i)-dum3(1)
      hlp2 = y(i)-dum3(2)
      hlp3 = z(i)-dum3(3)
      hlp4 = hlp1*hlp1
      hlp5 = hlp2*hlp2
      hlp6 = hlp3*hlp3
      thr_rutil(4,tpi2) = thr_rutil(4,tpi2) + hlp1*hlp1
      thr_rutil(5,tpi2) = thr_rutil(5,tpi2) + hlp2*hlp2
      thr_rutil(6,tpi2) = thr_rutil(6,tpi2) + hlp3*hlp3
      thr_rutil(7,tpi2) = thr_rutil(7,tpi2) + hlp1*hlp2
      thr_rutil(8,tpi2) = thr_rutil(8,tpi2) + hlp1*hlp3
      thr_rutil(9,tpi2) = thr_rutil(9,tpi2) + hlp2*hlp3
    end do
  end if
!
  if (job_flags(4).EQV..true.) then
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        tvec(:) = -bnd_params(1:3)*anint((com(imol,:)-(bnd_params(4:6)+0.5*bnd_params(1:3)))/bnd_params(1:3))
        x(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) = x(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) + tvec(1)
        y(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) = y(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) + tvec(2)
        z(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) = z(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) + tvec(3)
      else if (bnd_shape.eq.3) then
        tvec(1:2) = 0.0
        tvec(3) = -bnd_params(6)*anint((com(imol,3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6))
        z(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) = z(mlg_limits(idx,1,tpi2):mlg_limits(idx,2,tpi2)) + tvec(3)
      end if
    end if
  end if
!
  if ((job_flags(1).EQV..true.).OR.((job_flags(4).EQV..true.).AND.(bnd_type.eq.1))) then
    if (job_flags(1).EQV..true.) then
!$OMP BARRIER
!$OMP SINGLE
      rgv(imol) = sum(thr_rutil(4,1:thrdat%maxn)+thr_rutil(5,1:thrdat%maxn)+thr_rutil(6,1:thrdat%maxn))
      rgpcs(imol,1,1) = sum(thr_rutil(4,1:thrdat%maxn))
      rgpcs(imol,2,2) = sum(thr_rutil(5,1:thrdat%maxn))
      rgpcs(imol,3,3) = sum(thr_rutil(6,1:thrdat%maxn))
      rgpcs(imol,1,2) = sum(thr_rutil(7,1:thrdat%maxn))
      rgpcs(imol,2,1) = rgpcs(imol,1,2)
      rgpcs(imol,1,3) = sum(thr_rutil(8,1:thrdat%maxn))
      rgpcs(imol,3,1) = rgpcs(imol,1,3)
      rgpcs(imol,2,3) = sum(thr_rutil(9,1:thrdat%maxn))
      rgpcs(imol,3,2) = rgpcs(imol,2,3)
!$OMP END SINGLE
    end if
!$OMP SINGLE
    if (job_flags(1).EQV..true.) then
      rgv(imol) = sqrt(immol*rgv(imol))
      rgpcs(imol,:,:) = immol*rgpcs(imol,:,:)
    end if
    if ((job_flags(4).EQV..true.).AND.(bnd_type.eq.1)) then 
      com(imol,:) = com(imol,:) + tvec(:)
      comm(imol,:) = comm(imol,:) + tvec(:)
    end if
    if ((job_flags(1).EQV..true.).AND.(job_flags(3).EQV..true.)) then
      j = 3
      rgten(:,:) = rgpcs(imol,:,:)
      call mat_diag(j,rgten,dum1,evmat)
      rgpcs(imol,:,:) = transpose(evmat(:,:))
      rgevs(imol,:) = dum1(:)
    end if
!$OMP END SINGLE NOWAIT
  end if
!$OMP BARRIER
!
end
!
!---------------------------------------------------------------------------------------------------------
!
! makexyz has limited parallelizability (atom level) due to strong data dependency (worst case scenario of zero parallelizability)
! this means that many threads may have nothing to do at the various levels (depends on block size and # threads)
!
! this subroutine takes 3 alternative routes: single thread, all threads with barriers per hierarchical unit, subteam 
! the subteam route creates overhead and is sensitive to affinity settings, but needed to avoid unreasonable slowdown
! at large thread numbers
!
subroutine makexyz_threads(idx,tpi)
!
  use threads
  use zmatrix
!
  implicit none
!
  integer, INTENT(IN):: idx,tpi
!
  integer i,j,imol,chiral,i2,ii,i3,i4
  RTYPE bl,ba,baodi
  integer tpi2,npt,OMP_GET_THREAD_NUM
!  integer(KIND=8) tts(2)
!
!  if (tpi.eq.1) call System_Clock(count=tts(1))
  imol = mlg_limits(idx,5,tpi)
!
  if (mlg_limits(idx,12,1).eq.1) then
!$OMP SINGLE
    tpi2 = 1
    do i=1,mlg_limits(idx,11,1) ! must be the same for all threads
      do j=1,mlg_ixs(1,i,idx,tpi2)
        ii = mlg_ixs(j+1,i,idx,tpi2)
        i2 = iz(1,ii)
        i3 = iz(2,ii)
        i4 = iz(3,ii)
        chiral = iz(4,ii)
        bl = blen(ii)
        ba = bang(ii)
        baodi = ztor(ii)
        call genxyz(ii,i2,bl,i3,ba,i4,baodi,chiral)
      end do
    end do
!$OMP END SINGLE
  else if (mlg_limits(idx,12,1).le.thrdat%maxn) then
    do i=1,mlg_limits(idx,11,1) ! must be the same for all threads
      do j=1,mlg_ixs(1,i,idx,tpi)
        ii = mlg_ixs(j+1,i,idx,tpi)
        i2 = iz(1,ii)
        i3 = iz(2,ii)
        i4 = iz(3,ii)
        chiral = iz(4,ii)
        bl = blen(ii)
        ba = bang(ii)
        baodi = ztor(ii)
        call genxyz(ii,i2,bl,i3,ba,i4,baodi,chiral)
      end do
!$OMP BARRIER
    end do
  else
    npt = mlg_limits(idx,12,1)
!$OMP SINGLE
#ifdef DISABLE_OPENMP4
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(npt) PRIVATE(tpi2,i,j,i2,i3,i4,chiral,bl,ba,baodi,ii)
#else
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(npt) PRIVATE(tpi2,i,j,i2,i3,i4,chiral,bl,ba,baodi,ii) PROC_BIND(CLOSE)
#endif
    tpi2 = OMP_GET_THREAD_NUM()+1
!
    do i=1,mlg_limits(idx,11,1) ! must be the same for all threads
      do j=1,mlg_ixs(1,i,idx,tpi2)
        ii = mlg_ixs(j+1,i,idx,tpi2)
        i2 = iz(1,ii)
        i3 = iz(2,ii)
        i4 = iz(3,ii)
        chiral = iz(4,ii)
        bl = blen(ii)
        ba = bang(ii)
        baodi = ztor(ii)
        call genxyz(ii,i2,bl,i3,ba,i4,baodi,chiral)
      end do
!$OMP BARRIER
    end do
!$OMP END PARALLEL
!$OMP END SINGLE
  end if
!  if (tpi.eq.1) call System_Clock(count=tts(2))
!  if (tpi.eq.1) write(*,*) 1.0e3*(tts(2)-tts(1))/thrdat%rate
!
end
!
!---------------------------------------------------------------------------------------------------------
!
subroutine zmatfyc_threads(tpi,mode)
!
  use fyoc
  use sequen
  use zmatrix
  use threads
!
  implicit none
!
  integer, INTENT(IN):: tpi,mode
!
  integer ff,ll,yy,ww
  integer i,j
  integer(KIND=8) ttimer
!     
  if (thr_dlb(9,1).gt.0) then
    if (tpi.eq.1) thr_dlb(9,2) = thr_dlb(9,2) + 1
    call System_Clock(count=ttimer)
    thr_timings(21,tpi) = thr_timings(21,tpi) + ttimer
  end if
  do i=thr_limits(45,tpi),thr_limits(46,tpi)
!
    if (mode.eq.2) then
      if ((fline(i).gt.0).AND.(fline2(i).gt.0)) then 
        phish(i) = ztor(fline2(i)) - ztor(fline(i))
        if (phish(i).gt.180.0) phish(i) = phish(i) - 360.0
        if (phish(i).le.-180.0) phish(i) = phish(i) + 360.0
      end if
      if ((yline(i).gt.0).AND.(yline2(i).gt.0)) then 
        if (ztor(yline(i)).lt.0.0) then
          psish(i) = ztor(yline2(i)) - (ztor(yline(i)) + 180.0)
        else
          psish(i) = ztor(yline2(i)) - (ztor(yline(i)) - 180.0)
        end if
        if (psish(i).gt.180.0) psish(i) = psish(i) - 360.0
        if (psish(i).le.-180.0) psish(i) = psish(i) + 360.0
      end if
    end if

! exclude inactive residues
    if (notors(i).EQV..true.) cycle
!
! first phi
    ff = fline(i)
    if (ff.gt.0) then
      phi(i) = ztor(ff)
      if (phi(i).gt.180.0) phi(i) = phi(i) - 360.0
      if (phi(i).lt.-180.0) phi(i) = phi(i) + 360.0
    end if
!
! next set psi
    yy = yline(i)
    if (yy.gt.0) then
      psi(i) = ztor(yy)
    end if
!
! reset the psi-angle according to zmatrix format
    if (yy.gt.0) then
      if (psi(i) .lt. 0.0d0) then
        psi(i) = psi(i) + 180.0d0
      else
        psi(i) = psi(i) - 180.0d0
      end if
      if (psi(i).gt.180.0) psi(i) = psi(i) - 360.0
      if (psi(i).lt.-180.0) psi(i) = psi(i) + 360.0
    end if
!
! get nucleic acid backbone angles from z-matrix
    do j = 1,nnucs(i)
      ll = nucsline(j,i)
      nucs(j,i) = ztor(ll)
      if (nucs(j,i).gt.180.0) nucs(j,i) = nucs(j,i) - 360.0
      if (nucs(j,i).lt.-180.0) nucs(j,i) = nucs(j,i) + 360.0
    end do
!
! get chi angles from z-matrix
    do j = 1,nchi(i)
      ll = chiline(j,i)
      chi(j,i) = ztor(ll)
      if (chi(j,i).gt.180.0) chi(j,i) = chi(j,i) - 360.0
      if (chi(j,i).lt.-180.0) chi(j,i) = chi(j,i) + 360.0
    end do
!
! omega angles
    ww = wline(i)
    if (ww.gt.0) then
      omega(i) = ztor(ww)
      if (omega(i).gt.180.0) omega(i) = omega(i) - 360.0
      if (omega(i).lt.-180.0) omega(i) = omega(i) + 360.0
    end if
!
  end do
  if (thr_dlb(9,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(22,tpi) = thr_timings(22,tpi) + ttimer
  end if
!
end
!
!---------------------------------------------------------------------
!
subroutine fyczmat_threads(tpi)
!
  use fyoc
  use sequen
  use zmatrix
  use threads
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer ff,ll,yy,ww,f2,y2
  integer i,j
  integer(KIND=8) ttimer
!     
  if (thr_dlb(9,1).gt.0) then
    if (tpi.eq.1) thr_dlb(9,2) = thr_dlb(9,2) + 1
    call System_Clock(count=ttimer)
    thr_timings(21,tpi) = thr_timings(21,tpi) + ttimer
  end if
  do i=thr_limits(45,tpi),thr_limits(46,tpi)
!
! exclude inactive residues
    if (notors(i).EQV..true.) cycle
!
! first phi
    ff = fline(i)
    f2 = fline2(i)
    if (ff.gt.0) then
      ztor(ff) = phi(i)
    end if
!
! set other torsions that depend on phi
    if ((f2.gt.0).AND.(ff.gt.0)) then
      ztor(f2) = ztor(ff) + phish(i)
      if (ztor(f2).gt.180.0) ztor(f2) = ztor(f2) - 360.0
      if (ztor(f2).lt.-180.0) ztor(f2) = ztor(f2) + 360.0
    end if
!
! next set psi
    yy= yline(i)
    y2 = yline2(i)
    if (yy.gt.0) then
      ztor(yy) = psi(i)
    end if
!
! followed by torsions that depend on psi
    if ((y2.gt.0).AND.(yy.gt.0)) then
      ztor(y2) = ztor(yy) + psish(i)
      if (ztor(y2).gt.180.0) ztor(y2) = ztor(y2) - 360.0
      if (ztor(y2).lt.-180.0) ztor(y2) = ztor(y2) + 360.0
    end if
!
! reset the psi-angle according to zmatrix format
    if (yy.gt.0) then
      if(ztor(yy) .lt. 0.0d0) then
        ztor(yy) = ztor(yy) + 180.0d0
      else
        ztor(yy) = ztor(yy) - 180.0d0
      end if
    end if
!
! put nucleic acid backbone angles into the z-matrix
    do j = 1,nnucs(i)
      ll = nucsline(j,i)
      ztor(ll) = nucs(j,i)
    end do
!
! copy the chi angles into the z-matrix
    do j = 1,nchi(i)
      ll = chiline(j,i)
      ztor(ll) = chi(j,i)
    end do
!
! omega angles
    ww = wline(i)
    if (ww.gt.0) then
      ztor(ww) = omega(i)
    end if
!
  end do
!
  if (thr_dlb(9,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(22,tpi) = thr_timings(22,tpi) + ttimer
  end if
!
end
!
!-----------------------------------------------------------------------
!
! shorten the list, calculate parameters, find center, and sort
!
subroutine clusters_postproc_threads(it,nba,nsnps,upornot,tpi,thresh)
!
  use clusters
  use threads
  use interfaces
!
  implicit none
!
  integer, INTENT(IN):: tpi,nsnps
  logical, INTENT(IN):: upornot
  RTYPE, INTENT(IN):: thresh
!
  integer nba
  type(t_scluster) it(nba)
  integer, ALLOCATABLE:: fixlims(:,:)          ! boundaries for threads
  integer tpn,i,k,kk,ixx,klo,khi,centeri,buffi(2),l,ll,j,jj,lastc
  integer OMP_GET_NUM_THREADS
  RTYPE mind,rdv,frac(2)
  type(t_scluster), ALLOCATABLE:: copyit(:)
  logical lastcomp,firstcomp,lhlp
!
  tpn = OMP_GET_NUM_THREADS()
  allocate(fixlims(6,tpn))
  frac(1) = (1.0*nba)/(1.0*tpn)
  frac(2) = (1.0*nsnps)/(1.0*tpn)
  buffi(1:2) = 0
  do i=1,tpn
    fixlims(1,i) = min(nba,buffi(1)) + 1
    buffi(1) = nint(i*frac(1))
    fixlims(2,i) = max(1,buffi(1))
    fixlims(3,i) = min(nsnps,buffi(2)) + 1
    buffi(2) = nint(i*frac(2))
    fixlims(4,i) = max(1,buffi(2))
  end do
!
  do k=1,2
    ixx = 1
    do i=2,tpn
      if (fixlims(2*k,i).ge.fixlims(2*k-1,i)) then
        if (fixlims(2*k-1,i).le.fixlims(2*k,ixx)) then
          fixlims(2*k-1,i) = fixlims(2*k-1,i) + 1
        else
          ixx = i
        end if
      end if
    end do
  end do
!
  do i=fixlims(1,tpi),fixlims(2,tpi)
    if (it(i)%nmbrs.gt.0) then
      call cluster_calc_params(it(i),thresh) ! only modifies fields radius, diam, and quality in it(i)
    end if
  end do
!
!$OMP SINGLE
  if (nba.gt.1) allocate(thr_hlper(nba,5))
!$OMP END SINGLE
  if (cmode.ne.7) then
    do i=fixlims(1,tpi),fixlims(2,tpi)
      if ((it(i)%nmbrs.eq.1).OR.(it(i)%nmbrs.eq.2)) it(i)%center = minval(it(i)%snaps(1:it(i)%nmbrs))
      it(i)%nodewt(2) = HUGE(it(i)%nodewt(2)) ! hi-jack 
    end do
!$OMP BARRIER
    klo = 1
    fixlims(5,tpi) = 0
    fixlims(6,tpi) = nba
    thr_rutil(20,tpi) = 0.0
    firstcomp = .false.
    do i=1,nba
      if (klo.gt.fixlims(4,tpi)) then
        fixlims(6,tpi) = i - 1
        exit
      end if
      khi = klo + it(i)%nmbrs
      if (khi.le.fixlims(3,tpi)) then
        klo = khi
        cycle
      end if
      lastcomp = .false.
      mind = HUGE(mind)
      if (fixlims(5,tpi).eq.0) fixlims(5,tpi) = i ! first cluster to be touched by this thread

      if ((fixlims(3,tpi)-klo+1).le.1) then ! tpi has the first snapshot in this cluster -> use native arrays
        if (fixlims(5,tpi).eq.i) then ! this was also tpi's first cluster overall
          firstcomp = .true. ! to indicate that tpi's lower bound in cluster number is OK as it is
          thr_rutil(19,tpi) = HUGE(mind)
          thr_rutil(20,tpi) = 1.0*i
        end if
        if (it(i)%nmbrs.gt.2) then
          do k=max(1,fixlims(3,tpi)-klo+1),min(it(i)%nmbrs,fixlims(4,tpi)-klo+1)
            call snap_to_cluster_d(rdv,it(i),it(i)%snaps(k))
            if (rdv.lt.mind) then
              it(i)%center = it(i)%snaps(k)
              mind = rdv
              it(i)%nodewt(2) = rdv
            end if
          end do
        end if
      else ! tpi does not have the first snapshot in this cluster -> use temp arrays (there can only be on of those per thread)
        if (fixlims(5,tpi).eq.i) then ! this was also tpi's first cluster overall
          thr_rutil(19,tpi) = HUGE(mind)
          thr_rutil(20,tpi) = 1.0*i
        end if
        if (it(i)%nmbrs.gt.2) then
          do k=max(1,fixlims(3,tpi)-klo+1),min(it(i)%nmbrs,fixlims(4,tpi)-klo+1)
            call snap_to_cluster_d(rdv,it(i),it(i)%snaps(k))
            if (rdv.lt.thr_rutil(19,tpi)) then
              thr_rutil(19,tpi) = rdv
              thr_rutil(21,tpi) = 1.0*it(i)%snaps(k)
            end if
          end do
        end if
      end if
      if (it(i)%nmbrs.le.fixlims(4,tpi)-klo+1) lastcomp = .true. ! i ended this cluster
      klo = khi
    end do
!$OMP BARRIER
!$OMP SINGLE
    do i=1,thrdat%maxn
      lastc = nint(thr_rutil(20,i))
      if ((lastc.gt.0).AND.(lastc.le.nba)) then
        if ((it(lastc)%nmbrs.gt.2).AND.(thr_rutil(19,i).lt.it(lastc)%nodewt(2))) then
          it(lastc)%nodewt(2) = thr_rutil(19,i)
          it(lastc)%center = nint(thr_rutil(21,i))
        end if
      end if
    end do
!$OMP END SINGLE NOWAIT
  else
    klo = 1
    fixlims(5,tpi) = 0
    fixlims(6,tpi) = nba
    firstcomp = .false.
    do i=1,nba
      if (klo.gt.fixlims(4,tpi)) then
        fixlims(6,tpi) = i - 1
        exit
      end if
      khi = klo + it(i)%nmbrs
      if (khi.le.fixlims(3,tpi)) then
        klo = khi
        cycle
      end if
      lastcomp = .false.
      if (fixlims(5,tpi).eq.0) fixlims(5,tpi) = i ! first cluster to be touched by this thread
      if (fixlims(5,tpi).eq.i) then
        if (it(i)%nmbrs.gt.1) then
          if ((fixlims(3,tpi)-klo+1).le.1) firstcomp = .true. ! i started this cluster
        else
          firstcomp = .true.
        end if
      end if
      if (it(i)%nmbrs.le.fixlims(4,tpi)-klo+1) lastcomp = .true. ! i ended this cluster
      klo = khi
    end do
  end if ! whether in cmode 7
!
  if (tpi.eq.1) fixlims(5,tpi) = 1
  if ((fixlims(5,tpi).le.0).OR.(fixlims(6,tpi).le.0)) then
    fixlims(5,tpi) = 1
    fixlims(6,tpi) = 0
  end if
  if (nba.gt.1) then
!   sort according to size (a threaded NlogN library sort would come in handy)
    do i=fixlims(1,tpi),fixlims(2,tpi)
      thr_hlper(i,3) = i
      thr_hlper(i,1) = it(i)%nmbrs
    end do
    klo = 1
    khi = nba
    centeri = nba
    if (fixlims(6,tpi).ge.fixlims(5,tpi)) allocate(copyit(fixlims(6,tpi)-fixlims(5,tpi)+1))
  end if
!$OMP BARRIER
  if (nba.le.1) return
!$OMP SINGLE
  lhlp = upornot
  call merge_sort(centeri,lhlp,thr_hlper(:,1),thr_hlper(:,2),klo,khi,thr_hlper(:,3),thr_hlper(:,4))
!$OMP END SINGLE
  if ((firstcomp.EQV..false.).AND.(lastcomp.EQV..false.).AND.(fixlims(5,tpi).eq.fixlims(6,tpi))) then
    fixlims(5,tpi) = 1
    fixlims(6,tpi) = 0
  else if (firstcomp.EQV..false.) then
    fixlims(5,tpi) = fixlims(5,tpi) + 1
  end if
!
  do i=fixlims(5,tpi),fixlims(6,tpi)
    k = i - fixlims(5,tpi) + 1
    call copy_cluster(it(thr_hlper(i,4)),copyit(k))
  end do
!$OMP BARRIER
  do i=fixlims(5,tpi),fixlims(6,tpi)
    k = i - fixlims(5,tpi) + 1
    call copy_cluster(copyit(k),it(i))
  end do
  if (allocated(copyit).EQV..true.) deallocate(copyit)
!$OMP BARRIER
  if (resort_clustering.EQV..true.) then
    allocate(copyit(10))
    i = nba
    kk = 0
    do while (i.ge.1)
      k = i-1
      khi = i
      do while (it(i)%nmbrs.eq.it(k)%nmbrs)
        k = k - 1
        if (k.le.0) exit
      end do
      klo = k+1
      if (khi.gt.(klo+1)) then
        kk = kk + 1
        if (mod(kk,tpn).eq.(tpi-1)) then
          if ((khi-klo+1).gt.size(copyit)) then
            deallocate(copyit)
            allocate(copyit(khi-klo+1))
          end if
!         resort size ties by centroid (increasing)
          do j=klo,khi
            thr_hlper(j,3) = j
            thr_hlper(j,1) = it(j)%center
          end do
          l = 1
          ll = khi-klo+1
          centeri = ll
          lhlp = .true.
          call merge_sort(centeri,lhlp,thr_hlper(klo:khi,1),thr_hlper(klo:khi,2),l,ll,thr_hlper(klo:khi,3),thr_hlper(klo:khi,4))
          do j=1,centeri
            jj = thr_hlper(j+klo-1,4)
            call copy_cluster(it(jj),copyit(j))
          end do
          do j=1,centeri
            jj = j + klo - 1
            call copy_cluster(copyit(j),it(jj))
          end do
        end if
      else if (khi.eq.(klo+1)) then
        kk = kk + 1
        if (mod(kk,tpn).eq.(tpi-1)) then
          if (it(khi)%center.lt.it(klo)%center) then
            call copy_cluster(it(klo),copyit(1))
            call copy_cluster(it(khi),it(klo))
            call copy_cluster(copyit(1),it(khi))
          end if
        end if
      end if
      i = k
      if (i.le.1) exit
    end do
    if (allocated(copyit).EQV..true.) deallocate(copyit)
  end if
!$OMP BARRIER
!$OMP SINGLE
  deallocate(thr_hlper)
!$OMP END SINGLE
!
end
!
!------------------------------------------------------------------------
!
! test some threaded implementations in isolated fashion
!
subroutine threads_test()
!
  use atoms
  use molecule
  use threads
  use iounit
  use zmatrix
  use energies
  use forces
  use system
  use movesets
  use cutoffs
  use ewalds
#ifdef LINK_FFTW
  use, intrinsic :: ISO_C_BINDING
#endif
!
  implicit none
!
#ifdef LINK_FFTW
#ifdef ORACLE_FORTRAN
#include "fftw3_SUN.f03"
#else
#include "fftw3.f03"
#endif
#endif
  integer azero,tpi,tpn,i,j,k,OMP_GET_THREAD_NUM,imol,mm,tpnbu,bu1,lo,hi,mred,bigsz,smallsz
  RTYPE random,force3,fr1bu,fr2bu,fbu1,fr1,fr2
  logical jobflags(4),atrue,afalse
  RTYPE, ALLOCATABLE:: refvals(:,:),evals(:,:)
  integer, ALLOCATABLE:: refivals(:,:)
  integer(KIND=8) ttt(thrdat%maxn,4)
!
  write(ilog,*) 
  write(ilog,*) '*** BEGINNING OF TEST OF MULTI-THREADED ROUTINES ***'
  write(ilog,*)
  write(ilog,*) 'CAMPARI is always going to use ',thrdat%maxn,' threads. It is recommended to vary this number through the &
 &key-file. Errors and race conditions are best identified in oversubscribed usage.'
  write(ilog,*) 
!
 55 format(2x,'Threaded computation by ',a,' of ',a,' yielded a MUD of ',g16.8,a,' and a max error of ',g16.8,a,'.')
  azero = 0
  atrue = .true.
  afalse = .false.
  tpnbu = thrdat%maxn
  bu1 = ens%sysfrz
  fbu1 = ens%insK
  mred = 0
  bigsz = 100000
  smallsz = 100
!
  if (allocated(thr_rhlper).EQV..true.) deallocate(thr_rhlper)
!
  allocate(thr_rhlper(bigsz,3))
  allocate(refvals(bigsz,thrdat%maxn+1))
  do j=1,thrdat%maxn
    do i=1,smallsz
      refvals(i,j) = random()
    end do 
  end do
  refvals(1:smallsz,thrdat%maxn+1) = sum(refvals(1:smallsz,1:thrdat%maxn),dim=2)
  call OMP_SET_NUM_THREADS(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,j,k)
  mm = smallsz
  j = 1
  tpi = OMP_GET_THREAD_NUM()+1
  call thr_freduction(mred,refvals(1:mm,tpi),mm,j,tpi)
!$OMP END PARALLEL

  mm = smallsz
  fr1 = sum(abs(refvals(1:mm,thrdat%maxn+1)-refvals(1:mm,1)))/(1.0*mm)
  fr2 = maxval(abs(refvals(1:mm,thrdat%maxn+1)-refvals(1:mm,1)))
  write(ilog,55) 'thr_freduction','a small FP vector',fr1,' a.u.',fr2,' a.u'
!
  do j=1,thrdat%maxn
    do i=1,bigsz
      refvals(i,j) = random()
    end do 
  end do
  refvals(:,thrdat%maxn+1) = sum(refvals(:,1:thrdat%maxn),dim=2)
!
  if (allocated(thr_rhlper).EQV..true.) deallocate(thr_rhlper)
!
  allocate(thr_rhlper(bigsz,3))

  thr_rhlper(:,2) = 0.0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,j,k)
  mm = bigsz
  i = 67
  j = 1
  k = 2
  tpi = OMP_GET_THREAD_NUM()+1
  call thr_combinevec(refvals(1:mm,tpi),mm,j,k,i,tpi)
!$OMP END PARALLEL
  mm = bigsz
  thr_limits(67,:) = 1
  thr_limits(68,:) = 0
  fr1 = sum(abs(refvals(1:mm,thrdat%maxn+1)-thr_rhlper(1:mm,2)))/(1.0*mm)
  fr2 = maxval(abs(refvals(1:mm,thrdat%maxn+1)-thr_rhlper(1:mm,2)))
  write(ilog,55) 'thr_combinevec','a big FP vector',fr1,' a.u.',fr2,' a.u'
  deallocate(refvals)
  deallocate(thr_rhlper)
!
  if (allocated(thr_hlper).EQV..true.) deallocate(thr_hlper)
!
  allocate(thr_hlper(bigsz,3))
  allocate(refivals(bigsz,thrdat%maxn+1))
  do j=1,thrdat%maxn
    do i=1,bigsz
      refivals(i,j) = nint(random()*10.0)
    end do 
  end do
  refivals(:,thrdat%maxn+1) = sum(refivals(:,1:thrdat%maxn),dim=2)
  thr_hlper(:,2) = 0
  call OMP_SET_NUM_THREADS(thrdat%maxn)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,j,k)
  mm = smallsz
  i = 67
  j = 1
  k = 2
  tpi = OMP_GET_THREAD_NUM()+1
  call thr_combineivec(refivals(1:mm,tpi),mm,j,k,i,tpi)
!$OMP END PARALLEL
  mm = smallsz
  thr_limits(67,:) = 1
  thr_limits(68,:) = 0
  fr1 = 1.0*sum(abs(refivals(1:mm,thrdat%maxn+1)-thr_hlper(1:mm,2)))/(1.0*mm)
  fr2 = 1.0*maxval(abs(refivals(1:mm,thrdat%maxn+1)-thr_hlper(1:mm,2)))
  write(ilog,55) 'thr_combineivec','a small integer vector',fr1,' a.u.',fr2,' a.u'
!
  thr_hlper(:,2) = 0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,j,k)
  mm = bigsz
  i = 67
  j = 1
  k = 2
  tpi = OMP_GET_THREAD_NUM()+1
  call thr_combineivec(refivals(1:mm,tpi),mm,j,k,i,tpi)
!$OMP END PARALLEL
  mm = bigsz
  thr_limits(67,:) = 1
  thr_limits(68,:) = 0
  fr1 = sum(abs(refivals(1:mm,thrdat%maxn+1)-thr_hlper(1:mm,2)))/(1.0*mm)
  fr2 = maxval(abs(refivals(1:mm,thrdat%maxn+1)-thr_hlper(1:mm,2)))
  write(ilog,55) 'thr_combineivec','a big integer vector',fr1,' a.u.',fr2,' a.u'
  deallocate(refivals)
  deallocate(thr_hlper)
!
  if ((use_dyn.EQV..true.).AND.(fycxyz.eq.2)) then
    xref(1:n) = cart_v(1:n,1)
    yref(1:n) = cart_v(1:n,2)
    zref(1:n) = cart_v(1:n,3)
    do i=1,n
      cart_v(i,1) = 0.5*(random()-0.5)
      cart_v(i,2) = 0.5*(random()-0.5)
      cart_v(i,3) = 0.5*(random()-0.5)
    end do
    do imol=1,nmol
      call update_comv(imol)
    end do
    allocate(refvals(n,9))
    if (nmlgs.gt.0) then
      do mm=1,nmlgs
        imol = mlg_limits(mm,5,1)
        refvals(mm,4) = dc_di(imol)%v(1)
        refvals(mm,5) = dc_di(imol)%v(2)
        refvals(mm,6) = dc_di(imol)%v(3)
      end do
      jobflags(:) = .false.
      jobflags(3) = .true.
      call OMP_SET_NUM_THREADS(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i)
      tpi = OMP_GET_THREAD_NUM()+1
      do mm=1,nmlgs
        call molops_threads(mm,tpi,jobflags)
        imol = mlg_limits(mm,5,tpi)
!$OMP SINGLE
        refvals(mm,7:9) = abs(dc_di(imol)%v(1:3)-refvals(mm,4:6))
!$OMP END SINGLE
      end do
!$OMP END PARALLEL
      write(ilog,55) 'molops_threads','center of mass velocities',&
 &                    sum(refvals(1:nmlgs,7:9))/(3.0*nmlgs),' A/ps',maxval(refvals(1:nmlgs,7:9)),' A/ps'
      do mm=1,nmlgs
        imol = mlg_limits(mm,5,1)
        dc_di(imol)%v(1) = refvals(mm,4)
        dc_di(imol)%v(2) = refvals(mm,5)
        dc_di(imol)%v(3) = refvals(mm,6)
      end do
      write(ilog,*)
    end if
    ens%insK = 10.
    do tpn=1,3
      ens%sysfrz = tpn
      write(ilog,*) 'Now testing drift_removal(...) with FMCSC_SYSFRZ ',ens%sysfrz,':'
      refvals(1:n,4:6) = cart_v(1:n,1:3)
      call drift_removal(fr1,fr2,azero)
      refvals(1:n,1:3) = cart_v(1:n,1:3)
      cart_v(1:n,1:3) = refvals(1:n,4:6)
      fr1bu = fr1
      fr2bu = fr2
      call OMP_SET_NUM_THREADS(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,fr1,fr2,i)
      tpi = OMP_GET_THREAD_NUM()+1
      call drift_removal(fr1,fr2,tpi)
!$OMP BARRIER
!$OMP SINGLE
      write(ilog,55) 'drift_removal','translational and angular drift',&
 &                    (abs(fr1-fr1bu)+abs(fr2-fr2bu))/2.0,' %',max(abs(fr1-fr1bu),abs(fr2-fr2bu)),' %'
      refvals(1:n,7:9) = abs(cart_v(1:n,1:3)-refvals(1:n,1:3))
      write(ilog,55) 'drift_removal','Cartesian velocities',&
 &                    sum(refvals(1:n,7:9))/(3.0*n),' A/ps',maxval(refvals(1:n,7:9)),' A/ps'
!$OMP END SINGLE
!$OMP END PARALLEL
      cart_v(1:n,1:3) = refvals(1:n,4:6)
    end do
    cart_v(1:n,1) = xref(1:n)
    cart_v(1:n,2) = yref(1:n)
    cart_v(1:n,3) = zref(1:n)
    deallocate(refvals)
    write(ilog,*)
  else if ((use_dyn.EQV..true.).AND.(fycxyz.eq.1)) then
    do imol=1,nmol
      xref(imol) = dc_di(imol)%v(1)
      yref(imol) = dc_di(imol)%v(2)
      zref(imol) = dc_di(imol)%v(3)
    end do
    do imol=1,nmol
      dc_di(imol)%v(1) = 0.5*(random()-0.5)
      dc_di(imol)%v(2) = 0.5*(random()-0.5)
      dc_di(imol)%v(3) = 0.5*(random()-0.5)
    end do
    allocate(refvals(nmol,9))
    ens%insK = 10.
    do tpn=1,3
      ens%sysfrz = tpn
      write(ilog,*) 'Now testing drift_removal(...) with FMCSC_SYSFRZ ',ens%sysfrz,':'
      do imol=1,nmol
        refvals(imol,4) = dc_di(imol)%v(1)
        refvals(imol,5) = dc_di(imol)%v(2)
        refvals(imol,6) = dc_di(imol)%v(3)
      end do
      call drift_removal(fr1,fr2,azero)
      do imol=1,nmol
        refvals(imol,1) = dc_di(imol)%v(1)
        refvals(imol,2) = dc_di(imol)%v(2)
        refvals(imol,3) = dc_di(imol)%v(3)
        dc_di(imol)%v(1) = refvals(imol,4)
        dc_di(imol)%v(2) = refvals(imol,5)
        dc_di(imol)%v(3) = refvals(imol,6)
      end do
      fr1bu = fr1
      fr2bu = fr2
      call OMP_SET_NUM_THREADS(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,fr1,fr2,i)
      tpi = OMP_GET_THREAD_NUM()+1
      call drift_removal(fr1,fr2,tpi)
!$OMP BARRIER
!$OMP SINGLE
      write(ilog,55) 'drift_removal','translational and angular drift',&
 &                    (abs(fr1-fr1bu)+abs(fr2-fr2bu))/2.0,' %',max(abs(fr1-fr1bu),abs(fr2-fr2bu)),' %'
      do imol=1,nmol
        refvals(imol,7:9) = abs(dc_di(imol)%v(1:3)-refvals(imol,1:3))
      end do
      write(ilog,55) 'drift_removal','center of mass velocities',&
 &                    sum(refvals(1:nmol,7:9))/(3.0*nmol),' A/ps',maxval(refvals(1:nmol,7:9)),' A/ps'
!$OMP END SINGLE
!$OMP END PARALLEL
      do imol=1,nmol
        dc_di(imol)%v(1) = refvals(imol,4)
        dc_di(imol)%v(2) = refvals(imol,5)
        dc_di(imol)%v(3) = refvals(imol,6)
      end do
    end do
    do imol=1,nmol
      dc_di(imol)%v(1) = xref(imol)
      dc_di(imol)%v(2) = yref(imol)
      dc_di(imol)%v(3) = zref(imol)
    end do
    deallocate(refvals)
    write(ilog,*)
  end if
!
  ens%sysfrz = bu1
  ens%insK = fbu1
!
! forces and energies
  if (use_dyn.EQV..true.) then
    allocate(refvals(n,9))
    allocate(evals(MAXENERGYTERMS,3))
    fr1 = force3(evals(:,1),evals(:,2),evals(:,3),atrue)
    refvals(1:n,1:3) = cart_f(1:n,1:3)
    refvals(1:n,4:6) = cart_f_tr(1:n,1:3)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    k = MAXENERGYTERMS
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,atrue)
    tpi = OMP_GET_THREAD_NUM()+1
    call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
!$OMP END PARALLEL
    write(ilog,55) 'force3_threads with NBL update','net Cartesian forces',&
 &  sum(abs(cart_f(1:n,1:3)-refvals(1:n,1:3))/(3.0*n)),' kcal/mol/A',maxval(abs(cart_f(1:n,1:3)-refvals(1:n,1:3))),' kcal/mol/A'
    write(ilog,55) 'force3_threads with NBL update','TR Cartesian forces',&
 &  sum(abs(cart_f_tr(1:n,1:3)-refvals(1:n,4:6))/(3.0*n)),' kcal/mol/A',maxval(abs(cart_f_tr(1:n,1:3)-refvals(1:n,4:6))),&
 &' kcal/mol/A'
    write(ilog,55) 'force3_threads with NBL update','net energy terms',&
 &  sum(abs(esterms(1:k)-evals(1:k,1))/(1.0*k)),' kcal/mol',maxval(abs(esterms(1:k)-evals(1:k,1))),' kcal/mol'
    write(ilog,55) 'force3_threads with NBL update','TR energy terms',&
 &  sum(abs(esterms_tr(1:k)-evals(1:k,2))/(1.0*k)),' kcal/mol',maxval(abs(esterms_tr(1:k)-evals(1:k,2))),' kcal/mol'
    write(ilog,55) 'force3_threads with NBL update','LR energy terms',&
 &  sum(abs(esterms_lr(1:k)-evals(1:k,3))/(1.0*k)),' kcal/mol',maxval(abs(esterms_lr(1:k)-evals(1:k,3))),' kcal/mol'
    bu1 = nbl_up
    nbl_up = 1000
    fr1 = force3(evals(:,1),evals(:,2),evals(:,3),afalse)
    refvals(1:n,1:3) = cart_f(1:n,1:3)
    refvals(1:n,4:6) = cart_f_tr(1:n,1:3)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    k = MAXENERGYTERMS
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,atrue)
    tpi = OMP_GET_THREAD_NUM()+1
    call force3_threads(esterms,evals(:,2),evals(:,3),afalse,esave)
!$OMP END PARALLEL
    write(ilog,55) 'force3_threads without NBL update','net Cartesian forces',&
 &  sum(abs(cart_f(1:n,1:3)-refvals(1:n,1:3))/(3.0*n)),' kcal/mol/A',maxval(abs(cart_f(1:n,1:3)-refvals(1:n,1:3))),' kcal/mol/A'
    write(ilog,55) 'force3_threads without NBL update','TR Cartesian forces',&
 &  sum(abs(cart_f_tr(1:n,1:3)-refvals(1:n,4:6))/(3.0*n)),' kcal/mol/A',maxval(abs(cart_f_tr(1:n,1:3)-refvals(1:n,4:6))),&
 &' kcal/mol/A'
    write(ilog,55) 'force3_threads without NBL update','net energy terms',&
 &  sum(abs(esterms(1:k)-evals(1:k,1))/(1.0*k)),' kcal/mol',maxval(abs(esterms(1:k)-evals(1:k,1))),' kcal/mol'
    deallocate(evals)
    deallocate(refvals)
    nbl_up = bu1
    write(ilog,*)
  end if

! MC energies
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    allocate(evals(MAXENERGYTERMS,3))
    call energy(evals(:,1),fr1,azero)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    k = MAXENERGYTERMS
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,atrue)
    tpi = OMP_GET_THREAD_NUM()+1
    call energy(evals(:,2),fr2,tpi)
!$OMP END PARALLEL
    write(ilog,55) 'energy (no cutoffs)','net energy terms',&
 &  sum(abs(evals(1:k,2)-evals(1:k,1))/(1.0*k)),' kcal/mol',maxval(abs(evals(1:k,2)-evals(1:k,1))),' kcal/mol'
    call energy3(evals(:,1),atrue,fr1,azero)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    k = MAXENERGYTERMS
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,atrue)
    atrue = .true.
    tpi = OMP_GET_THREAD_NUM()+1
    call energy3(evals(:,2),atrue,fr2,tpi)
!$OMP END PARALLEL
    atrue = .true.
    write(ilog,55) 'energy3 (with NBL update)','net energy terms',&
 &  sum(abs(evals(1:k,2)-evals(1:k,1))/(1.0*k)),' kcal/mol',maxval(abs(evals(1:k,2)-evals(1:k,1))),' kcal/mol'
    call energy3(evals(:,1),afalse,fr1,azero)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    k = MAXENERGYTERMS
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,atrue)
    atrue = .false.
    tpi = OMP_GET_THREAD_NUM()+1
    call energy3(evals(:,2),atrue,fr2,tpi)
!$OMP END PARALLEL
    atrue = .true.
    write(ilog,55) 'energy3 (without NBL update)','net energy terms',&
 &  sum(abs(evals(1:k,2)-evals(1:k,1))/(1.0*k)),' kcal/mol',maxval(abs(evals(1:k,2)-evals(1:k,1))),' kcal/mol'
    deallocate(evals)
    write(ilog,*)
  end if


  if (nmlgs.gt.0) then
!   coordinates generation
    xref(1:n) = x(1:n)
    yref(1:n) = y(1:n)
    zref(1:n) = z(1:n)
    allocate(refvals(n,9))
    refvals(1:n,1) = blen(1:n)
    refvals(1:n,2) = bang(1:n)
    refvals(1:n,3) = ztor(1:n)
    do mm=1,nmlgs
      imol = mlg_limits(mm,5,1)
      do i=atmol(imol,1)+3,atmol(imol,2)
        blen(i) = blen(i) + 0.02*(random()-0.5)
        bang(i) = max(1.0,min(179.0,bang(i) + 1.5*(random()-0.5)))
        ztor(i) = ztor(i) + 5.5*(random()-0.5)
      end do
      call makexyz_formol(imol)
    end do
    refvals(1:n,4) = x(1:n)
    refvals(1:n,5) = y(1:n)
    refvals(1:n,6) = z(1:n)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    lo = 1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,hi)
    tpi = OMP_GET_THREAD_NUM()+1
    do mm=1,nmlgs
      call makexyz_threads(mm,tpi)
      imol = mlg_limits(mm,5,tpi)
      hi = atmol(imol,2)-atmol(imol,1)
!$OMP BARRIER
!$OMP SINGLE
      refvals(lo:(lo+hi),7) = abs(x(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),4))
      refvals(lo:(lo+hi),8) = abs(y(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),5))
      refvals(lo:(lo+hi),9) = abs(z(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),6))
      lo = lo + hi + 1
!$OMP END SINGLE
    end do
!$OMP END PARALLEL
    write(ilog,55) 'makexyz_threads','Cartesian coordinates from Z matrix',&
 &                    sum(refvals(1:(lo-1),7:9))/(3.0*(lo-1)),' A',maxval(refvals(1:(lo-1),7:9)),' A'
    blen(1:n) = refvals(1:n,1)
    bang(1:n) = refvals(1:n,2)
    ztor(1:n) = refvals(1:n,3)
    x(1:n) = xref(1:n)
    y(1:n) = yref(1:n)
    z(1:n) = zref(1:n)
    deallocate(refvals)
    write(ilog,*)
!
!   polymeric parameters
    allocate(refvals(nmlgs*17,3))
    refvals(:,:) = 0.0
    jobflags(:) = .false.
    jobflags(1) = .true.
    do mm=1,nmlgs
      imol = mlg_limits(mm,5,1)
      call update_rigidm(imol)
      refvals(mm*17-16,1) = rgvm(imol)
      refvals(mm*17-15:mm*17-13,1) = comm(imol,1:3)
      call update_rigid(imol)
      refvals(mm*17-12,1) = rgv(imol)
      refvals(mm*17-11:mm*17-9,1) = com(imol,1:3)
      refvals(mm*17-8:mm*17-6,1) = rgpcs(imol,1,1:3)
      refvals(mm*17-5:mm*17-3,1) = rgpcs(imol,2,1:3)
      refvals(mm*17-2:mm*17,1) =   rgpcs(imol,3,1:3)
    end do
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    lo = 1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,hi)
    tpi = OMP_GET_THREAD_NUM()+1
    do mm=1,nmlgs
      call molops_threads(mm,tpi,jobflags)
      imol = mlg_limits(mm,5,tpi)
!$OMP BARRIER
!$OMP SINGLE
      jobflags(3) = .true.
      refvals(mm*17-16,2) = abs(rgvm(imol)-refvals(mm*17-16,1))
      refvals(mm*17-15:mm*17-13,2) = abs(comm(imol,:)-refvals(mm*17-15:mm*17-13,1))
!$OMP END SINGLE
      call molops_threads_geo(mm,tpi,jobflags)
      imol = mlg_limits(mm,5,tpi)
!$OMP BARRIER
!$OMP SINGLE
      jobflags(3) = .false.
      refvals(mm*17-12,3) = abs(rgv(imol)-refvals(mm*17-12,1))
      refvals(mm*17-11:mm*17-9,3) = abs(com(imol,:)-refvals(mm*17-11:mm*17-9,1))
      refvals(mm*17-8:mm*17-6,3) = abs(rgpcs(imol,1,1:3)-refvals(mm*17-8:mm*17-6,1))
      refvals(mm*17-5:mm*17-3,3) = abs(rgpcs(imol,2,1:3)-refvals(mm*17-5:mm*17-3,1))
      refvals(mm*17-2:mm*17,3) =   abs(rgpcs(imol,3,1:3)-refvals(mm*17-2:mm*17,1))
!$OMP END SINGLE  
    end do
!$OMP END PARALLEL
    write(ilog,55) 'molops_threads','polymeric descriptors',&
 &                    sum(refvals(1:17*nmlgs,2))/(4*nmlgs),' (v.u.)',maxval(refvals(1:17*nmlgs,2)),' (v.u.)'
    write(ilog,55) 'molops_threads_geo','polymeric descriptors',&
 &                    sum(refvals(1:17*nmlgs,3))/(10*nmlgs),' (v.u.)',maxval(refvals(1:17*nmlgs,3)),' (v.u.)'
    do imol=1,nmol
      call update_rigid(imol)
      call update_rigidm(imol)
    end do
    write(ilog,*) 
    deallocate(refvals)
!
!
!   Z-matrix generation
    allocate(refvals(n,9))
    jobflags(:) = .false.
    jobflags(2) = .true.
    xref(1:n) = x(1:n)
    yref(1:n) = y(1:n)
    zref(1:n) = z(1:n)
    do i=1,n
      x(i) = x(i) + 0.5*(random()-0.5)
      y(i) = y(i) + 0.5*(random()-0.5)
      z(i) = z(i) + 0.5*(random()-0.5)
    end do
    do imol=1,nmol
      call genzmat(imol)
    end do
    refvals(1:n,1) = blen(1:n)
    refvals(1:n,2) = bang(1:n)
    refvals(1:n,3) = ztor(1:n)
    call OMP_SET_NUM_THREADS(thrdat%maxn)
    lo = 1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i,hi)
    tpi = OMP_GET_THREAD_NUM()+1
    do mm=1,nmlgs
      call molops_threads(mm,tpi,jobflags)
      imol = mlg_limits(mm,5,tpi)
      hi = atmol(imol,2)-atmol(imol,1)
!$OMP BARRIER
!$OMP SINGLE
      refvals(lo:(lo+hi),4) = abs(blen(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),1))
      refvals(lo:(lo+hi),5) = abs(bang(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),2))
      refvals(lo:(lo+hi),6) = abs(ztor(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),3))
!$OMP END SINGLE
      call molops_threads_geo(mm,tpi,jobflags)
      imol = mlg_limits(mm,5,tpi)
!$OMP BARRIER
!$OMP SINGLE
      refvals(lo:(lo+hi),7) = abs(blen(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),1))
      refvals(lo:(lo+hi),8) = abs(bang(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),2))
      refvals(lo:(lo+hi),9) = abs(ztor(atmol(imol,1):atmol(imol,2))-refvals(atmol(imol,1):atmol(imol,2),3))
      lo = lo + hi + 1
!$OMP END SINGLE  
    end do
!$OMP END PARALLEL
    write(ilog,55) 'molops_threads','Z-matrix coordinates',&
 &                    sum(refvals(1:(lo-1),4:6))/(3.0*(lo-1)),' (v.u.)',maxval(refvals(1:(lo-1),4:6)),' (v.u.)'
    write(ilog,55) 'molops_threads_geo','Z-matrix coordinates',&
 &                    sum(refvals(1:(lo-1),7:9))/(3.0*(lo-1)),' (v.u.)',maxval(refvals(1:(lo-1),7:9)),' (v.u.)'
    x(1:n) = xref(1:n)
    y(1:n) = yref(1:n)
    z(1:n) = zref(1:n)
    do imol=1,nmol
      call genzmat(imol)
    end do
    write(ilog,*) 
    deallocate(refvals)
!
!   internal coordinate force et al.
    if ((use_dyn.EQV..true.).AND.(fycxyz.eq.1)) then
      allocate(refvals(3*n,9))
      do imol=1,nmol
        dc_di(imol)%f(:) = 0.
      end do
      esave = force3(esterms,esterms_tr,esterms_lr,atrue)
      tpn = 1
      do imol=1,1!nmol
        refvals(tpn:(tpn+size(dc_di(imol)%v)-1),7) = dc_di(imol)%f(:) ! this part may not be zero
        tpn = tpn + size(dc_di(imol)%v)
      end do
      do imol=1,1!nmol
        do i=1,size(dc_di(imol)%v)
          dc_di(imol)%v(i) = 0.05*(random()-0.5)
        end do
      end do
      call cart2int_f(skip_frz,azero)
      call cart2int_I(azero)
      refvals(1:n,1:3) = cart_v(1:n,1:3)
      tpn = 1
      do imol=1,1!nmol
        refvals(tpn:(tpn+size(dc_di(imol)%v)-1),4) = dc_di(imol)%f(:)
        refvals(tpn:(tpn+size(dc_di(imol)%v)-1),5) = dc_di(imol)%im(:)
        refvals(tpn:(tpn+size(dc_di(imol)%v)-1),6) = dc_di(imol)%olddat(:,2)
        dc_di(imol)%f(:) = refvals(tpn:(tpn+size(dc_di(imol)%v)-1),7)
        tpn = tpn + size(dc_di(imol)%v)
      end do
      call OMP_SET_NUM_THREADS(thrdat%maxn)
      lo = 1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,mm,imol,i)
      tpi = OMP_GET_THREAD_NUM()+1
      call cart2int_f(skip_frz,tpi)
      call cart2int_I(tpi)
!$OMP BARRIER
!$OMP SINGLE
      do imol=1,1!nmol
        k=size(dc_di(imol)%v) - 1
        refvals(lo:(lo+k),7) = abs(dc_di(imol)%f(:)-refvals(lo:(lo+k),4))
        refvals(lo:(lo+k),8) = abs(dc_di(imol)%im(:)-refvals(lo:(lo+k),5))
        refvals(lo:(lo+k),9) = abs(dc_di(imol)%olddat(:,2)-refvals(lo:(lo+k),6))
        lo = lo + k + 1  
      end do
!$OMP END SINGLE
!$OMP END PARALLEL
      write(ilog,55) 'cart2int_f','internal coordinate forces',&
 &                    sum(refvals(1:(lo-1),7))/(1.0*(lo-1)),' (v.u.)',maxval(refvals(1:(lo-1),7)),' (v.u.)'
      write(ilog,55) 'cart2int_f','diagonal elements of mass matrix',&
 &                    sum(refvals(1:(lo-1),8))/(1.0*(lo-1)),' (v.u.)',maxval(refvals(1:(lo-1),8)),' (v.u.)'
      write(ilog,55) 'cart2int_I','backup of diagonal elements of mass matrix',&
 &                    sum(refvals(1:(lo-1),9))/(1.0*(lo-1)),' (v.u.)',maxval(refvals(1:(lo-1),9)),' (v.u.)'
      write(ilog,55) 'cart2int_f','derived Cartesian velocities',&
 &  sum(abs(cart_v(1:n,1:3)-refvals(1:n,1:3))/(3.0*n)),' A/ps',maxval(abs(cart_v(1:n,1:3)-refvals(1:n,1:3))),' A/ps'
      do imol=1,nmol
        dc_di(imol)%v(:) = 0.
        dc_di(imol)%f(:) = 0.
      end do
      call cart2int_I(azero)
      deallocate(refvals)
      write(ilog,*) 
    end if
  end if
!
#ifdef LINK_FFTW
  if ((use_dyn.EQV..true.).AND.(lrel_md.eq.2).AND.(ewald_mode.eq.1)) then
    write(ilog,*) 'Performance (not correctness) tests of linked FFTW library'
    write(ilog,*) 'From outside parallel code:'
!   test fftw plan
    do lo=1,thrdat%maxn
!
      call pme_plan_threads(lo)
!
      do i=1,kdims(1,2)
        do j=1,kdims(2,2)
          do k=1,kdims(3,2)
            Qew1(i,j,k) = (random()-0.5)*0.2
          end do
        end do
      end do
!
      call System_Clock(count=ttt(lo,1))      
      call fftw_execute_dft(fftplanb,Qew1,Qew2)
      call System_Clock(count=ttt(lo,2))
      
    end do
!
 66 format(10000(g12.5,1x))
    write(ilog,*) 'Time (ms) for maximum number of 1 ... ',thrdat%maxn,' threads used by FFTW: '
    write(ilog,66) 1000.0*(ttt(1:thrdat%maxn,2)-ttt(1:thrdat%maxn,1))/thrdat%rate 
    write(ilog,*) 
    write(ilog,*) 'From within parallel (nested) code:'
    do lo=1,thrdat%maxn
!
      call pme_plan_threads(lo)
!
      do i=1,kdims(1,2)
        do j=1,kdims(2,2)
          do k=1,kdims(3,2)
            Qew1(i,j,k) = (random()-0.5)*0.2
          end do
        end do
      end do
!
      call OMP_SET_NUM_THREADS(lo)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi)
      call OMP_SET_NESTED(.true.)
      tpi = OMP_GET_THREAD_NUM()+1
!$OMP SINGLE
      call System_Clock(count=ttt(lo,3))      
      call fftw_execute_dft(fftplanb,Qew1,Qew2)
      call System_Clock(count=ttt(lo,4))
!$OMP END SINGLE
!$OMP END PARALLEL      
    end do
    write(ilog,*) 'Time (ms) for maximum number of 1 ... ',thrdat%maxn,' threads used by FFTW: '
    write(ilog,66) 1000.0*(ttt(1:thrdat%maxn,4)-ttt(1:thrdat%maxn,3))/thrdat%rate 
!
    Qew1(:,:,:) = 0.0
    write(ilog,*) 
  end if
#endif
!
  write(ilog,*) '*** END OF TEST OF MULTI-THREADED ROUTINES ***'
  write(ilog,*) 
  write(ilog,*) 'CAMPARI will now terminate with an error. This is irrespective of the test results.'
  write(ilog,*) 
  call fexit()
!
end
!
#else
!
subroutine threads_sanity_checks()
!
  use iounit
!
  write(ilog,*) 'Fatal. Called threads_sanity_checks without having compiled support for threads. This &
 &is a bug.'
  call fexit()
!
end
!
#endif
!
!--------------------------------------------------------------------------------------------
!
