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
! CONTRIBUTIONS: Hoang Tran, Rohit Pappu, Adam Steffen                     !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!----------------------------------------------------------------------
!
! FIRST: wrapper routines for mcmove()
!
!----------------------------------------------------------------------
!
subroutine mcmove_chi(rs,mode,tpi)
!
  use iounit
  use cutoffs
  use sequen
  use energies
#ifdef ENABLE_THREADS
  use threads, ONLY: thr_mlgix
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi
!
  integer imol
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL,jfl(4)
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, mcmove_chi(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
#endif  
!
  imol = molofrs(rs)
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call makeref_poly(imol)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call makeref_forsc(rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
!$OMP SINGLE
#endif
    call sample_sc(rs,mode)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
!$OMP SINGLE
#endif
    call makexyz_forsc(rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call updateresgp(rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    end if
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
      if ((tpi.gt.0).AND.(thr_mlgix(imol).gt.0)) then
        jfl(1) = .true.
        jfl(2:4) = .false.
        jfl(3) = .true.
        call molops_threads_geo(thr_mlgix(imol),tpi,jfl)
      else
!$OMP SINGLE
#endif
      call update_rigid(imol)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end if
  else if (mode.eq.2) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call restore_forsc(rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
!$OMP SINGLE
#endif
    call getref_forsc(rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call updateresgp(rs)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_chi(...) with unknown mode.&
 & Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_genpiv(rs,mode,ct,movty,tpi)
!
  use iounit
  use cutoffs
  use molecule
  use energies
  use sequen
  use atoms
  use fyoc, ONLY: nchi
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat,thr_mlgix
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi,movty
  logical, INTENT(IN):: ct
!
  integer j,imol,azero,aone
#ifdef ENABLE_THREADS
  integer stx(2),tpn,sta,incr
  logical OMP_IN_PARALLEL,jfl(4)
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, mcmove_genpiv(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  tpn = thrdat%maxn
  if (tpi.gt.0) then
    sta = tpi-1
    incr = tpn
  else
    sta = 0
    incr = 1
  end if
#endif  
!
  if ((movty.le.0).OR.(movty.gt.4)) then
    write(ilog,*) 'Fatal. Called mcmove_genpiv(...) with unsupported or illegal move type. This is a bug.'
    call fexit()
  end if
  aone = 1
  azero = 0
  imol = molofrs(rs)
  if ((mode.eq.0).OR.(mode.eq.1)) then
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call makeref_poly(imol)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if ((movty.eq.1).OR.(movty.eq.4)) then
      call sample_bb(rs,mode)  ! only Z matrix
      if ((movty.eq.4).AND.(nchi(rs).gt.0)) call sample_sc(rs,mode)
    else if (movty.eq.2) then
      call sample_w(rs,mode)   ! only Z matrix
    else if (movty.eq.3) then
      call sample_nuc(rs,mode) ! only Z matrix
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    if (ct.EQV..true.) then
#ifdef ENABLE_THREADS
      call threads_bounds(atmol(imol,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        xref(stx(1):stx(2)) = x(stx(1):stx(2))
        yref(stx(1):stx(2)) = y(stx(1):stx(2))
        zref(stx(1):stx(2)) = z(stx(1):stx(2))
      end if
!$OMP BARRIER
#else
      call makeref_forbb(imol,rsmol(imol,1))
#endif
      if ((movty.eq.1).OR.(movty.eq.4)) then
        call quatrot_pivot(aone,rs,tpi)
        call quatrot_pivot(azero,rs,tpi)
      else if (movty.eq.2) then
        call quatrot_omega(aone,rs,tpi)
        call quatrot_omega(azero,rs,tpi)
      else if (movty.eq.3) then
        call quatrot_nuc(aone,rs,tpi)
        call quatrot_nuc(azero,rs,tpi)
      end if
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
        do j=rsmol(imol,1)+sta,rs,incr
#else
        do j=rsmol(imol,1),rs ! rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    else
#ifdef ENABLE_THREADS
      call threads_bounds(rsinfo(rs,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        xref(stx(1):stx(2)) = x(stx(1):stx(2))
        yref(stx(1):stx(2)) = y(stx(1):stx(2))
        zref(stx(1):stx(2)) = z(stx(1):stx(2))
      end if
!$OMP BARRIER
#else
      call makeref_forbb(imol,rs)
#endif
      if ((movty.eq.1).OR.(movty.eq.4)) then
        call quatrot_pivot(aone,rs,tpi)
      else if (movty.eq.2) then
        call quatrot_omega(aone,rs,tpi)
      else if (movty.eq.3) then
        call quatrot_nuc(aone,rs,tpi)
      end if
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
        do j=rs+sta,rsmol(imol,2),incr
#else
        do j=rs,rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    end if
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
      if ((tpi.gt.0).AND.(thr_mlgix(imol).gt.0)) then
        jfl(1) = .true.
        jfl(2:4) = .false.
        jfl(3) = .true.
        call molops_threads_geo(thr_mlgix(imol),tpi,jfl)
      else
!$OMP SINGLE
#endif
      call update_rigid(imol)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end if
  else if (mode.eq.2) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if ((movty.eq.1).OR.(movty.eq.4)) then
      call restore_forfy(rs) ! only Z matrix
      if ((movty.eq.4).AND.(nchi(rs).gt.0)) call restore_forsc(rs) ! only Z matrix
    else if (movty.eq.2) then
      call restore_forw(rs) ! only Z matrix
    else if (movty.eq.3) then
      call restore_fornuc(rs) ! only Z matrix
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    if (ct.EQV..true.) then
#ifdef ENABLE_THREADS
      call threads_bounds(atmol(imol,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        x(stx(1):stx(2)) = xref(stx(1):stx(2))
        y(stx(1):stx(2)) = yref(stx(1):stx(2))
        z(stx(1):stx(2)) = zref(stx(1):stx(2))
      end if
#else
      call getref_forbb(imol,rsmol(imol,1))
#endif
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
        do j=rsmol(imol,1)+sta,rs,incr
#else
        do j=rsmol(imol,1),rs ! rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    else
#ifdef ENABLE_THREADS
      call threads_bounds(rsinfo(rs,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        x(stx(1):stx(2)) = xref(stx(1):stx(2))
        y(stx(1):stx(2)) = yref(stx(1):stx(2))
        z(stx(1):stx(2)) = zref(stx(1):stx(2))
      end if
#else
      call getref_forbb(imol,rs)
#endif
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
        do j=rs+sta,rsmol(imol,2),incr
#else
        do j=rs,rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_genpiv(...) with unknown mode. Offending mode # is ',mode,'.'
    call fexit()
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_genpuck(rs,puckdofs,mode,ct,movty,tpi)
!
  use iounit
  use cutoffs
  use molecule
  use energies
  use sequen
  use atoms
  use zmatrix, ONLY: ztor
  use fyoc, ONLY: cur_pucks,phi,pucline
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat,thr_mlgix
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs,mode,tpi,movty
  logical, INTENT(IN):: ct
!
  RTYPE puckdofs(7)
  integer j,imol,azero,aone
#ifdef ENABLE_THREADS
  integer stx(2),tpn,sta,incr
  logical OMP_IN_PARALLEL,jfl(4)
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, mcmove_genpuck(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  tpn = thrdat%maxn
  if (tpi.gt.0) then
    sta = tpi-1
    incr = tpn
  else
    sta = 0
    incr = 1
  end if
#endif  
!
  if ((movty.le.0).OR.(movty.gt.3)) then
    write(ilog,*) 'Fatal. Called mcmove_genpuck(...) with unsupported or illegal move type. This is a bug.'
    call fexit()
  end if
  aone = 1
  azero = 0
  imol = molofrs(rs)
  if ((mode.eq.0).OR.(mode.eq.1)) then
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      call makeref_poly(imol)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
    if (ct.EQV..true.) then
#ifdef ENABLE_THREADS
      call threads_bounds(atmol(imol,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        xref(stx(1):stx(2)) = x(stx(1):stx(2))
        yref(stx(1):stx(2)) = y(stx(1):stx(2))
        zref(stx(1):stx(2)) = z(stx(1):stx(2))
      end if
#else
      call makeref_forbb(imol,rsmol(imol,1))
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
      if (movty.eq.1) then
        call sample_pucker(rs,puckdofs,cur_pucks,mode)  ! rebuilds some coordinates as well
      else if (movty.eq.2) then
        call sample_sugar(rs,puckdofs,mode)   ! rebuilds some coordinates as well
      end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      if (movty.eq.1) then
        call quatrot_pucker(aone,rs,tpi)
        call quatrot_pucker(azero,rs,tpi)
      else if (movty.eq.2) then
        call quatrot_sugar(aone,rs,tpi)
        call quatrot_sugar(azero,rs,tpi)
      end if
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
        do j=rsmol(imol,1)+sta,rs,incr
#else
        do j=rsmol(imol,1),rs ! rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    else
#ifdef ENABLE_THREADS
      call threads_bounds(rsinfo(rs,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        xref(stx(1):stx(2)) = x(stx(1):stx(2))
        yref(stx(1):stx(2)) = y(stx(1):stx(2))
        zref(stx(1):stx(2)) = z(stx(1):stx(2))
      end if
#else
      call makeref_forbb(imol,rs)
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
      if (movty.eq.1) then
        call sample_pucker(rs,puckdofs,cur_pucks,mode)  ! rebuilds some coordinates as well
      else if (movty.eq.2) then
        call sample_sugar(rs,puckdofs,mode)   ! rebuilds some coordinates as well
      end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      if (movty.eq.1) then
        call quatrot_pucker(aone,rs,tpi)
      else if (movty.eq.2) then
        call quatrot_sugar(aone,rs,tpi)
      end if
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
        do j=rs+sta,rsmol(imol,2),incr
#else
        do j=rs,rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    end if
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
      if ((tpi.gt.0).AND.(thr_mlgix(imol).gt.0)) then
        jfl(1) = .true.
        jfl(2:4) = .false.
        jfl(3) = .true.
        call molops_threads_geo(thr_mlgix(imol),tpi,jfl)
      else
!$OMP SINGLE
#endif
      call update_rigid(imol)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end if
  else if (mode.eq.2) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call restore_forring(rs) ! only Z matrix
    if (movty.eq.1) phi(rs) = ztor(pucline(rs))
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    if (ct.EQV..true.) then
#ifdef ENABLE_THREADS
      call threads_bounds(atmol(imol,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        x(stx(1):stx(2)) = xref(stx(1):stx(2))
        y(stx(1):stx(2)) = yref(stx(1):stx(2))
        z(stx(1):stx(2)) = zref(stx(1):stx(2))
      end if
#else
      call getref_forbb(imol,rsmol(imol,1))
#endif
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
        do j=rsmol(imol,1)+sta,rs,incr
#else
        do j=rsmol(imol,1),rs ! rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    else
#ifdef ENABLE_THREADS
      call threads_bounds(rsinfo(rs,1),atmol(imol,2),tpi,tpn,stx)
      if (stx(2).ge.stx(1)) then
        x(stx(1):stx(2)) = xref(stx(1):stx(2))
        y(stx(1):stx(2)) = yref(stx(1):stx(2))
        z(stx(1):stx(2)) = zref(stx(1):stx(2))
      end if
#else
      call getref_forbb(imol,rs)
#endif
      if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
        do j=rs+sta,rsmol(imol,2),incr
#else
        do j=rs,rsmol(imol,2)
#endif
          call updateresgp(j)
        end do
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_genpuck(...) with unknown mode. Offending mode # is ',mode,'.'
    call fexit()
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!----------------------------------------------------------------------
!
! there is a complication here with C-terminal alignment:
! the ct flag is meant to indicate that C-termini align, but this may or may not
! already be reflected in the rotation lists -> check and pass appropriately to getref/setref
!
subroutine mcmove_other(ati,mode,ct,setalC,isnat,tpi)
!
  use iounit
  use sequen
  use atoms
  use zmatrix
  use energies
  use cutoffs
  use math
  use movesets, ONLY: mvcur
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat,thr_mlgix
#endif
!
  implicit none
!
  integer, INTENT(IN):: mode,ati,tpi
  logical, INTENT(IN):: ct,isnat,setalC
!
  integer ix,rs
#ifdef ENABLE_THREADS
  integer sta,incr
  logical OMP_IN_PARALLEL,jfl(4)
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, mcmove_genpuck(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  if (tpi.gt.0) then
    sta = tpi-1
    incr = thrdat%maxn
  else
    sta = 0
    incr = 1
  end if
#endif 
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
  ix = 1
  if (ct.EQV..true.) ix = 2
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    call makeref_forrotlst(ati,setalC,tpi)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
      call makeref_poly(molofrs(atmres(ati)))
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call sample_other(ati,mode,mvcur%hlp(1),isnat)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    call quatxyz_forrotlst(ati,ct,mvcur%hlp(1),tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
      do rs=izrot(ati)%rsbnds(1,ix)+sta,izrot(ati)%rsbnds(2,ix),incr
#else
      do rs=izrot(ati)%rsbnds(1,ix),izrot(ati)%rsbnds(2,ix)
#endif
        call updateresgp(rs)
      end do
    end if
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
#ifdef ENABLE_THREADS
      if ((tpi.gt.0).AND.(thr_mlgix(molofrs(atmres(ati))).gt.0)) then
        jfl(1) = .true.
        jfl(2:4) = .false.
        jfl(3) = .true.
        call molops_threads_geo(thr_mlgix(molofrs(atmres(ati))),tpi,jfl)
      else
!$OMP SINGLE
#endif
      call update_rigid(molofrs(atmres(ati)))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end if
  else if (mode.eq.2) then
    call getref_forrotlst(ati,setalC,tpi)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    call restore_forother(ati,isnat)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
      do rs=izrot(ati)%rsbnds(1,ix)+sta,izrot(ati)%rsbnds(2,ix),incr
#else
      do rs=izrot(ati)%rsbnds(1,ix),izrot(ati)%rsbnds(2,ix)
#endif
        call updateresgp(rs)
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_other(...) with unknown mode&
 &. Offending mode # is ',mode,'.'
    call fexit()
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_cr(rsi,rsf,dfy,mode,ct)
!
  use iounit
  use molecule
  use cutoffs
  use sequen
  use movesets
  use energies
!
  implicit none
!
  integer rsi,rsf,i,mode,j
  logical ct
  RTYPE dfy(MAXCRDOF)
!
  if ((mode.eq.0).OR.(mode.eq.1)) then
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
      call makeref_poly(molofrs(rsi))
    end if
    if (ct.EQV..true.) then
      call makeref_forbb(molofrs(rsi),rsmol(molofrs(rsi),1))
      call sample_bb_cr(rsi,rsf,dfy)
      call alignC_CR(rsi,rsf,dfy)
      if (use_mcgrid.EQV..true.) then
        do j=rsmol(molofrs(rsi),1),rsmol(molofrs(rsi),2)
          call updateresgp(j)
        end do
      end if
    else
      call makeref_forbb(molofrs(rsi),rsi)
      call sample_bb_cr(rsi,rsf,dfy)
      if (use_mcgrid.EQV..true.) then
        do j=rsi,rsmol(molofrs(rsi),2)
          call updateresgp(j)
        end do
      end if
    end if
    if ((mode.eq.0).AND.(use_POLY.EQV..true.)) then
      call update_rigid(molofrs(rsi))
    end if
  else if (mode.eq.2) then
    do i=rsi,rsf-1
      call restore_forfy(i)
    end do
    if (ct.EQV..true.) then
      call getref_forbb(molofrs(rsi),rsmol(molofrs(rsi),1))
      if (use_mcgrid.EQV..true.) then
        do j=rsmol(molofrs(rsi),1),rsmol(molofrs(rsi),2)
          call updateresgp(j)
        end do
      end if
    else
      call getref_forbb(molofrs(rsi),rsi)
      if (use_mcgrid.EQV..true.) then
        do j=rsi,rsmol(molofrs(rsi),2)
          call updateresgp(j)
        end do
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_cr(...) with unknown mode. O&
 &ffending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! for the following three:
!
! mode 0: sample prerotation segment (+ swivelling end), do necessary backup work
!         dfy are determined on the outside
! mode 1: sample (or re-sample) everything including chain closure
! mode 2: restore to original state
! mode 3: restore just the prerotation without worrying about energy-associated things just yet
! mode 4: restore Cartesian coordinates of post-prerotation segment (this is tricky due to residue
!         being broken up)
! mode 5: restore to original state but do not worry about energetics
!
!
!-----------------------------------------------------------------------
!
subroutine mcmove_nuccr(rsi,rsf,dof,dfy,mode)
!
  use iounit
  use molecule
  use cutoffs
  use sequen
  use movesets
  use energies
  use fyoc
  use zmatrix
  use ujglobals
!
  implicit none
!
  integer rsi,rsf,mode,j,azero,aone,dof
  RTYPE dfy(MAXUJDOF+6)
  azero = 0
  aone = 1
!
  if (mode.eq.0) then
    call makeref_forbb(molofrs(rsi),rsi)
    call sample_bb_nuccrpr(rsi,rsf,dof,dfy,azero)
  else if (mode.eq.1) then
    call sample_bb_nuccrfull(rsi,rsf,dof,dfy,azero)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    do j=rsi,rsf
      call restore_fornuc(j)
    end do
    call getref_forbb_nuccr(molofrs(rsi),rsi)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.3) then
    do j=rsi,rsf-1 ! has to exclude rsf since rsf is never sampled via sample_bb_nuccrpr
      call restore_fornuc(j)
    end do
    call getref_forbb_nuccr(molofrs(rsi),rsi)
  else if (mode.eq.4) then
    call getref_forbb_nuccr_pre(molofrs(rsf),rsf-1)
  else if (mode.eq.5) then
    do j=rsi,rsf
      call restore_fornuc(j)
    end do
    call getref_forbb_nuccr(molofrs(rsi),rsi)
  else
    write(ilog,*) 'Fatal. Called mcmove_nuccr(...) with unknown mode.&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine mcmove_torcr_dj(rsi,rsf,dof,dfy,mode)
!
  use iounit
  use molecule
  use cutoffs
  use sequen
  use movesets
  use energies
  use fyoc
  use zmatrix
  use ujglobals
!
  implicit none
!
  integer rsi,rsf,mode,j,dof
  RTYPE dfy(MAXUJDOF+6)
!
  if (mode.eq.0) then
    call makeref_forbb(molofrs(rsi),rsi)
    call sample_bb_torcrpr2(rsi,rsf,dof,dfy)
  else if (mode.eq.1) then
    call sample_bb_torcrfull2(rsi,rsf,dof,dfy)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    do j=rsi,rsf
      call restore_forfy(j)
      if (j.le.rsf-2) then
        call restore_forw(j)
      end if
    end do
    call getref_forbb(molofrs(rsi),rsi)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.3) then
    call restore_forw(rsf-2)
    do j=rsi,rsf-3
      call restore_forfy(j)
      call restore_forw(j)
    end do
    call getref_forbb(molofrs(rsi),rsi)
  else if (mode.eq.4) then
    call getref_fordjo(molofrs(rsf),rsf-2)
  else if (mode.eq.5) then
    do j=rsi,rsf
      call restore_forfy(j)
      if (j.le.rsf-2) then
       call restore_forw(j)
      end if
    end do
    call getref_forbb(molofrs(rsi),rsi)
  else
    write(ilog,*) 'Fatal. Called mcmove_torcr_dj(...) with unknown mode.&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine mcmove_torcr_aux(rsi,rsf,which,curpks,mode,mode2)
!
  use sequen
  use movesets
  use ujglobals
  use iounit
!
  implicit none
!
  integer j,rsf,which,mode,k,atwo,mode2,rsoff,rsi
  RTYPE curpks(MAXSOLU*4,7,3),curpt(25),pdofs(7)
!
  atwo = 2
  if (mode2.eq.1) then
    rsoff = 2
  else if (mode2.eq.2) then
    rsoff = 1
  else
    write(ilog,*) 'Fatal. Called mcmove_torcr_aux(...) with unknown mode2.&
 &Offending mode2 # is ',mode2,'.'
    call fexit()
  end if
!
  if (mode.eq.1) then
    k = 0
    do j=rsf-rsoff,rsf
      if (seqflag(j).eq.5) then
        k = k + 1
        curpt(1:7) = curpks(which,1:7,k)
        call sample_pucker(j,pdofs,curpt,atwo)
      end if
    end do
  else if (mode.eq.2) then
    do j=rsf-rsoff,rsf
      if (seqflag(j).eq.5) then
        call restore_forring(j)
      end if
    end do
  else
    write(ilog,*) 'Fatal. Called mcmove_torcr_aux(...) with unknown mode.&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine mcmove_torcr_do(rsi,rsf,dof,dfy,mode)
!
  use iounit
  use molecule
  use cutoffs
  use sequen
  use movesets
  use energies
  use fyoc
  use zmatrix
  use ujglobals
!
  implicit none
!
  integer rsi,rsf,mode,j,dof
  RTYPE dfy(MAXUJDOF+6)
!
  if (mode.eq.0) then
    call makeref_forbb(molofrs(rsi),rsi)
    call sample_bb_torcrpr(rsi,rsf,dof,dfy)
  else if (mode.eq.1) then
    call sample_bb_torcrfull(rsi,rsf,dof,dfy)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    do j=rsi,rsf
      call restore_forfy(j)
      call restore_forw(j)
    end do
    call getref_forbb(molofrs(rsi),rsi)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.3) then
    do j=rsi,rsf-2
      call restore_forfy(j)
      call restore_forw(j)
    end do
    call getref_forbb(molofrs(rsi),rsi)
  else if (mode.eq.4) then
    call getref_fordo(molofrs(rsf),rsf-1)
  else if (mode.eq.5) then
    do j=rsi,rsf
      call restore_forfy(j)
      call restore_forw(j)
    end do
    call getref_forbb(molofrs(rsi),rsi)
  else
    write(ilog,*) 'Fatal. Called mcmove_torcr_do(...) with unknown mode.&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!
!-----------------------------------------------------------------------
!
! mode 0: sample prerotation segment (+ swivelling end), do necessary backup work
!         dfy are determined on the outside
! mode 1: sample (or re-sample) everything including chain closure
! mode 2: restore to original state
! mode 3: restore just the prerotation without worrying about energy-associated things just yet
! mode 4: restore Cartesian coordinates of post-prerotation segment (this is tricky due to residue
!         being broken up)
!         
!
subroutine mcmove_uj(rsi,rsf,dfy,mode)
!
  use iounit
  use molecule
  use cutoffs
  use sequen
  use movesets
  use energies
  use fyoc
  use zmatrix
  use ujglobals
!
  implicit none
!
  integer rsi,rsf,mode,j
  RTYPE dfy(MAXUJDOF+6)
!
  if (mode.eq.0) then
    if (use_POLY.EQV..true.) then
      call makeref_poly(molofrs(rsi))
    end if
    call makeref_forbb(molofrs(rsi),rsi)
    call sample_bb_ujpr(rsi,rsf-2,dfy)
  else if (mode.eq.1) then
    call sample_bb_ujcc(rsf,dfy)
    call sample_bb_ujpr(rsi,rsf-2,dfy)
    if (use_POLY.EQV..true.) then
      call update_rigid(molofrs(rsi))
    end if
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    do j=rsi,rsf-1
      call restore_forfy(j)
      call restore_forabg(j)
    end do
    call restore_forfy(rsf)
    call getref_forbb(molofrs(rsi),rsi)
    if (use_mcgrid.EQV..true.) then
      do j=rsi,rsf
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.3) then
    do j=rsi,rsf-2
      call restore_forfy(j)
      call restore_forabg(j)
    end do
    call getref_forbb(molofrs(rsi),rsi)
  else if (mode.eq.4) then
!    do j=rsf-1,rsmol(molofrs(rsf),2)
!      call restore_forfy(j)
!      call restore_forabg(j)
!    end do
    call getref_foruj(molofrs(rsf),rsf-1)
  else
    write(ilog,*) 'Fatal. Called mcmove_fyc(...) with unknown mode.&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------
!
subroutine mcmove_lct(ff,mode)
!
  use iounit
  use cutoffs
  use sequen
!
  implicit none
!
  integer ff,mode,j,aone
!
  aone = 1
!
  if (mode.eq.0) then
!    WARNING: totally broken!!!!
    call makeref_forbb(aone,aone)
    call sample_lct(ff)
    if (use_mcgrid.EQV..true.) then
      do j=1,nseq
        call updateresgp(j)
      end do
    end if
  else if (mode.eq.2) then
    call restore_forlct(ff)
!    WARNING: totally broken!!!!
    call getref_forbb(aone,aone)
    if (use_mcgrid.EQV..true.) then
      do j=1,nseq
        call updateresgp(j)
      end do
    end if
  else
    write(ilog,*) 'Fatal. Called mcmove_lct(...) with unknown mode. &
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
! wrapper function for setting (transferring) increments globally in torsional dynamics
! mode 0: override ztorpr
! mode 1: do not override ztorpr
!
subroutine tmdmove(imol,mode)
!
  use fyoc
  use molecule
  use forces
  use zmatrix
!
  implicit none
!
  integer, INTENT(IN):: imol,mode
!
  integer rs,ttc,azero,nt,i
  logical afalse
  RTYPE local_other,bufy(4),local_chis(MAXCHI),local_omega,local_nucs(6),local_phi,local_psi
  logical local_nucflag(6),local_chiflag(MAXCHI)
!
  afalse = .false.
  azero = 0
!
  do rs=rsmol(imol,1),rsmol(imol,2)
    if (wnr(rs).gt.0) then
      ttc = wnr(rs)
      local_omega = omega(rs) + dc_di(imol)%incr(ttc)
      if (local_omega.gt.180.0) local_omega = local_omega - 360.0
      if (local_omega.lt.-180.0) local_omega = local_omega + 360.0
      call setw(rs,local_omega,mode)
    end if
    local_phi = phi(rs) ! may be nonsensical
    local_psi = psi(rs) ! may be nonsensical
    if (fnr(rs).gt.0) then
      ttc = fnr(rs)
      local_phi = phi(rs) + dc_di(imol)%incr(ttc)
      if (local_phi.gt.180.0) local_phi = local_phi - 360.0
      if (local_phi.lt.-180.0) local_phi = local_phi + 360.0
    end if
    if (ynr(rs).gt.0) then
      ttc = ynr(rs)
      local_psi = psi(rs) + dc_di(imol)%incr(ttc)
      if (local_psi.gt.180.0) local_psi = local_psi - 360.0
      if (local_psi.lt.-180.0) local_psi = local_psi + 360.0
    end if
    if ((fnr(rs).gt.0).OR.(ynr(rs).gt.0)) then
      if (mode.eq.1) then
        if (fline(rs).gt.0)  bufy(1) = ztorpr(fline(rs))
        if (fline2(rs).gt.0) bufy(2) = ztorpr(fline2(rs))
        if (yline(rs).gt.0)  bufy(3) = ztorpr(yline(rs))
        if (yline2(rs).gt.0) bufy(4) = ztorpr(yline2(rs))
      end if
      call setfy(rs,local_phi,local_psi,azero)
      if (mode.eq.1) then
        if (fline(rs).gt.0)  ztorpr(fline(rs)) = bufy(1)
        if (fline2(rs).gt.0) ztorpr(fline2(rs)) = bufy(2)
        if (yline(rs).gt.0)  ztorpr(yline(rs)) = bufy(3)
        if (yline2(rs).gt.0) ztorpr(yline2(rs)) = bufy(4)
      end if
    end if
    do nt=1,nnucs(rs)
      if (nucsnr(nt,rs).le.0) then
        local_nucs(nt) = nucs(nt,rs)
        cycle
      end if
      ttc = nucsnr(nt,rs)
      local_nucs(nt) = nucs(nt,rs) + dc_di(imol)%incr(ttc)
      if (local_nucs(nt).gt.180.0) local_nucs(nt) = local_nucs(nt) - 360.0
      if (local_nucs(nt).lt.-180.0) local_nucs(nt) = local_nucs(nt) + 360.0
      local_nucflag(nt) = .true.
    end do
    if (nnucs(rs).gt.0) then
      call setnucs(rs,local_nucs,mode)
    end if
    do nt=1,nchi(rs)
      if (chinr(nt,rs).le.0) then
        local_chis(nt) = chi(nt,rs)
        cycle
      end if
      ttc = chinr(nt,rs)
      local_chis(nt) = chi(nt,rs) + dc_di(imol)%incr(ttc)
      if (local_chis(nt).gt.180.0) local_chis(nt) = local_chis(nt) - 360.0
      if (local_chis(nt).lt.-180.0) local_chis(nt) = local_chis(nt) + 360.0
      local_chiflag(nt) = .true.
    end do
    if (nchi(rs).gt.0) then
      call setchi(rs,local_chis,mode)
    end if
  end do
  if (othidxmol(moltypid(imol)).gt.0) then
    do i=1,dc_di(imol)%maxntor
      ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
      if (ttc.lt.othidxmol(moltypid(imol))) cycle
      local_other = ztor(dc_di(imol)%recurs(i,1)) + dc_di(imol)%incr(ttc)
      if (local_other.gt.180.0) local_other = local_other - 360.0
      if (local_other.lt.-180.0) local_other = local_other + 360.0
      call setother(dc_di(imol)%recurs(i,1),local_other,mode,afalse)
    end do
  end if
!
end
!
!----------------------------------------------------------------------
!
! SECOND: sampling routines for the different movetypes
!
!----------------------------------------------------------------------
!
! sample_sc just varies a single sidechain chis by completely
! random re-assignment
!
subroutine sample_sc(rs,mode)
!
  use iounit
  use sequen
  use fyoc
  use movesets
!
  implicit none
!
  integer rs,mode,i
  RTYPE random,cc(MAXCHI)
!
! mode=0: sample fresh, mode=1: apply current sampling move
  if (mode.eq.0) then
!   proline is special
!    if (seqtyp(rs).eq.9) then
!!     we're going to temporarily abuse cur_chis(4) and cc(4) to store
!!     the pucker indicator
!      if (random().lt.0.5) then
!        cc(4) = 0.0
!      else
!        cc(4) = 1.0
!      end if
!      cur_chis(4) = cc(4)
!      call pucker(rs,cc(4),mode)
!    else
    if (random().ge.chi_randfreq) then
      do i=1,nchi(rs)
        if (cur_chiflag(i).EQV..true.) then
          cc(i) = chi(i,rs) + (random()-0.5)*chi_stepsz
          if (cc(i).gt.180.0) cc(i) = cc(i) - 360.0
          if (cc(i).lt.-180.0) cc(i) = cc(i) + 360.0
          cur_chis(i) = cc(i)
        else
          cc(i) = chi(i,rs)
          cur_chis(i) = chi(i,rs)
        end if
      end do
      call setchi(rs,cc,mode)
    else
      do i=1,nchi(rs)
        if (cur_chiflag(i).EQV..true.) then
          cc(i) = random()*360.0 - 180.0
          cur_chis(i) = cc(i)
        else
          cc(i) = chi(i,rs)
          cur_chis(i) = chi(i,rs)
        end if
      end do
      call setchi(rs,cc,mode)
    end if
!      call setchi(seqtyp(rs),mode,rs)
!    end if
  else if (mode.eq.1) then
!    if (seqtyp(rs).eq.9) then
!      call pucker(rs,cur_chis(4),mode)
!    else
    call setchi(rs,cur_chis,mode)
!    end if
!    call setchi2(seqtyp(rs),rs,cur_chis)
  else
    write(ilog,*) 'Fatal. Called sample_sc(...) with unknown mode.&
 & Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------
!
! sample_sugar creates a new coupled backbone+sidechain state along the C4*-C3* in ribo-
! and deoxyriboses, i.e., for the sampling of polynucleotides
!
subroutine sample_sugar(rs,puckdofs,mode)
!
  use iounit
  use sequen
  use fyoc
  use zmatrix
  use polypep
  use movesets
  use atoms
  use aminos
  use molecule
  use math
  use ujglobals
  use cutoffs, ONLY: rsinfo
!
  implicit none
!
  integer rs,mode,aone,atwo,chosenmode
  RTYPE random,refvals(2),puckdofs(7)
  logical findsolu
!
  aone = 1
  atwo = 2
!
! mode=0: sample fresh, mode=1: apply current sampling move
  if (mode.eq.0) then
    call pucker_refvals(rs,refvals)
    if (random().lt.nucpuckrdfreq) then
      call freepucker(rs,puckdofs,aone)
      chosenmode = 1
    else
      call freepucker(rs,puckdofs,atwo)
      chosenmode = 2
    end if
!   backup and modify C5*-C4*-C3*-C2*, check wrap-around
    ztorpr(at(rs)%sc(1)) = ztor(at(rs)%sc(1))
    ztor(at(rs)%sc(1)) = ztor(at(rs)%sc(1)) + puckdofs(2)
    if (ztor(at(rs)%sc(1)).gt.180.0) ztor(at(rs)%sc(1)) = ztor(at(rs)%sc(1)) - 360.0
    if (ztor(at(rs)%sc(1)).lt.-180.0) ztor(at(rs)%sc(1)) = ztor(at(rs)%sc(1)) + 360.0
    cur_pucks(1) = ztor(at(rs)%sc(1))
!   backup and modify C4*-C3*-C2*-C1*
    ztorpr(at(rs)%sc(2)) = ztor(at(rs)%sc(2))
    ztor(at(rs)%sc(2)) = ztor(at(rs)%sc(2)) + puckdofs(1)
    if (ztor(at(rs)%sc(2)).gt.180.0) ztor(at(rs)%sc(2)) = ztor(at(rs)%sc(2)) - 360.0
    if (ztor(at(rs)%sc(2)).lt.-180.0) ztor(at(rs)%sc(2)) = ztor(at(rs)%sc(2)) + 360.0
    cur_pucks(2) = ztor(at(rs)%sc(2))
!   backup and modify C2*-C3*-C4*-O4*
    ztorpr(at(rs)%sc(3)) = ztor(at(rs)%sc(3))
    if (chosenmode.eq.2) then
      ztor(at(rs)%sc(3)) = ztor(at(rs)%sc(3)) + puckdofs(4)
    else
      ztor(at(rs)%sc(3)) = ztor(at(rs)%sc(3)) + puckdofs(2)
    end if
    if (ztor(at(rs)%sc(3)).gt.180.0) ztor(at(rs)%sc(3)) = ztor(at(rs)%sc(3)) - 360.0
    if (ztor(at(rs)%sc(3)).lt.-180.0) ztor(at(rs)%sc(3)) = ztor(at(rs)%sc(3)) + 360.0
    cur_pucks(3) = ztor(at(rs)%sc(3))
!   backup and modify angle C4*-C3*-C2*
    bangpr(at(rs)%sc(1)) = bang(at(rs)%sc(1))
    if (chosenmode.eq.2) bang(at(rs)%sc(1)) = bang(at(rs)%sc(1)) + puckdofs(5)
    cur_pucks(4) = bang(at(rs)%sc(1))
!   backup and modify angle C3*-C2*-C1*
    bangpr(at(rs)%sc(2)) = bang(at(rs)%sc(2))
    if (chosenmode.eq.2) bang(at(rs)%sc(2)) = bang(at(rs)%sc(2)) + puckdofs(6)
    cur_pucks(5) = bang(at(rs)%sc(2))
!   backup angle C3*-C4*-O4*
    bangpr(at(rs)%sc(3)) = bang(at(rs)%sc(3))
!   backup and modify remaining Z-matrix-used torsions around C3*-C4* and check wrap-arounds
    if (rs.eq.rsmol(molofrs(rs),2)) then
      if (moltermid(molofrs(rs),2).ne.1) then
        write(ilog,*) 'Fatal. Encountered unknown terminus t&
 &ype(s) for molecule ',molofrs(rs),' in sample_sugar(...).'
        call fexit()
      end if
    end if
!   backup and modify C5*-C4*-C3*-O3* (3'-terminal or +1)
    ztorpr(nucsline(6,rs)) = ztor(nucsline(6,rs))
    if (chosenmode.eq.2) then
      ztor(nucsline(6,rs)) =  ztor(nucsline(6,rs)) + puckdofs(3)
    else
      ztor(nucsline(6,rs)) =  ztor(nucsline(6,rs)) + puckdofs(2)
    end if
    if (ztor(nucsline(6,rs)).gt.180.0) ztor(nucsline(6,rs)) = ztor(nucsline(6,rs)) - 360.0
    if (ztor(nucsline(6,rs)).lt.-180.0) ztor(nucsline(6,rs)) = ztor(nucsline(6,rs)) + 360.0
    cur_pucks(7) = ztor(nucsline(6,rs))
!   re-build xyz
    call makexyz_forset(rsinfo(rs,1),rsinfo(rs,1)+rsinfo(rs,2))
!   now close the 5-ring (this call populates puckdofs(7))
    if (chosenmode.eq.2) then
      call ringclose(rs,refvals,puckdofs,ztor(at(rs)%sc(3))/RADIAN,findsolu)
    else
      puckdofs(7) = bang(at(rs)%sc(3))
      findsolu = .true.
    end if
    if (findsolu.EQV..true.) then
!     lastly modify angle C3*-C4*-O4*
      if (chosenmode.eq.2) bang(at(rs)%sc(3)) = puckdofs(7)
      cur_pucks(6) = bang(at(rs)%sc(3))
    else
      pc_wrncnt(1) = pc_wrncnt(1) + 1
      if (pc_wrncnt(1).eq.pc_wrnlmt(1)) then
        write(ilog,*) 'Warning. Failed to find a new ring conformation in sample_sugar(...). This &
 &is indicative of too large step sizes or a near-fatal conformation (if no bond angle potentials&
 & are used). This simulation may crash due to bond angle exceptions.'
        write(ilog,*) 'This was warning #',pc_wrncnt(1),' of this type not all of which may be displayed.'
       if (10.0*pc_wrnlmt(1).gt.0.5*HUGE(pc_wrnlmt(1))) then
          pc_wrncnt(1) = 0 ! reset
        else
          pc_wrnlmt(1) = pc_wrnlmt(1)*10
        end if
      end if
 !     nucpuckrdfreq = 1.0
      bang(at(rs)%sc(3)) = puckdofs(7)
      cur_pucks(6) = bangpr(at(rs)%sc(3))
      cur_pucks(1) = ztorpr(at(rs)%sc(1))
      cur_pucks(2) = ztorpr(at(rs)%sc(2))
      cur_pucks(3) = ztorpr(at(rs)%sc(3))
      cur_pucks(4) = bangpr(at(rs)%sc(1))
      cur_pucks(5) = bangpr(at(rs)%sc(2))
      cur_pucks(7) = ztorpr(nucsline(6,rs))
      ztor(at(rs)%sc(1)) = cur_pucks(1)
      ztor(at(rs)%sc(2)) = cur_pucks(2)
      ztor(at(rs)%sc(3)) = cur_pucks(3)
      bang(at(rs)%sc(1)) = cur_pucks(4)
      bang(at(rs)%sc(2)) = cur_pucks(5)
      ztor(nucsline(6,rs)) = cur_pucks(7)
    end if
  else if (mode.eq.1) then
!   re-assign C5*-C4*-C3*-C2*
    ztor(at(rs)%sc(1)) =  cur_pucks(1)
!   re-assign C4*-C3*-C2*-C1*
    ztor(at(rs)%sc(2)) =  cur_pucks(2)
!   re-assign C2*-C3*-C4*-O4*
    ztor(at(rs)%sc(3)) = cur_pucks(3)
!   re-assign angle C4*-C3*-C2*
    bang(at(rs)%sc(1)) = cur_pucks(4)
!   re-assign angle C3*-C2*-C1*
    bang(at(rs)%sc(2)) = cur_pucks(5)
!   re-assign angle C3*-C4*-O4*
    bang(at(rs)%sc(3)) = cur_pucks(6)
!   re-assign remaining Z-matrix-used torsions around C3*-C4*
    if (rs.eq.rsmol(molofrs(rs),2)) then
      if (moltermid(molofrs(rs),2).ne.1) then
        write(ilog,*) 'Fatal. Encountered unknown terminus t&
 &ype(s) for molecule ',molofrs(rs),' in sample_sugar(...).'
        call fexit()
      end if
    end if
!   re-assign C5*-C4*-C3*-O3*(3'-terminal or +1)
    ztor(nucsline(6,rs)) =  cur_pucks(7)
  else
    write(ilog,*) 'Fatal. Called sample_sugar(...) with unknown mode.&
 & Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------
!
! the same for cyclic peptide residues, i.e., PRO, HYP, PCA
!
subroutine sample_pucker(rs,puckdofs,curpks,mode)
!
  use iounit
  use sequen
  use fyoc
  use zmatrix
  use polypep
  use movesets
  use atoms
  use aminos
  use molecule
  use math
  use ujglobals
  use system
  use cutoffs, ONLY: rsinfo
!
  implicit none
!
  integer rs,mode,aone,atwo,shf,chosenmode
  RTYPE random,refvals(2),puckdofs(7),curpks(25)
  logical findsolu
!
  aone = 1
  atwo = 2
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! mode=0: sample fresh, mode=1: apply current sampling move
  if (mode.eq.0) then
    call pucker_refvals(rs,refvals)
    if (random().lt.puckerrdfreq) then
      call freepucker(rs,puckdofs,aone)
      chosenmode = 1
    else
      call freepucker(rs,puckdofs,atwo)
      chosenmode = 2
    end if
!   backup and modify C(-1)-N-CA-C, check wrap-around (this also works N-terminal for C-CA-N-HN1)
    ztorpr(pucline(rs)) = ztor(pucline(rs))
    if (chosenmode.eq.2) then
      ztor(pucline(rs)) = ztor(pucline(rs)) + puckdofs(3)
    else
      ztor(pucline(rs)) = ztor(pucline(rs)) + puckdofs(2)
    end if
    if (ztor(pucline(rs)).gt.180.0) ztor(pucline(rs)) = ztor(pucline(rs)) - 360.0
    if (ztor(pucline(rs)).lt.-180.0) ztor(pucline(rs)) = ztor(pucline(rs)) + 360.0
    curpks(1) = ztor(pucline(rs))
!   backup and modify C(-1)-N-CA-CB, check wrap-around
    ztorpr(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf))
    ztor(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf)) + puckdofs(2)
    if (ztor(at(rs)%sc(2-shf)).gt.180.0) ztor(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf)) - 360.0
    if (ztor(at(rs)%sc(2-shf)).lt.-180.0) ztor(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf)) + 360.0
    curpks(2) = ztor(at(rs)%sc(2-shf))
!   backup and modify N-CA-CB-CG
    ztorpr(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf))
    ztor(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf)) + puckdofs(1)
    if (ztor(at(rs)%sc(3-shf)).gt.180.0) ztor(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf)) - 360.0
    if (ztor(at(rs)%sc(3-shf)).lt.-180.0) ztor(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf)) + 360.0
    curpks(3) = ztor(at(rs)%sc(3-shf))
!   backup and modify CB-CA-N-CD 
    ztorpr(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf))
    if (chosenmode.eq.2) then
      ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) + puckdofs(4)
    else
      ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) + puckdofs(2) 
    end if
    if (ztor(at(rs)%sc(4-shf)).gt.180.0) ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) - 360.0
    if (ztor(at(rs)%sc(4-shf)).lt.-180.0) ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) + 360.0
    curpks(4) = ztor(at(rs)%sc(4-shf))
!   backup and modify angle N-CA-CB
    bangpr(at(rs)%sc(2-shf)) = bang(at(rs)%sc(2-shf))
    if (chosenmode.eq.2) bang(at(rs)%sc(2-shf)) = bang(at(rs)%sc(2-shf)) + puckdofs(5)
    curpks(5) = bang(at(rs)%sc(2-shf))
!   backup and modify angle CA-CB-CG
    bangpr(at(rs)%sc(3-shf)) = bang(at(rs)%sc(3-shf))
    if (chosenmode.eq.2) bang(at(rs)%sc(3-shf)) = bang(at(rs)%sc(3-shf)) + puckdofs(6)
    curpks(6) = bang(at(rs)%sc(3-shf))
!   backup angle CA-N-CD
    bangpr(at(rs)%sc(4-shf)) = bang(at(rs)%sc(4-shf))
!   re-build xyz
    call makexyz_forset(rsinfo(rs,1),rsinfo(rs,1)+rsinfo(rs,2))
!   now close the 5-ring (this call populates puckdofs(7))
    if (chosenmode.eq.2) then
      call ringclose(rs,refvals,puckdofs,ztor(at(rs)%sc(4-shf))/RADIAN,findsolu)
    else
      findsolu = .true.
      puckdofs(7) = bang(at(rs)%sc(4-shf))
    end if
    if (findsolu.EQV..true.) then
!     lastly modify angle CA-N-CD
      if (chosenmode.eq.2) bang(at(rs)%sc(4-shf)) = puckdofs(7)
      curpks(7) = bang(at(rs)%sc(4-shf))
    else
      pc_wrncnt(4) = pc_wrncnt(4) + 1
      if (pc_wrncnt(4).eq.pc_wrnlmt(4)) then
        write(ilog,*) 'Warning. Failed to find a new ring conformation in sample_pucker(...). This &
 &is indicative of too large step sizes or a near-fatal conformation (if no bond angle potentials&
 & are used). This simulation may crash due to bond angle exceptions.'
        write(ilog,*) 'This was warning #',pc_wrncnt(4),' of this type not all of which may be displayed.'
        if (10.0*pc_wrnlmt(4).gt.0.5*HUGE(pc_wrnlmt(4))) then
          pc_wrncnt(4) = 0 ! reset
        else
          pc_wrnlmt(4) = pc_wrnlmt(4)*10
        end if
      end if
      bang(at(rs)%sc(4-shf)) = puckdofs(7)
      curpks(7) = bang(at(rs)%sc(4-shf))
      curpks(1) = ztorpr(pucline(rs))
      curpks(2) = ztorpr(at(rs)%sc(2-shf)) 
      curpks(3) = ztorpr(at(rs)%sc(3-shf)) 
      curpks(4) = ztorpr(at(rs)%sc(4-shf))
      curpks(5) = bangpr(at(rs)%sc(2-shf))
      curpks(6) = bangpr(at(rs)%sc(3-shf)) 
      ztor(pucline(rs)) =  curpks(1)
      ztor(at(rs)%sc(2-shf)) = curpks(2)
      ztor(at(rs)%sc(3-shf)) = curpks(3)
      ztor(at(rs)%sc(4-shf)) = curpks(4)
      bang(at(rs)%sc(2-shf)) = curpks(5)
      bang(at(rs)%sc(3-shf)) = curpks(6)
    end if
!   keep phi pointer current
    phi(rs) = ztor(pucline(rs))
  else if (mode.eq.1) then
!   re-assign C(-1)-N-CA-C
    ztor(pucline(rs)) =  curpks(1)
!   re-assign C(-1)-N-CA-CB
    ztor(at(rs)%sc(2-shf)) = curpks(2)
!   re-assign N-CA-CB-CG
    ztor(at(rs)%sc(3-shf)) = curpks(3)
!   re-assign CA-CB-CG-CD
    ztor(at(rs)%sc(4-shf)) = curpks(4)
!   re-assign angle N-CA-CB
    bang(at(rs)%sc(2-shf)) = curpks(5)
!   re-assign angle CA-CB-CG
    bang(at(rs)%sc(3-shf)) = curpks(6)
!   re-assign angle CA-N-CD
    bang(at(rs)%sc(4-shf)) = curpks(7)
    phi(rs) = ztor(pucline(rs))
  else if (mode.eq.2) then
!   re-assign C(-1)-N-CA-CB
    ztor(at(rs)%sc(2-shf)) = curpks(2)
!   re-assign N-CA-CB-CG
    ztor(at(rs)%sc(3-shf)) = curpks(3)
!   re-assign CA-CB-CG-CD
    ztor(at(rs)%sc(4-shf)) = curpks(4)
!   re-assign angle N-CA-CB
    bang(at(rs)%sc(2-shf)) = curpks(5)
!   re-assign angle CA-CB-CG
    bang(at(rs)%sc(3-shf)) = curpks(6)
!   re-assign angle CA-N-CD
    bang(at(rs)%sc(4-shf)) = curpks(7)
!   re-build xyz
    call makexyz_forsc(rs)
  else
    write(ilog,*) 'Fatal. Called sample_pucker(...) with unknown mode.&
 & Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------
!
! sample_bb just varies a single set of phipsi by completely
! random re-assignment, steric grid-bias, or hc-bias 
!
subroutine sample_bb(rs,mode)
!
  use iounit
  use sequen
  use fyoc
  use grids
  use movesets
!
  implicit none
!
  integer rs,mode,ig,ng
  RTYPE random,yc,fc
!
  if (use_stericgrids.EQV..true.) then
!   mode=0: sample fresh, mode=1: apply current sampling move
    if (mode.eq.0) then
!     first gather new (biased) random variables from grids
      ng = stgr%ngr(seqtyp(rs))
      ig = floor(ng*random()) + 1
      if (seqflag(rs).eq.5) then ! do not sample phi for cyclic residue
        fc = phi(rs)
      else
        fc = stgr%it(seqtyp(rs),ig,1)
        cur_phi = fc
      end if
      yc = stgr%it(seqtyp(rs),ig,2)
      cur_psi = yc
!     then set those new values (mode=1 -> grids) for the chosen phi,psi pair
      call setfy(rs,fc,yc,1)
    else if (mode.eq.1) then
      call setfy(rs,cur_phi,cur_psi,0)
    else
      write(ilog,*) 'Fatal. Called sample_bb(...) with unknown mode.&
 & Offending mode # is ',mode,'.'
      call fexit()
    end if
  else
!   mode=0: sample fresh, mode=1: apply current sampling move
    if (mode.eq.0) then
      if (random().lt.pivot_randfreq) then
        yc = random()*360.0 - 180.0
        cur_psi = yc
        if (seqflag(rs).eq.5) then ! do not sample phi for cyclic residue
          fc = phi(rs)
          cur_phi = fc
        else
          fc = random()*360.0 - 180.0
          cur_phi = fc
        end if
      else
        yc = psi(rs) + (random()-0.5)*pivot_stepsz
        if (yc.gt.180.0) yc = yc - 360.0
        if (yc.lt.-180.0) yc = yc + 360.0
        cur_psi = yc
        if (seqflag(rs).eq.5) then ! do not sample phi for cyclic residue
          fc = phi(rs)
          cur_phi = fc
        else
          fc = phi(rs) + (random()-0.5)*pivot_stepsz
          if (fc.gt.180.0) fc = fc - 360.0
          if (fc.lt.-180.0) fc = fc + 360.0
          cur_phi = fc
        end if
      end if
      call setfy(rs,fc,yc,0)
    else if (mode.eq.1) then
       call setfy(rs,cur_phi,cur_psi,0) 
    else
      write(ilog,*) 'Fatal. Called sample_bb(...) with unknown mode.&
 & Offending mode # is ',mode,'.'
      call fexit()
    end if
  end if
!   
end
!c
!-------------------------------------------------------------------------
!
! sample_w just varies a single omega angle by bounded
! random re-assignment or stepwise perturbation
!
subroutine sample_w(rs,mode)
!
  use iounit
  use sequen
  use fyoc
  use movesets
!
  implicit none
!
  integer rs,mode
  RTYPE random,wc
!
! mode=0: sample fresh, mode=1: apply current sampling move
  if (mode.eq.0) then
    if (random().lt.omega_randfreq) then
      wc = random()*360.0 - 180.0
      cur_omega = wc
!      if (random().gt.0.5) then
!        wc = random()*20.0 - 10.0
!        cur_omega = wc
!      else
!        wc = random()*20.0 - 190.0
!        if (wc.lt.-180.0) wc = wc + 360.0
!        cur_omega = wc
!      end if
    else
      wc = omega(rs) + (random()-0.5)*omega_stepsz
      if (wc.gt.180.0) wc = wc - 360.0
      if (wc.lt.-180.0) wc = wc + 360.0
      cur_omega = wc
    end if
    call setw(rs,wc,mode)
  else if (mode.eq.1) then
    call setw(rs,cur_omega,mode) 
  else
    write(ilog,*) 'Fatal. Called sample_w(...) with unknown mode. Of&
 &fending mode # is ',mode,'.'
    call fexit()
  end if
!   
end
!
!-------------------------------------------------------------------------
!
! sample_nuc varies a selected number of nucleic acid backbone angles
! within a single residue
!
subroutine sample_nuc(rs,mode)
!
  use iounit
  use sequen
  use fyoc
  use movesets
!
  implicit none
!
  integer rs,mode,i
  RTYPE random,nnc(6)
!
! mode=0: sample fresh, mode=1: apply current sampling move
  if (mode.eq.0) then
    if (random().ge.nuc_randfreq) then
      do i=1,nnucs(rs)
        if (cur_nucflag(i).EQV..true.) then
          nnc(i) = nucs(i,rs) + (random()-0.5)*nuc_stepsz
          if (nnc(i).gt.180.0) nnc(i) = nnc(i) - 360.0
          if (nnc(i).lt.-180.0) nnc(i) = nnc(i) + 360.0
          cur_nucs(i) = nnc(i)
        else
          nnc(i) = nucs(i,rs)
          cur_nucs(i) = nucs(i,rs)
        end if
      end do
      call setnucs(rs,nnc,mode)
    else
      do i=1,nnucs(rs)
        if (cur_nucflag(i).EQV..true.) then
          nnc(i) = random()*360.0 - 180.0
          cur_nucs(i) = nnc(i)
        else
          nnc(i) = nucs(i,rs)
          cur_nucs(i) = nucs(i,rs)
        end if
      end do
      call setnucs(rs,nnc,mode)
    end if
  else if (mode.eq.1) then
    call setnucs(rs,cur_nucs,mode)
  else
    write(ilog,*) 'Fatal. Called sample_nuc(...) with unknown mode. &
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------
!
! sample_other varies a single dihedral angle directly addressing the Z matrix
!
subroutine sample_other(ati,mode,oth,isnat)
!
  use iounit
  use zmatrix
  use movesets
  use math
!
  implicit none
!
  integer, INTENT(IN):: ati,mode
  logical, INTENT(IN):: isnat
!
  RTYPE random,oth,doth
!
! mode=0: sample fresh, mode=1: apply current sampling move
  if (mode.eq.0) then
    if (random().ge.other_randfreq) then
      doth = (random()-0.5)*other_stepsz
      oth = ztor(ati) + doth
      if (oth.gt.180.0) oth = oth - 360.0
      if (oth.lt.-180.0) oth = oth + 360.0
      cur_other = oth
      call setother(ati,oth,mode,isnat)
    else
      oth = random()*360.0 - 180.0
      doth = oth - ztor(ati)
      if (doth.gt.180.0) doth = doth - 360.0
      if (doth.lt.-180.0) doth = doth + 360.0
      cur_other = oth
      call setother(ati,oth,mode,isnat)
    end if
  else if (mode.eq.1) then
    doth = cur_other - ztor(ati)
    if (doth.gt.180.0) doth = doth - 360.0
    if (doth.lt.-180.0) doth = doth + 360.0
    call setother(ati,cur_other,mode,isnat)
  else
    write(ilog,*) 'Fatal. Called sample_other(...) with unknown mode. &
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
  oth = doth/RADIAN
!
end
!
!----------------------------------------------------------------------
!
! sample_lct just varies a single LCT
!
subroutine sample_lct(ff)
!
  use iounit
  use sequen
  use fyoc
  use movesets
  use energies
!
  implicit none
!
  integer ff
  RTYPE random,fc
  integer aone
!
  aone = 1
!
! make this a flexible interval asap
  fc = random()*(lct_bnd(ff,2)-lct_bnd(ff,1)) + lct_bnd(ff,1)
  cur_lct = fc
  lctold(ff) = lct(ff)
  lct(ff) = fc
! this command will compute the corresponding native torsions,
! using the native set-commands (setfy and setchi2), which
! will do additional backups
  call lct2fyc()
! finally, here we'll build the correct Cartesian geometry from
! the now updated native internals
! WARNING: totally broken!!!!! 
  call makexyz_forbb(aone)
!   
end
!
!----------------------------------------------------------------------
!
! sample_bb_cr varies all the backbone angles included in a concerted 
! rotation move (rsi,...,rsf-1), this is more of a perturbation move
!
subroutine sample_bb_cr(rsi,rsf,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
!
  implicit none
!
  integer rsi,rsf,i
  RTYPE yc,fc,dphi(MAXCRDOF)
!
  do i=rsi,rsf-1
    fc = phi(i) + radian*dphi(2*(i-rsi)+1)
!    write(ilog,*) 'res.phi.:',i,radian*dphi(2*(i-rsi)+1)
    yc = psi(i) + radian*dphi(2*(i-rsi)+2)
!    write(ilog,*) 'res.psi.:',i,radian*dphi(2*(i-rsi)+2)
    call setfy(i,fc,yc,0)
  end do
!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!   
end
!
!-----------------------------------------------------------------------
!
subroutine sample_bb_torcrpr(rsi,rsf,dof,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
  use ujglobals
  use atoms
  use polypep
  use sequen
!
  implicit none
!
  integer rsi,rsf,dof,dofcnt,azero,aone,rs
  RTYPE yc,fc,wc,dphi(MAXUJDOF+6)
!
  dofcnt = 0
  azero = 0
  aone = 1
  do rs=rsf-2,rsi,-1
    fc = phi(rs)
    yc = psi(rs)
    wc = omega(rs)
!   backup
    call setw(rs,wc,azero)
    yc = psi(rs) + radian*dphi(dof-dofcnt)
    if (yc.lt.-180) yc = yc + 360.0
    if (yc.gt.180) yc = yc - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setfy(rs,fc,yc,azero)
      exit
    end if
    if (seqflag(rs).ne.5) then
      fc = phi(rs) + radian*dphi(dof-dofcnt)
      if (fc.lt.-180) fc = fc + 360.0
      if (fc.gt.180) fc = fc - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setfy(rs,fc,yc,azero)
        exit
      end if
    end if
    call setfy(rs,fc,yc,azero)
    wc = omega(rs) + radian*dphi(dof-dofcnt)
    if (wc.lt.-180) wc = wc + 360.0
    if (wc.gt.180) wc = wc - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setw(rs,wc,aone)
      exit
    end if
    call setw(rs,wc,aone)
  end do
!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!   
end
!
!-----------------------------------------------------------------------
!
subroutine sample_bb_torcrpr2(rsi,rsf,dof,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
  use ujglobals
  use atoms
  use polypep
  use sequen
!
  implicit none
!
  integer rsi,rsf,dof,dofcnt,azero,rs,aone
  RTYPE yc,fc,wc,dphi(MAXUJDOF+6)
!
  dofcnt = 0
  azero = 0
  aone = 1
  do rs=rsf-2,rsi,-1
    fc = phi(rs)
    yc = psi(rs)
    wc = omega(rs)
    call setw(rs,wc,azero)
    if (rs.lt.rsf-2) then
      yc = psi(rs) + radian*dphi(dof-dofcnt)
      if (yc.lt.-180) yc = yc + 360.0
      if (yc.gt.180) yc = yc - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setfy(rs,fc,yc,azero)
        exit
      end if
      if (seqflag(rs).ne.5) then
        fc = phi(rs) + radian*dphi(dof-dofcnt)
        if (fc.lt.-180) fc = fc + 360.0
        if (fc.gt.180) fc = fc - 360.0
        dofcnt = dofcnt + 1
        if (dofcnt.eq.dof) then
          call setfy(rs,fc,yc,azero)
          exit
        end if
      end if
      call setfy(rs,fc,yc,azero)
    end if
    wc = omega(rs) + radian*dphi(dof-dofcnt)
    if (wc.lt.-180) wc = wc + 360.0
    if (wc.gt.180) wc = wc - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setw(rs,wc,aone)
      exit
    end if
    call setw(rs,wc,aone)
  end do
!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!   
end
!
!-----------------------------------------------------------------------
!
subroutine sample_bb_torcrfull(rsi,rsf,dof,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
  use ujglobals
  use atoms
  use polypep
  use sequen
!
  implicit none
!
  integer rsi,rsf,dof,dofcnt,azero,rs
  RTYPE yc,fc,wc,dphi(MAXUJDOF+6)
!
  dofcnt = 0
  azero = 0
! closure segment
  do rs=rsf,rsf-1,-1
    yc = psi(rs) + radian*dphi(dof+(rs-rsf+1)*3+3)
    if (yc.lt.-180) yc = yc + 360.0
    if (yc.gt.180) yc = yc - 360.0
    fc = phi(rs) + radian*dphi(dof+(rs-rsf+1)*3+2)
    if (fc.lt.-180) fc = fc + 360.0
    if (fc.gt.180) fc = fc - 360.0
    call setfy(rs,fc,yc,azero)
    wc = omega(rs) + radian*dphi(dof+(rs-rsf+1)*3+1)
    if (wc.lt.-180) wc = wc + 360.0
    if (wc.gt.180) wc = wc - 360.0
    call setw(rs,wc,azero)    
  end do
!
! pre-rotation segment
  do rs=rsf-2,rsi,-1
    fc = phi(rs)
    yc = psi(rs)
    wc = omega(rs)
    yc = psi(rs) + radian*dphi(dof-dofcnt)
    if (yc.lt.-180) yc = yc + 360.0
    if (yc.gt.180) yc = yc - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setfy(rs,fc,yc,azero)
      exit
    end if
    if (seqflag(rs).ne.5) then
      fc = phi(rs) + radian*dphi(dof-dofcnt)
      if (fc.lt.-180) fc = fc + 360.0
      if (fc.gt.180) fc = fc - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setfy(rs,fc,yc,azero)
        exit
      end if
    end if
    call setfy(rs,fc,yc,azero)
    wc = omega(rs) + radian*dphi(dof-dofcnt)
    if (wc.lt.-180) wc = wc + 360.0
    if (wc.gt.180) wc = wc - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setw(rs,wc,azero)
      exit
    end if
    call setw(rs,wc,azero)
  end do
!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!   
end
!
!-----------------------------------------------------------------------
!
subroutine sample_bb_torcrfull2(rsi,rsf,dof,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
  use ujglobals
  use atoms
  use polypep
  use sequen
!
  implicit none
!
  integer rsi,rsf,dof,dofcnt,azero,rs
  RTYPE yc,fc,wc,dphi(MAXUJDOF+6)
!
  dofcnt = 0
  azero = 0
  do rs=rsf,rsf-2,-1
    yc = psi(rs) + radian*dphi(dof+(rs-rsf+2)*2+2)
    if (yc.lt.-180) yc = yc + 360.0
    if (yc.gt.180) yc = yc - 360.0
    fc = phi(rs) + radian*dphi(dof+(rs-rsf+2)*2+1)
    if (fc.lt.-180) fc = fc + 360.0
    if (fc.gt.180) fc = fc - 360.0
    call setfy(rs,fc,yc,azero)
    if (rs.gt.rsf-2) then
      wc = omega(rs)
      call setw(rs,wc,azero)    
    end if
  end do
!
  do rs=rsf-2,rsi,-1
    fc = phi(rs)
    yc = psi(rs)
    wc = omega(rs)
    if (rs.lt.rsf-2) then
      yc = psi(rs) + radian*dphi(dof-dofcnt)
      if (yc.lt.-180) yc = yc + 360.0
      if (yc.gt.180) yc = yc - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setfy(rs,fc,yc,azero)
        exit
      end if
      if (seqflag(rs).ne.5) then
        fc = phi(rs) + radian*dphi(dof-dofcnt)
        if (fc.lt.-180) fc = fc + 360.0
        if (fc.gt.180) fc = fc - 360.0
        dofcnt = dofcnt + 1
        if (dofcnt.eq.dof) then
          call setfy(rs,fc,yc,azero)
          exit
        end if
      end if
      call setfy(rs,fc,yc,azero)
    end if
    wc = omega(rs) + radian*dphi(dof-dofcnt)
    if (wc.lt.-180) wc = wc + 360.0
    if (wc.gt.180) wc = wc - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setw(rs,wc,azero)
      exit
    end if
    call setw(rs,wc,azero)
  end do

!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!   
end
!
!-----------------------------------------------------------------------
!
subroutine sample_bb_nuccrpr(rsi,rsf,dof,dphi,mode)
!
  use iounit
  use math
  use fyoc
  use movesets
  use ujglobals
  use sequen
  use molecule
!
  implicit none
!
  integer rsi,rsf,mode,dof,dofcnt,rs,aone,azero
  RTYPE newnucs(6),dphi(MAXUJDOF+6)
!
  dofcnt = 0
  aone = 1
  azero = 0
  do rs=rsf-1,rsi,-1
    newnucs(1:6) = nucs(1:6,rs)
    call setnucs(rs,newnucs,azero)
!   skip 4,5 for residue immediately preceding closure (used) and for caps (don't exist)
    if ((rs.lt.rsf-1).AND.(nucsline(4,rs).gt.0).AND.(nucsline(5,rs).gt.0)) then
      newnucs(5) = nucs(5,rs) + radian*dphi(dof-dofcnt)
      if (newnucs(5).lt.-180) newnucs(5) = newnucs(5) + 360.0
      if (newnucs(5).gt.180) newnucs(5) = newnucs(5) - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setnucs(rs,newnucs,aone)
        exit
      end if
      newnucs(4) = nucs(4,rs) + radian*dphi(dof-dofcnt)
      if (newnucs(4).lt.-180) newnucs(4) = newnucs(4) + 360.0
      if (newnucs(4).gt.180) newnucs(4) = newnucs(4) - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setnucs(rs,newnucs,aone)
        exit
      end if
    end if
    newnucs(3) = nucs(3,rs) + radian*dphi(dof-dofcnt)
    if (newnucs(3).lt.-180) newnucs(3) = newnucs(3) + 360.0
    if (newnucs(3).gt.180) newnucs(3) = newnucs(3) - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setnucs(rs,newnucs,aone)
      exit
    end if
    newnucs(2) = nucs(2,rs) + radian*dphi(dof-dofcnt)
    if (newnucs(2).lt.-180) newnucs(2) = newnucs(2) + 360.0
    if (newnucs(2).gt.180) newnucs(2) = newnucs(2) - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setnucs(rs,newnucs,aone)
      exit
    end if
!   note that 1 is never a suitable d.o.f if rsi is terminal, but should be skipped automatically
    newnucs(1) = nucs(1,rs) + radian*dphi(dof-dofcnt)
    if (newnucs(1).lt.-180) newnucs(1) = newnucs(1) + 360.0
    if (newnucs(1).gt.180) newnucs(1) = newnucs(1) - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setnucs(rs,newnucs,aone)
      exit
    end if
  end do
!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!   
end
!
!-----------------------------------------------------------------------
!
subroutine sample_bb_nuccrfull(rsi,rsf,dof,dphi,mode)
!
  use iounit
  use math
  use fyoc
  use movesets
  use ujglobals
  use sequen
!
  implicit none
!
  integer rsi,rsf,mode,dof,dofcnt,rs,azero,aone
  RTYPE newnucs(6),dphi(MAXUJDOF+6)
!
  dofcnt = 0
  aone = 1
  azero = 0
!
  newnucs(1:6) = nucs(1:6,rsf)
  newnucs(1) = nucs(1,rsf) + radian*dphi(dof+3)
  if (newnucs(1).lt.-180) newnucs(1) = newnucs(1) + 360.0
  if (newnucs(1).gt.180) newnucs(1) = newnucs(1) - 360.0
  newnucs(2) = nucs(2,rsf) + radian*dphi(dof+4)
  if (newnucs(2).lt.-180) newnucs(2) = newnucs(2) + 360.0
  if (newnucs(2).gt.180) newnucs(2) = newnucs(2) - 360.0
  newnucs(3) = nucs(3,rsf) + radian*dphi(dof+5)
  if (newnucs(3).lt.-180) newnucs(3) = newnucs(3) + 360.0
  if (newnucs(3).gt.180) newnucs(3) = newnucs(3) - 360.0
  newnucs(4) = nucs(4,rsf) + radian*dphi(dof+6)
  if (newnucs(4).lt.-180) newnucs(4) = newnucs(4) + 360.0
  if (newnucs(4).gt.180) newnucs(4) = newnucs(4) - 360.0
  call setnucs(rsf,newnucs,azero)
  newnucs(1:6) = nucs(1:6,rsf-1)
  newnucs(4) = nucs(4,rsf-1) + radian*dphi(dof+1)
  if (newnucs(4).lt.-180) newnucs(4) = newnucs(4) + 360.0
  if (newnucs(4).gt.180) newnucs(4) = newnucs(4) - 360.0
  newnucs(5) = nucs(5,rsf-1) + radian*dphi(dof+2)
  if (newnucs(5).lt.-180) newnucs(5) = newnucs(5) + 360.0
  if (newnucs(5).gt.180) newnucs(5) = newnucs(5) - 360.0
  call setnucs(rsf-1,newnucs,aone)
!
  do rs=rsf-1,rsi,-1
    newnucs(1:6) = nucs(1:6,rs)
!   skip 4,5 for residue immediately preceding closure (used) and for caps (don't exist)
    if ((rs.lt.rsf-1).AND.(nucsline(4,rs).gt.0).AND.(nucsline(5,rs).gt.0)) then
      newnucs(5) = nucs(5,rs) + radian*dphi(dof-dofcnt)
      if (newnucs(5).lt.-180) newnucs(5) = newnucs(5) + 360.0
      if (newnucs(5).gt.180) newnucs(5) = newnucs(5) - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setnucs(rs,newnucs,aone)
        exit
      end if
      newnucs(4) = nucs(4,rs) + radian*dphi(dof-dofcnt)
      if (newnucs(4).lt.-180) newnucs(4) = newnucs(4) + 360.0
      if (newnucs(4).gt.180) newnucs(4) = newnucs(4) - 360.0
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) then
        call setnucs(rs,newnucs,aone)
        exit
      end if
    end if
    newnucs(3) = nucs(3,rs) + radian*dphi(dof-dofcnt)
    if (newnucs(3).lt.-180) newnucs(3) = newnucs(3) + 360.0
    if (newnucs(3).gt.180) newnucs(3) = newnucs(3) - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setnucs(rs,newnucs,aone)
      exit
    end if
    newnucs(2) = nucs(2,rs) + radian*dphi(dof-dofcnt)
    if (newnucs(2).lt.-180) newnucs(2) = newnucs(2) + 360.0
    if (newnucs(2).gt.180) newnucs(2) = newnucs(2) - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setnucs(rs,newnucs,aone)
      exit
    end if
    newnucs(1) = nucs(1,rs) + radian*dphi(dof-dofcnt)
    if (newnucs(1).lt.-180) newnucs(1) = newnucs(1) + 360.0
    if (newnucs(1).gt.180) newnucs(1) = newnucs(1) - 360.0
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) then
      call setnucs(rs,newnucs,aone)
      exit
    end if
  end do
!
! remember that this will re-do the complete chain beyond and including rsi 
  call makexyz_forbb(rsi)
!
end
!
!-------------------------------------------------------------------------
!
subroutine sample_bb_ujpr(rsi,rsf,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
  use atoms
  use zmatrix
  use polypep
  use ujglobals
!
  implicit none
!
  integer rsi,rsf,i
  RTYPE yc,fc,dphi(MAXUJDOF+6),ac,bc,gc
!
! increment angles and dihedrals by passed on vector
  do i=rsi,rsf
    fc = phi(i) + dphi(2*(i-rsi)+1)
    yc = psi(i) + dphi(2*(i-rsi)+2)
! 
    ac = bang(ci(i)) + dphi(3*(i-rsi) + cur_ujsz*2 + 1)
    bc = bang(ni(i+1)) + dphi(3*(i-rsi) + cur_ujsz*2 + 2)
    gc = bang(cai(i+1)) + dphi(3*(i-rsi) + cur_ujsz*2 + 3)
!
    call setfy(i,fc,yc,0)
    call setabg(i,ac,bc,gc)
  end do
! 
  call makexyz_forbb(rsi)
!
end
!
!-----------------------------------------------------------------------
!
! sets internals for chain closure segment
! relies on extraneous call to makexyz!!
!
subroutine sample_bb_ujcc(rsf,dphi)
!
  use iounit
  use math
  use fyoc
  use movesets
  use atoms
  use zmatrix
  use polypep
  use ujglobals
!
  implicit none
!
  integer rsf
  RTYPE dphi(MAXUJDOF+6),ac,bc,gc
  RTYPE fc1,fc2,yc1,yc2
!
! increment angles and dihedrals by passed on vector: this is tricky
! because of the messed up ordering within dphi
  fc1 = phi(rsf-1) + dphi(5*cur_ujsz+1)
  yc1 = psi(rsf-1) + dphi(5*cur_ujsz+2)
  fc2 = phi(rsf) + dphi(5*cur_ujsz+3)
  yc2 = psi(rsf)
! 
  ac = bang(ci(rsf-1)) + dphi(5*cur_ujsz+4)
  bc = bang(ni(rsf)) + dphi(5*cur_ujsz+5)
  gc = bang(cai(rsf)) + dphi(5*cur_ujsz+6)
!
  call setfy(rsf-1,fc1,yc1,0)
  call setfy(rsf,fc2,yc2,0)
  call setabg(rsf-1,ac,bc,gc)
! 
! note that makexyz must be executed later in a different fxn!!
!
end
!
!----------------------------------------------------------------------
!
! THIRD: set routines for the different degrees freedom
!
!----------------------------------------------------------------------
!
! "setfy" sets the phi-psi values from within the allowed regions
! of a specified mesostate box
!     
! res  -  Residue ID (integer)
! fc   -  phi angle at the center of the chosen grid in phi,psi space
! yc   -  psi angle at the center of the chosen grid in phi,psi space
! mode -  if mode != 1 then set to exactly fc, yc, 
!         if mode == 1, then use grids (deprecated)
!
subroutine setfy(res,fc,yc,mode)
!
  use fyoc
  use grids
  use zmatrix
  use sequen
!
  implicit none
!
  integer ff,res,yy,y2,f2,mode,styp
  RTYPE displ,fc,random,yc
!
  styp = seqtyp(res)
!
! z-matrix line numbers for backbone torsions
  ff = fline(res)
  f2 = fline2(res)
  yy = yline(res)
  y2 = yline2(res)
!
! generate random displacements and assign desired phi,psi values
! first phi
  if (ff.gt.0) then
    ztorpr(ff) = ztor(ff)
    if (mode.eq.1) then
      displ = 2.0*(random()-0.5)
      ztor(ff) = stgr%sfyc(styp)*displ + fc
      cur_phi = ztor(ff)
    else
      ztor(ff) = fc
    end if
    call wrapper(ff)
!   update pointer array
    phi(res) = ztor(ff)
! set other torsions that depend on phi
    if (f2.gt.0) then
      ztorpr(f2) = ztor(f2)
      ztor(f2) = ztor(ff) + phish(res)
      if (ztor(f2).gt.180.0) ztor(f2) = ztor(f2) - 360.0
      if (ztor(f2).lt.-180.0) ztor(f2) = ztor(f2) + 360.0
    end if
  end if
!
! then psi
!
  if (yy.gt.0) then
    ztorpr(yy) = ztor(yy)
    if (mode.eq.1) then
      displ = 2.0*(random()-0.5)
      ztor(yy) = stgr%sfyc(styp)*displ + yc
      cur_psi = ztor(yy)
    else
      ztor(yy) = yc
    end if
    call wrapper(yy)
!   update pointer array
    psi(res) = ztor(yy)
!   set other torsions that depend on psi
    if (y2.gt.0) then
      ztorpr(y2) = ztor(y2)
      ztor(y2) = ztor(yy) + psish(res)
      if (ztor(y2).gt.180.0) ztor(y2) = ztor(y2) - 360.0
      if (ztor(y2).lt.-180.0) ztor(y2) = ztor(y2) + 360.0
    end if
!   reset the psi-angle according to zmatrix format
    if(ztor(yy).lt.0.0d0) then
      ztor(yy) = ztor(yy) + 180.0d0
    else
      ztor(yy) = ztor(yy) - 180.0d0
    end if
  end if
!
end
!
!----------------------------------------------------------------------------
!
! rs  -  Residue ID
! wc - omega angle 
! mode 0: sample freshly
! mode 1: re-sample (no back-up)
!
subroutine setw(rs,wc,mode)
!
  use fyoc
  use iounit
  use zmatrix
!
  implicit none
!
  integer ll,rs,mode
  RTYPE wc
!
! fresh sampling
  if (mode.eq.0) then
!   back-up current torsion and set new one
    ll = wline(rs)
    if (ll.gt.0) then
      ztorpr(ll) = ztor(ll)
      omega(rs) = wc
      ztor(ll) = wc
    end if
! re-sampling
  else if (mode.eq.1) then
!   just re-assign
    ll = wline(rs)
    if (ll.gt.0) then
      omega(rs) = wc
      ztor(ll) = wc
    end if
  else
    write(ilog,*) 'Fatal. Called setw(...) with unknown mode.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------
!
subroutine setchi(rs,cc,mode)
!
  use iounit
  use fyoc
  use zmatrix
!
  implicit none
!
  integer i,ll,rs,mode
  RTYPE cc(MAXCHI)
!
! fresh sampling
  if (mode.eq.0) then
!   back-up current torsions and set new ones
    do i=1,nchi(rs)
      ll = chiline(i,rs)
      ztorpr(ll) = ztor(ll)
      chi(i,rs) = cc(i)
      ztor(ll) = cc(i)
    end do
! re-sampling
  else if (mode.eq.1) then
!   just re-assign
    do i=1,nchi(rs)
      ll = chiline(i,rs)
      chi(i,rs) = cc(i)
      ztor(ll) = cc(i)
    end do
  else
    write(ilog,*) 'Fatal. Called setchi(...) with unknown mode.'
  end if
!
end
!
!----------------------------------------------------------------------------
!
subroutine setother(ati,torval,mode,isnat)
!
  use iounit
  use zmatrix
  use sequen
  use molecule
  use fyoc
  use atoms
!
  implicit none
!
  integer ati,mode,rsl(3),nrs,j,rs
  RTYPE torval
  logical isnat
!
! fresh sampling
  if (mode.eq.0) then
!   back-up current torsion and set new one
    ztorpr(ati) = ztor(ati)
    ztor(ati) = torval
! re-sampling
  else if (mode.eq.1) then
!   just re-assign
    ztor(ati) = torval
  else
    write(ilog,*) 'Fatal. Called setother(...) with unknown mode.'
    call fexit()
  end if
!
  if (isnat.EQV..true.) then
    rsl(1) = atmres(ati)
    nrs = 1
    if (atmres(ati).gt.rsmol(molofrs(atmres(ati)),1)) then 
      nrs = nrs + 1
      rsl(nrs) = atmres(ati)-1
    end if
    if (atmres(ati).lt.rsmol(molofrs(atmres(ati)),2)) then
      nrs = nrs + 1
      rsl(nrs) = atmres(ati)+1
    end if
    do rs=1,nrs
      if (ati.eq.wline(rsl(rs))) omega(rsl(rs)) = ztor(ati)
      if (ati.eq.fline(rsl(rs))) phi(rsl(rs)) = ztor(ati)
      if (ati.eq.yline(rsl(rs)))  then
        if(ztor(ati).lt.0.0) then
          psi(rsl(rs)) = ztor(ati) + 180.0
        else
          psi(rsl(rs)) = ztor(ati) - 180.0
        end if
      end if
      do j=1,nnucs(rsl(rs))
        if (ati.eq.nucsline(j,rsl(rs))) nucs(j,rsl(rs)) = ztor(ati)
      end do
      do j=1,nchi(rsl(rs))
        if (ati.eq.chiline(j,rsl(rs))) chi(j,rsl(rs)) = ztor(ati)
      end do
    end do
  end if
!
end
!
!----------------------------------------------------------------------------
!
subroutine setnucs(rs,nnc,mode)
!
  use iounit
  use fyoc
  use zmatrix
!
  implicit none
!
  integer i,ll,rs,mode
  RTYPE nnc(6)
!
! fresh sampling
  if (mode.eq.0) then
!   back-up current torsions and set new ones
    do i=1,nnucs(rs)
      ll = nucsline(i,rs)
      ztorpr(ll) = ztor(ll)
      nucs(i,rs) = nnc(i)
      ztor(ll) = nnc(i)
    end do
! re-sampling
  else if (mode.eq.1) then
!   just re-assign
    do i=1,nnucs(rs)
      ll = nucsline(i,rs)
      nucs(i,rs) = nnc(i)
      ztor(ll) = nnc(i)
    end do
  else
    write(ilog,*) 'Fatal. Called setnucs(...) with unknown mode.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a routine to set the three bond angles for a residue starting at its
! N, CA, and C-atoms respectively down the backbone chain
! spills into next residue -> sanity checks ...
!
subroutine setabg(res,ac,bc,gc)
!
  use zmatrix
  use atoms
  use polypep
  use iounit
  use sequen
  use molecule
!
  implicit none
!
  integer res
  integer aa,bb,gg
  RTYPE ac,bc,gc
!
  if ((res.eq.rsmol(molofrs(res),1)).AND.(seqflag(res).eq.10)) then
    write(ilog,*) 'Fatal. Called setabg(...) with N-terminal cap res&
 &idue. This is a bug. Please report.'
    call fexit()
  else if ((res.eq.rsmol(molofrs(res),2)).AND.(seqtyp(res).eq.30)) then
    write(ilog,*) 'Fatal. Called setabg(...) with residue immediatel&
 &y ahead of C-terminal NH2-cap. Please report this bug.'
    call fexit()
  end if
!
! ! z-matrix line numbers for backbone torsions
  aa = ci(res)
  bb = ni(res+1)
  gg = cai(res+1)
!
  bangpr(aa) = bang(aa)
  bangpr(bb) = bang(bb)
  bangpr(gg) = bang(gg)
!
  bang(aa) = ac
  bang(bb) = bc
  bang(gg) = gc
!
end
!!
!!-----------------------------------------------------------------------
!!
!!! the subroutine to set the pucker angle for proline
!!
!!! the rr is the identifier which conformation to adopt (>0.5 vs. <0.5)
!!
!subroutine! pucker(rs,rr,mode)
!!
!  use fyoc
!  use zmatrix
!  use iounit
!!
!  implicit none
!!
!  integer mode,rs
!  integer l1,l2
!  RTYPE rr
!!
!!! set ring geometry for proline to be exo or endo, 
!!! following the structure database distribution
!  if ((mode.lt.0).OR.(mode.gt.1)) then
!    write(ilog,*) 'Fatal. Called pucker with unknown mode (',&
! &mode,').'
!    call fexit()
!  end if
!!!     
!!! store old values, set chi angles based on chosen pucker mode
!  l1 = chiline(1,rs)
!  l2 = chiline(2,rs)
!  if (mode.eq.0) then
!    ztorpr(l1) = ztor(l1)
!    ztorpr(l2) = ztor(l2)
!  end if
!!
!!! WARNING PHI MISSING!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if (rr .lt. 0.5d0) then
!    if (chiral(rs).eq.1) then
!      ztor(l1) = -28.0d0
!      ztor(l2) = 4.3d0  
!    else
!      ztor(l1) = 31.0d0 !! 28.0d0
!      ztor(l2) = -11.9d0
!    end if
!!!    ztor(l2) = -1.3503d0*ztor(l1) + 3.2392d0
!    chi(1,rs) = ztor(l1)
!    chi(2,rs) = ztor(l2)
!  else
!    if (chiral(rs).eq.1) then
!      ztor(l1) =  31.0d0 !! 32.0d0
!      ztor(l2) = -11.9d0 
!    else
!      ztor(l1) =  -29.0d0 !! -32.0d0
!      ztor(l2) = 6.4d0 
!    end if
!!!    ztor(l2) = -1.3503d0*ztor(l1) + 3.2392d0
!    chi(1,rs) = ztor(l1)
!    chi(2,rs) = ztor(l2)
!  end if
!!
!  if (mode.eq.0) then
!    cur_chis(1) = chi(1,rs)
!    cur_chis(2) = chi(2,rs)
!  end if
!!
!end
!!
!!
!----------------------------------------------------------------------
!
! FOURTH: unset (restore) routines for the different degrees freedom
!
!----------------------------------------------------------------------
!    
! recover internals
!
subroutine restore(res)
!
  use fyoc
  use sequen
  use zmatrix
  use movesets
!
  implicit none
!
  integer i,res
  integer cc,ff,yy,y2,f2,oo
!
! restore phi and psi
  ff = fline(res)
  yy = yline(res)
  y2 = yline2(res)
  f2 = fline2(res)
  if (ff.gt.0) then
    ztor(ff) = ztorpr(ff)
    phi(res) = ztor(ff)
  end if
  if (yy.gt.0) then
    ztor(yy) = ztorpr(yy)
    if(ztor(yy).lt.0.0) then
     psi(res) = ztor(yy) + 180.0
    else
     psi(res) = ztor(yy) - 180.0
    end if
  end if
  if (y2.gt.0) ztor(y2) = ztorpr(y2)
  if (f2.gt.0) ztor(f2) = ztorpr(f2)     
!
! restore chi angles
  do i=1,nchi(res)
    cc = chiline(i,res)
    ztor(cc) = ztorpr(cc)
    chi(i,res) = ztor(cc)
  end do
!
! restore omega angles
  oo = wline(res)
  if (oo.gt.0) then
    ztor(oo) = ztorpr(oo)
    omega(res) = ztor(oo)
  end if
!
end
!
!------------------------------------------------------------------
!
subroutine restore_forfy(res)
!
  use fyoc
  use zmatrix
!
  implicit none
!
  integer res
  integer ff,yy,y2,f2
!
! restore phi and psi
  ff = fline(res)
  yy = yline(res)
  y2 = yline2(res)
  f2 = fline2(res)
  if (ff.gt.0) then
    ztor(ff) = ztorpr(ff)
    phi(res) = ztor(ff)
  end if
  if (yy.gt.0) then
    ztor(yy) = ztorpr(yy)
    if(ztor(yy).lt.0.0) then
      psi(res) = ztor(yy) + 180.0
    else
      psi(res) = ztor(yy) - 180.0
    end if
  end if
  if (y2.gt.0) ztor(y2) = ztorpr(y2)
  if (f2.gt.0) ztor(f2) = ztorpr(f2)    
!
end
!
!------------------------------------------------------------------
!
subroutine restore_forw(rs)
!
  use fyoc
  use zmatrix
!
  implicit none
!
  integer oo,rs
!
  oo = wline(rs)
  if (oo.gt.0) then
    ztor(oo) = ztorpr(oo)
    omega(rs) = ztor(oo)
  end if
!
end
!
!
!------------------------------------------------------------------
!
subroutine restore_fornuc(rs)
!
  use fyoc
  use zmatrix
!
  implicit none
!
  integer rs,i,cc
!
  do i=1,nnucs(rs)
    cc = nucsline(i,rs)
    ztor(cc) = ztorpr(cc)
    nucs(i,rs) = ztor(cc)
  end do
!
end
!
!------------------------------------------------------------------
!
subroutine restore_forsc(rs)
!
  use fyoc
  use zmatrix
!
  implicit none
!
  integer rs,i,cc
!
  do i=1,nchi(rs)
    cc = chiline(i,rs)
    ztor(cc) = ztorpr(cc)
    chi(i,rs) = ztor(cc)
  end do
!
end
!
!------------------------------------------------------------------
!
subroutine restore_forother(ati,isnat)
!
  use zmatrix
  use molecule
  use sequen
  use atoms
  use fyoc
!
  implicit none
!
  integer ati,rs,rsl(3),nrs,j
  logical isnat
!
  ztor(ati) = ztorpr(ati)
!
  if (isnat.EQV..true.) then
    rsl(1) = atmres(ati)
    nrs = 1
    if (atmres(ati).gt.rsmol(molofrs(atmres(ati)),1)) then 
      nrs = nrs + 1
      rsl(nrs) = atmres(ati)-1
    end if
    if (atmres(ati).lt.rsmol(molofrs(atmres(ati)),2)) then
      nrs = nrs + 1
      rsl(nrs) = atmres(ati)+1
    end if
    do rs=1,nrs
      if (ati.eq.wline(rsl(rs))) omega(rsl(rs)) = ztor(ati)
      if (ati.eq.fline(rsl(rs))) phi(rsl(rs)) = ztor(ati)
      if (ati.eq.yline(rsl(rs)))  then
        if(ztor(ati).lt.0.0) then
          psi(rsl(rs)) = ztor(ati) + 180.0
        else
          psi(rsl(rs)) = ztor(ati) - 180.0
        end if
      end if
      do j=1,nnucs(rsl(rs))
        if (ati.eq.nucsline(j,rsl(rs))) nucs(j,rsl(rs)) = ztor(ati)
      end do
      do j=1,nchi(rsl(rs))
        if (ati.eq.chiline(j,rsl(rs))) chi(j,rsl(rs)) = ztor(ati)
      end do
    end do
  end if
!
end
!
!---------------------------------------------------------------
!
subroutine restore_poly(imol)
!
  use molecule
!
  implicit none
!
  integer i,imol,j
!
  rgv(imol) = rgvref(imol)
  do i=1,3
    rgevs(imol,i) = rgevsref(imol,i)
    com(imol,i) = comref(imol,i)
  end do
  do i=1,3
    do j=1,3
      rgpcs(imol,i,j) = rgpcsref(imol,i,j)
    end do
  end do
!
end
!
!----------------------------------------------------------------
!
subroutine restore_forring(rs)
!
  use fyoc
  use zmatrix
  use polypep
  use sequen
  use molecule
  use iounit
  use system
  use aminos
!
  implicit none
!
  integer rs,shf
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! cyclic peptide residues
  if (pucline(rs).gt.0) then
    ztor(at(rs)%sc(2-shf)) = ztorpr(at(rs)%sc(2-shf)) 
    ztor(at(rs)%sc(3-shf)) = ztorpr(at(rs)%sc(3-shf))
    ztor(at(rs)%sc(4-shf)) = ztorpr(at(rs)%sc(4-shf))
    bang(at(rs)%sc(2-shf)) = bangpr(at(rs)%sc(2-shf))
    bang(at(rs)%sc(3-shf)) = bangpr(at(rs)%sc(3-shf))
    ztor(pucline(rs)) = ztorpr(pucline(rs))
    bang(at(rs)%sc(4-shf)) = bangpr(at(rs)%sc(4-shf))
! nucleotides
  else if (nucsline(6,rs).gt.0) then
    ztor(at(rs)%sc(1)) = ztorpr(at(rs)%sc(1)) 
    ztor(at(rs)%sc(2)) = ztorpr(at(rs)%sc(2))
    ztor(at(rs)%sc(3)) = ztorpr(at(rs)%sc(3))
    bang(at(rs)%sc(1)) = bangpr(at(rs)%sc(1))
    bang(at(rs)%sc(2)) = bangpr(at(rs)%sc(2))
    bang(at(rs)%sc(3)) = bangpr(at(rs)%sc(3))
    ztor(nucsline(6,rs)) = ztorpr(nucsline(6,rs))
  else
    write(ilog,*) 'Fatal. Called restore_forring with unsupported residue (',rs,' of&
 & type ',seqtyp(rs),': ',amino(seqtyp(rs)),').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------
!
subroutine restore_forabg(res)
!
  use fyoc
  use zmatrix
  use polypep
  use iounit
  use sequen
  use molecule
!
  implicit none
!
  integer res
  integer aa,bb,gg
!
  if ((res.eq.rsmol(molofrs(res),1)).AND.(seqflag(res).eq.10)) then
    write(ilog,*) 'Fatal. Called restore_forabg(...) with N-terminal&
 & cap residue. This is a bug. Please report.'
    call fexit()
  else if ((res.eq.rsmol(molofrs(res),2)).AND.(seqtyp(res).eq.30)) then
    write(ilog,*) 'Fatal. Called restore_forabg(...) with residue im&
 &mediately ahaed of C-terminal NH2-cap. Please report this bug.'
    call fexit()
  end if
!      
! restore the N-CA-C (aa), CA-C-N+1 (bb) and C-N+1-CA+1 (gg) angles
  aa = ci(res)
  bb = ni(res+1)
  gg = cai(res+1)
!
  bang(aa) = bangpr(aa)
  bang(bb) = bangpr(bb)
  bang(gg) = bangpr(gg)
!
end
!
!-----------------------------------------------------------------------------
!
! WARNING: not supported at the moment
!
subroutine restore_forlct(res)
!
  use fyoc
  use aminos
  use sequen
!
  implicit none
!
  integer res,i
!
! restore current LCT valuue
  lct(res) = lctold(res)
! now restore the corresponding internals
  do i=1,nseq
    if ((seqflag(i).eq.10).OR.(seqflag(i).eq.11).OR.(seqflag(i).eq.13)) then 
      cycle
    end if
    call restore(i)
  end do
! the xyz-restoration is handled outside (as usual)
!
end
!
!------------------------------------------------------------------------------
!
! "wrapper" corrects wrap around artifacts for a specific zmatrix
! torsion
!
subroutine wrapper(line)
!
  use zmatrix
!
  implicit none
!
  integer line
!
! correct wrap around artifiacts: interval is always [-180,180]
  if(ztor(line) .lt. -180.0d0 ) then 
    ztor(line) = ztor(line) + 360.0d0
  else if(ztor(line) .gt. 180.0d0) then
    ztor(line) = ztor(line) - 360.0d0
  end if
!
end
!
!--------------------------------------------------------------------------
!
