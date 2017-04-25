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
! MAIN AUTHOR:   Nicholas Lyle                                             !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!     
!
#include "macros.i"
!
!-----------------------------------------------------------------------------
! 
! this initialization routine performs a few checks, initializes values, and allocates memory
!
subroutine wl_init()
!
  use wl
  use molecule
  use iounit
!  
  implicit none
!  
  integer i,t1,t2
!
  if ((wld%wl_mode.gt.3).OR.(wld%wl_mode.lt.1)) then
    write(ilog,*) 'Fatal. Encountered unsupported mode for Wang-Landau sampling in wl_init(). This is a bug.'
    call fexit()
  end if
  if ((wld%wl_mode.eq.2).OR.((wld%wl_mode.eq.3).AND.(wld%dimensionality.eq.2))) then
    if ((wld%wl_mol(1).lt.1).OR.(wld%wl_mol(1).gt.nmol)) then
      write(ilog,*) 'Fatal. Selected molecule for Wang-Landau sampling (mode 2 or 3) does not exist. Please check input.'
      call fexit()
    else if ((atmol(wld%wl_mol(1),1).eq.atmol(wld%wl_mol(1),2)).OR.(ntormol(moltypid(wld%wl_mol(1))).le.0)) then
      write(ilog,*) 'Warning. Selected molecule for Wang-Landau method (mode 2 or 3) does or &
 &may lack internal flexibility (reaction coordinate is potentially invariant).'
    end if
  end if
  if ((wld%wl_mode.eq.2).AND.(wld%dimensionality.eq.2)) then
    if ((wld%wl_mol(2).lt.1).OR.(wld%wl_mol(2).gt.nmol)) then
      write(ilog,*) 'Fatal. Second selected molecule for 2D Wang-Landau sampling (mode 2) does not exist. Please check input.'
      call fexit()
    else if ((wld%wl_mol(2).eq.wld%wl_mol(1)).AND.(wld%wl_rc(2).eq.wld%wl_rc(1))) then
      write(ilog,*) 'Fatal. Selected identical molecules and reaction coordinates for 2D Wang-Landau sampling (mode 2). This &
 &is fatal (use 1D instead).'
      call fexit()
    else if ((atmol(wld%wl_mol(2),1).eq.atmol(wld%wl_mol(2),2)).OR.(ntormol(moltypid(wld%wl_mol(2))).le.0)) then
      write(ilog,*) 'Warning. Second selected molecule for 2D Wang-Landau method (mode 2) does or &
 &may lack internal flexibility (reaction coordinate is potentially invariant).'
    end if
  end if
  if ((wld%wl_mode.eq.3).AND.(wld%dimensionality.eq.2)) then
    wld%wl_mol(2) = wld%wl_mol(1)
    wld%wl_rc(2) = wld%wl_rc(1) ! 1 is Rg, 2 is ZA, 3 is ZB
    wld%wl_rc(1) = 0 ! 0 means use global(!) energy instead
  end if
!
! set up WL related histograms (note that g_max / g_max2d and g_binsz / g_binsz2d are provided by user input
  t1 = 100
  t2 = 100
  call wl_realloc_histos(t1,t2)
  if (wld%dimensionality.eq.1) then
    wld%g_min = wld%g_max - (wld%nbins-1)*wld%g_binsz
    wld%minv = wld%g_min - 0.5*wld%g_binsz
    wld%bctr(1) = wld%g_min
    do i=2,wld%nbins
     wld%bctr(i) = wld%bctr(i-1)+wld%g_binsz
    end do
    wld%minb = wld%nbins
  else if (wld%dimensionality.eq.2) then
    wld%g_min2d(:) = wld%g_max2d(:) - (wld%nbins2d(:)-1)*wld%g_binsz2d(:)
    wld%minv2d(:) = wld%g_min2d(:) - 0.5*wld%g_binsz2d(:)
    wld%bctr2d1(1) = wld%g_min2d(1)
    wld%bctr2d2(1) = wld%g_min2d(2)
    do i=2,wld%nbins2d(1)
       wld%bctr2d1(i) = wld%bctr2d1(i-1)+wld%g_binsz2d(1)
    end do
    do i=2,wld%nbins2d(2)
       wld%bctr2d2(i) = wld%bctr2d2(i-1)+wld%g_binsz2d(2)
    end do
    wld%minb2d(:) = wld%nbins2d(:)
  end if
!
! read in initialization files if requested (note that this may resize histograms)
  if (wld%use_ginitfile) then
    call strlims(wld%ginitfile,t1,t2)
    call wl_parse_ginitfile(wld%ginitfile(t1:t2))
  end if
!
end subroutine wl_init
!
!------------------------------------------------------------------------------------------
!
subroutine wl_update_histos(istep,idx_new,idx_new2)
!
  use wl
  use system
  use iounit
#ifdef ENABLE_MPI
  use mpistuff
  use interfaces
#endif
!
  implicit none
!
  integer idx_new,istep,idx_new2
  character(MAXSTRLEN) bname
  RTYPE critv
#ifdef ENABLE_MPI
  integer itmp
  integer, ALLOCATABLE:: gtmp2d(:,:,:),gtmp(:,:)
#endif
!
  if ((idx_new.le.0).OR.((wld%dimensionality.eq.1).AND.(idx_new.gt.wld%nbins)).OR.&
 & ((wld%dimensionality.eq.2).AND.((idx_new2.le.0).OR.(idx_new2.gt.wld%nbins2d(2)).OR.&
 &                                 (idx_new.gt.wld%nbins2d(1))))) then
    write(ilog,*) 'Fatal. Subroutine wl_update_histos(...) was called with corrupt bin indices. &
 &This is either a bug or it could indicate a numerically corrupt simulation.'
    call fexit()
  end if
!
! update the histograms every wld%hufreq MC steps
  if (mod(wld%stepnum,wld%hufreq).eq.0) then
!
!    Update histograms
    if (wld%dimensionality.eq.1) then
      wld%g(idx_new) =  wld%g(idx_new) + wld%fval  !fval is LOG(f) now
      wld%gh(idx_new) =  wld%gh(idx_new) + 1
      wld%ghbu(idx_new) =  wld%ghbu(idx_new) + 1
      wld%ghtot(idx_new) =  wld%ghtot(idx_new) + 1
    else if (wld%dimensionality.eq.2) then
      wld%g2d(idx_new,idx_new2) =  wld%g2d(idx_new,idx_new2) + wld%fval  !fval is LOG(f) now
      wld%gh2d(idx_new,idx_new2) =  wld%gh2d(idx_new,idx_new2) + 1
      wld%ghbu2d(idx_new,idx_new2) =  wld%ghbu2d(idx_new,idx_new2) + 1
      wld%ghtot2d(idx_new,idx_new2) =  wld%ghtot2d(idx_new,idx_new2) + 1
    end if
  end if
!
! Check energy histogram and update f value if needed
! valid only for first part of the algorithm
  if (wld%t1on.EQV..false.) then
    if (mod(wld%stepnum,wld%gh_flatcheck_freq).eq.0) then
#ifdef ENABLE_MPI
      if ((use_MPIAVG.EQV..true.).AND.((wld%freeze.EQV..false.).OR.(wld%flevel.eq.0))) then
        call wl_parallel_syncbins(idx_new,idx_new2)
      end if
#endif
      call wl_update_convergence(istep)
    end if
  end if
!
! Check if we switch to the 1/t portion of the algorithm.
! once active we never switch back
  if (wld%dimensionality.eq.1) then
    critv = dble(wld%nbins)
  else if (wld%dimensionality.eq.2) then
    critv = dble(wld%nbins2d(1)*wld%nbins2d(2))
  end if
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    itmp = mpi_nodes*wld%stepnum ! note that for hybrid MC/MD, the cycles are synchronized
    if (wld%t1on.EQV..true.) then
      if (mod(wld%stepnum,wld%gh_flatcheck_freq).eq.0) then
        if (wld%dimensionality.eq.1) then
          allocate(gtmp(wld%nbins,4))
          call wl_parallel_combine(gtmp=gtmp)
          if ((wld%ghtot(wld%minb).le.0).OR.(wld%ghtot(wld%maxb).le.0)) then
            write(ilog,*) 'Fatal. Bin range synchronization failed in parallel Wang-Landau method. This is a bug.'
            call fexit()
          end if
        else if (wld%dimensionality.eq.2) then
          allocate(gtmp2d(wld%nbins2d(1),wld%nbins2d(2),4))
          call wl_parallel_combine(gtmp2d=gtmp2d)
          if ((sum(wld%ghtot2d(wld%minb2d(1),:)).le.0).OR.(sum(wld%ghtot2d(wld%maxb2d(1),:)).le.0).OR.&
 &            (sum(wld%ghtot2d(:,wld%minb2d(2))).le.0).OR.(sum(wld%ghtot2d(:,wld%maxb2d(2))).le.0)) then
            write(ilog,*) 'Fatal. Bin range synchronization failed in parallel Wang-Landau method (2D). This is a bug.'
            call fexit()
          end if
        end if
      end if
      wld%fval = critv/dble(itmp)
    else if (wld%stepnum.gt.wld%buffer) then
      if ((wld%flevel.gt.0).AND.(wld%fval.le.(critv/dble(itmp)))) then
        wld%t1on = .true.
        wld%fval = critv/dble(itmp)
      end if
    end if
  else ! non-communicating
    if (wld%t1on.EQV..false.) then
      if (wld%stepnum.gt.wld%buffer) then
        if ((wld%flevel.gt.0).AND.(wld%fval.le.(critv/dble(wld%stepnum)))) then
          wld%t1on = .true.
          wld%fval = critv/dble(wld%stepnum)
        end if
      end if
    else
      wld%fval = critv/dble(wld%stepnum)
    end if
  end if
#else
  if (wld%t1on.EQV..false.) then
    if (wld%stepnum.gt.wld%buffer) then
      if ((wld%flevel.gt.0).AND.(wld%fval.le.(critv/dble(wld%stepnum)))) then
        wld%t1on = .true.
        wld%fval = critv/dble(wld%stepnum)
      end if
    end if
  else
    wld%fval = critv/dble(wld%stepnum)
  end if
#endif
!
! if in 1/t regime, output the histograms every 0.5M steps
  if ((debug_wanglandau.EQV..true.).AND.(mod(wld%stepnum,500000).eq.0).AND.(wld%t1on.EQV..true.)) then
    bname = '_TMP_WL_G'
#ifdef ENABLE_MPI
    if (use_MPIAVG.EQV..true.) then
      if (myrank.eq.0) then
!        this is not optimal in performance, but temporarily write sum into wld%gh
        if (wld%dimensionality.eq.1) then
          wld%gh(:) = gtmp(:,3)
        else if (wld%dimensionality.eq.2) then
          wld%gh2d(:,:) = gtmp2d(:,:,3)
        end if
      end if
    end if
#endif
    call wl_write_temphists(bname,wld%stepnum/500000)
!   reset histo
    if (wld%dimensionality.eq.1) then
      wld%gh(:) = 0
    else if (wld%dimensionality.eq.2) then
      wld%gh2d(:,:) = 0
    end if
  end if
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (wld%dimensionality.eq.1) then
      if (allocated(gtmp).EQV..true.) deallocate(gtmp)
    else if (wld%dimensionality.eq.2) then
      if (allocated(gtmp2d).EQV..true.) deallocate(gtmp2d)
    end if
  end if
#endif
!
!  write(*,*) 'S=',istep,' O=',g_coord_current,'N=',g_coord_proposed,moveok
!
end subroutine wl_update_histos
!
!-------------------------------------------------------------------------------------------------------
!
subroutine wl_get_binpos(inrange,ecurr,diff,idx_old,idx_new,idx_old2,idx_new2)
!
  use wl
  use iounit
  use torsn
  use molecule
!
  implicit none
!
  RTYPE ecurr,diff,g_coord_current(2),g_coord_proposed(2),zabsa(2)
  integer idx_new,idx_old,idxv(2),exby,exbv(2),idx_old2,idx_new2
  logical inrange
!
  inrange = .true.
  if ((wld%wl_mode.eq.1).OR.((wld%wl_mode.eq.3).AND.(wld%dimensionality.eq.1))) then
    if ((ecurr+diff).gt.(wld%g_max+wld%g_binsz)) then
      idx_new = wld%nbins + 1
    else
      idx_new = floor((ecurr + diff - wld%minv)/wld%g_binsz)+1 
    end if
    idx_old = floor((ecurr  - wld%minv)/wld%g_binsz)+1
  else if ((wld%wl_mode.eq.3).AND.(wld%dimensionality.eq.2)) then
    if ((ecurr+diff).gt.(wld%g_max2d(1)+wld%g_binsz2d(1))) then
      idx_new = wld%nbins2d(1) + 1
    else
      idx_new = floor((ecurr + diff - wld%minv2d(1))/wld%g_binsz2d(1))+1 
    end if
    idx_old = floor((ecurr  - wld%minv2d(1))/wld%g_binsz2d(1))+1
  end if
!
  if ((wld%wl_mode.eq.2).AND.(wld%dimensionality.eq.1)) then
!   the "current" RC value is what was already computed on the previous step
    if (wld%wl_rc(1).eq.1) then
      g_coord_current(1) = rgv(wld%wl_mol(1))
      call rg_internal(wld%wl_mol(1),g_coord_proposed(1))
    else if (wld%wl_rc(1).eq.2) then
      g_coord_current(1) = z_alpha(wld%wl_mol(1))
      call z_secondary(wld%wl_mol(1),g_coord_proposed(1),zabsa(2))
    else if (wld%wl_rc(1).eq.3) then
      g_coord_current(1) = z_beta(wld%wl_mol(1))
      call z_secondary(wld%wl_mol(1),zabsa(1),g_coord_proposed(1))
    end if
    if (g_coord_proposed(1).gt.(wld%g_max+wld%g_binsz)) then
      idx_new = wld%nbins + 1
    else
      idx_new = floor((g_coord_proposed(1) - wld%minv)/wld%g_binsz)+1 
    end if
    idx_old = floor((g_coord_current(1) - wld%minv)/wld%g_binsz)+1
  else if ((wld%wl_mode.eq.2).AND.(wld%dimensionality.eq.2)) then
!   the "current" RC value is what was already computed on the previous step
    if (wld%wl_rc(1).eq.1) then
      g_coord_current(1) = rgv(wld%wl_mol(1))
      call rg_internal(wld%wl_mol(1),g_coord_proposed(1))
    else if (wld%wl_rc(1).eq.2) then
      g_coord_current(1) = z_alpha(wld%wl_mol(1))
      call z_secondary(wld%wl_mol(1),g_coord_proposed(1),g_coord_proposed(2))
    else if (wld%wl_rc(1).eq.3) then
      g_coord_current(1) = z_beta(wld%wl_mol(1))
      call z_secondary(wld%wl_mol(1),g_coord_proposed(2),g_coord_proposed(1))
    end if
    if (g_coord_proposed(1).gt.(wld%g_max2d(1)+wld%g_binsz2d(1))) then
      idx_new = wld%nbins2d(1) + 1
    else
      idx_new = floor((g_coord_proposed(1) - wld%minv2d(1))/wld%g_binsz2d(1))+1 
    end if
    idx_old = floor((g_coord_current(1) - wld%minv2d(1))/wld%g_binsz2d(1))+1
  end if
  if (wld%dimensionality.eq.2) then ! second D is always RC
!   the "current" RC value is what was already computed on the previous step
    if (wld%wl_rc(2).eq.1) then
      g_coord_current(2) = rgv(wld%wl_mol(2))
      call rg_internal(wld%wl_mol(2),g_coord_proposed(2))
    else
!     avoid calling z_secondary again
      if ((wld%wl_mode.eq.2).AND.(wld%wl_mol(1).eq.wld%wl_mol(2)).AND.((wld%wl_rc(1).eq.2).OR.(wld%wl_rc(1).eq.3))) then
        if (wld%wl_rc(2).eq.2) then
          g_coord_current(2) = z_alpha(wld%wl_mol(2))
        else if (wld%wl_rc(2).eq.3) then
          g_coord_current(2) = z_beta(wld%wl_mol(2))
        end if
      else
        if (wld%wl_rc(2).eq.2) then
          g_coord_current = z_alpha(wld%wl_mol(2))
          call z_secondary(wld%wl_mol(2),g_coord_proposed(2),zabsa(2))
       else if (wld%wl_rc(2).eq.3) then
          g_coord_current = z_beta(wld%wl_mol(2))
          call z_secondary(wld%wl_mol(2),zabsa(1),g_coord_proposed(2))
        end if
      end if
    end if
    if (g_coord_proposed(2).gt.(wld%g_max2d(2)+wld%g_binsz2d(2))) then
      idx_new2 = wld%nbins2d(2) + 1
    else
      idx_new2 = floor((g_coord_proposed(2) - wld%minv2d(2))/wld%g_binsz2d(2))+1 
    end if
    idx_old2 = floor((g_coord_current(2)  - wld%minv2d(2))/wld%g_binsz2d(2))+1
  end if
!
! the rest are simple sanity checks: WL must start in the correct g bin range.  
!   terminate if this is not true
  if (wld%dimensionality.eq.1) then
    if (idx_old.gt.wld%nbins) then
      if ((wld%exrule.eq.3).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        call wl_extend_histos(idx_old,exby) 
!     for energy exit fatally, if upper limit exceeded on first step
      else if ((wld%stepnum.eq.1).AND.((wld%wl_mode.eq.1).OR.(wld%wl_mode.eq.3))) then
        write(ilog,*) 'Wang-Landau simulation is not properly equilibrated for data collection. Change upper limit &
 &on energy histogram criterion or increase equilibration time.'
        call fexit()
      else
!       use last (currently spanned) bin instead (wrong and ugly)
        if ((wld%flevel.gt.0).OR.(wld%freeze.EQV..true.)) then
          idx_old = wld%maxb
        else
          idx_old = wld%nbins
        end if
      end if
    end if
    if (idx_old.lt.1) then
      if ((wld%exrule.ge.2).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        call wl_extend_histos(idx_old,exby) 
        idx_new = idx_new + exby
        idx_old = idx_old + exby
      else
!       use first (currently spanned) bin instead (wrong and ugly)
        if ((wld%flevel.gt.0).OR.(wld%freeze.EQV..true.)) then
          idx_old = wld%minb
        else
          idx_old = 1
        end if
      end if
    end if
    if (idx_new.gt.wld%nbins) then
      if ((wld%exrule.eq.3).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        call wl_extend_histos(idx_new,exby)
      else
        inrange = .false.
        wld%overflow = wld%overflow + 1
      end if
    else if ((idx_new.gt.wld%maxb).AND.((wld%flevel.gt.0).AND.(wld%freeze.EQV..true.))) then
      inrange = .false.
      wld%overflow = wld%overflow + 1
    end if
    if (idx_new.lt.1) then
      if ((wld%exrule.ge.2).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        call wl_extend_histos(idx_new,exby) 
        idx_new = idx_new + exby
        idx_old = idx_old + exby
      else
        inrange = .false.
        wld%underflow = wld%underflow + 1
      end if
    else if ((idx_new.lt.wld%minb).AND.((wld%flevel.gt.0).AND.(wld%freeze.EQV..true.))) then
      inrange = .false.
      wld%underflow = wld%underflow + 1
    end if
  else if (wld%dimensionality.eq.2) then
    if ((idx_old.gt.wld%nbins2d(1)).OR.(idx_old2.gt.wld%nbins2d(2))) then
      if ((wld%exrule.eq.3).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        idxv(1) = idx_old
        idxv(2) = idx_old2
        call wl_extend_histos2d(idxv,exbv)
!     for energy exit fatally, if upper limit exceeded on first step
      else if ((wld%wl_mode.eq.3).AND.(idx_old.gt.wld%nbins2d(1)).AND.(wld%stepnum.eq.1)) then
        write(ilog,*) '2D Wang-Landau simulation is not properly equilibrated for data collection. Change upper limit &
 &on energy histogram criterion or increase equilibration time.'
        call fexit()
      else
!       use last (currently spanned) bin instead (wrong and ugly)
        if ((wld%flevel.gt.0).OR.(wld%freeze.EQV..true.)) then
          if (idx_old.gt.wld%nbins2d(1)) idx_old = wld%maxb2d(1)
          if (idx_old2.gt.wld%nbins2d(2)) idx_old2 = wld%maxb2d(2)
        else
          if (idx_old.gt.wld%nbins2d(1)) idx_old = wld%nbins2d(1)
          if (idx_old2.gt.wld%nbins2d(2)) idx_old2 = wld%nbins2d(2)
        end if
      end if
    end if
    if ((idx_old.lt.1).OR.(idx_old2.lt.1)) then
      if ((wld%exrule.ge.2).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        idxv(1) = idx_old
        idxv(2) = idx_old2
        call wl_extend_histos2d(idxv,exbv) 
        idx_new = idx_new + exbv(1)
        idx_old = idx_old + exbv(1)
        idx_new2 = idx_new2 + exbv(2)
        idx_old2 = idx_old2 + exbv(2)
      else
!       use first (currently spanned) bin instead (wrong and ugly)
        if ((wld%flevel.gt.0).OR.(wld%freeze.EQV..true.)) then
          if (idx_old.lt.1) idx_old = wld%minb2d(1)
          if (idx_old2.lt.1) idx_old2 = wld%minb2d(2)
        else
          if (idx_old.lt.1) idx_old = 1
          if (idx_old2.lt.1) idx_old2 = 1
        end if
      end if
    end if
    if ((idx_new.gt.wld%nbins2d(1)).OR.(idx_new2.gt.wld%nbins2d(2))) then
      if ((wld%exrule.eq.3).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        idxv(1) = idx_new
        idxv(2) = idx_new2
        call wl_extend_histos2d(idxv,exbv)
      else
        inrange = .false.
        wld%overflow = wld%overflow + 1
      end if
    else if (((idx_new.gt.wld%maxb2d(1)).OR.(idx_new2.gt.wld%maxb2d(2))).AND.&
 &           ((wld%flevel.gt.0).AND.(wld%freeze.EQV..true.))) then
      inrange = .false.
      wld%overflow = wld%overflow + 1
    end if
    if ((idx_new.lt.1).OR.(idx_new2.lt.1)) then
      if ((wld%exrule.ge.2).AND.((wld%flevel.eq.0).OR.(wld%freeze.EQV..false.))) then
        idxv(1) = idx_new
        idxv(2) = idx_new2
        call wl_extend_histos2d(idxv,exbv)
        idx_new = idx_new + exbv(1)
        idx_old = idx_old + exbv(1)
        idx_new2 = idx_new2 + exbv(2)
        idx_old2 = idx_old2 + exbv(2)
      else
        inrange = .false.
        wld%underflow = wld%underflow + 1
      end if
    else if (((idx_new.lt.wld%minb2d(1)).OR.(idx_new2.lt.wld%minb2d(2))).AND.&
 &           ((wld%flevel.gt.0).AND.(wld%freeze.EQV..true.))) then
      inrange = .false.
      wld%underflow = wld%underflow + 1
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
! this helper keeps the global RC arrays for the relevant molecules up-to-date (called in mcstat())
!
subroutine wl_keep_rc_updated()
!
  use wl
  use molecule
  use torsn
!
  implicit none
!
  if ((wld%dimensionality.eq.1).AND.(wld%wl_mode.eq.2)) then
    if (wld%wl_rc(1).eq.1) then
      call rg_internal(wld%wl_mol(1),rgv(wld%wl_mol(1)))
    else if ((wld%wl_rc(1).eq.2).OR.(wld%wl_rc(1).eq.3)) then
      call z_secondary(wld%wl_mol(1),z_alpha(wld%wl_mol(1)),z_beta(wld%wl_mol(1)))
    end if
  else if (wld%dimensionality.eq.2) then
    if (wld%wl_rc(2).eq.1) then
      call rg_internal(wld%wl_mol(2),rgv(wld%wl_mol(2)))
    else if ((wld%wl_rc(2).eq.2).OR.(wld%wl_rc(2).eq.3)) then
      call z_secondary(wld%wl_mol(2),z_alpha(wld%wl_mol(2)),z_beta(wld%wl_mol(2)))
    end if
  end if
!
end
!
!--------------------------------------------------------------------------
!
! an independent vectorized Rg calculation
!
subroutine rg_internal(imol,g_coord_proposed)
!
  use iounit
  use molecule
  use atoms
!
  implicit none
!
  integer i,j,imol
  RTYPE rgtentmp(3,3),xc,yc,zc,g_coord_proposed,fatm
!
  xc = 0.0d0
  yc = 0.0d0
  zc = 0.0d0
  rgtentmp(:,:) = 0.0
  fatm = dble(atmol(imol,2)-atmol(imol,1)+1)
!
  i=atmol(imol,1)
  j=atmol(imol,2)
  xc = sum(x(i:j))
  yc = sum(y(i:j))
  zc = sum(z(i:j))
  xc = xc / fatm
  yc = yc / fatm
  zc = zc / fatm
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
!    write(ilog,*) 'Error: Rg of single atom in WL'
    return
  end if
! update gyration tensor: calc in same way as
! external routines to maintain numerical accuracy  
  rgtentmp(1,1) = sum((x(i:j)-xc)**2)
  rgtentmp(1,2) = sum((x(i:j)-xc)*(y(i:j)-yc))
  rgtentmp(1,3) = sum((x(i:j)-xc)*(z(i:j)-zc))
  rgtentmp(2,2) = sum((y(i:j)-yc)**2)
  rgtentmp(2,3) = sum((y(i:j)-yc)*(z(i:j)-zc))
  rgtentmp(3,3) = sum((z(i:j)-zc)**2)
  rgtentmp(2,1) = rgtentmp(1,2)
  rgtentmp(3,1) = rgtentmp(1,3)
  rgtentmp(3,2) = rgtentmp(2,3)
  rgtentmp(:,:) = rgtentmp(:,:)/fatm
  g_coord_proposed = sqrt(rgtentmp(1,1)+rgtentmp(2,2)+rgtentmp(3,3))
!
end subroutine rg_internal
!
!---------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_MPI
!
! this routine sums up gh and ghbu over all nodes, increments ghtot with the non-self ghbu, and resets
! ghbu to zero for the next increment cycle
!
subroutine wl_parallel_combine(gtmp,gtmp2d)
!
  use wl
  use iounit
  use mpistuff
  use mpi
!
  implicit none
!
  integer ierr,msgtag(4),aone,msgsz(3)
  logical bcast
! the assumed size for gtmp is (1:wld%nbins,4)
! the assumed size for gtmp2d is (1:wld%nbins2d(1),1:wld%nbins2d(2),4)
  integer, OPTIONAL, ALLOCATABLE, INTENT(IN OUT):: gtmp(:,:),gtmp2d(:,:,:)
!
  msgtag(1) = 2000+mod(4*wld%stepnum/wld%gh_flatcheck_freq+30000,30000)
  msgtag(2) = msgtag(1) + 1
  msgtag(3) = msgtag(1) + 2
  msgtag(4) = msgtag(1) + 3
  bcast = .true. ! emulates MPI_ALLREDUCE
  aone = 1 ! 1 is summation
  msgsz(3) = 4
!
  if (wld%dimensionality.eq.1) then
    if (present(gtmp).EQV..false.) then
      write(ilog,*) 'Fatal. Called wl_parallel_combine(...) with optional argument missing for a case when it &
 &is required. This is a bug.'
      call fexit()
    end if
    msgsz(2) = wld%nbins
    gtmp(:,2) = wld%ghbu(:)
    gtmp(:,1) = wld%gh(:)
    call MPI_ALLINTS2D(msgtag(1:4),aone,ierr,bcast,msgsz(2:3),gtmp,use_MPIcolls)
    if (ierr.ne.0) then
      write(ilog,*) 'Fatal. MPI_ALLINTS2D(...) exited with non-zero error in wl_parallel_combine(...) on &
 &node ',myrank+1,'.'
      call fexit()
    end if
!   increment ghtot and g by total difference across all replicas
    wld%ghtot(:) = wld%ghtot(:) - wld%ghbu(:) + gtmp(:,4)
    wld%g(:) = wld%g(:) + (gtmp(:,4)-wld%ghbu(:))*wld%fval
    wld%ghbu(:) = 0
  else if (wld%dimensionality.eq.2) then
    if (present(gtmp2d).EQV..false.) then
      write(ilog,*) 'Fatal. Called wl_parallel_combine(...) with optional argument missing for a case when it &
 &is required. This is a bug.'
      call fexit()
    end if
    msgsz(1:2) = wld%nbins2d(:)
    gtmp2d(:,:,2) = wld%ghbu2d(:,:)
    gtmp2d(:,:,1) = wld%gh2d(:,:)
    call MPI_ALLINTS3D(msgtag(1:4),aone,ierr,bcast,msgsz(1:3),gtmp2d,use_MPIcolls)
    if (ierr.ne.0) then
      write(ilog,*) 'Fatal. MPI_ALLINTS3D(...) exited with non-zero error in wl_parallel_combine(...) on &
 &node ',myrank+1,'.'
      call fexit()
    end if
!   increment ghtot and g by total difference across all replicas
    wld%ghtot2d(:,:) = wld%ghtot2d(:,:) - wld%ghbu2d(:,:) + gtmp2d(:,:,4)
    wld%g2d(:,:) = wld%g2d(:,:) + (gtmp2d(:,:,4)-wld%ghbu2d(:,:))*wld%fval
    wld%ghbu2d(:,:) = 0
  end if
!
end
!
!-------------------------------------------------------------------------------------------------
!
subroutine wl_parallel_syncbins(idx_new,idx_new2)
!
  use wl
  use iounit
  use mpistuff
  use mpi
!
  implicit none
!
  integer idx_new,ierr,ij,ik,ikv(2),idxv(2),idx_new2,msgsz,atwo,msgtag(4)
  RTYPE mima(16)
  logical bcast
!
  msgtag(1) = 1001+mod(4*wld%stepnum/wld%gh_flatcheck_freq+30000,30000)
  msgtag(2) = msgtag(1) + 1
  msgtag(3) = msgtag(1) + 2
  msgtag(4) = msgtag(1) + 3
  bcast = .true.
  atwo = 2 ! 2 is max of
!
! synchronize bin ranges
  if (wld%dimensionality.eq.1) then
    mima(1) = -wld%bctr(wld%minb)
    mima(2) = wld%bctr(wld%maxb)
    mima(3) = -wld%g_min
    mima(4) = wld%g_max
    msgsz = 8
    call MPI_ALLREALS1D(msgtag,atwo,ierr,bcast,msgsz,mima(1:8),use_MPIcolls)
    if (ierr.ne.0) then
      write(ilog,*) 'Fatal. MPI_ALLREALS1D(...) exited with non-zero error in wl_parallel_syncbins(...) on &
 &node ',myrank+1,'.'
      call fexit()
    end if
    if ((mima(8)-wld%g_max).gt.0.5*wld%g_binsz) then
      ij = wld%nbins + nint((mima(8)-wld%g_max)/wld%g_binsz)
      call wl_extend_histos(ij,ik)
    end if
    if ((wld%g_min+mima(7)).gt.0.5*wld%g_binsz) then
      ij = 1 - nint((wld%g_min+mima(7))/wld%g_binsz)
      call wl_extend_histos(ij,ik)
      idx_new = idx_new + ik
    end if
    wld%minb = ceiling((-mima(5) - wld%minv)/wld%g_binsz)
    wld%maxb = ceiling((mima(6) - wld%minv)/wld%g_binsz)
  else if (wld%dimensionality.eq.2) then
    mima(1) = -wld%bctr2d1(wld%minb2d(1))
    mima(2) = -wld%bctr2d2(wld%minb2d(2))
    mima(3) = wld%bctr2d1(wld%maxb2d(1))
    mima(4) = wld%bctr2d2(wld%maxb2d(2))
    mima(5:6) = -wld%g_min2d(:)
    mima(7:8) = wld%g_max2d(:)
    msgsz = 16
    call MPI_ALLREALS1D(msgtag,atwo,ierr,bcast,msgsz,mima(1:16),use_MPIcolls)
    if (ierr.ne.0) then
      write(ilog,*) 'Fatal. MPI_ALLREALS1D(...) exited with non-zero error in wl_parallel_syncbins(...) on &
 &node ',myrank+1,'.'
      call fexit()
    end if
    if (((mima(15)-wld%g_max2d(1)).gt.0.5*wld%g_binsz2d(1)).OR.((mima(16)-wld%g_max2d(2)).gt.0.5*wld%g_binsz2d(2))) then
      idxv(1:2) = wld%nbins2d(1:2) + nint((mima(15:16)-wld%g_max2d(1:2))/wld%g_binsz2d(1:2))
      call wl_extend_histos2d(idxv,ikv)
    end if
    if (((mima(13)+wld%g_min2d(1)).gt.0.5*wld%g_binsz2d(1)).OR.((mima(14)+wld%g_max2d(2)).gt.0.5*wld%g_binsz2d(2))) then
      idxv(1:2) = 1 - nint((wld%g_min2d(1:2)+mima(13:14))/wld%g_binsz2d(1:2))
      call wl_extend_histos2d(idxv,ikv)
      idx_new = idx_new + ikv(1)
      idx_new2 = idx_new2 + ikv(2)
    end if
    wld%minb2d(1:2) = ceiling((-mima(9:10)  - wld%minv2d(1:2))/wld%g_binsz2d(1:2))
    wld%maxb2d(1:2) = ceiling((mima(11:12)  - wld%minv2d(1:2))/wld%g_binsz2d(1:2))
  end if
!
end
!
#endif
!
!--------------------------------------------------------------------------
!
subroutine wl_extend_histos(idxn,exby)
!
  use wl
  use mpistuff
!
  implicit none
!
  integer exby,idxn,i,k
  RTYPE, ALLOCATABLE:: gh1(:,:)
  integer, ALLOCATABLE:: ih1(:,:)
!
  if (idxn.eq.0) then ! only 1 extra bin would be required
    exby = 10
  else if (idxn.gt.wld%nbins) then
    exby = 10*ceiling(dble(idxn-wld%nbins)/10.)
  else
    exby = 10*ceiling((abs(idxn)+1.0)/10.)
  end if
!
  allocate(gh1(2,wld%nbins))
  allocate(ih1(3,wld%nbins))
  gh1(1,:) = wld%g(:)
  gh1(2,:) = wld%bctr(:)
  ih1(1,:) = wld%gh(:)
  ih1(2,:) = wld%ghtot(:)
  ih1(3,:) = wld%ghbu(:)
  deallocate(wld%g)
  deallocate(wld%gh)
  deallocate(wld%ghbu)
  deallocate(wld%ghtot)
  deallocate(wld%bctr)
  wld%nbins = wld%nbins + exby
  allocate(wld%g(wld%nbins))
  allocate(wld%gh(wld%nbins))
  allocate(wld%ghbu(wld%nbins))
  allocate(wld%ghtot(wld%nbins))
  allocate(wld%bctr(wld%nbins))
  if (idxn.le.0) then
    wld%minv = wld%minv - exby*wld%g_binsz
    do i=1,exby
      wld%bctr(i) = wld%minv + (i-0.5)*wld%g_binsz
      wld%ghtot(i) = 0
      wld%gh(i) = 0
      wld%ghbu(i) = 0
      wld%g(i) = 0.0
    end do
    k = exby + 1
    wld%g(k:wld%nbins) = gh1(1,:)
    wld%bctr(k:wld%nbins) = gh1(2,:)
    wld%gh(k:wld%nbins) = ih1(1,:)
    wld%ghtot(k:wld%nbins) = ih1(2,:)
    wld%ghbu(k:wld%nbins) = ih1(3,:)
    wld%g_min = wld%bctr(1)
    wld%minb = wld%minb + exby
    wld%maxb = wld%maxb + exby
  else
    do i=wld%nbins-exby+1,wld%nbins
      wld%bctr(i) = wld%minv + (i-0.5)*wld%g_binsz
      wld%ghtot(i) = 0
      wld%gh(i) = 0
      wld%ghbu(i) = 0
      wld%g(i) = 0.0
    end do
    wld%g(1:(wld%nbins-exby)) = gh1(1,:)
    wld%bctr(1:(wld%nbins-exby)) = gh1(2,:)
    wld%gh(1:(wld%nbins-exby)) = ih1(1,:)
    wld%ghtot(1:(wld%nbins-exby)) = ih1(2,:)
    wld%ghbu(1:(wld%nbins-exby)) = ih1(3,:)
    wld%g_max = wld%bctr(wld%nbins)
  end if
  deallocate(gh1)
  deallocate(ih1)
!
end subroutine wl_extend_histos
!
!--------------------------------------------------------------------------
!
! we may have the case of one dimension overflowing and the other underflowing, but
! extension may not be allowed to both sides -> screen
! this routine is never called if wl_exrule is 1 and/or histograms are frozen
!
subroutine wl_extend_histos2d(idxn,exby)
!
  use wl
  use mpistuff
!
  implicit none
!
  integer exby(2),idxn(2),i
  RTYPE, ALLOCATABLE:: gh1(:,:),ghb2(:),ghb1(:)
  integer, ALLOCATABLE:: ih1(:,:,:)
!
  do i=1,2
    if ((idxn(i).ge.1).AND.(idxn(i).le.wld%nbins2d(i))) then
      exby(i) = 0
    else if (idxn(i).le.0) then
      if (idxn(i).eq.0) then
        exby(i) = 10
      else
        exby(i) = 10*ceiling((abs(idxn(i))+1.0)/10.)
      end if
    else if ((idxn(i).gt.wld%nbins2d(i)).AND.(wld%exrule.eq.3)) then
      exby(i) = 10*ceiling(dble(idxn(i)-wld%nbins2d(i))/10.)
    else
      exby(i) = 0
    end if
  end do
!
  allocate(gh1(wld%nbins2d(1),wld%nbins2d(2)))
  if (exby(1).gt.0) allocate(ghb1(wld%nbins2d(1)))
  if (exby(2).gt.0) allocate(ghb2(wld%nbins2d(2)))
  allocate(ih1(3,wld%nbins2d(1),wld%nbins2d(2)))
  gh1(:,:) = wld%g2d(:,:)
  if (exby(1).gt.0) ghb1(:) = wld%bctr2d1(:)
  if (exby(2).gt.0) ghb2(:) = wld%bctr2d2(:)
  ih1(1,:,:) = wld%gh2d(:,:)
  ih1(2,:,:) = wld%ghtot2d(:,:)
  ih1(3,:,:) = wld%ghbu2d(:,:)
  deallocate(wld%g2d)
  deallocate(wld%gh2d)
  deallocate(wld%ghbu2d)
  deallocate(wld%ghtot2d)
  if (exby(1).gt.0) deallocate(wld%bctr2d1)
  if (exby(2).gt.0) deallocate(wld%bctr2d2)
  wld%nbins2d(:) = wld%nbins2d(:) + exby(:)
  allocate(wld%g2d(wld%nbins2d(1),wld%nbins2d(2)))
  allocate(wld%gh2d(wld%nbins2d(1),wld%nbins2d(2)))
  allocate(wld%ghbu2d(wld%nbins2d(1),wld%nbins2d(2)))
  allocate(wld%ghtot2d(wld%nbins2d(1),wld%nbins2d(2)))
  if (exby(1).gt.0) allocate(wld%bctr2d1(wld%nbins2d(1)))
  if (exby(2).gt.0) allocate(wld%bctr2d2(wld%nbins2d(2)))
!
  if (exby(1).gt.0) then
    if (idxn(1).gt.0) then
      wld%ghtot2d((wld%nbins2d(1)-exby(1)+1):wld%nbins2d(1),:) = 0
      wld%gh2d((wld%nbins2d(1)-exby(1)+1):wld%nbins2d(1),:) = 0
      wld%ghbu2d((wld%nbins2d(1)-exby(1)+1):wld%nbins2d(1),:) = 0
      wld%g2d((wld%nbins2d(1)-exby(1)+1):wld%nbins2d(1),:) = 0.0
      do i=(wld%nbins2d(1)-exby(1)+1),wld%nbins2d(1)
        wld%bctr2d1(i) = wld%minv2d(1) + (i-0.5)*wld%g_binsz2d(1)
      end do
      wld%bctr2d1(1:(wld%nbins2d(1)-exby(1))) = ghb1(:)
      wld%g_max2d(1) =  wld%bctr2d1(wld%nbins2d(1))
    else
      wld%minv2d(1) = wld%minv2d(1) - exby(1)*wld%g_binsz2d(1)
      wld%ghtot2d(1:exby(1),:) = 0
      wld%gh2d(1:exby(1),:) = 0
      wld%ghbu2d(1:exby(1),:) = 0
      wld%g2d(1:exby(1),:) = 0.0
      do i=1,exby(1)
        wld%bctr2d1(i) = wld%minv2d(1) + (i-0.5)*wld%g_binsz2d(1)
      end do
      wld%bctr2d1((exby(1)+1):wld%nbins2d(1)) = ghb1(:)
      wld%g_min2d(1) = wld%bctr2d1(1)
      wld%minb2d(1) = wld%minb2d(1) + exby(1)
      wld%maxb2d(1) = wld%maxb2d(1) + exby(1)
    end if
  end if
  if (exby(2).gt.0) then
    if (idxn(2).gt.0) then
      wld%ghtot2d(:,(wld%nbins2d(2)-exby(2)+1):wld%nbins2d(2)) = 0
      wld%gh2d(:,(wld%nbins2d(2)-exby(2)+1):wld%nbins2d(2)) = 0
      wld%ghbu2d(:,(wld%nbins2d(2)-exby(2)+1):wld%nbins2d(2)) = 0
      wld%g2d(:,(wld%nbins2d(2)-exby(2)+1):wld%nbins2d(2)) = 0.0
      do i=(wld%nbins2d(2)-exby(2)+1),wld%nbins2d(2)
        wld%bctr2d2(i) = wld%minv2d(2) + (i-0.5)*wld%g_binsz2d(2)
      end do
      wld%bctr2d2(1:(wld%nbins2d(2)-exby(2))) = ghb2(:)
      wld%g_max2d(2) =  wld%bctr2d2(wld%nbins2d(2))
    else
      wld%minv2d(2) = wld%minv2d(2) - exby(2)*wld%g_binsz2d(2)
      wld%ghtot2d(:,1:exby(2)) = 0
      wld%gh2d(:,1:exby(2)) = 0
      wld%ghbu2d(:,1:exby(2)) = 0
      wld%g2d(:,1:exby(2)) = 0.0
      do i=1,exby(2)
        wld%bctr2d2(i) = wld%minv2d(2) + (i-0.5)*wld%g_binsz2d(2)
      end do
      wld%bctr2d2((exby(2)+1):wld%nbins2d(2)) = ghb2(:)
      wld%g_min2d(2) = wld%bctr2d2(1)
      wld%minb2d(2) = wld%minb2d(2) + exby(2)
      wld%maxb2d(2) = wld%maxb2d(2) + exby(2)
    end if
  end if
  if ((idxn(1).gt.0).AND.(idxn(2).gt.0)) then
    wld%g2d(1:(wld%nbins2d(1)-exby(1)),1:(wld%nbins2d(2)-exby(2))) = gh1(:,:)
    wld%gh2d(1:(wld%nbins2d(1)-exby(1)),1:(wld%nbins2d(2)-exby(2))) = ih1(1,:,:)
    wld%ghtot2d(1:(wld%nbins2d(1)-exby(1)),1:(wld%nbins2d(2)-exby(2))) = ih1(2,:,:)
    wld%ghbu2d(1:(wld%nbins2d(1)-exby(1)),1:(wld%nbins2d(2)-exby(2))) = ih1(3,:,:)
  else if ((idxn(1).le.0).AND.(idxn(2).le.0)) then
    wld%g2d((exby(1)+1):wld%nbins2d(1),(exby(2)+1):wld%nbins2d(2)) = gh1(:,:)
    wld%gh2d((exby(1)+1):wld%nbins2d(1),(exby(2)+1):wld%nbins2d(2)) = ih1(1,:,:)
    wld%ghtot2d((exby(1)+1):wld%nbins2d(1),(exby(2)+1):wld%nbins2d(2)) = ih1(2,:,:)
    wld%ghbu2d((exby(1)+1):wld%nbins2d(1),(exby(2)+1):wld%nbins2d(2)) = ih1(3,:,:)
  else if ((idxn(1).gt.0).AND.(idxn(2).le.0)) then
    wld%g2d(1:(wld%nbins2d(1)-exby(1)),(exby(2)+1):wld%nbins2d(2)) = gh1(:,:)
    wld%gh2d(1:(wld%nbins2d(1)-exby(1)),(exby(2)+1):wld%nbins2d(2)) = ih1(1,:,:)
    wld%ghtot2d(1:(wld%nbins2d(1)-exby(1)),(exby(2)+1):wld%nbins2d(2)) = ih1(2,:,:)
    wld%ghbu2d(1:(wld%nbins2d(1)-exby(1)),(exby(2)+1):wld%nbins2d(2)) = ih1(3,:,:)
  else if ((idxn(1).le.0).AND.(idxn(2).gt.0)) then
    wld%g2d((exby(1)+1):wld%nbins2d(1),1:(wld%nbins2d(2)-exby(2))) = gh1(:,:)
    wld%gh2d((exby(1)+1):wld%nbins2d(1),1:(wld%nbins2d(2)-exby(2))) = ih1(1,:,:)
    wld%ghtot2d((exby(1)+1):wld%nbins2d(1),1:(wld%nbins2d(2)-exby(2))) = ih1(2,:,:)
    wld%ghbu2d((exby(1)+1):wld%nbins2d(1),1:(wld%nbins2d(2)-exby(2))) = ih1(3,:,:)
  end if
  deallocate(gh1)
  deallocate(ih1)
  if (exby(1).gt.0) deallocate(ghb1)
  if (exby(2).gt.0) deallocate(ghb2)
!
end subroutine wl_extend_histos2d
!
!-----------------------------------------------------------------------------------------
!
! simply resize histograms, losing all content and all information in bctr or bctr2d
!
subroutine wl_realloc_histos(s1,s2)
!
  use wl
  use iounit
!
  implicit none
!
  integer s1,s2
!
  if (wld%dimensionality.eq.1) then
    if (s1.le.1) then
      write(ilog,*) 'Fatal. Called wl_realloc_histos(...) with bad size request. This could be &
 &caused by bad (restart) file input, or be indicative of a bug.'
      call fexit()
    end if
    if (allocated(wld%gh).EQV..true.) deallocate(wld%gh)
    if (allocated(wld%ghbu).EQV..true.) deallocate(wld%ghbu)
    if (allocated(wld%ghtot).EQV..true.) deallocate(wld%ghtot)
    if (allocated(wld%g).EQV..true.) deallocate(wld%g)
    if (allocated(wld%bctr).EQV..true.) deallocate(wld%bctr)
    wld%nbins = s1
    allocate(wld%gh(wld%nbins))
    allocate(wld%ghbu(wld%nbins))
    allocate(wld%ghtot(wld%nbins))
    allocate(wld%g(wld%nbins))
    allocate(wld%bctr(wld%nbins))
    wld%gh(:) = 0
    wld%ghbu(:) = 0
    wld%ghtot(:) = 0
    wld%g(:) = 0
  else if (wld%dimensionality.eq.2) then
    if ((s1.le.1).OR.(s2.le.1)) then
      write(ilog,*) 'Fatal. Called wl_realloc_histos(...) with bad size request(s). This could be &
 &caused by bad (restart) file input, or be indicative of a bug.'
      call fexit()
    end if
    if (allocated(wld%gh2d).EQV..true.) deallocate(wld%gh2d)
    if (allocated(wld%ghbu2d).EQV..true.) deallocate(wld%ghbu2d)
    if (allocated(wld%ghtot2d).EQV..true.) deallocate(wld%ghtot2d)
    if (allocated(wld%g2d).EQV..true.) deallocate(wld%g2d)
    if (allocated(wld%bctr2d1).EQV..true.) deallocate(wld%bctr2d1)
    if (allocated(wld%bctr2d2).EQV..true.) deallocate(wld%bctr2d2)
    wld%nbins2d(1) = s1
    wld%nbins2d(2) = s2
    allocate(wld%gh2d(wld%nbins2d(1),wld%nbins2d(2)))
    allocate(wld%ghbu2d(wld%nbins2d(1),wld%nbins2d(2)))
    allocate(wld%ghtot2d(wld%nbins2d(1),wld%nbins2d(2)))
    allocate(wld%g2d(wld%nbins2d(1),wld%nbins2d(2)))
    allocate(wld%bctr2d1(wld%nbins2d(1)))
    allocate(wld%bctr2d2(wld%nbins2d(2)))
    wld%gh2d(:,:) = 0
    wld%ghbu2d(:,:) = 0
    wld%ghtot2d(:,:) = 0
    wld%g2d(:,:) = 0
  end if
!
end
!
!------------------------------------------------------------------------------------
!
subroutine wl_write_temphists(bname,danum)
!
  use wl
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer i,t1,t2,slen,itmp,danum
  character(MAXSTRLEN) bname,fn,fn2
  character(2) xpont
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
!
  call strlims(bname,t1,t2)
!
! write temp hists to track progress (these recycle due to two digits used)
  slen = 2
  call int2str(danum,xpont,slen)
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (myrank.ne.0) then
      return
    end if
    fn =   bname(t1:t2)//'_'//xpont(1:slen)
    fn2 =  bname(t1:t2)//'H_'//xpont(1:slen)
  else if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn =  'N_'//nod(1:tl)//bname(t1:t2)//'_'//xpont(1:slen)
    fn2 = 'N_'//nod(1:tl)//bname(t1:t2)//'H_'//xpont(1:slen)
  end if
#else
  fn =   bname(t1:t2)//'_'//xpont(1:slen)
  fn2 =  bname(t1:t2)//'H_'//xpont(1:slen)
#endif
  if (wld%dimensionality.eq.1) then
    call delete_then_openfile(itmp, fn)
    do i=wld%minb,wld%maxb
      write(itmp,*) wld%bctr(i),'  ',wld%g(i)
    end do
    call close_filehandle(itmp)
!
!   gh
    call delete_then_openfile(itmp, fn2)
    do i=wld%minb,wld%maxb
      write(itmp,*) wld%bctr(i),'  ',wld%gh(i)
    end do
    call close_filehandle(itmp)
  end if
!
end subroutine wl_write_temphists
!
!--------------------------------------------------------------------------
!
! Updates F using the first part of the algorithm
!
subroutine wl_update_convergence(istep)
!
  use wl
  use iounit
#ifdef ENABLE_MPI
  use mpistuff
  use interfaces
#endif
!
  implicit none
!
  integer istep,i,j
  character(MAXSTRLEN) bname
  RTYPE vtimes,minbinv
#ifdef ENABLE_MPI
  integer, ALLOCATABLE:: gtmp(:,:)
  integer, ALLOCATABLE:: gtmp2d(:,:,:)
#endif
!
! how many times do we need to visit each bin ?
  if (wld%hvmode.eq.1) then
    vtimes = 1.0d0/sqrt(wld%fval)
  else if (wld%hvmode.eq.2) then
    vtimes = 1
  else
    write(ilog,*) 'Fatal. Encountered an unsupported mode for WL convergence updates. This is a bug.'
    call fexit()
  end if
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (wld%dimensionality.eq.1) then
      allocate(gtmp(wld%nbins,4))
      call wl_parallel_combine(gtmp=gtmp)
    else if (wld%dimensionality.eq.2) then
      allocate(gtmp2d(wld%nbins2d(1),wld%nbins2d(2),4))
      call wl_parallel_combine(gtmp2d=gtmp2d)
    end if
  end if
#endif
!
! If each bin is visited enough, update the f value
! because the histogram is potentially too large, we discard consecutive, empty bins to either side
! in 2D this means only bins that are bracketed by finite counts in at least one of the dimensions are counted
  if (wld%dimensionality.eq.2) then
#ifdef ENABLE_MPI 
    minbinv = HUGE(minbinv)
    do i=1,wld%nbins2d(1)
      do j=1,wld%nbins2d(2)
        if (gtmp2d(i,j,3).lt.minbinv) then
          if (gtmp2d(i,j,3).le.0.0) then
            if ((i.gt.1).AND.(i.lt.wld%nbins2d(1)).AND.(j.gt.1).AND.(j.lt.wld%nbins2d(2))) then
              if (((sum(gtmp2d(1:i,j,3)).gt.0.0).AND.(sum(gtmp2d((i+1):wld%nbins2d(1),j,3)).gt.0.0)).OR.&
 &                ((sum(gtmp2d(i,1:j,3)).gt.0.0).AND.(sum(gtmp2d(i,(j+1):wld%nbins2d(2),3)).gt.0.0))) then
                minbinv = gtmp2d(i,j,3)
              end if
            end if
          else
            minbinv = gtmp2d(i,j,3)
          end if
        end if
        if (minbinv.le.0.0) exit
      end do
      if (minbinv.le.0.0) exit
    end do
  else if (wld%dimensionality.eq.1) then
    minbinv = minval(gtmp(wld%minb:wld%maxb,3))
  end if
#else
    minbinv = HUGE(minbinv)
    do i=1,wld%nbins2d(1)
      do j=1,wld%nbins2d(2)
        if (wld%gh2d(i,j).lt.minbinv) then
          if (wld%gh2d(i,j).le.0.0) then
            if ((i.gt.1).AND.(i.lt.wld%nbins2d(1)).AND.(j.gt.1).AND.(j.lt.wld%nbins2d(2))) then
              if (((sum(wld%gh2d(1:i,j)).gt.0.0).AND.(sum(wld%gh2d((i+1):wld%nbins2d(1),j)).gt.0.0)).OR.&
 &                ((sum(wld%gh2d(i,1:j)).gt.0.0).AND.(sum(wld%gh2d(i,(j+1):wld%nbins2d(2))).gt.0.0))) then
                minbinv = wld%gh2d(i,j)
              end if
            end if
          else
            minbinv = wld%gh2d(i,j)
          end if
        end if
        if (minbinv.le.0.0) exit
      end do
      if (minbinv.le.0.0) exit
    end do
  else if (wld%dimensionality.eq.1) then
    minbinv = minval(wld%gh(wld%minb:wld%maxb))
  end if
#endif
  if ((minbinv.ge.vtimes).AND.(wld%stepnum.gt.wld%buffer)) then
    wld%fval = wld%fval/2.0
    wld%flevel = wld%flevel + 1
    if (debug_wanglandau.EQV..true.) then
      write(ilog,*) '   '
      if ((wld%flevel.eq.1).OR.(wld%freeze.EQV..false.)) then
        if (wld%dimensionality.eq.1) then
          write(ilog,*) 'Minimum: ',wld%bctr(wld%minb)-0.5*wld%g_binsz
          write(ilog,*) 'Maximum: ',wld%bctr(wld%maxb)+0.5*wld%g_binsz
        else if (wld%dimensionality.eq.2) then
          write(ilog,*) 'Minimum in Dim 1: ',wld%bctr2d1(wld%minb2d(1))-0.5*wld%g_binsz2d(1)
          write(ilog,*) 'Maximum in Dim 1: ',wld%bctr2d1(wld%maxb2d(1))+0.5*wld%g_binsz2d(1)
          write(ilog,*) 'Minimum in Dim 2: ',wld%bctr2d2(wld%minb2d(2))-0.5*wld%g_binsz2d(2)
          write(ilog,*) 'Maximum in Dim 2: ',wld%bctr2d2(wld%maxb2d(2))+0.5*wld%g_binsz2d(2)
        end if
      end if
      write(ilog,*) 'F-level: ', wld%flevel
      write(ilog,*) 'New F-VALUE: ', wld%fval
      write(ilog,*) '1/t: ', dble(wld%nbins)/dble(wld%stepnum )
      if (wld%hvmode.eq.1) write(ilog,*) '# Per Bin: ', 1.0d0/sqrt(wld%fval)
      write(ilog,*) 'MC STEP: ',istep
      bname = '_TMP_WL_G'
#ifdef ENABLE_MPI
      if (use_MPIAVG.EQV..true.) then
        if (myrank.eq.0) then
!         this is not optimal in performance, but temporarily write sum into wld%gh
          if (wld%dimensionality.eq.1) then
            wld%gh(:) = gtmp(:,3)
          else if (wld%dimensionality.eq.2) then
            wld%gh2d(:,:) = gtmp2d(:,:,3)
          end if
        end if
      end if
#endif
      call wl_write_temphists(bname,wld%flevel-1)
    end if
!
!   Reset histo
    if (wld%dimensionality.eq.1) then
      wld%gh(:) = 0
    else if (wld%dimensionality.eq.2) then
      wld%gh2d(:,:) = 0
    end if
  else
#ifdef ENABLE_MPI
    if (use_MPIAVG.EQV..true.) then
      if (wld%dimensionality.eq.1) then
        deallocate(gtmp)
      else if (wld%dimensionality.eq.2) then
        deallocate(gtmp2d)
      end if
    end if
#endif
!   non-MPI do nothing
  end if
!
end subroutine wl_update_convergence
!
!--------------------------------------------------------------------------
!
subroutine wl_sim_finished
!
  use wl
  use iounit
#ifdef ENABLE_MPI
  use mpi
  use mpistuff
  use interfaces
#endif
!
  implicit none
!
  integer i,ij,ik,ijk(4)
  character(MAXSTRLEN) fn,fn2
  logical doprt
  integer, ALLOCATABLE:: veci(:)
  RTYPE, ALLOCATABLE:: vecr(:)
#ifdef ENABLE_MPI
  integer tl,k1,k2
  integer, ALLOCATABLE:: gtmp2d(:,:,:),gtmp(:,:)
  character(3) nod
#endif
!
 56   format(2x,a,i10)
 57   format(2x,a,g14.7)
 58   format(g24.14,1x,g25.15)
 59   format(g24.14,1x,i12)
 60   format(1000001(g25.15))
 61   format(' #                       ',1000000(g25.15))
 62   format(g24.14,1x,1000000(i12))
  ik = wld%stepnum-wld%stepnumbuf
  write(ilog,*)
  write(ilog,*) 'Wang Landau (WL) statistics:'
  write(ilog,56) 'TOTAL NUMBER OF ATTEMPTED WL MOVES     : ',ik
  write(ilog,57) 'FRACTION OF ACCEPTED WL MOVES          : ',dble(wld%accepted)/dble(ik)
!  write(ilog,57) 'Energy underflow fraction: ', dble(wld%underflow)/dble(wld%stepnum)
  if ((wld%exrule.ne.3).OR.(wld%freeze.EQV..true.)) then
    write(ilog,57) 'FRACTION OF MOVES EXCEEDING MAX. CRIT. : ',dble(wld%overflow)/dble(ik)
  end if
  if ((wld%exrule.eq.1).OR.(wld%freeze.EQV..true.)) then
    write(ilog,57) 'FRACTION OF MOVES BELOW MIN. CRIT.     : ',dble(wld%underflow)/dble(ik)
  end if
  write(ilog,56) 'NUMBER OF F-FACTOR REDUCTIONS          : ',wld%flevel
  write(ilog,*)
!
! if in MPIAVG, do one last round of syncing/combining data
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    k1 = 1 ! dummy
    k2 = 1 ! dummy
    call wl_parallel_syncbins(k1,k2)
    if (wld%dimensionality.eq.1) then
      allocate(gtmp(wld%nbins,4))
      call wl_parallel_combine(gtmp=gtmp)
    else if (wld%dimensionality.eq.2) then
      allocate(gtmp2d(wld%nbins2d(1),wld%nbins2d(2),4))
      call wl_parallel_combine(gtmp2d=gtmp2d)
      deallocate(gtmp2d)
    end if
  end if
#endif
!
! Get the effective # of bins
  doprt = .true.
  if (wld%dimensionality.eq.1) then
    if (sum(wld%ghtot).le.0) then
      write(ilog,*) 'Warning. Empty histogram of total counts in &
 &Wang-Landau sampling suggests corrupt calculation.'
      doprt = .false.
    end if
    i = wld%nbins
    do while (wld%ghtot(i).eq.0)
     i = i - 1
    end do
    ik = i
    i = 1
    do while (wld%ghtot(i).eq.0)
      i = i + 1
    end do
    ij = i
    if ((ij.ne.wld%minb).OR.(ik.ne.wld%maxb)) then
      write(ilog,*) 'Warning. Mismatch between visited range and possible range &
 &in Wang-Landau sampling suggests results that lack convergence.'
      write(*,*) wld%minb,wld%maxb,ij,ik
    end if
    ij = 1
  else if (wld%dimensionality.eq.2) then
    ijk(1) = 1
    do while (sum(wld%ghtot2d(ijk(1),:)).le.0) 
      ijk(1) = ijk(1) + 1
    end do
    ijk(2) = wld%nbins2d(1)
    do while (sum(wld%ghtot2d(ijk(2),:)).le.0) 
      ijk(2) = ijk(2) - 1
    end do
    ijk(3) = 1
    do while (sum(wld%ghtot2d(:,ijk(3))).le.0) 
      ijk(3) = ijk(3) + 1
    end do
    ijk(4) = wld%nbins2d(2)
    do while (sum(wld%ghtot2d(:,ijk(4))).le.0) 
      ijk(4) = ijk(4) - 1
    end do
    if ((ijk(1).ne.wld%minb2d(1)).OR.(ijk(2).ne.wld%maxb2d(1)).OR.&
 &      (ijk(3).ne.wld%minb2d(2)).OR.(ijk(4).ne.wld%maxb2d(2))) then
      write(ilog,*) 'Warning. Mismatch between visited range and possible range &
 &in 2D Wang-Landau sampling suggests results that lack convergence.'
    end if
  end if
!  
! Write final G to file
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (myrank.ne.0) then
      if (allocated(wld%g).EQV..true.) deallocate(wld%g)
      if (allocated(wld%gh).EQV..true.) deallocate(wld%gh)
      if (allocated(wld%ghtot).EQV..true.) deallocate(wld%ghtot)
      if (allocated(wld%ghbu).EQV..true.) deallocate(wld%ghbu)
      if (allocated(wld%bctr).EQV..true.) deallocate(wld%bctr)
      if (allocated(wld%g2d).EQV..true.) deallocate(wld%g2d)
      if (allocated(wld%gh2d).EQV..true.) deallocate(wld%gh2d)
      if (allocated(wld%ghtot2d).EQV..true.) deallocate(wld%ghtot2d)
      if (allocated(wld%ghbu2d).EQV..true.) deallocate(wld%ghbu2d)
      if (allocated(wld%bctr2d1).EQV..true.) deallocate(wld%bctr2d1)
      if (allocated(wld%bctr2d2).EQV..true.) deallocate(wld%bctr2d2)
      return
    end if
    fn =  'WANGLANDAU_G.dat'
    fn2 =  'WANGLANDAU_GH.dat'
  else if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn =  'N_'//nod(1:tl)//'_WANGLANDAU_G.dat'
    fn2 = 'N_'//nod(1:tl)//'_WANGLANDAU_GH.dat'
  end if
#else
  fn =  'WANGLANDAU_G.dat'
  fn2 = 'WANGLANDAU_GH.dat'
#endif 
  if (wld%dimensionality.eq.1) then
    call delete_then_openfile(wld%ig, fn)
    do i=ij,ik
      write(wld%ig,58) wld%bctr(i),wld%g(i)
    end do
    call close_filehandle(wld%ig)
!
!   Output G histogram (never reset) to a file 
    call delete_then_openfile(wld%ightot, fn2)
    do i=ij,ik
      write(wld%ightot,59) wld%bctr(i),wld%ghtot(i)
    end do
    call close_filehandle(wld%ightot)
  else if (wld%dimensionality.eq.2) then
    call delete_then_openfile(wld%ig, fn)
    write(wld%ig,61) wld%bctr2d2(ijk(3):ijk(4))
    allocate(vecr(ijk(4)-ijk(3)+1))
    do i=ijk(1),ijk(2)
      vecr = wld%g2d(i,ijk(3):ijk(4))
      write(wld%ig,60) wld%bctr2d1(i),vecr
    end do
    deallocate(vecr)
    call close_filehandle(wld%ig)
!
!   Output G histogram (never reset) to a file 
    call delete_then_openfile(wld%ightot, fn2)
    write(wld%ightot,61) wld%bctr2d2(ijk(3):ijk(4))
    allocate(veci(ijk(4)-ijk(3)+1))
    do i=ijk(1),ijk(2)
      veci = wld%ghtot2d(i,ijk(3):ijk(4))
      write(wld%ightot,62) wld%bctr2d1(i),veci
    end do
    deallocate(veci)
    call close_filehandle(wld%ightot)
  end if
!
  if (allocated(wld%g).EQV..true.) deallocate(wld%g)
  if (allocated(wld%gh).EQV..true.) deallocate(wld%gh)
  if (allocated(wld%ghtot).EQV..true.) deallocate(wld%ghtot)
  if (allocated(wld%ghbu).EQV..true.) deallocate(wld%ghbu)
  if (allocated(wld%bctr).EQV..true.) deallocate(wld%bctr)
  if (allocated(wld%g2d).EQV..true.) deallocate(wld%g2d)
  if (allocated(wld%gh2d).EQV..true.) deallocate(wld%gh2d)
  if (allocated(wld%ghtot2d).EQV..true.) deallocate(wld%ghtot2d)
  if (allocated(wld%ghbu2d).EQV..true.) deallocate(wld%ghbu2d)
  if (allocated(wld%bctr2d1).EQV..true.) deallocate(wld%bctr2d1)
  if (allocated(wld%bctr2d2).EQV..true.) deallocate(wld%bctr2d2)
!
end subroutine wl_sim_finished
!
!--------------------------------------------------------------------------
!
subroutine wl_parse_ginitfile(fpath)

  use wl
  use iounit
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  character(len=*) fpath
  character(MAXSTRLEN) str2
  integer iomess, i, k, st1, st2
  logical openfile,resetbins(2)
  RTYPE bd1, bd2
  RTYPE, ALLOCATABLE:: vecr(:)
!
  if (openfile(wld%igin,fpath).EQV..false.) then
    write(ilog,*) 'Warning. Initialization file for density histogram for Wang-Landau sampling (',fpath,') could &
 &not be opened. Histogram will remain initialized to zero.'
    return
  end if
!
 79 format(FORM_MAXSTRLEN)
!
  k = 0
  if (wld%dimensionality.eq.2) then
!   the first line in 2D case contains an integer with the number of columns first, followed by the bin centers
    read(wld%igin,*,iostat=iomess) str2
    call strlims_quiet(str2,st1,st2)
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Warning. Initialization file for density histogram for Wang-Landau sampling (',fpath,') is empty. &
 &Histogram will assume default settings and remain initialized to zero.'
      return
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing histogram initialization input for Wang-Landau sampling &
 &(got: ',fpath,').'
      call fexit()
    end if
    i = 1
    call extract_int(str2(st1:st2),wld%nbins2d(2),i)
  end if
!
  do while(.true.)
!   all other lines contain bin centers first, then actual value(s)
    read(wld%igin,79,iostat=iomess) str2
    call strlims_quiet(str2,st1,st2)
    if (st2.lt.(st1+2)) then
      write(ilog,*) 'Fatal. Incomplete or empty line while processing histogram initialization input for Wang-Landau sampling &
 &(got: ',fpath,'). Please remove all unnecessary lines.'
      call fexit()
    end if
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing histogram initialization input for Wang-Landau sampling &
 &(got: ',fpath,').'
      call fexit()
    end if
    k = k + 1
  end do
!
  if ((k.le.1).OR.((wld%dimensionality.eq.2).AND.(wld%nbins2d(2).le.1))) then
    write(ilog,*) 'Warning. Initialization file for density histogram for Wang-Landau sampling (',fpath,') is empty (at most &
 &one bin) or corrupt. Histogram will assume default settings and remain initialized to zero.'
    return
  end if
!
  st1 = k
  st2 = wld%nbins2d(2)
  call wl_realloc_histos(st1,st2)
!
  k = 0
  rewind(unit=wld%igin)
  resetbins(:) = .false.
!
  if (wld%dimensionality.eq.2) then
    read(wld%igin,*,iostat=iomess) st1,wld%bctr2d2(1:wld%nbins2d(2))
    bd2 = 1.0/(wld%bctr2d2(wld%nbins2d(2)) - wld%bctr2d2(1))
    do i=2,wld%nbins2d(2)
      if ((wld%bctr2d2(i)-wld%bctr2d2(i-1)).le.0.0) then
        write(ilog,*) 'Fatal. Bins have to be listed in ascending order (left-to-right) in histogram initialization input for &
 &Wang-Landau sampling and spacing has to be finite (dim. 2).'
        call fexit()
      end if
      if (abs(bd2*(wld%bctr2d2(i)-wld%bctr2d2(i-1))-1.0/dble(wld%nbins2d(2)-1)).gt.1.0e-2) then
        write(ilog,*) 'Fatal. Bin size in histogram initialization input for Wang-Landau sampling is not &
 &constant (even at low assumed precision).' 
        call fexit()
      else if ((abs(bd2*(wld%bctr2d2(i)-wld%bctr2d2(i-1))-1.0/dble(wld%nbins2d(2)-1)).gt.1.0e-9)&
 &              .AND.(resetbins(2).EQV..false.)) then
        write(ilog,*) 'Warning. Bin size in histogram initialization input for Wang-Landau sampling is not &
 &exactly constant. Resetting based on min/max range.' 
        resetbins(2) = .true.
      end if
    end do
    wld%g_binsz2d(2) = 1.0/(dble(wld%nbins2d(2)-1)*bd2)
    if (resetbins(2).EQV..true.) then
      do i=2,wld%nbins2d(2)
        wld%bctr2d2(i) = wld%bctr2d2(1) + (i-1)*wld%g_binsz2d(2)
      end do
    end if
  end if
!
! now read in actual data (along with bin centers for dim 1)
  if (wld%dimensionality.eq.2) allocate(vecr(wld%nbins2d(2)))
  do while(.true.)
    if (wld%dimensionality.eq.1) then
      read(wld%igin,*,iostat=iomess) wld%bctr(k+1),wld%g(k+1)
      if ((k+2).gt.wld%nbins) exit
    else if (wld%dimensionality.eq.2) then
      read(wld%igin,*,iostat=iomess) wld%bctr2d1(k+1),vecr(:)
      wld%g2d(k+1,1:wld%nbins2d(2)) = vecr(:)
      if ((k+2).gt.wld%nbins2d(1)) exit
    end if
    if (iomess.eq.IOSTAT_END) then
      if (wld%dimensionality.eq.2) then
        if ((k+1).lt.wld%nbins2d(1)) then
          write(*,*) 'Fatal. Initialization file for density histogram for Wang-Landau sampling (',fpath,') appears to be &
 &incomplete (unexpected end-of-file).'
          call fexit()
        end if
      end if
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing histogram initialization input for Wang-Landau sampling &
 &(got: ',fpath,').'
      call fexit()
    end if
    k = k + 1
  end do
  if (wld%dimensionality.eq.2) deallocate(vecr)
!
! sanity check bin ranges
  if (wld%dimensionality.eq.1) then
    bd1 = 1.0/(wld%bctr(wld%nbins) - wld%bctr(1))
    do i=2,wld%nbins
      if ((wld%bctr(i)-wld%bctr(i-1)).le.0.0) then
        write(ilog,*) 'Fatal. Bins have to be listed in ascending order (top-to-bottom) in histogram initialization input for &
 &Wang-Landau sampling and spacing has to be finite (dim. 1).'
        call fexit()
      end if
      if (abs(bd1*(wld%bctr(i)-wld%bctr(i-1))-1.0/dble(wld%nbins-1)).gt.1.0e-2) then
        write(ilog,*) 'Fatal. Bin size in histogram initialization input for Wang-Landau sampling is not &
 &constant (even at low assumed precision).' 
        call fexit()
      else if ((abs(bd1*(wld%bctr(i)-wld%bctr(i-1))-1.0/dble(wld%nbins-1)).gt.1.0e-9)&
 &              .AND.(resetbins(1).EQV..false.)) then
        write(ilog,*) 'Warning. Bin size in histogram initialization input for Wang-Landau sampling is not &
 &exactly constant. Resetting based on min/max range.'
        resetbins(1) = .true.
      end if
    end do
    wld%g_binsz = 1.0/(dble(wld%nbins-1)*bd1)
    if (resetbins(1).EQV..true.) then
      do i=2,wld%nbins
        wld%bctr(i) = wld%bctr(1) + (i-1)*wld%g_binsz
      end do
    end if
  else if (wld%dimensionality.eq.2) then
    bd1 = 1.0/(wld%bctr2d1(wld%nbins2d(1)) - wld%bctr2d1(1))
    do i=2,wld%nbins2d(1)
      if ((wld%bctr2d1(i)-wld%bctr2d1(i-1)).le.0.0) then
        write(ilog,*) 'Fatal. Bins have to be listed in ascending order (top-to-bottom) in histogram initialization input for &
 &Wang-Landau sampling and spacing has to be finite (dim. 1).'
        call fexit()
      end if
      if (abs(bd1*(wld%bctr2d1(i)-wld%bctr2d1(i-1))-1.0/dble(wld%nbins2d(1)-1)).gt.1.0e-2) then
        write(ilog,*) 'Fatal. Bin size in histogram initialization input for Wang-Landau sampling is not &
 &constant (even at low assumed precision).' 
        call fexit()
      else if ((abs(bd1*(wld%bctr2d1(i)-wld%bctr2d1(i-1))-1.0/dble(wld%nbins2d(1)-1)).gt.1.0e-9)&
 &              .AND.(resetbins(1).EQV..false.)) then
        write(ilog,*) 'Warning. Bin size in histogram initialization input for Wang-Landau sampling is not &
 &exactly constant. Resetting based on min/max range.'
        resetbins(1) = .true.
      end if
    end do
    wld%g_binsz2d(1) = 1.0/(dble(wld%nbins2d(1)-1)*bd1)
    if (resetbins(1).EQV..true.) then
      do i=2,wld%nbins2d(1)
        wld%bctr2d1(i) = wld%bctr2d1(1) + (i-1)*wld%g_binsz2d(1)
      end do
    end if
  end if
! 
  close(unit=wld%igin)
!
  if (wld%dimensionality.eq.1) then
    if (abs(wld%bctr(wld%nbins)-wld%g_max).gt.1.0e-9) then
!      write(ilog,*) 'Warning. Mismatch of last bin center while processing histogram initialization input for Wang-Landau&
! & sampling (input: file -> ',wld%bctr(wld%nbins),'; keyword -> ',wld%g_max,').'
      wld%g_max = wld%bctr(wld%nbins)
    end if
!
    wld%g_min = wld%g_max - (wld%nbins-1)*wld%g_binsz
    wld%minv = wld%g_min - 0.5*wld%g_binsz
    wld%minb = wld%nbins
  else if (wld%dimensionality.eq.2) then
    if (abs(wld%bctr2d1(wld%nbins2d(1))-wld%g_max2d(1)).gt.1.0e-9) then
!      write(ilog,*) 'Warning. Mismatch of last bin center while processing histogram initialization input for Wang-Landau&
! & sampling (dim. 1, input: file -> ',wld%bctr2d1(wld%nbins2d(1)),'; keyword -> ',wld%g_max2d(1),').'
      wld%g_max2d(1) = wld%bctr2d1(wld%nbins2d(1))
    end if
    if (abs(wld%bctr2d2(wld%nbins2d(2))-wld%g_max2d(2)).gt.1.0e-9) then
!      write(ilog,*) 'Warning. Mismatch of last bin center while processing histogram initialization input for Wang-Landau&
! & sampling (dim. 2, input: file -> ',wld%bctr2d2(wld%nbins2d(2)),'; keyword -> ',wld%g_max2d(2),').'
      wld%g_max2d(2) = wld%bctr2d2(wld%nbins2d(2))
    end if
!
    wld%g_min2d(:) = wld%g_max2d(:) - (wld%nbins2d(:)-1)*wld%g_binsz2d(:)
    wld%minv2d(:) = wld%g_min2d(:) - 0.5*wld%g_binsz2d(:)
    wld%minb2d(:) = wld%nbins2d(:)
  end if
!
  finit_wanglandau = .true.
!
end subroutine wl_parse_ginitfile
!
!--------------------------------------------------------------------------
!
