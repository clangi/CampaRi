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
! #######################################################
! ##                                                   ##
! ## subroutine ldmove -- Stochastic Dynamics moves    ##
! ##                                                   ##
! #######################################################
!     
! "ldmove" calls the LD integrator
!
!
subroutine cart_ldmove(istep,ndump,forcefirst)
!
  use iounit
  use math
  use atoms
  use molecule
  use forces
  use energies
  use movesets
  use system
  use units
  use mcsums
  use movesets
  use fyoc
  use torsn
  use cutoffs
  use shakeetal
!
  implicit none
!
  integer imol,j,istep,ndump,azero,aone,ttc,rs,i,incs,modstep
  RTYPE force3,fr1,fr2
  RTYPE ompl,ommi,bold,rndnew,normal
  RTYPE exgdt2,exgdt,dsplxyz(3),boxovol
  integer(KIND=8) t1,t2
  logical afalse,forcefirst
!
  afalse = .false.
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
! set the constants for the Langevin
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/(fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
!
! missing initialization
!
  azero = 0
  aone = 1
  nstep = nstep + 1
  mvcnt%nld = mvcnt%nld + 1
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!    write(*,*) (1.0 - exgdt)/(fric_ga*exgdt2),exgdt2,ompl
    do imol=1,nmol
!     sanity check in first step against (partially) unsupported 2-atom molecules
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
      end if
!     make sure rigid body coordinates are up-to-date
      call update_rigidm(imol)
    end do
!   get the deterministic force
    call System_Clock(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    if (no_shake.EQV..false.) call cart2cart_shake(azero)
    ens%insU = esave
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
!
    ttc = 0
!   set the frictional parameters
    do i=1,n
      if (mass(i).le.0.0) cycle
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        cart_ldp(i,j,1) = (u_dyn_kb*kelvin*ivms(i))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
        cart_ldp(i,j,2) = (u_dyn_kb*kelvin*ivms(i))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
        cart_ldp(i,j,3) = (u_dyn_kb*kelvin*ivms(i))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      end do
    end do
!
    call randomize_cart_velocities(azero)
    cart_a(:,:) = cart_v(:,:) ! velocity back-up
!
!   now do the first step integration
!   now loop over all molecules
    do imol=1,nmol
!
      call makeref_formol(imol)
!
      do i=atmol(imol,1),atmol(imol,2)
        if (mass(i).le.0.0) cycle
!       gather the displacements in atomic coordinate space, update velocities    
        do j=1,3
          if (cart_frz(i,j).EQV..true.) then
            cart_v(i,j) = 0.0
            dsplxyz(j) = 0.0
            cycle
          end if
!         "increment" frictional/stochastic parameter to starting value
          cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1))
          cart_ldp(i,j,5) = normal()
          rndnew = cart_ldp(i,j,4)*cart_ldp(i,j,5)
!         increment velocity from time zero to 1/2dt
          cart_v(i,j) = exgdt2*(cart_v(i,j) + ompl*u_dyn_fconv*dyn_dt*cart_f(i,j)*ivms(i) + rndnew)
!         increment position t to t + dt with staggered velocity at t + 1/2dt
          dsplxyz(j) = dyn_dt*cart_v(i,j)
        end do
        x(i) = x(i) + dsplxyz(1)
        y(i) = y(i) + dsplxyz(2)
        z(i) = z(i) + dsplxyz(3)
      end do
!
    end do
!
    if (no_shake.EQV..false.) call shake_wrap(azero)
!
    do imol=1,nmol
!
!     from fresh xyz, we need to recompute internals (rigid-body and torsions)
      call update_rigidm(imol)
      call update_comv(imol)
      call genzmat(imol)
!     and shift molecule into central cell if necessary
      call update_image(imol)
!     update grid association if necessary
      if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(rs)
        end do
      end if
!
    end do
!
!   update pointer arrays such that analysis routines can work properly 
    call zmatfyc2()
!
  else
!
    cart_a(:,:) = cart_v(:,:) ! velocity back-up
!
!   get the deterministic force
    call System_Clock(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    if (no_shake.EQV..false.) call cart2cart_shake(azero)
    ens%insU = esave
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
!
!   set the frictional parameters
    do i=1,n
      if (mass(i).le.0.0) cycle
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        cart_ldp(i,j,1) = (u_dyn_kb*kelvin*ivms(i))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
        cart_ldp(i,j,2) = (u_dyn_kb*kelvin*ivms(i))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
        cart_ldp(i,j,3) = (u_dyn_kb*kelvin*ivms(i))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      end do
    end do
!
!   now loop over all molecules
    do imol=1,nmol
      call makeref_formol(imol)
!
      do i=atmol(imol,1),atmol(imol,2)
        if (mass(i).le.0.0) cycle
!       gather the displacements in atomic coordinate space, update velocities    
        do j=1,3
          if (cart_frz(i,j).EQV..true.) then
            cart_v(i,j) = 0.0
            dsplxyz(j) = 0.0
            cycle
          end if
!         increment frictional/stochastic parameters to current values (trailing)
          bold = cart_ldp(i,j,2)/cart_ldp(i,j,4)
          cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1) + cart_ldp(i,j,3) - bold*bold)
          rndnew = bold*cart_ldp(i,j,5)
          cart_ldp(i,j,5) = normal()
          rndnew = rndnew + cart_ldp(i,j,4)*cart_ldp(i,j,5)
!         now increment velocity at time t - 1/2dt to t + 1/2dt
          cart_v(i,j) = exgdt2*(exgdt2*cart_v(i,j) + u_dyn_fconv*dyn_dt*cart_f(i,j)*ivms(i) + rndnew)
!         and finally increment position to t + dt
          dsplxyz(j) = dyn_dt*cart_v(i,j)
        end do
        x(i) = x(i) + dsplxyz(1)
        y(i) = y(i) + dsplxyz(2)
        z(i) = z(i) + dsplxyz(3)
      end do
!
    end do
!
    if (no_shake.EQV..false.) call shake_wrap(azero)
!
    incs = 0 
    ens%insR(6:7) = 0.
    modstep = 1
    if (istep.gt.nequil) modstep = mod(istep,nsancheck)
    do imol=1,nmol
!
!     from fresh xyz, we need to recompute internals (rigid-body and torsions)
      if ((ntormol(moltypid(imol)).eq.0).AND.(modstep.eq.0)) call compute_rotvel(imol,incs,ens%insR(6),ens%insR(7))
      call update_rigidm(imol)
      call update_comv(imol)
      call genzmat(imol)
!     and shift molecule into central cell if necessary
      call update_image(imol)
!     update grid association if necessary
      if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(rs)
        end do
      end if
!
    end do
    if (incs.gt.0) ens%avgcnt2 = ens%avgcnt2 + 1
    if (incs.gt.0) then
      ens%insR(5) = dble(incs)
    else
      ens%insR(5) = 0.0
    end if
!
!   update pointer arrays such that analysis routines can work properly 
    call zmatfyc2()
!
  end if
!
  cart_a(1:n,:) = (cart_v(1:n,:)-cart_a(1:n,:))/dyn_dt ! finite diff. accelerations
  ens%insR(1) = 0.0
  do i=1,n
    if (mass(i).gt.0.0) then
      ens%insR(1) = ens%insR(1) + sum((mass(i)*cart_a(i,:)/u_dyn_fconv - cart_f(i,:))**2)*ivms(i)
    end if
  end do
  call get_cart_ensv(boxovol,azero)
  call drift_removal(fr1,fr2,azero)
  call prt_curens(istep,boxovol,fr1,fr2)
!
end
!
!-------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine cart_ldmove_threads(istep,ndump,forcefirst,tpi)
!
  use iounit
  use math
  use atoms
  use sequen
  use molecule
  use forces
  use energies
  use movesets
  use system
  use units
  use mcsums
  use movesets
  use fyoc
  use torsn
  use cutoffs
  use shakeetal
  use threads
  use keys, ONLY: rndstat
!
  implicit none
!
  integer, INTENT(IN):: istep,ndump,tpi
  logical, INTENT(IN):: forcefirst
!
  integer imol,j,i,rs,incs,azero,ixx,modstep,atwo
  RTYPE dsplxyz(3),fr1,fr2,ompl,ommi,exgdt,exgdt2,rndnew,bold,temps(2)
  RTYPE boxovol,hlp1,ihlp1,hlp2,hlp3,hlp4
  logical ldummy,moljflags(4)
  integer(KIND=8) ttimer,t1,t2
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
! set the constants for the Langevin
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/(fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
  hlp1 = u_dyn_fconv*dyn_dt
  ihlp1 = 1.0/hlp1
  hlp2 = u_dyn_kb*kelvin*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
  hlp3 = u_dyn_kb*kelvin*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
  hlp4 = u_dyn_kb*kelvin*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
!
  azero = 0
  atwo = 2
  moljflags(:) = .true.
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
    moljflags(2:4) = .false.
    ixx = 1
    if (nmlgs.gt.0) then
      do while ((mlg_limits(ixx,5,tpi).lt.thr_limits(35,tpi)).AND.(ixx.lt.nmlgs))
        ixx = min(ixx+1,nmlgs)
      end do
    end if
    do imol=thr_limits(35,tpi),thr_limits(36,tpi)
      if (nmlgs.gt.0) then
        if (mlg_limits(ixx,5,tpi).eq.imol) then
          ixx = min(ixx+1,nmlgs)
          cycle
        end if
      end if
      call update_rigidm(imol)
    end do
    do i=1,nmlgs
!     covers update_rigidm, update_comv, genzmat, update_image
      call molops_threads(i,tpi,moljflags)
    end do
    moljflags(:) = .true.
!   note that the initialized velocities are exactly matched to the temperature
!   request -> no need to use thermostat in first step
!$OMP BARRIER
    call randomize_cart_velocities(tpi)
!$OMP BARRIER
  end if
!
!$OMP SINGLE
  nstep = nstep + 1
  mvcnt%nld = mvcnt%nld + 1
  ens%insR(1) = 0.0
  ens%insR(5:7) = 0.
!$OMP END SINGLE NOWAIT
! velocity backup
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),1) = cart_v(thr_limits(1,tpi):thr_limits(2,tpi),1)
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),2) = cart_v(thr_limits(1,tpi):thr_limits(2,tpi),2)
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),3) = cart_v(thr_limits(1,tpi):thr_limits(2,tpi),3)
! coordinate backup
  xref(thr_limits(1,tpi):thr_limits(2,tpi)) = x(thr_limits(1,tpi):thr_limits(2,tpi))
  yref(thr_limits(1,tpi):thr_limits(2,tpi)) = y(thr_limits(1,tpi):thr_limits(2,tpi))
  zref(thr_limits(1,tpi):thr_limits(2,tpi)) = z(thr_limits(1,tpi):thr_limits(2,tpi))
!
!$OMP BARRIER
  if (tpi.eq.1) call System_Clock(t1)
  ldummy = forcefirst
  call force3_threads(esterms,esterms_tr,esterms_lr,ldummy,esave)
!
  if (have_virtuals.EQV..true.) then
!$OMP BARRIER
    call cart2cart_shake(tpi)
  end if
!
  if (tpi.eq.1) then
    ens%insR(10) = ens%insU
    ens%insU = esave
    boxovol = ens%insV
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
  end if
!$OMP SINGLE
  call gen_rd_for_ld()
!$OMP END SINGLE NOWAIT 
!
  do i=thr_limits(1,tpi),thr_limits(2,tpi)
    if (mass(i).le.0.0) cycle
    do j=1,3
      if (cart_frz(i,j).EQV..true.) cycle
      cart_ldp(i,j,1) = hlp2*ivms(i)
      cart_ldp(i,j,2) = hlp3*ivms(i) 
      cart_ldp(i,j,3) = hlp4*ivms(i)
    end do
  end do
!
!$OMP BARRIER
!
! now loop over all atoms
  ixx = thr_limits(100,tpi)
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!
    do i=thr_limits(1,tpi),thr_limits(2,tpi)
!
      if (mass(i).le.0.0) cycle
!     gather the displacements in atomic coordinate space, update velocities
      do j=1,3
        if (cart_frz(i,j).EQV..true.) then
          cart_v(i,j) = 0.0
          dsplxyz(j) = 0.0
          cycle
        end if
        ixx = ixx + 1
!       "increment" frictional/stochastic parameter to starting value
        cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1))
        cart_ldp(i,j,5) = rndstat%rbuf(ixx) 
        rndnew = cart_ldp(i,j,4)*cart_ldp(i,j,5)
!       increment velocity from time zero to 1/2dt
        cart_v(i,j) = exgdt2*(cart_v(i,j) + ompl*hlp1*cart_f(i,j)*ivms(i) + rndnew)
!       increment position t to t + dt with staggered velocity at t + 1/2dt
        dsplxyz(j) = dyn_dt*cart_v(i,j)
      end do
      x(i) = x(i) + dsplxyz(1)
      y(i) = y(i) + dsplxyz(2)
      z(i) = z(i) + dsplxyz(3)
    end do
!
  else
!
    do i=thr_limits(1,tpi),thr_limits(2,tpi)
!
      if (mass(i).le.0.0) cycle
!     gather the displacements in atomic coordinate space, update velocities
      do j=1,3
        if (cart_frz(i,j).EQV..true.) then
          cart_v(i,j) = 0.0
          dsplxyz(j) = 0.0
          cycle
        end if
        ixx = ixx + 1
!       increment frictional/stochastic parameters to current values (trailing)
        bold = cart_ldp(i,j,2)/cart_ldp(i,j,4)
        cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1) + cart_ldp(i,j,3) - bold*bold)
        rndnew = bold*cart_ldp(i,j,5)
        cart_ldp(i,j,5) = rndstat%rbuf(ixx) ! normal()
        rndnew = rndnew + cart_ldp(i,j,4)*cart_ldp(i,j,5)
!       now increment velocity at time t - 1/2dt to t + 1/2dt
        cart_v(i,j) = exgdt2*(exgdt2*cart_v(i,j) + hlp1*cart_f(i,j)*ivms(i) + rndnew)
!       and finally increment position to t + dt
        dsplxyz(j) = dyn_dt*cart_v(i,j)
      end do
      x(i) = x(i) + dsplxyz(1)
      y(i) = y(i) + dsplxyz(2)
      z(i) = z(i) + dsplxyz(3)
    end do
!
  end if
!
!$OMP BARRIER
  if (no_shake.EQV..false.) call shake_wrap(tpi)
!
!  in NPT(E), update box
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!$OMP BARRIER
!$OMP SINGLE
    call manostat(boxovol)
!$OMP END SINGLE
  end if
!
!$OMP BARRIER
  if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
!   update grid association if needed 
    do rs=thr_limits(3,tpi),thr_limits(4,tpi)
      call updateresgp(rs)
    end do
  end if
!
! loop over all molecules again
  incs = 0 
  temps(1:2) = 0.0
  if (thr_dlb(10,1).gt.0) then
    if (tpi.eq.1) thr_dlb(10,2) = thr_dlb(10,2) + 1
    call System_Clock(count=ttimer)
    thr_timings(23,tpi) = thr_timings(23,tpi) + ttimer
  end if
  modstep = 1
  if (istep.gt.nequil) modstep = mod(istep,nsancheck)
  ixx = 1
  if (nmlgs.gt.0) then
    do while ((mlg_limits(ixx,5,tpi).lt.thr_limits(35,tpi)).AND.(ixx.lt.nmlgs))
      ixx = min(ixx+1,nmlgs)
    end do
  end if
  do imol=thr_limits(35,tpi),thr_limits(36,tpi)
    if ((ntormol(moltypid(imol)).eq.0).AND.(modstep.eq.0)) call compute_rotvel(imol,incs,temps(1),temps(2))
    if (nmlgs.gt.0) then
      if (mlg_limits(ixx,5,tpi).eq.imol) then
        ixx = min(ixx+1,nmlgs)
        cycle
      end if
    end if
!   we need to recompute internals
    call update_rigidm(imol)
    call update_comv(imol)
    call genzmat(imol)
!   and shift molecule into central cell if necessary
    call update_image(imol)
  end do
  if (thr_dlb(10,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(24,tpi) = thr_timings(24,tpi) + ttimer
  end if
!
  do i=1,nmlgs
!   covers update_rigidm, update_comv, genzmat, update_image
    call molops_threads(i,tpi,moljflags)
  end do
!
! update pointer arrays such that analysis routines can work properly
!$OMP BARRIER
  call zmatfyc_threads(tpi,atwo)
!
  bold = 0.0
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),1) = ihlp1*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),1) - &
 &                                                       cart_a(thr_limits(1,tpi):thr_limits(2,tpi),1))
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),2) = ihlp1*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),2) - &
 &                                                       cart_a(thr_limits(1,tpi):thr_limits(2,tpi),2))
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),3) = ihlp1*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),3) - &
 &                                                       cart_a(thr_limits(1,tpi):thr_limits(2,tpi),3))
  do i=thr_limits(1,tpi),thr_limits(2,tpi)
    if (mass(i).gt.0.0) then
      bold = bold + sum((mass(i)*cart_a(i,:) - cart_f(i,:))**2)*ivms(i)
    end if
  end do
  thr_rutil(5,tpi) = dble(incs)
  thr_rutil(1,tpi) = bold
  thr_rutil(6:7,tpi) = temps(1:2)
!$OMP BARRIER
!$OMP SINGLE
  ens%insR(1) = ens%insR(1) + sum(thr_rutil(1,1:thrdat%maxn))
  incs = nint(sum(thr_rutil(5,1:thrdat%maxn)))
  if (incs.gt.0) then
    ens%avgcnt2 = ens%avgcnt2 + 1
    ens%insR(5) = dble(incs)
    ens%insR(6:7) = sum(thr_rutil(6:7,1:thrdat%maxn),dim=2)
  end if
!$OMP END SINGLE
!
!  call System_Clock(thr_timings(66,tpi))
  call get_cart_ensv(boxovol,tpi)
  call drift_removal(fr1,fr2,tpi)
!  call System_Clock(thr_timings(67,tpi))
!  write(*,*) tpi,1.0e-3*(thr_timings(67,tpi)-thr_timings(66,tpi))
!$OMP BARRIER
!$OMP SINGLE
  call prt_curens(istep,boxovol,fr1,fr2)
!$OMP END SINGLE
!
end
!
#endif
!
!-----------------------------------------------------------------------
!
subroutine gen_rd_for_ld()
!
  use iounit
  use system
  use molecule
  use forces
  use atoms, ONLY: mass,n
  use keys, ONLY: rndstat
!
  implicit none
!
  RTYPE normal
  integer i,j,k
!
  if (allocated(rndstat%rbuf).EQV..false.) allocate(rndstat%rbuf(3*n))
!
  if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
    if (fycxyz.eq.1) then
      call fexit()
    else if (fycxyz.eq.2) then
      k = 0 
      do i=1,n
        if (mass(i).le.0.0) cycle
!       gather the displacements in atomic coordinate space, update velocities    
        do j=1,3
          if (cart_frz(i,j).EQV..true.) cycle
          k = k + 1
          rndstat%rbuf(k) = normal()
        end do
      end do
    end if
  else
    write(ilog,*) 'Encountered unsupported dynamics flag in gen_rd_for_ld(...). Offending code is ',dyn_mode,'. &
 &Please report this problem!'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------
!
! this routine is only used to ensure that a run which was restarted from
! an MC restart file. this is necessary because the non-first step LD
! algorithm requires ldp(i,4) and ldp(i,5)
! of course, it required inertial masses to be set and up-to-date
!
subroutine init_cart_ldps()
!
  use iounit
  use forces
  use atoms
  use system
  use units
!
  implicit none
!
  integer i,j
  RTYPE ompl,ommi,normal
  RTYPE exgdt2,exgdt
!
! set the constants
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/&
 &    (fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
  exgdt = exp(-fric_ga*dyn_dt)
  do i=1,n
    if (mass(i).le.0.0) cycle
    do j=1,3
      cart_ldp(i,j,1) = (u_dyn_kb*kelvin*ivms(i))*&
 &                  (2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
      cart_ldp(i,j,2) = (u_dyn_kb*kelvin*ivms(i))*&
 &                  (2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
      cart_ldp(i,j,3) = (u_dyn_kb*kelvin*ivms(i))*&
 &                  (2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      cart_ldp(i,j,4) = sqrt(cart_ldp(i,j,1))
      cart_ldp(i,j,5) = normal()
    end do
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a routine that does nothing but assemble CLD parameters into provided arrays
!
subroutine manage_clds(vars,mode,tpi)
!
  use iounit
  use forces
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: mode,tpi
!
  RTYPE vars(n,9)
#ifdef ENABLE_THREADS
  integer sta,sto
!
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
#endif
!
  if (mode.eq.1) then
!   store in vars
#ifdef ENABLE_THREADS
    vars(sta:sto,1) = cart_v(sta:sto,1)
    vars(sta:sto,2) = cart_v(sta:sto,2)
    vars(sta:sto,3) = cart_v(sta:sto,3)
    vars(sta:sto,4) = cart_ldp(sta:sto,1,4)
    vars(sta:sto,5) = cart_ldp(sta:sto,1,5)
    vars(sta:sto,6) = cart_ldp(sta:sto,2,4)
    vars(sta:sto,7) = cart_ldp(sta:sto,2,5)
    vars(sta:sto,8) = cart_ldp(sta:sto,3,4)
    vars(sta:sto,9) = cart_ldp(sta:sto,3,5)
#else
    vars(1:n,1) = cart_v(1:n,1)
    vars(1:n,2) = cart_v(1:n,2)
    vars(1:n,3) = cart_v(1:n,3)
    vars(1:n,4) = cart_ldp(1:n,1,4)
    vars(1:n,5) = cart_ldp(1:n,1,5)
    vars(1:n,6) = cart_ldp(1:n,2,4)
    vars(1:n,7) = cart_ldp(1:n,2,5)
    vars(1:n,8) = cart_ldp(1:n,3,4)
    vars(1:n,9) = cart_ldp(1:n,3,5)
#endif
  else if (mode.eq.2) then
#ifdef ENABLE_THREADS
    cart_v(sta:sto,1) = vars(sta:sto,1)
    cart_v(sta:sto,2) = vars(sta:sto,2)
    cart_v(sta:sto,3) = vars(sta:sto,3)
    cart_ldp(sta:sto,1,4) = vars(sta:sto,4)
    cart_ldp(sta:sto,1,5) = vars(sta:sto,5)
    cart_ldp(sta:sto,2,4) = vars(sta:sto,6)
    cart_ldp(sta:sto,2,5) = vars(sta:sto,7)
    cart_ldp(sta:sto,3,4) = vars(sta:sto,8)
    cart_ldp(sta:sto,3,5) = vars(sta:sto,9)
#else
    cart_v(1:n,1) = vars(1:n,1)
    cart_v(1:n,2) = vars(1:n,2)
    cart_v(1:n,3) = vars(1:n,3)
    cart_ldp(1:n,1,4) = vars(1:n,4)
    cart_ldp(1:n,1,5) = vars(1:n,5)
    cart_ldp(1:n,2,4) = vars(1:n,6)
    cart_ldp(1:n,2,5) = vars(1:n,7)
    cart_ldp(1:n,3,4) = vars(1:n,8)
    cart_ldp(1:n,3,5) = vars(1:n,9)
#endif
  else
    write(ilog,*) 'Fatal. Called manage_clds(...) with unknown mode (identifier is ',&
  &mode,'). This is a bug.'
    call fexit()
  end if
!
  end
!
!-----------------------------------------------------------------------
!
