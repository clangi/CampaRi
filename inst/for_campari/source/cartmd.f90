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
subroutine cart_mdmove(istep,ndump,forcefirst)
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
!
  implicit none
!
  integer, INTENT(IN):: istep,ndump
  logical, INTENT(IN):: forcefirst
!
  integer imol,j,i,rs,incs,azero,modstep
  RTYPE force3,dsplxyz(3),hlp1,fr1,fr2
  RTYPE tscs(tstat%n_tgrps),tsc,boxovol
  integer(KIND=8) t1,t2
!
  if ((tstat%flag.eq.3).AND.(pstat%flag.eq.3)) then
    call cart_mdmove_npt(istep,ndump,forcefirst)
    return
  end if
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
  nstep = nstep + 1
  mvcnt%nmd = mvcnt%nmd + 1
  azero = 0 
!
!
! initialize velocities
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!
    hlp1 = 0.5*u_dyn_fconv*dyn_dt
    do imol=1,nmol
!     sanity check in first step against (partially) unsupported 2-atom molecules
!      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
!        write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
! &ed in dynamics. Check back later.'
!        call fexit()
!      end if
      call update_rigidm(imol)
    end do
    call System_Clock(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    if (have_virtuals.EQV..true.) call cart2cart_shake(azero)
    ens%insR(10) = ens%insU
    ens%insU = esave
    boxovol = ens%insV
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
!
!   note that the initialized velocities are exactly matched to the temperature
!   request -> no need to use thermostat in first step
    call randomize_cart_velocities(azero)
    cart_a(:,:) = cart_v(:,:) ! velocity back-up
!
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
!         increment velocity at time t - 1/2dt to t + 1/2dt
          cart_v(i,j) = cart_v(i,j) + hlp1*cart_f(i,j)/mass(i)
!         increment position t to t + dt with staggered velocity at t + 1/2dt
          dsplxyz(j) = dyn_dt*cart_v(i,j)
        end do
        x(i) = x(i) + dsplxyz(1)
        y(i) = y(i) + dsplxyz(2)
        z(i) = z(i) + dsplxyz(3)
!
      end do
!
    end do
!
    if (no_shake.EQV..false.) call shake_wrap(azero)
!
!   in NPT(E), update box
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      call manostat(boxovol)
    end if
!
!   now loop over all molecules again
    do imol=1,nmol
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
    end do
!
!   update pointer arrays such that analysis routines can work properly 
    call zmatfyc2()
!
  else
!
    hlp1 = u_dyn_fconv*dyn_dt
    cart_a(:,:) = cart_v(:,:) ! velocity back-up
    call System_Clock(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    if (have_virtuals.EQV..true.) call cart2cart_shake(azero)
    ens%insR(10) = ens%insU
    ens%insU = esave
    boxovol = ens%insV
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
!
    if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
!     get the temperature re-scaling factor from thermostat
      call thermostat(tscs)
    else
      tscs(:) = 1.0
    end if
!
!   now loop over all molecules
    do imol=1,nmol
!
!     T-group specific T-rescaling
      tsc = tscs(tstat%molgrp(imol))
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
!         increment velocity at time t - 1/2dt to t + 1/2dt
          cart_v(i,j) = tsc*cart_v(i,j) + hlp1*cart_f(i,j)/mass(i)
!         increment position t to t + dt with staggered velocity at t + 1/2dt
          dsplxyz(j) = dyn_dt*cart_v(i,j)
        end do
        x(i) = x(i) + dsplxyz(1)
        y(i) = y(i) + dsplxyz(2)
        z(i) = z(i) + dsplxyz(3)
!
      end do
!
    end do
!
    if (no_shake.EQV..false.) call shake_wrap(azero)
!
!   in NPT(E), update box
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      call manostat(boxovol)
    end if
!
!   loop over all molecules again
    incs = 0 
    ens%insR(6:7) = 0.
    modstep = 1
    if (istep.gt.nequil) modstep = mod(istep,nsancheck)
    do imol=1,nmol 
!     we need to recompute internals
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
  ens%insR(1) = 0.0 ! GLCF
  do i=1,n
    if (mass(i).gt.0.0) then
      ens%insR(1) = ens%insR(1) + sum((mass(i)*cart_a(i,:)/u_dyn_fconv - cart_f(i,:))**2)/mass(i)
    end if
  end do
  call get_cart_ensv(boxovol,azero)
  call drift_removal(fr1,fr2,azero)
  call prt_curens(istep,boxovol,fr1,fr2) ! print
!
end
!
!-------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine cart_mdmove_threads(istep,ndump,forcefirst,tpi)
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
!
  implicit none
!
  integer, INTENT(IN):: istep,ndump,tpi
  logical, INTENT(IN):: forcefirst
!
  integer imol,j,i,rs,incs,azero,modstep,atwo
  RTYPE dsplxyz(3),fr1,fr2
  RTYPE tscs(tstat%n_tgrps),tsc,boxovol,hlp1,ihlp1,temps(2)
  logical ldummy,moljflags(4)
  integer(KIND=8) ttimer,t1,t2
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
  azero = 0
  atwo = 2
  moljflags(:) = .true.
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
    moljflags(2:4) = .false.
    hlp1 = 0.5*u_dyn_fconv*dyn_dt
    do imol=thr_limits(35,tpi),thr_limits(36,tpi)
      if (thr_mlgix(imol).gt.0) cycle
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
  else
    hlp1 = u_dyn_fconv*dyn_dt
  end if
  ihlp1 = 1.0/(u_dyn_fconv*dyn_dt)! hlp1
!
!$OMP SINGLE
  nstep = nstep + 1
  mvcnt%nmd = mvcnt%nmd + 1
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
  if ((istep.gt.1).AND.(forcefirst.EQV..false.)) then
!$OMP SINGLE 
    if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
!     get the temperature re-scaling factor from thermostat
      call thermostat(tscs)
    else
      tscs(:) = 1.0
    end if
!$OMP END SINGLE COPYPRIVATE(tscs)
  else
    tscs(:) = 1.0
  end if
!
!$OMP BARRIER
!
! now loop over all molecules
  do i=thr_limits(1,tpi),thr_limits(2,tpi)
!
!   T-group specific T-rescaling
    tsc = tscs(tstat%molgrp(molofrs(atmres(i)))) ! tscs(tstat%molgrp(imol))
    if (mass(i).le.0.0) cycle
!   gather the displacements in atomic coordinate space, update velocities
    do j=1,3
      if (cart_frz(i,j).EQV..true.) then
        cart_v(i,j) = 0.0
        dsplxyz(j) = 0.0
        cycle
      end if
!     increment velocity at time t - 1/2dt to t + 1/2dt
      cart_v(i,j) = tsc*cart_v(i,j) + hlp1*cart_f(i,j)/mass(i)
!     increment position t to t + dt with staggered velocity at t + 1/2dt
      dsplxyz(j) = dyn_dt*cart_v(i,j)
    end do
    x(i) = x(i) + dsplxyz(1)
    y(i) = y(i) + dsplxyz(2)
    z(i) = z(i) + dsplxyz(3)
!
  end do
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
!   upgrade grid association if needed 
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
  do imol=thr_limits(35,tpi),thr_limits(36,tpi)
    if ((ntormol(moltypid(imol)).eq.0).AND.(modstep.eq.0)) call compute_rotvel(imol,incs,temps(1),temps(2))
    if (thr_mlgix(imol).gt.0) cycle
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
  tsc = 0.0
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),1) = ihlp1*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),1) - &
 &                                                       cart_a(thr_limits(1,tpi):thr_limits(2,tpi),1))
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),2) = ihlp1*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),2) - &
 &                                                       cart_a(thr_limits(1,tpi):thr_limits(2,tpi),2))
  cart_a(thr_limits(1,tpi):thr_limits(2,tpi),3) = ihlp1*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),3) - &
 &                                                       cart_a(thr_limits(1,tpi):thr_limits(2,tpi),3))
  do i=thr_limits(1,tpi),thr_limits(2,tpi)
    if (mass(i).gt.0.0) then
      tsc = tsc + sum((mass(i)*cart_a(i,:) - cart_f(i,:))**2)/mass(i)
    end if
  end do
  thr_rutil(5,tpi) = dble(incs)
  thr_rutil(1,tpi) = tsc
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
subroutine cart_mdmove_npt(istep,ndump,forcefirst)
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
!
  implicit none
!
  integer imol,j,istep,ndump,i
  integer boxdof,alldof,alldof2,azero
  RTYPE force3,fr1,fr2,expt2,boxovol
  integer(KIND=8) t1,t2
  RTYPE mfst2,tfst2,mtst2,expt3,expt4,expt1
  logical afalse,forcefirst
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
  azero = 0
  afalse = .false.
  nstep = nstep + 1
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      boxdof = 1
    else if (pstat%flag.eq.3) then
      alldof = 3*(n)-n_constraints
      alldof2 = 3*(n+1) - n_constraints
      boxdof = 0
    end if
  else if (ens%flag.eq.1) then
    if (tstat%flag.eq.3) then
      alldof = 3*(n)-n_constraints
      alldof2 = 3*(n+1) - n_constraints
      boxdof = 0
    end if
  else if (ens%flag.eq.2) then
    alldof = 3*(n)-n_constraints
    alldof2 = alldof
    boxdof = 0
  end if
  if ((tstat%flag.eq.3).AND.(pstat%flag.eq.3)) then
    mfst2 = 1.0/(pstat%params(1)*pstat%params(1))
    tfst2 = 1.0/(tstat%params(1)*tstat%params(1))
    mtst2 = mfst2/tfst2
  end if
!
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!
    if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
      tstat%params(4) = 0.0 ! Jacobian for T-stat particle
      tstat%params(2) = 0.0 ! initial velocity for T-stat particle
    else
      tstat%params(4) = 0.0 ! Jacobian for T-stat particle
      tstat%params(2) = 0.0
    end if
    if ((ens%flag.eq.4).OR.(ens%flag.eq.3)) then
      pstat%params(2) = 0.0 ! initial velocity for P-stat particle
    else
      pstat%params(2) = 0.0
    end if
!
!   note that the initialized velocities are exactly matched to the temperature
!   request -> no need to use thermostat in first step
    do imol=1,nmol
!     sanity check in first step against (partially) unsupported 2-atom molecules
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
      end if
      call update_rigidm(imol)
    end do
!   initialize velocities
    call randomize_cart_velocities(azero)
    ens%insT = kelvin
!
  end if
!
! gather deterministic force
  call System_Clock(t1)
  esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
  ens%insR(10) = ens%insU
  ens%insU = esave
  call System_Clock(t2)
  time_energy = time_energy + 1.0*(t2 - t1)
  boxovol = ens%insV
! get unitless force on manostat particle
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (bnd_type.eq.1) then
      if (bnd_shape.eq.1) then
        bnd_fr = (invtemp/alldof)*&
 &   (1.0*(ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) - &
 &    3.0*(ens%insV*extpress/u_dyn_virconv))
      end if
    else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
      if (bnd_shape.eq.2) then
        bnd_fr = (invtemp/alldof)*&
 &   (1.0*(ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) - &
 &    3.0*(ens%insV*extpress/u_dyn_virconv) - bnd_params(4)*bnd_fr)
      end if
    end if
  end if
!
  if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
!   accumulate thermostat Jacobian 
    tstat%params(4) = tstat%params(4) + tstat%params(2)*dyn_dt*0.5
!   get velocity re-scaling factors (might use a series expansion to avoid singularity)
    if (abs(tstat%params(2)*dyn_dt*0.5).le.1.0e-6) then
      expt1 = -tstat%params(2)*dyn_dt*0.5
      expt1 = 1.0 + expt1 + expt1*expt1/2.0 + expt1*expt1*expt1/6.0
      expt2 = -dyn_dt*0.5
      expt2 = expt2 + expt2*expt2*tstat%params(2)/2.0 + &
 &  expt2*expt2*expt2*tstat%params(2)*tstat%params(2)/6.0 + &
 & expt2*expt2*expt2*expt2*&
 &  tstat%params(2)*tstat%params(2)*tstat%params(2)/12.0
      expt2 = -u_dyn_fconv*expt2
    else
      expt1 = exp(-tstat%params(2)*dyn_dt*0.5)
      expt2 = u_dyn_fconv*&
 &       (1.0 - exp(-tstat%params(2)*dyn_dt*0.5))/tstat%params(2)
    end if
  else
    expt1 = 1.0
    expt2 = 0.5*u_dyn_fconv*dyn_dt
  end if
!
! now loop over all molecules
  do imol=1,nmol
!
    call makeref_formol(imol)
!
    do i=atmol(imol,1),atmol(imol,2)
      if (mass(i).le.0.0) cycle
!     gather the displacements in atomic coordinate space, update velocities
      do j=1,3
        cart_v(i,j) = expt1*cart_v(i,j)+expt2*cart_f(i,j)/mass(i)
      end do
    end do
  end do
! 
! increment manostat velocity   
  if (ens%flag.eq.3) then
    pstat%params(2) = pstat%params(2) + 0.5*dyn_dt*&
 &   (mfst2*bnd_fr + alldof2*mtst2*tstat%params(2)*tstat%params(2))
  else if (ens%flag.eq.4) then
    pstat%params(2) = pstat%params(2) + 0.5*dyn_dt*mfst2*(bnd_fr&
 &                         + alldof2)
  end if
! calculate unitless force on thermostat particle
  if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
    ens%insK = 0.0
    do j=1,3
      ens%insK = ens%insK + (1.0/u_dyn_fconv)* &
 &           0.5*sum(cart_v(1:n,j)*cart_v(1:n,j)*mass(1:n))
    end do
    tstat%params(3) = -1.0 + 2.0*ens%insK*invtemp/alldof
!   increment thermostat velocity
    tstat%params(2) = tstat%params(2) + 0.5*dyn_dt*&
 &                                   tfst2*tstat%params(3)
    tstat%params(5) = 0.5*tstat%params(2)*tstat%params(2)/tfst2
    tstat%params(6) = tstat%params(4)
  end if
!
! calculate kinetic energy of manostat and thermostat, store Jacobian
  if ((ens%flag.eq.4).OR.(ens%flag.eq.3)) then
    pstat%params(5) = 0.5*pstat%params(2)*pstat%params(2)/mfst2
  end if
!
! calculate ensemble variables at proper place during integration
  call get_cart_ensv(boxovol,azero)
!
! now do second half-step
!
  if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
!   increment thermostat velocity again
    tstat%params(2) = tstat%params(2) + 0.5*dyn_dt*&
 &                                   tfst2*tstat%params(3)
!   accumulate jacobian again
    tstat%params(4) = tstat%params(4) + tstat%params(2)*dyn_dt*0.5
!   re-set variables for velocity increment
    if (abs(tstat%params(2)*dyn_dt*0.5).le.1.0e-6) then
      expt1 = -tstat%params(2)*dyn_dt*0.5
      expt1 = 1.0 + expt1 + expt1*expt1/2.0 + expt1*expt1*expt1/6.0
      expt2 = -dyn_dt*0.5
      expt2 = expt2 + expt2*expt2*tstat%params(2)/2.0 + &
 &  expt2*expt2*expt2*tstat%params(2)*tstat%params(2)/6.0 + &
 & expt2*expt2*expt2*expt2*&
 &  tstat%params(2)*tstat%params(2)*tstat%params(2)/12.0
      expt2 = -u_dyn_fconv*expt2
    else
      expt1 = exp(-tstat%params(2)*dyn_dt*0.5)
      expt2 = u_dyn_fconv*&
 &       (1.0 - exp(-tstat%params(2)*dyn_dt*0.5))/tstat%params(2)
    end if
  else
    expt1 = 1.0
    expt2 = 0.5*u_dyn_fconv*dyn_dt
  end if
!
! loop over all molecules to finalize velocities
  do imol=1,nmol
!   T-group specific T-rescaling
    do i=atmol(imol,1),atmol(imol,2)
      if (mass(i).le.0.0) cycle
!     gather the displacements in atomic coordinate space, update velocities
      do j=1,3
!       increment velocity at time t - 1/2dt to t + 1/2dt
        cart_v(i,j) = expt1*cart_v(i,j)+expt2*cart_f(i,j)/mass(i)
      end do
    end do
  end do
! increment manostat velocity again
  if (ens%flag.eq.3) then
    pstat%params(2) = pstat%params(2) + 0.5*dyn_dt*&
 &   (mfst2*bnd_fr + alldof2*mtst2*tstat%params(2)*tstat%params(2))
  else if (ens%flag.eq.4) then
    pstat%params(2) = pstat%params(2) + 0.5*dyn_dt*mfst2*(bnd_fr&
 &                         + alldof2)
  end if
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!   increment volume (isotropic)
    if (bnd_shape.eq.1) then
      bnd_params(1:3)=exp(pstat%params(2)*dyn_dt)*bnd_params(1:3)
    else if (bnd_shape.eq.2) then
      bnd_params(4) = exp(pstat%params(2)*dyn_dt)*bnd_params(4)
    end if
    call update_bound(afalse)
  end if
! increment thermostat velocity for the final time
  if (ens%flag.eq.3) then
    tstat%params(2) = exp(-alldof2*pstat%params(2)*dyn_dt)*&
 &                          tstat%params(2)
  end if
!
! and eventually atomic positions
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (abs(pstat%params(2)*dyn_dt).le.1.0e-6) then
      expt3 = pstat%params(2)*dyn_dt
      expt3 = 1.0 + expt3 + expt3*expt3/2.0+expt3*expt3*expt3/6.0
      expt4 = dyn_dt
      expt4 = expt4 + expt4*expt4*pstat%params(2)/2.0 + &
 &  expt4*expt4*expt4*pstat%params(2)*pstat%params(2)/6.0 + &
 & expt4*expt4*expt4*expt4*&
 &  pstat%params(2)*pstat%params(2)*pstat%params(2)/12.0
    else
      expt3 = exp(pstat%params(2)*dyn_dt)
      expt4 = (expt3 - 1.0)/pstat%params(2)
    end if
  else
    expt3 = 1.0
    expt4 = dyn_dt
  end if
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      x(i) = expt3*x(i) + expt4*cart_v(i,1)
      y(i) = expt3*y(i) + expt4*cart_v(i,2)
      z(i) = expt3*z(i) + expt4*cart_v(i,3)
    end do
!   we need to recompute internals
    call update_rigidm(imol)
    call update_comv(imol)
    call genzmat(imol)
!   and shift molecule into central cell if necessary
    call update_image(imol)
  end do
!  call gen_cartv(mshfs)
!
! update pointer arrays such that analysis routines can work properly 
  call zmatfyc2()
!
! remove c.o.m drift (in the right place?)
  call drift_removal(fr1,fr2,azero)
!
 997  format(' Kinetic E   | Potential E | Total E     | Te&
 &mperature | Drift-xyz   | Drift-Euler |')
 998  format(' Kinetic E   | Potential E | Enthalpy    | Te&
 &mperature | Pressure    | Box volume  |')
 999  format('---------------------------------------------&
 &---------------------------------------')
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!   write header line
    if (nsim.ge.nsancheck) then
      if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
        write(ilog,997)
        write(ilog,999) 
      else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
        write(ilog,998)
        write(ilog,999)
      end if
    end if
  end if
  if (mod(nstep,nsancheck).eq.0) then
    if (ens%flag.eq.2) then
      write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU,ens%insT,&
 &fr1,fr2
    else if (ens%flag.eq.1) then
      write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU + &
 &   (alldof/invtemp)*(tstat%params(6) + tstat%params(5)),ens%insT,&
 &fr1,fr2
    else if (ens%flag.eq.3) then
      write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU + &
 &extpress*ens%insV/u_dyn_virconv+(alldof/invtemp)*&
 &     (tstat%params(6) + tstat%params(5) + pstat%params(5)),&
 &ens%insT,ens%insP,boxovol/1.0e6
    else if (ens%flag.eq.4) then
      write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU + &
 &extpress*ens%insV/u_dyn_virconv + &
 & (alldof/invtemp)*pstat%params(5),ens%insT,ens%insP,boxovol/1.0e6
     end if
  end if
!
   56 format(50(g13.6,1x))
!
end
!
!-----------------------------------------------------------------------
!
subroutine update_comv(imol)
!
  use molecule
  use forces
  use atoms
!
  implicit none
!
  integer i,j,imol
!
! a trivial subroutine to populate dc_di(imol)%v(1..3) with center of mass
! velocity based on Cartesian velocities cart_v
!
  do j=1,3
    dc_di(imol)%v(j) = 0.0
  end do
  do i=atmol(imol,1),atmol(imol,2)
    do j=1,3
      dc_di(imol)%v(j) = dc_di(imol)%v(j) + mass(i)*cart_v(i,j)
    end do
  end do
  do j=1,3
    if (movingmass(imol,j).gt.0.0) dc_di(imol)%v(j) = dc_di(imol)%v(j)/movingmass(imol,j)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this function can be used to infer rigid-body rotational velocities for molecules rigidified with constraints
!
subroutine compute_rotvel(imol,incs,sumT1,sumT2)
!
  use molecule
  use atoms
  use units
  use system
!
  implicit none
!
  integer i,j,imol,incs
  RTYPE cds(3*(atmol(imol,2)-atmol(imol,1)+1),2),tvec(3),qrot(4),co(3),cn(3),rgtenm(3,3,2),hlp1(2),hlp2(2),hlp3(2),hlp4(2)
  RTYPE hlp5(2),hlp6(2),sumT1,sumT2,vvec(3),rotvec(3)
!
  if ((atmol(imol,2)-atmol(imol,1)).eq.1) return ! does not work for linear molecules
!
! a trivial subroutine to populate dc_di(imol)%v(1..3) with center of mass
! velocity based on Cartesian velocities cart_v
!
  do i=atmol(imol,1),atmol(imol,2)
    j = i-atmol(imol,1)+1
    cds(3*j-2,1) = xref(i)
    cds(3*j-1,1) = yref(i)
    cds(3*j,1) = zref(i)
    cds(3*j-2,2) = x(i)
    cds(3*j-1,2) = y(i)
    cds(3*j,2) = z(i)
  end do
!
  j = atmol(imol,2)-atmol(imol,1)+1
  call align_3D_wt(j,cds(:,1),cds(:,2),mass(atmol(imol,1):atmol(imol,2)),tvec,qrot,co,cn)
!
! from the optimal transformation, infer RB velocities
  vvec(1:3) = tvec(:)/dyn_dt
  rotvec(1:3) = 2.0*asin(qrot(2:4))/dyn_dt
!
! get inertia tensors
  rgtenm(:,:,:) = 0.0
  do i=atmol(imol,1),atmol(imol,2)
    hlp1(1) = xref(i)-co(1)
    hlp2(1) = yref(i)-co(2)
    hlp3(1) = zref(i)-co(3)
    hlp1(2) = x(i)-cn(1)
    hlp2(2) = y(i)-cn(2)
    hlp3(2) = z(i)-cn(3)
    hlp4(1:2) = hlp1(1:2)*hlp1(1:2)
    hlp5(1:2) = hlp2(1:2)*hlp2(1:2)
    hlp6(1:2) = hlp3(1:2)*hlp3(1:2)
    rgtenm(1,1,1:2) = rgtenm(1,1,1:2) + mass(i)*(hlp5(1:2)+hlp6(1:2))
    rgtenm(1,2,1:2) = rgtenm(1,2,1:2) - mass(i)*hlp1(1:2)*hlp2(1:2)
    rgtenm(1,3,1:2) = rgtenm(1,3,1:2) - mass(i)*hlp1(1:2)*hlp3(1:2)
    rgtenm(2,2,1:2) = rgtenm(2,2,1:2) + mass(i)*(hlp4(1:2)+hlp6(1:2))
    rgtenm(2,3,1:2) = rgtenm(2,3,1:2) - mass(i)*hlp2(1:2)*hlp3(1:2)
    rgtenm(3,3,1:2) = rgtenm(3,3,1:2) + mass(i)*(hlp4(1:2)+hlp5(1:2))
  end do
  rgtenm(2,1,1:2) = rgtenm(1,2,1:2)
  rgtenm(3,1,1:2) = rgtenm(1,3,1:2)
  rgtenm(3,2,1:2) = rgtenm(2,3,1:2)
!
  sumT2 = sumT2 + dot_product(rotvec(1:3),matmul((rgtenm(:,:,1)+rgtenm(:,:,2))*0.5,rotvec(1:3)))/u_dyn_kb
  sumT1 = sumT1 + molmass(moltypid(imol))*dot_product(vvec(1:3),vvec(1:3))/u_dyn_kb
  incs = incs + 3
!
end
!
!-----------------------------------------------------------------------
!
! a routine that does nothing but assemble CMD parameters into provided arrays
!
subroutine manage_cmds(vars,mode,tpi)
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
  RTYPE vars(n,3)
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
#else
    vars(1:n,1) = cart_v(1:n,1)
    vars(1:n,2) = cart_v(1:n,2)
    vars(1:n,3) = cart_v(1:n,3)
#endif
  else if (mode.eq.2) then
#ifdef ENABLE_THREADS
    cart_v(sta:sto,1) = vars(sta:sto,1)
    cart_v(sta:sto,2) = vars(sta:sto,2)
    cart_v(sta:sto,3) = vars(sta:sto,3)
#else
    cart_v(1:n,1) = vars(1:n,1)
    cart_v(1:n,2) = vars(1:n,2)
    cart_v(1:n,3) = vars(1:n,3)
#endif
  else
    write(ilog,*) 'Fatal. Called manage_cmds(...) with unknown mode (identifier is ',&
  &mode,'). This is a bug.'
    call fexit()
  end if
!
  end
!
!-----------------------------------------------------------------------
!
