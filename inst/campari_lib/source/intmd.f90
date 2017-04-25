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
! MAIN AUTHOR:   Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! #######################################################
! ##                                                   ##
! ## this subroutine propagates internal degrees of    ##
! ## freedom with a leapfrog integrator assuming       ##
! ## Newtonain (fully ballistic) dynamics              ##
! ##                                                   ##
! #######################################################
!     
!
subroutine int_mdmove(istep,ndump,forcefirst)
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
  use zmatrix
!
  implicit none
!
  integer imol,j,istep,ndump,aone,azero,ttc,rs,i
  RTYPE force3
  RTYPE t1,t2,t0,ilow,ihigh,frac1,frac2
  RTYPE netf,tscs(tstat%n_tgrps),boxovol,tsc
  RTYPE, ALLOCATABLE:: mshfs(:,:)
  logical afalse,forcefirst
  integer(KIND=8) tt1,tt2
!
  allocate(mshfs(3,nmol))
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
  afalse = .false.
  azero = 0
  aone = 1
  nstep = nstep + 1
  mvcnt%nmd = mvcnt%nmd + 1
  cur_trans(:) = 0.0
!
! initialize velocities
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!
    do imol=1,nmol
!     sanity check in first step against (partially) unsupported 2-atom molecules
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
      end if
      call update_rigidm(imol)
!     zero out internal coordinate forces
      dc_di(imol)%f(:) = 0.0
      dc_di(imol)%olddat(:,1) = dc_di(imol)%im(:)
      dc_di(imol)%olddat(:,2) = dc_di(imol)%im(:) ! safe side
      dc_di(imol)%v(:) = 0.0
    end do
!
    call System_Clock(tt1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    ens%insR(10) = ens%insU ! backup
    ens%insU = esave
    boxovol = ens%insV
    call System_Clock(tt2)
    call cart2int_f(skip_frz,azero)
    if (fudge_mass.EQV..true.) call fudge_masses()
    time_energy = time_energy + 1.0*(tt2 - tt1)
!
!   note that the initialized velocities are exactly matched to the temperature
!   request -> no need to use thermostat in first step
!   prior call to cart2int_f or similar is necessary to populate inertia initially
    call randomize_velocities(azero)
!
!   the velocities just need to be incremented by a half-step without having access to any
!   past information, coordinates are advanced a full step
    do imol=1,nmol
!
!     first: center-of-mass translation (linear motion)
      do j=1,3
!       increment velocity at time t_0 to t_0 + 1/2dt
        if (dc_di(imol)%frz(j).EQV..true.) then
          cur_trans(j) = 0.0
          cycle
        end if
        dc_di(imol)%v(j) = dc_di(imol)%v(j) + u_dyn_fconv*0.5*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j)
        cur_trans(j) = dyn_dt*dc_di(imol)%v(j)
      end do
!
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%v(ttc) = 0.0
          else
            t0 = 0.5*dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
            dc_di(imol)%v(ttc) = dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
          end if
          dc_di(imol)%incr(ttc) = dyn_dt*dc_di(imol)%v(ttc) ! in degrees
        end do
!       now transfer
        cur_rot(1:3) = dc_di(imol)%incr(4:6)/RADIAN
!       note that this will edit the Z-matrix and populate ztorpr
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
      end if
!
      call makeref_formol(imol)
!      call makeref_polym(imol)
      if (dc_di(imol)%maxntor.gt.0) then
        call IMD_prealign(imol,skip_frz)
        call makexyz_formol(imol)
      end if
!
!     finally: increment coordinates
!              note rotation is done first as it relies on comm
      if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
        call rotxyzm(imol,aone)
      end if
      call transxyzm(imol,aone,mshfs(:,imol))
!     due to torsional moves we need to recompute rigid-body coordinates
      call update_rigidm(imol)
!     update grid association if necessary
      if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(rs)
        end do
      end if
    end do
!
!   finally, in NPT(E), update box
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (bnd_shape.eq.2) then
        if (pstat%flag.eq.1) then
          netf= bnd_fr - 4.0*PI*bnd_params(5)*extpress/u_dyn_virconv
          bnd_v= bnd_v + 0.5*u_dyn_fconv*dyn_dt*netf/pstat%params(1)
          bnd_params(4) = bnd_params(4) + dyn_dt*bnd_v
          call update_bound(afalse)
        else
          write(ilog,*) 'Fatal. Encountered unsupported manostat for&
 & chosen box and ensemble (manostat-flag: ',pstat%flag,').'
          call fexit()
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary shape&
 & for chosen ensemble (Flag: ',ens%flag,').'
        call fexit()
      end if
    end if
!
  else
!
    do imol=1,nmol
!     zero out internal coordinate forces and back up inertia at t_1
      dc_di(imol)%f(:) = 0.0
      dc_di(imol)%olddat(:,1) = dc_di(imol)%im(:)
    end do
!
    call System_Clock(tt1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    ens%insR(10) = ens%insU ! backup
    ens%insU = esave
    boxovol = ens%insV
    call System_Clock(tt2)
    call cart2int_f(skip_frz,azero)
    if (fudge_mass.EQV..true.) call fudge_masses()
    time_energy = time_energy + 1.0*(tt2 - tt1)
!
    if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
!     get the temperature re-scaling factor from thermostat
      call thermostat(tscs)
    else
      tscs(:) = 1.0
    end if
!
!   if inertia are variable and we want to guess I(t1.5), the following is used
    if (dyn_integrator_ops(1).gt.0) then
!     backup entire Z matrix
      call makeref_zmat_gl()
!     now loop over all molecules processing everything but translational motion
      do imol=1,nmol
        tsc = tscs(tstat%molgrp(imol))
        dc_di(imol)%v(1:3) = tsc*dc_di(imol)%v(1:3)
!
        if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
          do i=-2,dc_di(imol)%maxntor
            if (i.le.0) then
              ttc = i+6 ! RB rotation
            else
              ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
            end if
            if (ttc.le.0) cycle
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              dc_di(imol)%incr(ttc) = 0.0
            else
              dc_di(imol)%v(ttc) = tsc*dc_di(imol)%v(ttc)
              t0 = 0.5*dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
!             use inertia at t1 and previous guess at t0.5
              t1 = sqrt(dc_di(imol)%olddat(ttc,2)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
              if (dyn_integrator.eq.2) then
                dc_di(imol)%incr(ttc) = 0.5*dyn_dt*t1 ! approximate
              else
                t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*dc_di(imol)%im(ttc)*dc_di(imol)%olddat(ttc,2) + &
 & 4.0*t0*dc_di(imol)%v(ttc)*dc_di(imol)%im(ttc)
                if (t2.gt.0.0) then ! solve quadratic equation exactly, decide on appropriate solution via similarity to approx.
                  t2 = sqrt(t2)
                  if (abs(t1-(t0 - t2)/(2.0*dc_di(imol)%im(ttc))).gt.abs(t1-(t0 + t2)/(2.0*dc_di(imol)%im(ttc)))) then
                    dc_di(imol)%incr(ttc) = 0.5*dyn_dt*(t0 + t2)/(2.0*dc_di(imol)%im(ttc))
                  else
                    dc_di(imol)%incr(ttc) = 0.5*dyn_dt*(t0 - t2)/(2.0*dc_di(imol)%im(ttc))
                  end if
                else
                  dc_di(imol)%incr(ttc) =  0.5*dyn_dt*t1
                end if
              end if
            end if
          end do
!         now transfer
          cur_rot(1:3) = dc_di(imol)%incr(4:6)/RADIAN
!         note that this will edit the Z-matrix, but not alter ztorpr
          if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,aone)
        end if
        call makeref_formol(imol)
        if (dc_di(imol)%maxntor.gt.0) then
          call IMD_prealign(imol,skip_frz)
          call makexyz_formol(imol)
        end if
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rotxyzm(imol,aone)
        end if
      end do
!
!     back up I(t0.5) and recompute inertia for t1.5
      do imol=1,nmol
        dc_di(imol)%olddat(:,3) = dc_di(imol)%olddat(:,2)
      end do
      call cart2int_I(azero) ! must use xyz only
!     restore Z matrix and transfer
      call getref_zmat_gl()
      call zmatfyc() ! zmatfyc2() would also update phish/psish
!     restore coordinates (cheaper than rebuilding of course)
      do imol=1,nmol
        call getref_formol(imol)
      end do
    end if
!
!   now loop over all molecules
    do imol=1,nmol
!
!     T-group specific T-rescaling
      if (dyn_integrator_ops(1).eq.0) then
        tsc = tscs(tstat%molgrp(imol))
      else
        tsc = 1.0
      end if
!
!     first: center-of-mass translation (linear motion)
      do j=1,3
!       increment velocity at time t - 1/2dt to t + 1/2dt
        if (dc_di(imol)%frz(j).EQV..true.) then
          cur_trans(j) = 0.0
          cycle
        end if
        dc_di(imol)%v(j) = tsc*dc_di(imol)%v(j) + u_dyn_fconv*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j)
!       increment position t to t + dt with staggered velocity at t + 1/2dt
        cur_trans(j) = dyn_dt*dc_di(imol)%v(j)
      end do
!     now all rotational motion (including RB) 
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
            dc_di(imol)%v(ttc) = 0.0
            cycle
          end if
!
!         thermostat (tsc is 1.0 if already performed)
          dc_di(imol)%v(ttc) = tsc*dc_di(imol)%v(ttc)
!
!         the inertia correction is lagging in time which causes drainage of Cartesian K.E.
          t0 = dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r          
          if (dyn_integrator_ops(1).le.0) then
            t1 = sqrt(dc_di(imol)%olddat(ttc,1)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
            if (dyn_integrator.eq.2) then
              dc_di(imol)%v(ttc) = t1
            else
              t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*dc_di(imol)%im(ttc)*dc_di(imol)%olddat(ttc,1) + &
 & 4.0*t0*dc_di(imol)%v(ttc)*dc_di(imol)%im(ttc)
              if (t2.gt.0.0) then
                t2 = sqrt(t2)
                if (abs(t1-(t0 - t2)/(2.0*dc_di(imol)%im(ttc))).gt.abs(t1-(t0 + t2)/(2.0*dc_di(imol)%im(ttc)))) then
                  dc_di(imol)%v(ttc) = (t0 + t2)/(2.0*dc_di(imol)%im(ttc))
                else
                  dc_di(imol)%v(ttc) = (t0 - t2)/(2.0*dc_di(imol)%im(ttc))
                end if
              else
                dc_di(imol)%v(ttc) = t1
              end if
            end if
          else
!           with the guessed correction, we still have the option to increment velocity in stepwise fashion
!           to better account for rapidly changing inertia
            t0 = t0/(1.0*dyn_integrator_ops(1))
            do j=1,dyn_integrator_ops(1)
              frac1 = (1.0*(j-1))/(1.0*dyn_integrator_ops(1))
              frac2 = (1.0*j)/(1.0*dyn_integrator_ops(1))
              if (frac1.le.0.5) then
                ilow = dc_di(imol)%olddat(ttc,3) + 2.0*frac1*(dc_di(imol)%im(ttc)-dc_di(imol)%olddat(ttc,3))
              else
                ilow = dc_di(imol)%im(ttc) + 2.0*(frac1-0.5)*(dc_di(imol)%olddat(ttc,2)-dc_di(imol)%im(ttc))
              end if
              if (frac2.le.0.5) then
                ihigh = dc_di(imol)%olddat(ttc,3) + 2.0*frac2*(dc_di(imol)%im(ttc)-dc_di(imol)%olddat(ttc,3))
              else
                ihigh = dc_di(imol)%im(ttc) + 2.0*(frac2-0.5)*(dc_di(imol)%olddat(ttc,2)-dc_di(imol)%im(ttc))
              end if
              t1 = sqrt(ilow/ihigh)*dc_di(imol)%v(ttc) + t0/ihigh
              if (dyn_integrator.eq.2) then
                dc_di(imol)%v(ttc) = t1 
              else
                t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*ihigh*ilow + 4.0*t0*dc_di(imol)%v(ttc)*ihigh
                if (t2.gt.0.0) then
                  t2 = sqrt(t2)
                  frac1 = dc_di(imol)%v(ttc)
                  if (abs(t1-(t0 - t2)/(2.0*ihigh)).gt.abs(t1-(t0 + t2)/(2.0*ihigh))) then
                    dc_di(imol)%v(ttc) = (t0 + t2)/(2.0*ihigh)
                  else
                    dc_di(imol)%v(ttc) = (t0 - t2)/(2.0*ihigh)
                  end if
                else
!                  frac1 = 0.5*(t1*t1*ihigh - (dc_di(imol)%v(ttc)**2)*ilow - (t1+dc_di(imol)%v(ttc))*t0)/(u_dyn_fconv_r)
!                  ens%avgR(10) = ens%avgR(10) + frac1
!                   write(*,*) ens%avgR(10)/(nstep*dyn_dt)
                  dc_di(imol)%v(ttc) = t1
                end if
              end if
            end do
          end if
!
          dc_di(imol)%incr(ttc) = dyn_dt*dc_di(imol)%v(ttc) ! in degrees
          if (abs(dc_di(imol)%incr(ttc)).gt.36.0) then
            fo_wrncnt(10) = fo_wrncnt(10) + 1
            if (fo_wrncnt(10).eq.fo_wrnlmt(10)) then
              if (i.le.0) then
                write(ilog,*) 'Warning. Extreme velocity for rigid rotation of molecule ',imol,' suggests unstable simulation.'
              else
                write(ilog,*) 'Warning. Extreme velocity for torsion defined by atom ',dc_di(imol)%recurs(i,1),' &
 &suggests unstable simulation.'
              end if
              write(ilog,*) 'This was warning #',fo_wrncnt(10),' of this type not all of which may be displayed.'
              if (10.0*fo_wrnlmt(10).gt.0.5*HUGE(fo_wrnlmt(10))) then
                fo_wrncnt(10) = 0
              else
                fo_wrnlmt(10) = fo_wrnlmt(10)*10
              end if
            end if
          end if
        end do
!       now transfer
        cur_rot(1:3) = dc_di(imol)%incr(4:6)/RADIAN
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
      end if
!
      call makeref_formol(imol)
!      call makeref_polym(imol)
      if (dc_di(imol)%maxntor.gt.0) then
        call IMD_prealign(imol,skip_frz)
        call makexyz_formol(imol)
      end if
!
!     finally: increment coordinates
!              note rotation is done first as it relies on comm
      if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
        call rotxyzm(imol,aone)
      end if
      call transxyzm(imol,aone,mshfs(:,imol))
!     due to torsional moves we need to recompute rigid-body coordinates
      call update_rigidm(imol)
!     update grid association if necessary
      if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
        do rs=rsmol(imol,1),rsmol(imol,2)
          call updateresgp(rs)
        end do
      end if
!
    end do
!
!   finally, in NPT(E), update box
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      tsc = tscs(tstat%molgrp(nmol+1))
      if (bnd_shape.eq.2) then
        if (pstat%flag.eq.1) then
          netf= bnd_fr - 4.0*PI*bnd_params(5)*extpress/u_dyn_virconv
          bnd_v= tsc*bnd_v + u_dyn_fconv*dyn_dt*netf/pstat%params(1)
          bnd_params(4) = bnd_params(4) + dyn_dt*bnd_v
          call update_bound(afalse)
        else
          write(ilog,*) 'Fatal. Encountered unsupported manostat for&
 & chosen box and ensemble (manostat-flag: ',pstat%flag,').'
          call fexit()
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary shape&
 & for chosen ensemble (Flag: ',ens%flag,').'
        call fexit()
      end if
    end if
!
  end if
!
  call get_ensv(boxovol,azero)
  call drift_removal(frac1,frac2,azero)
  call prt_curens(istep,boxovol,frac1,frac2)
!
  deallocate(mshfs)
!
  end
!
!------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine int_mdmove_threads(istep,ndump,forcefirst,tpi)
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
  use zmatrix
  use threads
!
  implicit none
!
  integer, INTENT(IN):: tpi,istep,ndump
  logical, INTENT(IN):: forcefirst
!
  integer imol,j,aone,azero,ttc,rs,i,ixx
  RTYPE t1,t2,t0,ilow,ihigh,frac1,frac2
  RTYPE netf,tscs(tstat%n_tgrps),boxovol,tsc
  RTYPE, ALLOCATABLE:: mshfs(:,:)
  logical atrue,afalse,jobflags(4)
  integer(KIND=8) tt1,tt2,ttimer
!
  allocate(mshfs(3,nmol))
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
  atrue = .true.
  afalse = .false.
  azero = 0
  aone = 1
!$OMP SINGLE
  nstep = nstep + 1
  mvcnt%nmd = mvcnt%nmd + 1
!$OMP END SINGLE
!
! initialize velocities
!
  if ((istep.eq.1).OR.(forcefirst.EQV..true.)) then
!
!$OMP SINGLE
    do imol=1,nmol
!     sanity check in first step against (partially) unsupported 2-atom molecules
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
      end if
      call update_rigidm(imol)
!     zero out internal coordinate forces
      dc_di(imol)%f(:) = 0.0
      dc_di(imol)%olddat(:,1) = dc_di(imol)%im(:)
      dc_di(imol)%olddat(:,2) = dc_di(imol)%im(:) ! safe side
      dc_di(imol)%v(:) = 0.0
    end do
!$OMP END SINGLE
!
    if (tpi.eq.1) call System_Clock(tt1)
    call force3_threads(esterms,esterms_tr,esterms_lr,forcefirst,esave)
    if (tpi.eq.1) then
      call System_Clock(tt2)
      time_energy = time_energy + 1.0*(tt2 - tt1)
    end if
    call cart2int_f(skip_frz,tpi)
!$OMP BARRIER
!$OMP SINGLE
    ens%insR(10) = ens%insU ! backup
    ens%insU = esave
    boxovol = ens%insV
    if (fudge_mass.EQV..true.) call fudge_masses()
!
!   note that the initialized velocities are exactly matched to the temperature
!   request -> no need to use thermostat in first step
!   prior call to cart2int_f or similar is necessary to populate inertia initially
!$OMP END SINGLE
    call randomize_velocities(tpi)
!$OMP BARRIER
!$OMP SINGLE
!
!   the velocities just need to be incremented by a half-step without having access to any
!   past information, coordinates are advanced a full step
    do imol=1,nmol
!
!     first: center-of-mass translation (linear motion)
      do j=1,3
!       increment velocity at time t_0 to t_0 + 1/2dt
        if (dc_di(imol)%frz(j).EQV..true.) then
          dc_di(imol)%incr(j) = 0.0
          cycle
        end if
        dc_di(imol)%v(j) = dc_di(imol)%v(j) + u_dyn_fconv*0.5*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j)
        dc_di(imol)%incr(j) = dyn_dt*dc_di(imol)%v(j)
      end do
!
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%v(ttc) = 0.0
          else
            t0 = 0.5*dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
            dc_di(imol)%v(ttc) = dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
          end if
          dc_di(imol)%incr(ttc) = dyn_dt*dc_di(imol)%v(ttc) ! in degrees
        end do
!       note that this will edit the Z-matrix and populate ztorpr
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
      end if
!
      call makeref_formol(imol)
!      call makeref_polym(imol)
      if (dc_di(imol)%maxntor.gt.0) then
        call IMD_prealign(imol,skip_frz)
        call makexyz_formol(imol)
      end if
!
!     finally: increment coordinates
!              note rotation is done first as it relies on comm
      call rbcxyzm(imol,azero,atrue)
!     due to torsional moves we need to recompute rigid-body coordinates
      call update_rigidm(imol)
    end do
!
!   finally, in NPT(E), update box
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (bnd_shape.eq.2) then
        if (pstat%flag.eq.1) then
          netf= bnd_fr - 4.0*PI*bnd_params(5)*extpress/u_dyn_virconv
          bnd_v= bnd_v + 0.5*u_dyn_fconv*dyn_dt*netf/pstat%params(1)
          bnd_params(4) = bnd_params(4) + dyn_dt*bnd_v
          call update_bound(afalse)
        else
          write(ilog,*) 'Fatal. Encountered unsupported manostat for&
 & chosen box and ensemble (manostat-flag: ',pstat%flag,').'
          call fexit()
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary shape&
 & for chosen ensemble (Flag: ',ens%flag,').'
        call fexit()
      end if
    end if
!$OMP END SINGLE
!
  else
!
    do imol=thr_limits(61,tpi),thr_limits(62,tpi)
!     zero out internal coordinate forces and back up inertia at t_1
      dc_di(imol)%f(:) = 0.0
      dc_di(imol)%olddat(:,1) = dc_di(imol)%im(:)
    end do
!   coordinate and Z matrix backup
    xref(thr_limits(1,tpi):thr_limits(2,tpi)) = x(thr_limits(1,tpi):thr_limits(2,tpi))
    yref(thr_limits(1,tpi):thr_limits(2,tpi)) = y(thr_limits(1,tpi):thr_limits(2,tpi))
    zref(thr_limits(1,tpi):thr_limits(2,tpi)) = z(thr_limits(1,tpi):thr_limits(2,tpi))
    if (dyn_integrator_ops(1).gt.0) then
      ztorpr(thr_limits(1,tpi):thr_limits(2,tpi)) = ztor(thr_limits(1,tpi):thr_limits(2,tpi))
      bangpr(thr_limits(1,tpi):thr_limits(2,tpi)) = bang(thr_limits(1,tpi):thr_limits(2,tpi))
    end if
!$OMP BARRIER
!
!
    if (tpi.eq.1) call System_Clock(tt1)
    call force3_threads(esterms,esterms_tr,esterms_lr,forcefirst,esave)
    if (tpi.eq.1) then
      call System_Clock(tt2)
      time_energy = time_energy + 1.0*(tt2 - tt1)
    end if
!    call System_Clock(count=thr_timings(66,tpi))
    call cart2int_f(skip_frz,tpi)
!    call System_Clock(count=thr_timings(67,tpi))
!    write(*,*) tpi,1.0e-3*(thr_timings(67,tpi)-thr_timings(66,tpi))
!$OMP BARRIER
!$OMP SINGLE
    ens%insR(10) = ens%insU ! backup
    ens%insU = esave
    boxovol = ens%insV
    if (fudge_mass.EQV..true.) call fudge_masses()
!
    if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
!     get the temperature re-scaling factor from thermostat
      call thermostat(tscs)
    else
      tscs(:) = 1.0
    end if
!$OMP END SINGLE COPYPRIVATE(tscs)
!
!   if inertia are variable and we want to guess I(t1.5), the following is used
    if (dyn_integrator_ops(1).gt.0) then
!     now loop over all molecules processing everything but translational motion
!     this first loop is only for molecules processed by only a single thread
      do imol=thr_limits(61,tpi),thr_limits(62,tpi)
        if (thr_mlgix(imol).gt.0) cycle
        tsc = tscs(tstat%molgrp(imol))
        dc_di(imol)%v(1:3) = tsc*dc_di(imol)%v(1:3)
        dc_di(imol)%incr(1:3) = 0.0
!
        if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
          do i=-2,dc_di(imol)%maxntor
            if (i.le.0) then
              ttc = i+6 ! RB rotation
            else
              ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
            end if
            if (ttc.le.0) cycle
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              dc_di(imol)%incr(ttc) = 0.0
            else
              dc_di(imol)%v(ttc) = tsc*dc_di(imol)%v(ttc)
              t0 = 0.5*dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
!             use inertia at t1 and previous guess at t0.5
              t1 = sqrt(dc_di(imol)%olddat(ttc,2)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
              if (dyn_integrator.eq.2) then
                dc_di(imol)%incr(ttc) = 0.5*dyn_dt*t1 ! approximate
              else
                t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*dc_di(imol)%im(ttc)*dc_di(imol)%olddat(ttc,2) + &
 & 4.0*t0*dc_di(imol)%v(ttc)*dc_di(imol)%im(ttc)
                if (t2.gt.0.0) then ! solve quadratic equation exactly, decide on appropriate solution via similarity to approx.
                  t2 = sqrt(t2)
                  if (abs(t1-(t0 - t2)/(2.0*dc_di(imol)%im(ttc))).gt.abs(t1-(t0 + t2)/(2.0*dc_di(imol)%im(ttc)))) then
                    dc_di(imol)%incr(ttc) = 0.5*dyn_dt*(t0 + t2)/(2.0*dc_di(imol)%im(ttc))
                  else
                    dc_di(imol)%incr(ttc) = 0.5*dyn_dt*(t0 - t2)/(2.0*dc_di(imol)%im(ttc))
                  end if
                else
                  dc_di(imol)%incr(ttc) =  0.5*dyn_dt*t1
                end if
              end if
            end if
          end do
!         note that this will edit the Z-matrix, but not alter ztorpr
          if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,aone)
        end if
        if (dc_di(imol)%maxntor.gt.0) then
          call IMD_prealign(imol,skip_frz)
          call makexyz_formol(imol)
        end if
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rbcxyzm(imol,azero,afalse)
        end if
      end do
!     the same for internally parallelized molecules
      do ixx=1,nmlgs
        imol = mlg_limits(ixx,5,tpi)
        tsc = tscs(tstat%molgrp(imol))
!
        if (tpi.eq.1) then
          dc_di(imol)%v(1:3) = tsc*dc_di(imol)%v(1:3)
          dc_di(imol)%incr(1:3) = 0.0
        end if
!
        if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
          do i=mlg_limits(ixx,13,tpi),mlg_limits(ixx,14,tpi) ! -2,dc_di(imol)%maxntor
            if (i.le.0) then
              ttc = i+6 ! RB rotation
            else
              ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
            end if
            if (ttc.le.0) cycle
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              dc_di(imol)%incr(ttc) = 0.0
            else
              dc_di(imol)%v(ttc) = tsc*dc_di(imol)%v(ttc)
              t0 = 0.5*dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
!             use inertia at t1 and previous guess at t0.5
              t1 = sqrt(dc_di(imol)%olddat(ttc,2)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
              if (dyn_integrator.eq.2) then
                dc_di(imol)%incr(ttc) = 0.5*dyn_dt*t1 ! approximate
              else
                t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*dc_di(imol)%im(ttc)*dc_di(imol)%olddat(ttc,2) + &
 & 4.0*t0*dc_di(imol)%v(ttc)*dc_di(imol)%im(ttc)
                if (t2.gt.0.0) then ! solve quadratic equation exactly, decide on appropriate solution via similarity to approx.
                  t2 = sqrt(t2)
                  if (abs(t1-(t0 - t2)/(2.0*dc_di(imol)%im(ttc))).gt.abs(t1-(t0 + t2)/(2.0*dc_di(imol)%im(ttc)))) then
                    dc_di(imol)%incr(ttc) = 0.5*dyn_dt*(t0 + t2)/(2.0*dc_di(imol)%im(ttc))
                  else
                    dc_di(imol)%incr(ttc) = 0.5*dyn_dt*(t0 - t2)/(2.0*dc_di(imol)%im(ttc))
                  end if
                else
                  dc_di(imol)%incr(ttc) =  0.5*dyn_dt*t1
                end if
              end if
            end if
          end do
!         note that this will edit the Z-matrix, but not alter ztorpr
!$OMP BARRIER
!$OMP SINGLE
          if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,aone) ! bottleneck
!$OMP END SINGLE
        end if
        if (dc_di(imol)%maxntor.gt.0) then
!$OMP SINGLE
          call IMD_prealign(imol,skip_frz) ! bottleneck
!$OMP END SINGLE
          call makexyz_threads(ixx,tpi)
        end if
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rbcxyzm(ixx,tpi,afalse) ! ixx is translated inside rbcxyzm if tpi is > 0
        end if
      end do
!$OMP BARRIER
!
!     back up I(t0.5) and recompute inertia for t1.5
      do imol=thr_limits(61,tpi),thr_limits(62,tpi)
        dc_di(imol)%olddat(:,3) = dc_di(imol)%olddat(:,2)
      end do
!$OMP BARRIER
      call cart2int_I(tpi) ! must use xyz only
!$OMP BARRIER
      ztor(thr_limits(1,tpi):thr_limits(2,tpi)) = ztorpr(thr_limits(1,tpi):thr_limits(2,tpi))
      bang(thr_limits(1,tpi):thr_limits(2,tpi)) = bangpr(thr_limits(1,tpi):thr_limits(2,tpi))
!$OMP BARRIER
      call zmatfyc_threads(tpi,aone)
      x(thr_limits(1,tpi):thr_limits(2,tpi)) = xref(thr_limits(1,tpi):thr_limits(2,tpi))
      y(thr_limits(1,tpi):thr_limits(2,tpi)) = yref(thr_limits(1,tpi):thr_limits(2,tpi))
      z(thr_limits(1,tpi):thr_limits(2,tpi)) = zref(thr_limits(1,tpi):thr_limits(2,tpi))
!$OMP BARRIER
    end if
!
!   first loop over all molecules that are treated by one thread only
    if (thr_dlb(12,1).gt.0) then
      if (tpi.eq.1) thr_dlb(12,2) = thr_dlb(12,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(27,tpi) = thr_timings(27,tpi) + ttimer
    end if
    do imol=thr_limits(61,tpi),thr_limits(62,tpi)
      if (thr_mlgix(imol).gt.0) cycle
!
!     T-group specific T-rescaling
      if (dyn_integrator_ops(1).eq.0) then
        tsc = tscs(tstat%molgrp(imol))
      else
        tsc = 1.0
      end if
!
!     first: center-of-mass translation (linear motion)
      do j=1,3
!       increment velocity at time t - 1/2dt to t + 1/2dt
        if (dc_di(imol)%frz(j).EQV..true.) then
          dc_di(imol)%incr(j) = 0.0
          cycle
        end if
        dc_di(imol)%v(j) = tsc*dc_di(imol)%v(j) + u_dyn_fconv*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j)
!       increment position t to t + dt with staggered velocity at t + 1/2dt
        dc_di(imol)%incr(j) = dyn_dt*dc_di(imol)%v(j)
      end do
!     now all rotational motion (including RB) 
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
            dc_di(imol)%v(ttc) = 0.0
            cycle
          end if
!
!         thermostat (tsc is 1.0 if already performed)
          dc_di(imol)%v(ttc) = tsc*dc_di(imol)%v(ttc)
!
!         the inertia correction is lagging in time which causes drainage of Cartesian K.E.
          t0 = dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r          
          if (dyn_integrator_ops(1).le.0) then
            t1 = sqrt(dc_di(imol)%olddat(ttc,1)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
            if (dyn_integrator.eq.2) then
              dc_di(imol)%v(ttc) = t1
            else
              t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*dc_di(imol)%im(ttc)*dc_di(imol)%olddat(ttc,1) + &
 & 4.0*t0*dc_di(imol)%v(ttc)*dc_di(imol)%im(ttc)
              if (t2.gt.0.0) then
                t2 = sqrt(t2)
                if (abs(t1-(t0 - t2)/(2.0*dc_di(imol)%im(ttc))).gt.abs(t1-(t0 + t2)/(2.0*dc_di(imol)%im(ttc)))) then
                  dc_di(imol)%v(ttc) = (t0 + t2)/(2.0*dc_di(imol)%im(ttc))
                else
                  dc_di(imol)%v(ttc) = (t0 - t2)/(2.0*dc_di(imol)%im(ttc))
                end if
              else
                dc_di(imol)%v(ttc) = t1
              end if
            end if
          else
!           with the guessed correction, we still have the option to increment velocity in stepwise fashion
!           to better account for rapidly changing inertia
            t0 = t0/(1.0*dyn_integrator_ops(1))
            do j=1,dyn_integrator_ops(1)
              frac1 = (1.0*(j-1))/(1.0*dyn_integrator_ops(1))
              frac2 = (1.0*j)/(1.0*dyn_integrator_ops(1))
              if (frac1.le.0.5) then
                ilow = dc_di(imol)%olddat(ttc,3) + 2.0*frac1*(dc_di(imol)%im(ttc)-dc_di(imol)%olddat(ttc,3))
              else
                ilow = dc_di(imol)%im(ttc) + 2.0*(frac1-0.5)*(dc_di(imol)%olddat(ttc,2)-dc_di(imol)%im(ttc))
              end if
              if (frac2.le.0.5) then
                ihigh = dc_di(imol)%olddat(ttc,3) + 2.0*frac2*(dc_di(imol)%im(ttc)-dc_di(imol)%olddat(ttc,3))
              else
                ihigh = dc_di(imol)%im(ttc) + 2.0*(frac2-0.5)*(dc_di(imol)%olddat(ttc,2)-dc_di(imol)%im(ttc))
              end if
              t1 = sqrt(ilow/ihigh)*dc_di(imol)%v(ttc) + t0/ihigh
              if (dyn_integrator.eq.2) then
                dc_di(imol)%v(ttc) = t1 
              else
                t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*ihigh*ilow + 4.0*t0*dc_di(imol)%v(ttc)*ihigh
                if (t2.gt.0.0) then
                  t2 = sqrt(t2)
                  frac1 = dc_di(imol)%v(ttc)
                  if (abs(t1-(t0 - t2)/(2.0*ihigh)).gt.abs(t1-(t0 + t2)/(2.0*ihigh))) then
                    dc_di(imol)%v(ttc) = (t0 + t2)/(2.0*ihigh)
                  else
                    dc_di(imol)%v(ttc) = (t0 - t2)/(2.0*ihigh)
                  end if
                else
!                  frac1 = 0.5*(t1*t1*ihigh - (dc_di(imol)%v(ttc)**2)*ilow - (t1+dc_di(imol)%v(ttc))*t0)/(u_dyn_fconv_r)
!                  ens%avgR(10) = ens%avgR(10) + frac1
!                   write(*,*) ens%avgR(10)/(nstep*dyn_dt)
                  dc_di(imol)%v(ttc) = t1
                end if
              end if
            end do
          end if
!
          dc_di(imol)%incr(ttc) = dyn_dt*dc_di(imol)%v(ttc) ! in degrees
          if (abs(dc_di(imol)%incr(ttc)).gt.36.0) then
!$OMP CRITICAL(IMD_WARNING)
            fo_wrncnt(10) = fo_wrncnt(10) + 1
            if (fo_wrncnt(10).eq.fo_wrnlmt(10)) then
              if (i.le.0) then
                write(ilog,*) 'Warning. Extreme velocity for rigid rotation of molecule ',imol,' suggests unstable simulation.'
              else
                write(ilog,*) 'Warning. Extreme velocity for torsion defined by atom ',dc_di(imol)%recurs(i,1),' &
 &suggests unstable simulation.'
              end if
              write(ilog,*) 'This was warning #',fo_wrncnt(10),' of this type not all of which may be displayed.'
              if (10.0*fo_wrnlmt(10).gt.0.5*HUGE(fo_wrnlmt(10))) then
                fo_wrncnt(10) = 0
              else
                fo_wrnlmt(10) = fo_wrnlmt(10)*10
              end if
            end if
!$OMP END CRITICAL(IMD_WARNING)
          end if
        end do
!       now transfer
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
      end if
!
      if (dc_di(imol)%maxntor.gt.0) then
        call IMD_prealign(imol,skip_frz)
        call makexyz_formol(imol)
      end if
      call rbcxyzm(imol,azero,atrue)
!     due to torsional moves we need to recompute rigid-body coordinates
      call update_rigidm(imol)
    end do
    if (thr_dlb(12,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(28,tpi) = thr_timings(28,tpi) + ttimer
    end if
!
!   now the same for any internally parallelized molecules
    if (nmlgs.gt.0) then
      jobflags(:) = .false.
      jobflags(1) = .true.
    end if
    do ixx=1,nmlgs
      imol = mlg_limits(ixx,5,tpi)
!
!     T-group specific T-rescaling
      if (dyn_integrator_ops(1).eq.0) then
        tsc = tscs(tstat%molgrp(imol))
      else
        tsc = 1.0
      end if
!
!     first: center-of-mass translation (linear motion)
      if (tpi.eq.1) then
        do j=1,3
!         increment velocity at time t - 1/2dt to t + 1/2dt
          if (dc_di(imol)%frz(j).EQV..true.) then
            dc_di(imol)%incr(j) = 0.0
            cycle
          end if
          dc_di(imol)%v(j) = tsc*dc_di(imol)%v(j) + u_dyn_fconv*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j)
!         increment position t to t + dt with staggered velocity at t + 1/2dt
          dc_di(imol)%incr(j) = dyn_dt*dc_di(imol)%v(j)
        end do
      end if
!     now all rotational motion (including RB) 
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=mlg_limits(ixx,13,tpi),mlg_limits(ixx,14,tpi) ! -2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
            dc_di(imol)%v(ttc) = 0.0
            cycle
          end if
!
!         thermostat (tsc is 1.0 if already performed)
          dc_di(imol)%v(ttc) = tsc*dc_di(imol)%v(ttc)
!
!         the inertia correction is lagging in time which causes drainage of Cartesian K.E.
          t0 = dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r          
          if (dyn_integrator_ops(1).le.0) then
            t1 = sqrt(dc_di(imol)%olddat(ttc,1)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + t0/dc_di(imol)%im(ttc)
            if (dyn_integrator.eq.2) then
              dc_di(imol)%v(ttc) = t1
            else
              t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*dc_di(imol)%im(ttc)*dc_di(imol)%olddat(ttc,1) + &
 & 4.0*t0*dc_di(imol)%v(ttc)*dc_di(imol)%im(ttc)
              if (t2.gt.0.0) then
                t2 = sqrt(t2)
                if (abs(t1-(t0 - t2)/(2.0*dc_di(imol)%im(ttc))).gt.abs(t1-(t0 + t2)/(2.0*dc_di(imol)%im(ttc)))) then
                  dc_di(imol)%v(ttc) = (t0 + t2)/(2.0*dc_di(imol)%im(ttc))
                else
                  dc_di(imol)%v(ttc) = (t0 - t2)/(2.0*dc_di(imol)%im(ttc))
                end if
              else
                dc_di(imol)%v(ttc) = t1
              end if
            end if
          else
!           with the guessed correction, we still have the option to increment velocity in stepwise fashion
!           to better account for rapidly changing inertia
            t0 = t0/(1.0*dyn_integrator_ops(1))
            do j=1,dyn_integrator_ops(1)
              frac1 = (1.0*(j-1))/(1.0*dyn_integrator_ops(1))
              frac2 = (1.0*j)/(1.0*dyn_integrator_ops(1))
              if (frac1.le.0.5) then
                ilow = dc_di(imol)%olddat(ttc,3) + 2.0*frac1*(dc_di(imol)%im(ttc)-dc_di(imol)%olddat(ttc,3))
              else
                ilow = dc_di(imol)%im(ttc) + 2.0*(frac1-0.5)*(dc_di(imol)%olddat(ttc,2)-dc_di(imol)%im(ttc))
              end if
              if (frac2.le.0.5) then
                ihigh = dc_di(imol)%olddat(ttc,3) + 2.0*frac2*(dc_di(imol)%im(ttc)-dc_di(imol)%olddat(ttc,3))
              else
                ihigh = dc_di(imol)%im(ttc) + 2.0*(frac2-0.5)*(dc_di(imol)%olddat(ttc,2)-dc_di(imol)%im(ttc))
              end if
              t1 = sqrt(ilow/ihigh)*dc_di(imol)%v(ttc) + t0/ihigh
              if (dyn_integrator.eq.2) then
                dc_di(imol)%v(ttc) = t1 
              else
                t2 = t0**2 + 4.0*(dc_di(imol)%v(ttc)**2)*ihigh*ilow + 4.0*t0*dc_di(imol)%v(ttc)*ihigh
                if (t2.gt.0.0) then
                  t2 = sqrt(t2)
                  frac1 = dc_di(imol)%v(ttc)
                  if (abs(t1-(t0 - t2)/(2.0*ihigh)).gt.abs(t1-(t0 + t2)/(2.0*ihigh))) then
                    dc_di(imol)%v(ttc) = (t0 + t2)/(2.0*ihigh)
                  else
                    dc_di(imol)%v(ttc) = (t0 - t2)/(2.0*ihigh)
                  end if
                else
!                  frac1 = 0.5*(t1*t1*ihigh - (dc_di(imol)%v(ttc)**2)*ilow - (t1+dc_di(imol)%v(ttc))*t0)/(u_dyn_fconv_r)
!                  ens%avgR(10) = ens%avgR(10) + frac1
!                   write(*,*) ens%avgR(10)/(nstep*dyn_dt)
                  dc_di(imol)%v(ttc) = t1
                end if
              end if
            end do
          end if
!
          dc_di(imol)%incr(ttc) = dyn_dt*dc_di(imol)%v(ttc) ! in degrees
          if (abs(dc_di(imol)%incr(ttc)).gt.36.0) then
!$OMP CRITICAL(IMD_WARNING)
            fo_wrncnt(10) = fo_wrncnt(10) + 1
            if (fo_wrncnt(10).eq.fo_wrnlmt(10)) then
              if (i.le.0) then
                write(ilog,*) 'Warning. Extreme velocity for rigid rotation of molecule ',imol,' suggests unstable simulation.'
              else
                write(ilog,*) 'Warning. Extreme velocity for torsion defined by atom ',dc_di(imol)%recurs(i,1),' &
 &suggests unstable simulation.'
              end if
              write(ilog,*) 'This was warning #',fo_wrncnt(10),' of this type not all of which may be displayed.'
              if (10.0*fo_wrnlmt(10).gt.0.5*HUGE(fo_wrnlmt(10))) then
                fo_wrncnt(10) = 0
              else
                fo_wrnlmt(10) = fo_wrnlmt(10)*10
              end if
            end if
!$OMP END CRITICAL(IMD_WARNING)
          end if
        end do
!$OMP BARRIER
!$OMP SINGLE
!       now transfer
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
!$OMP END SINGLE
      end if
!
      if (dc_di(imol)%maxntor.gt.0) then
!$OMP SINGLE
        call IMD_prealign(imol,skip_frz)
!$OMP END SINGLE
        call makexyz_threads(ixx,tpi) ! ends with a barrier
      end if
!
      call rbcxyzm(ixx,tpi,atrue) ! ends with an implied barrier
!     due to torsional moves we need to recompute rigid-body coordinates
      call molops_threads(ixx,tpi,jobflags)
!
    end do
!
!$OMP BARRIER
!   finally, in NPT(E), update box
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!$OMP SINGLE
      tsc = tscs(tstat%molgrp(nmol+1))
      if (bnd_shape.eq.2) then
        if (pstat%flag.eq.1) then
          netf= bnd_fr - 4.0*PI*bnd_params(5)*extpress/u_dyn_virconv
          bnd_v= tsc*bnd_v + u_dyn_fconv*dyn_dt*netf/pstat%params(1)
          bnd_params(4) = bnd_params(4) + dyn_dt*bnd_v
          call update_bound(afalse)
        else
          write(ilog,*) 'Fatal. Encountered unsupported manostat for&
 & chosen box and ensemble (manostat-flag: ',pstat%flag,').'
          call fexit()
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary shape&
 & for chosen ensemble (Flag: ',ens%flag,').'
        call fexit()
      end if
!$OMP END SINGLE
    end if
!
  end if
!
  if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
!$OMP BARRIER
!   upgrade grid association if needed 
    do rs=thr_limits(3,tpi),thr_limits(4,tpi)
      call updateresgp(rs)
    end do
  end if
!
  call get_ensv(boxovol,tpi)
  call drift_removal(frac1,frac2,tpi)
!$OMP SINGLE
  call prt_curens(istep,boxovol,frac1,frac2)
!$OMP END SINGLE
!
  deallocate(mshfs)
!
  end
!
#endif
!
!-----------------------------------------------------------------------
!
! right now this function initializes internal coordinate gradients itself
! therefore care has to be taken in making assumptions about the current
! values of dc_di%f and dc_di%im
!
subroutine cart2int_f(cycle_frz,tpi)
!
  use iounit
  use forces
  use molecule
  use polypep
  use atoms
  use fyoc
  use math
  use zmatrix
  use sequen
  use aminos
  use mcsums
  use system
  use units
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
  logical, INTENT(IN):: cycle_frz
!
  integer imol,i,j,k,aone,atwo
  RTYPE pos2(3),pos3(3),bv(22),glcfv
  logical atrue
#ifdef ENABLE_THREADS
  integer ixx,sta,sto
!
  if (tpi.gt.0) then
    sta = thr_limits(35,tpi)
    sto = thr_limits(36,tpi)
  else
    sta = 1
    sto = nmol
  end if
#endif
!
  glcfv = 0.
  atrue = .true.
  aone = 1
  atwo = 2
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
  ens%insR(1) = 0.0
!$OMP END SINGLE
!
  do imol=sta,sto
    if (tpi.gt.0) then
      if (thr_mlgix(imol).gt.0) cycle
    end if
#else
  ens%insR(1) = 0.0
!
  do imol=1,nmol
#endif
!
!   cycle out if molecule entirely frozen
    if ((molfrzidx(imol).eq.4).AND.(cycle_frz.EQV..true.)) cycle
!
    dc_di(imol)%im(:) = 0.0
!       
!   WARNING!!!!! treatment of rotation missing for diatomic molecules or other linear molecules
    if ((atmol(imol,2)-atmol(imol,1)).le.1) then
      dc_di(imol)%im(1:3) = molmass(moltypid(imol))
      dc_di(imol)%f(1) = sum(cart_f(atmol(imol,1):atmol(imol,2),1))
      dc_di(imol)%f(2) = sum(cart_f(atmol(imol,1):atmol(imol,2),2))
      dc_di(imol)%f(3) = sum(cart_f(atmol(imol,1):atmol(imol,2),3))
      cart_v(atmol(imol,1):atmol(imol,2),1) = dc_di(imol)%v(1)
      cart_v(atmol(imol,1):atmol(imol,2),2) = dc_di(imol)%v(2)
      cart_v(atmol(imol,1):atmol(imol,2),3) = dc_di(imol)%v(3)
      cycle
    end if
!
 66 format(g14.7,1x,g14.7,10(i8,1x))
    cart_a(atmol(imol,1):atmol(imol,2),:) = cart_v(atmol(imol,1):atmol(imol,2),:)
!
    if (allocated(dc_di(imol)%valrecurs).EQV..true.) dc_di(imol)%valrecurs(:,:) = 0.0
    do i=1,dc_di(imol)%maxntor
      j = dc_di(imol)%recurs(i,1)
      pos2(1) = x(iz(1,j))
      pos2(2) = y(iz(1,j))
      pos2(3) = z(iz(1,j))
      pos3(1) = x(iz(2,j))
      pos3(2) = y(iz(2,j))
      pos3(3) = z(iz(2,j))
      do k=1,izrot(j)%treevs(3) ! transfer for branch merger
        dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(2)) = dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(2)) + &
 &                                                    dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(5+k))
        dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(5+k)) = 0.0
      end do
      call Vrecurse_dtor1(j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)),atrue,atrue)
      j = dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
    end do
!   rigid-body rotation degrees of freedom by diff-set
    if (allocated(dc_di(imol)%valrecurs).EQV..true.) then
      bv(1:17) = sum(dc_di(imol)%valrecurs(1:17,:),dim=2)
      dc_di(imol)%valrecurs(1:17,:) = 0.0
    else
      bv(:) = 0.0
    end if
    call Vrecurse_drot1(atmol(imol,1),imol,bv,atrue)
    do k=1,izrot(atmol(imol,1))%treevs(3) ! transfer RB result to all outgoing branches
      dc_di(imol)%valrecurs(7:22,izrot(atmol(imol,1))%treevs(5+k)) = bv(7:22)
    end do
!
    do i=dc_di(imol)%maxntor,1,-1
      j = dc_di(imol)%recurs(i,1)
      pos2(1) = x(iz(1,j))
      pos2(2) = y(iz(1,j))
      pos2(3) = z(iz(1,j))
      pos3(1) = x(iz(2,j))
      pos3(2) = y(iz(2,j))
      pos3(3) = z(iz(2,j))
      call Vrecurse_dtor2(j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)),atrue,atrue)
!      if (izrot(j)%treevs(2).eq.5) write(*,*) j,dc_di(imol)%valrecurs(17,izrot(j)%treevs(2))
      do k=1,izrot(j)%treevs(3) ! transfer to all outgoing branches
!         if (izrot(j)%treevs(5+k).eq.5) write(*,*) dc_di(imol)%valrecurs(17,izrot(j)%treevs(2))
        dc_di(imol)%valrecurs(7:22,izrot(j)%treevs(5+k)) = dc_di(imol)%valrecurs(7:22,izrot(j)%treevs(2))
      end do
    end do
!  do i=1,size(dc_di(imol)%valrecurs(1,:))
!    write(*,44) i,dc_di(imol)%valrecurs(17:18,i)
!  end do
 44 format(i5,3(g13.6,1x))
!
    cart_a(atmol(imol,1):atmol(imol,2),:)=(cart_v(atmol(imol,1):atmol(imol,2),:)-cart_a(atmol(imol,1):atmol(imol,2),:))/dyn_dt
    do i=atmol(imol,1),atmol(imol,2)
      if (mass(i).gt.0.0) then
        glcfv = glcfv + sum((mass(i)*cart_a(i,:)/u_dyn_fconv - cart_f(i,:))**2)/mass(i)
      end if
    end do
!
!    write(*,*) ens%insK,ens%insK+dyn_dt*dot_product(dc_di(imol)%btq(:,7),dc_di(imol)%v(:))/u_dyn_fconv_r
!    write(*,*) ens%insK,ens%insK+dyn_dt*dot_product(dc_di(imol)%btq(:,8),dc_di(imol)%v(:))/u_dyn_fconv_r
!
!
!!   first outward recursion to populate cart_a
!    call Vrecurse_drot2(atwo,atmol(imol,1),imol,bv)
!    do k=1,izrot(atmol(imol,1))%treevs(3) ! transfer RB result to all outgoing branches
!      dc_di(imol)%valrecurs(:,izrot(atmol(imol,1))%treevs(5+k)) = bv(:)
!    end do
!    do i=dc_di(imol)%maxntor,1,-1
!      j = dc_di(imol)%recurs(i,1)
!      pos2(1) = x(iz(1,j))
!      pos2(2) = y(iz(1,j))
!      pos2(3) = z(iz(1,j))
!      pos3(1) = x(iz(2,j))
!      pos3(2) = y(iz(2,j))
!      pos3(3) = z(iz(2,j))
!      if (izrot(j)%treevs(5).eq.0) then
!        call Vrecurse_dtor2b(aone,j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)))
!      else
!        call Vrecurse_dtor2b(atwo,j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)))
!      end if
!      do k=1,izrot(j)%treevs(3) ! transfer to all outgoing branches
!        dc_di(imol)%valrecurs(:,izrot(j)%treevs(5+k)) = dc_di(imol)%valrecurs(:,izrot(j)%treevs(2))
!      end do
!    end do
!!   final inward recursion to populate dc_di(imol)%btq
!    do i=1,dc_di(imol)%maxntor
!      j = dc_di(imol)%recurs(i,1)
!      pos2(1) = x(iz(1,j))
!      pos2(2) = y(iz(1,j))
!      pos2(3) = z(iz(1,j))
!      pos3(1) = x(iz(2,j))
!      pos3(2) = y(iz(2,j))
!      pos3(3) = z(iz(2,j))
!      do k=1,izrot(j)%treevs(3) ! transfer for branch merger
!        dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)) = dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)) + &
! &                                                    dc_di(imol)%valrecurs(:,izrot(j)%treevs(5+k))
!        dc_di(imol)%valrecurs(:,izrot(j)%treevs(5+k)) = 0.0
!      end do
!      if (izrot(j)%treevs(5).eq.0) then
!        call Vrecurse_dtor3(aone,j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)))
!      else
!        call Vrecurse_dtor3(atwo,j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)))
!      end if
!    end do
!!   rigid-body rotation degrees of freedom by diff-set
!    if (allocated(dc_di(imol)%valrecurs).EQV..true.) then
!      bv(:) = sum(dc_di(imol)%valrecurs(:,:),dim=2)
!      dc_di(imol)%valrecurs(:,:) = 0.0
!    else
!      bv(:) = 0.0
!    end if
!    call Vrecurse_drot3(atwo,atmol(imol,1),imol,bv)
!    write(*,*) ens%insK2,ens%insK2+dyn_dt*dot_product(dc_di(imol)%btq(:,1),dc_di(imol)%v(:))/u_dyn_fconv_r
!
!    write(*,*) 0.5*((dc_di(imol)%v(9)/RADIAN)**2)*dc_di(imol)%im(9)/u_dyn_fconv
!    write(*,*) ens%insK,' to ',ens%insK + dyn_dt*dc_di(imol)%v(10)*dc_di(imol)%btq(10)/u_dyn_fconv_r
!    write(*,*) dc_di(imol)%im(9),' to2 ',dc_di(imol)%im(9) + 2.0*dyn_dt*dc_di(imol)%v(10)*dc_di(imol)%btq(10)/RADIAN
!
!    do i=atmol(imol,2),atmol(imol,2)
!      write(*,*) (x(i)-xref(i))/dyn_dt,(y(i)-yref(i))/dyn_dt,(z(i)-zref(i))/dyn_dt
!      write(*,*) cart_v(i,:)
!      write(*,*) (cart_v(i,:)-cart_ldp(i,:,2))/dyn_dt
!      write(*,*) 'AC',cart_a(i,:)
!    end do
!    dc_di(imol)%f(:) =  dc_di(imol)%f(:) - dc_di(imol)%btq(:,7)/u_dyn_fconv
  end do
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(19,tpi) = glcfv
!$OMP BARRIER
!$OMP SINGLE
    ens%insR(1) = ens%insR(1) + sum(thr_rutil(19,1:thrdat%maxn))
!$OMP END SINGLE NOWAIT       
  else
    ens%insR(1) = ens%insR(1) + glcfv
  end if
#else
  ens%insR(1) = ens%insR(1) + glcfv
#endif
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    do ixx=1,nmlgs
      j = mlg_limits(ixx,5,tpi)
      call cart2int_threads(ixx,j,cycle_frz,tpi,aone)
    end do
  end if
#endif
!
end
!
!------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine cart2int_threads(idx,imol,cycle_frz,tpi,mode)
!
  use iounit
  use forces
  use molecule
  use polypep
  use atoms
  use fyoc
  use math
  use zmatrix
  use sequen
  use aminos
  use mcsums
  use system
  use units
  use threads
  use cutoffs
!
  implicit none
!
  integer, INTENT(IN):: tpi,idx,imol,mode
  logical, INTENT(IN):: cycle_frz
!
  integer i,ii,j,k,aone,atwo,maxbr,cntbr(maxval(mlg_limits(idx,8,:))),OMP_GET_NUM_THREADS,bnds(2)
  RTYPE pos2(3),pos3(3),bv(22),idt
  RTYPE bvt(22,mlg_limits(idx,10,tpi)-mlg_limits(idx,9,tpi)+1)
  RTYPE bvb2(22,maxval(mlg_limits(idx,8,:)))
  logical touchbr(maxval(mlg_limits(idx,8,:)))
  logical afalse,atrue
!
  aone = 1
  atwo = 2
  afalse = .false.
  atrue = .true.
  idt = 1.0/dyn_dt
  if ((mode.ne.1).AND.(mode.ne.2)) then
    write(ilog,*) 'Fatal. Called cart2int_threads(...) with unsupported mode (got ',mode,'). This is a bug.'
  end if
!
! cycle out if molecule entirely frozen
  if ((molfrzidx(imol).eq.4).AND.(cycle_frz.EQV..true.)) return
!
! WARNING!!!!! treatment of rotation missing for diatomic molecules or other linear molecules
  if ((atmol(imol,2)-atmol(imol,1)).le.1) then
!$OMP SINGLE
    if (mode.eq.1) then
      dc_di(imol)%im(1:3) = molmass(moltypid(imol))
      dc_di(imol)%f(1) = sum(cart_f(atmol(imol,1):atmol(imol,2),1))
      dc_di(imol)%f(2) = sum(cart_f(atmol(imol,1):atmol(imol,2),2))
      dc_di(imol)%f(3) = sum(cart_f(atmol(imol,1):atmol(imol,2),3))
      cart_v(atmol(imol,1):atmol(imol,2),1) = dc_di(imol)%v(1)
      cart_v(atmol(imol,1):atmol(imol,2),2) = dc_di(imol)%v(2)
      cart_v(atmol(imol,1):atmol(imol,2),3) = dc_di(imol)%v(3)
    else
      dc_di(imol)%olddat(1:3,2) = molmass(moltypid(imol))
    end if
!$OMP END SINGLE
    return
  end if
!
  maxbr =  maxval(mlg_limits(idx,8,:)) ! highest branch indicator
  touchbr(:) = .false.
!$OMP SINGLE
  vhlper(1:22,1) = 0.
  if (mode.eq.1) then
    dc_di(imol)%im(:) = 0.0
  else if (mode.eq.2) then
    dc_di(imol)%olddat(:,2) = 0.0
  end if
!$OMP END SINGLE NOWAIT
!
  if (mode.eq.1) then
    cart_a(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),1) = cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),1)
    cart_a(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),2) = cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),2)
    cart_a(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),3) = cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),3)
    cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),1) = dc_di(imol)%v(1)
    cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),2) = dc_di(imol)%v(2)
    cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),3) = dc_di(imol)%v(3)
  end if
!
  if (allocated(dc_di(imol)%valrecurs).EQV..true.) then
    dc_di(imol)%valrecurs(:,mlg_limits(idx,7,tpi):mlg_limits(idx,8,tpi)) = 0.0
  end if
  bv(:) = 0.0
  cntbr(:) = 0
  bvt(:,:) = 0.0
  thr_cart2f_hlp(1:17,:,tpi) = 0.0
!$OMP BARRIER
!
  if (mode.eq.1) then
    do i=mlg_limits(idx,9,tpi),mlg_limits(idx,10,tpi) ! 1,dc_di(imol)%maxntor
      ii = i - mlg_limits(idx,9,tpi) + 1
      j = dc_di(imol)%recurs(i,1)
      pos2(1) = x(iz(1,j))
      pos2(2) = y(iz(1,j))
      pos2(3) = z(iz(1,j))
      pos3(1) = x(iz(2,j))
      pos3(2) = y(iz(2,j))
      pos3(3) = z(iz(2,j))
      if (izrot(j)%alsz2.le.0) then ! tips have no data dependency and can be incremented (finalized) right away
        call Vrecurse_dtor1(j,imol,pos2(:),pos3(:),bvt(:,ii),atrue,atrue)
      else
        call Vrecurse_dtor1(j,imol,pos2(:),pos3(:),bvt(:,ii),afalse,atrue)
      end if
      thr_cart2f_hlp(1:17,izrot(j)%treevs(2),tpi) = thr_cart2f_hlp(1:17,izrot(j)%treevs(2),tpi) + bvt(1:17,ii)
      cntbr(izrot(j)%treevs(2)) = cntbr(izrot(j)%treevs(2)) + 1
    end do
  else if (mode.eq.2) then
    do i=mlg_limits(idx,9,tpi),mlg_limits(idx,10,tpi) ! 1,dc_di(imol)%maxntor
      ii = i - mlg_limits(idx,9,tpi) + 1
      j = dc_di(imol)%recurs(i,1)
      pos2(1) = x(iz(1,j))
      pos2(2) = y(iz(1,j))
      pos2(3) = z(iz(1,j))
      pos3(1) = x(iz(2,j))
      pos3(2) = y(iz(2,j))
      pos3(3) = z(iz(2,j))
      if (izrot(j)%alsz2.le.0) then ! tips have no data dependency and can be incremented (finalized) right away
        call Vrecurse_dtor0(j,imol,pos2(:),pos3(:),bvt(:,ii),atrue,atrue)
      else
        call Vrecurse_dtor0(j,imol,pos2(:),pos3(:),bvt(:,ii),afalse,atrue)
      end if
      thr_cart2f_hlp(1:17,izrot(j)%treevs(2),tpi) = thr_cart2f_hlp(1:17,izrot(j)%treevs(2),tpi) + bvt(1:17,ii)
      cntbr(izrot(j)%treevs(2)) = cntbr(izrot(j)%treevs(2)) + 1
    end do
  end if
!$OMP BARRIER
! this collects the final values for each branch (to be transferred to parent branch)
  call threads_bounds(1,maxbr,tpi,thrdat%maxn,bnds(1:2))
  do i=bnds(1),bnds(2)
    dc_di(imol)%valrecurs(1:17,i) = sum(thr_cart2f_hlp(1:17,i,1:thrdat%maxn),dim=2)
  end do
!$OMP BARRIER
  thr_rutil(1:17,tpi) = 0.0
  do i=mlg_limits(idx,7,tpi),mlg_limits(idx,8,tpi)
    thr_rutil(1:17,tpi) = thr_rutil(1:17,tpi) + dc_di(imol)%valrecurs(1:17,i)
    dc_di(imol)%valrecurs(:,i) = 0.0
  end do
! next we make a thread-local copy of valrecurs and transfer endpoint values from other threads
  do i=1,maxbr
    bvb2(1:17,i) = 0.0
  end do
  do j=1,tpi-1
    do i=1,maxbr
      bvb2(1:17,i) = bvb2(1:17,i) + thr_cart2f_hlp(1:17,i,j)
    end do
  end do
!
! now we can populate %im (or %olddat(:,2)) and %f: the work for merging branches is done by all threads (wasteful)
  if (mode.eq.1) then
    do i=1,mlg_limits(idx,10,tpi)
      ii = i - mlg_limits(idx,9,tpi) + 1
      j = dc_di(imol)%recurs(i,1)
      do k=1,izrot(j)%treevs(3) ! transfer for branch merger
        bvb2(1:17,izrot(j)%treevs(2)) = bvb2(1:17,izrot(j)%treevs(2)) + bvb2(1:17,izrot(j)%treevs(5+k))
      end do
      if (i.ge.mlg_limits(idx,9,tpi)) then
        if (izrot(j)%alsz2.le.0) then
          bvb2(1:17,izrot(j)%treevs(2)) = bvb2(1:17,izrot(j)%treevs(2)) + bvt(1:17,ii)
          cycle
        end if
        pos2(1) = x(iz(1,j))
        pos2(2) = y(iz(1,j))
        pos2(3) = z(iz(1,j))
        pos3(1) = x(iz(2,j))
        pos3(2) = y(iz(2,j))
        pos3(3) = z(iz(2,j))
        bvb2(1:17,izrot(j)%treevs(2)) = bvb2(1:17,izrot(j)%treevs(2)) + bvt(1:17,ii)
        call Vrecurse_dtor1(j,imol,pos2(:),pos3(:),bvb2(:,izrot(j)%treevs(2)),atrue,afalse)
      end if
    end do
  else if (mode.eq.2) then
    do i=1,mlg_limits(idx,10,tpi)
      ii = i - mlg_limits(idx,9,tpi) + 1
      j = dc_di(imol)%recurs(i,1)
      do k=1,izrot(j)%treevs(3) ! transfer for branch merger
        bvb2(1:17,izrot(j)%treevs(2)) = bvb2(1:17,izrot(j)%treevs(2)) + bvb2(1:17,izrot(j)%treevs(5+k))
      end do
      if (i.ge.mlg_limits(idx,9,tpi)) then
        if (izrot(j)%alsz2.le.0) then
          bvb2(1:17,izrot(j)%treevs(2)) = bvb2(1:17,izrot(j)%treevs(2)) + bvt(1:17,ii)
          cycle
        end if
        pos2(1) = x(iz(1,j))
        pos2(2) = y(iz(1,j))
        pos2(3) = z(iz(1,j))
        pos3(1) = x(iz(2,j))
        pos3(2) = y(iz(2,j))
        pos3(3) = z(iz(2,j))
        bvb2(1:17,izrot(j)%treevs(2)) = bvb2(1:17,izrot(j)%treevs(2)) + bvt(1:17,ii)
        call Vrecurse_dtor0(j,imol,pos2(:),pos3(:),bvb2(:,izrot(j)%treevs(2)),atrue,afalse)
      end if
    end do
  end if
  cntbr(:) = 0
  bvt(:,:) = 0.0
!$OMP BARRIER
!$OMP SINGLE
  vhlper(1:17,1) = sum(thr_rutil(1:17,1:thrdat%maxn),dim=2)
! rigid body motion is handled by a single thread
  if (mode.eq.1) then
    call Vrecurse_drot1(atmol(imol,1),imol,vhlper(1:22,1),afalse)
  else if (mode.eq.2) then
    call Vrecurse_drot0(atmol(imol,1),imol,vhlper(1:22,1))
  end if
!$OMP END SINGLE
!
  if (mode.eq.1) then
    thr_cart2f_hlp(17:22,:,tpi) = 0.
    do i=mlg_limits(idx,10,tpi),1,-1 ! dc_di(imol)%maxntor,1,-1
      ii = i - mlg_limits(idx,9,tpi) + 1
      j = dc_di(imol)%recurs(i,1)
      if (i.ge.mlg_limits(idx,9,tpi)) then
        pos2(1) = x(iz(1,j))
        pos2(2) = y(iz(1,j))
        pos2(3) = z(iz(1,j))
        pos3(1) = x(iz(2,j))
        pos3(2) = y(iz(2,j))
        pos3(3) = z(iz(2,j))
        if (cntbr(izrot(j)%treevs(2)).eq.0) then
          if (izrot(j)%treevs(1).le.0) then
!            write(*,*) 'root branch ',izrot(j)%treevs(2)
            thr_cart2f_hlp(17:22,izrot(j)%treevs(2),tpi) = vhlper(17:22,1)
          end if
        end if
        call Vrecurse_dtor2(j,imol,pos2(:),pos3(:),bvt(:,ii),afalse,atrue)
!        if (izrot(j)%treevs(2).eq.5) write(*,*) j,tpi,bvt(17,ii)
        thr_cart2f_hlp(17:22,izrot(j)%treevs(2),tpi) = thr_cart2f_hlp(17:22,izrot(j)%treevs(2),tpi) + bvt(17:22,ii)
        cntbr(izrot(j)%treevs(2)) = cntbr(izrot(j)%treevs(2)) + 1
      end if
      do k=1,izrot(j)%treevs(3) ! assign for all outgoing branches
!          write(*,*) 'setting ',izrot(j)%treevs(5+k),' from ',izrot(j)%treevs(2),thr_cart2f_hlp(17,izrot(j)%treevs(2),tpi),tpi
          thr_cart2f_hlp(17:22,izrot(j)%treevs(5+k),tpi) = thr_cart2f_hlp(17:22,izrot(j)%treevs(2),tpi)
      end do
    end do
!
!$OMP BARRIER
!
!   again, we make a thread-local copy of valrecurs and transfer endpoint values from other threads
    do i=1,maxbr
      bvb2(17:22,i) = 0.0
    end do
    do j=tpi,OMP_GET_NUM_THREADS()
      do i=1,maxbr
        bvb2(17:22,i) = bvb2(17:22,i) + thr_cart2f_hlp(17:22,i,j)
      end do
    end do
    if (tpi.eq.1) then
!    do i=1,maxbr
!      write(*,44) i,bvb2(17:18,i)
!    end do
    end if
    do i=mlg_limits(idx,9,tpi),mlg_limits(idx,10,tpi)
      ii = i - mlg_limits(idx,9,tpi) + 1
      j = dc_di(imol)%recurs(i,1)
      pos2(1) = x(iz(1,j))
      pos2(2) = y(iz(1,j))
      pos2(3) = z(iz(1,j))
      pos3(1) = x(iz(2,j))
      pos3(2) = y(iz(2,j))
      pos3(3) = z(iz(2,j))
      call Vrecurse_dtor2(j,imol,pos2(:),pos3(:),bvb2(:,izrot(j)%treevs(2)),atrue,afalse)
      bvb2(17:22,izrot(j)%treevs(2)) = bvb2(17:22,izrot(j)%treevs(2)) - bvt(17:22,ii)
    end do
   44 format(i5,3(g13.6,1x))
!$OMP BARRIER
!
    cart_a(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),:)=idt*(cart_v(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),:) - &
                                                               cart_a(mlg_limits(idx,1,tpi):mlg_limits(idx,2,tpi),:))
    thr_rutil(18,tpi) = 0.
    do i=mlg_limits(idx,1,tpi),mlg_limits(idx,2,tpi)
      if (mass(i).gt.0.0) then
        thr_rutil(18,tpi) = thr_rutil(18,tpi) + sum((mass(i)*cart_a(i,:)/u_dyn_fconv - cart_f(i,:))**2)/mass(i)
      end if
    end do
!$OMP BARRIER
!$OMP SINGLE
    ens%insR(1) = ens%insR(1) + sum(thr_rutil(18,1:thrdat%maxn))
!$OMP END SINGLE
  end if
!
end
!
#endif
!
!-----------------------------------------------------------------------
!
subroutine cart2int_I(tpi)
!
  use iounit
  use forces
  use molecule
  use polypep
  use atoms
  use fyoc
  use math
  use zmatrix
  use sequen
  use aminos
  use mcsums
  use system
  use units
  use movesets
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer imol,i,j,k,aone,atwo
  RTYPE pos2(3),pos3(3),bv(22)
  logical atrue
#ifdef ENABLE_THREADS
  integer ixx,sta,sto
!
  if (tpi.gt.0) then
    sta = thr_limits(35,tpi)
    sto = thr_limits(36,tpi)
  else
    sta = 1
    sto = nmol
  end if
#endif
!
  atrue = .true.
  aone = 1
  atwo = 2
!
#ifdef ENABLE_THREADS
  do imol=sta,sto
    if (tpi.gt.0) then
      if (thr_mlgix(imol).gt.0) cycle
    end if
#else
  do imol=1,nmol
#endif
!
!   cycle out if molecule entirely frozen
    if ((molfrzidx(imol).eq.4).AND.(skip_frz.EQV..true.)) cycle
!
!   WARNING!!!!! treatment of rotation missing for diatomic molecules or other linear molecules
    if ((atmol(imol,2)-atmol(imol,1)).le.1) then
      dc_di(imol)%olddat(:,2) = molmass(moltypid(imol))
      cycle
    end if
!
    dc_di(imol)%olddat(:,2) = 0.0
!
    if (allocated(dc_di(imol)%valrecurs).EQV..true.) dc_di(imol)%valrecurs(:,:) = 0.0
    do i=1,dc_di(imol)%maxntor
      j = dc_di(imol)%recurs(i,1)
      pos2(1) = x(iz(1,j))
      pos2(2) = y(iz(1,j))
      pos2(3) = z(iz(1,j))
      pos3(1) = x(iz(2,j))
      pos3(2) = y(iz(2,j))
      pos3(3) = z(iz(2,j))
     do k=1,izrot(j)%treevs(3) ! transfer for branch merger
        dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(2)) = dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(2)) + &
 &                                                    dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(5+k))
        dc_di(imol)%valrecurs(1:17,izrot(j)%treevs(5+k)) = 0.0
      end do
      call Vrecurse_dtor0(j,imol,pos2(:),pos3(:),dc_di(imol)%valrecurs(:,izrot(j)%treevs(2)),atrue,atrue)
    end do
!   rigid-body rotation degrees of freedom by diff-set
    if (allocated(dc_di(imol)%valrecurs).EQV..true.) then
      bv(1:17) = sum(dc_di(imol)%valrecurs(1:17,:),dim=2)
      dc_di(imol)%valrecurs(1:17,:) = 0.0
    else
      bv(:) = 0.0
    end if
    call Vrecurse_drot0(atmol(imol,1),imol,bv)
  end do
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    do ixx=1,nmlgs
      j = mlg_limits(ixx,5,tpi)
      call cart2int_threads(ixx,j,skip_frz,tpi,atwo)
    end do
  end if
#endif
!
end
!
!----------------------------------------------------------------------------------------------------------------
!
! to give a single example of a recursion formula, consider the following:
! when using difference sets, we can generate Cartesian velocities in effectively a single O(N) loop:
! v_i = v_t + r_i x Sum_m( ax_m * omega_m ) - Sum_m [( p_m x ax_m )*omega_m]
! r_i is position vector, ax_m are all relevant rotation axes, p_m reference points, omega_m angular velocities, and v_t
! is the translational velocity of the molecule;
! note that v_t is set before and that base atoms are already incremented for RB rotation via Vrecurse_drot(...)
!
subroutine Vrecurse_dtor2(ati,imol,pos2,pos3,branchv,do_up,do_inc)
!
  use iounit
  use forces
  use atoms
  use zmatrix
  use math
  use units
  use system
!
  implicit none
!
  integer, INTENT(IN):: ati,imol
  logical, INTENT(IN):: do_up,do_inc
  RTYPE, INTENT(IN):: pos2(3),pos3(3)
!
  integer ttc,k,ii,ik,ilim
! integer i,ati2,kk,j
  RTYPE con1(3),con3(3),cp1(3),c2,c1,branchv(22),con0(3),veff
! RTYPE pos1(3),cp2(3),bv2(12)
!
  ilim = 0
  if (izrot(ati)%alsz2.eq.0) then
    ilim = izrot(ati)%alsz
    if (allocated(izrot(ati)%rotis).EQV..false.) then
      write(ilog,*) 'Fatal. Called Vrecurse_dtor2(...) for a line, for which the required atom set is not &
 &allocated. This is a bug.'
      call fexit()
    end if
  else 
    ilim = izrot(ati)%alsz2
  end if
!
  ttc = dc_di(imol)%recurs(izrot(ati)%treevs(4),3)
!
!! this block can be used to verify that the recursion scheme is working properly for dK.E./dphi_l based on
!! the active form (phi_l altering I_m further toward base)
!  hlp = 0.
!  do k=1,izrot(ati)%alsz
!    ii = izrot(ati)%rotis(k,1)
!    ik = izrot(ati)%rotis(k,2)
!    do i=ii,ik
!!     net accel for torsions not including self
!      cp2(1) =  (branchv(10)-branchv(7))*x(i) - branchv(11)*y(i) - branchv(13)*z(i) + branchv(14)
!      cp2(2) = -branchv(11)*x(i) + (branchv(10)-branchv(8))*y(i) - branchv(12)*z(i) + branchv(15)
!      cp2(3) = -branchv(13)*x(i) - branchv(12)*y(i) + (branchv(10)-branchv(9))*z(i) + branchv(16)
!      pos1(1) = x(i) - pos2(1)
!      pos1(2) = y(i) - pos2(2)
!      pos1(3) = z(i) - pos2(3)
!      call crossprod3(con1(:),pos1(:),cp1(:))
!      hlp = hlp + mass(i)*dot_product(cp1(:),cp2(:))
!    end do
!  end do
!
  if ((ttc.gt.0).AND.(do_inc.EQV..true.)) then
    con0(:) = pos3(:)-pos2(:)
    if (dc_di(imol)%recurs(izrot(ati)%treevs(4),2).eq.1) con0(:) = -con0(:)
    c2 = 1.0/(dot_product(con0(:),con0(:)))
    c1 = sqrt(c2)
    con1(:) = -con0(:)*c1
    veff = dc_di(imol)%v(ttc)/RADIAN
    con3(:) = veff*con1(:)
!    branchv(7:9) = branchv(7:9) + con3(:)*con3(:)
!    branchv(11) = branchv(11) + con3(1)*con3(2)
!    branchv(12) = branchv(12) + con3(2)*con3(3)
!    branchv(13) = branchv(13) + con3(3)*con3(1)
    call crossprod3(pos2(:),con3(:),cp1(:))
!    call crossprod3(cp1(:),con3(:),cp2(:))
!    branchv(14:16) = branchv(14:16) + cp2(:)
!    branchv(10) = branchv(10) + veff*veff
    branchv(17:19) = branchv(17:19) + con3(:)
    branchv(20:22) = branchv(20:22) + cp1(:)
!!   set bias term
!   call crossprod3(con1(:),pos2(:),cp1(:))
!   dc_di(imol)%btq(ttc,7) = dot_product(con1(:),branchv(4:6)) - dot_product(cp1(:),branchv(1:3))
!!    dc_di(imol)%btq(ttc,8) = hlp
  else
    veff = 0.0 ! note that the above increments are all zero with veff = 0.0
  end if
!  if (ttc.gt.0) write(*,*) 'VER',hlp,izrot(ati)%alsz2,ati
!ens%insK,ens%insK+dyn_dt*hlp*veff/u_dyn_fconv

!  if (ttc.gt.0) write(*,*) 'REC',dc_di(imol)%btq(ttc,7)
!ens%insK,ens%insK+dyn_dt*dc_di(imol)%btq(ttc,7)*veff/u_dyn_fconv
!
  if (do_up.EQV..true.) then
    do k=1,ilim
      if (izrot(ati)%alsz2.eq.0) then
        ii = izrot(ati)%rotis(k,1)
        ik = izrot(ati)%rotis(k,2)
      else
        ii = izrot(ati)%diffis(k,1)
        ik = izrot(ati)%diffis(k,2)
      end if
!      cart_a(ii:ik,1) = cart_a(ii:ik,1) + (branchv(10)-branchv(7))*x(ii:ik) - branchv(11)*y(ii:ik) - branchv(13)*z(ii:ik) + &
!   &branchv(14)
!      cart_a(ii:ik,2) = cart_a(ii:ik,2) - branchv(11)*x(ii:ik) + (branchv(10)-branchv(8))*y(ii:ik) - branchv(12)*z(ii:ik) + &
!   &branchv(15)
!      cart_a(ii:ik,3) = cart_a(ii:ik,3) - branchv(13)*x(ii:ik) - branchv(12)*y(ii:ik) + (branchv(10)-branchv(9))*z(ii:ik) + &
!   &branchv(16)
      cart_v(ii:ik,1) = cart_v(ii:ik,1) + branchv(18)*z(ii:ik) - branchv(19)*y(ii:ik) + branchv(20)
      cart_v(ii:ik,2) = cart_v(ii:ik,2) + branchv(19)*x(ii:ik) - branchv(17)*z(ii:ik) + branchv(21)
      cart_v(ii:ik,3) = cart_v(ii:ik,3) + branchv(17)*y(ii:ik) - branchv(18)*x(ii:ik) + branchv(22)
    end do
  end if
!
!  if (ttc.gt.0) then
!!   this is first incremented by the entire set of acceleration terms for just this torsion (but all atoms toward tip)
!    branchv(1:6) = branchv(1:6) + dc_di(imol)%btq(ttc,1:6)
!  end if
!  do k=1,ilim
!    if (izrot(ati)%alsz2.eq.0) then
!      ii = izrot(ati)%rotis(k,1)
!      ik = izrot(ati)%rotis(k,2)
!    else
!      ii = izrot(ati)%diffis(k,1)
!      ik = izrot(ati)%diffis(k,2)
!    end if
!  ! and now we remove the net accelerations for diff set only
!    do i=ii,ik
!      pos1(1) = x(i)
!      pos1(2) = y(i)
!      pos1(3) = z(i)
!      branchv(1:3) = branchv(1:3) - mass(i)*cart_a(i,:)
!      call crossprod3(pos1,cart_a(i,:),cp1(:))
!      branchv(4:6) = branchv(4:6) - mass(i)*cp1(:)
!    end do
!  end do
!  if (izrot(ati)%treevs(3).gt.0) then
!!   we cannot transfer the entire set to every branch, instead use branch-specific terms
!!   currently, these need to be recomputed
!    bv2(7:12) = 0.0
!    do kk=1,izrot(ati)%treevs(3)
!      do j=izrot(ati)%treevs(4),1,-1
!        ii = dc_di(imol)%recurs(j,1)
!        if ((izrot(ii)%treevs(1).eq.ati).AND.(izrot(ii)%treevs(2).ne.izrot(ati)%treevs(2)).AND.&
! &(izrot(ii)%treevs(2).eq.izrot(ati)%treevs(5+kk))) then
!          ati2 = ii
!          bv2(1:6) = 0.0
!          do k=1,izrot(ati2)%alsz
!            ii = izrot(ati2)%rotis(k,1)
!            ik = izrot(ati2)%rotis(k,2)
!            do i=ii,ik
!!             net accel for torsions not including self
!              cp2(1) =  (branchv(10)-branchv(7))*x(i) - branchv(11)*y(i) - branchv(13)*z(i) + branchv(14)
!              cp2(2) = -branchv(11)*x(i) + (branchv(10)-branchv(8))*y(i) - branchv(12)*z(i) + branchv(15)
!              cp2(3) = -branchv(13)*x(i) - branchv(12)*y(i) + (branchv(10)-branchv(9))*z(i) + branchv(16)
!              cp1(1) = y(i)*cp2(3) - z(i)*cp2(2)
!              cp1(2) = z(i)*cp2(1) - x(i)*cp2(3)
!              cp1(3) = x(i)*cp2(2) - y(i)*cp2(1)
!              bv2(1:3) = bv2(1:3) + mass(i)*cp2(:)
!              bv2(4:6) = bv2(4:6) + mass(i)*cp1(:)
!            end do
!          end do
!          dc_di(imol)%valrecurs(1:6,izrot(ati)%treevs(5+kk)) = bv2(1:6)
!          bv2(7:12) = bv2(7:12) + bv2(1:6)
!          exit
!        end if
!      end do
!    end do
!    branchv(1:6) = branchv(1:6) - bv2(7:12)
!  end if
!
end
!
!------------------------------------------------------------------------------------------------------
!
! just inertia
!
subroutine Vrecurse_dtor0(ati,imol,pos2,pos3,branchv,do_up,do_inc)
!
  use iounit
  use forces
  use atoms
  use zmatrix
  use math
!
  implicit none
!
  integer, INTENT(IN):: imol,ati
  logical, INTENT(IN):: do_up,do_inc
  RTYPE, INTENT(IN):: pos2(3),pos3(3)
!
  integer ttc,i,k,ii,ik,ilim
  RTYPE con1(3),c2,c1,branchv(22),hlp2(3),con0(3),hlp3(3),dp1
!
  if (do_inc.EQV..true.) then
    ilim = 0
    if (izrot(ati)%alsz2.eq.0) then
      ilim = izrot(ati)%alsz
      if (allocated(izrot(ati)%rotis).EQV..false.) then
        write(ilog,*) 'Fatal. Called Vrecurse_dtor0(...) for a line, for which the required atom set is not &
 &allocated. This is a bug.'
        call fexit()
      end if
      branchv(:) = 0.0
    else
      ilim = izrot(ati)%alsz2
    end if
    do k=1,ilim
      if (izrot(ati)%alsz2.eq.0) then
        ii = izrot(ati)%rotis(k,1)
        ik = izrot(ati)%rotis(k,2)
      else
        ii = izrot(ati)%diffis(k,1)
        ik = izrot(ati)%diffis(k,2)
      end if
      do i=ii,ik
        hlp2(1) = x(i)
        hlp2(2) = y(i)
        hlp2(3) = z(i)
        branchv(17) = branchv(17) + mass(i)*dot_product(hlp2(:),hlp2(:))
      end do
      branchv(1) = branchv(1) + sum(mass(ii:ik)*x(ii:ik))
      branchv(2) = branchv(2) + sum(mass(ii:ik)*y(ii:ik))
      branchv(3) = branchv(3) + sum(mass(ii:ik)*z(ii:ik))
      branchv(4) = branchv(4) + sum(mass(ii:ik))
      branchv(5) = branchv(5) + sum(mass(ii:ik)*x(ii:ik)*x(ii:ik))
      branchv(6) = branchv(6) + sum(mass(ii:ik)*y(ii:ik)*y(ii:ik))
      branchv(7) = branchv(7) + sum(mass(ii:ik)*z(ii:ik)*z(ii:ik))
      branchv(8) = branchv(8) + sum(mass(ii:ik)*x(ii:ik)*y(ii:ik))
      branchv(9) = branchv(9) + sum(mass(ii:ik)*y(ii:ik)*z(ii:ik))
      branchv(10) = branchv(10) + sum(mass(ii:ik)*z(ii:ik)*x(ii:ik))
    end do
  end if
!
  ttc = dc_di(imol)%recurs(izrot(ati)%treevs(4),3)
  if ((ttc.gt.0).AND.(do_up.EQV..true.)) then
    
    con0(:) = pos3(:)-pos2(:)
    if (dc_di(imol)%recurs(izrot(ati)%treevs(4),2).eq.1) con0(:) = -con0(:)
    c2 = 1.0/(dot_product(con0(:),con0(:)))
    c1 = sqrt(c2)
    con1(:) = -con0(:)*c1
    dp1 = dot_product(pos2(:),con0(:))
    hlp2(:) = con0(:)*con0(:)
    hlp3(1) = con0(1)*con0(2)
    hlp3(2) = con0(2)*con0(3)
    hlp3(3) = con0(3)*con0(1)
!
    dc_di(imol)%olddat(ttc,2) = dc_di(imol)%olddat(ttc,2) + branchv(4)*dot_product(pos2(:),pos2(:)) - &
 &                        2.0*dot_product(pos2(:),branchv(1:3)) &
 &                        + branchv(17) - c2*(dp1*dp1*branchv(4) - 2.0*dp1*dot_product(con0(:),branchv(1:3)) + &
 &                        dot_product(hlp2(:),branchv(5:7)) + 2.0*dot_product(hlp3(:),branchv(8:10)))
  end if
!
end
!
!------------------------------------------------------------------------------------------------------
!
! this populates forces and inertia for torsional degrees of freedom during the first inward recursion
! function can be optimized much further (ii,ik loop)
!
subroutine Vrecurse_dtor1(ati,imol,pos2,pos3,branchv,do_up,do_inc)
!
  use iounit
  use forces
  use atoms
  use zmatrix
  use math
!
  implicit none
!
  logical, INTENT(IN):: do_up,do_inc
  integer, INTENT(IN):: imol,ati
  RTYPE, INTENT(IN):: pos2(3),pos3(3)
!
  integer ttc,i,k,ii,ik,ilim
  RTYPE con1(3),cp1(3),c2,c1,dp1,branchv(22),hlp2(3),con0(3),hlp3(3)
! RTYPE veff,con3(3),cp2(3)
!
  if (do_inc.EQV..true.) then
    ilim = 0
    if (izrot(ati)%alsz2.eq.0) then
      ilim = izrot(ati)%alsz
      if (allocated(izrot(ati)%rotis).EQV..false.) then
        write(ilog,*) 'Fatal. Called Vrecurse_dtor1(...) for a line, for which the required atom set is not &
   &allocated. This is a bug.'
        call fexit()
      end if
      branchv(:) = 0.0
    else
      ilim = izrot(ati)%alsz2
    end if
!
    do k=1,ilim
      if (izrot(ati)%alsz2.eq.0) then
        ii = izrot(ati)%rotis(k,1)
        ik = izrot(ati)%rotis(k,2)
      else
        ii = izrot(ati)%diffis(k,1)
        ik = izrot(ati)%diffis(k,2)
      end if
      do i=ii,ik
        hlp2(1) = x(i)
        hlp2(2) = y(i)
        hlp2(3) = z(i)
        hlp3(:) = cart_f(i,:)
        call crossprod3(hlp2(:),hlp3(:),cp1(:))
        branchv(11:13) = branchv(11:13) + cp1(:)
        branchv(17) = branchv(17) + mass(i)*dot_product(hlp2(:),hlp2(:))
      end do
      branchv(1) = branchv(1) + sum(mass(ii:ik)*x(ii:ik))
      branchv(2) = branchv(2) + sum(mass(ii:ik)*y(ii:ik))
      branchv(3) = branchv(3) + sum(mass(ii:ik)*z(ii:ik))
      branchv(4) = branchv(4) + sum(mass(ii:ik))
      branchv(5) = branchv(5) + sum(mass(ii:ik)*x(ii:ik)*x(ii:ik))
      branchv(6) = branchv(6) + sum(mass(ii:ik)*y(ii:ik)*y(ii:ik))
      branchv(7) = branchv(7) + sum(mass(ii:ik)*z(ii:ik)*z(ii:ik))
      branchv(8) = branchv(8) + sum(mass(ii:ik)*x(ii:ik)*y(ii:ik))
      branchv(9) = branchv(9) + sum(mass(ii:ik)*y(ii:ik)*z(ii:ik))
      branchv(10) = branchv(10) + sum(mass(ii:ik)*z(ii:ik)*x(ii:ik))
      branchv(14) = branchv(14) + sum(cart_f(ii:ik,1))
      branchv(15) = branchv(15) + sum(cart_f(ii:ik,2))
      branchv(16) = branchv(16) + sum(cart_f(ii:ik,3))
    end do
  end if
!
  ttc = dc_di(imol)%recurs(izrot(ati)%treevs(4),3)
  if ((ttc.gt.0).AND.(do_up.EQV..true.)) then
    con0(:) = pos3(:)-pos2(:)
    if (dc_di(imol)%recurs(izrot(ati)%treevs(4),2).eq.1) con0(:) = -con0(:)
    c2 = 1.0/(dot_product(con0(:),con0(:)))
    c1 = sqrt(c2)
    con1(:) = -con0(:)*c1
!    veff = dc_di(imol)%v(ttc)/RADIAN
!    con3(:) = con1(:)*veff
!    call crossprod3(con3(:),pos2(:),cp1(:))
!    call crossprod3(cp1(:),con3(:),cp2(:))
!    hlp2(:) = con3(:)*con3(:)
!    hlp3(1) = con3(1)*con3(2)
!    hlp3(2) = con3(2)*con3(3)
!    hlp3(3) = con3(3)*con3(1)
!    dc_di(imol)%btq(ttc,1) = veff*veff*branchv(1) - branchv(4)*cp2(1) - hlp2(1)*branchv(1)-hlp3(1)*branchv(2)-hlp3(3)*branchv(3)
!    dc_di(imol)%btq(ttc,2) = veff*veff*branchv(2) - branchv(4)*cp2(2) - hlp3(1)*branchv(1)-hlp2(2)*branchv(2)-hlp3(2)*branchv(3)
!    dc_di(imol)%btq(ttc,3) = veff*veff*branchv(3) - branchv(4)*cp2(3) - hlp3(3)*branchv(1)-hlp3(2)*branchv(2)-hlp2(3)*branchv(3)
!    dc_di(imol)%btq(ttc,4) = hlp3(1)*branchv(10) + (hlp2(2)-hlp2(3))*branchv(9) + hlp3(2)*(branchv(7)-branchv(6)) - &
! &                           hlp3(3)*branchv(8) - branchv(2)*cp2(3) + branchv(3)*cp2(2)
!    dc_di(imol)%btq(ttc,5) = hlp3(3)*(branchv(5)-branchv(7)) + hlp3(2)*branchv(8) + (hlp2(3)-hlp2(1))*branchv(10) - &
! &                           hlp3(1)*branchv(9) - branchv(3)*cp2(1) + branchv(1)*cp2(3)
!    dc_di(imol)%btq(ttc,6) = (hlp2(1)-hlp2(2))*branchv(8) + hlp3(1)*(branchv(6)-branchv(5)) + hlp3(3)*branchv(9) - &
! &                           hlp3(2)*branchv(10) - branchv(1)*cp2(2) + branchv(2)*cp2(1)
    call crossprod3(pos3(:),branchv(14:16),cp1(:))
    dc_di(imol)%f(ttc) = dc_di(imol)%f(ttc) + dot_product(con1(:),branchv(11:13)) - dot_product(con1(:),cp1(:))
    dp1 = dot_product(pos2(:),con0(:))
    hlp2(:) = con0(:)*con0(:)
    hlp3(1) = con0(1)*con0(2)
    hlp3(2) = con0(2)*con0(3)
    hlp3(3) = con0(3)*con0(1)
!
    dc_di(imol)%im(ttc) = dc_di(imol)%im(ttc) + branchv(4)*dot_product(pos2(:),pos2(:)) - 2.0*dot_product(pos2(:),branchv(1:3)) &
 &                        + branchv(17) - c2*(dp1*dp1*branchv(4) - 2.0*dp1*dot_product(con0(:),branchv(1:3)) + &
 &                        dot_product(hlp2(:),branchv(5:7)) + 2.0*dot_product(hlp3(:),branchv(8:10)))
  end if
!
end
!
!--------------------------------------------------------------------------------------------------
!
subroutine Vrecurse_drot1(ati,imol,branchv,do_init)
!
  use iounit
  use forces
  use atoms
  use zmatrix
  use math
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: imol,ati
  logical, INTENT(IN):: do_init
!
  integer ttc,i,k,ii,ik,ilim,rrr
! integer j,kk,ati2
  RTYPE pos2(3,3),pos3(3),con1(3,3),cp1(3),c2,c1,branchv(22),hlp2(3),con3(3),con0(3),veff,hlp3(3),dp1
! RTYPE pos1(3),cp2(3),bv2(12)
!
  pos3(1:3) = comm(imol,1:3)
  do rrr=1,3
    pos2(1:3,rrr) = comm(imol,1:3) + rgpcsm(imol,rrr,1:3)
    con0(:) = pos3(:)-pos2(:,rrr)
    c2 = 1.0/(dot_product(con0(:),con0(:)))
    c1 = sqrt(c2)
    con1(:,rrr) = -con0(:)*c1
  end do
!
  ilim = 0
  if (izrot(ati)%alsz2.eq.0) then
    ilim = izrot(ati)%alsz
    if (allocated(izrot(ati)%rotis).EQV..false.) then
      write(ilog,*) 'Fatal. Called Vrecurse_drot1(...) for a line, for which the required atom set is not &
 &allocated. This is a bug.'
      call fexit()
    end if
    branchv(:) = 0.0
  else
    ilim = izrot(ati)%alsz2
  end if
!
  do k=1,ilim
    if (izrot(ati)%alsz2.eq.0) then
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
    else
      ii = izrot(ati)%diffis(k,1)
      ik = izrot(ati)%diffis(k,2)
    end if
    do i=ii,ik
      hlp2(1) = x(i)
      hlp2(2) = y(i)
      hlp2(3) = z(i)
      hlp3(:) = cart_f(i,:)
      call crossprod3(hlp2(:),hlp3(:),cp1(:))
      branchv(11:13) = branchv(11:13) + cp1(:)
      branchv(17) = branchv(17) + mass(i)*dot_product(hlp2(:),hlp2(:))
    end do
    branchv(1) = branchv(1) + sum(mass(ii:ik)*x(ii:ik))
    branchv(2) = branchv(2) + sum(mass(ii:ik)*y(ii:ik))
    branchv(3) = branchv(3) + sum(mass(ii:ik)*z(ii:ik))
    branchv(4) = branchv(4) + sum(mass(ii:ik))
    branchv(5) = branchv(5) + sum(mass(ii:ik)*x(ii:ik)*x(ii:ik))
    branchv(6) = branchv(6) + sum(mass(ii:ik)*y(ii:ik)*y(ii:ik))
    branchv(7) = branchv(7) + sum(mass(ii:ik)*z(ii:ik)*z(ii:ik))
    branchv(8) = branchv(8) + sum(mass(ii:ik)*x(ii:ik)*y(ii:ik))
    branchv(9) = branchv(9) + sum(mass(ii:ik)*y(ii:ik)*z(ii:ik))
    branchv(10) = branchv(10) + sum(mass(ii:ik)*z(ii:ik)*x(ii:ik))
    branchv(14) = branchv(14) + sum(cart_f(ii:ik,1))
    branchv(15) = branchv(15) + sum(cart_f(ii:ik,2))
    branchv(16) = branchv(16) + sum(cart_f(ii:ik,3))
  end do
!
!! the set of atoms is identical and complete for all three rigid rotations
!  do rrr=1,3
!    ttc = rrr + 3
!    veff = dc_di(imol)%v(ttc)/RADIAN
!    con3(:) = con1(:,rrr)*veff
!    call crossprod3(con3(:),pos2(:,rrr),cp1(:))
!    call crossprod3(cp1(:),con3(:),cp2(:))
!    hlp2(:) = con3(:)*con3(:)
!    hlp3(1) = con3(1)*con3(2)
!    hlp3(2) = con3(2)*con3(3)
!    hlp3(3) = con3(3)*con3(1)
!    dc_di(imol)%btq(ttc,1) = veff*veff*branchv(1) - branchv(4)*cp2(1) - hlp2(1)*branchv(1)-hlp3(1)*branchv(2)-hlp3(3)*branchv(3)
!    dc_di(imol)%btq(ttc,2) = veff*veff*branchv(2) - branchv(4)*cp2(2) - hlp3(1)*branchv(1)-hlp2(2)*branchv(2)-hlp3(2)*branchv(3)
!    dc_di(imol)%btq(ttc,3) = veff*veff*branchv(3) - branchv(4)*cp2(3) - hlp3(3)*branchv(1)-hlp3(2)*branchv(2)-hlp2(3)*branchv(3)
!    dc_di(imol)%btq(ttc,4) = hlp3(1)*branchv(10) + (hlp2(2)-hlp2(3))*branchv(9) + hlp3(2)*(branchv(7)-branchv(6)) - &
! &                           hlp3(3)*branchv(8) - branchv(2)*cp2(3) + branchv(3)*cp2(2)
!    dc_di(imol)%btq(ttc,5) = hlp3(3)*(branchv(5)-branchv(7)) + hlp3(2)*branchv(8) + (hlp2(3)-hlp2(1))*branchv(10) - &
! &                           hlp3(1)*branchv(9) - branchv(3)*cp2(1) + branchv(1)*cp2(3)
!    dc_di(imol)%btq(ttc,6) = (hlp2(1)-hlp2(2))*branchv(8) + hlp3(1)*(branchv(6)-branchv(5)) + hlp3(3)*branchv(9) - &
! &                           hlp3(2)*branchv(10) - branchv(1)*cp2(2) + branchv(2)*cp2(1)
!  end do
  call crossprod3(pos3(:),branchv(14:16),cp1(:))
  do rrr=1,3
    ttc = rrr + 3
    dc_di(imol)%f(ttc) = dc_di(imol)%f(ttc) + dot_product(con1(:,rrr),branchv(11:13)) - dot_product(con1(:,rrr),cp1(:))
    dp1 = dot_product(pos2(:,rrr),con1(:,rrr))
    hlp2(:) = con1(:,rrr)*con1(:,rrr)
    hlp3(1) = con1(1,rrr)*con1(2,rrr)
    hlp3(2) = con1(2,rrr)*con1(3,rrr)
    hlp3(3) = con1(3,rrr)*con1(1,rrr)
    dc_di(imol)%im(ttc) = dc_di(imol)%im(ttc) + branchv(4)*dot_product(pos2(:,rrr),pos2(:,rrr)) - &
 & 2.0*dot_product(pos2(:,rrr),branchv(1:3)) + branchv(17) - dp1*dp1*branchv(4) + &
 & 2.0*dp1*dot_product(con1(:,rrr),branchv(1:3)) - dot_product(hlp2(:),branchv(5:7)) - 2.0*dot_product(hlp3(:),branchv(8:10))
  end do
  dc_di(imol)%f(1:3) = branchv(14:16)
  dc_di(imol)%im(1:3) = molmass(moltypid(imol))
!
  branchv(:) = 0.0
! sum contributions from all RB rotations
!  branchv(1:6) = sum(dc_di(imol)%btq(4:6,1:6),dim=1)
!
! now accumulate the angular velocity-related terms for forward recursion starting with RB rot
  do rrr=1,3
    ttc = rrr + 3
    veff = dc_di(imol)%v(ttc)/RADIAN
    con3(:) = con1(:,rrr)*veff
!    branchv(7:9) = branchv(7:9) + con3(:)*con3(:)
!    branchv(11) = branchv(11) + con3(1)*con3(2)
!    branchv(12) = branchv(12) + con3(2)*con3(3)
!    branchv(13) = branchv(13) + con3(3)*con3(1)
    call crossprod3(pos2(:,rrr),con3(:),cp1(:))
!    call crossprod3(cp1(:),con3(:),cp2(:))
!    branchv(14:16) = branchv(14:16) + cp2(:)
!    branchv(10) = branchv(10) + veff*veff
    branchv(17:19) = branchv(17:19) + con3(:)
    branchv(20:22) = branchv(20:22) + cp1(:)
!!   to set the bias term, we subtract out contribution of this specific rotation (note that this is independent
!!   of any of the other oerpations performed in the latter part of this function
!    branchv(1:6) = branchv(1:6) - dc_di(imol)%btq(ttc,1:6)
!    call crossprod3(con1(:,rrr),pos2(:,rrr),cp1(:))
!    dc_di(imol)%btq(ttc,7) = dot_product(con1(:,rrr),branchv(4:6)) - dot_product(cp1(:),branchv(1:3))
!!   and add it back in
!    branchv(1:6) = branchv(1:6) + dc_di(imol)%btq(ttc,1:6)
  end do
!
! increment particle acc.s
  if (do_init.EQV..true.) then
    cart_v(atmol(imol,1):atmol(imol,2),1) = dc_di(imol)%v(1)
    cart_v(atmol(imol,1):atmol(imol,2),2) = dc_di(imol)%v(2)
    cart_v(atmol(imol,1):atmol(imol,2),3) = dc_di(imol)%v(3)
  end if
  do k=1,ilim
    if (izrot(ati)%alsz2.eq.0) then
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
    else
      ii = izrot(ati)%diffis(k,1)
      ik = izrot(ati)%diffis(k,2)
    end if
!    cart_a(ii:ik,1) = cart_a(ii:ik,1) + (branchv(10)-branchv(7))*x(ii:ik) - branchv(11)*y(ii:ik) - branchv(13)*z(ii:ik) + &
! &branchv(14)
!    cart_a(ii:ik,2) = cart_a(ii:ik,2) - branchv(11)*x(ii:ik) + (branchv(10)-branchv(8))*y(ii:ik) - branchv(12)*z(ii:ik) + &
! &branchv(15)
!    cart_a(ii:ik,3) = cart_a(ii:ik,3) - branchv(13)*x(ii:ik) - branchv(12)*y(ii:ik) + (branchv(10)-branchv(9))*z(ii:ik) + &
! &branchv(16)
    cart_v(ii:ik,1) = cart_v(ii:ik,1) + branchv(18)*z(ii:ik) - branchv(19)*y(ii:ik) + branchv(20)
    cart_v(ii:ik,2) = cart_v(ii:ik,2) + branchv(19)*x(ii:ik) - branchv(17)*z(ii:ik) + branchv(21)
    cart_v(ii:ik,3) = cart_v(ii:ik,3) + branchv(17)*y(ii:ik) - branchv(18)*x(ii:ik) + branchv(22)
  end do
!
!  do k=1,izrot(ati)%alsz2
!    ii = izrot(ati)%diffis(k,1)
!    ik = izrot(ati)%diffis(k,2)
!!   and now we remove the net accelerations for diff set only
!    do i=ii,ik
!      pos1(1) = x(i)
!      pos1(2) = y(i)
!      pos1(3) = z(i)
!      branchv(1:3) = branchv(1:3) - mass(i)*cart_a(i,:)
!      call crossprod3(pos1,cart_a(i,:),cp1(:))
!      branchv(4:6) = branchv(4:6) - mass(i)*cp1(:)
!    end do
!  end do
!
!  if (izrot(ati)%treevs(3).gt.0) then
!!   we cannot transfer the entire set to every branch, instead use branch-specific terms
!!   currently, these need to be recomputed (except main branch, hopefully the largest)
!    bv2(7:12) = 0.0
!    do kk=1,izrot(ati)%treevs(3)
!      do j=dc_di(imol)%maxntor,1,-1
!        ii = dc_di(imol)%recurs(j,1)
!        if ((izrot(ii)%treevs(1).eq.0).AND.(izrot(ii)%treevs(2).ne.izrot(ati)%treevs(2)).AND.&
! &(izrot(ii)%treevs(2).eq.izrot(ati)%treevs(5+kk))) then
!          ati2 = ii
!          bv2(1:6) = 0.0
!          do k=1,izrot(ati2)%alsz
!            ii = izrot(ati2)%rotis(k,1)
!            ik = izrot(ati2)%rotis(k,2)
!            do i=ii,ik
!!             net accel for torsions not including self
!              cp2(1) =  (branchv(10)-branchv(7))*x(i) - branchv(11)*y(i) - branchv(13)*z(i) + branchv(14)
!              cp2(2) = -branchv(11)*x(i) + (branchv(10)-branchv(8))*y(i) - branchv(12)*z(i) + branchv(15)
!              cp2(3) = -branchv(13)*x(i) - branchv(12)*y(i) + (branchv(10)-branchv(9))*z(i) + branchv(16)
!              cp1(1) = y(i)*cp2(3) - z(i)*cp2(2)
!              cp1(2) = z(i)*cp2(1) - x(i)*cp2(3)
!              cp1(3) = x(i)*cp2(2) - y(i)*cp2(1)
!              bv2(1:3) = bv2(1:3) + mass(i)*cp2(:)
!              bv2(4:6) = bv2(4:6) + mass(i)*cp1(:)
!            end do
!          end do
!          dc_di(imol)%valrecurs(1:6,izrot(ati)%treevs(5+kk)) = bv2(1:6)
!          bv2(7:12) = bv2(7:12) + bv2(1:6)
!          exit
!        else if ((izrot(ii)%treevs(1).eq.0).AND.(izrot(ii)%treevs(2).eq.izrot(ati)%treevs(2)).AND.&
! &(izrot(ii)%treevs(2).eq.izrot(ati)%treevs(5+kk))) then
!!         do nothing (added by diff)
!        end if
!      end do
!    end do
!    dc_di(imol)%valrecurs(1:6,izrot(ati)%treevs(2)) = branchv(1:6) - bv2(7:12)
!  end if
!
end
!
!--------------------------------------------------------------------------------------------------
!
! just inertia
!
subroutine Vrecurse_drot0(ati,imol,branchv)
!
  use iounit
  use forces
  use atoms
  use zmatrix
  use math
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: imol,ati
!
  integer ttc,i,k,ii,ik,ilim,rrr
  RTYPE dp1,pos2(3,3),pos3(3),con1(3,3),cp1(3),c2,c1,branchv(22),hlp2(3),con0(3),hlp3(3)
!
  pos3(1:3) = comm(imol,1:3)
  do rrr=1,3
    pos2(1:3,rrr) = comm(imol,1:3) + rgpcsm(imol,rrr,1:3)
    con0(:) = pos3(:)-pos2(:,rrr)
    c2 = 1.0/(dot_product(con0(:),con0(:)))
    c1 = sqrt(c2)
    con1(:,rrr) = -con0(:)*c1
  end do
!
  ilim = 0
  if (izrot(ati)%alsz2.eq.0) then
    ilim = izrot(ati)%alsz
    if (allocated(izrot(ati)%rotis).EQV..false.) then
      write(ilog,*) 'Fatal. Called Vrecurse_drot0(...) for a line, for which the required atom set is not &
 &allocated. This is a bug.'
      call fexit()
    end if
    branchv(:) = 0.0
  else
    ilim = izrot(ati)%alsz2
  end if
!
  do k=1,ilim
    if (izrot(ati)%alsz2.eq.0) then
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
    else
      ii = izrot(ati)%diffis(k,1)
      ik = izrot(ati)%diffis(k,2)
    end if
    do i=ii,ik
      hlp2(1) = x(i)
      hlp2(2) = y(i)
      hlp2(3) = z(i)
      branchv(17) = branchv(17) + mass(i)*dot_product(hlp2(:),hlp2(:))
    end do
    branchv(1) = branchv(1) + sum(mass(ii:ik)*x(ii:ik))
    branchv(2) = branchv(2) + sum(mass(ii:ik)*y(ii:ik))
    branchv(3) = branchv(3) + sum(mass(ii:ik)*z(ii:ik))
    branchv(4) = branchv(4) + sum(mass(ii:ik))
    branchv(5) = branchv(5) + sum(mass(ii:ik)*x(ii:ik)*x(ii:ik))
    branchv(6) = branchv(6) + sum(mass(ii:ik)*y(ii:ik)*y(ii:ik))
    branchv(7) = branchv(7) + sum(mass(ii:ik)*z(ii:ik)*z(ii:ik))
    branchv(8) = branchv(8) + sum(mass(ii:ik)*x(ii:ik)*y(ii:ik))
    branchv(9) = branchv(9) + sum(mass(ii:ik)*y(ii:ik)*z(ii:ik))
    branchv(10) = branchv(10) + sum(mass(ii:ik)*z(ii:ik)*x(ii:ik))
  end do
!
  call crossprod3(pos3(:),branchv(14:16),cp1(:))
  do rrr=1,3
    ttc = rrr + 3
    dc_di(imol)%f(ttc) = dc_di(imol)%f(ttc) + dot_product(con1(:,rrr),branchv(11:13)) - dot_product(con1(:,rrr),cp1(:))
    dp1 = dot_product(pos2(:,rrr),con1(:,rrr))
    hlp2(:) = con1(:,rrr)*con1(:,rrr)
    hlp3(1) = con1(1,rrr)*con1(2,rrr)
    hlp3(2) = con1(2,rrr)*con1(3,rrr)
    hlp3(3) = con1(3,rrr)*con1(1,rrr)
    dc_di(imol)%olddat(ttc,2) = dc_di(imol)%olddat(ttc,2) + branchv(4)*dot_product(pos2(:,rrr),pos2(:,rrr)) - &
 & 2.0*dot_product(pos2(:,rrr),branchv(1:3)) + branchv(17) - dp1*dp1*branchv(4) + &
 & 2.0*dp1*dot_product(con1(:,rrr),branchv(1:3)) - dot_product(hlp2(:),branchv(5:7)) - 2.0*dot_product(hlp3(:),branchv(8:10))
  end do
  dc_di(imol)%olddat(1:3,2) = molmass(moltypid(imol))
!
end
!
!-----------------------------------------------------------------------
!
subroutine fudge_masses()
!
  use forces
  use molecule
  use system
  use iounit
!
  implicit none
!
  integer imol,j
!
  if (fycxyz.ne.1) then
    write(ilog,*) 'Fatal. Called fudge_masses() while choice of degrees of freedom&
 & is not equivalent to torsional / rigid-body ones. This is a bug.'
    call fexit()
  end if
  do imol=1,nmol
    do j=1,size(dc_di(imol)%im)
      dc_di(imol)%im(j) = dc_di(imol)%im(j)/(1.0*dc_di(imol)%nattor(j))
    end do
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine set_IMDalign()
!
  use molecule
  use forces
  use iounit
  use sequen
  use atoms
  use movesets
  use zmatrix
  use fyoc
  use cutoffs, ONLY: molinfo
  use system, ONLY: fycxyz
!
  implicit none
!
  integer imol,tnat,nelgs,i,j,rs,kk,trat,kid,othst
  logical fndit,isbrnch
  integer, ALLOCATABLE:: tmpl(:,:)
  RTYPE rmass
!
  if ((align_NC.gt.4).OR.(align_NC.lt.1)) then
    write(ilog,*) 'Fatal. Called set_IMDalign() with unsupported chain alignment mode (',align_NC,&
   &'). This is a bug.'
    call fexit()
  end if
!
  do imol=1,nmol
    nelgs = 0
    tnat = atmol(imol,2) - atmol(imol,1) + 1
    if (tnat.eq.1) cycle ! no allocation for single atom species
    do i=atmol(imol,1)+3,atmol(imol,2)
      if (allocated(izrot(i)%rotis).EQV..true.) then
        nelgs = nelgs + 1
      end if
    end do
    dc_di(imol)%maxntor = nelgs
    dc_di(imol)%nattor(:) = tnat
    if (nelgs.le.0) cycle
    allocate(dc_di(imol)%recurs(nelgs,3))
    dc_di(imol)%recurs(:,:) = 0
    allocate(tmpl(maxval(izrot(atmol(imol,1):atmol(imol,2))%alsz)+4,2))
    fndit = .false.
    nelgs = 0  
!   populate pointer from recurs into izrot
    do i=atmol(imol,1)+3,atmol(imol,2)
      if (allocated(izrot(i)%rotis).EQV..true.) then
        nelgs = nelgs + 1
        dc_di(imol)%recurs(izrot(i)%treevs(4),1) = i
      end if
    end do
!   now check possibility of alignment, flag, and recompute rotation lists
    nelgs = 0
    do i=atmol(imol,1)+3,atmol(imol,2)
      if (allocated(izrot(i)%rotis).EQV..true.) then
        nelgs = nelgs + 1
        trat = 0
        rmass = 0.0
        do kk=1,izrot(i)%alsz
          trat = trat + izrot(i)%rotis(kk,2) - izrot(i)%rotis(kk,1) + 1
          rmass = rmass + sum(mass(izrot(i)%rotis(kk,1):izrot(i)%rotis(kk,2)))
        end do
        if ((izrot(i)%treevs(2).gt.1).AND.(izrot(i)%treevs(1).gt.0)) then
!         this is a branch
        else if ((izrot(i)%treevs(1).eq.0).AND.(izrot(i)%treevs(5).eq.0)) then
!         these should be all other terminal units not suitable for any alignment
        else
!         these should be on the main branch and are relevant for alignment
          isbrnch = .false.
          do j=izrot(i)%treevs(4),1,-1
            if ((izrot(dc_di(imol)%recurs(j,1))%treevs(1).eq.i).AND.(izrot(dc_di(imol)%recurs(j,1))%treevs(2).gt.1).AND.&
 &              (izrot(i)%treevs(1).eq.0).AND.(izrot(dc_di(imol)%recurs(j,1))%treevs(2).eq.izrot(i)%treevs(2))) then
              isbrnch = .true.
              exit
            end if
          end do
!         unless they are longer segments merging directly into the base filtered by this loop above
!         (note that all of this only works as long as izrot reflects the N-terminal building direction)
          if (isbrnch.EQV..true.) cycle
          if (rmass.le.11.5) cycle ! never align H-only cases
          if ((((tnat-trat-2).lt.trat).AND.(align_NC.ge.3)).OR.(align_NC.eq.2)) then
            dc_di(imol)%recurs(nelgs,2) = 1
            fndit = .true.
            kid = 0
            if (izrot(i)%rotis(1,1).gt.atmol(imol,1)) then
              kid = kid + 1
              tmpl(kid,1) = atmol(imol,1)
              tmpl(kid,2) = izrot(i)%rotis(1,1) - 1
            end if
            do kk=2,izrot(i)%alsz
              kid = kid + 1
              tmpl(kid,1) = izrot(i)%rotis(kk-1,2) + 1
              tmpl(kid,2) = izrot(i)%rotis(kk,1) - 1
            end do
            if (izrot(i)%rotis(izrot(i)%alsz,2).lt.atmol(imol,2)) then
              kid = kid + 1
              tmpl(kid,1) = izrot(i)%rotis(izrot(i)%alsz,2)+1
              tmpl(kid,2) = atmol(imol,2)
            end if
!           remove axis atoms
            do kk=1,kid
              if ((iz(1,i).gt.tmpl(kk,1)).AND.(iz(1,i).lt.tmpl(kk,2))) then
                kid = kid + 1
                tmpl(kid,1) = iz(1,i) + 1
                tmpl(kid,2) = tmpl(kk,2)
                tmpl(kk,2) = iz(1,i) - 1
              else if (iz(1,i).eq.tmpl(kk,1)) then
                tmpl(kk,1) = iz(1,i) + 1
              else if (iz(1,i).eq.tmpl(kk,2)) then
                tmpl(kk,2) = iz(1,i) - 1
              end if
              if ((iz(2,i).gt.tmpl(kk,1)).AND.(iz(2,i).lt.tmpl(kk,2))) then
                kid = kid + 1
                tmpl(kid,1) = iz(2,i) + 1
                tmpl(kid,2) = tmpl(kk,2)
                tmpl(kk,2) = iz(2,i) - 1
              else if (iz(2,i).eq.tmpl(kk,1)) then
                tmpl(kk,1) = iz(2,i) + 1
              else if (iz(2,i).eq.tmpl(kk,2)) then
                tmpl(kk,2) = iz(2,i) - 1
              end if
            end do            
!           rewrite list
            deallocate(izrot(i)%rotis)
            allocate(izrot(i)%rotis(kid,2))
            izrot(i)%alsz = kid
            izrot(i)%rotis(1:kid,:) = tmpl(1:kid,:)
          end if
        end if
      end if
    end do
    deallocate(tmpl)
    allocate(tmpl(nelgs,2))
    tmpl(:,:) = dc_di(imol)%recurs(:,1:2)
!   if necessary, redetermine rotation list structure
    if (fndit.EQV..true.) call parse_rotlsts(imol)
!   populate nattor and pointers from recurs structure to d.o.f structure
    nelgs = 0
    kid = 0
    dc_di(imol)%recurs(:,2) = 0 
    othst = othidxmol(moltypid(imol))
    do i=atmol(imol,1)+3,atmol(imol,2)
      if (allocated(izrot(i)%rotis).EQV..true.) then
        if (izrot(i)%treevs(2).gt.kid) kid = izrot(i)%treevs(2)
        nelgs = nelgs + 1
        dc_di(imol)%recurs(izrot(i)%treevs(4),1) = i
        do rs=max(rsmol(imol,1),atmres(i)-2),min(rsmol(imol,2),atmres(i)+2)
          if (i.eq.wline(rs)) dc_di(imol)%recurs(izrot(i)%treevs(4),3) = wnr(rs)
          if (i.eq.fline(rs)) dc_di(imol)%recurs(izrot(i)%treevs(4),3) = fnr(rs)
          if (i.eq.yline(rs)) dc_di(imol)%recurs(izrot(i)%treevs(4),3) = ynr(rs)
          do j=1,nchi(rs)
            if (i.eq.chiline(j,rs)) then
              dc_di(imol)%recurs(izrot(i)%treevs(4),3) = chinr(j,rs)
              exit
            end if 
          end do
          do j=1,nnucs(rs)
            if (i.eq.nucsline(j,rs)) then
              dc_di(imol)%recurs(izrot(i)%treevs(4),3) = nucsnr(j,rs)
              exit
            end if
          end do
        end do
        if (dc_di(imol)%recurs(izrot(i)%treevs(4),3).le.0) then ! must be an unsupported
          if (othst.le.0) then 
            write(ilog,*) 'Fatal. Setup of degrees of freedom not natively supported is inconsistent for &
 &molecule ',imol,'. This is a bug.'
            call fexit()
          end if
          dc_di(imol)%recurs(izrot(i)%treevs(4),3) = othst
          if (fycxyz.eq.1) then
            if (((dyn_integrator_ops(10).eq.0).OR.((dyn_integrator_ops(10).eq.1).AND.(seqtyp(atmres(iz(1,i))).ne.26)).OR.&
 &             ((dyn_integrator_ops(10).eq.2).AND.(seqtyp(atmres(iz(1,i))).eq.26))).AND.&
 &            (dc_di(imol)%frz(dc_di(imol)%recurs(izrot(i)%treevs(4),3)).EQV..false.)) then
              dc_di(imol)%frz(dc_di(imol)%recurs(izrot(i)%treevs(4),3)) = .true.
              n_constraints = n_constraints + 1
            end if
          end if
          othst = othst + 1
        end if
        if (dc_di(imol)%recurs(izrot(i)%treevs(4),3).gt.0) then
          if (tmpl(nelgs,2).eq.1) then
            dc_di(imol)%align(dc_di(imol)%recurs(izrot(i)%treevs(4),3)) = .true.
            dc_di(imol)%recurs(izrot(i)%treevs(4),2) = 1
          end if
          trat = 0
          do kk=1,izrot(i)%alsz
            trat = trat + izrot(i)%rotis(kk,2) - izrot(i)%rotis(kk,1) + 1
          end do
          dc_di(imol)%nattor(dc_di(imol)%recurs(izrot(i)%treevs(4),3)) = trat
        else
          if (tmpl(nelgs,2).eq.1) then
            dc_di(imol)%recurs(izrot(i)%treevs(4),2) = 1
          end if
        end if
      end if
    end do
    deallocate(tmpl)
    allocate(dc_di(imol)%valrecurs(22,kid))
  end do
!
  kk = 0 
  do imol=1,nmol
    molinfo(imol,2) = kk + 1
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
      kk = kk + 3
    else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
      kk = kk + 5
    else
      kk = kk + 6 + ntormol(moltypid(imol))
    end if
  end do
!
end
!
!---------------------------------------------------------------------
!
! this subroutine prealigns the torsions we wish to deal with in C-terminal alignment
! i.e., which generate a flailing arm on the N-terminal side
! the result is new xyz for the first three atoms in molecule imol, which a subsequent
! call to makexyz_formol then takes advantage of
!
subroutine IMD_prealign(imol,cycle_frz)
!
  use iounit
  use sequen
  use molecule
  use forces
  use polypep
  use movesets
  use fyoc
  use atoms
  use aminos
  use math
  use zmatrix
!
  implicit none
!
  integer i,j,imol,ati,atj,ttc
  RTYPE axl(3),refp(3)
  logical cycle_frz
!
  if ((molfrzidx(imol).gt.1).AND.(cycle_frz.EQV..true.)) return
!
  if (align_NC.eq.1) return
!
  ati = atmol(imol,1)
  atj = atmol(imol,1) + 2
!
  do i=1,dc_di(imol)%maxntor
    if (dc_di(imol)%recurs(i,3).le.0) cycle ! this torsion is not set up as d.o.f
    ttc = dc_di(imol)%recurs(i,3)
    if (dc_di(imol)%align(ttc).EQV..false.) cycle
    j = iz(1,dc_di(imol)%recurs(i,1))
    axl(1) = x(j) - x(iz(2,dc_di(imol)%recurs(i,1)))
    axl(2) = y(j) - y(iz(2,dc_di(imol)%recurs(i,1)))
    axl(3) = z(j) - z(iz(2,dc_di(imol)%recurs(i,1)))
    refp(1) = x(j)
    refp(2) = y(j)
    refp(3) = z(j)
    call alignC_onetor(ati,atj,refp,axl,dc_di(imol)%incr(ttc)/RADIAN)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a routine that does nothing but assemble TMD parameters into provided arrays
!
subroutine manage_imds(vars,mode,tpi)
!
  use iounit
  use forces
  use molecule
  use system
  use torsn
#ifdef ENABLE_THREADS
  use threads
  use cutoffs, ONLY: molinfo
#endif
!
  implicit none
!
  integer, INTENT(IN):: mode,tpi
!
  integer imol,j,ttc
  RTYPE vars(ndyntorsn+totrbd,3)
#ifdef ENABLE_THREADS
  integer sta,sto
!
  if (tpi.gt.0) then
    sta = thr_limits(61,tpi)
    sto = thr_limits(62,tpi)
  else
    sta =  1
    sto = nmol
  end if
#endif
!
!  allocate(vars(totrbd+ntorpuck,2))
!
#ifdef ENABLE_THREADS
  if ((sta.ge.1).AND.(sta.le.nmol)) then
    ttc = molinfo(sta,2) - 1
  else if (sto.ge.sta) then
    write(ilog,*) 'Fatal. Inconsistent thread bounds in manage_imds(...). This is a bug.'
    call fexit()
  end if
  do imol=sta,sto
#else
  ttc = 0
  do imol=1,nmol
#endif
    if (mode.eq.1) then
!     store in vars
      do j=1,3
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(j)
        vars(ttc,2) = dc_di(imol)%im(j)
        vars(ttc,3) = dc_di(imol)%olddat(j,2)
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(4)
        vars(ttc,2) = dc_di(imol)%im(4)
        vars(ttc,3) = dc_di(imol)%olddat(4,2)
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(5)
        vars(ttc,2) = dc_di(imol)%im(5)
        vars(ttc,3) = dc_di(imol)%olddat(5,2)
      else if ((atmol(imol,2)-atmol(imol,1)).gt.1) then
        do j=4,6
          ttc = ttc + 1
          vars(ttc,1) = dc_di(imol)%v(j)
          vars(ttc,2) = dc_di(imol)%im(j)
          vars(ttc,3) = dc_di(imol)%olddat(j,2)
        end do
      end if
      do j=1,ntormol(moltypid(imol))
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(6+j)
        vars(ttc,2) = dc_di(imol)%im(6+j)
        vars(ttc,3) = dc_di(imol)%olddat(6+j,2)
      end do
    else if (mode.eq.2) then
!     assign from vars
      do j=1,3
        ttc = ttc + 1
        dc_di(imol)%v(j) = vars(ttc,1)
        dc_di(imol)%im(j) = vars(ttc,2)
        dc_di(imol)%olddat(j,2) = vars(ttc,3)
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = ttc + 1
        dc_di(imol)%v(4) = vars(ttc,1)
        dc_di(imol)%im(4) = vars(ttc,2)
        dc_di(imol)%olddat(4,2) = vars(ttc,3)
        ttc = ttc + 1
        dc_di(imol)%v(5) = vars(ttc,1)
        dc_di(imol)%im(5) = vars(ttc,2)
        dc_di(imol)%olddat(5,2) = vars(ttc,3)
      else if ((atmol(imol,2)-atmol(imol,1)).gt.1) then
        do j=4,6
          ttc = ttc + 1
          dc_di(imol)%v(j) = vars(ttc,1)
          dc_di(imol)%im(j) = vars(ttc,2)
          dc_di(imol)%olddat(j,2) = vars(ttc,3)
        end do
      end if
      do j=1,ntormol(moltypid(imol))
        ttc = ttc + 1
        dc_di(imol)%v(6+j) = vars(ttc,1)
        dc_di(imol)%im(6+j) = vars(ttc,2)
        dc_di(imol)%olddat(6+j,2) = vars(ttc,3)
      end do
    else
      write(ilog,*) 'Fatal. Called manage_imds(...) with unknown mode (identifier is ',&
  &mode,'). This is a bug.'
      call fexit()
    end if
  end do
!
!  deallocate(vars)
!
  end
!
!------------------------------------------------------------------------------------------
!
subroutine rotlst_report()
!
  use molecule
  use zmatrix
  use iounit
  use atoms
  use params
  use system
  use forces
!
  implicit none
!
  integer imol,i,k,ixx,ixx2
  RTYPE tmpm
  logical have_tmd
!
  have_tmd = .false.
  if ((dyn_mode.ne.1).AND.(fycxyz.eq.1)) have_tmd = .true.
!
 32 format(' Mol.  # Ref.   Atoms   Mass      Frozen')
 33 format('   Atom # Rotat.  Mass    Unique Rank Order  Parent Ele. Frozen')
 34 format(2x,i7,1x,i6,1x,g10.4,i6,1x,i3,1x,i5,1x,i7,1x,a1,a1,a1,a1,2x,l1)
 35 format(2x,i7,1x,i6,1x,g10.4,i6,1x,i3,1x,i5,1x,i7,1x,a1,a1,a1,a1)
 36 format(i7,1x,i7,1x,i7,1x,g10.4,1x,6(l1))
!
  write(ilog,*)
  write(ilog,*) '--- Summary of Rotation List Setup by Molecule ---'

  do imol=1,nmol
    write(ilog,32)
    if (atmol(imol,1).eq.atmol(imol,2)) then
      if (have_tmd.EQV..true.) then
        write(ilog,36) imol,atmol(imol,1),atmol(imol,2)-atmol(imol,1)+1,sum(mass(atmol(imol,1):atmol(imol,2))),&
 &dc_di(imol)%frz(1:3)
      else
        write(ilog,36) imol,atmol(imol,1),atmol(imol,2)-atmol(imol,1)+1,sum(mass(atmol(imol,1):atmol(imol,2))) 
      end if
      cycle
    else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
      if (have_tmd.EQV..true.) then
        write(ilog,36) imol,atmol(imol,1),atmol(imol,2)-atmol(imol,1)+1,sum(mass(atmol(imol,1):atmol(imol,2))),&
 &dc_di(imol)%frz(1:5)
      else
        write(ilog,36) imol,atmol(imol,1),atmol(imol,2)-atmol(imol,1)+1,sum(mass(atmol(imol,1):atmol(imol,2))) 
      end if
      cycle
    else
      if (have_tmd.EQV..true.) then
        write(ilog,36) imol,atmol(imol,1),atmol(imol,2)-atmol(imol,1)+1,sum(mass(atmol(imol,1):atmol(imol,2))),&
 &dc_di(imol)%frz(1:6)
      else
        write(ilog,36) imol,atmol(imol,1),atmol(imol,2)-atmol(imol,1)+1,sum(mass(atmol(imol,1):atmol(imol,2))) 
      end if
    end if
!    format(i7,1x,i6,1x,
    if (maxval(izrot((atmol(imol,1)+1):atmol(imol,2))%alsz).gt.0) write(ilog,33)
    do i=atmol(imol,1)+1,atmol(imol,2)
      if (izrot(i)%alsz.gt.0) then
        ixx = 0 
        tmpm = 0.0
        do k=1,izrot(i)%alsz
          ixx = ixx + izrot(i)%rotis(k,2) - izrot(i)%rotis(k,1) + 1
          tmpm = tmpm + sum(mass(izrot(i)%rotis(k,1):izrot(i)%rotis(k,2)))
        end do
        if (izrot(i)%alsz2.gt.0) then
          ixx2 = 0
          do k=1,izrot(i)%alsz2
            ixx2 = ixx2 + izrot(i)%diffis(k,2) - izrot(i)%diffis(k,1) + 1
          end do
        else
          ixx2 = ixx
        end if
        if (have_tmd.EQV..true.) then
          write(ilog,34) i,ixx,tmpm,ixx2,izrot(i)%treevs(5),izrot(i)%treevs(4),izrot(i)%treevs(1),bio_code(b_type(iz(3,i)))(1:1),&
 &bio_code(b_type(iz(2,i)))(1:1),bio_code(b_type(iz(1,i)))(1:1),bio_code(b_type(i))(1:1),&
 &dc_di(imol)%frz(dc_di(imol)%recurs(izrot(i)%treevs(4),3))
        else
          write(ilog,35) i,ixx,tmpm,ixx2,izrot(i)%treevs(5),izrot(i)%treevs(4),izrot(i)%treevs(1),bio_code(b_type(iz(3,i)))(1:1),&
 &bio_code(b_type(iz(2,i)))(1:1),bio_code(b_type(iz(1,i)))(1:1),bio_code(b_type(i))(1:1)
        end if
      end if
    end do
  end do
  write(ilog,*) '---    End of Summary of Rotation List Setup   ---'
  write(ilog,*)
!
end
!
!------------------------------------------------------------------------------------------
!
