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
subroutine int_ldmove(istep,ndump,forcefirst)
!
  use iounit
  use math
  use molecule
  use forces
  use energies
  use movesets
  use system
  use units
  use mcsums
  use cutoffs
  use zmatrix
!
  implicit none
!
  integer imol,j,istep,ndump,aone,azero,ttc,rs,i
  RTYPE force3,t0,fr1,fr2
  RTYPE ompl,ommi,intS,bold,rndnew,normal
  RTYPE exgdt2,exgdt
  RTYPE, ALLOCATABLE:: mshfs(:,:)
  logical afalse,forcefirst
  integer(KIND=8) t1,t2
!
  afalse = .false.
!
 57   format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58   format('Mol. ',i6,' Res. ',i8,' :',g14.7)
!
  allocate(mshfs(3,nmol))
!
! set the constants for the Langevin
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/(fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
  intS = 0.0
!
! missing initialization
!
  aone = 1
  azero = 0
  nstep = nstep + 1
  mvcnt%nld = mvcnt%nld + 1
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
      dc_di(imol)%v(:) = 0.0
    end do
!
!   get the deterministic force
    call System_Clock(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    ens%insR(10) = ens%insU ! backup
    ens%insU = esave
    call System_Clock(t2)
    call cart2int_f(skip_frz,azero)
    time_energy = time_energy + 1.0*(t2 - t1)
!
!   set the frictional parameters
    do imol=1,nmol
      do j=1,size(dc_di(imol)%im)
        if ((dc_di(imol)%frz(j).EQV..true.)) cycle
        dc_di(imol)%ldp(j,1) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
        dc_di(imol)%ldp(j,2) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
        dc_di(imol)%ldp(j,3) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      end do
    end do
!
    call randomize_velocities(azero)
!
!   now loop over all molecules
    do imol=1,nmol
!
!     first: center-of-mass translation (linear motion)
      do j=1,3
!       increment velocity at time t - 1/2dt to t + 1/2dt
        if (dc_di(imol)%frz(j).EQV..true.) then
          cur_trans(j) = 0.0
          cycle
        end if
!       "increment" frictional/stochastic parameter to starting value
        dc_di(imol)%ldp(j,4) = sqrt(dc_di(imol)%ldp(j,1))
        dc_di(imol)%ldp(j,5) = normal()
        rndnew = dc_di(imol)%ldp(j,4)*dc_di(imol)%ldp(j,5)
!       increment velocity from time zero to 1/2dt
        dc_di(imol)%v(j) = exgdt2*(dc_di(imol)%v(j) + ompl*u_dyn_fconv*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j) + rndnew)
!       and finally increment position to dt
        cur_trans(j)=dc_di(imol)%v(j)*(1.0 - exgdt)/(fric_ga*exgdt2)
      end do
!     now all rotational motion (including RB) 
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc = dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
            dc_di(imol)%v(ttc) = 0.0
            cycle
          end if
!
          t0 = dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
!         "increment" frictional/stochastic parameter to starting value
          dc_di(imol)%ldp(ttc,4) = sqrt(dc_di(imol)%ldp(ttc,1))
          dc_di(imol)%ldp(ttc,5) = normal()
          rndnew = RADIAN*dc_di(imol)%ldp(ttc,4)*dc_di(imol)%ldp(ttc,5)
!         increment velocity from time zero to 1/2dt
          dc_di(imol)%v(ttc) = exgdt2*(dc_di(imol)%v(ttc) + ompl*t0/dc_di(imol)%im(ttc) + rndnew)
!         and propose torsion
          dc_di(imol)%incr(ttc) = dc_di(imol)%v(ttc)*(1.0 - exgdt)/(fric_ga*exgdt2)
        end do
!       now transfer
        cur_rot(1:3) = dc_di(imol)%incr(4:6)/RADIAN
!       note that this will edit the Z-matrix and populate ztorpr
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
      end if
!
      call makeref_formol(imol)
      call IMD_prealign(imol,skip_frz)
      call makexyz_formol(imol)
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
  else
!
    do imol=1,nmol
!     zero out internal coordinate forces and back up inertia at t_1
      dc_di(imol)%f(:) = 0.0
      dc_di(imol)%olddat(:,1) = dc_di(imol)%im(:)
    end do
!
!   get the deterministic force
    call System_Clock(t1)
    esave = force3(esterms,esterms_tr,esterms_lr,forcefirst)
    ens%insR(10) = ens%insU ! backup
    ens%insU = esave
    call System_Clock(t2)
    call cart2int_f(skip_frz,azero)
!    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
!
!   set the frictional parameters
    do imol=1,nmol
      do j=1,size(dc_di(imol)%im)
        if ((dc_di(imol)%frz(j).EQV..true.)) cycle
        dc_di(imol)%ldp(j,1) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
        dc_di(imol)%ldp(j,2) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
        dc_di(imol)%ldp(j,3) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      end do
    end do
!
!   now loop over all molecules
    do imol=1,nmol
!
!     first: center-of-mass translation (linear motion)
      do j=1,3
!       increment velocity at time t - 1/2dt to t + 1/2dt
        if (dc_di(imol)%frz(j).EQV..true.) then
          cur_trans(j) = 0.0
          cycle
        end if
!       increment frictional/stochastic parameters to current values (trailing)
        bold = dc_di(imol)%ldp(j,2)/dc_di(imol)%ldp(j,4)
        dc_di(imol)%ldp(j,4) = sqrt(dc_di(imol)%ldp(j,1) + dc_di(imol)%ldp(j,3) - bold*bold)
        rndnew = bold*dc_di(imol)%ldp(j,5)
        dc_di(imol)%ldp(j,5) = normal()
        rndnew = rndnew + dc_di(imol)%ldp(j,4)*dc_di(imol)%ldp(j,5)
!       now increment velocity at time t - 1/2dt to t + 1/2dt
        dc_di(imol)%v(j) = exgdt2*(exgdt2*dc_di(imol)%v(j) + u_dyn_fconv*dyn_dt*dc_di(imol)%f(j)/dc_di(imol)%im(j) + rndnew)
!       and finally increment position to t + dt
        cur_trans(j) = dc_di(imol)%v(j)*(1.0 - exgdt)/(fric_ga*exgdt2)
      end do
!     now all rotational motion (including RB) 
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc = dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
            dc_di(imol)%v(ttc) = 0.0
            cycle
          end if
!
!         dealing with fluc/diss and adhoc inertia correction
          t0 = dc_di(imol)%f(ttc)*dyn_dt*u_dyn_fconv_r
          bold = dc_di(imol)%ldp(ttc,2)/dc_di(imol)%ldp(ttc,4)
          dc_di(imol)%ldp(ttc,4) = sqrt(dc_di(imol)%ldp(ttc,1) + dc_di(imol)%ldp(ttc,3)-bold*bold)
          rndnew = RADIAN*bold*dc_di(imol)%ldp(ttc,5)
          dc_di(imol)%ldp(ttc,5) = normal()
          rndnew= rndnew + RADIAN*dc_di(imol)%ldp(ttc,4)*dc_di(imol)%ldp(ttc,5)
!         increment velocity at time t - 1/2dt to t + 1/2dt
          dc_di(imol)%v(ttc) = exgdt2*(exgdt2*sqrt(dc_di(imol)%olddat(ttc,1)/dc_di(imol)%im(ttc))*dc_di(imol)%v(ttc) + &
 &                                     t0/dc_di(imol)%im(ttc) + rndnew)
!         and propose torsion
          dc_di(imol)%incr(ttc) = dc_di(imol)%v(ttc)*(1.0 - exgdt)/(fric_ga*exgdt2)
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
      call IMD_prealign(imol,skip_frz)
      call makexyz_formol(imol)
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
  end if
!
  call get_ensv(ens%insV,azero)
  call drift_removal(fr1,fr2,azero)
  call prt_curens(istep,ens%insV,fr1,fr2)
!
  deallocate(mshfs)
!
end
!
!---------------------------------------------------------------------------
!
! this routine is only used to ensure that a run which was restarted from
! a will work properly. this is necessary because the non-first step LD
! algorithm requires ldp(i,4) and ldp(i,5)
!
subroutine init_ldps()
!
  use iounit
  use forces
  use molecule
  use system
  use units
  use movesets
!
  implicit none
!
  integer imol,j
  RTYPE ompl,ommi,normal
  RTYPE exgdt2,exgdt
!
! set the constants for the Langevin
  exgdt = exp(-fric_ga*dyn_dt)
  exgdt2 = exp(-fric_ga*0.5*dyn_dt)
  ompl = (exgdt - 1.0 + fric_ga*dyn_dt)/(fric_ga*dyn_dt*(1.0 - exgdt))
  ommi = 1.0 - ompl
  do imol=1,nmol
    do j=1,size(dc_di(imol)%im)
      if ((dc_di(imol)%frz(j).EQV..true.)) cycle
      dc_di(imol)%ldp(j,1) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ompl*ompl*fric_ga*dyn_dt + ompl - ommi)
      dc_di(imol)%ldp(j,2) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ompl*ommi*fric_ga*dyn_dt - ompl + ommi)
      dc_di(imol)%ldp(j,3) = (u_dyn_kb*kelvin/dc_di(imol)%im(j))*(2.0*ommi*ommi*fric_ga*dyn_dt + ompl - ommi)
      dc_di(imol)%ldp(j,4) = sqrt(dc_di(imol)%ldp(j,1))
      dc_di(imol)%ldp(j,5) = normal()
    end do
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a routine that does nothing but assemble TLD parameters into provided arrays
!
subroutine manage_ilds(vars,mode,tpi)
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
  RTYPE vars(totrbd+ndyntorsn,4)
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
!  allocate(vars(totrbd+ndyntorsn,4))
!
#ifdef ENABLE_THREADS
 if ((sta.ge.1).AND.(sta.le.nmol)) then
    ttc = molinfo(sta,2) - 1
  else if (sto.ge.sta) then
    write(ilog,*) 'Fatal. Inconsistent thread bounds in manage_ilds(...). This is a bug.'
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
        vars(ttc,3) = dc_di(imol)%ldp(j,4)
        vars(ttc,4) = dc_di(imol)%ldp(j,5)
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(4)
        vars(ttc,2) = dc_di(imol)%im(4)
        vars(ttc,3) = dc_di(imol)%ldp(4,4)
        vars(ttc,4) = dc_di(imol)%ldp(4,5)
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(5)
        vars(ttc,2) = dc_di(imol)%im(5)
        vars(ttc,3) = dc_di(imol)%ldp(5,4)
        vars(ttc,4) = dc_di(imol)%ldp(5,5)
      else if ((atmol(imol,2)-atmol(imol,1)).gt.1) then
        do j=4,6
          ttc = ttc + 1
          vars(ttc,1) = dc_di(imol)%v(j)
          vars(ttc,2) = dc_di(imol)%im(j)
          vars(ttc,3) = dc_di(imol)%ldp(j,4)
          vars(ttc,4) = dc_di(imol)%ldp(j,5)
        end do
      end if
      do j=1,ntormol(moltypid(imol))
        ttc = ttc + 1
        vars(ttc,1) = dc_di(imol)%v(6+j)
        vars(ttc,2) = dc_di(imol)%im(6+j)
        vars(ttc,3) = dc_di(imol)%ldp(6+j,4)
        vars(ttc,4) = dc_di(imol)%ldp(6+j,5)
      end do
    else if (mode.eq.2) then
!     assign from vars
      do j=1,3
        ttc = ttc + 1
        dc_di(imol)%v(j) = vars(ttc,1)
        dc_di(imol)%im(j) = vars(ttc,2)
        dc_di(imol)%ldp(j,4) = vars(ttc,3)
        dc_di(imol)%ldp(j,5) = vars(ttc,4)
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = ttc + 1
        dc_di(imol)%v(4) = vars(ttc,1)
        dc_di(imol)%im(4) = vars(ttc,2)
        dc_di(imol)%ldp(4,4) = vars(ttc,3)
        dc_di(imol)%ldp(4,5) = vars(ttc,4)
        ttc = ttc + 1
        dc_di(imol)%v(5) = vars(ttc,1)
        dc_di(imol)%im(5) = vars(ttc,2)
        dc_di(imol)%ldp(5,4) = vars(ttc,3)
        dc_di(imol)%ldp(5,5) = vars(ttc,4)
      else if ((atmol(imol,2)-atmol(imol,1)).gt.1) then
        do j=4,6
          ttc = ttc + 1
          dc_di(imol)%v(j) = vars(ttc,1)
          dc_di(imol)%im(j) = vars(ttc,2)
          dc_di(imol)%ldp(j,4) = vars(ttc,3)
          dc_di(imol)%ldp(j,5) = vars(ttc,4)
        end do
      end if
      do j=1,ntormol(moltypid(imol))
        ttc = ttc + 1
        dc_di(imol)%v(j+6) = vars(ttc,1)
        dc_di(imol)%im(j+6) = vars(ttc,2)
        dc_di(imol)%ldp(j+6,4) = vars(ttc,3)
        dc_di(imol)%ldp(j+6,5) = vars(ttc,4)
      end do
    else
      write(ilog,*) 'Fatal. Called manage_ilds(...) with unknown mode (identifier is ',&
  &mode,'). This is a bug.'
      call fexit()
    end if
  end do
!
!  deallocate(vars)
!
  end
!
!--------------------------------------------------------------------------------------
!
