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
! the strategy is the following:
! there are several types of actual boundary conditions identified through
! the integer-global bnd_type (periodic, hard)
! secondly, there are different shapes identified through the integer-global
! bnd_shape (cube, sphere)
! obviously, not all of them are compatible
! the routines herein handle:
! 1) boundary conditions for any type of com-displacement moves (involves internals!)
! 2) boundary conditions for distance evaluation
! 3) boundary conditions for randomization, image-finding, ...
! 4) boundary energy term
!
!----------------------------------------------------------------------
!
! a routine ensuring consistency in all boundary parameters
! we assume that the radius or the side-lengths have been updated for
! spherical and rectangular systems, respectively
!
subroutine update_bound(init)
!
  use iounit
  use system
  use math
  use units
  use forces
!
  implicit none
!
  RTYPE ovol
  logical init
!
  if (bnd_shape.eq.1) then
    if (init.EQV..true.) then
      ens%insV = bnd_params(1)*bnd_params(2)*bnd_params(3)
    else
      ovol = ens%insV
      ens%insV = bnd_params(1)*bnd_params(2)*bnd_params(3)
      bnd_pV = bnd_pV + (ens%insV-ovol)*(extpress/u_dyn_virconv)
    end if
  else if (bnd_shape.eq.2) then
    if (init.EQV..true.) then
      bnd_params(5) = bnd_params(4)*bnd_params(4)
      ens%insV = (4./3.)*PI*bnd_params(5)*bnd_params(4)
    else
      ovol = ens%insV
      bnd_params(5) = bnd_params(4)*bnd_params(4)
      ens%insV = (4./3.)*PI*bnd_params(5)*bnd_params(4)
      bnd_pV = bnd_pV + (ens%insV-ovol)*(extpress/u_dyn_virconv)
    end if
  else if (bnd_shape.eq.3) then
    if (init.EQV..true.) then
      bnd_params(5) = bnd_params(4)*bnd_params(4)
      ens%insV = PI*bnd_params(5)*bnd_params(6)
    else
      ovol = ens%insV
      bnd_params(5) = bnd_params(4)*bnd_params(4)
      ens%insV = PI*bnd_params(5)*bnd_params(6)
      bnd_pV = bnd_pV + (ens%insV-ovol)*(extpress/u_dyn_virconv)
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in updat&
 &e_bound() (code # ',bnd_shape,').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------
!
! a routine to get a valid translational displacement vector
!
subroutine trans_bound(imol,tvec)
!
  use iounit
  use molecule
  use system
  use movesets
!
  implicit none
!
  integer i,imol,k
  RTYPE tvec(3),como(3),dis2,random,tlen
  logical done
!
  do i=1,3
    como(i) = com(imol,i)
  end do
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      call ranvec(tvec)
      tlen = random()*trans_stepsz
      do i=1,3
        com(imol,i) = com(imol,i) + tlen*tvec(i)
      end do
      do i=1,3
        if (com(imol,i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-com(imol,i))/bnd_params(i)) + 1
          com(imol,i) = com(imol,i) + k*bnd_params(i)
        else if(com(imol,i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k = floor((com(imol,i)-(bnd_params(3+i)+bnd_params(i)))/bnd_params(i)) + 1
          com(imol,i) = com(imol,i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      call ranvec(tvec)
      tlen = random()*trans_stepsz
      do i=1,3
        com(imol,i) = com(imol,i) + tlen*tvec(i)
      end do
      if (com(imol,3).lt.bnd_params(3)-0.5*bnd_params(6)) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-com(imol,3))/bnd_params(6)) + 1
        com(imol,3) = com(imol,3) + k*bnd_params(6)
      else if (com(imol,3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((com(imol,3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        com(imol,3) = com(imol,3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in trans_bound() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
!   spherical system
    if (bnd_shape.eq.2) then
!     by filtering out moves that would leave the (hard) system boundary,
!     we might break the Markov chain, but hope to minimize boundary artifacts
      done = .false.
      do while (done.EQV..false.)
        call ranvec(tvec)
        tlen = random()*trans_stepsz
        do i=1,3
          tvec(i) = tvec(i)*tlen
        end do
        do i=1,3
          com(imol,i) = com(imol,i) + tvec(i)
        end do
        dis2 = 0.0
        do i=1,3
          dis2 = dis2 + (com(imol,i)-bnd_params(i))*&
 &                      (com(imol,i)-bnd_params(i))
        end do
        if (dis2.le.bnd_params(5)) then
          done = .true.
        else
          do i=1,3
            com(imol,i) = com(imol,i) - tvec(i)
          end do
        end if
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in trans_bound() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! soft wall BCs (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   any system: no constraints apply (all boundarys must be using a finite potential for nonspherical systems!)
    call ranvec(tvec)
    tlen = random()*trans_stepsz
    do i=1,3
      com(imol,i) = com(imol,i) + tlen*tvec(i)
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in trans_bound() (code # ',bnd_type,').'
    call fexit()
  end if
!
  do i=1,3
    tvec(i) = com(imol,i) - como(i)
  end do
!
end
!
!----------------------------------------------------------------------
!
! a routine to correct a translational displacement vector for BC
!
subroutine trans_bound2(imol,tvec,shiftv)
!
  use iounit
  use molecule
  use system
  use movesets
!
  implicit none
!
  integer i,imol,k
  RTYPE tvec(3),como(3),comn(3),shiftv(3)
!
  do i=1,3
    como(i) = comm(imol,i)
  end do
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        comn(i) = comm(imol,i) + tvec(i)
      end do
      do i=1,3
        if (comn(i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-comn(i))/bnd_params(i)) + 1
          shiftv(i) = k*bnd_params(i)
          comn(i) = comn(i) + shiftv(i)
        else if(comn(i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k = floor((comn(i)-(bnd_params(3+i)+bnd_params(i)))/bnd_params(i)) + 1
          shiftv(i) = -k*bnd_params(i)
          comn(i) = comn(i) + shiftv(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      do i=1,3
        comn(i) = comm(imol,i) + tvec(i)
      end do
      if (comn(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-comn(3))/bnd_params(6)) + 1
        shiftv(3)= k*bnd_params(6)
        comn(3) = comn(3) + shiftv(3)
      else if (comn(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((comn(3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        shiftv(3) = -k*bnd_params(6)
        comn(3) = comn(3) + shiftv(3)
      end if

    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in trans_bound2() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Hard-wall boundary condition is currently not supported. This is a bug.'
    call fexit()
! soft wall BC (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   any system: no constraints apply
    do i=1,3
      comn(i) = comm(imol,i) + tvec(i)
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in trans_bound2() (code # ',bnd_type,').'
    call fexit()
  end if
!
  do i=1,3
    tvec(i) = comn(i) - como(i)
  end do
!
end
!
!----------------------------------------------------------------------
!
! same as trans_bound(...) operating on a vector instead
!
subroutine trans_bound3(comv,tvec)
!
  use iounit
  use molecule
  use system
  use movesets
!
  implicit none
!
  integer i,k
  RTYPE tvec(3),como(3),comv(3),random,tlen
!
  do i=1,3
    como(i) = comv(i)
  end do
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      call ranvec(tvec)
      tlen = random()*trans_stepsz
      do i=1,3
        comv(i) = comv(i) + tlen*tvec(i)
      end do
      do i=1,3
        if (comv(i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-comv(i))/bnd_params(i)) + 1
          comv(i) = comv(i) + k*bnd_params(i)
        else if(comv(i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k = floor((comv(i)-(bnd_params(3+i)+bnd_params(i)))/&
 &                                             bnd_params(i)) + 1
          comv(i) = comv(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      call ranvec(tvec)
      tlen = random()*trans_stepsz
      do i=1,3
        comv(i) = comv(i) + tlen*tvec(i)
      end do
      if (comv(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-comv(3))/bnd_params(6)) + 1
        comv(3) = comv(3) + k*bnd_params(6)
      else if (comv(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((comv(3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        comv(3) = comv(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in trans_bound3() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Hard-wall boundary condition is not yet &
 &supported for cluster-based rigid body moves. This is an omission bug.'
    call fexit()
! soft wall BC (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   spherical system: no constraints apply
    call ranvec(tvec)
    tlen = random()*trans_stepsz
    do i=1,3
      comv(i) = comv(i) + tlen*tvec(i)
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in trans_bound3() (code # ',bnd_type,').'
    call fexit()
  end if
!
  do i=1,3
    tvec(i) = comv(i) - como(i)
  end do
!
end
!
!----------------------------------------------------------------------
!
! same as trans_bound2 operating on com, not comm, but not exporting shift vector
!
subroutine trans_bound4(imol,tvec)
!
  use iounit
  use molecule
  use system
  use movesets
!
  implicit none
!
  integer i,imol,k
  RTYPE tvec(3),como(3),comn(3)
!
  do i=1,3
    como(i) = com(imol,i)
  end do
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        comn(i) = com(imol,i) + tvec(i)
      end do
      do i=1,3
        if (comn(i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-comn(i))/bnd_params(i)) + 1
          comn(i) = comn(i) + k*bnd_params(i)
        else if(comn(i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k = floor((comn(i)-(bnd_params(3+i)+bnd_params(i)))/bnd_params(i)) + 1
          comn(i) = comn(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      do i=1,3
        comn(i) = com(imol,i) + tvec(i)
      end do
      if (comn(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-comn(3))/bnd_params(6)) + 1
        comn(3) = comn(3) + k*bnd_params(6)
      else if (comn(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((comn(3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        comn(3) = comn(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in trans_bound4() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Hard-wall boundary condition is currently not supported. This is a bug.'
    call fexit()
! soft wall BC (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   spherical system: no constraints apply
    do i=1,3
      comn(i) = com(imol,i) + tvec(i)
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in trans_bound4() (code # ',bnd_type,').'
    call fexit()
  end if
!
  do i=1,3
    tvec(i) = comn(i) - como(i)
  end do
!
end
!
!----------------------------------------------------------------------
!
! a routine to correct an out-of-box molecule to its proper image (based on comm, com shifted by the same)
!
subroutine update_image(imol)
!
  use iounit
  use molecule
  use system
  use movesets
  use atoms
!
  implicit none
!
  integer i,imol,k
  RTYPE tvec(3)
!
  tvec(:) = 0.0
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (comm(imol,i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-comm(imol,i))/bnd_params(i)) + 1
          tvec(i) = k*bnd_params(i)
        else if(comm(imol,i).gt.(bnd_params(3+i)+bnd_params(i)))then
          k = floor((comm(imol,i)-(bnd_params(3+i)+bnd_params(i)))/bnd_params(i)) + 1
          tvec(i) = -k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (comm(imol,3).lt.(bnd_params(3)-0.5*bnd_params(6))) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-comm(imol,3))/bnd_params(6)) + 1
        tvec(3) = k*bnd_params(6)
      else if(comm(imol,3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((comm(imol,3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        tvec(3) = -k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in update_image() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Hard-wall boundary condition is currently not supported. This is a bug.'
    call fexit()
! soft wall BC (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in update_image() (code # ',bnd_type,').'
    call fexit()
  end if
!
  comm(imol,:) = comm(imol,:) + tvec(:)
  com(imol,:) = com(imol,:) + tvec(:)
  do i=atmol(imol,1),atmol(imol,2)
    x(i) = x(i) + tvec(1)
    y(i) = y(i) + tvec(2)
    z(i) = z(i) + tvec(3)
  end do
!
end
!
!----------------------------------------------------------------------
!
! a routine to correct an out-of-box molecule to its proper image (based on com, comm shifted by the same)
!
subroutine update_image2(imol)
!
  use iounit
  use molecule
  use system
  use movesets
  use atoms
!
  implicit none
!
  integer i,imol,k
  RTYPE tvec(3)
!
  tvec(:) = 0.0
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (com(imol,i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-com(imol,i))/bnd_params(i)) + 1
          tvec(i) = k*bnd_params(i)
        else if(com(imol,i).gt.(bnd_params(3+i)+bnd_params(i)))then
          k = floor((com(imol,i)-(bnd_params(3+i)+bnd_params(i)))/bnd_params(i)) + 1
          tvec(i) = -k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (com(imol,3).lt.(bnd_params(3)-0.5*bnd_params(6))) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-com(imol,3))/bnd_params(6)) + 1
        tvec(3) = k*bnd_params(6)
      else if(com(imol,3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((com(imol,3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        tvec(3) = -k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in update_image2() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Hard-wall boundary condition is currently not supported. This is a bug.'
    call fexit()
! soft wall BC (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in update_image() (code # ',bnd_type,').'
    call fexit()
  end if
!
  com(imol,:) = com(imol,:) + tvec(:)
  comm(imol,:) = comm(imol,:) + tvec(:)
  do i=atmol(imol,1),atmol(imol,2)
    x(i) = x(i) + tvec(1)
    y(i) = y(i) + tvec(2)
    z(i) = z(i) + tvec(3)
  end do
!
end
!
!----------------------------------------------------------------------
!
! a routine to update RB coordinates after an internal (torsional) move,
! which correctly accounts for BC
!
subroutine internal_bound(imol,tpi)
!
  use iounit
  use molecule
  use system
  use atoms
  use energies, ONLY: use_POLY
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: imol,tpi
!
  integer i,k,stx(2)
  RTYPE tvec(3)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL,jfl(4)
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, internal_bound(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if ((tpi.gt.0).AND.(thr_mlgix(imol).gt.0)) then
    stx(1:2) = mlg_limits(thr_mlgix(imol),1:2,tpi)
  else if (tpi.le.1) then
    stx(1:2) = atmol(imol,1:2)
  else
    stx(1) = 1
    stx(2) = 0
  end if
#else
  stx(1:2) = atmol(imol,1:2)
#endif
!
  tvec(:) = 0.0
#ifdef ENABLE_THREADS
  if ((use_POLY.EQV..false.).AND.(tpi.gt.0).AND.(thr_mlgix(imol).gt.0)) then
    jfl(1) = .true.
    jfl(2:4) = .false.
    jfl(3) = .true.
    call molops_threads_geo(thr_mlgix(imol),tpi,jfl)
  else if (use_POLY.EQV..false.) then
!$OMP SINGLE
    call update_rigid(imol)
!$OMP END SINGLE NOWAIT
  end if
!$OMP BARRIER
#else
  if (use_POLY.EQV..false.) call update_rigid(imol)
#endif
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (com(imol,i).lt.bnd_params(3+i)) then
          k = floor((bnd_params(3+i)-com(imol,i))/bnd_params(i)) + 1
          tvec(i) = k*bnd_params(i)
        else if(com(imol,i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k = floor((com(imol,i)-(bnd_params(3+i)+bnd_params(i)))/bnd_params(i)) + 1
          tvec(i) = -k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (com(imol,3).lt.(bnd_params(3)-0.5*bnd_params(6))) then
        k = floor((bnd_params(3)-0.5*bnd_params(6)-com(imol,3))/bnd_params(6)) + 1
        tvec(3) = k*bnd_params(6)
      else if(com(imol,3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k = floor((com(imol,3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        tvec(3) = -k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in internal_bound() (code # ',bnd_shape,'). This is an omission bug.'
      call fexit()
    end if
! hard wall BC (HWBC)
  else if (bnd_type.eq.2) then
    write(ilog,*) 'Fatal. Hard-wall boundary condition not yet supported for torsional moves.'
    call fexit()
! soft wall BC (SWBC)
  else if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in internal_bound() (code # ',bnd_type,').'
    call fexit()
  end if
!
  if ((bnd_type.eq.1).AND.(stx(2).ge.stx(1))) then
    x(stx(1):stx(2)) = x(stx(1):stx(2)) + tvec(1)
    y(stx(1):stx(2)) = y(stx(1):stx(2)) + tvec(2)
    z(stx(1):stx(2)) = z(stx(1):stx(2)) + tvec(3)
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
  com(imol,:) = com(imol,:) + tvec(:)
  comm(imol,:) = comm(imol,:) + tvec(:)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!------------------------------------------------------------------------------------
!
! a routine to place a molecule somewhere in the box
!
subroutine randomize_trans(imol,tvec)
!
  use iounit
  use molecule
  use system
  use atoms
!
  implicit none
!
  integer i,imol
  RTYPE dis2,como(3),random,tvec(3)
  logical done
!
  do i=1,3
    como(i) = com(imol,i)
  end do
!
  if (bnd_shape.eq.1) then
    do i=1,3
      com(imol,i) = random()*bnd_params(i)
      com(imol,i) = com(imol,i) + bnd_params(3+i)
    end do
  else if (bnd_shape.eq.2) then ! sphere
    done = .false.
    do while (done.EQV..false.)
      do i=1,3
        com(imol,i) = random()*2.0*bnd_params(4) - bnd_params(4)
        com(imol,i) = com(imol,i) + bnd_params(i)
      end do
      dis2 = 0.0
      do i=1,3
        dis2 = dis2 + (com(imol,i)-bnd_params(i))*(com(imol,i)-bnd_params(i))
      end do
      if (dis2.le.bnd_params(5)) then
        done = .true.
      end if
    end do
  else if (bnd_shape.eq.3) then ! cylinder
    done = .false.
    do while (done.EQV..false.)
      com(imol,1) = bnd_params(4)*(2.0*random()-1.0)
      com(imol,2) = bnd_params(4)*(2.0*random()-1.0)
      com(imol,3) = bnd_params(6)*(random()-0.5)
      com(imol,:) = com(imol,:) + bnd_params(1:3)
      dis2 = 0.0
      do i=1,2
        dis2 = dis2 + (com(imol,i)-bnd_params(i))*(com(imol,i)-bnd_params(i))
      end do
      if (dis2.le.bnd_params(5)) then
        done = .true.
      end if
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in ran&
 &domize_trans() (code # ',bnd_shape,').'
    call fexit()
  end if
!
  do i=1,3
    tvec(i) = com(imol,i) - como(i)
  end do
!
end
!
!------------------------------------------------------------------------------------
!
! the same routine operating on a supplied vector
!
subroutine randomize_trans2(comv,tvec)
!
  use iounit
  use molecule
  use system
  use atoms
!
  implicit none
!
  integer i
  RTYPE dis2,como(3),comv(3),random,tvec(3)
  logical done
!
  do i=1,3
    como(i) = comv(i)
  end do
!
! rectangular box
  if (bnd_shape.eq.1) then
    do i=1,3
      comv(i) = random()*bnd_params(i)
      comv(i) = comv(i) + bnd_params(3+i)
    end do
! spherical system
  else if (bnd_shape.eq.2) then
    done = .false.
    do while (done.EQV..false.)
      do i=1,3
        comv(i) = random()*2.0*bnd_params(4) - bnd_params(4)
        comv(i) = comv(i) + bnd_params(i)
      end do
      dis2 = 0.0
      do i=1,3
        dis2 = dis2 + (comv(i)-bnd_params(i))*(comv(i)-bnd_params(i))
      end do
      if (dis2.le.bnd_params(5)) then
        done = .true.
      end if
    end do
! cylinder
  else if (bnd_shape.eq.3) then
    done = .false.
    do while (done.EQV..false.)
      comv(1) = bnd_params(4)*(2.0*random()-1.0)
      comv(2) = bnd_params(4)*(2.0*random()-1.0)
      comv(3) = bnd_params(6)*(random()-0.5)
      comv(:) = comv(:) + bnd_params(1:3)
      dis2 = 0.0
      do i=1,2
        dis2 = dis2 + (comv(i)-bnd_params(i))*(comv(i)-bnd_params(i))
      end do
      if (dis2.le.bnd_params(5)) then
        done = .true.
      end if
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in randomize_trans2() (code # ',bnd_shape,').'
    call fexit()
  end if
!
  do i=1,3
    tvec(i) = comv(i) - como(i)
  end do
!
end
!
!--------------------------------------------------------------------------
!
! this function takes the reference atoms for residues rs1,rs2 and computes
! the correct image distance vector between the two residues
! this can in turn be applied (using dis_bound4) to generate consistent image
! image interactions for residues rs1,rs2
!
! note that it IS ABSOLUTELY NECESSARY TO PRESERVE THE ORDER OF TAKING THE
! DISTANCE IN SUBSEQUENT CALLS TO dis_bound4/5, i.e., if atom ii comes from res.
! rs1, and atom kk from rs2, then the call to dis_bound4 must be
! call dis_bound4(ii,kk,svec,dis2), otherwise distances are nonsensical
!
! values are corrected for intramolecular terms (not shifted) including possible warnings
!
subroutine dis_bound_rs(rs1,rs2,svec)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer i,k,ii,kk
  RTYPE dvec(3),dis2,ddd(3),svec(3)
  logical viol
!
 33   format('Boundary violation for residues ',i6,' and ',i6,' in dis_bound_rs(...):')
 34   format('Image distance (',g14.6,') is less than long-range cutoff (',g14.6,').')
!
! init.
  svec(:) = 0.0
!
  ii = refat(rs1)
  kk = refat(rs2)
!
! distance vector
  viol = .false.
  dvec(1) = x(kk) - x(ii)
  dvec(2) = y(kk) - y(ii)
  dvec(3) = z(kk) - z(ii)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
            if (viol.EQV..false.) viol = .true.
          else
            svec(i) = k*bnd_params(i)
          end if
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
            if (viol.EQV..false.) viol = .true.
          else
            svec(i) = -k*bnd_params(i)
          end if
        end if
      end do
      if ((viol.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
        do i=1,3
          if (dvec(i).lt.(-0.5*bnd_params(i))) then
            k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i))+1
            ddd(i) = dvec(i) + k*bnd_params(i)
          else if(dvec(i).gt.0.5*bnd_params(i)) then
            k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i))+1
            ddd(i) = dvec(i) - k*bnd_params(i)
          else
            ddd(i) = dvec(i)
          end if
        end do
        dis2 = dot_product(ddd(:),ddd(:))
        if (dis2.lt.mcel_cutoff2) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(BND_WARNING)
#endif
          bnd_wrncnt(2) = bnd_wrncnt(2) + 1
          if (bnd_wrncnt(2).eq.bnd_wrnlmt(2)) then
            write(ilog,33) rs1,rs2
            write(ilog,34) sqrt(dis2),mcel_cutoff
            write(ilog,*) 'This was warning #',bnd_wrncnt(2),' of this type not all of which may be displayed.'
            if (10.0*bnd_wrnlmt(2).gt.0.5*HUGE(bnd_wrnlmt(2))) then
              bnd_wrncnt(2) = 0
            else
              bnd_wrnlmt(2) = bnd_wrnlmt(2)*10
            end if
          end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(BND_WARNING)
#endif
        end if
      end if
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
          if (viol.EQV..false.) viol = .true.
        else
          svec(3) = k*bnd_params(6)
        end if
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
          if (viol.EQV..false.) viol = .true.
        else
          svec(3) = -k*bnd_params(6)
        end if
      end if
      if ((viol.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
        ddd(:) = dvec(:)
        if (dvec(3).lt.(-0.5*bnd_params(6))) then
          k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6))+1
          ddd(3) = dvec(3) + k*bnd_params(6)
        else if(dvec(3).gt.0.5*bnd_params(6)) then
          k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6))+1
          ddd(3) = dvec(3) - k*bnd_params(6)
        end if
        dis2 = dot_product(ddd(:),ddd(:))
        if (dis2.lt.mcel_cutoff2) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(BND_WARNING)
#endif
          bnd_wrncnt(2) = bnd_wrncnt(2) + 1
          if (bnd_wrncnt(2).eq.bnd_wrnlmt(2)) then
            write(ilog,33) rs1,rs2
            write(ilog,34) sqrt(dis2),mcel_cutoff
            write(ilog,*) 'This was warning #',bnd_wrncnt(2),' of this type not all of which may be displayed.'
            if (10.0*bnd_wrnlmt(2).gt.0.5*HUGE(bnd_wrnlmt(2))) then
              bnd_wrncnt(2) = 0
            else
              bnd_wrnlmt(2) = bnd_wrnlmt(2)*10
            end if
          end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(BND_WARNING)
#endif
        end if
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound_rs() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound_rs() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this is the same routine, only trimmed down to the bare minimum
! correct for intramolecular terms, but no warnings produced
!
subroutine dis_bound_rs2(rs1,rs2,svec)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer i,k,ii,kk
  RTYPE dvec(3),svec(3)
!
  svec(:) = 0.0
!
  ii = refat(rs1)
  kk = refat(rs2)
!
! distance vector
  dvec(1) = x(kk) - x(ii)
  dvec(2) = y(kk) - y(ii)
  dvec(3) = z(kk) - z(ii)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
            svec(i) = k*bnd_params(i)
          end if
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
            svec(i) = -k*bnd_params(i)
          end if
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
          svec(3) = k*bnd_params(6)
        end if
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
          svec(3) = -k*bnd_params(6)
        end if
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound_rs2() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound_rs2() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! again the same (minimal) routine, which reports the actual distance as well
! again, correct for intramolecular terms, but no warnings produced
!
subroutine dis_bound_rs3(rs1,rs2,svec,dis)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none
!
  integer, INTENT(IN):: rs1,rs2
!
  integer i,k,ii,kk
  RTYPE dvec(3),svec(3),dis
!
  svec(:) = 0.0
!
  ii = refat(rs1)
  kk = refat(rs2)
!
! distance vector
  dvec(1) = x(kk) - x(ii)
  dvec(2) = y(kk) - y(ii)
  dvec(3) = z(kk) - z(ii)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
            svec(i) = k*bnd_params(i)
          end if
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
            svec(i) = -k*bnd_params(i)
          end if
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
          svec(3) = k*bnd_params(6)
        end if
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
          svec(3) = -k*bnd_params(6)
        end if
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound_rs3() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound_rs3() (code # ',bnd_type,').'
    call fexit()
  end if
!
  dvec(:) = dvec(:) + svec(:)
  dis = sqrt(sum(dvec(:)*dvec(:)))
!
end
!
!--------------------------------------------------------------------------
!
! functions for individual atom pairs: corrected for intramolecular terms with possible warnings
!
subroutine dis_bound(ii,kk,dis2)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none
!
  integer, INTENT(IN):: ii,kk
!
  integer i,k
  RTYPE dvec(3),dis2,ddd(3)
  logical viol
!
 33   format('Boundary violation for atoms ',i7,' and ',i7,' in dis_bound(...):')
 34   format('Image distance (',g14.6,') is less than long-range cutoff (',g14.6,').')
! distance vector
  viol = .false.
  dvec(1) = x(kk) - x(ii)
  dvec(2) = y(kk) - y(ii)
  dvec(3) = z(kk) - z(ii)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
            if (viol.EQV..false.) viol = .true.
          else
            dvec(i) = dvec(i) + k*bnd_params(i)
          end if
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
            if (viol.EQV..false.) viol = .true.
          else
            dvec(i) = dvec(i) - k*bnd_params(i)
          end if
        end if
      end do
      if ((viol.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
        do i=1,3
          if (dvec(i).lt.(-0.5*bnd_params(i))) then
            k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i))+1
            ddd(i) = dvec(i) + k*bnd_params(i)
          else if(dvec(i).gt.0.5*bnd_params(i)) then
            k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i))+1
            ddd(i) = dvec(i) - k*bnd_params(i)
          else
            ddd(i) = dvec(i)
          end if
        end do
        dis2 = dot_product(ddd(:),ddd(:))
        if (dis2.lt.mcel_cutoff2) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(BND_WARNING)
#endif
          bnd_wrncnt(1) = bnd_wrncnt(1) + 1
          if (bnd_wrncnt(1).eq.bnd_wrnlmt(1)) then
            write(ilog,33) ii,kk
            write(ilog,34) sqrt(dis2),mcel_cutoff
            write(ilog,*) 'This was warning #',bnd_wrncnt(1),' of this type not all of which may be displayed.'
            if (10.0*bnd_wrnlmt(1).gt.0.5*HUGE(bnd_wrnlmt(1))) then
              bnd_wrncnt(1) = 0
            else
              bnd_wrnlmt(1) = bnd_wrnlmt(1)*10
            end if
          end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(BND_WARNING)
#endif
        end if
      end if
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
          if (viol.EQV..false.) viol = .true.
        else
          dvec(3) = dvec(3) + k*bnd_params(6)
        end if
      else if (dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        if (molofrs(atmres(ii)).eq.molofrs(atmres(kk))) then
          if (viol.EQV..false.) viol = .true.
        else
          dvec(3) = dvec(3) - k*bnd_params(6)
        end if
      end if
      if ((viol.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
        ddd(:) = dvec(:)
        if (dvec(3).lt.(-0.5*bnd_params(6))) then
          k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6))+1
          ddd(3) = dvec(3) + k*bnd_params(6)
        else if(dvec(3).gt.0.5*bnd_params(6)) then
          k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6))+1
          ddd(3) = dvec(3) - k*bnd_params(6)
        end if
        dis2 = dot_product(ddd(:),ddd(:))
        if (dis2.lt.mcel_cutoff2) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(BND_WARNING)
#endif
          bnd_wrncnt(1) = bnd_wrncnt(1) + 1
          if (bnd_wrncnt(1).eq.bnd_wrnlmt(1)) then
            write(ilog,33) ii,kk
            write(ilog,34) sqrt(dis2),mcel_cutoff
            write(ilog,*) 'This was warning #',bnd_wrncnt(1),' of this type not all of which may be displayed.'
            if (10.0*bnd_wrnlmt(1).gt.0.5*HUGE(bnd_wrnlmt(1))) then
              bnd_wrncnt(1) = 0
            else
              bnd_wrnlmt(1) = bnd_wrnlmt(1)*10
            end if
          end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(BND_WARNING)
#endif
        end if
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound() (code # ',bnd_type,').'
    call fexit()
  end if
!
  dis2 = 0.0
  do i=1,3
    dis2 = dis2 + dvec(i)*dvec(i)
  end do
!
end
!
!--------------------------------------------------------------------------
!
! function for vectors: never corrected for intramolecular exceptions
!
subroutine dis_bound2(v1,v2,dis2)
!
  use iounit
  use system
!
  implicit none
!
  RTYPE, INTENT(IN):: v1(3),v2(3)
!
  integer i,k
  RTYPE dvec(3),dis2
!
! distance vector
  dvec(1) = v1(1) - v2(1)
  dvec(2) = v1(2) - v2(2)
  dvec(3) = v1(3) - v2(3)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) + k*bnd_params(i)
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(3))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        dvec(3) = dvec(3) + k*bnd_params(6)
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        dvec(3) = dvec(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound2() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound2() (code # ',bnd_type,').'
    call fexit()
  end if
!
  dis2 = 0.0
  do i=1,3
    dis2 = dis2 + dvec(i)*dvec(i)
  end do
!
end
!
!
!--------------------------------------------------------------------------
!
! function for atoms: never corrected for intramolecular exceptions, also returns dvec
!
subroutine dis_bound3(ii,kk,dvec,dis2)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none 
!
  integer, INTENT(IN):: ii,kk
!
  integer i,k
  RTYPE dvec(3),dis2
!
! distance vector
  dvec(1) = x(kk) - x(ii)
  dvec(2) = y(kk) - y(ii)
  dvec(3) = z(kk) - z(ii)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) + k*bnd_params(i)
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        dvec(3) = dvec(3) + k*bnd_params(6)
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        dvec(3) = dvec(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound3() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound3() (code # ',bnd_type,').'
    call fexit()
  end if
!
  dis2 = 0.0
  do i=1,3
    dis2 = dis2 + dvec(i)*dvec(i)
  end do
!
end
!
!--------------------------------------------------------------------------
!
! the same for reference coordinate set
!
subroutine dis_bound3r(ii,kk,dvec,dis2)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none 
!
  integer, INTENT(IN):: ii,kk
!
  integer i,k
  RTYPE dvec(3),dis2
!
! distance vector
  dvec(1) = xref(kk) - xref(ii)
  dvec(2) = yref(kk) - yref(ii)
  dvec(3) = zref(kk) - zref(ii)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) + k*bnd_params(i)
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        dvec(3) = dvec(3) + k*bnd_params(6)
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        dvec(3) = dvec(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound3r() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound3r() (code # ',bnd_type,').'
    call fexit()
  end if
!
  dis2 = 0.0
  do i=1,3
    dis2 = dis2 + dvec(i)*dvec(i)
  end do
!
end
!
!--------------------------------------------------------------------------
!
! same as dis_bound3 but also returns shift vector instead
!
subroutine dis_bound3t(ii,kk,tvec,dis2)
!
  use iounit
  use atoms
  use system
  use sequen
  use cutoffs
!
  implicit none 
!
  integer, INTENT(IN):: ii,kk
!
  integer i,k
  RTYPE tvec(3),dis2,dvec(3)
!
! distance vector
  dvec(1) = x(kk) - x(ii)
  dvec(2) = y(kk) - y(ii)
  dvec(3) = z(kk) - z(ii)
  tvec(:) = 0.0
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          tvec(i) = k*bnd_params(i)
        else if (dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          tvec(i) = -k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        tvec(3) = k*bnd_params(6)
      else if (dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        tvec(3) = -k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in dis_bound3t() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in dis_bound3t() (code # ',bnd_type,').'
    call fexit()
  end if
!
  dvec(:) = dvec(:) + tvec(:)
  dis2 = dot_product(dvec,dvec)
!
end
!
!--------------------------------------------------------------------------
!
! a simple helper routine to be called after dis_bound_rs or dis_bound_rs2
!
subroutine dis_bound4(ii,kk,svec,dis2)
!
  use atoms
!
  implicit none
!
  integer, INTENT(IN):: ii,kk
!
  RTYPE svec(3),dis2
!
  dis2 = ((x(kk) - x(ii) + svec(1))*(x(kk) - x(ii) + svec(1)))&
 &     + ((y(kk) - y(ii) + svec(2))*(y(kk) - y(ii) + svec(2)))&
 &     + ((z(kk) - z(ii) + svec(3))*(z(kk) - z(ii) + svec(3)))
!
end

!--------------------------------------------------------------------------
!
! another simple helper routine to be called after dis_bound_rs or dis_bound_rs2
!
subroutine dis_bound5(ii,kk,svec,dis2,dvec)
!
  use atoms
!
  implicit none
!
  integer, INTENT(IN):: ii,kk
!
  RTYPE svec(3),dvec(3),dis2
!
  dvec(1) = (x(kk) - x(ii) + svec(1))
  dvec(2) = (y(kk) - y(ii) + svec(2))
  dvec(3) = (z(kk) - z(ii) + svec(3))
  dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3)
!
end
!
!--------------------------------------------------------------------------
!
! a routine to bring molecule jmol into the right reference frame as
! molecule imol (periodic image), but only via comn
!
subroutine shift_bound(imol,jmol,comn)
!
  use iounit
  use system
  use molecule
!
  implicit none
!
  integer i,k,imol,jmol
  RTYPE comn(3),dvec(3)
!
! distance vector
  dvec(1) = com(jmol,1) - com(imol,1)
  dvec(2) = com(jmol,2) - com(imol,2)
  dvec(3) = com(jmol,3) - com(imol,3)
  comn(:) = com(jmol,:)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          comn(i) = com(jmol,i) + k*bnd_params(i)
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) - k*bnd_params(i)
          comn(i) = com(jmol,i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        comn(3) = com(jmol,3) + k*bnd_params(6)
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        dvec(3) = dvec(3) - k*bnd_params(6)
        comn(3) = com(jmol,3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in shift_bound() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in shift_bound() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a routine to compute the PBC-corrected distance vector between
! two molecules c.o.m.es
!
subroutine shift_bound2(imol,jmol,dvec)
!
  use iounit
  use system
  use molecule
!
  implicit none
!
  integer i,k,imol,jmol
  RTYPE dvec(3)
!
! distance vector
  dvec(1) = comm(jmol,1) - comm(imol,1)
  dvec(2) = comm(jmol,2) - comm(imol,2)
  dvec(3) = comm(jmol,3) - comm(imol,3)
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) + k*bnd_params(i)
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          dvec(i) = dvec(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        dvec(3) = dvec(3) + k*bnd_params(6)
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        dvec(3) = dvec(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in shift_bound2() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in shift_bound2() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!--------------------------------------------------------------------------
!
! a routine to compute the PBC-correction to the distance vector between
! two molecules' geometric centers
!
subroutine shift_bound3(imol,jmol,tvec)
!
  use iounit
  use system
  use molecule
!
  implicit none
!
  integer i,k,imol,jmol
  RTYPE dvec(3),tvec(3)
!
! distance vector
  dvec(1) = com(jmol,1) - com(imol,1)
  dvec(2) = com(jmol,2) - com(imol,2)
  dvec(3) = com(jmol,3) - com(imol,3)
  tvec(:) = 0.0
!
! periodic BC (PBC)
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        if (dvec(i).lt.(-0.5*bnd_params(i))) then
          k = floor((-0.5*bnd_params(i)-dvec(i))/bnd_params(i)) + 1
          tvec(i) = tvec(i) + k*bnd_params(i)
        else if(dvec(i).gt.0.5*bnd_params(i)) then
          k = floor((dvec(i)-0.5*bnd_params(i))/bnd_params(i)) + 1
          tvec(i) = tvec(i) - k*bnd_params(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      if (dvec(3).lt.(-0.5*bnd_params(6))) then
        k = floor((-0.5*bnd_params(6)-dvec(3))/bnd_params(6)) + 1
        tvec(3) = tvec(3) + k*bnd_params(6)
      else if(dvec(3).gt.0.5*bnd_params(6)) then
        k = floor((dvec(3)-0.5*bnd_params(6))/bnd_params(6)) + 1
        tvec(3) = tvec(3) - k*bnd_params(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in shift_bound3() (code # ',bnd_shape,').'
      call fexit()
    end if
! soft- or hard wall BC (S/HWBC)
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition in shift_bound3() (code # ',bnd_type,').'
    call fexit()
  end if
!
end
!
!------------------------------------------------------
!
subroutine pt_obox(vec)
!
  use system
  use iounit
!
  implicit none
!
  integer j
  RTYPE vec(3)
!
  if (bnd_shape.eq.1) then
    do j=1,3
      vec(j) = bnd_params(j+3) + 3.0*bnd_params(j)
    end do
  else if (bnd_shape.eq.2) then
    do j=1,3
      vec(j) = bnd_params(j) + 3.0*bnd_params(4)
    end do        
  else if (bnd_shape.eq.3) then
    vec(1:2) = bnd_params(1:2) + 3.0*bnd_params(4)
    vec(3) = vec(3) + 1.5*bnd_params(6)
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in pt_ob&
 &ox() (code # ',bnd_shape,').'
    call fexit()
  end if
!
end
!
!--------------------------------------------------------------------------------
!
! atom-wise Cartesian re-scaling associated with volume changes
!
subroutine rescale_xyz(scf,update)
!
  use atoms
  use system
  use forces
  use molecule
  use iounit
!
  implicit none
!
  integer imol
  RTYPE scf
  logical update
!
  if (bnd_shape.eq.1) then
    x(1:n) = scf*(x(1:n)-bnd_params(4)) + bnd_params(4)
    y(1:n) = scf*(y(1:n)-bnd_params(5)) + bnd_params(5)
    z(1:n) = scf*(z(1:n)-bnd_params(6)) + bnd_params(6)
  else if (bnd_shape.eq.2) then
    x(1:n) = scf*(x(1:n)-bnd_params(1)) + bnd_params(1)
    y(1:n) = scf*(y(1:n)-bnd_params(2)) + bnd_params(2)
    z(1:n) = scf*(z(1:n)-bnd_params(3)) + bnd_params(3)
  else if (bnd_shape.eq.3) then
    x(1:n) = scf*(x(1:n)-bnd_params(1)) + bnd_params(1)
    y(1:n) = scf*(y(1:n)-bnd_params(2)) + bnd_params(2)
    z(1:n) = scf*(z(1:n)-bnd_params(3)) + bnd_params(3)
  else 
    write(ilog,*) 'Fatal. Encountered unsupported box shape in resca&
 &le_xyz(...) (code # ',bnd_shape,').'
    call fexit()
  end if
!
  if (update.EQV..true.) then
!   we need to recompute internals
    do imol=1,nmol
      call update_rigidm(imol)
      call update_comv(imol)
      call genzmat(imol)
    end do
!   velocity adjustment based on box-scaling is a terrible idea
!    if (bnd_shape.eq.1) then
!      cart_v(1:n,1) = cart_v(1:n,1)+((scf-1.0)*x(1:n) + 
! &                                   (1.0-scf)*bnd_params(4))/dyn_dt
!      cart_v(1:n,2) = cart_v(1:n,2)+((scf-1.0)*y(1:n) + 
! &                                   (1.0-scf)*bnd_params(5))/dyn_dt
!      cart_v(1:n,3) = cart_v(1:n,3)+((scf-1.0)*z(1:n) + 
! &                                   (1.0-scf)*bnd_params(6))/dyn_dt
!    else
!      write(ilog,*) 'Fatal. Encountered unsupported box shape in res
! &cale_xyz(...) (code # ',bnd_shape,').'
!      call fexit()
!    end if
  end if
!
end
!
!-----------------------------------------------------------------------
