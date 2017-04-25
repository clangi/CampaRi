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
!------------------------------------------------------------------------
!
subroutine prt_restart(eee)
!
  use iounit
  use system
  use atoms
  use molecule
  use zmatrix
  use mcsums
  use energies
  use mpistuff
  use forces
  use grandensembles
  use movesets
  use wl
!
  implicit none
!
  RTYPE, INTENT(IN):: eee
!
  integer i,imol,irst,freeunit,t1,t2,intdim,k,aone,azero,mt,iomess
  integer iomess2
  character(60) dumpfile
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  logical exists,ismember,isone,fileerror,atrue
  data isone/.true./
!
  save isone
!
  fileerror = .false.
  atrue = .true.
!
#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  if (isone.EQV..true.) then
    dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_1.rst'
    isone = .false.
  else
    dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'_2.rst'
    isone = .true.
  end if
  call strlims(dumpfile,t1,t2)
  inquire(file=dumpfile(t1:t2),exist=exists)
  if (exists.EQV..true.) then
    irst = freeunit()
    open(unit=irst,file=dumpfile(t1:t2),status='replace',iostat=iomess)
    if (iomess.ne.0) then
      write(ilog,*) 'An error has occurred while replacing the old restart file. NFS/Fortran incompatibility bug?'
      write(ilog,*) 'The unit number is ',irst,' (filename: ',dumpfile(t1:t2),'). Check Fortran error message.'
      open(unit=irst,file='N_'//nod(1:tl)//'_RESTART.tmp',status='replace',iostat=iomess2)
      if (iomess2.ne.0) then
        write(ilog,*) 'Fatal. File system error. Unable to open temporary file. Check status and permissions on file system.'
        call fexit()
      else
        fileerror = .true.
        write(ilog,*) 'Warning. Writing restart file to temporary file due to problem described above.'
      end if
    end if
  else
    irst = freeunit()
    open(unit=irst,file=dumpfile(t1:t2),status='new',iostat=iomess)
    if (iomess.ne.0) then
      write(ilog,*) 'Fatal. An error has occurred while opening the new restart file. NFS/Fortran incompatibility bug?'
      write(ilog,*) 'The unit number is ',irst,' (filename: ',dumpfile(t1:t2),'). Check Fortran error message.'
      call fexit()
    end if
  end if 
#else
  if (isone.EQV..true.) then
    dumpfile = basename(1:bleng)//'_1.rst'
    isone = .false.
  else
    dumpfile = basename(1:bleng)//'_2.rst'
    isone = .true.
  end if
  call strlims(dumpfile,t1,t2)
  inquire(file=dumpfile(t1:t2),exist=exists)
  if (exists.EQV..true.) then
    irst = freeunit()
    open(unit=irst,file=dumpfile(t1:t2),status='replace',iostat=iomess)
    if (iomess.ne.0) then
      write(ilog,*) 'An error has occurred while replacing the old restart file. NFS/Fortran incompatibility bug?'
      write(ilog,*) 'The unit number is ',irst,' (filename: ',dumpfile(t1:t2),'). Check Fortran error message.'
      open(unit=irst,file='RESTART.tmp',status='replace',iostat=iomess2)
      if (iomess2.ne.0) then
        write(ilog,*) 'Fatal. File system error. Unable to open temporary file. Check status and permissions on file system.'
        call fexit()
      else
        fileerror = .true.
        write(ilog,*) 'Warning. Writing restart file to temporary file due to problem described above.'
      end if
    end if
  else
    irst = freeunit()
    open(unit=irst,file=dumpfile(t1:t2),status='new',iostat=iomess)
    if (iomess.ne.0) then
      write(ilog,*) 'Fatal. An error has occurred while opening the new restart file. NFS/Fortran incompatibility bug?'
      write(ilog,*) 'The unit number is ',irst,' (filename: ',dumpfile(t1:t2),'). Check Fortran error message.'
      call fexit()
    end if
  end if
#endif
!
#ifdef DISABLE_FLOAT
   16 format(i1,1x,i20,1x,i20,1x,g20.12)
   17 format(6(1x,i10))
   21 format(1(g20.12,1x))
   22 format(2(g20.12,1x))
   23 format(3(g20.12,1x))
   24 format(2(g20.12,1x),3(i20,1x))
   25 format(1(g20.12,1x),3(i20,1x))
#else
   16 format(i1,1x,i20,1x,i20,1x,g20.12)
   17 format(6(1x,i10))
   21 format(1(g20.12,1x))
   22 format(2(g20.12,1x))
   23 format(3(g20.12,1x))
   24 format(2(g20.12,1x),3(i20,1x))
   25 format(1(g20.12,1x),3(i20,1x))
#endif
   20 format(i20)
   19 format(2(i20,1x))
!
! step number
  write(irst,20) nstep
!
! box size and origin (not really needed yet)
  do i=1,7
    write(irst,21) bnd_params(i)
  end do
!
! first xyz
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      write(irst,23) x(i),y(i),z(i)
    end do
  end do
!
! now internals (should be removed)
  do imol=1,nmol
    do i=atmol(imol,1)+1,atmol(imol,2)
      if (i-atmol(imol,1).eq.1) then
        write(irst,21) blen(i)
      else if (i-atmol(imol,1).eq.2) then
        write(irst,22) blen(i),bang(i)
      else
        write(irst,23) blen(i),bang(i),ztor(i)
      end if
    end do
  end do
!
! if dynamics, we'll need velocities
  if (use_dyn.EQV..true.) then
!   if hybrid, write out current cycle
    if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      if (in_dyncyc.EQV..true.) then
        write(irst,20) 1
        write(irst,19) curcyc_start,curcyc_end
      else
        write(irst,20) 0
        write(irst,19) curcyc_start,curcyc_end
      end if
    end if
!   for internal space sampling we also need the inertial masses belonging
!   to the previous configuration
    if (fycxyz.eq.1) then
      do imol=1,nmol
        if (atmol(imol,1).eq.atmol(imol,2)) then
          intdim = 3
        else if (atmol(imol,1).eq.(atmol(imol,2)-1)) then
          intdim = 5
        else
          intdim = ntormol(moltypid(imol))+6
        end if
        if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
          do k=1,intdim 
            write(irst,22) dc_di(imol)%v(k),dc_di(imol)%im(k)
            write(irst,22) dc_di(imol)%ldp(k,4),dc_di(imol)%ldp(k,5)
          end do
        else if (dyn_integrator_ops(1).le.0) then
          do k=1,intdim 
            write(irst,22) dc_di(imol)%v(k),dc_di(imol)%im(k)
          end do
        else
          do k=1,intdim 
            write(irst,23) dc_di(imol)%v(k),dc_di(imol)%im(k),dc_di(imol)%olddat(k,2)
          end do
        end if
      end do
!   for Cartesian space sampling masses are constant -> simple
    else if (fycxyz.eq.2) then
      do imol=1,nmol
        do i=atmol(imol,1),atmol(imol,2)
          write(irst,23) cart_v(i,1),cart_v(i,2),cart_v(i,3)
          if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
            write(irst,23) cart_ldp(i,1,4),cart_ldp(i,2,4),cart_ldp(i,3,4)
            write(irst,23) cart_ldp(i,1,5),cart_ldp(i,2,5),cart_ldp(i,3,5)
          end if
        end do
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported choice of degree&
 &s of freedom in prt_restart(...). This is an omission bug.'
      call fexit()
    end if
!   in NPT/NPE, in certain cases, we need the velocity of the boundary (particle)
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (bnd_shape.eq.2) then
        if (pstat%flag.eq.1) then
          write(irst,21) bnd_v
        else
          write(ilog,*) 'Fatal. Encountered unsupported manostat for&
 & chosen box and ensemble (manostat-flag: ',pstat%flag,') in prt_re&
 &start(...).'
          call fexit()
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary shape&
 & for chosen ensemble (Flag: ',ens%flag,') in prt_restart(...).'
        call fexit()
      end if
    end if
  end if
!
! in (S)GCMC we need particle existences
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    aone = 1
    azero = 0
    do imol=1,nmol
      mt = moltypid(imol)
      if (ismember(fluctypes,mt).EQV..false.) cycle
!
      if (ismember(ispresent,imol).EQV..true.) then
        write(irst,19) imol,aone
      else
        write(irst,19) imol,azero 
      end if
    end do
  end if
!
! in WL, we need histograms, etc.
  if (mc_acc_crit.eq.3) then
    imol = 0
    if (wld%t1on.EQV..true.) imol = 1
    write(irst,16) imol,wld%stepnum,wld%flevel,wld%fval
    if (wld%dimensionality.eq.1) then
      write(irst,17) wld%nbins,wld%minb,wld%maxb
      do imol=1,wld%nbins
        write(irst,24) wld%bctr(imol),wld%g(imol),wld%gh(imol),wld%ghtot(imol),wld%ghbu(imol)
      end do
    else if (wld%dimensionality.eq.2) then
      write(irst,17) wld%nbins2d(:),wld%minb2d(:),wld%maxb2d(:)
      do imol=1,wld%nbins2d(1)
        write(irst,21) wld%bctr2d1(imol)
      end do
      do imol=1,wld%nbins2d(2)
        write(irst,21) wld%bctr2d2(imol)
      end do
      do imol=1,wld%nbins2d(1)
        do i=1,wld%nbins2d(2)
          write(irst,25) wld%g2d(imol,i),wld%gh2d(imol,i),wld%ghtot2d(imol,i),wld%ghbu2d(imol,i)
        end do
      end do
    end if
  end if
!
  write(irst,21) eee ! potential energy for current conformation
  if (dyn_mode.ne.1) then
    write(irst,21) ens%insK
  end if
!
  close(unit=irst)
!
  if (fileerror.EQV..true.) then
    write(ilog,*) 'Successfully completed writing restart-file to temporary file.'
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_restart()
!
  use iounit
  use system
  use atoms
  use molecule
  use zmatrix
  use mcsums
  use energies
  use mpistuff
  use torsn
  use polyavg
  use forces
  use units
  use math
  use grandensembles
  use movesets
  use wl
!
  implicit none
!
  integer i,imol,irst,freeunit,t1,t2,k,intdim,ttc,azero
  integer mt,dummy
  RTYPE eee,eee2,eeek
  RTYPE force3,random
#ifdef ENABLE_MPI
  integer modstep,tl
  character(3) nod
#endif
  character(60) dumpfile
  logical exists,ismember,atrue
!
  azero = 0
  atrue = .true.
!
#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  dumpfile = 'N_'//nod(1:tl)//'_'//basename(1:bleng)//'.rst'
  call strlims(dumpfile,t1,t2)
  inquire(file=dumpfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
      write(ilog,*) 'Restart run fails because restart-file is missi&
 &ng or corrupt (',dumpfile(t1:t2),'). Fatal exit.'
    call fexit()
  end if
  irst = freeunit()
  open(unit=irst,file=dumpfile(t1:t2),status='old')
#else
  dumpfile = basename(1:bleng)//'.rst'
  call strlims(dumpfile,t1,t2)
  inquire(file=dumpfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Restart run fails because restart-file is missing&
 & or corrupt (',dumpfile(t1:t2),'). Fatal exit.'
    call fexit()
  end if
  irst = freeunit()
  open(unit=irst,file=dumpfile(t1:t2),status='old')
#endif
!
#ifdef DISABLE_FLOAT
   16 format(i1,1x,i20,1x,i20,1x,g20.12)
   17 format(6(1x,i10))
   21 format(1(g20.12,1x))
   22 format(2(g20.12,1x))
   23 format(3(g20.12,1x))
   24 format(2(g20.12,1x),3(i20,1x))
   25 format(1(g20.12,1x),3(i20,1x))
#else
   16 format(i1,1x,i20,1x,i20,1x,g20.12)
   17 format(6(1x,i10))
   21 format(1(g20.12,1x))
   22 format(2(g20.12,1x))
   23 format(3(g20.12,1x))
   24 format(2(g20.12,1x),3(i20,1x))
   25 format(1(g20.12,1x),3(i20,1x))
#endif
   20 format(i20)
   19 format(2(i20,1x))
   18 format(1x,'File: ',g15.7,' Re-computed: ',g15.7)
!
! step number
  read(irst,20) nstep
!
! box size and origin (not really needed yet)
  do i=1,7
    read(irst,21) bnd_params(i)
  end do
  call update_bound(atrue)
!
! first xyz
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      read(irst,23) x(i),y(i),z(i)
    end do
  end do
!
! now internals (should be removed)
  do imol=1,nmol
    do i=atmol(imol,1)+1,atmol(imol,2)
      if (i-atmol(imol,1).eq.1) then
        read(irst,21) blen(i)
      else if (i-atmol(imol,1).eq.2) then
        read(irst,22) blen(i),bang(i)
      else
        read(irst,23) blen(i),bang(i),ztor(i)
      end if
    end do
    call genzmat(imol)
  end do
  call zmatfyc2()
!
! if dynamics, we'll also need velocities on the internals (for now) and the
! inertial masses belonging to the previous configuration (of which there is
! no other account!)
  if ((use_dyn.EQV..true.).AND.(mc_compat_flag.ne.1)) then
!   if hybrid, read-in current cycle
    if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      read(irst,20) dummy
      if (dummy.eq.1) then
        in_dyncyc = .true.
      else
        in_dyncyc = .false.
      end if
      read(irst,19) curcyc_start,curcyc_end
    end if
    if (fycxyz.eq.1) then
      do imol=1,nmol
        if (atmol(imol,1).eq.atmol(imol,2)) then
          intdim = 3
        else if (atmol(imol,1).eq.(atmol(imol,2)-1)) then
          intdim = 5
        else
          intdim = ntormol(moltypid(imol))+6
        end if
        if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
          do k=1,intdim 
            read(irst,22) dc_di(imol)%v(k),dc_di(imol)%im(k)
            read(irst,22) dc_di(imol)%ldp(k,4),dc_di(imol)%ldp(k,5)
          end do
        else if (dyn_integrator_ops(1).le.0) then
          do k=1,intdim 
            read(irst,22) dc_di(imol)%v(k),dc_di(imol)%im(k)
          end do
        else
          dc_di(imol)%olddat(:,2) = -1.0
          do k=1,intdim
            read(irst,23) dc_di(imol)%v(k),dc_di(imol)%im(k),dc_di(imol)%olddat(k,2)
            if (dc_di(imol)%olddat(k,2).le.0.0) dc_di(imol)%olddat(k,2) = dc_di(imol)%im(k)
          end do
        end if
      end do
    else if (fycxyz.eq.2) then
      do imol=1,nmol
        do i=atmol(imol,1),atmol(imol,2)
          read(irst,23) cart_v(i,1),cart_v(i,2),cart_v(i,3)
          if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
            read(irst,23) cart_ldp(i,1,4),cart_ldp(i,2,4),cart_ldp(i,3,4)
            read(irst,23) cart_ldp(i,1,5),cart_ldp(i,2,5),cart_ldp(i,3,5)
          end if
        end do
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported choice of degree&
 &s of freedom in read_restart(). This is an omission bug.'
      call fexit()
    end if
!   in NPT/NPE, in certain cases, we need the velocity of the boundary (particle)
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (bnd_shape.eq.2) then
        if (pstat%flag.eq.1) then
          read(irst,21) bnd_v
        else
          write(ilog,*) 'Fatal. Encountered unsupported manostat for&
 & chosen box and ensemble (manostat-flag: ',pstat%flag,') in read_r&
 &estart().'
          call fexit()
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary shape&
 & for chosen ensemble (Flag: ',ens%flag,') in read_restart().'
        call fexit()
      end if
    end if
    do imol=1,nmol
      call update_rigidm(imol)
      if ((dyn_mode.ne.1).AND.(fycxyz.eq.2)) call update_comv(imol)
    end do
  else if ((dyn_mode.ne.1).AND.(mc_compat_flag.eq.1)) then
!   for support of restarts from MC equilibration runs
    do imol=1,nmol
      call update_rigidm(imol)
      if ((dyn_mode.ne.1).AND.(fycxyz.eq.2)) call update_comv(imol)
    end do
  end if
!
! in (S)GCMC we need particle existences
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    do imol=1,nmol
      mt = moltypid(imol)
      if (ismember(fluctypes,mt).EQV..false.) cycle
      read(irst,19) ttc,i
      if (ttc.ne.imol) then
        write(ilog,*) 'Fatal. Set up for (S)GCMC is inconsistent wit&
 &h restart file. This implies either that it was modified, that the&
 & (S)GCMC input file was modified, that the program version changed&
 &, or simply a bug.'
        call fexit()
      end if
      if ((i.eq.1).AND.(ismember(ispresent,imol).EQV..false.)) then
        call insertmol(imol)
      else if ((i.eq.0).AND.(ismember(ispresent,imol).EQV..true.)) then
        call deletemol(imol)
      end if
    end do
  end if
!
! in WL, we need histograms, etc.
! note that - contrary to other variables of this type - we do read in
! wld%stepnum since it is required in the WL switching
  if (mc_acc_crit.eq.3) then
    read(irst,16) k,wld%stepnum,wld%flevel,wld%fval
    wld%stepnumbuf = wld%stepnum
    if (k.eq.1) wld%t1on = .true.
    if (wld%dimensionality.eq.1) then
      read(irst,17) wld%nbins,wld%minb,wld%maxb
      i = wld%nbins
      call wl_realloc_histos(i,imol)
      do imol=1,wld%nbins
        read(irst,24) wld%bctr(imol),wld%g(imol),wld%gh(imol),wld%ghtot(imol),wld%ghbu(imol)
      end do
      wld%g_binsz = (wld%bctr(wld%nbins)-wld%bctr(1))/dble(wld%nbins-1)
      wld%g_max = wld%bctr(wld%nbins)
      wld%g_min = wld%g_max - (wld%nbins-1)*wld%g_binsz
      wld%minv = wld%g_min - 0.5*wld%g_binsz
    else if (wld%dimensionality.eq.2) then
      read(irst,17) wld%nbins2d(:),wld%minb2d(:),wld%maxb2d(:)
      i = wld%nbins2d(1)
      k = wld%nbins2d(2)
      call wl_realloc_histos(i,k)
      do imol=1,wld%nbins2d(1)
        read(irst,21) wld%bctr2d1(imol)
      end do
      do imol=1,wld%nbins2d(2)
        read(irst,21) wld%bctr2d2(imol)
      end do
      do imol=1,wld%nbins2d(1)
        do i=1,wld%nbins2d(2)
          read(irst,25) wld%g2d(imol,i),wld%gh2d(imol,i),wld%ghtot2d(imol,i),wld%ghbu2d(imol,i)
        end do
      end do
      wld%g_binsz2d(1) = (wld%bctr2d1(wld%nbins2d(1))-wld%bctr2d1(1))/dble(wld%nbins2d(1)-1)
      wld%g_binsz2d(2) = (wld%bctr2d2(wld%nbins2d(2))-wld%bctr2d2(1))/dble(wld%nbins2d(2)-1)
      wld%g_max2d(1) = wld%bctr2d1(wld%nbins2d(1))
      wld%g_max2d(2) = wld%bctr2d2(wld%nbins2d(2))
      wld%g_min2d(:) = wld%g_max2d(:) - (wld%nbins2d(:)-1)*wld%g_binsz2d(:)
      wld%minv2d(:) = wld%g_min2d(:) - 0.5*wld%g_binsz2d(:)
    end if
  end if
!
!
! this should conclude all relevant recovery (i.e., the recovery
! of all quantities allowed to change throughout the simulation)
!
! now energy to check sanity
  read(irst,21) eee2
  if (dyn_mode.eq.1) then
    call energy(esterms,eee,azero)
  else
    eee = force3(esterms,esterms_tr,esterms_lr,atrue)
  end if
  if (abs(eee-eee2).gt.0.001) then
    write(ilog,*) 'WARNING. Energies do not match after reading in r&
 &estart-file. This implies either that the restart-file was modifie&
 &d, that the key-file was modified, that the program version change&
 &d, or simply a bug.'
    write(ilog,18) eee2,eee
  end if
  esave = eee
  if ((dyn_mode.ne.1).AND.(mc_compat_flag.ne.1)) then
    read(irst,21) eeek
    if ((in_dyncyc.EQV..false.).AND.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
!     do nothing since dynamics variables during MC are meaningless
    else
      call ensv_helper()
      if (abs(eeek-ens%insK).gt.0.001) then
        write(ilog,*) 'WARNING. Kinetic energies do not match after re&
 &ading in restart-file. This implies either that it was modified, t&
 &hat the key-file was modified, that the program version changed, o&
 &r simply a bug.'
        write(ilog,18) eeek,ens%insK
      end if
    end if
!   to support starting MD runs from MC equilibration, here we'll randomize velocities
  else if (((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.4).OR.(dyn_mode.eq.6)).AND.&
          &(mc_compat_flag.eq.1)) then
    if (fycxyz.eq.1) then
      call cart2int_I(azero)
      do i=1,nmol
        dc_di(i)%im(:) = dc_di(i)%olddat(:,2)
      end do
      call randomize_velocities(azero)
      ens%insT = kelvin
      if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
        tstat%grpT(:) = kelvin
      else if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
        call init_ldps()
      end if
    else if (fycxyz.eq.2) then
      call randomize_cart_velocities(azero)
      ens%insT = kelvin
      if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
        tstat%grpT(:) = kelvin
      else if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
        call init_cart_ldps()
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported choice of degree&
 &s of freedom in read_restart(). This is an omission bug.'
      call fexit()
    end if
!   always start a fake-restarted hybrid run with an MC portion
  else if (((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(mc_compat_flag.eq.1)) then
    in_dyncyc = .false.
    curcyc_start = nstep + 1
    curcyc_end = nstep + min_mccyclen + floor(random()*(max_mccyclen-min_mccyclen) + 0.5)
    curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
!     important: re_freq is never less than 2
      if (mod(curcyc_start,re_freq).eq.0) then
        curcyc_start = curcyc_start + 1
      end if
      if (mod(curcyc_end,re_freq).eq.0) then
        curcyc_end = curcyc_end + 1
      end if
!     publish the new cycle to all nodes and receive on slaves (note that this will override
!     the cycle settings just obtained on the slave node)
      call MPI_SyncHybridCycle(nstep)
    end if
#endif
  end if
!
! if in MPI, collect where every node is and reduce to the slowest one
! (yes, somewhat inaccurate, but we have no way to recover the buffered
! messages, and hence the calculation would be out of sync and eventually
! crash)
#ifdef ENABLE_MPI
  call MPI_Synchronize()
#endif
!
! fix broken hybrid cycle assignments and resync if needed
  if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    if (curcyc_start.le.nstep) then
      curcyc_start = nstep + 1
    else if ((curcyc_start.gt.nstep).AND.(curcyc_end.ge.curcyc_start)) then ! preserve length of stretch in restart file
      curcyc_end = curcyc_end - (curcyc_start - (nstep+1))
      curcyc_start = nstep + 1
    end if
    if (curcyc_end.lt.curcyc_start) then
      if (in_dyncyc.EQV..false.) then
        curcyc_end = curcyc_start - 1 + min_mccyclen + floor(random()*(max_mccyclen-min_mccyclen) + 0.5)
      else
        curcyc_end = curcyc_start - 1 + min_dyncyclen + floor(random()*(max_dyncyclen-min_dyncyclen) + 0.5)
      end if
      curcyc_end = min(curcyc_end,nsim)
    end if
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
!     important: re_freq is never less than 2
      if (mod(curcyc_start,re_freq).eq.0) then
        curcyc_start = curcyc_start + 1
      end if
      if (mod(curcyc_end,re_freq).eq.0) then
        curcyc_end = curcyc_end + 1
      end if
!     publish the new cycle to all nodes and receive on slaves (note that this will override
!     the cycle settings just obtained on the slave node)
      call MPI_SyncHybridCycle(nstep)
    end if
#endif
  end if
!
! fix certain necessary counters
  do i=1,nstep
#ifdef ENABLE_MPI
    if (use_MPIAVG.EQV..true.) then
      modstep = mod(i,enout)
      if (modstep.eq.0) then
        if (mpi_cnt_en.eq.(mpi_nodes-1)) then
          mpi_cnt_en = 0
        else
          mpi_cnt_en = mpi_cnt_en + 1
        end if
      end if
      modstep = mod(i,torout)
      if (modstep.eq.0) then
        if (mpi_cnt_tor.eq.(mpi_nodes-1)) then
          mpi_cnt_tor = 0
        else
          mpi_cnt_tor = mpi_cnt_tor + 1
        end if
      end if
      modstep = mod(i,savcalc)
      if ((i.gt.nequil).AND.(modstep.eq.0)) then
        if (mpi_cnt_sav.eq.(mpi_nodes-1)) then
          mpi_cnt_sav = 0
        else
          mpi_cnt_sav = mpi_cnt_sav + 1
        end if
      end if
      modstep = mod(i,polout)
      if (modstep.eq.0) then
        if (mpi_cnt_pol.eq.(mpi_nodes-1)) then
          mpi_cnt_pol = 0
        else
          mpi_cnt_pol = mpi_cnt_pol + 1
        end if
      end if
      modstep = mod(i,xyzout)
      if ((i.gt.nequil).AND.(modstep.eq.0)) then
        if (mpi_cnt_xyz.eq.(mpi_nodes-1)) then
          mpi_cnt_xyz = 0
        else
          mpi_cnt_xyz = mpi_cnt_xyz + 1
        end if
      end if
      modstep = mod(i,covcalc)
      if ((i.gt.nequil).AND.(modstep.eq.0)) then
       if (mpi_cnt_trcv.eq.(mpi_nodes-1)) then
          mpi_cnt_trcv = 0
        else
          mpi_cnt_trcv = mpi_cnt_trcv + 1
        end if
      end if
    end if
#endif
  end do
!
  close(unit=irst)
!
end
!
!-----------------------------------------------------------------------
! 

