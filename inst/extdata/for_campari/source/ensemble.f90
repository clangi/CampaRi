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
! CONTRIBUTIONS: Jose Pulido                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!     
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
! A collection of subroutines related to computing quantities
! relevant for characterizing the ensemble as a whole
!
!-----------------------------------------------------------------------
!
subroutine ensv_helper()
!
  use system
  use torsn
  use forces
  use atoms
  use units
  use molecule
  use math
!
  implicit none
!
  integer imol,i,j,ttc
  RTYPE ttt,temp,ttg(tstat%n_tgrps)
!
  ttt = 0.0
  ttg(:) = 0.0
!
  if (fycxyz.eq.1) then
    if (dyn_integrator_ops(1).gt.0) then ! we have estimated inertia at time t1.0
      do imol=1,nmol
        do j=1,3
          temp = dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
            temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
            ttt = ttt + temp
            ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          end do    
        end if
        do j=ttc+1,ttc+ntormol(moltypid(imol))
          temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    else
      do imol=1,nmol
        do j=1,3
          temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
            temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
            ttt = ttt + temp
            ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          end do    
        end if
        do j=ttc+1,ttc+ntormol(moltypid(imol))
          temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    end if
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        temp = pstat%params(1)*bnd_v*bnd_v
        ttt = ttt + temp
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      else if (pstat%flag.eq.3) then
!       do nothing (extended ensemble particles not considered)
      end if
    end if
    tstat%grpT(:) = ttg(:) / (tstat%grpdof(:)*u_dyn_kb)
    ens%insR(8) = ens%insK
    ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
    ens%insT = ttt / ((totrbd+ndyntorsn+ens%boxdof-n_constraints)*u_dyn_kb)
    ens%insR(9) = ens%insK2
    ttt = sum(mass(1:n)*(cart_v(1:n,1)*cart_v(1:n,1) + cart_v(1:n,2)*cart_v(1:n,2) + cart_v(1:n,3)*cart_v(1:n,3)))
    ens%insK2 = (1.0/u_dyn_fconv)*0.5*ttt
!
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
      call virial(ens%insVirT)
      call ekinten(ens%insKinT)
      ens%insVirT(:,:) = 0.0
      ens%insP = 0.0
      ens%insP = (1./3.)*(u_dyn_virconv/ens%insV)*((1.0/u_dyn_fconv)*&
   &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
   &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
    end if
!
  else if (fycxyz.eq.2) then
!
    do imol=1,nmol
      do i=atmol(imol,1),atmol(imol,2)
        do j=1,3
          temp = mass(i)*cart_v(i,j)*cart_v(i,j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        temp = pstat%params(1)*bnd_v*bnd_v
        ttt = ttt + temp
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      else if (pstat%flag.eq.2) then
!       do nothing (Markov chain not considered)
      else if (pstat%flag.eq.3) then
!       do nothing (extended ensemble particles not considered)
      end if
    end if
    tstat%grpT(:) = ttg(:) / (tstat%grpdof(:)*u_dyn_kb)
    ens%insR(8) = ens%insK
    ens%insK = (1.0/u_dyn_fconv)*0.5*ttt
    ens%insT = ttt / (dble(3*n+ens%boxdof-n_constraints)*u_dyn_kb)
!
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
      call virial(ens%insVirT)
      call ekinten(ens%insKinT)
      ens%insVirT(:,:) = 0.0
      ens%insP = 0.0
      ens%insP = (1./3.)*(u_dyn_virconv/ens%insV)*((1.0/u_dyn_fconv)*&
   &        dble(3.0*n+ens%boxdof)/dble(3*n+ens%boxdof-n_constraints)*&
   &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
   &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
    end if
  end if
end
!
!----------------------------------------------------------------------------------
!
subroutine get_ensv(boxovol,tpi)
!
  use iounit
  use system
  use molecule
  use forces
  use units
  use mcsums
  use math
  use torsn
  use energies
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer imol,i,j,ttc
  RTYPE ttt,temp,boxovol,ttg(tstat%n_tgrps)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL,doavg,usebu
  integer sta,sto
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, get_ensv(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
  if (tpi.le.1) then
    ens%insR(9) = ens%insK2
    ens%insK2 = 0.0
  end if
#else
  ens%insR(9) = ens%insK2
#endif
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (nstep.gt.nequil) then
      doavg = .true.
    else
      doavg = .false.
    end if
    if (dyn_integrator_ops(1).gt.0) then
      usebu = .true.
    else
      usebu = .false.
    end if
    call get_Ekin_threads(tpi,doavg,usebu)
    ttt = sum(mass(thr_limits(1,tpi):thr_limits(2,tpi))*(cart_v(thr_limits(1,tpi):thr_limits(2,tpi),1)**2 + &
 &                                                       cart_v(thr_limits(1,tpi):thr_limits(2,tpi),2)**2 + &
 &                                                       cart_v(thr_limits(1,tpi):thr_limits(2,tpi),3)**2))
    thr_rutil(20,tpi) = ttt
!$OMP BARRIER
!$OMP SINGLE
    ens%insK2 = sum(thr_rutil(20,1:thrdat%maxn))
!$OMP END SINGLE
  else
!$OMP SINGLE
#endif
  ttt = 0.0
  ttg(:) = 0.0
  if (dyn_integrator_ops(1).gt.0) then ! we have estimated inertia at time t1.0
    do imol=1,nmol
      do j=1,3
        temp = dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
!       WARNING: support missing
      else
        ttc = 6
        do j=4,6
          temp = (1.0/(RADIAN*RADIAN))*&
 &             dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (nstep.gt.nequil) then
            dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
          end if
        end do    
      end if
      do j=ttc+1,ttc+ntormol(moltypid(imol))
        temp = (1.0/(RADIAN*RADIAN))*&
 &            dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
    end do
  else
    do imol=1,nmol
      do j=1,3
        temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
!       WARNING: support missing
      else
        ttc = 6
        do j=4,6
          temp = (1.0/(RADIAN*RADIAN))*&
 &             dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (nstep.gt.nequil) then
            dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
          end if
        end do    
      end if
      do j=ttc+1,ttc+ntormol(moltypid(imol))
        temp = (1.0/(RADIAN*RADIAN))*&
 &            dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        if (nstep.gt.nequil) then
          dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp/u_dyn_kb
        end if
      end do
    end do
  end if
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      if (nstep.gt.nequil) then
        bnd_avgT = bnd_avgT + temp/u_dyn_kb
      end if
    else if (pstat%flag.eq.3) then
!     do nothing (extended ensemble particles not considered)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in get_ensv(...). T&
 &his is most likely an omission bug.'
      call fexit()
    end if
  end if
  ens%insT = ttt
  tstat%grpT(1:tstat%n_tgrps) = ttg(1:tstat%n_tgrps)
  ens%insK2 = sum(mass(1:n)*(cart_v(1:n,1)*cart_v(1:n,1) + cart_v(1:n,2)*cart_v(1:n,2) + cart_v(1:n,3)*cart_v(1:n,3)))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
  end if 
#endif
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  ens%insK2 = (1.0/u_dyn_fconv)*0.5*ens%insK2
  tstat%grpT(:) = tstat%grpT(:) / (tstat%grpdof(:)*u_dyn_kb)
  ens%insR(8) = ens%insK
  ens%insK = (1.0/u_dyn_fconv)*0.5*ens%insT
  ens%insT = ens%insT / ((totrbd+ndyntorsn+ens%boxdof-n_constraints)*u_dyn_kb)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
    call virial(ens%insVirT)
    call ekinten(ens%insKinT)
    ens%insVirT(:,:) = 0.0
    ens%insP = 0.0
    ens%insP = (1./3.)*(u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*&
 &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
 &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
  end if
!  
  if (nstep.gt.nequil) then
    if (nstep.eq.(nequil+1)) then
      if (ens%flag.eq.2) then
        ens%insR(2) = ens%insK + ens%insU
        ens%insR(3) = ens%insK2 + ens%insU
      end if
    end if
    ens%avgcnt = ens%avgcnt + 1
    ens%avgK = ens%avgK + ens%insK
    ens%avgK2 = ens%avgK2 + ens%insK2
    ens%avgU = ens%avgU + ens%insU
    ens%avgT = ens%avgT + ens%insT
    ens%avgR(4) = ens%avgR(4) + abs(ens%insK2-ens%insK)
    ens%avgR(1) = ens%avgR(1) + ens%insR(1)
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!      ens%avgV = ens%avgV + (nmol+1)*ens%insT*gasconst/
! &                              (boxovol*ens%insP/u_dyn_virconv)
      ens%avgV = ens%avgV + boxovol/1.0e6
      do i=1,3
        do j=1,3
          ens%avgVirT(i,j) = ens%avgVirT(i,j) +(u_dyn_virconv/boxovol)*ens%insVirT(i,j)
          ens%avgPT(i,j) = ens%avgPT(i,j) + (u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*ens%insKinT(i,j) + ens%insVirT(i,j))
        end do
      end do
      ens%avgP = ens%avgP + ens%insP
    else
      ens%avgR(2) = ens%avgR(2) + (ens%insU + ens%insK)**2
      ens%avgR(3) = ens%avgR(3) + (ens%insU + ens%insK2)**2
      ens%avgR(5) = ens%avgR(5) + ens%insU**2
      ens%avgR(6) = ens%avgR(6) + ens%insK**2
      ens%avgR(7) = ens%avgR(7) + ens%insK2**2
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!-----------------------------------------------------------------------
!
! same as get_ensv for straight Cartesian dynamics
!
subroutine get_cart_ensv(boxovol,tpi)
!
  use iounit
  use system
  use molecule
  use forces
  use units
  use mcsums
  use atoms
  use energies
#ifdef ENABLE_THREADS
  use threads
  use sequen
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer imol,i,j
  RTYPE ttt,temp,boxovol,ttg(tstat%n_tgrps)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL,doavg,afalse
  integer sta,sto
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, get_cart_ensv(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
#endif
!
  ttg(:) = 0.0
  ttt = 0.0
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (nstep.gt.nequil) then
      doavg = .true.
    else
      doavg = .false.
    end if
    afalse = .false.
    call get_Ekin_threads(tpi,doavg,afalse)
  else
!$OMP SINGLE
#endif
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      do j=1,3
        temp = mass(i)*cart_v(i,j)*cart_v(i,j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
!        if (nstep.gt.nequil) then
!          cart_avgT(i,j) = cart_avgT(i,j) + temp/u_dyn_kb
!        end if
      end do
    end do
  end do
! currently not supported
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      tstat%grpT(tstat%molgrp(nmol+1)) = tstat%grpT(tstat%molgrp(nmol+1)) + temp
      if (nstep.gt.nequil) then
        bnd_avgT = bnd_avgT + temp/u_dyn_kb
      end if
    else if (pstat%flag.eq.2) then
!     do nothing (Markov chain not considered)
    else if (pstat%flag.eq.3) then
!     do nothing (extended ensemble particles not considered)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in get_cart_ensv(..&
 &.). This is most likely an omission bug.'
      call fexit()
    end if
  end if
  ens%insT = ttt
  tstat%grpT(1:tstat%n_tgrps) = ttg(1:tstat%n_tgrps)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
  end if
!
!$OMP SINGLE
#endif
!
  tstat%grpT(:) =  tstat%grpT(:) / (tstat%grpdof(:)*u_dyn_kb)
  ens%insR(8) = ens%insK
  ens%insK = (1.0/u_dyn_fconv)*0.5*ens%insT
  ens%insT = ens%insT / (dble(3*n+ens%boxdof-n_constraints)*u_dyn_kb)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then ! nonfunctional
    call virial(ens%insVirT)
    call ekinten(ens%insKinT)
    ens%insVirT(:,:) = 0.0
    ens%insP = 0.0
    ens%insP = (1./3.)*(u_dyn_virconv/boxovol)*((1.0/u_dyn_fconv)*&
 &        dble(3.0*n+ens%boxdof)/dble(3*n+ens%boxdof-n_constraints)*&
 &            (ens%insKinT(1,1)+ens%insKinT(2,2)+ens%insKinT(3,3)) +&
 &            (ens%insVirT(1,1)+ens%insVirT(2,2)+ens%insVirT(3,3)) )
  end if
!  
  if (nstep.gt.nequil) then
    if (nstep.eq.(nequil+1)) then
      if (ens%flag.eq.2) then
        ens%insR(2) = ens%insK + ens%insU
      end if
    end if
    ens%avgcnt = ens%avgcnt + 1
    ens%avgK = ens%avgK + ens%insK
    ens%avgU = ens%avgU + ens%insU
    ens%avgT = ens%avgT + ens%insT
    ens%avgR(1) = ens%avgR(1) + ens%insR(1)
    if (ens%insR(5).gt.0.0) then
      ens%avgR(8) = ens%avgR(8) + ens%insR(6)/ens%insR(5)
      ens%avgR(9) = ens%avgR(9) + ens%insR(7)/ens%insR(5)
    end if
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!      ens%avgV = ens%avgV + (nmol+1)*ens%insT*gasconst/
! &                              (boxovol*ens%insP/u_dyn_virconv)
      ens%avgV = ens%avgV + boxovol/1.0e6
      do i=1,3
        do j=1,3
       ens%avgVirT(i,j) = ens%avgVirT(i,j) +(u_dyn_virconv/boxovol)*ens%insVirT(i,j)
!        ens%avgPT(i,j) = ens%avgPT(i,j) + (u_dyn_virconv/boxovol)*(
! &                            (1.0/u_dyn_fconv)*ens%insKinT(i,j) +
! &                                              ens%insVirT(i,j) )
        end do
      end do
      ens%avgP = ens%avgP + ens%insP
    else
      ens%avgR(2) = ens%avgR(2) + (ens%insU + ens%insK)**2
      ens%avgR(5) = ens%avgR(5) + ens%insU**2
      ens%avgR(6) = ens%avgR(6) + ens%insK**2
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!----------------------------------------------------------------------------------
!
  subroutine prt_curens(istep,boxovol,fr1,fr2)
!
  use system
  use iounit
  use mcsums
  use cutoffs
  use movesets
  use units
!
  implicit none
!
  RTYPE, INTENT(IN):: fr1,fr2,boxovol
  integer, INTENT(IN):: istep
!
 555 format(1000(g14.6,1x))

  if ((fycxyz.ne.1).AND.(fycxyz.ne.2)) then
    write(ilog,*) 'Fatal. Unsupported choice of degrees of freedom in prt_curens(...).&
  & This is definitely an omission bug.'
    call fexit()
  end if
!
 997  format(' Kinetic E   | Potential E | Total E     | Temperature | Drift-xyz   | Drift-Euler |')
 998  format(' Kinetic E   | Potential E | Enthalpy    | Temperature | Pressure    | Box volume  |')
 995  format(' Kinetic E   | Cart. K.E.  | Potential E | Total E     | Temperature | Drift-xyz   | Drift-Euler |')
 996  format(' Kinetic E   | Cart. K.E.  | Potential E | Enthalpy    | Temperature | Pressure    | Box volume  |')

 999  format('---------------------------------------------&
 &---------------------------------------')
  if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
!   write header line
    if ((istep.eq.curcyc_start).AND.(nsim.ge.nsancheck).AND.&
   &(floor(1.0*curcyc_end/nsancheck).gt.floor(1.0*istep/nsancheck))) then
      if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
        if (fycxyz.eq.1) write(ilog,995)
        if (fycxyz.eq.2) write(ilog,997)
      else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
        if (fycxyz.eq.1) write(ilog,996)
        if (fycxyz.eq.2) write(ilog,998)
      end if
    end if
  else if (istep.eq.1) then
!   write header line
    if (nsim.ge.nsancheck) then
      if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
        if (fycxyz.eq.1) write(ilog,995)
        if (fycxyz.eq.2) write(ilog,997)
      else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
        if (fycxyz.eq.1) write(ilog,996)
        if (fycxyz.eq.2) write(ilog,998)
      end if
      write(ilog,999)
    end if
  end if
  if (mod(nstep,nsancheck).eq.0) then
    if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
      if (fycxyz.eq.1) write(ilog,56) ens%insK,ens%insK2,ens%insU,ens%insK + ens%insU,ens%insT,fr1,fr2
      if (fycxyz.eq.2) write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU,ens%insT,fr1,fr2
    else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (fycxyz.eq.1) write(ilog,56) ens%insK,ens%insK2,ens%insU,ens%insK + ens%insU + extpress*ens%insV/u_dyn_virconv,&
 &                                    ens%insT,ens%insP,0.001*boxovol
      if (fycxyz.eq.2) write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU + extpress*ens%insV/u_dyn_virconv,&
 &                                    ens%insT,ens%insP,0.001*boxovol
      write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU + extpress*ens%insV/u_dyn_virconv,ens%insT,ens%insP,boxovol/1.0e3
     end if
  end if
!
   56 format(50(g13.6,1x))
!
end
!
!--------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine get_Ekin_threads(tpi,doavg,usebu)
!
  use iounit
  use threads
  use forces
  use system
  use molecule
  use math
  use sequen
  use atoms
  use units
  use mcsums
!
  implicit none
!
  integer, INTENT(IN):: tpi
  logical, INTENT(IN):: doavg,usebu
!
  integer i,j,ttc,imol,mollo,molhi
  RTYPE ttt,ttg(tstat%n_tgrps),temp,ihlp
!
  if (tpi.le.0) then
    write(ilog,*) 'Fatal. Subroutine get_Ekin_threads(...) must be called exclusively from within an OpenMP parallel &
 &region with correct thread identifiers. This is a bug.'
    call fexit()
  end if
  ttt = 0.
  ttg(:) = 0.0
!$OMP SINGLE
  ens%insT = 0.
  tstat%grpT(1:tstat%n_tgrps) = 0.0
!$OMP END SINGLE NOWAIT
!
  if (fycxyz.eq.1) then
    ihlp = 1.0/u_dyn_kb
    if (usebu.EQV..true.) then ! we have estimated inertia at time t1.0
      do imol=thr_limits(61,tpi),thr_limits(62,tpi)
        do j=1,3
          temp = dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (doavg.EQV..true.) dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp*ihlp
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5 
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
            temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
            ttt = ttt + temp
            ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
            if (doavg.EQV..true.) dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp*ihlp
          end do    
        end if
        do j=ttc+1,ttc+ntormol(moltypid(imol))
          temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%olddat(j,2)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (doavg.EQV..true.) dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp*ihlp
        end do
      end do
    else
      do imol=thr_limits(61,tpi),thr_limits(62,tpi)
        do j=1,3
          temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (doavg.EQV..true.) dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp*ihlp
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5 
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
            temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
            ttt = ttt + temp
            ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
            if (doavg.EQV..true.) dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp*ihlp
          end do    
        end if
        do j=ttc+1,ttc+ntormol(moltypid(imol))
          temp = (1.0/(RADIAN*RADIAN))*dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
          if (doavg.EQV..true.) dc_di(imol)%avgT(j) = dc_di(imol)%avgT(j) + temp*ihlp
        end do
      end do
    end if
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!$OMP SINGLE
      if (pstat%flag.eq.1) then
        temp = pstat%params(1)*bnd_v*bnd_v
        ttt = ttt + temp
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
        if (doavg.EQV..true.) bnd_avgT = bnd_avgT + temp*ihlp
      else if (pstat%flag.eq.3) then
!       do nothing
      else
        write(ilog,*) 'Fatal. Unsupported manostat in get_Ekin_threads(...). This is most likely an omission bug.'
        call fexit()
      end if
!$OMP END SINGLE
    end if
  else if (fycxyz.eq.2) then
    if (thr_limits(2,tpi).ge.thr_limits(1,tpi)) then
      mollo = molofrs(atmres(thr_limits(1,tpi)))
      molhi = molofrs(atmres(thr_limits(2,tpi)))
    else
      mollo = 1
      molhi = 0
    end if
    do imol=mollo,molhi
      do i=max(thr_limits(1,tpi),atmol(imol,1)),min(thr_limits(2,tpi),atmol(imol,2))
        do j=1,3
          temp = mass(i)*cart_v(i,j)*cart_v(i,j)
          ttt = ttt + temp
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
        end do
      end do
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!$OMP SINGLE
      if (pstat%flag.eq.1) then
        temp = pstat%params(1)*bnd_v*bnd_v
        ttt = ttt + temp
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
      else
        write(ilog,*) 'Fatal. Unsupported manostat in get_Ekin_threads(...). This is most likely an omission bug.'
        call fexit()
      end if
!$OMP END SINGLE
    end if
  else
    write(ilog,*) 'Fatal. Unsupported choice of degrees of freedom in get_Ekin_threads(...).&
  & This is definitely an omission bug.'
    call fexit()
  end if
!
  thr_rutil(1,tpi) = ttt
  thr_rutil(2:(tstat%n_tgrps+1),tpi) = ttg(1:tstat%n_tgrps)
!$OMP BARRIER
!$OMP SINGLE
  ens%insT = sum(thr_rutil(1,1:thrdat%maxn))
  tstat%grpT(1:tstat%n_tgrps) = sum(thr_rutil(2:(tstat%n_tgrps+1),1:thrdat%maxn),dim=2)
!$OMP END SINGLE
!
end 
!
#endif
!
!---------------------------------------------------------------------------
!
! the routine to initialize dof-velocities to a Boltzmann-distributed
! ensemble using target the target temperature and (current) masses
!
subroutine randomize_velocities(tpi)
!
  use system
  use forces
  use atoms
  use molecule
  use units
  use math
  use torsn
  use iounit
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer imol,j,k,i,intdim,ttc
  RTYPE vind,normal,ttt,fr1,fr2,temp
  RTYPE ttg(tstat%n_tgrps),corrfac(tstat%n_tgrps)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL,afalse
  integer sta,sto
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, randomize_velocities(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(61,tpi)
    sto = thr_limits(62,tpi)
  else
    sta = 1
    sto = nmol
  end if
  if (tpi.le.1) ens%insK = 10.0
#else
  ens%insK = 10.0 ! fake set
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  do imol=1,nmol
!
!   first: center-of-mass translation
!
    do j=1,3
!     velocities in A/ps
      if (dc_di(imol)%frz(j).EQV..true.) then
        dc_di(imol)%v(j) = 0.0
        cycle
      end if
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
      dc_di(imol)%v(j) = vind
    end do
!
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
!
!   second: rigid-body rotation

    if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
!
      ttc = 5
      write(ilog,*) 'Fatal. Diatomic molecules currently not support&
 &ed by dynamics implementation.'
      call fexit()
!     WARNING: support missing
!
    else
!
      ttc = 6
      do j=4,6
!       velocities in deg/ps
        if (dc_di(imol)%frz(j).EQV..true.) then
          dc_di(imol)%v(j) = 0.0
          cycle
        end if
        vind = normal()
        vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
        dc_di(imol)%v(j) = RADIAN*vind
      end do
!
    end if
!
!   third: torsional degrees of freedom
!
    do j=1,ntormol(moltypid(imol))
      if (dc_di(imol)%frz(ttc+j).EQV..true.) then
        dc_di(imol)%v(ttc+j) = 0.0
        cycle
      end if
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(ttc+j))
      dc_di(imol)%v(ttc+j) = RADIAN*vind
    end do
!
  end do
!
! fourth: box degrees of freedom (if applicable)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/pstat%params(1))
      bnd_v = vind
    else if (pstat%flag.eq.3) then
      pstat%params(2) = 0.0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_veloci&
 &ties(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
  call drift_removal(fr1,fr2,tpi)
!
! now re-scale to match exact temperature request
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
  if (tpi.gt.0) then
   afalse = .false.
   call get_Ekin_threads(tpi,afalse,afalse) ! sets preliminary values into ens%insT and tstat%grpT(:)
!   write(*,*) tpi,ens%insT,tstat%grpT(:)
  else
!$OMP SINGLE
#endif
  ttt = 0.0
  ttg(:) = 0.0
  do imol=1,nmol
    do j=1,3
      temp = dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
      ttt = ttt + temp
      ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
    end do
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
      ttc = 3
    else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
      ttc = 5
!     WARNING: support missing
    else
      ttc = 6
      do j=4,6
        temp = (1.0/(RADIAN*RADIAN))*&
 &           dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
      end do    
    end if
    do j=ttc+1,ttc+ntormol(moltypid(imol))
      temp = (1.0/(RADIAN*RADIAN))*&
 &             dc_di(imol)%im(j)*dc_di(imol)%v(j)*dc_di(imol)%v(j)
      ttt = ttt + temp
      ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
    else if (pstat%flag.eq.3) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_veloci&
 &ties(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
  ens%insT = ttt
  tstat%grpT(1:tstat%n_tgrps) = ttg(1:tstat%n_tgrps)
!
#ifdef ENABLE_THREADS
!$OMP END SINGLE
  end if
!$OMP SINGLE
#endif
!
  if ((totrbd+ndyntorsn+ens%boxdof-n_constraints).le.0) then
    write(ilog,*) 'Attempting a calculation with no remaining degree&
 &s of freedom. This is fatal.'
    call fexit()
  end if
  if (ens%insT.le.0.0) then
    write(ilog,*) 'System is initially void of kinetic energy even t&
 &hough there are unconstrained degrees of freedom. This is either a&
 & bug or a numerical precision error due to an extremely small temp&
 &erature request. Fatal exit.'
    call fexit()
  end if
  do i=1,tstat%n_tgrps
    if (tstat%grpT(i).le.0.0) then
      write(ilog,*) 'At least on of the T-coupling groups is initial&
 &ly void of kinetic energy even though there are unconstrained degr&
 &ees of freedom. This is either a bug or a numerical precision erro&
 &r due to an extremely small temperature request. Fatal exit.'
      call fexit()
    end if
  end do
  ens%insK = (1.0/u_dyn_fconv)*0.5*ens%insT
  ens%insT = ens%insT / (dble(totrbd+ndyntorsn+ens%boxdof-n_constraints)*u_dyn_kb)
  tstat%grpT(:) = tstat%grpT(:)/(u_dyn_kb*tstat%grpdof(:))
  corrfac(:) = sqrt(kelvin/tstat%grpT(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE COPYPRIVATE(corrfac)
!
  do imol=sta,sto
#else
  do imol=1,nmol
#endif
    if (atmol(imol,1).eq.atmol(imol,2)) then
      intdim = 3
    else if (atmol(imol,1).eq.(atmol(imol,2)-1)) then
      intdim = 5
    else
      intdim = ntormol(moltypid(imol))+6
    end if
    do k=1,intdim
      dc_di(imol)%v(k) =corrfac(tstat%molgrp(imol))*dc_di(imol)%v(k)
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
!$OMP SINGLE
    if (pstat%flag.eq.1) then
      bnd_v = corrfac(tstat%molgrp(nmol+1))*bnd_v
    end if
!$OMP END SINGLE
  end if
! 
end
!
!---------------------------------------------------------------------------
!
! the same for straight cartesian dynamics
!
subroutine randomize_cart_velocities(tpi)
!
  use system
  use forces
  use atoms
  use molecule
  use units
  use math
  use torsn
  use iounit
#ifdef ENABLE_THREADS
  use threads
  use sequen
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer imol,j,i
  RTYPE vind,ttg(tstat%n_tgrps),normal,ttt,fr1,fr2
  RTYPE temp,corrfac(tstat%n_tgrps)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL,afalse
  integer mollo,molhi,sta2,sto2
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, randomize_velocities(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta2 = thr_limits(1,tpi)
    sto2 = thr_limits(2,tpi)
  else
    sta2 = 1
    sto2 = n
  end if
  if (tpi.le.1) ens%insK = 10.0
#else
  ens%insK = 10.0 ! fake set
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  do imol=1,nmol
!
!   straightforward
!
    do i=atmol(imol,1),atmol(imol,2)
      if (mass(i).le.0.0) cycle
      do j=1,3
        if (cart_frz(i,j).EQV..true.) then
          cart_v(i,j) = 0.0
          cycle
        end if
!       velocities in A/ps
        vind = normal()
        vind = vind*sqrt(u_dyn_kb*kelvin/mass(i))
        cart_v(i,j) = vind
      end do
    end do
!
    call update_comv(imol)
!
  end do
!
! box degrees of freedom (if applicable)
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      vind = normal()
      vind = vind*sqrt(u_dyn_kb*kelvin/pstat%params(1))
      bnd_v = vind
    else if (pstat%flag.eq.2) then
!     do nothing for Markov chain
    else if (pstat%flag.eq.3) then
      pstat%params(2) = 0.0
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_cart&
 &_velocities(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
  call drift_removal(fr1,fr2,tpi)
!
! now re-scale to match exact temperature request
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
  if (tpi.gt.0) then
    afalse = .false.
    call get_Ekin_threads(tpi,afalse,afalse) ! sets preliminary values into ens%insT and tstat%grpT(:)
!   write(*,*) tpi,ens%insT,tstat%grpT(:)
  else
!$OMP SINGLE
#endif
  ttt = 0.0
  ttg(:) = 0.0
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
      do j=1,3
        temp = mass(i)*cart_v(i,j)*cart_v(i,j)
        ttt = ttt + temp
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + temp
      end do
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      temp = pstat%params(1)*bnd_v*bnd_v
      ttt = ttt + temp
      ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + temp
    else if (pstat%flag.eq.2) then
!     do nothing
    else if (pstat%flag.eq.3) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Unsupported manostat in randomize_cart&
 &_velocities(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
  ens%insT = ttt
  tstat%grpT(1:tstat%n_tgrps) = ttg(1:tstat%n_tgrps)
!
#ifdef ENABLE_THREADS
!$OMP END SINGLE
  end if
!
!$OMP SINGLE
#endif
  if ((3*n+ens%boxdof-n_constraints).lt.0) then
    write(ilog,*) 'Attempting a calculation with no remaining degree&
 &s of freedom. This is fatal.'
    call fexit()
  end if
  if (ens%insT.le.0.0) then
    write(ilog,*) 'System is initially void of kinetic energy even t&
 &hough there are unconstrained degrees of freedom. This is either a&
 & bug or a numerical precision error due to an extremely small temp&
 &erature request. Fatal exit.'
    call fexit()
  end if
  do i=1,tstat%n_tgrps
    if (tstat%grpT(i).le.0.0) then
      write(ilog,*) 'At least on of the T-coupling groups is initial&
 &ly void of kinetic energy even though there are unconstrained degr&
 &ees of freedom. This is either a bug or a numerical precision erro&
 &r due to an extremely small temperature request. Fatal exit.'
      call fexit()
    end if
  end do
  ens%insK = (1.0/u_dyn_fconv)*0.5*ens%insT
  ens%insT = ens%insT / (dble(3*n+ens%boxdof-n_constraints)*u_dyn_kb)
  tstat%grpT(:) = tstat%grpT(:)/(u_dyn_kb*tstat%grpdof(:))
  corrfac(:) = sqrt(kelvin/tstat%grpT(:))
#ifdef ENABLE_THREADS
!$OMP END SINGLE COPYPRIVATE(corrfac)
!
  if (sto2.ge.sta2) then
    mollo = molofrs(atmres(sta2))
    molhi = molofrs(atmres(sto2))
  else
    mollo = 1
    molhi = 0
  end if
  do imol=mollo,molhi
    do i=max(sta2,atmol(imol,1)),min(sto2,atmol(imol,2))
#else
  do imol=1,nmol
    do i=atmol(imol,1),atmol(imol,2)
#endif
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        cart_v(i,j) = corrfac(tstat%molgrp(imol))*cart_v(i,j)
      end do
    end do
  end do
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if (pstat%flag.eq.1) then
      bnd_v = corrfac(tstat%molgrp(nmol+1))*bnd_v
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
  end if
! 
end
!
!---------------------------------------------------------------------------
!
!  a subroutine to simply re-scale velocities by a T-increment
!
subroutine rescale_velocities(foreign_T,tpi)
!
  use system
  use forces
  use molecule
  use iounit
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
  RTYPE, INTENT(IN):: foreign_T
!
  integer imol,j,ttc
  RTYPE corrfac
#ifdef ENABLE_THREADS
  integer sta,sto
#endif
!
  corrfac = sqrt(kelvin/foreign_T)
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    sta = thr_limits(61,tpi)
    sto = thr_limits(62,tpi)
  else
    sta = 1
    sto = nmol
  end if
  do imol=sta,sto
#else
  do imol=1,nmol
#endif
!
!   first: center-of-mass translation
!
    do j=1,3
!     velocities in A/ps
      dc_di(imol)%v(j) = corrfac*dc_di(imol)%v(j)
    end do
!
    if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
!
!   second: rigid-body rotation

    if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
!
      ttc = 5
      write(ilog,*) 'Fatal. Diatomic molecules currently not support&
 &ed by dynamics implementation.'
      call fexit()
!     WARNING: support missing
!
    else
!
      ttc = 6
      do j=4,6
!       velocities in deg/ps
        dc_di(imol)%v(j) = corrfac*dc_di(imol)%v(j)
      end do
!
    end if
!
!   third: torsional degrees of freedom
!
    do j=1,ntormol(moltypid(imol))
      dc_di(imol)%v(ttc+j) = corrfac*dc_di(imol)%v(ttc+j) 
    end do
!
  end do
!
! fourth: box degrees of freedom (if applicable)
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      bnd_v = corrfac*bnd_v
    else if (pstat%flag.eq.3) then
      pstat%params(2) = corrfac*pstat%params(2)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in rescale_veloci&
 &ties(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!-----------------------------------------------------------------------
!
! the same for straight cartesian dynamics
!
subroutine rescale_cart_velocities(foreign_T,tpi)
!
  use system
  use forces
  use molecule
  use iounit
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
  RTYPE, INTENT(IN):: foreign_T

  integer imol
  RTYPE corrfac
#ifdef ENABLE_THREADS
  integer ixx
  logical moljflags(4)
#endif
!
  corrfac = sqrt(kelvin/foreign_T)
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    cart_v(thr_limits(1,tpi):thr_limits(2,tpi),1) = corrfac*cart_v(thr_limits(1,tpi):thr_limits(2,tpi),1)
    cart_v(thr_limits(1,tpi):thr_limits(2,tpi),2) = corrfac*cart_v(thr_limits(1,tpi):thr_limits(2,tpi),2)
    cart_v(thr_limits(1,tpi):thr_limits(2,tpi),3) = corrfac*cart_v(thr_limits(1,tpi):thr_limits(2,tpi),3)
  else
    cart_v(1:n,1:3) = corrfac*cart_v(1:n,1:3)
  end if
#else
  cart_v(1:n,1:3) = corrfac*cart_v(1:n,1:3)
#endif
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
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
      call update_comv(imol)
    end do
    moljflags(:) = .false.
    moljflags(2) = .true.
    do imol=1,nmlgs
      call molops_threads(imol,tpi,moljflags)
    end do
  else
#endif
  do imol=1,nmol
    call update_comv(imol)
  end do
#ifdef ENABLE_THREADS
  end if
#endif
!
! box degrees of freedom (if applicable)
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    if (pstat%flag.eq.1) then
      bnd_v = corrfac*bnd_v
    else if (pstat%flag.eq.2) then
!     do nothing for Markov chain
    else if (pstat%flag.eq.3) then
      pstat%params(2) = corrfac*pstat%params(2)
    else
      write(ilog,*) 'Fatal. Unsupported manostat in rescale_cart&
 &_velocities(...). This is most likely an omission bug.'
      call fexit()
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end
!
!-----------------------------------------------------------------------
!
! the subroutine to remove drift of the center of mass of the whole
! system (either as translation or rotation)
! not that drift is only computed over the intermolecular degrees of
! freedom
! in non-frictional dynamics, drift can easily absorb all of the kinetic
! energy of the system while effectively freezing the "real" degrees
! of freedom
!
subroutine drift_removal(fr1,fr2,tpi)
!
  use iounit
  use system
  use molecule
  use forces
  use mcsums
  use units
  use atoms
#ifdef ENABLE_THREADS
  use threads
  use cutoffs
  use sequen
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer j,imol,i,sta,sto
  RTYPE vavg(4),fr1,fr2,mte(3),ivmm(3)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer mollo,molhi,sta2,sto2
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, drift_removal(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(59,tpi)
    sto = thr_limits(60,tpi)
    if (thr_limits(2,tpi).ge.thr_limits(1,tpi)) then
      mollo = molofrs(atmres(thr_limits(1,tpi)))
      molhi = molofrs(atmres(thr_limits(2,tpi)))
    else
      mollo = 1
      molhi = 0
    end if
    sta2 = thr_limits(1,tpi)
    sto2 = thr_limits(2,tpi)
  else
    sta = 1
    sto = nmol
    mollo = 1
    molhi = nmol
    sta2 = 1
    sto2 = n
  end if
#else
!
  sta = 1
  sto = nmol
#endif
!
  if (n.eq.1) then
    fr1 = 0.0
    fr2 = 0.0
    return
  end if
  if ((nmol.eq.1).AND.(fycxyz.eq.1)) then
    fr1 = 0.0
    fr2 = 0.0
    return
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:4,tpi) = 0.0
    do imol=sta,sto
      thr_rutil(1:3,tpi) = thr_rutil(1:3,tpi) + movingmass(imol,1:3)*dc_di(imol)%v(1:3)
    end do
!$OMP BARRIER
    vavg(1:3) = sum(thr_rutil(1:3,1:thrdat%maxn),dim=2)
  else
#endif
  vavg(1:4) = 0.0
  do imol=sta,sto
    vavg(1:3) = vavg(1:3) + movingmass(imol,1:3)*dc_di(imol)%v(1:3)
  end do
#ifdef ENABLE_THREADS
  end if
#endif
!
  do j=1,3
    if (movingmass(nmol+1,j).gt.0.) then
      vavg(j) = vavg(j)/movingmass(nmol+1,j)
    else
      vavg(j) = 0.
    end if
  end do
!
  vavg(4) = movingmass(nmol+1,1)*vavg(1)*vavg(1) + movingmass(nmol+1,2)*vavg(2)*vavg(2) + movingmass(nmol+1,3)*vavg(3)*vavg(3)
  fr1 = 0.5*vavg(4)/(u_dyn_fconv*ens%insK)
  if (ens%sysfrz.ge.2) then
    if (fycxyz.eq.1) then
      do imol=sta,sto
        do j=1,3
          if (dc_di(imol)%frz(j).EQV..false.) dc_di(imol)%v(j) = dc_di(imol)%v(j) - vavg(j)
        end do
      end do
    else if (fycxyz.eq.2) then
#ifdef ENABLE_THREADS
      do imol=mollo,molhi
        ivmm(:) = 1.0/movingmass(imol,:)
        do i=max(sta2,atmol(imol,1)),min(sto2,atmol(imol,2))
          mte(:) = mass(i)*ivmm(:)
          do j=1,3
            if (cart_frz(i,j).EQV..false.) cart_v(i,j) = cart_v(i,j) - mte(j)*vavg(j)
          end do
        end do
      end do
#else
      do imol=1,nmol
        ivmm(:) = 1.0/movingmass(imol,:)
        do i=atmol(imol,1),atmol(imol,2)
          mte(:) = mass(i)*ivmm(:)
          do j=1,3
            if (cart_frz(i,j).EQV..false.) cart_v(i,j) = cart_v(i,j) - mte(j)*vavg(j)
          end do
        end do
      end do
#endif       
    end if
  end if
!
  fr2 = 0.
  if ((fycxyz.eq.1).AND.(nmol.gt.1)) then
    call rot_system_int(fr2,tpi)
  else if ((fycxyz.eq.2).AND.(n.gt.1)) then
    call rot_system_cart(fr2,tpi)
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice for degrees&
 & of freedom in drift_removal(...). This is a bug.'
    call fexit()
  end if
  fr2 = fr2/ens%insK
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!-----------------------------------------------------------------------
!
! note that this function cannot remove "internal" rotational drift
! potentially dangerous --> WARNING!!!!!!!!!!!
!
subroutine rot_system_int(erot,tpi)
!
  use iounit
  use forces
  use atoms
  use molecule
  use sequen
  use math
  use mcsums
  use system
  use units
#ifdef ENABLE_THREADS
  use threads
  use cutoffs
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer imol,i,j,athree,sta,sto
  RTYPE pos2(3),pos3(3),pos1(3),xr,yr,zr,tem
  RTYPE or_pl(3),nopl,erot
  RTYPE jte(3,3),jti(3,3),vel(3),vel2(3)
  logical afalse
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, rot_system_int(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(59,tpi)
    sto = thr_limits(60,tpi)
  else
    sta = 1
    sto = nmol
  end if
!$OMP BARRIER
#else
!
  sta = 1
  sto = nmol
#endif
!
  afalse = .false.
  athree = 3
!
  pos3(:) = 0.0
  vel(:) = 0.0
  jte(:,:) = 0.0
!
  do imol=sta,sto
    do j=1,3
      pos3(j) = pos3(j) + movingmass(imol,j)*comm(imol,j)
    end do
  end do
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:3,tpi) = pos3(:)
!$OMP BARRIER
    do j=1,3
      if (movingmass(nmol+1,j).gt.0.0) pos3(j) = sum(thr_rutil(j,1:thrdat%maxn))/movingmass(nmol+1,j)
    end do
  else
#endif
  do j=1,3
    if (movingmass(nmol+1,j).gt.0.0) pos3(j) = pos3(j)/movingmass(nmol+1,j)
  end do
#ifdef ENABLE_THREADS
  end if
#endif
!
! get net angular momentum in A^2 * g / (mol * ps)
  do imol=sta,sto
    do j=1,3
      pos2(j) = comm(imol,j)-pos3(j)
    end do
    do j=1,3
      pos1(j) = movingmass(imol,j)*dc_di(imol)%v(j)
    end do
    call crossprod(pos2,pos1,or_pl,nopl)
    do j=1,3
      vel(j) = vel(j) + or_pl(j)
    end do
  end do
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(4:6,tpi) = vel(1:3)
!$OMP BARRIER
    vel(:) = sum(thr_rutil(4:6,1:thrdat%maxn),dim=2)
  end if
#endif
! now get inertial tensor in A^2 * g / mol
  if (nmol.gt.2) then
    do imol=sta,sto
      xr = comm(imol,1)-pos3(1)
      yr = comm(imol,2)-pos3(2)
      zr = comm(imol,3)-pos3(3)
      jte(1,1) = jte(1,1) + molmass(moltypid(imol))*(yr*yr + zr*zr)
      jte(1,2) = jte(1,2) - molmass(moltypid(imol))*xr*yr
      jte(1,3) = jte(1,3) - molmass(moltypid(imol))*xr*zr
      jte(2,2) = jte(2,2) + molmass(moltypid(imol))*(xr*xr + zr*zr)
      jte(2,3) = jte(2,3) - molmass(moltypid(imol))*yr*zr
      jte(3,3) = jte(3,3) + molmass(moltypid(imol))*(xr*xr + yr*yr)
    end do
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      thr_rutil(7:9,tpi) = jte(:,1)
      thr_rutil(10:12,tpi) = jte(:,2)
      thr_rutil(13:15,tpi) = jte(:,3)
!$OMP BARRIER
      jte(:,1) = sum(thr_rutil(7:9,1:thrdat%maxn),dim=2)
      jte(:,2) = sum(thr_rutil(10:12,1:thrdat%maxn),dim=2)
      jte(:,3) = sum(thr_rutil(13:15,1:thrdat%maxn),dim=2)
    end if
#endif
    jte(3,2) = jte(2,3)
    jte(3,1) = jte(1,3)
    jte(2,1) = jte(1,2)
!   invert it
    call invert_inertial_tensor(jte,jti)
!   derive angular velocity in 1 / ps
    do i=1,3
      vel2(i) = 0.0
      do j=1,3
        vel2(i) = vel2(i) + jti(i,j)*vel(j)
      end do
    end do
  else if (nmol.eq.2) then
    tem = (comm(1,1) - comm(2,1))**2 +&
 &             (comm(1,2) - comm(2,2))**2 +&
 &             (comm(1,3) - comm(2,3))**2 
    tem = tem*molmass(moltypid(1))*molmass(moltypid(2))/(molmass(moltypid(1))+molmass(moltypid(2)))
    do j=1,3
      vel2(j) = vel(j)/tem
    end do
    if (movingmass(nmol+1,1).le.0.0) vel2(2:3) = 0.
    if (movingmass(nmol+1,3).le.0.0) vel2(1:2) = 0.
    if (movingmass(nmol+1,2).le.0.0) then
      vel2(1) = 0.
      vel2(3) = 0.
    end if
  end if
  erot = 0.0
  do i=1,3
     erot = erot + vel2(i)*vel(i)
  end do
  erot = 0.5*erot/u_dyn_fconv
! and finally remove by adjusting c.o.m.-velocties
  if (ens%sysfrz.eq.3) then
    do imol=sta,sto
      xr = comm(imol,1)-pos3(1)
      yr = comm(imol,2)-pos3(2)
      zr = comm(imol,3)-pos3(3)
      dc_di(imol)%v(1) = dc_di(imol)%v(1) - vel2(2)*zr + vel2(3)*yr
      dc_di(imol)%v(2) = dc_di(imol)%v(2) - vel2(3)*xr + vel2(1)*zr
      dc_di(imol)%v(3) = dc_di(imol)%v(3) - vel2(1)*yr + vel2(2)*xr
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same for full Cartesian description (simpler and slower)
!
subroutine rot_system_cart(erot,tpi)
!
  use iounit
  use forces
  use atoms
  use molecule
  use sequen
  use math
  use mcsums
  use system
  use units
#ifdef ENABLE_THREADS
  use threads
  use cutoffs
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer i,j,athree,sta,sto
  RTYPE pos2(3),pos3(3),pos1(3),xr,yr,zr,tem
  RTYPE or_pl(3),nopl,erot
  RTYPE jte(3,3),jti(3,3),vel(3),vel2(3)
  logical afalse
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, rot_system_cart(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
!$OMP BARRIER
#else
!
  sta = 1
  sto = n
#endif
!
  afalse = .false.
  athree = 3
!
  pos3(:) = 0.0
  vel(:) = 0.0
  jte(:,:) = 0.0
!
  do i=sta,sto
    pos3(1) = pos3(1) + mass(i)*x(i)
    pos3(2) = pos3(2) + mass(i)*y(i)
    pos3(3) = pos3(3) + mass(i)*z(i)
  end do
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    thr_rutil(1:3,tpi) = pos3(1:3)
!$OMP BARRIER
    do j=1,3
      if (movingmass(nmol+1,j).gt.0.0) then
        pos3(j) = sum(thr_rutil(j,1:thrdat%maxn))/movingmass(nmol+1,j)
      else
        pos3(j) = sum(thr_rutil(j,1:thrdat%maxn))
      end if
    end do
  else
#endif
  do j=1,3
    if (movingmass(nmol+1,j).gt.0.0) pos3(j) = pos3(j)/movingmass(nmol+1,j)
  end do
#ifdef ENABLE_THREADS
  end if
#endif
!
! now get inertial tensor in A^2 * g / mol and net angular momentum in A^2 * g / (mol * ps)
  if (n.gt.2) then
    do i=sta,sto
      xr = x(i)-pos3(1)
      yr = y(i)-pos3(2)
      zr = z(i)-pos3(3)
      do j=1,3
        if (cart_frz(i,j).EQV..false.) then
          pos1(j) = cart_v(i,j)*mass(i)
        else
          pos1(j) =  0.
        end if
      end do
      jte(1,1) = jte(1,1) + mass(i)*(yr*yr + zr*zr)
      jte(1,2) = jte(1,2) - mass(i)*xr*yr
      jte(1,3) = jte(1,3) - mass(i)*xr*zr
      jte(2,2) = jte(2,2) + mass(i)*(xr*xr + zr*zr)
      jte(2,3) = jte(2,3) - mass(i)*yr*zr
      jte(3,3) = jte(3,3) + mass(i)*(xr*xr + yr*yr)
      pos2(1) = xr
      pos2(2) = yr
      pos2(3) = zr
      call crossprod(pos2,pos1,or_pl,nopl)
      vel(:) = vel(:) + or_pl(:)
    end do
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      thr_rutil(4:6,tpi) = vel(:)
      thr_rutil(7:9,tpi) = jte(:,1)
      thr_rutil(10:12,tpi) = jte(:,2)
      thr_rutil(13:15,tpi) = jte(:,3)
!$OMP BARRIER
      vel(1:3) = sum(thr_rutil(4:6,1:thrdat%maxn),dim=2)
      jte(:,1) = sum(thr_rutil(7:9,1:thrdat%maxn),dim=2)
      jte(:,2) = sum(thr_rutil(10:12,1:thrdat%maxn),dim=2)
      jte(:,3) = sum(thr_rutil(13:15,1:thrdat%maxn),dim=2)
    end if
#endif
    jte(3,2) = jte(2,3)
    jte(3,1) = jte(1,3)
    jte(2,1) = jte(1,2)
!   invert it
    call invert_inertial_tensor(jte,jti)
!   derive angular velocity in 1 / ps
    do i=1,3
      vel2(i) = 0.0
      do j=1,3
        vel2(i) = vel2(i) + jti(i,j)*vel(j)
      end do
    end do
  else if (n.eq.2) then
    do i=1,2
      pos2(1) = mass(i)*(x(i)-pos3(1))
      pos2(2) = mass(i)*(y(i)-pos3(2))
      pos2(3) = mass(i)*(z(i)-pos3(3))
      pos1(:) = cart_v(i,:)
      call crossprod(pos2,pos1,or_pl,nopl)
      vel(:) = vel(:) + or_pl(:)
    end do
    tem = (x(1) - x(2))**2 +&
 &             (y(1) - y(2))**2 +&
 &             (z(1) - z(2))**2 
    tem = tem*mass(1)*mass(2)/(mass(1)+mass(2))
    do j=1,3
      vel2(j) = vel(j)/tem
    end do
    if (movingmass(nmol+1,1).le.0.0) vel2(2:3) = 0.
    if (movingmass(nmol+1,3).le.0.0) vel2(1:2) = 0.
    if (movingmass(nmol+1,2).le.0.0) then
      vel2(1) = 0.
      vel2(3) = 0.
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  end if
  erot = 0.0
  do i=1,3
     erot = erot + vel2(i)*vel(i)
  end do
  erot = 0.5*erot/u_dyn_fconv
! and finally remove by adjusting c.o.m.-velocties
  if (ens%sysfrz.eq.3) then
    do i=sta,sto
      xr = x(i)-pos3(1)
      yr = y(i)-pos3(2)
      zr = z(i)-pos3(3)
      cart_v(i,1) = cart_v(i,1) - vel2(2)*zr + vel2(3)*yr
      cart_v(i,2) = cart_v(i,2) - vel2(3)*xr + vel2(1)*zr
      cart_v(i,3) = cart_v(i,3) - vel2(1)*yr + vel2(2)*xr
    end do
!$OMP BARRIER
  end if
!
end
!
!-----------------------------------------------------------------------
!
! here we compute the (internal, potential) pressure tensor from the virial     
!
subroutine virial(virpot)
!
  use iounit
  use system
  use forces
  use atoms
  use molecule
!
  implicit none
!
  integer i,j,k
  RTYPE virpot(3,3)
!
! basically, every force which acts along a distance or directly on a positional
! coordinate will affect the virial
!
  if ((bnd_type.eq.3).OR.(bnd_type.eq.4)) then
    if (bnd_shape.ne.2) then
      write(ilog,*) 'Fatal. Encountered unsupported box shape in virial calculation. This is an omission bug.'
      call fexit()
    end if
!   this is really simple assuming we have all the relevant forces collected in cart_f
!   note that this includes boundary forces
    do j=1,3
      do k=1,3
        virpot(k,j) = 0.0
      end do
    end do
    do i=1,nmol
      do j=1,3
!       note that we compute the virial over the previous configuration, since all forces
!       are computed for that
        virpot(j,1) = virpot(j,1) + dc_di(i)%f(j)*&
 & (commref(i,j)-bnd_params(1))
        virpot(j,2) = virpot(j,2) + dc_di(i)%f(j)*&
 & (commref(i,j)-bnd_params(2))
        virpot(j,3) = virpot(j,3) + dc_di(i)%f(j)*&
 & (commref(i,j)-bnd_params(3))
      end do
    end do
    do i=1,3
      do j=1,3
!       subtract out the boundary term (which is what maintains mechanical equilibrium and
!       averages the pressure to be zero)
        virpot(i,j) = virpot(i,j) + bnd_f(i,j)
      end do
    end do
  else if (bnd_type.eq.1) then
!   with periodic images it is much more complicated
!   if we can adopt the single-sum implementation as shown in GROMACS it might be ok;
!   otherwise collecting on-the-fly becomes inevitable (position vectors need to be
!   minimum image-corrected ...)
    write(ilog,*) 'Fatal. Pressure calculation currently does not support periodic boundary conditions. Check back later.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! currently not in use
!
subroutine gen_cartv(mshfs)
!
  use iounit
  use system
  use atoms
  use molecule
  use mcsums
  use forces
!
  implicit none
!
  integer imol,k
  RTYPE mshfs(3,nmol),tst
!
! PBC
  if (bnd_type.eq.1) then
    if (nstep.eq.1) then
!     reconstructed velocities are techniqually for the half-step t + 1/2dt and not for time t
!     (as would be desired by most ensemble routines)
      cart_v(1:n,1) = (x(1:n) - xref(1:n))/dyn_dt
      cart_v(1:n,2) = (y(1:n) - yref(1:n))/dyn_dt
      cart_v(1:n,3) = (z(1:n) - zref(1:n))/dyn_dt
      do imol=1,nmol
        do k=1,3
          tst = mshfs(k,imol)/dyn_dt
          cart_v(atmol(imol,1):atmol(imol,2),k) = cart_v(atmol(imol,1):atmol(imol,2),k) - tst
        end do
      end do
    else
!     if we have an old cart_v, we can estimate the velocity at timestep t by averaging
!     over the velocities at times t - 1/2dt and t + 1/2dt
      cart_v(1:n,1) = 0.5*(cart_v(1:n,1)+(x(1:n)-xref(1:n))/dyn_dt)
      cart_v(1:n,2) = 0.5*(cart_v(1:n,2)+(y(1:n)-yref(1:n))/dyn_dt)
      cart_v(1:n,3) = 0.5*(cart_v(1:n,3)+(z(1:n)-zref(1:n))/dyn_dt)
      do imol=1,nmol
        do k=1,3
          tst = 0.5*mshfs(k,imol)/dyn_dt
          cart_v(atmol(imol,1):atmol(imol,2),k) = cart_v(atmol(imol,1):atmol(imol,2),k) - tst
        end do
      end do
    end if
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
    if (nstep.eq.1) then
!     reconstructed velocities are techniqually for the half-step t + 1/2dt and not for time t
!     (as would be desired by most ensemble routines)
      cart_v(1:n,1) = (x(1:n) - xref(1:n))/dyn_dt
      cart_v(1:n,2) = (y(1:n) - yref(1:n))/dyn_dt
      cart_v(1:n,3) = (z(1:n) - zref(1:n))/dyn_dt
    else
!     if we have an old cart_v, we can estimate the velocity at timestep t by averaging
!     over the velocities at times t - 1/2dt and t + 1/2dt
      cart_v(1:n,1) = 0.5*(cart_v(1:n,1)+(x(1:n)-xref(1:n))/dyn_dt)
      cart_v(1:n,2) = 0.5*(cart_v(1:n,2)+(y(1:n)-yref(1:n))/dyn_dt)
      cart_v(1:n,3) = 0.5*(cart_v(1:n,3)+(z(1:n)-zref(1:n))/dyn_dt)
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & and/or box shape in gen_cartv(...). Please check back later.'
    call fexit()
  end if 
!
end
!
!-----------------------------------------------------------------------
!
! here we compute the (ideal, kinetic) pressure tensor from the (Cartesian) velocities
!     
subroutine ekinten(virkin)
!
  use iounit
  use system
  use forces
  use atoms
  use mcsums
  use molecule
!
  implicit none
!
  integer j,imol
  RTYPE virkin(3,3)
!
  if (fycxyz.eq.1) then
    virkin(:,:) = 0.0
    do j=1,3
      do imol=1,nmol
        virkin(j,1)=virkin(j,1) + dc_di(imol)%v(j)*dc_di(imol)%v(1)*movingmass(imol,1)
        virkin(j,2)=virkin(j,2) + dc_di(imol)%v(j)*dc_di(imol)%v(2)*movingmass(imol,2)
        virkin(j,3)=virkin(j,3) + dc_di(imol)%v(j)*dc_di(imol)%v(3)*movingmass(imol,3)
      end do
    end do
  else if (fycxyz.eq.2) then
!   these are half-step velocities and hence slightly incorrect
    do j=1,3
      virkin(j,1) = sum(cart_v(1:n,j)*cart_v(1:n,1)*mass(1:n))
      virkin(j,2) = sum(cart_v(1:n,j)*cart_v(1:n,2)*mass(1:n))
      virkin(j,3) = sum(cart_v(1:n,j)*cart_v(1:n,3)*mass(1:n))
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice of degrees &
 &of freedom in ekinten(...). This is an omission bug.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine thermostat(tsc)
!
  use system
  use iounit
  use units
  use forces
  use system
  use atoms
  use molecule
  use math
  use mcsums
!
  implicit none
!
  RTYPE tsc(tstat%n_tgrps),random,normal,vind,decayrate,expt,exp2t,tgtt,rn1,gn1,rndgamma
  integer j,imol,i,ttc
  logical didup
!
  decayrate = dyn_dt/tstat%params(1)
!
  if (tstat%flag.eq.1) then
!   Berendsen is really simple
    tsc(:) = sqrt( 1.0 + decayrate*((kelvin/tstat%grpT(:)) - 1.0) )
  else if (tstat%flag.eq.2) then
!   Andersen is also simple: note that we do not use tsc() here, but modify velocities
!   directly, and that we ignore T-groups (since they are irrelevant, but can still be used
!   for analysis!)
    tsc(:) = 1.0
    if (fycxyz.eq.1) then
!     in non-Cartesian dynamics we have the problem of vastly inhomogeneous degrees of freedom,
!     and we'll therefore couple each d.o.f. independently 
!     we also assume the inertial masses are xyz up-to-date (T-stat called soon after cart2int()!)
      do imol=1,nmol
!       first: center-of-mass translation
        do j=1,3
!         velocities in A/ps
          if (dc_di(imol)%frz(j).EQV..true.) then
            cycle
          end if
          if (random().lt.decayrate) then
            vind = normal()
            vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
            dc_di(imol)%v(j) = vind
          end if
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
!       second: rigid-body rotation
        if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          call fexit()
!         WARNING: support missing
        else
          ttc = 6
          do j=4,6
!           velocities in deg/ps
            if (dc_di(imol)%frz(j).EQV..true.) then
              cycle
            end if
            if (random().lt.decayrate) then
              vind = normal()
              vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(j))
              dc_di(imol)%v(j) = RADIAN*vind
            end if
          end do
        end if
!       third: torsional degrees of freedom
        do j=1,ntormol(moltypid(imol))
          if (dc_di(imol)%frz(ttc+j).EQV..true.) then
            cycle
          end if
          if (random().lt.decayrate) then
            vind = normal()
            vind = vind*sqrt(u_dyn_kb*kelvin/dc_di(imol)%im(ttc+j))
            dc_di(imol)%v(ttc+j) = RADIAN*vind
          end if
        end do
!       fourth: box degrees of freedom (if applicable)
        if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
          if (random().lt.decayrate) then
            vind = normal()
            vind = vind*sqrt(u_dyn_kb*kelvin/pstat%params(1))
            bnd_v = vind
          end if
        end if
      end do
    else if (fycxyz.eq.2) then
!     in Cartesian dynamics we'll simply couple each atom individually
      do imol=1,nmol
        didup = .false.
        do i=atmol(imol,1),atmol(imol,2)
          if (mass(i).le.0.0) cycle
          if (random().lt.decayrate) then
            didup = .true.
            do j=1,3
!             velocities in A/ps
              vind = normal()
              vind = vind*sqrt(u_dyn_kb*kelvin/mass(i))
              cart_v(i,j) = vind
            end do
          end if
        end do
        if (didup.EQV..true.) then
          call update_comv(imol)
        end if
      end do
    end if
! for the Nose-Hoover analog method (by Stern), we'll simply set this to what is needed for the integrator
  else if (tstat%flag.eq.3) then
    tsc(:) = exp(-tstat%params(2)*0.5*dyn_dt)
! Bussi is really simple as well - note they strongly recommend using the exact (analytical) form derived in the appendix
! if all Gaussian random numbers are zero, and the gamma random # is ~ grpdof, analysis of leading terms recovers Berendsen
! conversely, if decayrate -> 0.0, NVE is recovered
  else if (tstat%flag.eq.4) then
    expt = exp(-decayrate)
    exp2t = 2.0*exp(-0.5*decayrate)
    do i=1,tstat%n_tgrps
      tgtt = kelvin/(tstat%grpT(i)*tstat%grpdof(i))
!     note that the factor of 2.0 in front of rndgamma is not immediately apparent from the reference
!     but is crucial for relating sums of squared Gaussian RNs to values pulled from the gamma-dist.
      if (tstat%grpdof(i).eq.1) then
        gn1 = 0.0
      else if (tstat%grpdof(i).eq.2) then
        rn1 = normal()
        gn1 = rn1*rn1
      else if (mod(floor(tstat%grpdof(i)+1.0e-7)-1,2).eq.0) then
        gn1 = 2.0*rndgamma((1.0*tstat%grpdof(i)-1.0)/2.0)
!     this split may not be necessary
      else
        rn1 = normal()
        gn1 = 2.0*rndgamma((1.0*tstat%grpdof(i)-2.0)/2.0) + rn1*rn1
      end if
      rn1 = normal()
      tsc(i) = sqrt(expt + tgtt*(1.0-expt)*(rn1*rn1 + gn1) + exp2t*rn1*sqrt(tgtt*(1.0-expt)))
    end do
  else
    write(ilog,*) 'Encountered unsupported thermostat flag. Offendin&
 &g code is ',tstat%flag,'. Please report this problem!'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine gen_rd_for_tstat()
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
  RTYPE random,normal,decayrate
  integer j,imol,i,ttc
!
  decayrate = dyn_dt/tstat%params(1)
!
  if (tstat%flag.eq.2) then
    if (fycxyz.eq.1) then
      i = 0
      do imol=1,nmol
        do j=1,3
          if (dc_di(imol)%frz(j).EQV..true.) then
            cycle
          end if
          if (random().lt.decayrate) then
            i = i + 1 
            rndstat%rbuf(i) = normal()
          end if
        end do
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) cycle
        if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          call fexit()
        else
          ttc = 6
          do j=4,6
            if (dc_di(imol)%frz(j).EQV..true.) then
              cycle
            end if
            if (random().lt.decayrate) then
              i = i + 1 
              rndstat%rbuf(i) = normal()
            end if
          end do
        end if
        do j=1,ntormol(moltypid(imol))
          if (dc_di(imol)%frz(ttc+j).EQV..true.) then
            cycle
          end if
          if (random().lt.decayrate) then
            i = i + 1 
            rndstat%rbuf(i) = normal()
          end if
        end do
        if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
          if (random().lt.decayrate) then
            i = i + 1 
            rndstat%rbuf(i) = normal()
          end if
        end if
      end do
    else if (fycxyz.eq.2) then
      j = 0 
      do i=1,n
        if (mass(i).le.0.0) cycle
        if (random().lt.decayrate) then
          j = j + 3
          rndstat%rbuf(j-2) = normal()
          rndstat%rbuf(j-1) = normal()
          rndstat%rbuf(j) = normal()
        end if
      end do
    end if
  else
    write(ilog,*) 'Encountered unsupported thermostat flag in gen_rd_for_tstat(...). Offending code is ',tstat%flag,'. &
 &Please report this problem!'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
subroutine init_thermostat()
!
  use system
  use iounit
  use molecule
  use forces
  use torsn
  use atoms
  use shakeetal
  use sequen
!
  implicit none
!
  integer ttg(tstat%n_tgrps),i,ttf(tstat%n_tgrps),ttc,imol,j
  integer drift_c,alldof
  logical samxyz(3)
!
  ttg(:) = 0
  ttf(:) = 0
  drift_c = 0
! note that we have to spread out the drift-correction constraints onto
! the T-coupling groups
  if ((nmol.gt.1).AND.(fycxyz.eq.1)) then
    samxyz(:) = .false.
    do i=1,nmol
      if (dc_di(i)%frz(1).EQV..false.) samxyz(1) = .true.
      if (dc_di(i)%frz(2).EQV..false.) samxyz(2) = .true.
      if (dc_di(i)%frz(3).EQV..false.) samxyz(3) = .true.
    end do
    j = 0
    do i=1,3
      if (samxyz(i).EQV..false.) j = j + 1 ! j can never be 3
    end do
    if (ens%sysfrz.eq.3) then
!     if there's exactly two molecules one rotational dof is missing
      if (nmol.eq.2) then
        drift_c = 5 - j - min(j,1)*j
      else
        drift_c = 6 - j - min(j,1)*(j+1)
      end if
    else if (ens%sysfrz.eq.2) then
      drift_c = 3 - j
    else
      drift_c = 0
    end if
  else if (fycxyz.eq.1) then
    drift_c = 0 ! drift_c are lumped in dc_di%frz
  else
    samxyz(:) = .false.
    do i=1,n
      if (cart_frz(i,1).EQV..false.) samxyz(1) = .true.
      if (cart_frz(i,2).EQV..false.) samxyz(2) = .true.
      if (cart_frz(i,3).EQV..false.) samxyz(3) = .true.
    end do
    j = 0
    do i=1,3
      if (samxyz(i).EQV..false.) j = j + 1 ! j can never be 3
    end do
    if (ens%sysfrz.eq.3) then
!     if there's exactly two molecules one rotational dof is missing
      if (n.eq.2) then
        drift_c = 5 - j - min(j,1)*j
      else if (n.eq.1) then
        drift_c = 3 - j
      else
        drift_c = 6 - j - min(j,1)*(j+1)
      end if
    else if (ens%sysfrz.eq.2) then
      drift_c = 3 - j
    end if
  end if
  if (fycxyz.eq.1) then
    do imol=1,nmol
      do j=1,3
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
        if (dc_di(imol)%frz(j).EQV..true.) then
          ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + 1
        end if
      end do
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
        write(ilog,*) 'Fatal. Diatomic molecules are not yet support&
 &ed in dynamics. Check back later.'
        call fexit()
!       WARNING: support missing
      else
        ttc = 6
        do j=4,6
          ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
          if (dc_di(imol)%frz(j).EQV..true.) then
            ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + 1
          end if
        end do    
      end if
      do j=ttc+1,ttc+ntormol(moltypid(imol))
        ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
        if (dc_di(imol)%frz(j).EQV..true.) then
          ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + 1
        end if
      end do
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + 1
      end if
    end if
    alldof = totrbd+ndyntorsn-n_constraints+ens%boxdof
    if ((sum(ttg)-sum(ttf)-drift_c).ne.alldof) then
      write(ilog,*) 'Fatal. Setup of degrees of freedom is inconsist&
 &ent. This is most certainly a bug.'
       call fexit()
    end if
  else if (fycxyz.eq.2) then
    do imol=1,nmol
      do i=atmol(imol,1),atmol(imol,2)
        do j=1,3
          if (cart_frz(i,j).EQV..false.) ttg(tstat%molgrp(imol)) = ttg(tstat%molgrp(imol)) + 1
        end do
      end do
    end do
    do i=1,cart_cons_grps
      if (constraints(i)%nr.gt.0) then
        imol = molofrs(atmres(constraints(i)%idx(1,1)))
        ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + constraints(i)%nr
      end if
      if (constraints(i)%nr3.gt.0) then
        imol = molofrs(atmres(constraints(i)%idx3(1,1))) ! the same
        ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + constraints(i)%nr3
      end if
      if (constraints(i)%nr4.gt.0) then
        imol = molofrs(atmres(constraints(i)%idx4(1,1))) ! the same
        ttf(tstat%molgrp(imol)) = ttf(tstat%molgrp(imol)) + constraints(i)%nr4
      end if
    end do
    if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (pstat%flag.eq.1) then
        ttg(tstat%molgrp(nmol+1)) = ttg(tstat%molgrp(nmol+1)) + 1
      end if
    end if
    alldof = 3*n-n_constraints+ens%boxdof
    if ((sum(ttg)-sum(ttf)-drift_c).ne.alldof) then
      write(ilog,*) 'Fatal. Setup of degrees of freedom is inconsist&
 &ent. This is most certainly a bug.'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice of degrees &
 &of freedom in init_thermostat(). This is an omission bug.'
    call fexit()
  end if
!
  do i=1,tstat%n_tgrps
    tstat%grpT(i) = 0.0
    tstat%grpdof(i) = ttg(i) - ttf(i) - &
 &    dble(drift_c*(ttg(i) - ttf(i)))/dble(alldof+drift_c)
    if (tstat%grpdof(i).le.0.0) then
      write(ilog,*) 'Fatal. T-coupling group #',i,' has no active de&
 &grees of freedom. Check input files for unwanted constraints.'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine manostat(boxovol)
!
  use system
  use iounit
  use units
  use forces
  use system
  use energies
  use molecule
  use math
  use atoms
!
  implicit none
!
  RTYPE netf,deltaV,boxovol,oldbnd(3),evn(MAXENERGYTERMS)
  RTYPE random,scf,diff
  logical atrue,afalse
!
  atrue = .true.
  afalse = .false.
  bnd_params(8) = 1000.0
!
! in constant pressure simulations update box
  if (bnd_shape.eq.2) then
    if (pstat%flag.eq.1) then
      netf= bnd_fr - 4.0*PI*bnd_params(5)*extpress/u_dyn_virconv
      bnd_v= bnd_v + 0.5*u_dyn_fconv*dyn_dt*netf/pstat%params(1)
      bnd_params(4) = bnd_params(4) + dyn_dt*bnd_v
      call update_bound(afalse)
    else if (pstat%flag.eq.2) then
!     backup old radius
      oldbnd(1) = bnd_params(4)
      deltaV = (random()-0.5)*bnd_params(8)
      bnd_params(4) = ((3./(4.0*PI))*(boxovol+deltaV))**(1./3.)
      call update_bound(afalse)
      call e_boundary(evn)
      diff = 3.0*(evn(12)-esterms(12)) + &
 &                       3.0*deltaV*extpress/u_dyn_virconv
      invtemp = 1./(gasconst*ens%insT)
      if ((diff.le.0.0).OR.(random().lt.exp(-diff*invtemp))) then
        esterms(12) = evn(12)
      else
        bnd_params(4) = oldbnd(1)
        call update_bound(afalse)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported manostat for cho&
 &sen box and ensemble (manostat-flag: ',pstat%flag,').'
      call fexit()
    end if
  else if (bnd_shape.eq.1) then
    if (pstat%flag.eq.2) then
!     backup old radius
      oldbnd(1:3) = bnd_params(1:3)
      deltaV = (random()-0.5)*bnd_params(8)
      if ((boxovol+deltaV).le.0.0) then
        write(ilog,*) 'Fatal. Volume fluctuations yielded zero or ne&
 &gative volume in Markov chain. This indicates a poorly set up simu&
 &lation. Please check settings!'
        call fexit()
      end if
!     isotropic
      scf = ((boxovol+deltaV)**(1./3.))/bnd_params(1)
      bnd_params(1:3) = scf*bnd_params(1:3)
      call update_bound(afalse)
      diff = -3.0*deltaV*extpress/u_dyn_virconv
      invtemp = 1./(gasconst*ens%insT)
      if ((diff.le.0.0).OR.(random().lt.exp(-diff*invtemp))) then
        call rescale_xyz(scf,atrue)
      else
        bnd_params(1:3) = oldbnd(1:3)
        call update_bound(afalse)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported manostat for cho&
 &sen box and ensemble (manostat-flag: ',pstat%flag,').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary shape for&
 & chosen ensemble (Flag: ',ens%flag,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine set_movingmass()
!
  use molecule
  use system
  use forces
  use atoms
!
  implicit none
!
  integer imol,k
!
  movingmass(:,:) = 0.0
  do imol=1,nmol
    if (fycxyz.eq.2) then
      do k=atmol(imol,1),atmol(imol,2)
        if (cart_frz(k,1).EQV..false.) movingmass(imol,1) = movingmass(imol,1) + mass(k)
        if (cart_frz(k,2).EQV..false.) movingmass(imol,2) = movingmass(imol,2) + mass(k)
        if (cart_frz(k,3).EQV..false.) movingmass(imol,3) = movingmass(imol,3) + mass(k)
      end do
    else
      if (dc_di(imol)%frz(1).EQV..false.) movingmass(imol,1) = molmass(moltypid(imol))
      if (dc_di(imol)%frz(2).EQV..false.) movingmass(imol,2) = molmass(moltypid(imol))
      if (dc_di(imol)%frz(3).EQV..false.) movingmass(imol,3) = molmass(moltypid(imol))
    end if
  end do
  do k=1,3
    movingmass(nmol+1,k) = sum(movingmass(1:nmol,k))
  end do
!
end
!
!-----------------------------------------------------------------------
!
