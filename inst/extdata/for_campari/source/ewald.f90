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
!
subroutine Straight_PBCSum(evec)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use system
  use energies
  use math
  use ewalds
  use iounit
  use sequen
#ifdef PGI_FORTRAN
  use libm
#endif
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
  integer i,k,ii,kk,xlo,ylo,zlo,xhi,yhi,zhi,xi,yi,zi,j,fp,azero
  integer rs,rs2
  integer(KIND=8) t1,t2
  RTYPE dvec(3),svec(3),evec(MAXENERGYTERMS),term0,dis
  RTYPE cut2,voli,dis2
  RTYPE tem0
  RTYPE ca1(n,3),ca2(3,n),ca3(3,n)
!
  azero = 0
  fp = 46
!
  if ((bnd_type.ne.1).OR.(bnd_shape.ne.1)) call fexit()
!
  xlo = -1
  ylo = -1
  zlo = -1
  xhi = 1
  yhi = 1
  zhi = 1
  cut2 = mcnb_cutoff2
  voli = bnd_params(1)*bnd_params(2)*bnd_params(3)
!
  call System_Clock(t1)
! central unit cell
  evec(6) = 0.0
  cart_f(:,:) = 0.0
  do rs=1,nseq
   do i=1,at(rs)%npol
    ii = at(rs)%pol(i)
    do rs2=rs+1,nseq
     do k=1,at(rs2)%npol
      kk = at(rs2)%pol(k)
!      call dis_bound3(ii,kk,dvec,dis2)
      dvec(1) = x(kk) - x(ii)! + svec(1)
      dvec(2) = y(kk) - y(ii)! + svec(2)
      dvec(3) = z(kk) - z(ii)! + svec(3)
      dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2)&
 &                                     + dvec(3)*dvec(3)
      dis = sqrt(dis2)
      term0 = electric*atq(kk)*atq(ii)/dis
      evec(6) = evec(6) + scale_POLAR*term0
!      do j=1,3
!        cart_f(ii,j) = cart_f(ii,j) -
! &            scale_POLAR*(term0/dis2)*dvec(j)
!        cart_f(kk,j) = cart_f(kk,j) +
! &            scale_POLAR*(term0/dis2)*dvec(j)
!      end do
     end do
    end do
   end do
  end do
  tem0 = 0.0
  do xi=xlo,xhi
    do yi=ylo,yhi
      do zi=zlo,zhi
        if ((xi.eq.0).AND.(yi.eq.0).AND.(zi.eq.0)) cycle
        svec(1) = dble(xi)*bnd_params(1)
        svec(2) = dble(yi)*bnd_params(2)
        svec(3) = dble(zi)*bnd_params(3)
        do rs=1,nseq
         do i=1,at(rs)%npol
          ii = at(rs)%pol(i)
          do rs2=1,nseq
           do k=1,at(rs2)%npol
            kk = at(rs2)%pol(k)
            dvec(1) = x(kk) - x(ii) + svec(1)
            dvec(2) = y(kk) - y(ii) + svec(2)
            dvec(3) = z(kk) - z(ii) + svec(3)
            dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2)&
 &                                     + dvec(3)*dvec(3)
            dis = sqrt(dis2)
            tem0 = tem0 + electric*atq(ii)*atq(kk)
            term0 = electric*atq(ii)*atq(kk)/dis
            evec(6) = evec(6) + scale_POLAR*term0
!            do j=1,3
!              cart_f(ii,j) = cart_f(ii,j) -
! &               scale_POLAR*(term0/dis2)*dvec(j)
!              cart_f(kk,j) = cart_f(kk,j) +
! &               scale_POLAR*(term0/dis2)*dvec(j)
!            end do
           end do
          end do
         end do
        end do
      end do
    end do
  end do
  call System_Clock(t2)
!
  write(*,*) 'DUMB: ',t2-t1,evec(6)
  evec(:) = 0.0
  return
!
  write(*,*) cart_f(cglst%it(fp),:)
!
  call System_Clock(t1)
! central unit cell with minimum image and no cutoffs
  evec(6) = 0.0
  ca1(:,:) = 0.0
  do i=1,cglst%ncs
    ii = cglst%it(i)
    do k=i+1,cglst%ncs
      kk = cglst%it(k)
      call dis_bound3(ii,kk,dvec,dis2)
      dis = sqrt(dis2)
      term0 = electric*cglst%tc(k)*cglst%tc(i)/dis
      evec(6) = evec(6) + scale_POLAR*term0
      do j=1,3
        ca1(ii,j) = ca1(ii,j) -&
 &            scale_POLAR*(term0/dis2)*dvec(j)
        ca1(kk,j) = ca1(kk,j) +&
 &            scale_POLAR*(term0/dis2)*dvec(j)
      end do
    end do
  end do
  call System_Clock(t2)
  write(*,*) 'MINI: ',t2-t1,evec(6)
  write(*,*) ca1(cglst%it(fp),:)
!
! standard Ewald
  ca2(:,:) = 0.0
  call setup_ewald_constV()
  call System_Clock(t1)
  call force_ewald(evec,ca2,azero)
  call System_Clock(t2)
!
  write(*,*) 'EWAL: ',t2-t1,evec(6)
  write(*,*) ca2(:,cglst%it(fp))
!
!
#ifdef LINK_FFTW
  ca3(:,:) = 0.0
  call setup_pme_constV()
  call System_Clock(t1)
  call force_pme(evec,ca3,azero)
  call System_Clock(t2)
  write(*,*) 'PME : ',t2-t1,evec(6)
  write(*,*) ca3(:,cglst%it(fp))
!
#endif
!
end
!
!-----------------------------------------------------------------------
!
subroutine force_ewald(evec,ca_f,tpi2)
!
  use cutoffs
  use ewalds
  use math
  use units
  use system
  use atoms
  use forces
  use sequen
  use polypep
  use energies
#ifdef ENABLE_THREADS
  use threads
  use iounit
#endif
  use mcsums
!
  implicit none
!
  integer,INTENT(IN):: tpi2
!
  integer xi,yi,zi,i,ii,rs,sta,sto,shfv(3),tpi
  RTYPE voli,ivec(3),mvec(3),bterm,term0p,ewrat,terms
  RTYPE evec(MAXENERGYTERMS),recsum,m2
  RTYPE ca_f(3,n)
  complex(KIND=8) exptm,exptms,sm1,twopii,term0(3),SFac,SFacs
  PARAMETER (sm1=(0,1))
  integer(KIND=8) t1,t2
#ifdef ENABLE_THREADS
  integer tpn
  integer(KIND=8) ttimer
  logical OMP_IN_PARALLEL
!
  if ((tpi2.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, force_ewald(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi2.gt.0) then
    if (thrdat%subnrs(1).gt.0) then
      tpi = tpi2 + thrdat%maxn - thrdat%subnrs(1)
      tpn = thrdat%subnrs(1)
    else
      tpi = tpi2
      tpn = thrdat%maxn
    end if
    sta = thr_limits(27,tpi)
    sto = thr_limits(28,tpi)
    if (tpi2.eq.1) call System_Clock(t1)
  else
    sta = 1
    sto = nseq
    tpn = 1
    tpi = 1
    call System_Clock(t1)
  end if
#else
  sta = 1
  sto = nseq
  tpi = 1
  call System_Clock(t1)
#endif

!
! some initialization
  voli = bnd_params(1)*bnd_params(2)*bnd_params(3)
  bterm = scale_POLAR*electric/(2.0*PI*voli)
  ewrat = PI*PI/(ewpm*ewpm)
  twopii = 2.0*PI*sm1
  shfv(:) = -kdims(1:3,1) + 1
!
! the constant term is ... well ... constant
  if (tpi2.le.1) evec(6) = evec(6) + ewcnst
!
! reciprocal space setup and sum
  recsum = 0.0
#ifdef ENABLE_THREADS
!  if (sta.le.sto) then
!    ewSvd(:,at(sta)%bb(1):(at(sto)%bb(1)+at(sto)%na-1)) = 0.0
!    ewSvds(:,at(sta)%bb(1):(at(sto)%bb(1)+at(sto)%na-1)) = 0.0
!  end if
  ewSfac(:,:,:,tpi) = 0.0
  ewSfacs(:,:,:,tpi) = 0.0
  if (thr_dlb(5,1).gt.0) then
    if (tpi2.eq.1) thr_dlb(5,2) = thr_dlb(5,2) + 1
    call System_Clock(count=ttimer)
    thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
  end if
#else
!  ewSvd(:,:) = 0.0
!  ewSvds(:,:) = 0.0
  ewSfac(:,:,:,1) = 0.0
  ewSfacs(:,:,:,1) = 0.0
#endif
  do xi=kdims(1,1),kdims(1,2)
    do yi=kdims(2,1),kdims(2,2)
      do zi=kdims(3,1),kdims(3,2)
        if ((xi.eq.0).AND.(yi.eq.0).AND.(zi.eq.0)) cycle
        ivec(1) = dble(xi)/bnd_params(1)
        ivec(2) = dble(yi)/bnd_params(2)
        ivec(3) = dble(zi)/bnd_params(3)
!       get the structure factors
        do rs=sta,sto
          do i=1,at(rs)%npol
            ii = at(rs)%pol(i)
            mvec(1) = (x(ii)-bnd_params(4))*ivec(1)
            mvec(2) = (y(ii)-bnd_params(5))*ivec(2)
            mvec(3) = (z(ii)-bnd_params(6))*ivec(3)
            exptm = atq(ii)*exp(twopii*sum(mvec))
            exptms = atq(ii)*exp(-twopii*sum(mvec))
            ewSfac(xi+shfv(1),yi+shfv(2),zi+shfv(3),tpi) = ewSfac(xi+shfv(1),yi+shfv(2),zi+shfv(3),tpi) + exptm
            ewSfacs(xi+shfv(1),yi+shfv(2),zi+shfv(3),tpi) = ewSfacs(xi+shfv(1),yi+shfv(2),zi+shfv(3),tpi) + exptms
          end do
        end do
      end do
    end do
  end do
#ifdef ENABLE_THREADS
  if (thr_dlb(5,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
  end if
!$OMP BARRIER
  if (thr_dlb(5,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
  end if
#endif
  do xi=kdims(1,1),kdims(1,2)
    do yi=kdims(2,1),kdims(2,2)
      do zi=kdims(3,1),kdims(3,2)
        if ((xi.eq.0).AND.(yi.eq.0).AND.(zi.eq.0)) cycle
        ivec(1) = dble(xi)/bnd_params(1)
        ivec(2) = dble(yi)/bnd_params(2)
        ivec(3) = dble(zi)/bnd_params(3)
        m2 = sum(ivec(:)*ivec(:))
        term0p = bterm*(1.0/m2)*exp(-ewrat*m2)
        term0(:) = twopii*ivec(:)*term0p
#ifdef ENABLE_THREADS
        Sfac = sum(ewSfac(xi+shfv(1),yi+shfv(2),zi+shfv(3),1:tpn))
        Sfacs = sum(ewSfacs(xi+shfv(1),yi+shfv(2),zi+shfv(3),1:tpn))
        terms = REAL(Sfac*Sfacs)
!$OMP MASTER
        recsum = recsum + terms*term0p
!$OMP END MASTER
#else
        Sfac = ewSfac(xi+shfv(1),yi+shfv(2),zi+shfv(3),1)
        Sfacs = ewSfacs(xi+shfv(1),yi+shfv(2),zi+shfv(3),1)
        terms = REAL(Sfac*Sfacs)
        recsum = recsum + terms*term0p
#endif
        do rs=sta,sto
          do i=1,at(rs)%npol
            ii = at(rs)%pol(i)
!           use partials of the structure factor and structure factor itself to get
!           Cartesian force
            mvec(1) = (x(ii)-bnd_params(4))*ivec(1)
            mvec(2) = (y(ii)-bnd_params(5))*ivec(2)
            mvec(3) = (z(ii)-bnd_params(6))*ivec(3)
            exptm = atq(ii)*exp(twopii*sum(mvec))
            exptms = atq(ii)*exp(-twopii*sum(mvec))
            ca_f(:,ii) = ca_f(:,ii) - term0*(Sfacs*exptm - Sfac*exptms)
          end do
        end do
      end do
    end do
  end do
#ifdef ENABLE_THREADS
  if (thr_dlb(5,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
  end if
!$OMP MASTER
  evec(6) = evec(6) + recsum
!$OMP END MASTER
#else
  evec(6) = evec(6) + recsum
#endif
!  write(*,*) evec(6)
!
#ifdef ENABLE_THREADS
  if (tpi2.le.1) then
    call System_Clock(t2)
    time_pme = time_pme + 1.0*(t2-t1)
    time_energy = time_energy + 1.0*(t1-t2)
  end if
#else
  call System_Clock(t2)
  time_pme = time_pme + 1.0*(t2-t1)
  time_energy = time_energy + 1.0*(t1-t2)
#endif
!
end
!
!-----------------------------------------------------------------------
!
! WARNING!!!! this routine does not work for asymmetric grids and odd spline
! order. clarify why (real-to-complex FFT???)
!
! when using threads, we have to face the inconvenience that FFTW will spawn
! its own team of threads of unknown size (scaling is generally not good and it may often use fewer threads than possible)
! in order to keep this procedure reasonably efficient, OMP_PROC_BIND must be false OR a MASTER construct
! must be used instead of the more appopriate SINGLE construct around the fftw_execute calls
! for small and medium systems, the parallel performance is likely to be limited by the needless overhead
!
subroutine force_pme(evec,ca_f,tpi2)
!
#ifdef LINK_FFTW
!
  use cutoffs
  use ewalds
  use sequen
  use polypep
  use math
  use units
  use system
  use atoms
  use forces
  use energies
  use mcsums
#ifdef ENABLE_THREADS
  use threads
  use iounit
#endif
  use, intrinsic :: ISO_C_BINDING
!
!
  implicit none
#ifdef ORACLE_FORTRAN
#include "fftw3_SUN.f03"
#else
#include "fftw3.f03"
#endif
!
  integer, INTENT(IN):: tpi2
!
  integer xi,yi,zi,i,j,k,l,ii,cc,gdz(3),sh,px,py,pz
  integer rs,sta,sto,sta2,sto2
  RTYPE mvec(3),ewrat,term0,sqrtel,ts(3)
  RTYPE evec(MAXENERGYTERMS),recsum,gdzv(3),ewtmpf(3),catmpf(3)
  RTYPE ca_f(3,n),temj,temk,teml,tem,bsplbu(splor)
  integer xmap(3*kdims(1,2)),ymap(3*kdims(2,2)),zmap(3*kdims(3,2))
  logical zyzero
  integer(KIND=8):: t1,t2
#ifdef ENABLE_THREADS
  integer tpi
  integer(KIND=8) ttimer
  logical OMP_IN_PARALLEL
!
  if ((tpi2.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, force_pme(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi2.gt.0) then
    if (thrdat%subnrs(1).gt.0) then
      tpi = tpi2 + thrdat%maxn - thrdat%subnrs(1)
    else
      tpi = tpi2
    end if
    sta = thr_limits(27,tpi)
    sto = thr_limits(28,tpi)
    sta2 = thr_limits(29,tpi)
    sto2 = thr_limits(30,tpi)
    if (tpi2.eq.1) then
      call System_Clock(t1)
      ewnetf(:) = 0.0
    end if
  else
    sta = 1
    sto = nseq
    sta2 = 1
    sto2 = kdims(3,2)
    call System_Clock(t1)
    tpi = 0
  end if
#else
  sta = 1
  sto = nseq
  sta2 = 1
  sto2 = kdims(3,2)
  call System_Clock(t1)
#endif
!
  sqrtel = sqrt(electric)
  sh = splor - 1! + mod(splor/2)
  ewrat = PI*PI/(ewpm*ewpm)
  do j=1,3
    gdz(j) = kdims(j,2)
    gdzv(j) = dble(gdz(j))/bnd_params(j)
  end do
! move to global / setup
  do i=1,3*gdz(1)
    xi = i - sh
    if (xi.gt.gdz(1)) xi = xi - gdz(1)
    if (xi.gt.gdz(1)) xi = xi - gdz(1)
    if (xi.gt.gdz(1)) xi = xi - gdz(1)
    if (xi.lt.1) xi = xi + gdz(1)
    if (xi.lt.1) xi = xi + gdz(1)
    if (xi.lt.1) xi = xi + gdz(1)
    xmap(i) = xi
  end do
  do i=1,3*gdz(2)
    yi = i  - sh
    if (yi.gt.gdz(2)) yi = yi - gdz(2)
    if (yi.gt.gdz(2)) yi = yi - gdz(2)
    if (yi.gt.gdz(2)) yi = yi - gdz(2)
    if (yi.lt.1) yi = yi + gdz(2)
    if (yi.lt.1) yi = yi + gdz(2)
    if (yi.lt.1) yi = yi + gdz(2)
    ymap(i) = yi
  end do
  do i=1,3*gdz(3)
    zi = i  - sh
    if (zi.gt.gdz(3)) zi = zi - gdz(3)
    if (zi.gt.gdz(3)) zi = zi - gdz(3)
    if (zi.gt.gdz(3)) zi = zi - gdz(3)
    if (zi.lt.1) zi = zi + gdz(3)
    if (zi.lt.1) zi = zi + gdz(3)
    if (zi.lt.1) zi = zi + gdz(3)
    zmap(i) = zi
  end do

! initialize increment variables
  recsum = 0.0
  if (tpi2.le.1) evec(6) = evec(6) + ewcnst
  ewtmpf(:) = 0.0
  Qew1(:,:,sta2:sto2) = 0.0
!
! spread scaled coordinates on grid by incrementing the spline array
#ifdef ENABLE_THREADS
  if ((thr_dlb(5,1).gt.0).AND.(tpi.gt.0)) then
    if (tpi2.eq.1) thr_dlb(5,2) = thr_dlb(5,2) + 1
    call System_Clock(count=ttimer)
    thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
  end if
! split into two loops, the second adds overhead
  do rs=sta,sto !1,nseq
    do i=1,at(rs)%npol
      ii = at(rs)%pol(i)
      mvec(1:3) = gdzv(1:3)*(atinfo(1:3,ii)-bnd_params(4:6))
      term0 = atinfo(4,ii)*sqrtel
      ewgfls(1:3,ii) = floor(mvec(1:3))
      ts(1:3) = mvec(1:3) - ewgfls(1:3,ii)
!     first get the Mn(xi), i.e., the relevant B-spline values for the grid-points
!     in the vicinity of atom ii: this can be done separately for all axes
!     get the d/dxi (Mn(xi)): note that the grid axes are aligned with the coordinate axes,
!     so only diagonal terms contribute, i.e., d/dyi (Mn(xi)) is obviously zero
!     WARNING: check how this is behaved in generalized (triclinic) boxes
      call cardBspline(ts(1),splor,bspl(:,1,ii),bspld(:,1,ii))
      call cardBspline(ts(2),splor,bspl(:,2,ii),bspld(:,2,ii))
      call cardBspline(ts(3),splor,bspl(:,3,ii),bspld(:,3,ii))
      bspld(1:splor,1,ii) = bspld(1:splor,1,ii)*gdzv(1)
      bspld(1:splor,2,ii) = bspld(1:splor,2,ii)*gdzv(2)
      bspld(1:splor,3,ii) = bspld(1:splor,3,ii)*gdzv(3)
!
      ewgfls(1:3,ii) = ewgfls(1:3,ii) + gdz(1:3) ! shift to have positive integers
    end do
  end do
  if ((thr_dlb(5,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
  end if
  if (tpi2.gt.0) then
!$OMP BARRIER
  end if
  if ((thr_dlb(5,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
  end if
  do rs=1,nseq
    do i=1,at(rs)%npol
      ii = at(rs)%pol(i)
      term0 = atinfo(4,ii)*sqrtel
      do l=1,splor
        teml = term0*bspl(l,3,ii)
        zi = zmap(ewgfls(3,ii) + l)
        if (zi.lt.sta2) cycle
        if (zi.gt.sto2) cycle
        do k=1,splor
          yi = ymap(ewgfls(2,ii) + k)
          temk = teml*bspl(k,2,ii)
          do j=1,splor
            xi = xmap(ewgfls(1,ii) + j)
            temj = temk*bspl(j,1,ii)
!           intrinsic conversion
            Qew1(xi,yi,zi) = Qew1(xi,yi,zi) + temj
          end do
        end do
      end do
    end do
  end do
#else
  do rs=1,nseq
    do i=1,at(rs)%npol
      ii = at(rs)%pol(i)
      mvec(1:3) = gdzv(1:3)*(atinfo(1:3,ii)-bnd_params(4:6))
      term0 = atinfo(4,ii)*sqrtel
      ewgfls(1:3,ii) = floor(mvec(1:3))
      ts(1:3) = mvec(1:3) - ewgfls(1:3,ii)
!     first get the Mn(xi), i.e., the relevant B-spline values for the grid-points
!     in the vicinity of atom ii: this can be done separately for all axes
!     get the d/dxi (Mn(xi)): note that the grid axes are aligned with the coordinate axes,
!     so only diagonal terms contribute, i.e., d/dyi (Mn(xi)) is obviously zero
!     WARNING: check how this is behaved in generalized (triclinic) boxes
      call cardBspline(ts(1),splor,bspl(:,1,ii),bspld(:,1,ii))
      call cardBspline(ts(2),splor,bspl(:,2,ii),bspld(:,2,ii))
      call cardBspline(ts(3),splor,bspl(:,3,ii),bspld(:,3,ii))
      bspld(1:splor,1,ii) = bspld(1:splor,1,ii)*gdzv(1)
      bspld(1:splor,2,ii) = bspld(1:splor,2,ii)*gdzv(2)
      bspld(1:splor,3,ii) = bspld(1:splor,3,ii)*gdzv(3)
!
      ewgfls(1:3,ii) = ewgfls(1:3,ii) + gdz(1:3) ! shift to have positive integers
      do l=1,splor
        teml = term0*bspl(l,3,ii)
        zi = zmap(ewgfls(3,ii) + l)
        do k=1,splor
          yi = ymap(ewgfls(2,ii) + k)
          temk = teml*bspl(k,2,ii)
          do j=1,splor
            xi = xmap(ewgfls(1,ii) + j)
            temj = temk*bspl(j,1,ii)
!           intrinsic conversion
            Qew1(xi,yi,zi) = Qew1(xi,yi,zi) + temj
          end do
        end do
      end do
    end do
  end do
#endif
#ifdef ENABLE_THREADS
  if ((thr_dlb(5,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
  end if
  if (tpi2.gt.0) then
!$OMP BARRIER
!$OMP MASTER
! now inverse FFT Qarr
!    call System_Clock(count=ttt(1))
    call fftw_execute_dft(fftplanb,Qew1,Qew2)
!    call System_Clock(count=ttt(2))
!    write(*,*) 'Plan B',1.0e-3*(ttt(2)-ttt(1))
!$OMP END MASTER
!$OMP BARRIER
  else
    call fftw_execute_dft(fftplanb,Qew1,Qew2)
  end if
#else
  call fftw_execute_dft(fftplanb,Qew1,Qew2)
#endif
!
! now accumulate energy
! note that the mapping between real and reciprocal space is as follows:
! the edges in the transformed, reciprocal matrix correspond to the nearest images
! in real-space
  do pz=sta2,sto2
    zi = pz - 1
    if (zi.gt.(gdz(3)/2)) zi = zi - gdz(3)
    do py=1,gdz(2)
      yi = py - 1
      if (yi.gt.(gdz(2)/2)) yi = yi - gdz(2)
      zyzero = .false.
      if ((zi.eq.0).AND.(yi.eq.0)) zyzero = .true.
      do px=1,gdz(1)
        xi = px - 1
        if (xi.gt.(gdz(1)/2)) xi = xi - gdz(1)
        if ((zyzero.EQV..true.).AND.(xi.eq.0)) cycle
!        ivec(1) = xi/bnd_params(1)
!        ivec(2) = yi/bnd_params(2)
!        ivec(3) = zi/bnd_params(3)
!        m2 = sum(ivec(:)*ivec(:))
!        temj = bsmx(px)*bsmy(py)*bsmz(pz)
!        term0 = (1.0/m2)*exp(-ewrat*m2)*scale_POLAR
! &                            /(2.0*temj*PI*voli)
        tem =         REAL(Qew2(px,py,pz))*REAL(Qew2(px,py,pz))&
 &                  + AIMAG(Qew2(px,py,pz))*AIMAG(Qew2(px,py,pz))
        recsum = recsum + QewBC(px,py,pz)*tem
        Qew2(px,py,pz) = Qew2(px,py,pz)*2.0*QewBC(px,py,pz)
      end do
    end do
  end do
#ifdef ENABLE_THREADS
  if (tpi2.gt.0) then
!$OMP BARRIER
!$OMP MASTER
! FFT Qew2 into Qew1 (convolution of the approx. structure factors with quasi-array term0)
!    call System_Clock(count=ttt(1))
    call fftw_execute_dft(fftplanf,Qew2,Qew1)
!    call System_Clock(count=ttt(2))
!    write(*,*) 'Plan F',1.0e-3*(ttt(2)-ttt(1))
!$OMP END MASTER
!$OMP BARRIER
  else
    call fftw_execute_dft(fftplanf,Qew2,Qew1)
  end if
#else
  call fftw_execute_dft(fftplanf,Qew2,Qew1)
#endif
!
! now loop over all charges, identify relevant points on dQ/drij and multiply with Qew1
! note that energies are interpolated, not charges, so Newton's 2nd law is broken beyond
! machine precision: the easiest fix is to accumulate the net force and subtract it out
#ifdef ENABLE_THREADS
  if ((thr_dlb(5,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
  end if
#endif
  do rs=sta,sto !1,nseq
    do i=1,at(rs)%npol
      ii = at(rs)%pol(i)
      term0 = atq(ii)*sqrtel
!     get the d/dxi (Mn(xi)): note that the grid axes are aligned with the coordinate axes,
!     so only diagonal terms contribute, i.e., d/dyi (Mn(xi)) is obviously zero
!     WARNING: check how this is behaved in generalized (triclinic) boxes
      catmpf(:) = 0.0
      do cc=1,3
        bsplbu(1:splor) = bspl(1:splor,cc,ii)
        bspl(1:splor,cc,ii) = bspld(1:splor,cc,ii)
!       also note, however, that the values in Q originally are products of spline values, hence
!       the derivative array dQ/dxi has as many terms in it as the array Q itself with just atom
!       i populating it
        do l=1,splor
          temj = term0*bspl(l,3,ii)
          zi = zmap(ewgfls(3,ii) + l)
          do k=1,splor
            yi = ymap(ewgfls(2,ii) + k)
            temk = temj*bspl(k,2,ii)
            do j=1,splor
              xi = xmap(ewgfls(1,ii) + j)
              teml = Qew1(xi,yi,zi)*temk*bspl(j,1,ii)
!             intrinsic conversion
              ewtmpf(cc) = ewtmpf(cc) - teml
              catmpf(cc) = catmpf(cc) - teml
!              write(*,*) AIMAG(Qarr(xi,yi,zi)*teml),
! &                         REAL(Qarr(xi,yi,zi)*teml)
            end do
          end do
        end do
        bspl(1:splor,cc,ii) = bsplbu(1:splor)
      end do
      ca_f(:,ii) = ca_f(:,ii) + catmpf(:)
    end do
  end do
#ifdef ENABLE_THREADS
  if ((thr_dlb(5,1).gt.0).AND.(tpi.gt.0)) then
    call System_Clock(count=ttimer)
    thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
  end if
  if (tpi2.gt.0) then
!$OMP CRITICAL(EWF_UPDATE)
    ewnetf(:) = ewnetf(:) + ewtmpf(:)
!$OMP END CRITICAL(EWF_UPDATE)
!$OMP BARRIER
  else
    ewnetf(:) = ewtmpf(:)
  end if
#else
  ewnetf(:) = ewtmpf(:)
#endif
!
  do rs=sta,sto !1,nseq
    do i=1,at(rs)%npol
      ii = at(rs)%pol(i)
      ca_f(:,ii) = ca_f(:,ii) - ewnetf(:)*ewinvm
    end do
  end do
!
  evec(6) = evec(6) + recsum
!
#ifdef ENABLE_THREADS
  if (tpi2.le.1) then
    call System_Clock(t2)
    time_pme = time_pme + 1.0*(t2-t1)
    time_energy = time_energy + 1.0*(t1-t2)
  end if
#else
  call System_Clock(t2)
  time_pme = time_pme + 1.0*(t2-t1)
  time_energy = time_energy + 1.0*(t1-t2)
#endif
!  write(*,*) 'PME: ',t2-t1
#else
!
  use energies
  use atoms
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: tpi2
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  write(ilog,*) 'Fatal. Code is not compiled against FFTW. PME metho&
 &d unavailable.'
  call fexit()
#endif
!
end
!
!-----------------------------------------------------------------------
!
! this routine shall identify the 1D-spline values for a real number
! between zero and the interpolation order, splor
!
! splor=1 is the square hat (binning) function (M_1)
! splor=2 is linear interpolation: M_2(rrr) = 1 - |rrr - 1| over 0 < rrr < 2 only (linear hat), can be written as
!         the sum of two terms that are products of M_1 functions with the fractional coordinate
! splor=3 is a quadratic function available recursively as the sum of two terms that are products of M_2 functions with
!         the fractional coordinate
! splor=4 is a cubic function defined recursively as the sum of two terms that are products of M_3 functions with
!         the fractional coordinate
! splor=5+ continues analogously
!
subroutine cardBspline(rrr,splor,res,dres)
!
  implicit none
!
  integer bi,splor,k
  RTYPE rrr,u,res(splor),inv,dres(splor)
!
  u = rrr
  if ((u.lt.0.0).OR.(u.gt.splor))  then
    res(:) = 0.0
    dres(:) = 0.0
    return
  end if
!
  if (splor.eq.1) then
    res(1) = 1.0
    dres(1) = 0.0
    return
  end if
!
  res(1) = 1.0 - abs(u)
  res(2) = 1.0 - abs(u-1.0)
  do k=3,splor-1
    inv = 1.0/dble(1.0*k-1.0)
    res(k) = inv*u*res(k-1)
    do bi=1,k-2
      res(k-bi) = inv*&
 & ((u+dble(bi))*res(k-bi-1) + (dble(k-bi)-u)*res(k-bi))
    end do
    res(1) = inv*(1.0-u)*res(1)
  end do
!
  if (splor.gt.2) then
    res(splor) = 0.0
    dres(1) = -res(1)
    do k=2,splor
      dres(k) = res(k-1) - res(k)
    end do
    do k=splor,splor
      inv = 1.0/dble(1.0*k-1.0)
      res(k) = inv*u*res(k-1)
      do bi=1,k-2
        res(k-bi) = inv*&
 & ((u+dble(bi))*res(k-bi-1) + (dble(k-bi)-u)*res(k-bi))
      end do
      res(1) = inv*(1.0-u)*res(1)
    end do
  else
    dres(1) = -1.0
    dres(2) = 1.0
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine get_bsplmods(splor,xd,yd,zd,bsmx,bsmy,bsmz)
!
  use math
!
  implicit none
!
  integer xi,yi,zi,xd,yd,zd,splor,k
  RTYPE bsmx(xd),bsmy(yd),bsmz(zd),splx(splor),azero,dum(splor)
  complex(KIND=8) prefac,sm1,acc
  PARAMETER (sm1 = (0,1))
!
  azero = 0.0
  call cardBspline(azero,splor,splx,dum)
!
  prefac = 2.0*PI*sm1/dble(xd)
  do xi=1,xd
    acc = 0.0
    do k=0,splor-2
      acc = acc + splx(k+1)*exp(prefac*dble((xi-1)*k))
    end do
    bsmx(xi) = REAL(acc)*REAL(acc) + AIMAG(acc)*AIMAG(acc)
  end do
!
  prefac = 2.0*PI*sm1/dble(yd)
  do yi=1,yd
    acc = 0.0
    do k=0,splor-2
      acc = acc + splx(k+1)*exp(prefac*dble((yi-1)*k))
    end do
    bsmy(yi) = REAL(acc)*REAL(acc) + AIMAG(acc)*AIMAG(acc)
  end do
!
  prefac = 2.0*PI*sm1/dble(zd)
  do zi=1,zd
    acc = 0.0
    do k=0,splor-2
      acc = acc + splx(k+1)*exp(prefac*dble((zi-1)*k))
    end do
    bsmz(zi) = REAL(acc)*REAL(acc) + AIMAG(acc)*AIMAG(acc)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this routine serves the following tasks:
! 1) get a reasonable grid size for the system at hand
! 2) initialize and store the constant term
! 3) allocate and initialize the grid itself
! 4) allocate and store the B-spline moduli for the fixed grid
!
subroutine setup_pme_constV()
!
  use iounit
  use atoms
  use ewalds
  use cutoffs
  use polypep
  use system
  use sequen
  use units
  use math
  use energies
#ifdef ENABLE_THREADS
  use threads
#endif
#ifdef LINK_FFTW
  use, intrinsic :: ISO_C_BINDING
#endif
!
  implicit none
!
#ifdef LINK_FFTW
!
#ifdef ORACLE_FORTRAN
#include "fftw3_SUN.f03"
#else
#include "fftw3.f03"
#endif
!
  integer(C_INT):: errcode
  character(KIND=C_CHAR):: stringy(0:MAXSTRLEN)
  integer i,j,k,rs,gdz(3),xi,yi,zi,px,py,pz,t1,t2
  RTYPE netQ,voli,gsz,der,def,ivec(3),temj,m2,ewrat
  RTYPE memreq,phiterm,Hterm,spfac,slen,area,QN
#ifdef ENABLE_THREADS
  integer threadsok
#endif
!
  if ((bnd_type.ne.1).OR.(bnd_shape.ne.1)) then
    write(ilog,*) 'Fatal. Ewald sums are not supported for non-3D-periodic systems.'
    call fexit()
  end if
!
  do j=1,3
    kdims(j,1) = nint(bnd_params(j)/ewfspac)
    kdims(j,2) = 2*kdims(j,1)
  end do
  gdz(:) = kdims(:,2)
  gsz = min(gdz(1),min(gdz(2),gdz(3)))
!
  if (gsz.le.splor) then
    write(ilog,*) 'Warning. Decreasing Fourier spacing for PME due to conflict with B-spline order.'
  end if
  do while (gsz.le.splor)
    ewfspac = ewfspac/2.0
    do j=1,3
      kdims(j,1) = nint(bnd_params(j)/ewfspac)
      kdims(j,2) = 2*kdims(j,1)
    end do
    gdz(:) = kdims(:,2)
    gsz = min(gdz(1),min(gdz(2),gdz(3)))
  end do
  Hterm = gsz
  k = 1
  do while (Hterm.gt.2.0)
    k = k + 1
    Hterm = Hterm/2.0
  end do
! not clear what is the optimal form of Hterm to use
!  Hterm = 2.0**(-1.0*k)
!
  memreq = dble(gdz(1))*dble(gdz(2))*dble(gdz(3))
  if (memreq.gt.1.0e8) then
    write(ilog,*) 'Warning. Mesh in PME is exceptionally large. Cons&
 &ider increasing both FMCSC_EWFSPAC and FMCSC_BSPLINE. This run mig&
 &ht swap and/or crash.'
  end if
!
  voli = bnd_params(1)*bnd_params(2)*bnd_params(3)
  slen = voli**(1.0/3.0)
  area = voli**(2.0/3.0)
! simplified but may work better (more conservative)
  Hterm = ewfspac/slen
!
  netQ = 0.0
  QN = 0.0
  do rs=1,nseq
    do i=1,at(rs)%npol
      QN = QN + 1.0
      netQ = netQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
    end do
  end do
!
! WARNING: the search procedure is broken!!!!!!!!!!!!
!          (bad formulas for accuracy/tolerance)
  spfac = dble(splor+1)
  do i=splor,1,-1
    spfac = spfac*dble(i)
  end do
  phiterm = 1.0/spfac
! use user-provided Ewald-parameter
  if (ewpm_pre.gt.0.0) then
    ewpm = ewpm_pre
    ewetol = huge(ewetol)
    der = (1.0/sqrt(QN*mcnb_cutoff*voli))*&
 &    exp(-ewpm*ewpm*mcnb_cutoff2)*netQ*electric*2.0
!   reciprocal-space part
    def = 2.0*(PI**0.25)*(netQ/area)*phiterm*((2.0*ewpm*slen*Hterm)**(splor+1.0))*&
 &         exp(0.5*dble(splor+1)*(log(splor+1.0) - log(2.0) + 1.0))*&
 &         sqrt(6.0*ewpm*slen/(QN*(2.0*splor+3.0)))
    ewetol = QN*sqrt(der*der+def*def)
! scan reasonable regime for first minimum
  else
    ewpm = 0.0
    ewetol = huge(ewetol)
    do while (ewpm.lt.1.0)
      ewpm = ewpm + 0.01
!     real-space part
      der = (1.0/sqrt(mcnb_cutoff*voli))*&
 &    exp(-ewpm*ewpm*mcnb_cutoff2)*netQ*electric*2.0
!     reciprocal-space part
      def = 2.0*(PI**0.25)*(netQ/area)*phiterm*((2.0*ewpm*slen*Hterm)**(splor+1.0))*&
 &         exp(0.5*dble(splor+1)*(log(splor+1.0) - log(2.0) + 1.0))*&
 &         sqrt(6.0*ewpm*slen/(QN*(2.0*splor+3.0)))
      if (QN*sqrt(der*der+def*def).gt.ewetol) then
        ewpm = ewpm - 0.01
        exit
      end if
      ewetol = QN*sqrt(der*der+def*def)
    end do
  end if
!
!  ewpm = 0.1
  ewpm2 = ewpm*ewpm
  ewpite = 2.0*ewpm/sqrt(PI)
  ewcnst = -scale_POLAR*electric*(ewpm/sqrt(PI))*netQ
!
  allocate(bsmx(gdz(1)))
  allocate(bsmy(gdz(2)))
  allocate(bsmz(gdz(3)))
  call get_bsplmods(splor,gdz(1),gdz(2),gdz(3),bsmx,bsmy,bsmz)
!
  nqew1 = gdz(1)*gdz(2)*gdz(3)
  nqew2 = nqew1
  pqew1 = fftw_alloc_complex(nqew1)
  call C_F_POINTER(pqew1, Qew1, gdz)
  pqew2 = fftw_alloc_complex(nqew2)
  call C_F_POINTER(pqew2, Qew2, gdz)
!  allocate(Qew1(gdz(1),gdz(2),gdz(3)))
!  allocate(Qew2(gdz(1),gdz(2),gdz(3)))
!
#ifdef ENABLE_THREADS
! first do thread initialization for FFTW
  threadsok = fftw_init_threads()
  if (threadsok.eq.0) then
    write(ilog,*) 'Fatal. Unresolved error occurred during thread initialization for FFTW.'
    call fexit()
  end if
! then provide max number of concurrent threads for later use
! note that the only threadsafe function are variants of fftw_execute
  call fftw_plan_with_nthreads(thrdat%maxn)
#endif
!
! should yield the same plan in terms of algorithm, but different pointers ...
! order of array indeces must be reversed
  if (fftlevel.eq.1) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_ESTIMATE)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_ESTIMATE)
  else if (fftlevel.eq.2) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_MEASURE)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_MEASURE)
  else if (fftlevel.eq.3) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_PATIENT)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_PATIENT)
  else if (fftlevel.eq.4) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  else if (fftlevel.eq.5) then
    i = MAXSTRLEN
    call strlims(fftfname,t1,t2)
    do i=t1,t2
      stringy(i-t1) = fftfname(i:i)
    end do
    stringy((t2-t1+1)) = C_NULL_CHAR
    errcode = fftw_import_wisdom_from_filename(stringy);
    if (errcode.eq.0) then
      write(ilog,*) 'Fatal. FFTW wisdom could not be read (got ',fftfname(t1:t2),' as filename).'
      call fexit()
    end if
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_WISDOM_ONLY)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_WISDOM_ONLY)
  end if
  if ((fftprint.eq.1).AND.(fftlevel.ne.5)) then
    i = MAXSTRLEN
    call strlims(fftfname,t1,t2)
    do i=t1,t2
      stringy(i-t1) = fftfname(i:i)
    end do
    stringy((t2-t1+1)) = C_NULL_CHAR
    errcode = fftw_export_wisdom_to_filename(stringy);
  end if
!
! populate static grid variables: perfect for constant volume sim.s
  allocate(QewBC(gdz(1),gdz(2),gdz(3)))
  ewrat = PI*PI/(ewpm*ewpm)
  do px=1,gdz(1)
    do py=1,gdz(2)
      do pz=1,gdz(3)
        xi = px - 1
        yi = py - 1
        zi = pz - 1
        if (xi.gt.(gdz(1)/2)) xi = xi - gdz(1)
        if (yi.gt.(gdz(2)/2)) yi = yi - gdz(2)
        if (zi.gt.(gdz(3)/2)) zi = zi - gdz(3)
        if ((xi.eq.0).AND.(yi.eq.0).AND.(zi.eq.0)) then
          QewBC(px,py,pz) = 0.0
        end if
        ivec(1) = dble(xi)/bnd_params(1)
        ivec(2) = dble(yi)/bnd_params(2)
        ivec(3) = dble(zi)/bnd_params(3)
        m2 = sum(ivec(:)*ivec(:))
        temj = bsmx(px)*bsmy(py)*bsmz(pz)
        if ((xi.eq.0).AND.(yi.eq.0).AND.(zi.eq.0)) cycle
        QewBC(px,py,pz) = (1.0/m2)*exp(-ewrat*m2)*scale_POLAR&
 &                            /(2.0*temj*PI*voli)
      end do
    end do
  end do
!
! allocate B-spline parameters and grid associations
  allocate(bspl(splor,3,n))
  allocate(bspld(splor,3,n))
  allocate(ewgfls(3,n))
!
! get inverse of total number of charges
  k = 0
  do rs=1,nseq
    do i=1,at(rs)%npol
      k = k + 1
    end do
  end do
  ewinvm = 1.0/dble(k)
!
#else
  write(ilog,*) 'Fatal. Code is not compiled against FFTW. PME metho&
 &d unavailable.'
  call fexit()
#endif
!
end
!
!-------------------------------------------------------------------------------------
!
#ifdef LINK_FFTW
#ifdef ENABLE_THREADS
!
! redo FFTW plans with a specified max number of threads
!
subroutine pme_plan_threads(ntt)
!
  use ewalds
  use iounit
  use, intrinsic :: ISO_C_BINDING
!
  implicit none
!
#ifdef ORACLE_FORTRAN
#include "fftw3_SUN.f03"
#else
#include "fftw3.f03"
#endif
!
  integer, INTENT(IN):: ntt
!
  integer gdz(3)
!
  gdz(:) = kdims(:,2)
!
  if (fftlevel.eq.5) then
    write(ilog,*) 'Warning. FFTW plan creation using wisdom is disable in tests of thread efficiency &
 &or when changing numbers of dedicated nodes. Using FFTW_PATIENT instead. This may take a while.'
    fftlevel = 3
  end if
!
  call fftw_destroy_plan(fftplanf)
  call fftw_destroy_plan(fftplanb)
  call fftw_cleanup_threads()
  call fftw_plan_with_nthreads(ntt)
!
  if (fftlevel.eq.1) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_ESTIMATE)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_ESTIMATE)
  else if (fftlevel.eq.2) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_MEASURE)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_MEASURE)
  else if (fftlevel.eq.3) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_PATIENT)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_PATIENT)
  else if (fftlevel.eq.4) then
    fftplanb = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew1,Qew2,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
    fftplanf = fftw_plan_dft_3d(gdz(3),gdz(2),gdz(1),Qew2,Qew1,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  end if
!
end
!
#endif
#endif
!
!-----------------------------------------------------------------------
!
! this routine serves the following tasks:
! 1) get a reasonable parameter set for the system at hand
! 2) initialize and store the constant term
!
! note that in the interest of the user the short-range (real-space) cutoff needs to be
! fixed, which leaves three options:
! i)   user-fix reciprocal cutoff as well and determine ewpm for max. accuracy
! ii)  user-fix ewpm and determine reciprocal cutoff for max. accuracy
! iii) user-fix accuracy and determine reciprocal cutoff as well as ewpm optimally
!
subroutine setup_ewald_constV()
!
  use iounit
  use atoms
  use ewalds
  use cutoffs
  use polypep
  use system
  use sequen
  use units
  use math
  use energies
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer i,k,rs,QN,dimv(3)
  RTYPE netQ,voli,gsz,der,def,area,lend
!
  if ((bnd_type.ne.1).OR.(bnd_shape.ne.1)) then
    write(ilog,*) 'Fatal. Ewald sums are not supported for non-3D-periodic systems.'
    call fexit()
  end if
!
! set the reciprocal space cutoffs according to desired spacing
  do k=1,3
    kdims(k,1) = -ceiling(bnd_params(k)/ewfspac)
    kdims(k,2) = ceiling(bnd_params(k)/ewfspac)
  end do
  dimv(:) = kdims(1:3,2)-kdims(1:3,1)+1
#ifdef ENABLE_THREADS
  allocate(ewSfac(dimv(1),dimv(2),dimv(3),thrdat%maxn))
  allocate(ewSfacs(dimv(1),dimv(2),dimv(3),thrdat%maxn))
#else
  allocate(ewSfac(dimv(1),dimv(2),dimv(3),1))
  allocate(ewSfacs(dimv(1),dimv(2),dimv(3),1))
#endif
!  allocate(ewSvd(3,n))
!  allocate(ewSvds(3,n))
! get net squared charge
  netQ = 0.0
  QN = 0
  do rs=1,nseq
    do i=1,at(rs)%npol
      QN = QN + 1
      netQ = netQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
    end do
  end do
! get system volume and find dimension of lowest resolution
  voli = bnd_params(1)*bnd_params(2)*bnd_params(3)
  gsz = huge(gsz)
  do k=1,3
    if ((kdims(k,2)-kdims(k,1)+1).lt.gsz) then
      gsz = kdims(k,2)-kdims(k,1)+1
    end if
  end do
  gsz = (gsz-1)/2
  area = voli**(2./3.)
  lend = voli**(1./3.)
! use user-provided Ewald parameter
  if (ewpm_pre.gt.0.0) then
    ewpm = ewpm_pre
    der = (1.0/sqrt(QN*mcnb_cutoff*voli))*&
 &    exp(-ewpm*ewpm*mcnb_cutoff2)*netQ*electric*2.0
    def= sqrt(8.0/gsz)*exp((-(PI*gsz/ewpm)**2)/area)*netQ*electric/(PI*lend)
    ewetol = QN*sqrt(der*der+def*def) !der+def
! determine Ewald parameter by scanning reasonable regime
  else
    ewpm = 0.0
    ewetol = HUGE(ewetol)
    do while (ewpm.lt.1.0)
      ewpm = ewpm + 0.01
      der = (1.0/sqrt(QN*mcnb_cutoff*voli))*&
 &    exp(-ewpm*ewpm*mcnb_cutoff2)*netQ*electric*2.0
      def= sqrt(8.0/(1.0*QN*gsz))*exp((-(1.0/area)*(PI*gsz/ewpm)**2))*netQ*electric*ewpm/(PI*lend)
      if (QN*sqrt(der*der+def*def).gt.ewetol) then
        ewpm = ewpm - 0.01
        exit
      end if
      ewetol = QN*sqrt(der*der+def*def)
    end do
  end if
!
  ewpm2 = ewpm*ewpm
  ewpite = 2.0*ewpm/sqrt(PI)
  ewcnst = -scale_POLAR*electric*(ewpm/sqrt(PI))*netQ
!  write(*,*) netQ,ewcnst,ewetol
!
end
!
!-----------------------------------------------------------------------
!
! right now this routine eliminates calculations rather indiscriminantly
! adjust if support increases
!
subroutine setup_ewald()
!
  use cutoffs
  use ewalds
  use iounit
  use system
  use energies
  use sequen
!
  implicit none
!
  RTYPE minsdl,maxrrd
  integer i
!
  if (use_POLAR.EQV..false.) then
    lrel_md = 1
    return
  end if
!
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    write(ilog,*) 'Fatal. Ewald summation is currently not supported&
 & for ensembles with fluctuating volumes. Check back later.'
    call fexit()
  end if
!
  if (use_FEG.EQV..true.) then
    write(ilog,*) 'Fatal. Ewald summation is currently not supported&
 & for free energy growth calculations. Check back later.'
    call fexit()
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    write(ilog,*) 'Fatal. Ewald summation is currently not supported&
 & for the ABSINTH implicit solvent model. Check back later.'
    call fexit()
  end if
!
  if (is_pewlj.EQV..false.) then
    write(ilog,*) 'Fatal. Current Hamiltonian does not support Ewald&
 & summation. Please check back later.'
    call fexit()
  end if
!
  if ((bnd_type.ne.1).OR.(bnd_shape.ne.1)) then
    write(ilog,*) 'Fatal. Ewald sums are not supported for non-3D-periodic systems.'
    call fexit()
  end if
!
  if (bnd_shape.ne.1) then
    write(ilog,*) 'Fatal. Ewald summation is currently not supported&
 & for non-rectangular, periodic boxes. Check back later.'
    call fexit()
  end if
!
  if (use_cutoffs.EQV..false.) then
    write(ilog,*) 'Fatal. Ewald summation requires using and setting&
 &a real-space finite cutoff.'
    call fexit()
  else
    mcel_cutoff = mcnb_cutoff
    mcel_cutoff2 = mcnb_cutoff2
    minsdl = huge(minsdl)
    do i=1,3
      if (bnd_params(i).lt.minsdl) then
        minsdl = bnd_params(i)
      end if
    end do
    maxrrd = tiny(maxrrd)
    do i=1,nseq
      if (resrad(i).gt.maxrrd) maxrrd = resrad(i)
    end do
    if ((mcnb_cutoff+2.0*maxrrd).gt.(0.5*minsdl)) then
      write(ilog,*) 'Fatal. Ewald sums are technically incompatible &
 &with the given cutoff size. Residue-wise shift vectors have to ap&
 &ply for all atom pairs in a residue-pair. Increase system size or&
 & decrease cutoff.'
      call fexit()
    end if
  end if
!
  if (ewald_mode.eq.1) then
    call setup_pme_constV()
  else if (ewald_mode.eq.2) then
    call setup_ewald_constV()
  else
    write(ilog,*) 'Fatal. Encountered unsupported mode for Ewald sum&
 &s (offending mode is ',ewald_mode,'). Please report this bug.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
