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
! MAIN AUTHOR:   Hoang Tran                                                !
! CONTRIBUTIONS: Rohit Pappu, Andreas Vitalis                              !
!                                                                          !
!--------------------------------------------------------------------------!
!
!
#include "macros.i"
#define VC_SIZE 8
#define VC_SIZEM1 7
!
! ###################################################
! ##                                               ##
! ##  subroutine readgrid -- read in phi-psi grids ##
! ##                                               ##
! ###################################################
!
! "readgrid" reads in a predetermined grid of important phi,psi
! values into memory for later use in Markov chain sampling
!
subroutine readgrid()
!
  use grids
  use iounit
  use aminos
!
  implicit none
!
  integer funit,freeunit,next,i,t1,t2
  RTYPE ff,yy
  character(60) fyfile
  character(80) record
  character(3) aa3lc
  logical exists
!
  call strlims(griddir,t1,t2)
  write(ilog,*)
  write(ilog,*) 'Using grids in ',griddir(t1:t2)
!
  do i = 1,MAXAMINO
     aa3lc = amino(i)
     call tolower(aa3lc)
     fyfile=griddir(t1:t2)//aa3lc//'_grid.dat'
     inquire(file=fyfile,exist=exists)
     if (exists.EQV..false.) then
!
!       read in default grid
        fyfile=griddir(t1:t2)//'def_grid.dat'
     end if
     next = 1
!    
! read in grid centers from fyfile
     inquire(file=fyfile,exist=exists)
     if(exists.EQV..false.) then 
        call fexit()
     else
        funit=freeunit()
        open(unit=funit,file=fyfile,status='old')
        do while (.true.)
           read(funit,10,end=30) record
 10            format(a80)
           read(record,*,err=20,end=20)ff,yy
 20            continue
           stgr%it(i,next,1) = ff
           stgr%it(i,next,2) = yy
           next = next + 1
        end do
 30         continue
        next = next - 1 
        stgr%ngr(i) = next
        close(unit=funit)
     end if
      
     write(ilog,35) stgr%ngr(i),aa3lc//'_grid.dat',amino(i)
 35      format('     Read in [',i5,'] grid centers from: ',a20,&
 &          ' for ',a3)
  end do
!
  stgr%sfyc(1) = stgr%halfwindow
  write(ilog,*)
  do i = 1,MAXAMINO
     stgr%sfyc(i) = stgr%sfyc(1)
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_example(evec,rs1,ca_f,hira,hira2,which)
!
  use energies
  use polypep
  use forces
  use atoms
  use params
  use units
  use cutoffs
  use sequen
  use system
  use molecule
  use ewalds, ONLY: ewpite,ewpm,ewpm2
  use tabpot, ONLY: p_terfcoverr,i_terfcoverr,terfcoverr
!
  implicit none
!
  integer, INTENT(IN):: rs1,hira,hira2,which
!
  integer hira3,ii,j,hi,lo,k,i,jj,imol,chnksz,chnk,kkk
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvecx(VC_SIZE),dvecy(VC_SIZE),dvecz(VC_SIZE)
  RTYPE foix(VC_SIZE),foiy(VC_SIZE),foiz(VC_SIZE)
  RTYPE dvecix(hira2),dveciy(hira2),dveciz(hira2)
  RTYPE for_x(hira2),for_y(hira2),for_z(hira2)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na)
  RTYPE term0(hira2)
  integer k1(hira),k2(hira),ka(hira)
  RTYPE drax(hira),dray(hira),draz(hira),srax(hira),sray(hira),sraz(hira)
  RTYPE lop1,lop2,lop3,lop4,epsik,radik
!  
! initialize
  imol = molofrs(rs1)
  if (which.eq.1) then
    hira3 = 3*rs_nbl(rs1)%nwnbs
  else if (which.eq.2) then
    hira3 = 3*rs_nbl(rs1)%nwnbtrs
  else
    return
  end if
  for_x(1:hira3) = 0.0
  for_y(1:hira3) = 0.0
  for_z(1:hira3) = 0.0
!
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+at(rs1)%na-1
    j = ii-rsinfo(rs1,1)+1
    k0(j) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
  end do
!
  if (which.eq.1) then ! standard range, recompute shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnb(1:hira)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    drax(1:hira) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1))
    dray(1:hira) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1))
    draz(1:hira) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1))
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        lop1 = 1.0/bnd_params(1)
        lop3 = -bnd_params(1)
        srax(1:hira) = lop3*anint(drax(1:hira)*lop1)
        lop1 = 1.0/bnd_params(2)
        lop3 = -bnd_params(2)
        sray(1:hira) = lop3*anint(dray(1:hira)*lop1)
        lop1 = 1.0/bnd_params(3)
        lop3 = -bnd_params(3)
        sraz(1:hira) = lop3*anint(draz(1:hira)*lop1)
      else if (bnd_shape.eq.3) then
        srax(1:hira) = 0.0
        sray(1:hira) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        sraz(1:hira) = lop3*anint(draz(1:hira)*lop1)
      end if
    else 
      srax(:) = 0.0
      sray(:) = 0.0
      sraz(:) = 0.0
    end if
    drax(1:hira) = drax(1:hira) + srax(1:hira)
    dray(1:hira) = dray(1:hira) + sray(1:hira)
    draz(1:hira) = draz(1:hira) + sraz(1:hira)
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(1:hira)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    srax(1:hira) = rs_nbl(rs1)%wtrsvec(1:hira,1)
    sray(1:hira) = rs_nbl(rs1)%wtrsvec(1:hira,2)
    sraz(1:hira) = rs_nbl(rs1)%wtrsvec(1:hira,3)
    drax(1:hira) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1)) + srax(1:hira)
    dray(1:hira) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1)) + sray(1:hira)
    draz(1:hira) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1)) + sraz(1:hira)
  end if
!
! we can skip shift vector checks here
  ka(1) = 1
  do k=2,hira
    ka(k) = ka(k-1)+3
  end do
  do k=1,hira
    do j=k1(k)+1,k1(k)+3
      jj = ka(k) + j - k1(k) - 1
      dvecix(jj) = srax(k) + atinfo(1,j)
      dveciy(jj) = sray(k) + atinfo(2,j)
      dveciz(jj) = sraz(k) + atinfo(3,j)
      term0(jj) = atinfo(4,j)
    end do
  end do
  if (jj.ne.hira3) call fexit()
!
  kkk = b_type(rsinfo(rs1,1))
  radik = lj_sig(bio_ljtyp(kkk),bio_ljtyp(kkk))
  epsik = lj_eps(bio_ljtyp(kkk),bio_ljtyp(kkk))
  lop1 = 4.0*scale_IPP*epsik
  lop2 = 4.0*scale_attLJ*epsik
  lop3 = electric*scale_POLAR
!
! LJ
  do chnk=1,hira,VC_SIZE
    hi = min(chnk+VC_SIZEM1,hira)
    chnksz = hi-chnk+1
    lo = chnk + hira3
    ii = hi + hira3
    foix(1:VC_SIZE) = 1.0/((drax(chnk:hi))**2 + (dray(chnk:hi))**2 + (draz(chnk:hi))**2)
    term6(1:VC_SIZE) = (radik*foix(1:VC_SIZE))**3
    foiz(1:VC_SIZE) = lop1*term6(1:VC_SIZE)*term6(1:VC_SIZE)
    foiy(1:VC_SIZE) = lop2*term6(1:VC_SIZE)
    evec(1) = evec(1) + sum(foiz(1:VC_SIZE))
    evec(3) = evec(3) - sum(foiy(1:VC_SIZE))
    term6(1:VC_SIZE) = foix(1:VC_SIZE)*(12.0*foiz(1:VC_SIZE)-6.0*foiy(1:VC_SIZE))
!   force
    for_x(lo:ii) = drax(chnk:hi)*term6(1:VC_SIZE)
    for_y(lo:ii) = dray(chnk:hi)*term6(1:VC_SIZE)
    for_z(lo:ii) = draz(chnk:hi)*term6(1:VC_SIZE)
    ca_f(1,k0(1)) = ca_f(1,k0(1)) - sum(for_x(lo:ii))
    ca_f(2,k0(1)) = ca_f(2,k0(1)) - sum(for_y(lo:ii))
    ca_f(3,k0(1)) = ca_f(3,k0(1)) - sum(for_z(lo:ii))
  end do
!
! Coulomb
  if ((is_plj.EQV..true.).OR.(is_fegplj.EQV..true.)) then
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvecx(1:chnksz) = dvecix(chnk:hi)
      dvecy(1:chnksz) = dveciy(chnk:hi)
      dvecz(1:chnksz) = dveciz(chnk:hi)
      do i=2,at(rs1)%na
        lop4 = lop3*ixis(4,i)
        foix(1:chnksz) = 1.0/((dvecx(1:chnksz)-ixis(1,i))**2 + &
 &                            (dvecy(1:chnksz)-ixis(2,i))**2 + &
 &                            (dvecz(1:chnksz)-ixis(3,i))**2)
!       Coulomb potential
        foiy(1:chnksz) = lop4*xt0(1:chnksz)*sqrt(foix(1:chnksz))
        evec(6) = evec(6) + sum(foiy(1:chnksz))
        term6(1:chnksz) = foix(1:chnksz)*foiy(1:chnksz)
!       force 
        foix(1:chnksz) = (dvecx(1:chnksz)-ixis(1,i))*term6(1:chnksz)
        foiy(1:chnksz) = (dvecy(1:chnksz)-ixis(2,i))*term6(1:chnksz)
        foiz(1:chnksz) = (dvecz(1:chnksz)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foix(1:chnksz))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foiy(1:chnksz))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foiz(1:chnksz))
        for_x(chnk:hi) = for_x(chnk:hi) + foix(1:chnksz)
        for_y(chnk:hi) = for_y(chnk:hi) + foiy(1:chnksz)
        for_z(chnk:hi) = for_z(chnk:hi) + foiz(1:chnksz)


!       not faster with GNU or Intel or Sun on SSE4
!        lop4 = lop3*ixis(4,i)
!        xt1(1:chnksz) = dvec(1:chnksz,1)-ixis(1,i)
!        xt2(1:chnksz) = dvec(1:chnksz,2)-ixis(2,i)
!        xt3(1:chnksz) = dvec(1:chnksz,3)-ixis(3,i)

!        xt4(1:chnksz) = xt1(1:chnksz)**2 + xt2(1:chnksz)**2 + xt3(1:chnksz)**2 
!        foin(1:chnksz,1) = sqrt(1.0/xt4(1:chnksz))
!        evec(6) = evec(6) + sum(lop4*xt0(1:chnksz)*foin(1:chnksz,1))
!        term6(1:chnksz) = lop4*xt0(1:chnksz)*foin(1:chnksz,1)*foin(1:chnksz,1)*foin(1:chnksz,1)
!        xt1(1:chnksz) = xt1(1:chnksz)*term6(1:chnksz)
!        xt2(1:chnksz) = xt2(1:chnksz)*term6(1:chnksz)
!        xt3(1:chnksz) = xt3(1:chnksz)*term6(1:chnksz)
!        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(xt1(1:chnksz))
!        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(xt2(1:chnksz))
!        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(xt3(1:chnksz))
!        for_k(chnk:hi,1) = for_k(chnk:hi,1) + xt1(1:chnksz)
!        for_k(chnk:hi,2) = for_k(chnk:hi,2) + xt2(1:chnksz)
!        for_k(chnk:hi,3) = for_k(chnk:hi,3) + xt3(1:chnksz)


      end do
    end do
  end if
!
  k1(1:hira) = k1(1:hira) + 1
  do k=1,hira
    ca_f(1,k1(k)-1) = ca_f(1,k1(k)-1) + for_x(k+hira3)
    ca_f(2,k1(k)-1) = ca_f(2,k1(k)-1) + for_y(k+hira3)
    ca_f(3,k1(k)-1) = ca_f(3,k1(k)-1) + for_z(k+hira3)
    ca_f(1,k1(k):(k1(k)+2)) = ca_f(1,k1(k):(k1(k)+2)) + for_x(ka(k):(ka(k)+2))
    ca_f(2,k1(k):(k1(k)+2)) = ca_f(2,k1(k):(k1(k)+2)) + for_y(ka(k):(ka(k)+2))
    ca_f(3,k1(k):(k1(k)+2)) = ca_f(3,k1(k):(k1(k)+2)) + for_z(ka(k):(ka(k)+2))
  end do
!
end

