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
#define VC_SIZE 48
#define VC_SIZEM1 47
!
!-----------------------------------------------------------------------
!
!
!             ####################################
!             #                                  #
!             #  A SET OF FORCE/ENERGY ROUTINES  #
!             #  WHICH WILL USUALLY TAKE UP THE  #
!             #  BULK OF CPU TIME FOR "VACUUM"   #
!             #  (CONDENSED PHASE) CALCULATIONS  #
!             #                                  #
!             ####################################
!
!
!-----------------------------------------------------------------------
!
! the first two routines do not use cutoffs, and consequently no neighbor lists
! (Vforce_PLJ, Vforce_LJ)
!
subroutine Vforce_PLJ(evec,rs1,ca_f)
!
  use energies
  use polypep
  use atoms
  use params, ONLY: lj_sig,lj_eps
  use units
  use fyoc
  use cutoffs
  use sequen
  use system
!
  implicit none
!
  integer, INTENT(IN):: rs1
!
  integer ii,j,hi,i,shii,ki,kj,imol,chnksz,chnk,rs2
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE foin(VC_SIZE,3),xt1(VC_SIZE)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  RTYPE lop1,lop2,lop3,lop4,lopb(3),lopc(3)
!  
  if (rs1.ge.(nseq-1)) return
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
  end do
  ii = refat(rs1)
!
! reuse parameters
  lop1 = 4.0*scale_IPP
  lop2 = 4.0*scale_attLJ
  lop3 = electric*scale_POLAR
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
      do j=1,3
        lopb(j) = 1.0/bnd_params(j)
        lopc(j) = -bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      lopb(3) = 1.0/bnd_params(6)
      lopc(3) = -bnd_params(6)
    end if
  end if 
!
! gas phase (this routine is not generally safe to use with ghosting)
  if (is_plj.EQV..true.) then
    do chnk=at(rs1+2)%bb(1),n,VC_SIZE
      hi = min(chnk + VC_SIZEM1,n)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = atinfo(1,chnk:hi) 
      dvec(1:chnksz,2) = atinfo(2,chnk:hi)
      dvec(1:chnksz,3) = atinfo(3,chnk:hi)
      xt0(1:chnksz) = atinfo(4,chnk:hi)
      xt1(1:chnksz) = 1.0
      do rs2=atmres(chnk),atmres(hi)
        ki = max(chnk,at(rs2)%bb(1))-chnk+1
        kj = min(hi,at(rs2)%bb(1)+at(rs2)%na-1)-chnk+1
        if (disulf(rs1).eq.rs2) xt1(ki:kj) = 0.0
        if (molofrs(rs2).ne.imol) then
          if (bnd_type.eq.1) then
!           cubic box
            if (bnd_shape.eq.1) then
              do j=1,3
                dvec(ki:kj,j) = dvec(ki:kj,j) + lopc(j)*anint((atinfo(j,refat(rs2))-atinfo(j,refat(rs1)))*lopb(j))
              end do
            else if (bnd_shape.eq.3) then
              dvec(ki:kj,3) = dvec(ki:kj,3) + lopc(3)*anint((atinfo(3,refat(rs2))-atinfo(3,refat(rs1)))*lopb(3))
            end if
          end if
        end if
      end do
      do i=1,at(rs1)%na
!       inverse distance squared and LJ potential
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
        term6(1:chnksz) = xt1(1:chnksz)*(lj_sig(attyp(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
        foin(1:chnksz,2) = lj_eps(attyp(chnk:hi),k0(i,1))
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
        term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        evec(3) = evec(3) - sum(term6(1:chnksz))
!       Coulomb potential
        foin(1:chnksz,2) = lop4*xt1(1:chnksz)*xt0(1:chnksz)*sqrt(foin(1:chnksz,1))
        evec(6) = evec(6) + sum(foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*(foin(1:chnksz,2)-6.0*term6(1:chnksz)+12.0*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        ca_f(:,chnk:hi) = ca_f(:,chnk:hi) + transpose(foin(1:chnksz,1:3))
      end do
    end do
  else
    return
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine Vforce_LJ(evec,rs1,ca_f)
!
  use energies
  use polypep
  use atoms
  use params, ONLY: lj_sig,lj_eps
  use units
  use fyoc
  use cutoffs
  use sequen
  use system
!
  implicit none
!
  integer, INTENT(IN):: rs1
!
  integer ii,j,hi,i,shii,ki,kj,imol,chnksz,chnk,rs2
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE foin(VC_SIZE,3),xt1(VC_SIZE)
  RTYPE ixis(3,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  RTYPE lop1,lop2,lopb(3),lopc(3)
!  
  if (rs1.ge.(nseq-1)) return
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
  end do
  ii = refat(rs1)
!
! reuse parameters
  lop1 = 4.0*scale_IPP
  lop2 = 4.0*scale_attLJ
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
      do j=1,3
        lopb(j) = 1.0/bnd_params(j)
        lopc(j) = -bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      lopb(3) = 1.0/bnd_params(6)
      lopc(3) = -bnd_params(6)
    end if
  end if 
!
! gas phase (this routine is not generally safe to use with ghosting)
  if ((is_lj.EQV..true.).OR.(is_ev.EQV..true.)) then
    do chnk=at(rs1+2)%bb(1),n,VC_SIZE
      hi = min(chnk + VC_SIZEM1,n)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = atinfo(1,chnk:hi) 
      dvec(1:chnksz,2) = atinfo(2,chnk:hi)
      dvec(1:chnksz,3) = atinfo(3,chnk:hi)
      xt0(1:chnksz) = atinfo(4,chnk:hi)
      xt1(1:chnksz) = 1.0
      do rs2=atmres(chnk),atmres(hi)
        ki = max(chnk,at(rs2)%bb(1))-chnk+1
        kj = min(hi,at(rs2)%bb(1)+at(rs2)%na-1)-chnk+1
        if (disulf(rs1).eq.rs2) xt1(ki:kj) = 0.0
        if (molofrs(rs2).ne.imol) then
          if (bnd_type.eq.1) then
!           cubic box
            if (bnd_shape.eq.1) then
              do j=1,3
                dvec(ki:kj,j) = dvec(ki:kj,j) + lopc(j)*anint((atinfo(j,refat(rs2))-atinfo(j,refat(rs1)))*lopb(j))
              end do
            else if (bnd_shape.eq.3) then
              dvec(ki:kj,3) = dvec(ki:kj,3) + lopc(3)*anint((atinfo(3,refat(rs2))-atinfo(3,refat(rs1)))*lopb(3))
            end if
          end if
        end if
      end do
      if (use_attLJ.EQV..true.) then
        do i=1,at(rs1)%na
!         inverse distance squared and LJ potential
          foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
     &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
     &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
          term6(1:chnksz) = xt1(1:chnksz)*(lj_sig(attyp(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
          foin(1:chnksz,2) = lj_eps(attyp(chnk:hi),k0(i,1))
          foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
          term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
          evec(1) = evec(1) + sum(foin(1:chnksz,3))
          evec(3) = evec(3) - sum(term6(1:chnksz))
          term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3)-6.0*term6(1:chnksz))
!         force 
          foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
          foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
          foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
          ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
          ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
          ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
          ca_f(:,chnk:hi) = ca_f(:,chnk:hi) + transpose(foin(1:chnksz,1:3))
        end do
      else
        do i=1,at(rs1)%na
!         inverse distance squared and LJ potential
          foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
     &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
     &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
          term6(1:chnksz) = xt1(1:chnksz)*(lj_sig(attyp(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
          foin(1:chnksz,2) = lj_eps(attyp(chnk:hi),k0(i,1))
          foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
          evec(1) = evec(1) + sum(foin(1:chnksz,3))
          term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3))
!         force 
          foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
          foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
          foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
          ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
          ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
          ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
          ca_f(:,chnk:hi) = ca_f(:,chnk:hi) + transpose(foin(1:chnksz,1:3))
        end do
      end if
    end do
  else
    return
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the following are the most important routines for standard MD/LD calculations
! these neighbor-list based fxns must be made entirely consistent with the
! residue-based functions (in particular Vforce_rsp) by disabling atom-based cutoff
! truncation there
! otherwise, the set of interactions considered for LJ/WCA will differ slightly
!
! for Coulomb interactions, this is currently not an issue, since
! the residue-based functions do this consistently regardless
!
! the which flag controls the type of neighbor list to use
! different functional forms for electrostatic interactions are picked by is_* type flags
!
subroutine Vforce_PLJ_C(evec,rs1,ca_f,lora,hira3,hira2,which)
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
#ifdef DISABLE_ERFTAB
#ifdef PGI_FORTRAN
  use libm
#endif
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,lo,k,i,shii,jj,ki,kj,imol,chnksz,chnk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)

!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE foin(VC_SIZE,3),xt1(VC_SIZE),xt2(VC_SIZE),xt3(VC_SIZE)
  RTYPE tfc(VC_SIZE,4),ts(VC_SIZE)
  integer bins(VC_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  integer ka(hira2),kb(hira2)
  RTYPE term0(hira2)
  integer sh(hira3-lora+1),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8
!  
! probably seg-faulted already if hira is 0
  hira = hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
  end do
  ii = refat(rs1)
!
  if ((which.eq.1).OR.(which.eq.3)) then ! standard range, recompute shift vectors
    if (which.eq.1) then
      k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    else
      k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira3)
    end if
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else if (which.eq.4) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
! reuse parameters
  lop1 = 4.0*scale_IPP
  lop2 = 4.0*scale_attLJ
  lop3 = electric*scale_POLAR
!
! this loop has data dependency and a conditional -> as few operations as possible in here
  lo = 1
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
  do k=1,hira
    do j=k1(k),k1(k)+sh(k)
      jj = ka(k) + j - k1(k)
      kb(jj) = attyp(j)
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
      term0(jj) = atinfo(4,j)
    end do
  end do
!
! truncated Coulomb
  if ((is_plj.EQV..true.).OR.(is_fegplj.EQV..true.)) then
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       inverse distance squared and LJ potential
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
        term6(1:chnksz) = (lj_sig(kb(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
        term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        evec(3) = evec(3) - sum(term6(1:chnksz))
!       Coulomb potential
        foin(1:chnksz,2) = lop4*xt0(1:chnksz)*sqrt(foin(1:chnksz,1))
        evec(6) = evec(6) + sum(foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*(foin(1:chnksz,2)-6.0*term6(1:chnksz)+12.0*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! Generalized reaction-field (cutoffs must be forced explicitly)
  else if ((is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
    lop5 = mcel_cutoff2
    lop6 = -2.0*par_POLAR(1)
    lop7 = par_POLAR(1)
    lop8 = par_POLAR(2)
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       inverse distance squared and LJ potential
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                         (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                         (dvec(1:chnksz,3)-ixis(3,i))**2
        xt1(1:chnksz) = lop7*foin(1:chnksz,1)-lop8
        foin(1:chnksz,1) = 1.0/foin(1:chnksz,1)
        term6(1:chnksz) = (lj_sig(kb(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
        term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        evec(3) = evec(3) - sum(term6(1:chnksz))
!       RF Coulomb potential ( U_ij = f_0*(1.0/r_ij + f_1*r_ij^2 + f_2) )
        foin(1:chnksz,2) = min(1,int(lop5*foin(1:chnksz,1)))*lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(foin(1:chnksz,2)*(sqrt(foin(1:chnksz,1))+xt1(1:chnksz)))
        term6(1:chnksz) = foin(1:chnksz,2)*(foin(1:chnksz,1)*sqrt(foin(1:chnksz,1))+lop6) + &
 &                         (12.0*foin(1:chnksz,3)-6.0*term6(1:chnksz))*foin(1:chnksz,1)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! transformed Coulomb for Ewald summation
  else if (is_pewlj.EQV..true.) then
    lop5 = ewpm
    lop6 = -ewpm2
    lop7 = ewpite
    lop8 = 1.0/p_terfcoverr(1)
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       inverse distance squared and LJ potential
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                              (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                              (dvec(1:chnksz,3)-ixis(3,i))**2)
        xt1(1:chnksz) = sqrt(foin(1:chnksz,1))
!        foin(1:chnksz,1) = 1.0/foin(1:chnksz,1)
        term6(1:chnksz) = (lj_sig(kb(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        xt2(1:chnksz) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
        term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
        evec(1) = evec(1) + sum(xt2(1:chnksz))
        evec(3) = evec(3) - sum(term6(1:chnksz))
!
#ifdef DISABLE_ERFTAB
        foin(1:chnksz,3) = 1.0/xt1(1:chnksz)
        foin(1:chnksz,2) = xt1(1:chnksz)*erfc(lop5*foin(1:chnksz,3))
        foin(1:chnksz,3) = (lop7*exp(lop6*foin(1:chnksz,3)*foin(1:chnksz,3)) + foin(1:chnksz,2))*foin(1:chnksz,3)
#else
!       Ewald Coulomb potential ( U_ij = f_0*(1.0-erf(r_ij*f_1))/r_ij ) -> via inlined interpolation fxn
        bins(1:chnksz) = min(i_terfcoverr,int(lop8*xt1(1:chnksz)) )
!       this assignment is slow
        tfc(1:chnksz,1) = terfcoverr(bins(1:chnksz),1)
        tfc(1:chnksz,2) = terfcoverr(bins(1:chnksz),2)
        tfc(1:chnksz,3) = terfcoverr(bins(1:chnksz),3)
        tfc(1:chnksz,4) = terfcoverr(bins(1:chnksz),4)
!       compute the cubic spline
        ts(1:chnksz) = lop8*xt1(1:chnksz) - real(bins(1:chnksz))
        foin(1:chnksz,2) = tfc(1:chnksz,1) + ts(1:chnksz)*(tfc(1:chnksz,2) + tfc(1:chnksz,3)*ts(1:chnksz)*ts(1:chnksz) +&
 &                                                         tfc(1:chnksz,4)*ts(1:chnksz))
        foin(1:chnksz,3) = lop8*(tfc(1:chnksz,2) + ts(1:chnksz)*(3.0*tfc(1:chnksz,3)*ts(1:chnksz) + 2.0*tfc(1:chnksz,4)))
#endif
        xt3(1:chnksz) = lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(xt3(1:chnksz)*foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*(12.0*xt2(1:chnksz)-6.0*term6(1:chnksz) + xt3(1:chnksz)*xt1(1:chnksz)*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  end if
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
  end do
!
end
!
!---------------------------------------------------------------------------------------------
!
! this routine employs modified (buffered) forms of LJ and Cb potentials
! it also supports RF electrostatics with a linearly scaled Cb
! otherwise, functionality is equivalent to Vforce_PLJ_C
! the evaluation of the buffered potentials contains the standard potentials as
! limiting cases but is significantly more expensive
!
subroutine Vforce_FEGPLJ_C(evec,rs1,ca_f,lora,hira3,hira2,which)
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
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,lo,k,i,shii,jj,ki,kj,imol,chnksz,chnk,hira
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term0(hira2)
  RTYPE xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE term6(VC_SIZE),foin(VC_SIZE,3),xt1(VC_SIZE),xt2(VC_SIZE),xt3(VC_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na,2)
  integer sh(hira3-lora+1)
  integer ka(hira2),kb(hira2),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE ca_f(3,n)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8,lop9,lop10,lop11
!  
! probably seg-faulted already if hira is 0
  hira = hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
    ixis(4,j) = atq(ii)
  end do
  ii = refat(rs1)
!
  if (which.eq.1) then ! standard range, recompute shift vectors
    k2(1:hira) = rs_nbl(rs1)%gnb(lora:hira3)
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%gnbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%gtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%gtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%gtrsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
! reuse parameters
  lop1 = par_FEG2(1)*4.0
  lop2 = par_FEG2(5)*4.0
  lop3 = electric*par_FEG2(9)
  lop9 = par_FEG2(2)
  lop10 = par_FEG2(6)
  lop11 = par_FEG2(10)
!
! this loop has data dependency and a conditional -> as few operations as possible in here
  lo = 1
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
  do k=1,hira
    do j=k1(k),k1(k)+sh(k)
      jj = ka(k) + j - k1(k)
      kb(jj) = attyp(j)
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
      term0(jj) = atinfo(4,j)
    end do
  end do
!
! ghosted, all with simple truncation
  if (is_fegplj.EQV..true.) then
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       the buffered ghost potentials are quite costly
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                         (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                         (dvec(1:chnksz,3)-ixis(3,i))**2
        xt3(1:chnksz) = 1.0/foin(1:chnksz,1)
        term6(1:chnksz) = (foin(1:chnksz,1)/lj_sig(kb(chnk:hi),k0(i,1)))**3
        xt1(1:chnksz) = 1.0/(term6(1:chnksz)+lop9)
        xt2(1:chnksz) = 1.0/(term6(1:chnksz)+lop10)
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*xt1(1:chnksz)*xt1(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        foin(1:chnksz,3) = foin(1:chnksz,3)*12.0*xt1(1:chnksz)*term6(1:chnksz)*xt3(1:chnksz)
        foin(1:chnksz,2) = lop2*foin(1:chnksz,2)*xt2(1:chnksz)
        evec(3) = evec(3) - sum(foin(1:chnksz,2))
        foin(1:chnksz,2) = foin(1:chnksz,2)*6.0*xt2(1:chnksz)*term6(1:chnksz)*xt3(1:chnksz)
!       ghosted Coulomb potential (this is needlessly expensive if the potential is unbuffered and just scaled)
        xt2(1:chnksz) = sqrt(foin(1:chnksz,1))
        xt1(1:chnksz) = 1.0/(xt2(1:chnksz)+lop11)
        foin(1:chnksz,1) = lop4*xt0(1:chnksz)*xt1(1:chnksz)
        evec(6) = evec(6) + sum(foin(1:chnksz,1))
        term6(1:chnksz) = foin(1:chnksz,1)*xt1(1:chnksz)*sqrt(xt3(1:chnksz))-foin(1:chnksz,2)+foin(1:chnksz,3)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! ghosted, Coulomb with generalized reaction-field (cutoffs must be forced explicitly)
  else if (is_fegprflj.EQV..true.) then
    lop5 = mcel_cutoff2
    lop6 = -2.0*par_POLAR(1)
    lop7 = par_POLAR(1)
    lop8 = par_POLAR(2)
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       the buffered ghost potentials are quite costly
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                         (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                         (dvec(1:chnksz,3)-ixis(3,i))**2
        xt3(1:chnksz) = 1.0/foin(1:chnksz,1)
        term6(1:chnksz) = (foin(1:chnksz,1)/lj_sig(kb(chnk:hi),k0(i,1)))**3
        xt1(1:chnksz) = 1.0/(term6(1:chnksz)+lop9)
        xt2(1:chnksz) = 1.0/(term6(1:chnksz)+lop10)
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*xt1(1:chnksz)*xt1(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        foin(1:chnksz,3) = foin(1:chnksz,3)*12.0*xt1(1:chnksz)*term6(1:chnksz)*xt3(1:chnksz)
        foin(1:chnksz,2) = lop2*foin(1:chnksz,2)*xt2(1:chnksz)
        evec(3) = evec(3) - sum(foin(1:chnksz,2))
        foin(1:chnksz,2) = foin(1:chnksz,2)*6.0*xt2(1:chnksz)*term6(1:chnksz)*xt3(1:chnksz)
!       RF Coulomb potential ( U_ij = f_0*(1.0/r_ij + f_1*r_ij^2 + f_2) ) - only allowed option
!       is to have this simply scaled -> no buffering
        xt2(1:chnksz) = min(1,int(lop5*xt3(1:chnksz)))*lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(xt2(1:chnksz)*(sqrt(xt3(1:chnksz))+lop7*foin(1:chnksz,1)-lop8))
        term6(1:chnksz) = xt2(1:chnksz)*(xt3(1:chnksz)*sqrt(xt3(1:chnksz))+lop6)-foin(1:chnksz,2)+foin(1:chnksz,3)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  end if
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
  end do
!
end
!
!------------------------------------------------------------------------------------------------------
!
! the same as Vforce_PLJ_C with attempted optimizations for 3-site water models
!
! all molecules in the neighbor list have 3 atoms
!
! atom 1 is always the reference atom
!
! atom 1 is the only LJ interaction particle, and the parameters are identical for all sites

! atoms 1-3 carry arbitrary partial charges
!
! this fxn sometimes runs minimally slower than Vforce_PLJ_C, so it is currently not in use
! (standard loops are always used for 3-site water models)
! the theoretical savings come from:
! - compute 1 instead of 9 LJ interactions per molecule pair
! - no need to assign LJ parameters from global arrays for all interactions
!
subroutine Vforce_PLJ_C_T3P(evec,rs1,ca_f,lora,hira3,hira2,which)
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
#ifdef DISABLE_ERFTAB
#ifdef PGI_FORTRAN
  use libm
#endif
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,k,i,jj,imol,chnksz,chnk,kkk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE foin(VC_SIZE,3),xt1(VC_SIZE),xt3(VC_SIZE)
  RTYPE tfc(VC_SIZE,4)
  integer bins(VC_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na)
  RTYPE term0(hira2)
  integer k1(hira3-lora+1),k2(hira3-lora+1),ka(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8,epsik,radik
!  
  hira = hira3-lora+1
!
! initialize
  imol = molofrs(rs1)
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
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
    k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira3)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1))
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1))
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1))
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
    drav(1:hira,1) = drav(1:hira,1) + srav(1:hira,1)
    drav(1:hira,2) = drav(1:hira,2) + srav(1:hira,2)
    drav(1:hira,3) = drav(1:hira,3) + srav(1:hira,3)
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira3)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira3,3)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1)) + srav(1:hira,1)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1)) + srav(1:hira,2)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1)) + srav(1:hira,3)
  else
    return
  end if
!
! we can skip shift vector checks here
  ka(1) = 1
  do k=2,hira
    ka(k) = ka(k-1)+3
  end do
  do k=1,hira
    do j=k1(k),k1(k)+2
      jj = ka(k) + j - k1(k)
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
      term0(jj) = atinfo(4,j)
    end do
  end do
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
    foin(1:chnksz,1) = 1.0/((drav(chnk:hi,1))**2 + (drav(chnk:hi,2))**2 + (drav(chnk:hi,3))**2)
    term6(1:chnksz) = (radik*foin(1:chnksz,1))**3
    foin(1:chnksz,3) = lop1*term6(1:chnksz)*term6(1:chnksz)
    foin(1:chnksz,2) = lop2*term6(1:chnksz)
    evec(1) = evec(1) + sum(foin(1:chnksz,3))
    evec(3) = evec(3) - sum(foin(1:chnksz,2))
    term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3)-6.0*foin(1:chnksz,2))
!   force
    for_k(ka(chnk:hi),1) = drav(chnk:hi,1)*term6(1:chnksz)
    for_k(ka(chnk:hi),2) = drav(chnk:hi,2)*term6(1:chnksz)
    for_k(ka(chnk:hi),3) = drav(chnk:hi,3)*term6(1:chnksz)
    ca_f(1,k0(1)) = ca_f(1,k0(1)) - sum(for_k(ka(chnk:hi),1))
    ca_f(2,k0(1)) = ca_f(2,k0(1)) - sum(for_k(ka(chnk:hi),2))
    ca_f(3,k0(1)) = ca_f(3,k0(1)) - sum(for_k(ka(chnk:hi),3))
  end do
!
! Coulomb
  if ((is_plj.EQV..true.).OR.(is_fegplj.EQV..true.)) then
    do chnk=1,hira2-VC_SIZE,VC_SIZE
      hi = chnk + VC_SIZEM1
      xt0(1:VC_SIZE) = term0(chnk:hi)
      dvec(1:VC_SIZE,1) = dveci(chnk:hi,1)
      dvec(1:VC_SIZE,2) = dveci(chnk:hi,2)
      dvec(1:VC_SIZE,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
        lop4 = lop3*ixis(4,i)
        foin(1:VC_SIZE,1) = 1.0/((dvec(1:VC_SIZE,1)-ixis(1,i))**2 + &
 &                              (dvec(1:VC_SIZE,2)-ixis(2,i))**2 + &
 &                              (dvec(1:VC_SIZE,3)-ixis(3,i))**2)
!       Coulomb potential
        foin(1:VC_SIZE,2) = lop4*xt0(1:VC_SIZE)*sqrt(foin(1:VC_SIZE,1))
        evec(6) = evec(6) + sum(foin(1:VC_SIZE,2))
        term6(1:VC_SIZE) = foin(1:VC_SIZE,1)*foin(1:VC_SIZE,2)
!       force 
        foin(1:VC_SIZE,1) = (dvec(1:VC_SIZE,1)-ixis(1,i))*term6(1:VC_SIZE)
        foin(1:VC_SIZE,2) = (dvec(1:VC_SIZE,2)-ixis(2,i))*term6(1:VC_SIZE)
        foin(1:VC_SIZE,3) = (dvec(1:VC_SIZE,3)-ixis(3,i))*term6(1:VC_SIZE)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:VC_SIZE,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:VC_SIZE,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:VC_SIZE,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:VC_SIZE,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:VC_SIZE,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:VC_SIZE,3)
      end do
    end do
!    chnk = chnk - VC_SIZE 
    hi = min(chnk + VC_SIZEM1,hira2)
    chnksz = hi-chnk+1
    if (chnksz.gt.0) then
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
!       Coulomb potential
        foin(1:chnksz,2) = lop4*xt0(1:chnksz)*sqrt(foin(1:chnksz,1))
        evec(6) = evec(6) + sum(foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*foin(1:chnksz,2)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end if
! Generalized reaction-field (cutoffs must be forced explicitly)
  else if ((is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
    lop5 = mcel_cutoff2
    lop6 = -2.0*par_POLAR(1)
    lop7 = par_POLAR(1)
    lop8 = par_POLAR(2)
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       inverse distance squared
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                         (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                         (dvec(1:chnksz,3)-ixis(3,i))**2
        foin(1:chnksz,2) = lop7*foin(1:chnksz,1)-lop8
        foin(1:chnksz,1) = 1.0/foin(1:chnksz,1)
!       RF Coulomb potential ( U_ij = f_0*(1.0/r_ij + f_1*r_ij^2 + f_2) )
        foin(1:chnksz,3) = min(1,int(lop5*foin(1:chnksz,1)))*lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(foin(1:chnksz,3)*(sqrt(foin(1:chnksz,1))+foin(1:chnksz,2)))
        term6(1:chnksz) = foin(1:chnksz,3)*(foin(1:chnksz,1)*sqrt(foin(1:chnksz,1))+lop6)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! transformed Coulomb for Ewald summation
  else if (is_pewlj.EQV..true.) then
    lop5 = ewpm
    lop6 = -ewpm2
    lop7 = ewpite
    lop8 = 1.0/p_terfcoverr(1)
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
!       inverse distance squared
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                              (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                              (dvec(1:chnksz,3)-ixis(3,i))**2)
        xt1(1:chnksz) = sqrt(foin(1:chnksz,1))
#ifdef DISABLE_ERFTAB
        foin(1:chnksz,3) = 1.0/xt1(1:chnksz)
        foin(1:chnksz,2) = xt1(1:chnksz)*erfc(lop5*foin(1:chnksz,3))
        foin(1:chnksz,3) = (lop7*exp(lop6*foin(1:chnksz,3)*foin(1:chnksz,3)) + foin(1:chnksz,2))*foin(1:chnksz,3)
#else
!       Ewald Coulomb potential ( U_ij = f_0*(1.0-erf(r_ij*f_1))/r_ij ) -> via inlined interpolation fxn
        bins(1:chnksz) = min(i_terfcoverr,int(lop8*xt1(1:chnksz)) )
!       this assignment is slow
        tfc(1:chnksz,1) = terfcoverr(bins(1:chnksz),1)
        tfc(1:chnksz,2) = terfcoverr(bins(1:chnksz),2)
        tfc(1:chnksz,3) = terfcoverr(bins(1:chnksz),3)
        tfc(1:chnksz,4) = terfcoverr(bins(1:chnksz),4)
!       compute the cubic spline
        term6(1:chnksz) = lop8*xt1(1:chnksz) - real(bins(1:chnksz))
        foin(1:chnksz,2) = tfc(1:chnksz,1) + term6(1:chnksz)*(tfc(1:chnksz,2) + &
 &                         tfc(1:chnksz,3)*term6(1:chnksz)*term6(1:chnksz) + tfc(1:chnksz,4)*term6(1:chnksz))
        foin(1:chnksz,3) = lop8*(tfc(1:chnksz,2) + term6(1:chnksz)*(3.0*tfc(1:chnksz,3)*term6(1:chnksz) + 2.0*tfc(1:chnksz,4)))
#endif
        xt3(1:chnksz) = lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(xt3(1:chnksz)*foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*(xt3(1:chnksz)*xt1(1:chnksz)*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  end if
!
  do k=1,hira
    ca_f(:,k1(k):(k1(k)+2)) = ca_f(:,k1(k):(k1(k)+2)) + transpose(for_k(ka(k):ka(k)+2,1:3))
  end do
!
end
!
!-------------------------------------------------------------------------------------------
!
! the same as Vforce_PLJ_C for interactions between standard 4-site water molecules (model: TIP4P)
!
! atom 1 must be reference atom and site of LJ interaction, and all molecules must use same parameters
!
! atoms 2-4 can carry arbitrary charges but must not have LJ interactions
! 
subroutine Vforce_PLJ_C_T4P(evec,rs1,ca_f,lora,hira4,hira2,which)
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
#ifdef DISABLE_ERFTAB
#ifdef PGI_FORTRAN
  use libm
#endif
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira4,hira2,which
!
  integer hira3,ii,j,hi,lo,k,i,jj,imol,chnksz,chnk,kkk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE foin(VC_SIZE,3),xt1(VC_SIZE),xt3(VC_SIZE)
  RTYPE tfc(VC_SIZE,4)
  integer bins(VC_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na)
  RTYPE term0(hira2)
  integer k1(hira4-lora+1),k2(hira4-lora+1),ka(hira4-lora+1)
  RTYPE drav(hira4-lora+1,3),srav(hira4-lora+1,3)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8,epsik,radik
!
  hira = hira4-lora+1
!
! initialize
  imol = molofrs(rs1)
  if (which.eq.1) then
    hira3 = 3*hira ! rs_nbl(rs1)%nwnbs
  else if (which.eq.2) then
    hira3 = 3*hira ! rs_nbl(rs1)%nwnbtrs
  else
    return
  end if
  for_k(1:hira3,1) = 0.0
  for_k(1:hira3,2) = 0.0
  for_k(1:hira3,3) = 0.0
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
    k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira4)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1))
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1))
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1))
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
    drav(1:hira,1) = drav(1:hira,1) + srav(1:hira,1)
    drav(1:hira,2) = drav(1:hira,2) + srav(1:hira,2)
    drav(1:hira,3) = drav(1:hira,3) + srav(1:hira,3)
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira4)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira4,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira4,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira4,3)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1)) + srav(1:hira,1)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1)) + srav(1:hira,2)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1)) + srav(1:hira,3)
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
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
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
    foin(1:chnksz,1) = 1.0/((drav(chnk:hi,1))**2 + (drav(chnk:hi,2))**2 + (drav(chnk:hi,3))**2)
    term6(1:chnksz) = (radik*foin(1:chnksz,1))**3
    foin(1:chnksz,3) = lop1*term6(1:chnksz)*term6(1:chnksz)
    foin(1:chnksz,2) = lop2*term6(1:chnksz)
    evec(1) = evec(1) + sum(foin(1:chnksz,3))
    evec(3) = evec(3) - sum(foin(1:chnksz,2))
    term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3)-6.0*foin(1:chnksz,2))
!   force
    for_k(lo:ii,1) = drav(chnk:hi,1)*term6(1:chnksz)
    for_k(lo:ii,2) = drav(chnk:hi,2)*term6(1:chnksz)
    for_k(lo:ii,3) = drav(chnk:hi,3)*term6(1:chnksz)
    ca_f(1,k0(1)) = ca_f(1,k0(1)) - sum(for_k(lo:ii,1))
    ca_f(2,k0(1)) = ca_f(2,k0(1)) - sum(for_k(lo:ii,2))
    ca_f(3,k0(1)) = ca_f(3,k0(1)) - sum(for_k(lo:ii,3))
  end do
!
! Coulomb
  if ((is_plj.EQV..true.).OR.(is_fegplj.EQV..true.)) then
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=2,at(rs1)%na
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                              (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                              (dvec(1:chnksz,3)-ixis(3,i))**2)
!       Coulomb potential
        foin(1:chnksz,2) = lop4*xt0(1:chnksz)*sqrt(foin(1:chnksz,1))
        evec(6) = evec(6) + sum(foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*foin(1:chnksz,2)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! Generalized reaction-field (cutoffs must be forced explicitly)
  else if ((is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
    lop5 = mcel_cutoff2
    lop6 = -2.0*par_POLAR(1)
    lop7 = par_POLAR(1)
    lop8 = par_POLAR(2)
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=2,at(rs1)%na
!       inverse distance squared
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                         (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                         (dvec(1:chnksz,3)-ixis(3,i))**2
        foin(1:chnksz,2) = lop7*foin(1:chnksz,1)-lop8
        foin(1:chnksz,1) = 1.0/foin(1:chnksz,1)
!       RF Coulomb potential ( U_ij = f_0*(1.0/r_ij + f_1*r_ij^2 + f_2) )
        foin(1:chnksz,3) = min(1,int(lop5*foin(1:chnksz,1)))*lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(foin(1:chnksz,3)*(sqrt(foin(1:chnksz,1))+foin(1:chnksz,2)))
        term6(1:chnksz) = foin(1:chnksz,3)*(foin(1:chnksz,1)*sqrt(foin(1:chnksz,1))+lop6)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! transformed Coulomb for Ewald summation
  else if (is_pewlj.EQV..true.) then
    lop5 = ewpm
    lop6 = -ewpm2
    lop7 = ewpite
    lop8 = 1.0/p_terfcoverr(1)
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=2,at(rs1)%na
!       inverse distance squared
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                              (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                              (dvec(1:chnksz,3)-ixis(3,i))**2)
        xt1(1:chnksz) = sqrt(foin(1:chnksz,1))
#ifdef DISABLE_ERFTAB
        foin(1:chnksz,3) = 1.0/xt1(1:chnksz)
        foin(1:chnksz,2) = xt1(1:chnksz)*erfc(lop5*foin(1:chnksz,3))
        foin(1:chnksz,3) = (lop7*exp(lop6*foin(1:chnksz,3)*foin(1:chnksz,3)) + foin(1:chnksz,2))*foin(1:chnksz,3)
#else
!       Ewald Coulomb potential ( U_ij = f_0*(1.0-erf(r_ij*f_1))/r_ij ) -> via inlined interpolation fxn
        bins(1:chnksz) = min(i_terfcoverr,int(lop8*xt1(1:chnksz)) )
!       this assignment is slow
        tfc(1:chnksz,1) = terfcoverr(bins(1:chnksz),1)
        tfc(1:chnksz,2) = terfcoverr(bins(1:chnksz),2)
        tfc(1:chnksz,3) = terfcoverr(bins(1:chnksz),3)
        tfc(1:chnksz,4) = terfcoverr(bins(1:chnksz),4)
!       compute the cubic spline
        term6(1:chnksz) = lop8*xt1(1:chnksz) - real(bins(1:chnksz))
        foin(1:chnksz,2) = tfc(1:chnksz,1) + term6(1:chnksz)*(tfc(1:chnksz,2) + &
 &                         tfc(1:chnksz,3)*term6(1:chnksz)*term6(1:chnksz) + tfc(1:chnksz,4)*term6(1:chnksz))
        foin(1:chnksz,3) = lop8*(tfc(1:chnksz,2) + term6(1:chnksz)*(3.0*tfc(1:chnksz,3)*term6(1:chnksz) + 2.0*tfc(1:chnksz,4)))
#endif
        xt3(1:chnksz) = lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(xt3(1:chnksz)*foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*(xt3(1:chnksz)*xt1(1:chnksz)*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  end if
!
  k1(1:hira) = k1(1:hira) + 1
  do k=1,hira
    ca_f(:,k1(k)-1) = ca_f(:,k1(k)-1) + for_k(k+hira3,1:3)
    ca_f(:,k1(k):(k1(k)+2)) = ca_f(:,k1(k):(k1(k)+2)) + transpose(for_k(ka(k):(ka(k)+2),1:3))
  end do
!
end
!
!-------------------------------------------------------------------------------------------
!
! the same as Vforce_PLJ_C for interactions between standard five site water molecules (model: TIP5P)
!
! atom 1 must be reference atom and site of LJ interaction, and all molecules must use same parameters
!
! atoms 2-5 can carry arbitrary charges but must not have LJ interactions
! 
subroutine Vforce_PLJ_C_T5P(evec,rs1,ca_f,lora,hira4,hira2,which)
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
#ifdef DISABLE_ERFTAB
#ifdef PGI_FORTRAN
  use libm
#endif
#endif
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira4,hira2,which
!
  integer hira3,ii,j,hi,lo,k,i,jj,ki,kj,imol,chnksz,chnk,kkk,hira
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
!
  RTYPE term6(VC_SIZE),xt0(VC_SIZE),dvec(VC_SIZE,3)
  RTYPE foin(VC_SIZE,3),xt1(VC_SIZE),xt3(VC_SIZE)
  RTYPE tfc(VC_SIZE,4)
  integer bins(VC_SIZE)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(4,at(rs1)%na)
  integer k0(at(rs1)%na)
  RTYPE term0(hira2)
  integer k1(hira4-lora+1),k2(hira4-lora+1),ka(hira4-lora+1)
  RTYPE drav(hira4-lora+1,3),srav(hira4-lora+1,3)
  RTYPE lop1,lop2,lop3,lop4,lop5,lop6,lop7,lop8,epsik,radik
!  
  hira = hira4-lora+1
!
! initialize
  imol = molofrs(rs1)
  if (which.eq.1) then
    hira3 = 4*hira ! rs_nbl(rs1)%nwnbs
  else if (which.eq.2) then
    hira3 = 4*hira ! rs_nbl(rs1)%nwnbtrs
  else
    return
  end if
  for_k(1:hira3,1) = 0.0
  for_k(1:hira3,2) = 0.0
  for_k(1:hira3,3) = 0.0
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
    k2(1:hira) = rs_nbl(rs1)%wnb(lora:hira4)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1))
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1))
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1))
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
    drav(1:hira,1) = drav(1:hira,1) + srav(1:hira,1)
    drav(1:hira,2) = drav(1:hira,2) + srav(1:hira,2)
    drav(1:hira,3) = drav(1:hira,3) + srav(1:hira,3)
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%wnbtr(lora:hira4)
    k1(1:hira) = rsinfo(k2(1:hira),1)
    srav(1:hira,1) = rs_nbl(rs1)%wtrsvec(lora:hira4,1)
    srav(1:hira,2) = rs_nbl(rs1)%wtrsvec(lora:hira4,2)
    srav(1:hira,3) = rs_nbl(rs1)%wtrsvec(lora:hira4,3)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,k0(1)) + srav(1:hira,1)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,k0(1)) + srav(1:hira,2)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,k0(1)) + srav(1:hira,3)
  end if
!
! we can skip shift vector checks here
  ka(1) = 1
  do k=2,hira
    ka(k) = ka(k-1)+4
  end do
  do k=1,hira
    do j=k1(k)+1,k1(k)+4
      jj = ka(k) + j - k1(k) - 1
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
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
    foin(1:chnksz,1) = 1.0/((drav(chnk:hi,1))**2 + (drav(chnk:hi,2))**2 + (drav(chnk:hi,3))**2)
    term6(1:chnksz) = (radik*foin(1:chnksz,1))**3
    foin(1:chnksz,3) = lop1*term6(1:chnksz)*term6(1:chnksz)
    foin(1:chnksz,2) = lop2*term6(1:chnksz)
    evec(1) = evec(1) + sum(foin(1:chnksz,3))
    evec(3) = evec(3) - sum(foin(1:chnksz,2))
    term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3)-6.0*foin(1:chnksz,2))
!   force
    for_k(lo:ii,1) = drav(chnk:hi,1)*term6(1:chnksz)
    for_k(lo:ii,2) = drav(chnk:hi,2)*term6(1:chnksz)
    for_k(lo:ii,3) = drav(chnk:hi,3)*term6(1:chnksz)
    ca_f(1,k0(1)) = ca_f(1,k0(1)) - sum(for_k(lo:ii,1))
    ca_f(2,k0(1)) = ca_f(2,k0(1)) - sum(for_k(lo:ii,2))
    ca_f(3,k0(1)) = ca_f(3,k0(1)) - sum(for_k(lo:ii,3))
  end do
!
! Coulomb
  if ((is_plj.EQV..true.).OR.(is_fegplj.EQV..true.)) then
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=2,at(rs1)%na
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                              (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                              (dvec(1:chnksz,3)-ixis(3,i))**2)
!       Coulomb potential
        foin(1:chnksz,2) = lop4*xt0(1:chnksz)*sqrt(foin(1:chnksz,1))
        evec(6) = evec(6) + sum(foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*foin(1:chnksz,2)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! Generalized reaction-field (cutoffs must be forced explicitly)
  else if ((is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
    lop5 = mcel_cutoff2
    lop6 = -2.0*par_POLAR(1)
    lop7 = par_POLAR(1)
    lop8 = par_POLAR(2)
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=2,at(rs1)%na
!       inverse distance squared
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                         (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                         (dvec(1:chnksz,3)-ixis(3,i))**2
        foin(1:chnksz,2) = lop7*foin(1:chnksz,1)-lop8
        foin(1:chnksz,1) = 1.0/foin(1:chnksz,1)
!       RF Coulomb potential ( U_ij = f_0*(1.0/r_ij + f_1*r_ij^2 + f_2) )
        foin(1:chnksz,3) = min(1,int(lop5*foin(1:chnksz,1)))*lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(foin(1:chnksz,3)*(sqrt(foin(1:chnksz,1))+foin(1:chnksz,2)))
        term6(1:chnksz) = foin(1:chnksz,3)*(foin(1:chnksz,1)*sqrt(foin(1:chnksz,1))+lop6)
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
! transformed Coulomb for Ewald summation
  else if (is_pewlj.EQV..true.) then
    lop5 = ewpm
    lop6 = -ewpm2
    lop7 = ewpite
    lop8 = 1.0/p_terfcoverr(1)
    do chnk=1,hira3,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira3)
      chnksz = hi-chnk+1
      xt0(1:chnksz) = term0(chnk:hi)
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=2,at(rs1)%na
!       inverse distance squared
        lop4 = lop3*ixis(4,i)
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
 &                              (dvec(1:chnksz,2)-ixis(2,i))**2 + &
 &                              (dvec(1:chnksz,3)-ixis(3,i))**2)
        xt1(1:chnksz) = sqrt(foin(1:chnksz,1))
#ifdef DISABLE_ERFTAB
        foin(1:chnksz,3) = 1.0/xt1(1:chnksz)
        foin(1:chnksz,2) = xt1(1:chnksz)*erfc(lop5*foin(1:chnksz,3))
        foin(1:chnksz,3) = (lop7*exp(lop6*foin(1:chnksz,3)*foin(1:chnksz,3)) + foin(1:chnksz,2))*foin(1:chnksz,3)
#else
!       Ewald Coulomb potential ( U_ij = f_0*(1.0-erf(r_ij*f_1))/r_ij ) -> via inlined interpolation fxn
        bins(1:chnksz) = min(i_terfcoverr,int(lop8*xt1(1:chnksz)) )
!       this assignment is slow
        tfc(1:chnksz,1) = terfcoverr(bins(1:chnksz),1)
        tfc(1:chnksz,2) = terfcoverr(bins(1:chnksz),2)
        tfc(1:chnksz,3) = terfcoverr(bins(1:chnksz),3)
        tfc(1:chnksz,4) = terfcoverr(bins(1:chnksz),4)
!       compute the cubic spline
        term6(1:chnksz) = lop8*xt1(1:chnksz) - real(bins(1:chnksz))
        foin(1:chnksz,2) = tfc(1:chnksz,1) + term6(1:chnksz)*(tfc(1:chnksz,2) + &
 &                         tfc(1:chnksz,3)*term6(1:chnksz)*term6(1:chnksz) + tfc(1:chnksz,4)*term6(1:chnksz))
        foin(1:chnksz,3) = lop8*(tfc(1:chnksz,2) + term6(1:chnksz)*(3.0*tfc(1:chnksz,3)*term6(1:chnksz) + 2.0*tfc(1:chnksz,4)))
#endif
        xt3(1:chnksz) = lop4*xt0(1:chnksz)
        evec(6) = evec(6) + sum(xt3(1:chnksz)*foin(1:chnksz,2))
        term6(1:chnksz) = foin(1:chnksz,1)*(xt3(1:chnksz)*xt1(1:chnksz)*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i)) = ca_f(1,k0(i)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i)) = ca_f(2,k0(i)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i)) = ca_f(3,k0(i)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  end if
!
  k1(1:hira) = k1(1:hira) + 1
  do k=1,hira
    ki = k1(k)+3
    kj = ka(k)+3
    ca_f(:,k1(k)-1) = ca_f(:,k1(k)-1) + for_k(k+hira3,1:3)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! just the LJ potential
!
subroutine Vforce_LJ_C(evec,rs1,ca_f,lora,hira3,hira2,which)
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
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,j,hi,lo,k,i,shii,jj,ki,kj,imol,chnksz,chnk,hira
  RTYPE evec(MAXENERGYTERMS)
  RTYPE term6(VC_SIZE),dvec(VC_SIZE,3)!,dvex(VC_SIZE),dvey(VC_SIZE),dvez(VC_SIZE),t0(VC_SIZE)
  RTYPE foin(VC_SIZE,3)
  RTYPE dveci(hira2,3)
  RTYPE for_k(hira2,3)
  RTYPE ixis(3,at(rs1)%na)
  integer k0(at(rs1)%na,3)
  integer sh(hira3-lora+1)
  integer ka(hira2),kb(hira2),k1(hira3-lora+1),k2(hira3-lora+1)
  RTYPE drav(hira3-lora+1,3),srav(hira3-lora+1,3)
  RTYPE ca_f(3,n)
  RTYPE lop1,lop2,lop3
!  
! probably seg-faulted already if hira is 0
  hira = hira3-lora+1
  if ((hira.le.0).OR.(hira2.lt.hira)) return
!
! initialize
  for_k(:,1) = 0.0
  for_k(:,2) = 0.0
  for_k(:,3) = 0.0
!
  imol = molofrs(rs1)
  shii = at(rs1)%na - 1
  k0(:,3) = 1
  do ii=rsinfo(rs1,1),rsinfo(rs1,1)+shii
    j = ii-rsinfo(rs1,1)+1
    k0(j,1) = attyp(ii)
    if (maxval(lj_eps(:,attyp(ii))).le.0.0) k0(j,3) = 0
    k0(j,2) = ii
    ixis(1,j) = x(ii)
    ixis(2,j) = y(ii)
    ixis(3,j) = z(ii)
  end do
  ii = refat(rs1)
!
  if (which.eq.1) then ! standard range, recompute shift vectors
    k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    k1(1:hira) = refat(k2(1:hira))
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else
    return
  end if
  sh(1:hira) = rsinfo(k2(1:hira),2)
!
! reuse parameters
  lop1 = 4.0*scale_IPP
  lop2 = 4.0*scale_attLJ
!
! this loop has data dependency and a conditional -> as few operations as possible in here
  lo = 1
  if (molofrs(k2(1)).eq.imol) srav(1,:) = 0.0
  ka(1) = 1
  k1(1) = rsinfo(k2(1),1)
  do k=2,hira
    k1(k) =  rsinfo(k2(k),1)
    if (molofrs(k2(k)).eq.imol) srav(k,:) = 0.0
    ka(k) = ka(k-1)+sh(k-1)+1
  end do
  do k=1,hira
    do j=k1(k),k1(k)+sh(k)
      jj = ka(k) + j - k1(k)
      kb(jj) = attyp(j)
      dveci(jj,1) = srav(k,1) + atinfo(1,j)
      dveci(jj,2) = srav(k,2) + atinfo(2,j)
      dveci(jj,3) = srav(k,3) + atinfo(3,j)
    end do
  end do
!
  if (is_lj.EQV..true.) then
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
        if (k0(i,3).eq.0) cycle
!       use vector forms as much as possible
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
        term6(1:chnksz) = (lj_sig(kb(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
        term6(1:chnksz) = lop2*foin(1:chnksz,2)*term6(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        evec(3) = evec(3) - sum(term6(1:chnksz))
        term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3)-6.0*term6(1:chnksz))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  else ! must be is_ev 
    do chnk=1,hira2,VC_SIZE
      hi = min(chnk + VC_SIZEM1,hira2)
      chnksz = hi-chnk+1
      dvec(1:chnksz,1) = dveci(chnk:hi,1)
      dvec(1:chnksz,2) = dveci(chnk:hi,2)
      dvec(1:chnksz,3) = dveci(chnk:hi,3)
      do i=1,at(rs1)%na
        if (k0(i,3).eq.0) cycle
!       use vector forms as much as possible
        foin(1:chnksz,2) = lj_eps(kb(chnk:hi),k0(i,1))
        foin(1:chnksz,1) = 1.0/((dvec(1:chnksz,1)-ixis(1,i))**2 + &
   &                            (dvec(1:chnksz,2)-ixis(2,i))**2 + &
   &                            (dvec(1:chnksz,3)-ixis(3,i))**2)
        term6(1:chnksz) = (lj_sig(kb(chnk:hi),k0(i,1))*foin(1:chnksz,1))**3
        foin(1:chnksz,3) = lop1*foin(1:chnksz,2)*term6(1:chnksz)*term6(1:chnksz)
        evec(1) = evec(1) + sum(foin(1:chnksz,3))
        term6(1:chnksz) = foin(1:chnksz,1)*(12.0*foin(1:chnksz,3))
!       force 
        foin(1:chnksz,1) = (dvec(1:chnksz,1)-ixis(1,i))*term6(1:chnksz)
        foin(1:chnksz,2) = (dvec(1:chnksz,2)-ixis(2,i))*term6(1:chnksz)
        foin(1:chnksz,3) = (dvec(1:chnksz,3)-ixis(3,i))*term6(1:chnksz)
        ca_f(1,k0(i,2)) = ca_f(1,k0(i,2)) - sum(foin(1:chnksz,1))
        ca_f(2,k0(i,2)) = ca_f(2,k0(i,2)) - sum(foin(1:chnksz,2))
        ca_f(3,k0(i,2)) = ca_f(3,k0(i,2)) - sum(foin(1:chnksz,3))
        for_k(chnk:hi,1) = for_k(chnk:hi,1) + foin(1:chnksz,1)
        for_k(chnk:hi,2) = for_k(chnk:hi,2) + foin(1:chnksz,2)
        for_k(chnk:hi,3) = for_k(chnk:hi,3) + foin(1:chnksz,3)
      end do
    end do
  end if
!
  do k=1,hira
    ki = k1(k)+sh(k)
    kj = ka(k)+sh(k)
    ca_f(:,k1(k):ki) = ca_f(:,k1(k):ki) + transpose(for_k(ka(k):kj,1:3))
  end do
!
end
!
!-------------------------------------------------------------------------------
!
! a specialized routine suitable for sparse tabulated potentials
! note that this is not crosslink-adjusted  and does not require 
! calls to the correction fxn in the wrapper construct (force_wrap)
!
subroutine Vforce_TAB(evec,rs1,ca_f)
!
  use energies
  use polypep
  use forces
  use atoms
  use molecule
  use sequen
  use params
  use tabpot
  use units
  use math
  use cutoffs
!
  implicit none
!
  integer, INTENT(IN):: rs1
!
  integer rs2,ii,kk,hi,lo,alcsz,k,i,jj,k2
  RTYPE sv(3),pfac,pfac2,hlp1,hlp2,hlp3,hlp4,dhlp1,dhlp2,dhlp3,evec(MAXENERGYTERMS)
  RTYPE dvec(rs_nbl(rs1)%ntabias,3)
  RTYPE d2(rs_nbl(rs1)%ntabias)
  RTYPE term0(rs_nbl(rs1)%ntabias),term1(rs_nbl(rs1)%ntabias),term2(rs_nbl(rs1)%ntabias),term3(rs_nbl(rs1)%ntabias)
  RTYPE id1(rs_nbl(rs1)%ntabias)
  RTYPE d1(rs_nbl(rs1)%ntabias),ca_f(3,n)
  integer tbin(rs_nbl(rs1)%ntabias)
!  
  lo = 1 
  do jj=1,rs_nbl(rs1)%ntabnbs
    rs2 = rs_nbl(rs1)%tabnb(jj)
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    call dis_bound_rs2(rs1,rs2,sv)
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
      kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
      dvec(lo+k,1) = x(kk) - x(ii)
      dvec(lo+k,2) = y(kk) - y(ii)
      dvec(lo+k,3) = z(kk) - z(ii)
      k = k + 1
    end do
    dvec(lo:hi,1) =  dvec(lo:hi,1) + sv(1)
    dvec(lo:hi,2) =  dvec(lo:hi,2) + sv(2)
    dvec(lo:hi,3) =  dvec(lo:hi,3) + sv(3)
    lo = hi + 1
  end do
!
  alcsz = hi
  d2(1:alcsz) = dvec(1:alcsz,1)**2 + dvec(1:alcsz,2)**2 + dvec(1:alcsz,3)**2
  d1(1:alcsz) = sqrt(d2(1:alcsz))
  pfac = 1.0/tbp%res
  pfac2 = pfac*scale_TABUL
  term0(1:alcsz) = pfac*(d1(1:alcsz)-tbp%dis(1))
  tbin(1:alcsz) = min(tbp%bins,max(1,ceiling(term0(1:alcsz))))
  term1(1:alcsz) = pfac*(d1(1:alcsz) - tbp%dis(tbin(1:alcsz))) ! linear
  term2(1:alcsz) = term1(1:alcsz)*term1(1:alcsz) ! quadratic
  term3(1:alcsz) = term2(1:alcsz)*term1(1:alcsz) ! cubic
  id1(1:alcsz) = 1.0/d1(1:alcsz)
!  tbin(1:alcsz) = floor((d1(1:alcsz)-tbp%dis(1))/tbp%res) + 1
!  pfac = scale_TABUL*(1.0/tbp%res)
!
  lo = 1
  do jj=1,rs_nbl(rs1)%ntabnbs
    rs2 = rs_nbl(rs1)%tabnb(jj) ! always > rs1
    if ((rs2-rs1).eq.1) then
!     essentially cycle since tabnb still contains sequence neighbors
      hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
      lo = hi + 1
      cycle
    end if
    hi = lo + (tbp%rsmat(rs2,rs1)-tbp%rsmat(rs1,rs2))
    k = 0
    do i=tbp%rsmat(rs1,rs2),tbp%rsmat(rs2,rs1)
      if (tbin(lo+k).ge.tbp%bins) then
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),tbp%bins)
      else if (tbin(lo+k).ge.1) then
        k2 = lo+k
        ii = atmol(molofrs(rs1),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1) 
        kk = atmol(molofrs(rs2),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
        hlp1 = 2.0*term3(k2) - 3.0*term2(k2)
        hlp2 = term3(k2) - term2(k2)
        hlp3 = hlp2 - term2(k2) + term1(k2)
        dhlp1 = 6.0*term2(k2) - 6.0*term1(k2)
        dhlp2 = 3.0*term2(k2) - 2.0*term1(k2)
        dhlp3 = dhlp2 - 2.0*term1(k2) + 1.0
        hlp4 = (hlp1+1.0)*tbp%pot(tbp%lst(i,3),tbin(k2)) - hlp1*tbp%pot(tbp%lst(i,3),tbin(k2)+1) + &
 &               tbp%res*(hlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + hlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1))
!        term1 = (tbp%dis(tbin(k2)+1) - d1(k2))*tbp%pot(tbp%lst(i,3),tbin(k2))
!        term2 = (d1(k2) - tbp%dis(tbin(k2)))*tbp%pot(tbp%lst(i,3),tbin(k2)+1)
!        term3 = pfac*(term1+term2)
        evec(9) = evec(9) + scale_TABUL*hlp4
        hlp4 = pfac2*id1(k2)*(dhlp1*(tbp%pot(tbp%lst(i,3),tbin(k2)) - tbp%pot(tbp%lst(i,3),tbin(k2)+1)) + &
 &               tbp%res*(dhlp3*tbp%tang(tbp%lst(i,3),tbin(k2)) + dhlp2*tbp%tang(tbp%lst(i,3),tbin(k2)+1)))
!        term4 = id1(k2)*pfac*(tbp%pot(tbp%lst(i,3),tbin(k2)+1)-tbp%pot(tbp%lst(i,3),tbin(k2)))
        ca_f(:,ii) = ca_f(:,ii) + hlp4*dvec(k2,1:3)
        ca_f(:,kk) = ca_f(:,kk) - hlp4*dvec(k2,1:3)
      else
        evec(9) = evec(9) + scale_TABUL*tbp%pot(tbp%lst(i,3),1)
      end if
      k = k + 1
    end do
    lo = hi + 1
  end do
!
  end
!
!-------------------------------------------------------------------------------
!
! the same thing for standard NB-list
!
subroutine Vforce_TAB_C(evec,rs1,ca_f,lora,hira3,hira2,which)
!
  use energies
  use polypep
  use forces
  use atoms
  use molecule
  use sequen
  use params
  use tabpot
  use units
  use math
  use cutoffs
  use system
!
  implicit none
!
  integer, INTENT(IN):: rs1,lora,hira3,hira2,which
!
  integer ii,hi,lo,k,i,j,jj,chnk,chnksz,rslo,rshi,hira
!
  RTYPE srav(hira3-lora+1,3),drav(hira3-lora+1,3)
  integer k2(hira3-lora+1),k1(hira3-lora+1)
!
  RTYPE dvec(hira2,3)
  integer kx(hira2,2),ky(hira2)
!
  RTYPE id1(VC_SIZE),d1(VC_SIZE),term0(VC_SIZE),term1(VC_SIZE),term2(VC_SIZE),term3(VC_SIZE),term4(VC_SIZE)
  RTYPE term5(VC_SIZE),term6(VC_SIZE),term7(VC_SIZE),term8(VC_SIZE)
  integer tbin(VC_SIZE)
!
  RTYPE ca_f(3,n)
  RTYPE lop1,lop2,lop3,hlp4,evec(MAXENERGYTERMS)
!
  hira =  hira3-lora+1
!
  if (which.eq.1) then ! standard range, recompute shift vectors
    k2(1:hira) = rs_nbl(rs1)%nb(lora:hira3)
    k1(1:hira) = refat(k2(1:hira))
    ii = refat(rs1)
    drav(1:hira,1) = atinfo(1,k1(1:hira)) - atinfo(1,ii)
    drav(1:hira,2) = atinfo(2,k1(1:hira)) - atinfo(2,ii)
    drav(1:hira,3) = atinfo(3,k1(1:hira)) - atinfo(3,ii)
!   PBC
    if (bnd_type.eq.1) then
!     cubic box
      if (bnd_shape.eq.1) then
        do j=1,3
          lop1 = 1.0/bnd_params(j)
          lop3 = -bnd_params(j)
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      else if (bnd_shape.eq.3) then
        srav(1:hira,1:2) = 0.0
        lop1 = 1.0/bnd_params(6)
        lop3 = -bnd_params(6)
        do j=3,3
          srav(1:hira,j) = lop3*anint(drav(1:hira,j)*lop1)
        end do
      end if
    else 
      srav(:,:) = 0.0
    end if
  else if (which.eq.2) then ! twin-range, use precomputed shift vectors
    k2(1:hira) = rs_nbl(rs1)%nbtr(lora:hira3)
    srav(1:hira,1) = rs_nbl(rs1)%trsvec(lora:hira3,1)
    srav(1:hira,2) = rs_nbl(rs1)%trsvec(lora:hira3,2)
    srav(1:hira,3) = rs_nbl(rs1)%trsvec(lora:hira3,3)
  else
    return
  end if
!     
  lo = 1 
  do jj=1,hira
    rshi = max(rs1,k2(jj))
    rslo = min(rs1,k2(jj)) 
    hi = lo + (tbp%rsmat(rshi,rslo)-tbp%rsmat(rslo,rshi))
    k = 0
    do i=tbp%rsmat(rslo,rshi),tbp%rsmat(rshi,rslo)
      kx(lo+k,1) = atmol(molofrs(rslo),1) + tbp%lst(i,1) - atmol(molofrs(atmres(tbp%lst(i,1))),1)
      kx(lo+k,2) = atmol(molofrs(rshi),1) + tbp%lst(i,2) - atmol(molofrs(atmres(tbp%lst(i,2))),1)
      dvec(lo+k,1) = x(kx(lo+k,2)) - x(kx(lo+k,1))
      dvec(lo+k,2) = y(kx(lo+k,2)) - y(kx(lo+k,1))
      dvec(lo+k,3) = z(kx(lo+k,2)) - z(kx(lo+k,1))
      k = k + 1
    end do
    ky(lo:hi) = tbp%lst(tbp%rsmat(rslo,rshi):tbp%rsmat(rshi,rslo),3)
    if (rslo.eq.rs1) then
      dvec(lo:hi,1) =  dvec(lo:hi,1) + srav(jj,1)
      dvec(lo:hi,2) =  dvec(lo:hi,2) + srav(jj,2)
      dvec(lo:hi,3) =  dvec(lo:hi,3) + srav(jj,3)
    else
      dvec(lo:hi,1) =  dvec(lo:hi,1) - srav(jj,1)
      dvec(lo:hi,2) =  dvec(lo:hi,2) - srav(jj,2)
      dvec(lo:hi,3) =  dvec(lo:hi,3) - srav(jj,3)
    end if
    lo = hi + 1
  end do
!
  lop1 = 1.0/tbp%res
  lop2 = scale_TABUL
  lop3 = lop1*scale_TABUL
!
  do chnk=1,hira2,VC_SIZE
    hi = min(chnk + VC_SIZEM1,hira2)
    chnksz = hi-chnk+1
!
    d1(1:chnksz) = sqrt(dvec(chnk:hi,1)**2 + dvec(chnk:hi,2)**2 + dvec(chnk:hi,3)**2)
    id1(1:chnksz) = 1.0/d1(1:chnksz)
    term0(1:chnksz) = lop1*(d1(1:chnksz)-tbp%dis(1))
    tbin(1:chnksz) = min(tbp%bins,max(1,ceiling(term0(1:chnksz))))
    term1(1:chnksz) = lop1*(d1(1:chnksz) - tbp%dis(tbin(1:chnksz))) ! linear
    term2(1:chnksz) = term1(1:chnksz)*term1(1:chnksz) ! quadratic
    term3(1:chnksz) = term2(1:chnksz)*term1(1:chnksz) ! cubic
!
    term0(1:chnksz) = 2.0*term3(1:chnksz) - 3.0*term2(1:chnksz)
    term4(1:chnksz) = term3(1:chnksz) - term2(1:chnksz)
    term5(1:chnksz) = term4(1:chnksz) - term2(1:chnksz) + term1(1:chnksz)
    term6(1:chnksz) = 6.0*term2(1:chnksz) - 6.0*term1(1:chnksz)
    term7(1:chnksz) = 3.0*term2(1:chnksz) - 2.0*term1(1:chnksz)
    term8(1:chnksz) = term7(1:chnksz) - 2.0*term1(1:chnksz) + 1.0
!
    do i=1,chnksz
      k = chnk + i - 1
      if (tbin(i).ge.tbp%bins) then
        evec(9) = evec(9) + lop2*tbp%pot(ky(k),tbp%bins)
      else if (tbin(i).ge.1) then
        hlp4 = (term0(i)+1.0)*tbp%pot(ky(k),tbin(i)) - term0(i)*tbp%pot(ky(k),tbin(i)+1) + &
 &            tbp%res*(term5(i)*tbp%tang(ky(k),tbin(i)) + term4(i)*tbp%tang(ky(k),tbin(i)+1))
        evec(9) = evec(9) + lop2*hlp4
        hlp4 = lop3*id1(i)*(term6(i)*(tbp%pot(ky(k),tbin(i)) - tbp%pot(ky(k),tbin(i)+1)) + &
 &               tbp%res*(term8(i)*tbp%tang(ky(k),tbin(i)) + term7(i)*tbp%tang(ky(k),tbin(i)+1)))
        ca_f(:,kx(k,1)) = ca_f(:,kx(k,1)) + hlp4*dvec(k,1:3)
        ca_f(:,kx(k,2)) = ca_f(:,kx(k,2)) - hlp4*dvec(k,1:3)
      else
        evec(9) = evec(9) + lop2*tbp%pot(ky(k),1)
      end if
    end do
  end do
!
end
!
!--------------------------------------------------------------------------------------------------------
!
