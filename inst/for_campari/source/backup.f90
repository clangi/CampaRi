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
! CONTRIBUTIONS: Rohit Pappu, Adam Steffen, Hoang Tran                     !
!                                                                          !
!--------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------
!
! a set of trivial routines which copy subsets of atomic coordinates
! from or to backup arrays
!
!-----------------------------------------------------------------------
!
!
subroutine makeref_forrotlst(ati,alC,tpi)
!
  use atoms
  use zmatrix
  use sequen
  use molecule
  use atoms
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,ati
  logical, INTENT(IN):: alC
!
  integer ii,ik,k,imol,shf,shf2
#ifdef ENABLE_THREADS
  integer stx(2),tpn
!
  tpn = thrdat%maxn
#endif
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
!
  if (alC.EQV..true.) then
    shf = 1
    shf2 = 1
    imol = molofrs(atmres(ati))
    if (izrot(ati)%rotis(1,1).eq.atmol(imol,1)) shf = 2
    if (izrot(ati)%rotis(izrot(ati)%alsz,2).eq.atmol(imol,2)) shf2 = 0
    do k=shf,izrot(ati)%alsz+shf2
      if (k.eq.1) then 
        ii = atmol(imol,1)
      else
        ii = izrot(ati)%rotis(k-1,2)+1
      end if
      if (k.gt.izrot(ati)%alsz) then 
        ik = atmol(imol,2)
      else
        ik = izrot(ati)%rotis(k,1)-1
      end if
#ifdef ENABLE_THREADS
      if ((ik-ii+1).gt.10*thrdat%maxn) then
        call threads_bounds(ii,ik,tpi,tpn,stx)
        xref(stx(1):stx(2)) = x(stx(1):stx(2))
        yref(stx(1):stx(2)) = y(stx(1):stx(2))
        zref(stx(1):stx(2)) = z(stx(1):stx(2))
      else
!$OMP SINGLE
#endif
      xref(ii:ik) = x(ii:ik)
      yref(ii:ik) = y(ii:ik)
      zref(ii:ik) = z(ii:ik)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end do
  else
    do k=1,izrot(ati)%alsz
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
#ifdef ENABLE_THREADS
      if ((ik-ii+1).gt.10*thrdat%maxn) then
        call threads_bounds(ii,ik,tpi,tpn,stx)
        xref(stx(1):stx(2)) = x(stx(1):stx(2))
        yref(stx(1):stx(2)) = y(stx(1):stx(2))
        zref(stx(1):stx(2)) = z(stx(1):stx(2))
      else
!$OMP SINGLE
#endif
      xref(ii:ik) = x(ii:ik)
      yref(ii:ik) = y(ii:ik)
      zref(ii:ik) = z(ii:ik)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine makeref_forsc(rs)
!
  use atoms
  use polypep
!
  implicit none
!
  integer i,rs,ii
!
  do i=1,at(rs)%nbb+at(rs)%nsc
    if (i.le.at(rs)%nbb) then
      ii = at(rs)%bb(i)
    else
      ii = at(rs)%sc(i-at(rs)%nbb)
    end if
    xref(ii) = x(ii)
    yref(ii) = y(ii)
    zref(ii) = z(ii)
  end do
!
end
!
!---------------------------------------------------------------
!
!
subroutine makeref_forbb(imol,rs)
!
  use atoms
  use polypep
  use molecule
!
  implicit none
!
  integer i,rs,imol
!
! for a pivot move all atoms 'behind' active residue are moved as well
! note that this assumes that the first atom in the residue is a "backbone" atom 
  do i=at(rs)%bb(1),atmol(imol,2)
    xref(i) = x(i)
    yref(i) = y(i)
    zref(i) = z(i)
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_formol(imol)
!
  use atoms
  use molecule
!
  implicit none
!
  integer i,imol
!
! all atoms in the molecule are backed up
  do i=atmol(imol,1),atmol(imol,2)
    xref(i) = x(i)
    yref(i) = y(i)
    zref(i) = z(i)
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_poly(imol)
!
  use molecule
!
  implicit none
!
  integer i,imol,j
!
  rgvref(imol) = rgv(imol)
  do i=1,3
    rgevsref(imol,i) = rgevs(imol,i)
    comref(imol,i) = com(imol,i)
  end do
  do i=1,3
    do j=1,3
      rgpcsref(imol,i,j) = rgpcs(imol,i,j)
    end do
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_polym(imol)
!
  use molecule
!
  implicit none
!
  integer imol
!
!  rgvmref(imol) = rgvm(imol)
!  rgevsmref(imol,1:3) = rgevsm(imol,1:3)
  commref(imol,1:3) = comm(imol,1:3)
!  do i=1,3
!    do j=1,3
!      rgpcsmref(imol,i,j) = rgpcsm(imol,i,j)
!    end do
!  end do
!
end
!
!---------------------------------------------------------------
!
subroutine makeref_zmat_gl()
!
  use atoms
  use zmatrix
!
  implicit none
!
  ztorpr(1:n) = ztor(1:n)
  bangpr(1:n) = bang(1:n)
!
end
!
!-----------------------------------------------------------------------
!
subroutine getref_forrotlst(ati,alC,tpi)
!
  use atoms
  use zmatrix
  use sequen
  use molecule
  use atoms
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,ati
  logical, INTENT(IN):: alC
!
  integer ii,ik,k,imol,shf,shf2
#ifdef ENABLE_THREADS
  integer stx(2),tpn
!
  tpn = thrdat%maxn
#endif
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
!
  if (alC.EQV..true.) then
    shf = 1
    shf2 = 1
    imol = molofrs(atmres(ati))
    if (izrot(ati)%rotis(1,1).eq.atmol(imol,1)) shf = 2
    if (izrot(ati)%rotis(izrot(ati)%alsz,2).eq.atmol(imol,2)) shf2 = 0
    do k=shf,izrot(ati)%alsz+shf2
      if (k.eq.1) then 
        ii = atmol(imol,1)
      else
        ii = izrot(ati)%rotis(k-1,2)+1
      end if
      if (k.gt.izrot(ati)%alsz) then 
        ik = atmol(imol,2)
      else
        ik = izrot(ati)%rotis(k,1)-1
      end if
#ifdef ENABLE_THREADS
      if ((ik-ii+1).gt.10*thrdat%maxn) then
        call threads_bounds(ii,ik,tpi,tpn,stx)
        x(stx(1):stx(2)) = xref(stx(1):stx(2))
        y(stx(1):stx(2)) = yref(stx(1):stx(2))
        z(stx(1):stx(2)) = zref(stx(1):stx(2))
      else
!$OMP SINGLE
#endif
      x(ii:ik) = xref(ii:ik)
      y(ii:ik) = yref(ii:ik)
      z(ii:ik) = zref(ii:ik)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end do
  else
    do k=1,izrot(ati)%alsz
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
#ifdef ENABLE_THREADS
      if ((ik-ii+1).gt.10*thrdat%maxn) then
        call threads_bounds(ii,ik,tpi,tpn,stx)
        x(stx(1):stx(2)) = xref(stx(1):stx(2))
        y(stx(1):stx(2)) = yref(stx(1):stx(2))
        z(stx(1):stx(2)) = zref(stx(1):stx(2))
      else
!$OMP SINGLE
#endif
      x(ii:ik) = xref(ii:ik)
      y(ii:ik) = yref(ii:ik)
      z(ii:ik) = zref(ii:ik)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
      end if
#endif
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine getref_forsc(rs)
!
  use atoms
  use polypep
!
  implicit none
!
  integer i,rs,ii
!
  do i=1,at(rs)%nbb+at(rs)%nsc
    if (i.le.at(rs)%nbb) then
      ii = at(rs)%bb(i)
    else
      ii = at(rs)%sc(i-at(rs)%nbb)
    end if
    x(ii) = xref(ii)
    y(ii) = yref(ii)
    z(ii) = zref(ii)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_forbb(imol,rs)
!
  use atoms
  use polypep
  use molecule
!
  implicit none
!
  integer i,rs,imol
!
  do i=at(rs)%bb(1),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_forbb_nuccr(imol,rs)
!
  use atoms
  use polypep
  use molecule
!
  implicit none
!
  integer i,rs,imol
!
  do i=nuci(rs,1),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_forbb_nuccr_pre(imol,rs)
!
  use atoms
  use polypep
  use molecule
  use system
!
  implicit none
!
  integer i,rs,imol
!
! rsf-1 should never be a 5'-cap residue
  do i=nuci(rs,6),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_formol(imol)
!
  use atoms
  use molecule
!
  implicit none
!
  integer i,imol
!
! all atoms in the molecule are restored
  do i=atmol(imol,1),atmol(imol,2)
    x(i) = xref(i)
    y(i) = yref(i)
    z(i) = zref(i)
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine getref_poly(imol)
!
  use molecule
!
  implicit none
!
  integer i,imol,j
!
  rgv(imol) = rgvref(imol)
  do i=1,3
    rgevs(imol,i) = rgevsref(imol,i)
    com(imol,i) = comref(imol,i)
  end do
  do i=1,3
    do j=1,3
      rgpcs(imol,i,j) = rgpcsref(imol,i,j)
    end do
  end do
!
end
!
!---------------------------------------------------------------
!
subroutine getref_zmat_gl()
!
  use atoms
  use zmatrix
!
  implicit none
!
  ztor(1:n) = ztorpr(1:n)
  bang(1:n) = bangpr(1:n)
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_foruj(imol,rs)
!
  use atoms
  use iounit
  use polypep
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,imol
!
  if (seqpolty(rs).eq.'P') then
    do i=at(rs)%bb(3),atmol(imol,2)
      x(i) = xref(i)
      y(i) = yref(i)
      z(i) = zref(i)
    end do
  else
    write(ilog,*) 'Fatal. Called getref_foruj(...) with unsupported &
 &polymer or cap type. Please report this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_fordjo(imol,rs)
!
  use atoms
  use iounit
  use polypep
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,imol
!
  if (seqpolty(rs).eq.'P') then
    x(at(rs)%bb(4)) = xref(at(rs)%bb(4))
    y(at(rs)%bb(4)) = yref(at(rs)%bb(4))
    z(at(rs)%bb(4)) = zref(at(rs)%bb(4))
    do i=at(rs+1)%bb(1),atmol(imol,2)
      x(i) = xref(i)
      y(i) = yref(i)
      z(i) = zref(i)
    end do
  else
    write(ilog,*) 'Fatal. Called getref_fordjo(...) with unsupported &
 &polymer or cap type. Please report this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------
!
subroutine getref_fordo(imol,rs)
!
  use atoms
  use iounit
  use polypep
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,imol
!
  if (seqpolty(rs).eq.'P') then
    do i=at(rs)%bb(2),atmol(imol,2)
      x(i) = xref(i)
      y(i) = yref(i)
      z(i) = zref(i)
    end do
  else
    write(ilog,*) 'Fatal. Called getref_fordo(...) with unsupported &
 &polymer or cap type. Please report this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------
!

