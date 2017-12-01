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
! CONTRIBUTIONS: Xiaoling Wang, Rohit Pappu                                !
!                                                                          !
!--------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------
!
! this subroutine updates the real array (Z-coords) based on the settings
! in the pointer array (phi,psi,chi,...)
! note that phi and psi essentially control two internal degrees of freedom,
! i.e., two dihedrals, and hence their value must be imposed on both of them
! (even though one is strictly dependent)
! the fxn skips all dihedrals in the pointer array, for which no pointer is defined
! (i.e., for which ?line is zero)
!
subroutine fyczmat()
!
  use fyoc
  use sequen
  use zmatrix
!
  implicit none
!
  integer ff,f2,ll,yy,y2,ww
  integer i,j
!
!     
  do i = 1,nseq
!
! exclude inactive residues
    if (notors(i).EQV..true.) cycle
!
! first phi
    ff = fline(i)
    f2 = fline2(i)
    if (ff.gt.0) then
      ztor(ff) = phi(i)
    end if
!
! set other torsions that depend on phi
    if ((f2.gt.0).AND.(ff.gt.0)) then
      ztor(f2) = ztor(ff) + phish(i)
!      if (ztor(ff) .gt. 0.0d0) then
!        ztor(f2) = ztor(ff)-180.0d0
!      else
!        ztor(f2) = ztor(ff)+180.0d0
!      end if
      if (ztor(f2).gt.180.0) ztor(f2) = ztor(f2) - 360.0
      if (ztor(f2).lt.-180.0) ztor(f2) = ztor(f2) + 360.0
    end if
!
! next set psi
    yy= yline(i)
    y2 = yline2(i)
    if (yy.gt.0) then
      ztor(yy) = psi(i)
    end if
!
! followed by torsions that depend on psi
    if ((y2.gt.0).AND.(yy.gt.0)) then
!      ztor(y2) = ztor(yy)
      ztor(y2) = ztor(yy) + psish(i)
      if (ztor(y2).gt.180.0) ztor(y2) = ztor(y2) - 360.0
      if (ztor(y2).lt.-180.0) ztor(y2) = ztor(y2) + 360.0
    end if
!
! reset the psi-angle according to zmatrix format
    if (yy.gt.0) then
      if(ztor(yy) .lt. 0.0d0) then
        ztor(yy) = ztor(yy) + 180.0d0
      else
        ztor(yy) = ztor(yy) - 180.0d0
      end if
    end if
!
! put nucleic acid backbone angles into the z-matrix
    do j = 1,nnucs(i)
      ll = nucsline(j,i)
      ztor(ll) = nucs(j,i)
    end do
!
! copy the chi angles into the z-matrix
    do j = 1,nchi(i)
      ll = chiline(j,i)
      ztor(ll) = chi(j,i)
    end do
!
! omega angles
    ww = wline(i)
    if (ww.gt.0) then
      ztor(ww) = omega(i)
    end if
!
  end do
!
end
!
!-------------------------------------------------------------------------------
!
!
! zmatfyc is the inverse routine, it updates the pointer array based on the
! real array
! note that we essentially have a choice of how to assign phi and psi because of the
! 1:2-mapping (see above)
!
subroutine zmatfyc()
!
  use fyoc
  use sequen
  use zmatrix
!
  implicit none
!
  integer ff,ll,yy,ww
  integer i,j
!
!     
  do i = 1,nseq
!
! exclude inactive residues
    if (notors(i).EQV..true.) cycle
!
! first phi
    ff = fline(i)
    if (ff.gt.0) then
      phi(i) = ztor(ff)
      if (phi(i).gt.180.0) phi(i) = phi(i) - 360.0
      if (phi(i).lt.-180.0) phi(i) = phi(i) + 360.0
    end if
!
! next set psi
    yy = yline(i)
    if (yy.gt.0) then
      psi(i) = ztor(yy)
    end if
!
! reset the psi-angle according to zmatrix format
    if (yy.gt.0) then
      if (psi(i) .lt. 0.0d0) then
        psi(i) = psi(i) + 180.0d0
      else
        psi(i) = psi(i) - 180.0d0
      end if
      if (psi(i).gt.180.0) psi(i) = psi(i) - 360.0
      if (psi(i).lt.-180.0) psi(i) = psi(i) + 360.0
    end if
!
! get nucleic acid backbone angles from z-matrix
    do j = 1,nnucs(i)
      ll = nucsline(j,i)
      nucs(j,i) = ztor(ll)
      if (nucs(j,i).gt.180.0) nucs(j,i) = nucs(j,i) - 360.0
      if (nucs(j,i).lt.-180.0) nucs(j,i) = nucs(j,i) + 360.0
    end do
!
! get chi angles from z-matrix
    do j = 1,nchi(i)
      ll = chiline(j,i)
      chi(j,i) = ztor(ll)
      if (chi(j,i).gt.180.0) chi(j,i) = chi(j,i) - 360.0
      if (chi(j,i).lt.-180.0) chi(j,i) = chi(j,i) + 360.0
    end do
!
! omega angles
    ww = wline(i)
    if (ww.gt.0) then
      omega(i) = ztor(ww)
      if (omega(i).gt.180.0) omega(i) = omega(i) - 360.0
      if (omega(i).lt.-180.0) omega(i) = omega(i) + 360.0
    end if
!
  end do
!
end
!
!---------------------------------------------------------------------------------
!
! a wrapper updating psish and phish beforehand(!)
!
subroutine zmatfyc2()
!
  use molecule
  use sequen
  use fyoc
  use zmatrix
!
  implicit none
!
  integer imol,rs
!
  do imol=1,nmol
    do rs=rsmol(imol,1),rsmol(imol,2)
!     because of the way, phi- and psi-angles are coded up, there are two associated
!     ztor()-entries for each
!     in case the geometry comes from pdb we cannot assume what their (static) relationship
!     is and need to determine it explicitly
      if ((fline(rs).gt.0).AND.(fline2(rs).gt.0)) then 
        phish(rs) = ztor(fline2(rs)) - ztor(fline(rs))
        if (phish(rs).gt.180.0) phish(rs) = phish(rs) - 360.0
        if (phish(rs).le.-180.0) phish(rs) = phish(rs) + 360.0
      end if
      if ((yline(rs).gt.0).AND.(yline2(rs).gt.0)) then 
        if (ztor(yline(rs)).lt.0.0) then
          psish(rs) = ztor(yline2(rs)) - (ztor(yline(rs)) + 180.0)
        else
          psish(rs) = ztor(yline2(rs)) - (ztor(yline(rs)) - 180.0)
        end if
        if (psish(rs).gt.180.0) psish(rs) = psish(rs) - 360.0
        if (psish(rs).le.-180.0) psish(rs) = psish(rs) + 360.0
      end if
    end do
  end do
!
  call zmatfyc()
!
end
!
!-----------------------------------------------------------------------
!
subroutine zmatfyc_rs(i)
!
  use fyoc
  use sequen
  use zmatrix
!
  implicit none
!
  integer, INTENT(IN):: i
!
  integer ff,ll,yy,ww
  integer j     
!
! exclude inactive residues
  if (notors(i).EQV..true.) return
!
! first phi
  ff = fline(i)
  if (ff.gt.0) then
    phi(i) = ztor(ff)
    if (phi(i).gt.180.0) phi(i) = phi(i) - 360.0
    if (phi(i).lt.-180.0) phi(i) = phi(i) + 360.0
  end if
!
! next set psi
  yy = yline(i)
  if (yy.gt.0) then
    psi(i) = ztor(yy)
  end if
!
! reset the psi-angle according to zmatrix format
  if (yy.gt.0) then
    if (psi(i) .lt. 0.0d0) then
      psi(i) = psi(i) + 180.0d0
    else
      psi(i) = psi(i) - 180.0d0
    end if
    if (psi(i).gt.180.0) psi(i) = psi(i) - 360.0
    if (psi(i).lt.-180.0) psi(i) = psi(i) + 360.0
  end if
!
! get nucleic acid backbone angles from z-matrix
  do j = 1,nnucs(i)
    ll = nucsline(j,i)
    nucs(j,i) = ztor(ll)
    if (nucs(j,i).gt.180.0) nucs(j,i) = nucs(j,i) - 360.0
    if (nucs(j,i).lt.-180.0) nucs(j,i) = nucs(j,i) + 360.0
  end do
!
! get chi angles from z-matrix
  do j = 1,nchi(i)
    ll = chiline(j,i)
    chi(j,i) = ztor(ll)
    if (chi(j,i).gt.180.0) chi(j,i) = chi(j,i) - 360.0
    if (chi(j,i).lt.-180.0) chi(j,i) = chi(j,i) + 360.0
  end do
!
! omega angles
  ww = wline(i)
  if (ww.gt.0) then
    omega(i) = ztor(ww)
    if (omega(i).gt.180.0) omega(i) = omega(i) - 360.0
    if (omega(i).lt.-180.0) omega(i) = omega(i) + 360.0
  end if
!
end
!
!-----------------------------------------------------------------------------------------------
!

