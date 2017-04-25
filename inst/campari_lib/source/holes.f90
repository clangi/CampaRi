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
! MAIN AUTHOR:   Hoang Tran                                                !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
!
#include "macros.i"
!
! Hoang Tran, 22 Aug., 2005
! calculate if a sphere of size 'da' placed at a specific 
! position will intersect with the chain
!
!
! wrapper routine for void space analysis
!
!
subroutine holes()
!
  use mcsums
  use sequen
  use molecule
!
  implicit none
!
  logical empty
  integer isempty
  RTYPE estrg,xcm,ycm,zcm,rsx,rsy,rsz,ru(3),rbegin,abegin
  RTYPE rstep,rstop,rnow,astep,astop,anow
  integer numa,numr

!
! estimate IPP radius of gyration
  estrg = (0.8*molcontlen(1))**0.6
!
! estimate a value to use for the upper limits (rstop,astop),
! based on IPP good solvent scaling
  rstop = estrg*2.0
  astop = estrg/2.0
  rstep = 5.0d0
  astep = 2.5d0
  rbegin = 0.0d0
  abegin = astep

  numa = 0
  numr = 0

  if (firsthole) then
     anow = astep
     astop = max(astop,astep)
     do while (anow .le. astop)
        numa = numa + 1
        anow = anow + astep
     end do
     rnow = 0.0d0
     do while (rnow .le. rstop)
        numr = numr + 1
        rnow = rnow + rstep
     end do
     write(iholes,10)rbegin,rstep,numr
 10      format('Rstart, Rstep, Numr = ',f6.2,f6.2,i4)
     write(iholes,20)abegin,astep,numa
 20      format('Astart, Astep, Numa = ',f6.2,f6.2,i4)
     firsthole = .false.
  end if

!cccccccccc
! get the center of mass
  call gcom(xcm,ycm,zcm)
! get a random unit vector
  call ranvec(ru)
!cccccccccccc
!
! initialize 
  anow = abegin
  do while (anow .le. astop)
     rnow = rbegin
     do while (rnow .le. rstop)
!       set position of sphere
        rsx = ru(1)*rnow + xcm
        rsy = ru(2)*rnow + ycm
        rsz = ru(3)*rnow + zcm
        call spherefits(rsx,rsy,rsz,anow,empty)
        if (empty.EQV..true.) then
           isempty = 1
        else
           isempty = 0
        end if
        write(iholes,30,advance="no") isempty
 30         format(i2)
        rnow = rnow + rstep
     end do
     write(iholes,40,advance="no") '  '
 40      format(a2)
     anow = anow + astep
  end do
  write(iholes,*)

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccc
!     
! compute (global) center of mass of system
subroutine gcom(xcm,ycm,zcm)
!
  use atoms
!
  implicit none
!
  integer i
  RTYPE total,xcm,ycm,zcm,weigh
!
! compute the position of the center of mass
!
  total = 0.0d0
  xcm = 0.0d0
  ycm = 0.0d0
  zcm = 0.0d0
  do i = 1, n
    weigh = mass(i)
    total = total + weigh
    xcm = xcm + x(i)*weigh
    ycm = ycm + y(i)*weigh
    zcm = zcm + z(i)*weigh
  end do
  xcm = xcm / total
  ycm = ycm / total
  zcm = zcm / total

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! computes whether a sphere of radius 'da' can fit at the
! specified position

subroutine spherefits(rsx,rsy,rsz,da,fits)
!
  use atoms
  use polypep
  use sequen
  use molecule
  use system
  use grandensembles
!
  implicit none
!  
  integer i,ii,j,jj,imol
  RTYPE rsx,rsy,rsz,da,dxr,dyr,dzr
  RTYPE dx,dy,dz,resdist,dist,room
  logical fits,ismember

  fits = .true.

!
! loop over all residues
  do imol=1,nmol
!   cycle non-presents (should not apply)
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
    do i=rsmol(imol,1),rsmol(imol,2)
     ii = refat(i)
     dxr = x(ii) - rsx
     dyr = y(ii) - rsy
     dzr = z(ii) - rsz
     resdist = sqrt(dxr*dxr + dyr*dyr + dzr*dzr)
     if (resdist .lt. (resrad(i) + da + 10.0d0)) then
!
! loop over the atoms within this residue
        do j=1,at(i)%nbb+at(i)%nsc
           if (j .le. at(i)%nbb) then
              jj = at(i)%bb(j)
           else
              jj = at(i)%sc(j-at(i)%nbb)
           end if
           dx = x(jj) - rsx
           dy = y(jj) - rsy
           dz = z(jj) - rsz
           dist = sqrt(dx*dx + dy*dy + dz*dz)
           room = da + atr(jj)
           if (dist .ge. room) then
              continue
           else
              fits = .false.
              exit
           end if
        end do
     end if
     if (.not. fits) exit
    end do
  end do
!
end
!
