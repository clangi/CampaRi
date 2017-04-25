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
! CONTRIBUTIONS: Xiaoling Wang                                             !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! "readfyc" read in torsions to re-build a single polymer chain
! from the proper torsional output file
! this is terribly non-general and heavily discouraged
!
subroutine readfyc()
!
  use fyoc
  use sequen
  use iounit
  use torsn
  use system
  use molecule
  use aminos
  use polypep
  use zmatrix
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer fyc,njunk
  integer i,freeunit,iomess,t1,t2,j,k,cnt,shf,shfbu
  RTYPE ejunk
  RTYPE, ALLOCATABLE:: tvec(:)
  character(10*ntorpuck+50), ALLOCATABLE:: fycl(:)
  logical exists
!
  if (nmol.ne.1) then
    write(ilog,*) 'Fatal. Reading in data for multiple molecules from an FYC-file&
 & is impossible. Please use a different way to specify structural input.'
    call fexit()
  end if
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  allocate(fycl(2))
  allocate(tvec(ntorpuck))
!
  call strlims(fycfile,t1,t2)
  inquire(file=fycfile(t1:t2),exist=exists)
  fyc = freeunit()
  if(exists.EQV..true.) then
    open (unit=fyc,file=fycfile(t1:t2),status='old')
    rewind(unit=fyc)
  else
    write(ilog,*) 'Fatal. Could not open FYC-file (',fycfile(t1:t2),&
 &'). Maybe it was (re)moved?'
    call fexit()
  end if
!
! read in the FYC structure file
! The structure specified in the last line of the file is used
 11   format(i12,1x,g15.6,1x,2000000(f9.3,1x))
  cnt = 0
  do while (.true.)
 10     format(a)
    read(fyc,10,iostat=iomess) fycl(1)
    if (iomess.eq.IOSTAT_END) exit
    if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing FYC-inpu&
 &t (got: ',fycfile(t1:t2),').'
      call fexit()
    end if
    cnt = cnt + 1
    if (cnt.eq.1) cycle ! discard header
    read(fycl(1),11,iostat=iomess) njunk,ejunk,(tvec(j),j=1,ntorpuck)
  end do
!
! now transfer the torsions from the temporary array to the pointer arrays
  j = 0
  do i=1,nseq
    if (wline(i).gt.0) then
      j = j + 1
      omega(i) = tvec(j)
    end if
    if ((fline(i).gt.0).AND.(seqflag(i).ne.5)) then
      j = j + 1
      phi(i) = tvec(j)
    end if
    if (yline(i).gt.0) then
      j = j + 1
      psi(i) = tvec(j)
    end if
    do k=1,nnucs(i)
      j = j + 1
      nucs(k,i) = tvec(j)
    end do
    do k=1,nchi(i)
      j = j + 1
      chi(k,i) = tvec(j)
    end do
!   direct Z-matrix edit for pucker
    if ((pucline(i).gt.0).OR.((seqpolty(i).eq.'N').AND.(nucsline(6,i).gt.0))) then
      if (seqpolty(i).eq.'N') then
        shfbu = shf
        shf = 1
      end if
      if (seqpolty(i).eq.'N') then
        j = j + 1
        ztor(nucsline(6,i)) = tvec(j)
      else
        j = j + 1
        ztor(fline(i)) = tvec(j)
        phi(i) = tvec(j) ! necessary, otherwise fyczmat will overwrite
      end if
      j = j + 1
      ztor(at(i)%sc(2-shf)) = tvec(j)
      j = j + 1
      ztor(at(i)%sc(3-shf)) = tvec(j)
      j = j + 1
      ztor(at(i)%sc(4-shf)) = tvec(j)
      j = j + 1
      bang(at(i)%sc(2-shf)) = tvec(j)
      j = j + 1
      bang(at(i)%sc(3-shf)) = tvec(j)
      j = j + 1
      bang(at(i)%sc(4-shf)) = tvec(j)
      if (seqpolty(i).eq.'N') shf = shfbu
    end if
  end do
! transfer pointer arrays to Z-matrix 
  call fyczmat()
! build Cartesian coordinates
  call makexyz2()
!
  deallocate(fycl)
  deallocate(tvec)
!
end
!
!----------------------------------------------------------------------
