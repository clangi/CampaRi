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


! a tool to correctly print very large hexidecimal numbers from binary masks with "sz" elements

program mask_gen

  use, INTRINSIC:: ISO_FORTRAN_ENV
!
  implicit none
!
  REAL(KIND=REAL128) deint,tmp,deint2
  integer i,k,kk,sz,lv,st,shf1,shf2,sta
  logical notdone,started
  character(100) hexn
  character hexchar
!
  sz = 72
!
  notdone = .true.
  lv=1
  do while (notdone.EQV..true.)
    if (mod(sz,lv).ne.0) then
      lv=lv+1
      if (lv.gt.sz) exit
      cycle
    end if
    if (lv.gt.sz) exit
    st=sz/lv
    do kk=1,lv
      deint = 0.
      do k=(kk-1)*st,kk*st-1
        tmp = 2.0_REAL128**k
        deint = deint + tmp
      end do
      hexn(:) = ' '
      tmp = deint+0.005 ! buffer the value so that the floating point comparisons all become safe except the last one
      started = .false.
      do i=ceiling(sz/4.0),0,-1
        if (i.eq.0) tmp = tmp - 0.005 ! remove buffer
        if (log(tmp).lt.(i*log(16.0_REAL128))) then
          if (started.EQV..true.) hexn((ceiling(sz/4.0)-i+1):(ceiling(sz/4.0)-i+1)) = '0'
          cycle
        end if
        started = .true.
        if (i.gt.0) then
          deint2 = exp(log(tmp)-(i*log(16.0_REAL128)))
          k = floor(deint2)
        else
          k = nint(tmp)
        end if
        if (k.gt.0) tmp = tmp - k*(16.0_REAL128**i)
        hexn((ceiling(sz/4.0)-i+1):(ceiling(sz/4.0)-i+1)) = hexchar(k)
      end do
      write(*,'(1x,i4,1x,i4,1x,a,g32.26)') (kk-1)*st,kk*st-1,hexn(1:(ceiling(sz/4.0)+2)),deint
    end do
    if ((sz/lv).ge.2) then
      do kk=1,lv
        deint = 0.
        if (mod(sz,2*lv).eq.0) then        
          st=sz/lv/2
          do k=(kk-1)*st,kk*st-1
            tmp = 2.0_REAL128**k
            deint = deint + tmp
          end do
          do k=sz/2+(kk-1)*st,sz/2+kk*st-1
            tmp = 2.0_REAL128**k
            deint = deint + tmp
          end do
        else
          sta = 1
          if (mod(kk,2).eq.0) then
            shf1 = 0
            shf2 = -1
          else
            shf1 = -1
            shf2 = 0
          end if
          st = floor(sz/lv/2.0)
          do k=nint((kk-1)*(st+0.5))+shf2,nint((kk-1)*(st+0.5))+st+shf1+shf2
            tmp = 2.0_REAL128**k
            deint = deint + tmp
          end do
          do k=sz/2+nint((kk-1)*(st+0.5)),sz/2+nint((kk-1)*(st+0.5))+st+shf2
            tmp = 2.0_REAL128**k
            deint = deint + tmp
          end do
        end if
        hexn(:) = ' '
        tmp = deint+0.1 ! buffer the value so that the floating point comparisons all become safe except the last one
        started = .false.
        do i=ceiling(sz/4.0),0,-1
          if (i.eq.0) tmp = tmp - 0.1 ! remove buffer
          if (log(tmp).lt.(i*log(16.0_REAL128))) then
            if (started.EQV..true.) hexn((ceiling(sz/4.0)-i+1):(ceiling(sz/4.0)-i+1)) = '0'
            cycle
          end if
          started = .true.
          if (i.gt.0) then
            deint2 = exp(log(tmp)-(i*log(16.0_REAL128)))
            k = floor(deint2)
          else
            k = nint(tmp)
          end if
          if (k.gt.0) tmp = tmp - k*(16.0_REAL128**i)
          hexn((ceiling(sz/4.0)-i+1):(ceiling(sz/4.0)-i+1)) = hexchar(k)
        end do
        if (mod(sz,2*lv).eq.0) then
          write(*,'(1x,i4,1x,i4,1x,i4,1x,i4,1x,a,g32.26)') (kk-1)*st,kk*st-1,sz/2+(kk-1)*st,sz/2+kk*st-1,&
 &hexn(1:(ceiling(sz/4.0)+2)),deint
        else
          write(*,'(1x,f6.1,1x,f6.1,1x,f6.1,1x,f6.1,1x,a,g32.26)') (kk-1)*(st+0.5),kk*(st+0.5)-1,sz/2+(kk-1)*(st+0.5),&
 &sz/2+kk*(st+0.5)-1,hexn(1:(ceiling(sz/4.0)+2)),deint
        end if
      end do
    end if
    lv=lv+1
  end do
!
end program
!
function hexchar(i)
!
  integer i
  character hexchar
!
  if (i.eq.0) hexchar = '0'
  if (i.eq.1) hexchar = '1'
  if (i.eq.2) hexchar = '2'
  if (i.eq.3) hexchar = '3'
  if (i.eq.4) hexchar = '4'
  if (i.eq.5) hexchar = '5'
  if (i.eq.6) hexchar = '6'
  if (i.eq.7) hexchar = '7'
  if (i.eq.8) hexchar = '8'
  if (i.eq.9) hexchar = '9'
  if (i.eq.10) hexchar = 'A'
  if (i.eq.11) hexchar = 'B'
  if (i.eq.12) hexchar = 'C'
  if (i.eq.13) hexchar = 'D'
  if (i.eq.14) hexchar = 'E'
  if (i.eq.15) hexchar = 'F'
!
end function hexchar 
