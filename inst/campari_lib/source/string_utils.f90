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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!     
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
!
!             ###############################
!             #                             #
!             #  A SET OF SUBROUTINES TO    #
!             #  PARSE AND EDIT STRINGS     #
!             #                             #
!             ###############################
!
!
!-------------------------------------------------------------------
!
! a subroutine to determine the limits for non-blank characters in
! an input string (f = first, l = last)
!
subroutine strlims(str,f,l)
!
  implicit none
!
  integer f,l,last,i,leng
  character(*) str
  character nix
!
  l = 0
  nix = char(0)
  leng = len(str)
  last = leng
  f = last+1
  do i=1,leng
    if (str(i:i).gt.' ') then
      f = i
      exit
    end if
  end do
  do i=1,leng
    if (str(i:i).eq.nix) then
      last = i - 1
      exit
    end if
  end do
!
  do i=last,1,-1
    if (str(i:i).gt.' ') then
      l = i
      exit
    end if
  end do
!
  if (f.gt.i) then
    write(*,*) 'WARNING: Got empty string in strlims(...).'
    f = 1
    i = 1
  end if
!
end
!
!-------------------------------------------------------------------
!
! the same, only silent
!
subroutine strlims_quiet(str,f,l)
!
  implicit none
!
  integer f,l,last,i,leng
  character(*) str
  character nix
!
  l = 0
  nix = char(0)
  leng = len(str)
  last = leng
  f = last+1
  do i=1,leng
    if (str(i:i).gt.' ') then
      f = i
      exit
    end if
  end do
  do i=1,leng
    if (str(i:i).eq.nix) then
      last = i - 1
      exit
    end if
  end do
!
  do i=last,1,-1
    if (str(i:i).gt.' ') then
      l = i
      exit
    end if
  end do
!
  if (f.gt.i) then
    f = 1
    i = 1
  end if
!
end
!
!-------------------------------------------------------------------
!
! similar routine operating on a left-truncated substring
!
subroutine strlims_ll(str,f,l,k)
!
  implicit none
!
  integer f,l,last,i,leng,k
  character(*) str
  character nix
!
  l = 0
  nix = char(0)
  leng = len(str)
  last = leng
  f = last+1
  do i=k,leng
    if (str(i:i).gt.' ') then
      f = i
      exit
    end if
  end do
  do i=k,leng
    if (str(i:i).eq.nix) then
      last = i - 1
      exit
    end if
  end do
!
  do i=last,k,-1
    if (str(i:i).gt.' ') then
      l = i
      exit
    end if
  end do
!
  if (f.gt.i) then
    write(*,*) 'WARNING: Got empty string in strlims(...).'
    f = 1
    i = 1
  end if
!
end
!
!---------------------------------------------------------------------
!
! this one is only used for getting keywords
! warning omitted
!
subroutine strlims_k(str,f,l)
!
  implicit none
!
  integer f,l,last,i,leng
  character str(200)
  character nix
!
!
  l = 0
  nix = char(0)
  leng = size(str)
  last = leng
  f = last+1
  do i=1,leng
    if (str(i).gt.' ') then
      f = i
      exit
    end if
  end do
  do i=1,leng
    if (str(i).eq.nix) then
      last = i - 1
      exit
    end if
  end do
!
  do i=last,1,-1
    if (str(i).gt.' ') then
      l = i
      exit
    end if
  end do
!
  if (f.gt.i) then
    f = 1
    i = 1
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this is a standard string operation
! note that the intrinsic ichar() uses numerical character values
! which are NOT necessarily ASCII codes (for that iachar())
!
subroutine toupper(string)
!
  implicit none
!
  integer i,leng,nasc
  character(*) string
!
  leng=len(string)
  do i=1,leng
    nasc=ichar(string(i:i))
!   the numerical range of standard lower-case letters
    if ((nasc.ge.97).AND.(nasc.le.122)) then
      string(i:i)=char(nasc-32)
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! to upper first in word, only blank/tab supported as separator
!
subroutine toupper_first(string)
!
  implicit none
!
  integer i,leng,nasc
  character(*) string
!
  leng=len(string)
  do i=1,leng
    nasc=ichar(string(i:i))
!   the numerical range of standard lower-case letters
    if (i.gt.1) then
      if ((nasc.ge.97).AND.(nasc.le.122).AND.((string((i-1):(i-1)).eq.' ').OR.(string((i-1):(i-1)).eq.'\t'))) then
        string(i:i)=char(nasc-32)
      end if
    else
      if ((nasc.ge.97).AND.(nasc.le.122)) then
        string(i:i)=char(nasc-32)
      end if
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! analogous
!
subroutine tolower(string)
!
  implicit none
!
  integer i,leng,nasc
  character(*) string
!
  leng=len(string)
  do i=1,leng
    nasc=ichar(string(i:i))
!   the numerical range of standard upper-case letters
    if ((nasc.ge.65).AND.(nasc.le.90)) then
      string(i:i)=char(nasc+32)
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this routine searches for an integer number in an input string
! makes - like all other string manipulation functions - assumptions
! about ASCII-character ordering
! input "cont" is the index where to start processing
! output "cont" is the index where a non-numerical sequence resumes
!
subroutine extract_int(string,vali,cont)
!
  implicit none
!
  integer i,j,cont,t1,t2,s1,s2
  integer vali,tenex,ten2s(10)
  logical isminus,isnum
  character letter
  character(*) string
  data ten2s  /1,10,100,1000,10000,100000,1000000,10000000,&
 &              100000000,1000000000/
!
  isminus = .false.
  isnum = .false.
  call strlims_ll(string,t1,t2,cont)
!
  s1 = 1
  s2 = 0
  do i=cont,t2
    letter = string(i:i)
    if ((letter.ge.'0').AND.(letter.le.'9')) then
      if (isnum.EQV..false.) then
        isnum = .true.
        s1 = i
      end if
      if (i.eq.t2) then
        s2 = t2
        cont = i+1
      end if
    else if ((letter.eq.'-').AND.(isminus.EQV..false.)) then
      isminus = .true.
    else if (isnum.EQV..true.) then
      if ((letter.le.' ').OR.(letter.eq.',').OR.&
 &         (letter.eq.';').OR.(letter.eq.'_')) then
        s2 = i - 1
        cont = i
        exit
      end if
    else if ((letter.gt.' ').OR.(isminus.EQV..true.)) then
      exit
    end if
  end do
!
  if (s2.gt.(s1+9)) then
    write(*,*) 'Warning. Processing integer which is too large in&
 & extract_int(...). Truncating ...'
    s2 = s1+9
  end if
!
  j = 0
  tenex = 0
  do i=s2,s1,-1
    j = j + 1
    if (string(i:i).eq.'0') then
      tenex = 0
    else if (string(i:i).eq.'1') then
      tenex = 1
    else if (string(i:i).eq.'2') then
      tenex = 2
    else if (string(i:i).eq.'3') then
      tenex = 3
    else if (string(i:i).eq.'4') then
      tenex = 4
    else if (string(i:i).eq.'5') then
      tenex = 5
    else if (string(i:i).eq.'6') then
      tenex = 6
    else if (string(i:i).eq.'7') then
      tenex = 7
    else if (string(i:i).eq.'8') then
      tenex = 8
    else if (string(i:i).eq.'9') then
      tenex = 9
    else
      write(*,*) 'Fatal. Found non-numerical string in extract_int(...) after string filtering. The string in question is:'
      write(*,*) string(t1:t2)
      write(*,*) 'This is usually indicative of a bad input file or an error in the key-file (mismatched data type).'
      call fexit()
    end if
    if (i.eq.s2) then
      vali = tenex*ten2s(j)
    else
      vali = vali + tenex*ten2s(j)
    end if
  end do
!
  if (s2.lt.s1) then
    if ((t2.gt.t1).OR.(string(t1:t1).gt.' ')) then
      write(*,*) 'Warning. Failed to extract integer from first entry in string in extract_int(..). The string in question is:'
      write(*,*) string(t1:t2)
    end if
  else if (isminus.EQV..true.) then
    vali = -vali
  end if
!
end
!
!---------------------------------------------------------------------
!
! read consecutive entries of numbers
!
subroutine get_ints_from_str(string,nvals,valv,cont,nvalsread,forcer)
!
  implicit none
!
  character(*) string
  character(len(string)) string2
  integer nvals,cont,cont2,kk,iomess,nvalsread,ll
  integer valv(nvals)
  logical done,isnum,forcer
!
  nvalsread = 0
  done = .false.
  if (nvals.le.0) done = .true.
  cont2 = cont
  do while (done.EQV..false.)
    call extract_str(string,string2,cont2)
    if (cont2.eq.cont) then ! this means string is over 
      done = .true.
    else
      isnum = .false.
      ll = 1
      if ((string2(1:1).eq.'-').AND.(cont.lt.(cont2-1))) then
        ll = 2
      end if
      if ((string2(ll:ll).eq.'0').OR.(string2(ll:ll).eq.'1').OR.(string2(ll:ll).eq.'2').OR.(string2(ll:ll).eq.'3').OR.&
 &        (string2(ll:ll).eq.'4').OR.(string2(ll:ll).eq.'5').OR.(string2(ll:ll).eq.'6').OR.(string2(ll:ll).eq.'7').OR.&
 &        (string2(ll:ll).eq.'8').OR.(string2(ll:ll).eq.'9')) isnum = .true.
      if (isnum.EQV..true.) then
        read(string2,*,iostat=iomess) kk
        if (iomess.eq.0) then
          nvalsread = nvalsread + 1
          valv(nvalsread) = kk
          if (nvalsread.eq.nvals) done = .true.
        else ! not a valid integer even though it looks like one in the beginning
          if (forcer.EQV..false.) done = .true. ! exit unless we are force-reading
        end if
      else ! not a valid integer
        if (forcer.EQV..false.) done = .true. ! exit unless we are force-reading
      end if
    end if
    cont = cont2
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine get_reals_from_str(string,nvals,valv,cont,nvalsread,forcer)
!
  implicit none
!
  character(*) string
  character(len(string)) string2
  integer nvals,cont,cont2,iomess,nvalsread,ll
  RTYPE valv(nvals)
  RTYPE rr
  logical done,isnum,forcer
!
  nvalsread = 0
  done = .false.
  if (nvals.le.0) done = .true.
  cont2 = cont
  do while (done.EQV..false.)
    call extract_str(string,string2,cont2)
    if (cont2.eq.cont) then ! this means string is over 
      done = .true.
    else
      isnum = .false.
      ll = 1
      if ((string2(1:1).eq.'-').AND.(cont.lt.(cont2-1))) then
        ll = 2
      end if
!     negative number support 
      if ((string2(ll:ll).eq.'0').OR.(string2(ll:ll).eq.'1').OR.(string2(ll:ll).eq.'2').OR.(string2(ll:ll).eq.'3').OR.&
 &        (string2(ll:ll).eq.'4').OR.(string2(ll:ll).eq.'5').OR.(string2(ll:ll).eq.'6').OR.(string2(ll:ll).eq.'7').OR.&
 &        (string2(ll:ll).eq.'8').OR.(string2(ll:ll).eq.'9')) isnum = .true.
      if (string2(ll:ll).eq.'.') then
!       support for pos. or neg. numbers starting with fractional dot (0 omitted)
        if (cont.lt.(cont2-2)) then
          ll = ll + 1
          if ((string2(ll:ll).eq.'0').OR.(string2(ll:ll).eq.'1').OR.(string2(ll:ll).eq.'2').OR.(string2(ll:ll).eq.'3').OR.&
 &            (string2(ll:ll).eq.'4').OR.(string2(ll:ll).eq.'5').OR.(string2(ll:ll).eq.'6').OR.(string2(ll:ll).eq.'7').OR.&
 &            (string2(ll:ll).eq.'8').OR.(string2(ll:ll).eq.'9')) isnum = .true.
        end if
      end if
      if (isnum.EQV..true.) then
        read(string2,*,iostat=iomess) rr
        if (iomess.eq.0) then
          nvalsread = nvalsread + 1
          valv(nvalsread) = rr
          if (nvalsread.eq.nvals) done = .true.
        else ! not a valid number even though it looks like one in the beginning
          if (forcer.EQV..false.) done = .true. ! exit unless we are force-reading
        end if
      else ! definitely not a number
        if (forcer.EQV..false.) done = .true. ! exit unless we are force-reading
      end if
    end if
    cont = cont2
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine get_floats_from_str(string,nvals,valv,cont,nvalsread,forcer)
!
  implicit none
!
  character(*) string
  character(len(string)) string2
  integer nvals,cont,cont2,iomess,nvalsread,ll
  real(KIND=4) valv(nvals)
  real(KIND=4) rr
  logical done,isnum,forcer
!
  nvalsread = 0
  done = .false.
  if (nvals.le.0) done = .true.
  cont2 = cont
  do while (done.EQV..false.)
    call extract_str(string,string2,cont2)
    if (cont2.eq.cont) then ! this means string is over 
      done = .true.
    else
      isnum = .false.
      ll = 1
      if ((string2(1:1).eq.'-').AND.(cont.lt.(cont2-1))) then
        ll = 2
      end if
      if ((string2(ll:ll).eq.'0').OR.(string2(ll:ll).eq.'1').OR.(string2(ll:ll).eq.'2').OR.(string2(ll:ll).eq.'3').OR.&
 &        (string2(ll:ll).eq.'4').OR.(string2(ll:ll).eq.'5').OR.(string2(ll:ll).eq.'6').OR.(string2(ll:ll).eq.'7').OR.&
 &        (string2(ll:ll).eq.'8').OR.(string2(ll:ll).eq.'9')) isnum = .true.
      if (string2(ll:ll).eq.'.') then
!       support for pos. or neg. numbers starting with fractional dot (0 omitted)
        if (cont.lt.(cont2-2)) then
          ll = ll + 1
          if ((string2(ll:ll).eq.'0').OR.(string2(ll:ll).eq.'1').OR.(string2(ll:ll).eq.'2').OR.(string2(ll:ll).eq.'3').OR.&
 &            (string2(ll:ll).eq.'4').OR.(string2(ll:ll).eq.'5').OR.(string2(ll:ll).eq.'6').OR.(string2(ll:ll).eq.'7').OR.&
 &            (string2(ll:ll).eq.'8').OR.(string2(ll:ll).eq.'9')) isnum = .true.
        end if
      end if
      if (isnum.EQV..true.) then
        read(string2,*,iostat=iomess) rr
        if (iomess.eq.0) then
          nvalsread = nvalsread + 1
          valv(nvalsread) = rr
          if (nvalsread.eq.nvals) done = .true.
        else
          if (forcer.EQV..false.) done = .true.
        end if
      else
        if (forcer.EQV..false.) done = .true.
      end if
    end if
    cont = cont2
  end do
!
end
!
!------------------------------------------------------------------------
!
! a slightly modified version of extract_int that gives a little more information
! and reads through non-blank characters as well, it doesn't accept
! negative numbers
! (used for automatic parsing of systematically numbered input files only)
!
subroutine extract_int2(string,vali,cont,s12,leadzero)
!
  implicit none
!
  integer i,j,cont,t1,t2,s1,s2,s12,vali,leadzero
  integer tenex,ten2s(10)
  logical isnum
  character letter
  character(*) string
  data ten2s /1,10,100,1000,10000,100000,1000000,10000000,&
 &              100000000,1000000000/
!
  isnum = .false.
  call strlims_ll(string,t1,t2,cont)
  s1 = 1
  s2 = 0
  do i=cont,t2
    letter = string(i:i)
    if ((letter.ge.'0').AND.(letter.le.'9')) then
      if (isnum.EQV..false.) then
        isnum = .true.
        s1 = i
      end if
      if (i.eq.t2) then
        s2 = t2
        cont = i+1
      end if
    else if (isnum.EQV..true.) then
      s2 = i - 1
      cont = i
      exit
    end if
  end do
!
  if (s2.gt.(s1+9)) then
    write(*,*) 'Warning. Processing integer which is too large in&
 & extract_int2(...). Truncating ...'
    s2 = s1+9
  end if

  s12 = s1
  if ((string(s1:s1).eq.'0').AND.(s2.gt.s1)) then
    leadzero = s2-s1+1
  else
    leadzero = 0
  end if
!
  j = 0
  vali = 0
  tenex = 0
  do i=s2,s1,-1
    j = j + 1
    if (string(i:i).eq.'0') then
      tenex = 0
    else if (string(i:i).eq.'1') then
      tenex = 1
    else if (string(i:i).eq.'2') then
      tenex = 2
    else if (string(i:i).eq.'3') then
      tenex = 3
    else if (string(i:i).eq.'4') then
      tenex = 4
    else if (string(i:i).eq.'5') then
      tenex = 5
    else if (string(i:i).eq.'6') then
      tenex = 6
    else if (string(i:i).eq.'7') then
      tenex = 7
    else if (string(i:i).eq.'8') then
      tenex = 8
    else if (string(i:i).eq.'9') then
      tenex = 9
    else
      write(*,*) 'Fatal. Found non-numerical string in extract_int2(.&
 &..) after string filtering. This is a bug.'
      call fexit()
    end if
    vali = vali + tenex*ten2s(j)
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine extract_int3(string,t1,t2,toout)
!
  implicit none
!
  integer i,j,t1,t2,neg
  character letter
  character(*) string,toout
  logical innum
!
  j = len(toout)
  do i=1,j
    toout(i:i) = ' '
  end do
!
  neg = 0
  innum = .false.
  do i=t1,t2
    letter = string(i:i)
    if (letter.ge.'0' .and. letter.le.'9') then
      if (innum.EQV..false.) then
        if (neg.gt.0) then
          if (neg.eq.(i-1)) then
            innum = .true.
          end if
        else
          innum = .true.
        end if
      end if
      toout(i:i) = letter
    else if (letter.eq.'-') then
      neg = i
      toout(i:i) = letter
    else if (letter.eq.'.') then

    else
      if (innum.EQV..true.) then
        innum = .false.
        exit
      else
        if (neg.gt.0) then
          if (neg.eq.(i-1)) then
            toout(i-1:i-1) = ' '
          end if
          neg = 0
        end if
        toout(i:i) = ' '
      end if
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this routine simply copies the next non-blank string
! from the input ("string") starting to read at "cont".
! the isolated string is returned in "strout" and the "cont"-
! counter is set to the next blank
!
subroutine extract_str(string,strout,cont)
!
  implicit none
!
  integer insz,outsz,fi,st,cont,i
  character(*) string,strout
!
  insz = len(string(cont:))
  outsz = len(strout)
  if (outsz.le.0) then
    write(*,*) 'Fatal. Called extract_str(...) with unallocated stri&
 &ng input. This is a bug.'
    call fexit()
  end if
  if (insz.le.0) then
!   nothing to do
    do i=1,outsz
      strout(i:i) = ' '
    end do
    return
  end if
!
  st = len(string) + 1
  fi = len(string)
  do i=cont,len(string)
    if (string(i:i).gt.' ') then
      st = i
      exit
    end if
  end do
  do i=st+1,len(string)
    if (string(i:i).le.' ') then
      fi = i-1
      exit
    end if
  end do
!
  if (fi.lt.st) then
    do i=1,outsz
      strout(i:i) = ' '
    end do
  else if ((fi-st+1).gt.outsz) then
    cont = fi + 1
    fi = st + outsz - 1
    strout(1:outsz) = string(st:fi)
  else
    strout(1:(fi-st+1)) = string(st:fi)
    do i=fi-st+2,outsz
      strout(i:i) = ' '
    end do
    cont = fi + 1
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this routine simply copies the next a??-character string
! from the input ("string") starting to read at "cont".
! the isolated string is returned in "strout" and the "cont"-
! counter is set to the next blank
! (note that the target string must start with a standard letter,
!  but can otherwise contain any alpha-numerical character)
!
subroutine extract_abc(string,strout,cont)
!
  implicit none
!
  integer insz,outsz,cont,st,fi,nasc,i
  character(*) string,strout
!
  insz = len(string(cont:))
  outsz = len(strout)
  if (outsz.le.0) then
    write(*,*) 'Fatal. Called extract_abc(...) with unallocated stri&
 &ng input. This is a bug.'
    call fexit()
  end if
  if (insz.le.0) then
!   nothing to do
    do i=1,outsz
      strout(i:i) = ' '
    end do
    return
  end if
!
  st = len(string) + 1
  fi = len(string)
  do i=cont,len(string)
    nasc = ichar(string(i:i))
    if (((nasc.ge.65).AND.(nasc.le.90)).OR.&
 &      ((nasc.ge.97).AND.(nasc.le.122))) then
      st = i
      exit
    end if
  end do
  do i=st+1,len(string)
    if (string(i:i).le.' ') then
      fi = i-1
      exit
    end if
  end do
!
  if (fi.lt.st) then
    do i=1,outsz
      strout(i:i) = ' '
    end do
  else if ((fi-st+1).gt.outsz) then
    cont = fi + 1
    fi = st + outsz - 1
    strout(1:outsz) = string(st:fi)
  else
    strout(1:(fi-st+1)) = string(st:fi)
    do i=fi-st+2,outsz
      strout(i:i) = ' '
    end do
    cont = fi + 1
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this routine simply copies the next double-quoted string
! from the input ("string") starting to read at "cont".
! the isolated string is returned in "strout" and the "cont"-
! counter is set to the next blank
!
subroutine extract_quo(string,strout,cont)
!
  implicit none
!
  integer insz,outsz,cont,st,fi,i
  character(*) string,strout
!
  insz = len(string(cont:))
  outsz = len(strout)
  if (outsz.le.0) then
    write(*,*) 'Fatal. Called extract_quo(...) with unallocated stri&
 &ng input. This is a bug.'
    call fexit()
  end if
  if (insz.le.0) then
!   nothing to do
    do i=1,outsz
      strout(i:i) = ' '
    end do
    return
  end if
!
  st = len(string) + 1
  fi = len(string)
  do i=cont,len(string)
    if (string(i:i).eq.'"') then
      st = i+1
      exit
    end if
  end do
  do i=st+1,len(string)
    if (string(i:i).eq.'"') then
      fi = i-1
      exit
    end if
  end do
!
  if (fi.lt.st) then
    do i=1,outsz
      strout(i:i) = ' '
    end do
  else if ((fi-st+1).gt.outsz) then
    cont = fi + 2
    fi = st + outsz - 1
    strout(1:outsz) = string(st:fi)
  else
    strout(1:(fi-st+1)) = string(st:fi)
    do i=fi-st+2,outsz
      strout(i:i) = ' '
    end do
    cont = fi + 2
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this routine converts the integer in "ii" to a string (in "string")
! range exception handling is provided
! if the fxn is called with a string of size "strsz" left-padding with zeros
! is obtained
! if the fxn is called with "strsz" = 0 right-justification withe left-padding
! blanks is obtained
! if the fxn is called with a string of size larger than "strsz" right-padding
! with blanks is obtained
!
subroutine int2str(ii,string,strsz)
!
  implicit none
!
  integer mx(0:16),leng,i,ii,strsz,minsz,cpy,sumi
  character it(0:9)
  character(*) string
  logical just_r
  data it /'0','1','2','3','4','5','6','7','8','9'/
!
  if (strsz.le.0) then
    just_r = .true.
    strsz = 1
  else
    just_r = .false.
  end if
  minsz = strsz
  leng = len(string)
!
! invert sign if necessary (technically this routines only
! processes unsigned integers)
  if (ii .lt. 0) then
    i = -ii
  end if
!
! this is standard: the mx(...) are integers of course
!
  sumi = 0
  do i=16,0,-1
    mx(i) = (ii-sumi)/(10.0**i)
    if (i.gt.0) then
      sumi = sumi + (10.0**i) * mx(i)
    end if
  end do
!
! now we have the digits of our string in integer form (0 through 9)
!
! check for integer-size
!
  if ((mx(10).gt.0).OR.(mx(11).gt.0).OR.(mx(12).gt.0).OR.&
 &    (mx(13).gt.0).OR.&
 &    (mx(14).gt.0).OR.(mx(15).gt.0).OR.(mx(16).gt.0)) then
    write(*,*) 'Warning. Possibly bad result from int2str(...) (inte&
 &ger overflow).'
  end if
!
! find the final string length
!
  strsz = 1
  do i=16,1,-1
    if (mx(i).ne.0) then
      strsz = i+1
      exit
    end if
  end do
!
! correct if too large or small
  strsz = min(strsz,leng)
  strsz = max(strsz,minsz)
!
  if (strsz.gt.17) then
    write(*,*) 'Warning. Possibly bad result from int2str(...) (digi&
 &t overflow).'
  end if
!
! convert individual digits to a string of numeric characters
!
  do i=1,strsz
    string(i:i) = it(mx(strsz-i))
  end do
!
! left/right-justification
!
  if (just_r.EQV..true.) then
    do i=strsz,1,-1
      cpy = leng-strsz+i
      string(cpy:cpy) = string(i:i)
    end do
    do i=1,leng-strsz
      string(i:i) = ' '
    end do
  else
    do i = strsz+1,leng
      string(i:i) = ' '
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!


