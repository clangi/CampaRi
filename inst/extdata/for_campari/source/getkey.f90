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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------
!
! this subroutine scans a text-file and stores its contents in a dynamic
! data structure (see keys.i)
!
#include"macros.i"

subroutine getkey()
!
  use commline
  use iounit
  use keys
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer i,iu,t1,t2,t3,t4,j
  integer freeunit,iomess
  character(MAXARGLEN) keyfile
  character(MAXARGLEN) argstr
  character(MAXKEYLEN) linec
  character linea(MAXKEYLEN)
  logical exists
!
! scan command line arguments for key-file
!
  exists = .false.
  do i=1,nargs
    argstr = args(i)%it
    if ((argstr(1:2).eq.'-K').OR.(argstr(1:2).eq.'-k').OR.&
 &      (argstr(1:5).eq.'-KEY').OR.(argstr(1:5).eq.'-key').OR.&
 &      (argstr(1:5).eq.'-Key')) then
      keyfile = args(i+1)%it
      call strlims(keyfile,t1,t2)
      inquire (file=keyfile(t1:t2),exist=exists)
      if (exists.EQV..false.) then
        write (ilog,*) 'Fatal. Could not find or open requested keyf&
 &ile (got: ',keyfile(t1:t2),').'
        call fexit()
      else
        exit
      end if
    else if (argstr(1:1).eq.'-') then
      call strlims(argstr,t3,t4)
      write(ilog,'(3(a))') ' Fatal. Unrecognized argument/option (got ',argstr(t3:t4),').'
      write(ilog,*)
      write(ilog,'(a)') ' USAGE: CAMPARI requires a mandatory key-file with input options &
 &that specify the calculation to be attempted.'
      write(ilog,*)
      write(ilog,'(a)') ' EXAMPLE: ${PATH_TO_CAMPARI}/bin/${ARCH}/campari -k ${PATH_TO_CAMPARI}/examples/tutorial1/TEMPLATE.key'
      call fexit()
    end if
  end do
!
! this is just for braindead cases (vanilla key-file in local directory) ...
!
  if (exists.EQV..false.) then
    keyfile = 'campari.key'
    call strlims(keyfile,t1,t2)
    inquire (file=keyfile(t1:t2),exist=exists)
  end if
!
! now obtain and store keyfile in a set of lines (see keys.i) in two steps
!
  nkey = 0
! note this format spec. (#20) has to be manually synched with MAXKEYLEN
 20   format (a200)
  if (exists.EQV..true.) then
!
    iu = freeunit()
    open(unit=iu,file=keyfile(t1:t2),status='old')
    rewind(unit=iu)
!
    do while (.true.)
!
      read(iu,20,iostat=iomess) linec
      do j=1,len(linec)
        linea(j) = linec(j:j)
      end do
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then 
        write(ilog,*) 'Fatal. I/O-error with specified keyfile (got:&
 & ',keyfile(t1:t2),').'
        call fexit()
      end if
!
      call strlims_k(linea,t3,t4)
      if (t4.ge.t3) then
        nkey = nkey + 1
      end if
    end do
!
!   now allocate memory for keyfile, remember that too long lines will always be 
!   chopped to MAXKEYLEN
!
    allocate(key(nkey))
!
    rewind (unit=iu)
    nkey = 0
!
    do while (.true.)
      read(iu,20,iostat=iomess) linec
      do j=1,len(linec)
        linea(j) = linec(j:j)
      end do
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then 
        write(ilog,*) 'Fatal. I/O-error with specified keyfile (got:&
 & ',keyfile(t1:t2),').'
        call fexit()
      end if
!
      call strlims_k(linea,t3,t4)
      if (t4.ge.t3) then
        nkey = nkey + 1
        if ((t4-t3+1).ge.MAXKEYLEN) then
          allocate(key(nkey)%line(MAXKEYLEN))
          do j=1,MAXKEYLEN
            key(nkey)%line(j) = linea(j+t3-1)
          end do
        else
          allocate(key(nkey)%line(t4-t3+1))
          do j=1,(t4-t3+1)
            key(nkey)%line(j) = linea(j+t3-1)
          end do
        end if
      end if
    end do
!
    close (unit=iu)
!
  end if
!
end
!
!------------------------------------------------------------------------
