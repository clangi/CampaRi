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
!
!             ##################################
!             #                                #
!             #  A SET OF SUBROUTINES WHICH    #
!             #  DEAL WITH PROGRAM FLOW AND    #
!             #  FILE INTERFACES               #
!             #                                #
!             ##################################
!
!
!-----------------------------------------------------------------------------
!
! fatal exit: call when encountering an internal problem
!
subroutine fexit()
!
#ifdef ENABLE_MPI
  use mpi
#endif
  use iounit
!
  implicit none
!
  integer k
  character(500) intfile
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  character(500) tnum
  integer lext
#endif
!
#ifdef ENABLE_MPI
  integer ierr
!
  call makelogio(2)
  call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#endif
!
#ifdef ENABLE_THREADS
  if (OMP_IN_PARALLEL().EQV..true.) then
    lext = int(log10(1.0*OMP_GET_NUM_THREADS()) + 1)
    call int2str(OMP_GET_THREAD_NUM,tnum,lext)
    intfile = 'CRASH_T'//tnum(1:lext)//' '
  else
    intfile = 'CRASH '
  end if
  if (OMP_IN_PARALLEL().EQV..true.) then
!$OMP CRITICAL(CRASHDUMP)
    call FMCSC_dump(intfile)
!$OMP END CRITICAL(CRASHDUMP)
  else
    call FMCSC_dump(intfile)
  end if
#else
  intfile = 'CRASH '
  call FMCSC_dump(intfile)
#endif
!
  k = std_f_out
  write(k,*)
  write(k,*) '------------------------------------------------>'
  write(k,*)
  write(k,*) 'CAMPARI CRASHED UNEXPECTEDLY. PLEASE RECORD ANY '
  write(k,*) ' INFORMATION ABOUT THE PROBLEM PROVIDED ABOVE.'
  write(k,*)
  write(k,*) '<------------------------------------------------'
!
  stop 98
!
end
!
!-----------------------------------------------------------------------
!
! a trivial fxn to obtain command line arguments
!
subroutine grab_args()
!
  use commline
  use iounit
!
  implicit none
!
  integer argstat,i,k
!
! the FORTRAN intrinsics can handle this well in modern Fortran
! note the zero-index argument is the program call itself and is
! not counted in COMMAND_ARGUMENT_COUNT
!
  nargs = COMMAND_ARGUMENT_COUNT()
!
  if (nargs.eq.0) then
    k = std_f_out
    write(k,*) 'USAGE: CAMPARI requires a mandatory key-file with input options &
 &that specify the calculation to be attempted.'
    write(k,*)
    write(k,*) 'EXAMPLE: ${PATH_TO_CAMPARI}/bin/${ARCH}/campari -k ${PATH_TO_CAMPARI}/examples/tutorial1/TEMPLATE.key'
    call fexit()
  end if
!
  allocate(args(nargs))     
!
  do i=1,nargs
    call GET_COMMAND_ARGUMENT(i,args(i)%it,args(i)%leng,argstat)
  end do
     
!
end
!
!-----------------------------------------------------------------------
!
! a function which provides an empty file handle
! it just loops up until it finds an open one (tests explicitly)
! 
! this routine is not thread-safe and should always only be called by the master thread
!
function freeunit()
!
  use iounit
  use mcsums, ONLY: itrajhlp,itraj
!
  implicit none
!
  integer freeunit
  logical inuse
!
! try each logical unit until an unopened one is found
!
  freeunit = 0
  inuse = .true.
  do while (inuse)
    freeunit = freeunit + 1
    if ((freeunit.ne.std_f_in).AND.(freeunit.ne.std_f_out).AND.(freeunit.ne.itrajhlp).AND.(freeunit.ne.itraj)) then
      inquire (unit=freeunit,opened=inuse)
    end if
!
    if (freeunit.eq.huge(freeunit)) then
      write(ilog,*) 'Fatal. No free I/O unit could be found. This mo&
 &st likely reports on a more severe underlying problem with the fil&
 &e system.'
      call fexit()
    end if
  end do
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! the same for a number of different unopened units at the same time
!
! also not thread-safe (master only)
!
subroutine get_freeunits(howmany,funits)
!
  use iounit
  use mcsums, ONLY: itrajhlp,itraj
!
  implicit none
!
  integer, INTENT(IN):: howmany
  integer, INTENT(OUT):: funits(howmany)
!
  integer i,tryu
  logical inuse
!
! try each logical unit until an unopened one is found
!
  if (howmany.le.0) then
    write(ilog,*) 'Fatal. Called get_freeunits(...) with a nonsensical argument for the requested number. This is a bug.'
    call fexit()
  end if
  i = 1
  funits(:) = 0
  tryu = 0
  do while (i.le.howmany)
    tryu = tryu + 1
    if ((tryu.ne.std_f_in).AND.(tryu.ne.std_f_out).AND.(tryu.ne.itrajhlp).AND.(tryu.ne.itraj)) then
      inquire (unit=tryu,opened=inuse)
    end if
    if (inuse.EQV..false.) then
      funits(i) = tryu
      i = i + 1
    end if
!
    if (tryu.eq.huge(tryu)) then
      write(ilog,*) 'Fatal. No free I/O unit could be found in get_freeunits(...). This mo&
 &st likely reports on a more severe underlying problem with the fil&
 &e system or a bug.'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this fxn saves about one line of code and checks whether a given
! unit is associated with an open file
!
function is_filehandle_open(fh)
!
  implicit none
!  
  integer fh
  logical is_filehandle_open
!
  inquire (unit=fh,opened=is_filehandle_open)
!
  return
!
end function is_filehandle_open
!
!-----------------------------------------------------------------------
!
! this routine is a shortcut for deleting a file, and reopening it 
! to the handle provided
! it is important that the passed string is properly stripped of blanks on the calling side
!
subroutine delete_then_openfile(fh,fpath)
!  
  implicit none
!  
  integer fh,freeunit
  character(len=*) fpath
  logical exists
!
  inquire(file=fpath,exist=exists)
  if (exists.EQV..true.) then
    fh = freeunit()
    open (unit=fh,file=fpath,status='old')
    close(unit=fh,status='delete')
  end if
!
  fh = freeunit()
  open (unit=fh,file=fpath,status='new') 
!
end subroutine delete_then_openfile
!
!-----------------------------------------------------------------------
!
! this routine is a shortcut for opening an existing file 
! it throws a warning upon failure to find such a file, but does not
! open a new file
! it is important that the passed string is properly stripped of blanks on the calling side
!
function openfile(fh,fpath)
!
  implicit none
!
  integer fh, freeunit
  character(len=*) fpath
  logical exists,openfile
!
  inquire(file=fpath,exist=exists)
  if (exists.EQV..true.) then
    fh = freeunit()
    open (unit=fh,file=fpath,status='old')
    openfile = .true.
  else
    write(*,*) 'WARNING: Cannot open input file (',fpath,') in openfile(...).'
    openfile = .false.
  end if
!
  return
!  
end function openfile
!
!-----------------------------------------------------------------------
!
! this is a redundant wrapper that saves zero lines of code
!
subroutine close_filehandle(fh)
!    
  implicit none
!
  integer fh
  logical is_filehandle_open
!
  if(is_filehandle_open(fh))  close(unit=fh)
!
end subroutine close_filehandle
!
!-----------------------------------------------------------------------
!
! further functions may be added
!
!   To be added. Common file manipulation support is lacking in fortran
!   subroutine copy_file(ffrom, fto)
!   
!     implicit none
! 
! !     character(len=*) :: fpath
! !     integer isrc, idst
! ! 
! !     OPEN(UNIT=ISRC, FILE='', ACCESS='DIRECT', STATUS='OLD', ACTION='READ', IOSTAT=IERR, RECL=1)
! ! OPEN(UNIT=IDST, FILE='', ACCESS='DIRECT', STATUS='REPLACE', ACTION='WRITE', IOSTATE=IERR, RE)
! ! IREC = 1
! ! DO
! !   READ(UNIT=ISRC, REC=IREC, IOSTAT=IERR) CHAR
! !   IF (IERR.NE.0) EXIT
! !   WRITE(UNIT=IDST, REC=I) CHAR
! !   IREC = IREC + 1
! ! END DO
! ! 
! ! 
!   end subroutine copy_file
!
