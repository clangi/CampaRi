
function get_unused_unit()

  use iounit

  implicit none
!
  integer get_unused_unit
  logical inuse
!
! try each logical unit until an unopened one is found
!
  get_unused_unit = 0
  inuse = .true.
  do while (inuse)
    get_unused_unit = get_unused_unit + 1
    if ((get_unused_unit.ne.std_f_in).AND.(get_unused_unit.ne.std_f_out)) then
      inquire (unit=get_unused_unit,opened=inuse)
    end if
!
    if (get_unused_unit.eq.huge(get_unused_unit)) then
      write(*,*) 'Fatal. No free I/O unit could be found. This mo&
&st likely reports on a more severe underlying problem with the fil&
&e system.'
      call fexit()
    end if
  end do
end function get_unused_unit

function is_filehandle_open(fh)

  implicit none
  
  integer fh
  logical is_filehandle_open

  inquire (unit=fh,opened=is_filehandle_open)
end function is_filehandle_open

!return a valid filehandle for the path fpath
!if this file exists, delete the existing file
subroutine delete_then_openfile(fh,fpath)
  
  implicit none
  
  integer fh, get_unused_unit
  character(len=*) :: fpath
  logical exists

  inquire(file=fpath,exist=exists)
  if (exists.EQV..true.) then
    fh = get_unused_unit()
    open (unit=fh,file=fpath,status='old')
    close(unit=fh,status='delete')
  end if

  fh = get_unused_unit()
  open (unit=fh,file=fpath,status='new') 
  
end subroutine delete_then_openfile

!return a valid filehandle for the path fpath
!file must exist on open, otherwise error
function openfile(fh,fpath)

  implicit none

  integer fh, get_unused_unit
  character(len=*) :: fpath
  logical exists, openfile

  inquire(file=fpath,exist=exists)
  if (exists.EQV..true.) then
    fh = get_unused_unit()
    open (unit=fh,file=fpath,status='old')
    openfile = .true.
  else
    write(*,*) 'WARNING: Cannot open input file (',fpath,') in openfile()'
    openfile = .false.
  end if
  
end function openfile

subroutine close_filehandle(fh)
    
   implicit none

   integer fh
   logical is_filehandle_open

    if(is_filehandle_open(fh))  close(unit=fh)

end subroutine close_filehandle

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
! 
! 
!   end subroutine copy_file


