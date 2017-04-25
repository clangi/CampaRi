! ****** NOTE: This has been incorporated into a modern Fortran library ******
! ****** See https://github.com/wesbarnett/libgmxfort ******
! ****** 14 Dec 2016 *****

! XDR Fortran Interface Example Program
! 2014 (c) James W. Barnett <jbarnet4@tulane.edu>
! https://github.com/wesbarnett/

program read_xtc_prog

    use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
    use xtc

    implicit none
    character (len=1024) :: filename
    real, allocatable :: pos(:,:)
    integer :: NATOMS, STEP, STAT
    real :: box(3,3), prec, time, box_trans(3,3)
    type(C_PTR) :: xd_c
    type(xdrfile), pointer :: xd
    logical :: ex

    ! Set the file name for C.
    filename = "traj.xtc"//C_NULL_CHAR

    inquire(file=trim(filename),exist=ex)

    if (ex .eqv. .false.) then
        write(0,*)
        write(0,'(a)') " Error: "//trim(filename)//" does not exist."
        write(0,*)
        stop
    end if

    STAT = read_xtc_natoms(filename,NATOMS)
    allocate(pos(3,NATOMS))

    ! Open the file for reading. Convert C pointer to Fortran pointer.
    xd_c = xdrfile_open(filename,"r")
    call c_f_pointer(xd_c,xd)

    STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)

    do while ( STAT == 0 )

        ! C is row-major, whereas Fortran is column major. Hence the following.
        box = transpose(box_trans)

        ! Just an example to show what was read in
        write(*,'(a,f12.6,a,i0)') " Time (ps): ", time, "  Step: ", STEP
        write(*,'(a,f12.6,a,i0)') " Precision: ", prec, "  No. Atoms: ", NATOMS
        write(*,'(3f9.3)') pos

        ! This is the same order as found in the GRO format fyi
        write(*,'(11f9.5)') box(1,1), box(2,2), box(3,3), &
                            box(1,2), box(1,3), & 
                            box(2,1), box(2,3), &
                            box(3,1), box(3,2) 

        STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)

    end do

    STAT = xdrfile_close(xd)
    deallocate(pos)

end program read_xtc_prog
