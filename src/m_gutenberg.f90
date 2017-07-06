module gutenberg
  logical mute !global verbose var
  ! functions for printing variables
contains
  ! string print
  subroutine spr(str)
    implicit none
    character(*) :: str
    if(.not.mute) call intpr(str,-1,0,0)
  end subroutine
  ! integer print
  subroutine ipr(int)
    implicit none
    integer :: int
    if(.not.mute) call intpr('',-1,int,1)
  end subroutine
  ! real print
  subroutine rpr(re)
    implicit none
    real :: re
    if(.not.mute) call realpr('',-1,re,1)
  end subroutine
  ! string + real print -> TODO: making it on the same line
  subroutine srpr(str, re)
    implicit none
    character(*) :: str
    real :: re
    if(.not.mute) call realpr(str,-1,re,1)
  end subroutine
  ! real vector print
  subroutine rvpr(re)
    implicit none
    real, intent(in) :: re(:)
    if(.not.mute) call realpr('',-1,re,size(re))
  end subroutine
  ! print a string like "gimaica: ",int
  subroutine sipr(str, int)
    implicit none
    character(*), intent(in) :: str
    integer, intent(in) :: int
    character(255) :: strint
    integer :: counted_digits
    ! let's count the digits
    call count_digits_int(int,counted_digits)
    ! let's make the integer a string
    call int2str(int,strint,counted_digits)
    ! print it
    if(.not.mute) call spr(trim(str)//" "//trim(strint))
  end subroutine
  ! skip line
  subroutine sl()
    implicit none
    if(.not.mute) call intpr(' ',-1,0,0)
  end subroutine
  ! counts the digits in int
  subroutine count_digits_int(int,dig)
    integer, intent(in) :: int
    integer, intent(in out) :: dig
    real :: to_floor
    integer :: to_flooi
    do dig = 1,300
      to_floor = to_flooi/10
      if(dig.eq.1) to_floor = int/10
      to_flooi = floor(to_floor)
      if(to_flooi.eq.0) exit
    end do
  end subroutine

  ! original struncture of cluster summary
  ! 67 format(i5,3x,i10,4x,g14.4,1x,i11,4x,i12)
  ! 68 format(i5,3x,i10,7x,a7,5x,i11,4x,i12)
  ! call spr('Level    #Clusters      TotalSnaps    Tot Children')
  ! write(ilog,68) 1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
  ! sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
  ! do i=2,c_nhier+1
  !   write(ilog,67) i,birchtree(i)%ncls,scrcts(i),sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
  !   sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren)
  ! end do

  ! SPECIAL prints: cluster summary
  subroutine clu_summary_linepr(int5, int10, int11, int12)
    implicit none
    integer, intent(in) :: int5, int10, int11, int12
    character(255) :: strint
    character(255) :: str
    integer :: counted_digits, l

    ! int5
    call count_digits_int(int5,counted_digits)
    call int2str(int5,strint,counted_digits)
    str = ''
    do l=1,(5-counted_digits) !adding white space in front
      str = str//' '
    end do
    str = str//trim(strint)//'   '

    ! int10
    call count_digits_int(int10,counted_digits)
    call int2str(int10,strint,counted_digits)
    do l=1,(10-counted_digits) !adding white space in front
      str = str//' '
    end do
    str = str//trim(strint)//'     '

    ! int11
    call count_digits_int(int11,counted_digits)
    call int2str(int11,strint,counted_digits)
    do l=1,(11-counted_digits) !adding white space in front
      str = str//' '
    end do
    str = str//trim(strint)//'    '

    ! int12
    call count_digits_int(int12,counted_digits)
    call int2str(int12,strint,counted_digits)
    do l=1,(11-counted_digits) !adding white space in front
      str = str//' '
    end do
    str = str//trim(strint)

    ! final call to string print
    if(.not.mute) call spr(str)
  end subroutine
end module
