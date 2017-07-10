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
  ! integer vector print
  subroutine ivpr(intr)
    implicit none
    integer :: intr(:)
    if(.not.mute) call intpr('',-1,intr,size(intr))
  end subroutine
  ! real print
  subroutine rpr(re)
    implicit none
    real, intent(in) :: re
    real :: re_help
    character(255) :: strR, strL, str_fin
    integer int, n_digits, counted_digits
    n_digits = 4
    int = floor(re)
    call count_digits_int(int, counted_digits)
    call int2str(int, strR, counted_digits)
    re_help = re - int*1.0
    re_help = re_help*(10**n_digits)
    int = floor(re_help)
    call count_digits_int(int, counted_digits)
    call int2str(int, strL, counted_digits)
    str_fin = trim(strR)//'.'//trim(strL)
    call spr(trim(str_fin))
  end subroutine
  ! string + real print
  subroutine srpr(str, re)
    implicit none
    character(*), intent(in) :: str
    real, intent(in) :: re
    character(255) :: str_h, str_fin
    call re2str(re, 4, str_h)
    str_fin = trim(str)//' '//trim(str_h)
    call spr(trim(str_fin))
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
  subroutine re2str(re, n_digits, str)
    character(*), intent(in) :: str
    real, intent(in) :: re
    real :: re_help
    character(255) :: strR, strL
    integer int, n_digits, counted_digits
    ! thi routine will split the real in 2 integers (per side of .)
    int = floor(re)
    call count_digits_int(int, counted_digits)
    call int2str(int, strR, counted_digits)
    re_help = re - int*1.0
    re_help = re_help*(10**n_digits)
    int = floor(re_help)
    call count_digits_int(int, counted_digits)
    call int2str(int, strL, counted_digits)
    str = trim(strR)//'.'//trim(strL)
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
    character(255) :: str_tmp
    integer :: counted_digits
    str = ''
    ! int5
    call count_digits_int(int5,counted_digits)
    call int2str(int5,strint,counted_digits)
    str_tmp = trim(strint)
    call add_spaces_in_front(str_tmp, 5-counted_digits)
    str = trim(trim(str)//trim(str_tmp))
    ! int10
    call count_digits_int(int10,counted_digits)
    call int2str(int10,strint,counted_digits)
    str_tmp = trim(strint)
    call add_spaces_in_front(str_tmp, 13-counted_digits)
    str = trim(trim(str)//trim(str_tmp))
    ! real7 + 7spaces
    ! int11
    call count_digits_int(int11,counted_digits)
    call int2str(int11,strint,counted_digits)
    str_tmp = trim(strint)
    call add_spaces_in_front(str_tmp, 16-counted_digits)
    str = trim(trim(str)//trim(str_tmp))
    ! int12
    call count_digits_int(int12,counted_digits)
    call int2str(int12,strint,counted_digits)
    str_tmp = trim(strint)
    call add_spaces_in_front(str_tmp, 16-counted_digits)
    str = trim(trim(str)//trim(str_tmp))
    ! final call to string print
    call spr(trim(str))
  end subroutine
  ! add space in front of a string
  subroutine add_spaces_in_front(str, n_spaces)
    implicit none
    character(*) :: str
    integer :: l,n_spaces
    do l=1,n_spaces !adding white space in front
      str = ' '//str
    end do
    str = trim(str)
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
  !---------------------------------------------------------------------------------------
  !
  subroutine int2str(ii,string,strsz)
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
    call spr('Warning. Possibly bad result from int2str(...) (inte&
   &ger overflow).')
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
      call spr('Warning. Possibly bad result from int2str(...) (digi&
   &t overflow).')
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
end module
