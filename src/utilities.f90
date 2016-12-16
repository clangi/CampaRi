subroutine Merge(A,Aix,NA,B,Bix,NB,C,Cix,NC,or) !(T,T2,NA,A(NA+1),ix(NA+1),NB,A,ix,N)

   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   real, intent(in out) :: A(NA) !T       ! B overlays C(NA+1:NC)
   real, intent(in)     :: B(NB) !A(NA+1)
   real, intent(in out) :: C(NC) !A
   integer, intent(in out) :: Aix(NA) !Tix      ! B overlays C(NA+1:NC)
   integer, intent(in)     :: Bix(NB) !ix(NA+1)
   integer, intent(in out) :: Cix(NC) !ix
   logical, intent(in) :: or

   integer :: I,J,K
   logical :: cond2

   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)

    if(or) then
      cond2 = (A(I) <= B(J))
    else
      cond2 = (A(I) >= B(J))
    end if
    if (cond2) then
      C(K) = A(I)
      Cix(K) = Aix(I)
      I = I+1
    else
      C(K) = B(J)
      Cix(K) = Bix(J)
      J = J+1
    endif
    K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      Cix(K) = Aix(I)
      I = I + 1
      K = K + 1
   enddo
   return
end subroutine merge

recursive subroutine MergeSort(A,ix,N,T,Tix,order_min)

   integer, intent(in) :: N !number of elements to sort
   real, dimension(N), intent(in out) :: A !vector to sort
   integer, dimension(N), intent(in out) :: ix !indexes of the vector to sort
   real, dimension((N+1)/2), intent (out) :: T !helper for recursion
   integer, dimension((N+1)/2), intent (out) :: Tix !helper for indexes recursion
   logical, intent(in) :: order_min !if true the order is growing in values

   integer :: NA,NB,V2
   real :: V
   logical :: cond



  if (N < 2) return
  if (N == 2) then
    !setting the order
    if(order_min) then
      cond =  A(1) > A(2)
    else
      cond = A(1) < A(2)
    end if
    if (cond) then
       V = A(1)
       A(1) = A(2)
       A(2) = V

       V2 = ix(1)
       ix(1) = ix(2)
       ix(2) = V2
    endif
    return
  endif
   NA=(N+1)/2
   NB=N-NA

   !recursing
   call MergeSort(A,ix,NA,T,Tix,order_min)
   call MergeSort(A(NA+1),ix(NA+1),NB,T,Tix,order_min)

   !setting the order
   if(order_min) then
     cond =  A(NA) > A(NA+1)
   else
     cond = A(NA) < A(NA+1)
   end if

   !merging
   if (cond) then
      T(1:NA)=A(1:NA)
      Tix(1:NA)=ix(1:NA)
      call Merge(T,Tix,NA,A(NA+1),ix(NA+1),NB,A,ix,N,order_min)
   endif

   return

end subroutine MergeSort

! program TestMergeSort
!
!    integer, parameter :: N = 8
!    real, dimension(N) :: A = (/ 1.1, 5.0, 2.1, 7.4, 3.8, 9.1, 4.1, 6.3 /)
!    integer, dimension(N) :: ix
!    real, dimension ((N+1)/2) :: T
!    integer, dimension ((N+1)/2) :: Tix
!
!    integer i
!
!    do i=1,N
!       ix(i) = i
!    end do
!
!    call MergeSort(A,ix,N,T,Tix)
!
!
!    write(ilog,*)'Sorted array :',A
!    write(ilog,'(A,/,10I3)')'Sorted indexes :', ix
!
! end program TestMergeSort


!-----------------------------------------------------------------------
! a function which provides an empty file handle
! it just loops up until it finds an open one (tests explicitly)
!
function freeunit()
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
    if ((freeunit.ne.5).AND.(freeunit.ne.6)) then
      inquire (unit=freeunit,opened=inuse)
    end if
  end do
!
end
!
!-------------------------------------------------------------------
!
! binary search for a target interval confinement in a sorted array of integers
! equivalence is allowed with left interval bound and match for series of identical values is rightmost
!
subroutine ibinary_search(arsz,arvec,keyv,ifound)
!
  implicit none
!
  integer arsz,ifound,k,ivsz
  integer keyv,arvec(*)
  logical notdone
!
  notdone = .true.
  if (keyv.lt.arvec(1)) then
    ifound = 0
    return
  else if (keyv.ge.arvec(arsz)) then
    ifound = arsz+1
    return
  end if
!
  k = 1
  ivsz = max(1,nint(dble(arsz)/2.0))
  do while (k.lt.arsz)
    if (arvec(min(arsz,k+ivsz)).gt.keyv) then
!     stay in lower half
    else
      k = min(arsz,k+ivsz)
    end if
    if (ivsz.eq.1) exit
    ivsz = max(1,nint(dble(ivsz)/2.0))
  end do
!
  ifound = k
!
end
!
!-----------------------------------------------------------------------
!
! the PRNG from the references below (ideally, multiple options
! to be offered here)
!
! literature references:
!
! P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)
!
! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
! Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
! University Press, 1992, Section 7-1
!
!
function random_or()
!
  use m_variables_gen
  implicit none
  !
  integer MAXKEYLEN,MAXKWLEN
  parameter (MAXKEYLEN=200)
  parameter (MAXKWLEN=25)
  !
  integer nkey
  logical have_legit_key
  !
  type t_charmat
    character, ALLOCATABLE:: line(:)
    logical legit(3)
  end type t_charmat
  type(t_charmat), ALLOCATABLE:: key(:)
  integer im1,ia1,iq1,ir1
  integer im2,ia2,iq2,ir2
  integer big,ntable
  integer imm1,ndiv,j
  real factor,random_or
  parameter (im1=2147483563)
  parameter (ia1=40014)
  parameter (iq1=53668)
  parameter (ir1=12211)
  parameter (im2=2147483399)
  parameter (ia2=40692)
  parameter (iq2=52774)
  parameter (ir2=3791)
  parameter (big=141803398)
  parameter (ntable=32)
  parameter (imm1=im1-1)
  parameter (ndiv=1+imm1/ntable)
  parameter (factor=1.0d0/im1)
  integer i,k
  integer(KIND=8) dts(8)
  character(12) rc(3)
  integer cont,iomess
  character(MAXKWLEN) seedstr
  character(MAXKEYLEN) kline,readfrom
!
!
! the default seed is set to a large number and then time-incremented
!
  if (firstcall.EQV..true.) then
    firstcall = .false.
    seed = big

    call Date_and_Time(rc(1),rc(2),rc(3),dts)
    dts(1) = mod(int(dts(1)),10)-1
    seed = seed + 32140800*dts(1) + 2678400*(dts(2)-1)
    seed = seed + 86400*(dts(3)-1) + 3600*dts(5)
    seed = seed + 60*dts(6) + dts(7) + dts(8)
!
!   provide further increment by process ID
    ! seed = abs(seed + handwrapped_getpid()) NOT IMPLEMENTED
!
!   get a user-specified seed
    if(rand_seed.gt.0) seed = rand_seed
!
    write(ilog,*)
    write(ilog,*) 'Initialized PRNG with Seed of ',seed
    write(ilog,*)
!
!   initialize shuffling table
!
    seed2 = seed
    do i=ntable+8,1,-1
      k = seed / iq1
      seed = ia1 * (seed-k*iq1) - k*ir1
      if (seed.lt.0)  seed = seed + im1
      if (i.le.ntable)  itable(i) = seed
    end do
    iy = itable(1)
  end if

!
! get a new random number value each call
!
  k = seed/iq1
  seed = ia1*(seed-k*iq1) - k*ir1
  if (seed.lt.0)  seed = seed + im1
  k = seed2/iq2
  seed2 = ia2*(seed2-k*iq2) - k*ir2
  if (seed2.lt.0)  seed2 = seed2 + im2
  i = 1 + iy/ndiv
  iy = itable(i) - seed2
  itable(i) = seed
  if (iy.lt.1)  iy = iy + imm1
  random_or = factor*iy
!
end

!-------------------------------------------------------------------------------
!
! Random generator based on the standard random_number gnu subroutine.
! This is a wrapper to have a function.
!
function random_st()
  use m_variables_gen
  implicit none
  integer,allocatable :: seed_a(:)
  integer n
  real random_st
  if(rand_seed.gt.0) seed = rand_seed
  if(seed.gt.0) then
    call random_seed(size = n)
    allocate(seed_a(n))
    seed_a(:) = seed
    call random_seed(put = seed_a)
  end if

  call random_number(random_st)
!
end

!-----------------------------------------------------------------------------
!
! fatal exit: call when encountering an internal problem
!
subroutine fexit()
!
!
  implicit none
!
  integer k
!
!
  k = 6
  write(k,*)
  write(k,*) '------------------------------------------------>'
  write(k,*)
  write(k,*) 'CAMPARI CRASHED UNEXPECTEDLY. PLEASE RECORD ANY '
  write(k,*) ' INFORMATION ABOUT THE PROBLEM PROVIDED ABOVE.'
  write(k,*)
  write(k,*) '<------------------------------------------------'
!
  ! stop 98
!
end
!
!-----------------------------------------------------------------------
