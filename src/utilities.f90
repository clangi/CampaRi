subroutine Merge(A,Aix,NA,B,Bix,NB,C,Cix,NC) !(T,T2,NA,A(NA+1),ix(NA+1),NB,A,ix,N)

   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   real, intent(in out) :: A(NA) !T       ! B overlays C(NA+1:NC)
   real, intent(in)     :: B(NB) !A(NA+1)
   real, intent(in out) :: C(NC) !A
   integer, intent(in out) :: Aix(NA) !Tix      ! B overlays C(NA+1:NC)
   integer, intent(in)     :: Bix(NB) !ix(NA+1)
   integer, intent(in out) :: Cix(NC) !ix

   integer :: I,J,K

   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
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

recursive subroutine MergeSort(A,ix,N,T,Tix)

   integer, intent(in) :: N
   real, dimension(N), intent(in out) :: A
   integer, dimension(N), intent(in out) :: ix
   real, dimension((N+1)/2), intent (out) :: T
   integer, dimension((N+1)/2), intent (out) :: Tix

   integer :: NA,NB,V2
   real :: V

   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
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

   call MergeSort(A,ix,NA,T,Tix)
   call MergeSort(A(NA+1),ix(NA+1),NB,T,Tix)

   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      Tix(1:NA)=ix(1:NA)
      call Merge(T,Tix,NA,A(NA+1),ix(NA+1),NB,A,ix,N)
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
