program to_del
  ! use optimization_maha_metric
  implicit none
  integer, parameter :: n_elements = 1000
  integer, parameter :: n_feature = 5
  logical :: mute_it
  real :: clu1e1(n_feature,n_elements), clu1e2(n_feature,n_elements)
  real :: clu2(n_feature,n_elements)
  real :: X_MA(n_feature,n_feature)
  real :: shift = 1.2

  ! ! testing the where() function
  ! real :: testm(2,2)
  ! real :: testm2(2,2)
  ! testm(1,:) = (/1,2/)
  ! testm(2,:) = (/-1,2/)
  ! testm2(1,:) = (/5,6/)
  ! testm2(2,:) = (/234,21/)
  ! where(testm .gt. 0.0) testm2 = testm2*0.11111
  ! print *, testm2
  ! ! end of testing where

  ! ! Stuff for testing the matrix gen from vector multiplication
  ! real :: testv1(2,1)
  ! real :: testv2(2,1)
  ! real :: onesv(2,2)
  ! integer :: uia
  ! testv1(:,1) = (/1,2/)
  ! testv2(:,1) = (/3,4/)
  ! onesv = matmul(testv1, transpose(testv2))
  ! do uia=1,2
  !   print *, onesv(uia,:)
  ! end do
  ! mute_it = .false.
  ! ! end of the testing for outer product


  X_MA = -1 ! Just to be sure that is actually changing it
  mute_it = .false.

  call random_number(clu1e1)
  call random_number(clu1e2)
  call random_number(clu2)
  ! print *, clu2

  clu1e1(2:3,:) = clu1e1(2:3,:) + shift
  clu1e2(2:3,:) = clu1e2(2:3,:) + shift
  clu2(2:3,:) = clu2(2:3,:) - shift

  ! clu1e1(3,:) = clu1e1(3,:) + shift
  ! clu1e2(3,:) = clu1e2(3,:) + shift
  ! clu2(3,:) = clu2(3,:) - shift
  !
  ! clu1e1(:,:) = clu1e1(:,:) + shift
  ! clu1e2(:,:) = clu1e2(:,:) + shift
  ! clu2(:,:) = clu2(:,:) - shift

  ! normalization
  clu1e1 = clu1e1 / maxval(abs(clu1e1))
  clu1e2 = clu1e2 / maxval(abs(clu1e2))
  clu2 = clu2 / maxval(abs(clu2))
  ! print *, clu2

  ! The matrix multiplication mistery
  ! print *, clu1e1(1:2,1)
  ! print *, clu1e1(1:2,2)
  ! print *,
  ! print *, clu1e2(1:2,1)
  ! print *,
  ! clu1e1(1:2,1) = matmul(clu1e1(1:2,1:2), clu1e2(1:2,1))
  ! print *, clu1e1(1:2,1)
  ! print *, clu1e1(1,1)*clu1e2(1,1)
  ! print *, clu1e1(2,1)*clu1e2(2,1)
  ! print *,
  ! print *, clu2(1:2,1)
  ! print *,
  ! print *, dot_product(clu2(1:2,1), clu1e1(1:2,1))
  ! Mystery still not solved


  call find_mahalanobis_F(clu1e1, clu1e2, clu2, &
  & n_feature, n_elements, X_MA, .true., mute_it)


end program
