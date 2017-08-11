subroutine find_mahalanobis_F(clu1e1, clu1e2, clu2, n_feature, n_elements, X_MA)
  use optimization_maha_metric
  use arpack_eig
  implicit none
  integer, intent(in) :: n_feature, n_elements
  real, dimension(n_feature, n_elements), intent(in) :: clu1e1, clu1e2
  real, dimension(n_feature, n_elements), intent(in) :: clu2
  real, dimension(n_feature, n_feature), intent(inout) :: X_MA
  real, dimension(:,:), allocatable :: grad_obj_fu
  real epsilon
  integer :: K=10, J=10  ! max number of iterations
  integer :: ki, ji, i
  ! real :: rho(K)
  real :: rho(1)
  real :: rho_out_tmp, f_out_tmp
  real :: positive_rho(n_elements)
  real :: functio_out_pos_rho(n_elements)
  integer :: kkk
  real :: of1, of2, of3, val1, val2
  real :: alpha ! to be found by line search
  real :: eig_i(1,3) ! eigenvalue
  real, allocatable :: eigv_i(:,:)
  real, allocatable :: search_dir(:,:)
  integer :: index_of_max
  integer AllocateStatus, DeAllocateStatus
  integer AllocateStatus2, DeAllocateStatus2
  integer AllocateStatus3, DeAllocateStatus3
  integer AllocateStatus4, DeAllocateStatus4
  real :: time1, time2, time3

  ! assignment to the module variables
  n_triplets = n_elements
  n_features = n_feature

  allocate(X_k(n_features, n_features), STAT = AllocateStatus)
  if(AllocateStatus /= 0) stop "*** Not enough memory ***"
  allocate(grad_obj_fu(n_features, n_features), STAT = AllocateStatus2)
  if(AllocateStatus /= 0) stop "*** Not enough memory ***"
  allocate(X(n_features, n_features), STAT = AllocateStatus3)
  if(AllocateStatus /= 0) stop "*** Not enough memory ***"
  allocate(X_tmp(n_features, n_features), STAT = AllocateStatus4)
  if(AllocateStatus /= 0) stop "*** Not enough memory ***"
  allocate(i_j(n_features, n_triplets))
  allocate(i_k(n_features, n_triplets))
  allocate(A(n_features, n_features))
  allocate(search_dir(n_features, n_features))
  allocate(eigv_i(1, n_features))

  i_j = clu1e1 - clu1e2
  i_k = clu1e1 - clu2

  if(all(i_j .eq. 0) .or. all(i_k .eq. 0)) then
    print *, "Problem: you inserted identical triplets (internal or external snapshots are identical)."
    stop
  end if

  ! variables initialization
  K = 10
  J = 10

  X_k = 0
  do k=1,n_features
    X_k(k,k) = 1.0/n_features
  end do
  ! X_k(1,:) = 1
  epsilon = 0.0001
  positive_rho = 0
  functio_out_pos_rho = 0
  print_rho_components = .false.
  call cpu_time(time1)
  do ki=1,K
    ! print *, "k = ", ki
    !compute rho by solving the argmax problem
    kkk = 1
    do i=1, n_triplets
      call rho_f(i, rho_out_tmp)
      call obj_f(rho_out_tmp, f_out_tmp)
      if(rho_out_tmp .gt. 0.001) then
        positive_rho(kkk) = rho_out_tmp
        functio_out_pos_rho(kkk) = f_out_tmp
        kkk = kkk + 1
      end if
    end do
    kkk = kkk - 1
    if(any(positive_rho(1:kkk) .le. 0)) then
      print *, "Found non strictly positive rho"
      stop
    end if
    if(kkk .eq. 0) then
      print *, "No positive rho found. The introduced triplets are not suited for this learning."
      stop
    end if
    index_of_max = maxloc(functio_out_pos_rho(1:kkk), dim = 1)
    rho(k) = positive_rho(index_of_max)
    ! print *, "selected: ", index_of_max
    ! print *, "-----------------"
    ! print *, functio_out_pos_rho
    ! print *, kkk
    ! print *, "-----------------"
    ! print *, positive_rho
    ! print *, "-----------------"
    ! compute X as argmax of the obj_function
    X_tmp = X_k ! X_tmp is the X BEFORE optimization (aka X k-1)
    ! print *, "X_k: ", X_k
    ! print *, ' '
    call cpu_time(time2)
    print *, "Time to find Rho: ", time2-time1
    do ji=1,J
      ! find the grad
      call grad_obj_f(rho(k), grad_obj_fu)
      ! compute eig_i max eigenvalue with eigv_i its vector from grad_obj
      A = grad_obj_fu
      ! print *, "A: ", A
      nx = n_features
      nev = 1
      ncv = 2*nev + 1
      call find_eigens(eig_i, eigv_i)
      ! print *, '=------------------------'
      ! print *, 'eigenvalue', eig_i(1,1)
      ! print *, '=------------------------'
      search_dir = matmul(transpose(eigv_i(1:1,:)),eigv_i(1:1,:))
      ! if(abs(eig_i(1,1)) .lt. epsilon .or. abs(dot_product(eigv_i(1,:),eigv_i(1,:))) .lt. epsilon) then
      if(eig_i(1,1) .lt. epsilon) then ! ATTENTION! originally without abs
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!egenvalues converged!!!!!!!!!!!!!!!!!!!!!!!'
        exit ! converged
      end if
      ! print *, '=------------------------'
      ! print *, "Direction: ", search_dir
      ! print *, '=------------------------'
      ! search alpha using backtracking line search
      call line_search_alpha(alpha, grad_obj_fu, rho(k), search_dir)
      ! alpha = 1.5
      X_k = X_k + alpha*(search_dir - X_k) ! alpha must be found with linear searching
      ! print *, "X_ksad: ", X_k
    end do
    call cpu_time(time1)
    print *, "Time to find the new X_k: ", time1-time2
    X = X_k ! this is the new X (X_k)
    ! condition to break it
    if(k.gt.1) then
      call obj_f(rho(k), of1) ! it is using already the new value (X)
      X_k = X_tmp ! this is from before!
      call obj_f(rho(k), of2) ! it is using X_k-1 (the one from before)
      call obj_f(rho(k-1), of3)
      val1 = abs(of1 - of2)
      val2 = abs(of2 - of3)
      if(val1 .lt. epsilon .and. val2 .lt. epsilon) then
        print *, ' !!!!!!!!!!!!!!!!!!!!!!!!X_k converged to X!!!!!!!!!!!!!!!!!!!!!!!'
        exit ! converged
      end if
      X_k = X ! it did not converge -> continue iterations
    end if
  end do
  ! print *, rho
  ! print *, "X_k: ", X_k ! it must be the one from before
  print *, " "
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print *, "                        ------X-----"
  print *, ""
  do k=1,n_features
    print *, X(k,:)
  end do
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  ! print *,
  print_rho_components = .true.
  X_k = X
  call rho_f(1,rho(1))
  print *, ' ----------------------------------- '
  X_k = 0
  do k=1,n_features
    X_k(k,k) = 1
  end do
  call rho_f(1,rho(1))
  print *, ' ----------------------------------- '
  X_k = 0
  X_k(3,3) = 1
  X_k(2,2) = 1
  call rho_f(1,rho(1))
  X_MA = X
  ! deallocation
  deallocate(X_k, STAT = DeAllocateStatus)
  deallocate(grad_obj_fu, STAT = DeAllocateStatus2)
  deallocate(X, STAT = DeAllocateStatus3)
  deallocate(X_tmp, STAT = DeAllocateStatus4)
  deallocate(i_j)
  deallocate(i_k)
  deallocate(A)
end
