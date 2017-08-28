subroutine find_mahalanobis_F(clu1e1, clu1e2, clu2, &
  &n_feature, n_elements, X_MA, must_be_positive, mute_it)

  ! imports
  use gutenberg
  use optimization_maha_metric
  use arpack_eig
  implicit none

  ! variables in input
  integer, intent(in) :: n_feature, n_elements
  real, dimension(n_feature, n_elements), intent(in) :: clu1e1, clu1e2
  real, dimension(n_feature, n_elements), intent(in) :: clu2
  real, dimension(n_feature, n_feature), intent(inout) :: X_MA
  logical, intent(in) :: mute_it
  logical, intent(in) :: must_be_positive ! the distance function must or not be positive

  ! stuff that will be in input (maybe)
  ! logical mute ! to comment when in R

  ! internal variables
  real, dimension(:,:), allocatable :: grad_obj_fu
  real epsilon
  integer :: K=10, J=10  ! max number of iterations
  integer :: ki, ji, i, kkk ! looping vars
  real :: rho(2) ! 2 is the now rho, 1 is the last iteration rho
  real :: rho_out_tmp, f_out_tmp ! tmp
  real :: chosen_rho, functio_out_pos_rho
  real :: of1, of2, of3, val1, val2 ! values for the last convergence condition
  real :: eig_i(1,3) ! eigenvalue
  real, allocatable :: eigv_i(:,:) ! eigenvector
  real, allocatable :: search_dir(:,:)
  real, allocatable :: alpha_dir(:,:)

  ! small garbage
  integer AllocateStatus, DeAllocateStatus
  integer AllocateStatus2, DeAllocateStatus2
  integer AllocateStatus3, DeAllocateStatus3
  integer AllocateStatus4, DeAllocateStatus4
  real :: time1, time2, time3, tot_time

  ! assignment to the module variables
  n_triplets = n_elements
  n_features = n_feature

  ! Total timing
  call cpu_time(tot_time)

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
  allocate(alpha_dir(n_features, n_features))
  allocate(eigv_i(1, n_features))

  i_j = clu1e1 - clu1e2
  i_k = clu1e1 - clu2

  if(all(i_j .eq. 0.0) .or. all(i_k .eq. 0.0)) then
    call fexit("Problem: you inserted identical triplets &
    &(internal or external snapshots are identical).")
  end if

  call cpu_time(time1)

  ! variables initialization
  K = 10 ! big loop
  J = 10 ! small loop

  ! X_0
  X_k = 0
  ! do k=1,n_features
  !   X_k(k,k) = 1.0/n_features ! not correct
  ! end do
  X_k(1,:) = 1 ! I can use only this for this calculations (must be rank=1)
  ! in order to calculate the eigenvalues correctly I need a symmetric matrix
  X_k(:,1) = 1 ! this is not hampering the rank discussion as it is still rank=1

  ! sensibility
  epsilon = 0.00001
  print_rho_components = .false.
  detailed_print_rho = .false.
  mute = mute_it

  ! K LOOP ---------------------------------------------------------------------
  do ki=1,K

    ! print k
    call sipr('k (over 10): ', ki)

    !compute rho by solving the argmax problem
    if(ki .gt. 1) rho(1) = rho(2)
    kkk = 0
    do i=1, n_triplets
      call rho_f(i, rho_out_tmp)
      call obj_f(rho_out_tmp, f_out_tmp)
      ! select the positive rho which maximize obj_f
      if(rho_out_tmp .gt. 0.0) then
        if(kkk .eq. 0 .or. f_out_tmp .gt. functio_out_pos_rho) then
          functio_out_pos_rho = f_out_tmp
          chosen_rho = rho_out_tmp
          kkk = kkk + 1
        end if
      end if
    end do

    ! some check on the found rho
    if(chosen_rho .le. 0) then
      call fexit("Found non strictly positive rho")
    end if
    if(kkk .eq. 0) then
      call fexit("No positive rho found. The introduced triplets &
      & are not suited for this learning.")
    end if

    ! finding max value of the objfunction and its rho
    rho(2) = chosen_rho
    call cpu_time(time2)
    call srpr("Time to find Rho[s]: ", time2-time1)

    ! printing the results
    ! print *, "-----------------"
    ! print *, "Output of obj_f: ", functio_out_pos_rho
    ! print *, "-----------------"
    ! print *, "rho: ", rho(2)
    ! print *, "-----------------"

    ! Compute X as argmax of the obj_function
    X_tmp = X_k ! X_tmp is the X BEFORE optimization (aka X k - 1)

    ! J LOOP -------------------------------------------------------------------
    do ji=1,J
      ! print j
      call sipr('j (over 10): ', ji)

      ! Calculate the gradient matrix (it depends on the matrix X_k only)
      call grad_obj_f(rho(2), grad_obj_fu)

      ! Find first eigenvalue (and first eigenvector)
      A = grad_obj_fu
      ! n_features = Dimension of the eigen-problem
      ! 1 = number of eigenvalues to calculate
      ! 'LM' = largest magnitude
      call spr('=------------------------')
      call find_eigens(eig_i, eigv_i, n_features, 1, 'LM')
      call srpr('eigenvalue ', eig_i(1,1))
      call spr('=------------------------')
      ! print *, 'eigenvector ', eigv_i(1,:)
      ! print *, '=------------------------'
      if(eig_i(1,1) .lt. epsilon) then ! without abs is much fuster
        call spr('---------------------------------------&
        &-------------------- Egenvalues converged')
        exit ! eigenvalue < epsilon  ==  converged
      end if

      ! compute search direction p_i
      search_dir = matmul(transpose(eigv_i(1:1,:)),eigv_i(1:1,:)) - X_k
      ! print *, '=------------------------'
      ! do kkk=1,n_features
      !   print *, "Direction: ", search_dir(kkk,:)
      ! end do
      ! print *, '=------------------------'

      ! search alpha using backtracking line search
      call line_search_alpha(alpha_dir, grad_obj_fu, rho(2), search_dir, &
      & must_be_positive)
      ! line search alpha_dir convergence
      ! do kkk=1,n_features
      !   print *, 'alpha*search_dir: ', alpha_dir(kkk,:)
      ! end do
      if(all(alpha_dir .lt. 0.01)) then
        call spr('-------------------------------&
        &-------------------------- Alpha made convergence')
        exit
      end if

      X_k = X_k + alpha_dir ! alpha must be found with linear searching
    end do


    ! setting it for the next iteration
    call cpu_time(time1)
    call srpr("Time to find the new X_k[s]: ", time1-time2)
    X = X_k ! this is the new X (X_k) (saved)

    ! condition to break J loop
    if(ki .gt. 1) then
      call obj_f(rho(2), of1) ! it is using already the new value (X)
      X_k = X_tmp ! this is from before!
      call obj_f(rho(2), of2) ! it is using X_k-1 (the one from before)
      call obj_f(rho(1), of3)
      val1 = abs(of1 - of2)
      val2 = abs(of2 - of3)

      if(val1 .lt. epsilon .and. val2 .lt. epsilon) then
        call spr('-----------------------------------&
        &---------------------- Distance mat converged')
        exit
      end if
      X_k = X ! it did not converge -> continue iterations
    end if
  end do

  ! print summary
  call cpu_time(time3)
  call sl()
  call srpr("Time needed to find the Mahalanobis distance[s]:", time3-tot_time)
  call spr('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~&
  &~~~~~~~~~~~~~')
  call spr("                                      ------X-----")
  call sl()
  kkk = n_features
  do k=1,n_features
    if(n_features .gt. kkk) exit
    kkk = n_features
    if(n_features .gt. 10) kkk = 10
    call rvpr(X(k,1:kkk))
  end do
  call sl()
  call spr('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~&
  &~~~~~~~~~~~~~')
  call spr("Resulting summed learned distance for Internal/External elements: ")
  print_rho_components = .true.
  detailed_print_rho = .true. ! the not detailed one was WRONG
  X_k = X
  call rho_f(1,rho(1))
  call spr(' ----------------------------------- ')
  call spr("Resulting summed Euclidean distance for &
  &Internal/External elements: ")
  X_k = 0
  do k=1,n_features
    X_k(k,k) = 1.0/n_features
  end do
  call rho_f(1,rho(1))
  call spr(' ----------------------------------- ')

  ! main return of the function
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
