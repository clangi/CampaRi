module optimization_maha_metric
  ! main variables
  logical :: print_rho_components, detailed_print_rho
  integer n_features
  integer n_triplets
  real, allocatable :: i_j(:,:), i_k(:,:) ! i-j internal vectors, i-k external
  real, dimension(:,:), allocatable :: X_k, X, X_tmp ! current and previous maha matrix

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
! rho function: it is the function that has to be maximised and it is composed
! by the inter-cluster distance MINUS the intra-cluster distance
!
  subroutine rho_f(which, rho_out)
    use gutenberg
    implicit none
    integer, intent(in) :: which
    real, intent(out) :: rho_out
    real :: help_v(n_features, 1)
    real :: help_v2(n_features, 1)
    real :: ave_distance(n_triplets)
    integer iter, iter2

    ! help_v3(:,1) = i_k(:,which) ! deprecated to i_k(:,which:which)
    help_v = matmul(X_k, i_k(:, which:which))
    help_v2 = matmul(X_k, i_j(:,which:which))
    rho_out = dot_product(i_k(:, which), help_v(:,1)) - dot_product(i_j(:, which), help_v2(:,1))

    ! I want to print the components (sum of the distances)
    ! ------------------------------ I_K (external clusters)
    if (print_rho_components) then
      do iter=1, n_triplets
        if(detailed_print_rho) then
          do iter2=1, n_triplets
            help_v = matmul(X_k, i_j(:, iter:iter) - i_k(:,iter2:iter2))
            ave_distance(iter) = ave_distance(iter) + &
            & dot_product(i_j(:, iter) - i_k(:,iter2), help_v(:,1))
          end do
          ave_distance(iter) = ave_distance(iter) / (n_triplets*1.0)
        else
          help_v = matmul(X_k, i_k(:, iter:iter))
          ave_distance(iter) = dot_product(i_k(:, iter), help_v(:,1))
        end if
      end do
      ! print it
      call srpr(" : different clusters(i_k) distance(mean): ", &
      & sum(ave_distance)/(n_triplets*1.0))
      ! ------------------------------ I_J (internal clusters)
      do iter=1, n_triplets
        if(detailed_print_rho) then
          do iter2=1, n_triplets
            help_v = matmul(X_k, i_j(:,iter:iter) - i_j(:,iter:iter))
            ave_distance(iter) = ave_distance(iter) + &
            & dot_product(i_j(:, iter) - i_j(:, iter2), help_v(:,1))
          end do
          ave_distance(iter) = ave_distance(iter) / (n_triplets*1.0)
        else
          help_v = matmul(X_k, i_j(:,iter:iter))
          ave_distance(iter) = dot_product(i_j(:, iter), help_v(:,1))
        end if
      end do
      ! print it
      call srpr(" : same cluster(i_j) distance(mean): ", &
      & sum(ave_distance)/(n_triplets*1.0))
    end if
  end subroutine

!-------------------------------------------------------------------------------
! Objective function
!
  subroutine obj_f(rho_k, obj_f_out)
    implicit none
    real, intent(out) :: obj_f_out
    real, intent(in) :: rho_k
    real :: loss_ou
    call huber_loss(rho_k, loss_ou)
    obj_f_out = rho_k - loss_ou
  end subroutine

!-------------------------------------------------------------------------------
! Huber loss soubroutine
!
  subroutine huber_loss(rho_k, loss_out) !S, n_features
    implicit none
    real, intent(in) :: rho_k
    real, dimension(n_triplets) :: loss_out_v
    real, intent(out) :: loss_out
    real C ! algorithm parameter tha balaces the margin violation
    real Ar_dot_X, delta_rho
    real h
    real helper
    integer ii
    C = 1
    h = 0.5
    do ii=1, n_triplets
      call rho_f(ii, Ar_dot_X)
      delta_rho = Ar_dot_X - rho_k
      if(delta_rho.ge.h) then
        helper = 0.0
      else if(delta_rho.gt.(-h) .and. delta_rho.lt.h) then
        helper = (h - delta_rho**2)/(4*h)
      else if(delta_rho.le.(-h)) then
        helper = -delta_rho
      end if
      loss_out_v(ii) = helper
    end do
    loss_out = C * sum(loss_out_v)
  end subroutine

!-------------------------------------------------------------------------------
! gradient of the objective function (is a matrix it depends only on X_k)
!
  subroutine grad_obj_f(rho_k, grad_obj_f_out)
    implicit none
    real, intent(in) :: rho_k
    real, dimension(n_features, n_features), intent(out) :: grad_obj_f_out
    call huber_loss_deriv(rho_k, grad_obj_f_out)
    grad_obj_f_out = - grad_obj_f_out
  end subroutine


!-------------------------------------------------------------------------------
! HUBER loss derivative
! In comparison to the python code here the derivative of the linear matrix is
! not simply 1 but it depends on the constant matrix in front (e.g. i_k*i_k)
  subroutine huber_loss_deriv(rho_k, deriv_loss_out)
    implicit none
    real, intent(in) :: rho_k
    real, dimension(n_features, n_features), intent(out) :: deriv_loss_out
    real, dimension(n_features, n_features) :: dlo_tmp
    real, dimension(n_features, n_features) :: A_r_mat
    real :: helper
    real :: delta_rho
    real :: Ar_dot_X_deriv
    real :: C, h
    integer ii
    C = 1
    h = 0.5
    do ii=1, n_triplets
      ! print *, "ii", ii
      A_r_mat = matmul(i_k(:,ii:ii), transpose(i_k(:,ii:ii))) - &
      & matmul(i_j(:,ii:ii), transpose(i_j(:,ii:ii)))
      ! print *, "ARMAT", A_r_mat
      call rho_f(ii, Ar_dot_X_deriv)
      delta_rho = Ar_dot_X_deriv - rho_k
      ! this condition must be calculated GENERALLY
      ! Ar_dot_X_deriv = i_k(gi1, ii) * X_k(gi1, gj1) * i_k(gj1, ii) - &
      !                 & i_j(gi1, ii) * X_k(gi1, gj1) * i_j(gj1, ii)
      ! delta_rho = Ar_dot_X_deriv - rho_k
      if(delta_rho .ge. h) then
        helper = 0.0
      else if(delta_rho .gt. (-h) .and. delta_rho .lt. h) then
        helper = (-delta_rho)/(2*h)
      else if(delta_rho .le. (-h)) then
        ! helper = - (i_k(gi1, ii) * 1 * i_k(gj1, ii) - &
        !           & i_j(gi1, ii) * 1 * i_j(gj1, ii))
        helper = -1
      end if
      if(ii .eq. 1) then
        dlo_tmp = helper*A_r_mat
      else
        dlo_tmp = deriv_loss_out + helper*A_r_mat
      end if
      deriv_loss_out = dlo_tmp
      ! print *, 'helper', helper
      ! print *, 'deriv_loss_out', deriv_loss_out
    end do
    deriv_loss_out = C * deriv_loss_out
  end subroutine

!-------------------------------------------------------------------------------
! Backtracking line search algorithm to find the best alpha
! The wolfe condition are indicated by cond1 and cond2
!
! For this implementation (as from the paper) the Wolfe conditions seems not
! to work correctly, i.e. they never stop the loop.
! For this reason I added conditions on alpha that will respect the positiveness
! of the X_k matrix, and they will converge faster (at a certain point even
! only 2 iterations). Even with these modification the Wolfe conditions are
! not effectively exiting the loop.
!
!
  subroutine line_search_alpha(alpha_dir, gradient_f, rho_in, search_dir, &
    & must_positive)
    implicit none
    real, intent(in) :: gradient_f(n_features, n_features)
    real, intent(in) :: search_dir(n_features, n_features)
    real, intent(in) :: rho_in
    ! alpha dir (output) is alpha*search_dir
    real, dimension(n_features, n_features), intent(out) :: alpha_dir
    real  x_supertemp(n_features, n_features)
    real  cond1(n_features, n_features)
    real  cond2(n_features, n_features)
    real  m1(n_features, n_features)
    real  step, c1, c2
    integer, dimension(n_features, n_features) :: a
    real  f1, f2
    integer count, count_max
    logical must_positive
    count_max = 50
    if(.not. must_positive) count_max = 1

    ! Some parameter is inspired by scipy implementation
    c1 = 0.0001 ! 0 < c1 < c2 < 1
    c2 = 0.9
    ! amax = 1.0 ! amax maximum step size (if >1.0 it is exploding here)
    step = 0.9
    alpha_dir = 1.0/step ! std for final alpha
    a = 0
    cond1 = -1 ! default to enter the condition
    cond2 = -1
    count = 1

    ! Calculating X_i and saving X_k for the end
    call obj_f(rho_in, f1) ! X_i, rho_i
    x_supertemp = X_k

    ! main Wolfe condition loop
    do while (any(cond1 .lt. 0.0) .or. any(cond2 .lt. 0.0))
      ! X_i + alpha*p_i
      where(a .eq. 0) alpha_dir = alpha_dir*step
      X_k = x_supertemp + alpha_dir*search_dir

      ! Blocking a from going down more than 0.01 (useless)
      if(any(alpha_dir .lt. 0.01)) then
        where(alpha_dir .lt. 0.01) alpha_dir = 0.0
        where(alpha_dir .lt. 0.01) a = 1
      end if

      ! if the X_k is positive after the first alpha -> found
      where(X_k .gt. 0.0) a = 1

      ! Check to avoid useless loops (original x small and negative direction)
      if(any(x_supertemp .lt. 0.01 .and. search_dir .lt. 0.0) ) then
        where(x_supertemp .lt. 0.01 .and. search_dir .lt. 0.0) a = 1
        where(x_supertemp .lt. 0.01 .and. search_dir .lt. 0.0) alpha_dir = 0.0
      end if

      ! faster convergence of a to 0 if it is next to 0
      if(any(a .eq. 0 .and. alpha_dir .lt. 0.25)) step = 0.6
      ! if(any(a .eq. 0 .and. alpha_dir .lt. 0.25)) print *, "step to 0.6"
      if(any(a .eq. 0 .and. alpha_dir .lt. 0.1)) step = 0.3
      ! if(any(a .eq. 0 .and. alpha_dir .lt. 0.1)) print *, "step to 0.3"

      ! Calculate the new X_k to use (with new direction)
      call obj_f(rho_in, f2)
      call grad_obj_f(rho_in, m1)

      ! Calculate the condition (on one single side)
      cond1 = c1*alpha_dir*matmul(transpose(search_dir), gradient_f) +&
      & f1 - f2
      cond2 = c2*abs(matmul(transpose(search_dir), gradient_f)) - &
      & abs(matmul(transpose(search_dir), m1))
      ! print *, "cond1:", cond1
      ! print *, "cond2:", cond2
      ! print *, count, ": ", any(cond1 .lt. 0.0), any(cond2 .lt. 0.0)

      ! Exit condition for a positive X (a == 1 means that alpha found is right)
      if(all(X_k .gt. 0.0) .and. all(a .eq. 1)) exit

      ! Counting the number of iterations
      if(count .eq. count_max) then
        exit
      end if
      count = count + 1
    end do

    ! If we do not need the metric to be stricly positive we can do this
    if(.not. must_positive) alpha_dir = 0.2

    ! putting back the matrix and defining the alpha_dir for return
    X_k = x_supertemp
    alpha_dir = alpha_dir*search_dir

    ! Eventual print
    ! print *, "======= a_fin:", alpha_dir, "(count=", count, ")"
  end subroutine
end module
