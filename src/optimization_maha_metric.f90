module optimization_maha_metric
  ! main variables
  logical :: print_rho_components
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
    implicit none
    integer, intent(in) :: which
    real, intent(out) :: rho_out
    real :: help_v(n_features, 1)
    real :: help_v2(n_features, 1)
    real :: help_v3(n_features, 1)
    real :: ave_distance(n_triplets)
    integer iter
    help_v3(:,1) = i_k(:,which)
    help_v = matmul(X_k, help_v3) ! to check the which
    help_v3(:,1) = i_j(:,which)
    help_v2 = matmul(X_k, help_v3)
    rho_out = dot_product(i_k(:, which), help_v(:,1)) - dot_product(i_j(:, which), help_v2(:,1))
    if (print_rho_components)then
      do iter=1, n_triplets
        help_v3(:,1) = i_k(:,iter)
        help_v = matmul(X_k, help_v3) ! to check the which
        ave_distance(iter) = dot_product(i_k(:, iter), help_v(:,1))
      end do
      print *, "i_k distance(mean): ", sum(ave_distance)/(n_triplets*1.0)
      do iter=1, n_triplets
        help_v3(:,1) = i_j(:,iter)
        help_v = matmul(X_k, help_v3) ! to check the which
        ave_distance(iter) = dot_product(i_j(:, iter), help_v(:,1))
      end do
      print *, "i_j distance(mean): ", sum(ave_distance)/(n_triplets*1.0)
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
! gradient of the objective function (is a matrix it depends only on X_k)
!
  subroutine grad_obj_f(rho_k, grad_obj_f_out)
    implicit none
    real, intent(in) :: rho_k
    real, dimension(n_features, n_features), intent(out) :: grad_obj_f_out
    integer gi, gj !iterations for the matrix
    real deriv_loss_ou
    do gi=1,n_features
      do gj=1,n_features
        call huber_loss_deriv(rho_k, gi, gj, deriv_loss_ou)
        grad_obj_f_out(gi,gj) = - deriv_loss_ou
      end do
    end do
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
        helper = (h-delta_rho**2)/(4*h)
      else if(delta_rho.le.(-h)) then
        helper = -delta_rho
      end if
      loss_out_v(ii) = helper
    end do
    loss_out = C * sum(loss_out_v)
  end subroutine

!-------------------------------------------------------------------------------
! HUBER loss derivative
! In comparison to the python code here the derivative of the linear matrix is
! not simply 1 but it depends on the constant matrix in front (e.g. i_k*i_k)
  subroutine huber_loss_deriv(rho_k, gi1, gj1, deriv_loss_out)
    implicit none
    real, intent(in) :: rho_k
    integer, intent(in) :: gi1, gj1
    real, dimension(n_triplets) :: deriv_loss_out_v
    real, intent(out) :: deriv_loss_out
    real :: delta_rho
    real :: Ar_dot_X_deriv
    real :: C, h
    real helper
    integer sum1, sum2
    integer ii
    C = 1
    h = 0.5
    do ii=1, n_triplets
      Ar_dot_X_deriv = i_k(gi1, ii) * X_k(gi1, gj1) * i_k(gj1, ii) - &
                      & i_j(gi1, ii) * X_k(gi1, gj1) * i_j(gj1, ii)
      delta_rho = Ar_dot_X_deriv - rho_k
      if(delta_rho.ge.h) then
        helper = 0.0
      else if(delta_rho.gt.(-h) .and. delta_rho.lt.h) then
        helper = (-delta_rho)/(2*h)
      else if(delta_rho.le.(-h)) then
        helper = - (i_k(gi1, ii) * 1 * i_k(gj1, ii) - &
                  & i_j(gi1, ii) * 1 * i_j(gj1, ii))
        ! print *, helper
        ! helper = -1 ! POSSIBLE MODIFICATION
      end if
      deriv_loss_out_v(ii) = helper
    end do
    deriv_loss_out = C * sum(deriv_loss_out_v)
  end subroutine

!-------------------------------------------------------------------------------
! Backtracking line search algorithm to find the best alpha
! The wolfe condition are indicated by cond1 and cond2
!
  subroutine line_search_alpha(alpha_fin, gradient_f, rho_in, search_dir)
    implicit none
    real, intent(out) :: alpha_fin
    real, intent(in) :: gradient_f(n_features, n_features)
    real, intent(in) :: search_dir(n_features, n_features)
    real, intent(in) :: rho_in
    real  x_supertemp(n_features, n_features)
    real  cond1(n_features, n_features)
    real  cond2(n_features, n_features)
    real  m1(n_features, n_features)
    real  step, c1, c2
    real  a
    real  f1, f2
    integer count
    c1 = 0.0001!0 < c1 < c2 < 1
    c2 = 0.9999
    step = 0.8
    alpha_fin = 1
    a = 1.0/step
    cond1 = -1
    cond2 = -1
    count = 1
    call obj_f(rho_in, f1) ! X_i, rho_i
    x_supertemp = X_k
    do while (any(cond1 .lt. 0.0) .or. any(cond2 .lt. 0.0))
      ! X_i + alpha*p_i
      a = a*step
      ! print *, "----- ----- ----- -----"
      ! print *, "----- a:", a
      ! print *, "----- ----- ----- -----"
      X_k = x_supertemp + a*search_dir
      call obj_f(rho_in, f2)
      call grad_obj_f(rho_in,m1)
      cond1 = c1*a*matmul(transpose(search_dir), gradient_f) + f1 - f2
      cond2 = c2*abs(matmul(transpose(search_dir), gradient_f)) - &
      & abs(matmul(transpose(search_dir), m1))
      ! print *, "cond1:", cond1
      ! print *, count, ": ", any(cond1 .lt. 0.0), any(cond2 .lt. 0.0)
      ! print *, "cond2:", cond2
      count = count + 1
      if(count .eq. 50) then
        ! print *, "max count 50 reached"
        exit
      end if
    end do
    X_k = x_supertemp
    alpha_fin = a
    alpha_fin = 0.8
    ! print *, "======= a_fin:", alpha_fin, "(count=",count,")"
  end subroutine
end module
