module lapack_eig
  character, dimension(1) :: jobz = 'V' ! give back also eigenvectors (vs 'N')
  character, dimension(1) :: range = 'I' ! I is index mode of selected eigenvals
  integer ndim ! order of the matrix
  real :: diag ! diagonal
  integer :: iu = 1 ! largest eigenval
  integer :: il ! to set ndim - 1
  integer :: M = 1 ! iu - il + 1
  integer :: ldz = ndim
  ! On what these lapack functions do:
  ! http://www.netlib.org/lapack/lug/node30.html
  ! The one function needed for computing selected eigenvals
  ! http://www.netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_ga323734560b8bd052fbc474dc2f0b5605.html#ga323734560b8bd052fbc474dc2f0b5605
  ! The scipy source of the call to lapack
  ! https://github.com/scipy/scipy/blob/v0.19.1/scipy/linalg/decomp.py#L241-L445
  integer :: info

contains
  subroutine find_eigs_l(eval, evec)
    implicit none
    real, intent(inout) :: eval(1) ! OUT eigenval
    real, intent(inout) :: evec(ndim, M) ! OUT eigenvec

    call dstevr(jobz, range, ndim)


  end subroutine

end module
