!*******************************************************************************
!
!! SNSIMP is a simple program to call ARPACK for a symmetric eigenproblem.
!
!     This example program is intended to illustrate the
!     simplest case of using ARPACK in considerable detail.
!     This code may be used to understand basic usage of ARPACK
!     and as a template for creating an interface to ARPACK.
!
!     This code shows how to use ARPACK to find a few eigenvalues
!     (lambda) and corresponding eigenvectors (x) for the standard
!     eigenvalue problem:
!
!                        A*x = lambda*x
!
!     where A is a n by n real nonsymmetric matrix.
!
!     The main points illustrated here are
!
!        1) How to declare sufficient memory to find NEV
!           eigenvalues of largest magnitude.  Other options
!           are available.
!
!        2) Illustration of the reverse communication interface
!           needed to utilize the top level ARPACK routine SNAUPD
!           that computes the quantities needed to construct
!           the desired eigenvalues and eigenvectors(if requested).
!
!        3) How to extract the desired eigenvalues and eigenvectors
!           using the ARPACK routine SNEUPD.
!
!     The only thing that must be supplied in order to use this
!     routine on your problem is to change the array dimensions
!     appropriately, to specify WHICH eigenvalues you want to compute
!     and to supply a matrix-vector product
!
!                         w <-  Av
!
!     in place of the call to AV( )  below.
!
!     Once usage of this routine is understood, you may wish to explore
!     the other available options to improve convergence, to solve generalized
!     problems, etc.  Look at the file ex-nonsym.doc in DOCUMENTS directory.
!     This codes implements
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is obtained from the standard central difference
!         discretization of the convection-diffusion operator
!                 (Laplacian u) + rho*(du / dx)
!         on the unit square, with zero Dirichlet boundary condition.
!
!     ... OP = A  and  B = I.
!     ... Assume "call av (nx,x,y)" computes y = A*x
!     ... Use mode 1 of SNAUPD.
!
!\BeginLib
!
!\Routines called:
!     snaupd  ARPACK reverse communication interface routine.
!     sneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     tv      Matrix vector multiplication routine that computes T*x,
!             where T is a tridiagonal matrix.  It is used in routine
!             av.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: nsimp.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | The maximum dimensions for all arrays are            |
!     | set here to accommodate a problem size of            |
!     | N <= MAXN                                            |
!     |                                                      |
!     | NEV is the number of eigenvalues requested.          |
!     |     See specifications for ARPACK usage below.       |
!     |                                                      |
!     | NCV is the largest number of basis vectors that will |
!     |     be used in the Implicitly Restarted Arnoldi      |
!     |     Process.  Work per major iteration is            |
!     |     proportional to N*NCV*NCV.                       |
!     |                                                      |
!     | You must set:                                        |
!     |                                                      |
!     | MAXN:   Maximum dimension of the A allowed.          |
!     | MAXNEV: Maximum NEV allowed.                         |
!     | MAXNCV: Maximum NCV allowed.                         |
!     |                                                      |
!     | Specifications for ARPACK usage are set              |
!     | below:                                               |
!     |------------------------------------------------------|
!     |    1) NEV = 4  asks for 4 eigenvalues to be          |
!     |       computed.                                      |
!     |                                                      |
!     |    2) NCV = 20 sets the length of the Arnoldi        |
!     |       factorization.                                 |
!     |                                                      |
!     |    3) This is a standard problem.                    |
!     |         (indicated by bmat  = 'I')                   |
!     |                                                      |
!     |    4) Ask for the NEV eigenvalues of                 |
!     |       largest magnitude.                             |
!     |         (indicated by which = 'LM')                  |
!     |       See documentation in SNAUPD for the            |
!     |       other options SM, LR, SR, LI, SI.              |
!     |                                                      |
!     | Note: NEV and NCV must satisfy the following         |
!     | conditions:                                          |
!     |              NEV <= MAXNEV                           |
!     |          NEV + 2 <= NCV <= MAXNCV                    |
!     |                                                      |
!     |------------------------------------------------------|
!     |                                                      |
!     | Specification of stopping rules and initial          |
!     | conditions before calling SNAUPD                     |
!     |                                                      |
!     | TOL  determines the stopping criterion.              |
!     |                                                      |
!     |      Expect                                          |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC)  |
!     |               computed   true                        |
!     |                                                      |
!     |      If TOL <= 0,  then TOL <- macheps               |
!     |           (machine precision) is used.               |
!     |                                                      |
!     | IDO  is the REVERSE COMMUNICATION parameter          |
!     |      used to specify actions to be taken on return   |
!     |      from SNAUPD. (see usage below)                  |
!     |                                                      |
!     |      It MUST initially be set to 0 before the first  |
!     |      call to SNAUPD.                                 |
!     |                                                      |
!     | INFO on entry specifies starting vector information  |
!     |      and on return indicates error codes             |
!     |                                                      |
!     |      Initially, setting INFO=0 indicates that a      |
!     |      random starting vector is requested to          |
!     |      start the ARNOLDI iteration.  Setting INFO to   |
!     |      a nonzero value on the initial call is used     |
!     |      if you want to specify your own starting        |
!     |      vector (This vector must be placed in RESID).   |
!     |                                                      |
!     | The work array WORKL is used in SNAUPD as            |
!     | workspace.  Its dimension LWORKL is set as           |
!     | illustrated below.                                   |
!     |                                                      |
!     %------------------------------------------------------%
!
module arpack_eig
  logical :: superverbose = .false.
  logical :: verbose = .false.
  real, allocatable :: A(:,:)
  integer nx ! dimension of matrix A
  integer maxn ! maximum number of dimensions for A
  integer maxnev ! maximum number of requested eigenvalues
  integer maxncv ! maximum number of requested eigenvectors
  integer ldv ! maxn
  parameter (maxn=25600, maxnev=1200, maxncv=3000, ldv=maxn) !256,12,30
  integer nev ! number of eigenvalues
  integer ncv ! length of the Arnoldi factorization
  real zero
  real tol ! tolerance for convergence. If <=0 machine precision
  parameter (zero = 0.0E+0)
  character bmat*1, which*2 ! another way to initialize it!
  ! 'LM' -> want the NEV eigenvalues of largest magnitude.
  ! 'SM' -> want the NEV eigenvalues of smallest magnitude.
  ! 'LR' -> want the NEV eigenvalues of largest real part.
  ! 'SR' -> want the NEV eigenvalues of smallest real part.
  ! 'LI' -> want the NEV eigenvalues of largest imaginary part.
  ! 'SI' -> want the NEV eigenvalues of smallest imaginary part.

  ! At present there is no a-priori analysis to guide the selection
  ! of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
  ! However, it is recommended that NCV .ge. 2*NEV+1.

  contains
    subroutine find_eigens(eigval, eigvec)
    ! implicit none !it is craching with some arpack definition

    !     %--------------%
    !     | input variab |
    !     %--------------%
          real eigval(nev,3), eigvec(nev, nx)

    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
          integer iparam(11), ipntr(14)
          logical select(maxncv)
          Real ax(maxn), resid(maxn), workd(3*maxn), workev(3*maxncv), &
          workl(3*maxncv*maxncv+6*maxncv)
          real v(ldv, maxncv) ! max number of vars (nx*nx) ,  max num Arnold factorizations
          real d(maxncv,3)
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
          integer ido, n, lworkl, info, ierr, &
          j, ishfts, maxitr, mode1, nconv
          Real sigmar, sigmai
          logical first, rvec
    !     %-----------------------------%
    !     | BLAS & LAPACK routines used |
    !     %-----------------------------%
          Real slapy2, snrm2
          external slapy2, snrm2, saxpy
    !     %--------------------%
    !     | Intrinsic function |
    !     %--------------------%
          intrinsic abs
    !     %-------------------------------------------------%
    !     | The following sets dimensions for this problem. |
    !     %-------------------------------------------------%
          n = nx*nx
    !     %-------------------------%
    !     | Some check for the vars |
    !     %-------------------------%

          bmat  = 'I'
          which = 'LM' !order? largest major
          tol    = zero
          if ( n > maxn ) then
             print *, ' ERROR with _NSIMP: N is greater than MAXN '
             stop
          else if ( nev > maxnev ) then
             print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
             stop
          else if ( ncv > maxncv ) then
             print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
             stop
          end if
          if (ncv .lt. (2*nev) .or. ncv .gt. nx*nx) then
            print *, 'ERROR: ncv should be higher than 2*nev (+ 1) and greater than nx*nx'
            stop
          end if
    !     %-------------------------------------------------%
    !     | The following include statement and assignments |
    !     | initiate trace output from the internal         |
    !     | actions of ARPACK.  See debug.doc in the        |
    !     | DOCUMENTS directory for usage.  Initially, the  |
    !     | most useful information will be a breakdown of  |
    !     | time spent in the various stages of computation |
    !     | given by setting mnaupd = 1.                    |
    !     %-------------------------------------------------%
          ndigit = -3
          logfil = 6
          mnaitr = 0
          mnapps = 0
          mnaupd = 1
          mnaup2 = 0
          mneigh = 0
          mneupd = 0
          ! other different strange things
          lworkl  = 3*ncv**2+6*ncv
          ido     = 0
          info    = 0
    !     %---------------------------------------------------%
    !     | Specification of Algorithm Mode:                  |
    !     |                                                   |
    !     | This program uses the exact shift strategy        |
    !     | (indicated by setting IPARAM(1) = 1).             |
    !     | IPARAM(3) specifies the maximum number of Arnoldi |
    !     | iterations allowed.  Mode 1 of SNAUPD is used     |
    !     | (IPARAM(7) = 1). All these options can be changed |
    !     | by the user. For details see the documentation in |
    !     | SNAUPD.                                           |
    !     %---------------------------------------------------%
          ishfts = 1
          maxitr = 300
          mode1 = 1

          iparam(1) = ishfts
          iparam(3) = maxitr
          iparam(7) = mode1
    !     %-------------------------------------------%
    !     | M A I N   L O O P (Reverse communication) |
    !     %-------------------------------------------%
     10   continue
    !        %---------------------------------------------%
    !        | Repeatedly call the routine SNAUPD and take |
    !        | actions indicated by parameter IDO until    |
    !        | either convergence is indicated or maxitr   |
    !        | has been exceeded.                          |
    !        %---------------------------------------------%
             call snaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                           v, ldv, iparam, ipntr, workd, workl, lworkl, &
                           info )
             if (ido == -1 .or. ido == 1) then
    !           %-------------------------------------------%
    !           | Perform matrix vector multiplication      |
    !           |                y <--- Op*x                |
    !           | The user should supply his/her own        |
    !           | matrix vector multiplication routine here |
    !           | that takes workd(ipntr(1)) as the input   |
    !           | vector, and return the matrix vector      |
    !           | product to workd(ipntr(2)).               |
    !           %-------------------------------------------%
                call av (nx, workd(ipntr(1)), workd(ipntr(2)))
    !           %-----------------------------------------%
    !           | L O O P   B A C K to call SNAUPD again. |
    !           %-----------------------------------------%
                go to 10
             endif
    !     %----------------------------------------%
    !     | Either we have convergence or there is |
    !     | an error.                              |
    !     %----------------------------------------%
          if ( info < 0 ) then
    !        %--------------------------%
    !        | Error message, check the |
    !        | documentation in SNAUPD. |
    !        %--------------------------%
             print *, ' Error with SNAUPD, info = ', info
             print *, ' Check the documentation of snaupd.'
             print *, ' This is the main Ritz vectorization.'
             if (info .eq. (-3)) then
               if(superverbose) print *, "Remember: NCV-NEV >= 2 and less than or equal to N."
               if(superverbose) print *, 'ncv: ', ncv
               if(superverbose) print *, 'nev: ', nev
               if(superverbose) print *, 'n: ', n ! nx*nx
             end if
             stop
          else
    !        %-------------------------------------------%
    !        | No fatal errors occurred.                 |
    !        | Post-Process using SNEUPD.                |
    !        |                                           |
    !        | Computed eigenvalues may be extracted.    |
    !        |                                           |
    !        | Eigenvectors may be also computed now if  |
    !        | desired.  (indicated by rvec = .true.)    |
    !        |                                           |
    !        | The routine SNEUPD now called to do this  |
    !        | post processing (Other modes may require  |
    !        | more complicated post processing than     |
    !        | mode1,)                                   |
    !        %-------------------------------------------%
             rvec = .true.
            ! before the sneupd the v vector is populated with the Arnoldi basis vector
             call sneupd ( rvec, 'A', select, d, d(1,2), v, ldv, &
                  sigmar, sigmai, workev, bmat, n, which, nev, tol, &
                  resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
                  lworkl, ierr )
                  ! if(superverbose) print *, "Real parts of eigenvalues:", d(:,1)
             if(superverbose) print *, "The corresponding eigenvectors are &
             returned in the first ", iparam(5), " columns (iparam(5))."
    !        %------------------------------------------------%
    !        | The real parts of the eigenvalues are returned |
    !        | in the first column of the two dimensional     |
    !        | array D, and the IMAGINARY part are returned   |
    !        | in the second column of D.  The corresponding  |
    !        | eigenvectors are returned in the first         |
    !        | NCONV (= IPARAM(5)) columns of the two         |
    !        | dimensional array V if requested.  Otherwise,  |
    !        | an orthogonal basis for the invariant subspace |
    !        | corresponding to the eigenvalues in D is       |
    !        | returned in V.                                 |
    !        %------------------------------------------------%
            !  do i=1, iparam(5)
            !    if(superverbose) print *, i,' eigenvector:   ', v(1:nx,i)
            !  end do
             if ( ierr /= 0) then
    !           %------------------------------------%
    !           | Error condition:                   |
    !           | Check the documentation of SNEUPD. |
    !           %------------------------------------%
                print *, ' Error with SNEUPD, info (ierr) = ', ierr
                print *, ' Check the documentation of sneupd. '
                print *, ' This is the egenval/vec extractor. '
                if (info .eq. (-3)) then
                  if(superverbose) print *, "Remember: NCV-NEV >= 2 and less than or equal to N."
                  if(superverbose) print *, 'ncv: ', ncv
                  if(superverbose) print *, 'nev: ', nev
                  if(superverbose) print *, 'n: ', n ! nx*nx
                end if
                stop
             else
                first = .true.
                nconv =  iparam(5)
                do 20 j=1, nconv
    !              %---------------------------%
    !              | Compute the residual norm |
    !              |                           |
    !              |   ||  A*x - lambda*x ||   |
    !              |                           |
    !              | for the NCONV accurately  |
    !              | computed eigenvalues and  |
    !              | eigenvectors.  (IPARAM(5) |
    !              | indicates how many are    |
    !              | accurate to the requested |
    !              | tolerance)                |
    !              %---------------------------%
                   if (d(j,2) == zero)  then ! Ritz value is real
                      call av(nx, v(1,j), ax)
                      call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                      d(j,3) = snrm2(n, ax, 1)
                      d(j,3) = d(j,3) / abs(d(j,1))
                   else if (first) then
    !
    !                 %------------------------%
    !                 | Ritz value is complex. |
    !                 | Residual of one Ritz   |
    !                 | value of the conjugate |
    !                 | pair is computed.      |
    !                 %------------------------%
    !
                      call av(nx, v(1,j), ax)
                      call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                      call saxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                      d(j,3) = snrm2(n, ax, 1)
                      call av(nx, v(1,j+1), ax)
                      call saxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                      call saxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                      d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                      d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
                      d(j+1,3) = d(j,3)
                      first = .false.
                   else
                      first = .true.
                   end if
     20         continue
                !  Display computed residuals.
                ! call smout(6, nconv, 3, d, maxncv, -6, &
                !      'Ritz values (Real, Imag) and residual residuals')
             end if
             !  Print additional convergence information.
             if ( info == 1) then
               if(superverbose) print *, ' Maximum number of iterations reached (info == 1).'
             else if ( info == 3) then
               if(superverbose) print *, ' No shifts could be applied during implicit', &
                                         ' Arnoldi update, try increasing NCV.'
             end if
             if(superverbose) print *, '% ======= %'
             if(superverbose) print *, '| Summary |'
             if(superverbose) print *, '% ======= %'
             if(superverbose) print *, 'Size of the matrix is ', n
             if(superverbose) print *, 'The number of Ritz values (eigenval) requested is ', nev
             if(superverbose) print *, 'The number of Arnoldi vectors (eigenvec) generated', &
                      ' (NCV) is ', ncv
             if(superverbose) print *, 'What portion of the spectrum: ', which
             if(superverbose) print *, 'The number of converged Ritz values is ', nconv
             if(superverbose) print *, 'The number of Implicit Arnoldi update', &
                      ' iterations taken is ', iparam(3)
             if(superverbose) print *, 'The number of A*x (matrix multiplication) is ', iparam(9)
             if(superverbose) print *, 'The convergence criterion is ', tol
             if(superverbose) print *, 'The eigenvectors v are present in the shape of: ', size(v(:,1)), size(v(1,:))
            !  do i=1, iparam(5)
            !    if(superverbose) print *, i,' eigenvector:   ', v(1:(nx),i)
            !  end do
             do i=1,nev
               if(verbose) print *, i,' eigenvalue (Re, Im ,Residual):   ', d(i,1), d(i,2), d(i,3)
               if(verbose) print *, i,' eigenvector:   ', v(1:(nx),i)
               eigvec(i,:) = v(1:nx,i)
               eigval(i, 1) = d(i,1) ! Re
               eigval(i, 2) = d(i,2) ! Im
               eigval(i, 3) = d(i,3) ! Residual
             end do
             if(superverbose) print *, ' '
          end if
     9000 continue
          end

    !*******************************************************************************
    ! This function is the main input to our problems
    ! You can define a new product or a particolar relationship
    subroutine av (nx, x, y)
      implicit none
      integer nx, j
      Real x(nx), y(nx), help(nx,1)
      Real one
      parameter (one = 1.0E+0)
      help(:,1) = x
      help = matmul(A, help) ! simple matrix A * x multiplication
      y = help(:,1)
      return
    end
end module
