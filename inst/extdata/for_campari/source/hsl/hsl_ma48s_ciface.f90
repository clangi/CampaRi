module hsl_ma48_single_ciface
   use iso_c_binding
   use hsl_ma48_single, only:                                           &
      f_ma48_factors                => ma48_factors,                    &
      f_ma48_control                => ma48_control,                    &
      f_ma48_ainfo                  => ma48_ainfo,                      &
      f_ma48_finfo                  => ma48_finfo,                      &
      f_ma48_sinfo                  => ma48_sinfo,                      &
      f_ma48_initialize             => ma48_initialize,                 &
      f_ma48_analyse                => ma48_analyse,                    &
      f_ma48_get_perm               => ma48_get_perm,                   &
      f_ma48_factorize              => ma48_factorize,                  &
      f_ma48_solve                  => ma48_solve,                      &
      f_ma48_finalize               => ma48_finalize,                   &
      f_ma48_special_rows_and_cols  => ma48_special_rows_and_cols,      &
      f_ma48_determinant            => ma48_determinant
   use hsl_zd11_single
   implicit none

   integer, parameter :: wp = C_FLOAT

   type, bind(C) :: ma48_control
      integer(C_INT) :: f_arrays
      real(wp) :: multiplier
      real(wp) :: u
      real(wp) :: switch
      real(wp) :: drop
      real(wp) :: tolerance
      real(wp) :: cgce
      integer(C_INT) :: lp
      integer(C_INT) :: wp
      integer(C_INT) :: mp
      integer(C_INT) :: ldiag
      integer(C_INT) :: btf
      integer(C_INT) :: struct
      integer(C_INT) :: maxit
      integer(C_INT) :: factor_blocking
      integer(C_INT) :: solve_blas
      integer(C_INT) :: pivoting
      integer(C_INT) :: diagonal_pivoting
      integer(C_INT) :: fill_in
      integer(C_INT) :: switch_mode
   end type ma48_control

   type, bind(C) :: ma48_ainfo
      real(wp) :: ops
      integer(C_INT) :: flag
      integer(C_INT) :: more
      integer(C_LONG) :: lena_analyse
      integer(C_LONG) :: lenj_analyse
      integer(C_LONG) :: lena_factorize
      integer(C_LONG) :: leni_factorize
      integer(C_INT) :: ncmpa
      integer(C_INT) :: rank
      integer(C_LONG) :: drop
      integer(C_INT) :: struc_rank
      integer(C_LONG) :: oor
      integer(C_LONG) :: dup
      integer(C_INT) :: stat
      integer(C_INT) :: lblock
      integer(C_INT) :: sblock
      integer(C_LONG) :: tblock
   end type ma48_ainfo

   type, bind(C) :: ma48_finfo
      real(wp) :: ops
      integer(C_INT) :: flag
      integer(C_INT) :: more
      integer(C_LONG) :: size_factor
      integer(C_LONG) :: lena_factorize
      integer(C_LONG) :: leni_factorize
      integer(C_LONG) :: drop
      integer(C_INT) :: rank
      integer(C_INT) :: stat
   end type ma48_finfo

   type, bind(C) :: ma48_sinfo
      integer(C_INT) :: flag
      integer(C_INT) :: more
      integer(C_INT) :: stat
   end type ma48_sinfo
contains
   subroutine copy_control_in(ccontrol, fcontrol, f_arrays)
      type(ma48_control), intent(in) :: ccontrol
      type(f_ma48_control), intent(out) :: fcontrol
      logical, intent(out) :: f_arrays

      f_arrays = (ccontrol%f_arrays.ne.0)
      fcontrol%multiplier = ccontrol%multiplier
      fcontrol%u = ccontrol%u
      fcontrol%switch = ccontrol%switch
      fcontrol%drop = ccontrol%drop
      fcontrol%tolerance = ccontrol%tolerance
      fcontrol%cgce = ccontrol%cgce
      fcontrol%lp = ccontrol%lp
      fcontrol%wp = ccontrol%wp
      fcontrol%mp = ccontrol%mp
      fcontrol%ldiag = ccontrol%ldiag
      fcontrol%btf = ccontrol%btf
      fcontrol%struct = (ccontrol%struct.ne.0)
      fcontrol%maxit = ccontrol%maxit
      fcontrol%factor_blocking = ccontrol%factor_blocking
      fcontrol%solve_blas = ccontrol%solve_blas
      fcontrol%pivoting = ccontrol%pivoting
      fcontrol%diagonal_pivoting = (ccontrol%diagonal_pivoting.ne.0)
      fcontrol%fill_in = ccontrol%fill_in
      fcontrol%switch_mode = (ccontrol%switch_mode.ne.0)
   end subroutine copy_control_in

   subroutine copy_ainfo_out(fainfo, cainfo)
      type(f_ma48_ainfo), intent(in) :: fainfo
      type(ma48_ainfo), intent(out) :: cainfo

      cainfo%ops = fainfo%ops
      cainfo%flag = fainfo%flag
      cainfo%more = fainfo%more
      cainfo%lena_analyse = fainfo%lena_analyse
      cainfo%lenj_analyse = fainfo%lenj_analyse
      cainfo%lena_factorize = fainfo%lena_factorize
      cainfo%leni_factorize = fainfo%leni_factorize
      cainfo%ncmpa = fainfo%ncmpa
      cainfo%rank = fainfo%rank
      cainfo%drop = fainfo%drop
      cainfo%struc_rank = fainfo%struc_rank
      cainfo%oor = fainfo%oor
      cainfo%dup = fainfo%dup
      cainfo%stat = fainfo%stat
      cainfo%lblock = fainfo%lblock
      cainfo%sblock = fainfo%sblock
      cainfo%tblock = fainfo%tblock
   end subroutine copy_ainfo_out

   subroutine copy_finfo_out(ffinfo, cfinfo)
      type(f_ma48_finfo), intent(in) :: ffinfo
      type(ma48_finfo), intent(out) :: cfinfo

      cfinfo%ops = ffinfo%ops
      cfinfo%flag = ffinfo%flag
      cfinfo%more = ffinfo%more
      cfinfo%size_factor = ffinfo%size_factor
      cfinfo%lena_factorize = ffinfo%lena_factorize
      cfinfo%leni_factorize = ffinfo%leni_factorize
      cfinfo%drop = ffinfo%drop
      cfinfo%rank = ffinfo%rank
      cfinfo%stat = ffinfo%stat
   end subroutine copy_finfo_out

   subroutine copy_sinfo_out(fsinfo, csinfo)
      type(f_ma48_sinfo), intent(in) :: fsinfo
      type(ma48_sinfo), intent(out) :: csinfo

      csinfo%flag = fsinfo%flag
      csinfo%more = fsinfo%more
      csinfo%stat = fsinfo%stat
   end subroutine copy_sinfo_out

end module hsl_ma48_single_ciface

subroutine ma48_default_control_s(ccontrol) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   type(ma48_control), intent(out) :: ccontrol

   type(f_ma48_control) :: fcontrol

   ccontrol%f_arrays = 0 ! false

   call f_ma48_initialize(control=fcontrol)
   ccontrol%multiplier = fcontrol%multiplier
   ccontrol%u = fcontrol%u
   ccontrol%switch = fcontrol%switch
   ccontrol%drop = fcontrol%drop
   ccontrol%tolerance = fcontrol%tolerance
   ccontrol%cgce = fcontrol%cgce
   ccontrol%lp = fcontrol%lp
   ccontrol%wp = fcontrol%wp
   ccontrol%mp = fcontrol%mp
   ccontrol%ldiag = fcontrol%ldiag
   ccontrol%btf = fcontrol%btf
   ccontrol%struct = 0 ! false
   if(fcontrol%struct) ccontrol%struct = 1 ! true
   ccontrol%maxit = fcontrol%maxit
   ccontrol%factor_blocking = fcontrol%factor_blocking
   ccontrol%solve_blas = fcontrol%solve_blas
   ccontrol%pivoting = fcontrol%pivoting
   ccontrol%diagonal_pivoting = 0 ! false
   if(fcontrol%diagonal_pivoting) ccontrol%diagonal_pivoting = 1 ! true
   ccontrol%fill_in = fcontrol%fill_in
   ccontrol%switch_mode = 0 ! false
   if(fcontrol%switch_mode) ccontrol%switch_mode = 1 ! true
end subroutine ma48_default_control_s

subroutine ma48_initialize_s (cfactors) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   type(C_PTR), intent(out) :: cfactors

   type(f_ma48_factors), pointer :: ffactors

   ! Copy data in and associate pointers correctly
   allocate(ffactors)
   cfactors = c_loc(ffactors)

   ! Call the Fortran routine
   call f_ma48_initialize(factors=ffactors)
end subroutine ma48_initialize_s

subroutine ma48_analyse_s (m, n, ne, crow, ccol, cval, cfactors, &
      ccontrol, cainfo, cpfinfo, cperm, cendcol) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   integer(C_INT), value :: m
   integer(C_INT), value :: n
   integer(C_LONG), value :: ne
   type(C_PTR), value :: crow
   type(C_PTR), value :: ccol
   type(C_PTR), value :: cval
   type(C_PTR), value :: cfactors
   type(ma48_control), intent(in) :: ccontrol
   type(ma48_ainfo), intent(out) :: cainfo
   type(C_PTR), value :: cpfinfo
   type(C_PTR), value :: cperm
   type(C_PTR), value :: cendcol

   logical :: f_arrays
   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), pointer :: fcol
   real(wp), dimension(:), pointer :: fval
   type(zd11_type) :: matrix
   type(f_ma48_factors), pointer :: ffactors
   type(f_ma48_control) :: fcontrol
   type(f_ma48_ainfo) :: fainfo
   type(f_ma48_finfo) :: ffinfo
   integer, dimension(:), pointer :: fperm
   integer, dimension(:), allocatable, target :: fperm_alloc
   integer, dimension(:), pointer :: fendcol
   integer, dimension(:), allocatable, target :: fendcol_alloc
   type(ma48_finfo), pointer :: cfinfo

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   matrix%m = m
   matrix%n = n
   if(ne .gt. huge(matrix%ne)) then
      fainfo%flag = -3 ! ne out of range
      if(fcontrol%lp.gt.0) write(fcontrol%lp, "(a)") &
         "HSL_MA48 Error -3: ne must be representable as 32-bit"
      call copy_ainfo_out(fainfo, cainfo)
      return
   endif
   matrix%ne = int(ne, kind=kind(0))
   call C_F_POINTER(crow, frow, shape = (/ ne /))
   call C_F_POINTER(ccol, fcol, shape = (/ ne /))
   call C_F_POINTER(cval, fval, shape = (/ ne /))
   allocate(matrix%row(ne), matrix%col(ne), matrix%val(ne))
   if(f_arrays) then
      matrix%row(:) = frow(:)
      matrix%col(:) = fcol(:)
   else
      matrix%row(:) = frow(:) + 1
      matrix%col(:) = fcol(:) + 1
   endif
   matrix%val(:) = fval(:)
   call C_F_POINTER(cfactors, ffactors)
   if(C_ASSOCIATED(cperm)) then
      call C_F_POINTER(cperm, fperm, shape = (/ m+n /))
      if(.not.f_arrays) then
         allocate(fperm_alloc(m+n))
         fperm_alloc(:) = fperm(:) + 1
         fperm => fperm_alloc
      endif
   else
      nullify(fperm)
   endif
   if(C_ASSOCIATED(cendcol)) then
      call C_F_POINTER(cendcol, fendcol, shape = (/ m+n /))
      if(.not.f_arrays) then
         allocate(fendcol_alloc(m+n))
         fendcol_alloc(:) = fendcol(:) + 1
         fendcol => fendcol_alloc
      endif
   else
      nullify(fendcol)
   endif
   if(C_ASSOCIATED(cpfinfo)) then
      call C_F_POINTER(cpfinfo, cfinfo)
   else
      nullify(cfinfo)
   endif

   ! Call the Fortran routine
   if(associated(cfinfo)) then
      if(associated(fperm)) then
         if(associated(fendcol)) then
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               finfo=ffinfo, perm=fperm, endcol=fendcol)
         else
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               finfo=ffinfo, perm=fperm)
         endif
      else ! .not. associated(fperm)
         if(associated(fendcol)) then
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               finfo=ffinfo, endcol=fendcol)
         else
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               finfo=ffinfo)
         endif
      endif
      ! Copy data out
      call copy_ainfo_out(fainfo, cainfo)
      call copy_finfo_out(ffinfo, cfinfo)
   else ! .not. associated(cfinfo)
      if(associated(fperm)) then
         if(associated(fendcol)) then
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               perm=fperm, endcol=fendcol)
         else
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               perm=fperm)
         endif
      else ! .not. associated(fperm)
         if(associated(fendcol)) then
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo, &
               endcol=fendcol)
         else
            call f_ma48_analyse(matrix, ffactors, fcontrol, fainfo)
         endif
      endif
      ! Copy data out
      call copy_ainfo_out(fainfo, cainfo)
   endif
end subroutine ma48_analyse_s

subroutine ma48_factorize_s(m, n, ne, crow, ccol, cval, cfactors, &
      ccontrol, cfinfo, fast, partial) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   integer(C_INT), value :: m
   integer(C_INT), value :: n
   integer(C_LONG), value :: ne
   type(C_PTR), value :: crow
   type(C_PTR), value :: ccol
   type(C_PTR), value :: cval
   type(C_PTR), value :: cfactors
   type(ma48_control), intent(in) :: ccontrol
   type(ma48_finfo), intent(out) :: cfinfo
   integer(C_INT), value:: fast
   integer(C_INT), value:: partial

   type(zd11_type) :: matrix
   logical :: f_arrays
   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), pointer :: fcol
   real(wp), dimension(:), pointer :: fval
   type(f_ma48_factors), pointer :: ffactors
   type(f_ma48_control) :: fcontrol
   type(f_ma48_finfo) :: ffinfo

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   matrix%m = m
   matrix%n = n
   if(ne .gt. huge(matrix%ne)) then
      ffinfo%flag = -3 ! ne out of range
      if(fcontrol%lp.gt.0) write(fcontrol%lp, "(a)") &
         "HSL_MA48 Error -3: ne must be representable as 32-bit"
      call copy_finfo_out(ffinfo, cfinfo)
      return
   endif
   matrix%ne = int(ne, kind=kind(0))
   call C_F_POINTER(crow, frow, shape = (/ ne /))
   call C_F_POINTER(ccol, fcol, shape = (/ ne /))
   call C_F_POINTER(cval, fval, shape = (/ ne /))
   allocate(matrix%row(ne), matrix%col(ne), matrix%val(ne))
   if(f_arrays) then
      matrix%row(:) = frow(:)
      matrix%col(:) = fcol(:)
   else
      matrix%row(:) = frow(:) + 1
      matrix%col(:) = fcol(:) + 1
   endif
   matrix%val(:) = fval(:)
   call C_F_POINTER(cfactors, ffactors)

   ! Call the Fortran routine
   if(fast.ne.0) then ! fast is true
      if(partial.ne.0) then ! partial is true
         call f_ma48_factorize(matrix, ffactors, fcontrol, ffinfo, &
            fast=1, partial=1)
      else ! partial is false
         call f_ma48_factorize(matrix, ffactors, fcontrol, ffinfo, &
            fast=1)
      endif
   else ! fast is false
      if(partial.ne.0) then ! partial is true
         call f_ma48_factorize(matrix, ffactors, fcontrol, ffinfo, &
            partial=1)
      else ! partial is false
         call f_ma48_factorize(matrix, ffactors, fcontrol, ffinfo)
      endif
   endif

   ! Copy data out
   call copy_finfo_out(ffinfo, cfinfo)

end subroutine ma48_factorize_s

subroutine ma48_solve_s (m, n, ne, crow, ccol, cval, cfactors, crhs, &
      cx, ccontrol, csinfo, trans, cresid, cerror) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   integer(C_INT), value :: m
   integer(C_INT), value :: n
   integer(C_LONG), value :: ne
   type(C_PTR), value :: crow
   type(C_PTR), value :: ccol
   type(C_PTR), value :: cval
   type(C_PTR), value :: cfactors
   type(C_PTR), value:: crhs
   type(C_PTR), value :: cx
   type(ma48_control), intent(in) :: ccontrol
   type(ma48_sinfo), intent(out) :: csinfo
   integer(C_INT), value:: trans
   type(C_PTR), value :: cresid
   type(C_PTR), value :: cerror

   logical :: f_arrays
   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), pointer :: fcol
   real(wp), dimension(:), pointer :: fval
   type(f_ma48_factors), pointer :: ffactors
   type(f_ma48_control) :: fcontrol
   type(zd11_type) :: matrix
   real(wp), dimension(:), pointer :: frhs
   real(wp), dimension(:), pointer :: fx
   type(f_ma48_sinfo) :: fsinfo
   real(wp), dimension(:), pointer :: fresid
   real(wp), pointer :: ferror

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   matrix%m = m
   matrix%n = n
   if(ne .gt. huge(matrix%ne)) then
      fsinfo%flag = -3 ! ne out of range
      if(fcontrol%lp.gt.0) write(fcontrol%lp, "(a)") &
         "HSL_MA48 Error -3: ne must be representable as 32-bit"
      call copy_sinfo_out(fsinfo, csinfo)
      return
   endif
   matrix%ne = int(ne, kind=kind(0))
   call C_F_POINTER(crow, frow, shape = (/ ne /))
   call C_F_POINTER(ccol, fcol, shape = (/ ne /))
   call C_F_POINTER(cval, fval, shape = (/ ne /))
   allocate(matrix%row(ne), matrix%col(ne), matrix%val(ne))
   if(f_arrays) then
      matrix%row(:) = frow(:)
      matrix%col(:) = fcol(:)
   else
      matrix%row(:) = frow(:) + 1
      matrix%col(:) = fcol(:) + 1
   endif
   matrix%val(:) = fval(:)
   call C_F_POINTER(cfactors, ffactors)
   call C_F_POINTER(crhs, frhs, shape = (/ n /))
   call C_F_POINTER(cx, fx, shape = (/ n /))
   if(C_ASSOCIATED(cresid)) then
      call C_F_POINTER(cresid, fresid, shape = (/ 2 /))
   else
      nullify(fresid)
   endif
   if(C_ASSOCIATED(cerror)) then
      call C_F_POINTER(cerror, ferror)
   else
      nullify(ferror)
   endif

   ! Call the Fortran routine
   if(trans.ne.0) then ! trans is true
      if(associated(fresid)) then
         if(associated(ferror)) then
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, trans=1, resid=fresid, error=ferror)
         else ! .not.associated(ferror)
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, trans=1, resid=fresid)
         endif
      else ! .not.associated(fresid)
         if(associated(ferror)) then
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, trans=1, error=ferror)
         else ! .not.associated(ferror)
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, trans=1)
         endif
      endif
   else ! trans is false
      if(associated(fresid)) then
         if(associated(ferror)) then
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, resid=fresid, error=ferror)
         else ! .not.associated(ferror)
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, resid=fresid)
         endif
      else ! .not.associated(fresid)
         if(associated(ferror)) then
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo, error=ferror)
         else ! .not.associated(ferror)
            call f_ma48_solve(matrix, ffactors, frhs, fx, fcontrol, &
               fsinfo)
         endif
      endif
   endif

   ! Copy data out
   call copy_sinfo_out(fsinfo, csinfo)
end subroutine ma48_solve_s

integer(C_INT) function ma48_finalize_s (cfactors, ccontrol) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   type(C_PTR) :: cfactors
   type(ma48_control), intent(in) :: ccontrol

   type(f_ma48_factors), pointer :: ffactors
   type(f_ma48_control) :: fcontrol
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cfactors, ffactors)

   ! Call the Fortran routine
   call f_ma48_finalize(ffactors, fcontrol, ma48_finalize_s)

   ! Free memory
   deallocate(ffactors)
   cfactors = C_NULL_PTR

end function ma48_finalize_s

integer(C_INT) function ma48_determinant_s (cfactors, sgndet, logdet, &
      ccontrol) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   type(C_PTR), value:: cfactors
   integer(C_INT), intent(out) :: sgndet
   real(wp), intent(out) :: logdet
   type(ma48_control), intent(in) :: ccontrol

   type(f_ma48_factors), pointer :: ffactors
   type(f_ma48_control) :: fcontrol
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cfactors, ffactors)

   ! Call the Fortran routine
   call f_ma48_determinant(ffactors, sgndet, logdet, fcontrol, &
      ma48_determinant_s)

end function ma48_determinant_s

integer(C_INT) function ma48_special_rows_and_cols_s (cfactors, rank, crows, &
      ccols, ccontrol) bind(C)
   use hsl_ma48_single_ciface
   implicit none

   type(C_PTR), value:: cfactors
   integer(C_INT), intent(out) :: rank
   type(C_PTR), value :: crows
   type(C_PTR), value :: ccols
   type(ma48_control), intent(in) :: ccontrol

   type(f_ma48_factors), pointer :: ffactors
   type(f_ma48_control) :: fcontrol
   logical :: f_arrays
   integer(C_INT), dimension(:), pointer :: frows
   integer(C_INT), dimension(:), pointer :: fcols

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cfactors, ffactors)
   call C_F_POINTER(crows, frows, shape = (/ ffactors%m /))
   call C_F_POINTER(ccols, fcols, shape = (/ ffactors%n /))

   ! Call the Fortran routine
   call f_ma48_special_rows_and_cols(ffactors, rank, frows, fcols, &
      fcontrol, ma48_special_rows_and_cols_s)

   ! Copy data out
   if(.not.f_arrays) then
      frows(:) = frows(:) - 1
      fcols(:) = fcols(:) - 1
   endif

end function ma48_special_rows_and_cols_s
