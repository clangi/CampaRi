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
  use gutenberg
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
    call sl()
    call sipr('Initialized PRNG with Seed of ',seed)
    call sl()
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
  use gutenberg
  implicit none
!
  call sl()
  call spr('------------------------------------------------>')
  call sl()
  call spr('CAMPARI CRASHED UNEXPECTEDLY. PLEASE RECORD ANY ')
  call spr(' INFORMATION ABOUT THE PROBLEM PROVIDED ABOVE.')
  call sl()
  call spr('<------------------------------------------------')
!
  ! stop 98
!
end
!
!-----------------------------------------------------------------------


#ifdef LINK_NETCDF

!------------------------------------------------------------------------------
!
! #ifdef LINK_NETCDF
!
subroutine dump_nbl_nc()
!
  use m_mst_dumping
  use m_variables_gen
  use netcdf
!
  implicit none
!
  integer i,ncid,ii,jj,freeunit,xlen,istart
  integer, ALLOCATABLE :: helper(:)
  logical exists
  character(len = 100) attstring, dumpfile
  real, ALLOCATABLE :: prthlp(:)
  integer :: cnc_ids(10)
!
  dumpfile = 'MST_DUMPLING.nc' !TODO custom filename and more attributes
  inquire(file=dumpfile,exist=exists)
  if (exists.EQV..true.) then !manual delete of the file if it exists
    call spr('Dumping file MST_DUMPLING.nc already existent. It will be overwritten')
    ncid = freeunit()
    open(unit=ncid,file=dumpfile,status='old')
    close(unit=ncid,status='delete')
  end if
  ncid = freeunit()
  call check_err( nf90_create(path=dumpfile, cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
!
! enable definition
  do i=1,100
    attstring(i:i) = " "
  end do
  attstring(1:7) = "CAMPARI"
  call check_err( nf90_put_att(ncid, NF90_GLOBAL, "program", attstring(1:7)) )
! define dimensions
  attstring(1:10) = "framepairs"
  call check_err( nf90_def_dim(ncid, attstring(1:10), sum(approxmst(1:n_snaps)%deg), cnc_ids(1)) )
  attstring(1:10) = "     "
! define (not set) variables to hold distance and type of distance information
  attstring(1:9) = "snapshots"
  call check_err( nf90_def_var(ncid, attstring(1:9), NF90_INT, cnc_ids(1), cnc_ids(2)) )
  attstring(1:9) = "neighbors"
  call check_err( nf90_def_var(ncid, attstring(1:9), NF90_INT, cnc_ids(1), cnc_ids(3)) )
  attstring(1:9) = "distances"
  call check_err( nf90_def_var(ncid, attstring(1:9), NF90_FLOAT, cnc_ids(1), cnc_ids(4)) )
  attstring(1:9) = "         "
  if (dis_method.eq.1) then
    xlen = 14
    attstring(1:xlen) = "torsional RMSD"
  else if (dis_method.eq.5) then
    xlen = 21
    attstring(1:xlen) = "unaligned atomic RMSD"
  else
    call spr('Fatal. Unsupported distance criterion in dump_nbl_nc(...). This is an omission bug.')
    call fexit()
  end if
  call check_err( nf90_put_att(ncid, cnc_ids(1), "type", attstring(1:xlen)) )
!  call check_err( nf90_def_var(ncid, attstring(1:xlen), NF90_FLOAT, cnc_ids(2), cnc_ids(5)) )
  attstring(1:5) = "units"
  if (dis_method.le.2) then
    attstring(6:12) = "degrees"
    call check_err( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:12)) )
  else if (dis_method.eq.5) then
    attstring(6:13) = "angstrom"
    call check_err( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:13)) )
  end if
  attstring(1:13) = "               "
! quit define mode
  call check_err( nf90_enddef(ncid) )
  call check_err( nf90_sync(ncid) )
  ! write(ilog,*) "Definition completed"
!
! put the data
  allocate(helper(n_snaps))
  allocate(prthlp(n_snaps))

  istart = 1
  do i=1,n_snaps
    if (approxmst(i)%deg.le.0) cycle
    helper(1:approxmst(i)%deg) = i
    prthlp(1:approxmst(i)%deg) = real(approxmst(i)%dist(1:approxmst(i)%deg),KIND=4)
    call check_err( nf90_put_var(ncid, cnc_ids(2), helper(1:approxmst(i)%deg), &
 &                                       start = (/ istart /), count = (/ approxmst(i)%deg /)) )
    call check_err( nf90_put_var(ncid, cnc_ids(3), approxmst(i)%adj(1:approxmst(i)%deg), &
 &                                       start = (/ istart /), count = (/ approxmst(i)%deg /)) )
    call check_err( nf90_put_var(ncid, cnc_ids(4), real(prthlp(1:approxmst(i)%deg),KIND=4), &
 &                                       start = (/ istart /), count = (/ approxmst(i)%deg /)) )
    istart = istart + approxmst(i)%deg
  end do
!
  deallocate(prthlp)
  deallocate(helper)
!   do i=1,1
!     write(ilog,*) "DEG",approxmst(i)%deg
!     write(ilog,*) "ADJ",approxmst(i)%adj(:)
!     write(ilog,*) "DIST",approxmst(i)%dist(:)
!   end do
! ! close
  call check_err( nf90_close(ncid) )
  call sipr('Saved a total number of snapshots of: ', sum(approxmst(1:n_snaps)%deg))
!
end
!
!-------------------------------------------------------------------------------
!
! Always check the return code of every netCDF function call. In
! this example program, wrapping netCDF calls with "call check()"
! makes sure that any return which is not equal to nf90_noerr (0)
! will print a netCDF error message and exit.
!
!
subroutine check_err(ret)
!
  use m_variables_gen
  use netcdf
!
  implicit none
!
  integer ret
!
  if (ret.ne.NF90_NOERR) then
    call sipr('Fatal. File I/O error in NetCDF operation. Returned error status is: ', ret)
    call spr(nf90_strerror(ret))
    call spr('By the following NetCDF library linked to CAMPARI:')
    call spr(nf90_inq_libvers())
! ' Reading from xtc-file (',netcdfinfile(t1:t2),')&
! &exited with an error (got ',ret,'). Sequence mismatch? Number of s&
! &napshots incorrect?'
!    ret2 = nf90_close(ncid)
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------------
!
!
subroutine read_nbl_nc()
!
  use gutenberg
  use m_variables_gen
  use m_mst_dumping
  use netcdf
!
  implicit none
!
  integer ncid,t1,t2,fndds,dimlen,nframes,i,ret,ilast,istart
  logical exists
  character(len = 100) ucstr, trystr, dumpfile
  integer, ALLOCATABLE :: vdeg(:), vsnp(:)
  real, ALLOCATABLE :: vdist(:)
  integer :: cnc_ids(10)
!
  dumpfile = 'MST_DUMPLING.nc'
  inquire(file=dumpfile,exist=exists)
!
  if (exists.EQV..false.) then
    call spr('Fatal. Cannot open spec.d input file for reading &
 &from NetCDF in setup_netcdftraj().')
    call fexit()
  end if
! open
  call check_err( nf90_open(dumpfile, NF90_NOwrite, ncid) )
!
! find the necessary dimensions: three are required
  fndds = 0
  nframes = 0
  do i=1,NF90_MAX_DIMS
    ret = nf90_inquire_dimension(ncid,i,trystr,dimlen)
    if (ret.eq.NF90_NOERR) then
      ucstr(1:15) = trystr(1:15)
      call toupper(ucstr(1:15))
      if (ucstr(1:10).eq.'FRAMEPAIRS') then
        if (fndds.eq.1) then
          call spr('Warning. Ambiguous dimensions in NetCDF file. Encountered '//trim(ucstr(1:10))//' twice but keeping&
          &only first.')
        else
          nframes = dimlen
          fndds = 1
          cnc_ids(1) = i
        end if
      end if
    else if (ret.eq.NF90_EBADDIM) then
!     do nothing
    else ! get us out of here
      call check_err( nf90_inquire_dimension(ncid,i,trystr,dimlen) )
    end if
  end do
!
  if (nframes.lt.1) then
      call spr('Fatal. NetCDF-file ('//trim(dumpfile)//') has no neighbor &
      &data (empty containers).')
    call fexit()
  end if
!
! now find the necessary variables, only coordinates are required
  ret = nf90_inq_varid(ncid,"snapshots",cnc_ids(2))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
      call spr('Fatal. Variable "snapshots" not found in NetCDF-file ('//trim(dumpfile)//'&
      ). Use NetCDFs ncdump utility to check file.')
    call fexit()
  else
    call check_err( nf90_inq_varid(ncid,"snapshots",cnc_ids(2)) )
  end if
  ret = nf90_inq_varid(ncid,"neighbors",cnc_ids(3))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    call spr('Fatal. Variable "neighbors" not found in NetCDF-file ('//trim(dumpfile)//'&
    ). Use NetCDFs ncdump utility to check file.')
    call fexit()
  else
    call check_err( nf90_inq_varid(ncid,"neighbors",cnc_ids(3)) )
  end if
  ret = nf90_inq_varid(ncid,"distances",cnc_ids(4))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    call spr('Fatal. Variable "distances" not found in NetCDF-file ('//trim(dumpfile)//'&
    ). Use NetCDFs ncdump utility to check file.')
    call fexit()
  else
    call check_err( nf90_inq_varid(ncid,"distances",cnc_ids(4)) )
  end if
!
  allocate(vdeg(nframes))
  allocate(vdist(nframes))
  allocate(vsnp(nframes))
!
! for some reason nf90_get_var chokes if the increment (count) is very large
  istart = 1
  ilast = min(nframes,10000)
  do while (istart.le.nframes)
    call check_err( nf90_get_var(ncid, cnc_ids(2), vsnp(istart:ilast), &
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    call check_err( nf90_get_var(ncid, cnc_ids(3), vdeg(istart:ilast),&
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    call check_err( nf90_get_var(ncid, cnc_ids(4), vdist(istart:ilast),&
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    istart = ilast + 1
    ilast = min(nframes,istart+10000)
  end do
!
  allocate(approxmst(n_snaps))
  approxmst(:)%deg = 0
  i = 1
  ilast = 0
  do istart=vsnp(1),vsnp(nframes)
    if (istart.ne.vsnp(i)) cycle
    do while (vsnp(i).eq.istart)
      if (i.eq.nframes) exit
      if (vsnp(i+1).ne.vsnp(i)) exit
      i = i + 1
    end do
    approxmst(vsnp(i))%deg = i - ilast
    approxmst(vsnp(i))%alsz = approxmst(vsnp(i))%deg
    if (approxmst(vsnp(i))%deg.gt.0) then
      allocate(approxmst(vsnp(i))%adj(approxmst(vsnp(i))%deg))
      allocate(approxmst(vsnp(i))%dist(approxmst(vsnp(i))%deg))
      approxmst(vsnp(i))%adj(:) = vdeg(ilast+1:i)
      approxmst(vsnp(i))%dist(:) = real(vdist(ilast+1:i),KIND=4)
    else
      approxmst(vsnp(i))%deg = 0
      approxmst(vsnp(i))%alsz = 2
      allocate(approxmst(vsnp(i))%adj(approxmst(vsnp(i))%alsz))
      allocate(approxmst(vsnp(i))%dist(approxmst(vsnp(i))%alsz))
    end if
    if (i.eq.nframes) exit
    ilast = i
    i = i  +1
  end do
  call sipr('Read a total number of snapshots of: ', sum(approxmst(1:n_snaps)%deg))
  deallocate(vsnp)
  ! deallocate(vnbs)
  deallocate(vdist)
!
end
!
!-----------------------------------------------------------------------
!
! this is a standard string operation
! note that the intrinsic ichar() uses numerical character values
! which are NOT necessarily ASCII codes (for that iachar())
!
subroutine toupper(string)
!
  implicit none
!
  integer i,leng,nasc
  character(*) string
!
  leng=len(string)
  do i=1,leng
    nasc=ichar(string(i:i))
!   the numerical range of standard lower-case letters
    if ((nasc.ge.97).AND.(nasc.le.122)) then
      string(i:i)=char(nasc-32)
    end if
  end do
!
end

#endif
