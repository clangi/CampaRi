
#ifdef LINK_NETCDF

!------------------------------------------------------------------------------
!
! #ifdef LINK_NETCDF
!
subroutine dump_nbl_nc()
!
  use gutenberg
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
  use gutenberg
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
