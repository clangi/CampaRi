!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 3.0                                                           !
!                                                                          !
!    Copyright (C) 2017, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  !
!                        Davide Garolini, Jiri Vymetal                     !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    CAMPARI is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Andreas Vitalis, Marco Bacci                              !
!                                                                          !
!--------------------------------------------------------------------------!
!     
!
#include "macros.i"
!
!
! #################################################
! ##                                             ##
! ## THIS IS THE MASTER SOURCE FILE FOR THE      ##
! ## LIGHT DATA MINING EXECUTION MODE OF CAMPARI ##
! ##                                             ##
! #################################################
!
!
program datasaw
!
  use iounit
  use, INTRINSIC:: ISO_C_BINDING
  use threads
  use ncdm
  use interfaces
  use clusters
  use mcsums
  use, INTRINSIC:: ISO_FORTRAN_ENV
!
  implicit none
!
  integer dt(8)
  integer(KIND=8) ee,ee2,ee1,t1,t2
  character (len=12) rc(3)
  integer rs,aone,tpi,atwo
#ifdef ENABLE_THREADS
  integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,tpn
#endif
  logical afalse
#ifdef ENABLE_MPI
  integer masterrank,ierr
#else
  integer st1,st2
  character(len=1,kind=c_char), dimension(MAXSTRLEN+1):: honam(0:MAXSTRLEN)
  character(MAXSTRLEN) hnam
  integer(c_int) hnerr
  integer(c_size_t) hnml
#endif
!
#ifdef ENABLE_MPI
  masterrank = 0
#endif
!
  afalse = .false.
  aone = 1
  atwo = 2
#ifdef ENABLE_MPI
!
! MPI: obviously some things need to be different ...
! this routine will have the master distribute the simulation information
! and run the setup procedures for all CPUs
  time_comm = 0.0
  call System_Clock(ee1)
  call MPI_STARTMC()
  call System_Clock(t1)
  time_comm = time_comm + 1.0*(t2 - t1)
  call System_Clock(t2)
#else
!
! when we start
  hnml = MAXSTRLEN
  honam(:) = ' '
  hnerr = handwrapped_gethostname(honam(0:(MAXSTRLEN-1)),hnml)
  do rs=1,MAXSTRLEN
    hnam(rs:rs) = honam(rs-1)
  end do
  call strlims(hnam,st1,st2)
  if (st2.gt.st1) then
    if (hnam(st2:st2).eq.C_NULL_CHAR) st2 = st2 - 1
  end if
  call Date_and_Time(rc(1),rc(2),rc(3),dt)
  write(*,*)
  call write_license_nc()
  write(*,*) 'Execution started on host ',hnam(st1:st2),' at ',rc(2)(1:2),':',rc(2)(3:4),&
 &' on ',rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
  call System_Clock(ee1)
  write(*,*)
!     
! all initializations for start up
  call initial()
!
! we need to read (not yet parse) the key (input) file
  call getkey()
  call init_campprng(aone)
!
! level 1 input reader
  write(ilog,*)
  write(ilog,*) '---   Now parsing keywords ...            ---'
  write(ilog,*)
  call parsekey(1)
! level 2 input reader
  call parsekey(2)
! level 3 input reader
  call parsekey(3)
  write(ilog,*)
  write(ilog,*) '---   ... finished parsing keywords.      ---'
  write(ilog,*)
#endif
!
! ------------------- NET-CDF-BASED DATA MINING ----------------------------------------------
!
#ifdef LINK_NETCDF
!
  call ncdm_sanity_checks_1()
  call ncdm_summary_1()
  call ncdm_makeio(aone)
!
  call System_Clock(t1)
!
! --- set ini ---
  call ncdm_initial()
!
! --- input files ---
  if (ncdm_isthere_framesfl.EQV..true.) then
    call ncdm_read_framesfl()
  else
    call ncdm_manage_whichframes(afalse)
  end if
  if (ncdm_isthere_cfile.EQV..true.) then 
    call ncdm_read_cfile()
  else
    call ncdm_manage_whichfeats(afalse)
  end if
! --- end of input files ---
!
  call ncdm_manage_ccollect(afalse) 
!
  if (ncdm_doas.EQV..true.) then
    call ncdm_convert_ascii()
  end if
  if (ncdm_donc.EQV..true.) then  !mutually exclusive with ncdm_doas (see ncdm_sanity_checks_1)
    call ncdm_read_ncfl()
  end if
  call System_Clock(t2)
  time_ph = time_ph + 1.0*(t2 - t1) ! hi-jacked
!
! --- data mining ---
!
  if ((ncdm_donc.EQV..true.).OR.(ncdm_isdmonas.EQV..true.)) then    !otherwise only conversion was requested
    call ncdm_fill_auxcludata()           !fill cludata related quantities
    call ncdm_sanity_checks_2()           !all the sanity checks that require the knowledge of cludata params
    call ncdm_summary_2()                 !all the summary that requires the knowledge of cludata params
    call build_sconnect(clagt_msm,aone)   !build connectivity map
    if (cdis_crit.le.4) then
      call ncdm_scale_periodicity_andsincos()
    end if
!
#ifdef ENABLE_THREADS
    allocate(thr_limits(110,thrdat%maxn))
    allocate(thr_timings(200,thrdat%maxn))
    allocate(thr_rutil(100,thrdat%maxn))
    thr_timings(:,:) = 0
    allocate(thr_dlb(20,3))
    call omp_set_num_threads(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpi,tpn)
    tpn = omp_get_num_threads()
    tpi = omp_get_thread_num() + 1
#else
    tpi = 0
#endif
    if (tpi.le.1) call System_Clock(t1)
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call flushopen()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
    call do_clustering(tpi)
    if (tpi.le.1) then
      call System_Clock(t2)
      time_analysis = time_analysis + 1.0*(t2 - t1)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
  end if !((ncdm_donc.EQV..true.).OR.(ncdm_isdmonas.EQV..true.)) then
!
#else
!
  write(ilog,*) "Fatal. CAMPARI's NetCDF data mining executable was compiled without NetCDF support. This is meaningless."
  call fexit()
!
#endif
! ----------------------------------------
!
! finalize I/O (close files)
  call ncdm_makeio(atwo)
!
! summarize performance information
  call System_Clock(ee2)
  ee = 1.0*(ee2 - ee1)
  write(ilog,*)
  write(ilog,*) 'Total CPU time elapsed [s]     : ',ee/thrdat%rate
  write(ilog,*) 'Fraction for analysis routines : ',100.0*time_analysis/ee,'%'
  write(ilog,*) 'Fraction for input file proc.  : ',100.0*time_ph/ee,'%'
#ifdef ENABLE_MPI
  write(ilog,*) 'Fraction for communication     : ',100.0*time_comm/ee,'%'
  write(ilog,*) 'Fraction of remainder          : ',100.0*(1.0 - time_analysis/ee - time_comm/ee - time_ph/ee),'%'
#else
  write(ilog,*) 'Fraction of remainder          : ',100.0*(1.0 - time_analysis/ee - time_ph/ee),'%'
#endif
  write(ilog,*)
!
  call Date_and_Time(rc(1),rc(2),rc(3),dt)
  write(ilog,*) 'Execution ended at ',rc(2)(1:2),':',rc(2)(3:4),&
 &' on ',rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
!
#ifdef ENABLE_MPI
  call makelogio(2)
#endif
!
! in parallel mode, need to shutdown MPI universe
#ifdef ENABLE_MPI
  call MPI_StopMC()
  call MPI_FINALIZE(ierr)
#endif
!
end program datasaw
!
!
! license function repeated for potential separability
!
!------------------------------------------------------------------------------
!
subroutine write_license_nc()
!
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
  write(*,*) '                  ---        Welcome to CAMPARI       ---'
  write(*,*)
  write(*,*) '                  ---       Released Version 3.0      ---'
  write(*,*)
  write(*,*) '    Copyright (C) 2017, CAMPARI Development Team'
  write(*,*)
  write(*,*) '    Websites: http://campari.sourceforge.net'
  write(*,*) '              http://sourceforge.net/projects/campari/'
  write(*,*) '              http://pappulab.wustl.edu/campari'
  write(*,*)
  write(*,*) '    CAMPARI is free software: you can redistribute it and/or modify'
  write(*,*) '    it under the terms of the GNU General Public License as published by'
  write(*,*) '    the Free Software Foundation, either version 3 of the License, or'
  write(*,*) '    (at your option) any later version.'
  write(*,*)
  write(*,*) '    CAMPARI is distributed in the hope that it will be useful,'
  write(*,*) '    but WITHOUT ANY WARRANTY; without even the implied warranty of'
  write(*,*) '    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
  write(*,*) '    GNU General Public License for more details. '
  write(*,*)
  write(*,*) '    You should have received a copy of the GNU General Public License' 
  write(*,*) '    along with CAMPARI. If not, see <http://www.gnu.org/licenses/>.'
  write(*,*)
  write(*,*)
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
  write(*,*)
!
end
!
!------------------------------------------------------------------------------
!

