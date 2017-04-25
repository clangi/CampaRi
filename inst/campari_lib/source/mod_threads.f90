!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 2.0                                                           !
!                                                                          !
!    Copyright (C) 2014, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger                !
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
! MAIN AUTHOR:   Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module threads
!
  integer MAXTHREADS
  parameter (MAXTHREADS=128)
!
  type t_threads
    integer maxn                      ! (maximum) number of threads to use
    integer verbosity                 ! how verbose to be in output
    integer dlbfreq                   ! frequency to re-enable DLB
    integer dlbscale                  ! scale to perform slower but potentially more accurate DLB
    integer dlblen                    ! maximum length of a DLB segment
    logical test                      ! whether to run a test of threaded routines instead of a simulation
    logical suppflag(5)               ! flags whether to suppress internal parallelization of larger entities
    integer subnrs(5)                 ! team sizes for spawning of subteams
    RTYPE rate                        ! clock ticks per second for system clock
  end type t_threads
!
  type t_thr_rsp_nbl
    integer, ALLOCATABLE:: sr(:,:),tr(:,:),lr(:,:),tmpl(:)
    integer sralsz,tralsz,lralsz,srnrs,trnrs,lrnrs
  end type t_thr_rsp_nbl
!
  type(t_threads):: thrdat
!
  integer ithread                        ! file handle for logging
!
  integer, ALLOCATABLE:: thr_limits(:,:)   ! operation limits for various functionalities
  RTYPE, ALLOCATABLE:: thr_ca_f(:,:,:)     ! store Cartesian force per atom in transpose 
  RTYPE, ALLOCATABLE:: thr_ca_f_tr(:,:,:)  ! store mid-range Cartesian force per atom in transpose
  RTYPE, ALLOCATABLE:: thr_svte(:,:)       ! for ABSINTH, we need to replicate some more arrays (see also thr_sisa in forces)
  RTYPE, ALLOCATABLE:: thr_sum_scrcbs(:,:)  ! ditto
  RTYPE, ALLOCATABLE:: thr_sum_scrcbs_tr(:,:)  ! ditto
  RTYPE, ALLOCATABLE:: thr_cart2f_hlp(:,:,:)   ! temporary array for TMD
  type(t_thr_rsp_nbl), ALLOCATABLE:: thr_rsp_nbl(:)  ! temporary array for cutoffs in MC
  integer(KIND=8), ALLOCATABLE:: thr_timings(:,:)    ! timing information for DLB
  integer, ALLOCATABLE:: thr_dlb(:,:)                ! flag to store which load limits are currently adjustable + counter
  integer, ALLOCATABLE:: mlg_limits(:,:,:)           ! for parallelization within molecules
  integer, ALLOCATABLE:: mlg_ixs(:,:,:,:)            ! ditto
  integer, ALLOCATABLE:: thr_mlgix(:)                ! ditto (reverse array)
  integer nmlgs                                      ! number of parallelized molecules
  integer, ALLOCATABLE:: ccg_limits(:,:,:)           ! for parallelization within holo-constraints groups
  integer nccgs                                      ! number of parallelized holo-groups
  RTYPE, ALLOCATABLE:: thr_rutil(:,:)                ! utility vector for reduce ops (permanently allocated) 
  integer, ALLOCATABLE:: thr_hlper(:,:)              ! helper array for unknown sizes (realloc)
  RTYPE, ALLOCATABLE:: thr_rhlper(:,:)               ! ditto
  logical, ALLOCATABLE:: thr_lhlper(:,:)             ! ditto
!
!  integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!
end module threads
!

