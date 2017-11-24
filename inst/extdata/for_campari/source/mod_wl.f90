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
! MAIN AUTHOR:   Nicholas Lyle                                             !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module wl
!
  implicit none
!
  type t_wl
    integer wl_mode    ! 1 = generate g(E), 2 = generate PMF(RC), 3 = generate PMF(E) 
    integer dimensionality               ! 1 or 2 
    RTYPE g_min,g_min2d(2)               !minimum value of the DOS bin
    RTYPE g_max,g_max2d(2)               !maximum value of the DOS bin
    RTYPE g_binsz,g_binsz2d(2)           !DOS bin size
    RTYPE fval         !WL simulation convergence value initial value (LOG(f))
    RTYPE minv,minv2d(2)                 !current, left-most boundary of histograms
    integer gh_flatcheck_freq            !how often to check if reference energy histogram is flat
    integer buffer     !steps to delay f-updates by
    logical freeze     !whether to restrict histograms to bins seen at point buffer
    integer wl_mol(2)  ! over which molecule(s) to traverse p(RC)
    integer wl_rc(2)   ! what molecular reaction coordinates to use
    RTYPE, ALLOCATABLE:: g(:),g2d(:,:)   !density of states (g is stored as ln(g)) as function of E or arbitrary coordinate
    integer, ALLOCATABLE:: gh(:),ghbu(:),gh2d(:,:),ghbu2d(:,:) !g visitation histogram. used as for flatness criterion
    integer, ALLOCATABLE:: ghtot(:),ghtot2d(:,:)               !g visitation histogram that is never reset
    RTYPE, ALLOCATABLE:: bctr(:),bctr2d1(:),bctr2d2(:)         !gh and g bin centers
    integer nbins,nbins2d(2)                                   !current total bins in histograms
    integer hvmode     !how many times to visit each energy bin in order to update the F value 
    integer exrule     !whether to grow histograms dynamically or not
    logical t1on       !activate 1/t stage of the algorithm
    integer stepnum    !total WL iterations (not = istep) differs because 1 chi move has several iterations
    integer stepnumbuf ! a buffer used in restarted runs
    integer flevel     ! number of "stages" (gh resets) before using 1/t algorithm
    integer hufreq     !Decorrelation: update histograms with this frequency
!   Various WL statistics
    integer overflow   !counts that went over max bin window
    integer underflow  !counts that went under bin window
    integer accepted   !Total WL moves accepted
    integer maxb,minb,maxb2d(2),minb2d(2)!current bin ranges to use for convergence checks
!   fileIO
    integer ightot     !g visitation histogram that is never reset
    integer ig         !converged g output
    integer igin       !initializer for g(e)
    logical use_ginitfile
    character(MAXSTRLEN) ginitfile !Initialize g(E) with these values instead of zeros
  end type t_wl
  type(t_wl):: wld
!
  logical do_wanglandau,debug_wanglandau,finit_wanglandau,do_accelsim
!
! data structure for the weakly related accelerated molecular dynamics of McCammon and colleagues
  type t_hmjam
    integer nlst                    !number of boosted energy terms
    logical fill                    !whether to fill basins
    logical shave                   !whether to shave barriers
    logical, ALLOCATABLE:: isin(:)  !boost dimensions by energy term 
    integer, ALLOCATABLE:: lst(:)   !boost dimensions (temp)
    RTYPE, ALLOCATABLE:: alpha(:,:) !boost smoothifier
    RTYPE, ALLOCATABLE:: thresh(:,:)!boost thresholds
    RTYPE, ALLOCATABLE:: boosts(:,:)!boost energies
    RTYPE, ALLOCATABLE:: ca_f(:,:)  !temporary force array
    RTYPE, ALLOCATABLE:: t_ca_f(:,:)   !temporary force array in transpose
    RTYPE, ALLOCATABLE:: ca_f_tr(:,:)  !temporary force array #2
    RTYPE, ALLOCATABLE:: ca_f_bu(:,:)  !temporary force array #3
    RTYPE, ALLOCATABLE:: evec_thr(:,:) !temporay energy array for threads
    integer iwtsfile                !output file handle for frames with non-unity weights
    integer prtfrmwts               !whether to write frames requiring reweighting right away
    RTYPE threshwt                  !threshold weight for printing step number and associated weight
  end type t_hmjam
  type(t_hmjam):: hmjam
!
end module wl
!
