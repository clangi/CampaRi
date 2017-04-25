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
module paircorr
!
  integer MAXPCBINS
  parameter (MAXPCBINS=2500)
!
  type t_gpc
    integer mode,nos,atms,nlst,all_lst,all_atm
    logical, ALLOCATABLE:: lty(:)
    integer(8), ALLOCATABLE:: navg(:)
    integer, ALLOCATABLE:: lst(:,:)
    integer, ALLOCATABLE:: mat(:,:),rsmat(:,:)
    RTYPE, ALLOCATABLE:: pc(:,:)
  end type t_gpc
!
  type t_ampc
    RTYPE, ALLOCATABLE:: bb(:,:),ss(:,:),bs(:,:)
    logical do_ss,do_bb,do_bs
  end type t_ampc
!
  type (t_ampc) ampc
  type (t_gpc) gpc
  integer pccalc
  integer(8), ALLOCATABLE:: npcavg(:,:)
  logical gpc_report,do_amidepc
  character(MAXSTRLEN) pccodefile
!
  RTYPE, ALLOCATABLE:: rbc_pc(:,:,:),voli_pc(:)
  RTYPE pc_binsz
  integer inst_gpc   ! whether to write instantaneous distances 
  integer n_inst_gpc ! the total number of unique distances requested
  integer n_pccalls  ! total number of times do_general_pc has been called
!
end module paircorr
!
