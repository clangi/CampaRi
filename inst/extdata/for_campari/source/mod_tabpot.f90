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
! MAIN AUTHOR:   Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module tabpot
!
  type t_tbp
    integer mode,all_lst,nlst,nlstgrp,nos,atms,bins,tangmode
    integer, ALLOCATABLE:: rsmat(:,:),rsvec(:)
    integer, ALLOCATABLE:: lst(:,:)
    logical, ALLOCATABLE:: lty(:)
    RTYPE, ALLOCATABLE:: pot(:,:),dis(:),tang(:,:)
    RTYPE res,ipprm(10)
  end type t_tbp
!
  type (t_tbp) tbp
!
  character(MAXSTRLEN) tabcodefile,tabpotfile,tabtangfile
!
  logical tabul_report
!
  RTYPE, ALLOCATABLE:: terfcoverr(:,:) ! tabulated version of erfc(a*sqrt(r2))/sqrt(r2)
  RTYPE p_terfcoverr(4)                ! parameters (bin size, upper bound, max error, max gradient error) for terfcoverr
  integer i_terfcoverr                 ! number of bins for terfcoverr
!
end module tabpot
!
