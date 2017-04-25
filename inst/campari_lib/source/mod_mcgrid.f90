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
module mcgrid
!
!
  type t_mcgrid
    integer dim(3),bu(3),reset,chkrad,chl,pnts,maxgpnbs,maxnbs,ntmpl
    RTYPE origin(3),deltas(3),diags(2,2,2),sz(3)
    integer, ALLOCATABLE:: nblist(:,:),nbnr(:),resgp(:),tgrpgp(:),gpnblist(:,:),gpnbnr(:),gpnbmd(:,:),tmpl(:,:)
    RTYPE, ALLOCATABLE:: maxr(:)
    integer, ALLOCATABLE:: b_gpnbnr(:),b_gpnblist(:,:),b_gpnbmd(:,:)
    logical report
  end type t_mcgrid
  type(t_mcgrid) grid
!
end module mcgrid
!
