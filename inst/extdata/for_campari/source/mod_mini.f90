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
! CONTRIBUTIONS: Nicholas Lyle, Adam Steffen                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module mini
!  
  implicit none
!  
  integer mini_mode,mini_mem,mini_sc_sdsteps
  RTYPE mini_econv,mini_stepsize,mini_xyzstep,mini_rotstep,mini_intstep,mini_uphill
  RTYPE mini_sc_tbath,mini_sc_heat
!
  type t_minivs
    RTYPE, ALLOCATABLE:: cgrv(:),pgrv(:),dposv(:),peakg(:),aheadg(:),cgn(:),avstp(:),pcgn(:)
    RTYPE, ALLOCATABLE:: b_a(:),b_q(:),b_r(:),b_s(:,:),b_y(:,:),mshfs(:,:)
    RTYPE grms,hlp(10)
  end type t_minivs
!
  type(t_minivs):: minivs
!  
end module mini
!       
