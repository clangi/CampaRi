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
module shakeetal
!
  type t_constraints 
    RTYPE, ALLOCATABLE:: eqd2(:),eqa3(:),eqi(:),amat(:,:),ccoeff(:,:)
    integer, ALLOCATABLE:: idx(:,:),idx3(:,:),uidx(:),mapidx(:,:),cidx(:,:),ncis(:)
    integer, ALLOCATABLE:: idx4(:,:)
    integer nr,nr3,nr4,nats,maxcpl,itercnt,iters,lorder
    logical isrigid,hasmassless
    logical, ALLOCATABLE:: massless(:)
  end type t_constraints
  type(t_constraints), ALLOCATABLE:: constraints(:)
!
  RTYPE shake_tol,shake_atol
  integer cart_cons_grps,cart_cons_mode,shake_maxiter,cart_cons_method,cart_cons_source,lincs_order,lincs_iter
  integer, ALLOCATABLE:: settlelst1(:),settlelst2(:),settlelst3(:),nonsettlelst(:),settlelst4(:),settlelst5(:),settlelst6(:)
  integer settle_tip3ps,settle_spcs,settle_rest,shake_cnt,settle_tip4ps,settle_tip5ps,settle_tip4pes
  integer cs_wrnlmt(20),cs_wrncnt(20)
  character(MAXSTRLEN) cart_cons_file
  logical no_shake,add_settle,have_virtuals
!
end module shakeetal
!
