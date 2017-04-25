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
module forces
!
  integer MAXNB
  parameter (MAXNB=300)
!
  integer n_constraints,dyn_integrator,dyn_integrator_ops(10)
  RTYPE, ALLOCATABLE:: cart_f(:,:),sav_dr(:,:),scrq_dr(:,:)
  RTYPE, ALLOCATABLE:: cart_f_tr(:,:),cart_f_lr(:,:),tmpcf(:,:,:,:)
  RTYPE, ALLOCATABLE:: cart_a(:,:),cart_v(:,:)
  RTYPE, ALLOCATABLE:: cart_ldp(:,:,:),sum_scrcbs(:),sum_scrcbs_tr(:),sum_scrcbs_lr(:)
  RTYPE bnd_f(3,3),bnd_v,bnd_avgT,bnd_pV,bnd_fr
  integer fo_wrncnt(14),fo_wrnlmt(14)
  logical grad_check
  logical, ALLOCATABLE:: cart_frz(:,:)
  type t_dc_di
    integer maxntor
    integer, ALLOCATABLE:: nattor(:),recurs(:,:)
    RTYPE, ALLOCATABLE:: ldp(:,:),incr(:),f(:),im(:),olddat(:,:),v(:),avgT(:),valrecurs(:,:)
    LOGICAL, ALLOCATABLE:: frz(:),align(:)
  end type t_dc_di
  type(t_dc_di), ALLOCATABLE:: dc_di(:)
  type t_sisa
    integer nix,alsz 
    integer, ALLOCATABLE:: ix(:)
    RTYPE, ALLOCATABLE:: dr(:,:),qr(:,:)
  end type t_sisa
  type(t_sisa), ALLOCATABLE:: sisa(:),thr_sisa(:,:)
!
end module forces
!
