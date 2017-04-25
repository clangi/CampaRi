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
! module for average polymer properties
!
module polyavg
!
  integer DLSTBINS,TPBINS
  parameter (DLSTBINS=100)
  parameter (TPBINS=100)
!
 type t_polynew
    RTYPE, ALLOCATABLE:: rg(:),rg2(:),rgt1(:),rgt2(:),rgt12(:),vt(:),ree(:),rhr(:)
    RTYPE, ALLOCATABLE:: rgev(:,:),cost(:),cosa(:,:),rh(:,:),psc(:,:),dep(:,:)
    integer, ALLOCATABLE:: map(:),nm(:,:),nmr(:)
  end type t_polynew
!
  type(t_polynew) polynew

  RTYPE, ALLOCATABLE:: rgavg(:),rg2avg(:),rhavg(:,:),reteavg(:),persavg(:,:)
  RTYPE, ALLOCATABLE:: rgterm1avg(:),rgterm2avg(:),rgterm12avg(:)
  RTYPE, ALLOCATABLE:: dlavg(:),asphavg(:),acylavg(:),dlstavg(:)
  RTYPE, ALLOCATABLE:: rgpcavg(:,:),rdhist(:,:,:),retehist(:,:)
  RTYPE, ALLOCATABLE:: turns_rs(:,:),pscavg(:,:),densproavg(:,:)
  RTYPE, ALLOCATABLE:: vtavg(:),rghist(:,:),rhravg(:)
  RTYPE rg_binsz,qv_res,tp_exp,dlst_binsz,tp_binsz,ndepravg
  integer polcalc,polout,sctcalc,nqv,holescalc,rhcalc
  integer maxrgbins
  integer, ALLOCATABLE:: nsctavg(:),npolavg(:),npersavg(:)
  integer, ALLOCATABLE:: nrhavg(:)
  logical inst_pers
!
end module polyavg
!

