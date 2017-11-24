!!--------------------------------------------------------------------------!
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

#include "macros.i"
!
module system
!
  type t_ensemble
    RTYPE avgT,avgK,avgK2,avgU,avgP,avgV,avgPT(3,3),avgVirT(3,3),avgKinT(3,3),avgR(10)
    RTYPE insT,insK,insK2,insU,insP,insV,insPT(3,3),insVirT(3,3),insKinT(3,3),insR(10)
    integer flag,avgcnt,sysfrz,boxdof,avgcnt2
  end type t_ensemble
!
  type t_thermostat
     integer flag
     RTYPE params(20)
     integer n_tgrps
     character(MAXSTRLEN) fnam
     integer, ALLOCATABLE:: molgrp(:)
     RTYPE, ALLOCATABLE:: grpT(:),grpdof(:)
  end type t_thermostat
!
  type t_manostat
     integer flag
     RTYPE params(20)
  end type t_manostat
!
! names
  character(MAXSTRLEN) basename
  integer bleng
!
! ensemble-related variables
  RTYPE invtemp,kelvin,extpress,invcompress,fric_ga
  type(t_ensemble):: ens
  type(t_thermostat):: tstat
  type(t_manostat):: pstat
!
! box variables
  RTYPE bnd_params(20)
  integer bnd_shape,bnd_type,bnd_wrncnt(2),bnd_wrnlmt(2)
!
! simulation dimensions
  integer nsim,nequil,dyn_mode,fycxyz,mc_compat_flag
  logical pdb_analyze,do_restart
  RTYPE dyn_dt
  logical use_dyn,dyn_report
  RTYPE sysmass(2),sysvol(3)
!
! misc.
  integer fudge_mass
  integer ua_model
!
end module system
!
