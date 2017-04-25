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
module dssps
!
  integer MAXHBS
  parameter (MAXHBS=10)
!
  type t_hbl
    integer nahbs,ndhbs
    integer ahb(MAXHBS),dhb(MAXHBS)
    RTYPE enahb(MAXHBS),endhb(MAXHBS)
    RTYPE denahb(6,MAXHBS),dendhb(6,MAXHBS)
  end type t_hbl
!
  type t_bdgtb
    integer nrs
    integer, ALLOCATABLE:: lims(:,:)
    integer, ALLOCATABLE:: ty(:)
    RTYPE, ALLOCATABLE:: phbe(:)
  end type t_bdgtb
!
  type (t_hbl), ALLOCATABLE:: hbl(:)
  type (t_bdgtb) bdgtb
!
  integer, ALLOCATABLE:: hb_mat(:,:),pepmol(:),molpep(:,:),dssp_ass(:,:)
  integer, ALLOCATABLE:: peprs(:),rspep(:),pepmolsz(:),dssp_cnt(:)
  RTYPE, ALLOCATABLE:: dssp_2dhist(:,:),dssp_hists(:,:),dssp_perrs(:,:,:)
  RTYPE, ALLOCATABLE:: whsc_der(:,:,:),wesc_der(:,:,:)
  integer pep_sz,dsspcalc,par_DSSP2(10),npepmol
  RTYPE par_DSSP(10)
  logical inst_dssp
  character dssp_map(8,11)
  data dssp_map / ' ',' ',' ',' ',' ',' ',' ',' ',&
                 &'E','3','4','5',' ','-',' ',' ',&
                 &'E','<','<','<',' ','+',' ',' ',&
                 &'E','>','>','>',' ',' ',' ',' ',&
                 &'E','X','X','X',' ',' ',' ',' ',&
                 &'E',' ',' ',' ',' ',' ',' ',' ',&
                 &'H',' ',' ',' ',' ',' ',' ',' ',&
                 &'G',' ',' ',' ',' ',' ',' ',' ',&
                 &'I',' ',' ',' ',' ',' ',' ',' ',&
                 &'T',' ',' ',' ',' ',' ',' ',' ',&
                 &'S',' ',' ',' ','S',' ',' ',' '/
!
end module dssps
!
