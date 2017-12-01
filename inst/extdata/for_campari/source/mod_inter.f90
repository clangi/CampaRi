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
module inter
!
  type t_fudge
    RTYPE, ALLOCATABLE:: rsin(:)
    RTYPE, ALLOCATABLE:: rsnb(:)
    RTYPE, ALLOCATABLE:: rsin_ljs(:)
    RTYPE, ALLOCATABLE:: rsnb_ljs(:)
    RTYPE, ALLOCATABLE:: rsin_lje(:)
    RTYPE, ALLOCATABLE:: rsnb_lje(:)
    RTYPE, ALLOCATABLE:: elin(:)
    RTYPE, ALLOCATABLE:: elnb(:)
  end type t_fudge
!
  type t_iaa
    integer, ALLOCATABLE:: atin(:,:)
    integer, ALLOCATABLE:: atnb(:,:)
    integer, ALLOCATABLE:: polin(:,:)
    integer, ALLOCATABLE:: polnb(:,:)
    integer, ALLOCATABLE:: expolin(:,:)
    integer, ALLOCATABLE:: expolnb(:,:)
    integer, ALLOCATABLE:: bl(:,:)
    integer, ALLOCATABLE:: ba(:,:)
    integer, ALLOCATABLE:: di(:,:)
    integer, ALLOCATABLE:: impt(:,:)
    integer, ALLOCATABLE:: cm(:,:)
    RTYPE, ALLOCATABLE:: par_bl(:,:)
    RTYPE, ALLOCATABLE:: par_ba(:,:)
    RTYPE, ALLOCATABLE:: par_di(:,:)
    RTYPE, ALLOCATABLE:: par_impt(:,:)
    integer, ALLOCATABLE:: typ_bl(:)
    integer, ALLOCATABLE:: typ_ba(:)
    integer, ALLOCATABLE:: typ_di(:)
    integer, ALLOCATABLE:: typ_impt(:)
    integer, ALLOCATABLE:: typ_cm(:)
  end type t_iaa
!
  type(t_fudge), ALLOCATABLE :: fudge(:)
  type(t_iaa), ALLOCATABLE :: iaa(:)
  integer, ALLOCATABLE:: nrsintra(:),nrsnb(:)
  integer, ALLOCATABLE:: nrpolintra(:),nrpolnb(:)
  integer, ALLOCATABLE:: nrexpolin(:),nrexpolnb(:),nrsbaeff(:)
  integer, ALLOCATABLE:: nrsbl(:),nrsba(:),nrsdi(:),nrsdieff(:),nrscmeff(:)
  integer, ALLOCATABLE:: nrsimpt(:),nrsimpteff(:),nrsbleff(:),nrscm(:)
  RTYPE fudge_st_14,fudge_el_14
  integer eps_type,mode_14,nbsr_model,elec_model
  logical use_14,ia_report
!
end module inter
!
