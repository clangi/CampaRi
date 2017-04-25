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
module molecule
!
  integer, ALLOCATABLE:: atmol(:,:),rsmol(:,:),moltermid(:,:)
  integer, ALLOCATABLE:: moltyp(:,:),moltypid(:),molangr(:,:)
  integer, ALLOCATABLE:: ntormol(:),an_grp_mol(:),molfrzidx(:),solutes(:),othidxmol(:)
  integer nmol,nmoltyp,totrbd,nangrps,longestmol,nsolutes
  RTYPE, ALLOCATABLE:: com(:,:),rgpcs(:,:,:),rgevs(:,:),rgv(:)
  RTYPE, ALLOCATABLE:: comm(:,:),rgpcsm(:,:,:),rgevsm(:,:),rgvm(:)
  RTYPE, ALLOCATABLE:: rgpcsref(:,:,:),rgevsref(:,:),commref(:,:)
  RTYPE, ALLOCATABLE:: rgvref(:),comref(:,:),molvol(:,:)
  RTYPE, ALLOCATABLE:: molmass(:),molcontlen(:),movingmass(:,:)
  RTYPE rstime,rstime2,rstime3,rstime4
  logical, ALLOCATABLE:: do_tors(:),do_pol(:),do_pers(:),molischg(:),is_solvent(:)
  character(MAXSTRLEN) angrpfile
!
end module molecule
!
