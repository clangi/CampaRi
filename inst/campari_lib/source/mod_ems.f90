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
module ems
!
  type t_emgrid
    integer dim(3),cnt,idx(3,2),shf(3),dimc(3)
    integer, ALLOCATABLE:: blklms(:,:,:)
    integer, ALLOCATABLE:: xmap(:),ymap(:),zmap(:),xmap2(:),ymap2(:),zmap2(:)
    RTYPE origin(3),deltas(3)
    RTYPE, ALLOCATABLE:: mass(:,:,:),avgmass(:,:,:),dens(:,:,:),deltadens(:,:,:)
    RTYPE, ALLOCATABLE:: mass2(:,:,:),optdens(:,:,:),rslc(:),rlin(:,:),rblk(:,:,:)
    logical, ALLOCATABLE:: lslc(:,:),llin(:,:,:),lblk(:,:,:,:)
  end type t_emgrid
  type(t_emgrid) emgrid,emingrid,emcoarsegrid
!
  character(MAXSTRLEN) emmapfile,emdensunitstr,emspatunitstr,emoriunitstr,empropfile
  integer emmapsz(3),emnetcdf_ids(100),emcalc,emsplor,empotmode,emmcacc,emheuristic,empotprop
  RTYPE embgdensity,emthreshdensity,emtotalmass,emavgdensity,emflatval,emthreshdensity2
  RTYPE fdlts(3),emtruncate,emoptbg,emiw,embuf
  logical curmassm
  integer, ALLOCATABLE:: emgfls(:,:)
  RTYPE, ALLOCATABLE:: embspld(:,:,:),embspl(:,:,:),emcustom(:,:)
!
end module ems

