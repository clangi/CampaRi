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
module torsn
!
  integer NSEG,MAXFYMAPS,MAXTORRES
  parameter (NSEG=8)
  parameter (MAXFYMAPS=1000)
  parameter (MAXTORRES=12)
!
  integer ntorsn,ntorpuck,ndyntorsn
!
  RTYPE, ALLOCATABLE::  z_alpha(:),z_beta(:),z_hist(:,:,:)
  RTYPE, ALLOCATABLE::  z_2dhist(:,:,:),genfymap(:,:),z_hist2(:,:,:)
  RTYPE, ALLOCATABLE::  specfymap(:,:,:),molfymap(:,:,:)
  RTYPE, ALLOCATABLE::  torlc_coeff(:,:),torlc_coeff2(:,:)
  RTYPE, ALLOCATABLE::  torlc_hist(:,:),lct_weight(:),jc_res(:)
  RTYPE, ALLOCATABLE::  dih_hist(:,:),imp_mids(:),imp_hist(:,:)
  RTYPE, ALLOCATABLE::  bln_hist(:,:),ang_hist(:,:),ang_mids(:)
  RTYPE, ALLOCATABLE::  bln_mids(:),curtvec(:),seg_perrs(:,:,:)
  RTYPE fyres,torlc_params(2),intres(4)
!
  integer lof(NSEG),loy(NSEG),hif(NSEG),hiy(NSEG),segfymap(36,36)
  integer, ALLOCATABLE:: jccnt(:),sgprscnt(:),torzmlst(:)
  integer molrama(MAXFYMAPS),specrama(MAXFYMAPS),seg_alsz,nzmlst
  integer ntorlcs,torlccalc,torout,toroutmode,segcalc,angcalc,nrmolrama,rama_alsz(2)
  integer nrspecrama,maxseglen,torbins,intcalc,intszs(4,2)
  character(MAXSTRLEN) torlcfile,bbsegfile
  logical do_ints(4)
!
! covariance stuff
  integer, ALLOCATABLE::  ncovaravg(:),itrcv(:)
  integer covcalc,covmode
  logical, ALLOCATABLE:: trcvopen(:)
!
end module torsn
!
