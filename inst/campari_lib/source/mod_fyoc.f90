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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! phi,psi,omega,chi,nucs     : value arrays for various degrees of freedom per residue 
! fline,yline,wline,...      : pointer arrays into Z-matrix for various degrees of freedom per residue
! fline2,yline2              : pointer arrays for dependent torsions (on phi,psi) per residue
! phish,psish                : angular relation arrays for dependent torsions (on phi,psi) per residue
! nchi,nnucs                 : net number arrays for degrees of freedom of identical type per residue
! cur_phi,cur_psi,...        : current proposed moves for various degrees of freedom
! fnr,ynr,...                : pointer arrays into molecular force arrays for various degrees of freedom per residue
! chiral                     : chirality flag per residue
! disulf                     : overlaps with crlk_idx (clean up eventually) in mod_sequen.f90
! lct,lctold                 : disregard
! 
module fyoc
!
  integer MAXCHI
  parameter (MAXCHI=6)
!
  integer, ALLOCATABLE:: disulf(:),chiral(:),nchi(:),chiline(:,:)
  integer, ALLOCATABLE:: fline(:),fline2(:),yline(:),yline2(:)
  integer, ALLOCATABLE:: wline(:),nucsline(:,:),nnucs(:),pucline(:)
  integer, ALLOCATABLE:: fnr(:),ynr(:),wnr(:),chinr(:,:),nucsnr(:,:)
  RTYPE, ALLOCATABLE:: phi(:),psi(:),omega(:),nucs(:,:),phish(:)
  RTYPE, ALLOCATABLE:: chi(:,:),lct(:),lctold(:),psish(:)
  RTYPE cur_phi,cur_psi,cur_chis(MAXCHI),cur_omega,cur_lct
  RTYPE cur_nucs(6),cur_pucks(25)
  logical harappa,cur_chiflag(MAXCHI),cur_nucflag(6)
!
end module fyoc
!
