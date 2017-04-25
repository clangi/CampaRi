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
! n12,i12,... : topology arrays and variables: number and indeces of of 1-X-partners
! x,y,z       : Cartesian coordinate arrays
! xref,...    : Cartesian coordinate backup arrays
! n           : total atom number
! mass        : mass per atom
! ivms        : inverse mass per atom
! atmres      : residue-# pointer per atom
! attyp       : LJ atom type number per atom
! atnam       : LJ atom type name
! b_type      : biotype number per atom
! atr         : radius for solvation model (should be Pauling-like, not MMFF-like) per atom
! atvol       : naked atomic volume for solvation model per atom
! atbvol      : solvation shell volume for solvation model per atom
! atsav       : SAV for solvation model per atom
! atsavred
! atsavprm
! atsavmaxfr
! svte
! svbu
! atstv       : various solvation model parameters per atom
! atq         : partial charge per atom
! scrq        : partial charge screening factor per atom
! eqatm       : rigid topology equivalence atom pointer per atom (usually self)
! which_fosg  : solvation group index per atom (within residue -> atmres)
! which_dipg  : dipole/charge group index per atom (within residue -> atmres)
!
module atoms
!
  integer MAXVLNC
  parameter (MAXVLNC=4)
!
  integer n
  integer, ALLOCATABLE:: attyp(:),b_type(:)
  integer, ALLOCATABLE:: n12(:),i12(:,:),n13(:),i13(:,:)
  integer, ALLOCATABLE:: n14(:),i14(:,:),atmres(:),eqatm(:)
  integer, ALLOCATABLE:: which_fosg(:),which_dpg(:)
  character(3), ALLOCATABLE:: atnam(:)
!
  RTYPE, ALLOCATABLE:: x(:),y(:),z(:),xref(:),yref(:),zref(:)
  RTYPE, ALLOCATABLE:: atr(:),atvol(:),atsav(:),atstv(:),mass(:),ivms(:)
  RTYPE, ALLOCATABLE:: svte(:),svbu(:),atbvol(:),atsavmaxfr(:)
  RTYPE, ALLOCATABLE:: atsavred(:),atsavprm(:,:),atq(:),scrq(:)
!
end module atoms
!
