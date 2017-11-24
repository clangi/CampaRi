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
! MAIN AUTHOR:   Hoang Tran                                                !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! phs       : target pH
! ionres    : simple list ionizable residues
! nionres   : number of ionizable residues (size of ionres)
! ionstate  : ionization state per ionres
! pkai      : "intrinsic" pKas per ionres
! titsite   : pointer array into atom numbers of titratable sites per ionres
! ionicstr  : ionic strength (M)
! invdhlen  : 1/debye-huckel screening length (1/Angst.)
!
module ionize
!
  integer nionres
  integer, ALLOCATABLE:: ionres(:),ionstate(:),titsite(:),prtionst(:)
  RTYPE phs,ionicstr,invdhlen
  RTYPE, ALLOCATABLE:: pkai(:)
!
end module ionize
!
