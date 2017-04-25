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
! MAIN AUTHOR:   Albert Mao                                                !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module grandensembles
    implicit none
!
!    ! We either use absolute chemical potentials or prescribe a set
!    ! of equilibrium concentrations and use excess chemical potentials
    integer gc_mode
!    
!    ! This type and the subroutines and functions from
!    ! particlefluctuation.f are supposed to be grouped together into
!    ! one module that encapsulates the concept of an integer set.
!    ! Unfortunately, in order to fit within the framework Andreas
!    ! designed for modules in CAMPARI, the derived type definition
!    ! and the methods that operate on it must be separated because
!    ! modules cannot use other modules.
    type integer_set
      integer universesize, nummembers
      integer, ALLOCATABLE:: indexof(:), array(:)
    end type integer_set
!    
    character(MAXSTRLEN):: particleflucfile
    type(integer_set) fluctypes
    type(integer_set) ispresent
!    ! Second index is 1 for present molecules, 2 for nonpresent
!    ! molecules:
    type(integer_set), ALLOCATABLE:: typesets(:,:)
    RTYPE, ALLOCATABLE:: chempot(:),eqnum(:)
    RTYPE, ALLOCATABLE:: thermalvolume(:) ! not currently used
    logical sgc_report
!    
!    ! Variables used for particle number analysis, which is
!    ! performed in particleflucstat.
    integer particlenumcalc
    integer, ALLOCATABLE:: numberhistogram(:,:)
!
!   ! Warning counters
    integer gc_wrncnt(5),gc_wrnlmt(5)
!
end module grandensembles
!
