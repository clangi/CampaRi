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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ni          :  pointer array into atomic arrays for backbone amid(n)e nitrogen per residue
! cai         :  pointer array into atomic arrays for backbone carbon-alpha per residue
! ci          :  pointer array into atomic arrays for backbone carbon(x)yl carbon per residue
! oi          :  pointer array into atomic arrays for backbone carbon(x)yl oxygen per residue
! hni         :  pointer array into atomic arrays for backbone amid(n)e hydrogen per residue
! nuci        :  pointer array into atomic arrays for the six minimal nucleic acid backbone atoms (O3-P-O5-C5-C4-C3)
! at%bb       :  pointer array into atomic arrays per backbone atom per residue
! at%nbb      :  number of backbone atoms per residue
! at%sc       :  pointer array into atomic arrays per sidechain atom per residue
! at%nsc      :  number of sidechain atoms per residue
! at%pol      :  pointer array into atomic arrays per polar atom per residue
! at%npol     :  number of polar atoms per residue
! at%ndpgrps  :  number of dipole (charge) groups (=CG) per residue
!                at%dpgrp%nats    :  number of atoms per CG per residue
!                at%dpgrp%cgn     :  pointer array into cglst per CG per residue
!                at%dpgrp%nc      :  net (integer) charge per CG per residue (exception handling: set to next lowest
!                                    integer for groups with negative fractional charge and to next highest for pos. charge)
!                at%dpgrp%qnm     :  sum of absolute value of all charges per CG per residue
!                at%dpgrp%tc      :  net (floating point) charge per CG per residue
!                at%dpgrp%ats     :  pointer array into atomic arrays per CG per residue
! at%nfosgrps :  number of solvation groups (=SG) per residue
!                at%fosgrp%nats   :  number of atoms per SG per residue
!                at%fosgrp%val    :  free energy of solvation (at ref. T), enthalpy of solvation,
!                                    and heat capacity change per SG per residue
!                at%fosgrp%wts    :  weight factor per atom per SG per residue
!                at%fosgrp%ats    :  pointer array into atomic arrays per SG per residue
!
module polypep
!
  type t_dpgrp
    integer, ALLOCATABLE:: ats(:)
    integer nats,nc,cgn
    RTYPE qnm,tc,tsavq
  end type t_dpgrp
  type t_fosgrp
    integer, ALLOCATABLE:: ats(:)
    RTYPE, ALLOCATABLE:: wts(:)
    integer nats
    RTYPE val(3)
  end type t_fosgrp
  type t_at
    integer, ALLOCATABLE:: bb(:),sc(:),pol(:)
    integer nbb,nsc,npol,ndpgrps,nfosgrps,na,nncgrps
    type (t_dpgrp), ALLOCATABLE:: dpgrp(:)
    type (t_fosgrp), ALLOCATABLE:: fosgrp(:)
  end type t_at
!
  type(t_at), ALLOCATABLE:: at(:)
  integer, ALLOCATABLE:: ni(:),cai(:),ci(:),oi(:),hni(:)
  integer, ALLOCATABLE:: nuci(:,:)
!
end module polypep
!
!
