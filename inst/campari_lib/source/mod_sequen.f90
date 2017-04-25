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
! nseq        : number of residues
! seqtyp      : integer code for residue type per residue
! seqpolty    : character code for polymer type per residue
! seqflag     : another residue class identifier independent of seqtyp and seqpolty
! resrad      : effective maximum radius per residue
! natres      : upper limit for number of atoms per residue
! molofrs     : pointer array into molecular arrays per residue
! refat       : pointer array for "central" atom into atomic arrays per residue
! resfy,...   : temporary list of residues with a particular set of degrees of freedom
! notors      : flag for rigid residues (usually ions and small molecules)
! chgflag     : flag whether a residue has at least one charge group with a non-zero net charge per residue
! netchg      : integer giving the (typically) full charge the residue is carrying: for the few allowed
!               exceptions with non-integer net charges on a residue, this is set to the next lowest integer
!               for negative charges, and the next highest integer for positive charges (!!!)
! *infile     : structural input files mostly for trajectory analysis
! n_crosslinks: number of chemical crosslinks
! crlk_idx    : crosslink pointer per residue
! crosslink%exclin      : exclusion list for short-range, steric interactions created by crosslink
! crosslink%nrsin       : corresponding length
! crosslink%exclpol     : exclusion list for short-range, Coulombic interactions created by crosslink
! crosslink%nrspol      : corresponding length
! crosslink%itstype     : chemical type identifier
! crosslink%rsnrs       : constituent residues by number
! crosslink%ljpars      : custom LJ parameters for short-range, steric interactions created by crosslink
! crosslink%cbpars      : custom point charge (product) parameters for short-range, Coulombic interactions created by crosslink
! crosslink%is14in      : corresponding 14-indicator (steric)
! crosslink%is14pol     : corresponding 14-indicator (Coulombic)
! crosslink%olks        : coupling set for crosslink
! crosslink%nolks       : corresponding length
!  
!
module sequen
!
  type t_crosslink
    integer, ALLOCATABLE:: olks(:),exclin(:,:),exclpol(:,:)
    integer nolks,itstype,rsnrs(2),nrsin,nrspol
    RTYPE, ALLOCATABLE:: ljpars(:,:),cbpars(:)
    logical, ALLOCATABLE:: is14in(:),is14pol(:)
  end type t_crosslink
!
  type(t_crosslink), ALLOCATABLE:: crosslink(:)
  integer nseq,nressolute,n_crosslinks,crlk_mode,globrandomize,globrdmatts
  integer, ALLOCATABLE:: seqtyp(:),natres(:),molofrs(:),seqflag(:)
  integer, ALLOCATABLE:: an_grp_res(:),respuc(:),crlk_idx(:)
  integer, ALLOCATABLE:: resw(:),resfy(:),netchg(:),reschi(:)
  integer, ALLOCATABLE:: rs_vec(:),rs_vec2(:),refat(:),resnuc(:),resnucpuc(:)
  RTYPE, ALLOCATABLE:: resrad(:),resvol(:,:)
  RTYPE mrrd,globrdmthresh
  logical globcyclic,fycinput,pdbinput,seq_report
  logical, ALLOCATABLE:: chgflag(:),notors(:)
  character(1), ALLOCATABLE:: seqpolty(:)
  character(MAXSTRLEN) seqfile,fycfile,pdbinfile,xtcinfile,dcdinfile,netcdfinfile
!
end module sequen
!
