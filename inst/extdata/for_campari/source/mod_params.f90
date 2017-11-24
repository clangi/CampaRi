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
! n_ljtyp      : number of unique LJ types in parameter file
! n_ctyp       : number of unique charge types in parameter file
! n_biotyp     : number of unique biotypes in parameter file
! n_botyp      : number of unique bonded atom types in parameter file
! n_bondtyp    : number of unique bond definitions in parameter file
! n_angltyp    : number of unique angle definitions in parameter file
! n_torstyp    : number of unique dihedral definitions in parameter file
! n_cmstyp     : number of unique CMAP definitions in parameter file
! reduce       : size reduction factor (universal for all atoms)
! sigrule      : combination rule for LJ sigmas (1 = arithmetic, 2 = geometric)
! epsrule      : combination rule for LJ epsilons (1 = arithmetic, 2 = geometric)
! bio_ljtyp    : assigned LJ type per biotype
! bio_ctyp     : assigned charge type per biotype
! bio_botyp    : assigned bonded type per biotype
! lj_atnum     : assigned proton number (nucleus) per LJ type
! lj_val       : assigned valence per LJ type
! lj_describe  : assigned description per LJ type
! lj_symbol    : assigned atomic symbol per LJ type
! lj_weight    : assigned mass per LJ type
! bio_code     : assigned atomic symbol per biotype
! bio_res      : assigned residue string per biotype
! c_charge     : assigned partial charge (in e) per charge type
! c_note       : assigned description per charge type
! lj_rad       : assigned self-term radius per LJ type
! lj_eps       : assigned pairwise LJ interaction parameter per pair of LJ types
! lj_sig       : assigned pairwise LJ size parameter per pair of LJ types
! bo_typ       : type of bond length potential per bond definition
! ba_typ       : type of bond angle potential per angle definition
! di_typ       : type of dihedral angle potential per torsion definition
! cm_typ       : type and bin number of CMAP potential per CMAP definition
! bo_par       : parameters for bond length potential per bond definition
! ba_par       : parameters for bond angle potential per angle definition
! di_par       : parameters for dihedral angle potential per torsion definition
! cm_par       : grid parameters for CMAP potential per CMAP definition
! cm_par2      : bin width for CMAP potential per CMAP definition
! cm_file      : filename for location of grid data relative to FMCSC_CMAPDIR
! cm_dir       : directory name for location of all CMAP files
! bo_lst       : code of bond length potential per bonded type pair
! ba_lst       : code of bond angle potential per bonded type triple
! di_lst       : code of dihedral angle potential per bonded type quadruple
! impt_lst     : code of improper dihedral angle potential per bonded type quadruple
! cm_lst       : code of CMAP potential per bonded type quintuple
! improper_conv: integers used to switch between AMBER/OPLS and CHARMM/GROMOS improper conventions
! biotpatchfile: string to input file for applying patches to biotype assignments (has dep. changes)
! cpatchfile   : string to input file for applying charge patches (override of biotype-derived partial charges)
! bpatchfile   : string to input file for applying patches to bonded parameters
! ljpatchfile  : string to input file for applying patches to LJ parameters
! masspatchfile: string to input file for applying patches to atomic masses
! radpatchfile : string to input file for applying patches to atomic radii
! ncpatchfile  : string to input file for applying patches to residue-level net charge flags
! guess_bonded : whether or not to make up bonded parameters from Engh-Huber or input geometries
! be_unsafe    : whether or not to disable a few sanity checks and turn them into warnings
! *_patched    : whether override parameters where read from specific input files
! xxx_report   : report flags
!
module params
!
  integer MAXTYP,MAXBOPAR,MAXBAPAR,MAXDIPAR
  parameter (MAXTYP=1000)
  parameter (MAXBOPAR=3)
  parameter (MAXBAPAR=4)
  parameter (MAXDIPAR=9)
!
  character(MAXSTRLEN) paramfile
!
  integer n_ljtyp,n_ctyp,n_biotyp,n_bondtyp,n_angltyp,n_torstyp,n_cmstyp,improper_conv(2)
  integer n_botyp,ba_lstsz,bo_lstsz,di_lstsz,impt_lstsz,cm_lstsz,maxcmpar,cm_splor
  integer, ALLOCATABLE:: bio_ljtyp(:),bio_ctyp(:),bio_botyp(:)
  integer, ALLOCATABLE:: lj_atnum(:),lj_val(:)
  character(3), ALLOCATABLE:: lj_symbol(:),bio_code(:)
  character(20), ALLOCATABLE:: lj_describe(:),bio_res(:),c_note(:)
  RTYPE, ALLOCATABLE:: lj_weight(:),c_charge(:)
  RTYPE, ALLOCATABLE:: ba_par(:,:),bo_par(:,:),di_par(:,:)
  RTYPE, ALLOCATABLE:: cm_par(:,:,:,:),cm_par2(:)
  character(MAXSTRLEN), ALLOCATABLE:: cm_file(:)
  character(MAXSTRLEN) cm_dir,cpatchfile,bpatchfile,ljpatchfile,ncpatchfile,masspatchfile,radpatchfile,biotpatchfile
  integer, ALLOCATABLE:: ba_typ(:),bo_typ(:),di_typ(:),cm_typ(:,:)
  integer, ALLOCATABLE:: ba_lst(:,:),bo_lst(:,:),di_lst(:,:) 
  integer, ALLOCATABLE:: impt_lst(:,:),cm_lst(:,:)
!
  integer epsrule,sigrule,guess_bonded
  logical vdW_report,elec_report,dip_report,bonded_report,be_unsafe,charge_patched,bonded_patched
  logical lj_patched,nc_patched,rad_patched,mass_patched,bt_patched
  RTYPE, ALLOCATABLE:: lj_eps(:,:),lj_sig(:,:),lj_rad(:),lj_eps_14(:,:),lj_sig_14(:,:)
  RTYPE reduce
!
end module params
!
