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
! pdbxyzs            : data structure for storing specific atoms by name per residue
! pdb_fileformat     : what type of file format to process in trajectory analysis 
! pdb_readmode       : how to process info in PDB (1: infer dihedrals; 2: read xyz; 3: read xyz with more rebuilding poss.)
! pdb_hmode          : whether to let CAMPARI construct hydrogens completely afresh
! pdb_nucmode        : assumed convention for which residue O3* belongs to in nucleic acids
! pdb_convention     : input/output naming conventions for PDB files (AMBER/CHARMM/etc)
! pdb_formstr        : input/output format strings (allows extension of coordinate precision -> nonstandard PDB)
! pdb_tolerance      : on read, tolerances for default (hard-coded) bond angles/lengths -> otherwise triggers rebuild
! pdb_mpimany        : on read, whether different MPI replicas use different inputs or not
! pdb_suff/pdb_pref  : suffix and prefix when dealing with single PDB files on read
! pdb_format         : helper array for single PDB file name construction
! dcd_withbox,xtc_prec,netcdf_ids : special flags/helpers for specific binary input formats
! use_pdb_template   : whether PDB map is required and has been read for binary input trajectories (also unsupp. res.)
! pdbtmplfile        : the file name
! pdbmap             : map for reordering on read
! pdboutmap          : map for reordering on write
! pdb_writemode      : whether to write several files or single trajectory (usually tied rigidly to format)
! pdb_force_box      : on both read/write, whether to assume/put everything in the unit cell or leave mol. intact (default)
! just_solutes       : on write, whether to skip solvent molecules
! use_trajidx        : on write, whether user-selected list of atoms is written 
! trajidxfile        : the file name
! pdbeffn            : effective number of atoms for structural output (FMCSC_TRAJIDXFILE)
! pdboutlst,pdboutlog: list and "present-in-output" vectors (FMCSC_TRAJIDXFILE)
! n_pdbunk           : number of unsupported residues in sequence file (PDB input required)
! pdb_unknonws       : their names
! pdb_unkbnd         : data structure containing their bounds in input lists
! pdb_rmol           : reference molecule for image selection for structural output (and other things)
! pdbxyzs            : data structure for storing specific atoms by name per residue
! xpdb,ypdb,zpdb     : PDB coordinate arrays
! select_frames      : whether to anaylze only user-selected frames in traj. analysis (FMCSC_FRAMESFILE)
! frameidxfile       : the file name
! framelst,framelst2 : list of user-selected frames
! framecnt           : size of list
! curframe,netcdf_fr1: helpers for list management
! framewts,framewts2 : FP weights per user-selected frames
!
module pdb
!
  type t_pdbxyzs
    RTYPE nl(3),cal(3),cl(3),ol(3),no4(3),oz(3),oe2(3),od2(3)
    RTYPE n(3),ca(3),c(3),cb(3),cg(3),cg1(3),hnz1(3),hd2(3),he2(3)
    RTYPE og(3),sg(3),hg(3),hg1(3),og1(3),ch3(3),ch1(3),ch(3),ci(3)
    RTYPE cd(3),cd1(3),nd1(3),od1(3),sd(3)
    RTYPE ce(3),oe1(3),ne(3),cz(3),nz(3),n1(3)
    RTYPE o(3),ow(3),hw1(3),hw2(3),no1(3),h(3)
    RTYPE ht2(3),ht1(3),ct(3),c1(3),c21(3),ho(3)
    RTYPE hnt1(3),hct1(3),cct(3),cnt(3),c31(3)
    RTYPE c4(3),h21(3),h11(3),c2(3),hbt1(3),cbt(3)
    RTYPE oxt(3),ce2(3),hh(3),oh(3),h1(3),nd(3)
    RTYPE hn1(3),hn2(3),oxt1(3),n2(3),c6(3),n9(3)
    RTYPE c8(3),no5(3),np(3),nc5(3),nc4(3),nc3(3)
    RTYPE no3(3),nc3l(3),nc4l(3),nhop(3),no32(3)
    RTYPE nho3(3),no3l(3),no32l(3),nc1(3),nc2(3)
    RTYPE no2(3),no1s(3),nho2(3),nho1(3),nc5l(3)
    RTYPE nho5(3),thyh1(3),c5m(3),c5(3),s(3),c3(3),c12(3)
    RTYPE ct1(3),ct2(3),cb1(3),cb2(3),hs(3),ht11(3),c11(3)
    RTYPE ht12(3),ht13(3),ct3(3),o1(3),o2(3),chl(3)
  end type t_pdbxyzs
!
  integer pdb_format(4),pdb_fileformat,pdb_readmode,pdb_writemode,n_pdbunk,pdb_ihlp,pdb_rmol
  integer pdb_hmode,pdb_nucmode,pdb_convention(2),pdb_force_box,pdbeffn,framecnt,curframe,netcdf_fr1
#ifdef LINK_NETCDF
  integer netcdf_ids(40)
#endif
  logical use_pdb_template,just_solutes,use_trajidx,select_frames,use_frame_weights
  RTYPE, ALLOCATABLE:: xpdb(:),ypdb(:),zpdb(:),framewts(:),framewts2(:)
  character(3), ALLOCATABLE:: pdb_unknowns(:)
  integer, ALLOCATABLE:: pdbmap(:),pdboutlst(:),framelst(:),framelst2(:),pdb_unkbnd(:,:),pdboutmap(:),pdb_didread(:,:),pdbmap2(:)
  logical, ALLOCATABLE:: pdboutlog(:)
  RTYPE pdb_tolerance(3)
  logical dcd_withbox,pdb_mpimany
#ifdef LINK_XDR
  RTYPE xtc_prec
#endif
  character(MAXSTRLEN) pdb_suff,pdb_pref,pdbtmplfile,trajidxfile,frameidxfile,pdb_formstr(2)
!
end module pdb
!
