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
module fos
!
!
! parameters (all in kcal/mol) and global variables for implicit solvent model
! we keep those in capitalized pseudo-globals rather than an array structure,
! because none of the organizational units in the program (residues, molecules, ...)
! matches the model compound parsing exactly
!
  type t_savreq
    integer nats,instfreq
    RTYPE, ALLOCATABLE:: hists(:,:,:)
    integer, ALLOCATABLE:: idx(:)
    character(MAXSTRLEN) filen
  end type t_savreq
!
  type(t_savreq) savreq
  logical fos_report,ard_patched,asm_patched,fos_patched,use_Tdepfosvals
  RTYPE fos_Tdepref
  character(MAXSTRLEN) fospatchfile,asmpatchfile,ardpatchfile
!
  RTYPE, ALLOCATABLE:: freesolv(:,:)
! polypeptide backbone
  RTYPE FOS_PEP_BB(3),FOS_PEP_PROBB(3),FOS_PEP_BB_FOR(3)
  RTYPE FOS_PEP_BB_NH2(3),FOS_PEP_PROBB_FOR(3),FOS_PEP_CCT(3)
  RTYPE FOS_PEP_UCT(3),FOS_PEP_CNT(3),FOS_PEP_PROCNT(3)
  RTYPE FOS_PEP_UNT(3),FOS_PEP_PROUNT(3)
! polypeptide sidechains
  RTYPE FOS_ALA(3),FOS_VAL(3),FOS_LEU(3),FOS_ILE(3),FOS_HIE(3),FOS_HID(3)
  RTYPE FOS_MET(3),FOS_SER(3),FOS_THR(3),FOS_CYS(3),FOS_NVA(3),FOS_NLE(3)
  RTYPE FOS_PHE(3),FOS_TYR(3),FOS_TRP(3),FOS_ASN(3),FOS_PRO(3),FOS_GAM(3)
  RTYPE FOS_GLN(3),FOS_ASP(3),FOS_GLU(3),FOS_ARG(3),FOS_ORN(3),FOS_DAB(3)
  RTYPE FOS_LYS(3),FOS_HYP(3),FOS_ABA(3),FOS_PCA(3),FOS_AIB(3),FOS_HIP(3)
  RTYPE FOS_GLH(3),FOS_ASH(3),FOS_TYO(3),FOS_LYD(3),FOS_CYX(3),FOS_PTR(3)
  RTYPE FOS_TPO(3),FOS_SEP(3),FOS_KAC(3),FOS_KM1(3),FOS_KM2(3),FOS_KM3(3)
! polypeptide sidechain crosslinks
  RTYPE FOS_LINK_CC(3)
! small model compounds
  RTYPE FOS_ACA(3),FOS_NMF(3),FOS_NMA(3),FOS_URE(3),FOS_DMA(3),FOS_FOA(3),FOS_PPA(3)
  RTYPE FOS_SPC(3),FOS_T3P(3),FOS_T4P(3),FOS_T5P(3),FOS_T4E(3),FOS_CH4(3),FOS_MOH(3),FOS_PCR(3)
  RTYPE FOS_NA(3),FOS_CL(3),FOS_K(3),FOS_BR(3),FOS_I(3),FOS_CS(3),FOS_NH4(3),FOS_O2(3)
  RTYPE FOS_AC(3),FOS_GDN(3),FOS_CYT(3),FOS_URA(3),FOS_THY(3),FOS_GUA(3),FOS_ADE(3),FOS_PUR(3)
  RTYPE FOS_PRP(3),FOS_IBU(3),FOS_NBU(3),FOS_MSH(3),FOS_EOH(3),FOS_EMT(3),FOS_BEN(3),FOS_NAP(3)
  RTYPE FOS_TOL(3),FOS_IMD(3),FOS_IME(3),FOS_MIN(3),FOS_1MN(3),FOS_2MN(3),FOS_LCP(3),FOS_NO3(3)
! nucleic acid backbone
  RTYPE FOS_NUC_PO4M(3),FOS_NUC_PO4D(3),FOS_NUC_RIBO3(3),FOS_NUC_RIBO2(3)
  RTYPE FOS_NUC_RIBO1(3),FOS_NUC_MTHF(3)
! nucleic acid sidechain
  RTYPE FOS_XPU(3),FOS_XPT(3),FOS_XPC(3),FOS_XPG(3),FOS_XPA(3)
!
end module fos
