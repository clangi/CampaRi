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
! CONTRIBUTIONS: Albert Mao, Adam Steffen                                  !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! mc_acc_crit     : what acceptance criterion is in use
! xxx_move        : what type of MC move is current
! xxx_randfreq    : fraction of certain movetype to be fully randomizing
! xxx_freq        : tree-based fraction of move types (see doc.)
! xxx_stepsz      : step size of certain movetype for local perturbation
! cur_trans       : current rigid-body translation
! cur_rot         : current rigid-body rotation move
! trans_xyzbuf    : buffer (scale) for fully randomizing translation for non-periodic boundaries
! cur_clrgpcs     : current gyration tensor EVs for molecular cluster
! cur_clcom       : current center of mass for molecular cluster
! clurb_clus      : list of molecule clusters according to last call to molclusters2(...)
! clurb_cluszs    : vector of molecule cluster sizes according to last call to molclusters2(...)
! clurb_clun      : number of molecule clusters according to last call to molclusters2(...)
! align_NC        : how to align the chain for pivot-type moves (includes CR, omega)
! use_xxx         : certain special flags for individual movetypes (see mcmove.f)
! have_xxx        : indicators used to globally indicate whether a move type has any candidate entities
! Concerted rotation variables: Sjunesson et al. (implementation by Xiaoling Wang)
! roa    : rotation matrix with alph (N-Cai-C bond angle)
! rob    : rotation matrix with beta (Cai-C-N bond angle)
! bnc    : nitrogen alph carbon bond length
! bcc    : alpha carbon carbonyl carbon bond length
! bcn    : carbonyl carbon nitrogen carbon bond length
! bco    : carbonyl carbon oxygen bond length
! cco    : Cai-C-O bond angle
! mcn    : constant rotation matrix about peptide bond
!
module movesets
!
  integer MAXCRDOF
  parameter (MAXCRDOF=20)
!
! general frequencies, flags, and files
  RTYPE chifreq,crfreq,pivot_stepsz
  RTYPE nuc_stepsz,nuc_randfreq,nucfreq,nuccrfreq,nucpuckfreq,nucpuckrdfreq
  RTYPE omega_stepsz,omega_randfreq,omegafreq,puckerfreq,puckerrdfreq,pucker_anstp
  RTYPE pivot_randfreq,chi_stepsz,chi_randfreq,rigid_randfreq,pucker_distp
  RTYPE trans_stepsz,rot_stepsz,rigidfreq,rotfreq,phfreq,angcrfreq,otherfreq,other_randfreq,other_stepsz
  RTYPE particleflucfreq,clurb_freq,clurb_rdfreq,torcrfreq,torcrfreq_omega,other_natfreq,other_unkfreq
  logical use_globmoves,chi_move,pivot_move,SJconrot_move,nuccr_move
  logical fyc_move,lct_move,nuc_move,torcr_move_omega
  logical use_lctmoves,rigid_move,UJconrot_move,pucker_move,torcr_move
  logical rot_move,trans_move,use_ph,ph_move,omega_move,nucpuck_move,other_move
  logical particleinsertion_move, particledeletion_move
  logical particleidentity_move,clurb_move,in_dyncyc
  logical have_omega,have_chi,have_nuc,have_pucker,have_nucpuck,have_torcr
  logical have_djcr,have_docr,have_nuccr,have_sjcr,have_ujcr,have_rigid,have_clurb
  logical have_ph,have_particlefluc,have_pivot,have_cr,psw_report,have_other,have_unkother,have_natother,have_unsother
  integer mc_acc_crit,align_NC
  character(MAXSTRLEN) prefsamplingfile
!
! current move structure
  type t_mvcur
    integer rs,rsi,rsf,imol,dofix,nbllen,solc,dof,dof2
    RTYPE diff,bias,jac_a,jac_b,det_L,det_L_star,dsquare,dsquare_star,hlp(10)
    RTYPE, ALLOCATABLE:: evecp(:),eveca(:),evecd(:),covec(:)
    logical ok,ct,flip,isnat,solf
  end type t_mvcur
  type(t_mvcur):: mvcur
!
! counter structure
  type t_mvcnt
    integer nchi,nfy,nujcr,nsjcr,ndjcr,ndocr,nre,nomega,nlct,nclurb,nrb,ntrans,nrot,nnuc,nnuccr
    integer ninsert,ndelete,nidentitychange,nld,nmd,nbd,nph,nfyc,npucker,nother
    integer ndynseg,nmcseg
    RTYPE avgdynseglen,avgmcseglen
  end type t_mvcnt
  type(t_mvcnt):: mvcnt
!
! list of eligible entities for implementing picking probabilities
  type t_entlst
    integer nr ! number of eligible entities
    integer, ALLOCATABLE:: idx(:),idx2(:) ! indices of eligible entities, idx2 is auxiliary for CR
    RTYPE, ALLOCATABLE:: wt(:) ! corresponding FP weights for pref. sampling
  end type t_entlst
!
! concerted rotation parameters
  integer nr_crdof,nr_crfit,nr_crres,cr_mode,sj_wrncnt(10),sj_wrnlmt(10)
  RTYPE cr_a,cr_b,roa(3,3),rob(3,3),mcn(3,3),cco,bco,bcc,bcn,bnc
  type(t_entlst):: ujcrlst,djcrlst,sjcrlst,docrlst,nuccrlst
!
! pivot move parameters
  logical use_stericgrids
  type(t_entlst):: fylst,unklst,unslst,natlst
  RTYPE cur_other
!
! chi-move parameters
  integer nrchis
  type(t_entlst):: chilst
!
! omega move parameters
  type(t_entlst):: wlst
!
! nuc-move parameters
  integer nrnucim
  type(t_entlst):: nuclst
!
! nuc-pucker move parameters
  type(t_entlst):: nucpuclst
!
! pucker move parameters
  type(t_entlst):: puclst
!
! rigid-body move parameters
  integer clurb_maxsz,clurb_clun,clurb_incr,clurb_reset
  integer, ALLOCATABLE:: clurb_clus(:,:),clurb_cluszs(:)
  integer, ALLOCATABLE:: clurbf(:),clurbi(:)
  RTYPE cur_trans(3),cur_rot(3),cur_clcom(3),cur_clrgpcs(3,3)
  RTYPE clurb_strbail,clurb_discrit,trans_xyzbuf
  RTYPE, ALLOCATABLE:: cur_clurbt(:,:)
  logical use_coupledrigid
  type(t_entlst):: rblst,clurblst
!
! GC move parameters
  type(t_entlst):: fluclst
!
! hybrid dynamics move parameters
  integer curcyc_end,curcyc_start,max_dyncyclen,min_dyncyclen,max_mccyclen,min_mccyclen
  integer first_mccyclen
!
! constraint requests
  character(MAXSTRLEN) frzfile
  logical do_frz,frz_report,skip_frz
!
end module movesets
!
