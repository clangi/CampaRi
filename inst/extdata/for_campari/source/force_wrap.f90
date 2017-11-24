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
!------------------------------------------------------------------------
!
!
!             #################################
!             #                               #
!             # FIRST: WRAPPER ROUTINES       #
!             #        FOR CARTESIAN FORCES   #
!             #        AND GLOBAL ENERGIES    #
!             #                               #
!             #################################
!
!
!------------------------------------------------------------------------
!
function force1(evec)
!
  use sequen
  use energies
  use forces
  use molecule
  use atoms
  use fyoc
  use system
  use iounit
  use polypep
  use wl
  use cutoffs, ONLY: t_ca_f,atinfo
  use atoms
!
  implicit none
!
  integer i,j,imol,rs,azero,aone,atwo,afour
  RTYPE evec(MAXENERGYTERMS),force1,edum
  logical sayno,sayyes,pflag
!
  atinfo(1,1:n) = x(1:n)
  atinfo(2,1:n) = y(1:n)
  atinfo(3,1:n) = z(1:n)
!
  azero = 0
  aone = 1
  atwo = 2
  afour = 4
  evec(:) = 0.0
  bnd_f(:,:) = 0.0
  bnd_fr = 0.0
  sayno = .false.
  sayyes = .true.
! for later use: pass on to inner_loops fxns (virial is non-trivial for PBC!)
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    pflag = .true.
  else
    pflag = .false.
  end if
  ens%insVirT(:,:) = 0.0
  cart_f(:,:) = 0.0
  t_ca_f(:,:) = 0.0
  if (do_accelsim.EQV..true.) then
    hmjam%ca_f(:,:) = 0.0
    hmjam%t_ca_f(:,:) = 0.0
  end if
  if (use_IMPSOLV.EQV..true.) then
    sisa(:)%nix = 0
    sum_scrcbs(:) = 0.0
    sum_scrcbs_tr(:) = 0.0
    sum_scrcbs_lr(:) = 0.0
    sav_dr(:,:) = 0.0
  end if
!
  if (use_IMPSOLV.EQV..true.) call init_svte(3)
!
  do i=1,nseq
    do j=i,nseq
      if (disulf(i).eq.j) cycle
      call Vforce_rsp(evec,i,j,sayno,cart_f,sisa)
    end do
  end do
! crosslinks will always be treated as short-range 
  do i=1,n_crosslinks
    call Vforce_rsp2(evec,crosslink(i)%rsnrs(1),crosslink(i)%rsnrs(2),sayno,t_ca_f,svte,sisa)
  end do
!
  if (use_IMPSOLV.EQV..true.) then
    call init_svte(6)
    call Vforce_dsavdr(sav_dr,azero)
    if (use_POLAR.EQV..true.) then
      call force_setup_scrqs2(azero)
    end if
  end if
!
  do i=1,nseq
    do j=i,nseq
      if (disulf(i).eq.j) cycle
      call Vforce_rsp_long(evec,i,j,cart_f,sum_scrcbs)
    end do
  end do
! crosslinks will always be treated as short-range 
  do i=1,n_crosslinks
    call Vforce_rsp_long2(evec,crosslink(i)%rsnrs(1),crosslink(i)%rsnrs(2),t_ca_f,sum_scrcbs)
  end do
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    do i=1,nseq
      call Vforce_dcbdscrq(i,t_ca_f)
    end do
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(13).EQV..true.)) then
    i = 13
    hmjam%evec_thr(1,1) = sum(evec)
    call els_manage_one_t(i,hmjam%evec_thr(1,1),hmjam%t_ca_f,afour) ! everything so far is the total nonbonded force
    evec(13) = hmjam%evec_thr(1,1)
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(4).EQV..true.)) then
      do i=1,nseq
        call force_freesolv(i,evec,hmjam%t_ca_f)
      end do
      i = 4
      call els_manage_one_t(i,evec(4),t_ca_f,aone)
    else
      do i=1,nseq
        call force_freesolv(i,evec,t_ca_f)
      end do
    end if
  end if
!
  if (use_BOND(1).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(15).EQV..true.)) then
      do rs=1,nseq
        call force_bonds(rs,evec,hmjam%t_ca_f)
      end do
      i = 15
      call els_manage_one_t(i,evec(15),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_bonds(rs,evec,t_ca_f)
      end do
    end if
  end if
  if (use_BOND(2).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(16).EQV..true.)) then
      do rs=1,nseq
        call force_angles(rs,evec,hmjam%t_ca_f)
      end do
      i = 16
      call els_manage_one_t(i,evec(16),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_angles(rs,evec,t_ca_f)
      end do
    end if
  end if
  if (use_BOND(3).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(17).EQV..true.)) then
      do rs=1,nseq
        call force_impropers(rs,evec,hmjam%t_ca_f)
      end do
      i = 17
      call els_manage_one_t(i,evec(17),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_impropers(rs,evec,t_ca_f)
      end do
    end if
  end if
  if (use_BOND(4).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(18).EQV..true.)) then
      do rs=1,nseq
        call force_torsions(rs,evec,hmjam%t_ca_f)
      end do
      i = 18
      call els_manage_one_t(i,evec(18),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_torsions(rs,evec,t_ca_f)
      end do
    end if
  end if
  if (use_BOND(5).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(20).EQV..true.)) then
      do rs=1,nseq
        call force_cmap(rs,evec,hmjam%t_ca_f)
      end do
      i = 20
      call els_manage_one_t(i,evec(20),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_cmap(rs,evec,t_ca_f)
      end do
    end if
  end if
  if (use_CORR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(2).EQV..true.)) then
      do rs=1,nseq
        call force_corrector(rs,evec,hmjam%ca_f)
      end do
      i = 2
      call els_manage_one(i,evec(2),cart_f,aone)
    else
      do rs=1,nseq
        call force_corrector(rs,evec,cart_f)
      end do
    end if
  end if
!
  if (use_TOR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(7).EQV..true.)) then
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec,i,hmjam%t_ca_f)
        end if
      end do
      i = 7
      call els_manage_one_t(i,evec(7),t_ca_f,aone)
    else
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec,i,t_ca_f)
        end if
      end do
    end if
  end if
!
  if ((do_accelsim.EQV..true.).AND.((hmjam%isin(12).EQV..true.).OR.(hmjam%isin(2).EQV..true.))) then
    do rs=1,nseq
      call force_boundary_rs(rs,evec,aone,hmjam%t_ca_f)
    end do
    if (hmjam%isin(12).EQV..true.) then
      i = 12
      call els_manage_one_t(i,evec(12),t_ca_f,aone)
    end if
    if (hmjam%isin(2).EQV..true.) then
      i = 2
      call els_manage_one_t(i,evec(11),t_ca_f,aone)
    end if
  else
    do rs=1,nseq
      call force_boundary_rs(rs,evec,aone,t_ca_f)
    end do
  end if
!
  if (use_ZSEC.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(8).EQV..true.)) then
      do imol=1,nmol
        call force_zsec_gl(imol,evec,hmjam%t_ca_f)
      end do
      i = 8
      call els_manage_one_t(i,evec(8),t_ca_f,aone)
    else
      do imol=1,nmol
        call force_zsec_gl(imol,evec,t_ca_f)
      end do
    end if
  end if
!
  if (use_DSSP.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(19).EQV..true.)) then
      call force_dssp_gl(evec,hmjam%t_ca_f)
      i = 19
      call els_manage_one_t(i,evec(19),t_ca_f,aone)
    else
      call force_dssp_gl(evec,t_ca_f)
    end if
  end if
!
  if (use_EMICRO.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(21).EQV..true.)) then
      call force_emicro_gl(evec,hmjam%t_ca_f,azero)
      i = 21
      call els_manage_one_t(i,evec(21),t_ca_f,aone)
    else
      call force_emicro_gl(evec,t_ca_f,azero)
    end if
  end if
!
  if (use_POLY.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
      do imol=1,nmol
        call force_poly_gl(imol,evec,hmjam%t_ca_f)
      end do
      i = 14
      call els_manage_one_t(i,evec(14),t_ca_f,aone)
    else
      do imol=1,nmol
        call force_poly_gl(imol,evec,t_ca_f)
      end do
    end if
  end if
!
  if (use_DREST.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(10).EQV..true.)) then
      call force_drest(evec,hmjam%t_ca_f,azero)
      i = 10
      call els_manage_one_t(i,evec(10),t_ca_f,aone)
    else
      call force_drest(evec,t_ca_f,azero)
    end if
  end if
!
  force1 = sum(evec(:))
  cart_f(:,:) = cart_f(:,:) + transpose(t_ca_f(:,:))
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(11).EQV..true.)) then
    hmjam%ca_f(:,:) = cart_f(:,:)
    edum = force1
    i = 11
    call els_manage_one(i,edum,cart_f,atwo)
    evec(11) = hmjam%boosts(11,1)
    force1 = force1 + evec(11)
  end if
!
end
!
!-----------------------------------------------------------------------
!
function force3(evec,evec_tr,evec_lr,forceflag)
!
  use sequen
  use energies
  use forces
  use molecule
  use atoms
  use iounit
  use cutoffs
  use mcsums
  use inter
  use polypep
  use units
  use fyoc
  use system
  use movesets
#ifdef ENABLE_MPI
  use mpistuff
#endif
  use wl
  use tabpot
!
  implicit none
!
  integer i,j,imol,rs,tpi,azero,aone,atwo,athree,afour
  RTYPE evec(MAXENERGYTERMS),force3
  RTYPE evec_tr(MAXENERGYTERMS),evec_lr(MAXENERGYTERMS)
  logical sayno,sayyes,forceflag,freshnbl,pflag
  RTYPE evec_thr(MAXENERGYTERMS,3),edum
  integer(KIND=8) t1,t2
!
  atinfo(1,1:n) = x(1:n)
  atinfo(2,1:n) = y(1:n)
  atinfo(3,1:n) = z(1:n)
  if (use_POLAR.EQV..true.) atinfo(4,1:n) = atq(1:n)
  if (use_IMPSOLV.EQV..true.) then
    if (par_IMPSOLV2(1).eq.1) then
      atsavinfo(1,1:n) = atr(1:n)
      atsavinfo(2,1:n) = atsavred(1:n)*atvol(1:n)
    else if (par_IMPSOLV2(1).eq.2) then
      atsavinfo(1,1:n) = atr(1:n)
      atsavinfo(2,1:n) = atsavred(1:n)
      atsavinfo(3,1:n) = atvol(1:n)
    end if
  end if
!
! first determine whether we need to update the NB-list
  freshnbl = .false.
  if ((mod(nstep,nbl_up).eq.0).OR.(nstep.le.1).OR.&
 &    (forceflag.EQV..true.)) then
    freshnbl = .true.
  end if
!
! init
  azero = 0
  aone = 1
  atwo = 2
  athree = 3
  afour = 4
  evec_thr(:,:) = 0.0
  tpi = 0
  evec(:) = 0.0
  bnd_fr = 0.0
  sayno = .false.
  sayyes = .true.
! for later use: pass on to inner_loops fxns (virial is non-trivial for PBC!)
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    pflag = .true.
  else
    pflag = .false.
  end if
  ens%insVirT(:,:) = 0.0
  bnd_f(:,:) = 0.0
  cart_f(:,:) = 0.0
  t_ca_f(:,:) = 0.0
! the mid- and long-range interaction arrays (energies and cartesian forces) are only updated
! every so many steps
  if ((freshnbl.EQV..true.)) then
    if ((use_mcgrid.EQV..true.).AND.(use_cutoffs.EQV..true.)) rs_nbl(:)%ntmpanb = 0
    evec_tr(:) = 0.0
    cart_f_tr(:,:) = 0.0
    t_ca_f_tr(:,:) = 0.0
    if (use_IMPSOLV.EQV..true.) sum_scrcbs_tr(:) = 0.0
  else
    if (lrel_md.eq.2) then
      cart_f_tr(:,:) = 0.0
      evec_lr(:) = 0.0
    end if
  end if
  if ((freshnbl.EQV..true.)) then
    evec_lr(:) = 0.0
    if (use_IMPSOLV.EQV..true.) sum_scrcbs_lr(:) = 0.0
  end if
  if (do_accelsim.EQV..true.) then
    hmjam%ca_f(:,:) = 0.0
    hmjam%t_ca_f(:,:) = 0.0
  end if
  if (use_IMPSOLV.EQV..true.) then
    sisa(:)%nix = 0
    sum_scrcbs(:) = 0.0
    sav_dr(:,:) = 0.0
  end if
!
  if ((use_cutoffs.EQV..true.).AND.(freshnbl.EQV..true.).AND.(ideal_run.EQV..false.)) then
    call System_Clock(t1)
    call all_respairs_nbl(skip_frz,is_tab,tpi)
!    if (forceflag.EQV..true.) call tryit()
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t1 - t2)
    time_nbl = time_nbl + 1.0*(t2 - t1)
  end if
!
  if (is_plj.EQV..true.) then
!
!   first calculations which use the standard NB gas phase potential (Coulomb+LJ)
!   note that many additional potentials can still be handled including EXTRA (CORR)
!   these routines are strictly incompatible with IMPSOLV, TABUL, and WCA, but will
!   also seg-fault when called WITHOUT the Coulomb term
!
    if (use_cutoffs.EQV..true.) then
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,afour)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,athree)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              if (is_3site.EQV..true.) then
                call Vforce_PLJ_C_T3P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_4site.EQV..true.) then
                call Vforce_PLJ_C_T4P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_5site.EQV..true.) then
                call Vforce_PLJ_C_T5P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              end if
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            if (is_3site.EQV..true.) then
              call Vforce_PLJ_C_T3P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_4site.EQV..true.) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_5site.EQV..true.) then
              call Vforce_PLJ_C_T5P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            end if
          end if
        end do
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
               call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
               call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
        end do
      end if
!     add long-range electrostatics corrections
      call force_P_LR(evec_thr(:,3),t_ca_f_tr,sum_scrcbs_lr,freshnbl,azero)
!
    else
!     no cutoffs and no simplifying assumptions
      do i=1,nseq
        call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
        call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
        if (i.lt.nseq) then
          call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
        end if
        if (i.lt.(nseq-1)) then
          call Vforce_PLJ(evec_thr(:,1),i,t_ca_f)
        end if
      end do
    end if
!   crosslinks will always be treated as short-range 
    if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1),sayno,t_ca_f,svte,sisa,sum_scrcbs)
!
  else if (is_pewlj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle Ewald sums
!
    if (use_cutoffs.EQV..true.) then
!     note that we enforce during parsing that there is no twin-regime, so we
!     don't have to check it here.
!     also note that this is potentially a weak assumption for LJ interactions,
!     but that Ewald is generally safer if more density is in the real-space 
!     component using a slightly larger SR cutoff anyway
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call force_rsp_long_mod(evec_thr(:,1),i,i+1,sayno,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,athree)
          end if
        end do
        do i=rsw1,nseq
          call force_rsp_long_mod(evec_thr(:,1),i,i,sayno,cart_f)
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            if (is_3site.EQV..true.) then
              call Vforce_PLJ_C_T3P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_4site.EQV..true.) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_5site.EQV..true.) then
              call Vforce_PLJ_C_T5P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            end if
          end if
        end do
      else
        do i=1,nseq
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,sayno,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call force_rsp_long_mod(evec_thr(:,1),i,i+1,sayno,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
        end do
      end if
!     add reciprocal sum/forces and constant term to energy
      t_ca_f_tr(:,:) = 0.0
      call force_P_LR(evec_thr(:,3),t_ca_f_tr,sum_scrcbs_lr,freshnbl,azero)
!     crosslinks will always be treated as short-range 
      if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1),sayno,t_ca_f,svte,sisa,sum_scrcbs)
!
    else
!     Ewald must have a real-space cutoff
      write(ilog,*) 'Fatal. Ewald summation requires the use of a cutoff.'
      call fexit()
    end if
!
  else if (is_prflj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle (G)RF electrostatics
!
    if (use_cutoffs.EQV..true.) then
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call force_rsp_long_mod(evec_thr(:,2),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,afour)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call force_rsp_long_mod(evec_thr(:,1),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,athree)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              if (is_3site.EQV..true.) then
                call Vforce_PLJ_C_T3P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_4site.EQV..true.) then
                call Vforce_PLJ_C_T4P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_5site.EQV..true.) then
                call Vforce_PLJ_C_T5P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              end if
            end if
          end if
          call force_rsp_long_mod(evec_thr(:,1),i,i,use_cutoffs,cart_f)
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            if (is_3site.EQV..true.) then
              call Vforce_PLJ_C_T3P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_4site.EQV..true.) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_5site.EQV..true.) then
              call Vforce_PLJ_C_T5P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            end if
          end if
        end do
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call force_rsp_long_mod(evec_thr(:,2),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call force_rsp_long_mod(evec_thr(:,1),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
        end do
      end if
!     just increments energy
      call force_P_LR(evec_thr(:,3),t_ca_f_tr,sum_scrcbs_lr,freshnbl,azero)
!     crosslinks will always be treated as short-range 
      if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1),use_cutoffs,t_ca_f,svte,sisa,sum_scrcbs)
!
    else
!     RF must have a real-space cutoff
      write(ilog,*) 'Fatal. (G)RF corrections require the use of a cutoff.'
      call fexit()
    end if
!
  else if (is_fegplj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle ghosted Hamiltonians
!   essentially, separate neighbor lists are set up for ghosted interactions
!
    if (use_cutoffs.EQV..true.) then
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,afour)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%ngnbtrs,rs_nbl(i)%ngnbtrats,atwo)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,athree)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%ngnbs,rs_nbl(i)%ngnbats,aone)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              if (is_3site.EQV..true.) then
                call Vforce_PLJ_C_T3P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_4site.EQV..true.) then
                call Vforce_PLJ_C_T4P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_5site.EQV..true.) then
                call Vforce_PLJ_C_T5P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              end if
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%ngnbtrs,rs_nbl(i)%ngnbtrats,atwo)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            if (is_3site.EQV..true.) then
              call Vforce_PLJ_C_T3P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_4site.EQV..true.) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_5site.EQV..true.) then
              call Vforce_PLJ_C_T5P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            end if
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%ngnbs,rs_nbl(i)%ngnbats,aone)
          end if
        end do
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%ngnbtrs,rs_nbl(i)%ngnbtrats,atwo)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%ngnbs,rs_nbl(i)%ngnbats,aone)
          end if
        end do
      end if
!     add long-range electrostatics corrections
      call force_P_LR(evec_thr(:,3),t_ca_f_tr,sum_scrcbs_lr,freshnbl,azero)
!
    else
!     no cutoffs, no simplifying assumptions, and no speed-up
      do i=1,nseq
        do j=i,nseq
          if (disulf(i).eq.j) cycle
          call Vforce_rsp(evec_thr(:,1),i,j,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,j,cart_f,sum_scrcbs)
        end do
      end do
    end if
!   crosslinks will always be treated as short-range 
    if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1),sayno,t_ca_f,svte,sisa,sum_scrcbs)
!
  else if (is_fegprflj.EQV..true.) then
!
!   this is the modified tree with all routines being able to handle ghosted Hamiltonians + RF
!   for simple Cb-scaling (linear) only
!
!   cutoffs are mandatory due to RF
!
!   uses Vforce_PLJ_C, Vforce_FEGPLJ_C, force_rsp_long_mod, and force_P_LR
!
    if (use_cutoffs.EQV..true.) then
!
      if (use_waterloops.EQV..true.) then
        do i=1,rsw1-1
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call force_rsp_long_mod(evec_thr(:,2),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,afour)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%ngnbtrs,rs_nbl(i)%ngnbtrats,atwo)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call force_rsp_long_mod(evec_thr(:,1),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,athree)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%ngnbs,rs_nbl(i)%ngnbats,aone)
          end if
        end do
        do i=rsw1,nseq
          if ((freshnbl.EQV..true.)) then 
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%nwnbtrs.gt.0) then
              if (is_3site.EQV..true.) then
                call Vforce_PLJ_C_T3P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_4site.EQV..true.) then
                call Vforce_PLJ_C_T4P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              else if (is_5site.EQV..true.) then
                call Vforce_PLJ_C_T5P(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
              end if
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%ngnbtrs,rs_nbl(i)%ngnbtrats,atwo)
            end if
          end if
          call force_rsp_long_mod(evec_thr(:,1),i,i,use_cutoffs,cart_f)
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%nwnbs.gt.0) then
            if (is_3site.EQV..true.) then
              call Vforce_PLJ_C_T3P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_4site.EQV..true.) then
              call Vforce_PLJ_C_T4P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            else if (is_5site.EQV..true.) then
              call Vforce_PLJ_C_T5P(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
            end if
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%ngnbs,rs_nbl(i)%ngnbats,aone)
          end if
        end do
      else
        do i=1,nseq
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
                call force_rsp_long_mod(evec_thr(:,2),i,i+1,use_cutoffs,cart_f_tr)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
            if (rs_nbl(i)%ngnbtrs.gt.0) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%ngnbtrs,rs_nbl(i)%ngnbtrats,atwo)
            end if
          end if
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,use_cutoffs,cart_f)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
              call force_rsp_long_mod(evec_thr(:,1),i,i+1,use_cutoffs,cart_f)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_PLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
          if (rs_nbl(i)%ngnbs.gt.0) then
            call Vforce_FEGPLJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%ngnbs,rs_nbl(i)%ngnbats,aone)
          end if
        end do
      end if
!     add long-range electrostatics corrections
      call force_P_LR(evec_thr(:,3),t_ca_f_tr,sum_scrcbs_lr,freshnbl,azero)
!     crosslinks will always be treated as short-range 
      if (n_crosslinks.gt.0) call all_crosslink_corrs(evec_thr(:,1),sayno,t_ca_f,svte,sisa,sum_scrcbs)
    else
!     RF must have a real-space cutoff
      write(ilog,*) 'Fatal. (G)RF corrections require the use of a cutoff.'
      call fexit()
    end if
!
  else if ((is_lj.EQV..true.).OR.(is_ev.EQV..true.)) then
!
!   a tree for just LJ models: uses Vforce_LJ_C and Vforce_LJ
!
    if (use_cutoffs.EQV..true.) then
!
      do i=1,nseq
        if ((freshnbl.EQV..true.)) then
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.2) then
              call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
            end if
          end if
          if (rs_nbl(i)%nnbtrs.gt.0) then
            call Vforce_LJ_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
          end if
        end if
        call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
        if (i.lt.nseq) then
          if (rsp_vec(i).eq.1) then
            call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
          end if
        end if
        if (rs_nbl(i)%nnbs.gt.0) then
          call Vforce_LJ_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
        end if
      end do
!
    else
!     no cutoffs and no simplifying assumptions
      do i=1,nseq
        call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f) !cart_f,sisa)
        if (i.lt.nseq) then
          call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f) !cart_f,sisa)
        end if
        if (i.lt.(nseq-1)) then
          call Vforce_LJ(evec_thr(:,1),i,t_ca_f)
        end if
      end do
    end if
!   crosslinks will always be treated as short-range 
    do i=1,n_crosslinks
      call Vforce_rsp2(evec_thr(:,1),crosslink(i)%rsnrs(1),crosslink(i)%rsnrs(2),sayno,t_ca_f,svte,sisa)
    end do
!
!
  else if (is_tab.EQV..true.) then
! 
!   a modified tree designed to handle calculations with just tabulated potentials
!   including very sparse assignments
!
    if (use_cutoffs.EQV..true.) then
!
      do i=1,nseq
        if ((freshnbl.EQV..true.)) then
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.2) then
!              call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
              call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
            end if
          end if
          if (rs_nbl(i)%nnbtrs.gt.0) then
            imol = 0
            do j=1,rs_nbl(i)%nnbtrs
              imol = imol + tbp%rsmat(rs_nbl(i)%nbtr(j),i) - tbp%rsmat(i,rs_nbl(i)%nbtr(j)) + 1
            end do
            call Vforce_TAB_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,imol,atwo)
          end if
        end if
!        call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
        call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
        if (i.lt.nseq) then
          if (rsp_vec(i).eq.1) then
!            call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
            call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
          end if
        end if
        if (rs_nbl(i)%nnbs.gt.0) then
          imol = 0
          do j=1,rs_nbl(i)%nnbs
            imol = imol + tbp%rsmat(rs_nbl(i)%nb(j),i) - tbp%rsmat(i,rs_nbl(i)%nb(j)) + 1
          end do
          if (imol.gt.0) call Vforce_TAB_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,imol,aone)
        end if
      end do
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_rsp_long2(evec_thr(:,1),crosslink(i)%rsnrs(1),crosslink(i)%rsnrs(2),t_ca_f,sum_scrcbs)
      end do
!
    else ! no cutoffs
!
      do i=1,nseq
!        call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
        call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
        if (i.lt.nseq) then
!          call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
          call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
        end if
        if (rs_nbl(i)%ntabnbs.gt.0) then
          call Vforce_TAB(evec_thr(:,1),i,t_ca_f)
        end if
      end do
!     note that Vforce_TAB does not skip crrosslinked interactions and that no
!     corrections are needed
    end if
!
!
!   missing
!
  else if (ideal_run.EQV..true.) then
!
!   do nothing
!
  else
!
!   the standard calculation with support for arbitrary Hamiltonians
!   we'll do some subdistinctions for specialized routines, but keep the framework general 
!   any one-shot Hamiltonians should probably be separated out into a construct as above,
!   since distances will be calculated twice ...
!
    if (use_IMPSOLV.EQV..true.) call init_svte(3)
!
    if (use_cutoffs.EQV..true.) then
!
!     first the short-range terms (ATTLJ, IPP, WCA, FOS)
      if ((is_impljp.EQV..true.).OR.(is_implj.EQV..true.)) then
!       this is the specialized set of fxns for LJ 6/12 with ABSINTH implicit solvent model
        do i=1,nseq
!         if a nb-list update is due we recompute all short-range terms within second cutoff
!         note that this should never yield any extra terms in sisa%dr
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
              end if
            end if
            if (rs_nbl(i)%nnbtrs.gt.0) then
              call Vforce_LJIMP_C(evec_thr(:,2),i,t_ca_f_tr,svte,sisa,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
            end if
          end if
!         now compute true short-range terms
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
            end if
          end if
          if (rs_nbl(i)%nnbs.gt.0) then
            call Vforce_LJIMP_C(evec_thr(:,1),i,t_ca_f,svte,sisa,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
          end if
        end do
      else
!       the general case
        do i=1,nseq
!         if a nb-list update is due we recompute all short-range terms within second cutoff
!         note that this should never yield any extra terms in sisa%dr
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
              end if
            end if
            do j=1,rs_nbl(i)%nnbtrs
              call Vforce_rsp(evec_thr(:,2),i,rs_nbl(i)%nbtr(j),sayno,cart_f_tr,sisa)
            end do
            if (use_FEG.EQV..true.) then
              do j=1,rs_nbl(i)%ngnbtrs
                call Vforce_rsp(evec_thr(:,2),i,rs_nbl(i)%gnbtr(j),sayno,cart_f_tr,sisa)
              end do
            end if
          end if
!         now compute true short-range terms
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
            end if
          end if
          do j=1,rs_nbl(i)%nnbs
            call Vforce_rsp(evec_thr(:,1),i,rs_nbl(i)%nb(j),sayno,cart_f,sisa)
          end do
          if (use_FEG.EQV..true.) then
            do j=1,rs_nbl(i)%ngnbs
              call Vforce_rsp(evec_thr(:,1),i,rs_nbl(i)%gnb(j),sayno,cart_f,sisa)
            end do
          end if
        end do
      end if
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_rsp2(evec_thr(:,1),crosslink(i)%rsnrs(1),crosslink(i)%rsnrs(2),sayno,t_ca_f,svte,sisa)
      end do
!
!     setup the screened charge terms (which have the same short-range
!     as the sisa%dr terms unless large charge groups are used with coupled models:
!     in that latter case larger range is always accounted for fully(!))
      if (use_IMPSOLV.EQV..true.) then
        call init_svte(6)
        call Vforce_dsavdr(sav_dr,azero)
        if (use_POLAR.EQV..true.) then
          call force_setup_scrqs2(azero)
        end if
      end if
!
      if (is_implj.EQV..true.) then
!       do nothing here
      else if (is_impljp.EQV..true.) then
!       this is the specialized set of fxns for POLAR with ABSINTH implicit solvent model
        if ((scrq_model.le.2).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
          do i=1,nseq
!           compute the long-range terms within short cutoff
            call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
              end if
            end if
!
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM1256_C(evec_thr(:,1),i,t_ca_f,sum_scrcbs,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
            end if
!           and now add the long-range terms within twin-range
!           note that this has the funny effect that the terms in sum_scrcbs are
!           or are not incremented further on account of charge-charge interactions
!           between particles relatively far apart, i.e., it's a pseudo-cutoff
!           since the atom will still be close to one of the two charges 
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
!
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM1256_C(evec_thr(:,2),i,t_ca_f_tr,sum_scrcbs_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
              end if
            end if
          end do
        else if (scrq_model.eq.4) then
          do i=1,nseq
            call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM4_C(evec_thr(:,1),i,t_ca_f,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
            end if
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM4_C(evec_thr(:,2),i,t_ca_f_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
              end if
            end if
          end do
        else if ((scrq_model.eq.3).OR.((scrq_model.ge.7).AND.(scrq_model.le.9))) then
          do i=1,nseq
            call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
              end if
            end if
!
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PSCRM3789_C(evec_thr(:,1),i,t_ca_f,sum_scrcbs,aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
            end if
            if ((freshnbl.EQV..true.)) then
              if (i.lt.nseq) then
                if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
                end if
              end if
!
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PSCRM3789_C(evec_thr(:,2),i,t_ca_f_tr,sum_scrcbs_tr,aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
              end if
            end if
          end do
        else
          write(ilog,*) 'Fatal. Encountered unsupported screening model in force3(...).&
                       & This is most certainly an omission bug.'
          call fexit()
        end if
      else
!       the general case
        do i=1,nseq
!         now compute the long-range terms within short cutoff
          call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
          if (i.lt.nseq) then
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
            end if
          end if
!
          do j=1,rs_nbl(i)%nnbs
            call Vforce_rsp_long(evec_thr(:,1),i,rs_nbl(i)%nb(j),cart_f,sum_scrcbs)
          end do
          if (use_FEG.EQV..true.) then
            do j=1,rs_nbl(i)%ngnbs
              call Vforce_rsp_long(evec_thr(:,1),i,rs_nbl(i)%gnb(j),cart_f,sum_scrcbs)
            end do
          end if
!         and now add the long-range terms within twin-range
!         note that this has the funny effect that the terms in sum_scrcbs are
!         or are not incremented further on account of charge-charge interactions
!         between particles relatively far apart, i.e., it's a pseudo-cutoff
!         since the atom will still be close to one of the two charges 
          if ((freshnbl.EQV..true.)) then
            if (i.lt.nseq) then
              if (rsp_vec(i).eq.2) then
                call Vforce_rsp_long(evec_thr(:,2),i,i+1,cart_f_tr,sum_scrcbs_tr)
              end if
            end if
!
            do j=1,rs_nbl(i)%nnbtrs
              call Vforce_rsp_long(evec_thr(:,2),i,rs_nbl(i)%nbtr(j),cart_f_tr,sum_scrcbs_tr)
            end do
            if (use_FEG.EQV..true.) then
              do j=1,rs_nbl(i)%ngnbtrs
                call Vforce_rsp_long(evec_thr(:,2),i,rs_nbl(i)%gnbtr(j),cart_f_tr,sum_scrcbs_tr)
              end do
            end if
          end if
        end do
      end if
!     finally, add long-range electrostatics corrections (note that these include
!     the next neighbor problem explicitly)
!     remember that IMPSOLV and FEG are mutually exclusive, if that ever changes
!     the construct inside the fxn below has to be changed as well!
      if (use_POLAR.EQV..true.) then
        call force_P_LR(evec_thr(:,3),t_ca_f_tr,sum_scrcbs_lr,freshnbl,azero)
      end if
!     crosslinks will always be treated as short-range 
      do i=1,n_crosslinks
        call Vforce_rsp_long2(evec_thr(:,1),crosslink(i)%rsnrs(1),crosslink(i)%rsnrs(2),t_ca_f,sum_scrcbs)
      end do
!
    else
!
      do i=1,nseq
        do j=i+2,nseq
          call Vforce_rsp2(evec_thr(:,1),i,j,sayno,t_ca_f,svte,sisa)
        end do
        call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
        if ((i.lt.nseq).AND.(disulf(i).ne.(i+1))) then
          call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa)
        else if (disulf(i).eq.(i+1)) then
          call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,t_ca_f,svte,sisa)
        end if
      end do
!
      if (use_IMPSOLV.EQV..true.) then
        call init_svte(6)
        call Vforce_dsavdr(sav_dr,azero)
        if (use_POLAR.EQV..true.) then
          call force_setup_scrqs2(azero)
        end if
      end if
!
      do i=1,nseq
        do j=i+2,nseq
          call Vforce_rsp_long2(evec_thr(:,1),i,j,t_ca_f,sum_scrcbs)
        end do
        call Vforce_rsp_long(evec_thr(:,1),i,i,cart_f,sum_scrcbs)
        if ((i.lt.nseq).AND.(disulf(i).ne.(i+1))) then
          call Vforce_rsp_long(evec_thr(:,1),i,i+1,cart_f,sum_scrcbs)
        else if (disulf(i).eq.(i+1)) then
          call Vforce_rsp_long2(evec_thr(:,1),i,i+1,t_ca_f,sum_scrcbs)
        end if
      end do
!
    end if
!
  end if
!
! the rest are cutoff-independent terms
!
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.)) then
    do i=1,nseq
      call Vforce_dcbdscrq(i,t_ca_f)
    end do
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(13).EQV..true.)) then
    i = 13
    if (freshnbl.EQV..true.) then
      hmjam%evec_thr(1,1) = sum(evec_thr)
    else
      hmjam%evec_thr(1,1) = sum(evec_thr) + evec_tr(1) + evec_tr(3) + evec_tr(5) + evec_tr(6) + evec_tr(9) + evec_lr(6) + evec_lr(9)
    end if
    call els_manage_one_t(i,hmjam%evec_thr(1,1),hmjam%t_ca_f,afour) ! everything so far is the total nonbonded force
    evec_thr(13,1) = hmjam%evec_thr(1,1)
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(4).EQV..true.)) then
      do i=1,nseq
        call force_freesolv(i,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 4
      call els_manage_one_t(i,evec_thr(4,1),t_ca_f,aone)
    else
      do i=1,nseq
        call force_freesolv(i,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
!
  if (use_BOND(1).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(15).EQV..true.)) then
      do rs=1,nseq
        call force_bonds(rs,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 15
      call els_manage_one_t(i,evec_thr(15,1),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_bonds(rs,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
  if (use_BOND(2).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(16).EQV..true.)) then
      do rs=1,nseq
        call force_angles(rs,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 16
      call els_manage_one_t(i,evec_thr(16,1),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_angles(rs,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
  if (use_BOND(3).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(17).EQV..true.)) then
      do rs=1,nseq
        call force_impropers(rs,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 17
      call els_manage_one_t(i,evec_thr(17,1),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_impropers(rs,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
  if (use_BOND(4).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(18).EQV..true.)) then
      do rs=1,nseq
        call force_torsions(rs,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 18
      call els_manage_one_t(i,evec_thr(18,1),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_torsions(rs,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
  if (use_BOND(5).EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(20).EQV..true.)) then
      do rs=1,nseq
        call force_cmap(rs,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 20
      call els_manage_one_t(i,evec_thr(20,1),t_ca_f,aone)
    else
      do rs=1,nseq
        call force_cmap(rs,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
  if (use_CORR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(2).EQV..true.)) then
      do rs=1,nseq
        call force_corrector(rs,evec_thr(:,1),hmjam%ca_f)
      end do
      i = 2
      call els_manage_one(i,evec_thr(2,1),cart_f,aone)
    else
      do rs=1,nseq
        call force_corrector(rs,evec_thr(:,1),cart_f)
      end do
    end if
  end if
!
  if (use_TOR.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(7).EQV..true.)) then
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec_thr(:,1),i,hmjam%t_ca_f)
        end if
      end do
      i = 7
      call els_manage_one_t(i,evec_thr(7,1),t_ca_f,aone)
    else
      do i=1,nseq
        if (par_TOR2(i).gt.0) then
          call force_torrs(evec_thr(:,1),i,t_ca_f)
        end if
      end do
    end if
  end if
!
  if ((do_accelsim.EQV..true.).AND.((hmjam%isin(12).EQV..true.).OR.(hmjam%isin(2).EQV..true.))) then
    do rs=1,nseq
      call force_boundary_rs(rs,evec_thr(:,1),1,hmjam%t_ca_f)
    end do
    if (hmjam%isin(12).EQV..true.) then
      i = 12
      call els_manage_one_t(i,evec_thr(12,1),t_ca_f,aone)
    end if
    if (hmjam%isin(2).EQV..true.) then
      i = 2
      call els_manage_one_t(i,evec_thr(11,1),t_ca_f,aone)
    end if
  else
    do rs=1,nseq
      call force_boundary_rs(rs,evec_thr(:,1),1,t_ca_f)
    end do
  end if
!
  if (use_ZSEC.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(8).EQV..true.)) then
      do imol=1,nmol
        call force_zsec_gl(imol,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 8
      call els_manage_one_t(i,evec_thr(8,1),t_ca_f,aone)
    else
      do imol=1,nmol
        call force_zsec_gl(imol,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
!       call System_Clock(t2)
!        write(*,*) 'boundary, zsec, tor',t2-t1
!
!       call System_Clock(t1)
!
!
  if (use_DSSP.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(19).EQV..true.)) then
      call force_dssp_gl(evec_thr(:,1),hmjam%t_ca_f)
      i = 19
      call els_manage_one_t(i,evec_thr(19,1),t_ca_f,aone)
    else
      call force_dssp_gl(evec_thr(:,1),t_ca_f)
    end if
  end if
!
  if (use_EMICRO.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(21).EQV..true.)) then
      call force_emicro_gl(evec_thr(:,1),hmjam%t_ca_f,azero)
      i = 21
      call els_manage_one_t(i,evec_thr(21,1),t_ca_f,aone)
    else
      call force_emicro_gl(evec_thr(:,1),t_ca_f,azero)
    end if
  end if
!
  if (use_POLY.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
      do imol=1,nmol
        call force_poly_gl(imol,evec_thr(:,1),hmjam%t_ca_f)
      end do
      i = 14
      call els_manage_one_t(i,evec_thr(14,1),t_ca_f,aone)
    else
      do imol=1,nmol
        call force_poly_gl(imol,evec_thr(:,1),t_ca_f)
      end do
    end if
  end if
!
  if (use_DREST.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(10).EQV..true.)) then
      call force_drest(evec_thr(:,1),hmjam%t_ca_f,azero)
      i = 10
      call els_manage_one_t(i,evec_thr(10,1),t_ca_f,aone)
    else
      call force_drest(evec_thr(:,1),t_ca_f,azero)
    end if
  end if
!
  cart_f(1:n,1) = cart_f(1:n,1) + t_ca_f(1,1:n)
  cart_f(1:n,2) = cart_f(1:n,2) + t_ca_f(2,1:n)
  cart_f(1:n,3) = cart_f(1:n,3) + t_ca_f(3,1:n)
  if ((freshnbl.EQV..true.).OR.(lrel_md.eq.2)) then
    cart_f_tr(1:n,:) = cart_f_tr(1:n,:) + transpose(t_ca_f_tr(:,1:n))
    evec_lr(:) = evec_thr(:,3)
  end if
!  write(*,*) cart_f_tr(120,:)
!
!  write(*,*) invtemp*sum(evec_thr(:,1))
  if ((freshnbl.EQV..true.)) evec_tr(:) = evec_thr(:,2)
  if (freshnbl.EQV..true.) nbl_meannbs = nint(1.0*sum(rs_nbl(1:nseq)%nnbs)/(1.0*nseq))
  evec(:) = evec_thr(:,1) + evec_tr(:) + evec_lr(:)
  cart_f(:,:) = cart_f(:,:) + cart_f_tr(:,:)
!  write(*,*) cart_f(120,:)
!
  force3 = sum(evec(:))
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(11).EQV..true.)) then
    hmjam%ca_f(:,:) = cart_f(:,:)
    edum = force3
    i = 11
    call els_manage_one(i,edum,cart_f,atwo)
    evec(11) = hmjam%boosts(11,1)
    force3 = force3 + evec(11)
  end if
!
  return
!
end
!
!------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
! the same as force3 for threaded execution
!
! the routine starts with an initialization and transfer block that uses static limits
! for the large arrays
!
! this is followed by (roughly) the computation of nonbonded energies and forces (using dynamic
! load balancing) and computation of other terms (bias and bonded potentials)
!
! the routine finishes by assembling net forces and energies from the thread-specific arrays
!
! we must enter the subroutine from inside a parallel region, which implies that stack
! variables are going to be thread-private by default
! 
subroutine force3_threads(evec,evec_tr,evec_lr,forceflag,totale) !,tpi,tpn)
!
  use system
  use energies
  use mcsums
  use forces
  use cutoffs
  use atoms
  use iounit
  use sequen
  use fyoc
  use threads
  use movesets
  use molecule
  use tabpot
  use wl
  use inter, ONLY: nrsintra,nrpolintra
!
  implicit none
!
  logical, INTENT(IN):: forceflag
  RTYPE, INTENT(INOUT):: evec(MAXENERGYTERMS),evec_tr(MAXENERGYTERMS),evec_lr(MAXENERGYTERMS)
  RTYPE, INTENT(OUT):: totale
!
  integer azero,aone,atwo,athree,afour,asix,rs,i,j,tpi,tpi2,OMP_GET_THREAD_NUM,sta,sto,imol,rslo,rshi,tpp,bdums(19),ixes(6)
  logical sayno,sayyes,freshnbl,pflag,fnbl2,have_rf,need_mod,didelse
  integer(KIND=8) ttimer,ttimer2
  RTYPE evec_thr(MAXENERGYTERMS,3),rdum
!
  tpi = omp_get_thread_num() + 1
  evec_thr(:,:) = 0.0
  azero = 0
  aone = 1
  atwo = 2
  athree = 3
  afour = 4
  asix = 6
  sayno = .false.
  sayyes = .true.
  freshnbl = .false.
  didelse = .false.
  if ((mod(nstep,nbl_up).eq.0).OR.(nstep.le.1).OR.(forceflag.EQV..true.)) then
    freshnbl = .true.
  end if
  if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    pflag = .true.
  else
    pflag = .false.
  end if
  call OMP_SET_NESTED(.true.)
!
!$OMP SINGLE
  totale = 0.0
  bnd_fr = 0.0
  bnd_f(:,:) = 0.0
  evec(:) = 0.0
  if (freshnbl.EQV..true.) then
    evec_tr(:) = 0.0
    evec_lr(:) = 0.0
  else
    if (lrel_md.eq.2) evec_lr(:) = 0.0
    evec_thr(:,3) = evec_lr(:)
    if (lrel_md.ne.2) evec_thr(:,2) = evec_tr(:)
  end if
!$OMP END SINGLE NOWAIT
!
  do i=1,thrdat%maxn
    thr_ca_f(:,thr_limits(1,tpi):thr_limits(2,tpi),i) = 0.0
  end do
  cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) = 0.0
  if (do_accelsim.EQV..true.) then
    t_ca_f(:,thr_limits(1,tpi):thr_limits(2,tpi)) = 0.0
    hmjam%evec_thr(:,tpi) = 0.0
  end if
! the mid- and long-range interaction arrays (energies and cartesian forces) are only updated
! every so many steps
  if ((freshnbl.EQV..true.)) then
    if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) rs_nbl(thr_limits(3,tpi):thr_limits(4,tpi))%ntmpanb = 0
    cart_f_tr(thr_limits(1,tpi):thr_limits(2,tpi),:) = 0.0
    do i=1,thrdat%maxn
      thr_ca_f_tr(:,thr_limits(1,tpi):thr_limits(2,tpi),i) = 0.0
    end do
    if ((use_POLAR.EQV..true.).AND.(use_IMPSOLV.EQV..true.)) then
      do i=1,thrdat%maxn
        thr_sum_scrcbs_tr(thr_limits(1,tpi):thr_limits(2,tpi),i) = 0.0
      end do
      sum_scrcbs_tr(thr_limits(1,tpi):thr_limits(2,tpi)) = 0.0
      sum_scrcbs_lr(thr_limits(1,tpi):thr_limits(2,tpi)) = 0.0
    end if
  else
    if (lrel_md.eq.2) then
      cart_f_tr(thr_limits(1,tpi):thr_limits(2,tpi),:) = 0.0
      do i=1,thrdat%maxn
        thr_ca_f_tr(:,thr_limits(1,tpi):thr_limits(2,tpi),i) = 0.0
      end do
    end if
  end if

  if (use_IMPSOLV.EQV..true.) then
    sisa(thr_limits(1,tpi):thr_limits(2,tpi))%nix = 0
    sav_dr(:,thr_limits(1,tpi):thr_limits(2,tpi)) = 0.0
    do i=1,thrdat%maxn
      thr_sisa(thr_limits(1,tpi):thr_limits(2,tpi),i)%nix = 0
      thr_svte(thr_limits(1,tpi):thr_limits(2,tpi),i) = 0.0
    end do
    if (use_POLAR.EQV..true.) then
      sum_scrcbs(thr_limits(1,tpi):thr_limits(2,tpi)) = 0.0
      do i=1,thrdat%maxn
        thr_sum_scrcbs(thr_limits(1,tpi):thr_limits(2,tpi),i) = 0.0
      end do
    end if
  end if
! transfer
  atinfo(1,thr_limits(1,tpi):thr_limits(2,tpi)) = x(thr_limits(1,tpi):thr_limits(2,tpi))
  atinfo(2,thr_limits(1,tpi):thr_limits(2,tpi)) = y(thr_limits(1,tpi):thr_limits(2,tpi))
  atinfo(3,thr_limits(1,tpi):thr_limits(2,tpi)) = z(thr_limits(1,tpi):thr_limits(2,tpi))
  if (use_POLAR.EQV..true.) atinfo(4,thr_limits(1,tpi):thr_limits(2,tpi)) = atq(thr_limits(1,tpi):thr_limits(2,tpi))
  if (use_IMPSOLV.EQV..true.) then
    if (par_IMPSOLV2(1).eq.1) then
      atsavinfo(1,thr_limits(1,tpi):thr_limits(2,tpi)) = atr(thr_limits(1,tpi):thr_limits(2,tpi))
      atsavinfo(2,thr_limits(1,tpi):thr_limits(2,tpi)) = atsavred(thr_limits(1,tpi):thr_limits(2,tpi))*&
 &                                                     atvol(thr_limits(1,tpi):thr_limits(2,tpi))
    else if (par_IMPSOLV2(1).eq.2) then
      atsavinfo(1,thr_limits(1,tpi):thr_limits(2,tpi)) = atr(thr_limits(1,tpi):thr_limits(2,tpi))
      atsavinfo(2,thr_limits(1,tpi):thr_limits(2,tpi)) = atsavred(thr_limits(1,tpi):thr_limits(2,tpi))
      atsavinfo(3,thr_limits(1,tpi):thr_limits(2,tpi)) = atvol(thr_limits(1,tpi):thr_limits(2,tpi))
    end if
  end if
!
!$OMP BARRIER
!
! all_respairs_nbl performs its own load balancing and ends with a barrier
  if ((use_cutoffs.EQV..true.).AND.(ideal_run.EQV..false.).AND.(freshnbl.EQV..true.)) then
    if (tpi.eq.1) call System_Clock(count=ttimer)
    call all_respairs_nbl(skip_frz,is_tab,tpi)
    if (tpi.eq.1) then
      call System_Clock(count=ttimer2)
      time_energy = time_energy + 1.0*(ttimer-ttimer2)
      time_nbl = time_nbl + 1.0*(ttimer2-ttimer)
    end if
  end if
!
! general note on crosslinks: because crosslinks generally require short-range adjustments, crosslinked
! pairs do not show up in neighbor lists (rs_nbl or rs_vec)
! consequently, there must always be separate calls for crosslinks when relying on rs_nbl and rs_vec
! to obtain interacting pairs
!
  if ((is_plj.EQV..true.).OR.(is_prflj.EQV..true.).OR.(is_pewlj.EQV..true.).OR.&
 &    (is_fegplj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
!
    if (is_pewlj.EQV..true.) then
      have_rf = .false.
    else 
      have_rf = .true. ! only used if either Ewald or RF are present
    end if
    if ((is_prflj.EQV..true.).OR.(is_pewlj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
      need_mod = .true.
    else
      need_mod = .false.
    end if
!
    if (use_cutoffs.EQV..true.) then
!
      if ((is_pewlj.EQV..true.).AND.(thrdat%subnrs(1).gt.0)) then
!       create a sub-team
!$OMP MASTER
        write(*,*) 'PME by ',tpi
#ifdef DISABLE_OPENMP4
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(thrdat%subnrs(1)) PRIVATE(tpi2,fnbl2) REDUCTION(+:evec_thr)
#else
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(thrdat%subnrs(1)) PRIVATE(tpi2,fnbl2) REDUCTION(+:evec_thr) PROC_BIND(CLOSE)
#endif
        fnbl2 = freshnbl
        tpi2 = OMP_GET_THREAD_NUM()+1
!       only evec_thr and thr_ca_f_tr will be modified, sum_scrcbs is not used in any way
        call force_P_LR(evec_thr(:,3),thr_ca_f_tr(:,:,tpi2),sum_scrcbs,fnbl2,tpi2)
!$OMP END PARALLEL
!$OMP END MASTER ! SINGLE NOWAIT
!       and another
!$OMP SINGLE
        tpp = thrdat%maxn-thrdat%subnrs(1)
#ifdef DISABLE_OPENMP4
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(tpp) PRIVATE(i,tpi2,sta,sto) REDUCTION(+:evec_thr)
#else
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(tpp) PRIVATE(i,tpi2,sta,sto) REDUCTION(+:evec_thr) PROC_BIND(SPREAD)
#endif
        tpi2 = OMP_GET_THREAD_NUM()+1
!       load balancing for the nonbonded loops
        if (thr_dlb(3,1).gt.0) then
          if (tpi2.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
          call System_Clock(count=ttimer)
          thr_timings(5,tpi2) = thr_timings(5,tpi2) + ttimer
        end if
        if (use_waterloops.EQV..true.) then
          do i=thr_limits(31,tpi2),thr_limits(32,tpi2) ! 1,rsw1-1
            if ((freshnbl.EQV..true.)) then
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
              end if
              if (rs_nbl(i)%nwnbtrs.gt.0) then
                call Vforce_PLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,afour)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
            end if
            if (rs_nbl(i)%nwnbs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,athree)
            end if
          end do
          do i=thr_limits(33,tpi2),thr_limits(34,tpi2) ! rsw1,nseq
            if ((freshnbl.EQV..true.)) then
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
              end if
              if (rs_nbl(i)%nwnbtrs.gt.0) then
                if (is_3site.EQV..true.) then
                  call Vforce_PLJ_C_T3P(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
                else if (is_4site.EQV..true.) then
                  call Vforce_PLJ_C_T4P(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
                else if (is_5site.EQV..true.) then
                  call Vforce_PLJ_C_T5P(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nwnbtrs,rs_nbl(i)%nwnbtrats,atwo)
                end if
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
            end if
            if (rs_nbl(i)%nwnbs.gt.0) then
              if (is_3site.EQV..true.) then
                call Vforce_PLJ_C_T3P(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
              else if (is_4site.EQV..true.) then
                call Vforce_PLJ_C_T4P(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
              else if (is_5site.EQV..true.) then
                call Vforce_PLJ_C_T5P(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nwnbs,rs_nbl(i)%nwnbats,aone)
              end if
            end if
          end do
        else
          do i=thr_limits(9,tpi2),thr_limits(10,tpi2) ! 1,nseq
            if ((freshnbl.EQV..true.)) then
              if (rs_nbl(i)%nnbtrs.gt.0) then
                call Vforce_PLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi2),aone,rs_nbl(i)%nnbtrs,rs_nbl(i)%nnbtrats,atwo)
              end if
            end if
            if (rs_nbl(i)%nnbs.gt.0) then
              call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi2),aone,rs_nbl(i)%nnbs,rs_nbl(i)%nnbats,aone)
            end if
          end do
        end if
!       residue self-interactions and sequence-neighbor interactions
        if (use_waterloops.EQV..true.) then
!         the loop over water molecules itself (limits 33:34) is mostly taken care of already since:
!            1) there are no intra-residue interactions (this was checked) -> only EW/RF corrections missing
!            2) sequence neighbors were added to NB list
          sta = thr_limits(31,tpi2)
          sto = thr_limits(32,tpi2)-1
        else
          sta = thr_limits(9,tpi2)
          sto = thr_limits(10,tpi2)-1
        end if
        do i=sta,sto
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa) ! thread-safe since i=i
          call force_rsp_long_mod(evec_thr(:,1),i,i,sayno,cart_f)
          if (rsp_vec(i).eq.1) then
            call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa) ! thread-safe since i+1 always <= thr_limits(10,tpi2)
            call force_rsp_long_mod(evec_thr(:,1),i,i+1,sayno,cart_f)
          else if ((rsp_vec(i).eq.2).AND.(freshnbl.EQV..true.)) then
            call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
            call force_rsp_long_mod(evec_thr(:,2),i,i+1,sayno,cart_f_tr)
          end if
          if (disulf(i).gt.i) then
            call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
            call force_rsp_crlk_mod(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi))
          end if
        end do
        if (use_waterloops.EQV..true.) then
          do i=thr_limits(33,tpi2),thr_limits(34,tpi2)
            call force_rsp_long_mod(evec_thr(:,1),i,i,sayno,cart_f)
          end do
        end if
        if ((sto+1).ge.sta) then
          i = sto+1
          call Vforce_rsp(evec_thr(:,1),i,i,sayno,cart_f,sisa)
          call force_rsp_long_mod(evec_thr(:,1),i,i,sayno,cart_f)
          if (disulf(i).gt.i) then
            call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
            call force_rsp_crlk_mod(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi))
          end if
        end if
        if (thr_dlb(3,1).gt.0) then
          call System_Clock(count=ttimer)
          thr_timings(6,tpi2) = thr_timings(6,tpi2) + ttimer  
        end if
!$OMP BARRIER
        if ((sto+1).ge.sta) then
!$OMP CRITICAL(NB_FORCE_SNBOVER)
          if (rsp_vec(i).eq.1) then
            call Vforce_rsp(evec_thr(:,1),i,i+1,sayno,cart_f,sisa) 
            call force_rsp_long_mod(evec_thr(:,1),i,i+1,sayno,cart_f)
          else if ((rsp_vec(i).eq.2).AND.(freshnbl.EQV..true.)) then
            call Vforce_rsp(evec_thr(:,2),i,i+1,sayno,cart_f_tr,sisa)
            call force_rsp_long_mod(evec_thr(:,2),i,i+1,sayno,cart_f_tr)
          end if
!$OMP END CRITICAL(NB_FORCE_SNBOVER)
        end if
!$OMP END PARALLEL
!$OMP END SINGLE NOWAIT
!
      else
        if (freshnbl.EQV..true.) then
          ixes(1) = 31
          ixes(2) = 32
          ixes(3) = 73
          ixes(4) = 74
          ixes(5) = 75
          ixes(6) = 76
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
!         load balancing for the nonbonded loops with TR update
          if (thr_dlb(13,1).gt.0) then
            if (tpi.eq.1) thr_dlb(13,2) = thr_dlb(13,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(29,tpi) = thr_timings(29,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(16) = rs_nbl(i)%nnbats
              bdums(9) = 1
              bdums(10) = rs_nbl(i)%ngnbs
              bdums(17) = rs_nbl(i)%ngnbats
              bdums(11) = 1
              bdums(12) = rs_nbl(i)%nwnbs
              bdums(18) = rs_nbl(i)%nwnbats
              bdums(1) = 1
              bdums(2) = rs_nbl(i)%nnbtrs
              bdums(13) = rs_nbl(i)%nnbtrats
              bdums(3) = 1
              bdums(4) = rs_nbl(i)%ngnbtrs
              bdums(14) = rs_nbl(i)%ngnbtrats
              bdums(5) = 1
              bdums(6) = rs_nbl(i)%nwnbtrs
              bdums(15) = rs_nbl(i)%nwnbtrats
              bdums(19) = 1
            end if
            if (bdums(2).ge.bdums(1)) then
              call Vforce_PLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(1),bdums(2),bdums(13),atwo)
            end if
            if (bdums(4).ge.bdums(3)) then
              call Vforce_FEGPLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(3),bdums(4),bdums(14),atwo)
            end if
            if ((use_waterloops.EQV..true.).AND.(i.ge.rsw1)) then
              if (bdums(6).ge.bdums(5)) then
                if (is_3site.EQV..true.) then
                  call Vforce_PLJ_C_T3P(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(5),bdums(6),bdums(15),atwo)
                else if (is_4site.EQV..true.) then
                  call Vforce_PLJ_C_T4P(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(5),bdums(6),bdums(15),atwo)
                else if (is_5site.EQV..true.) then
                  call Vforce_PLJ_C_T5P(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(5),bdums(6),bdums(15),atwo)
                end if
              end if
            else
              if (bdums(6).ge.bdums(5)) then
                call Vforce_PLJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(5),bdums(6),bdums(15),afour)
              end if
            end if
            if (bdums(8).ge.bdums(7)) then
              call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(7),bdums(8),bdums(16),aone)
            end if
              if (bdums(10).ge.bdums(9)) then
              call Vforce_FEGPLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(9),bdums(10),bdums(17),aone)
            end if
            if ((use_waterloops.EQV..true.).AND.(i.ge.rsw1)) then
              if ((bdums(19).eq.1).AND.(need_mod.EQV..true.)) call force_rsp_long_mod2(evec_thr(:,1),i,i,have_rf,thr_ca_f(:,:,tpi))
              if (bdums(12).ge.bdums(11)) then
                if (is_3site.EQV..true.) then
                  call Vforce_PLJ_C_T3P(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),aone)
                else if (is_4site.EQV..true.) then
                  call Vforce_PLJ_C_T4P(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),aone)
                else if (is_5site.EQV..true.) then
                  call Vforce_PLJ_C_T5P(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),aone)
                end if
              end if
            else
              if (bdums(12).ge.bdums(11)) then
                call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),athree)
              end if
            end if
            if ((i.lt.rsw1).AND.(bdums(19).eq.1)) then
              if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
              else if (rsp_vec(i).eq.2) then
                call Vforce_rsp2(evec_thr(:,2),i,i+1,sayno,thr_ca_f_tr(:,:,tpi),svte,sisa)
              end if
              if (need_mod.EQV..true.) then
                call force_rsp_long_mod2(evec_thr(:,1),i,i,have_rf,thr_ca_f(:,:,tpi))
                if (rsp_vec(i).eq.1) then
                  call force_rsp_long_mod2(evec_thr(:,1),i,i+1,have_rf,thr_ca_f(:,:,tpi))
                else if (rsp_vec(i).eq.2) then
                  call force_rsp_long_mod2(evec_thr(:,2),i,i+1,have_rf,thr_ca_f_tr(:,:,tpi))
                end if
              else
!               the check against nrpolintra is legitimate only because we know use_TABUL to be false in this branch
                if (nrpolintra(i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
                else if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long2(evec_thr(:,2),i,i+1,thr_ca_f_tr(:,:,tpi),sum_scrcbs)
                end if
              end if
              if (disulf(i).gt.i) then
                call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
                if (need_mod.EQV..true.) then
                  call force_rsp_crlk_mod(evec_thr(:,1),i,disulf(i),have_rf,thr_ca_f(:,:,tpi))
                else
                  call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),sum_scrcbs)
                end if
              end if
            end if
          end do
          if (thr_dlb(13,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(30,tpi) = thr_timings(30,tpi) + ttimer  
          end if
        else
          ixes(1) = 9
          ixes(2) = 10
          ixes(3) = 69
          ixes(4) = 70
          ixes(5) = 71
          ixes(6) = 72
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
!         load balancing for the nonbonded loops without TR update
          if (thr_dlb(3,1).gt.0) then
            if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(16) = rs_nbl(i)%nnbats
              bdums(9) = 1
              bdums(10) = rs_nbl(i)%ngnbs
              bdums(17) = rs_nbl(i)%ngnbats
              bdums(11) = 1
              bdums(12) = rs_nbl(i)%nwnbs
              bdums(18) = rs_nbl(i)%nwnbats
              bdums(19) = 1
            end if
            if (bdums(8).ge.bdums(7)) then
              call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(7),bdums(8),bdums(16),aone)
            end if
            if (bdums(10).ge.bdums(9)) then
              call Vforce_FEGPLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(9),bdums(10),bdums(17),aone)
            end if
            if ((use_waterloops.EQV..true.).AND.(i.ge.rsw1)) then
              if ((bdums(19).eq.1).AND.(need_mod.EQV..true.)) call force_rsp_long_mod2(evec_thr(:,1),i,i,have_rf,thr_ca_f(:,:,tpi))
              if (bdums(12).ge.bdums(11)) then
                if (is_3site.EQV..true.) then
                  call Vforce_PLJ_C_T3P(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),aone)
                else if (is_4site.EQV..true.) then
                  call Vforce_PLJ_C_T4P(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),aone)
                else if (is_5site.EQV..true.) then
                  call Vforce_PLJ_C_T5P(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),aone)
                end if
              end if
            else
              if (bdums(12).ge.bdums(11)) then
                call Vforce_PLJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(11),bdums(12),bdums(18),athree)
              end if
            end if
            if ((bdums(19).eq.1).AND.(i.lt.rsw1)) then
              if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
              if (rsp_vec(i).eq.1) call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
              if (need_mod.EQV..true.) then
                call force_rsp_long_mod2(evec_thr(:,1),i,i,have_rf,thr_ca_f(:,:,tpi))
                if (rsp_vec(i).eq.1) call force_rsp_long_mod2(evec_thr(:,1),i,i+1,have_rf,thr_ca_f(:,:,tpi))
              else
!               the check against nrpolintra is legal only because we know use_TABUL is false in this branch
                if (nrpolintra(i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
                if (rsp_vec(i).eq.1) call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
              end if
              if (disulf(i).gt.i) then
                call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
                if (need_mod.EQV..true.) then
                  call force_rsp_crlk_mod(evec_thr(:,1),i,disulf(i),have_rf,thr_ca_f(:,:,tpi))
                else
                  call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),sum_scrcbs)
                end if
              end if
            end if
          end do
          if (thr_dlb(3,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer  
          end if
        end if ! whether freshnbl is true or not (twin-range update is performed or not)
        if (mod(nstep,HUGE(nstep)).eq.0) then ! activate for debugging by replacing HUGE(nstep) with something reasonable
!$OMP BARRIER
!$OMP MASTER
          call check_thread_loop_bounds2(freshnbl,ixes)
!$OMP END MASTER
        end if
      end if ! whether a separate team is doing reciprocal sum or not
!     long-range electrostatics treatment, this uses its own load balancing counter
!     we need a barrier for lrel_md being 4 or 5 since the ranges are not the same and thr_ca_f_tr may already have been used above
!     (for Ewald, twin-range is always empty, and for RF nothing else is incremented)
      if (is_pewlj.EQV..true.) then
!       only evec_thr and thr_ca_f_tr will be modified, sum_scrcbs is not used in any way
        if (thrdat%subnrs(1).eq.0) call force_P_LR(evec_thr(:,3),thr_ca_f_tr(:,:,tpi),sum_scrcbs,freshnbl,tpi)
      else if ((is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
!       nothing happens except that evec_thr may be modified by RF constant
        call force_P_LR(evec_thr(:,3),thr_ca_f_tr(:,:,tpi),sum_scrcbs,freshnbl,tpi)
      else
        if (lrel_md.ne.1) then
!$OMP BARRIER
        end if
!       sum_scrcbs is still not used
        call force_P_LR(evec_thr(:,3),thr_ca_f_tr(:,:,tpi),sum_scrcbs,freshnbl,tpi)
      end if

    else
!
!     no cutoffs: this is a debugging or, most likely, nonsensical calculation
      if (thr_dlb(3,1).gt.0) then
        if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
        call System_Clock(count=ttimer)
        thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
      end if
      do i=thr_limits(9,tpi),thr_limits(10,tpi)
        if (i.lt.(nseq-1)) call Vforce_PLJ(evec_thr(:,1),i,thr_ca_f(:,:,tpi))
        if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
!       the check against nrpolintra is legal only because we know use_TABUL is false in this branch
        if (nrpolintra(i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
        if (i.lt.nseq) then
          call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
          call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
        end if
        if (disulf(i).gt.i) then
          call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
          call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),sum_scrcbs) ! sum_scrcbs unused
        end if
      end do
      if (thr_dlb(3,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer  
      end if
    end if
!
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(13).EQV..true.)) then
      i = 13
      hmjam%evec_thr(i,tpi) = sum(evec_thr)
!     implies initial barrier
      call els_manage_one_t_threads(i,atwo,evec_thr(13,1),tpi) ! everything so far is the total nonbonded force
    end if
!
! the block for pure 6/12 interactions
!
  else if ((is_lj.EQV..true.).OR.(is_ev.EQV..true.)) then
!
    if (use_cutoffs.EQV..true.) then
!
      if (freshnbl.EQV..true.) then
        ixes(1) = 31
        ixes(2) = 32
        ixes(3) = 73
        ixes(4) = 74
        ixes(5) = 75
        ixes(6) = 76
        if (thr_limits(ixes(4),tpi).le.0) then
          call get_thread_loop_bounds2(freshnbl,ixes,tpi)
        else
          call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
        end if
!       load balancing for the nonbonded loops with TR update
        if (thr_dlb(13,1).gt.0) then
          if (tpi.eq.1) thr_dlb(13,2) = thr_dlb(13,2) + 1
          call System_Clock(count=ttimer)
          thr_timings(29,tpi) = thr_timings(29,tpi) + ttimer
        end if
        do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
          if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
            call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
          else
            bdums(7) = 1
            bdums(8) = rs_nbl(i)%nnbs
            bdums(16) = rs_nbl(i)%nnbats
            bdums(1) = 1
            bdums(2) = rs_nbl(i)%nnbtrs
            bdums(13) = rs_nbl(i)%nnbtrats
            bdums(19) = 1
          end if
          if (bdums(2).ge.bdums(1)) then
            call Vforce_LJ_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(1),bdums(2),bdums(13),atwo)
          end if
          if (bdums(8).ge.bdums(7)) then
            call Vforce_LJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(7),bdums(8),bdums(16),aone)
          end if
          if (bdums(19).eq.1) then
            if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
            if (rsp_vec(i).eq.1) then
              call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
            else if (rsp_vec(i).eq.2) then
              call Vforce_rsp2(evec_thr(:,2),i,i+1,sayno,thr_ca_f_tr(:,:,tpi),svte,sisa)
            end if
            if (disulf(i).gt.i) then
              call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
            end if
          end if
        end do
        if (thr_dlb(13,1).gt.0) then
          call System_Clock(count=ttimer)
          thr_timings(30,tpi) = thr_timings(30,tpi) + ttimer  
        end if
      else
        ixes(1) = 9
        ixes(2) = 10
        ixes(3) = 69
        ixes(4) = 70
        ixes(5) = 71
        ixes(6) = 72
        if (thr_limits(ixes(4),tpi).le.0) then
          call get_thread_loop_bounds2(freshnbl,ixes,tpi)
        else
          call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
        end if
!       load balancing for the nonbonded loops without TR update
        if (thr_dlb(3,1).gt.0) then
          if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
          call System_Clock(count=ttimer)
          thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
        end if
        do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
          if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
            call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
          else
            bdums(7) = 1
            bdums(8) = rs_nbl(i)%nnbs
            bdums(16) = rs_nbl(i)%nnbats
            bdums(19) = 1
          end if
          if (bdums(8).ge.bdums(7)) then
            call Vforce_LJ_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(7),bdums(8),bdums(16),aone)
          end if
           if (bdums(19).eq.1) then
            if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
            if (rsp_vec(i).eq.1) call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
            if (disulf(i).gt.i) then
              call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
            end if
          end if
        end do
        if (thr_dlb(3,1).gt.0) then
          call System_Clock(count=ttimer)
          thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer  
        end if
      end if ! whether freshnbl is true or not (twin-range update is performed or not)
    else
!     no cutoffs: this is a debugging or, most likely, nonsensical calculation
      if (thr_dlb(3,1).gt.0) then
        if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
        call System_Clock(count=ttimer)
        thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
      end if
      do i=thr_limits(9,tpi),thr_limits(10,tpi)
        if (i.lt.(nseq-1)) call Vforce_LJ(evec_thr(:,1),i,thr_ca_f(:,:,tpi))
        if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
        if (i.lt.nseq) call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
        if (disulf(i).gt.i) then
          call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa) ! args unused mostly
        end if
      end do
      if (thr_dlb(3,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer  
      end if
    end if
!
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(13).EQV..true.)) then
      i = 13
      hmjam%evec_thr(i,tpi) = sum(evec_thr)
!     implies initial barrier
      call els_manage_one_t_threads(i,atwo,evec_thr(13,1),tpi) ! everything so far is the total nonbonded force
    end if
!
! the block for pure tabulated interactions
!
  else if (is_tab.EQV..true.) then
!
    if (use_cutoffs.EQV..true.) then
!
      if (thr_dlb(3,1).gt.0) then
        if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
        call System_Clock(count=ttimer)
        thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
      end if
      do i=thr_limits(9,tpi),thr_limits(10,tpi) ! 1,nseq
        if ((freshnbl.EQV..true.)) then
          if (rs_nbl(i)%nnbtrs.gt.0) then
            imol = 0
            do j=1,rs_nbl(i)%nnbtrs
              rslo = min(rs_nbl(i)%nbtr(j),i)
              rshi = max(rs_nbl(i)%nbtr(j),i)
              imol = imol + tbp%rsmat(rshi,rslo) - tbp%rsmat(rslo,rshi) + 1
            end do
            if (imol.gt.0) call Vforce_TAB_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),aone,rs_nbl(i)%nnbtrs,imol,atwo)
          end if
        end if
        if (rs_nbl(i)%nnbs.gt.0) then
          imol = 0
          do j=1,rs_nbl(i)%nnbs
            rslo = min(rs_nbl(i)%nb(j),i)
            rshi = max(rs_nbl(i)%nb(j),i)
            imol = imol + tbp%rsmat(rshi,rslo) - tbp%rsmat(rslo,rshi) + 1
          end do
          if (imol.gt.0) call Vforce_TAB_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),aone,rs_nbl(i)%nnbs,imol,aone)
        end if
      end do
      do i=thr_limits(9,tpi),thr_limits(10,tpi)
!       the check against tbp is legal only because we know use_POLAR is false in this branch
        if (tbp%rsmat(i,i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
        if (rsp_vec(i).eq.1) then
          call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
        else if ((rsp_vec(i).eq.2).AND.(freshnbl.EQV..true.)) then
          call Vforce_rsp_long2(evec_thr(:,2),i,i+1,thr_ca_f_tr(:,:,tpi),sum_scrcbs)
        end if
        if (disulf(i).gt.i) then
          call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),sum_scrcbs) ! sum_scrcbs unused
        end if
      end do
      if (thr_dlb(3,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer  
      end if
!
    else
!     no cutoffs: this is a debugging or, most likely, nonsensical calculation
      if (thr_dlb(3,1).gt.0) then
        if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
        call System_Clock(count=ttimer)
        thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
      end if
      do i=thr_limits(9,tpi),thr_limits(10,tpi)
        if (rs_nbl(i)%ntabnbs.gt.0) then
          call Vforce_TAB(evec_thr(:,1),i,thr_ca_f(:,:,tpi))
        end if
!       the check against tbp is legal only because we know use_POLAR is false in this branch
        if (tbp%rsmat(i,i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
        if (i.lt.nseq) call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
      end do
      if (thr_dlb(3,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer  
      end if
    end if
!
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(13).EQV..true.)) then
      i = 13
      hmjam%evec_thr(i,tpi) = sum(evec_thr)
      call els_manage_one_t_threads(i,atwo,evec_thr(13,1),tpi) ! everything so far is the total nonbonded force
    end if
!
!
  else if (ideal_run.EQV..true.) then
!
!   do nothing
!
  else
!
!   the standard calculation with support for arbitrary Hamiltonians
!   we'll do some subdistinctions for specialized routines, but keep the framework general 
!   any one-shot Hamiltonians should probably be separated out into a construct as above,
!   since distances will be calculated twice ...
!
    didelse = .true.
!
    if (use_IMPSOLV.EQV..true.) then
      call init_svte_threads(athree,tpi)
!$OMP BARRIER
    end if
!
    if (use_cutoffs.EQV..true.) then
!
!     first the short-range terms (ATTLJ, IPP, WCA, FOS)
      if ((is_impljp.EQV..true.).OR.(is_implj.EQV..true.)) then
!       this is the specialized set of fxns for LJ 6/12 with ABSINTH implicit solvent model
        if (freshnbl.EQV..true.) then
          ixes(1) = 31
          ixes(2) = 32
          ixes(3) = 73
          ixes(4) = 74
          ixes(5) = 75
          ixes(6) = 76
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
!         load balancing for the nonbonded loops with TR update
          if (thr_dlb(13,1).gt.0) then
            if (tpi.eq.1) thr_dlb(13,2) = thr_dlb(13,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(29,tpi) = thr_timings(29,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(16) = rs_nbl(i)%nnbats
              bdums(1) = 1
              bdums(2) = rs_nbl(i)%nnbtrs
              bdums(13) = rs_nbl(i)%nnbtrats
              bdums(19) = 1
            end if
            if (bdums(2).ge.bdums(1)) then
              call Vforce_LJIMP_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi),bdums(1),bdums(2),&
 &bdums(13),atwo)
            end if
            if (bdums(8).ge.bdums(7)) then
              call Vforce_LJIMP_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi),bdums(7),bdums(8),&
 &bdums(16),aone)
            end if
            if (bdums(19).eq.1) then
              if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              else if (rsp_vec(i).eq.2) then
                call Vforce_rsp2(evec_thr(:,2),i,i+1,sayno,thr_ca_f_tr(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end if
              if (disulf(i).gt.i) then
                call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end if
            end if
          end do
          if (thr_dlb(13,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(30,tpi) = thr_timings(30,tpi) + ttimer  
          end if
        else
          ixes(1) = 9
          ixes(2) = 10
          ixes(3) = 69
          ixes(4) = 70
          ixes(5) = 71
          ixes(6) = 72
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
          if (thr_dlb(3,1).gt.0) then
            if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(16) = rs_nbl(i)%nnbats
              bdums(19) = 1
            end if
            if (bdums(8).ge.bdums(7)) then
              call Vforce_LJIMP_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi),bdums(7),bdums(8),&
 &bdums(16),aone)
            end if
            if (bdums(19).eq.1) then
              if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end if
              if (disulf(i).gt.i) then
                call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end if
            end if
          end do
          if (thr_dlb(3,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer
          end if
        end if
!
      else if ((use_IPP.EQV..false.).AND.(use_ATTLJ.EQV..false.).AND.(use_IMPSOLV.EQV..false.).AND.(use_WCA.EQV..false.)) then
!
!       do nothing
!
      else
!
!       the general case
        if (freshnbl.EQV..true.) then
          ixes(1) = 31
          ixes(2) = 32
          ixes(3) = 73
          ixes(4) = 74
          ixes(5) = 75
          ixes(6) = 76
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
!         load balancing for the nonbonded loops with TR update
          if (thr_dlb(13,1).gt.0) then
            if (tpi.eq.1) thr_dlb(13,2) = thr_dlb(13,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(29,tpi) = thr_timings(29,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(1) = 1
              bdums(2) = rs_nbl(i)%nnbtrs
              bdums(3) = 1
              bdums(4) = rs_nbl(i)%ngnbtrs
              bdums(5) = 1
              bdums(6) = rs_nbl(i)%nwnbtrs
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(9) = 1
              bdums(10) = rs_nbl(i)%ngnbs
              bdums(11) = 1
              bdums(12) = rs_nbl(i)%nwnbs
              bdums(19) = 1
            end if
            if (use_IMPSOLV.EQV..true.) then
              do j=bdums(1),bdums(2)
                call Vforce_rsp2(evec_thr(:,2),i,rs_nbl(i)%nbtr(j),sayno,thr_ca_f_tr(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(3),bdums(4)
                call Vforce_rsp2(evec_thr(:,2),i,rs_nbl(i)%gnbtr(j),sayno,thr_ca_f_tr(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(5),bdums(6)
                call Vforce_rsp2(evec_thr(:,2),i,rs_nbl(i)%wnbtr(j),sayno,thr_ca_f_tr(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(7),bdums(8)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%nb(j),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              if (bdums(19).eq.1) then
                if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                else if (rsp_vec(i).eq.2) then
                  call Vforce_rsp2(evec_thr(:,2),i,i+1,sayno,thr_ca_f_tr(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                end if
              end if
            else
              do j=bdums(1),bdums(2)
                call Vforce_rsp2(evec_thr(:,2),i,rs_nbl(i)%nbtr(j),sayno,thr_ca_f_tr(:,:,tpi),svte,sisa)
              end do
              do j=bdums(3),bdums(4)
                call Vforce_rsp2(evec_thr(:,2),i,rs_nbl(i)%gnbtr(j),sayno,thr_ca_f_tr(:,:,tpi),svte,sisa)
              end do
              do j=bdums(5),bdums(6)
                call Vforce_rsp2(evec_thr(:,2),i,rs_nbl(i)%wnbtr(j),sayno,thr_ca_f_tr(:,:,tpi),svte,sisa)
              end do
              do j=bdums(7),bdums(8)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%nb(j),sayno,thr_ca_f(:,:,tpi),svte,sisa)
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),sayno,thr_ca_f(:,:,tpi),svte,sisa)
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),sayno,thr_ca_f(:,:,tpi),svte,sisa)
              end do
              if (bdums(19).eq.1) then
                if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
                else if (rsp_vec(i).eq.2) then
                  call Vforce_rsp2(evec_thr(:,2),i,i+1,sayno,thr_ca_f_tr(:,:,tpi),svte,sisa)
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa)
                end if
              end if
            end if
          end do
          if (thr_dlb(13,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(30,tpi) = thr_timings(30,tpi) + ttimer  
          end if
        else
          ixes(1) = 9
          ixes(2) = 10
          ixes(3) = 69
          ixes(4) = 70
          ixes(5) = 71
          ixes(6) = 72
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
          if (thr_dlb(3,1).gt.0) then
            if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(9) = 1
              bdums(10) = rs_nbl(i)%ngnbs
              bdums(11) = 1
              bdums(12) = rs_nbl(i)%nwnbs
              bdums(19) = 1
            end if
            if (use_IMPSOLV.EQV..true.) then
              do j=bdums(7),bdums(8)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%nb(j),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
              end do
              if (bdums(19).eq.1) then
                if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
                end if
              end if
            else
              do j=bdums(7),bdums(8)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%nb(j),sayno,thr_ca_f(:,:,tpi),svte,sisa)
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),sayno,thr_ca_f(:,:,tpi),svte,sisa)
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),sayno,thr_ca_f(:,:,tpi),svte,sisa)
              end do
              if (bdums(19).eq.1) then
                if (nrsintra(i).gt.0) call Vforce_rsp2(evec_thr(:,1),i,i,sayno,thr_ca_f(:,:,tpi),svte,sisa)
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp2(evec_thr(:,1),i,i+1,sayno,thr_ca_f(:,:,tpi),svte,sisa)
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp2(evec_thr(:,1),i,disulf(i),sayno,thr_ca_f(:,:,tpi),svte,sisa)
                end if
              end if
            end if
          end do
          if (thr_dlb(3,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer
          end if
        end if
!
      end if
!
!     setup the screened charge terms (which have the same short-range
!     as the sisa%dr terms unless large charge groups are used with coupled models:
!     in that latter case larger range is always accounted for fully(!))
!$OMP BARRIER
      if (use_IMPSOLV.EQV..true.) then
        call init_svte_threads(asix,tpi)
        call Vforce_dsavdr(sav_dr,tpi)
        if (use_POLAR.EQV..true.) then
!$OMP BARRIER
          call force_setup_scrqs2(tpi)
!$OMP BARRIER
        end if
      end if
!
      if (is_implj.EQV..true.) then
!
!       do nothing here
!
      else if (is_impljp.EQV..true.) then
!
!       this is the specialized set of fxns for POLAR with ABSINTH implicit solvent model
!       note that in the twin-range increment of sum_scrcbs are
!       not incremented further on account of charge-charge interactions
!       between particles relatively far apart, i.e., it's a pseudo-cutoff
!       since the atom will still be close to one of the two charges 
        if (freshnbl.EQV..true.) then
          ixes(1) = 33
          ixes(2) = 34
          ixes(3) = 77
          ixes(4) = 78
          ixes(5) = 79
          ixes(6) = 80
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
!         load balancing for the nonbonded loops with TR update
          if (thr_dlb(14,1).gt.0) then
            if (tpi.eq.1) thr_dlb(14,2) = thr_dlb(14,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(31,tpi) = thr_timings(31,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(16) = rs_nbl(i)%nnbats
              bdums(1) = 1
              bdums(2) = rs_nbl(i)%nnbtrs
              bdums(13) = rs_nbl(i)%nnbtrats
              bdums(19) = 1
            end if
            if (bdums(2).ge.bdums(1)) then
              if ((scrq_model.le.2).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
                call Vforce_PSCRM1256_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi),bdums(1),bdums(2),bdums(13),&
 &atwo)
              else if (scrq_model.eq.4) then
                call Vforce_PSCRM4_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),bdums(1),bdums(2),bdums(13),atwo)
              else if ((scrq_model.eq.3).OR.((scrq_model.ge.7).AND.(scrq_model.le.9))) then
                call Vforce_PSCRM3789_C(evec_thr(:,2),i,thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi),bdums(1),bdums(2),bdums(13),&
 &atwo)
              end if
            end if
            if (bdums(8).ge.bdums(7)) then
              if ((scrq_model.le.2).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
                call Vforce_PSCRM1256_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi),bdums(7),bdums(8),bdums(16),aone)
              else if (scrq_model.eq.4) then
                call Vforce_PSCRM4_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(7),bdums(8),bdums(16),aone)
              else if ((scrq_model.eq.3).OR.((scrq_model.ge.7).AND.(scrq_model.le.9))) then
                call Vforce_PSCRM3789_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi),bdums(7),bdums(8),bdums(16),aone)
              end if
            end if
            if (bdums(19).eq.1) then
!             the check against nrpolintra is legal only because we know use_TABUL is false in this branch
              if (nrpolintra(i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              else if (rsp_vec(i).eq.2) then
                call Vforce_rsp_long2(evec_thr(:,2),i,i+1,thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi))
              end if
              if (disulf(i).gt.i) then
                call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end if
            end if
          end do
          if (thr_dlb(14,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(32,tpi) = thr_timings(32,tpi) + ttimer  
          end if
        else
          ixes(1) = 41
          ixes(2) = 42
          ixes(3) = 81
          ixes(4) = 82
          ixes(5) = 83
          ixes(6) = 84
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
          if (thr_dlb(7,1).gt.0) then
            if (tpi.eq.1) thr_dlb(7,2) = thr_dlb(7,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(17,tpi) = thr_timings(17,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(16) = rs_nbl(i)%nnbats
              bdums(19) = 1
            end if
            if (bdums(8).ge.bdums(7)) then
              if ((scrq_model.le.2).OR.((scrq_model.ge.5).AND.(scrq_model.le.6))) then
                call Vforce_PSCRM1256_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi),bdums(7),bdums(8),bdums(16),aone)
              else if (scrq_model.eq.4) then
                call Vforce_PSCRM4_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),bdums(7),bdums(8),bdums(16),aone)
              else if ((scrq_model.eq.3).OR.((scrq_model.ge.7).AND.(scrq_model.le.9))) then
                call Vforce_PSCRM3789_C(evec_thr(:,1),i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi),bdums(7),bdums(8),bdums(16),aone)
              end if
            end if
            if (bdums(19).eq.1) then
!             the check against nrpolintra is legal only because we know use_TABUL is false in this branch
              if (nrpolintra(i).gt.0) call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              if (rsp_vec(i).eq.1) then
                call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end if
              if (disulf(i).gt.i) then
                call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end if
            end if
          end do
          if (thr_dlb(7,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(18,tpi) = thr_timings(18,tpi) + ttimer
          end if
        end if
!
      else if ((use_TABUL.EQV..false.).AND.(use_POLAR.EQV..false.)) then
!
!       do nothing
!
      else
!
!       the general case
        if (freshnbl.EQV..true.) then
          ixes(1) = 33
          ixes(2) = 34
          ixes(3) = 77
          ixes(4) = 78
          ixes(5) = 79
          ixes(6) = 80
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
!         load balancing for the nonbonded loops with TR update
          if (thr_dlb(14,1).gt.0) then
            if (tpi.eq.1) thr_dlb(14,2) = thr_dlb(14,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(31,tpi) = thr_timings(31,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(1) = 1
              bdums(2) = rs_nbl(i)%nnbtrs
              bdums(3) = 1
              bdums(4) = rs_nbl(i)%ngnbtrs
              bdums(5) = 1
              bdums(6) = rs_nbl(i)%nwnbtrs
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(9) = 1
              bdums(10) = rs_nbl(i)%ngnbs
              bdums(11) = 1
              bdums(12) = rs_nbl(i)%nwnbs
              bdums(19) = 1
            end if
            if ((use_POLAR.EQV..true.).AND.(use_IMPSOLV.EQV..true.)) then
              do j=bdums(1),bdums(2)
                call Vforce_rsp_long2(evec_thr(:,2),i,rs_nbl(i)%nbtr(j),thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi))
              end do
              do j=bdums(3),bdums(4)
                call Vforce_rsp_long2(evec_thr(:,2),i,rs_nbl(i)%gnbtr(j),thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi))
              end do
              do j=bdums(5),bdums(6)
                call Vforce_rsp_long2(evec_thr(:,2),i,rs_nbl(i)%wnbtr(j),thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi))
              end do
              do j=bdums(7),bdums(8)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%nb(j),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end do
              if (bdums(19).eq.1) then
                call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
                else if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long2(evec_thr(:,2),i,i+1,thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi))
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
                end if
              end if
            else
              do j=bdums(1),bdums(2)
                call Vforce_rsp_long2(evec_thr(:,2),i,rs_nbl(i)%nbtr(j),thr_ca_f_tr(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(3),bdums(4)
                call Vforce_rsp_long2(evec_thr(:,2),i,rs_nbl(i)%gnbtr(j),thr_ca_f_tr(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(5),bdums(6)
                call Vforce_rsp_long2(evec_thr(:,2),i,rs_nbl(i)%wnbtr(j),thr_ca_f_tr(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(7),bdums(8)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%nb(j),thr_ca_f(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),thr_ca_f(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),thr_ca_f(:,:,tpi),sum_scrcbs)
              end do
              if (bdums(19).eq.1) then
                call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
                else if (rsp_vec(i).eq.2) then
                  call Vforce_rsp_long2(evec_thr(:,2),i,i+1,thr_ca_f_tr(:,:,tpi),sum_scrcbs)
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),sum_scrcbs)
                end if
              end if
            end if
          end do
          if (thr_dlb(14,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(32,tpi) = thr_timings(32,tpi) + ttimer  
          end if
        else
          ixes(1) = 41
          ixes(2) = 42
          ixes(3) = 81
          ixes(4) = 82
          ixes(5) = 83
          ixes(6) = 84
          if (thr_limits(ixes(4),tpi).le.0) then
            call get_thread_loop_bounds2(freshnbl,ixes,tpi)
          else
            call fix_thread_loop_bounds2(freshnbl,ixes,tpi)
          end if
          if (thr_dlb(7,1).gt.0) then
            if (tpi.eq.1) thr_dlb(7,2) = thr_dlb(7,2) + 1
            call System_Clock(count=ttimer)
            thr_timings(17,tpi) = thr_timings(17,tpi) + ttimer
          end if
          do i=thr_limits(ixes(1),tpi),thr_limits(ixes(2),tpi)
            if ((i.eq.thr_limits(ixes(1),tpi)).OR.(i.eq.thr_limits(ixes(2),tpi))) then
              call get_thread_loop_inst(tpi,freshnbl,i,ixes,bdums)
            else
              bdums(7) = 1
              bdums(8) = rs_nbl(i)%nnbs
              bdums(9) = 1
              bdums(10) = rs_nbl(i)%ngnbs
              bdums(11) = 1
              bdums(12) = rs_nbl(i)%nwnbs
              bdums(19) = 1
            end if
            if ((use_POLAR.EQV..true.).AND.(use_IMPSOLV.EQV..true.)) then
              do j=bdums(7),bdums(8)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%nb(j),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
              end do
              if (bdums(19).eq.1) then
                call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
                end if
              end if
            else
              do j=bdums(7),bdums(8)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%nb(j),thr_ca_f(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(9),bdums(10)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%gnb(j),thr_ca_f(:,:,tpi),sum_scrcbs)
              end do
              do j=bdums(11),bdums(12)
                call Vforce_rsp_long2(evec_thr(:,1),i,rs_nbl(i)%wnb(j),thr_ca_f(:,:,tpi),sum_scrcbs)
              end do
              if (bdums(19).eq.1) then
                call Vforce_rsp_long2(evec_thr(:,1),i,i,thr_ca_f(:,:,tpi),sum_scrcbs)
                if (rsp_vec(i).eq.1) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,i+1,thr_ca_f(:,:,tpi),sum_scrcbs)
                end if
                if (disulf(i).gt.i) then
                  call Vforce_rsp_long2(evec_thr(:,1),i,disulf(i),thr_ca_f(:,:,tpi),sum_scrcbs)
                end if
              end if
            end if
          end do
          if (thr_dlb(7,1).gt.0) then
            call System_Clock(count=ttimer)
            thr_timings(18,tpi) = thr_timings(18,tpi) + ttimer
          end if
        end if
      end if
!
!     finally, add long-range electrostatics corrections (note that these include
!     the next neighbor problem explicitly)
!     remember that IMPSOLV and FEG are mutually exclusive, if that ever changes
!     the construct inside the fxn below has to be changed as well!
      if (use_POLAR.EQV..true.) then
        call force_P_LR(evec_thr(:,3),thr_ca_f_tr(:,:,tpi),thr_sum_scrcbs_tr(:,tpi),freshnbl,tpi)
      end if
!
!   this the final drop-through scenario that operates similar to force1(), it offers only coarse load balancing
!   at the residue level (mostly for debugging)
!
    else
!
      if ((use_IMPSOLV.EQV..true.).OR.(use_IPP.EQV..true.).OR.(use_attLJ.EQV..true.).OR.(use_WCA.EQV..true.)) then
        if (thr_dlb(3,1).gt.0) then
          if (tpi.eq.1) thr_dlb(3,2) = thr_dlb(3,2) + 1
          call System_Clock(count=ttimer)
          thr_timings(5,tpi) = thr_timings(5,tpi) + ttimer
        end if
        if (use_IMPSOLV.EQV..true.) then
          do i=thr_limits(9,tpi),thr_limits(10,tpi) ! 1,nseq
            do j=i,nseq
              call Vforce_rsp2(evec_thr(:,1),i,j,sayno,thr_ca_f(:,:,tpi),thr_svte(:,tpi),thr_sisa(:,tpi))
            end do
          end do
        else
          do i=thr_limits(9,tpi),thr_limits(10,tpi) ! 1,nseq
            do j=i,nseq
              call Vforce_rsp2(evec_thr(:,1),i,j,sayno,thr_ca_f(:,:,tpi),svte,sisa)
            end do
          end do
        end if
        if (thr_dlb(3,1).gt.0) then
          call System_Clock(count=ttimer)
          thr_timings(6,tpi) = thr_timings(6,tpi) + ttimer
        end if
      end if
!
!$OMP BARRIER
      if (use_IMPSOLV.EQV..true.) then
        call init_svte_threads(asix,tpi)
        call Vforce_dsavdr(sav_dr,tpi)
        if (use_POLAR.EQV..true.) then
!$OMP BARRIER
          call force_setup_scrqs2(tpi)
!$OMP BARRIER
        end if
      end if
!
      if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
        if (thr_dlb(7,1).gt.0) then
          if (tpi.eq.1) thr_dlb(7,2) = thr_dlb(7,2) + 1
          call System_Clock(count=ttimer)
          thr_timings(17,tpi) = thr_timings(17,tpi) + ttimer
        end if
        if ((use_POLAR.EQV..true.).AND.(use_IMPSOLV.EQV..true.)) then
          do i=thr_limits(41,tpi),thr_limits(42,tpi) ! 1,nseq
            do j=i,nseq
              call Vforce_rsp_long2(evec_thr(:,1),i,j,thr_ca_f(:,:,tpi),thr_sum_scrcbs(:,tpi))
            end do
          end do
        else
          do i=thr_limits(41,tpi),thr_limits(42,tpi) ! 1,nseq
            do j=i,nseq
              call Vforce_rsp_long2(evec_thr(:,1),i,j,thr_ca_f(:,:,tpi),sum_scrcbs)
            end do
          end do
        end if
        if (thr_dlb(7,1).gt.0) then
          call System_Clock(count=ttimer)
          thr_timings(18,tpi) = thr_timings(18,tpi) + ttimer
        end if
      end if
!
    end if
!
  end if
!
!$OMP BARRIER
  if ((use_IMPSOLV.EQV..true.).AND.(use_POLAR.EQV..true.).AND.(scrq_model.ne.4)) then
    do i=1,thrdat%maxn
      sum_scrcbs(thr_limits(1,tpi):thr_limits(2,tpi)) = sum_scrcbs(thr_limits(1,tpi):thr_limits(2,tpi)) + &
 &                                                      thr_sum_scrcbs(thr_limits(1,tpi):thr_limits(2,tpi),i)
      if (freshnbl.EQV..true.) sum_scrcbs_tr(thr_limits(1,tpi):thr_limits(2,tpi)) = &
 & sum_scrcbs_tr(thr_limits(1,tpi):thr_limits(2,tpi)) + thr_sum_scrcbs_tr(thr_limits(1,tpi):thr_limits(2,tpi),i)
    end do
!$OMP BARRIER
  end if
!
  if ((do_accelsim.EQV..true.).AND.((hmjam%isin(4).EQV..true.).OR.(hmjam%isin(7).EQV..true.).OR.(hmjam%isin(12).EQV..true.).OR.&
 &    (hmjam%isin(13).EQV..true.).OR.(hmjam%isin(15).EQV..true.).OR.(hmjam%isin(16).EQV..true.).OR.(hmjam%isin(17).EQV..true.).OR.&
 &    (hmjam%isin(18).EQV..true.).OR.(hmjam%isin(20).EQV..true.).OR.(hmjam%isin(2).EQV..true.))) then
!   rest of polar force must be dealt with first
    if ((use_POLAR.EQV..true.).AND.(use_IMPSOLV.EQV..true.)) then
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call Vforce_dcbdscrq(rs,thr_ca_f(:,:,tpi)) ! just updates force, requires sisq_dr and all sum_scrbs
      end do
    end if
    if ((didelse.EQV..true.).AND.(hmjam%isin(13).EQV..true.)) then
      i = 13
      hmjam%evec_thr(i,tpi) = sum(evec_thr)
      call els_manage_one_t_threads(i,atwo,evec_thr(13,1),tpi) ! everything so far is the total nonbonded force
    end if
!   load balancing for residue-based potentials with no data dependency
!   forces can overlap, therefore threadprivate arrays must be used
    if (thr_dlb(6,1).gt.0) then
      if (tpi.eq.1) thr_dlb(6,2) = thr_dlb(6,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(13,tpi) = thr_timings(13,tpi) + ttimer
    end if
    do rs=thr_limits(37,tpi),thr_limits(38,tpi)
      if ((hmjam%isin(15).EQV..false.).AND.(use_BOND(1).EQV..true.)) call force_bonds(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if ((hmjam%isin(16).EQV..false.).AND.(use_BOND(2).EQV..true.)) call force_angles(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if ((hmjam%isin(17).EQV..false.).AND.(use_BOND(3).EQV..true.)) call force_impropers(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if ((hmjam%isin(18).EQV..false.).AND.(use_BOND(4).EQV..true.)) call force_torsions(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if ((hmjam%isin(20).EQV..false.).AND.(use_BOND(5).EQV..true.)) call force_cmap(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if ((hmjam%isin(4).EQV..false.).AND.(use_IMPSOLV.EQV..true.)) then
        call force_freesolv(rs,evec_thr(:,1),thr_ca_f(:,:,tpi)) ! updates f/e, requires sav_dr/sisa%dr
      end if
      if ((hmjam%isin(12).EQV..false.).AND.(hmjam%isin(2).EQV..false.)) &
 &       call force_boundary_rs(rs,evec_thr(:,1),aone,thr_ca_f(:,:,tpi))
      if ((hmjam%isin(7).EQV..false.).AND.(use_TOR.EQV..true.)) then
        if (par_TOR2(rs).gt.0) call force_torrs(evec_thr(:,1),rs,thr_ca_f(:,:,tpi))
      end if
    end do
    if (thr_dlb(6,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(14,tpi) = thr_timings(14,tpi) + ttimer
    end if
    i = 0
    call els_manage_one_t_threads(i,azero,rdum,tpi) ! move entries in thr_ca_f to safety
    if ((use_IMPSOLV.EQV..true.).AND.(hmjam%isin(4).EQV..true.)) then
      i = 4
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_freesolv(rs,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(4,1),tpi)
    end if
    if (hmjam%isin(7).EQV..true.) then
      i = 7
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        if (par_TOR2(rs).gt.0) call force_torrs(hmjam%evec_thr(:,tpi),rs,thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(7,1),tpi)
    end if
    if ((hmjam%isin(12).EQV..true.).OR.(hmjam%isin(2).EQV..true.)) then
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_boundary_rs(rs,hmjam%evec_thr(:,tpi),aone,thr_ca_f(:,:,tpi))
      end do
      if (hmjam%isin(2).EQV..false.) then
        if (use_OSMO.EQV..true.) evec_thr(11,1) = evec_thr(11,1) + hmjam%evec_thr(11,tpi)
        i = 12
        call els_manage_one_t_threads(i,aone,evec_thr(12,1),tpi)
      else if (hmjam%isin(12).EQV..false.) then
        evec_thr(12,1) = evec_thr(12,1) + hmjam%evec_thr(12,tpi)
        i = 2
        call els_manage_one_t_threads(i,aone,evec_thr(11,1),tpi)
      else
        i = 2
        call els_manage_one_t_threads(i,aone,evec_thr(11,1),tpi)
        i = 12
        call els_manage_one_t_threads(i,aone,evec_thr(12,1),tpi)
      end if
    end if
    if ((use_BOND(1).EQV..true.).AND.(hmjam%isin(15).EQV..true.)) then
      i = 15
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_bonds(rs,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(15,1),tpi)
    end if
    if ((use_BOND(2).EQV..true.).AND.(hmjam%isin(16).EQV..true.)) then
      i = 16
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_angles(rs,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(16,1),tpi)
    end if
    if ((use_BOND(3).EQV..true.).AND.(hmjam%isin(17).EQV..true.)) then
      i = 17
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_impropers(rs,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(17,1),tpi)
    end if
    if ((use_BOND(4).EQV..true.).AND.(hmjam%isin(18).EQV..true.)) then
      i = 18
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_torsions(rs,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(18,1),tpi)
    end if
    if ((use_BOND(5).EQV..true.).AND.(hmjam%isin(20).EQV..true.)) then
      i = 20
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call force_cmap(rs,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(20,1),tpi)
    end if
  else
!   load balancing for residue-based potentials with no data dependency
!   forces can overlap, therefore threadprivate arrays must be used
    if (thr_dlb(6,1).gt.0) then
      if (tpi.eq.1) thr_dlb(6,2) = thr_dlb(6,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(13,tpi) = thr_timings(13,tpi) + ttimer
    end if
    do rs=thr_limits(37,tpi),thr_limits(38,tpi)
      if (use_BOND(1).EQV..true.) call force_bonds(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if (use_BOND(2).EQV..true.) call force_angles(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if (use_BOND(3).EQV..true.) call force_impropers(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if (use_BOND(4).EQV..true.) call force_torsions(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if (use_BOND(5).EQV..true.) call force_cmap(rs,evec_thr(:,1),thr_ca_f(:,:,tpi))
      if (use_IMPSOLV.EQV..true.) then
        call force_freesolv(rs,evec_thr(:,1),thr_ca_f(:,:,tpi)) ! updates f/e, requires sav_dr/sisa%dr
        if (use_POLAR.EQV..true.) then
          call Vforce_dcbdscrq(rs,thr_ca_f(:,:,tpi)) ! just updates force, requires sisq_dr and 
        end if
      end if
      call force_boundary_rs(rs,evec_thr(:,1),aone,thr_ca_f(:,:,tpi))
      if (use_TOR.EQV..true.) then
        if (par_TOR2(rs).gt.0) call force_torrs(evec_thr(:,1),rs,thr_ca_f(:,:,tpi))
      end if
    end do
    if (thr_dlb(6,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(14,tpi) = thr_timings(14,tpi) + ttimer
    end if
  end if
!
  if ((do_accelsim.EQV..true.).AND.((hmjam%isin(19).EQV..true.).OR.(hmjam%isin(8).EQV..true.))) then
    if ((hmjam%isin(19).EQV..true.).AND.(use_DSSP.EQV..true.)) then
      i = 19
      call els_manage_one_t_threads(i,azero,rdum,tpi)
!$OMP SINGLE
      call force_dssp_gl(hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
!$OMP END SINGLE
      call els_manage_one_t_threads(i,aone,evec_thr(19,1),tpi)
    else if (use_DSSP.EQV..true.) then
!$OMP SINGLE
      call force_dssp_gl(evec_thr(:,1),thr_ca_f(:,:,tpi))
!$OMP END SINGLE
    end if
    if ((use_ZSEC.EQV..true.).AND.(use_DSSP.EQV..false.)) then
!$OMP BARRIER
    end if
    if ((hmjam%isin(8).EQV..true.).AND.(use_ZSEC.EQV..true.)) then
      i = 8
      call els_manage_one_t_threads(i,azero,rdum,tpi)
      do imol=thr_limits(35,tpi),thr_limits(36,tpi)
        call force_zsec_gl(imol,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      end do
      call els_manage_one_t_threads(i,aone,evec_thr(8,1),tpi)
    else if (use_ZSEC.EQV..true.) then
      do imol=thr_limits(35,tpi),thr_limits(36,tpi)
        call force_zsec_gl(imol,evec_thr(:,1),thr_ca_f(:,:,tpi))
      end do
    end if
  else if ((use_ZSEC.EQV..true.).OR.(use_DSSP.EQV..true.)) then
    if (use_ZSEC.EQV..true.) then
!$OMP BARRIER
    end if
!$OMP SECTIONS
!$OMP SECTION
    if (use_ZSEC.EQV..true.) then
      do imol=1,nmol
        call force_zsec_gl(imol,evec_thr(:,1),thr_ca_f(:,:,tpi))
      end do
    end if
!$OMP SECTION
    if (use_DSSP.EQV..true.) then
      call force_dssp_gl(evec_thr(:,1),thr_ca_f(:,:,tpi))
    end if
!$OMP END SECTIONS
  end if
!
! the rest are cutoff-independent terms
  if (use_DREST.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(10).EQV..true.)) then
      i = 10
      call els_manage_one_t_threads(i,azero,rdum,tpi)
      call force_drest(hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi),tpi)
      call els_manage_one_t_threads(i,aone,evec_thr(10,1),tpi)
    else
      call force_drest(evec_thr(:,1),thr_ca_f(:,:,tpi),tpi)
    end if
  end if
!
  if (use_EMICRO.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(21).EQV..true.)) then
      i = 21
      call els_manage_one_t_threads(i,azero,rdum,tpi)
      call force_emicro_gl(hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi),tpi)
      call els_manage_one_t_threads(i,aone,evec_thr(21,1),tpi)
    else
      call force_emicro_gl(evec_thr(:,1),thr_ca_f(:,:,tpi),tpi)
    end if
  end if
!
  if (use_POLY.EQV..true.) then
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
      i = 14
      call els_manage_one_t_threads(i,azero,rdum,tpi)
    end if
    i = 1
    if (nmlgs.gt.0) then
      do while ((mlg_limits(i,5,tpi).lt.thr_limits(35,tpi)).AND.(i.lt.nmlgs))
        i = min(i+1,nmlgs)
      end do
    end if
    do imol=thr_limits(35,tpi),thr_limits(36,tpi)
      if (nmlgs.gt.0) then
        if (mlg_limits(i,5,tpi).eq.imol) then
          i = min(i+1,nmlgs)
          cycle
        end if
      end if
      if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
        call force_poly_gl(imol,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi))
      else
        call force_poly_gl(imol,evec_thr(:,1),thr_ca_f(:,:,tpi))
      end if
    end do
    do i=1,nmlgs
      if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
        call force_poly_threads(i,hmjam%evec_thr(:,tpi),thr_ca_f(:,:,tpi),tpi)
      else
        call force_poly_threads(i,evec_thr(:,1),thr_ca_f(:,:,tpi),tpi)
      end if
    end do
    if ((do_accelsim.EQV..true.).AND.(hmjam%isin(14).EQV..true.)) then
      i = 14
      call els_manage_one_t_threads(i,aone,evec_thr(14,1),tpi)
    end if
  end if
  thr_rutil(1:MAXENERGYTERMS,tpi) = evec_thr(:,1)
  if (freshnbl.EQV..true.) thr_rutil((MAXENERGYTERMS+1):(2*MAXENERGYTERMS),tpi) = evec_thr(:,2)
  if ((lrel_md.eq.2).OR.(freshnbl.EQV..true.)) thr_rutil((2*MAXENERGYTERMS+1):(3*MAXENERGYTERMS),tpi) = evec_thr(:,3)
!
!$OMP BARRIER
!
  do i=1,thrdat%maxn
    cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) = cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) + &
 & transpose(thr_ca_f(:,thr_limits(1,tpi):thr_limits(2,tpi),i))
    if ((freshnbl.EQV..true.).OR.(lrel_md.eq.2)) then
      cart_f_tr(thr_limits(1,tpi):thr_limits(2,tpi),:) = cart_f_tr(thr_limits(1,tpi):thr_limits(2,tpi),:) + &
 & transpose(thr_ca_f_tr(:,thr_limits(1,tpi):thr_limits(2,tpi),i))
    end if
  end do
!
!$OMP SINGLE
  if ((freshnbl.EQV..true.).OR.(lrel_md.eq.2)) then
    evec_lr(:) = sum(thr_rutil((2*MAXENERGYTERMS+1):(3*MAXENERGYTERMS),1:thrdat%maxn),dim=2)
    if (freshnbl.EQV..true.) evec_tr(:) = sum(thr_rutil((MAXENERGYTERMS+1):(2*MAXENERGYTERMS),1:thrdat%maxn),dim=2)
  end if
  evec(:) = sum(thr_rutil(1:MAXENERGYTERMS,1:thrdat%maxn),dim=2) + evec_tr(:) + evec_lr(:)
  totale = sum(evec)
!$OMP END SINGLE NOWAIT
!
  if (use_IMPSOLV.EQV..true.) then
    if (mod(nstep,100).eq.0) then
      do i=thr_limits(1,tpi),thr_limits(2,tpi)
        if (sisa(i)%alsz.gt.(sisa(i)%nix+10)) call Vsisa_resize(i,sisa,azero)
      end do
      do i=1,n
        if (thr_sisa(i,tpi)%alsz.gt.(thr_sisa(i,tpi)%nix+10)) call Vsisa_resize(i,thr_sisa(:,tpi),azero)
      end do
    end if
  end if
!  write(*,*) cart_f_tr(120,:)
!
  if (use_cutoffs.EQV..true.) then
    cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) = cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) + &
 &                                             cart_f_tr(thr_limits(1,tpi):thr_limits(2,tpi),:)
  end if
  if (do_accelsim.EQV..true.) then
    cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) = cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) + &
 &                                        transpose(t_ca_f(:,thr_limits(1,tpi):thr_limits(2,tpi)))
  end if
!
  if ((do_accelsim.EQV..true.).AND.(hmjam%isin(11).EQV..true.)) then
    i = 11
    hmjam%evec_thr(:,tpi) = 0.0
    if (tpi.le.1) hmjam%evec_thr(i,tpi) = totale
    call els_manage_one_t_threads(i,athree,evec_thr(11,1),tpi) ! everything so far is the total nonbonded force
    if (tpi.le.1) totale = evec_thr(11,1)
  end if
!$OMP BARRIER
!  if (tpi.le.1) write(*,*) sum(thr_sisa(:,:)%alsz),sum(sisa(:)%alsz)
!  write(*,*) cart_f(120,:)
!
  return
!
end
!
#endif 
!
!------------------------------------------------------------------------
!
! a wrapper which calls the "correct" long-range electrostatics subroutine
! note that "long-range" really stands for multiple, fundamentally different
! things
! note that evec and ca_f are not initialized 
!
subroutine force_P_LR(evec,ca_f,sum_s,freshnbl,tpi)
!
  use energies
  use atoms
  use cutoffs
  use ewalds
  use iounit
  use mcsums
  use movesets
  use sequen
!
  implicit none
!
  integer, INTENT(IN):: tpi
  logical, INTENT(IN):: freshnbl
!
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n)
  RTYPE sum_s(n)
!
  if (use_POLAR.EQV..false.) return
!
  if ((lrel_md.eq.4).OR.(lrel_md.eq.5)) then
    if (freshnbl.EQV..true.) then
      call force_P_LR_NBL(evec,ca_f,skip_frz,sum_s,tpi)
    end if
  else if (lrel_md.eq.2) then
!   the reciprocal space sum must be calculated at every step
    if (ewald_mode.eq.1) then
      call force_pme(evec,ca_f,tpi)
    else if (ewald_mode.eq.2) then
      call force_ewald(evec,ca_f,tpi)
    end if
  else if (lrel_md.eq.3) then
!   do nothing except increment energy by self-term (RF has no LR or quasi-LR terms)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if (freshnbl.EQV..true.) then
      evec(6) = evec(6) + rfcnst
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  else if (lrel_md.eq.1) then
!   really do nothing
  else
    write(ilog,*) 'Fatal. Called force_P_LR(...) with unsupported &
 &LR electrostatics model (offending mode is ',lrel_md,'). Please re&
 &port this bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------------
!
! wrapper
!
subroutine all_crosslink_corrs(evec,cut,ca_f,s_sv,si,sum_s)
!
  use sequen
  use energies
  use forces
  use atoms, ONLY: n
!
  implicit none
!
  logical cut
  RTYPE evec(MAXENERGYTERMS),ca_f(3,n),sum_s(n),s_sv(n)
  type(t_sisa) si(n)
  integer lk
  logical cutt
!
  cutt = cut
  if ((is_pewlj.EQV..true.).OR.(is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
    if ((is_prflj.EQV..true.).OR.(is_fegprflj.EQV..true.)) then
      cutt = .true.
    else
      cutt = .false.
    end if
    do lk=1,n_crosslinks
      call Vforce_rsp2(evec,crosslink(lk)%rsnrs(1),crosslink(lk)%rsnrs(2),cutt,ca_f,s_sv,si)
      call force_rsp_crlk_mod(evec,crosslink(lk)%rsnrs(1),crosslink(lk)%rsnrs(2),cutt,ca_f)
    end do
  else
    do lk=1,n_crosslinks
      call Vforce_rsp2(evec,crosslink(lk)%rsnrs(1),crosslink(lk)%rsnrs(2),cutt,ca_f,s_sv,si)
      call Vforce_rsp_long2(evec,crosslink(lk)%rsnrs(1),crosslink(lk)%rsnrs(2),ca_f,sum_s)
    end do
  end if
!
end
!
!------------------------------------------------------------------------------------
!
! this is a subroutine providing support for foreign energy calculations
! in REMD runs
! care has to be taken such that forces from the original Hamiltonian
! are retained in principle (this routine does NOT deal with swaps per se)
!
subroutine lamforce(rve,fve,tpi)
!
  use energies
  use atoms
  use iounit
  use molecule
  use mpistuff
  use system
  use forces
  use units
  use cutoffs
  use dssps
  use ems
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  RTYPE rve(mpi_nodes),fve(mpi_nodes)
#ifdef ENABLE_MPI
#ifdef ENABLE_THREADS
  integer sta2,sto2
#endif
  integer i,j,jj,k,imol,which,i_start,i_end,sta,sto
  RTYPE force3,evec(MAXENERGYTERMS),evec_tr(MAXENERGYTERMS),evec_lr(MAXENERGYTERMS)
  RTYPE vbu(MAXREDIMS),eva,dum,dum2,evaref
  logical needforcebu,needsavup,needemup,badflg(MAXREDIMS),atrue
  integer vbui(MAXREDIMS)
  RTYPE, ALLOCATABLE:: ifbu(:)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, lamforce(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
#endif
!
#ifdef ENABLE_MPI
!
  needforcebu = .false.
  needsavup = .false.
  needemup = .false.
  atrue = .true.
  dum2 = 0.0
  evaref = esave
  evec_lr(:) = esterms_lr(:)
  evec_tr(:) = esterms_tr(:)
  evec(:) = esterms(:)
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
!$OMP BARRIER
!$OMP MASTER
#else
  sta = 1
  sto = n
#endif
  rve(:) = 0.0
  fve(:) = 0.0
  do i=1,re_conddim
!    if ((re_types(i).eq.7)).OR.&
! &      ((re_types(i).ge.12).AND.(re_types(i).le.19))) then
!      write(ilog,*) 'Fatal. REMD with the chosen condition (',&
! &re_types(i),') is not yet supported. Check back later.'
!      call fexit() 
!    end if
    if (re_types(i).eq.1) vbu(i)=kelvin
    if (re_types(i).eq.2) vbu(i)=scale_IPP
    if (re_types(i).eq.3) vbu(i)=scale_attLJ
    if (re_types(i).eq.4) vbu(i)=scale_WCA
    if (re_types(i).eq.5) vbu(i)=scale_POLAR
    if (re_types(i).eq.6) vbu(i)=scale_IMPSOLV
    if (re_types(i).eq.7) vbu(i)=par_IMPSOLV(2)
    if (re_types(i).eq.8) vbu(i)=scale_TOR
    if (re_types(i).eq.9) vbu(i)=scale_ZSEC
    if (re_types(i).eq.10) vbu(i)=par_ZSEC(1)
    if (re_types(i).eq.11) vbu(i)=par_ZSEC(3)
    if (re_types(i).eq.12) vbui(i)=scrq_model
    if (re_types(i).eq.13) vbu(i)=par_IMPSOLV(3)
    if (re_types(i).eq.14) vbu(i)=par_IMPSOLV(4)
    if (re_types(i).eq.15) vbu(i)=par_IMPSOLV(6)
    if (re_types(i).eq.16) vbu(i)=par_IMPSOLV(7)
    if (re_types(i).eq.17) vbu(i)=par_IMPSOLV(8)
    if (re_types(i).eq.18) vbui(i)=i_sqm
    if (re_types(i).eq.19) vbu(i)=par_IMPSOLV(9)
    if (re_types(i).eq.20) vbu(i)=scale_FEGS(1)
    if (re_types(i).eq.21) vbu(i)=scale_FEGS(3)
    if (re_types(i).eq.22) vbu(i)=scale_FEGS(6)
    if (re_types(i).eq.23) vbu(i)=scale_TABUL
    if (re_types(i).eq.24) vbu(i)=scale_POLY
    if (re_types(i).eq.25) vbu(i)=scale_DREST
    if (re_types(i).eq.26) vbu(i)=scale_FEGS(15)
    if (re_types(i).eq.27) vbu(i)=scale_FEGS(16)
    if (re_types(i).eq.28) vbu(i)=scale_FEGS(17)
    if (re_types(i).eq.29) vbu(i)=scale_FEGS(18)
    if (re_types(i).eq.30) vbu(i)=par_DSSP(9)
    if (re_types(i).eq.31) vbu(i)=par_DSSP(7)
!   if 32 do nothing
    if (re_types(i).eq.33) vbu(i)=scale_EMICRO
    if (re_types(i).eq.34) vbu(i)=emthreshdensity
  end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  if (needforcebu.EQV..true.) then ! this conditional is currently not trigger-able -> untested code
    if (fycxyz.eq.2) then
      cart_a(sta:sto,:) = cart_f(sta:sto,:) ! hi-jack of cart_a
    else if (fycxyz.eq.1) then
      cart_a(sta:sto,:) = cart_f(sta:sto,:)
      k = 0
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        sta2 = thr_limits(61,tpi)
        sto2 = thr_limits(62,tpi)
      else
        sta2 = 1
        sto2 = nmol
      end if
      do imol=sta2,sto2
#else
      do imol=1,nmol
#endif
        k = k + size(dc_di(imol)%f)
      end do
      allocate(ifbu(k))
      jj = 0
#ifdef ENABLE_THREADS
      do imol=sta2,sto2
#else
      do imol=1,nmol
#endif
        do j=1,size(dc_di(imol)%f)
          jj = jj + 1
          ifbu(jj) = dc_di(imol)%f(j)
        end do
      end do
    end if
  end if
!
  if (reol_all.EQV..true.) then
    i_start = 1
    i_end = re_conditions
  else
    i_start = max((myrank+1)-1,1)
    i_end = min((myrank+1)+1,re_conditions)
  end if
  do i=i_start,i_end
!  do i=1,re_conditions
!   set the other conditions (as requested: either neighbors only or all)
!   remember that we don't change the use_XX-flags
!   this implies, however, that for scale_XX = 0.0, we compute
!   a bunch of terms all multiplied by 0.0. in order to preserve
!   this information (to compute derivatives with respect to scale_XX),
!   we therefore use a little detour (set to 1.0, subtract out)
    eva = 0.0
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    do j=1,re_conddim
      badflg(j) = .false.
      if (re_types(j).eq.1) then
        kelvin = re_mat(i,j)
        invtemp = 1.0/(gasconst*kelvin)
      else if (re_types(j).eq.2) then
        scale_IPP   = re_mat(i,j)
        if (scale_IPP.le.0.0) then
          badflg(j) = .true.
          scale_IPP = 1.0
        end if
      else if (re_types(j).eq.3) then
        scale_attLJ = re_mat(i,j)
        if (scale_attLJ.le.0.0) then
          badflg(j) = .true.
          scale_attLJ = 1.0
        end if
      else if (re_types(j).eq.4) then
        scale_WCA   = re_mat(i,j)
        if (scale_WCA.le.0.0) then
          badflg(j) = .true.
          scale_WCA = 1.0
        end if
      else if (re_types(j).eq.5) then
        scale_POLAR = re_mat(i,j)
        if (scale_POLAR.le.0.0) then
          badflg(j) = .true.
          scale_POLAR = 1.0
        end if
      else if (re_types(j).eq.6) then
        scale_IMPSOLV = re_mat(i,j)
        if (scale_IMPSOLV.le.0.0) then
          badflg(j) = .true.
          scale_IMPSOLV = 1.0
        end if
      else if (re_types(j).eq.7) then
!       note that the relevance of this relies entirely on whether use_IMPSOLV is true
        par_IMPSOLV(2) = re_mat(i,j)
        if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
          coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
        else
          coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
        end if
        if (lrel_md.eq.3) then
          call setup_rfcnst2()
        end if
      else if (re_types(j).eq.8) then
        scale_TOR = re_mat(i,j)
        if (scale_TOR.le.0.0) then
          badflg(j) = .true.
          scale_TOR = 1.0
        end if
      else if (re_types(j).eq.9) then
        scale_ZSEC = re_mat(i,j)
        if (scale_ZSEC.le.0.0) then
          badflg(j) = .true.
          scale_ZSEC = 1.0
        end if
      else if (re_types(j).eq.10) then
!       note that the relevance of this hinges on scale_ZSEC
        par_ZSEC(1) = re_mat(i,j)
      else if (re_types(j).eq.11) then
!       note that the relevance of this hinges on scale_ZSEC
        par_ZSEC(3) = re_mat(i,j)
      else if (re_types(j).eq.12) then
!       note that the relevance of this hinges on scale_IMPSOLV
        scrq_model = nint(re_mat(i,j))
        if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
          coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
        else
          coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
        end if
      else if (re_types(j).eq.13) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(3) = re_mat(i,j)
      else if (re_types(j).eq.14) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(4) = re_mat(i,j)
      else if (re_types(j).eq.15) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(6) = re_mat(i,j)
      else if (re_types(j).eq.16) then
!       note that the relevance of this hinges on scale_IMPSOLV
        needsavup = .true.
        par_IMPSOLV(7) = re_mat(i,j)
      else if (re_types(j).eq.17) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        par_IMPSOLV(8) = 1./re_mat(i,j)
      else if (re_types(j).eq.18) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        i_sqm = nint(re_mat(i,j))
      else if (re_types(j).eq.19) then
!       note that the relevance of this hinges on scale_IMPSOLV and scrq_model
        par_IMPSOLV(9) = re_mat(i,j)
      else if (re_types(j).eq.20) then
        scale_FEGS(1) = re_mat(i,j)
        call setup_parFEG(1)
      else if (re_types(j).eq.21) then
        scale_FEGS(3) = re_mat(i,j)
        call setup_parFEG(3)
      else if (re_types(j).eq.22) then
        scale_FEGS(6) = re_mat(i,j)
        call setup_parFEG(6)
        if (lrel_md.eq.3) then
          call setup_rfcnst()
        end if
      else if (re_types(j).eq.23) then
        scale_TABUL = re_mat(i,j)
        if (scale_TABUL.le.0.0) then
          badflg(j) = .true.
          scale_TABUL = 1.0
        end if
      else if (re_types(j).eq.24) then
        scale_POLY = re_mat(i,j)
        if (scale_POLY.le.0.0) then
          badflg(j) = .true.
          scale_POLY = 1.0
        end if
      else if (re_types(j).eq.25) then
        scale_DREST = re_mat(i,j)
        if (scale_DREST.le.0.0) then
          badflg(j) = .true.
          scale_DREST = 1.0
        end if
      else if (re_types(j).eq.26) then
        scale_FEGS(15) = re_mat(i,j)
      else if (re_types(j).eq.27) then
        scale_FEGS(16) = re_mat(i,j)
      else if (re_types(j).eq.28) then
        scale_FEGS(17) = re_mat(i,j)
      else if (re_types(j).eq.29) then
        scale_FEGS(18) = re_mat(i,j)
      else if (re_types(j).eq.30) then
        par_DSSP(9) = re_mat(i,j)
      else if (re_types(j).eq.31) then
        par_DSSP(7) = re_mat(i,j)
!     if 32 do nothing
      else if (re_types(j).eq.33) then
        scale_EMICRO = re_mat(i,j)
        if (scale_EMICRO.le.0.0) then
          badflg(j) = .true.
          scale_EMICRO = 1.0
        end if
      else if (re_types(j).eq.34) then
        emthreshdensity = re_mat(i,j)
        needemup = .true.
      end if
    end do
    if (needsavup.EQV..true.) call absinth_savprm()
    if (needemup.EQV..true.) then
      call scale_emmap(emthreshdensity,dum,dum2)
      call precompute_diff_emmap()
    end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!   now compute the current energy at this condition with
!   our structure
!   note that we correct for energy/structure mismatch for our home node
!   (artifact of the leapfrog algorithm) by re-computing ourselves as well(!)
    if (Tonly.EQV..true.) then
      eva = evaref
    else
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
      else
        esave = force3(esterms,esterms_tr,esterms_lr,atrue)
      end if
#else
      esave = force3(esterms,esterms_tr,esterms_lr,atrue)
#endif
      eva = esave
    end if
!   the computation of derivatives is really only meaningful when
!   noTI is false
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    fve(i) = 0.0
    do j=1,re_conddim
      if (badflg(j).EQV..true.) then
        if (re_types(j).eq.2) then
          fve(i) = fve(i) + esterms(1)
          esterms(1) = 0.0
        else if (re_types(j).eq.3) then
          fve(i) = fve(i) + esterms(3)
          esterms(3) = 0.0
        else if (re_types(j).eq.4) then
          fve(i) = fve(i) + esterms(5)
          esterms(5) = 0.0
        else if (re_types(j).eq.5) then
          fve(i) = fve(i) + esterms(6)
          esterms(6) = 0.0
        else if (re_types(j).eq.6) then
          fve(i) = fve(i) + esterms(4)
          esterms(4) = 0.0
        else if (re_types(j).eq.8) then
          fve(i) = fve(i) + esterms(7)
          esterms(7) = 0.0
        else if (re_types(j).eq.9) then
          fve(i) = fve(i) + esterms(8)
          esterms(8) = 0.0
        else if (re_types(j).eq.23) then
          fve(i) = fve(i) + esterms(9)
          esterms(9) = 0.0
        else if (re_types(j).eq.24) then
          fve(i) = fve(i) + esterms(14)
          esterms(14) = 0.0
        else if (re_types(j).eq.25) then
          fve(i) = fve(i) + esterms(10)
          esterms(10) = 0.0
        else if (re_types(j).eq.33) then
          fve(i) = fve(i) + esterms(21)
          esterms(21) = 0.0
        end if
      else
        if (re_types(j).eq.2) then
          fve(i) = fve(i) + esterms(1)/scale_IPP
        else if (re_types(j).eq.3) then
          fve(i) = fve(i) + esterms(3)/scale_attLJ
        else if (re_types(j).eq.4) then
          fve(i) = fve(i) + esterms(5)/scale_WCA
        else if (re_types(j).eq.5) then
          fve(i) = fve(i) + esterms(6)/scale_POLAR
        else if (re_types(j).eq.6) then
          fve(i) = fve(i) + esterms(4)/scale_IMPSOLV
        else if (re_types(j).eq.8) then
          fve(i) = fve(i) + esterms(7)/scale_TOR
        else if (re_types(j).eq.9) then
          fve(i) = fve(i) + esterms(8)/scale_ZSEC
        else if (re_types(j).eq.10) then
          which = 1
          do imol=1,nmol
            call der_zsec_gl(imol,fve(i),which)
          end do
        else if (re_types(j).eq.11) then
          which = 2
          do imol=1,nmol
            call der_zsec_gl(imol,fve(i),which)
          end do
        else if (re_types(j).eq.23) then
          fve(i) = fve(i) + esterms(9)/scale_TABUL
        else if (re_types(j).eq.24) then
          fve(i) = fve(i) + esterms(14)/scale_POLY
        else if (re_types(j).eq.25) then
          fve(i) = fve(i) + esterms(10)/scale_DREST
!        else if (re_types(j).eq.30) then
!          which = 2
!          do imol=1,nmol
!            call der_dssp_gl(imol,fve(i),which)
!          end do
!        else if (re_types(j).eq.31) then
!          which = 1
!          do imol=1,nmol
!            call der_dssp_gl(imol,fve(i),which)
!          end do
        else if (re_types(j).eq.33) then
          fve(i) = fve(i) + esterms(21)/scale_EMICRO
        end if
      end if
    end do
!   now recover the actual energy at this condition (see above)
    rve(i) = eva*invtemp
    needsavup = .false.
    needemup = .false.
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end do
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  do i=1,re_conddim
    if (re_types(i).eq.1) then 
      kelvin = vbu(i)
      invtemp = 1.0/(gasconst*kelvin)
    else if (re_types(i).eq.2) then
      scale_IPP   = vbu(i)
    else if (re_types(i).eq.3) then
      scale_attLJ = vbu(i)
    else if (re_types(i).eq.4) then
      scale_WCA   = vbu(i)
    else if (re_types(i).eq.5) then
      scale_POLAR = vbu(i)
    else if (re_types(i).eq.6) then
      scale_IMPSOLV = vbu(i)
    else if (re_types(i).eq.7) then
      par_IMPSOLV(2) = vbu(i)
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
      else
        coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
      end if
      if (lrel_md.eq.3) then
        call setup_rfcnst2()
      end if
    else if (re_types(i).eq.8) then
      scale_TOR = vbu(i)
    else if (re_types(i).eq.9) then
      scale_ZSEC = vbu(i)
    else if (re_types(i).eq.10) then
      par_ZSEC(1) = vbu(i)
    else if (re_types(i).eq.11) then
      par_ZSEC(3) = vbu(i)
    else if (re_types(i).eq.12) then
      scrq_model = vbui(i)
      if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
        coul_scr = 1.0 - 1.0/par_IMPSOLV(2)
      else
        coul_scr = 1.0 - 1.0/sqrt(par_IMPSOLV(2))
      end if
    else if (re_types(i).eq.13) then
      needsavup = .true.
      par_IMPSOLV(3) = vbu(i)
    else if (re_types(i).eq.14) then
      needsavup = .true.
      par_IMPSOLV(4) = vbu(i)
    else if (re_types(i).eq.15) then
      needsavup = .true.
      par_IMPSOLV(6) = vbu(i)
    else if (re_types(i).eq.16) then
      needsavup = .true.
      par_IMPSOLV(7) = vbu(i)
    else if (re_types(i).eq.17) then
      par_IMPSOLV(8) = vbu(i)
    else if (re_types(i).eq.18) then
      i_sqm = vbui(i)
    else if (re_types(i).eq.19) then
      par_IMPSOLV(9) = vbu(i)
    else if (re_types(i).eq.20) then
      scale_FEGS(1) = vbu(i)
      call setup_parFEG(1) 
    else if (re_types(i).eq.21) then
      scale_FEGS(3) = vbu(i)
      call setup_parFEG(3)
    else if (re_types(i).eq.22) then
      scale_FEGS(6) = vbu(i)
      call setup_parFEG(6)
      if (lrel_md.eq.3) then
        call setup_rfcnst()
      end if
    else if (re_types(i).eq.23) then
      scale_TABUL = vbu(i)
    else if (re_types(i).eq.24) then
      scale_POLY  = vbu(i)
    else if (re_types(i).eq.25) then
      scale_DREST = vbu(i)
    else if (re_types(i).eq.26) then
      scale_FEGS(15) = vbu(i)
    else if (re_types(i).eq.27) then
      scale_FEGS(16) = vbu(i)
    else if (re_types(i).eq.28) then
      scale_FEGS(17) = vbu(i)
    else if (re_types(i).eq.29) then
      scale_FEGS(18) = vbu(i)
    else if (re_types(i).eq.30) then
      par_DSSP(9) = vbu(i)
    else if (re_types(i).eq.31) then
      par_DSSP(7) = vbu(i)
!   if 32 do nothing
    else if (re_types(i).eq.33) then
      scale_EMICRO = vbu(i)
    else if (re_types(i).eq.34) then
      emthreshdensity = vbu(i)
      needemup = .true.
    end if
  end do
  if (needsavup.EQV..true.) call absinth_savprm()
  if (needemup.EQV..true.) then
    call scale_emmap(emthreshdensity,dum,dum2)
    call precompute_diff_emmap()
  end if
  esave = evaref
  esterms_lr(:) = evec_lr(:)
  esterms_tr(:) = evec_tr(:)
  esterms(:) = evec(:)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  if (needforcebu.EQV..true.) then
    if (fycxyz.eq.2) then
      cart_f(sta:sto,:) = cart_a(sta:sto,:)
    else if (fycxyz.eq.1) then
      cart_f(sta:sto,:) = cart_a(sta:sto,:)
      jj = 0
#ifdef ENABLE_THREADS
      do imol=sta2,sto2
#else
      do imol=1,nmol
#endif
        do j=1,size(dc_di(imol)%f)
          jj = jj + 1
          dc_di(imol)%f(j) = ifbu(jj)
        end do
      end do
      deallocate(ifbu)
    end if
  end if
#else
!
  write(ilog,*) 'Fatal. Called lamforce(...) in non MPI-calculation.&
 & This is most definitely a bug.'
  call fexit()
#endif
!
end
!
!---------------------------------------------------------------------cc
!
! this routine assumes that there is only a single ghosted
! residue with effective scaled charges of par_FEG2(9)
! we also assume this residue is a dipole residue, not a group
! with a net charge
!
subroutine setup_rfcnst()
!
  use energies
  use iounit
  use cutoffs
  use sequen
  use polypep
  use atoms
  use system
  use math
  use units
!
  implicit none
!
  integer rs,i,k
  RTYPE resQ,netQ,istr,t1,t2,kap,epsr
!
  if (use_FEG.EQV..false.) return
!
  istr= 0.0
!
! we can only assume these are homogeneously distributed (see comments in
! grf_setup() in polar.f)
  netQ = 0.0
  k = 0
  do i=1,cglst%ncs
    if (par_FEG(atmres(cglst%it(i))).EQV..true.) then
      k = k + cglst%nc(i)
      netQ = netQ + cglst%tc(i)
    end if
    istr = istr + cglst%tc(i)*cglst%tc(i)
  end do
  if (k.ne.0) then
    write(ilog,*) 'Fatal. FEG in conjunction with (G)RF treatment &
 & is only allowed for net-neutral (sets of) residues.'
    call fexit()
  else if (abs(netQ).gt.1.0e-3) then
    write(ilog,*) 'Warning. Using a (set of) group(s) in FEG that has a net absolute fractional charge  &
 &exceeding 0.001. Check input parameters.'
  end if
  epsr = par_IMPSOLV(2)
  netQ = 0.0
  k = 0
  do rs=1,nseq
    resQ = 0.0
    if (par_FEG3(rs).EQV..true.) then
!     remember that fully de-coupled residues are treated as if they had zero charge
      cycle
    else if (par_FEG(rs).EQV..true.) then
      k = k + 1
      do i=1,at(rs)%npol
        resQ = resQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
      end do
      netQ = netQ + par_FEG2(9)*par_FEG2(9)*resQ
    else
      do i=1,at(rs)%npol
        resQ = resQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
      end do
      netQ = netQ + scale_POLAR*resQ
    end if
  end do
!
! standard RF (GRF limiting case with zero ionic strength)
  if (rf_mode.eq.2) then
    par_POLAR(1) = (epsr-1.0)/((2.0*epsr+1.0)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
! GRF
  else
    kap = sqrt(4.0*PI*electric*invtemp*istr/ens%insV)
    t1 = 1.0 + kap*mcel_cutoff
    t2 = kap*kap*mcel_cutoff2
    par_POLAR(1) = (t1*(epsr-1.0) + 0.5*epsr*t2)/&
 &        (((2.0*epsr+1.0)*t1 + epsr*t2)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
  end if
  rfcnst = -0.5*electric*netQ*par_POLAR(2)
!
end
!
!-----------------------------------------------------------------------
!
subroutine setup_rfcnst2()
!
  use energies
  use iounit
  use cutoffs
  use sequen
  use polypep
  use atoms
  use system
  use math
  use units
!
  implicit none
!
  integer rs,i
  RTYPE netQ,istr,t1,t2,kap,epsr
!
  if (use_FEG.EQV..true.) then
    call setup_rfcnst()
    return
  end if
! 
  istr= 0.0
!
! we can only assume these are homogeneously distributed (see comments in
! grf_setup() in polar.f)
  do i=1,cglst%ncs
    istr = istr + cglst%tc(i)*cglst%tc(i)
  end do
  epsr = par_IMPSOLV(2)
  netQ = 0.0
  do rs=1,nseq
    do i=1,at(rs)%npol
      netQ = netQ + atq(at(rs)%pol(i))*atq(at(rs)%pol(i))
    end do
  end do
!
! standard RF (GRF limiting case with zero ionic strength)
  if (rf_mode.eq.2) then
    par_POLAR(1) = (epsr-1.0)/((2.0*epsr+1.0)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
! GRF
  else
    kap = sqrt(4.0*PI*electric*invtemp*istr/ens%insV)
    t1 = 1.0 + kap*mcel_cutoff
    t2 = kap*kap*mcel_cutoff2
    par_POLAR(1) = (t1*(epsr-1.0) + 0.5*epsr*t2)/&
 &        (((2.0*epsr+1.0)*t1 + epsr*t2)*(mcel_cutoff**3.0))
    par_POLAR(2) = 1.0/mcel_cutoff + mcel_cutoff2*par_POLAR(1)
  end if
  rfcnst = -0.5*electric*scale_POLAR*netQ*par_POLAR(2)
!
end
!
!-----------------------------------------------------------------------
