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
!--------------------------------------------------------------------
!
subroutine allocate_rest()
!
  use mpistuff
  use energies
!
  implicit none
!
  integer mode
!
  mode = 1
  call allocate_accept(mode)
  call allocate_polyavg(mode)
  call allocate_contacts(mode)
  call allocate_pcs(mode)
  call allocate_dipolavg(mode)
  call allocate_diffrac(mode)
  if (use_EMICRO.EQV..false.) call allocate_ems(mode)
#ifdef ENABLE_MPI
  if ((use_REMC.EQV..true.).OR.((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.))) then
    call allocate_mpire(mode)
  end if
#endif
!
end
!
!--------------------------------------------------------------------
!
subroutine deallocate_all()
!
  use iounit
  use sequen
  use energies
  use inter
  use tabpot
  use polypep
  use paircorr
  use cutoffs
  use params
  use system
  use mpistuff
  use clusters
  use pdb
  use keys
  use commline
  use fos
  use wl
!
  implicit none
!
  integer mode,rs,i
!
  mode = 2
!
  if (pdb_analyze.EQV..true.) then
    if (align%yes.EQV..true.) then
      if (allocated(align%set).EQV..true.) deallocate(align%set)
      if (allocated(align%refxyz).EQV..true.) deallocate(align%refxyz)
      if (allocated(align%curxyz).EQV..true.) deallocate(align%curxyz)
      if (allocated(align%hlpxyz).EQV..true.) deallocate(align%hlpxyz)
      if (allocated(align%diffset).EQV..true.) deallocate(align%diffset)
    end if
    if (allocated(framelst).EQV..true.) deallocate(framelst)
    if (allocated(framewts).EQV..true.) deallocate(framewts)
    if (allocated(framelst2).EQV..true.) deallocate(framelst2)
    if (allocated(framewts2).EQV..true.) deallocate(framewts2)
  end if
  if (allocated(pdbmap).EQV..true.) deallocate(pdbmap)
  if (allocated(pdbmap2).EQV..true.) deallocate(pdbmap2)
  if (use_trajidx.EQV..true.) then
    if (allocated(pdboutlog).EQV..true.) deallocate(pdboutlog)
    if (allocated(pdboutlst).EQV..true.) deallocate(pdboutlst)
  end if
  if (savreq%nats.gt.0) then
    deallocate(savreq%hists)
    deallocate(savreq%idx)
  end if
  if (allocated(hmjam%isin).EQV..true.) deallocate(hmjam%isin)
  if (allocated(hmjam%alpha).EQV..true.) deallocate(hmjam%alpha)
  if (allocated(hmjam%thresh).EQV..true.) deallocate(hmjam%thresh)
  if (allocated(hmjam%boosts).EQV..true.) deallocate(hmjam%boosts)
  if (allocated(hmjam%ca_f).EQV..true.) deallocate(hmjam%ca_f)
  if (allocated(hmjam%ca_f_tr).EQV..true.) deallocate(hmjam%ca_f_tr)
  if (allocated(hmjam%ca_f_bu).EQV..true.) deallocate(hmjam%ca_f_bu)

#ifdef ENABLE_MPI
  if ((use_REMC.EQV..true.).OR.((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.))) then
    call allocate_mpire(mode)
  end if
#endif
!  write(*,*) 'done mpire'
  call allocate_ems(mode)
!  write(*,*) 'done ems'
  call allocate_clusters(mode)
!  write(*,*) 'done clusters'
  call allocate_diffrac(mode)
!  write(*,*) 'done diffrac'
  call allocate_dipolavg(mode)
!  write(*,*) 'done dipolavg'
  call allocate_contacts(mode)
!  write(*,*) 'done contacts'
  call allocate_pcs(mode)
!  write(*,*) 'done pcs'
  call allocate_polyavg(mode)
!  write(*,*) 'done polyavg'
#ifdef ENABLE_THREADS
  call allocate_threads(mode)
!  write(*,*) 'done threads'
#endif
  call allocate_ewalds(mode)
!  write(*,*) 'done ewalds'
  call allocate_accept(mode)
!  write(*,*) 'done accept'
  if (gpc%nos.gt.0) then
    call allocate_gpcpre(mode)
  end if
!  write(*,*) 'done gpcpre'
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    call allocate_moveset(mode)
  end if
!  write(*,*) 'done moveset'
  call allocate_mcgrid(mode)
!  write(*,*) 'done mcgrid'
  call allocate_dssps(mode)
!  write(*,*) 'done dssps'
  call allocate_torsn(mode)
!  write(*,*) 'done torsn'
  call allocate_particlefluc(mode)
!  write(*,*) 'done particlefluc'
  if (use_FEG.EQV..true.) then
    if ((allocated(par_FEG).EQV..false.).OR.&
 &      (allocated(par_FEG3).EQV..false.)) then
    else
      deallocate(par_FEG)
      deallocate(par_FEG3)
    end if
  end if
  call allocate_stericgrid(mode)
!  write(*,*) 'done stericgrids'
  call allocate_distrest(mode)
!  write(*,*) 'done distrest'
  if ((allocated(par_TOR2).EQV..false.).OR.&
 &    (allocated(par_TOR).EQV..false.)) then
  else
    deallocate(par_TOR)
    deallocate(par_TOR2)
  end if
  call allocate_ionize(mode)
!  write(*,*) 'done ionize'
  if ((allocated(par_POLY2).EQV..false.).OR.&
 &    (allocated(par_POLY).EQV..false.)) then
  else
    deallocate(par_POLY)
    deallocate(par_POLY2)
  end if
  if (use_TABUL.EQV..true.) then
    call allocate_tabpot(mode)
    call allocate_tbppre(mode)
  end if
!  write(*,*) 'done tabpot'
  if (use_POLAR.EQV..true.) then
    if (allocated(cglst%it).EQV..true.) then
      deallocate(cglst%it)
      deallocate(cglst%nc)
      deallocate(cglst%tc)
      deallocate(cglst%rsl)
      deallocate(cglst%irsl)
    end if
  end if
!  write(*,*) 'done cglst'
  call allocate_molecule_type(mode)
!  write(*,*) 'done molecule type'
  call allocate_atomprms(mode)
!  write(*,*) 'done atomprms'
  call allocate_forces(mode)
!  write(*,*) 'done forces'
  call deallocate_iaa()
!  write(*,*) 'done iaas/fudges'
  do rs=1,nseq
    do i=1,at(rs)%ndpgrps
      if (allocated(at(rs)%dpgrp(i)%ats).EQV..true.) then
        deallocate(at(rs)%dpgrp(i)%ats)
      end if
    end do
    if (at(rs)%ndpgrps.gt.0) then
      deallocate(at(rs)%dpgrp)
    end if
  end do
!  write(*,*) 'done dpgrps'
  do rs=1,nseq
    do i=1,at(rs)%nfosgrps
      if (allocated(at(rs)%fosgrp(i)%ats).EQV..true.) then
        deallocate(at(rs)%fosgrp(i)%ats)
      end if
      if (allocated(at(rs)%fosgrp(i)%wts).EQV..true.) then
        deallocate(at(rs)%fosgrp(i)%wts)
      end if
    end do
    if (at(rs)%nfosgrps.gt.0) then
      deallocate(at(rs)%fosgrp)
    end if
  end do
!  write(*,*) 'done fosgrps'
  call allocate_inter(mode)
!  write(*,*) 'done rest inter'
  call allocate_atomestls(mode)
!  write(*,*) 'done atom-xyzs'
  call allocate_polypep(mode)
!  write(*,*) 'done polypep'
  call allocate_fyoc(mode)
!  write(*,*) 'done fyoc'
  call allocate_molecule(mode)
!  write(*,*) 'done molecule'
  call allocate_sequen(mode)
!  write(*,*) 'done sequen'
  if (allocated(pdb_unknowns).EQV..true.) deallocate(pdb_unknowns)
  if (allocated(pdb_unkbnd).EQV..true.) deallocate(pdb_unkbnd)
  if (allocated(cm_lst).EQV..true.) deallocate(cm_lst)
  if (allocated(impt_lst).EQV..true.) deallocate(impt_lst)
  if (allocated(di_lst).EQV..true.) deallocate(di_lst)
  if (allocated(ba_lst).EQV..true.) deallocate(ba_lst)
  if (allocated(bo_lst).EQV..true.) deallocate(bo_lst)
  call allocate_params(mode)
!  write(*,*) 'done params'
  if (allocated(key).EQV..true.) then
    do i=1,nkey
      if (allocated(key(i)%line).EQV..true.) deallocate(key(i)%line)
    end do
    deallocate(key)
  end if
  if (allocated(rndstat%itable).EQV..true.) deallocate(rndstat%itable)
  if (allocated(rndstat%rbuf).EQV..true.) deallocate(rndstat%rbuf)
!  write(*,*) 'done keys'
  call allocate_mpistuff(mode)
!  write(*,*) 'done mpistuff'
  if (allocated(args).EQV..true.) deallocate(args)
!  write(*,*) 'done args'
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_molecule(mode)
!
  use iounit
  use sequen
  use atoms
  use molecule
  use movesets
  use cutoffs
!
  implicit none
!
  integer i,mode
!
  if (mode.eq.1) then
    allocate(atmol(nmol,2))
    allocate(rsmol(nmol,2))
    allocate(moltermid(nmol,2))
    allocate(moltypid(nmol))
    allocate(an_grp_mol(nmol))
    allocate(is_solvent(nmol))
    allocate(com(nmol,3))
    allocate(rgpcs(nmol,3,3))
    allocate(rgevs(nmol,3))
    allocate(rgv(nmol))
    allocate(comm(nmol,3))
    allocate(rgpcsm(nmol,3,3))
    allocate(rgevsm(nmol,3))
    allocate(rgvm(nmol))
    allocate(rgevsref(nmol,3))
    allocate(rgpcsref(nmol,3,3))
    allocate(rgvref(nmol))
    allocate(comref(nmol,3))
    allocate(commref(nmol,3))
    allocate(molfrzidx(nmol))
    allocate(molinfo(nmol,3))
    allocate(movingmass(nmol+1,3))
!
! molecular termini flags and default freeze index
!
    do i=1,nmol
      moltermid(i,1) = 0
      moltermid(i,2) = 0
      molfrzidx(i) = 1
      is_solvent(i) = .false.
    end do
!
  else if (mode.eq.2) then
    if ((allocated(atmol).EQV..false.).OR.&
 &      (allocated(rsmol).EQV..false.).OR.&
 &      (allocated(moltermid).EQV..false.).OR.&
 &      (allocated(moltypid).EQV..false.).OR.&
 &      (allocated(an_grp_mol).EQV..false.).OR.&
 &      (allocated(is_solvent).EQV..false.).OR.&
 &      (allocated(com).EQV..false.).OR.&
 &      (allocated(rgpcs).EQV..false.).OR.&
 &      (allocated(rgevs).EQV..false.).OR.&
 &      (allocated(rgv).EQV..false.).OR.&
 &      (allocated(comm).EQV..false.).OR.&
 &      (allocated(rgpcsm).EQV..false.).OR.&
 &      (allocated(rgevsm).EQV..false.).OR.&
 &      (allocated(rgvm).EQV..false.).OR.&
 &      (allocated(rgpcsref).EQV..false.).OR.&
 &      (allocated(rgevsref).EQV..false.).OR.&
 &      (allocated(rgvref).EQV..false.).OR.&
 &      (allocated(comref).EQV..false.).OR.&
 &      (allocated(molinfo).EQV..false.).OR.&
 &      (allocated(commref).EQV..false.).OR.&
 &      (allocated(molfrzidx).EQV..false.).OR.&
 &      (allocated(movingmass).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all molecular arr&
 &ays were not properly allocated, or prematurely freed. Please fix &
 &or report this problem.'
      call fexit()
    else
      deallocate(atmol)
      deallocate(rsmol)
      deallocate(moltermid)
      deallocate(moltypid)
      deallocate(an_grp_mol)
      deallocate(is_solvent)
      deallocate(com)
      deallocate(rgpcs)
      deallocate(rgevs)
      deallocate(rgv)
      deallocate(comm)
      deallocate(rgpcsm)
      deallocate(rgevsm)
      deallocate(rgvm)
      deallocate(rgevsref)
      deallocate(rgpcsref)
      deallocate(rgvref)
      deallocate(comref)
      deallocate(commref)
      deallocate(molfrzidx)
      deallocate(molinfo)
      deallocate(movingmass)
      if (allocated(molangr)) deallocate(molangr)
      if (allocated(solutes)) deallocate(solutes)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_molecule(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_molecule_type(mode)
!
  use iounit
  use molecule
  use shakeetal
!
  implicit none
!
  integer i,mode
!
  if (mode.eq.1) then
    allocate(moltyp(nmoltyp,2))
    allocate(molvol(nmoltyp,2))
    allocate(molmass(nmoltyp))
    allocate(molcontlen(nmoltyp))
    allocate(molischg(nmoltyp))
    allocate(do_pol(nmoltyp))
    allocate(do_tors(nmoltyp))
    allocate(do_pers(nmoltyp))
    allocate(ntormol(nmoltyp))
    allocate(othidxmol(nmoltyp))
    allocate(constraints(nmoltyp))
    do i=1,nmoltyp
      molmass(i) = 0.0
      molcontlen(i) = 0.0
      molvol(i,1) = 0.0
      molvol(i,2) = 0.0
      molischg(i) = .false.
      do_tors(i) = .false.
      do_pol(i) = .false.
      do_pers(i) = .false.
      ntormol(i) = 0
      othidxmol(i) = 0
      constraints(i)%nr = 0
    end do
  else if (mode.eq.2) then
    if (allocated(settlelst1).EQV..true.) deallocate(settlelst1)
    if (allocated(settlelst2).EQV..true.) deallocate(settlelst2)
    if (allocated(settlelst3).EQV..true.) deallocate(settlelst3)
    if (allocated(settlelst4).EQV..true.) deallocate(settlelst4)
    if (allocated(settlelst5).EQV..true.) deallocate(settlelst5)
    if (allocated(settlelst6).EQV..true.) deallocate(settlelst6)
    if (allocated(nonsettlelst).EQV..true.) deallocate(nonsettlelst)
    do i=1,cart_cons_grps
      if (allocated(constraints(i)%idx).EQV..true.) deallocate(constraints(i)%idx)
      if (allocated(constraints(i)%idx3).EQV..true.) deallocate(constraints(i)%idx3)
      if (allocated(constraints(i)%idx4).EQV..true.) deallocate(constraints(i)%idx4)
      if (allocated(constraints(i)%uidx).EQV..true.) deallocate(constraints(i)%uidx)
      if (allocated(constraints(i)%cidx).EQV..true.) deallocate(constraints(i)%cidx)
      if (allocated(constraints(i)%ncis).EQV..true.) deallocate(constraints(i)%ncis)
      if (allocated(constraints(i)%ccoeff).EQV..true.) deallocate(constraints(i)%ccoeff)
      if (allocated(constraints(i)%mapidx).EQV..true.) deallocate(constraints(i)%mapidx)
      if (allocated(constraints(i)%eqd2).EQV..true.) deallocate(constraints(i)%eqd2)
      if (allocated(constraints(i)%eqa3).EQV..true.) deallocate(constraints(i)%eqa3)
      if (allocated(constraints(i)%eqi).EQV..true.) deallocate(constraints(i)%eqi)
      if (allocated(constraints(i)%amat).EQV..true.) deallocate(constraints(i)%amat)
      if (allocated(constraints(i)%massless).EQV..true.) deallocate(constraints(i)%massless)
      if (allocated(constraints(i)%isinter).EQV..true.) deallocate(constraints(i)%isinter)
    end do
    if ((allocated(molvol).EQV..false.).OR.&
 &      (allocated(moltyp).EQV..false.).OR.&
 &      (allocated(molmass).EQV..false.).OR.&
 &      (allocated(molcontlen).EQV..false.).OR.&
 &      (allocated(molischg).EQV..false.).OR.&
 &      (allocated(do_pol).EQV..false.).OR.&
 &      (allocated(do_tors).EQV..false.).OR.&
 &      (allocated(do_pers).EQV..false.).OR.&
 &      (allocated(ntormol).EQV..false.).OR.&
 &      (allocated(othidxmol).EQV..false.).OR.&
 &      (allocated(constraints).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all molecule type&
 & arrays were not properly allocated, or prematurely freed. Please &
 &fix or report this problem.'
      call fexit()
    else
      deallocate(moltyp)
      deallocate(molvol)
      deallocate(molmass)
      deallocate(molcontlen)
      deallocate(molischg)
      deallocate(do_pol)
      deallocate(do_tors)
      deallocate(do_pers)
      deallocate(ntormol)
      deallocate(othidxmol)
      deallocate(constraints)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_molecule_type(...) with un&
 &known mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_sequen(mode)
!
  use iounit
  use sequen
  use energies
  use system
  use fos
  use contacts
  use cutoffs
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer i,j,mode
!
  if (mode.eq.1) then
    allocate(seqtyp(nseq))
    allocate(seqflag(nseq))
    allocate(natres(nseq))
    allocate(refat(nseq))
    allocate(resrad(nseq))
    allocate(resvol(nseq,3))
    allocate(chgflag(nseq))
    allocate(netchg(nseq))
    allocate(notors(nseq))
    allocate(molofrs(nseq))
    allocate(resfy(nseq))
    allocate(resnuc(nseq))
    allocate(resnucpuc(nseq))
    allocate(resw(nseq))
    allocate(reschi(nseq))
    allocate(respuc(nseq))
    allocate(seqpolty(nseq))
    if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.&
 &  ((nequil.lt.nsim).AND.(contactcalc.le.nsim))) then
      allocate(rsp_mat(nseq,nseq))
    end if
    rsp_nbl(:)%sralsz = 100
    rsp_nbl(:)%tralsz = 100
    rsp_nbl(:)%lralsz = 100
    do i=1,2
      allocate(rsp_nbl(i)%sr(rsp_nbl(i)%sralsz,2))
      allocate(rsp_nbl(i)%tr(rsp_nbl(i)%tralsz,2))
      allocate(rsp_nbl(i)%lr(rsp_nbl(i)%lralsz,2))
    end do
    allocate(rsp_vec(nseq))
    allocate(rs_nbl(nseq))
    if (use_IMPSOLV.EQV..true.) then
      allocate(freesolv(nseq,3))
      allocate(rs_vec(nseq))
      allocate(rs_vec2(nseq))
    end if
    allocate(rsinfo(nseq,4))
    allocate(vhlper(max(25,nseq),3))
!
    do i=1,nseq
      seqpolty(i) = 'X'
      notors(i) = .false.
      chgflag(i) = .false.
      netchg(i) = 0
      resrad(i) = 15.0
      natres(i) = -1000  ! trigger seg-fault if not set properly
      seqtyp(i) = -1000  ! --
      seqflag(i) = -1000
      refat(i) = -1000   ! --
      molofrs(i) = -1000 ! --
      rs_nbl(i)%ntmpanb = 0
      rs_nbl(i)%nnbs = 0
      rs_nbl(i)%nnbats = 0
      rs_nbl(i)%nnbtrs = 0
      rs_nbl(i)%nnbtrats = 0
      rs_nbl(i)%nnblrs = 0
      rs_nbl(i)%nnblrats = 0
      rs_nbl(i)%ngnbs = 0
      rs_nbl(i)%ngnbats = 0
      rs_nbl(i)%ngnbtrs = 0
      rs_nbl(i)%ngnbtrats = 0
      rs_nbl(i)%nwnbs = 0
      rs_nbl(i)%nwnbats = 0
      rs_nbl(i)%nwnbtrs = 0
      rs_nbl(i)%nwnbtrats = 0
      rs_nbl(i)%nwnblrs = 0
      rs_nbl(i)%nwnblrats = 0
      rs_nbl(i)%ntabnbs = 0
      rs_nbl(i)%ntabias = 0
      rs_nbl(i)%nbalsz = 10
      rs_nbl(i)%tralsz = 10
      rs_nbl(i)%lralsz = 10
      rs_nbl(i)%gnbalsz = 10
      rs_nbl(i)%gtralsz = 10
      rs_nbl(i)%glralsz = 10
      rs_nbl(i)%wnbalsz = 10
      rs_nbl(i)%wtralsz = 10
      rs_nbl(i)%wlralsz = 10
      rs_nbl(i)%tmpalsz = 10
      allocate(rs_nbl(i)%nb(rs_nbl(i)%nbalsz))
      allocate(rs_nbl(i)%nbtr(rs_nbl(i)%tralsz))
      allocate(rs_nbl(i)%trsvec(rs_nbl(i)%tralsz,3))
      allocate(rs_nbl(i)%nblr(rs_nbl(i)%lralsz))
      allocate(rs_nbl(i)%gnb(rs_nbl(i)%gnbalsz))
      allocate(rs_nbl(i)%gnbtr(rs_nbl(i)%gtralsz))
      allocate(rs_nbl(i)%gtrsvec(rs_nbl(i)%gtralsz,3))
      allocate(rs_nbl(i)%tmpanb(rs_nbl(i)%tmpalsz))
      if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.&
 &  ((nequil.lt.nsim).AND.(contactcalc.le.nsim))) then
        do j=1,nseq
          rsp_mat(i,j) = 0
        end do
      end if
      rsp_vec(i) = 0
      if (use_IMPSOLV.EQV..true.) then
        rs_vec(i) = 0
        rs_vec2(i) = 0
      end if
    end do
    if (n_crosslinks.gt.0) then
      allocate(crlk_idx(nseq))
      allocate(crosslink(n_crosslinks))
      do i=1,n_crosslinks
        crosslink(i)%nolks = 0
        crosslink(i)%itstype = 0
        allocate(crosslink(i)%olks(n_crosslinks-1))
      end do
    end if
!
  else if (mode.eq.2) then
    if ((allocated(seqtyp).EQV..false.).OR.&
 &      (allocated(seqflag).EQV..false.).OR.&
 &      (allocated(natres).EQV..false.).OR.&
 &      (allocated(refat).EQV..false.).OR.&
 &      (allocated(seqpolty).EQV..false.).OR.&
 &      (allocated(resrad).EQV..false.).OR.&
 &      (allocated(resvol).EQV..false.).OR.&
 &      (allocated(chgflag).EQV..false.).OR.&
 &      (allocated(netchg).EQV..false.).OR.&
 &      (allocated(notors).EQV..false.).OR.&
 &      (allocated(molofrs).EQV..false.).OR.&
 &      (allocated(rs_nbl).EQV..false.).OR.&
 &      (allocated(rsp_vec).EQV..false.).OR.&
 &      (allocated(rsinfo).EQV..false.).OR.&
 &      (allocated(vhlper).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all sequence&
 & arrays were not properly allocated, or prematurely freed. Please &
 &fix or report this problem.'
      call fexit()
    else
      deallocate(seqtyp)
      deallocate(seqflag)
      deallocate(natres)
      deallocate(seqpolty)
      deallocate(refat)
      deallocate(resrad)
      deallocate(molofrs)
      deallocate(resvol)
      deallocate(chgflag)
      deallocate(netchg)
      deallocate(notors)
      deallocate(rsinfo)
      deallocate(vhlper)
      deallocate(rsp_vec)
      if (allocated(rsp_mat)) deallocate(rsp_mat)
      do i=1,2
        if (allocated(rsp_nbl(i)%sr).EQV..true.) deallocate(rsp_nbl(i)%sr)
        if (allocated(rsp_nbl(i)%tr).EQV..true.) deallocate(rsp_nbl(i)%tr)
        if (allocated(rsp_nbl(i)%lr).EQV..true.) deallocate(rsp_nbl(i)%lr)
      end do
      if (allocated(rs_nbl).EQV..true.) then
        do i=1,nseq
          if (allocated(rs_nbl(i)%tmpanb).EQV..true.) then
            deallocate(rs_nbl(i)%tmpanb)
          end if
          if (allocated(rs_nbl(i)%nb).EQV..true.) then
            deallocate(rs_nbl(i)%nb)
          end if
          if (allocated(rs_nbl(i)%gnb).EQV..true.) then
            deallocate(rs_nbl(i)%gnb)
          end if
          if (allocated(rs_nbl(i)%wnb).EQV..true.) then
            deallocate(rs_nbl(i)%wnb)
          end if
          if (allocated(rs_nbl(i)%nbtr).EQV..true.) then
            deallocate(rs_nbl(i)%nbtr)
            deallocate(rs_nbl(i)%trsvec)
          end if
          if (allocated(rs_nbl(i)%gnbtr).EQV..true.) then
            deallocate(rs_nbl(i)%gnbtr)
            deallocate(rs_nbl(i)%gtrsvec)
          end if
          if (allocated(rs_nbl(i)%wnbtr).EQV..true.) then
            deallocate(rs_nbl(i)%wnbtr)
            deallocate(rs_nbl(i)%wtrsvec)
          end if
          if (allocated(rs_nbl(i)%nblr).EQV..true.) then
            deallocate(rs_nbl(i)%nblr)
          end if
          if (allocated(rs_nbl(i)%wnblr).EQV..true.) then
            deallocate(rs_nbl(i)%wnblr)
          end if
          if (allocated(rs_nbl(i)%tabnb).EQV..true.) then
            deallocate(rs_nbl(i)%tabnb)
          end if
        end do
      end if
      deallocate(rs_nbl)
      if (allocated(topgrplst).EQV..true.) deallocate(topgrplst)
      if (allocated(topgrpxyz).EQV..true.) deallocate(topgrpxyz)
      if (allocated(which_topgrp).EQV..true.) deallocate(which_topgrp)
    end if
    if (use_IMPSOLV.EQV..true.) then
      if ((allocated(rs_vec).EQV..false.).OR.&
 &        (allocated(rs_vec2).EQV..false.).OR.&
 &        (allocated(freesolv).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Rs_vec{2} was not prope&
 &rly allocated, or prematurely freed. Please fix or report this pro&
 &blem.'
        call fexit()
      else
        deallocate(rs_vec)
        deallocate(rs_vec2)
        deallocate(freesolv)
      end if
    end if
    if (n_crosslinks.gt.0) then
      do i=1,n_crosslinks
        if (allocated(crosslink(i)%olks).EQV..true.) deallocate(crosslink(i)%olks)
        if (allocated(crosslink(i)%exclin).EQV..true.) deallocate(crosslink(i)%exclin)
        if (allocated(crosslink(i)%is14in).EQV..true.) deallocate(crosslink(i)%is14in)
        if (allocated(crosslink(i)%ljpars).EQV..true.) deallocate(crosslink(i)%ljpars)
        if (allocated(crosslink(i)%cbpars).EQV..true.) deallocate(crosslink(i)%cbpars)
        if (allocated(crosslink(i)%exclpol).EQV..true.) deallocate(crosslink(i)%exclpol)
        if (allocated(crosslink(i)%is14pol).EQV..true.) deallocate(crosslink(i)%is14pol)
      end do
      if ((allocated(crosslink).EQV..false.).OR.&
 &        (allocated(crlk_idx).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Crosslink data structures were not prope&
 &rly allocated, or prematurely freed. Please fix or report this problem.'
        call fexit()
      else
        deallocate(crosslink)
        deallocate(crlk_idx)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_sequen(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_polypep(mode)
!
  use iounit
  use sequen
  use polypep
  use energies
!
  implicit none
!
  integer i,mode
!
  if (mode.eq.1) then
    allocate(ni(nseq))
    allocate(cai(nseq))
    allocate(ci(nseq))
    allocate(oi(nseq))
    allocate(hni(nseq))
    allocate(nuci(nseq,6))
    allocate(at(nseq))
    do i=1,nseq
!     this could be done more carefully: we net overallocate by at least a factor of 2
      at(i)%nbb = 0
      at(i)%nsc = 0
      at(i)%npol = 0
      at(i)%ndpgrps = 0
      at(i)%nfosgrps = 0
      allocate(at(i)%bb(natres(i)))
      allocate(at(i)%sc(natres(i)))
      if (use_POLAR.EQV..true.) then
        allocate(at(i)%pol(natres(i)))
      end if
    end do
  else if (mode.eq.2) then
    if ((allocated(ni).EQV..false.).OR.&
 &      (allocated(cai).EQV..false.).OR.&
 &      (allocated(ci).EQV..false.).OR.&
 &      (allocated(oi).EQV..false.).OR.&
 &      (allocated(hni).EQV..false.).OR.&
 &      (allocated(nuci).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some atomic pointer array&
 &s were not properly allocated, or prematurely freed. Please fix or&
 & report this problem.'
      call fexit()
    else
      deallocate(ni)
      deallocate(cai)
      deallocate(ci)
      deallocate(oi)
      deallocate(hni)
      deallocate(nuci)
    end if
    if (allocated(at).EQV..false.) then
      write(ilog,*) 'Fatal. This is a bug. BB/SC atomic pointer arra&
 &y was not properly allocated, or prematurely freed. Please fix or &
 &report this problem.'
      call fexit()
    else
      do i=1,nseq
        if ((allocated(at(i)%bb).EQV..false.).OR.&
 &          (allocated(at(i)%sc).EQV..false.)) then
          write(ilog,*) 'Fatal. This is a bug. BB/SC atomic pointer &
 & arrays were not properly allocated, or prematurely freed for resi&
 &due ',i,'. Please fix or report this problem.'
          call fexit()
        else
          deallocate(at(i)%bb)
          deallocate(at(i)%sc)
        end if
        if (use_POLAR.EQV..true.) then
          if (allocated(at(i)%pol).EQV..false.) then
            write(ilog,*) 'Fatal. This is a bug. Atomic charge point&
 &er arrays were not properly allocated, or prematurely freed for re&
 &sidue ',i,'. Please fix or report this problem.'
          call fexit()
          else
            deallocate(at(i)%pol)
          end if
        end if
      end do
      deallocate(at)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_polypep(...) with unknown &
 &mode (offending mode is ',mode,').'
    call fexit()
  end if
!
  end     
!
!---------------------------------------------------------------------------------
!
subroutine allocate_fyoc(mode)
!
  use iounit
  use sequen
  use fyoc
  use torsn
  use movesets
!
  implicit none
!
  integer i,j,mode
!
  if (mode.eq.1) then
!   dihedral pointer arrays, as long as MAXCHI is small the static allocation is OK
    allocate(phi(nseq))
    allocate(phish(nseq))
    allocate(fline(nseq))
    allocate(fline2(nseq))
    allocate(pucline(nseq))
    allocate(fnr(nseq))
    allocate(psi(nseq))
    allocate(psish(nseq))
    allocate(yline(nseq))
    allocate(yline2(nseq))
    allocate(ynr(nseq))
    allocate(nucs(6,nseq))
    allocate(nnucs(nseq))
    allocate(nucsline(6,nseq))
    allocate(nucsnr(6,nseq))
    allocate(omega(nseq))
    allocate(wline(nseq))
    allocate(wnr(nseq))
    allocate(nchi(nseq))
    allocate(chi(MAXCHI,nseq))
    allocate(chiline(MAXCHI,nseq))
    allocate(chinr(MAXCHI,nseq))
    allocate(chiral(nseq))
    allocate(disulf(nseq))
!   LCT stuff
    if (use_lctmoves.EQV..true.) then
      allocate(lct(ntorlcs))
      allocate(lctold(ntorlcs))
    end if
!   initialize (these must be by default <0, as we often check "if (?line.gt.0)")
!   also, the spacings are set to what they need to be for standard geometry
    do i=1,nseq
      fline(i) = -1000
      fline2(i) = -1000
      yline(i) = -1000
      yline2(i) = -1000
      wline(i) = -1000
      wnr(i) = -1000
      fnr(i) = -1000
      ynr(i) = -1000
      do j=1,MAXCHI
        chiline(j,i) = -1000
        chinr(j,i) = -1000
      end do
      do j=1,6
        nucsline(j,i) = -1000
        nucsnr(j,i) = -1000
      end do
      pucline(i) = -1000
      disulf(i) = 0
      phish(i) = 180.0d0
      psish(i) = 0.0d0
      psi(i) = 180.0
      phi(i) = 180.0
      omega(i) = 180.0
      nucs(:,i) = 180.0
      chi(:,i) = 180.0
    end do
!
  else if (mode.eq.2) then
    if ((allocated(phi).EQV..false.).OR.&
 &      (allocated(phish).EQV..false.).OR.&
 &      (allocated(fline).EQV..false.).OR.&
 &      (allocated(pucline).EQV..false.).OR.&
 &      (allocated(fline2).EQV..false.).OR.&
 &      (allocated(fnr).EQV..false.).OR.&
 &      (allocated(psi).EQV..false.).OR.&
 &      (allocated(psish).EQV..false.).OR.&
 &      (allocated(yline).EQV..false.).OR.&
 &      (allocated(yline2).EQV..false.).OR.&
 &      (allocated(ynr).EQV..false.).OR.&
 &      (allocated(nucs).EQV..false.).OR.&
 &      (allocated(nnucs).EQV..false.).OR.&
 &      (allocated(nucsline).EQV..false.).OR.&
 &      (allocated(nucsnr).EQV..false.).OR.&
 &      (allocated(nchi).EQV..false.).OR.&
 &      (allocated(chi).EQV..false.).OR.&
 &      (allocated(chiline).EQV..false.).OR.&
 &      (allocated(chinr).EQV..false.).OR.&
 &      (allocated(omega).EQV..false.).OR.&
 &      (allocated(wline).EQV..false.).OR.&
 &      (allocated(wnr).EQV..false.).OR.&
 &      (allocated(chiral).EQV..false.).OR.&
 &      (allocated(disulf).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all dihedral poin&
 &ter arrays were not properly allocated, or prematurely freed. Plea&
 &se fix or report this problem.'
      call fexit()
    else
      deallocate(phi)
      deallocate(phish)
      deallocate(fline)
      deallocate(fline2)
      deallocate(pucline)
      deallocate(fnr)
      deallocate(psi)
      deallocate(psish)
      deallocate(yline)
      deallocate(yline2)
      deallocate(ynr)
      deallocate(omega)
      deallocate(wline)
      deallocate(wnr)
      deallocate(nchi)
      deallocate(chi)
      deallocate(chiline)
      deallocate(chinr)
      deallocate(chiral)
      deallocate(disulf)
      deallocate(nucs)
      deallocate(nnucs)
      deallocate(nucsline)
      deallocate(nucsnr)
    end if
    if (use_lctmoves.EQV..true.) then
      if ((allocated(lct).EQV..false.).OR.&
 &      (allocated(lctold).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Some or all dihedral po&
 &inter arrays were not properly allocated, or prematurely freed. Pl&
 &ease fix or report this problem.'
        call fexit()
      else
        deallocate(lct)
        deallocate(lctold)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_fyoc(...) with unknown mod&
 &e (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_inter(mode)
!
  use iounit
  use sequen
  use inter
  use energies
!
  implicit none
!
  integer mode
!
  if (mode.eq.1) then
!   the rest of this is done in unbond/polar_groups 
    allocate(nrsintra(nseq))
    allocate(nrsnb(nseq))
    allocate(nrsbl(nseq))
    allocate(nrsba(nseq))
    allocate(nrsbleff(nseq))
    allocate(nrsbaeff(nseq))
    allocate(nrsdi(nseq))
    allocate(nrsdieff(nseq))
    allocate(nrsimpt(nseq))
    allocate(nrscm(nseq))
    allocate(nrscmeff(nseq))
    allocate(nrsimpteff(nseq))
!   charge arrays
    if (use_POLAR.EQV..true.) then
      allocate(nrpolintra(nseq))
      allocate(nrpolnb(nseq))
      allocate(nrexpolin(nseq))
      allocate(nrexpolnb(nseq))
    end if
  else if (mode.eq.2) then
    if ((allocated(nrsintra).EQV..false.).OR.&
 &      (allocated(nrsnb).EQV..false.).OR.&
 &      (allocated(nrsbl).EQV..false.).OR.&
 &      (allocated(nrsba).EQV..false.).OR.&
 &      (allocated(nrsdi).EQV..false.).OR.&
 &      (allocated(nrsbaeff).EQV..false.).OR.&
 &      (allocated(nrsdieff).EQV..false.).OR.&
 &      (allocated(nrsdieff).EQV..false.).OR.&
 &      (allocated(nrsimpt).EQV..false.).OR.&
 &      (allocated(nrsimpteff).EQV..false.).OR.&
 &      (allocated(nrscm).EQV..false.).OR.&
 &      (allocated(nrscmeff).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all dihedral poin&
 &ter arrays were not properly allocated, or prematurely freed. Plea&
 &se fix or report this problem.'
      call fexit()
    else
      deallocate(nrsintra)
      deallocate(nrsnb)
      deallocate(nrsbl)
      deallocate(nrsba)
      deallocate(nrsbleff)
      deallocate(nrsbaeff)
      deallocate(nrsdi)
      deallocate(nrsdieff)
      deallocate(nrsimpt)
      deallocate(nrsimpteff)
      deallocate(nrscm)
      deallocate(nrscmeff)
    end if
    if (use_POLAR.EQV..true.) then
      if ((allocated(nrpolintra).EQV..false.).OR.&
 &      (allocated(nrpolnb).EQV..false.).OR.&
 &      (allocated(nrexpolin).EQV..false.).OR.&
 &      (allocated(nrexpolnb).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Some or all polar inter&
 &action arrays were not properly allocated, or prematurely freed. P&
 &lease fix or report this problem.'
        call fexit()
      else
        deallocate(nrpolintra)
        deallocate(nrpolnb)
        deallocate(nrexpolin)
        deallocate(nrexpolnb)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_inter(...) with unknown mo&
 &de (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine deallocate_iaa()
!
  use inter
  use iounit
  use sequen
!
  implicit none
!
  integer rs
!
  if (allocated(iaa).EQV..false.) then
    write(ilog,*) 'Fatal. This is a bug. The interaction array struc&
 &ture was not properly allocated, or prematurely freed. Please fix &
 &or report this problem.'
    call fexit()
  else
    do rs=1,nseq
      if (allocated(iaa(rs)%atin)) then
        deallocate(iaa(rs)%atin)
      end if
      if (allocated(iaa(rs)%atnb)) then
        deallocate(iaa(rs)%atnb)
      end if
      if (allocated(iaa(rs)%polin)) then
        deallocate(iaa(rs)%polin)
      end if
      if (allocated(iaa(rs)%polnb)) then
        deallocate(iaa(rs)%polnb)
      end if
      if (allocated(iaa(rs)%expolin)) then
        deallocate(iaa(rs)%expolin)
      end if
      if (allocated(iaa(rs)%expolnb)) then
        deallocate(iaa(rs)%expolnb)
      end if
      if (allocated(iaa(rs)%bl)) then
        deallocate(iaa(rs)%bl)
      end if
      if (allocated(iaa(rs)%typ_bl)) then
        deallocate(iaa(rs)%typ_bl)
      end if
      if (allocated(iaa(rs)%par_bl)) then
        deallocate(iaa(rs)%par_bl)
      end if
      if (allocated(iaa(rs)%ba)) then
        deallocate(iaa(rs)%ba)
      end if
      if (allocated(iaa(rs)%typ_ba)) then
        deallocate(iaa(rs)%typ_ba)
      end if
      if (allocated(iaa(rs)%par_ba)) then
        deallocate(iaa(rs)%par_ba)
      end if
      if (allocated(iaa(rs)%di)) then
        deallocate(iaa(rs)%di)
      end if
      if (allocated(iaa(rs)%typ_di)) then
        deallocate(iaa(rs)%typ_di)
      end if
      if (allocated(iaa(rs)%par_di)) then
        deallocate(iaa(rs)%par_di)
      end if
      if (allocated(iaa(rs)%impt)) then
        deallocate(iaa(rs)%impt)
      end if
      if (allocated(iaa(rs)%typ_impt)) then
        deallocate(iaa(rs)%typ_impt)
      end if
      if (allocated(iaa(rs)%par_impt)) then
        deallocate(iaa(rs)%par_impt)
      end if
      if (allocated(iaa(rs)%cm)) then
        deallocate(iaa(rs)%cm)
      end if
      if (allocated(iaa(rs)%typ_cm)) then
        deallocate(iaa(rs)%typ_cm)
      end if
    end do
    deallocate(iaa)
  end if
  if (allocated(fudge).EQV..false.) then
    write(ilog,*) 'Fatal. This is a bug. The fudge factor structure &
 &was not properly allocated, or prematurely freed. Please fix or re&
 &port this problem.'
    call fexit()
  else
    do rs=1,nseq
      if (allocated(fudge(rs)%rsin)) then
        deallocate(fudge(rs)%rsin)
      end if
      if (allocated(fudge(rs)%rsnb)) then
        deallocate(fudge(rs)%rsnb)
      end if
      if (allocated(fudge(rs)%elin)) then
        deallocate(fudge(rs)%elin)
      end if
      if (allocated(fudge(rs)%elnb)) then
        deallocate(fudge(rs)%elnb)
      end if
      if (allocated(fudge(rs)%rsin_ljs)) then
        deallocate(fudge(rs)%rsin_ljs)
      end if
      if (allocated(fudge(rs)%rsnb_ljs)) then
        deallocate(fudge(rs)%rsnb_ljs)
      end if
      if (allocated(fudge(rs)%rsin_lje)) then
        deallocate(fudge(rs)%rsin_lje)
      end if
      if (allocated(fudge(rs)%rsnb_lje)) then
        deallocate(fudge(rs)%rsnb_lje)
      end if
    end do
    deallocate(fudge)
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_atomestls(mode)
!
  use iounit
  use sequen
  use atoms
  use zmatrix
!
  implicit none
!
  integer i,mode,natts
!
! we'll overallocate (the +20 is irrelevant, but remember that the natres-values have
! a little buffer-space, too)
  natts = 0
  do i=1,nseq
    natts = natts + natres(i)
  end do
  natts = natts + 20
  if (mode.eq.1) then
    allocate(x(natts))
    allocate(y(natts))
    allocate(z(natts))
    allocate(xref(natts))
    allocate(yref(natts))
    allocate(zref(natts))
    allocate(atnam(natts))
    allocate(attyp(natts))
    allocate(b_type(natts))
    allocate(mass(natts))
    allocate(ivms(natts))
    allocate(atmres(natts))
    allocate(eqatm(natts))
    allocate(n12(natts))
    allocate(i12(MAXVLNC,natts))
    allocate(n13(natts))
    allocate(i13(MAXVLNC*(MAXVLNC-1),natts))
    allocate(n14(natts))
    allocate(i14(MAXVLNC*(MAXVLNC-1)*(MAXVLNC-1),natts))
    allocate(iz(4,natts))
    allocate(blen(natts))
    allocate(bang(natts))
    allocate(ztor(natts))
    allocate(izrot(natts))
    allocate(ztorpr(natts))
    allocate(bangpr(natts))
    allocate(blenpr(natts))
    allocate(iadd(2,natts))
!   hopefully if through a bug an atom slips through our attention, this will seg-fault it
    atmres(:) = -1000
    eqatm(:) = -1000
    attyp(:) = 0
    b_type(:) = 0
    x(:) = 0.0
    y(:) = 0.0
    z(:) = 0.0
    iz(:,:) = 0
   else if (mode.eq.2) then
    if ((allocated(x).EQV..false.).OR.&
 &      (allocated(y).EQV..false.).OR.&
 &      (allocated(z).EQV..false.).OR.&
 &      (allocated(xref).EQV..false.).OR.&
 &      (allocated(yref).EQV..false.).OR.&
 &      (allocated(zref).EQV..false.).OR.&
 &      (allocated(atnam).EQV..false.).OR.&
 &      (allocated(attyp).EQV..false.).OR.&
 &      (allocated(b_type).EQV..false.).OR.&
 &      (allocated(mass).EQV..false.).OR.&
 &      (allocated(ivms).EQV..false.).OR.&
 &      (allocated(n12).EQV..false.).OR.&
 &      (allocated(i12).EQV..false.).OR.&
 &      (allocated(n13).EQV..false.).OR.&
 &      (allocated(i13).EQV..false.).OR.&
 &      (allocated(n14).EQV..false.).OR.&
 &      (allocated(i14).EQV..false.).OR.&
 &      (allocated(iz).EQV..false.).OR.&
 &      (allocated(izrot).EQV..false.).OR.&
 &      (allocated(blen).EQV..false.).OR.&
 &      (allocated(bang).EQV..false.).OR.&
 &      (allocated(ztor).EQV..false.).OR.&
 &      (allocated(ztorpr).EQV..false.).OR.&
 &      (allocated(bangpr).EQV..false.).OR.&
 &      (allocated(blenpr).EQV..false.).OR.&
 &      (allocated(iadd).EQV..false.).OR.&
 &      (allocated(atmres).EQV..false.).OR.&
 &      (allocated(eqatm).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all atom arrays w&
 &ere not properly allocated, or prematurely freed. Please fix or re&
 &port this problem.'
      call fexit()
    else
      do i=1,n
        if (allocated(izrot(i)%rotis).EQV..true.) deallocate(izrot(i)%rotis)
        if (allocated(izrot(i)%diffis).EQV..true.) deallocate(izrot(i)%diffis)
        if (allocated(izrot(i)%treevs).EQV..true.) deallocate(izrot(i)%treevs)
      end do
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(xref)
      deallocate(yref)
      deallocate(zref)
      deallocate(attyp)
      deallocate(b_type)
      deallocate(ivms)
      deallocate(mass)
      deallocate(atnam)
      deallocate(atmres)
      deallocate(eqatm)
      deallocate(n12)
      deallocate(i12)
      deallocate(n13)
      deallocate(i13)
      deallocate(n14)
      deallocate(i14)
      deallocate(iz)
      deallocate(izrot)
      deallocate(blen)
      deallocate(bang)
      deallocate(ztor)
      deallocate(ztorpr)
      deallocate(bangpr)
      deallocate(blenpr)
      deallocate(iadd)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_atomestls(...) with unknow&
 &n mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!
!---------------------------------------------------------------------------------
!
subroutine allocate_atomprms(mode)
!
  use iounit
  use atoms
  use energies
  use molecule
  use system
  use cutoffs
!
  implicit none
!
  integer mode
!
  if (mode.eq.1) then
    allocate(atr(n+20))
    allocate(atvol(n+20))
    allocate(atsav(n+20))
    allocate(atstv(n+20))
    allocate(atbvol(n+20))
    allocate(atsavred(n+20))
    allocate(atsavmaxfr(n+20))
    allocate(atinfo(4,n))
    if (par_IMPSOLV2(1).eq.2) then
      allocate(atsavinfo(3,n))
    else
      allocate(atsavinfo(2,n))
    end if
    allocate(t_ca_f(3,n))
    allocate(t_ca_f_tr(3,n))
!
    if ((use_IMPSOLV.EQV..true.).OR.(savcalc.le.nsim)) then
      allocate(atsavprm(n,8))
      if (savcalc.le.nsim) then
        allocate(atsavavg(2,n))
        allocate(natsavavg(nmol))
        atsavavg(:,:) = 0.0
        natsavavg(:) = 0
      end if
    end if
!
    allocate(svte(n))
    svte(:) = 0.0
    if (use_IMPSOLV.EQV..true.) then
      allocate(svbu(n))
      allocate(which_fosg(n))
      svte(:) = 0.0
      svbu(:) = 0.0
      which_fosg(:) = 0
    end if
!
    if (use_POLAR.EQV..true.) then
      allocate(atq(n))
      allocate(which_dpg(n))
      atq(:) = 0.0
      which_dpg(:) = 0
      if (use_IMPSOLV.EQV..true.) then
        allocate(scrq(n))
        scrq(:) = 0.0
      end if
    end if
!
  else if (mode.eq.2) then
    if ((allocated(atr).EQV..false.).OR.&
 &      (allocated(atvol).EQV..false.).OR.&
 &      (allocated(svte).EQV..false.).OR.&
 &      (allocated(atsav).EQV..false.).OR.&
 &      (allocated(atstv).EQV..false.).OR.&
 &      (allocated(atbvol).EQV..false.).OR.&
 &      (allocated(atinfo).EQV..false.).OR.&
 &      (allocated(atsavinfo).EQV..false.).OR.&
 &      (allocated(t_ca_f).EQV..false.).OR.&
 &      (allocated(t_ca_f_tr).EQV..false.).OR.&
 &      (allocated(atsavred).EQV..false.).OR.&
 &      (allocated(atsavmaxfr).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all atom paramete&
 &rs were not properly allocated, or prematurely freed. Please fix o&
 &r report this problem.'
      call fexit()
    else
      deallocate(atr)
      deallocate(atvol)
      deallocate(svte)
      deallocate(atsav)
      deallocate(atstv)
      deallocate(atbvol)
      deallocate(atinfo)
      deallocate(atsavinfo)
      deallocate(t_ca_f)
      deallocate(t_ca_f_tr)
      deallocate(atsavred)
      deallocate(atsavmaxfr)
    end if
    if ((use_IMPSOLV.EQV..true.).OR.(savcalc.le.nsim)) then
      if (allocated(atsavprm).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. Some atom parameter (at&
 &savprm) was not properly allocated, or prematurely freed. Please f&
 &ix or report this problem.'
      else 
        deallocate(atsavprm)
      end if
      if (savcalc.le.nsim) then
        if ((allocated(atsavavg).EQV..false.).OR.&
 &         (allocated(natsavavg).EQV..false.)) then
          write(ilog,*) 'Fatal. This is a bug. Some atom parameter (&
 &(n)atsavavg) was not properly allocated, or prematurely freed. Please&
 & fix or report this problem.'
        else
          deallocate(atsavavg)
          deallocate(natsavavg)
        end if
      end if
    end if
    if (use_IMPSOLV.EQV..true.) then
      if ((allocated(svbu).EQV..false.).OR.&
 &        (allocated(which_fosg).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Some or all atom parame&
 &ters were not properly allocated, or prematurely freed. Please fix&
 & or report this problem (IMPSOLV-dependent).'
      else 
        deallocate(svbu)
        deallocate(which_fosg)
      end if
    end if
    if (use_POLAR.EQV..true.) then
      if ((allocated(atq).EQV..false.).OR.&
 &        (allocated(which_dpg).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Atomic charge array was&
 & not properly allocated, or prematurely freed. Please fix or repor&
 &t this problem.'
      else 
        deallocate(atq)
        deallocate(which_dpg)
      end if
      if (use_IMPSOLV.EQV..true.) then
        if (allocated(scrq).EQV..false.) then
          write(ilog,*) 'Fatal. This is a bug. Screened charge array&
 & was not properly allocated, or prematurely freed. Please fix or r&
 &eport this problem.'
        else 
          deallocate(scrq)
        end if
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_atomprms(...) with unknown&
 & mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_threads(mode)
!
  use iounit
  use atoms
  use forces
  use threads
  use system
  use energies
!
  implicit none
!
  integer mode,i,tpi
!
  if (mode.eq.1) then
#ifdef ENABLE_THREADS
    allocate(thr_limits(110,thrdat%maxn))
    allocate(thr_timings(200,thrdat%maxn))
    allocate(thr_rutil(max(tstat%n_tgrps+1,max(MAXENERGYTERMS*3,100)),thrdat%maxn))
    thr_timings(:,:) = 0
    allocate(thr_dlb(20,3))
    if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      allocate(thr_rsp_nbl(2*thrdat%maxn))
      thr_rsp_nbl(:)%sralsz = 100
      thr_rsp_nbl(:)%tralsz = 100
      thr_rsp_nbl(:)%lralsz = 100
      do i=1,2*thrdat%maxn
        allocate(thr_rsp_nbl(i)%sr(thr_rsp_nbl(i)%sralsz,3))
        allocate(thr_rsp_nbl(i)%tr(thr_rsp_nbl(i)%tralsz,3))
        allocate(thr_rsp_nbl(i)%lr(thr_rsp_nbl(i)%lralsz,3))
      end do
      thr_rsp_nbl(:)%srnrs = 0
      thr_rsp_nbl(:)%trnrs = 0
      thr_rsp_nbl(:)%lrnrs = 0
    end if
    allocate(thr_svte(n,thrdat%maxn))
    thr_svte(:,:) = 0.0
    if (use_dyn.EQV..true.) then
      allocate(thr_ca_f(3,n,thrdat%maxn))
      allocate(thr_ca_f_tr(3,n,thrdat%maxn))
      if (use_IMPSOLV.EQV..true.) then
        allocate(thr_sisa(n,thrdat%maxn))
        thr_sisa(:,:)%nix = 0
        thr_sisa(:,:)%alsz = 0
        if (use_POLAR.EQV..true.) then
          allocate(thr_sum_scrcbs(n,thrdat%maxn))
          allocate(thr_sum_scrcbs_tr(n,thrdat%maxn))
        end if
      end if
    end if
#else
    allocate(thr_timings(200,1))
#endif
!
  else if (mode.eq.2) then
    if ((allocated(thr_dlb).EQV..false.).OR.&
 &      (allocated(thr_limits).EQV..false.).OR.&
 &      (allocated(thr_timings).EQV..false.).OR.&
 &      (allocated(thr_rutil).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some or all basic arrays for multi-threaded use &
 &of CAMPARI were not properly allocated or prematurely freed. Please fix or report this problem.'
      call fexit()
    else
      deallocate(thr_dlb)
      deallocate(thr_timings)
      deallocate(thr_limits)
      deallocate(thr_rutil)
      if (allocated(thr_rsp_nbl).EQV..true.) then
        do i=1,thrdat%maxn
          if (allocated(thr_rsp_nbl(i)%sr).EQV..true.) deallocate(thr_rsp_nbl(i)%sr)
          if (allocated(thr_rsp_nbl(i)%tr).EQV..true.) deallocate(thr_rsp_nbl(i)%tr)
          if (allocated(thr_rsp_nbl(i)%lr).EQV..true.) deallocate(thr_rsp_nbl(i)%lr)
          if (allocated(thr_rsp_nbl(i)%tmpl).EQV..true.) deallocate(thr_rsp_nbl(i)%tmpl)
        end do
        deallocate(thr_rsp_nbl)
      end if
      if (allocated(ccg_limits).EQV..true.) deallocate(ccg_limits)
      if (allocated(mlg_limits).EQV..true.) deallocate(mlg_limits)
      if (allocated(mlg_ixs).EQV..true.) deallocate(mlg_ixs)
      if (allocated(thr_mlgix).EQV..true.) deallocate(thr_mlgix)
      if (allocated(thr_ca_f_tr).EQV..true.) deallocate(thr_ca_f_tr)
      if (allocated(thr_ca_f).EQV..true.) deallocate(thr_ca_f)
      if (allocated(thr_sisa).EQV..true.) then
        do tpi=1,thrdat%maxn
          do i=1,n
            if (allocated(thr_sisa(i,tpi)%ix).EQV..true.) deallocate(thr_sisa(i,tpi)%ix)
            if (allocated(thr_sisa(i,tpi)%dr).EQV..true.) deallocate(thr_sisa(i,tpi)%dr)
            if (allocated(thr_sisa(i,tpi)%qr).EQV..true.) deallocate(thr_sisa(i,tpi)%qr)
          end do
        end do
        deallocate(thr_sisa)
      end if
      if (allocated(thr_svte).EQV..true.) deallocate(thr_svte)
      if (allocated(thr_lhlper).EQV..true.) deallocate(thr_lhlper)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_threads(...) with unknown&
 & mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!------------------------------------------------------------------------------------
!
subroutine allocate_moveset(mode)
!
  use iounit
  use sequen
  use accept
  use molecule
  use movesets
  use energies, ONLY: MAXENERGYTERMS
!
  implicit none
!
  integer mode
!
  if (mode.eq.1) then
    if ((rigidfreq.gt.0.0).AND.(clurb_freq.gt.0.0)) then
      allocate(cur_clurbt(nmol,3))
      allocate(clurbf(nmol))
      allocate(clurbi(nmol))
      if (clurb_rdfreq.lt.1.0) then
        allocate(clurb_clus(nmol,nmol))
        allocate(clurb_cluszs(nmol))
      end if
    end if
    allocate(mvcur%evecd(MAXENERGYTERMS))
    allocate(mvcur%evecp(MAXENERGYTERMS))
    allocate(mvcur%eveca(MAXENERGYTERMS))
    allocate(mvcur%covec(longestmol*10))
  else if (mode.eq.2) then
    if (allocated(mvcur%evecd).EQV..true.) deallocate(mvcur%evecd)
    if (allocated(mvcur%eveca).EQV..true.) deallocate(mvcur%eveca)
    if (allocated(mvcur%evecp).EQV..true.) deallocate(mvcur%evecp)
    if (allocated(mvcur%covec).EQV..true.) deallocate(mvcur%covec)
    if (allocated(chilst%idx).EQV..true.) deallocate(chilst%idx)
    if (allocated(chilst%idx2).EQV..true.) deallocate(chilst%idx2)
    if (allocated(chilst%wt).EQV..true.) deallocate(chilst%wt)
    if (allocated(fylst%idx).EQV..true.) deallocate(fylst%idx)
    if (allocated(fylst%idx2).EQV..true.) deallocate(fylst%idx2)
    if (allocated(fylst%wt).EQV..true.) deallocate(fylst%wt)
    if (allocated(nuclst%idx).EQV..true.) deallocate(nuclst%idx)
    if (allocated(nuclst%idx2).EQV..true.) deallocate(nuclst%idx2)
    if (allocated(nuclst%wt).EQV..true.) deallocate(nuclst%wt)
    if (allocated(wlst%idx).EQV..true.) deallocate(wlst%idx)
    if (allocated(wlst%idx2).EQV..true.) deallocate(wlst%idx2)
    if (allocated(wlst%wt).EQV..true.) deallocate(wlst%wt)
    if (allocated(puclst%idx).EQV..true.) deallocate(puclst%idx)
    if (allocated(puclst%idx2).EQV..true.) deallocate(puclst%idx2)
    if (allocated(puclst%wt).EQV..true.) deallocate(puclst%wt)
    if (allocated(nucpuclst%idx).EQV..true.) deallocate(nucpuclst%idx)
    if (allocated(nucpuclst%idx2).EQV..true.) deallocate(nucpuclst%idx2)
    if (allocated(nucpuclst%wt).EQV..true.) deallocate(nucpuclst%wt)
    if (allocated(djcrlst%idx).EQV..true.) deallocate(djcrlst%idx)
    if (allocated(djcrlst%idx2).EQV..true.) deallocate(djcrlst%idx2)
    if (allocated(djcrlst%wt).EQV..true.) deallocate(djcrlst%wt)
    if (allocated(docrlst%idx).EQV..true.) deallocate(docrlst%idx)
    if (allocated(docrlst%idx2).EQV..true.) deallocate(docrlst%idx2)
    if (allocated(docrlst%wt).EQV..true.) deallocate(docrlst%wt)
    if (allocated(nuccrlst%idx).EQV..true.) deallocate(nuccrlst%idx)
    if (allocated(nuccrlst%idx2).EQV..true.) deallocate(nuccrlst%idx2)
    if (allocated(nuccrlst%wt).EQV..true.) deallocate(nuccrlst%wt)
    if (allocated(ujcrlst%idx).EQV..true.) deallocate(ujcrlst%idx)
    if (allocated(ujcrlst%idx2).EQV..true.) deallocate(ujcrlst%idx2)
    if (allocated(ujcrlst%wt).EQV..true.) deallocate(ujcrlst%wt)
    if (allocated(sjcrlst%idx).EQV..true.) deallocate(sjcrlst%idx)
    if (allocated(sjcrlst%idx2).EQV..true.) deallocate(sjcrlst%idx2)
    if (allocated(sjcrlst%wt).EQV..true.) deallocate(sjcrlst%wt)
    if (allocated(rblst%idx).EQV..true.) deallocate(rblst%idx)
    if (allocated(rblst%idx2).EQV..true.) deallocate(rblst%idx2)
    if (allocated(rblst%wt).EQV..true.) deallocate(rblst%wt)
    if (allocated(clurblst%idx).EQV..true.) deallocate(clurblst%idx)
    if (allocated(clurblst%idx2).EQV..true.) deallocate(clurblst%idx2)
    if (allocated(clurblst%wt).EQV..true.) deallocate(clurblst%wt)
    if (allocated(fluclst%idx).EQV..true.) deallocate(fluclst%idx)
    if (allocated(fluclst%idx2).EQV..true.) deallocate(fluclst%idx2)
    if (allocated(fluclst%wt).EQV..true.) deallocate(fluclst%wt)
    if (allocated(unklst%idx).EQV..true.) deallocate(unklst%idx)
    if (allocated(unklst%idx2).EQV..true.) deallocate(unklst%idx2)
    if (allocated(unklst%wt).EQV..true.) deallocate(unklst%wt)
    if (allocated(unslst%idx).EQV..true.) deallocate(unslst%idx)
    if (allocated(unslst%idx2).EQV..true.) deallocate(unslst%idx2)
    if (allocated(unslst%wt).EQV..true.) deallocate(unslst%wt)
    if (allocated(natlst%idx).EQV..true.) deallocate(natlst%idx)
    if (allocated(natlst%idx2).EQV..true.) deallocate(natlst%idx2)
    if (allocated(natlst%wt).EQV..true.) deallocate(natlst%wt)
    if ((rigidfreq.gt.0.0).AND.(clurb_freq.gt.0.0)) then
      deallocate(cur_clurbt)
      deallocate(clurbf)
      deallocate(clurbi)
      if (clurb_rdfreq.lt.1.0) then
        deallocate(clurb_clus)
        deallocate(clurb_cluszs)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_moveset(...) with unknown &
 &mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine allocate_accept(mode)
!
  use iounit
  use sequen
  use accept
  use molecule
  use movesets
!
  implicit none
!
  integer i,mode
!
  if (mode.eq.1) then
    if (have_rigid.EQV..true.) then
      allocate(acc%rigid(nmol))
      do i=1,nmol
        acc%rigid(i) = 0
      end do
    end if
    if (have_omega.EQV..true.) then
      allocate(acc%omega(nseq))
      do i=1,nseq
        acc%omega(i) = 0
      end do
    end if
    if ((have_nuccr.EQV..true.).OR.(have_cr.EQV..true.)) then
      allocate(acc%cr(nseq))
      do i=1,nseq
        acc%cr(i) = 0
      end do
    end if
    if ((have_chi.EQV..true.).AND.(phfreq.lt.1.0)) then
      allocate(acc%chi(nseq))
      do i=1,nseq
        acc%chi(i) = 0
      end do
    end if
    if ((have_pucker.EQV..true.).OR.(have_nucpuck.EQV..true.)) then
      allocate(acc%pucker(nseq))
      do i=1,nseq
        acc%pucker(i) = 0
      end do
    end if
    if (have_pivot.EQV..true.) then
      allocate(acc%fy(nseq))
      do i=1,nseq
        acc%fy(i) = 0
      end do
      if ((.NOT.(allocated(acc%chi))).AND.(use_globmoves.EQV..true.)) then
        allocate(acc%chi(nseq))
        acc%chi(1:nseq) = 0
      end if
    end if
    if ((have_nuc.EQV..true.).AND.(nucpuckfreq.lt.1.0).AND.(nuccrfreq.lt.1.0)) then
      allocate(acc%nuc(nseq))
      do i=1,nseq
        acc%nuc(i) = 0
      end do
    end if
    if (have_other.EQV..true.) then
      allocate(acc%other(nseq))
      do i=1,nseq
        acc%other(i) = 0
      end do
    end if
    if (have_particlefluc.EQV..true.) then
      allocate(acc%insert(nmol))
      allocate(acc%delete(nmol))
      allocate(acc%permute(nmol))
      do i=1,nmol
        acc%insert(i) = 0
        acc%delete(i) = 0
        acc%permute(i) = 0
      end do
    end if
  else if (mode.eq.2) then
    if (allocated(acc%rigid)) then
      deallocate(acc%rigid)
    end if
    if (allocated(acc%omega)) then
      deallocate(acc%omega)
    end if
    if (allocated(acc%fy)) then
      deallocate(acc%fy)
    end if
    if (allocated(acc%cr)) then
      deallocate(acc%cr)
    end if
    if (allocated(acc%chi)) then
      deallocate(acc%chi)
    end if
    if (allocated(acc%pucker)) then
      deallocate(acc%pucker)
    end if
    if (allocated(acc%nuc)) then
      deallocate(acc%nuc)
    end if
    if (allocated(acc%delete)) then
      deallocate(acc%delete)
    end if
    if (allocated(acc%insert)) then
      deallocate(acc%insert)
    end if
    if (allocated(acc%permute)) then
      deallocate(acc%permute)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_accept(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_mcgrid(mode)
!
  use iounit
  use mcgrid
  use cutoffs
  use sequen
  use energies
  use system
!
  implicit none
!
  integer i,mode
!
  if (use_mcgrid.EQV..false.) return
!
  if (mode.eq.1) then
    if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      allocate(grid%tmpl(nseq,1))
      if ((use_cutoffs.EQV..true.).AND.(use_POLAR.EQV..true.).AND.(cglst%ncrs.gt.0)) then
        if (lrel_mc.eq.1) then
          allocate(rsp_lmat(cglst%ncrs,nseq))
        else if ((lrel_mc.eq.2).OR.(lrel_mc.eq.3)) then
          allocate(rsp_lmat(cglst%ncrs,cglst%ncrs))
        end if
      end if
    end if
    allocate(grid%nbnr(grid%dim(1)*grid%dim(2)*grid%dim(3)+1))
    allocate(grid%maxr(grid%dim(1)*grid%dim(2)*grid%dim(3)+1))
    allocate(grid%gpnbnr(grid%dim(1)*grid%dim(2)*grid%dim(3)+1))
    allocate(grid%b_gpnbnr(grid%dim(1)*grid%dim(2)*grid%dim(3)+1))
    do i=1,grid%dim(1)*grid%dim(2)*grid%dim(3)+1
      grid%nbnr(i) = 0
      grid%gpnbnr(i) = 0
      grid%b_gpnbnr(i) = 0
    end do
    allocate(grid%resgp(nseq))
  else if (mode.eq.2) then
    if ((allocated(grid%nbnr).EQV..false.).OR.&
 &      (allocated(grid%maxr).EQV..false.).OR.&
 &      (allocated(grid%gpnbnr).EQV..false.).OR.&
 &      (allocated(grid%b_gpnbnr).EQV..false.).OR.&
 &      (allocated(grid%resgp).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. MC-grid parameters were n&
 &ot properly allocated, or prematurely freed. Please fix or report &
 &this problem.'
      call fexit()
    else
      deallocate(grid%nbnr)
      deallocate(grid%maxr)
      deallocate(grid%gpnbnr)
      deallocate(grid%b_gpnbnr)
      deallocate(grid%resgp)
      if (allocated(grid%nblist).EQV..true.) deallocate(grid%nblist)
      if (allocated(grid%gpnblist).EQV..true.) deallocate(grid%gpnblist)
      if (allocated(grid%b_gpnblist).EQV..true.) deallocate(grid%b_gpnblist)
      if (allocated(grid%tgrpgp).EQV..true.) deallocate(grid%tgrpgp)
      if (allocated(grid%tmpl).EQV..true.) deallocate(grid%tmpl)
      if (allocated(rsp_lmat).EQV..true.) deallocate(rsp_lmat)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_mcgrid(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_ems(mode)
!
  use iounit
  use ems
  use system
  use energies
  use atoms
!
  implicit none
!
  integer mode,k,t1,t2,ints(3)
  logical atrue
!
  atrue = .true.
!
  if ((emcalc.gt.nsim).AND.(use_EMICRO.EQV..false.)) return
!
  if (mode.eq.1) then
 46 format('Warning. In dimension ',i2,', the evaluation grid is not at least of the &
 &resolution of the provided density map (',a,') (',g10.4,' vs. ',g10.4,' A). Adjusting evaluation grid.')
 47 format('Warning. In dimension ',i2,', the resolution of the evaluation grid is not an integer multiple &
 &of that of the provided density map (',a,') (',g10.4,' vs. ',g10.4,' A). Adjusting evaluation grid.')
 477 format('Warning. In dimension ',i2,', the density analysis grid does not line up with the &
 &system dimensions. Adjusting resolution slightly.')
 478 format('Warning. In dimension ',i2,', the density analysis grid does not line up with the &
 &system dimensions. Adjusting system size slightly.')
    if (use_EMICRO.EQV..true.) then
      call strlims(emmapfile,t1,t2)
      if ((emgrid%deltas(1).eq.0.0).AND.(emgrid%deltas(2).eq.0.0).AND.(emgrid%deltas(3).eq.0.0)) then
        emgrid%deltas(:) = emingrid%deltas(:)
      else
        do k=1,3
          if (emingrid%deltas(k).lt.emgrid%deltas(k)) then
            emgrid%deltas(k) = emingrid%deltas(k)
            write(ilog,46) k,emmapfile(t1:t2),emgrid%deltas(k),emingrid%deltas(k)
          else if (emingrid%deltas(k).gt.emgrid%deltas(k)) then
            emgrid%deltas(k) = emingrid%deltas(k)/(1.0*nint(emingrid%deltas(k)/emgrid%deltas(k)))
            if (mod(emingrid%deltas(k),emgrid%deltas(k)).gt.1.0e-7) then
              write(ilog,47) k,emmapfile(t1:t2),emgrid%deltas(k),emingrid%deltas(k)
            end if
          end if
        end do
      end if
    else
      if ((emgrid%deltas(1).eq.0.0).AND.(emgrid%deltas(2).eq.0.0).AND.(emgrid%deltas(3).eq.0.0)) then
        write(ilog,*) 'Warning. Specification of grid resolution is required for spatial density &
 &analysis. Disabling.'
        emcalc = nsim + 1
      else if ((emgrid%deltas(1).le.0.0).OR.(emgrid%deltas(2).le.0.0).OR.(emgrid%deltas(3).le.0.0)) then
        write(ilog,*) 'Warning. Grid resolution (FMCSC_EMDELTAS) must be positive, real numbers for all three &
 &dimensions. Disabling spatial density analysis.'
        emcalc = nsim + 1
      end if
      if ((empotprop.eq.3).AND.(use_POLAR.EQV..false.)) then
        write(ilog,*) 'Warning. Use of the polar potential is required for spatial charge density analysis. Disabling.'
        emcalc = nsim + 1
      end if
      if (emcalc.gt.nsim) return
    end if
    if ((bnd_shape.eq.1).AND.(bnd_type.eq.1)) then
      emgrid%dim(:) = ceiling((bnd_params(1:3)-0.1*emgrid%deltas(:))/emgrid%deltas(:))
      emgrid%origin(:) = bnd_params(4:6)
      if (use_EMICRO.EQV..false.) then
        do k=1,3
          if (abs(bnd_params(k)-emgrid%dim(k)*emgrid%deltas(k)).gt.1.0e-5) then
            write(ilog,477) k
          end if
          emgrid%deltas(k) = bnd_params(k)/(1.0*emgrid%dim(k))
        end do
      else
        do k=1,3
          if (abs(bnd_params(k)-emgrid%dim(k)*emgrid%deltas(k)).gt.1.0e-5) then
            write(ilog,478) k
          end if
          bnd_params(k) = emgrid%dim(k)*emgrid%deltas(k)
        end do
        call update_bound(atrue)
      end if
    else if ((bnd_shape.eq.1).AND.((bnd_type.ge.2).AND.(bnd_type.le.4))) then
      emgrid%dim(:) = ceiling((embuf*bnd_params(1:3)-0.1*emgrid%deltas(:))/emgrid%deltas(:))
      emgrid%origin(:) = bnd_params(4:6) - 0.5*(embuf-1.0)*bnd_params(1:3)
      emgrid%origin(:) = emgrid%origin(:) + mod(abs(emingrid%origin(:)-emgrid%origin(:)),emgrid%deltas(:))
    else if (bnd_shape.eq.2) then ! necessarily non-PBC
      emgrid%dim(:) = ceiling(2.0*embuf*bnd_params(4)/emgrid%deltas(:))
      emgrid%origin(:) = bnd_params(1:3) - embuf*bnd_params(4)
      emgrid%origin(:) = emgrid%origin(:) + mod(abs(emingrid%origin(:)-emgrid%origin(:)),emgrid%deltas(:))
    else if ((bnd_shape.eq.3).AND.((bnd_type.ge.2).AND.(bnd_type.le.4))) then
      emgrid%dim(1:2) = ceiling(2.0*embuf*bnd_params(4)/emgrid%deltas(1:2))
      emgrid%dim(3) =  ceiling((embuf*bnd_params(6)-0.1*emgrid%deltas(3))/emgrid%deltas(3))
      emgrid%origin(1:2) = bnd_params(1:2) - embuf*bnd_params(4)
      emgrid%origin(1:2) = emgrid%origin(1:2) + mod(abs(emingrid%origin(1:2)-emgrid%origin(1:2)),emgrid%deltas(1:2))
      emgrid%origin(3) = bnd_params(3) - 0.5*embuf*bnd_params(6)
      emgrid%origin(3) = emgrid%origin(3) + mod(abs(emingrid%origin(3)-emgrid%origin(3)),emgrid%deltas(3))
    else if ((bnd_shape.eq.3).AND.(bnd_type.eq.1)) then ! axis-PBC
      emgrid%dim(1:2) = ceiling(2.0*embuf*bnd_params(4)/emgrid%deltas(1:2))
      emgrid%dim(3) =  ceiling((bnd_params(6)-0.1*emgrid%deltas(3))/emgrid%deltas(3))
      emgrid%origin(1:2) = bnd_params(1:2) - embuf*bnd_params(4)
      emgrid%origin(1:2) = emgrid%origin(1:2) + mod(abs(emingrid%origin(1:2)-emgrid%origin(1:2)),emgrid%deltas(1:2))
      emgrid%origin(3) = bnd_params(3) - 0.5*bnd_params(6)
      emgrid%origin(3) = emgrid%origin(3) + mod(abs(emingrid%origin(3)-emgrid%origin(3)),emgrid%deltas(3))
    end if
 48 format('Warning. In dimension ',i2,', the provided density map (',a,') does not line up with the &
 &grid cells of the evaluation grid. Adjusting origin of provided map slightly (from ',g13.6,' to ',g13.6,').')
 49 format('Fatal. In dimension ',i2,', the provided density map (',a,') does not fit inside the dimensions &
 &of the evaluation grid on the left (lower) side (',g20.12,' vs. ',g20.12,').')
 499 format('Fatal. In dimension ',i2,', the provided density map (',a,') does not fit inside the dimensions &
 &of the evaluation grid on the right (upper) side (',g20.12,' vs. ',g20.12,').')
 50 format('Fatal. In dimension ',i2,', the internal EM map does not have high enough resolution &
 &to support the requested order of B-splines (',i5,' vs. ',i5,').')
    if (use_EMICRO.EQV..true.) then
      call strlims(emmapfile,t1,t2)
      do k=1,3
        if ((bnd_type.eq.1).AND.((bnd_shape.eq.1).OR.((bnd_shape.eq.3).AND.(k.eq.3)))) then
          if (mod(abs(emingrid%origin(k)-emgrid%origin(k)),emgrid%deltas(k)).gt.1.0e-7) then
            if (emingrid%origin(k).gt.emgrid%origin(k)) then
              write(ilog,48) k,emmapfile(t1:t2),emingrid%origin(k),&
 &                                            emingrid%origin(k)-mod(emingrid%origin(k)-emgrid%origin(k),emgrid%deltas(k))
              emingrid%origin(k) = emingrid%origin(k) - mod(emingrid%origin(k)-emgrid%origin(k),emgrid%deltas(k))
            else
              write(ilog,48) k,emmapfile(t1:t2),emingrid%origin(k),&
 &                                            emingrid%origin(k)+mod(emgrid%origin(k)-emingrid%origin(k),emgrid%deltas(k))
              emingrid%origin(k) = emingrid%origin(k) + mod(emgrid%origin(k)-emingrid%origin(k),emgrid%deltas(k))
            end if
          end if
        end if
        if (emgrid%origin(k)-emingrid%origin(k).gt.1.0e-7) then
          write(ilog,49) k,emmapfile(t1:t2),emingrid%origin(k),emgrid%origin(k)
          call fexit()
        else if ((emingrid%origin(k)+emingrid%dim(k)*emingrid%deltas(k)) - &
 &               (emgrid%origin(k)+emgrid%dim(k)*emgrid%deltas(k)).gt.1.0e-7) then
          write(ilog,499) k,emmapfile(t1:t2),emingrid%origin(k)+emingrid%dim(k)*emingrid%deltas(k),&
 &                        emgrid%origin(k)+emgrid%dim(k)*emgrid%deltas(k)
          call fexit()
        end if
        emgrid%idx(k,1) = nint((emingrid%origin(k)-emgrid%origin(k))/emgrid%deltas(k)) + 1
        emgrid%idx(k,2) = nint((emingrid%origin(k)+emingrid%dim(k)*emingrid%deltas(k)-emgrid%origin(k))/emgrid%deltas(k))
        emingrid%shf(k) = 1 - emgrid%idx(k,1)
      end do
      ints(:) = nint(emingrid%deltas(:)/(emgrid%deltas(:)))
      if (bnd_shape.eq.1) then
        if (bnd_type.eq.1) then
          emcoarsegrid%dim(:) = ceiling((bnd_params(1:3)-0.5*emingrid%deltas(:))/emingrid%deltas(:))
          do k=1,3
            if (mod(emgrid%dim(k),emcoarsegrid%dim(k)).ne.0) then
              write(ilog,*) 'Fatal. When using mismatched grids for input density and evaluation, finer number &
 &of cells in dimension ',k,' has to be an exact multiple of coarser number (cuboid with 3D PBC).'
              call fexit()
            end if
          end do
          emcoarsegrid%origin(:) = bnd_params(4:6)
        else
          emcoarsegrid%dim(:) = ceiling((embuf*bnd_params(1:3)-0.5*emingrid%deltas(:))/emingrid%deltas(:))
          emcoarsegrid%origin(:) = emgrid%origin(:)
        end if
      else if (bnd_shape.eq.2) then
        emcoarsegrid%dim(:) = ceiling((2.0*embuf*bnd_params(4)-0.5*emingrid%deltas(:))/emingrid%deltas(:))
        emcoarsegrid%origin(:) = emgrid%origin(:)
      else if ((bnd_shape.eq.3).AND.((bnd_type.ge.2).AND.(bnd_type.le.4))) then
        emcoarsegrid%dim(1:2) = ceiling((2.0*embuf*bnd_params(4)-0.5*emingrid%deltas(1:2))/emingrid%deltas(1:2))
        emcoarsegrid%dim(3) =  ceiling((embuf*bnd_params(6)-0.5*emingrid%deltas(3))/emingrid%deltas(3))
        emcoarsegrid%origin(1:3) = emgrid%origin(1:3)
      else if ((bnd_shape.eq.3).AND.(bnd_type.eq.1)) then
        emcoarsegrid%dim(1:2) = ceiling((2.0*embuf*bnd_params(4)-0.5*emingrid%deltas(1:2))/emingrid%deltas(1:2))
        emcoarsegrid%dim(3) =  ceiling((bnd_params(6)-0.5*emingrid%deltas(3))/emingrid%deltas(3))
        do k=3,3
          if (mod(emgrid%dim(k),emcoarsegrid%dim(k)).ne.0) then
            write(ilog,*) 'Fatal. When using mismatched grids for input density and evaluation, finer number &
 &of cells in dimension',k,' has to be an exact multiple of coarser number.'
            call fexit()
          end if
        end do
        emcoarsegrid%origin(1:3) = emgrid%origin(1:3)
      end if
      emcoarsegrid%deltas(:) = emingrid%deltas(:)
      do k=1,3
        if (ints(k).gt.1) then
          if (mod(ceiling((emingrid%origin(k)-emcoarsegrid%origin(k))/emgrid%deltas(k)+0.5),ints(k)).ne.1) then
            write(ilog,*) 'Fatal. For the mismatched grids for input density and evaluation, registry is off for &
 &dimension ',k,'. Please adjust system origin (and maybe size).'
            call fexit()
          end if
        end if
      end do
    else
      emcoarsegrid%deltas(:) = emgrid%deltas(:)
    end if
!
    do k=1,3
      if (emsplor.ge.floor((emgrid%dim(k)*1.0)/2.0)) then
        write(ilog,50) k,emgrid%dim(k),emsplor
        call fexit()
      end if
    end do
!
    allocate(emgrid%mass(emgrid%dim(1),emgrid%dim(2),emgrid%dim(3)))
    allocate(emgrid%mass2(emgrid%dim(1),emgrid%dim(2),emgrid%dim(3)))
    allocate(emgrid%avgmass(emgrid%dim(1),emgrid%dim(2),emgrid%dim(3)))
    allocate(emgrid%dens(emgrid%dim(1),emgrid%dim(2),emgrid%dim(3)))
!
    allocate(embspl(emsplor,3,n))
    allocate(embspld(emsplor,3,n))
    allocate(emgfls(n,3))
!
    if (use_EMICRO.EQV..true.) then
      allocate(emcoarsegrid%dens(emcoarsegrid%dim(1),emcoarsegrid%dim(2),emcoarsegrid%dim(3)))
      allocate(emcoarsegrid%deltadens(emcoarsegrid%dim(1),emcoarsegrid%dim(2),emcoarsegrid%dim(3)))
    end if
!
    emgrid%mass(:,:,:) = 0.0
    emgrid%mass2(:,:,:) = 0.0
    emgrid%avgmass(:,:,:) = 0.0
    emgrid%dens(:,:,:) = 0.0
!
    call em_getintmaps()
!
  else if (mode.eq.2) then
    if ((allocated(emgrid%mass).EQV..false.).OR.&
 &      (allocated(emgrid%mass2).EQV..false.).OR.&
 &      (allocated(emgrid%dens).EQV..false.).OR.&
 &      (allocated(emgrid%avgmass).EQV..false.).OR.&
 &      (allocated(embspl).EQV..false.).OR.&
 &      (allocated(embspld).EQV..false.).OR.&
 &      (allocated(emgfls).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. EM-grid or spline parameters were n&
 &ot properly allocated, or prematurely freed. Please fix or report &
 &this problem.'
      call fexit()
    else
      if (allocated(emingrid%dens).EQV..true.) deallocate(emingrid%dens)
      if (allocated(emingrid%optdens).EQV..true.) deallocate(emingrid%optdens)
      if (allocated(emcoarsegrid%deltadens).EQV..true.) deallocate(emcoarsegrid%deltadens)
      if (allocated(emcoarsegrid%dens).EQV..true.) deallocate(emcoarsegrid%dens)
      deallocate(emgrid%mass)
      deallocate(emgrid%mass2)
      deallocate(emgrid%avgmass)
      deallocate(emgrid%dens)
      if (allocated(emgrid%xmap).EQV..true.) deallocate(emgrid%xmap)
      if (allocated(emgrid%ymap).EQV..true.) deallocate(emgrid%ymap)
      if (allocated(emgrid%zmap).EQV..true.) deallocate(emgrid%zmap)
      if (allocated(emgrid%xmap2).EQV..true.) deallocate(emgrid%xmap2)
      if (allocated(emgrid%ymap2).EQV..true.) deallocate(emgrid%ymap2)
      if (allocated(emgrid%zmap2).EQV..true.) deallocate(emgrid%zmap2)
      if (allocated(emcoarsegrid%lslc).EQV..true.) deallocate(emcoarsegrid%lslc)
      if (allocated(emcoarsegrid%llin).EQV..true.) deallocate(emcoarsegrid%llin)
      if (allocated(emcoarsegrid%lblk).EQV..true.) deallocate(emcoarsegrid%lblk)
      if (allocated(emcoarsegrid%rslc).EQV..true.) deallocate(emcoarsegrid%rslc)
      if (allocated(emcoarsegrid%rlin).EQV..true.) deallocate(emcoarsegrid%rlin)
      if (allocated(emcoarsegrid%rblk).EQV..true.) deallocate(emcoarsegrid%rblk)
      if (allocated(emcoarsegrid%blklms).EQV..true.) deallocate(emcoarsegrid%blklms)
      deallocate(emgfls)
      deallocate(embspl)
      deallocate(embspld)
      if (allocated(emcustom).EQV..true.) deallocate(emcustom)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_ems(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_stericgrid(mode)
!
! this as of yet pseudo-dynamic: instead of MAXGRID it should be the actual ngrid
!
  use iounit
  use grids
  use movesets
  use aminos
!
  implicit none
!
  integer mode
!
  if (use_stericgrids.EQV..false.) return
!
  if (mode.eq.1) then
    allocate(stgr%it(MAXAMINO,MAXGRID,2))
    allocate(stgr%sfyc(MAXAMINO))
    allocate(stgr%ngr(MAXAMINO))
  else if (mode.eq.2) then
    if ((allocated(stgr%sfyc).EQV..false.).OR.&
 &      (allocated(stgr%it).EQV..false.).OR.&
 &      (allocated(stgr%ngr).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Steric grid parameters we&
 &re not properly allocated, or prematurely freed. Please fix or rep&
 &ort this problem.'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_stericgrid(...) with unkno&
 &wn mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_contacts(mode)
!
  use iounit
  use sequen
  use contacts
  use system
  use molecule
!
  implicit none
!
  integer i,j,mode
!
  if (contactcalc.gt.nsim) return
!
  if (mode.eq.1) then
    allocate(contact_maps(nressolute,nressolute))
    allocate(contact_hists(2,20*nressolute))
    allocate(ncontactavg(nsolutes,nsolutes))
    allocate(clumol(nsolutes,nsolutes))
    allocate(clusys(nsolutes))
    allocate(clucoo(nsolutes,100))
!    allocate(clu_nbl(nsolutes))
!
    do i=1,nressolute
      do j=1,nressolute
        contact_maps(i,j) = 0.0
      end do
    end do
    do i=1,20*nressolute
      contact_hists(1,i) = 0.0
      contact_hists(2,i) = 0.0
    end do
    do i=1,nsolutes
      do j=1,nsolutes
        clumol(i,j) = 0.0
        ncontactavg(i,j) = 0
      end do
      clusys(i) = 0.0
      do j=1,100
        clucoo(i,j) = 0.0
      end do
!      clu_nbl(i)%nnbs = 0
!      clu_nbl(i)%alsz = 10
    end do
  else if (mode.eq.2) then
    if ((allocated(contact_maps).EQV..false.).OR.&
 &      (allocated(contact_hists).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Contact maps and/or histo&
 &grams were not properly allocated, or prematurely freed. Please fi&
 &x or report this problem.'
      call fexit()
    else
      deallocate(contact_maps)
      deallocate(contact_hists)
    end if
    if ((allocated(clumol).EQV..false.).OR.&
 &      (allocated(clusys).EQV..false.).OR.&
 &      (allocated(ncontactavg).EQV..false.).OR.&
! &      (allocated(clu_nbl).EQV..false.).OR.&
 &      (allocated(clucoo).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Some cluster analysis var&
 &iables were not properly allocated, or prematurely freed. Please f&
 &ix or report this problem.'
      call fexit()
    else
      deallocate(clumol)
      deallocate(clusys)
      deallocate(ncontactavg)
      deallocate(clucoo)
    end if
!    do i=1,nsolutes
!      if (allocated(clu_nbl(i)%nb)) deallocate(clu_nbl(i)%nb)
!      if (allocated(clu_nbl(i)%link)) deallocate(clu_nbl(i)%link)
!    end do
!    deallocate(clu_nbl)
  else
    write(ilog,*) 'Fatal. Called allocate_contacts(...) with unknown&
 & mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_polyavg(mode)
!
  use iounit
  use sequen
  use system
  use molecule
  use polyavg
!
  implicit none
!
  integer i,j,mode,k
!
  if ((polcalc.gt.nsim).AND.(rhcalc.gt.nsim)) return
!
  if (mode.eq.1) then
!
    if (polcalc.le.nsim) then
      allocate(rgavg(nangrps))
      allocate(rg2avg(nangrps))
      allocate(vtavg(nangrps))
      allocate(reteavg(nangrps))
      allocate(rgpcavg(nangrps,3))
      allocate(rgterm1avg(nangrps))
      allocate(rgterm2avg(nangrps))
      allocate(rgterm12avg(nangrps))
      allocate(dlavg(nangrps))
      allocate(asphavg(nangrps))
      allocate(acylavg(nangrps))
      allocate(dlstavg(nangrps))
      allocate(rdhist(nangrps,TPBINS,DLSTBINS))
      allocate(rghist(nangrps,maxrgbins))
      allocate(retehist(nangrps,maxrgbins))
      allocate(densproavg(nangrps,maxrgbins))
      allocate(persavg(nangrps,longestmol))
      allocate(turns_rs(nseq,9))
      allocate(npolavg(nangrps))
      allocate(npersavg(nangrps))
      allocate(polynew%rg(nangrps))
      allocate(polynew%rg2(nangrps))
      allocate(polynew%rgt1(nangrps))
      allocate(polynew%rgt2(nangrps))
      allocate(polynew%rgt12(nangrps))
      allocate(polynew%rgev(nangrps,3))
      allocate(polynew%vt(nangrps))
      allocate(polynew%ree(nangrps))
      allocate(polynew%cost(longestmol))
      allocate(polynew%cosa(nangrps,longestmol))
      allocate(polynew%map(nangrps))
      allocate(polynew%dep(nangrps,maxrgbins))
!     initialize
      do i=1,nangrps
        npolavg(i) = 0
        npersavg(i) = 0
        rgavg(i) = 0.0
        rg2avg(i) = 0.0
        vtavg(i) = 0.0
        do j=1,3
          rgpcavg(i,j) = 0.0
        end do
        asphavg(i) = 0.0
        acylavg(i) = 0.0
        dlavg(i) = 0.0
        dlstavg(i) = 0.0
        rgterm1avg(i) = 0.0
        rgterm12avg(i) = 0.0
        rgterm2avg(i) = 0.0
        reteavg(i) = 0.0
        do j=1,TPBINS
          do k=1,DLSTBINS
            rdhist(i,j,k) = 0.0
          end do
        end do
        do j=1,longestmol
          persavg(i,j) = 0.0
        end do
        do j=1,maxrgbins
          retehist(i,j) = 0.0
          rghist(i,j) = 0.0
          densproavg(i,j) = 0.0
        end do
      end do
!
      do i=1,nseq
        do j=1,9
          turns_rs(i,j) = 0.0
        end do
      end do
    end if
!
    if (rhcalc.le.nsim) then
      allocate(rhavg(nangrps,longestmol))
      allocate(pscavg(nangrps,nqv))
      allocate(nrhavg(nangrps))
      allocate(rhravg(nangrps))
      allocate(nsctavg(nangrps))
      allocate(polynew%rhr(nangrps))
      allocate(polynew%rh(nangrps,longestmol))
      allocate(polynew%psc(nangrps,nqv))
      allocate(polynew%nm(nangrps,longestmol))
      allocate(polynew%nmr(nangrps))
!     initialize
      do i=1,nangrps
        nrhavg(i) = 0
        rhravg(i) = 0.0
        nsctavg(i) = 0
        do j=1,longestmol
          rhavg(i,j) = 0.0
        end do
        do j=1,nqv
          pscavg(i,j) = 0.0
        end do
      end do
!
    end if
!
  else if (mode.eq.2) then
!
    if (polcalc.le.nsim) then
      if ((allocated(rgavg).EQV..false.).OR.&
 &        (allocated(rg2avg).EQV..false.).OR.&
 &        (allocated(vtavg).EQV..false.).OR.&
 &        (allocated(reteavg).EQV..false.).OR.&
 &        (allocated(persavg).EQV..false.).OR.&
 &        (allocated(rgterm1avg).EQV..false.).OR.&
 &        (allocated(rgterm2avg).EQV..false.).OR.&
 &        (allocated(rgterm12avg).EQV..false.).OR.&
 &        (allocated(dlavg).EQV..false.).OR.&
 &        (allocated(asphavg).EQV..false.).OR.&
 &        (allocated(acylavg).EQV..false.).OR.&
 &        (allocated(dlstavg).EQV..false.).OR.&
 &        (allocated(rgpcavg).EQV..false.).OR.&
 &        (allocated(rdhist).EQV..false.).OR.&
 &        (allocated(retehist).EQV..false.).OR.&
 &        (allocated(rghist).EQV..false.).OR.&
 &        (allocated(turns_rs).EQV..false.).OR.&
 &        (allocated(densproavg).EQV..false.).OR.&
 &        (allocated(npolavg).EQV..false.).OR.&
 &        (allocated(npersavg).EQV..false.).OR.&
 &        (allocated(polynew%rg).EQV..false.).OR.&
 &        (allocated(polynew%rg2).EQV..false.).OR.&
 &        (allocated(polynew%rgt1).EQV..false.).OR.&
 &        (allocated(polynew%rgt2).EQV..false.).OR.&
 &        (allocated(polynew%rgt12).EQV..false.).OR.&
 &        (allocated(polynew%vt).EQV..false.).OR.&
 &        (allocated(polynew%rgev).EQV..false.).OR.&
 &        (allocated(polynew%ree).EQV..false.).OR.&
 &        (allocated(polynew%dep).EQV..false.).OR.&
 &        (allocated(polynew%cosa).EQV..false.).OR.&
 &        (allocated(polynew%cost).EQV..false.).OR.&
 &        (allocated(polynew%map).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Standard polymer analys&
 &is arrays were not properly allocated, or prematurely freed. Pleas&
 &e fix or report this problem.'
        call fexit()
      else
        deallocate(rgavg)
        deallocate(rg2avg)
        deallocate(vtavg)
        deallocate(reteavg)
        deallocate(persavg)
        deallocate(rgterm1avg)
        deallocate(rgterm2avg)
        deallocate(rgterm12avg)
        deallocate(dlavg)
        deallocate(asphavg)
        deallocate(acylavg)
        deallocate(dlstavg)
        deallocate(rgpcavg)
        deallocate(rdhist)
        deallocate(retehist)
        deallocate(rghist)
        deallocate(turns_rs)
        deallocate(densproavg)
        deallocate(npolavg)
        deallocate(npersavg)
        deallocate(polynew%rg)
        deallocate(polynew%rg2)
        deallocate(polynew%rgt1)
        deallocate(polynew%rgt2)
        deallocate(polynew%rgt12)
        deallocate(polynew%vt)
        deallocate(polynew%rgev)
        deallocate(polynew%ree)
        deallocate(polynew%dep)
        deallocate(polynew%cosa)
        deallocate(polynew%cost)
        deallocate(polynew%map)
      end if
    end if
!
    if (rhcalc.le.nsim) then
      if ((allocated(rhavg).EQV..false.).OR.&
 &        (allocated(pscavg).EQV..false.).OR.&
 &        (allocated(nrhavg).EQV..false.).OR.&
 &        (allocated(rhravg).EQV..false.).OR.&
 &        (allocated(nsctavg).EQV..false.).OR.&
 &        (allocated(polynew%rhr).EQV..false.).OR.&
 &        (allocated(polynew%rh).EQV..false.).OR.&
 &        (allocated(polynew%psc).EQV..false.).OR.&
 &        (allocated(polynew%nm).EQV..false.).OR.&
 &        (allocated(polynew%nmr).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Expensive polymer analy&
 &sis arrays were not properly allocated, or prematurely freed. Plea&
 &se fix or report this problem.'
        call fexit()
      else
        deallocate(rhavg)
        deallocate(pscavg)
        deallocate(nrhavg)
        deallocate(rhravg)
        deallocate(nsctavg)
        deallocate(polynew%rhr)
        deallocate(polynew%rh)
        deallocate(polynew%psc)
        deallocate(polynew%nm)
        deallocate(polynew%nmr)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_polyavg(...) with unknown &
 &mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_torsn(mode)
!
  use iounit
  use sequen
  use system
  use molecule
  use torsn
  use polypep
  use inter
  use polyavg
  use energies
  use atoms
  use zmatrix
  use wl
!
  implicit none
!
  integer j,mode,k,imol,rs
!
!  if ((segcalc.gt.nsim).AND.(angcalc.gt.nsim).AND.&
! &    ((torlccalc.gt.nsim).OR.(ntorlcs.le.0)).AND.&
! &    (use_ZSEC.EQV..false.).AND.(intcalc.gt.nsim)) return
!
  if (mode.eq.1) then
!
    j = ntorpuck
    do k=1,n_crosslinks
      if (crosslink(k)%itstype.le.2) j = j + 3
    end do
    allocate(curtvec(j))
!    
    if (angcalc.le.nsim) then
      torbins = floor(360.0/fyres)+1
      allocate(genfymap(torbins,torbins))
      allocate(jc_res(nseq))
      j = 0
      do imol=1,nmol
        if (do_tors(moltypid(imol)).EQV..false.) cycle
        do rs=rsmol(imol,1)+1,rsmol(imol,2)-1
          if (rs.eq.(rsmol(imol,1)+1)) j = j + 1
          if (at(rs)%nsc.le.0) cycle
          if ((hni(rs).le.0).OR.(cai(rs).le.0).OR.(ni(rs).le.0)) cycle
          if (at(rs)%sc(1).gt.0) then
            if ((mass(at(rs)%sc(1)).gt.4.0).OR.(iz(1,at(rs)%sc(1)).ne.cai(rs)).OR.(n12(at(rs)%sc(1)).gt.1)) cycle
          end if
        end do
      end do
      allocate(jccnt(max(j,2)))
      jccnt(:) = 0
      if (nrspecrama.gt.0) then
        allocate(specfymap(nrspecrama,torbins,torbins))
        specfymap(:,:,:) = 0.0
      end if
      if (nrmolrama.gt.0) then
        allocate(molfymap(nrmolrama,torbins,torbins))
        molfymap(:,:,:) = 0.0
      end if
!     initialize, note we store original size in rama_alsz(1), since nrmolrama, and nrspecrama may be reduced later
      rama_alsz(1) = nrspecrama
      rama_alsz(2) = nrmolrama
      genfymap(:,:) = 0.0
      jc_res(:) = 0.0
    end if
!
    if (intcalc.le.nsim) then
      intszs(1:3,2) = 101
      if (do_ints(1).EQV..true.) then
        intszs(1,1) = sum(nrsbl(1:nseq))
        if (intszs(1,1).gt.0) then
          allocate(bln_mids(intszs(1,1)))
          allocate(bln_hist(intszs(1,1),intszs(1,2)))
        else
          do_ints(1) = .false.
        end if
      end if
      if (do_ints(2).EQV..true.) then
        intszs(2,1) = sum(nrsba(1:nseq))
        if (intszs(2,1).gt.0) then
          allocate(ang_mids(intszs(2,1)))
          allocate(ang_hist(intszs(2,1),intszs(2,2)))
        else
          do_ints(2) = .false.
        end if
      end if
      if (do_ints(3).EQV..true.) then
        intszs(3,1) = sum(nrsimpt(1:nseq))
        if (intszs(3,1).gt.0) then
          allocate(imp_mids(intszs(3,1)))
          allocate(imp_hist(intszs(3,1),intszs(3,2)))
        else
          do_ints(3) = .false.
        end if
      end if
      intszs(4,2) = 100
      if (do_ints(4).EQV..true.) then
        intszs(4,1) = sum(nrsdi(1:nseq))
        if (intszs(4,1).gt.0) then
          allocate(dih_hist(intszs(4,1),intszs(4,2)))
        else
          do_ints(4) = .false.
        end if
      end if
    end if
    if ((do_ints(1).EQV..false.).AND.(do_ints(2).EQV..false.).AND.&
 &      (do_ints(3).EQV..false.).AND.(do_ints(4).EQV..false.)) then
      intcalc = nsim + 1
    end if
!
    if ((torlccalc.le.nsim).AND.(ntorlcs.gt.0)) then
      allocate(torlc_coeff(ntorlcs,ntorsn))
      allocate(torlc_coeff2(ntorlcs,ntorsn))
      allocate(torlc_hist(ntorlcs,MAXTORLCBINS))
      allocate(lct_weight(ntorlcs))
      lct_weight(:) = 0.0
      torlc_hist(:,:) = 0.0
      torlc_coeff(:,:) = 0.0
      torlc_coeff2(:,:) = 0.0
    end if
!
    if (covcalc.le.nsim) then
      allocate(itrcv(nangrps))
      allocate(ncovaravg(nangrps))
      allocate(trcvopen(nangrps))
      ncovaravg(:) = 0
      trcvopen(:) = .false.
    end if
!
    if (intcalc.le.nsim) then
      if (do_ints(1).EQV..true.) then
        bln_hist(:,:) = 0.0 
      end if
      if (do_ints(2).EQV..true.) then
        ang_hist(:,:) = 0.0
      end if
      if (do_ints(3).EQV..true.) then
        imp_hist(:,:) = 0.0
      end if
      if (do_ints(4).EQV..true.) then
        dih_hist(:,:) = 0.0
      end if
    end if
!
    if ((segcalc.le.nsim).OR.(use_ZSEC.EQV..true.).OR.(do_wanglandau.EQV..true.)) then
      allocate(z_alpha(nmol))
      allocate(z_beta(nmol))
      if (segcalc.le.nsim) then
        seg_alsz = 0
        do k=1,nangrps
          if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
          if (rsmol(molangr(k,1),2).gt.seg_alsz) seg_alsz = rsmol(molangr(k,1),2)
        end do
        if (seg_alsz.le.0) then
          write(ilog,*) 'No polypeptides found with suitable residues to perform secondary structure&
 & segment or global content analyses. Turning off analyses.'
          segcalc = nsim + 1
        end if
      end if
!
      if (segcalc.le.nsim) then
        allocate(sgprscnt(nangrps))
        allocate(seg_perrs(NSEG,seg_alsz,maxseglen))
        allocate(z_hist(nangrps,2,100))
        allocate(z_hist2(nangrps,100,100))
        allocate(z_2dhist(nangrps,TPBINS,100))
!       initialize
        z_hist(:,:,:) = 0.0
        z_hist2(:,:,:) = 0.0
        z_2dhist(:,:,:) = 0.0
        sgprscnt(:) = 0
        seg_perrs(:,:,:) = 0.0
      end if
    end if
!
  else if (mode.eq.2) then
!
    if (allocated(torzmlst).EQV..true.) deallocate(torzmlst)
!
    if (allocated(curtvec).EQV..false.) then
      write(ilog,*) 'Fatal. This is a bug. Temporary torsion arr&
 &ay was not properly allocated, or prematurely freed. Please fix or&
 & report this problem.'
      call fexit()
    else
      deallocate(curtvec)
    end if
    if (angcalc.le.nsim) then
      if (allocated(jc_res).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. J-coupling constant arr&
 &ay was not properly allocated, or prematurely freed. Please fix or&
 & report this problem.'
        call fexit()
      else
        deallocate(jc_res)
      end if
      if (allocated(jccnt).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. J-coupling constant counter &
 &array was not properly allocated, or prematurely freed. Please fix or&
 & report this problem.'
        call fexit()
      else
        deallocate(jccnt)
      end if
      if (allocated(genfymap).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. Global Ramachandran map&
 & was not properly allocated, or prematurely freed. Please fix or r&
 &eport this problem.'
        call fexit()
      else
        deallocate(genfymap)
      end if
      if (rama_alsz(1).gt.0) then
        if (allocated(specfymap).EQV..false.) then
          write(ilog,*) 'Fatal. This is a bug. Some Ramachandran map&
 & was not properly allocated, or prematurely freed. Please fix or r&
 &eport this problem.'
          call fexit()
        else
          deallocate(specfymap)
        end if
      end if
      if (rama_alsz(2).gt.0) then
        if (allocated(molfymap).EQV..false.) then
          write(ilog,*) 'Fatal. This is a bug. Some Ramachandran map&
 & was not properly allocated, or prematurely freed. Please fix or r&
 &eport this problem.'
          call fexit()
        else
          deallocate(molfymap)
        end if
      end if
    end if
!
    if (intcalc.le.nsim) then
      if (do_ints(1).EQV..true.) then
        if ((allocated(bln_mids).EQV..false.).OR.&
 &          (allocated(bln_hist).EQV..false.)) then
          write(ilog,*) 'Fatal. This is a bug. Some internal coordin&
 &ate histogram was not properly allocated or prematurely freed. Ple&
 &ase fix or report this problem.'
          call fexit()
        else
          deallocate(bln_mids)
          deallocate(bln_hist)
        end if
      end if
      if (do_ints(2).EQV..true.) then
        if ((allocated(ang_mids).EQV..false.).OR.&
 &          (allocated(ang_hist).EQV..false.)) then
          write(ilog,*) 'Fatal. This is a bug. Some internal coordin&
 &ate histogram was not properly allocated or prematurely freed. Ple&
 &ase fix or report this problem.'
          call fexit()
        else
          deallocate(ang_mids)
          deallocate(ang_hist)
        end if
      end if
      if (do_ints(3).EQV..true.) then
        if ((allocated(imp_mids).EQV..false.).OR.&
 &          (allocated(imp_hist).EQV..false.)) then
          write(ilog,*) 'Fatal. This is a bug. Some internal coordin&
 &ate histogram was not properly allocated or prematurely freed. Ple&
 &ase fix or report this problem.'
          call fexit()
        else
          deallocate(imp_mids)
          deallocate(imp_hist)
        end if
      end if
      if (do_ints(4).EQV..true.) then
        if (allocated(dih_hist).EQV..false.) then
          write(ilog,*) 'Fatal. This is a bug. Some internal coordin&
 &ate histogram was not properly allocated or prematurely freed. Ple&
 &ase fix or report this problem.'
          call fexit()
        else
          deallocate(dih_hist)
        end if
      end if
    end if
!
    if ((use_ZSEC.EQV..true.).OR.(segcalc.le.nsim).OR.(do_wanglandau.EQV..true.)) then
      if ((allocated(z_alpha).EQV..false.).OR.&
 &        (allocated(z_beta).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Z-sec structures were n&
 &ot properly allocated, or prematurely freed. Please fix or report &
 &this problem.'
        call fexit()
      else
        deallocate(z_alpha)
        deallocate(z_beta) 
      end if
    end if
    if (segcalc.le.nsim) then            
      if ((allocated(z_hist).EQV..false.).OR.&
 &        (allocated(z_hist2).EQV..false.).OR.&
 &        (allocated(z_2dhist).EQV..false.).OR.&
 &        (allocated(sgprscnt).EQV..false.).OR.&
 &        (allocated(seg_perrs).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Segment distribution ar&
 &rays were not properly allocated, or prematurely freed. Please fix&
 & or report this problem.'
        call fexit()
      else
        deallocate(z_hist)
        deallocate(z_hist2)
        deallocate(z_2dhist)
        deallocate(seg_perrs)
        deallocate(sgprscnt)
      end if
    end if
!
    if ((torlccalc.le.nsim).AND.(ntorlcs.gt.0)) then
      if ((allocated(torlc_coeff).EQV..false.).OR.&
 &        (allocated(torlc_coeff2).EQV..false.).OR.&
 &        (allocated(torlc_hist).EQV..false.).OR.&
 &        (allocated(lct_weight).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Torsional LC analysis a&
 &rrays were not properly allocated, or prematurely freed. Please fi&
 &x or report this problem.'
        call fexit()
      else
        deallocate(torlc_coeff)
        deallocate(torlc_coeff2)
        deallocate(torlc_hist)
        deallocate(lct_weight)
      end if
    end if
!
    if (covcalc.le.nsim) then
      if ((allocated(itrcv).EQV..false.).OR.&
 &        (allocated(trcvopen).EQV..false.).OR.&
 &        (allocated(ncovaravg).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. Covariance analysis arr&
 &ays were not properly allocated, or prematurely freed. Please fix &
 &or report this problem.'
        call fexit()
      else
        deallocate(itrcv)
        deallocate(ncovaravg)
        deallocate(trcvopen)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_torsn(...) with unknown mo&
 &de (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------------
!
subroutine allocate_gpcpre(mode)
!
  use iounit
  use paircorr
  use system
!
  implicit none
!
  integer mode
!
  if (pccalc.gt.nsim) return
!
! here we just allocate the data storage structure for the GPC requests
  if (mode.eq.1) then
    allocate(gpc%lst(gpc%all_lst,3))
    allocate(gpc%lty(gpc%all_lst))
    allocate(gpc%navg(gpc%all_lst))
    gpc%navg(:) = 0
  else if (mode.eq.2) then
    if ((allocated(gpc%lst).EQV..false.).OR.&
 &      (allocated(gpc%lty).EQV..false.).OR.&
 &      (allocated(gpc%navg).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. GPC request storage struc&
 &ture was not properly allocated, or prematurely freed. Please fix &
 &or report this problem.'
      call fexit()
    else
      deallocate(gpc%lst)
      deallocate(gpc%lty)
      deallocate(gpc%navg)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_gpcpre(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------------
!
subroutine allocate_pcs(mode)
!
  use iounit
  use sequen
  use system
  use molecule
  use paircorr
!
  implicit none
!
  integer i,j,mode,bba,imol,sca
!
  if (pccalc.gt.nsim) return
!
! first set up for whether specific amide-PC applies
  ampc%do_bb = .false.
  ampc%do_ss = .false.
  ampc%do_bs = .false.
  do imol=1,nmol
    if (do_amidepc.EQV..false.) exit
    sca = 0
    bba = 0
    do i=rsmol(imol,1),rsmol(imol,2)
      if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.22)).OR.&
 &         (seqtyp(i).eq.24).OR.&
 &        ((seqtyp(i).ge.27).AND.(seqtyp(i).le.30)).OR.&
 &        ((seqtyp(i).ge.33).AND.(seqtyp(i).le.35))) then
        bba = bba + 1
      end if
      if ((seqtyp(i).eq.17).OR.(seqtyp(i).eq.19)) then
        sca = sca + 1
      end if
    end do
    if (sca.ge.2) then
      ampc%do_ss = .true.
      ampc%do_bs = .true.
    else if (sca.eq.1) then
      ampc%do_bs = .true.
    end if
    if (bba.ge.2) then
      ampc%do_bb = .true.
    end if
  end do
!
! now (de)allocate
  if (mode.eq.1) then
!   amide-PC
    if (ampc%do_bb.EQV..true.) then
      allocate(ampc%bb(2,MAXPCBINS))
      do i=1,MAXPCBINS
        ampc%bb(1,i) = 0.0
        ampc%bb(2,i) = 0.0
      end do
    end if
    if (ampc%do_ss.EQV..true.) then
      allocate(ampc%ss(3,MAXPCBINS))
      do i=1,MAXPCBINS
        ampc%ss(1,i) = 0.0
        ampc%ss(2,i) = 0.0
        ampc%ss(3,i) = 0.0
      end do
    end if
    if (ampc%do_bs.EQV..true.) then
      allocate(ampc%bs(4,MAXPCBINS))
      do i=1,MAXPCBINS
        ampc%bs(1,i) = 0.0
        ampc%bs(2,i) = 0.0
        ampc%bs(3,i) = 0.0
        ampc%bs(4,i) = 0.0
      end do
    end if
!   general-PC
    if (gpc%nos.gt.0) then
      allocate(gpc%pc(gpc%nos,MAXPCBINS))
      do i=1,MAXPCBINS
        do j=1,gpc%nos
          gpc%pc(j,i) = 0.0
        end do
      end do
    end if
!   RBC-PC
    if (nsolutes.gt.1) then
      allocate(rbc_pc(nangrps,nangrps,MAXPCBINS))
      allocate(npcavg(nangrps,nangrps))
      npcavg(:,:) = 0
      rbc_pc(:,:,:) = 0.0
    end if
!   volume element
    allocate(voli_pc(MAXPCBINS))
    voli_pc(:) = 0.0
  else if (mode.eq.2) then
    if (ampc%do_bb.EQV..true.) then
      if (allocated(ampc%bb).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. Amide-BBBB array was no&
 &t properly allocated, or prematurely freed. Please fix or report t&
 &his problem.'
        call fexit()
      else
        deallocate(ampc%bb)
      end if
    end if
    if (ampc%do_bs.EQV..true.) then
      if (allocated(ampc%bs).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. Amide-BBSC array was no&
 &t properly allocated, or prematurely freed. Please fix or report t&
 &his problem.'
        call fexit()
      else
        deallocate(ampc%bs)
      end if
    end if
    if (ampc%do_ss.EQV..true.) then
      if (allocated(ampc%ss).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. Amide-SCSC array was no&
 &t properly allocated, or prematurely freed. Please fix or report t&
 &his problem.'
        call fexit()
      else
        deallocate(ampc%ss)
      end if
    end if
    if (gpc%nos.gt.0) then
      if (allocated(gpc%pc).EQV..false.) then
        write(ilog,*) 'Fatal. This is a bug. General-PC array was no&
 &t properly allocated, or prematurely freed. Please fix or report t&
 &his problem.'
        call fexit()
      else
        deallocate(gpc%pc)
      end if
    end if
    if (nsolutes.gt.1) then
      if ((allocated(rbc_pc).EQV..false.).OR.&
 &        (allocated(npcavg).EQV..false.)) then
        write(ilog,*) 'Fatal. This is a bug. RBC-PC arrays were not pr&
 &operly allocated, or prematurely freed. Please fix or report this &
 &problem.'
        call fexit()
      else
        deallocate(rbc_pc)
        deallocate(npcavg)
      end if
    end if
    if (allocated(voli_pc).EQV..false.) then
      write(ilog,*) 'Fatal. This is a bug. PC volume element array was not pr&
 &operly allocated, or prematurely freed. Please fix or report this &
 &problem.'
      call fexit()
    else
      deallocate(voli_pc)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_pcs(...) with unknown mode&
 & (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------------
!
subroutine allocate_tbppre(mode)
!
  use iounit
  use tabpot
  use system
  use energies
!
  implicit none
!
  integer mode
!
  if (use_TABUL.EQV..false.) return
!
! here we just allocate the data storage structure for the GPC requests
  if (mode.eq.1) then
    allocate(tbp%lst(tbp%all_lst,3))
    allocate(tbp%lty(tbp%all_lst))
  else if (mode.eq.2) then
    if ((allocated(tbp%lst).EQV..false.).OR.&
 &      (allocated(tbp%lty).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. TBP request storage struc&
 &ture was not properly allocated, or prematurely freed. Please fix &
 &or report this problem.'
      call fexit()
    else
      deallocate(tbp%lst)
      deallocate(tbp%lty)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_tbppre(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------------
!
subroutine allocate_tabpot(mode)
!
  use iounit
  use tabpot
  use system
  use energies
  use sequen
!
  implicit none
!
  integer mode
!
  if (use_TABUL.EQV..false.) return
!
! here we just allocate the data storage structure for the actual TBPs
  if (mode.eq.1) then
    allocate(tbp%pot(tbp%nos,tbp%bins))
    allocate(tbp%tang(tbp%nos,tbp%bins))
    allocate(tbp%dis(tbp%bins))
    allocate(tbp%rsmat(nseq,nseq))
    allocate(tbp%rsvec(nseq))
  else if (mode.eq.2) then
    if ((allocated(tbp%pot).EQV..false.).OR.&
 &      (allocated(tbp%dis).EQV..false.).OR.&
 &      (allocated(tbp%tang).EQV..false.).OR.&
 &      (allocated(tbp%rsmat).EQV..false.).OR.&
 &      (allocated(tbp%rsvec).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. TBP data structure was no&
 &t properly allocated, or prematurely freed. Please fix or report t&
 &his problem.'
      call fexit()
    else
      deallocate(tbp%pot)
      deallocate(tbp%tang)
      deallocate(tbp%dis)
      deallocate(tbp%rsmat)
      deallocate(tbp%rsvec)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_tabpot(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_ionize(mode)
!
  use iounit
  use ionize
  use movesets
!
  implicit none
!
  integer mode
!
  if (phfreq.le.0.0) return
!
  if (mode.eq.1) then
    allocate(ionres(nionres))
    allocate(titsite(nionres))
    allocate(ionstate(nionres))
    allocate(pkai(nionres))
    allocate(prtionst(nionres))
!
  else if (mode.eq.2) then
    if ((allocated(ionres).EQV..false.).OR.&
 &      (allocated(titsite).EQV..false.).OR.&
 &      (allocated(ionstate).EQV..false.).OR.&
 &      (allocated(prtionst).EQV..false.).OR.&
 &      (allocated(pkai).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Titration move arrays wer&
 &e not properly allocated, or prematurely freed. Please fix or repo&
 &rt this problem.'
      call fexit()
    else
      deallocate(ionres)
      deallocate(titsite)
      deallocate(ionstate) 
      deallocate(prtionst)
      deallocate(pkai)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_ionize(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_distrest(mode)
!
  use iounit
  use distrest
!
  implicit none
!
  integer mode
!
  if (allndr.le.0) return
!
  if (mode.eq.1) then
    allocate(dresat(allndr,2))
    allocate(dresprm(allndr,2))
    allocate(drestyp(allndr))
!
  else if (mode.eq.2) then
    if ((allocated(dresat).EQV..false.).OR.&
 &      (allocated(dresprm).EQV..false.).OR.&
 &      (allocated(drestyp).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Distance restraint terms &
 &were not properly allocated, or prematurely freed. Please fix or r&
 &eport this problem.'
      call fexit()
    else
      deallocate(dresat)
      deallocate(dresprm)
      deallocate(drestyp)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_distrest(...) with unknown&
 & mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_mpistuff(mode)
!
  use iounit
  use mpistuff
  use system
!
  implicit none
!
  integer mode
!
#ifdef ENABLE_MPI
  integer i,j,masterrank
!
  masterrank = 0
!
  if (mode.eq.1) then
    if (use_REMC.EQV..false.) return
    allocate(re_trans(mpi_nodes,mpi_nodes))
    allocate(re_probs(mpi_nodes,mpi_nodes))
    allocate(re_olap(mpi_nodes,4))
    allocate(re_mat(mpi_nodes,MAXREDIMS))
    allocate(mpi_lmap(2*mpi_nodes,2))
    do i=1,mpi_nodes
      mpi_lmap(i,2) = i
    end do
    do i=1,mpi_nodes
      do j=1,MAXREDIMS
        re_mat(i,j) = 273.0 + (i-1)*5.0
      end do
    end do
    do i=1,mpi_nodes
      do j=1,mpi_nodes
        re_trans(i,j) = 0
        re_probs(i,j) = 0.0
      end do
      do j=1,4
        re_olap(i,j) = 0.0
      end do
    end do
!
  else if (mode.eq.2) then
    if ((use_REMC.EQV..true.).AND.((allocated(re_trans).EQV..false.).OR.&
 &      (allocated(re_probs).EQV..false.).OR.&
 &      (allocated(re_olap).EQV..false.).OR.&
 &      ((allocated(mpi_lmap).EQV..false.).AND.(myrank.eq.masterrank)).OR.&
 &      (allocated(re_mat).EQV..false.))) then
      write(ilog,*) 'Fatal. This is a bug. REMC-arrays were not prop&
 &erly allocated, or prematurely freed. Please fix or report this pr&
 &oblem.'
      call fexit()
    else
      if (allocated(mpi_granullst%sends).EQV..true.) deallocate(mpi_granullst%sends)
      if (allocated(mpi_granullst%sends2).EQV..true.) deallocate(mpi_granullst%sends2)
      if (allocated(mpi_granullst%recvs).EQV..true.) deallocate(mpi_granullst%recvs)
      if (allocated(mpi_granullst%recvs2).EQV..true.) deallocate(mpi_granullst%recvs2)
      if (allocated(mpi_globlst%sends).EQV..true.) deallocate(mpi_globlst%sends)
      if (allocated(mpi_globlst%sends2).EQV..true.) deallocate(mpi_globlst%sends2)
      if (allocated(mpi_globlst%recvs).EQV..true.) deallocate(mpi_globlst%recvs)
      if (allocated(mpi_globlst%recvs2).EQV..true.) deallocate(mpi_globlst%recvs2)
      if (use_REMC.EQV..true.) then
        deallocate(re_trans)
        deallocate(re_probs)
        deallocate(re_mat)
        deallocate(re_olap)
        deallocate(mpi_lmap)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_mpistuff(...) with unknown&
 & mode (offending mode is ',mode,').'
    call fexit()
  end if
!
#endif
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_params(mode)
!
  use iounit
  use params
!
  implicit none
!
  integer i,j,mode
!
  if (mode.eq.1) then
    allocate(lj_atnum(n_ljtyp))
    allocate(lj_val(n_ljtyp))
    allocate(lj_weight(n_ljtyp))
    allocate(lj_symbol(n_ljtyp))
    allocate(lj_describe(n_ljtyp))
    allocate(lj_eps(n_ljtyp,n_ljtyp))
    allocate(lj_sig(n_ljtyp,n_ljtyp))
    allocate(lj_eps_14(n_ljtyp,n_ljtyp))
    allocate(lj_sig_14(n_ljtyp,n_ljtyp))
    allocate(lj_rad(n_ljtyp))
    do i=1,n_ljtyp
      lj_symbol(i) = "   "
      lj_describe(i) = "                    "
      lj_val(i) = -1
      lj_atnum(i) = -1
      lj_weight(i) = -1.0d0
      do j=1,n_ljtyp
        lj_eps(i,j) = -1.0
        lj_sig(i,j) = -1.0
        lj_eps_14(i,j) = -1.0
        lj_sig_14(i,j) = -1.0
      end do
    end do
!
    allocate(bio_ljtyp(n_biotyp))
    allocate(bio_ctyp(n_biotyp))
    allocate(bio_botyp(n_biotyp))
    allocate(bio_code(n_biotyp))
    allocate(bio_res(n_biotyp))
    do i=1,n_biotyp
      bio_code(i) = "   "
      bio_res(i) = "                    "
      bio_ljtyp(i) = -1
      bio_ctyp(i) = -1
      bio_botyp(i) = -1
    end do
!
    allocate(c_charge(n_ctyp))
    allocate(c_note(n_ctyp))
    do i=1,n_ctyp
      c_charge(i) = 0.0d0
      c_note(i) = "                    "
    end do
!
    allocate(bo_typ(n_bondtyp))
    allocate(bo_par(n_bondtyp,MAXBOPAR))
    do i=1,n_bondtyp
      bo_typ(i) = -1
    end do
    bo_par(:,:) = 0.0
!
    allocate(ba_typ(n_angltyp))
    allocate(ba_par(n_angltyp,MAXBAPAR))
    do i=1,n_angltyp
      ba_typ(i) = -1
    end do
    ba_par(:,:) = 0.0
!
    allocate(di_typ(n_torstyp))
    allocate(di_par(n_torstyp,MAXDIPAR))
    do i=1,n_torstyp
      di_typ(i) = -1
    end do
    di_par(:,:) = 0.0
!
    allocate(cm_typ(n_cmstyp,2))
    allocate(cm_par(n_cmstyp,maxcmpar,maxcmpar,17))
    allocate(cm_par2(n_cmstyp))
    allocate(cm_file(n_cmstyp))
    do i=1,n_cmstyp
      cm_typ(i,1) = -1
      cm_typ(i,2) = 4
      cm_file(i) = "CMAP_01.dat"
    end do
    cm_par(:,:,:,:) = 0.0
    cm_par2(:) = 0.0
!
  else if (mode.eq.2) then
    if ((allocated(lj_atnum).EQV..false.).OR.&
 &      (allocated(lj_val).EQV..false.).OR.&
 &      (allocated(lj_weight).EQV..false.).OR.&
 &      (allocated(lj_symbol).EQV..false.).OR.&
 &      (allocated(lj_describe).EQV..false.).OR.&
 &      (allocated(lj_eps).EQV..false.).OR.&
 &      (allocated(lj_sig).EQV..false.).OR.&
 &      (allocated(lj_eps_14).EQV..false.).OR.&
 &      (allocated(lj_sig_14).EQV..false.).OR.&
 &      (allocated(lj_rad).EQV..false.).OR.&
 &      (allocated(bio_code).EQV..false.).OR.&
 &      (allocated(bio_res).EQV..false.).OR.&
 &      (allocated(bio_ljtyp).EQV..false.).OR.&
 &      (allocated(bio_ctyp).EQV..false.).OR.&
 &      (allocated(bio_botyp).EQV..false.).OR.&
 &      (allocated(c_charge).EQV..false.).OR.&
 &      (allocated(c_note).EQV..false.).OR.&
 &      (allocated(bo_typ).EQV..false.).OR.&
 &      (allocated(ba_typ).EQV..false.).OR.&
 &      (allocated(di_typ).EQV..false.).OR.&
 &      (allocated(cm_typ).EQV..false.).OR.&
 &      (allocated(bo_par).EQV..false.).OR.&
 &      (allocated(ba_par).EQV..false.).OR.&
 &      (allocated(di_par).EQV..false.).OR.&
 &      (allocated(cm_par).EQV..false.).OR.&
 &      (allocated(cm_file).EQV..false.).OR.&
 &      (allocated(cm_par2).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Parameter storage arrays &
 &were not properly allocated, or prematurely freed. Please fix or r&
 &eport this problem.'
      call fexit()
    else
      deallocate(lj_atnum)
      deallocate(lj_val)
      deallocate(lj_symbol)
      deallocate(lj_weight)
      deallocate(lj_describe)
      deallocate(lj_eps)
      deallocate(lj_sig)
      deallocate(lj_eps_14)
      deallocate(lj_sig_14)
      deallocate(lj_rad)
      deallocate(bio_code)
      deallocate(bio_res)
      deallocate(bio_ljtyp)
      deallocate(bio_ctyp)
      deallocate(bio_botyp)
      deallocate(c_charge)
      deallocate(c_note)
      deallocate(bo_typ)
      deallocate(ba_typ)
      deallocate(di_typ)
      deallocate(cm_typ)
      deallocate(bo_par)
      deallocate(ba_par)
      deallocate(di_par)
      deallocate(cm_par2)
      deallocate(cm_par)
      deallocate(cm_file)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_params(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_dipolavg(mode)
!
  use iounit
  use dipolavg
  use molecule
  use sequen
  use system
!
  implicit none
!
  integer mode
!
  if (dipcalc.gt.nsim) return
!
  if (mode.eq.1) then
    allocate(nmtdavg(nangrps))
    allocate(nrsdavg(nseq))
    allocate(mtdavg(nangrps,DIPDIM))
    allocate(rsdavg(nseq,DIPDIM))
    rsdavg(:,:) = 0.0
    mtdavg(:,:) = 0.0
    nmtdavg(:) = 0
    nrsdavg(:) = 0
!
  else if (mode.eq.2) then
    if ((allocated(nmtdavg).EQV..false.).OR.&
 &      (allocated(nrsdavg).EQV..false.).OR.&
 &      (allocated(mtdavg).EQV..false.).OR.&
 &      (allocated(rsdavg).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Dipole analysis arrays we&
 &re not properly allocated, or prematurely freed. Please fix or rep&
 &ort this problem.'
      call fexit()
    else
      deallocate(nmtdavg)
      deallocate(nrsdavg)
      deallocate(mtdavg)
      deallocate(rsdavg)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_dipolavg(...) with unknown&
 & mode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_forces(mode)
!
  use iounit
  use forces
  use atoms
  use molecule
  use system
!
  implicit none
!
  integer alcsz,mode,imol
!
!  if (use_FORCES.EQV..false.) return
!
  if (mode.eq.1) then
    allocate(cart_f(n,3))
    allocate(cart_f_tr(n,3))
    allocate(cart_v(n,3))
    cart_v(:,:) = 0.0
    allocate(cart_a(n,3))
    allocate(cart_frz(n,3))
    cart_frz(:,:) = .false.
    allocate(sav_dr(3,n))
    allocate(scrq_dr(3,n))
    scrq_dr(:,:) = 0.0
    allocate(sisa(n))
    sisa(1:n)%nix = 0
    sisa(1:n)%alsz = 0
    allocate(sum_scrcbs(n))
    allocate(sum_scrcbs_tr(n))
    allocate(sum_scrcbs_lr(n))
    allocate(cart_ldp(n,3,5))
    allocate(tstat%molgrp(nmol+1)) ! +1 to accomodate boundary particle in some NPT
    allocate(tstat%grpT(tstat%n_tgrps)) ! note that this is de-/re-allocated later (see parsefiles.f)
    allocate(tstat%grpdof(tstat%n_tgrps)) ! note that this is de-/re-allocated later (see parsefiles.f)
    tstat%molgrp(:) = 1
    allocate(dc_di(nmol))
    do imol=1,nmol
      if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
        alcsz = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        alcsz = 5
      else
        alcsz = ntormol(moltypid(imol))+6
      end if
      allocate(dc_di(imol)%f(alcsz))
      allocate(dc_di(imol)%im(alcsz))
      allocate(dc_di(imol)%olddat(alcsz,3))
      allocate(dc_di(imol)%v(alcsz))
      allocate(dc_di(imol)%avgT(alcsz))
      allocate(dc_di(imol)%frz(alcsz))
      allocate(dc_di(imol)%align(alcsz))
      allocate(dc_di(imol)%incr(alcsz))
      allocate(dc_di(imol)%nattor(alcsz))
      allocate(dc_di(imol)%ldp(alcsz,5))
      if (fudge_mass.gt.0) then
        allocate(dc_di(imol)%fudge(alcsz))
        dc_di(imol)%fudge(:) = 1.0
      end if
      dc_di(imol)%maxntor = 0
      dc_di(imol)%avgT(:) = 0.0
      dc_di(imol)%frz(:) = .false.
      dc_di(imol)%align(:) = .false.
      dc_di(imol)%nattor(:) = 0
      dc_di(imol)%incr(:) = 0.0
      dc_di(imol)%ldp(:,:) = 0.0
      dc_di(imol)%im(:) = 0.0
      dc_di(imol)%olddat(:,:) = 0.0
      dc_di(imol)%v(:) = 0.0
      dc_di(imol)%f(:) = 0.0
    end do
!
  else if (mode.eq.2) then
    if ((allocated(cart_f).EQV..false.).OR.&
 &      (allocated(cart_f_tr).EQV..false.).OR.&
 &      (allocated(cart_v).EQV..false.).OR.&
 &      (allocated(cart_a).EQV..false.).OR.&
 &      (allocated(cart_frz).EQV..false.).OR.&
 &      (allocated(sav_dr).EQV..false.).OR.&
 &      (allocated(scrq_dr).EQV..false.).OR.&
 &      (allocated(sisa).EQV..false.).OR.&
 &      (allocated(sum_scrcbs).EQV..false.).OR.&
 &      (allocated(sum_scrcbs_tr).EQV..false.).OR.&
 &      (allocated(sum_scrcbs_lr).EQV..false.).OR.&
 &      (allocated(cart_ldp).EQV..false.).OR.&
 &      (allocated(dc_di).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. Force arrays were not pro&
 &perly allocated, or prematurely freed. Please fix or report or thi&
 &s problem.'
      call fexit()
    else
      deallocate(cart_f)
      deallocate(cart_f_tr)
      deallocate(cart_v)
      deallocate(cart_a)
      deallocate(cart_frz)
      deallocate(sav_dr)
      deallocate(scrq_dr)
      do imol=1,n
        if (allocated(sisa(imol)%ix).EQV..true.) deallocate(sisa(imol)%ix)
        if (allocated(sisa(imol)%dr).EQV..true.) deallocate(sisa(imol)%dr)
        if (allocated(sisa(imol)%qr).EQV..true.) deallocate(sisa(imol)%qr)
      end do
      deallocate(sisa)
      deallocate(sum_scrcbs)
      deallocate(sum_scrcbs_tr)
      deallocate(sum_scrcbs_lr)
      deallocate(cart_ldp)
      do imol=1,nmol
        if (allocated(dc_di(imol)%f).EQV..true.) deallocate(dc_di(imol)%f)
        if (allocated(dc_di(imol)%im).EQV..true.) deallocate(dc_di(imol)%im)
        if (allocated(dc_di(imol)%olddat).EQV..true.) deallocate(dc_di(imol)%olddat)
        if (allocated(dc_di(imol)%v).EQV..true.) deallocate(dc_di(imol)%v)
        if (allocated(dc_di(imol)%avgT).EQV..true.) deallocate(dc_di(imol)%avgT)
        if (allocated(dc_di(imol)%ldp).EQV..true.) deallocate(dc_di(imol)%ldp)
        if (allocated(dc_di(imol)%fudge).EQV..true.) deallocate(dc_di(imol)%fudge)
        if (allocated(dc_di(imol)%frz).EQV..true.) deallocate(dc_di(imol)%frz)
        if (allocated(dc_di(imol)%align).EQV..true.) deallocate(dc_di(imol)%align)
        if (allocated(dc_di(imol)%incr).EQV..true.) deallocate(dc_di(imol)%incr)
        if (allocated(dc_di(imol)%nattor).EQV..true.) deallocate(dc_di(imol)%nattor)
        if (allocated(dc_di(imol)%recurs).EQV..true.) deallocate(dc_di(imol)%recurs)
        if (allocated(dc_di(imol)%valrecurs).EQV..true.) deallocate(dc_di(imol)%valrecurs)
      end do
      deallocate(dc_di)
      if (allocated(tstat%molgrp).EQV..true.) deallocate(tstat%molgrp)
      if (allocated(tstat%grpT).EQV..true.) deallocate(tstat%grpT)
      if (allocated(tstat%grpdof).EQV..true.) deallocate(tstat%grpdof)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_forces(...) with unknown m&
 &ode (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!------------------------------------------------------------------------------
!
subroutine allocate_diffrac(mode)
!
  use iounit
  use system
  use diffrac
!
  implicit none
!
  integer i,j,mode
!
  if (diffrcalc.gt.nsim) return
!
  if (mode.eq.1) then
    allocate(am_diffr_map(rr_max,zz_max))
    do i=1,rr_max
      do j=1,zz_max
        am_diffr_map(i,j) = 0.0
      end do
    end do
!
  else if (mode.eq.2) then
    if (allocated(am_diffr_map).EQV..false.) then
      write(ilog,*) 'Fatal. This is a bug. Diffraction analysis arra&
 &ys were not properly allocated, or prematurely freed. Please fix o&
 &r report this problem.'
      call fexit()
    else
      deallocate(am_diffr_map)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_diffr(...) with unknown mo&
 &de (offending mode is ',mode,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------------
!
subroutine allocate_particlefluc(mode)
!
  use molecule
  use grandensembles
  use iounit
!
  implicit none
!
  integer mode, typid
!  
  if (mode.eq.1) then
    call constructset(fluctypes, nmoltyp)
    call constructset(ispresent, nmol)
    allocate(typesets(nmoltyp,2))
    do typid = 1, nmoltyp
      call constructset(typesets(typid,1), nmol)
      call constructset(typesets(typid,2), nmol)
    end do
    allocate(chempot(nmoltyp))
    allocate(thermalvolume(nmoltyp))
    allocate(eqnum(nmoltyp))
    allocate(numberhistogram(nmoltyp*3,0:nmol))
    numberhistogram(:,:) = 0
  else if (mode.eq.2) then
    if (allocated(numberhistogram).EQV..true.) then
      deallocate(numberhistogram)
    end if
    if (allocated(thermalvolume).EQV..true.) then
      deallocate(thermalvolume)
    end if
    if (allocated(chempot).EQV..true.) then
      deallocate(chempot)
    end if
    if (allocated(eqnum).EQV..true.) then
      deallocate(eqnum)
    end if
    if (allocated(typesets).EQV..true.) then
      do typid = 1, nmoltyp
        call destroyset(typesets(typid,1))
        call destroyset(typesets(typid,2))
      end do
      deallocate(typesets)
    end if
    call destroyset(ispresent)
    call destroyset(fluctypes)
  else
    write(ilog,*) 'Fatal. Called allocate_particlefluc with unknown &
 &mode (offending mode is ', mode, ').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine allocate_dssps(mode)
!
  use iounit
  use sequen
  use dssps
  use system
  use molecule
  use energies
!
  implicit none
!
  integer mode
!
  if ((dsspcalc.gt.nsim).AND.(use_DSSP.EQV..false.)) return
!
  if (mode.eq.1) then
    allocate(peprs(nseq))
    allocate(rspep(pep_sz))
    allocate(hb_mat(pep_sz,pep_sz))
    allocate(hbl(pep_sz))
    allocate(pepmolsz(npepmol))
    allocate(dssp_cnt(npepmol))
    allocate(pepmol(nmol))
    allocate(molpep(npepmol,2))
    allocate(bdgtb%lims(2*nseq,4))
    allocate(bdgtb%ty(2*nseq))
    allocate(bdgtb%phbe(2*nseq))
    allocate(dssp_ass(pep_sz,8))
    allocate(dssp_perrs(pep_sz,11,pep_sz))
    allocate(dssp_2dhist(100,100))
    allocate(dssp_hists(4,100))
    allocate(whsc_der(3,4,pep_sz))
    allocate(wesc_der(3,4,pep_sz))
    peprs(:) = 0
    pepmolsz(:) = 0
    pepmol(:) = 0
    hb_mat(:,:) = 0
    hbl(:)%nahbs = 0
    hbl(:)%ndhbs = 0
    dssp_ass(:,:) = 0
    dssp_cnt(:) = 0
    dssp_perrs(:,:,:) = 0.0
    dssp_hists(:,:) = 0.0
    dssp_2dhist(:,:) = 0.0
  else if (mode.eq.2) then
    if ((allocated(peprs).EQV..false.).OR.&
 &      (allocated(hb_mat).EQV..false.).OR.&
 &      (allocated(hbl).EQV..false.).OR.&
 &      (allocated(rspep).EQV..false.).OR.&
 &      (allocated(pepmolsz).EQV..false.).OR.&
 &      (allocated(pepmol).EQV..false.).OR.&
 &      (allocated(molpep).EQV..false.).OR.&
 &      (allocated(dssp_ass).EQV..false.).OR.&
 &      (allocated(dssp_cnt).EQV..false.).OR.&
 &      (allocated(dssp_hists).EQV..false.).OR.&
 &      (allocated(dssp_2dhist).EQV..false.).OR.&
 &      (allocated(dssp_perrs).EQV..false.).OR.&
 &      (allocated(whsc_der).EQV..false.).OR.&
 &      (allocated(wesc_der).EQV..false.).OR.&
 &      (allocated(bdgtb%lims).EQV..false.).OR.&
 &      (allocated(bdgtb%phbe).EQV..false.).OR.&
 &      (allocated(bdgtb%ty).EQV..false.)) then
      write(ilog,*) 'Fatal. This is a bug. DSSP analysis arrays were&
 & not properly allocated, or prematurely freed. Please fix or repor&
 &t this problem.'
      call fexit()
    else
      deallocate(peprs)
      deallocate(hb_mat)
      deallocate(hbl)
      deallocate(rspep)
      deallocate(pepmolsz)
      deallocate(pepmol)
      deallocate(molpep)
      deallocate(dssp_ass)
      deallocate(dssp_cnt)
      deallocate(dssp_2dhist)
      deallocate(dssp_hists)
      deallocate(dssp_perrs)
      deallocate(wesc_der)
      deallocate(whsc_der)
      deallocate(bdgtb%lims)
      deallocate(bdgtb%phbe)
      deallocate(bdgtb%ty)
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_dssps with unknown mode (o&
 &ffending mode is ',mode, ').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine allocate_clusters(mode)
!
  use iounit
  use system
  use clusters
!
  implicit none
!
  integer mode,i,j
!
  if (cstorecalc.gt.nsim) return
!
  if (mode.eq.1) then
!   do nothing
  else if (mode.eq.2) then
    if (allocated(cdofset).EQV..false.) then
      write(ilog,*) 'Fatal. This is a bug. Structural clustering arrays were&
 & not properly allocated, or prematurely freed. Please fix or repor&
 &t this problem.'
      call fexit()
    else
      if (allocated(scluster).EQV..true.) then
        do i=1,size(scluster)
          if (allocated(scluster(i)%snaps).EQV..true.) deallocate(scluster(i)%snaps)
          if (allocated(scluster(i)%tmpsnaps).EQV..true.) deallocate(scluster(i)%tmpsnaps)
          if (allocated(scluster(i)%sums).EQV..true.) deallocate(scluster(i)%sums)
          if (allocated(scluster(i)%map).EQV..true.) deallocate(scluster(i)%map)
          if (allocated(scluster(i)%children).EQV..true.) deallocate(scluster(i)%children)
          if (allocated(scluster(i)%wghtsnb).EQV..true.) deallocate(scluster(i)%wghtsnb)
          if (allocated(scluster(i)%lstnb).EQV..true.) deallocate(scluster(i)%lstnb)
          if (allocated(scluster(i)%flwnb).EQV..true.) deallocate(scluster(i)%flwnb)
          if (allocated(scluster(i)%lensnb).EQV..true.) deallocate(scluster(i)%lensnb)
          if (allocated(scluster(i)%fewtsnb).EQV..true.) deallocate(scluster(i)%fewtsnb)
        end do
        deallocate(scluster)
      end if
      if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
      if (allocated(cl_mwvec).EQV..true.) deallocate(cl_mwvec)
      if (allocated(cludata).EQV..true.) deallocate(cludata)
      if (allocated(clufold).EQV..true.) deallocate(clufold)
      if (allocated(cluunfold).EQV..true.) deallocate(cluunfold)
      if (allocated(erp_msm_fwr).EQV..true.) deallocate(erp_msm_fwr)
      if (allocated(eip_msm_bwr).EQV..true.) deallocate(eip_msm_bwr)
      if (allocated(eip_msm_fwr).EQV..true.) deallocate(eip_msm_fwr)
      if (allocated(eigvec_msm_fwr).EQV..true.) deallocate(eigvec_msm_fwr)
      if (allocated(eigvec_msm_bwr).EQV..true.) deallocate(eigvec_msm_bwr)
      if (allocated(pfolds).EQV..true.) deallocate(pfolds)
      if (allocated(cnblst).EQV..true.) then
        do i=1,size(cnblst)
          if (allocated(cnblst(i)%dis).EQV..true.) deallocate(cnblst(i)%dis)
          if (allocated(cnblst(i)%idx).EQV..true.) deallocate(cnblst(i)%idx)
          if (allocated(cnblst(i)%tagged).EQV..true.) deallocate(cnblst(i)%tagged)
        end do
        deallocate(cnblst)
      end if
      if (allocated(birchtree).EQV..true.) then
        do i=1,c_nhier+1
          do j=1,size(birchtree(i)%cls) ! birchtree(i)%ncls
            if (allocated(birchtree(i)%cls(j)%snaps).EQV..true.) deallocate(birchtree(i)%cls(j)%snaps)
            if (allocated(birchtree(i)%cls(j)%tmpsnaps).EQV..true.) deallocate(birchtree(i)%cls(j)%tmpsnaps)
            if (allocated(birchtree(i)%cls(j)%sums).EQV..true.) deallocate(birchtree(i)%cls(j)%sums)
            if (allocated(birchtree(i)%cls(j)%map).EQV..true.) deallocate(birchtree(i)%cls(j)%map)
            if (allocated(birchtree(i)%cls(j)%children).EQV..true.) deallocate(birchtree(i)%cls(j)%children)
            if (allocated(birchtree(i)%cls(j)%wghtsnb).EQV..true.) deallocate(birchtree(i)%cls(j)%wghtsnb)
            if (allocated(birchtree(i)%cls(j)%lstnb).EQV..true.) deallocate(birchtree(i)%cls(j)%lstnb)
            if (allocated(birchtree(i)%cls(j)%flwnb).EQV..true.) deallocate(birchtree(i)%cls(j)%flwnb)
            if (allocated(birchtree(i)%cls(j)%lensnb).EQV..true.) deallocate(birchtree(i)%cls(j)%lensnb)
            if (allocated(birchtree(i)%cls(j)%fewtsnb).EQV..true.) deallocate(birchtree(i)%cls(j)%fewtsnb)
          end do
          deallocate(birchtree(i)%cls)
        end do
        deallocate(birchtree)
      end if
      if (allocated(approxmst).EQV..true.) then
        do i=1,size(approxmst)
          if (allocated(approxmst(i)%adj).EQV..true.) deallocate(approxmst(i)%adj)
          if (allocated(approxmst(i)%dist).EQV..true.) deallocate(approxmst(i)%dist)
        end do
        deallocate(approxmst)
      end if
    end if
  else
    write(ilog,*) 'Fatal. Called allocate_clusters with unknown mode (o&
 &ffending mode is ',mode, ').'
    call fexit()
  end if
!
end
!
#ifdef ENABLE_MPI
!
!-----------------------------------------------------------------------
!
subroutine allocate_mpire(mode)
!
  use iounit
  use mpistuff
  use system
  use atoms
  use torsn
  use molecule
!
  implicit none
!
  integer mode,masterrank
!
  masterrank = 0
!
  if ((use_REMC.EQV..false.).AND.(use_MPIMultiSeed.EQV..false.)) return
!
  if (mode.eq.1) then
!
    allocate(mpi_evec(mpi_nodes))
    allocate(mpi_fvec(mpi_nodes))
!
    if (myrank.eq.masterrank) allocate(mpi_emat(mpi_nodes,mpi_nodes))
    if ((fycxyz.eq.2).OR.(force_rexyz.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      allocate(mpi_cdsvec(n,3,2))
    end if
    if (((dyn_mode.eq.1).AND.(fycxyz.eq.1).AND.(force_rexyz.EQV..false.)).OR.&
 &        ((use_dyn.EQV..true.).AND.(fycxyz.eq.1))) then
      allocate(mpi_tvec(ntorpuck,2))
      allocate(mpi_rbvec(nmol*9,2))
    end if
!     
    if (dyn_mode.eq.1) then
!     do nothing 
    else if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
      if (fycxyz.eq.2) then
        allocate(mpi_dvec(n,3,2))
      else
        allocate(mpi_dvec(totrbd+ndyntorsn,3,2))
      end if
    else if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
      if (fycxyz.eq.2) then
        allocate(mpi_dvec(n,9,2))
      else
        allocate(mpi_dvec(totrbd+ndyntorsn,4,2))
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported dynamics method (',dyn_mode,') in&
 & allocate_mpire(...). This is most certainly a bug.'
      call fexit()
    end if
  else if (mode.eq.2) then
    if (allocated(mpi_evec).EQV..true.) deallocate(mpi_evec)
    if (allocated(mpi_emat).EQV..true.) deallocate(mpi_emat)
    if (allocated(mpi_fvec).EQV..true.) deallocate(mpi_fvec)
    if (allocated(mpi_tvec).EQV..true.) deallocate(mpi_tvec)
    if (allocated(mpi_rbvec).EQV..true.) deallocate(mpi_rbvec)
    if (allocated(mpi_cdsvec).EQV..true.) deallocate(mpi_cdsvec)
    if (allocated(mpi_dvec).EQV..true.) deallocate(mpi_dvec)
  else
    write(ilog,*) 'Fatal. Called allocate_mpire with unknown mode (o&
 &ffending mode is ',mode, ').'
    call fexit()
  end if
!
end
!
#endif
!
!-----------------------------------------------------------------------------
!
subroutine allocate_ewalds(mode)
!
  use iounit
  use ewalds
#ifdef LINK_FFTW
  use, intrinsic :: ISO_C_BINDING
#endif
!
  implicit none
!
#ifdef LINK_FFTW
#ifdef ORACLE_FORTRAN
#include "fftw3_SUN.f03"
#else
#include "fftw3.f03"
#endif
#endif
!
  integer mode
!
  if (mode.eq.1) return
!
  if (mode.eq.2) then
!
    if (allocated(bsmx).EQV..true.) deallocate(bsmx)
    if (allocated(bsmy).EQV..true.) deallocate(bsmy)
    if (allocated(bsmz).EQV..true.) deallocate(bsmz)
!    if (allocated(Qew1).EQV..true.) deallocate(Qew1)
!    if (allocated(Qew2).EQV..true.) deallocate(Qew2)
#ifdef LINK_FFTW
    call fftw_free(pqew1)
    call fftw_free(pqew2)
#endif
    if (allocated(QewBC).EQV..true.) deallocate(QewBC)
    if (allocated(bspl).EQV..true.) deallocate(bspl)
    if (allocated(bspld).EQV..true.) deallocate(bspld)
    if (allocated(ewgfls).EQV..true.) deallocate(ewgfls)
    if (allocated(ewSfac).EQV..true.) deallocate(ewSfac)
    if (allocated(ewSfacs).EQV..true.) deallocate(ewSfacs)
    if (allocated(ewSvd).EQV..true.) deallocate(ewSvd)
    if (allocated(ewSvds).EQV..true.) deallocate(ewSvds)
!
  else
    write(ilog,*) 'Fatal. Called allocate_ewalds(...) with unknown mode (o&
 &ffending mode is ',mode, ').'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------------------------
!
