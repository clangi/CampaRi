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
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!--------------------------------------------------------------------------------------------
!
subroutine gen_graph_from_clusters(it,nbasins,prtclu,tpi)
!
  use interfaces
  use clusters
  use mpistuff
  use iounit
  use system
  use threads
!
  implicit none
!
  integer, INTENT(IN):: nbasins,tpi
  logical, INTENT(IN):: prtclu
!
  integer i,j,k,iu,freeunit,ii,jj,ll,kk,ik,jk,inissnap,endssnap,athree,azero,aone,ssix
  integer, ALLOCATABLE:: ivx(:,:)
  character(MAXSTRLEN) fn
  logical foundit,exists,dofwr
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
  type(t_scluster) it(nbasins)
  RTYPE dis,dval
  character(max(1,max(1,ceiling(log10(1.0*ceiling(log10(1.0*nbasins+0.5))))))) bnrs
  character(max(1,ceiling(log10(1.0*(1+c_multires)+0.5)))) bnrs2
  character(MAXSTRLEN) fmtstr
  character(3) tmsfx
  integer(KIND=8) tcnt(20)
!
! init
  azero = 0
  aone = 1
  athree = 3
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  call System_Clock(count=tcnt(1))
  write(ilog,*)
  write(ilog,'(1x,a)') 'Now constructing graph from clustering ...'
  if ((prtclu.EQV..true.).AND.(cmode.ne.6).AND.(cmode.ne.7)) then
!
    ik = 0
    jk = max(1,ceiling(log10(1.0*nbasins+0.5)))
    ii = max(1,ceiling(log10(1.0*jk+0.5)))
    call int2str(jk,bnrs,ik)
    ik = 0
    jk = 1+c_multires
    jj = max(1,ceiling(log10(1.0*jk+0.5)))
    call int2str(jk,bnrs2,ik)
    fmtstr='('//bnrs2(1:jj)//'(i'//bnrs(1:ii)//',1x))'
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'_STRUCT_CLUSTERING.clu'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'STRUCT_CLUSTERING.clu'
    end if
#else
    fn = 'STRUCT_CLUSTERING.clu'
#endif
!
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if (exists.EQV..true.) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old',position='append')
      close(unit=iu,status='delete')
    end if
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
  end if
!
  ll = 0
  do j=1,nbasins
    do k=1,it(j)%nmbrs
      if ((it(j)%snaps(k).le.0).OR.(it(j)%snaps(k).gt.cstored)) then
        write(ilog,*) 'Fatal. Encountered nonexistent snapshot in cluster member list. This is a bug.'
        call fexit()
      end if
      sconnect(it(j)%snaps(k),4) = j
    end do
    it(j)%nodewt(1) = (1.0*it(j)%nmbrs)/(1.0*cstored)
    it(j)%nodewt(2:5) = 0.0
  end do
  if ((c_multires.gt.0).AND.((cmode.eq.5).OR.((cmode.eq.4).AND.(cprogindex.eq.2)))) then
    allocate(ivx(cstored,c_multires))
    do ll=c_nhier,c_nhier-c_multires+1,-1
      do j=1,birchtree(ll)%ncls
        do k=1,birchtree(ll)%cls(j)%nmbrs
          if ((birchtree(ll)%cls(j)%snaps(k).le.0).OR.(birchtree(ll)%cls(j)%snaps(k).gt.cstored)) then
            write(ilog,*) 'Fatal. Encountered nonexistent snapshot in cluster member list. This is a bug.'
            call fexit()
          end if
          ivx(birchtree(ll)%cls(j)%snaps(k),c_nhier-ll+1) = j
        end do
      end do
    end do
  else if (c_multires.gt.0) then
    c_multires = 0
  end if
  if ((prtclu.EQV..true.).AND.(cmode.ne.6).AND.(cmode.ne.7)) then
    if (c_multires.gt.0) then
      do i=1,cstored
        write(iu,fmtstr) sconnect(i,4),(ivx(i,ii),ii=1,c_multires)
      end do
    else
      do i=1,cstored
        write(iu,fmtstr) sconnect(i,4)
      end do
    end if
    close(unit=iu)
  end if
  if (allocated(ivx).EQV..true.) deallocate(ivx)
!
  do i=1,nbasins
    if (allocated(it(i)%wghtsnb).EQV..true.) deallocate(it(i)%wghtsnb)
    if (allocated(it(i)%lensnb).EQV..true.) deallocate(it(i)%lensnb)
    if (allocated(it(i)%fewtsnb).EQV..true.) deallocate(it(i)%fewtsnb)
    if (allocated(it(i)%map).EQV..true.) deallocate(it(i)%map)
    if (allocated(it(i)%lstnb).EQV..true.) deallocate(it(i)%lstnb)
    if (allocated(it(i)%flwnb).EQV..true.) deallocate(it(i)%flwnb)
    it(i)%nbalsz = 2
    allocate(it(i)%wghtsnb(it(i)%nbalsz,2))
    allocate(it(i)%lensnb(it(i)%nbalsz,2))
    allocate(it(i)%fewtsnb(it(i)%nbalsz,4))
    it(i)%wghtsnb(:,:) = 0
    it(i)%lensnb(:,:) = 0.0
    it(i)%fewtsnb(:,:) = 0.0
    allocate(it(i)%lstnb(it(i)%nbalsz))
    allocate(it(i)%map(it(i)%nbalsz))
    it(i)%nb = 0
  end do
!
  if ((brklnk_report.ge.1).AND.(cmode.ne.7)) then
 456 format(2(i10,1x),g14.6)
    write(ilog,*)
    write(ilog,'(a)') '--------------- REPORT OF SNAPSHOT LINK DISTANCES ------------------'
    if (brklnk_report.eq.1) then
      write(ilog,'(a,i7,a)') ' Only geometric distances for links with a lag time different from ',clagt_msm,' stored snapshots &
 &are printed.'
    end if
  end if
  do i=1,sconnectsz
    if ((sconnect(i,1).le.0).OR.(sconnect(i,2).le.0)) cycle
    if ((sconnect(i,1).gt.cstored).OR.(sconnect(i,2).gt.cstored)) then
      write(ilog,*) 'Warning. Skipping all graph-related analysis due to mismatch between snapshot connectivity matrix &
 &and data size. This is most likely caused by a truncated trajectory analysis run.'
      sconnectsz = -sconnectsz
      exit
    end if
    ll = sconnect(sconnect(i,1),4)
    ii = sconnect(sconnect(i,2),4)
    if (cmode.ne.7) then
      call snap_to_snap_d(dis,sconnect(i,1),sconnect(i,2))
      if (brklnk_report.ge.1) then
        if ((brklnk_report.eq.2).OR.(sconnect(i,2)-sconnect(i,1).ne.clagt_msm)) then
          write(ilog,456) (sconnect(i,j),j=1,2),dis
        end if
      end if
      dval = dis*dis
    else
      dval = 1.0
    end if
    foundit = .false.
    do j=1,it(ii)%nb
      if (it(ii)%lstnb(j).eq.ll) then
        it(ii)%wghtsnb(j,2) = it(ii)%wghtsnb(j,2) + 1
        it(ii)%lensnb(j,2) = it(ii)%lensnb(j,2) + dval
        it(ii)%fewtsnb(j,2) = it(ii)%fewtsnb(j,2) + 1.0
        foundit = .true.
        ik = j
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      it(ii)%nb = it(ii)%nb + 1
      ik = it(ii)%nb
      if (it(ii)%nb.gt.it(ii)%nbalsz) then
        call scluster_resizenb(it(ii))
      end if
!      it(ii)%map(ll) = it(ii)%nb
      it(ii)%lstnb(it(ii)%nb) = ll
      it(ii)%wghtsnb(it(ii)%nb,2) = 1
      it(ii)%wghtsnb(it(ii)%nb,1) = 0
      it(ii)%lensnb(it(ii)%nb,2) = dval
      it(ii)%lensnb(it(ii)%nb,1) = 0.0
      it(ii)%fewtsnb(it(ii)%nb,2) = 1.0
      it(ii)%fewtsnb(it(ii)%nb,1) = 0.0
    end if
    foundit = .false.
    do j=1,it(ll)%nb
      if (it(ll)%lstnb(j).eq.ii) then
        it(ll)%wghtsnb(j,1) = it(ll)%wghtsnb(j,1) + 1
        it(ll)%lensnb(j,1) = it(ll)%lensnb(j,1) + dval
        it(ll)%fewtsnb(j,1) = it(ll)%fewtsnb(j,1) + 1.0
        foundit = .true.
        jk = j
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      it(ll)%nb = it(ll)%nb + 1
      jk = it(ll)%nb
      if (it(ll)%nb.gt.it(ll)%nbalsz) then
        call scluster_resizenb(it(ll))
      end if
!      it(ll)%map(ii) = it(ll)%nb
      it(ll)%lstnb(it(ll)%nb) = ii
      it(ll)%wghtsnb(it(ll)%nb,1) = 1
      it(ll)%wghtsnb(it(ll)%nb,2) = 0
      it(ll)%lensnb(it(ll)%nb,1) = dval
      it(ll)%lensnb(it(ll)%nb,2) = 0.0
      it(ll)%fewtsnb(it(ll)%nb,1) = 1.0
      it(ll)%fewtsnb(it(ll)%nb,2) = 0.0
    end if
!   we store the index for the list of node 2 (ll ; jk) in node 1 at the position of the neighbor (ii ; ik)
    it(ii)%map(ik) = jk
!   we store the index for the list of node 1 (ii ; ik) in node 2 at the position of the neighbor (ll ; jk)
    it(ll)%map(jk) = ik
  end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  if (sconnectsz.le.0) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    if (tpi.le.1) sconnectsz = -sconnectsz
    return
  end if
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  if ((brklnk_report.ge.1).AND.(cmode.ne.7)) then
    write(ilog,'(a)') '----------- END OF REPORT OF SNAPSHOT LINK DISTANCES ---------------'
    write(ilog,*)
  end if
!
! allocate flow double vector and normalize MSD
  do i=1,nbasins
    if (it(i)%nb.gt.0) then
      allocate(it(i)%flwnb(it(i)%nbalsz,2))
      it(i)%flwnb(:,:) = 0
    end if
    do j=1,it(i)%nb
      if (it(i)%wghtsnb(j,1).gt.0) then
        it(i)%lensnb(j,1) = it(i)%lensnb(j,1)/(1.0*it(i)%wghtsnb(j,1))
      end if
      if (it(i)%wghtsnb(j,2).gt.0) then
        it(i)%lensnb(j,2) = it(i)%lensnb(j,2)/(1.0*it(i)%wghtsnb(j,2))
      end if
    end do
  end do
  k = 1
  call System_Clock(count=tcnt(2))
  write(ilog,'(1x,a,g11.3,a)') 'Time elapsed for graph construction: ',(tcnt(2)-tcnt(1))/thrdat%rate,' [s]'
  write(ilog,*)
!
! if requested, augment network
  call augment_network(it,nbasins,caddlkmd) ! implies calling tarjan_SCC -> adds links and computes strongly connected components
!    call augment_network(it,nbasins,k,thr_svte(:,1)) ! also calls tarjan_SCC
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  if ((tpi.le.1).AND.((cequil.gt.0).OR.(dopfold.EQV..true.).OR.(synmode.ne.0).OR.(eigvalmd.ne.0).OR.(ccfepmode.gt.0))) then
    call System_Clock(count=tcnt(19))
    write(ilog,*)
    write(ilog,'(1x,a)') 'Now processing (editing, analyzing) graph ...'
  end if

! if requested, perform weight extraction for bias removal (reflected in graphml file)
  if (cequil.gt.0) then  
    dval = 1.0e-10
    call bias_removal(it,nbasins,dval,tpi)
  end if
!
  if ((dopfold.EQV..true.).OR.(synmode.ne.0).OR.(eigvalmd.ne.0).OR.(ccfepmode.gt.0)) then
!   set refscc, inissnap, and endssnap to legal values 
    call sel_ssnap(it,nbasins,sconnect,sconnectsz,inissnap,endssnap,tpi)
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  if ((dopfold.EQV..true.).OR.(synmode.ne.0).OR.(eigvalmd.ne.0).OR.(ccfepmode.gt.0)) then
!   set refscc, inissnap, and endssnap to legal values 
    if (tmat_report.EQV..true.) then
      i = 1
      j = 1
      call edit_tmat(it,nbasins,i,j)
      i = refscc
      if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
        j = 3
        tmsfx = 'FWD'
        call prt_tmat(it,nbasins,j,tmsfx,i)
      end if
      if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
        j = 4
        tmsfx = 'BWD'
        call prt_tmat(it,nbasins,j,tmsfx,i)
      end if
    end if
  else
    if (tmat_report.EQV..true.) then
      i = 1
      j = 0
      call edit_tmat(it,nbasins,i,j)
      i = 0
      if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
        j = 3
        tmsfx = 'FWD'
        call prt_tmat(it,nbasins,j,tmsfx,i)
      end if
      if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
        j = 4
        tmsfx = 'BWD'
        call prt_tmat(it,nbasins,j,tmsfx,i)
      end if
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  if ((tpi.le.1).AND.((cequil.gt.0).OR.(dopfold.EQV..true.).OR.(synmode.ne.0).OR.(eigvalmd.ne.0).OR.(ccfepmode.gt.0))) then
    call System_Clock(count=tcnt(20))
    write(ilog,'(1x,a,g11.3,a)') 'Time elapsed for graph processing: ',(tcnt(20)-tcnt(19))/thrdat%rate,' [s]'
    write(ilog,*)
  end if
!
#ifdef LINK_HSL
! Eigenvalues and eigenvectors (timescales and projections) 
  if (eigvalmd.ne.0) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call System_Clock(count=tcnt(3))
    write(ilog,*) 
    write(ilog,'(1x,a)') 'Now computing requested spectral properties from MSM transition matrix(ces) ...'
!   Fix the initial and end snapshots first and reference component
    if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
      dofwr = .true.
      call sparse_eigen(it,nbasins,refscc,dofwr,tpi,inissnap) 
    end if
    if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
      dofwr = .false.
      call sparse_eigen(it,nbasins,refscc,dofwr,tpi,inissnap) 
    end if
    call System_Clock(count=tcnt(4))
    write(ilog,'(1x,a,g11.3,a)') 'Time elapsed for spectral property solvers: ',(tcnt(4)-tcnt(3))/thrdat%rate,' [s]'
    write(ilog,*)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end if
!
  if (dopfold.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call System_Clock(count=tcnt(5))
    write(ilog,*)
    write(ilog,'(1x,a)') 'Now computing requested committor probabilities from MSM transition matrix(ces) ...'
    if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
      dofwr = .true.
      call solve_pfold_lsys(it,nbasins,sconnect,sconnectsz,dofwr,tpi)
    end if
    if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
      dofwr = .false.
      call solve_pfold_lsys(it,nbasins,sconnect,sconnectsz,dofwr,tpi)
    end if
    call System_Clock(count=tcnt(6))
    write(ilog,*)
    write(ilog,'(1x,a,g11.3,a)') 'Time elapsed for committor probability solvers: ',(tcnt(6)-tcnt(5))/thrdat%rate,' [s]'
    write(ilog,*)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end if
#endif
!
! Random walker on network: synthetic trajectories: since it does not use matvals_*wr it's fine even after pfold
  if (synmode.ne.0) then
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
!   Fix the initial and end snapshots first and reference component
    if (synmode.eq.1) then
      if (it(sconnect(inissnap,4))%inscc.ne.it(sconnect(endssnap,4))%inscc) then
        write(ilog,*) 'Warning. Initial and target snapshots for the generation of synthetic trajectories do not belong to the & 
 &same strongly connected component. This implies irreversibility and is not allowed. Disabling.'
        synmode = 0
      else if (maxval(it(1:nbasins)%inscc).ne.minval(it(1:nbasins)%inscc)) then
        write(ilog,*) 'Warning. Multiple strongly connected components were detected. For the current settings of SYNTRAJ_MD, & 
 &the search is limited to that component the initial and target cluster belong to.'
      end if
    else if (maxval(it(1:nbasins)%inscc).ne.minval(it(1:nbasins)%inscc)) then
      write(ilog,*) 'Warning. Multiple strongly connected components detected. For the current settings of SYNTRAJ_MD, & 
 &this implies that the random walker is restricted to the component the initial cluster belongs to.'
    end if
    if (synmode.ne.0) then
!     populate fewtsnb(:,3:4) with the cumulative sums for picking via binary_search
      call edit_tmat(it,nbasins,athree,aone)
    end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    if (synmode.ne.0) then
      if (tpi.le.1) then
        call System_Clock(count=tcnt(7))
        write(ilog,*)
        write(ilog,'(1x,a)') 'Now performing requested random walks (synthetic trajectories -> FMCSC_SYNTRAJ_MD) ...'
      end if
      if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
        dofwr = .true.
        call synth_traj_MSM(it,nbasins,sconnect(1:cstored,4),inissnap,endssnap,dofwr,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      end if
      if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
        dofwr = .false.
        call synth_traj_MSM(it,nbasins,sconnect(1:cstored,4),inissnap,endssnap,dofwr,tpi)
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    if (tpi.le.1) then
      call System_Clock(count=tcnt(8))
      write(ilog,*)
      write(ilog,'(1x,a,g11.3,a)') 'Time elapsed for random walkers: ',(tcnt(8)-tcnt(7))/thrdat%rate,' [s]'
      write(ilog,*)
    end if
  end if
!
  if (ccfepmode.gt.0) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call System_Clock(count=tcnt(9))
    write(ilog,*)
    write(ilog,'(1x,a)') 'Now computing cut-based free energy profiles ...'
    dval = 1.0e-10
    if (size(csccwts).gt.1) then
      if ((inissnap_usrslct.gt.0).OR.((ccfepmode.gt.7).AND.(ccfepmode.le.10))) then
        write(ilog,*) 'Warning. The computed cut-based free energy profiles will only include that strongly connected &
 &component containing the chosen reference snapshot (FMCSC_INISYNSNAP).'
      else if (inissnap_usrslct.eq.0) then
        write(ilog,*) 'Warning. The computed cut-based free energy profiles will only include that strongly connected &
 &component containing the cluster with the largest sampling weight.'
      else
        write(ilog,*) 'Remark. Potentially computing separate cut-based free energy profiles for individual strongly connected &
 &components from the respective clusters with the largest sampling weights.'
      end if
    end if
    do kk=1,size(csccwts) ! loop over components
!     use only target component unless requested otherwise
      if (((inissnap_usrslct.ge.0).OR.((ccfepmode.gt.7).AND.(ccfepmode.le.10))).AND.(refscc.ne.kk)) cycle
      foundit = .false.
      do j=1,2
        if ((which_timeflow.eq.1).AND.(j.eq.2)) cycle
        if ((which_timeflow.eq.2).AND.(j.eq.1)) cycle
!       always use raw sampling weight as reference probability vector (even though this may not truly be the "largest" in ss)
        ssix = 1
        if (j.eq.1) then
          dofwr = .true.
        else
          dofwr = .false.
        end if
!       for the weight in ssix, find the largest cluster for the component in question
        dis = 0.
        ll = 0
        do i=1,nbasins
          if (it(i)%inscc.eq.kk) then
            ll = ll + 1
            if (it(i)%nodewt(ssix).gt.dis) then
              k = i
              dis = it(i)%nodewt(ssix)
            end if
          end if
        end do
        if (ll.gt.1) then
          if (inissnap_usrslct.gt.0) k = sconnect(inissnap,4)
          if ((ccfepmode.le.7).OR.(ccfepmode.gt.10)) then
            call cfep(it,nbasins,k,dval,ccfepmode,dofwr)
          else
            if ((ccfepmode.eq.8).OR.(ccfepmode.eq.10)) then ! (+) committor
              ii = 8
              call cfep(it,nbasins,k,dval,ii,dofwr)
            end if
            if ((ccfepmode.eq.9).OR.(ccfepmode.eq.10)) then ! (-) committor
              ii = 9
              call cfep(it,nbasins,k,dval,ii,dofwr)
            end if
          end if
        else if (foundit.EQV..false.) then
          write(ilog,*) 'Warning. Skipping calculation of cut-based free energy profile(s) for component ',kk,' as it only &
 &comprises a single cluster.'
          foundit = .true.
        end if
      end do
    end do
    call System_Clock(count=tcnt(10))
    write(ilog,*)
    write(ilog,'(1x,a,g11.3,a)') 'Time elapsed for cut-based free energy profile(s): ',(tcnt(10)-tcnt(9))/thrdat%rate,' [s]'
    write(ilog,*)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!--------------------------------------------------------------------------------
!
! get the number and list of unassigned snapshots that space transitions between basins
! WARNING: incomplete
!
subroutine get_graph_unassigned(it,nbasins,ivm)
!
  use clusters
!
  implicit none
!
  integer i,j,k,ll,nbasins,sep
  integer, ALLOCATABLE:: iv1(:)
  integer ivm(nbasins,nbasins)
  logical foundit
  type(t_scluster) it(nbasins)

  allocate(iv1(cstored))
  iv1(:) = 0
  ll = 0
  do i=1,cstored
    foundit = .false.
    do j=1,nbasins
      do k=1,it(j)%nmbrs
        if (i.eq.it(j)%snaps(k)) then
          foundit = .true.
          iv1(i) = j
          exit
        end if
      end do
      if (foundit.EQV..true.) exit
    end do
  end do
!
  k = 1
  do while (iv1(k).eq.0)
    k = k + 1
  end do
  ll = iv1(k)
  sep = 0
  do i=k+1,cstored
    if (iv1(i).eq.0) then
      sep = sep + 1
      cycle
    end if
!
    if ((iv1(i).gt.0).AND.(ll.gt.0)) then
      ivm(ll,iv1(i)) = ivm(ll,iv1(i)) + sep
      sep = 0
    end if
    ll = iv1(i)
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine graphml_helper_for_clustering(it,nclusters)
!
  use iounit
  use mpistuff
  use clusters
  use molecule
  use sequen
  use aminos
  use system
  use mcsums
  use pdb
!
  implicit none
!
  integer i,j,k,ii,jj,ipdbh,freeunit,nclusters,ilen,jlen,klen
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
  character(MAXSTRLEN) istr,jstr,kstr,iistr
  logical exists,havems(4)
  type(t_scluster) it(nclusters)
!

#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  if (use_REMC.EQV..true.) then
    fn = 'N_'//nod(1:tl)//'_STRUCT_CLUSTERING.graphml'
  else
    fn = 'STRUCT_CLUSTERING.graphml'
  end if
#else
  fn = 'STRUCT_CLUSTERING.graphml'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ipdbh = freeunit()
    open(unit=ipdbh,file=fn(ii:jj),status='old',position='append')
    close(unit=ipdbh,status='delete')
  end if
  ipdbh=freeunit()
  open(unit=ipdbh,file=fn(ii:jj),status='new')
  havems(:) = .false.
!
 77  format(a)
 78 format(a,es11.5,a)
  write(ipdbh,77)  '<?xml version="1.0" encoding="UTF-8"?>'
  write(ipdbh,77) '<graphml xmlns="http://graphml.graphdrawing.org/xmlns"'
  write(ipdbh,*) '    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
  write(ipdbh,*) '    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">'
  write(ipdbh,*) ' <key id="d0" for="node" attr.name="nsampwt" attr.type="double"/>'
  if (cmode.ne.7) write(ipdbh,*) ' <key id="d1" for="node" attr.name="nradius" attr.type="double"/>'
  if (cmode.ne.7) write(ipdbh,*) ' <key id="d2" for="node" attr.name="ndiamet" attr.type="double"/>'
  if (cmode.ne.7) write(ipdbh,*) ' <key id="d3" for="node" attr.name="nqualit" attr.type="double"/>'
  if (sum(it(1:nclusters)%nodewt(2)).gt.0.0) then
    havems(1) = .true.
    write(ipdbh,*) ' <key id="d6" for="node" attr.name="nitfwwt" attr.type="double"/>'
  end if
  if (sum(it(1:nclusters)%nodewt(3)).gt.0.0) then
    havems(2) = .true.
    write(ipdbh,*) ' <key id="d7" for="node" attr.name="nitbwwt" attr.type="double"/>'
  end if
  if (sum(it(1:nclusters)%nodewt(4)).gt.0.0) then
    havems(3) = .true.
    write(ipdbh,*) ' <key id="d8" for="node" attr.name="nlafwwt" attr.type="double"/>'
  end if
  if (sum(it(1:nclusters)%nodewt(5)).gt.0.0) then
    havems(4) = .true.
    write(ipdbh,*) ' <key id="d9" for="node" attr.name="nlabwwt" attr.type="double"/>'
  end if
  write(ipdbh,*) ' <key id="i0" for="node" attr.name="ncenter" attr.type="int"/>'
  write(ipdbh,*) ' <key id="i2" for="node" attr.name="nsccidx" attr.type="int"/>'
  write(ipdbh,*) ' <key id="i1" for="edge" attr.name="eweight" attr.type="int"/>'
  if (cmode.ne.7) write(ipdbh,*) ' <key id="d4" for="edge" attr.name="egeomsd" attr.type="double"/>'
  write(ipdbh,*) ' <key id="d5" for="edge" attr.name="efloatw" attr.type="double"/>'
  write(ipdbh,*)
  write(ipdbh,*) ' <graph id="G" edgedefault="directed">'
!
  do i=1,nclusters
    ilen = log10(1.0*i) + 1
    write(istr,*) ilen
! 33 format(a,i<ilen>,a)
    write(ipdbh,'("    <node id=""n",i' // ADJUSTL(istr) // ',""">")') i
    write(ipdbh,78) '      <data key="d0">',it(i)%nodewt(1),'</data>'
    if (cmode.ne.7) write(ipdbh,78) '      <data key="d1">',it(i)%radius,'</data>'
    if (cmode.ne.7) write(ipdbh,78) '      <data key="d2">',it(i)%diam,'</data>'
    if (cmode.ne.7) write(ipdbh,78) '      <data key="d3">',max(0.0,it(i)%quality),'</data>'
    if (havems(1).EQV..true.) write(ipdbh,78) '      <data key="d6">',it(i)%nodewt(2),'</data>'
    if (havems(2).EQV..true.) write(ipdbh,78) '      <data key="d7">',it(i)%nodewt(3),'</data>'
    if (havems(3).EQV..true.) write(ipdbh,78) '      <data key="d8">',it(i)%nodewt(4),'</data>'
    if (havems(4).EQV..true.) write(ipdbh,78) '      <data key="d9">',it(i)%nodewt(5),'</data>'
    write(iistr,*) int(log10(1.0*max(it(i)%center,1)) + 1)
    write(ipdbh,'("      <data key=""i0"">",i' // ADJUSTL(iistr) // ',"</data>")') it(i)%center
    write(iistr,*) int(log10(1.0*max(it(i)%inscc,1)) + 1)
    write(ipdbh,'("      <data key=""i2"">",i' // ADJUSTL(iistr) // ',"</data>")') it(i)%inscc
    write(ipdbh,*) '    </node>'
  end do
  k = 1
  do i=1,nclusters
    ilen = log10(1.0*i) + 1
    write(istr,*) ilen
    do j=1,it(i)%nb
      if (it(i)%wghtsnb(j,1).le.0) cycle
      if (k.le.0) then
        klen = 1
      else
        klen = log10(1.0*k) + 1
      end if
      jlen = log10(1.0*max(it(i)%lstnb(j),1)) + 1
      write(jstr,*) jlen
      write(kstr,*) klen
      write(iistr,*) int(log10(1.0*max(it(i)%wghtsnb(j,1),1)) + 1)
! 34   format(a,i<klen>,a,i<ilen>,a,i<jlen>,a)
!      write(ipdbh,34) '     <edge id="e',k,'" source="n',i,'" target="n',it(i)%lstnb(j),'">'
      write(ipdbh,'("     <edge id=""e",i' // ADJUSTL(kstr) // ',""" source=""n",i' // ADJUSTL(istr) // &
 &'""" target=""n",i' // ADJUSTL(jstr) // '""">")') k,i,it(i)%lstnb(j)
      write(ipdbh,'("      <data key=""i1"">",i' // ADJUSTL(iistr) // ',"</data>")') it(i)%wghtsnb(j,1)
      if (cmode.ne.7) write(ipdbh,78) '      <data key="d4">',it(i)%lensnb(j,1),'</data>'
      write(ipdbh,78) '      <data key="d5">',it(i)%fewtsnb(j,1),'</data>'
!      write(ipdbh,*) '      <data key="i1">',it(i)%wghtsnb(j,1),'</data>'
      write(ipdbh,*) '    </edge>'
      k = k + 1
    end do 
  end do
  write(ipdbh,*) '  </graph>'
  write(ipdbh,77) '</graphml>'
  close(unit=ipdbh)
!
end
!
!----------------------------------------------------------------------------------------------------
!
! a subroutine to add links in controlled fashion
!
subroutine augment_network(it,nbasins,mode)
!
  use clusters
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: nbasins,mode
  type(t_scluster):: it(nbasins)
!
  integer i,j,k,ncycles,nedges,nsccs
  logical haveself
!
  k = 1
  call tarjan_SCC(it,nbasins,k)
  nsccs = maxval(it(1:nbasins)%inscc)
!
  if (mode.le.0) return
!
  write(ilog,*) 'Now augmenting network ...'
!
  ncycles = 0 
  nedges = 0
! fix missing self transitions and reverse edges
  do i=1,nbasins
    haveself = .false.
    do j=1,it(i)%nb
      if (it(i)%lstnb(j).eq.i) then
        haveself = .true.
      end if
      if ((it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc).AND.((mode.eq.1).OR.(mode.eq.3))) then
        if (it(i)%wghtsnb(j,1).le.0) then
          nedges = nedges + 1
          it(i)%wghtsnb(j,1) = 1
          it(i)%fewtsnb(j,1) = caddlinkwt
          it(i)%lensnb(j,1) = it(it(i)%lstnb(j))%lensnb(it(i)%map(j),1) ! use same value as existing edge
          it(it(i)%lstnb(j))%wghtsnb(it(i)%map(j),2) = 1
          it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),2) = caddlinkwt
        else if (it(i)%wghtsnb(j,2).le.0) then
          nedges = nedges + 1
          it(i)%wghtsnb(j,2) = 1
          it(i)%fewtsnb(j,2) = caddlinkwt
          it(i)%lensnb(j,2) = it(it(i)%lstnb(j))%lensnb(it(i)%map(j),2) ! use same value as existing edge
          it(it(i)%lstnb(j))%wghtsnb(it(i)%map(j),1) = 1
          it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),1) = caddlinkwt
        end if
      end if
      if ((mode.eq.4).OR.(mode.eq.5)) then
        if (it(i)%wghtsnb(j,2).gt.it(i)%wghtsnb(j,1)) then
          k = it(i)%wghtsnb(j,2) - it(i)%wghtsnb(j,1)
          nedges = nedges + k
          it(i)%lensnb(j,1) = (it(i)%wghtsnb(j,1)*it(i)%lensnb(j,1) + k*it(i)%lensnb(j,2))/(1.0*it(i)%wghtsnb(j,2))
          it(i)%fewtsnb(j,1) = it(i)%fewtsnb(j,1) + 1.0*k
          it(i)%wghtsnb(j,1) = it(i)%wghtsnb(j,1) + k
          it(it(i)%lstnb(j))%lensnb(it(i)%map(j),2) = it(i)%lensnb(j,1)
          it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),2) = it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),2) + 1.0*k
          it(it(i)%lstnb(j))%wghtsnb(it(i)%map(j),2) = it(it(i)%lstnb(j))%wghtsnb(it(i)%map(j),2) + k
        else if (it(i)%wghtsnb(j,1).gt.it(i)%wghtsnb(j,2)) then
          k = it(i)%wghtsnb(j,1) - it(i)%wghtsnb(j,2)
          nedges = nedges + k
          it(i)%lensnb(j,2) = (it(i)%wghtsnb(j,2)*it(i)%lensnb(j,2) + k*it(i)%lensnb(j,1))/(1.0*it(i)%wghtsnb(j,1))
          it(i)%fewtsnb(j,2) = it(i)%fewtsnb(j,2) + 1.0*k
          it(i)%wghtsnb(j,2) = it(i)%wghtsnb(j,2) + k
          it(it(i)%lstnb(j))%lensnb(it(i)%map(j),1) = it(i)%lensnb(j,2)
          it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),1) = it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),1) + 1.0*k
          it(it(i)%lstnb(j))%wghtsnb(it(i)%map(j),1) = it(it(i)%lstnb(j))%wghtsnb(it(i)%map(j),1) + k
        end if
      end if
    end do
    if ((haveself.EQV..false.).AND.((mode.eq.2).OR.(mode.eq.3).OR.(mode.eq.5))) then
      ncycles = ncycles + 1
      it(i)%nb = it(i)%nb + 1
      if (it(i)%nb.gt.it(i)%nbalsz) then
        call scluster_resizenb(it(i))
      end if
      it(i)%map(it(i)%nb) = it(i)%nb
      it(i)%lstnb(it(i)%nb) = i
      it(i)%wghtsnb(it(i)%nb,1:2) = 1
      it(i)%lensnb(it(i)%nb,1:2) = (0.5*cradius)**2
      it(i)%fewtsnb(it(i)%nb,1:2) = caddlinkwt
    end if
  end do
!
! recompute SCCs
  k = 2
  call tarjan_SCC(it,nbasins,k)
  k = maxval(it(1:nbasins)%inscc)
 560 format('Done. Augmentation used a threshold distance of ',g14.5,' (various units) based on the quantile of ',g12.5,'.')
 561 format('Added ',i7,' self transitions and ',i7,' missing reverse links for connectnedness.',//,&
 &        'The recalculated number of strongly connected components is ',i7,' (from ',i7,').')
 562 format('Imposition of detailed balance added a total link weight of ',i7,', and ',i7,' self transitions were added.',//,&
 &        'The recalculated number of strongly connected components is ',i7,' (from ',i7,').')
  if (mode.le.3) then
    write(ilog,561) ncycles,nedges,k,nsccs
  else
    write(ilog,562) nedges,ncycles,k,nsccs
  end if
  nsccs = k
!
end

!
!-------------------------------------------------------------------------
!
! snode: reference node (used in 2-7 and always to identify ref. component)
! dofwr: direction of time in transition matrix
! nmode: 1 MFPT; 2-7: exp.; 8-11: committors
!
subroutine cfep(it,nbasins,snode,thresher,nmode,dofwr)
!
  use clusters
  use iounit
  use mpistuff
  use interfaces
  use system
!
  implicit none
!
  logical, INTENT(IN):: dofwr
  integer nbasins,snode,i,j,k,ipdbh,ii,jj,freeunit,tx,nmode,t1,t2,ssix,fbix
  RTYPE thresher,totew,progger
  type(t_scluster) it(nbasins)
  RTYPE, ALLOCATABLE:: mfptv(:),incr(:),hlp1(:),rcutter(:)!,mfpth(:),mfptmsd(:)
  integer, ALLOCATABLE:: cutter(:),iv3(:),iv4(:)
  logical upornot,exists,reqreal,ldum
  character(8) nor
  character(3) tstr
  character(MAXSTRLEN) fn,bitt
#ifdef ENABLE_MPI
  character(3) nod
  integer tl
#endif
!
  if (dofwr.EQV..true.) then
    tstr = 'FWD'
    fbix = 1
    ssix = 0
    do k=4,0,-2
      do i=1,nbasins
        if (it(i)%inscc.ne.it(snode)%inscc) cycle
        if (it(i)%nodewt(max(1,k)).gt.0.) then
          ssix = max(1,k)
          exit
        end if
      end do
      if (ssix.gt.0) exit
    end do
  else
    tstr = 'BWD'
    fbix = 2
    ssix = 0
    do k=5,1,-2
      do i=1,nbasins
        if (it(i)%inscc.ne.it(snode)%inscc) cycle
        if (it(i)%nodewt(k).gt.0.) then
          ssix = k
          exit
        end if
      end do
      if (ssix.gt.0) exit
    end do
  end if
!
! for the cut function, we use the raw link weights, which generally are integer, but may be not 
  reqreal = .false.
  if ((caddlinkwt.ne.1.0).AND.(caddlkmd.ne.0).AND.(caddlkmd.ne.4)) reqreal = .true.
  if (cequil.eq.2) reqreal = .true. ! scaled links
!
! first determine or copy the ordering parameter
  allocate(mfptv(nbasins))
! 
  if ((nmode.le.2).OR.(nmode.ge.11)) then
    i = 1
    j = 1
    call edit_tmat(it,nbasins,i,j)
  end if
  if (nmode.eq.1) then
    call iterative_mfpt(it,nbasins,snode,mfptv,thresher,dofwr,maxtime_eqmsm)
  else if (nmode.eq.2) then ! redundant mode right now
    reqreal = .true.
    call iterative_mfpt(it,nbasins,snode,mfptv,thresher,dofwr,maxtime_eqmsm) 
  else if (nmode.le.7) then
    k = nmode - 2
    call dijkstra_sp(it,nbasins,snode,mfptv,k)
  else if ((nmode.eq.8).AND.(dofwr.EQV..true.)) then
    if ((dopfold.EQV..false.).OR.(maxval(pfolds(:,1)).le.0)) then
      write(ilog,*) 'Warning. Aborting generation of pfold forward time (+) cut profile because the required committor &
 &probabilities are missing or corrupt.'
      return
    end if
    do i=1,nbasins
      mfptv(i) = pfolds(i,1)  !,1 is pfold fwr +
    end do
  else if ((nmode.eq.8).AND.(dofwr.EQV..false.)) then
    if ((dopfold.EQV..false.).OR.(maxval(pfolds(:,2)).le.0)) then
      write(ilog,*) 'Warning. Aborting generation of pfold backward time (+) cut profile because the required committor &
 &probabilities are missing or corrupt.'
      return
    end if
    do i=1,nbasins
      mfptv(i) = pfolds(i,2)  !,2 is pfold bwr +
    end do
  else if ((nmode.eq.9).AND.(dofwr.EQV..true.)) then
    if ((dopfoldminus.EQV..false.).OR.(maxval(pfolds(:,3)).le.0)) then
      write(ilog,*) 'Warning. Aborting generation of pfold forward time (-) cut profile because the required committor &
 &probabilities are missing or corrupt.'
      return
    end if
    do i=1,nbasins
      mfptv(i) = pfolds(i,3)  !,3 is pfold fwr -
    end do
  else if ((nmode.eq.9).AND.(dofwr.EQV..false.)) then
    if ((dopfoldminus.EQV..false.).OR.(maxval(pfolds(:,4)).le.0)) then
      write(ilog,*) 'Warning. Aborting generation of pfold backward time (-) cut profile because the required committor &
 &probabilities are missing or corrupt.'
      return
    end if
    do i=1,nbasins
      mfptv(i) = pfolds(i,4)  !,4 is pfold bwr -
    end do
  else if (nmode.eq.11) then
    ldum = .false.
    call network_order_by_flux(it,nbasins,snode,mfptv,dofwr,ldum)
  else if (nmode.eq.12) then
    ldum = .true.
    call network_order_by_flux(it,nbasins,snode,mfptv,dofwr,ldum)
  end if
!
  totew = 0.0
  do j=1,nbasins
    if (it(j)%nmbrs.gt.0) then
      totew = totew + sum(it(j)%fewtsnb(1:it(j)%nb,fbix))
    end if
  end do
!
  tx = 8
  call int2str(snode,nor,tx)
  bitt(:) = ' '
  bitt(1:5) = 'MFPT '
  if (nmode.gt.1) then
    if (nmode.eq.2) bitt(1:10) = 'DIFFUSION '
    if ((nmode.gt.2).AND.(nmode.le.7)) bitt(1:10) = 'SHORTPATH '
    if (nmode.eq.8) bitt(1:11) = 'PFOLD_PLUS '
    if (nmode.eq.9) bitt(1:12) = 'PFOLD_MINUS '
    if (nmode.eq.11) bitt(1:13) = 'FLUX_PIX_OUT '
    if (nmode.eq.12) bitt(1:12) = 'FLUX_PIX_IN '
  end if
  call strlims(bitt,t1,t2)
#ifdef ENABLE_MPI
  tl = 3
  call int2str(myrank,nod,tl)
  if (use_REMC.EQV..true.) then
    fn = 'N_'//nod(1:tl)//'_'//bitt(t1:t2)//'_'//tstr//'_CFEP_'//nor(1:tx)//'.dat'
  else
    fn = bitt(t1:t2)//'_'//tstr//'_CFEP_'//nor(1:tx)//'.dat'
  end if
#else
  fn = bitt(t1:t2)//'_'//tstr//'_CFEP_'//nor(1:tx)//'.dat'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ipdbh = freeunit()
    open(unit=ipdbh,file=fn(ii:jj),status='old',position='append')
    close(unit=ipdbh,status='delete')
  end if
  ipdbh=freeunit()
  open(unit=ipdbh,file=fn(ii:jj),status='new')
!
  allocate(incr(nbasins))
  allocate(iv3(nbasins))
  allocate(iv4(nbasins))
  allocate(hlp1(nbasins))
  incr(:) = mfptv(:)
  do i=1,nbasins
    iv3(i) = i
  end do
  ii = 1
  jj = nbasins
  k = nbasins
  if (nmode.eq.8) then 
    upornot = .false.
    call merge_sort(k,upornot,incr,hlp1,ii,jj,iv3,iv4)
  else  !pfep (+)
    upornot = .true.
    call merge_sort(k,upornot,incr,hlp1,ii,jj,iv3,iv4)
  end if
!
  deallocate(incr)
  deallocate(iv3)
!
!  allocate(mfpth(100))
  if (reqreal.EQV..true.) then
    allocate(rcutter(nbasins))
    rcutter(:) = 0.0
  end if
  allocate(cutter(nbasins))
  cutter(:) = 0
!  allocate(mfptmsd(nbasins))
!
! for all edges 
  hlp1(:) = 0.0
!  mfpth(:) = 0.0
!  mfptmsd(:) = 0.0
!  binsz = (maxval(mfptv(:))-minval(mfptv(:)))/100.
!  mfptmsd(1) = 1.0
  do i=1,nbasins
    if (it(iv4(i))%inscc.ne.it(snode)%inscc) cycle
    do j=1,it(iv4(i))%nb
      jj = it(iv4(i))%lstnb(j) ! loop over all edges
      if (it(jj)%inscc.ne.it(snode)%inscc) cycle
      if (nmode.eq.8) then ! (+) committor
        do k=i,nbasins ! all clusters with smaller values of (+) committor prob. than the source of the current edge, iv4(i)
          if (mfptv(jj).le.mfptv(iv4(k))) then
!           the total cut weight as the sum over all snap-to-snap transitions
            if (reqreal.EQV..true.) then
              rcutter(k) = rcutter(k) + it(iv4(i))%fewtsnb(j,fbix)
            end if
            cutter(k) = cutter(k) + it(iv4(i))%wghtsnb(j,fbix)
!           the mean square distance across the cut as a linear average over all snap-to-snap transitions 
            hlp1(k) = hlp1(k) + it(iv4(i))%lensnb(j,fbix)*it(iv4(i))%wghtsnb(j,fbix)
!           the MSD in the variable used for ordering across the cut
!           mfptmsd(k) = mfptmsd(k) + it(iv4(i))%wghtsnb(j,fbix)*(mfptv(jj)-mfptv(iv4(k)))**2
          else
            exit ! condition cannot hold for any further cluster due to sorting -> move to next edge
          end if
        end do
      else
        do k=i,nbasins ! all clusters with larger order parameter value than the source of the current edge, iv4(i) 
          if (mfptv(jj).ge.mfptv(iv4(k))) then
!           the total cut weight as the sum over all snap-to-snap transitions
            if (reqreal.EQV..true.) then
              rcutter(k) = rcutter(k) + it(iv4(i))%fewtsnb(j,fbix)
            end if
            cutter(k) = cutter(k) + it(iv4(i))%wghtsnb(j,fbix)
!           the mean square distance across the cut as a linear average over all snap-to-snap transitions 
            hlp1(k) = hlp1(k) + it(iv4(i))%lensnb(j,fbix)*it(iv4(i))%wghtsnb(j,fbix)
!           the MSD in the variable used for ordering across the cut
!           mfptmsd(k) = mfptmsd(k) + it(iv4(i))%wghtsnb(j,fbix)*(mfptv(jj)-mfptv(iv4(k)))**2
          else
            exit ! condition cannot hold for any further cluster due to sorting -> move to next edge
          end if
        end do
      end if
    end do
!    k = max(1,min(100,ceiling((mfptv(iv4(i))-minval(mfptv(:)))/binsz)))
!    mfpth(k) = mfpth(k) + 1.0*it(iv4(i))%nmbrs
  end do
!
 14 format(i8,1x,i8,1x,g15.8,1x,i8,1x,g11.4,1x,g13.6,1x,g13.6,1x,g13.6,1x,g13.6)
 15 format(i8,1x,i8,1x,g15.8,1x,g13.6,1x,g11.4,1x,g13.6,1x,g13.6,1x,g13.6,1x,g13.6)
  jj = 0
  progger = 0.
!
  do i=1,nbasins
    if (it(iv4(i))%inscc.ne.it(snode)%inscc) cycle
    jj = jj + 1
    if (reqreal.EQV..true.) then
      rcutter(i) = max(min(caddlinkwt,1.0),rcutter(i)) ! we may have missed a single transition in opposite direction
      write(ipdbh,15) jj,it(iv4(i))%center,progger,rcutter(i),-(1.0/invtemp)*log(rcutter(i)/totew),mfptv(iv4(i)),&
 &                    hlp1(i)/(2.0*cstorecalc*cutter(i))
    else
      cutter(i) = max(1,cutter(i)) ! we may have missed a single transition in opposite direction
      write(ipdbh,14) jj,it(iv4(i))%center,progger,cutter(i),-(1.0/invtemp)*log(cutter(i)/totew),mfptv(iv4(i)),&
 &                    hlp1(i)/(2.0*cstorecalc*cutter(i))!,&
! &                  -(1.0/invtemp)*log(mfptmsd(i)/normerr),mfptmsd(i)/(2.0*cstorecalc*cutter(i))
!  do i=1,100
    end if
!    write(ipdbh,14) i,101-i,minval(mfptv(:))+(i-0.5)*binsz,int(mfpth(i)),-(1.0/invtemp)*log((mfpth(i)+1.)/(1.0*cstored))
    progger = progger + it(iv4(i))%nodewt(ssix)
  end do
!
  close(unit=ipdbh)
!
  if (allocated(cutter).EQV..true.) deallocate(cutter)
  if (allocated(rcutter).EQV..true.) deallocate(rcutter)
  deallocate(iv4)
  deallocate(hlp1)
  deallocate(mfptv)
!  deallocate(mfptmsd)
!  deallocate(mfpth)
!
end
!
!--------------------------------------------------------------------------
!
! WARNING: incomplete
!
subroutine reweight_edges(it,nbasins)
!
  use clusters
  use iounit
!
  implicit none
!
  integer nbasins,i,j
  type(t_scluster) it(nbasins)
!
  do i=1,nbasins
    do j=1,it(i)%nb
      if (it(i)%lensnb(j,1).gt.0.0) it(i)%fewtsnb(j,1) = it(i)%fewtsnb(j,1)*it(i)%lensnb(j,1)/(it(i)%wghtsnb(j,1)*cradius)
      if (it(i)%lensnb(j,2).gt.0.0) it(i)%fewtsnb(j,2) = it(i)%fewtsnb(j,2)*it(i)%lensnb(j,2)/(it(i)%wghtsnb(j,2)*cradius)
    end do
  end do
!
end
!
!-------------------------------------------------------------------------------------
!
! assuming some changes have been made to the transition matrix, equilibrate a MSM to convergence (thresh)
!
subroutine bias_removal(it,nbasins,thresh,tpi)
!
  use clusters
  use iounit
  use system, ONLY: basename,bleng
#ifdef ENABLE_MPI
  use mpistuff, ONLY: myrank
#endif
  use threads, ONLY: thrdat,thr_hlper,thr_rhlper
!
  implicit none
!
  integer, INTENT(IN):: nbasins,tpi
  RTYPE, INTENT(IN):: thresh
!
  integer iterations,i,j,k,nsccs,ii,cbnds(2)!,nfails,ncycles,didtb
  type(t_scluster) it(nbasins)
  RTYPE maxerr,selffac!,normerr,normeo,fraccer,lasterr,normer
!  RTYPE, ALLOCATABLE:: wvec(:),pbase(:,:)
!  logical, ALLOCATABLE:: hit(:)
  logical staterr,augment_network,write_weights,dofwr
  character(MAXSTRLEN) intfile
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
#ifdef ENABLE_THREADS
  integer OMP_GET_NUM_THREADS
!
  if (tpi.gt.0) then
    i = 1
    j = nbasins
    k = OMP_GET_NUM_THREADS()
    call threads_bounds(i,j,tpi,k,cbnds)
  else
#else
  integer freeunit
#endif
!
  cbnds(1) = 1
  cbnds(2) = nbasins  
#ifdef ENABLE_THREADS
  end if
#endif
!
  nsccs = size(csccwts)
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  if ((cequil.eq.2).OR.(cequil.eq.4)) then
    call reweight_edges(it,nbasins)
  end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
! set flag to identify whether network was augmented (which implies that it(:)%fewtsnb is populated)
  if (caddlinks.gt.0.0) then
    augment_network = .true.
  else
    augment_network = .false.
  end if
! whether to directly write frame weights
  write_weights = .true.
!
! simple MSM reweighting
  if ((cequil.eq.1).OR.(cequil.eq.2)) then
    k = 2
    if ((cequil.eq.2).OR.(augment_network.EQV..true.)) k = 3
    if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
      allocate(thr_hlper(nbasins,1))
      allocate(thr_rhlper(nbasins,5))
      write(ilog,*) 'Now reequilibrating forward time network ...'
#ifdef ENABLE_THREADS
!      call System_Clock(count=thr_timings(121,1))
!$OMP END MASTER
!$OMP BARRIER
#endif
      thr_rhlper(cbnds(1):cbnds(2),:) = 0.0
      if (k.eq.3) then
        thr_rhlper(cbnds(1):cbnds(2),1) = it(cbnds(1):cbnds(2))%nodewt(1)
      end if
      selffac = thresh
      j = 2 ! rely on fewtsnb
      dofwr = .true.
      call equilibrate_MSM(it,nbasins,thr_rhlper(:,1),j,k,selffac,iterations,dofwr,maxtime_eqmsm,tpi)
      it(cbnds(1):cbnds(2))%nodewt(2) = thr_rhlper(cbnds(1):cbnds(2),2)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
!      call System_Clock(count=thr_timings(122,1))
!      write(*,*) 'Time for equilibrate_MSM: ',(thr_timings(122,1)-thr_timings(121,1))/thrdat%rate
#endif
      write(ilog,*) 'Done after ',iterations,' iterations.'
      if (which_timeflow.eq.1) then
        deallocate(thr_rhlper)
        deallocate(thr_hlper)
      end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    end if
    if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
      if (which_timeflow.eq.2) then
        allocate(thr_hlper(nbasins,1))
        allocate(thr_rhlper(nbasins,5))
      end if
      write(ilog,*) 'Now reequilibrating backward time network ...'
#ifdef ENABLE_THREADS
!      call System_Clock(count=thr_timings(121,1))
!$OMP END MASTER
!$OMP BARRIER
#endif
      thr_rhlper(cbnds(1):cbnds(2),:) = 0.0
      if (k.eq.3) then
        thr_rhlper(cbnds(1):cbnds(2),1) = it(cbnds(1):cbnds(2))%nodewt(1)
      end if
      selffac = thresh
      j = 2 ! rely on fewtsnb
      dofwr = .false.
      call equilibrate_MSM(it,nbasins,thr_rhlper(:,1),j,k,selffac,iterations,dofwr,maxtime_eqmsm,tpi)
      it(cbnds(1):cbnds(2))%nodewt(3) = thr_rhlper(cbnds(1):cbnds(2),2)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
!      call System_Clock(count=thr_timings(122,1))
!      write(*,*) 'Time for equilibrate_MSM: ',(thr_timings(122,1)-thr_timings(121,1))/thrdat%rate
#endif
      write(ilog,*) 'Done after ',iterations,' iterations.'
      deallocate(thr_rhlper)
      deallocate(thr_hlper)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    end if
!
!  else if ((cequil.eq.3).OR.(cequil.eq.4)) then
!    allocate(pbase(nsccs,2))
!    pbase(:,:) = 0.0
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!    allocate(thr_hlper(cstored,2))
!    allocate(thr_rhlper(nbasins,6))
!    write(ilog,*) 'Now reequilibrating network ...'
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    j = 2
!    k = 3
!    call equilibrate_MSM(it,nbasins,thr_rhlper(:,1),j,k,thresh,iterations,dofwr,maxtime_eqmsm,tpi)
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!    thr_rhlper(1:nsccs,3:4) = 0.0
!    it(1:nbasins)%nodewt = thr_rhlper(:,2)
!    write(ilog,*) 'Done after ',iterations,' iterations.'
!!   guess the underlying flat propagator
!    do i=1,nbasins
!      thr_rhlper(it(i)%inscc,3) = thr_rhlper(it(i)%inscc,3) + 1.0
!      thr_rhlper(it(i)%inscc,4) = thr_rhlper(it(i)%inscc,4) + it(i)%nodewt
!      normerr = 0.0
!      nfails = 0
!      do j=1,it(i)%nb
!        if (it(i)%fewtsnb(j,1).gt.0.0) then
!          it(i)%fewtsnb(j,1) = it(i)%fewtsnb(j,1)*(it(it(i)%lstnb(j))%nodewt + it(i)%nodewt)/it(it(i)%lstnb(j))%nodewt
!          if (it(i)%lstnb(j).ne.i) normerr = normerr + it(i)%fewtsnb(j,1)
!          if (it(i)%lstnb(j).eq.i) k = j
!          nfails = nfails + it(i)%wghtsnb(j,1)
!        end if
!      end do
!      if (nfails.gt.it(i)%wghtsnb(k,1)) then
!!        it(i)%fewtsnb(k,1) = normerr*it(i)%wghtsnb(k,1)/(1.0*(nfails - it(i)%wghtsnb(k,1)))
!      end if
!    end do
!
!    thr_rhlper(1:nsccs,4) = thr_rhlper(1:nsccs,4)/thr_rhlper(1:nsccs,3)
!    do i=1,nbasins
!      thr_rhlper(i,1) = thr_rhlper(it(i)%inscc,4)
!    end do
!!    thr_rhlper(:,1) = 1.0/(1.0*nbasins)
!    thr_hlper(:,2) = 0
!    do j=1,nbasins
!      do k=1,it(j)%nmbrs
!        thr_hlper(it(j)%snaps(k),2) = j
!      end do
!    end do
!    write(ilog,*) 'Now iterating matrix to give uniform steady state ...'
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    pbase(1:nsccs,1:2) = thr_rhlper(1:nsccs,3:4)
!!   now iterate this matrix to yield a uniform steady state
! 563 format('Total error is ',g20.10,' for steady state obtained by ',i12,' iterations in cycle # ',i5,'.')
! 564 format('Min/max ratios to perfect distribution are ',g16.6,' and ',g16.6,', respectively.')
! 565 format('Self-transition probabilities are ',g16.6,' vs. ',g16.6,' on average (raw vs. uniform).')
!    fraccer = 1.0
!    nfails = 1000
!    do ncycles=1,nfails
!      k = 6
!      j = 2
!      if (ncycles.le.nfails/10.0) then
!        fraccer = thresh
!      else
!        fraccer = max(1.0e-14,thresh-(ncycles-nfails/10.0)*(thresh-1.0e-14)/(0.8*nfails))
!      end if
!      call equilibrate_MSM(it,nbasins,thr_rhlper(:,1),j,k,fraccer,iterations,dofwr,maxtime_eqmsm,tpi)
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!      thr_rhlper(1,5) = 0.0
!      do i=1,nbasins
!        thr_rhlper(1,5) = thr_rhlper(1,5) + abs(thr_rhlper(i,2) - pbase(it(i)%inscc,2))
!        do j=1,it(i)%nb
!          if ((it(i)%fewtsnb(j,1).gt.0.0)) then
!            it(i)%fewtsnb(j,1) = it(i)%fewtsnb(j,1)*0.5*(thr_rhlper(it(i)%lstnb(j),2) + &
! &                            thr_rhlper(i,2))/(thr_rhlper(it(i)%lstnb(j),2))
!          end if
!        end do
!      end do
!      do i=1,nbasins
!        thr_rhlper(i,1) =  pbase(it(i)%inscc,2)
!      end do
!      thr_rhlper(:,1) = thr_rhlper(:,1)/sum(thr_rhlper(:,1))
!!      write(ilog,563) thr_rhlper(1,5),iterations,ncycles
!!      write(ilog,564) minval(thr_rhlper(:,2))*nbasins,maxval(thr_rhlper(:,2))*nbasins
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    end do
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!!   the resultant matrix may not map the kinetics well: reset by adding increment to self-transition
!!   (this does not change the steady state)
!    thr_rhlper(1,4) = 0.0
!    thr_rhlper(1,5) = 0.0
!    fraccer = 0.0
!    do i=1,nbasins
!      do j=1,it(i)%nb
!        if (it(i)%lstnb(j).eq.i) then
!          thr_rhlper(1,4) = thr_rhlper(1,4) + (1.0*it(i)%wghtsnb(j,1)) ! /(1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))
!          thr_rhlper(1,5) = thr_rhlper(1,5) + it(i)%fewtsnb(j,1)
!        end if
!      end do
!      fraccer = fraccer + sum(it(i)%fewtsnb(1:it(i)%nb,1))
!    end do
!    thr_rhlper(1,4) = thr_rhlper(1,4)/(1.0*cstored)
!    thr_rhlper(1,5) = thr_rhlper(1,5)/fraccer
!    write(ilog,565) thr_rhlper(1,4),thr_rhlper(1,5)
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    if (thr_rhlper(1,5).lt.thr_rhlper(1,4)) then
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!      do i=1,nbasins
!        do j=1,it(i)%nb
!          if (it(i)%lstnb(j).eq.i) then
!            it(i)%fewtsnb(j,1) =  it(i)%fewtsnb(j,1) + (thr_rhlper(1,4)-thr_rhlper(1,5))*sum(it(i)%fewtsnb(1:it(i)%nb,1))
!          end if
!        end do
!      end do
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!      j = 2
!      k = 6
!      fraccer = 1.0e-14
!      call equilibrate_MSM(it,nbasins,thr_rhlper(:,1),j,k,fraccer,iterations,dofwr,maxtime_eqmsm,tpi)
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!      normeo = 0.0
!      do i=1,nbasins
!        normeo = normeo + abs(thr_rhlper(i,2) - pbase(it(i)%inscc,2))
!      end do
!      write(ilog,563) normeo,iterations,ncycles
!!      write(ilog,564) minval(thr_rhlper(:,2))*nbasins,maxval(thr_rhlper(:,2))*nbasins
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    end if
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!!    write(*,*) thr_rhlper(:,2)
!    write(ilog,*) 'Done after ',ncycles,' cycles.'
!!   now we have a flat propagator matrix obtained as a perturbation
!!   by emulating the set of stretches actually sampled, we get an idea of the initial state bias 
! 555 format(i10,3(g16.8,1x))
!    write(ilog,*) 'Now computing initial state bias based on uniform propagator ...'
!    thr_rhlper(1:nbasins,1) = 1.0/(1.0*nbasins)
!    didtb = 0
!    if (allocated(trbrkslst).EQV..false.) then
!      didtb = 1
!      allocate(trbrkslst(2))
!      trbrkslst(:) = cstored
!      ntbrks2 = 1
!    end if
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    ii = 3
!    j = 50
!    k = 2
!    if (write_weights.EQV..true.) then
!      i = 3
!      allocate(wvec(cstored))
!!      call propagate_stretches(it,nbasins,k,cstored,trbrkslst(1:ntbrks2),ntbrks2,thr_rhlper(1:nbasins,1),j,wvec,cstored,i,ii,tpi)
!      call propagate_stretches2(it,nbasins,k,cstored,trbrkslst(1:ntbrks2),ntbrks2,j,thr_rhlper(1:nbasins,1),wvec,cstored,i,ii,tpi)
!    else
!      i = 0
!      allocate(wvec(aone))
!      call propagate_stretches(it,nbasins,k,cstored,trbrkslst(1:ntbrks2),ntbrks2,thr_rhlper(1:nbasins,1),j,wvec,aone,i,ii,tpi)
!    end if
#ifdef ENABLE_THREADS
!! !$OMP BARRIER
!! !$OMP MASTER
#endif
!    if (didtb.eq.1) then
!      ntbrks2 = 0
!      deallocate(trbrkslst)
!    end if
!    write(ilog,564) minval(thr_rhlper(1:nbasins,2))*nbasins,maxval(thr_rhlper(1:nbasins,2))*nbasins
!    it(1:nbasins)%nodewt = it(1:nbasins)%nmbrs/(1.0*cstored)/thr_rhlper(1:nbasins,2)
!    write(ilog,*) 'Done.'
!    deallocate(thr_rhlper)
!    deallocate(thr_hlper)
#ifdef ENABLE_THREADS
!! !$OMP END MASTER
!! !$OMP BARRIER
#endif
!    
!!    do i=1,nbasins
!!      write(*,555) i,it(i)%sums(:,1)/(1.0*it(i)%nmbrs),thr_rhlper(i,1:2)
!!    end do
!
!!    fraccer = 300.0/(1.0*cstored)
!!    thr_rhlper(:,3) =  1.0/(1.0*nbasins) ! it(1:nbasins)%nodewt
!!    thr_rhlper(:,4) = thr_rhlper(:,3)
!!    normeo = HUGE(normeo)
!!    do ncycles=1,10
!!      call propagate_stretches(it,nbasins,iv1,cstored,trbrkslst(1:ntbrks2),ntbrks2,thr_rhlper(:,3),thr_rhlper(:,2),aone)
!!      do i=1,nbasins
!!        do j=1,it(i)%nb
!!          if (abs(it(i)%fewtsnb(j,2)-1.0*it(i)%wghtsnb(j,1)).gt.20.) then
!!            write(*,*) i,j,it(i)%fewtsnb(j,2),it(i)%wghtsnb(j,1)
!!          end if
!!        end do
!!        if (i.le.10) write(*,*) i,thr_rhlper(i,3)
!!      end do
!!      thr_rhlper(:,4) = thr_rhlper(:,3)
!!      normerr = 0.0
!!      do i=1,nbasins
!!        do j=1,it(i)%nb
!!          maxerr = it(i)%fewtsnb(j,2) - 1.0*it(i)%wghtsnb(j,1)
!!          normerr = normerr + abs(maxerr)
!!!          thr_rhlper(i,4) = thr_rhlper(i,4)*(1.0 - fraccer*maxerr)
!!          thr_rhlper(it(i)%lstnb(j),4) = thr_rhlper(it(i)%lstnb(j),4)*(1.0 - fraccer*maxerr)
!!        end do
!!      end do
!!!      thr_rhlper(:,4) = thr_rhlper(:,4) + minval(thr_rhlper(:,4)) + 1.0e-9 ! max(thr_rhlper(:,4),1.0e-9)
!!      if (normerr.lt.normeo) then
!!        thr_rhlper(:,3) = thr_rhlper(:,4)/sum(thr_rhlper(:,4))
!!        normeo = normerr
!!        write(*,*) 'Productive',ncycles,normerr
!!      else
!!        fraccer = fraccer/2.0
!!        write(*,*) 'Reduction',normerr,fraccer
!!      end if
!!    end do

!!    it(:)%nodewt = thr_rhlper(:,3)
!!    deallocate(pbase)
!!    deallocate(hit)
!!    deallocate(iv1)
!
!  else if (cequil.eq.5) then
!
!    allocate(thr_rhlper(nbasins,4))
!    write(ilog,*) 'Now reequilibrating network ...'
!    k = 3
!    j = 2
!    call equilibrate_MSM(it,nbasins,thr_rhlper(:,1),j,k,thresh,iterations,dofwr,maxtime_eqmsm,tpi)
!    it(1:nbasins)%nodewt = thr_rhlper(:,2)
!    write(ilog,*) 'Done after ',iterations,' iterations.'

!    allocate(hit(nbasins))
!    allocate(pbase(nsccs,2))
!    fraccer = 0.4
!!    selffac = 0.0
!!    do while (1.eq.1)
!!      if (abs(erfc(selffac*diffth)-0.5).le.1.0e-6) exit
!!      selffac = selffac + 0.1*(erfcnormffac*diffth) - 0.5)
!!    end do
!!    do j=1,nbasins
!!      it(j)%nodewt = it(j)%nodewt*(1.0 + 0.2 - 0.4*random())
!!    end do
!!    it(1:nbasins)%nodewt = it(1:nbasins)%nodewt/sum(it(1:nbasins)%nodewt)
!    thr_rhlper(:,1) = it(1:nbasins)%nodewt
!    normer = 1.0/(sum(dhist)+distbins/(1.0*it(1)%nmbrs))
!    lasterr = HUGE(lasterr)
!    selffac = 1.0/(1.5*maxval(it(1:nbasins)%nb))
!
!    do ncycles=1,1 ! 250
!!     generate new T
!      do i=1,nbasins
!        normeo = 0.0
!        do j=1,it(i)%nb
!          if (it(i)%fewtsnb(j,1).gt.0.0) then
!!            call clusters_mean_diff(it(i),it(it(i)%lstnb(j)),maxerr)
!!            k = min(max(1,ceiling(maxerr*2.0/cradius)),distbins)
!!            if (k.eq.1) then
!!              lipval = dhist(1) + 1.0 + (dhist(2)-dhist(1))*(maxerr - 0.25*cradius)/(0.5*cradius)
!!            else if (k.eq.distbins) then
!!              lipval = dhist(distbins) + 1
!!            else if (maxerr.gt.(k*0.5*cradius - 0.25*cradius)) then
!!              lipval = dhist(k)+1.0/(1.0*it(1)%nmbrs)+(dhist(k+1)-dhist(k))*(maxerr - k*0.5*cradius+0.25*cradius)/(0.5*cradius)
!!            else
!!              lipval = dhist(k)+1.0/(1.0*it(1)%nmbrs)-(dhist(k-1)-dhist(k))*(maxerr - k*0.5*cradius+0.25*cradius)/(0.5*cradius)
!!            end if
!!             write(*,*) maxerr,erfc(selffac*maxerr),lipval
!!            it(i)%fewtsnb(j,1) = it(it(i)%lstnb(j))%nodewt*(lipval*normer + (it(i)%wghtsnb(j,1)/(1.0*it(i)%nmbrs))) ! &
!! & erfc(selffac*maxerr)
!            if (i.eq.it(i)%lstnb(j)) then
!!              it(i)%fewtsnb(j,1) = it(i)%fewtsnb(j,1)*2.0 
!              it(i)%fewtsnb(j,1) = (1.0 - selffac)*0.5
!            else
!!              it(i)%fewtsnb(j,1) = it(i)%fewtsnb(j,1)*(it(it(i)%lstnb(j))%nodewt + it(i)%nodewt)/it(it(i)%lstnb(j))%nodewt
!              it(i)%fewtsnb(j,1) = selffac*it(it(i)%lstnb(j))%nodewt/(it(it(i)%lstnb(j))%nodewt + it(i)%nodewt)
!            end if
!          end if
!        end do
!      end do
!
!      k = 6
!      thr_rhlper(:,2) = it(1:nbasins)%nodewt
!      j = 1
!      call equilibrate_MSM(it,nbasins,thr_rhlper(:,2),j,k,thresh,iterations,dofwr,maxtime_eqmsm,tpi)
!! 
!!     check consistency and adjust
!      normerr = sum(abs(it(1:nbasins)%nodewt-thr_rhlper(:,1)))/(1.0*nbasins)
!      if ((normerr.gt.lasterr).AND.(ncycles.gt.5)) exit
!      lasterr = normerr
!      it(1:nbasins)%nodewt = thr_rhlper(:,1) ! it(1:nbasins)%nodewt - fraccer*(thr_rhlper(1:nbasins,1) - it(1:nbasins)%nodewt)
!      write(*,*) ncycles,normerr
!
!    end do

!    deallocate(pbase)
!    deallocate(hit)
!
!  else if (cequil.eq.4) then
!
!    allocate(thr_rhlper(nbasins,4))
!    allocate(pbase(nsccs,2))
!    fraccer = 1.0
!
!    normerr = HUGE(normerr)
!    do ncycles=1,100000
!      thr_rhlper(:,2) = it(:)%nodewt
!!     compute mutual consistency of T matrix
!      lasterr = normerr
!      normerr = 0.0
!      do i=1,nbasins
!        normeo = 0.0
!        do j=1,it(i)%nb
!          if (it(i)%fewtsnb(j,1).gt.0.0) then
!            normeo = normeo + it(it(i)%lstnb(j))%nodewt*min((2.0/cradius),sqrt(1.0/it(i)%lensnb(j,1)))
!          end if
!        end do
!        normeo = 1.0/normeo
!        do j=1,it(i)%nb
!          if (it(i)%fewtsnb(j,1).gt.0.0) then
!            maxerr = it(i)%fewtsnb(j,1)*thr_rhlper(i,3) - &
! &    normeo*it(it(i)%lstnb(j))%nodewt*min((2.0/cradius),sqrt(1.0/it(i)%lensnb(j,1)))
!            thr_rhlper(it(i)%lstnb(j),2) = max(0.0, thr_rhlper(it(i)%lstnb(j),2)*&
! &         (1.0 + fraccer*maxerr/(normeo*min((2.0/cradius),sqrt(1.0/it(i)%lensnb(j,1))))))
!            normerr = normerr + abs(maxerr)
!          end if
!        end do
!      end do
!      if ((normerr.gt.lasterr).AND.(ncycles.gt.1)) then
!        write(*,*) ncycles,lasterr
!        exit
!      end if
!      if (mod(ncycles,1000).eq.0) write(*,*) ncycles,normerr
!      maxerr = sum(thr_rhlper(:,2))
!      it(:)%nodewt = thr_rhlper(:,2)/maxerr
!    end do
!    deallocate(pbase)
!
  end if
!
! renorm the total(!) in case something did not converge properly
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    maxerr = 1.0/sum(it(1:nbasins)%nodewt(2))
    it(1:nbasins)%nodewt(2) = maxerr*it(1:nbasins)%nodewt(2)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
  if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    maxerr = 1.0/sum(it(1:nbasins)%nodewt(3))
    it(1:nbasins)%nodewt(3) = maxerr*it(1:nbasins)%nodewt(3)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
!
  if (write_weights.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP MASTER
    allocate(thr_hlper(2,1))
    i = 2
    call get_freeunits(i,thr_hlper(:,1))
!$OMP END MASTER
!$OMP BARRIER
#endif
    if ((which_timeflow.eq.1).OR.(which_timeflow.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
#ifdef ENABLE_MPI
      tl = 3
      call int2str(myrank,nod,tl)
      intfile = 'N_'//nod//'_'//basename(1:bleng)//'_FWD.fwts'
#else
      intfile = basename(1:bleng)//'_FWD.fwts'
#endif
      call strlims(intfile,i,j)
      inquire(file=intfile(i:j),exist=staterr)
      if (staterr.EQV..true.) then
#ifdef ENABLE_THREADS
        ii = thr_hlper(1,1)
#else
        ii = freeunit()
#endif
        open (unit=ii,file=intfile(i:j),status='old')
        close(unit=ii,status='delete')
      end if
#ifdef ENABLE_THREADS
      ii = thr_hlper(1,1)
#else
      ii = freeunit()
#endif
      open (unit=ii,file=intfile(i:j),status='new')
 566 format(i14,1x,g19.11E3)
      do j=1,cstored
        write(ii,566) j,cstored*it(sconnect(j,4))%nodewt(2)/(1.0*it(sconnect(j,4))%nmbrs)
      end do
      close(unit=ii)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
    if ((which_timeflow.eq.2).OR.(which_timeflow.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
#ifdef ENABLE_MPI
      tl = 3
      call int2str(myrank,nod,tl)
      intfile = 'N_'//nod//'_'//basename(1:bleng)//'_BWD.fwts'
#else
      intfile = basename(1:bleng)//'_BWD.fwts'
#endif
      call strlims(intfile,i,j)
      inquire(file=intfile(i:j),exist=staterr)
      if (staterr.EQV..true.) then
#ifdef ENABLE_THREADS
        ii = thr_hlper(2,1)
#else
        ii = freeunit()
#endif
        open (unit=ii,file=intfile(i:j),status='old')
        close(unit=ii,status='delete')
      end if
#ifdef ENABLE_THREADS
      ii = thr_hlper(2,1)
#else
      ii = freeunit()
#endif
      open (unit=ii,file=intfile(i:j),status='new')
      do j=1,cstored
        write(ii,566) j,cstored*it(sconnect(j,4))%nodewt(3)/(1.0*it(sconnect(j,4))%nmbrs)
      end do
      close(unit=ii)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
  if (write_weights.EQV..true.) then
!$OMP SINGLE
    deallocate(thr_hlper)
!$OMP END SINGLE 
  end if
#endif
!
end
!
!----------------------------------------------------------------------------------------------------
!
! a subroutine to add links in controlled fashion
!
! dhist(distbins) is the histogram of pairwise distances for the native transitions
!
subroutine augment_network_dev(it,nbasins,distbins,dhist)
!
  use clusters
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: nbasins,distbins
  type(t_scluster):: it(nbasins)
  RTYPE, INTENT(IN):: dhist(distbins)
!
  integer i,j,k,ii,ncycles,nfails,iterations,nsccs
  logical staterr
  RTYPE diffth,norm,lipval,rval
!
  write(ilog,*) 'Now augmenting network ...'
  k = 1
  call tarjan_SCC(it,nbasins,k)
  nsccs = maxval(it(1:nbasins)%inscc)
  ncycles = 0 
  nfails = 0
! fix missing self transitions and reverse edges
  do i=1,nbasins
    staterr = .false.
    do j=1,it(i)%nb
      if (it(i)%lstnb(j).eq.i) then
        staterr = .true.
      else if (it(i)%wghtsnb(j,1).le.0) then
!        write(*,*) 'Missing link from ',it(i)%sums(:,1)/(1.0*it(i)%nmbrs),' to ',&
! &it(it(i)%lstnb(j))%sums(:,1)/(1.0*it(it(i)%lstnb(j))%nmbrs)
        it(i)%fewtsnb(j,1) = caddlinkwt
        it(i)%lensnb(j,1) = it(it(i)%lstnb(j))%lensnb(it(i)%map(j),1)
        nfails = nfails + 1
!        write(*,*) 'Aug',i,it(i)%lstnb(j),it(i)%lensnb(j,1)
      else
!        write(*,*) 'Exi',i,it(i)%lstnb(j),it(i)%lensnb(j,1)
      end if
    end do
    if (staterr.EQV..false.) then
      ncycles = ncycles + 1
      it(i)%nb = it(i)%nb + 1
      if (it(i)%nb.gt.it(i)%nbalsz) then
        call scluster_resizenb(it(i))
      end if
      it(i)%map(it(i)%nb) = it(i)%nb
      it(i)%lstnb(it(i)%nb) = i
      it(i)%wghtsnb(it(i)%nb,1:2) = 0
      it(i)%lensnb(it(i)%nb,1:2) = 0.5*cradius
      it(i)%fewtsnb(it(i)%nb,1:2) = caddlinkwt
    end if
  end do
!
! scan the histogram to find the width cutoff
  norm = sum(dhist)
  lipval = 0
  do i=1,distbins
    lipval = lipval + dhist(i)
    if (lipval/norm.gt.caddlinks) then
      diffth = 0.5*(i-1)*cradius*(lipval - norm*caddlinks)/dhist(i) + 0.5*i*cradius*(1.0 - (lipval - norm*caddlinks)/dhist(i))
      exit
    end if
  end do
!
  iterations = 0
  do i=1,nbasins
    do j=1,it(i)%nb
!                  write(*,*) 'dicomp',sqrt(it(i)%lensnb(j,1)),rval!,it(i)%sums(:,1)/(1.0*it(i)%nmbrs)-&
!  &it(it(i)%lstnb(j))%sums(:,1)/(1.0*it(it(i)%lstnb(j))%nmbrs)
      if (it(i)%lstnb(j).eq.i) cycle
      do k=1,it(it(i)%lstnb(j))%nb
        if (it(it(i)%lstnb(j))%lstnb(k).eq.i) cycle
        staterr = .false.
        do ii=1,it(i)%nb
          if (it(it(i)%lstnb(j))%lstnb(k).eq.it(i)%lstnb(ii)) then
            staterr = .true.
            exit
          end if
        end do
        if (staterr.EQV..false.) then
          call clusters_mean_diff(it(i),it(it(it(i)%lstnb(j))%lstnb(k)),rval)
          if (rval.le.diffth) then
            iterations = iterations + 1
            it(i)%nb = it(i)%nb + 1
            if (it(i)%nb.gt.it(i)%nbalsz) then
              call scluster_resizenb(it(i))
            end if
!            it(i)%map(it(i)%nb) = it(i)%nb
            it(i)%lstnb(it(i)%nb) = it(it(i)%lstnb(j))%lstnb(k)
            it(i)%wghtsnb(it(i)%nb,1:2) = 0
            it(i)%lensnb(it(i)%nb,1:2) = rval
            it(i)%fewtsnb(it(i)%nb,1:2) = caddlinkwt
!          write(*,*) 'Suggest from ',it(i)%sums(:,1)/(1.0*it(i)%nmbrs),' to ',&
! &it(it(it(i)%lstnb(j))%lstnb(k))%sums(:,1)/(1.0*it(it(it(i)%lstnb(j))%lstnb(k))%nmbrs)
          end if
        end if
      end do
    end do
  end do
! recompute SCCs
  k = 2
  call tarjan_SCC(it,nbasins,k)
  k = maxval(it(1:nbasins)%inscc)
 560 format('Done. Augmentation used a threshold distance of ',g14.5,' (various units) based on the quantile of ',g12.5,'.')
 561 format('Added ',i7,' new links and ',i7,' self transitions, and added ',i7,' missing reverse links.',//,&
 &        'The recalculated number of strongly connected components is ',i7,' (from ',i7,').')
  write(ilog,560) diffth,caddlinks
  write(ilog,561) iterations,ncycles,nfails,k,nsccs
  nsccs = k
!
end
!
!------------------------------------------------------------------------------------------------------
!
! this subroutine takes the set of simulation stretches given by the coarse-grained trajectory (thr_hlper(1:tlen,tlstx)) and
! supplied breaks list (bklst,blen) and emulates them exactly as probability propagations on the network defined by: 
! Tij = Z*it(i)%fewtsnb(k(j),1), where Z is the normalizer (per i) and mode2 is 3
! Tij = Z*it(i)%wghtsnb(k(j),1), where Z is the normalizer (per i) and mode2 is 2
! Tij = Z*it(i)%fewtsnb(k(j),1)*p(j)/(p(j)+p(i)), where Z is the normalizer (per i), p is provided as input (inss), and mode2 is 1
!
! thr_rhlper(1:nbasins,outpx) is a resultant probability distribution weighted correctly in terms of stretch lengths
! if mode is 1 or 2, it(:)fewtsnb(:,2) is additionally populated as a predicted count matrix
! if mode is 2 or 3, wvals is additionally populated as a weight vector per snapshot
! 
subroutine propagate_stretches(it,nbasins,tlstx,tlen,bklst,blen,inss,outpx,wvals,wlen,mode,mode2,tpi)
!
  use clusters
  use iounit
  use threads, ONLY: thr_hlper,thr_rhlper
!
  implicit none
!
  integer, INTENT(IN):: tlen,tlstx,blen,bklst(blen),outpx,mode,mode2,nbasins,tpi,wlen
  RTYPE, INTENT(IN):: inss(nbasins)
  type(t_scluster):: it(nbasins)
  RTYPE wvals(wlen)
!
  integer i,j,iii,ilen,ilo,ihi,stretchi,cbnd(2),sts(4)
  RTYPE incr,msd,msdmax
  RTYPE, ALLOCATABLE:: ptmp(:)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer OMP_GET_NUM_THREADS,k
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, propagate_stretches(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  i = 1
  j = nbasins
  if (tpi.gt.0) then
    k = OMP_GET_NUM_THREADS()
    call threads_bounds(i,j,tpi,k,cbnd)
    j = 4
    call get_thread_loop_bounds(i,j,sts(1:2),sts(3:4),tpi)
    sts(2) = sts(3)
  else
    cbnd(1) = 1
    cbnd(2) = nbasins
    sts(1:2) = cbnd(1:2)
  end if
!$OMP BARRIER
#else
  cbnd(1) = 1
  cbnd(2) = nbasins
  sts(1:2) = cbnd(1:2)
#endif
  allocate(ptmp(nbasins))
  do i=sts(1),sts(2)
    if (mode2.eq.1) then
      ptmp(i) = 0.0
      do j=1,it(i)%nb
        ptmp(i) = ptmp(i) + it(i)%fewtsnb(j,1)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if (mode2.eq.2) then
      ptmp(i) = 1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))
    else if (mode2.eq.3) then
      ptmp(i) = sum(it(i)%fewtsnb(1:it(i)%nb,1))
    else
      write(ilog,*) 'Called propagate_stretches(...) with unsupported mode2 (got ',mode2,'). This is a bug.'
      call fexit()
    end if
    if (ptmp(i).gt.0.0) ptmp(i) = 1.0/ptmp(i) ! normalizer
!
    if (mode2.eq.1) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i)*it(i)%fewtsnb(j,1)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if (mode2.eq.2) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i)*it(i)%wghtsnb(j,1)
      end do
    else if (mode2.eq.3) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i)*it(i)%fewtsnb(j,1)
      end do
    end if
  end do
!
  thr_rhlper(cbnd(1):cbnd(2),outpx) = 0.0
  do iii=0,blen
    if (iii.eq.0) then
      ilen = bklst(1)
      ilo = 1
      ihi = bklst(1)
    else if (iii.eq.blen) then
      ilen = cstored - bklst(iii)
      ilo = bklst(iii) + 1
      ihi = cstored
    else
      ilen = bklst(iii+1) - bklst(iii)
      ilo = bklst(iii) + 1
      ihi = bklst(iii+1)
    end if
    if (ilen.eq.0) cycle
    thr_rhlper(cbnd(1):cbnd(2),4) = 0.0
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
    thr_rhlper(thr_hlper(ilo,tlstx),4) = 1.0
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    do stretchi=1,ilen
      ptmp(:) = 0.0 
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      do i=sts(1),sts(2)
        if (thr_rhlper(i,4).le.0.0) cycle
        do j=1,it(i)%nb
          incr = thr_rhlper(i,4)*it(i)%fewtsnb(j,3)
          ptmp(it(i)%lstnb(j)) = ptmp(it(i)%lstnb(j)) + incr
          if (mode.le.2) it(i)%fewtsnb(j,2) = it(i)%fewtsnb(j,2) + incr
        end do
      end do
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then 
        j = 1
        k = 4
        i = 67
!$OMP BARRIER
        thr_rhlper(cbnd(1):cbnd(2),k) = 0.0
        call thr_combinevec(ptmp,nbasins,j,k,i,tpi)
      else
        thr_rhlper(1:nbasins,4) = ptmp(:)
      end if
#else
      thr_rhlper(1:nbasins,4) = ptmp(:)
#endif
    end do
!    if ((mode.eq.2).OR.(mode.eq.3)) then
!      do stretchi=ilo,ihi
!        wvals(stretchi) = (1.0/(1.0*nbasins) + wbuf)/(phlp2(thr_hlper(stretchi,tlstx)) + wbuf)
!      end do
!    end if
    msd = 0.0
    msdmax = 0.0
    do i=cbnd(1),cbnd(2)
      call snap_to_snap_d(incr,it(thr_hlper(ilo,tlstx))%center,it(i)%center)
      msd = msd + thr_rhlper(i,4)*incr*incr
      msdmax = msdmax + incr*incr/(1.0*nbasins)
    end do
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      i = 0
      j = 1
      k = 5
      call thr_reduction(i,msd,j,k)
      msd = thr_rhlper(j,k)
      k = 3
      call thr_reduction(i,msdmax,j,k)
      msdmax = thr_rhlper(j,k)
    end if
!$OMP BARRIER
#endif
    if (tpi.le.1) write(*,*) ilo,msd,msdmax,msd/ilen,-msdmax*log(max(1.0e-6,1.0 - msd/msdmax))/(1.0*ilen)
    if ((mode.eq.2).OR.(mode.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      thr_rhlper(1,3) = 0.0
      thr_rhlper(1,5) = 0.0
      incr = -msdmax*log(max(1.0e-6,1.0 - msd/msdmax))/(1.0*ilen) ! leading term is msd/ilen
      wvals(ilo:ihi) = incr
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    end if
    thr_rhlper(cbnd(1):cbnd(2),4) = ilen*thr_rhlper(cbnd(1):cbnd(2),4)/(1.0*tlen)
    thr_rhlper(cbnd(1):cbnd(2),outpx) = thr_rhlper(cbnd(1):cbnd(2),outpx) + thr_rhlper(cbnd(1):cbnd(2),4)
  end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  thr_rhlper(1:nbasins,outpx) = thr_rhlper(1:nbasins,outpx)/sum(thr_rhlper(1:nbasins,outpx))
  do stretchi=1,cstored
!    write(*,*) wvals(stretchi),(1.0/(1.0*nbasins))/(thr_rhlper(thr_hlper(stretchi,tlstx),outpx))
!    wvals(stretchi) = wvals(stretchi)*(1.0/(1.0*nbasins))/(thr_rhlper(thr_hlper(stretchi,tlstx),outpx))
  end do
!  write(*,*) minval(wvals),maxval(wvals)
!  do k=1,nbasins
!    do j=1,it(k)%nb
!      if (it(k)%lstnb(j).eq.k) write(*,*) k,it(k)%fewtsnb(j,1)/sum(it(k)%fewtsnb(:,1)),thr_rhlper(k,outpx)
!    end do
!  end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  deallocate(ptmp)
!
end
!
!---------------------------------------------------------------------------------------------------------
!
! 
subroutine propagate_stretches2(it,nbasins,tlstx,tlen,bklst,blen,difflen,inss,wvals,wlen,mode,mode2,tpi)
!
  use clusters
  use iounit
  use threads, ONLY: thr_hlper,thr_rhlper
!
  implicit none
!
  integer, INTENT(IN):: tlen,tlstx,blen,bklst(blen),mode,mode2,nbasins,tpi,wlen,difflen
  RTYPE, INTENT(IN):: inss(nbasins)
  type(t_scluster):: it(nbasins)
  RTYPE wvals(wlen)
!
  integer i,j,k,iii,ilen,ilo,ihi,stretchi,cbnd(2),sts(4),flen(2)
  RTYPE incr,msd(3),tsh,sfac(2)
  RTYPE, ALLOCATABLE:: ptmp(:,:)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer OMP_GET_NUM_THREADS
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, propagate_stretches2(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  i = 1
  j = nbasins
  if (tpi.gt.0) then
    k = OMP_GET_NUM_THREADS()
    call threads_bounds(i,j,tpi,k,cbnd)
    j = 4
    call get_thread_loop_bounds(i,j,sts(1:2),sts(3:4),tpi)
    sts(2) = sts(3)
  else
    cbnd(1) = 1
    cbnd(2) = nbasins
    sts(1:2) = cbnd(1:2)
  end if
!$OMP BARRIER
#else
  cbnd(1) = 1
  cbnd(2) = nbasins
  sts(1:2) = cbnd(1:2)
#endif
  allocate(ptmp(nbasins,2))
  do i=sts(1),sts(2)
    if (mode2.eq.1) then
      ptmp(i,1) = 0.0
      do j=1,it(i)%nb
        ptmp(i,1) = ptmp(i,1) + it(i)%fewtsnb(j,1)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if (mode2.eq.2) then
      ptmp(i,1) = 1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))
    else if (mode2.eq.3) then
      ptmp(i,1) = sum(it(i)%fewtsnb(1:it(i)%nb,1))
    else
      write(ilog,*) 'Called propagate_stretches(...) with unsupported mode2 (got ',mode2,'). This is a bug.'
      call fexit()
    end if
    if (ptmp(i,1).gt.0.0) ptmp(i,1) = 1.0/ptmp(i,1) ! normalizer
!
    if (mode2.eq.1) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i,1)*it(i)%fewtsnb(j,1)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if (mode2.eq.2) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i,1)*it(i)%wghtsnb(j,1)
      end do
    else if (mode2.eq.3) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i,1)*it(i)%fewtsnb(j,1)
      end do
    end if
  end do
!
  thr_rhlper(cbnd(1):cbnd(2),2:3) = 0.0
  do iii=0,blen
    if (iii.eq.0) then
      ilen = bklst(1)
      ilo = 1
      ihi = bklst(1)
    else if (iii.eq.blen) then
      ilen = cstored - bklst(iii)
      ilo = bklst(iii) + 1
      ihi = cstored
    else
      ilen = bklst(iii+1) - bklst(iii)
      ilo = bklst(iii) + 1
      ihi = bklst(iii+1)
    end if
    if (ilen.eq.0) cycle
    tsh = 1.0*(1.0 - exp(-(1.0*ilen)/cequilbuf))
    j = 4
    k = mode2
    thr_rhlper(cbnd(1):cbnd(2),6) = 0.0
    flen(:) = 0
    if (thr_hlper(ilo,tlstx).eq.thr_hlper(ihi,tlstx)) cycle
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
!    thr_rhlper(thr_hlper(ilo,tlstx),6) = 1.0
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!    call flow_MSM(it,nbasins,thr_hlper(ihi,tlstx),thr_rhlper(:,6),inss,j,k,tsh,flen(1),ilen,tpi)
!    thr_rhlper(cbnd(1):cbnd(2),6) = 0.0
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
!    thr_rhlper(thr_hlper(ihi,tlstx),6) = 1.0
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!    call flow_MSM(it,nbasins,thr_hlper(ilo,tlstx),thr_rhlper(:,6),inss,j,k,tsh,flen(2),ilen,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    thr_rhlper(cbnd(1):cbnd(2),4:5) = 0.0
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
!    write(*,*) ilen,flen(1:2)
    thr_rhlper(thr_hlper(ilo,tlstx),4) = 1.0
    thr_rhlper(thr_hlper(ihi,tlstx),5) = 1.0
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    sfac(1:2) = 1.0 ! (1.0*ilen)/(1.0*flen(1:2))
    do stretchi=1,ilen ! maxval(flen)
      ptmp(:,:) = 0.0 
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!      if (stretchi.le.flen(1)) then
        do i=sts(1),sts(2)
          if (thr_rhlper(i,4).le.0.0) cycle
          do j=1,it(i)%nb
            incr = thr_rhlper(i,4)*it(i)%fewtsnb(j,3)
            ptmp(it(i)%lstnb(j),1) = ptmp(it(i)%lstnb(j),1) + incr
          end do
        end do
!      end if
!      if (stretchi.le.flen(2)) then
        do i=sts(1),sts(2)
          if (thr_rhlper(i,5).le.0.0) cycle
          do j=1,it(i)%nb
            incr = thr_rhlper(i,5)*it(i)%fewtsnb(j,3)
            ptmp(it(i)%lstnb(j),2) = ptmp(it(i)%lstnb(j),2) + incr
          end do
        end do
!      end if
!
      if (stretchi.eq.min(difflen,ilen)) then ! maxval(flen))) then
        msd(:) = 0.0
        do i=cbnd(1),cbnd(2)
          call snap_to_snap_d(incr,it(thr_hlper(ilo,tlstx))%center,it(i)%center)
          msd(1) = msd(1) + thr_rhlper(i,4)*incr*incr
          msd(3) = msd(3) + incr*incr/(1.0*nbasins)
          call snap_to_snap_d(incr,it(thr_hlper(ihi,tlstx))%center,it(i)%center)
          msd(2) = msd(2) + thr_rhlper(i,5)*incr*incr
        end do
      end if
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then 
        j = 1
        i = 67
!$OMP BARRIER
!        if (stretchi.le.flen(1)) then
          k = 4
          thr_rhlper(cbnd(1):cbnd(2),k) = 0.0
          call thr_combinevec(ptmp(:,1),nbasins,j,k,i,tpi)
!        end if
!        if (stretchi.le.flen(2)) then
          k = 5
          thr_rhlper(cbnd(1):cbnd(2),k) = 0.0
          call thr_combinevec(ptmp(:,2),nbasins,j,k,i,tpi)
!        end if
      else
!        if (stretchi.le.flen(1)) thr_rhlper(1:nbasins,4) = ptmp(:,1)
!        if (stretchi.le.flen(2)) thr_rhlper(1:nbasins,5) = ptmp(:,2)
        thr_rhlper(1:nbasins,4) = ptmp(:,1)
        thr_rhlper(1:nbasins,5) = ptmp(:,2)
      end if
#else
!      if (stretchi.le.flen(1)) thr_rhlper(1:nbasins,4) = ptmp(:,1)
!      if (stretchi.le.flen(2)) thr_rhlper(1:nbasins,5) = ptmp(:,2)
      thr_rhlper(1:nbasins,4) = ptmp(:,1)
      thr_rhlper(1:nbasins,5) = ptmp(:,2)
#endif
!      if (stretchi.le.flen(1)) then
        thr_rhlper(cbnd(1):cbnd(2),2) = thr_rhlper(cbnd(1):cbnd(2),2) + sfac(1)*thr_rhlper(cbnd(1):cbnd(2),4)
!      end if
!      if (stretchi.le.flen(2)) then
        thr_rhlper(cbnd(1):cbnd(2),3) = thr_rhlper(cbnd(1):cbnd(2),3) + sfac(2)*thr_rhlper(cbnd(1):cbnd(2),5)
!      end if
    end do
!
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
!$OMP BARRIER
      i = 0
      j = 1
      k = 5
      if (tpi.le.1) thr_rhlper(1,5) = 0.0
      call thr_reduction(i,msd(1),j,k)
      msd(1) = thr_rhlper(j,k)
      if (tpi.le.1) thr_rhlper(1,5) = 0.0
      call thr_reduction(i,msd(2),j,k)
      msd(2) = thr_rhlper(j,k)
      if (tpi.le.1) thr_rhlper(1,5) = 0.0
      call thr_reduction(i,msd(3),j,k)
      msd(3) = thr_rhlper(j,k)
    end if
!$OMP BARRIER
#endif
 5762 format(i8,5(g14.6,1x))
    if (tpi.le.1) write(*,5762) ilo,msd(:),-msd(3)*log(max(1.0e-6,1.0 - msd(1:2)/msd(3)))/(1.0*min(ilen,difflen))
    if ((mode.eq.2).OR.(mode.eq.3)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      incr = log(max(1.0e-6,1.0 - msd(2)/msd(3)))/log(max(1.0e-6,1.0 - msd(1)/msd(3))) ! leading term is msd/ilen
      do i=ilo,ihi
        wvals(ilo:ihi) = 1.0 ! incr
      end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  thr_rhlper(1:nbasins,2) = thr_rhlper(1:nbasins,2)/sum(thr_rhlper(1:nbasins,2))
  thr_rhlper(1:nbasins,3) = thr_rhlper(1:nbasins,3)/sum(thr_rhlper(1:nbasins,3))
!  do i=1,nbasins
!    write(*,5762) i,it(i)%sums(1,1)/(1.0*it(i)%nmbrs),thr_rhlper(i,2:3),thr_rhlper(i,3)/thr_rhlper(i,2)
!  end do
  do stretchi=1,cstored
!    write(*,*) wvals(stretchi),(1.0/(1.0*nbasins))/(thr_rhlper(thr_hlper(stretchi,tlstx),outpx))
    wvals(stretchi) = wvals(stretchi)*((thr_rhlper(thr_hlper(stretchi,tlstx),3))/(thr_rhlper(thr_hlper(stretchi,tlstx),2)))**2
  end do
!  write(*,*) minval(wvals),maxval(wvals)
!  do k=1,nbasins
!    do j=1,it(k)%nb
!      if (it(k)%lstnb(j).eq.k) write(*,*) k,it(k)%fewtsnb(j,1)/sum(it(k)%fewtsnb(:,1)),thr_rhlper(k,outpx)
!    end do
!  end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  deallocate(ptmp)
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! this subroutine computes the steady-state of the MSM defined as:
! modes 1,4:
! Tij = Z*it(i)%fewtsnb(k(j),1)*p(j)/(p(j)+p(i)), where Z is the normalizer (per i) and p is provided as input (inss)
! modes 2,5:
! Tij = Z*it(i)%wghtsnb(k(j),1), where Z is the normalizer
! modes 3,6:
! Tij = Z*it(i)%fewtsnb(k(j),1), where Z is the normalizer
!
! in modes 4-6, inss is taken as initial guess , otherwise nodewt(1)
! if there is more than one SCC, it assumes that it(:)%nodewt(1) respects their relative weights
!
! dofwr controls whether to use forward or backward time matrix
!
! thr_hlper and thr_rhlper must be allocated with sizes (nbasins,1) and (nbasins,5) at least
!
subroutine equilibrate_MSM(it,nbasins,inss,outpx,mode,thresh,iterations,dofwr,maxtime,tpi)
!
  use clusters
  use iounit
#ifdef ENABLE_THREADS
  use threads
#else
  use threads, ONLY: thr_rhlper,thrdat
#endif
!
  implicit none
!
  integer, INTENT(IN):: mode,nbasins,tpi,outpx
  RTYPE, INTENT(IN):: thresh,inss(nbasins)
  logical, INTENT(IN):: dofwr
  integer (KIND=8), INTENT(IN):: maxtime
!
  type(t_scluster) it(nbasins)
  integer i,j,iterations,nfails,maxscc,cbnd(2),sts(4),fbix
  RTYPE normeo,normerr,maxerr
  logical notdone,staterr
  integer(KIND=8) tts(2)
!
  RTYPE, ALLOCATABLE:: ptmp(:)
  integer, ALLOCATABLE:: hit(:)
  RTYPE, ALLOCATABLE:: pbase(:,:)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer OMP_GET_NUM_THREADS,k
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, equilibrate_MSM(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  i = 1
  j = nbasins
  if (tpi.gt.0) then
    k = OMP_GET_NUM_THREADS()
    call threads_bounds(i,j,tpi,k,cbnd)
    j = 4
    call get_thread_loop_bounds(i,j,sts(1:2),sts(3:4),tpi)
    sts(2) = sts(3)
  else
    cbnd(1) = 1
    cbnd(2) = nbasins
    sts(1:2) = cbnd(1:2)
  end if
!$OMP MASTER
  thr_hlper(1:nbasins,1) = 0
  call System_Clock(count=tts(1))
!$OMP END MASTER
!$OMP BARRIER
#else
  cbnd(1) = 1
  cbnd(2) = nbasins
  sts(1:2) = cbnd(1:2)
  call System_Clock(count=tts(1))
#endif
!
  if (dofwr.EQV..true.) then
    fbix = 3
  else
    fbix = 4
  end if
  maxscc = size(csccwts)
  allocate(ptmp(nbasins))
  allocate(hit(nbasins))
  allocate(pbase(maxscc,2))
  pbase(:,:) = 0.0
  do i=sts(1),sts(2)
    pbase(it(i)%inscc,1) = pbase(it(i)%inscc,1) + it(i)%nodewt(1)
    if ((mode.eq.1).OR.(mode.eq.4)) then
      ptmp(i) = 0.0
      do j=1,it(i)%nb
        if (it(it(i)%lstnb(j))%inscc.eq.it(i)%inscc) then
          ptmp(i) = ptmp(i) + it(i)%fewtsnb(j,fbix-2)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
        end if
      end do
    else if ((mode.eq.2).OR.(mode.eq.5)) then
      ptmp(i) = 0.0
      do j=1,it(i)%nb
        if (it(it(i)%lstnb(j))%inscc.eq.it(i)%inscc) then
          ptmp(i) = ptmp(i) + 1.0*it(i)%wghtsnb(j,fbix-2)
        end if
      end do
    else if ((mode.eq.3).OR.(mode.eq.6)) then
      ptmp(i) = 0.0
      do j=1,it(i)%nb
        if (it(it(i)%lstnb(j))%inscc.eq.it(i)%inscc) then
          ptmp(i) = ptmp(i) + it(i)%fewtsnb(j,fbix-2)
        end if
      end do
    else
      write(ilog,*) 'Called equilibrate_MSM(...) with unsupported mode (got ',mode,'). This is a bug.'
      call fexit()
    end if
    if (ptmp(i).gt.0.0) ptmp(i) = 1.0/ptmp(i) ! normalizer
!
    if ((mode.eq.1).OR.(mode.eq.4)) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,fbix) = ptmp(i)*it(i)%fewtsnb(j,fbix-2)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if ((mode.eq.2).OR.(mode.eq.5)) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,fbix) = ptmp(i)*it(i)%wghtsnb(j,fbix-2)
      end do
    else if ((mode.eq.3).OR.(mode.eq.6)) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,fbix) = ptmp(i)*it(i)%fewtsnb(j,fbix-2)
      end do
    end if
  end do
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    j = 1
    k = 4
    i = 69
!$OMP BARRIER
!$OMP SINGLE
    thr_rhlper(1:maxscc,4) = 0.0
!$OMP END SINGLE
    call thr_combinevec(pbase(:,1),maxscc,j,k,i,tpi)
  else
    thr_rhlper(1:maxscc,4) = pbase(:,1)
  end if
#else
  thr_rhlper(1:maxscc,4) = pbase(:,1)
#endif
!
! initial guess
  if (mode.le.3) then
    thr_rhlper(cbnd(1):cbnd(2),outpx) = it(cbnd(1):cbnd(2))%nodewt(1) ! raw sampling weight
  else
    thr_rhlper(cbnd(1):cbnd(2),outpx) = inss(cbnd(1):cbnd(2))
  end if
  notdone = .true.
  iterations = 0
  nfails = 0
  normeo = 1.0
  staterr = .false.
  do while (notdone.EQV..true.)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif 
    iterations = iterations + 1
    ptmp(:) = 0.0
    pbase(:,2) = 0.0
    hit(:) = 0
    do i=sts(1),sts(2)
      do j=1,it(i)%nb
        if ((it(it(i)%lstnb(j))%inscc.eq.it(i)%inscc).AND.(it(i)%wghtsnb(j,fbix-2).gt.0)) then
          ptmp(it(i)%lstnb(j)) = ptmp(it(i)%lstnb(j)) + thr_rhlper(i,outpx)*it(i)%fewtsnb(j,fbix)
          hit(it(i)%lstnb(j)) = 1
        end if
      end do
    end do
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then 
      j = 1
      k = 3
      i = 67
      thr_rhlper(cbnd(1):cbnd(2),3) = 0.0
      thr_hlper(cbnd(1):cbnd(2),1) = 0   
      call thr_combinevec(ptmp,nbasins,j,k,i,tpi)
      k = 1
      call thr_combineivec(hit,nbasins,j,k,i,tpi)
    else
      thr_hlper(1:nbasins,1) = hit(:)
      thr_rhlper(1:nbasins,3) = ptmp(:)
    end if
#else
    thr_rhlper(1:nbasins,3) = ptmp(:)
#endif
!
    if (maxscc.gt.1) then
      do i=cbnd(1),cbnd(2)
#ifdef ENABLE_THREADS
        if (thr_hlper(i,1).eq.0) then
#else
        if (hit(i).eq.0) then ! SCCs of size 1
#endif
          thr_rhlper(i,3) = it(i)%nodewt(1) ! unreachable
        end if
        pbase(it(i)%inscc,2) = pbase(it(i)%inscc,2) + thr_rhlper(i,3)
        thr_rhlper(i,3) = thr_rhlper(it(i)%inscc,4)*thr_rhlper(i,3)
      end do
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        j = 1
        k = 4
        i = 69
!$OMP BARRIER
!$OMP SINGLE
        thr_rhlper(1:maxscc,4) = 0.0
!$OMP END SINGLE
        call thr_combinevec(pbase(:,2),maxscc,j,k,i,tpi)
        pbase(:,2) = thr_rhlper(1:maxscc,4)
      end if
#endif
!
      do i=1,maxscc
        if (pbase(i,2).gt.0.0) pbase(i,2) = 1.0/pbase(i,2)
      end do
      do i=cbnd(1),cbnd(2)
        thr_rhlper(i,3) = thr_rhlper(i,3)*pbase(it(i)%inscc,2)
      end do
    end if
!
    if (mod(iterations,10).eq.0) then
      maxerr = maxval(abs(thr_rhlper(cbnd(1):cbnd(2),outpx)-thr_rhlper(cbnd(1):cbnd(2),3)))
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        i = 1
        j = 1
        k = 5
!$OMP BARRIER
!$OMP SINGLE
        thr_rhlper(j,k) = 0.0
!$OMP END SINGLE
        call thr_reduction(j,maxerr,j,k)
        maxerr = thr_rhlper(j,k)
      end if
#endif
      if (maxerr.le.thresh) notdone = .false.
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif 
    if ((mod(iterations,100).eq.0).AND.(notdone.EQV..true.)) then
      staterr = .false.
      normerr = sum(abs(thr_rhlper(cbnd(1):cbnd(2),3)-thr_rhlper(cbnd(1):cbnd(2),outpx)))/(1.0*nbasins)
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      call System_Clock(count=tts(2))
      if ((tts(2)-tts(1)).gt.maxtime) then
        write(ilog,*) 'Aborting after ',iterations,' iterations due to maximum time of ',(1.0*maxtime)/thrdat%rate,&
 &' [s] being exceeded.'
#ifdef ENABLE_THREADS
!       set a persistent mark for all threads to exit
        thr_hlper(1,1) = -100 ! overwrites hit(1) 
#else
        exit
#endif     
      end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
      if (tpi.gt.0) then
        i = 0
        j = 1
        k = 5
!$OMP BARRIER
!$OMP SINGLE
        thr_rhlper(j,k) = 0.0
!$OMP END SINGLE
        call thr_reduction(i,normerr,j,k)
        normerr = thr_rhlper(j,k)
      end if
!$OMP BARRIER
#endif
      thr_rhlper(cbnd(1):cbnd(2),5) = abs(thr_rhlper(cbnd(1):cbnd(2),3)-thr_rhlper(cbnd(1):cbnd(2),outpx))
!     max a prediction for when the error would reach 0
      if (normerr.ge.normeo) then ! error increases
        nfails = nfails + 1
      else if (normerr/(0.01*(normeo-normerr)).gt.1.0e8) then ! error is not decreasing significantly
        nfails = nfails + 1
      end if 
      if (nfails.gt.100) then
        staterr = .true.
        exit
      end if
      normeo = normerr
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    thr_rhlper(cbnd(1):cbnd(2),outpx) = thr_rhlper(cbnd(1):cbnd(2),3)
#ifdef ENABLE_THREADS
    if (thr_hlper(1,1).eq.-100) exit
#endif
    if (notdone.EQV..false.) exit
  end do
! 
 457 format('Node (cluster) ',i10,' has unresolvable increment in probability of ',g14.6,' with original population of ',g14.6,'.')
  if (staterr.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
 458 format('Some nodes could not be reequilibrated, which is likely due to periodicity in the network (list follows).')
    write(ilog,458)
    do i=1,nbasins
      if (abs(thr_rhlper(i,5)).gt.thresh) write(ilog,457) i,abs(thr_rhlper(i,5)),it(i)%nodewt(1)
    end do
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
!
 459 format('Network imbalance measure (violation of detailed balance) is ',g14.6,' (range 0.0 to 1.0).')
! measure network imbalance
  normeo = 0.0
  normerr = 0.0
  do i=sts(1),sts(2)
    do j=1,it(i)%nb
      if (it(i)%lstnb(j).eq.i) then
        normeo = normeo + thr_rhlper(i,outpx)*it(i)%fewtsnb(j,fbix)
      else if (it(i)%wghtsnb(j,fbix-2).gt.0) then
        normerr = normerr + abs(thr_rhlper(i,outpx)*it(i)%fewtsnb(j,fbix) - &
 &                thr_rhlper(it(i)%lstnb(j),outpx)*it(it(i)%lstnb(j))%fewtsnb(it(i)%map(j),fbix) )
      end if
    end do
  end do
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    i = 0
    j = 1
    k = 5
!$OMP BARRIER
!$OMP SINGLE
    thr_rhlper(j,k) = 0.0
!$OMP END SINGLE
    call thr_reduction(i,normerr,j,k)
    normerr = thr_rhlper(j,k)
!$OMP BARRIER
!$OMP SINGLE
    thr_rhlper(j,k) = 0.0
!$OMP END SINGLE
    call thr_reduction(i,normeo,j,k)
    normeo = thr_rhlper(j,k)
  end if
!$OMP BARRIER
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  write(ilog,459) normerr/(1.0-normeo)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  deallocate(pbase)
  deallocate(hit)
  deallocate(ptmp)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! this subroutine computes a non-equilibrium on a network defined as:
! mode 1:
! Tij = Z*it(i)%fewtsnb(k(j),1)*p(j)/(p(j)+p(i)), where Z is the normalizer (per i) and p is provided as input (inss)
! mode 2:
! Tij = Z*it(i)%wghtsnb(k(j),1), where Z is the normalizer
! mode 3:
! Tij = Z*it(i)%fewtsnb(k(j),1), where Z is the normalizer
! if tnode is not in the same SCC as any entry with finite sdist, a fatal exit is triggered
! the routine returns the number of steps used to achieve a probability greater than thresh in tnode and the final distribution
!
! thr_rhlper must be allocated with size (nbasins,5) at least
!
subroutine flow_MSM(it,nbasins,tnode,sdist,inss,outpx,mode,thresh,iterations,stopit,tpi)
!
  use clusters
  use iounit
#ifdef ENABLE_THREADS
  use threads
#else
  use threads, ONLY: thr_rhlper
#endif
!
  implicit none
!
  integer, INTENT(IN):: mode,nbasins,tpi,outpx,tnode,stopit
  RTYPE, INTENT(IN):: thresh,inss(nbasins),sdist(nbasins)
!
  type(t_scluster) it(nbasins)
  integer i,j,iterations,nfails,maxscc,cbnd(2),sts(4)
  logical notdone
!
  RTYPE, ALLOCATABLE:: ptmp(:)
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer OMP_GET_NUM_THREADS,k
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, flow_MSM(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  i = 1
  j = nbasins
  if (tpi.gt.0) then
    k = OMP_GET_NUM_THREADS()
    call threads_bounds(i,j,tpi,k,cbnd)
    j = 4
    call get_thread_loop_bounds(i,j,sts(1:2),sts(3:4),tpi)
    sts(2) = sts(3)
  else
    cbnd(1) = 1
    cbnd(2) = nbasins
    sts(1:2) = cbnd(1:2)
  end if
!$OMP MASTER
  do i=1,nbasins
    if ((sdist(i).gt.0.0).AND.(it(i)%inscc.ne.it(tnode)%inscc)) then
      write(ilog,*) 'Fatal. Called flow_MSM(...) with incompatible source distribution. This is a bug.'
      call fexit()
    end if
  end do
!$OMP END MASTER
!$OMP BARRIER
#else
  cbnd(1) = 1
  cbnd(2) = nbasins
  sts(1:2) = cbnd(1:2)
  do i=1,nbasins
    if ((sdist(i).gt.0.0).AND.(it(i)%inscc.ne.it(tnode)%inscc)) then
      write(ilog,*) 'Fatal. Called flow_MSM(...) with incompatible source distribution. This is a bug.'
      call fexit()
    end if
  end do
#endif
  maxscc = maxval(it(1:nbasins)%inscc)
  allocate(ptmp(nbasins))
  do i=sts(1),sts(2)
    if (mode.eq.1) then
      ptmp(i) = 0.0
      do j=1,it(i)%nb
        ptmp(i) = ptmp(i) + it(i)%fewtsnb(j,1)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if (mode.eq.2) then
      ptmp(i) = 1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))
    else if (mode.eq.3) then
      ptmp(i) = sum(it(i)%fewtsnb(1:it(i)%nb,1))
    else
      write(ilog,*) 'Called flow_MSM(...) with unsupported mode (got ',mode,'). This is a bug.'
      call fexit()
    end if
    if (ptmp(i).gt.0.0) ptmp(i) = 1.0/ptmp(i) ! normalizer
!
    if (mode.eq.1) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i)*it(i)%fewtsnb(j,1)*inss(it(i)%lstnb(j))/(inss(it(i)%lstnb(j)) + inss(i))
      end do
    else if (mode.eq.2) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i)*it(i)%wghtsnb(j,1)
      end do
    else if (mode.eq.3) then
      do j=1,it(i)%nb
        it(i)%fewtsnb(j,3) = ptmp(i)*it(i)%fewtsnb(j,1)
      end do
    end if
  end do
!
! source distribution
  thr_rhlper(cbnd(1):cbnd(2),outpx) = sdist(cbnd(1):cbnd(2))
  notdone = .true.
  iterations = 0
  nfails = 0
!
  do while ((notdone.EQV..true.).AND.(iterations.lt.stopit))
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif 
    iterations = iterations + 1
    ptmp(:) = 0.0
    do i=sts(1),sts(2)
      if (i.eq.tnode) then
        ptmp(i) = thr_rhlper(i,outpx)
        cycle
      end if
      do j=1,it(i)%nb
        if ((it(it(i)%lstnb(j))%inscc.eq.it(i)%inscc).AND.(it(i)%fewtsnb(j,1).gt.0.0)) then
          ptmp(it(i)%lstnb(j)) = ptmp(it(i)%lstnb(j)) + thr_rhlper(i,outpx)*it(i)%fewtsnb(j,3)
        end if
      end do
    end do
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then 
      j = 1
      k = 3
      i = 67
      thr_rhlper(cbnd(1):cbnd(2),5) = 0.0
      call thr_combinevec(ptmp,nbasins,j,k,i,tpi)
    else
      thr_rhlper(1:nbasins,5) = ptmp(:)
    end if
#else
    thr_rhlper(1:nbasins,5) = ptmp(:)
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    thr_rhlper(cbnd(1):cbnd(2),outpx) = thr_rhlper(cbnd(1):cbnd(2),5)
    if (thr_rhlper(tnode,5).ge.thresh) notdone = .false.
    if (notdone.EQV..false.) exit
  end do
! 
  deallocate(ptmp)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!-------------------------------------------------------------------------
!
! this subroutine computes different variants of shortest path measures for
! reaching arbitrary nodes from a source (snode)
! mode = 1: shortest path by number of edges
! mode = 2: shortest path by total geometric length of edges
! mode = 3: most probable path (minimal in negative logarithm) - this penalizes paths with low likelihood nodes on them
! mode = 4: most probable path assuming that all nodes are populated equally
! mode = 5: most probable path assuming equal capacity for all existing edges
! 
! name of the subroutine may be a misnomer
!
subroutine dijkstra_sp(it,nbasins,snode,plen,mode)
!
  use clusters
  use iounit
!
  implicit none
!
  integer nbasins,mode,snode,iterations,processed,i,j,hits,nohits,ssix
  RTYPE netlen,minlen,plen(nbasins),mtmp
  type(t_scluster) it(nbasins)
  integer, ALLOCATABLE:: touched(:)
  logical notdone
!
  write(ilog,*) 'Now computing shortest path (',mode,') from node ',snode,' ...'
!
  allocate(touched(nbasins))
!
  ssix = 1 ! only raw sampling weight at the moment (for mode 3 only)
  iterations = 0
  plen(:) = HUGE(plen)
  plen(snode) = 0.0
  if ((mode.eq.3).OR.(mode.eq.5)) plen(snode) = -log(it(snode)%nodewt(ssix))
  touched(:) = 0
  touched(snode) = 1
  notdone = .true.
  processed = nbasins - 1
  nohits = 0
  minlen = 0.0
!
  do while (notdone.EQV..true.)
!
    iterations = iterations + 1
    hits = 0
    mtmp = 0.001*HUGE(mtmp)
    do i=1,nbasins
      if (touched(i).ne.1) cycle
      do j=1,it(i)%nb
        if (it(i)%wghtsnb(j,1).le.0) cycle
        if (mode.eq.1) then
          netlen = 1.0
        else if (mode.eq.2) then
          netlen = it(i)%lensnb(j,1)
        else if (mode.eq.3) then
          netlen = -log(it(i)%nodewt(ssix)*it(i)%wghtsnb(j,1)/(1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))))
        else if (mode.eq.4) then
          netlen = -log(1.0*it(i)%wghtsnb(j,1)/(1.0*sum(it(i)%wghtsnb(1:it(i)%nb,1))))
        else if (mode.eq.5) then
          netlen = -log(it(it(i)%lstnb(j))%nodewt(ssix)/(1.0*it(i)%nb))
        end if
        mtmp = min(mtmp,netlen)
        if ((plen(i)+netlen).lt.plen(it(i)%lstnb(j))) then
          hits = hits + 1
          plen(it(i)%lstnb(j)) = plen(i) + netlen
          if (touched(it(i)%lstnb(j)).eq.0) processed = processed - 1
          touched(it(i)%lstnb(j)) = 1
        end if
      end do
      touched(i) = 2
    end do
    if (hits.eq.0) then
      nohits = nohits + 1
    else
      nohits = 0
    end if
    minlen = minlen + mtmp
    if ((maxval(plen(:)).lt.(minlen)).AND.(processed.eq.0)) exit
    if (nohits.eq.1000) exit ! emergency
  end do
!
! fix non-reachable
  netlen = 0.0
  do i=1,nbasins
    if (plen(i).eq.HUGE(plen(1))) cycle
    netlen = max(netlen,plen(i))
  end do
  j = 0
  do i=1,nbasins
    if (plen(i).eq.HUGE(plen(1))) then
      j = j + 1
      plen(i) = netlen + 1.0*j
      if ((j.eq.1).AND.(cequil.le.0)) then
        write(ilog,*) 'Warning. Path length values larger than ',netlen,' correspond to &
 &unreachable nodes (arbitrary values given).'
      end if
    end if
  end do
!
  deallocate(touched)
!
  write(ilog,*) 'Done after ',iterations,' iterations.'
!
end
!
!-------------------------------------------------------------------------
!
! iterative determination of the MFPT to snode
! 
! this routine is not thread-safe
!
! this routine relies on fewts(:,3 or 4) and edit_tmat should be called according to choice before(!) entering it
!
subroutine iterative_mfpt(it,nbasins,snode,mfptv,thresher,dofwr,maxtime)
!
  use clusters
  use iounit
  use threads, ONLY: thrdat
!
  implicit none
!
  integer, INTENT(IN):: nbasins,snode
  logical, INTENT(IN):: dofwr
  RTYPE, INTENT(IN):: thresher
  integer(KIND=8), INTENT(IN):: maxtime
!
  type(t_scluster) it(nbasins)
  RTYPE mfptv(nbasins)
!
  integer iterations,i,j,nsccs,fbix,tonext,lastt
  RTYPE normerr,kr
  logical notdone
  RTYPE, ALLOCATABLE:: incr(:)
  logical, ALLOCATABLE:: hit(:)
  integer(KIND=8) tts(2)
  character(20) timename
!
  call System_Clock(count=tts(1))
!
  if (dofwr.EQV..true.) then
    fbix = 3
    timename = 'forward '
  else
    fbix = 4
    timename = 'backward '
  end if
  call strlims(timename,i,j)
  write(ilog,'(a,i8,a,a,a)') 'Now iteratively computing mean first passage time to node ',snode,' from ',timename(i:j),&
 &' time matrix ...'
  nsccs = size(csccwts)
  allocate(incr(nbasins))
  allocate(hit(nbasins))
!
  notdone = .true.
  iterations = 0
  hit(:) = .false.
  mfptv(:) = 1.0
  mfptv(snode) = 0.0
  tonext = 10
  lastt = 0
  do while (notdone.EQV..true.)
    iterations = iterations + 1
    incr(:) = 1.0
    do i=1,nbasins
      if (i.eq.snode) cycle
      if (it(i)%inscc.ne.it(snode)%inscc) cycle
      do j=1,it(i)%nb
        if (it(it(i)%lstnb(j))%inscc.ne.it(snode)%inscc) cycle
        if (it(i)%fewtsnb(j,fbix).gt.0.) then
          incr(i) = incr(i) + it(i)%fewtsnb(j,fbix)*mfptv(it(i)%lstnb(j))
          hit(i) = .true.
        end if
      end do
    end do
    incr(snode) = 0.0
    normerr = sum(abs(incr(:)-mfptv(:)))/(1.0*nbasins)
    mfptv(:) = incr(:)
    if (normerr.lt.thresher) exit
!    if (maxval(mfptv(:)).ge.cstored) exit ! emergency exit
    if (iterations.eq.tonext) then
      call System_Clock(count=tts(2))
      if ((tts(2)-tts(1)).gt.maxtime) then
        write(ilog,*) 'Aborting after ',iterations,' iterations due to maximum time of ',(1.0*maxtime)/thrdat%rate,&
 &' [s] being exceeded.'
        exit
      else
        kr = (1.0*iterations)/((tts(2)-tts(1))/thrdat%rate) ! iterations per second
        lastt = tonext
        i = nint(kr*0.1*maxtime/thrdat%rate) ! 10% into the future of maximum allowed time
        tonext = iterations + max(1,i)
      end if
    end if
  end do
!
  kr = maxval(mfptv(:))
  j = 0
  do i=1,nbasins
    if ((hit(i).EQV..false.).AND.(i.ne.snode)) then
      kr = kr + 1.0
      mfptv(i) = kr
      j = j + 1
    end if
  end do
!
  deallocate(hit)
  deallocate(incr)
!
  write(ilog,*) 'Done after ',iterations,' iterations.'
!
end
!
!--------------------------------------------------------------------------------------------------
!
! order the network in a progress index-like manner by flux
!
subroutine network_order_by_flux(it,nbasins,snode,mfptv,dofwr,inward)
!
  use iounit
  use clusters
!
  implicit none
!
  integer, INTENT(IN):: nbasins,snode
  logical, INTENT(IN):: dofwr,inward
  type(t_scluster), INTENT(IN):: it(nbasins)
  RTYPE, INTENT(OUT):: mfptv(nbasins)
!
!  character(len=:), ALLOCATABLE:: timename
  integer i,ii,j,k,plstsz,rscc,tmpix,fbix,ssix,nclus_inrscc
  integer, ALLOCATABLE:: plst(:)
  RTYPE tmpr
  RTYPE, ALLOCATABLE:: helpvec(:)
  logical notdone
!
  if (dofwr.EQV..true.) then
    fbix = 3
    ssix = 4
!    timename = 'forward '
  else
    fbix = 4
    ssix = 5
!    timename = 'backward '
  end if
!
  rscc = it(snode)%inscc
  nclus_inrscc = 0
  do j=ssix,ssix-2,-2
    tmpr = 0.
    do i=1,nbasins
      if (it(i)%inscc.eq.rscc) then
        tmpr = tmpr + it(i)%nodewt(j)
        if (j.eq.ssix) then
          nclus_inrscc = nclus_inrscc + 1
        end if
      end if
    end do
    if (tmpr.gt.0.0) then
      exit
    end if
  end do
  ssix = max(1,j) ! loop ends on 0 or 1 if unsuccessful -> resort to sampling weight
!
  allocate(plst(nbasins))
  allocate(helpvec(nbasins))
  plst(1) = snode
  mfptv(:) = 0
  mfptv(snode) = 1.
  plstsz = 1
  notdone = .true.
  helpvec(:) = 0
  if (inward.EQV..false.) then ! find the largest increment from set A to B under ss-conditions
    do while (plstsz.lt.nclus_inrscc)
      tmpix = 0
      i = plstsz                        !consider only the last added one as the ref. val. for the others is already there
      ii = plst(i)                      !cluster in set A in position i (last added)
      do j=1,it(ii)%nb
        k = it(ii)%lstnb(j)
        if (it(k)%inscc.ne.rscc) cycle  !not in the same component
        if (mfptv(k).gt.0) cycle        !already added in set A
        helpvec(k) = helpvec(k) + it(ii)%nodewt(ssix)*it(ii)%fewtsnb(j,fbix)  !cumulating
      end do
      tmpix = maxloc(helpvec,dim=1)  !if I do not put dim=1, then tmpix has to be an array of rank 1
      plstsz = plstsz + 1
      mfptv(tmpix) = 1.0*plstsz
      plst(plstsz) = tmpix
      helpvec(tmpix) = 0. !the just added guy won't be detected by maxloc anymore, while the others will keep on growing
    end do
  else ! here we find the largest incoming increment into set A from B under ss-conditions
    do while (plstsz.lt.nclus_inrscc)
      tmpix = 0
      i = plstsz
      ii = plst(i)
      do j=1,it(ii)%nb
        k = it(ii)%lstnb(j)
        if (it(k)%inscc.ne.rscc) cycle
        if (mfptv(k).gt.0) cycle
        !helpvec(k) = helpvec(k) + it(k)%nodewt(ssix)*it(k)%fewtsnb(it(ii)%map(j),fbix) !candidate std state, singles disadvantaged
        helpvec(k) = helpvec(k) + it(ii)%nodewt(ssix)*it(k)%fewtsnb(it(ii)%map(j),fbix) !parent std state, singles advantaged
      end do
      tmpix = maxloc(helpvec,dim=1)
      plstsz = plstsz + 1
      mfptv(tmpix) = 1.0*plstsz
      plst(plstsz) = tmpix
      helpvec(tmpix) = 0.
    end do
  end if
  deallocate(plst)
  deallocate(helpvec)
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine min_st_cut(it,nbasins,snode,tnode,cutmode,cutval,spart,spartsz)
!
  use clusters
  use iounit
!
  implicit none
!
  integer ii,i,j,k,tnodei
  integer nbasins,snode,tnode,cutval,spart(nbasins),spartsz
  integer startnode,compcut,setincr,oldsetsz,curnode,getexflow,setsz
  integer, ALLOCATABLE:: setmap(:)
  logical cutmode,notdone
  type(t_scluster) it(nbasins)
!
  if (snode.eq.tnode) then
    write(ilog,*) 'Fatal. Called min_st_cut(...) with identical source and target nodes. This is a bug.'
    call fexit()
  end if
!
  do i=1,nbasins
    it(i)%ldis = 0
    it(i)%active = 0
    it(i)%flwnb(:,:) = 0
    it(i)%rflw = 0
  end do
!
  
!
  startnode = -1
!  write(*,*) snode,tnode,it(snode)%nb,it(tnode)%nb
!  write(*,*) 'SL'
!  write(*,*) it(snode)%lstnb(1:it(snode)%nb)
!   write(*,*) 'SW'
!  write(*,*) it(snode)%wghtsnb(1:it(snode)%nb,:)
!  write(*,*) 'TL'
!  write(*,*) it(tnode)%lstnb(1:it(tnode)%nb)
!  write(*,*) 'TW'
!  write(*,*) it(tnode)%wghtsnb(1:it(tnode)%nb,:)


  do i=1,it(snode)%nb
    if (it(snode)%lstnb(i).eq.snode) cycle
    j = it(snode)%lstnb(i)
    if (j.eq.tnode) tnodei = i
    it(snode)%flwnb(i,1) = it(snode)%flwnb(i,1) + it(snode)%wghtsnb(i,1)
    k = it(snode)%map(i) ! k = it(j)%map(snode)
    it(j)%flwnb(k,2) = it(j)%flwnb(k,2) + it(snode)%wghtsnb(i,1)
    if ((it(snode)%wghtsnb(i,1).gt.0).AND.(j.ne.tnode)) then
      startnode = j
      it(j)%active = 1
    end if
  end do
  it(snode)%ldis = nbasins
!
  if (startnode.eq.-1) then
    cutval = it(snode)%wghtsnb(it(tnode)%map(tnodei),1) ! it(snode)%wghtsnb(it(snode)%map(tnode),1)
!   wait
  else
    notdone = .true.
    do while (notdone.EQV..true.) 
      call process_pushpreflow(it,startnode,snode,tnode,cutmode,curnode,nbasins)
      if (curnode.eq.-1) then
        notdone = .false.
        do i=1,nbasins
          if ((i.ne.snode).AND.(i.ne.tnode).AND.(it(i)%active.gt.0)) then
            startnode = i
            notdone = .true.
          end if
        end do
      else
        startnode = curnode
      end if
    end do
  end if
!
  setincr = 1
  it(tnode)%rflw = 1
  setsz = 1
  allocate(setmap(nbasins))
  setmap(:) = 0
  setmap(tnode) = 1
  do while (setincr.gt.0) 
    oldsetsz = setsz
    do i=1,nbasins
      if (it(i)%rflw.eq.0) cycle
      call process_pushresflow(it,i,nbasins,setmap,setsz)
    end do
    setincr = setsz - oldsetsz
  end do
  spartsz = 0
  do i=1,nbasins
    if (setmap(i).eq.0) then
      spartsz = spartsz + 1
      spart(spartsz) = i
    end if
  end do
!  write(*,*) 'SPART: ',spart(1:spartsz)
!
  cutval = 0
  do ii=1,nbasins
    do i=1,it(ii)%nb
      if (it(ii)%lstnb(i).eq.ii) cycle
      j = it(ii)%lstnb(i)
      if (j.lt.ii) then
        if ((setmap(ii).eq.0).AND.(setmap(j).eq.1)) then
          cutval = cutval + it(ii)%wghtsnb(i,1)
        else if ((setmap(j).eq.0).AND.(setmap(ii).eq.1)) then
          cutval = cutval + it(ii)%wghtsnb(i,2)
        end if
      end if
    end do
  end do
  compcut = getexflow(it,tnode)
  if (cutval.ne.compcut) then
    write(ilog,*) 'Fatal. Push-relabel algorithm failure in min_st_cut(...). &
 &Sum of cuts is ',cutval,'  but excess flow at target is ',compcut,'. This is a bug.'
    call fexit()
  end if
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine process_pushpreflow(it,inode,snode,tnode,cutmode,curnode,nbasins)
!  
  use clusters
  use iounit
!
  implicit none
!
  integer inode,snode,tnode,curnode,nbasins
  integer i,j,k,minldis,maxfl,avail,exflow
  logical cutmode
  integer block(nbasins)
  type(t_scluster) it(nbasins)
!
  curnode = -1
  exflow = 0
  avail = it(inode)%nb
!
  do i=1,it(inode)%nb
    if (it(inode)%lstnb(i).eq.inode) then
      block(i) = 0
      avail = avail - 1
      cycle
    end if
    j = it(inode)%lstnb(i) 
    block(i) = 0
    exflow = exflow + it(inode)%flwnb(i,2) - it(inode)%flwnb(i,1)
    if (it(j)%ldis.ne.(it(inode)%ldis-1)) then
      block(i) = 1
      avail = avail - 1
    end if
  end do
!  write(*,*) 'Excess flow for ',inode,' is ',exflow
!
  do while (exflow.ne.0)
    if (avail.gt.0) then
      do i=1,it(inode)%nb
        if (it(inode)%lstnb(i).eq.inode) cycle
        j = it(inode)%lstnb(i)
        k = it(inode)%map(i) ! it(j)%map(inode)
        if (block(i).eq.0) then
          maxfl = it(inode)%wghtsnb(i,1) - it(inode)%flwnb(i,1) + it(j)%flwnb(k,1)
!          write(*,*) 'Pushing ',maxfl,' / ',exflow,' down ',i,' ( ',j,').'
          if ((curnode.eq.-1).AND.(j.ne.tnode).AND.(j.ne.snode)) curnode = j
          if (j.ne.tnode) then
            if (cutmode.EQV..true.) then
              if (it(j)%ldis.lt.nbasins) it(j)%active = 1
            else
              it(j)%active = 1
            end if
          end if
          if (maxfl.ge.exflow) then
            it(j)%flwnb(k,2) = it(j)%flwnb(k,2) + exflow
            it(inode)%flwnb(i,1) = it(inode)%flwnb(i,1) + exflow
            exflow = 0
          else
            it(j)%flwnb(k,2) = it(j)%flwnb(k,2) + maxfl
            it(inode)%flwnb(i,1) = it(inode)%flwnb(i,1) + maxfl
            exflow = exflow - maxfl
            block(i) = 1
            avail = avail - 1
          end if
          if (exflow.eq.0) exit
        end if
      end do
    else
      if (exflow.eq.0) exit
      minldis = HUGE(minldis)
      do i=1,it(inode)%nb
        if (it(inode)%lstnb(i).eq.inode) cycle
        j = it(inode)%lstnb(i) 
        k = it(inode)%map(i) ! it(j)%map(inode)
        maxfl = it(inode)%wghtsnb(i,1) - it(inode)%flwnb(i,1) + it(j)%flwnb(k,1)
        if (maxfl.gt.0) then
!          write(*,*) 'eligible: ',j, '(',inode,').'
          if (it(j)%ldis.lt.minldis) minldis = it(j)%ldis
        end if
      end do
!       write(*,*) 'Relabeling ',inode,' to ',minldis + 1,'.'
      it(inode)%ldis = minldis + 1
      do i=1,it(inode)%nb
        if (it(inode)%lstnb(i).eq.inode) cycle
        j = it(inode)%lstnb(i)
        k = it(inode)%map(i) ! it(j)%map(inode)
        maxfl = it(inode)%wghtsnb(i,1) - it(inode)%flwnb(i,1) + it(j)%flwnb(k,1)
!        write(*,*) 'Max flow for ',j,' is ',maxfl,'.'
        if ((maxfl.gt.0).AND.(it(j)%ldis.eq.(it(inode)%ldis-1))) then
          block(i) = 0
          avail =  avail + 1
        end if
      end do
    end if
  end do
!
  it(inode)%active = 0
!
end
!
!-------------------------------------------------------------------------------------
!
subroutine process_pushresflow(it,inode,nbasins,setmap,setsz)
!  
  use clusters
  use iounit
!
  implicit none
!
  integer i,j,inode,nbasins,setmap(nbasins),setsz
  type(t_scluster) it(nbasins)
!
  do i=1,it(inode)%nb
    if (it(inode)%lstnb(i).eq.inode) cycle
    j = it(inode)%lstnb(i)
    if (it(inode)%wghtsnb(i,2).gt.(it(inode)%flwnb(i,2)-it(inode)%flwnb(i,1))) then
      it(j)%rflw = 1
      if (setmap(j).eq.0) then
        setsz = setsz + 1
        setmap(j) = 1
      end if
    end if
  end do
!
end
!
!--------------------------------------------------------------------------------------------
!
function getexflow(it,inode)
!
  use clusters
!
  implicit none
!
  integer getexflow,inode,i
  type(t_scluster) it(*)
!
  getexflow = 0
  do i=1,it(inode)%nb
    if (it(inode)%lstnb(i).eq.inode) cycle
    getexflow = getexflow + it(inode)%flwnb(i,2) - it(inode)%flwnb(i,1)
  end do
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine tarjan_SCC(it,nbasins,mode)
!
  use clusters
!
  implicit none
!
  integer, INTENT(IN):: mode,nbasins
!
  integer cursccid,glidx,i,idxcnt,mcopy
  integer, ALLOCATABLE:: idxlst(:)
  type (t_scluster) it(nbasins)
  RTYPE norm
!
  mcopy = mode
  allocate(idxlst(nbasins+10))
!
  it(1:nbasins)%active = 0
  it(1:nbasins)%rflw = -1
  glidx = 0
  cursccid = 1
  it(1:nbasins)%inscc = 0
  idxcnt = 0
!
  do i=1,nbasins
    if (it(i)%inscc.gt.0) cycle
    call graph_dfs(it,nbasins,i,glidx,cursccid,idxlst,idxcnt,mcopy)
  end do
!  
  cursccid = maxval(it(1:nbasins)%inscc)
  if (allocated(csccwts).EQV..true.) deallocate(csccwts)
  allocate(csccwts(cursccid))
  csccwts(:) = 0.
  idxlst(:) = 0
  do i=1,nbasins
    idxlst(it(i)%inscc) = idxlst(it(i)%inscc) + 1
    csccwts(it(i)%inscc) = csccwts(it(i)%inscc) + it(i)%nodewt(1)
    it(i)%ix(2) = idxlst(it(i)%inscc) ! cluster numbering relative to component
  end do
  norm = 1.0/sum(csccwts)
  csccwts(:) = norm*csccwts(:)
!
  deallocate(idxlst)
!
end
!
!-------------------------------------------------------------------------------------------
!
recursive subroutine graph_dfs(it,nbasins,which,glidx,cursccid,idxlst,idxcnt,mcopy)
!
  use clusters
!
  implicit none
!
  integer glidx,i,nbasins,cursccid,which,idxcnt,idxlst(nbasins+10),mcopy
  type (t_scluster) it(nbasins)
!
  it(which)%ldis = glidx
  it(which)%active = glidx
  it(which)%rflw = 1
  glidx = glidx + 1
  idxcnt = idxcnt + 1
  idxlst(idxcnt) = which
  it(which)%inscc = cursccid
!
!  if (it(it(which)%lstnb(1))%nb.eq.1) then
!    write(*,*) it(which)%lstnb(1),' is a leaf'
!  else if ((it(it(which)%lstnb(1))%nb.eq.2).AND.&
! &    ((it(it(which)%lstnb(1))%wghtsnb(1,1).eq.0).AND.(it(it(which)%lstnb(1))%wghtsnb(2,2).eq.0)).OR.&
! &    ((it(it(which)%lstnb(1))%wghtsnb(1,2).eq.0).AND.(it(it(which)%lstnb(1))%wghtsnb(2,1).eq.0))) then
!    write(*,*) it(which)%lstnb(1),' is on a loop'
!  end if
  do i=1,it(which)%nb
    if (it(which)%lstnb(i).eq.which) cycle
    if (((mcopy.eq.1).AND.(it(which)%wghtsnb(i,1).gt.0)).OR.&
 &      ((mcopy.eq.2).AND.(it(which)%fewtsnb(i,1).gt.0.))) then
      if (it(it(which)%lstnb(i))%active.le.0) then
        call graph_dfs(it,nbasins,it(which)%lstnb(i),glidx,cursccid,idxlst,idxcnt,mcopy)
        it(which)%ldis = min(it(which)%ldis,it(it(which)%lstnb(i))%ldis)
      else if (it(it(which)%lstnb(i))%rflw.eq.1) then
        it(which)%ldis = min(it(which)%ldis,it(it(which)%lstnb(i))%active)
      end if
    end if
  end do
!
  if (it(which)%ldis.eq.it(which)%active) then
    do i=idxcnt,1,-1
      it(idxlst(i))%inscc = cursccid
      it(idxlst(i))%rflw = -1
      if (idxlst(i).eq.which) then
        idxcnt = i-1
        exit
      end if
    end do
    cursccid = cursccid + 1
  end if
!
end
!
!--------------------------------------------------------------------------------
!
! This subroutine generates synthetic trajectories through random walks on a complex network
! It is assumed that it(:)%fewtsnb(:,3:4) are set to the cumulative sums already (see edit_tmat(...))
!
! iv1: cluster assignment per snapshot
! inissnap, endssnap: snaposhots to read start (modes 1-2) and target (mode 1) clusters from
! dofwr: whether we are using forward or backward time matrix
!
! for threads: each random walk has 100% data dependency on the prior state (Markovian), i.e., 
! parallelization is only across walks and not within 
!
subroutine synth_traj_MSM(it,nbasins,iv1,inissnap,endssnap,dofwr,tpi)
!
  use clusters
  use iounit
  use interfaces
  use threads, ONLY: thr_hlper
!
  implicit none
!
  integer, INTENT(IN):: nbasins,iv1(cstored),inissnap,endssnap,tpi
  logical, INTENT(IN):: dofwr
!
  type(t_scluster) it(nbasins)
!
  integer i,curclu,endclu,istraj,iiter,nextclu,fbix,ssix!,nextsnap
  integer tname,lastini,iu,iline,tl,k,iniclu_s12,stao(1:2),rscc,t1,t2
  integer(KIND=8) tlens(2)
  RTYPE randnum,hlp1(maxval(it(1:nbasins)%nb))
  integer, ALLOCATABLE:: straj(:)  !synthetic trajectory
  integer, ALLOCATABLE:: tstats(:,:)
  RTYPE, ALLOCATABLE:: ptmp(:),ptmp2(:)
  character(MAXSTRLEN) filename,timename
  character(5) nod
  character(3) ttype
  logical notdone,norm !norm: if .true., the cumulative sum have to end in 1 (avoid numerical problems).
#ifdef ENABLE_THREADS
  RTYPE, ALLOCATABLE:: rdbuf(:)
  logical OMP_IN_PARALLEL
  integer OMP_GET_NUM_THREADS,tpn,rdbufcnt,rdbufsize,rdbufthresh,kix,kk
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, synth_traj_MSM(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  if (tpi.gt.0) then
    tpn = OMP_GET_NUM_THREADS()
    i = 1
    k = nstraj
    call threads_bounds(i,k,tpi,tpn,stao(1:2))
  else
    tpn = 1
    stao(1) = 1
    stao(2) = nstraj
  end if
  rdbufsize = nssnap + 1
  rdbufthresh = nssnap + 1
  allocate(rdbuf(rdbufsize))
  rdbufcnt = rdbufthresh
#else
  RTYPE random
  integer freeunit
!
  stao(1) = 1
  stao(2) = nstraj
#endif
!
#ifdef ENABLE_THREADS
!$OMP MASTER
  allocate(thr_hlper(tpn,1))
  call get_freeunits(tpn,thr_hlper(:,1))
!$OMP END MASTER
!$OMP BARRIER
  iu = thr_hlper(tpi,1)
#else
  iu = freeunit()
#endif
! initializations
  rscc = it(iv1(inissnap))%inscc
  istraj = 0
  norm = .true. !all cumulative sums have to sum up to 1 in this subroutine
  if (dofwr.EQV..true.) then
    ttype = 'FWD'
    fbix = 3
    timename = 'forward '
  else
    ttype = 'BWD'
    fbix = 4
    timename = 'backward '
  end if 
!
! allocate
  allocate(ptmp(nbasins))
  allocate(ptmp2(nbasins))
  allocate(straj(nssnap))
  allocate(tstats(nstraj,3))
  tstats(:,:) = 0
!
! Where appropriate, define initial and target (end) clusters
  if ((synmode.eq.1).or.(synmode.eq.2)) then !Ini clus is defined in mode 3 only to specify target component.
    iniclu_s12 = iv1(inissnap)
    if ((iniclu_s12.le.0).OR.(iniclu_s12.gt.nbasins)) then
      write(ilog,*) 'Fatal. Unable to identify initial cluster from initial snapshot (out of range). This is a bug.'
      call fexit()
    end if
    if (tpi.le.1) write(ilog,*) 'Mapping the initial snapshot to cluster ',iniclu_s12,'.'
  else if (synmode.eq.3) then
    ssix = 1
    do i=1,nbasins
      if (it(i)%inscc.eq.rscc) then
        ptmp(i) = it(i)%nodewt(ssix)/csccwts(rscc)
      else
        ptmp(i) = 0.
      end if
    end do
    k = nbasins
    call cumsum(ptmp(1:k),k,norm,ptmp2(1:k)) ! ptmp2 for reuse later -> pick via binary search only from rscc
  end if
!
! few checks
  if (synmode.eq.1) then                     !Target (end cluster) is defined only in mode 1
    endclu = iv1(endssnap)                   !endssnap defined in absolute terms
    if (tpi.le.1)  write(ilog,*) 'Mapping the target snapshot to cluster ',endclu,'.'
    if (iv1(inissnap).eq.iv1(endssnap)) then
      write(ilog,*)
      if (tpi.le.1) write(ilog,*) 'Warning. Initial and target snapshots for the generation of synthetic &
 &trajectories belong to the same cluster. Disabling.'
      write(ilog,*)
      synmode = 0
      return 
    end if
    if ((endclu.le.0).OR.(endclu.gt.nbasins)) then
      write(ilog,*) 'Fatal. Unable to identify target cluster from target snapshot (out of range). This is a bug.'
      call fexit()
    end if
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  tname = 0
  tlens(:) = 0
! Generation of synthetic trajectories
  do istraj=stao(1),stao(2)
!   replenish random numbers if needed
#ifdef ENABLE_THREADS
    if (((rdbufcnt+nssnap+2).gt.rdbufthresh).AND.(rdbufcnt.le.rdbufthresh)) then
      kix = rdbufthresh+1
      if (kix.gt.rdbufsize) kix = 1
      kk = nssnap + 1
!$OMP CRITICAL(RDBUF_REPL)
      call get_nrandoms(kk,rdbuf(kix:(kix+nssnap)))
!$OMP END CRITICAL(RDBUF_REPL)
      rdbufthresh = rdbufthresh + nssnap + 1
      if (rdbufthresh.gt.rdbufsize) rdbufthresh = nssnap + 1
    end if
#endif
    if (synmode.eq.3) then
#ifdef ENABLE_THREADS
      rdbufcnt = rdbufcnt + 1
      if (rdbufcnt.gt.rdbufsize) rdbufcnt = 1
      randnum = rdbuf(rdbufcnt)
#else
      randnum = random()
#endif
      k = nbasins
      call binary_search(k,ptmp2(1:k),randnum,i) !select a random starting cluster in rscc
      i = i + 1  !the return of binary search starts from 0 here as my input vector starts at finite weight value
      curclu = i
      if (it(curclu)%inscc.ne.rscc) then !This check explicitly refers to the whole data set
        write(ilog,*) 'Fatal. Selected initial cluster does not belong to the chosen reference component. This is a bug.'
        call fexit()
      end if
    else !establish the initial cluster for synmode 1 or 2.
      curclu = iniclu_s12
    end if
    straj(:) = 0
    iiter = 1
    straj(iiter) = curclu
    notdone = .true.
    lastini = 1
    do iiter=2,nssnap
#ifdef ENABLE_THREADS
      rdbufcnt = rdbufcnt + 1
      if (rdbufcnt.gt.rdbufsize) rdbufcnt = 1
      randnum = rdbuf(rdbufcnt)
#else
      randnum = random()
#endif
      hlp1(1:it(curclu)%nb) = it(curclu)%fewtsnb(1:it(curclu)%nb,fbix)
      call binary_search(it(curclu)%nb,hlp1(1:it(curclu)%nb),randnum,i) !where the walker is going next
      i = i + 1                         !the return of binary search starts from 0 as above 
      nextclu = it(curclu)%lstnb(i) !enumeration on the ref. comp.
      if (nextclu.eq.iniclu_s12) lastini = iiter
      if (rscc.ne.it(nextclu)%inscc) then !This check explicitly refers to the whole data set
        write(ilog,*) 'Fatal. Leaving the reference strongly connected component in synth_traj_MSM(...). This is a bug.'
        call fexit()
      end if
!     for picking snapshots from clusters instead
!      randnum = random()
!      nextsnap = tmat(curclu)%snaps(ceiling(tmat(curclu)%nmbrs*randnum)) ! it(curclu)%center
      straj(iiter) = nextclu
      if ((synmode.eq.1).AND.(nextclu.eq.endclu)) then
        tlens(1) = tlens(1) + iiter
        tlens(2) = tlens(2) + iiter - lastini
        tname = tname + 1
        exit
      end if
      if (iiter.eq.nssnap) exit
      curclu = nextclu
    end do
!
    tstats(istraj-stao(1)+1,1) = straj(1)
    tstats(istraj-stao(1)+1,2) = straj(iiter)
    tstats(istraj-stao(1)+1,3) = iiter
!
    if ((nsskip.gt.0).AND.(((synmode.eq.1).AND.(straj(iiter).eq.endclu)).OR.(synmode.ne.1))) then ! write 
      tl =  5
      call int2str(istraj,nod,tl)
      filename(:) = ' '
      filename = 'MSM_SYN_TRAJ_'//nod(1:tl)//'_'//ttype//'.frames'
      open(unit=iu,file=trim(filename),status='unknown')
      write(iu,'("# from ",i9," to ",i9," in ",i10," steps")') (tstats(istraj-stao(1)+1,i),i=1,3)
      do iline=nsskip,iiter,nsskip
        write(iu,'(I9)') straj(iline)
      end do
      close(unit=iu)
    end if
  end do
!
  if (synmode.eq.1) then
    call strlims(timename,t1,t2)
 41 format(' In ',i9,' attempts, thread ',i5,' produced no reactive trajectories (random walks) from ',a,' time matrix.')
 42 format(' In ',i9,' attempts, produced no reactive trajectories (random walks) from ',a,' time matrix.')
 43 format(' In ',i9,' attempts, thread ',i5,' produced ',i9,' reactive trajectories (random walks) with an average length of '&
 &,g13.6,' steps (',g13.6,' steps since last visit to starting point) from ',a,' time matrix.')
 44 format(' In ',i9,' attempts, produced ',i9,' reactive trajectories (random walks) with an average length of ',g13.6,' steps&
 & (',g13.6,' steps since last visit to starting point) from ',a,' time matrix.')
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SYNTH_TRAJ_REPORT)
    if (tname.gt.0) then
      write(ilog,43) stao(2)-stao(1)+1,tpi,tname,(1.0*tlens(1:2))/(1.0*tname),timename(t1:t2)
    else if (stao(2).ge.stao(1)) then
      write(ilog,41) stao(2)-stao(1)+1,tpi,timename(t1:t2)
    end if
!$OMP END CRITICAL(SYNTH_TRAJ_REPORT)
  deallocate(rdbuf)
#else
    if (tname.gt.0) then
      write(ilog,44) stao(2)-stao(1)+1,tname,(1.0*tlens(1:2))/(1.0*tname),timename(t1:t2)
    else if (stao(2).ge.stao(1)) then
      write(ilog,42) stao(2)-stao(1)+1,timename(t1:t2)
    end if
#endif
  end if
  deallocate(ptmp2)
  deallocate(ptmp)
  deallocate(straj)
  deallocate(tstats)
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
  deallocate(thr_hlper)
!$OMP END SINGLE
#endif
!
end subroutine synth_traj_MSM
!
!--------------------------------------------------------------------------------
!
! a subroutine to select eligible snapshots/clusters for graph operations
! requiring start/end (F/U) states
!
  subroutine sel_ssnap(it,nbasins,ivx,ivxsz,inissnap,endssnap,tpi)
!
  use clusters
  use iounit
  use pdb, ONLY: select_frames
!
  implicit none
!
  integer, INTENT(IN)::  nbasins,ivxsz,ivx(ivxsz,5),tpi
  integer, INTENT(OUT):: inissnap,endssnap
!
  type(t_scluster) it(nbasins)
!
  integer j,k
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, sel_ssnap(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  inissnap = 0
  endssnap = 0
  if ((dopfold.EQV..false.).AND.(synmode.eq.0).AND.(eigvalmd.eq.0).AND.(ccfepmode.eq.0)) return
!
  k = cstored
!
!  what if 0: set to centroid of largest cluster (always true regardless settings for framesfile, ccollect, ...)
  if ((inissnap_usrslct.eq.0).OR.(inissnap_usrslct.eq.-1)) then   !default
    inissnap = it(1)%center         !centroid of the largest cluster
    do j=2,nbasins
      if (it(j)%nmbrs.gt.it(1)%nmbrs) then
        write(ilog,*) 'Fatal. Clusters must be sorted by decreasing size in sel_ssnap(...). This is a bug.'
        call fexit()
      end if
    end do
  else ! otherwise, we have already remapped the input to 1:cstored (now in inissnap_usrslct)
    inissnap = inissnap_usrslct
  end if
!
  if (dopfold.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if (allocated(clufold).EQV..false.) then
      if (inissnap.gt.0) then
        allocate(clufold(1))
        clufold(1) = inissnap
      else
        write(ilog,*) 'Fatal. Received an illegal value for the initial snapshot from FMCSC_CLUFOLDFILE assignment &
 &in sel_ssnap(...). This is a bug.'
        call fexit()
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
!
  if ((synmode.eq.1).OR.(dopfold.EQV..true.)) then  !cases where end snapshot is important
    if (endssnap_usrslct.eq.0) then   !default
      endssnap = cstored              !last stored snapshot
    else ! otherwise, we have already remapped the input to 1:cstored (now in inissnap_usrslct)
      endssnap = endssnap_usrslct
    end if
!
    if (dopfold.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      if (allocated(cluunfold).EQV..false.) then
        if (endssnap.gt.0) then
          allocate(cluunfold(1))
          cluunfold(1) = endssnap
        else
          write(ilog,*) 'Fatal. Received an illegal value for the target snapshot from FMCSC_CLUUNFOLDFILE assignment &
 &in sel_ssnap(...). This is a bug.'
          call fexit()
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
    end if
  end if
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
! print final selection (this may be redundant)
  refscc = it(ivx(inissnap,4))%inscc 
  write(ilog,*) 
  write(ilog,'(a)') '------------------------------------------------------------------------------------' 
  write(ilog,'(A,I9,A,I9,A)') 'Final result from FMCSC_INISYNSNAP selection: ',inissnap,' (component: ',refscc,')'
  if ((synmode.eq.1).OR.(dopfold.EQV..true.)) then
    write(ilog,'(A,I9,A,I9,A)') 'Final result from FMCSC_ENDSYNSNAP selection: ',endssnap,' (component: ',&
 &it(ivx(endssnap,4))%inscc,')'
  end if
  write(ilog,'(a)') '------------------------------------------------------------------------------------'
  write(ilog,*) 
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
end subroutine sel_ssnap
!
!-------------------------------------------------------------------------------------------
!
subroutine mergetrees(ntrees,a,b,ib,snap2tree,lsnap2tree)
!
  use iounit
  use clusters
!
  implicit none
!
  integer ntrees
  type(t_progindextree) a,b ! component a is added to component b
  integer ib ! index of component b
  integer lsnap2tree
  integer snap2tree(lsnap2tree)
  integer, ALLOCATABLE:: snaps2(:)
  integer i ! local snapshot index, loop variable
  integer nsnaps2 ! new number of snapshots
!
  nsnaps2 = b%nsnaps + a%nsnaps
  allocate(snaps2(nsnaps2))
!
  ntrees = ntrees - 1
! update information for component b
  do i=1,b%nsnaps
    snaps2(i) = b%snaps(i)
  end do
  do i=(b%nsnaps+1),nsnaps2
    snaps2(i) = a%snaps(i)
  end do
  b%nsnaps = nsnaps2
  deallocate(b%snaps)
  allocate(b%snaps(nsnaps2))
  do i=1,nsnaps2
    b%snaps(i) = snaps2(i)
  end do
! UPDATE SNAPS2COMP
  do i=1,a%nsnaps
    snap2tree(a%snaps(i)) = ib
  end do
! indicate component a as empty
  a%nsnaps = -1
  deallocate(a%snaps)
end
!
!-----------------------------------------------------------------------------------
!
subroutine heapify(heap,n,key,source,a)
!
  implicit none
!
  integer n,a,i,x,hleft,hright
  integer heap(n),source(n)
  real(KIND=4) key(n)
  logical isheap
!
  isheap = .false.
  i = a
  do while(isheap.EQV..false.)
    x = i
    if (hleft(i).le.n) then
      if (key(hleft(i)).lt.key(x)) then
        x = hleft(i)
      end if
    end if
    if (hright(i).le.n) then
      if (key(hright(i)).lt.key(x)) then
        x = hright(i)
      end if
    end if
    if (x.eq.i) then
      isheap = .true.
    else
      call hswap(heap,n,key,source,i,x)
      i = x
    end if
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hswap(heap,n,key,source,i,x)
!
  implicit none
!
  integer n,i,x
  integer heap(n),source(n)
  real(KIND=4) key(n)
  integer temp,tempsource
  real(KIND=4) tempkey
!
  temp = heap(i)
  tempkey = key(i)
  tempsource = source(i)
  heap(i) = heap(x)
  source(i) = source(x)
  key(i) = key(x)
  heap(x) = temp
  key(x) = tempkey
  source(x) = tempsource
end
!
!-----------------------------------------------------------------------------------
!
subroutine hbuild(heap,n,key,source)
!
  implicit none
!
  integer n,hparent,a
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  do a=hparent(n),1,-1
    call heapify(heap,n,key,source,a)
  end do
end
!
!-----------------------------------------------------------------------------------
!
function hparent(i)
!
  implicit none
!
  integer hparent,i
!
  if (i.eq.1) then
!    write(ilog,*) 'NICO: Bug! Parent of root of heap requested.'
!    call fexit()
  end if
  hparent = i/2
!
end
!
!-----------------------------------------------------------------------------------
!
function hleft(i)
!
  implicit none
!
  integer hleft,i
!
  hleft = 2*i
end
!
!-----------------------------------------------------------------------------------
!
function hright(i)
!
  implicit none
!
  integer hright,i
!
  hright = 2*i+1
end
!
!-----------------------------------------------------------------------------------
!
subroutine hprint(heap,n,key,source)
!
  use iounit
!
  implicit none
!
  integer n,i,hparent,hleft,hright
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  do i=1,n
    write(ilog,*) heap(i),key(i),source(i),hparent(i),hleft(i),hright(i)
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hinsert(heap,n,key,source,new,newkey,newsource)
!
  use clusters
  use iounit
!
  implicit none
!
  integer n,new,hparent,i,newsource
  integer heap(n+1),source(n+1)
  real(KIND=4) key(n+1)
  real(KIND=4) newkey
!
  if (n.ge.cstored) then
    write(ilog,*) 'NICO: Bug! Max heap size reached.'
    call fexit()
  end if
  n = n+1
  heap(n) = new
  key(n) = newkey
  source(n) = newsource
  i = n
  do while (i.gt.1)
    if (key(i).lt.key(hparent(i))) then
      call hswap(heap,n,key,source,i,hparent(i))
      i = hparent(i)
    else
      exit
    end if
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hremovemin(heap,n,key,source)
!
  implicit none
!
  integer n
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  heap(1) = heap(n)
  key(1) = key(n)
  source(1) = source(n)
  n = n-1
  call heapify(heap,n,key,source,1)
end
!
!-----------------------------------------------------------------------------------
!
subroutine pidxtree_growsiblings(t1)
!
  use clusters
!
  type(t_progindextree) t1
  integer tmp(max(t1%nsibalsz,1)),k
!
  if (t1%nsibalsz.eq.0) then
    t1%nsibalsz = 5
    allocate(t1%siblings(t1%nsibalsz))
  else
    tmp(:) = t1%siblings(1:t1%nsibalsz)
    deallocate(t1%siblings)
    k = t1%nsibalsz
    t1%nsibalsz = t1%nsibalsz*2
    allocate(t1%siblings(t1%nsibalsz))
    t1%siblings(1:k) = tmp(:)
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
! in the PIGS trace, the reseeding step is masked (zero increment)
! the trajectories contain the information prior to reseeding, and the trace provides posterior mappings
!
! the subroutine rewrites the (simple) input map (tlst) to reroute transitions (snapshot level)
! according to what happened in PIGS
!
! this is an extremely annoying algorithm to do anything with due to there being so many different
! parameters to keep track; it is currently assumed general, however
!
subroutine PIGS_ReadTrace(tlst,tlstsz)
!
  use clusters
  use mpistuff
  use iounit
  use interfaces
  use system, ONLY: nsim
  use pdb, ONLY: select_frames
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer, INTENT(IN):: tlstsz
!
  integer tlst(tlstsz,5) ! a preassembled list of all default i:i+1 transition pairs: to be altered here columns 1 and 5
!
  integer istep,kstep,lstep,sstep,xstep,iomess,spert,i,j,jj,k,kk,kkk,ii,which,klines,oldwhich,recwhich,minstep
  integer iu,freeunit,t1,t2,lmap(re_conditions),history
  integer, ALLOCATABLE:: lmapmat(:,:)
  logical exists,allout
!  RTYPE dis
!
  call strlims(re_traceinfile,t1,t2)
  inquire(file=re_traceinfile(t1:t2),exist=exists)
  if (exists.EQV..true.) then
    iu = freeunit()
    open(unit=iu,file=re_traceinfile(t1:t2),status='old',action='read')
  else
    return
  end if
!
  if ((re_conditions.le.0).OR.(re_aux(4).le.0).OR.(re_aux(8).le.0)) then
    write(ilog,*) 'Fatal. The processing of a PIGS trace in clustering graph-based analyses &
 &requires setting FMCSC_RE_TRAJOUT, FMCSC_RE_TRAJSKIP, FMCSC_RE_TRAJTOTAL, and FMCSC_REPLICAS to appropriate values.'
    call fexit()
!  else if (nequil.gt.0) then
!    write(ilog,*) 'Fatal. The processing of a PIGS trace in clustering graph-based analyses is incompatible with setting &
! &FMCSC_EQUIL to something other than 0.'
!    call fexit()
  else if (select_frames.EQV..true.) then
!    if (pdb_
!    if (re_aux(8)*re_conditions.gt.abs(curframe)) then ! curframe holds original nsim as backup
!      write(ilog,*) 'Fatal. The processing of a PIGS trace in clustering graph-based analyses&
! & requires setting FMCSC_RE_TRAJTOTAL, FMCSC_REPLICAS, and FMCSC_NRSTEPS to compatible values.'
!      call fexit()
!    end if
#ifdef ENABLE_MPI
  else if (re_aux(8).ne.nsim) then
    write(ilog,*) 'Warning. The processing of a PIGS trace in clustering graph-based analyses&
 & in parallel analysis mode requires setting FMCSC_RE_TRAJTOTAL and FMCSC_NRSTEPS to the same values. Adjusted.'
    re_aux(8) = nsim
  end if
#else
  else if (re_aux(8)*re_conditions.ne.nsim) then
    write(ilog,*) 'Fatal. The processing of a PIGS trace in clustering graph-based analyses&
 & requires setting FMCSC_RE_TRAJTOTAL, FMCSC_REPLICAS, and FMCSC_NRSTEPS to compatible values.'
    call fexit()
  end if
#endif
!
  if ((re_conditions.le.0).OR.(re_aux(4).le.0).OR.(re_aux(4).ge.cstored)) then
    write(ilog,*) 'Fatal. Inconsistent parameters in PIGS_ReadTrace(...). Inferring connectivity changes from &
 &PIGS trace file requires setting FMCSC_RE_TRAJOUT, FMCSC_RE_TRAJSKIP, and FMCSC_REPLICAS appropriately.'
    call fexit()
  end if
  if (mod(cstored,re_conditions).ne.0) then
    write(ilog,*) 'Fatal. The number of processed snapshots is not an exact multiple of the number of assumed replicas &
 &in PIGS_ReadTrace(...). This is likely caused by FMCSC_FRAMESFILE, FMCSC_NRSTEPS, or FMCSC_CCOLLECT.'
    call fexit()
  end if
  spert = cstored/re_conditions
  sstep = re_aux(8) !  - re_aux(5)/re_aux(4)
  xstep = re_aux(5)-mod(re_aux(5),re_aux(4)) ! the last snapshot not to have been written in original (trace) numbering
  k = sstep
  kk = 0
  do istep=1,cstored-1
    if (tlst(istep,3).gt.k) then
      if (kk.ne.spert) then
        write(ilog,*) 'Fatal. The number of processed snapshots per replica is not constant. This is likely caused by &
 &FMCSC_FRAMESFILE, FMCSC_NRSTEPS, or FMCSC_CCOLLECT.'
        call fexit()
      end if
      k = k + sstep
      kk = 1
    else
      kk = kk + 1
    end if
  end do
  kk = kk + 1
  if (kk.ne.spert) then
    write(ilog,*) 'Fatal. The number of processed snapshots per replica is not constant. This is likely caused by &
 &FMCSC_FRAMESFILE, FMCSC_NRSTEPS, or FMCSC_CCOLLECT.'
    call fexit()
  end if
!
  do istep=spert,cstored-1,spert
    tlst(istep,2) = 0 ! breaks between entire copies
    tlst(istep+1,5) = 0
  end do
!
! first sanity check and read
  kstep = 0
  kk = 0
  do while (kstep.le.(re_aux(8)*re_aux(4)+xstep))
    read(iu,*,iostat=iomess) kstep,lmap(1:re_conditions)
    if (kstep.gt.re_aux(8)*re_aux(4)+xstep) exit
    if (iomess.eq.IOSTAT_END) then
      allout = .true.
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while reading PIGS trace file (...).'
      call fexit()
    end if
    kk = kk + 1
    do k=1,re_conditions
      if ((lmap(k).le.0).OR.(lmap(k).gt.re_conditions)) then
        write(ilog,*) 'Warning. Out-of-range replica number while reading PIGS trace file (...). &
 &FMCSC_RE_TRAJSKIP, FMCSC_RE_TRAJOUT, and FMCSC_REPLICAS must match the trace file. Assuming no reseeding occurred.'
        lmap(k) = k
      end if
    end do
  end do
  rewind(unit=iu)
  klines = kk
  allocate(lmapmat(re_conditions+1,klines))
  do i=1,klines
    read(iu,*,iostat=iomess) lmapmat(:,i)
    do k=1,re_conditions
      if ((lmapmat(k+1,i).le.0).OR.(lmapmat(k+1,i).gt.re_conditions)) then
        lmapmat(k+1,i) = k
      end if
    end do
  end do
  kk = 1
  oldwhich = 0
  do istep=1,cstored-1
    if (tlst(istep,2).eq.0) cycle
    which = tlst(istep,3)/sstep + 1
    if (which.ne.oldwhich) then
      kk = 1
      oldwhich = which
    end if
    kstep = lmapmat(1,kk) + (which-1)*sstep*re_aux(4) - xstep
    do while ((tlst(istep,3)*re_aux(4).gt.kstep).AND.(kk.lt.klines))
      kk = kk + 1
      kstep = lmapmat(1,kk) + (which-1)*sstep*re_aux(4) - xstep
      if (kk.eq.klines) exit
    end do
    if (mod(istep,spert).eq.0) cycle
    if ((tlst(istep,3)*re_aux(4).le.kstep).AND.(tlst(istep+1,3)*re_aux(4).gt.kstep)) then
!     track back reseedings in this interval
      history = 0
      kkk = kstep
      if (kk.lt.klines) then
        do while (tlst(istep+1,3)*re_aux(4).gt.kkk)
          if (lmapmat(which+1,kk).ne.which) then
            history = history + 1
          end if
          kk = kk + 1
          if (kk.gt.klines) exit
          kkk = lmapmat(1,kk) + (which-1)*sstep*re_aux(4) - xstep
        end do
        kk = kk - 1
      end if
!      write(*,*) kstep,history,kk,kkk,xstep
!     if there were any, we have work to do
      if (history.gt.0) then
        lstep = lmapmat(which+1,kk)
!        write(*,*) 'hit ',istep,tlst(istep,3),tlst(istep+1,3),which,lstep
        tlst(istep,1) = 0
        kkk = kk
        recwhich = which
        ii = istep
        minstep = tlst(istep,3) - (which-1)*sstep
        do while ((tlst(istep,1).eq.0).AND.(lstep.gt.0))
          i = lstep*spert
          do while (((i.gt.(lstep-1)*spert).OR.(kkk.gt.0)).AND.(tlst(istep,1).eq.0))
!           this block traces back the trajectory until a relevant step number is reached (redone for every track switch)
            if (i.gt.(lstep-1)*spert) then
              if ((tlst(i,3)-(lstep-1)*sstep).gt.minstep) then
!                write(*,*) 'skipping ',tlst(i,3)-(lstep-1)*sstep,lstep,recwhich,minstep
                i = i - 1
                cycle
              end if
            end if
!           if our current relevant snapshot is a break point, exit
            if ((ntbrks.gt.0).AND.(i.gt.(lstep-1)*spert)) then
              j = i ! tlst(i,3)
              call binary_search(ntbrks,trbrkslst(1:ntbrks),j,jj)
              if (trbrkslst(max(1,min(ntbrks,jj))).eq.i) then
!                write(*,*) 'Killing ',istep,' recovery due to break crossing (',trbrkslst(max(1,min(ntbrks,jj))),').'
                tlst(istep,1) = istep
                tlst(tlst(istep,2),5) = 0
                tlst(istep,2) = 0
                exit
              end if
            end if
!           if we ran out of trajectory in the current track, we may still be able to switch and recover
            if (i.eq.(lstep-1)*spert) then
!              write(*,*) 'have to decrease from ',kkk,lmapmat(1,kkk)
              if (kkk.gt.0) then
                if (lmapmat(lstep+1,kkk).ne.lstep) then
                  recwhich = lstep
                  lstep = lmapmat(lstep+1,kkk)
                  minstep = min(minstep,(lmapmat(1,kkk)-xstep)/re_aux(4))
!                  write(*,*) 'have to switch track from ',recwhich,' to ',lstep,istep
                  exit
                else
                  kkk = kkk - 1
!                  minstep = min(minstep,(lmapmat(1,kkk)-xstep)/re_aux(4)+1) ! integer division
!                  write(*,*) 'could decrease kkk to ',kkk
                end if
              end if
!           no more reseedings to handle means we accept
            else if (kkk.eq.0) then
              tlst(istep,1) = istep
              tlst(tlst(istep,2),5) = i
!              write(*,*) 'Would take2 ',istep,i!,tlst(i,3), ' to ',tlst(tlst(istep,2),3)
              i = i - 1
              if (i.eq.(lstep-1)*spert) exit
              minstep = min(minstep,tlst(i,3)-(lstep-1)*sstep)
!           the last two blocks deal with tracking back either trace (first) or trajectory (second) and possibly accepting
            else if (re_aux(4)*tlst(i,3).le.(lmapmat(1,kkk)+re_aux(4)*(lstep-1)*sstep-xstep)) then
!              write(*,*) 'have to decrease from ',kkk,lmapmat(1,kkk)
              if (kkk.gt.0) then
                if (lmapmat(lstep+1,kkk).ne.lstep) then
                  recwhich = lstep
                  lstep = lmapmat(lstep+1,kkk)
                  ii = i
                  minstep = min(minstep,(lmapmat(1,kkk)-xstep)/re_aux(4))
!                  write(*,*) 'have to switch track from ',recwhich,' to ',lstep,minstep
                  exit
                else
                  kkk = kkk - 1
!                 the minimum step should be the one that is within 
!                  minstep = min(minstep,(lmapmat(1,kkk)-xstep)/re_aux(4)+1) ! integer division
!                  write(*,*) 'could decrease kkk to ',kkk
                end if
              end if
            else if (re_aux(4)*tlst(i,3).ge.(lmapmat(1,kkk)+re_aux(4)*(lstep-1)*sstep-xstep)) then
              tlst(istep,1) = istep
              tlst(tlst(istep,2),5) = i
!              write(*,*) 'Would take ',istep,i!,tlst(i,3), ' to ',tlst(tlst(istep,2),3)
              i = i - 1
              if (i.eq.(lstep-1)*spert) exit
              minstep = min(minstep,tlst(i,3)-(lstep-1)*sstep)
            end if
!            write(*,*) 'searching ',lstep,tlst(i,3),kkk
          end do
          if (kkk.gt.0) then
            if (lmapmat(1,kkk).le.re_aux(5)) then
!              write(*,*) 'Did run out ',istep
              tlst(istep,1) = istep
              tlst(tlst(istep,2),5) = 0
              tlst(istep,2) = 0
              exit
            end if
          end if
          if (i.eq.(lstep-1)*spert) then
            lstep = 0
!            write(*,*) 'would have run out',i
          end if
        end do
      end if
      if (kk.lt.klines) kk = kk + 1
    end if
!    if (minval(tlst(istep,1:2)).gt.0) then
!      call snap_to_snap_d(dis,tlst(istep,1),tlst(istep,2))
!      write(*,*) 'HOP ',istep,dis,mod(tlst(tlst(istep,2),3),sstep)-mod(tlst(tlst(istep,1),3),sstep)
!    end if
  end do
!
  close(unit=iu)
  deallocate(lmapmat)
!
end
!
!-----------------------------------------------------------------------------------------------------------------------
!
! a subroutine to build a default connectivity map based on a lag time (SW by default)
!
! structure of the map after completion:
! sconnect(:,1) is index of "from" relative to cstored (cludata)
! sconnect(:,2) is index of "to" relative to cstored (cludata)
! sconnect(:,3) is map of the index itself (not necessarily "from" or "to") to original index reference
!               (lost in random-access frames file)
! sconnect(:,4:5) are used as temporaries here
! sconnect(:,4) will later hold the CG trajectory 
!
subroutine build_sconnect(lagt,mode)
!
  use clusters
  use mpistuff
  use iounit
  use interfaces
  use pdb, ONLY: framecnt,framelst,select_frames,pdb_fileformat,use_frame_weights
  use system, ONLY: nsim,nequil
!
  implicit none
!
  integer, INTENT(IN):: lagt,mode
!
  integer i,j,k,ilag,t1,t2,spert,nrep,bhlp,jj1,jj2,lkcnt,cst_bu,csttmp
  logical have_trace,have_warned
!
! if still relevant, read breaks and links files
  call read_trajbrksfile()
  call read_trajlinksfile()
!
! when called after data are read and cstored is populated
  if (mode.eq.1) then
    csttmp = cstored
! when called after parameters are set but before data are read (cstored is 0)
  else if (mode.eq.2) then
    csttmp = 0
    if ((select_frames.EQV..true.).AND.(((pdb_fileformat.ne.5).AND.(pdb_fileformat.ne.2)).OR.(use_frame_weights.EQV..true.))) then
      k = 1
      do i=cstorecalc,nsim,cstorecalc
        do while (k.le.framecnt)
          if (framelst(k).eq.i) csttmp = csttmp + 1
          if (framelst(k).ge.i) exit
          k = k + 1
        end do
      end do
    else
      k = nsim
      if (select_frames.EQV..true.) k = framecnt
      do i=cstorecalc,k,cstorecalc
        if (i.gt.nequil) then
          csttmp = csttmp + 1
        end if
      end do
#ifdef ENABLE_MPI
      if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiseed.EQV..false.)) then
        csttmp = csttmp*re_conditions
      end if
#endif
    end if
  else
    write(ilog,*) 'Called build_sconnect(...) with unsupported mode (got ',mode,'). This is a bug.'
    call fexit()
  end if
!
! max allocation size is cstored + ntlnks
  sconnectsz = csttmp+ntlnks
  allocate(sconnect(sconnectsz,5))
  sconnect(:,:) = 0
!
  call strlims(re_traceinfile,t1,t2)
  inquire(file=re_traceinfile(t1:t2),exist=have_trace)
#ifdef ENABLE_MPI
  if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiseed.EQV..false.)) then
!   do nothing    
  else
    have_trace = .false. ! PIGS analysis mode or REX analysis modes never read trace 
  end if
#endif
!
! create a map
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    i = cstorecalc*floor((re_conditions*nsim+0.1)/cstorecalc)
  else
    i = cstorecalc*floor((nsim+0.1)/cstorecalc)
  end if
#else
  i = cstorecalc*floor((nsim+0.1)/cstorecalc)
#endif
! the default map with lag time 1
  j = csttmp
  if ((select_frames.EQV..true.).AND.(((pdb_fileformat.ne.5).AND.(pdb_fileformat.ne.2)).OR.(use_frame_weights.EQV..true.))) then
    k = framecnt
    do while ((i.ge.1).AND.(j.ge.1).AND.(k.ge.1))
      do while (i.gt.framelst(k))
        i = i - cstorecalc
      end do
      if (i.eq.framelst(k)) then
        sconnect(j,3) = i
        i = i - cstorecalc
        j = j - 1
      end if
      k = k - 1
    end do      
  else
    if (select_frames.EQV..true.) i = framecnt
    do while ((i.ge.1).AND.(j.ge.1))
      sconnect(j,3) = i
      i = i - cstorecalc
      j = j - 1
    end do
  end if
!
  do i=1,csttmp-1
    sconnect(i,1) = i
    sconnect(i,2) = i+1
    sconnect(i+1,5) = i
  end do
  sconnect(csttmp,1) = csttmp
!
  if (have_trace.EQV..true.) then
    nrep = re_conditions
    spert = csttmp/re_conditions
    cst_bu = cstored
    cstored = csttmp
    call PIGS_ReadTrace(sconnect,sconnectsz) ! remaps sconnect(:,1)
  else
    nrep = 1
    spert = csttmp
    cstored = csttmp
  end if
!
  if (lagt.ge.spert) then
    write(ilog,*) 'Fatal. Called build_sconnect with an ineligible lag time. This is a bug.'
    call fexit()
  end if
!
! add custom breaks (still assuming lag time 1)
  if (ntbrks.gt.0) then
    bhlp = csttmp
    do i=1,ntbrks
      j = trbrkslst(i) 
      call binary_search(bhlp,sconnect(1:csttmp,3),j,k)
      if (sconnect(max(1,min(bhlp,k)),3).eq.trbrkslst(i)) then
        if (sconnect(max(1,min(bhlp,k)),2).gt.0) then
          if (sconnect(sconnect(max(1,min(bhlp,k)),2),5).eq.max(1,min(bhlp,k))) sconnect(sconnect(max(1,min(bhlp,k)),2),5) = 0
          sconnect(max(1,min(bhlp,k)),2) = 0
        end if
      else if ((k.le.1).OR.(k.gt.bhlp)) then
!       do nothing (out of range)
      else
        if (sconnect(sconnect(k,2),5).eq.k) sconnect(sconnect(k,2),5) = 0
        sconnect(k,2) = 0
      end if
    end do
  end if
!
! the same for links
  lkcnt = 0
  if (ntlnks.gt.0) then
    bhlp = csttmp
    have_warned = .false.
    do i=1,ntlnks
      j = trlnkslst(i,1)
      call binary_search(bhlp,sconnect(1:csttmp,3),j,k)
      if (sconnect(max(1,min(bhlp,k)),3).eq.trlnkslst(i,1)) then
        jj1 = max(1,min(bhlp,k))
      else if (((k.le.0).OR.(k.gt.bhlp)).OR.(remap_trajlnks.EQV..false.)) then
        jj1 = 0
      else
        jj1 = k
      end if
      j = trlnkslst(i,2)
      call binary_search(bhlp,sconnect(1:csttmp,3),j,k)
      if (sconnect(max(1,min(bhlp,k)),3).eq.trlnkslst(i,2)) then
        jj2 = max(1,min(bhlp,k))
      else if (((k.le.0).OR.(k.ge.bhlp)).OR.(remap_trajlnks.EQV..false.)) then
        jj2 = 0
      else
        jj2 = k+1
      end if
      if (jj1.gt.0) then
        lkcnt = lkcnt + 1
        sconnect(csttmp+lkcnt,1) = jj1
        sconnect(csttmp+lkcnt,3) = sconnect(jj1,3)
        sconnect(csttmp+lkcnt,5) = sconnect(jj1,5)
        if (jj2.gt.0) then
          sconnect(csttmp+lkcnt,2) = jj2
        end if
      end if
      if (((jj1.le.0).OR.(jj2.le.0)).AND.(have_warned.EQV..false.)) then
        write(ilog,*) 'Warning. At least one custom link requested through FMCSC_TRAJLINKSFILE could not be used because frames &
 &are not present (by action of FMCSC_CCOLLECT or FMCSC_FRAMESFILE).'
        have_warned = .true.
      end if
!      if ((jj1.le.0).OR.(jj2.le.0)) then
!        write(*,*) trlnkslst(i,1:2),jj1,jj2
!      end if
    end do
!   if remapping of links is on, check that through trace or breaks none of the remapped connections could be compromised
!   if so, disable link
    if (remap_trajlnks.EQV..true.) then
      have_warned = .true.
      do i=csttmp+1,csttmp+ntlnks
        if ((sconnect(i,1).gt.0).AND.(sconnect(i,2).gt.0)) then
          if (trlnkslst(i-csttmp,1).ne.sconnect(i,3)) then
            if ((sconnect(sconnect(i,1),1).ne.sconnect(i,1)).OR.(sconnect(sconnect(i,1),2).le.0)) then
              sconnect(i,2) = 0
              have_warned = .false.
              cycle
            end if
          end if
          if (trlnkslst(i-csttmp,2).ne.sconnect(sconnect(i,2),3)) then
            if (sconnect(i,2).gt.1) then
              if ((sconnect(sconnect(i,2)-1,1).ne.sconnect(i,2)-1).OR.(sconnect(sconnect(i,2)-1,2).le.0)) then
                sconnect(i,2) = 0
                have_warned = .false.
                cycle
              end if
            end if
          end if
        end if
      end do
      if (have_warned.EQV..false.) then
        write(ilog,*) 'Warning. At least one custom link requested through FMCSC_TRAJLINKSFILE could not be used because breaks &
 &interfere (FMCSC_TRAJBREAKSFILE or FMCSC_TRACEFILE).'
      end if
    end if
  end if
!
! now translate everything (including custom links) to final lag time by simply following the reverse map
  sconnect(:,4) = sconnect(:,5)
  do ilag=1,lagt-1
    do i=1,csttmp
      if (sconnect(i,4).gt.0) then
        sconnect(i,4) = sconnect(sconnect(i,4),5)
      else
        sconnect(i,4) = 0
      end if
    end do
  end do
!
! and rewrite
  k = 0
  do i=1,csttmp
    if (sconnect(i,4).gt.0) then
      k = k + 1
      sconnect(k,1) = sconnect(i,4)
      sconnect(k,2) = sconnect(i,1)
    end if
  end do
  do i=k+1,csttmp
    sconnect(i,1:2) = 0
  end do
!
  if (brklnk_report.ge.1) then
  456 format(100(i10,1x))
    write(ilog,*)
    write(ilog,'(a)') '--------------- REPORT OF SNAPSHOT LINK STRUCTURE ------------------'
    if (brklnk_report.eq.1) then
      write(ilog,'(a,i7,a)') ' Only links with a lag time different from ',clagt_msm,' stored snapshots are printed.'
    end if
    write(ilog,'(a)') '  "FROM"  |   "TO"   | INDEX "F"| INDEX "T"|'
    do i=1,sconnectsz
      if ((sconnect(i,1).le.0).OR.(sconnect(i,2).le.0)) cycle
      if ((brklnk_report.ne.1).OR.((sconnect(i,2)-sconnect(i,1)).ne.clagt_msm)) then
        k = k -1
        write(ilog,456) (sconnect(i,j),j=1,2),(sconnect(sconnect(i,j),3),j=1,2)
      end if
    end do
    if (brklnk_report.eq.1) then
      write(ilog,'(a,i9,a)') ' In addition, there are ',k,' links with the expected time lag.'
    end if
    write(ilog,'(a)') '----------- END OF REPORT OF SNAPSHOT LINK STRUCTURE ---------------'
    write(ilog,*)
  end if
  sconnect(:,4:5) = 0
!
! lastly, remap additional input snapshots where needed
  if ((eigvalmd.ne.0).OR.(dopfold.EQV..true.).OR.(synmode.ne.0)) then
!   Initial snapshot 
    if (inissnap_usrslct.gt.0) then !what if we have a frames file or ccollect is not 1 or both: remapping through sconnect(:,3)
      bhlp = csttmp
      j = inissnap_usrslct
      call binary_search(bhlp,sconnect(1:bhlp,3),j,k)
      if (sconnect(max(1,min(bhlp,k)),3).eq.inissnap_usrslct) then
        write(ilog,*) 'Identified initial snapshot selected by FMCSC_INISYNSNAP (',inissnap_usrslct,') as ',j,'th stored snapshot.'
        inissnap_usrslct = max(1,min(bhlp,k))
      else
        write(ilog,*) 'Fatal. Initial snapshot selected by FMCSC_INISYNSNAP could not be found amongst stored snapshots.'
        call fexit() ! to help the user
      end if
    end if
  end if
!
! end snapshot 
  if ((synmode.eq.1).OR.(dopfold.EQV..true.)) then  !cases where end snapshot is important
!   Initial snapshot 
    if (endssnap_usrslct.gt.0) then !what if we have a frames file or ccollect is not 1 or both: remapping through sconnect(:,3)
      bhlp = csttmp
      j = endssnap_usrslct
      call binary_search(bhlp,sconnect(1:bhlp,3),j,k)
      if (sconnect(max(1,min(bhlp,k)),3).eq.endssnap_usrslct) then
        write(ilog,*) 'Identified target snapshot selected by FMCSC_ENDSYNSNAP (',endssnap_usrslct,') as ',j,'th stored snapshot.'
        endssnap_usrslct = max(1,min(bhlp,k))
      else
        write(ilog,*) 'Fatal. Target snapshot selected by FMCSC_ENDSYNSNAP could not be found amongst stored snapshots.'
        call fexit() ! to help the user
      end if
    end if
  end if
!
  if (dopfold.EQV..true.) then
    k = 1
    call read_snapsetfile(sconnect(:,1:3),sconnectsz,k) ! reads F-set to mapped snapshots
    if (allocated(clufold).EQV..false.) then
      if (inissnap_usrslct.gt.0) then
        allocate(clufold(1))
        clufold(1) = inissnap_usrslct
      end if
    end if
    k = 2
    call read_snapsetfile(sconnect(:,1:3),sconnectsz,k) ! reads U-set to mapped snapshots
    if (allocated(cluunfold).EQV..false.) then
      if (endssnap_usrslct.gt.0) then
        allocate(cluunfold(1))
        cluunfold(1) = endssnap_usrslct
      end if
    end if
  end if
!
!
end
!
!-------------------------------------------------------------------------------------------------------------------------
!
subroutine build_graph_from_file()
!
  use clusters
!
  implicit none
!
  integer iunit,freeunit,i,nba,j,k,ik,jk,ii,ll,iu,whichscc
  logical foundit,dofwr
  integer, ALLOCATABLE:: remap(:)
!
  type(t_scluster), ALLOCATABLE:: it(:)  ! sparse transition matrix
  RTYPE ptmp
  RTYPE, ALLOCATABLE:: p1(:)
 
!
  iunit = freeunit()
  open(unit=iunit,file="tmat_abeta_add0_lag20/STRUCT_CLUSTERING.clu",status='old')
!
  do i=1,cstored
    read(iunit,*) sconnect(i,4)
  end do
!
  close(unit=iunit)
!
  nba = maxval(sconnect(1:cstored,4))
  write(*,*) nba
!
  allocate(it(nba))
!
  do i=1,nba
    if (allocated(it(i)%wghtsnb).EQV..true.) deallocate(it(i)%wghtsnb)
    if (allocated(it(i)%lensnb).EQV..true.) deallocate(it(i)%lensnb)
    if (allocated(it(i)%fewtsnb).EQV..true.) deallocate(it(i)%fewtsnb)
    if (allocated(it(i)%map).EQV..true.) deallocate(it(i)%map)
    if (allocated(it(i)%lstnb).EQV..true.) deallocate(it(i)%lstnb)
    if (allocated(it(i)%flwnb).EQV..true.) deallocate(it(i)%flwnb)
    it(i)%nbalsz = 2
    allocate(it(i)%wghtsnb(it(i)%nbalsz,2))
    allocate(it(i)%lensnb(it(i)%nbalsz,2))
    allocate(it(i)%fewtsnb(it(i)%nbalsz,4))
    it(i)%wghtsnb(:,:) = 0
    it(i)%lensnb(:,:) = 0.0
    it(i)%fewtsnb(:,:) = 0.0
    allocate(it(i)%lstnb(it(i)%nbalsz))
    allocate(it(i)%map(it(i)%nbalsz))
    it(i)%nb = 0
    it(i)%nmbrs = 0
    it(i)%nodewt(:) = 0.0
  end do
!
  do i=1,cstored
    j = sconnect(i,4)
    it(j)%nodewt(1) = it(j)%nodewt(1) + 1.0
  end do
  ptmp = 1.0/(1.0*cstored)
  it(:)%nodewt(1) = ptmp*it(:)%nodewt(1)
!
  do i=1,sconnectsz
    if ((sconnect(i,1).le.0).OR.(sconnect(i,2).le.0)) cycle
    ll = sconnect(sconnect(i,1),4)
    ii = sconnect(sconnect(i,2),4)
    
    foundit = .false.
    do j=1,it(ii)%nb
      if (it(ii)%lstnb(j).eq.ll) then
        it(ii)%wghtsnb(j,2) = it(ii)%wghtsnb(j,2) + 1
!        it(ii)%lensnb(j,2) = it(ii)%lensnb(j,2) + dval
        it(ii)%fewtsnb(j,2) = it(ii)%fewtsnb(j,2) + 1.0
        foundit = .true.
        ik = j
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      it(ii)%nb = it(ii)%nb + 1
      ik = it(ii)%nb
      if (it(ii)%nb.gt.it(ii)%nbalsz) then
        call scluster_resizenb(it(ii))
      end if
!      it(ii)%map(ll) = it(ii)%nb
      it(ii)%lstnb(it(ii)%nb) = ll
      it(ii)%wghtsnb(it(ii)%nb,2) = 1
      it(ii)%wghtsnb(it(ii)%nb,1) = 0
!      it(ii)%lensnb(it(ii)%nb,2) = dval
!      it(ii)%lensnb(it(ii)%nb,1) = 0.0
      it(ii)%fewtsnb(it(ii)%nb,2) = 1.0
      it(ii)%fewtsnb(it(ii)%nb,1) = 0.0
    end if
    foundit = .false.
    do j=1,it(ll)%nb
      if (it(ll)%lstnb(j).eq.ii) then
        it(ll)%wghtsnb(j,1) = it(ll)%wghtsnb(j,1) + 1
!        it(ll)%lensnb(j,1) = it(ll)%lensnb(j,1) + dval
        it(ll)%fewtsnb(j,1) = it(ll)%fewtsnb(j,1) + 1.0
        foundit = .true.
        jk = j
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      it(ll)%nb = it(ll)%nb + 1
      jk = it(ll)%nb
      if (it(ll)%nb.gt.it(ll)%nbalsz) then
        call scluster_resizenb(it(ll))
      end if
!      it(ll)%map(ii) = it(ll)%nb
      it(ll)%lstnb(it(ll)%nb) = ii
      it(ll)%wghtsnb(it(ll)%nb,1) = 1
      it(ll)%wghtsnb(it(ll)%nb,2) = 0
!      it(ll)%lensnb(it(ll)%nb,1) = dval
!      it(ll)%lensnb(it(ll)%nb,2) = 0.0
      it(ll)%fewtsnb(it(ll)%nb,1) = 1.0
      it(ll)%fewtsnb(it(ll)%nb,2) = 0.0
    end if
!   we store the index for the list of node 2 (ll ; jk) in node 1 at the position of the neighbor (ii ; ik)
    it(ii)%map(ik) = jk
!   we store the index for the list of node 1 (ii ; ik) in node 2 at the position of the neighbor (ll ; jk)
    it(ll)%map(jk) = ik
  end do
!
  call augment_network(it,nba,caddlkmd)
  k=MAXVAL(it(:)%inscc)
  allocate(p1(k))
  p1(:) = 0.
  do i=1,nba
    j=it(i)%inscc
    p1(j) = p1(j) + it(i)%nodewt(1)
  end do
  whichscc = MAXLOC(p1,dim=1)
  deallocate(p1)
!
  dofwr = .true.
  if (dofwr.EQV..true.) then
    jk = 1
  else
    jk = 2
  end if
!
  allocate(remap(nba))
  remap(:) = 0
  ik = 0
  do i=1,nba
    if (it(i)%inscc.eq.whichscc) then
      ik = ik + 1
      remap(i) = ik
    end if
    ptmp = 0.
    do j=1,it(i)%nb
      if (it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc) cycle
      ptmp = ptmp + it(i)%fewtsnb(j,jk)
    end do
    if (ptmp.gt.0.0) then
      ptmp = 1.0/ptmp
      it(i)%fewtsnb(1:it(i)%nb,3) = ptmp*it(i)%fewtsnb(1:it(i)%nb,jk)
      do j=1,it(i)%nb
        if (it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc) it(i)%fewtsnb(j,3) = 0.0
      end do
    end if
  end do
!
  if (tmat_report.EQV..true.) then
    iu = freeunit()
    if (dofwr.EQV..true.) then
      open(unit=iu,file='map_tmat_fwr.dat',status='unknown')
    else
      open(unit=iu,file='map_tmat_bwr.dat',status='unknown')
    end if
    do i=1,nba
      if (it(i)%inscc.ne.whichscc) cycle
      do j=1,it(i)%nb
        if (it(i)%fewtsnb(j,3).gt.0.0) then
          write(iu,'(I6,2x,I6,2x,I6,2x,I6,2x,G19.11,2x,G14.6)') remap(i),remap(it(i)%lstnb(j)),i,it(i)%lstnb(j),&
 &it(i)%fewtsnb(j,3),it(i)%fewtsnb(j,1)
        end if
      end do
    end do
    close(unit=iu)
  end if
!
end
!
!-------------------------------------------------------------------------------------------------------------
!
! mode 0: row-normalized transition matrix in columns 3 and 4 of fewtsnb based on wghtsnb(:,1:2)
! mode 1: row-normalized transition matrix in columns 3 and 4 of fewtsnb based on fewtsnb(:,1:2)
! mode 2: same as 0 but putting the cumulative sum
! mode 3: same as 1 but putting the cumulative sum
!
! mode2 1: cycle for edges bridging separate SCCs
! mode2 else: do not cycle
!
subroutine edit_tmat(it,nbasins,mode,mode2)
!
  use clusters
!
  implicit none
!
  integer, INTENT(IN):: nbasins,mode,mode2
!
  type(t_scluster) it(nbasins)
!
  integer i,j,jk
  RTYPE ptmp
  RTYPE, ALLOCATABLE:: pvec(:)
  logical atrue
!
  atrue = .true.
!
  if ((mode.eq.2).OR.(mode.eq.3)) allocate(pvec(maxval(it(1:nbasins)%nb)))
!
  do i=1,nbasins
    do jk=1,2
      ptmp = 0.
      if ((mode.eq.0).OR.(mode.eq.2)) then
        do j=1,it(i)%nb
          if ((it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc).AND.(mode2.eq.1)) cycle
          ptmp = ptmp + 1.0*it(i)%wghtsnb(j,jk)
        end do
        if (ptmp.gt.0.0) then
          ptmp = 1.0/ptmp
          it(i)%fewtsnb(1:it(i)%nb,jk+2) = ptmp*it(i)%wghtsnb(1:it(i)%nb,jk)
          if (mode2.eq.1) then
            do j=1,it(i)%nb
              if (it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc) it(i)%fewtsnb(j,jk+2) = 0.0
            end do
          end if
        end if
      else if ((mode.eq.1).OR.(mode.eq.3)) then
        do j=1,it(i)%nb
          if ((it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc).AND.(mode2.eq.1)) cycle
          ptmp = ptmp + it(i)%fewtsnb(j,jk)
        end do
        if (ptmp.gt.0.0) then
          ptmp = 1.0/ptmp
          it(i)%fewtsnb(1:it(i)%nb,jk+2) = ptmp*it(i)%fewtsnb(1:it(i)%nb,jk)
          if (mode2.eq.1) then
            do j=1,it(i)%nb
              if (it(i)%inscc.ne.it(it(i)%lstnb(j))%inscc) it(i)%fewtsnb(j,jk+2) = 0.0
            end do
          end if
        end if
      end if
      if ((mode.eq.2).OR.(mode.eq.3)) then
        call cumsum(it(i)%fewtsnb(1:it(i)%nb,jk+2),it(i)%nb,atrue,pvec(1:it(i)%nb))
        it(i)%fewtsnb(1:it(i)%nb,jk+2) = pvec(1:it(i)%nb)
      end if
    end do
  end do
!
  if (allocated(pvec).EQV..true.) deallocate(pvec)
!
end
!
!--------------------------------------------------------------------------------
!
! a routine to print all nonzero elements in it(:)fewtsnb(:,fbix), either for all (sccix is zero)
! nbasins, or remapped to an SCC (sccix is nonzero)
!
! the values in fewtsnb and of suffix must be set meaningfully on the outside
!
subroutine prt_tmat(it,nbasins,fbix,suffix,sccix)
!
  use iounit
  use clusters
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer, INTENT(IN):: fbix,nbasins,sccix
  type (t_scluster), INTENT(IN):: it(nbasins)
  character(3), INTENT(IN):: suffix
!
  integer i,j,iu,freeunit,t1,t2,tl,k
  character(MAXSTRLEN) fn
  character(6) compo
#ifdef ENABLE_MPI
  character(3) nod
#endif
!
  iu = freeunit()
  k = 2 - mod(fbix,2)
!
  if (sccix.gt.0) then
    tl = 6
    call int2str(sccix,compo,tl)
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'TMAT_'//compo//'_'//suffix//'.dat'
    else if (use_MPIAVG.EQV..true.) then
      if (myrank.ne.0) return
      fn = 'TMAT_'//compo//'_'//suffix//'.dat'
    end if
#else
    fn = 'TMAT_'//compo//'_'//suffix//'.dat'
#endif
  else
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      tl = 3
      call int2str(myrank,nod,tl)
      fn = 'N_'//nod(1:tl)//'TMAT_'//suffix//'.dat'
    else if (use_MPIAVG.EQV..true.) then
      if (myrank.ne.0) return
      fn = 'TMAT_'//suffix//'.dat'
    end if
#else
    fn = 'TMAT_'//suffix//'.dat'
#endif
  end if
!
  call strlims(fn,t1,t2)
  open(unit=iu,file=fn(t1:t2),status='unknown')
!
  do i=1,nbasins
    if ((it(i)%inscc.ne.sccix).AND.(sccix.gt.0)) cycle
    do j=1,it(i)%nb
      if (it(i)%fewtsnb(j,fbix).gt.0.0) then
        if (sccix.gt.0) then
          write(iu,'(I7,1x,I7,1x,I7,1x,I7,1x,G19.11,2x,I9)') &
 & it(i)%ix(2),it(it(i)%lstnb(j))%ix(2),i,it(i)%lstnb(j),it(i)%fewtsnb(j,fbix),it(i)%wghtsnb(j,k)
        else
          write(iu,'(I7,1x,I7,1x,G19.11,2x,I9)') i,it(i)%lstnb(j),it(i)%fewtsnb(j,fbix),it(i)%wghtsnb(j,k)
        end if
      end if
    end do
  end do
  close(unit=iu)
!
end
!
!--------------------------------------------------------------------------------
!
! check a vector for whether it makes a possible steady-state for the SCC defined by sccix
!
! specifically: 1) check for sign change; 2) check for zeros; 3) check for max ratio
! it returns safe_ss as -1 if criteria 1 or 2 are violated, 0 if 3 is violated, and 1 otherwise
!
subroutine check_msm_ss(nsynclu,eigvec_msm,safe_ss,it,nbasins,ssixin,sccix)
!
  use iounit
  use clusters, ONLY: t_scluster,csccwts
!
  implicit none
!
  integer, INTENT(IN):: nsynclu,nbasins,sccix,ssixin
  integer, INTENT(OUT):: safe_ss
  type(t_scluster), INTENT(IN):: it(nbasins)
  RTYPE eigvec_msm(1:nsynclu)
!
  RTYPE, ALLOCATABLE:: tmpvr(:)
  integer i,k,counter,ssix
  RTYPE maxthresh,minthresh,eigsum
!
  safe_ss = 1
!
  if ((ssixin.lt.2).OR.(ssixin.gt.5)) then
    write(ilog,*) 'Warning. Called check_msm_ss(...) with a nonsensical setting for the input type. This is a(n omission) bug.'
  end if
!
  maxthresh = maxval(eigvec_msm(:))
  if (maxthresh.le.0) then                !reverse sign
    do i=1,nsynclu
      eigvec_msm(i) = -1.*eigvec_msm(i)
    end do  
  end if
  maxthresh = maxval(eigvec_msm(:))
  !write(ilog,*) 'Maxthresh',maxthresh
  minthresh = 1.*maxthresh/1.E06          !6 orders of magnitude smaller, heuristic for possible problems detection
  !write(ilog,*) 'Minthresh',minthresh
!
  do i=1,nsynclu
    if (eigvec_msm(i).le.0.0) then ! presence of exact zeroes
      write(ilog,*) 'Warning. Found an entry that is strictly zero in an assumed steady state population vector. This &
 &indicates a numerically ill-defined graph, failure of the routine computing the vector to converge, or a bug.'
      safe_ss = -1
      return
    end if
    if (eigvec_msm(i)*eigvec_msm(1).lt.0.) then !change in sign
      write(ilog,*) 'Warning. Found a change in sign in the entries in an assumed steady state population vector. This &
 &indicates a numerically ill-defined input graph, failure of the routine computing the vector to converge, or a bug.'
      safe_ss = -1
      return
    end if
  end do
!
  counter = 0
  do i=1,nsynclu
    !write(ilog,*) 'eigvec_msm(i),minthresh',eigvec_msm(i),minthresh
    if (eigvec_msm(i).lt.minthresh) then
      safe_ss = 0
      counter = counter + 1
    end if
  end do
!
  eigsum = sum(eigvec_msm(:))
  eigvec_msm(:) = eigvec_msm(:)/eigsum    !normalization
!
  if (safe_ss.eq.0) then
    write(ilog,69) counter
  else
    write(ilog,*) 'An assumed steady state population vector spans maximum and minimum &
 &probabilities from ',minval(eigvec_msm),' to ',maxval(eigvec_msm),'.'
  end if
!
 69 format('WARNING. An assumed steady state population vector has ',i6,' nodes that are more than 6 orders of magnitude lighter &
 &than the heaviest node.',/,'This might create (or have created) severe numerical stability and convergence issues for any &
 &algorithm that relies on or implies this equilibrium distribution.')
!
  if (ssixin.gt.3) then
    ssix = ssixin - 2 ! cross-check HSL against iterative
  else if (ssixin.gt.1) then
    ssix = ssixin + 2 ! cross-check iterative against HSL
  else
    return
  end if
  counter = 0
  if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
  allocate(tmpvr(nsynclu))
  k = 0
! the values in nodewt always preserve the component raw weight -> adjust for 1:1 comparison 
! based on check/return structure above, we know that eigvec_msm is positive definite at least
  do i=1,nbasins
    if (it(i)%inscc.eq.sccix) then
      k = k + 1
      if (it(i)%nodewt(ssix).gt.0.0) then
        tmpvr(k) = 100.*abs(it(i)%nodewt(ssix)/csccwts(sccix)-eigvec_msm(k))/min(eigvec_msm(k),(it(i)%nodewt(ssix)/csccwts(sccix)))
        counter = counter + 1
      end if
    end if
  end do
  if (counter.gt.0) then
    write(ilog,*) 'The maximum difference between the assumed steady state vector and an available algorithm steady state &
 &probabilities is [%]: ',maxval(tmpvr(:))
  end if
  if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
!
end subroutine check_msm_ss
!
!----------------------------------------------------------------------------------------------------------
!
! It computes a user-defined set of eigenvalues and potentially eigenvectors by 
! using an external library that handles sparse-matrix algebra.
!
! This is the outer function that takes CAMPARI-data structures as input and converts them
! to the required vectors
! A sub-matrix is used based on sccix, the user choice for SCC
!
! This outer routine could also provide an interface to the corresponding LAPACK solver
!
subroutine sparse_eigen(it,nbasins,sccix,dofwr,tpi,inissnap) 
!
  use interfaces
  use clusters
  use iounit 
!
  implicit none
!
  integer exstat  
  logical, INTENT(IN):: dofwr
  integer, INTENT(IN):: nbasins,tpi,sccix,inissnap
!
  type (t_scluster) it(nbasins)
!
  integer i,j,k,which_eigvals,iu,freeunit,actedim,fbix,safe_ss
!
  character(MAXSTRLEN) filename
  character(len=3) fnsffx
  character(len=17) eigfmt
!
  RTYPE tmpr,random
  RTYPE, ALLOCATABLE:: tmpvr(:,:)
  RTYPE, ALLOCATABLE:: erp_msm(:),eip_msm(:),eigvec_msm(:,:)
!
  integer nz                        ! how many non-zero values in our forward or backward matrix
  RTYPE, ALLOCATABLE:: matvals(:)   ! Non-zero values of the matrix (forward or backward time transitions)
  integer, ALLOCATABLE:: irn(:)     ! Row indexes of non-zero values of matrix
  integer, ALLOCATABLE:: jcn(:)     ! Col indexes of non-zero values of matrix
  integer matord                    ! >= 3 ; Order of matrix A
  integer LN                        ! LN >= matord. First dimension of the arrays XX, UU, and WW. Not altered by the routine
!
!
  which_eigvals = eigvalmd - 1 !for consistency with the external routines: 1 -> 0 and so on. This is IND in the orgnl hsl
!
! initializations
  if (refscc.ne.it(sconnect(inissnap,4))%inscc) then
    write(ilog,*) 'Fatal. Gotten refscc.ne.it(sconnect(inissnap,4))%insc in sparse_eigen. This is a bug.'
    call fexit()
  end if
  if (dofwr.EQV..true.) then
    fnsffx = 'FWD'
    fbix = 3
  else 
    fnsffx = 'BWD'
    fbix = 4
  end if
!
  i = 1
  k = 1
  call edit_tmat(it,nbasins,i,k)
  k = 0
  nz = 0
  do i=1,nbasins
    if (it(i)%inscc.eq.sccix) then
      k = k + 1
      do j=1,it(i)%nb
        if (it(i)%fewtsnb(j,fbix).gt.0.0) nz = nz + 1
      end do
    end if
  end do
  matord = k
  if (k.le.0) then
    write(ilog,*) 'Fatal. Called sparse_eigen(...) with a component index that has no cluster members. This is a bug.'
    call fexit()
  end if  
!
  LN = matord !same dimension of tmat(1:nsynclu) set by transition_matrix
  call check_Arnoldis(matord,numeig,arnstps,arnblks)
!
! populate input vectors
  allocate(irn(nz))
  allocate(jcn(nz))
  allocate(matvals(nz))
  allocate(tmpvr(LN,1))
!
  nz = 0
  k = 0
  do i=1,nbasins
    if (it(i)%inscc.eq.sccix) then
      k = k + 1
      if (eigvalmd.eq.2) tmpvr(k,1) = 1.0*max(sum(it(i)%wghtsnb(1:it(i)%nb,1)),sum(it(i)%wghtsnb(1:it(i)%nb,2)))
      do j=1,it(i)%nb
        if (it(i)%fewtsnb(j,fbix).gt.0.0) then
          nz = nz + 1
          irn(nz) = it(i)%ix(2)
          jcn(nz) = it(it(i)%lstnb(j))%ix(2)
          matvals(nz) = it(i)%fewtsnb(j,fbix)
        end if
      end do
    end if
  end do
  if (eigvalmd.eq.2) then 
    tmpr = 1.0/sum(tmpvr(1:LN,1))
    tmpvr(1:LN,1) = tmpr*tmpvr(1:LN,1)
  else
    do i=1,LN
      tmpvr(i,1) = 2.0*(random()-0.5)
    end do
  end if
!
! allocate output arrays
  allocate(erp_msm(numeig+1))
  allocate(eip_msm(numeig+1))
  allocate(eigvec_msm(LN,numeig+1))
  eip_msm(:) = 0.
  erp_msm(:) = 0.
  eigvec_msm(:,:) = 0.
! solve the spectral problem
  which_eigvals = eigvalmd - 1 ! local copy of eigvalmd: options 0-2 
#ifdef LINK_HSL
  call calc_eigs_msm(which_eigvals,numeig,matord,LN,nz,irn,jcn,matvals,arnstps,arnblks,arnmaxitr,&
 &doeigvec_msm,erp_msm,eip_msm,eigvec_msm,actedim,arntol,tmpvr(1:LN,1),hsl_verbosity,exstat)
#else
  write(ilog,*) 'Fatal. Executing sparse_eigen(...) despite HSL not being linked. This is a bug.'
  call fexit()
#endif
!
  actnumeig = actedim  !this is the global
!
  if ((exstat.lt.0).AND.(exstat.ne.-8).AND.(exstat.ne.-10)) then
    write(ilog,*) 'Warning. Due to previous errors from HSL routines, spectral decomposition is aborted.'
    if (actedim.ne.0) then
      write(ilog,*) 'Fatal. Gotten wrong actedim / exstat combo. This is a bug.',actedim,exstat
      call fexit()
    end if
    return
  else if (actedim.eq.0) then
    write(ilog,*) 'Warning. Due to no available converged eigenvalues, spectral decomposition is aborted.'
    return
  end if
! if possible, check steady state and normalize it
  safe_ss = 0
  if ((eigvalmd.eq.2).AND.(doeigvec_msm.EQV..true.)) then 
    write(ilog,*)
    write(ilog,'(1x,a)') 'Checking the steady state from spectral decomposition (FMCSC_EIGVAL_MD):'
    i = fbix+1
    call check_msm_ss(LN,eigvec_msm(1:LN,1),safe_ss,it,nbasins,i,refscc)  !check whether the HSL steady state is reliable
  end if
  if ((safe_ss.eq.1).AND.(exstat.eq.-10)) then
    write(ilog,*) 'Fatal. Gotten wrong combo of safe_ss and exstat. This is a bug.'
    call fexit()
  end if
!
! write to file the results
  iu = freeunit()
  filename = 'EIGENS_'//fnsffx//'.dat'
  open(unit=iu,file=trim(filename),status='unknown')
  if (doeigvec_msm.EQV..false.) then  
    do i=1,actedim
      write(iu,'(G19.11,2x,G19.11)') erp_msm(i),eip_msm(i)
    end do
  end if
  if (doeigvec_msm.EQV..true.) then
    write(eigfmt,'(A1,I9,"G",I2,".",I2,A1)') '(',matord+2,19,11,')'
    j = 1
    do while (j.le.actedim)
      if (eip_msm(j).eq.0) then
        write(iu,fmt=eigfmt) erp_msm(j),eip_msm(j),eigvec_msm(:,j)
        j = j + 1
      else if (eip_msm(j).gt.0) then
        write(iu,fmt=eigfmt) erp_msm(j),eip_msm(j),eigvec_msm(:,j)
        write(iu,fmt=eigfmt) erp_msm(j),eip_msm(j),eigvec_msm(:,j+1)
        j = j + 2
      end if
    end do
  end if
  close(unit=iu)
  if (actedim.lt.numeig) then
    iu = freeunit()
    filename = 'EIGVALS_NOTCONVRGD_'//fnsffx//'.dat'
    open(unit=iu,file=trim(filename),status='unknown')
    do i=actedim+1,numeig
      write(iu,'(G19.11,2x,G19.11)') erp_msm(i),eip_msm(i)
    end do 
    close(unit=iu)
  end if
!
  if (safe_ss.eq.1) then
    do i=1,nbasins
      if (it(i)%inscc.eq.sccix) it(i)%nodewt(fbix+1) = csccwts(refscc)*eigvec_msm(it(i)%ix(2),1)
    end do
  else if (eigvalmd.eq.2) then
    write(ilog,'(1x,a)') 'Steady state from spectral decomposition (FMCSC_EIGVAL_MD) was deemed unreliable. &
 &This means it was only written to file, but is not used otherwise.'
  end if
  write(ilog,*)
!
! propagate to common variables
  if (dofwr.EQV..true.) then  
    if (allocated(erp_msm_fwr).EQV..true.) deallocate(erp_msm_fwr)
    allocate(erp_msm_fwr(size(erp_msm)))
    erp_msm_fwr = erp_msm
    if (allocated(eip_msm_fwr).EQV..true.) deallocate(eip_msm_fwr)
    allocate(eip_msm_fwr(size(eip_msm)))
    eip_msm_fwr = eip_msm
    if (doeigvec_msm.EQV..true.) then
      if (allocated(eigvec_msm_fwr).EQV..true.) deallocate(eigvec_msm_fwr)
      allocate(eigvec_msm_fwr(size(eigvec_msm,dim=1),size(eigvec_msm,dim=2)))
    end if
  else
    if (allocated(erp_msm_bwr).EQV..true.) deallocate(erp_msm_bwr)
    allocate(erp_msm_bwr(size(erp_msm)))
    erp_msm_bwr = erp_msm
    if (allocated(eip_msm_bwr).EQV..true.) deallocate(eip_msm_bwr)
    allocate(eip_msm_bwr(size(eip_msm)))
    eip_msm_bwr = eip_msm
    if (doeigvec_msm.EQV..true.) then
      if (allocated(eigvec_msm_bwr).EQV..true.) deallocate(eigvec_msm_bwr)
      allocate(eigvec_msm_bwr(size(eigvec_msm,dim=1),size(eigvec_msm,dim=2)))
    end if
  end if
!
  if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
  if (allocated(erp_msm).EQV..true.) deallocate(erp_msm)
  if (allocated(eip_msm).EQV..true.) deallocate(eip_msm)
  if (allocated(eigvec_msm).EQV..true.) deallocate(eigvec_msm)
!
end subroutine sparse_eigen
!
!--------------------------------------------------------------------------------------------------
!
! a simple routine to check and adjust parameters for the Arnoldi method
!
subroutine check_Arnoldis(matord,neig,stps,blks)
!
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: matord,neig
  integer stps,blks
!
! checks for Arnoldi method
  if (stps.gt.floor(1.0*matord/blks)) then
    stps = max(2,floor(1.0*matord/blks))   !controlled by a keyword
    write(ilog,*) 'WARNING. The number of Arnoldi steps is too large for the actual transition matrix. &
 & Decreasing to the maximum possible allowed value, ',stps
  end if
  !
  ! checks for consistency with external routins
  if (matord.lt.3) then
    write(ilog,*) 'WARNING. The order of the transition matrix is lower than 3. External routines will crash.'
  end if
  if ((neig.lt.1).OR.(neig.gt.(matord - 2))) then
    write(ilog,*) 'WARNING. The number of requested eigenvalues is either less than 1 or bigger than the &
 &order of the transition matrix. External routines will crash.'
  end if
  if (blks.lt.1) then
    write(ilog,*) 'WARNING. The number of blocks for the Arnoldi method is smaller than 1. &
 &. External routines will crash. This is a bug.'
  end if
  if (stps.lt.2) then
    write(ilog,*) 'WARNING. The number of steps for the Arnoldi method is smaller than 2. &
 &. External routines will crash. This is a bug.'
  end if
  if ((stps*blks).gt.matord) then
    blks = max(1,floor(1.0*matord/stps))
    write(ilog,*) 'WARNING. The number steps*nblocks for the Arnoldi method is larger than the order of the matrix &
 &. Decreasing blocks to maximum possible allowe value.',blks
  end if
  if ((stps*blks).lt.min(matord,neig + 2)) then
    write(ilog,*) 'WARNING. The number steps*nblocks for the Arnoldi method is too small. External routines will crash.'
  end if
!
end
!
!-------------------------------------------------------------------------------------------------------------------------
!
#ifdef LINK_HSL
!
! the inner function for sparse matrix spectral (eigen)decomposition
! it requires passing as vectors the "i", "j", and "t_ij" values as irn, jcn, and matvals, respectively (for all nonzero t_ij)
!
! it populates OER, OEI (real and imaginary parts of eigenvalues, complex conjugates are expected to follow each other)
! and evec (for the paired complex conjugates, the complex vector of just that eigenvalue with positive imaginary part is
! printed, which is why the max im must be neigv+1; the same is true for OER OEI if the last selected one is complex)
!
subroutine calc_eigs_msm(which_eigvals,neigv,matord,LN,nz,irn,jcn,matvals,nstps,nblks,nmaxitr,doevec,OER,OEI,evec,NEV,etol,guess,&
 &vb,exstat)
!
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: neigv                    !number of eigenvalues requested, can vary to account for partial convergence
  logical, INTENT(IN):: doevec                   !whether to record eigenvectors as well
  integer, INTENT(IN):: LN                       !number of rows and columns
  integer, INTENT(IN):: nz                       !number of nonzero elements in matrix
  integer, INTENT(IN):: which_eigvals            !criterion for sorting (descending in REAL, NORM, COMPLEX)
  integer, INTENT(IN):: matord                   !order of matrix
  integer, INTENT(IN):: nstps,nblks,nmaxitr      !controls
  integer, INTENT(IN):: vb                       !verbosity level
  RTYPE, INTENT(OUT):: evec(LN,neigv+1)          !to hold eigenvectors, neigv+1 is maximum size
  RTYPE, INTENT(OUT):: OER(neigv+1),OEI(neigv+1) !arrays to hold real and imaginary part of eigenvalues
  RTYPE, INTENT(IN):: guess(LN)                  !initial guess

  integer i,j
  integer irn(nz),jcn(nz)
  RTYPE matvals(nz)
  integer exstat                       ! tell the status of INFO(1)
!
  RTYPE etol,random,eigsum             ! etol: tolerance for convergence
  RTYPE, ALLOCATABLE:: tmpvr(:)
  integer iarn
  logical notdone
  RTYPE, ALLOCATABLE:: ER(:),EI(:)     ! arrays of length nstps*nblks that will hold real and imaginary part of eigvals 
  integer ICNTL(11)                    ! On return it contains default values
  RTYPE CNTL(1)                        ! On return it contains a default value
  RTYPE, ALLOCATABLE:: XX(:,:)         ! dimensions (LN, nstps*nblks). 
                                       ! the first nblks columns of X must contain an initial estimate of the nblks basis vectors.
  RTYPE, ALLOCATABLE:: UU(:,:),WW(:,:) ! dimensions (LN,nstps*nblks).
  integer, ALLOCATABLE:: IKEEP(:)      ! length 11+nstps*nblks. 
  RTYPE, ALLOCATABLE:: RKEEP(:)        ! length 13+nstps*nblks*(2*nstps*nblks+ 10). 
  integer INFO(6)                      ! length 6 that need not be set by the user. 
  RTYPE RINFO(5)                       ! length 5 that need not be set by the user. 
  integer NEV                          ! number of eigenvalues successfully computed by EB13A/AD.
  RTYPE, ALLOCATABLE:: RES(:)          ! array of length NEV that need not be set by the user. 
!
! finish allocating 
  allocate(ER(nstps*nblks))
  allocate(EI(nstps*nblks))
  allocate(XX(1:LN,nstps*nblks))
  allocate(UU(1:LN,nstps*nblks))
  allocate(WW(1:LN,nstps*nblks))
  allocate(IKEEP(11 + nstps*nblks))
  allocate(RKEEP(13 + nstps*nblks*(2*nstps*nblks + 10)))
!
! Initialization
  notdone = .true.
  call EB13ID(ICNTL,CNTL)
  ICNTL(6) = vb     !verbose level, 6 on TEST!
  IKEEP(1) = 0
  ICNTL(11) = 999   !maximum number of iterations. I found no importance as long is > 0.
  ICNTL(5) = 99999  !maximum number of vector-matrix products. Now dynamically set in case of problems (unlikely).
  ICNTL(9) = 2      !2: the only safe option for the moment is (blocked) Arnoldi's with acceleration of starting vectors
  ICNTL(10) = 1     !1: we supply an estimate of the basis vectors (only of first one and generate the others as they do) 
  if (etol.gt.0) then
    CNTL(1) = etol  !otherwise it is computed automatically by HSL
  end if
  exstat = 0
  ICNTL(7) = 0      !0: Use norm and supply it
  RKEEP(1) = 0      !compute Frobenius norm
  do i=1,nz
    RKEEP(1) = RKEEP(1) + matvals(i)*matvals(i)
  end do
  RKEEP(1) = sqrt(RKEEP(1))  !Frobenius norm
!
! basis vectors
  XX(1:LN,1) = guess(1:LN)
  if (nblks.ge.2) then
    do j=2,nblks
      do i=1,LN
        XX(i,j) = 2.*(random()-0.5)
      end do
    end do
  end if
!
!!! First solve eigenvalue problem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  iarn = 1
  do while (notdone.EQV..true.) 
    call EB13AD(which_eigvals,matord,neigv,nblks,nstps,ER,EI,LN,XX,UU,WW,IKEEP,RKEEP,ICNTL,CNTL,INFO,RINFO)
    exstat = INFO(1) 
    if (INFO(1).eq.-9) then
      if (iarn.le.nmaxitr) then  !nmaxitr is max number of restarts set with NEIGRST
        iarn = iarn + 1
        write(ilog,'(A102)') ' Warning. EB13A failed to converge within due to a too small maximum number of iterations. &
 &Increasing.'
        ICNTL(11) = 2*ICNTL(11)  !just increase by a factor of two
        INFO(1) = 0
      else
        notdone = .false.
        write(ilog,'(A68)') ' Warning. External call to EB13A failed. No results will be written.'
        NEV = 0
      end if
    else if (INFO(1).eq.-8) then !not all is converged
      if (iarn.le.nmaxitr) then  
        iarn = iarn + 1
        write(ilog,'(A42,I9,A28,I6,A3,I6,A12,G14.6,A14)') ' Warning. EB13A failed to converge within ',nstps,' steps. &
 &Restarting attempt: ',iarn,' on',nmaxitr,'. Tolerance:',CNTL(1)
        ICNTL(10) = 1 !just restart with the estimated solution as I do not want to resize nstps
        IKEEP(1) = 0   
        ICNTL(7) = 0   
        INFO(1) = 0
      else
        notdone = .false.
        write(ilog,*) 'Accepted convergence tolerance: ',RINFO(4)
        write(ilog,'(A50,I4,A5,I6)') ' Number of eigenvalues within accepted tolerance: ',INFO(4),' on: ',neigv
        NEV = INFO(4)
      end if
    else if (INFO(1).eq.-7) then
      if (iarn.le.nmaxitr) then
        iarn = iarn + 1
        write(ilog,'(A105)') ' Warning. EB13A failed to converge within due to a too small max number of mat-vect. multipl. &
 &Increasing.'
        ICNTL(5) = 2*ICNTL(5)  
        INFO(1) = 0
        call sparse_matmul(matord,IKEEP(2),IKEEP(3),nz,jcn,irn,matvals,UU,WW,LN)  !here as written in the doc
      else
        notdone = .false.
        NEV = 0
      end if
    else if (INFO(1).lt.0) then
      notdone = .false.
      write(ilog,'(A68)') ' Warning. External call to EB13A failed. No results will be written.'
      NEV = 0
    else if (INFO(1).eq.2) then
      if (iarn.le.nmaxitr) then
        iarn = iarn + 1
        write(ilog,'(A42,I9,A28,I6,A3,I6,A12,G14.6)') ' Warning. EB13A failed to converge within ',nstps,' steps. &
 &Restarting attempt: ',iarn,' on',nmaxitr,'. Tolerance:',CNTL(1)
        ICNTL(10) = 1  
        IKEEP(1) = 0   
        ICNTL(7) = 0   
        INFO(1) = 0
      else
        notdone = .false.
        write(ilog,'(A32,I6)') 'Accepted convergence tolerance: ',RINFO(4)
        write(ilog,'(A50,I4,A5,I6)') ' Number of eigenvalues within accepted tolerance: ',INFO(4),' on: ',neigv
        NEV = INFO(4)
      end if
    else if (IKEEP(1).eq.0) then
      notdone = .false.
      NEV = INFO(4)
      exstat = INFO(1) 
    else
      call sparse_matmul(matord,IKEEP(2),IKEEP(3),nz,jcn,irn,matvals,UU,WW,LN)
    end if
  end do !do while (notdone.EQV..true.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! some checks on obtained eigenvalues, possible early return
  if (exstat.ne.INFO(1)) then
    write(ilog,*) 'Fatal. exstat.ne.INFO(1). This is a bug.',exstat,INFO(1)
    call fexit()
  end if
  write(ilog,*) 'Finally gotten ',NEV,' converged eigenvalues. User requested: ',neigv
  if (NEV.gt.neigv) then
    write(ilog,*) 'Warning. The last computed eigenvalue (at least) is complex.'
  end if
  if (NEV.eq.0) then
    write(ilog,*) 'Warning. Due to previous errors (the last one of which is reported below), the execution of the innermost &
 &routine for spectral decompostion is terminated here and no analysis that relies on these results will be perfomed.'
    call errmess_EB13AD(INFO(1))
    deallocate(ER)
    deallocate(EI)
    deallocate(XX)
    deallocate(UU)
    deallocate(WW)
    deallocate(IKEEP)
    deallocate(RKEEP)
    if (allocated(RES).EQV..true.) deallocate(RES)
    if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
    return
  end if
  if (NEV.lt.neigv) then
    write(ilog,*) 'Warning. Not all eigenvalues are converged. Extra output for unconverged eigenvalues will be written.'
    if ((exstat.ne.-8).AND.(exstat.ne.2)) then
      write(ilog,*) 'Fatal. Found some unconverged eigenvalues but with wrong exit status. This is a bug.'
      call fexit()
    end if
  end if
  do i=1,NEV !neigv
    if (EI(i).ne.0) then  
      write(ilog,*) 'Warning. Found complex eigenvalue(s). This normally indicates that detailed balance &
 &is not satisfied by the input data, which can be imposed by setting FMCSC_CADDLINKMODE to 3 or 4 if desired.'
      exit
    end if
  end do
  do i=1,NEV-1 !neigv-1
    do j=i+1,NEV
      if ((ER(i).eq.ER(j)).AND.(EI(i).eq.ER(j))) then
        write(ilog,69) j,i  
      end if
    end do
    if (EI(i).gt.0.) then ! complex eigenvalue -> complex vector across two rows
      if (abs(EI(i)+(EI(i+1))).gt.1.0e-9) then
        write(ilog,*) 'Warning. Complex eigenvalues are wrong or not in expected order in calc_eigs_msm(...). ',i,'th eigenvector&
 & may not be meaningful.'
      end if
    end if
  end do
! transfer everything (also possibly unconverged ones) to output
  do i=1,neigv
    OER(i) = ER(i)
    OEI(i) = EI(i)
  end do
!
!!! Then move to eigenvectors if required !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (doevec.EQV..true.) then
    write(ilog,*) 'Computing the eigenvectors corresponding to the converged eigenvalues...'
    if (allocated(RES).EQV..true.) deallocate(RES)
    if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
    allocate(RES(1:NEV))
    allocate(tmpvr(1:LN))
!   To compute the corresponding eigenvectors
    call EB13BD(matord,nblks,nstps,ER,EI,LN,XX,UU,WW,NEV,evec(1:LN,1:NEV),RES,IKEEP,RKEEP)
    call sparse_matmul(matord,IKEEP(2),IKEEP(3),nz,jcn,irn,matvals,UU,WW,LN)      !take the transpose, irn,jcn -> jcn,irn
    call EB13BD(matord,nblks,nstps,ER,EI,LN,XX,UU,WW,NEV,evec(1:LN,1:NEV),RES,IKEEP,RKEEP)
!   check first eigenvector
    if (which_eigvals.eq.1) then
      write(ilog,*) 'Checking first eigenvector...'
      do i=1,LN
        if (evec(i,1)*evec(1,1).lt.0) then  !change in sign
          if (etol.le.epsilon(1.d0)*10) then
            write(ilog,*) 'WARNING. Discovered a change in sign in the entries of the first eigenvector. This is an &
 &unrecoverable error. Results should be ignored.' !,sum(tmpvr)
            exstat = -10  !my hack to recognize this situation
          end if
        end if
      end do
    end if
    i = 1
    do while (i.le.NEV)
      if (EI(i).gt.0) then
        eigsum = sqrt(1./sum(evec(:,i)*evec(:,i)+evec(:,i+1)*evec(:,i+1))) 
        evec(:,i) = evec(:,i)*eigsum
        evec(:,i+1) = evec(:,i+1)*eigsum
        i = i + 2
      else
        eigsum = sqrt(1./dot_product(evec(:,i),evec(:,i)))
        evec(:,i) = evec(:,i)*eigsum
        i = i + 1
      end if
    end do
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
69 format('WARNING. Eigenvalue ',i6,' is equal to eigenvalue ',i6,'. This is likely indicative of numerical problems.')
  deallocate(ER)
  deallocate(EI)
  deallocate(XX)
  deallocate(UU)
  deallocate(WW)
  deallocate(IKEEP)
  deallocate(RKEEP)
  if (allocated(RES).EQV..true.) deallocate(RES)
  if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
!
end subroutine calc_eigs_msm
!
!----------------------------------------------------------------------------------------------------------
!
subroutine errmess_EB13AD(exstat)
!
  use iounit
!
  implicit none 
!
  integer, INTENT(IN):: exstat
!
  if (exstat.eq.0) return !nothing to say
!
  if (exstat.eq.-1) then
    write(ilog,*) 'Value of matrix rank out of range. N < 3.'
  else if (exstat.eq.-2) then
    write(ilog,*) 'Value of NBLOCK out of range. NBLOCK < 1.'
  else if (exstat.eq.-3) then
    write(ilog,*) 'Value of NUMEIG out of range.'
  else if (exstat.eq.-4) then
    write(ilog,*) 'Value of NSTEPS out of range.'
  else if (exstat.eq.-5) then
    write(ilog,*) 'Value of LN out of range.'
  else if (exstat.eq.-6) then
    write(ilog,*) 'Too small norm of input matrix.'
  else if (exstat.eq.-7) then
    write(ilog,*) 'Algorithm required more matrix-vector multiplications than specified by the control parameter ICNTL(5).'
  else if (exstat.eq.-8) then
    write(ilog,*) 'The algorithm was unable to find all the NUMEIG wanted eigenvalues with the user-supplied parameters because, &
 &at some stage, there were no points on which to construct the Chebychev ellipse. This could happen if, for example, the matrix &
 &has multiple eigenvalues.'
  else if (exstat.eq.-9) then
    write(ilog,*) 'Algorithm required more iterations than specified by the control parameter ICNTL(11).'
  else if (exstat.eq.1) then
    write(ilog,*) 'The user-supplied convergence tolerance lies outside the interval (u,1.0). Default used.'
  else if (exstat.eq.2) then
    write(ilog,*) 'Tolerance not achieved.'
  else if (exstat.eq.3) then
    write(ilog,*) 'The Arnoldi iteration stopped in fewer than NSTEPS. Some spurious zero (or almost zero) eigenvalues may be &
 &returned. The user may try choosing a smaller value for NSTEPS.'
  end if
!
end subroutine errmess_EB13AD
!
!----------------------------------------------------------------------------------------------------------
!
subroutine solve_pfold_lsys(it,nbasins,ivx,ivxsz,dofwr,tpi)
!
  use iounit
  use interfaces
  use clusters        
  use hsl_zd11_double !without this compilation/linking error
  use hsl_ma48_double !HSL_MA48_DOUBLE
  use threads, ONLY: thrdat
!
  implicit none
!
  logical, INTENT(IN):: dofwr
  integer, INTENT(IN):: nbasins,ivxsz,ivx(ivxsz,5),tpi
  type(t_scluster) it(nbasins)
!
  logical atrue,foundit
  integer iu,freeunit,i,j,k,ii,kk,iii,kkk,info,counter,nextras,refsccbu
  integer exstat
  RTYPE res(2),errr,tmpr,maxdev,tev(2,2)
  integer nclus(3)                            !how many clusters make up the folded, intermediate, and unfolded sets
  integer fbix                                !column of fewtsnb to pick dependent on dofwr
  integer ssix                                !index of nodewt to pick dependent on dofwr and prior solutions
  integer, ALLOCATABLE:: clusets(:,:)         !cluster sets: F (:,1), I (:,2), U (:,3), reverse map for I (:,4)
  integer, ALLOCATABLE:: tmpvi(:,:)           !temporary
  RTYPE, ALLOCATABLE:: tmpvr(:,:),tmpvr2(:,:),b(:),xx(:)
  integer have_ss                   ! for pfold(-): whether steady state solution already exists and is usable
  integer matord                    ! >= 3 ; Order of matrix A
  integer LN                        ! LN >= matord. First dimension of the arrays XX, UU, and WW. Not altered by the routine
  integer which_eigvals,neig        ! needed for the call to calc_eigs_msm
  integer actedim                   ! return dimensionality from spectral solver
  integer stps,blks,nmaxitr         ! num. of steps, blocks and max iter. in the Arnoldi's method for solution of spectr. prob.
  RTYPE etol                        ! tolerance for Arnoldi method
  type(ZD11_TYPE) mat               ! Pfold solver variables, mat is the vector-ized sub-matrix 
  type(MA48_CONTROL) control
  type(MA48_AINFO) ainfo
  type(MA48_FINFO) finfo
  type(MA48_SINFO) sinfo
  type(MA48_FACTORS) factors
  character(len=3) fnsffx           ! names and strings for output control
  character(len=17) eigfmt
  character(MAXSTRLEN) filename
!
! initialize 
  atrue = .true.
  nclus(1) = size(clufold)
  nclus(2) = 0
  nclus(3) = size(cluunfold)
  if (dofwr.EQV..true.) then
    fnsffx = 'fwr'
    fbix = 3
    ssix = 4
  else
    fnsffx = 'bwr'
    fbix = 4
    ssix = 5
  end if
!
  write(ilog,*) 
  if (dofwr.EQV..true.) then
    write(ilog,'(1x,a)') '-------- Setting up and solving system for (+) committor (FWD-time) -----------'
  else
    write(ilog,'(1x,a)') '-------- Setting up and solving system for (+) committor (BWD-time) -----------'
  end if
!
! some checks
  if ((nclus(1).lt.1).OR.(nclus(3).lt.1)) then
    write(ilog,*) 'Fatal. There are no eligible snapshots in lists defining folded and unfolded sets for Pfold &
 &analysis. This is a bug (check for calls to sel_ssnap(...) and read_snapsetfile(...)).'
    call fexit()
  end if
  allocate(clusets(nbasins,4))
  clusets(:,:) = 0
!
  refsccbu = refscc
  refscc = it(ivx(clufold(1),4))%inscc ! first cluster in the folded set
  write(ilog,*) 'Reference component (from first eligible snapshot constituting folded set):',refscc
!
  kk = 0
  do i=1,nclus(1)
    if (it(ivx(clufold(i),4))%inscc.ne.refscc) then
      write(ilog,*) 'Warning. Snapshot ',clufold(i),' is part of a cluster that does not belong to the reference component. &
 &Discarded a member of the folded set (in solve_pfold_lsys(...)).'
    else
      foundit = .false.
      do k=1,i-1
        if (ivx(clufold(k),4).eq.ivx(clufold(i),4)) then
          write(ilog,*) 'Warning. Snapshot ',clufold(i),' maps to the same cluster as another member of the folded set. &
 &Discarded a member of the folded set (in solve_pfold_lsys(...)).'
          foundit = .true.
          exit
        end if
      end do
      if (foundit.EQV..false.) then
        clusets(ivx(clufold(i),4),2) = -1  
      end if
    end if
  end do
  do i=1,nclus(3)
    if (it(ivx(cluunfold(i),4))%inscc.ne.refscc) then
      write(ilog,*) 'Warning. Snapshot ',cluunfold(i),' is part of a cluster that does not belong to the reference component. &
 &Discarded a member of the unfolded set (in solve_pfold_lsys(...)).'
    else
      foundit = .false.
      do k=1,i-1
        if (ivx(cluunfold(k),4).eq.ivx(cluunfold(i),4)) then
          write(ilog,*) 'Warning. Snapshot ',cluunfold(i),' maps to the same cluster as another member of the unfolded set. &
 &Discarded a member of the unfolded set (in solve_pfold_lsys(...)).'
          foundit = .true.
          exit
        end if
      end do
      if (foundit.EQV..false.) then
        if (clusets(ivx(cluunfold(i),4),2).eq.-1) then
          write(ilog,*) 'Warning. Snapshot ',cluunfold(i),' is part of a cluster that belongs to the folded set. &
 &Discarded a member of the unfolded set (in solve_pfold_lsys(...)).'
        else
          clusets(ivx(cluunfold(i),4),2) = 1
        end if
      end if
    end if
  end do
!
  nclus(:) = 0
  do i=1,nbasins
    if (clusets(i,2).eq.-1) then
      nclus(1) = nclus(1) + 1
      clusets(nclus(1),1) = i
    else if (clusets(i,2).eq.0) then
      if (it(i)%inscc.eq.refscc) then
        nclus(2) = nclus(2) + 1
        clusets(nclus(2),2) = i ! can't overtake
        clusets(i,4) = nclus(2)
      end if
    else if (clusets(i,2).eq.1) then
      nclus(3) = nclus(3) + 1
      clusets(nclus(3),3) = i
    end if
  end do
!
  if ((nclus(1).le.0).OR.(nclus(3).le.0)) then
    write(ilog,'(1x,a,2(i7,"|"),i7,a)') 'Warning. Computation of committor probabilities aborted because at least one of the 3 &
 &sets (F, I, U: ',nclus(1:3),') is empty.'
    deallocate(clusets)
    dopfold = .false.
    dopfoldminus = .false.
    if (ccfepmode.ge.8) then
      write(ilog,*) 'Consequently, aborting the generation of committor probability-based cut profiles.'
      ccfepmode = 0
    end if
    refscc = refsccbu
    return
  end if
!
!!! First solve the problem for the forward (+) committors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (pfold_report.EQV..true.) then
    iu = freeunit()
    open(unit=iu,file="fold_clus.out",status='unknown')
    do i=1,nclus(1)
      write(iu,'(I9,2x,I9)') clusets(i,1),it(clusets(i,1))%center
    end do
    close(unit=iu)
    iu = freeunit()
    open(unit=iu,file="unfold_clus.out",status='unknown')
    do i=1,nclus(3)
      write(iu,'(I9,2x,I9)') clusets(i,3),it(clusets(i,3))%center
    end do
    close(unit=iu)
  end if
!
! setting the order of the system to the number of intermediate clusters
  mat%m = nclus(2)                   !number of rows, the linear system is limited to the intermediate clusters
  mat%n = nclus(2)                   !square matrix
!
! determine the size of the sparse matrix "vector"
  counter = 0
  j = 1 ! work with fewtsnb(:,1:2)
  k = 1 ! zero out SCC-spanning transitions and normalize to resultant sub-matrices
  call edit_tmat(it,nbasins,j,k)
!
  nextras = 0
  do ii=1,nclus(2)
    i = clusets(ii,2)
    foundit = .false.
    do j=1,it(i)%nb
      if (clusets(it(i)%lstnb(j),4).eq.0) cycle ! not in set I
      if (it(i)%fewtsnb(j,fbix).gt.0.0) then
        counter = counter + 1
        if (i.eq.it(i)%lstnb(j)) then
          foundit = .true.
        end if
      end if
    end do
    if (foundit.EQV..false.) nextras = nextras + 1
  end do
!
  mat%ne = counter + nextras !all nonzero elements of matrix after(!) subtraction of identity matrix
  if (allocated(mat%val)) deallocate(mat%val)
  if (allocated(mat%row)) deallocate(mat%row)
  if (allocated(mat%col)) deallocate(mat%col)
  allocate(mat%val(mat%ne),mat%row(mat%ne),mat%col(mat%ne))
  if (allocated(b)) deallocate(b)
  allocate(b(nclus(2)))
  if (allocated(xx)) deallocate(xx) !this vector will contain tha pfold values from linear system solution
  allocate(xx(nclus(2)))
!
! after the proper allocation we perform the same cycle to assign the values
  counter = 0                       !how many non-null transitions among intermediate states
  do ii=1,nclus(2)
    i = clusets(ii,2)
    foundit = .false.
    do j=1,it(i)%nb
      if (clusets(it(i)%lstnb(j),4).eq.0) cycle ! not in set I
      if (it(i)%fewtsnb(j,fbix).gt.0.0) then
        counter = counter + 1
        mat%row(counter) = ii
        mat%col(counter) = clusets(it(i)%lstnb(j),4)
        mat%val(counter) = it(i)%fewtsnb(j,fbix) ! normalized to sub-matrix spanning refscc
        if (i.eq.it(i)%lstnb(j)) then
          mat%val(counter) =  mat%val(counter) - 1.0 ! subtract identity matrix from diagonal
          foundit = .true.
        end if
      end if
    end do
    if (foundit.EQV..false.) then
      counter = counter + 1
      mat%row(counter) = ii
      mat%col(counter) = ii
      mat%val(counter) = -1.0 ! subtract identity matrix from diagonal for missing diagonal element
    end if
  end do
!
! outputting
  if (pfold_report.EQV..true.) then
    iu = freeunit()
    filename = 'mat_pfold_'//fnsffx//'_plus.dat'
    open(unit=iu,file=trim(filename),status='unknown')
    do i=1,mat%ne
      write(iu,'(I9,2x,I9,2x,I9,2x,I9,2x,F12.9)') mat%row(i),mat%col(i), &
 &clusets(mat%row(i),2),clusets(mat%col(i),2),mat%val(i)
    end do
    close(unit=iu)
  end if
!
! right hand side of the linear system: we use the value in the transition matrix Tij only when i is intermediate and j folded
  b(:) = 0.
  do ii=1,nclus(2)
    i = clusets(ii,2)
    do j=1,it(i)%nb
      if (it(i)%fewtsnb(j,fbix).le.0.0) cycle ! nothing to subtract in any case
      if (clusets(it(i)%lstnb(j),4).eq.0) then ! not in set I
        do k=1,nclus(1)
          if (clusets(k,1).eq.it(i)%lstnb(j)) then ! search in the set F
            b(ii) = b(ii) - it(i)%fewtsnb(j,fbix) ! this edge goes into F
            exit
          end if
        end do
      end if
    end do
  end do
!
! outputting
  if (pfold_report.EQV..true.) then
    iu = freeunit()
    filename = 'rhs_pfold_'//fnsffx//'_plus.dat'
    open(unit=iu,file=trim(filename),status='unknown')
    do i=1,nclus(2)
      write(iu,'(F12.9)') b(i)
    end do
    close(unit=iu)
  end if
!
! Initialize the structures
  call MA48_INITIALIZE(factors,control)
  control%lp = ilog
  control%wp = ilog
  control%mp = ilog
  control%ldiag = hsl_verbosity                 !verbose level
  if (control%u < 1.0) control%u = 1.0          !0.01 default, increase to improve numerical accuracy they say...
  if (control%u < 1000) control%maxit = 1000    !10 default. Max number of refinement iterations
  if (control%cgce > 0.1) control%cgce = 0.1    !0.5 default. Smaller if you accept slower convergence rate
!
! Analyse and factorize
  call MA48_ANALYSE(mat,factors,control,ainfo,finfo)
  if (ainfo%flag.lt.0) then
    write(ilog,'(1x,A,I2,A)') 'Fatal. Failure of MA48_ANALYSE with error code ',ainfo%flag,'.'
    call fexit()
  else if (ainfo%flag.gt.0) then
    write(ilog,'(1x,A,I2,A)') 'Warning issued by MA48_ANALYSE with code ',ainfo%flag,'.'
  end if
  if (hsl_verbosity.ge.3) then
    write(ilog,'(1x,A60,I8)') 'Number of dropped entries in matrix analysis:',ainfo%drop
    write(ilog,'(1x,A60,I8)') 'Estimation of matrix rank from analysis:',ainfo%rank
    write(ilog,'(1x,A60,I8)') 'Structural rank from analysis:',ainfo%struc_rank
    write(ilog,'(1x,A60,I8)') 'Number of compresses in analysis:',ainfo%ncmpa
  end if
!
! Factorize
  call MA48_FACTORIZE(mat,factors,control,finfo) !fast
  if (finfo%flag.lt.0) then
    write(ilog,'(1x,A,I2,A)') 'Fatal. Failure of MA48_FACTORIZE with error code ',finfo%flag,'.'
    call fexit()
  else if (finfo%flag.gt.0) then
    write(ilog,'(1x,A,I2,A)') 'Warning issued by MA48_FACTORIZE with code ',finfo%flag,'.'
  end if
  if (hsl_verbosity.ge.3) then
    write(ilog,'(1x,A60,I8)') 'Number of dropped entries in matrix factorization:',finfo%drop
    write(ilog,'(1x,A60,I8)') 'Estimation of matrix rank from factorization:',finfo%rank
  end if
!
! Solve with iterative refinement
  call MA48_SOLVE(mat,factors,b,xx,control,sinfo,resid=res,error=errr)
  if (sinfo%flag.ne.0) then
    write(ilog,'(1x,A,I2,A)') 'Fatal. Failure of MA48_SOLVE (error in solution of linear system) with error code ',sinfo%flag,'.'
    call fexit()
  end if
  iu = freeunit()
  if (dofwr.EQV..true.) then
    if (hsl_verbosity.ge.1) then
      write(ilog,*) 'Forward-time, forward (+) committor has been computed.' 
      write(ilog,*) 'Residuals: ',res(1),res(2)
      write(ilog,*) 'Error estimation: ',errr
    end if
    open(unit=iu,file='PFOLD_PLUS_FWD.dat')
  else
    if (hsl_verbosity.ge.1) then
      write(ilog,*) 'Backward-time, forward (+) committor has been computed.'
      write(ilog,*) 'Residuals: ',res(1),res(2)
      write(ilog,*) 'Error estimation: ',errr
    end if
    open(unit=iu,file='PFOLD_PLUS_BWD.dat')
  end if
  if (allocated(pfolds).EQV..false.) then 
    allocate(pfolds(1:nbasins,4))
    pfolds = 0.
  end if
  do i=1,nclus(1)
    if (pfolds(clusets(i,1),fbix-2).ne.0.) then
      write(ilog,*) 'Fatal. Found non-zero entry in pfold for folded set. This is a bug.'
      call fexit()
    end if
    pfolds(clusets(i,1),fbix-2) = 1.
    write(iu,'(I9,2x,F12.6)') clusets(i,1),pfolds(clusets(i,1),fbix-2)
  end do
  do i=1,nclus(2)
    if (pfolds(clusets(i,2),fbix-2).ne.0.) then
      write(ilog,*) 'Fatal. Found non-zero entry in pfold for intermediate set. This is a bug.'
      call fexit()
    end if
    pfolds(clusets(i,2),fbix-2) = xx(i)
    write(iu,'(I9,2x,F12.6)') clusets(i,2),pfolds(clusets(i,2),fbix-2)
    if ( ( xx(i).gt.(1.0+1.0E-06) ).OR.( xx(i).lt.(-1.0E-06) ) ) then
      write(ilog,*) 'Warning. Pfold (+) out of range. This is either a numerical problem or a bug.', &
 &i,clusets(i,2),xx(i)
    end if
  end do
  do i=1,nclus(3)
    if (pfolds(clusets(i,3),fbix-2).ne.0.) then
      write(ilog,*) 'Fatal. Found non-zero entry in pfold (+) for unfolded set. This is a bug.'
      call fexit()
    end if
    pfolds(clusets(i,3),fbix-2) = 0.
    write(iu,'(I9,2x,F12.6)') clusets(i,3),pfolds(clusets(i,3),fbix-2)
  end do
  close(unit=iu)
!
! Clean up
  call MA48_FINALIZE(factors,control,info)
  write(ilog,'(1x,a)') '-------------------------------------------------------------------------------'
  write(ilog,*) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!! if the user requests the computation of the backward pfold too, i.e. minus (-) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (dopfoldminus.EQV..true.) then
    if (dofwr.EQV..true.) then
      write(ilog,'(1x,a)') '-------- Setting up and solving system for (-) committor (FWD-time) -----------'
    else
      write(ilog,'(1x,a)') '-------- Setting up and solving system for (-) committor (BWD-time) -----------'
    end if
!
!   the implied "inverse" propagator needs to be computed from the same transition matrix using the steady-state
!   --> we may need to recompute it
    have_ss = 0
    matord = sum(nclus(1:3))
    LN = matord
    if (refscc.eq.refsccbu) then
      do j=ssix,ssix-2,-2
        tmpr = 0.
        do i=1,nbasins
          if (it(i)%inscc.eq.refscc) tmpr = tmpr + it(i)%nodewt(j)
        end do
        if (tmpr.gt.0.0) then
          have_ss = 1
          exit
        end if  
      end do
      if (have_ss.eq.1) ssix = j
    end if
    if (have_ss.eq.1) then
      if (allocated(tmpvr2).EQV..true.) deallocate(tmpvr2)
      allocate(tmpvr2(LN,3))
      do i=1,nbasins
        if (it(i)%inscc.eq.refscc) then
          tmpvr2(it(i)%ix(2),1) = it(i)%nodewt(ssix)
        end if
      end do
!     check it
      call check_msm_ss(LN,tmpvr2(:,1),have_ss,it,nbasins,ssix,refscc)
      if (have_ss.le.0) write(*,'(1x,a)') 'Based on check, existing steady-state vector was rejected in HSL solver &
 &(FMCSC_DOPFOLD). This indicates poor settings for convergence, a numerically ill-behaved input graph, or a program &
 &inconsistency. Do not trust subsequent results unless the reasons for this warning are clear.'
    end if
    if (have_ss.le.0) then
      k = 0
      do i=1,nbasins
        if (it(i)%inscc.eq.refscc) then
          do j=1,it(i)%nb
            if (it(i)%fewtsnb(j,fbix).gt.0.0) k = k + 1
          end do
        end if
      end do
 !    use HSL Arnoldi method to retrieve it
      blks = 3
      stps = 3 
      nmaxitr = arnmaxitr
      etol = arntol
      neig = 1
      call check_Arnoldis(matord,neig,stps,blks)
      which_eigvals = 1        !hard coded, we always want rightmost ones
      if (allocated(tmpvi).EQV..true.) deallocate(tmpvi)
      allocate(tmpvi(k,2))
      if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
      allocate(tmpvr(k,1))
      if (allocated(tmpvr2).EQV..true.) deallocate(tmpvr2)
      allocate(tmpvr2(LN,3))
      k = 0
      do i=1,nbasins
        if (it(i)%inscc.eq.refscc) then
          tmpvr2(it(i)%ix(2),1) = it(i)%nodewt(1)  !nodeweight always ok here as which_eigvals hard-coded to 1
          do j=1,it(i)%nb
            if (it(i)%fewtsnb(j,fbix).gt.0.0) then
              k = k + 1
              tmpvi(k,1) = it(i)%ix(2)
              tmpvi(k,2) = it(it(i)%lstnb(j))%ix(2)
              tmpvr(k,1) = it(i)%fewtsnb(j,fbix)
            end if
          end do
        end if
      end do
      call calc_eigs_msm(which_eigvals,neig,matord,LN,k,tmpvi(:,1),tmpvi(:,2),tmpvr(:,1),stps,blks,nmaxitr,atrue,&
 &tev(:,1),tev(:,2),tmpvr2(:,1:2),actedim,etol,tmpvr2(:,1),hsl_verbosity,exstat)
      if ((exstat.lt.0).AND.(exstat.ne.-8).AND.(exstat.ne.-10)) then
        write(ilog,*) 'Warning. Due to previous errors from HSL routines, spectral decomposition is aborted.'
        write(ilog,*) 'Consequently, all pfold (-) analyses will be skipped.'
        dopfoldminus = .false.
        deallocate(clusets)
        if (ccfepmode.eq.9) then
          write(ilog,*) 'Consequently, aborting the generation of pfold(-)-based cut profiles.'
          ccfepmode = 0
        else if (ccfepmode.eq.10) then
          write(ilog,*) 'Consequently, aborting the generation of pfold(-)-based cut profiles.'
          ccfepmode = 8
        end if
        refscc = refsccbu
        if (allocated(tmpvi).EQV..true.) deallocate(tmpvi)
        if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
        if (allocated(tmpvr2).EQV..true.) deallocate(tmpvr2)
        return
      else if (actedim.eq.0) then
        write(ilog,*) 'Warning. Due to no available converged eigenvalues, spectral decomposition is aborted.'
                write(ilog,*) 'Consequently, all pfold (-) analyses will be skipped.'
        dopfoldminus = .false.
        deallocate(clusets)
        if (ccfepmode.eq.9) then
          write(ilog,*) 'Consequently, aborting the generation of pfold(-)-based cut profiles.'
          ccfepmode = 0
        else if (ccfepmode.eq.10) then
          write(ilog,*) 'Consequently, aborting the generation of pfold(-)-based cut profiles.'
          ccfepmode = 8
        end if
        refscc = refsccbu
        if (allocated(tmpvi).EQV..true.) deallocate(tmpvi)
        if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
        if (allocated(tmpvr2).EQV..true.) deallocate(tmpvr2)
        return
      end if
      if (pfold_report.EQV..true.) then
        iu = freeunit()
        filename = 'ss_pfold_'//fnsffx//'.dat'
        open(unit=iu,file=trim(filename),status='unknown')
        write(eigfmt,'(A1,I9,"G",I2,".",I2,A1)') '(',matord+2,14,6,')'
        do j=1,actedim
          write(iu,fmt=eigfmt) tev(j,1),tev(j,2),tmpvr2(:,j)
        end do
        close(unit=iu)
      end if
      deallocate(tmpvi)
      deallocate(tmpvr)
    end if
!     
    if (have_ss.le.0) then
!     check it
      call check_msm_ss(LN,tmpvr2(:,1),have_ss,it,nbasins,ssix,refscc)
      if (have_ss.le.0) write(*,'(1x,a)') 'Based on check, newly computed steady-state vector was rejected in HSL solver &
 &(FMCSC_DOPFOLD). This indicates poor settings for the Arnoldi method or a numerically ill-behaved input graph. Do &
 &not trust subsequent results unless the reasons for this warning are clear.'
    end if
    if ((have_ss.eq.1).AND.(ssix.gt.3)) then
!     a new LA solution: store it
      do i=1,nbasins
        if (it(i)%inscc.eq.refscc) it(i)%nodewt(ssix) = csccwts(refscc)*tmpvr2(it(i)%ix(2),1)
      end do
    end if
!
!   in the setup of data structures for the linear system, the sizes remain the same (F,I,U,fbix unchanged)!
    counter = 0
    do ii=1,nclus(2)
      i = clusets(ii,2) ! only set I
      iii = it(i)%ix(2)
      foundit = .false.
      do j=1,it(i)%nb
        if (clusets(it(i)%lstnb(j),4).eq.0) cycle ! not in set I
        kkk = it(it(i)%lstnb(j))%ix(2)
        if (it(i)%fewtsnb(j,fbix).gt.0.0) then
          counter = counter + 1
          mat%row(counter) = clusets(it(i)%lstnb(j),4) ! transposed relative to (+) case
          mat%col(counter) = ii
          mat%val(counter) = tmpvr2(iii,1)*it(i)%fewtsnb(j,fbix)/tmpvr2(kkk,1) ! implied "backward" propagator
          if (i.eq.it(i)%lstnb(j)) then
            foundit = .true.
            mat%val(counter) = mat%val(counter) - 1.0
          end if
        end if
      end do
      if (foundit.EQV..false.) then
        counter = counter + 1
        mat%row(counter) = ii
        mat%col(counter) = ii
        mat%val(counter) = -1.0 ! subtract identity matrix from diagonal for missing diagonal element
      end if
    end do
!   right hand side of the linear system: we use the value in the transition matrix Tij only when j is intermediate and i unfolded
    b(:) = 0.
    do ii=1,nclus(3)
      i = clusets(ii,3) ! just set U
      iii = it(i)%ix(2) 
      do j=1,it(i)%nb
        if (it(i)%fewtsnb(j,fbix).le.0.0) cycle ! nothing to subtract in any case
        if (clusets(it(i)%lstnb(j),4).gt.0) then ! is in set I
          kkk = it(it(i)%lstnb(j))%ix(2)
          b(clusets(it(i)%lstnb(j),4)) = b(clusets(it(i)%lstnb(j),4)) - tmpvr2(iii,1)*it(i)%fewtsnb(j,fbix)/tmpvr2(kkk,1)
        end if
      end do
    end do
!   check transition matrix integrity
    tmpvr2(:,3) = 0.
    res(1) = HUGE(res(1))
    res(2) = 0.
    do i=1,nbasins
      if (it(i)%inscc.ne.refscc) cycle
      iii = it(i)%ix(2)
      do j=1,it(i)%nb
        if (it(i)%fewtsnb(j,fbix).gt.0.) then
          kkk = it(it(i)%lstnb(j))%ix(2)
          tmpr = tmpvr2(iii,1)*it(i)%fewtsnb(j,fbix)/tmpvr2(kkk,1)
          tmpvr2(kkk,3) = tmpvr2(kkk,3) + tmpr
          res(1) = min(res(1),tmpr)
          res(2) = max(res(2),tmpr)
        end if
      end do
    end do
    maxdev = 0.
    do i=1,nbasins
      if (it(i)%inscc.ne.refscc) cycle
      maxdev = max(abs(1.0-tmpvr2(it(i)%ix(2),3)),maxdev)
      if (abs(1.0-tmpvr2(it(i)%ix(2),3)).gt.1E-03) then
        write(ilog,69) i,tmpvr2(it(i)%ix(2),3)
      end if
    end do
 69 format( 'WARNING. Row ',i6,' of (-) transition matrix does not sum up to 1 by far, but to: ',g14.6,/,& 
 &'. This can likely generate arbitrarily large errors in backward (-) committor, for example pfold (-) > 1 or < 0.')
!
!   write some info
    write(ilog,'(A63,g16.8)') 'Largest deviation from 1 of row sum of (-) transition matrix: ',maxdev
    write(ilog,'(A63,g16.8)') 'Minimum not-null value in (-) transition matrix: ',res(1)
    write(ilog,'(A63,g16.8)') 'Maximum value in (-) transition matrix: ',res(2)
!
!   outputting
    if (pfold_report.EQV..true.) then
      iu = freeunit()
      filename = 'mat_pfold_'//fnsffx//'_minus.dat'
      open(unit=iu,file=trim(filename),status='unknown')
      do i=1,mat%ne
        write(iu,'(I9,2x,I9,2x,I9,2x,I9,2x,F12.9)') mat%row(i),mat%col(i), &
   &clusets(mat%row(i),2),clusets(mat%col(i),2),mat%val(i)
      end do
      close(unit=iu)
    end if
!
!   outputting
    if (pfold_report.EQV..true.) then
      iu = freeunit()
      filename = 'rhs_pfold_'//fnsffx//'_minus.dat'
      open(unit=iu,file=trim(filename),status='unknown')
      do i=1,nclus(2)
        write(iu,'(F12.9)') b(i)
      end do
      close(unit=iu)
    end if
!
!   Initialize the structures
    call MA48_INITIALIZE(factors,control)
    control%lp = ilog
    control%wp = ilog
    control%mp = ilog
    control%ldiag = 2 !verbose level
    if (control%u < 1.0) control%u = 1.0      
    if (control%u < 1000) control%maxit = 1000
    if (control%cgce > 0.1) control%cgce = 0.1  
!
!   Analyse and factorize
!   Analyse and factorize
    call MA48_ANALYSE(mat,factors,control,ainfo,finfo)
    if (ainfo%flag.lt.0) then
      write(ilog,'(1x,A,I2,A)') 'Failure of MA48_ANALYSE with error code ',ainfo%flag,'.'
      call fexit()
    else if (ainfo%flag.gt.0) then
      write(ilog,'(1x,A,I2,A)') 'Warning issued by MA48_ANALYSE with code ',ainfo%flag,'.'
    end if
    if (hsl_verbosity.ge.3) then
      write(ilog,'(1x,A60,I8)') 'Number of dropped entries in matrix analysis:',ainfo%drop
      write(ilog,'(1x,A60,I8)') 'Estimation of matrix rank from analysis:',ainfo%rank
      write(ilog,'(1x,A60,I8)') 'Structural rank from analysis:',ainfo%struc_rank
      write(ilog,'(1x,A60,I8)') 'Number of compresses in analysis:',ainfo%ncmpa
    end if
!
!   Factorize
    call MA48_FACTORIZE(mat,factors,control,finfo) !fast
    if (finfo%flag.lt.0) then
      write(ilog,'(1x,A,I2,A)') 'Failure of MA48_FACTORIZE with error code ',finfo%flag,'.'
      call fexit()
    else if (finfo%flag.gt.0) then
      write(ilog,'(1x,A,I2,A)') 'Warning issued by MA48_FACTORIZE with code ',finfo%flag,'.'
    end if
    if (hsl_verbosity.ge.3) then
      write(ilog,'(1x,A60,I8)') 'Number of dropped entries in matrix factorization:',finfo%drop
      write(ilog,'(1x,A60,I8)') 'Estimation of matrix rank from factorization:',finfo%rank
    end if
!
!   Solve with iterative refinement
    call MA48_SOLVE(mat,factors,b,xx,control,sinfo,resid=res,error=errr)
    if (sinfo%flag.ne.0) then
      write(ilog,'(1x,A,I2,A)') 'Fatal. Failure of MA48_SOLVE (error in solution of linear system) with error code ',sinfo%flag,'.'
      call fexit()
    end if
    iu = freeunit()
    if (dofwr.EQV..true.) then
      if (hsl_verbosity.ge.1) then
        write(ilog,*) 'Forward-time, backward (-) committor has been computed.'
        write(ilog,*) 'Residuals: ',res(1),res(2)
        write(ilog,*) 'Error estimation: ',errr
      end if
      open(unit=iu,file='PFOLD_MINUS_FWD.dat')
    else
      if (hsl_verbosity.ge.1) then
        write(ilog,*) 'Backward-time, backward (-) committor has been computed.'
        write(ilog,*) 'Residuals: ',res(1),res(2)
        write(ilog,*) 'Error estimation: ',errr
      end if
      open(unit=iu,file='PFOLD_MINUS_BWD.dat')
    end if
    do i=1,nclus(1)
      if (pfolds(clusets(i,1),fbix).ne.0.) then
        write(ilog,*) 'Fatal. Found non-zero entry in pfold for folded set. This is a bug.'
         call fexit()
      end if
      pfolds(clusets(i,1),fbix) = 0.
      write(iu,'(I9,2x,F12.6)') clusets(i,1),pfolds(clusets(i,1),fbix)
    end do
    do i=1,nclus(2)
      if (pfolds(clusets(i,2),fbix).ne.0.) then
        write(ilog,*) 'Fatal. Found non-zero entry in pfold (-) for intermediate set. This is a bug.'
        call fexit()
      end if
      pfolds(clusets(i,2),fbix) = xx(i)
      write(iu,'(I9,2x,F12.6)') clusets(i,2),pfolds(clusets(i,2),fbix)
      if ( ( xx(i).gt.(1.0+1.0E-06) ).OR.( xx(i).lt.(-1.0E-06) ) ) then
        write(ilog,*) 'Warning. Pfold (-) out of range. This is either a numerical problem or a bug.', &
 &i,clusets(i,2),xx(i)
      end if
    end do
    do i=1,nclus(3)
      if (pfolds(clusets(i,3),fbix).ne.0.) then
        write(ilog,*) 'Fatal. Found non-zero entry in pfold for folded set. This is a bug.'
        call fexit()
      end if
      pfolds(clusets(i,3),fbix) = 1.
      write(iu,'(I9,2x,F12.6)') clusets(i,3),pfolds(clusets(i,3),fbix)
    end do
    close(unit=iu)
!
!   Clean up
    call MA48_FINALIZE(factors,control,info)
  end if !dopfoldminus
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  write(ilog,*) '-------------------------------------------------------------------------------'
!
  if (allocated(mat%val).EQV..true.) deallocate(mat%val)
  if (allocated(mat%row).EQV..true.) deallocate(mat%row)
  if (allocated(mat%col).EQV..true.) deallocate(mat%col)
  if (allocated(b).EQV..true.) deallocate(b)
  if (allocated(xx).EQV..true.) deallocate(xx)
  if (allocated(tmpvi).EQV..true.) deallocate(tmpvi)
  if (allocated(tmpvr).EQV..true.) deallocate(tmpvr)
  if (allocated(tmpvr2).EQV..true.) deallocate(tmpvr2)
!
  refscc = refsccbu
!
end subroutine solve_pfold_lsys
!
#endif
!--------------------------------------------------------------------------------------------
!
