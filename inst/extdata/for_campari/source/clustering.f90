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
! CONTRIBUTIONS: Nicolas Bloechliger, Jiri Vymetal, Marco Bacci            !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-------------------------------------------------------------------------
!
subroutine store_for_clustering(tpi)
!
  use sequen
  use forces
  use clusters
  use fyoc
  use polypep
  use atoms
  use system
  use aminos
  use molecule
  use movesets
  use iounit
  use zmatrix
  use math
  use mcsums
  use interfaces
#ifdef ENABLE_MPI
  use mpistuff
  use mpi
#endif
  use pdb, ONLY: pdb_rmol
  use threads, ONLY: thrdat
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer kk,ttc,ttc2,i,rs,imol,lmol,kidx,cdofsz,shf
  logical afalse
  RTYPE d2v,mmtw,clutmp(calcsz),tmpzt,getztor,tvec(3)
  RTYPE, ALLOCATABLE:: testvals(:,:)
#ifdef ENABLE_MPI
  integer clumpiavgtag,mstatus(MPI_STATUS_SIZE),masterrank,ierr
#endif
  integer(KIND=8) tts(3)
!
  afalse = .false.
!
#ifdef ENABLE_MPI
  masterrank = 0
  clumpiavgtag = 355
#endif
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  cstored = cstored + 1
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
! torsional
  if ((cdis_crit.ge.1).AND.(cdis_crit.le.4)) then
    if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4)) then
      if (cchangeweights.eq.0) then ! weights are not replaced later -> inertial masses are relevant
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
!       may require backup 
        if ((dyn_integrator_ops(1).gt.0).AND.(use_dyn.EQV..true.).AND.(pdb_analyze.EQV..false.).AND.(fycxyz.eq.1)) then
          do imol=1,nmol
            dc_di(imol)%olddat(:,3) = dc_di(imol)%olddat(:,2)
          end do
        end if
        ! update com and populate olddat(:,2)
        do imol=1,nmol
          call update_rigidm(imol)
        end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call cart2int_I(tpi)
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    cdofsz = calcsz/2
    if (cdis_crit.eq.1) cdofsz = calcsz
    if (cdis_crit.eq.4) cdofsz = calcsz/3
!
    clutmp(:) = 0.0
!   all torsions
    ttc = 0
    ttc2 = 1
    do rs=1,nseq
      imol = molofrs(rs)
!     omega
      if (wline(rs).gt.0) then
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(wnr(rs).gt.0)) mmtw = dc_di(imol)%olddat(wnr(rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = omega(rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = omega(rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(omega(rs)/RADIAN)
            clutmp(2*ttc2) = cos(omega(rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(omega(rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(omega(rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     phi
      if (fline(rs).gt.0) then
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(fnr(rs).gt.0)) mmtw = dc_di(imol)%olddat(fnr(rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = phi(rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = phi(rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(phi(rs)/RADIAN)
            clutmp(2*ttc2) = cos(phi(rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(phi(rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(phi(rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     psi
      if (yline(rs).gt.0) then
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(ynr(rs).gt.0)) mmtw = dc_di(imol)%olddat(ynr(rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (yline2(rs).gt.0) then
            tmpzt = ztor(yline2(rs))
          else
            tmpzt = psi(rs)
          end if 
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = tmpzt
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = tmpzt
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(tmpzt/RADIAN)
            clutmp(2*ttc2) = cos(tmpzt/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(tmpzt/RADIAN)
            clutmp(3*ttc2-1) = cos(tmpzt/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     nucleic acid bb angles
      do kk = 1,nnucs(rs)
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(nucsnr(kk,rs).gt.0)) mmtw = dc_di(imol)%olddat(nucsnr(kk,rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = nucs(kk,rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = nucs(kk,rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(nucs(kk,rs)/RADIAN)
            clutmp(2*ttc2) = cos(nucs(kk,rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(nucs(kk,rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(nucs(kk,rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end do
      if (ttc2.gt.cdofsz) exit
!     sugar bond in nucleotides
      if (seqpolty(rs).eq.'N') then
        kidx = nucsline(6,rs)
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(nnucs(rs).gt.1)) then
          if (nucsnr(nnucs(rs)-1,rs).gt.0) mmtw = dc_di(imol)%olddat(nucsnr(nnucs(rs)-1,rs),2)
        end if
!       note the cheating on the inertial mass for the sugar bond
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = ztor(kidx)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = ztor(kidx)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(ztor(kidx)/RADIAN)
            clutmp(2*ttc2) = cos(ztor(kidx)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(ztor(kidx)/RADIAN)
            clutmp(3*ttc2-1) = cos(ztor(kidx)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end if
!     chi angles
      do kk = 1,nchi(rs)
        ttc = ttc + 1
        mmtw = 1.0
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).AND.(chinr(kk,rs).gt.0)) mmtw = dc_di(imol)%olddat(chinr(kk,rs),2)
        if (cdofset(ttc2,1).eq.ttc) then
          if (cdis_crit.eq.1) then
            clutmp(ttc2) = chi(kk,rs)
          else if (cdis_crit.eq.2) then
            clutmp(2*ttc2-1) = chi(kk,rs)
            clutmp(2*ttc2) = mmtw
          else if (cdis_crit.eq.3) then
            clutmp(2*ttc2-1) = sin(chi(kk,rs)/RADIAN)
            clutmp(2*ttc2) = cos(chi(kk,rs)/RADIAN)
          else if (cdis_crit.eq.4) then
            clutmp(3*ttc2-2) = sin(chi(kk,rs)/RADIAN)
            clutmp(3*ttc2-1) = cos(chi(kk,rs)/RADIAN)
            clutmp(3*ttc2) = mmtw
          end if
          ttc2 = ttc2 + 1
          if (ttc2.gt.cdofsz) exit
        end if
      end do
      if (ttc2.gt.cdofsz) exit
!     completion for crosslinks
      if (disulf(rs).gt.0) then
        if (crosslink(crlk_idx(rs))%itstype.le.2) then ! disulfide
          ttc = ttc + 1
          if (cdofset(ttc2,1).eq.ttc) then
            tmpzt = getztor(cai(rs),at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf))
            if (cdis_crit.eq.1) then
              clutmp(ttc2) = tmpzt
            else if (cdis_crit.eq.2) then
              clutmp(2*ttc2-1) = tmpzt
              clutmp(2*ttc2) = 1.0 ! meaningless value
            else if (cdis_crit.eq.3) then
              clutmp(2*ttc2-1) = sin(tmpzt/RADIAN)
              clutmp(2*ttc2) = cos(tmpzt/RADIAN)
            else if (cdis_crit.eq.4) then
              clutmp(3*ttc2-2) = sin(tmpzt/RADIAN)
              clutmp(3*ttc2-1) = cos(tmpzt/RADIAN)
              clutmp(3*ttc2) = 1.0 ! meaningless value
            end if
            ttc2 = ttc2 + 1
            if (ttc2.gt.cdofsz) exit
          end if
          if (disulf(rs).gt.rs) then
            ttc = ttc + 1
            if (cdofset(ttc2,1).eq.ttc) then
              tmpzt = getztor(at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf),at(disulf(rs))%sc(2-shf))
              if (cdis_crit.eq.1) then
                clutmp(ttc2) = tmpzt
              else if (cdis_crit.eq.2) then
                clutmp(2*ttc2-1) = tmpzt
                clutmp(2*ttc2) = 1.0 ! meaningless value
              else if (cdis_crit.eq.3) then
                clutmp(2*ttc2-1) = sin(tmpzt/RADIAN)
                clutmp(2*ttc2) = cos(tmpzt/RADIAN)
              else if (cdis_crit.eq.4) then
                clutmp(3*ttc2-2) = sin(tmpzt/RADIAN)
                clutmp(3*ttc2-1) = cos(tmpzt/RADIAN)
                clutmp(3*ttc2) = 1.0 ! meaningless value
              end if
              ttc2 = ttc2 + 1
              if (ttc2.gt.cdofsz) exit
            end if
          end if
        end if
      end if
!     possibly completion for unsupported
      if (((dyn_integrator_ops(10).eq.3).OR.((dyn_integrator_ops(10).eq.2).AND.(seqtyp(rs).ne.26)).OR.&
                                            ((dyn_integrator_ops(10).eq.1).AND.(seqtyp(rs).eq.26))).AND.&
 &        ((unslst%nr.gt.0).OR.(unklst%nr.gt.0))) then
        do kk=max(atmol(molofrs(rs),1)+3,at(rs)%bb(1)),at(rs)%bb(1)+at(rs)%na-1
          if (izrot(kk)%alsz.gt.0) then
            if (natlst%nr.gt.0) then
              call binary_search(natlst%nr,natlst%idx(1:natlst%nr),kk,kidx)
              if (natlst%idx(max(1,min(natlst%nr,kidx))).eq.kk) cycle
            end if
            if ((kk.eq.fline(rs)).OR.(kk.eq.yline(rs)).OR.(kk.eq.fline2(rs)).OR.(kk.eq.yline2(rs))) cycle
            if (((dyn_integrator_ops(10).ge.2).AND.(seqtyp(rs).ne.26)).OR.& ! must be in unslst
                ((dyn_integrator_ops(10).ne.2).AND.(seqtyp(rs).eq.26))) then ! must be in unklst
              ttc = ttc + 1 
              if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4)) then
                if (dc_di(imol)%recurs(izrot(kk)%treevs(4),3).gt.0) then
                  mmtw = dc_di(imol)%olddat(dc_di(imol)%recurs(izrot(kk)%treevs(4),3),2)
                end if
              end if
              if (cdofset(ttc2,1).eq.ttc) then
                if (cdis_crit.eq.1) then
                  clutmp(ttc2) = ztor(kk)
                else if (cdis_crit.eq.2) then
                  clutmp(2*ttc2-1) = ztor(kk)
                  clutmp(2*ttc2) = mmtw
                else if (cdis_crit.eq.3) then
                  clutmp(2*ttc2-1) = sin(ztor(kk)/RADIAN)
                  clutmp(2*ttc2) = cos(ztor(kk)/RADIAN)
                else if (cdis_crit.eq.4) then
                  clutmp(3*ttc2-2) = sin(ztor(kk)/RADIAN)
                  clutmp(3*ttc2-1) = cos(ztor(kk)/RADIAN)
                  clutmp(3*ttc2) = mmtw
                end if
                ttc2 = ttc2 + 1
                if (ttc2.gt.cdofsz) exit
              end if
            end if
          end if
        end do
        if (ttc2.gt.cdofsz) exit
      end if
    end do
    if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4)) then
      if (cchangeweights.gt.0) then
!       may require restore
        if ((dyn_integrator_ops(1).gt.0).AND.(use_dyn.EQV..true.).AND.(pdb_analyze.EQV..false.).AND.(fycxyz.eq.1)) then
          do imol=1,nmol
            dc_di(imol)%olddat(:,2) = dc_di(imol)%olddat(:,3)
          end do
        end if
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
! xyz atom set
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
!   respect ref mol shifts
    lmol = 0
    tvec(:) = 0.0
    do i=1,clstsz/3
      imol = molofrs(atmres(cdofset(i,1)))
      if (pdb_rmol.gt.0) then
        if (lmol.ne.imol) then
          call shift_bound3(pdb_rmol,imol,tvec)
        end if
        lmol = imol
      end if
      clutmp(3*i-2) = x(cdofset(i,1)) + tvec(1)
      clutmp(3*i-1) = y(cdofset(i,1)) + tvec(2)
      clutmp(3*i) = z(cdofset(i,1)) + tvec(3)
    end do
    if (cdis_crit.eq.10) clutmp((clstsz+1):calcsz) = 1.0 ! currently no locally adaptive weight is collected live
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
! internal distance vector
  else if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
    if (cdistransform.eq.3) then
      tmpzt = 0.5*PI/cdistrans_params(1)
    end if
    do i=1,clstsz
      call dis_bound(cdofset(i,1),cdofset(i,2),d2v)
      clutmp(i) = sqrt(d2v)
      if (cdistransform.eq.0) then
!       linear
        clutmp(i) = sqrt(d2v)
      else if (cdistransform.eq.1) then
!       sigmoidal
        clutmp(i) = 1.0 - 1.0/(1.0 + exp(-(sqrt(d2v)-cdistrans_params(1))/cdistrans_params(2)))
      else if (cdistransform.eq.2) then
!       hyperbolic
        clutmp(i) = (sqrt(d2v)+cdistrans_params(1))**(-1.0/cdistrans_params(2))
      else if (cdistransform.eq.3) then
!       sine
        d2v = sqrt(d2v)
        if (d2v.le.cdistrans_params(1)) then
          clutmp(i) = sin(d2v*tmpzt)
        else
          clutmp(i) = 1.0
        end if
      end if
    end do
    if (cdis_crit.eq.9) clutmp((clstsz+1):calcsz) = 1.0 ! currently no locally adaptive weight is collected live
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end if
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) then
    if (myrank.eq.masterrank) then
      cludata(:,cstored) = clutmp(:)
      do i=1,mpi_nodes-1
        call MPI_RECV(clutmp,calcsz,MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.clumpiavgtag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected ',clumpiavgtag,'.'
          call fexit()
        end if
        cludata(:,i*cmaxsnaps+cstored) = clutmp(:)
      end do
    else
      call MPI_SEND(clutmp,calcsz,MPI_RTYPE,masterrank,clumpiavgtag,MPI_COMM_WORLD,ierr)
    end if
  else
    cludata(:,cstored) = clutmp(:)
  end if
#else
  cludata(:,cstored) = clutmp(:)
#endif
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  if ((cstored.eq.1000).AND.(cmode.ge.3)) then
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.((use_MPIMultiSeed.EQV..true.).OR.(myrank.ne.0))) return
#endif
#ifdef ENABLE_THREADS
    allocate(testvals(1000,ceiling(1000./(1.0*thrdat%maxn))+1))
    if (tpi.le.1) call System_Clock(count=tts(1))
    ttc2 = 0
    do i=max(tpi,1),1000,max(thrdat%maxn,1)
      ttc2 = ttc2 + 1
      do ttc=1,1000
        call snap_to_snap_d(testvals(ttc,ttc2),i,ttc)
      end do 
    end do
!$OMP BARRIER
#else
    allocate(testvals(1000,1000))
    call System_Clock(count=tts(1))
    ttc2 = 0
    do i=1,1000
      ttc2 = ttc2 + 1
      do ttc=1,1000
        call snap_to_snap_d(testvals(ttc,ttc2),i,ttc)
      end do
    end do
#endif
    if (tpi.le.1) then
      call System_Clock(count=tts(2))
 345 format('Performance test on distance evaluations yielded a rate of ',g12.4,' millions/sec.')
 346 format('This suggests a low estimate for the post-processing time of ',g12.4,' s (',g12.4,' h, ',g12.4,' d)')
      tts(3) = 0
      if (cmode.eq.3) then
        tts(3) = 2*cmaxsnaps*(cmaxsnaps-1)
      else if (cmode.eq.5) then
        tts(3) = 100*cmaxsnaps
      else if ((cmode.eq.4).AND.(cprogindex.eq.1)) then
        tts(3) = 3*cmaxsnaps*(cmaxsnaps-1)
      else if ((cmode.eq.4).AND.(cprogindex.eq.2)) then
        tts(3) = cmaxsnaps*(100 + cprogindrmax*max(2,nint(log(1.0*cmaxsnaps)))) ! all are very approximate of course
      end if
      write(ilog,345) thrdat%rate/(1.0*(tts(2)-tts(1)))
      if (tts(3).gt.0) then
        d2v = (1.0e-6*tts(3))/(thrdat%rate/(1.0*(tts(2)-tts(1))))
        write(ilog,346) d2v,d2v/3600,d2v/86400
      end if
    end if
  end if
!
end
!
!--------------------------------------------------------------------------------------------------
!
! allow centering, variance normalization, and smoothing prior to clustering
!
subroutine preprocess_cludata(modei)
!
  use clusters
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: modei
!
  RTYPE cmean,icst,cvar
  integer effdd,dd
  RTYPE, ALLOCATABLE:: cvx(:),cvx2(:)
!
  if (modei.gt.0) then
    if (csmoothie.ge.(cstored/2-1)) then
      write(ilog,*) 'Warning. Selected smoothing window size is too large (FMCSC_CSMOOTHORDER ',csmoothie,').'
      csmoothie = cstored/2-1
      write(ilog,*) 'Value is adjusted to (roughly) half the data set length (',csmoothie,').'
    end if
  end if
!
  if (cprepmode.le.0) return
!
  if ((cdis_crit.eq.1).OR.(cdis_crit.eq.2)) then
    if (modei.gt.0) then
      write(ilog,*) 'Warning. Disabling data preprocessing for structural clustering due to use of &
 &explicit dihedral angles (periodic quantities with fixed range centered at 0.0).'
    end if
    return
  end if
!
  if (modei.gt.0) then
    if ((cdis_crit.eq.3).OR.(cdis_crit.eq.4)) then
      write(ilog,*) 'Warning. Data preprocessing for structural clustering may not be particularly &
 &meaningful when performed on sine and cosine terms of dihedral angles.'
    end if
    if ((align_for_clustering.EQV..true.).AND.((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10))) then
      if (cprepmode.gt.1) then
        write(ilog,*) 'Warning. Disabling data preprocessing for structural clustering except centering due to use of &
 &coordinate RMSD with alignment. Use a prealigned trajectory instead.'
        cprepmode = 1
      end if
    end if
  end if
!
  icst = 1.0/(1.0*cstored)
!
  effdd = clstsz
  if (cdis_crit.eq.4) effdd = 3*effdd/2
  if (cprepmode.ge.3) then
    allocate(cvx(cstored))
    allocate(cvx2(cstored))
  end if
  do dd=1,effdd
    if ((cdis_crit.eq.4).AND.(mod(dd,3).eq.0)) cycle ! skip weight dimension
    if (cprepmode.ne.3) then
      cmean = icst*sum(cludata(dd,1:cstored))
      cludata(dd,1:cstored) = cludata(dd,1:cstored) - cmean
    end if
    if ((cprepmode.eq.2).OR.(cprepmode.eq.5)) then
      cvar = 1.0/sqrt(sum(cludata(dd,1:cstored)**2)/(1.0*(cstored-1)))
      cludata(dd,1:cstored) = cvar*cludata(dd,1:cstored)
    end if
    if (cprepmode.ge.3) then
      cvx(:) = cludata(dd,1:cstored)
      call vsmooth(cvx,cstored,csmoothie,cvx2)
      cludata(dd,1:cstored) = cvx2(:)
    end if
  end do
!
  if (allocated(cvx).EQV..true.) deallocate(cvx)
  if (allocated(cvx2).EQV..true.) deallocate(cvx2)
!
end
!
!---------------------------------------------------------------------------------------------------
!
! allows weights to be repopulated with different types of information
! 0: leave as is
! 1: inverse standard deviation 
! 2: static weights using ACF at fixed lag
! 3: square root of inverse std dev and ACF at fixed lag
! 4: rate of transitions across global mean
! 5: same as 4 with smoothing for weights generation
! 6: sqrt(2*4)
! 7: sqrt(2*5)
! 8: rate of transitions across separators defined by histogram minima
! 9: same as 8 with smoothing for weights generation
!
subroutine repopulate_weights()
!
  use clusters
  use iounit
  use math
#ifdef ENABLE_MPI
  use mpistuff, ONLY: use_MPIAVG,mpi_nodes
#endif
!
  implicit none
!
  integer dd,k,i,j,ii,ntrans,ncvs,nbins,nbinsm10,nbinsblk,effdd,wpos,lastk
  RTYPE cmean,cvar,icst,cacf,cm1,cm2
  RTYPE, ALLOCATABLE:: cvx(:),cvx2(:),cvw(:),cvlst(:)
#ifdef ENABLE_MPI
!  RTYPE, ALLOCATABLE:: piggsy(:,:)
#endif
  integer, ALLOCATABLE:: cvi(:),ixlst(:)
!
 677 format('Mean weights per dimension: ',1000000(g12.5,1x))
!
  if (cchangeweights.eq.0) return
!
  if (cstored.le.4) then
    write(ilog,*) 'Warning. Refusing to adjust weights for clustering due to insufficient data length.'
    return
  end if
!
  if (cchangeweights.eq.2) then
    if (cwacftau.ge.(cstored/2-1)) then
      write(ilog,*) 'Warning. Selected lag time for ACF weights is too large (FMCSC_CLAGTIME ',cwacftau,').'
      cwacftau = cstored/2-1
      write(ilog,*) 'Value is adjusted to roughly half the data set length (',cwacftau,').'
    end if
  end if
!
  nbins = 100
  nbinsm10 = nbins - 10
  nbinsblk = 3
  if (cchangeweights.ge.8) then
    allocate(cvi(nbins))
    allocate(cvlst(nbinsm10))
    allocate(ixlst(nbinsm10))
  end if
  allocate(cvx(cstored))
!
  if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
    allocate(cvw(clstsz))
    cvw(:) = 0.0
  end if
  if ((cchangeweights.eq.1).OR.(cchangeweights.eq.3).OR.&
 &    (cchangeweights.eq.5).OR.(cchangeweights.eq.7).OR.(cchangeweights.eq.9)) allocate(cvx2(cstored))
!
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    if (cwwindowsz.ge.cstored) then
      write(ilog,*) 'Warning. Selected window size for locally adaptive weights is too large (FMCSC_CWINDOWSIZE ',cwwindowsz,').'
      cwwindowsz = cstored/2
      write(ilog,*) 'Value is adjusted to roughly half the data set length (',cwwindowsz,').'
    end if
  end if
  if (cchangeweights.eq.0) return
!
! locally adaptive weights of all types
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    lastk = 2
    if (cdis_crit.ge.9) lastk = 1
!   to avoid having to rewrite the same code several times, we pull the selector inside
    ii = cwwindowsz/2
    icst = 1.0/(1.0*cstored)
    effdd = clstsz
    if (cdis_crit.eq.4) effdd = effdd/2
    do dd=1,effdd
      do k=1,lastk
        if (cdis_crit.eq.4) then
          cvx(:) = cludata(3*dd-3+k,1:cstored)
          wpos = 3*dd
        else if (cdis_crit.eq.2) then
          if (k.eq.1) cvx(:) = sin(cludata(2*dd-1,1:cstored)/RADIAN)
          if (k.eq.2) cvx(:) = cos(cludata(2*dd-1,1:cstored)/RADIAN)
          wpos = 2*dd
        else
          cvx(:) = cludata(dd,1:cstored) 
          wpos = dd + clstsz
        end if
!       if so desired, populate temp vector of static weights with ACF values
        if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
#ifdef ENABLE_MPI
! concatenated trajectories imply a blockiness that can be analyzed
!          if (use_MPIAVG.EQV..true.) then 
!            allocate(piggsy(mpi_nodes,2))
!            nbinsblk = mpi_nodes
!            call vmsf_block(cvx,cstored,nbinsblk,piggsy(:,1),piggsy(:,2))
! 333 format(80(f9.3,1x))
!             write(*,333) piggsy(:,2),sqrt(sum((cvx(401:600)-(sum(cvx(401:600))/200.))**2)/200.)
!            cm1 = sum(piggsy(:,1))/(1.0*mpi_nodes) 
!            cvw(dd) = sum(piggsy(:,2)) / sqrt(sum((piggsy(:,1)-cm1)**2))
!            deallocate(piggsy)
!          else
            call vacf_fixtau(cvx,cstored,cwacftau,cacf)
            cvw(dd) = max(cvw(dd),cacf)
!          end if
#else
          call vacf_fixtau(cvx,cstored,cwacftau,cacf)
          cvw(dd) = max(cvw(dd),cacf)
#endif
        end if
!       do nothing more for static ACF weights
        if (cchangeweights.eq.2) then
!       locally adaptive measure based on variance possibly combined with ACF
        else if ((cchangeweights.eq.1).OR.(cchangeweights.eq.3)) then
          i = cwwindowsz
          call vmsf_window(cvx,cstored,i,cvx2)
          do i=1,cstored
            if (k.eq.1) cludata(wpos,i) = max(0.0,cvx2(i))
            if (k.eq.2) cludata(wpos,i) = cludata(wpos,i) + max(0.0,cvx2(i))
            if ((k.eq.lastk).AND.(cludata(wpos,i).gt.0.0)) cludata(wpos,i) = sqrt(1.0/(cludata(wpos,i)))
          end do
!       locally adaptive measure based on crossings of the global mean, possibly estimated from smoothed data,
!       possibly combined with ACF
        else if ((cchangeweights.ge.4).AND.(cchangeweights.le.7)) then
          if ((cchangeweights.eq.5).OR.(cchangeweights.eq.7)) then
            call vsmooth(cvx,cstored,csmoothie,cvx2)
            cvx(:) = cvx2(:)
          end if
          cmean = icst*sum(cvx)
          ntrans = 0
          do i=2,cwwindowsz
            if (((cvx(i-1).gt.cmean).AND.(cvx(i).lt.cmean)).OR.((cvx(i-1).lt.cmean).AND.(cvx(i).gt.cmean))) then
              ntrans = ntrans + 1
            end if
          end do
!          write(*,*) dd,cmean,ntrans
          if (k.eq.1) cludata(wpos,ii) = 1.0*ntrans+cdynwbuf
          if (k.eq.2) cludata(wpos,ii) = 1.0*min(cludata(wpos,ii),ntrans+cdynwbuf)
          do i=cwwindowsz+1,cstored
            if (((cvx(i-1).gt.cmean).AND.(cvx(i).lt.cmean)).OR.((cvx(i-1).lt.cmean).AND.(cvx(i).gt.cmean))) ntrans = ntrans + 1
            if (((cvx(i-cwwindowsz).gt.cmean).AND.(cvx(i-cwwindowsz+1).lt.cmean)).OR.&
 &              ((cvx(i-cwwindowsz).lt.cmean).AND.(cvx(i-cwwindowsz+1).gt.cmean))) ntrans = ntrans - 1
!
            if (k.eq.1) cludata(wpos,i-ii) = 1.0*ntrans+cdynwbuf
            if (k.eq.2) cludata(wpos,i-ii) = 1.0*min(cludata(wpos,i-ii),ntrans+cdynwbuf)
          end do
!       locally adaptive measure based on passes of separators, possibly of smoothed data
        else if ((cchangeweights.eq.8).OR.(cchangeweights.eq.9)) then
          cvi(:) = 0
          call vautohist(cvx,cstored,nbins,cvi,cm1,cm2)
          ncvs = 0
          ixlst(:) = 0
          call vminima_int(cvi,nbins,ncvs,ixlst,nbinsm10,nbinsblk)
          do i=1,ncvs
            cvlst(i) = cm1 + (ixlst(i)-0.5)*cm2
          end do
          if (ncvs.le.0) then
            ncvs = 1
            cvlst(ncvs) = icst*sum(cvx)
          end if
          if (cchangeweights.eq.9) then
            call vsmooth(cvx,cstored,csmoothie,cvx2)
            cvx(:) = cvx2(:)
          end if
!         we need to penalize low transition rates based on separators that are within a low likelihood tail
!         cm1 should correspond to the maximal value of the smaller fractional population for any considered separator
          cm1 = 0.0
          do j=1,ncvs
            if (ixlst(j).le.0) then
              cm1 = max(0.5,cm1)
            else
              cm1 = max(cm1,min(sum(cvi(1:ixlst(j)))/(1.0*cstored),sum(cvi((ixlst(j)+1):100))/(1.0*cstored)))
            end if
          end do
          ntrans = 0
          do i=2,cwwindowsz
            do j=1,ncvs
              if (((cvx(i-1).gt.cvlst(j)).AND.(cvx(i).lt.cvlst(j))).OR.((cvx(i-1).lt.cvlst(j)).AND.(cvx(i).gt.cvlst(j)))) then
                ntrans = ntrans + 1
                exit
              end if
            end do
          end do
          if (k.eq.1) cludata(wpos,ii) = (1.0*ntrans+cdynwbuf)/cm1
          if (k.eq.2) cludata(wpos,ii) = min(cludata(wpos,ii),(1.0*ntrans+cdynwbuf)/cm1)
          do i=cwwindowsz+1,cstored
            do j=1,ncvs
              if (((cvx(i-1).gt.cvlst(j)).AND.(cvx(i).lt.cvlst(j))).OR.((cvx(i-1).lt.cvlst(j)).AND.(cvx(i).gt.cvlst(j)))) then
                ntrans = ntrans + 1
                exit
              end if
            end do
            do j=1,ncvs
              if (((cvx(i-cwwindowsz).gt.cvlst(j)).AND.(cvx(i-cwwindowsz+1).lt.cvlst(j))).OR.&
 &                ((cvx(i-cwwindowsz).lt.cvlst(j)).AND.(cvx(i-cwwindowsz+1).gt.cvlst(j)))) then
                ntrans = ntrans - 1
                exit
              end if
            end do
            if (k.eq.1) cludata(wpos,i-ii) = (1.0*ntrans+cdynwbuf)/cm1
            if (k.eq.2) cludata(wpos,i-ii) = min(cludata(wpos,i-ii),(1.0*ntrans+cdynwbuf)/cm1)
          end do
        else
          write(ilog,*) 'Fatal. Weights of type ',cchangeweights,' are not supported at the moment. &
 &This is likely to be an omission bug.'
          call fexit()
        end if
      end do
    end do
    if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
      if (sum(cvw(1:effdd)).le.0.0) then
        cvw(1:effdd) = 1.0 ! if all weights are negligible, turn to flat weights
      end if
    end if
    do dd=1,effdd
      if (cdis_crit.eq.4) then
        wpos = 3*dd
      else if (cdis_crit.eq.2) then
        wpos = 2*dd
      else
        wpos = dd+clstsz
      end if
!     complete those data with access to incomplete windows
      if (cchangeweights.gt.3) then  
        cludata(wpos,1:(ii-1)) = cludata(wpos,ii)
        cludata(wpos,(cstored-ii+1):cstored) = cludata(wpos,(cstored-ii))
      end if
      if ((cchangeweights.ge.4).AND.(cchangeweights.le.9)) then
        cludata(wpos,1:cstored) = 1.0/cludata(wpos,1:cstored)
      end if
      if ((cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
        cludata(wpos,1:cstored) = sqrt(cvw(dd)*cludata(wpos,1:cstored))
      else if (cchangeweights.eq.2) then
        cludata(wpos,1:cstored) = cvw(dd)
      end if
!      write(*,*) dd,icst*sum(cludata(wpos,1:cstored)),&
! &icst*sum((cludata(wpos,1:cstored)-icst*sum(cludata(wpos,1:cstored)))**2)
    end do
!    do i=1,cstored
!      write(0,677) cludata((clstsz+1):(2*clstsz),i)!,cludata(calcsz,i)
!    end do
! static weights are only for interatomic distances at the moment
  else if (cdis_crit.eq.8) then
    icst = 1.0/(1.0*cstored)
    if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        call vacf_fixtau(cvx,cstored,cwacftau,cacf)
        cvw(dd) = max(cvw(dd),cacf)
      end do
      if (sum(cvw(1:clstsz)).le.0.0) then
        cvw(1:clstsz) = 1.0 ! if all weights are negligible, turn to flat weights
      end if
    end if
    if ((cchangeweights.eq.1).OR.(cchangeweights.eq.3)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        cmean = icst*sum(cvx)
        cvar = icst*sum((cvx(:)- cmean)**2)
        if (cvar.le.0.0) then
          cl_imvec(dd) = 0.0
          cycle
        end if
        cl_imvec(dd) = sqrt(1.0/cvar)
        if (cchangeweights.eq.3) cl_imvec(dd) = sqrt(cl_imvec(dd)*cvw(dd))
      end do
!   based on ACF
    else if (cchangeweights.eq.2) then
      cl_imvec(1:clstsz) = cvw(1:clstsz)
!   based on passes of global mean
    else if ((cchangeweights.ge.4).AND.(cchangeweights.le.7)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        if ((cchangeweights.eq.5).OR.(cchangeweights.eq.7)) then
          call vsmooth(cvx,cstored,csmoothie,cvx2)
          cvx(:) = cvx2(:)
        end if
        cmean = icst*sum(cvx)
        ntrans = 0
        do i=2,cstored
          if (((cvx(i-1).gt.cmean).AND.(cvx(i).lt.cmean)).OR.((cvx(i-1).lt.cmean).AND.(cvx(i).gt.cmean))) then
            ntrans = ntrans + 1
          end if
        end do
        cl_imvec(dd) = 1.0/(1.0*(ntrans+cdynwbuf))
        if ((cchangeweights.eq.6).OR.(cchangeweights.eq.7)) cl_imvec(dd) = sqrt(cl_imvec(dd)*cvw(dd))
      end do
!   based on passes of separators, possibly of smoothed data
    else if ((cchangeweights.eq.8).OR.(cchangeweights.eq.9)) then
      do dd=1,clstsz
        cvx(:) = cludata(dd,1:cstored)
        cvi(:) = 0
        call vautohist(cvx,cstored,nbins,cvi,cm1,cm2)
        ncvs = 0
        ixlst(:) = 0
        call vminima_int(cvi,nbins,ncvs,ixlst,nbinsm10,nbinsblk)
        do i=1,ncvs
          cvlst(i) = cm1 + (ixlst(i)-0.5)*cm2
        end do
        if (ncvs.le.0) then
          ncvs = 1
          cvlst(ncvs) = icst*sum(cvx)
        end if
        if (cchangeweights.eq.9) then
          call vsmooth(cvx,cstored,csmoothie,cvx2)
          cvx(:) = cvx2(:)
        end if
        ntrans = 0
        do i=2,cstored
          do k=1,ncvs
            if (((cvx(i-1).gt.cvlst(k)).AND.(cvx(i).lt.cvlst(k))).OR.((cvx(i-1).lt.cvlst(k)).AND.(cvx(i).gt.cvlst(k)))) then
              ntrans = ntrans + 1
              exit
            end if
          end do
          if (i.eq.cwwindowsz) ntrans = 0
        end do
!       we need to penalize low transition rates based on separators that are within a low likelihood tail
!       cm1 should correspond to the maximal value of the smaller fractional population for any considered separator
        cm1 = 0.0
        do k=1,ncvs
          if (ixlst(k).le.0) then
            cm1 = max(0.5,cm1)
          else
            cm1 = max(cm1,min(sum(cvi(1:ixlst(k)))/(1.0*cstored),sum(cvi((ixlst(k)+1):100))/(1.0*cstored)))
          end if
        end do
        cl_imvec(dd) = cm1/(1.0*(ntrans+cdynwbuf))
      end do
    else
      write(ilog,*) 'Fatal. Weights of type ',cchangeweights,' are not supported for interatomic distances. &
 &This is an omission bug.'
      call fexit()
    end if
  end if
!
  if (allocated(cvi).EQV..true.) deallocate(cvi)
  if (allocated(cvx).EQV..true.) deallocate(cvx)
  if (allocated(cvx2).EQV..true.) deallocate(cvx2)
  if (allocated(ixlst).EQV..true.) deallocate(ixlst)
  if (allocated(cvlst).EQV..true.) deallocate(cvlst)
  if (allocated(cvw).EQV..true.) deallocate(cvw)
!
end
!
!------------------------------------------------------------------------------------
!
! dynamic weights should be normalized per snapshot so that different snapshots 
! have equivalent "net" sizes
! this function will also print out summary information regarding weights (incl. static ones)
! it should be called after any calls to "reduce_cludimensionality"
!
subroutine normalize_weights()
!
  use clusters
  use iounit
!
  implicit none
!
  integer i,dd,wpos,effdd
  RTYPE normer,icst
  RTYPE, ALLOCATABLE:: cvw(:)
  logical have_warned
!
  have_warned = .false.
  icst = 1.0/(1.0*cstored)
!
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    effdd = clstsz
    if (cdis_crit.eq.4) effdd = effdd/2
    do i=1,cstored
      normer = 0.
      do dd=1,effdd
        if (cdis_crit.eq.4) then
          wpos = 3*dd
        else if (cdis_crit.eq.2) then
          wpos = 2*dd
        else
          wpos = dd+clstsz
        end if
        normer = normer + cludata(wpos,i)
      end do
      if (normer.le.0.) then
        if (have_warned.EQV..false.) then
          write(ilog,*) 'Warning. At least one snapshot has dynamic weights that are all zero. Adjusted to flat weights.'
          have_warned = .true.
        end if
        do dd=1,effdd
          if (cdis_crit.eq.4) then
            wpos = 3*dd
          else if (cdis_crit.eq.2) then
            wpos = 2*dd
          else
            wpos = dd+clstsz
          end if
          cludata(wpos,i) = 1.0/(1.0*effdd)
        end do
      else  
!      write(*,*) i,normer
        do dd=1,effdd
          if (cdis_crit.eq.4) then
            wpos = 3*dd
          else if (cdis_crit.eq.2) then
            wpos = 2*dd
          else
            wpos = dd+clstsz
          end if
          cludata(wpos,i) = cludata(wpos,i)/normer
        end do
      end if
    end do
  end if
!
 677 format('Mean weights per dimension: ',1000000(g12.5,1x))
 678 format('Static weights per dimension: ',1000000(g12.5,1x))
!
! print out information about mean weights
  if (cdis_crit.eq.8) then
    normer = sum(cl_imvec(:))
    if (normer.le.0.0) then ! no variance in data
      cl_imvec(:) = sqrt(1.0/(1.0*clstsz))
    else
      cl_imvec(:) = sqrt(cl_imvec(:)/normer)
    end if
    write(ilog,678) cl_imvec(1:clstsz)**2 ! print the square as weights are already adjusted (square root) for use in prod. form
  else if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    if (allocated(cvw).EQV..false.) allocate(cvw(effdd))
    do dd=1,effdd
      if (cdis_crit.eq.4) then
        wpos = 3*dd
      else if (cdis_crit.eq.2) then
        wpos = 2*dd
      else
        wpos = dd+clstsz
      end if
      cvw(dd) = icst*sum(cludata(wpos,1:cstored))
    end do
    write(ilog,677) cvw(1:effdd)/sum(cvw(1:effdd))
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
subroutine do_clustering(tpi)
!
  use clusters
  use iounit
  use threads, ONLY: thrdat
#ifdef ENABLE_MPI
  use mpistuff, ONLY: use_MPIAVG,use_MPIMultiSeed
  use mcsums, ONLY: nstep
  use system, ONLY: pdb_analyze
#endif
#ifdef LINK_NETCDF
  use ncdm
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer exitcode,azero,atwo,aone
  integer(KIND=8) t1,t2
  logical atrue
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
#endif
!
  exitcode = 0
  azero = 0
  atrue = .true.
  atwo = 2
  aone = 1
!
  if (cstored.le.0) return

#ifdef ENABLE_THREADS
  if (tpi.le.0) then
    if (OMP_IN_PARALLEL().EQV..true.) then
      write(ilog,*) 'Fatal. When using multi-threaded code, do_clustering(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
      call fexit()
    end if
  end if
#endif
!
! for a PIGS analysis run, we drop out here and simply call ASMaster
#ifdef ENABLE_MPI
  if ((pdb_analyze.EQV..true.).AND.(use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
    exitcode = nstep
    call MPI_ASMaster(exitcode,tpi)
    return
!  for a PIGS simulation run, the last stretch of data is simply ignored
  else if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
    return
  end if
#endif  
!
! preprocess data if so desired
  if (((cchangeweights.gt.0).OR.(cprepmode.gt.0)).AND.(tpi.le.1)) write(ilog,*)
  call System_Clock(t1)
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  call preprocess_cludata(cmode)
  call repopulate_weights()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  call System_Clock(t2)
  if ((cchangeweights.gt.0).OR.(cprepmode.gt.0)) then
    if (tpi.le.1) write(ilog,*) 'Time elapsed for data preprocessing for structural clustering: ',(t2-t1)/thrdat%rate, ' [s]'
  end if
!
  if (pcamode.gt.1) then
#ifdef LINK_LAPACK
    if ((cdis_crit.ge.3)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      if ((pcamode.ne.3).AND.(pcamode.ne.5)) then
        call reduce_cludimensionality()
      end if
      call normalize_weights()
      call get_principal_components(calcsz,cstored,clstsz)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    else if (reduced_dim_clustering.gt.0) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      call reduce_cludimensionality()
      call normalize_weights()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    end if
#endif
  else
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
    call reduce_cludimensionality()
    call normalize_weights()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end if
#ifdef LINK_NETCDF
  if (cfeature_dump.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.ge.9)) then
      call ncdm_set_cludatamasks(atrue,cdis_crit,ncdm_featvals_includata,clstsz,ncdm_dynwghts_includata) 
    else
      call ncdm_set_cludatamasks(atrue,cdis_crit,ncdm_featvals_includata,clstsz)
    end if
    call ncdm_write_ncfl(aone,atrue) !write down all the features (sin and cos are different and no extra dim for them in cdf)
#ifdef ENABLE_THREADS
!$OMP END SINGLE NOWAIT
#endif
  end if
#endif
!
  if ((cmode.eq.1).OR.(cmode.eq.2)) then
    call System_Clock(t1)
    call leader_clustering(cmode,nstruccls,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    call System_Clock(t2)
    if (tpi.le.1) write(ilog,*) 'Time elapsed for clustering and graph-based analyses: ',(t2-t1)/thrdat%rate, ' [s]'
  else if (cmode.eq.5) then
    call System_Clock(t1)
#ifdef ENABLE_THREADS
    call birch_clustering_threads(cmode,nstruccls,tpi)
!$OMP BARRIER
#else
    call birch_clustering(cmode,nstruccls)
#endif
    call System_Clock(t2)
    if (tpi.le.1) write(ilog,*) 'Time elapsed for clustering and graph-based analyses: ',(t2-t1)/thrdat%rate, ' [s]'
  else if (cmode.eq.6) then
    call System_Clock(t1)
    call file_clustering(cmode,nstruccls,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    call System_Clock(t2)
    if (tpi.le.1) write(ilog,*) 'Time elapsed for clustering and graph-based analyses: ',(t2-t1)/thrdat%rate, ' [s]'
  else
    if (cmode.eq.4) then
      call System_Clock(t1)
      if (cprogindex.eq.1) then
#ifdef LINK_NETCDF
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        if (read_nbl_from_nc.EQV..true.) then
          call read_nbl_nc()
        else
          call gennbl_for_clustering()
        end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#else
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        call gennbl_for_clustering()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#endif
        call System_Clock(t2)
        if (tpi.le.1) write(ilog,*) 'Time elapsed for neighbor list creation/reading: ',(t2-t1)/thrdat%rate, ' [s]'
        call System_Clock(t1)
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        call gen_MST_from_nbl()
        call do_prog_index()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call System_Clock(t2)
        if (tpi.le.1) write(ilog,*) 'Time elapsed for generating progress index data: ',(t2-t1)/thrdat%rate, ' [s]'
      else
#ifdef ENABLE_THREADS
        call birch_clustering_threads(cmode,nstruccls,tpi)
!$OMP BARRIER
#else
        call birch_clustering(cmode,nstruccls)
#endif
        call System_Clock(t2)
        if (tpi.le.1) write(ilog,*) 'Time elapsed for clustering and graph-based analyses: ',(t2-t1)/thrdat%rate, ' [s]'
        call System_Clock(t1)
        call gen_MST_from_treeclustering(tpi)
        call System_Clock(t2)
        if (tpi.le.1) write(ilog,*) 'Time elapsed for MST building: ',(t2-t1)/thrdat%rate, ' [s]'
        call System_Clock(t1)
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        call do_prog_index()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call System_Clock(t2)
        if (tpi.le.1) write(ilog,*) 'Time elapsed for generating progress index data: ',(t2-t1)/thrdat%rate, ' [s]'
      end if
    else if (cmode.eq.3) then
      call System_Clock(t1)
#ifdef LINK_NETCDF
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      if (read_nbl_from_nc.EQV..true.) then
        call read_nbl_nc()
      else
        call gennbl_for_clustering()
      end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#else
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      call gennbl_for_clustering()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#endif
      call System_Clock(t2)
      if (tpi.le.1) write(ilog,*) 'Time elapsed for neighbor list creation/reading: ',(t2-t1)/thrdat%rate, ' [s]'
      call System_Clock(t1)
      call hierarchical_clustering(nstruccls,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      call System_Clock(t2)
      if (tpi.le.1) write(ilog,*) 'Time elapsed for clustering and graph-based analyses: ',(t2-t1)/thrdat%rate, ' [s]'
    else
      if (tpi.le.1) write(ilog,*) 'Fatal. Encountered unknown clustering mode in do_clustering(...). This&
 & is a bug.'
      call fexit()
    end if
    if (exitcode.ne.0) then
      if (tpi.le.1) write(ilog,*) 'Fatal. Clustering algorithm in do_clustering(...) exited with an error. This&
 & is either a bug or a corrupt setup (mode ',cmode,').'
      call fexit()
    end if
  end if
!
end
!
!--------------------------------------------------------------------------------------
!
subroutine gennbl_for_clustering()
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer i,j,ii,jj,kk,ll,k,l,mi,mj,maxalcsz,nsets,nzeros,nclalcsz
  integer(KIND=8) testcnt,testcnt2
  RTYPE rmsdtmp,rmsdctr,bucr,rmsdclcl
  RTYPE, ALLOCATABLE:: buf1(:),buf2(:),maxrads(:)
  integer, ALLOCATABLE:: bui1(:),bui2(:)
  logical atrue
!
  atrue = .true.
!
  allocate(cnblst(cstored))
  cnblst(:)%nbs = 0
  cnblst(:)%alsz = 4
  maxalcsz = 4
  do i=1,cstored
    allocate(cnblst(i)%idx(cnblst(i)%alsz))
    allocate(cnblst(i)%dis(cnblst(i)%alsz))
  end do
!
! local screening
  bucr = cradius 
  cradius = cmaxrad
  write(ilog,*)
  write(ilog,*) 'Now using truncated leader algorithm for pre-screening in neighbor list generation ...'
  nclalcsz = 10
  allocate(scluster(nclalcsz))
  scluster(:)%nmbrs = 0
  scluster(:)%alsz = 0
  scluster(:)%nb = 0
  scluster(:)%nchildren = 0
  scluster(:)%nbalsz = 0
  scluster(:)%chalsz = 0
  k = 0
  do i=1,cstored
    ii = -1
    do j=k,max(k-500,1),-1
      call snap_to_cluster_d(rmsdtmp,scluster(j),i)
      if (rmsdtmp.lt.cradius) then
        call cluster_addsnap(scluster(j),i,rmsdtmp)
        ii = j
        exit
      end if
    end do
    if (ii.eq.-1) then
      k = k + 1
      if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
      call cluster_addsnap(scluster(k),i,rmsdtmp)
    end if
  end do
  nsets = k
!  call leader_clustering(1,nsets)
  write(ilog,*) '... done with initial cluster generation for neighbor list.'
  write(ilog,*)
!
  write(ilog,*) 'Now finding maximum radius of all ',nsets,' identified clusters ...'
  allocate(maxrads(nsets))
  do i=1,nsets
    maxrads(i) = 0.0
    do j=1,scluster(i)%nmbrs
      call snap_to_cluster_d(rmsdtmp,scluster(i),scluster(i)%snaps(j))
      if (rmsdtmp.gt.maxrads(i)) then
        maxrads(i) = rmsdtmp
      end if
    end do
  end do
  write(ilog,*) '... done.'
  write(ilog,*)
!
! now compare all blocks to each other (the slowest part) taking advantage of information
! generated previously (otherwise intractable) 
  ii = 0
  testcnt = 0
  testcnt2 = 0
  write(ilog,*) 'Now computing cutoff-assisted neighbor list ...'
  do i=1,nsets
    ii = ii + 1
    jj = ii 
    do kk=1,scluster(i)%nmbrs
      k = scluster(i)%snaps(kk)
      do ll=kk+1,scluster(i)%nmbrs
        l = scluster(i)%snaps(ll)
        call snap_to_snap_d(rmsdtmp,k,l)
        testcnt = testcnt + 1
        if (rmsdtmp.lt.chardcut) then
!         we'll store these redundantly
          testcnt2 = testcnt2 + 1
          cnblst(k)%nbs = cnblst(k)%nbs + 1
          if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(k)
          if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
          cnblst(k)%idx(cnblst(k)%nbs) = l
          cnblst(k)%dis(cnblst(k)%nbs) = rmsdtmp
          cnblst(l)%nbs = cnblst(l)%nbs + 1
          if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(l)
          if (cnblst(l)%alsz.gt.maxalcsz) maxalcsz = cnblst(l)%alsz
          cnblst(l)%idx(cnblst(l)%nbs) = k
          cnblst(l)%dis(cnblst(l)%nbs) = rmsdtmp
        end if
      end do
    end do
    do j=i+1,nsets
      jj = jj + 1
      call cluster_to_cluster_d(rmsdclcl,scluster(i),scluster(j))
      testcnt = testcnt + 1
      if ((rmsdclcl - maxrads(i) - maxrads(j)).gt.chardcut) cycle
      if ((scluster(i)%nmbrs.eq.1).AND.(scluster(j)%nmbrs.eq.1)) then
        if (rmsdclcl.lt.chardcut) then
          k = scluster(i)%snaps(1)
          l = scluster(j)%snaps(1)
!         we'll store these redundantly
          testcnt2 = testcnt2 + 1
          cnblst(k)%nbs = cnblst(k)%nbs + 1
          if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(k)
          if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
          cnblst(k)%idx(cnblst(k)%nbs) = l
          cnblst(k)%dis(cnblst(k)%nbs) = rmsdclcl
          cnblst(l)%nbs = cnblst(l)%nbs + 1
          if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(l)
          if (cnblst(l)%alsz.gt.maxalcsz) maxalcsz = cnblst(l)%alsz
          cnblst(l)%idx(cnblst(l)%nbs) = k
          cnblst(l)%dis(cnblst(l)%nbs) = rmsdclcl
        end if        
        cycle
      end if
!      if (maxrads(j).le.maxrads(i)) then
      if (scluster(i)%nmbrs.le.scluster(j)%nmbrs) then
        mj = j
        mi = i
      else
        mj = i
        mi = j
      end if
      do kk=1,scluster(mi)%nmbrs
        k = scluster(mi)%snaps(kk)
        call snap_to_cluster_d(rmsdctr,scluster(mj),k)
        testcnt = testcnt + 1
        if ((rmsdctr - maxrads(mj)).gt.chardcut) cycle
        do ll=1,scluster(mj)%nmbrs
          l = scluster(mj)%snaps(ll)
          if (scluster(mj)%nmbrs.eq.1) then
            rmsdtmp = rmsdctr
          else
            call snap_to_snap_d(rmsdtmp,k,l)
            testcnt = testcnt + 1
          end if
          if (rmsdtmp.lt.chardcut) then
!           we'll store these redundantly
            testcnt2 = testcnt2 + 1
            cnblst(k)%nbs = cnblst(k)%nbs + 1
            if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(k)
            if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
            cnblst(k)%idx(cnblst(k)%nbs) = l
            cnblst(k)%dis(cnblst(k)%nbs) = rmsdtmp
            cnblst(l)%nbs = cnblst(l)%nbs + 1
            if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(l)
            if (cnblst(l)%alsz.gt.maxalcsz) maxalcsz = cnblst(l)%alsz
            cnblst(l)%idx(cnblst(l)%nbs) = k
            cnblst(l)%dis(cnblst(l)%nbs) = rmsdtmp
          end if
        end do
      end do
    end do
  end do
  write(ilog,*) '... done after computing ',(100.0*testcnt)/(0.5*cstored*(cstored-1)),'% of &
 &possible terms with ',(100.0*testcnt2)/(1.0*testcnt),'% successful.'
  write(ilog,*)
!
  write(ilog,*) 'Now sorting neighbor lists ...'
 33 format(200000000(i6,1x))
 34 format(200000000(f10.3,1x))
  allocate(buf1(maxalcsz))
  allocate(buf2(maxalcsz))
  allocate(bui1(maxalcsz))
  allocate(bui2(maxalcsz))
  nzeros = 0
  do i=1,cstored
    if (cnblst(i)%nbs.le.0) then
      nzeros = nzeros + 1
!      write(ilog,*) 'Warning. Snapshot # ',i,' is without a neighbor (similar) structure. This &
! &may cause the clustering algorithm to crash or misbehave otherwise.'
    else
      if (cnblst(i)%nbs.gt.1) then
        buf1(1:cnblst(i)%nbs) = cnblst(i)%dis(1:cnblst(i)%nbs)
        do ii=1,cnblst(i)%nbs
          bui1(ii) = ii
        end do
        ii = 1
        jj = cnblst(i)%nbs
        call merge_sort(ldim=cnblst(i)%nbs,up=atrue,list=buf1(1:cnblst(i)%nbs),olist=buf2(1:cnblst(i)%nbs),ilo=ii,ihi=jj,&
 &                  idxmap=bui1(1:cnblst(i)%nbs),olist2=bui2(1:cnblst(i)%nbs))
        buf1(1:cnblst(i)%nbs) = cnblst(i)%dis(1:cnblst(i)%nbs)
        bui1(1:cnblst(i)%nbs) = cnblst(i)%idx(1:cnblst(i)%nbs)
        do j=1,cnblst(i)%nbs
          cnblst(i)%dis(j) = buf1(bui2(j))
          cnblst(i)%idx(j) = bui1(bui2(j))
        end do
      end if
      allocate(cnblst(i)%tagged(cnblst(i)%nbs))
      cnblst(i)%tagged(:) = .false.
    end if
  end do
  if (nzeros.gt.0) then
    write(ilog,*) 'Warning. ',nzeros,' snapshots are without a neighbor (similar) structure. This &
 &may in some cases cause the clustering algorithm to misbehave.'
  end if
  write(ilog,*) '... done.'
  write(ilog,*) 
!
#ifdef LINK_NETCDF
  write(ilog,*) 'Dumping to binary NetCDF-file ...'
    call dump_nbl_nc()
  write(ilog,*) '... done.'
  write(ilog,*)
#endif
!
  cradius = bucr
!
  deallocate(bui2)
  deallocate(bui1)
  deallocate(buf2)
  deallocate(buf1)
  deallocate(maxrads)
  do i=1,nsets
    if (allocated(scluster(i)%snaps).EQV..true.) deallocate(scluster(i)%snaps)
    if (allocated(scluster(i)%tmpsnaps).EQV..true.) deallocate(scluster(i)%tmpsnaps)
    if (allocated(scluster(i)%sums).EQV..true.) deallocate(scluster(i)%sums)
    if (allocated(scluster(i)%map).EQV..true.) deallocate(scluster(i)%map)
    if (allocated(scluster(i)%children).EQV..true.) deallocate(scluster(i)%children)
    if (allocated(scluster(i)%wghtsnb).EQV..true.) deallocate(scluster(i)%wghtsnb)
    if (allocated(scluster(i)%lstnb).EQV..true.) deallocate(scluster(i)%lstnb)
    if (allocated(scluster(i)%flwnb).EQV..true.) deallocate(scluster(i)%flwnb)
  end do
  deallocate(scluster)
!
end
!
!------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
subroutine dump_nbl_nc()
!
  use iounit
  use system
  use mcsums
  use clusters
  use netcdf
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer i,ncid,ii,jj,freeunit,xlen,istart
  integer, ALLOCATABLE:: helper(:)
!  integer nf90_redef,nf90_put_att,nf90_def_dim
  logical exists
  character(MAXSTRLEN) attstring,dumpfile
  real(KIND=4), ALLOCATABLE:: prthlp(:)
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
!
#ifdef ENABLE_MPI
  call int2str(myrank,nod,re_aux(10))
  if (use_REMC.EQV..true.) then
    dumpfile = 'N_'//nod(1:re_aux(10))//'_FRAMES_NBL.nc'
  else
    dumpfile = 'FRAMES_NBL.nc'
  end if
#else
  dumpfile = 'FRAMES_NBL.nc'
#endif
  call strlims(dumpfile,ii,jj)
  inquire(file=dumpfile(ii:jj),exist=exists)
  if (exists.EQV..true.) then
    ncid = freeunit()
    open(unit=ncid,file=dumpfile(ii:jj),status='old')
    close(unit=ncid,status='delete')
  end if
  ncid = freeunit()
  call check_fileionetcdf( nf90_create(path=dumpfile(ii:jj), cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
!
! enable definition
  do i=1,MAXSTRLEN
    attstring(i:i) = " "
  end do
  attstring(1:7) = "CAMPARI"
  call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "program", attstring(1:7)) )
  attstring(1:7) = "       "
  attstring(1:3) = "XXX"
  call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "programVersion", attstring(1:3)) )
  attstring(1:3) = "   "
  call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "title", basename(1:bleng)) )
! define dimensions
  attstring(1:10) = "framepairs"
  call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:10), sum(cnblst(1:cstored)%nbs), cnc_ids(1)) )
  attstring(1:10) = "     " 
! define (not set) variables to hold distance and type of distance information
  attstring(1:9) = "snapshots"
  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_INT, cnc_ids(1), cnc_ids(2)) )
  attstring(1:9) = "neighbors"
  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_INT, cnc_ids(1), cnc_ids(3)) )
  attstring(1:9) = "distances"
  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_FLOAT, cnc_ids(1), cnc_ids(4)) )
  attstring(1:9) = "         "
  if (cdis_crit.eq.1) then
    xlen = 14
    attstring(1:xlen) = "torsional RMSD"
  else if (cdis_crit.eq.2) then
    xlen = 37
    attstring(1:xlen) = "inertial-mass-weighted torsional RMSD"
  else if (cdis_crit.eq.3) then
    xlen = 31
    attstring(1:xlen) = "RMSD of torsional Fourier terms"
  else if (cdis_crit.eq.4) then
    xlen = 54
    attstring(1:xlen) = "inertial-mass-weighted RMSD of torsional Fourier terms"
  else if (cdis_crit.eq.5) then
    if (align_for_clustering.EQV..true.) then
      xlen = 19
      attstring(1:xlen) = "aligned atomic RMSD"
    else
      xlen = 21
      attstring(1:xlen) = "unaligned atomic RMSD"
    end if
  else if (cdis_crit.eq.6) then
    if (align_for_clustering.EQV..true.) then
      xlen = 42
      attstring(1:xlen) = "aligned atomic RMSD (may be separate sets)"
    else
      xlen = 21
      attstring(1:xlen) = "unaligned atomic RMSD"
    end if
  else if (cdis_crit.eq.7) then
    xlen = 29
    attstring(1:xlen) = "internal distance vector RMSD"
  else if (cdis_crit.eq.8) then
    xlen = 43
    attstring(1:xlen) = "mass-weighted internal distance vector RMSD"
  else if (cdis_crit.eq.9) then
    xlen = 46
    attstring(1:xlen) = "locally weighted internal distance vector RMSD" 
  else if (cdis_crit.eq.10) then
    xlen = 43
    attstring(1:xlen) = "locally weighted, unaligned coordinate RMSD"  
  else
    write(ilog,*) 'Fatal. Unsupported distance criterion in dump_nbl_nc(...). This is an omission bug.'
    call fexit()
  end if
  call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), "type", attstring(1:xlen)) )
!  call check_fileionetcdf( nf90_def_var(ncid, attstring(1:xlen), NF90_FLOAT, cnc_ids(2), cnc_ids(5)) )
  attstring(1:5) = "units"
  if (cdis_crit.le.2) then
    attstring(6:12) = "degrees"
    call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:12)) )
  else if ((cdis_crit.ge.3).AND.(cdis_crit.le.4)) then
    attstring(6:13) = "unitless"
    call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:13)) )
  else
    attstring(6:13) = "angstrom"
    call check_fileionetcdf( nf90_put_att(ncid, cnc_ids(1), attstring(1:5), attstring(6:13)) )
  end if
  attstring(1:13) = "               "
! quit define mode
  call check_fileionetcdf( nf90_enddef(ncid) )
  call check_fileionetcdf( nf90_sync(ncid) )
!
! put the data
  allocate(helper(cstored))
  allocate(prthlp(cstored))
  istart = 1
  do i=1,cstored
    if (cnblst(i)%nbs.le.0) cycle
    helper(1:cnblst(i)%nbs) = i
    prthlp(1:cnblst(i)%nbs) = REAL(cnblst(i)%dis(1:cnblst(i)%nbs),KIND=4)
    call check_fileionetcdf( nf90_put_var(ncid, cnc_ids(2), helper(1:cnblst(i)%nbs), &
 &                                       start = (/ istart /), count = (/ cnblst(i)%nbs /)) )
    call check_fileionetcdf( nf90_put_var(ncid, cnc_ids(3), cnblst(i)%idx(1:cnblst(i)%nbs), &
 &                                       start = (/ istart /), count = (/ cnblst(i)%nbs /)) )
    call check_fileionetcdf( nf90_put_var(ncid, cnc_ids(4), REAL(prthlp(1:cnblst(i)%nbs),KIND=8), &
 &                                       start = (/ istart /), count = (/ cnblst(i)%nbs /)) )
    istart = istart + cnblst(i)%nbs
  end do
!
  deallocate(prthlp)
  deallocate(helper)
!
! close
  call check_fileionetcdf( nf90_close(ncid) )
! 
end
!
!---------------------------------------------
!
subroutine read_nbl_nc()
!
  use netcdf
  use clusters
  use iounit
!
  implicit none
!
  integer ncid,t1,t2,fndds,dimlen,nframes,i,ret,ilast,istart
  logical exists
  character(MAXSTRLEN) ucstr,trystr
  integer, ALLOCATABLE:: vnbs(:),vsnp(:)
  real(KIND=4), ALLOCATABLE:: vdis(:)
!
  call strlims(nblfilen,t1,t2)
  inquire (file=nblfilen(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from NetCDF (',nblfilen(t1:t2),') in setup_netcdftraj().'
    call fexit()
  end if
!
! open
 44 format('Warning. Ambiguous dimensions in NetCDF file. Encountered ',a,' twice but keeping&
 &only first.')
  call check_fileionetcdf( nf90_open(nblfilen(t1:t2), NF90_NOWRITE, ncid) )
!
! find the necessary dimensions: three are required
  fndds = 0
  nframes = 0
  do i=1,NF90_MAX_DIMS
    ret = nf90_inquire_dimension(ncid,i,trystr,dimlen)
    if (ret.eq.NF90_NOERR) then
      ucstr(1:15) = trystr(1:15)
      call toupper(ucstr(1:15))
      if (ucstr(1:10).eq.'FRAMEPAIRS') then
        if (fndds.eq.1) then 
          write(ilog,44) ucstr(1:10)
        else
          nframes = dimlen
          fndds = 1
          cnc_ids(1) = i
        end if
      end if
    else if (ret.eq.NF90_EBADDIM) then
!     do nothing
    else ! get us out of here
      call check_fileionetcdf( nf90_inquire_dimension(ncid,i,trystr,dimlen) )
    end if
  end do
!
  if (nframes.lt.1) then
    write(ilog,*) 'Fatal. NetCDF-file (',nblfilen(t1:t2),') has no neighbor &
 &data (empty containers).'
    call fexit()
  end if
!
! now find the necessary variables, only coordinates are required
  ret = nf90_inq_varid(ncid,"snapshots",cnc_ids(2))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "snapshots" not found in NetCDF-file (',&
 &nblfilen(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ncid,"snapshots",cnc_ids(2)) )
  end if
  ret = nf90_inq_varid(ncid,"neighbors",cnc_ids(3))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "neighbors" not found in NetCDF-file (',&
 &nblfilen(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ncid,"neighbors",cnc_ids(3)) )
  end if
  ret = nf90_inq_varid(ncid,"distances",cnc_ids(4))
  if (ret.eq.NF90_NOERR) then
!   do nothing
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "distances" not found in NetCDF-file (',&
 &nblfilen(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(ncid,"distances",cnc_ids(4)) )
  end if
!
  allocate(vnbs(nframes))
  allocate(vdis(nframes))
  allocate(vsnp(nframes))
!
! for some reason nf90_get_var chokes if the increment (count) is very large
  istart = 1
  ilast = min(nframes,10000)
  do while (istart.le.nframes)
    call check_fileionetcdf( nf90_get_var(ncid, cnc_ids(2), vsnp(istart:ilast), &
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    call check_fileionetcdf( nf90_get_var(ncid, cnc_ids(3), vnbs(istart:ilast),&
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    call check_fileionetcdf( nf90_get_var(ncid, cnc_ids(4), vdis(istart:ilast),&
 &                                        start = (/ istart /) , count = (/ ilast-istart+1 /)) )
    istart = ilast + 1
    ilast = min(nframes,istart+10000)
  end do
! 
  allocate(cnblst(cstored))
  cnblst(:)%nbs = 0
  i = 1
  ilast = 0
  do istart=vsnp(1),vsnp(nframes)
    if (istart.ne.vsnp(i)) cycle
    do while (vsnp(i).eq.istart)
      if (i.eq.nframes) exit
      if (vsnp(i+1).ne.vsnp(i)) exit
      i = i + 1
    end do
    cnblst(vsnp(i))%nbs = i - ilast
    cnblst(vsnp(i))%alsz = cnblst(vsnp(i))%nbs
    if (cnblst(vsnp(i))%nbs.gt.0) then
      allocate(cnblst(vsnp(i))%idx(cnblst(vsnp(i))%nbs))
      allocate(cnblst(vsnp(i))%dis(cnblst(vsnp(i))%nbs))
      allocate(cnblst(vsnp(i))%tagged(cnblst(vsnp(i))%nbs))
      cnblst(vsnp(i))%idx(:) = vnbs(ilast+1:i)
      cnblst(vsnp(i))%dis(:) = REAL(vdis(ilast+1:i),KIND=4)
      cnblst(vsnp(i))%tagged(:) = .false.
    else
      cnblst(vsnp(i))%nbs = 0
      cnblst(vsnp(i))%alsz = 2
      allocate(cnblst(vsnp(i))%idx(cnblst(vsnp(i))%alsz))
      allocate(cnblst(vsnp(i))%dis(cnblst(vsnp(i))%alsz))
      allocate(cnblst(vsnp(i))%tagged(cnblst(vsnp(i))%alsz))
    end if
    if (i.eq.nframes) exit
    ilast = i
    i = i  +1 
  end do
  deallocate(vsnp)
  deallocate(vnbs)
  deallocate(vdis)
!
end
!
#endif
!
!----------------------------------------------------------------------------
!
! this routine has a threads-specific variant below, birch_clustering() itself is not thread-safe
!
! set modei to 0 to obtain a silent version without post-processing
! set modei to -1 to obtain a silent version without post-processing and without populating scluster
! set modei to -2 to obtain a slightly less silent version of -1
!
subroutine birch_clustering(modei,nnodes)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer i,j,k,ii,jj,kk,ll,mm,atwo,kkf,thekk,nlst1,snapstart,snapend,snapinc,nnodes,modei
  integer, ALLOCATABLE:: kklst(:,:),kkhistory(:,:)
  integer(KIND=8) cnt1,cnt2
  RTYPE, ALLOCATABLE:: scrcts(:)
  RTYPE rdv,mind,helper,maxd(2),normer(4),qualmet(4)
  logical atrue,afalse
!
 33 format('ERROR B: ',i6,' could be part of ',i5,' at ',g14.7,' (last d',g14.6,')')
 34 format('NEXT HIGHER: ',i6,' could have been part of ',i5,' at ',g14.7,'.')
 35 format('ERROR C: ',i6,' is meant to a child of ',i5,' at level ',i5,' but has distance ',g14.7,'.')
!
  atrue = .true.
  afalse = .false.
  atwo = 2
  if (cmaxrad.le.cradius) cmaxrad = 2.0*cradius
!
  if (cleadermode.le.2) then
    snapstart = 1
    snapend = cstored
    snapinc = 1
  else
    snapstart = cstored
    snapend = 1
    snapinc = -1
  end if
!
  allocate(kkhistory(cstored,c_nhier+1))
  allocate(birchtree(c_nhier+1))
  do ii=1,c_nhier+1
    allocate(birchtree(ii)%cls(10))
    birchtree(ii)%ncls = 0
    birchtree(ii)%nclsalsz = 10
  end do
  allocate(scrcts(c_nhier+1))
  birchtree(1)%ncls = 1
  scrcts(c_nhier+1) = cradius
  if (c_nhier.gt.1) then
    scrcts(2) = cmaxrad
    do i=3,c_nhier
      scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - cradius) ! linear so far
    end do
  end if
  do i=1,c_nhier+1
    birchtree(i)%cls(:)%nmbrs = 0
    birchtree(i)%cls(:)%nb = 0
    birchtree(i)%cls(:)%nchildren = 0
    birchtree(i)%cls(:)%alsz = 0
    birchtree(i)%cls(:)%nbalsz = 0
    birchtree(i)%cls(:)%chalsz = 0
    birchtree(i)%cls(:)%parent = 0
  end do
!
  i = 1
  rdv = 1.0
  call cluster_addsnap(birchtree(1)%cls(1),i,rdv) ! init arrays etc.
  if (allocated(birchtree(1)%cls(1)%snaps).EQV..true.) deallocate(birchtree(1)%cls(1)%snaps)
  allocate(birchtree(1)%cls(1)%snaps(cstored))
  birchtree(1)%cls(1)%alsz = cstored
  birchtree(1)%cls(1)%nmbrs = cstored
  do i=1,cstored
    birchtree(1)%cls(1)%snaps(i) = i
  end do
  allocate(birchtree(1)%cls(1)%children(2))
  birchtree(1)%cls(1)%chalsz = 2
!
  if (modei.gt.0) then
    write(ilog,*)
    write(ilog,*) 'Now performing tree-based clustering ...'
  end if
!
  cnt1 = 0
  cnt2 = 0
!
! loop from coarsest to finest level
  do ii=2,c_nhier
    kk = 0 
    do mm=1,birchtree(ii-1)%ncls ! split
      if (ii.eq.2) then
        do k=snapstart,snapend,snapinc
          i = k
          mind = HUGE(mind)
          do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over clusters generated thus far by the snapshots in parent cluster
            ll = birchtree(ii-1)%cls(mm)%children(j)
            cnt1 = cnt1 + 1
            if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
            call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
            if (rdv.lt.mind) then
              mind = rdv
              jj = j
            end if
          end do
!         make a new cluster or not?
          if (mind.lt.scrcts(ii)) then
            kk = birchtree(ii-1)%cls(mm)%children(jj)
            call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
          else
            birchtree(ii)%ncls = birchtree(ii)%ncls + 1
            if (birchtree(ii)%ncls.gt.birchtree(ii)%nclsalsz) call scluster_resizelst(birchtree(ii)%nclsalsz,birchtree(ii)%cls)
            call cluster_addsnap(birchtree(ii)%cls(birchtree(ii)%ncls),i,rdv)
            call cluster_addchild(birchtree(ii-1)%cls(mm),mm,birchtree(ii)%cls(birchtree(ii)%ncls),birchtree(ii)%ncls)
            cnt2 = cnt2 + 1
          end if
        end do 
      else
        do k=1,birchtree(ii-1)%cls(mm)%nmbrs
          i = birchtree(ii-1)%cls(mm)%snaps(k)
          mind = HUGE(mind)
          do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over clusters generated thus far by the snapshots in parent cluster
            ll = birchtree(ii-1)%cls(mm)%children(j)
            cnt1 = cnt1 + 1
            if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
            call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
            if (rdv.lt.mind) then
              mind = rdv
              jj = j
            end if
          end do
!         make a new cluster or not?
          if (mind.lt.scrcts(ii)) then
            kk = birchtree(ii-1)%cls(mm)%children(jj)
            call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
          else
            birchtree(ii)%ncls = birchtree(ii)%ncls + 1
            if (birchtree(ii)%ncls.gt.birchtree(ii)%nclsalsz) call scluster_resizelst(birchtree(ii)%nclsalsz,birchtree(ii)%cls)
            call cluster_addsnap(birchtree(ii)%cls(birchtree(ii)%ncls),i,rdv)
            call cluster_addchild(birchtree(ii-1)%cls(mm),mm,birchtree(ii)%cls(birchtree(ii)%ncls),birchtree(ii)%ncls)
            cnt2 = cnt2 + 1
          end if
        end do
      end if
    end do
  end do
!  write(*,*) cnt1,' after clustering1'
!
! now refilter through static tree and store filter path
  kkhistory(:,:) = 0
  kkhistory(:,1) = 1
  thekk = 0
  do i=snapstart,snapend,snapinc
    do ii=2,c_nhier
      mind = HUGE(mind)
      mm = kkhistory(i,ii-1)
      if (mm.eq.0) cycle
      do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over children only
        ll = birchtree(ii-1)%cls(mm)%children(j)
        cnt1 = cnt1 + 1
        if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
        call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
        if (rdv.lt.mind) then
          mind = rdv
          jj = j
        end if
      end do
!     we will borrow the search space from the nearest available cluster, which affects the 2nd pass clusterings
!     only indirectly
      kkhistory(i,ii) = birchtree(ii-1)%cls(mm)%children(jj)
    end do
  end do
!  write(*,*) cnt1,' after clustering2'
!
! now recluster for one or more levels based on the filter paths (2-pass)
  do ii=c_nhier+1,min(c_nhier+1,max(3,c_nhier+1-c_multires)),-1
    kk = 0
    if (ii.le.c_nhier) then
      birchtree(ii)%ncls = 0
      birchtree(ii)%cls(:)%nmbrs = 0
      birchtree(ii)%cls(:)%parent = 0
      birchtree(ii-1)%cls(:)%nchildren = 0
    end if
!
    do i=snapstart,snapend,snapinc
      mm = kkhistory(i,ii-1)
      mind = HUGE(mind)
!     split by mm
      do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over clusters generated thus far by the snapshots in parent cluster
        ll = birchtree(ii-1)%cls(mm)%children(j)
        cnt1 = cnt1 + 1
        if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
        call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
        if (rdv.lt.mind) then
          mind = rdv
          jj = j
        end if
      end do
!     make a new cluster or not?
      if (mind.lt.scrcts(ii)) then
        kk = birchtree(ii-1)%cls(mm)%children(jj)
        call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
      else
        birchtree(ii)%ncls = birchtree(ii)%ncls + 1
        if (birchtree(ii)%ncls.gt.birchtree(ii)%nclsalsz) call scluster_resizelst(birchtree(ii)%nclsalsz,birchtree(ii)%cls)
        call cluster_addsnap(birchtree(ii)%cls(birchtree(ii)%ncls),i,rdv)
        call cluster_addchild(birchtree(ii-1)%cls(mm),mm,birchtree(ii)%cls(birchtree(ii)%ncls),birchtree(ii)%ncls)
        cnt2 = cnt2 + 1
        kk = birchtree(ii)%ncls
      end if
    end do
  end do
!  write(*,*) cnt1,' after clustering3'
!
  deallocate(kkhistory)
!
 66 format('Level    # Clusters     Threshold     Total Snaps    Total Children')
 67 format(i9,i10,1x,g14.4,4x,i12,4x,i12)
 68 format(i9,i10,5x,a7,7x,i12,4x,i12)
  if ((modei.gt.0).OR.(modei.eq.-2)) then
    write(ilog,66)
    write(ilog,68) c_nhier+1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
 &               sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
    do i=2,c_nhier+1
      write(ilog,67) c_nhier+2-i,birchtree(i)%ncls,scrcts(i),sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
 &               sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren)
    end do
    write(ilog,*) '---------------------------------------------------------------------'
    write(ilog,*)
    write(ilog,*) '... done after a total of ',cnt1,' distance evaluations.'
    write(ilog,*)
  end if
!
  if (refine_clustering.EQV..true.) then
!
    allocate(kklst(cstored,2))
!
    do j=1,birchtree(c_nhier+1)%ncls
      call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
    end do
    call quality_of_clustering(birchtree(c_nhier+1)%ncls,birchtree(c_nhier+1)%cls,scrcts(c_nhier+1),qualmet)
    if (modei.gt.0) then
      write(ilog,*) 'Now merging clusters that yield joint reduced average intracluster distance ...'
    end if
 77 format('Would join ',i5,' (',i6,') and ',i5,'(',i6,') from: ',/,'Diam: ',g10.4,' / ',g10.4,'; Rad.: ',&
 &g10.4,' / ',g10.4,' to ',g10.4,' / ',g10.4,'.')
    cnt1 = 0
    cnt2 = 0
    do i=1,birchtree(c_nhier)%ncls
      nlst1 = 0
      do kkf=1,birchtree(c_nhier)%ncls
        if (i.eq.kkf) cycle
        call cluster_to_cluster_d(rdv,birchtree(c_nhier)%cls(i),birchtree(c_nhier)%cls(kkf))
        cnt1 = cnt1 + 1
        if (rdv.lt.scrcts(c_nhier)) then
          nlst1 = nlst1 + 1
          kklst(nlst1,1) = kkf
        end if
      end do
      do j=1,birchtree(c_nhier)%cls(i)%nchildren
        jj = birchtree(c_nhier)%cls(i)%children(j)
        if (birchtree(c_nhier+1)%cls(jj)%nmbrs.le.0) cycle
        do kkf=1,nlst1
          do k=1,birchtree(c_nhier)%cls(kklst(kkf,1))%nchildren
            kk = birchtree(c_nhier)%cls(kklst(kkf,1))%children(k)
            if (birchtree(c_nhier+1)%cls(kk)%nmbrs.le.0) cycle
            if (jj.eq.kk) cycle
            call clusters_joint_diam(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk),rdv,helper)
            cnt1 = cnt1 + 1
            normer(1) = 0.5*birchtree(c_nhier+1)%cls(jj)%nmbrs*(birchtree(c_nhier+1)%cls(jj)%nmbrs-1.0)
            normer(2) = 0.5*birchtree(c_nhier+1)%cls(kk)%nmbrs*(birchtree(c_nhier+1)%cls(kk)%nmbrs-1.0)
            normer(3) = 1.0/(1.0*(birchtree(c_nhier+1)%cls(jj)%nmbrs + birchtree(c_nhier+1)%cls(kk)%nmbrs))
            normer(4) = 0.0
            if (sum(normer(1:2)).gt.0) normer(4) = 1.0/sum(normer(1:2))
            maxd(1) =  normer(4)*(normer(1)*birchtree(c_nhier+1)%cls(jj)%diam + &
              &                   normer(2)*birchtree(c_nhier+1)%cls(kk)%diam)
            maxd(2) =  normer(3)*(birchtree(c_nhier+1)%cls(jj)%nmbrs*birchtree(c_nhier+1)%cls(jj)%radius + &
   &                              birchtree(c_nhier+1)%cls(kk)%nmbrs*birchtree(c_nhier+1)%cls(kk)%radius)
            if ((helper.le.maxd(2)).OR.(rdv.le.maxd(1))) then
              if (birchtree(c_nhier+1)%cls(jj)%nmbrs.gt.birchtree(c_nhier+1)%cls(kk)%nmbrs) then
                call join_clusters(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk))
                birchtree(c_nhier+1)%cls(jj)%diam = rdv
                birchtree(c_nhier+1)%cls(jj)%radius = helper
                cnt2 = cnt2 + 1
              else
                call join_clusters(birchtree(c_nhier+1)%cls(kk),birchtree(c_nhier+1)%cls(jj))
                birchtree(c_nhier+1)%cls(kk)%diam = rdv
                birchtree(c_nhier+1)%cls(kk)%radius = helper
                cnt2 = cnt2 + 1
                exit
              end if
            end if
          end do
          if (birchtree(c_nhier+1)%cls(jj)%nmbrs.eq.0) exit
        end do
      end do
    end do
    if (modei.gt.0) then
      write(ilog,*) '... done after a total of ',cnt2,' merges requiring ',cnt1,' additional &
 &distance or joint size evaluations.'
      write(ilog,*)
    end if
    deallocate(kklst)
  end if
!
! now shorten list and resort
  do k=c_nhier+1,2,-1
    call clusters_shorten(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls)
    do j=1,birchtree(k)%ncls
      call cluster_calc_params(birchtree(k)%cls(j),scrcts(k))
    end do
    do i=1,birchtree(k)%ncls
      call cluster_getcenter(birchtree(k)%cls(i))
    end do
    call clusters_sort(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls,afalse)
  end do
!
  if (modei.ge.0) then
!   lastly, copy into global cluster array
    allocate(scluster(birchtree(c_nhier+1)%ncls))
    do i=1,birchtree(c_nhier+1)%ncls
      call copy_cluster(birchtree(c_nhier+1)%cls(i),scluster(i))
    end do
  end if
!
  nnodes = birchtree(c_nhier+1)%ncls
!
 63 format(i7,1x,i7,1x,i8,1000(1x,g12.5))
 64 format(1000(g12.5,1x))
 69 format('------------- CLUSTER SUMMARY (THRESHOLD OF ',g12.5,')')
!
  if (modei.gt.0) then
    do k=c_nhier+1,max(2,c_nhier+1-c_multires),-1
      write(ilog,69) scrcts(k)
      write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
!    do i=1,birchtree(c_nhier+1)%ncls
!      write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
!    end do
      do i=1,birchtree(k)%ncls
        write(ilog,63) i,birchtree(k)%cls(i)%nmbrs,birchtree(k)%cls(i)%center,birchtree(k)%cls(i)%diam,birchtree(k)%cls(i)%radius
      end do
      write(ilog,*) '----------------------------------------------------'
      write(ilog,*)
      call quality_of_clustering(birchtree(k)%ncls,birchtree(k)%cls(1:birchtree(k)%ncls),scrcts(k),qualmet)
    end do
!
    atrue = .true.
    kk = 0
    call gen_graph_from_clusters(scluster,nnodes,atrue,kk)
    call graphml_helper_for_clustering(scluster,nnodes)
    call vmd_helper_for_clustering(scluster,nnodes)
  end if
!
  deallocate(scrcts)
!
end
!
!----------------------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
! set modei to 0 to obtain a silent version without post-processing
! set modei to -1 to obtain a silent version without post-processing and without populating scluster
! set modei to -2 to obtain a slightly less silent version of -1
!
! have a look at birch_clustering(...) to understand the decomposition of the algorithm first
!
subroutine birch_clustering_threads(modei,nnodes,tpi2)
!
  use clusters
  use iounit
  use interfaces
  use threads, ONLY: thr_hlper,thr_rhlper,thrdat,thr_limits,ithread
!
  implicit none
!
  integer, INTENT(IN):: tpi2
!
  integer i,j,k,ii,jj,kk,ll,mm,atwo,ithresh,kkf,nlst1,snapstart,snapend,snapinc,nnodes,modei,tpi
  integer(KIND=8) cnt1,cnt2,tmvars(8)
  RTYPE, ALLOCATABLE:: scrcts(:)
  RTYPE rdv,mind,helper,maxd(2),normer(4),qualmet(4),rdv2,rdv3,rdv4
  logical atrue,afalse
  integer tpx,tpn,OMP_GET_NUM_THREADS,tbnds(6),bla,bla2,bla3,bla4,adds1,adds2
!
 33 format('ERROR B: ',i6,' could be part of ',i5,' at ',g14.7,' (last d',g14.6,')')
 34 format('NEXT HIGHER: ',i6,' could have been part of ',i5,' at ',g14.7,'.')
 35 format('ERROR C: ',i6,' is meant to a child of ',i5,' at level ',i5,' but has distance ',g14.7,'.')
!
  if (tpi2.le.0) then
    write(ilog,*) 'Fatal. Birch_clustering_threads(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
!
  tpn = OMP_GET_NUM_THREADS()
  tpi = tpi2
  snapstart = 1
  snapend = cstored
  call threads_bounds(snapstart,snapend,tpi,tpn,tbnds(1:2))
!$OMP SINGLE
  allocate(thr_hlper(cstored,tpn))
  allocate(thr_rhlper(2,tpn))
  allocate(csnap2clus(cstored,c_nhier+1))
!$OMP END SINGLE
!
  if (tpi.le.1) call System_Clock(count=tmvars(1))
!
  atrue = .true.
  afalse = .false.
  atwo = 2
  if (cmaxrad.le.cradius) cmaxrad = 2.0*cradius
!
! initialize the global tree structure first (could be simplified)
!
!$OMP SINGLE
  allocate(birchtree(c_nhier+1))
  do ii=1,c_nhier+1
    allocate(birchtree(ii)%cls(10))
    birchtree(ii)%ncls = 0
    birchtree(ii)%nclsalsz = 10
  end do
  birchtree(1)%ncls = 1
!
  do i=1,c_nhier+1
    birchtree(i)%cls(:)%nmbrs = 0
    birchtree(i)%cls(:)%nb = 0
    birchtree(i)%cls(:)%nchildren = 0
    birchtree(i)%cls(:)%alsz = 0
    birchtree(i)%cls(:)%nbalsz = 0
    birchtree(i)%cls(:)%chalsz = 0
    birchtree(i)%cls(:)%parent = 0
  end do
!
  allocate(birchtree(1)%cls(1)%snaps(cstored))
  birchtree(1)%cls(1)%alsz = cstored
  birchtree(1)%cls(1)%nmbrs = cstored
  allocate(birchtree(1)%cls(1)%children(2))
  birchtree(1)%cls(1)%chalsz = 2
  cdevalcnt = 0
  if (allocated(thr_btree).EQV..true.) deallocate(thr_btree)
  allocate(thr_btree(c_nhier+1,tpn))
!$OMP END SINGLE
  if (cleadermode.le.2) then
    do i=tbnds(1),tbnds(2)
      birchtree(1)%cls(1)%snaps(i) = i
    end do 
  else
    do i=tbnds(1),tbnds(2)
      birchtree(1)%cls(1)%snaps(i) = cstored-i+1
    end do 
  end if
!
! setup the global but thread-locally-used tree object next
  allocate(thr_btree(1,tpi)%cls(1))
  thr_btree(1,tpi)%ncls = 1
  thr_btree(1,tpi)%nclsalsz = 1
  do i=2,c_nhier+1
    allocate(thr_btree(i,tpi)%cls(10))
    thr_btree(i,tpi)%ncls = 0
    thr_btree(i,tpi)%nclsalsz = 10
    thr_btree(i,tpi)%cls(:)%nmbrs = 0
    thr_btree(i,tpi)%cls(:)%nb = 0
    thr_btree(i,tpi)%cls(:)%nchildren = 0
    thr_btree(i,tpi)%cls(:)%alsz = 0
    thr_btree(i,tpi)%cls(:)%nbalsz = 0
    thr_btree(i,tpi)%cls(:)%chalsz = 0
    thr_btree(i,tpi)%cls(:)%parent = 0
    thr_btree(i,tpi)%cls(:)%geni = 0
  end do
!
! the schedule of thresholds
  allocate(scrcts(c_nhier+1))
  scrcts(c_nhier+1) = cradius
  if (c_nhier.gt.1) then
    scrcts(2) = cmaxrad
    do i=3,c_nhier
      scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - cradius) ! linear so far
    end do
  end if
!
  if ((modei.gt.0).AND.(tpi.le.1)) then
    write(ilog,*)
    write(ilog,*) 'Now performing tree-based clustering ...'
  end if
!
! counters for cost and balance reporting
  cnt1 = 0
  cnt2 = 0
!
! the algorithm is split into 3 parts: 1) initial tree creation
!                                      2) static tree path generation
!                                      3) reclustering from leafs inward
!
! step 1) is the most difficult to parallelize and involves a loop from coarsest to finest level with an embedded loop
! across snapshots
! parallelization eventually becomes simple by splitting the load according to the parent clusters at the next level closer
! to the root (no cross-interference possible within algorithm)
! the problem is that for levels near the root (where there are 1 or just few clusters), this is not possible
! we therefore use an alternative strategy requiring a subsequent and heuristic merge
! note that it is infeasible to try to let multiple threads operate on the birchtree structure simultaneously
! (the synchronization cost, even with possible delayed updates, will always be large compared to the extremely cheap
! cost per snapshot); this means that all synchronization operations are located outside of the snapshot loops
  ithresh = cstored/(1.0*cbirchbatch)
!$OMP BARRIER
  do ii=2,c_nhier
    kkf = 0
    adds1 = 0
    adds2 = 0
    helper = scrcts(ii)
!
!   reshape thread-specific tree object for prior level
    if (allocated(thr_btree(ii-1,tpi)%cls).EQV..true.) deallocate(thr_btree(ii-1,tpi)%cls)
    allocate(thr_btree(ii-1,tpi)%cls(birchtree(ii-1)%ncls))
    thr_btree(ii-1,tpi)%ncls = birchtree(ii-1)%ncls
    thr_btree(ii-1,tpi)%nclsalsz = birchtree(ii-1)%ncls
    thr_btree(ii-1,tpi)%cls(:)%nmbrs = 0
    thr_btree(ii-1,tpi)%cls(:)%nchildren = 0
    thr_btree(ii-1,tpi)%cls(:)%alsz = 0
    thr_btree(ii-1,tpi)%cls(:)%chalsz = 0
    thr_btree(ii-1,tpi)%cls(:)%parent = 0
    thr_btree(ii-1,tpi)%cls(:)%geni = 0
!   first pass focuses on those cases requiring merging
    do mm=1,birchtree(ii-1)%ncls
      if (birchtree(ii-1)%cls(mm)%nmbrs.le.ithresh) cycle ! only go over the internally parallelized cases
      kkf = kkf + 1
       k = 1
       i = birchtree(ii-1)%cls(mm)%nmbrs
      call threads_bounds(k,i,tpi,tpn,tbnds(3:4))
      do k=tbnds(3),tbnds(4) ! tpi,birchtree(ii-1)%cls(mm)%nmbrs,tpn
        i = birchtree(ii-1)%cls(mm)%snaps(k)
        mind = HUGE(mind)
        do j=1,thr_btree(ii-1,tpi)%cls(mm)%nchildren  ! these children clusters are only the thread-local ones
          ll = thr_btree(ii-1,tpi)%cls(mm)%children(j)
          cnt1 = cnt1 + 1
          if ((thr_btree(ii,tpi)%cls(ll)%center.le.0).OR.(thr_btree(ii,tpi)%cls(ll)%center.gt.cstored)) call fexit()
          call snap_to_cluster_d(rdv,thr_btree(ii,tpi)%cls(ll),i)
          if (rdv.lt.mind) then
            mind = rdv
            jj = j
          end if
        end do
!       make a new thread-local cluster or not?
        if (mind.lt.scrcts(ii)) then
          kk = thr_btree(ii-1,tpi)%cls(mm)%children(jj)
          call cluster_addsnap(thr_btree(ii,tpi)%cls(kk),i,rdv)
        else
          thr_btree(ii,tpi)%ncls = thr_btree(ii,tpi)%ncls + 1
          if (thr_btree(ii,tpi)%ncls.gt.thr_btree(ii,tpi)%nclsalsz) then
            call scluster_resizelst(thr_btree(ii,tpi)%nclsalsz,thr_btree(ii,tpi)%cls)
          end if
          thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls)%geni = mm
          call cluster_addsnap(thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls),i,rdv)
          call cluster_addchild(thr_btree(ii-1,tpi)%cls(mm),mm,thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls),thr_btree(ii,tpi)%ncls)
          cnt2 = cnt2 + 1
        end if
      end do
    end do
    adds1 = thr_btree(ii,tpi)%ncls
!   now deal with the remaining entries; note that the schedule requires constant synchronization of sorts
    if (kkf.lt.birchtree(ii-1)%ncls) then
!$OMP DO SCHEDULE(DYNAMIC,1)
      do mm=1,birchtree(ii-1)%ncls
        if (birchtree(ii-1)%cls(mm)%nmbrs.gt.ithresh) cycle ! only go over the trivially parallelized cases
        do k=1,birchtree(ii-1)%cls(mm)%nmbrs
          i = birchtree(ii-1)%cls(mm)%snaps(k)
          mind = HUGE(mind)
          do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over clusters generated thus far by the snapshots in parent cluster
            ll = birchtree(ii-1)%cls(mm)%children(j)
            cnt1 = cnt1 + 1
            if ((thr_btree(ii,tpi)%cls(ll)%center.le.0).OR.(thr_btree(ii,tpi)%cls(ll)%center.gt.cstored)) call fexit()
            call snap_to_cluster_d(rdv,thr_btree(ii,tpi)%cls(ll),i)
            if (rdv.lt.mind) then
              mind = rdv
              jj = j
            end if
          end do
!         make a new cluster or not?
          if (mind.lt.scrcts(ii)) then
            kk = birchtree(ii-1)%cls(mm)%children(jj)
            call cluster_addsnap(thr_btree(ii,tpi)%cls(kk),i,rdv)
          else
            thr_btree(ii,tpi)%ncls = thr_btree(ii,tpi)%ncls + 1
            if (thr_btree(ii,tpi)%ncls.gt.thr_btree(ii,tpi)%nclsalsz) then
              call scluster_resizelst(thr_btree(ii,tpi)%nclsalsz,thr_btree(ii,tpi)%cls)
            end if
            thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls)%geni = mm
            call cluster_addsnap(thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls),i,rdv)
!           here, we can modify the mother cluster directly
            call cluster_addchild(birchtree(ii-1)%cls(mm),mm,thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls),thr_btree(ii,tpi)%ncls)
            cnt2 = cnt2 + 1
          end if
        end do
!        if (ii.le.5) write(*,*) 17-ii,tpi,birchtree(ii-1)%cls(mm)%nmbrs,birchtree(ii-1)%cls(mm)%nchildren
      end do
!$OMP END DO
    end if
!    if ((kkf.le.0).AND.(tpi.le.1)) write(*,*) 'no mor merges at ',ii
!   the rest of stage 1 is concerned with assembling an appropriate birchtree(ii) structure by copying/merging
    adds2 = thr_btree(ii,tpi)%ncls
    thr_btree(ii,tpi)%ncls = adds1 ! fake it
!$OMP BARRIER
!    if (tpi.le.1) write(*,*) ii,'merging ',sum(thr_btree(ii,1:tpn)%ncls),' down ...'
    bla = 0
    bla2 = 0
    if (kkf.gt.0) then ! this large block is for creating the results in birchtree of the large ii-1 clusters scanned above
!   learn a merging threshold from the data obtained thus far (the problem is that the most appropriate value
!   depends on tree height, tree level, and data dimensionality); default is 1.0*scrcts(ii) as above
      thr_rhlper(1:2,tpi) = 0.0  
      do i=1,thr_btree(ii,tpi)%ncls
        if (thr_btree(ii,tpi)%cls(i)%nmbrs.gt.50) then
          call cluster_calc_params(thr_btree(ii,tpi)%cls(i),rdv)
          thr_rhlper(1,tpi) = thr_rhlper(1,tpi) + thr_btree(ii,tpi)%cls(i)%nmbrs*thr_btree(ii,tpi)%cls(i)%radius
          thr_rhlper(2,tpi) = thr_rhlper(2,tpi) + thr_btree(ii,tpi)%cls(i)%nmbrs
        end if
      end do
!$OMP BARRIER
      rdv = sum(thr_rhlper(2,1:tpn))
      if (rdv.gt.0.0) then
        helper = min(scrcts(ii),sum(thr_rhlper(1,1:tpn))/rdv)
      end if
!      if (tpi.le.1) write(*,*) tpi,helper,scrcts(ii),helper/scrcts(ii)
!     to merge smaller into larger clusters, create some helper data structures
      if (tpi.eq.1) then
        bla = 1
      else
        bla = sum(thr_btree(ii,1:(tpi-1))%ncls) + 1
      end if
      bla2 = sum(thr_btree(ii,1:tpn)%ncls)
      do mm=1,thr_btree(ii,tpi)%ncls
        thr_btree(ii,tpi)%cls(mm)%parent = bla + mm - 1
        thr_hlper(bla+mm-1,1) = thr_btree(ii,tpi)%cls(mm)%nmbrs
      end do
!$OMP BARRIER
!     now determine a merging map, hi-jack parent for this purpose
      mm = 0
      do tpx=1,tpn
        if (tpx.eq.tpi) then
          mm = mm + thr_btree(ii,tpi)%ncls
          cycle
        end if
        do jj=1,thr_btree(ii,tpx)%ncls
          mm = mm + 1
          do kk=1,thr_btree(ii,tpi)%ncls
            if (thr_btree(ii,tpi)%cls(kk)%geni.ne.thr_btree(ii,tpx)%cls(jj)%geni) cycle ! do not merge across parents
!           calculate the radius (and diameter) of the joint cluster efficiently
!           also calculate the centroid-to-centroid and the mean snap-to-snap distances (all possibly approximate)
            call clusters_mean_diff(thr_btree(ii,tpi)%cls(kk),thr_btree(ii,tpx)%cls(jj),rdv4)
            call cluster_to_cluster_d(rdv3,thr_btree(ii,tpi)%cls(kk),thr_btree(ii,tpx)%cls(jj))
            if (thr_btree(ii,tpi)%cls(kk)%nmbrs.gt.thr_btree(ii,tpx)%cls(jj)%nmbrs) then
              call clusters_joint_radius(thr_btree(ii,tpi)%cls(kk),thr_btree(ii,tpx)%cls(jj),rdv,rdv2)
            else
              call clusters_joint_radius(thr_btree(ii,tpx)%cls(jj),thr_btree(ii,tpi)%cls(kk),rdv,rdv2)
            end if
 3 format(i6,i6,10(g11.5,1x))
            if (thr_hlper(thr_btree(ii,tpi)%cls(kk)%parent,1).le.thr_hlper(mm,1)) then
              if ((thr_hlper(thr_btree(ii,tpi)%cls(kk)%parent,1).eq.thr_hlper(mm,1)).AND.&
 &                (thr_btree(ii,tpi)%cls(kk)%parent.lt.mm)) cycle
              if ((rdv3.le.cmergecrit*helper).AND.((rdv4-rdv2)/rdv2.le.(cmergecrit-1.0))) then
!            write(*,3) thr_btree(ii,tpi)%cls(kk)%nmbrs,thr_btree(ii,tpx)%cls(jj)%nmbrs,scrcts(ii),helper,rdv,rdv2,rdv3,rdv4
                thr_btree(ii,tpi)%cls(kk)%parent = mm
              end if
            end if
          end do
        end do
      end do
!$OMP SINGLE
      mm = 0
      do tpx=1,tpn
        do jj=1,thr_btree(ii,tpx)%ncls
          mm = mm + 1
          csnap2clus(mm,1) = tpx
          csnap2clus(mm,2) = jj
        end do
      end do
!$OMP END SINGLE
      bla2 = 0
      mm = 0
      do tpx=1,tpn
        do jj=1,thr_btree(ii,tpx)%ncls
          mm = mm + 1
          if (thr_btree(ii,tpx)%cls(jj)%parent.eq.mm) then
            bla2 = bla2 + 1
            if (tpi.eq.tpx) thr_hlper(mm,1) = bla2 ! map for indexing birchtree(ii)%cls(:)
          else if (tpi.eq.tpx) then
!           possibly find our true parent
            kk = thr_btree(ii,tpx)%cls(jj)%parent
            do while (kk.ne.thr_btree(ii,csnap2clus(kk,1))%cls(csnap2clus(kk,2))%parent)
              kk = thr_btree(ii,csnap2clus(kk,1))%cls(csnap2clus(kk,2))%parent
            end do
            thr_btree(ii,tpx)%cls(jj)%parent = kk
            thr_hlper(mm,1) = 0
          end if
        end do
      end do
!$OMP BARRIER
    end if
!    if (tpi.le.1) write(*,*) 'to ...',bla2
    thr_btree(ii,tpi)%ncls = adds2 - adds1
!$OMP BARRIER
    if (tpi.le.1) then
      bla3 = bla2
    else
      bla3 = sum(thr_btree(ii,1:(tpi-1))%ncls) + bla2
    end if
    bla4 = sum(thr_btree(ii,1:tpn)%ncls)
!$OMP SINGLE
!   recreate birchtree data structures as needed
    if (allocated(birchtree(ii)%cls).EQV..true.) deallocate(birchtree(ii)%cls)
    allocate(birchtree(ii)%cls(bla2+bla4))
    birchtree(ii)%ncls = bla2+bla4
    birchtree(ii)%nclsalsz = bla2+bla4
    do mm=1,birchtree(ii-1)%ncls
!      if (allocated(birchtree(ii-1)%cls(mm)%children).EQV..true.) deallocate(birchtree(ii-1)%cls(mm)%children)
      birchtree(ii-1)%cls(mm)%nchildren = 0
    end do
!$OMP END SINGLE
    thr_btree(ii,tpi)%ncls = adds1
    if (bla2.gt.0) then ! globally the same
!$OMP BARRIER
      do mm=bla,bla+thr_btree(ii,tpi)%ncls-1
        if (thr_hlper(mm,1).gt.0) then ! let the true owner of the mapped cluster create the superset
          call copy_cluster(thr_btree(ii,tpi)%cls(mm-bla+1),birchtree(ii)%cls(thr_hlper(mm,1)))
!          birchtree(ii)%cls(thr_hlper(mm,1))%parent = 1
          do tpx=1,tpn
            do jj=1,thr_btree(ii,tpx)%ncls
              if ((tpi.eq.tpx).AND.(jj.eq.(mm-bla+1))) cycle
              if (thr_btree(ii,tpx)%cls(jj)%parent.eq.mm) then
                call join_clusters(birchtree(ii)%cls(thr_hlper(mm,1)),thr_btree(ii,tpx)%cls(jj))
              end if
            end do
          end do
        end if
      end do
!$OMP BARRIER
    end if
    if (bla4.gt.0) then ! globally the same
!     simple copy back
      do i=adds1+1,adds2
        call copy_cluster(thr_btree(ii,tpi)%cls(i),birchtree(ii)%cls(bla3+i-adds1))
!       birchtree(ii-1)%cls(xx) was exclusively operated on by this thread and has the right number of children
!       thus, simply overwrite index
        kk = thr_btree(ii,tpi)%cls(i)%geni
        birchtree(ii-1)%cls(kk)%nchildren = birchtree(ii-1)%cls(kk)%nchildren + 1
        birchtree(ii-1)%cls(kk)%children(birchtree(ii-1)%cls(kk)%nchildren) = bla3+i-adds1
      end do
!$OMP BARRIER
    end if
    if (bla2.gt.0) then ! globally the same
!$OMP SINGLE
!     fix parent-child relations for the merged clusters
      mm = 0
      do tpx=1,tpn
        do jj=1,thr_btree(ii,tpx)%ncls
          mm = mm + 1
          if (thr_hlper(mm,1).gt.0) then
            kk = thr_btree(ii,tpx)%cls(jj)%geni
            call cluster_addchild(birchtree(ii-1)%cls(kk),kk,birchtree(ii)%cls(thr_hlper(mm,1)),thr_hlper(mm,1))
          end if
        end do
      end do
!$OMP END SINGLE NOWAIT
    end if
!$OMP BARRIER
  end do
  if (thrdat%verbosity.gt.0) then
!$OMP CRITICAL(DEBUG_PRINT)
    if (tpi.le.1) then
      call System_Clock(count=tmvars(2))
      if (thrdat%verbosity.gt.0) write(ithread,*) (tmvars(2)-tmvars(1))/thrdat%rate,' s for Clustering #1'
      tmvars(1) = tmvars(2)
    end if
    if (thrdat%verbosity.gt.1) write(ithread,*) tpi,' so far did ',cnt1
!$OMP END CRITICAL(DEBUG_PRINT)
  end if
!
! step 2 is the simplest as it uses a static shared data structures and has no competing writes
! it consists of looping over all snapshots and the H-1 levels again to store a path through the static tree
! this refiltering is needed for appropriate cluster consistency
  do ii=2,c_nhier
    if (allocated(thr_btree(ii,tpi)%cls).EQV..true.) deallocate(thr_btree(ii,tpi)%cls)
    allocate(thr_btree(ii,tpi)%cls(birchtree(ii)%ncls))
    thr_btree(ii,tpi)%ncls = birchtree(ii)%ncls
    thr_btree(ii,tpi)%nclsalsz = birchtree(ii)%ncls
    thr_btree(ii,tpi)%cls(:)%alsz = 0
    thr_btree(ii,tpi)%cls(:)%nmbrs = 0
    thr_btree(ii,tpi)%cls(:)%parent = 0
  end do
  do ii=2,c_nhier
    csnap2clus(tbnds(1):tbnds(2),ii) = 0
  end do
  csnap2clus(tbnds(1):tbnds(2),1) = 1
!
  if (cleadermode.le.2) then
    snapstart = tbnds(1)
    snapend = tbnds(2)
    snapinc = 1 ! tpn
  else
    snapstart = -cstored
    snapend = -1
    call threads_bounds(snapstart,snapend,tpi,tpn,tbnds(5:6))
    snapstart = -tbnds(5)
    snapend = -tbnds(6)
    snapinc = -1
  end if
!
!$OMP BARRIER
  do ii=2,c_nhier
    do i=snapstart,snapend,snapinc ! hard to know whether this is going to be balanced
!    do ii=2,c_nhier
      mind = HUGE(mind)
      mm = csnap2clus(i,ii-1)
      if (mm.eq.0) call fexit()
      jj = 0
      do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over children only
        ll = birchtree(ii-1)%cls(mm)%children(j)
        cnt1 = cnt1 + 1
        if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
        call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
        if (rdv.lt.mind) then
          mind = rdv
          jj = j
        end if
      end do
!     we will borrow the search space from the nearest available cluster, which affects the 2nd pass clusterings
!     only indirectly
      kk = birchtree(ii-1)%cls(mm)%children(jj)
      call cluster_addjustsnap(thr_btree(ii,tpi)%cls(kk),i) ! thread-safe and simplified
      csnap2clus(i,ii) = birchtree(ii-1)%cls(mm)%children(jj) ! thread-safe
    end do
  end do
!
!$OMP BARRIER
  do ii=2,c_nhier
    csnap2clus(tbnds(1):tbnds(2),ii) = 0
  end do
  if (tpi.eq.1) then
    csnap2clus(1,1) = 1
    thr_hlper(1,1) = 1
  end if
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC,1)
  do ii=2,c_nhier
    birchtree(ii-1)%cls(1:birchtree(ii-1)%ncls)%nchildren = 0
    thr_hlper(ii,1) = 0
    do j=1,birchtree(ii)%ncls
      do tpx=1,tpn
        if (thr_btree(ii,tpx)%cls(j)%nmbrs.gt.0) then
          thr_hlper(ii,1) =  thr_hlper(ii,1) + 1
          csnap2clus(thr_hlper(ii,1),ii) = j  
          exit
        end if
      end do
    end do
  end do
!$OMP END DO
!
  if (thrdat%verbosity.gt.0) then
!$OMP CRITICAL(DEBUG_PRINT)
    if (thrdat%verbosity.gt.1) write(ithread,*) tpi,' now did ',cnt1
    if (tpi.le.1) then
      call System_Clock(count=tmvars(2))
      if (thrdat%verbosity.gt.0) write(ithread,*) (tmvars(2)-tmvars(1))/thrdat%rate,' s for Clustering #2'
      tmvars(1) = tmvars(2)
    end if
!$OMP END CRITICAL(DEBUG_PRINT)
  end if
!
! the third and last part involves creating the actual final clusterings for 1 or more levels.
! the cluster merges across threads implemented for step 1 obviously affect the results, and
! here we should therefore avoid all such merges; this unfortunately means that load balance may
! not be achievable; to help this somewhat, we parallelize across both levels and parent clusters
! simultaneously, which costs some in memory footprint (since thr_btree objects are appended)
!
! with a minimal tree, we are missing everything from above -> populate single root cluster in thr_btree
  if (c_nhier.eq.1) then
    thr_btree(1,tpi)%cls(1)%nmbrs = 0
    thr_btree(1,tpi)%cls(1)%alsz = max(tbnds(2)-tbnds(1)+1,1)
    allocate(thr_btree(1,tpi)%cls(1)%snaps(thr_btree(1,tpi)%cls(1)%alsz))
    do i=tbnds(1),tbnds(2)
      thr_btree(1,tpi)%cls(1)%nmbrs = thr_btree(1,tpi)%cls(1)%nmbrs + 1 
      thr_btree(1,tpi)%cls(1)%snaps(thr_btree(1,tpi)%cls(1)%nmbrs) = i
    end do 
  end if
  bla = sum(thr_hlper((min(c_nhier+1,max(3,c_nhier+1-c_multires))-1):c_nhier,1))
!$OMP SINGLE
  allocate(approxmst(bla))
!$OMP END SINGLE
!$OMP DO SCHEDULE(DYNAMIC,1)
  do bla2=1,bla
    bla3 = 0
    do ii=min(c_nhier+1,max(3,c_nhier+1-c_multires)),c_nhier+1
      bla3 = bla3 + thr_hlper(ii-1,1)
      if (bla3.ge.bla2) then
        mm = bla2 + thr_hlper(ii-1,1) - bla3
        exit
      end if
    end do
    mm = csnap2clus(mm,ii-1)
    approxmst(bla2)%deg = 0
    do tpx=1,tpn
      approxmst(bla2)%deg = approxmst(bla2)%deg + thr_btree(ii-1,tpx)%cls(mm)%nmbrs
    end do
    if (approxmst(bla2)%deg.gt.0) then
      allocate(approxmst(bla2)%adj(approxmst(bla2)%deg))
      jj = 0
      do tpx=1,tpn
        if (thr_btree(ii-1,tpx)%cls(mm)%nmbrs.gt.0) then
          approxmst(bla2)%adj((jj+1):(jj+thr_btree(ii-1,tpx)%cls(mm)%nmbrs)) = &
 &                  thr_btree(ii-1,tpx)%cls(mm)%snaps(1:thr_btree(ii-1,tpx)%cls(mm)%nmbrs)
          jj = jj + thr_btree(ii-1,tpx)%cls(mm)%nmbrs
        end if
      end do
    end if
  end do 
!$OMP BARRIER
  do ii=min(c_nhier+1,max(3,c_nhier+1-c_multires)),c_nhier+1
    deallocate(thr_btree(ii,tpi)%cls)
    thr_btree(ii,tpi)%ncls = 0
    thr_btree(ii,tpi)%nclsalsz = 0
  end do
!$OMP DO SCHEDULE(DYNAMIC,1)
  do bla2=1,bla
    bla3 = 0
    do ii=min(c_nhier+1,max(3,c_nhier+1-c_multires)),c_nhier+1
      bla3 = bla3 + thr_hlper(ii-1,1)
      if (bla3.ge.bla2) then
        mm = bla2 + thr_hlper(ii-1,1) - bla3
        exit
      end if
    end do
    mm = csnap2clus(mm,ii-1)
    do k=1,approxmst(bla2)%deg
!    if (birchtree(ii-1)%cls(mm)%nchildren.gt.10) write(*,*) tpi,birchtree(ii-1)%cls(mm)%nchildren,mm
      i = approxmst(bla2)%adj(k)
      mind = HUGE(mind)
      do j=1,birchtree(ii-1)%cls(mm)%nchildren  ! loop over clusters generated thus far by the snapshots in parent cluster
        ll = birchtree(ii-1)%cls(mm)%children(j)
        cnt1 = cnt1 + 1
        if ((thr_btree(ii,tpi)%cls(ll)%center.le.0).OR.(thr_btree(ii,tpi)%cls(ll)%center.gt.cstored)) call fexit()
        call snap_to_cluster_d(rdv,thr_btree(ii,tpi)%cls(ll),i)
        if (rdv.lt.mind) then
          mind = rdv
          jj = j
        end if
      end do
!     make a new cluster or not?
      if (mind.lt.scrcts(ii)) then
        kk = birchtree(ii-1)%cls(mm)%children(jj)
        call cluster_addsnap(thr_btree(ii,tpi)%cls(kk),i,rdv)
      else
        thr_btree(ii,tpi)%ncls = thr_btree(ii,tpi)%ncls + 1
        if (thr_btree(ii,tpi)%ncls.gt.thr_btree(ii,tpi)%nclsalsz) then
          call scluster_resizelst(thr_btree(ii,tpi)%nclsalsz,thr_btree(ii,tpi)%cls)
        end if
        call cluster_addsnap(thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls),i,rdv)
        call cluster_addchild(birchtree(ii-1)%cls(mm),mm,thr_btree(ii,tpi)%cls(thr_btree(ii,tpi)%ncls),thr_btree(ii,tpi)%ncls)
        cnt2 = cnt2 + 1
        kk = thr_btree(ii,tpi)%ncls
      end if
    end do
  end do
!$OMP BARRIER
!$OMP SINGLE
  deallocate(csnap2clus)
  deallocate(thr_hlper)
  deallocate(thr_rhlper)
  deallocate(approxmst)
!$OMP END SINGLE NOWAIT
!
  do ii=c_nhier+1,min(c_nhier+1,max(3,c_nhier+1-c_multires)),-1
    if (tpi.le.1) then
      bla = 1
    else
      bla = sum(thr_btree(ii,1:(tpi-1))%ncls) + 1
    end if
    bla2 = sum(thr_btree(ii,1:tpn)%ncls)
!$OMP BARRIER
!$OMP SINGLE
    if (allocated(birchtree(ii)%cls).EQV..true.) deallocate(birchtree(ii)%cls)
    allocate(birchtree(ii)%cls(bla2))
    birchtree(ii)%ncls = bla2
    birchtree(ii)%nclsalsz = bla2
!$OMP END SINGLE
    do i=bla,bla+thr_btree(ii,tpi)%ncls-1
      call copy_cluster(thr_btree(ii,tpi)%cls(i-bla+1),birchtree(ii)%cls(i))
    end do
    if (allocated(thr_btree(ii,tpi)%cls).EQV..true.) deallocate(thr_btree(ii,tpi)%cls)
  end do    
!
!$OMP BARRIER
!$OMP CRITICAL(DEVALCNT_UP)
  if (thrdat%verbosity.gt.1) write(ithread,*) tpi,' did ',cnt1
  cdevalcnt = cdevalcnt + cnt1
!$OMP END CRITICAL(DEVALCNT_UP)
  if (tpi.eq.1) then
    deallocate(thr_btree)
  end if
!$OMP BARRIER
  if (tpi.le.1) then
    call System_Clock(count=tmvars(2))
    if (thrdat%verbosity.gt.0) write(ithread,*) (tmvars(2)-tmvars(1))/thrdat%rate,' s for Clustering #3'
    tmvars(1) = tmvars(2)
  end if

 66 format('Level    # Clusters     Threshold     Total Snaps    Total Children')
 67 format(i9,i10,1x,g14.4,4x,i12,4x,i12)
 68 format(i9,i10,5x,a7,7x,i12,4x,i12)
  if (((modei.gt.0).OR.(modei.eq.-2)).AND.(tpi.le.1)) then
    write(ilog,66)
    write(ilog,68) c_nhier+1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
 &               sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
    do i=2,c_nhier+1
      write(ilog,67) c_nhier+2-i,birchtree(i)%ncls,scrcts(i),sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
 &               sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren)
    end do
    write(ilog,*) '---------------------------------------------------------------------'
    write(ilog,*)
    write(ilog,*) '... done after a total of ',cdevalcnt,' distance evaluations.'
    write(ilog,*)
  end if
!
  if (refine_clustering.EQV..true.) then
!$OMP BARRIER
!$OMP MASTER
    allocate(csnap2tree(maxval(birchtree(2:(c_nhier+1))%ncls)))
    do j=1,birchtree(c_nhier+1)%ncls
      call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
    end do
    call quality_of_clustering(birchtree(c_nhier+1)%ncls,birchtree(c_nhier+1)%cls,scrcts(c_nhier+1),qualmet)
    if (modei.gt.0) then
      write(ilog,*) 'Now merging clusters that yield joint reduced average intracluster distance ...'
    end if
 77 format('Would join ',i5,' (',i6,') and ',i5,'(',i6,') from: ',/,'Diam: ',g10.4,' / ',g10.4,'; Rad.: ',&
 &g10.4,' / ',g10.4,' to ',g10.4,' / ',g10.4,'.')
    cnt1 = 0
    cnt2 = 0
    do i=1,birchtree(c_nhier)%ncls
      nlst1 = 0
      do kkf=1,birchtree(c_nhier)%ncls
        if (i.eq.kkf) cycle
        call cluster_to_cluster_d(rdv,birchtree(c_nhier)%cls(i),birchtree(c_nhier)%cls(kkf))
        cnt1 = cnt1 + 1
        if (rdv.lt.scrcts(c_nhier)) then
          nlst1 = nlst1 + 1
          csnap2tree(nlst1) = kkf
        end if
      end do
      do j=1,birchtree(c_nhier)%cls(i)%nchildren
        jj = birchtree(c_nhier)%cls(i)%children(j)
        if (birchtree(c_nhier+1)%cls(jj)%nmbrs.le.0) cycle
        do kkf=1,nlst1
          do k=1,birchtree(c_nhier)%cls(csnap2tree(kkf))%nchildren
            kk = birchtree(c_nhier)%cls(csnap2tree(kkf))%children(k)
            if (birchtree(c_nhier+1)%cls(kk)%nmbrs.le.0) cycle
            if (jj.eq.kk) cycle
            call clusters_joint_diam(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk),rdv,helper)
            cnt1 = cnt1 + 1
            normer(1) = 0.5*birchtree(c_nhier+1)%cls(jj)%nmbrs*(birchtree(c_nhier+1)%cls(jj)%nmbrs-1.0)
            normer(2) = 0.5*birchtree(c_nhier+1)%cls(kk)%nmbrs*(birchtree(c_nhier+1)%cls(kk)%nmbrs-1.0)
            normer(3) = 1.0/(1.0*(birchtree(c_nhier+1)%cls(jj)%nmbrs + birchtree(c_nhier+1)%cls(kk)%nmbrs))
            normer(4) = 0.0
            if (sum(normer(1:2)).gt.0) normer(4) = 1.0/sum(normer(1:2))
            maxd(1) =  normer(4)*(normer(1)*birchtree(c_nhier+1)%cls(jj)%diam + &
              &                   normer(2)*birchtree(c_nhier+1)%cls(kk)%diam)
            maxd(2) =  normer(3)*(birchtree(c_nhier+1)%cls(jj)%nmbrs*birchtree(c_nhier+1)%cls(jj)%radius + &
   &                              birchtree(c_nhier+1)%cls(kk)%nmbrs*birchtree(c_nhier+1)%cls(kk)%radius)
            if ((helper.le.maxd(2)).OR.(rdv.le.maxd(1))) then
              if (birchtree(c_nhier+1)%cls(jj)%nmbrs.gt.birchtree(c_nhier+1)%cls(kk)%nmbrs) then
                call join_clusters(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk))
                birchtree(c_nhier+1)%cls(jj)%diam = rdv
                birchtree(c_nhier+1)%cls(jj)%radius = helper
                cnt2 = cnt2 + 1
              else
                call join_clusters(birchtree(c_nhier+1)%cls(kk),birchtree(c_nhier+1)%cls(jj))
                birchtree(c_nhier+1)%cls(kk)%diam = rdv
                birchtree(c_nhier+1)%cls(kk)%radius = helper
                cnt2 = cnt2 + 1
                exit
              end if
            end if
          end do
          if (birchtree(c_nhier+1)%cls(jj)%nmbrs.eq.0) exit
        end do
      end do
    end do
    if (modei.gt.0) then
      write(ilog,*) '... done after a total of ',cnt2,' merges requiring ',cnt1,' additional &
 &distance or joint size evaluations.'
      write(ilog,*)
    end if
    deallocate(csnap2tree)
!$OMP END MASTER
!$OMP BARRIER
    if (tpi.le.1) then
      call System_Clock(count=tmvars(2))
      if (thrdat%verbosity.gt.0) write(ithread,*) (tmvars(2)-tmvars(1))/thrdat%rate,' s for Clustering Refinement'
    end if
  end if
!
  if (tpi.gt.0) then
    do k=c_nhier+1,2,-1
      jj = sum(birchtree(k)%cls(1:birchtree(k)%ncls)%nmbrs)
      call clusters_postproc_threads(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls,jj,afalse,tpi,scrcts(k))
!$OMP SINGLE
      do while (birchtree(k)%cls(birchtree(k)%ncls)%nmbrs.le.0)
        birchtree(k)%ncls = birchtree(k)%ncls - 1
      end do
!$OMP END SINGLE
    end do
  else
    do k=c_nhier+1,2,-1
      call clusters_shorten(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls)
      do j=1,birchtree(k)%ncls
        call cluster_calc_params(birchtree(k)%cls(j),scrcts(k))
      end do
      do i=1,birchtree(k)%ncls
        call cluster_getcenter(birchtree(k)%cls(i))
      end do
      call clusters_sort(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls,afalse)
    end do
  end if
!$OMP BARRIER
  if (tpi.le.1) then
    call System_Clock(count=tmvars(2))
    if (thrdat%verbosity.gt.0) write(ithread,*) (tmvars(2)-tmvars(1))/thrdat%rate,' s for Clustering Post-Processing'
  end if
!
  if (modei.ge.0) then
!$OMP SINGLE
!   lastly, copy into global cluster array
    allocate(scluster(birchtree(c_nhier+1)%ncls))
!$OMP END SINGLE
    do i=max(1,tpi),birchtree(c_nhier+1)%ncls,tpn
      call copy_cluster(birchtree(c_nhier+1)%cls(i),scluster(i))
    end do
  end if
!$OMP CRITICAL(NNODESSET)
  nnodes = birchtree(c_nhier+1)%ncls
!$OMP END CRITICAL(NNODESSET)
!
 63 format(i7,1x,i7,1x,i8,1000(1x,g12.5))
 64 format(1000(g12.5,1x))
 69 format('------------- CLUSTER SUMMARY (THRESHOLD OF ',g12.5,')')
!
  if (modei.gt.0) then
!$OMP BARRIER
!$OMP MASTER
    do k=c_nhier+1,max(2,c_nhier+1-c_multires),-1
      write(ilog,69) scrcts(k)
      write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
!    do i=1,birchtree(c_nhier+1)%ncls
!      write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
!    end do
      do i=1,birchtree(k)%ncls
        write(ilog,63) i,birchtree(k)%cls(i)%nmbrs,birchtree(k)%cls(i)%center,birchtree(k)%cls(i)%diam,birchtree(k)%cls(i)%radius
      end do
      write(ilog,*) '----------------------------------------------------'
      write(ilog,*)
      call quality_of_clustering(birchtree(k)%ncls,birchtree(k)%cls(1:birchtree(k)%ncls),scrcts(k),qualmet)
    end do
!
!$OMP END MASTER
!$OMP BARRIER
    atrue = .true.
    call gen_graph_from_clusters(scluster,nnodes,atrue,tpi)
!$OMP BARRIER
!$OMP MASTER
    flush(ilog)
    call graphml_helper_for_clustering(scluster,nnodes)
    call vmd_helper_for_clustering(scluster,nnodes)
    flush(ilog)
!$OMP END MASTER
!$OMP BARRIER
  end if
!
  deallocate(scrcts)
!$OMP SINGLE
  cdevalcnt = 0
!$OMP END SINGLE
!
  if (tpi.le.1) then
    call System_Clock(count=tmvars(2))
    if (thrdat%verbosity.gt.0) write(ithread,*) (tmvars(2)-tmvars(1))/thrdat%rate,' s for Clustering Rest'
  end if
!
end
!
#endif
!
!----------------------------------------------------------------------------
!
! !!!!!!!!!!!!! WARNING !!!!!!!!!!!!
! !!! DO NOT USE THIS SUBROUTINE !!!
! !!!!!!!!!!!!! WARNING !!!!!!!!!!!!
!
! the old workflow, which has different exception handling and therefore gives slightly different
! results; safe for multi-threaded execution but no actual parallelism inside except what may
! be implemented in gen_graph_from_clusters(...) 
!
! set modei to 0 to obtain a silent version without post-processing
! set modei to -1 to obtain a silent version without post-processing and without populating scluster
!
subroutine birch_clustering_old(modei,nnodes,tpi)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer, INTENT(IN):: tpi,modei
!
  integer i,j,k,ii,jj,kk,ll,mm,atwo,fail,kkf,thekk,nlst1,nlst2,snapstart,snapend,snapinc,nnodes,cnhmax
  integer, ALLOCATABLE:: kklst(:,:),errcnt(:),kkhistory(:)
  integer(KIND=8) cnt1,cnt2
  RTYPE, ALLOCATABLE:: scrcts(:)
  RTYPE rdv,mind,helper,maxd(2),normer(4),qualmet(4)
  logical atrue,afalse,notdone,cycit
!
 33 format('ERROR B: ',i6,' could be part of ',i5,' at ',g14.7,' (last d',g14.6,')')
 34 format('NEXT HIGHER: ',i6,' could have been part of ',i5,' at ',g14.7,'.')
 35 format('ERROR C: ',i6,' is meant to a child of ',i5,' at level ',i5,' but has distance ',g14.7,'.')
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  atrue = .true.
  afalse = .false.
  atwo = 2
  if (cmaxrad.le.cradius) cmaxrad = 2.0*cradius
!
  if (cleadermode.le.2) then
    snapstart = 1
    snapend = cstored
    snapinc = 1
  else
    snapstart = cstored
    snapend = 1
    snapinc = -1
  end if
!
  allocate(kklst(cstored,2))
  allocate(kkhistory(c_nhier+1))
  allocate(errcnt(c_nhier+1))
  allocate(birchtree(c_nhier+1))
  do ii=1,c_nhier+1
    allocate(birchtree(ii)%cls(10))
    birchtree(ii)%ncls = 0
    birchtree(ii)%nclsalsz = 10
  end do
  allocate(scrcts(c_nhier+1))
  errcnt(:) = 0
  birchtree(1)%ncls = 1
  scrcts(c_nhier+1) = cradius
  if (c_nhier.gt.1) then
    scrcts(2) = cmaxrad
    do i=3,c_nhier
      scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - cradius) ! linear so far
    end do
  end if
  do i=1,c_nhier+1
    birchtree(i)%cls(:)%nmbrs = 0
    birchtree(i)%cls(:)%nb = 0
    birchtree(i)%cls(:)%nchildren = 0
    birchtree(i)%cls(:)%alsz = 0
    birchtree(i)%cls(:)%nbalsz = 0
    birchtree(i)%cls(:)%chalsz = 0
    birchtree(i)%cls(:)%parent = 0
  end do

  allocate(birchtree(1)%cls(1)%snaps(2))
  birchtree(1)%cls(1)%alsz = 2
  allocate(birchtree(1)%cls(1)%children(2))
  birchtree(1)%cls(1)%chalsz = 2
!
  if (modei.gt.0) then
    write(ilog,*)
    write(ilog,*) 'Now performing tree-based clustering ...'
  end if
  cycit = .false.
  cnt1 = 0
  cnt2 = 0
  do i=snapstart,snapend,snapinc
    kk = 1
    notdone = .true.
    fail = -1
    nlst1 = 1
    kklst(1,1) = kk
!    cnt1 = 0
    do ii=2,c_nhier
      jj = -1
      mind = HUGE(mind)
      nlst2 = 0
      do mm=1,nlst1
        kk = kklst(mm,1)
        do j=1,birchtree(ii-1)%cls(kk)%nchildren
          ll = birchtree(ii-1)%cls(kk)%children(j)
          cnt1 = cnt1 + 1
          if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
          call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
          if (rdv.lt.mind) then
            mind = rdv
            jj = j
            thekk = kk
          end if
          if (rdv.lt.scrcts(ii)) then
            nlst2 = nlst2 + 1
            kklst(nlst2,2) = ll
          end if
        end do
      end do
!     store the path
      if (jj.eq.-1) then
        kkhistory(ii-1) = 1
      else
        kkhistory(ii-1) = thekk
      end if
      if ((ii.le.(c_nhier-1)).AND.(jj.gt.0)) then
        if (nlst2.le.0) then ! absolutely nothing nearby 
          nlst1 = 1
          kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
          if (fail.eq.-1) then
!            write(*,*) i,' failed at ',ii,' w/ ',mind
            fail = ii
            kkf = thekk
          end if
        else
          if (fail.gt.0) then
            do j=fail,ii-1
              birchtree(j)%ncls = birchtree(j)%ncls + 1
              if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
              kkhistory(j) = birchtree(j)%ncls
              if (j.gt.fail) then
                call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
 &                                    birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                cnt2 = cnt2 + 1
              end if
            end do
            call cluster_addchild(birchtree(fail-1)%cls(kkf),kkf,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
            call cluster_addchild(birchtree(ii-1)%cls(birchtree(ii-1)%ncls),&
 &                                birchtree(ii)%cls(jj)%parent,birchtree(ii)%cls(jj),jj)
            fail = -1
          end if
          nlst1 = 1
          kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
        end if
        cycit = .true.
!     leaf
      else if ((ii.eq.c_nhier).AND.(mind.lt.scrcts(ii))) then
        kk = birchtree(ii-1)%cls(thekk)%children(jj)
        call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
        do j=2,c_nhier
          call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),i,rdv)
        end do
        notdone = .false.
      else
        if (fail.eq.-1) then
          fail = ii
          kkf = kk
        end if
        kk = kkf
        do j=fail,c_nhier
          birchtree(j)%ncls = birchtree(j)%ncls + 1
          if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
          call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
          if (j.gt.fail) then
            call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
 &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
            cnt2 = cnt2 + 1
          end if
        end do
        call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
        cnt2 = cnt2 + 1
        do j=2,fail
          call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),i,rdv)
        end do
        call snap_to_cluster_d(maxd(1),birchtree(fail-1)%cls(kk),i)
        notdone = .false.
      end if
      if (cycit.EQV..true.) then
        cycit = .false.
        cycle
      end if
      if (notdone.EQV..false.) exit
    end do
  end do
!
!  multi-pass 
  do cnhmax=c_nhier+1,min(c_nhier+1,max(3,c_nhier+1-c_multires)),-1
    if (cnhmax.lt.(c_nhier+1)) then
      birchtree(cnhmax)%ncls = 0
      birchtree(cnhmax)%cls(:)%nmbrs = 0
      birchtree(cnhmax)%cls(:)%parent = 0
      birchtree(cnhmax-1)%cls(:)%nchildren = 0
    end if
!
    cycit = .false.
    do i=snapstart,snapend,snapinc
      kk = 1
      notdone = .true.
      fail = -1
      nlst1 = 1
      kklst(1,1) = kk
      do ii=2,cnhmax
        jj = -1
        mind = HUGE(mind)
        nlst2 = 0
        do mm=1,nlst1
          kk = kklst(mm,1)
          do j=1,birchtree(ii-1)%cls(kk)%nchildren
            ll = birchtree(ii-1)%cls(kk)%children(j)
            cnt1 = cnt1 + 1
            if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
            call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
            if (rdv.lt.mind) then
              mind = rdv
              jj = j
              thekk = kk
            end if
            if (rdv.lt.scrcts(ii)) then
              nlst2 = nlst2 + 1
              kklst(nlst2,2) = ll
            end if
          end do
        end do
!       store the path
        if (jj.eq.-1) then
          kkhistory(ii-1) = 1
        else
          kkhistory(ii-1) = thekk
        end if
        if ((ii.lt.cnhmax).AND.(jj.gt.0)) then
          if (nlst2.le.0) then ! absolutely nothing nearby 
            nlst1 = 1
            kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
            if (fail.eq.-1) then
!              write(*,*) i,' failed at ',ii,' w/ ',mind
              fail = ii
              kkf = thekk
            end if
          else 
            nlst1 = 1
            kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
          end if
          cycit = .true.
!       leaf
        else if ((ii.eq.cnhmax).AND.(mind.lt.scrcts(ii))) then
          kk = birchtree(ii-1)%cls(thekk)%children(jj)
          call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
          notdone = .false.
        else if (fail.gt.0) then
          kk = kkf
          do j=fail,cnhmax
            birchtree(j)%ncls = birchtree(j)%ncls + 1
            if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
            call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
            if (j.gt.fail) then
              call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
   &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                       cnt2 = cnt2 + 1
            end if
          end do
          call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
          cnt2 = cnt2 + 1
          notdone = .false.
        else
          fail = ii
          j = cnhmax
          birchtree(j)%ncls = birchtree(j)%ncls + 1
          if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
          call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
          call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
          notdone = .false.
        end if
        if (cycit.EQV..true.) then
          cycit = .false.
          cycle
        end if
        if (notdone.EQV..false.) exit
      end do
    end do
  end do

 66 format('Level    # Clusters     Threshold     Total Snaps    Total Children')
 67 format(i9,i10,1x,g14.4,4x,i12,4x,i12)
 68 format(i9,i10,5x,a7,7x,i12,4x,i12)
  if (modei.gt.0) then
    write(ilog,66)
    write(ilog,68) c_nhier+1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
 &               sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
    do i=2,c_nhier+1
      write(ilog,67) c_nhier+2-i,birchtree(i)%ncls,scrcts(i),sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
 &               sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren)
    end do
    write(ilog,*) '---------------------------------------------------------------------'
    write(ilog,*)
    write(ilog,*) '... done after a total of ',cnt1,' distance evaluations.'
    write(ilog,*)
  end if
!
  if (refine_clustering.EQV..true.) then
!
    do j=1,birchtree(c_nhier+1)%ncls
      call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
    end do
    call quality_of_clustering(birchtree(c_nhier+1)%ncls,birchtree(c_nhier+1)%cls,scrcts(c_nhier+1),qualmet)
    if (modei.gt.0) then
      write(ilog,*) 'Now merging clusters that yield joint reduced average intracluster distance ...'
    end if
 77 format('Would join ',i5,' (',i6,') and ',i5,'(',i6,') from: ',/,'Diam: ',g10.4,' / ',g10.4,'; Rad.: ',&
 &g10.4,' / ',g10.4,' to ',g10.4,' / ',g10.4,'.')
    cnt1 = 0
    cnt2 = 0
    do i=1,birchtree(c_nhier)%ncls
      nlst1 = 0
      do kkf=1,birchtree(c_nhier)%ncls
        if (i.eq.kkf) cycle
        call cluster_to_cluster_d(rdv,birchtree(c_nhier)%cls(i),birchtree(c_nhier)%cls(kkf))
        cnt1 = cnt1 + 1
        if (rdv.lt.scrcts(c_nhier)) then
          nlst1 = nlst1 + 1
          kklst(nlst1,1) = kkf
        end if
      end do
      do j=1,birchtree(c_nhier)%cls(i)%nchildren
        jj = birchtree(c_nhier)%cls(i)%children(j)
        if (birchtree(c_nhier+1)%cls(jj)%nmbrs.le.0) cycle
        do kkf=1,nlst1
          do k=1,birchtree(c_nhier)%cls(kklst(kkf,1))%nchildren
            kk = birchtree(c_nhier)%cls(kklst(kkf,1))%children(k)
            if (birchtree(c_nhier+1)%cls(kk)%nmbrs.le.0) cycle
            if (jj.eq.kk) cycle
            call clusters_joint_diam(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk),rdv,helper)
            cnt1 = cnt1 + 1
            normer(1) = 0.5*birchtree(c_nhier+1)%cls(jj)%nmbrs*(birchtree(c_nhier+1)%cls(jj)%nmbrs-1.0)
            normer(2) = 0.5*birchtree(c_nhier+1)%cls(kk)%nmbrs*(birchtree(c_nhier+1)%cls(kk)%nmbrs-1.0)
            normer(3) = 1.0/(1.0*(birchtree(c_nhier+1)%cls(jj)%nmbrs + birchtree(c_nhier+1)%cls(kk)%nmbrs))
            normer(4) = 0.0
            if (sum(normer(1:2)).gt.0) normer(4) = 1.0/sum(normer(1:2))
            maxd(1) =  normer(4)*(normer(1)*birchtree(c_nhier+1)%cls(jj)%diam + &
              &                   normer(2)*birchtree(c_nhier+1)%cls(kk)%diam)
            maxd(2) =  normer(3)*(birchtree(c_nhier+1)%cls(jj)%nmbrs*birchtree(c_nhier+1)%cls(jj)%radius + &
   &                              birchtree(c_nhier+1)%cls(kk)%nmbrs*birchtree(c_nhier+1)%cls(kk)%radius)
            if ((helper.le.maxd(2)).OR.(rdv.le.maxd(1))) then
              if (birchtree(c_nhier+1)%cls(jj)%nmbrs.gt.birchtree(c_nhier+1)%cls(kk)%nmbrs) then
                call join_clusters(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk))
                birchtree(c_nhier+1)%cls(jj)%diam = rdv
                birchtree(c_nhier+1)%cls(jj)%radius = helper
                cnt2 = cnt2 + 1
              else
                call join_clusters(birchtree(c_nhier+1)%cls(kk),birchtree(c_nhier+1)%cls(jj))
                birchtree(c_nhier+1)%cls(kk)%diam = rdv
                birchtree(c_nhier+1)%cls(kk)%radius = helper
                cnt2 = cnt2 + 1
                exit
              end if
            end if
          end do
          if (birchtree(c_nhier+1)%cls(jj)%nmbrs.eq.0) exit
        end do
      end do
    end do
    if (modei.gt.0) then
      write(ilog,*) '... done after a total of ',cnt2,' merges requiring ',cnt1,' additional &
 &distance or joint size evaluations.'
      write(ilog,*)
    end if
  end if
!
! now shorten list and resort
  do k=c_nhier+1,2,-1
    call clusters_shorten(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls)
    do j=1,birchtree(k)%ncls
      call cluster_calc_params(birchtree(k)%cls(j),scrcts(k))
    end do
    do i=1,birchtree(k)%ncls
      call cluster_getcenter(birchtree(k)%cls(i))
    end do
    call clusters_sort(birchtree(k)%cls(1:birchtree(k)%ncls),birchtree(k)%ncls,afalse)
  end do
!
  if (modei.ge.0) then
!   lastly, copy into global cluster array
    allocate(scluster(birchtree(c_nhier+1)%ncls))
    do i=1,birchtree(c_nhier+1)%ncls
      call copy_cluster(birchtree(c_nhier+1)%cls(i),scluster(i))
    end do
  end if
!
  nnodes = birchtree(c_nhier+1)%ncls
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
 63 format(i7,1x,i7,1x,i8,1000(1x,g12.5))
 64 format(1000(g12.5,1x))
 69 format('------------- CLUSTER SUMMARY (THRESHOLD OF ',g12.5,')')
!
  if (modei.gt.0) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
    do k=c_nhier+1,c_nhier+1-c_multires,-1
      write(ilog,69) scrcts(k)
      write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
!    do i=1,birchtree(c_nhier+1)%ncls
!      write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
!    end do
      do i=1,birchtree(k)%ncls
        write(ilog,63) i,birchtree(k)%cls(i)%nmbrs,birchtree(k)%cls(i)%center,birchtree(k)%cls(i)%diam,birchtree(k)%cls(i)%radius
      end do
      write(ilog,*) '----------------------------------------------------'
      write(ilog,*)
      call quality_of_clustering(birchtree(k)%ncls,birchtree(k)%cls(1:birchtree(k)%ncls),scrcts(k),qualmet)
    end do
!
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    atrue = .true.
    call gen_graph_from_clusters(scluster,nnodes,atrue,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call graphml_helper_for_clustering(scluster,nnodes)
    call vmd_helper_for_clustering(scluster,nnodes)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  deallocate(kklst)
  deallocate(kkhistory)
  deallocate(scrcts)
  deallocate(errcnt)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif 
!
end subroutine birch_clustering_old
!
!----------------------------------------------------------------------------
!
! not thread-parallelized but can be called by multiple threads (relevant for subcalls) if nnodes is globally visible
!
subroutine leader_clustering(mode,nnodes,tpi)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer, INTENT(IN):: mode,tpi
!
  integer i,j,ii,jj,k,kk,l,nsets,atwo,nsets_old,globi,nclalcsz,azero
  integer snapstart,snapend,snapinc,nnodes
  integer(KIND=8) cnt1,cnt2
  integer, ALLOCATABLE:: iv1(:)
  RTYPE rdv,coreval,qualmet(4)
  logical atrue,afalse
!
  atrue = .true.
  afalse = .false.
  azero = 0
  atwo = 2
  cnt1 = 0
  if (cleadermode.le.2) then
    snapstart = 1
    snapend = cstored
    snapinc = 1
  else
    snapstart = cstored
    snapend = 1
    snapinc = -1
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  write(ilog,*)
!
  if (mode.eq.1) then
    write(ilog,*) 'Now performing LEADER-clustering ...'
    nclalcsz = 10
    allocate(scluster(nclalcsz))
    scluster(:)%nmbrs = 0
    scluster(:)%alsz = 0
    scluster(:)%nb = 0
    scluster(:)%nchildren = 0
    scluster(:)%nbalsz = 0
    scluster(:)%chalsz = 0
    scluster(:)%parent = 0
    k = 0 
    if ((cleadermode.eq.1).OR.(cleadermode.eq.3)) then
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=k,1,-1
          cnt1 = cnt1 + 1
          call snap_to_snap_d(rdv,i,scluster(j)%center)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
   else
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=1,k
          cnt1 = cnt1 + 1
          call snap_to_snap_d(rdv,i,scluster(j)%center)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
    end if
    write(ilog,*) 'Done after evaluating ',cnt1,' pairwise distances.'
    write(ilog,*)
!
    nsets = k
!
!   now sort
    call clusters_sort(scluster,nsets,afalse)
    do i=1,nsets
      call cluster_sortsnaps(scluster(i))
      call cluster_calc_params(scluster(i),cradius)
    end do
!
  else if (mode.eq.2) then
!
    write(ilog,*) 'Now performing first stage of modified LEADER-clustering ...'
    nclalcsz = 10
    allocate(scluster(nclalcsz))
    scluster(:)%nmbrs = 0
    scluster(:)%alsz = 0
    scluster(:)%nb = 0
    scluster(:)%nchildren = 0
    scluster(:)%nbalsz = 0
    scluster(:)%chalsz = 0
    scluster(:)%parent = 0
    k = 0
    if ((cleadermode.eq.1).OR.(cleadermode.eq.3)) then
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=k,1,-1 !max(k-100,1),-1
          cnt1 = cnt1 + 1
          call snap_to_cluster_d(rdv,scluster(j),i)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
    else
      do i=snapstart,snapend,snapinc
        ii = -1
        do j=1,k
          cnt1 = cnt1 + 1
          call snap_to_cluster_d(rdv,scluster(j),i)
          if (rdv.lt.cradius) then
            call cluster_addsnap(scluster(j),i,rdv)
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          k = k + 1
          if (k.gt.nclalcsz) call scluster_resizelst(nclalcsz,scluster)
          call cluster_addsnap(scluster(k),i,rdv)
        end if
      end do
    end if
    write(ilog,*) 'Done after evaluating ',cnt1,' pairwise distances.'
    write(ilog,*)
!
    nsets = k
    nsets_old = k
!
    call clusters_sort(scluster,nsets,afalse)
    do i=1,nsets
      call cluster_sortsnaps(scluster(i))
      call cluster_calc_params(scluster(i),cradius)
    end do
!
    cnt1 = 0
    cnt2 = 0
!
    if (refine_clustering.EQV..true.) then
      allocate(iv1(cstored))
!     the first block is to remove overlap
      write(ilog,*) 'Now removing cluster overlap by reassigning ambiguous frames to larger clusters ...'
!     we're doing this exactly twice - iteration would not be guaranteed to converge due to moving centers
      do globi=1,2
        do i=nsets,1,-1
          if (allocated(scluster(i)%snaps).EQV..false.) cycle
          do j=1,nsets
            if (i.eq.j) cycle
            if (allocated(scluster(j)%snaps).EQV..false.) cycle
            k = i
            l = j
            if (scluster(i)%nmbrs.gt.scluster(j)%nmbrs) then
              k = j
              l = i
            end if
            kk = 0
            call cluster_to_cluster_d(rdv,scluster(j),scluster(i))
            cnt1 = cnt1 + 1
            if (rdv.le.(scluster(i)%diam+scluster(j)%diam)) then
              do jj=1,scluster(k)%nmbrs
                if ((scluster(k)%snaps(jj).eq.scluster(k)%center).AND.(scluster(k)%nmbrs.gt.1)) cycle
                cnt2 = cnt2 + 1
                call snap_to_cluster_d(coreval,scluster(l),scluster(k)%snaps(jj))
                if (coreval.lt.cradius) then
                  kk = kk + 1
                  iv1(kk) = scluster(k)%snaps(jj)
                  if (coreval.gt.scluster(l)%diam) scluster(l)%diam = coreval
                end if
              end do
              if (kk.gt.0) then
                call cluster_transferframes(scluster(k),scluster(l),iv1(1:kk),kk)
              end if
            end if
            if (allocated(scluster(i)%snaps).EQV..false.) exit
          end do
        end do
      end do
!     now shorten list and resort
      call clusters_shorten(scluster(1:nsets),nsets)
      call clusters_sort(scluster(1:nsets),nsets,afalse)
      do i=nsets,1,-1
        call cluster_calc_params(scluster(i),cradius)
      end do
      deallocate(iv1)
      write(ilog,*) 'Done with ',nsets,' remaining of ',nsets_old,' original clusters using ',&
   & cnt1+cnt2,' additional pairwise distance evaluations.'
      write(ilog,*)
    end if
  end if
!
!
 66 format(i7,1x,i7,1x,i8,1x,g12.5,1x,g12.5,1x,g12.5)
 64 format(1000(g12.5,1x))
!
  write(ilog,*) '------------- CLUSTER SUMMARY ------------------'
  write(ilog,*) ' #       No.     Origin    Diameter     Radius      '
  do i=1,nsets
    write(ilog,66) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
  end do
  write(ilog,*) '------------------------------------------------'
  write(ilog,*)
  nnodes = nsets
  call quality_of_clustering(nnodes,scluster(1:nnodes),cradius,qualmet)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  call gen_graph_from_clusters(scluster,nnodes,atrue,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  call graphml_helper_for_clustering(scluster,nnodes)
  call vmd_helper_for_clustering(scluster,nnodes)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  end
!
!---------------------------------------------------------------------------------------------
!
! not thread-parallelized but can be called by multiple threads (relevant for subcalls) if nnodes is globally visible
!
subroutine hierarchical_clustering(nnodes,tpi)
!
  use clusters
  use mpistuff
  use iounit
  use math
  use interfaces
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer i,j,k,l,chosei,chosej,nclusters,nclalsz,aone,allnbs,globi,nnodes,azero
  integer, ALLOCATABLE:: tryat(:),inset(:),nbsnr(:),alllnks(:,:),iv1(:),iv2(:,:),iv3(:)
  real(KIND=4), ALLOCATABLE:: alldiss(:),tmpv(:)
  logical notdone,candid,atrue,afalse
  RTYPE rdv,cdiameter,qualmet(4),rdvtmp
!
  atrue = .true.
  afalse = .false.
  aone = 1
  azero = 0
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  allocate(tryat(cstored))
  allocate(inset(cstored))
  allocate(nbsnr(cstored))
  nclalsz = 10
  allocate(scluster(nclalsz))
  scluster(:)%nmbrs = 0
  scluster(:)%alsz = 0
  scluster(:)%nb = 0
  scluster(:)%nchildren = 0
  scluster(:)%nbalsz = 0
  scluster(:)%chalsz = 0
  scluster(:)%parent = 0
!
  tryat(:) = 1
  inset(:) = 0
  notdone = .true.
  nclusters = 0
  cdiameter = 2.0*cradius
  candid = .false.
  nbsnr(:) = cnblst(1:cstored)%nbs
!
! this first step prunes the nb-list to our desired size
  do i=1,cstored
    do j=1,cnblst(i)%nbs
      if (cnblst(i)%dis(j).gt.cdiameter) exit
    end do
    cnblst(i)%nbs = j-1
  end do
!
  allnbs = sum(cnblst(1:cstored)%nbs)
!
  if (allnbs.gt.0) then
    write(ilog,*)
    write(ilog,*) 'Now creating global sorted list of neighbor pairs ...'
    allocate(iv2(allnbs,2))
    allocate(alldiss(allnbs))
    allocate(tmpv(allnbs))
    allocate(iv1(allnbs))
    allocate(iv3(allnbs))
    j = 1
    do i=1,cstored
      if (cnblst(i)%nbs.gt.0) then
        iv2(j:j+cnblst(i)%nbs-1,1) = i
        iv2(j:j+cnblst(i)%nbs-1,2) = cnblst(i)%idx(1:cnblst(i)%nbs)
        tmpv(j:j+cnblst(i)%nbs-1) = cnblst(i)%dis(1:cnblst(i)%nbs)
        j = j + cnblst(i)%nbs
      end if
    end do
    do j=1,allnbs
      iv1(j) = j
    end do
    call merge_sort(ldim=allnbs,up=atrue,list=tmpv(1:allnbs),olist=alldiss(1:allnbs),&
 &                ilo=aone,ihi=allnbs,idxmap=iv1(1:allnbs),olist2=iv3(1:allnbs))
    deallocate(tmpv)
    deallocate(iv1)
    allocate(alllnks(allnbs,2))
    do i=1,allnbs
      alllnks(i,:) = iv2(iv3(i),:)
    end do
    deallocate(iv2)
    deallocate(iv3)
    write(ilog,*) '... done.'
    write(ilog,*)
  else
    notdone = .false.
  end if
!
  write(ilog,*) 'Now performing hierarchical clustering by considering shortest remaining link ...'
  globi = 1
  do while (notdone.EQV..true.)
    rdv = alldiss(globi)
    chosei = alllnks(globi,1)
    chosej = alllnks(globi,2)
    globi = globi + 1
    if (globi.gt.allnbs) exit
    if (rdv.gt.cdiameter) exit
    if ((inset(chosei).gt.0).AND.(inset(chosej).gt.0)) then
!     merge inset(chosei) and inset(chosej) if permissible
      if (inset(chosei).ne.inset(chosej)) then
        if (clinkage.eq.1) then ! maximum linkage
          candid = .true.
          do i=1,scluster(inset(chosei))%nmbrs
            do j=1,scluster(inset(chosej))%nmbrs
              k = 0
              do l=1,cnblst(scluster(inset(chosej))%snaps(j))%nbs
                if (inset(cnblst(scluster(inset(chosej))%snaps(j))%idx(l)).eq.inset(chosei)) k = k + 1
                if (k.eq.scluster(inset(chosei))%nmbrs) exit
              end do
              if (k.lt.scluster(inset(chosei))%nmbrs) then
                candid = .false.
                exit
              end if
            end do
            if (candid.EQV..false.) exit
          end do
        else if (clinkage.eq.2) then ! minimum linkage
          candid = .false.
          do i=1,scluster(inset(chosei))%nmbrs
            do j=1,scluster(inset(chosej))%nmbrs
              do l=1,cnblst(scluster(inset(chosej))%snaps(j))%nbs
                if (inset(cnblst(scluster(inset(chosej))%snaps(j))%idx(l)).eq.inset(chosei)) then
                  candid = .true.
                  exit
                end if
              end do
              if (candid.EQV..true.) exit
            end do
            if (candid.EQV..true.) exit
          end do
        else if (clinkage.eq.3) then ! mean linkage
          candid = .false.
          call cluster_to_cluster_d(rdvtmp,scluster(inset(chosei)),scluster(inset(chosej)))
          if (rdvtmp.le.cradius) candid = .true.
        end if
        if (candid.EQV..true.) then
          if (scluster(inset(chosei))%nmbrs.gt.scluster(inset(chosej))%nmbrs) then
            k = inset(chosei)
            l = inset(chosej)
          else
            k = inset(chosej)
            l = inset(chosei)
          end if
          inset(scluster(l)%snaps(1:scluster(l)%nmbrs)) = k
          call join_clusters(scluster(k),scluster(l))
!         transfer last to l
          if (l.lt.nclusters) then
            call copy_cluster(scluster(nclusters),scluster(l))
            deallocate(scluster(nclusters)%snaps)
            scluster(nclusters)%alsz = 0
            scluster(nclusters)%nmbrs = 0
            nclusters = nclusters - 1
            do j=1,scluster(l)%nmbrs
              if (inset(scluster(l)%snaps(j)).ne.nclusters+1) call fexit()
              inset(scluster(l)%snaps(j)) = l
            end do
          else if (nclusters.eq.l) then
            nclusters = nclusters - 1
          end if
        end if
      end if
!   append inset(chosei) if permissible
    else if (inset(chosei).gt.0) then
      if (clinkage.eq.1) then ! maximum linkage
        k = 0
        candid = .false.
        do i=1,cnblst(chosej)%nbs
          if (inset(cnblst(chosej)%idx(i)).eq.inset(chosei)) k = k + 1
          if (k.eq.scluster(inset(chosei))%nmbrs) exit
        end do
        if (k.eq.scluster(inset(chosei))%nmbrs) candid = .true.
      else if (clinkage.eq.2) then ! minimum
        candid = .false.
!       note that tryat(chosej) has to be 1, otherwise it inset(chosej) cannot be zero
        do i=1,cnblst(chosej)%nbs
          if (inset(cnblst(chosej)%idx(i)).eq.inset(chosei)) then
            candid = .true.
            exit
          end if
        end do
      else if (clinkage.eq.3) then ! mean linkage
        candid = .false.
        call snap_to_cluster_d(rdvtmp,scluster(inset(chosei)),chosej)
        if (rdvtmp.le.cradius) candid = .true.
      end if
      if (candid.EQV..true.) then
        inset(chosej) = inset(chosei)
        call cluster_addsnap(scluster(inset(chosei)),chosej,rdv)
      end if
!   append inset(chosej) if permissible
    else if (inset(chosej).gt.0) then
      if (clinkage.eq.1) then ! maximum linkage
        k = 0
        candid = .false.
        do i=1,cnblst(chosei)%nbs
          if (inset(cnblst(chosei)%idx(i)).eq.inset(chosej)) k = k + 1
          if (k.eq.scluster(inset(chosej))%nmbrs) exit
        end do
        if (k.eq.scluster(inset(chosej))%nmbrs) candid = .true.
      else if (clinkage.eq.2) then ! minimum
        candid = .false.
!       note that tryat(chosei) has to be 1, otherwise it inset(chosei) cannot be zero
        do i=1,cnblst(chosei)%nbs
          if (inset(cnblst(chosei)%idx(i)).eq.inset(chosej)) then
            candid = .true.
            exit
          end if
        end do
      else if (clinkage.eq.3) then ! mean linkage
        candid = .false.
        call snap_to_cluster_d(rdvtmp,scluster(inset(chosej)),chosei)
        if (rdvtmp.le.cradius) candid = .true.
      end if
      if (candid.EQV..true.) then
        inset(chosei) = inset(chosej)
        call cluster_addsnap(scluster(inset(chosej)),chosei,rdv)
      end if 
!   create new cluster of size 2
    else
      nclusters = nclusters + 1
      if (nclusters.gt.nclalsz) call scluster_resizelst(nclalsz,scluster)
      call cluster_addsnap(scluster(nclusters),chosei,rdv)
      call cluster_addsnap(scluster(nclusters),chosej,rdv)
      inset(chosei) = nclusters
      inset(chosej) = nclusters
    end if
!
    tryat(chosei) = tryat(chosei) + 1
!   update nb-list counter
    if (tryat(chosei).le.cnblst(chosei)%nbs) then
      do while ((inset(chosei).eq.inset(cnblst(chosei)%idx(tryat(chosei)))).AND.(inset(chosei).gt.0))
        tryat(chosei) = tryat(chosei) + 1
        if (tryat(chosei).gt.cnblst(chosei)%nbs) exit
      end do
    end if
    if (tryat(chosej).le.cnblst(chosej)%nbs) then
      if (cnblst(chosej)%idx(tryat(chosej)).eq.chosei) tryat(chosej) = tryat(chosej) + 1
      if (tryat(chosej).le.cnblst(chosej)%nbs) then
        do while ((inset(chosej).eq.inset(cnblst(chosej)%idx(tryat(chosej))).AND.inset(chosej).gt.0))
          tryat(chosej) = tryat(chosej) + 1
          if (tryat(chosej).gt.cnblst(chosej)%nbs) exit
        end do
      end if
    end if
  end do
  write(ilog,*) '... done.'
  if (allnbs.gt.0) then
    deallocate(alldiss)
    deallocate(alllnks)
  end if
  do i=1,cstored
    if (inset(i).le.0) then
      nclusters = nclusters + 1
      if (nclusters.gt.nclalsz) call scluster_resizelst(nclalsz,scluster)
      call cluster_addsnap(scluster(nclusters),i,rdv)
    end if
  end do
!
! shorten list
  call clusters_shorten(scluster(1:nclusters),nclusters)
  do i=1,nclusters
!   get centroid representative
    call cluster_getcenter(scluster(i))
  end do
  call clusters_sort(scluster(1:nclusters),nclusters,afalse)
!
 66 format(i7,1x,i7,1x,i8,1x,g12.5,1x,g12.5)
 65 format(1000(g12.5,1x))
!
  write(ilog,*) '------------- CLUSTER SUMMARY ------------------'
  write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
  do i=1,nclusters
!   get radius etc.
    call cluster_calc_params(scluster(i),cradius)
    write(ilog,66) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
  end do
  write(ilog,*) '------------------------------------------------'
  write(ilog,*)
!
  call quality_of_clustering(nclusters,scluster(1:nclusters),cradius,qualmet)
!
  nnodes = nclusters
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  if (nnodes.gt.0) call gen_graph_from_clusters(scluster,nnodes,atrue,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  if (nnodes.gt.0) call vmd_helper_for_clustering(scluster,nclusters)
  if (nnodes.gt.0) call graphml_helper_for_clustering(scluster,nclusters)
!
  cnblst(1:cstored)%nbs = nbsnr(:)
!
  deallocate(nbsnr)
  deallocate(tryat)
  deallocate(inset)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
end
!
!------------------------------------------------------------------------------------------
!
! this simply regenerates a clustering from cludata and the data in sconnect(:,4)
! cstored has to be correct
! it populates scluster(1:nba) and sets nba itself to the number of clusters
!
subroutine file_clustering(modei,nba,tpi)
!
  use clusters
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: tpi,modei
!
  integer i,k,nba,tpi2
  RTYPE rdv,qualmet(4)
  logical atrue,afalse
  integer, ALLOCATABLE:: fixlims(:,:)
#ifdef ENABLE_THREADS
  integer tpn,ixx
  integer OMP_GET_NUM_THREADS
!
  atrue = .true.
  afalse = .false.
!
!$OMP MASTER
#else
  atrue = .true.
  afalse = .false.
!
#endif
  if (allocated(sconnect).EQV..false.) then
    write(ilog,*) 'Fatal. Entered file_clustering(...) without the snapshot map ("sconnect") being allocated. This is a bug.'
    call fexit()
  end if
  nba = MAXVAL(sconnect(1:cstored,4))
  nstruccls = nba
!
  allocate(scluster(nba))
!
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
  if (tpi.gt.0) then
    tpn = OMP_GET_NUM_THREADS()
    allocate(fixlims(6,tpn))
    do i=1,tpn
      fixlims(1,i) = min(nba,nint((i-1)*(1.0*nba/(1.0*tpn)))) + 1
      fixlims(2,i) = max(1,nint(i*(1.0*nba/(1.0*tpn))))
      fixlims(3,i) = min(cstored,nint((i-1)*(1.0*cstored/(1.0*tpn)))) + 1
      fixlims(4,i) = max(1,nint(i*(1.0*cstored/(1.0*tpn))))
    end do
    do k=1,2
      ixx = 1
      do i=2,tpn
        if (fixlims(2*k,i).ge.fixlims(2*k-1,i)) then
          if (fixlims(2*k-1,i).le.fixlims(2*k,ixx)) then
            fixlims(2*k-1,i) = fixlims(2*k-1,i) + 1
          else
            ixx = i
          end if
        end if
      end do
    end do
    tpi2 = tpi
  else
#endif
  allocate(fixlims(6,1))
  tpi2 = 1
  fixlims(:,1) = 1
  fixlims(2,1) = nba
  fixlims(4,1) = cstored
#ifdef ENABLE_THREADS
  end if
#endif
!
  do i=fixlims(1,tpi2),fixlims(2,tpi2)
    scluster(i)%center = 0
    scluster(i)%geni = 0
    scluster(i)%genidx = 0
    scluster(i)%nchildren = 0
    scluster(i)%parent = 0
    scluster(i)%nbalsz = 0
    scluster(i)%nb = 0
    scluster(i)%nmbrs = 0
    scluster(i)%alsz = 0
  end do
!
  if (cmode.eq.6) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
    if (tpi.gt.0) then
      do i=1,cstored
        if (mod(sconnect(i,4),tpn).ne.(tpi-1)) cycle
        rdv = 1.0
        call cluster_addsnap(scluster(sconnect(i,4)),i,rdv)
      end do
    else
#endif
    do i=1,cstored
      rdv = 1.0
      call cluster_addsnap(scluster(sconnect(i,4)),i,rdv)
    end do
#ifdef ENABLE_THREADS
    end if
!$OMP BARRIER
#endif
  else if (cmode.eq.7) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
    if (tpi.gt.0) then
      do i=1,cstored
        if (mod(sconnect(i,4),tpn).ne.(tpi-1)) cycle
        rdv = 1.0
        call cluster_addjustsnap(scluster(sconnect(i,4)),i)
      end do
    else
#endif
    do i=1,cstored
      rdv = 1.0
      call cluster_addjustsnap(scluster(sconnect(i,4)),i)
    end do
#ifdef ENABLE_THREADS
    end if
!$OMP BARRIER
#endif
  end if

!
  k = cstored
  rdv = cradius
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    call clusters_postproc_threads(scluster,nba,k,afalse,tpi,rdv)
  else
#endif
  call clusters_shorten(scluster(1:nba),nba)
  do i=1,nba
    call cluster_calc_params(scluster(i),rdv)
  end do
  if (cmode.ne.7) then
    do i=1,nba
      call cluster_getcenter(scluster(i))
    end do
  end if
  call clusters_sort(scluster(1:nba),nba,afalse)
#ifdef ENABLE_THREADS
  end if
#endif
!
 63 format(i7,1x,i7,1x,i8,1000(1x,g12.5))
 64 format(1000(g12.5,1x))
 69 format('------------- CLUSTER SUMMARY (ASSUMED THRESHOLD OF ',g12.5,')')
!
  if (modei.gt.0) then
!$OMP BARRIER
!$OMP MASTER
    write(ilog,69) rdv
    write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
    do i=1,nba
      write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
    end do
    write(ilog,*) '----------------------------------------------------'
    write(ilog,*)
    call quality_of_clustering(nba,scluster(1:nba),rdv,qualmet)
!$OMP END MASTER
!$OMP BARRIER
  end if
  call gen_graph_from_clusters(scluster,nba,atrue,tpi)
!$OMP BARRIER
!$OMP MASTER
  call graphml_helper_for_clustering(scluster,nba)
  call vmd_helper_for_clustering(scluster,nba)
!$OMP END MASTER
!$OMP BARRIER
! 
end
!
!---------------------------------------------------------------------------------
!
! this subroutine generates an exact MST assuming it is provided with a nb-list
! object that holds all the necessary edges
! this routine is very memory-intensive due to the duplication of the already large nb-list object
! it is highly related to hierarchical clustering with minimum linkage and max threshold
!
subroutine gen_MST_from_nbl()
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer allnbs,i,j,k,aone,globi,chosei,chosej,ntrees,nlnks,forsz
  real(KIND=4), ALLOCATABLE:: alldiss(:),tmpv(:)
  RTYPE rdv
  integer, ALLOCATABLE:: iv2(:,:),iv1(:),iv3(:),alllnks(:,:)
  type(t_scluster), ALLOCATABLE:: it(:)
  logical atrue,notdone
!
  write(ilog,*)
  write(ilog,*) 'Now creating global sorted list of neighbor pairs ...'
!
  aone = 1
  atrue = .true.
  notdone = .true.
  allnbs = sum(cnblst(1:cstored)%nbs)/2
  allocate(iv2(allnbs,2))
  allocate(alldiss(allnbs))
  allocate(tmpv(allnbs))
  allocate(iv1(allnbs))
  allocate(iv3(allnbs))
  j = 0
  do i=1,cstored
    do k=1,cnblst(i)%nbs
      if (cnblst(i)%idx(k).gt.i) then
        j = j + 1
        iv2(j,1) = i
        iv2(j,2) = cnblst(i)%idx(k)
        tmpv(j) = cnblst(i)%dis(k)
      end if
    end do
  end do
  do j=1,allnbs
    iv1(j) = j
  end do
  call merge_sort(ldim=allnbs,up=atrue,list=tmpv(1:allnbs),olist=alldiss(1:allnbs),&
 &                ilo=aone,ihi=allnbs,idxmap=iv1(1:allnbs),olist2=iv3(1:allnbs))
  deallocate(tmpv)
  deallocate(iv1)
  allocate(alllnks(2,allnbs))
  do i=1,allnbs
    alllnks(:,i) = iv2(iv3(i),:)
!    write(*,*) alllnks(1:2,i),alldiss(i)
  end do
  deallocate(iv2)
  deallocate(iv3)
  write(ilog,*) '... done.'
  write(ilog,*)
!
  allocate(iv1(cstored))
  iv1(:) = 0
  write(ilog,*) 'Now generating MST by considering shortest remaining link and merging ...'
  globi = 1
  nlnks = 0
  ntrees = 0
  forsz = 10
  allocate(it(forsz))
  do i=1,forsz
    it(i)%alsz = 0
    it(i)%nmbrs = 0
  end do
  do while (notdone.EQV..true.)
    rdv = alldiss(globi)
    chosei = alllnks(1,globi)
    chosej = alllnks(2,globi)
    globi = globi + 1
    if (globi.gt.allnbs) exit
    if ((iv1(chosei).le.0).AND.(iv1(chosej).le.0)) then
      nlnks = nlnks + 1
      ntrees = ntrees + 1
      iv1(chosei) = ntrees
      iv1(chosej) = ntrees
      if (ntrees.gt.forsz) call scluster_resizelst(forsz,it)
      call cluster_addsnap(it(iv1(chosei)),chosei,0.0)
      call cluster_addsnap(it(iv1(chosei)),chosej,0.0)
    else if ((iv1(chosei).gt.0).AND.(iv1(chosej).gt.0)) then
      if (iv1(chosei).eq.iv1(chosej)) cycle
      nlnks = nlnks + 1
      if (it(iv1(chosei))%nmbrs.gt.it(iv1(chosej))%nmbrs) then
        j = iv1(chosej)
        do i=1,it(j)%nmbrs
          iv1(it(j)%snaps(i)) = iv1(chosei)
        end do
        call join_clusters(it(iv1(chosei)),it(j))
      else
        j = iv1(chosei)
        do i=1,it(j)%nmbrs
          iv1(it(j)%snaps(i)) = iv1(chosej)
        end do
        call join_clusters(it(iv1(chosej)),it(j))
      end if
    else if (iv1(chosei).gt.0) then
      nlnks = nlnks + 1
      iv1(chosej) = iv1(chosei)
      call cluster_addsnap(it(iv1(chosei)),chosej,0.0)
    else
      nlnks = nlnks + 1
      iv1(chosei) = iv1(chosej)
      call cluster_addsnap(it(iv1(chosej)),chosei,0.0)
    end if
    alllnks(1,nlnks) = chosei
    alllnks(2,nlnks) = chosej
    alldiss(nlnks) = rdv
    if (nlnks.eq.(cstored-1)) exit
  end do
!
  if (nlnks.ne.(cstored-1)) then
    write(ilog,*) 'Fatal. Neighbor list is insufficient to create minimum spanning tree. &
 &Increase relevant thresholds.'
    call fexit()
  end if
!
  deallocate(iv1)
  do i=1,forsz
    if (allocated(it(i)%snaps).EQV..true.) deallocate(it(i)%snaps)
    if (allocated(it(i)%tmpsnaps).EQV..true.) deallocate(it(i)%tmpsnaps)
    if (allocated(it(i)%sums).EQV..true.) deallocate(it(i)%sums)
    if (allocated(it(i)%map).EQV..true.) deallocate(it(i)%map)
    if (allocated(it(i)%children).EQV..true.) deallocate(it(i)%children)
    if (allocated(it(i)%wghtsnb).EQV..true.) deallocate(it(i)%wghtsnb)
    if (allocated(it(i)%lstnb).EQV..true.) deallocate(it(i)%lstnb)
    if (allocated(it(i)%flwnb).EQV..true.) deallocate(it(i)%flwnb)
  end do
  deallocate(it)
!
  allocate(approxmst(cstored))
  approxmst(1:cstored)%deg = 0
 587 format(' Weight of Minimum Spanning Tree: ',1x,g12.5,a)
  write(ilog,*)
  if (cdis_crit.le.2) then
    write(ilog,587) sum(alldiss(1:cstored-1)),' degrees'
  else if (cdis_crit.le.4) then
    write(ilog,587) sum(alldiss(1:cstored-1)),' '
  else if (cdis_crit.le.10) then
    write(ilog,587) sum(alldiss(1:cstored-1)),' Angstrom'
  end if
  write(ilog,*)
!
  call gen_MST(alllnks(:,1:(cstored-1)),alldiss(1:(cstored-1)),approxmst)
!
  deallocate(alldiss)
  deallocate(alllnks)
!
  write(ilog,*) '... done.'
  write(ilog,*)
!
end
!
!-----------------------------------------------------------------------------------------
!
! this subroutine generates a set of links and their lengths that constitute an approximate
! minimum spanning tree (SST) based on results from tree-based clustering in the birchtree object
! it then transcribes this into an adjacency list object
! the accuracy of the SST depends on the number of guesses (cprogindrmax), the properties
! of the clustering, and the auxiliary search depth being utilized (cprogrdepth)
!
subroutine gen_MST_from_treeclustering(tpi2)
!
  use clusters
  use iounit
  use interfaces
#ifdef ENABLE_THREADS
  use threads
  use keys
#else
  use threads, ONLY: thrdat
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi2
!
  integer N_NEARS
  parameter(N_NEARS=5)  ! number of guesses per snapshot to keep in memory
!
  integer i,j,k,ki,e,i1,i2,l,m,mm,ixx,ixx2,cursnp,thismode,b1,b2,mycl,tpi,cbnds(2)
  integer tlstsz,ishfx,kshfx
  integer(KIND=8) testcnt,trcnt,tdcnt
  integer(KIND=8) tmvars(2)!,tmm(12)
  logical atrue
  integer ntrees,oldntrees ! number of active trees
  integer, ALLOCATABLE:: tmpix(:,:)
  RTYPE, ALLOCATABLE:: tmpdis(:,:)
  RTYPE tmpcld(calcsz)!,tps(8)
  RTYPE jkdist ! distance between two snapshots
  integer boruvkasteps! counts how many boruvka steps are needed
  integer nmstedges   ! current number of edges in the growing trees
  integer, ALLOCATABLE:: testlst(:,:) ! temporary list of eligible candidates for random search with annotations
  integer, ALLOCATABLE:: snplst(:)    ! a temp array to hold lists of snapshots to scan for eligible distances
  integer a           ! numbers of snapshots, used for sophisticated nearest neighbor guessing
  integer kk,kix,kixi ! local snapshot index within cluster, used for sophisticated nearest neighbor guessing
  integer, ALLOCATABLE:: satisfied(:) ! array to indicate whether enough distances have been probed
  integer, ALLOCATABLE:: mstedges(:,:)! contains all the edges in the SST
                                      ! Indices: edgenumber (1-nedges), endpoint (1-2)
  real(KIND=4), ALLOCATABLE:: lmstedges(:) ! contains the lengths of all the edges in the SST
#ifdef ENABLE_THREADS
  integer, ALLOCATABLE:: fixlims(:,:)    ! boundaries for threads
  integer tpn                            ! thread index and number
  integer OMP_GET_NUM_THREADS
  logical OMP_IN_PARALLEL
  RTYPE, ALLOCATABLE:: rdbuf(:)          ! random number buffer
  integer rdbufsize,rdbufcnt,rdbufthresh ! associated parameters
!  type(t_rndstat) rngt
#else
  RTYPE random                           ! random number
#endif
!
!  tps(:) = 0.0
  atrue = .true.
  ntrees = cstored
  boruvkasteps = 0
  tlstsz = cstored
  testcnt = 0
  trcnt = 0
  tdcnt = 0
  nmstedges = 0
#ifdef ENABLE_THREADS
! the fact that threads "compete" for access to the PRNG means that one can no longer use
! a fixed seed to guarantee getting the same result
  rdbufsize = 4*cprogindrmax ! min(100000,cprogindrmax*cstored/max(1,thrdat%maxn))
  rdbufthresh = 4*cprogindrmax
#endif
!
#ifdef ENABLE_THREADS
  if (tpi2.le.0) then
    if (OMP_IN_PARALLEL().EQV..true.) then
      write(ilog,*) 'Fatal. When using multi-threaded code, gen_MST_from_treeclustering(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
      call fexit()
    end if
    tpn = 1
    allocate(fixlims(2,1))
    fixlims(2,1) = cstored
    fixlims(1,1) = 1
    tpi = 1
  else
    tpn = OMP_GET_NUM_THREADS()
    tpi = tpi2
    allocate(fixlims(2,tpn))
    do i=1,tpn
      fixlims(1,i) = min(cstored,nint((i-1)*(1.0*cstored/(1.0*tpn)))) + 1
      fixlims(2,i) = max(1,nint(i*(1.0*cstored/(1.0*tpn))))
    end do
    ixx = 1
    do i=2,tpn
      k=1
      if (fixlims(2*k,i).ge.fixlims(2*k-1,i)) then
        if (fixlims(2*k-1,i).le.fixlims(2*k,ixx)) then
          fixlims(2*k-1,i) = fixlims(2*k-1,i) + 1
        else
          ixx = i
        end if
      end if
    end do
    tlstsz = max(1,fixlims(2,tpi) - fixlims(1,tpi) + 1)
    i = -1
    k = 1
!    call init_anyprng(rngt,k,i)
  end if
!
!$OMP MASTER
  allocate(tmptree(ntrees))
  allocate(csnap2tree(cstored))
  allocate(csnap2clus(cstored,c_nhier+1))
  allocate(mstedges(2,cstored-1))
  allocate(lmstedges(cstored-1))

#else
!
  tpi = 1
  allocate(tmptree(ntrees))
  allocate(csnap2tree(cstored))
  allocate(csnap2clus(cstored,c_nhier+1))
  allocate(mstedges(2,cstored-1))
  allocate(lmstedges(cstored-1))
!
#endif
  write(ilog,*)
  write(ilog,*) 'Now generating approximate MST based on tree-based clustering ...'
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
!
  do i=fixlims(1,tpi),fixlims(2,tpi)
    allocate(tmptree(i)%snaps(1))
    tmptree(i)%nsnaps = 1
    tmptree(i)%mine(1) = -1
    tmptree(i)%mine(2) = -1
    tmptree(i)%mind = HUGE(jkdist)
    tmptree(i)%nsibalsz = 0
    tmptree(i)%snaps(1) = i 
    csnap2tree(i) = i
  end do
#else
  do i=1,ntrees
    allocate(tmptree(i)%snaps(1))
    tmptree(i)%snaps(1) = i 
    csnap2tree(i) = i
  end do
  tmptree(:)%nsnaps = 1
  tmptree(:)%mine(1) = -1
  tmptree(:)%mine(2) = -1
  tmptree(:)%mind = HUGE(jkdist)
  tmptree(:)%nsibalsz = 0
#endif
!
! Generate map (snapshot, tree level) -> cluster based on birchtree;
! in the process remove double entries for non-terminal levels and populate persistent tree index vector per cluster
#ifdef ENABLE_THREADS
!$OMP DO SCHEDULE(STATIC, 1)
#endif
  do l=c_nhier+1,1,-1
    csnap2clus(:,l) = 0
    do k=1,birchtree(l)%ncls
      if (allocated(birchtree(l)%cls(k)%tmpsnaps).EQV..true.) deallocate(birchtree(l)%cls(k)%tmpsnaps)
      allocate(birchtree(l)%cls(k)%tmpsnaps(birchtree(l)%cls(k)%nmbrs))
      do j=1,birchtree(l)%cls(k)%nmbrs
        birchtree(l)%cls(k)%tmpsnaps(j) = csnap2tree(birchtree(l)%cls(k)%snaps(j))
        if (csnap2clus(birchtree(l)%cls(k)%snaps(j),l).gt.0) then
          if (birchtree(l)%cls(csnap2clus(birchtree(l)%cls(k)%snaps(j),l))%nmbrs.gt.birchtree(l)%cls(k)%nmbrs) then
            call cluster_removesnap(birchtree(l)%cls(k),birchtree(l)%cls(k)%snaps(j))
          else
            call cluster_removesnap(birchtree(l)%cls(csnap2clus(birchtree(l)%cls(k)%snaps(j),l)),birchtree(l)%cls(k)%snaps(j))
            csnap2clus(birchtree(l)%cls(k)%snaps(j),l) = k
          end if
        else
          csnap2clus(birchtree(l)%cls(k)%snaps(j),l) = k
        end if
      end do
    end do
  end do
#ifdef ENABLE_THREADS
!$OMP END DO
#endif
!
! temporary working variables
  allocate(satisfied(tlstsz))
  satisfied(:) = 0
  allocate(snplst(tlstsz))
  snplst(:) = c_nhier + 1
  allocate(tmpix(N_NEARS,tlstsz))
  allocate(tmpdis(N_NEARS,tlstsz)) ! N_NEARS+1))
  tmpdis(:,:) = HUGE(tmpdis(1,1)) ! initialize
  tmpix(:,:) = 0                  ! initialize
  allocate(testlst(cstored,4))
  testlst(:,4) = 0
#ifdef ENABLE_THREADS
  allocate(rdbuf(rdbufsize))
  rdbufcnt = rdbufthresh
#endif
!
!
  if (tpi.le.1) call System_Clock(count=tmvars(1))
!
! allocation and initialization for temporary working variables
!
  do while (ntrees.ge.2)
!
    if (boruvkasteps.gt.0) then
!      call System_Clock(count=tmm(1))
!     sort snapshots in clusters according to tree membership (this may destroy the integrity of some aspects of the tree structure
!      -> don't use these aspects hereafter) and update the persistent sorted list itself
#ifdef ENABLE_THREADS
!$OMP DO SCHEDULE(STATIC, 1)
#endif
      do l=c_nhier+1,1,-1 
        do i=1,birchtree(l)%ncls
          if (birchtree(l)%cls(i)%nmbrs.le.1) cycle
          do j=1,birchtree(l)%cls(i)%nmbrs
            testlst(j,4) = csnap2tree(birchtree(l)%cls(i)%snaps(j))
            testlst(j,1) = j
          end do
          ixx = 1
          ixx2 = birchtree(l)%cls(i)%nmbrs
          mm = birchtree(l)%cls(i)%nmbrs
          call merge_sort(ldim=mm,up=atrue,list=testlst(1:mm,4),olist=testlst(1:mm,3),ilo=ixx,ihi=ixx2,&
 &                        idxmap=testlst(1:mm,1),olist2=testlst(1:mm,2))
!          write(*,*) l,i,':'
!          write(*,555) testlst(1:mm,3) 
          do j=1,mm
            testlst(j,4) = birchtree(l)%cls(i)%snaps(testlst(j,2))
          end do
          birchtree(l)%cls(i)%snaps(1:mm) = testlst(1:mm,4)
          birchtree(l)%cls(i)%tmpsnaps(1:mm) = testlst(1:mm,3)
!          write(*,555) birchtree(l)%cls(i)%snaps(1:mm) ! testlst(1:mm,3) 
!          write(*,555) snap2tree(birchtree(l)%cls(i)%snaps(1:mm))
        end do
      end do
#ifdef ENABLE_THREADS
!$OMP END DO
#endif
!      call System_Clock(count=tmm(2))
!      tps(1) = tps(1) + (tmm(2)-tmm(1))/thrdat%rate
    end if
!
!   first get a set of guesses based on fixed trees, clustering data structure, and prior guesses 
!   for simpler parallelizability, we do this in snapshot space
#ifdef ENABLE_THREADS
    do i=fixlims(1,tpi),fixlims(2,tpi)
      ishfx = i - fixlims(1,tpi) + 1
#else
    do i=1,cstored
      ishfx = i 
#endif
      l = snplst(ishfx) ! this is remembered across Boruvka stages
      do while (satisfied(ishfx).lt.cprogindrmax)
!      
        mycl = csnap2clus(i,l)
        tmpcld(1:calcsz) = cludata(1:calcsz,i)
!        call System_Clock(count=tmm(3)) 
        a = birchtree(l)%cls(mycl)%nmbrs
        ixx = a
!       we can quickly test whether search space is large enough to trigger random method
        ixx2 = 0
        i1 = -1
        i2 = -1
        if ((cprogindrmax-satisfied(ishfx)).gt.a) then
          thismode = 1 ! determ required
          cbnds(:) = 0
        else if (boruvkasteps.gt.0) then
!          if (csnap2tree(birchtree(l)%cls(mycl)%snaps(ixx)).lt.csnap2tree(i)) then
!            i1 = ixx
!            if (ixx.lt.a) then
!              if (csnap2tree(birchtree(l)%cls(mycl)%snaps(ixx+1)).lt.csnap2tree(i)) i1 = ixx + 1
!            end if
!          else if (csnap2tree(birchtree(l)%cls(mycl)%snaps(1)).lt.csnap2tree(i)) then
          if (csnap2tree(birchtree(l)%cls(mycl)%snaps(1)).lt.csnap2tree(i)) then
!            testlst(1:ixx,4) = csnap2tree(birchtree(l)%cls(mycl)%snaps(1:ixx))
            call binary_search(ixx,birchtree(l)%cls(mycl)%tmpsnaps(1:ixx),csnap2tree(i)-1,ixx2)
            cbnds(1) = ixx2
            i1 = ixx2 
          else
            i1 = 0
            cbnds(1) = a
          end if
!          if (csnap2tree(birchtree(l)%cls(mycl)%snaps(a-ixx+1)).gt.csnap2tree(i)) then
!            i2 = ixx
!            if ((a-ixx+1).gt.1) then
!              if (csnap2tree(birchtree(l)%cls(mycl)%snaps(a-ixx)).gt.csnap2tree(i)) i2 = ixx + 1
!            end if
!          else if (csnap2tree(birchtree(l)%cls(mycl)%snaps(a)).gt.csnap2tree(i)) then
          if (csnap2tree(birchtree(l)%cls(mycl)%snaps(ixx)).gt.csnap2tree(i)) then
!            testlst(1:ixx,4) = csnap2tree(birchtree(l)%cls(mycl)%snaps(1:ixx))
            call binary_search(ixx,birchtree(l)%cls(mycl)%tmpsnaps(1:ixx),csnap2tree(i),ixx2)
            cbnds(2) = ixx2+1
            i2 = ixx - ixx2
          else
            i2 = 0
            cbnds(2) = 1
          end if
          if ((i1+i2).gt.(cprogindrmax-satisfied(ishfx))) then
            thismode = 0
          else
            thismode = 1
          end if
        else if (boruvkasteps.eq.0) then
          cbnds(:) = 0
          if ((a-1).gt.(cprogindrmax-satisfied(ishfx))) then
            thismode = 0
          else
            thismode = 1
          end if
        end if
!        call System_Clock(count=tmm(4))
!        tps(2) = tps(2) + (tmm(4)-tmm(3))/thrdat%rate
        if (thismode.eq.1) then ! deterministic
!          call System_Clock(count=tmm(5))
          do j=1,a
            k = birchtree(l)%cls(mycl)%snaps(j)
            if (csnap2tree(k).eq.csnap2tree(i)) exit
            testcnt = testcnt + 1
            tdcnt = tdcnt + 1
            call clustering_distance(jkdist,tmpcld(1:calcsz),cludata(1:calcsz,k))
!   update shortest eligible edges per snap
            if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
              do ixx=1,N_NEARS
                if (tmpix(ixx,ishfx).eq.k) exit
                if (jkdist.lt.tmpdis(ixx,ishfx)) then
                  do ixx2=N_NEARS,ixx+1,-1 ! shift
                    tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                    tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                  end do
                  tmpix(ixx,ishfx) = k
                  tmpdis(ixx,ishfx) = jkdist 
                  exit
                end if
              end do
            end if
#ifdef ENABLE_THREADS
            if ((k.ge.fixlims(1,tpi)).AND.(k.le.fixlims(2,tpi))) then
            kshfx = k - fixlims(1,tpi) + 1
#else
            kshfx = k
#endif
            if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
              do ixx=1,N_NEARS
                if (tmpix(ixx,kshfx).eq.i) exit
                if (jkdist.lt.tmpdis(ixx,kshfx)) then
                  do ixx2=N_NEARS,ixx+1,-1 ! shift
                    tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                    tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                  end do
                  tmpix(ixx,kshfx) = i
                  tmpdis(ixx,kshfx) = jkdist 
                  exit
                end if
              end do
            end if
#ifdef ENABLE_THREADS
            end if
#endif
            satisfied(ishfx) = satisfied(ishfx) + 1
          end do
          do j=a,1,-1
            k = birchtree(l)%cls(mycl)%snaps(j)
            if (csnap2tree(k).eq.csnap2tree(i)) exit
            tdcnt = tdcnt + 1
            testcnt = testcnt + 1
            call clustering_distance(jkdist,tmpcld(1:calcsz),cludata(1:calcsz,k))
!   update shortest eligible edges per snap
            if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
              do ixx=1,N_NEARS
                if (tmpix(ixx,ishfx).eq.k) exit
                if (jkdist.lt.tmpdis(ixx,ishfx)) then
                  do ixx2=N_NEARS,ixx+1,-1 ! shift
                    tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                    tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                  end do
                  tmpix(ixx,ishfx) = k
                  tmpdis(ixx,ishfx) = jkdist 
                  exit
                end if
              end do
            end if
#ifdef ENABLE_THREADS
            if ((k.ge.fixlims(1,tpi)).AND.(k.le.fixlims(2,tpi))) then
            kshfx = k - fixlims(1,tpi) + 1
#else
            kshfx = k
#endif
            if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
              do ixx=1,N_NEARS
                if (tmpix(ixx,kshfx).eq.i) exit
                if (jkdist.lt.tmpdis(ixx,kshfx)) then
                  do ixx2=N_NEARS,ixx+1,-1 ! shift
                    tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                    tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                  end do
                  tmpix(ixx,kshfx) = i
                  tmpdis(ixx,kshfx) = jkdist 
                  exit
                end if
              end do
            end if
#ifdef ENABLE_THREADS
            end if
#endif
            satisfied(ishfx) = satisfied(ishfx) + 1
          end do
!          call System_Clock(count=tmm(6))
!          tps(3) = tps(3) + (tmm(6)-tmm(5))/thrdat%rate
        else ! random
!          call System_Clock(count=tmm(7))
          m = a
          b1 = 0 ! max(1,a/2)
          b2 = 0 ! b1
          if (cbnds(1).gt.0) then
            if (cbnds(2).eq.1) then
              b2 = 1
              m = cbnds(1)
            else if (cbnds(1).gt.cbnds(2)) then
              b2 = cbnds(2)
              m = cbnds(1) - cbnds(2) + 1
            else if (cbnds(2).gt.cbnds(1)) then
              b2 = cbnds(2)
              m = a - (cbnds(2) - cbnds(1) - 1)
            end if
          end if
#ifdef ENABLE_THREADS
          if (((rdbufcnt+cprogindrmax+1).gt.rdbufthresh).AND.(rdbufcnt.le.rdbufthresh)) then
            kix = rdbufthresh+1
            if (kix.gt.rdbufsize) kix = 1
            kk = 2*cprogindrmax
!$OMP CRITICAL(RDBUF_REPL)
            call get_nrandoms(kk,rdbuf(kix:(kix+2*cprogindrmax-1)))
!$OMP END CRITICAL(RDBUF_REPL)
            rdbufthresh = rdbufthresh + 2*cprogindrmax
            if (rdbufthresh.gt.rdbufsize) rdbufthresh = 2*cprogindrmax
          end if
#endif
          if (boruvkasteps.gt.0) then
            do while (satisfied(ishfx).lt.cprogindrmax)
#ifdef ENABLE_THREADS
              rdbufcnt = rdbufcnt + 1
              if (rdbufcnt.gt.rdbufsize) rdbufcnt = 1 
              kk = int(rdbuf(rdbufcnt)*m) + b2
#else
              kk = int(random()*m) + b2
#endif
              do ki=1,cprogbatchsz
                kix = kk
                if (kix.gt.a) kix = kix - a
                if (kix.le.0) kix = kix + a
                k = birchtree(l)%cls(mycl)%snaps(kix)
                if (csnap2tree(i).eq.csnap2tree(k)) then
                  call fexit()
                end if
                trcnt = trcnt + 1
                testcnt = testcnt + 1
                call clustering_distance(jkdist,tmpcld(1:calcsz),cludata(1:calcsz,k))
!                update shortest eligible edges per snap
                if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,ishfx).eq.k) exit
                    if (jkdist.lt.tmpdis(ixx,ishfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                        tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                      end do
                      tmpix(ixx,ishfx) = k
                      tmpdis(ixx,ishfx) = jkdist 
                      exit
                    end if
                  end do
                end if
#ifdef ENABLE_THREADS
                if ((k.ge.fixlims(1,tpi)).AND.(k.le.fixlims(2,tpi))) then
                kshfx = k - fixlims(1,tpi) + 1
#else
                kshfx = k
#endif
                if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,kshfx).eq.i) exit
                    if (jkdist.lt.tmpdis(ixx,kshfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                        tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                      end do
                      tmpix(ixx,kshfx) = i
                      tmpdis(ixx,kshfx) = jkdist 
                      exit
                    end if
                  end do
                end if
#ifdef ENABLE_THREADS
                end if
#endif
                kk = kk + 1
                if (kk.gt.(m+b2-1)) kk = b2
                satisfied(ishfx) = satisfied(ishfx) + 1
                if (satisfied(ishfx).eq.cprogindrmax) exit
              end do
            end do
          else
            kixi = max(1,a/cprogbatchsz)
            do while (satisfied(ishfx).lt.cprogindrmax)
#ifdef ENABLE_THREADS
              rdbufcnt = rdbufcnt + 1
              if (rdbufcnt.gt.rdbufsize) rdbufcnt = 1
              kk = int(rdbuf(rdbufcnt)*a)
#else
              kk = int(random()*a)
#endif
              do ki=1,cprogbatchsz
                kix = kk
                if (kix.gt.a) kix = kix - a
                if (kix.le.0) kix = kix + a
                k = birchtree(l)%cls(mycl)%snaps(kix)
                if (csnap2tree(i).eq.csnap2tree(k)) then
                  cycle ! in first stage, only a single snapshot per cluster needs to be cycled
                end if
                call clustering_distance(jkdist,tmpcld(1:calcsz),cludata(1:calcsz,k))
                trcnt = trcnt + 1
                testcnt = testcnt + 1
!   update shortest eligible edges per snap
                if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,ishfx).eq.k) exit
                    if (jkdist.lt.tmpdis(ixx,ishfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                        tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                      end do
                      tmpix(ixx,ishfx) = k
                      tmpdis(ixx,ishfx) = jkdist
                      exit
                    end if
                  end do
                end if
#ifdef ENABLE_THREADS
                if ((k.ge.fixlims(1,tpi)).AND.(k.le.fixlims(2,tpi))) then
                kshfx = k - fixlims(1,tpi) + 1
#else
                kshfx = k
#endif
                if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,kshfx).eq.i) exit
                    if (jkdist.lt.tmpdis(ixx,kshfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                        tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                      end do
                      tmpix(ixx,kshfx) = i
                      tmpdis(ixx,kshfx) = jkdist
                      exit
                    end if
                  end do
                end if
#ifdef ENABLE_THREADS
                end if
#endif
                kk = kk + kixi
                if (kk.gt.a) kk = kk - a
                satisfied(ishfx) = satisfied(ishfx) + 1
                if (satisfied(ishfx).eq.cprogindrmax) exit
              end do
            end do
          end if
!          call System_Clock(count=tmm(8))
!          tps(4) = tps(4) + (tmm(8)-tmm(7))/thrdat%rate
        end if
        l = l - 1
        if (satisfied(ishfx).eq.0) snplst(ishfx) = l
        if ((snplst(ishfx)-(l+1)).eq.cprogrdepth) exit
        if (l.eq.0) exit
      end do
    end do
!
!    call System_Clock(count=tmm(9))
    testlst(1:ntrees,3) = 0
!
#ifdef ENABLE_THREADS
    do cursnp=fixlims(1,tpi),fixlims(2,tpi)
#else
    do cursnp=1,cstored
#endif
#ifdef ENABLE_THREADS
      ishfx = cursnp - fixlims(1,tpi) + 1
#else
      ishfx = cursnp
#endif
      j = csnap2tree(cursnp)
      if (testlst(j,3).eq.0) then ! first edge for that tree
        testlst(j,1) = cursnp
        testlst(j,2) = tmpix(1,ishfx)
        testlst(j,3) = ishfx
        testlst(j,4) = 1
      else if (tmpdis(1,ishfx).lt.tmpdis(testlst(j,4),testlst(j,3))) then
        testlst(j,1) = cursnp
        testlst(j,2) = tmpix(1,ishfx)
        testlst(j,3) = ishfx
        testlst(j,4) = 1
      end if
      do ixx=1,N_NEARS
        if (tmpix(ixx,ishfx).le.0) cycle
        k = csnap2tree(tmpix(ixx,ishfx))
        if (testlst(k,3).eq.0) then ! first edge for that tree
          testlst(k,1) = tmpix(ixx,ishfx)
          testlst(k,2) = cursnp
          testlst(k,3) = ishfx
          testlst(k,4) = ixx
        else if (tmpdis(ixx,ishfx).lt.tmpdis(testlst(k,4),testlst(k,3))) then
          testlst(k,1) = tmpix(ixx,ishfx)
          testlst(k,2) = cursnp
          testlst(k,3) = ishfx
          testlst(k,4) = ixx
        end if
      end do
    end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(COLLECT_EDGES)
#endif
    do k=1,ntrees ! could be made faster for first stage at least
      if (testlst(k,3).gt.0) then
        if (tmpdis(testlst(k,4),testlst(k,3)).lt.tmptree(k)%mind) then
          tmptree(k)%mind = tmpdis(testlst(k,4),testlst(k,3))
          tmptree(k)%mine(1) = testlst(k,1)
          tmptree(k)%mine(2) = testlst(k,2)
        end if
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(COLLECT_EDGES)
#endif
!    call System_Clock(count=tmm(10))
!    tps(5) = tps(5) + (tmm(10)-tmm(9))/thrdat%rate
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!    call System_Clock(count=tmm(11))
!$OMP MASTER
#endif
! merge trees (& update approximate minimum spanning tree):
    tmptree(1:ntrees)%nsiblings = 0
    do k=1,ntrees
      tmptree(k)%ptr = k
    end do
!
    do e=1,ntrees ! nedges
      i1 = tmptree(e)%ptr                               ! tree index of first edge endpoint
      i2 = tmptree(csnap2tree(tmptree(e)%mine(2)))%ptr  ! tree index of second edge endpoint
      if (i1.ne.i2) then ! edge does not introduce cycle
! add corresponding edge to approximate minimum spanning tree:
        nmstedges = nmstedges+1
        if (nmstedges.ge.cstored) then
          write(ilog,*) 'Fatal. The number of edges in the approximate minimum spanning tree is not correct. Please report &
 &this bug.'
          call fexit()
        end if
        mstedges(1,nmstedges) = tmptree(e)%mine(1) ! edgelst(e,1)
        mstedges(2,nmstedges) = tmptree(e)%mine(2) ! edgelst(e,2)
        lmstedges(nmstedges) = tmptree(e)%mind     ! ledgelst(e)
! update pointer for tree and all its siblings obtained through prior merge operations
        tmptree(i1)%ptr = i2
        tmptree(i2)%nsiblings = tmptree(i2)%nsiblings + 1
        if (tmptree(i2)%nsiblings.gt.tmptree(i2)%nsibalsz) call pidxtree_growsiblings(tmptree(i2))
        tmptree(i2)%siblings(tmptree(i2)%nsiblings) = i1
        do j=1,tmptree(i1)%nsiblings
          mm = tmptree(i1)%siblings(j)
          tmptree(i2)%nsiblings = tmptree(i2)%nsiblings + 1
          if (tmptree(i2)%nsiblings.gt.tmptree(i2)%nsibalsz) call pidxtree_growsiblings(tmptree(i2))
          tmptree(i2)%siblings(tmptree(i2)%nsiblings) = mm
          tmptree(mm)%ptr = i2
        end do
      end if
    end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
! reassign csnap2tree based on pointer
#ifdef ENABLE_THREADS
    do j=fixlims(1,tpi),fixlims(2,tpi)
#else
    do j=1,cstored
#endif
      csnap2tree(j) = tmptree(csnap2tree(j))%ptr
    end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
! now construct a new set of sorted from csnap2tree and readjust csnap2tree (this does not scale)
    mm = 0
    oldntrees = ntrees
    ntrees = 0
    do j=1,cstored
      k = csnap2tree(j)
      if (testlst(k,3).ge.0) then
        ntrees = ntrees + 1
        testlst(k,3) = -ntrees
      end if
    end do
!  
#ifdef ENABLE_THREADS
!$OMP BARRIER
    do j=fixlims(1,tpi),fixlims(2,tpi)
#else
    do j=1,cstored
#endif
      csnap2tree(j) = -testlst(csnap2tree(j),3)
    end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    do j=ntrees+1,oldntrees
      if (allocated(tmptree(j)%siblings).EQV..true.) deallocate(tmptree(j)%siblings)
    end do
    tmptree(1:ntrees)%mind = HUGE(jkdist)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
    boruvkasteps = boruvkasteps+1
! manage stored guesses
#ifdef ENABLE_THREADS
    do j=fixlims(1,tpi),fixlims(2,tpi)
      ishfx = j - fixlims(1,tpi) + 1
#else
    do j=1,cstored
      ishfx = j
#endif
      satisfied(ishfx) = 0
      ixx2 = 0
      do ixx=1,N_NEARS
        if (tmpix(ixx,ishfx).gt.0) then
          if (csnap2tree(j).ne.csnap2tree(tmpix(ixx,ishfx))) then
            ixx2 = ixx2 + 1
            tmpix(ixx2,ishfx) = tmpix(ixx,ishfx)
            tmpdis(ixx2,ishfx) = tmpdis(ixx,ishfx)
          end if
        end if
      end do
      do ixx=max(1,ixx2),N_NEARS
        tmpix(ixx,ishfx) = 0
        tmpdis(ixx,ishfx) = HUGE(tmpdis(ixx,ishfx))
      end do
    end do
!    call System_Clock(count=tmm(12))
!    tps(6) = tps(6) + (tmm(12)-tmm(11))/thrdat%rate
!
 567 format('... time for Boruvka stage ',i4,': ',g11.3,'s ...')
 456 format(i4,8(g12.4,1x),i10,i10,i6,i6,1x,g12.4,g12.4,g12.4)
#ifdef ENABLE_THREADS
    if (tpi.le.1) then
      call System_Clock(count=tmvars(2))
      write(*,567) boruvkasteps,(tmvars(2)-tmvars(1))/thrdat%rate
      tmvars(1) = tmvars(2)
      flush(ilog)
    end if
!    write(*,456) tpi,tps(3:5),tps(3)+tps(4),sum(tps(3:5)),1.0e6*tps(3)/tdcnt,1.0e6*tps(4)/trcnt,&
! &               (tps(4)/tps(3))/(1.0*trcnt/(1.0*tdcnt)),tdcnt,trcnt
!    tps(:) = 0.0
!    trcnt = 0
!    tdcnt = 0
!$OMP BARRIER
#else
    call System_Clock(count=tmvars(2))
    write(*,567) boruvkasteps,(tmvars(2)-tmvars(1))/thrdat%rate
    tmvars(1) = tmvars(2)
#endif
!
  end do
!
#ifdef ENABLE_THREADS
!  deallocate(rdlst)
!$OMP CRITICAL(TOTCNTER)
  cdevalcnt = cdevalcnt + testcnt
!$OMP END CRITICAL(TOTCNTER)
#else
  cdevalcnt = testcnt
#endif
!
  deallocate(testlst)
  deallocate(snplst)
  deallocate(satisfied)
  deallocate(tmpdis)
  deallocate(tmpix)
#ifdef ENABLE_THREADS
  deallocate(rdbuf)
#endif
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  if (cstored-1.ne.nmstedges) then
    write(ilog,*) 'Fatal. The number of edges in the approximate minimum spanning tree is not correct. Please report &
 &this bug.'
    call fexit()
  end if
!
  allocate(approxmst(cstored))
  approxmst(1:cstored)%deg = 0
 587 format(' Weight of Short Spanning Tree: ',1x,g12.5,a)
  write(ilog,*)
  if (cdis_crit.le.2) then
    write(ilog,587) sum(lmstedges(1:cstored-1)),' degrees'
  else if (cdis_crit.le.4) then
    write(ilog,587) sum(lmstedges(1:cstored-1)),' '
  else if (cdis_crit.le.10) then
    write(ilog,587) sum(lmstedges(1:cstored-1)),' Angstrom'
  end if
  write(ilog,*)
  call gen_MST(mstedges,lmstedges,approxmst)
  do i=1,ntrees
    if (allocated(tmptree(i)%siblings).EQV..true.) deallocate(tmptree(i)%siblings)
  end do
  deallocate(tmptree)
  deallocate(csnap2tree)
  deallocate(csnap2clus)
  deallocate(mstedges)
  deallocate(lmstedges)
!
 77 format(a,20(i18,1x))
  write(ilog,*) '... done after ',cdevalcnt,' additional distance evaluations.'
  write(ilog,*)
#ifdef ENABLE_THREADS
!$OMP END MASTER
  deallocate(fixlims)
!$OMP BARRIER
#endif
!
end
!
!--------------------------------------------------------------------------------
!
! this routines transcribes a list of edges with their distances, that effectively 
! describe the MST, into an array of adjacency list objects
!
subroutine gen_MST(mstedges,lmstedges,mst)
!
  use clusters
!
  implicit none
!
  integer e,v
  type(t_adjlist) mst(cstored)
  integer mstedges(2,cstored-1)
  REAL(kind=4) lmstedges(cstored-1)
!
! transform edgelist of MST (mstedges) to adjacencylist:
  do e=1,cstored
    allocate(mst(e)%adj(1))
    allocate(mst(e)%dist(1))
  end do
  do e=1,cstored-1
    do v=1,2
      if (mst(mstedges(v,e))%deg.gt.0) call extend_adjlst_byone(mst(mstedges(v,e)))
      mst(mstedges(v,e))%deg = mst(mstedges(v,e))%deg + 1
      mst(mstedges(v,e))%adj(mst(mstedges(v,e))%deg) = mstedges(3-v,e)
      mst(mstedges(v,e))%dist(mst(mstedges(v,e))%deg) = lmstedges(e)
    end do
  end do
!
end
!
!------------------------------------------------------------------------------------------------------
!
! a simple helper to grow as conservatively (and slowly) as possible an adjacency list object
!
subroutine extend_adjlst_byone(mstnode)
!
  use clusters
!
  implicit none
!
  integer, ALLOCATABLE:: itmp1(:)
  real(KIND=4), ALLOCATABLE:: rtmp1(:)
  type(t_adjlist) mstnode
!
  allocate(itmp1(mstnode%deg))
  allocate(rtmp1(mstnode%deg))
!
  itmp1(1:mstnode%deg) = mstnode%adj(1:mstnode%deg)
  rtmp1(1:mstnode%deg) = mstnode%dist(1:mstnode%deg)
  deallocate(mstnode%adj)
  deallocate(mstnode%dist)
!
  allocate(mstnode%adj(mstnode%deg+1))
  allocate(mstnode%dist(mstnode%deg+1))
  mstnode%adj(1:mstnode%deg) = itmp1(1:mstnode%deg)
  mstnode%dist(1:mstnode%deg) = rtmp1(1:mstnode%deg)
!
  deallocate(itmp1)
  deallocate(rtmp1)
!
end
!
!-----------------------------------------------------------------------------------------------------
!
! a subroutine to contract terminal vertices into their parent by artificially setting the distance 
! to a negative value
!
subroutine contract_mst(alst,nrnds)
!
  use clusters
  use iounit, ONLY: ilog
!
  implicit none
!
  integer, INTENT(IN):: nrnds
!
  integer i,j,kk,thej,ll,mm
  type(t_adjlist) alst(cstored)
!
  logical, ALLOCATABLE:: terminal(:)
!
  allocate(terminal(cstored))
  terminal(:) = .false.
!
  mm = 0
  ll = 0
  do i=1,cstored
    if (alst(i)%deg.eq.1) then
      terminal(i) = .true.
      ll = ll + 1
    end if
  end do
!
  mm = 0
  do kk=1,nrnds
    do i=1,cstored
      if (terminal(i).EQV..true.) then
        do j=1,alst(i)%deg
          if (alst(i)%dist(j).ge.0.0) then
            thej = alst(i)%adj(j)
            if (alst(i)%dist(j).eq.0.0) then
              alst(i)%dist(j) = -10.0*TINY(alst(i)%dist(j))
            else
              alst(i)%dist(j) = -1.0/alst(i)%dist(j) 
            end if
            mm = mm + 1
            exit
          end if
        end do
        do ll=1,alst(thej)%deg
          if (alst(thej)%adj(ll).eq.i) then
            alst(thej)%dist(ll) = alst(i)%dist(j)
            exit
          end if
        end do
      end if
    end do
!   reset terminal status
    terminal(:) = .false.
    ll = 0
    do i=1,cstored
      thej = 0
      do j=1,alst(i)%deg
        if (alst(i)%dist(j).ge.0.0) thej = thej + 1
      end do
      if (thej.eq.1) then
        terminal(i) = .true.
        ll = ll + 1
      end if
    end do
    if (ll.eq.0) exit
  end do
!
  write(ilog,*) 
  write(ilog,104) mm,100.0*(mm)/(1.0*cstored)
  write(ilog,105) ll,100.0*(ll)/(1.0*cstored)
 104 format('The number of promoted (folded) edges in the ST is ',i10,' (',f8.3,' %).')
 105 format('The number of remaining leaf vertices in the ST is ',i10,' (',f8.3,' %).')
  deallocate(terminal)
!
end
!
!
!------------------------------------------------------------------------------------------------------
!
! using PRIM's algorithm, this routine will trace the MST adjacency list to derive the progress index per se
! (along with the minimal distance to the current set (distv), the inverse map (invvec) and the parent/source vector (iv2))
!
subroutine gen_progind_from_adjlst(alst,starter,progind,distv,invvec,iv2)
!
  use clusters
  use iounit
!
  implicit none
!
  type(t_adjlist) alst(cstored)
  logical, ALLOCATABLE:: added(:), inprogind(:)
  integer, ALLOCATABLE:: heap(:),hsource(:)
  integer heapsize,lprogind,invvec(cstored+2),progind(cstored),iv2(cstored),j,starter
  real(KIND=4), ALLOCATABLE:: key(:)
  real(KIND=4) distv(cstored)
!
  allocate(inprogind(cstored))
  allocate(added(cstored))
  allocate(key(cstored))
  allocate(heap(cstored))
  allocate(hsource(cstored))
  added(:) = .false.
  inprogind(:) = .false.
  distv(:) = 0.0
  key(:) = 0.0
  hsource(:) = 0
!
! add first snapshot to progress index:
  if (starter.le.cstored) then
    progind(1) = starter
  else
    progind(1) = 1
    write(ilog,*) 'Warning. The snapshot index requested is not available (there are ',cstored,' snapshots in memory). &
 &Using first one instead.'
  end if
  lprogind = 1
  distv(lprogind) = 0.0
  added(progind(lprogind)) = .true.
  inprogind(progind(lprogind)) = .true.
  invvec(progind(lprogind)+1) = 1
  iv2(lprogind) = progind(lprogind)
! build heap:
  heapsize = alst(progind(lprogind))%deg
  heap(1:heapsize) = alst(progind(lprogind))%adj(:)
  key(1:heapsize) = alst(progind(lprogind))%dist(:)
  hsource(1:heapsize) = progind(lprogind)
  call hbuild(heap(1:heapsize),heapsize,key(1:heapsize),hsource(1:heapsize))
  do j=1,alst(progind(lprogind))%deg
    added(alst(progind(lprogind))%adj(j)) = .true.
  end do
  do while (lprogind.lt.cstored)
    do j=1,alst(progind(lprogind))%deg
! add neighbors of last snapshot in progress index to heap
      if (added(alst(progind(lprogind))%adj(j)).EQV..false.) then
        call hinsert(heap(1:(heapsize+1)),heapsize,key(1:(heapsize+1)),hsource(1:(heapsize+1)),&
 &                   alst(progind(lprogind))%adj(j),alst(progind(lprogind))%dist(j),progind(lprogind))
        added(alst(progind(lprogind))%adj(j)) = .true.
      end if
    end do
! append next snapshot to progind():
    lprogind = lprogind+1
    progind(lprogind) = heap(1)
    iv2(lprogind) = hsource(1)
    distv(lprogind) = key(1)
    invvec(heap(1)+1) = lprogind
    inprogind(progind(lprogind)) = .true.
! remove added snapshot from heap:
    call hremovemin(heap(1:heapsize),heapsize,key(1:heapsize),hsource(1:heapsize))
  end do
!
  do j=1,cstored
    if (distv(j).lt.0.0) then
      if (distv(j).ge.-10.0*TINY(distv(j))) then
        distv(j) = 0.0
      else
        distv(j) = -1.0/distv(j)
      end if
    end if
  end do
!
  deallocate(key)
  deallocate(added)
  deallocate(inprogind)
  deallocate(heap)
  deallocate(hsource)
!
end
!
!------------------------------------------------------------------------------------------------------
!
subroutine do_prog_index()
!
  use clusters
  use iounit
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer refi,i,kk,aone,idxc,ii,jj
  RTYPE random
  integer, ALLOCATABLE:: progind(:),invvec(:),cutv(:),iv2(:),proflst(:)
  real(KIND=4), ALLOCATABLE:: distv(:)
  logical candid,profile_min,exists
  integer tl2,iu,freeunit
  character(12) nod2
  character(MAXSTRLEN) fn
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
!
  write(ilog,*)
  write(ilog,*) 'Now deriving progress index and computing MFPT-related quantities ...'
!
  allocate(distv(cstored))
  allocate(progind(cstored))
  allocate(invvec(cstored+2))
  allocate(iv2(cstored))
  allocate(proflst(max(1,cstored/csivmin)))
  if (csivmin.ge.cstored) then
    write(ilog,*) 'Warning. Setting for FMCSC_CBASINMAX is inappropriate (too large). Adjusting to '&
 &,max(cstored/10,1),'.'
    csivmin = max(1,cstored/10)
  end if
  if ((cprogindstart.eq.-1).AND.(cprogindex.eq.2)) then
    write(ilog,*) 'Using representative snapshot of largest cluster (#',birchtree(c_nhier+1)%cls(1)%center,&
 &') as reference for progress index (setting -1).'
    cprogindstart = birchtree(c_nhier+1)%cls(1)%center ! center of largest cluster in hierarchical tree
  end if
  if (cprogindstart.gt.cstored) then
    write(ilog,*) 'Warning. Setting for FMCSC_CPROGINDSTART requests a snapshot that does not exist. Remember &
 &that numbering is with respect to the snapshots currently in memory. Using first snapshot instead.'
    cprogindstart = 1
  end if
  if (cprogpwidth.ge.cstored/2) then
    write(ilog,*) 'Warning. Setting for FMCSC_CPROGINDWIDTH is inappropriate (too large). Adjusting to '&
 &,max(cstored/10,1),'.'
    cprogpwidth = max(1,cstored/10)
  end if
!
  aone = 1
!
  if (cprogfold.gt.0) then
    call contract_mst(approxmst,cprogfold)
  end if
!
  if (cprogindstart.eq.0) then
    allocate(cutv(cstored))
    invvec(1) = 3*cstored
    invvec(cstored+2) = 6*cstored
    refi = ceiling(random()*cstored)
!
    call gen_progind_from_adjlst(approxmst,refi,progind,distv,invvec,iv2)
!
!   generate N(A->B) + N(B->A) with boundary condition allB
    cutv(1) = 2 
    do i=2,cstored
      kk = progind(i) + 1
      if ((invvec(kk+1).lt.i).AND.(invvec(kk-1).lt.i)) then
        cutv(i) = cutv(i-1) - 2
      else if ((invvec(kk+1).gt.i).AND.(invvec(kk-1).gt.i)) then
        cutv(i) = cutv(i-1) + 2
      else
        cutv(i) = cutv(i-1)
      end if
    end do
    idxc = 0
    do i=csivmin+1,cstored-csivmin
      candid = profile_min(cutv(1:cstored),cstored,i,csivmin,csivmax,aone,kk)
      if (candid.EQV..true.) then
        idxc = idxc + 1
        proflst(idxc) = i
      end if
    end do
    if (idxc.le.0) then
      write(ilog,*) 'Warning. Automatic identification of starting snapshots for profile &
 &generation failed. Creating a single profile from snapshot #1.'
      idxc = 1
      proflst(idxc) = 1
    end if
  else if (cprogindstart.gt.0) then
    idxc = 1
    proflst(idxc) = cprogindstart
  end if
!
  do i=1,idxc
    invvec(1) = 3*cstored
    invvec(cstored+2) = 6*cstored
!
    call gen_progind_from_adjlst(approxmst,proflst(i),progind,distv,invvec,iv2)
    tl2 = 12
    call int2str(proflst(i),nod2,tl2)
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,nod,re_aux(10))
      fn = 'N_'//nod(1:re_aux(10))//'_PROGIDX_'//nod2(1:tl2)//'.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'PROGIDX_'//nod2(1:tl2)//'.dat'
    end if
#else
    fn = 'PROGIDX_'//nod2(1:tl2)//'.dat'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old',position='append')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
    call gen_manycuts(progind,distv,invvec,iv2,iu,cprogpwidth)
    close(unit=iu)
  end do
!
  deallocate(proflst)
  deallocate(iv2)
  deallocate(invvec)
  deallocate(distv)
  deallocate(progind)
!
  write(ilog,*) '... done.'
  write(ilog,*)
!
end
!
!----------------------------------------------------------------------------------------------------
!
subroutine gen_manycuts(setis,distv,invvec,ivec2,iu,pwidth)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer tmat(3,3),vv1(3),vv2(3),vv3(3)
  integer iu,i,j,k,pwidth,ii,jj,kk,ik,iik,jk
  integer invvec(cstored+2),ivec2(cstored),setis(cstored)
  real(KIND=4) distv(cstored)
  integer, ALLOCATABLE:: cutv2(:),statev(:)
  logical atrue,docyc
!
  type t_lkmap
    integer, ALLOCATABLE:: lkix(:)
    integer nlk,alsz
  end type t_lkmap
!
  type(t_lkmap), ALLOCATABLE:: lkmap(:)
!
  atrue = .true.
!
! create a data structure that allows immediate access to existing links (sconnect) from snapshot index
  allocate(lkmap(cstored))
  lkmap(:)%nlk = 0
  lkmap(:)%alsz = 2
  do i=1,cstored
    allocate(lkmap(i)%lkix(lkmap(i)%alsz))
  end do
!
  do i=1,sconnectsz
    if ((sconnect(i,1).le.0).OR.(sconnect(i,2).le.0)) cycle
    if ((sconnect(i,1).gt.cstored).OR.(sconnect(i,2).gt.cstored)) then
      write(ilog,*) 'Warning. Skipping progress index analysis due to mismatch between snapshot connectivity matrix &
 &and data size. This is most likely caused by a truncated trajectory analysis run.'
      return
    end if
    ii = sconnect(i,1)
    jj = sconnect(i,2)
    lkmap(ii)%nlk = lkmap(ii)%nlk + 1
    if (lkmap(ii)%nlk.gt.lkmap(ii)%alsz) then
      kk = lkmap(ii)%alsz
      lkmap(ii)%alsz = lkmap(ii)%alsz + 10
      call resize_vec(kk,lkmap(ii)%lkix,lkmap(ii)%alsz,atrue)
    end if
    lkmap(ii)%lkix(lkmap(ii)%nlk) = i
    lkmap(jj)%nlk = lkmap(jj)%nlk + 1
    if (lkmap(jj)%nlk.gt.lkmap(jj)%alsz) then
      kk = lkmap(jj)%alsz
      lkmap(jj)%alsz = lkmap(jj)%alsz + 10
      call resize_vec(kk,lkmap(jj)%lkix,lkmap(jj)%alsz,atrue)
    end if
    lkmap(jj)%lkix(lkmap(jj)%nlk) = i
  end do
!
  allocate(cutv2(cstored))
  cutv2(:) = 0
!
! global cut for 2-state model (c(AB)+c(BA)) is stored in cutv
  do i=1,cstored
    if (i.gt.1) cutv2(i) = cutv2(i-1)
    if ((setis(i).le.clagt_msm).OR.(setis(i).gt.(cstored-clagt_msm))) then  ! boundary condition is BX...XB 
      cutv2(i) = cutv2(i) + 1
    end if
    do j=1,lkmap(setis(i))%nlk
      jj = lkmap(setis(i))%lkix(j)
      kk = sconnect(jj,2) + 1
      if (sconnect(jj,2).eq.setis(i)) kk = sconnect(jj,1) + 1
      if (invvec(kk).lt.i) then
        cutv2(i) = cutv2(i) - 1
      else
        cutv2(i) = cutv2(i) + 1
      end if
    end do 
  end do
!
! in the 3-state model, the boundary condition is all-C
  allocate(statev(cstored))
  statev(:) = 3
  tmat(:,:) = 0
  tmat(3,3) = sum(lkmap(1:cstored)%nlk)/2 + 2*clagt_msm ! all map-based links + transitions with lag time into BC (state C)
! first populate tmat to state first pwidth snapshots in PIX order being in state B
  do i=1,pwidth
    ii = setis(i)
    if (ii.le.clagt_msm) then
      tmat(3,statev(ii)) = tmat(3,statev(ii)) - 1
    end if
    if (ii.gt.(cstored-clagt_msm)) then  ! boundary condition is CX...XC
      tmat(statev(ii),3) = tmat(statev(ii),3) - 1
    end if
    do j=1,lkmap(ii)%nlk
      jj = sconnect(lkmap(ii)%lkix(j),1)
      kk = sconnect(lkmap(ii)%lkix(j),2)
      tmat(statev(jj),statev(kk)) = tmat(statev(jj),statev(kk)) - 1
    end do
    statev(ii) = 2
    if (ii.le.clagt_msm) then
      tmat(3,statev(ii)) = tmat(3,statev(ii)) + 1
    end if
    if (ii.gt.(cstored-clagt_msm)) then  ! boundary condition is CX...XC
      tmat(statev(ii),3) = tmat(statev(ii),3) + 1
    end if
    do j=1,lkmap(ii)%nlk
      jj = sconnect(lkmap(ii)%lkix(j),1)
      kk = sconnect(lkmap(ii)%lkix(j),2)
      tmat(statev(jj),statev(kk)) = tmat(statev(jj),statev(kk)) + 1
    end do
  end do
! at each subsequent step, one snap becomes A, one becomes B, and one becomes C except at the boundaries of (1:cstored)
  do i=1,cstored
    do k=-1,1
      if ((i+k*pwidth.le.0).OR.((i+k*pwidth).gt.cstored)) cycle
      ii = setis(i+k*pwidth)
      if (ii.le.clagt_msm) then
        tmat(3,statev(ii)) = tmat(3,statev(ii)) - 1
      end if
      if (ii.gt.(cstored-clagt_msm)) then  ! boundary condition is CX...XC
        tmat(statev(ii),3) = tmat(statev(ii),3) - 1
      end if
      do j=1,lkmap(ii)%nlk
!       links will be duplicated if 2 or 3 of the actual snaps at the 3 positions are indeed linked in sconnect -> screen
        docyc = .false.
        do ik=-1,k-1
          if ((i+ik*pwidth.le.0).OR.((i+ik*pwidth).gt.cstored)) cycle
          iik = setis(i+ik*pwidth)
          do jk=1,lkmap(iik)%nlk
            if (lkmap(iik)%lkix(jk).eq.lkmap(ii)%lkix(j)) docyc = .true.
          end do
        end do
        if (docyc.EQV..true.) cycle
        jj = sconnect(lkmap(ii)%lkix(j),1)
        kk = sconnect(lkmap(ii)%lkix(j),2)
        tmat(statev(jj),statev(kk)) = tmat(statev(jj),statev(kk)) - 1
      end do
    end do
    if ((cstored-i).ge.pwidth) statev(setis(i+pwidth)) = 2 ! from 3
    if (i.gt.pwidth) statev(setis(i-pwidth)) = 3 ! from 1
    statev(setis(i)) = 1 ! from 2
    do k=-1,1
      if ((i+k*pwidth.le.0).OR.((i+k*pwidth).gt.cstored)) cycle
      ii = setis(i+k*pwidth)
      if (ii.le.clagt_msm) then
        tmat(3,statev(ii)) = tmat(3,statev(ii)) + 1
      end if
      if (ii.gt.(cstored-clagt_msm)) then  ! boundary condition is CX...XC
        tmat(statev(ii),3) = tmat(statev(ii),3) + 1
      end if
      do j=1,lkmap(ii)%nlk
!       links will be duplicated if 2 or 3 of the actual snaps at the 3 positions are indeed linked in sconnect -> screen
        docyc = .false.
        do ik=-1,k-1
          if ((i+ik*pwidth.le.0).OR.((i+ik*pwidth).gt.cstored)) cycle
          iik = setis(i+ik*pwidth)
          do jk=1,lkmap(iik)%nlk
            if (lkmap(iik)%lkix(jk).eq.lkmap(ii)%lkix(j)) docyc = .true.
          end do
        end do
        if (docyc.EQV..true.) cycle
        jj = sconnect(lkmap(ii)%lkix(j),1)
        kk = sconnect(lkmap(ii)%lkix(j),2)
        tmat(statev(jj),statev(kk)) = tmat(statev(jj),statev(kk)) + 1
      end do
    end do
    ii = min(i,pwidth)
    jj = min(cstored-i,pwidth)
    kk = max(0,cstored - ii - jj)
    vv1(:) = tmat(1,:)
    vv2(:) = tmat(2,:)
    vv3(:) = tmat(3,:)
    write(iu,666) i,cstored-i,setis(i),cutv2(i),distv(i),ivec2(i),invvec(ivec2(i)+1),cutv2(invvec(ivec2(i)+1)),vv1,vv2,vv3,ii,jj,kk
  end do
 666 format(4(i10,1x),g12.5,1x,1000(i10,1x))
!
  deallocate(statev)
  deallocate(cutv2)
  do i=1,cstored
    if (allocated(lkmap(i)%lkix).EQV..true.) deallocate(lkmap(i)%lkix)
  end do
  deallocate(lkmap)
!
end
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
