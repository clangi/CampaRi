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
! CONTRIBUTIONS: Albert Mao                                                !
!                                                                          !
!--------------------------------------------------------------------------!
!-----------------------------------------------------------------------
!
!    some routines that read-in auxiliary files
!    mostly pulled out of parsekey.f
!
!-----------------------------------------------------------------------
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
! a generic function to read an atom index list into array tmpvec (at least size n)
!
subroutine read_atmidx(iunit,nlst,tmpvec,tvsz,lsort,checkbounds)
!
  use iounit
  use atoms
  use interfaces
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,iomess,nlst,ii,mn,tvsz,i,j,jj
  integer tmpvec(tvsz)
  integer, ALLOCATABLE:: iv1(:)
  logical lsort,checkbounds,notdone
  character(MAXSTRLEN) str2
!
  nlst = 0
!
 79   format(FORM_MAXSTRLEN)
!
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_atmidx(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if (checkbounds.EQV..true.) then
      if ((mn.gt.n).OR.(mn.le.0)) then
        write(ilog,*) 'Requested set member ',&
 &mn,' for (atom) index list, which does not exist. Ignored ...'
        cycle
      end if
    end if
!
    nlst = nlst + 1
    if (checkbounds.EQV..true.) then
      if (nlst.gt.n) then
        write(ilog,*) 'Fatal. Set list provided to read_atmidx(...) is longer &
 &than the number of atoms in the system. This is fatal.'
        call fexit()
      end if
    end if
    if (nlst.gt.tvsz) then
      write(ilog,*) 'Fatal. List size violation in read_atmidx(...). Check for &
 &double entries in index input files.'
      call fexit()
    end if
    tmpvec(nlst) = mn
!
  end do
!
  if (nlst.le.1) return
!
  if (lsort.EQV..true.) then
!
    allocate(iv1(nlst))
    iv1(:) = tmpvec(1:nlst)
    ii = 1
    jj = nlst
    call merge_sort(ldim=nlst,up=lsort,list=iv1,olist=tmpvec(1:nlst),ilo=ii,ihi=jj)
    deallocate(iv1)
!   remove double entries
    notdone = .true.
    i = 2
    do while (notdone.EQV..true.)
      if (tmpvec(i).eq.tmpvec(i-1)) then
        write(ilog,*) 'Warning. Removing double entry from index list in read_atmidx(...).'
        tmpvec((i-1):(nlst-1)) = tmpvec(i:nlst)
        nlst = nlst - 1
      else
        i = i + 1
      end if
      if (i.gt.nlst) notdone = .false.
    end do
  else
!   remove double entries
    notdone = .true.
    i = 2
    do while (notdone.EQV..true.)
      do j=i,nlst
        if (tmpvec(j).eq.tmpvec(i-1)) then
          write(ilog,*) 'Warning. Removing double entry from index list in read_atmidx(...).'
          tmpvec((i-1):(nlst-1)) = tmpvec(i:nlst)
          nlst = nlst - 1
          exit
        end if
        if (j.eq.nlst) i = i + 1
      end do
      if (i.gt.nlst) notdone = .false.
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a generic function to read an pair index list into array tmpvec
!
subroutine read_atmpairidx(iunit,nlst,tmpvec,tvsz,lsort,checkbounds,swapem)
!
  use iounit
  use atoms
  use interfaces
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,iomess,nlst,ii,mn,mn2,tvsz,i,j,ista,kk
  integer tmpvec(tvsz,2)
  integer, ALLOCATABLE:: iv1(:,:)
  logical lsort,checkbounds,notdone,atrue,swapem
  character(MAXSTRLEN) str2,str3
  integer comp1
!
  nlst = 0
  atrue = .true.
  if (n.le.40000) then
    comp1 = n*(n-1)/2
  else
    comp1 = 40000*39999/2
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_atmpairidx(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((ii.ne.1).AND.(checkbounds.EQV..true.)) then
      if ((mn.gt.n).OR.(mn.le.0)) then
        write(ilog,*) 'Requested set member ',mn,' for (atom) pair index list, which does not exist. Ignored ...'
        cycle
      end if
    end if
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    ii = 1
    mn2 = -1
    call extract_int(str3,mn2,ii)
    if ((ii.ne.1).AND.(checkbounds.EQV..true.)) then
      if ((mn2.gt.n).OR.(mn2.le.0)) then
        write(ilog,*) 'Requested set member ',mn2,' for (atom) pair index list, which does not exist. Ignored ...'
        cycle
      end if
      if (mn.eq.mn2) then
        write(ilog,*) 'Requested identical indices (',mn,') for (atom) pair index list. Ignored ...'
        cycle
      end if
    end if
    if (ii.eq.1) then
      write(ilog,*) 'Fatal. Provided list does not (always) possess two entries in &
 &read_atmpairidx(...). Check input files.'
      call fexit()
    end if
!
    nlst = nlst + 1
    if (checkbounds.EQV..true.) then
      if (nlst.gt.comp1) then
        write(ilog,*) 'Fatal. Set list provided to read_atmpairidx(...) is longer &
 &than the number of unique atom pairs in the system or overflows otherwise. This is fatal.'
        call fexit()
      end if
    end if
    if (nlst.gt.tvsz) then
      write(ilog,*) 'Fatal. List size violation in read_atmpairidx (...). Check for &
 &double entries in index input files.'
      call fexit()
    end if
    if (swapem.EQV..true.) then
      tmpvec(nlst,1) = min(mn,mn2)
      tmpvec(nlst,2) = max(mn,mn2)
    else
      tmpvec(nlst,1) = mn
      tmpvec(nlst,2) = mn2
    end if
!
  end do
!
  if (nlst.le.1) return
!
  if (lsort.EQV..true.) then
!
!   sort first column
    allocate(iv1(nlst,4))
    iv1(:,1:2) = tmpvec(1:nlst,1:2)
    ii = 1
    kk = nlst
    do i=1,nlst
      iv1(i,3) = i
    end do
    call merge_sort(ldim=nlst,up=lsort,list=iv1(:,1),olist=tmpvec(1:nlst,1),ilo=ii,ihi=kk,idxmap=iv1(:,3),olist2=iv1(:,4))
    do i=1,nlst
      if (i.eq.iv1(i,4)) cycle
      tmpvec(i,2) = iv1(iv1(i,4),2)
    end do 
    deallocate(iv1)
!   resort second column for identical values in first
    i = 1
    ista = i
    do while (1.eq.1)
      if (i.lt.nlst) then
        if (tmpvec(i,1).eq.tmpvec(i+1,1)) then 
          i = i + 1
          cycle
        end if
      end if
      kk = i-ista+1
      if (kk.gt.1) then
        call isort(ldim=kk,up=atrue,ilist=tmpvec(ista:i,2),ilist2=tmpvec(ista:i,1))
      end if
      i = i + 1
      ista = i
      if (i.gt.nlst) exit
    end do
!
!   check for duplicates
    notdone = .true.
    i = 2
    do while (notdone.EQV..true.)
      if ((tmpvec(i,1).eq.tmpvec(i-1,1)).AND.(tmpvec(i,2).eq.tmpvec(i-1,2))) then
        write(ilog,*) 'Warning. Removing double entry from index list in read_atmpairidx(...).'
        tmpvec((i-1):(nlst-1),:) = tmpvec(i:nlst,:)
        nlst = nlst - 1
      else
        i = i + 1
      end if
      if (i.gt.nlst) notdone = .false.
    end do
  else
!
!   check for duplicates
    notdone = .true.
    i = 2
    do while (notdone.EQV..true.)
      do j=i,nlst
        if (((tmpvec(j,1).eq.tmpvec(i-1,1)).AND.(tmpvec(j,2).eq.tmpvec(i-1,2))).OR.&
 &          ((tmpvec(j,2).eq.tmpvec(i-1,1)).AND.(tmpvec(j,1).eq.tmpvec(i-1,2)))) then
          write(ilog,*) 'Warning. Removing double entry from index list in read_atmpairidx(...).'
          tmpvec((i-1):(nlst-1),:) = tmpvec(i:nlst,:)
          nlst = nlst - 1
          exit
        end if
        if (j.eq.nlst) i = i + 1
      end do
      if (i.gt.nlst) notdone = .false.
    end do
  end if
!
end
!
!----------------------------------------------------------------------------
!
subroutine read_alignfile(algstruc)
!
  use iounit
  use atoms
  use system
  use clusters
  use pdb, ONLY: pdb_rmol
  use sequen, ONLY: molofrs
  use molecule, ONLY: nmol
!
  implicit none
!
  integer iunit,freeunit,t1,t2,i,lmol
  logical exists,atrue
  integer, ALLOCATABLE:: tmpvec(:),cnts(:)
  type(t_align) algstruc
!
  call strlims(algstruc%filen,t1,t2)
  inquire(file=algstruc%filen(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &algstruc%filen(t1:t2),') for FMCSC_ALIGNFILE. Turning off alignment for trajectory analysis.'
    algstruc%yes = .false.
    algstruc%calc = nsim + 1
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  atrue = .true.
  open (unit=iunit,file=algstruc%filen(t1:t2),status='old')
!
  allocate(tmpvec(n))
!
! get the unique, sorted set of atom indices
  call read_atmidx(iunit,algstruc%nr,tmpvec,n,atrue,atrue)
  close(unit=iunit)
  if (algstruc%nr.lt.3) then
    write(ilog,*) 'Set of atom indices is too small (',algstruc%nr,') for successful alignment. &
 &Disabling this functionality.'
    algstruc%yes = .false.
    algstruc%calc = nsim + 1 
    deallocate(tmpvec)
    return
  end if
  allocate(algstruc%set(algstruc%nr))
  allocate(algstruc%refxyz(3*algstruc%nr))
  allocate(algstruc%curxyz(3*algstruc%nr))
  algstruc%set(1:algstruc%nr) = tmpvec(1:algstruc%nr)
!
  if ((algstruc%mmol.eq.0).AND.(pdb_rmol.gt.0).AND.(pdb_rmol.le.nmol)) then
    algstruc%mmol = pdb_rmol
    if ((nmol.gt.2).AND.(bnd_type.eq.1)) write(ilog,*) 'The reference molecule for alignment (FMCSC_ALIGNCALC) is &
 &provided by the setting for FMCSC_XYZ_REFMOL (',pdb_rmol,').'
  else if (algstruc%mmol.eq.0) then ! find the molecule that contributes most to the alignment set
    allocate(cnts(nmol))
    cnts(:) = 0
    do i=1,algstruc%nr
      cnts(molofrs(atmres(algstruc%set(i)))) = cnts(molofrs(atmres(algstruc%set(i)))) + 1
    end do
    lmol = 0
    do i=1,nmol
      if (cnts(i).gt.lmol) then
        lmol = cnts(i)
        algstruc%mmol = i
      end if
    end do
    deallocate(cnts)
    if ((nmol.gt.2).AND.(bnd_type.eq.1)) write(ilog,*) 'The reference molecule for alignment (FMCSC_ALIGNCALC) was &
 &determined heuristically and is # ',algstruc%mmol,' by occurrence in the sequence file.'
  end if
  deallocate(tmpvec)
!
end
!
!----------------------------------------------------------------------------
!
! this is the routine for reading the dimension selection input file for a clustering (FMCSC_CFILE)
! it performs a number of setup tasks and sanity checks
!
! it is somewhat confusingly named as there is also a routine to read a file with an actual 
! clustering result (read_clufile(...) below)
!
subroutine read_clusteringfile()
!
  use iounit
  use atoms
  use system
  use clusters
  use sequen
  use aminos
  use fyoc
  use params
  use interfaces
  use polypep
  use zmatrix, ONLY: izrot
  use movesets
  use molecule, ONLY: atmol
  use forces, ONLY: dyn_integrator_ops
#ifdef ENABLE_MPI
  use mpistuff
  use mcsums
#endif
!
  implicit none
!
  integer i,j,k,rs,ttc,ttc2,iunit,freeunit,t1,t2,kk,k1,k2,modstep,maxstp,minstp,kidx
  logical exists,atrue
  RTYPE random
  integer, ALLOCATABLE:: tmpvec(:),tmpvec2(:,:)
  type(t_align) algi
#ifdef ENABLE_MPI
  integer masterrank
#endif
!
  algi%yes = .false.
  algi%nr = 0
!
  cmaxsnaps = 0
!
! parsefiles: test
  maxstp = nsim
  minstp = min(nsim,nequil-1)
  do i=minstp,maxstp
    modstep = mod(i,cstorecalc)
    if ((i.gt.nequil).AND.(modstep.eq.0)) then
      cmaxsnaps = cmaxsnaps + 1
    end if
  end do
#ifdef ENABLE_MPI
  if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.).AND.(pdb_analyze.EQV..false.)) then
    if (re_freq.lt.cstorecalc) then
      write(ilog,*) 'Fatal. With the specified settings, some or all of the reseeding intervals for PIGS collect &
 &no data whatsoever. Please check settings for FMCSC_CCOLLECT and FMCSC_REFREQ.'
      call fexit()
    end if
    if ((re_freq.lt.10*cstorecalc).AND.(re_freq.ge.cstorecalc).AND.(re_freq.lt.nsim)) then
      write(ilog,*) 'Warning. With the specified settings, the reseeding intervals for PIGS collect less than 10 snapshots &
 &per replica. Are the values for FMCSC_CCOLLECT and FMCSC_REFREQ as intended?'
    end if
    if ((re_freq.lt.nsim).AND.(rstout.le.nsim).AND.(mod(rstout,re_freq).ne.0)) then
      write(ilog,*) 'Warning. To avoid loss of information on restarts, the interval for writing restart files should be &
 &an integer multiple of the reseeding interval. Consider changing this.'
    end if
    cmaxsnaps = 0
    do i=cstorecalc,re_freq+cstorecalc-1
      modstep = mod(i,cstorecalc)
      if (modstep.eq.0) then
        cmaxsnaps = cmaxsnaps + 1
      end if
    end do
  end if
#endif
#ifdef ENABLE_MPI
  if (.NOT.((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.))) then
    if (cmaxsnaps.le.3) then
      write(ilog,*) 'Disabling structural clustering due to insufficient number &
 &of snapshots being analyzed.'
      cstorecalc = nsim + 1
      return
    end if
  end if
#else
  if (cmaxsnaps.le.3) then
    write(ilog,*) 'Disabling structural clustering due to insufficient number &
 &of snapshots being analyzed.'
    cstorecalc = nsim + 1
    return
  end if
#endif
#ifdef ENABLE_MPI
  masterrank = 0
  if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..false.)) then
    if (do_restart.EQV..true.) then
      write(ilog,*) 'Warning. Structural clustering analysis is not supported in restarted MPI averaging runs. Disabled.'
      cstorecalc = nsim + 1
      return
    else
      write(ilog,*) 'Warning. For structural clustering analysis in MPI averaging calculations, memory for trajectory storage &
 &is allocated exclusively on the head node. This may lead to unexpected crashes.'
    end if
  end if
#endif
  if (do_restart.EQV..true.) then
    write(ilog,*) 'Warning. When using structural clustering analysis for a restarted simulation run, memory may be &
 &drastically overallocated. It is generally not recommended to use these features together.'
  end if
! out of place
  if ((cdis_crit.eq.4).OR.(cdis_crit.eq.2).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    sumssz = 5
  else
    sumssz = 2
  end if
  if ((cdis_crit.eq.6).AND.(align_for_clustering.EQV..false.)) then
    write(ilog,*) 'Warning. It is pointless to request mode 6 for FMCSC_CDISTANCE if alignment is not performed.'
  end if
  if ((cdis_crit.eq.10).AND.(align_for_clustering.EQV..true.)) then
    write(ilog,*) 'Fatal. Locally adaptive weights on coordinate RMSD (10 for FMCSC_CDISTANCE) are not compatible with &
 &alignment at the moment.'
    call fexit()
  end if
  if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then 
    if (n.le.1) then
      write(ilog,*) 'Fatal. It is nonsensical to use interatomic distances (7-9 for FMCSC_CDISTANCE) for clustering with &
 &just a single atom in the system.'
      call fexit()
    end if
  end if
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4)) then 
    if (n_crosslinks.gt.0) then
      write(ilog,*) 'Warning. Values used for rotational masses are always nonsensical for residues participating in stretches &
 &containing crosslinks. If any of these are selected, it is strongly recommended to use other weights (see FMCSC_CMODWEIGHTS).'
    end if
  end if
!
  call strlims(cfilen,t1,t2)
  inquire(file=cfilen(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &cfilen(t1:t2),') for FMCSC_CFILE. Assuming full-size request (all-atom, all-torsion, etc. ...).'
    if ((cdis_crit.ge.1).AND.(cdis_crit.le.4)) then
      ttc = 0
      do rs=1,nseq
        if (wline(rs).gt.0) ttc = ttc + 1
        if (fline(rs).gt.0) ttc = ttc + 1
        if (yline(rs).gt.0) ttc = ttc + 1
        ttc = ttc + nchi(rs) + nnucs(rs)
        if (seqpolty(rs).eq.'N') then
          ttc = ttc + 1
        end if
        if (disulf(rs).gt.0) then
          if (crosslink(crlk_idx(rs))%itstype.le.2) then
            ttc = ttc + 1
            if (disulf(rs).gt.rs) ttc = ttc + 1
          end if
        end if
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
              end if
            end if
          end do
        end if
      end do
      if (ttc.le.0) then
        write(ilog,*) 'No torsions selected for structural clustering.'
        if (pdb_analyze.EQV..true.) then
          write(ilog,*) 'This is fatal in trajectory analysis runs.'
          call fexit()
        else
          write(ilog,*) 'Ignoring request to perform structural clustering.'
          cstorecalc = nsim + 1
          return
        end if
      end if
      if (cdis_crit.eq.1) calcsz = ttc ! fync
      if (cdis_crit.eq.2) calcsz = 2*ttc ! LAW fync
      if (cdis_crit.eq.3) calcsz = 2*ttc ! sin/cos
      if (cdis_crit.eq.4) calcsz = 3*ttc ! LAW sin/cos
    else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then ! RMSD
      if ((cdis_crit.eq.6).AND.(align_for_clustering.EQV..true.)) then
        write(ilog,*) 'Warning. It is not possible to split alignment and distance computation atom sets for &
 &FMCSC_CDISTANCE 6 if FMCSC_CFILE is not specified or corrupt. Using complete atom set for both.'
      end if
      ttc = n
      cdofsbnds(1) = 1
      cdofsbnds(3) = 1
      cdofsbnds(2) = 3*n
      cdofsbnds(4) = 3*n
      calcsz = 3*n
      if (cdis_crit.eq.10) calcsz = 2*calcsz
    else if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then ! DRMS
      if (n.le.10) then
        ttc = min(3*n,n*(n-1)/2)
      else
        ttc = 3*n
      end if
      write(ilog,*) 'Due to missing input for FMCSC_CFILE for distance criteria 7-9, a set of ',ttc,' interatomic &
 &distances with randomly selected atoms will be used.'
      calcsz = ttc
      if (cdis_crit.eq.9) calcsz = 2*calcsz
    else
      write(ilog,*) 'Fatal. Encountered unsupported mode for data collection for &
 &structural clustering in read_clusteringfile(). This is an omission bug.'
      call fexit()
    end if
    if (cdis_crit.eq.4) then
      clstsz =  2*ttc
    else if ((cdis_crit.ge.5).AND.(cdis_crit.le.6)) then
      clstsz = calcsz
    else if ((cdis_crit.eq.2).OR.(cdis_crit.eq.9)) then
      clstsz = ttc
    else if (cdis_crit.eq.3) then
      clstsz = calcsz
    else if (cdis_crit.eq.10) then
      clstsz = 3*n
    else
      clstsz = calcsz
    end if
    if (((cdis_crit.eq.5).OR.(cdis_crit.eq.6)).AND.(calcsz.lt.9).AND.(align_for_clustering.EQV..true.)) then
      write(ilog,*) 'Set of atom indices is too small (',calcsz/3,') for successful alignment required &
 &for Cartesian coordinate-based, structural clustering.'
      if (pdb_analyze.EQV..true.) then
        write(ilog,*) 'This is fatal in trajectory analysis runs.'
        call fexit()
      else
        write(ilog,*) 'Ignoring request to perform structural clustering.'
        cstorecalc = nsim + 1
        return
      end if
    end if
#ifdef ENABLE_MPI
    if (use_MPIAVG.EQV..true.) then
      if (myrank.eq.masterrank) then
        allocate(cludata(calcsz,mpi_nodes*cmaxsnaps))
        allocate(cl_mwvec(calcsz))
        cl_mwvec(:) = 0.0
      end if
    else
      allocate(cludata(calcsz,cmaxsnaps))
      allocate(cl_mwvec(calcsz))
      cl_mwvec(:) = 0.0
    end if
#else
    allocate(cludata(calcsz,cmaxsnaps))
    allocate(cl_mwvec(calcsz))
    cl_mwvec(:) = 0.0
#endif
    allocate(cdofset(ttc,2))
    if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then ! DRMS
      do k=1,clstsz
        ttc = 1
        do while (ttc.ne.0)
          cdofset(k,1) = max(1,min(n,ceiling(random()*n)))
          j = cdofset(k,1) 
          do while (j.eq.cdofset(k,1))
            j = max(1,min(n,ceiling(random()*n)))
          end do
          cdofset(k,2) = j
          if (cdofset(k,1).gt.cdofset(k,2)) then
            kk = cdofset(k,1)
            cdofset(k,1) = cdofset(k,2)
            cdofset(k,2) = kk
          end if
          ttc = 0
          do kk=1,k-1
            if ((cdofset(k,1).eq.cdofset(kk,1)).AND.(cdofset(k,2).eq.cdofset(kk,2))) then
              ttc = kk
              exit
            end if
          end do 
        end do 
      end do
 98   format('Atoms ',i6,' (',a4,') in Res. ',i6,' (',a3,') and ',i6,' (',a4,') in Res. ',i6,' (',a3,')')
      write(ilog,*)
      write(ilog,*) 'Summary of randomly chosen interatomic distances for structural clustering:'
      write(ilog,*) '---------------------------------------------------------------------------'
      do k=1,clstsz
        write(ilog,98) cdofset(k,1),bio_code(b_type(cdofset(k,1))),atmres(cdofset(k,1)),amino(seqtyp(atmres(cdofset(k,1)))),&
 &                     cdofset(k,2),bio_code(b_type(cdofset(k,2))),atmres(cdofset(k,2)),amino(seqtyp(atmres(cdofset(k,2))))
      end do
      write(ilog,*)
      if (cdis_crit.eq.8) then
        allocate(cl_imvec(calcsz))
        cl_imvec(:) = mass(cdofset(:,1))+mass(cdofset(:,2))
      end if
    else
      do i=1,ttc
        cdofset(i,1) = i
      end do
    end if
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  atrue = .true.
  open (unit=iunit,file=cfilen(t1:t2),status='old')
!
! get the unique, sorted set of indices
  if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
    if (n.le.1000) then
      kk = n*(n-1)/2
    else
      kk = 1000*999/2
    end if
    allocate(tmpvec2(kk,2))
    calcsz = 0
    call read_atmpairidx(iunit,calcsz,tmpvec2,kk,atrue,atrue,atrue)
    close(unit=iunit)
!
    if (calcsz.le.0) then
      write(ilog,*) 'No distances selected for structural clustering.'
      if (pdb_analyze.EQV..true.) then
        write(ilog,*) 'This is fatal in trajectory analysis runs.'
        call fexit()
      else
        write(ilog,*) 'Ignoring request to perform structural clustering.'
        cstorecalc = nsim + 1
        deallocate(tmpvec2)
        return
      end if
    end if
    allocate(cdofset(calcsz,2))
    cdofset(1:calcsz,:) = tmpvec2(1:calcsz,:)
    clstsz = calcsz
    if (cdis_crit.eq.9) calcsz = 2*calcsz
#ifdef ENABLE_MPI
    if (use_MPIAVG.EQV..true.) then
      if (myrank.eq.masterrank) then
        allocate(cludata(calcsz,mpi_nodes*cmaxsnaps))
        allocate(cl_mwvec(calcsz))
        cl_mwvec(:) = 0.0
      end if
    else
      allocate(cludata(calcsz,cmaxsnaps))
      allocate(cl_mwvec(calcsz))
      cl_mwvec(:) = 0.0
    end if
#else
    allocate(cludata(calcsz,cmaxsnaps))
    allocate(cl_mwvec(calcsz))
    cl_mwvec(:) = 0.0
#endif
    deallocate(tmpvec2)
    if (cdis_crit.eq.8) then
      allocate(cl_imvec(calcsz))
      cl_imvec(:) = mass(cdofset(:,1))+mass(cdofset(:,2))
    end if
!
  else
!   for aligned coordinates, offer option to separate alignment and distance sets via FMCSC_ALIGNFILE
!   we'll merge these sets later and set up corresponding weights
    if ((cdis_crit.eq.6).AND.(align_for_clustering.EQV..true.)) then
      call strlims(align%filen,t1,t2)
      inquire(file=align%filen(t1:t2),exist=exists)
      if (exists.EQV..true.) then
        algi%filen = align%filen
        algi%yes = .true.
        call read_alignfile(algi)
        if (allocated(algi%refxyz).EQV..true.) deallocate(algi%refxyz)
        if (allocated(algi%curxyz).EQV..true.) deallocate(algi%curxyz)
      end if
    end if
    allocate(tmpvec(n))
    calcsz = 0
    call read_atmidx(iunit,calcsz,tmpvec,n,atrue,atrue)
    close(unit=iunit)
!
!   sanity checks
    if ((cdis_crit.ge.1).AND.(cdis_crit.le.4)) then
      ttc = 0
      do rs=1,nseq
        if (wline(rs).gt.0) ttc = ttc + 1
        if (fline(rs).gt.0) ttc = ttc + 1
        if (yline(rs).gt.0) ttc = ttc + 1
        ttc = ttc + nchi(rs) + nnucs(rs)
        if (seqpolty(rs).eq.'N') then
          ttc = ttc + 1
        end if
        if (disulf(rs).gt.0) then
          if (crosslink(crlk_idx(rs))%itstype.le.2) then
            ttc = ttc + 1
            if (disulf(rs).gt.rs) ttc = ttc + 1
          end if
        end if
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
              end if
            end if
          end do
        end if
      end do
      i = 1
      do while (i.le.calcsz)
        if (tmpvec(i).gt.ttc) then
          write(ilog,*) 'Warning. Requested illegal torsional index in read_clusteringfile(). Please &
 &consult the documentation on how to access specific torsions by index request. Ignored.'
          do j=i+1,calcsz
            tmpvec(j-1) = tmpvec(j)
          end do
          calcsz = calcsz - 1
        else
          i = i + 1
        end if
      end do
      ttc2 = 1
      ttc = 0
 66   format('# ',i5,' is ',a7,'   angle of residue # ',i7,' (',a3,')')
 69   format('# ',i5,' is ',a7,'   angle of residues # ',i7,' (',a3,') and ',i7,' (',a3,')')
 67   format('# ',i5,' is ',a7,'-',i1,' angle of residue # ',i7,' (',a3,')')
 68   format('# ',i5,' is ',a7,' (',i10,') angle of residue # ',i7,' (',a3,')')

      write(ilog,*)
      write(ilog,*) 'Summary of chosen dihedral angles for structural clustering:'
      write(ilog,*) '------------------------------------------------------------'
      do rs=1,nseq
        if (wline(rs).gt.0) ttc = ttc + 1
        if (tmpvec(ttc2).eq.ttc) then
          write(ilog,66) ttc2,'  OMEGA',rs,amino(seqtyp(rs))
          ttc2 = ttc2 + 1
        end if
        if (fline(rs).gt.0) ttc = ttc + 1
        if (tmpvec(ttc2).eq.ttc) then
          write(ilog,66) ttc2,'    PHI',rs,amino(seqtyp(rs))
          ttc2 = ttc2 + 1
        end if
        if (yline(rs).gt.0) ttc = ttc + 1
        if (tmpvec(ttc2).eq.ttc) then
          write(ilog,66) ttc2,'    PSI',rs,amino(seqtyp(rs))
          ttc2 = ttc2 + 1
        end if
        do j=1,nnucs(rs)
          ttc = ttc + 1
          if (tmpvec(ttc2).eq.ttc) then
            write(ilog,67) ttc2,'NUCLEIC',j,rs,amino(seqtyp(rs))
            ttc2 = ttc2 + 1
          end if
        end do
        if (seqpolty(rs).eq.'N') then
          ttc = ttc + 1
          if (tmpvec(ttc2).eq.ttc) then
            write(ilog,66) ttc2,'  SUGAR',rs,amino(seqtyp(rs))
            ttc2 = ttc2 + 1
          end if
        end if
        do j=1,nchi(rs)
          ttc = ttc + 1
          if (tmpvec(ttc2).eq.ttc) then
            write(ilog,67) ttc2,'    CHI',j,rs,amino(seqtyp(rs))
            ttc2 = ttc2 + 1
          end if
        end do
        if (disulf(rs).gt.0) then
          if (crosslink(crlk_idx(rs))%itstype.le.2) then
            ttc = ttc + 1
            if (tmpvec(ttc2).eq.ttc) then
              j = 2
              write(ilog,67) ttc2,'    CHI',j,rs,amino(seqtyp(rs))
              ttc2 = ttc2 + 1
            end if
            if (disulf(rs).gt.rs) then
              ttc = ttc + 1
              if (tmpvec(ttc2).eq.ttc) then
                write(ilog,69) ttc2,'   CSSC',rs,amino(seqtyp(rs)),disulf(rs),amino(seqtyp(disulf(rs)))
                ttc2 = ttc2 + 1
              end if
            end if
          end if
        end if
        if (((dyn_integrator_ops(10).eq.3).OR.((dyn_integrator_ops(10).eq.2).AND.(seqtyp(rs).ne.26)).OR.&
                                              ((dyn_integrator_ops(10).eq.1).AND.(seqtyp(rs).eq.26))).AND.&
   &        ((unslst%nr.gt.0).OR.(unklst%nr.gt.0))) then
          j = 0
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
                j = j + 1
                if (tmpvec(ttc2).eq.ttc) then
                  write(ilog,68) ttc2,'   UNS',kk,rs,amino(seqtyp(rs))
                  ttc2 = ttc2 + 1
                end if
              end if
            end if
          end do
        end if
      end do
      write(ilog,*) '------------------------------------------------------------'
      write(ilog,*)
    end if
    if (calcsz.le.0) then
      write(ilog,*) 'No torsions/atoms selected for structural clustering.'
      if (pdb_analyze.EQV..true.) then
        write(ilog,*) 'This is fatal in trajectory analysis runs.'
        call fexit()
      else
        write(ilog,*) 'Ignoring request to perform structural clustering.'
        cstorecalc = nsim + 1
        deallocate(tmpvec)
        return
      end if
    end if
!
    if (((cdis_crit.eq.5).OR.(cdis_crit.eq.6)).AND.(calcsz.lt.3)) then
      if (align_for_clustering.EQV..true.) then
        write(ilog,*) 'Set of atom indices is too small (',calcsz,') for successful alignment required &
 &for Cartesian coordinate-based, structural clustering.'
        if (pdb_analyze.EQV..true.) then
          write(ilog,*) 'This is fatal in trajectory analysis runs.'
          call fexit()
        else
          write(ilog,*) 'Ignoring request to perform structural clustering.'
          cstorecalc = nsim + 1
          deallocate(tmpvec)
          return
        end if
      end if
    end if
!
    if ((cdis_crit.eq.6).AND.(algi%yes.EQV..true.)) then
      k1 = 1
      k2 = 1
      ttc = 0
!     create a joint set with the overlaps in the middle, low end for distance, high end for alignment set
      cdofsbnds(1) = 1
      allocate(cdofset(calcsz+algi%nr,2))
      do while ((k1.le.calcsz).OR.(k2.le.algi%nr))
        if (k2.gt.algi%nr) then
          ttc = ttc + 1
          cdofset(ttc,1) = tmpvec(k1)
          k1 = k1 + 1
        else if (k1.gt.calcsz) then
          exit
        else if (tmpvec(k1).eq.algi%set(k2)) then
          k1 = k1 + 1
          k2 = k2 + 1
        else if (tmpvec(k1).lt.algi%set(k2)) then
          ttc = ttc + 1
          cdofset(ttc,1) = tmpvec(k1)
          k1 = k1 + 1
        else
          k2 = k2 + 1
        end if
      end do
      k1 = 1
      k2 = 1
      cdofsbnds(3) = 3*ttc + 1
      do while ((k1.le.calcsz).OR.(k2.le.algi%nr))
        if ((k2.gt.algi%nr).OR.(k1.gt.calcsz)) then
          exit
        else if (tmpvec(k1).eq.algi%set(k2)) then
          ttc = ttc + 1
          cdofset(ttc,1) = tmpvec(k1)
          k1 = k1 + 1
          k2 = k2 + 1
        else if (tmpvec(k1).lt.algi%set(k2)) then
          k1 = k1 + 1
        else
          k2 = k2 + 1
        end if
      end do
      k1 = 1
      k2 = 1
      cdofsbnds(2) = 3*ttc
      do while ((k1.le.calcsz).OR.(k2.le.algi%nr))
        if (k2.gt.algi%nr) then
          exit
        else if (k1.gt.calcsz) then
          ttc = ttc + 1
          cdofset(ttc,1) = algi%set(k2)
          k2 = k2 + 1
        else if (tmpvec(k1).eq.algi%set(k2)) then
          k1 = k1 + 1
          k2 = k2 + 1
        else if (tmpvec(k1).lt.algi%set(k2)) then
          k1 = k1 + 1
        else
          ttc = ttc + 1
          cdofset(ttc,1) = algi%set(k2)
          k2 = k2 + 1
        end if
      end do
      cdofsbnds(4) = 3*ttc
      deallocate(tmpvec)
      allocate(tmpvec(ttc))
      tmpvec(:) = cdofset(1:ttc,1)
      deallocate(cdofset)
      allocate(cdofset(ttc,2))
      cdofset(:,1) = tmpvec(:)
      calcsz = ttc
      deallocate(algi%set)
    else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
      if ((cdis_crit.eq.6).AND.(align_for_clustering.EQV..true.)) then
        write(ilog,*) 'Warning. Failed to read in a separate atom index set for alignment to use with &
 &FMCSC_CDISTANCE 6. Using same set for alignment and distance computation.'
      end if
      allocate(cdofset(calcsz,2))
      cdofset(1:calcsz,1) = tmpvec(1:calcsz)
      cdofsbnds(1) = 1
      cdofsbnds(2) = 3*calcsz
      cdofsbnds(3:4) = cdofsbnds(1:2)
    else
      allocate(cdofset(calcsz,2))
      cdofset(1:calcsz,1) = tmpvec(1:calcsz)
    end if
    if (cdis_crit.eq.4) then
      calcsz = calcsz*3
      clstsz = (calcsz*2)/3
    else if ((cdis_crit.ge.5).AND.(cdis_crit.le.6)) then
      calcsz = calcsz*3
      clstsz = calcsz
    else if (cdis_crit.eq.10) then
      calcsz = calcsz*3
      clstsz = calcsz
      calcsz = 2*calcsz
    else if (cdis_crit.eq.2) then
      clstsz = calcsz
      calcsz = calcsz*2
    else if (cdis_crit.eq.3) then
      calcsz = calcsz*2
      clstsz = calcsz
    else ! option 1
      clstsz = calcsz
    end if
#ifdef ENABLE_MPI
    if (use_MPIAVG.EQV..true.) then
      if (myrank.eq.masterrank) then
        allocate(cludata(calcsz,mpi_nodes*cmaxsnaps))
        allocate(cl_mwvec(calcsz))
        cl_mwvec(:) = 0.0
      end if
    else
      allocate(cludata(calcsz,cmaxsnaps))
      allocate(cl_mwvec(calcsz))
      cl_mwvec(:) = 0.0
    end if
#else
    allocate(cludata(calcsz,cmaxsnaps))
    allocate(cl_mwvec(calcsz))
    cl_mwvec(:) = 0.0
#endif
    deallocate(tmpvec)
!
  end if
!
end
!
!------------------------------------------------------------------------------------------
!
! this is the routine for reading back an existing clustering (such as in STRUCT_CLUSTERING.clu)
! the first two arguments specify the variable to hold the integer (coarse-grained) trajectory
! which refers to the column in the file with name cfn
!
subroutine read_clufile(cgtrj,cgtrjsz,cfn,which)
!
  use iounit
  use interfaces
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer, INTENT(IN):: which,cgtrjsz
  character(MAXSTRLEN), INTENT(IN):: cfn
!
  integer i,s1,s2,iunit,freeunit,tmplst(which),vread,vreq,j,cgtrj(cgtrjsz),iomess,ii,jj,kk
  character(which*MAXSTRLEN) line
  logical afalse,atrue,hw,exists
  integer, ALLOCATABLE:: iv1(:,:),iv2(:)
!
  afalse = .false.
  atrue = .true.
  vreq = which
!
  call strlims(cfn,s1,s2)
  inquire(file=cfn(s1:s2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file (',cfn(s1:s2),') for FMCSC_CLUFILE while in mode 6/7 for FMCSC_CMODE.'
    call fexit()
  end if
!
  iunit = freeunit()
  open(unit=iunit,file=cfn(s1:s2),status='old')
!
  cgtrj(:) = -1
  do i=1,cgtrjsz
    read(iunit,'(a)',iostat=iomess) line
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Fatal. Failed in reading line ',i,' of  clustering input file (',cfn(s1:s2),'). Expected ',&
 &cgtrjsz,' lines to be present.'
      call fexit()
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_clufile(...).'
      call fexit()
    end if
    j = 1
    call get_ints_from_str(line,vreq,tmplst(:),j,vread,afalse)
    if (vread.lt.which) then
      write(ilog,*) 'Fatal. Line ',i,' of clustering input file (',cfn(s1:s2),') does not have the required number &
 &of values for the setting of FMCSC_CLUFILECOL.'
      call fexit()
    end if
    cgtrj(i) = tmplst(which)
  end do
!
  if (MINVAL(cgtrj(1:cgtrjsz)).le.0) then
    write(ilog,*) 'Fatal. Cluster indices in input file ',cfn(s1:s2),' (column ',which,') must be positive integers.'
    call fexit()
  end if
!
  allocate(iv1(cgtrjsz,2))
  iv1(:,1) = cgtrj(1:cgtrjsz)
  ii = 1
  jj = cgtrjsz
  j = cgtrjsz
  call merge_sort(ldim=j,up=atrue,list=iv1(:,1),olist=iv1(:,2),ilo=ii,ihi=jj)
!
  kk = 1
  do i=2,cgtrjsz
    if ((iv1(i,2)-iv1(i-1,2)).ge.1) kk = kk + 1
  end do
  if (MAXVAL(cgtrj(1:cgtrjsz)).gt.cgtrjsz) then
    write(ilog,*) 'Warning. The highest cluster number referred to in input file ',cfn(s1:s2),' (column ',which,') exceeds &
 &the number of snapshots. If input is arbitrary, this could lead to a memory-related crash.'
  end if
  allocate(iv2(MAXVAL(cgtrj(1:cgtrjsz))))
  kk = 1
  do i=2,cgtrjsz
    if ((iv1(i,2)-iv1(i-1,2)).ge.1) then
      iv2(iv1(i-1,2)) = kk
      kk = kk + 1
    end if
  end do
  iv2(iv1(cgtrjsz,2)) = kk
!
  hw = .false.
  do i=1,cgtrjsz
    if ((hw.EQV..false.).AND.(iv2(cgtrj(i)).ne.cgtrj(i))) then
      write(ilog,*) 'Warning. Cluster numbers have been remapped (non-contiguous input in ',cfn(s1:s2),', column ',which,'). &
 &This can be confusing.'
      hw = .true.
    end if
    cgtrj(i) = iv2(cgtrj(i))
  end do
!
  deallocate(iv2)
  deallocate(iv1)
!
  close(unit=iunit)
! 
end
!
!----------------------------------------------------------------------------
!
subroutine read_trajbrksfile()
!
  use clusters
  use iounit
  use system
#ifdef ENABLE_MPI
  use mpistuff, ONLY: re_conditions,use_MPIAVG
#else
  use pdb, ONLY: framecnt,select_frames,pdb_fileformat,use_frame_weights
#endif
!
  implicit none
!
  integer k,k2,t1,t2,iunit,i,j,freeunit,vmax
  integer, ALLOCATABLE:: tmpvec(:)
  logical exists,atrue,afalse
!
  atrue = .true.
  afalse = .false.
  vmax = nsim
#ifdef ENABLE_MPI  
  if (use_MPIAVG.EQV..true.) then
    vmax = re_conditions*nsim
  end if
#else
 if ((select_frames.EQV..true.).AND.(((pdb_fileformat.eq.5).OR.(pdb_fileformat.eq.2)).AND.(use_frame_weights.EQV..false.))) then
    vmax = framecnt
  end if
#endif
  ntbrks = 0
!
  call strlims(tbrkfilen,t1,t2)
  inquire(file=tbrkfilen(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &tbrkfilen(t1:t2),') for FMCSC_TRAJBREAKSFILE. All transitions will be kept (this may include replica-exchange swaps,&
 & transitions caused by trajectory concatenation, etc).'
    return
  end if
  iunit = freeunit()
  open(unit=iunit,file=tbrkfilen(t1:t2),status='old')
!
  allocate(tmpvec(vmax))
!
  call read_atmidx(iunit,ntbrks,tmpvec,vmax,atrue,afalse)
  close(unit=iunit)
!
  k2 = ntbrks
  k = ntbrks
  do i=1,k
    if (tmpvec(i).ge.1) then
      if (i.gt.1) then
        do j=1,k-(i-1)
          tmpvec(j) = tmpvec(j+i-1)
        end do
        ntbrks = ntbrks - (i-1)
      end if
      exit
    end if
  end do
  k = ntbrks
  do i=1,k
    if (tmpvec(i).gt.vmax) then
      ntbrks = i-1
      exit
    end if
  end do
!
  if (ntbrks.lt.k2) then
    write(ilog,*) 'Warning. ',k2-ntbrks,' trajectory breaks in ',tbrkfilen(t1:t2),' have been discarded (out of range).'
  end if
!
  if (ntbrks.le.0) then
    deallocate(tmpvec)
    write(ilog,*) 'Warning. Set of user-requested trajectory breaks for network-related analyses is empty.'
    return
  end if
!
  allocate(trbrkslst(ntbrks))
  trbrkslst(1:ntbrks) = tmpvec(1:ntbrks)
  deallocate(tmpvec)
!
end
!
!----------------------------------------------------------------------------
!
subroutine read_trajlinksfile()
!
  use clusters
  use iounit
  use system
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
#ifdef ENABLE_MPI
  use mpistuff, ONLY: re_conditions,use_MPIAVG
#else
  use pdb, ONLY: framecnt,select_frames,pdb_fileformat,use_frame_weights
#endif
!
  implicit none
!
  integer k,t1,t2,iunit,i,freeunit,iomess,ii,vmax
  integer, ALLOCATABLE:: tmpvec(:,:)
  logical exists,atrue,afalse
  character(MAXSTRLEN) str2,abc
!
  atrue = .true.
  afalse = .false.
  vmax = nsim
#ifdef ENABLE_MPI  
  if (use_MPIAVG.EQV..true.) then
    vmax = re_conditions*nsim
  end if
#else
  if ((select_frames.EQV..true.).AND.(((pdb_fileformat.eq.5).OR.(pdb_fileformat.eq.2)).AND.(use_frame_weights.EQV..false.))) then
    vmax = framecnt
  end if
#endif
  ntlnks = 0
!
  call strlims(tlnkfilen,t1,t2)
  inquire(file=tlnkfilen(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &tlnkfilen(t1:t2),') for FMCSC_TRAJLINKSFILE. No transitions will be added.'
    return
  end if
  iunit = freeunit()
  open(unit=iunit,file=tlnkfilen(t1:t2),status='old')
!
  read(iunit,'(a)',iostat=iomess) str2
!
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file (',&
 &tlnkfilen(t1:t2),') for FMCSC_TRAJLINKSFILE. No transitions will be added.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing input file for FMCSC_TRAJLINKSFILE (got: ',tlnkfilen(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  abc(:) = ' '
  call extract_str(str2,abc,ii)
  if ((abc(1:1).eq.'Y').OR.(abc(1:1).eq.'y')) then
    remap_trajlnks = .true.
  else if ((abc(1:1).eq.'N').OR.(abc(1:1).eq.'n')) then
    remap_trajlnks = .false.
  else
    write(ilog,*) 'Warning. Got bad mode identifier in first line of FMCSC_TRAJLINKSFILE (',tlnkfilen(t1:t2),'). The first &
 &line should be either "Y" or "N" (without quotes).'
    call fexit()
  end if
!
  allocate(tmpvec(vmax,4))
!
  call read_atmpairidx(iunit,ntlnks,tmpvec(1:vmax,1:2),vmax,atrue,afalse,afalse)
  close(unit=iunit)
!
! discard definitely illegal ones
  k = ntlnks
  ntlnks = 0
  do i=1,k
    if ((tmpvec(i,1).ge.1).AND.(tmpvec(i,1).le.vmax).AND.(tmpvec(i,2).ge.1).AND.(tmpvec(i,2).le.vmax)) then
      ntlnks = ntlnks + 1
      tmpvec(ntlnks,3:4) = tmpvec(i,1:2)
    end if
  end do
!
  if (ntlnks.lt.k) then
    write(ilog,*) 'Warning. ',k-ntlnks,' trajectory links in ',tlnkfilen(t1:t2),' have been discarded (out of range).'
  end if
!
  if (ntlnks.le.0) then
    deallocate(tmpvec)
    write(ilog,*) 'Warning. Set of user-requested trajectory links to add for network-related analyses is empty.'
    return
  end if
!
  allocate(trlnkslst(ntlnks,2))
  trlnkslst(1:ntlnks,1:2) = tmpvec(1:ntlnks,3:4)
  deallocate(tmpvec)
!
end
!
!----------------------------------------------------------------------------
!
subroutine read_trajidxfile()
!
  use iounit
  use atoms
  use system
  use pdb
!
  implicit none
!
  integer i,iunit,freeunit,t1,t2,buff
  logical exists,atrue
  integer, ALLOCATABLE:: tmpvec(:)
!
  buff = pdbeffn
!
  call strlims(trajidxfile,t1,t2)
  inquire(file=trajidxfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &trajidxfile(t1:t2),') for FMCSC_TRAJIDXFILE. Ignoring request ...'
    use_trajidx = .false.
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  atrue = .true.
  open (unit=iunit,file=trajidxfile(t1:t2),status='old')
!
  allocate(tmpvec(n))
!
! get the unique, sorted set of atom indices
  pdbeffn = 0
  call read_atmidx(iunit,pdbeffn,tmpvec,n,atrue,atrue)
  close(unit=iunit)
  if (pdbeffn.lt.1) then
    write(ilog,*) 'Set of atom indices for trajectory output has zero size. &
 &Adjust FMCSC_XYZOUT to disable trajectory output completely (request ignored).'
    use_trajidx = .false.
    pdbeffn = buff
    deallocate(tmpvec)
    return
  end if
!
  allocate(pdboutlst(pdbeffn))
  allocate(pdboutlog(n))
  pdboutlst(1:pdbeffn) = tmpvec(1:pdbeffn)
  pdboutlog(:) = .false.
  do i=1,pdbeffn
    pdboutlog(pdboutlst(i)) = .true.
  end do
  deallocate(tmpvec)
!
end
!
!----------------------------------------------------------------------------
!
subroutine read_framesfile()
!
  use iounit
  use atoms
  use system
  use pdb
  use interfaces
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,t1,t2,buff,lstsz,iomess,ii,kk,nlst,mn
  logical exists,notdone,skipsort
  character(MAXSTRLEN) str2,str3
  integer, ALLOCATABLE:: tmpvec(:,:)
  RTYPE, ALLOCATABLE:: tmpvec2(:),tmpvec3(:)
  RTYPE fpv
!
  buff = nsim
!
  call strlims(frameidxfile,t1,t2)
  inquire(file=frameidxfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &frameidxfile(t1:t2),') for FMCSC_FRAMESFILE. Ignoring request ...'
    select_frames = .false.
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=frameidxfile(t1:t2),status='old')
!
  lstsz = 0
  iomess = 0
  do while (iomess.eq.0)
    read(iunit,79,iostat=iomess) str2
    lstsz = lstsz + 1
  end do
  rewind(unit=iunit)
  allocate(tmpvec(lstsz,5))
  allocate(tmpvec2(lstsz))
!
  nlst = 0
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_framesfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((mn.gt.nsim).OR.(mn.le.0)) then
      write(ilog,*) 'Requested frame number ',mn,' for trajectory frame list,&
 & which is invalid or inconsistent with the provided simulation length. Ignored ...'
      cycle
    end if
    nlst = nlst + 1
    if (nlst.gt.nsim) then
      if (.NOT.((pdb_analyze.EQV..true.).AND.((pdb_fileformat.eq.5).OR.(pdb_fileformat.eq.2)).AND.&
 &                                           (use_frame_weights.EQV..false.))) then
        write(ilog,*) 'Fatal. Frame list provided to read_framesfile(...) is longer &
 &than the specified simulation length. This is only allowed for frame lists processed "as is" (&
 &using input from individual PDB files or a NetCDF trajectory and not using frame weights).'
        call fexit()
      end if
    end if
!
    tmpvec(nlst,3) = nlst
    tmpvec(nlst,1) = mn
    kk = 1
    call extract_str(str2(ii:MAXSTRLEN),str3,kk)
    if (kk.eq.1) then
!     assume weight of unity
      tmpvec2(nlst) = 1.0
    else
      read(str3,*,iostat=iomess) fpv
      if ((iomess.eq.IOSTAT_END).OR.(iomess.eq.1)) then
        tmpvec2(nlst) = 1.0
      else
        tmpvec2(nlst) = fpv
        use_frame_weights = .true.
      end if
    end if
!
  end do
  if (nlst.lt.1) then
    write(ilog,*) 'Fatal. Set of frame indices for PDB analysis mode has zero size (in ',&
 &frameidxfile(t1:t2),').'
    select_frames = .false.
    deallocate(tmpvec)
    deallocate(tmpvec2)
    call fexit()
  end if
!
  close(unit=iunit)
!
  skipsort = .false.
  if ((pdb_analyze.EQV..true.).AND.((pdb_fileformat.eq.5).OR.(pdb_fileformat.eq.2)).AND.(use_frame_weights.EQV..false.)) then
    write(ilog,*) 'Skipping the sorting of input for FMCSC_FRAMESFILE due to compatible file format. &
 &This means that duplicates are allowed, snapshots are processed exactly in input order, and all step &
 &numbering from now on refers to the line numbers in ',frameidxfile(t1:t2),'.' 
    skipsort = .true.
    tmpvec(1:nlst,4) = tmpvec(1:nlst,1)
  else
    if ((pdb_analyze.EQV..true.).AND.((pdb_fileformat.eq.5).OR.(pdb_fileformat.eq.2))) then
      write(ilog,*) 'Warning. Despite using a file format compatible with processing trajectory frames in arbitrary order, &
 &the input for FMCSC_FRAMESFILE (',frameidxfile(t1:t2),') is sorted, and duplicates are removed. This is because of the use &
 &of floating point weights.'
    end if
!   sort
    ii = 1
    kk = nlst
    mn = nlst
    call merge_sort(ldim=mn,up=.TRUE.,list=tmpvec(:,1),olist=tmpvec(:,4),ilo=ii,ihi=kk,idxmap=tmpvec(:,3),olist2=tmpvec(:,2))
!   in-place sorting of associated array (tmpvec2) requires use of reverse index
    do ii=1,nlst
      tmpvec(tmpvec(ii,2),5) = ii
    end do
    do ii=1,nlst
      if (ii.eq.tmpvec(ii,2)) cycle
      fpv = tmpvec2(ii)
      tmpvec2(ii) = tmpvec2(tmpvec(ii,2))
      tmpvec2(tmpvec(ii,2)) = fpv
      kk = tmpvec(ii,5)
      mn = tmpvec(ii,2)
      tmpvec(ii,5) = tmpvec(mn,5)
      tmpvec(mn,5) = kk
      tmpvec(kk,2) = tmpvec(ii,2)
    end do

!   remove double entries
    ii = 2
    notdone = .true.
    do while (notdone.EQV..true.)
      if (tmpvec(ii,4).eq.tmpvec(ii-1,4)) then
        write(ilog,*) 'Warning. Removing double entry from index list in read_framesfile(...).'
        tmpvec((ii-1):(nlst-1),4) = tmpvec(ii:nlst,4)
        tmpvec2((ii-1):(nlst-1)) = tmpvec2(ii:nlst)
        nlst = nlst - 1
      else
        ii = ii + 1
      end if
      if (ii.gt.nlst) notdone = .false.
    end do
  end if
!
  buff = nlst
  do lstsz=1,nlst
    if (tmpvec(lstsz,4).gt.nsim) then
      buff = lstsz-1
      exit
    end if
  end do
  if (buff.le.0) then
    write(ilog,*) 'Fatal. Setting for trajectory file length (',nsim,') and selected frames &
 &in file "',frameidxfile(t1:t2),'" have no overlap. Please adjust input files.'
    deallocate(tmpvec)
    deallocate(tmpvec2)
    call fexit()
  end if
  framecnt = buff
  allocate(framelst(framecnt))
  if (use_frame_weights.EQV..true.) then
    allocate(framewts(framecnt))
    framewts(1:framecnt) = tmpvec2(1:framecnt)
  end if
  framelst(1:framecnt) = tmpvec(1:framecnt,4)
!
  if (use_frame_weights.EQV..true.) then
!   eventually create a list sorted by increasing weight
    allocate(tmpvec3(framecnt))
    mn = framecnt
    tmpvec(1:framecnt,1) = framelst(1:framecnt)
    tmpvec2(1:framecnt) = framewts(1:framecnt)
    do ii=1,framecnt
      tmpvec(ii,3) = ii
    end do
    ii = 1
    kk = framecnt
    mn = framecnt
    call merge_sort(ldim=mn,up=.TRUE.,list=tmpvec2(1:mn),olist=tmpvec3(1:mn),ilo=ii,ihi=kk,idxmap=tmpvec(1:mn,3),&
   &olist2=tmpvec(1:mn,2))
    do ii=1,framecnt
      tmpvec(tmpvec(ii,2),5) = ii
    end do
    do ii=1,framecnt
      if (ii.eq.tmpvec(ii,2)) cycle
      iomess = tmpvec(ii,1)
      tmpvec(ii,1) = tmpvec(tmpvec(ii,2),1)
      tmpvec(tmpvec(ii,2),1) = iomess
      kk = tmpvec(ii,5)
      mn = tmpvec(ii,2)
      tmpvec(ii,5) = tmpvec(mn,5)
      tmpvec(mn,5) = kk
      tmpvec(kk,2) = tmpvec(ii,2)
    end do
!
    allocate(framelst2(framecnt))
    allocate(framewts2(framecnt))
    framelst2(1:framecnt) = tmpvec(1:framecnt,1)
    framewts2(1:framecnt) = tmpvec3(1:framecnt)
!
    deallocate(tmpvec3)
  end if
!
  deallocate(tmpvec)
  deallocate(tmpvec2)
  if ((skipsort.EQV..false.).AND.(nsim.gt.framelst(framecnt))) then
    write(ilog,*) 'Warning. Adjusting run length to match last eligible frame in provided set (',&
 &framelst(framecnt),').'
    curframe = -nsim ! hi-jack for backup
    nsim = framelst(framecnt)
  else if ((skipsort.EQV..true.).AND.(framecnt.gt.nsim)) then
    curframe = nsim ! hi-jack for backup
    nsim = framecnt
  else
    curframe = -nsim ! hi-jack for consistency
  end if
!
end
!
!----------------------------------------------------------------------------
!
subroutine read_savreqfile()
!
  use iounit
  use atoms
  use system
  use fos
!
  implicit none
!
  integer iunit,freeunit,t1,t2
  logical exists,atrue
  integer, ALLOCATABLE:: tmpvec(:)
!
  call strlims(savreq%filen,t1,t2)
  inquire(file=savreq%filen(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &savreq%filen(t1:t2),') for FMCSC_SAVATOMFILE. Ignoring request ...'
    savreq%nats = 0
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  atrue = .true.
  open (unit=iunit,file=savreq%filen(t1:t2),status='old')
!
  allocate(tmpvec(n))
!
! get the unique, sorted set of atom indices
  savreq%nats = 0
  call read_atmidx(iunit,savreq%nats,tmpvec,n,atrue,atrue)
  close(unit=iunit)
  if (savreq%nats.lt.1) then
    write(ilog,*) 'Set of atom indices for individual solvent-accessibility histograms has zero size. &
 &Check input in ',savreq%filen(t1:t2),' (request ignored).'
    savreq%nats = 0
    deallocate(tmpvec)
    return
  end if
!
  allocate(savreq%idx(savreq%nats))
  allocate(savreq%hists(savreq%nats,2,100))
  savreq%idx(1:savreq%nats) = tmpvec(1:savreq%nats)
  savreq%hists(:,:,:) = 0.0
!
  deallocate(tmpvec)
!
end
!
!----------------------------------------------------------------------------------------------
!
subroutine read_torfile()
!
  use energies
  use iounit
  use sequen
  use molecule
  use aminos
  use fyoc
  use torsn
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
  use polypep
  use interfaces
  use system, ONLY: ua_model
  use zmatrix, ONLY: izrot,ztor
  use movesets, ONLY: natlst
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,i,j,ii,jj,k,dc,lastin,tl,shf
  character(2) rsnma
  character(3) resname
  character(10) toty,abc
  character(10000) str2,str3
  character(14), ALLOCATABLE:: stringy(:)
  logical exists,badflag
  RTYPE getpuckertor,getztor,tmpv(MAXTORRES)
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  call strlims(torfile,t1,t2)
  inquire(file=torfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &torfile(t1:t2),') for FMCSC_TORFILE. Turning off torsional terms.'
    use_TOR = .false.
    scale_TOR = 0.0
    return
  end if
!
 79   format(a10000)
!
  allocate(par_TOR(nseq,2*MAXTORRES))
  allocate(par_TOR2(nseq))
  par_TOR2(:) = 0
  par_TOR(:,:) = 0.0
!
  iunit = freeunit()
  call strlims(torfile,t1,t2)
  open (unit=iunit,file=torfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
!
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for torsional &
 &input (',torfile(t1:t2),'). Turning off torsional bias term.'
    use_TOR = .false.
    scale_TOR = 0.0
    close(unit=iunit)
    deallocate(par_TOR2)
    deallocate(par_TOR)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing torsional &
 &input (got: ',torfile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_str(str2,abc,ii)
  if ((abc(1:1).eq.'G').OR.(abc(1:1).eq.'g')) then
    read(iunit,79) str2
    ii = 1
    par_TOR2(1) = -1
    call extract_int(str2,par_TOR2(1),ii)
    str3 = str2(ii:10000)
!
!   first harmonic + global
!
    if (par_TOR2(1).eq.1) then
      read(str3,*,iostat=iomess) (par_TOR(1,i),i=1,2*MAXTORRES)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for torsio&
 &nal input (',torfile(t1:t2),'). Turning off torsional bias term.'
        use_TOR = .false.
        scale_TOR = 0.0
        close(unit=iunit)
        deallocate(par_TOR2)
        deallocate(par_TOR)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing torsio&
 &nal input (got: ',torfile(t1:t2),').'
        call fexit()
      end if
!
      do i=1,MAXTORRES
        if ((par_TOR(1,i).gt.180.0).OR.(par_TOR(1,i).lt.-180.0)) then
          par_TOR2(1) = 0
          scale_TOR = 0.0
          use_TOR = .false.
          write(ilog,*) 'Bad minimum position for torsional bias poten&
 &tial (',par_TOR(1,i),'). Turning bias term off.'
        end if
      end do
      do i=MAXTORRES+1,2*MAXTORRES
        if (par_TOR(1,i).lt.0.0) then
          par_TOR(1,i) = -par_TOR(1,i)
          write(ilog,*) 'Warning. Negative force const. for harmonic t&
 &orsional bias potential. Assuming ',par_TOR(1,i),'.'
        end if
      end do
      do i=2,nseq
        par_TOR2(i) = 1
        par_TOR(i,:) = par_TOR(1,:)
      end do 
!
!   now Gaussian and global
!
    else if (par_TOR2(1).eq.2) then
      read(str3,*,iostat=iomess)  (par_TOR(1,i),i=1,2*MAXTORRES)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for torsio&
 &nal input (',torfile(t1:t2),'). Turning off torsional bias term.'
        use_TOR = .false.
        scale_TOR = 0.0
        close(unit=iunit)
        deallocate(par_TOR2)
        deallocate(par_TOR)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing torsio&
 &nal input (got: ',torfile(t1:t2),').'
        call fexit()
      end if
!
      do i=1,MAXTORRES
        if ((par_TOR(1,i).gt.180.0).OR.(par_TOR(1,i).lt.-180.0)) then
          par_TOR2(1) = 0
          scale_TOR = 0.0
          use_TOR = .false.
          write(ilog,*) 'Bad minimum position for torsional bias poten&
 &tial (',par_TOR(1,i),'). Turning bias term off.'
        end if
      end do
      do i=MAXTORRES+1,2*MAXTORRES
        if (par_TOR(1,i).lt.0.0) then
          par_TOR(1,i) = -par_TOR(1,i)
          write(ilog,*) 'Warning. Negative std. dev. for Gaussian tors&
 &ional bias potential. Assuming ',par_TOR(1,i),'.'
        end if
        if (par_TOR(1,i).gt.0.0) then
          par_TOR(1,i) = 1.0/(2.0*par_TOR(1,i)**2)
        end if
      end do
      do i=2,nseq
        par_TOR2(i) = 2
        par_TOR(i,:) = par_TOR(1,:)
      end do
!
!   now harmonic and each residue according to structural input
!
    else if (par_TOR2(1).eq.3) then
!     this mode has no minimum position
      read(str3,*,iostat=iomess) (par_TOR(1,i),i=MAXTORRES+1,2*MAXTORRES)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for torsio&
 &nal input (',torfile(t1:t2),'). Turning off torsional bias term.'
        use_TOR = .false.
        scale_TOR = 0.0
        close(unit=iunit)
        deallocate(par_TOR2)
        deallocate(par_TOR)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing torsio&
 &nal input (got: ',torfile(t1:t2),').'
        call fexit()
      end if
!
      do i=MAXTORRES+1,2*MAXTORRES
        if (par_TOR(1,i).lt.0.0) then
          par_TOR(1,i) = -par_TOR(1,i)
          write(ilog,*) 'Warning. Negative force const. for harmonic i&
 &nit. tors. bias potential. Assuming ',par_TOR(1,i),'.'
        end if
      end do
      do i=2,nseq
        par_TOR2(i) = 3
        par_TOR(i,MAXTORRES+1:2*MAXTORRES) = par_TOR(1,MAXTORRES+1:2*MAXTORRES)
      end do
!
!   now Gaussian and each residue according to structural input
!
    else if (par_TOR2(1).eq.4) then
!     this mode has no minimum position
      read(str3,*,iostat=iomess) (par_TOR(1,i),i=MAXTORRES+1,2*MAXTORRES)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for torsio&
 &nal input (',torfile(t1:t2),'). Turning off torsional bias term.'
        use_TOR = .false.
        scale_TOR = 0.0
        close(unit=iunit)
        deallocate(par_TOR2)
        deallocate(par_TOR)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing torsio&
 &nal input (got: ',torfile(t1:t2),').'
        call fexit()
      end if
!
      do i=MAXTORRES+1,2*MAXTORRES
        if (par_TOR(1,i).lt.0.0) then
          par_TOR(1,i) = -par_TOR(1,i)
          write(ilog,*) 'Warning. Negative std. dev. for Gaussian init&
 &. tors. bias potential. Assuming ',par_TOR(1,i),'.'
        end if
        if (par_TOR(1,i).gt.0.0) then
          par_TOR(1,i) = 1.0/(2.0*par_TOR(1,i)**2)
        end if
      end do
      do i=2,nseq
        par_TOR2(i) = 4
        par_TOR(i,MAXTORRES+1:2*MAXTORRES) = par_TOR(1,MAXTORRES+1:2*MAXTORRES)
      end do
    else
      write(ilog,*) 'Unknown torsional mode while reading&
 & file ',torfile(1:t2),'. Fatal exit.'
      close(unit=iunit)
      call fexit()
    end if
!
  else if ((abc(1:1).eq.'R').OR.(abc(1:1).eq.'r')) then
    do while (.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing torsio&
 &nal input (got: ',torfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      k = -1
      call extract_int(str2,k,ii)
!
      if (ii.eq.1) exit
      if ((k.gt.nseq).OR.(k.le.0)) then
        write(ilog,*) 'Warning. Requested torsional bias potential for re&
 &sidue ',k,', which does not exist. Skipping ...'
        cycle
      end if
!
      par_TOR2(k) = -1
      call extract_int(str2,par_TOR2(k),ii)
      str3 = str2(ii:10000)
      if (par_TOR2(k).eq.0) then
        cycle
!
!     now harmonic for specific residues
!
      else if (par_TOR2(k).eq.1) then
        read(str3,*,iostat=iomess) (par_TOR(k,i),i=1,2*MAXTORRES)
        if (iomess.eq.IOSTAT_END) then
          write(ilog,*) 'Warning. Got empty/incomplete file for tors&
 &ional input (',torfile(t1:t2),'). Turning off torsional bias term.'
          use_TOR = .false.
          scale_TOR = 0.0
          close(unit=iunit)
          deallocate(par_TOR2)
          deallocate(par_TOR)
          return
        else if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing tors&
 &ional input (got: ',torfile(t1:t2),').'
          call fexit()
        end if
!
        badflag = .false.
        do i=1,MAXTORRES
          if ((par_TOR(k,i).gt.180.0).OR.(par_TOR(k,i).lt.-180.0)) then
            par_TOR2(k) = 0
            write(ilog,*) 'Bad minimum position for torsional bias poten&
 &tial for residue ',k,'  (',par_TOR(k,i),'). Turning off.'
            badflag = .true.
          end if
        end do
        if (badflag.EQV..true.) cycle
!
        do i=MAXTORRES+1,2*MAXTORRES
          if (par_TOR(k,i).lt.0.0) then
            par_TOR(k,i) = -par_TOR(k,i)
            write(ilog,*) 'Warning. Negative force const. for torsiona&
 &l pot. for residue ',k,'. Assuming ',par_TOR(k,i),'.'
          end if
        end do
!
!     Gaussian for specific residues
!
      else if (par_TOR2(k).eq.2) then
        read(str3,*,iostat=iomess) (par_TOR(k,i),i=1,2*MAXTORRES)
        if (iomess.eq.IOSTAT_END) then
          write(ilog,*) 'Warning. Got empty/incomplete file for tors&
 &ional input (',torfile(t1:t2),'). Turning off torsional bias term.'
          use_TOR = .false.
          scale_TOR = 0.0
          close(unit=iunit)
          deallocate(par_TOR2)
          deallocate(par_TOR)
          return
        else if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing tors&
 &ional input (got: ',torfile(t1:t2),').'
          call fexit()
        end if
!
        badflag = .false.
        do i=1,MAXTORRES
          if ((par_TOR(k,i).gt.180.0).OR.(par_TOR(k,i).lt.-180.0)) then
            par_TOR2(k) = 0
            write(ilog,*) 'Bad minimum position for torsion&
 &al pot. for residue ',k,' (',par_TOR(k,i),').'
            write(ilog,*) 'Turning bias term off.'
            badflag = .true.
          end if
        end do
        if (badflag.EQV..true.) cycle
!
        do i=MAXTORRES+1,2*MAXTORRES
          if (par_TOR(k,i).lt.0.0) then
            par_TOR(k,i) = -par_TOR(k,i)
            write(ilog,*) 'Warning. Negative std. dev. for init. tors.&
 & pot. for residue ',k,'. Assuming ',par_TOR(k,i),'.'
          end if
          if (par_TOR(k,i).gt.0.0) then
            par_TOR(k,i) = 1.0/(2.0*par_TOR(k,i)**2)
          end if
        end do
!
!     harmonic for specific residues with minimum from initial structure
!
      else if (par_TOR2(k).eq.3) then
        read(str3,*,iostat=iomess) (par_TOR(k,i),i=MAXTORRES+1,2*MAXTORRES)
        if (iomess.eq.IOSTAT_END) then
          write(ilog,*) 'Warning. Got empty/incomplete file for tors&
 &ional input (',torfile(t1:t2),'). Turning off torsional bias term.'
          use_TOR = .false.
          scale_TOR = 0.0
          close(unit=iunit)
          deallocate(par_TOR2)
          deallocate(par_TOR)
          return
        else if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing tors&
 &ional input (got: ',torfile(t1:t2),').'
          call fexit()
        end if
!
        do i=MAXTORRES+1,2*MAXTORRES
          if (par_TOR(k,i).lt.0.0) then
            par_TOR(k,i) = -par_TOR(k,i)
            write(ilog,*) 'Warning. Negative force const. for torsiona&
 &l pot. for residue ',k,'. Assuming ',par_TOR(k,i),'.'
          end if
        end do
!
!     and finally residue-specific with Gaussian bias and minimum from initial structure
!
      else if (par_TOR2(k).eq.4) then
        read(str3,*,iostat=iomess) (par_TOR(k,i),i=MAXTORRES+1,2*MAXTORRES)
        if (iomess.eq.IOSTAT_END) then
          write(ilog,*) 'Warning. Got empty/incomplete file for tors&
 &ional input (',torfile(t1:t2),'). Turning off torsional bias term.'
          use_TOR = .false.
          scale_TOR = 0.0
          close(unit=iunit)
          deallocate(par_TOR2)
          deallocate(par_TOR)
          return
        else if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing tors&
 &ional input (got: ',torfile(t1:t2),').'
          call fexit()
        end if
!
        do i=MAXTORRES+1,2*MAXTORRES
          if (par_TOR(k,i).lt.0.0) then
            par_TOR(k,i) = -par_TOR(k,i)
            write(ilog,*) 'Warning. Negative std. dev. for init. tors.&
 & pot. for residue ',k,'. Assuming ',par_TOR(k,i),'.'
          end if
          if (par_TOR(k,i).gt.0.0) then
            par_TOR(k,i) = 1.0/(2.0*par_TOR(k,i)**2)
          end if
        end do
      else
        write(ilog,*) 'Unknown torsional mode while reading file ',&
 &torfile(1:t2),'. Fatal exit.'
        close(unit=iunit)
        call fexit()
      end if
    end do
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',torfile(1:t2),': ',abc,'! Fatal exit.'
    close(unit=iunit)
    call fexit()
  end if
  close(unit=iunit)
!
!
 24   format(i4,' (',a3,'): ',a8)
 25   format(12x,'Phi0: ',f7.2,' Psi0: ',f7.2,' KPhi: ',f7.2,' KPsi: '&
 &,f7.2)
 26   format(12x,'Phi0: ',f7.2,' Psi0: ',f7.2,' DPhi: ',f7.2,' DPsi: '&
 &,f7.2)
 27   format(200(f6.1,',',f6.1,'|'))
 37   format(i10,1x,i1,1x,240(g16.9,' ',g16.9,' '))
 28   format(14x)
 29   format(1000(a14))
!
  j = 0
  do i=1,nseq
    dc = 0
    if (par_TOR2(i).gt.0) then
      if (wline(i).gt.0) then
        dc = dc + 1
        if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) par_TOR(i,1) = omega(i)
      else
        if (nnucs(i).eq.0) par_TOR(i,MAXTORRES+1) = 0.0
      end if
      if (fline(i).gt.0) then
        dc = dc + 1
        if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) par_TOR(i,2) = phi(i)
      else
        if (nnucs(i).eq.0) par_TOR(i,MAXTORRES+2) = 0.0
      end if
      if (yline(i).gt.0) then
        dc = dc + 1
        if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) par_TOR(i,3) = psi(i)
      else
        if (nnucs(i).eq.0) par_TOR(i,MAXTORRES+3) = 0.0
      end if
      if (dc.gt.0) par_TOR(i,MAXTORRES+4:MAXTORRES+6) = 0.0
      if ((dc.gt.0).AND.(nnucs(i).gt.0)) then
        write(ilog,*) 'Fatal. Torsional bias potentials do not currently support polymers which&
 & utilize both fyo-d.o.f.s and nucleic acid d.o.f.s. This is most certainly an omission bug.'
        call fexit()
      end if
      do k=1,nnucs(i) 
        dc = dc + 1
        if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) par_TOR(i,k) = nucs(k,i)
      end do
!     add the sugar
      if (seqpolty(i).eq.'N') then
        dc = dc + 1
        if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) par_TOR(i,dc) = getpuckertor(i,tl)
      end if
      if ((nnucs(i).gt.0).AND.(nnucs(i).lt.5)) par_TOR(i,MAXTORRES+nnucs(i)+2:MAXTORRES+6) = 0.0
      if (dc.gt.6) then
        write(ilog,*) 'Fatal. Torsional bias potentials do not currently support polymers which&
 & have more than six nucleic acid d.o.f.s. This is most certainly an omission bug.'
        call fexit()
      end if
      do k=1,nchi(i)
        dc = dc + 1
        if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) par_TOR(i,6+k) = chi(k,i)
      end do
      k = 6+nchi(i)
!     add crosslink dihedrals
      if (disulf(i).gt.0) then
        if (crosslink(crlk_idx(i))%itstype.le.2) then
          dc = dc + 1
          k = 6+nchi(i)+1
          if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) then
            par_TOR(i,6+nchi(i)+1) = getztor(cai(i),at(i)%sc(2-shf),at(i)%sc(3-shf),at(disulf(i))%sc(3-shf))
          end if
          if (disulf(i).gt.i) then
            dc = dc + 1
            k = 6+nchi(i)+2
            if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) then
              par_TOR(i,6+nchi(i)+2) = getztor(at(i)%sc(2-shf),at(i)%sc(3-shf),at(disulf(i))%sc(3-shf),at(disulf(i))%sc(2-shf))
            end if
          end if
        end if
      end if
!     add nonnative dihedrals in unsupported residues 
      if (seqtyp(i).eq.26) then
        if ((seqpolty(i).ne.'P').AND.(seqpolty(i).ne.'N')) then
         jj  = 0
        else
         jj  = 6
        end if
        do k=at(i)%bb(1),at(i)%bb(1)+at(i)%na-1
          if (izrot(k)%alsz.gt.0) then
            if (natlst%nr.gt.0) then
              call binary_search(natlst%nr,natlst%idx(1:natlst%nr),k,ii)
              if (natlst%idx(max(1,min(natlst%nr,ii))).eq.k) cycle
            end if
            if ((k.eq.fline(i)).OR.(k.eq.yline(i)).OR.(k.eq.fline2(i)).OR.(k.eq.yline2(i))) cycle ! these are so annoying
            dc = dc + 1
            jj = jj + 1
            if (jj.gt.12) exit
            if ((par_TOR2(i).eq.3).OR.(par_TOR2(i).eq.4)) then
              par_TOR(i,jj) = ztor(k)
            end if
          end if
        end do
      end if
!
      if (k.lt.MAXTORRES) par_TOR(i,MAXTORRES+k+1:2*MAXTORRES) = 0.0
!
      do k=MAXTORRES+1,2*MAXTORRES
        if (par_TOR(i,k).gt.0.0) then
          j = j + 1
          exit
        end if
      end do
    end if
    jj = 0
    do k=MAXTORRES+1,2*MAXTORRES
      if (par_TOR(i,k).gt.0.0) then
        jj = jj + 1
      end if
    end do
    if (jj.le.0) par_TOR2(i) = 0
  end do
  if (j.eq.0) then
    write(ilog,*) 
    write(ilog,*) 'WARNING: No torsional bias potentials (FMCSC_SC_TOR) requested (?). Disabling term.'
    write(ilog,*)
    use_TOR = .false.
    scale_TOR = 0.0
    return
  end if
!
  allocate(stringy(MAXTORRES))
!
  if (tor_report.EQV..true.) then
    tl = 2
    write(ilog,*)
    write(ilog,*) '--- Summary of Torsional Bias Terms in the System ---'
    write(ilog,*)
    do i=1,nseq
      stringy(:) = '              '
      resname = amino(seqtyp(i))
      lastin = 0
      if (wline(i).gt.0) then
        stringy(1) = '     OME      '
        lastin = 1
      end if
      if (fline(i).gt.0) then
        stringy(2) = '     PHI      '
        lastin = 2
      end if
      if (yline(i).gt.0) then
        stringy(3) = '     PSI      '
        lastin = 3
      end if 
      do j=1,nnucs(i)
        call int2str(j,rsnma,tl)
        stringy(j) = '    NUC '//rsnma(1:tl)//'    '
        lastin = j
      end do
      if (seqpolty(i).eq.'N') stringy(6) = '    SUGAR     '
      do j=1,nchi(i)
        call int2str(j,rsnma,tl)
        stringy(6+j) = '    CHI '//rsnma(1:tl)//'    '
        lastin = 6+j
      end do
      if (disulf(i).gt.0) then
        if (crosslink(crlk_idx(i))%itstype.le.2) then
          lastin = lastin + 1
          stringy(lastin) = '    CHI 02    '
          if (disulf(i).gt.i) then
            lastin = lastin + 1
            stringy(lastin) = '    CSSC      '
          end if
        end if
      end if
      if (seqtyp(i).eq.26) then
        if ((seqpolty(i).ne.'P').AND.(seqpolty(i).ne.'N')) then
          lastin = 0
        else
          lastin = 6
        end if
        jj = 0 
        do k=at(i)%bb(1),at(i)%bb(1)+at(i)%na-1
          if (izrot(k)%alsz.gt.0) then
            if (natlst%idx(max(1,min(natlst%nr,ii))).eq.k) cycle
            if ((k.eq.fline(i)).OR.(k.eq.yline(i)).OR.(k.eq.fline2(i)).OR.(k.eq.yline2(i))) cycle ! these are so annoying
            lastin = lastin + 1
            jj = jj + 1
            call int2str(jj,rsnma,tl)
            stringy(lastin) = '    UNS '//rsnma(1:tl)//'    '
            if (lastin.eq.12) exit
          end if
        end do
      end if
!
      if (lastin.eq.0) cycle
      if (par_TOR2(i).eq.0) then
        toty = 'NONE      '
        write(ilog,24) i,resname,toty
      else if ((par_TOR2(i).eq.1).OR.(par_TOR2(i).eq.3)) then
        if (par_TOR2(i).eq.1) toty = 'HARMONIC  '
        if (par_TOR2(i).eq.3) toty = 'HARMONIC_I'
        write(ilog,24) i,resname,toty
        write(ilog,29) (stringy(k),k=1,lastin)
        do k=1,lastin
          if (par_TOR(i,MAXTORRES+k).gt.0.0) then 
            write(ilog,27,advance="no") par_TOR(i,k),par_TOR(i,MAXTORRES+k)
          else
            write(ilog,28,advance="no")
          end if
        end do
        write(ilog,*)
      else if ((par_TOR2(i).eq.2).OR.(par_TOR2(i).eq.4)) then
        if (par_TOR2(i).eq.2) toty = 'GAUSSIAN  '
        if (par_TOR2(i).eq.4) toty = 'GAUSSIAN_I'
        write(ilog,24) i,resname,toty
        write(ilog,29) (stringy(k),k=1,lastin)
        do k=1,lastin
          if (par_TOR(i,MAXTORRES+k).gt.0.0) then 
            write(ilog,27,advance="no") par_TOR(i,k),1.0/sqrt(2.0*par_TOR(i,MAXTORRES+k))
          else
            write(ilog,28,advance="no")
          end if
        end do 
        write(ilog,*)
!        write(ilog,26) par_TOR(i,1),par_TOR(i,2),&
! &       1.0/sqrt(2.0*par_TOR(i,3)),1.0/sqrt(2.0*par_TOR(i,4))
      else
        toty = 'UNKNOWN   '
        write(ilog,24) i,resname,toty
      end if
    end do
    write(ilog,*)
!
#ifdef ENABLE_MPI
    j = 3
    call int2str(myrank,resname,j)
    str3 ='N_'//resname(1:j)//'_SAMPLE_TORFILE.dat'
#else
    str3 = 'SAMPLE_TORFILE.dat'
#endif
    call strlims(str3,t1,t2)
    inquire(file=str3(t1:t2),exist=exists)
    if(exists) then
      iunit = freeunit()
      open(unit=iunit,file=str3(t1:t2),status='old')
      close(unit=iunit,status='delete')
    end if
    iunit=freeunit()
    open(unit=iunit,file=str3(t1:t2),status='new')
!
    write(iunit,*) 'R'
    do i=1,nseq
      if ((par_TOR2(i).eq.1).OR.(par_TOR2(i).eq.3)) then
        write(iunit,37) i,1,par_TOR(i,1:MAXTORRES),par_TOR(i,(MAXTORRES+1):(2*MAXTORRES))
      else if ((par_TOR2(i).eq.2).OR.(par_TOR2(i).eq.4)) then
        do j=1,MAXTORRES
          if (par_TOR(i,j+MAXTORRES).gt.0.0) then
            tmpv(j) = 1.0/sqrt(2.0*par_TOR(i,j+MAXTORRES))
          else
            tmpv(j) = 0.0 
          end if
        end do
        write(iunit,37) i,2,par_TOR(i,1:MAXTORRES),tmpv(1:MAXTORRES)
      end if      
    end do
    close(unit=iunit)
  end if
!
  deallocate(stringy)
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_polfile()
!
  use energies
  use iounit
  use molecule
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer NRPOLYPARAMS
  parameter (NRPOLYPARAMS=4)
!
  integer iunit,freeunit,iomess,t1,t2,ii,i,j,k,polyp2(nmol),imol,mn
  integer mt
  character(10) toty,abc
  character(MAXSTRLEN) str2,str3
  logical exists
  RTYPE polyp(nmol,NRPOLYPARAMS)
!
! the only thing we already checked for is the existence of the actual input file
! UNLESS mr.user elected to use the default -> check again
  call strlims(polfile,t1,t2)
  inquire(file=polfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &polfile(t1:t2),') for FMCSC_POLYFILE. Turning off polymeric terms.'
    use_POLY = .false.
    scale_POLY = 0.0
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  allocate(par_POLY(nmol,4))
  allocate(par_POLY2(nmol))
  do i=1,nmol
    par_POLY2(i) = 0        !mode of polymeric pot.
    do j=1,4
      par_POLY(i,j) = 0.0   !specific parameters (like harmonic spring constant)
    end do
  end do
!
  iunit = freeunit()
  open (unit=iunit,file=polfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for polymeric &
 &input (',polfile(t1:t2),'). Turning off polymeric term.'
    use_POLY = .false.
    scale_POLY = 0.0
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing polymeric &
 &input (got: ',polfile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_str(str2,abc,ii)
!
! first molecule-based
!
  if ((abc(1:1).eq.'M').OR.(abc(1:1).eq.'m')) then
    polbiasmode = 1
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing polyme&
 &ric input (got: ',polfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      mn = -1
      call extract_int(str2,mn,ii)
      if ((mn.gt.nmol).OR.(mn.le.0)) then
        write(ilog,*) 'Requested polymeric potential for molecule ',&
 &mn,', which does not exist. Skipping ...'
        cycle
      end if
      if (ii.eq.1) exit
!
      par_POLY2(mn) = 1
      str3 = str2(ii:MAXSTRLEN)
!     for both std. harmonic and exponential torsional potentials there are four parameters:
!     target-Rg, target-delta, spring-Rg, spring-delta
      read(str3,*,iostat=iomess) par_POLY(mn,1),par_POLY(mn,2),&
 &                  par_POLY(mn,3),par_POLY(mn,4)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for polyme&
 &ric input (',polfile(t1:t2),'). Turning off polymeric term.'
        use_POLY = .false.
        scale_POLY = 0.0
        close(unit=iunit)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing polyme&
 &ric input (got: ',polfile(t1:t2),').'
          call fexit()
        end if
!
      if ((par_POLY(mn,1).lt.0.0).OR.(par_POLY(mn,1).gt.1.0).OR.&
 &     (par_POLY(mn,2).gt.1.0).OR.(par_POLY(mn,2).lt.0.0)) then
        par_POLY2(mn) = 0
        write(ilog,*) 'Bad minimum position for polymeric bias poten&
 &tial (',par_POLY(mn,1),',',par_POLY(mn,2),').'
        write(ilog,*) 'Turning bias term off for mol. #',&
 & mn,'.'
      end if
      if (par_POLY(mn,3).lt.0.0) then
        par_POLY(mn,3) = -par_POLY(mn,3)
        write(ilog,*) 'Warning. Negative force const. for harmonic t&
 & bias potential. Assuming ',par_POLY(mn,3),' for mol. #',mn,'.'
      end if
      if (par_POLY(mn,4).lt.0.0) then
        par_POLY(mn,4) = -par_POLY(mn,4)
        write(ilog,*) 'Warning. Negative force const. for harmonic a&
 &spher. bias pot.. Assuming ',par_POLY(mn,4),' for mol. #',mn,'.'
      end if 
    end do
!
! and now molecule-type-based
!
  else if ((abc(1:1).eq.'T').OR.(abc(1:1).eq.'t')) then
    polbiasmode = 2
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing polyme&
 &ric input (got: ',polfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      mt = -1
      call extract_int(str2,mt,ii)
      if ((mt.gt.nmoltyp).OR.(mt.le.0)) then
        write(ilog,*) 'Fatal. Requested polymeric potential for mole&
 &cule type ',mt,', which does not exist. Skipping ...'
        cycle
      end if
      if (ii.eq.1) exit
      par_POLY2(mt) = 1
      str3 = str2(ii:MAXSTRLEN)
!     for both std. harmonic and exponential torsional potentials there are four parameters:
!     target-Rg, target-delta, spring-Rg, spring-delta
      read(str3,*,iostat=iomess) par_POLY(mt,1),par_POLY(mt,2),&
 &                            par_POLY(mt,3),par_POLY(mt,4)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for polyme&
 &ric input (',polfile(t1:t2),'). Turning off polymeric term.'
        use_POLY = .false.
        scale_POLY = 0.0
        close(unit=iunit)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing polyme&
 &ric input (got: ',polfile(t1:t2),').'
        call fexit()
      end if
!
      if ((par_POLY(mt,1).lt.0.0).OR.(par_POLY(mt,1).gt.1.0).OR.&
 &       (par_POLY(mt,2).gt.1.0).OR.(par_POLY(mt,2).lt.0.0)) then
        par_POLY2(mt) = 0
        write(ilog,*) 'Bad minimum position for polymeric bias poten&
 &tial (',par_POLY(mt,1),',',par_POLY(mt,2),').'
        write(ilog,*) 'Turning bias term off for mol. type #',&
 & mt,'.'
        cycle
      end if
      if (par_POLY(mt,3).lt.0.0) then
        par_POLY(mt,3) = -par_POLY(mt,3)
        write(ilog,*) 'Warning. Negative force const. for harmonic t&
 & bias potential. Assuming ',par_POLY(mt,3),' for mol. type #',mt,'.'
      end if
      if (par_POLY(mt,4).lt.0.0) then
          par_POLY(mt,4) = -par_POLY(mt,4)
        write(ilog,*) 'Warning. Negative force const. for harmonic a&
 &spher. bias pot.. Assuming ',par_POLY(mt,4),' for mol. type #',mt,'.'
      end if
    end do
!
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',polfile(1:t2),': ',abc,'! Fatal exit.'
      close(unit=iunit)
      call fexit()
  end if
  close(unit=iunit)
!
!
 24   format('Mol. #',i4,': ',a8)
 25   format(12x,'T_0: ',f7.2,' Dl_0: ',f7.2,' K_T: ',f7.2,' K_Dl: '&
 &,f7.2)
!
  j = 0
!
! per molecule setup
  if (polbiasmode.eq.1) then
    do imol=1,nmol
      mt = moltypid(imol)
      if (par_POLY2(imol).eq.1) then
        j = j + 1
      end if
    end do
! per molecule-type setup
  else
!   first save mt-dependent info in polyp and then reset par_POLY
    do i=1,nmoltyp
      do k=1,NRPOLYPARAMS
        polyp(i,k) = par_POLY(i,k)
      end do
      polyp2(i) = par_POLY2(i)
    end do
    do i=1,nmol
      do k=1,NRPOLYPARAMS
        par_POLY(i,k) = 0.0
      end do
      par_POLY2(i) = 0
    end do
!   now assign par_POLY based on values in polyp
    do mt=1,nmoltyp
      if (polyp2(mt).eq.1) then
        do imol=1,nmol
          if (moltypid(imol).eq.mt) then
            j = j + 1
            par_POLY2(imol) = 1
            do k=1,NRPOLYPARAMS
              par_POLY(imol,k) = polyp(mt,k)
            end do
          end if
        end do
      end if
    end do
  end if
  if (j.eq.0) then
    write(ilog,*) 
    write(ilog,*) 'WARNING: No polymeric biasing potentials requeste&
 &d (?).'
    write(ilog,*)
    use_POLY = .false.
    scale_POLY = 0.0
    return
  end if
!
  if (poly_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '--- Summary of Polymer Bias Terms in the System -&
 &--'
    write(ilog,*)
    do imol=1,nmol
      if (par_POLY2(imol).eq.0) then
        toty = 'NONE      '
        write(ilog,24) imol,toty
      else if (par_POLY2(imol).eq.1) then
        toty = 'HARMONIC  '
        write(ilog,24) imol,toty
        write(ilog,25) par_POLY(imol,1),par_POLY(imol,2),&
 &                     par_POLY(imol,3),par_POLY(imol,4)
      end if
    end do
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_atpatchfile(mode)
!
  use energies
  use iounit
  use atoms
  use fos
  use params
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer i,iunit,freeunit,iomess,t1,t2,t3,t4,ii,nlst,mn,mode
  character(MAXSTRLEN) str2,str3,fneff,str4,str5
  logical exists
  RTYPE mnrt
  RTYPE, ALLOCATABLE:: newasm(:)
  integer, ALLOCATABLE:: tmpvec(:)
  logical, ALLOCATABLE:: fdvec(:)
!
  nlst = 0
!
  if (mode.eq.1) then 
    call strlims(asmpatchfile,t1,t2)
    fneff = asmpatchfile
    str4 = 'maximum SAV fraction'
  else if (mode.eq.2) then
    call strlims(ardpatchfile,t1,t2)
    fneff = ardpatchfile
    str4 = 'volume reduction factor'
  else if (mode.eq.3) then
    call strlims(masspatchfile,t1,t2)
    fneff = masspatchfile
    str4 = 'mass'
  else if (mode.eq.4) then
    call strlims(radpatchfile,t1,t2)
    fneff = radpatchfile
    str4 = 'radius'
  else
    write(ilog,*) 'Fatal. Called read_atpatchfile with unsupported mode (got ',mode,'). This is a bug.'
    call fexit()
  end if
  call strlims(str4,t3,t4)
  str5(t3:t4) = str4(t3:t4)
  call toupper_first(str5(t3:t4))
!
  inquire(file=fneff(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No patches to atomic ',str4(t3:t4),' applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()

  open (unit=iunit,file=fneff(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for atomic ',str4(t3:t4),' &
 &patch input (',fneff(t1:t2),'). Using default values (from parameters and/or topology).'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing atomic ',str4(t3:t4),' patch &
 &input (',fneff(t1:t2),').'
    call fexit()
  end if
  rewind(unit=iunit)
!
! now just read present entries
!
  allocate(fdvec(n))
  allocate(newasm(n))
  allocate(tmpvec(n))
  fdvec(:) = .false.
!
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_atpatchfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((mn.gt.n).OR.(mn.le.0)) then
      write(ilog,*) 'Requested atom number ',mn,' for applying a patch to its ',str4(t3:t4),'.&
 & It does not exist, however, and the request will be ignored ...'
      cycle
    end if
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    read(str3,*) mnrt
    if ((((mnrt.lt.0.0).OR.(mnrt.gt.1.0)).AND.(mode.le.2)).OR.((mnrt.le.0.9).AND.(mode.eq.3)).OR.&
 &((mnrt.lt.0.0).AND.(mode.eq.4))) then
      write(ilog,*) 'Requested an illegal value for patching ',str4(t3:t4),' of atom ',mn,' &
 &(',mnrt,'). Request will be ignored ...'
      cycle
    end if
!
    nlst = nlst + 1
    if (nlst.gt.n) then
      write(ilog,*) 'Fatal. More entries in file supplied to read_atpatchfile(...) &
 &than there are atoms in the system. This is fatal.'
      deallocate(fdvec)
      deallocate(newasm)
      deallocate(tmpvec)
      call fexit()
    end if
    if (fdvec(mn).EQV..false.) then
      tmpvec(nlst) = mn
      newasm(nlst) = mnrt
      fdvec(mn) = .true.
    else
      nlst = nlst - 1
      write(ilog,*) 'Warning. List provided to read_atpatchfile(...) contains &
 &double entry ( atom # ',mn,'). Ignored (first instance used) ...'
    end if
!
  end do
!
  close(unit=iunit)
!
  if ((nlst.gt.0).AND.(mode.eq.1)) asm_patched = .true.
  if ((nlst.gt.0).AND.(mode.eq.2)) ard_patched = .true.
  if ((nlst.gt.0).AND.(mode.eq.3)) mass_patched = .true.
  if ((nlst.gt.0).AND.(mode.eq.4)) rad_patched = .true.
!
 24   format(' Atom #',i7,' was patched from a value of ',g12.5,' to ',g12.5,a,' for its ',a,'.')
!
  write(ilog,*)
  write(ilog,*) '-- Summary of Applied Patches to Atomic ',str5(t3:t4),' --'
  write(ilog,*)
  do i=1,nlst
    if (mode.eq.1) then
      write(ilog,24) tmpvec(i),atsavmaxfr(tmpvec(i)),newasm(i),'',str4(t3:t4)
      atsavmaxfr(tmpvec(i)) = newasm(i)
    else if (mode.eq.2) then
      write(ilog,24) tmpvec(i),atsavred(tmpvec(i)),newasm(i),'',str4(t3:t4)
      atsavred(tmpvec(i)) = newasm(i)
    else if (mode.eq.3) then
      write(ilog,24) tmpvec(i),mass(tmpvec(i)),newasm(i),'g/mol',str4(t3:t4)
      mass(tmpvec(i)) = newasm(i)
    else if (mode.eq.4) then
      write(ilog,24) tmpvec(i),atr(tmpvec(i)),newasm(i),'Angstrom',str4(t3:t4)
      atr(tmpvec(i)) = newasm(i)
    end if
  end do
  write(ilog,*)
  deallocate(fdvec)
  deallocate(newasm)
  deallocate(tmpvec)
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_chargepatchfile(tmpvec,nlst)
!
  use energies
  use iounit
  use atoms
  use params
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,t3,t4,i,ii,nlst,mn
  character(MAXSTRLEN) str2,str3
  logical exists
  RTYPE mnrt
  RTYPE, ALLOCATABLE:: newatq(:)
  integer tmpvec(n)
  logical, ALLOCATABLE:: fdvec(:)
!
  nlst = 0
  if (use_POLAR.EQV..false.) return
!
! the only thing we already checked for is the existence of the actual input file,
! UNLESS mr.user elected to use the default -> check again
  call strlims(cpatchfile,t1,t2)
  inquire(file=cpatchfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No charge patches applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()

  open (unit=iunit,file=cpatchfile(t1:t2),status='old')
  call strlims(paramfile,t3,t4)
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for atomic charge patch &
 &input (',cpatchfile(t1:t2),'). Using biotype-based charges from parameter file ',paramfile(t3:t4),'.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing atomic charge patch &
 &input (',cpatchfile(t1:t2),').'
    call fexit()
  end if
  rewind(unit=iunit)
!
! now just read present entries
!
  allocate(fdvec(n))
  allocate(newatq(n))
  fdvec(:) = .false.
!
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_chargepatchfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((mn.gt.n).OR.(mn.le.0)) then
      write(ilog,*) 'Requested atom number ',&
 &mn,' for applying a charge patch. It does not exist, however, and the request will be ignored ...'
      cycle
    end if
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    read(str3,*) mnrt
!
    nlst = nlst + 1
    if (nlst.gt.n) then
      write(ilog,*) 'Fatal. More entries in file supplied to read_chargepatchfile(...) &
 &than there are atoms in the system. This is fatal.'
      deallocate(fdvec)
      deallocate(newatq)
      call fexit()
    end if
    if (fdvec(mn).EQV..false.) then
      tmpvec(nlst) = mn
      newatq(nlst) = mnrt
      fdvec(mn) = .true.
    else
      nlst = nlst - 1
      write(ilog,*) 'Warning. List provided to read_chargepatchfile(...) contains &
 &double entry ( atom # ',mn,'). Ignored (first instance used) ...'
    end if
!
  end do
!
  close(unit=iunit)
!
  if (nlst.gt.0) charge_patched = .true.
!
 24   format(' Atom #',i7,' was patched from a point charge of ',g12.5,'e to ',g12.5,'e.')
!
  write(ilog,*)
  write(ilog,*) '-- Summary of Applied Charge Patches --'
  write(ilog,*)
  do i=1,nlst
    write(ilog,24) tmpvec(i),atq(tmpvec(i)),newatq(i)
    atq(tmpvec(i)) = newatq(i)
  end do
  write(ilog,*)
  deallocate(fdvec)
  deallocate(newatq)
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_ncpatchfile(rsncp)
!
  use iounit
  use atoms
  use sequen
  use params
  use clusters
  use polypep
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,t3,t4,ii,mn,mnrt,i
  character(MAXSTRLEN) str2,str3
  logical exists,afalse
  type(t_cnblst) rsncp(nseq)
  RTYPE tc,ceps
!
  afalse = .false.
!
  call strlims(ncpatchfile,t1,t2)
  inquire(file=ncpatchfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No patches to net charge flags applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  str2(:) = ' '
  str3(:) = ' '
  iunit = freeunit()

  open (unit=iunit,file=ncpatchfile(t1:t2),status='old')
  call strlims(paramfile,t3,t4)
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for net charge flag patch input (',ncpatchfile(t1:t2),'). &
 & Using default assignments.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing net charge flag patch input (',ncpatchfile(t1:t2),').'
    call fexit()
  end if
  rewind(unit=iunit)
!
! now just read present entries
!
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_ncpatchfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((mn.gt.nseq).OR.(mn.le.0)) then
      write(ilog,*) 'Requested residue number ',mn,' for applying a patch to its net charge group parsing.&
 & It does not exist, however, and the request will be ignored ...'
      cycle
    end if
    if (at(mn)%npol.le.0) then
      write(ilog,*) 'Warning. Requested patch to net charge group status for a residue with &
 &polar atoms. Ignored ...'
      cycle
    end if
    if (rsncp(mn)%alsz.ge.0) then
      write(ilog,*) 'Warning. List provided to read_ncpatchfile(...) contains &
 &double entry (residue # ',mn,'). Ignored (first instance used) ...'
      cycle
    end if
!
    tc = sum(atq(at(mn)%pol(1:at(mn)%npol)))
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    ii = 1
    mnrt = -1
    call extract_int(str3,mnrt,ii)
    if ((mnrt.ne.1).AND.(mnrt.ne.0)) then
      write(ilog,*) 'Encountered uninterpretable entry in net charge patch input (',mnrt,'). First number after residue &
 &index has to be either 0 or 1. Ignoring entire request for residue ',mn,'.'
      cycle
    end if
    rsncp(mn)%alsz = mnrt
    if (mnrt.gt.0) then
      call get_reals_from_str(str3,at(mn)%na,rsncp(mn)%dis(:),ii,rsncp(mn)%nbs,afalse)
      if (abs(tc-sum(rsncp(mn)%dis(1:rsncp(mn)%nbs))).gt.(rsncp(mn)%nbs*(10.0**(-precision(ceps))))) then
        write(ilog,*) 'Fatal. The requested set of charge group targets for residue ',mn,' is inconsistent with &
 &its total charge (after eventual patches via FMCSC_CPATCHFILE). Please check input and/or parameters.'
        call fexit()
      end if
    end if
!
  end do
!
  close(unit=iunit)
!
  if (maxval(rsncp(1:nseq)%nbs).gt.0) nc_patched = .true.
!
 26   format(' Residue #',i9,' was patched to now use charge group targets of ',10000(f10.6,1x),'.')
!
  if (nc_patched.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '-- Summary of Applied Patches to Residue-Level Charge Group Targets --'
    write(ilog,*)
    do i=1,nseq
      if (rsncp(i)%nbs.gt.0) then
        write(ilog,26) i,rsncp(i)%dis(1:rsncp(i)%nbs)
      end if
    end do
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_biotpatchfile()
!
  use iounit
  use atoms
  use params
  use sequen
  use molecule
  use pdb
  use aminos
  use polypep
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,t3,t4,i,ii,k,nlst,mn,mnrt,mode,imol
  character(MAXSTRLEN) str2,str3
  character(10) abc
  logical exists,pwarn
  integer, ALLOCATABLE:: newbtyp(:),tmpvec(:),rstlst(:),seqvec(:)
  logical, ALLOCATABLE:: fdvec(:)
!
  nlst = 0
!
  call strlims(biotpatchfile,t1,t2)
  inquire(file=biotpatchfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No patches to biotype (re)assignments applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()

  open (unit=iunit,file=biotpatchfile(t1:t2),status='old')
  call strlims(paramfile,t3,t4)
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for patch input (',biotpatchfile(t1:t2),') for biotype &
 &(re)assignments. No changes applied.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing patch input (',biotpatchfile(t1:t2),'). for biotype (re)assignments.'
    call fexit()
  end if
  ii = 1
  call extract_str(str2,abc,ii)
  mode = 1 ! atom-based
  if ((abc(1:1).eq.'r').OR.(abc(1:1).eq.'R').OR.(abc(1:1).eq.'t').OR.(abc(1:1).eq.'T')) then
    mode = 2
  else if ((abc(1:1).eq.'a').OR.(abc(1:1).eq.'A')) then
    mode = 1
  else
    write(ilog,*) 'Warning. Unknown mode identifier in file for patch input (',biotpatchfile(t1:t2),') for biotype &
 &(re)assignments. Assuming pure index-based reassignment.'
  end if
!
  allocate(fdvec(n))
  allocate(newbtyp(n))
  allocate(tmpvec(n))
  fdvec(:) = .false.
!
! for residue type assignments, first set up an auxiliary array to filter out inappropriate entries
  if (mode.eq.2) then
    allocate(seqvec(nseq))
    if (n_pdbunk.gt.0) then
      allocate(rstlst(4*MAXAMINO+pdb_unkbnd(n_pdbunk,3)))
    else
      allocate(rstlst(4*MAXAMINO))
    end if
    rstlst(:) = 0
    ii = 0
    do i=1,nseq
      imol = molofrs(i)
      if (seqtyp(i).eq.26) then
        ii = ii + 1
        if (pdb_unkbnd(ii,3).ne.i) then
          write(ilog,*) 'Fatal. Setup of unsupported residues is inconsistent in read_biotpatchfile(...). This is a bug.'
          call fexit()
        end if
        if (pdb_unkbnd(ii,4).eq.0) then
          seqvec(i) = MAXAMINO + pdb_unkbnd(ii,3)
        else
          seqvec(i) = MAXAMINO + pdb_unkbnd(ii,4)
        end if
      else
        if ((i.eq.rsmol(imol,1)).AND.(i.eq.rsmol(imol,2))) then
          seqvec(i) = 4*seqtyp(i)
        else if (i.eq.rsmol(imol,2)) then
          seqvec(i) = 4*seqtyp(i) - 1
        else if (i.eq.rsmol(imol,1)) then
          seqvec(i) = 4*seqtyp(i) - 2
        else
          seqvec(i) = 4*seqtyp(i) - 3
        end if
      end if
      if (rstlst(seqvec(i)).eq.0) rstlst(seqvec(i)) = i   
    end do
  end if
!
  pwarn = .false.
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file input in read_biotpatchfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
!
    if (ii.eq.1) cycle
!
    if ((mn.gt.n).OR.(mn.le.0)) then
      write(ilog,*) 'Requested atom number ',mn,' for applying a patch to biotype assignment.&
 & It does not exist, however, and the request will be ignored ...'
      cycle
    end if
!
    if (mode.eq.2) then
      if (atmres(mn).ne.rstlst(seqvec(atmres(mn)))) then
        write(ilog,*) 'Requested atom number ',mn,' for applying a patch to residue type-based biotype assignment.&
 & However, this atom does not correspond to the first instance of a residue of this type, and the request will be ignored ...'
        cycle
      end if
    end if
!
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    read(str3,*) mnrt
    if ((mnrt.gt.n_biotyp).OR.(mnrt.le.0)) then
      write(ilog,*) 'Requested biotype number ',mnrt,' for an atom to be patched, but the request is illegal &
 &(corrupt entry or lack of support in parameter file). Request will be ignored ...'
      cycle
    end if
    if (mnrt.eq.b_type(mn)) then
      if (pwarn.EQV..false.) then
        pwarn = .true.
        write(ilog,*) 'A requested patch for biotype reassignment is redundant, since the atom &
 &already has been assigned that type (in ',paramfile(t3:t4),'). Request will be ignored (further warnings skipped) ...'
      end if
      cycle
    end if
!
    nlst = nlst + 1
    if (nlst.gt.n) then
      write(ilog,*) 'Fatal. More entries in file supplied to read_biotpatchfile(...) &
 &than there are number of atoms in the system. This is fatal.'
      call fexit()
    end if
    if (fdvec(mn).EQV..false.) then
      tmpvec(nlst) = mn
      newbtyp(nlst) = mnrt
      fdvec(mn) = .true.
    else
      nlst = nlst - 1
      write(ilog,*) 'Warning. List provided to read_biotpatchfile(...) contains &
 &double entry (atom # ',mn,'). Ignored (first instance used) ...'
    end if
!
  end do
!
  close(unit=iunit)
!
  if (nlst.gt.0) bt_patched = .true.
!
 24   format(' Atom #',i9,' was patched from biotype ',i6,' to type ',i6,' (in ',a,').')
!
  if (bt_patched.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '-- Summary of Applied Patches to Biotype Assignments --'
    write(ilog,*)
    pwarn = .false.
    do i=1,nlst
      write(ilog,24) tmpvec(i),b_type(tmpvec(i)),newbtyp(i),paramfile(t3:t4)
      if ((pwarn.EQV..false.).AND.(seqtyp(atmres(tmpvec(i))).ne.26)) then
        if ((lj_val(attyp(tmpvec(i))).ne.lj_val(bio_ljtyp(newbtyp(i)))).OR.&
 &          (lj_atnum(attyp(tmpvec(i))).ne.lj_atnum(bio_ljtyp(newbtyp(i)))).OR.&
 &          (lj_weight(attyp(tmpvec(i))).ne.lj_weight(bio_ljtyp(newbtyp(i))))) then
          pwarn = .true.
          write(ilog,*) 'Warning. Atomic parameters such as mass, valence, and proton number may ALL be altered by an applied &
 &patch to biotype assignments.'
        end if
      end if
      b_type(tmpvec(i)) = newbtyp(i)
      attyp(tmpvec(i)) = bio_ljtyp(b_type(tmpvec(i)))
      if (attyp(tmpvec(i)).eq.0) then
        write(ilog,*) 'Fatal. Parameter file lacks support for biotype ',b_type(tmpvec(i)),' selected via patch for atom #',n,').'
        call fexit()
      else
        atnam(tmpvec(i))(1:2) = lj_symbol(attyp(tmpvec(i)))(1:2)
        mass(tmpvec(i)) = lj_weight(attyp(tmpvec(i)))
      end if
      if (mode.eq.2) then
        do k=atmres(tmpvec(i))+1,nseq
          if (seqvec(k).eq.seqvec(atmres(tmpvec(i)))) then
            ii = tmpvec(i) - at(atmres(tmpvec(i)))%bb(1) + at(k)%bb(1)
            write(ilog,24) ii,b_type(ii),newbtyp(i),paramfile(t3:t4)
            if (fdvec(ii).EQV..true.) then
              write(ilog,*) 'Fatal. Residue type-based biotype reassignment is inconsistent in read_biotpatchfile(...). &
 &This is a bug.'
              call fexit()
            end if
            b_type(ii) = newbtyp(i)
            attyp(ii) = bio_ljtyp(b_type(ii))
            atnam(ii)(1:2) = lj_symbol(attyp(ii))(1:2)
            mass(ii) = lj_weight(attyp(ii))
          end if
        end do
      end if
    end do
    write(ilog,*)
  end if
!
  if (mode.eq.2) then
    deallocate(rstlst)
    deallocate(seqvec)
  end if
  deallocate(fdvec)
  deallocate(newbtyp)
  deallocate(tmpvec)
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_ljpatchfile()
!
  use iounit
  use atoms
  use params
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,t3,t4,i,ii,nlst,mn,mnrt
  character(MAXSTRLEN) str2,str3
  logical exists,pwarn
  integer, ALLOCATABLE:: newattyp(:),tmpvec(:)
  logical, ALLOCATABLE:: fdvec(:)
!
  nlst = 0
!
  call strlims(ljpatchfile,t1,t2)
  inquire(file=ljpatchfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No patches to steric and dispersion parameters applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()

  open (unit=iunit,file=ljpatchfile(t1:t2),status='old')
  call strlims(paramfile,t3,t4)
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for patch input (',ljpatchfile(t1:t2),') for steric and dispersion &
 &parameters. Using biotype-based assignments from parameter file ',paramfile(t3:t4),'.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing patch input (',ljpatchfile(t1:t2),'). for steric and dispersion &
 &parameters.'
    call fexit()
  end if
  rewind(unit=iunit)
!
! now just read present entries
!
  allocate(fdvec(n))
  allocate(newattyp(n))
  allocate(tmpvec(n))
  fdvec(:) = .false.
!
  pwarn = .false.
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_ljpatchfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((mn.gt.n).OR.(mn.le.0)) then
      write(ilog,*) 'Requested atom number ',mn,' for applying a patch to steric and dispersion parameters.&
 & It does not exist, however, and the request will be ignored ...'
      cycle
    end if
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    read(str3,*) mnrt
    if ((mnrt.gt.n_ljtyp).OR.(mnrt.le.0)) then
      write(ilog,*) 'Requested LJ type number ',mnrt,' for an atom to be patched, but the request is illegal &
 &(corrupt entry or lack of support in parameter file). Request will be ignored ...'
      cycle
    end if
    if (mnrt.eq.attyp(mn)) then
      if (pwarn.EQV..false.) then
        pwarn = .true.
        write(ilog,*) 'A requested patch for steric and dispersion parameters is redundant, since the atom &
 &already has been assigned that type (in ',paramfile(t3:t4),'). Request will be ignored (further warnings skipped) ...'
      end if
      cycle
    end if
!
    nlst = nlst + 1
    if (nlst.gt.n) then
      write(ilog,*) 'Fatal. More entries in file supplied to read_ljpatchfile(...) &
 &than there are number of atoms in the system. This is fatal.'
      deallocate(fdvec)
      deallocate(newattyp)
      call fexit()
    end if
    if (fdvec(mn).EQV..false.) then
      tmpvec(nlst) = mn
      newattyp(nlst) = mnrt
      fdvec(mn) = .true.
    else
      nlst = nlst - 1
      write(ilog,*) 'Warning. List provided to read_ljpatchfile(...) contains &
 &double entry (atom # ',mn,'). Ignored (first instance used) ...'
    end if
!
  end do
!
  close(unit=iunit)
!
  if (nlst.gt.0) lj_patched = .true.
!
 24   format(' Atom #',i9,' was patched from LJ type ',i6,' to type ',i6,' (in ',a,').')
!
  if (lj_patched.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '-- Summary of Applied Patches to Steric and Dispersion Parameters --'
    write(ilog,*)
    pwarn = .false.
    do i=1,nlst
      write(ilog,24) tmpvec(i),attyp(tmpvec(i)),newattyp(i),paramfile(t3:t4)
      if ((pwarn.EQV..false.).AND.((lj_val(attyp(tmpvec(i))).ne.lj_val(newattyp(i))).OR.&
 & (lj_atnum(attyp(tmpvec(i))).ne.lj_atnum(newattyp(i))).OR.(lj_weight(attyp(tmpvec(i))).ne.lj_weight(newattyp(i))))) then
        pwarn = .true.
        write(ilog,*) 'Warning. Atomic parameters such as mass, valence, and proton number are not modified by an applied &
 &patch of steric and dispersion parameters, even though the underlying types in the parameter file suggest so.'
      end if
      attyp(tmpvec(i)) = newattyp(i)
    end do
    write(ilog,*)
  end if
!
  deallocate(fdvec)
  deallocate(newattyp)
  deallocate(tmpvec)
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_bondedpatchfile(allm)
!
  use energies
  use iounit
  use atoms
  use params
  use inter
  use math
  use aminos
  use sequen
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,t3,t4,i,ii,j,foundit,ix(6),rs,next,nlst
  character(MAXSTRLEN) str2,str3,indix
  logical exists,allm(5)
!
! the only thing we already checked for is the existence of the actual input file,
! UNLESS mr.user elected to use the default -> check again
  call strlims(bpatchfile,t1,t2)
  inquire(file=bpatchfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No patches to bonded parameters applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
 76   format('Replacing bond length potential for atoms ',&
 &i7,' and ',i7,' in residue # ',i6, ' with #',i4,'.')
 77   format('Adding missing bond length potential (#',i4,') for atoms ',&
 &i7,' and ',i7,' in residue # ',i6, '.')
 75   format('Warning. A patch requested for a bond length potential is invalid &
 &(atoms ',i7,' and ',i7,'; type ',i4,') and will be ignored.')
 78   format('--> Biotypes ',i4,' (',a4,') and ',i4,' (',a4,') in residu&
 &e type ',a3,'.')
 66   format('Replacing bond angle potential for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6,' with #',i4,'.')
 67   format('Adding missing bond angle potential (#',i4,') for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6,'.')
 65   format('Warning. A patch requested for a bond angle potential is invalid &
 &(atoms ',i7,', ',i7,', and ,',i7,'; type ',i4,') and will be ignored.')
 68   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), and ',i4,' &
 &(',a4,') in residue type ',a3,'.')
 56   format('Replacing torsional potential for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6, ' with #',i4,'.')
 57   format('Adding missing torsional potential (#',i4,')  for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6,'.')
 55   format('Warning. A patch requested for a torsional potential is invalid &
 &(atoms ',i7,', ',i7,', ',i7,', and ',i7,'; type ',i4,') and will be ignored.')
 58   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), and ',i4,' (',a4,') in residue type ',a3,'.')
 46   format('Replacing improper dihedral angle potential for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6, ' with #',i4,'.')
 47   format('Adding missing improper dihedral angle potential (#',i4,')  for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6,'.')
 45   format('Warning. A patch requested for an improper dihedral angle potential is invalid &
 &(atoms ',i7,', ',i7,', ',i7,', and ',i7,'; type ',i4,') and will be ignored.')
 48   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), and ',i4,' (',a4,') in residue type ',a3,'.')
 36   format('Replacing CMAP potential for atoms ',&
 &i7,',',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6, ' with #',i4,'.')
 37   format('Adding missing CMAP potential (#',i4,') for atoms ',&
 &i7,',',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6,'.')
 35   format('Warning. A patch requested for a CMAP potential is invalid &
 &(atoms ',i7,', ',i7,', ',i7,', ',i7,', and ',i7,'; type ',i4,') and will be ignored.')
 38   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), ',i4,' (',a4,'), and ',i4,' (',a4,') -> residue type of middle atom: ',a3,'.')
!
  iunit = freeunit()

  open (unit=iunit,file=bpatchfile(t1:t2),status='old')
  call strlims(paramfile,t3,t4)
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for bonded parameters patch &
 &input (',bpatchfile(t1:t2),'). Using assignments from parameter file ',paramfile(t3:t4),'.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing bonded parameters patch &
 &input (',bpatchfile(t1:t2),').'
    call fexit()
  end if
  rewind(unit=iunit)
!
  nlst = 0
  if (bonded_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '--------    Patched Terms (',bpatchfile(t1:t2),')    --------'
  else
    write(ilog,*)
    write(ilog,*) '-------- List of Patched Bonded Interactions ( from ',bpatchfile(t1:t2),') --------'
  end if
!
  do while(.true.)
!
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_bondedpatchfile(...).'
      call fexit()
    end if
    next = 1
    call extract_str(str2,indix,next)
    if (next.eq.1) cycle ! blank line
    call toupper(indix)
!
    if (indix(1:11) .eq. 'PATCH_BOND ') then
      if (use_BOND(1).EQV..false.) cycle
      str3 = str2(next:MAXSTRLEN)
      read (str3,*,iostat=iomess) ix(1),ix(2),ix(3)
      if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading bonded parameters patch file&
 & (got: ',bpatchfile(t1:t2),').'
        call fexit()
      end if
      foundit = 0
      do i=1,2
        if ((ix(i).le.0).OR.(ix(i).gt.n)) cycle
        do j=1,nrsbl(atmres(ix(i)))
          if (((iaa(atmres(ix(i)))%bl(j,1).eq.ix(1)).AND.(iaa(atmres(ix(i)))%bl(j,2).eq.ix(2))).OR.&
 &            ((iaa(atmres(ix(i)))%bl(j,1).eq.ix(2)).AND.(iaa(atmres(ix(i)))%bl(j,2).eq.ix(1)))) then
            foundit = j
            rs = atmres(ix(i))
            ix(1:2) = iaa(atmres(ix(i)))%bl(j,1:2)
            exit
          end if
        end do
        if (foundit.gt.0) exit
      end do    
      if ((foundit.le.0).OR.(ix(3).le.0).OR.(ix(3).gt.n_bondtyp)) then
        write(ilog,75) ix(1),ix(2),ix(3)
        cycle
      else
        nlst = nlst + 1
        if (iaa(rs)%typ_bl(j).le.0) then
          allm(1) = .false.
          write(ilog,77) ix(3),ix(1),ix(2),rs
        else
          write(ilog,76) ix(1),ix(2),rs,ix(3)
        end if
        write(ilog,78) b_type(ix(1)),bio_code(b_type(ix(1))),b_type(ix(2)),bio_code(b_type(ix(2))),amino(seqtyp(rs))
        iaa(rs)%typ_bl(j) = bo_typ(ix(3))
        iaa(rs)%par_bl(j,1:MAXBOPAR) = bo_par(ix(3),1:MAXBOPAR)
        if (iaa(rs)%typ_bl(j).eq.3) then ! GROMOS quartic: store squared equ length
          iaa(rs)%par_bl(j,3) = iaa(rs)%par_bl(j,2)*iaa(rs)%par_bl(j,2)
        end if
      end if
    else if (indix(1:12) .eq. 'PATCH_ANGLE ') then
      if (use_BOND(2).EQV..false.) cycle
      str3 = str2(next:MAXSTRLEN)
      read (str3,*,iostat=iomess) ix(1),ix(2),ix(3),ix(4)
      if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading bonded patch file&
 & (got: ',bpatchfile(t1:t2),').'
        call fexit()
      end if
      foundit = 0
      do i=1,3
        if ((ix(i).le.0).OR.(ix(i).gt.n)) cycle
        ii = atmres(ix(i))
        do j=1,nrsba(ii)
          if (((iaa(ii)%ba(j,1).eq.ix(1)).AND.(iaa(ii)%ba(j,2).eq.ix(2)).AND.(iaa(ii)%ba(j,3).eq.ix(3))).OR.&
 &            ((iaa(ii)%ba(j,1).eq.ix(3)).AND.(iaa(ii)%ba(j,2).eq.ix(2)).AND.(iaa(ii)%ba(j,3).eq.ix(1)))) then
            foundit = j
            rs = atmres(ix(i))
             ix(1:3) = iaa(ii)%ba(j,1:3)
            exit
          end if
        end do
        if (foundit.gt.0) exit
      end do    
      if ((foundit.le.0).OR.(ix(4).le.0).OR.(ix(4).gt.n_angltyp)) then
        write(ilog,65) ix(1),ix(2),ix(3),ix(4)
        cycle
      else
        nlst = nlst + 1
        if (iaa(rs)%typ_ba(j).le.0) then
          allm(2) = .false.
          write(ilog,67) ix(4),ix(1),ix(2),ix(3),rs
        else
          write(ilog,66) ix(1),ix(2),ix(3),rs,ix(4)
        end if
        write(ilog,68) b_type(ix(1)),bio_code(b_type(ix(1))),b_type(ix(2)),bio_code(b_type(ix(2))),b_type(ix(3)),&
 &bio_code(b_type(ix(3))),amino(seqtyp(rs))
        iaa(rs)%typ_ba(j) = ba_typ(ix(4))
        iaa(rs)%par_ba(j,1:MAXBAPAR) = ba_par(ix(4),1:MAXBAPAR)
        if ((iaa(rs)%typ_ba(j).ge.1).AND.(iaa(rs)%typ_ba(j).le.2)) then 
          iaa(rs)%par_ba(j,1)=iaa(rs)%par_ba(j,1)/(RADIAN*RADIAN)
        end if
        if (iaa(rs)%typ_ba(j).eq.3) then ! GROMOS cos-harmonic: store cosine of equ angle
          iaa(rs)%par_ba(j,3) = cos(iaa(rs)%par_ba(j,2)/RADIAN)
        end if
      end if
    else if (indix(1:15) .eq. 'PATCH_IMPROPER ') then
      if (use_BOND(3).EQV..false.) cycle
      str3 = str2(next:MAXSTRLEN)
      read (str3,*,iostat=iomess) ix(1),ix(2),ix(3),ix(4),ix(5)
      if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading bonded patch file&
 & (got: ',bpatchfile(t1:t2),').'
        call fexit()
      end if
      foundit = 0
      do i=1,4
        if ((ix(i).le.0).OR.(ix(i).gt.n)) cycle
        ii = atmres(ix(i))
        do j=1,nrsimpt(ii)
          if (((iaa(ii)%impt(j,1).eq.ix(1)).AND.(iaa(ii)%impt(j,2).eq.ix(2)).AND.&
 &             (iaa(ii)%impt(j,3).eq.ix(3)).AND.(iaa(ii)%impt(j,4).eq.ix(4))).OR.&
 &            ((iaa(ii)%impt(j,1).eq.ix(1)).AND.(iaa(ii)%impt(j,2).eq.ix(3)).AND.&
 &             (iaa(ii)%impt(j,3).eq.ix(2)).AND.(iaa(ii)%impt(j,4).eq.ix(4)))) then
            foundit = j
            rs = atmres(ix(i))
            ix(1:4) = iaa(ii)%impt(j,1:4)
            exit
          end if
        end do
        if (foundit.gt.0) exit
      end do    
      if ((foundit.le.0).OR.(ix(5).le.0).OR.(ix(5).gt.n_torstyp)) then
        write(ilog,45) ix(1),ix(2),ix(3),ix(4),ix(5)
        cycle
      else
        nlst = nlst + 1
        if (iaa(rs)%typ_impt(j).le.0) then
          allm(4) = .false.
          write(ilog,47) ix(5),ix(1),ix(2),ix(3),ix(4),rs
        else
          write(ilog,46) ix(1),ix(2),ix(3),ix(4),rs,ix(5)
        end if
        write(ilog,48) b_type(ix(1)),bio_code(b_type(ix(1))),b_type(ix(2)),bio_code(b_type(ix(2))),b_type(ix(3)),&
 &bio_code(b_type(ix(3))),b_type(ix(4)),bio_code(b_type(ix(4))),amino(seqtyp(rs))
        iaa(rs)%typ_impt(j) = di_typ(ix(5))
        iaa(rs)%par_impt(j,1:MAXDIPAR) = di_par(ix(5),1:MAXDIPAR)
!       convert Ryckaert-Bellemans from assumed polymer convention to internal peptide convention
        if (di_typ(ix(5)).eq.3) then
          iaa(rs)%par_impt(j,2) = -iaa(rs)%par_impt(j,2)
          iaa(rs)%par_impt(j,4) = -iaa(rs)%par_impt(j,4)
          iaa(rs)%par_impt(j,6) = -iaa(rs)%par_impt(j,6)
          iaa(rs)%par_impt(j,8) = -iaa(rs)%par_impt(j,8)
!       convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
        else if (di_typ(ix(5)).eq.2) then
          iaa(rs)%par_impt(j,1)=iaa(rs)%par_impt(j,1)/(RADIAN*RADIAN)
        end if
!       re-flag polymer convention Ryckaert-Bellemans
        if (di_typ(ix(5)).eq.3) then
          iaa(rs)%typ_impt(j) = 1
        end if
      end if
    else if (indix(1:14) .eq. 'PATCH_TORSION ') then
      if (use_BOND(4).EQV..false.) cycle
      str3 = str2(next:MAXSTRLEN)
      read (str3,*,iostat=iomess) ix(1),ix(2),ix(3),ix(4),ix(5)
      if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading bonded patch file&
 & (got: ',bpatchfile(t1:t2),').'
        call fexit()
      end if
      foundit = 0
      do i=1,4
        if ((ix(i).le.0).OR.(ix(i).gt.n)) cycle
        ii = atmres(ix(i))
        do j=1,nrsdi(ii)
          if (((iaa(ii)%di(j,1).eq.ix(1)).AND.(iaa(ii)%di(j,2).eq.ix(2)).AND.&
 &             (iaa(ii)%di(j,3).eq.ix(3)).AND.(iaa(ii)%di(j,4).eq.ix(4))).OR.&
 &            ((iaa(ii)%di(j,1).eq.ix(4)).AND.(iaa(ii)%di(j,2).eq.ix(3)).AND.&
 &             (iaa(ii)%di(j,3).eq.ix(2)).AND.(iaa(ii)%di(j,4).eq.ix(1)))) then
            foundit = j
            rs = atmres(ix(i))
            ix(1:4) = iaa(ii)%di(j,1:4)
            exit
          end if
        end do
        if (foundit.gt.0) exit
      end do    
      if ((foundit.le.0).OR.(ix(5).le.0).OR.(ix(5).gt.n_torstyp)) then
        write(ilog,55) ix(1),ix(2),ix(3),ix(4),ix(5)
        cycle
      else
        nlst = nlst + 1
        if (iaa(rs)%typ_di(j).le.0) then
          allm(4) = .false.
          write(ilog,57) ix(5),ix(1),ix(2),ix(3),ix(4),rs
        else
          write(ilog,56) ix(1),ix(2),ix(3),ix(4),rs,ix(5)
        end if
        write(ilog,58) b_type(ix(1)),bio_code(b_type(ix(1))),b_type(ix(2)),bio_code(b_type(ix(2))),b_type(ix(3)),&
 &bio_code(b_type(ix(3))),b_type(ix(4)),bio_code(b_type(ix(4))),amino(seqtyp(rs))
        iaa(rs)%typ_di(j) = di_typ(ix(5))
        iaa(rs)%par_di(j,1:MAXDIPAR) = di_par(ix(5),1:MAXDIPAR)
!       convert Ryckaert-Bellemans from assumed polymer convention to internal peptide convention
        if (di_typ(ix(5)).eq.3) then
          iaa(rs)%par_di(j,2) = -iaa(rs)%par_di(j,2)
          iaa(rs)%par_di(j,4) = -iaa(rs)%par_di(j,4)
          iaa(rs)%par_di(j,6) = -iaa(rs)%par_di(j,6)
          iaa(rs)%par_di(j,8) = -iaa(rs)%par_di(j,8)
!       convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
        else if (di_typ(ix(5)).eq.2) then
          iaa(rs)%par_di(j,1)=iaa(rs)%par_di(j,1)/(RADIAN*RADIAN)
        end if
!       re-flag polymer convention Ryckaert-Bellemans
        if (di_typ(ix(5)).eq.3) then
          iaa(rs)%typ_di(j) = 1
        end if
      end if
    else if (indix(1:11) .eq. 'PATCH_CMAP ') then
      if (use_BOND(5).EQV..false.) cycle
      str3 = str2(next:MAXSTRLEN)
      read (str3,*,iostat=iomess) ix(1),ix(2),ix(3),ix(4),ix(5),ix(6)
      if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while reading bonded patch file&
 & (got: ',bpatchfile(t1:t2),').'
        call fexit()
      end if
      foundit = 0
      do i=1,5
        if ((ix(i).le.0).OR.(ix(i).gt.n)) cycle
        ii = atmres(ix(i))
        do j=1,nrscm(ii)
          if ((iaa(ii)%cm(j,1).eq.ix(1)).AND.(iaa(ii)%cm(j,2).eq.ix(2)).AND.&
 &            (iaa(ii)%cm(j,3).eq.ix(3)).AND.(iaa(ii)%cm(j,4).eq.ix(4)).AND.(iaa(ii)%cm(j,5).eq.ix(5))) then
            foundit = j
            rs = atmres(ix(i))
            exit
          end if
        end do
        if (foundit.gt.0) exit
      end do    
      if ((foundit.le.0).OR.(ix(6).le.0).OR.(ix(6).gt.n_cmstyp)) then
        write(ilog,35) ix(1),ix(2),ix(3),ix(4),ix(5),ix(6)
        cycle
      else
        nlst = nlst + 1
        if (iaa(rs)%typ_cm(j).le.0) then
          allm(5) = .false.
          write(ilog,37) ix(6),ix(1),ix(2),ix(3),ix(4),ix(5),rs
        else
          write(ilog,36) ix(1),ix(2),ix(3),ix(4),ix(5),rs,ix(6)
        end if
        write(ilog,38) b_type(ix(1)),bio_code(b_type(ix(1))),b_type(ix(2)),bio_code(b_type(ix(2))),b_type(ix(3)),&
 &bio_code(b_type(ix(3))),b_type(ix(4)),bio_code(b_type(ix(4))),b_type(ix(5)),bio_code(b_type(ix(5))),&
 &amino(seqtyp(atmres(ix(3))))
        iaa(rs)%typ_cm(j) = ix(6)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported keyword in bonded parameters patch file&
 & ',bpatchfile(t1:t2),'. Use only PATCH_BOND, PATCH_ANGLE, PATCH_IMPROPER, PATCH_TORSION, or PATCH_CMAP.'
      call fexit()
    end if
!
  end do
!
  close(unit=iunit)
!
  if (nlst.gt.0) then
    bonded_patched = .true.
    write(ilog,*)
  else
    write(ilog,'(1x,a)') 'Failed to read in any eligible patches (check input).'
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_fospatchfile(tmpmat,tmpmat2,nlst)
!
  use energies
  use iounit
  use atoms
  use params
  use fos
  use interfaces
  use sequen
  use molecule
  use polypep
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,t3,t4,i,j,k,ii,kk,nlst,mn,mng,maxid,iref,lk,readrfos
  character(MAXSTRLEN) str2,str3,str4
  character(10) abc
  logical exists,wrnd,afalse
  RTYPE mnwt,mnrfos(3),tmpmat2(n,4)
  integer tmpmat(n,2)
  integer, ALLOCATABLE:: hlp(:,:)
  logical, ALLOCATABLE:: fdvec(:)
!
  afalse = .false.
  nlst = 0
  if (use_IMPSOLV.EQV..false.) return
!
! the only thing we already checked for is the existence of the actual input file,
! UNLESS mr.user elected to use the default -> check again
  call strlims(fospatchfile,t1,t2)
  inquire(file=fospatchfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'No FOS patches applied (no file specified or file corrupt/unreadable).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()

  open (unit=iunit,file=fospatchfile(t1:t2),status='old')
  call strlims(paramfile,t3,t4)
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for FOS patch &
 &input (',fospatchfile(t1:t2),'). Using default group-based FOS from parameter file ',paramfile(t3:t4),'.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing FOS patch input (',fospatchfile(t1:t2),').'
    call fexit()
  end if
  ii = 1
  call extract_str(str2,abc,ii)
!
! now just read present entries
!
  allocate(fdvec(n))
  fdvec(:) = .false.
!
  do while(.true.)
    read(iunit,79,iostat=iomess) str2
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_fospatchfile(...).'
      call fexit()
    end if
!
    ii = 1
    mn = -1
    call extract_int(str2,mn,ii)
    if (ii.eq.1) cycle
    if ((mn.gt.n).OR.(mn.le.0)) then
      write(ilog,*) 'Requested atom number ',&
 &mn,' for applying a FOS patch. It does not exist, however, and the request will be ignored ...'
      cycle
    end if
    str3(:) = ' '
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    read(str3,*,iostat=iomess) mng
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Fatal. Incomplete line while processing file &
 &input in read_fospatchfile(...).'
      call fexit()
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_fospatchfile(...).'
      call fexit()
    end if
    ii = 1
    k = -56879214
    call extract_int(str3,k,ii)
    str2(:) = ' '
    str2(1:MAXSTRLEN-ii+1) = str3(ii:MAXSTRLEN)
    read(str2,*,iostat=iomess) mnwt
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Fatal. Incomplete line while processing file &
 &input in read_fospatchfile(...).'
      call fexit()
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing file &
 &input in read_fospatchfile(...).'
      call fexit()
    end if
    if (mnwt.lt.0.0) then
      write(ilog,*) 'Requested negative atomic weight (',&
 &mnwt,') in applying a FOS patch. Request will be ignored (most likely fatal later on) ...'
      cycle
    end if
    ii = 1
    call extract_str(str2,str4,ii)
    str3(:) = ' '
    str3(1:MAXSTRLEN-ii+1) = str2(ii:MAXSTRLEN)
    i = 1
    mnrfos(:) = 0.0
    readrfos = 0
    lk = 3
    call get_reals_from_str(str3,lk,mnrfos(1:3),i,readrfos,afalse)
!
    nlst = nlst + 1
    if (fdvec(mn).EQV..false.) then
      tmpmat(nlst,1) = mn
      tmpmat(nlst,2) = mng
      tmpmat2(nlst,1) = mnwt
      if (readrfos.gt.0) then
        tmpmat2(nlst,2:4) = mnrfos(1:3)
        if (readrfos.eq.1) then
          tmpmat2(nlst,3) = mnrfos(1)
          tmpmat2(nlst,4) = 0.0
        else if (readrfos.eq.2) then
          tmpmat2(nlst,4) = 0.0
        end if
      else
        readrfos = 0
        do j=1,nlst-1
          if (mng.eq.tmpmat(j,2)) then
            tmpmat2(nlst,2:4) = tmpmat2(j,2:4)
            readrfos = 1
            exit
          end if
        end do
        if (readrfos.eq.0) then
          write(ilog,*) 'Fatal. First atom to a newly defined FOS group must specify reference &
 &free energy of solvation in applying a FOS patch.'
          call fexit()
        end if
      end if
      fdvec(mn) = .true.
    else
      nlst = nlst - 1
      write(ilog,*) 'Warning. List provided to read_fospatchfile(...) contains &
 &double entry ( atom # ',mn,'). Ignored (first instance used) ...'
    end if
!
  end do
!
  close(unit=iunit)
!
! sort by index of group
  allocate(hlp(nlst,5))
  do i=1,nlst
    hlp(i,3) = i
  end do
  wrnd = .TRUE.
  i = 1
  j = nlst
  k = nlst
  hlp(:,4) = tmpmat(1:nlst,2)
  call merge_sort(ldim=k,up=wrnd,list=hlp(:,4),olist=hlp(:,1),ilo=i,ihi=j,idxmap=hlp(:,3),olist2=hlp(:,2))
!  
  do i=1,nlst
    hlp(hlp(i,2),5) = i
  end do
  do i=1,nlst
    if (i.eq.hlp(i,2)) cycle
    do k=1,4
      mnwt = tmpmat2(i,k)
      tmpmat2(i,k) = tmpmat2(hlp(i,2),k)
      tmpmat2(hlp(i,2),k) = mnwt
    end do
    do k=1,2
      mn = tmpmat(i,k)
      tmpmat(i,k) = tmpmat(hlp(i,2),k)
      tmpmat(hlp(i,2),k) = mn
    end do
    k = hlp(i,5)
    mn = hlp(i,2)
    hlp(i,5) = hlp(mn,5)
    hlp(mn,5) = k
    hlp(k,2) = hlp(i,2)
  end do
  deallocate(hlp)
!
  23 format(i6,1x,i6,1x,g12.5,1x,g12.5)
  22 format(a,i6,a,g10.5,a)
  21 format(a,i6,a,i8,a)
  i = 1
  maxid = 0
  wrnd = .false.
  do while (i.le.nlst)
    mng = tmpmat(i,2)
    k = atmres(tmpmat(i,1))
    if ((abc(1:1).eq.'t').OR.(abc(1:1).eq.'T')) then
      if (molofrs(k).ne.moltyp(moltypid(molofrs(k)),1)) then
        write(ilog,22) 'Fatal. Atom starting patched FOS group labeled ',mng,' is not part of the first &
 &molecule of that type. This is fatal in conjunction with input mode "T".'
        call fexit()
      end if
    end if
    j = i
    mnwt = 0.0
    if (mng.gt.maxid) maxid = mng
    do while (tmpmat(i,2).eq.mng)
      mnwt = mnwt + tmpmat2(i,1)
      if (i.gt.j) tmpmat2(i,2:4) = tmpmat2(j,2:4)
      if (molofrs(k).ne.molofrs(atmres(tmpmat(i,1)))) then
        write(ilog,22) 'Fatal. Atoms constituting patched FOS group labeled ',mng,' are not all part of &
 &the same molecule (first is ',molofrs(k),').'
        call fexit()
      end if
      if (k.ne.atmres(tmpmat(i,1))) then
        if ((abc(1:1).eq.'R').OR.(abc(1:1).eq.'r')) then
          write(ilog,22) 'Fatal. Atoms constituting patched FOS group labeled ',mng,' are not all part of &
 &the same residue (first is ',k,').'
          write(ilog,*) 'In mode R, this is currently not supported. Please adjust input file.'
          call fexit()
        else
          if (wrnd.EQV..false.) then
            wrnd = .true.
            write(ilog,22) 'Warning. Atoms constituting patched FOS group labeled ',mng,' are not all part of &
 &the same residue (first is ',k,').'
            write(ilog,*) 'There is only very limited support for such requests at the moment (further warnings omitted).'
          end if
        end if
      end if
      i = i + 1
      if (i.gt.nlst) exit
    end do
    if (abs(mnwt-1.0).gt.1.0e-5) then
      write(ilog,22) 'Fatal. Sum of atomic weights for patched FOS group labeled ',mng,' is not close enough to &
 &unity (',mnwt,'). Please correct input file.'
      call fexit()
    end if
  end do
!
  if ((abc(1:1).eq.'t').OR.(abc(1:1).eq.'T')) then
    i = 1
    mn = nlst
    do while (i.le.mn)
      k = molofrs(atmres(tmpmat(i,1)))
      mng = tmpmat(i,2)
      iref = i
      do while (tmpmat(i,2).eq.mng)
        i = i + 1
        if (i.gt.mn) exit
      end do
      do j=1,nmol
        if (moltypid(j).ne.moltypid(k)) cycle
        if (moltyp(moltypid(k),1).eq.j) cycle
        maxid = maxid + 1
        do ii=1,(i-iref)
          nlst = nlst + 1
          tmpmat(nlst,1) = tmpmat(iref+ii-1,1) + atmol(j,1) - atmol(k,1)
          tmpmat(nlst,2) = maxid
          tmpmat2(nlst,:) = tmpmat2(iref+ii-1,:)
        end do
      end do
    end do
  else if ((abc(1:1).eq.'r').OR.(abc(1:1).eq.'R')) then
    i = 1
    mn = nlst
    do while (i.le.mn)
      k = seqtyp(atmres(tmpmat(i,1)))
      mng = tmpmat(i,2)
      iref = i
      do while (tmpmat(i,2).eq.mng)
        i = i + 1
        if (i.gt.mn) exit
      end do
      mng = atmres(tmpmat(iref,1))
      kk = 0
      do j=1,nseq
        if (seqtyp(j).ne.k) cycle
        if ((j.eq.rsmol(molofrs(j),1)).AND.(j.eq.rsmol(molofrs(j),2)).AND.&
 &          ((mng.ne.rsmol(molofrs(mng),1)).OR.(mng.ne.rsmol(molofrs(mng),2)))) cycle
        if ((mng.eq.rsmol(molofrs(mng),1)).AND.(mng.eq.rsmol(molofrs(mng),2)).AND.&
 &          ((j.ne.rsmol(molofrs(j),1)).OR.(j.ne.rsmol(molofrs(j),2)))) cycle
        if ((j.eq.rsmol(molofrs(j),1)).AND.(mng.ne.rsmol(molofrs(mng),1))) cycle
        if ((j.ne.rsmol(molofrs(j),1)).AND.(mng.eq.rsmol(molofrs(mng),1))) cycle
        if ((j.eq.rsmol(molofrs(j),2)).AND.(mng.ne.rsmol(molofrs(mng),2))) cycle
        if ((j.ne.rsmol(molofrs(j),2)).AND.(mng.eq.rsmol(molofrs(mng),2))) cycle
        exists = .false.
        do lk=1,n_crosslinks
          if ((crosslink(lk)%rsnrs(1).eq.j).OR.(crosslink(lk)%rsnrs(2).eq.j).OR.&
 &            (crosslink(lk)%rsnrs(1).eq.mng).OR.(crosslink(lk)%rsnrs(2).eq.mng)) then
            exists = .true.
            exit
          end if
        end do
        if (exists.EQV..true.) cycle        
        kk = kk + 1
        if (kk.eq.1) then
          if (j.ne.mng) then
            write(ilog,22) 'Fatal. Atom starting patched FOS group labeled ',tmpmat(iref,2),' is not part of the first &
 &residue that type. This is fatal in conjunction with input mode "R".'
            call fexit()
          end if
        end if
        if (j.eq.mng) cycle
        maxid = maxid + 1
        do ii=1,(i-iref)
          nlst = nlst + 1
          tmpmat(nlst,1) = tmpmat(iref+ii-1,1) + at(j)%bb(1) - at(mng)%bb(1)
          tmpmat(nlst,2) = maxid
          tmpmat2(nlst,:) = tmpmat2(iref+ii-1,:)
        end do
      end do
    end do
  else if ((abc(1:1).eq.'a').OR.(abc(1:1).eq.'A')) then
!    do nothing
  else
    write(ilog,*) 'Fatal. Unrecognized or unsupported input mode in read_fospatchfile(...).'
    call fexit()
  end if
!  do i=1,nlst
!    write(ilog,23) tmpmat(i,1),tmpmat(i,2),tmpmat2(i,1),tmpmat2(i,2)
!  end do
!  write(ilog,*)
  deallocate(fdvec)
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_drestfile()
!
  use energies
  use iounit
  use sequen
  use distrest
  use aminos
  use atoms
  use params
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,i,ii,datt(2),aone,atwo,daty
  character(MAXSTRLEN) str2,str3
  logical exists
  RTYPE daprm(2)
!
  ndrest = 0
!
  call strlims(drestfile,t1,t2)
  inquire(file=drestfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (&
 &',drestfile(t1:t2),') for FMCSC_DRESTFILE. Turning off DR terms.'
    use_DREST = .false.
    scale_DREST = 0.0
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  aone = 1
  atwo = 2
!
  iunit = freeunit()

  open (unit=iunit,file=drestfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for distance r&
 &estraint input (',drestfile(t1:t2),'). Turning off this term.'
    use_DREST = .false.
    scale_DREST = 0.0
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing distance r&
 &estraint input (got: ',drestfile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  allndr = -1
  call extract_int(str2,allndr,ii)
!
  if (allndr.le.0) then
    write(ilog,*) 'Warning. The first line of the distance restraint input file (',drestfile(t1:t2),') must &
 &give the number of terms to be read. Turning term off due to uninterpretable file.'
    allndr = 0
    use_DREST = .false.
    scale_DREST = 0.0
    close(unit=iunit)
    return
  end if
!
  call allocate_distrest(aone)
!
! now just read the requested number of entries
!
  do i=1,allndr
    read(iunit,79,iostat=iomess) str3
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Warning. Got empty/incomplete file for distance&
 & restraint input (',drestfile(t1:t2),'). Turning off this term.'
      use_DREST = .false.
      scale_DREST = 0.0
      ndrest = 0
      allndr = 0
      close(unit=iunit)
      call allocate_distrest(atwo)
      return
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing distance&
 & restraint input (got: ',drestfile(t1:t2),').'
      call fexit()
    end if
!
    daty = 0
    datt(1:2) = 0
    daprm(:) = 0.
    read(str3,*,iostat=iomess) datt(1),datt(2),daty,daprm(1),daprm(2)
    if ((datt(1).le.0).OR.(datt(2).le.0).OR.&
 &      (datt(1).gt.n).OR.(datt(2).gt.n)) then
      write(ilog,*) 'Illegal entry in distance restraint input (atom&
 & numbers out of range in ',drestfile(t1:t2),').'
      cycle
    else if ((daty.lt.1).OR.(daty.gt.3)) then
      write(ilog,*) 'Illegal entry in distance restraint input (unsupported&
 & restraint type in ',drestfile(t1:t2),').'
      cycle
    else if (daprm(1).lt.0.0) then
      write(ilog,*) 'Illegal entry in distance restraint input (nega&
 &tive minimum distance in ',drestfile(t1:t2),').'
      cycle
    else if ((daty.le.3).AND.(daprm(2).lt.0.0)) then
      write(ilog,*) 'Got negative force constant for (half-)harmonic restraint in distance restrai&
 &nt input (',drestfile(t1:t2),'). Assuming ',-daprm(2),' instead.'
      daprm(2) = -daprm(2)
    end if
!
    ndrest = ndrest + 1
    drestyp(ndrest) = daty
    dresat(ndrest,:) = datt(:)
    dresprm(ndrest,:) = daprm(:)
!
  end do
!
  close(unit=iunit)
!
  if (ndrest.eq.0) then
    write(ilog,*) 'No valid entries in distance restraint input (',&
 &drestfile(t1:t2),'). Turning term off.'
    use_DREST = .false.
    scale_DREST = 0.0
    allndr = 0
    return
  end if    
!
 24   format('Atoms #',i6,' (',a4,' in ',a3,') and '&
 &                 ,i6,' (',a4,' in ',a3,'))')
 25   format(8x,'Type: ',i2,' D_0: ',f7.2,' K_d: ',f7.2)
!
  if (drest_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '-- Summary of Distance Restraint Terms in the Sys&
 &tem --'
    write(ilog,*)
    do i=1,ndrest
      write(ilog,24) dresat(i,1),bio_code(b_type(dresat(i,1))),&
 &                   amino(seqtyp(atmres(dresat(i,1)))),&
 &                   dresat(i,2),bio_code(b_type(dresat(i,2))),&
 &                   amino(seqtyp(atmres(dresat(i,2))))
      write(ilog,25) drestyp(i),dresprm(i,1),dresprm(i,2)
    end do
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_gpcfile()
!
  use iounit
  use sequen
  use molecule
  use aminos
  use atoms
  use paircorr
  use params
  use mpistuff
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,i,j,k,szz,aone,atwo,dum(3)
  integer, ALLOCATABLE:: dummat(:,:),pclst(:)
  integer ii,jj,nm1,nm2,typ1,typ2,imol,jmol,kmol,ati,atj,pcnlst
  character(10) abc
  character(MAXSTRLEN) str2,str3
  logical exists,novel,havetogo
!
! the only thing we already checked for is the existence of the actual input file,
! UNLESS mr.user elected to use the default -> check again
  call strlims(pccodefile,t1,t2)
  inquire(file=pccodefile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (&
 &',pccodefile(t1:t2),') for FMCSC_PCCODEFILE. Turning off general P&
 &C analysis.'
    gpc%nos = 0
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
 78   format('Skipping #1:',i6,' #2:',i6,' w/ PC#:',i6)
 77   format('Atoms ',i6,' and ',i6,' are in mol. ',i6,', not ',i6,'.')
 76   format('Atom ',i6,' is part of mol. ',i6,', not ',i6,'.')
 75   format('Transferring input # ',i6,' to output number #',i6,'.')
 74   format('Skipping mode: ',a1,' w/ #1:',i6,' and #2:',i6)
 73   format('Skipping mode: ',a1,' w/ type:',i6)
!
  aone = 1
  atwo = 2
!
  iunit = freeunit()

  open (unit=iunit,file=pccodefile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for general PC&
 &input (',pccodefile(t1:t2),'). Turning off this analysis.'
    gpc%nos = 0
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing general PC&
 & input: ',pccodefile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_int(str2,gpc%mode,ii) 
  str3 = str2(ii:MAXSTRLEN)
  szz = -1
  call extract_int(str3,szz,ii)
!
! first list-based formats
  if ((gpc%mode.eq.1).OR.(gpc%mode.eq.2)) then
    if (szz.le.0) then
      write(ilog,*) 'Warning. Second number on header line has to provide number of entries for list-based general&
 & PC input (',pccodefile(t1:t2),'). Turning off analysis.'
      gpc%nos = 0
      close(unit=iunit)
      return
    end if
    gpc%all_lst = szz
    gpc%all_atm = 0
    call allocate_gpcpre(aone)
!
    gpc%nlst = 0
    do i=1,gpc%all_lst
      read(iunit,79,iostat=iomess) str3
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for genera&
 &l PC input (',pccodefile(t1:t2),'). Turning off analysis.'
        gpc%nos = 0
        close(unit=iunit)
        call allocate_gpcpre(atwo)
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing genera&
 &l PC input (got: ',pccodefile(t1:t2),').'
        call fexit()
      end if
!
      dum(1) = -1
      dum(2) = -1
      dum(3) = -1
      read(str3,*,iostat=iomess) dum(1),dum(2),dum(3)
!
      if (gpc%mode.eq.1) then
!       the atom-specific format is simple enough
        if ((dum(1).le.0).OR.(dum(1).gt.n).OR.&
 &          (dum(2).le.0).OR.(dum(2).gt.n).OR.(dum(3).le.0)) then
          write(ilog,*) 'Range exception occurred while reading GPC &
 &input:'
          write(ilog,78) dum(1),dum(2),dum(3)
          cycle
        end if
        gpc%nlst = gpc%nlst + 1
        gpc%lst(gpc%nlst,1) = dum(1)
        gpc%lst(gpc%nlst,2) = dum(2)
        gpc%lst(gpc%nlst,3) = dum(3)
        gpc%lty(gpc%nlst) = .false.
!       more sanity checks are needed for the molecule-type-generalized format
      else if (gpc%mode.eq.2) then
        if ((dum(1).le.0).OR.(dum(1).gt.n).OR.&
 &          (dum(2).le.0).OR.(dum(2).gt.n).OR.(dum(3).le.0)) then
          write(ilog,*) 'Range exception occurred while reading GPC &
 &input:'
          write(ilog,78) dum(1),dum(2),dum(3)
          cycle
        end if
        imol = molofrs(atmres(dum(1)))
        jmol = molofrs(atmres(dum(2)))
!       the first case is an intramolecular term:
!       this has to be the first mol. of its type
        if (imol.eq.jmol) then
          if (imol.ne.molangr(an_grp_mol(imol),1)) then
            write(ilog,*) 'Atoms for intramolecular terms have to co&
 &rrespond to the first molecule in its analysis group in sequence:'
            write(ilog,77) dum(1),dum(2),imol,&
 &                           molangr(an_grp_mol(imol),1)
            cycle
          end if
          gpc%nlst = gpc%nlst + 1
          gpc%lst(gpc%nlst,1) = dum(1)
          gpc%lst(gpc%nlst,2) = dum(2)
          gpc%lst(gpc%nlst,3) = dum(3)
          gpc%lty(gpc%nlst) = .true.
!       the second case is an intermolecular term between different analysis groups:
!       both mol.s have to be the first occurrence of their respective group
        else if (an_grp_mol(imol).ne.an_grp_mol(jmol)) then
          if ((imol.ne.molangr(an_grp_mol(imol),1)).OR.&
 &            (jmol.ne.molangr(an_grp_mol(jmol),1))) then
            write(ilog,*) 'Atoms for intermolecular terms have&
 & to correspond to the first molecules in the two (different) analysis groups in sequence:'
            if (imol.ne.molangr(an_grp_mol(imol),1)) then
              write(ilog,76) dum(1),imol,molangr(an_grp_mol(imol),1)
            else
              write(ilog,76) dum(2),jmol,molangr(an_grp_mol(jmol),1)
            end if
            cycle
          end if
          gpc%nlst = gpc%nlst + 1
          gpc%lst(gpc%nlst,1) = dum(1)
          gpc%lst(gpc%nlst,2) = dum(2)
          gpc%lst(gpc%nlst,3) = dum(3)
          gpc%lty(gpc%nlst) = .true.
!       the third and final case is an intermolecular term between the same analysis group:
!       the first mol has to be the first in its group, the second the second
        else
          if (imol.ne.molangr(an_grp_mol(imol),1)) then
            write(ilog,*) 'Atoms for group-internal intermolecular terms have &
 &to correspond to the first two molecules in their analysis group in sequence:'
            write(ilog,76) dum(1),imol,molangr(an_grp_mol(imol),1)
            cycle
          end if
          do kmol=1,nmol
            if ((an_grp_mol(kmol).eq.an_grp_mol(imol)).AND.(kmol.ne.imol)) exit
          end do
          if (jmol.ne.kmol) then
            write(ilog,*) 'Atoms for group-internal intermolecular terms have &
 &to correspond to the first two molecules in their analysis group in sequence:'
            write(ilog,76) dum(2),jmol,kmol
            cycle
          end if
          gpc%nlst = gpc%nlst + 1
          gpc%lst(gpc%nlst,1) = dum(1)
          gpc%lst(gpc%nlst,2) = dum(2)
          gpc%lst(gpc%nlst,3) = dum(3)
          gpc%lty(gpc%nlst) = .true.
        end if
      end if
    end do
!
! second matrix-based format
  else if (gpc%mode.eq.3) then
!
!   first do a dry run to print out all the warnings and to determine the needed memory
    do while (.true.)
      havetogo = .false.
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while processing GPC input.'
        call fexit()
      end if
      abc(1:1) = ' '
      typ1 = -1
      typ2 = -1
      read(str2,*,iostat=iomess) abc,typ1,typ2
!     intermolecular type-wise requests
      if (abc(1:1).eq.'T') then
        if ((typ1.le.0).OR.(typ1.gt.nangrps).OR.&
 &          (typ2.le.0).OR.(typ2.gt.nangrps)) then
          write(ilog,*) 'Range exception occurred while processing GPC input:'
          write(ilog,74) abc(1:1),typ1,typ2
          cycle
        end if
        if (typ1.eq.typ2) then
          if (molangr(typ1,2).eq.1) then
            write(ilog,*) 'Intermolecular term for analysis group ',&
 &typ1,' is ignored, as there is only one molecule in this group.'
            cycle
          end if
          do kmol=1,nmol
            if ((an_grp_mol(kmol).eq.typ1).AND.(kmol.ne.molangr(typ1,1))) exit
          end do
          nm1 = atmol(molangr(typ1,1),2)-atmol(molangr(typ1,1),1)+1
          allocate(dummat(nm1,nm1))
!         this is interpreted to be an intermolecular term, but the
!         matrix is symmetric and only half+diag of it is used (has to be
!         fully there, though!)
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
            if (iomess.eq.IOSTAT_END) then
              write(ilog,*) 'Warning. Incomplete matrix found while processing GPC input:'
              write(ilog,74) abc(1:1),typ1,typ2
              havetogo = .true.
              exit
            else if (iomess.eq.2) then
              write(ilog,*) 'Fatal. I/O error while processing GPC input.'
              call fexit()
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=ii,nm1
              if (dummat(ii,jj).gt.0) then
                gpc%all_lst = gpc%all_lst + 1
              end if
            end do
          end do
          deallocate(dummat)
        else
          nm1 = atmol(molangr(typ1,1),2)-atmol(molangr(typ1,1),1)+1
          nm2 = atmol(molangr(typ2,1),2)-atmol(molangr(typ2,1),1)+1
          allocate(dummat(nm1,nm2))
!         this is by necessity an intermolecular term, with a
!         non-symmetric matrix
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm2)
            if (iomess.eq.IOSTAT_END) then
              write(ilog,*) 'Warning. Incomplete matrix found while processing GPC input:'
              write(ilog,74) abc(1:1),typ1,typ2
              havetogo = .true.
              exit
            else if (iomess.eq.2) then
              write(ilog,*) 'Fatal. I/O error while processing GPC input.'
              call fexit()
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=1,nm2
              if (dummat(ii,jj).gt.0) then
                gpc%all_lst = gpc%all_lst + 1
              end if
            end do
          end do
          deallocate(dummat)
        end if
!     intramolecular type-wise requests
      else if (abc(1:1).eq.'t') then
        if ((typ1.le.0).OR.(typ1.gt.nangrps)) then
          write(ilog,*) 'Range exception occurred while processing GPC input:'
          write(ilog,73) abc(1:1),typ1
          cycle
        end if
        nm1 = atmol(molangr(typ1,1),2)-atmol(molangr(typ1,1),1)+1
        allocate(dummat(nm1,nm1))
!       as an intramolecular term, the matrix is symmetric and only half sans diag
!       of it is used (has to be fully there, though!)
        do ii=1,nm1
          read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
          if (iomess.eq.IOSTAT_END) then
            write(ilog,*) 'Warning. Incomplete matrix found while processing GPC input:'
            write(ilog,73) abc(1:1),typ1
            havetogo = .true.
            exit
          else if (iomess.eq.2) then
            write(ilog,*) 'Fatal. I/O error while processing GPC input.'
            call fexit()
          end if
        end do
        if (havetogo.EQV..true.) then
          deallocate(dummat)
          exit
        end if
        do ii=1,nm1
          do jj=ii+1,nm1
            if (dummat(ii,jj).gt.0) then
              gpc%all_lst = gpc%all_lst + 1
            end if
          end do
        end do
        deallocate(dummat)
      end if
    end do
!
!   now allocate memory, rewind file, and read-off header line
    call allocate_gpcpre(aone)
    rewind(unit=iunit)
    read(iunit,79,iostat=iomess) str3
    do while (.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      end if
      abc = ' '
      typ1 = -1
      typ2 = -1
      read(str2,*,iostat=iomess) abc,typ1,typ2
!     intermolecular type-wise requests
      if (abc.eq.'T') then
        if ((typ1.le.0).OR.(typ1.gt.nangrps).OR.(typ2.le.0).OR.(typ2.gt.nangrps)) then
          cycle
        end if
        if (typ1.eq.typ2) then
          if (molangr(typ1,2).eq.1) then
            cycle
          end if
          do kmol=1,nmol
            if ((an_grp_mol(kmol).eq.typ1).AND.(kmol.ne.molangr(typ1,1))) exit
          end do
          nm1 = atmol(molangr(typ1,1),2)-atmol(molangr(typ1,1),1)+1
          allocate(dummat(nm1,nm1))
!         this is interpreted to be an intermolecular term, but the
!         matrix is symmetric and only half+diag of it is used (has to be
!         fully there, though!)
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
            if (iomess.eq.IOSTAT_END) then
              havetogo = .true.
              exit
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=ii,nm1
              if (dummat(ii,jj).gt.0) then
                gpc%nlst = gpc%nlst + 1
                gpc%lst(gpc%nlst,1) = atmol(molangr(typ1,1),1)+ii-1
                gpc%lst(gpc%nlst,2) = atmol(kmol,1)+jj-1
                gpc%lst(gpc%nlst,3) = dummat(ii,jj)
                gpc%lty(gpc%nlst) = .true.
              end if
            end do
          end do
          deallocate(dummat)
        else
          nm1 = atmol(molangr(typ1,1),2)-atmol(molangr(typ1,1),1)+1
          nm2 = atmol(molangr(typ2,1),2)-atmol(molangr(typ2,1),1)+1
          allocate(dummat(nm1,nm2))
!         this is by necessity an intermolecular term, with a
!         non-symmetric matrix
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm2)
            if (iomess.eq.IOSTAT_END) then
              havetogo = .true.
              exit
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=1,nm2
              if (dummat(ii,jj).gt.0) then
                gpc%nlst = gpc%nlst + 1
                gpc%lst(gpc%nlst,1) = atmol(molangr(typ1,1),1)+ii-1
                gpc%lst(gpc%nlst,2) = atmol(molangr(typ2,1),1)+jj-1
                gpc%lst(gpc%nlst,3) = dummat(ii,jj)
                gpc%lty(gpc%nlst) = .true.
              end if
            end do
          end do
          deallocate(dummat)
        end if
!     intramolecular type-wise requests
      else if (abc(1:1).eq.'t') then
        if ((typ1.le.0).OR.(typ1.gt.nangrps)) then
          cycle
        end if
        nm1 = atmol(molangr(typ1,1),2)-atmol(molangr(typ1,1),1)+1
        allocate(dummat(nm1,nm1))
        do ii=1,nm1
          read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
          if (iomess.eq.IOSTAT_END) then
            havetogo = .true.
            exit
          end if
        end do
        if (havetogo.EQV..true.) then
          deallocate(dummat)
          exit
        end if
        do ii=1,nm1
          do jj=ii+1,nm1
            if (dummat(ii,jj).gt.0) then
              gpc%nlst = gpc%nlst + 1
              gpc%lst(gpc%nlst,1) = atmol(molangr(typ1,1),1)+ii-1
              gpc%lst(gpc%nlst,2) = atmol(molangr(typ1,1),1)+jj-1
              gpc%lst(gpc%nlst,3) = dummat(ii,jj)
              gpc%lty(gpc%nlst) = .true.
            end if
          end do
        end do
        deallocate(dummat)
      end if
    end do
!
  else
    write(ilog,*) 'Unknown mode for general PC analysis (got: ', &
 &gpc%mode,'). Turning feature off.'
    gpc%nos = 0
    close(unit=iunit)
    return
  end if
!
  close(unit=iunit)
!
  gpc%nos = 0
  allocate(pclst(gpc%all_lst))
  do i=1,gpc%nlst
    novel = .true.
    do j=1,gpc%nos
      if (gpc%lst(i,3).eq.pclst(j)) then
        novel = .false.
        exit
      end if
    end do
    if (novel.EQV..true.) then
      gpc%nos = gpc%nos + 1
      pclst(gpc%nos) = gpc%lst(i,3)
    end if
  end do
  if (gpc%nos.eq.0) then
    write(ilog,*) 'No general PCs requested in ',pccodefile(t1:t2),'&
 &? Turning feature off.'
    deallocate(pclst)
    call allocate_gpcpre(atwo)
    return
  end if
  call sort(gpc%nos,pclst)
!
  do j=1,gpc%nos
    if (pclst(j).ne.j) then
      write(ilog,*) 'Discovered discontinous numbering of PCs while &
 &processing GPC input.'
      do jj=j,gpc%nos
        write(ilog,75) pclst(jj),jj
        do i=1,gpc%nlst
          if (gpc%lst(i,3).eq.pclst(jj)) gpc%lst(i,3) = jj
        end do
      end do
      exit
    end if
  end do
!
  deallocate(pclst)
!
! lastly a sanity check against the illegal mixing of intra- and intermolecular terms
  do i=1,gpc%nlst
    imol = molofrs(atmres(gpc%lst(i,1)))
    jmol = molofrs(atmres(gpc%lst(i,2)))
    if (imol.eq.jmol) then
      do j=i+1,gpc%nlst
        if (gpc%lst(j,3).ne.gpc%lst(i,3)) cycle
        if (molofrs(atmres(gpc%lst(j,1))).ne.molofrs(atmres(gpc%lst(j,2)))) then
          write(ilog,*) 'Fatal. Discovered attempt to unify intra- and intermolecular&
 & terms under the same index while processing GPC input. This is not allowed due to&
 & the different types of averaging performed for both.'
          call fexit()
        end if
      end do
    else
      do j=i+1,gpc%nlst
        if (gpc%lst(j,3).ne.gpc%lst(i,3)) cycle
        if (molofrs(atmres(gpc%lst(j,1))).eq.molofrs(atmres(gpc%lst(j,2)))) then
          write(ilog,*) 'Fatal. Discovered attempt to unify intra- and intermolecular&
 & terms under the same index while processing GPC input. This is not allowed due to&
 & the different types of averaging performed for both.'
          call fexit()
        end if
      end do
    end if
  end do
!
  if (gpc_report.EQV..true.) then
#ifdef ENABLE_MPI
    if ((use_REMC.EQV..true.).AND.(myrank.eq.0)) then
      str3 = 'N_000_GENERAL_PC.idx'
    else if ((use_MPIAVG.EQV..true.).AND.(myrank.eq.0)) then
      str3 = 'GENERAL_PC.idx'
    else
      return
    end if
#else
    str3 = 'GENERAL_PC.idx'
#endif
    call strlims(str3,t1,t2)
!
   99   format('--- Generalized Distance Distribution Component #',i5)
   98   format('Specific: Atom ',i6,' (',a4,') in Res. ',i6,' (',a3,')')
   97   format('          Atom ',i6,' (',a4,') in Res. ',i6,' (',a3,')')
   96   format('Group-wise and|',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &                                 ') of An. Group ',i4)
   95   format('Intramolecular|',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &')')
   94   format('Group-wise and|',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &                                 ') of An. Group ',i4)
   93   format('Intermolecular|',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &                                 ') of An. Group ',i4)
!
    inquire(file=str3(t1:t2),exist=exists)
    if(exists) then
      iunit = freeunit()
      open(unit=iunit,file=str3(t1:t2),status='old')
      close(unit=iunit,status='delete')
    end if
    iunit=freeunit()
    open(unit=iunit,file=str3(t1:t2),status='new')
!
    write(iunit,*) 'List of Terms for Columns in GENERAL_PC.dat'
    write(iunit,*) '-------------------------------------------'
    write(iunit,*)
    allocate(pclst(gpc%nlst))
    pcnlst = 0
    do i=1,gpc%nos
      do j=1,gpc%nlst
        if (gpc%lst(j,3).eq.i) then
          pcnlst = pcnlst + 1
          pclst(pcnlst) = j
        end if
      end do
!
      write(iunit,99) i
      do j=1,pcnlst
        k = pclst(j)
        ati = gpc%lst(k,1)
        atj = gpc%lst(k,2)
        if (gpc%lty(k).EQV..false.) then
          write(iunit,98) ati,bio_code(b_type(ati)),atmres(ati),&
 &                     amino(seqtyp(atmres(ati)))
          write(iunit,97) atj,bio_code(b_type(atj)),atmres(atj),&
 &                     amino(seqtyp(atmres(atj)))
        else
          imol = molofrs(atmres(ati))
          jmol = molofrs(atmres(atj))
          ii = atmres(ati)
          jj = atmres(atj)
          if (imol.eq.jmol) then
            write(iunit,96)bio_code(b_type(ati)),ii-rsmol(imol,1)+1,&
 &amino(seqtyp(ii)),an_grp_mol(imol)
            write(iunit,95)bio_code(b_type(atj)),jj-rsmol(jmol,1)+1,&
 &amino(seqtyp(jj))
          else
            write(iunit,94)bio_code(b_type(ati)),ii-rsmol(imol,1)+1,&
 &amino(seqtyp(ii)),an_grp_mol(imol)
            write(iunit,93)bio_code(b_type(atj)),jj-rsmol(jmol,1)+1,&
 &amino(seqtyp(jj)),an_grp_mol(jmol)
          end if
        end if
      end do
      write(iunit,*)
      pcnlst = 0
    end do
    deallocate(pclst)
    close(unit=iunit)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_tabfiles()
!
  use iounit
  use sequen
  use molecule
  use aminos
  use atoms
  use tabpot
  use params
  use mpistuff
  use energies
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,i,j,k,szz,aone,atwo,dum(3)
  integer, ALLOCATABLE:: dummat(:,:),pclst(:)
  RTYPE, ALLOCATABLE:: potrs(:)
  RTYPE disres,disold,eps,fac(2)
  integer ii,jj,nm1,nm2,typ1,typ2,imol,jmol,kmol,ati,atj,pcnlst
  integer rs1,rs2,rs1o,rs2o,off1,off2,oo1,oo2
  character(10) abc
  character(MAXSTRLEN) str2,str3
  logical exists,novel,havetogo
!
! the only thing we already checked for is the existence of the actual input file,
! UNLESS mr.user elected to use the default -> check again
  call strlims(tabcodefile,t1,t2)
  inquire(file=tabcodefile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open chosen/default input file (&
 &',tabcodefile(t1:t2),') for FMCSC_TABCODEFILE. Turning off tabulated potentials.'
    tbp%nos = 0
    use_TABUL = .false.
    scale_TABUL = 0.0
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
 78   format('Skipping #1:',i6,' #2:',i6,' w/ TBP#:',i6)
 77   format('Atoms ',i6,' and ',i6,' are in mol. ',i6,', not ',i6,'.')
 76   format('Atom ',i6,' is part of mol. ',i6,', not ',i6,'.')
 74   format('Skipping mode: ',a1,' w/ #1:',i6,' and #2:',i6)
 73   format('Skipping mode: ',a1,' w/ type:',i6)
!
!
  eps = 1000.0*(10.0**(-precision(eps)))
  aone = 1
  atwo = 2
!
  iunit = freeunit()

  open (unit=iunit,file=tabcodefile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for tabulated potentials &
 &input (',tabcodefile(t1:t2),'). Turning off tabulated potential.'
    tbp%nos = 0
    use_TABUL = .false.
    scale_TABUL = 0.0
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing tabulated potentials input: ',tabcodefile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  tbp%mode = -1
  call extract_int(str2,tbp%mode,ii) 
  str3 = str2(ii:MAXSTRLEN)
  szz = -1
  call extract_int(str3,szz,ii)
!
! first list-based formats
  if ((tbp%mode.eq.1).OR.(tbp%mode.eq.2)) then
    if (szz.le.0) then
      write(ilog,*) 'Warning. Second number on header line has to provide number of entries for list-based general&
 & tabulated potentials input (',tabcodefile(t1:t2),'). Turning off term.'
      tbp%nos = 0
      use_TABUL = .false.
      scale_TABUL = 0.0
      close(unit=iunit)
      return
    end if
    tbp%all_lst = szz
    call allocate_tbppre(aone)
!
    tbp%nlst = 0
    do i=1,tbp%all_lst
      read(iunit,79,iostat=iomess) str3
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got empty/incomplete file for tabulated potentials &
 &input (',tabcodefile(t1:t2),'). Turning off term.'
        close(unit=iunit)
        call allocate_tbppre(atwo)
        tbp%nos = 0
        use_TABUL = .false.
        scale_TABUL = 0.0
        return
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing tabulated potentials input (got: ',tabcodefile(t1:t2),').'
        call fexit()
      end if
!
      dum(1) = -1
      dum(2) = -1
      dum(3) = -1
      read(str3,*,iostat=iomess) dum(1),dum(2),dum(3)
!
      if (tbp%mode.eq.1) then
!       the atom-specific format is simple enough
        if ((dum(1).le.0).OR.(dum(1).gt.n).OR.&
 &          (dum(2).le.0).OR.(dum(2).gt.n).OR.(dum(3).le.0)) then
          write(ilog,*) 'Warning. Range exception occurred while reading tabulated potentials input:'
          write(ilog,78) dum(1),dum(2),dum(3)
          cycle
        end if
        tbp%nlst = tbp%nlst + 1
        tbp%lst(tbp%nlst,1) = dum(1)
        tbp%lst(tbp%nlst,2) = dum(2)
        tbp%lst(tbp%nlst,3) = dum(3)
        tbp%lty(tbp%nlst) = .false.
!       more sanity checks are needed for the molecule-type-generalized format
      else if (tbp%mode.eq.2) then
        if ((dum(1).le.0).OR.(dum(1).gt.n).OR.&
 &          (dum(2).le.0).OR.(dum(2).gt.n).OR.(dum(3).le.0)) then
          write(ilog,*) 'Warning. Range exception occurred while reading tabulated potentials input:'
          write(ilog,78) dum(1),dum(2),dum(3)
          cycle
        end if
        imol = molofrs(atmres(dum(1)))
        jmol = molofrs(atmres(dum(2)))
!       the first case is an intramolecular term:
!       this has to be the first mol. of its type
        if (imol.eq.jmol) then
          if (imol.ne.moltyp(moltypid(imol),1)) then
            write(ilog,*) 'Warning. In tabulated potentials, atoms for intramolecular terms have to co&
 &rrespond to the first molecule of their type in sequence:'
            write(ilog,77) dum(1),dum(2),imol,&
 &                           moltyp(moltypid(imol),1)
            cycle
          end if
          tbp%nlst = tbp%nlst + 1
          tbp%lst(tbp%nlst,1) = dum(1)
          tbp%lst(tbp%nlst,2) = dum(2)
          tbp%lst(tbp%nlst,3) = dum(3)
          tbp%lty(tbp%nlst) = .true.
!       the second case is an intermolecular term between different types:
!       both mol.s have to be the first of their respective type
        else if (moltypid(imol).ne.moltypid(jmol)) then
          if ((imol.ne.moltyp(moltypid(imol),1)).OR.&
 &            (jmol.ne.moltyp(moltypid(jmol),1))) then
            write(ilog,*) 'Warning. In tabulated potentials, atoms for cross-intermolecular terms have&
 & to correspond to the first molecules of their types in sequence:'
            if (imol.ne.moltyp(moltypid(imol),1)) then
              write(ilog,76) dum(1),imol,moltyp(moltypid(imol),1)
            else
              write(ilog,76) dum(2),jmol,moltyp(moltypid(jmol),1)
            end if
            cycle
          end if
          tbp%nlst = tbp%nlst + 1
          tbp%lst(tbp%nlst,1) = dum(1)
          tbp%lst(tbp%nlst,2) = dum(2)
          tbp%lst(tbp%nlst,3) = dum(3)
          tbp%lty(tbp%nlst) = .true.
!       the third and final case is an intermolecular term between the same type:
!       the first mol has to be the first of its type, the second the second
        else
          if (imol.ne.moltyp(moltypid(imol),1)) then
            write(ilog,*) 'Warning. In tabulated potentials, atoms for self-intermolecular terms have &
 &to correspond to the first two molecules of their type in sequence:'
            write(ilog,76) dum(1),imol,moltyp(moltypid(imol),1)
            cycle
          end if
          do kmol=1,nmol
            if ((moltypid(kmol).eq.moltypid(imol))&
 &                            .AND.(kmol.ne.imol)) exit
          end do
          if (jmol.ne.kmol) then
            write(ilog,*) 'Warning. In tabulated potentials, atoms for self-intermolecular terms have &
 &to correspond to the first two molecules of their type in sequence:'
            write(ilog,76) dum(2),jmol,kmol
            cycle
          end if
          tbp%nlst = tbp%nlst + 1
          tbp%lst(tbp%nlst,1) = dum(1)
          tbp%lst(tbp%nlst,2) = dum(2)
          tbp%lst(tbp%nlst,3) = dum(3)
          tbp%lty(tbp%nlst) = .true.
        end if
      end if
    end do
!
! second matrix-based format
  else if (tbp%mode.eq.3) then
!
!   first do a dry run to print out all the warnings and to determine the needed memory
    do while (.true.)
      havetogo = .false.
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while processing tabulated potentials input.'
        call fexit()
      end if
      abc(1:1) = ' '
      typ1 = -1
      typ2 = -1
      read(str2,*,iostat=iomess) abc,typ1,typ2
!     intermolecular type-wise requests
      if (abc(1:1).eq.'T') then
        if ((typ1.le.0).OR.(typ1.gt.nmoltyp).OR.&
 &          (typ2.le.0).OR.(typ2.gt.nmoltyp)) then
          write(ilog,*) 'Warning. Range exception occurred while processing tabulated potentials input:'
          write(ilog,74) abc(1:1),typ1,typ2
          cycle
        end if
        if (typ1.eq.typ2) then
          if (moltyp(typ1,2).eq.1) then
            write(ilog,*) 'Warning. In selection of tabulated potentials, intermolecular term for molecule type ',&
 &typ1,' is ignored, as there is only one molecule of this type.'
            cycle
          end if
          do kmol=1,nmol
            if ((moltypid(kmol).eq.typ1).AND.(kmol.ne.moltyp(typ1,1))) exit
          end do
          nm1 = atmol(moltyp(typ1,1),2)-atmol(moltyp(typ1,1),1)+1
          allocate(dummat(nm1,nm1))
!         this is interpreted to be an intermolecular term, but the
!         matrix is symmetric and only half+diag of it is used (has to be
!         fully there, though!)
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
            if (iomess.eq.IOSTAT_END) then
              write(ilog,*) 'Warning. Incomplete matrix found while processing tabulated potentials input:'
              write(ilog,74) abc(1:1),typ1,typ2
              havetogo = .true.
              exit
            else if (iomess.eq.2) then
              write(ilog,*) 'Fatal. I/O error while processing tabulated potentials input.'
              call fexit()
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
!         mirror the matrix
          do ii=1,nm1
            do jj=ii+1,nm1
              dummat(jj,ii) = dummat(ii,jj)
            end do
          end do
          do ii=1,nm1
            do jj=1,nm1
              if (dummat(ii,jj).gt.0) then
                tbp%all_lst = tbp%all_lst + 1
              end if
            end do
          end do
          deallocate(dummat)
        else
          nm1 = atmol(moltyp(typ1,1),2)-atmol(moltyp(typ1,1),1)+1
          nm2 = atmol(moltyp(typ2,1),2)-atmol(moltyp(typ2,1),1)+1
          allocate(dummat(nm1,nm2))
!         this is by necessity an intermolecular term, with a
!         non-symmetric matrix
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm2)
            if (iomess.eq.IOSTAT_END) then
              write(ilog,*) 'Warning. Incomplete matrix found while processing tabulated potentials input:'
              write(ilog,74) abc(1:1),typ1,typ2
              havetogo = .true.
              exit
            else if (iomess.eq.2) then
              write(ilog,*) 'Fatal. I/O error while processing tabulated potentials input.'
              call fexit()
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=1,nm2
              if (dummat(ii,jj).gt.0) then
                tbp%all_lst = tbp%all_lst + 1
              end if
            end do
          end do
          deallocate(dummat)
        end if
!     intramolecular type-wise requests
      else if (abc(1:1).eq.'t') then
        if ((typ1.le.0).OR.(typ1.gt.nmoltyp)) then
          write(ilog,*) 'Warning. Range exception occurred while processing tabulated potentials input:'
          write(ilog,73) abc(1:1),typ1
          cycle
        end if
        nm1 = atmol(moltyp(typ1,1),2)-atmol(moltyp(typ1,1),1)+1
        allocate(dummat(nm1,nm1))
!       as an intramolecular term, the matrix is symmetric and only half sans diag
!       of it is used (has to be fully there, though!)
        do ii=1,nm1
          read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
          if (iomess.eq.IOSTAT_END) then
            write(ilog,*) 'Warning. Incomplete matrix found while processing tabulated potentials input:'
            write(ilog,73) abc(1:1),typ1
            havetogo = .true.
            exit
          else if (iomess.eq.2) then
            write(ilog,*) 'Fatal. I/O error while processing tabulated potentials input.'
            call fexit()
          end if
        end do
        if (havetogo.EQV..true.) then
          deallocate(dummat)
          exit
        end if
        do ii=1,nm1
          do jj=ii+1,nm1
            if (dummat(ii,jj).gt.0) then
              tbp%all_lst = tbp%all_lst + 1
            end if
          end do
        end do
        deallocate(dummat)
      end if
    end do
!
!   now allocate memory, rewind file, and read-off header line
    call allocate_tbppre(aone)
    rewind(unit=iunit)
    read(iunit,79,iostat=iomess) str3
    do while (.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      end if
      abc = ' '
      typ1 = -1
      typ2 = -1
      read(str2,*,iostat=iomess) abc,typ1,typ2
!     intermolecular type-wise requests
      if (abc.eq.'T') then
        if ((typ1.le.0).OR.(typ1.gt.nmoltyp).OR.&
 &          (typ2.le.0).OR.(typ2.gt.nmoltyp)) then
          cycle
        end if
        if (typ1.eq.typ2) then
          if (moltyp(typ1,2).eq.1) then
            cycle
          end if
          do kmol=1,nmol
            if ((moltypid(kmol).eq.typ1).AND.(kmol.ne.moltyp(typ1,1))) exit
          end do
          nm1 = atmol(moltyp(typ1,1),2)-atmol(moltyp(typ1,1),1)+1
          allocate(dummat(nm1,nm1))
!         this is interpreted to be an intermolecular term, but the
!         matrix is symmetric and only half+diag of it is used (has to be
!         fully there, though!)
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
            if (iomess.eq.IOSTAT_END) then
              havetogo = .true.
              exit
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=ii+1,nm1
              dummat(jj,ii) = dummat(ii,jj)
            end do
          end do
          do ii=1,nm1
            do jj=1,nm1
              if (dummat(ii,jj).gt.0) then
                tbp%nlst = tbp%nlst + 1
                tbp%lst(tbp%nlst,1) = atmol(moltyp(typ1,1),1)+ii-1
                tbp%lst(tbp%nlst,2) = atmol(kmol,1)+jj-1
                tbp%lst(tbp%nlst,3) = dummat(ii,jj)
                tbp%lty(tbp%nlst) = .true.
              end if
            end do
          end do
          deallocate(dummat)
        else
          nm1 = atmol(moltyp(typ1,1),2)-atmol(moltyp(typ1,1),1)+1
          nm2 = atmol(moltyp(typ2,1),2)-atmol(moltyp(typ2,1),1)+1
          allocate(dummat(nm1,nm2))
!         this is by necessity an intermolecular term, with a
!         non-symmetric matrix
          do ii=1,nm1
            read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm2)
            if (iomess.eq.IOSTAT_END) then
              havetogo = .true.
              exit
            end if
          end do
          if (havetogo.EQV..true.) then
            deallocate(dummat)
            exit
          end if
          do ii=1,nm1
            do jj=1,nm2
              if (dummat(ii,jj).gt.0) then
                tbp%nlst = tbp%nlst + 1
                tbp%lst(tbp%nlst,1) = atmol(moltyp(typ1,1),1)+ii-1
                tbp%lst(tbp%nlst,2) = atmol(moltyp(typ2,1),1)+jj-1
                tbp%lst(tbp%nlst,3) = dummat(ii,jj)
                tbp%lty(tbp%nlst) = .true.
              end if
            end do
          end do
          deallocate(dummat)
        end if
!     intramolecular type-wise requests
      else if (abc(1:1).eq.'t') then
        if ((typ1.le.0).OR.(typ1.gt.nmoltyp)) then
          cycle
        end if
        nm1 = atmol(moltyp(typ1,1),2)-atmol(moltyp(typ1,1),1)+1
        allocate(dummat(nm1,nm1))
        do ii=1,nm1
          read(iunit,*,iostat=iomess) (dummat(ii,jj),jj=1,nm1)
          if (iomess.eq.IOSTAT_END) then
            havetogo = .true.
            exit
          end if
        end do
        if (havetogo.EQV..true.) then
          deallocate(dummat)
          exit
        end if
        do ii=1,nm1
          do jj=ii+1,nm1
            if (dummat(ii,jj).gt.0) then
              tbp%nlst = tbp%nlst + 1
              tbp%lst(tbp%nlst,1) = atmol(moltyp(typ1,1),1)+ii-1
              tbp%lst(tbp%nlst,2) = atmol(moltyp(typ1,1),1)+jj-1
              tbp%lst(tbp%nlst,3) = dummat(ii,jj)
              tbp%lty(tbp%nlst) = .true.
            end if
          end do
        end do
        deallocate(dummat)
      end if
    end do
!
  else
    write(ilog,*) 'Warning. Unknown mode for tabulated potentials input (identifier: ', &
 &tbp%mode,'). Turning term off.'
    tbp%nos = 0
    use_TABUL = .false.
    scale_TABUL = 0.0
    close(unit=iunit)
    return
  end if
!
  tbp%nos = 0
  allocate(pclst(tbp%all_lst))
  do i=1,tbp%nlst
    novel = .true.
    do j=1,tbp%nos
      if (tbp%lst(i,3).eq.pclst(j)) then
        novel = .false.
        exit
      end if
    end do
    if (novel.EQV..true.) then
      tbp%nos = tbp%nos + 1
      pclst(tbp%nos) = tbp%lst(i,3)
    end if
  end do
  if (tbp%nos.eq.0) then
    write(ilog,*) 'Warning. No applicable tabulated potentials appear to be requested in ',&
 &tabcodefile(t1:t2),'. Turning energy term off.'
    deallocate(pclst)
    call allocate_tbppre(atwo)
    use_TABUL = .false.
    scale_TABUL = 0.0
    close(unit=iunit)
    return
  end if
  call sort(tbp%nos,pclst)
!
  do j=1,tbp%nos
    if (pclst(j).ne.j) then
      write(ilog,*) 'Warning. Discovered discontinous numbering of tabulated &
 &potentials while processing corresponding input. Turning off term.'
      deallocate(pclst)
      tbp%nos = 0
      call allocate_tbppre(atwo)
      use_TABUL = .false.
      scale_TABUL = 0.0
      close(unit=iunit)
      return
    end if
  end do
!
  deallocate(pclst)
  close(unit=iunit)
!
! and now the potential input
!
! first allocate a dummy array for the dry run
  allocate(potrs(tbp%nos+1))
! open the file
  iunit = freeunit()
  call strlims(tabpotfile,t1,t2)

  inquire(file=tabpotfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open chosen input file (&
 &',tabpotfile(t1:t2),') for FMCSC_TABPOTFILE. Turning off tabulat&
 &ed potentials.'
    tbp%nos = 0
    use_TABUL = .false.
    scale_TABUL = 0.0
    return
  end if
!
  open (unit=iunit,file=tabpotfile(t1:t2),status='old')
!
  i = 0
  do while (1.eq.1)
    i = i + 1
    read(iunit,*,iostat=iomess) potrs(1),(potrs(k+1),k=1,tbp%nos)
    if (iomess.eq.IOSTAT_END) then
      i = i - 1
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. I/O error while processing tabulated potentials input.'
      call fexit()
    end if
    if (potrs(1).lt.0.0) then
      write(ilog,*) 'Warning. Tabulated potentials have negative &
 &distances (in ',tabpotfile(t1:t2),'). Turning off term.'
      deallocate(potrs)
      tbp%nos = 0
      call allocate_tbppre(atwo)
      use_TABUL = .false.
      scale_TABUL = 0.0
      return
    end if
    if (i.eq.2) disres = potrs(1) - disold
    if (i.gt.2) then
      if (abs(disres-potrs(1)+disold).gt.eps) then
        write(ilog,*) 'Warning. Tabulated potentials have variable d&
 &istance intervals (in ',tabpotfile(t1:t2),'). Turning off term.'
        deallocate(potrs)
        tbp%nos = 0
        call allocate_tbppre(atwo)
        use_TABUL = .false.
        scale_TABUL = 0.0
        return
      end if
    end if
    if (i.gt.1) then
      if (disres.le.0.0) then
        write(ilog,*) 'Warning. Tabulated potentials have disallowed d&
 &istance column (in ',tabpotfile(t1:t2),'). Turning off term.'
        deallocate(potrs)
        tbp%nos = 0
        call allocate_tbppre(atwo)
        use_TABUL = .false.
        scale_TABUL = 0.0
        return
      end if
    end if
    disold = potrs(1)
  end do
!
  if (i.le.1) then
    write(ilog,*) 'Warning. Got empty/incomplete file for tabulated poten&
 &tials (',tabpotfile(t1:t2),'). Turning off term.'
    tbp%nos = 0
    call allocate_tbppre(atwo)
    use_TABUL = .false.
    scale_TABUL = 0.0
    deallocate(potrs)
    return
  end if
!
! finally set the parameters necessary for allocation
  tbp%bins = i
  tbp%res = disres
  call allocate_tabpot(aone)
!
! rewind file and read in potentials
  rewind(unit=iunit)
  do i=1,tbp%bins
    read(iunit,*,iostat=iomess) tbp%dis(i),&
 &                      (tbp%pot(k,i),k=1,tbp%nos)
  end do
!
  close(unit=iunit)
!
! see if the user also provided the tangents 
  iunit = freeunit()
  call strlims(tabtangfile,t1,t2)

  inquire(file=tabtangfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Warning. Cannot open optional input file (&
 &',tabtangfile(t1:t2),') for FMCSC_TABTANGFILE. This means tangents are &
 &constructed from potential data according to spline settings.'
    tbp%tangmode = 0
    j = 1
!
  else
!
    open (unit=iunit,file=tabtangfile(t1:t2),status='old')
!
    i = 0
    do while (1.eq.1)
      i = i + 1
      read(iunit,*,iostat=iomess) (potrs(k),k=1,tbp%nos)
      if (iomess.eq.IOSTAT_END) then
        i = i - 1
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. I/O error while processing tangents input for tabulated potentials.'
        call fexit()
      end if
      tbp%tang(1:tbp%nos,i) = potrs(1:tbp%nos)
    end do
!
    tbp%tangmode = 1
    if (i.lt.tbp%bins) then
      write(ilog,*) 'Warning. Got empty/incomplete file for tangents for tabulated poten&
 &tials (',tabtangfile(t1:t2),'). Remaining values will be constructed from spline settings.'
      tbp%tangmode = 2
    end if
!
    j = i + 1
!
    close(unit=iunit)
!
  end if
!
  deallocate(potrs)
!
  if (j.le.tbp%bins) then
    fac(1) = (1.0-tbp%ipprm(1))*(1.0+tbp%ipprm(2))
    fac(2) = (1.0-tbp%ipprm(1))*(1.0-tbp%ipprm(2))
!
!   set the (remaining) tangents based on interpolation scheme (these are simplified from Kochanek-Bartels splines lacking
!   the harmful continuity parameter)
    if (j.eq.1) then
      do k=1,tbp%nos
!       approximate end-point to decay to flatness (alt. as line extension)
        tbp%tang(k,1) = 0.0 ! (tbp%pot(k,2)-tbp%pot(k,1))/(tbp%dis(2)-tbp%dis(1))
      end do
      j = 2
    end if
    do k=1,tbp%nos
!     approximate end-point to decay to flatness (alt. as line extension)
      tbp%tang(k,tbp%bins) = 0.0 ! (tbp%pot(k,tbp%bins)-tbp%pot(k,tbp%bins-1))/(tbp%dis(tbp%bins)-tbp%dis(tbp%bins-1))
    end do
!   this doesn't respond correctly to irregularly spaced bins
    if (j.lt.tbp%bins) then
      do i=j,tbp%bins-1
        do k=1,tbp%nos
          tbp%tang(k,i) = (fac(1)*(tbp%pot(k,i)-tbp%pot(k,i-1)) + fac(2)*(tbp%pot(k,i+1)-tbp%pot(k,i)))/&
 &                         (tbp%dis(i+1)-tbp%dis(i-1))
        end do
      end do
    end if
  end if
!
! initialize res.-res. matrix
  do rs1=1,nseq
    do rs2=1,nseq
      tbp%rsmat(rs1,rs2) = 0
    end do
    tbp%rsvec(rs1) = 0
  end do
!
! make sure first residue is always first
  do i=1,tbp%nlst
    if (atmres(tbp%lst(i,1)).gt.atmres(tbp%lst(i,2))) then
      k = tbp%lst(i,1)
      tbp%lst(i,1) = tbp%lst(i,2)
      tbp%lst(i,2) = k
    end if
  end do
!
! now sort first residue column
  allocate(potrs(3))
  do i=1,tbp%nlst
    rs1 = atmres(tbp%lst(i,1))
    ii = i
    do j=i+1,tbp%nlst
      rs2 = atmres(tbp%lst(j,1))
      if (rs2.lt.rs1) then
        ii = j
        rs1 = atmres(tbp%lst(ii,1))
      end if
    end do
    if (ii.ne.i) then
      do k=1,3
        potrs(k) = tbp%lst(i,k)
        tbp%lst(i,k) = tbp%lst(ii,k)
        tbp%lst(ii,k) = potrs(k)
      end do
    end if
  end do
! and second one (conditionally)
  do i=1,tbp%nlst
    rs1 = atmres(tbp%lst(i,2))
    ii = i
    do j=i+1,tbp%nlst
      rs2 = atmres(tbp%lst(j,2))
      if (atmres(tbp%lst(j,1)).eq.atmres(tbp%lst(i,1))) then
        if (rs2.lt.rs1) then
          ii = j
          rs1 = atmres(tbp%lst(ii,2))
        end if
      end if
    end do
    if (ii.ne.i) then
      do k=1,3
        potrs(k) = tbp%lst(i,k)
        tbp%lst(i,k) = tbp%lst(ii,k)
        tbp%lst(ii,k) = potrs(k)
      end do
    end if
  end do
!
! 444 format(10(i8,1x))
!  do i=1,tbp%nlst
!    write(*,444) i,(tbp%lst(i,k),k=1,3),atmres(tbp%lst(i,1)),
! &           atmres(tbp%lst(i,2))
!  end do
  deallocate(potrs)           
!
! now we can set the bounds for each residue-residue pair
  rs1o = 0
  rs2o = 0
  do i=1,tbp%nlst
    rs1 = atmres(tbp%lst(i,1))
    rs2 = atmres(tbp%lst(i,2))
    if (tbp%lty(i).EQV..false.) then
      if ((rs2.ne.rs2o).OR.(rs1.ne.rs1o)) then
        tbp%rsmat(rs1,rs2) = i
        if ((rs2o.gt.0).AND.(rs1o.gt.0)) then
          if (rs1o.ne.rs2o) then
            tbp%rsmat(rs2o,rs1o) = i-1
          else
            tbp%rsvec(rs1o) = i-1
          end if
        end if
      end if
      if (i.eq.tbp%nlst) then
        if (rs2.ne.rs1) then
          tbp%rsmat(rs2,rs1) = i
        else
          tbp%rsvec(rs1) = i
        end if
      end if
    else
      typ1 = moltypid(molofrs(rs1))
      typ2 = moltypid(molofrs(rs2))
      off1 = rs1 - rsmol(molofrs(rs1),1)
      off2 = rs2 - rsmol(molofrs(rs2),1)
      if ((rs2.ne.rs2o).OR.(rs1.ne.rs1o)) then
        do imol=1,nmol
          if (molofrs(rs1).eq.molofrs(rs2)) then
            if (moltypid(imol).eq.typ1) then
              tbp%rsmat(rsmol(imol,1)+off1,rsmol(imol,1)+off2) = i
            end if
          else
            do jmol=imol+1,nmol
              if ((moltypid(imol).eq.typ1).AND.&
 &                   (moltypid(jmol).eq.typ2)) then
                tbp%rsmat(rsmol(imol,1)+off1,rsmol(jmol,1)+off2) = i
              else if ((moltypid(jmol).eq.typ1).AND.&
 &                   (moltypid(imol).eq.typ2)) then
                tbp%rsmat(rsmol(jmol,1)+off2,rsmol(imol,1)+off1) = i
              end if
            end do
          end if
        end do
        if ((rs2o.gt.0).AND.(rs1o.gt.0)) then
          oo1 = rs1o - rsmol(molofrs(rs1o),1)
          oo2 = rs2o - rsmol(molofrs(rs2o),1)
          do imol=1,nmol
            if (molofrs(rs1o).eq.molofrs(rs2o)) then
              if (moltypid(imol).eq.moltypid(molofrs(rs1o))) then
                if (rs1o.ne.rs2o) then
                  tbp%rsmat(rsmol(imol,1)+oo2,rsmol(imol,1)+oo1)=i-1
                else
                  tbp%rsvec(rsmol(imol,1)+oo1) = i-1
                end if
              end if
            else
              do jmol=imol+1,nmol
                if ((moltypid(imol).eq.moltypid(molofrs(rs1o))).AND.&
 &              (moltypid(jmol).eq.moltypid(molofrs(rs2o)))) then
                tbp%rsmat(rsmol(jmol,1)+oo2,rsmol(imol,1)+oo1) = i-1
                else if ((moltypid(jmol).eq.moltypid(molofrs(rs2o)))&
 &      .AND.(moltypid(imol).eq.moltypid(molofrs(rs1o)))) then
                  tbp%rsmat(rsmol(imol,1)+oo1,rsmol(jmol,1)+oo2)=i-1
                end if
              end do
            end if
          end do
        end if
      end if
      if (i.eq.tbp%nlst) then
        do imol=1,nmol
          if (molofrs(rs1).eq.molofrs(rs2)) then
            if (moltypid(imol).eq.typ1) then
              if (rs1.ne.rs2) then
                tbp%rsmat(rsmol(imol,1)+off2,rsmol(imol,1)+off1)= i
              else
                tbp%rsvec(rsmol(imol,1)+off1) = i
              end if
            end if
          else
            do jmol=imol+1,nmol
              if ((moltypid(imol).eq.typ1).AND.&
 &                   (moltypid(jmol).eq.typ2)) then
                tbp%rsmat(rsmol(jmol,1)+off2,rsmol(imol,1)+off1)= i
              else if ((moltypid(jmol).eq.typ1).AND.&
 &                   (moltypid(imol).eq.typ2)) then
                tbp%rsmat(rsmol(imol,1)+off1,rsmol(jmol,1)+off2)= i
              end if
            end do
          end if
        end do
      end if
    end if
    rs1o = rs1
    rs2o = rs2
  end do
!
!  do rs1=1,nseq
!    if (tbp%rsmat(rs1,rs1).gt.0) then
!        write(*,*) rs1,rs1,'start ',tbp%rsmat(rs1,rs1),' and end',tbp%rsvec(rs1)
!    end if
!    do rs2=rs1+1,nseq
!      if (tbp%rsmat(rs1,rs2).gt.0) then
!        write(*,*) rs1,rs2,'start ',tbp%rsmat(rs1,rs2),' and end',tbp%rsmat(rs2,rs1)
!      end if
!      if (((tbp%rsmat(rs1,rs2).gt.0).AND.(tbp%rsmat(rs2,rs1).eq.0)).OR.&
! &((tbp%rsmat(rs1,rs2).eq.0).AND.(tbp%rsmat(rs2,rs1).gt.0)))then
!        write(*,*) 'Fatal for ',rs1,rs2
!      end if
!    end do
!  end do
!
!
  if (tabul_report.EQV..true.) then
#ifdef ENABLE_MPI
    if ((use_REMC.EQV..true.).AND.(myrank.eq.0)) then
      str3 = 'N_000_TABULATED_POT.idx'
    else if ((use_MPIAVG.EQV..true.).AND.(myrank.eq.0)) then
      str3 = 'TABULATED_POT.idx'
    else
      return
    end if
#else
    str3 = 'TABULATED_POT.idx'
#endif
    call strlims(str3,t1,t2)
!
   99   format('--- Tabulated Distance-Dependent Potential #',i5)
   98   format('Specific: Atom ',i6,' (',a4,') in Res. ',i6,' (',a3,')')
   97   format('          Atom ',i6,' (',a4,') in Res. ',i6,' (',a3,')')
   96   format('Type-wise and |',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &                                 ') of Mol. Type ',i4)
   95   format('Intramolecular|',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &')')
   94   format('Type-wise and |',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &                                 ') of Mol. Type ',i4)
   93   format('Intermolecular|',6x,' (',a4,') in Res. ',i6,' (',a3,&
 &                                 ') of Mol. Type ',i4)
!
    inquire(file=str3(t1:t2),exist=exists)
    if(exists) then
      iunit = freeunit()
      open(unit=iunit,file=str3(t1:t2),status='old')
      close(unit=iunit,status='delete')
    end if
    iunit=freeunit()
    open(unit=iunit,file=str3(t1:t2),status='new')
!
    write(iunit,*) 'List of Terms for Columns in Tabulated Potential&
 & Input'
    write(iunit,*) '------------------------------------------------&
 &------'
    write(iunit,*)
    allocate(pclst(tbp%nlst))
    pcnlst = 0
    do i=1,tbp%nos
      do j=1,tbp%nlst
        if (tbp%lst(j,3).eq.i) then
          pcnlst = pcnlst + 1
          pclst(pcnlst) = j
        end if
      end do
!
      write(iunit,99) i
      do j=1,pcnlst
        k = pclst(j)
        ati = tbp%lst(k,1)
        atj = tbp%lst(k,2)
        if (tbp%lty(k).EQV..false.) then
          write(iunit,98) ati,bio_code(b_type(ati)),atmres(ati),&
 &                     amino(seqtyp(atmres(ati)))
          write(iunit,97) atj,bio_code(b_type(atj)),atmres(atj),&
 &                     amino(seqtyp(atmres(atj)))
        else
          imol = molofrs(atmres(ati))
          jmol = molofrs(atmres(atj))
          ii = atmres(ati)
          jj = atmres(atj)
          if (imol.eq.jmol) then
            write(iunit,96)bio_code(b_type(ati)),ii-rsmol(imol,1)+1,&
 &amino(seqtyp(ii)),moltypid(imol)
            write(iunit,95)bio_code(b_type(atj)),jj-rsmol(jmol,1)+1,&
 &amino(seqtyp(jj))
          else
            write(iunit,94)bio_code(b_type(ati)),ii-rsmol(imol,1)+1,&
 &amino(seqtyp(ii)),moltypid(imol)
            write(iunit,93)bio_code(b_type(atj)),jj-rsmol(jmol,1)+1,&
 &amino(seqtyp(jj)),moltypid(jmol)
          end if
        end if
      end do
      write(iunit,*)
      pcnlst = 0
    end do
    deallocate(pclst)
    close(unit=iunit)
  end if
!
!
!
!  if (tabul_report.EQV..true.) then
!   write the rs-rs-matrix
! 77     format(1000(i1,1x))
! 78     format('The last informative distance bin seems to be #',i3)
! 80     format('This is a total number of ',i8,' interactions.')
! 79     format('Atoms ',i6,' and ',i6,' through Pot. # ',i5)
!    write(ilog,*)
!    write(ilog,*) '--- Summary of Tabulated Terms in the System ---'
!    write(ilog,*)
!    write(ilog,*) '--- Residue-Residue Matrix (any pair at all) ---'
!    write(ilog,*) 
!    do i=1,nseq
!      write(ilog,77) (tab_rsmat(i,j),j=1,nseq)
!    end do
!    write(ilog,*)
!    write(ilog,78) dummy
!    tt = 0 
!    do ii=1,par_TABUL(3)
!      do kk=ii+1,par_TABUL(3)
!        if (map_tab(ii,kk).gt.0) then
!          tt = tt + 1
!          write(ilog,79) ii,kk,map_tab(ii,kk)
!        end if
!      end do
!    end do
!    write(ilog,*)
!    write(ilog,80) tt 
!    write(ilog,*)
!  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine read_frzfile()
!
  use energies
  use iounit
  use molecule
  use sequen
  use movesets
  use aminos
  use fyoc
  use system
  use forces
  use interfaces
  use zmatrix
  use atoms
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,ii,j,jj,imol,mn,rs
  integer mt,aone,atwo,athree,frzidx,ncons(50)
  integer, ALLOCATABLE:: rslst(:)
  character(3) dofs
  character(10) abc
  character(MAXSTRLEN) str2
  logical exists,sayno,sayyes
!
  aone = 1
  atwo = 2
  athree = 3
  sayno = .false.
  sayyes = .true.
!
  call strlims(frzfile,t1,t2)
  inquire(file=frzfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &frzfile(t1:t2),') for FMCSC_FRZFILE. All eligible d.o.f. will be sampled.'
    do_frz = .false.
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=frzfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for constraint&
 & input (',frzfile(t1:t2),'). All eligible d.o.f. will be sampled.'
    do_frz = .false.
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing constraint&
 & input (got: ',frzfile(t1:t2),').'
    call fexit()
  end if
!
  allocate(rslst(nseq))
!
  ii = 1
  call extract_str(str2,abc,ii)
!
! first bail-out for Cartesian dynamics -> such constraints are not supported via this path
  if ((fycxyz.eq.2).AND.(.NOT.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))).AND.&
 &    (abc(1:1).ne.'a').AND.(abc(1:1).ne.'A')) then
    write(ilog,*) 'Warning. Only atom-specific constraint requests (mode "A") are supported in FMCSC_FRZFILE if&
  & the run uses a pure Cartesian dynamics treatment (also consider holonomic constraints -> FMCSC_SHAKESET).'
    do_frz = .false.
    return
  else if ((fycxyz.eq.1).AND.(dyn_mode.eq.1).AND.((abc(1:1).eq.'a').OR.(abc(1:1).eq.'A'))) then
    write(ilog,*) 'Warning. Atom-specific constraint requests (mode "A") via FMCSC_FRZFILE are not supported if&
  & the run is a pure Monte Carlo run in internal coordinate space.'
    do_frz = .false.
    return
  else if ((fycxyz.eq.2).AND.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
    write(ilog,*) 'Warning. In hybrid runs with Cartesian space sampling for the dynamics, the request of any type of &
 &constraints via FMCSC_FRZFILE implies that these constraints will apply exclusively to either internal coordinates &
 &(MC portion) or atomic positions (dynamics portion). This is generally a BAD idea (inconsistencies). The only workaround &
 &would be the presence of additional, holonomic constraints (-> FMCSC_SHAKESET) to match the MC constraint set exactly.'
  else if ((fycxyz.eq.1).AND.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
    if ((abc(1:1).eq.'a').OR.(abc(1:1).eq.'A')) then
      write(ilog,*) 'In hybrid internal coordinate space simulations, atom-specific constraint requests (mode "A") will only &
 &be active for the dynamics portion (use a different mode to create a consistent set of constraints).'
    end if
  end if
!
  ncons(:) = 0 
!
! first molecule-based
!
  if ((abc(1:1).eq.'M').OR.(abc(1:1).eq.'m')) then
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing constr&
 &aint input (got: ',frzfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      mn = -1
      call extract_int(str2,mn,ii)
      if ((mn.le.0).OR.(mn.gt.nmol)) then
        write(ilog,*) 'Requested constraints for molecule ',&
 &mn,', which does not exist. Skipping ...'
        cycle
      end if
      if (ii.eq.1) exit
      call extract_abc(str2,dofs,ii)
      call toupper(dofs)
!
      if (dofs.eq.'INT') then
        if (do_pol(moltypid(mn)).EQV..false.) then
          write(ilog,*) 'Warning: Molecule ',mn,' has no internal de&
 &grees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(aone,mn,sayyes,sayyes,sayyes,sayyes,sayyes,rslst,ncons)
        end if
      else if (dofs.eq.'RBC') then
        call remove_rbcdof(aone,mn,ncons)
      else if (dofs.eq.'ALL') then
        if (do_pol(moltypid(mn)).EQV..false.) then
          write(ilog,*) 'Warning: Molecule ',mn,' has no internal de&
 &grees of freedom -> just constraining RBC.'
        else
          call remove_intdof(aone,mn,sayyes,sayyes,sayyes,sayyes,sayyes,rslst,ncons)
        end if
        call remove_rbcdof(aone,mn,ncons)
      else if (dofs.eq.'CHI') then
        if (do_pol(moltypid(mn)).EQV..false.) then
          write(ilog,*) 'Warning: Molecule ',mn,' has no internal de&
 &grees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(aone,mn,sayno,sayno,sayyes,sayno,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'FYO') then
        if (do_pol(moltypid(mn)).EQV..false.) then
          write(ilog,*) 'Warning: Molecule ',mn,' has no internal de&
 &grees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(aone,mn,sayyes,sayyes,sayno,sayno,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'PUC') then
        if (do_pol(moltypid(mn)).EQV..false.) then
          write(ilog,*) 'Warning: Molecule ',mn,' has no internal de&
 &grees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(aone,mn,sayno,sayno,sayno,sayyes,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'OTH') then
        if (do_pol(moltypid(mn)).EQV..false.) then
          write(ilog,*) 'Warning: Molecule ',mn,' has no internal de&
 &grees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(aone,mn,sayno,sayno,sayno,sayno,sayyes,rslst,ncons)
        end if
      else
        write(ilog,*) 'Unknown selection for constraint request for &
 &molecule ',mn,' in ',frzfile(t1:t2),': ',dofs,' - ignored.'
      end if
!
    end do
!
! and now molecule-type-based
!
  else if ((abc(1:1).eq.'T').OR.(abc(1:1).eq.'t')) then
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing constr&
 &aint input (got: ',frzfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      mt = -1
      call extract_int(str2,mt,ii)
      if ((mt.le.0).OR.(mt.gt.nmoltyp)) then
        write(ilog,*) 'Fatal. Requested constraints for molecule typ&
 &e ',mt,', which does not exist. Skipping ...'
        cycle
      end if
      if (ii.eq.1) exit
      call extract_abc(str2,dofs,ii)
      call toupper(dofs)
!
      if (dofs.eq.'INT') then
        if (do_pol(mt).EQV..false.) then
          write(ilog,*) 'Warning: Molecule type ',mt,' has no intern&
 &al degrees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(atwo,mt,sayyes,sayyes,sayyes,sayyes,sayyes,rslst,ncons)
        end if
      else if (dofs.eq.'RBC') then
        call remove_rbcdof(atwo,mt,ncons)
      else if (dofs.eq.'ALL') then
        if (do_pol(mt).EQV..false.) then
          write(ilog,*) 'Warning: Molecule type ',mt,' has no intern&
 &al degrees of freedom -> just constraining RBC.'
        else
          call remove_intdof(atwo,mt,sayyes,sayyes,sayyes,sayyes,sayyes,rslst,ncons)
        end if
        call remove_rbcdof(atwo,mt,ncons)
      else if (dofs.eq.'CHI') then
        if (do_pol(mt).EQV..false.) then
          write(ilog,*) 'Warning: Molecule type ',mt,' has no intern&
 &al degrees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(atwo,mt,sayno,sayno,sayyes,sayno,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'FYO') then
        if (do_pol(mt).EQV..false.) then
          write(ilog,*) 'Warning: Molecule type ',mt,' has no intern&
 &al degrees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(atwo,mt,sayyes,sayyes,sayno,sayno,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'PUC') then
        if (do_pol(mt).EQV..false.) then
          write(ilog,*) 'Warning: Molecule type ',mt,' has no intern&
 &al degrees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(atwo,mt,sayno,sayno,sayno,sayyes,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'OTH') then
        if (do_pol(mt).EQV..false.) then
          write(ilog,*) 'Warning: Molecule type ',mt,' has no intern&
 &al degrees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(atwo,mt,sayno,sayno,sayno,sayno,sayyes,rslst,ncons)
        end if
      else
        write(ilog,*) 'Unknown selection for constraint request for &
 &molecule type ',mt,' in ',frzfile(t1:t2),': ',dofs,' - ignored.'
      end if
!
    end do
!
! residue-based (ugly for RBC)
  else if ((abc(1:1).eq.'R').OR.(abc(1:1).eq.'r')) then
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing constr&
 &aint input (got: ',frzfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      rs = -1
      call extract_int(str2,rs,ii)
      if ((rs.le.0).OR.(rs.gt.nseq)) then
        write(ilog,*) 'Requested constraints for residue ',&
 &rs,', which does not exist. Skipping ...'
        cycle
      end if
      if (ii.eq.1) exit
      call extract_abc(str2,dofs,ii)
      call toupper(dofs)
!
      if (dofs.eq.'INT') then
        if (notors(rs).EQV..true.) then
          write(ilog,*) 'Warning: Residue ',rs,' has no internal deg&
 &rees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(athree,rs,sayyes,sayyes,sayyes,sayyes,sayyes,rslst,ncons)
        end if
      else if (dofs.eq.'RBC') then
        write(ilog,*) 'Warning: Rigid-body constraints for residues &
 &will fixate the whole molecule they belong to(!).'
        call remove_rbcdof(aone,molofrs(rs),ncons)
      else if (dofs.eq.'ALL') then
        if (notors(rs).EQV..true.) then
          write(ilog,*) 'Warning: Residue ',rs,' has no internal deg&
 &rees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(athree,rs,sayyes,sayyes,sayyes,sayyes,sayyes,rslst,ncons)
        end if
        write(ilog,*) 'Warning: Rigid-body constraints for residues &
 &will fixate the whole molecule they belong to(!).'
        call remove_rbcdof(aone,molofrs(rs),ncons)
      else if (dofs.eq.'CHI') then
        if (nchi(rs).le.0) then
          write(ilog,*) 'Warning: Residue ',rs,' has no sidechain deg&
 &rees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(athree,rs,sayno,sayno,sayyes,sayno,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'FYO') then
        if (notors(rs).EQV..true.) then
          write(ilog,*) 'Warning: Residue ',rs,' has no internal deg&
 &rees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(athree,rs,sayyes,sayyes,sayno,sayno,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'PUC') then
        if (.NOT.((pucline(rs).gt.0).OR.(seqflag(rs).eq.5).OR.(seqpolty(rs).eq.'N'))) then
          write(ilog,*) 'Warning: Residue ',rs,' has no pucker deg&
 &rees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(athree,rs,sayno,sayno,sayno,sayyes,sayno,rslst,ncons)
        end if
      else if (dofs.eq.'OTH') then
        if (notors(rs).EQV..true.) then
          write(ilog,*) 'Warning: Residue ',rs,' has no internal deg&
 &rees of freedom, which could be constrained -> ignored.'
        else
          call remove_intdof(athree,rs,sayno,sayno,sayno,sayno,sayyes,rslst,ncons)
        end if
      else
        write(ilog,*) 'Unknown selection for constraint request for &
 &residue',rs,' in ',frzfile(t1:t2),': ',dofs,' - ignored.'
      end if
!
    end do
!
! atom-based (the xyz option is only available if we have Cartesian space sampling, the INT option only for
! TMD-containing runs)
  else if ((abc(1:1).eq.'A').OR.(abc(1:1).eq.'a')) then
    do while (.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing constr&
 &aint input (got: ',frzfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      rs = -1
      call extract_int(str2,rs,ii)
      if ((rs.le.0).OR.(rs.gt.n)) then
        write(ilog,*) 'Requested constraints for atom ',rs,', which does not exist. Skipping ...'
        cycle
      end if
      if (ii.eq.1) exit
      call extract_abc(str2,dofs,ii)
      call toupper(dofs)
!
      if (fycxyz.eq.2) then
        if (dofs.eq.'XYZ') then
          do j=1,3
            if (cart_frz(rs,j).EQV..false.) then 
              cart_frz(rs,j) = .true.
              n_constraints = n_constraints + 1
              ncons(23+j) = ncons(23+j) + 1
            end if
          end do
        else if (dofs(1:2).eq.'XY') then
          do j=1,2
            if (cart_frz(rs,j).EQV..false.) then 
              cart_frz(rs,j) = .true.
              n_constraints = n_constraints + 1
              ncons(23+j) = ncons(23+j) + 1
            end if
          end do
        else if (dofs(1:2).eq.'XZ') then
          do jj=3,4
            j = jj
            if (jj.eq.4) j = 1
            if (cart_frz(rs,j).EQV..false.) then 
              cart_frz(rs,j) = .true.
              n_constraints = n_constraints + 1
              ncons(23+j) = ncons(23+j) + 1
            end if
          end do
        else if (dofs(1:2).eq.'YZ') then
          do j=2,3
            if (cart_frz(rs,j).EQV..false.) then 
              cart_frz(rs,j) = .true.
              n_constraints = n_constraints + 1
              ncons(23+j) = ncons(23+j) + 1
            end if
          end do
        else if (dofs(1:1).eq.'X') then
          if (cart_frz(rs,1).EQV..false.) then
            cart_frz(rs,1) = .true.
            n_constraints = n_constraints + 1
            ncons(24) = ncons(24) + 1
          end if
        else if (dofs(1:1).eq.'Y') then
          if (cart_frz(rs,2).EQV..false.) then
            cart_frz(rs,2) = .true.
            n_constraints = n_constraints + 1
            ncons(25) = ncons(25) + 1
          end if
        else if (dofs(1:1).eq.'Z') then
          if (cart_frz(rs,3).EQV..false.) then
            cart_frz(rs,3) = .true.
            n_constraints = n_constraints + 1
            ncons(26) = ncons(26) + 1
          end if
        else
          write(ilog,*) 'Unknown selection for Cartesian space constraint request for &
   &atom ',rs,' in ',frzfile(t1:t2),': ',dofs,' - ignored.'
        end if
      else if (use_dyn.EQV..true.) then
        if (dofs.eq.'INT') then
          imol = molofrs(atmres(rs))
          if ((izrot(rs)%alsz.le.0).OR.(rs.eq.atmol(imol,1))) then 
            write(ilog,*) 'Warning. Requested dihedral angle constraint for Z-matrix line ',rs,&
 &', which does not correspond to a rotatable dihedral angle. Ignored.'
            cycle
          else if ((izrot(rs)%treevs(4).le.0).OR.(izrot(rs)%treevs(4).gt.dc_di(imol)%maxntor))  then
            write(ilog,*) 'Fatal. Requested dihedral angle constraint for Z-matrix line ',rs,&
 &', which is eligible, but has no rank. This is most likely a bug.'
            call fexit()
          else if (dc_di(imol)%recurs(izrot(rs)%treevs(4),3).le.0)  then
            write(ilog,*) 'Fatal. Requested dihedral angle constraint for Z-matrix line ',rs,&
 &', which is eligible, but not set up in the recursive structure. This is most likely a bug.'
            call fexit()
          end if
          if (dc_di(imol)%frz(dc_di(imol)%recurs(izrot(rs)%treevs(4),3)).EQV..false.) then
            dc_di(imol)%frz(dc_di(imol)%recurs(izrot(rs)%treevs(4),3)) = .true.
            n_constraints = n_constraints + 1
            ncons(27) = ncons(27) + 1
          end if
        else if ((dofs.eq.'XYZ').OR.(dofs.eq.'ABC').OR.(dofs(1:2).eq.'XY').OR.(dofs(1:2).eq.'XZ').OR.&
 &   (dofs(1:2).eq.'YZ').OR.(dofs(1:2).eq.'AB').OR.(dofs(1:2).eq.'BC').OR.(dofs(1:2).eq.'AC').OR.(dofs(1:1).eq.'X').OR.&
 &   (dofs(1:1).eq.'Y').OR.(dofs(1:1).eq.'Z').OR.((dofs(1:1).eq.'A').AND.(dofs.ne.'ALL')).OR.(dofs(1:1).eq.'B').OR.&
 &((dofs(1:1).eq.'C').AND.(dofs.ne.'CHI'))) then
          imol = molofrs(atmres(rs))
          if (rs.ne.atmol(imol,1)) then
            write(ilog,*) 'Warning. Requested rigid-body motion constraint for atom ',rs,&
 &', which does not correspond to the first atom in its molecule. Ignored.'
            cycle
          else if (izrot(rs)%alsz.le.0) then
            write(ilog,*) 'Fatal. Requested rigid-body motion constraint for atom ',rs,&
 &', but its moving atoms list is not set up. This is most likely a bug.'
            call fexit()
          end if
          if (dofs.eq.'XYZ') then
            do j=1,3
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else if (dofs(1:2).eq.'XY') then
            do j=1,2
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else if (dofs(1:2).eq.'XZ') then
            do jj=3,4
              j = jj
              if (jj.eq.4) j = 1
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else if (dofs(1:2).eq.'YZ') then
            do j=2,3
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else if (dofs(1:1).eq.'X') then
            do j=1,1
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else if (dofs(1:1).eq.'Y') then
            do j=2,2
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else if (dofs(1:1).eq.'Z') then
            do j=3,3
              if (dc_di(imol)%frz(j).EQV..false.) then
                dc_di(imol)%frz(j) = .true.
                n_constraints = n_constraints + 1
                ncons(27+j) = ncons(27+j) + 1
              end if
            end do
          else
            if (atmol(imol,2).eq.atmol(imol,1)) then
              write(ilog,*) 'Warning. Requested rigid-body rotational constraint for atom ',rs,&
 &', which is the only atom in a single-atom molecule. Ignored.'
              cycle
            else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
              write(ilog,*) 'Fatal. Two-atom molecules are not yet support&
 &ed in torsional/rigid-body dynamics. Check back later.'
              call fexit()
            end if
            if (dofs.eq.'ABC') then
              do j=4,6
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            else if (dofs(1:2).eq.'AB') then
              do j=4,5
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            else if (dofs(1:2).eq.'BC') then
              do j=5,6
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            else if (dofs(1:2).eq.'AC') then
              do jj=6,7
                j = jj
                if (jj.eq.7) j = 4
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            else if ((dofs(1:1).eq.'A').AND.(dofs.ne.'ALL')) then
              do j=4,4
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            else if (dofs(1:1).eq.'B') then
              do j=5,5
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            else if ((dofs(1:1).eq.'C').AND.(dofs.ne.'CHI')) then
              do j=6,6
                if (dc_di(imol)%frz(j).EQV..false.) then
                  dc_di(imol)%frz(j) = .true.
                  n_constraints = n_constraints + 1
                  ncons(27+j) = ncons(27+j) + 1
                end if
              end do
            end if
          end if
        else
          write(ilog,*) 'Unknown selection for internal coordinate space constraint request for &
   &Z-matrix line ',rs,' in ',frzfile(t1:t2),': ',dofs,' - ignored.'
        end if
      end if
!
    end do
!
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',frzfile(t1:t2),': ',abc(1:1),'! Fatal exit.'
    close(unit=iunit)
    call fexit()
  end if
  close(unit=iunit)
!
  deallocate(rslst)
!
! spend some time to pre-organize constraints to map constrained interactions to proper hierarchy level
! (only for dynamics): per molecule: 4) all frozen, 3) all torsions frozen, 2) all backbone and OTHER torsions frozen
  if (use_dyn.EQV..true.) then
    do imol=1,nmol
      frzidx = 4
      do j=1,min(size(dc_di(imol)%frz(:)),6)
        if ((dc_di(imol)%frz(j).EQV..false.).AND.(frzidx.gt.3)) then
          frzidx = 3
        end if
      end do
      do rs=rsmol(imol,1),rsmol(imol,2)
        if (wline(rs).gt.0) then
          if (wnr(rs).gt.0) then
            if ((dc_di(imol)%frz(wnr(rs)).EQV..false.).AND.(frzidx.gt.1)) then
              frzidx = 1
            end if
          end if
        end if
        if (fline(rs).gt.0) then
          if (fnr(rs).gt.0) then
            if ((dc_di(imol)%frz(fnr(rs)).EQV..false.).AND.(frzidx.gt.1)) then
              frzidx = 1
            end if
          end if
        end if
        if (yline(rs).gt.0) then
          if (ynr(rs).gt.0) then
            if ((dc_di(imol)%frz(ynr(rs)).EQV..false.).AND.(frzidx.gt.1)) then
              frzidx = 1
            end if
          end if
        end if
        do j=1,nnucs(rs)
          if (nucsnr(j,rs).gt.0) then
            if ((dc_di(imol)%frz(nucsnr(j,rs)).EQV..false.).AND.(frzidx.gt.1)) then
              frzidx = 1
            end if
          end if
        end do
        do j=1,nchi(rs)
          if (chinr(j,rs).gt.0) then
            if ((dc_di(imol)%frz(chinr(j,rs)).EQV..false.).AND.(frzidx.gt.2)) then
              frzidx = 2
            end if
          end if
        end do
      end do
      if (othidxmol(moltypid(imol)).gt.0) then
        do j=othidxmol(moltypid(imol)),ntormol(moltypid(imol))+6
          if ((dc_di(imol)%frz(j).EQV..false.).AND.(frzidx.gt.1)) then
            frzidx = 1 ! note that frzidx=2 is of little use anyway, so it's okay to be overeager here 
          end if
        end do
      end if
      molfrzidx(imol) = frzidx
    end do
  end if
!
! for Cartesian constraints: set total mass of atoms frozen completely
  sysmass(2) = 0.
  if (fycxyz.eq.2) then
    do j=1,n
      if ((cart_frz(j,1).EQV..true.).AND.(cart_frz(j,2).EQV..true.).AND.(cart_frz(j,3).EQV..true.)) then
        sysmass(2) = sysmass(2) + mass(j)
      end if
    end do
  end if
!
 45   format(4x,'Molecule #',i6)
 46   format(4x,'Residue  #',i6,' (',a3,')')
 52   format(i10,' degrees of freedom have been removed from ',a)
 47   format(i10,' residues have been removed from ',a)
 48   format(i10,' molecules have been removed from ',a)
 49   format(i10,' ',a,' angles have been frozen in internal coordinate space dynamics.')
 50   format(i10,' rigid body degrees of freedom have been frozen in internal coordinate space dynamics.')
 51   format(i10,' atomic ',a,'-coordinates have been explicitly frozen in Cartesian space dynamics.')
 53   format(i10,' molecular rigid body rotation axes have been frozen in internal coordinate space dynamics.')
 54   format(i10,' molecular rigid body translation axes have been frozen in internal coordinate space dynamics.')


!
  if (frz_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '----- Summary of constraints  -----'
    write(ilog,*) 
    if ((fycxyz.eq.1).OR.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
      write(ilog,*) '----- MC constraints  -----'
      write(ilog,*) '(lists only those that would be sampled by current move set)'
      write(ilog,*)
      if ((ncons(21).gt.0).AND.(have_rigid.EQV..true.).AND.(clurb_freq.lt.1.0)) then
        write(ilog,48) ncons(21),'single molecule rigid-body sampling.'
      end if
      if ((ncons(22).gt.0).AND.(have_clurb.EQV..true.)) then
        write(ilog,48) ncons(22),'multi-molecule rigid-body sampling.'
      end if
      if ((ncons(1).gt.0).AND.(have_pivot.EQV..true.)) then
        write(ilog,47) ncons(1),'polypeptide phi/psi pivot sampling.'
      end if
      if ((ncons(2).gt.0).AND.(have_omega.EQV..true.)) then
        write(ilog,47) ncons(2),'omega sampling (FMCSC_OMEGAFREQ).'
      end if
      if ((ncons(3).gt.0).AND.(have_chi.EQV..true.)) then
        write(ilog,47) ncons(3),'sidechain sampling (FMCSC_CHIFREQ).'
      end if
      if ((ncons(4).gt.0).AND.(have_nuc.EQV..true.).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0)) then
        write(ilog,47) ncons(4),'nucleic avid pivot sampling (FMCSC_NUCFREQ).'
      end if
      if ((ncons(5).gt.0).AND.(have_pucker.EQV..true.)) then
        write(ilog,47) ncons(5),'polypeptide pucker sampling (FMCSC_PKRFREQ).'
      end if
      if ((ncons(6).gt.0).AND.(have_nucpuck.EQV..true.)) then
        write(ilog,47) ncons(6),'sugar pucker sampling (FMCSC_SUGARFREQ).'
      end if
      if ((ncons(7).gt.0).AND.(have_natother.EQV..true.)) then
        write(ilog,52) ncons(7),'single dihedral angle pivot sampling of native CAMPARI d.o.f.s (FMCSC_OTHERNATFREQ).'
      end if
      if ((ncons(9).gt.0).AND.(have_unkother.EQV..true.)) then
        write(ilog,52) ncons(9),'single dihedral angle pivot sampling of unsupported residues (FMCSC_OTHERUNKFREQ).'
      end if
      if ((ncons(8).gt.0).AND.(have_unsother.EQV..true.)) then
        write(ilog,52) ncons(8),'single dihedral angle pivot sampling of nonnative CAMPARI d.o.f.s in supported residues &
 &(FMCSC_OTHERFREQ).'
      end if
      if ((ncons(10).gt.0).AND.(have_djcr.EQV..true.)) then
        write(ilog,47) ncons(10),'Dinner-Ulmschneider CR moves without omega sampling (FMCSC_TORCRFREQ).'
      end if
      if ((ncons(11).gt.0).AND.(have_docr.EQV..true.)) then
        write(ilog,47) ncons(11),'Dinner-Ulmschneider CR moves with omega sampling (FMCSC_TORCROFREQ).'
      end if
      if ((ncons(12).gt.0).AND.(have_nuccr.EQV..true.)) then
        write(ilog,47) ncons(12),'Dinner-Ulmschneider CR moves for nucleic acids (FMCSC_NUCCRFREQ).'
      end if
      if ((ncons(13).gt.0).AND.(have_ujcr.EQV..true.)) then
        write(ilog,47) ncons(13),'Ulmschneider-Jorgensen CR moves with bond angle sampling (FMCSC_ANGCRFREQ).'
      end if
      if ((ncons(14).gt.0).AND.(have_sjcr.EQV..true.)) then
        write(ilog,47) ncons(14),'Sjunnesson-Favrin approximate CR sampling (FMCSC_CRFREQ).'
      end if
      write(ilog,*)
    end if
    if ((fycxyz.eq.1).AND.(use_dyn.EQV..true.)) then
      write(ilog,*) '----- TMD constraints  -----'
      write(ilog,*)
      if (ncons(23).gt.0) then
        write(ilog,50) ncons(23)
      end if
      if (ncons(16).gt.0) then
        write(ilog,49) ncons(16),'polypeptide phi'
      end if
      if (ncons(17).gt.0) then
        write(ilog,49) ncons(17),'polypeptide psi'
      end if
      if (ncons(15).gt.0) then
        write(ilog,49) ncons(15),'peptide (amide) omega'
      end if
      if (ncons(19).gt.0) then
        write(ilog,49) ncons(19),'sidechain'
      end if
      if (ncons(18).gt.0) then
        write(ilog,49) ncons(18),'nucleic acid backbone'
      end if
      if (ncons(20).gt.0) then
        write(ilog,49) ncons(20),'nonnative CAMPARI (OTHER)'
      end if
      if (ncons(27).gt.0) then ! mutually exclusive
        write(ilog,49) ncons(27),'individual degrees of freedom (mode "A"), i.e., dihedral'
      end if
      if (sum(ncons(28:30)).gt.0) then
        write(ilog,54) sum(ncons(28:30))
      end if
      if (sum(ncons(31:33)).gt.0) then
        write(ilog,53) sum(ncons(31:33))
      end if
      write(ilog,*)
    end if
    if (fycxyz.eq.2) then
      write(ilog,*) '----- CMD constraints  -----'
      write(ilog,*) 
      if (ncons(24).gt.0) then
        write(ilog,51) ncons(24),'x'
      end if
      if (ncons(25).gt.0) then
        write(ilog,51) ncons(25),'y'
      end if
      if (ncons(26).gt.0) then
        write(ilog,51) ncons(26),'z'
      end if
      write(ilog,*)
    end if
    write(ilog,*)
    write(ilog,*) '------------------------------------'
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine remove_intdof(mode,idx,fyfl,omfl,chifl,pucfl,othfl,rslst,ncons)
!
  use iounit
  use molecule
  use sequen
  use movesets
  use forces
  use fyoc
  use system
  use interfaces
  use ujglobals
  use zmatrix
  use atoms
  use polypep
!
  implicit none
!
  integer mode,idx,k,imol,rs,j,jj,i,ii,nlst,nlst2,rslst(nseq),ncons(50)
  logical fyfl,omfl,chifl,pucfl,othfl,eliminate
!
!
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
!
! some minor consistency checks (may actually be safer to remove these eventually)
    if ((dyn_mode.eq.1).AND.((rigidfreq.ge.1.0).OR.(particleflucfreq.ge.1.0))) then
      write(ilog,*) 'Warning. Constraining internal degrees of freed&
 &om makes no sense if those are not actually sampled. Ignored ...'
      return
    end if
!
    if (((have_other.EQV..false.).OR.((other_unkfreq.le.0.0).AND.(other_natfreq.ge.1.0))).AND.(othfl.EQV..true.).AND.&
 &(fyfl.EQV..false.).AND.(omfl.EQV..false.).AND.(chifl.EQV..false.).AND.(pucfl.EQV..false.).AND.(dyn_mode.eq.1)) then
      write(ilog,*) 'Warning. Constraining unsupported degrees of free&
 &dom makes no sense if those are not actually sampled (via OTHER). Ignored ...'
      return
    end if
!
!   assemble a list of residues to operate on
    if (mode.eq.3) then
      nlst = 1
      rslst(nlst) = idx
    else if (mode.eq.1) then
      imol = idx 
      nlst = 0
      do rs=rsmol(imol,1),rsmol(imol,2)
        nlst = nlst + 1
        rslst(nlst) = rs
      end do
    else if (mode.eq.2) then
      nlst = 0
      nlst2 = 0
      do imol=moltyp(idx,1),nmol
        if (moltypid(imol).eq.idx) then
          nlst2 = nlst2 + 1
          do rs=rsmol(imol,1),rsmol(imol,2)
            nlst = nlst + 1
            rslst(nlst) = rs
          end do
        end if
        if (nlst2.eq.moltyp(idx,2)) exit
      end do
    else
      write(ilog,*) 'Called remove_intdof(...) with unknown mode (of&
 &fending mode is ',mode,'). This is fatal.'
      call fexit()
    end if
!
!   now scan each residue in the list 
    do nlst2=1,nlst
      rs = rslst(nlst2)
      if ((fylst%nr.gt.0).AND.(fyfl.EQV..true.)) then
        call binary_search(fylst%nr,fylst%idx(1:fylst%nr),rs,k)
        if (fylst%idx(max(1,min(fylst%nr,k))).eq.rs) then
          do j=max(1,min(fylst%nr,k)),fylst%nr-1
            fylst%idx(j) = fylst%idx(j+1)
          end do
          fylst%nr = fylst%nr - 1
          ncons(1) = ncons(1) + 1
        end if
      end if
      if ((wlst%nr.gt.0).AND.(omfl.EQV..true.)) then
        call binary_search(wlst%nr,wlst%idx(1:wlst%nr),rs,k)
        if (wlst%idx(max(1,min(wlst%nr,k))).eq.rs) then
          do j=max(1,min(wlst%nr,k)),wlst%nr-1
            wlst%idx(j) = wlst%idx(j+1)
          end do
          wlst%nr = wlst%nr - 1
          ncons(2) = ncons(2) + 1
        end if
      end if
      if ((chilst%nr.gt.0).AND.(chifl.EQV..true.)) then
        call binary_search(chilst%nr,chilst%idx(1:chilst%nr),rs,k)
        if (chilst%idx(max(1,min(chilst%nr,k))).eq.rs) then
          do j=max(1,min(chilst%nr,k)),chilst%nr-1
            chilst%idx(j) = chilst%idx(j+1)
          end do
          chilst%nr = chilst%nr - 1
          ncons(3) = ncons(3) + 1
        end if
      end if
      if ((nuclst%nr.gt.0).AND.(fyfl.EQV..true.)) then
        call binary_search(nuclst%nr,nuclst%idx(1:nuclst%nr),rs,k)
        if (nuclst%idx(max(1,min(nuclst%nr,k))).eq.rs) then
          do j=max(1,min(nuclst%nr,k)),nuclst%nr-1
            nuclst%idx(j) = nuclst%idx(j+1)
          end do
          nuclst%nr = nuclst%nr - 1
          ncons(4) = ncons(4) + 1
        end if
      end if
      if ((puclst%nr.gt.0).AND.(pucfl.EQV..true.)) then
        call binary_search(puclst%nr,puclst%idx(1:puclst%nr),rs,k)
        if (puclst%idx(max(1,min(puclst%nr,k))).eq.rs) then
          do j=max(1,min(puclst%nr,k)),puclst%nr-1
            puclst%idx(j) = puclst%idx(j+1)
          end do
          puclst%nr = puclst%nr - 1
          ncons(5) = ncons(5) + 1
        end if
      end if
      if ((nucpuclst%nr.gt.0).AND.(pucfl.EQV..true.)) then
        call binary_search(nucpuclst%nr,nucpuclst%idx(1:nucpuclst%nr),rs,k)
        if (nucpuclst%idx(max(1,min(nucpuclst%nr,k))).eq.rs) then
          do j=max(1,min(nucpuclst%nr,k)),nucpuclst%nr-1
            nucpuclst%idx(j) = nucpuclst%idx(j+1)
          end do
          nucpuclst%nr = nucpuclst%nr - 1
          ncons(6) = ncons(6) + 1
        end if
      end if
      if ((natlst%nr.gt.0).AND.((fyfl.EQV..true.).OR.(omfl.EQV..true.).OR.(chifl.EQV..true.))) then
        call binary_search(natlst%nr,natlst%idx(1:natlst%nr),at(rs)%bb(1),k)
        i = max(1,min(natlst%nr,k))
        if (i.gt.1) then
          i = i-1
        end if
        jj = i
        do while (jj.le.natlst%nr)
          eliminate = .false.
          if ((omfl.EQV..true.).AND.(natlst%idx(jj).eq.wline(rs))) eliminate = .true.
          if (fyfl.EQV..true.) then
            if ((natlst%idx(jj).eq.fline(rs)).OR.(natlst%idx(jj).eq.fline2(rs)).OR.(natlst%idx(jj).eq.yline(rs)).OR.&
                (natlst%idx(jj).eq.yline2(rs))) eliminate = .true.
            do ii=1,nnucs(rs)
              if (natlst%idx(jj).eq.nucsline(ii,rs)) eliminate = .true.
            end do
          end if
          if (chifl.EQV..true.) then
            do ii=1,nchi(rs)
              if (natlst%idx(jj).eq.chiline(ii,rs)) eliminate = .true.
            end do
          end if
          if (eliminate.EQV..true.) then
            do j=jj,natlst%nr-1
              natlst%idx(j) = natlst%idx(j+1)
            end do
            natlst%nr = natlst%nr - 1
            ncons(7) = ncons(7) + 1
          else
            jj = jj + 1
          end if
          if (jj.gt.natlst%nr) exit
          if (atmres(natlst%idx(jj)).gt.(rs+1)) exit
        end do
      end if
      if ((unslst%nr.gt.0).AND.(othfl.EQV..true.)) then
        call binary_search(unslst%nr,unslst%idx(1:unslst%nr),at(rs)%bb(1),k)
        i = max(1,min(unslst%nr,k))
        if (i.gt.1) then
          i = i-1
        end if
        jj = i
        do while (jj.le.unslst%nr)
          if ((atmres(iz(1,unslst%idx(jj))).eq.rs)) then
            do j=jj,unslst%nr-1
              unslst%idx(j) = unslst%idx(j+1)
            end do
            unslst%nr = unslst%nr - 1
            ncons(8) = ncons(8) + 1
          else
            jj = jj + 1
          end if
          if (jj.gt.unslst%nr) exit
          if (atmres(unslst%idx(jj)).gt.(rs+1)) exit
        end do
      end if
      if ((unklst%nr.gt.0).AND.(othfl.EQV..true.)) then
        call binary_search(unklst%nr,unklst%idx(1:unklst%nr),at(rs)%bb(1),k)
        i = max(1,min(unklst%nr,k))
        if (i.gt.1) then
          i = i-1
        end if
        jj = i
        do while (jj.le.unklst%nr)
          if ((atmres(iz(1,unklst%idx(jj))).eq.rs)) then
            do j=jj,unklst%nr-1
              unklst%idx(j) = unklst%idx(j+1)
            end do
            unklst%nr = unklst%nr - 1
            ncons(9) = ncons(9) + 1
          else
            jj = jj + 1
          end if
          if (jj.gt.unklst%nr) exit
          if (atmres(unklst%idx(jj)).gt.(rs+1)) exit
        end do
      end if
      if ((djcrlst%nr.gt.0).AND.((fyfl.EQV..true.).OR.(omfl.EQV..true.).OR.(pucfl.EQV..true.))) then
        j = 1
        do while (j.le.djcrlst%nr)
          if ((molofrs(rs).eq.molofrs(djcrlst%idx(j))).AND.(djcrlst%idx(j).ge.rs).AND.(djcrlst%idx2(j).le.rs)) then
            if ((djcrlst%idx(j)-rs).ge.4) then
              djcrlst%idx2(j) = rs + 1
              j = j + 1
            else
              djcrlst%nr = djcrlst%nr - 1
              ncons(10) = ncons(10) + 1
              do k=j,djcrlst%nr
                djcrlst%idx(k) = djcrlst%idx(k+1)
                djcrlst%idx2(k) = djcrlst%idx2(k+1)
              end do
            end if
          else
            j = j + 1
          end if
        end do
      end if
      if ((docrlst%nr.gt.0).AND.((fyfl.EQV..true.).OR.(omfl.EQV..true.).OR.(pucfl.EQV..true.))) then
        j = 1
        do while (j.le.docrlst%nr)
          if ((molofrs(rs).eq.molofrs(docrlst%idx(j))).AND.(docrlst%idx(j).ge.rs).AND.(docrlst%idx2(j).le.rs)) then
            if ((docrlst%idx(j)-rs).ge.3) then
              docrlst%idx2(j) = rs + 1
              j = j + 1
            else
              docrlst%nr = docrlst%nr - 1
              ncons(11) = ncons(11) + 1
              do k=j,docrlst%nr
                docrlst%idx(k) = docrlst%idx(k+1)
                docrlst%idx2(k) = docrlst%idx2(k+1)
              end do
            end if
          else
            j = j + 1
          end if
        end do
      end if
      if ((nuccrlst%nr.gt.0).AND.(fyfl.EQV..true.)) then
        j = 1
        do while (j.le.nuccrlst%nr)
          if ((molofrs(rs).eq.molofrs(nuccrlst%idx(j))).AND.(nuccrlst%idx(j).ge.rs).AND.(nuccrlst%idx2(j).le.rs)) then
            if ((nuccrlst%idx(j)-rs).ge.2) then
              nuccrlst%idx2(j) = rs + 1
              j = j + 1
            else
              nuccrlst%nr = nuccrlst%nr - 1
              ncons(12) = ncons(12) + 1
              do k=j,nuccrlst%nr
                nuccrlst%idx(k) = nuccrlst%idx(k+1)
                nuccrlst%idx2(k) = nuccrlst%idx2(k+1)
              end do
            end if
          else
            j = j + 1
          end if
        end do
      end if
      if ((ujcrlst%nr.gt.0).AND.(fyfl.EQV..true.)) then
        j = 1
        do while (j.le.ujcrlst%nr)
          if ((molofrs(rs).eq.molofrs(ujcrlst%idx(j))).AND.(ujcrlst%idx(j).ge.rs).AND.(ujcrlst%idx2(j).le.rs)) then
            if ((ujcrlst%idx(j)-rs).ge.ujminsz) then
              ujcrlst%idx2(j) = rs + 1
              j = j + 1
            else
              ujcrlst%nr = ujcrlst%nr - 1
              ncons(13) = ncons(13) + 1
              do k=j,ujcrlst%nr
                ujcrlst%idx(k) = ujcrlst%idx(k+1)
                ujcrlst%idx2(k) = ujcrlst%idx2(k+1)
              end do
            end if
          else
            j = j + 1
          end if
        end do
      end if
      if ((sjcrlst%nr.gt.0).AND.(fyfl.EQV..true.)) then
        j = 1
        do while (j.le.sjcrlst%nr)
          if ((molofrs(rs).eq.molofrs(sjcrlst%idx(j))).AND.(sjcrlst%idx(j).ge.rs).AND.(sjcrlst%idx2(j).le.rs)) then
            if ((sjcrlst%idx(j)-rs).ge.(nr_crres)) then
              sjcrlst%idx2(j) = rs + 1
              j = j + 1
            else
              sjcrlst%nr = sjcrlst%nr - 1
              ncons(14) = ncons(14) + 1
              do k=j,sjcrlst%nr
                sjcrlst%idx(k) = sjcrlst%idx(k+1)
                sjcrlst%idx2(k) = sjcrlst%idx2(k+1)
              end do
            end if
          else
            j = j + 1
          end if
        end do
      end if
    end do
!
!   let us simplify by not allowing any originally active type of move to be removed entirely by constraints
    if (pdb_analyze.EQV..false.) then
      if ((have_pivot.EQV..true.).AND.(fylst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for polypeptide backbone pivot &
   &phi/psi moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_omega.EQV..true.).AND.(wlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for polypeptide backbone omega &
   &moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_chi.EQV..true.).AND.(chilst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for sidechain &
   &moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_nuc.EQV..true.).AND.(nuclst%nr.eq.0).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for NUC-type pivot &
   &moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_pucker.EQV..true.).AND.(puclst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for polypeptide pucker &
   &moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_nucpuck.EQV..true.).AND.(nucpuclst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for nucleic acid pucker &
   &moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_djcr.EQV..true.).AND.(djcrlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for Dinner-Ulmschneider backbone pivot &
   &concerted rotation moves without omega sampling entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_docr.EQV..true.).AND.(docrlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for Dinner-Ulmschneider polypeptide &
   &concerted rotation moves with omega sampling entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_nuccr.EQV..true.).AND.(nuccrlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for Dinner-Ulmschneider nucleic acid &
   &concerted rotation moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_ujcr.EQV..true.).AND.(ujcrlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for Ulmschneider-Jorgensen polypeptide &
   &concerted rotation moves with bond angle changes entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_sjcr.EQV..true.).AND.(sjcrlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for Sjunnesson-Favrin polypeptide &
   &concerted rotation moves entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_unkother.EQV..true.).AND.(have_other.EQV..true.).AND.(unklst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for single torsion pivot moves &
   &in unsupported residues entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_natother.EQV..true.).AND.(have_other.EQV..true.).AND.(natlst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for single torsion pivot moves &
   &(OTHER) on natively supported dihedral angles entirely. Please adjust the move set instead.'
        call fexit()
      end if
      if ((have_unsother.EQV..true.).AND.(have_other.EQV..true.).AND.(unslst%nr.eq.0)) then
        write(ilog,*) 'Fatal. The constraints requested deplete the eligible residues for single torsion pivot moves &
   &(OTHER) on unsupported dihedral angles in supported residues entirely. Please adjust the move set instead.'
        call fexit()
      end if
    end if
  end if ! if using MC moves
!
  if ((use_dyn.EQV..true.).AND.(fycxyz.eq.1)) then
!   (re)assemble a list of residues to operate on
    if (mode.eq.3) then
      nlst = 1
      rslst(nlst) = idx
    else if (mode.eq.1) then
      imol = idx 
      nlst = 0
      do rs=rsmol(imol,1),rsmol(imol,2)
        nlst = nlst + 1
        rslst(nlst) = rs
      end do
    else if (mode.eq.2) then
      nlst = 0
      nlst2 = 0
      do imol=moltyp(idx,1),nmol
        if (moltypid(imol).eq.idx) then
          nlst2 = nlst2 + 1
          do rs=rsmol(imol,1),rsmol(imol,2)
            nlst = nlst + 1
            rslst(nlst) = rs
          end do
        end if
        if (nlst2.eq.moltyp(idx,2)) exit
      end do
    else
      write(ilog,*) 'Called remove_intdof(...) with unknown mode (of&
 &fending mode is ',mode,'). This is fatal.'
      call fexit()
    end if
!
!   now scan each residue in the list 
    do nlst2=1,nlst
      rs = rslst(nlst2)
      imol = molofrs(rs)
      if (wnr(rs).gt.0) then
        if (omfl.EQV..true.) then
          if (dc_di(imol)%frz(wnr(rs)).EQV..false.) then
            dc_di(imol)%frz(wnr(rs)) = .true.
            n_constraints = n_constraints + 1
            ncons(15) = ncons(15) + 1
          end if
        end if
      end if
      if (fnr(rs).gt.0) then
        if (fyfl.EQV..true.) then
          if (dc_di(imol)%frz(fnr(rs)).EQV..false.) then
            dc_di(imol)%frz(fnr(rs)) = .true.
            n_constraints = n_constraints + 1
            ncons(16) = ncons(16) + 1
          end if
        end if
      end if
      if (ynr(rs).gt.0) then
        if (fyfl.EQV..true.) then
          if (dc_di(imol)%frz(ynr(rs)).EQV..false.) then
            dc_di(imol)%frz(ynr(rs)) = .true.
            n_constraints = n_constraints + 1
            ncons(17) = ncons(17) + 1
          end if
        end if
      end if
      do j=1,nnucs(rs)
        if (nucsnr(j,rs).le.0) cycle
        if (fyfl.EQV..true.) then
          if (dc_di(imol)%frz(nucsnr(j,rs)).EQV..false.) then
            dc_di(imol)%frz(nucsnr(j,rs)) = .true.
            n_constraints = n_constraints + 1
            ncons(18) = ncons(18) + 1
          end if
        end if
      end do
      do j=1,nchi(rs)
        if (chinr(j,rs).le.0) cycle
        if (chifl.EQV..true.) then
          if (dc_di(imol)%frz(chinr(j,rs)).EQV..false.) then
            dc_di(imol)%frz(chinr(j,rs)) = .true.
            n_constraints = n_constraints + 1
            ncons(19) = ncons(19) + 1
          end if
        end if
      end do
      if (othfl.EQV..true.) then
        do j=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1
          if (izrot(j)%alsz.le.0) cycle
          if (j.eq.at(rsmol(imol,1))%bb(1)) cycle
          if (dc_di(imol)%recurs(izrot(j)%treevs(4),3).ge.othidxmol(moltypid(imol))) then
            if (dc_di(imol)%frz(dc_di(imol)%recurs(izrot(j)%treevs(4),3)).EQV..false.) then
              dc_di(imol)%frz(dc_di(imol)%recurs(izrot(j)%treevs(4),3)) = .true.
              n_constraints = n_constraints + 1
              ncons(20) = ncons(20) + 1
            end if
          end if
        end do
     end if
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine remove_rbcdof(mode,idx,ncons)
!
  use iounit
  use molecule
  use sequen
  use movesets
  use system
  use forces
  use interfaces
!
  implicit none
!
  integer mode,idx,k,imol,j,ttc,jmol,ncons(50)
  integer, ALLOCATABLE:: mollst(:)
!
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
!
    if (((particleflucfreq.ge.1.0).OR.(rigidfreq.le.0.0)).AND.(dyn_mode.eq.1)) then
      write(ilog,*) 'Warning. Constraining rigid-body degrees of fre&
 &edom makes no sense if those are not actually sampled. Ignored.'
      return
    end if
!
!   molecule-dep. is trivial
    if (mode.eq.1) then
      ttc = 1
      allocate(mollst(ttc))
      mollst(1) = idx
!   type-dep. too
    else if (mode.eq.2) then
      ttc = 0
      allocate(mollst(moltyp(idx,2)))
      do imol=1,nmol
        if (moltypid(imol).eq.idx) then
          ttc = ttc + 1
          mollst(ttc) = imol
        end if
      end do
    else
      write(ilog,*) 'Called remove_rbcdof(...) with unknown mode (of&
 &fending mode is ',mode,'). This is fatal.'
      call fexit()
    end if
!
    do jmol=1,ttc
      imol = mollst(jmol)
      if (rblst%nr.gt.0) then
        call binary_search(rblst%nr,rblst%idx(1:rblst%nr),imol,k)
        if (rblst%idx(max(1,min(rblst%nr,k))).eq.imol) then
          ncons(21) = ncons(21) + 1
          do j=max(1,min(rblst%nr,k)),rblst%nr-1
            rblst%idx(j) = rblst%idx(j+1)
          end do
          rblst%nr = rblst%nr - 1
        end if
      end if
      if (clurblst%nr.gt.0) then
        call binary_search(clurblst%nr,clurblst%idx(1:clurblst%nr),imol,k)
        if (clurblst%idx(max(1,min(clurblst%nr,k))).eq.imol) then
          ncons(22) = ncons(22) + 1
          do j=max(1,min(clurblst%nr,k)),clurblst%nr-1
            clurblst%idx(j) = clurblst%idx(j+1)
          end do
          clurblst%nr = clurblst%nr - 1
        end if
      end if
    end do

!
    deallocate(mollst)
!
    if ((have_rigid.EQV..true.).AND.(clurb_freq.lt.1.0).AND.(rblst%nr.eq.0)) then
      write(ilog,*) 'Fatal. The constraints requested deplete the eligible molecules for single molecule, rigid-body &
 &moves entirely. Please adjust the move set instead.'
      call fexit()
    end if
    if ((have_clurb.EQV..true.).AND.(clurblst%nr.le.1)) then
      write(ilog,*) 'Fatal. The constraints requested deplete the eligible molecules for cluster rigid-body &
 &moves. Please adjust the move set instead.'
      call fexit()
    end if

  end if
!
! post-consistency check
  if ((have_clurb.EQV..true.).AND.(clurblst%nr.lt.clurb_maxsz)) then
    write(ilog,*)' Warning. After application of constraints, maximum cluster size for cluster rigid &
 &body moves is reduced to the remaining number of mobile molecules (',clurblst%nr,').'
    clurb_maxsz = clurblst%nr
  end if
!
  if ((use_dyn.EQV..true.).AND.(fycxyz.eq.1)) then
!   molecule-based
    if (mode.eq.1) then
      imol = idx
      if (atmol(imol,2).eq.atmol(imol,1)) then
        ttc = 3
      else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        ttc = 5
      else
        ttc = 6
      end if
      do j=1,ttc
        if (dc_di(imol)%frz(j).EQV..false.) then
          dc_di(imol)%frz(j) = .true.
          n_constraints = n_constraints + 1
          ncons(23) = ncons(23) + 1
        end if
      end do
    else if (mode.eq.2) then
      do imol=1,nmol
        if (moltypid(imol).ne.idx) cycle
        if (atmol(imol,2).eq.atmol(imol,1)) then
          ttc = 3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ttc = 5
        else
          ttc = 6
        end if
        do j=1,ttc
          if (dc_di(imol)%frz(j).EQV..false.) then
            dc_di(imol)%frz(j) = .true.
            n_constraints = n_constraints + 1
            ncons(23) = ncons(23) + 1
          end if
        end do
      end do
    else
      write(ilog,*) 'Called remove_rbcdof(...) with unknown mode (of&
 &fending mode is ',mode,'). This is fatal.'
      call fexit()
    end if
  end if
!
end
!
!----------------------------------------------------------------------------
!
! this function reads the (optional) input file for preferential sampling
! (mostly residue-based), and sets the desired sampling weights
! it is called after all other move set operations (frzfile, etc.) have been
! handled!
!
subroutine preferential_sampling_setup()
!
  use iounit
  use movesets
  use interfaces
  use sequen
  use system
  use molecule
  use zmatrix
  use atoms
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,st1,st2,ii,i,k,kk
  character(MAXSTRLEN) str2,str3
  logical exists,inrslp,inmllp,inmtlp,touched(18),inatlp,inother
  RTYPE ivwt,wt_vec(11)
!
! first initialize weight lists for default case
  if (fylst%nr.gt.0) then
    ivwt = 1.0/dble(fylst%nr)
    do i=1,fylst%nr-1
      fylst%wt(i) = dble(i)*ivwt
    end do
    fylst%wt(fylst%nr) = 1.0
  end if
  if (wlst%nr.gt.0) then
    ivwt = 1.0/dble(wlst%nr)
    do i=1,wlst%nr-1
      wlst%wt(i) = dble(i)*ivwt
    end do
    wlst%wt(wlst%nr) = 1.0
  end if
  if (chilst%nr.gt.0) then
    ivwt = 1.0/dble(chilst%nr)
    do i=1,chilst%nr-1
      chilst%wt(i) = dble(i)*ivwt
    end do
    chilst%wt(chilst%nr) = 1.0
  end if
  if (nuclst%nr.gt.0) then
    ivwt = 1.0/dble(nuclst%nr)
    do i=1,nuclst%nr-1
      nuclst%wt(i) = dble(i)*ivwt
    end do
    nuclst%wt(nuclst%nr) = 1.0
  end if
  if (puclst%nr.gt.0) then
    ivwt = 1.0/dble(puclst%nr)
    do i=1,puclst%nr-1
      puclst%wt(i) = dble(i)*ivwt
    end do
    puclst%wt(puclst%nr) = 1.0
  end if
  if (nucpuclst%nr.gt.0) then
    ivwt = 1.0/dble(nucpuclst%nr)
    do i=1,nucpuclst%nr-1
      nucpuclst%wt(i) = dble(i)*ivwt
    end do
    nucpuclst%wt(nucpuclst%nr) = 1.0
  end if
  if (djcrlst%nr.gt.0) then
    allocate(djcrlst%wt(djcrlst%nr))
    ivwt = 1.0/dble(djcrlst%nr)
    do i=1,djcrlst%nr-1
      djcrlst%wt(i) = dble(i)*ivwt
    end do
    djcrlst%wt(djcrlst%nr) = 1.0
  end if
  if (docrlst%nr.gt.0) then
    allocate(docrlst%wt(docrlst%nr))
    ivwt = 1.0/dble(docrlst%nr)
    do i=1,docrlst%nr-1
      docrlst%wt(i) = dble(i)*ivwt
    end do
    docrlst%wt(docrlst%nr) = 1.0
  end if
  if (nuccrlst%nr.gt.0) then
    allocate(nuccrlst%wt(nuccrlst%nr))
    ivwt = 1.0/dble(nuccrlst%nr)
    do i=1,nuccrlst%nr-1
      nuccrlst%wt(i) = dble(i)*ivwt
    end do
    nuccrlst%wt(nuccrlst%nr) = 1.0
  end if
  if (sjcrlst%nr.gt.0) then
    allocate(sjcrlst%wt(sjcrlst%nr))
    ivwt = 1.0/dble(sjcrlst%nr)
    do i=1,sjcrlst%nr-1
      sjcrlst%wt(i) = dble(i)*ivwt
    end do
    sjcrlst%wt(sjcrlst%nr) = 1.0
  end if
  if (ujcrlst%nr.gt.0) then
    allocate(ujcrlst%wt(ujcrlst%nr))
    ivwt = 1.0/dble(ujcrlst%nr)
    do i=1,ujcrlst%nr-1
      ujcrlst%wt(i) = dble(i)*ivwt
    end do
    ujcrlst%wt(ujcrlst%nr) = 1.0
  end if
  if (natlst%nr.gt.0) then
    ivwt = 1.0/dble(natlst%nr)
    do i=1,natlst%nr-1
      natlst%wt(i) = dble(i)*ivwt
    end do
    natlst%wt(natlst%nr) = 1.0
  end if
  if (unklst%nr.gt.0) then
    ivwt = 1.0/dble(unklst%nr)
    do i=1,unklst%nr-1
      unklst%wt(i) = dble(i)*ivwt
    end do
    unklst%wt(unklst%nr) = 1.0
  end if
  if (unslst%nr.gt.0) then
    ivwt = 1.0/dble(unslst%nr)
    do i=1,unslst%nr-1
      unslst%wt(i) = dble(i)*ivwt
    end do
    unslst%wt(unslst%nr) = 1.0
  end if
  if (rblst%nr.gt.0) then
    allocate(rblst%wt(rblst%nr))
    ivwt = 1.0/dble(rblst%nr)
    do i=1,rblst%nr-1
      rblst%wt(i) = dble(i)*ivwt
    end do
    rblst%wt(rblst%nr) = 1.0
  end if
  if (clurblst%nr.gt.0) then
    allocate(clurblst%wt(clurblst%nr))
    ivwt = 1.0/dble(clurblst%nr)
    do i=1,clurblst%nr-1
      clurblst%wt(i) = dble(i)*ivwt
    end do
    clurblst%wt(clurblst%nr) = 1.0
  end if
  if (fluclst%nr.gt.0) then
    allocate(fluclst%wt(fluclst%nr))
    ivwt = 1.0/dble(fluclst%nr)
    do i=1,fluclst%nr-1
      fluclst%wt(i) = dble(i)*ivwt
    end do
    fluclst%wt(fluclst%nr) = 1.0
  end if
!
  call strlims(prefsamplingfile,t1,t2)
  inquire(file=prefsamplingfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Warning. Cannot open assumed input file (',&
 &prefsamplingfile(t1:t2),') for FMCSC_PSWFILE. Assuming default weights (see documentation).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=prefsamplingfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for preferential sampling input&
 & (',prefsamplingfile(t1:t2),'). Assuming default weights.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing preferential sampling&
 & input (got: ',prefsamplingfile(t1:t2),').'
    call fexit()
  end if
!
! if we're here, exit via "return" is no longer allowed -> reset weights to unity
  if (fylst%nr.gt.0) fylst%wt(:) = 1.0
  if (wlst%nr.gt.0) wlst%wt(:) = 1.0
  if (chilst%nr.gt.0) chilst%wt(:) = 1.0
  if (nuclst%nr.gt.0) nuclst%wt(:) = 1.0
  if (puclst%nr.gt.0) puclst%wt(:) = 1.0
  if (nucpuclst%nr.gt.0) nucpuclst%wt(:) = 1.0
  if (djcrlst%nr.gt.0) djcrlst%wt(:) = 1.0
  if (docrlst%nr.gt.0) docrlst%wt(:) = 1.0
  if (nuccrlst%nr.gt.0) nuccrlst%wt(:) = 1.0
  if (ujcrlst%nr.gt.0) ujcrlst%wt(:) = 1.0
  if (sjcrlst%nr.gt.0) sjcrlst%wt(:) = 1.0
  if (rblst%nr.gt.0) rblst%wt(:) = 1.0
  if (clurblst%nr.gt.0) clurblst%wt(:) = 1.0
  if (fluclst%nr.gt.0) fluclst%wt(:) = 1.0
  if (unklst%nr.gt.0) unklst%wt(:) = 1.0
  if (unslst%nr.gt.0) unslst%wt(:) = 1.0
  if (natlst%nr.gt.0) natlst%wt(:) = 1.0
!
  inrslp = .false.
  inmllp = .false.
  inmtlp = .false.
  inatlp = .false.
  touched(:) = .false.
!
  call strlims_quiet(str2,st1,st2)
  do while (.true.)
    if ((str2(st1:st1).eq.'r').OR.(str2(st1:st1).eq.'R')) then
      if (inmllp.EQV..true.) then
        inmllp = .false.
      end if
      if (inmtlp.EQV..true.) then
        inmtlp = .false.
      end if
      if (inatlp.EQV..true.) then
        inatlp = .false.
      end if
      inrslp = .true.
    else if ((str2(st1:st1).eq.'M').OR.(str2(st1:st1).eq.'m')) then
      if (inrslp.EQV..true.) then
        inrslp = .false.
      end if
      if (inmtlp.EQV..true.) then
        inmtlp = .false.
      end if
      if (inatlp.EQV..true.) then
        inatlp = .false.
      end if
      inmllp = .true.
    else if ((str2(st1:st1).eq.'T').OR.(str2(st1:st1).eq.'t')) then
      if (inrslp.EQV..true.) then
        inrslp = .false.
      end if
      if (inmllp.EQV..true.) then
        inmllp = .false.
      end if
      if (inatlp.EQV..true.) then
        inatlp = .false.
      end if
      inmtlp = .true.
    else if ((str2(st1:st1).eq.'A').OR.(str2(st1:st1).eq.'a')) then
      if (inrslp.EQV..true.) then
        inrslp = .false.
      end if
      if (inmllp.EQV..true.) then
        inmllp = .false.
      end if
      if (inmtlp.EQV..true.) then
        inmtlp = .false.
      end if
      inatlp = .true.
    else if (str2(st1:st1).eq.' ') then
!     do nothing for emtpy line
    else
      ii = 1
      k = -1
      call extract_int(str2,k,ii)
      if ((inrslp.EQV..true.).AND.((k.le.nseq).AND.(k.ge.1))) then
        str3 = str2(ii:MAXSTRLEN)
        wt_vec(:) = 1.0
        read(str3,*,iostat=iomess) (wt_vec(i), i=1,11) ! truncation is tolerated 
        if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing preferential sampling&
 & input (got: ',prefsamplingfile(t1:t2),').'
          call fexit()
        end if
!       now go through each of the 11 residue-based lists, locate correct entry, and set weight to target
        if ((wt_vec(1).ge.0.0).AND.(fylst%nr.gt.0)) then
          call binary_search(fylst%nr,fylst%idx(1:fylst%nr),k,kk)
          kk = min(max(1,kk),fylst%nr)
          if (k.ne.fylst%idx(kk)) then
            if (wt_vec(1).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for phi/psi sampling for &
 &residue ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(1).ne.fylst%wt(kk)) then
            touched(1) = .true.
            fylst%wt(kk) = wt_vec(1)
          end if
        end if
        if ((wt_vec(2).ge.0.0).AND.(wlst%nr.gt.0)) then
          call binary_search(wlst%nr,wlst%idx(1:wlst%nr),k,kk)
          kk = min(max(1,kk),wlst%nr)
          if (k.ne.wlst%idx(kk)) then
            if (wt_vec(2).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for omega sampling for &
 &residue ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(2).ne.wlst%wt(kk)) then
            touched(2) = .true.
            wlst%wt(kk) = wt_vec(2)
          end if
        end if
        if ((wt_vec(3).ge.0.0).AND.(chilst%nr.gt.0)) then
          call binary_search(chilst%nr,chilst%idx(1:chilst%nr),k,kk)
          kk = min(max(1,kk),chilst%nr)
          if (k.ne.chilst%idx(kk)) then
            if (wt_vec(3).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for chi sampling for &
 &residue ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(3).ne.chilst%wt(kk)) then
            touched(3) = .true.
            chilst%wt(kk) = wt_vec(3)
          end if
        end if
        if ((wt_vec(4).ge.0.0).AND.(nuclst%nr.gt.0)) then
          call binary_search(nuclst%nr,nuclst%idx(1:nuclst%nr),k,kk)
          kk = min(max(1,kk),nuclst%nr)
          if (k.ne.nuclst%idx(kk)) then
            if (wt_vec(4).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for pivot NUC sampling for &
 &residue ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(4).ne.nuclst%wt(kk)) then
            touched(4) = .true.
            nuclst%wt(kk) = wt_vec(4)
          end if
        end if
        if ((wt_vec(5).ge.0.0).AND.(puclst%nr.gt.0)) then
          call binary_search(puclst%nr,puclst%idx(1:puclst%nr),k,kk)
          kk = min(max(1,kk),puclst%nr)
          if (k.ne.puclst%idx(kk)) then
            if (wt_vec(5).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for polypeptide pucker sampling for &
 &residue ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(5).ne.puclst%wt(kk)) then
            touched(5) = .true.
            puclst%wt(kk) = wt_vec(5)
          end if
        end if
        if ((wt_vec(6).ge.0.0).AND.(nucpuclst%nr.gt.0)) then
          call binary_search(nucpuclst%nr,nucpuclst%idx(1:nucpuclst%nr),k,kk)
          kk = min(max(1,kk),nucpuclst%nr)
          if (k.ne.nucpuclst%idx(kk)) then
            if (wt_vec(6).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for nucleic acid pucker sampling for &
 &residue ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(6).ne.nucpuclst%wt(kk)) then
            touched(6) = .true.
            nucpuclst%wt(kk) = wt_vec(6)
          end if
        end if
        if ((wt_vec(7).ge.0.0).AND.(djcrlst%nr.gt.0)) then
          call binary_search(djcrlst%nr,djcrlst%idx(1:djcrlst%nr),k,kk)
          kk = min(max(1,kk),djcrlst%nr)
          if (k.ne.djcrlst%idx(kk)) then
            if (wt_vec(7).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for Dinner-Ulmschneider polypeptide concerted rotations without &
 &omega sampling for residue ',k,', which is not eligible as a final residue in a stretch. Ignored.'
            end if
          else if (wt_vec(7).ne.djcrlst%wt(kk)) then
            touched(7) = .true.
            djcrlst%wt(kk) = wt_vec(7)
          end if
        end if
        if ((wt_vec(8).ge.0.0).AND.(docrlst%nr.gt.0)) then
          call binary_search(docrlst%nr,docrlst%idx(1:docrlst%nr),k,kk)
          kk = min(max(1,kk),docrlst%nr)
          if (k.ne.docrlst%idx(kk)) then
            if (wt_vec(8).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for Dinner-Ulmschneider polypeptide concerted rotations with &
 &omega sampling for residue ',k,', which is not eligible as a final residue in a stretch. Ignored.'
            end if
          else if (wt_vec(8).ne.docrlst%wt(kk)) then
            touched(8) = .true.
            docrlst%wt(kk) = wt_vec(8)
          end if
        end if
        if ((wt_vec(9).ge.0.0).AND.(nuccrlst%nr.gt.0)) then
          call binary_search(nuccrlst%nr,nuccrlst%idx(1:nuccrlst%nr),k,kk)
          kk = min(max(1,kk),nuccrlst%nr)
          if (k.ne.nuccrlst%idx(kk)) then
            if (wt_vec(9).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for Dinner-Ulmschneider nucleic acid concerted rotations for &
 &residue ',k,', which is not eligible as a final residue in a stretch. Ignored.'
            end if
          else if (wt_vec(9).ne.nuccrlst%wt(kk)) then
            touched(9) = .true.
            nuccrlst%wt(kk) = wt_vec(9)
          end if
        end if
        if ((wt_vec(10).ge.0.0).AND.(ujcrlst%nr.gt.0)) then
          call binary_search(ujcrlst%nr,ujcrlst%idx(1:ujcrlst%nr),k,kk)
          kk = min(max(1,kk),ujcrlst%nr)
          if (k.ne.ujcrlst%idx(kk)) then
            if (wt_vec(10).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for Ulmschneider-Jorgensen polypeptide concerted rotations for &
 &residue ',k,', which is not eligible as a final residue in a stretch. Ignored.'
            end if
          else if (wt_vec(10).ne.ujcrlst%wt(kk)) then
            touched(10) = .true.
            ujcrlst%wt(kk) = wt_vec(10)
          end if
        end if
        if ((wt_vec(11).ge.0.0).AND.(sjcrlst%nr.gt.0)) then
          call binary_search(sjcrlst%nr,sjcrlst%idx(1:sjcrlst%nr),k,kk)
          kk = min(max(1,kk),sjcrlst%nr)
          if (k.ne.sjcrlst%idx(kk)) then
            if (wt_vec(11).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for Sjunnesson-Favrin polypeptide concerted rotations for &
 &residue ',k,', which is not eligible as a final residue in a stretch. Ignored.'
            end if
          else if (wt_vec(11).ne.sjcrlst%wt(kk)) then
            touched(11) = .true.
            sjcrlst%wt(kk) = wt_vec(11)
          end if
        end if
      else if ((inmllp.EQV..true.).AND.((k.le.nmol).AND.(k.ge.1))) then
        str3 = str2(ii:MAXSTRLEN)
        wt_vec(:) = 1.0
        read(str3,*,iostat=iomess) (wt_vec(i), i=1,2) ! truncation is tolerated 
        if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing preferential sampling&
 & input (got: ',prefsamplingfile(t1:t2),').'
          call fexit()
        end if
!       go through each of the 2 molecule lists, locate correct entry, and set weight to target
        if ((wt_vec(1).ge.0.0).AND.(rblst%nr.gt.0)) then
          call binary_search(rblst%nr,rblst%idx(1:rblst%nr),k,kk)
          kk = min(max(1,kk),rblst%nr)
          if (k.ne.rblst%idx(kk)) then
            if (wt_vec(1).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for rigid-body sampling for &
 &molecule ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(1).ne.rblst%wt(kk)) then
            touched(12) = .true.
            rblst%wt(kk) = wt_vec(1)
          end if
        end if
        if ((wt_vec(2).ge.0.0).AND.(clurblst%nr.gt.0)) then
          call binary_search(clurblst%nr,clurblst%idx(1:clurblst%nr),k,kk)
          kk = min(max(1,kk),clurblst%nr)
          if (k.ne.clurblst%idx(kk)) then
            if (wt_vec(2).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for cluster rigid-body sampling for &
 &molecule ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(2).ne.clurblst%wt(kk)) then
            touched(13) = .true.
            clurblst%wt(kk) = wt_vec(2)
          end if
        end if
      else if ((inmtlp.EQV..true.).AND.((k.le.nmoltyp).AND.(k.ge.1))) then
        str3 = str2(ii:MAXSTRLEN)
        wt_vec(:) = 1.0
        read(str3,*,iostat=iomess) wt_vec(1)  
        if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing preferential sampling&
 & input (got: ',prefsamplingfile(t1:t2),').'
          call fexit()
        end if
!       go through the list, locate correct entry, and set weight to target
        if ((wt_vec(1).ge.0.0).AND.(fluclst%nr.gt.0)) then
          call binary_search(fluclst%nr,fluclst%idx(1:fluclst%nr),k,kk)
          kk = min(max(1,kk),fluclst%nr)
          if (k.ne.fluclst%idx(kk)) then
            if (wt_vec(1).ne.1.0) then
              write(ilog,*) 'Warning. Requested altered weight for particle insertion/deletion sampling for &
 &molecule type ',k,', which is not eligible for this move type. Ignored.'
            end if
          else if (wt_vec(1).ne.fluclst%wt(kk)) then
            touched(14) = .true.
            fluclst%wt(kk) = wt_vec(1)
          end if
        end if

      else if ((inatlp.EQV..true.).AND.((k.le.n).AND.(k.ge.1))) then
        str3 = str2(ii:MAXSTRLEN)
        wt_vec(:) = 1.0
        read(str3,*,iostat=iomess) wt_vec(1)
        if (iomess.eq.2) then
          write(ilog,*) 'Fatal. File I/O error while processing preferential sampling&
 & input (got: ',prefsamplingfile(t1:t2),').'
          call fexit()
        end if
!       go through the list, locate correct entry, and set weight to target
        if ((izrot(k)%alsz.le.0).OR.(k.eq.atmol(molofrs(atmres(k)),1))) then
          write(ilog,*) 'Warning. Requested altered weight for OTHER dihedral angle sampling for &
 &Z-matrix coordinate defined by atom ',k,', which is not rotatable or otherwise ineligible. Ignored.'
        else
          inother = .false.
          if (wt_vec(1).ge.0.0) then
!           any given dihedral should occur in at most one of these lists
            if (unklst%nr.gt.0) then
              call binary_search(unklst%nr,unklst%idx(1:unklst%nr),k,kk)
              kk = min(max(1,kk),unklst%nr)
              if (k.eq.unklst%idx(kk)) inother = .true.
              if ((k.eq.unklst%idx(kk)).AND.(wt_vec(1).ne.unklst%wt(kk))) then
                touched(15) = .true.
                unklst%wt(kk) = wt_vec(1)
              end if
            end if
            if (unslst%nr.gt.0) then
              call binary_search(unslst%nr,unslst%idx(1:unslst%nr),k,kk)
              kk = min(max(1,kk),unslst%nr)
              if (k.eq.unslst%idx(kk)) inother = .true.
              if ((k.eq.unslst%idx(kk)).AND.(wt_vec(1).ne.unslst%wt(kk))) then
                touched(16) = .true.
                unslst%wt(kk) = wt_vec(1)
              end if
            end if
            if (natlst%nr.gt.0) then
              call binary_search(natlst%nr,natlst%idx(1:natlst%nr),k,kk)
              kk = min(max(1,kk),natlst%nr)
              if (k.eq.natlst%idx(kk)) inother = .true.
              if ((k.eq.natlst%idx(kk)).AND.(wt_vec(1).ne.natlst%wt(kk))) then
                touched(17) = .true.
                natlst%wt(kk) = wt_vec(1)
              end if
            end if
          end if
          if ((inother.EQV..false.).AND.(wt_vec(1).ne.1.0)) then
            write(ilog,*) 'Warning. Requested altered weight for OTHER dihedral angle sampling for &
 &Z-matrix coordinate defined by atom ',k,', which is not eligible for this move type. Ignored.'
          end if
        end if
      else
        write(ilog,*) 'Warning. Unprocessed entry in preferential sampling input file. &
 &Check for unwanted or incorrect entries or missing section headers.'
      end if
    end if
    str2(:) = ' '
    read(iunit,79,iostat=iomess) str2
    call strlims_quiet(str2,st1,st2)
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing preferential sampling&
 & input (got: ',prefsamplingfile(t1:t2),').'
      call fexit()
    end if
  end do
!
! now re-compute the weight boundaries
!
  if (fylst%nr.gt.0) then
    if (sum(fylst%wt(1:fylst%nr)).le.0.0) then
      if (have_pivot.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, phi/psi sampling has no eligible &
 &residues with finite weight left. Adjust input file or disable phi/psi moves altogether.'
        call fexit()
      else
        fylst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(fylst%wt(1:fylst%nr))
    fylst%wt(1) = ivwt*fylst%wt(1)
    do i=2,fylst%nr-1
      fylst%wt(i) = fylst%wt(i-1)+fylst%wt(i)*ivwt
    end do
    fylst%wt(fylst%nr) = 1.0
  end if
  if (wlst%nr.gt.0) then
    if (sum(wlst%wt(1:wlst%nr)).le.0.0) then
      if (have_omega.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, omega sampling has no eligible &
 &residues with finite weight left. Adjust input file or disable omega moves altogether.'
        call fexit()
      else
        wlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(wlst%wt(1:wlst%nr))
    wlst%wt(1) = ivwt*wlst%wt(1)
    do i=2,wlst%nr-1
      wlst%wt(i) = wlst%wt(i-1)+wlst%wt(i)*ivwt
    end do
    wlst%wt(wlst%nr) = 1.0
  end if
  if (chilst%nr.gt.0) then
    if (sum(chilst%wt(1:chilst%nr)).le.0.0) then
      if ((have_chi.EQV..true.).AND.(phfreq.lt.1.0)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, chi sampling has no eligible &
 &residues with finite weight left. Adjust input file or disable chi moves altogether.'
        call fexit()
      else
        chilst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(chilst%wt(1:chilst%nr))
    chilst%wt(1) = ivwt*chilst%wt(1)
    do i=2,chilst%nr-1
      chilst%wt(i) = chilst%wt(i-1)+chilst%wt(i)*ivwt
    end do
    chilst%wt(chilst%nr) = 1.0
  end if
  if (nuclst%nr.gt.0) then
    if (sum(nuclst%wt(1:nuclst%nr)).le.0.0) then
      if ((have_nuc.EQV..true.).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, NUC pivot sampling has no eligible &
 &residues with finite weight left. Adjust input file or disable NUC pivot moves altogether.'
        call fexit()
      else
        nuclst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(nuclst%wt(1:nuclst%nr))
    nuclst%wt(1) = ivwt*nuclst%wt(1)
    do i=2,nuclst%nr-1
      nuclst%wt(i) = nuclst%wt(i-1)+nuclst%wt(i)*ivwt
    end do
    nuclst%wt(nuclst%nr) = 1.0
  end if
  if (puclst%nr.gt.0) then
    if (sum(puclst%wt(1:puclst%nr)).le.0.0) then
      if (have_pucker.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, polypeptide pucker sampling has no eligible &
 &residues with finite weight left. Adjust input file or disable polypeptide pucker moves altogether.'
        call fexit()
      else
        puclst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(puclst%wt(1:puclst%nr))
    puclst%wt(1) = ivwt*puclst%wt(1)
    do i=2,puclst%nr-1
      puclst%wt(i) = puclst%wt(i-1)+puclst%wt(i)*ivwt
    end do
    puclst%wt(puclst%nr) = 1.0
  end if
  if (nucpuclst%nr.gt.0) then
    if (sum(nucpuclst%wt(1:nucpuclst%nr)).le.0.0) then
      if (have_nucpuck.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, nucleic acid pucker sampling has no eligible &
 &residues with finite weight left. Adjust input file or disable nucleic acid pucker moves altogether.'
        call fexit()
      else
        nucpuclst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(nucpuclst%wt(1:nucpuclst%nr))
    nucpuclst%wt(1) = ivwt*nucpuclst%wt(1)
    do i=2,nucpuclst%nr-1
      nucpuclst%wt(i) = nucpuclst%wt(i-1)+nucpuclst%wt(i)*ivwt
    end do
    nucpuclst%wt(nucpuclst%nr) = 1.0
  end if
  if (djcrlst%nr.gt.0) then
    if (sum(djcrlst%wt(1:djcrlst%nr)).le.0.0) then
      if (have_djcr.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights,  Dinner-Ulmschneider polypeptide concerted rotation &
 &moves without omega sampling have no residues with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        djcrlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(djcrlst%wt(1:djcrlst%nr))
    djcrlst%wt(1) = ivwt*djcrlst%wt(1)
    do i=2,djcrlst%nr-1
      djcrlst%wt(i) = djcrlst%wt(i-1)+djcrlst%wt(i)*ivwt
    end do
    djcrlst%wt(djcrlst%nr) = 1.0
  end if
  if (docrlst%nr.gt.0) then
    if (sum(docrlst%wt(1:docrlst%nr)).le.0.0) then
      if (have_docr.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, Dinner-Ulmschneider polypeptide concerted rotation &
 &moves with omega sampling have no residues with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        docrlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(docrlst%wt(1:docrlst%nr))
    docrlst%wt(1) = ivwt*docrlst%wt(1)
    do i=2,docrlst%nr-1
      docrlst%wt(i) = docrlst%wt(i-1)+docrlst%wt(i)*ivwt
    end do
    docrlst%wt(docrlst%nr) = 1.0
  end if
  if (nuccrlst%nr.gt.0) then
    if (sum(nuccrlst%wt(1:nuccrlst%nr)).le.0.0) then
      if (have_nuccr.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, Dinner-Ulmschneider nucleic acid concerted rotation &
 &moves have no residues with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        nuccrlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(nuccrlst%wt(1:nuccrlst%nr))
    nuccrlst%wt(1) = ivwt*nuccrlst%wt(1)
    do i=2,nuccrlst%nr-1
      nuccrlst%wt(i) = nuccrlst%wt(i-1)+nuccrlst%wt(i)*ivwt
    end do
    nuccrlst%wt(nuccrlst%nr) = 1.0
  end if
  if (ujcrlst%nr.gt.0) then
    if (sum(ujcrlst%wt(1:ujcrlst%nr)).le.0.0) then
      if (have_ujcr.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, Ulmschneider-Jorgensen bond angle concerted rotation &
 &moves for polypeptides have no residues with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        ujcrlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(ujcrlst%wt(1:ujcrlst%nr))
    ujcrlst%wt(1) = ivwt*ujcrlst%wt(1)
    do i=2,ujcrlst%nr-1
      ujcrlst%wt(i) = ujcrlst%wt(i-1)+ujcrlst%wt(i)*ivwt
    end do
    ujcrlst%wt(ujcrlst%nr) = 1.0
  end if
  if (sjcrlst%nr.gt.0) then
    if (sum(sjcrlst%wt(1:sjcrlst%nr)).le.0.0) then
      if (have_sjcr.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, Sjunnesson-Favrin polypeptide concerted rotation &
 &moves have no eligible residues with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        sjcrlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(sjcrlst%wt(1:sjcrlst%nr))
    sjcrlst%wt(1) = ivwt*sjcrlst%wt(1)
    do i=2,sjcrlst%nr-1
      sjcrlst%wt(i) = sjcrlst%wt(i-1)+sjcrlst%wt(i)*ivwt
    end do
    sjcrlst%wt(sjcrlst%nr) = 1.0
  end if
  if (rblst%nr.gt.0) then
    if (sum(rblst%wt(1:rblst%nr)).le.0.0) then
      if ((have_rigid.EQV..true.).AND.(clurb_freq.lt.1.0)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, single molecule, rigid-body moves &
 &have no eligible molecules with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        rblst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(rblst%wt(1:rblst%nr))
    rblst%wt(1) = ivwt*rblst%wt(1)
    do i=2,rblst%nr-1
      rblst%wt(i) = rblst%wt(i-1)+rblst%wt(i)*ivwt
    end do
    rblst%wt(rblst%nr) = 1.0
  end if
  if (clurblst%nr.gt.0) then
    if (sum(clurblst%wt(1:clurblst%nr)).le.0.0) then
      if (have_clurb.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, cluster rigid-body moves &
 &have no eligible molecules with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        clurblst%wt(:) = 1.0
      end if
    end if
    if (sum(clurblst%wt(1:clurblst%nr)).eq.maxval(clurblst%wt(1:clurblst%nr))) then
      if (have_clurb.EQV..true.) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, cluster rigid-body moves &
 &have only a single molecule with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        clurblst%wt(:) = 1.0
      end if
    end if
    k = 0
    do i=1,clurblst%nr
      if (clurblst%wt(i).le.0.0) k = k + 1
    end do
    if ((clurblst%nr-k).lt.clurb_maxsz) then
      write(ilog,*)' Warning. After application of preferential sampling weights, maximum cluster size for cluster rigid &
 &body moves is reduced to the remaining number of eligible molecules (',clurblst%nr-k,').'
      clurb_maxsz = clurblst%nr - k
    end if
    ivwt = 1.0/sum(clurblst%wt(1:clurblst%nr))
    clurblst%wt(1) = ivwt*clurblst%wt(1)
    do i=2,clurblst%nr-1
      clurblst%wt(i) = clurblst%wt(i-1)+clurblst%wt(i)*ivwt
    end do
    clurblst%wt(clurblst%nr) = 1.0
  end if
  if (fluclst%nr.gt.0) then
    if (sum(fluclst%wt(1:fluclst%nr)).le.0.0) then
      if ((have_particlefluc.EQV..true.).AND.(ens%flag.eq.5)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, particle insertion/deletion moves &
 &have no eligible molecule types with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        fluclst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(fluclst%wt(1:fluclst%nr))
    fluclst%wt(1) = ivwt*fluclst%wt(1)
    do i=2,fluclst%nr-1
      fluclst%wt(i) = fluclst%wt(i-1)+fluclst%wt(i)*ivwt
    end do
    fluclst%wt(fluclst%nr) = 1.0
  end if
  if (unklst%nr.gt.0) then
    if (sum(unklst%wt(1:unklst%nr)).le.0.0) then
      if ((have_other.EQV..true.).AND.(have_unkother.EQV..true.)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, OTHER sampling for unsupported residues has no &
 &eligible degrees of freedom with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        unklst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(unklst%wt(1:unklst%nr))
    unklst%wt(1) = ivwt*unklst%wt(1)
    do i=2,unklst%nr-1
      unklst%wt(i) = unklst%wt(i-1)+unklst%wt(i)*ivwt
    end do
    unklst%wt(unklst%nr) = 1.0
  end if
  if (unslst%nr.gt.0) then
    if (sum(unslst%wt(1:unslst%nr)).le.0.0) then
      if ((have_other.EQV..true.).AND.(have_unsother.EQV..true.)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, OTHER sampling for unsupported degrees of freedom in &
 &supported residues has no eligible angles with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        unslst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(unslst%wt(1:unslst%nr))
    unslst%wt(1) = ivwt*unslst%wt(1)
    do i=2,unslst%nr-1
      unslst%wt(i) = unslst%wt(i-1)+unslst%wt(i)*ivwt
    end do
    unslst%wt(unslst%nr) = 1.0
  end if
  if (natlst%nr.gt.0) then
    if (sum(natlst%wt(1:natlst%nr)).le.0.0) then
      if ((have_other.EQV..true.).AND.(have_natother.EQV..true.)) then
        write(ilog,*) 'Fatal. After applying preferential sampling weights, OTHER sampling for natively supported degrees of &
 &freedom has no dihedral angles with finite weight left. Adjust input file or disable these moves altogether.'
        call fexit()
      else
        natlst%wt(:) = 1.0
      end if
    end if
    ivwt = 1.0/sum(natlst%wt(1:natlst%nr))
    natlst%wt(1) = ivwt*natlst%wt(1)
    do i=2,natlst%nr-1
      natlst%wt(i) = natlst%wt(i-1)+natlst%wt(i)*ivwt
    end do
    natlst%wt(natlst%nr) = 1.0
  end if
!
  close(unit=iunit)
!
 555 format(i12,': ',g12.5,'%')
 556 format('Residue        Sampling Weight')
 557 format('Molecule       Sampling Weight')
 558 format('Molecule Type  Sampling Weight')
 560 format('Z-matrix Entry Sampling Weight')


  if (psw_report.EQV..true.) then
    if ((have_particlefluc.EQV..true.).AND.(touched(14).EQV..true.).AND.(ens%flag.ne.6)) then
      write(ilog,*) '--------------------------------------- --------'
      write(ilog,*) '-- Preferential Sampling Weights for GC moves --'
      write(ilog,*) '------------------------------------------------'
      write(ilog,558)
      write(ilog,555) fluclst%idx(1),fluclst%wt(1)*100.0
      do i=2,fluclst%nr
        write(ilog,555) fluclst%idx(i),(fluclst%wt(i)-fluclst%wt(i-1))*100.0
      end do
    end if
    if ((have_rigid.EQV..true.).AND.(clurb_freq.lt.1.0).AND.(touched(12).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for RIGID moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,557)
      write(ilog,555) rblst%idx(1),rblst%wt(1)*100.0
      do i=2,rblst%nr
        write(ilog,555) rblst%idx(i),(rblst%wt(i)-rblst%wt(i-1))*100.0
      end do
    end if
    if ((have_clurb.EQV..true.).AND.(touched(13).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for CLURB moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,557)
      write(ilog,555) clurblst%idx(1),clurblst%wt(1)*100.0
      do i=2,clurblst%nr
        write(ilog,555) clurblst%idx(i),(clurblst%wt(i)-clurblst%wt(i-1))*100.0
      end do
    end if
    if ((have_pivot.EQV..true.).AND.(touched(1).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for PIVOT moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) fylst%idx(1),fylst%wt(1)*100.0
      do i=2,fylst%nr
        write(ilog,555) fylst%idx(i),(fylst%wt(i)-fylst%wt(i-1))*100.0
      end do
    end if
    if ((have_omega.EQV..true.).AND.(touched(2).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for OMEGA moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) wlst%idx(1),wlst%wt(1)*100.0
      do i=2,wlst%nr
        write(ilog,555) wlst%idx(i),(wlst%wt(i)-wlst%wt(i-1))*100.0
      end do
    end if
    if ((have_chi.EQV..true.).AND.(phfreq.lt.1.0).AND.(touched(3).EQV..true.)) then
      write(ilog,*) '-------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for CHI moves --'
      write(ilog,*) '-------------------------------------------------'
      write(ilog,556)
      write(ilog,555) chilst%idx(1),chilst%wt(1)*100.0
      do i=2,chilst%nr
        write(ilog,555) chilst%idx(i),(chilst%wt(i)-chilst%wt(i-1))*100.0
      end do
    end if
    if ((have_nuc.EQV..true.).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0).AND.(touched(4).EQV..true.)) then
      write(ilog,*) '-------------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for NUC-PIVOT moves --'
      write(ilog,*) '-------------------------------------------------------'
      write(ilog,556)
      write(ilog,555) nuclst%idx(1),nuclst%wt(1)*100.0
      do i=2,nuclst%nr
        write(ilog,555) nuclst%idx(i),(nuclst%wt(i)-nuclst%wt(i-1))*100.0
      end do
    end if
    if ((have_pucker.EQV..true.).AND.(touched(5).EQV..true.)) then
      write(ilog,*) '----------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for PUCKER moves --'
      write(ilog,*) '----------------------------------------------------'
      write(ilog,556)
      write(ilog,555) puclst%idx(1),puclst%wt(1)*100.0
      do i=2,puclst%nr
        write(ilog,555) puclst%idx(i),(puclst%wt(i)-puclst%wt(i-1))*100.0
      end do
    end if
    if ((have_nucpuck.EQV..true.).AND.(touched(6).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for SUGAR moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) nucpuclst%idx(1),nucpuclst%wt(1)*100.0
      do i=2,nucpuclst%nr
        write(ilog,555) nucpuclst%idx(i),(nucpuclst%wt(i)-nucpuclst%wt(i-1))*100.0
      end do
    end if
    if ((have_djcr.EQV..true.).AND.(touched(7).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for DJ-CR moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) djcrlst%idx(1),djcrlst%wt(1)*100.0
      do i=2,djcrlst%nr
        write(ilog,555) djcrlst%idx(i),(djcrlst%wt(i)-djcrlst%wt(i-1))*100.0
      end do
    end if
    if ((have_docr.EQV..true.).AND.(touched(8).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for DO-CR moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) docrlst%idx(1),docrlst%wt(1)*100.0
      do i=2,docrlst%nr
        write(ilog,555) docrlst%idx(i),(docrlst%wt(i)-docrlst%wt(i-1))*100.0
      end do
    end if
    if ((have_nuccr.EQV..true.).AND.(touched(9).EQV..true.)) then
      write(ilog,*) '----------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for NUC-CR moves --'
      write(ilog,*) '----------------------------------------------------'
      write(ilog,556)
      write(ilog,555) nuccrlst%idx(1),nuccrlst%wt(1)*100.0
      do i=2,nuccrlst%nr
        write(ilog,555) nuccrlst%idx(i),(nuccrlst%wt(i)-nuccrlst%wt(i-1))*100.0
      end do
    end if
    if ((have_ujcr.EQV..true.).AND.(touched(10).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for UJ-CR moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) ujcrlst%idx(1),ujcrlst%wt(1)*100.0
      do i=2,ujcrlst%nr
        write(ilog,555) ujcrlst%idx(i),(ujcrlst%wt(i)-ujcrlst%wt(i-1))*100.0
      end do
    end if
    if ((have_sjcr.EQV..true.).AND.(touched(11).EQV..true.)) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for SJ-CR moves --'
      write(ilog,*) '---------------------------------------------------'
      write(ilog,556)
      write(ilog,555) sjcrlst%idx(1),sjcrlst%wt(1)*100.0
      do i=2,sjcrlst%nr
        write(ilog,555) sjcrlst%idx(i),(sjcrlst%wt(i)-sjcrlst%wt(i-1))*100.0
      end do
    end if
    if ((have_other.EQV..true.).AND.((touched(15).EQV..true.).OR.(touched(16).EQV..true.).OR.(touched(17).EQV..true.))) then
      write(ilog,*) '---------------------------------------------------'
      write(ilog,*) '-- Preferential Sampling Weights for OTHER moves --'
      write(ilog,*) '---------------------------------------------------'
      if ((have_unkother.EQV..true.).AND.(touched(15).EQV..true.)) then
        write(ilog,*)'-- List for unsupported residues:'
        write(ilog,560)
        write(ilog,555) unklst%idx(1),unklst%wt(1)*100.0
        do i=2,unklst%nr
          write(ilog,555) unklst%idx(i),(unklst%wt(i)-unklst%wt(i-1))*100.0
        end do
      end if
      if ((have_unsother.EQV..true.).AND.(touched(16).EQV..true.)) then
        write(ilog,*)'-- List for unsupported degrees of freedom in supported residues:'
        write(ilog,560)
        write(ilog,555) unslst%idx(1),unslst%wt(1)*100.0
        do i=2,unslst%nr
          write(ilog,555) unslst%idx(i),(unslst%wt(i)-unslst%wt(i-1))*100.0
        end do
      end if
      if ((have_natother.EQV..true.).AND.(touched(17).EQV..true.)) then
        write(ilog,*)'-- List for natively supported degrees of freedom:'
        write(ilog,560)
        write(ilog,555) natlst%idx(1),natlst%wt(1)*100.0
        do i=2,natlst%nr
          write(ilog,555) natlst%idx(i),(natlst%wt(i)-natlst%wt(i-1))*100.0
        end do
      end if
    end if
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_fegfile()
!
  use energies
  use iounit
  use molecule
  use sequen
  use params
  use aminos
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,ii,i,imol,mn,rs,mt,decfl
  logical feglst(nseq),declst(nseq)
  character(10) abc
  character(MAXSTRLEN) str2,str3
  logical exists,havem
!
! some sanity checks before entering the parameter setup routines 
!
  if ((use_hardsphere.EQV..true.).AND.((use_FEGS(1).EQV..true.).OR.(use_FEGS(3).EQV..true.))) then
    write(ilog,*) 'Fatal. Ghosting supports only soft-core excluded &
 &volume. Hard-sphere potentials will never work for FEG.'
    call fexit()
  else if ((nindx.ne.12).AND.((use_FEGS(1).EQV..true.).OR.(use_FEGS(3).EQV..true.))) then
    write(ilog,*) 'Fatal. Ghosting requires the IPP-exponent to be e&
 &xactly 12.'
    call fexit()
  end if
!
  if (use_FEGS(1).EQV..true.) then
    call setup_parFEG(1)
  end if
  if (use_FEGS(3).EQV..true.) then
    call setup_parFEG(3)
  end if
  if (use_FEGS(6).EQV..true.) then
    call setup_parFEG(6)
  end if
!
  feglst(:) = .false.
  declst(:) = .false.
!
  call strlims(fegfile,t1,t2)
  inquire(file=fegfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &fegfile(t1:t2),') for FMCSC_FEGFILE. Standard potentials used.'
    use_FEG = .false.
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=fegfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for ghost part&
 &icle input (',fegfile(t1:t2),'). Turning off terms.'
    use_FEG = .false.
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing ghost part&
 &icle input (got: ',fegfile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_str(str2,abc,ii)
!
! first molecule-based
!
  if ((abc(1:1).eq.'M').OR.(abc(1:1).eq.'m')) then
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing ghost &
 &particle input (got: ',fegfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      mn = -1
      call extract_int(str2,mn,ii)
      if ((mn.gt.nmol).OR.(mn.le.0)) then
        write(ilog,*) 'Requested ghosting for molecule ',&
 &mn,', which does not exist. Ignored ...'
        cycle
      end if
      if (ii.eq.1) exit
!
      do rs=rsmol(mn,1),rsmol(mn,2)
        feglst(rs) = .true.
      end do
!
!     attempt to read flag for full decoupling overwrite
      str3 = str2(ii:MAXSTRLEN)
      ii = 1
      decfl = -1
      call extract_int(str3,decfl,ii)
      if (ii.eq.1) cycle ! no extra integer found
      if (decfl.eq.0) then
       do rs=rsmol(mn,1),rsmol(mn,2)
          declst(rs) = .true.
        end do
      else
        write(ilog,*) 'Incomprehensible input for extra flag in ghos&
 &t particle input. Ignored ...'
      end if
!
    end do
!
! and now molecule-type-based
!
  else if ((abc(1:1).eq.'T').OR.(abc(1:1).eq.'t')) then
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing ghost &
 &particle input (got: ',fegfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      mt = -1
      call extract_int(str2,mt,ii)
      if ((mt.gt.nmoltyp).OR.(mt.le.0)) then
        write(ilog,*) 'Requested ghosting for molecule type',&
 &mt,', which does not exist. Ignored ...'
        cycle
      end if
      if (ii.eq.1) exit
!
      do imol=1,nmol
        if (moltypid(imol).eq.mt) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            feglst(rs) = .true.
          end do
        end if
      end do
!
!     attempt to read flag for full decoupling overwrite
      str3 = str2(ii:MAXSTRLEN)
      ii = 1
      decfl = -1
      call extract_int(str3,decfl,ii)
      if (ii.eq.1) cycle ! no extra integer found
      if (decfl.eq.0) then
        do imol=1,nmol
          if (moltypid(imol).eq.mt) then
            do rs=rsmol(imol,1),rsmol(imol,2)
              declst(rs) = .true.
            end do
          end if
        end do
      else
        write(ilog,*) 'Incomprehensible input for extra flag in ghos&
 &t particle input. Ignored ...'
      end if
!
    end do
!
  else if ((abc(1:1).eq.'R').OR.(abc(1:1).eq.'r')) then
    do while(.true.)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing ghost &
 &particle input (got: ',fegfile(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      rs = -1
      call extract_int(str2,rs,ii)
      if ((rs.gt.nseq).OR.(rs.le.0)) then
        write(ilog,*) 'Requested ghosting for residue',&
 &rs,', which does not exist. Ignored ...'
        cycle
      end if
      if (ii.eq.1) exit
!
      feglst(rs) = .true.
!
!     attempt to read flag for full decoupling overwrite
      str3 = str2(ii:MAXSTRLEN)
      ii = 1
      decfl = -1
      call extract_int(str3,decfl,ii)
      if (ii.eq.1) cycle ! no extra integer found
      if (decfl.eq.0) then
        declst(rs) = .true.
      else
        write(ilog,*) 'Incomprehensible input for extra flag in ghos&
 &t particle input. Ignored ...'
      end if
!
    end do
!
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',fegfile(t1:t2),': ',abc(1:1),'! Fatal exit.'
    close(unit=iunit)
    call fexit()
  end if
  close(unit=iunit)
!
  havem = .false.
  do rs=1,nseq
    if (feglst(rs).EQV..true.) then
      havem = .true.
      exit
    end if
  end do
!
  if (havem.EQV..false.) then
    write(ilog,*) 'Warning. No ghost molecules or residues requested&
 & in ',fegfile(t1:t2),'? Turning off potentials ...'
    use_FEG = .false.
    return 
  end if
!
! now we have finalized that we do indeed use FEG
  allocate(par_FEG(nseq))
  allocate(par_FEG3(nseq))
  do i=1,nseq
    par_FEG(i) = feglst(i)
    par_FEG3(i) = declst(i)
  end do
!
 45   format(4x,'Residue #',i6,' (',a3,')')
 46   format(4x,'Residue #',i6,' (',a3,') <-- always fully de-coupled')
!
  if (feg_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '-- List of ghost residues -----------------------&
 &---'
    do rs=1,nseq
      if (par_FEG3(rs).EQV..true.) then
        if (par_FEG(rs).EQV..true.) then
          write(ilog,46) rs,amino(seqtyp(rs)) 
        end if
      else
        if (par_FEG(rs).EQV..true.) then
          write(ilog,45) rs,amino(seqtyp(rs)) 
        end if
      end if
    end do
    write(ilog,*) '-------------------------------------------------&
 &---'
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_tgrpfile()
!
  use energies
  use iounit
  use molecule
  use sequen
  use params
  use aminos
  use system
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,ii,i,j,jj,kk,mn
  character(10) abc
  character(MAXSTRLEN) str2
  integer, ALLOCATABLE:: tmpv(:),tmpv2(:)
  logical exists,novel
!
! some sanity checks
!
  if (dyn_mode.eq.1) then
    write(ilog,*) 'Fatal. Called read_tgrpfile() while not performin&
 &g a dynamics calculation. This is a bug.'
    call fexit()
  end if
!
! for Andersen, Langevin, BD don't bother (invidual coupling anyway)
  if ((tstat%flag.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.4).OR.(dyn_mode.eq.8)) then
    return
  end if
!
  call strlims(tstat%fnam,t1,t2)
  inquire(file=tstat%fnam(t1:t2),exist=exists)
  if ((exists.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.5))&
 &.AND.(ens%flag.eq.1).AND.((tstat%flag.eq.1).OR.(tstat%flag.eq.4))) then
!   only throw warning for Berendsen and Newtonian MD
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &tstat%fnam(t1:t2),') for FMCSC_TSTAT_FILE. Assuming a single tempe&
 &rature coupling group encompassing all molecules.'
!   note that in allocate the Tstat is set up such that by default the
!   complete system will form a single T-coupling group
    return
  else if (exists.EQV..false.) then
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=tstat%fnam(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for T-coupling&
 & group input (',tstat%fnam(t1:t2),'). Assuming a single group enco&
 &mpassing all molecules.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing T-coupling&
 & group input (got: ',tstat%fnam(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_str(str2,abc,ii)
!
! first molecule-based (this file has a hard format: nmol+1 lines are required after the header!
!
  allocate(tmpv(nmol+1))
  if ((abc(1:1).eq.'M').OR.(abc(1:1).eq.'m')) then
    do i=1,nmol+1
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Fatal. Input file for T-coupling groups (mode&
 &s M or m) is incomplete. Make sure there are nmol+1 rows after th&
 &e header line (',tstat%fnam(t1:t2),').'
        close(unit=iunit)
        deallocate(tmpv)
        call fexit()
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing T-coupling&
 & group input (got: ',tstat%fnam(t1:t2),').'
        deallocate(tmpv)
        call fexit()
      end if
!
      ii = 1
      tmpv(i) = -48571903
      call extract_int(str2,tmpv(i),ii)
!
      if (ii.eq.1) then
        write(ilog,*) 'Fatal. Input file for T-coupling groups (mode&
 &s M or m) is incomplete. Make sure there are nmol+1 rows after th&
 &e header line (',tstat%fnam(t1:t2),').'
        close(unit=iunit)
        deallocate(tmpv)
        call fexit()
      end if
!
    end do
!
! and now molecule-type-based
!
  else if ((abc(1:1).eq.'T').OR.(abc(1:1).eq.'t')) then
    do i=1,nmoltyp+1
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Fatal. Input file for T-coupling groups (mode&
 &s T or t) is incomplete. Make sure there are nmoltyp+1 rows after&
 & the header line (',tstat%fnam(t1:t2),').'
        close(unit=iunit)
        deallocate(tmpv)
        call fexit()
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing T-coupling&
 & group input (got: ',tstat%fnam(t1:t2),').'
        deallocate(tmpv)
        call fexit()
      end if
!
      ii = 1
      mn = -48571903
      call extract_int(str2,mn,ii)
!
      if (ii.eq.1) then
        write(ilog,*) 'Fatal. Input file for T-coupling groups (mode&
 &s T or t) is incomplete. Make sure there are nmoltyp+1 rows after&
 & the header line (',tstat%fnam(t1:t2),').'
        close(unit=iunit)
        deallocate(tmpv)
        call fexit()
      end if
!
      if (i.le.nmoltyp) then
        do j=1,nmol
          if (moltypid(j).eq.i) then
            tmpv(j) = mn
          end if
        end do
      else
        tmpv(nmol+1) = mn
      end if
!
    end do
!
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',tstat%fnam(t1:t2),': ',abc(1:1),'! Fatal exit.'
    close(unit=iunit)
    deallocate(tmpv)
    call fexit()
  end if
  close(unit=iunit)
!
! now we need to check the consistency of the request
  do i=1,nmol+1
    if (tmpv(i).le.0) then
      write(ilog,*) 'Fatal. Illegal request for T-coupling group enc&
 &ountered first at molecule #',i,'. T-coupling group numbers have t&
 &o be positive integers.'
      deallocate(tmpv)
      call fexit()
    end if
  end do
  kk = 0
  allocate(tmpv2(nmol+1))
  tmpv2(:) = 0
  do i=1,nmol+1
    novel = .true.
    do j=1,kk
      if (tmpv(i).eq.tmpv2(j)) then
        novel = .false.
        exit
      end if
    end do
    if (novel.EQV..true.) then
      kk = kk + 1
      tmpv2(kk) = tmpv(i)
    end if
  end do
  tstat%n_tgrps = kk
!
  if (tstat%n_tgrps.eq.1) then
!   WARNING:: this might change in the future!!!
    if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
      write(ilog,*) 'The entire system forms a single T-coupling gro&
 &up.'
      write(ilog,*)
    else
      write(ilog,*) 'The entire system forms a single T-coupling gro&
 &up for analysis purposes.'
      write(ilog,*)
    end if
    deallocate(tmpv2)
    deallocate(tmpv)
    return
  end if
!
! if we have multiple groups, check and correct numbering
  call sort(tstat%n_tgrps,tmpv2)
!
  do j=1,tstat%n_tgrps
    if (tmpv2(j).ne.j) then
      write(ilog,*) 'Discovered discontinous numbering of T-coupling&
 & groups. Renumbering ...'
      do jj=j,tstat%n_tgrps
        do i=1,nmol+1
          if (tmpv(i).eq.tmpv2(j)) tmpv(i) = jj
        end do
      end do
      exit
    end if
  end do
!
! now transfer
  deallocate(tstat%grpdof)
  deallocate(tstat%grpT)
  allocate(tstat%grpT(tstat%n_tgrps))
  allocate(tstat%grpdof(tstat%n_tgrps))
  tstat%molgrp(:) = tmpv(:)
  deallocate(tmpv)
  deallocate(tmpv2)
!
 45   format(4x,'Molecule #',i6,' : ',i6)
 46   format(4x,'Boundary particle: ',i6)
!
  if (dyn_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '---List of molecules and their T-coupling group n&
 &umber ---'
    do i=1,nmol
      write(ilog,45) i,tstat%molgrp(i)
    end do
    if (ens%flag.eq.3) then
      write(ilog,46) tstat%molgrp(nmol+1)
    end if
    write(ilog,*) '-------------------------------------------------&
 &---------'
    write(ilog,*)
  end if
!
end
!
!-----------------------------------------------------------------------
!
!
! a generic function to read an atom index list into array tmpvec (at least size n)
!
subroutine read_cartconsfile(idxlst,nidx)
!
  use iounit
  use atoms
  use molecule
  use shakeetal
  use system
  use sequen
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,iomess,ii,jj,kk,imol,jmol,kmol,t1,t2,nidx
  integer idxlst(max(6*n-6,1),4),freeunit
  logical exists
  character(MAXSTRLEN) str2
  character(3) abc
!
  if (cart_cons_mode.ne.6) return
!
! some sanity checks
!
  if ((dyn_mode.eq.1).OR.(fycxyz.ne.2)) then
    write(ilog,*) 'Fatal. Called read_cartconsfile() while not performin&
 &g a Cartesian dynamics calculation. This is a bug.'
    call fexit()
  end if
!
  call strlims(cart_cons_file,t1,t2)
  inquire(file=cart_cons_file(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file (',&
 &cart_cons_file(t1:t2),') for FMCSC_SHAKEFILE despite requesting custom constraint &
 & set. Please check key-file and/or location/path/accessibility of input file.'
    call fexit()
  end if
!
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=cart_cons_file(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for custom constraints (SHAKE/LINCS)&
 & input (',cart_cons_file(t1:t2),'). Turning off constraints.'
    close(unit=iunit)
    cart_cons_mode = 1
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing custom constraints (SHAKE/LINCS)&
 & input (got: ',cart_cons_file(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_str(str2,abc,ii)
!
! first molecule-based (this file has a hard format: nmol+1 lines are required after the header!
!
  if ((abc(1:1).eq.'T').OR.(abc(1:1).eq.'t')) then
    do while (1.eq.1)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing SHAKE/LINCS&
 & input (got: ',cart_cons_file(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      kk = -1
      call extract_int(str2,kk,ii)
      if (ii.eq.1) exit
      if ((kk.gt.n).OR.(kk.le.0)) then
        write(ilog,*) 'Requested constraint for atom ',&
 &kk,', which does not exist. Ignored ...'
        cycle
      end if
      imol = molofrs(atmres(kk))
      if (imol.ne.moltyp(moltypid(imol),1)) then
        write(ilog,*) 'Fatal. Constraints in molecule type-resolved mode have to be&
 & provided only for the instances of atoms in the first molecule of each type.'
        call fexit()
      end if
      jj = -1
      call extract_int(str2,jj,ii)
      if (ii.eq.1) exit
      if ((jj.gt.n).OR.(jj.le.0)) then
        write(ilog,*) 'Requested constraint for atom ',&
 &jj,', which does not exist. Ignored ...'
        cycle
      end if
      jmol = molofrs(atmres(jj))
      if (jmol.ne.moltyp(moltypid(jmol),1)) then
        write(ilog,*) 'Fatal. Constraints in molecule type-resolved mode have to be&
 & provided only for the instances of atoms in the first molecule of each type.'
        call fexit()
      end if
      if (add_settle.EQV..true.) then
        if ((atmres(kk).eq.atmres(jj)).AND.((seqtyp(atmres(kk)).eq.39).OR.(seqtyp(atmres(kk)).eq.40))) then
          write(ilog,*) 'Requested constraint for (exemplary) atoms ',&
 &kk,' and ',jj,' which will be covered by SETTLE. Ignoring all such constraints ...'
          cycle
        end if
      end if
      if (imol.eq.jmol) then
        do kmol=1,nmol
          if (moltypid(kmol).eq.moltypid(imol)) then
            nidx = nidx + 1
            idxlst(nidx,1) = kk-atmol(imol,1)+atmol(kmol,1)
            idxlst(nidx,2) = jj-atmol(jmol,1)+atmol(kmol,1)
            idxlst(nidx,3) = 0
          end if
        end do
      else
        if ((moltyp(moltypid(imol),1).gt.1).OR.(moltyp(moltypid(jmol),1).gt.1)) then
          write(ilog,*) 'Fatal. Constraints spanning different molecules are ambiguous in &
 &molecule type-resolved input if more than molecule is present of either of the respective&
 & types. Use atom-resolved input mode instead.'
          call fexit()
        end if
        nidx = nidx + 1
        idxlst(nidx,1) = kk
        idxlst(nidx,2) = jj
        idxlst(nidx,3) = 0
      end if
!
    end do
 !
  else if ((abc(1:1).eq.'A').OR.(abc(1:1).eq.'a')) then
    do while (1.eq.1)
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        exit
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing SHAKE/LINCS&
 & input (got: ',cart_cons_file(t1:t2),').'
        call fexit()
      end if
!
      ii = 1
      kk = -1
      call extract_int(str2,kk,ii)
      if (ii.eq.1) exit
      if ((kk.gt.n).OR.(kk.le.0)) then
        write(ilog,*) 'Requested constraint for atom ',&
 &kk,', which does not exist. Ignored ...'
        cycle
      end if
      jj = -1
      call extract_int(str2,jj,ii)
      if (ii.eq.1) exit
      if ((jj.gt.n).OR.(jj.le.0)) then
        write(ilog,*) 'Requested constraint for atom ',&
 &jj,', which does not exist. Ignored ...'
        cycle
      end if
      if (add_settle.EQV..true.) then
        if ((atmres(kk).eq.atmres(jj)).AND.((seqtyp(atmres(kk)).eq.39).OR.(seqtyp(atmres(kk)).eq.40))) then
          write(ilog,*) 'Requested constraint for atoms ',&
 &kk,' and ',jj,' which will be covered by SETTLE. Ignored ...'
          cycle
        end if
      end if
      nidx = nidx + 1
      idxlst(nidx,1) = kk
      idxlst(nidx,2) = jj
      idxlst(nidx,3) = 0
!
    end do
!
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',cart_cons_file(t1:t2),': ',abc(1:1),'! Fatal exit.'
    close(unit=iunit)
    call fexit()
  end if
!
  close(unit=iunit)
!
end
!
!----------------------------------------------------------------------------
!
subroutine read_analysisgrpfile()
!
  use energies
  use iounit
  use molecule
  use sequen
  use params
  use aminos
  use system
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,ii,i,j,k,mn,jmol,refgrp,ngrps
  character(10) abc
  character(MAXSTRLEN) str2
  integer, ALLOCATABLE:: tmpv(:),mapgrp(:)
  logical exists,foundnew
!
  call strlims(angrpfile,t1,t2)
  inquire(file=angrpfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &angrpfile(t1:t2),') for FMCSC_ANGRPFILE. Assuming default assignment (&
 &per molecule type; all single-residue molecules are solvent).'
    return
  end if
!
 79   format(FORM_MAXSTRLEN)
!
  iunit = freeunit()
  open (unit=iunit,file=angrpfile(t1:t2),status='old')
!
  read(iunit,79,iostat=iomess) str2
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Warning. Got empty/incomplete file for analysis&
 & group input (',angrpfile(t1:t2),'). Assuming default assignment.'
    close(unit=iunit)
    return
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing analysis&
 & group input (got: ',angrpfile(t1:t2),').'
    call fexit()
  end if
!
  ii = 1
  call extract_str(str2,abc,ii)
!
! first molecule-based (this file has a hard format: nmol lines are required after the header!
!
  allocate(tmpv(nmol))
  allocate(mapgrp(nmol))
  ngrps = 0
  if ((abc(1:1).eq.'M').OR.(abc(1:1).eq.'m')) then
    do j=1,nmol
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Fatal. Input file for analysis groups (mode&
 &s M or m) is incomplete. Make sure there are as many rows as molecules after&
 & the header line (',angrpfile(t1:t2),').'
        close(unit=iunit)
        deallocate(mapgrp)
        deallocate(tmpv)
        call fexit()
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing analysis&
 & group input (got: ',angrpfile(t1:t2),').'
        deallocate(tmpv)
        deallocate(mapgrp)
        call fexit()
      end if
!
      ii = 1
      mn = -87981233
      call extract_int(str2,mn,ii)
!
      if (ii.eq.1) then
        write(ilog,*) 'Fatal. Input file for analysis groups (mode&
 &s M or m) is incomplete. Make sure there are as many rows as molecules after&
 & the header line (',angrpfile(t1:t2),').'
        close(unit=iunit)
        deallocate(mapgrp)
        deallocate(tmpv)
        call fexit()
      end if
!
      if (mn.ne.0) then
        tmpv(j) = mn
        foundnew = .true.
        do k=1,j-1
          if (mn.eq.tmpv(k)) then
            foundnew = .false.
            exit
          end if
        end do
        if (foundnew.EQV..true.) then
          ngrps = ngrps + 1
          mapgrp(ngrps) = mn
        end if
        if (mn.lt.0) then
          is_solvent(j) = .true.
        else
          is_solvent(j) = .false.
        end if
      else
        write(ilog,*) 'Fatal. Invalid entry for analysis groups (mode M or m)&
 &. Analysis group numbers have to be positive (solute) or negative (solvent) integers.'
        call fexit()
      end if
    end do
!
!   parse what we found
    do j=1,ngrps
      refgrp = mapgrp(j)
      jmol = 0
      do i=1,nmol
        if (tmpv(i).eq.refgrp) then
          if (jmol.eq.0) then
            jmol = moltypid(i)
          else
            if (moltypid(i).ne.jmol) then
              write(ilog,*) 'Fatal. Invalid entries for analysis groups (mode M or m)&
 &. Analysis groups have to be subsets of molecule types, i.e., cannot encompass molecules&
 & of different types:'
              write(ilog,*) 'Molecule #',i,' of type #',moltypid(i),' has group #',tmpv(i),' but that group also&
 & encompasses a molecule of type #',jmol,'.'
              call fexit()
            end if
          end if
        end if
      end do
    end do
!
!   assign: note we completely ignore the numbering the user provided
    do i=1,nmol
      do j=1,ngrps
        if (mapgrp(j).eq.tmpv(i)) then
          refgrp = j
          exit
        end if
      end do
      an_grp_mol(i) = refgrp
    end do
!
! this mode is only to support labeling of certain molecule types as solvent
  else if ((abc(1:1).eq.'T').OR.(abc(1:1).eq.'t')) then
    ngrps = nmoltyp
    do j=1,nmoltyp
      read(iunit,79,iostat=iomess) str2
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Fatal. Input file for analysis groups (mode&
 &s T or t) is incomplete. Make sure there are as many rows as molecule types after&
 & the header line (',angrpfile(t1:t2),').'
        close(unit=iunit)
        deallocate(mapgrp)
        deallocate(tmpv)
        call fexit()
      else if (iomess.eq.2) then
        write(ilog,*) 'Fatal. File I/O error while processing analysis&
 & group input (got: ',angrpfile(t1:t2),').'
        deallocate(mapgrp)
        deallocate(tmpv)
        call fexit()
      end if
!
      ii = 1
      mn = -87981233 
      call extract_int(str2,mn,ii)
!
      if (ii.eq.1) then
        write(ilog,*) 'Fatal. Input file for analysis groups (mode&
 &s T or t) is incomplete. Make sure there are as many rows as molecule types after&
 & the header line (',angrpfile(t1:t2),').'
        close(unit=iunit)
        deallocate(mapgrp)
        deallocate(tmpv)
        call fexit()
      end if
!
      if (mn.ne.0) then
        if (mn.lt.0) then
          do i=moltyp(j,1),nmol
            if (moltypid(i).eq.j) then
              is_solvent(i) = .true.
            end if
          end do
        else
         do i=moltyp(j,1),nmol
            if (moltypid(i).eq.j) then
              is_solvent(i) = .false.
            end if
          end do
        end if
      else
        write(ilog,*) 'Fatal. Invalid entry for analysis groups (mode T or t)&
 &. Analysis group numbers have to be positive (solute) or negative (solvent) integers.'
        call fexit()
      end if
    end do
!
!   assign: note we completely ignore the numbering the user provided
    do i=1,nmoltyp
      do j=1,nmol
        if (moltypid(j).eq.i) then
          an_grp_mol(j) = i
        end if
      end do
    end do
!
  else
    write(ilog,*) 'Unknown mode identifier while reading first line &
 &of file ',angrpfile(t1:t2),': ',abc(1:1),'! Fatal exit.'
    close(unit=iunit)
    deallocate(mapgrp)
    deallocate(tmpv)
    call fexit()
  end if
  close(unit=iunit)
!
  deallocate(tmpv)
  deallocate(mapgrp)
!
  nangrps = ngrps
!
  end
!
!--------------------------------------------------------------------------------------------
!
subroutine alloc_angrps(anidx,solvidx)
!
  use system
  use molecule
  use sequen
!
  implicit none
!
  integer i,j
  logical solvidx(nmol)
  integer anidx(nmol)
!
! allocate and transfer (note deallocation of molangr/solutes is handled in allocate_molecule())
  allocate(molangr(nangrps,2))
  molangr(:,1) = 0
  molangr(:,2) = 0
  do i=1,nangrps
    do j=1,nmol
      if ((anidx(j).eq.i).AND.(molangr(i,1).eq.0)) then
        molangr(i,1) = j
        molangr(i,2) = molangr(i,2) + 1
      else if (anidx(j).eq.i) then
        molangr(i,2) = molangr(i,2) + 1
      end if
    end do
  end do
  nsolutes = 0
  do j=1,nmol
    if (solvidx(j).EQV..false.) then
      nsolutes = nsolutes + 1
    end if
  end do
  if (nsolutes.gt.0) allocate(solutes(nsolutes))
  i = 0
  nressolute = 0
  do j=1,nmol
    if (solvidx(j).EQV..false.) then
      i = i + 1
      solutes(i) = j
      nressolute = nressolute + rsmol(j,2) - rsmol(j,1) + 1
    end if
  end do
!
end
!
!------------------------------------------------------------------------------------
!
subroutine read_particleflucfile()
!
  use grandensembles
  use molecule
  use iounit
  use system
  use units
  use math
  use movesets
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer cardinality, getmember, freeunit
  logical ismember
  integer t1,t2,i,imol,linecount,typid,initialnumber,iunit
  integer iomess, pos, aone
  logical exists
  character(MAXSTRLEN) line
!  
 79   format(FORM_MAXSTRLEN)
!
  aone = 1
!
! Makes sure the file exists and opens it
  call strlims(particleflucfile,t1,t2)
  inquire(file=particleflucfile(t1:t2), exist=exists)
  if (.not. exists) then
    write(ilog,*) 'FATAL: Cannot open specified input file (',&
 &particleflucfile(t1:t2), ') for FMCSC_PARTICLEFLUCFILE.'
    call fexit()
  end if
  iunit = freeunit()
  open (unit=iunit,file=particleflucfile(t1:t2),status='old')
!  
! Initialize data structures
  call allocate_particlefluc(aone)
  do typid = 1, nmoltyp
    thermalvolume(typid) = (planckconstant * sqrt(invtemp/2/PI/molma&
 &ss(typid))) ** 3
  end do
  do imol = 1, nmol
    call insertmol(imol)
  end do
!  
!  write(ilog,*) 'Parsing particle fluctuation file:'
! Reads the number of lines that specify fluctuating types
  read(iunit,79,iostat=iomess) line
  pos = 1
  linecount = -1
  call extract_int(line, linecount, pos)
  if (linecount.lt.1) then
    write(ilog,*) 'FATAL: In GCMC-input, an invalid number was speci&
 &fied for the number of lines describing fluctuating types.'
    call fexit()
  end if
  if ((ens%flag.eq.6).AND.(linecount.lt.2)) then
    write(ilog,*) 'FATAL: In SGCMC-input fewer than two fluctuating &
 &types were specified for a semigrand ensemble simulation.'
    call fexit()
  end if
!  
! Reads the fluctuating moltyps, initial number of particles, and
! their chemical potentials:
  fluclst%nr = 0
  do i = 1, linecount
    read(iunit, 79, iostat=iomess) line
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Fatal. Input file for (S)GCMC ensemble does not contain &
 &as many lines as indicated on the first line. Make sure there are as many rows &
 &as fluctuating molecule types after the header line (',particleflucfile(t1:t2),').'
      close(unit=iunit)
      call fexit()
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing (S)GCMC&
 & input (got: ',particleflucfile(t1:t2),').'
      call fexit()
    end if
!
    pos = 1
    typid = -1
    call extract_int(line, typid, pos)
!    call extract_int(line, initialnumber, pos)
    read (line(pos:),*) eqnum(typid),chempot(typid)
    initialnumber = nint(eqnum(typid))
    if ((typid.lt.1).OR.(typid.gt.nmoltyp)) then
      write(ilog,*) 'FATAL: In (S)GCMC-input, an invalid number was &
 &specified for the TypeID.'
      call fexit()
    end if
    if (ismember(fluctypes, typid).EQV..true.) then
      write(ilog,*) 'FATAL: In (S)GCMC-input, a duplicate TypeID was&
 & specified.'
      call fexit()
    end if
    if (initialnumber.gt.cardinality(typesets(typid,1))) then
      write(ilog,*) 'FATAL: In (S)GCMC-input, the initial number of &
 &particles exceeds the number of reserve particles available for th&
 &at type.'
      call fexit()
    end if
!    write(ilog,*) 'TypeID: ', typid, ' Initial number: ',&
! &initialnumber, ' mu: ', chempot(typid), 'kcal/mol'
    call addtoset(fluctypes, typid)
    fluclst%nr = fluclst%nr + 1 
    do while (cardinality(typesets(typid,1)).gt.initialnumber)
      call deletemol(getmember(typesets(typid,1), 1))
    end do
  end do
!
  allocate(fluclst%idx(fluclst%nr))
  do i=1,fluclst%nr
    fluclst%idx(i) = fluctypes%array(i)
  end do
! 
  close(unit=iunit)
!
 77   format('Mol. type #:',i4,' (example molecule is ',i7,') has a maxi&
 &mum number of ',i7,' molecules.')
 76   format('Its starting number is ',i7,' molecules and its total chemical p&
 &otential is ',g12.5,' kcal/mol.')
 75   format('Its expected average as well as starting number is ',i7,' molecules and its excess chemical p&
 &otential is ',g12.5,' kcal/mol.')
  if (sgc_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '- Summary of Molecule Types for (S)GC-Ensemble -'
    do i=1,fluclst%nr
      write(ilog,77) fluctypes%array(i),&
 &   moltyp(fluctypes%array(i),1),moltyp(fluctypes%array(i),2)
      if (gc_mode.eq.1) then
        write(ilog,76) cardinality(typesets(fluctypes%array(i),1)),&
 &chempot(fluctypes%array(i))
      else if (gc_mode.eq.2) then
        write(ilog,75) cardinality(typesets(fluctypes%array(i),1)),&
 &chempot(fluctypes%array(i))
      else
        write(ilog,*) 'Fatal. Encountered unknown mode for GC implementation in &
 &read_particleflucfile(). This is a bug.'
        call fexit()
      end if
    end do
    write(ilog,*) '-       End of Summary for (S)GC-Ensemble      -'
    write(ilog,*) 
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_refile()
!
  use iounit
  use mpistuff
  use energies
  use system
  use units
  use dssps
  use ems
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer iunit,freeunit,iomess,t1,t2,i,k
  logical exists,chgd,reqimps(4)
!
  reqimps(:) = .false.
!
  call strlims(re_infile,t1,t2)
  inquire(file=re_infile(t1:t2),exist=exists)
  if (exists.EQV..true.) then
    iunit = freeunit()
    open (unit=iunit,file=re_infile(t1:t2),status='old')
  else
    write(ilog,*) 'FATAL: Cannot open default input file (', &
 &re_infile(t1:t2),') for replica exchange spec.s.'
    call fexit()
  end if
!
  read(iunit,*,iostat=iomess) (re_types(k),k=1,re_conddim)
  if (iomess.eq.IOSTAT_END) then
    write(ilog,*) 'Fatal. Not enough dimensions in first line of RE &
 &input file (got: ',re_infile(t1:t2),').'
    call fexit()
  else if (iomess.eq.2) then
    write(ilog,*) 'Fatal. File I/O error while processing replica ex&
 &change input (got: ',re_infile(t1:t2),').'
    call fexit()
  end if
!
  do i=1,re_conddim
    do k=i+1,re_conddim
      if (re_types(i).eq.re_types(k)) then
        write(ilog,*) 'Fatal. Specified condition twice in file ',&
 &re_infile(t1:t2),' for RE input.'
        write(ilog,*) 'Positions ',i,' and ',k,' with ',re_types(i),&
 &'!'
        close(unit=iunit)
        call fexit()
      end if
    end do
  end do
!
  do i=1,re_conditions
    read(iunit,*,iostat=iomess) (re_mat(i,k),k=1,re_conddim)
!    write(*,*) (re_mat(i,k),k=1,re_conddim)
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Fatal. Not enough dimensions or conditions in R&
 &E input file (line ',i+1,' in ',re_infile(t1:t2),').'
      close(unit=iunit)
      call fexit()
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while processing replica &
 &exchange input (got: ',re_infile(t1:t2),').'
      close(unit=iunit)
      call fexit()
    end if
  end do
  close(unit=iunit)
!
! now parse the matrix and set the appropriate parameters
! note that we won't tolerate ANY out-of-range parameters to avoid
! confusion
! alse note that in RE scale_XX=0.0 does not have to equal use_XX = .false.
! the reason is that in order to compute foreign energies some arrays need to
! be allocated and some parameters set, which depend on the use_XX-flags  
  do k=1,re_conddim
    if (re_types(k).eq.1) then
      kelvin = re_mat(myrank+1,k)
      if (kelvin.le.0.0d0) then
        write(ilog,*) 'Fatal. Requested illegal value for temperatur&
 &e in REMC.'
        call fexit()
      end if
      invtemp = 1.0/(gasconst*kelvin)
    else if (re_types(k).eq.2) then
      scale_IPP = re_mat(myrank+1,k)
      use_IPP = .true.
      if (scale_IPP.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for IPP scalin&
 &g factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.3) then
      scale_attLJ = re_mat(myrank+1,k)
      use_attLJ = .true.
      if (scale_attLJ.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for attLJ scal&
 &ing factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.4) then
      scale_WCA = re_mat(myrank+1,k)
      use_WCA = .true.
      if (scale_WCA.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for WCA scalin&
 &g factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.5) then
      scale_POLAR = re_mat(myrank+1,k)
      use_POLAR = .true.
      if (scale_POLAR.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for POLAR scal&
 &ing factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.6) then
      scale_IMPSOLV = re_mat(myrank+1,k)
      use_IMPSOLV = .true.
      if (scale_IMPSOLV.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for IMPSOLV sc&
 &aling factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.7) then
      par_IMPSOLV(2)=re_mat(myrank+1,k)
      if (par_IMPSOLV(2).lt.1.0) then
        write(ilog,*) 'Fatal. Requested illegal value for implicit s&
 &olvent dielectric in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.8) then
      scale_TOR=re_mat(myrank+1,k)
      use_TOR = .true.
      if (scale_TOR.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for TOR scalin&
 &g factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.9) then
      scale_ZSEC=re_mat(myrank+1,k)
      use_ZSEC = .true.
      if (scale_ZSEC.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for ZSEC scali&
 &ng factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.10) then
      par_ZSEC(1)=re_mat(myrank+1,k)
      if ((par_ZSEC(1).lt.0.0).OR.&
 &      (par_ZSEC(1).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for target alp&
 &ha-content in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.11) then
      par_ZSEC(3)=re_mat(myrank+1,k)
      if ((par_ZSEC(3).lt.0.0).OR.&
 &            (par_ZSEC(3).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for target bet&
 &a-content in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.12) then
      scrq_model=nint(re_mat(myrank+1,k))
      reqimps(4) = .true.
      if ((scrq_model.lt.1).OR.(scrq_model.gt.9)) then
        write(ilog,*) 'Fatal. Requested illegal value for screening &
 &model in REMC.'
        call fexit()
      end if
    else  if (re_types(k).eq.13) then
      par_IMPSOLV(3)=re_mat(myrank+1,k)
      if ((par_IMPSOLV(3).le.0.0).OR.&
 &               (par_IMPSOLV(3).gt.100.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for FOS-tau in&
 & REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.14) then
      par_IMPSOLV(4)=re_mat(myrank+1,k)
      if ((par_IMPSOLV(4).le.0.0).OR.&
 &            (par_IMPSOLV(4).gt.100.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for SCR-tau in&
 & REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.15) then
      par_IMPSOLV(6)=re_mat(myrank+1,k)
      if ((par_IMPSOLV(6).lt.0.0).OR.&
 &            (par_IMPSOLV(6).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for FOS-midpoi&
 &nt in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.16) then
      par_IMPSOLV(7)=re_mat(myrank+1,k)
      if ((par_IMPSOLV(7).lt.0.0).OR.&
 &            (par_IMPSOLV(7).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for SCR-midpoi&
 &nt in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.17) then
      par_IMPSOLV(8)=re_mat(myrank+1,k)
      if (par_IMPSOLV(8).lt.1.0) then
        write(ilog,*) 'Fatal. Requested illegal value for contact di&
 &electric in REMC.'
        call fexit()
      end if
      par_IMPSOLV(8) = 1./par_IMPSOLV(8)
    else if (re_types(k).eq.18) then
      i_sqm=nint(re_mat(myrank+1,k))
      if ((i_sqm.lt.-10).OR.(i_sqm.gt.10)) then
        write(ilog,*) 'Fatal. Requested illegal value for generalize&
 &d mean in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.19) then
      par_IMPSOLV(9)=re_mat(myrank+1,k)
      if ((par_IMPSOLV(9).lt.0.0).OR.&
 &              (par_IMPSOLV(9).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for SCRMIX in &
 &REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.20) then
      use_FEGS(1) = .true.
      scale_FEGS(1) = re_mat(myrank+1,k)
      if (scale_FEGS(1).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_IPP in&
 & REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.21) then
      use_FEGS(3) = .true.
      scale_FEGS(3) = re_mat(myrank+1,k)
      if (scale_FEGS(3).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_attLJ &
 &in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.22) then
      use_FEGS(6) = .true.
      scale_FEGS(6) = re_mat(myrank+1,k)
      if (scale_FEGS(6).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_POLAR &
 &in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.23) then
      scale_TABUL = re_mat(myrank+1,k)
      use_TABUL = .true.
      if (scale_TABUL.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for TABUL scal&
 &ing factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.24) then
      scale_POLY = re_mat(myrank+1,k)
      use_POLY = .true.
      if (scale_POLY.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for POLY scali&
 &ng factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.25) then
      scale_DREST = re_mat(myrank+1,k)
      use_DREST = .true.
      if (scale_DREST.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for DREST scal&
 &ing factor in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.26) then
      use_FEGS(15) = .true.
      scale_FEGS(15) = re_mat(myrank+1,k)
      if (scale_FEGS(15).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_BONDED&
 &_B in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.27) then
      use_FEGS(16) = .true.
      scale_FEGS(16) = re_mat(myrank+1,k)
      if (scale_FEGS(16).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_BONDED&
 &_A in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.28) then
      use_FEGS(17) = .true.
      scale_FEGS(17) = re_mat(myrank+1,k)
      if (scale_FEGS(17).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_BONDED&
 &_I in REMC.'
        call fexit()
      end if
    else if (re_types(k).eq.29) then
      use_FEGS(18) = .true.
      scale_FEGS(18) = re_mat(myrank+1,k)
      if (scale_FEGS(18).lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for FEG_BONDED&
 &_T in REMC.'
        call fexit()
      end if
    else if ((re_types(k).eq.30)) then
      par_DSSP(9) = re_mat(myrank+1,k)
      if ((par_DSSP(9).lt.0.0d0).OR.(par_DSSP(9).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for FMCSC_DSSP_ESC.'
        call fexit()
      end if
      write(ilog,*) 'Warning. Overlap/TI calculation with DSSP restraint windows&
  & is not yet implemented.  Please write the necessary der_dssp_gl function.'
    else if ((re_types(k).eq.31)) then
      par_DSSP(7) = re_mat(myrank+1,k)
      if ((par_DSSP(7).lt.0.0d0).OR.(par_DSSP(7).gt.1.0)) then
        write(ilog,*) 'Fatal. Requested illegal value for FMCSC_DSSP_HSC.'
        call fexit()
      end if
      write(ilog,*) 'Warning. Overlap/TI calculation with DSSP restraint windows&
  & is not yet implemented.  Please write the necessary der_dssp_gl function.'
!   for MC this is a dummy dimension
    else if (re_types(k).eq.32) then
      dyn_dt = re_mat(myrank+1,k)
      if (dyn_dt.le.0.0d0) then
        dyn_dt = 0.002 ! 2fs
        write(ilog,*) 'Fatal. Requested illegal value for FMCSC_TIMESTEP in REMC.'
        call fexit()
      end if
    else  if (re_types(k).eq.33) then
      scale_EMICRO = re_mat(myrank+1,k)
      use_EMICRO = .true.
      if (scale_EMICRO.lt.0.0) then
        write(ilog,*) 'Fatal. Requested illegal value for EMICRO scaling factor in&
 & REMC.'
        call fexit()
      end if
#ifdef LINK_NETCDF
!     do nothing further
#else
      write(ilog,*) 'Fatal. Spatial density restraints (EMICRO) are available only if &
 &CAMPARI is linked against NetCDF (pass LINK_NETCDF to the compiler).'
      call fexit
#endif
    else  if (re_types(k).eq.34) then
      emthreshdensity = re_mat(myrank+1,k)
#ifdef LINK_NETCDF
!     do nothing further
#else
      write(ilog,*) 'Fatal. Spatial density restraints (EMICRO) are available only if &
 &CAMPARI is linked against NetCDF (pass LINK_NETCDF to the compiler).'
      call fexit
#endif
    end if
  end do
!
! some sanity checks for RE-dimensions depending on other
! parameters, which might or might not be set properly
  if (reqimps(4).EQV..false.) then
    if ((scrq_model.ge.5).AND.(scrq_model.le.8)) reqimps(1) = .true.
    if ((scrq_model.eq.3).OR.((scrq_model.ge.7).AND.(scrq_model.le.9))) reqimps(2:3) = .true.
    if (scrq_model.eq.4) reqimps(2) = .true.
  else
    do k=1,re_conddim
      if (re_types(k).eq.12) then
        do i=1,re_conditions
          if ((re_mat(i,k).ge.5).AND.(re_mat(i,k).le.8)) reqimps(1) = .true.
          if ((re_mat(i,k).eq.3).OR.((re_mat(i,k).ge.7).AND.(re_mat(i,k).le.9))) reqimps(2:3) = .true.
          if (re_mat(i,k).eq.4) reqimps(2) = .true.
        end do
      end if
    end do
  end if
  do k=1,re_conddim
    if (((re_types(k).eq.12).OR.&
 &     ((re_types(k).ge.13).AND.(re_types(k).le.19)))&
 &        .AND.(use_IMPSOLV.EQV..false.)) then
      write(ilog,*) 'Fatal. Specified implicit solvent parameter(s) &
 &as RE dimension(s), but did not turn on IMPSOLV. Check input.'
      call fexit()
    end if
    if (((re_types(k).eq.17).AND.(reqimps(2).EQV..false.)).OR.((re_types(k).eq.19).AND.(reqimps(3).EQV..false.)).OR.&
        ((re_types(k).eq.18).AND.(reqimps(1).EQV..false.))) then
      write(ilog,*) 'Fatal. Specified implicit solvent parameter(s) &
 &as RE dimension(s), but the screening model(s) in use does/do not require these parameters.'
      call fexit()
    end if
    if ((reqimps(4).EQV..false.).AND.(scrq_model.eq.4).AND.((re_types(k).eq.14).OR.(re_types(k).eq.16).OR.&
 &      (re_types(k).eq.18).OR.(re_types(k).eq.19))) then
      write(ilog,*) 'Fatal. Specified implicit solvent parameter(s) &
 &as RE dimension(s), but did not turn on IMPSOLV. Check input.'
      call fexit()
    end if
    if (((re_types(k).eq.20).OR.(re_types(k).eq.21).OR.&
 &       (re_types(k).eq.22).OR.(re_types(k).eq.26).OR.&
 &       (re_types(k).eq.27).OR.(re_types(k).eq.28).OR.&
 &       (re_types(k).eq.29)).AND.&
 &           (use_FEG.EQV..false.)) then
      write(ilog,*) 'Fatal. Specified FEG scaling factors but did no&
 &t request ghosting. Check input.'
      call fexit()
    end if
    if (((re_types(k).eq.10).OR.(re_types(k).eq.11)).AND.&
 &            (use_ZSEC.EQV..false.)) then
      write(ilog,*) 'Fatal. Specified target content for 2-structure&
 & biasing potential as RE dimension, but set scaling factor to 0?'
      call fexit()
    end if
    if (((re_types(k).eq.30).OR.(re_types(k).eq.31)).AND.&
 &            (use_DSSP.EQV..false.)) then
      write(ilog,*) 'Fatal. Specified target content for DSSP&
 & biasing potential as RE dimension, but set scaling factor to 0?'
      call fexit()
    end if
    if ((re_types(k).eq.34).AND.(use_EMICRO.EQV..false.)) then
      write(ilog,*) 'Fatal. Specified threshold density for EM &
 &restraint potential as RE dimension, but set scaling factor to 0?'
      call fexit()
    end if
  end do
!
! some more checks for invariant RE-dimensions
  do k=1,re_conddim
    chgd = .false.
    do i=2,re_conditions
      if (re_mat(i,k).ne.re_mat(i-1,k)) then
        chgd = .true.
        exit
      end if
    end do
    if ((chgd.EQV..false.).AND.(re_types(k).ne.32)) then
      write(ilog,*) 'Fatal. Dimension #',k,' (type #',re_types(k),')&
 & is invariant across replicas. Please eliminate!'
      call fexit()
    end if
  end do
!
! since we don't tolerate exceptions, we know the other replicas will be
! fine if the program is still running -> no need to check more on them
!
! now do some setup work for overlap measures
!
  Tonly = .false.
  Tthere = .false.
  noTI = .false.
  do k=1,re_conddim
    if (re_types(k).eq.1) then
      noTI = .true.
      Tthere = .true.
      if (re_conddim.eq.1) then
        Tonly = .true.
      end if
      Tdim = k
    else if ((re_types(k).eq.7).OR.&
 &  ((re_types(k).ge.12).AND.(re_types(k).le.22)).OR.&
 &  ((re_types(k).ge.30).AND.(re_types(k).le.31)).OR.&
 &  ((re_types(k).eq.34))) then
      noTI = .true.
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine read_snapsetfile(ivx,ivxsz,which)
!
  use clusters
  use iounit
  use interfaces
!
  implicit none
!
  integer, INTENT(IN):: ivxsz,ivx(ivxsz,3),which
!
  integer iunit,i,j,k,freeunit,nfold,counter,tmp,t1,t2,jj
  integer, ALLOCATABLE:: tmpvec(:),tmpvec_srt(:)
  logical exists,atrue,afalse
  character(MAXSTRLEN) fn
!
  atrue = .true.
  afalse = .false.
!
  fn(:) = ' '
  if (which.eq.1) then
    fn = clufoldfile
  else if (which.eq.2) then
    fn = cluunfoldfile
  else
    write(ilog,*) 'Fatal. Called read_snapsetfile(...) with unsupported mode (got ',which,'). This is a bug'
    call fexit()
  end if
  call strlims(fn,t1,t2)

  inquire(file=fn(t1:t2),exist=exists)
  if (exists.EQV..false.) then  !clufold does not get allocated and this is important in graph_algorithms
    if (which.eq.1) then
      write(ilog,*) 'Warning. ',fn(t1:t2),' not found (FMCSC_CLUFOLDFILE). Input from &
 &FMCSC_INISYNSNAP or centroid of largest cluster will be used instead.'
      return
    else if (which.eq.2) then
      write(ilog,*) 'Warning. ',fn(t1:t2),' not found (FMCSC_CLUUNFOLDFILE). Input from &
 &FMCSC_ENDSYNSNAP or cluster containing last stored snapshot will be used instead.'
      return
    end if
  end if
  iunit = freeunit()
  open(unit=iunit,file=fn(t1:t2),status='old')
!
  allocate(tmpvec(cstored))
  allocate(tmpvec_srt(cstored))
!
  nfold = 0 ! transparency
  call read_atmidx(iunit,nfold,tmpvec,cstored,afalse,afalse)
  close(unit=iunit) 
!
! we try to find (exactly) the read snapshots numbers in the actual stored data
  k = cstored
  counter = 0
  tmp = maxval(ivx(:,3))
  do i=1,nfold
    if ((tmpvec(i).le.0).OR.(tmpvec(i).gt.tmp)) then
      write(ilog,*) 'Warning. Found a negative or too large number in CLUFOLDFILE at position: ',i,'. Ignoring.'
      tmpvec(i) = 0
      cycle
    end if
    counter = counter + 1
    j = tmpvec(i)
    call binary_search(k,ivx(1:k,3),j,jj)
    if (ivx(max(1,min(k,jj)),3).ne.tmpvec(i)) then
      write(ilog,*) 'Fatal. Unable to find snapshot ',tmpvec(i),' in the data from the original having actually been stored &
 &(FRAMESFILE, CCOLLECT, ...) in read_snapsetfile(...).'
      call fexit()
    end if
  end do
!
  if ((nfold.le.0).OR.(counter.le.0)) then
    deallocate(tmpvec)
    if (which.eq.1) then
      write(ilog,*) 'Warning. ',fn(t1:t2),' (FMCSC_CLUFOLDFILE) is empty or contains no interpretable numbers. Input from &
 &FMCSC_INISYNSNAP or centroid of largest cluster will be used instead.'
      return
    else if (which.eq.2) then
      write(ilog,*) 'Warning. ',fn(t1:t2),' (FMCSC_CLUUNFOLDFILE) is empty or contains no interpretable numbers. Input from &
 &FMCSC_ENDSYNSNAP or cluster containing last stored snapshot will be used instead.'
      return
    end if
  end if
!
! sort the temporary array
  k = 1
  call merge_sort(ldim=nfold,up=atrue,list=tmpvec,olist=tmpvec_srt,ilo=k,ihi=nfold)
  deallocate(tmpvec)
!
! remove duplicates and zeroes
  if (nfold.ge.2) then
    exists = .true.
    i = 2
    do while (exists.EQV..true.)
      if (tmpvec_srt(i).eq.tmpvec_srt(i-1)) then
        write(ilog,*) 'Warning. Removing duplicated snapshot from snapshot list in ',fn(t1:t2),'.'
        tmpvec_srt((i-1):(nfold-1)) = tmpvec_srt(i:nfold)
        nfold = nfold - 1
      else if (tmpvec_srt(i-1).eq.0) then
        tmpvec_srt((i-1):(nfold-1)) = tmpvec_srt(i:nfold)
        nfold = nfold - 1
      else 
        i = i + 1
      end if
      if (i.gt.nfold) exists = .false.
    end do
  end if
!
! at this point, we still have snapshot lists, which later need to be translated to cluster lists (additional duplicates poss.!)
  if (which.eq.1) then
    allocate(clufold(nfold))
    clufold(1:nfold) = tmpvec_srt(1:nfold)
  else if (which.eq.2) then
    allocate(cluunfold(nfold))
    cluunfold(1:nfold) = tmpvec_srt(1:nfold)
  end if
  deallocate(tmpvec_srt)
!
end
!
!---------------------------------------------------------------------------------------------------------------------------------
!
