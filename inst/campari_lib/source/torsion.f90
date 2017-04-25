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
! #################################################################
! ##                                                             ##
! ## subroutine torsions -- periodic dumpout of peptide torsions ##
! ##                                                             ##
! #################################################################
!
! "torsions" dumps out the backbone and sidechain torsions
!
!
subroutine torsion()
!
  use mcsums
  use mpistuff
  use grandensembles
  use system
  use molecule
  use zmatrix
  use atoms
  use torsn
  use sequen
  use polypep
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer j,ttc,imol,shf,rs
  RTYPE offset,getztor
  logical ismember
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
 10   format(i12,1x,20000000(f9.3,1x))
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    imol = 0
    do ttc=1,nzmlst
      if (molofrs(atmres(abs(torzmlst(ttc)))).ne.imol) then
        imol = molofrs(atmres(abs(torzmlst(ttc))))
        offset = 0.0
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
            if (ismember(ispresent,imol).EQV..false.) offset = 720.0
          end if
        end if
      end if
      if (torzmlst(ttc).gt.0) then
        curtvec(ttc) = ztor(torzmlst(ttc)) + offset
      else
        rs = atmres(abs(torzmlst(ttc)))
        if (disulf(rs).gt.0) then
          if (crosslink(crlk_idx(rs))%itstype.le.2) then
            if (abs(torzmlst(ttc)).eq.at(rs)%sc(2-shf)) then
              curtvec(ttc) = getztor(at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf),at(disulf(rs))%sc(2-shf))
            else if (abs(torzmlst(ttc)).eq.at(rs)%sc(3-shf)) then
              curtvec(ttc) = getztor(cai(disulf(rs)),at(disulf(rs))%sc(2-shf),at(disulf(rs))%sc(3-shf),at(rs)%sc(3-shf))
            else
              curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
            end if
          else
            curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
          end if
        else
          curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
        end if
      end if
    end do
!
    write(idihed,10) nstep,(curtvec(j),j=1,nzmlst)
!
  else if (use_MPIAVG.EQV..true.) then
    if (myrank.eq.0) then
      call MPI_AVGWriteTor()
    else
      call MPI_AVGSendTor()
    end if
    if (mpi_cnt_tor.eq.(mpi_nodes-1)) then
      mpi_cnt_tor = 0
    else
      mpi_cnt_tor = mpi_cnt_tor + 1
    end if

  end if
#else
!
  imol = 0
  do ttc=1,nzmlst
    if (molofrs(atmres(abs(torzmlst(ttc)))).ne.imol) then
      imol = molofrs(atmres(abs(torzmlst(ttc))))
      offset = 0.0
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
          if (ismember(ispresent,imol).EQV..false.) offset = 720.0
        end if
      end if
    end if
    if (torzmlst(ttc).gt.0) then
      curtvec(ttc) = ztor(torzmlst(ttc)) + offset
    else
      rs = atmres(abs(torzmlst(ttc)))
      if (disulf(rs).gt.0) then
        if (crosslink(crlk_idx(rs))%itstype.le.2) then
          if (abs(torzmlst(ttc)).eq.at(rs)%sc(2-shf)) then
            curtvec(ttc) = getztor(at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf),at(disulf(rs))%sc(2-shf))
          else if (abs(torzmlst(ttc)).eq.at(rs)%sc(3-shf)) then
            curtvec(ttc) = getztor(cai(disulf(rs)),at(disulf(rs))%sc(2-shf),at(disulf(rs))%sc(3-shf),at(rs)%sc(3-shf))
          else
            curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
          end if
        else
          curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
        end if
      else
        curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
      end if
    end if
  end do
!
  write(idihed,10) nstep,(curtvec(j),j=1,nzmlst)
!
#endif
!
end
!
!----------------------------------------------------------------------
!
subroutine torsion_header()
!
  use fyoc
  use sequen
  use movesets
  use mcsums
  use aminos
  use torsn
  use mpistuff
  use iounit
  use polypep
  use system
  use atoms
!
  implicit none
!
  integer i,j,tl,ttc,shf,shfbu
  character(2) rsnma
  character(3) resname
  character(10), ALLOCATABLE:: stringy(:)
  integer, ALLOCATABLE:: tmplst(:)
  logical ressbeg
!
  tl = 2
  ttc = 0
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  call torsion_assemble()
!
  if (toroutmode.eq.1) then ! the traditional version
    allocate(stringy(size(curtvec)))
    if (allocated(torzmlst).EQV..true.) deallocate(torzmlst)
    allocate(tmplst(n))
    nzmlst = 0
    do i=1,size(curtvec)
      stringy(i) = '          '
    end do
!
    do i=1,nseq
      ressbeg = .false.
      resname = amino(seqtyp(i))
      if (wline(i).gt.0) then
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = wline(i)
        if (ressbeg.EQV..false.) then
          stringy(nzmlst) = '| '//resname//' OME '
          ressbeg = .true.
        else
          stringy(nzmlst) = '   OME    '
        end if
      end if
      if ((fline(i).gt.0).AND.(seqflag(i).ne.5)) then
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = fline(i)
        if (ressbeg.EQV..false.) then
          stringy(nzmlst) = '| '//resname//' PHI '
          ressbeg = .true.
        else
          stringy(nzmlst) = '   PHI    '
        end if
      end if
      if (yline2(i).gt.0) then
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = yline2(i)
        if (ressbeg.EQV..false.) then
          stringy(nzmlst) = '| '//resname//' PSI '
          ressbeg = .true.
        else
          stringy(nzmlst) = '   PSI    '
        end if
      else if (yline(i).gt.0) then
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = yline(i)
        if (ressbeg.EQV..false.) then
          stringy(nzmlst) = '| '//resname//' PSI '
          ressbeg = .true.
        else
          stringy(nzmlst) = '   PSI    '
        end if
      end if
      do j=1,nnucs(i)
        call int2str(j,rsnma,tl)
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = nucsline(j,i)
        if (ressbeg.EQV..false.) then
          stringy(nzmlst) = '|'//resname//' NUC'//rsnma(1:tl)
          ressbeg = .true.
        else
          stringy(nzmlst) = '  NUC '//rsnma(1:tl)//'  '
        end if
      end do
      do j=1,nchi(i)
        call int2str(j,rsnma,tl)
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = chiline(j,i)
        if (ressbeg.EQV..false.) then
          stringy(nzmlst) = '|'//resname//' CHI'//rsnma(1:tl)
          ressbeg = .true.
        else
          stringy(nzmlst) = '  CHI '//rsnma(1:tl)//'  '
        end if
      end do
      if (disulf(i).gt.0) then
        call int2str(j,rsnma,tl)
        if (crosslink(crlk_idx(i))%itstype.le.2) then ! disulfide
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = -at(disulf(i))%sc(3-shf)
          if (ressbeg.EQV..false.) then
            stringy(nzmlst) = '|'//resname//' CHI02'
            ressbeg = .true.
          else
            stringy(nzmlst) = '   CHI 02 '
          end if
          if (disulf(i).gt.i) then
            nzmlst = nzmlst + 1
            tmplst(nzmlst) = -at(disulf(i))%sc(2-shf)
            if (ressbeg.EQV..false.) then
              stringy(nzmlst) = '|'//resname//' CSSC '
              ressbeg = .true.
            else
              stringy(nzmlst) = '   CSSC   '
            end if
          end if
        end if
      end if
      if ((pucline(i).gt.0).OR.((seqpolty(i).eq.'N').AND.(nucsline(6,i).gt.0))) then
        do j=1,7
          call int2str(j,rsnma,tl)
          nzmlst = nzmlst + 1
          if (ressbeg.EQV..false.) then
            stringy(nzmlst) = '|'//resname//' PUC'//rsnma(1:tl)
            ressbeg = .true.
          else
            stringy(nzmlst) = '  PUC '//rsnma(1:tl)//'  '
          end if
        end do
        if (seqpolty(i).eq.'N') then
          shfbu = shf
          shf = 1
        end if
        if (seqpolty(i).eq.'N') then
          tmplst(nzmlst-6) = nucsline(6,i)
        else
          tmplst(nzmlst-6) = fline(i)
        end if
        tmplst(nzmlst-5) = at(i)%sc(2-shf)
        tmplst(nzmlst-4) = at(i)%sc(3-shf)
        tmplst(nzmlst-3) = at(i)%sc(4-shf)
        tmplst(nzmlst-2) = -at(i)%sc(2-shf)
        tmplst(nzmlst-1) = -at(i)%sc(3-shf)
        tmplst(nzmlst) = -at(i)%sc(4-shf)
        if (seqpolty(i).eq.'N') shf = shfbu
      end if
    end do
    if (nzmlst.gt.0) then
      allocate(torzmlst(nzmlst))
      torzmlst(:) = tmplst(1:nzmlst)
    end if
    deallocate(tmplst)
  end if
!
  if (nzmlst.le.0) then
    write(ilog,*) 'Warning. No eligible dihedral angles to print to FYC.dat. Disabling &
 &instantaneous output of torsions (FMCSC_TOROUT).'
    torout = nsim+1
    return
  end if
!
!  write(idihed,*) '# Step Number      ',stringy(1:ttc*8) 
 34 format('#_Step_No.__|',20000000(a10))
 35 format('#_Step_No.__|',20000000(i9,'|'))
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    if (toroutmode.eq.1) then
      write(idihed,34) (stringy(j),j=1,nzmlst)
    else if (toroutmode.eq.2) then
      write(idihed,35) torzmlst(1:nzmlst)
    else 
      write(ilog,*) 'Fatal. Called torsion_header(...) with unknown output mode. This is most likely an omission bug.'
      call fexit()
    end if
  else if (use_MPIAVG.EQV..true.) then
    if (myrank.eq.0) then
      if (toroutmode.eq.1) then
        write(idihed,34) (stringy(j),j=1,nzmlst)
      else if (toroutmode.eq.2) then
        write(idihed,35) torzmlst(1:nzmlst)
      else 
        write(ilog,*) 'Fatal. Called torsion_header(...) with unknown output mode. This is most likely an omission bug.'
        call fexit()
      end if
    end if
  end if
#else
  if (toroutmode.eq.1) then
    write(idihed,34) (stringy(j),j=1,nzmlst)
  else if (toroutmode.eq.2) then
    write(idihed,35) torzmlst(1:nzmlst)
  else 
    write(ilog,*) 'Fatal. Called torsion_header(...) with unknown output mode. This is most likely an omission bug.'
    call fexit()
  end if
#endif
!
  if (allocated(stringy).EQV..true.) deallocate(stringy)
!
end
!
!---------------------------------------------------------------------
!
! a routine to assemble lists of and quasi flexible Z-matrix dihedral angles written to FYC.dat
! excludes all(!) bond angles
!
subroutine torsion_assemble()
!
  use atoms
  use molecule
  use system
  use movesets
  use forces
  use zmatrix
  use fyoc
  use torsn
  use interfaces
  use polypep
  use sequen
  use ujglobals
!
  implicit none
!
  integer i,k,ttc,imol,rs,ii,jj,shf,hangover
  integer, ALLOCATABLE:: tmplst(:),iv(:,:)
  logical atrue
  logical, ALLOCATABLE:: isin(:)
!
  nzmlst = 0
  if (toroutmode.ne.2) return
  allocate(tmplst(n))
  allocate(isin(n))
  isin(:) = .false.
  atrue = .true.
  shf = 0
  if (ua_model.gt.0) shf = 1
!
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    if ((have_omega.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,wlst%nr
        if (((wlst%wt(i).gt.wlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(wlst%wt(i).gt.0.0))).AND.&
 &          (isin(wline(wlst%idx(i))).EQV..false.)) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = wline(wlst%idx(i))
          isin(wline(wlst%idx(i))) = .true.
        end if          
      end do
    end if
    if ((have_pivot.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,fylst%nr
        if ((fylst%wt(i).gt.fylst%wt(max(1,i-1))).OR.((i.eq.1).AND.(fylst%wt(i).gt.0.0))) then
          if ((yline2(fylst%idx(i)).gt.0).AND.(yline(fylst%idx(i)).gt.0)) then
            if ((isin(yline(fylst%idx(i))).EQV..false.).AND.(izrot(yline(fylst%idx(i)))%alsz.gt.0)) then
              nzmlst = nzmlst + 1
              tmplst(nzmlst) = yline2(fylst%idx(i))
              isin(yline(fylst%idx(i))) = .true.
            end if
          else if (yline(fylst%idx(i)).gt.0) then
            if ((isin(yline(fylst%idx(i))).EQV..false.).AND.(izrot(yline(fylst%idx(i)))%alsz.gt.0)) then
              nzmlst = nzmlst + 1
              tmplst(nzmlst) = yline(fylst%idx(i))
              isin(yline(fylst%idx(i))) = .true.
            end if
          end if
          if (fline(fylst%idx(i)).gt.0) then
            if ((isin(fline(fylst%idx(i))).EQV..false.).AND.(izrot(fline(fylst%idx(i)))%alsz.gt.0)) then
              nzmlst = nzmlst + 1
              tmplst(nzmlst) = fline(fylst%idx(i))
              isin(fline(fylst%idx(i))) = .true.
            end if
          end if
        end if          
      end do
    end if
    if (((have_nuc.EQV..true.).AND.(nucpuckfreq.lt.1.0).AND.(nuccrfreq.lt.1.0)).OR.(pdb_analyze.EQV..true.)) then
      do i=1,nuclst%nr
        if ((nuclst%wt(i).gt.nuclst%wt(max(1,i-1))).OR.((i.eq.1).AND.(nuclst%wt(i).gt.0.0))) then
          do k=1,nnucs(nuclst%idx(i))
            if (isin(nucsline(k,nuclst%idx(i))).EQV..false.) then
              nzmlst = nzmlst + 1
              tmplst(nzmlst) = nucsline(k,nuclst%idx(i))
              isin(nucsline(k,nuclst%idx(i))) = .true.
            end if
          end do
        end if          
      end do
    end if
    if ((have_nucpuck.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,nucpuclst%nr
        if (((nucpuclst%wt(i).gt.nucpuclst%wt(max(1,i-1))).OR.((i.eq.1).AND.(nucpuclst%wt(i).gt.0.0))).AND.&
 &          (isin(nucsline(6,nucpuclst%idx(i))).EQV..false.)) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = nucsline(6,nucpuclst%idx(i))
          isin(nucsline(6,nucpuclst%idx(i))) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(nucpuclst%idx(i))%sc(1)
          isin(at(nucpuclst%idx(i))%sc(1)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(nucpuclst%idx(i))%sc(2)
          isin(at(nucpuclst%idx(i))%sc(2)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(nucpuclst%idx(i))%sc(3)
          isin(at(nucpuclst%idx(i))%sc(3)) = .true.
        end if          
      end do
    end if
    if ((have_pucker.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,puclst%nr
        if (((puclst%wt(i).gt.puclst%wt(max(1,i-1))).OR.((i.eq.1).AND.(puclst%wt(i).gt.0.0))).AND.&
 &          (isin(fline(puclst%idx(i))).EQV..false.)) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = fline(puclst%idx(i))
          isin(fline(puclst%idx(i))) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(puclst%idx(i))%sc(2-shf)
          isin(at(puclst%idx(i))%sc(2-shf)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(puclst%idx(i))%sc(3-shf)
          isin(at(puclst%idx(i))%sc(3-shf)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(puclst%idx(i))%sc(4-shf)
          isin(at(puclst%idx(i))%sc(4-shf)) = .true.
        end if          
      end do
    end if
    if ((have_chi.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,chilst%nr
        if ((chilst%wt(i).gt.chilst%wt(max(1,i-1))).OR.((i.eq.1).AND.(chilst%wt(i).gt.0.0))) then
          do k=1,nchi(chilst%idx(i))
            if (isin(chiline(k,chilst%idx(i))).EQV..false.) then
              nzmlst = nzmlst + 1
              tmplst(nzmlst) = chiline(k,chilst%idx(i))
              isin(chiline(k,chilst%idx(i))) = .true.
            end if
          end do
        end if          
      end do
    end if
    if ((have_unsother.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,unslst%nr
        if ((unslst%wt(i).gt.unslst%wt(max(1,i-1))).OR.((i.eq.1).AND.(unslst%wt(i).gt.0.0))) then
          if (isin(unslst%idx(i)).EQV..false.) then
            nzmlst = nzmlst + 1
            tmplst(nzmlst) = unslst%idx(i)
            isin(unslst%idx(i)) = .true.
          end if
        end if          
      end do
    end if
    if ((have_unkother.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,unklst%nr
        if ((unklst%wt(i).gt.unklst%wt(max(1,i-1))).OR.((i.eq.1).AND.(unklst%wt(i).gt.0.0))) then
          if (isin(unklst%idx(i)).EQV..false.) then
            nzmlst = nzmlst + 1
            tmplst(nzmlst) = unklst%idx(i)
            isin(unklst%idx(i)) = .true.
          end if
        end if          
      end do
    end if
    if ((have_natother.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,natlst%nr
        if ((natlst%wt(i).gt.natlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(natlst%wt(i).gt.0.0))) then
          if (isin(natlst%idx(i)).EQV..false.) then
            nzmlst = nzmlst + 1
            tmplst(nzmlst) = natlst%idx(i)
            isin(natlst%idx(i)) = .true.
          end if
        end if          
      end do
    end if
    if ((have_nuccr.EQV..true.).OR.(pdb_analyze.EQV..true.)) then
      do i=1,nuccrlst%nr
        if ((nuccrlst%wt(i).gt.nuccrlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(nuccrlst%wt(i).gt.0.0))) then
          do rs=nuccrlst%idx2(i),nuccrlst%idx(i)
            do k=1,nnucs(rs)
              if ((isin(nucsline(k,rs)).EQV..false.).AND.(izrot(nucsline(k,rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = nucsline(k,rs)
                isin(nucsline(k,nuclst%idx(i))) = .true.
              end if
            end do
          end do
        end if          
      end do
    end if
    if (have_djcr.EQV..true.) then
      do i=1,djcrlst%nr
        if ((djcrlst%wt(i).gt.djcrlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(djcrlst%wt(i).gt.0.0))) then
          do rs=djcrlst%idx2(i),djcrlst%idx(i)
            imol = molofrs(rs)
            hangover = 0
            if (rs.lt.(djcrlst%idx(i)-2-((torcrmaxsz2+1)/3))) cycle
            if (rs.eq.(djcrlst%idx(i)-2-((torcrmaxsz2+1)/3))) hangover = mod(torcrmaxsz2+2,3)
            if ((yline2(rs).gt.0).AND.(yline(rs).gt.0)) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline2(rs)
                isin(yline(rs)) = .true.
              end if
            else if (yline(rs).gt.0) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline(rs)
                isin(yline(rs)) = .true.
              end if
            end if
            if ((fline(rs).gt.0).AND.(rs.gt.rsmol(imol,1)).AND.(hangover.ne.1)) then
              if ((isin(fline(rs)).EQV..false.).AND.(izrot(fline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = fline(rs)
                isin(fline(rs)) = .true.
              end if
            end if
            if ((fline(rs).gt.0).AND.(rs.ge.(djcrlst%idx(i)-2)).AND.(seqflag(rs).eq.5)) then ! pucker
              if (isin(fline(rs)).EQV..false.) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = fline(rs)
                isin(fline(rs)) = .true.
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = at(rs)%sc(2-shf)
                isin(at(rs)%sc(2-shf)) = .true.
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = at(rs)%sc(3-shf)
                isin(at(rs)%sc(3-shf)) = .true.
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = at(rs)%sc(4-shf)
                isin(at(rs)%sc(4-shf)) = .true.
              end if
            end if
            if ((wline(rs).gt.0).AND.(rs.lt.(djcrlst%idx(i)-1)).AND.(hangover.ne.1).AND.(hangover.ne.2)) then
              if ((isin(wline(rs)).EQV..false.).AND.(izrot(wline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = wline(rs)
                isin(wline(rs)) = .true.
              end if
            end if
          end do
        end if          
      end do
    end if
    if (have_docr.EQV..true.) then
      do i=1,docrlst%nr
        if ((docrlst%wt(i).gt.docrlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(docrlst%wt(i).gt.0.0))) then
          do rs=docrlst%idx2(i),docrlst%idx(i)
            imol = molofrs(rs)
            hangover = 0
            if (rs.lt.(docrlst%idx(i)-2-((torcrmaxsz-1)/3))) cycle
            if (rs.eq.(docrlst%idx(i)-2-((torcrmaxsz-1)/3))) hangover = mod(torcrmaxsz,3)
            if ((yline2(rs).gt.0).AND.(yline(rs).gt.0)) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline2(rs)
                isin(yline(rs)) = .true.
              end if
            else if (yline(rs).gt.0) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline(rs)
                isin(yline(rs)) = .true.
              end if
            end if
            if ((fline(rs).gt.0).AND.(rs.gt.rsmol(imol,1)).AND.(hangover.ne.1)) then
              if ((isin(fline(rs)).EQV..false.).AND.(izrot(fline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = fline(rs)
                isin(fline(rs)) = .true.
              end if
            end if
            if ((fline(rs).gt.0).AND.(rs.ge.(docrlst%idx(i)-1)).AND.(seqflag(rs).eq.5)) then ! pucker
              if (isin(fline(rs)).EQV..false.) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = fline(rs)
                isin(fline(rs)) = .true.
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = at(rs)%sc(2-shf)
                isin(at(rs)%sc(2-shf)) = .true.
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = at(rs)%sc(3-shf)
                isin(at(rs)%sc(3-shf)) = .true.
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = at(rs)%sc(4-shf)
                isin(at(rs)%sc(4-shf)) = .true.
              end if
            end if
            if ((wline(rs).gt.0).AND.(hangover.ne.1).AND.(hangover.ne.2)) then
              if ((isin(wline(rs)).EQV..false.).AND.(izrot(wline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = wline(rs)
                isin(wline(rs)) = .true.
              end if
            end if
          end do
        end if          
      end do
    end if
    if (have_sjcr.EQV..true.) then
      do i=1,sjcrlst%nr
        if ((sjcrlst%wt(i).gt.sjcrlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(sjcrlst%wt(i).gt.0.0))) then
          do rs=sjcrlst%idx2(i),sjcrlst%idx(i)
            imol = molofrs(rs)
            if (rs.lt.(sjcrlst%idx(i)-nr_crres+1)) cycle
            if ((yline2(rs).gt.0).AND.(yline(rs).gt.0)) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline2(rs)
                isin(yline(rs)) = .true.
              end if
            else if (yline(rs).gt.0) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline(rs)
                isin(yline(rs)) = .true.
              end if
            end if
            if ((fline(rs).gt.0).AND.(rs.gt.rsmol(imol,1))) then
              if ((isin(fline(rs)).EQV..false.).AND.(izrot(fline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = fline(rs)
                isin(fline(rs)) = .true.
              end if
            end if
          end do
        end if          
      end do
    end if
    if (have_ujcr.EQV..true.) then
      do i=1,ujcrlst%nr
        if ((ujcrlst%wt(i).gt.ujcrlst%wt(max(1,i-1))).OR.((i.eq.1).AND.(ujcrlst%wt(i).gt.0.0))) then
          do rs=ujcrlst%idx2(i),ujcrlst%idx(i)
            imol = molofrs(rs)
            if (rs.lt.(ujcrlst%idx(i)-ujmaxsz+1)) cycle
            if ((yline2(rs).gt.0).AND.(yline(rs).gt.0)) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline2(rs)
                isin(yline(rs)) = .true.
              end if
            else if (yline(rs).gt.0) then
              if ((isin(yline(rs)).EQV..false.).AND.(izrot(yline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = yline(rs)
                isin(yline(rs)) = .true.
              end if
            end if
            if ((fline(rs).gt.0).AND.(rs.gt.rsmol(imol,1))) then
              if ((isin(fline(rs)).EQV..false.).AND.(izrot(fline(rs))%alsz.gt.0)) then
                nzmlst = nzmlst + 1
                tmplst(nzmlst) = fline(rs)
                isin(fline(rs)) = .true.
              end if
            end if
          end do
        end if          
      end do
    end if
  end if
  if ((use_dyn.EQV..true.).AND.(fycxyz.eq.1)) then
    do imol=1,nmol
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc = dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.6) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) cycle
!
          rs = atmres(dc_di(imol)%recurs(i,1))
          if ((dc_di(imol)%recurs(i,1).eq.yline(rs)).AND.(yline2(rs).gt.0)) then
            if (isin(yline(rs)).EQV..false.) then
              nzmlst = nzmlst + 1
              tmplst(nzmlst) = yline2(rs)
              isin(yline(rs)) = .true.
            end if
          else if (isin(dc_di(imol)%recurs(i,1)).EQV..false.) then
            nzmlst = nzmlst + 1
            tmplst(nzmlst) = dc_di(imol)%recurs(i,1)
            isin(dc_di(imol)%recurs(i,1)) = .true.
          end if
        end do
      end if
    end do
  else if (fycxyz.eq.2) then
    do i=1,n
      rs = atmres(i)
      if ((i.eq.yline(rs)).AND.(yline2(rs).gt.0)) then
        if (isin(yline(rs)).EQV..false.) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = yline2(rs)
          isin(i) = .true.
        end if
      else if ((izrot(i)%alsz.gt.0).OR.(i.eq.fline(rs)).OR.(i.eq.nucsline(6,rs))) then
        if (isin(i).EQV..true.) cycle
        if (i.eq.atmol(molofrs(rs),1)) cycle
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = i
        isin(i) = .true.
        if ((i.eq.fline(rs)).AND.(izrot(i)%alsz.le.0)) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(rs)%sc(2-shf)
          isin(at(rs)%sc(2-shf)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(rs)%sc(3-shf)
          isin(at(rs)%sc(3-shf)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(rs)%sc(4-shf)
          isin(at(rs)%sc(4-shf)) = .true.
        else if (i.eq.nucsline(6,rs)) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(rs)%sc(1)
          isin(at(rs)%sc(1)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(rs)%sc(2)
          isin(at(rs)%sc(2)) = .true.
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = at(rs)%sc(3)
          isin(at(rs)%sc(3)) = .true.
        end if
      end if
    end do
  end if
  do rs=1,nseq
    if (disulf(rs).gt.0) then
      if (crosslink(crlk_idx(rs))%itstype.le.2) then ! disulfide
        nzmlst = nzmlst + 1
        tmplst(nzmlst) = -at(disulf(rs))%sc(3-shf)
        if (disulf(rs).gt.rs) then
          nzmlst = nzmlst + 1
          tmplst(nzmlst) = -at(disulf(rs))%sc(2-shf)
        end if
      end if
    end if
  end do
!
  if (nzmlst.gt.0) then
!   sort
    allocate(iv(nzmlst,5))
    do i=1,nzmlst
      iv(i,3) = i
    end do
    iv(:,1) = tmplst(1:nzmlst)
    ii = 1
    jj = nzmlst
    call merge_sort(nzmlst,atrue,iv(:,1),iv(:,2),ii,jj)
    allocate(torzmlst(nzmlst))
    torzmlst(:) = iv(:,2)
    deallocate(iv)
  end if
  deallocate(tmplst)
  deallocate(isin)
!
end
!
!----------------------------------------------------------------------------------------------------
!
subroutine do_rama()
!
  use iounit
  use fyoc
  use molecule
  use torsn
  use system
  use grandensembles
  use sequen
  use pdb
!
  implicit none
!
  integer rs,i,bins,iphi,ipsi,imol,moli
  logical ismember
  RTYPE incr
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
  bins = torbins
! add current phi/psi for residue rs to map
!
  do imol=1,nmol
    if (do_tors(moltypid(imol)).EQV..false.) cycle
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    moli = -1
    do i=1,nrmolrama
      if (molrama(i).eq.an_grp_mol(imol)) then
        moli = i
        exit
      end if
    end do
!   terminal residues can never be part of meaningful peptide torsional analysis
    do rs=rsmol(imol,1)+1,rsmol(imol,2)-1
      if ((fline(rs).gt.0).AND.(yline(rs).gt.0)) then
!       overstock terminal bins
        iphi = max(1,min(bins,floor((phi(rs)+180.0)/fyres)+1))
        ipsi = max(1,min(bins,floor((psi(rs)+180.0)/fyres)+1))
        genfymap(iphi,ipsi) = genfymap(iphi,ipsi) + incr
        if (moli.gt.0) then
          molfymap(moli,iphi,ipsi) = molfymap(moli,iphi,ipsi) + incr
        end if
      end if
    end do
  end do
!
  do i=1,nrspecrama ! tested previously for fline/yline
    rs= specrama(i)
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(molofrs(rs))).EQV..true.) then
        if (ismember(ispresent,molofrs(rs)).EQV..false.) cycle
      end if
    end if
!   overstock terminal bins
    iphi = max(1,min(bins,floor((phi(rs)+180.0)/fyres)+1))
    ipsi = max(1,min(bins,floor((psi(rs)+180.0)/fyres)+1))
    if ((iphi.gt.bins).OR.(ipsi.gt.bins)) then
      write(ilog,*) 'Fatal. This means a range exception in th&
 &e pointer arrays (phi,psi) has occurred for residue ',rs,'.'
      call fexit()
    end if
    specfymap(i,iphi,ipsi) = specfymap(i,iphi,ipsi) + incr
  end do
!
! eventually do HN-Halpha J-coupling analysis
  if (ua_model.eq.0) then
    call do_jcoupl()
  end if
!
end
!
!--------------------------------------------------------------------------
!
! we use the average of four different (but similar) parametrizations of 
! the Karplus equation: see
!  Avbelj, F. and Baldwin, R.L., PNAS 100, 10, p5742-5747 (2003)
!
subroutine do_jcoupl()
!
  use iounit
  use torsn
  use molecule
  use polypep
  use fyoc
  use math
  use sequen
  use grandensembles
  use system
  use atoms
  use zmatrix
  use pdb
!
  implicit none
!
  integer imol,rs,cmol
  RTYPE tt,tt2,jc,jc2,getztor,scal
  logical gotone,ismember
!
  if (use_frame_weights.EQV..true.) then
    scal = 1.0*framewts(curframe)
  else
    scal = 1.0
  end if
!
  gotone = .false.
  cmol = 0
  do imol=1,nmol
    if (do_tors(moltypid(imol)).EQV..false.) cycle
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!   terminal residues can never be part of meaningful torsional analysis
    do rs=rsmol(imol,1)+1,rsmol(imol,2)-1
      if (rs.eq.(rsmol(imol,1)+1)) then
        cmol = cmol + 1
        jccnt(cmol) = jccnt(cmol) + 1
      end if
      if (at(rs)%nsc.le.0) cycle
      if ((hni(rs).le.0).OR.(cai(rs).le.0).OR.(ni(rs).le.0)) cycle
      if (at(rs)%sc(1).gt.0) then
        if ((mass(at(rs)%sc(1)).gt.4.0).OR.(iz(1,at(rs)%sc(1)).ne.cai(rs)).OR.(n12(at(rs)%sc(1)).gt.1)) cycle
      end if
      if ((n12(at(rs)%sc(2)).gt.1).OR.(mass(at(rs)%sc(2)).gt.4.0)) then
        tt = getztor(hni(rs),ni(rs),cai(rs),at(rs)%sc(1))/radian
        jc = 0.25*( (6.40*cos(tt)*cos(tt) - 1.40*cos(tt) + 1.90) +&
 &                  (6.60*cos(tt)*cos(tt) - 1.30*cos(tt) + 1.50) +&
 &                  (6.51*cos(tt)*cos(tt) - 1.76*cos(tt) + 1.60) +&
 &                  (7.90*cos(tt)*cos(tt) - 1.05*cos(tt) + 0.65) )
        jc_res(rs) = jc_res(rs) + scal*jc
      else
        tt = getztor(hni(rs),ni(rs),cai(rs),at(rs)%sc(1))/radian
        tt2 = getztor(hni(rs),ni(rs),cai(rs),at(rs)%sc(2))/radian
        jc = 0.25*( (6.40*cos(tt)*cos(tt) - 1.40*cos(tt) + 1.90) +&
 &                  (6.60*cos(tt)*cos(tt) - 1.30*cos(tt) + 1.50) +&
 &                  (6.51*cos(tt)*cos(tt) - 1.76*cos(tt) + 1.60) +&
 &                  (7.90*cos(tt)*cos(tt) - 1.05*cos(tt) + 0.65) )
        jc2 = 0.25*( (6.40*cos(tt2)*cos(tt2) - 1.40*cos(tt2) + 1.90) +&
 &                   (6.60*cos(tt2)*cos(tt2) - 1.30*cos(tt2) + 1.50) +&
 &                   (6.51*cos(tt2)*cos(tt2) - 1.76*cos(tt2) + 1.60) +&
 &                   (7.90*cos(tt2)*cos(tt2) - 1.05*cos(tt2) + 0.65) )
        jc_res(rs) = jc_res(rs) + scal*0.5*(jc + jc2)
      end if
    end do
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a setup routine for 1D-histograms on all torsions, bond angles and lengths
!
subroutine setup_inthists()
!
  use torsn
  use iounit
  use polypep
  use inter
  use sequen
  use params
  use atoms
  use mpistuff
!
  implicit none
!
  integer rs,i,kbl,kba,kdi,kim,freeunit,ir,ii,jj,kk,ll
  RTYPE bli,bai,imi,getztor,getbang,getblen
  character(MAXSTRLEN) fnam
  character(4) b1,b2,b3,b4
#ifdef ENABLE_MPI
  character(3) xpont
  integer lext2
#endif
  logical exists
!
  kbl = 0
  kba = 0
  kim = 0
  do rs=1,nseq
!
    if (do_ints(1).EQV..true.) then
      do i=1,nrsbl(rs)
        kbl = kbl + 1
        bli = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
        bln_mids(kbl) = bli
      end do
    end if
!
    if (do_ints(2).EQV..true.) then
      do i=1,nrsba(rs)
        kba = kba + 1
        bai=getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
        ang_mids(kba) = bai
      end do
    end if
!
    if (do_ints(3).EQV..true.) then
      do i=1,nrsimpt(rs)
        kim = kim + 1
        imi = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),&
 &                    iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
        imp_mids(kim) = imi
      end do
    end if
!
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    lext2 = 3
    if (myrank.eq.0) then
      call int2str(myrank,xpont,lext2)
      fnam =  'N_'//xpont(1:lext2)//'_INTERNAL_COORDS.idx'
    else
      return
    end if
  else if (use_MPIAVG.EQV..true.) then
    if (myrank.ne.0) return
    fnam = 'INTERNAL_COORDS.idx'
  end if
#else
  fnam = 'INTERNAL_COORDS.idx'
#endif
  call strlims(fnam,ii,jj)
  inquire(file=fnam(ii:jj),exist=exists)
  if(exists) then
    ir = freeunit()
    open(unit=ir,file=fnam(ii:jj),status='old')
    close(unit=ir,status='delete')
  end if
  ir=freeunit()
  open(unit=ir,file=fnam(ii:jj),status='new')
!
 66   format(i6,':',1x,2(a4,' '),1x,'(',2(i7,' '),'): ',g12.5)
 67   format(i6,':',1x,3(a4,' '),1x,'(',3(i7,' '),'): ',g12.5)
 68   format(i6,':',1x,4(a4,' '),1x,'(',4(i7,' '),'): ',g12.5)
 69   format(i6,':',1x,4(a4,' '),1x,'(',4(i7,' '),')')
!
  if (do_ints(1).EQV..true.) then
    kbl = 0
    write(ir,*)
    write(ir,*) '----------     Bond Lengths     -----------'
    write(ir,*)
    do rs=1,nseq
      do i=1,nrsbl(rs)
        kbl = kbl + 1
        ii = iaa(rs)%bl(i,1)
        jj = iaa(rs)%bl(i,2)
        b1 = bio_code(b_type(ii))
        b2 = bio_code(b_type(jj))
        write(ir,66) kbl,b1,b2,ii,jj,bln_mids(kbl)
      end do
    end do
  end if
  if (do_ints(2).EQV..true.) then
    kba = 0
    write(ir,*)
    write(ir,*) '----------      Bond Angles     -----------'
    write(ir,*)
    do rs=1,nseq
      do i=1,nrsba(rs)
        kba = kba + 1
        ii = iaa(rs)%ba(i,1)
        jj = iaa(rs)%ba(i,2)
        kk = iaa(rs)%ba(i,3)
        b1 = bio_code(b_type(ii))
        b2 = bio_code(b_type(jj))
        b3 = bio_code(b_type(kk))
        write(ir,67) kba,b1,b2,b3,ii,jj,kk,ang_mids(kba)
      end do
    end do
  end if
  if (do_ints(3).EQV..true.) then
    kim = 0
    write(ir,*)
    write(ir,*) '----------  Improper Torsions   -----------'
    write(ir,*)
    do rs=1,nseq
      do i=1,nrsimpt(rs)
        kim = kim + 1
        ii = iaa(rs)%impt(i,1)
        jj = iaa(rs)%impt(i,2)
        kk = iaa(rs)%impt(i,3)
        ll = iaa(rs)%impt(i,4)
        b1 = bio_code(b_type(ii))
        b2 = bio_code(b_type(jj))
        b3 = bio_code(b_type(kk))
        b4 = bio_code(b_type(ll))
        write(ir,68) kim,b1,b2,b3,b4,ii,jj,kk,ll,imp_mids(kim)
      end do
    end do
  end if
  if (do_ints(4).EQV..true.) then
    kdi = 0
    write(ir,*)
    write(ir,*) '----------   Proper Dihedrals   -----------'
    write(ir,*)
    do rs=1,nseq
      do i=1,nrsdi(rs)
        kdi = kdi + 1
        ii = iaa(rs)%di(i,1)
        jj = iaa(rs)%di(i,2)
        kk = iaa(rs)%di(i,3)
        ll = iaa(rs)%di(i,4)
        b1 = bio_code(b_type(ii))
        b2 = bio_code(b_type(jj))
        b3 = bio_code(b_type(kk))
        b4 = bio_code(b_type(ll))
        write(ir,69) kdi,b1,b2,b3,b4,ii,jj,kk,ll
      end do
    end do
  end if
!
  close(unit=ir)
!
end
!
!-----------------------------------------------------------------------
!
! a simple routine to get 1D-histograms on all torsions, bond angles and bond lengths
!
subroutine do_inthists()
!
  use torsn
  use iounit
  use polypep
  use inter
  use sequen
  use system
  use grandensembles
  use molecule
  use pdb
!
  implicit none
!
  integer rs,i,kbl,kba,kdi,kim,sz(4),it,imol
  RTYPE bli,bai,dii,imi,getztor,getbang,getblen,os(3),effd,incr
  logical ismember
!
  kbl = 0
  kba = 0
  kdi = 0
  kim = 0
  sz(1:4) = intszs(1:4,2)
  os(1:3) = ((sz(1:3)-1)/2.0) + 0.5
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!   
  do imol=1,nmol
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      if (do_ints(1).EQV..true.) then
        do i=1,nrsbl(rs)
          kbl = kbl + 1
          bli = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
          it= floor((bli-(bln_mids(kbl)-os(1)*intres(1)))/intres(1))+1
          if (it.lt.1) it = 1
          if (it.gt.sz(1)) it = sz(1)
          bln_hist(kbl,it) = bln_hist(kbl,it) + incr
        end do
      end if
!
      if (do_ints(2).EQV..true.) then
        do i=1,nrsba(rs)
          kba = kba + 1
          bai=getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
          it= floor((bai-(ang_mids(kba)-os(2)*intres(2)))/intres(2))+1
          if (it.lt.1) it = 1
          if (it.gt.sz(2)) it = sz(2)
          ang_hist(kba,it) = ang_hist(kba,it) + incr
        end do
      end if
!
      if (do_ints(3).EQV..true.) then
        do i=1,nrsimpt(rs)
          kim = kim + 1
          imi = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),&
 &                      iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
          effd = imi-imp_mids(kim)
          if (effd.ge.180.0) effd = effd - 360.0
          if (effd.lt.-180.0) effd = effd + 360.0
          it = floor((effd + os(3)*intres(3))/intres(3)) + 1
          if (it.lt.1) it = 1
          if (it.gt.sz(3)) it = sz(3)
          imp_hist(kim,it) = imp_hist(kim,it) + incr
        end do
      end if
!
      if (do_ints(4).EQV..true.) then
        do i=1,nrsdi(rs)
          kdi = kdi + 1
          dii = getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),&
 &                      iaa(rs)%di(i,3),iaa(rs)%di(i,4))
          it = floor((dii + 180.0 - 0.5*intres(4))/intres(4)) + 1
          if (it.eq.0) then
            it = sz(4)
!           the [0.0:intres(4)/2.0] interval is wrapped into the last bin
!           to avoid the case of 180.0 being an exact bin boundary
          else if (it.lt.1) then
            it = 1 ! should never happen
            write(ilog,*) 'Warning: Dihedral angle out of bounds in do&
 &_inthists().'
          else if (it.gt.sz(4)) then
            it = sz(4) ! should never happen
            write(ilog,*) 'Warning: Dihedral angle out of bounds in do&
 &_inthists().'
          end if
          dih_hist(kdi,it) = dih_hist(kdi,it) + incr
        end do
      end if
    end do
!
  end do
!
end
!
!---------------------------------------------------------------------------
!
subroutine printrama()
!
  use iounit
  use torsn
  use mpistuff
  use molecule
  use fyoc
  use aminos
  use sequen
  use zmatrix
  use polypep
  use atoms
  use pdb
  use system
!
  implicit none
!
  integer rs,i,j,k,rr,ir,lext,freeunit,imol,ii,jj,cmol,jclen
  logical exists,dospecr(nrspecrama),dospecm(nrmolrama)
  character(5) string
  character(60) raspec
  RTYPE totcnt
  RTYPE, ALLOCATABLE:: normer(:)
#ifdef ENABLE_MPI
  character(3) xpont
  integer lext2
#endif
!
 27   format(36000(g12.5,1x))
 28   format(i10,1x,g14.7)
!
  dospecr(:) = .true.
  dospecm(:) = .true.
!
  do k=1,nrspecrama
    totcnt = sum(specfymap(k,1:torbins,1:torbins))
    if (totcnt.gt.0.0) then
      specfymap(k,:,:) = specfymap(k,:,:)/totcnt
    else
      dospecr(k) = .false.
    end if
  end do
  do k=1,nrmolrama
    totcnt = sum(molfymap(k,1:torbins,1:torbins))
    if (totcnt.gt.0.0) then
      molfymap(k,:,:) = molfymap(k,:,:)/totcnt
    else
      dospecm(k) = .false.
    end if
  end do
  totcnt = sum(genfymap(1:torbins,1:torbins))
  if (totcnt.gt.0.0) genfymap(:,:) = genfymap(:,:)/totcnt
!
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    raspec =  'N_'//xpont(1:lext2)//'_RAMACHANDRAN.dat'
  else if (use_MPIAVG.EQV..true.) then
    raspec = 'RAMACHANDRAN.dat'
  end if
#else
  raspec = 'RAMACHANDRAN.dat'
#endif
  call strlims(raspec,ii,jj)
  inquire(file=raspec(ii:jj),exist=exists)
  if(exists) then
    ir = freeunit()
    open(unit=ir,file=raspec(ii:jj),status='old')
    close(unit=ir,status='delete')
  end if
  if (totcnt.gt.0.0) then
    ir=freeunit()
    open(unit=ir,file=raspec(ii:jj),status='new')
    do i=1,torbins
      write(ir,27) (genfymap(j,i),j=1,torbins)
    end do
    close(unit=ir)
  end if
!
  lext = 5
  do rr=1,nrspecrama
    call int2str(specrama(rr),string,lext)
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      raspec =  'N_'//xpont(1:lext2)//'_RESRAMA_'//&
 &string(1:lext)//'.dat'
    else if (use_MPIAVG.EQV..true.) then
      raspec = 'RESRAMA_'//string(1:lext)//'.dat'
    end if
#else
    raspec = 'RESRAMA_'//string(1:lext)//'.dat'
#endif
    call strlims(raspec,ii,jj)
    inquire(file=raspec(ii:jj),exist=exists)
    if(exists) then
      ir = freeunit()
      open(unit=ir,file=raspec(ii:jj),status='old')
      close(unit=ir,status='delete')
    end if
    if (dospecr(rr).EQV..false.) cycle
    ir=freeunit()
    open(unit=ir,file=raspec(ii:jj),status='new')
    do i=1,torbins
      write(ir,27) (specfymap(rr,j,i),j=1,torbins)
    end do
    close(unit=ir)
  end do
!
  do rr=1,nrmolrama
    call int2str(molrama(rr),string,lext)
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      raspec =  'N_'//xpont(1:lext2)//'_MOLRAMA_'//&
 &string(1:lext)//'.dat'
    else if (use_MPIAVG.EQV..true.) then
      raspec = 'MOLRAMA_'//string(1:lext)//'.dat'
    end if
#else
    raspec = 'MOLRAMA_'//string(1:lext)//'.dat'
#endif
    call strlims(raspec,ii,jj)
    inquire(file=raspec(ii:jj),exist=exists)
    if(exists) then
      ir = freeunit()
      open(unit=ir,file=raspec(ii:jj),status='old')
      close(unit=ir,status='delete')
    end if
    if (dospecm(rr).EQV..false.) cycle
    ir=freeunit()
    open(unit=ir,file=raspec(ii:jj),status='new')
    do i=1,torbins
      write(ir,27) (molfymap(rr,j,i),j=1,torbins)
    end do
    close(unit=ir)
  end do
!
  if (sum(jccnt).gt.0) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      raspec =  'N_'//xpont(1:lext2)//'_JCOUPLING.dat'
    else if (use_MPIAVG.EQV..true.) then
      raspec = 'JCOUPLING.dat'
    end if
#else
    raspec = 'JCOUPLING.dat'
#endif
    call strlims(raspec,ii,jj)
    inquire(file=raspec(ii:jj),exist=exists)
    if(exists) then
      ir = freeunit()
      open(unit=ir,file=raspec(ii:jj),status='old')
      close(unit=ir,status='delete')
    end if
    ir=freeunit()
    open(unit=ir,file=raspec(ii:jj),status='new')
!
    jclen = 0
    do imol=1,nmol
      if (do_tors(moltypid(imol)).EQV..false.) cycle
!     terminal residues can never be part of meaningful torsional analysis
      do rs=rsmol(imol,1)+1,rsmol(imol,2)-1
        if (rs.eq.(rsmol(imol,1)+1)) jclen = jclen + 1
        if (at(rs)%nsc.le.0) cycle
        if ((hni(rs).le.0).OR.(cai(rs).le.0).OR.(ni(rs).le.0)) cycle
        if (at(rs)%sc(1).gt.0) then
          if ((mass(at(rs)%sc(1)).gt.4.0).OR.(iz(1,at(rs)%sc(1)).ne.cai(rs)).OR.(n12(at(rs)%sc(1)).gt.1)) cycle
        end if
      end do
    end do
    allocate(normer(jclen))
!
    if (use_frame_weights.EQV..true.) then
      normer(:) = 0.0
      k = 0
      do i=1,framecnt
        if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),angcalc).eq.0)) then
          k = k + 1
          normer(1:jclen) = normer(1:jclen) + framewts2(i)
        end if
      end do
      if (k.ne.jccnt(1)) then
        write(ilog,*) 'Warning. NMR J-coupling constants have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in JCOUPLING.dat or analogous files.'
      end if
    else
      normer(1:jclen) = 1.0*jccnt(1:jclen)
    end if
!
    cmol = 0
    do imol=1,nmol
      if (do_tors(moltypid(imol)).EQV..false.) cycle
!     terminal residues can never be part of meaningful torsional analysis
      do rs=rsmol(imol,1)+1,rsmol(imol,2)-1
        if (rs.eq.(rsmol(imol,1)+1)) cmol = cmol + 1
        if (at(rs)%nsc.le.0) cycle
        if ((hni(rs).le.0).OR.(cai(rs).le.0).OR.(ni(rs).le.0)) cycle
        if (at(rs)%sc(1).gt.0) then
          if ((mass(at(rs)%sc(1)).gt.4.0).OR.(iz(1,at(rs)%sc(1)).ne.cai(rs)).OR.(n12(at(rs)%sc(1)).gt.1)) cycle
        end if
        write(ir,28) rs,jc_res(rs)/(normer(cmol))
      end do
    end do
    close(unit=ir)
    deallocate(normer)
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine prt_inthists()
!
  use torsn
  use iounit
  use inter
  use torsn
  use sequen
  use mpistuff
!
  implicit none
!
  integer i,k,ii,jj
  integer iba,freeunit
  logical exists
#ifdef ENABLE_MPI
  character(3) xpont
  integer lext2
#endif
  character(MAXSTRLEN) fnm
  RTYPE nmz
!
 77   format(100000000(g12.5))
!
#ifdef ENABLE_MPI
  lext2 = 3
#endif
  if (do_ints(1).EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      fnm =  'N_'//xpont(1:lext2)//'_INTHISTS_BL.dat'
    else if (use_MPIAVG.EQV..true.) then
      fnm = 'INTHISTS_BL.dat'
    end if
#else
    fnm = 'INTHISTS_BL.dat'
#endif
    call strlims(fnm,ii,jj)
    inquire(file=fnm(ii:jj),exist=exists)
    if(exists) then
      iba = freeunit()
      open(unit=iba,file=fnm(ii:jj),status='old')
      close(unit=iba,status='delete')
    end if
    iba=freeunit()
    open(unit=iba,file=fnm(ii:jj),status='new')
!
    do i=1,intszs(1,1)
      nmz = sum(bln_hist(i,:))
      if (nmz.gt.0.0) then
        bln_hist(i,:) = bln_hist(i,:)/nmz
      end if
    end do
    do k=1,intszs(1,2)
      write(iba,77) (bln_hist(i,k),i=1,intszs(1,1))
    end do
    close(unit=iba)
  end if
!
  if (do_ints(2).EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      fnm =  'N_'//xpont(1:lext2)//'_INTHISTS_BA.dat'
    else if (use_MPIAVG.EQV..true.) then
      fnm = 'INTHISTS_BA.dat'
    end if
#else
    fnm = 'INTHISTS_BA.dat'
#endif
    call strlims(fnm,ii,jj)
    inquire(file=fnm(ii:jj),exist=exists)
    if(exists) then
      iba = freeunit()
      open(unit=iba,file=fnm(ii:jj),status='old')
      close(unit=iba,status='delete')
    end if
    iba=freeunit()
    open(unit=iba,file=fnm(ii:jj),status='new')
!
    do i=1,intszs(2,1)
      nmz = sum(ang_hist(i,:))
      if (nmz.gt.0.0) then
        ang_hist(i,:) = ang_hist(i,:)/nmz
      end if
    end do
    do k=1,intszs(2,2)
      write(iba,77) (ang_hist(i,k),i=1,intszs(2,1))
    end do
    close(unit=iba)
  end if
!
  if (do_ints(3).EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      fnm =  'N_'//xpont(1:lext2)//'_INTHISTS_IM.dat'
    else if (use_MPIAVG.EQV..true.) then
      fnm = 'INTHISTS_IM.dat'
    end if
#else
    fnm = 'INTHISTS_IM.dat'
#endif
    call strlims(fnm,ii,jj)
    inquire(file=fnm(ii:jj),exist=exists)
    if(exists) then
      iba = freeunit()
      open(unit=iba,file=fnm(ii:jj),status='old')
      close(unit=iba,status='delete')
    end if
    iba=freeunit()
    open(unit=iba,file=fnm(ii:jj),status='new')
!
    do i=1,intszs(3,1)
      nmz = sum(imp_hist(i,:))
      if (nmz.gt.0.0) then
        imp_hist(i,:) = imp_hist(i,:)/nmz
      end if
    end do
    do k=1,intszs(3,2)
      write(iba,77) (imp_hist(i,k),i=1,intszs(3,1))
    end do
    close(unit=iba)
  end if
!
  if (do_ints(4).EQV..true.) then
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      fnm =  'N_'//xpont(1:lext2)//'_INTHISTS_DI.dat'
    else if (use_MPIAVG.EQV..true.) then
      fnm = 'INTHISTS_DI.dat'
    end if
#else
    fnm = 'INTHISTS_DI.dat'
#endif
    call strlims(fnm,ii,jj)
    inquire(file=fnm(ii:jj),exist=exists)
    if(exists) then
      iba = freeunit()
      open(unit=iba,file=fnm(ii:jj),status='old')
      close(unit=iba,status='delete')
    end if
    iba=freeunit()
    open(unit=iba,file=fnm(ii:jj),status='new')
!
    do i=1,intszs(4,1)
      nmz = sum(dih_hist(i,:))
      if (nmz.gt.0.0) then
        dih_hist(i,:) = dih_hist(i,:)/nmz
      end if
    end do
    do k=1,intszs(4,2)
      write(iba,77) (dih_hist(i,k),i=1,intszs(4,1))
    end do
    close(unit=iba)
  end if
!
end
!
!------------------------------------------------------------------------
!
subroutine bb_segment_reader()
!
  use iounit
  use torsn
!
  implicit none
!
  integer i,j,freeunit,it,ii,jj
  logical exists
!
 90   format(36i2)
!
  call strlims(bbsegfile,ii,jj)
  inquire(file=bbsegfile(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=bbsegfile(ii:jj),status='old')
!   note that the file is written such that is visually reproduces a standard Rama.
!   map, hence, the flipped i-loop
    do i=36,1,-1
      read(it,90) (segfymap(i,j),j=1,36)
    end do
    close(unit=it)
  else
    write(ilog,*) 'Fatal. File with backbone segment parsing of Ramachandran map is missing:'
    write(ilog,*) bbsegfile(ii:jj)
    write(ilog,*) 'If this is the desired file, check whether it exists and is readable? Otherwise, &
 &check choice for keyword FMCSC_BBSEGFILE in key-file.'
    call fexit()
  end if
!
end
!
!------------------------------------------------------------------------
!
subroutine bb_segments()
!
  use iounit
  use fyoc
  use molecule
  use polyavg
  use torsn
  use system
  use grandensembles
  use pdb
  use sequen
  use dssps
  use mcsums
!
  implicit none
!
  integer i,k,rs,fi,yi,z1,z2,imol,refrs,mt
  integer seg_len,cur_seg,z3,ptmin,ptmax
  logical ismember
!  integer, ALLOCATABLE:: segrs(:)
  RTYPE incr
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
  ptmax = 0
  ptmin = HUGE(ptmin)
!
! we want to parse the current backbone and identify the segment
! distribution for different basins defined in a datafile
!
!  if (inst_dssp.EQV..true.) then
!    allocate(segrs(nseq))
!    segrs(:) = 0
!  end if
!
  do imol=1,nmol
!
    if (do_tors(moltypid(imol)).EQV..false.) cycle
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
!   some initialization
    seg_len = 0
    cur_seg = 0
    mt = an_grp_mol(imol)
    sgprscnt(mt) = sgprscnt(mt) + 1
!
    do i=rsmol(imol,1)+1,rsmol(imol,2)-1
!
      refrs = rsmol(molangr(mt,1),1) + (i-rsmol(imol,1))
!
      if ((fline(i).gt.0).AND.(yline(i).gt.0)) then
        fi = max(1,min(floor((phi(i)+180.0)/10.0)+1,36))
        yi = max(1,min(floor((psi(i)+180.0)/10.0)+1,36))
        k = segfymap(yi,fi)
      else
        k = 0
      end if
!      if (inst_dssp.EQV..true.) then
!        ptmin = min(i,ptmin)
!        segrs(i) = k
!        ptmax = max(i,ptmax)
!      end if
!
      if (k.gt.0) then
!       for a dipeptide, distribution is trivial and not very useful
        if ((i.eq.rsmol(imol,1)+1).AND.&
 &          (i.eq.rsmol(imol,2)-1)) then
          seg_perrs(k,refrs,1) = seg_perrs(k,refrs,1) + incr
!       first eligible residue always makes a new segment
        else if (i.eq.rsmol(imol,1)+1) then
          seg_len = seg_len + 1
          cur_seg = k
!       last eligible residue always closes current segment
        else if (i.eq.rsmol(imol,2)-1) then
!         elongate one last time and store
          if (k.eq.cur_seg) then
            seg_len = seg_len + 1
            if (seg_len.gt.maxseglen) then
              write(ilog,*) 'WARNING: Segment of type ',cur_seg,' is&
 & too long (Residues ',i-seg_len+1,' to ',i,'). Omitting ...'
            else
              do rs=refrs-seg_len+1,refrs
                seg_perrs(cur_seg,rs,seg_len) =&
 &                   seg_perrs(cur_seg,rs,seg_len) + incr
              end do
            end if
!         terminate previous segment, increment and store new one of length 1
          else
            if ((cur_seg.gt.0).AND.(seg_len.gt.0)) then
              if (seg_len.gt.maxseglen) then
                write(ilog,*) 'WARNING: Segment of type ',cur_seg,' is&
 & too long (Residues ',i-seg_len,' to ',i-1,'). Omitting ...'
              else
                do rs=refrs-seg_len,refrs-1
                  seg_perrs(cur_seg,rs,seg_len) =&
 &                   seg_perrs(cur_seg,rs,seg_len) + incr
                end do
              end if
            end if
            seg_perrs(k,refrs,1) = seg_perrs(k,refrs,1) + incr
          end if
!       elongation of current segment
        else if (k.eq.cur_seg) then
          seg_len = seg_len + 1
!       store old segment (if there was one) and start new one
        else
          if ((cur_seg.gt.0).AND.(seg_len.gt.0)) then
            if (seg_len.gt.maxseglen) then
              write(ilog,*) 'WARNING: Segment of type ',cur_seg,' is t&
 &oo long (Residues ',i-seg_len,' to ',i-1,'). Omitting ...'
            else
              do rs=refrs-seg_len,refrs-1
                seg_perrs(cur_seg,rs,seg_len) =&
 &                 seg_perrs(cur_seg,rs,seg_len) + incr
              end do
            end if
          end if
          seg_len = 1
          cur_seg = k
        end if
      else
!       if we found no valid conformation terminate current segment (if there was one)
        if ((seg_len.gt.0).AND.(cur_seg.gt.0)) then
          do rs=refrs-seg_len,refrs-1
            seg_perrs(cur_seg,rs,seg_len) =&
 &               seg_perrs(cur_seg,rs,seg_len) + incr
          end do
        end if
        seg_len = 0
        cur_seg = 0
      end if
!
    end do
!
!   and secondary structure analysis according to global order parameters z_alpha, z_beta
    call z_secondary(imol,z_alpha(imol),z_beta(imol))
!
!   increment histograms for the latter
    k = an_grp_mol(imol)
    z1 = max(1,min(floor(z_alpha(imol)/0.01) + 1,100))
    z2 = max(1,min(floor(z_beta(imol)/0.01) + 1,100))
    z3 = max(1,min(floor(2.5*(1.75*rgv(imol)/molcontlen(moltypid(imol)))&
 &          **(4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)/tp_binsz) + 1,TPBINS))
    z_hist(k,1,z1) = z_hist(k,1,z1) + incr
    z_hist(k,2,z2) = z_hist(k,2,z2) + incr
    z_2dhist(k,z3,z2) = z_2dhist(k,z3,z2) + incr
    z_hist2(k,z1,z2) = z_hist2(k,z1,z2) + incr
  end do
!
!  if (inst_dssp.EQV..true.) then
! 5446 format(1000000(i2))
! 5447 format(1000000(i7))
!    seg_perrs(2:3,ptmin:ptmax,maxseglen) = 0
!    do i=ptmin,ptmax
!      if ((segrs(i).eq.1).OR.(segrs(i).eq.2)) then
!        if (seg_perrs(4,i,maxseglen).gt.seg_perrs(1,i,maxseglen)) then
!          seg_perrs(3,i,maxseglen) = nstep - seg_perrs(4,i,maxseglen)
!        end if
!        seg_perrs(1,i,maxseglen) = nstep
!      else if (segrs(i).eq.4) then
!        if (seg_perrs(1,i,maxseglen).gt.seg_perrs(4,i,maxseglen)) then
!          seg_perrs(2,i,maxseglen) = nstep - seg_perrs(1,i,maxseglen)
!        end if
!        seg_perrs(4,i,maxseglen) = nstep
!      end if
!      if (seg_perrs(2,i,maxseglen).gt.0) then 
!        write(0,*) i, '-1',nstep,nstep-nint(seg_perrs(2,i,maxseglen)),nint(seg_perrs(2,i,maxseglen))
!      else if (seg_perrs(3,i,maxseglen).gt.0) then 
!        write(0,*) i, '-2',nstep,nstep-nint(seg_perrs(3,i,maxseglen)),nint(seg_perrs(3,i,maxseglen))
!      end if
!    end do
!    if (ptmax.ge.ptmin) then
!      write(0,5446) segrs(ptmin:ptmax)
!    end if
!    deallocate(segrs)
!  end if
!
end
!
!--------------------------------------------------------------------------
!
subroutine prt_bb_segments()
!
  use iounit
  use sequen
  use torsn
  use molecule
  use polyavg
  use mpistuff
  use pdb
  use system
!
  implicit none
!
  integer i,j,k,maxi,freeunit,it,it2,ii,jj,effmoltyp,mt
  logical exists,fly
  character(60) fn
#ifdef ENABLE_MPI
  character(3) xpont
  integer lext2
#endif
  RTYPE, ALLOCATABLE:: segmol(:,:)
  RTYPE, ALLOCATABLE:: normer(:),tot(:,:)
!
 24   format('# Analysis group ',i3,' (Ref.: Mol. ',i5,' / Res. ',i5,'-',&
 &i5,'):')
 46   format('# Molecule ',i5,': Residues ,',i7,'-',i7)
 27   format(8i10)
 28   format(8(g14.6,1x))
!
  effmoltyp = 0
  do k=1,nangrps
    if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
    effmoltyp = effmoltyp + 1
  end do
!
  if (effmoltyp.eq.0) return
!
  allocate(normer(nangrps))
!
  if (use_frame_weights.EQV..false.) then
#ifdef ENABLE_MPI
    lext2 = 3
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      fn =  'N_'//xpont(1:lext2)//'_BB_SEGMENTS_RES.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'BB_SEGMENTS_RES.dat'
    end if
#else
    fn = 'BB_SEGMENTS_RES.dat'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      it = freeunit()
      open(unit=it,file=fn(ii:jj),status='old')
      close(unit=it,status='delete')
    end if
    it=freeunit()
    open(unit=it,file=fn(ii:jj),status='new')
    maxi = 1
    do i=2,maxseglen
      fly = .false.
      do j=1,seg_alsz
        do k=1,NSEG
          if (seg_perrs(k,j,i).gt.0.0) then
            maxi = i
            fly = .true.
            exit
          end if
        end do
        if (fly.EQV..true.) exit
      end do
    end do
   902  format(i8,1000000(1x,i10))
    do k=1,nangrps
      if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
      write(it,24) k,molangr(k,1),rsmol(molangr(k,1),1),&
   &rsmol(molangr(k,1),2)
      do i=rsmol(molangr(k,1),1)+1,rsmol(molangr(k,1),2)-1
        write(it,902) i,((nint(seg_perrs(j,i,mt)), j=1,NSEG), mt=1,maxi)
      end do
    end do
    close(unit=it)
  end if
!
! provide normalized values
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    fn =  'N_'//xpont(1:lext2)//'_BB_SEGMENTS_NORM_RES.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'BB_SEGMENTS_NORM_RES.dat'
  end if
#else
  fn = 'BB_SEGMENTS_NORM_RES.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
  maxi = 1
  do i=2,maxseglen
    fly = .false.
    do j=1,seg_alsz
      do k=1,NSEG
        if (seg_perrs(k,j,i).gt.0.0) then
          maxi = i
          fly = .true.
          exit
        end if
      end do
      if (fly.EQV..true.) exit
    end do
  end do
 903  format(i8,1000000(1x,g12.5))
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),segcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)*molangr(:,2)
      end if
    end do
    if (k.ne.sgprscnt(1)/molangr(1,2)) then
      write(ilog,*) 'Warning. Segment distributions have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output related to segment distributions (BB_SEGMENTS.dat and &
 &related files).'
    end if
  else
    normer(:) = 1.0*sgprscnt(:)
  end if
  do k=1,nangrps
    if (normer(k).le.0.0) cycle
    if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
    write(it,24) k,molangr(k,1),rsmol(molangr(k,1),1),&
 &rsmol(molangr(k,1),2)
    do i=rsmol(molangr(k,1),1)+1,rsmol(molangr(k,1),2)-1
      write(it,903) i,((seg_perrs(j,i,mt)/normer(k), j=1,NSEG), mt=1,maxi)
    end do
  end do
  close(unit=it)
!
! integrate res.-specific information to get per molecule type information
!
  if (use_frame_weights.EQV..false.) then
#ifdef ENABLE_MPI
    lext2 = 3
    if (use_REMC.EQV..true.) then
      call int2str(myrank,xpont,lext2)
      fn =  'N_'//xpont(1:lext2)//'_BB_SEGMENTS.dat'
    else if (use_MPIAVG.EQV..true.) then
      fn = 'BB_SEGMENTS.dat'
    end if
#else
    fn = 'BB_SEGMENTS.dat'
#endif
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      it = freeunit()
      open(unit=it,file=fn(ii:jj),status='old')
      close(unit=it,status='delete')
    end if
    it=freeunit()
    open(unit=it,file=fn(ii:jj),status='new')
  end if
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    fn =  'N_'//xpont(1:lext2)//'_BB_SEGMENTS_NORM.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'BB_SEGMENTS_NORM.dat'
  end if
#else
  fn = 'BB_SEGMENTS_NORM.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it2 = freeunit()
    open(unit=it2,file=fn(ii:jj),status='old')
    close(unit=it2,status='delete')
  end if
  it2=freeunit()
  open(unit=it2,file=fn(ii:jj),status='new')
!
  allocate(segmol(NSEG,maxseglen))
!
  do mt=1,nangrps
    if (normer(mt).le.0.0) cycle
    if (do_tors(moltypid(molangr(mt,1))).EQV..false.) cycle
    maxi = 1
    segmol(:,:) = 0.0
    do k=1,NSEG
      do j=1,maxseglen
        do i=rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
          if (seg_perrs(k,i,j).gt.0.0) then
            if (j.gt.maxi) maxi = j
            segmol(k,j) = segmol(k,j) + seg_perrs(k,i,j)
          end if
        end do
      end do
    end do
    do j=1,maxseglen
      do k=1,NSEG
        segmol(k,j) = segmol(k,j)/(1.0*j)
      end do
    end do
    if (do_tors(moltypid(molangr(mt,1))).EQV..true.) then
      if (use_frame_weights.EQV..false.) then
        write(it,24) mt,molangr(mt,1),rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
        do j=1,maxi
          write(it,27) (nint(segmol(i,j)),i=1,NSEG)
        end do
      end if
      write(it2,24) mt,molangr(mt,1),rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
      do j=1,maxi
        write(it2,28) (segmol(i,j)/normer(mt),i=1,NSEG)
      end do
    end if
  end do
  if (use_frame_weights.EQV..false.) close(unit=it)
  close(unit=it2)
  deallocate(segmol)
!
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    fn =  'N_'//xpont(1:lext2)//'_ZSEC_HIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'ZSEC_HIST.dat'
  end if
#else
  fn = 'ZSEC_HIST.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
 91   format(f7.3,1x,g14.5,1x,g14.5)
!
  allocate(tot(nangrps,4))
  do k=1,nangrps
    tot(k,1) = sum(z_hist(k,1,1:100))
    if (tot(k,1).gt.0.0) z_hist(k,1,:) = z_hist(k,1,:)/tot(k,1)
    tot(k,2) = sum(z_hist(k,2,1:100))
    if (tot(k,2).gt.0.0) z_hist(k,2,:) = z_hist(k,2,:)/tot(k,2)
    tot(k,3) = sum(z_hist2(k,1:100,1:100))
    if (tot(k,3).gt.0.0) z_hist2(k,:,:) = z_hist2(k,:,:)/tot(k,3)
  end do
!
  do k=1,nangrps
    if ((tot(k,1).le.0.0).AND.(tot(k,2).le.0.0)) cycle
    if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
    write(it,24) k,molangr(k,1),rsmol(molangr(k,1),1),&
 &rsmol(molangr(k,1),2)
    do i=1,100
      write(it,91) (i-0.5)*0.01,z_hist(k,1,i),z_hist(k,2,i)
    end do
  end do
  close(unit=it)
!
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    fn =  'N_'//xpont(1:lext2)//'_ZAB_2DHIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'ZAB_2DHIST.dat'
  end if
#else
  fn = 'ZAB_2DHIST.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
 336  format(100(g14.7,1x))
  do k=1,nangrps
    if (tot(k,3).le.0.0) cycle
    if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
    write(it,24) k,molangr(k,1),rsmol(molangr(k,1),1),&
 &rsmol(molangr(k,1),2)
    do i=1,100
      write(it,336) (z_hist2(k,j,i), j=1,100)
    end do
  end do
  close(unit=it)
!
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    fn =  'N_'//xpont(1:lext2)//'_ZBETA_RG.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'ZBETA_RG.dat'
  end if
#else
  fn = 'ZBETA_RG.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=fn(ii:jj),status='old')
    close(unit=it,status='delete')
  end if
  it=freeunit()
  open(unit=it,file=fn(ii:jj),status='new')
 145  format(100(g14.7,1x))
!
  do k=1,nangrps
    tot(k,4) = sum(z_2dhist(k,1:TPBINS,1:100))
    if (tot(k,4).gt.0.0) z_2dhist(k,:,:) = z_2dhist(k,:,:)/tot(k,4)
  end do
!
  do k=1,nangrps
    if (tot(k,4).le.0.0) cycle
    if (do_tors(moltypid(molangr(k,1))).EQV..false.) cycle
    write(it,24) k,molangr(k,1),rsmol(molangr(k,1),1),&
 &rsmol(molangr(k,1),2)
    do j=1,100
      write(it,145) (z_2dhist(k,i,j),i=1,TPBINS)
    end do
  end do
  close(unit=it)
!
  deallocate(tot)
  deallocate(normer)
!
end
!
!----------------------------------------------------------------------
!
! this routine is needed by the global biasing potentials acting
! on z_alpha and z_beta
!
subroutine z_secondary(imol,zav,zbv)
!
  use iounit
  use sequen
  use energies
  use fyoc
  use molecule
  use torsn
  !use mcsums
!
  implicit none
!
  integer i,j,imol,nress !,refmol,freeunit,iu
  RTYPE dis,vv(2),ww(2),uu(2),zav,zbv
!
  zav = 0.0
  zbv = 0.0
  !refmol = 1
!
  if (do_tors(moltypid(imol)).EQV..false.) return
  nress = rsmol(imol,2)-rsmol(imol,1)-1
!
  do i=rsmol(imol,1)+1,rsmol(imol,2)-1
!
    if ((fline(i).le.0).OR.(yline(i).le.0)) cycle
!
!   alpha distance
    vv(1) = phi(i) - par_ZSEC2(1) 
    vv(2) = psi(i) - par_ZSEC2(2)
!   periodic boundary conditions (so to speak)
    do j=1,2
      if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
      if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
    end do
!   now the distance
    dis = sqrt(vv(1)**2 + vv(2)**2)
    if (dis.le.par_ZSEC2(3)) then
!     this residue is fully alpha
      zav = zav + 1.0
    else
!     this residue is contributing less than a full count (most likely ~0.0)
      ww(1) = par_ZSEC2(1) + (par_ZSEC2(3)/dis)*vv(1)
      ww(2) = par_ZSEC2(2) + (par_ZSEC2(3)/dis)*vv(2)
!     PBC #2
      do j=1,2
        if (ww(j).lt.-180.0) ww(j) = ww(j) + 360.0
        if (ww(j).gt.180.0) ww(j) = ww(j) - 360.0
      end do
      uu(1) = phi(i) - ww(1)
      uu(2) = psi(i) - ww(2)
!     PBC #3
      do j=1,2
        if (uu(j).lt.-180.0) uu(j) = uu(j) + 360.0
        if (uu(j).gt.180.0) uu(j) = uu(j) - 360.0
      end do
      zav = zav + exp(-par_ZSEC2(4)*(uu(1)**2+uu(2)**2))
    end if
!   now beta
    vv(1) = phi(i) - par_ZSEC2(5) 
    vv(2) = psi(i) - par_ZSEC2(6)
    do j=1,2
      if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
      if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
    end do
    dis = sqrt(vv(1)**2 + vv(2)**2)
    if (dis.le.par_ZSEC2(7)) then
      zbv = zbv + 1.0
    else
      ww(1) = par_ZSEC2(5) + (par_ZSEC2(7)/dis)*vv(1)
      ww(2) = par_ZSEC2(6) + (par_ZSEC2(7)/dis)*vv(2)
      do j=1,2
        if (ww(j).lt.-180.0) ww(j) = ww(j) + 360.0
        if (ww(j).gt.180.0) ww(j) = ww(j) - 360.0
      end do
      uu(1) = phi(i) - ww(1)
      uu(2) = psi(i) - ww(2)
      do j=1,2
        if (uu(j).lt.-180.0) uu(j) = uu(j) + 360.0
        if (uu(j).gt.180.0) uu(j) = uu(j) - 360.0
      end do
      zbv = zbv + exp(-par_ZSEC2(8)*(uu(1)**2+uu(2)**2))
    end if
  end do
!
  zav = zav/(1.0*nress)
  zbv = zbv/(1.0*nress)
!
  !if (imol.eq.refmol) then
  !  iu = freeunit()
  !  open(unit=iu,file='ZSEC_inst.dat',status='unknown',position='append')  
  !  write(iu,'(I9,2x,F12.6,2x,F12.6)') nstep,zav,zbv
  !  close(unit=iu)
  !end if
!
end
!
!-----------------------------------------------------------------------------
!
subroutine setup_lc_tor()
!
  use iounit
  use energies
  use torsn
  use fyoc
  use movesets
  use system
!
  implicit none
!
  integer i,j,iu,freeunit,t1,t2,hh1,hh2
  logical exists
!
! first read in the coefficients
  call strlims(torlcfile,t1,t2)
  inquire(file=torlcfile(t1:t2),exist=exists)
  if (exists.EQV..false.) then
    write(ilog,*) 'Coefficent file for torsional LCs does not seem t&
 &o exist. Turning off analysis/pot. (',torlcfile(t1:t2),').'
    torlccalc = nsim + 1
    return
  end if
  iu = freeunit()
  open(unit=iu,file=torlcfile(t1:t2),status='old')
  do i=1,ntorlcs
    if (torlcmode.eq.1) then
      read(iu,*) lct_weight(i),(torlc_coeff(i,j),j=1,2*ntorsn)
    else if (torlcmode.eq.2) then
      read(iu,*) lct_weight(i),(torlc_coeff(i,j),j=1,ntorsn)
    end if
  end do 
  close(unit=iu)
!
! if LCT moves are requested, we'll have to read in the inverse, too
  if (use_lctmoves.EQV..true.) then
    call strlims(torlcfile2,t1,t2)
    inquire(file=torlcfile2(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*) 'Coefficent matrix for inverse LCT transform doe&
 &s not exist. Turning off LCT moves (',torlcfile2(t1:t2),') and LCT&
 &potential.'
      use_lctmoves = .false.
      use_LCTOR = .false.
      scale_LCTOR = 0.0
    else
      iu = freeunit()
      open(unit=iu,file=torlcfile2(t1:t2),status='old')
      if (torlcmode.eq.1) then
        if (2*ntorsn.ne.ntorlcs) then
          write(ilog,*) 'Fatal error in setting up LCT moves.'
          write(ilog,*) 'Full transform has to be provided (i.e., nu&
 &mber of LCTs needs to equal 2*N_tor (',torlcfile2(t1:t2),').'
          call fexit()
        else
          do i=1,2*ntorsn
            read(iu,*) (torlc_coeff2(i,j),j=1,ntorlcs)
          end do
        end if
      else if (torlcmode.eq.2) then
        if (ntorsn.ne.ntorlcs) then
          write(ilog,*) 'Fatal error in setting up LCT moves.'
          write(ilog,*) 'Full transform has to be provided (i.e., nu&
 &mber of LCTs needs to equal N_tor (',torlcfile2(t1:t2),').'
          call fexit()
        else
          do i=1,ntorsn
            read(iu,*) (torlc_coeff2(i,j),j=1,ntorlcs)
          end do
        end if
      end if
      close(unit=iu)
    end if
  end if
!
! another thing to do is to set up the natural bounds for the individual LCTs
  if (use_lctmoves.EQV..true.) then
    if (use_LCTOR.EQV..true.) then
      do i=1,ntorlcs
        hh1 = -1
        hh2 = -1
        do j=1,par_LCTOR2(1)
          if (lct_pot(i,j).le.par_LCTOR(3)) then
            if (hh1.eq.(-1)) then
              hh1 = j
            end if
            hh2 = j
          end if
        end do
        if ((hh1.le.(-1)).OR.(hh2.le.(-1)).OR.&
 &           (hh2.lt.hh1)) then
          write(ilog,*) 'Fatal error in setting up LCT moves. Energy&
 & for LCT ',i,' never below threshold (',par_LCTOR(3),'). Turning o&
 &ff LCT moves and energy term.'
          use_LCTOR = .false.
          scale_LCTOR = 0.0
          use_lctmoves = .false.
          exit
        end if
        if (hh1.gt.1) hh1 = hh1 - 1
        if (hh2.lt.par_LCTOR2(1)) hh2 = hh2 + 1
        lct_bnd(i,1) = par_LCTOR(2) + (hh1-1.0)*par_LCTOR(1)
        lct_bnd(i,2) = par_LCTOR(2) + (hh2-1.0)*par_LCTOR(1)
!        write(*,*) i,lct_bnd(i,1),lct_bnd(i,2)
      end do
    end if
  end if
!
! now set boundaries for tor-lc values
!  maxdum = -1.0
!  do i=1,ntorlcs
!    do j=1,ntorsn
!      dum = dum + abs(torlc_coeff(i,2*j-1) + torlc_coeff(i,2*j))
!    end do
!    if (dum.gt.maxdum) maxdum = dum
!  end do
!
!
  do while ((floor(2.0*abs(torlc_params(1))/torlc_params(2))+1).gt.&
 &                                                 MAXTORLCBINS)
    torlc_params(2) = torlc_params(2)*2.0
    write(ilog,*) 'WARNING: Adjusting resolution for torsional LCs. &
 &Increase MAXTORLCBIN and recompile to overcome this problem.'
    if ((floor(2.0*abs(torlc_params(1))/torlc_params(2))+1).le.&
 &                                            MAXTORLCBINS) then
      write(ilog,*)
      write(ilog,*) 'Adjustment complete. Final value: ',&
 &torlc_params(2)
    end if
  end do
!
end
!
!-----------------------------------------------------------------------------
!
subroutine collect_lc_tor()
!
  use iounit
  use torsn
  use aminos
  use sequen
  use math
  use fyoc
  use movesets
  use energies
  use pdb
!
  implicit none
!
  integer i,j,dc,bin,rs
  RTYPE, ALLOCATABLE:: ccc(:)
  RTYPE incr
!
  allocate(ccc(ntorlcs))
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  dc = 0
  do rs=1,nseq
    if (wline(rs).gt.0) then
      dc = dc + 1
      curtvec(dc) = omega(rs)/RADIAN
    end if
    if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
      dc = dc + 1
      curtvec(dc) = phi(rs)/RADIAN
    end if
    if (yline(rs).gt.0) then
      dc = dc + 1
      curtvec(dc) = psi(rs)/RADIAN
    end if
    do i=1,nnucs(rs)
      dc = dc + 1
      curtvec(dc) = nucs(i,rs)/RADIAN
    end do
    do i=1,nchi(rs)
      dc = dc + 1
      curtvec(dc) = chi(i,rs)/RADIAN
    end do
  end do
!
  do i=1,ntorlcs
    if (use_lctmoves.EQV..false.) then
      ccc(i) = 0.0
      do j=1,ntorsn
        if (torlcmode.eq.1) then
          ccc(i) = ccc(i) + torlc_coeff(i,2*j-1)*cos(curtvec(j))&
 &                        + torlc_coeff(i,2*j)*sin(curtvec(j))
        else if (torlcmode.eq.2) then
          ccc(i) = ccc(i) + torlc_coeff(i,j)*curtvec(j)
        end if
      end do
    else
      ccc(i) = lct(i)
    end if
    bin = max(1,min(floor((ccc(i) - torlc_params(1))/torlc_params(2)) + 1,MAXTORLCBINS))
    torlc_hist(i,bin) = torlc_hist(i,bin) + incr
  end do
!
  deallocate(ccc)
!
end
!
!---------------------------------------------------------------------------------
!
subroutine prt_lc_tor()
!
  use iounit
  use torsn
  use mpistuff
  use energies
!
  implicit none
!
  integer i,j,iu,freeunit,ii,jj
  RTYPE dum
  character(100) fn
#ifdef ENABLE_MPI
  character(3) nod
  integer lext
#endif
  logical exists
!
! normalize
  do i=1,ntorlcs
    dum = sum(torlc_hist(i,1:MAXTORLCBINS))
    if (dum.gt.0.0) torlc_hist(i,:) = torlc_hist(i,:)/dum
  end do
  if (dum.le.0.0) return
!
! print out
 44   format(g12.4,1x,20000000(g14.6))
#ifdef ENABLE_MPI
  lext = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,lext)
    fn = 'N_'//nod//'_LC_TORSIONS.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'LC_TORSIONS.dat'
  end if
#else
  fn = 'LC_TORSIONS.dat'
#endif 
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
  do i=1,MAXTORLCBINS
    write(iu,44) (i-0.5)*torlc_params(2)+torlc_params(1),&
 &                            (torlc_hist(j,i),j=1,ntorlcs)
  end do
  close(unit=iu)
!
end
!
!----------------------------------------------------------------------------
!
subroutine lct2fyc()
!
  use iounit
  use torsn
  use fyoc
  use aminos
  use sequen
  use math
!
  implicit none
!
  integer i,j,rs,dc,azero
  RTYPE fc,yc,xc(MAXCHI+6),wc
!
  azero = 0
  dc = 0
! 
  do rs=1,nseq
!   omega
    if (wline(rs).gt.0) then
      dc = dc + 1
      wc = 0.0
      do j=1,ntorlcs
        wc = wc + torlc_coeff2(dc,j)*lct(j)
      end do
      wc = wc*RADIAN
      if ((wc.gt.180.0).OR.(wc.lt.-180.0)) harappa = .true.
      do while (wc.gt.180.0)
        wc = wc - 360.0
      end do
      do while (wc.lt.-180.0)
        wc = wc + 360.0
      end do
!     write(*,*) fc,wc
      call setw(rs,wc,azero)
    end if
!   phi
    if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
      dc = dc + 1
      fc = 0.0
      do j=1,ntorlcs
        fc = fc + torlc_coeff2(dc,j)*lct(j)
      end do
      fc = fc*RADIAN
      if ((fc.gt.180.0).OR.(fc.lt.-180.0)) harappa = .true.
      do while (fc.gt.180.0)
        fc = fc - 360.0
      end do
      do while (fc.lt.-180.0)
        fc = fc + 360.0
      end do
    else if (fline(rs).gt.0) then
      fc = phi(rs)
    end if
    if (yline(rs).gt.0) then
!     psi
      dc = dc + 1
      yc = 0.0
      do j=1,ntorlcs
        yc = yc + torlc_coeff2(dc,j)*lct(j)
      end do
      yc = yc*RADIAN
      if ((yc.gt.180.0).OR.(yc.lt.-180.0)) harappa = .true.
      do while (yc.gt.180.0)
        yc = yc - 360.0
      end do
      do while (yc.lt.-180.0)
        yc = yc + 360.0
      end do
!     write(*,*) fc,yc
    end if
    call setfy(rs,fc,yc,azero)
!   nucs
    do i=1,nnucs(rs)
      dc = dc + 1
      xc(i) = 0.0
      do j=1,ntorlcs
        xc(i) = xc(i) + torlc_coeff2(dc,j)*lct(j)
      end do
      xc(i) = xc(i)*RADIAN
      if ((xc(i).gt.180.0).OR.(xc(i).lt.-180.0)) harappa = .true.
      do while (xc(i).gt.180.0)
        xc(i) = xc(i) - 360.0
      end do
      do while (xc(i).lt.-180.0)
        xc(i) = xc(i) + 360.0
      end do
    end do
    if (nnucs(rs).gt.0) then
!     update the nucs
      call setnucs(rs,xc(1:6),azero)
    end if
!   chis
    do i=1,nchi(rs)
      dc = dc + 1
      xc(i) = 0.0
      do j=1,ntorlcs
        xc(i) = xc(i) + torlc_coeff2(dc,j)*lct(j)
      end do
      xc(i) = xc(i)*RADIAN
      if ((xc(i).gt.180.0).OR.(xc(i).lt.-180.0)) harappa = .true.
      do while (xc(i).gt.180.0)
        xc(i) = xc(i) - 360.0
      end do
      do while (xc(i).lt.-180.0)
        xc(i) = xc(i) + 360.0
      end do
    end do
    if (nchi(rs).gt.0) then
!     update the chi angles
      call setchi(rs,xc(1:MAXCHI),azero)
    end if
  end do
!  write(*,*) phi(2),psi(2)
!
end
!
!-------------------------------------------------------------------------
!
subroutine fyc2lct()
!
  use iounit
  use fyoc
  use aminos
  use sequen
  use torsn
  use math
  use energies
!
  implicit none
!
  integer i,j,dc,rs
  RTYPE dum
!
  dc = 0
  do rs=1,nseq
    if (wline(rs).gt.0) then
      dc = dc + 1
      curtvec(dc) = omega(rs)/RADIAN
    end if
    if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
      dc = dc + 1
      curtvec(dc) = phi(rs)/RADIAN
    end if
    if (yline(rs).gt.0) then
      dc = dc + 1
      curtvec(dc) = psi(rs)/RADIAN
    end if
    do i=1,nnucs(rs)
      dc = dc + 1
      curtvec(dc) = nucs(i,rs)/RADIAN
    end do
    do i=1,nchi(rs)
      dc = dc + 1
      curtvec(dc) = chi(i,rs)/RADIAN
    end do
  end do
!
  do i=1,ntorlcs
    dum = 0.0
    do j=1,ntorsn
      if (torlcmode.eq.1) then
        dum = dum + torlc_coeff(i,2*j-1)*cos(curtvec(j))&
 &                + torlc_coeff(i,2*j)*sin(curtvec(j))
      else if (torlcmode.eq.2) then
        dum = dum + torlc_coeff(i,j)*curtvec(j)
      end if
    end do
    lct(i) = dum
  end do
!
end
! 

