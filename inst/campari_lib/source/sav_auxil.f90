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
!---------------------------------------------------------------------------
!
! this has nothing but the setup routine for the model compound-based limitations
! on the maximum expected SAV in it
! this currently very crudely deals with unknown residues - the entire residue is
! considered a model compound and no dummy hydrogen are created (should be improvable
! [with considerable effort] for supported based topologies, e.g., polypeptides)
!
! the only other functionality is for determining actual SAVs during the course
! of the simulations (analysis tool, required when simulation does not use
! ABSINTH implicit solvent model)
!
! changes to the fundamental way how overlap volumes are calculated have to be made
! in two places: solv_incr_helper and the analysis loop
!
!---------------------------------------------------------------------------------
!
subroutine setup_savol(toponly,totsav)
!
  use atoms
  use aminos
  use polypep
  use iounit
  use inter
  use sequen
  use molecule
  use energies
  use fyoc
  use system
  use params
  use grandensembles
!
  implicit none
!
  integer i,k,ii,kk,rs,rs2,bblst(maxval(at(1:nseq)%nbb+at(1:nseq)%nsc)+50),bbnr,aone,athree
  integer ptlst(20),ptnr,sclst(maxval(at(1:nseq)%nbb+at(1:nseq)%nsc)),scnr,k_12
  RTYPE dis,maxd,totsav,rkk,rii,dis2,svec(3),vs(50),dum,olaps(3)
  logical toponly,ismember
  character(3) resname,resname2
!
  aone = 1
  athree = 3
!
! a few internals so we don't pass real numbers
  vs(1) = 1.09d0
  vs(2) = 109.47d0 
  vs(3) = 119.85d0
  vs(4) = 119.5d0
  vs(5) = 119.35d0
  vs(6) = 127.1d0
  vs(7) = 126.8d0
  vs(8) = 1.01d0
  vs(9) = 180.0d0
  vs(10) = 0.0d0
  vs(11) = 110.3d0
  vs(12) = 145.2d0
  vs(13) = 108.0d0
  vs(14) = -142.3d0
  vs(15) = 120.9d0
!
! this first branch (toponly is true) is the major setup function for initializing 
! i)  the atstv to their exclusion list-derived maxima (to account for all interactions
!     which are not actually computed as NB-interactions during energy evaluations)
! ii) the atsavmaxfr to their model compound-derived constant value
!
  if (toponly.EQV..true.) then
!
!   initialize to atbvol (== volume of atom's solvation shell) 
!
    do i=1,n
      atstv(i) = atbvol(i)
    end do
!
!   account for non-computed interactions (see above)
!
    do rs=1,nseq
!
      resname = amino(seqtyp(rs))
      if (rs.gt.1) resname2 = amino(seqtyp(rs-1))
!
      call calc_bonded_reduc_rs(rs)
!
    end do
!
!   now set up the maximum atom-based acc.-fraction based on the underlying model compounds
!   used for the solvation model
!   these preferably rigid small molecules are more or less fully represented in the
!   the different amino acids, the small differences are overcome by creating and 
!   adding dummy atoms. other than that this part of the code just assembles the
!   necessary atoms for each model compound
!
    do rs=1,nseq
!
!     first prepare by setting up arrays with the relevant atoms (incl. dummies)
!
      resname = amino(seqtyp(rs))
      if (rs.gt.1) then 
        resname2 = amino(seqtyp(rs-1))
      else
        resname2 = '   '
      end if
!
      bbnr = 0
      ptnr = 0
!
!     we'll split this into three separate processes:
!     1) treatment of peptide backbone ((rs-1,rs)-unit assigned to res. rs) OR
!        nucleic acid sugar (all units localized on residue rs: the parsing strategy
!           is more complex since the model compounds have overlap in important atoms
!           (like the glyocsyl bond) -> use of unusual sugar model compounds) OR
!        treatment of full molecule (for simplistic residues)
!     2) treatment of sidechains
!     3) treatment of termini (i.e., exceptions that do not replace one of the
!        previous two groups, but rather have to be considered in addition to them
!
!     1) backbone and simplistic residues
!
      call gen_dummyH_bb(rs,resname,resname2,ptnr,ptlst,vs)
!
!     now that we have generated ptnr dummy atoms in the array ptlst, that have all relevant properties
!     set up, we can go on and add the actually present atoms to make up the full model compound
!
      call gen_complst_bb(rs,resname,resname2,ptnr,ptlst,&
 &                         bbnr,bblst)
!
!     we have assembled all atoms belonging to the model compound in two arrays, ptlst and bblst
!     now finally, initialize the actual array
!
      do i=1,bbnr
        ii = bblst(i)
        atsavmaxfr(ii) = atbvol(ii)
      end do
!
! 211      format('Nr. ',i8,' with radius ',f7.4)
!      write(*,*) 'BBLST for res. ',rs
!      do i=1,bbnr
!        write(*,211) bblst(i),atr(bblst(i))
!      end do
!      write(*,*) 'PTLST for res. ',rs
!      do i=1,ptnr
!        write(*,211) ptlst(i),atr(ptlst(i))
!      end do
!
!     determine SAV as usual
!
      do i=1,bbnr
        ii = bblst(i)
        do k=i+1,bbnr
          kk = bblst(k)
          rkk = atsavred(kk)
          rii = atsavred(ii)
          do k_12=1,n12(ii)
            if (kk.eq.i12(k_12,ii)) then
              rkk = 1.0
              rii = 1.0
              exit
            end if
          end do
          dis = sqrt((x(ii) - x(kk))**2 +&
 &                   (y(ii) - y(kk))**2 +&
 &                   (z(ii) - z(kk))**2)
          maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
          if (dis.lt.maxd) then
            call solv_incr_helper(dis,atsavmaxfr(ii),atsavmaxfr(kk),atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),athree)
          end if
        end do
        do k=1,ptnr
          kk = ptlst(k)
          rkk = atsavred(kk)
          rii = atsavred(ii)
          do k_12=1,n12(kk)
            if (ii.eq.i12(k_12,kk)) then
              rkk = 1.0
              rii = 1.0
              exit
            end if
          end do
          dis = sqrt((x(ii) - x(kk))**2 +&
 &                   (y(ii) - y(kk))**2 +&
 &                   (z(ii) - z(kk))**2)
          maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
          if (dis.lt.maxd) then
            call solv_incr_helper(dis,atsavmaxfr(ii),dum,atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),aone)
          end if
        end do
      end do
! 
!
!     2) sidechains
!
      call gen_all_sc(rs,resname,ptnr,ptlst,&
 &                         scnr,sclst,vs)
!
!     initialization of actual array
!
      do i=1,scnr
        ii = sclst(i)
        atsavmaxfr(ii) = atbvol(ii)
      end do
!
!     determine SAV as usual
!
      do i=1,scnr
        ii = sclst(i)
!
        do k=i+1,scnr
          kk = sclst(k)
          rkk = atsavred(kk)
          rii = atsavred(ii)
          do k_12=1,n12(ii)
            if (kk.eq.i12(k_12,ii)) then
              rkk = 1.0
              rii = 1.0
              exit
            end if
          end do
          dis = sqrt((x(ii) - x(kk))**2 +&
 &                   (y(ii) - y(kk))**2 +&
 &                   (z(ii) - z(kk))**2)
          maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
          if (dis.lt.maxd) then
            call solv_incr_helper(dis,atsavmaxfr(ii),atsavmaxfr(kk),atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),athree)
          end if
        end do
!
        do k=1,ptnr
          kk = ptlst(k)
          rkk = atsavred(kk)
          rii = atsavred(ii)
          do k_12=1,n12(kk)
            if (ii.eq.i12(k_12,kk)) then
              rkk = 1.0
              rii = 1.0
              exit
            end if
          end do
          dis = sqrt((x(ii) - x(kk))**2 +&
 &                   (y(ii) - y(kk))**2 +&
 &                   (z(ii) - z(kk))**2)
          maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
          if (dis.lt.maxd) then
            call solv_incr_helper(dis,atsavmaxfr(ii),dum,atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),aone)
          end if
        end do
!
      end do
!
!     for AIB we actually did this only for one methyl-group (see
!     gen_all_sc(...)) --> now just copy results over
      if (resname.eq.'AIB') then
        atsavmaxfr(at(rs)%sc(2)) = atsavmaxfr(at(rs)%sc(1))
        if (ua_model.eq.0) then
          atsavmaxfr(at(rs)%sc(6)) = atsavmaxfr(at(rs)%sc(3))
          atsavmaxfr(at(rs)%sc(7)) = atsavmaxfr(at(rs)%sc(4))
          atsavmaxfr(at(rs)%sc(8)) = atsavmaxfr(at(rs)%sc(5))
        end if
      end if
!
!
!     cycle out for single-residue molecules or most caps (24 is explicitly defined as 5'-nucleoside, so only sugar
!     to cover)
!     note that unsupporteds are always dealt with in the small molecule vein
!
      if (((seqflag(rs).ge.101).AND.(seqflag(rs).le.103)).OR.(seqflag(rs).eq.24).OR.&
          (resname.eq.'FOR').OR.(resname.eq.'ACE').OR.&
 &        (resname.eq.'NH2').OR.(resname.eq.'NME').OR.(resname.eq.'UNK')) then
        cycle
      end if
!
!
!     3) exceptions (termini) and phosphate unit in polynucleotides
!
!     first phosphates and peptide N-termini, our reference compounds here are
!     dimethylphosphate ((CH3)2PO4-) and methylhydrogenphosphate (HCH3PO4-) for the nucleotides,
!     methylamine for the amino acids, and dimethylamine for the imino acids
!
      if ((rs.eq.rsmol(molofrs(rs),1)).OR.(aminopolty(seqtyp(rs)).eq.'N')) then
!
        call gen_all_Nter(rs,resname,ptnr,ptlst,bbnr,bblst,vs)
!
!       initialization of actual array
!
        do i=1,bbnr
          ii = bblst(i)
          atsavmaxfr(ii) = atbvol(ii)
        end do
!
!       determine SAV as usual
!
        do i=1,bbnr
          ii = bblst(i)
          do k=i+1,bbnr
            kk = bblst(k)
            rkk = atsavred(kk)
            rii = atsavred(ii)
            do k_12=1,n12(ii)
              if (kk.eq.i12(k_12,ii)) then
                rkk = 1.0
                rii = 1.0
                exit
              end if
            end do
            dis = sqrt((x(ii) - x(kk))**2 +&
 &                     (y(ii) - y(kk))**2 +&
 &                     (z(ii) - z(kk))**2)
            maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
            if (dis.lt.maxd) then
              call solv_incr_helper(dis,atsavmaxfr(ii),atsavmaxfr(kk),atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),athree)
            end if
          end do
          do k=1,ptnr
            kk = ptlst(k)
            rkk = atsavred(kk)
            rii = atsavred(ii)
            do k_12=1,n12(kk)
              if (ii.eq.i12(k_12,kk)) then
                rkk = 1.0
                rii = 1.0
                exit
              end if
            end do
            dis = sqrt((x(ii) - x(kk))**2 +&
 &                     (y(ii) - y(kk))**2 +&
 &                     (z(ii) - z(kk))**2)
            maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
            if (dis.lt.maxd) then
              call solv_incr_helper(dis,atsavmaxfr(ii),dum,atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),aone)
            end if
          end do
        end do
!
      end if !whether N-terminal residue
!
!     and now C-terminus, our model compound is acetic acid in all cases
!
      if (rs.eq.rsmol(molofrs(rs),2)) then
!
        call gen_all_Cter(rs,resname,ptnr,ptlst,bbnr,bblst,vs)
!
!       initialize the actual array
!
        do i=1,bbnr
          ii = bblst(i)
          atsavmaxfr(ii) = atbvol(ii)
        end do
!
!       determine SAV as usual
!
        do i=1,bbnr
          ii = bblst(i)
          do k=i+1,bbnr
            kk = bblst(k)
            rkk = atsavred(kk)
            rii = atsavred(ii)
            do k_12=1,n12(ii)
              if (kk.eq.i12(k_12,ii)) then
                rkk = 1.0
                rii = 1.0
                exit
              end if
            end do
            dis = sqrt((x(ii) - x(kk))**2 +&
 &                     (y(ii) - y(kk))**2 +&
 &                     (z(ii) - z(kk))**2)
            maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
            if (dis.lt.maxd) then
              call solv_incr_helper(dis,atsavmaxfr(ii),atsavmaxfr(kk),atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),athree)
            end if
          end do
          do k=1,ptnr
            kk = ptlst(k)
            rkk = atsavred(kk)
            rii = atsavred(ii)
            do k_12=1,n12(kk)
              if (ii.eq.i12(k_12,kk)) then
                rkk = 1.0
                rii = 1.0
                exit
              end if
            end do
            dis = sqrt((x(ii) - x(kk))**2 +&
 &                     (y(ii) - y(kk))**2 +&
 &                     (z(ii) - z(kk))**2)
            maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
            if (dis.lt.maxd) then
              call solv_incr_helper(dis,atsavmaxfr(ii),dum,atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),aone)
            end if
          end do
        end do
!
      end if !whether C-terminal residue
!
!     DEBUGGING ONLY
!      write(*,*) '2BBLST for res. ',rs
!      do i=1,bbnr
!        write(*,211) bblst(i),atr(bblst(i))
!      end do
!      write(*,*) '2PTLST for res. ',rs
!      do i=1,ptnr
!        write(*,211) ptlst(i),atr(ptlst(i))
!        write(*,212) x(ptlst(i)),y(ptlst(i)),z(ptlst(i))
!      end do
!      write(*,212) x(ni(rs)),y(ni(rs)),z(ni(rs))
!      write(*,212) x(ci(rs)),y(ci(rs)),z(ci(rs))
!
    end do
!
!   finally - normalization to give atsavmaxfr!
!
    do i=1,n
      atsavmaxfr(i) = max(0.0,atsavmaxfr(i)/atbvol(i))
!
!      write(*,*) i,': ',atstv(i)/atbvol(i),' ',atsavmaxfr(i)
!
    end do
!
!   after the topology-based assignment, allow overrides
    call read_atpatchfile(aone)
!
!
!
! finally we have the analysis loop to compute SAVs during the simulation if the implicit
! solvent isn't actually used
!
  else  !if toponly - now finally the loop to do SAV analysis during the simulation
!
    do i=1,n
      atsav(i) = atstv(i)
    end do
!
    do rs=1,nseq
!
!     within-residue relevant interactions
      do i=1,nrsintra(rs)
        ii = iaa(rs)%atin(i,1)
        kk = iaa(rs)%atin(i,2)
        rkk = atsavred(kk)
        rii = atsavred(ii)
        dis2 = (x(ii) - x(kk))**2&
 &           + (y(ii) - y(kk))**2&
 &           + (z(ii) - z(kk))**2
        dis = sqrt(dis2)
        maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
        if (dis.lt.maxd) then
          if (par_IMPSOLV2(1).eq.1) then
            if (dis.gt.(maxd-2*atr(kk))) then
              atsav(ii) = atsav(ii) - ((maxd-dis)/(2.0*atr(kk)))*rkk*atvol(kk)
            else if (dis.lt.(atr(ii)+atr(kk))) then
              atsav(ii) = atsav(ii) - (dis/(atr(ii)+atr(kk)))*rkk*atvol(kk)
            else
              atsav(ii) = atsav(ii) - rkk*atvol(kk)
            end if
            if (dis.gt.(maxd-2*atr(ii))) then
              atsav(kk) = atsav(kk) - ((maxd-dis)/(2.0*atr(ii)))*rii*atvol(ii)
            else if (dis.lt.(atr(ii)+atr(kk))) then
              atsav(kk) = atsav(kk) - (dis/(atr(kk)+atr(ii)))*rii*atvol(ii)
            else
             atsav(kk) = atsav(kk) - rii*atvol(ii)
            end if
          else if (par_IMPSOLV2(1).eq.2) then
            call sphere_overlaps(atr(ii),atr(kk),dis,par_IMPSOLV(1),atvol(ii),atvol(kk),olaps(1:3))
            if (olaps(2).gt.0.0) then
              atsav(ii) = atsav(ii) - (olaps(2)-olaps(1))*atsavred(kk)
            end if
            if (olaps(3).gt.0.0) then
              atsav(kk) = atsav(kk) - (olaps(3)-olaps(1))*atsavred(ii)
            end if
          end if
        end if
      end do
!
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(molofrs(rs))).EQV..true.) then
          if (ismember(ispresent,molofrs(rs)).EQV..false.) cycle
        end if
      end if
!
      do rs2=rs+1,nseq
!
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,moltypid(molofrs(rs2))).EQV..true.) then
            if (ismember(ispresent,molofrs(rs2)).EQV..false.) cycle
          end if
        end if
!
        call dis_bound_rs(rs,rs2,svec)
!
        if (rs2.eq.(rs+1)) then
!
!         nearest neighbor relevant interactions
          do i=1,nrsnb(rs)
            ii = iaa(rs)%atnb(i,1)
            kk = iaa(rs)%atnb(i,2)
            rkk = atsavred(kk)
            rii = atsavred(ii)
            call dis_bound4(ii,kk,svec,dis2)
            dis = sqrt(dis2)
            maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
            if (dis.lt.maxd) then
              if (par_IMPSOLV2(1).eq.1) then
                if (dis.gt.(maxd-2*atr(kk))) then
                  atsav(ii) = atsav(ii) - ((maxd-dis)/(2.0*atr(kk)))*rkk*atvol(kk)
                else if (dis.lt.(atr(ii)+atr(kk))) then
                  atsav(ii) = atsav(ii) - (dis/(atr(ii)+atr(kk)))*rkk*atvol(kk)
                else
                  atsav(ii) = atsav(ii) - rkk*atvol(kk)
                end if
                if (dis.gt.(maxd-2*atr(ii))) then
                  atsav(kk) = atsav(kk) - ((maxd-dis)/(2.0*atr(ii)))*rii*atvol(ii)
                else if (dis.lt.(atr(ii)+atr(kk))) then
                  atsav(kk) = atsav(kk) - (dis/(atr(kk)+atr(ii)))*rii*atvol(ii)
                else
                  atsav(kk) = atsav(kk) - rii*atvol(ii)
                end if
              else if (par_IMPSOLV2(1).eq.2) then
                call sphere_overlaps(atr(ii),atr(kk),dis,par_IMPSOLV(1),atvol(ii),atvol(kk),olaps(1:3))
                if (olaps(2).gt.0.0) then
                  atsav(ii) = atsav(ii) - (olaps(2)-olaps(1))*atsavred(kk)
                end if
                if (olaps(3).gt.0.0) then
                  atsav(kk) = atsav(kk) - (olaps(3)-olaps(1))*atsavred(ii)
                end if
              end if
            end if
          end do
! 
        else
!
!         all other interactions
          do i=1,at(rs)%nbb+at(rs)%nsc
            if (i.le.at(rs)%nbb) then
              ii = at(rs)%bb(i)
            else
              ii = at(rs)%sc(i-at(rs)%nbb)
            end if
!
            do k=1,at(rs2)%nbb+at(rs2)%nsc
              if (k.le.at(rs2)%nbb) then
                kk = at(rs2)%bb(k)
              else
               kk = at(rs2)%sc(k-at(rs2)%nbb)
              end if
!
              call dis_bound4(ii,kk,svec,dis2)
              dis = sqrt(dis2)
              maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
              if (dis.lt.maxd) then
                if (par_IMPSOLV2(1).eq.1) then
                  if (dis.gt.(maxd-2*atr(kk))) then
                    atsav(ii)=atsav(ii) - ((maxd-dis)/(2.0*atr(kk)))*atsavred(kk)*atvol(kk)
                  else if (dis.lt.(atr(ii)+atr(kk))) then
                    atsav(ii)=atsav(ii)- (dis/(atr(ii)+atr(kk)))*atsavred(kk)*atvol(kk)
                  else
                    atsav(ii) = atsav(ii) - atsavred(kk)*atvol(kk)
                  end if
                  if (dis.gt.(maxd-2*atr(ii))) then
                    atsav(kk)=atsav(kk) - ((maxd-dis)/(2.0*atr(ii)))*atsavred(ii)*atvol(ii)
                  else if (dis.lt.(atr(ii)+atr(kk))) then
                    atsav(kk)=atsav(kk) - (dis/(atr(kk)+atr(ii)))*atsavred(ii)*atvol(ii)
                  else
                    atsav(kk) = atsav(kk) - atsavred(ii)*atvol(ii)
                  end if
                else if (par_IMPSOLV2(1).eq.2) then
                  call sphere_overlaps(atr(ii),atr(kk),dis,par_IMPSOLV(1),atvol(ii),atvol(kk),olaps(1:3))
                  if (olaps(2).gt.0.0) then
                    atsav(ii) = atsav(ii) - (olaps(2)-olaps(1))*atsavred(kk)
                  end if
                  if (olaps(3).gt.0.0) then
                    atsav(kk) = atsav(kk) - (olaps(3)-olaps(1))*atsavred(ii)
                  end if
                end if
              end if
!
            end do
          end do
        end if
      end do ! over rs2
!
    end do ! over rs
!
  end if
!
  totsav = 0.0
  if (toponly.EQV..true.) then
    do i=1,n
      totsav = totsav + atstv(i)
    end do
  else
    do i=1,n
      totsav = totsav + atsav(i)
    end do
  end if
!
end
!
!-----------------------------------------------------------------------
!
! we're using the following strategy: because the actual SAVs are going
! to computed within the short-range energy evaluation, we need to set up
! the contribution from interactions, which are not computed, separately.
! this topology-based contribution could be found in the same way the non-changing
! interactions are filtered out in setup_srinter(), however, for better code
! clarity we're just going to use the results from setup_srinter() directly
! (--> NEVER change the order setup_savol() and setup_srinter() are called)
!
subroutine calc_bonded_reduc_rs(rs)
!
  use atoms
  use polypep
  use iounit
  use inter
  use sequen
  use molecule
  use energies
!
  implicit none
!
  integer i,k,ii,kk,jj,rs,k_12,athree
  RTYPE dis,maxd,rkk,rii
  logical proceed
!
  athree = 3
!
  do i=1,at(rs)%nbb+at(rs)%nsc
    if (i.le.at(rs)%nbb) then
      ii = at(rs)%bb(i)
    else
      ii = at(rs)%sc(i-at(rs)%nbb)
    end if
    do k=i+1,at(rs)%nbb+at(rs)%nsc
      if (k.le.at(rs)%nbb) then
        kk = at(rs)%bb(k)
      else
        kk = at(rs)%sc(k-at(rs)%nbb)
      end if
      proceed = .true.
      do jj=1,nrsintra(rs)
        if &
 &(((ii.eq.iaa(rs)%atin(jj,1)).AND.(kk.eq.iaa(rs)%atin(jj,2))).OR.&
 & ((kk.eq.iaa(rs)%atin(jj,1)).AND.(ii.eq.iaa(rs)%atin(jj,2)))) then
          proceed = .false.
          exit
        end if
      end do
!
      if (proceed.EQV..true.) then
        rkk = atsavred(kk)
        rii = atsavred(ii)
        do k_12=1,n12(ii)
          if (kk.eq.i12(k_12,ii)) then
            rkk = 1.0
            rii = 1.0
            exit
          end if
        end do
        dis = sqrt((x(ii) - x(kk))**2 +&
 &                 (y(ii) - y(kk))**2 +&
 &                 (z(ii) - z(kk))**2)
        maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
        if (dis.lt.maxd) then
          call solv_incr_helper(dis,atstv(ii),atstv(kk),atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),athree)
        end if
      end if
    end do
  end do
!
! and now through the next neighbor-loop
! note that the code always assumes that rigid connectivity never covers
! more than two adjacent residues
  do i=1,at(rs)%nbb+at(rs)%nsc
!   the last residue in a molecule has no topologically relevant forward next-neighbor
    if (rs.eq.rsmol(molofrs(rs),2)) exit
    if (i.le.at(rs)%nbb) then
      ii = at(rs)%bb(i)
    else
      ii = at(rs)%sc(i-at(rs)%nbb)
    end if
    do k=1,at(rs+1)%nbb+at(rs+1)%nsc
      if (k.le.at(rs+1)%nbb) then
        kk = at(rs+1)%bb(k)
      else
        kk = at(rs+1)%sc(k-at(rs+1)%nbb)
      end if
!
      proceed = .true.
!
      do jj=1,nrsnb(rs)
        if &
 &  (((ii.eq.iaa(rs)%atnb(jj,1)).AND.(kk.eq.iaa(rs)%atnb(jj,2))).OR.&
 &   ((kk.eq.iaa(rs)%atnb(jj,1)).AND.&
 &               (ii.eq.iaa(rs)%atnb(jj,2)))) then
          proceed = .false.
          exit
        end if
      end do
!
      if (proceed.EQV..true.) then
        rkk = atsavred(kk)
        rii = atsavred(ii)
        do k_12=1,n12(ii)
          if (kk.eq.i12(k_12,ii)) then
            rkk = 1.0
            rii = 1.0
            exit
          end if
        end do
        dis = sqrt((x(ii) - x(kk))**2 +&
 &                 (y(ii) - y(kk))**2 +&
 &                 (z(ii) - z(kk))**2)
        maxd = (atr(ii) + atr(kk)) + par_IMPSOLV(1)
        if (dis.lt.maxd) then
          call solv_incr_helper(dis,atstv(ii),atstv(kk),atr(ii),atr(kk),rii,rkk,atvol(ii),atvol(kk),athree)
        end if
      end if
    end do
  end do
!
end
!
!------------------------------------------------------------------------
!
! this routine computes the pairwse SAV volume increment for the supplied parameters
!
! not that option 1 is not consistent with the way SAVs are approximated during the
! actual simulation in case the atoms overlap directly and do not have the same radius
!
! also, the logic of the reduction factors in the fusion limit is different for cases 1 vs 2
!
subroutine solv_incr_helper(dis,valii,valkk,rii,rkk,redii,redkk,volii,volkk,mode)
!
  use iounit
  use energies, ONLY: par_IMPSOLV,par_IMPSOLV2
!
  implicit none
!
  integer, INTENT(IN):: mode
  RTYPE, INTENT(IN):: rii,rkk,volii,volkk,redii,redkk,dis
  RTYPE, INTENT(INOUT):: valii,valkk
!
  RTYPE maxd,olaps(3)
!
  maxd = par_IMPSOLV(1) + rii + rkk
  if (dis.ge.maxd) return
!
  if (par_IMPSOLV2(1).eq.1) then
    if ((mode.eq.1).OR.(mode.eq.3)) then
      if (dis.gt.(maxd-2*rkk)) then
        valii = valii - ((maxd-dis)/(2.0*rkk))*redkk*volkk
      else if (dis.lt.(rii+rkk)) then
        if (rii.ge.rkk) then
          if (dis.gt.(rii-rkk)) then
            valii = valii - ((dis-(rii-rkk))/(2.0*rkk))*redkk*volkk
          end if
        else
          if (dis.gt.(rkk-rii)) then
            valii = valii - (redkk*volkk - ((rkk+rii-dis)/(2.0*rii))*redii*volii)
          else
            valii = valii - (redkk*volkk - redii*volii)
          end if
        end if
      else
        valii = valii - redkk*volkk
      end if
    end if
    if ((mode.eq.2).OR.(mode.eq.3)) then
      if (dis.gt.(maxd-2*rii)) then
        valkk = valkk - ((maxd-dis)/(2.0*rii))*redii*volii
      else if (dis.lt.(rii+rkk)) then
        if (rkk.ge.rii) then
          if (dis.gt.(rkk-rii)) then
            valkk = valkk - ((dis-(rkk-rii))/(2.0*rii))*redii*volii
          end if
        else
          if (dis.gt.(rii-rkk)) then
            valkk = valkk - (redii*volii - ((rii+rkk-dis)/(2.0*rkk))*redkk*volkk)
          else
            valkk = valkk - (redii*volii - redkk*volkk)
          end if
        end if
      else
        valkk = valkk - redii*volii
      end if
    end if
  else if (par_IMPSOLV2(1).eq.2) then
    call sphere_overlaps(rii,rkk,dis,par_IMPSOLV(1),volii,volkk,olaps(1:3))
    if ((mode.eq.1).OR.(mode.eq.3)) then
      if (olaps(2).gt.0.0) then 
        valii = valii - (olaps(2)-olaps(1))*redkk
      end if
    end if
    if ((mode.eq.2).OR.(mode.eq.3)) then
      if (olaps(3).gt.0.0) then
        valkk = valkk - (olaps(3)-olaps(1))*redii
      end if
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported mode for sphere overlap calculation in solv_incr_helper(...). This &
 &is an omission bug.'
    call fexit()
  end if
end
!
!-----------------------------------------------------------------------
!
! the following subroutines cover the following systems in that order
! 1) treatment of peptide backbone ((rs-1,rs)-unit assigned to res. rs) OR
!    nucleic acid sugar (all units localized on residue rs: the parsing strategy
!       is more complex since the model compounds have overlap in important atoms
!       (like the glyocsyl bond) -> use of unusual sugar model compounds) OR
!    treatment of full molecule (for simplistic residues)
!    --> for this system, there are two subroutines, the first one generating
!        dummy hydrogens (gen_dummyH_bb), the second one generating the actual
!        atom lists (gen_complst_bb)
! 2) treatment of sidechains
!    --> for this system, there is only one subroutine (gen_all_sc)
! 3) treatment of termini (i.e., exceptions that do not replace one of the
!    previous two groups, but rather have to be considered in addition to them
!
!-----------------------------------------------------------------------
!
! dummy hydrogen generation for peptide backbone and nucleotide sugars 
!
subroutine gen_dummyH_bb(rs,resname,resname2,ptnr,ptlst,vs)
!
  use atoms
  use iounit
  use molecule
  use sequen
  use aminos
  use polypep 
  use energies
  use fyoc
  use system
!
  implicit none
!
  integer i,rs,ptlst(20),ptnr,refh,shf2,shf
  RTYPE vs(50)
  character(3) resname,resname2
!
! 1) backbone and simplistic residues
!
  ptnr = 0
!
! we'll be needing a lot of dummy hydrogens for the C-alphas of residues rs and rs-1
! since those are fairly generic we'll set them up first
  if ((ua_model.gt.0).OR.((rs.eq.rsmol(molofrs(rs),1)).AND.&
 &        ((aminopolty(seqtyp(rs)).eq.'P').OR.(resname.eq.'ACE').OR.&
 &         (resname.eq.'FOR').OR.(resname.eq.'UNK')))) then
!   do nothing
  else if (resname.eq.'GLY') then 
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-1
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &              refh,vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptnr = 1
  else if (resname.eq.'AIB') then
    refh = at(rs)%sc(at(rs)%nsc)
    ptlst(1) = n+20-3
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &              ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-2
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &              ci(rs),vs(2),-chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptlst(3) = n+20-1
    call genxyz(ptlst(3),cai(rs),vs(1),ni(rs),vs(2),&
 &                  ptlst(2),vs(2),chiral(rs))
    atr(ptlst(3)) = atr(refh)
    atvol(ptlst(3)) = atvol(refh)
    atsavred(ptlst(3)) = atsavred(refh)
    i12(1,ptlst(3)) = cai(rs)
    ptnr = 3
  else if (resname.eq.'PRO') then
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-3
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &              ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-2
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &              refh,vs(2),chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptlst(3) = n+20-1
    call genxyz(ptlst(3),at(rs)%sc(4),vs(1),ni(rs),vs(2),&
 &              at(rs)%sc(9),vs(2),1)
    atr(ptlst(3)) = atr(refh)
    atvol(ptlst(3)) = atvol(refh)
    atsavred(ptlst(3)) = atsavred(refh)
    i12(1,ptlst(3)) = at(rs)%sc(4)
    ptnr = 3
  else if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &      (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &      (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &      (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &      (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &      (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &      (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &      (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &      (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &      (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &      (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &      (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &      (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &      (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &      (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &      (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &      (resname.eq.'HIP').OR.(resname.eq.'KAC').OR.&
 &      (resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &      (resname.eq.'KM3')) then
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-2
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &              ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-1
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &              refh,vs(2),chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptnr = 2
  else if ((resname.eq.'NME').OR.(resname.eq.'ACE').OR.&
 &         (resname.eq.'NH2').OR.(resname.eq.'FOR').OR.&
 &         (resname.eq.'NA+').OR.(resname.eq.'CL-').OR.&
 &         (resname.eq.'K+ ').OR.(resname.eq.'BR-').OR.&
 &         (resname.eq.'CS+').OR.(resname.eq.'I- ').OR.&
 &         (resname.eq.'NH4').OR.(resname.eq.'AC-').OR.&
 &         (resname.eq.'1MN').OR.(resname.eq.'LCP').OR.&
 &         (resname.eq.'2MN').OR.(resname.eq.'NO3').OR.&
 &         (resname.eq.'SPC').OR.(resname.eq.'T3P').OR.&
 &         (resname.eq.'T4P').OR.(resname.eq.'T4E').OR.&
 &         (resname.eq.'T5P').OR.(resname.eq.'PUR').OR.&
 &         (resname.eq.'URE').OR.(resname.eq.'NMA').OR.&
 &         (resname.eq.'NMF').OR.(resname.eq.'ACA').OR.&
 &         (resname.eq.'PPA').OR.(resname.eq.'FOA').OR.&
 &         (resname.eq.'CH4').OR.(resname.eq.'MOH').OR.&
 &         (resname.eq.'PCR').OR.(resname.eq.'DMA').OR.&
 &         (resname.eq.'GDN').OR.(resname.eq.'URA').OR.&
 &         (resname.eq.'THY').OR.(resname.eq.'CYT').OR.&
 &         (resname.eq.'ADE').OR.(resname.eq.'GUA').OR.&
 &         (resname.eq.'PRP').OR.(resname.eq.'NBU').OR.&
 &         (resname.eq.'IBU').OR.(resname.eq.'TOL').OR.&
 &         (resname.eq.'BEN').OR.(resname.eq.'NAP').OR.&
 &         (resname.eq.'EOH').OR.(resname.eq.'MSH').OR.&
 &         (resname.eq.'EMT').OR.(resname.eq.'IMD').OR.&
 &         (resname.eq.'IME').OR.(resname.eq.'MIN').OR.&
 &         (resname.eq.'UNK').OR.(resname.eq.'O2 ')) then
    ptnr = 0
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &         (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &         (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &         (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &         (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &         (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
!   select the first 5*-H as reference
    shf = 0
    if (rs.eq.rsmol(molofrs(rs),2)) shf = 1
    shf2 = 0
    if ((resname.eq.'R5P').OR.(resname.eq.'D5P')) shf2 = 1
    refh = at(rs)%bb(9+shf)
!   now gather dummy atoms (not all are always needed)
    ptlst(1) = n+20-3+shf+shf2
    call genxyz(ptlst(1),nuci(rs,4),vs(1),nuci(rs,5),vs(2),&
 &              refh,vs(2),-1)
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = nuci(rs,4)
    ptnr = 1
!   only if the residue is not 3' we need to add an extra H3*
    if (rs.lt.rsmol(molofrs(rs),2)) then
      ptnr = ptnr + 1
      ptlst(ptnr) = n+20-2+shf+shf2
      call genxyz(ptlst(ptnr),nuci(rs,6),vs(1),nuci(rs,5),&
 &            vs(11),nuci(rs,4),vs(12),0)
      atr(ptlst(ptnr)) = atr(refh)
      atvol(ptlst(ptnr)) = atvol(refh)
      atsavred(ptlst(ptnr)) = atsavred(refh)
      i12(1,ptlst(ptnr)) = nuci(rs,6)
    end if
!   only if a base is present we need to add an extra H1*
    if ((resname.ne.'R5P').AND.(resname.ne.'D5P')) then
      ptnr = ptnr + 1
      ptlst(ptnr) = n+20-1
      call genxyz(ptlst(ptnr),at(rs)%sc(2),vs(1),at(rs)%sc(3),&
 &          vs(13),nuci(rs,5),vs(14),0)
      atr(ptlst(ptnr)) = atr(refh)
      atvol(ptlst(ptnr)) = atvol(refh)
      atsavred(ptlst(ptnr)) = atsavred(refh)
      i12(1,ptlst(ptnr)) = at(rs)%sc(2)
    end if
  else if ((resname.eq.'RIC').OR.(resname.eq.'DIC').OR.&
 &         (resname.eq.'RIU').OR.(resname.eq.'DIU').OR.&
 &         (resname.eq.'RIT').OR.(resname.eq.'DIT').OR.&
 &         (resname.eq.'RIA').OR.(resname.eq.'DIA').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'DIG').OR.&
 &         (resname.eq.'RIB').OR.(resname.eq.'DIB')) then
!   select the first 5*-H as reference
    shf2 = 0
    if ((resname.eq.'RIB').OR.(resname.eq.'DIB')) shf2 = 1
    refh = at(rs)%bb(5)
!   now gather dummy atoms (not all are always needed)
    ptnr =  1
    ptlst(ptnr) = n+20-2+shf2
    call genxyz(ptlst(ptnr),nuci(rs,4),vs(1),nuci(rs,3),&
 &            vs(11),nuci(rs,2),vs(12),0)
    atr(ptlst(ptnr)) = atr(refh)
    atvol(ptlst(ptnr)) = atvol(refh)
    atsavred(ptlst(ptnr)) = atsavred(refh)
    i12(1,ptlst(ptnr)) = nuci(rs,4)
!   only if a base is present we need to add an extra H1*
    if ((resname.ne.'RIB').AND.(resname.ne.'DIB')) then
      ptnr = ptnr + 1
      ptlst(ptnr) = n+20-1
      call genxyz(ptlst(ptnr),at(rs)%sc(2),vs(1),at(rs)%sc(3),&
 &          vs(13),nuci(rs,3),vs(14),0)
      atr(ptlst(ptnr)) = atr(refh)
      atvol(ptlst(ptnr)) = atvol(refh)
      atsavred(ptlst(ptnr)) = atsavred(refh)
      i12(1,ptlst(ptnr)) = at(rs)%sc(2)
    end if
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Implicit solvent does not support re&
 &sidue ',resname,' in pos. ',rs,' at the moment. Check back later.'
      call fexit()
    end if
  end if
!
  if ((ua_model.gt.0).OR.(resname2.eq.'UNK').OR.((rs.eq.rsmol(molofrs(rs),1)).OR.&
 &           (aminopolty(seqtyp(rs)).eq.'N'))) then
!   do nothing 
  else if (resname2.eq.'GLY') then
    refh = at(rs-1)%sc(1)
    ptlst(1+ptnr) = n+20-ptnr-1
    call genxyz(ptlst(1+ptnr),cai(rs-1),vs(1),ci(rs-1),vs(2)&
 &              ,refh,vs(2),-chiral(rs-1))
    atr(ptlst(1+ptnr)) = atr(refh)
    atvol(ptlst(1+ptnr)) = atvol(refh)
    atsavred(ptlst(1+ptnr)) = atsavred(refh)
    i12(1,ptlst(1+ptnr)) = cai(rs-1)
    ptnr = ptnr + 1
  else if (resname2.eq.'AIB') then
    refh = at(rs-1)%sc(at(rs-1)%nsc)
    ptlst(1+ptnr) = n+20-ptnr-3
    call genxyz(ptlst(1+ptnr),cai(rs-1),vs(1),ni(rs-1),vs(2)&
 &              ,ci(rs-1),vs(2),chiral(rs-1))
    atr(ptlst(1+ptnr)) = atr(refh)
    atvol(ptlst(1+ptnr)) = atvol(refh)
    atsavred(ptlst(1+ptnr)) = atsavred(refh)
    i12(1,ptlst(1+ptnr)) = cai(rs-1)
    ptlst(2+ptnr) = n+20-ptnr-2
    call genxyz(ptlst(2+ptnr),cai(rs-1),vs(1),ni(rs-1),vs(2)&
 &              ,ci(rs-1),vs(2),-chiral(rs-1))
    atr(ptlst(2+ptnr)) = atr(refh)
    atvol(ptlst(2+ptnr)) = atvol(refh)
    atsavred(ptlst(2+ptnr)) = atsavred(refh)
    i12(1,ptlst(2+ptnr)) = cai(rs-1)
    ptlst(3+ptnr) = n+20-ptnr-1
    call genxyz(ptlst(3+ptnr),cai(rs-1),vs(1),ci(rs-1),vs(2)&
 &                  ,ptlst(2+ptnr),vs(2),-chiral(rs-1))
    atr(ptlst(3+ptnr)) = atr(refh)
    atvol(ptlst(3+ptnr)) = atvol(refh)
    atsavred(ptlst(3+ptnr)) = atsavred(refh)
    i12(1,ptlst(3+ptnr)) = cai(rs-1)
    ptnr = ptnr + 3
  else if ((resname2.eq.'ALA').OR.(resname2.eq.'ABA').OR.&
 &      (resname2.eq.'NVA').OR.(resname2.eq.'VAL').OR.&
 &      (resname2.eq.'LEU').OR.(resname2.eq.'ILE').OR.&
 &      (resname2.eq.'MET').OR.(resname2.eq.'PHE').OR.&
 &      (resname2.eq.'NLE').OR.(resname2.eq.'SER').OR.&
 &      (resname2.eq.'THR').OR.(resname2.eq.'CYS').OR.&
 &      (resname2.eq.'HIE').OR.(resname2.eq.'HID').OR.&
 &      (resname2.eq.'TYR').OR.(resname2.eq.'TRP').OR.&
 &      (resname2.eq.'ASP').OR.(resname2.eq.'GLU').OR.&
 &      (resname2.eq.'LYS').OR.(resname2.eq.'ARG').OR.&
 &      (resname2.eq.'GLN').OR.(resname2.eq.'ASN').OR.&
 &      (resname2.eq.'GLH').OR.(resname2.eq.'ASH').OR.&
 &      (resname2.eq.'TYO').OR.(resname2.eq.'CYX').OR.&
 &      (resname2.eq.'LYD').OR.(resname2.eq.'PTR').OR.&
 &      (resname2.eq.'SEP').OR.(resname2.eq.'TPO').OR.&
 &      (resname2.eq.'KAC').OR.(resname2.eq.'KM1').OR.&
 &      (resname2.eq.'KM2').OR.(resname2.eq.'KM3').OR.&
 &      (resname2.eq.'PRO').OR.&
 &      (resname2.eq.'ORN').OR.(resname2.eq.'DAB').OR.&
 &      (resname2.eq.'HIP')) then
    refh = at(rs-1)%sc(1)
    ptlst(1+ptnr) = n+20-ptnr-2
    call genxyz(ptlst(1+ptnr),cai(rs-1),vs(1),ni(rs-1),vs(2)&
 &              ,ci(rs-1),vs(2),chiral(rs-1))
    atr(ptlst(1+ptnr)) = atr(refh)
    atvol(ptlst(1+ptnr)) = atvol(refh)
    atsavred(ptlst(1+ptnr)) = atsavred(refh)
    i12(1,ptlst(1+ptnr)) = cai(rs-1)
    ptlst(2+ptnr) = n+20-ptnr-1
    call genxyz(ptlst(2+ptnr),cai(rs-1),vs(1),ci(rs-1),vs(2)&
 &              ,refh,vs(2),-chiral(rs-1))
    atr(ptlst(2+ptnr)) = atr(refh)
    atvol(ptlst(2+ptnr)) = atvol(refh)
    atsavred(ptlst(2+ptnr)) = atsavred(refh)
    i12(1,ptlst(2+ptnr)) = cai(rs-1)
    ptnr = ptnr + 2
  else if ((resname2.eq.'ACE').OR.(resname2.eq.'FOR')) then
!   do nothing
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Implicit solvent does not currently &
 &support residue ',resname2,' in pos ',rs-1,' . Check back later.'
      call fexit()
    end if
  end if
  do i=1,ptnr
    if (ptlst(i).le.n) then
      write(ilog,*) 'Fatal. Abused real atom as holding place fo&
 &r dummy atom. This is most certainly a bug.'
      write(ilog,*) 'Atom ',ptlst(i),' in res.',rs,'.'
      call fexit()
    end if
    n12(ptlst(i)) = 1
  end do
!
! 266      format(i4,1x,i4,1x,3f8.3)
!  do i=1,ptnr
!    write(*,266) rs,ptnr,x(ptlst(i)),y(ptlst(i)),z(ptlst(i))
!  end do
!  write(*,212) x(ci(rs)),y(ci(rs)),z(ci(rs))
!  write(*,212) x(ni(rs-1)),y(ni(rs-1)),z(ni(rs-1))
!
end
!
!-----------------------------------------------------------------------
!
! assembly of atom list for backbone model compounds (incl. nucl. sugars) 
!
subroutine gen_complst_bb(rs,resname,resname2,ptnr,ptlst,bbnr,bblst)
!
  use atoms
  use iounit
  use molecule
  use sequen
  use aminos
  use polypep
  use system 
!
  implicit none
!
  integer i,rs,ptlst(20),ptnr,bblst(maxval(at(1:nseq)%nbb+at(1:nseq)%nsc)+50),bbnr,off,shf,shf2,shf3
  character(3) resname,resname2
!
! now that we have generated ptnr dummy atoms in the array ptlst, that have all relevant properties
! set up, we can go on and add the actually present atoms to make up the full model compound
!
  bbnr = 0
  shf = 0
  shf2 = 0
  shf3 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
  end if
!
! for small molecule residues it's outrageously simple
! note, however, that at least for PPA, EMT, and NBU we have significant internal degrees of freedom
! and therefore use a potentially weak approximation in modelling them all as maximally extended
! (in particular for NBU). note also that all atoms are assumed to be "backbone" for the model compounds
! here, even though a lot of the atoms might actually be in the sc-array.
  if (((seqflag(rs).ge.101).AND.(seqflag(rs).le.103)).OR.(resname.eq.'UNK')) then
    if ((resname.eq.'UNK').AND.(seqpolty(rs).eq.'N').OR.(seqpolty(rs).eq.'P')) then
      write(ilog,*)' Warning. Maximum solvent-accessible volume fractions are computed in simplified fashion for &
 &unsupported residue ',rs,', even though it is recognized as a supported polymer type (use FMCSC_SAVPATCHFILE to fix this).'
    else if ((resname.eq.'UNK').AND.((rsmol(molofrs(rs),1).ne.rs).OR.(rsmol(molofrs(rs),2).ne.rs))) then
      write(ilog,*)' Warning. Maximum solvent-accessible volume fractions are computed in simplified fashion for &
 &unsupported residue ',rs,', which is of an unsupported polymer type (use FMCSC_SAVPATCHFILE to fix this).'
    end if
    do i=1,at(rs)%nbb
      bbnr = bbnr + 1
      bblst(bbnr) = at(rs)%bb(i)
    end do
    do i=1,at(rs)%nsc
      bbnr = bbnr + 1
      bblst(bbnr) = at(rs)%sc(i)
    end do
! we'll treat a.a. termini later
  else if ((rs.eq.rsmol(molofrs(rs),1)).AND.&
 &    ((aminopolty(seqtyp(rs)).eq.'P').OR.(resname.eq.'ACE').OR.&
 &     (resname.eq.'FOR'))) then
!   do nothing
  else if (resname.eq.'GLY') then
    bblst(1) = ni(rs)
    bblst(2) = hni(rs)
    bblst(3) = cai(rs)
    if (ua_model.eq.0) then
      bblst(4) = at(rs)%sc(1)
      bblst(5) = at(rs)%sc(2)
      bbnr = 5
    else
      bbnr = 3
    end if
  else if ((resname.eq.'NME').OR.(resname.eq.'NH2')) then
    do i=1,at(rs)%nbb
      bblst(i) = at(rs)%bb(i)
    end do
    do i=1,at(rs)%nsc
      bblst(at(rs)%nbb+i) = at(rs)%sc(i)
    end do
    bbnr = at(rs)%nsc + at(rs)%nbb
  else if (resname.eq.'AIB') then
    bblst(1) = ni(rs)
    bblst(2) = hni(rs)
    bblst(3) = cai(rs)
    bbnr = 3
  else if (resname.eq.'PRO') then
    bblst(1) = ni(rs)
    bblst(2) = cai(rs)
    if (ua_model.eq.0) then
      bblst(3) = at(rs)%sc(1)
      bbnr = 3
    else
      bbnr = 2
    end if
    ptlst(ptnr+1) = at(rs)%sc(4-shf)
    if (ua_model.eq.0) then
      ptlst(ptnr+2) = at(rs)%sc(9)
      ptlst(ptnr+3) = at(rs)%sc(10)
      ptnr = ptnr + 3
    else
      ptnr = ptnr + 1
    end if
  else if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &      (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &      (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &      (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &      (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &      (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &      (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &      (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &      (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &      (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &      (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &      (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &      (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &      (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &      (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &      (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &      (resname.eq.'HIP').OR.(resname.eq.'KAC').OR.&
 &      (resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &      (resname.eq.'KM3')) then
    bblst(1) = ni(rs)
    bblst(2) = hni(rs)
    bblst(3) = cai(rs)
    if (ua_model.eq.0) then
      bblst(4) = at(rs)%sc(1)
      bbnr = 4
    else
      bbnr = 3
    end if
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &         (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &         (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &         (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &         (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &         (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
!   C5*, C4*, C3*, C2*, C1*, and O4*
    bblst(1) = nuci(rs,4)
    bblst(2) = nuci(rs,5)
    bblst(3) = nuci(rs,6)
    bblst(4) = at(rs)%sc(1)
    bblst(5) = at(rs)%sc(2)
    bblst(6) = at(rs)%sc(3)
    off = 0
    if (ua_model.eq.0) then
      if (rs.eq.rsmol(molofrs(rs),2)) off = 1
!     1H5*, 2H5*, H4*, H3*
      bblst(7) = at(rs)%bb(9+off)
      bblst(8) = at(rs)%bb(10+off)
      bblst(9) = at(rs)%bb(11+off)
      bblst(10) = at(rs)%bb(12+off)
      bbnr = 10
    else
      bbnr = 6
    end if
!   for 3'-residues: O3* and HO3*
    if (rs.eq.rsmol(molofrs(rs),2)) then
      if (moltermid(molofrs(rs),2).ne.1) then
        write(ilog,*) 'Fatal. Encountered unsupported C-terminus&
 & (RES: ',rs,') in setup_savol().'
        call fexit()
      end if     
      bblst(bbnr+1) = at(rs)%bb(9)
      bblst(bbnr+2) = at(rs)%bb(at(rs)%nbb)
      bbnr = bbnr + 2
    end if
!   and finally the rest of the sugar 
    if ((resname.eq.'D5P').OR.(resname.eq.'R5P')) then
      do i=4,at(rs)%nsc
        bbnr = bbnr + 1
        bblst(bbnr) = at(rs)%sc(i)
      end do
    else
      if ((resname.eq.'RPC').OR.(resname.eq.'RPT').OR.&
 &        (resname.eq.'RPA').OR.(resname.eq.'RPG').OR.&
 &        (resname.eq.'RPU')) then
!       O2*
        bblst(bbnr+1) = at(rs)%sc(4)
        bbnr = bbnr + 1
!       H2*, H1*
        if (ua_model.eq.0) then
          bblst(bbnr+1) = at(rs)%sc(5)
          bblst(bbnr+2) = at(rs)%sc(6)
          bbnr = bbnr + 2
        end if
!       HO2*
        bblst(bbnr+1) = at(rs)%sc(7-shf2)
        bbnr = bbnr + 1
      else
!       2xH2*, H1*
        if (ua_model.eq.0) then
          bblst(bbnr+1) = at(rs)%sc(4)
          bblst(bbnr+2) = at(rs)%sc(5)
          bblst(bbnr+3) = at(rs)%sc(6)
          bbnr = bbnr + 3
        end if
      end if
     end if
  else if ((resname.eq.'RIC').OR.(resname.eq.'DIC').OR.&
 &         (resname.eq.'RIU').OR.(resname.eq.'DIU').OR.&
 &         (resname.eq.'RIT').OR.(resname.eq.'DIT').OR.&
 &         (resname.eq.'RIA').OR.(resname.eq.'DIA').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'DIG').OR.&
 &         (resname.eq.'RIB').OR.(resname.eq.'DIB')) then
!   O5*, C5*, C4*, C3*, and all hydrogens
    bbnr = 0
    do i=1,at(rs)%nbb
      bbnr = bbnr + 1
      bblst(bbnr) = at(rs)%bb(i)
    end do
!   C2*, C1*, and O4*
    bblst(bbnr+1) = at(rs)%sc(1)
    bblst(bbnr+2) = at(rs)%sc(2)
    bblst(bbnr+3) = at(rs)%sc(3)
    bbnr = bbnr + 3
!   and finally the rest of the sugar 
    if ((resname.eq.'DIB').OR.(resname.eq.'RIB')) then
      do i=4,at(rs)%nsc
        bbnr = bbnr + 1
        bblst(bbnr) = at(rs)%sc(i)
      end do
    else
      if (ua_model.eq.0) then
        bblst(bbnr+1) = at(rs)%sc(4) ! O2* or 1H2*
        bblst(bbnr+2) = at(rs)%sc(5) ! H2* or 2H2*
        bblst(bbnr+3) = at(rs)%sc(6) ! H1*
        bbnr = bbnr + 3
      else if ((resname.eq.'RIC').OR.(resname.eq.'RIT').OR.&
 &        (resname.eq.'RIA').OR.(resname.eq.'RIG').OR.&
 &        (resname.eq.'RIU')) then
        bblst(bbnr+1) = at(rs)%sc(4) ! O2*
        bbnr = bbnr + 1
      end if
      if ((resname.eq.'RIC').OR.(resname.eq.'RIT').OR.&
 &        (resname.eq.'RIA').OR.(resname.eq.'RIG').OR.&
 &        (resname.eq.'RIU')) then
        bblst(bbnr+1) = at(rs)%sc(7-shf2) ! HO2**
        bbnr = bbnr + 1
      end if
    end if
  end if
!
! now again check on the N-terminal nearest neighbor (if applicable)
  if ((rs.eq.rsmol(molofrs(rs),1)).OR.&
 &           (aminopolty(seqtyp(rs)).eq.'N')) then
!   do nothing
  else if (resname2.eq.'GLY') then
    bblst(1+bbnr) = ci(rs-1)
    bblst(2+bbnr) = oi(rs-1)
    bbnr = bbnr + 2
!   despite this being of minute relevance, we always want to make sure
!   (if possible) that the CH-unit is considered part of a peptide model compound
    if ((rs-1).eq.rsmol(molofrs(rs),1)) then
      bblst(1+bbnr) = cai(rs-1)
      if (ua_model.eq.0) then
        bblst(2+bbnr) = at(rs-1)%sc(1)
        bblst(3+bbnr) = at(rs-1)%sc(2)
        bbnr = bbnr + 3
      else
        bbnr = bbnr + 1
      end if
    else
      ptlst(1+ptnr) = cai(rs-1)
      if (ua_model.eq.0) then
        ptlst(2+ptnr) = at(rs-1)%sc(1)
        ptlst(3+ptnr) = at(rs-1)%sc(2)
        ptnr = ptnr + 3
      else
        ptnr = ptnr + 1
      end if
    end if
  else if ((resname2.eq.'ACE').OR.(resname2.eq.'FOR')) then
    do i=1,at(rs-1)%nbb
      bblst(i+bbnr) = at(rs-1)%bb(i)
    end do
    do i=1,at(rs-1)%nsc
      bblst(at(rs-1)%nbb+i+bbnr) = at(rs-1)%sc(i)
    end do
    bbnr = bbnr + at(rs-1)%nsc + at(rs-1)%nbb
  else if (resname2.eq.'UNK') then
    if (ci(rs-1).gt.0) then
      ptlst(1+ptnr) = ci(rs-1)
      ptnr = ptnr + 1
    end if
    if (oi(rs-1).gt.0) then
      ptlst(1+ptnr) = oi(rs-1)
      ptnr = ptnr + 1
    end if
    if (cai(rs-1).gt.0) then
      ptlst(1+ptnr) = cai(rs-1)
      ptnr = ptnr + 1
    end if
  else if (resname2.eq.'AIB') then
    bblst(1+bbnr) = ci(rs-1)
    bblst(2+bbnr) = oi(rs-1)
    bbnr = bbnr + 2
    if ((rs-1).eq.rsmol(molofrs(rs),1)) then
      bblst(1+bbnr) = cai(rs-1)
      bbnr = bbnr + 1
    else
      ptlst(1+ptnr) = cai(rs-1)
      ptnr = ptnr + 1
    end if
  else if ((resname2.eq.'ALA').OR.(resname2.eq.'ABA').OR.&
 &      (resname2.eq.'NVA').OR.(resname2.eq.'VAL').OR.&
 &      (resname2.eq.'LEU').OR.(resname2.eq.'ILE').OR.&
 &      (resname2.eq.'MET').OR.(resname2.eq.'PHE').OR.&
 &      (resname2.eq.'NLE').OR.(resname2.eq.'SER').OR.&
 &      (resname2.eq.'THR').OR.(resname2.eq.'CYS').OR.&
 &      (resname2.eq.'HIE').OR.(resname2.eq.'HID').OR.&
 &      (resname2.eq.'TYR').OR.(resname2.eq.'TRP').OR.&
 &      (resname2.eq.'ASP').OR.(resname2.eq.'GLU').OR.&
 &      (resname2.eq.'LYS').OR.(resname2.eq.'ARG').OR.&
 &      (resname2.eq.'GLN').OR.(resname2.eq.'ASN').OR.&
 &      (resname2.eq.'GLH').OR.(resname2.eq.'ASH').OR.&
 &      (resname2.eq.'TYO').OR.(resname2.eq.'CYX').OR.&
 &      (resname2.eq.'LYD').OR.(resname2.eq.'PTR').OR.&
 &      (resname2.eq.'SEP').OR.(resname2.eq.'TPO').OR.&
 &      (resname2.eq.'KAC').OR.(resname2.eq.'KM1').OR.&
 &      (resname2.eq.'KM2').OR.(resname2.eq.'KM3').OR.&
 &      (resname2.eq.'PRO').OR.&
 &      (resname2.eq.'ORN').OR.(resname2.eq.'DAB').OR.&
 &      (resname2.eq.'HIP')) then
    bblst(1+bbnr) = ci(rs-1)
    bblst(2+bbnr) = oi(rs-1)
    bbnr = bbnr + 2
    if ((rs-1).eq.rsmol(molofrs(rs),1)) then
      bblst(1+bbnr) = cai(rs-1)
      if (ua_model.eq.0) then
        bblst(2+bbnr) = at(rs-1)%sc(1)
        bbnr = bbnr + 2
      else
        bbnr = bbnr + 1
      end if
    else
      ptlst(1+ptnr) = cai(rs-1)
      if (ua_model.eq.0) then
        ptlst(2+ptnr) = at(rs-1)%sc(1)
        ptnr = ptnr + 2
      else
        ptnr = ptnr + 1
      end if
    end if
  end if
!
! DEBUGGING ONLY
! 211      format('Nr. ',i8,' (',a5,') with radius ',f7.4)
! 213      format('Nr. ',i8,' with radius ',f7.4)
! 212      format('HETATM  ',i3,'  CL  CL- A ',i3,'    ',
! &                           f8.3,f8.3,f8.3)
!  write(*,*) 'BBLST for res. ',rs
!  do i=1,bbnr
!   write(*,211)bblst(i),bio_code(b_type(bblst(i))),atr(bblst(i))
!  end do
!  write(*,*) 'PTLST for res. ',rs
!  do i=1,ptnr
!    write(*,213) ptlst(i),atr(ptlst(i))
!    write(*,212) n+i,n+i,x(ptlst(i)),y(ptlst(i)),z(ptlst(i))
!  end do
!  write(*,212) x(ni(rs)),y(ni(rs)),z(ni(rs))
!  write(*,212) x(ci(rs)),y(ci(rs)),z(ci(rs))
!
end
!
!-----------------------------------------------------------------------
!
! dummy hydrogen generation and atom list generation are more straightforwardly
! dealt with in one shot for sidechains (no neighbor dependence, few dummies needed)
!
subroutine gen_all_sc(rs,resname,ptnr,ptlst,scnr,sclst,vs)
!
  use atoms
  use iounit
  use aminos
  use polypep 
  use energies
  use system
  use sequen
  use fyoc
!
  implicit none
!
  integer i,rs,ptlst(20),ptnr,sclst(maxval(at(1:nseq)%nbb+at(1:nseq)%nsc)),scnr,off,shf
  RTYPE vs(50)
  character(3) resname
!
! similar (less complicated, but a lot more cases) setup for sidechains
!
  scnr = 0
  ptnr = 0
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  if (disulf(rs).gt.0) then
    if ((resname.ne.'CYS').OR.(seqtyp(disulf(rs)).ne.8)) then
      write(ilog,*) 'Solvation setup is incomplete for crosslink between&
 & residues ',rs,' and ',disulf(rs),'.'
      if (use_IMPSOLV.EQV..true.) then
        write(ilog,*) 'This is a fatal omission bug if ABSINTH implicit solvent model is in use.'
        call fexit()
      else if (savcalc.le.nsim) then
        write(ilog,*) 'Disabling all SAV-related analysis.'
        savcalc = nsim+1
      end if
    end if
  end if
!
  if ((resname.eq.'ACE').OR.(resname.eq.'GLY').OR.&
 &    (resname.eq.'FOR').OR.(resname.eq.'NH2').OR.(resname.eq.'NME').OR.&
 &    (resname.eq.'DIB').OR.(resname.eq.'RIB').OR.&
 &    (resname.eq.'D5P').OR.(resname.eq.'R5P').OR.&
 &    (((seqflag(rs).ge.101).AND.(seqflag(rs).le.103)).OR.(resname.eq.'UNK'))) then
    scnr = 0
    ptnr = 0
  else if (resname.eq.'ALA') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(3),vs(2)&
 &              ,at(rs)%sc(4),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(3))
      atvol(ptlst(1)) = atvol(at(rs)%sc(3))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(3))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'ABA') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(4),vs(2)&
 &              ,at(rs)%sc(5),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(3))
      atvol(ptlst(1)) = atvol(at(rs)%sc(3))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(3))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
!    write(*,*) x(ptlst(1)),' ',y(ptlst(1)),' ',z(ptlst(1))
  else if (resname.eq.'VAL') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(5),vs(2)&
 &              ,at(rs)%sc(4),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(5))
      atvol(ptlst(1)) = atvol(at(rs)%sc(5))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(5))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'NVA') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(5),vs(2)&
 &              ,at(rs)%sc(6),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(5))
      atvol(ptlst(1)) = atvol(at(rs)%sc(5))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(5))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'LEU') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(7),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'ILE') then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to account for changes in accessibility. because it's a purely 
!            hydrophobic sidechain, the pre-built chain might not be the best
!            reference state. this is UNfixed at the moment.
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(3),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'NLE') then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to account for changes in accessibility. because it's a purely 
!            hydrophobic sidechain, the extended chain might not be a good
!            reference state. this is UNfixed at the moment.
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(7),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'MET') then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to account for the changing effect of the CH2CH3 on the
!            accessibility of the S-CH3-remaineder. chi-angle is to the
!            most likely value of 180.0 by default (makepept.f).
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(7),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'SER') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(4),vs(2)&
 &              ,at(rs)%sc(5),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(4))
      atvol(ptlst(1)) = atvol(at(rs)%sc(4))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(4))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'SEP') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8),vs(2)&
 &              ,at(rs)%sc(9),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(8))
      atvol(ptlst(1)) = atvol(at(rs)%sc(8))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(8))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'CYS').OR.(resname.eq.'CYX')) then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(4),vs(2)&
 &              ,at(rs)%sc(5),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(4))
      atvol(ptlst(1)) = atvol(at(rs)%sc(4))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(4))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
    if (disulf(rs).gt.0) then
      do i=2-shf,at(disulf(rs))%nsc
        ptnr = ptnr + 1
        ptlst(ptnr) = at(disulf(rs))%sc(i)
      end do
      if (ua_model.eq.0) then
        ptnr = ptnr + 1
        ptlst(ptnr) = n+20-3
        call genxyz(ptlst(ptnr),at(disulf(rs))%sc(2),vs(1),at(disulf(rs))%sc(4),vs(2)&
 &              ,at(disulf(rs))%sc(5),vs(2),1)
        atr(ptlst(ptnr)) = atr(at(disulf(rs))%sc(4))
        atvol(ptlst(ptnr)) = atvol(at(disulf(rs))%sc(4))
        atsavred(ptlst(ptnr)) = atsavred(at(disulf(rs))%sc(4))
        i12(1,ptlst(ptnr)) = at(disulf(rs))%sc(2)
      end if
    end if
  else if (resname.eq.'THR') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(5),vs(2)&
 &              ,at(rs)%sc(3),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(5))
      atvol(ptlst(1)) = atvol(at(rs)%sc(5))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(5))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'TPO') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(9),vs(2)&
 &              ,at(rs)%sc(3),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(9))
      atvol(ptlst(1)) = atvol(at(rs)%sc(9))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(9))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'ASN') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(7),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'GLN') then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to account for the changing effect of the terminal CH3 on 
!            the accessibility of the polar group, adjusted chi-angle in
!            makepept.f to 90.0 for now
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(7),vs(2)&
 &              ,at(rs)%sc(8),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(7))
      atvol(ptlst(1)) = atvol(at(rs)%sc(7))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(7))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'ASP').OR.(resname.eq.'ASH')) then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(7),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'GLU').OR.(resname.eq.'GLH')) then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to account for the changing effect of the terminal CH3 on 
!            the accessibility of the COO, adjusted chi-angle in
!            makepept.f to 90.0 for now
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
     call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(7),vs(2)&
 &              ,at(rs)%sc(8),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(7))
      atvol(ptlst(1)) = atvol(at(rs)%sc(7))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(7))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'LYS').OR.(resname.eq.'LYD')) then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to get account for the changing effect of the aliphatic tail.
!            however, the extended conformation seems the best reference state
!            at least for NH3.
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(7),vs(2)&
 &              ,at(rs)%sc(8),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(7))
      atvol(ptlst(1)) = atvol(at(rs)%sc(7))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(7))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'ORN') then
!   WARNING: see LYS
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(6),vs(2)&
 &              ,at(rs)%sc(7),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'DAB') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(5),vs(2)&
 &              ,at(rs)%sc(6),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(5))
      atvol(ptlst(1)) = atvol(at(rs)%sc(5))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(5))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'ARG') then
!   WARNING: this is one of the sidechains we would actually have to sample
!            to get account for the changing effect of the aliphatic tail.
!            however, the extended conformation is the best reference state
!            given the dominance of the charged guanidino group.
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(7),vs(2)&
 &              ,at(rs)%sc(8),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(7))
      atvol(ptlst(1)) = atvol(at(rs)%sc(7))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(7))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'PHE').OR.(resname.eq.'KM2')) then
!   WARNING: for KM2 case, see LYS
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(9),vs(2)&
 &              ,at(rs)%sc(10),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(9))
      atvol(ptlst(1)) = atvol(at(rs)%sc(9))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(9))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'TYR').OR.(resname.eq.'TYO').OR.(resname.eq.'KAC').OR.(resname.eq.'KM3')) then
!   WARNING: for KAC/KM3 cases, see LYS
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(10),&
 &                          vs(2),at(rs)%sc(11),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(10))
      atvol(ptlst(1)) = atvol(at(rs)%sc(10))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(10))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'PTR') then
!     WARNING: this is adjusted by default to have the P oop to avoid side asymmetry
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(14),&
 &                          vs(2),at(rs)%sc(15),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(14))
      atvol(ptlst(1)) = atvol(at(rs)%sc(14))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(14))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'TRP') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(12),&
 &                            vs(2),at(rs)%sc(13),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(12))
      atvol(ptlst(1)) = atvol(at(rs)%sc(12))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(12))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if ((resname.eq.'HIE').OR.(resname.eq.'HID').OR.(resname.eq.'HIP').OR.(resname.eq.'KM1')) then
!   WARNING: for KM1 case, see LYS
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8),&
 &                            vs(2),at(rs)%sc(9),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(8))
      atvol(ptlst(1)) = atvol(at(rs)%sc(8))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(8))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptnr = 1
    end if
  else if (resname.eq.'PRO') then
    do i=2-shf,at(rs)%nsc
      sclst(i-1+shf) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-1+shf
!   setup partner atoms (including dummies to reproduce model compound)
    if (ua_model.eq.0) then
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(5),vs(2)&
 &              ,at(rs)%sc(6),vs(2),1)
      atr(ptlst(1)) = atr(at(rs)%sc(5))
      atvol(ptlst(1)) = atvol(at(rs)%sc(5))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(5))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptlst(2) = n+20-3
      call genxyz(ptlst(2),at(rs)%sc(4),vs(1),at(rs)%sc(9),vs(2)&
 &              ,at(rs)%sc(10),vs(2),-1)
      atr(ptlst(2)) = atr(at(rs)%sc(9))
      atvol(ptlst(2)) = atvol(at(rs)%sc(9))
      atsavred(ptlst(2)) = atsavred(at(rs)%sc(9))
      i12(1,ptlst(2)) = at(rs)%sc(4)
      ptnr = 2
    end if
  else if (resname.eq.'AIB') then
    sclst(1) = at(rs)%sc(1)
    if (ua_model.eq.0) then
      sclst(2) = at(rs)%sc(3)
      sclst(3) = at(rs)%sc(4)
      sclst(4) = at(rs)%sc(5)
      scnr = 4
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),sclst(1),vs(1),sclst(2),vs(2)&
 &              ,sclst(3),vs(2),1)
      atr(ptlst(1)) = atr(sclst(2))
      atvol(ptlst(1)) = atvol(sclst(2))
      atsavred(ptlst(1)) = atsavred(sclst(2))
      i12(1,ptlst(1)) = sclst(1)
      ptnr = 1
    else
      scnr = 1
    end if
!   remember that for all bases, we use N-methylated forms and that there are no
!   significant degrees freedom within the model compound
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &         (resname.eq.'RIC').OR.(resname.eq.'DIC')) then
    off = 0
    if ((resname.eq.'DPC').OR.(resname.eq.'DIC')) off = 1
    if (ua_model.gt.0) then
      off = 2
      if ((resname.eq.'DPC').OR.(resname.eq.'DIC')) off = 4
    end if
    do i=8-off,at(rs)%nsc
      sclst(i-(7-off)) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-(7-off)
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),-1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptlst(2) = n+20-3
      call genxyz(ptlst(2),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),1)
      atr(ptlst(2)) = atr(at(rs)%sc(6))
      atvol(ptlst(2)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(2)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(2)) = at(rs)%sc(2)
      ptnr = 2
    end if
  else if ((resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &         (resname.eq.'RIU').OR.(resname.eq.'DIU')) then
    off = 0
    if ((resname.eq.'DPU').OR.(resname.eq.'DIU')) off = 1
    if (ua_model.gt.0) then
      off = 2
      if ((resname.eq.'DPU').OR.(resname.eq.'DIU')) off = 4
    end if
    do i=8-off,at(rs)%nsc
      sclst(i-(7-off)) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-(7-off)
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),-1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptlst(2) = n+20-3
      call genxyz(ptlst(2),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),1)
      atr(ptlst(2)) = atr(at(rs)%sc(6))
      atvol(ptlst(2)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(2)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(2)) = at(rs)%sc(2)
      ptnr = 2
    end if
  else if ((resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &         (resname.eq.'RIT').OR.(resname.eq.'DIT')) then
    off = 0
    if ((resname.eq.'DPT').OR.(resname.eq.'DIT')) off = 1
    if (ua_model.gt.0) then
      off = 2
      if ((resname.eq.'DPT').OR.(resname.eq.'DIT')) off = 4
    end if
    do i=8-off,at(rs)%nsc
      sclst(i-(7-off)) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-(7-off)
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),-1)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptlst(2) = n+20-3
      call genxyz(ptlst(2),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),1)
      atr(ptlst(2)) = atr(at(rs)%sc(6))
      atvol(ptlst(2)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(2)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(2)) = at(rs)%sc(2)
      ptnr = 2
    end if
  else if ((resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &         (resname.eq.'RIA').OR.(resname.eq.'DIA')) then
    off = 0
    if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) off = 1
    if (ua_model.gt.0) then
      off = 2
      if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) off = 4
    end if
    do i=8-off,at(rs)%nsc
      sclst(i-(7-off)) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-(7-off)
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),0)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptlst(2) = n+20-3
      call genxyz(ptlst(2),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),1)
      atr(ptlst(2)) = atr(at(rs)%sc(6))
      atvol(ptlst(2)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(2)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(2)) = at(rs)%sc(2)
      ptnr = 2
    end if
  else if ((resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &         (resname.eq.'RIG').OR.(resname.eq.'DIG')) then
    off = 0
    if ((resname.eq.'DPG').OR.(resname.eq.'DIG')) off = 1
    if (ua_model.gt.0) then
      off = 2
      if ((resname.eq.'DPG').OR.(resname.eq.'DIG')) off = 4
    end if
    do i=8-off,at(rs)%nsc
      sclst(i-(7-off)) = at(rs)%sc(i)
    end do
    scnr = at(rs)%nsc-(7-off)
    if (ua_model.eq.0) then
!     setup partner atoms (including dummies to reproduce model compound)
      ptlst(1) = n+20-4
      call genxyz(ptlst(1),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),0)
      atr(ptlst(1)) = atr(at(rs)%sc(6))
      atvol(ptlst(1)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(1)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(1)) = at(rs)%sc(2)
      ptlst(2) = n+20-3
      call genxyz(ptlst(2),at(rs)%sc(2),vs(1),at(rs)%sc(8-off)&
 &                            ,vs(2),at(rs)%sc(6),vs(2),1)
      atr(ptlst(2)) = atr(at(rs)%sc(6))
      atvol(ptlst(2)) = atvol(at(rs)%sc(6))
      atsavred(ptlst(2)) = atsavred(at(rs)%sc(6))
      i12(1,ptlst(2)) = at(rs)%sc(2)
      ptnr = 2
    end if
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*) 'Fatal. Implicit solvent does not support re&
 &sidue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
!   WARNING!!! this check we're only going to do once
!   make sure implementation is COMPLETE before removing it 
    if (savcalc.le.nsim) then
      write(ilog,*)
      write(ilog,*) 'Warning. SAV-analysis not yet supported for&
 & residue ',resname,' in pos. ',rs,' at the moment. Turned off.'
      savcalc = nsim+1
    end if
  end if
!
  do i=1,ptnr
    if ((disulf(rs).eq.0).AND.(ptlst(i).le.n)) then
      write(ilog,*) 'Fatal. Abused real atom as holding place fo&
 &r dummy atom. This is most certainly a bug.'
      write(ilog,*) 'Atom ',ptlst(i),' in res.',rs,'.'
      call fexit()
    else if ((disulf(rs).gt.0).AND.(ptlst(i).le.n)) then
!     do nothing
    else
      n12(ptlst(i)) = 1
    end if
  end do
!
  if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RIC').OR.(resname.eq.'DIC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RIU').OR.(resname.eq.'DIU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RIT').OR.(resname.eq.'DIT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RIA').OR.(resname.eq.'DIA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'RIG').OR.(resname.eq.'DIG')) then
    ptlst(ptnr+1) = at(rs)%sc(2) ! the C1*
    if (ua_model.eq.0) then
      ptlst(ptnr+2) = at(rs)%sc(6) ! the H1*
      ptnr = ptnr + 2
    else
      ptnr = ptnr + 1
    end if
  else
!   do nothing (sanity checked multiple times elsewhere)
  end if
!
end
!
!-----------------------------------------------------------------------
!
! assembly of dummy hydrogen and atom list for N-terminal peptide backbone
! model compounds as well as nucleotide phosphates
!
subroutine gen_all_Nter(rs,resname,ptnr,ptlst,bbnr,bblst,vs)
!
  use atoms
  use iounit
  use molecule
  use sequen
  use aminos
  use energies
  use polypep
  use fyoc 
  use system
  use params
!
  implicit none
!
  integer i,rs,ptlst(20),ptnr,bblst(maxval(at(1:nseq)%nbb+at(1:nseq)%nsc)+50),bbnr,off,off2,refh,shf
  RTYPE vs(50)
  character(3) resname
!
  ptnr = 0
  bbnr = 0
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  if (ua_model.gt.0) then
!   do nothing
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
!
!   let's start by making 1-3 dummy hydrogens for the esterified (soon-to-be) methyl carbons
    off = 0
    if (rs.eq.rsmol(molofrs(rs),1)) off = 2
    off2 = 0
    if (rs.eq.rsmol(molofrs(rs),2)) off2 = 1
    ptnr = 1
    refh = at(rs)%bb(9+off2) ! first H5*
!   the extra H5*
    ptlst(ptnr) = n+20-3+off
    call genxyz(ptlst(ptnr),nuci(rs,4),vs(1),&
 & at(rs)%bb(9+off2),vs(2),at(rs)%bb(10+off2),vs(2),1)
    atr(ptlst(ptnr)) = atr(refh)
    atvol(ptlst(ptnr)) = atvol(refh)
    atsavred(ptlst(ptnr)) = atsavred(refh)
    i12(1,ptlst(ptnr)) = nuci(rs,4)
!   the two extra 3*H
    if (rs.gt.rsmol(molofrs(rs),1)) then
      if (seqpolty(rs-1).eq.'N') then
        off = 0
        if ((nuci(rs,5).le.0).AND.(nuci(rs,6).le.0).AND.(nuci(rs,4).gt.0)) off = 4
        ptnr = ptnr + 1
        ptlst(ptnr) = n+20-2
        call genxyz(ptlst(ptnr),at(rs-1)%bb(4+off),vs(1),&
 &      at(rs)%bb(1),vs(2),at(rs-1)%bb(8+off),vs(2),1)
        atr(ptlst(ptnr)) = atr(refh)
        atvol(ptlst(ptnr)) = atvol(refh)
        atsavred(ptlst(ptnr)) = atsavred(refh)
        i12(1,ptlst(ptnr)) = at(rs-1)%bb(4+off)
        ptnr = ptnr + 1
        ptlst(ptnr) = n+20-1
        call genxyz(ptlst(ptnr),at(rs-1)%bb(4+off),vs(1),&
 &      at(rs-1)%bb(3+off),vs(2),at(rs-1)%bb(8+off),vs(2),1)
        atr(ptlst(ptnr)) = atr(refh)
        atvol(ptlst(ptnr)) = atvol(refh)
        atsavred(ptlst(ptnr)) = atsavred(refh)
        i12(1,ptlst(ptnr)) = at(rs-1)%bb(4+off)
      end if
    end if
  else if (resname.eq.'GLY') then 
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-1
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &            refh,vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptnr = 1
  else if (resname.eq.'AIB') then
    refh = at(rs)%sc(at(rs)%nsc)
    ptlst(1) = n+20-3
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-2
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),-chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptlst(3) = n+20-1
    call genxyz(ptlst(3),cai(rs),vs(1),ni(rs),vs(2),&
 &                ptlst(2),vs(2),chiral(rs))
    atr(ptlst(3)) = atr(refh)
    atvol(ptlst(3)) = atvol(refh)
    atsavred(ptlst(3)) = atsavred(refh)
    i12(1,ptlst(3)) = cai(rs)
    ptnr = 3
  else if (resname.eq.'PRO') then
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-3
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-2
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &            refh,vs(2),chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptlst(3) = n+20-1
    call genxyz(ptlst(3),at(rs)%sc(4),vs(1),ni(rs),vs(2),&
 &            at(rs)%sc(9),vs(2),1)
    atr(ptlst(3)) = atr(refh)
    atvol(ptlst(3)) = atvol(refh)
    atsavred(ptlst(3)) = atsavred(refh)
    i12(1,ptlst(3)) = at(rs)%sc(4)
    ptnr = 3
  else if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &  (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &  (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &  (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &  (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &  (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &  (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &  (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &  (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &  (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &  (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &  (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &  (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &  (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &  (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &  (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &  (resname.eq.'HIP').OR.(resname.eq.'KAC').OR.&
 &  (resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &  (resname.eq.'KM3')) then
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-2
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-1
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &            refh,vs(2),chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptnr = 2
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &N-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end if
!
  do i=1,ptnr
    if (ptlst(i).le.n) then
      write(ilog,*) 'Fatal. Abused real atom as holding place &
 &for dummy atom. This is most certainly a bug.'
      write(ilog,*) 'Atom ',ptlst(i),' in res.',rs,'.'
      call fexit()
    end if
    n12(ptlst(i)) = 1
  end do
!
! nucleotides:
! phosphate will always be covered (caps already cycled out for)
! peptides:
! only the NH2/3-group will be directly covered, unless it's a free
! amino acid, in which case NHxCHx will be covered
  if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
    bblst(1) = nuci(rs,1)
    bblst(2) = nuci(rs,2)
    bblst(3) = nuci(rs,3)
    bblst(4) = at(rs)%bb(4)
    bblst(5) = at(rs)%bb(5)
    bbnr = 5
    if (rs.eq.rsmol(molofrs(rs),1)) then
      off = 0
      if (rs.eq.rsmol(molofrs(rs),2)) off = 1
      if (moltermid(molofrs(rs),1).ne.1) then
        write(ilog,*) 'Fatal. Encountered unsupported N-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
      if (moltermid(molofrs(rs),2).ne.1) then
        write(ilog,*) 'Fatal. Encountered unsupported C-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
      bblst(6) = at(rs)%bb(at(rs)%nbb-off)
      bbnr = 6
    end if
  else if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &    (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &    (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &    (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &    (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &    (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &    (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &    (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &    (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &    (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &    (resname.eq.'GLY').OR.(resname.eq.'AIB').OR.&
 &    (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &    (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &    (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &    (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &    (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &    (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &    (resname.eq.'HIP').OR.(resname.eq.'KAC').OR.&
 &    (resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &    (resname.eq.'KM3')) then
    bblst(1) = ni(rs)
    if (rs.eq.rsmol(molofrs(rs),2)) then
      bblst(2) = cai(rs)
      bblst(3) = at(rs)%bb(6)
      bblst(4) = at(rs)%bb(7)
      bbnr = 4
      if (ua_model.eq.0) then
        if (resname.eq.'GLY') then
          bblst(5) = at(rs)%sc(1)
          bblst(6) = at(rs)%sc(2)
          bbnr = 6
        else if (resname.ne.'AIB') then
          bblst(5) = at(rs)%sc(1)
          bbnr = 5
        end if
      end if
      if (moltermid(molofrs(rs),1).eq.1) then
        bblst(bbnr+1) = at(rs)%bb(8)
        bbnr = bbnr + 1
      else if (moltermid(molofrs(rs),2).eq.2) then
!       do nothing
      else
        write(ilog,*) 'Fatal. Encountered unsupported N-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
    else
      bblst(2) = at(rs)%bb(5)
      bblst(3) = at(rs)%bb(6)
      if (moltermid(molofrs(rs),1).eq.1) then
        bblst(4) = at(rs)%bb(7)
        bbnr = 4
      else if (moltermid(molofrs(rs),1).eq.2) then
        bbnr = 3
      else
        write(ilog,*) 'Fatal. Encountered unsupported N-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
    end if
  else if (resname.eq.'PRO') then
    bblst(1) = ni(rs)
    if (rs.eq.rsmol(molofrs(rs),2)) then
      bblst(2) = at(rs)%bb(6)
      bblst(3) = cai(rs)
      bbnr = 3
      if (ua_model.eq.0) then
        bblst(bbnr+1) = at(rs)%sc(1)
        bbnr = bbnr + 1
      end if
      if (moltermid(molofrs(rs),1).eq.1) then
        bblst(bbnr+1) = at(rs)%bb(7)
        bbnr = bbnr+1
      else if (moltermid(molofrs(rs),1).eq.2) then
!       do nothing else
      else
        write(ilog,*) 'Fatal. Encountered unsupported N-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
      if (moltermid(molofrs(rs),2).ne.1) then
        write(ilog,*) 'Fatal. Encountered unsupported C-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
    else
      bblst(2) = at(rs)%bb(5)
      bbnr = 2
      if (moltermid(molofrs(rs),1).eq.1) then
        bblst(3) = at(rs)%bb(6)
        bbnr = 3
      else if (moltermid(molofrs(rs),2).eq.2) then
        bbnr = 2
      else
        write(ilog,*) 'Fatal. Encountered unsupported N-termin&
 &us (RES: ',rs,') in setup_savol().'
        call fexit()
      end if
    end if
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &N-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end if
! and now let's assemble the full partner atom list
  if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &    (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &    (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &    (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &    (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &    (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &    (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &    (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &    (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &    (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &    (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &    (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &    (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &    (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &    (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &    (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &    (resname.eq.'HIP').OR.(resname.eq.'KAC').OR.&
 &    (resname.eq.'KM1').OR.(resname.eq.'KM2').OR.&
 &    (resname.eq.'KM3')) then
    if (rs.eq.rsmol(molofrs(rs),2)) then
      ptlst(1+ptnr) = cai(rs)
      if (ua_model.eq.0) then
        ptlst(2+ptnr) = at(rs)%sc(1)
        ptnr = ptnr + 2
      else
        ptnr = ptnr + 1
      end if
    end if
  else if (resname.eq.'GLY') then
    if (rs.ne.rsmol(molofrs(rs),2)) then
      ptlst(1+ptnr) = cai(rs)
      if (ua_model.eq.0) then
        ptlst(2+ptnr) = at(rs)%sc(1)
        ptlst(3+ptnr) = at(rs)%sc(2)
        ptnr = ptnr + 3
      else
        ptnr = ptnr + 1
      end if
    end if
  else if (resname.eq.'AIB') then
    if (rs.eq.rsmol(molofrs(rs),2)) then
      ptlst(1+ptnr) = cai(rs)
      ptnr = ptnr + 1
    end if
  else if (resname.eq.'PRO') then
    ptlst(1+ptnr) = at(rs)%sc(4-shf)
    if (ua_model.eq.0) then
      ptlst(2+ptnr) = at(rs)%sc(9)
      ptlst(3+ptnr) = at(rs)%sc(10)
      if (rs.eq.rsmol(molofrs(rs),2)) then
        ptnr = ptnr + 3
      else
        ptlst(4+ptnr) = cai(rs)
        ptlst(5+ptnr) = at(rs)%sc(1)
        ptnr = ptnr + 5
      end if
    else
      if (rs.eq.rsmol(molofrs(rs),2)) then
        ptnr = ptnr + 1
      else
        ptlst(2+ptnr) = cai(rs)
        ptnr = ptnr + 2
      end if
    end if
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
    off = 0
!   always add the C5* and the two H5*
    if (rs.eq.rsmol(molofrs(rs),2)) off = 1
    ptlst(1+ptnr) = at(rs)%bb(6)
    if (ua_model.eq.0) then
      ptlst(2+ptnr) = at(rs)%bb(9+off)
      ptlst(3+ptnr) = at(rs)%bb(10+off)
      if (rs.eq.rsmol(molofrs(rs),1)) then
        ptnr = ptnr + 3
      else if (seqpolty(rs-1).eq.'N') then
!       for non-5'-residue we need the previous methyl fragment (C3* and H3*)
        off2 = 0
        if ((nuci(rs,5).le.0).AND.(nuci(rs,6).le.0).AND.(nuci(rs,4).gt.0)) off2 = 4
        ptlst(4+ptnr) = at(rs-1)%bb(4+off2)
        ptlst(5+ptnr) = at(rs-1)%bb(8+off2)
        ptnr = ptnr + 5
      end if
    else
      if (rs.eq.rsmol(molofrs(rs),1)) then
        ptnr = ptnr + 1
      else if ((seqflag(rs-1).eq.22).OR.(seqflag(rs-1).eq.24)) then
!       for non-5'-residue we need the previous methyl fragment (C3* and H3*)
        off2 = 0
        if (seqflag(rs-1).eq.22) off2 = 4
        ptlst(2+ptnr) = at(rs-1)%bb(4+off2)
        ptnr = ptnr + 2
      end if
    end if
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &N-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! assembly of dummy hydrogen and atom list for C-terminal peptide backbone
! model compounds 
!
subroutine gen_all_Cter(rs,resname,ptnr,ptlst,bbnr,bblst,vs)
!
  use atoms
  use iounit
  use molecule
  use sequen
  use aminos
  use energies
  use polypep
  use fyoc
  use system 
!
  implicit none
!
  integer i,rs,ptlst(20),ptnr,bblst(maxval(at(1:nseq)%nbb+at(1:nseq)%nsc)+50),bbnr,refh
  RTYPE vs(50)
  character(3) resname
!
  ptnr = 0
  bbnr = 0
!
  if (ua_model.gt.0) then
!   do nothing
  else if (resname.eq.'GLY') then 
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-1
    call genxyz(ptlst(1),cai(rs),vs(1),ci(rs),vs(2),&
 &            refh,vs(2),-chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptnr = 1
  else if (resname.eq.'AIB') then
    refh = at(rs)%sc(at(rs)%nsc)
    ptlst(1) = n+20-3
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-2
    call genxyz(ptlst(2),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),-chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptlst(3) = n+20-1
    call genxyz(ptlst(3),cai(rs),vs(1),ci(rs),vs(2),&
 &                ptlst(2),vs(2),-chiral(rs))
    atr(ptlst(3)) = atr(refh)
    atvol(ptlst(3)) = atvol(refh)
    atsavred(ptlst(3)) = atsavred(refh)
    i12(1,ptlst(3)) = cai(rs)
    ptnr = 3
  else if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &  (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &  (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &  (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &  (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &  (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &  (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &  (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &  (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &  (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &  (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &  (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &  (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &  (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &  (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &  (resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &  (resname.eq.'KM2').OR.(resname.eq.'KM3').OR.&
 &  (resname.eq.'PRO').OR.&
 &  (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &  (resname.eq.'HIP')) then
    refh = at(rs)%sc(1)
    ptlst(1) = n+20-2
    call genxyz(ptlst(1),cai(rs),vs(1),ni(rs),vs(2),&
 &            ci(rs),vs(2),chiral(rs))
    atr(ptlst(1)) = atr(refh)
    atvol(ptlst(1)) = atvol(refh)
    atsavred(ptlst(1)) = atsavred(refh)
    i12(1,ptlst(1)) = cai(rs)
    ptlst(2) = n+20-1
    call genxyz(ptlst(2),cai(rs),vs(1),ci(rs),vs(2),&
 &            refh,vs(2),-chiral(rs))
    atr(ptlst(2)) = atr(refh)
    atvol(ptlst(2)) = atvol(refh)
    atsavred(ptlst(2)) = atsavred(refh)
    i12(1,ptlst(2)) = cai(rs)
    ptnr = 2
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
    if (moltermid(molofrs(rs),2).ne.1) then
      write(ilog,*) 'Fatal. Encountered unsupported C-terminus&
 & (RES: ',rs,') in setup_savol().'
      call fexit()
    end if
!   do nothing
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &C-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end if
!
  do i=1,ptnr
    if (ptlst(i).le.n) then
      write(ilog,*) 'Fatal. Abused real atom as holding place &
 &for dummy atom. This is most certainly a bug.'
      write(ilog,*) 'Atom ',ptlst(i),' in res.',rs,'.'
      call fexit()
    end if
    n12(ptlst(i)) = 1
  end do
!
! only the COO-group will be directly covered
  if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &    (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &    (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &    (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &    (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &    (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &    (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &    (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &    (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &    (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &    (resname.eq.'GLY').OR.(resname.eq.'AIB').OR.&
 &    (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &    (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &    (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &    (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &    (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &    (resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &    (resname.eq.'KM2').OR.(resname.eq.'KM3').OR.&
 &    (resname.eq.'PRO').OR.&
 &    (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &    (resname.eq.'HIP')) then
    bblst(1) = ci(rs)
    bblst(2) = at(rs)%bb(4)
    bblst(3) = at(rs)%bb(5)
    if (moltermid(molofrs(rs),2).eq.1) then
      bbnr = 3
    else
      write(ilog,*) 'Fatal. Encountered unsupported C-terminus&
 & (RES: ',rs,') in setup_savol().'
      call fexit()
    end if
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
    if (moltermid(molofrs(rs),2).ne.1) then
      write(ilog,*) 'Fatal. Encountered unsupported C-terminus&
 & (RES: ',rs,') in setup_savol().'
      call fexit()
    end if
!   do nothing
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &C-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end if
! and now let's assemble the full partner atom list
  if ((resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &    (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &    (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &    (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &    (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &    (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &    (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &    (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &    (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &    (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &    (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &    (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &    (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &    (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &    (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &    (resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &    (resname.eq.'KM2').OR.(resname.eq.'KM3').OR.&
 &    (resname.eq.'PRO').OR.&
 &    (resname.eq.'ORN').OR.(resname.eq.'DAB').OR.&
 &    (resname.eq.'HIP')) then
    ptlst(1+ptnr) = cai(rs)
    if (ua_model.eq.0) then
      ptlst(2+ptnr) = at(rs)%sc(1)
      ptnr = ptnr + 2
    else
      ptnr = ptnr + 1
    end if
  else if (resname.eq.'GLY') then
    ptlst(1+ptnr) = cai(rs)
    if (ua_model.eq.0) then
      ptlst(2+ptnr) = at(rs)%sc(1)
      ptlst(3+ptnr) = at(rs)%sc(2)
      ptnr = ptnr + 3
    else
      ptnr = ptnr + 1
    end if
  else if (resname.eq.'AIB') then
    ptlst(1+ptnr) = cai(rs)
    ptnr = ptnr + 1
  else if ((resname.eq.'RPC').OR.(resname.eq.'DPC').OR.&
 &    (resname.eq.'RPU').OR.(resname.eq.'DPU').OR.&
 &    (resname.eq.'RPT').OR.(resname.eq.'DPT').OR.&
 &    (resname.eq.'RPA').OR.(resname.eq.'DPA').OR.&
 &    (resname.eq.'RPG').OR.(resname.eq.'DPG').OR.&
 &    (resname.eq.'R5P').OR.(resname.eq.'D5P')) then
    if (moltermid(molofrs(rs),2).ne.1) then
      write(ilog,*) 'Fatal. Encountered unsupported C-terminus&
 & (RES: ',rs,') in setup_savol().'
      call fexit()
    end if
!   do nothing
  else
    ptnr = 0
    if (use_IMPSOLV.EQV..true.) then
      write(ilog,*) 'Fatal. Implicit solvent does not support &
 &C-terminal residue ',resname,' at the moment. Check back later.'
      call fexit()
    end if
  end if
!
end
!
!----------------------------------------------------------------------
