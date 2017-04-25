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
!-----------------------------------------------------------------------
!
! "setup_srinter" generates interaction lists on a per-residue basis
! for active non-bonded interactions, this includes intra-residue
! and nearest neighbor contributions (mapped to first of the two residues
! in sequence). all pairs over residues with a seq. separation of at least 2
! are assumed to be contributing completely (i.e., this subroutine does
! not support the case in which one desires to exclude interactions spanning
! more than the neighboring residues!!!)
! the strategy of using equivalence atoms (eqatm) to determine quasi-
! rigid relationships is explained elsewhere
! for bonded interactions the lists are simply set up over residues and
! do not worry about sequence separation, but use 12/13/14 arrays instead
!
!     
subroutine setup_srinter()
!
  use atoms
  use inter
  use aminos
  use iounit
  use params
  use polypep
  use sequen
  use molecule
  use energies
  use system
  use math
  use fyoc
  use zmatrix
!
  implicit none
!
  integer i,j,k,rs,ii,kk,iii,jjj,jj,q,qq,shf
  integer rsi,rsj,l,m,ll,ils(4),bbc
  integer, ALLOCATABLE:: tempa(:,:),tempb(:,:),tempc(:,:)
  integer, ALLOCATABLE:: tempd(:,:),tempm(:,:)
  integer imol,skipper,alcsz
  logical proceed,its14,delay_report
  RTYPE getblen,getbang,getztor,getblen_ref,getbang_ref,getztor_ref
!
  delay_report = .false.
  do i=1,n
    if (eqatm(i).lt.0) eqatm(i) = i
  end do
!
  skipper= MAXVLNC + MAXVLNC*(MAXVLNC-1) + MAXVLNC*(MAXVLNC-1)*(MAXVLNC-1)
  allocate(fudge(nseq))
  allocate(iaa(nseq))
!
  do imol=1,nmol
! 
    do rs=rsmol(imol,1),rsmol(imol,2)
!
      if (natres(rs).gt.1) then
        allocate(fudge(rs)%rsin((natres(rs)*(natres(rs)-1))/2))
        allocate(fudge(rs)%rsin_ljs((natres(rs)*(natres(rs)-1))/2))
        allocate(fudge(rs)%rsin_lje((natres(rs)*(natres(rs)-1))/2))
        allocate(iaa(rs)%atin((natres(rs)*(natres(rs)-1))/2,2))
        do k=1,(natres(rs)*(natres(rs)-1))/2
          fudge(rs)%rsin(k) = 1.0
        end do
      end if
      if (rs.lt.nseq) then
        alcsz = natres(rs)*natres(rs+1)
        if (disulf(rs).gt.0) then
          alcsz = max(alcsz,natres(rs)*natres(disulf(rs)))
        end if
        allocate(fudge(rs)%rsnb(alcsz))
        allocate(fudge(rs)%rsnb_ljs(alcsz))
        allocate(fudge(rs)%rsnb_lje(alcsz))
        allocate(iaa(rs)%atnb(alcsz,2))
        do k=1,alcsz
          fudge(rs)%rsnb(k) = 1.0
        end do
      end if
!
!     get relevant intra-residue interactions for fixed distances,angles
!
      nrsintra(rs) = 0
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
!       dummy atoms are still in those arrays! (dangerous ...)
!       skip mercilessly
        if (attyp(ii) .le. 0)  cycle
        do k=i+1,at(rs)%nbb+at(rs)%nsc
          if (k.le.at(rs)%nbb) then
            kk = at(rs)%bb(k)
          else
            kk = at(rs)%sc(k-at(rs)%nbb)
          end if
          if (attyp(kk) .le. 0)  cycle
          proceed = .true.
          its14 = .false.
!
!         in Cartesian space, connectivity is the only thing that matters (even though
!         a lot of interactions will show very little fluctuations)
          if ((nbsr_model.eq.2).OR.(nbsr_model.eq.3)) then
            do j=1,n12(ii)
              if (i12(j,ii).eq.kk) proceed = .false.
            end do
            do j=1,n13(ii)
              if (i13(j,ii).eq.kk) proceed = .false.
            end do
            do j=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(j,ii).eq.kk) proceed = .false.
              else
                if (i14(j,ii).eq.kk) its14 = .true.
              end if
            end do
            if (proceed.EQV..true.) then
              nrsintra(rs) = nrsintra(rs) + 1 
              iaa(rs)%atin(nrsintra(rs),1) = ii
              iaa(rs)%atin(nrsintra(rs),2) = kk
              fudge(rs)%rsin_ljs(nrsintra(rs)) = lj_sig(attyp(ii),attyp(kk))
              fudge(rs)%rsin_lje(nrsintra(rs)) = lj_eps(attyp(ii),attyp(kk))
              if (its14.EQV..true.) then
                fudge(rs)%rsin(nrsintra(rs)) = fudge_st_14
                fudge(rs)%rsin_ljs(nrsintra(rs)) = lj_sig_14(attyp(ii),attyp(kk))
                fudge(rs)%rsin_lje(nrsintra(rs)) = lj_eps_14(attyp(ii),attyp(kk))
              end if
            end if
!
!         with hard-coded rigid constraints, a lot more interactions can be excluded
!         however, setup work also becomes equally more complicated
          else if (nbsr_model.eq.1) then
            if (seqtyp(rs).eq.26) then
              delay_report = .true.
              cycle ! handled separately
            end if
            if ((eqatm(ii).eq.eqatm(kk)).OR.&
 &   (eqatm(ii).eq.kk).OR.(eqatm(kk).eq.ii)) proceed = .false.
!
            do j=1,n12(ii)
              if (i12(j,ii).eq.kk) proceed = .false.
              if (i12(j,ii).eq.eqatm(kk)) proceed = .false.
            end do
            if (eqatm(ii).ne.ii) then
              do j=1,n12(eqatm(ii))
                if (i12(j,eqatm(ii)).eq.kk) proceed = .false.
                if (i12(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
              end do
            end if
!
            do j=1,n13(ii)
              if (i13(j,ii).eq.kk) proceed = .false.
              if (i13(j,ii).eq.eqatm(kk)) proceed = .false.
            end do
            if (eqatm(ii).ne.ii) then
              do j=1,n13(eqatm(ii))
                if (i13(j,eqatm(ii)).eq.kk) proceed = .false.
                if (i13(j,eqatm(ii)).eq.eqatm(kk)) proceed = .false.
              end do
            end if
!
            do j=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(j,ii).eq.kk) proceed = .false.
                if (i14(j,ii).eq.eqatm(kk)) proceed = .false.
              else
                if (i14(j,ii).eq.kk) its14 = .true.
                if (mode_14.eq.2) then
                  if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                end if
              end if
            end do
            if (eqatm(ii).ne.ii) then
              do j=1,n14(eqatm(ii))
                if (use_14.EQV..false.) then
                  if (i14(j,eqatm(ii)).eq.kk) proceed = .false.
                  if (i14(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                else
                  if (mode_14.eq.2) then
                    if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                    if(i14(j,eqatm(ii)).eq.eqatm(kk)) its14 = .true.
                   end if
                 end if
               end do
            end if
!
            if (proceed.EQV..true.) then
              nrsintra(rs) = nrsintra(rs) + 1 
              iaa(rs)%atin(nrsintra(rs),1) = ii
              iaa(rs)%atin(nrsintra(rs),2) = kk
              fudge(rs)%rsin_ljs(nrsintra(rs)) = lj_sig(attyp(ii),attyp(kk))
              fudge(rs)%rsin_lje(nrsintra(rs)) = lj_eps(attyp(ii),attyp(kk))
              if (its14.EQV..true.) then
                fudge(rs)%rsin(nrsintra(rs)) = fudge_st_14
                fudge(rs)%rsin_ljs(nrsintra(rs)) = lj_sig_14(attyp(ii),attyp(kk))
                fudge(rs)%rsin_lje(nrsintra(rs)) = lj_eps_14(attyp(ii),attyp(kk))
              end if
            end if
!
          end if ! which value of nbsr_model
!
        end do
      end do
!
!     cycle out if last residue in sequence
      nrsnb(rs) = 0
      if (rs.eq.nseq) cycle
!
!     get relevant next-residue in chain contacts for fixed distances/angles
!
!     this will work whether the residue is the last residue in a molecule or not
!     (checked inside). in the former case, all possible interaction to the next
!     residue are of course added
!
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
!       dummy atoms might still be in those arrays! (dangerous ...)
!       skip mercilessly
        if (attyp(ii) .le. 0)  cycle
        do k=1,at(rs+1)%nbb+at(rs+1)%nsc
          if (k.le.at(rs+1)%nbb) then
            kk = at(rs+1)%bb(k)
          else
            kk = at(rs+1)%sc(k-at(rs+1)%nbb)
          end if
          if (attyp(kk) .le. 0)  cycle
          proceed = .true.
          its14 = .false.
!
!         Cartesian space, see above
          if ((nbsr_model.eq.2).OR.(nbsr_model.eq.3)) then
            do j=1,n12(ii)
              if (i12(j,ii).eq.kk) proceed = .false.
            end do
            do j=1,n13(ii)
              if (i13(j,ii).eq.kk) proceed = .false.
            end do
            do j=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(j,ii).eq.kk) proceed = .false.
              else
                if (i14(j,ii).eq.kk) its14 = .true.
              end if
            end do
            if (proceed.EQV..true.) then
              nrsnb(rs) = nrsnb(rs) + 1 
              iaa(rs)%atnb(nrsnb(rs),1) = ii
              iaa(rs)%atnb(nrsnb(rs),2) = kk
              fudge(rs)%rsnb_ljs(nrsnb(rs)) = lj_sig(attyp(ii),attyp(kk))
              fudge(rs)%rsnb_lje(nrsnb(rs)) = lj_eps(attyp(ii),attyp(kk))
              if (its14.EQV..true.) then
                fudge(rs)%rsnb(nrsnb(rs)) = fudge_st_14
                fudge(rs)%rsnb_ljs(nrsnb(rs)) = lj_sig_14(attyp(ii),attyp(kk))
                fudge(rs)%rsnb_lje(nrsnb(rs)) = lj_eps_14(attyp(ii),attyp(kk))
              end if
            end if
!
!         with rigid constraints, see above
          else if (nbsr_model.eq.1) then
            if (((seqtyp(rs).eq.26).OR.(seqtyp(rs+1).eq.26)).AND.(molofrs(rs).eq.molofrs(rs+1))) then
              delay_report = .true.
              cycle ! handled separately
            end if
            if ((eqatm(ii).eq.eqatm(kk)).OR.(eqatm(ii).eq.kk).OR.(eqatm(kk).eq.ii)) proceed = .false.
!
            do j=1,n12(ii)
              if (i12(j,ii).eq.kk) proceed = .false.
              if (i12(j,ii).eq.eqatm(kk)) proceed = .false.
            end do
            if (eqatm(ii).ne.ii) then
              do j=1,n12(eqatm(ii))
                if (i12(j,eqatm(ii)).eq.kk) proceed = .false.
                if (i12(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
              end do
            end if
!
            do j=1,n13(ii)
              if (i13(j,ii).eq.kk) proceed = .false.
              if (i13(j,ii).eq.eqatm(kk)) proceed = .false.
            end do
            if (eqatm(ii).ne.ii) then
              do j=1,n13(eqatm(ii))
                if (i13(j,eqatm(ii)).eq.kk) proceed = .false.
                if (i13(j,eqatm(ii)).eq.eqatm(kk)) proceed= .false.
              end do
            end if
!
            do j=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(j,ii).eq.kk) proceed = .false.
                if (i14(j,ii).eq.eqatm(kk)) proceed = .false.
              else
                if (i14(j,ii).eq.kk) its14 = .true.
                if (mode_14.eq.2) then
                  if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                end if
              end if
            end do
!
            if (eqatm(ii).ne.ii) then
              do j=1,n14(eqatm(ii))
                if (use_14.EQV..false.) then
                  if (i14(j,eqatm(ii)).eq.kk) proceed = .false.
                  if (i14(j,eqatm(ii)).eq.eqatm(kk)) proceed=.false.
                else
                  if (mode_14.eq.2) then
                    if (i14(j,ii).eq.eqatm(kk)) its14 = .true.
                    if (i14(j,eqatm(ii)).eq.eqatm(kk)) its14= .true.
                  end if
                end if
              end do
            end if
!
            if (proceed.EQV..true.) then
              nrsnb(rs) = nrsnb(rs) + 1 
              iaa(rs)%atnb(nrsnb(rs),1) = ii
              iaa(rs)%atnb(nrsnb(rs),2) = kk
              fudge(rs)%rsnb_ljs(nrsnb(rs)) = lj_sig(attyp(ii),attyp(kk))
              fudge(rs)%rsnb_lje(nrsnb(rs)) = lj_eps(attyp(ii),attyp(kk))
              if (its14.EQV..true.) then
                fudge(rs)%rsnb(nrsnb(rs)) = fudge_st_14
                fudge(rs)%rsnb_ljs(nrsnb(rs)) = lj_sig_14(attyp(ii),attyp(kk))
                fudge(rs)%rsnb_lje(nrsnb(rs)) = lj_eps_14(attyp(ii),attyp(kk))
              end if
            end if
! 
          end if ! which value of nbsr_model
!
        end do
      end do
!
    end do !loop over residues
!
  end do !loop over molecules
!
  if (nbsr_model.eq.3) then
    call GROMOS_excludes()
  end if
!
!
! ok, now on to the bonded interactions
!
  shf = 0
  if (ua_model.gt.0) shf = 1
!
  do rs=1,nseq
!
    nrsbl(rs) = 0
    nrsbleff(rs) = 0
    nrsba(rs) = 0
    nrsbaeff(rs) = 0
    nrsdi(rs) = 0
    nrscm(rs) = 0
    nrsdieff(rs) = 0
    nrsimpt(rs) = 0
    nrsimpteff(rs) = 0
    nrscmeff(rs) = 0
!
    skipper = natres(rs)*MAXVLNC
    allocate(tempa(skipper,2))
    skipper = natres(rs)*MAXVLNC*(MAXVLNC-1)
    allocate(tempb(skipper,3))
    skipper = natres(rs)*MAXVLNC*(MAXVLNC-1)*(MAXVLNC-1)
    allocate(tempc(skipper,4))
    allocate(tempd(skipper,4))
    allocate(tempm(skipper,5))
!
    do i=1,at(rs)%nbb+at(rs)%nsc
      if (i.le.at(rs)%nbb) then
        ii = at(rs)%bb(i)
      else
        ii = at(rs)%sc(i-at(rs)%nbb)
      end if
!
!     still make sure to cycle away dummy atoms
      if (attyp(ii).le.0)  cycle
!
!     bonds
      do j=1,n12(ii)
        jj = i12(j,ii)
        proceed = .true.
        rsj = atmres(jj)
        rsj = min(rsj,rs)
        do rsi=rsj,rs-1
          do k=1,nrsbl(rsi)
            if (&
 &    ((iaa(rsi)%bl(k,1).eq.ii).AND.(iaa(rsi)%bl(k,2).eq.jj))&
 &.OR.((iaa(rsi)%bl(k,1).eq.jj).AND.(iaa(rsi)%bl(k,2).eq.ii))) then
              proceed = .false.
              exit
            end if
          end do
        end do
        if (proceed.EQV..true.) then
          do k=1,nrsbl(rs)
            if (&
 &    ((tempa(k,1).eq.ii).AND.(tempa(k,2).eq.jj))&
 &.OR.((tempa(k,1).eq.jj).AND.(tempa(k,2).eq.ii))) then
              proceed = .false.
              exit
            end if
          end do
        end if
        if (proceed.EQV..false.) cycle
        nrsbl(rs) = nrsbl(rs) + 1
        tempa(nrsbl(rs),1) = ii
        tempa(nrsbl(rs),2) = jj
      end do
!
!     angles
      do j=1,n13(ii)
        jj = i13(j,ii)
!       find the middle atom
        do k=1,n12(ii)
          do l=1,n12(jj)
            if (i12(k,ii).eq.i12(l,jj)) then
!             kk is our current candidate to complet the set ii,jj,kk
              kk = i12(k,ii)
              proceed = .true.
              rsj = min(atmres(jj),atmres(kk),rs)
              do rsi=rsj,rs-1
                do m=1,nrsba(rsi)
                  do iii=1,3
                    if ((iaa(rsi)%ba(m,iii).ne.ii).AND.&
 &                      (iaa(rsi)%ba(m,iii).ne.jj).AND.&
 &                      (iaa(rsi)%ba(m,iii).ne.kk)) then
                     proceed = .true.
                      exit
                    else
                      proceed = .false.
                    end if
                  end do
                  if (proceed.EQV..false.) exit
                end do
                if (proceed.EQV..false.) exit
              end do
              if (proceed.EQV..true.) then
                do m=1,nrsba(rs)
                  do iii=1,3
                    if ((tempb(m,iii).ne.ii).AND.&
 &                      (tempb(m,iii).ne.jj).AND.&
 &                      (tempb(m,iii).ne.kk)) then
                      proceed = .true.
                      exit
                    else
                      proceed = .false.
                    end if
                  end do
                  if (proceed.EQV..false.) exit
                end do
              end if
              if (proceed.EQV..true.) then
                nrsba(rs) = nrsba(rs) + 1
                tempb(nrsba(rs),1) = ii
                tempb(nrsba(rs),2) = kk
                tempb(nrsba(rs),3) = jj
              end if
            end if
          end do
        end do
      end do
!
!     improper dihedrals (typically planar trigonal centers or chiral tetrahedral centers)
      if ((n12(ii).eq.3).OR.(n12(ii).eq.4)) then
        if (n12(ii).eq.3) then
!         heaviest atom (mass) defines default last atom (only relevant for cases with identity)
!         if two heavier substituents are of same mass, third defaults to last (unique case)
!         if all atoms are equivalent, order is "as is" (see permutations below)
          nrsimpt(rs) = nrsimpt(rs) + 1
          tempd(nrsimpt(rs),1) = ii
          if (mass(i12(1,ii)).gt.mass(i12(2,ii))) then
            if (mass(i12(1,ii)).gt.mass(i12(3,ii))) then
              tempd(nrsimpt(rs),4) = i12(1,ii)
              tempd(nrsimpt(rs),2) = i12(2,ii)
              tempd(nrsimpt(rs),3) = i12(3,ii)
            else if (mass(i12(1,ii)).eq.mass(i12(3,ii))) then
              tempd(nrsimpt(rs),4) = i12(2,ii)
              tempd(nrsimpt(rs),2) = i12(1,ii)
              tempd(nrsimpt(rs),3) = i12(3,ii)
            else
              tempd(nrsimpt(rs),4) = i12(3,ii)
              tempd(nrsimpt(rs),2) = i12(1,ii)
              tempd(nrsimpt(rs),3) = i12(2,ii)
            end if
          else if (mass(i12(1,ii)).eq.mass(i12(2,ii))) then
            if (mass(i12(1,ii)).eq.mass(i12(3,ii))) then
              tempd(nrsimpt(rs),4) = i12(1,ii)
              tempd(nrsimpt(rs),2) = i12(2,ii)
              tempd(nrsimpt(rs),3) = i12(3,ii)
            else
              tempd(nrsimpt(rs),4) = i12(3,ii)
              tempd(nrsimpt(rs),2) = i12(1,ii)
              tempd(nrsimpt(rs),3) = i12(2,ii)
            end if
          else
            if (mass(i12(2,ii)).gt.mass(i12(3,ii))) then
              tempd(nrsimpt(rs),4) = i12(2,ii)
              tempd(nrsimpt(rs),2) = i12(1,ii)
              tempd(nrsimpt(rs),3) = i12(3,ii)
            else if (mass(i12(2,ii)).eq.mass(i12(3,ii))) then
              tempd(nrsimpt(rs),4) = i12(1,ii)
              tempd(nrsimpt(rs),2) = i12(2,ii)
              tempd(nrsimpt(rs),3) = i12(3,ii)
            else
              tempd(nrsimpt(rs),4) = i12(3,ii)
              tempd(nrsimpt(rs),2) = i12(1,ii)
              tempd(nrsimpt(rs),3) = i12(2,ii)
            end if
          end if
        else if (n12(ii).eq.4) then
!         check if any one of the four atoms is completely transparent
          iii = 0
          do j=1,n12(ii)
            if (bio_ctyp(b_type(i12(j,ii))).gt.0) then
              if ((lj_eps(bio_ljtyp(b_type(i12(j,ii))),bio_ljtyp(b_type(i12(j,ii)))).ne.0.0).OR.&
 &              (c_charge(bio_ctyp(b_type(i12(j,ii)))).ne.0.0)) then
                iii = iii + 1
                ils(iii) = i12(j,ii)
              end if
            else
              if (lj_eps(bio_ljtyp(b_type(i12(j,ii))),bio_ljtyp(b_type(i12(j,ii)))).ne.0.0) then
                iii = iii + 1
                ils(iii) = i12(j,ii)
              end if
            end if
          end do
          if (iii.eq.3) then
            nrsimpt(rs) = nrsimpt(rs) + 1
            tempd(nrsimpt(rs),1) = ii
            if (mass(ils(1)).gt.mass(ils(2))) then
              if (mass(ils(1)).gt.mass(ils(3))) then
                tempd(nrsimpt(rs),4) = ils(1)
                tempd(nrsimpt(rs),2) = ils(2)
                tempd(nrsimpt(rs),3) = ils(3)
              else if (mass(ils(1)).eq.mass(ils(3))) then
                tempd(nrsimpt(rs),4) = ils(2)
                tempd(nrsimpt(rs),2) = ils(1)
                tempd(nrsimpt(rs),3) = ils(3)
              else
                tempd(nrsimpt(rs),4) = ils(3)
                tempd(nrsimpt(rs),2) = ils(1)
                tempd(nrsimpt(rs),3) = ils(2)
              end if
            else if (mass(ils(1)).eq.mass(ils(2))) then
              if (mass(ils(1)).eq.mass(ils(3))) then
                tempd(nrsimpt(rs),4) = ils(1)
                tempd(nrsimpt(rs),2) = ils(2)
                tempd(nrsimpt(rs),3) = ils(3)
              else
                tempd(nrsimpt(rs),4) = ils(3)
                tempd(nrsimpt(rs),2) = ils(1)
                tempd(nrsimpt(rs),3) = ils(2)
              end if
            else
              if (mass(ils(2)).gt.mass(ils(3))) then
                tempd(nrsimpt(rs),4) = ils(2)
                tempd(nrsimpt(rs),2) = ils(1)
                tempd(nrsimpt(rs),3) = ils(3)
              else if (mass(ils(2)).eq.mass(ils(3))) then
                tempd(nrsimpt(rs),4) = ils(1)
                tempd(nrsimpt(rs),2) = ils(2)
                tempd(nrsimpt(rs),3) = ils(3)
              else
                tempd(nrsimpt(rs),4) = ils(3)
                tempd(nrsimpt(rs),2) = ils(1)
                tempd(nrsimpt(rs),3) = ils(2)
              end if
            end if
          end if
!
        end if
      end if
!
!     dihedrals
      do j=1,n12(ii)
        jj = i12(j,ii)
!       find a non-reversing chain of connected atoms
        do k=1,n12(jj)
          kk = i12(k,jj)
          do l=1,n12(kk)
            ll = i12(l,kk)
            if ((ll.eq.ii).OR.(ll.eq.jj).OR.&
 &              (kk.eq.ii)) cycle
            proceed = .true.
!           kk and ll are our current axis candidates to complete the set ii,jj,kk,ll
            rsj = min(atmres(jj),atmres(kk),atmres(ll),rs)
            do rsi=rsj,rs-1
              do jjj=1,nrsdi(rsi)
                do iii=1,4
                  if ((iaa(rsi)%di(jjj,iii).ne.ii).AND.&
 &                    (iaa(rsi)%di(jjj,iii).ne.jj).AND.&
 &                    (iaa(rsi)%di(jjj,iii).ne.kk).AND.&
 &                    (iaa(rsi)%di(jjj,iii).ne.ll)) then
                    proceed = .true.
                    exit
                  else
                    proceed = .false.
                  end if
                end do
                if (proceed.EQV..false.) exit
              end do
              if (proceed.EQV..false.) exit
            end do
            if (proceed.EQV..true.) then
              do jjj=1,nrsdi(rs)
                do iii=1,4
                  if ((tempc(jjj,iii).ne.ii).AND.&
 &                    (tempc(jjj,iii).ne.jj).AND.&
 &                    (tempc(jjj,iii).ne.kk).AND.&
 &                    (tempc(jjj,iii).ne.ll)) then
                    proceed = .true.
                    exit
                  else
                    proceed = .false.
                  end if
                end do
                if (proceed.EQV..false.) exit
              end do
            end if
            if (proceed.EQV..true.) then
              nrsdi(rs) = nrsdi(rs) + 1
              tempc(nrsdi(rs),1) = ii
              tempc(nrsdi(rs),2) = jj
              tempc(nrsdi(rs),3) = kk
              tempc(nrsdi(rs),4) = ll
            end if
          end do
        end do
      end do
!
!     last and least: CMAP (two immediately connected all-backbone dihedrals)
      if (i.gt.at(rs)%nbb) cycle
      if (n12(ii).eq.1) cycle
      do j=1,n12(ii)
        jj = i12(j,ii)
        proceed = .false.
        do bbc=1,at(atmres(jj))%nbb
          if (at(atmres(jj))%bb(bbc).eq.jj) proceed = .true.
        end do
        if (proceed.EQV..false.) cycle
!       find a non-reversing chain of connected atoms
        do k=1,n12(jj)
          kk = i12(k,jj)
          if (kk.eq.ii) cycle
!          if (atmres(kk).ne.rs) cycle
          proceed = .false.
          do bbc=1,at(atmres(kk))%nbb
            if (at(atmres(kk))%bb(bbc).eq.kk) proceed = .true.
          end do
          if (proceed.EQV..false.) cycle
          do l=1,n12(kk)
            ll = i12(l,kk)
            if ((ll.eq.ii).OR.(ll.eq.jj)) cycle
            proceed = .false.
            do bbc=1,at(atmres(ll))%nbb
              if (at(atmres(ll))%bb(bbc).eq.ll) proceed = .true.
            end do
            if (proceed.EQV..false.) cycle
            do q=1,n12(ll)
              qq = i12(q,ll)
              if (n12(qq).eq.1) cycle
              if ((qq.eq.ii).OR.(qq.eq.jj).OR.(qq.eq.kk)) cycle
              proceed = .false.
              do bbc=1,at(atmres(qq))%nbb
                if (at(atmres(qq))%bb(bbc).eq.qq) proceed = .true.
              end do
              if (proceed.EQV..false.) cycle
              rsj = min(atmres(jj),atmres(kk),atmres(ll),atmres(qq),rs)
              do rsi=rsj,rs-1
                do jjj=1,nrscm(rsi)
                  do iii=1,5
                    if ((iaa(rsi)%cm(jjj,iii).ne.ii).AND.&
 &                      (iaa(rsi)%cm(jjj,iii).ne.jj).AND.&
 &                      (iaa(rsi)%cm(jjj,iii).ne.kk).AND.&
 &                      (iaa(rsi)%cm(jjj,iii).ne.ll).AND.&
 &                      (iaa(rsi)%cm(jjj,iii).ne.qq)) then
                      proceed = .true.
                      exit
                    else
                      proceed = .false.
                    end if
                  end do
                  if (proceed.EQV..false.) exit
                end do
                if (proceed.EQV..false.) exit
              end do
              if (proceed.EQV..true.) then
                do jjj=1,nrscm(rs)
                  do iii=1,5
                    if ((tempm(jjj,iii).ne.ii).AND.&
 &                      (tempm(jjj,iii).ne.jj).AND.&
 &                      (tempm(jjj,iii).ne.kk).AND.&
 &                      (tempm(jjj,iii).ne.ll).AND.&
 &                      (tempm(jjj,iii).ne.qq)) then
                      proceed = .true.
                      exit
                    else
                      proceed = .false.
                    end if
                  end do
                  if (proceed.EQV..false.) exit
                end do
              end if
              if (proceed.EQV..true.) then
                nrscm(rs) = nrscm(rs) + 1
                tempm(nrscm(rs),1) = ii
                tempm(nrscm(rs),2) = jj
                tempm(nrscm(rs),3) = kk
                tempm(nrscm(rs),4) = ll
                tempm(nrscm(rs),5) = qq
                nrscm(rs) = nrscm(rs) + 1
                tempm(nrscm(rs),1) = qq
                tempm(nrscm(rs),2) = ll
                tempm(nrscm(rs),3) = kk
                tempm(nrscm(rs),4) = jj
                tempm(nrscm(rs),5) = ii
              end if
            end do
          end do
        end do
      end do
!
    end do
!
!   now allocate and transfer
  77    format('biotype_bond ',3(1x,i5))
!    write(*,*) rs,':'
    allocate(iaa(rs)%bl(nrsbl(rs),2))
    allocate(iaa(rs)%par_bl(nrsbl(rs),MAXBOPAR+1))
    allocate(iaa(rs)%typ_bl(nrsbl(rs)))
    iaa(rs)%typ_bl(:) = 0
    iaa(rs)%par_bl(:,:) = 0.
    do i=1,nrsbl(rs)
      do j=1,2
        iaa(rs)%bl(i,j) = tempa(i,j)
      end do
!     we have to skip those bonded terms coming form crosslinks
      proceed = .true.
      do k=1,n_crosslinks
        if (((atmres(tempa(i,1)).eq.crosslink(k)%rsnrs(1)).AND.(atmres(tempa(i,2)).eq.crosslink(k)%rsnrs(2))).OR.&
 &          ((atmres(tempa(i,1)).eq.crosslink(k)%rsnrs(2)).AND.(atmres(tempa(i,2)).eq.crosslink(k)%rsnrs(1)))) then
          if (crosslink(k)%itstype.le.2) then
            if ((at(atmres(tempa(i,1)))%sc(3-shf).eq.tempa(i,1)).AND.(at(atmres(tempa(i,2)))%sc(3-shf).eq.tempa(i,2))) then
              iaa(rs)%par_bl(i,MAXBOPAR+1)  = 2.03
              proceed = .false.
              exit
            end if
          else
            write(ilog,*) 'Fatal. Encountered unsupported crosslink type in setup_srinter(...). This is an omission bug.'
            call fexit()
          end if
        end if
      end do
      if (proceed.EQV..false.) cycle
      if ((atmres(tempa(i,1)).ne.atmres(tempa(i,2))).AND.&
 &    ((seqtyp(atmres(tempa(i,1))).eq.26).OR.(seqtyp(atmres(tempa(i,2))).eq.26))) then
        iaa(rs)%par_bl(i,2) = getblen_ref(tempa(i,1),tempa(i,2))
      else
        iaa(rs)%par_bl(i,2) = getblen(tempa(i,1),tempa(i,2))
      end if
      iaa(rs)%par_bl(i,MAXBOPAR+1) = iaa(rs)%par_bl(i,2)
      if (guess_bonded.ge.1) then
        iaa(rs)%par_bl(i,1) = 300.0
        iaa(rs)%typ_bl(i) = 1
      end if
    end do
    allocate(iaa(rs)%ba(nrsba(rs),3))
    allocate(iaa(rs)%par_ba(nrsba(rs),MAXBAPAR+1))
    allocate(iaa(rs)%typ_ba(nrsba(rs)))
    iaa(rs)%typ_ba(:) = 0
    iaa(rs)%par_ba(:,:) = 0.
    do i=1,nrsba(rs)
      do j=1,3
        iaa(rs)%ba(i,j) = tempb(i,j)
      end do
!     we have to skip those bonded terms coming form crosslinks
      proceed = .true.
      do k=1,n_crosslinks
        if (((maxval(atmres(tempb(i,1:3))).eq.crosslink(k)%rsnrs(1)).AND.(minval(atmres(tempb(i,1:3))).eq.crosslink(k)%rsnrs(2)))&
 &.OR.((maxval(atmres(tempb(i,1:3))).eq.crosslink(k)%rsnrs(2)).AND.(minval(atmres(tempb(i,1:3))).eq.crosslink(k)%rsnrs(1)))) then
          if (crosslink(k)%itstype.le.2) then
            if (((at(atmres(tempb(i,1)))%sc(2-shf).eq.tempb(i,1)).AND.(at(atmres(tempb(i,2)))%sc(3-shf).eq.tempb(i,2)).AND.&
 &               (at(atmres(tempb(i,3)))%sc(3-shf).eq.tempb(i,3))).OR.&
 &              ((at(atmres(tempb(i,1)))%sc(3-shf).eq.tempb(i,1)).AND.(at(atmres(tempb(i,2)))%sc(3-shf).eq.tempb(i,2)).AND.&
 &               (at(atmres(tempb(i,3)))%sc(2-shf).eq.tempb(i,3)))) then
              iaa(rs)%par_ba(i,MAXBAPAR+1) = 103.0
              proceed = .false.
              exit
            end if
          else
            write(ilog,*) 'Fatal. Encountered unsupported crosslink type in setup_srinter(...). This is an omission bug.'
            call fexit()
          end if

            proceed = .false.
            exit
        end if
      end do
      if (proceed.EQV..false.) cycle
      iaa(rs)%par_ba(i,1) = 80.0/(RADIAN*RADIAN)
      if ((maxval(atmres(tempb(i,1:3))).ne.minval(atmres(tempb(i,1:3)))).AND.&
 &  ((seqtyp(atmres(tempb(i,1))).eq.26).OR.(seqtyp(atmres(tempb(i,2))).eq.26).OR.(seqtyp(atmres(tempb(i,3))).eq.26))) then
        iaa(rs)%par_ba(i,2) = getbang_ref(tempb(i,1),tempb(i,2),tempb(i,3))
      else
        iaa(rs)%par_ba(i,2) = getbang(tempb(i,1),tempb(i,2),tempb(i,3))
      end if
      iaa(rs)%par_ba(i,MAXBAPAR+1) = iaa(rs)%par_ba(i,2)
      if (guess_bonded.ge.1) then
        iaa(rs)%typ_ba(i) = 1
        iaa(rs)%par_ba(i,1) = 80.0/(RADIAN*RADIAN)
      end if
    end do
    allocate(iaa(rs)%di(nrsdi(rs),4))
    allocate(iaa(rs)%par_di(nrsdi(rs),MAXDIPAR))
    allocate(iaa(rs)%typ_di(nrsdi(rs)))
    iaa(rs)%typ_di(:) = 0
!    write(*,*) 'dihedrals'
    do i=1,nrsdi(rs)
      do j=1,4
        iaa(rs)%di(i,j) = tempc(i,j)
      end do
    end do
!
!   we will triplicate impropers to allow for the various possible definitions of the reference plane
!   holding the central atom
    allocate(iaa(rs)%impt(3*nrsimpt(rs),4))
    allocate(iaa(rs)%par_impt(3*nrsimpt(rs),MAXDIPAR))
    allocate(iaa(rs)%typ_impt(3*nrsimpt(rs)))
    iaa(rs)%typ_impt(:) = 0
    iaa(rs)%par_impt(:,:) = 0.0
!    write(*,*) 'improper dihedrals'
!   note that in case of a single degeneracy (two identical biotypes), the unique atom
!   is in tempd(i,4): it is important that we choose a set of permutations which places it in
!   each position once -> remember that the redundant permutation for planar impropers is flipping 2 and 3
!   (so alternatively we could do 1324, 1423, and 1342 for those)
!   otherwise, the comparisons via biotype to the parameters fail later on (incomplete)
    do i=1,nrsimpt(rs)
      iaa(rs)%impt(3*i-2,1) = tempd(i,1)
      iaa(rs)%impt(3*i-2,2) = tempd(i,2)
      iaa(rs)%impt(3*i-2,3) = tempd(i,3)
      iaa(rs)%impt(3*i-2,4) = tempd(i,4)
!     2nd permutation
      iaa(rs)%impt(3*i-1,1) = tempd(i,1)
      iaa(rs)%impt(3*i-1,2) = tempd(i,2)
      iaa(rs)%impt(3*i-1,3) = tempd(i,4)
      iaa(rs)%impt(3*i-1,4) = tempd(i,3)
!     3rd permutation
      iaa(rs)%impt(3*i,1) = tempd(i,1)
      iaa(rs)%impt(3*i,2) = tempd(i,4)
      iaa(rs)%impt(3*i,3) = tempd(i,3)
      iaa(rs)%impt(3*i,4) = tempd(i,2)
!      write(*,*) iaa(rs)%impt(i,:)
    end do
    nrsimpt(rs) = nrsimpt(rs)*3
    if (guess_bonded.ge.2) then
      do i=1,nrsimpt(rs)
        do k=1,n_crosslinks
         if (((maxval(atmres(iaa(rs)%impt(i,1:4))).eq.crosslink(k)%rsnrs(1)).AND.&
 &            (minval(atmres(iaa(rs)%impt(i,1:4))).eq.crosslink(k)%rsnrs(2))).OR.&
 &           ((maxval(atmres(iaa(rs)%impt(i,1:4))).eq.crosslink(k)%rsnrs(2)).AND.&
 &            (minval(atmres(iaa(rs)%impt(i,1:4))).eq.crosslink(k)%rsnrs(1)))) then
            if (crosslink(k)%itstype.le.2) then
              if ((at(atmres(iaa(rs)%impt(i,1)))%sc(3-shf).eq.iaa(rs)%impt(i,1)).OR.&
 &                (at(atmres(iaa(rs)%impt(i,2)))%sc(3-shf).eq.iaa(rs)%impt(i,2)).OR.&
 &                (at(atmres(iaa(rs)%impt(i,3)))%sc(3-shf).eq.iaa(rs)%impt(i,3)).OR.&
 &                (at(atmres(iaa(rs)%impt(i,4)))%sc(3-shf).eq.iaa(rs)%impt(i,4))) then
                write(ilog,*) 'Fatal. There should be no impropers on the -S-S- atoms in a disulfide. This is a bug.' 
                call fexit()
              end if
            end if
         end if
        end do
        if ((maxval(atmres(iaa(rs)%impt(i,1:4))).ne.minval(atmres(iaa(rs)%impt(i,1:4)))).AND.&
 &  ((seqtyp(atmres(iaa(rs)%impt(i,1))).eq.26).OR.(seqtyp(atmres(iaa(rs)%impt(i,2))).eq.26).OR.&
 &   (seqtyp(atmres(iaa(rs)%impt(i,3))).eq.26).OR.(seqtyp(atmres(iaa(rs)%impt(i,4))).eq.26))) then
          iaa(rs)%par_impt(i,2) = getztor_ref(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
        else
          iaa(rs)%par_impt(i,2) = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
        end if
        iaa(rs)%typ_impt(i) = 2
        iaa(rs)%par_impt(i,1) = 40.0/(RADIAN*RADIAN)
      end do
    end if
!
    allocate(iaa(rs)%cm(nrscm(rs),5))
!    allocate(iaa(rs)%par_cm(nrscm(rs),MAXCMPAR,MAXCMPAR))
    allocate(iaa(rs)%typ_cm(nrscm(rs)))
    iaa(rs)%typ_cm(:) = 0
!    write(*,*) 'dihedrals'
    do i=1,nrscm(rs)
      do j=1,5
        iaa(rs)%cm(i,j) = tempm(i,j)
      end do
    end do
!
!   clean up temporaries
    deallocate(tempa)
    deallocate(tempb)
    deallocate(tempc)
    deallocate(tempd)
    deallocate(tempm)
!    write(*,*) rs,': ',nrsba(rs)
!    call force_torsions(rs,esterms)
!
  end do
!
! finally, give option to write report (mostly for debugging)
!
  if ((ia_report.EQV..true.).AND.(delay_report.EQV..false.)) then
    call interreport()
  end if
!
end
!
!-----------------------------------------------------------------------------------------------------------
!
subroutine interreport()
!
  use iounit
  use inter
  use polypep
  use atoms
  use sequen
  use aminos
!
  implicit none
!
  integer MAXNBATM
  parameter (MAXNBATM=500)
!
  integer, ALLOCATABLE:: mm(:,:)
  integer rs,rs1,rs2,i,zz,ii,kk,j,k
  character(3) resname
!
  allocate(mm(MAXNBATM,MAXNBATM))
  mm(:,:) = 0
!
 44   format('Residue ',i3,' (',a3,') has ',i3,' total atoms:')
 46   format(i7,' ',i7,' ',f6.3)
 47   format('There are ',i4,' intra-residue interactions:')
 48   format('There are ',i4,' to-next-neighbor interactions:')
 49   format('Atom #1 Atom #2 Fudge-factor')
 50   format(5000(i6,1x))
!
  write(ilog,*)
  write(ilog,*) '--- Summary of NB interactions ---'
  write(ilog,*)
  if (n.le.MAXNBATM) then
    write(ilog,*) 'Atom-atom matrix of unique NB interactions (does not reflect crosslink excludes)'
    write(ilog,*)
    zz = 0
    do rs1=1,nseq
      do rs2=rs1,nseq
        if (rs1.eq.rs2) then
          do i=1,nrsintra(rs1)
            zz = zz + 1
            ii = iaa(rs1)%atin(i,1)
            kk = iaa(rs1)%atin(i,2)
            mm(ii,kk) = zz
            mm(kk,ii) = zz
          end do
        else if (rs2-rs1.eq.1) then
          do i=1,nrsnb(rs1)
            zz = zz + 1
            ii = iaa(rs1)%atnb(i,1)
            kk = iaa(rs1)%atnb(i,2)
            mm(ii,kk) = zz
            mm(kk,ii) = zz
          end do
        else
          do i=1,at(rs1)%nsc+at(rs1)%nbb
            if (i.le.at(rs1)%nbb) then
              ii = at(rs1)%bb(i)
            else
              ii = at(rs1)%sc(i-at(rs1)%nbb)
            end if
            do k=1,at(rs2)%nbb+at(rs2)%nsc
              if (k.le.at(rs2)%nbb) then
                kk = at(rs2)%bb(k)
              else
                kk = at(rs2)%sc(k-at(rs2)%nbb)
              end if
              zz = zz + 1
              mm(ii,kk) = zz
              mm(kk,ii) = zz
            end do
          end do
        end if
      end do
    end do
    do i=1,n
      write(ilog,50) (mm(i,j),j=1,n)
    end do
  end if
  deallocate(mm)
  write(ilog,*)
  do rs=1,nseq
    resname = amino(seqtyp(rs))
    write(ilog,44) rs,resname,at(rs)%nbb+at(rs)%nsc
    write(ilog,*)
    write(ilog,47) nrsintra(rs)
    if (nrsintra(rs).gt.0) then
      write(ilog,49)
      do i=1,nrsintra(rs)
        write(ilog,46) iaa(rs)%atin(i,1),iaa(rs)%atin(i,2),&
 &fudge(rs)%rsin(i)
      end do
    end if
    write(ilog,*)
    if (rs.lt.nseq) then
    write(ilog,48) nrsnb(rs)
    if (nrsnb(rs).gt.0) then
      write(ilog,49)
      do i=1,nrsnb(rs)
        write(ilog,46) iaa(rs)%atnb(i,1),iaa(rs)%atnb(i,2),&
 &fudge(rs)%rsnb(i)
      end do
    end if
    write(ilog,*)
    end if
  end do
!
end
!
!---------------------------------------------------------------------
!
! GROMOS demands NB-excludes similar but not identical to nbsr_model 1, hence a separate
! and non-general fxn
!
  subroutine GROMOS_excludes()
!
  use iounit
  use inter
  use sequen
  use aminos
  use atoms
  use polypep
  use molecule
  use system
!
  implicit none
! 
  integer i,j,rs,ii,kk,hb1,shf2,shf3
  character(3) resname
  logical its14
!
  shf2 = 0
  shf3 = 0
  if (ua_model.gt.0) then 
    shf2 = 2
    shf3 = 3
    if (ua_model.eq.2) then
      write(ilog,*) 'Fatal. The GROMOS exclusion rules are largely specific for hydrogens &
 &attached to aromatic rings and make no sense if those are not explicitly represented (choice &
 &of FMCSC_UAMODEL). Adjust setting for FMCSC_INTERMODEL.'
      call fexit()
    end if
  end if
!
  do rs=1,nseq
    do i=1,nrsintra(rs)
!
      resname = amino(seqtyp(rs))
      ii = iaa(rs)%atin(i,1)
      kk = iaa(rs)%atin(i,2)
      if ((attyp(ii).eq.0).OR.(attyp(kk).eq.0)) cycle
!
      its14 = .false.
!
      do j=1,n14(ii)
        if (i14(j,ii).eq.kk) its14 = .true.
      end do
      if (its14.EQV..true.) then
        if ((resname.eq.'PHE').OR.(resname.eq.'TRP').OR.(resname.eq.'TYO').OR.&
 &          (resname.eq.'HIE').OR.(resname.eq.'HID').OR.(resname.eq.'HIP')) then
          if (ua_model.eq.0) then
            hb1 = 8
            if (resname.eq.'PHE') hb1 = 9
            if (resname.eq.'TYO') hb1 = 10
            if (resname.eq.'TRP') hb1 = 12
!           both are sidechain atoms but not beta-hydrogens
            if (((ii.ge.at(rs)%sc(2)).AND.(ii.ne.at(rs)%sc(hb1)).AND.(ii.ne.at(rs)%sc(hb1+1))).AND.&
 &            ((kk.ge.at(rs)%sc(2)).AND.(kk.ne.at(rs)%sc(hb1)).AND.(kk.ne.at(rs)%sc(hb1+1)))) then
              fudge(rs)%rsin(i) = 0.0
            end if
          else
!           both are sidechain atoms
            if ((ii.ge.at(rs)%sc(1)).AND.(kk.ge.at(rs)%sc(1))) then
              fudge(rs)%rsin(i) = 0.0
            end if
          end if
        else if (resname.eq.'PTR') then
          write(ilog,*) 'Fatal. GROMOS exclusion rules are currently not supported for residue type PTR. Please &
 &check back later.'
          call fexit()
        else if (resname.eq.'TYR') then
          hb1 = 10
!         both are sidechain atoms but not beta-hydrogens and not phenolic hydrogen
          if (ua_model.eq.0) then
            if (((ii.ge.at(rs)%sc(2)).AND.(ii.ne.at(rs)%sc(hb1)).AND.&
 &             (ii.ne.at(rs)%sc(hb1+1)).AND.(ii.ne.at(rs)%sc(16))).AND.&
 &            ((kk.ge.at(rs)%sc(2)).AND.(kk.ne.at(rs)%sc(hb1)).AND.&
 &             (kk.ne.at(rs)%sc(hb1+1)).AND.(kk.ne.at(rs)%sc(16)))) then
              fudge(rs)%rsin(i) = 0.0
            end if
          else
!           both are sidechain atoms and not phenolic hydrogen
            if (((ii.ge.at(rs)%sc(1)).AND.(ii.ne.at(rs)%sc(at(rs)%nsc))).AND.&
 &              ((kk.ge.at(rs)%sc(1)).AND.(kk.ne.at(rs)%sc(at(rs)%nsc)))) then
              fudge(rs)%rsin(i) = 0.0
            end if
          end if
        else if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(19-shf2)).AND.&
 &            (kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(19-shf2))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(19-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(19-shf2)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'RPU').OR.(resname.eq.'RIU')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(kk.ge.at(rs)%sc(8-shf2))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'RPC').OR.(resname.eq.'RIC')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(17-shf2)).AND.&
 &            (kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(17-shf2))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(17-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(17-shf2)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'RPG').OR.(resname.eq.'RIG')) then
          if ((ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(20-shf2)).AND.&
 &            (kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(20-shf2))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(8-shf2)).AND.(kk.le.at(rs)%sc(20-shf2))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(8-shf2)).AND.(ii.le.at(rs)%sc(20-shf2)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(18-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(18-shf3))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(18-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(18-shf3)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'DPT').OR.(resname.eq.'DIT')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(17-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(17-shf3))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(17-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(17-shf3)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'DPC').OR.(resname.eq.'DIC')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(16-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(16-shf3))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(16-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(16-shf3)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        else if ((resname.eq.'DPG').OR.(resname.eq.'DIG')) then
          if ((ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(19-shf3)).AND.&
 &            (kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(19-shf3))) then
            fudge(rs)%rsin(i) = 0.0
          else if (((ii.eq.at(rs)%sc(2)).AND.(kk.ge.at(rs)%sc(7-shf3)).AND.(kk.le.at(rs)%sc(19-shf3))).OR.&
 &                 ((kk.eq.at(rs)%sc(2)).AND.(ii.ge.at(rs)%sc(7-shf3)).AND.(ii.le.at(rs)%sc(19-shf3)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        end if
      end if
      if ((resname.eq.'RPA').OR.(resname.eq.'RIA')) then
        if (((ii.eq.at(rs)%sc(12-shf2)).AND.((kk.eq.at(rs)%sc(20-shf2)).OR.(kk.eq.at(rs)%sc(21-shf2)))).OR.&
            ((kk.eq.at(rs)%sc(12-shf2)).AND.((ii.eq.at(rs)%sc(20-shf2)).OR.(ii.eq.at(rs)%sc(21-shf2))))) then
          fudge(rs)%rsin(i) = 0.0
        end if
      else if ((resname.eq.'DPA').OR.(resname.eq.'DIA')) then
        if (((ii.eq.at(rs)%sc(11-shf3)).AND.((kk.eq.at(rs)%sc(19-shf3)).OR.(kk.eq.at(rs)%sc(20-shf3)))).OR.&
            ((kk.eq.at(rs)%sc(11-shf3)).AND.((ii.eq.at(rs)%sc(19-shf3)).OR.(ii.eq.at(rs)%sc(20-shf3))))) then
          fudge(rs)%rsin(i) = 0.0
        end if
      else if (resname.eq.'RPU') then
        if (((ii.eq.at(rs)%bb(1)).AND.(kk.eq.at(rs)%sc(18-shf2))).OR.&
            ((kk.eq.at(rs)%bb(1)).AND.(ii.eq.at(rs)%sc(18-shf2)))) then
          fudge(rs)%rsin(i) = 0.0
        end if
        if (rs.eq.rsmol(molofrs(rs),2)) then ! O3* vs. O2* exclusion
          if (((ii.eq.at(rs)%sc(4)).AND.(kk.eq.at(rs)%bb(9))).OR.&
            ((kk.eq.at(rs)%sc(4)).AND.(ii.eq.at(rs)%bb(9)))) then
            fudge(rs)%rsin(i) = 0.0
          end if
        end if
      end if
!      if (resname.eq.'TRP') then
!!       15-interactions H(N) to CE3 and CH2 should be gone but aren't --> ):
!        if (((ii.eq.at(rs)%sc(15)).AND.((kk.eq.at(rs)%sc(8)).OR.(kk.eq.at(rs)%sc(11)))).OR.&
! &              ((kk.eq.at(rs)%sc(15)).AND.((ii.eq.at(rs)%sc(8)).OR.(ii.eq.at(rs)%sc(11))))) then
!              fudge(rs)%rsin(i) = 0.0
!        end if
!      end if
    end do
    do i=1,nrsnb(rs)
!
      resname = amino(seqtyp(rs))
      ii = iaa(rs)%atnb(i,1)
      kk = iaa(rs)%atnb(i,2)
      if ((attyp(ii).eq.0).OR.(attyp(kk).eq.0)) cycle
!
      if (resname.eq.'RPU') then
        if (((ii.eq.at(rs)%sc(4)).AND.(kk.eq.at(rs+1)%bb(1))).OR.&
            ((kk.eq.at(rs)%sc(4)).AND.(ii.eq.at(rs+1)%bb(1)))) then
          fudge(rs)%rsnb(i) = 0.0
        end if
      end if
    end do
  end do
!
end
!
!--------------------------------------------------------------------------------------------
!
! get exclusion lists for cross-links
!
subroutine crosslink_excludes()
!
  use iounit
  use sequen
  use polypep
  use atoms
  use inter
  use params
  use energies
!
  implicit none
!
  integer i,j,k,l,rs1,rs2,ii,kk
  logical proceed,its14
!
  do j=1,n_crosslinks
    crosslink(j)%nrsin = 0
    if ((crosslink(j)%itstype.eq.1).OR.(crosslink(j)%itstype.eq.2)) then
      rs1 = crosslink(j)%rsnrs(1)
      rs2 = crosslink(j)%rsnrs(2)
      allocate(crosslink(j)%exclin(at(rs1)%na*at(rs2)%na,2))
      allocate(crosslink(j)%ljpars(at(rs1)%na*at(rs2)%na,2))
      allocate(crosslink(j)%is14in(at(rs1)%na*at(rs2)%na))
!
      do i=1,at(rs1)%nbb+at(rs1)%nsc
        if (i.le.at(rs1)%nbb) then
          ii = at(rs1)%bb(i)
        else
          ii = at(rs1)%sc(i-at(rs1)%nbb)
        end if
        if (attyp(ii) .le. 0)  cycle
        do k=1,at(rs2)%nbb+at(rs2)%nsc
          if (k.le.at(rs2)%nbb) then
            kk = at(rs2)%bb(k)
          else
            kk = at(rs2)%sc(k-at(rs2)%nbb)
          end if
          if (attyp(kk) .le. 0)  cycle
          proceed = .false.
          its14 = .false.
!
!         this assumes there cannot be any non-canonical excludes between rs1,rs2
          do l=1,n12(ii)
            if (i12(l,ii).eq.kk) proceed = .true.
          end do
          if (proceed.EQV..false.) then
            do l=1,n13(ii)
              if (i13(l,ii).eq.kk) proceed = .true.
            end do
          end if
          if (proceed.EQV..false.) then
            do l=1,n14(ii)
              if (use_14.EQV..false.) then
                if (i14(l,ii).eq.kk) proceed = .true.
              else
                if (i14(l,ii).eq.kk) then
                  its14 = .true.
                  proceed = .true.
                end if
              end if
            end do
          end if
          if (proceed.EQV..true.) then
            crosslink(j)%nrsin = crosslink(j)%nrsin + 1
            crosslink(j)%exclin(crosslink(j)%nrsin,1) = ii
            crosslink(j)%exclin(crosslink(j)%nrsin,2) = kk
            if (its14.EQV..true.) then 
              crosslink(j)%ljpars(crosslink(j)%nrsin,1) = lj_sig_14(attyp(ii),attyp(kk))
              crosslink(j)%ljpars(crosslink(j)%nrsin,2) = (1.0-fudge_st_14)*lj_eps_14(attyp(ii),attyp(kk))
              crosslink(j)%is14in(crosslink(j)%nrsin) = .true.
            else
              crosslink(j)%ljpars(crosslink(j)%nrsin,1) = lj_sig(attyp(ii),attyp(kk))
              crosslink(j)%ljpars(crosslink(j)%nrsin,2) = lj_eps(attyp(ii),attyp(kk))
              crosslink(j)%is14in(crosslink(j)%nrsin) = .false.
            end if
          end if
        end do
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported crosslink type in crosslink_excludes(...).&
 & This is most likely an omission bug.'
      call fexit()
    end if
  end do
!
end
!
!------------------------------------------------------------------------------------
!
! a subroutine specifically to correct the missing exclusion setup for unsupported residues 
! (requires rotation lists to be done)
!
subroutine correct_srinter()
!
  use atoms
  use inter
  use zmatrix
  use sequen
  use molecule
  use polypep
  use params
!
  implicit none
!
  integer rs,ii,kk,i,k
  logical proceed,its14,write_report
!
  write_report = .false.
  if (nbsr_model.ne.1) return
!
  do rs=1,nseq
!
    if (seqtyp(rs).ne.26) cycle
!
    write_report = .true.
    nrsintra(rs) = 0
    do i=1,at(rs)%nbb+at(rs)%nsc
      if (i.le.at(rs)%nbb) then
        ii = at(rs)%bb(i)
      else
        ii = at(rs)%sc(i-at(rs)%nbb)
      end if
      if (attyp(ii).le.0)  cycle
      do k=i+1,at(rs)%nbb+at(rs)%nsc
        if (k.le.at(rs)%nbb) then
          kk = at(rs)%bb(k)
        else
          kk = at(rs)%sc(k-at(rs)%nbb)
        end if
        if (attyp(kk).le.0)  cycle
        proceed = .true.
        its14 = .false.
!
!       standard topology checks and rotation list-based checks
        call ia_rotlsts(ii,kk,proceed,its14)
!
        if ((its14.EQV..true.).AND.(use_14.EQV..false.)) proceed = .false.
!
        if (proceed.EQV..true.) then
          nrsintra(rs) = nrsintra(rs) + 1 
          iaa(rs)%atin(nrsintra(rs),1) = ii
          iaa(rs)%atin(nrsintra(rs),2) = kk
          fudge(rs)%rsin_ljs(nrsintra(rs)) = lj_sig(attyp(ii),attyp(kk))
          fudge(rs)%rsin_lje(nrsintra(rs)) = lj_eps(attyp(ii),attyp(kk))
          if (its14.EQV..true.) then
            fudge(rs)%rsin(nrsintra(rs)) = fudge_st_14
            fudge(rs)%rsin_ljs(nrsintra(rs)) = lj_sig_14(attyp(ii),attyp(kk))
            fudge(rs)%rsin_lje(nrsintra(rs)) = lj_eps_14(attyp(ii),attyp(kk))
          end if
        end if
!
      end do
    end do
  end do
!
  do rs=1,nseq-1
!
    if ((seqtyp(rs).ne.26).AND.(seqtyp(rs+1).ne.26)) cycle
!
    if (molofrs(rs).ne.molofrs(rs+1)) cycle
!
    write_report = .true.
    nrsnb(rs) = 0
!
!   get relevant next-residue in chain contacts for fixed distances/angles
!
    do i=1,at(rs)%nbb+at(rs)%nsc
      if (i.le.at(rs)%nbb) then
        ii = at(rs)%bb(i)
      else
        ii = at(rs)%sc(i-at(rs)%nbb)
      end if
      if (attyp(ii).le.0)  cycle
      do k=1,at(rs+1)%nbb+at(rs+1)%nsc
        if (k.le.at(rs+1)%nbb) then
          kk = at(rs+1)%bb(k)
        else
          kk = at(rs+1)%sc(k-at(rs+1)%nbb)
        end if
        if (attyp(kk).le.0)  cycle
        proceed = .true.
        its14 = .false.
!
!       standard topology checks and rotation list-based checks
        call ia_rotlsts(ii,kk,proceed,its14)
!
        if ((its14.EQV..true.).AND.(use_14.EQV..false.)) proceed = .false.
!
        if (proceed.EQV..true.) then
          nrsnb(rs) = nrsnb(rs) + 1 
          iaa(rs)%atnb(nrsnb(rs),1) = ii
          iaa(rs)%atnb(nrsnb(rs),2) = kk
          fudge(rs)%rsnb_ljs(nrsnb(rs)) = lj_sig(attyp(ii),attyp(kk))
          fudge(rs)%rsnb_lje(nrsnb(rs)) = lj_eps(attyp(ii),attyp(kk))
          if (its14.EQV..true.) then
            fudge(rs)%rsnb(nrsnb(rs)) = fudge_st_14
            fudge(rs)%rsnb_ljs(nrsnb(rs)) = lj_sig_14(attyp(ii),attyp(kk))
            fudge(rs)%rsnb_lje(nrsnb(rs)) = lj_eps_14(attyp(ii),attyp(kk))
          end if
        end if
!
      end do
    end do
!
  end do !loop over residues
!
!
  if ((write_report.EQV..true.).AND.(ia_report.EQV..true.)) then
    call interreport()
  end if
!
end
!
!--------------------------------------------------------------------------------------------------------------
!
! revert residue order in crosslinks that was introduced for hierarchy reasons in randomization
!
subroutine correct_crosslinks()
!
  use sequen
  use polypep
!
  implicit none
!
  integer i,j,k,l,rs1,rs2,kk,ranky(maxval(natres)**2),bus(maxval(natres)**2,3)
  RTYPE bus2(maxval(natres)**2,2)
!
  do i=1,n_crosslinks
!   itstype, olks/nolks, nrsin, nrspol are not affected (neither is crlk_idx)
    if (crosslink(i)%rsnrs(1).lt.crosslink(i)%rsnrs(2)) cycle
    if (crosslink(i)%nrsin.gt.0) then
      do j=1,crosslink(i)%nrsin ! is14in and ljpars are not affected
        kk = crosslink(i)%exclin(j,1)
        crosslink(i)%exclin(j,1) = crosslink(i)%exclin(j,2)
        crosslink(i)%exclin(j,2) = kk
      end do
    end if
    if (crosslink(i)%nrsin.gt.0) then
      do j=1,crosslink(i)%nrspol ! is14pol and cbpars are not affected
        kk = crosslink(i)%exclpol(j,1)
        crosslink(i)%exclpol(j,1) = crosslink(i)%exclpol(j,2)
        crosslink(i)%exclpol(j,2) = kk
      end do
    end if
    kk = crosslink(i)%rsnrs(1)
    crosslink(i)%rsnrs(1) =  crosslink(i)%rsnrs(2)
    crosslink(i)%rsnrs(2) = kk
  end do
!
  do i=1,n_crosslinks
    rs1 = crosslink(i)%rsnrs(1)
    rs2 = crosslink(i)%rsnrs(2)
    if (crosslink(i)%nrsin.gt.0) then
      bus(:,:) = 0
      ranky(1:(at(rs1)%na*at(rs2)%na)) = 0
      kk = 0
      do j=1,at(rs1)%na
        do k=1,at(rs2)%na
          do l=1,crosslink(i)%nrsin
            if ((crosslink(i)%exclin(l,1).eq.(at(rs1)%bb(1)+j-1)).AND.&
 &              (crosslink(i)%exclin(l,2).eq.(at(rs2)%bb(1)+k-1))) then
              kk = kk + 1
              ranky(l) = kk
              bus2(kk,1:2) =  crosslink(i)%ljpars(l,1:2)
              bus(kk,1:2) = crosslink(i)%exclin(l,1:2)
              if (crosslink(i)%is14in(l).EQV..true.) then
                bus(kk,3) = 1
              end if
            end if
          end do
        end do
      end do
      do l=1,crosslink(i)%nrsin
        crosslink(i)%exclin(l,1:2) = bus(l,1:2)
        if (bus(l,3).gt.0) then
          crosslink(i)%is14in(l) = .true.
        else
          crosslink(i)%is14in(l) = .false.
        end if
        crosslink(i)%ljpars(l,1:2) = bus2(l,1:2)
      end do
    end if
    if (crosslink(i)%nrspol.gt.0) then
      bus(:,:) = 0
      ranky(1:(at(rs1)%npol*at(rs2)%npol)) = 0
      kk = 0
      do j=1,at(rs1)%npol
        do k=1,at(rs2)%npol
          do l=1,crosslink(i)%nrspol
            if ((crosslink(i)%exclpol(l,1).eq.at(rs1)%pol(j)).AND.&
 &              (crosslink(i)%exclpol(l,2).eq.at(rs2)%pol(k))) then
              kk = kk + 1
              ranky(l) = kk
              bus2(kk,1) = crosslink(i)%cbpars(l)
              bus(kk,1:2) = crosslink(i)%exclpol(l,1:2)
              if (crosslink(i)%is14pol(l).EQV..true.) then
                bus(kk,3) = 1
              end if
            end if
          end do
        end do
      end do
      do l=1,crosslink(i)%nrspol
        crosslink(i)%exclpol(l,1:2) = bus(l,1:2)
        if (bus(l,3).gt.0) then
          crosslink(i)%is14pol(l) = .true.
        else
          crosslink(i)%is14pol(l) = .false.
        end if
        crosslink(i)%cbpars(l) = bus2(l,1)
      end do
    end if
  end do


end
!
!----------------------------------------------------------------------------------------------------------------
!
