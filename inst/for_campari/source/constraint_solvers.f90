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
!------------------------------------------------------------------------------------------------------
!
subroutine shake_setup()
!
  use shakeetal
  use molecule
  use atoms
  use iounit
  use sequen
  use forces
  use math
  use interfaces
  use zmatrix
  use system
  use mini
  use params, ONLY: MAXBAPAR,MAXBOPAR,be_unsafe
  use energies, ONLY: use_BOND
  use inter
!
  implicit none
!
  integer i,j,k,kk,ii,ii3,ii4,jj,mt,mt2,imol,hits(MAXVLNC*MAXVLNC),nhits,t1,t2,nidx,nangls,azero
  integer ii1,ii2,jj1,jj2,candii,candi,effats,effcns,nidx2,rs,rs1,rs2,rslst(3),blix,baix
  integer rnglst(10)
  RTYPE getblen,getbang,l1,l2,l3,l4,lhh,loh,tmpr
  integer, ALLOCATABLE:: idxlst(:,:),taglst(:),fats(:),rngnrs(:),idxlst2(:,:)
  RTYPE, ALLOCATABLE:: eqdlst(:)
  logical foundit,atoms_areinset
  character(MAXSTRLEN) stringy
!
  azero = 0
  cart_cons_grps = 0
  if ((cart_cons_mode.eq.1).AND.(add_settle.EQV..false.)) then
    rs1 = 0
    do i=1,n
      if (mass(i).le.0.0) rs1 = rs1 + 1
    end do
    if (rs1.gt.0) then
      if ((dyn_mode.ne.6).OR.(mini_mode.eq.4)) then
        write(ilog,*) 'Fatal. Massless (virtual) particles (sites) are fatal if in Cartesian dynamics they &
 &are not rigorously connected to a constraint group.'
        call fexit()
      end if
    end if
    no_shake = .true.
    return
  end if
!
! WARNING: using hard-coded parameters for TipNp and SPC geometry
  if (add_settle.EQV..true.) then
    l1 = 0.9572d0
    l2 = 1.0d0
    l3 = l1*sin(104.52/RADIAN)/sin(0.5*(180.0-104.52)/RADIAN) ! 1.513901
    l4 = l2*sin(109.47/RADIAN)/sin(0.5*(180.0-109.47)/RADIAN) ! 1.632981
!    l5 = 0.15d0 ! O-MW in TIP4P
!    l6 = sqrt(0.25*l3*l3 + (cos(0.5*104.52/RADIAN)*l1-l5)**2) ! 0.873480
!    l7 = 0.15d0 ! O-MW in TIP4P-Ew
!    l8 = sqrt(0.25*l3*l3 + (cos(0.5*104.52/RADIAN)*l1-l6)**2)
  end if
!
  allocate(idxlst(max(6*n-6,1),4))
  allocate(eqdlst(max(6*n-6,1)))
  allocate(taglst(max(6*n-6,1)))
  allocate(rngnrs(n))
!
  do i=1,n
    call determine_upto6ring(i,rngnrs(i),rnglst)
  end do
!
  nidx = 0
  eqdlst(:) = 0.0
  idxlst(:,:) = 0
  do imol=1,nmol
    if (cart_cons_mode.eq.6) exit
    do i=atmol(imol,1),atmol(imol,2)
      if (cart_cons_mode.eq.1) exit
      if (n12(i).eq.0) cycle
      do j=1,n12(i)
        candii = -1
        if (cart_cons_mode.eq.2) then ! all terminal bonds with "hydrogen"
          if ((atnam(i12(j,i))(1:1).eq.'H').AND.(mass(i12(j,i)).le.3.5).AND.&
 &            (n12(i12(j,i)).eq.1).AND.(i12(j,i).gt.i)) candii = i12(j,i)
        else if ((cart_cons_mode.ge.3).AND.(cart_cons_mode.le.5)) then ! at least all bonds
          if ((i12(j,i).gt.i)) candii = i12(j,i)
        else
          write(ilog,*) 'Fatal. Encountered unsupported constraint set selection in shake_setup(). This&
 & is an omission bug.'
          call fexit()
        end if
        if (candii.gt.0) then
          nidx = nidx + 1
          idxlst(nidx,1) = i
          idxlst(nidx,2) = candii
          idxlst(nidx,3:4) = 0
          eqdlst(nidx) = getblen(i,candii)
          cart_cons_grps = cart_cons_grps + 1
          taglst(nidx) = cart_cons_grps
          nhits = 0
          do k=1,nidx-1
            if ((i.eq.idxlst(k,1)).OR.(candii.eq.idxlst(k,1)).OR.(i.eq.idxlst(k,2)).OR.(candii.eq.idxlst(k,2))) then
              if (nhits.eq.2) then
                if ((taglst(k).ne.taglst(hits(1))).AND.(taglst(k).ne.taglst(hits(2)))) then
                  write(ilog,*) 'Fatal. Inconsistent setup of SHAKE units. This is a bug.'
                  call fexit()
                end if
              else if (nhits.eq.1) then
                if (taglst(k).ne.taglst(hits(1))) then
                  nhits = nhits + 1
                  hits(nhits) = k
                end if
              else
                nhits = nhits + 1
                hits(nhits) = k
              end if
            end if
          end do
          if (nhits.eq.0) then ! do nothing
          else if (nhits.eq.1) then ! merge into existing group
            taglst(nidx) = taglst(hits(1))
            cart_cons_grps = cart_cons_grps - 1
          else if (nhits.eq.2) then ! merge two newly joined groups
            t1 = min(taglst(hits(1)),taglst(hits(2)))
            t2 = max(taglst(hits(1)),taglst(hits(2)))
            do k=1,nidx-1
              if (taglst(k).eq.t2) taglst(k) = t1
            end do
            taglst(nidx) = t1
            cart_cons_grps = cart_cons_grps - 2
            do k=1,nidx-1
              if (taglst(k).gt.t2) taglst(k) = taglst(k) - 1
            end do
          else 
            write(ilog,*) 'Fatal. Inconsistent setup of SHAKE units. This is a bug.'
            call fexit()
          end if
        end if
      end do
      if ((cart_cons_mode.eq.4).OR.(cart_cons_mode.eq.5)) then
        nangls = 0
        do j=1,n12(i)
          do jj=j+1,n12(i)
            candi = i12(j,i)
            candii = i12(jj,i)
!           we cycle out intra-ring angles since otherwise system accumulates
!           even more redundant constraints
            if ((n12(i).gt.2).AND.(rngnrs(candi).gt.0).AND.(rngnrs(candii).gt.0).AND.(rngnrs(i).gt.0)) cycle
            nangls = nangls + 1
            nidx = nidx + 1
            if (cart_cons_mode.eq.4) then
              idxlst(nidx,1) = candi
              idxlst(nidx,2) = i
              idxlst(nidx,3) = candii
              idxlst(nidx,4) = 0
              eqdlst(nidx) = cos(getbang(candi,i,candii)/RADIAN)
            else
              idxlst(nidx,1) = candi
              idxlst(nidx,2) = candii
              idxlst(nidx,3:4) = 0
              eqdlst(nidx) = getblen(candi,candii)
            end if
            cart_cons_grps = cart_cons_grps + 1
            taglst(nidx) = cart_cons_grps
            nhits = 0
            do k=1,nidx-1
              if ((candi.eq.idxlst(k,1)).OR.(candi.eq.idxlst(k,2)).OR.&
 &                (candii.eq.idxlst(k,1)).OR.(candii.eq.idxlst(k,2)).OR.&
 &       (((i.eq.idxlst(k,1)).OR.(i.eq.idxlst(k,2))).AND.(cart_cons_mode.eq.4))) then
                if (nhits.eq.2) then
                  if ((taglst(k).ne.taglst(hits(1))).AND.(taglst(k).ne.taglst(hits(2)))) then
                    write(ilog,*) 'Fatal. Inconsistent setup of SHAKE units. This is a bug.'
                    call fexit()
                  end if
                else if (nhits.eq.1) then
                  if (taglst(k).ne.taglst(hits(1))) then
                    nhits = nhits + 1
                    hits(nhits) = k
                  end if
                else
                  nhits = nhits + 1
                  hits(nhits) = k
                end if
              end if
            end do
            if (nhits.eq.0) then ! do nothing
            else if (nhits.eq.1) then ! merge into existing group
              taglst(nidx) = taglst(hits(1))
              cart_cons_grps = cart_cons_grps - 1
            else if (nhits.eq.2) then ! merge two newly joined groups
              t1 = min(taglst(hits(1)),taglst(hits(2)))
              t2 = max(taglst(hits(1)),taglst(hits(2)))
              do k=1,nidx-1
                if (taglst(k).eq.t2) taglst(k) = t1
              end do
              taglst(nidx) = t1
              cart_cons_grps = cart_cons_grps - 2
              do k=1,nidx-1
                if (taglst(k).gt.t2) taglst(k) = taglst(k) - 1
              end do
            else 
              write(ilog,*) 'Fatal. Inconsistent setup of SHAKE units. This is a bug.'
              call fexit()
            end if
            if (nangls.eq.(3*(n12(i)+1)-6-n12(i))) exit
!            if (nangls.eq.(n12(i)-1)) exit
          end do
          if (nangls.eq.(3*(n12(i)+1)-6-n12(i))) exit
!          if (nangls.eq.(n12(i)-1)) exit
        end do
      end if
    end do   
    if (add_settle.EQV..true.) then
      if ((rsmol(imol,2).eq.rsmol(imol,1)).AND.((seqtyp(rsmol(imol,1)).eq.40).OR.(seqtyp(rsmol(imol,1)).eq.102).OR.&
 &         (seqtyp(rsmol(imol,1)).eq.45).OR.(seqtyp(rsmol(imol,1)).eq.39).OR.(seqtyp(rsmol(imol,1)).eq.103))) then
        if (seqtyp(rsmol(imol,1)).eq.39) then ! SPC
          lhh = l4
          loh = l2
        else
          lhh = l3
          loh = l1
        end if
        k = nidx
        do while (nidx.gt.0)
          if     ((idxlst(nidx,1).ge.atmol(imol,1)).AND.(idxlst(nidx,2).ge.atmol(imol,1)).AND.&
 &                (idxlst(nidx,1).le.atmol(imol,2)).AND.(idxlst(nidx,2).le.atmol(imol,2)).AND.&
 &        (((idxlst(nidx,3).le.atmol(imol,2)).AND.(idxlst(nidx,3).ge.atmol(imol,1))).OR.(idxlst(nidx,3).eq.0))) then
!           do nothing
          else 
            exit
          end if
          nidx = nidx - 1
        end do
        if (k.eq.nidx) cart_cons_grps = cart_cons_grps + 1
        nidx = nidx + 1
        idxlst(nidx,1) = atmol(imol,1)
        idxlst(nidx,2) = atmol(imol,1)+1
        eqdlst(nidx) = loh
        nidx = nidx + 1
        idxlst(nidx,1) = atmol(imol,1)
        idxlst(nidx,2) = atmol(imol,1)+2
        eqdlst(nidx) = loh
        nidx = nidx + 1
        idxlst(nidx,1) = atmol(imol,1)+1
        idxlst(nidx,2) = atmol(imol,1)+2
        eqdlst(nidx) = lhh
        if ((seqtyp(rsmol(imol,1)).eq.45).OR.(seqtyp(rsmol(imol,1)).eq.102).OR.(seqtyp(rsmol(imol,1)).eq.103)) then
!         note that the constraint lengths involving virtual sites are never actually used
          nidx = nidx + 1
          idxlst(nidx,1) = atmol(imol,1)
          idxlst(nidx,2) = atmol(imol,1)+3
          eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
          nidx = nidx + 1
          idxlst(nidx,1) = atmol(imol,1)+1
          idxlst(nidx,2) = atmol(imol,1)+3
          eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
          nidx = nidx + 1
          idxlst(nidx,1) = atmol(imol,1)+2
          idxlst(nidx,2) = atmol(imol,1)+3
          eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
          if (seqtyp(rsmol(imol,1)).eq.102) then
            nidx = nidx + 1
            idxlst(nidx,1) = atmol(imol,1)
            idxlst(nidx,2) = atmol(imol,1)+4
            eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
            nidx = nidx + 1
            idxlst(nidx,1) = atmol(imol,1)+1
            idxlst(nidx,2) = atmol(imol,1)+4
            eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
            nidx = nidx + 1
            idxlst(nidx,1) = atmol(imol,1)+2
            idxlst(nidx,2) = atmol(imol,1)+4
            eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
            taglst(nidx-8:nidx) = cart_cons_grps
            idxlst(nidx-8:nidx,3:4) = 0
          else
            taglst(nidx-5:nidx) = cart_cons_grps
            idxlst(nidx-5:nidx,3:4) = 0
          end if
        else
          taglst(nidx-2:nidx) = cart_cons_grps
          idxlst(nidx-2:nidx,3:4) = 0
        end if    
      end if
    end if
  end do
  if (cart_cons_mode.eq.6) then
    call read_cartconsfile(idxlst,nidx)
    if (add_settle.EQV..true.) then
      allocate(idxlst2(nidx,4))
      idxlst2(1:nidx,1:4) = idxlst(1:nidx,1:4)
      nidx2 = nidx
!     manually remove preexisting constraints on water molecules
      foundit = .false.
      nidx = 0
      do i=1,nidx2
        rs1 = atmres(idxlst(i,1))
        rs2 = atmres(idxlst(i,2))
        if ((rs1.eq.rs2).AND.((seqtyp(rs1).eq.40).OR.(seqtyp(rs1).eq.45).OR.(seqtyp(rs1).eq.39).OR.&
 &                            (seqtyp(rs1).eq.102).OR.(seqtyp(rs1).eq.103))) then
          foundit = .true.
        else
          nidx = nidx + 1
          idxlst(nidx,:) = idxlst2(i,:)
        end if
      end do
      if (foundit.EQV..true.) then
        write(ilog,*) 'Warning. Removed one or more constraints from file input (FMCSC_SHAKESET 6) that &
 &were already covered by the choice for FMCSC_SETTLEH2O.'
      end if
      deallocate(idxlst2)
      do imol=1,nmol
        if ((rsmol(imol,2).eq.rsmol(imol,1)).AND.((seqtyp(rsmol(imol,1)).eq.40).OR.(seqtyp(rsmol(imol,1)).eq.102).OR.&
 &         (seqtyp(rsmol(imol,1)).eq.45).OR.(seqtyp(rsmol(imol,1)).eq.39).OR.(seqtyp(rsmol(imol,1)).eq.103))) then
          if (seqtyp(rsmol(imol,1)).eq.39) then
            lhh = l4
            loh = l2
          else
            lhh = l3
            loh = l1
          end if
          nidx = nidx + 1
          idxlst(nidx,1) = atmol(imol,1)
          idxlst(nidx,2) = atmol(imol,1)+1
          eqdlst(nidx) = loh
          nidx = nidx + 1
          idxlst(nidx,1) = atmol(imol,1)
          idxlst(nidx,2) = atmol(imol,1)+2
          eqdlst(nidx) = loh
          nidx = nidx + 1
          idxlst(nidx,1) = atmol(imol,1)+1
          idxlst(nidx,2) = atmol(imol,1)+2
          eqdlst(nidx) = lhh
          if ((seqtyp(rsmol(imol,1)).eq.45).OR.(seqtyp(rsmol(imol,1)).eq.102).OR.(seqtyp(rsmol(imol,1)).eq.103)) then
            nidx = nidx + 1
            idxlst(nidx,1) = atmol(imol,1)
            idxlst(nidx,2) = atmol(imol,1)+3
            eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
            nidx = nidx + 1
            idxlst(nidx,1) = atmol(imol,1)+1
            idxlst(nidx,2) = atmol(imol,1)+3
            eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
            nidx = nidx + 1
            idxlst(nidx,1) = atmol(imol,1)+2
            idxlst(nidx,2) = atmol(imol,1)+3
            eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
            if (seqtyp(rsmol(imol,1)).eq.102) then
              nidx = nidx + 1
              idxlst(nidx,1) = atmol(imol,1)
              idxlst(nidx,2) = atmol(imol,1)+4
              eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
              nidx = nidx + 1
              idxlst(nidx,1) = atmol(imol,1)+1
              idxlst(nidx,2) = atmol(imol,1)+4
              eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
              nidx = nidx + 1
              idxlst(nidx,1) = atmol(imol,1)+2
              idxlst(nidx,2) = atmol(imol,1)+4
              eqdlst(nidx) = getblen(idxlst(nidx,1),idxlst(nidx,2))
              idxlst(nidx-8:nidx,3:4) = 0
            else
              idxlst(nidx-5:nidx,3:4) = 0
            end if
          else
            idxlst(nidx-2:nidx,3:4) = 0
          end if  
!          cart_cons_grps = cart_cons_grps + 1
!          taglst(nidx-2:nidx) = cart_cons_grps
        end if
      end do
    end if
    do i=1,nidx
      if (eqdlst(i).le.0.0) eqdlst(i) = getblen(idxlst(i,1),idxlst(i,2))
      cart_cons_grps = cart_cons_grps + 1
      taglst(i) = cart_cons_grps
      nhits = 0
      do k=1,i-1
        if ((idxlst(i,1).eq.idxlst(k,1)).OR.(idxlst(i,2).eq.idxlst(k,1)).OR.&
 &          (idxlst(i,1).eq.idxlst(k,2)).OR.(idxlst(i,2).eq.idxlst(k,2))) then
          if (nhits.eq.2) then
            if ((taglst(k).ne.taglst(hits(1))).AND.(taglst(k).ne.taglst(hits(2)))) then
              write(ilog,*) 'Fatal. Inconsistent setup of SHAKE units. This is a bug.'
              call fexit()
            end if
          else if (nhits.eq.1) then
            if (taglst(k).ne.taglst(hits(1))) then
              nhits = nhits + 1
              hits(nhits) = k
            end if
          else
            nhits = nhits + 1
            hits(nhits) = k
          end if
        end if
      end do
      if (nhits.eq.0) then ! do nothing
      else if (nhits.eq.1) then ! merge into existing group
        taglst(i) = taglst(hits(1))
        cart_cons_grps = cart_cons_grps - 1
      else if (nhits.eq.2) then ! merge two newly joined groups
        t1 = min(taglst(hits(1)),taglst(hits(2)))
        t2 = max(taglst(hits(1)),taglst(hits(2)))
        do k=1,i-1
          if (taglst(k).eq.t2) taglst(k) = t1
        end do
        taglst(i) = t1
        cart_cons_grps = cart_cons_grps - 2
        do k=1,i-1
          if (taglst(k).gt.t2) taglst(k) = taglst(k) - 1
        end do
      else 
        write(ilog,*) 'Fatal. Inconsistent setup of SHAKE units. This is a bug.'
        call fexit()
      end if
    end do
  end if
!
  deallocate(constraints)
  allocate(constraints(cart_cons_grps))
!
  do i=1,cart_cons_grps
    constraints(i)%nr = 0
    constraints(i)%nr3 = 0
    constraints(i)%nr4 = 0
    constraints(i)%lorder = lincs_order
    constraints(i)%iters = 0
    constraints(i)%itercnt = 0
    constraints(i)%isrigid = .false.
    do k=1,nidx
      if (taglst(k).eq.i) then
        if (idxlst(k,3).eq.0) then
          constraints(i)%nr = constraints(i)%nr + 1
        else if (idxlst(k,4).eq.0) then
          constraints(i)%nr3 = constraints(i)%nr3 + 1
        else if (idxlst(k,4).gt.0) then
          constraints(i)%nr4 = constraints(i)%nr4 + 1
        end if
      end if
    end do
!
    if (((constraints(i)%nr3.gt.0).OR.(constraints(i)%nr4.gt.0)).AND.(cart_cons_method.ne.1)) then
      write(ilog,*) 'Fatal. Explicit bond angle or improper dihedral angle constraints are&
 & currently only supported with the simple SHAKE procedure. This is a bug.'
      call fexit()
    end if
!
    if (constraints(i)%nr.gt.0) then
      allocate(constraints(i)%idx(constraints(i)%nr+min(3,constraints(i)%nr3),2)) ! annoying extra storage (SETTLE transfer)
      allocate(constraints(i)%eqd2(constraints(i)%nr+min(3,constraints(i)%nr3)))
      constraints(i)%eqd2(:) = 1.0
      if ((cart_cons_method.eq.2).OR.(cart_cons_method.eq.3)) then
!       this may be badly memory-inefficient (worst case = (3N-6)*(3N-6) reals extra storage)
        allocate(constraints(i)%amat(constraints(i)%nr+constraints(i)%nr3,constraints(i)%nr+constraints(i)%nr3))
        allocate(constraints(i)%mapidx(constraints(i)%nr,2))
      end if
    end if
    if (constraints(i)%nr3.gt.0) then
      allocate(constraints(i)%idx3(constraints(i)%nr3,3))
      allocate(constraints(i)%eqa3(constraints(i)%nr3))
      constraints(i)%eqa3(:) = 0.5
    end if
    if (constraints(i)%nr4.gt.0) then
      allocate(constraints(i)%idx4(constraints(i)%nr4,4))
      allocate(constraints(i)%eqi(constraints(i)%nr4))
      constraints(i)%eqi(:) = 0.0
    end if
!
    ii = 0
    ii3 = 0
    ii4 = 0
    do k=1,nidx
      if ((taglst(k).eq.i).AND.(idxlst(k,3).eq.0)) then
        ii = ii + 1
        constraints(i)%idx(ii,1:2) = idxlst(k,1:2)
        constraints(i)%eqd2(ii) = eqdlst(k)
      else if ((taglst(k).eq.i).AND.(idxlst(k,4).eq.0)) then
        ii3 = ii3 + 1
        constraints(i)%idx3(ii3,1:3) = idxlst(k,1:3)
        constraints(i)%eqa3(ii3) = eqdlst(k)
      else if ((taglst(k).eq.i).AND.(idxlst(k,4).gt.0)) then
        ii4 = ii4 + 1
        constraints(i)%idx4(ii4,1:4) = idxlst(k,1:4)
        constraints(i)%eqi(ii4) = eqdlst(k)
      end if
    end do
    if (constraints(i)%nr.gt.0) constraints(i)%eqd2(:) = constraints(i)%eqd2(:)*constraints(i)%eqd2(:)
    n_constraints = n_constraints + constraints(i)%nr + constraints(i)%nr3 + constraints(i)%nr4
!
    allocate(fats(2*constraints(i)%nr+3*constraints(i)%nr3 + 4*constraints(i)%nr4))
    mt = 0
    allocate(constraints(i)%massless(constraints(i)%nr + constraints(i)%nr3 + constraints(i)%nr4))
    allocate(constraints(i)%isinter(constraints(i)%nr + constraints(i)%nr3 + constraints(i)%nr4))
    constraints(i)%massless(:) = .false.
    constraints(i)%hasmassless = .false.
    constraints(i)%isinter(:) = .false.
!
    do k=1,constraints(i)%nr
      if ((mass(constraints(i)%idx(k,1)).le.0.0).OR.(mass(constraints(i)%idx(k,2)).le.0.0)) then
        constraints(i)%massless(k) = .true.
        constraints(i)%hasmassless = .true.
      end if
      if (molofrs(atmres(constraints(i)%idx(k,1))).ne.molofrs(atmres(constraints(i)%idx(k,2)))) then
        constraints(i)%isinter(k) = .true.
      end if
      do ii=1,2
        if ((cart_frz(constraints(i)%idx(k,ii),1).EQV..true.).OR.(cart_frz(constraints(i)%idx(k,ii),2).EQV..true.).OR.&
 &          (cart_frz(constraints(i)%idx(k,ii),3).EQV..true.)) then
          write(ilog,*) 'Fatal. Concurrent use of explicit atom (FMCSC_FRZFILE) and holonomic constraints (FMCSC_SHAKESET) &
 &is only allowed if no trivially constrained atom is part of any holonomic constraint group.'
          call fexit()
        end if
        foundit = .false.
        do j=1,mt
          if (fats(j).eq.constraints(i)%idx(k,ii)) foundit = .true.
          if (foundit.EQV..true.) exit
        end do
        if (foundit.EQV..false.) then
          mt = mt + 1
          fats(mt) = constraints(i)%idx(k,ii)
          if ((cart_cons_method.eq.2).OR.(cart_cons_method.eq.3)) constraints(i)%mapidx(k,ii) = mt
        else
          if ((cart_cons_method.eq.2).OR.(cart_cons_method.eq.3)) constraints(i)%mapidx(k,ii) = j
        end if
      end do
    end do
    do k=1,constraints(i)%nr3
      if ((mass(constraints(i)%idx3(k,1)).le.0.0).OR.(mass(constraints(i)%idx3(k,2)).le.0.0).OR.&
 &        (mass(constraints(i)%idx3(k,3)).le.0.0)) then
        constraints(i)%massless(constraints(i)%nr+k) = .true.
        constraints(i)%hasmassless = .true.
      end if
      if ((molofrs(atmres(constraints(i)%idx3(k,1))).ne.molofrs(atmres(constraints(i)%idx3(k,2)))).OR.&
 &        (molofrs(atmres(constraints(i)%idx3(k,1))).ne.molofrs(atmres(constraints(i)%idx3(k,3))))) then
        constraints(i)%isinter(constraints(i)%nr+k) = .true.
      end if
      do ii=1,3
        if ((cart_frz(constraints(i)%idx3(k,ii),1).EQV..true.).OR.(cart_frz(constraints(i)%idx3(k,ii),2).EQV..true.).OR.&
 &          (cart_frz(constraints(i)%idx3(k,ii),3).EQV..true.)) then
          write(ilog,*) 'Fatal. Concurrent use of explicit atom (FMCSC_FRZFILE) and holonomic constraints (FMCSC_SHAKESET) &
 &is only allowed if no trivially constrained atom is part of any holonomic constraint group.'
          call fexit()
        end if
        foundit = .false.
        do j=1,mt
          if (fats(j).eq.constraints(i)%idx3(k,ii)) foundit = .true.
          if (foundit.EQV..true.) exit
        end do
        if (foundit.EQV..false.) then
          mt = mt + 1
          fats(mt) = constraints(i)%idx3(k,ii)
        end if
      end do
    end do
!   these have to be intramolecular (isinter is always false)
    do k=1,constraints(i)%nr4
      if ((mass(constraints(i)%idx4(k,1)).le.0.0).OR.(mass(constraints(i)%idx4(k,2)).le.0.0).OR.&
 &        (mass(constraints(i)%idx4(k,3)).le.0.0).OR.(mass(constraints(i)%idx4(k,4)).le.0.0)) then
        constraints(i)%massless(constraints(i)%nr+constraints(i)%nr3+k) = .true.
        constraints(i)%hasmassless = .true.
      end if
      ii1 = molofrs(atmres(constraints(i)%idx4(k,1)))
      do ii=1,4
        if ((cart_frz(constraints(i)%idx4(k,ii),1).EQV..true.).OR.(cart_frz(constraints(i)%idx4(k,ii),2).EQV..true.).OR.&
 &          (cart_frz(constraints(i)%idx4(k,ii),3).EQV..true.)) then
          write(ilog,*) 'Fatal. Concurrent use of explicit atom (FMCSC_FRZFILE) and holonomic constraints (FMCSC_SHAKESET) &
 &is only allowed if no trivially constrained atom is part of any holonomic constraint group.'
          call fexit()
        end if
        if (molofrs(atmres(constraints(i)%idx4(k,ii))).ne.ii1) then
          write(ilog,*) 'Fatal. Improper dihedral angle constraints have to consist of atoms that are all in the &
 &same molecule. This is an omission bug.'
          call fexit()
        end if
        foundit = .false.
        do j=1,mt
          if (fats(j).eq.constraints(i)%idx4(k,ii)) foundit = .true.
          if (foundit.EQV..true.) exit
        end do
        if (foundit.EQV..false.) then
          mt = mt + 1
          fats(mt) = constraints(i)%idx4(k,ii)
        end if
      end do
    end do
!
    if (mt.eq.2) then
      if (constraints(i)%nr+constraints(i)%nr3+constraints(i)%nr4+5.gt.3*mt) then
        write(ilog,*) 'Fatal. Constraints group has too few remaining degrees of freedom. This is &
 &either a bug or the result of bad user-input when using mode 6 for FMCSC_SHAKESET.'
        call fexit()
      else if (constraints(i)%nr+constraints(i)%nr3+constraints(i)%nr4+5.eq.3*mt) then
        constraints(i)%isrigid = .true.
      end if
    else
      if (constraints(i)%nr+constraints(i)%nr3+constraints(i)%nr4+6.gt.3*mt) then
        write(ilog,*) 'Fatal. Constraints group has too few remaining degrees of freedom. This is &
 &either a bug or the result of bad user-input when using mode 6 for FMCSC_SHAKESET.'
        call fexit()
      else if (constraints(i)%nr+constraints(i)%nr3+constraints(i)%nr4+6.eq.3*mt) then
        constraints(i)%isrigid = .true.
      end if
    end if
    constraints(i)%nats = mt
    allocate(constraints(i)%uidx(mt))
    constraints(i)%uidx(:) = fats(1:mt)
    deallocate(fats)
  end do 
!
  deallocate(idxlst)
  deallocate(eqdlst)
  deallocate(taglst)
  deallocate(rngnrs)
!
! identify those constraint groups suitable for settle, assemble into 3 lists 
  settle_tip3ps = 0
  settle_tip4ps = 0
  settle_tip4pes = 0
  settle_tip5ps = 0
  settle_spcs = 0
  settle_rest = 0
  shake_cnt = 0
  do mt=1,cart_cons_grps
    shake_cnt = shake_cnt + 1
    effcns = 0
    effats = 0
    do i=1,constraints(mt)%nr+constraints(mt)%nr3+constraints(mt)%nr4
      if (constraints(mt)%massless(i).EQV..false.) effcns = effcns + 1
    end do
    do i=1,constraints(mt)%nats
      if (mass(constraints(mt)%uidx(i)).gt.0.0) effats = effats + 1
    end do
    if ((effcns.eq.3).AND.(effats.eq.3)) then
      if (b_type(constraints(mt)%uidx(2)).ne.b_type(constraints(mt)%uidx(3))) cycle
      if (b_type(constraints(mt)%idx(1,1)).ne.b_type(constraints(mt)%uidx(1))) cycle
      if (b_type(constraints(mt)%idx(2,1)).ne.b_type(constraints(mt)%uidx(1))) cycle
      k = 0
      do i=1,n12(constraints(mt)%uidx(1))
        if (i12(i,constraints(mt)%uidx(1)).eq.constraints(mt)%uidx(2)) k = k + 1
        if (i12(i,constraints(mt)%uidx(1)).eq.constraints(mt)%uidx(3)) k = k + 1
      end do
      if (k.ne.2) cycle
      if ((atmres(constraints(mt)%uidx(1)).eq.atmres(constraints(mt)%uidx(2))).AND.&
 &        (atmres(constraints(mt)%uidx(1)).eq.atmres(constraints(mt)%uidx(3)))) then
!       note that the TIPNP models have the same HOH geometry
        if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.40) then
          settle_tip3ps = settle_tip3ps + 1
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.45) then
          settle_tip4ps = settle_tip4ps + 1
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.102) then
          settle_tip5ps = settle_tip5ps + 1
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.103) then
          settle_tip4pes = settle_tip4pes + 1
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.39) then ! SPC
          settle_spcs = settle_spcs + 1
        else
          settle_rest = settle_rest + 1
        end if
      else
        if (molofrs(atmres(constraints(mt)%uidx(1))).ne.molofrs(atmres(constraints(mt)%uidx(3)))) then
!         boundary wraparound not supported, so leave at SHAKE
          cycle
        end if
        settle_rest = settle_rest + 1
      end if
    else
      cycle
    end if
    shake_cnt = shake_cnt - 1
  end do
  if (settle_tip3ps.gt.0) allocate(settlelst1(settle_tip3ps))
  if (settle_spcs.gt.0) allocate(settlelst2(settle_spcs))
  if (settle_rest.gt.0) allocate(settlelst3(settle_rest))
  if (settle_tip4ps.gt.0) allocate(settlelst4(settle_tip4ps))
  if (settle_tip4pes.gt.0) allocate(settlelst6(settle_tip4pes))
  if (settle_tip5ps.gt.0) allocate(settlelst5(settle_tip5ps))

  allocate(nonsettlelst(shake_cnt+1))
  settle_tip3ps = 0
  settle_tip4ps = 0
  settle_tip4pes = 0
  settle_tip5ps = 0
  settle_spcs = 0
  settle_rest = 0
  shake_cnt = 0
  do mt=1,cart_cons_grps
    shake_cnt = shake_cnt + 1
    effcns = 0
    effats = 0
    do i=1,constraints(mt)%nr+constraints(mt)%nr3+constraints(mt)%nr4
      if (constraints(mt)%massless(i).EQV..false.) effcns = effcns + 1
    end do
    do i=1,constraints(mt)%nats
      if (mass(constraints(mt)%uidx(i)).gt.0.0) effats = effats + 1
    end do
    nonsettlelst(shake_cnt) = mt
    if ((effcns.eq.3).AND.(effats.eq.3)) then
      if (b_type(constraints(mt)%uidx(2)).ne.b_type(constraints(mt)%uidx(3))) cycle
      if (b_type(constraints(mt)%idx(1,1)).ne.b_type(constraints(mt)%uidx(1))) cycle
      if (b_type(constraints(mt)%idx(2,1)).ne.b_type(constraints(mt)%uidx(1))) cycle
      k = 0
      do i=1,n12(constraints(mt)%uidx(1))
        if (i12(i,constraints(mt)%uidx(1)).eq.constraints(mt)%uidx(2)) k = k + 1
        if (i12(i,constraints(mt)%uidx(1)).eq.constraints(mt)%uidx(3)) k = k + 1
      end do
      if (k.ne.2) cycle
      if (constraints(mt)%eqd2(1).ne.constraints(mt)%eqd2(2)) then
        tmpr = 0.5*(sqrt(constraints(mt)%eqd2(1)) + sqrt(constraints(mt)%eqd2(2)))
        constraints(mt)%eqd2(1:2) = tmpr*tmpr
      end if
      if ((atmres(constraints(mt)%uidx(1)).eq.atmres(constraints(mt)%uidx(2))).AND.&
 &        (atmres(constraints(mt)%uidx(1)).eq.atmres(constraints(mt)%uidx(3)))) then
        if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.40) then
          settle_tip3ps = settle_tip3ps + 1
          settlelst1(settle_tip3ps) = mt
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.45) then
          settle_tip4ps = settle_tip4ps + 1
          settlelst4(settle_tip4ps) = mt
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.102) then
          settle_tip5ps = settle_tip5ps + 1
          settlelst5(settle_tip5ps) = mt
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.103) then
          settle_tip4pes = settle_tip4pes + 1
          settlelst6(settle_tip4pes) = mt
        else if (seqtyp(atmres(constraints(mt)%uidx(1))).eq.39) then ! SPC
          settle_spcs = settle_spcs + 1
          settlelst2(settle_spcs) = mt
        else if (cart_cons_source.eq.1) then
          settle_rest = settle_rest + 1
          settlelst3(settle_rest) = mt
        else
          cycle
        end if
      else if (cart_cons_source.eq.1) then
        if (molofrs(atmres(constraints(mt)%uidx(1))).ne.molofrs(atmres(constraints(mt)%uidx(3)))) then
!         boundary wraparound not supported, so leave at SHAKE
          cycle
        end if
        settle_rest = settle_rest + 1
        settlelst3(settle_rest) = mt
      else
        cycle
      end if
    else
      cycle
    end if
    shake_cnt = shake_cnt - 1
  end do
!
! disable P-SHAKE categorically unless UNSAFE is used, switch anyway if non-SETTLE groups with massless particles exist
  if ((shake_cnt.gt.0).AND.((cart_cons_method.eq.3).OR.(cart_cons_method.eq.2))) then
    write(ilog,*) 'Fatal. The P-SHAKE algorithm is not properly supported. You can skip this error using &
 &FMCSC_UNSAFE, but crashes are likely.'
    if (be_unsafe.EQV..false.) call fexit()
  end if
  do mt=1,shake_cnt
    i = nonsettlelst(mt)
    if ((cart_cons_method.eq.3).OR.((cart_cons_method.eq.2).AND.(constraints(i)%isrigid.EQV..true.))) then
      if (constraints(i)%hasmassless.EQV..true.) then
        write(ilog,*) 'Warning. Constraint groups with massless particles that are not &
 &covered by SETTLE are incompatible with P-SHAKE. Switching to standard SHAKE.'
        cart_cons_method = 1
        exit
      end if
    end if
    if (cart_cons_method.eq.1) exit
  end do
!
  effats = 0
  do i=1,cart_cons_grps
!   sort atomic indices and migrate massless atoms to end
!   this is because in the presence of multiple massless sites, their rebuilding may be hierarchical
    allocate(idxlst2(constraints(i)%nats,2))
    idxlst2(:,1) = constraints(i)%uidx(1:constraints(i)%nats)
    rs1 = 1
    rs2 = constraints(i)%nats
    ii = constraints(i)%nats
    foundit = .true.
    call merge_sort(ldim=ii,up=foundit,ilo=rs1,ihi=rs2,olist=idxlst2(:,2),list=idxlst2(:,1))
    rs1 = 0
    do k=1,ii
      if (mass(idxlst2(k,2)).gt.0.0) then
        rs1 = rs1 + 1
        constraints(i)%uidx(rs1) = idxlst2(k,2)
      end if
    end do
    do k=1,ii
      if (mass(idxlst2(k,2)).le.0.0) then
        effats = effats + 1
        rs1 = rs1 + 1
        constraints(i)%uidx(rs1) = idxlst2(k,2)
      end if
    end do
    deallocate(idxlst2)
    if (constraints(i)%hasmassless.EQV..true.) then
      if (constraints(i)%isrigid.EQV..false.) then
        write(ilog,*) 'Fatal. Even though it may in some cases be feasible, CAMPARI currently does not support &
 &massless (virtual) particles (sites) as part of constraint groups that are not completely rigid &
 &internally.'
        call fexit()
      end if
      k = constraints(i)%nats
      ii = 3
      do while (mass(constraints(i)%uidx(k)).le.0.0)
        if (atoms_areinset(iz(1:3,constraints(i)%uidx(k)),ii,constraints(i)%uidx(1:(k-1)),k-1).EQV..false.) then
          write(ilog,*) 'Fatal. The position of a massless (virtual) particle (site) does &
 &not depend solely and hierarchically on the position of atoms in the same constraint group.'
          call fexit()
        end if
        k = k - 1
      end do
    end if
  end do
!
! with massless constraints in generic SETTLE groups, we want to reorder the constraints as well
  do ii=1,settle_rest
    mt = settlelst3(ii)
    i = 1
    do while (i.le.constraints(mt)%nr3)
      if (constraints(mt)%massless(constraints(mt)%nr+i).EQV..false.) then
        constraints(mt)%nr = constraints(mt)%nr + 1
        k = constraints(mt)%nr
        constraints(mt)%idx(k,1) = constraints(mt)%idx3(i,1)
        constraints(mt)%idx(k,2) = constraints(mt)%idx3(i,3)
        constraints(mt)%eqd2(k) = getblen(constraints(mt)%idx3(i,1),constraints(mt)%idx3(i,3))**2
        foundit = constraints(mt)%massless(k+constraints(mt)%nr3-1)
        constraints(mt)%massless((k+1):(k+constraints(mt)%nr3-1)) = constraints(mt)%massless(k:(k+constraints(mt)%nr3-2))
        constraints(mt)%massless(k) = .false.
        if (i.lt.constraints(mt)%nr3) then
          constraints(mt)%idx3(i,:) = constraints(mt)%idx3(constraints(mt)%nr3,:)
          constraints(mt)%eqa3(i) = constraints(mt)%eqa3(constraints(mt)%nr3)
          constraints(mt)%massless(k+i) = foundit
        end if
        constraints(mt)%nr3 = constraints(mt)%nr3 - 1
      end if
      i = i + 1
    end do
    do i=1,constraints(mt)%nr
      if (constraints(mt)%massless(i).EQV..true.) then
        do k=constraints(mt)%nr,i+1,-1
          if (constraints(mt)%massless(k).EQV..false.) then
            rs1 = constraints(mt)%idx(k,1)
            rs2 = constraints(mt)%idx(k,2)
            tmpr = constraints(mt)%eqd2(k)
            constraints(mt)%idx(k,:) = constraints(mt)%idx(i,:)
            constraints(mt)%eqd2(k) = constraints(mt)%eqd2(i)
            constraints(mt)%idx(i,1) = rs1
            constraints(mt)%idx(i,2) = rs2
            constraints(mt)%eqd2(i) = tmpr
            constraints(mt)%massless(k) = .true.
            constraints(mt)%massless(i) = .false.
          end if
        end do
      end if
    end do
    if (constraints(mt)%eqd2(1).ne.constraints(mt)%eqd2(2)) then
      tmpr = 0.5*(sqrt(constraints(mt)%eqd2(1)) + sqrt(constraints(mt)%eqd2(2)))
      constraints(mt)%eqd2(1:2) = tmpr*tmpr
    end if
  end do
!
! if we have massless particles in unconstrained groups, exit as well (unless minimizer)
  rs1 = 0
  do i=1,n
    if (mass(i).le.0.0) rs1 = rs1 + 1
  end do
  if (effats.gt.rs1) then
    write(ilog,*) 'Fatal. Inconsistent setup of constraint groups. This is a bug.'
    call fexit()
  else if (effats.lt.rs1) then
    if ((dyn_mode.ne.6).OR.(mini_mode.eq.4)) then
      write(ilog,*) 'Fatal. Massless (virtual) particles (sites) are fatal if in Cartesian dynamics they &
 &are not rigidly connected to a constraint group.'
      call fexit()
    end if
  end if
  have_virtuals = .false.
  if (effats.gt.0) have_virtuals = .true.
!
! now reset targets to user choice
  if (cart_cons_source.le.2) then
    blix = MAXBOPAR + 1
    baix = MAXBAPAR + 1
    if (cart_cons_source.eq.2) blix = 2
    if (cart_cons_source.eq.2) baix = 2
!   note that the par_bl variables in question hold either the force field value (if available) or a 
!   value based on initial(!) geometry
    do mt2=1,shake_cnt+settle_rest
      if (mt2.le.shake_cnt) then
        mt = nonsettlelst(mt2)
      else
        mt = settlelst3(mt2-shake_cnt)
      end if
      if ((use_BOND(1).EQV..false.).AND.(constraints(mt)%nr.gt.0).AND.(cart_cons_source.eq.2)) then
        write(ilog,*) 'Fatal. Please enable bond length potentials in order to use force field-derived values &
 &for holonomic constraints.'
        call fexit()
      end if
      if ((use_BOND(2).EQV..false.).AND.(constraints(mt)%nr3.gt.0).AND.(cart_cons_source.eq.2)) then
        write(ilog,*) 'Fatal. Please enable bond angle potentials in order to use force field-derived values &
 &for holonomic (angle) constraints.'
        call fexit()
      end if
      kk = 0
      do i=1,constraints(mt)%nr
        foundit = .false.
        rs = min(atmres(constraints(mt)%idx(i,1)),atmres(constraints(mt)%idx(i,2)))
        do j=1,nrsbl(rs)
          if (((iaa(rs)%bl(j,1).eq.constraints(mt)%idx(i,1)).AND.(iaa(rs)%bl(j,2).eq.constraints(mt)%idx(i,2))).OR.&
 &            ((iaa(rs)%bl(j,2).eq.constraints(mt)%idx(i,1)).AND.(iaa(rs)%bl(j,1).eq.constraints(mt)%idx(i,2)))) then
            foundit = .true.
            constraints(mt)%eqd2(i) = iaa(rs)%par_bl(j,blix)*iaa(rs)%par_bl(j,blix)
          end if
        end do
        if ((foundit.EQV..false.).AND.(rs.ne.max(atmres(constraints(mt)%idx(i,1)),atmres(constraints(mt)%idx(i,2))))) then
          rs = max(atmres(constraints(mt)%idx(i,1)),atmres(constraints(mt)%idx(i,2)))
          do j=1,nrsbl(rs)
            if (((iaa(rs)%bl(j,1).eq.constraints(mt)%idx(i,1)).AND.(iaa(rs)%bl(j,2).eq.constraints(mt)%idx(i,2))).OR.&
 &              ((iaa(rs)%bl(j,2).eq.constraints(mt)%idx(i,1)).AND.(iaa(rs)%bl(j,1).eq.constraints(mt)%idx(i,2)))) then
              foundit = .true.
              constraints(mt)%eqd2(i) = iaa(rs)%par_bl(j,blix)*iaa(rs)%par_bl(j,blix)
            end if
          end do
        end if
        if (foundit.EQV..false.) then
          ii3 = 0
          do j=1,n13(constraints(mt)%idx(i,1))
            if (i13(j,constraints(mt)%idx(i,1)).eq.constraints(mt)%idx(i,2)) then
              do k=1,n12(constraints(mt)%idx(i,1))
                do kk=1,n12(constraints(mt)%idx(i,2))
                  if (i12(k,constraints(mt)%idx(i,1)).eq.i12(kk,constraints(mt)%idx(i,2))) then
                    ii3 = i12(k,constraints(mt)%idx(i,1))
                    exit
                  end if
                end do
                if (ii3.gt.0) exit
              end do
            end if
            if (ii3.gt.0) exit
          end do
          if (ii3.gt.0) then
            rslst(1) = atmres(constraints(mt)%idx(i,1))
            rslst(3) = atmres(constraints(mt)%idx(i,2))
            rslst(2) = atmres(ii3)
            do k=1,3
              do j=1,nrsbl(rslst(k))
                if (((iaa(rslst(k))%bl(j,1).eq.ii3).AND.(iaa(rslst(k))%bl(j,2).eq.constraints(mt)%idx(i,2))).OR.&
 &                  ((iaa(rslst(k))%bl(j,2).eq.ii3).AND.(iaa(rslst(k))%bl(j,1).eq.constraints(mt)%idx(i,2)))) then
                  l2 = iaa(rslst(k))%par_bl(j,blix)
                end if
                if (((iaa(rslst(k))%bl(j,1).eq.ii3).AND.(iaa(rslst(k))%bl(j,2).eq.constraints(mt)%idx(i,1))).OR.&
 &                  ((iaa(rslst(k))%bl(j,2).eq.ii3).AND.(iaa(rslst(k))%bl(j,1).eq.constraints(mt)%idx(i,1)))) then
                  l1 = iaa(rslst(k))%par_bl(j,blix)
                end if
              end do
            end do
            do k=1,3
              do j=1,nrsba(rslst(k))
                if (((iaa(rslst(k))%ba(j,1).eq.constraints(mt)%idx(i,1)).AND.&
 &                   (iaa(rslst(k))%ba(j,3).eq.constraints(mt)%idx(i,2)).AND.(iaa(rslst(k))%ba(j,2).eq.ii3)).OR.&
 &                  ((iaa(rslst(k))%ba(j,3).eq.constraints(mt)%idx(i,1)).AND.&
 &                   (iaa(rslst(k))%ba(j,1).eq.constraints(mt)%idx(i,2)).AND.(iaa(rslst(k))%ba(j,2).eq.ii3))) then
                  foundit = .true.
                  if ((iaa(rslst(k))%typ_ba(j).eq.2).AND.(cart_cons_source.eq.2)) then
                    constraints(mt)%eqd2(i) = iaa(rslst(k))%par_ba(j,4)*iaa(rslst(k))%par_ba(j,4)
                  else if ((iaa(rslst(k))%typ_ba(j).eq.1).OR.(iaa(rslst(k))%typ_ba(j).eq.3).OR.&
 &                         ((iaa(rslst(k))%typ_ba(j).eq.2).AND.(cart_cons_source.eq.1))) then              
                    constraints(mt)%eqd2(i) = l1*l1 + l2*l2 - 2*l1*l2*cos(iaa(rslst(k))%par_ba(j,baix)/RADIAN)
                  end if
                end if
                if (foundit.EQV..true.) exit
              end do
              if (foundit.EQV..true.) exit
            end do
          end if
        end if
      end do
!     easier for explicit angle constraints (we know they are topology-consistent)
      do i=1,constraints(mt)%nr3
        foundit = .false.
        rslst(1) = atmres(constraints(mt)%idx3(i,1))
        rslst(2) = atmres(constraints(mt)%idx3(i,2))
        rslst(3) = atmres(constraints(mt)%idx3(i,3))
        do k=1,3
          do j=1,nrsba(rslst(k))
            if (((iaa(rslst(k))%ba(j,1).eq.constraints(mt)%idx3(i,1)).AND.&
 &        (iaa(rslst(k))%ba(j,3).eq.constraints(mt)%idx3(i,3)).AND.(iaa(rslst(k))%ba(j,2).eq.constraints(mt)%idx3(i,2))).OR.&
 &              ((iaa(rslst(k))%ba(j,3).eq.constraints(mt)%idx3(i,1)).AND.&
 &        (iaa(rslst(k))%ba(j,1).eq.constraints(mt)%idx3(i,3)).AND.(iaa(rslst(k))%ba(j,2).eq.constraints(mt)%idx3(i,2)))) then
              foundit = .true.
              constraints(mt)%eqa3(i) = cos(iaa(rslst(k))%par_ba(j,baix)/RADIAN)
            end if
            if (foundit.EQV..true.) exit
          end do
          if (foundit.EQV..true.) exit
        end do
      end do
      if (mt2.gt.shake_cnt) then
        if (constraints(mt)%eqd2(1).ne.constraints(mt)%eqd2(2)) then
          tmpr = 0.5*(sqrt(constraints(mt)%eqd2(1)) + sqrt(constraints(mt)%eqd2(2)))
          constraints(mt)%eqd2(1:2) = tmpr*tmpr
        end if
      end if
    end do
  end if
!
! setup constraint coupling arrays (explicit coupling only)
  if (cart_cons_method.eq.4) then
    do mt2=1,shake_cnt
      mt = nonsettlelst(mt2)
      kk = 0
      do i=1,constraints(mt)%nr
        if (constraints(mt)%massless(i).EQV..true.) then
          write(ilog,*) 'Fatal. LINCS implementation currently does not support massless &
 &particles as part of LINCS constraint groups. Use SETTLE (if applicable) or SHAKE.'
          call fexit()
        end if
        ii1 = constraints(mt)%idx(i,1)
        ii2 = constraints(mt)%idx(i,2)
        do j=i+1,constraints(mt)%nr
          jj1 = constraints(mt)%idx(j,1)
          jj2 = constraints(mt)%idx(j,2)
          if ((ii1.eq.jj1).OR.(ii1.eq.jj2).OR.(ii2.eq.jj1).OR.(ii2.eq.jj2)) kk = kk + 1
        end do
      end do
      constraints(mt)%maxcpl = kk
    end do
!
    do mt2=1,shake_cnt
      mt = nonsettlelst(mt2)
      allocate(constraints(mt)%ncis(constraints(mt)%nr))
      constraints(mt)%ncis(:) = 0
      if (constraints(mt)%maxcpl.eq.0) cycle
!
      allocate(constraints(mt)%cidx(constraints(mt)%nr,constraints(mt)%maxcpl))
      allocate(constraints(mt)%ccoeff(constraints(mt)%nr,constraints(mt)%maxcpl))
      do i=1,constraints(mt)%nr
        ii1 = constraints(mt)%idx(i,1)
        ii2 = constraints(mt)%idx(i,2)
        do j=1,constraints(mt)%nr
          if (i.eq.j) cycle
          jj1 = constraints(mt)%idx(j,1)
          jj2 = constraints(mt)%idx(j,2)
          if (ii1.eq.jj1) then
            constraints(mt)%ncis(i) = constraints(mt)%ncis(i) + 1
            constraints(mt)%cidx(i,constraints(mt)%ncis(i)) = j
            constraints(mt)%ccoeff(i,constraints(mt)%ncis(i)) = (1.0/sqrt(1.0/mass(ii1) + 1.0/mass(ii2)))*&
 &(1.0/sqrt(1.0/mass(jj1) + 1.0/mass(jj2)))*(1.0/mass(ii1))*(-1.0)
          else if (ii1.eq.jj2) then
            constraints(mt)%ncis(i) = constraints(mt)%ncis(i) + 1
            constraints(mt)%cidx(i,constraints(mt)%ncis(i)) = j
            constraints(mt)%ccoeff(i,constraints(mt)%ncis(i)) = (1.0/sqrt(1.0/mass(ii1) + 1.0/mass(ii2)))*&
 &(1.0/sqrt(1.0/mass(jj1) + 1.0/mass(jj2)))*(1.0/mass(ii1))*(1.0)
          else if (ii2.eq.jj1) then
            constraints(mt)%ncis(i) = constraints(mt)%ncis(i) + 1
            constraints(mt)%cidx(i,constraints(mt)%ncis(i)) = j
            constraints(mt)%ccoeff(i,constraints(mt)%ncis(i)) = (1.0/sqrt(1.0/mass(ii1) + 1.0/mass(ii2)))*&
 &(1.0/sqrt(1.0/mass(jj1) + 1.0/mass(jj2)))*(1.0/mass(ii2))*(1.0)
          else if (ii2.eq.jj2) then
            constraints(mt)%ncis(i) = constraints(mt)%ncis(i) + 1
            constraints(mt)%cidx(i,constraints(mt)%ncis(i)) = j
            constraints(mt)%ccoeff(i,constraints(mt)%ncis(i)) = (1.0/sqrt(1.0/mass(ii1) + 1.0/mass(ii2)))*&
 &(1.0/sqrt(1.0/mass(jj1) + 1.0/mass(jj2)))*(1.0/mass(ii2))*(-1.0)
          end if
        end do
      end do
    end do
  end if
!
! disable mixed mode if pointless
  foundit = .false.
  if (cart_cons_method.eq.2) then
    do i=1,shake_cnt
      mt = nonsettlelst(i)
      if (constraints(mt)%isrigid.EQV..true.) then
        foundit = .true.
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      write(ilog,*) 'Warning. There are no internally rigid constraint groups to handle for P-SHAKE. Switching to &
 &standard SHAKE.'
      cart_cons_method = 1
    end if
  end if
!
  ii = 0
  do mt=1,cart_cons_grps
    ii = ii + constraints(mt)%nr + constraints(mt)%nr3 + constraints(mt)%nr4
  end do
  if (ii.eq.0) then
    no_shake = .true.
    return
  end if
!
  do i=1,shake_cnt
    mt = nonsettlelst(i)
    if ((constraints(mt)%nr + constraints(mt)%nr3 + constraints(mt)%nr4).gt.constraints(mt)%nats-1) then
      write(ilog,*) 'Warning. At least one constraint group has more holonomic constraints than&
 & the minimum possible number of bonds within the group. This MAY lead to unexpected crashes &
 &and generally poor convergence of constraint solvers. This warning may occur repeatedly &
 &for systems containing rings or if higher-order constraints are requested.'
      exit
    end if
  end do
#ifdef ENABLE_THREADS
! swap some constraints for potentially better partitioning (SHAKE with distances only)
  if (cart_cons_method.eq.1) then
    do i=1,cart_cons_grps
      k = 1
      do while (k.lt.constraints(i)%nr)
        ii = constraints(i)%idx(k,1)
        kk = constraints(i)%idx(k,2)
        if (atmres(ii).ne.atmres(kk)) then
          j = k + 1
          do while ((atmres(ii).eq.atmres(constraints(i)%idx(j,1))).AND.(atmres(ii).eq.atmres(constraints(i)%idx(j,2))))
            j = j + 1
          end do
          if ((atmres(kk).eq.atmres(constraints(i)%idx(j,1))).AND.(atmres(kk).eq.atmres(constraints(i)%idx(j,2))).AND.&
   &          (j.ne.(k+1))) then
!           swap j-1 and k
            tmpr = constraints(i)%eqd2(k)
            ii3 = constraints(i)%idx(k,1)
            ii4 = constraints(i)%idx(k,2)
            foundit = constraints(i)%massless(k)
            constraints(i)%eqd2(k) = constraints(i)%eqd2(j-1)
            constraints(i)%idx(k,:) = constraints(i)%idx(j-1,:)
            constraints(i)%massless(k) = constraints(i)%massless(j-1)
            constraints(i)%eqd2(j-1) = tmpr
            constraints(i)%idx(j-1,1) = ii3
            constraints(i)%idx(j-1,2) = ii4
            constraints(i)%massless(j-1) = foundit
          end if
          k = j
        else
          k = k + 1
        end if
      end do
    end do
  end if
  call threads_setup_shake()
#endif
!  do i=1,cart_cons_grps
!    write(*,*) 'IDX'
!    write(*,*) constraints(i)%idx(1:constraints(i)%nr,1)
!    write(*,*) constraints(i)%idx(1:constraints(i)%nr,2)
!    write(*,*) constraints(i)%massless(1:constraints(i)%nr)
!    write(*,*) sqrt(constraints(i)%eqd2(1:constraints(i)%nr))
!    write(*,*) constraints(i)%idx3(1:constraints(i)%nr3,1)
!    write(*,*) constraints(i)%idx3(1:constraints(i)%nr3,2)
!    write(*,*) constraints(i)%idx3(1:constraints(i)%nr3,3)
!    write(*,*) RADIAN*acos(constraints(i)%eqa3(1:constraints(i)%nr3))
!    write(*,*) constraints(i)%massless((constraints(i)%nr+1):(constraints(i)%nr+constraints(i)%nr3))
!    write(*,*) constraints(i)%uidx(1:constraints(i)%nats)
!  end do
!
#ifdef LINK_LAPACK
  if (cart_cons_method.eq.3) then
    do i=1,shake_cnt
      mt = nonsettlelst(i)
      call shake_p_setA(mt)
    end do
  else if (cart_cons_method.eq.2) then
    do i=1,shake_cnt
      mt = nonsettlelst(i)
      if (constraints(mt)%isrigid.EQV..true.) call shake_p_setA(mt)
    end do
  end if
#endif
!
! finally, adjust constraints (this changes atomic positions!)
  if (cart_cons_source.le.2) then
    xref(1:n) = x(1:n)
    yref(1:n) = y(1:n)
    zref(1:n) = z(1:n)
    if ((pdb_analyze.EQV..false.).AND.(do_restart.EQV..false.).AND.(fycxyz.eq.2).AND.(use_dyn.EQV..true.)) then
      call shake_wrap(azero)
      do imol=1,nmol
        call genzmat(imol)
      end do
      call zmatfyc2()
      stringy = 'START '
      call FMCSC_dump(stringy)
    end if
  end if
!
end
!
! this small helper functon finds (inefficiently) whether all entries of testl are also part of refl
!
function atoms_areinset(testl,tln,refl,rln)
!
  integer i,k,tln,rln,testl(tln),refl(rln)
  logical foundit,atoms_areinset
!
  atoms_areinset = .true.
  do i=1,tln
    foundit = .false.
    do k=1,rln
      if (testl(i).eq.refl(k)) then
        foundit = .true.
        exit
      end if
    end do
    if (foundit.EQV..false.) then
      atoms_areinset = .false.
      exit
    end if
  end do
!
  return
!
end
!
!------------------------------------------------------------------------------------------------------
!
! the wrapper for constraint solvers when called from dynamics - could be made more efficient
!
! it can be parallelized in a thread-safe manner by distributing calls to any of the specific functions
! this is because every atom can at most be part of 1 constraint group and because functions operate
! strictly on groups 
!
subroutine shake_wrap(tpi)
!
  use shakeetal
  use iounit 
  use mcsums
  use molecule
  use sequen
  use atoms
  use math
  use zmatrix
#ifdef ENABLE_THREADS
  use threads
  use cutoffs
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
  integer i,k
  RTYPE mvals1(2),mvals2(2),mvals4(2),mvals5(2),mvals6(2),rvals(4),rvals2(4),mvals3(2),rvals3(4)
  RTYPE rvals1(4),rvals4(4),rvals5(4),rvals6(4)
  integer(KIND=8) t1,t2
#ifdef ENABLE_THREADS
  integer sta,sto,ixx
  integer(KIND=8) ttimer
  logical OMP_IN_PARALLEL
#endif
!
#ifdef ENABLE_THREADS
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, shake_wrap(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    atinfo(1,thr_limits(1,tpi):thr_limits(2,tpi)) = x(thr_limits(1,tpi):thr_limits(2,tpi))
    atinfo(2,thr_limits(1,tpi):thr_limits(2,tpi)) = y(thr_limits(1,tpi):thr_limits(2,tpi))
    atinfo(3,thr_limits(1,tpi):thr_limits(2,tpi)) = z(thr_limits(1,tpi):thr_limits(2,tpi))
!$OMP BARRIER
    if (thr_dlb(4,1).gt.0) then
      if (tpi.eq.1) thr_dlb(4,2) = thr_dlb(4,2) + 1
      call System_Clock(count=t1)
      thr_timings(9,tpi) = thr_timings(9,tpi) + t1
    end if
  end if
  if (tpi.le.1) call System_Clock(count=t1)
#else
  call System_Clock(t1)
#endif
!
! first get SETTLE out of the way
! WARNING: using hard-coded parameters for TipNp and SPC geometry
  if ((settle_tip3ps.gt.0).OR.(settle_tip4ps.gt.0).OR.(settle_tip5ps.gt.0).OR.(settle_tip4pes.gt.0)) then
    rvals(3) = 0.5*0.9572*sin(104.52/RADIAN)/sin(0.5*(180.0-104.52)/RADIAN) ! 1.513901
    rvals(1) = 2.0*sqrt(0.9572*0.9572 - rvals(3)*rvals(3))
    rvals(2) = sqrt(0.9572*0.9572 - rvals(3)*rvals(3))
    rvals(4) = 2.0*rvals(3)
  end if
  if (settle_tip3ps.gt.0) then
    mvals1(1:2) = mass(constraints(settlelst1(1))%uidx(1:2))/molmass(moltypid(molofrs(atmres(constraints(settlelst1(1))%uidx(1)))))
    rvals1(3:4) = rvals(3:4)
    rvals1(1) = mvals1(2)*rvals(1)
    rvals1(2) = rvals(2) - rvals1(1)
  end if
  if (settle_tip4ps.gt.0) then
    mvals4(1:2) = mass(constraints(settlelst4(1))%uidx(1:2))/molmass(moltypid(molofrs(atmres(constraints(settlelst4(1))%uidx(1)))))
    rvals4(3:4) = rvals(3:4)
    rvals4(1) = mvals4(2)*rvals(1)
    rvals4(2) = rvals(2) - rvals4(1)
  end if
  if (settle_tip4pes.gt.0) then
    mvals6(1:2) = mass(constraints(settlelst6(1))%uidx(1:2))/molmass(moltypid(molofrs(atmres(constraints(settlelst6(1))%uidx(1)))))
    rvals6(3:4) = rvals(3:4)
    rvals6(1) = mvals6(2)*rvals(1)
    rvals6(2) = rvals(2) - rvals6(1)
  end if
  if (settle_tip5ps.gt.0) then
    mvals5(1:2) = mass(constraints(settlelst5(1))%uidx(1:2))/molmass(moltypid(molofrs(atmres(constraints(settlelst5(1))%uidx(1)))))
    rvals5(3:4) = rvals(3:4)
    rvals5(1) = mvals5(2)*rvals(1)
    rvals5(2) = rvals(2) - rvals5(1)
  end if

  if (settle_spcs.gt.0) then
    mvals2(1:2) = mass(constraints(settlelst2(1))%uidx(1:2))/molmass(moltypid(molofrs(atmres(constraints(settlelst2(1))%uidx(1)))))
    rvals2(3) = 0.5*1.0*sin(109.47/RADIAN)/sin(0.5*(180.0-109.47)/RADIAN) ! 1.632981
    rvals2(1) = 2.0*mvals2(2)*sqrt(1.0 - rvals2(3)*rvals2(3))
    rvals2(2) = sqrt(1.0 - rvals2(3)*rvals2(3)) - rvals2(1)
    rvals2(4) = 2.0*rvals2(3)
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = settle_tip3ps
  else
    sta = thr_limits(11,tpi)
    sto = thr_limits(12,tpi)
  end if
  do i=sta,sto
#else
  do i=1,settle_tip3ps
#endif
    call settle(settlelst1(i),mvals1,rvals1)
  end do
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = settle_spcs
  else
    sta = thr_limits(13,tpi)
    sto = thr_limits(14,tpi)
  end if
  do i=sta,sto
#else
  do i=1,settle_spcs
#endif
    call settle(settlelst2(i),mvals2,rvals2)
  end do
! for the 4+ site models, the virtual sites' Z-matrix entries are reset
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = settle_tip4ps
  else
    sta = thr_limits(15,tpi)
    sto = thr_limits(16,tpi)
  end if
  do i=sta,sto
#else
  do i=1,settle_tip4ps
#endif
    blen(constraints(settlelst4(i))%uidx(4)) = 0.15
    bang(constraints(settlelst4(i))%uidx(4)) = 52.26
    ztor(constraints(settlelst4(i))%uidx(4)) = 0.0
    call settle(settlelst4(i),mvals4,rvals4)
  end do
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = settle_tip5ps
  else
    sta = thr_limits(17,tpi)
    sto = thr_limits(18,tpi)
  end if
  do i=sta,sto
#else
  do i=1,settle_tip5ps
#endif
    blen(constraints(settlelst5(i))%uidx(4:5)) = 0.7
    bang(constraints(settlelst5(i))%uidx(4:5)) = 110.69477
    ztor(constraints(settlelst5(i))%uidx(4:5)) = 110.69477
    call settle(settlelst5(i),mvals5,rvals5)
  end do
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = settle_tip4pes
  else
    sta = thr_limits(19,tpi)
    sto = thr_limits(20,tpi)
  end if
  do i=sta,sto
#else
  do i=1,settle_tip4pes
#endif
    blen(constraints(settlelst6(i))%uidx(4)) = 0.125
    bang(constraints(settlelst6(i))%uidx(4)) = 52.26
    ztor(constraints(settlelst6(i))%uidx(4)) = 0.0
    call settle(settlelst6(i),mvals6,rvals6)
  end do
! for the remaining groups, take values from constraint settings, and leave Z-matrix entries
! of eventual virtual sites untouched (this also happens in SHAKE)
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = settle_rest
  else
    sta = thr_limits(21,tpi)
    sto = thr_limits(22,tpi)
  end if
  do i=sta,sto
#else
  do i=1,settle_rest
#endif
    mvals3(1:2) = mass(constraints(settlelst3(i))%uidx(1:2))/sum(mass(constraints(settlelst3(i))%uidx(1:3)))
    rvals3(4) = sqrt(constraints(settlelst3(i))%eqd2(3))
    rvals3(3) = 0.5*rvals3(4)
    rvals3(1) = 2.0*mvals3(2)*sqrt(constraints(settlelst3(i))%eqd2(1) - rvals3(3)*rvals3(3))
    rvals3(2) = sqrt(constraints(settlelst3(i))%eqd2(1) - rvals3(3)*rvals3(3)) - rvals3(1)
    call settle(settlelst3(i),mvals3,rvals3)
  end do
!
! now deal with remaining constraint groups
  if (cart_cons_method.eq.1) then
#ifdef ENABLE_THREADS
    if (tpi.eq.0) then
      sta = 1
      sto = shake_cnt
    else
      sta = thr_limits(23,tpi)
      sto = thr_limits(24,tpi)
      ixx = 1
      if (nccgs.gt.0) then
        do while ((ccg_limits(ixx,6,tpi).lt.sta).AND.(ixx.lt.nccgs))
          ixx = min(ixx+1,nccgs)
        end do
      end if
    end if
    do k=sta,sto
      if ((nccgs.gt.0).AND.(tpi.gt.0)) then
        if (ccg_limits(ixx,6,tpi).eq.k) then
          ixx = min(ixx+1,nccgs)
          cycle
        end if
      end if
#else
    do k=1,shake_cnt
#endif
      i = nonsettlelst(k)
      if ((constraints(i)%nr3.eq.0).AND.(constraints(i)%nr4.eq.0)) then
        call shake_rcb_just2(i)
      else if (constraints(i)%nr4.eq.0) then
        call shake_rcb_just23(i)
      else
        call shake_rcb(i)
      end if
    end do
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      if (thr_dlb(4,1).gt.0) then
        call System_Clock(count=ttimer)
        thr_timings(10,tpi) = thr_timings(10,tpi) + ttimer
      end if
      do k=1,nccgs
        i = ccg_limits(k,5,tpi)
        call shake_rcb_just2_threads(i,tpi,k)
      end do
    end if
#endif
  else if (cart_cons_method.eq.2) then
#ifdef ENABLE_THREADS
    if (tpi.eq.0) then
      sta = 1
      sto = shake_cnt
    else
      sta = thr_limits(23,tpi)
      sto = thr_limits(24,tpi)
    end if
    do k=sta,sto
#else
    do k=1,shake_cnt
#endif
      i = nonsettlelst(k)
      if (constraints(i)%isrigid.EQV..true.) then
        call shake_p(i)
      else
        if (constraints(i)%nr3.eq.0) then
          call shake_rcb_just2(i)
        else if (constraints(i)%nr4.eq.0) then
          call shake_rcb_just23(i)
        else
          call shake_rcb(i)
        end if
      end if
    end do
  else if (cart_cons_method.eq.3) then
#ifdef ENABLE_THREADS
    if (tpi.eq.0) then
      sta = 1
      sto = shake_cnt
    else
      sta = thr_limits(23,tpi)
      sto = thr_limits(24,tpi)
    end if
    do k=sta,sto
#else
    do k=1,shake_cnt
#endif
      i = nonsettlelst(k)
      call shake_p(i)
    end do
#ifdef LINK_LAPACK
    if (mod(nstep,100).eq.0) then
      do k=1,shake_cnt
        i = nonsettlelst(k)
        if (constraints(i)%isrigid.EQV..false.) then
          call shake_p_setA(i)
        end if
      end do
    end if
#endif
  else if (cart_cons_method.eq.4) then
#ifdef ENABLE_THREADS
    if (tpi.eq.0) then
      sta = 1
      sto = shake_cnt
    else
      sta = thr_limits(23,tpi)
      sto = thr_limits(24,tpi)
    end if
    do k=sta,sto
#else
    do k=1,shake_cnt
#endif
      i = nonsettlelst(k)
      if (constraints(i)%maxcpl.gt.0) then
        call lincs(i)
      else
        call shake_rcb_just2(i)
      end if
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported method for holonomic constraints in &
 &shake_wrap(). This is a bug.'
    call fexit()
  end if
!
#ifdef ENABLE_THREADS
  if ((tpi.gt.0).AND.(cart_cons_method.ne.1)) then
    if (thr_dlb(4,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(10,tpi) = thr_timings(10,tpi) + ttimer
    end if
  end if
  if (tpi.eq.1) then 
    call System_Clock(t2)
    time_holo = time_holo + 1.0*(t2 - t1)
  end if
#else
  call System_Clock(t2)
  time_holo = time_holo + 1.0*(t2 - t1)
#endif
!
end
!
!---------------------------------------------------------------------------------------------------------
!
subroutine shake_rcb(grpid)
!
  use atoms
  use shakeetal
  use system
  use molecule
  use iounit
  use mcsums
  use forces
  use math
  use zmatrix
!
  implicit none
! 
  integer grpid,nc,k,iter
  RTYPE dt2,maxdev,idt,maxdev3,maxdev4,dumm
  RTYPE dvec(3),d2,xvt,dveca1(3),dveca2(3),da1,da2,denomt,dotpt
  RTYPE vv1(3),vv2(3),vv3(3),vv4(3),ctm,stm,smin,cmin,ct2,st2,si
  RTYPE dvecref(constraints(grpid)%nr,3)
  RTYPE vvec(constraints(grpid)%nr,3)
  RTYPE vvec1(constraints(grpid)%nr,3)
  RTYPE vvec2(constraints(grpid)%nr,3)
  RTYPE vveca1(constraints(grpid)%nr3,3)
  RTYPE vveca2(constraints(grpid)%nr3,3)
  RTYPE vveca3(constraints(grpid)%nr3,3)
  RTYPE vveci1(3,constraints(grpid)%nr4)
  RTYPE vveci2(3,constraints(grpid)%nr4)
  RTYPE vveci3(3,constraints(grpid)%nr4)
  RTYPE vveci4(3,constraints(grpid)%nr4)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
  RTYPE mterm(constraints(grpid)%nr3,3)
  RTYPE mterm4(constraints(grpid)%nr4,4)
  RTYPE dvecra1(constraints(grpid)%nr3,3)
  RTYPE dvecra2(constraints(grpid)%nr3,3)
  RTYPE dra1(constraints(grpid)%nr3)
  RTYPE dra2(constraints(grpid)%nr3)
  RTYPE dotp(constraints(grpid)%nr3)
  RTYPE denom(constraints(grpid)%nr3)
  RTYPE lambdas(constraints(grpid)%nr)
  RTYPE lambdas3(constraints(grpid)%nr3)
  RTYPE lambdas4(constraints(grpid)%nr4)
  RTYPE deviats(constraints(grpid)%nr)
  RTYPE deviats3(constraints(grpid)%nr3)
  RTYPE deviats4(constraints(grpid)%nr4)
!
  dt2 = dyn_dt*dyn_dt
  idt = 1.0/dyn_dt
  do nc=1,constraints(grpid)%nr
    if (constraints(grpid)%massless(nc).EQV..true.) then
      mterm1(nc) = 0.0
      mterm2(nc) = 0.0
      cycle
    end if
    mterm1(nc) = dt2/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = dt2/mass(constraints(grpid)%idx(nc,2))
    call dis_bound3r(constraints(grpid)%idx(nc,2),constraints(grpid)%idx(nc,1),dvec(:),dumm)
    dvecref(nc,:) = dvec(:)
  end do
  do k=1,3
    vvec1(:,k) = dvecref(:,k)*mterm1(:)
    vvec2(:,k) = -dvecref(:,k)*mterm2(:)
    vvec(:,k) = vvec1(:,k) - vvec2(:,k)
  end do
!
  do nc=1,constraints(grpid)%nr3
    if (constraints(grpid)%massless(nc+constraints(grpid)%nr).EQV..true.) then
      mterm(nc,1) = 0.0
      mterm(nc,2) = 0.0
      mterm(nc,3) = 0.0
      cycle
    end if
    mterm(nc,1) = 0.5*dt2/mass(constraints(grpid)%idx3(nc,1))
    mterm(nc,2) = 0.5*dt2/mass(constraints(grpid)%idx3(nc,2))
    mterm(nc,3) = 0.5*dt2/mass(constraints(grpid)%idx3(nc,3))
    call dis_bound3r(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,1),dvec(:),dumm)
    dvecra1(nc,:) = dvec(:)
    call dis_bound3r(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,3),dvec(:),dumm)
    dvecra2(nc,:) = dvec(:)
  end do
  if (constraints(grpid)%nr3.gt.0) then
    dra1(:) = dvecra1(:,1)**2 + dvecra1(:,2)**2 + dvecra1(:,3)**2
    dra2(:) = dvecra2(:,1)**2 + dvecra2(:,2)**2 + dvecra2(:,3)**2
    denom(:) = 1.0/(dra1(:)*dra2(:))
    dra1(:) = sqrt(dra1(:))
    dra2(:) = sqrt(dra2(:))
    dotp(:) = dvecra1(:,1)*dvecra2(:,1) + dvecra1(:,2)*dvecra2(:,2) + dvecra1(:,3)*dvecra2(:,3)
    do k=1,3
      vveca1(:,k) = denom(:)*(dra1(:)*dra2(:)*dvecra2(:,k) - dotp(:)*dra2(:)*dvecra1(:,k)/dra1(:))
      vveca3(:,k) = denom(:)*(dra1(:)*dra2(:)*dvecra1(:,k) - dotp(:)*dra1(:)*dvecra2(:,k)/dra2(:))
      vveca2(:,k) = - (vveca1(:,k) + vveca3(:,k))
      vveca1(:,k) = vveca1(:,k)*mterm(:,1)
      vveca2(:,k) = vveca2(:,k)*mterm(:,2)
      vveca3(:,k) = vveca3(:,k)*mterm(:,3)
    end do
  end if
!
  do nc=1,constraints(grpid)%nr4
    if (constraints(grpid)%massless(nc+constraints(grpid)%nr+constraints(grpid)%nr3).EQV..true.) then
      mterm4(nc,1) = 0.0
      mterm4(nc,2) = 0.0
      mterm4(nc,3) = 0.0
      mterm4(nc,4) = 0.0
      cycle
    end if
    mterm4(nc,1) = 0.5*dt2/mass(constraints(grpid)%idx4(nc,1))
    mterm4(nc,2) = 0.5*dt2/mass(constraints(grpid)%idx4(nc,2))
    mterm4(nc,3) = 0.5*dt2/mass(constraints(grpid)%idx4(nc,3))
    mterm4(nc,4) = 0.5*dt2/mass(constraints(grpid)%idx4(nc,4))
    call onetor_deriv_ref(constraints(grpid)%idx4(nc,1),constraints(grpid)%idx4(nc,2),&
 &constraints(grpid)%idx4(nc,3),constraints(grpid)%idx4(nc,4),ctm,stm,&
 &vveci1(:,nc),vveci2(:,nc),vveci3(:,nc),vveci4(:,nc))
    vveci1(:,nc) = vveci1(:,nc)*mterm4(nc,1)
    vveci2(:,nc) = vveci2(:,nc)*mterm4(nc,2)
    vveci3(:,nc) = vveci3(:,nc)*mterm4(nc,3)
    vveci4(:,nc) = vveci4(:,nc)*mterm4(nc,4)
  end do
!
  iter = 0
  maxdev = 2.0*(shake_tol) + 1.0
  maxdev3 = 2.0*(shake_atol) + 1.0
  maxdev4 = 2.0*(shake_tol) + 1.0
!
  do while ((maxdev.gt.shake_tol).OR.(maxdev3.gt.shake_atol).OR.(maxdev4.gt.shake_tol))
    iter = iter + 1 
    do nc=1,constraints(grpid)%nr
      if (constraints(grpid)%massless(nc).EQV..true.) cycle
      call dis_bound3(constraints(grpid)%idx(nc,2),constraints(grpid)%idx(nc,1),dvec(:),d2)
!
      xvt = dvec(1)*vvec(nc,1) + dvec(2)*vvec(nc,2) + dvec(3)*vvec(nc,3)
      if (abs(xvt).lt.(dt2*1.0e-8)) then
        deviats(nc) = 0.0
        lambdas(nc) = 0.0
      else
        xvt = 0.5/xvt
        deviats(nc) = constraints(grpid)%eqd2(nc) - d2
        lambdas(nc) = deviats(nc)*xvt
      end if
!
      x(constraints(grpid)%idx(nc,1)) = x(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,1)
      x(constraints(grpid)%idx(nc,2)) = x(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,1)
      y(constraints(grpid)%idx(nc,1)) = y(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,2)
      y(constraints(grpid)%idx(nc,2)) = y(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,2)
      z(constraints(grpid)%idx(nc,1)) = z(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,3)
      z(constraints(grpid)%idx(nc,2)) = z(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,3)
    end do
!
    do nc=1,constraints(grpid)%nr3
      if (constraints(grpid)%massless(constraints(grpid)%nr+nc).EQV..true.) cycle
      call dis_bound3(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,1),dveca1(:),da1)
      call dis_bound3(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,3),dveca2(:),da2)
      denomt = 1.0/(da1*da2)
      da1 = sqrt(da1)
      da2 = sqrt(da2)
      dotpt = dveca1(1)*dveca2(1) + dveca1(2)*dveca2(2) + dveca1(3)*dveca2(3)
      vv1(:) = denomt*(da1*da2*dveca2(:) - dotpt*da2*dveca1(:)/da1)
      vv3(:) = denomt*(da1*da2*dveca1(:) - dotpt*da1*dveca2(:)/da2)
      vv2(:) = - (vv1(:) + vv3(:))
      xvt = dot_product(vv1(:),vveca1(nc,:)) + dot_product(vv2(:),vveca2(nc,:)) + dot_product(vv3(:),vveca3(nc,:))
!      if (abs(xvt).lt.(dt2*1.0e-8)) then
!        deviats3(nc) = 0.0
!        lambdas3(nc) = 0.0
!      else
        xvt = 1.0/xvt
        deviats3(nc) = constraints(grpid)%eqa3(nc) - dotpt*sqrt(denomt)
        lambdas3(nc) = deviats3(nc)*xvt
!      end if
!
      x(constraints(grpid)%idx3(nc,1)) = x(constraints(grpid)%idx3(nc,1)) + lambdas3(nc)*vveca1(nc,1)
      x(constraints(grpid)%idx3(nc,2)) = x(constraints(grpid)%idx3(nc,2)) + lambdas3(nc)*vveca2(nc,1)
      x(constraints(grpid)%idx3(nc,3)) = x(constraints(grpid)%idx3(nc,3)) + lambdas3(nc)*vveca3(nc,1)
      y(constraints(grpid)%idx3(nc,1)) = y(constraints(grpid)%idx3(nc,1)) + lambdas3(nc)*vveca1(nc,2)
      y(constraints(grpid)%idx3(nc,2)) = y(constraints(grpid)%idx3(nc,2)) + lambdas3(nc)*vveca2(nc,2)
      y(constraints(grpid)%idx3(nc,3)) = y(constraints(grpid)%idx3(nc,3)) + lambdas3(nc)*vveca3(nc,2)
      z(constraints(grpid)%idx3(nc,1)) = z(constraints(grpid)%idx3(nc,1)) + lambdas3(nc)*vveca1(nc,3)
      z(constraints(grpid)%idx3(nc,2)) = z(constraints(grpid)%idx3(nc,2)) + lambdas3(nc)*vveca2(nc,3)
      z(constraints(grpid)%idx3(nc,3)) = z(constraints(grpid)%idx3(nc,3)) + lambdas3(nc)*vveca3(nc,3)
    end do
!
    do nc=1,constraints(grpid)%nr4
      if (constraints(grpid)%massless(constraints(grpid)%nr+constraints(grpid)%nr3+nc).EQV..true.) cycle
      call onetor_deriv(constraints(grpid)%idx4(nc,1),constraints(grpid)%idx4(nc,2),&
 &constraints(grpid)%idx4(nc,3),constraints(grpid)%idx4(nc,4),ctm,stm,vv1(:),vv2(:),vv3(:),vv4(:))
      cmin = cos(constraints(grpid)%eqi(nc)/RADIAN) 
      smin = sin(constraints(grpid)%eqi(nc)/RADIAN)
!     get cos(phi-phi_0) and sin(phi-phi_0)
      ct2 = cmin*ctm + smin*stm
      st2 = cmin*stm - smin*ctm
      if (st2.lt.0.0) then
        si = -1.0
      else
        si = 1.0
      end if
      if (abs(ct2).gt.0.1) then ! safe regime for asin
        deviats4(nc) = RADIAN*asin(st2)
        if (ct2.lt.0.0) deviats4(nc) = 180.0 - deviats4(nc) ! safe for asin but angle very far away 
      else ! switch to cosine (means angle is far away from target)
        deviats4(nc) = si*RADIAN*acos(ct2)
      end if
      if (deviats4(nc).lt.-180.0) deviats4(nc) = deviats4(nc) + 360.0
      if (deviats4(nc).gt.180.0) deviats4(nc) = deviats4(nc) - 360.0
      deviats4(nc) = -deviats4(nc)
      xvt = dot_product(vv1(:),vveci1(:,nc)) + dot_product(vv2(:),vveci2(:,nc)) + &
 &          dot_product(vv3(:),vveci3(:,nc)) + dot_product(vv4(:),vveci4(:,nc))
      xvt = 1.0/xvt
      lambdas4(nc) = deviats4(nc)*xvt
!
      x(constraints(grpid)%idx4(nc,1)) = x(constraints(grpid)%idx4(nc,1)) + lambdas4(nc)*vveci1(1,nc)
      x(constraints(grpid)%idx4(nc,2)) = x(constraints(grpid)%idx4(nc,2)) + lambdas4(nc)*vveci2(1,nc)
      x(constraints(grpid)%idx4(nc,3)) = x(constraints(grpid)%idx4(nc,3)) + lambdas4(nc)*vveci3(1,nc)
      x(constraints(grpid)%idx4(nc,4)) = x(constraints(grpid)%idx4(nc,4)) + lambdas4(nc)*vveci4(1,nc)
      y(constraints(grpid)%idx4(nc,1)) = y(constraints(grpid)%idx4(nc,1)) + lambdas4(nc)*vveci1(2,nc)
      y(constraints(grpid)%idx4(nc,2)) = y(constraints(grpid)%idx4(nc,2)) + lambdas4(nc)*vveci2(2,nc)
      y(constraints(grpid)%idx4(nc,3)) = y(constraints(grpid)%idx4(nc,3)) + lambdas4(nc)*vveci3(2,nc)
      y(constraints(grpid)%idx4(nc,4)) = y(constraints(grpid)%idx4(nc,4)) + lambdas4(nc)*vveci4(2,nc)
      z(constraints(grpid)%idx4(nc,1)) = z(constraints(grpid)%idx4(nc,1)) + lambdas4(nc)*vveci1(3,nc)
      z(constraints(grpid)%idx4(nc,2)) = z(constraints(grpid)%idx4(nc,2)) + lambdas4(nc)*vveci2(3,nc)
      z(constraints(grpid)%idx4(nc,3)) = z(constraints(grpid)%idx4(nc,3)) + lambdas4(nc)*vveci3(3,nc)
      z(constraints(grpid)%idx4(nc,4)) = z(constraints(grpid)%idx4(nc,4)) + lambdas4(nc)*vveci4(3,nc)
    end do
!
    maxdev = 0.0
    maxdev3 = 0.0
    maxdev4 = 0.0
    do nc=1,constraints(grpid)%nr
      if (constraints(grpid)%massless(nc).EQV..true.) cycle
      if (abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc)).gt.maxdev) then
        maxdev = abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc))
      end if
    end do
    do nc=1,constraints(grpid)%nr3
      if (constraints(grpid)%massless(constraints(grpid)%nr+nc).EQV..true.) cycle
      if (abs(deviats3(nc)).gt.maxdev3) then
        maxdev3 = abs(deviats3(nc))
      end if
    end do
    do nc=1,constraints(grpid)%nr4
      if (constraints(grpid)%massless(constraints(grpid)%nr+constraints(grpid)%nr3+nc).EQV..true.) cycle
      if (abs(deviats4(nc)).gt.maxdev4) then
        maxdev4 = abs(deviats4(nc))
      end if
    end do
!
    if (iter.gt.shake_maxiter) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_WARNINGS)
#endif
      cs_wrncnt(1) = cs_wrncnt(1) + 1
      if (cs_wrncnt(1).eq.cs_wrnlmt(1)) then
        write(ilog,*) 'Warning. SHAKE did not converge within ',shake_maxiter,' steps. Inspect&
 & simulation results with utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(1),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(1).gt.0.5*HUGE(cs_wrnlmt(1))) then
          cs_wrncnt(1) = 0
        else
          cs_wrnlmt(1) = cs_wrnlmt(1)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_WARNINGS)
#endif
      exit
    end if
  end do
!
! now readjust cartesian velocities
! reconstructed velocities (override!) for the half-step t + 1/2dt
  do nc=1,constraints(grpid)%nats
    if (mass(constraints(grpid)%uidx(nc)).le.0.0) then
!     atoms in group have to be in order for this rebuilding of massless sites to work
      k = constraints(grpid)%uidx(nc)
      call genxyz(k,iz(1,k),blen(k),iz(2,k),bang(k),iz(3,k),ztor(k),iz(4,k))
    end if
    cart_v(constraints(grpid)%uidx(nc),1) = (x(constraints(grpid)%uidx(nc)) - xref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),2) = (y(constraints(grpid)%uidx(nc)) - yref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),3) = (z(constraints(grpid)%uidx(nc)) - zref(constraints(grpid)%uidx(nc)))*idt
  end do
!
! if we timed out, we infer that our vectors were insufficient (probably planar) -> update targets to avoid this
! problem re-occurring
! if, however, deviation large -> print extra warning
  if (iter.ge.shake_maxiter) then
    if (maxdev.gt.min(0.01d0,100.0*shake_tol)) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BL_WARNINGS)
#endif
      cs_wrncnt(2) = cs_wrncnt(2) + 1
      if (cs_wrncnt(2).eq.cs_wrnlmt(2)) then
        write(ilog,*) 'Warning. Maximum bond length deviation remains large after premature exit of SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(2),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(2).gt.0.5*HUGE(cs_wrnlmt(2))) then
          cs_wrncnt(2) = 0
        else
          cs_wrnlmt(2) = cs_wrnlmt(2)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BL_WARNINGS)
#endif
    else
      do nc=1,constraints(grpid)%nr
        if ((deviats(nc)/abs(constraints(grpid)%eqd2(nc))).gt.shake_tol) then
          constraints(grpid)%eqd2(nc) = constraints(grpid)%eqd2(nc) - deviats(nc)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BLA_WARNINGS)
#endif
      cs_wrncnt(4) = cs_wrncnt(4) + 1
      if (cs_wrncnt(4).eq.cs_wrnlmt(4)) then
        write(ilog,*) 'Warning. Resetting bond length targets (<1%) for SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(4),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(4).gt.0.5*HUGE(cs_wrnlmt(4))) then
          cs_wrncnt(4) = 0
        else
          cs_wrnlmt(4) = cs_wrnlmt(4)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BLA_WARNINGS)
#endif
    end if
    if (maxdev3.gt.min(0.01d0,100.0*shake_atol)) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BA_WARNINGS)
#endif
      cs_wrncnt(3) = cs_wrncnt(3) + 1
      if (cs_wrncnt(3).eq.cs_wrnlmt(3)) then
        write(ilog,*) 'Warning. Maximum bond angle deviation remains large after premature exit of SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(3).gt.0.5*HUGE(cs_wrnlmt(3))) then
          cs_wrncnt(3) = 0
        else
          cs_wrnlmt(3) = cs_wrnlmt(3)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BA_WARNINGS)
#endif
    else
      do nc=1,constraints(grpid)%nr3
        if (deviats3(nc).gt.shake_atol) then
          constraints(grpid)%eqa3(nc) = constraints(grpid)%eqa3(nc) - deviats3(nc)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BAA_WARNINGS)
#endif
      cs_wrncnt(5) = cs_wrncnt(5) + 1
      if (cs_wrncnt(5).eq.cs_wrnlmt(5)) then
        write(ilog,*) 'Warning. Resetting bond angle targets (<1%) for SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(5),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(5).gt.0.5*HUGE(cs_wrnlmt(5))) then
          cs_wrncnt(5) = 0
        else
          cs_wrnlmt(5) = cs_wrnlmt(5)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BAA_WARNINGS)
#endif
    end if
  end if
  constraints(grpid)%iters = constraints(grpid)%iters + iter
  constraints(grpid)%itercnt = constraints(grpid)%itercnt + 1
!
end
!
!---------------------------------------------------------------------------------------------------------
!
subroutine shake_rcb_just23(grpid)
!
  use atoms
  use shakeetal
  use system
  use molecule
  use iounit
  use mcsums
  use forces
  use math
  use zmatrix
!
  implicit none
! 
  integer grpid,nc,k,iter
  RTYPE dt2,maxdev,idt,maxdev3,dumm
  RTYPE dvec(3),d2,xvt,dveca1(3),dveca2(3),da1,da2,denomt,dotpt,vv1(3),vv2(3),vv3(3)
  RTYPE dvecref(constraints(grpid)%nr,3)
  RTYPE vvec(constraints(grpid)%nr,3)
  RTYPE vvec1(constraints(grpid)%nr,3)
  RTYPE vvec2(constraints(grpid)%nr,3)
  RTYPE vveca1(constraints(grpid)%nr3,3)
  RTYPE vveca2(constraints(grpid)%nr3,3)
  RTYPE vveca3(constraints(grpid)%nr3,3)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
  RTYPE mterm(constraints(grpid)%nr3,3)
  RTYPE dvecra1(constraints(grpid)%nr3,3)
  RTYPE dvecra2(constraints(grpid)%nr3,3)
  RTYPE dra1(constraints(grpid)%nr3)
  RTYPE dra2(constraints(grpid)%nr3)
  RTYPE dotp(constraints(grpid)%nr3)
  RTYPE denom(constraints(grpid)%nr3)
  RTYPE lambdas(constraints(grpid)%nr)
  RTYPE lambdas3(constraints(grpid)%nr3)
  RTYPE deviats(constraints(grpid)%nr)
  RTYPE deviats3(constraints(grpid)%nr3)
!
  dt2 = dyn_dt*dyn_dt
  idt = 1.0/dyn_dt
  do nc=1,constraints(grpid)%nr
    if (constraints(grpid)%massless(nc).EQV..true.) then
      mterm1(nc) = 0.0
      mterm2(nc) = 0.0
      cycle
    end if
    mterm1(nc) = dt2/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = dt2/mass(constraints(grpid)%idx(nc,2))
    call dis_bound3r(constraints(grpid)%idx(nc,2),constraints(grpid)%idx(nc,1),dvec(:),dumm)
    dvecref(nc,:) = dvec(:)
  end do
  do k=1,3
    vvec1(:,k) = dvecref(:,k)*mterm1(:)
    vvec2(:,k) = -dvecref(:,k)*mterm2(:)
    vvec(:,k) = vvec1(:,k) - vvec2(:,k)
  end do
!
  do nc=1,constraints(grpid)%nr3
    if (constraints(grpid)%massless(constraints(grpid)%nr+nc).EQV..true.) then
      mterm(nc,1) = 0.0
      mterm(nc,2) = 0.0
      mterm(nc,3) = 0.0
      cycle
    end if
    mterm(nc,1) = 0.5*dt2/mass(constraints(grpid)%idx3(nc,1))
    mterm(nc,2) = 0.5*dt2/mass(constraints(grpid)%idx3(nc,2))
    mterm(nc,3) = 0.5*dt2/mass(constraints(grpid)%idx3(nc,3))
    call dis_bound3r(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,1),dvec(:),dumm)
    dvecra1(nc,:) = dvec(:)
    call dis_bound3r(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,3),dvec(:),dumm)
    dvecra2(nc,:) = dvec(:)
  end do
  if (constraints(grpid)%nr3.gt.0) then
    dra1(:) = dvecra1(:,1)**2 + dvecra1(:,2)**2 + dvecra1(:,3)**2
    dra2(:) = dvecra2(:,1)**2 + dvecra2(:,2)**2 + dvecra2(:,3)**2
    denom(:) = 1.0/(dra1(:)*dra2(:))
    dra1(:) = sqrt(dra1(:))
    dra2(:) = sqrt(dra2(:))
    dotp(:) = dvecra1(:,1)*dvecra2(:,1) + dvecra1(:,2)*dvecra2(:,2) + dvecra1(:,3)*dvecra2(:,3)
    do k=1,3
      vveca1(:,k) = denom(:)*(dra1(:)*dra2(:)*dvecra2(:,k) - dotp(:)*dra2(:)*dvecra1(:,k)/dra1(:))
      vveca3(:,k) = denom(:)*(dra1(:)*dra2(:)*dvecra1(:,k) - dotp(:)*dra1(:)*dvecra2(:,k)/dra2(:))
      vveca2(:,k) = - (vveca1(:,k) + vveca3(:,k))
      vveca1(:,k) = vveca1(:,k)*mterm(:,1)
      vveca2(:,k) = vveca2(:,k)*mterm(:,2)
      vveca3(:,k) = vveca3(:,k)*mterm(:,3)
    end do
  end if
!
  iter = 0
  maxdev = 2.0*(shake_tol) + 1.0
  maxdev3 = 2.0*(shake_atol) + 1.0
!
  do while ((maxdev.gt.shake_tol).OR.(maxdev3.gt.shake_atol))
    iter = iter + 1 
    do nc=1,constraints(grpid)%nr
      if (constraints(grpid)%massless(nc).EQV..true.) cycle
      call dis_bound3(constraints(grpid)%idx(nc,2),constraints(grpid)%idx(nc,1),dvec(:),d2)
!
      xvt = dvec(1)*vvec(nc,1) + dvec(2)*vvec(nc,2) + dvec(3)*vvec(nc,3)
      if (abs(xvt).lt.(dt2*1.0e-8)) then
        deviats(nc) = 0.0
        lambdas(nc) = 0.0
      else
        xvt = 0.5/xvt
        deviats(nc) = constraints(grpid)%eqd2(nc) - d2
        lambdas(nc) = deviats(nc)*xvt
      end if
!
      x(constraints(grpid)%idx(nc,1)) = x(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,1)
      x(constraints(grpid)%idx(nc,2)) = x(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,1)
      y(constraints(grpid)%idx(nc,1)) = y(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,2)
      y(constraints(grpid)%idx(nc,2)) = y(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,2)
      z(constraints(grpid)%idx(nc,1)) = z(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,3)
      z(constraints(grpid)%idx(nc,2)) = z(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,3)
    end do
!
    do nc=1,constraints(grpid)%nr3
      if (constraints(grpid)%massless(constraints(grpid)%nr+nc).EQV..true.) cycle
      call dis_bound3(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,1),dveca1(:),da1)
      call dis_bound3(constraints(grpid)%idx3(nc,2),constraints(grpid)%idx3(nc,3),dveca2(:),da2)
      denomt = 1.0/(da1*da2)
      da1 = sqrt(da1)
      da2 = sqrt(da2)
      dotpt = dveca1(1)*dveca2(1) + dveca1(2)*dveca2(2) + dveca1(3)*dveca2(3)
      vv1(:) = denomt*(da1*da2*dveca2(:) - dotpt*da2*dveca1(:)/da1)
      vv3(:) = denomt*(da1*da2*dveca1(:) - dotpt*da1*dveca2(:)/da2)
      vv2(:) = - (vv1(:) + vv3(:))
      xvt = dot_product(vv1(:),vveca1(nc,:)) + dot_product(vv2(:),vveca2(nc,:)) + dot_product(vv3(:),vveca3(nc,:))
!      if (abs(xvt).lt.(dt2*1.0e-8)) then
!        deviats3(nc) = 0.0
!        lambdas3(nc) = 0.0
!      else
        xvt = 1.0/xvt
        deviats3(nc) = constraints(grpid)%eqa3(nc) - dotpt*sqrt(denomt)
        lambdas3(nc) = deviats3(nc)*xvt
!      end if
!
      x(constraints(grpid)%idx3(nc,1)) = x(constraints(grpid)%idx3(nc,1)) + lambdas3(nc)*vveca1(nc,1)
      x(constraints(grpid)%idx3(nc,2)) = x(constraints(grpid)%idx3(nc,2)) + lambdas3(nc)*vveca2(nc,1)
      x(constraints(grpid)%idx3(nc,3)) = x(constraints(grpid)%idx3(nc,3)) + lambdas3(nc)*vveca3(nc,1)
      y(constraints(grpid)%idx3(nc,1)) = y(constraints(grpid)%idx3(nc,1)) + lambdas3(nc)*vveca1(nc,2)
      y(constraints(grpid)%idx3(nc,2)) = y(constraints(grpid)%idx3(nc,2)) + lambdas3(nc)*vveca2(nc,2)
      y(constraints(grpid)%idx3(nc,3)) = y(constraints(grpid)%idx3(nc,3)) + lambdas3(nc)*vveca3(nc,2)
      z(constraints(grpid)%idx3(nc,1)) = z(constraints(grpid)%idx3(nc,1)) + lambdas3(nc)*vveca1(nc,3)
      z(constraints(grpid)%idx3(nc,2)) = z(constraints(grpid)%idx3(nc,2)) + lambdas3(nc)*vveca2(nc,3)
      z(constraints(grpid)%idx3(nc,3)) = z(constraints(grpid)%idx3(nc,3)) + lambdas3(nc)*vveca3(nc,3)
    end do

    maxdev = 0.0
    maxdev3 = 0.0
    do nc=1,constraints(grpid)%nr
      if (constraints(grpid)%massless(nc).EQV..true.) cycle
      if (abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc)).gt.maxdev) then
        maxdev = abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc))
      end if
    end do
    do nc=1,constraints(grpid)%nr3
      if (constraints(grpid)%massless(constraints(grpid)%nr+nc).EQV..true.) cycle
      if (abs(deviats3(nc)).gt.maxdev3) then
        maxdev3 = abs(deviats3(nc))
      end if
    end do
!
    if (iter.gt.shake_maxiter) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_WARNINGS)
#endif
      cs_wrncnt(1) = cs_wrncnt(1) + 1
      if (cs_wrncnt(1).eq.cs_wrnlmt(1)) then
        write(ilog,*) 'Warning. SHAKE did not converge within ',shake_maxiter,' steps. Inspect&
 & simulation results with utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(1),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(1).gt.0.5*HUGE(cs_wrnlmt(1))) then
          cs_wrncnt(1) = 0
        else
          cs_wrnlmt(1) = cs_wrnlmt(1)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_WARNINGS)
#endif
      exit
    end if
  end do
!
! now readjust cartesian velocities
! reconstructed velocities (override!) for the half-step t + 1/2dt
  do nc=1,constraints(grpid)%nats
    if (mass(constraints(grpid)%uidx(nc)).le.0.0) then
!     atoms in group have to be in order for this rebuilding of massless sites to work
      k = constraints(grpid)%uidx(nc)
      call genxyz(k,iz(1,k),blen(k),iz(2,k),bang(k),iz(3,k),ztor(k),iz(4,k))
    end if
    cart_v(constraints(grpid)%uidx(nc),1) = (x(constraints(grpid)%uidx(nc)) - xref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),2) = (y(constraints(grpid)%uidx(nc)) - yref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),3) = (z(constraints(grpid)%uidx(nc)) - zref(constraints(grpid)%uidx(nc)))*idt
  end do
!
! if we timed out, we infer that our vectors were insufficient (probably planar) -> update targets to avoid this
! problem re-occurring
! if, however, deviation large -> print extra warning
  if (iter.ge.shake_maxiter) then
    if (maxdev.gt.min(0.01d0,100.0*shake_tol)) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BL_WARNINGS)
#endif
      cs_wrncnt(2) = cs_wrncnt(2) + 1
      if (cs_wrncnt(2).eq.cs_wrnlmt(2)) then
        write(ilog,*) 'Warning. Maximum bond length deviation remains large after premature exit of SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(2),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(2).gt.0.5*HUGE(cs_wrnlmt(2))) then
          cs_wrncnt(2) = 0
        else
          cs_wrnlmt(2) = cs_wrnlmt(2)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BL_WARNINGS)
#endif
    else
      do nc=1,constraints(grpid)%nr
        if ((deviats(nc)/abs(constraints(grpid)%eqd2(nc))).gt.shake_tol) then
          constraints(grpid)%eqd2(nc) = constraints(grpid)%eqd2(nc) - deviats(nc)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BLA_WARNINGS)
#endif
      cs_wrncnt(4) = cs_wrncnt(4) + 1
      if (cs_wrncnt(4).eq.cs_wrnlmt(4)) then
        write(ilog,*) 'Warning. Resetting bond length targets (<1%) for SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(4),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(4).gt.0.5*HUGE(cs_wrnlmt(4))) then
          cs_wrncnt(4) = 0
        else
          cs_wrnlmt(4) = cs_wrnlmt(4)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BLA_WARNINGS)
#endif
    end if
    if (maxdev3.gt.min(0.01d0,100.0*shake_atol)) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BA_WARNINGS)
#endif
      cs_wrncnt(3) = cs_wrncnt(3) + 1
      if (cs_wrncnt(3).eq.cs_wrnlmt(3)) then
        write(ilog,*) 'Warning. Maximum bond angle deviation remains large after premature exit of SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(3).gt.0.5*HUGE(cs_wrnlmt(3))) then
          cs_wrncnt(3) = 0
        else
          cs_wrnlmt(3) = cs_wrnlmt(3)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BA_WARNINGS)
#endif
    else
      do nc=1,constraints(grpid)%nr3
        if (deviats3(nc).gt.shake_atol) then
          constraints(grpid)%eqa3(nc) = constraints(grpid)%eqa3(nc) - deviats3(nc)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BAA_WARNINGS)
#endif
      cs_wrncnt(5) = cs_wrncnt(5) + 1
      if (cs_wrncnt(5).eq.cs_wrnlmt(5)) then
        write(ilog,*) 'Warning. Resetting bond angle targets (<1%) for SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(5),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(5).gt.0.5*HUGE(cs_wrnlmt(5))) then
          cs_wrncnt(5) = 0
        else
          cs_wrnlmt(5) = cs_wrnlmt(5)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BAA_WARNINGS)
#endif
    end if
  end if
  constraints(grpid)%iters = constraints(grpid)%iters + iter
  constraints(grpid)%itercnt = constraints(grpid)%itercnt + 1
!
end
!
!---------------------------------------------------------------------------------------------------------
!
subroutine shake_rcb_just2(grpid)
!
  use atoms
  use shakeetal
  use system
  use molecule
  use iounit
  use mcsums
  use forces
  use math
  use zmatrix
!
  implicit none
! 
  integer grpid,nc,k,iter
  RTYPE dt2,maxdev,idt,dumm
  RTYPE dvec(3),d2,xvt
  RTYPE dvecref(constraints(grpid)%nr,3)
  RTYPE vvec(constraints(grpid)%nr,3)
  RTYPE vvec1(constraints(grpid)%nr,3)
  RTYPE vvec2(constraints(grpid)%nr,3)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
  RTYPE lambdas(constraints(grpid)%nr)
  RTYPE deviats(constraints(grpid)%nr)
!
  dt2 = dyn_dt*dyn_dt
  idt = 1.0/dyn_dt
  do nc=1,constraints(grpid)%nr
    if (constraints(grpid)%massless(nc).EQV..true.) then
      mterm1(nc) = 0.0
      mterm2(nc) = 0.0
      cycle
    end if
    mterm1(nc) = dt2/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = dt2/mass(constraints(grpid)%idx(nc,2))
    call dis_bound3r(constraints(grpid)%idx(nc,2),constraints(grpid)%idx(nc,1),dvec(:),dumm)
    dvecref(nc,:) = dvec(:)
  end do
  do k=1,3
    vvec1(:,k) = dvecref(:,k)*mterm1(:)
    vvec2(:,k) = -dvecref(:,k)*mterm2(:)
    vvec(:,k) = vvec1(:,k) - vvec2(:,k)
  end do
!
  iter = 0
  maxdev = 2.0*(shake_tol) + 1.0
!
  do while (maxdev.gt.shake_tol)
    iter = iter + 1 
    do nc=1,constraints(grpid)%nr
      if (constraints(grpid)%massless(nc).EQV..true.) cycle
      call dis_bound3(constraints(grpid)%idx(nc,2),constraints(grpid)%idx(nc,1),dvec(:),d2)
!
      xvt = dvec(1)*vvec(nc,1) + dvec(2)*vvec(nc,2) + dvec(3)*vvec(nc,3)
      if (abs(xvt).lt.(dt2*1.0e-8)) then
        deviats(nc) = 0.0
        lambdas(nc) = 0.0
      else
        xvt = 0.5/xvt
        deviats(nc) = constraints(grpid)%eqd2(nc) - d2
        lambdas(nc) = deviats(nc)*xvt
      end if
!
      x(constraints(grpid)%idx(nc,1)) = x(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,1)
      x(constraints(grpid)%idx(nc,2)) = x(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,1)
      y(constraints(grpid)%idx(nc,1)) = y(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,2)
      y(constraints(grpid)%idx(nc,2)) = y(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,2)
      z(constraints(grpid)%idx(nc,1)) = z(constraints(grpid)%idx(nc,1)) + lambdas(nc)*vvec1(nc,3)
      z(constraints(grpid)%idx(nc,2)) = z(constraints(grpid)%idx(nc,2)) + lambdas(nc)*vvec2(nc,3)
    end do
!
    maxdev = 0.0
    do nc=1,constraints(grpid)%nr
      if (constraints(grpid)%massless(nc).EQV..true.) cycle
      if (abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc)).gt.maxdev) then
        maxdev = abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc))
      end if
    end do
!
    if (iter.gt.shake_maxiter) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_WARNINGS)
#endif
      cs_wrncnt(1) = cs_wrncnt(1) + 1
      if (cs_wrncnt(1).eq.cs_wrnlmt(1)) then
        write(ilog,*) 'Warning. SHAKE did not converge within ',shake_maxiter,' steps. Inspect&
 & simulation results with utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(1),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(1).gt.0.5*HUGE(cs_wrnlmt(1))) then
          cs_wrncnt(1) = 0
        else
          cs_wrnlmt(1) = cs_wrnlmt(1)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_WARNINGS)
#endif
      exit
    end if
  end do
!
! now readjust cartesian velocities
! reconstructed velocities (override!) for the half-step t + 1/2dt
  do nc=1,constraints(grpid)%nats
    if (mass(constraints(grpid)%uidx(nc)).le.0.0) then
!     atoms in group have to be in order for this rebuilding of massless sites to work
      k = constraints(grpid)%uidx(nc)
      call genxyz(k,iz(1,k),blen(k),iz(2,k),bang(k),iz(3,k),ztor(k),iz(4,k))
    end if
    cart_v(constraints(grpid)%uidx(nc),1) = (x(constraints(grpid)%uidx(nc)) - xref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),2) = (y(constraints(grpid)%uidx(nc)) - yref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),3) = (z(constraints(grpid)%uidx(nc)) - zref(constraints(grpid)%uidx(nc)))*idt
  end do
!
! if we timed out, we infer that our vectors were insufficient (probably planar) -> update targets to avoid this
! problem re-occurring
! if, however, deviation large -> print extra warning
  if (iter.ge.shake_maxiter) then
    if (maxdev.gt.min(0.01d0,100.0*shake_tol)) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BL_WARNINGS)
#endif
      cs_wrncnt(2) = cs_wrncnt(2) + 1
      if (cs_wrncnt(2).eq.cs_wrnlmt(2)) then
        write(ilog,*) 'Warning. Maximum bond length deviation remains large after premature exit of SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(2),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(2).gt.0.5*HUGE(cs_wrnlmt(2))) then
          cs_wrncnt(2) = 0
        else
          cs_wrnlmt(2) = cs_wrnlmt(2)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BL_WARNINGS)
#endif
    else
      do nc=1,constraints(grpid)%nr
        if ((deviats(nc)/abs(constraints(grpid)%eqd2(nc))).gt.shake_tol) then
          constraints(grpid)%eqd2(nc) = constraints(grpid)%eqd2(nc) - deviats(nc)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SHAKE_BLA_WARNINGS)
#endif
      cs_wrncnt(4) = cs_wrncnt(4) + 1
      if (cs_wrncnt(4).eq.cs_wrnlmt(4)) then
        write(ilog,*) 'Warning. Resetting bond length targets (<1%) for SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(4),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(4).gt.0.5*HUGE(cs_wrnlmt(4))) then
          cs_wrncnt(4) = 0
        else
          cs_wrnlmt(4) = cs_wrnlmt(4)*10
        end if
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SHAKE_BLA_WARNINGS)
#endif
    end if
  end if
  constraints(grpid)%iters = constraints(grpid)%iters + iter
  constraints(grpid)%itercnt = constraints(grpid)%itercnt + 1
!
end
!
!---------------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
! attempt to parallelize solver for larger groups, only supports distance constraints
! must be used only for short constraints when using PBC
! constraints must be preorganized to allow operation on contiguous blocks with minimal
! overlap
!
subroutine shake_rcb_just2_threads(grpid,tpi,idx)
!
  use atoms
  use shakeetal
  use system
  use iounit
  use mcsums
  use forces
  use zmatrix
  use threads
  use cutoffs
!
  implicit none
! 
  integer, INTENT(IN):: grpid,tpi,idx
!
!  integer(KIND=8) ttt(2)
  integer nc,nc2,k,iter,ki,kj,alsz,alsz2,j
  RTYPE dt2,maxdev,idt,minidt2,lambda
  RTYPE dvec(3),d2,xvt,lopb(3),lopc(3)
  RTYPE vvec(3,max(1,ccg_limits(idx,2,tpi)-ccg_limits(idx,1,tpi)+1))
  RTYPE vvec1(3,max(1,ccg_limits(idx,2,tpi)-ccg_limits(idx,1,tpi)+1))
  RTYPE vvec2(3,max(1,ccg_limits(idx,2,tpi)-ccg_limits(idx,1,tpi)+1))
  integer ivec(2,max(1,ccg_limits(idx,2,tpi)-ccg_limits(idx,1,tpi)+1))
  RTYPE deviats(max(1,ccg_limits(idx,2,tpi)-ccg_limits(idx,1,tpi)+1))
  RTYPE id2(max(1,ccg_limits(idx,2,tpi)-ccg_limits(idx,1,tpi)+1))
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).OR.((tpi.gt.1).AND.(OMP_IN_PARALLEL().EQV..false.))) then
    write(ilog,*) 'Fatal. Subroutine shake_rcb_just2_threads(...) can only be called &
 &by all threads of a team in a parallel binding region. This is a bug.'
    call fexit()
  end if
!  if (tpi.eq.1) then
!    call System_Clock(ttt(1))
!  end if
!
  ki = ccg_limits(idx,1,tpi)
  kj = ccg_limits(idx,2,tpi)
  alsz = kj-ki+1
  dt2 = dyn_dt*dyn_dt
  minidt2 = 1.0e-8*dt2
  idt = 1.0/dyn_dt
! deal with potential image shifts for intermolecular constraints
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
      do j=1,3
        lopb(j) = 1.0/bnd_params(j)
        lopc(j) = -bnd_params(j)
      end do
    else if (bnd_shape.eq.3) then
      lopb(3) = 1.0/bnd_params(6)
      lopc(3) = -bnd_params(6)
    end if
  end if
  do nc=ki,kj
    nc2 = nc-ki+1
    if (constraints(grpid)%massless(nc).EQV..true.) then
      cycle
    end if
    dvec(1) = xref(constraints(grpid)%idx(nc,1))-xref(constraints(grpid)%idx(nc,2))
    dvec(2) = yref(constraints(grpid)%idx(nc,1))-yref(constraints(grpid)%idx(nc,2))
    dvec(3) = zref(constraints(grpid)%idx(nc,1))-zref(constraints(grpid)%idx(nc,2))
    if (constraints(grpid)%isinter(nc).EQV..true.) then
      if (bnd_type.eq.1) then
!       cubic box
        if (bnd_shape.eq.1) then
          dvec(1) = dvec(1) + lopc(1)*anint((xref(constraints(grpid)%idx(nc,1))-xref(constraints(grpid)%idx(nc,2)))*lopb(1))
          dvec(2) = dvec(2) + lopc(2)*anint((yref(constraints(grpid)%idx(nc,1))-yref(constraints(grpid)%idx(nc,2)))*lopb(2))
          dvec(3) = dvec(3) + lopc(3)*anint((zref(constraints(grpid)%idx(nc,1))-zref(constraints(grpid)%idx(nc,2)))*lopb(3))
        else if (bnd_shape.eq.3) then
          dvec(3) = dvec(3) + lopc(3)*anint((zref(constraints(grpid)%idx(nc,1))-zref(constraints(grpid)%idx(nc,2)))*lopb(3))
        end if
      end if
    end if
    vvec1(:,nc2) = dvec(:)*dt2/mass(constraints(grpid)%idx(nc,1))
    vvec2(:,nc2) = -dvec(:)*dt2/mass(constraints(grpid)%idx(nc,2))
    ivec(:,nc2) = constraints(grpid)%idx(nc,1:2)
    id2(nc2) = 1.0/constraints(grpid)%eqd2(nc)
  end do
  if (alsz.gt.0) vvec(:,1:alsz) = vvec1(:,1:alsz) - vvec2(:,1:alsz)
!
  iter = 0
  maxdev = 2.0*(shake_tol) + 1.0
!$OMP BARRIER
!
  do while ((maxval(ccg_limits(idx,2,1:thrdat%maxn)).gt.0).AND.(iter.le.shake_maxiter))
    iter = iter + 1 
    if (ccg_limits(idx,2,tpi).eq.constraints(grpid)%nr) then
      alsz2 = alsz
    else
      alsz2 = alsz-1
    end if
    do nc=1,alsz2 ! the last constraint in the set is the one with overlap with another set (except for the last nonzero set)

      nc2 = nc + ccg_limits(idx,1,tpi) - 1
      if (constraints(grpid)%massless(nc2).EQV..true.) cycle
      dvec(:) = atinfo(1:3,ivec(1,nc))-atinfo(1:3,ivec(2,nc))
      if (constraints(grpid)%isinter(nc).EQV..true.) then
        if (bnd_shape.eq.1) then
          do j=1,3
            dvec(j) = dvec(j) + lopc(j)*anint((atinfo(j,ivec(1,nc))-atinfo(j,ivec(2,nc)))*lopb(j))
          end do
        else if (bnd_shape.eq.3) then
          dvec(3) = dvec(3) + lopc(3)*anint((atinfo(3,ivec(1,nc))-atinfo(3,ivec(2,nc)))*lopb(3))
        end if
      end if
      d2 = dot_product(dvec,dvec)
      xvt = dot_product(dvec,vvec(:,nc))
!
      if (abs(xvt).lt.(minidt2)) then
        deviats(nc) = 0.0
        lambda = 0.0
        cycle
      else
        xvt = 0.5/xvt
        deviats(nc) = constraints(grpid)%eqd2(nc2) - d2
        lambda = deviats(nc)*xvt
        atinfo(1:3,ivec(1,nc)) = atinfo(1:3,ivec(1,nc)) + lambda*vvec1(:,nc)
        atinfo(1:3,ivec(2,nc)) = atinfo(1:3,ivec(2,nc)) + lambda*vvec2(:,nc)
      end if
    end do
!
!$OMP BARRIER
    if ((alsz.gt.0).AND.(alsz2.lt.alsz)) then
      nc2 = ccg_limits(idx,1,tpi) + alsz - 1 ! ccg_limits(idx,2,tpi)
      dvec(:) = atinfo(1:3,ivec(1,alsz))-atinfo(1:3,ivec(2,alsz))
      d2 = dot_product(dvec,dvec)
      xvt = dot_product(dvec,vvec(:,alsz))
!
      if (abs(xvt).lt.(minidt2)) then
        deviats(alsz) = 0.0
        lambda = 0.0
      else
        xvt = 0.5/xvt
        deviats(alsz) = constraints(grpid)%eqd2(nc2) - d2
        lambda = deviats(alsz)*xvt
        atinfo(1:3,ivec(1,alsz)) = atinfo(1:3,ivec(1,alsz)) + lambda*vvec1(:,alsz)
        atinfo(1:3,ivec(2,alsz)) = atinfo(1:3,ivec(2,alsz)) + lambda*vvec2(:,alsz)
      end if
    end if
!
    maxdev = 0.0
    do nc=1,alsz
      nc2 = nc + ccg_limits(idx,1,tpi) - 1
      if (constraints(grpid)%massless(nc2).EQV..true.) cycle
      lambda = abs(deviats(nc))*id2(nc)
      if (lambda.gt.maxdev) maxdev = lambda
    end do
    if (maxdev.le.shake_tol) ccg_limits(idx,2,tpi) = 0
!
    if (iter.gt.shake_maxiter) then
!$OMP SINGLE
!$OMP CRITICAL(SHAKE_WARNINGS)
      cs_wrncnt(1) = cs_wrncnt(1) + 1
      if (cs_wrncnt(1).eq.cs_wrnlmt(1)) then
        write(ilog,*) 'Warning. SHAKE did not converge within ',shake_maxiter,' steps. Inspect&
 & simulation results with utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(1),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(1).gt.0.5*HUGE(cs_wrnlmt(1))) then
          cs_wrncnt(1) = 0
        else
          cs_wrnlmt(1) = cs_wrnlmt(1)*10
        end if
      end if
!$OMP END CRITICAL(SHAKE_WARNINGS)
!$OMP END SINGLE
    end if
!$OMP BARRIER
  end do
!
!
! transfer coordinates and readjust cartesian velocities
  do nc=ccg_limits(idx,3,tpi),ccg_limits(idx,4,tpi)
    k = constraints(grpid)%uidx(nc)
    x(k) = atinfo(1,k)
    y(k) = atinfo(2,k)
    z(k) = atinfo(3,k)
  end do
!$OMP BARRIER
! restore
  ccg_limits(idx,2,tpi) = ccg_limits(idx,1,tpi) + alsz - 1
  do nc=ccg_limits(idx,3,tpi),ccg_limits(idx,4,tpi)
    k = constraints(grpid)%uidx(nc)
    if (mass(constraints(grpid)%uidx(nc)).le.0.0) then
!     atoms in group have to be in order for this rebuilding of massless sites to work
      call genxyz(k,iz(1,k),blen(k),iz(2,k),bang(k),iz(3,k),ztor(k),iz(4,k))
    end if
    cart_v(k,1) = (x(k) - xref(k))*idt
    cart_v(k,2) = (y(k) - yref(k))*idt
    cart_v(k,3) = (z(k) - zref(k))*idt
  end do
!
! if we timed out, we infer that our vectors were insufficient (probably planar) -> update targets to avoid this
! problem re-occurring
! if, however, deviation large -> print extra warning
  if (iter.ge.shake_maxiter) then
    if (maxdev.gt.min(0.01d0,100.0*shake_tol)) then
!$OMP CRITICAL(SHAKE_BL_WARNINGS)
      cs_wrncnt(2) = cs_wrncnt(2) + 1
      if (cs_wrncnt(2).eq.cs_wrnlmt(2)) then
        write(ilog,*) 'Warning. Maximum bond length deviation remains large after premature exit of SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(2),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(2).gt.0.5*HUGE(cs_wrnlmt(2))) then
          cs_wrncnt(2) = 0
        else
          cs_wrnlmt(2) = cs_wrnlmt(2)*10
        end if
      end if
!$OMP END CRITICAL(SHAKE_BL_WARNINGS)
    else if (maxdev.gt.shake_tol) then
      do nc=ki,kj
        nc2 = nc-ki+1
        if ((abs(deviats(nc2))*id2(nc2)).gt.shake_tol) then
          constraints(grpid)%eqd2(nc) = constraints(grpid)%eqd2(nc) - deviats(nc2)
        end if
      end do
!$OMP CRITICAL(SHAKE_BLA_WARNINGS)
      cs_wrncnt(4) = cs_wrncnt(4) + 1
      if (cs_wrncnt(4).eq.cs_wrnlmt(4)) then
        write(ilog,*) 'Warning. Resetting bond length targets (<1%) for SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(4),' of this type not all of which may be displayed.'
        if (10.0*cs_wrnlmt(4).gt.0.5*HUGE(cs_wrnlmt(4))) then
          cs_wrncnt(4) = 0
        else
          cs_wrnlmt(4) = cs_wrnlmt(4)*10
        end if
      end if
!$OMP END CRITICAL(SHAKE_BLA_WARNINGS)
    end if
  end if
!$OMP SINGLE
  constraints(grpid)%iters = constraints(grpid)%iters + iter
  constraints(grpid)%itercnt = constraints(grpid)%itercnt + 1
!$OMP END SINGLE NOWAIT
!
!  if (tpi.eq.1) then
!    call System_Clock(ttt(2))
!    write(*,*) 1.0e-3*(ttt(2)-ttt(1))
!  end if
!
end
!
#endif
!
!---------------------------------------------------------------------------------------------------------
!

subroutine settle(grpid,mvals,rvals)
!
! mvals: 1-2 fractional masses
! rvals: 1: mass-weighted projection length of OH onto HH; 2: remainder of 1
! rvals: 3: half distance between terminal atoms, 4: distance between terminal atoms
!
  use atoms
  use shakeetal
  use system
  use iounit
  use mcsums
  use forces
  use zmatrix
!
  implicit none
!
  integer nc,k,grpid
  RTYPE mvals(2),rvals(4)
  RTYPE comn(3),rvec(3,3),svec(3,3),ovec(3,3),ivec(3,2),tvec(3,3),idt,ang1,ang2,ang3
  RTYPE tvn2(3),itvn(3),sine1,cose1,sine2,cose2,sine3,cose3,val1,val2,hlp(6),dp(11)
!
  idt = 1.0/dyn_dt
!
  ivec(1,1) = xref(constraints(grpid)%uidx(2)) - xref(constraints(grpid)%uidx(1))
  ivec(2,1) = yref(constraints(grpid)%uidx(2)) - yref(constraints(grpid)%uidx(1))
  ivec(3,1) = zref(constraints(grpid)%uidx(2)) - zref(constraints(grpid)%uidx(1))
  ivec(1,2) = xref(constraints(grpid)%uidx(3)) - xref(constraints(grpid)%uidx(1))
  ivec(2,2) = yref(constraints(grpid)%uidx(3)) - yref(constraints(grpid)%uidx(1))
  ivec(3,2) = zref(constraints(grpid)%uidx(3)) - zref(constraints(grpid)%uidx(1))
!
  comn(1) = mvals(1)*x(constraints(grpid)%uidx(1)) + mvals(2)*(x(constraints(grpid)%uidx(2)) + &
 &          x(constraints(grpid)%uidx(3)))
  comn(2) = mvals(1)*y(constraints(grpid)%uidx(1)) + mvals(2)*(y(constraints(grpid)%uidx(2)) + &
 &          y(constraints(grpid)%uidx(3)))
  comn(3) = mvals(1)*z(constraints(grpid)%uidx(1)) + mvals(2)*(z(constraints(grpid)%uidx(2)) + &
 &          z(constraints(grpid)%uidx(3)))
  rvec(1,1) = x(constraints(grpid)%uidx(1)) - comn(1)
  rvec(2,1) = y(constraints(grpid)%uidx(1)) - comn(2)
  rvec(3,1) = z(constraints(grpid)%uidx(1)) - comn(3)
  rvec(1,2) = x(constraints(grpid)%uidx(2)) - comn(1)
  rvec(2,2) = y(constraints(grpid)%uidx(2)) - comn(2)
  rvec(3,2) = z(constraints(grpid)%uidx(2)) - comn(3)
  rvec(1,3) = x(constraints(grpid)%uidx(3)) - comn(1)
  rvec(2,3) = y(constraints(grpid)%uidx(3)) - comn(2)
  rvec(3,3) = z(constraints(grpid)%uidx(3)) - comn(3)
  call crossprod2(ivec(:,1),ivec(:,2),tvec(:,3),tvn2(3))
  call crossprod2(rvec(:,1),tvec(:,3),tvec(:,1),tvn2(1))
  call crossprod2(tvec(:,3),tvec(:,1),tvec(:,2),tvn2(2))
  itvn(:) = sqrt(1.0/tvn2(:))
  do k=1,3
    tvec(:,k) = tvec(:,k)*itvn(k)
  end do
!
  dp(5) = dot_product(tvec(:,3),rvec(:,1))
  sine1 = dp(5)/rvals(1)
  val1 = (1.0 - sine1*sine1)
  if (val1.le.0.0) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SETTLE_WARNINGS)
#endif
    cs_wrncnt(9) = cs_wrncnt(9) + 1
    if (cs_wrncnt(9).eq.cs_wrnlmt(9)) then
      write(ilog,*) 'Warning. SETTLE failed due to incompatible geometry, and SHAKE will be attempted &
 &instead.  This simulation must be analyzed with the utmost care.'
      write(ilog,*) 'This was warning #',cs_wrncnt(9),' of this type not all of which may be displayed.'
      if (10.0*cs_wrnlmt(9).gt.0.5*HUGE(cs_wrnlmt(9))) then
        cs_wrncnt(9) = 0
      else
        cs_wrnlmt(9) = cs_wrnlmt(9)*10
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SETTLE_WARNINGS)
#endif
    call shake_rcb_just2(grpid)
    return
  end if
  cose1 = val1*sqrt(1.0/val1)
!
  dp(8) = dot_product(tvec(:,3),rvec(:,2))
  dp(11) = dot_product(tvec(:,3),rvec(:,3))
  sine2 = (dp(8) - dp(11)) / (rvals(4)*cose1)
  val2 = 1.0 - sine2*sine2
  if (val2.le.0.0) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(SETTLE_WARNINGS)
#endif
    cs_wrncnt(9) = cs_wrncnt(9) + 1
    if (cs_wrncnt(9).eq.cs_wrnlmt(9)) then
      write(ilog,*) 'Warning. SETTLE failed due to incompatible geometry, and SHAKE will be attempted &
 &instead.  This simulation must be analyzed with the utmost care.'
      write(ilog,*) 'This was warning #',cs_wrncnt(9),' of this type not all of which may be displayed.'
      if (10.0*cs_wrnlmt(9).gt.0.5*HUGE(cs_wrnlmt(9))) then
        cs_wrncnt(9) = 0
      else
        cs_wrnlmt(9) = cs_wrnlmt(9)*10
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(SETTLE_WARNINGS)
#endif
    call shake_rcb_just2(grpid)
    return
  end if
  cose2 = val2*(sqrt(1.0/val2))
!
  dp(1) = dot_product(tvec(:,1),ivec(:,1))
  dp(2) = dot_product(tvec(:,2),ivec(:,1))
  dp(3) = dot_product(tvec(:,1),ivec(:,2))
  dp(4) = dot_product(tvec(:,2),ivec(:,2))
  dp(6) = dot_product(tvec(:,1),rvec(:,2))
  dp(7) = dot_product(tvec(:,2),rvec(:,2))
  dp(9) = dot_product(tvec(:,1),rvec(:,3))
  dp(10) = dot_product(tvec(:,2),rvec(:,3))
  hlp(1) = rvals(1)*cose1
  hlp(2) = -rvals(3)*cose2
  hlp(3) = -rvals(2)*cose1
  hlp(4) = rvals(3)*sine1*sine2
  hlp(5) = hlp(3) - hlp(4)
  hlp(6) = hlp(3) + hlp(4)
!
  ang1 = hlp(2) * (dp(1) - dp(3)) + dp(2)*hlp(5) + dp(4)*hlp(6)
  ang2 = hlp(2) * (dp(4) - dp(2)) + dp(1)*hlp(5) + dp(3)*hlp(6)
  ang3 = dp(1)*dp(7) - dp(6)*dp(2) + dp(3)*dp(10) - dp(9)*dp(4)
  val2 = ang1*ang1 + ang2*ang2
  val1 = val2 - ang3*ang3
  sine3 = (ang1*ang3 - ang2*val1*sqrt(1.0/val1)) / val2
  val1 = 1.0 - sine3*sine3
  cose3 = val1*sqrt(1.0/val1)
!
  svec(1,1) = -sine3*hlp(1)
  svec(2,1) = cose3*hlp(1)
  svec(3,1) = dp(5)
  svec(1,2) = cose3*hlp(2) - sine3*hlp(5)
  svec(2,2) = sine3*hlp(2) + cose3*hlp(5)
  svec(3,2) = dp(8)
  svec(1,3) = -cose3*hlp(2) - sine3*hlp(6)
  svec(2,3) = -sine3*hlp(2) + cose3*hlp(6)
  svec(3,3) = dp(11)
  ovec(:,1) = tvec(:,1)*svec(1,1) + tvec(:,2)*svec(2,1) + tvec(:,3)*svec(3,1)
  ovec(:,2) = tvec(:,1)*svec(1,2) + tvec(:,2)*svec(2,2) + tvec(:,3)*svec(3,2)
  ovec(:,3) = tvec(:,1)*svec(1,3) + tvec(:,2)*svec(2,3) + tvec(:,3)*svec(3,3)
!
! lastly, update coordinates
  x(constraints(grpid)%uidx(1)) = comn(1) + ovec(1,1)
  y(constraints(grpid)%uidx(1)) = comn(2) + ovec(2,1)
  z(constraints(grpid)%uidx(1)) = comn(3) + ovec(3,1)
  x(constraints(grpid)%uidx(2)) = comn(1) + ovec(1,2)
  y(constraints(grpid)%uidx(2)) = comn(2) + ovec(2,2)
  z(constraints(grpid)%uidx(2)) = comn(3) + ovec(3,2)
  x(constraints(grpid)%uidx(3)) = comn(1) + ovec(1,3)
  y(constraints(grpid)%uidx(3)) = comn(2) + ovec(2,3)
  z(constraints(grpid)%uidx(3)) = comn(3) + ovec(3,3)
!
! for water molecules with massless sites, the assumed principle of least constraints
! solution can be obtained by simply rebuilding the massless points based on rigid geometry
!  if (at(atmres(constraints(grpid)%uidx(1)))%
!
! and fix velocities
  do nc=1,constraints(grpid)%nats
    if (mass(constraints(grpid)%uidx(nc)).le.0.0) then
!     atoms in group have to be in order for this rebuilding of massless sites to work
      k = constraints(grpid)%uidx(nc)
      call genxyz(k,iz(1,k),blen(k),iz(2,k),bang(k),iz(3,k),ztor(k),iz(4,k))
    end if
    cart_v(constraints(grpid)%uidx(nc),1) = (x(constraints(grpid)%uidx(nc)) - xref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),2) = (y(constraints(grpid)%uidx(nc)) - yref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),3) = (z(constraints(grpid)%uidx(nc)) - zref(constraints(grpid)%uidx(nc)))*idt
  end do
!
end
!
!---------------------------------------------------------------------------------------------------------
!
! this is pre-conditioned SHAKE (J. Comput. Phys., 220 (2007) 740-750)
! note that this method only works infrequently for systems of coupling
! that is too high for SHAKE  - possible implementation bug that causes just linearization to fail?
!
subroutine shake_p(grpid)
!
  use atoms
  use shakeetal
  use system
  use iounit
  use mcsums
  use forces
!
  implicit none
! 
  integer nc,iter,i,j,k,np,grpid
  RTYPE dt2,maxdev,idt
  RTYPE dvec(constraints(grpid)%nr,3)
  RTYPE dvecref(constraints(grpid)%nr,3)
  RTYPE vmatp(3*constraints(grpid)%nats,constraints(grpid)%nr)
  RTYPE vmat(3*constraints(grpid)%nats,constraints(grpid)%nr)
  RTYPE vvec(constraints(grpid)%nr,3)
  RTYPE d2(constraints(grpid)%nr)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
  RTYPE lambdas(constraints(grpid)%nr)
  RTYPE deviats(constraints(grpid)%nr)
  RTYPE xvt(constraints(grpid)%nr)
!
  dt2 = dyn_dt*dyn_dt
  idt = 1.0/dyn_dt
!
  np = 0
  vmat(:,:) = 0.0
  vmatp(:,:) = 0.0
  do nc=1,constraints(grpid)%nr
    mterm1(nc) = dt2/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = dt2/mass(constraints(grpid)%idx(nc,2))
    dvecref(nc,1) = xref(constraints(grpid)%idx(nc,1)) - xref(constraints(grpid)%idx(nc,2))
    dvecref(nc,2) = yref(constraints(grpid)%idx(nc,1)) - yref(constraints(grpid)%idx(nc,2))
    dvecref(nc,3) = zref(constraints(grpid)%idx(nc,1)) - zref(constraints(grpid)%idx(nc,2))
    vmatp((3*constraints(grpid)%mapidx(nc,1)-2):(3*constraints(grpid)%mapidx(nc,1)),nc) = dvecref(nc,:)*mterm1(nc)
    vmatp((3*constraints(grpid)%mapidx(nc,2)-2):(3*constraints(grpid)%mapidx(nc,2)),nc) = -dvecref(nc,:)*mterm2(nc)
  end do
!
! sparse matrix multiplication
  do nc=1,constraints(grpid)%nr
    i = 3*constraints(grpid)%mapidx(nc,1)-2
    j = 3*constraints(grpid)%mapidx(nc,1)
    do k=1,constraints(grpid)%nr
      vmat(i:j,k) = vmat(i:j,k) + vmatp(i:j,nc)*constraints(grpid)%amat(nc,k)
    end do
    i = 3*constraints(grpid)%mapidx(nc,2)-2
    j = 3*constraints(grpid)%mapidx(nc,2)
    do k=1,constraints(grpid)%nr
      vmat(i:j,k) = vmat(i:j,k) + vmatp(i:j,nc)*constraints(grpid)%amat(nc,k)
    end do
  end do
!  10 format(8(g12.5,1x))
!  write(*,10) transpose(constraints(grpid)%amat)
!  write(*,10) transpose(vmat)
  do nc=1,constraints(grpid)%nr
    vvec(nc,:) = vmat((3*constraints(grpid)%mapidx(nc,1)-2):(3*constraints(grpid)%mapidx(nc,1)),nc) - &
 &               vmat((3*constraints(grpid)%mapidx(nc,2)-2):(3*constraints(grpid)%mapidx(nc,2)),nc)
  end do
!
  iter = 0
  maxdev = (shake_tol+1.0)*2.0
  do while (maxdev.gt.shake_tol)
    iter = iter + 1 
    do nc=1,constraints(grpid)%nr
      dvec(nc,1) = x(constraints(grpid)%idx(nc,1)) - x(constraints(grpid)%idx(nc,2))
      dvec(nc,2) = y(constraints(grpid)%idx(nc,1)) - y(constraints(grpid)%idx(nc,2))
      dvec(nc,3) = z(constraints(grpid)%idx(nc,1)) - z(constraints(grpid)%idx(nc,2))
    end do
!
    d2(:) = dvec(:,1)*dvec(:,1) + dvec(:,2)*dvec(:,2) + dvec(:,3)*dvec(:,3)
!
    xvt(:) = 0.5/(dvec(:,1)*vvec(:,1) + dvec(:,2)*vvec(:,2) + dvec(:,3)*vvec(:,3))
    deviats(:) = constraints(grpid)%eqd2(:) - d2(:)
!
    maxdev = 0.0
    do nc=1,constraints(grpid)%nr
      if (abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc)).gt.maxdev) then
        maxdev = abs(deviats(nc))/abs(constraints(grpid)%eqd2(nc))
      end if
    end do
!
    lambdas(:) = deviats(:)*xvt(:)
    do i=1,constraints(grpid)%nats
      x(constraints(grpid)%uidx(i)) = x(constraints(grpid)%uidx(i)) + dot_product(lambdas(:),vmat(3*i-2,:))
      y(constraints(grpid)%uidx(i)) = y(constraints(grpid)%uidx(i)) + dot_product(lambdas(:),vmat(3*i-1,:))
      z(constraints(grpid)%uidx(i)) = z(constraints(grpid)%uidx(i)) + dot_product(lambdas(:),vmat(3*i,:))
    end do
!
    if (iter.gt.shake_maxiter) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(PSHAKE_WARNINGS)
#endif
      cs_wrncnt(11) = cs_wrncnt(11) + 1
      if (cs_wrncnt(11).eq.cs_wrnlmt(11)) then
        write(ilog,*) 'Warning. P-SHAKE did not converge within ',shake_maxiter,' steps. Inspect&
 & simulation results with utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(11),' of this type not all of which may be displayed.'
        cs_wrnlmt(11) = cs_wrnlmt(11)*10
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(PSHAKE_WARNINGS)
#endif
      exit
    end if
  end do
!
! now readjust cartesian velocities
! reconstructed velocities (override!) for the half-step t + 1/2dt
  do nc=1,constraints(grpid)%nats
    cart_v(constraints(grpid)%uidx(nc),1) = (x(constraints(grpid)%uidx(nc)) - xref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),2) = (y(constraints(grpid)%uidx(nc)) - yref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),3) = (z(constraints(grpid)%uidx(nc)) - zref(constraints(grpid)%uidx(nc)))*idt
  end do
!
! if we timed out, we infer that our vectors were insufficient (probably planar) -> update targets to avoid this
! problem re-occurring
! if, however, deviation large -> print extra warning
  if (iter.ge.shake_maxiter) then
    if (maxdev.gt.min(0.01d0,100.0*shake_tol)) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(PSHAKE_BL_WARNINGS)
#endif
      cs_wrncnt(12) = cs_wrncnt(12) + 1
      if (cs_wrncnt(12).eq.cs_wrnlmt(12)) then
        write(ilog,*) 'Warning. Maximum bond length deviation remains large after premature exit of P-SHAKE. This&
 & simulation is unlikely to be stable.'
        write(ilog,*) 'This was warning #',cs_wrncnt(12),' of this type not all of which may be displayed.'
        cs_wrnlmt(12) = cs_wrnlmt(12)*10
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(PSHAKE_BL_WARNINGS)
#endif
    else
      do nc=1,constraints(grpid)%nr
        if ((deviats(nc)/abs(constraints(grpid)%eqd2(nc))).gt.shake_tol) then
          constraints(grpid)%eqd2(nc) = constraints(grpid)%eqd2(nc) - deviats(nc)
        end if
      end do
#ifdef ENABLE_THREADS
!$OMP CRITICAL(PSHAKE_BLA_WARNINGS)
#endif
      cs_wrncnt(14) = cs_wrncnt(14) + 1
      if (cs_wrncnt(14).eq.cs_wrnlmt(14)) then
        write(ilog,*) 'Warning. Resetting bond length targets (<1%) for P-SHAKE to circumvent poor convergence. This&
 & simulation must be analyzed with the utmost care.'
        write(ilog,*) 'This was warning #',cs_wrncnt(14),' of this type not all of which may be displayed.'
        cs_wrnlmt(14) = cs_wrnlmt(14)*10
      end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(PSHAKE_BLA_WARNINGS)
#endif
    end if
  end if
!
  constraints(grpid)%iters = constraints(grpid)%iters + iter
  constraints(grpid)%itercnt = constraints(grpid)%itercnt + 1
!
end
!
!---------------------------------------------------------------------------------------------------------
!
#ifdef LINK_LAPACK
!
! this is the setup routine for the pre-conditioner in P-SHAKE (J. Comput. Phys., 220, 740-750 (2007))
!
subroutine shake_p_setA(grpid)
!
  use atoms
  use shakeetal
  use system
  use molecule
  use iounit
!
  implicit none
!
  integer nc,i1,j1,i2,j2,k,kk,iter,grpid,lapinf
  integer rowidx(constraints(grpid)%nr-1)
  RTYPE maxdev,dt2
  RTYPE dvec(constraints(grpid)%nr,3)
  RTYPE rmat(constraints(grpid)%nr-1,constraints(grpid)%nr-1)
  RTYPE rmatref(constraints(grpid)%nr,constraints(grpid)%nr)
  RTYPE hvec2(constraints(grpid)%nr-1)
  RTYPE rmat2(constraints(grpid)%nr-1,constraints(grpid)%nr-1)
  RTYPE amattmp(constraints(grpid)%nr,constraints(grpid)%nr)
  RTYPE vmat(3*constraints(grpid)%nats,constraints(grpid)%nr)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
!
  dt2 = dyn_dt*dyn_dt
!
  if (constraints(grpid)%nr.eq.1) then
    constraints(grpid)%amat(1,1) = 1.0
    return
  end if
!
  constraints(grpid)%amat(:,:) = 0.0
  do k=1,constraints(grpid)%nr
    constraints(grpid)%amat(k,k) = 1.0
  end do
!
  vmat(:,:) = 0.0
  do nc=1,constraints(grpid)%nr
    mterm1(nc) = dt2/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = dt2/mass(constraints(grpid)%idx(nc,2))
    dvec(nc,1) = x(constraints(grpid)%idx(nc,1)) - x(constraints(grpid)%idx(nc,2))
    dvec(nc,2) = y(constraints(grpid)%idx(nc,1)) - y(constraints(grpid)%idx(nc,2))
    dvec(nc,3) = z(constraints(grpid)%idx(nc,1)) - z(constraints(grpid)%idx(nc,2))
    vmat((3*constraints(grpid)%mapidx(nc,1)-2):(3*constraints(grpid)%mapidx(nc,1)),nc) = dvec(nc,:)*mterm1(nc)
    vmat((3*constraints(grpid)%mapidx(nc,2)-2):(3*constraints(grpid)%mapidx(nc,2)),nc) = -dvec(nc,:)*mterm2(nc)
  end do
!
  maxdev = 1.0
  iter = 0
  do while (maxdev.gt.1.0e-8)
!
    iter = iter + 1 
!
    amattmp(:,:) = 0.0
    do nc=1,constraints(grpid)%nr
      amattmp(nc,nc) = 1.0
    end do
 10  format(18(5(g13.6),/))
    do nc=1,constraints(grpid)%nr
      rmatref(:,nc) = (vmat(3*(constraints(grpid)%mapidx(nc,1))-2,:)-vmat(3*(constraints(grpid)%mapidx(nc,2))-2,:))*dvec(nc,1) + &
 &                    (vmat(3*(constraints(grpid)%mapidx(nc,1))-1,:)-vmat(3*(constraints(grpid)%mapidx(nc,2))-1,:))*dvec(nc,2) + &
 &                    (vmat(3*(constraints(grpid)%mapidx(nc,1)),:)-vmat(3*(constraints(grpid)%mapidx(nc,2)),:))*dvec(nc,3)
    end do
!
    do nc=1,constraints(grpid)%nr
!     assemble r
      i1 = 1
      j1 = nc-1
      i2 = nc+1
      j2 = constraints(grpid)%nr
      if (j1.ge.i1) hvec2(i1:j1) = rmatref(i1:j1,nc)
      if (j2.ge.i2) hvec2(i2-1:j2-1) = rmatref(i2:j2,nc)
      k = 0
      do kk=1,constraints(grpid)%nr
        if (kk.eq.nc) cycle
        k = k + 1
        if (j1.ge.i1) rmat(i1:j1,k) = rmatref(i1:j1,kk)
        if (j2.ge.i2) rmat(i2-1:j2-1,k) = rmatref(i2:j2,kk)
      end do
!     solve
      hvec2 = -hvec2
!      write(*,*) transpose(rmat)
!      write(*,*) hvec2
      rmat2 = rmat
      call dgesv(constraints(grpid)%nr-1,1,rmat,constraints(grpid)%nr-1,rowidx,hvec2,constraints(grpid)%nr-1,lapinf)
!      write(*,*) matmul(rmat2,hvec2)
      if (lapinf.ne.0) then
        write(ilog,*) 'Fatal. LAPACK routine dgesv returned with error status ',lapinf,' in&
 & shake_p_setA(...) for constraint group ',grpid,'.'
        call fexit() 
      end if
      if (j1.ge.i1) amattmp(i1:j1,nc) = 1.0*hvec2(i1:j1)
      if (j2.ge.i2) amattmp(i2:j2,nc) = 1.0*hvec2(i2-1:j2-1)
    end do
!    write(*,10) transpose(amattmp)
!
    vmat = matmul(vmat,amattmp)
    constraints(grpid)%amat = matmul(constraints(grpid)%amat,amattmp)
    maxdev = 0.0
    do nc=1,constraints(grpid)%nr
      do k=1,constraints(grpid)%nr
        if (nc.eq.k) cycle
        if (abs(amattmp(nc,k)).gt.maxdev) maxdev = abs(amattmp(nc,k))
      end do
    end do
!
  end do
!
end
!
#endif
!
!---------------------------------------------------------------------------------------------------------
!
! and the update routine for the pre-conditioner in P-SHAKE (J. Comput. Phys., 220, 740-750 (2007))
! this routine is not yet functional -> not sure of what the mistake is
!
subroutine shake_p_updateA(grpid)
!
  use atoms
  use shakeetal
  use system
  use molecule
!
  implicit none
!
  integer nc,grpid
  RTYPE dt2
  RTYPE dvec(constraints(grpid)%nr,3)
  RTYPE amattmp(constraints(grpid)%nr,constraints(grpid)%nr)
  RTYPE vmat(3*constraints(grpid)%nats,constraints(grpid)%nr)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
  RTYPE denom(3,constraints(grpid)%nr)
!
  dt2 = dyn_dt*dyn_dt
!
  if (constraints(grpid)%nr.eq.1) then
    constraints(grpid)%amat(1,1) = 1.0
    return
  end if
!
  vmat(:,:) = 0.0
  do nc=1,constraints(grpid)%nr
    mterm1(nc) = dt2/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = dt2/mass(constraints(grpid)%idx(nc,2))
    dvec(nc,1) = x(constraints(grpid)%idx(nc,1)) - x(constraints(grpid)%idx(nc,2))
    dvec(nc,2) = y(constraints(grpid)%idx(nc,1)) - y(constraints(grpid)%idx(nc,2))
    dvec(nc,3) = z(constraints(grpid)%idx(nc,1)) - z(constraints(grpid)%idx(nc,2))
    vmat((3*constraints(grpid)%mapidx(nc,1)-2):(3*constraints(grpid)%mapidx(nc,1)),nc) = dvec(nc,:)*mterm1(nc)
    vmat((3*constraints(grpid)%mapidx(nc,2)-2):(3*constraints(grpid)%mapidx(nc,2)),nc) = -dvec(nc,:)*mterm2(nc)
  end do
!
  do nc=1,constraints(grpid)%nr
    amattmp(:,nc) = (vmat(3*(constraints(grpid)%mapidx(nc,1))-2,:)-vmat(3*(constraints(grpid)%mapidx(nc,2))-2,:))*dvec(nc,1) + &
 &                  (vmat(3*(constraints(grpid)%mapidx(nc,1))-1,:)-vmat(3*(constraints(grpid)%mapidx(nc,2))-1,:))*dvec(nc,2) + &
 &                  (vmat(3*(constraints(grpid)%mapidx(nc,1)),:)-vmat(3*(constraints(grpid)%mapidx(nc,2)),:))*dvec(nc,3)
    denom(:,nc) = vmat(3*constraints(grpid)%mapidx(nc,1)-2:3*constraints(grpid)%mapidx(nc,1),nc) - &
 &                vmat(3*constraints(grpid)%mapidx(nc,2)-2:3*constraints(grpid)%mapidx(nc,2),nc)
  end do
  do nc=1,constraints(grpid)%nr
    amattmp(nc,:) = -amattmp(nc,:)/(denom(1,nc)*dvec(:,1) + denom(2,nc)*dvec(:,2) + denom(3,nc)*dvec(:,3))
    amattmp(nc,nc) = 1.0
  end do
  constraints(grpid)%amat = matmul(constraints(grpid)%amat,amattmp)
!
end
!
!-------------------------------------------------------------------------------------------------------------
!
! this is LINCS (J. Comput. Chem., 18 (12), 1463-1472 (1997))
!
subroutine lincs(grpid)
!
  use atoms
  use shakeetal
  use system
  use iounit
  use mcsums
  use forces
!
  implicit none
! 
  integer nc,k,np,grpid,jumpidx,loi,liter
  RTYPE maxdev,idt
  RTYPE dvec(constraints(grpid)%nr,3)
  RTYPE dvecref(constraints(grpid)%nr,3)
  RTYPE d2ref(constraints(grpid)%nr)
  RTYPE id1ref(constraints(grpid)%nr)
  RTYPE mterm(constraints(grpid)%nr)
  RTYPE mterm1(constraints(grpid)%nr)
  RTYPE mterm2(constraints(grpid)%nr)
  RTYPE amat(MAXVAL(constraints(grpid)%ncis),constraints(grpid)%nr)
  RTYPE rvec(constraints(grpid)%nr,3)
  RTYPE sol(constraints(grpid)%nr)
  RTYPE lterm(constraints(grpid)%nr)
!
  idt = 1.0/dyn_dt
  np = 0
  do nc=1,constraints(grpid)%nr
    mterm1(nc) = 1.0/mass(constraints(grpid)%idx(nc,1))
    mterm2(nc) = 1.0/mass(constraints(grpid)%idx(nc,2))
    dvecref(nc,1) = xref(constraints(grpid)%idx(nc,1)) - xref(constraints(grpid)%idx(nc,2))
    dvecref(nc,2) = yref(constraints(grpid)%idx(nc,1)) - yref(constraints(grpid)%idx(nc,2))
    dvecref(nc,3) = zref(constraints(grpid)%idx(nc,1)) - zref(constraints(grpid)%idx(nc,2))
    dvec(nc,1) = x(constraints(grpid)%idx(nc,1)) - x(constraints(grpid)%idx(nc,2))
    dvec(nc,2) = y(constraints(grpid)%idx(nc,1)) - y(constraints(grpid)%idx(nc,2))
    dvec(nc,3) = z(constraints(grpid)%idx(nc,1)) - z(constraints(grpid)%idx(nc,2))
    lterm(nc) = sqrt(constraints(grpid)%eqd2(nc))
  end do
  mterm(:) = sqrt(1.0/(mterm1(:) + mterm2(:)))
  d2ref(:) = dvecref(:,1)*dvecref(:,1) + dvecref(:,2)*dvecref(:,2) + dvecref(:,3)*dvecref(:,3)
  id1ref(:) = sqrt(1.0/d2ref(:))
  do k=1,3
    dvecref(:,k) = dvecref(:,k)*id1ref(:)
  end do
!
  do nc=1,constraints(grpid)%nr
    amat(1:constraints(grpid)%ncis(nc),nc) = constraints(grpid)%ccoeff(nc,1:constraints(grpid)%ncis(nc))* ( &
 &            dvecref(nc,1)*dvecref(constraints(grpid)%cidx(nc,1:constraints(grpid)%ncis(nc)),1) + &
 &            dvecref(nc,2)*dvecref(constraints(grpid)%cidx(nc,1:constraints(grpid)%ncis(nc)),2) + &
 &            dvecref(nc,3)*dvecref(constraints(grpid)%cidx(nc,1:constraints(grpid)%ncis(nc)),3) )
  end do
  rvec(:,1) = mterm(:) * ( (dvecref(:,1)*dvec(:,1) + dvecref(:,2)*dvec(:,2) + dvecref(:,3)*dvec(:,3)) - lterm(:) )
  sol(:) = rvec(:,1)
!
! solve (first time)
  jumpidx = 2
  do loi=1,constraints(grpid)%lorder
    rvec(:,jumpidx) = 0.0
    do nc=1,constraints(grpid)%nr
      rvec(nc,jumpidx) = rvec(nc,jumpidx) + sum( amat(1:constraints(grpid)%ncis(nc),nc) * &
 &                              rvec(constraints(grpid)%cidx(nc,1:constraints(grpid)%ncis(nc)),3-jumpidx) )
    end do
    sol(:) = sol(:) + rvec(:,jumpidx)
    jumpidx = 3 - jumpidx
  end do
! reuse
  do k=1,3
    dvec(:,k) = mterm(:)*sol(:)*dvecref(:,k)
  end do
  do nc=1,constraints(grpid)%nr
    x(constraints(grpid)%idx(nc,1)) = x(constraints(grpid)%idx(nc,1)) - mterm1(nc)*dvec(nc,1)
    y(constraints(grpid)%idx(nc,1)) = y(constraints(grpid)%idx(nc,1)) - mterm1(nc)*dvec(nc,2)
    z(constraints(grpid)%idx(nc,1)) = z(constraints(grpid)%idx(nc,1)) - mterm1(nc)*dvec(nc,3)
    x(constraints(grpid)%idx(nc,2)) = x(constraints(grpid)%idx(nc,2)) + mterm2(nc)*dvec(nc,1)
    y(constraints(grpid)%idx(nc,2)) = y(constraints(grpid)%idx(nc,2)) + mterm2(nc)*dvec(nc,2)
    z(constraints(grpid)%idx(nc,2)) = z(constraints(grpid)%idx(nc,2)) + mterm2(nc)*dvec(nc,3)
  end do
!
  do liter=1,lincs_iter
! correct rotational lengthening
  rvec(:,1) = mterm(:) * (lterm(:) - sqrt( 2.0*constraints(grpid)%eqd2(1:constraints(grpid)%nr) - &
 &             (x(constraints(grpid)%idx(:,1)) - x(constraints(grpid)%idx(:,2)))**2 - &
 &             (y(constraints(grpid)%idx(:,1)) - y(constraints(grpid)%idx(:,2)))**2 - &
 &             (z(constraints(grpid)%idx(:,1)) - z(constraints(grpid)%idx(:,2)))**2 ) )
  sol(:) = rvec(:,1)
!
! solve (second time)
  jumpidx = 2
  do loi=1,constraints(grpid)%lorder
    rvec(:,jumpidx) = 0.0
    do nc=1,constraints(grpid)%nr
      rvec(nc,jumpidx) = rvec(nc,jumpidx) + sum( amat(1:constraints(grpid)%ncis(nc),nc) * &
 &                              rvec(constraints(grpid)%cidx(nc,1:constraints(grpid)%ncis(nc)),3-jumpidx) )
    end do
    sol(:) = sol(:) + rvec(:,jumpidx)
    jumpidx = 3 - jumpidx
  end do
! reuse
  do k=1,3
    dvec(:,k) = mterm(:)*sol(:)*dvecref(:,k)
  end do
  do nc=1,constraints(grpid)%nr
    x(constraints(grpid)%idx(nc,1)) = x(constraints(grpid)%idx(nc,1)) - mterm1(nc)*dvec(nc,1)
    y(constraints(grpid)%idx(nc,1)) = y(constraints(grpid)%idx(nc,1)) - mterm1(nc)*dvec(nc,2)
    z(constraints(grpid)%idx(nc,1)) = z(constraints(grpid)%idx(nc,1)) - mterm1(nc)*dvec(nc,3)
    x(constraints(grpid)%idx(nc,2)) = x(constraints(grpid)%idx(nc,2)) + mterm2(nc)*dvec(nc,1)
    y(constraints(grpid)%idx(nc,2)) = y(constraints(grpid)%idx(nc,2)) + mterm2(nc)*dvec(nc,2)
    z(constraints(grpid)%idx(nc,2)) = z(constraints(grpid)%idx(nc,2)) + mterm2(nc)*dvec(nc,3)
  end do
  end do
!
! now readjust cartesian velocities
! reconstructed velocities (override!) for the half-step t + 1/2dt
  do nc=1,constraints(grpid)%nats
    cart_v(constraints(grpid)%uidx(nc),1) = (x(constraints(grpid)%uidx(nc)) - xref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),2) = (y(constraints(grpid)%uidx(nc)) - yref(constraints(grpid)%uidx(nc)))*idt
    cart_v(constraints(grpid)%uidx(nc),3) = (z(constraints(grpid)%uidx(nc)) - zref(constraints(grpid)%uidx(nc)))*idt
  end do
!
! check accuracy and eventually warn/adjust
  do nc=1,constraints(grpid)%nr
    dvec(nc,1) = x(constraints(grpid)%idx(nc,1)) - x(constraints(grpid)%idx(nc,2))
    dvec(nc,2) = y(constraints(grpid)%idx(nc,1)) - y(constraints(grpid)%idx(nc,2))
    dvec(nc,3) = z(constraints(grpid)%idx(nc,1)) - z(constraints(grpid)%idx(nc,2))
  end do
  d2ref(:) = abs(dvec(:,1)*dvec(:,1) + dvec(:,2)*dvec(:,2) + dvec(:,3)*dvec(:,3) - &
 &    constraints(grpid)%eqd2(1:constraints(grpid)%nr) ) / constraints(grpid)%eqd2(1:constraints(grpid)%nr)
!
  maxdev = 0.0
  do nc=1,constraints(grpid)%nr
    if (d2ref(nc).gt.maxdev) maxdev = d2ref(nc)
  end do
!
  constraints(grpid)%iters = constraints(grpid)%iters + constraints(grpid)%lorder
!
  if (maxdev.gt.shake_tol) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(LINCS_WARNINGS)
#endif
    cs_wrncnt(19) = cs_wrncnt(19) + 1
    if (cs_wrncnt(19).eq.cs_wrnlmt(19)) then
      write(ilog,*) 'Warning. Maximum bond length deviation remains large after LINCS termination. This &
 &may happen frequently if too small a tolerance was chosen or during rapid equilibration.'
      write(ilog,*) 'This was warning #',cs_wrncnt(19),' of this type not all of which may be displayed.'
      cs_wrnlmt(19) = cs_wrnlmt(19)*10
    end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(LINCS_WARNINGS)
#endif
    constraints(grpid)%lorder = min(constraints(grpid)%lorder + 1,16)
  else if (maxdev.lt.0.1*shake_tol) then
    constraints(grpid)%lorder = max(constraints(grpid)%lorder - 1,2)
  end if
!
  constraints(grpid)%itercnt = constraints(grpid)%itercnt + 1
!
end
!
!----------------------------------------------------------------------------------------------------
!
! when using massless (virtual) sites, the forces acting on them must be redistributed to those atoms
! the sites are rigidly connected to
! this routine is here, because it is only needed in Cartesian MD, and because Cartesian MD with virtual
! sites absolutely requires SHAKE et al.
! the strategy is to find compensating internal forces (bond stretch, angle bend, torsional) for the Z-matrix
! entries defining the virtual site, that - when added successively - exactly zero out the force on the site
! and add additional forces on iz(1,k)-iz(3,k)
! because of the dependence on Z-matrix structure, the solution is not unique, but hopefully reasonable
!
subroutine cart2cart_shake(tpi)
!
  use shakeetal
  use iounit
  use forces
  use atoms
  use zmatrix
  use math
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
! 
  integer, INTENT(IN):: tpi
  integer i,k,kk,jj,sta,sto
  RTYPE or_pl(3),nopl,in_pl(3),nipl,dax,dox,cosine,sine,alph,netm
  RTYPE torstar,con1(3),c1,con2(3),c2,con3(3),c3,getbang,getztor,trigs(2),tdvec(3,4)
!
#ifdef ENABLE_THREADS
  if (tpi.eq.0) then
    sta = 1
    sto = cart_cons_grps
  else
    sta = thr_limits(25,tpi)
    sto = thr_limits(26,tpi)
  end if
#else
  sta = 1
  sto = cart_cons_grps
#endif
  do i=sta,sto
    if (constraints(i)%hasmassless.EQV..false.) cycle
    do kk=constraints(i)%nats,1,-1
      k = constraints(i)%uidx(kk)
      if (mass(k).gt.0.0) cycle
!     first get the Z-matrix relevant bond (or pseudobond) vectors, their norms and in-plane and out-of-plane vec.s
      con2(1) = x(k)-x(iz(1,k))
      con2(2) = y(k)-y(iz(1,k))
      con2(3) = z(k)-z(iz(1,k))
      c2 = sqrt(dot_product(con2,con2))
      con2(:) = con2(:)/c2
      alph = getbang(k,iz(1,k),iz(2,k))/RADIAN
      dax = sin(alph)*c2
      dox = cos(alph)*c2
      torstar = getztor(iz(3,k),iz(2,k),iz(1,k),k)/RADIAN
      con1(1) = x(iz(1,k)) - x(iz(2,k))
      con1(2) = y(iz(1,k)) - y(iz(2,k))
      con1(3) = z(iz(1,k)) - z(iz(2,k))
      c1 = sqrt(dot_product(con1,con1))
      con1(:) = con1(:)/c1
      con3(1) = x(iz(2,k)) - x(iz(3,k))
      con3(2) = y(iz(2,k)) - y(iz(3,k))
      con3(3) = z(iz(2,k)) - z(iz(3,k))
      c3 = sqrt(dot_product(con3,con3))
      con3(:) = con3(:)/c3
      call crossprod(con3,con1,or_pl,nopl)
      cosine = dot_product(con3,con1)
      sine = sqrt(max(1.0d0-cosine**2,1.0d-7))
      if (abs(cosine) .ge. 1.0d0) then
        write(ilog,10)  i,cosine
 10     format (/,' Bad torsion setup for virtual site in constraint group',i6,': ',f12.4)
      end if
      or_pl(:) = or_pl(:)/sine
      call crossprod(or_pl,con1,in_pl,nipl)
!     now the actual corrections; the strategy is always to compute dx,y,z/dq, dot it with dU/dx,dy,dz (site) to
!     get a net force on q, and to then distribute a force of opposite sign via dU/dq * dq/dx,y,z onto the
!     max. 4 participating atoms
!     first bond stretch
!      netm = dot_product((dax*(sin(torstar)*or_pl(:) + cos(torstar)*in_pl(:)) - dox*con1(:))/c2,cart_f(k,:))
      netm = -dot_product(con2,cart_f(k,:))
!      call onebond_deriv(k,iz(1,k),tdvec(:,1),tdvec(:,2))
!      do jj=1,1
!        cart_f(iz(jj,k),:) = cart_f(iz(jj,k),:) + netm*(tdvec(:,jj+1))
!      end do
!      cart_f(k,:) = cart_f(k,:) - netm*(tdvec(:,1))
      cart_f(iz(1,k),:) = cart_f(iz(1,k),:) - netm*con2(:)
      cart_f(k,:) = cart_f(k,:) + netm*con2(:)
!     now angle bending
      netm = -dot_product(dox*(sin(torstar)*or_pl(:) + cos(torstar)*in_pl(:)) + dax*con1(:),cart_f(k,:))
      call oneang_deriv(k,iz(1,k),iz(2,k),trigs(1),tdvec(:,1),tdvec(:,2),tdvec(:,3))
      do jj=1,2
        cart_f(iz(jj,k),:) = cart_f(iz(jj,k),:) + netm*(tdvec(:,jj+1))/RADIAN
      end do
      cart_f(k,:) = cart_f(k,:) + netm*(tdvec(:,1))/RADIAN
!     and finally torsional
      netm = -dot_product(dax*(cos(torstar)*or_pl(:) - sin(torstar)*in_pl(:)),cart_f(k,:))
      call onetor_deriv(k,iz(1,k),iz(2,k),iz(3,k),trigs(1),trigs(2),tdvec(:,1),tdvec(:,2),tdvec(:,3),tdvec(:,4))
      do jj=1,3
        cart_f(iz(jj,k),:) = cart_f(iz(jj,k),:) + netm*(tdvec(:,jj+1))/RADIAN
      end do
      cart_f(k,:) = cart_f(k,:) + netm*(tdvec(:,1))/RADIAN
      if (dot_product(cart_f(k,:),cart_f(k,:)).gt.1.0e-9) then
#ifdef ENABLE_THREADS
!$OMP CRITICAL(VIRTUALSITE_WARNINGS)
#endif
        cs_wrncnt(18) = cs_wrncnt(18) + 1
        if (cs_wrncnt(18).eq.cs_wrnlmt(18)) then
          write(ilog,*) 'Warning. Net force from virtual site could not be remapped to atoms used to construct &
 &site. This can happen if the Z-matrix is improperly set up or geometries are simplisitc (colinear bonds). It may also &
 &indicate a bug.'
          write(ilog,*) 'This was warning #',cs_wrncnt(18),' of this type not all of which may be displayed.'
          cs_wrnlmt(18) = cs_wrnlmt(18)*10
        end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(VIRTUALSITE_WARNINGS)
#endif
      end if
    end do
  end do
!
end
!
!--------------------------------------------------------------------------------------------------------
!
