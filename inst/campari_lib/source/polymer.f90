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
! CONTRIBUTIONS: Rohit Pappu, Hoang Tran                                   !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!
! #######################################################
! ##                                                   ##
! ## subroutine polymer -- record polymeric properties ##
! ##                                                   ##
! #######################################################
!
!
! in this instant. write-out, we provide system-wide quantities which
! are only occassionally meaningful, but complement other
! metrics of the overall simulation (like MC acceptance, energies, dihedrals, ensemble var.s)
!
subroutine polymer()
!
  use mcsums
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  RTYPE rgs,rgsten(3,3),rgsev(3)
  RTYPE asphs,acyls,dlts,dltsts,tps
!
 10   format(i12,2x,g14.7,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5,1x,g14.7,&
 &1x,g14.7,1x,g14.7)
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
!
    call gyrate2(rgs,tps,rgsten,rgsev,asphs,acyls,dlts,dltsts)
!
    write(ipolmr,10) nstep,rgs,tps,asphs/(rgs**2),acyls/(rgs**2),&
 &dltsts,rgsev(3),rgsev(2),rgsev(1)
!
  else if (use_MPIAVG.EQV..true.) then
    if (myrank.eq.0) then
      call MPI_AVGWritePol()
    else
      call MPI_AVGSendPol()
    end if
    if (mpi_cnt_pol.eq.(mpi_nodes-1)) then
      mpi_cnt_pol = 0
    else
      mpi_cnt_pol = mpi_cnt_pol + 1
    end if
  end if
#else
!
  call gyrate2(rgs,tps,rgsten,rgsev,asphs,acyls,dlts,dltsts)
!
  write(ipolmr,10) nstep,rgs,tps,asphs/(rgs**2),acyls/(rgs**2),&
 &dltsts,rgsev(3),rgsev(2),rgsev(1)
#endif
!
end
!
! #######################################################
! ##                                                   ##
! ## subroutine do_poly -- average polymeric properties ##
! ##                                                   ##
! #######################################################
!
!
subroutine do_poly()
!
  use molecule
  use sequen
  use polyavg
  use iounit
  use mcsums
  use grandensembles
  use system
  use pdb
!
  implicit none
!
  RTYPE t1n,t2n,t12n,incr
  integer j,i,imol,k,d1,d2
  logical ismember,doperser
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  if (((dyn_mode.ge.2).AND.(dyn_mode.le.4))&
 &.OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    do imol=1,nmol
      call update_rigid(imol)
    end do
  end if
!
  polynew%vt(:) = 0.0
  polynew%rg(:) = 0.0
  polynew%rg2(:) = 0.0
  polynew%rgt1(:) = 0.0
  polynew%rgt2(:) = 0.0
  polynew%rgt12(:) = 0.0
  polynew%rgev(:,:) = 0.0
  polynew%map(:) = npolavg(:)
!
! all the gyration tensor-dependent statistics are automatically available
! through the rigid-body coordinates
  doperser = .false.
  do imol=1,nmol
!
!   set doperser (should be global)
    if (do_pers(moltypid(imol)).EQV..true.) doperser = .true.
!   cycle out for inappropriate molecules
    if (do_pol(moltypid(imol)).EQV..false.) cycle
!
    k = an_grp_mol(imol)
!
!  cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
!   increment counter
    npolavg(k) = npolavg(k) + 1
    npersavg(k) = npersavg(k) + 1
!
!   average Rg, principal components of Rg-tensor (ascending order), and some combinations
    polynew%rg(k) = polynew%rg(k) + rgv(imol)
    polynew%rg2(k) = polynew%rg2(k) + rgv(imol)*rgv(imol)
    polynew%vt(k) = polynew%vt(k) + 2.5*(1.75*rgv(imol)/&
 &                             molcontlen(moltypid(imol)))&
 &          **(4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)
    do i=1,3
      polynew%rgev(k,i) = polynew%rgev(k,i) + rgevs(imol,i)
    end do
    t1n = rgevs(imol,1)*rgevs(imol,2) + rgevs(imol,2)*rgevs(imol,3) &
 &                + rgevs(imol,1)*rgevs(imol,3)
    t2n = (rgevs(imol,1) + rgevs(imol,2) + rgevs(imol,3))**2
    t12n = t1n/t2n
    polynew%rgt1(k) = polynew%rgt1(k) + t1n
    polynew%rgt2(k) = polynew%rgt2(k) + t2n
    polynew%rgt12(k) = polynew%rgt12(k) + t12n
!
!   take care of the t vs. delta* histogram (which of course also contains the t-histogram)
!   note that t is pseudo-normalized [0:~2.0], and that the hard-coded array size might have to
!   be re-evaluated for lower-density systems
    d1 = max(1,min(floor(2.5*(1.75*rgv(imol)/molcontlen(moltypid(imol)))&
 &          **(4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)/tp_binsz) + 1,TPBINS))
    d2 = max(1,min(floor((1.0 - 3.0*(t12n))/dlst_binsz) + 1,DLSTBINS))
!   in case of a (hopefully rare) range exception, mercilessly overstock last bin
    rdhist(k,d1,d2) = rdhist(k,d1,d2) + incr
!   take care of the simple Rg-histogram
    d1 = max(1,min(floor(rgv(imol)/rg_binsz) + 1,maxrgbins))
!   in case of range exception, mercilessly overstock last bin
    rghist(k,d1) = rghist(k,d1) + incr
  end do
!
  if (use_frame_weights.EQV..true.) then
    polynew%vt(:) = incr*polynew%vt(:)
    polynew%rg(:) = incr*polynew%rg(:)
    polynew%rg2(:) = incr*polynew%rg2(:) 
    polynew%rgt1(:) = incr*polynew%rgt1(:)
    polynew%rgt2(:) = incr*polynew%rgt2(:)
    polynew%rgt12(:) = incr*polynew%rgt12(:)
    polynew%rgev(:,:) = incr*polynew%rgev(:,:)
  end if
!
! end-to-end distance (note that this takes care of Ree-histograms)
  call rete(polynew%ree)
!
! density profiles
  call do_density_profile()
!
! now we can construct the running averages 
  do i=1,nangrps
    if (npolavg(i).le.0) cycle
    rgavg(i) = (polynew%map(i)*rgavg(i) + polynew%rg(i))/(1.0*npolavg(i))
    rg2avg(i) = (polynew%map(i)*rg2avg(i) + polynew%rg2(i))/(1.0*npolavg(i))
    vtavg(i) = (polynew%map(i)*vtavg(i) + polynew%vt(i))/(1.0*npolavg(i))
    do j=1,3
      rgpcavg(i,j) = (polynew%map(i)*rgpcavg(i,j) + polynew%rgev(i,j))/(1.0*npolavg(i))
    end do
    rgterm1avg(i) = (polynew%map(i)*rgterm1avg(i) + polynew%rgt1(i))/(1.0*npolavg(i))
    rgterm2avg(i) = (polynew%map(i)*rgterm2avg(i) + polynew%rgt2(i))/(1.0*npolavg(i))
    rgterm12avg(i) = (polynew%map(i)*rgterm12avg(i) + polynew%rgt12(i))/(1.0*npolavg(i))
    dlstavg(i) = 1.0 - 3.0*rgterm12avg(i)
!   this condition needs to be checked in case weights are used and the very first weight is zero
    if (rgterm2avg(i).gt.0.0) dlavg(i) = 1.0 - 3.0*(rgterm1avg(i)/rgterm2avg(i))
    asphavg(i) = rgpcavg(i,3) - 0.5*(rgpcavg(i,2) + rgpcavg(i,1))
    acylavg(i) = rgpcavg(i,2) - rgpcavg(i,1)
    reteavg(i) = (polynew%map(i)*reteavg(i) + polynew%ree(i))/(1.0*npolavg(i))
  end do
!
  if (doperser.EQV..true.) then
!   cumulative average of angular correlation as fxn of distance in sequence
!   allows computation of persistence length
    call persis(polynew%cosa)
!    if (inst_pers) then
!       write(ipers,70)(polynew%cosa(i),i=1,nseq)
!   70      format(500f8.4)
!    end if
    do i=1,nangrps
      do j=1,rsmol(molangr(i,1),2)-rsmol(molangr(i,1),1)+1
        persavg(i,j)=persavg(i,j)+polynew%cosa(i,j)
      end do
    end do
  end if
!
end
!
!
!------------------------------------------------------------------------
!
! a mini-routine necessary because Rh-computation scales with N^2 and
! should not be incorporated into the standard polymer stuff, all of
! which is just ~N
!
subroutine hydroradavg()
!
  use iounit
  use polyavg
  use molecule
  use sequen
  use mcsums
  use system
  use grandensembles
  use pdb
!
  implicit none
!
  integer i,j,k,nrhc,imol
  logical dosct,ismember
  RTYPE incr
!
  nrhc = nstep/rhcalc
!
  if (mod(nrhc,sctcalc).eq.0) then
    dosct = .true.
  else
    dosct = .false.
  end if
!
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  do imol=1,nmol
!
!   cycle out for inappropriate molecules
    if (do_pol(moltypid(imol)).EQV..false.) cycle
!
    k = an_grp_mol(imol)
!
!  cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
!   increment counter
    nrhavg(k) = nrhavg(k) + 1
    if (dosct.EQV..true.) nsctavg(k) = nsctavg(k) + 1
!
  end do
!
  call hydroradius(polynew%rhr,polynew%rh,polynew%psc,dosct)
!
  do i=1,nangrps
    rhravg(i) = rhravg(i) + incr*polynew%rhr(i)
    do j=1,rsmol(molangr(i,1),2)-rsmol(molangr(i,1),1)+1
      rhavg(i,j) = rhavg(i,j) + incr*polynew%rh(i,j)
    end do
    if (dosct.EQV..true.) then
      do j=1,nqv
        pscavg(i,j) = pscavg(i,j) + incr*polynew%psc(i,j)
      end do
    end if
  end do
!
end
!
!-------------------------------------------------------------------------
!
!
! #############################################
! ##                                         ##
! ## subroutine rete -- end-to-end distances ##
! ##                                         ##
! #############################################
!
!
! "rete" computes the end-to-end distance for polypeptides
!
subroutine rete(ree) 
!
  use polypep
  use sequen
  use polyavg
  use molecule
  use aminos
  use iounit
  use grandensembles
  use system
  use pdb
!
  implicit none
!
  RTYPE getblen,ree(nangrps),reenew,incr
  integer refi,refj,imol,d1
  logical ismember
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  ree(:) = 0.0
!
  do imol=1,nmol
    if (rsmol(imol,1).eq.rsmol(imol,2)) cycle
    if (do_pol(moltypid(imol)).EQV..false.) cycle
!   cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
    if (seqpolty(rsmol(imol,1)).eq.'P') then
      refi = cai(rsmol(imol,1))
    else if (seqpolty(rsmol(imol,1)).eq.'N') then
      refi = nuci(rsmol(imol,1),2)
    else if ((amino(seqtyp(rsmol(imol,1))).eq.'ACE').OR.&
 &           (amino(seqtyp(rsmol(imol,1))).eq.'FOR')) then
      refi = at(rsmol(imol,1))%bb(1)
    else if (amino(seqtyp(rsmol(imol,1))).eq.'UNK') then
      refi = refat(rsmol(imol,1))
    else
      write(ilog,*) 'Fatal. Called rete(...) with unsupported residu&
 &e-type (first). This is most likely a minor bug during setup.'
      call fexit()
    end if
!
    if (seqpolty(rsmol(imol,2)).eq.'P') then
      refj = cai(rsmol(imol,2))
    else if (seqpolty(rsmol(imol,2)).eq.'N') then
      refj = nuci(rsmol(imol,2),6)
    else if ((amino(seqtyp(rsmol(imol,2))).eq.'NME').OR.&
 &           (amino(seqtyp(rsmol(imol,2))).eq.'NH2')) then
      refj = at(rsmol(imol,2))%bb(1)
    else if (amino(seqtyp(rsmol(imol,2))).eq.'UNK') then
      refj = refat(rsmol(imol,2))
    else
      write(ilog,*) 'Fatal. Called rete(...) with unsupported residu&
 &e-type (last). This is most likely a minor bug during setup.'
      call fexit()
    end if
!
    reenew = getblen(refi,refj)
!    chle = sqrt(dble(rsmol(imol,2)-rsmol(imol,1)+1))
    d1 = max(1,min(floor(reenew/rg_binsz) + 1,maxrgbins))
    if (d1.gt.maxrgbins) d1 = maxrgbins
    retehist(an_grp_mol(imol),d1) = retehist(an_grp_mol(imol),d1) + incr
    ree(an_grp_mol(imol)) = ree(an_grp_mol(imol)) + incr*reenew
  end do
!
end
!
!
!
! #############################################
! ##                                         ##
! ## subroutine persis -- persistence length ##
! ##                     via correlation fxn ##
! ##                                         ##
! #############################################
!
!
! just computes the average angle between the CA->CO
! vectors going through the sequence as a function of distance
! in sequence
!
subroutine persis(cosavg) 
!
  use iounit
  use atoms
  use math
  use polypep
  use sequen
  use polyavg
  use molecule
  use system
  use grandensembles
  use aminos
  use pdb
!
  implicit none
!
  integer i,j,k,lastrs,firstrs,imol,rsref,mt,rei1,rei2,rej1,rej2
  RTYPE cosine,cosavg(nangrps,longestmol),xva,yva,zva,xvb,yvb,zvb,rva2,rvb2,dot
  RTYPE ang,dih,getztor,incr
  logical ismember
!
  cosavg(:,:) = 0.0
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  do imol=1,nmol
!
!   cycle out if this molecule type does not support persistence length analysis
    if (do_pers(moltypid(imol)).EQV..false.) cycle
!
    mt = an_grp_mol(imol)
!
!   cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
    do i=1,rsmol(imol,2)-rsmol(imol,1)+1
      polynew%cost(i) = 0.0
    end do
!
    if ((ci(rsmol(imol,1)).gt.0).AND.(ni(rsmol(imol,1)).gt.0)) then
      firstrs = rsmol(imol,1)
    else if ((nuci(rsmol(imol,1),2).gt.0).AND.&
 &      (nuci(rsmol(imol,1),6).gt.0)) then
      firstrs = rsmol(imol,1)
    else
      firstrs = rsmol(imol,1) + 1
    end if
!    
    if ((ci(rsmol(imol,2)).gt.0).AND.(ni(rsmol(imol,2)).gt.0)) then
      lastrs = rsmol(imol,2)
    else if ((nuci(rsmol(imol,2),2).gt.0).AND.&
 &      (nuci(rsmol(imol,2),6).gt.0)) then
      lastrs = rsmol(imol,2)
    else
      lastrs = rsmol(imol,2)-1
    end if
!
    if (lastrs.le.firstrs) then
      write(ilog,*) 'Fatal. Angular correlation analysis was called with an &
 &incompatible/incompetent polymer. This is a setup bug in parse_sequence(...).'
      call fexit()
    end if
!
    do i=firstrs,lastrs-1
      if (firstrs.eq.rsmol(imol,1)) then
        rsref = rsmol(molangr(mt,1),1) + (i-firstrs)
      else
        rsref = rsmol(molangr(mt,1),1) + 1 + (i-firstrs)
      end if
      do j=(i+1),lastrs
!       get the two N->C-vectors and their norms
        if (ni(i).gt.0) then ! parse_sequence(...) makes sure that this is safe
          rei2 = ci(i)
          rei1 = ni(i)
          rej2 = ci(j)
          rej1 = ni(j)
        else if (nuci(i,2).gt.0) then
          rei2 = nuci(i,6)
          rei1 = nuci(i,2)
          rej2 = nuci(j,6)
          rej1 = nuci(j,2)
        end if
        xva = x(rei2) - x(rei1)
        yva = y(rei2) - y(rei1)
        zva = z(rei2) - z(rei1)
        rva2 = xva*xva + yva*yva + zva*zva
        xvb = x(rej2) - x(rej1)
        yvb = y(rej2) - y(rej1)
        zvb = z(rej2) - z(rej1)
        rvb2 = xvb*xvb + yvb*yvb + zvb*zvb
!       compute the cosine
        if ((rva2.gt.0.0).AND.(rvb2.gt.0.0)) then
          dot = xva*xvb + yva*yvb + zva*zvb
          cosine = dot/sqrt(rva2*rvb2)
          polynew%cost(j-i) = polynew%cost(j-i) + cosine
          ang = acos(max(-1.0,min(1.0,cosine)))
          dih = getztor(rei1,rei2,rej1,rej2)
          if (ang.lt.0.0) then
            ang = -ang
          end if
!         turn prop.s are accumulated residue-specific (first res. as ref. res.)
!         now specific per-molecule parsing necessary
          do k=1,3
            if ((j-i).eq.k) then
              turns_rs(rsref,3*k-2) = turns_rs(rsref,3*k-2) + incr*ang
              turns_rs(rsref,3*k-1) = turns_rs(rsref,3*k-1) + incr*cos(dih/RADIAN)
              turns_rs(rsref,3*k) = turns_rs(rsref,3*k) + incr*sin(dih/RADIAN)
            end if
          end do
        else
          write(ilog,*) 'Fatal: Backbone directional vector for residue ',i,' or ',j,' has zero length!'
          call fexit()
        end if
      end do
    end do
!   normalize
    do j=1,lastrs-firstrs
      polynew%cost(j) = polynew%cost(j)/dble(lastrs-firstrs+1-j)
    end do
!   now transfer into allocated molecule slot
!   value @first res. is value for seq. distance of zero
!   following values are straightforward, molecule always starts @res1
    cosavg(mt,1) = cosavg(mt,1) + incr*1.0
    do j=2,lastrs-firstrs+1
      cosavg(mt,j) = cosavg(mt,j) + incr*polynew%cost(j-1)
    end do
  end do
!
end
!
!
!-----------------------------------------------------------------------------
!
subroutine do_density_profile()
!
  use iounit
  use molecule
  use atoms
  use polyavg
  use sequen
  use system
  use grandensembles
  use pdb
!
  implicit none
!
  integer i,k,d2,imol,mt
  RTYPE cmass(3),dis,tmass,incr
  logical ismember
!
  polynew%dep(:,:) = 0.0
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
  do imol=1,nmol
!
!   cycle out if this molecule type does not support density profile analysis
    if (do_pol(moltypid(imol)).EQV..false.) cycle
!
    mt = an_grp_mol(imol)
!
!   cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
!   get the molecular center of mass
    do k=1,3
      cmass(k) = 0.0
    end do
    tmass = 0.0
    do i=atmol(imol,1),atmol(imol,2)
      tmass = tmass + mass(i)
      cmass(1) = cmass(1) + mass(i)*x(i)
      cmass(2) = cmass(2) + mass(i)*y(i)
      cmass(3) = cmass(3) + mass(i)*z(i)
    end do
    do k=1,3
      cmass(k) = cmass(k) / tmass
    end do
!
!   now loop through the atoms and bin them according to the distance to the center of mass
    do i=atmol(imol,1),atmol(imol,2)
      dis = sqrt( (x(i) - cmass(1))**2 + &
 &                (y(i) - cmass(2))**2 +&
 &                (z(i) - cmass(3))**2)
      d2 = max(1,min(floor(dis/rg_binsz) + 1,maxrgbins))
!     in case of range exception, mercilessly overstock last bin
      polynew%dep(mt,d2) = polynew%dep(mt,d2) + mass(i)
    end do
!   
  end do
!
  densproavg(:,:) = densproavg(:,:) + incr*polynew%dep(:,:)
!
end
!
!-----------------------------------------------------------------------------
!
subroutine prt_pers()
!
  use iounit
  use sequen
  use math
  use molecule
  use polyavg
  use mpistuff
  use mcsums
  use polypep
  use aminos
  use system
  use pdb
!
  implicit none
!
  integer i,k,it,freeunit,persdis,mt,ii,jj
  integer lastrs,firstrs,sss
  character(60) fn
#ifdef ENABLE_MPI
  character(3) xpont
  integer lext2
#endif
  logical exists
  RTYPE sssterm,dihi(3)
  RTYPE, ALLOCATABLE:: normer(:)
!
  exists = .false.
  do i=1,nmol
    if (do_pers(moltypid(i)).EQV..true.) then
      exists = .true.
      exit
    end if
  end do
  if (exists.EQV..false.) return
!
  if (sum(npersavg(:)).le.0) return
!
  allocate(normer(nangrps))
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),polcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)*molangr(:,2)
      end if
    end do
    if (k.ne.npersavg(1)/molangr(1,2)) then
      write(ilog,*) 'Warning. Generic topological analyses have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in PERSISTENCE.dat, TURNS_RES.dat, and &
 &related files.'
    end if
  else
    normer(:) = 1.0*npersavg(:)
  end if
!
  if (inst_pers) then
!   do nothing if instantaneous values have already been printed out
  else
!
 24   format('# Analysis group ',i5,' (Ref.: Mol. ',i6,' / Res. ',i6,'-',&
 &i6,'):')
 89   format(i10,' ',g14.7,1x,g14.7)
!
    do mt=1,nangrps
      if (do_pers(moltypid(molangr(mt,1))).EQV..false.) cycle
      write(ipers,24) mt,molangr(mt,1),rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
      persdis = rsmol(molangr(mt,1),2)-rsmol(molangr(mt,1),1)+1
!
      if (ni(rsmol(molangr(mt,1),1)+1).gt.0) then ! hierarchy of detection -> see parse_sequence(...)
        if ((ni(rsmol(molangr(mt,1),1)).le.0).OR.(ci(rsmol(molangr(mt,1),1)).le.0)) persdis = persdis - 1
        if ((ni(rsmol(molangr(mt,1),2)).le.0).OR.(ci(rsmol(molangr(mt,1),2)).le.0)) persdis = persdis - 1
      else if (nuci(rsmol(molangr(mt,1),1)+1,2).gt.0) then
        if ((nuci(rsmol(molangr(mt,1),1),2).le.0).OR.(nuci(rsmol(molangr(mt,1),1),6).le.0)) persdis = persdis - 1
        if ((nuci(rsmol(molangr(mt,1),2),2).le.0).OR.(nuci(rsmol(molangr(mt,1),2),6).le.0)) persdis = persdis - 1
      end if
!
      if (npersavg(mt).gt.0) then
        do i=1,persdis
          persavg(mt,i) = persavg(mt,i)/normer(mt)
        end do
      end if
      do i=1,persdis
        sss = i-1
        sssterm = persavg(mt,i)*((1.0*sss)**1.5)
        write(ipers,89) sss,persavg(mt,i),sssterm
      end do
    end do
  end if
!
 94   format(i10,6(1x,g14.7))
#ifdef ENABLE_MPI
  lext2 = 3
  if (use_REMC.EQV..true.) then
    call int2str(myrank,xpont,lext2)
    fn =  'N_'//xpont(1:lext2)//'_TURNS_RES.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'TURNS_RES.dat'
  end if
#else
  fn = 'TURNS_RES.dat'
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
!
  do mt=1,nangrps
    if (do_pers(moltypid(molangr(mt,1))).EQV..false.) cycle
!
    firstrs = rsmol(molangr(mt,1),1) 
    lastrs = rsmol(molangr(mt,1),2)
    if (ni(rsmol(molangr(mt,1),1)+1).gt.0) then ! hierarchy of detection -> see parse_sequence(...)
      if ((ni(rsmol(molangr(mt,1),1)).le.0).OR.(ci(rsmol(molangr(mt,1),1)).le.0)) firstrs = firstrs + 1
      if ((ni(rsmol(molangr(mt,1),2)).le.0).OR.(ci(rsmol(molangr(mt,1),2)).le.0)) lastrs = lastrs - 1
    else if (nuci(rsmol(molangr(mt,1),1)+1,2).gt.0) then
      if ((nuci(rsmol(molangr(mt,1),1),2).le.0).OR.(nuci(rsmol(molangr(mt,1),1),6).le.0)) firstrs = firstrs + 1
      if ((nuci(rsmol(molangr(mt,1),2),2).le.0).OR.(nuci(rsmol(molangr(mt,1),2),6).le.0)) lastrs = lastrs - 1
    end if
!
    if (npersavg(mt).gt.0) then
      write(it,24) mt,molangr(mt,1),rsmol(molangr(mt,1),1),rsmol(molangr(mt,1),2)
      do i=firstrs,lastrs-1
        do k=1,3
          turns_rs(i,3*k-2)=radian*turns_rs(i,3*k-2)/normer(mt)
          turns_rs(i,3*k-1)=turns_rs(i,3*k-1)/normer(mt)
          turns_rs(i,3*k) = turns_rs(i,3*k)/normer(mt)
        end do
      end do
      do i=firstrs,lastrs-1
        do k=1,3
          dihi(k) = 0.0
          if ((i+k).gt.lastrs) cycle
          if (turns_rs(i,3*k-1).ge.0.0) then
            dihi(k) = RADIAN*atan(turns_rs(i,3*k)/turns_rs(i,3*k-1))
          else if (turns_rs(i,3*k).gt.0.0) then
            dihi(k)=RADIAN*(atan(turns_rs(i,3*k)/turns_rs(i,3*k-1))+PI)
          else
            dihi(k)=RADIAN*(atan(turns_rs(i,3*k)/turns_rs(i,3*k-1))-PI)
          end if
        end do
        write(it,94) i,(turns_rs(i,3*k-2),dihi(k),k=1,3)
      end do
    end if
  end do
!
  deallocate(normer)
!
end
!
!---------------------------------------------------------------------
!
subroutine prt_rdhist()
!
  use iounit
  use molecule
  use math
  use polyavg
  use mpistuff
  use units
  use system
  use pdb
!
  implicit none
!
  integer i,j,k,iu,freeunit,maxi,effmoltyp,ii,jj
  RTYPE tot,rval
  character(100) fn
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  logical exists
  RTYPE, ALLOCATABLE:: normer(:),volis(:),normer2(:)
!
 24   format('# Molecule type ',i5,' (Ref.: Mol. ',i6,' / Res. ',i6,'-',i6,'):')
 25   format(100(g14.7,1x))
 26   format('# Mol. type: ',100(i7,8x))
 27   format(f12.5,1x,100(g14.7,1x))
 28   format(12(g14.7,1x))
 29   format('# Analy. grp:',100(i7,8x))
 30   format('# Analysis group ',i5,' (Ref.: Mol. ',i6,' / Res. ',i6,'-',i6,'):')
!
  allocate(normer(nangrps))
  allocate(normer2(nangrps))
  allocate(volis(maxrgbins))
!
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),polcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)*molangr(:,2)
      end if
    end do
    if (k.ne.npolavg(1)/molangr(1,2)) then
      write(ilog,*) 'Warning. Simple (linear) polymeric analyses have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in POLYAVG.dat, RDHIST.dat, DENSPROF.dat, and &
 &related files.'
    end if
  else
    normer(:) = 1.0*npolavg(:)
  end if
!
  if (allocated(nrhavg).EQV..true.) then
    if (use_frame_weights.EQV..true.) then
      normer2(:) = 0.0
      k = 0
      do i=1,framecnt
        if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),rhcalc).eq.0)) then
          k = k + 1
          normer2(:) = normer2(:) + framewts2(i)*molangr(:,2)
        end if
      end do
      if (k.ne.nrhavg(1)/molangr(1,2)) then
        write(ilog,*) 'Warning. Pairwise, polymeric analyses have inconsistent normalization. This is &
   &a bug. Please report this error and ignore all output in INTSCAL.dat, KRATKY.dat, and &
   &related files.'
      end if
    else
      normer2(:) = 1.0*nrhavg(:)
    end if
  end if
!
  maxi = 1
  effmoltyp = 0
  do k=1,nangrps
    if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
    tot = 0.0
    if (allocated(rdhist).EQV..true.) then
      tot = sum(rdhist(k,1:TPBINS,1:DLSTBINS))
    end if
    if (tot.gt.0.0) then
      effmoltyp = effmoltyp + 1 
      rdhist(k,1:TPBINS,1:DLSTBINS) = rdhist(k,1:TPBINS,1:DLSTBINS)/tot
    end if
  end do
  if (allocated(densproavg).EQV..true.) then
    do i=1,maxrgbins
      volis(i) = ((i*rg_binsz)**3 - ((i-1)*rg_binsz)**3)*4.*PI/3.
    end do
    do k=1,nangrps
      if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
      if (normer(k).le.0.0) cycle
      densproavg(k,:) = convdens*densproavg(k,:)/(normer(k)*volis(:))
    end do
  end if
  deallocate(volis)
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_POLYAVG.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'POLYAVG.dat'
  end if
#else
  fn = 'POLYAVG.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then   
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    do k=1,nangrps
!
      if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
!
      if (use_frame_weights.EQV..true.) then
        rgavg(k) = rgavg(k)*npolavg(k)/normer(k)
        rg2avg(k) = rg2avg(k)*npolavg(k)/normer(k)
        vtavg(k) = vtavg(k)*npolavg(k)/normer(k)
        asphavg(k) = asphavg(k)*npolavg(k)/normer(k)
        acylavg(k) = acylavg(k)*npolavg(k)/normer(k)
        rgterm12avg(k) = rgterm12avg(k)*npolavg(k)/normer(k)
        rgterm1avg(k) = rgterm1avg(k)*npolavg(k)/normer(k)
        rgterm2avg(k) = rgterm2avg(k)*npolavg(k)/normer(k)
        dlstavg(k) = 1.0 - 3.0*rgterm12avg(k)
        dlavg(k) = 1.0 - 3.0*(rgterm1avg(k)/rgterm2avg(k))
        reteavg(k) = reteavg(k)*npolavg(k)/normer(k)
        do ii=1,3
          rgpcavg(k,ii) = rgpcavg(k,ii)*npolavg(k)/normer(k)
        end do
      end if 
!
      exists = .false.
      if (allocated(nrhavg).EQV..true.) then
        if (nrhavg(k).gt.0) then
          rhravg(k) = normer2(k)/rhravg(k)
          exists = .true.
        end if
      end if
      write(iu,30) k,molangr(k,1),rsmol(molangr(k,1),1),rsmol(molangr(k,1),2)
      if (exists.EQV..true.) then
        write(iu,28) rgavg(k),rg2avg(k),vtavg(k),asphavg(k)/rg2avg(k),acylavg(k)/rg2avg(k),&
 &                   dlstavg(k),dlavg(k),reteavg(k),rhravg(k),rgpcavg(k,3),rgpcavg(k,2),rgpcavg(k,1)
      else
        write(iu,28) rgavg(k),rg2avg(k),vtavg(k),asphavg(k)/rg2avg(k),acylavg(k)/rg2avg(k),&
 &                   dlstavg(k),dlavg(k),reteavg(k),0.0,rgpcavg(k,3),rgpcavg(k,2),rgpcavg(k,1)
      end if
    end do
    close(unit=iu)
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_RDHIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'RDHIST.dat'
  end if
#else
  fn = 'RDHIST.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then   
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    do k=1,nangrps
      if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
      write(iu,30) k,molangr(k,1),rsmol(molangr(k,1),1),rsmol(molangr(k,1),2)
      do j=1,DLSTBINS
        write(iu,25) (rdhist(k,i,j),i=1,TPBINS)
      end do
    end do
    close(unit=iu)
  end if
!
! now Ree-histograms
  maxi = 0
  effmoltyp = 0
  do k=1,nangrps
    tot = 0.0
    if (allocated(retehist).EQV..true.) then
      do i=1,maxrgbins
        if ((retehist(k,i).gt.0.0).AND.(i.gt.maxi)) then
          maxi = i
        end if
      end do
      tot = sum(retehist(k,:))
    end if
    if (tot.gt.0.0) then
      effmoltyp = effmoltyp + 1
      polynew%map(effmoltyp) = k
      retehist(k,:) = retehist(k,:)/tot
    end if
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_RETEHIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'RETEHIST.dat'
  end if
#else
  fn = 'RETEHIST.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    write(iu,29) (polynew%map(k),k=1,effmoltyp)
    do i=1,maxi
      rval = (i-0.5)*rg_binsz
      write(iu,27) rval,(retehist(polynew%map(k),i),k=1,effmoltyp)
    end do
    close(unit=iu)
  end if
!
! and Rg-histograms
  maxi = 0
  effmoltyp = 0
  do k=1,nangrps
    tot = 0.0
    if (allocated(rghist).EQV..true.) then
      do i=1,maxrgbins
        if ((rghist(k,i).gt.0.0).AND.(i.gt.maxi)) then
          maxi = i
        end if
      end do
      tot = sum(rghist(k,:))
    end if
    if (tot.gt.0.0) then
      effmoltyp = effmoltyp + 1
      polynew%map(effmoltyp) = k
      rghist(k,:) = rghist(k,:)/tot
    end if
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_RGHIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'RGHIST.dat'
  end if
#else
  fn = 'RGHIST.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    write(iu,29) (polynew%map(k),k=1,effmoltyp)
    do i=1,maxi
      rval = (i-0.5)*rg_binsz
      write(iu,27) rval,(rghist(polynew%map(k),i),k=1,effmoltyp)
    end do
    close(unit=iu)
  end if
!
! and finally density profiles
  maxi = 0
  effmoltyp = 0
  do k=1,nangrps
    if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
    if (allocated(densproavg).EQV..true.) then
      do i=1,maxrgbins
        if ((densproavg(k,i).gt.0.0).AND.(i.gt.maxi)) then
          maxi = i
        end if
      end do
    end if
    if (maxi.gt.0) then
      effmoltyp = effmoltyp + 1
      polynew%map(effmoltyp) = k
    end if
  end do
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_DENSPROF.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'DENSPROF.dat'
  end if
#else
  fn = 'DENSPROF.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then   
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    write(iu,29) (polynew%map(k),k=1,effmoltyp)
    do i=1,maxi
      rval = (i-0.5)*rg_binsz
      write(iu,27) rval,(densproavg(polynew%map(k),i),k=1,effmoltyp)
    end do
    close(unit=iu)
  end if
!
  deallocate(normer2)
  deallocate(normer)
!
end
!
!----------------------------------------------------------------------
!
! #######################################################
! ##                                                   ##
! ## subroutine rtube -- thickness of a consensus tube ##
! ##                                                   ##
! #######################################################
!
!
! "rtube" computes the thickness of a Banavar (generalized Edwards) 
! tube to estimate the wiggle room around a specific conformation
! provides an estimate of how the local radius of curvature changes
! in the ensemble; outputs the largest and smallest triplet radii
!
!
subroutine rtube(rthick)
!
  use iounit
  use polypep
  use sequen
!
  implicit none
!
  integer i,ires,j,jres,k,kres,ica,jca,kca
  integer is,js,ks,nloop
  RTYPE rr,rthick
!
! set logical variable to be true
  rthick = 1.0d080
  nloop = nseq - 2
  do i = 1,nloop-2
     ires = i + 1
     ica = cai(ires)
     do j = i+1,nloop-1
        jres = j + 1
        jca = cai(jres)
        do k = j+1,nloop
           kres = k + 1
           kca = cai(kres)
           call triplet (ica,jca,kca,rr)
           if(rr .lt. rthick) then 
              rthick = rr
              is = i
              js = j
              ks = k
           end if
        end do
     end do
  end do
  return
end
!
!
!
! #################################################################
! ##                                                             ##
! ## subroutine triplet - computes the radius of a circumcircle  ##
! ##                      that passes through three chosen atoms ##
! ##                                                             ##
! #################################################################
!
!
! "triplet" calculates the radius of the circle that circumscribes
! a triplet of atoms; the chosen atom numbers are input to the routine
! execution is terminated if any two atom numbers are identical
!     
subroutine triplet (i,j,k,rad)
!
  use atoms
!
  implicit none
!
  integer i,j,k
  RTYPE a,b,c,denmtr,rad,s
  RTYPE xij,yij,zij,xik,yik,zik,xjk,yjk,zjk
!
! Literature reference
! See http://www.math.niu.edu/~rusin/known-math/99/area_circle
!
!
  xij = x(i) - x(j)
  yij = y(i) - y(j)
  zij = z(i) - z(j)
  a = sqrt(xij**2 + yij**2 + zij**2)
!     
  xik = x(i) - x(k)
  yik = y(i) - y(k)
  zik = z(i) - z(k)
  b = sqrt(xik**2 + yik**2 + zik**2)
!     
  xjk = x(j) - x(k)
  yjk = y(j) - y(k)
  zjk = z(j) - z(k)
  c = sqrt(xjk**2 + yjk**2 + zjk**2)
!     
  s = 0.5d0*(a + b + c)
  denmtr = sqrt(s*(s-a)*(s-b)*(s-c))
  if(denmtr .gt. 0.0d0) then
     rad = 0.25d0*a*b*c/denmtr
  else
     rad = huge(rad)
  end if
!
  return
!
end
!
!----------------------------------------------------------------------------
!
! system-wide radius of gyration tensor
! note that since this is intermolecular, boundary conditions have to be accounted
! for (whereas this is non-sensical for all purely intramolecular measures)
!
subroutine gyrate2(rgs,tps,rgsten,ev,asphs,acyls,dlts,dltsts)
!
  use iounit
  use atoms
  use molecule
  use polyavg
  use sequen
!
  implicit none
!
  integer i,j,imol
  RTYPE rgs,rgsten(3,3),asphs,acyls,dlts,dltsts,t1,t2,t12,tps
  RTYPE xc,yc,zc,evmat(3,3),ev(3),cln
  RTYPE shcom(nmol,3),sh(3)
!
  do i=1,3
    do j=1,3
      rgsten(i,j) = 0.0
    end do
  end do
  xc = 0.0
  yc = 0.0
  zc = 0.0
!
  do imol=1,nmol
    call update_rigid(imol)
    call shift_bound(1,imol,sh)
    do j=1,3
      shcom(imol,j) = sh(j)
    end do
    xc = xc + (atmol(imol,2)-atmol(imol,1)+1)*shcom(imol,1)
    yc = yc + (atmol(imol,2)-atmol(imol,1)+1)*shcom(imol,2)
    zc = zc + (atmol(imol,2)-atmol(imol,1)+1)*shcom(imol,3)
  end do
!
  xc = xc/dble(n)
  yc = yc/dble(n)
  zc = zc/dble(n)
!
  do imol=1,nmol
    do i=1,3
      sh(i) = shcom(imol,i) - com(imol,i)
    end do
    do i=atmol(imol,1),atmol(imol,2)
      rgsten(1,1) = rgsten(1,1) + (x(i)+sh(1)-xc)*(x(i)+sh(1)-xc)
      rgsten(1,2) = rgsten(1,2) + (x(i)+sh(1)-xc)*(y(i)+sh(2)-yc)
      rgsten(1,3) = rgsten(1,3) + (x(i)+sh(1)-xc)*(z(i)+sh(3)-zc)
      rgsten(2,2) = rgsten(2,2) + (y(i)+sh(2)-yc)*(y(i)+sh(2)-yc)
      rgsten(2,3) = rgsten(2,3) + (y(i)+sh(2)-yc)*(z(i)+sh(3)-zc)
      rgsten(3,3) = rgsten(3,3) + (z(i)+sh(3)-zc)*(z(i)+sh(3)-zc)
    end do
  end do
!
  rgsten(2,1) = rgsten(1,2)
  rgsten(3,1) = rgsten(1,3)
  rgsten(3,2) = rgsten(2,3)
  do i=1,3
    do j=1,3
      rgsten(i,j) = rgsten(i,j)/(dble(n))
    end do
  end do
!
  rgs = sqrt(rgsten(1,1)+rgsten(2,2)+rgsten(3,3))
  cln = 0.0
  do i=1,nmoltyp
    cln = cln + moltyp(i,2)*molcontlen(i)
  end do
  if (cln.gt.0.0) then 
    tps = 2.5*(1.75*rgs/cln)**(4.0/(dble(nseq))**tp_exp)
  else
    tps = 0.0
  end if
  call mat_diag(3,rgsten,ev,evmat)
  t1 = ev(1)*ev(2) + ev(2)*ev(3) + ev(1)*ev(3)
  t2 = (ev(1) + ev(2) + ev(3))**2
  t12 = t1/t2
  dltsts = 1.0 - 3.0*t12
  dlts = 1.0 - 3.0*(t1/t2)
  asphs = ev(3) - 0.5*(ev(2) + ev(1))
  acyls = ev(2) - ev(1)
!
end
!
!----------------------------------------------------------------------------
!
! "hydroradius" computes internal scaling behaviour through
! averaging of all unique internal distances r_ij as a function
! of sequence distance.
! the contribution of adding the self-distance (r_ii) is trivial
! to add analytically (this was claimed by people in the lab to be
! relevant / conventional)
!
subroutine hydroradius(rhr,rh,psc,dosct)
!
  use atoms
  use molecule
  use sequen
  use polypep
  use polyavg
  use system
  use grandensembles
!
  implicit none
!
  integer i,j,k,imol,rs,rs2,ii,kk,seqsep,mty,qv,nmt,incr
  RTYPE rh(nangrps,longestmol),xd,yd,zd,qrr,rvv
  RTYPE psc(nangrps,nqv),rhr(nangrps)
  logical dosct,ismember,refmoled(nangrps),foundit
!
  do i=1,nangrps
    polynew%nmr(i) = 0
    rhr(i) = 0.0
    do j=1,rsmol(molangr(i,1),2)-rsmol(molangr(i,1),1)+1
      polynew%nm(i,j) = 0 
      rh(i,j) = 0.0d0
    end do
    if (dosct.EQV..true.) then
      do j=1,nqv
        psc(i,j) = 0.0
      end do
    end if
  end do
!
  refmoled(:) = .false.
!
  do imol=1,nmol
!
    foundit = .true.
!
    if (do_pol(moltypid(imol)).EQV..false.) cycle
    mty = an_grp_mol(imol)
!   cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
!
    do rs=rsmol(imol,1),rsmol(imol,2)
      do rs2=rs,rsmol(imol,2)
!       we offset this by 1 to avoid writing to the zero-element
        seqsep = rs2-rs+1
!
!       get the normalization number (only for the first successful occurence)
        if (refmoled(mty).EQV..false.) then
          if (seqsep.eq.1) then
            incr = ((at(rs)%nbb+at(rs)%nsc)*(at(rs2)%nbb+at(rs2)%nsc-1))/2
            polynew%nm(mty,seqsep) = polynew%nm(mty,seqsep) + incr
            polynew%nmr(mty) = polynew%nmr(mty) + incr
          else
            incr = (at(rs)%nbb+at(rs)%nsc)*(at(rs2)%nbb+at(rs2)%nsc)
            polynew%nm(mty,seqsep) = polynew%nm(mty,seqsep) + incr
            polynew%nmr(mty) = polynew%nmr(mty) + incr
          end if
          if (foundit.EQV..false.) foundit = .true.
        end if
!
        do i=1,at(rs)%nbb+at(rs)%nsc
          if (i.le.at(rs)%nbb) then
            ii = at(rs)%bb(i)
          else
            ii = at(rs)%sc(i-at(rs)%nbb)
          end if
          if (seqsep.eq.1) then
            do k=i+1,at(rs2)%nbb+at(rs2)%nsc
              if (k.le.at(rs2)%nbb) then
                kk = at(rs2)%bb(k)
              else
                kk = at(rs2)%sc(k-at(rs2)%nbb)
              end if
              xd = x(ii) - x(kk)
              yd = y(ii) - y(kk)
              zd = z(ii) - z(kk)
              rvv = sqrt(xd*xd + yd*yd + zd*zd)
              rh(mty,seqsep) = rh(mty,seqsep) + rvv
              rhr(mty) = rhr(mty) + 1.0/rvv
              if (dosct.EQV..true.) then
                do qv=1,nqv
                  qrr = ((qv-0.5)*qv_res) *rvv
                  psc(mty,qv) = psc(mty,qv) + sin(qrr)/qrr
                end do
              end if
            end do
          else
            do k=1,at(rs2)%nbb+at(rs2)%nsc
              if (k.le.at(rs2)%nbb) then
                kk = at(rs2)%bb(k)
              else
                kk = at(rs2)%sc(k-at(rs2)%nbb)
              end if
              xd = x(ii) - x(kk)
              yd = y(ii) - y(kk)
              zd = z(ii) - z(kk)
              rvv = sqrt(xd*xd + yd*yd + zd*zd)
              rh(mty,seqsep) = rh(mty,seqsep) + rvv
              rhr(mty) = rhr(mty) + 1./rvv
              if (dosct.EQV..true.) then
                do qv=1,nqv
                  qrr = ((qv-0.5)*qv_res) *rvv
                  psc(mty,qv) = psc(mty,qv) + sin(qrr)/qrr
                end do
              end if
            end do
          end if
        end do
      end do
    end do
    if (foundit.EQV..true.) refmoled(mty) = .true.
  end do
!
  do i=1,nangrps
    if (polynew%nmr(i).gt.0) then
      rhr(i) = rhr(i)/dble(polynew%nmr(i))
    end if
    nmt = 0
    do j=1,rsmol(molangr(i,1),2)-rsmol(molangr(i,1),1)+1
      if (polynew%nm(i,j).gt.0) then
        rh(i,j) = rh(i,j)/dble(polynew%nm(i,j))
        nmt = nmt + polynew%nm(i,j)
      end if
    end do
    if (dosct.EQV..true.) then
      if (nmt.gt.0) then
        do j=1,nqv
          psc(i,j) = psc(i,j)/(1.0*nmt)
        end do
      end if
    end if
  end do
! 
end
!
!-----------------------------------------------------------------------------
!
subroutine prt_hydrorad()
!
  use iounit
  use molecule
  use polyavg
  use mpistuff
  use pdb
  use system
!
  implicit none
!
  integer i,k,ii,jj,iu,freeunit,effmoltyp
  character(100) fn
  RTYPE qvv
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  RTYPE, ALLOCATABLE:: normer(:)
  logical exists
!
 24   format('# Molecule type ',i3,' (Ref.: Mol. ',i5,' / Res. ',i5,'-',&
 &i5,'):')
 27   format(i6,1x,g14.7)
 28   format(f12.6,1x,g14.7,1x,g14.7)
!
  effmoltyp = 0
  do k=1,nangrps
    if (nrhavg(k).gt.0) effmoltyp = effmoltyp + 1
  end do
!
  allocate(normer(nangrps))
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),rhcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)*molangr(:,2)
      end if
    end do
    if (k.ne.nrhavg(1)/molangr(1,2)) then
      write(ilog,*) 'Warning. Pairwise, polymeric analyses have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in INTSCAL.dat and related files.'
    end if
  else
    normer(:) = 1.0*nrhavg(:)
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_INTSCAL.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'INTSCAL.dat'
  end if
#else
  fn = 'INTSCAL.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    do k=1,nangrps
      if (nrhavg(k).gt.0) then
        do i=1,rsmol(molangr(k,1),2)-rsmol(molangr(k,1),1)+1
          rhavg(k,i) = rhavg(k,i)/normer(k)
        end do
      end if
    end do
    do k=1,nangrps
      if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
      write(iu,24) k,molangr(k,1),rsmol(molangr(k,1),1),rsmol(molangr(k,1),2)
      do i=1,rsmol(molangr(k,1),2)-rsmol(molangr(k,1),1)+1
        write(iu,27) i-1,rhavg(k,i)
      end do
    end do
    close(unit=iu)
  end if
!
  effmoltyp = 0
  do k=1,nangrps
    if (nsctavg(k).gt.0) effmoltyp = effmoltyp + 1
  end do
!
  if (use_frame_weights.EQV..true.) then
    normer(:) = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),rhcalc).eq.0).AND.(mod(framelst2(i)/rhcalc,sctcalc).eq.0)) then
        k = k + 1
        normer(:) = normer(:) + framewts2(i)*molangr(:,2)
      end if
    end do
    if (k.ne.nsctavg(1)/molangr(1,2)) then
      write(ilog,*) 'Warning. Theoretical scattering analyses have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in KRATKY.dat and related files.'
    end if
  else
    normer(:) = 1.0*nsctavg(:)
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_KRATKY.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'KRATKY.dat'
  end if
#else
  fn = 'KRATKY.dat'
#endif
  call strlims(fn,ii,jj)
!
  if (effmoltyp.gt.0) then
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      iu = freeunit()
      open(unit=iu,file=fn(ii:jj),status='old')
      close(unit=iu,status='delete')
    end if
    iu=freeunit()
    open(unit=iu,file=fn(ii:jj),status='new')
!
    do k=1,nangrps
      if (nsctavg(k).gt.0) then
        do i=1,nqv
          pscavg(k,i) = pscavg(k,i)/normer(k)
        end do
      end if
    end do
    do k=1,nangrps
      if (do_pol(moltypid(molangr(k,1))).EQV..false.) cycle
      write(iu,24) k,molangr(k,1),rsmol(molangr(k,1),1),rsmol(molangr(k,1),2)
      do i=1,nqv
        qvv = (i-0.5)*qv_res
        write(iu,28) qvv,pscavg(k,i),qvv*qvv*pscavg(k,i)
      end do
    end do
    close(unit=iu)
  end if
!
  deallocate(normer)
!
end
!
!------------------------------------------------------------------------------
!
subroutine do_diffraction()
!
  use atoms
  use math
  use units
  use diffrac
  use iounit
  use polypep
  use molecule
  use sequen
  use system
  use grandensembles
  use pdb
!
  implicit none
!
  integer i,j,rr,zz,cyli,cylj,nscat,imol
  RTYPE rmax,rt,cylv(3),cylvn,v1(3),v2(3),rfv(3)
  RTYPE v1n,sine,cosine,rfr,rfz,rfvn,cosine2,rval,zval
  RTYPE fac,term0,term1,term2,bessel,eps,axon(3),incr
  RTYPE, ALLOCATABLE:: cyz(:),cyr(:),cyt(:),cym(:)
  logical ismember
!
#ifdef DISABLE_FLOAT
  eps = 0.0000001d0
#else
! the low-precision xyz-coordinates lead to some hysteresis error in converting back and forth
  eps = 0.02d0
#endif
!
   23 format(4(g14.6,1x))
!
  allocate(cyz(n))
  allocate(cyt(n))
  allocate(cyr(n))
  allocate(cym(n))
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if
!
! first convert the system to cylindrical coordinates
!
  if (diffr_uax.EQV..true.) then
!   we'll anchor the axis to the point provided by the user
    do j=1,3
      cylv(j) = diffr_axis(j)
      axon(j) = diffr_axon(j)
    end do
    cyli = 0
    cylj = 0
  else
!   find the axis unit vector and a normal which defines the zero-point for the theta-angle
    rmax = 0.0
    do i=1,n
      do j=i+1,n
        rt = (x(i)-x(j))*(x(i)-x(j)) + &
 &           (y(i)-y(j))*(y(i)-y(j)) +&
 &           (z(i)-z(j))*(z(i)-z(j))
        if (rt.gt.rmax) then
          rmax = rt
          cylv(1) = x(i)-x(j)
          cylv(2) = y(i)-y(j)
          cylv(3) = z(i)-z(j)
          cyli = i
          cylj = j
        end if
      end do
    end do
    axon(1) = x(cylj)
    axon(2) = y(cylj)
    axon(3) = z(cylj)
  end if
!  cylj = 2473
!  cyli = 503
!       cylv(1) = x(cyli)-x(cylj)
!        cylv(2) = y(cyli)-y(cylj)
!        cylv(3) = z(cyli)-z(cylj)
!  write(*,*) cyli,cylj
!  write(*,*) cylv(1),cylv(2),cylv(3)
  cylvn = sqrt(cylv(1)*cylv(1) + cylv(2)*cylv(2) + cylv(3)*cylv(3))
  do i=1,3
    cylv(i) = cylv(i)/cylvn
  end do
  v1(1) = -axon(1)
  v1(2) = -axon(2)
  v1(3) = -axon(3)
  v1n = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3))
  do j=1,3
    v1(j) = v1(j)/v1n
  end do
  call crossprod(cylv,v1,v2,sine)
  rfr = v1n*sine
  cosine = v1(1)*cylv(1) + v1(2)*cylv(2) + v1(3)*cylv(3)
  rfz = v1n*cosine
  rfv(1) = -(axon(1) + rfz*cylv(1))
  rfv(2) = -(axon(2) + rfz*cylv(2))
  rfv(3) = -(axon(3) + rfz*cylv(3))
  rfvn = sqrt(rfv(1)*rfv(1) + rfv(2)*rfv(2) + rfv(3)*rfv(3))
  do j=1,3
    rfv(j) = rfv(j)/rfvn
  end do
!  write(*,*) cylv(1)*rfv(1)+cylv(2)*rfv(2)+cylv(3)*rfv(3)
!
! now generate cylindrical coordinates for each (present) atom
  nscat = 0
  do imol=1,nmol
!   cycle out for non-present molecules
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    do i=atmol(imol,1),atmol(imol,2)
      if (mass(i).le.4.0) cycle
      nscat = nscat + 1
!     mass
      cym(nscat) = mass(i)
!     first get the distance from the axis
      if ((i.eq.cyli).OR.(i.eq.cylj)) then
        cyr(nscat) = 0.0
      else
        v1(1) = x(i)-axon(1)
        v1(2) = y(i)-axon(2)
        v1(3) = z(i)-axon(3)
        v1n = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3))
        do j=1,3
          v1(j) = v1(j)/v1n
        end do
        call crossprod(cylv,v1,v2,sine)
        cyr(nscat) = v1n*sine
      end if
!     now get the z-coordinate along the axis
      if (i.eq.cylj) then
        cyz(nscat) = 0.0
      else if (i.eq.cyli) then
        cyz(nscat) = cylvn
      else
        cosine = v1(1)*cylv(1) + v1(2)*cylv(2) + v1(3)*cylv(3)
        cyz(nscat) = v1n*cosine
      end if
!     and finally the angular coordinate
      if ((i.eq.cyli).OR.(i.eq.cylj)) then
        cyt(nscat) = PI
      else
        v1(1) = x(i) - (cyz(nscat)*cylv(1) + axon(1))
        v1(2) = y(i) - (cyz(nscat)*cylv(2) + axon(2))
        v1(3) = z(i) - (cyz(nscat)*cylv(3) + axon(3))
        v1n = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3))
        do j=1,3
          v1(j) = v1(j)/v1n
        end do
        call crossprod(rfv,v1,v2,sine)
        cosine =  v1(1)*rfv(1) + v1(2)*rfv(2) + v1(3)*rfv(3)
        cosine2 = v2(1)*cylv(1) + v2(2)*cylv(2) + v2(3)*cylv(3)
        if ((cosine.gt.(1.0+eps)).OR.(cosine.lt.-(1.0+eps))) then
          write(ilog,*) 'Fatal. This is a bug in the computation of cy&
 &lindrical coordinates for the purpose of diffraction analysis. Ple&
 &ase report this problem.'
          call fexit()
        else if (cosine.ge.0.0) then
          cosine = min(cosine,1.0d0)
          if (cosine2.ge.0.0) then
            cyt(nscat) = acos(cosine)
          else
            cyt(nscat) = 2.0*PI - acos(cosine)
          end if
        else
          cosine = max(cosine,-1.0d0)
          if (cosine2.ge.0.0) then
            cyt(nscat) = acos(cosine)
          else
            cyt(nscat) = 2.0*PI - acos(cosine)
          end if
        end if
      end if
    end do
  end do
!
! now apply the Bessel-Fourier transform over the 2D-map
! note that this can be an extremely expensive operation depending
! on the dimensions of the system and the resolution of the map
  do rr=1,rr_max
    rval = (rr-0.5)*rr_res - rr_res*(rr_max/2)
    do zz=1,zz_max
      zval = (zz-0.5)*zz_res - zz_res*(zz_max/2)
      do j=0,bes_max
        term1 = 0.0
        term2 = 0.0 
        do i=1,nscat
!         very crude mass proportionality (which should really be
!         the proper angle-dependent integral over the actual atomic
!         electron density -> ha!)
          fac = cym(i)/12.0108
          term0 = fac*bessel(j,rval*cyr(i))
          term1 = term1 + term0*cos(zval*cyz(i) - j*cyt(i))
          term2 = term2 + term0*sin(zval*cyz(i) - j*cyt(i))
        end do
        am_diffr_map(rr,zz) = am_diffr_map(rr,zz) + &
 &                                 incr*(term1*term1 + term2*term2)
      end do
    end do
  end do
!
  diffr_cnt = diffr_cnt + 1
!  
  deallocate(cyr)
  deallocate(cyt)
  deallocate(cyz)
  deallocate(cym)
!
end
!
!---------------------------------------------------------------------
!
subroutine prt_diffraction()
!
  use iounit
  use molecule
  use math
  use diffrac
  use mpistuff
  use units
  use system
  use pdb
!
  implicit none
!
  integer i,j,ii,jj,iu,freeunit,k
  character(100) fn
  logical exists
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  RTYPE normer
!
 25   format(10000000(g14.7,1x))
!
  if (use_frame_weights.EQV..true.) then
    normer = 0.0
    k = 0
    do i=1,framecnt
      if ((framelst2(i).gt.nequil).AND.(mod(framelst2(i),diffrcalc).eq.0)) then
        k = k + 1
        normer = normer + framewts2(i)
      end if
    end do
    if (k.ne.diffr_cnt) then
      write(ilog,*) 'Warning. Average diffraction maps have inconsistent normalization. This is &
 &a bug. Please report this error and ignore all output in DIFFRACTION.dat and related files.'
    end if
  else
    normer = 1.0*diffr_cnt
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_DIFFRACTION.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'DIFFRACTION.dat'
  end if
#else
  fn = 'DIFFRACTION.dat'
#endif
!  
  call strlims(fn,ii,jj) 
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
!
  do i=1,rr_max
    do j=1,zz_max
      am_diffr_map(i,j) = am_diffr_map(i,j)/normer
    end do
  end do
!
  do j=1,zz_max
    write(iu,25) (am_diffr_map(i,j),i=1,rr_max)
  end do
  close(unit=iu)
!
end
!
!-----------------------------------------------------------------------
!
subroutine bessel_read()
!
  use iounit
  use system
  use diffrac
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer i,j,freeunit,it,t1,t2,iomess
  logical exists
!
  call strlims(besselfile,t1,t2)
  inquire(file=besselfile(t1:t2),exist=exists)
  if(exists) then
    it = freeunit()
    open(unit=it,file=besselfile(t1:t2),status='old')
    read(it,*,iostat=iomess) bes_num,bes_bins,bes_res
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Warning. Got empty/incomplete file for Bessel f&
 &unctions (',besselfile(t1:t2),'). Turning off diffraction calc..'
      diffrcalc = nsim + 1
    end if
    allocate(bessels(bes_num,bes_bins))
    do i=1,bes_bins
      read(it,*,iostat=iomess) (bessels(j,i),j=1,bes_num)
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*) 'Warning. Got incomplete file for Bessel funct&
 &ions (',besselfile(t1:t2),'). Turning off diffraction calc..'
        diffrcalc = nsim + 1
        deallocate(bessels)
      end if
    end do
    close(unit=it)
  else
    write(ilog,*) 'File with tabulated Bessel Functions missing (',&
 &besselfile(t1:t2),'). Turning off diffraction calc..'
    diffrcalc = nsim + 1
  end if
!
! adjust max order if necessary
  if (bes_max.gt.bes_num) then
    write(ilog,*) 'Requested higher order of Bessel functions (',&
 &bes_max,') than are provided in file (',besselfile(t1:t2),'). Adju&
 &sting to maximum possible order (',bes_num,').'
    bes_max = bes_num
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine bessel_close()
!
  use iounit
  use diffrac
!
  implicit none
!
  if (allocated(bessels).EQV..false.) then
    write(ilog,*) 'Fatal. This is a bug. Bessel function array was n&
 &ot properly allocated, or prematurely freed. Please fix or report &
 &this problem.'
    call fexit()
  else
    deallocate(bessels)
  end if
!
end
!
!-----------------------------------------------------------------------
!
! this returns the real part for the Bessel fxn of order 'order'
! for the argument 'val'
! note that for negative arguments the solution is formally complex,
! but the imaginary components are always zero(?)
! Bessel fxn's of even order are axis-symmetric while those of odd
! order are point-symmetric 
!
function bessel(order,val)
!
  use diffrac
!
  implicit none
!
  integer order,bin
  RTYPE val,bessel,fr1
!
  bin = floor(abs(val)/bes_res)
!
  if (bin.lt.1) then
    bessel = bessels(order+1,1)
  else if (bin.ge.bes_bins) then
    bessel = bessels(order+1,bes_bins)
  else
    fr1 = (abs(val)/bes_res) - 1.0*bin
    bessel = fr1*bessels(order+1,bin) + &
 &                (1.0-fr1)*bessels(order+1,bin+1)
  end if
!
  if (val.lt.0.0) then
    if ((order.ne.0).AND.(mod(order,2).ne.0)) then
      bessel = -bessel
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
