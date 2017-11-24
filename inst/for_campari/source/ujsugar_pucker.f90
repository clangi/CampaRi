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
! MAIN AUTHOR:   Adam Steffen                                              !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
subroutine freepucker(rs,d_pucker,mode)
! ---------------------------------------------------------------------
!  Subroutine that varies the 4 free dofs of a 5 member ring
!
!  INPUT: rs - index of ring
!  OUTPUT: d_freepucker - set of change in DOFs 
! ---------------------------------------------------------------------
!
  use iounit
  use polypep
  use zmatrix
  use atoms
  use math
  use aminos
  use sequen
  use movesets
  use molecule
  use system
  use fyoc
!
  implicit none
!
  integer i,shf,idx3
  RTYPE random
  integer, INTENT (IN) :: rs, mode
  RTYPE, INTENT (OUT) :: d_pucker(7)
!
! has no memory of previous position.  Is this better than applying a deviation?
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! mode 1) reflect phi and chi1 2) random purturbation move
  if (mode.eq.1) then
!   cyclic peptide residues
    if (seqflag(rs).eq.5) then
      d_pucker(1) = -ztor(at(rs)%sc(3-shf))
      d_pucker(2) = -ztor(at(rs)%sc(4-shf))
      d_pucker(5) = bang(at(rs)%sc(2-shf)) 
      d_pucker(6) = bang(at(rs)%sc(3-shf))
!   nucleotides
    else if (seqpolty(rs).eq.'N') then
      d_pucker(1) = -ztor(at(rs)%sc(2))
      d_pucker(2) = -ztor(at(rs)%sc(3))
      d_pucker(5) = bang(at(rs)%sc(1))
      d_pucker(6) = bang(at(rs)%sc(2))
    end if
    d_pucker(7) = 0.0
  else if (mode.eq.2) then
!   cyclic peptide residues
    if (seqflag(rs).eq.5) then
      idx3 = pucline(rs)
      d_pucker(1) = ztor(at(rs)%sc(3-shf)) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(2) = ztor(at(rs)%sc(4-shf)) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(3) = ztor(idx3) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(4) = ztor(at(rs)%sc(2-shf)) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(5) = bang(at(rs)%sc(2-shf)) + random()*pucker_anstp - 0.5*pucker_anstp
      d_pucker(6) = bang(at(rs)%sc(3-shf)) + random()*pucker_anstp - 0.5*pucker_anstp
      d_pucker(7) = 0.0
!   nucleotides
    else if (seqpolty(rs).eq.'N') then
      if (rs.eq.rsmol(molofrs(rs),2)) then
        if (moltermid(molofrs(rs),2).ne.1) then
          write(ilog,*) 'Fatal. Encountered unknown terminus t&
 &ype(s) for molecule ',molofrs(rs),' in sample_sugar(...).'
          call fexit()
        end if
      end if
      idx3 = nucsline(6,rs)
      d_pucker(1) = ztor(at(rs)%sc(2)) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(2) = ztor(at(rs)%sc(3)) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(3) = ztor(idx3) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(4) = ztor(at(rs)%sc(1)) + random()*pucker_distp - 0.5*pucker_distp
      d_pucker(5) = bang(at(rs)%sc(1)) + random()*pucker_anstp - 0.5*pucker_anstp
      d_pucker(6) = bang(at(rs)%sc(2)) + random()*pucker_anstp - 0.5*pucker_anstp
      d_pucker(7) = 0.0
    else
      write(ilog,*) 'Fatal. Called freepucker(...) with unsupported residue (',rs,' of&
 & type ',seqtyp(rs),': ',amino(seqtyp(rs)),').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Called freepucker method with unsupported mode (',mode,').'
    call fexit()
  end if
!
! get differences
! cyclic peptide residues
  if (mode.eq.2) then
!   wrap-arounds are dealt with in setconf.f90- here check only for bond angle exceptions
    do i=5,6
      if ((d_pucker(i).gt.180.0).OR.(d_pucker(i).lt.0.0)) then
        write(ilog,*) 'Fatal. Bond angle out of bounds in freepucker(...).'
        call fexit()
      end if
    end do
    if (seqflag(rs).eq.5) then
      d_pucker(1) = d_pucker(1) - ztor(at(rs)%sc(3-shf))
      d_pucker(2) = d_pucker(2) - ztor(at(rs)%sc(4-shf))
      d_pucker(3) = d_pucker(3) - ztor(idx3)
      d_pucker(4) = d_pucker(4) - ztor(at(rs)%sc(2-shf))
      d_pucker(5) = d_pucker(5) - bang(at(rs)%sc(2-shf))
      d_pucker(6) = d_pucker(6) - bang(at(rs)%sc(3-shf))
!   nucleotides
    else if (seqpolty(rs).eq.'N') then
      d_pucker(1) = d_pucker(1) - ztor(at(rs)%sc(2))
      d_pucker(2) = d_pucker(2) - ztor(at(rs)%sc(3))
      d_pucker(3) = d_pucker(3) - ztor(idx3)
      d_pucker(4) = d_pucker(4) - ztor(at(rs)%sc(1))
      d_pucker(5) = d_pucker(5) - bang(at(rs)%sc(1))
      d_pucker(6) = d_pucker(6) - bang(at(rs)%sc(2)) 
    end if
  else
    d_pucker(3:4) = 0.0
    if (seqflag(rs).eq.5) then
      d_pucker(1) = d_pucker(1) - ztor(at(rs)%sc(3-shf))
      d_pucker(2) = d_pucker(2) - ztor(at(rs)%sc(4-shf))
      d_pucker(5) = d_pucker(5) - bang(at(rs)%sc(2-shf))
      d_pucker(6) = d_pucker(6) - bang(at(rs)%sc(3-shf))
!   nucleotides
    else if (seqpolty(rs).eq.'N') then
      d_pucker(1) = d_pucker(1) - ztor(at(rs)%sc(2))
      d_pucker(2) = d_pucker(2) - ztor(at(rs)%sc(3))
      d_pucker(5) = d_pucker(5) - bang(at(rs)%sc(1))
      d_pucker(6) = d_pucker(6) - bang(at(rs)%sc(2)) 
    end if
  end if
!
end subroutine freepucker
!
! !-----------------------------------------------------------------------
! ! Subroutine that finds the jacobian factor of the 5x5 matrix needed for
! ! chain closure.
! !
! ! INPUT:  rs - index of ring
! ! OUTPUT: jac_x - that jacobian value
! ! ----------------------------------------------------------------------
subroutine jacobian_pucker(rs,jac_x)
!
  use iounit
  use polypep
  use atoms
  use math
  use sequen
  use molecule
  use aminos
  use system
!
  implicit none
!
  integer k,indx(5),shf,shfn
  RTYPE c1,c2,c3,c4,u1(3),u2(3),u3(3),u4(3),s1(3),s2(3),s3(3)
  RTYPE s4(3),s5(3),w1(3),w2(3),w3(3),t0(3),cp,cpo(3)
  RTYPE r1(3),r2(3),v1(3),v2(3),v3(3),v4(3),transform(5,5),det_a
  integer, INTENT (IN) :: rs
  RTYPE, INTENT (OUT) :: jac_x
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! cyclic peptide residues
  if (seqflag(rs).eq.5) then
!   CA
    s5(1) = x(at(rs)%bb(2))
    s5(2) = y(at(rs)%bb(2))
    s5(3) = z(at(rs)%bb(2))
!   N
    s4(1) = x(at(rs)%bb(1))
    s4(2) = y(at(rs)%bb(1))
    s4(3) = z(at(rs)%bb(1))
!   CD
    s3(1) = x(at(rs)%sc(4-shf))
    s3(2) = y(at(rs)%sc(4-shf))
    s3(3) = z(at(rs)%sc(4-shf))
!   CG
    s2(1) = x(at(rs)%sc(3-shf))
    s2(2) = y(at(rs)%sc(3-shf))
    s2(3) = z(at(rs)%sc(3-shf))
!   CB
    s1(1) = x(at(rs)%sc(2-shf))
    s1(2) = y(at(rs)%sc(2-shf))
    s1(3) = z(at(rs)%sc(2-shf))
! nucleotides
  else if (seqpolty(rs).eq.'N') then
    if (seqflag(rs).eq.24) then ! DIB, RIB, RIX, or DIX
      shfn = 2
    else
      shfn = 0
    end if
!   C3*
    s5(1) = x(nuci(rs,6-shfn))
    s5(2) = y(nuci(rs,6-shfn))
    s5(3) = z(nuci(rs,6-shfn))
!   C4*
    s4(1) = x(nuci(rs,5-shfn))
    s4(2) = y(nuci(rs,5-shfn))
    s4(3) = z(nuci(rs,5-shfn))
!   O4*
    s3(1) = x(at(rs)%sc(3))
    s3(2) = y(at(rs)%sc(3))
    s3(3) = z(at(rs)%sc(3))
!   C1*
    s2(1) = x(at(rs)%sc(2))
    s2(2) = y(at(rs)%sc(2))
    s2(3) = z(at(rs)%sc(2))
!   C2*
    s1(1) = x(at(rs)%sc(1))
    s1(2) = y(at(rs)%sc(1))
    s1(3) = z(at(rs)%sc(1))
  else
    write(ilog,*) 'Fatal. Called jacobian_pucker(...) with unsupported residue (',rs,' of&
 & type ',seqtyp(rs),': ',amino(seqtyp(rs)),').'
    call fexit()
  end if
!
  do k=1,3
    u1(k) = s2(k) - s1(k)
    u2(k) = s3(k) - s2(k)
    u3(k) = s4(k) - s3(k)
    u4(k) = s5(k) - s4(k)  
  end do
! 
  c1 = 0.0
  c2 = 0.0
  c3 = 0.0
  c4 = 0.0
  do k=1,3
    c1 = c1 + u1(k)*u1(k)
    c2 = c2 + u2(k)*u2(k)
    c3 = c3 + u3(k)*u3(k)
    c4 = c4 + u4(k)*u4(k)
  end do
  do k=1,3
    u1(k) = u1(k)/sqrt(c1)
    u2(k) = u2(k)/sqrt(c2)
    u3(k) = u3(k)/sqrt(c3)
    u4(k) = u4(k)/sqrt(c4)
  end do
!
  call crossprod(u2,u1,cpo,cp)
  w1 = cpo/cp
  call crossprod(u3,u2,cpo,cp)
  w2 = cpo/cp
  call crossprod(u4,u3,cpo,cp)
  w3 = cpo/cp
!
  r1 = s4-s2
  r2 = s4-s3
!
  call crossprod(u1,r1,v1,cp)
  call crossprod(u2,r2,v2,cp)
  call crossprod(w1,r1,v3,cp)
  call crossprod(w2,r2,v4,cp)
  do k=1,3
    transform(k,1) = v1(k)
    transform(k,2) = v2(k)
    transform(k,3) = v3(k)
    transform(k,4) = v4(k)
    transform(k,5) = 0.0
  end do
  call crossprod(u1,u4,t0,cp)
  transform(4,1) = t0(1)
  transform(5,1) = t0(2)
  call crossprod(u2,u4,t0,cp)
  transform(4,2) = t0(1)
  transform(5,2) = t0(2)
  call crossprod(w1,u4,t0,cp)
  transform(4,3) = t0(1)
  transform(5,3) = t0(2)
  call crossprod(w2,u4,t0,cp)
  transform(4,4) = t0(1)
  transform(5,4) = t0(2)
  call crossprod(w3,u4,t0,cp)
  transform(4,5) = t0(1)
  transform(5,5) = t0(2)
!
  call dtrm(transform,5,det_a,indx)
!
  jac_x = 1.0/abs(det_a)
!
end subroutine jacobian_pucker
!---------------------------------------------------------------------
!
subroutine pucker_refvals(rs,refvals)
! !------------------------------------------------------------
! ! subroutine to store reference values needed for ring closure
! ! INPUTS:
! !     rs -  reference residue for chain closure
! ! MODIFIES:
! !     oldvals - {r45,r51}
! !------------------------------------------------------------
!
  use iounit
  use polypep
  use fyoc
  use atoms
  use molecule
  use sequen
  use aminos
  use system
!
  implicit none
!
  integer rs,shf,shfn
  RTYPE refvals(2),getblen
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! cyclic peptide residues
  if (seqflag(rs).eq.5) then
!   CD-N
    refvals(1) = getblen(at(rs)%bb(1),at(rs)%sc(4-shf))
!   CG-CD
    refvals(2) = getblen(at(rs)%sc(4-shf),at(rs)%sc(3-shf))
! nucleotides
  else if (seqpolty(rs).eq.'N') then
    if (seqflag(rs).eq.24) then ! DIB, RIB, RIX, or DIX
      shfn = 2
    else
      shfn = 0
    end if
!   O4*-C4*
    refvals(1) = getblen(nuci(rs,5-shfn),at(rs)%sc(3))
!   C1*-O4*
    refvals(2) = getblen(at(rs)%sc(3),at(rs)%sc(2))
  else
    write(ilog,*) 'Fatal. Called pucker_refvals(...) with unsupported residue (',rs,' of&
 & type ',seqtyp(rs),': ',amino(seqtyp(rs)),').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------
!
! mode: 1 for dj, 2 for do
!
subroutine proline_cr(rsf,puckdofs,solmat,jacv,solcnt,curpks,mode)
!
  use fyoc
  use aminos
  use sequen
  use ujglobals
  use movesets
  use zmatrix
  use polypep
  use math
  use iounit
  use system
!
  implicit none
!
  integer i,j,solcnt,rsf,rs,k,rsoff,mode,shf
  RTYPE random,newphi,puckdofs(7),solmat(MAXSOLU*4,6),jacv(MAXSOLU*4),jacv2(MAXSOLU*4)
  RTYPE jaca,ztorbu,refvals(2),curpks(MAXSOLU*4,7,3),curpks2(MAXSOLU*4,7,3),solmat2(MAXSOLU*4,6)
  RTYPE normal
  logical discard,findsolu,badsolu(MAXSOLU*4),stayold
!
  badsolu(:) = .false.
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  k = 0
  if (mode.eq.1) then
    rsoff = 2
  else if (mode.eq.2) then
    rsoff = 1
  else
    write(ilog,*) 'Fatal. Called proline_cr(...) with unknown mode.&
 &Offending mode # is ',mode,'.'
    call fexit()
  end if
!
  do rs=rsf-rsoff,rsf
    if (pucline(rs).gt.0) then
      k = k + 1
      do i=1,solcnt
!       if prior failure we can save some work
        if (badsolu(i).EQV..true.) cycle
        discard = .false.
        if (rs.eq.rsf) newphi = ztor(at(rs)%sc(4-shf)) + solmat(i,5)
        if ((rs.eq.rsf-1).AND.(mode.eq.1)) newphi = ztor(at(rs)%sc(4-shf)) + solmat(i,3)
        if ((rs.eq.rsf-1).AND.(mode.eq.2)) newphi = ztor(at(rs)%sc(4-shf)) + solmat(i,2)
        if (rs.eq.rsf-2) newphi = ztor(at(rs)%sc(4-shf)) + solmat(i,1)
        if (newphi.gt.180.0) newphi = newphi - 360.0
        if (newphi.le.-180.0) newphi = newphi + 360.0
        if ((newphi.lt.-55.0).OR.(newphi.gt.50.0)) then
          badsolu(i) = .true.
          cycle
        end if
!       biased but necessary
        if (newphi.ge.0.0) then
          puckdofs(1) = -40.0 + normal()*10.0
        else
          puckdofs(1) = 35.0 + normal()*10.0
        end if
        call pucker_refvals(rs,refvals)
        stayold = .false.
        if (rs.eq.rsf) then
          puckdofs(2) = ztor(at(rs)%sc(4-shf)) + solmat(i,5)
          puckdofs(3) = ztor(pucline(rs)) + solmat(i,5)
          puckdofs(4) = ztor(at(rs)%sc(2-shf)) + solmat(i,5) + (random()-0.5)*3.0 ! decouple slightly with hc stpsz
          if (abs(solmat(i,5)).le.1.0e-9) then
!           leave original
            stayold = .true.
          end if
        else if (rs.eq.rsf-1) then
          if (mode.eq.1) then
            puckdofs(2) = ztor(at(rs)%sc(4-shf)) + solmat(i,3)
            puckdofs(3) = ztor(pucline(rs)) + solmat(i,3)
            puckdofs(4) = ztor(at(rs)%sc(2-shf)) + solmat(i,3) + (random()-0.5)*3.0
            if (abs(solmat(i,3)).le.1.0e-9) then
!             leave original
              stayold = .true.
            end if
          else if (mode.eq.2) then
            puckdofs(2) = ztor(at(rs)%sc(4-shf)) + solmat(i,2)
            puckdofs(3) = ztor(pucline(rs)) + solmat(i,2)
            puckdofs(4) = ztor(at(rs)%sc(2-shf)) + solmat(i,2) + (random()-0.5)*3.0
            if (abs(solmat(i,2)).le.1.0e-9) then
!             leave original
              stayold = .true.
            end if
          end if
        else if (rs.eq.rsf-2) then
          puckdofs(2) = ztor(at(rs)%sc(4-shf)) + solmat(i,1)
          puckdofs(3) = ztor(pucline(rs)) + solmat(i,1)
          puckdofs(4) = ztor(at(rs)%sc(2-shf)) + solmat(i,1) + (random()-0.5)*3.0
          if (abs(solmat(i,1)).le.1.0e-9) then
!           leave original
            stayold = .true.
          end if
        end if
        if (stayold.EQV..true.) then
          call jacobian_pucker(rs,jaca)
          jacv(i) = jacv(i)*jaca
          curpks(i,7,k) = bang(at(rs)%sc(4-shf))
          curpks(i,1,k) = ztor(pucline(rs))
          curpks(i,2,k) = ztor(at(rs)%sc(2-shf)) 
          curpks(i,3,k) = ztor(at(rs)%sc(3-shf)) 
          curpks(i,4,k) = ztor(at(rs)%sc(4-shf))
          curpks(i,5,k) = bang(at(rs)%sc(2-shf))
          curpks(i,6,k) = bang(at(rs)%sc(3-shf))
          ztorpr(at(rs)%sc(2-shf)) = curpks(i,2,k)
          ztorpr(at(rs)%sc(3-shf)) = curpks(i,3,k)
          ztorpr(at(rs)%sc(4-shf)) = curpks(i,4,k)
          bangpr(at(rs)%sc(2-shf)) = curpks(i,5,k)
          bangpr(at(rs)%sc(3-shf)) = curpks(i,6,k)
          bangpr(at(rs)%sc(4-shf)) = curpks(i,7,k)
          cycle 
        end if
        puckdofs(5) = bang(at(rs)%sc(2-shf))
        puckdofs(6) = bang(at(rs)%sc(3-shf)) + (random()-0.5)*1.0 ! small hc step size
        puckdofs(7) = 0.0
!
!       check for wrap-arounds
        do j=1,4
          if (puckdofs(j).gt.180.0) puckdofs(j) = puckdofs(j) - 360.0
          if (puckdofs(j).lt.-180.0) puckdofs(j) = puckdofs(j) + 360.0
        end do
        do j=5,6
          if ((puckdofs(j).gt.180.0).OR.(puckdofs(j).lt.0.0)) then
            write(ilog,*) 'Fatal. Bond angle out of bounds in proline_cr(...).'
            call fexit()
          end if
        end do
!       get differences
        puckdofs(1) = puckdofs(1) - ztor(at(rs)%sc(3-shf))
        puckdofs(2) = puckdofs(2) - ztor(at(rs)%sc(4-shf))
        puckdofs(3) = puckdofs(3) - ztor(pucline(rs))
        puckdofs(4) = puckdofs(4) - ztor(at(rs)%sc(2-shf))
        puckdofs(5) = puckdofs(5) - bang(at(rs)%sc(2-shf))
        puckdofs(6) = puckdofs(6) - bang(at(rs)%sc(3-shf))
!       apply the new set
!       double backup ztorpr(fline), backup and modify C(-1)-N-CA-C, check wrap-around
        ztorbu = ztorpr(pucline(rs))
        ztorpr(pucline(rs)) = ztor(pucline(rs))
        ztor(pucline(rs)) = ztor(pucline(rs)) + puckdofs(3)
        if (ztor(pucline(rs)).gt.180.0) ztor(pucline(rs)) = ztor(pucline(rs)) - 360.0
        if (ztor(pucline(rs)).lt.-180.0) ztor(pucline(rs)) = ztor(pucline(rs)) + 360.0
        curpks(i,1,k) = ztor(pucline(rs))
!       backup and modify C(-1)-N-CA-CB, check wrap-around
        ztorpr(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf))
        ztor(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf)) + puckdofs(4)
        if (ztor(at(rs)%sc(2-shf)).gt.180.0) ztor(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf)) - 360.0
        if (ztor(at(rs)%sc(2-shf)).lt.-180.0) ztor(at(rs)%sc(2-shf)) = ztor(at(rs)%sc(2-shf)) + 360.0
        curpks(i,2,k) = ztor(at(rs)%sc(2-shf))
!       backup and modify N-CA-CB-CG
        ztorpr(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf))
        ztor(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf)) + puckdofs(1)
        if (ztor(at(rs)%sc(3-shf)).gt.180.0) ztor(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf)) - 360.0
        if (ztor(at(rs)%sc(3-shf)).lt.-180.0) ztor(at(rs)%sc(3-shf)) = ztor(at(rs)%sc(3-shf)) + 360.0
        curpks(i,3,k) = ztor(at(rs)%sc(3-shf))
!       backup and modify CB-CA-N-CD 
        ztorpr(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf))
        ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) + puckdofs(2)
        if (ztor(at(rs)%sc(4-shf)).gt.180.0) ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) - 360.0
        if (ztor(at(rs)%sc(4-shf)).lt.-180.0) ztor(at(rs)%sc(4-shf)) = ztor(at(rs)%sc(4-shf)) + 360.0
        curpks(i,4,k) = ztor(at(rs)%sc(4-shf))
!       backup and modify angle N-CA-CB
        bangpr(at(rs)%sc(2-shf)) = bang(at(rs)%sc(2-shf))
        bang(at(rs)%sc(2-shf)) = bang(at(rs)%sc(2-shf)) + puckdofs(5)
        curpks(i,5,k) = bang(at(rs)%sc(2-shf))
!       backup and modify angle CA-CB-CG
        bangpr(at(rs)%sc(3-shf)) = bang(at(rs)%sc(3-shf))
        bang(at(rs)%sc(3-shf)) = bang(at(rs)%sc(3-shf)) + puckdofs(6)
        curpks(i,6,k) = bang(at(rs)%sc(3-shf))
!       backup angle CB-CG-CD
        bangpr(at(rs)%sc(4-shf)) = bang(at(rs)%sc(4-shf))
!       re-build xyz
        call makexyz_forsc(rs)
!       now close the 5-ring (populates puckdofs(5))
        call ringclose(rs,refvals,puckdofs(1:7),ztor(at(rs)%sc(4-shf))/RADIAN,findsolu)
        if ((findsolu.EQV..true.).AND.(puckdofs(7).lt.130.0).AND.(puckdofs(7).gt.90.0)) then
!         lastly modify angle CA-N-CD
          bang(at(rs)%sc(4-shf)) = puckdofs(7)
          curpks(i,7,k) = bang(at(rs)%sc(4-shf))
        else
          if (findsolu.EQV..false.) then
            pc_wrncnt(5) = pc_wrncnt(5) + 1
            if (pc_wrncnt(5).eq.pc_wrnlmt(5)) then
              write(ilog,*) 'Warning. Failed to find a new ring conformation in proline_cr(...). This &
 &is usually indicative of an incompatible phi-angle perturbation, but may report on more serious &
 &problems for instance due to bond angle exceptions as well.'
              write(ilog,*) 'This was warning #',pc_wrncnt(5),' of this type not all of which may be displayed.'
              if (10.0*pc_wrnlmt(5).gt.0.5*HUGE(pc_wrnlmt(5))) then
                pc_wrncnt(5) = 0 ! reset
              else
                pc_wrnlmt(5) = pc_wrnlmt(5)*10
              end if
            end if
          end if
          bang(at(rs)%sc(4-shf)) = puckdofs(7)
          curpks(i,7,k) = bang(at(rs)%sc(4-shf))
          curpks(i,1,k) = ztorpr(pucline(rs))
          curpks(i,2,k) = ztorpr(at(rs)%sc(2-shf)) 
          curpks(i,3,k) = ztorpr(at(rs)%sc(3-shf)) 
          curpks(i,4,k) = ztorpr(at(rs)%sc(4-shf))
          curpks(i,5,k) = bangpr(at(rs)%sc(2-shf))
          curpks(i,6,k) = bangpr(at(rs)%sc(3-shf)) 
          ztor(pucline(rs)) =  curpks(i,1,k)
          ztor(at(rs)%sc(2-shf)) = curpks(i,2,k)
          ztor(at(rs)%sc(3-shf)) = curpks(i,3,k)
          ztor(at(rs)%sc(4-shf)) = curpks(i,4,k)
          bang(at(rs)%sc(2-shf)) = curpks(i,5,k)
          bang(at(rs)%sc(3-shf)) = curpks(i,6,k)
          call restore_forring(rs)
          ztorpr(pucline(rs)) = ztorbu
          call makexyz_forsc(rs)
          badsolu(i) = .true.
          cycle
        end if
!       re-build xyz
        call makexyz_forsc(rs)
!       the jacobian relies exclusively on xyz
        call jacobian_pucker(rs,jaca)
        jacv(i) = jacv(i)*jaca
        call restore_forring(rs)
        ztorpr(pucline(rs)) = ztorbu
        call makexyz_forsc(rs)  
      end do
    else if (seqflag(rs).eq.5) then
      pc_wrncnt(3) = pc_wrncnt(3) + 1
      if (pc_wrncnt(3).eq.pc_wrnlmt(3)) then
        write(ilog,*) 'Warning. Called proline_cr(...) with an apparent cyclic residue without proper &
 &puckering setup. This is likely to lead to topology violations and may indicate a bug.'
        write(ilog,*) 'This was warning #',pc_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*pc_wrnlmt(3).gt.0.5*HUGE(pc_wrnlmt(3))) then
          pc_wrncnt(3) = 0 ! reset
        else
          pc_wrnlmt(3) = pc_wrnlmt(3)*10
        end if
      end if
    end if
  end do
!
  jacv2(1:solcnt) = jacv(1:solcnt)
  curpks2(1:solcnt,:,:) = curpks(1:solcnt,:,:)
  solmat2(1:solcnt,:) = solmat(1:solcnt,:)
  k = 0
  do i=1,solcnt
    if (badsolu(i).EQV..false.) then
      k = k + 1
      jacv(k) = jacv2(i)
      solmat(k,:) = solmat2(i,:)
      curpks(k,:,:) = curpks2(i,:,:)
    end if
  end do
  solcnt = k
!
end subroutine proline_cr
!
!---------------------------------------------------------------------
!
subroutine ringclose(rs,refvals,d_pucker,w1val,findsolu)
! !------------------------------------------------------------
! ! Subroutine to calculate the constraints for chain closure.
! ! INPUTS:
! !     rs -  index of ring for chain closure
! !     oldvals - constant dof i.e. bond lengths
! ! OUTPUT:
! !     d_dependents - array of change in angles needed for closure
! !------------------------------------------------------------
!
  use iounit
  use polypep
  use atoms
  use math
  use fyoc
  use movesets
  use ujglobals
  use aminos
  use mcsums
  use molecule
  use sequen
  use fyoc
  use zmatrix
  use system
!
  implicit none
!
  integer rs,shf,shfn
  logical findsolu,logi
  RTYPE w1val,nopl,q0(3)
  RTYPE r(3),a3,ax,s(3),con1xy(3),a1,ay,Ra1(3,3),Ra2(3,3),Ra3(3,3),az
  RTYPE newxaxis(3),newyaxis(3),newzaxis(3),yaxis(3),zaxis(3),a2
  RTYPE xaxis(3),RaE(3,3),or_pl(3),or_pl2(3)
  RTYPE origin(3),s1(3),posC3(3),posC4(3),posC2(3),posC1(3)
  RTYPE con32(3),con43(3),d_pucker(7),refvals(2)
  RTYPE p43,p15,p54
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! get the current position of atoms used in ring closure
!
  origin(1) = 0.0
  origin(2) = 0.0
  origin(3) = 0.0
! cyclic peptide residues
  if (seqflag(rs).eq.5) then
!   CB
    posC2(1) = x(at(rs)%sc(2-shf))
    posC2(2) = y(at(rs)%sc(2-shf))
    posC2(3) = z(at(rs)%sc(2-shf))
!   CA
    posC3(1) = x(at(rs)%bb(2))
    posC3(2) = y(at(rs)%bb(2))
    posC3(3) = z(at(rs)%bb(2))
!   N
    posC4(1) = x(at(rs)%bb(1))
    posC4(2) = y(at(rs)%bb(1))
    posC4(3) = z(at(rs)%bb(1))
!   CG
    posC1(1) = x(at(rs)%sc(3-shf))
    posC1(2) = y(at(rs)%sc(3-shf))
    posC1(3) = z(at(rs)%sc(3-shf))
! nucleotides
  else if (seqpolty(rs).eq.'N') then
    if (seqflag(rs).eq.24) then ! DIB, RIB, RIX, or DIX
      shfn = 2
    else
      shfn = 0
    end if
!   C2*
    posC2(1) = x(at(rs)%sc(1))
    posC2(2) = y(at(rs)%sc(1))
    posC2(3) = z(at(rs)%sc(1))
!   C3*
    posC3(1) = x(nuci(rs,6-shfn))
    posC3(2) = y(nuci(rs,6-shfn))
    posC3(3) = z(nuci(rs,6-shfn))
!   C4*
    posC4(1) = x(nuci(rs,5-shfn))
    posC4(2) = y(nuci(rs,5-shfn))
    posC4(3) = z(nuci(rs,5-shfn))
!   C1*
    posC1(1) = x(at(rs)%sc(2))
    posC1(2) = y(at(rs)%sc(2))
    posC1(3) = z(at(rs)%sc(2))
  else
    write(ilog,*) 'Fatal. Called ringclose(...) with unsupported residue (',rs,' of&
 & type ',seqtyp(rs),': ',amino(seqtyp(rs)),').'
    call fexit()
  end if
!
  con32(:) = posC2(:)-posC3(:)
  con43(:) = posC4(:)-posC3(:)
!
! get bond lengths
  p43 = sqrt(sum((posC3(:)-posC4(:))**2))
  p54 = refvals(1)
  p15 = refvals(2)
!
  xaxis = (/ 1, 0, 0 /)
  yaxis = (/ 0, 1, 0 /)
  zaxis = (/ 0, 0, 1 /)
!
! ! projection of con1 on xy-plane
  con1xy(1) = posC4(1) - posC3(1)
  con1xy(2) = posC4(2) - posC3(2)
  con1xy(3) = 0.0
!
! ! rotate around Z axis
  if (abs(sum(con1xy(1:3))).le.1.0e-9) then
    a1 = 0.0
  else
    call bondang2(xaxis,origin,con1xy,a1)
! ! compare to y-axis to find if +/- rotation
    call bondang2(con1xy,origin,yaxis,ay)
    if(ay.gt.(PI/2.0)) then
      a1 = -a1
    end if
  end if
!  
! ! rotation matrix for z-rotation     
  Ra1(1,1) = cos(a1)
  Ra1(1,2) = sin(a1)
  Ra1(1,3) = 0.0
  Ra1(2,1) = -Ra1(1,2)
  Ra1(2,2) = Ra1(1,1)
  Ra1(2,3) = 0.0
  Ra1(3,1) = 0.0
  Ra1(3,2) = 0.0
  Ra1(3,3) = 1.0
!
! ! rotate around new Y axis (angle between new x-axis and con1)
  newxaxis(1) = Ra1(1,1)
  newxaxis(2) = Ra1(1,2)
  newxaxis(3) = Ra1(1,3)
  call bondang2(con43,origin,newxaxis,a2)
  call bondang2(con43,origin,zaxis,az)
  if(az.lt.(PI/2)) then
    a2 = -a2
  end if
!
! ! rotation matrix for y-rotation
  Ra2(1,1) = cos(a2)
  Ra2(1,2) = 0.0
  Ra2(1,3) = -sin(a2)
  Ra2(2,1) = 0.0
  Ra2(2,2) = 1.0
  Ra2(2,3) = 0.0
  Ra2(3,1) = -Ra2(1,3)
  Ra2(3,2) = 0.0
  Ra2(3,3) = Ra2(1,1)
!
! ! rotate around new X axis so that Y'' and p1 align
  RaE = Matmul(Ra2,Ra1)
  newzaxis(1) = RaE(3,1)
  newzaxis(2) = RaE(3,2)
  newzaxis(3) = RaE(3,3)
  newxaxis(1) = RaE(1,1)
  newxaxis(2) = RaE(1,2)
  newxaxis(3) = RaE(1,3)
  newyaxis(1) = Ra1(2,1)
  newyaxis(2) = Ra1(2,2)
  newyaxis(3) = Ra1(2,3)
  call crossprod(con32,con43,or_pl,nopl)
  call bondang2(newzaxis,origin,or_pl,a3)
  call crossprod(con32,con43,or_pl2,nopl)
  call bondang2(newyaxis,origin,or_pl2,ax)
  if(ax.lt.(PI/2)) then
    a3 = -a3
  end if
!
  Ra3(1,1) = 1.0
  Ra3(1,2) = 0.0
  Ra3(1,3) = 0.0
  Ra3(2,1) = 0.0
  Ra3(2,2) = cos(a3)
  Ra3(2,3) = sin(a3)
  Ra3(3,1) = 0.0
  Ra3(3,2) = -Ra3(2,3)
  Ra3(3,3) = Ra3(2,2)
!
! ! build full euler rotation matrix
  RaE = Matmul(Ra3,Matmul(Ra2,Ra1))
!
! ! set q0 = p43
  q0(1) = p43
  q0(2) = 0.0
  q0(3) = 0.0
!
! ! get pos S = RaE.posC1
  s1(1) = posC1(1)-posC3(1)
  s1(2) = posC1(2)-posC3(2)
  s1(3) = posC1(3)-posC3(3)
  s(1) = sum(RaE(1,:)*s1(:))
  s(2) = sum(RaE(2,:)*s1(:))
  s(3) = sum(RaE(3,:)*s1(:))
!
  findsolu = .false.
!
  call checkreal(w1val,logi)
!
  if (logi.EQV..true.) then
    call alpha1_of_w1(w1val)
    findsolu = .true.
  else
!   cyclic peptide residues
    if (seqflag(rs).eq.5) then
      d_pucker(7) = bang(at(rs)%sc(4-shf))
    else if (seqpolty(rs).eq.'N') then
      d_pucker(7) = bang(at(rs)%sc(3))
    end if
  end if
!    
  CONTAINS
! -------------------------------------------------------------------
! Following subroutines are helpers for the UJchainclosure method
!
! -------------------------------------------------------------------  
!
! check for real solutions
subroutine checkreal(w1val,isreal)
!
  RTYPE w1val, arg1, arg2
  logical isreal
!
  isreal = .false.
  call r_of_w1(r,w1val)
!
  arg1 = 4*(p54**2)*(r(1)**2+r(2)**2)
  arg2 = (r(1)**2 + r(2)**2 + r(3)**2 + p54**2 - p15**2)**2
!
  if(arg1.ge.arg2) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
end subroutine checkreal
! -------------------------------------------------------------------  
!
! subroutine to get r(w1) = inv(R1)*(s-q0)
subroutine r_of_w1(r,w1val)
!
  RTYPE r(3),w1val
!
  r(1) = s(1) - q0(1)
  r(2) = cos(w1val)*(s(2) - q0(2)) + sin(w1val)*(s(3) - q0(3))
  r(3) = -sin(w1val)*(s(2) - q0(2)) + cos(w1val)*(s(3) - q0(3))
!
end subroutine r_of_w1
! -------------------------------------------------------------------  
!
! subroutine to find alpha1
subroutine alpha1_of_w1(w1val)
!
  RTYPE w1val,w,arg1,arg2,arg1neg,arg2neg,neg
!
  call r_of_w1(r,w1val)
  w = r(1)**2 + r(2)**2 + r(3)**2 + p54**2 - p15**2
!     
    arg1neg=(w*r(1) + r(2)*sqrt(4*p54**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p54*(r(1)**2 + r(2)**2))
    arg2neg=(w*r(2) - r(1)*sqrt(4*p54**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p54*(r(1)**2 + r(2)**2))
!
  arg1=(w*r(1) - r(2)*sqrt(4*p54**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p54*(r(1)**2 + r(2)**2))
  arg2=(w*r(2) + r(1)*sqrt(4*p54**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p54*(r(1)**2 + r(2)**2))
!  
  d_pucker(7) = 180.0 + atan2(arg2,arg1)*RADIAN
  neg = 180.0 + atan2(arg2neg,arg1neg)*RADIAN
!
end subroutine alpha1_of_w1
! -------------------------------------------------------------------  
!
! ! subroutine to generate T(i) and its inverse (invTi)
subroutine genT(alpha,Ti,invTi)
!
  RTYPE Ti(3,3),invTi(3,3)
  RTYPE, INTENT(IN) :: alpha
!
  Ti(1,1) = cos(alpha)
  Ti(1,2) = -sin(alpha)
  Ti(1,3) = 0.0
  Ti(2,1) = -1*Ti(1,2)
  Ti(2,2) = Ti(1,1)
  Ti(2,3) = 0.0
  Ti(3,1) = 0.0
  Ti(3,2) = 0.0
  Ti(3,3) = 1.0
!
  invTi = Ti
  invTi(1,2) = Ti(2,1)
  invTi(2,1) = Ti(1,2)
!
end subroutine genT
!
! -------------------------------------------------------------------  
!
! ! subroutine to generate R(i) and its inverse (invRi)
subroutine genR(phi1,Ri,invRi)
!
  RTYPE Ri(3,3),invRi(3,3)
  RTYPE, INTENT(IN) :: phi1
!
  Ri(1,1) = 1.0
  Ri(1,2) = 0.0
  Ri(1,3) = 0.0
  Ri(2,1) = 0.0
  Ri(2,2) = cos(phi1)
  Ri(2,3) = -sin(phi1)
  Ri(3,1) = 0.0
  Ri(3,2) = -Ri(2,3)
  Ri(3,3) = Ri(2,2)
!
  invRi = Ri
  invRi(2,3) = Ri(3,2)
  invRi(3,2) = Ri(2,3)
!
end subroutine genR
! ------------------------------------------------------------------
end subroutine ringclose
