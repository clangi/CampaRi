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
! MAIN AUTHOR:   Adam Steffen                                              !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
subroutine buildUJdof(rsi,rsf,dof,ujmat)
! !------------------------------------------------------------
! ! Subroutine to calculate the matrix I = (dA/dphi)**2 
! ! for UJ prerotation.
! ! 
! ! INPUTS:     rsi, rsf - inital and final fixed residues number
! ! MODIFIES:   ujmat - five DOF per residue, except for last 2
! !------------------------------------------------------------
!
  use iounit
  use ujglobals
  use polypep
  use atoms
  use math
  use fyoc
  use movesets
  use zmatrix
  use sequen
!
  implicit none
!
  integer j,k,m,rsi,rsf,dof,rs
  RTYPE pos2(3),pos3(3),pos4(3)
  RTYPE ujmat(MAXUJDOF,MAXUJDOF),doflever(3,MAXUJDOF),dvec(3)
!
! the reference atom is the CA of residue rsf-1
!
! loop over all free residues of the prerotation
  do rs=rsi,rsf-2
!
!   phi angle
    pos4(1) = x(ci(rs-1))
    pos4(2) = y(ci(rs-1))
    pos4(3) = z(ci(rs-1))
    pos3(1) = x(ni(rs))
    pos3(2) = y(ni(rs))
    pos3(3) = z(ni(rs))
    pos2(1) = x(cai(rs))
    pos2(2) = y(cai(rs))
    pos2(3) = z(cai(rs))
    call get_tor_lever(molofrs(atmres(cai(rsf-1))),cai(rsf-1),&
 &pos2,pos3,pos4,dvec)
    doflever(:,2*(rs-rsi)+1) = dvec(:)/RADIAN
!
!   psi angle
    pos4(1) = x(ni(rs))
    pos4(2) = y(ni(rs))
    pos4(3) = z(ni(rs))
    pos3(1) = x(cai(rs))
    pos3(2) = y(cai(rs))
    pos3(3) = z(cai(rs))
    pos2(1) = x(ci(rs))
    pos2(2) = y(ci(rs))
    pos2(3) = z(ci(rs))
    call get_tor_lever(molofrs(atmres(cai(rsf-1))),cai(rsf-1),&
 &pos2,pos3,pos4,dvec)
    doflever(:,2*(rs-rsi)+2) = dvec(:)/RADIAN
!
!   N-CA-C bond angle
    pos4(1) = x(ci(rs-1))
    pos4(2) = y(ci(rs-1))
    pos4(3) = z(ci(rs-1))
    pos3(1) = x(ni(rs))
    pos3(2) = y(ni(rs))
    pos3(3) = z(ni(rs))
    pos2(1) = x(cai(rs))
    pos2(2) = y(cai(rs))
    pos2(3) = z(cai(rs))
    call get_ang_lever(molofrs(atmres(cai(rsf-1))),cai(rsf-1),&
 &ci(rs),pos2,pos3,pos4,dvec)
    doflever(:,3*(rs-rsi)+ 2*cur_ujsz + 1) = dvec(:)/RADIAN
!
!   CA-C-N(+1) bond angle
    pos4(1) = x(ni(rs))
    pos4(2) = y(ni(rs))
    pos4(3) = z(ni(rs))
    pos3(1) = x(cai(rs))
    pos3(2) = y(cai(rs))
    pos3(3) = z(cai(rs))
    pos2(1) = x(ci(rs))
    pos2(2) = y(ci(rs))
    pos2(3) = z(ci(rs))
    call get_ang_lever(molofrs(atmres(cai(rsf-1))),cai(rsf-1),&
 &ni(rs+1),pos2,pos3,pos4,dvec)
    doflever(:,3*(rs-rsi)+ 2*cur_ujsz + 2) = dvec(:)/RADIAN
!
!   C-N(+1)-CA(+1) bond angle
    pos4(1) = x(cai(rs))
    pos4(2) = y(cai(rs))
    pos4(3) = z(cai(rs))
    pos3(1) = x(ci(rs))
    pos3(2) = y(ci(rs))
    pos3(3) = z(ci(rs))
    pos2(1) = x(ni(rs+1))
    pos2(2) = y(ni(rs+1))
    pos2(3) = z(ni(rs+1))
    call get_ang_lever(molofrs(atmres(cai(rsf-1))),cai(rsf-1),&
 &cai(rs+1),pos2,pos3,pos4,dvec)
    doflever(:,3*(rs-rsi)+ 2*cur_ujsz + 3) = dvec(:)/RADIAN
!
  end do
!
! ! build symetric matrix I = dA/dphi*dA/dphi (outer_product)
  do j=1,dof
    do k=1,dof
!     ! inner product
      ujmat(j,k) = 0.0
      do m=1,3
        ujmat(j,k) = ujmat(j,k) + doflever(m,j)*doflever(m,k)
      end do
    end do
  end do
!
end subroutine buildUJdof
!
!
!-----------------------------------------------------------------------
!
! a simple subroutine to obtain the derivative dr/dphi for atom i (down
! the chain) with respect to a torsion earlier in the chain defined by
! positions pos4 (furthest), pos3, pos2 and an unneeded posx
!
subroutine get_tor_lever(imol,i,pos2,pos3,pos4,dvec)
!
  use iounit
  use forces
  use atoms
  use molecule
  use sequen
!
  implicit none
!
  integer j
!
  integer imol,i
  RTYPE pos1(3),pos2(3),pos3(3),pos4(3),con1(3),con2(3)
  RTYPE or_pl(3),in_pl(3),nopl,nipl,c1,c2,alph,torstar,effbond
  RTYPE cosine,sine,eps,dax,dvec(3)
!
#ifdef DISABLE_FLOAT
  eps = 0.0000001d0
#else
! the low-precision xyz-coordinates lead to some hysteresis error in converting back and forth
  eps = 0.02d0
#endif
!
  pos1(1) = x(i)
  pos1(2) = y(i)
  pos1(3) = z(i)
!
  call dihed(pos1,pos2,pos3,pos4,torstar)
!
! now the derivative is simply dref/dtorstar * dtorstar/dphi
! - the latter is always unity
! - the former can be computed from the analytical formula in genxyz
!
! first get normalized connection vectors 23 and 34 to set up frame
  c1 = 0.0
  c2 = 0.0
  do j=1,3
    con1(j) = pos2(j) - pos3(j)
    c1 = c1 + con1(j)*con1(j)
    con2(j) = pos3(j) - pos4(j)
    c2 = c2 + con2(j)*con2(j)
  end do
  c1 = sqrt(c1)
  c2 = sqrt(c2)
  do j=1,3
    con1(j) = con1(j)/c1
    con2(j) = con2(j)/c2
  end do
!
! a normal vector to the reference plane is given by 23x34
  call crossprod(con2,con1,or_pl,nopl)
! normalize it
  cosine = con1(1)*con2(1)+con1(2)*con2(2)+con1(3)*con2(3)
  sine = sqrt(max(1.0d0-cosine**2,eps))
  if (abs(cosine) .ge. 1.0d0) then
    write(ilog,10)  i,torstar,cosine
 10     format (/,' Bad torsion for residue',i6,': ',2f12.4)
  end if
  do j=1,3
    or_pl(j) = or_pl(j)/sine
  end do
! and the vector within the reference plane is given by (23x34)x34
  call crossprod(or_pl,con1,in_pl,nipl)
!  
! the effective bond angles and bond lengths are readily calculated
  do j=1,3
    con1(j) = pos2(j) - pos1(j)
  end do
  effbond = sqrt(con1(1)**2 + con1(2)**2 + con1(3)**2)
  call bondang2(pos1,pos2,pos3,alph)
  dax = sin(alph)*effbond
!
  do j=1,3
    dvec(j) = &
 & effbond*(-in_pl(j)*sin(alph)*sin(torstar) +&
 &         or_pl(j)*sin(alph)*cos(torstar))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! a slightly more complicated routine to obtain the derivative dr/dalph for atom ii (down
! the chain) with respect to a bond angle earlier in the chain.
! a constant torsion is defined by positions pos4 (furthest), pos3, pos2 and posx (atom kk),
! the bond angle (a1) is between pos3, pos2 and posx
! from this the derivative dposx/da1 is computed
! then (if necessary) the frame is shifted one up and the derivative of the same formula is
! used on pos3, pos2, posx, and pos1 now determining the (3x3)-derivative dpos1/dposx
! the latter is then matrix-multiplied with dposx/da1 to obtain the actual quantity of
! interest, i.e., dpos1/da1
!
subroutine get_ang_lever(imol,ii,kk,pos2,pos3,pos4,dvec)
!
  use iounit
  use forces
  use atoms
  use molecule
  use sequen
  use math
  use zmatrix
!
  implicit none
!
  integer j,k
  integer imol,ii,kk
  RTYPE pos1(3),pos2(3),pos3(3),pos4(3),con1(3),con2(3),conx(3)
  RTYPE or_pl(3),in_pl(3),nopl,nipl,c1,c2,alph,torstar,effbond
  RTYPE cosine,sine,eps,dvec(3),posx(3)
  RTYPE dor_pl(3,3),din_pl(3,3)
  RTYPE dc1dpx(3,3),dsindpx(3),dpxda(3),defbdpx(3),dvecp(3,3)
!
#ifdef DISABLE_FLOAT
  eps = 0.0000001d0
#else
! the low-precision xyz-coordinates lead to some hysteresis error in converting back and forth
  eps = 0.02d0
#endif
!
  posx(1) = x(kk)
  posx(2) = y(kk)
  posx(3) = z(kk)
!
  call dihed(posx,pos2,pos3,pos4,torstar)
!
! first get normalized connection vectors 23 and 34 to set up frame
  c1 = 0.0
  c2 = 0.0
  do j=1,3
    con1(j) = pos2(j) - pos3(j)
    c1 = c1 + con1(j)*con1(j)
    con2(j) = pos3(j) - pos4(j)
    c2 = c2 + con2(j)*con2(j)
  end do
  c1 = sqrt(c1)
  c2 = sqrt(c2)
  do j=1,3
    con1(j) = con1(j)/c1
    con2(j) = con2(j)/c2
  end do
!
! a normal vector to the reference plane is given by 23x34
  call crossprod(con2,con1,or_pl,nopl)
! normalize it
  cosine = con1(1)*con2(1)+con1(2)*con2(2)+con1(3)*con2(3)
  sine = sqrt(max(1.0d0-cosine**2,eps))
  if (abs(cosine) .ge. 1.0d0) then
    write(ilog,10)  kk,torstar,cosine
 10     format (/,' Bad torsion for residue',i6,': ',2f12.4)
  end if
  do j=1,3
    or_pl(j) = or_pl(j)/sine
  end do
! and the vector within the reference plane is given by (23x34)x34
  call crossprod(or_pl,con1,in_pl,nipl)
!  
! the effective bond angles and bond lengths are readily calculated
  call bondang2(pos3,pos2,posx,alph)
  do j=1,3
    conx(j) = posx(j) - pos2(j)
  end do
  effbond = sqrt(conx(1)**2 + conx(2)**2 + conx(3)**2)
!
  do j=1,3
    dvec(j) = effbond*(&
 &          in_pl(j)*cos(alph)*cos(torstar) +&
 &          or_pl(j)*cos(alph)*sin(torstar) +&
 &          con1(j)*sin(alph))
  end do
!
! if the original 4th atom is the same as the reference, we're done 
  if (ii.eq.kk) return
!
!
!
! if not, we have to repeat with a shifted frame and take the derivative
! with respect to the position of the (now) 3rd position whose partial
! is in dvec
!
  dpxda(:) = dvec(:)
  pos1(1) = x(ii)
  pos1(2) = y(ii)
  pos1(3) = z(ii)
!
  call dihed(pos1,posx,pos2,pos3,torstar)
!
! first get normalized connection vectors 23 and 34 to set up frame
  con1(:) = posx(:) - pos2(:)
  c1 = sqrt(sum(con1(:)*con1(:)))
  con2(:) = pos2(:) - pos3(:)
  c2 = sqrt(sum(con2(:)*con2(:)))
!  vector derivative dcon1(norm.)/dposx 
  dc1dpx(1,1) = -((con1(1)*con1(1))/(c1**3)) + 1.0/c1
  dc1dpx(1,2) = -con1(1)*con1(2)/(c1**3)
  dc1dpx(1,3) = -con1(1)*con1(3)/(c1**3)
  dc1dpx(2,1) = -con1(2)*con1(1)/(c1**3)
  dc1dpx(2,2) = -((con1(2)*con1(2))/(c1**3)) + 1.0/c1
  dc1dpx(2,3) = -con1(2)*con1(3)/(c1**3)
  dc1dpx(3,1) = -con1(3)*con1(1)/(c1**3)
  dc1dpx(3,2) = -con1(3)*con1(2)/(c1**3)
  dc1dpx(3,3) = -((con1(3)*con1(3))/(c1**3)) + 1.0/c1
  con1(:) = con1(:)/c1
  con2(:) = con2(:)/c2
!
! a normal vector to the reference plane is given by 23x34
  call crossprod(con2,con1,or_pl,nopl)
  dor_pl(1,1) = con2(2)*dc1dpx(3,1) - con2(3)*dc1dpx(2,1)
  dor_pl(1,2) = con2(2)*dc1dpx(3,2) - con2(3)*dc1dpx(2,2)
  dor_pl(1,3) = con2(2)*dc1dpx(3,3) - con2(3)*dc1dpx(2,3)
  dor_pl(2,1) = con2(3)*dc1dpx(1,1) - con2(1)*dc1dpx(3,1)
  dor_pl(2,2) = con2(3)*dc1dpx(1,2) - con2(1)*dc1dpx(3,2)
  dor_pl(2,3) = con2(3)*dc1dpx(1,3) - con2(1)*dc1dpx(3,3)
  dor_pl(3,1) = con2(1)*dc1dpx(2,1) - con2(2)*dc1dpx(1,1)
  dor_pl(3,2) = con2(1)*dc1dpx(2,2) - con2(2)*dc1dpx(1,2)
  dor_pl(3,3) = con2(1)*dc1dpx(2,3) - con2(2)*dc1dpx(1,3)
! normalize it
  cosine = con1(1)*con2(1)+con1(2)*con2(2)+con1(3)*con2(3)
  sine = sqrt(max(1.0d0-cosine**2,eps))
  dsindpx(1) = (-1.0/sine)*cosine*sum(dc1dpx(1,:)*con2(:))
  dsindpx(2) = (-1.0/sine)*cosine*sum(dc1dpx(2,:)*con2(:))
  dsindpx(3) = (-1.0/sine)*cosine*sum(dc1dpx(3,:)*con2(:))
  if (abs(cosine) .ge. 1.0d0) then
    write(ilog,10)  ii,torstar,cosine
  end if
  do j=1,3
    do k=1,3
      dor_pl(j,k) = dor_pl(j,k)/sine - or_pl(j)*dsindpx(k)/(sine**2)
    end do
    or_pl(j) = or_pl(j)/sine
  end do
!
! and the vector within the reference plane is given by (23x34)x34
  call crossprod(or_pl,con1,in_pl,nipl)
  din_pl(1,1) = or_pl(2)*dc1dpx(3,1) + dor_pl(2,1)*con1(3)&
 &           -  or_pl(3)*dc1dpx(2,1) - dor_pl(3,1)*con1(2)
  din_pl(1,2) = or_pl(2)*dc1dpx(3,2) + dor_pl(2,2)*con1(3)&
 &           -  or_pl(3)*dc1dpx(2,2) - dor_pl(3,2)*con1(2)
  din_pl(1,3) = or_pl(2)*dc1dpx(3,3) + dor_pl(2,3)*con1(3)&
 &           -  or_pl(3)*dc1dpx(2,3) - dor_pl(3,3)*con1(2)
  din_pl(2,1) = or_pl(3)*dc1dpx(1,1) + dor_pl(3,1)*con1(1)&
 &           -  or_pl(1)*dc1dpx(3,1) - dor_pl(1,1)*con1(3)
  din_pl(2,2) = or_pl(3)*dc1dpx(1,2) + dor_pl(3,2)*con1(1)&
 &           -  or_pl(1)*dc1dpx(3,2) - dor_pl(1,2)*con1(3)
  din_pl(2,3) = or_pl(3)*dc1dpx(1,3) + dor_pl(3,3)*con1(1)&
 &           -  or_pl(1)*dc1dpx(3,3) - dor_pl(1,3)*con1(3)
  din_pl(3,1) = or_pl(1)*dc1dpx(2,1) + dor_pl(1,1)*con1(2)&
 &           -  or_pl(2)*dc1dpx(1,1) - dor_pl(2,1)*con1(1)
  din_pl(3,2) = or_pl(1)*dc1dpx(2,2) + dor_pl(1,2)*con1(2)&
 &           -  or_pl(2)*dc1dpx(1,2) - dor_pl(2,2)*con1(1)
  din_pl(3,3) = or_pl(1)*dc1dpx(2,3) + dor_pl(1,3)*con1(2)&
 &           -  or_pl(2)*dc1dpx(1,3) - dor_pl(2,3)*con1(1)
!  
! the effective bond angles and bond lengths are readily calculated
  call bondang2(pos2,posx,pos1,alph)
  do j=1,3
    conx(j) = pos1(j) - posx(j)
  end do
  effbond = sqrt(conx(1)**2 + conx(2)**2 + conx(3)**2)
  defbdpx(1) = (-1.0/effbond)*conx(1)
  defbdpx(2) = (-1.0/effbond)*conx(2)
  defbdpx(3) = (-1.0/effbond)*conx(3)
!
! note that we assume the geometry to remain fixed which is why partials
! with respect to torstar, alph and effbond are omitted (all of which explicitly 
! depend on posx)
! the addition of the unity matrix is reflecting the translational component of
! the dpos1/dposx derivative while the rest is accoungting for rotation
  do j=1,3
    do k=1,3
      dvecp(j,k) = effbond*(&
 &       din_pl(j,k)*sin(alph)*cos(torstar) +&
! &       in_pl(j)*cos(alph)*fvx(k)*cos(torstar) -
! &       in_pl(j)*sin(alph)*sin(torstar)*tvx(k) +&
 &       dor_pl(j,k)*sin(alph)*sin(torstar)&! +
! &       or_pl(j)*cos(alph)*fvx(k)*sin(torstar) +
! &       or_pl(j)*sin(alph)*cos(torstar)*tvx(k) -&
 &       - dc1dpx(j,k)*cos(alph))
! &       + con1(j)*sin(alph)*fvx(k)) + 
! &            defbdpx(k)*(
! &       in_pl(j)*sin(alph)*cos(torstar) +
! &       or_pl(j)*sin(alph)*sin(torstar) -
! &       con1(j)*cos(alph))
    end do
    dvecp(j,j) = dvecp(j,j) + 1.0
  end do
!
! finally the last chain rule derivative dpos1/dposx*dposx/dalpha = dpos1/dalpha
! (where alpha is the original bond angle down the chain)
  do k=1,3
    dvec(k) = sum(dvecp(k,:)*dpxda(:))
  end do
!
end
!
! !-----------------------------------------------------------------------
! ! Subroutine that finds the jacobian factor of the 5x5 matrix needed for
! ! chain closure.
! !
! ! INPUT:  rsf - index of final residue, allows location of s and u
! ! OUTPUT: jac_x - that jacobian value
! ! ----------------------------------------------------------------------
subroutine jacobian(rsf,jac_x)
!
  use polypep
  use atoms
  use math
  use sequen
  use iounit
!
 implicit none
!
  integer k,indx(5)
  RTYPE c1,c2,c3,c4,u1(3),u2(3),u3(3),u4(3),s1(3),s2(3),s3(3)
  RTYPE s4(3),s5(3),w1(3),w2(3),w3(3),t0(3),cp,cpo(3)
  RTYPE r1(3),r2(3),v1(3),v2(3),v3(3),v4(3),transform(5,5),det_a
  integer, INTENT (IN) :: rsf
  RTYPE, INTENT (OUT) :: jac_x
!
    if(rsf.gt.nseq) then
      write(ilog,*) 'Fatal. Called uj_conrot chain closure &
 &jacobian. Residue index out of bounds.'
      write(ilog,*) 'nseq=',nseq,'rsf=',rsf
      call fexit()
    end if
!
    s5(1) = x(cai(rsf))
    s5(2) = y(cai(rsf))
    s5(3) = z(cai(rsf))
    s4(1) = x(ni(rsf))
    s4(2) = y(ni(rsf))
    s4(3) = z(ni(rsf))
    s3(1) = x(ci(rsf-1))
    s3(2) = y(ci(rsf-1))
    s3(3) = z(ci(rsf-1))
    s2(1) = x(cai(rsf-1))
    s2(2) = y(cai(rsf-1))
    s2(3) = z(cai(rsf-1))
    s1(1) = x(ni(rsf-1))
    s1(2) = y(ni(rsf-1))
    s1(3) = z(ni(rsf-1))
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
end subroutine jacobian
! -----------------------------------------------------------------------------
!
! -----------------------------------------------------------------
subroutine uj_refvals(rsf,ccmat,oldvals)
! !------------------------------------------------------------
! ! subroutine to store values needed for chain closure
! ! INPUTS:
! !     rsf -  reference residue for chain closure
! ! MODIFIES:
! !     oldvals - {omega1,p1,p2,p3}
! !     ccmat - reference angle values
! !------------------------------------------------------------
!
  use iounit
  use polypep
  use atoms
  use math
  use fyoc
  use movesets
!
  implicit none
!
  integer rsf
  RTYPE ccmat(6),oldvals(4),getblen,getztor,getbang
!
! ! get old dihedrals
  ccmat(1) = getztor(ci(rsf-2),ni(rsf-1),cai(rsf-1),ci(rsf-1))
  ccmat(2) = getztor(ni(rsf-1),cai(rsf-1),ci(rsf-1),ni(rsf))
  ccmat(3) = getztor(ci(rsf-1),ni(rsf),cai(rsf),ci(rsf))
  oldvals(1) = getztor(cai(rsf-1),ci(rsf-1),ni(rsf),cai(rsf)) ! omega
!
! ! get old bond angles
  ccmat(4) = getbang(ni(rsf-1),cai(rsf-1),ci(rsf-1))
  ccmat(5) = getbang(cai(rsf-1),ci(rsf-1),ni(rsf))
  ccmat(6) = getbang(ci(rsf-1),ni(rsf),cai(rsf))
!
  oldvals(2) = getblen(ni(rsf-1),cai(rsf-1))
  oldvals(3) = getblen(cai(rsf-1),ci(rsf-1))
  oldvals(4) = getblen(ci(rsf-1),ni(rsf))
!
end
!
!---------------------------------------------------------------------
subroutine UJchainclose(rsf,ccmat,dccmat,oldvals,findsolu)
! !------------------------------------------------------------
! ! Subroutine to calculate the constraints for chain closure
! ! for UJ prerotation.
! ! INPUTS:
! !     rsf -  reference residue for chain closure
! !     oldvals - {w1old,omega1,p1,p2,p3}
! !     ccmat - reference angle values
! ! OUTPUT:
! !     dccmat - array of change in angles needed for closure
! !------------------------------------------------------------
!
  use iounit
  use polypep
  use atoms
  use math
  use fyoc
  use movesets
  use ujglobals
!
  implicit none
!
  integer j,k,rsf,iright,c1count,c2count,nstps
  logical findsolu,logi,searchflag1,searchflag2,logi2
  RTYPE w1,ccmat(6),dccmat(6),oldvals(4)
  RTYPE p1,p2,p3,p4,p5,nopl,u(3),v(3),con1(3),posCI_0(3)
  RTYPE q0(3),posCA(3),posNI(3),posNF(3),posCAF(3),posCIF(3)
  RTYPE next
  RTYPE rots(3,3),T3(3,3),rots2(3,3),epsilon,midpoint
  RTYPE q1(3,2),R3(3,3),R2(3,3),R1(3,3),T2(3,3),invT3(3,3)
  RTYPE invR1(3,3),invR2(3,3),invR3(3,3),invT1(3,3),invT2(3,3)
  RTYPE T1(3,3),gn,g1,intercpt,r(3),newccmat(2,6)
  RTYPE a3,ax,s(3),scaniv
  RTYPE con1xy(3),a1,ay,Ra1(3,3),Ra2(3,3),Ra3(3,3),az
  RTYPE newxaxis(3),newyaxis(3),newzaxis(3),yaxis(3),zaxis(3),a2
  RTYPE xaxis(3),RaE(3,3),or_pl(3),or_pl2(3),con5(3),con4(3)
  RTYPE con0(3),origin(3),s1(3),lcon5(3),w1_0,searchstp
!
! ! get the current position of atoms used in chain closure
!
  origin(1) = 0.0
  origin(2) = 0.0
  origin(3) = 0.0
  posCI_0(1) = x(ci(rsf-2))
  posCI_0(2) = y(ci(rsf-2))
  posCI_0(3) = z(ci(rsf-2))
  posNI(1) = x(ni(rsf-1))
  posNI(2) = y(ni(rsf-1))
  posNI(3) = z(ni(rsf-1))
  posCA(1) = x(cai(rsf-1))
  posCA(2) = y(cai(rsf-1))
  posCA(3) = z(cai(rsf-1))
  posNF(1) = x(ni(rsf))
  posNF(2) = y(ni(rsf))
  posNF(3) = z(ni(rsf))
  posCAF(1) = x(cai(rsf))
  posCAF(2) = y(cai(rsf))
  posCAF(3) = z(cai(rsf))
  posCIF(1) = x(ci(rsf))
  posCIF(2) = y(ci(rsf))
  posCIF(3) = z(ci(rsf)) 
!
  con0(:) = posCI_0(:)-posNI(:)
  con1(:) = posCA(:)-posNI(:)
  con4(:) = posCAF(:)-posNF(:)
  con5(:) = posCIF(:)-posCAF(:)
!
! ! get bond lengths N-CA, CA-C, and C-N
  p1 = oldvals(2)
  p2 = oldvals(3)
  p3 = oldvals(4)
  p4 = sqrt(sum((posCAF(:)-posNF(:))**2))
  p5 = sqrt(sum((posCIF(:)-posCAF(:))**2))
!
!  write(*,*) 'pre cc'
!  write(*,*) 'atom1={',posCI_0,'};'
!  write(*,*) 'atom2={',posNI,'};'
!  write(*,*) 'atom3={',posCA,'};'
!  write(*,*) 'atom5={',posNF,'};'
!  write(*,*) 'atom6={',posCAF,'};'
!  write(*,*) 'atom7={',posCIF,'};'
!  write(*,*) 'p1={',p1,',0,0}'
!  write(*,*) 'p2={',p2,',0,0}'
!  write(*,*) 'p3={',p3,',0,0}'
!  write(*,*) 'p4={',p4,',0,0}'
!  write(*,*) 'p5={',p5,',0,0}'
! ! build Euler rotation matrix
!
  xaxis = (/ 1, 0, 0 /)
  yaxis = (/ 0, 1, 0 /)
  zaxis = (/ 0, 0, 1 /)
!
! ! projection of con1 on xy-plane
  con1xy(1) = posCA(1) - posNI(1)
  con1xy(2) = posCA(2) - posNI(2)
  con1xy(3) = 0.0
!
! ! rotate around Z axis
  call bondang2(xaxis,origin,con1xy,a1)
! ! compare to y-axis to find if +/- rotation
  call bondang2(con1xy,origin,yaxis,ay)
  if(ay.gt.(PI/2.0)) then
    a1 = -a1
  end if
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
  call bondang2(con1,origin,newxaxis,a2)
  call bondang2(con1,origin,zaxis,az)
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
  call crossprod(newxaxis,con0,or_pl,nopl)
  call bondang2(newzaxis,origin,or_pl,a3)
  call crossprod(newxaxis,con0,or_pl2,nopl)
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
! ! create rotation matrix R3 and its inverse
  call genR((omega(rsf))/RADIAN,R3,invR3)
!
! ! set q0 = p1
  q0(1) = p1
  q0(2) = 0.0
  q0(3) = 0.0
!
! ! get pos S = RaE.posNF
  s1(1) = posNF(1)-posNI(1)
  s1(2) = posNF(2)-posNI(2)
  s1(3) = posNF(3)-posNI(3)
  s(1) = sum(RaE(1,:)*s1(:))
  s(2) = sum(RaE(2,:)*s1(:))
  s(3) = sum(RaE(3,:)*s1(:))
!  write(*,*) 's=',s
!
! ! find u
  u(1) = sum(RaE(1,:)*con4(:))
  u(2) = sum(RaE(2,:)*con4(:))
  u(3) = sum(RaE(3,:)*con4(:))
!  write(*,*) 'u=',u
!
! ! find v
  lcon5(1) = sum(RaE(1,:)*con5(:))
  lcon5(2) = sum(RaE(2,:)*con5(:))
  lcon5(3) = sum(RaE(3,:)*con5(:))
  call crossprod(lcon5,u,or_pl,nopl)
  call crossprod(u,or_pl,v,nopl)
!  write(*,*) 'v=',v
!
! ! step through w1 ranging from -20 degrees to +20 degrees of original
  findsolu = .false.
  searchflag1 = .true.
  searchflag2 = .true.
  c1count = 0
  c2count = 0
!
  scaniv = UJ_params(5)
  searchstp = min((scaniv/20.0),UJ_params(4))
  w1_0 = (ccmat(1)-scaniv-searchstp)/radian
  nstps = 2.0*ceiling(scaniv/searchstp)
! PBC
  if (w1_0.lt.-PI) w1_0 = w1_0 + 2.0*PI
!
90    do j=1,nstps
    w1 = w1_0 + searchstp*j/radian
!   PBC
    if (w1.gt.PI) w1 = w1 - 2.0*PI
    next = w1 +searchstp/radian
!   PBC
    if (next.gt.PI) next = next - 2.0*PI
    call checkreal(w1,logi)
    call checkreal(next,logi2)
    if ((logi2.EQV..true.).AND.(logi.EQV..true.)) then
      if(searchflag1) then
        c1count = c1count +1
        call g_of_w1(next,gn,1)
        call g_of_w1(w1,g1,1)
        intercpt = gn*g1
        if(intercpt.lt.0.0d0) then 
          iright = 1
          call bisection(w1,next,1)
!          write(ilog,*) 'Bisection Method::found a solution 1'
!          write(*,*) 'w11 = ',RADIAN*w1
          findsolu = .true.
          searchflag1 = .false.
        end if
      end if
      w1 = w1_0 + searchstp*j/radian
      if (w1.gt.PI) w1 = w1 - 2.0*PI
      next = w1 + searchstp/radian
      if (next.gt.PI) next = next - 2.0*PI
      if(searchflag2) then
        c2count = c2count +1
        call g_of_w1(next,gn,2)
        call g_of_w1(w1,g1,2)
        intercpt = gn*g1
        if(intercpt.lt.0.0d0) then
          iright = 2
          call bisection(w1,next,2)
!          write(ilog,*) 'Bisection Method::found a solution 2'
!          write(*,*) 'w12 = ',RADIAN*w1
          findsolu = .true.
          searchflag2 = .false.
        end if
      end if
      if((.NOT.searchflag1).AND.(.NOT.searchflag2)) then
!        write(*,*)'both solutions found'
        iright = 3
        EXIT
      end if
    end if
  end do
!
  if(.NOT.findsolu) then    
!    write(ilog,*) 'Did not find a solution to the chain closure'
  else
    call get_newccmat(newccmat,iright)
    if (.NOT.findsolu) then
!     do nothing
    else
!      write(*,*) 'new values 1 = ', newccmat(1,:)
!      write(*,*) 'new values 2 = ', newccmat(2,:)
      do k=1,6
        dccmat(k) = RADIAN*newccmat(iright,k)-ccmat(k)
      end do
    end if
  end if
!
!  write(*,*) 'counts = ',c1count,c2count
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
  implicit none
!
  RTYPE w1val, arg1, arg2
  logical isreal
!
  isreal = .false.
  call r_of_w1(r,w1val)
!
  arg1 = 4*(p2**2)*(r(1)**2+r(2)**2)
  arg2 = (r(1)**2 + r(2)**2 + r(3)**2 + p2**2 - p3**2)**2
!
  if(arg1.gt.arg2) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
end subroutine checkreal
! -------------------------------------------------------------------  
!
! ! bisection(w1)
subroutine bisection(w1val,nextw1,t)
!
  implicit none
!
  integer t
  RTYPE w1val,nextw1,w1tmp,grr,grr2
!
  epsilon = 2.0*(10.0**(-15.0))
! ! start loop
  do while (abs(nextw1 - w1val).gt.(2.0*epsilon))
!
! ! calculate midpoint of domain
    if ((nextw1 - w1val).gt.PI) then
      w1tmp = w1val + 2.0*PI
      midpoint = (nextw1 + w1tmp) / 2.0
    else if ((nextw1 - w1val).lt.-PI) then
      w1tmp = w1val - 2.0*PI
      midpoint = (nextw1 + w1tmp) / 2.0
    else
      midpoint = (nextw1 + w1val) / 2.0
    end if
    if (midpoint.ge.PI) midpoint = midpoint - 2.0*PI
    if (midpoint.lt.-PI) midpoint = midpoint + 2.0*PI
!
!   ! find f(midpoint)
    call g_of_w1(w1val,grr,t)
    call g_of_w1(midpoint,grr2,t)
    intercpt = grr*grr2
    if (intercpt.gt. 0.0) then
!     ! throw away left half
      w1val = midpoint
    else
!     ! throw away right half
      nextw1 = midpoint
    end if 
  
  end do
! ! return final midpoint value
  w1val = (nextw1 + w1val) / 2.0
!
  newccmat(t,1) = w1val
!  write(*,*) 'in bis:',t,RADIAN*newccmat(t,1)
!
end subroutine bisection
! -------------------------------------------------------------------  
!
! ! subroutine to get r(w1) = inv(R1)*(s-q0)
subroutine r_of_w1(r,w1val)
!
  implicit none
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
! ! subroutine to find alpha1
subroutine alpha1_of_w1(newccmat,t)
!
  implicit none
!
  integer t
  RTYPE newccmat(2,6),w,arg1,arg2
!
  call r_of_w1(r,newccmat(t,1))
  w = r(1)**2 + r(2)**2 + r(3)**2 + p2**2 - p3**2
!     
  if (t.eq.1) then
    arg1=(w*r(1) + r(2)*sqrt(4*p2**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p2*(r(1)**2 + r(2)**2))
    arg2=(w*r(2) - r(1)*sqrt(4*p2**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p2*(r(1)**2 + r(2)**2))
  else if(t.eq.2) then
    arg1=(w*r(1) - r(2)*sqrt(4*p2**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p2*(r(1)**2 + r(2)**2))
    arg2=(w*r(2) + r(1)*sqrt(4*p2**2 * (r(1)**2 + r(2)**2) - w**2))/&
 &(2*p2*(r(1)**2 + r(2)**2))
  else
    write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
 &ain closure algorithm) must be either 1 or 2.'
    call fexit()
  end if
!
  newccmat(t,4) = atan2(arg2,arg1)
!  write(*,*) 'a1',t,arg1,arg2
!
end subroutine alpha1_of_w1
! -------------------------------------------------------------------  
!
! ! subroutine to find alpha2 and w2
subroutine alpha2_w2_of_w1(newccmat,t)
!
  implicit none
!
  integer t
  RTYPE newccmat(2,6),arg1,arg2
!
  call alpha1_of_w1(newccmat,t)
!
  q1(1,t) = cos(newccmat(t,4))*r(1) + sin(newccmat(t,4))*r(2)
  q1(2,t) = -sin(newccmat(t,4))*r(1) + cos(newccmat(t,4))*r(2)
  q1(3,t) = r(3)
!
  newccmat(t,5) = acos((q1(1,t)-p2)/p3)
!  write(*,*) 'a2',t,(q1(1,t)-p2)/p3
!
  arg1 = q1(2,t)/sin(newccmat(t,5))/p3
  arg2 = q1(3,t)/sin(newccmat(t,5))/p3
  newccmat(t,2) = atan2(arg2,arg1)
!  write(*,*) 'w2',t,arg1,arg2
!
end subroutine alpha2_w2_of_w1
! -------------------------------------------------------------------  
!
! ! subroutine to find alpha3 and w4
subroutine get_newccmat(newccmat,iright)
!
  implicit none
!
  integer iright,mode
  logical validsol(2)
  RTYPE dmat(3),dw,dalpha,newccmat(2,6),arg1,arg2,random
!
! ! choose the solution that is within +/-20 alpha and +/-50 w
  dw = 50.0*PI/180.0
  dalpha = 20.0*PI/180.0
  mode = iright
  iright = 0
  validsol(:) = .false.
!  write(*,*) 'in get_newccmat'
!
!  write(*,*) 'mode = ',mode
  do k=1,2
    if ((k.eq.mode).OR.(mode.eq.3)) then
      call alpha2_w2_of_w1(newccmat,k)
!
      call genR(newccmat(k,2),R2,invR2)
      call genR(newccmat(k,1),R1,invR1)
      call genT(newccmat(k,5),T2,invT2)
      call genT(newccmat(k,4),T1,invT1)
!
!     ! inv(R3)*inv(T2)*inv(R2)*inv(T1)*inv(R1)
      rots = matmul(invR3,matmul(invT2,matmul(invR2,matmul(invT1,inv&
 &R1))))
!     ! find alpha3
      arg1 = rots(1,1)*u(1) + rots(1,2)*u(2) + rots(1,3)*u(3)
      arg2 = rots(2,1)*u(1) + rots(2,2)*u(2) + rots(2,3)*u(3)
!      write(*,*) 'utest: ',rots(3,1)*u(1) + rots(3,2)*u(2) + 
! &rots(3,3)*u(3)
      newccmat(k,6) = atan2(arg2,arg1) 
!
!     ! now w4 using alpha3
      call genT(newccmat(k,6),T3,invT3)
      rots2 = matmul(invT3,rots)
!      write(*,*) 'vtest: ',rots2(1,1)*v(1) + rots2(1,2)*v(2) +
! &                         rots2(1,3)*v(3)
      arg1 = rots2(2,1)*v(1) + rots2(2,2)*v(2) + rots2(2,3)*v(3)
      arg2 = rots2(3,1)*v(1) + rots2(3,2)*v(2) + rots2(3,3)*v(3)
      newccmat(k,3) = atan2(arg2,arg1)
!
!      xx(:) = q0(:) + Matmul(Matmul(R1,T1),q1(:,k))
!      write(*,*) 'stest1: ',xx(1),s(1)
!      write(*,*) 'stest2: ',xx(2),s(2)
!      write(*,*) 'stest3: ',xx(3),s(3)
!      write(*,*) RADIAN*newccmat(k,3)
  if (newccmat(k,4).lt.PI) newccmat(k,4) = PI - abs(newccmat(k,4))
  if (newccmat(k,5).lt.PI) newccmat(k,5) = PI - abs(newccmat(k,5))
  if (newccmat(k,6).lt.PI) newccmat(k,6) = PI - abs(newccmat(k,6))
!  newccmat(k,3) = newccmat(k,3) + PI
!  if (newccmat(k,3).gt.PI) newccmat(k,3) = newccmat(k,3) - 2.0*PI
!
!      write(*,*) 'newccmat = ',newccmat(k,:)
!
      dmat(1) = abs(ccmat(2)/RADIAN-newccmat(k,2))
      dmat(2) = abs(ccmat(4)/RADIAN-newccmat(k,4))
      dmat(3) = abs(ccmat(5)/RADIAN-newccmat(k,5))
!      write(*,*) 'dmat = ',dmat
      if ((dmat(1).lt.dw).AND.(dmat(2).lt.dalpha)&
 &   .AND.(dmat(3).lt.dalpha)) then
        validsol(k) = .true.
!        write(ilog,*) 'valid solution in branch ',iright
        findsolu = .true.
      end if
    end if  
  end do
!  write(*,*) 'ccmat   : ',ccmat

!  if ((iright.eq.1).OR.(iright.eq.3)) then
!    write(*,*) 'newccmat: ',RADIAN*newccmat(1,:)
!  end if
!  if ((iright.eq.2).OR.(iright.eq.3)) then
!    write(*,*) 'newccmat: ',RADIAN*newccmat(2,:)
!  end if
! ! check that at least one valid solution was found
  if ((validsol(1).EQV..false.).AND.(validsol(2).EQV..false.)) then
!    write(ilog,*) 'Warning: neither branch gave valid solution'
    findsolu = .false.
  else if ((validsol(1).EQV..true.).AND.(validsol(2).EQV..true.))&
 &then
    iright = 1
    if (random().gt.0.5) iright = 2
!    write(*,*) 'picked random: ',iright
  else if (validsol(1).EQV..true.) then
    iright = 1
  else
    iright = 2
  end if
!
end subroutine get_newccmat
! -------------------------------------------------------------------  
!
! ! subroutine to generate T(i) and its inverse (invTi)
subroutine genT(alpha,Ti,invTi)
!
  implicit none
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
  implicit none
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
!
! -------------------------------------------------------------------  
!
! ! subroutine outputs g(w1)
!
subroutine g_of_w1(w1val,g,t)
!
  implicit none
!
  integer t
  RTYPE w1val,g
  RTYPE invTRmat(3,3)
!
! ! assign w1
  newccmat(t,1) = w1val
!
! ! using branch t get solution for g(w1), find 3 dofs needed
  call alpha2_w2_of_w1(newccmat,t)
!
! ! calculate g(w1)
  call genR(newccmat(t,1),R1,invR1)
  call genR(newccmat(t,2),R2,invR2)
  call genT(newccmat(t,4),T1,invT1)
  call genT(newccmat(t,5),T2,invT2)
  invTRmat = matmul(invR3,matmul(invT2,matmul(invR2,matmul(invT1,inv&
 &R1))))
  g = sum(invTRmat(3,:)*u(:))
!
end subroutine g_of_w1
!
! ------------------------------------------------------------------
end subroutine UJchainclose
