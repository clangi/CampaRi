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
subroutine buildnuccrdof(rsi,rsf,dof,nuccrmat)
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
  use molecule
!
  implicit none
!
  integer j,k,m,rsi,rsf,dof,rs,dofcnt,nucpos,shf
  RTYPE pos2(3),pos3(3),pos4(3)
  RTYPE nuccrmat(MAXUJDOF,MAXUJDOF),doflever(3,MAXUJDOF),dvec(3)
!
! dof has total DOF involved in current prerotation
! loop over all free degrees of freedom of the prerotation

  dofcnt = 0
  nucpos = 6
  do rs=rsf-1,rsi,-1
    shf = 0
    if (rs.eq.(rsmol(molofrs(rsi),1)+1)) then
      if (seqflag(rsmol(molofrs(rsi),1)).eq.24) then
        shf = 2
      end if
    end if
    pos4(1) = x(nuci(rs,2))
    pos4(2) = y(nuci(rs,2))
    pos4(3) = z(nuci(rs,2))
    pos3(1) = x(nuci(rs,3))
    pos3(2) = y(nuci(rs,3))
    pos3(3) = z(nuci(rs,3))
    pos2(1) = x(nuci(rs,4))
    pos2(2) = y(nuci(rs,4))
    pos2(3) = z(nuci(rs,4))
    call get_tor_lever(molofrs(atmres(nuci(rsf,5))),nuci(rsf,5),pos2,pos3,pos4,dvec)
    doflever(:,dof-dofcnt) = dvec(:)/RADIAN
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) exit
    pos4(1) = x(nuci(rs,1))
    pos4(2) = y(nuci(rs,1))
    pos4(3) = z(nuci(rs,1))
    pos3(1) = x(nuci(rs,2))
    pos3(2) = y(nuci(rs,2))
    pos3(3) = z(nuci(rs,2))
    pos2(1) = x(nuci(rs,3))
    pos2(2) = y(nuci(rs,3))
    pos2(3) = z(nuci(rs,3))
    call get_tor_lever(molofrs(atmres(nuci(rsf,5))),nuci(rsf,5),pos2,pos3,pos4,dvec)
    doflever(:,dof-dofcnt) = dvec(:)/RADIAN
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) exit

!    if (rs.lt.rsf-1) then
      pos4(1) = x(nuci(rs-1,nucpos-shf))
      pos4(2) = y(nuci(rs-1,nucpos-shf))
      pos4(3) = z(nuci(rs-1,nucpos-shf))
      pos3(1) = x(nuci(rs,1))
      pos3(2) = y(nuci(rs,1))
      pos3(3) = z(nuci(rs,1))
      pos2(1) = x(nuci(rs,2))
      pos2(2) = y(nuci(rs,2))
      pos2(3) = z(nuci(rs,2))
      call get_tor_lever(molofrs(atmres(nuci(rsf,5))),nuci(rsf,5),pos2,pos3,pos4,dvec)
      doflever(:,dof-dofcnt) = dvec(:)/RADIAN
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) exit
      pos4(1) = x(nuci(rs-1,nucpos-shf-1))
      pos4(2) = y(nuci(rs-1,nucpos-shf-1))
      pos4(3) = z(nuci(rs-1,nucpos-shf-1))
      pos3(1) = x(nuci(rs-1,nucpos-shf))
      pos3(2) = y(nuci(rs-1,nucpos-shf))
      pos3(3) = z(nuci(rs-1,nucpos-shf))
      pos2(1) = x(nuci(rs,1))
      pos2(2) = y(nuci(rs,1))
      pos2(3) = z(nuci(rs,1))
      call get_tor_lever(molofrs(atmres(nuci(rsf,5))),nuci(rsf,5),pos2,pos3,pos4,dvec)
      doflever(:,dof-dofcnt) = dvec(:)/RADIAN
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) exit
      pos4(1) = x(nuci(rs-1,nucpos-shf-3))
      pos4(2) = y(nuci(rs-1,nucpos-shf-3))
      pos4(3) = z(nuci(rs-1,nucpos-shf-3))
      pos3(1) = x(nuci(rs-1,nucpos-shf-2))
      pos3(2) = y(nuci(rs-1,nucpos-shf-2))
      pos3(3) = z(nuci(rs-1,nucpos-shf-2))
      pos2(1) = x(nuci(rs-1,nucpos-shf-1))
      pos2(2) = y(nuci(rs-1,nucpos-shf-1))
      pos2(3) = z(nuci(rs-1,nucpos-shf-1))
      call get_tor_lever(molofrs(atmres(nuci(rsf,5))),nuci(rsf,5),pos2,pos3,pos4,dvec)
      doflever(:,dof-dofcnt) = dvec(:)/RADIAN
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) exit
!    end if

  end do
!
! ! build symetric matrix I = dA/dphi*dA/dphi (outer_product)
  do j=1,dof
    do k=1,dof
!     ! inner product
      nuccrmat(j,k) = 0.0d0
      do m=1,3
        nuccrmat(j,k) = nuccrmat(j,k) + doflever(m,j)*doflever(m,k)
      end do
    end do
  end do
!
end subroutine buildnuccrdof
!
!-----------------------------------------------------------------------
! Subroutine that finds the jacobian factor of the 5x5 matrix needed for
! chain closure.
!
! INPUT:  rsf - index of final residue, allows location of s and u
! OUTPUT: jac_x - that jacobian value
! ----------------------------------------------------------------------
subroutine jacobian_nuccr(rsf,jac_x)
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
  RTYPE c6,c1,c3,c4,c5,c7,u6(3),u1(3),u3(3),u4(3),u5(3),s1(3),s2(3),u7(3)
  RTYPE s3(3),s4(3),s5(3),t0(3),cp,r1(3),r3(3),r4(3),r5(3),s0(3)
  RTYPE v1(3),v3(3),v4(3),v5(3),transform(5,5),det_a,s7(3),s6(3)
  integer, INTENT (IN) :: rsf
  RTYPE, INTENT (OUT) :: jac_x
!
    if(rsf.gt.nseq) then
      write(ilog,*) 'Fatal. Called torcr_conrot chain closure &
 &jacobian. Residue index out of bounds.'
      write(ilog,*) 'nseq=',nseq,'rsf=',rsf
      call fexit()
    end if
!
    s7(1) = x(nuci(rsf,5))
    s7(2) = y(nuci(rsf,5))
    s7(3) = z(nuci(rsf,5))
    s6(1) = x(nuci(rsf,4))
    s6(2) = y(nuci(rsf,4))
    s6(3) = z(nuci(rsf,4))
    s5(1) = x(nuci(rsf,3))
    s5(2) = y(nuci(rsf,3))
    s5(3) = z(nuci(rsf,3))
    s4(1) = x(nuci(rsf,2))
    s4(2) = y(nuci(rsf,2))
    s4(3) = z(nuci(rsf,2))
    s3(1) = x(nuci(rsf,1))
    s3(2) = y(nuci(rsf,1))
    s3(3) = z(nuci(rsf,1))
    s2(1) = x(nuci(rsf-1,6))
    s2(2) = y(nuci(rsf-1,6))
    s2(3) = z(nuci(rsf-1,6))
    s1(1) = x(nuci(rsf-1,5))
    s1(2) = y(nuci(rsf-1,5))
    s1(3) = z(nuci(rsf-1,5))
    s0(1) = x(nuci(rsf-1,4))
    s0(2) = y(nuci(rsf-1,4))
    s0(3) = z(nuci(rsf-1,4))
!
    do k=1,3
      u1(k) = s1(k) - s0(k)
      u3(k) = s3(k) - s2(k)
      u4(k) = s4(k) - s3(k)
      u5(k) = s5(k) - s4(k)  
      u6(k) = s6(k) - s5(k)
      u7(k) = s7(k) - s6(k)
    end do
!   
    c1 = 0.0
    c3 = 0.0
    c4 = 0.0
    c5 = 0.0
    c6 = 0.0
    c7 = 0.0
    do k=1,3
      c1 = c1 + u1(k)*u1(k)
      c3 = c3 + u3(k)*u3(k)
      c4 = c4 + u4(k)*u4(k)
      c5 = c5 + u5(k)*u5(k)
      c6 = c6 + u6(k)*u6(k)
      c7 = c7 + u7(k)*u7(k)
    end do
    u1(:) = u1(:)/sqrt(c1)
    u3(:) = u3(:)/sqrt(c3)
    u4(:) = u4(:)/sqrt(c4)
    u5(:) = u5(:)/sqrt(c5)
    u6(:) = u6(:)/sqrt(c6)
    u7(:) = u7(:)/sqrt(c7)
    r1(:) = s6(:) - s1(:)
    r3(:) = s6(:) - s3(:)
    r4(:) = s6(:) - s4(:)
    r5(:) = s6(:) - s5(:)
    call crossprod(u1,r1,v1,cp)
    call crossprod(u3,r3,v3,cp)
    call crossprod(u4,r4,v4,cp)
    call crossprod(u5,r5,v5,cp)
    do k=1,3
      transform(k,1) = v1(k)
      transform(k,2) = v3(k)
      transform(k,3) = v4(k)
      transform(k,4) = v5(k)
      transform(k,5) = 0.0
    end do
!
    call crossprod(u1,u7,t0,cp)
    transform(4,1) = t0(1)
    transform(5,1) = t0(2)
    call crossprod(u3,u7,t0,cp)
    transform(4,2) = t0(1)
    transform(5,2) = t0(2)
    call crossprod(u4,u7,t0,cp)
    transform(4,3) = t0(1)
    transform(5,3) = t0(2)
    call crossprod(u5,u7,t0,cp)
    transform(4,4) = t0(1)
    transform(5,4) = t0(2)
    call crossprod(u6,u7,t0,cp)
    transform(4,5) = t0(1)
    transform(5,5) = t0(2)
!
    call dtrm(transform,5,det_a,indx)
!
    jac_x = 1.0/abs(det_a)
!
end subroutine jacobian_nuccr
!
! -----------------------------------------------------------------
subroutine nuccr_refvals(rsf,ccmat_nuc,oldvals_nuc)
!------------------------------------------------------------
! subroutine to store values needed for chain closure
! INPUTS:
!     rsf -  reference residue for chain closure
! MODIFIES:
!     oldvals_nuc - {p1,p2,p3,p4,p5,p6,a1,a2,a3,a4,a5,a6}
!     ccmat_nuc - {w1,w2,w3,w4,w5,w6}
!------------------------------------------------------------
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
  RTYPE ccmat_nuc(7),oldvals_nuc(14),getblen,getztor,getbang
!
! get old bond lengths
  oldvals_nuc(1) = getblen(nuci(rsf-1,4),nuci(rsf-1,5))
  oldvals_nuc(2) = getblen(nuci(rsf-1,5),nuci(rsf-1,6))
  oldvals_nuc(3) = getblen(nuci(rsf-1,6),nuci(rsf,1))
  oldvals_nuc(4) = getblen(nuci(rsf,1),nuci(rsf,2))
  oldvals_nuc(5) = getblen(nuci(rsf,2),nuci(rsf,3))
  oldvals_nuc(6) = getblen(nuci(rsf,3),nuci(rsf,4))
  oldvals_nuc(7) = getblen(nuci(rsf,4),nuci(rsf,5))
!
! get old bond angles
  oldvals_nuc(8) =  PI - getbang(nuci(rsf-1,4),nuci(rsf-1,5),nuci(rsf-1,6))/RADIAN
  oldvals_nuc(9) =  PI - getbang(nuci(rsf-1,5),nuci(rsf-1,6),nuci(rsf,1))/RADIAN
  oldvals_nuc(10) = PI - getbang(nuci(rsf-1,6),nuci(rsf,1),nuci(rsf,2))/RADIAN
  oldvals_nuc(11) = PI - getbang(nuci(rsf,1),nuci(rsf,2),nuci(rsf,3))/RADIAN
  oldvals_nuc(12) = PI - getbang(nuci(rsf,2),nuci(rsf,3),nuci(rsf,4))/RADIAN
  oldvals_nuc(13) = PI - getbang(nuci(rsf,3),nuci(rsf,4),nuci(rsf,5))/RADIAN
  oldvals_nuc(14) = PI - getbang(nuci(rsf,4),nuci(rsf,5),nuci(rsf,6))/RADIAN
!
! get old torsions
  ccmat_nuc(1) = getztor(nuci(rsf-1,3),nuci(rsf-1,4),nuci(rsf-1,5),nuci(rsf-1,6))/RADIAN
  ccmat_nuc(7) = getztor(nuci(rsf-1,4),nuci(rsf-1,5),nuci(rsf-1,6),nuci(rsf,1))/RADIAN
  ccmat_nuc(2) = getztor(nuci(rsf-1,5),nuci(rsf-1,6),nuci(rsf,1),nuci(rsf,2))/RADIAN
  ccmat_nuc(3) = getztor(nuci(rsf-1,6),nuci(rsf,1),nuci(rsf,2),nuci(rsf,3))/RADIAN
  ccmat_nuc(4) = getztor(nuci(rsf,1),nuci(rsf,2),nuci(rsf,3),nuci(rsf,4))/RADIAN
  ccmat_nuc(5) = getztor(nuci(rsf,2),nuci(rsf,3),nuci(rsf,4),nuci(rsf,5))/RADIAN
  ccmat_nuc(6) = getztor(nuci(rsf,3),nuci(rsf,4),nuci(rsf,5),nuci(rsf,6))/RADIAN
!
end
!
!---------------------------------------------------------------------
subroutine picksolu_nc(rsi,rsf,dof,dphi,solmat,solucount,solcnt2,jac_b_vec,dpfmat,j,prob_move,prob_back)
!------------------------------------------------------------
!
  use iounit
  use ujglobals
!
  implicit none
!
  integer rsi,rsf,dof,solucount,afive,aone,j,solcnt2
  RTYPE dpfmat(MAXUJDOF+6),dphi(MAXUJDOF+6),random,jac_b_vec(4*MAXSOLU)
  RTYPE solmat(MAXSOLU*4,6),dummy(2),prob_move,prob_back
!
  afive = 5
  aone = 1
!
  if (torcrmode.eq.2) then
    j = random()*solucount + 1
    dpfmat(1:dof) = dphi(1:dof)
    dpfmat(dof+1:dof+6) = solmat(j,1:6)
!   restore, then move entire internal chain
    call mcmove_nuccr(rsi,rsf,dof,dpfmat,afive)
    call mcmove_nuccr(rsi,rsf,dof,dpfmat,aone)
    return
  end if
!
! we get into trouble if prerotation bias is numerically tiny but other
! set is empty (happens for large dof and large bias due to exp())
  if ((solucount.eq.0).AND.(prob_back.le.0.0)) then
    prob_back = 1.0
  else if ((solcnt2.eq.0).AND.(prob_move.le.0.0)) then
    prob_move = 1.0
  end if
  dummy(1:2) = 0.0
  if (solucount.gt.0) then
    dummy(1) = sum(prob_move*jac_b_vec(1:solucount))
  end if
  if (solcnt2.gt.0) then
    dummy(2) = sum(prob_back*jac_b_vec(solucount+1:solucount+solcnt2))
  end if
  dummy(1) = random()*(dummy(1)+dummy(2))
  if (1.gt.solucount) then
    dummy(2) = jac_b_vec(1)*prob_back
  else
    dummy(2) = jac_b_vec(1)*prob_move
  end if
  do j=1,solucount+solcnt2
    if (dummy(1).lt.dummy(2)) then
      if (j.gt.solucount) then
        dpfmat(1:dof) = 0.0
      else
        dpfmat(1:dof) = dphi(1:dof)
      end if
      dpfmat(dof+1:dof+6) = solmat(j,1:6)
!     restore, then move entire internal chain
      call mcmove_nuccr(rsi,rsf,dof,dpfmat,afive)
      call mcmove_nuccr(rsi,rsf,dof,dpfmat,aone)
      exit
    else
      if (j.eq.(solucount+solcnt2)) then
        write(ilog,*) 'Fatal. Jacobian weighting of solutions is inconsistent. This is a bug.'
        call fexit()
      end if
      if ((j+1).gt.solucount) then
        dummy(2) = dummy(2) + jac_b_vec(j+1)*prob_back
      else
        dummy(2) = dummy(2) + jac_b_vec(j+1)*prob_move
      end if
    end if
  end do
!
  end
!
!---------------------------------------------------------------------
subroutine nuccrchainclose(rsi,rsf,dof,dphi,ccmat_nuc,incrmat,oldvals_nuc,solucount,jac_b_vec)
!------------------------------------------------------------
! Subroutine to calculate the constraints for chain closure
! for UJ prerotation.
! INPUTS:
!     rsi -  reference residue for start of perturbation
!     rsf -  reference residue for chain closure
!     dof -  total number of pre-rotation dof.s
!     dphi - pre-rotation increments
!     oldvals_nuc - reference bond lengths and angles
!     ccmat_nuc - reference torsion angle values
! OUTPUT:
!     dpfmat - array of change in angles needed for entire perturbation
!     solucount - number of solutions found for closure
!     jac_b - Jacobian of selected closure solution
!     
!------------------------------------------------------------
!
  use iounit
  use polypep
  use atoms
  use math
  use fyoc
  use movesets
  use ujglobals
  use mcsums
!
  implicit none
!
  integer j,k,kk,kkk,rsi,rsf,solucount,npreintvs,jj,jjj
  integer findsolu(4),dof,athree,atwo,aone,afive,nintvs(4)
  logical logi1,logi1b,logi1a,goodsol,really,really2,really3,atrue!,slog(3)
  RTYPE w,w1,ccmat_nuc(7),dpfmat(MAXUJDOF+6),oldvals_nuc(14),dphi(MAXUJDOF+6)
  RTYPE p1,p2,p3,p4,p5,p6,p7,nopl,u(3),con1(3),posO5(3),gnpo,gnpr,maxincr
  RTYPE q0(3),posC4(3),posC5(3),posC3F(3),posC4F(3),posC5F(3),intvalbnds(4,MAXSOLU*4,2)
  RTYPE q1(3),q2(3),R3(3,3),R2(3,3),R1(3,3),T2(3,3),invT3(3,3)
  RTYPE invR1(3,3),invR2(3,3),invR3(3,3),invT1(3,3),invT2(3,3),T5(3,3)
  RTYPE R4(3,3),invR4(3,3),T4(3,3),invT4(3,3),R5(3,3),invR5(3,3),invT5(3,3)
  RTYPE T1(3,3),gn,g1,intercpt,newccmat2(4,6),incrmat(MAXSOLU*4,6)
  RTYPE a3,ax,s(3),T3(3,3),epsilon,midpoint,w1mat(4,MAXSOLU)
  RTYPE con1xy(3),a1,ay,Ra1(3,3),Ra2(3,3),Ra3(3,3),az,increment,curslope(4)
  RTYPE newxaxis(3),newyaxis(3),newzaxis(3),yaxis(3),zaxis(3),a2,ltestw1,htestw1,lastslope(4)
  RTYPE xaxis(3),RaE(3,3),or_pl(3),or_pl2(3),preintvbs(MAXSOLU*4,2)
  RTYPE con0(3),origin(3),s1(3),searchstp,con7(3),oldws(4),oldgs(4)
  RTYPE invmat(3,3),u1,invRaE(3,3),jac_b_vec(MAXSOLU*4),newt
  RTYPE sinov(14),cosov(14),sinccn7,cosccn7,q1sq,q2sq,nouse(20),w5denom
!
! helpers
  aone = 1
  atwo = 2
  athree = 3
  afive = 5
  atrue = .true.
!
! get the current position of atoms used in chain closure
  posO5(1) = x(nuci(rsf-1,3))
  posO5(2) = y(nuci(rsf-1,3))
  posO5(3) = z(nuci(rsf-1,3))
  posC5(1) = x(nuci(rsf-1,4))
  posC5(2) = y(nuci(rsf-1,4))
  posC5(3) = z(nuci(rsf-1,4))
  posC4(1) = x(nuci(rsf-1,5))
  posC4(2) = y(nuci(rsf-1,5))
  posC4(3) = z(nuci(rsf-1,5))
  posC5F(1) = x(nuci(rsf,4))
  posC5F(2) = y(nuci(rsf,4))
  posC5F(3) = z(nuci(rsf,4)) 
  posC4F(1) = x(nuci(rsf,5))
  posC4F(2) = y(nuci(rsf,5))
  posC4F(3) = z(nuci(rsf,5))
  posC3F(1) = x(nuci(rsf,6))
  posC3F(2) = y(nuci(rsf,6))
  posC3F(3) = z(nuci(rsf,6))
!
  con0(:) = posO5(:)-posC5(:)
  con1(:) = posC4(:)-posC5(:)
  con7(:) = posC4F(:)-posC5F(:)
!
! get bond lengths
  p1 = oldvals_nuc(1)
  p2 = oldvals_nuc(2)
  p3 = oldvals_nuc(3)
  p4 = oldvals_nuc(4)
  p5 = oldvals_nuc(5)
  p6 = oldvals_nuc(6)
  p7 = oldvals_nuc(7)
!
! trigonometric terms
  sinov(8) = sin(oldvals_nuc(8))
  cosov(8) = cos(oldvals_nuc(8))
  sinov(9) = sin(oldvals_nuc(9))
  cosov(9) = cos(oldvals_nuc(9))
  sinov(10) = sin(oldvals_nuc(10))
  cosov(10) = cos(oldvals_nuc(10))
  sinov(11) = sin(oldvals_nuc(11))
  cosov(11) = cos(oldvals_nuc(11))
  sinov(12) = sin(oldvals_nuc(12))
  cosov(12) = cos(oldvals_nuc(12))
  sinov(13) = sin(oldvals_nuc(13))
  cosov(13) = cos(oldvals_nuc(13))
  sinccn7 = sin(ccmat_nuc(7))
  cosccn7 = cos(ccmat_nuc(7))
!
! static matrices
  call genR(ccmat_nuc(7),R2,invR2)
  call genT(oldvals_nuc(8),T1,invT1)
  call genT(oldvals_nuc(9),T2,invT2)
  call genT(oldvals_nuc(10),T3,invT3)
  call genT(oldvals_nuc(11),T4,invT4)
  call genT(oldvals_nuc(12),T5,invT5)
!
! build Euler rotation matrix
  origin(1) = 0.0
  origin(2) = 0.0
  origin(3) = 0.0
  xaxis = (/ 1, 0, 0 /)
  yaxis = (/ 0, 1, 0 /)
  zaxis = (/ 0, 0, 1 /)

! projection of con1 on xy-plane
  con1xy(1) = posC4(1) - posC5(1)
  con1xy(2) = posC4(2) - posC5(2)
  con1xy(3) = 0.0

! rotate around Z axis
  call bondang2(xaxis,origin,con1xy,a1)
! compare to y-axis to find if +/- rotation
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

! rotate around new Y axis (angle between new x-axis and con1)
  newxaxis(1) = Ra1(1,1)
  newxaxis(2) = Ra1(1,2)
  newxaxis(3) = Ra1(1,3)
  call bondang2(con1,origin,newxaxis,a2)
  call bondang2(con1,origin,zaxis,az)
  if(az.lt.(PI/2)) then
    a2 = -a2
  end if

! rotation matrix for y-rotation
  Ra2(1,1) = cos(a2)
  Ra2(1,2) = 0.0
  Ra2(1,3) = -sin(a2)
  Ra2(2,1) = 0.0
  Ra2(2,2) = 1.0
  Ra2(2,3) = 0.0
  Ra2(3,1) = -Ra2(1,3)
  Ra2(3,2) = 0.0
  Ra2(3,3) = Ra2(1,1)

! rotate around new X axis so that Y'' and p1 align
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

! build full euler rotation matrix
  RaE = Matmul(Ra3,Matmul(Ra2,Ra1))
  invRaE(1,1) = RaE(1,1)
  invRaE(1,2) = RaE(2,1)
  invRaE(1,3) = RaE(3,1)
  invRaE(2,1) = RaE(1,2)
  invRaE(2,2) = RaE(2,2)
  invRaE(2,3) = RaE(3,2)
  invRaE(3,1) = RaE(1,3)
  invRaE(3,2) = RaE(2,3)
  invRaE(3,3) = RaE(3,3)

! get q1 an q2
  q1(1) = cosov(10)*p4 + p3
  q1(2) = sinov(10)*p4
  q1(3) = 0.0
  q1sq = sum(q1(:)*q1(:))
  q2(1) = cosov(12)*p6 + p5
  q2(2) = sinov(12)*p6
  q2(3) = 0.0
  q2sq = sum(q2(:)*q2(:))
  w5denom = (q2(2)*q2(2) + q2(3)*q2(3))*sinov(11)

! get pos S = RaE.posC4F
  s1(1) = posC5F(1)-posC5(1)
  s1(2) = posC5F(2)-posC5(2)
  s1(3) = posC5F(3)-posC5(3)
  s(1) = sum(RaE(1,:)*s1(:))
  s(2) = sum(RaE(2,:)*s1(:))
  s(3) = sum(RaE(3,:)*s1(:))

! find u
  u(1) = sum(RaE(1,:)*con7(:))
  u(2) = sum(RaE(2,:)*con7(:))
  u(3) = sum(RaE(3,:)*con7(:))
  u1 = sqrt(sum(u(:)*u(:)))
  u(:) = u(:)/u1
!
  findsolu(:) = 0
!
! first scan the interval for imaginary solutions in all of the four branches
!
  w1 = -PI
  g1 = 0.0
  npreintvs = 0
  increment = 0.0
  logi1b = .false.
  kk = 0
!
  do while (1.eq.1)
    oldws(1) = w1
    oldgs(1) = g1
    w1 = w1 + increment
    kk = kk + 1
    if (w1.gt.PI) exit
    
    call checkreal1_nc(w1,logi1,g1,curslope(1))
    if (curslope(1).ne.0.0) then
      newt = abs(g1/curslope(1))
    else if ((increment.ne.0.0).AND.(g1.ne.oldgs(1))) then
      newt = abs(g1*increment/(g1-oldgs(1)))
    else
      newt = 1.0e-8 ! bail
    end if
    if (g1.le.0.0) then
      if (oldgs(1).gt.0.0) then
!        write(*,*) 'accidental cross out',RADIAN*w1,RADIAN*oldws(1)
        preintvbs(npreintvs,2) = oldws(1)
        logi1b = .false.
        increment = 1.0d-9
        cycle
      end if
!     can only use g1-slope for guidance
      if (curslope(1).gt.0.0) then
!       converging case
        increment = min(UJ_params(4)/RADIAN,0.5*newt)
        increment = max(min(increment,0.5*(PI-w1)),1.0d-9)
        if ((increment+w1).eq.w1) then
!          write(*,*) 'convergent switch',w1*RADIAN
          increment = 1.0d-9
          cycle
        end if
      else
!       diverging case
        increment = min(UJ_params(4)/RADIAN,0.5*newt)
        increment = max(min(increment,0.5*(PI-w1)),1.0d-9)
        if ((increment+w1).eq.w1) then
!          write(*,*) 'divergent switch',w1*RADIAN!,g1,curslope(1)
          increment = 1.0d-9
          cycle
        end if
      end if
    else ! g1 greater zero
      if (oldgs(1).le.0.0) then
!        write(*,*) 'accidental cross in',RADIAN*w1,RADIAN*oldws(1)
        npreintvs = npreintvs + 1
        preintvbs(npreintvs,1) = w1
        logi1b = .true.
        increment = 1.0d-9
        cycle
      end if
      if (curslope(1).le.0.0) then
!       converging case
        increment = min(UJ_params(4)/RADIAN,0.5*newt)
        increment = max(min(increment,0.5*(PI-w1)),1.0d-9)
        if ((increment+w1).eq.w1) then
!          write(*,*) 'convergent switch',w1*RADIAN
          increment = 1.0d-9
          cycle
        end if
      else
!       diverging case
        increment = min(UJ_params(4)/RADIAN,0.5*newt)
        increment = max(min(increment,0.5*(PI-w1)),1.0d-9)
        if ((increment+w1).eq.w1) then
!          write(*,*) 'divergent switch',w1*RADIAN!,g1,curslope(1)
          increment = 1.0d-9
          cycle
        end if
      end if
    end if
  end do
!
  if (logi1b.EQV..true.) then
    preintvbs(npreintvs,2) = PI
  end if
!  write(*,*) 'PRE',npreintvs
!  do j=1,npreintvs
!    write(*,*) RADIAN*preintvbs(j,1),RADIAN*preintvbs(j,2)
!  end do
!
! first bail-out for no solution
  if (npreintvs.eq.0) then
    solucount = 0
    call mcmove_nuccr(rsi,rsf,dof,dphi,athree)
    return
  end if
!
  do kkk=2,1,-1
    nintvs(kkk) = 0
    do jj=1,npreintvs 
      increment = 0.0
      w1 = preintvbs(jj,1)
      kk = 0
      logi1b = .false.
      logi1a = .true.
      oldgs(:) = -1.0
      do while (1.eq.1)
        oldws(1) = w1
        w1 = w1 + increment
        kk = kk + 1
        if (w1.gt.preintvbs(jj,2)) exit
        maxincr = max(1.0d-9,0.5*(preintvbs(jj,2)-w1))
        call checkreal2_nc(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),2*kkk,atrue) ! branch 2 and 1 are identical, as are 3 and 4
        if (curslope(2).ne.0.0) then
          newt = abs(gn/curslope(2))
        else if ((increment.ne.0.0).AND.(gn.ne.oldgs(2))) then
          newt = abs(gn*increment/(gn-oldgs(2)))
        else
          newt = 1.0e-8 ! bail
        end if
!
!       initial segment
        if ((logi1a.EQV..true.).AND.(gn.gt.0.0)) then
          logi1a = .false.
!          write(*,*) 'first start',RADIAN*w1
          nintvs(kkk) = nintvs(kkk) + 1
          intvalbnds(kkk,nintvs(kkk),1) = w1
          logi1b = .true.
        else if ((logi1a.EQV..true.).AND.(gn.le.0.0)) then
          logi1a = .false.
!       accidental cross out
        else if ((oldgs(2).gt.0.0).AND.(gn.le.0.0)) then
!          write(*,*) 'accidental cross out',RADIAN*w1,RADIAN*oldws(1)
          intvalbnds(kkk,nintvs(kkk),2) = oldws(1)
          logi1b = .false.
          increment = 1.0d-9
          oldgs(2) = gn
          cycle
!       accidental cross in
        else if ((oldgs(2).le.0.0).AND.(gn.gt.0.0)) then
!          write(*,*) 'accidental cross in',RADIAN*w1
          nintvs(kkk) = nintvs(kkk) + 1
          intvalbnds(kkk,nintvs(kkk),1) = w1
          logi1b = .true.
          increment = 1.0d-9
          oldgs(2) = gn
          cycle
        end if
        if (gn*curslope(2).lt.0.0) then
!         converging case
          increment = min(UJ_params(4)/RADIAN,0.8*newt,maxincr)
          oldgs(2) = gn
          if ((increment+w1).eq.w1) then
!          write(*,*) 'convergent switch'
            if (gn.gt.0.0) then
              if (logi1b.EQV..false.) then
!                write(*,*) 'convergent start',RADIAN*w1
                nintvs(kkk) = nintvs(kkk) + 1
                intvalbnds(kkk,nintvs(kkk),1) = w1
                logi1b = .true.
                increment = 1.0d-9
                cycle
              else
!                write(*,*) 'convergent end',RADIAN*w1
                intvalbnds(kkk,nintvs(kkk),2) = w1
                logi1b = .false.
                increment = 1.0d-9
                cycle
              end if
            else
              increment = 1.0d-9
              cycle
            end if
          end if
        else
!          diverging case
          oldgs(2) = gn
          increment = min(UJ_params(4)/RADIAN,1.0*newt,maxincr)
!         gn may start extremely small 
          if ((increment+w1).eq.w1) then
            if (gn.gt.0.0) then
              if (logi1b.EQV..false.) then
!                write(*,*) 'divergent start',RADIAN*w1
                nintvs(kkk) = nintvs(kkk) + 1
                intvalbnds(kkk,nintvs(kkk),1) = w1
                logi1b = .true.
                increment = 1.0d-9
                cycle
              else
!                write(*,*) 'divergent end',RADIAN*w1
                intvalbnds(kkk,nintvs(kkk),2) = w1
                logi1b = .false.
                increment = 1.0d-9
                cycle
              end if
            else
              increment = 1.0d-9
              cycle
            end if
          end if
        end if
!      write(*,*) RADIAN*w1,g1,curslope(1),gn,curslope(2),logi1
      end do
      if (logi1b.EQV..true.) then
        intvalbnds(kkk,nintvs(kkk),2) = preintvbs(jj,2)
      end if
    end do
!
    if ((logi1b.EQV..false.).AND.(logi1.EQV..true.)) then
!     this shouldn't really happen
      call checkreal1_nc(oldws(1),logi1,g1,curslope(1))
      if (logi1.EQV..true.) then
        nintvs(kkk) = nintvs(kkk) + 1
        intvalbnds(kkk,nintvs(kkk),1) = oldws(1)
        intvalbnds(kkk,nintvs(kkk),2) = oldws(1)
      end if
    end if
!
!   now check the identified intervals
    if (nintvs(kkk).gt.0) then
!     check for vanishing intervals
      do j=1,nintvs(kkk)
        if (abs(intvalbnds(kkk,j,1)-intvalbnds(kkk,j,2)).lt.2.0d-9) then
!          write(*,*) 'vanishing interval:',j,intvalbnds(kkk,j,1)
          do jj=2*kkk-1,2*kkk
            w1 = 0.5*(intvalbnds(kkk,j,1)+intvalbnds(kkk,j,2))
            call checkreal2_nc(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2_nc(w1,newccmat2,g1,jj,really)
!             really-check is redundant since get_newccmat2_nc calls checkreal2_nc
              if ((really.EQV..true.).AND.(g1.eq.0.0)) then
                findsolu(jj) = findsolu(jj) + 1
                w1mat(jj,findsolu(jj)) = w1
              end if
            end if
          end do
          do jj=j,nintvs(kkk)-1
            intvalbnds(kkk,jj,1) = intvalbnds(kkk,jj+1,1)
            intvalbnds(kkk,jj,2) = intvalbnds(kkk,jj+1,2)
          end do
          nintvs(kkk) = nintvs(kkk) - 1
        end if
      end do
!     unify if necessary
      if (nintvs(kkk).gt.1) then
        if ((intvalbnds(kkk,1,1).le.-PI).AND.(intvalbnds(kkk,nintvs(kkk),2).ge.PI)) then
!        write(*,*) 'wrapped around'
          intvalbnds(kkk,1,1) = intvalbnds(kkk,nintvs(kkk),1)
          intvalbnds(kkk,1,2) = intvalbnds(kkk,1,2) + 2.0*PI
          nintvs(kkk) = nintvs(kkk) - 1
        end if
      end if
      do j=1,nintvs(kkk)
!       shorten the interval slightly to be numerically safe
        intvalbnds(kkk,j,1) = intvalbnds(kkk,j,1) + 2.0d-9
        intvalbnds(kkk,j,2) = intvalbnds(kkk,j,2) - 2.0d-9
!       check one more time
        if (abs(intvalbnds(kkk,j,1)-intvalbnds(kkk,j,2)).lt.4.0d-9) then
!         write(*,*) 'vanishing interval:',j,intvalbnds(kkk,j,1)
          do jj=2*kkk-1,2*kkk
            w1 = 0.5*(intvalbnds(kkk,j,1)+intvalbnds(kkk,j,2))
            call checkreal2_nc(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2_nc(w1,newccmat2,g1,jj,really)
!             see above
              if ((really.EQV..true.).AND.(g1.eq.0.0)) then
                findsolu(jj) = findsolu(jj) + 1
                w1mat(jj,findsolu(jj)) = w1
              end if
            end if
          end do
          do jj=j,nintvs(kkk)-1
            intvalbnds(kkk,jj,1) = intvalbnds(kkk,jj+1,1)
            intvalbnds(kkk,jj,2) = intvalbnds(kkk,jj+1,2)
          end do
          nintvs(kkk) = nintvs(kkk) - 1
        end if
      end do
    end if
!    write(*,*) 'NU',kkk,kk,' steps'
!    do jj=1,nintvs(kkk)
!    write(*,*) kkk,jj,RADIAN*intvalbnds(kkk,jj,1:2)
!    end do
  end do
!
! and lastly scan them
  kk = 0
  do kkk=1,4
    jj = 1
    if (kkk.ge.3) jj = 2
    do j=1,nintvs(jj)
      lastslope(:) = 0.0
      increment = 0.0
      oldws(:) = 0.0
      oldgs(:) = 0.0
      g1 = 0.0
      w1 = intvalbnds(jj,j,1)
!
      do while (1.eq.1)
        kk = kk + 1
        oldws(kkk) = w1
        if (w1.ge.intvalbnds(jj,j,2)) exit
        if ((w1.lt.intvalbnds(jj,j,2)).AND.(w1+increment.gt.intvalbnds(jj,j,2))) then
          w1 = intvalbnds(jj,j,2)
        else
          w1 = w1 + increment
        end if
        oldgs(1) = g1
        call get_newccmat2_nc(w1,newccmat2,g1,kkk,really)
!       bail out if we nail it
        if ((really.EQV..true.).AND.(g1.eq.0.0)) then
!          write(*,*) 'nailed it ',w1*RADIAN,' in ',kkk
          findsolu(kkk) = findsolu(kkk) + 1
          w1mat(kkk,findsolu(kkk)) = w1
          increment = 1.0d-9
          cycle
        end if
!       get approximation of slope 
        ltestw1 = w1 - 1.0d-9
        htestw1 = w1 + 1.0d-9
        if (ltestw1.lt.intvalbnds(jj,j,1)) then
          ltestw1 = intvalbnds(jj,j,1)
          searchstp = 1.0d-9
          gnpr = g1
          really3 = really
          call get_newccmat2_nc(htestw1,newccmat2,gnpo,kkk,really2)
!         bail out if we nail it
          if ((really2.EQV..true.).AND.(gnpo.eq.0.0)) then
!          write(*,*) 'nailed it fpo',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = htestw1
            increment = 1.0d-9
            cycle
          end if
        else if (htestw1.gt.intvalbnds(jj,j,2)) then
          htestw1 = intvalbnds(jj,j,2)
          searchstp = 1.0d-9
          really2 = really
          call get_newccmat2_nc(ltestw1,newccmat2,gnpr,kkk,really3)
!         bail out if we nail it
          if ((really3.EQV..true.).AND.(gnpr.eq.0.0)) then
!          write(*,*) 'nailed it fpr',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = ltestw1
            increment = 1.0d-9
            cycle
          end if
          gnpo = g1
        else
          searchstp = 2.0d-9
          call get_newccmat2_nc(ltestw1,newccmat2,gnpr,kkk,really3)
!         bail out if we nail it
          if ((really3.EQV..true.).AND.(gnpr.eq.0.0)) then
!          write(*,*) 'nailed it pr',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = ltestw1
            increment = 1.0d-9
            cycle
          end if
          call get_newccmat2_nc(htestw1,newccmat2,gnpo,kkk,really2)
          if ((really2.EQV..true.).AND.(gnpo.eq.0.0)) then
!          write(*,*) 'nailed it po',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = htestw1
            increment = 1.0d-9
            cycle
          end if
        end if
!       we build in a sanity check to make sure gnpo and gnpr are not NaN (other inifinite looping looms)
        if ((really.EQV..false.).OR.(really2.EQV..false.).OR.(really3.EQV..false.)) then
!         cancel this interval
          nc_wrncnt(2) = nc_wrncnt(2) + 1
          if (nc_wrncnt(2).eq.nc_wrnlmt(2)) then
            write(ilog,*) 'Warning. Encountered interval inconsistency in loop closure of Dinner-Ulmschneider&
 & algorithm for nucleic acids. This means a number was tested which was previously found to give a real solu&
 &tion and now produced NaN. This is ultimately a bug the causes of which are currently poorly understood.'
            write(ilog,*) 'This was warning #',nc_wrncnt(2),' of this type not all of which may be displayed.'
            if (10.0*nc_wrnlmt(2).gt.0.5*HUGE(nc_wrnlmt(2))) then
              nc_wrncnt(2) = 0 ! reset
            else
              nc_wrnlmt(2) = nc_wrnlmt(2)*10
            end if
          end if
          exit
        end if
!       we can bisect the search interval straight away
        if (((gnpr.lt.0.0).AND.(gnpo.gt.0.0)).OR.((gnpr.gt.0.0).AND.(gnpo.lt.0.0))) then
!          write(*,*) 'bisecting ',ltestw1*RADIAN,' and ',htestw1*RADIAN,' in ',kkk
          call bisection2_nuc(ltestw1,htestw1,kkk,findsolu(kkk))
          increment = 2.0d-9
          g1 = 0.0 ! sets oldgs to zero
          cycle
        end if
!       check for accidental crossover
        if (((oldgs(kkk).lt.0.0).AND.(g1.gt.0.0)).OR.((oldgs(kkk).gt.0.0).AND.(g1.lt.0.0))) then ! yes
!          write(*,*) 'crossover: ',oldws(kkk)*RADIAN,' and ',w1*RADIAN,' in ',kkk
          call bisection2_nuc(oldws(kkk),w1,kkk,findsolu(kkk))
          increment = 1.0d-9
          cycle
        end if
!       now use Newton method to advance
        curslope(kkk) = (gnpo-gnpr)/searchstp
        if (curslope(kkk).ne.0.0) then
          newt = abs(g1/curslope(kkk))
        else if (lastslope(kkk).ne.0.0) then
          newt = abs(g1/lastslope(kkk))
        else
          newt = 1.0e-8 ! bail
        end if
        if (g1*curslope(kkk).lt.0.0) then
          increment = min(UJ_params(4)/RADIAN,0.5*newt)
        else
          increment = min(UJ_params(4)/RADIAN,0.8*newt)
        end if
        increment = max(min(increment,0.5*(intvalbnds(jj,j,2)-w1)),1.0d-9)
!        increment = min(UJ_params(4)/RADIAN,0.5*abs(g1/curslope(kkk)))
        if ((w1+increment).eq.w1) then
!         this means we approach a root with a limiting slope of zero
!          write(*,*) 'convergent root: ',w1*RADIAN,' in ',kkk
          findsolu(kkk) = findsolu(kkk) + 1
          w1mat(kkk,findsolu(kkk)) = w1
          increment = 1.0d-9
          cycle
        end if
        lastslope(kkk) = curslope(kkk)
      end do
    end do
  end do
!
  solucount = sum(findsolu(1:4))
!
  if (solucount.eq.0) then    
!    just exit if we found no solution
!    write(ilog,*) 'Did not find a solution to the chain closure'
    call mcmove_nuccr(rsi,rsf,dof,dphi,athree)
    return
  else if (solucount.eq.1) then
!   for a single solution it's easy: assign, restore, perturb, get jacobian, exit
    dpfmat(1:dof) = dphi(1:dof)
    do j=1,4
      if(findsolu(j).gt.0) then
        call get_newccmat2_nc(w1mat(j,1),newccmat2,g1,j,really)
        goodsol = .true.
        do jjj=1,6
          if ((newccmat2(j,jjj).lt.-PI).OR.(newccmat2(j,jjj).gt.PI).OR.(really.EQV..false.)) then
            goodsol = .false.
          end if
        end do
        if (goodsol.EQV..false.) then
          solucount = solucount - 1
          call mcmove_nuccr(rsi,rsf,dof,dphi,athree)
          return
        end if
        incrmat(1,1:6) = newccmat2(j,1:6) - ccmat_nuc(1:6)
        exit
      end if
    end do
    dpfmat(dof+1:dof+6) = incrmat(1,1:6)
!
!   restore pre-rotation segment (initial state fully restored)
    call mcmove_nuccr(rsi,rsf,dof,dphi,athree)
!
!   move entire internal chain
    call mcmove_nuccr(rsi,rsf,dof,dpfmat,aone)
    if ((abs(x(nuci(rsf,6))-posC3F(1))+abs(y(nuci(rsf,6))-posC3F(2))+abs(z(nuci(rsf,6))-posC3F(3))).gt.1.0e-4) then
      solucount = solucount - 1
      call mcmove_nuccr(rsi,rsf,dof,dpfmat,atwo)
      return
    end if
!
!   get jacobian for new state 
    call jacobian_nuccr(rsf,jac_b_vec(1))
!   restore
    call mcmove_nuccr(rsi,rsf,dof,dpfmat,atwo)
!
  else
!
!   for multiple solutions it's more work: for all: assign, restore, perturb, get jacobian, then pick randomly
    dpfmat(1:dof) = dphi(1:dof)
!
    kk = 0
    logi1 = .false.
    do j=1,4
!
      do kkk=1,findsolu(j)
        kk = kk + 1
        call get_newccmat2_nc(w1mat(j,kkk),newccmat2,g1,j,really)
        goodsol = .true.
        do jjj=1,6
          if ((newccmat2(j,jjj).lt.-PI).OR.(newccmat2(j,jjj).gt.PI).OR.(really.EQV..false.)) then
            goodsol = .false.
            exit
          end if
        end do
        if (goodsol.EQV..false.) then
          solucount = solucount - 1
          do jjj=kkk,findsolu(j)-1
            w1mat(j,jjj) = w1mat(j,jjj+1) 
          end do
          findsolu(j) = findsolu(j) - 1
          kk = kk - 1
          cycle
        end if
        incrmat(kk,1:6) = newccmat2(j,1:6) - ccmat_nuc(1:6)
!
        dpfmat(dof+1:dof+6) = incrmat(kk,1:6)
!
!       restore everything appropriately
        if (logi1.EQV..false.) then
          call mcmove_nuccr(rsi,rsf,dof,dphi,athree)
        else
          call mcmove_nuccr(rsi,rsf,dof,dpfmat,afive)
        end if
        logi1 = .true.
!
!       move entire internal chain
        call mcmove_nuccr(rsi,rsf,dof,dpfmat,aone)
        if ((abs(x(nuci(rsf,6))-posC3F(1))+abs(y(nuci(rsf,6))-posC3F(2))+abs(z(nuci(rsf,6))-posC3F(3))).gt.1.0e-4) then
          solucount = solucount - 1
          do jjj=kkk,findsolu(j)-1
            w1mat(j,jjj) = w1mat(j,jjj+1) 
          end do
          findsolu(j) = findsolu(j) - 1
          kk = kk - 1
          cycle
        end if
!
!       get jacobian for new state 
        call jacobian_nuccr(rsf,jac_b_vec(kk))
!
      end do
      
    end do
!
    if (logi1.EQV..false.) then
      call mcmove_nuccr(rsi,rsf,dof,dphi,athree)
    else
      call mcmove_nuccr(rsi,rsf,dof,dpfmat,atwo)
    end if
!
  end if
!
  CONTAINS
! -------------------------------------------------------------------  
!
! r(w1) = invT1*invR1*(s-q0)
subroutine r_of_w1w2_nc(r,w1val,rnorm,svq)
!
  RTYPE, INTENT(OUT):: r(3),svq(3),rnorm
  RTYPE, INTENT(IN):: w1val
!
!
! set q0 = p1 + R1.T1.p2
  q0(1) = p1 + cosov(8)*p2
  q0(2) = cos(w1val)*sinov(8)*p2
  q0(3) = sin(w1val)*sinov(8)*p2
!
  svq(:) = s(:) - q0(:)
!  write(*,*) "q0",q0

  call genR(w1val,R1,invR1)
  invmat = matmul(invT2,matmul(invR2,matmul(invT1,invR1)))
!
  r(1) = sum(invmat(1,:)*svq(:))
  r(2) = sum(invmat(2,:)*svq(:))
  r(3) = sum(invmat(3,:)*svq(:))
!
  rnorm = sum(r(:)*r(:))
!
!
end subroutine r_of_w1w2_nc
!
! -------------------------------------------------------------------  
!
! get arg1, arg2, and w and (on request) some derivatives
! 
subroutine deriva_w1_nc(sw1,cw1,rnorm,r,svq,arg1,arg2,darg1dw1,darg2dw1,w,drdw1,dwndw1,do_deriv)
!
  RTYPE, INTENT(IN):: r(3),sw1,cw1,rnorm,svq(3)
  logical, INTENT(IN):: do_deriv
  RTYPE, INTENT(OUT):: arg1,arg2,darg1dw1,darg2dw1,w,drdw1(3),dwndw1
  RTYPE sa1,sa2,ca1,ca2,sw2,cw2
!
  sa1 = sinov(8)
  sa2 = sinov(9)
  ca1 = cosov(8)
  ca2 = cosov(9)
  sw2 = sinccn7
  cw2 = cosccn7
!
  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
  arg2 = w*w*q1(2)*q1(2)
!
  if (do_deriv.EQV..true.) then
    drdw1(1) = -svq(2)*(sw1*sa1*ca2 + sa2*(sw1*ca1*cw2 + cw1*sw2)) + &
              & svq(3)*(cw1*sa1*ca2 + sa2*(cw1*ca1*cw2 - sw1*sw2)) + &
             & (sa1*cw1*ca2 + sa2*(ca1*cw1*cw2 - sw1*sw2))*sw1*sa1*p2 - &
             & (sa1*sw1*ca2 + sa2*(ca1*sw1*cw2 + cw1*sw2))*cw1*sa1*p2
    drdw1(2) = svq(2)*(sw1*sa1*sa2 - ca2*(sw1*ca1*cw2 + cw1*sw2)) - &
             & svq(3)*(cw1*sa1*sa2 - ca2*(cw1*ca1*cw2 - sw1*sw2)) + &
           & (-sa1*cw1*sa2 + ca2*(ca1*cw1*cw2 - sw1*sw2))*sw1*sa1*p2 - &
           & (-sa1*sw1*sa2 + ca2*(ca1*sw1*cw2 + cw1*sw2))*cw1*sa1*p2
    drdw1(3) = svq(2)*(sw1*ca1*sw2 - cw1*cw2) - svq(3)*(cw1*ca1*sw2 + sw1*cw2) - &
            & (ca1*cw1*sw2 + sw1*cw2)*sw1*sa1*p2 - cw1*sa1*p2*(-ca1*sw1*sw2 + cw1*cw2)
    dwndw1 = (1.0/(2.0*q1(2)))*(-2.0*drdw1(1)*q1(1) + sum(2.0*r(:)*drdw1(:)))
    darg1dw1 = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw1(2) + r(3)*drdw1(3)))
    darg2dw1 = 2.0*q1(2)*q1(2)*w*dwndw1
  end if
!
end subroutine deriva_w1_nc
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal1_nc(w1val,isreal,howfaroff,testslope)
!
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: howfaroff,testslope
  RTYPE arg1,arg2,sw1,cw1
  RTYPE drdw(3),darg1dr,darg2dr,dwndr,svq(3),rnorm,r(3)
  logical, INTENT(OUT):: isreal
  logical atrue
!
  isreal = .false.
  atrue = .true.
!
  sw1 = sin(w1val)
  cw1 = cos(w1val)
!
  call r_of_w1w2_nc(r,w1val,rnorm,svq)
!
  call deriva_w1_nc(sw1,cw1,rnorm,r,svq,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,atrue)
!
!  drdw(1) = -svq(2)*(sw1*sa1*ca2 + sa2*(sw1*ca1*cw2 + cw1*sw2)) + &
!           & svq(3)*(cw1*sa1*ca2 + sa2*(cw1*ca1*cw2 - sw1*sw2)) + &
!           & (sa1*cw1*ca2 + sa2*(ca1*cw1*cw2 - sw1*sw2))*sw1*sa1*p2 - &
!           & (sa1*sw1*ca2 + sa2*(ca1*sw1*cw2 + cw1*sw2))*cw1*sa1*p2
!  drdw(2) = svq(2)*(sw1*sa1*sa2 - ca2*(sw1*ca1*cw2 + cw1*sw2)) - &
!           & svq(3)*(cw1*sa1*sa2 - ca2*(cw1*ca1*cw2 - sw1*sw2)) + &
!           & (-sa1*cw1*sa2 + ca2*(ca1*cw1*cw2 - sw1*sw2))*sw1*sa1*p2 - &
!           & (-sa1*sw1*sa2 + ca2*(ca1*sw1*cw2 + cw1*sw2))*cw1*sa1*p2
!  drdw(3) = svq(2)*(sw1*ca1*sw2 - cw1*cw2) - svq(3)*(cw1*ca1*sw2 + sw1*cw2) - &
!          & (ca1*cw1*sw2 + sw1*cw2)*sw1*sa1*p2 - cw1*sa1*p2*(-ca1*sw1*sw2 + cw1*cw2)
!!
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  dwnewdr = (1.0/(2.0*q1(2)))*(-2.0*drdw(1)*q1(1) + sum(2.0*r(:)*drdw(:)))
!!
!  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  darg1dr = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw(2) + r(3)*drdw(3)))
!  arg2 = w*w*q1(2)*q1(2)
!  darg2dr = (2.0*q1(2)*q1(2))*w*dwnewdr
!
  howfaroff = arg1 - arg2
  testslope = darg1dr - darg2dr
  if (howfaroff.gt.0.0) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
end subroutine checkreal1_nc
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal2_nc(w1val,isreal2,howfaroff,targ1,tslp1,tslp2,w3val,sw3,cw3,t1ex,r,t,do_deriv)
!
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: w3val,sw3,cw3,t1ex,r(3),targ1,tslp1,tslp2,howfaroff
  RTYPE arg1,arg2,darg1dr,darg2dr,dwndr,drdw(3)
  RTYPE denom,enum11,enum21,enum2,sarg1,sarg2,denum2dr,sarg21
  RTYPE dsarg1dr,dsarg2dr,denum21dr,denum11dr,dw3valdr,dt1exdr,enuma,enumb
  RTYPE denumadr,denumbdr,sw1,cw1,svq(3),rnorm
  logical, INTENT(OUT):: isreal2
  logical, INTENT(IN):: do_deriv
  logical isreal
  integer, INTENT(IN):: t
!
  tslp2 = 0.0
  targ1 = 0.0
  isreal = .false.
!
  call r_of_w1w2_nc(r,w1val,rnorm,svq)
!
  sw1 = sin(w1val)
  cw1 = cos(w1val)
!
  call deriva_w1_nc(sw1,cw1,rnorm,r,svq,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,do_deriv)
!
  howfaroff = arg1 - arg2
  if (do_deriv.EQV..true.) tslp1 = darg1dr - darg2dr
  if (howfaroff.gt.0.0) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
  isreal2 = .false.
  if (isreal.EQV..true.) then
    denom = arg1
    enum2 = sqrt(arg1 - arg2)
    enuma = r(2)*q1(2) + r(3)*q1(3) 
    enumb = r(3)*q1(2) - r(2)*q1(3)
    if ((t.eq.1).OR.(t.eq.2)) then
      enum21 = enumb*w*q1(2) - enuma*enum2
      enum11 = enuma*w*q1(2) + enumb*enum2
    else if  ((t.eq.3).OR.(t.eq.4)) then
      enum21 = enumb*w*q1(2) + enuma*enum2
      enum11 = enuma*w*q1(2) - enumb*enum2
    else
      write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
 &ain closure algorithm) must be between 1 and 4.'
      call fexit()
    end if
!
    sarg1 = enum11/denom
    sarg2 = enum21/denom
    sarg21 = sarg2/sarg1
!
    w3val = atan2(sarg2,sarg1) ! or use ccmat
    sw3 = sin(w3val)
    cw3 = cos(w3val)
    call genR(w3val,R3,invR3)
!
    t1ex = cosov(10)*(r(1)-q1(1)) + sinov(10)*(r(2)*cw3+r(3)*sw3-q1(2))
    targ1 = q2sq*sinov(11)*sinov(11) + 2.0*q2(1)*t1ex*cosov(11)- q2(1)*q2(1) - t1ex*t1ex
!
    if (do_deriv.EQV..true.) then
      denum2dr = tslp1/(2.0*enum2)
      denumadr = drdw(2)*q1(2) + drdw(3)*q1(3)
      denumbdr = drdw(3)*q1(2) - drdw(2)*q1(3)
      if ((t.eq.1).OR.(t.eq.2)) then
        denum21dr = -denumadr*enum2 - enuma*denum2dr + q1(2)*(dwndr*enumb + w*denumbdr)
        denum11dr = denumbdr*enum2 + enumb*denum2dr + q1(2)*(dwndr*enuma + w*denumadr)
      else if  ((t.eq.3).OR.(t.eq.4)) then
        denum21dr = denumadr*enum2 + enuma*denum2dr + q1(2)*(dwndr*enumb + w*denumbdr)
        denum11dr = -denumbdr*enum2 - enumb*denum2dr + q1(2)*(dwndr*enuma + w*denumadr)
      end if
      dsarg1dr = (denom*denum11dr - enum11*enum2*darg1dr)/(denom*denom)
      dsarg2dr = (denom*denum21dr - enum21*enum2*darg1dr)/(denom*denom)
      dw3valdr = (1.0/(1.0 + sarg21*sarg21))*((sarg1*dsarg2dr - dsarg1dr*sarg2)/(sarg1*sarg1))
      dt1exdr = cosov(10)*drdw(1) + sinov(10)*(-r(2)*sw3*dw3valdr + &
 & drdw(2)*cw3 + drdw(3)*sw3 + r(3)*cw3*dw3valdr)
      tslp2 = 2.0*q2(1)*cosov(11)*dt1exdr - 2.0*t1ex*dt1exdr
    end if
!
    if (targ1.gt.0.0) then
      isreal2 = .true.
    else
      isreal2 = .false.
    end if
  end if
!
end subroutine checkreal2_nc
!
!------------------------------------------------------------------------------------------------------
! 
subroutine bisection2_nuc(w1,next,t,nsolu)
!
  implicit none
!
  integer, INTENT(IN):: t
  integer nsolu
  RTYPE w1,next,w1val,nextw1,grr,grr2
  logical really,really2
!
  w1val = w1
  nextw1 = next
  epsilon = 2.0*(10.0**(-15.0))
! ! start loop
  do while (abs(nextw1 - w1val).gt.(2.0*epsilon))
!
    midpoint = (nextw1 + w1val) / 2.0
!
    call get_newccmat2_nc(w1val,newccmat2,grr,t,really)
    call get_newccmat2_nc(midpoint,newccmat2,grr2,t,really2)
    if ((really.EQV..false.).OR.(really2.EQV..false.)) then
      nc_wrncnt(3) = nc_wrncnt(3) + 1
      if (nc_wrncnt(3).eq.nc_wrnlmt(3)) then
        write(ilog,*) 'Bisection method failed for Dinner-Ulmschneider concerted rotation&
 & algorithm for polynucleotides. This indicates significantly too coarse search&
 & settings (interval inconsistencies).'
        write(ilog,*) 'This was warning #',nc_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*nc_wrnlmt(3).gt.0.5*HUGE(nc_wrnlmt(3))) then
          nc_wrncnt(3) = 0 ! reset
        else
          nc_wrnlmt(3) = nc_wrnlmt(3)*10
        end if
      end if
      return
    end if
 !
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
  if (w1val.gt.PI) w1val = w1val - 2.0*PI
  if (w1val.le.-PI) w1val = w1val + 2.0*PI
!
  nsolu = nsolu + 1
  w1mat(t,nsolu) = w1val
!
end subroutine bisection2_nuc
! -------------------------------------------------------------------  
!!
!! ! subroutine to find w3
!subroutine get_w3(w1val,w3val,t)
!!
!  integer t
!  RTYPE w,arg1,arg2,w1val,w3val,denom,enum2,enum11,enum12,rnorm,svq(3),r(3)
!!
!  call r_of_w1w2_nc(r,w1val,rnorm,svq)
!
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  denom = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  enum2 = sqrt(denom - w*w*q1(2)*q1(2))
!  enum11 = r(2)*q1(2) + r(3)*q1(3)
!  enum12 = r(3)*q1(2) - r(2)*q1(3)
!     
 ! if ((t.eq.1).OR.(t.eq.2)) then
 !   arg1 = (enum11*w*q1(2) + enum12*enum2)/denom
 !   arg2 = (enum12*w*q1(2) - enum11*enum2)/denom! r(2)*q1(3))*enum2)/denom
 ! else if ((t.eq.3).OR.(t.eq.4)) then
 !   arg1 = (enum11*w*q1(2) - enum12*enum2)/denom
!    arg2 = (enum12*w*q1(2) + enum11*enum2)/denom! r(2)*q1(3))*enum2)/denom
!  else
!    write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
! &ain closure algorithm) must be between 1 and 4.'
!    call fexit()
!  end if
!!
!  w3val = atan2(arg2,arg1) ! or use ccmat
!  call genR(w3val,R3,invR3)
!!
!end subroutine get_w3
! -------------------------------------------------------------------  
!
! subroutine to find w4
subroutine get_newccmat2_nc(w1val,newccmat2,g,t,isreal2)
!
  integer, INTENT(IN):: t
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: newccmat2(4,6),g
  logical, INTENT(OUT):: isreal2
  RTYPE arg1,arg2,w3val,t1ex,t1y,t1z,t2val(3),w4val,w5val,w6val,r(3),howfaroff,targ1
  RTYPE ay,pos7(3),TRTR(3),R1TR(3),T1R2x(3),T1R2(3,3),T2R3x(3),T2R3(3,3)
  RTYPE T4R5p5(3),T3R4(3,3),T3R4x(3),T4R5(3,3),ctw5,stw5,stw3,ctw3,enum1
  logical afalse
!
  afalse = .false.
!
! first get w3
  call checkreal2_nc(w1val,isreal2,howfaroff,targ1,nouse(1),nouse(2),w3val,stw3,ctw3,t1ex,r,t,afalse)
  call genR(w3val,R3,invR3)
!
  if (isreal2.EQV..false.) then
    return ! deal on the outside
  end if

!  t1x = cosov(10)*(r(1)-q1(1)) + sinov(10)*(r(2)*ctw3+r(3)*stw3-q1(2))
!  denom = (q2(2)*q2(2) + q2(3)*q2(3))*sinov(11)
  enum1 = sqrt(q2sq*sinov(11)*sinov(11) + 2.0*q2(1)*t1ex*cosov(11) - q2(1)*q2(1) - t1ex*t1ex)
!!
  if ((t.eq.1).OR.(t.eq.3)) then
    arg1 = (q2(2)*(q2(1)*cosov(11) - t1ex) + q2(3)*enum1)/w5denom
    arg2 = (q2(3)*(-q2(1)*cosov(11) + t1ex) - q2(2)*enum1)/w5denom
  else if ((t.eq.2).OR.(t.eq.4)) then
    arg1 = (q2(2)*(q2(1)*cosov(11) - t1ex) - q2(3)*enum1)/w5denom
    arg2 = (q2(3)*(-q2(1)*cosov(11) + t1ex) + q2(2)*enum1)/w5denom
  else
    write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
 &ain closure algorithm) must be between 1 and 4.'
    call fexit()
  end if
!
  w5val = atan2(arg2,arg1)
  ctw5 = cos(w5val)
  stw5 = sin(w5val)
!
! -------------------------------------------------------------------  
!
  t2val(1) = q2(1)*cosov(11) - q2(2)*sinov(11)*ctw5 + q2(3)*sinov(11)*stw5
  t2val(2) = q2(1)*sinov(11) + q2(2)*cosov(11)*ctw5 - q2(3)*cosov(11)*stw5
  t2val(3) = q2(2)*stw5 + q2(3)*ctw5
! 
  t1y = -sinov(10)*(r(1)-q1(1)) + cosov(10)*(r(2)*ctw3 + r(3)*stw3 - q1(2))
  t1z = -r(2)*stw3 + r(3)*ctw3 - q1(3)
!
  arg1 = (t1y*t2val(2) + t1z*t2val(3))/(t1y*t1y + t1z*t1z)
  arg2 = (t1z*t2val(2) - t1y*t2val(3))/(t1y*t1y + t1z*t1z)
  w4val = atan2(arg2,arg1)
!
! -------------------------------------------------------------------  
!
  call genR(w3val,R3,invR3)
  call genT(oldvals_nuc(10),T3,invT3)
  call genR(w4val,R4,invR4)
  call genT(oldvals_nuc(11),T4,invT4)
  call genR(w5val,R5,invR5)
  call genT(oldvals_nuc(12),T5,invT5)
  invmat = matmul(invT5,matmul(invR5,matmul(invT4,matmul(invR4,matmul&
 &(invT3,matmul(invR3,matmul(invT2,matmul(invR2,matmul(invT1,invR1)))))))))
!
  arg1 = sum(invmat(2,:)*u(:))/sinov(13)
  arg2 = sum(invmat(3,:)*u(:))/sinov(13)
!
  w6val = atan2(arg2,arg1)
  newccmat2(t,1) = w1val
  newccmat2(t,2) = w3val
  newccmat2(t,3) = w4val
  newccmat2(t,4) = w5val
  newccmat2(t,5) = w6val
!
  T1R2 = matmul(T1,R2)
  T2R3 = matmul(T2,R3)
  T3R4 = matmul(T3,R4)
  T4R5 = matmul(T4,R5)
  T4R5p5(1) = T4R5(1,1)*p5 
  T4R5p5(2) = T4R5(2,1)*p5
  T4R5p5(3) = T4R5(3,1)*p5
  T4R5p5(1) = T4R5p5(1) + p4
  do k=1,3
    T3R4x(k) = sum(T3R4(k,:)*T4R5p5(:))
  end do
  T3R4x(1) = T3R4x(1) + p3
  do k=1,3
    T2R3x(k) = sum(T2R3(k,:)*T3R4x(:))
  end do
  T2R3x(1) = T2R3x(1) + p2
  do k=1,3
    T1R2x(k) = sum(T1R2(k,:)*T2R3x(:))
  end do
  T1R2x(1) = T1R2x(1) + p1
  do k=1,3
    R1TR(k) = sum(R1(k,:)*T1R2x(:))
  end do
  do k=1,3
    TRTR(k) = sum(invRaE(k,:)*R1TR(:))
  end do
  pos7(:) = TRTR(:) + posC5(:)
  call dihed(pos7,posC5F,posC4F,posC3F,ay)
  newccmat2(t,6) = ay
!
  g = u(1)*invmat(1,1) + u(2)*invmat(1,2) + u(3)*invmat(1,3) - cosov(13)
!
end subroutine get_newccmat2_nc
!!
!! -------------------------------------------------------------------  
!!
!! subroutine outputs g(w1)
!!
!subroutine g_of_w1w3_nuc(w1val,g,t)
!!
!  integer t
!  RTYPE w1val,g
!  RTYPE TRmat(3,3)
!!
!! using branch t get solution for g(w1)
!  call get_newccmat2(w1val,newccmat2,t)
!!
! calculate g(w1)
!  TRmat = matmul(invT5,matmul(invR5,matmul(invT4,matmul(invR4,matmul(invT3,matmul&
! &(invR3,matmul(invT2,matmul(invR2,matmul(invT1,invR1)))))))))
!  g = u(1)*TRmat(1,1) + u(2)*TRmat(1,2) + u(3)*TRmat(1,3) - cosov(13)
!!
!end subroutine g_of_w1w3_nuc
!!
! ------------------------------------------------------------------
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
  Ti(2,1) = -1.0*Ti(1,2)
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
end subroutine nuccrchainclose
