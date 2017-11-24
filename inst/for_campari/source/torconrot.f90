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
subroutine buildtorcrdof_do(rsi,rsf,dof,torcrmat,omepos,nome)
!------------------------------------------------------------
! Subroutine to calculate the matrix I = (dA/dphi)**2 for tor_cr prerotation.
! 
! INPUTS:     rsi, rsf - inital and final fixed residues number
! MODIFIES:   torcrmat - mode1: 3 DOF per residue (torsions only)
!------------------------------------------------------------
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
  integer j,k,m,rsi,rsf,dof,rs,dofcnt,omepos(MAXUJDOF),nome
  RTYPE pos2(3),pos3(3),pos4(3)
  RTYPE torcrmat(MAXUJDOF,MAXUJDOF),doflever(3,MAXUJDOF),dvec(3)
!
! loop over all free residues of the prerotation
  dofcnt = 0
!  omepos = 0 ! (will default to the most conservative value of 3)
  nome = 0
!  whereome = .false.
  do rs=rsf-2,rsi,-1
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
    call get_tor_lever(molofrs(atmres(cai(rsf))),cai(rsf),pos2,pos3,pos4,dvec)
    doflever(:,dof-dofcnt) = dvec(:)/RADIAN
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) exit
!
!   phi angle
    if (seqflag(rs).ne.5) then
      pos4(1) = x(ci(rs-1))
      pos4(2) = y(ci(rs-1))
      pos4(3) = z(ci(rs-1))
      pos3(1) = x(ni(rs))
      pos3(2) = y(ni(rs))
      pos3(3) = z(ni(rs))
      pos2(1) = x(cai(rs))
      pos2(2) = y(cai(rs))
      pos2(3) = z(cai(rs))
      call get_tor_lever(molofrs(atmres(cai(rsf))),cai(rsf),pos2,pos3,pos4,dvec)
      doflever(:,dof-dofcnt) = dvec(:)/RADIAN
      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) exit
    end if
!
!   omega angle
    pos4(1) = x(iz(3,wline(rs)))
    pos4(2) = y(iz(3,wline(rs)))
    pos4(3) = z(iz(3,wline(rs)))
    pos3(1) = x(ci(rs-1))
    pos3(2) = y(ci(rs-1))
    pos3(3) = z(ci(rs-1))
    pos2(1) = x(ni(rs))
    pos2(2) = y(ni(rs))
    pos2(3) = z(ni(rs))
    call get_tor_lever(molofrs(atmres(cai(rsf))),cai(rsf),pos2,pos3,pos4,dvec)
    doflever(:,dof-dofcnt) = dvec(:)/RADIAN
!    if (whereome.EQV..false.) then
!      whereome = .true.
!      omepos = dof-dofcnt
!    end if
    nome = nome + 1
    omepos(nome) = dof-dofcnt
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) exit
!
  end do
!
! ! build symetric matrix I = dA/dphi*dA/dphi (outer_product)
  do j=1,dof
    do k=1,dof
!     ! inner product
      torcrmat(j,k) = 0.0d0
      do m=1,3
        torcrmat(j,k) = torcrmat(j,k) + doflever(m,j)*doflever(m,k)
      end do
    end do
  end do
!
!  do while (omepos.gt.0) 
!    omepos = omepos - 3
!  end do
!  omepos = omepos + 3
!
end subroutine buildtorcrdof_do
!
!-----------------------------------------------------------------------
! Subroutine that finds the jacobian factor of the 5x5 matrix needed for
! chain closure.
!
! INPUT:  rsf - index of final residue, allows location of s and u
! OUTPUT: jac_x - that jacobian value
! ----------------------------------------------------------------------
subroutine jacobian_torcr_do(rsf,jac_x)
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
  RTYPE c6,c1,c2,c3,c4,c5,u6(3),u1(3),u2(3),u3(3),u4(3),u5(3),s1(3),s2(3)
  RTYPE s3(3),s4(3),s5(3),t0(3),cp,r51(3),r52(3),r53(3),r54(3),s0(3)
  RTYPE v1(3),v2(3),v3(3),v4(3),transform(5,5),det_a,s6(3)
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
    s6(1) = x(ci(rsf))
    s6(2) = y(ci(rsf))
    s6(3) = z(ci(rsf))
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
    s0(1) = x(ci(rsf-2))
    s0(2) = y(ci(rsf-2))
    s0(3) = z(ci(rsf-2))
!
    do k=1,3
      u1(k) = s1(k) - s0(k)
      u2(k) = s2(k) - s1(k)
      u3(k) = s3(k) - s2(k)
      u4(k) = s4(k) - s3(k)  
      u5(k) = s5(k) - s4(k)
      u6(k) = s6(k) - s5(k)
    end do
!   
    c1 = 0.0
    c2 = 0.0
    c3 = 0.0
    c4 = 0.0
    c5 = 0.0
    c6 = 0.0
    do k=1,3
      c1 = c1 + u1(k)*u1(k)
      c2 = c2 + u2(k)*u2(k)
      c3 = c3 + u3(k)*u3(k)
      c4 = c4 + u4(k)*u4(k)
      c5 = c5 + u5(k)*u5(k)
      c6 = c6 + u6(k)*u6(k)
    end do
    u1(:) = u1(:)/sqrt(c1)
    u2(:) = u2(:)/sqrt(c2)
    u3(:) = u3(:)/sqrt(c3)
    u4(:) = u4(:)/sqrt(c4)
    u5(:) = u5(:)/sqrt(c5)
    u6(:) = u6(:)/sqrt(c6)
    r51(:) = s5(:) - s1(:)
    r52(:) = s5(:) - s2(:)
    r53(:) = s5(:) - s3(:)
    r54(:) = s5(:) - s4(:)
    call crossprod(u1,r51,v1,cp)
    call crossprod(u2,r52,v2,cp)
    call crossprod(u3,r53,v3,cp)
    call crossprod(u4,r54,v4,cp)
    do k=1,3
      transform(k,1) = v1(k)
      transform(k,2) = v2(k)
      transform(k,3) = v3(k)
      transform(k,4) = v4(k)
      transform(k,5) = 0.0
    end do
!
    call crossprod(u1,u6,t0,cp)
    transform(4,1) = t0(1)
    transform(5,1) = t0(2)
    call crossprod(u2,u6,t0,cp)
    transform(4,2) = t0(1)
    transform(5,2) = t0(2)
    call crossprod(u3,u6,t0,cp)
    transform(4,3) = t0(1)
    transform(5,3) = t0(2)
    call crossprod(u4,u6,t0,cp)
    transform(4,4) = t0(1)
    transform(5,4) = t0(2)
    call crossprod(u5,u6,t0,cp)
    transform(4,5) = t0(1)
    transform(5,5) = t0(2)
!
    call dtrm(transform,5,det_a,indx)
!
    jac_x = 1.0/abs(det_a)
!
end subroutine jacobian_torcr_do
!
! -----------------------------------------------------------------
subroutine torcr_refvals_do(rsf,ccmat,oldvals2)
!------------------------------------------------------------
! subroutine to store values needed for chain closure
! INPUTS:
!     rsf -  reference residue for chain closure
! MODIFIES:
!     oldvals2 - {p1,p2,p3,p4,p5,p6,a1,a2,a3,a4,a5,a6}
!     ccmat - {w1,w2,w3,w4,w5,w6}
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
  RTYPE ccmat(6),oldvals2(12),getblen,getztor,getbang
!
! get old bond lengths
  oldvals2(1) = getblen(ci(rsf-2),ni(rsf-1))
  oldvals2(2) = getblen(ni(rsf-1),cai(rsf-1))
  oldvals2(3) = getblen(cai(rsf-1),ci(rsf-1))
  oldvals2(4) = getblen(ci(rsf-1),ni(rsf))
  oldvals2(5) = getblen(ni(rsf),cai(rsf))
  oldvals2(6) = getblen(cai(rsf),ci(rsf))
!
! get old bond angles
  oldvals2(7) = PI - getbang(ci(rsf-2),ni(rsf-1),cai(rsf-1))/RADIAN
  oldvals2(8) = PI - getbang(ni(rsf-1),cai(rsf-1),ci(rsf-1))/RADIAN
  oldvals2(9) = PI - getbang(cai(rsf-1),ci(rsf-1),ni(rsf))/RADIAN
  oldvals2(10) = PI - getbang(ci(rsf-1),ni(rsf),cai(rsf))/RADIAN
  oldvals2(11) = PI - getbang(ni(rsf),cai(rsf),ci(rsf))/RADIAN
  oldvals2(12) = PI - getbang(cai(rsf),ci(rsf),oi(rsf))/RADIAN
!
! get old torsions
  ccmat(1) = getztor(cai(rsf-2),ci(rsf-2),ni(rsf-1),cai(rsf-1))/RADIAN
  ccmat(2) = getztor(ci(rsf-2),ni(rsf-1),cai(rsf-1),ci(rsf-1))/RADIAN
  ccmat(3) = getztor(ni(rsf-1),cai(rsf-1),ci(rsf-1),ni(rsf))/RADIAN
  ccmat(4) = getztor(cai(rsf-1),ci(rsf-1),ni(rsf),cai(rsf))/RADIAN
  ccmat(5) = getztor(ci(rsf-1),ni(rsf),cai(rsf),ci(rsf))/RADIAN
  ccmat(6) = getztor(ni(rsf),cai(rsf),ci(rsf),oi(rsf))/RADIAN
!
end
!
!---------------------------------------------------------------------
subroutine picksolu_do(rsi,rsf,dof,dphi,solmat,curpks,solucount,solcnt2,jac_b_vec,dpfmat,j,prob_move,prob_back)
!------------------------------------------------------------
!
  use iounit
  use ujglobals
  use sequen
!
  implicit none
!
  integer rsi,rsf,dof,solucount,afive,aone,j,solcnt2,rs,k,atwo
  RTYPE dpfmat(MAXUJDOF+6),dphi(MAXUJDOF+6),random,jac_b_vec(4*MAXSOLU),curpt(25)
  RTYPE solmat(MAXSOLU*4,6),dummy(2),prob_move,prob_back,pdofs(7),curpks(4*MAXSOLU,7,3)
!
  afive = 5
  aone = 1
  atwo = 2
!
  if (torcrmode.eq.2) then
    j = random()*solucount + 1
    dpfmat(1:dof) = dphi(1:dof)
    dpfmat(dof+1:dof+6) = solmat(j,1:6)
!   restore, then move entire internal chain
    call mcmove_torcr_do(rsi,rsf,dof,dpfmat,afive)
    call mcmove_torcr_do(rsi,rsf,dof,dpfmat,aone)
!   if necessary, add re-puckered sidechains
    k = 0
    do rs=rsf-1,rsf
      if (seqflag(rs).eq.5) then
        k = k + 1
        curpt(1:7) = curpks(j,1:7,k)
        call sample_pucker(rs,pdofs,curpt,atwo)
      end if
    end do
    return
  end if
!
! we get into trouble if prerotation bias is numerically tiny but other
! set is empty (should nonetheless never happen)
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
      call mcmove_torcr_do(rsi,rsf,dof,dpfmat,afive)
      call mcmove_torcr_do(rsi,rsf,dof,dpfmat,aone)
      k = 0
      do rs=rsf-1,rsf
        if (seqflag(rs).eq.5) then
          k = k + 1
          curpt(1:7) = curpks(j,1:7,k)
          call sample_pucker(rs,pdofs,curpt,atwo)
        end if
      end do
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
!  j = random()*(solucount+solcnt2) + 1
!  if (j.gt.solucount) then
!    dpfmat(1:dof) = 0.0
!  else
!    dpfmat(1:dof) = dphi(1:dof)
!  end if
!  dpfmat(dof+1:dof+6) = solmat(j,1:6)
!! restore, then move entire internal chain
!  call mcmove_torcr_do(rsi,rsf,dof,dpfmat,afive)
!  call mcmove_torcr_do(rsi,rsf,dof,dpfmat,aone)
!
  end
!
!
!---------------------------------------------------------------------
subroutine torcrchainclose_do(rsi,rsf,dof,dphi,ccmat,incrmat,oldvals2,solucount,jac_b_vec)
!------------------------------------------------------------
! Subroutine to calculate the constraints for chain closure
! for UJ prerotation.
! INPUTS:
!     rsi,rsf -  reference residues for chain closure
!     oldvals2 - reference bond lengths and angles
!     ccmat - reference torsion angle values
!     dphi - prerotation displacement vector
!     dof - degrees of freedom in pre-rotation
! OUTPUT:
!     solmat - array of change in torsions needed for closure (1:solucount,1:6)
!     solucount - number of solutions
!     jac_b_vec - vector of jacobians for solutions
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
  use zmatrix
!
  implicit none
!
  integer j,k,kk,rsi,rsf,solucount,dof,aone,athree,afive
  logical logi1,logi1a,logi1b,goodsol,atrue,really,really2,really3!slog(3)
  integer findsolu(4),kkk,jj,nintvs(4),jjj,npreintvs,atwo!,dums(20)
  RTYPE w,w1,ccmat(6),dpfmat(MAXUJDOF+6),oldvals2(12)
  RTYPE p1,p2,p3,p4,p5,p6,nopl,u(3),con1(3),posCI_0(3),dphi(MAXUJDOF+6)
  RTYPE q0(3),posCA(3),posNI(3),posNF(3),posCIF(3),posCA_0(3),w1mat(4,MAXSOLU)
  RTYPE q1(3),q2(3),R3(3,3),R2(3,3),R1(3,3),T2(3,3),invT3(3,3)
  RTYPE invR1(3,3),invR2(3,3),invR3(3,3),invT1(3,3),invT2(3,3)
  RTYPE R4(3,3),invR4(3,3),T4(3,3),invT4(3,3)
  RTYPE T1(3,3),gn,g1,intercpt,newccmat2(4,6)
  RTYPE a3,ax,s(3),T3(3,3),epsilon,midpoint,ltestw1,htestw1
  RTYPE con1xy(3),a1,ay,Ra1(3,3),Ra2(3,3),Ra3(3,3),az,increment,gnpr,gnpo,preintvbs(MAXSOLU*4,2)
  RTYPE newxaxis(3),newyaxis(3),newzaxis(3),yaxis(3),zaxis(3),a2,oldgs(4),oldws(4),maxincr
  RTYPE xaxis(3),RaE(3,3),or_pl(3),or_pl2(3),lastslope(4),curslope(4)
  RTYPE con0(3),origin(3),s1(3),searchstp,con6(3),jac_b_vec(MAXSOLU*4),blas(3)
  RTYPE invmat(3,3),u1,invRaE(3,3),incrmat(MAXSOLU*4,6),intvalbnds(4,MAXSOLU*4,2)
  RTYPE sinov(12),cosov(12),q1sq,q2sq,svq(3),w4denom,nouse(17),newt
!
! test for closure and prerotation
  blas(1) = x(ci(rsf))
  blas(2) = y(ci(rsf))
  blas(3) = z(ci(rsf))
!
! helpers
  aone = 1
  atwo = 2
  athree = 3
  afive = 5
  atrue = .true.
!
! get the current position of atoms used in chain closure
  posCA_0(1) = x(cai(rsf-2))
  posCA_0(2) = y(cai(rsf-2))
  posCA_0(3) = z(cai(rsf-2))
  posCI_0(1) = x(ci(rsf-2))
  posCI_0(2) = y(ci(rsf-2))
  posCI_0(3) = z(ci(rsf-2))
  posNI(1) = x(ni(rsf-1))
  posNI(2) = y(ni(rsf-1))
  posNI(3) = z(ni(rsf-1))
  posCA(1) = x(cai(rsf))
  posCA(2) = y(cai(rsf))
  posCA(3) = z(cai(rsf))
  posCIF(1) = x(ci(rsf))
  posCIF(2) = y(ci(rsf))
  posCIF(3) = z(ci(rsf)) 
  posNF(1) = x(oi(rsf))
  posNF(2) = y(oi(rsf))
  posNF(3) = z(oi(rsf))
!  posCAF(1) = x(cai(rsf))
!  posCAF(2) = y(cai(rsf))
!  posCAF(3) = z(cai(rsf))
!
  con0(:) = posCA_0(:)-posCI_0(:)
  con1(:) = posNI(:)-posCI_0(:)
  con6(:) = posCIF(:)-posCA(:)
!
! get bond lengths
  p1 = oldvals2(1)
  p2 = oldvals2(2)
  p3 = oldvals2(3)
  p4 = oldvals2(4)
  p5 = oldvals2(5)
  p6 = oldvals2(6)
  call genT(oldvals2(7),T1,invT1)
  call genT(oldvals2(9),T3,invT3)
  sinov(7) = sin(oldvals2(7))
  cosov(7) = cos(oldvals2(7))
  sinov(8) = sin(oldvals2(8))
  cosov(8) = cos(oldvals2(8))
  sinov(9) = sin(oldvals2(9))
  cosov(9) = cos(oldvals2(9))
  sinov(11) = sin(oldvals2(11))
  cosov(11) = cos(oldvals2(11))

! build Euler rotation matrix
  origin(1) = 0.0
  origin(2) = 0.0
  origin(3) = 0.0
  xaxis = (/ 1, 0, 0 /)
  yaxis = (/ 0, 1, 0 /)
  zaxis = (/ 0, 0, 1 /)

! projection of con1 on xy-plane
  con1xy(1) = posNI(1) - posCI_0(1)
  con1xy(2) = posNI(2) - posCI_0(2)
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
    
! set q0 = p1
  q0(1) = p1
  q0(2) = 0.0
  q0(3) = 0.0

! get q1 an q2
  call genT(oldvals2(8),T2,invT2)
  q1(1) = T2(1,1)*p3 + p2
  q1(2) = T2(2,1)*p3
  q1(3) = T2(3,1)*p3
  q1sq = sum(q1(:)*q1(:))
  call genT(oldvals2(10),T4,invT4)
  q2(1) = T4(1,1)*p5 + p4
  q2(2) = T4(2,1)*p5
  q2(3) = T4(3,1)*p5
  q2sq = sum(q2(:)*q2(:))
  w4denom = (q2(2)*q2(2)+q2(3)*q2(3))*sinov(9)

! get pos S = RaE.posNF
  s1(1) = posCA(1)-posCI_0(1)
  s1(2) = posCA(2)-posCI_0(2)
  s1(3) = posCA(3)-posCI_0(3)
  s(1) = sum(RaE(1,:)*s1(:))
  s(2) = sum(RaE(2,:)*s1(:))
  s(3) = sum(RaE(3,:)*s1(:))
  svq(:) = s(:) - q0(:)

! find u
  u(1) = sum(RaE(1,:)*con6(:))
  u(2) = sum(RaE(2,:)*con6(:))
  u(3) = sum(RaE(3,:)*con6(:))
  u1 = sqrt(sum(u(:)*u(:)))
  u(:) = u(:)/u1

! propagate search variable w1 via Newton to identify branch-independent R-solution space 
  w1mat(:,:) = 0.0
  findsolu(:) = 0
!
  kk = 0
  w1 = -PI
  g1 = 0.0
  increment = 0.0
  npreintvs = 0
  logi1b = .false.
  do while (1.eq.1)
    oldws(1) = w1
    oldgs(1) = g1
    w1 = w1 + increment
    kk = kk + 1
    if (w1.gt.PI) exit
!   curslope(1) is analytical derivative dg1/dw1: g1 > 0 means real solution
    call checkreal1_do(w1,logi1,g1,curslope(1))
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
        increment = min(UJ_params(4)/RADIAN,0.8*newt)
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
        increment = min(UJ_params(4)/RADIAN,0.8*newt)
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
    call mcmove_torcr_do(rsi,rsf,dof,dphi,athree)
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
        call checkreal3_do(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),&
 &nouse(4),nouse(5:7),2*kkk,atrue) ! branch 2 and 1 are identical, as are 3 and 4
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
          increment = min(UJ_params(4)/RADIAN,0.5*newt,maxincr)
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
          increment = min(UJ_params(4)/RADIAN,0.8*newt,maxincr)
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
!    write(*,*) kkk,nintvs(kkk),kk
    if ((logi1b.EQV..false.).AND.(logi1.EQV..true.)) then
!     this shouldn't really happen
      call checkreal1_do(oldws(1),logi1,g1,curslope(1))
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
            call checkreal3_do(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2_do(w1,newccmat2,g1,jj,really)
!             remember that checkreal3_do is called by get_newccmat2_do again -> really-check redundant
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
        intvalbnds(kkk,j,1) = intvalbnds(kkk,j,1) + 1.0d-9
        intvalbnds(kkk,j,2) = intvalbnds(kkk,j,2) - 1.0d-9
!       check one more time
        if (abs(intvalbnds(kkk,j,1)-intvalbnds(kkk,j,2)).lt.2.0d-9) then
!         write(*,*) 'vanishing interval:',j,intvalbnds(kkk,j,1)
          do jj=2*kkk-1,2*kkk
            w1 = 0.5*(intvalbnds(kkk,j,1)+intvalbnds(kkk,j,2))
            call checkreal3_do(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2_do(w1,newccmat2,g1,jj,really)
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
!    write(*,*) 'DO',kkk,kk,' steps'
!    do jj=1,nintvs(kkk)
!    write(*,*) kkk,jj,RADIAN*intvalbnds(kkk,jj,1:2)
!    end do
  end do
!
! and lastly scan them
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
        oldws(kkk) = w1
        if (w1.ge.intvalbnds(jj,j,2)) exit
        if ((w1.lt.intvalbnds(jj,j,2)).AND.(w1+increment.gt.intvalbnds(jj,j,2))) then
          w1 = intvalbnds(jj,j,2)
        else
          w1 = w1 + increment
        end if
        oldgs(1) = g1
        call get_newccmat2_do(w1,newccmat2,g1,kkk,really)
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
          call get_newccmat2_do(htestw1,newccmat2,gnpo,kkk,really2)
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
          call get_newccmat2_do(ltestw1,newccmat2,gnpr,kkk,really3)
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
          call get_newccmat2_do(ltestw1,newccmat2,gnpr,kkk,really3)
!         bail out if we nail it
          if ((really3.EQV..true.).AND.(gnpr.eq.0.0)) then
!          write(*,*) 'nailed it pr',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = ltestw1
            increment = 1.0d-9
            cycle
          end if
          call get_newccmat2_do(htestw1,newccmat2,gnpo,kkk,really2)
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
          do_wrncnt(2) = do_wrncnt(2) + 1
          if (do_wrncnt(2).eq.do_wrnlmt(2)) then
            write(ilog,*) 'Warning. Encountered interval inconsistency in loop closure of Dinner-Ulmschneider&
 & algorithm with omega sampling. This means a number was tested which was previously found to give a real solu&
 &tion and now produced NaN. This is ultimately a bug the causes of which are currently poorly understood.'
            write(ilog,*) 'This was warning #',do_wrncnt(2),' of this type not all of which may be displayed.'
            if (10.0*do_wrnlmt(2).gt.0.5*HUGE(do_wrnlmt(2))) then
              do_wrncnt(2) = 0 ! reset
            else
              do_wrnlmt(2) = do_wrnlmt(2)*10
            end if
          end if
          exit
        end if
!       we can bisect the search interval straight away
        if (((gnpr.lt.0.0).AND.(gnpo.gt.0.0)).OR.((gnpr.gt.0.0).AND.(gnpo.lt.0.0))) then
!          write(*,*) 'bisecting ',ltestw1*RADIAN,' and ',htestw1*RADIAN,' in ',kkk
          call bisection2_do(ltestw1,htestw1,kkk,findsolu(kkk))
          increment = 2.0d-9
          g1 = 0.0 ! sets oldgs to zero
          cycle
        end if
!       check for accidental crossover
        if (((oldgs(kkk).lt.0.0).AND.(g1.gt.0.0)).OR.((oldgs(kkk).gt.0.0).AND.(g1.lt.0.0))) then ! yes
!          write(*,*) 'crossover: ',oldws(kkk)*RADIAN,' and ',w1*RADIAN,' in ',kkk
          call bisection2_do(oldws(kkk),w1,kkk,findsolu(kkk))
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
! note that we filter solutions for two things: exactness of closure and NaN
! this incapacitates any outside checks of this type
!
  if (solucount.eq.0) then    
!    just exit if we found no solution
!    write(ilog,*) 'Did not find a solution to the chain closure'
    call mcmove_torcr_do(rsi,rsf,dof,dphi,athree)
    return
  else if (solucount.eq.1) then
!   for a single solution it's easy: assign, restore, perturb, get jacobian, exit
    dpfmat(1:dof) = dphi(1:dof)
    do j=1,4
      if(findsolu(j).gt.0) then
        call get_newccmat2_do(w1mat(j,1),newccmat2,nouse(1),j,really)
        goodsol = .true.
        do jjj=1,6
          if ((newccmat2(j,jjj).lt.-PI).OR.(newccmat2(j,jjj).gt.PI).OR.(really.EQV..false.)) then
            goodsol = .false.
          end if
        end do
        if (goodsol.EQV..false.) then
          solucount = solucount - 1
          call mcmove_torcr_do(rsi,rsf,dof,dphi,athree)
          return
        end if
        incrmat(1,1:6) = newccmat2(j,1:6) - ccmat(1:6)
        exit
      end if
    end do

    dpfmat(dof+1:dof+6) = incrmat(1,1:6)
!
!   restore pre-rotation segment (initial state fully restored)
    call mcmove_torcr_do(rsi,rsf,dof,dphi,athree)
!
!   move entire internal chain
    call mcmove_torcr_do(rsi,rsf,dof,dpfmat,aone)
    if ((abs(x(ci(rsf))-blas(1))+abs(y(ci(rsf))-blas(2))+abs(z(ci(rsf))-blas(3))).gt.1.0e-4) then
      solucount = solucount - 1
      call mcmove_torcr_do(rsi,rsf,dof,dpfmat,atwo)
      return
    end if
!
!   get jacobian for new state 
    call jacobian_torcr_do(rsf,jac_b_vec(1))
!    write(*,*) '1 sol. - jacobian: ',jac_b_vec(1)
    call mcmove_torcr_do(rsi,rsf,dof,dpfmat,atwo)
!
  else
!   for multiple solutions it's more work: for all: assign, restore, perturb, get jacobian, then pick according to latter
!
    dpfmat(1:dof) = dphi(1:dof)
!
    kk = 0
    logi1 = .false.
    do j=1,4
!
      do kkk=1,findsolu(j)
!
        kk = kk + 1
        call get_newccmat2_do(w1mat(j,kkk),newccmat2,nouse(1),j,really)
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
        incrmat(kk,1:6) = newccmat2(j,1:6) - ccmat(1:6)
        dpfmat(dof+1:dof+6) = incrmat(kk,1:6)
!
!       restore everything appropriately
        if (logi1.EQV..false.) then
          call mcmove_torcr_do(rsi,rsf,dof,dphi,athree)
        else
          call mcmove_torcr_do(rsi,rsf,dof,dpfmat,afive)
        end if
        logi1 = .true.
!
!       move entire internal chain
        call mcmove_torcr_do(rsi,rsf,dof,dpfmat,aone)
        if ((abs(x(ci(rsf))-blas(1))+abs(y(ci(rsf))-blas(2))+abs(z(ci(rsf))-blas(3))).gt.1.0e-4) then
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
        call jacobian_torcr_do(rsf,jac_b_vec(kk))
!
      end do
!  
    end do
!
    if (logi1.EQV..false.) then
      call mcmove_torcr_do(rsi,rsf,dof,dphi,athree)
    else
      call mcmove_torcr_do(rsi,rsf,dof,dpfmat,atwo)
    end if
!
  end if
!
  CONTAINS
!
! -------------------------------------------------------------------  
!
! r(w1) = invT1*invR1*(s-q0)
subroutine r_of_w1w2_do(r,w1val,rnorm)
!
  RTYPE, INTENT(OUT):: r(3),rnorm
  RTYPE, INTENT(IN):: w1val
!
  call genR(w1val,R1,invR1)
  invmat = matmul(invT1,invR1)
!
  r(1) = sum(invmat(1,:)*svq(:))
  r(2) = sum(invmat(2,:)*svq(:))
  r(3) = sum(invmat(3,:)*svq(:))
!
  rnorm = sum(r(:)*r(:))
!
end subroutine r_of_w1w2_do
!
!-------------------------------------------------------------------------
!
! get arg1, arg2, and w and (on request) some derivatives
! 
subroutine deriva_w1_do(sw1,cw1,rnorm,r,arg1,arg2,darg1dw1,darg2dw1,w,drdw1,dwndw1,do_deriv)
!
  RTYPE, INTENT(IN):: r(3),sw1,cw1,rnorm
  logical, INTENT(IN):: do_deriv
  RTYPE, INTENT(OUT):: arg1,arg2,darg1dw1,darg2dw1,w,drdw1(3),dwndw1
!
  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
  arg2 = w*w*q1(2)*q1(2)
!
  if (do_deriv.EQV..true.) then
    drdw1(1) = -svq(2)*sw1*sinov(7) + svq(3)*cw1*sinov(7)
    drdw1(2) = -svq(2)*sw1*cosov(7) + svq(3)*cw1*cosov(7)
    drdw1(3) = -svq(2)*cw1 - svq(3)*sw1
    dwndw1 = (1.0/(2.0*q1(2)))*(-2.0*drdw1(1)*q1(1) + sum(2.0*r(:)*drdw1(:)))
    darg1dw1 = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw1(2) + r(3)*drdw1(3)))
    darg2dw1 = 2.0*q1(2)*q1(2)*w*dwndw1
  end if
!
end subroutine deriva_w1_do
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal1_do(w1val,isreal,howfaroff,testslope)
!
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: howfaroff,testslope
  RTYPE arg1,arg2,darg1dr,darg2dr,drdw(3),rnorm,sw1,cw1,r(3),dwndr
  logical, INTENT(OUT):: isreal
  logical atrue
!
  isreal = .false.
  atrue = .true.
!
  call r_of_w1w2_do(r,w1val,rnorm)
  sw1 = sin(w1val)
  cw1 = cos(w1val)
  call deriva_w1_do(sw1,cw1,rnorm,r,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,atrue)
!
  howfaroff = arg1 - arg2
  testslope = darg1dr - darg2dr
  if (howfaroff.gt.0.0) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
end subroutine checkreal1_do
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal3_do(w1val,isreal2,howfaroff,targ1,tslp1,tslp2,w2val,sw2,cw2,t1ex,r,t,do_deriv)
!
  logical, INTENT(IN):: do_deriv
  logical, INTENT(OUT):: isreal2
  RTYPE, INTENT(OUT):: howfaroff,targ1,tslp1,tslp2,cw2,sw2
  RTYPE w1val,arg1,arg2,darg1dr,darg2dr,dwndr,drdw(3),r(3)
  RTYPE denom,enum11,enum21,enum2,sarg1,sarg2,w2val,t1ex,denum2dr,rnorm,sw1,cw1
  RTYPE dsarg1dr,dsarg2dr,denum21dr,denum11dr,dw2valdr,dt1exdr,sarg21
  logical isreal
  integer, INTENT(IN):: t
!
  tslp2 = 0.0
  targ1 = 0.0
  isreal = .false.
!
  call r_of_w1w2_do(r,w1val,rnorm)
  sw1 = sin(w1val)
  cw1 = cos(w1val)
  call deriva_w1_do(sw1,cw1,rnorm,r,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,do_deriv)
!  drdw(1) = -svq(2)*sw1*sinov(7) + svq(3)*cw1*sinov(7)
!  drdw(2) = -svq(2)*sw1*cosov(7) + svq(3)*cw1*cosov(7)
!  drdw(3) = -svq(2)*cw1 - svq(3)*sw1
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  dwnewdr = (1.0/(2.0*q1(2)))*(-2.0*drdw(1)*q1(1) + sum(2.0*r(:)*drdw(:)))
!  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  darg1dr = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw(2) + r(3)*drdw(3)))
!  arg2 = w*w*q1(2)*q1(2)
!  darg2dr = 2.0*q1(2)*q1(2)*w*dwnewdr
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

    if ((t.eq.1).OR.(t.eq.2)) then
      enum21 = (r(3)*q1(2) - r(2)*q1(3))*w*q1(2) - (r(2)*q1(2) + r(3)*q1(3))*enum2
      enum11 = (r(2)*q1(2) + r(3)*q1(3))*w*q1(2) + (r(3)*q1(2) - r(2)*q1(3))*enum2
    else if  ((t.eq.3).OR.(t.eq.4)) then
      enum21 = (r(3)*q1(2) - r(2)*q1(3))*w*q1(2) + (r(2)*q1(2) + r(3)*q1(3))*enum2
      enum11 = (r(2)*q1(2) + r(3)*q1(3))*w*q1(2) - (r(3)*q1(2) - r(2)*q1(3))*enum2
    else
      write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
 &ain closure algorithm) must be between 1 and 4.'
      call fexit()
    end if
    sarg1 = enum11/denom
    sarg2 = enum21/denom
    w2val = atan2(sarg2,sarg1) ! or use ccmat
    sw2 = sin(w2val)
    cw2 = cos(w2val)
    sarg21 = (sarg2/sarg1)
    t1ex = T2(1,1)*(r(1)-q1(1)) + T2(2,1)*(r(2)*cw2+r(3)*sw2-q1(2))
    targ1 = q2sq*sinov(9)*sinov(9) + 2.0*q2(1)*t1ex*cosov(9) - q2(1)*q2(1) - t1ex*t1ex
!
    if (do_deriv.EQV..true.) then
      denum2dr = tslp1/(2.0*enum2)
      if ((t.eq.1).OR.(t.eq.2)) then
        denum21dr = (-q1(2)*drdw(2) - q1(3)*drdw(3))*enum2 - (r(2)*q1(2) + r(3)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(3)*q1(2) - r(2)*q1(3)) + w*(drdw(3)*q1(2) - drdw(2)*q1(3)))
        denum11dr = (q1(2)*drdw(3) - q1(3)*drdw(2))*enum2 + (r(3)*q1(2) - r(2)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(2)*q1(2) + r(3)*q1(3)) + w*(drdw(2)*q1(2) + drdw(3)*q1(3)))
      else if  ((t.eq.3).OR.(t.eq.4)) then
        denum21dr = (q1(2)*drdw(2) + q1(3)*drdw(3))*enum2 + (r(2)*q1(2) + r(3)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(3)*q1(2) - r(2)*q1(3)) + w*(drdw(3)*q1(2) - drdw(2)*q1(3)))
        denum11dr = (-q1(2)*drdw(3) + q1(3)*drdw(2))*enum2 - (r(3)*q1(2) - r(2)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(2)*q1(2) + r(3)*q1(3)) + w*(drdw(2)*q1(2) + drdw(3)*q1(3)))
      end if
      dsarg1dr = (denom*denum11dr - enum11*enum2*darg1dr)/(denom*denom)
      dsarg2dr = (denom*denum21dr - enum21*enum2*darg1dr)/(denom*denom)
      dw2valdr = (1.0/(1.0 + sarg21*sarg21))*((sarg1*dsarg2dr - dsarg1dr*sarg2)/(sarg1*sarg1))
      dt1exdr = T2(1,1)*drdw(1) + T2(2,1)*(-r(2)*sw2*dw2valdr + &
 &           drdw(2)*cw2 + drdw(3)*sw2 + r(3)*cw2*dw2valdr)
      tslp2 = 2.0*q2(1)*cosov(9)*dt1exdr - 2.0*t1ex*dt1exdr
!
    end if
!
    if (targ1.gt.0.0) then
      isreal2 = .true.
    else
      isreal2 = .false.
    end if
  end if
!  write(*,*) RADIAN*w2val!RADIAN*w1val,t,isreal,isreal2
!
end subroutine checkreal3_do
!
! -------------------------------------------------------------------  
!
! bisection(w1): note that the w1-values are unwrapped on the outside
!
subroutine bisection2_do(w1,next,t,nsolu)
!
  implicit none
!
  integer t,nsolu
  logical really,really2
  RTYPE w1,next,w1val,nextw1,grr,grr2
!
  w1val = w1
  nextw1 = next
!
  epsilon = 2.0*(10.0**(-15.0))
! ! start loop
  do while (abs(nextw1 - w1val).gt.(2.0*epsilon))
!
    midpoint = (nextw1 + w1val) / 2.0
!
!   ! find f(midpoint)
    call get_newccmat2_do(w1val,newccmat2,grr,t,really)
    call get_newccmat2_do(midpoint,newccmat2,grr2,t,really2)
    if ((really.EQV..false.).OR.(really2.EQV..false.)) then
      do_wrncnt(3) = do_wrncnt(3) + 1
      if (do_wrncnt(3).eq.do_wrnlmt(3)) then
        write(ilog,*) 'Bisection method failed for Dinner-Ulmschneider concerted rotation&
 & algorithm with omega bond sampling. This indicates significantly too coarse search&
 & settings (interval inconsistencies).'
        write(ilog,*) 'This was warning #',do_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*do_wrnlmt(3).gt.0.5*HUGE(do_wrnlmt(3))) then
          do_wrncnt(3) = 0 ! reset
        else
          do_wrnlmt(3) = do_wrnlmt(3)*10
        end if
      end if
      return
    end if
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
!  write(*,*) t,UJ_params(4),w1val,newccmat2(t,2)
  
!
end subroutine bisection2_do
!!
!! -------------------------------------------------------------------  
!!
!! ! subroutine to find w2
!subroutine get_w2_do(w1val,w2val,t)
!!
!  integer t
!  RTYPE w,arg1,arg2,w1val,w2val,rnorm,enum2,denom,r(3)
!!
!  call r_of_w1w2_do(r,w1val,rnorm)
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  denom = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  enum2 = sqrt(denom - w*w*q1(2)*q1(2))
!!     
!  if ((t.eq.1).OR.(t.eq.2)) then
!    arg1 = ((r(2)*q1(2) + r(3)*q1(3))*w*q1(2) + (r(3)*q1(2) - r(2)*q1(3))*enum2)/denom
!    arg2 = ((r(3)*q1(2) - r(2)*q1(3))*w*q1(2) - (r(2)*q1(2) + r(3)*q1(3))*enum2)/denom
!  else if ((t.eq.3).OR.(t.eq.4)) then
!    arg1 = ((r(2)*q1(2) + r(3)*q1(3))*w*q1(2) - (r(3)*q1(2) - r(2)*q1(3))*enum2)/denom
!    arg2 = ((r(3)*q1(2) - r(2)*q1(3))*w*q1(2) + (r(2)*q1(2) + r(3)*q1(3))*enum2)/denom
!  else
!    write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
! &ain closure algorithm) must be between 1 and 4.'
!    call fexit()
!  end if
!!
!  w2val = atan2(arg2,arg1) ! or use ccmat
!  call genR(w2val,R2,invR2)
!!
!end subroutine get_w2_do
! -------------------------------------------------------------------  
!
! subroutine to find w4
subroutine get_newccmat2_do(w1val,newccmat2,g,t,isreal2)
!
  integer, INTENT(IN):: t
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: newccmat2(4,6),g
  RTYPE arg1,arg2,w2val,t1ex,t1y,t1z,t2val(3),w4val,w3val,w5val,howfaroff,targ1
  RTYPE ay,pos6(3),TRTR(3),R1TR(3),T1R2x(3),T1R2(3,3),T2R3x(3),T2R3(3,3),dums(3),r(3)
  RTYPE T3R4p4(3),T3R4(3,3),ctw4,stw4,stw2,ctw2
  logical afalse,isreal2
!
  afalse = .false.
!
! first get w2
  call checkreal3_do(w1val,isreal2,howfaroff,targ1,nouse(1),nouse(2),w2val,stw2,ctw2,t1ex,r,t,afalse)
  call genR(w2val,R2,invR2)
!
  if (isreal2.EQV..false.) then
    return ! deal on the outside
  end if
!  call get_w2_do(w1val,w2val,t)
!  stw2 = sin(w2val)
!  ctw2 = cos(w2val)
!  t1x = T2(1,1)*(r(1)-q1(1)) + T2(2,1)*(r(2)*cos(w2val)+r(3)*sin(w2val)-q1(2))
!  denom = (q2(2)*q2(2)+q2(3)*q2(3))*sinov(9)
!  enum21 = q2sq*sinov(9)*sinov(9) + 2.0*q2(1)*t1x*cosov(9)
!  enum22 =  q2(1)*q2(1) + t1x*t1x
!  dums(3) = sqrt(enum21-enum22)
  dums(3) = sqrt(targ1)
!!
  if ((t.eq.1).OR.(t.eq.3)) then
    arg1 = (q2(2)*(q2(1)*cosov(9) - t1ex) + q2(3)*dums(3))/w4denom
    arg2 = (q2(3)*(-q2(1)*cosov(9) + t1ex) - q2(2)*dums(3))/w4denom
  else if ((t.eq.2).OR.(t.eq.4)) then
    arg1 = (q2(2)*(q2(1)*cosov(9) - t1ex) - q2(3)*dums(3))/w4denom
    arg2 = (q2(3)*(-q2(1)*cosov(9) + t1ex) + q2(2)*dums(3))/w4denom
  end if
!
  w4val = atan2(arg2,arg1)
  ctw4 = cos(w4val)
  stw4 = sin(w4val)
!
! -------------------------------------------------------------------  
!
  t2val(1) = q2(1)*cosov(9) - q2(2)*sinov(9)*ctw4 + q2(3)*sinov(9)*stw4
  t2val(2) = q2(1)*sinov(9) + q2(2)*cosov(9)*ctw4 - q2(3)*cosov(9)*stw4
  t2val(3) = q2(2)*stw4 + q2(3)*ctw4
! 
  t1y = -sinov(8)*(r(1)-q1(1)) + cosov(8)*(r(2)*ctw2 + r(3)*stw2 - q1(2))
  t1z = -r(2)*stw2 + r(3)*ctw2 - q1(3)
!
  arg1 = (t1y*t2val(2) + t1z*t2val(3))/(t1y*t1y + t1z*t1z)
  arg2 = (t1z*t2val(2) - t1y*t2val(3))/(t1y*t1y + t1z*t1z)
  w3val = atan2(arg2,arg1)
!
! -------------------------------------------------------------------  
!
  call genR(w3val,R3,invR3)
  call genR(w4val,R4,invR4)
  call genR(w1val,R1,invR1)
  invmat = matmul(invT4,matmul(invR4,matmul(invT3,matmul(invR3,matmul(invT2,matmul(invR2,&
 &          matmul(invT1,invR1)))))))
!
  arg1 = sum(invmat(2,:)*u(:))/sinov(11)
  arg2 = sum(invmat(3,:)*u(:))/sinov(11)
!
  w5val = atan2(arg2,arg1)
  newccmat2(t,1) = w1val
  newccmat2(t,2) = w2val
  newccmat2(t,3) = w3val
  newccmat2(t,4) = w4val
  newccmat2(t,5) = w5val
!
  T1R2 = matmul(T1,R2)
  T2R3 = matmul(T2,R3)
  T3R4 = matmul(T3,R4)
  T3R4p4(1) = T3R4(1,1)*p4 
  T3R4p4(2) = T3R4(2,1)*p4
  T3R4p4(3) = T3R4(3,1)*p4
  T3R4p4(1) = T3R4p4(1) + p3
  do k=1,3
    T2R3x(k) = sum(T2R3(k,:)*T3R4p4(:))
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
  pos6(:) = TRTR(:) + posCI_0(:)
  call dihed(pos6,posCA,posCIF,posNF,ay)
  newccmat2(t,6) = ay
!
  g = u(1)*invmat(1,1) + u(2)*invmat(1,2) + u(3)*invmat(1,3) - cosov(11)
!
end subroutine get_newccmat2_do
!!
!! -------------------------------------------------------------------  
!!
!! subroutine outputs g(w1)
!!
!subroutine g_of_w1w2_do(w1val,g,t)
!!
!  integer t
!  RTYPE w1val,g
!  RTYPE TRmat(3,3)
!!
!! using branch t get solution for g(w1)
!  call get_newccmat2_do(w1val,newccmat2,t)
!!
!! calculate g(w1)
!  TRmat = matmul(invT4,matmul(invR4,matmul(invT3,matmul(invR3,matmul&
! &(invT2,matmul(invR2,matmul(invT1,invR1)))))))
!  g = u(1)*TRmat(1,1) + u(2)*TRmat(1,2) + u(3)*TRmat(1,3) - cosov(11)
!!
!end subroutine g_of_w1w2_do
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
end subroutine torcrchainclose_do
!
!
!
!
!
subroutine buildtorcrdof_dj(rsi,rsf,dof,torcrmat,omepos,nome)
!------------------------------------------------------------
! Subroutine to calculate the matrix I = (dA/dphi)**2 for tor_cr prerotation.
! 
! INPUTS:     rsi, rsf, dof - inital and final fixed residues number, and number of pre-rot. dof.s
! MODIFIES:   torcrmat
!------------------------------------------------------------
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
  integer j,k,m,rsi,rsf,dof,rs,dofcnt,omepos(MAXUJDOF),nome
  RTYPE pos2(3),pos3(3),pos4(3)
  RTYPE torcrmat(MAXUJDOF,MAXUJDOF),doflever(3,MAXUJDOF+1),dvec(3)
!
! loop over all free residues of the prerotation
  dofcnt = 0
!  whereome = .false.
!  omepos = 0
  nome = 0
  do rs=rsf-2,rsi,-1
!
    if (rs.lt.rsf-2) then
!     psi angle
      pos4(1) = x(ni(rs))
      pos4(2) = y(ni(rs))
      pos4(3) = z(ni(rs))
      pos3(1) = x(cai(rs))
      pos3(2) = y(cai(rs))
      pos3(3) = z(cai(rs))
      pos2(1) = x(ci(rs))
      pos2(2) = y(ci(rs))
      pos2(3) = z(ci(rs))
      call get_tor_lever(molofrs(atmres(cai(rsf))),cai(rsf),pos2,pos3,pos4,dvec)
      doflever(1:3,dof-dofcnt) = dvec(:)/RADIAN

      dofcnt = dofcnt + 1
      if (dofcnt.eq.dof) exit
!
      if (seqflag(rs).ne.5) then
!       phi angle
        pos4(1) = x(ci(rs-1))
        pos4(2) = y(ci(rs-1))
        pos4(3) = z(ci(rs-1))
        pos3(1) = x(ni(rs))
        pos3(2) = y(ni(rs))
        pos3(3) = z(ni(rs))
        pos2(1) = x(cai(rs))
        pos2(2) = y(cai(rs))
        pos2(3) = z(cai(rs))
        call get_tor_lever(molofrs(atmres(cai(rsf))),cai(rsf),pos2,pos3,pos4,dvec)
        doflever(1:3,dof-dofcnt) = dvec(:)/RADIAN
        dofcnt = dofcnt + 1
        if (dofcnt.eq.dof) exit
      end if
    end if
!
!   omega angle
    pos4(1) = x(iz(3,wline(rs)))
    pos4(2) = y(iz(3,wline(rs)))
    pos4(3) = z(iz(3,wline(rs)))
    pos3(1) = x(ci(rs-1))
    pos3(2) = y(ci(rs-1))
    pos3(3) = z(ci(rs-1))
    pos2(1) = x(ni(rs))
    pos2(2) = y(ni(rs))
    pos2(3) = z(ni(rs))
    call get_tor_lever(molofrs(atmres(cai(rsf))),cai(rsf),pos2,pos3,pos4,dvec)
    doflever(1:3,dof-dofcnt) = dvec(:)/RADIAN
!    if (whereome.EQV..false.) then
!      whereome = .true.
!      omepos = dof-dofcnt
!    end if
    nome = nome + 1
    omepos(nome) = dof-dofcnt
    dofcnt = dofcnt + 1
    if (dofcnt.eq.dof) exit
!
  end do
!
! ! build symetric matrix I = dA/dphi*dA/dphi (outer_product)
  do j=1,dof
    do k=1,dof
!     ! inner product
      torcrmat(j,k) = 0.0d0
      do m=1,3
        torcrmat(j,k) = torcrmat(j,k) + doflever(m,j)*doflever(m,k)
      end do
    end do
  end do
!
!  do while (omepos.gt.0) 
!    omepos = omepos - 3
!  end do
!  omepos = omepos + 3
!
end subroutine buildtorcrdof_dj
!
!-----------------------------------------------------------------------
! Subroutine that finds the jacobian factor of the 5x5 matrix needed for
! chain closure.
!
! INPUT:  rsf - index of final residue, allows location of s and u
! OUTPUT: jac_x - that jacobian value
! ----------------------------------------------------------------------
subroutine jacobian_torcr_dj(rsf,jac_x)
!
  use polypep
  use atoms
  use math
  use iounit
  use sequen
!
  implicit none
!
  integer k,indx(5)
  RTYPE c8,c1,c2,c7,c4,c5,u8(3),u1(3),u2(3),u7(3),u4(3),u5(3),s1(3),s2(3)
  RTYPE s3(3),s4(3),s5(3),t0(3),cp,r1(3),r2(3),r4(3),r5(3),s0(3)
  RTYPE v1(3),v2(3),v4(3),v5(3),transform(5,5),det_a,s7(3),s6(3),s8(3)
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
    s8(1) = x(ci(rsf))
    s8(2) = y(ci(rsf))
    s8(3) = z(ci(rsf))
    s7(1) = x(cai(rsf))
    s7(2) = y(cai(rsf))
    s7(3) = z(cai(rsf))
    s6(1) = x(ni(rsf))
    s6(2) = y(ni(rsf))
    s6(3) = z(ni(rsf))
    s5(1) = x(ci(rsf-1))
    s5(2) = y(ci(rsf-1))
    s5(3) = z(ci(rsf-1))
    s4(1) = x(cai(rsf-1))
    s4(2) = y(cai(rsf-1))
    s4(3) = z(cai(rsf-1))
    s3(1) = x(ni(rsf-1))
    s3(2) = y(ni(rsf-1))
    s3(3) = z(ni(rsf-1))
    s2(1) = x(ci(rsf-2))
    s2(2) = y(ci(rsf-2))
    s2(3) = z(ci(rsf-2))
    s1(1) = x(cai(rsf-2))
    s1(2) = y(cai(rsf-2))
    s1(3) = z(cai(rsf-2))
    s0(1) = x(ni(rsf-2))
    s0(2) = y(ni(rsf-2))
    s0(3) = z(ni(rsf-2))
!
    do k=1,3
      u1(k) = s1(k) - s0(k)
      u2(k) = s2(k) - s1(k)
      u4(k) = s4(k) - s3(k)  
      u5(k) = s5(k) - s4(k)
      u7(k) = s7(k) - s6(k)
      u8(k) = s8(k) - s7(k)
    end do
!   
    c1 = 0.0
    c2 = 0.0
    c4 = 0.0
    c5 = 0.0
    c7 = 0.0
    c8 = 0.0
    do k=1,3
      c1 = c1 + u1(k)*u1(k)
      c2 = c2 + u2(k)*u2(k)
      c4 = c4 + u4(k)*u4(k)
      c5 = c5 + u5(k)*u5(k)
      c7 = c7 + u7(k)*u7(k)
      c8 = c8 + u8(k)*u8(k)
    end do
    u1(:) = u1(:)/sqrt(c1)
    u2(:) = u2(:)/sqrt(c2)
    u4(:) = u4(:)/sqrt(c4)
    u5(:) = u5(:)/sqrt(c5)
    u7(:) = u7(:)/sqrt(c7)
    u8(:) = u8(:)/sqrt(c8)
    r1(:) = s7(:) - s1(:)
    r2(:) = s7(:) - s2(:)
    r4(:) = s7(:) - s4(:)
    r5(:) = s7(:) - s5(:)
    call crossprod(u1,r1,v1,cp)
    call crossprod(u2,r2,v2,cp)
    call crossprod(u4,r4,v4,cp)
    call crossprod(u5,r5,v5,cp)
    do k=1,3
      transform(k,1) = v1(k)
      transform(k,2) = v2(k)
      transform(k,3) = v4(k)
      transform(k,4) = v5(k)
      transform(k,5) = 0.0
    end do
!
    call crossprod(u1,u8,t0,cp)
    transform(4,1) = t0(1)
    transform(5,1) = t0(2)
    call crossprod(u2,u8,t0,cp)
    transform(4,2) = t0(1)
    transform(5,2) = t0(2)
    call crossprod(u4,u8,t0,cp)
    transform(4,3) = t0(1)
    transform(5,3) = t0(2)
    call crossprod(u5,u8,t0,cp)
    transform(4,4) = t0(1)
    transform(5,4) = t0(2)
    call crossprod(u7,u8,t0,cp)
    transform(4,5) = t0(1)
    transform(5,5) = t0(2)
!
    call dtrm(transform,5,det_a,indx)
!
    jac_x = 1.0/abs(det_a)
!
end subroutine jacobian_torcr_dj
!
! -----------------------------------------------------------------
subroutine torcr_refvals_dj(rsf,ccmat,oldvals3)
!------------------------------------------------------------
! subroutine to store values needed for chain closure
! INPUTS:
!     rsf -  reference residue for chain closure
! MODIFIES:
!     oldvals3 - {p1,p2,p3,p4,p5,p6,a1,a2,a3,a4,a5,a6}
!     ccmat - {w1,w2,w3,w4,w5,w6}
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
  RTYPE ccmat(6),oldvals3(18),getblen,getztor,getbang
!
! get old bond lengths
  oldvals3(1) = getblen(ni(rsf-2),cai(rsf-2))
  oldvals3(2) = getblen(cai(rsf-2),ci(rsf-2))
  oldvals3(3) = getblen(ci(rsf-2),ni(rsf-1))
  oldvals3(4) = getblen(ni(rsf-1),cai(rsf-1))
  oldvals3(5) = getblen(cai(rsf-1),ci(rsf-1))
  oldvals3(6) = getblen(ci(rsf-1),ni(rsf))
  oldvals3(7) = getblen(ni(rsf),cai(rsf))
  oldvals3(8) = getblen(cai(rsf),ci(rsf))
!
! get old bond angles
  oldvals3(9) = PI - getbang(ni(rsf-2),cai(rsf-2),ci(rsf-2))/RADIAN
  oldvals3(10) = PI - getbang(cai(rsf-2),ci(rsf-2),ni(rsf-1))/RADIAN
  oldvals3(11) = PI - getbang(ci(rsf-2),ni(rsf-1),cai(rsf-1))/RADIAN
  oldvals3(12) = PI - getbang(ni(rsf-1),cai(rsf-1),ci(rsf-1))/RADIAN
  oldvals3(13) = PI - getbang(cai(rsf-1),ci(rsf-1),ni(rsf))/RADIAN
  oldvals3(14) = PI - getbang(ci(rsf-1),ni(rsf),cai(rsf))/RADIAN
  oldvals3(15) = PI - getbang(ni(rsf),cai(rsf),ci(rsf))/RADIAN
  oldvals3(16) = PI - getbang(cai(rsf),ci(rsf),oi(rsf))/RADIAN
!
! get omega torsions
  oldvals3(17) = getztor(cai(rsf-2),ci(rsf-2),ni(rsf-1),cai(rsf-1))/RADIAN
  oldvals3(18) = getztor(cai(rsf-1),ci(rsf-1),ni(rsf),cai(rsf))/RADIAN
!
! get old torsions
  ccmat(1) = getztor(ci(rsf-3),ni(rsf-2),cai(rsf-2),ci(rsf-2))/RADIAN
  ccmat(2) = getztor(ni(rsf-2),cai(rsf-2),ci(rsf-2),ni(rsf-1))/RADIAN
  ccmat(3) = getztor(ci(rsf-2),ni(rsf-1),cai(rsf-1),ci(rsf-1))/RADIAN
  ccmat(4) = getztor(ni(rsf-1),cai(rsf-1),ci(rsf-1),ni(rsf))/RADIAN
  ccmat(5) = getztor(ci(rsf-1),ni(rsf),cai(rsf),ci(rsf))/RADIAN
  ccmat(6) = getztor(ni(rsf),cai(rsf),ci(rsf),oi(rsf))/RADIAN
!
end subroutine torcr_refvals_dj
!
!---------------------------------------------------------------------
subroutine picksolu_dj(rsi,rsf,dof,dphi,solmat,curpks,solucount,solcnt2,jac_b_vec,dpfmat,j,prob_move,prob_back)
!------------------------------------------------------------
!
  use iounit
  use ujglobals
  use sequen
!
  implicit none
!
  integer rsi,rsf,dof,solucount,afive,aone,j,solcnt2,rs,k,atwo
  RTYPE dpfmat(MAXUJDOF+6),dphi(MAXUJDOF+6),random,jac_b_vec(4*MAXSOLU),curpt(25)
  RTYPE solmat(MAXSOLU*4,6),dummy(2),prob_move,prob_back,pdofs(7),curpks(4*MAXSOLU,7,3)
!
  afive = 5
  aone = 1
  atwo = 2
!
  if (torcrmode.eq.2) then
    j = random()*solucount + 1
    dpfmat(1:dof) = dphi(1:dof)
    dpfmat(dof+1:dof+6) = solmat(j,1:6)
!   restore, then move entire internal chain
    call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,afive)
    call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,aone)
!   if necessary, add re-puckered sidechains
    k = 0
    do rs=rsf-2,rsf
      if (seqflag(rs).eq.5) then
        k = k + 1
        curpt(1:7) = curpks(j,1:7,k)
        call sample_pucker(rs,pdofs,curpt,atwo)
      end if
    end do
    return
  end if
!
! we get into trouble if prerotation bias is numerically tiny but other
! set is empty (should nonetheless never happen)
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
      call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,afive)
      call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,aone)
!     if necessary, add re-puckered sidechains
      k = 0
      do rs=rsf-2,rsf
        if (seqflag(rs).eq.5) then
          k = k + 1
          curpt(1:7) = curpks(j,1:7,k)
          call sample_pucker(rs,pdofs,curpt,atwo)
        end if
      end do
      exit
    else
      if (j.eq.(solucount+solcnt2)) then
        write(ilog,*) 'Fatal. Jacobian weighting of solutions is inconsistent. This is a bug.'
!        write(*,*) solucount,solcnt2
!        write(*,*) prob_back,prob_move,jac_b_vec(1:solucount+solcnt2)
!        write(*,*) sum(prob_move*jac_b_vec(1:solucount))
!        write(*,*) sum(prob_back*jac_b_vec(solucount+1:solucount+solcnt2))
!        write(*,*) dummy(1),dummy(2)
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
subroutine torcrchainclose_dj(rsi,rsf,dof,dphi,ccmat,incrmat,oldvals3,solucount,jac_b_vec)
!------------------------------------------------------------
! Subroutine to calculate the constraints for chain closure
! for UJ prerotation.
! INPUTS:
!     rsf -  reference residue for chain closure
!     oldvals3 - reference bond lengths and angles
!     ccmat - reference torsion angle values
! OUTPUT:
!     dccmat - array of change in angles needed for closure
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
  integer j,k,kk,kkk,rsi,rsf,solucount,findsolu(4)
  integer dof,athree,aone,afive,npreintvs,nintvs(2),jj,jjj,atwo
  logical logi1,logi1b,logi1a,goodsol,atrue,really,really2,really3!,slog(3)
  RTYPE w,w1,ccmat(6),dpfmat(MAXUJDOF+6),oldvals3(18),dphi(MAXUJDOF+6)
  RTYPE p1,p2,p3,p4,p5,p6,p7,p8,nopl,u(3),con1(3),posCI_0(3),incrmat(MAXSOLU*4,6)
  RTYPE q0(3),posCA(3),posNI(3),posNF(3),posCIF(3),posCAF(3),w1mat(4,MAXSOLU)
  RTYPE q1(3),q2(3),R3(3,3),R2(3,3),R1(3,3),T2(3,3),invT3(3,3),invT6(3,3),maxincr
  RTYPE invR1(3,3),invR2(3,3),invR3(3,3),invT1(3,3),invT2(3,3),invR6(3,3),htestw1,ltestw1
  RTYPE R4(3,3),invR4(3,3),T4(3,3),invT4(3,3),TR1(3,3),TR2(3,3),intvalbnds(4,MAXSOLU*4,2)
  RTYPE T1(3,3),gn,g1,intercpt,newccmat2(4,6),T5(3,3),R6(3,3),T6(3,3)
  RTYPE a3,ax,s(3),T3(3,3),epsilon,midpoint,invT5(3,3)
  RTYPE con1xy(3),a1,ay,Ra1(3,3),Ra2(3,3),Ra3(3,3),az,invR5(3,3),R5(3,3),gnpr,gnpo
  RTYPE newxaxis(3),newyaxis(3),newzaxis(3),yaxis(3),zaxis(3),a2,preintvbs(MAXSOLU*4,2)
  RTYPE xaxis(3),RaE(3,3),or_pl(3),or_pl2(3),oldws(4),curslope(4),lastslope(4)
  RTYPE con0(3),origin(3),s1(3),searchstp,con8(3),oldgs(4),increment
  RTYPE invmat(3,3),u1,invRaE(3,3),jac_b_vec(MAXSOLU*4)
  RTYPE sinov(18),cosov(18),svq(3),q1sq,q2sq,w4denom,nouse(20),newt
!
! helpers
  aone = 1
  atwo = 2
  athree = 3
  afive = 5
  atrue = .true.
!
! get the current position of atoms used in chain closure
  posCI_0(1) = x(ci(rsf-3))
  posCI_0(2) = y(ci(rsf-3))
  posCI_0(3) = z(ci(rsf-3))
  posNI(1) = x(ni(rsf-2))
  posNI(2) = y(ni(rsf-2))
  posNI(3) = z(ni(rsf-2))
  posCA(1) = x(cai(rsf-2))
  posCA(2) = y(cai(rsf-2))
  posCA(3) = z(cai(rsf-2))
  posCAF(1) = x(cai(rsf))
  posCAF(2) = y(cai(rsf))
  posCAF(3) = z(cai(rsf)) 
  posCIF(1) = x(ci(rsf))
  posCIF(2) = y(ci(rsf))
  posCIF(3) = z(ci(rsf))
  posNF(1) = x(oi(rsf))
  posNF(2) = y(oi(rsf))
  posNF(3) = z(oi(rsf))
!
  con0(:) = posCI_0(:)-posNI(:)
  con1(:) = posCA(:)-posNI(:)
  con8(:) = posCIF(:)-posCAF(:)
!
! get bond lengths
  p1 = oldvals3(1)
  p2 = oldvals3(2)
  p3 = oldvals3(3)
  p4 = oldvals3(4)
  p5 = oldvals3(5)
  p6 = oldvals3(6)
  p7 = oldvals3(7)
  p8 = oldvals3(8)
  sinov(9) = sin(oldvals3(9))
  cosov(9) = cos(oldvals3(9))
  sinov(12) = sin(oldvals3(12))
  cosov(12) = cos(oldvals3(12))
  sinov(15) = sin(oldvals3(15))
  cosov(15) = cos(oldvals3(15))

! build Euler rotation matrix
  origin(1) = 0.0
  origin(2) = 0.0
  origin(3) = 0.0
  xaxis = (/ 1, 0, 0 /)
  yaxis = (/ 0, 1, 0 /)
  zaxis = (/ 0, 0, 1 /)

! projection of con1 on xy-plane
  con1xy(1) = posCA(1) - posNI(1)
  con1xy(2) = posCA(2) - posNI(2)
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
!  write(*,*) "RaE",RaE
  invRaE(1,1) = RaE(1,1)
  invRaE(1,2) = RaE(2,1)
  invRaE(1,3) = RaE(3,1)
  invRaE(2,1) = RaE(1,2)
  invRaE(2,2) = RaE(2,2)
  invRaE(2,3) = RaE(3,2)
  invRaE(3,1) = RaE(1,3)
  invRaE(3,2) = RaE(2,3)
  invRaE(3,3) = RaE(3,3)
    
! set q0 = p1
  q0(1) = p1
  q0(2) = 0.0
  q0(3) = 0.0
!  write(*,*) "q0",q0
  call genT(oldvals3(9),T1,invT1)
  call genT(oldvals3(10),T2,invT2)
  call genT(oldvals3(11),T3,invT3)
  call genT(oldvals3(12),T4,invT4)
  call genT(oldvals3(13),T5,invT5)
  call genT(oldvals3(14),T6,invT6)
  call genR(oldvals3(17),R3,invR3)
  call genR(oldvals3(18),R6,invR6)
  TR1 = matmul(T2,matmul(R3,T3))
  TR2 = matmul(T5,matmul(R6,T6))
! get q1 an q2
  q1(1) = T2(1,1)*p3 + TR1(1,1)*p4 + p2
  q1(2) = T2(2,1)*p3 + TR1(2,1)*p4
  q1(3) = T2(3,1)*p3 + TR1(3,1)*p4
  q1sq = sum(q1(:)*q1(:))
  q2(1) = T5(1,1)*p6 + TR2(1,1)*p7 + p5
  q2(2) = T5(2,1)*p6 + TR2(2,1)*p7
  q2(3) = T5(3,1)*p6 + TR2(3,1)*p7
  q2sq = sum(q2(:)*q2(:))
  w4denom = ((q2(2)*q2(2)+q2(3)*q2(3))*sinov(12))
!  write(*,*) "q1",q1
!  write(*,*) "q2",q2

! get pos S = RaE.posCIF
  s1(1) = posCAF(1)-posNI(1)
  s1(2) = posCAF(2)-posNI(2)
  s1(3) = posCAF(3)-posNI(3)
  s(1) = sum(RaE(1,:)*s1(:))
  s(2) = sum(RaE(2,:)*s1(:))
  s(3) = sum(RaE(3,:)*s1(:))
  svq(:) = s(:) - q0(:)
!  write(*,*) "s",s

! find u
  u(1) = sum(RaE(1,:)*con8(:))
  u(2) = sum(RaE(2,:)*con8(:))
  u(3) = sum(RaE(3,:)*con8(:))
  u1 = sqrt(sum(u(:)*u(:)))
  u(:) = u(:)/u1
!  write(*,*) "u",u

! step through w1 ranging from -scaniv degrees to +scaniv degrees of original
  findsolu(:) = 0
  newccmat2(:,:) = 0.0
!
  w1 = -PI
  g1 = 0.0
  npreintvs = 0
  increment = 0.0
  logi1b = .false.
  kk = 0
!  write(*,*) 'setup'
!
  do while (1.eq.1)
    oldws(1) = w1
    oldgs(1) = g1
    w1 = w1 + increment
    kk = kk + 1
    if (w1.gt.PI) exit
    call checkreal1_dj(w1,logi1,g1,curslope(1))
    if (curslope(1).ne.0.0) then
      newt = abs(g1/curslope(1))
    else if ((increment.ne.0.0).AND.(g1.ne.oldgs(1))) then
      newt = abs(g1*increment/(g1-oldgs(1)))
    else
      newt = 1.0e-8 ! bail
    end if
!    call checkreal1_dj(w1-1.0d-9,logi1a,gnpr,curslope(2))
!    call checkreal1_dj(w1+1.0d-9,logi1b,gnpo,curslope(3))
!    write(*,*) (gnpo-gnpr)/2.0d-9,curslope(1)
!
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
    call mcmove_torcr_dj(rsi,rsf,dof,dphi,athree)
    return
  end if
!  write(*,*) kk
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
        call checkreal3_dj(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
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
      call checkreal1_dj(oldws(1),logi1,g1,curslope(1))
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
            call checkreal3_dj(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2_dj(w1,newccmat2,g1,jj,really)
!             note that get_newccmat2_dj calls checkreal3_dj, hence the really-check is redundant
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
            call checkreal3_dj(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2_dj(w1,newccmat2,g1,jj,really)
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
!    write(*,*) kk,' steps'
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
        call get_newccmat2_dj(w1,newccmat2,g1,kkk,really)
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
          call get_newccmat2_dj(htestw1,newccmat2,gnpo,kkk,really2)
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
          call get_newccmat2_dj(ltestw1,newccmat2,gnpr,kkk,really3)
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
          call get_newccmat2_dj(ltestw1,newccmat2,gnpr,kkk,really3)
!         bail out if we nail it
          if ((really3.EQV..true.).AND.(gnpr.eq.0.0)) then
!          write(*,*) 'nailed it pr',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = ltestw1
            increment = 1.0d-9
            cycle
          end if
          call get_newccmat2_dj(htestw1,newccmat2,gnpo,kkk,really2)
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
          dj_wrncnt(2) = dj_wrncnt(2) + 1
          if (dj_wrncnt(2).eq.dj_wrnlmt(2)) then
            write(ilog,*) 'Warning. Encountered interval inconsistency in loop closure of Dinner-Ulmschneider&
 & algorithm without omega sampling. This means a number was tested which was previously found to give a real &
 &solution and now produced NaN. This is ultimately a bug the causes of which are currently poorly understood.'
            write(ilog,*) 'This was warning #',dj_wrncnt(2),' of this type not all of which may be displayed.'
            if (10.0*dj_wrnlmt(2).gt.0.5*HUGE(dj_wrnlmt(2))) then
              dj_wrncnt(2) = 0 ! reset
            else
              dj_wrnlmt(2) = dj_wrnlmt(2)*10
            end if
          end if
          exit
        end if
!       we can bisect the search interval straight away
        if (((gnpr.lt.0.0).AND.(gnpo.gt.0.0)).OR.((gnpr.gt.0.0).AND.(gnpo.lt.0.0))) then
!          write(*,*) 'bisecting ',ltestw1*RADIAN,' and ',htestw1*RADIAN,' in ',kkk
          call bisection2_dj(ltestw1,htestw1,kkk,findsolu(kkk))
!          findsolu(kkk) = findsolu(kkk) + 1
          increment = 2.0d-9
          g1 = 0.0 ! sets oldgs to zero
          cycle
        end if
!       check for accidental crossover
        if (((oldgs(kkk).lt.0.0).AND.(g1.gt.0.0)).OR.((oldgs(kkk).gt.0.0).AND.(g1.lt.0.0))) then ! yes
!          write(*,*) 'crossover: ',oldws(kkk)*RADIAN,' and ',w1*RADIAN,' in ',kkk
          call bisection2_dj(oldws(kkk),w1,kkk,findsolu(kkk))
!          findsolu(kkk) = findsolu(kkk) + 1
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
! below is the old incremental step search with bisection
!
!  searchstp = UJ_params(4)
!  scaniv = UJ_params(5)
!  w1_0 = ccmat(1) - (scaniv+searchstp)/radian
!  nstps = 2.0*ceiling(scaniv/searchstp)
!  if (w1_0.lt.-PI) w1_0 = w1_0 + 2.0*PI
!
!  do j=1,nstps
!    w1 = w1_0 + searchstp*j/radian
!    if (w1.gt.PI) w1 = w1 - 2.0*PI
!    next = w1 +searchstp/radian
!    if (next.gt.PI) next = next - 2.0*PI
!!
!    call checkreal1(w1,logi1)
!    call checkreal1(next,logi2)
!    if ((logi1.EQV..true.).AND.(logi2.EQV..true.)) then
!      call checkreal2(w1,logi3,1)
!      call checkreal2(next,logi4,1)
!      if ((logi3.EQV..true.).AND.(logi4.EQV..true.)) then
!        call g_of_w1w2(w1,g1,1)
!        call g_of_w1w2(next,gn,1)
!        intercpt = gn*g1
!        if(intercpt.lt.0.0) then
!          call bisection2(w1,next,1,findsolu(1)+1)
!          findsolu(1) = findsolu(1) + 1
!        end if
!        call g_of_w1w2(w1,g1,2)
!        call g_of_w1w2(next,gn,2)
!        intercpt = gn*g1
!        if(intercpt.lt.0.0) then
!          call bisection2(w1,next,2,findsolu(2)+1)
!          findsolu(2) = findsolu(2) + 1
!        end if
!      end if
!      call checkreal2(w1,logi3,3)
!      call checkreal2(next,logi4,3)
!      if ((logi3.EQV..true.).AND.(logi4.EQV..true.)) then
!        call g_of_w1w2(w1,g1,3)
!        call g_of_w1w2(next,gn,3)
!        intercpt = gn*g1
!        if(intercpt.lt.0.0) then
!          call bisection2(w1,next,3,findsolu(3)+1)
!          findsolu(3) = findsolu(3) + 1
!        end if
!        call g_of_w1w2(w1,g1,4)
!        call g_of_w1w2(next,gn,4)
!        intercpt = gn*g1
!        if(intercpt.lt.0.0) then
!          call bisection2(w1,next,4,findsolu(4)+1)
!          findsolu(4) = findsolu(4) + 1
!        end if
!      end if
!    end if
!  end do
!
  solucount = sum(findsolu(1:4))
!
  if (solucount.eq.0) then    
!    just exit if we found no solution
     call mcmove_torcr_dj(rsi,rsf,dof,dphi,athree)
!    write(ilog,*) 'Did not find a solution to the chain closure'
  else if (solucount.eq.1) then
!   for a single solution it's easy: assign, restore, perturb, get jacobian, exit
    dpfmat(1:dof) = dphi(1:dof)
    do j=1,4
      if(findsolu(j).gt.0) then
        call get_newccmat2_dj(w1mat(j,1),newccmat2,g1,j,really)
        goodsol = .true.
        do jjj=1,6
          if ((newccmat2(j,jjj).lt.-PI).OR.(newccmat2(j,jjj).gt.PI).OR.(really.EQV..false.)) then
            goodsol = .false.
          end if
        end do
        if (goodsol.EQV..false.) then
          solucount = solucount - 1
          call mcmove_torcr_dj(rsi,rsf,dof,dphi,athree)
          return
        end if
        incrmat(1,1:6) = newccmat2(j,1:6) - ccmat(1:6)
        exit
      end if
    end do
    dpfmat(dof+1:dof+6) = incrmat(1,1:6)
!
!   restore pre-rotation segment (initial state fully restored)
    call mcmove_torcr_dj(rsi,rsf,dof,dphi,athree)
!
!   move entire internal chain
    call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,aone)
    if ((abs(x(ci(rsf))-posCIF(1))+abs(y(ci(rsf))-posCIF(2))+abs(z(ci(rsf))-posCIF(3))).gt.1.0e-4) then
      solucount = solucount - 1
      call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,atwo)
      return
    end if
!
!   get jacobian for new state 
    call jacobian_torcr_dj(rsf,jac_b_vec(1))
!   restore
    call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,atwo)
!
  else
!   for multiple solutions it's more work: for all: assign, restore, perturb, get jacobian, then pick according to latter
!
    dpfmat(1:dof) = dphi(1:dof)
!
    kk = 0
    logi1 = .false.
    do j=1,4
!
      do kkk=1,findsolu(j)
        kk = kk + 1
        call get_newccmat2_dj(w1mat(j,kkk),newccmat2,g1,j,really)
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
        incrmat(kk,1:6) = newccmat2(j,1:6) - ccmat(1:6)
!
        dpfmat(dof+1:dof+6) = incrmat(kk,1:6)
!
!       restore everything appropriately
        if (logi1.EQV..false.) then
          call mcmove_torcr_dj(rsi,rsf,dof,dphi,athree)
        else
          call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,afive)
        end if
!
        logi1 = .true.
!       move entire internal chain
        call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,aone)
        if ((abs(x(ci(rsf))-posCIF(1))+abs(y(ci(rsf))-posCIF(2))+abs(z(ci(rsf))-posCIF(3))).gt.1.0e-4) then
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
        call jacobian_torcr_dj(rsf,jac_b_vec(kk))
!
      end do
      
    end do
!
    if (logi1.EQV..false.) then
      call mcmove_torcr_dj(rsi,rsf,dof,dphi,athree)
    else
      call mcmove_torcr_dj(rsi,rsf,dof,dpfmat,atwo)
    end if
!
  end if
!
  CONTAINS
! -------------------------------------------------------------------  
!
! r(w1) = invT1*invR1*(s-q0)
subroutine r_of_w1w2_dj(r,w1val,rnorm)
!
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: r(3),rnorm
!
  call genR(w1val,R1,invR1)
  invmat = matmul(invT1,invR1)
!
  r(1) = sum(invmat(1,:)*svq(:))
  r(2) = sum(invmat(2,:)*svq(:))
  r(3) = sum(invmat(3,:)*svq(:))
!
  rnorm = sum(r(:)*r(:))
!
!
end subroutine r_of_w1w2_dj
!
! -------------------------------------------------------------------  
!
subroutine deriva_w1_dj(sw1,cw1,rnorm,r,arg1,arg2,darg1dw1,darg2dw1,w,drdw1,dwndw1,do_deriv)
!
  RTYPE, INTENT(IN):: r(3),sw1,cw1,rnorm
  logical, INTENT(IN):: do_deriv
  RTYPE, INTENT(OUT):: arg1,arg2,darg1dw1,darg2dw1,w,drdw1(3),dwndw1
!
  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
  arg2 = w*w*q1(2)*q1(2)
!
  if (do_deriv.EQV..true.) then
    drdw1(1) = -svq(2)*sw1*sinov(9) + svq(3)*cw1*sinov(9)
    drdw1(2) = -svq(2)*sw1*cosov(9) + svq(3)*cw1*cosov(9)
    drdw1(3) = -svq(2)*cw1 - svq(3)*sw1
    dwndw1 = (1.0/(2.0*q1(2)))*(-2.0*drdw1(1)*q1(1) + sum(2.0*r(:)*drdw1(:)))
    darg1dw1 = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw1(2) + r(3)*drdw1(3)))
    darg2dw1 = 2.0*q1(2)*q1(2)*w*dwndw1
  end if
!
end subroutine deriva_w1_dj
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal1_dj(w1val,isreal,howfaroff,testslope)
!
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: howfaroff,testslope
  RTYPE sw1,cw1,drdw(3),darg1dr,darg2dr,w,dwndr,rnorm,arg1,arg2,r(3)
  logical, INTENT(OUT):: isreal
  logical atrue
!
  isreal = .false.
  atrue = .true.
!
  call r_of_w1w2_dj(r,w1val,rnorm)
!
  sw1 = sin(w1val)
  cw1 = cos(w1val)
!
  call deriva_w1_dj(sw1,cw1,rnorm,r,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,atrue)
!  drdw(1) = -svq(2)*sw1*sinov(9) + svq(3)*cw1*sinov(9)
!  drdw(2) = -svq(2)*sw1*cosov(9) + svq(3)*cw1*cosov(9)
!  drdw(3) = -svq(2)*cw1 - svq(3)*sw1
!
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  dwnewdr = (1.0/(2.0*q1(2)))*(-2.0*drdw(1)*q1(1) + sum(2.0*r(:)*drdw(:)))
!
!  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  darg1dr = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw(2) + r(3)*drdw(3)))
!  arg2 = w*w*q1(2)*q1(2)
!  darg2dr = 2.0*q1(2)*q1(2)*w*dwnewdr
!
  howfaroff = arg1 - arg2
  testslope = darg1dr - darg2dr
  if (howfaroff.gt.0.0) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
end subroutine checkreal1_dj
!
! -------------------------------------------------------------------
!
subroutine checkreal3_dj(w1val,isreal2,howfaroff,targ1,tslp1,tslp2,w2val,sw2,cw2,t1ex,r,t,do_deriv)
!
  logical, INTENT(IN):: do_deriv
  logical, INTENT(OUT):: isreal2
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: howfaroff,targ1,tslp1,tslp2,cw2,sw2
  RTYPE arg1,arg2,darg1dr,darg2dr,dwndr,drdw(3),r(3),denumadr,denumbdr,drrdr(3)
  RTYPE denom,enum11,enum21,enum2,sarg1,sarg2,w2val,t1ex,denum2dr,rnorm,sw1,cw1,enuma,enumb
  RTYPE dsarg1dr,dsarg2dr,denum21dr,denum11dr,dw2valdr,dt1exdr,sarg21,t1vp(3,3),rr(3)
  logical isreal
  integer, INTENT(IN):: t
!
  tslp2 = 0.0
  targ1 = 0.0
  isreal = .false.
!
  call r_of_w1w2_dj(r,w1val,rnorm)
  sw1 = sin(w1val)
  cw1 = cos(w1val)
  call deriva_w1_dj(sw1,cw1,rnorm,r,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,do_deriv)
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
!
    denom = arg1
    enum2 = sqrt(denom - arg2)
    enuma = r(2)*q1(2) + r(3)*q1(3)
    enumb = r(3)*q1(2) - r(2)*q1(3)
!     
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
    w2val = atan2(sarg2,sarg1)
    sw2 = sin(w2val)
    cw2 = cos(w2val)
    sarg21 = (sarg2/sarg1)
    call genR(w2val,R2,invR2)
!
    t1vp = matmul(invT3,matmul(invR3,invT2))
    rr = matmul(invR2,r)
    t1ex = t1vp(1,1)*(rr(1)-q1(1)) + t1vp(1,2)*(rr(2)-q1(2)) + t1vp(1,3)*(rr(3)-q1(3))
!
    targ1 = q2sq*sinov(12)*sinov(12) + 2.0*q2(1)*t1ex*cosov(12) - q2(1)*q2(1) - t1ex*t1ex
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
      dw2valdr = (1.0/(1.0 + sarg21*sarg21))*((sarg1*dsarg2dr - dsarg1dr*sarg2)/(sarg1*sarg1))
      drrdr(1) = drdw(1)
      drrdr(2) = drdw(2)*cw2 + drdw(3)*sw2 + (-r(2)*sw2 + r(3)*cw2)*dw2valdr
      drrdr(3) = -drdw(2)*sw2 + drdw(3)*cw2 + (-r(2)*cw2 - r(3)*sw2)*dw2valdr
      dt1exdr = t1vp(1,1)*drrdr(1) + t1vp(1,2)*drrdr(2) + t1vp(1,3)*drrdr(3)
      tslp2 = 2.0*q2(1)*cosov(12)*dt1exdr - 2.0*t1ex*dt1exdr
    end if
!
    if (targ1.gt.0.0) then
      isreal2 = .true.
    else
      isreal2 = .false.
    end if
!
  end if
!
end subroutine checkreal3_dj
!
! -------------------------------------------------------------------  
! bisection(w1)
subroutine bisection2_dj(w1,next,t,nsolu)
!
  implicit none
!
  integer t,nsolu
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
!   ! find f(midpoint)
    call get_newccmat2_dj(w1val,newccmat2,grr,t,really)
    call get_newccmat2_dj(midpoint,newccmat2,grr2,t,really2)
    if ((really.EQV..false.).OR.(really2.EQV..false.)) then
      dj_wrncnt(3) = dj_wrncnt(3) + 1
      if (dj_wrncnt(3).eq.dj_wrnlmt(3)) then
        write(ilog,*) 'Bisection method failed for Dinner-Ulmschneider concerted rotation&
 & algorithm without omega bond sampling. This indicates significantly too coarse search&
 & settings (interval inconsistencies).'
        write(ilog,*) 'This was warning #',dj_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*dj_wrnlmt(3).gt.0.5*HUGE(dj_wrnlmt(3))) then
          dj_wrncnt(3) = 0 ! reset
        else
          dj_wrnlmt(3) = dj_wrnlmt(3)*10
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
end subroutine bisection2_dj
!
! -------------------------------------------------------------------  
!!
!! ! subroutine to find w2
!subroutine get_w2_dj(w1val,w2val,t)
!!
!  integer t
!  RTYPE w,arg1,arg2,w1val,w2val,denom,enum2,enum11,enum12,rnorm
!!
!  call r_of_w1w2_dj(r,w1val,rnorm)
!
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  denom = arg1
!  arg2 = w*w*q1(2)*q1(2)
!  enum2 = denom - arg2
!! 23  format(3(f22.16))
!!  write(*,23) denom,w,(sum(r(:)*r(:)) + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  enum2 = sqrt(enum2)
!  enum11 = r(2)*q1(2) + r(3)*q1(3)
!  enum12 = r(3)*q1(2) - r(2)*q1(3)
!!     
!  if ((t.eq.1).OR.(t.eq.2)) then
!    arg1 = (enum11*w*q1(2) + enum12*enum2)/denom
!    arg2 = (enum12*w*q1(2) - enum11*enum2)/denom
!  else if ((t.eq.3).OR.(t.eq.4)) then
!    arg1 = (enum11*w*q1(2) - enum12*enum2)/denom
!    arg2 = (enum12*w*q1(2) + enum11*enum2)/denom
!  else
!    write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
! &ain closure algorithm) must be between 1 and 4.'
!    call fexit()
!  end if
!!
!  w2val = atan2(arg2,arg1) ! or use ccmat
!  call genR(w2val,R2,invR2)
!!
!end subroutine get_w2_dj
!!
! -------------------------------------------------------------------  
!
subroutine get_newccmat2_dj(w1val,newccmat2,g,t,isreal2)
!
  integer, INTENT(IN):: t
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: newccmat2(4,6),g
  RTYPE arg1,arg2,w2val,t1ex,t1y,t1z,t2val(3),w4val,w3val,w5val,howfaroff,targ1,enum21,enum22
  RTYPE ay,pos8(3),TRTR(3),R1TR(3),T1R2x(3),T1R2(3,3),T2R3x(3),T2R3(3,3),dums(3),r(3)
  RTYPE T3R4x(3),T3R4(3,3),T4R5(3,3),T4R5x(3),T5R6(3,3),T5R6p6(3),ctw4,stw4,stw2,ctw2,t1vp(3,3),rr(3)
  logical afalse
  logical, INTENT(OUT):: isreal2
!
  afalse = .false.
!
! first get w2
  call checkreal3_dj(w1val,isreal2,howfaroff,targ1,nouse(1),nouse(2),w2val,stw2,ctw2,t1ex,r,t,afalse)
!
  if (isreal2.EQV..false.) then
    return ! deal on the outside
  end if
!
  t1vp = matmul(invT3,matmul(invR3,invT2))
  rr = matmul(invR2,r)
  t1ex = t1vp(1,1)*(rr(1)-q1(1)) + t1vp(1,2)*(rr(2)-q1(2)) + t1vp(1,3)*(rr(3)-q1(3))
  t1y = t1vp(2,1)*(rr(1)-q1(1)) + t1vp(2,2)*(rr(2)-q1(2)) + t1vp(2,3)*(rr(3)-q1(3))
  t1z = t1vp(3,1)*(rr(1)-q1(1)) + t1vp(3,2)*(rr(2)-q1(2)) + t1vp(3,3)*(rr(3)-q1(3))
  enum21 = q2sq*sinov(12)*sinov(12) + 2.0*q2(1)*t1ex*cosov(12)
  enum22 = q2(1)*q2(1) + t1ex*t1ex
!
  dums(3) = sqrt(enum21-enum22)
  if ((t.eq.1).OR.(t.eq.3)) then
    arg1 = (q2(2)*(q2(1)*cosov(12) - t1ex) + q2(3)*dums(3))/w4denom
    arg2 = (q2(3)*(-q2(1)*cosov(12) + t1ex) + q2(2)*dums(3))/w4denom
  else if ((t.eq.2).OR.(t.eq.4)) then
    arg1 = (q2(2)*(q2(1)*cosov(12) - t1ex) - q2(3)*dums(3))/w4denom
    arg2 = (q2(3)*(-q2(1)*cosov(12) + t1ex) - q2(2)*dums(3))/w4denom
  else
    write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
 &ain closure algorithm) must be between 1 and 4.'
    call fexit()
  end if
!
  w4val = atan2(arg2,arg1)
  ctw4 = cos(w4val)
  stw4 = sin(w4val)
!
! -------------------------------------------------------------------  
!
  t2val(1) = q2(1)*cosov(12) - q2(2)*sinov(12)*ctw4 + q2(3)*sinov(12)*stw4
  t2val(2) = q2(1)*sinov(12) + q2(2)*cosov(12)*ctw4 - q2(3)*cosov(12)*stw4
  t2val(3) = q2(2)*stw4 + q2(3)*ctw4
!
  arg1 = (t1y*t2val(2) + t1z*t2val(3))/(t1y*t1y + t1z*t1z)
  arg2 = (t1z*t2val(2) - t1y*t2val(3))/(t1y*t1y + t1z*t1z)
  w3val = atan2(arg2,arg1)
!
! -------------------------------------------------------------------  
!
  call genR(w3val,R4,invR4)
  call genR(w4val,R5,invR5)
  call genR(w1val,R1,invR1)
  invmat = matmul(invT6,matmul(invR6,matmul(invT5,matmul(invR5,matmul(invT4,matmul&
 &(invR4,matmul(invT3,matmul(invR3,matmul(invT2,matmul(invR2,matmul(invT1,invR1)))))))))))
!
  arg1 = sum(invmat(2,:)*u(:))/sinov(15)
  arg2 = sum(invmat(3,:)*u(:))/sinov(15)
!
  w5val = atan2(arg2,arg1)
!  write(*,*) "w7",w5val
  newccmat2(t,1) = w1val
  newccmat2(t,2) = w2val
  newccmat2(t,3) = w3val
  newccmat2(t,4) = w4val
  newccmat2(t,5) = w5val
!
  T1R2 = matmul(T1,R2)
  T2R3 = matmul(T2,R3)
  T3R4 = matmul(T3,R4)
  T4R5 = matmul(T4,R5)
  T5R6 = matmul(T5,R6)
  T5R6p6(1) = T5R6(1,1)*p6 
  T5R6p6(2) = T5R6(2,1)*p6
  T5R6p6(3) = T5R6(3,1)*p6
  T5R6p6(1) = T5R6p6(1) + p5
  do k=1,3
    T4R5x(k) = sum(T4R5(k,:)*T5R6p6(:))
  end do
  T4R5x(1) = T4R5x(1) + p4
  do k=1,3
    T3R4x(k) = sum(T3R4(k,:)*T4R5x(:))
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
  pos8(:) = TRTR(:) + posNI(:)
  call dihed(pos8,posCAF,posCIF,posNF,ay)
  newccmat2(t,6) = ay
!
  g = u(1)*invmat(1,1) + u(2)*invmat(1,2) + u(3)*invmat(1,3) - cosov(15)
!  write(*,*) "w8",ay
!
end subroutine get_newccmat2_dj
!
! -------------------------------------------------------------------  
!
! subroutine outputs g(w1)
!!
!subroutine g_of_w1w2_dj(w1val,g,t)
!!
!  integer t
!  RTYPE w1val,g
!  RTYPE TRmat(3,3)
!!
!! using branch t get solution for g(w1)
!  call get_newccmat2_dj(w1val,newccmat2,t)
!!
!! calculate g(w1)
!  TRmat = matmul(invT6,matmul(invR6,matmul(invT5,matmul(invR5,matmul(invT4,matmul&
! &(invR4,matmul(invT3,matmul(invR3,matmul(invT2,matmul(invR2,matmul(invT1,invR1)))))))))))
!  g = u(1)*TRmat(1,1) + u(2)*TRmat(1,2) + u(3)*TRmat(1,3) - cosov(15)
!!
!end subroutine g_of_w1w2_dj
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
end subroutine torcrchainclose_dj
!
!
!
! -----------------------------------------------------------------
subroutine refvals(i1,i2,i3,i4,i5,i6,i7,i8,i9,ccmat,oldvals2)
!------------------------------------------------------------
! subroutine to store values needed for chain closure
! INPUTS:
!     rsf -  reference residue for chain closure
! MODIFIES:
!     oldvals2 - {p1,p2,p3,p4,p5,p6,a1,a2,a3,a4,a5,a6}
!     ccmat - {w1,w2,w3,w4,w5,w6}
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
  integer i1,i2,i3,i4,i5,i6,i7,i8,i9
  RTYPE ccmat(6),oldvals2(12),getblen,getztor,getbang
!
! get old bond lengths
  oldvals2(1) = getblen(i2,i3)
  oldvals2(2) = getblen(i3,i4)
  oldvals2(3) = getblen(i4,i5)
  oldvals2(4) = getblen(i5,i6)
  oldvals2(5) = getblen(i6,i7)
  oldvals2(6) = getblen(i7,i8)
!
! get old bond angles
  oldvals2(7) = PI - getbang(i2,i3,i4)/RADIAN
  oldvals2(8) = PI - getbang(i3,i4,i5)/RADIAN
  oldvals2(9) = PI - getbang(i4,i5,i6)/RADIAN
  oldvals2(10) = PI - getbang(i5,i6,i7)/RADIAN
  oldvals2(11) = PI - getbang(i6,i7,i8)/RADIAN
  oldvals2(12) = PI - getbang(i7,i8,i9)/RADIAN
!
! get old torsions
  ccmat(1) = getztor(i1,i2,i3,i4)/RADIAN
  ccmat(2) = getztor(i2,i3,i4,i5)/RADIAN
  ccmat(3) = getztor(i3,i4,i5,i6)/RADIAN
  ccmat(4) = getztor(i4,i5,i6,i7)/RADIAN
  ccmat(5) = getztor(i5,i6,i7,i8)/RADIAN
  ccmat(6) = getztor(i6,i7,i8,i9)/RADIAN
!
end
!
!---------------------------------------------------------------------
subroutine torchainclose(i1,i2,i3,i7,i8,i9,ccmat,incrmat,oldvals2,solucount)
!------------------------------------------------------------
! Subroutine to calculate the constraints for chain closure
! for UJ prerotation.
! INPUTS:
!     rsi,rsf -  reference residues for chain closure
!     oldvals2 - reference bond lengths and angles
!     ccmat - reference torsion angle values
! OUTPUT:
!     solmat - array of change in torsions needed for closure (1:solucount,1:6)
!     solucount - number of solutions
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
  use zmatrix
!
  implicit none
!
  integer j,k,kk,solucount,aone,athree,afive
  integer i1,i2,i3,i7,i8,i9
  logical logi1,logi1a,logi1b,goodsol,atrue,really,really2,really3!slog(3)
  integer findsolu(4),kkk,jj,nintvs(4),jjj,npreintvs,atwo!,dums(20)
  RTYPE w,w1,ccmat(6),dpfmat(MAXUJDOF+6),oldvals2(12)
  RTYPE p1,p2,p3,p4,p5,p6,nopl,u(3),con1(3),posCI_0(3)
  RTYPE q0(3),posCA(3),posNI(3),posNF(3),posCIF(3),posCA_0(3),w1mat(4,MAXSOLU)
  RTYPE q1(3),q2(3),R3(3,3),R2(3,3),R1(3,3),T2(3,3),invT3(3,3)
  RTYPE invR1(3,3),invR2(3,3),invR3(3,3),invT1(3,3),invT2(3,3)
  RTYPE R4(3,3),invR4(3,3),T4(3,3),invT4(3,3)
  RTYPE T1(3,3),gn,g1,intercpt,newccmat2(4,6)
  RTYPE a3,ax,s(3),T3(3,3),epsilon,midpoint,ltestw1,htestw1
  RTYPE con1xy(3),a1,ay,Ra1(3,3),Ra2(3,3),Ra3(3,3),az,increment,gnpr,gnpo,preintvbs(MAXSOLU*4,2)
  RTYPE newxaxis(3),newyaxis(3),newzaxis(3),yaxis(3),zaxis(3),a2,oldgs(4),oldws(4),maxincr
  RTYPE xaxis(3),RaE(3,3),or_pl(3),or_pl2(3),lastslope(4),curslope(4)
  RTYPE con0(3),origin(3),s1(3),searchstp,con6(3),blas(3)
  RTYPE invmat(3,3),u1,invRaE(3,3),incrmat(MAXSOLU*4,6),intvalbnds(4,MAXSOLU*4,2)
  RTYPE sinov(12),cosov(12),q1sq,q2sq,svq(3),w4denom,nouse(17),newt
!
! test for closure and prerotation
  blas(1) = x(i8)
  blas(2) = y(i8)
  blas(3) = z(i8)
!
! helpers
  aone = 1
  atwo = 2
  athree = 3
  afive = 5
  atrue = .true.
!
! get the current position of atoms used in chain closure
  posCA_0(1) = x(i1)
  posCA_0(2) = y(i1)
  posCA_0(3) = z(i1)
  posCI_0(1) = x(i2)
  posCI_0(2) = y(i2)
  posCI_0(3) = z(i2)
  posNI(1) = x(i3)
  posNI(2) = y(i3)
  posNI(3) = z(i3)
  posCA(1) = x(i7)
  posCA(2) = y(i7)
  posCA(3) = z(i7)
  posCIF(1) = x(i8)
  posCIF(2) = y(i8)
  posCIF(3) = z(i8) 
  posNF(1) = x(i9)
  posNF(2) = y(i9)
  posNF(3) = z(i9)
!
  con0(:) = posCA_0(:)-posCI_0(:)
  con1(:) = posNI(:)-posCI_0(:)
  con6(:) = posCIF(:)-posCA(:)
!
! get bond lengths
  p1 = oldvals2(1)
  p2 = oldvals2(2)
  p3 = oldvals2(3)
  p4 = oldvals2(4)
  p5 = oldvals2(5)
  p6 = oldvals2(6)
  call genT(oldvals2(7),T1,invT1)
  call genT(oldvals2(9),T3,invT3)
  sinov(7) = sin(oldvals2(7))
  cosov(7) = cos(oldvals2(7))
  sinov(8) = sin(oldvals2(8))
  cosov(8) = cos(oldvals2(8))
  sinov(9) = sin(oldvals2(9))
  cosov(9) = cos(oldvals2(9))
  sinov(11) = sin(oldvals2(11))
  cosov(11) = cos(oldvals2(11))

! build Euler rotation matrix
  origin(1) = 0.0
  origin(2) = 0.0
  origin(3) = 0.0
  xaxis = (/ 1, 0, 0 /)
  yaxis = (/ 0, 1, 0 /)
  zaxis = (/ 0, 0, 1 /)

! projection of con1 on xy-plane
  con1xy(1) = posNI(1) - posCI_0(1)
  con1xy(2) = posNI(2) - posCI_0(2)
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
    
! set q0 = p1
  q0(1) = p1
  q0(2) = 0.0
  q0(3) = 0.0

! get q1 an q2
  call genT(oldvals2(8),T2,invT2)
  q1(1) = T2(1,1)*p3 + p2
  q1(2) = T2(2,1)*p3
  q1(3) = T2(3,1)*p3
  q1sq = sum(q1(:)*q1(:))
  call genT(oldvals2(10),T4,invT4)
  q2(1) = T4(1,1)*p5 + p4
  q2(2) = T4(2,1)*p5
  q2(3) = T4(3,1)*p5
  q2sq = sum(q2(:)*q2(:))
  w4denom = (q2(2)*q2(2)+q2(3)*q2(3))*sinov(9)

! get pos S = RaE.posNF
  s1(1) = posCA(1)-posCI_0(1)
  s1(2) = posCA(2)-posCI_0(2)
  s1(3) = posCA(3)-posCI_0(3)
  s(1) = sum(RaE(1,:)*s1(:))
  s(2) = sum(RaE(2,:)*s1(:))
  s(3) = sum(RaE(3,:)*s1(:))
  svq(:) = s(:) - q0(:)

! find u
  u(1) = sum(RaE(1,:)*con6(:))
  u(2) = sum(RaE(2,:)*con6(:))
  u(3) = sum(RaE(3,:)*con6(:))
  u1 = sqrt(sum(u(:)*u(:)))
  u(:) = u(:)/u1

! propagate search variable w1 via Newton to identify branch-independent R-solution space 
  w1mat(:,:) = 0.0
  findsolu(:) = 0
!
  kk = 0
  w1 = -PI
  g1 = 0.0
  increment = 0.0
  npreintvs = 0
  logi1b = .false.
  do while (1.eq.1)
    oldws(1) = w1
    oldgs(1) = g1
    w1 = w1 + increment
    kk = kk + 1
    if (w1.gt.PI) exit
!   curslope(1) is analytical derivative dg1/dw1: g1 > 0 means real solution
    call checkreal1(w1,logi1,g1,curslope(1))
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
        increment = min(UJ_params(4)/RADIAN,0.8*newt)
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
        increment = min(UJ_params(4)/RADIAN,0.8*newt)
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
        call checkreal3(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
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
          increment = min(UJ_params(4)/RADIAN,0.5*newt,maxincr)
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
          increment = min(UJ_params(4)/RADIAN,0.8*newt,maxincr)
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
!    write(*,*) kkk,nintvs(kkk),kk
    if ((logi1b.EQV..false.).AND.(logi1.EQV..true.)) then
!     this shouldn't really happen
      call checkreal1(oldws(1),logi1,g1,curslope(1))
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
            call checkreal3(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2(w1,newccmat2,g1,jj,really)
!             remember that checkreal3 is called by get_newccmat2 again -> really-check redundant
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
        intvalbnds(kkk,j,1) = intvalbnds(kkk,j,1) + 1.0d-9
        intvalbnds(kkk,j,2) = intvalbnds(kkk,j,2) - 1.0d-9
!       check one more time
        if (abs(intvalbnds(kkk,j,1)-intvalbnds(kkk,j,2)).lt.2.0d-9) then
!         write(*,*) 'vanishing interval:',j,intvalbnds(kkk,j,1)
          do jj=2*kkk-1,2*kkk
            w1 = 0.5*(intvalbnds(kkk,j,1)+intvalbnds(kkk,j,2))
            call checkreal3(w1,logi1,g1,gn,curslope(1),curslope(2),nouse(1),nouse(2),nouse(3),nouse(4),&
 &nouse(5:7),jj,atrue)
            if (logi1.EQV..true.) then
              call get_newccmat2(w1,newccmat2,g1,jj,really)
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
!    write(*,*) 'DO',kkk,kk,' steps'
!    do jj=1,nintvs(kkk)
!    write(*,*) kkk,jj,RADIAN*intvalbnds(kkk,jj,1:2)
!    end do
  end do
!
! and lastly scan them
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
        oldws(kkk) = w1
        if (w1.ge.intvalbnds(jj,j,2)) exit
        if ((w1.lt.intvalbnds(jj,j,2)).AND.(w1+increment.gt.intvalbnds(jj,j,2))) then
          w1 = intvalbnds(jj,j,2)
        else
          w1 = w1 + increment
        end if
        oldgs(1) = g1
        call get_newccmat2(w1,newccmat2,g1,kkk,really)
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
          call get_newccmat2(htestw1,newccmat2,gnpo,kkk,really2)
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
          call get_newccmat2(ltestw1,newccmat2,gnpr,kkk,really3)
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
          call get_newccmat2(ltestw1,newccmat2,gnpr,kkk,really3)
!         bail out if we nail it
          if ((really3.EQV..true.).AND.(gnpr.eq.0.0)) then
!          write(*,*) 'nailed it pr',w1*RADIAN,' in ',kkk
            findsolu(kkk) = findsolu(kkk) + 1
            w1mat(kkk,findsolu(kkk)) = ltestw1
            increment = 1.0d-9
            cycle
          end if
          call get_newccmat2(htestw1,newccmat2,gnpo,kkk,really2)
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
          do_wrncnt(2) = do_wrncnt(2) + 1
          if (do_wrncnt(2).eq.do_wrnlmt(2)) then
            write(ilog,*) 'Warning. Encountered interval inconsistency in loop closure of generic CR algorithm&
 & for six subsequent dihedrals. This means a number was tested which was previously found to give a real solu&
 &tion and now produced NaN. This is ultimately a bug the causes of which are currently poorly understood.'
            write(ilog,*) 'This was warning #',do_wrncnt(2),' of this type not all of which may be displayed.'
            if (10.0*do_wrnlmt(2).gt.0.5*HUGE(do_wrnlmt(2))) then
              do_wrncnt(2) = 0 ! reset
            else
              do_wrnlmt(2) = do_wrnlmt(2)*10
            end if
          end if
          exit
        end if
!       we can bisect the search interval straight away
        if (((gnpr.lt.0.0).AND.(gnpo.gt.0.0)).OR.((gnpr.gt.0.0).AND.(gnpo.lt.0.0))) then
!          write(*,*) 'bisecting ',ltestw1*RADIAN,' and ',htestw1*RADIAN,' in ',kkk
          call bisection2(ltestw1,htestw1,kkk,findsolu(kkk))
          increment = 2.0d-9
          g1 = 0.0 ! sets oldgs to zero
          cycle
        end if
!       check for accidental crossover
        if (((oldgs(kkk).lt.0.0).AND.(g1.gt.0.0)).OR.((oldgs(kkk).gt.0.0).AND.(g1.lt.0.0))) then ! yes
!          write(*,*) 'crossover: ',oldws(kkk)*RADIAN,' and ',w1*RADIAN,' in ',kkk
          call bisection2(oldws(kkk),w1,kkk,findsolu(kkk))
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
! note that we filter solutions for two things: exactness of closure and NaN
! this incapacitates any outside checks of this type
!
  if (solucount.eq.0) then    
!    just exit if we found no solution
!    write(ilog,*) 'Did not find a solution to the chain closure'
    return
  else if (solucount.eq.1) then
!   for a single solution it's easy: assign, restore, perturb, get jacobian, exit
    do j=1,4
      if(findsolu(j).gt.0) then
        call get_newccmat2(w1mat(j,1),newccmat2,nouse(1),j,really)
        goodsol = .true.
        do jjj=1,6
          if ((newccmat2(j,jjj).lt.-PI).OR.(newccmat2(j,jjj).gt.PI).OR.(really.EQV..false.)) then
            goodsol = .false.
          end if
        end do
        if (goodsol.EQV..false.) then
          solucount = solucount - 1
          return
        end if
        incrmat(1,1:6) = newccmat2(j,1:6) - ccmat(1:6)
        exit
      end if
    end do

    dpfmat(1:6) = incrmat(1,1:6)
!
  else
!
!
    kk = 0
    logi1 = .false.
    do j=1,4
!
      do kkk=1,findsolu(j)
!
        kk = kk + 1
        call get_newccmat2(w1mat(j,kkk),newccmat2,nouse(1),j,really)
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
        incrmat(kk,1:6) = newccmat2(j,1:6) - ccmat(1:6)
        dpfmat(1:6) = incrmat(kk,1:6)
!
      end do
!  
    end do
!
  end if
!
  CONTAINS
!
! -------------------------------------------------------------------  
!
! r(w1) = invT1*invR1*(s-q0)
subroutine r_of_w1w2(r,w1val,rnorm)
!
  RTYPE, INTENT(OUT):: r(3),rnorm
  RTYPE, INTENT(IN):: w1val
!
  call genR(w1val,R1,invR1)
  invmat = matmul(invT1,invR1)
!
  r(1) = sum(invmat(1,:)*svq(:))
  r(2) = sum(invmat(2,:)*svq(:))
  r(3) = sum(invmat(3,:)*svq(:))
!
  rnorm = sum(r(:)*r(:))
!
end subroutine r_of_w1w2
!
!-------------------------------------------------------------------------
!
! get arg1, arg2, and w and (on request) some derivatives
! 
subroutine deriva_w1(sw1,cw1,rnorm,r,arg1,arg2,darg1dw1,darg2dw1,w,drdw1,dwndw1,do_deriv)
!
  RTYPE, INTENT(IN):: r(3),sw1,cw1,rnorm
  logical, INTENT(IN):: do_deriv
  RTYPE, INTENT(OUT):: arg1,arg2,darg1dw1,darg2dw1,w,drdw1(3),dwndw1
!
  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
  arg2 = w*w*q1(2)*q1(2)
!
  if (do_deriv.EQV..true.) then
    drdw1(1) = -svq(2)*sw1*sinov(7) + svq(3)*cw1*sinov(7)
    drdw1(2) = -svq(2)*sw1*cosov(7) + svq(3)*cw1*cosov(7)
    drdw1(3) = -svq(2)*cw1 - svq(3)*sw1
    dwndw1 = (1.0/(2.0*q1(2)))*(-2.0*drdw1(1)*q1(1) + sum(2.0*r(:)*drdw1(:)))
    darg1dw1 = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw1(2) + r(3)*drdw1(3)))
    darg2dw1 = 2.0*q1(2)*q1(2)*w*dwndw1
  end if
!
end subroutine deriva_w1
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal1(w1val,isreal,howfaroff,testslope)
!
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: howfaroff,testslope
  RTYPE arg1,arg2,darg1dr,darg2dr,drdw(3),rnorm,sw1,cw1,r(3),dwndr
  logical, INTENT(OUT):: isreal
  logical atrue
!
  isreal = .false.
  atrue = .true.
!
  call r_of_w1w2(r,w1val,rnorm)
  sw1 = sin(w1val)
  cw1 = cos(w1val)
  call deriva_w1(sw1,cw1,rnorm,r,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,atrue)
!
  howfaroff = arg1 - arg2
  testslope = darg1dr - darg2dr
  if (howfaroff.gt.0.0) then
    isreal = .true.
  else
    isreal = .false.
  end if
!
end subroutine checkreal1
!
! -------------------------------------------------------------------  
!
! check for real solutions to w2
subroutine checkreal3(w1val,isreal2,howfaroff,targ1,tslp1,tslp2,w2val,sw2,cw2,t1ex,r,t,do_deriv)
!
  logical, INTENT(IN):: do_deriv
  logical, INTENT(OUT):: isreal2
  RTYPE, INTENT(OUT):: howfaroff,targ1,tslp1,tslp2,cw2,sw2
  RTYPE w1val,arg1,arg2,darg1dr,darg2dr,dwndr,drdw(3),r(3)
  RTYPE denom,enum11,enum21,enum2,sarg1,sarg2,w2val,t1ex,denum2dr,rnorm,sw1,cw1
  RTYPE dsarg1dr,dsarg2dr,denum21dr,denum11dr,dw2valdr,dt1exdr,sarg21
  logical isreal
  integer, INTENT(IN):: t
!
  tslp2 = 0.0
  targ1 = 0.0
  isreal = .false.
!
  call r_of_w1w2(r,w1val,rnorm)
  sw1 = sin(w1val)
  cw1 = cos(w1val)
  call deriva_w1(sw1,cw1,rnorm,r,arg1,arg2,darg1dr,darg2dr,w,drdw,dwndr,do_deriv)
!  drdw(1) = -svq(2)*sw1*sinov(7) + svq(3)*cw1*sinov(7)
!  drdw(2) = -svq(2)*sw1*cosov(7) + svq(3)*cw1*cosov(7)
!  drdw(3) = -svq(2)*cw1 - svq(3)*sw1
!  w = (rnorm + q1sq - q2sq - 2.0*r(1)*q1(1))/(2.0*q1(2))
!  dwnewdr = (1.0/(2.0*q1(2)))*(-2.0*drdw(1)*q1(1) + sum(2.0*r(:)*drdw(:)))
!  arg1 = (r(2)*r(2) + r(3)*r(3))*(q1(2)*q1(2) + q1(3)*q1(3))
!  darg1dr = (q1(2)*q1(2) + q1(3)*q1(3))*(2.0*(r(2)*drdw(2) + r(3)*drdw(3)))
!  arg2 = w*w*q1(2)*q1(2)
!  darg2dr = 2.0*q1(2)*q1(2)*w*dwnewdr
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

    if ((t.eq.1).OR.(t.eq.2)) then
      enum21 = (r(3)*q1(2) - r(2)*q1(3))*w*q1(2) - (r(2)*q1(2) + r(3)*q1(3))*enum2
      enum11 = (r(2)*q1(2) + r(3)*q1(3))*w*q1(2) + (r(3)*q1(2) - r(2)*q1(3))*enum2
    else if  ((t.eq.3).OR.(t.eq.4)) then
      enum21 = (r(3)*q1(2) - r(2)*q1(3))*w*q1(2) + (r(2)*q1(2) + r(3)*q1(3))*enum2
      enum11 = (r(2)*q1(2) + r(3)*q1(3))*w*q1(2) - (r(3)*q1(2) - r(2)*q1(3))*enum2
    else
      write(ilog,*) 'Fatal error. Assignment of t (branch number in ch&
 &ain closure algorithm) must be between 1 and 4.'
      call fexit()
    end if
    sarg1 = enum11/denom
    sarg2 = enum21/denom
    w2val = atan2(sarg2,sarg1) ! or use ccmat
    sw2 = sin(w2val)
    cw2 = cos(w2val)
    sarg21 = (sarg2/sarg1)
    t1ex = T2(1,1)*(r(1)-q1(1)) + T2(2,1)*(r(2)*cw2+r(3)*sw2-q1(2))
    targ1 = q2sq*sinov(9)*sinov(9) + 2.0*q2(1)*t1ex*cosov(9) - q2(1)*q2(1) - t1ex*t1ex
!
    if (do_deriv.EQV..true.) then
      denum2dr = tslp1/(2.0*enum2)
      if ((t.eq.1).OR.(t.eq.2)) then
        denum21dr = (-q1(2)*drdw(2) - q1(3)*drdw(3))*enum2 - (r(2)*q1(2) + r(3)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(3)*q1(2) - r(2)*q1(3)) + w*(drdw(3)*q1(2) - drdw(2)*q1(3)))
        denum11dr = (q1(2)*drdw(3) - q1(3)*drdw(2))*enum2 + (r(3)*q1(2) - r(2)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(2)*q1(2) + r(3)*q1(3)) + w*(drdw(2)*q1(2) + drdw(3)*q1(3)))
      else if  ((t.eq.3).OR.(t.eq.4)) then
        denum21dr = (q1(2)*drdw(2) + q1(3)*drdw(3))*enum2 + (r(2)*q1(2) + r(3)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(3)*q1(2) - r(2)*q1(3)) + w*(drdw(3)*q1(2) - drdw(2)*q1(3)))
        denum11dr = (-q1(2)*drdw(3) + q1(3)*drdw(2))*enum2 - (r(3)*q1(2) - r(2)*q1(3))*denum2dr + &
   & q1(2)*(dwndr*(r(2)*q1(2) + r(3)*q1(3)) + w*(drdw(2)*q1(2) + drdw(3)*q1(3)))
      end if
      dsarg1dr = (denom*denum11dr - enum11*enum2*darg1dr)/(denom*denom)
      dsarg2dr = (denom*denum21dr - enum21*enum2*darg1dr)/(denom*denom)
      dw2valdr = (1.0/(1.0 + sarg21*sarg21))*((sarg1*dsarg2dr - dsarg1dr*sarg2)/(sarg1*sarg1))
      dt1exdr = T2(1,1)*drdw(1) + T2(2,1)*(-r(2)*sw2*dw2valdr + &
 &           drdw(2)*cw2 + drdw(3)*sw2 + r(3)*cw2*dw2valdr)
      tslp2 = 2.0*q2(1)*cosov(9)*dt1exdr - 2.0*t1ex*dt1exdr
!
    end if
!
    if (targ1.gt.0.0) then
      isreal2 = .true.
    else
      isreal2 = .false.
    end if
  end if
!  write(*,*) RADIAN*w2val!RADIAN*w1val,t,isreal,isreal2
!
end subroutine checkreal3
!
! -------------------------------------------------------------------  
!
! bisection(w1): note that the w1-values are unwrapped on the outside
!
subroutine bisection2(w1,next,t,nsolu)
!
  implicit none
!
  integer t,nsolu
  logical really,really2
  RTYPE w1,next,w1val,nextw1,grr,grr2
!
  w1val = w1
  nextw1 = next
!
  epsilon = 2.0*(10.0**(-15.0))
! ! start loop
  do while (abs(nextw1 - w1val).gt.(2.0*epsilon))
!
    midpoint = (nextw1 + w1val) / 2.0
!
!   ! find f(midpoint)
    call get_newccmat2(w1val,newccmat2,grr,t,really)
    call get_newccmat2(midpoint,newccmat2,grr2,t,really2)
    if ((really.EQV..false.).OR.(really2.EQV..false.)) then
      do_wrncnt(3) = do_wrncnt(3) + 1
      if (do_wrncnt(3).eq.do_wrnlmt(3)) then
        write(ilog,*) 'Bisection method failed for generic CR algorithm for six subsequent &
 &dihedral angles. This indicates significantly too coarse search settings (interval inconsistencies).'
        write(ilog,*) 'This was warning #',do_wrncnt(3),' of this type not all of which may be displayed.'
        if (10.0*do_wrnlmt(3).gt.0.5*HUGE(do_wrnlmt(3))) then
          do_wrncnt(3) = 0 ! reset
        else
          do_wrnlmt(3) = do_wrnlmt(3)*10
        end if
      end if
      return
    end if
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
!  write(*,*) t,UJ_params(4),w1val,newccmat2(t,2)
  
!
end subroutine bisection2
! -------------------------------------------------------------------  
!
! subroutine to find w4
subroutine get_newccmat2(w1val,newccmat2,g,t,isreal2)
!
  integer, INTENT(IN):: t
  RTYPE, INTENT(IN):: w1val
  RTYPE, INTENT(OUT):: newccmat2(4,6),g
  RTYPE arg1,arg2,w2val,t1ex,t1y,t1z,t2val(3),w4val,w3val,w5val,howfaroff,targ1
  RTYPE ay,pos6(3),TRTR(3),R1TR(3),T1R2x(3),T1R2(3,3),T2R3x(3),T2R3(3,3),dums(3),r(3)
  RTYPE T3R4p4(3),T3R4(3,3),ctw4,stw4,stw2,ctw2
  logical afalse,isreal2
!
  afalse = .false.
!
! first get w2
  call checkreal3(w1val,isreal2,howfaroff,targ1,nouse(1),nouse(2),w2val,stw2,ctw2,t1ex,r,t,afalse)
  call genR(w2val,R2,invR2)
!
  if (isreal2.EQV..false.) then
    return ! deal on the outside
  end if
!  call get_w2(w1val,w2val,t)
!  stw2 = sin(w2val)
!  ctw2 = cos(w2val)
!  t1x = T2(1,1)*(r(1)-q1(1)) + T2(2,1)*(r(2)*cos(w2val)+r(3)*sin(w2val)-q1(2))
!  denom = (q2(2)*q2(2)+q2(3)*q2(3))*sinov(9)
!  enum21 = q2sq*sinov(9)*sinov(9) + 2.0*q2(1)*t1x*cosov(9)
!  enum22 =  q2(1)*q2(1) + t1x*t1x
!  dums(3) = sqrt(enum21-enum22)
  dums(3) = sqrt(targ1)
!!
  if ((t.eq.1).OR.(t.eq.3)) then
    arg1 = (q2(2)*(q2(1)*cosov(9) - t1ex) + q2(3)*dums(3))/w4denom
    arg2 = (q2(3)*(-q2(1)*cosov(9) + t1ex) - q2(2)*dums(3))/w4denom
  else if ((t.eq.2).OR.(t.eq.4)) then
    arg1 = (q2(2)*(q2(1)*cosov(9) - t1ex) - q2(3)*dums(3))/w4denom
    arg2 = (q2(3)*(-q2(1)*cosov(9) + t1ex) + q2(2)*dums(3))/w4denom
  end if
!
  w4val = atan2(arg2,arg1)
  ctw4 = cos(w4val)
  stw4 = sin(w4val)
!
! -------------------------------------------------------------------  
!
  t2val(1) = q2(1)*cosov(9) - q2(2)*sinov(9)*ctw4 + q2(3)*sinov(9)*stw4
  t2val(2) = q2(1)*sinov(9) + q2(2)*cosov(9)*ctw4 - q2(3)*cosov(9)*stw4
  t2val(3) = q2(2)*stw4 + q2(3)*ctw4
! 
  t1y = -sinov(8)*(r(1)-q1(1)) + cosov(8)*(r(2)*ctw2 + r(3)*stw2 - q1(2))
  t1z = -r(2)*stw2 + r(3)*ctw2 - q1(3)
!
  arg1 = (t1y*t2val(2) + t1z*t2val(3))/(t1y*t1y + t1z*t1z)
  arg2 = (t1z*t2val(2) - t1y*t2val(3))/(t1y*t1y + t1z*t1z)
  w3val = atan2(arg2,arg1)
!
! -------------------------------------------------------------------  
!
  call genR(w3val,R3,invR3)
  call genR(w4val,R4,invR4)
  call genR(w1val,R1,invR1)
  invmat = matmul(invT4,matmul(invR4,matmul(invT3,matmul(invR3,matmul(invT2,matmul(invR2,&
 &          matmul(invT1,invR1)))))))
!
  arg1 = sum(invmat(2,:)*u(:))/sinov(11)
  arg2 = sum(invmat(3,:)*u(:))/sinov(11)
!
  w5val = atan2(arg2,arg1)
  newccmat2(t,1) = w1val
  newccmat2(t,2) = w2val
  newccmat2(t,3) = w3val
  newccmat2(t,4) = w4val
  newccmat2(t,5) = w5val
!
  T1R2 = matmul(T1,R2)
  T2R3 = matmul(T2,R3)
  T3R4 = matmul(T3,R4)
  T3R4p4(1) = T3R4(1,1)*p4 
  T3R4p4(2) = T3R4(2,1)*p4
  T3R4p4(3) = T3R4(3,1)*p4
  T3R4p4(1) = T3R4p4(1) + p3
  do k=1,3
    T2R3x(k) = sum(T2R3(k,:)*T3R4p4(:))
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
  pos6(:) = TRTR(:) + posCI_0(:)
  call dihed(pos6,posCA,posCIF,posNF,ay)
  newccmat2(t,6) = ay
!
  g = u(1)*invmat(1,1) + u(2)*invmat(1,2) + u(3)*invmat(1,3) - cosov(11)
!
end subroutine get_newccmat2
! ------------------------------------------------------------------
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
end subroutine torchainclose
!
!-------------------------------------------------------------------------------------------
!
! find nine consecutive atoms defining 6 torsional degrees of freedom which we can use to close the chain
!
subroutine crosslink_indices(i,i1,i2,i3,i4,i5,i6,i7,i8,i9)
!
  use iounit
  use sequen
  use molecule
  use polypep
  use system
!
  implicit none
!
  integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9
  integer rs1,rs2,shf
!
  rs1 = crosslink(i)%rsnrs(1)
  rs2 = crosslink(i)%rsnrs(2)
  i1 = -1
  i2 = -1
  i3 = -1
  i4 = -1
  i5 = -1
  i6 = -1
  i7 = -1
  i8 = -1
  i9 = -1
!
  if (molofrs(rs1).ne.molofrs(rs2)) then
    write(ilog,*) 'Fatal. Encountered intermolecular crosslink in crosslink_indices(...).&
 & Multiple intermolecular links for the same pair of molecules are not yet supported.'
    call fexit()
  end if
!
  if (crosslink(i)%itstype.eq.1) then
    if ((rs2.lt.rsmol(molofrs(rs2),1)+2).OR.(rs2.lt.(rs1+2))) then
      write(ilog,*) 'Fatal. The random generation of structures satisfying a disulfide bridge between neighboring residues &
 &is currently not supported.'
      call fexit()
    end if
    if (ci(rs2-1).gt.0) i1 = ci(rs2-1)
    if (ni(rs2).gt.0) i2 = ni(rs2)
    if (cai(rs2).gt.0) i3 = cai(rs2)
    shf = 0
    if (ua_model.gt.0) shf = 1
    if (at(rs2)%sc(2-shf).gt.0) i4 = at(rs2)%sc(2-shf) ! CB on rs2
    if (at(rs2)%sc(3-shf).gt.0) i5 = at(rs2)%sc(3-shf) ! SG on rs2
    if (at(rs1)%sc(3-shf).gt.0) i6 = at(rs1)%sc(3-shf) ! SG on rs1
    if (at(rs1)%sc(2-shf).gt.0) i7 = at(rs1)%sc(2-shf) ! CB on rs1
    if (cai(rs1).gt.0) i8 = cai(rs1)
    if (ni(rs1).gt.0) i9 = ni(rs1)
  else
    write(ilog,*) 'Fatal. Encountered unsupported crosslink type in crosslink_indices(...).&
 & This is most likely an omission bug.'
    call fexit()
  end if
!
  if ((i1.eq.-1).OR.(i2.eq.-1).OR.(i3.eq.-1).OR.(i4.eq.-1).OR.(i5.eq.-1).OR.(i6.eq.-1)&
 &.OR.(i7.eq.-1).OR.(i8.eq.-1).OR.(i9.eq.-1)) then
    write(ilog,*) 'Fatal. Missing reference atoms in crosslink_indices(...). This is most certainly&
 & a bug.'
    call fexit()
  end if
!
end
!
!------------------------------------------------------------------------------------------------------
!
! find nine consecutive atoms defining 6 torsional degrees of freedom which we can use to close the chain
!
subroutine crosslink_refvals(i,i1,i2,i3,i4,i5,i6,i7,i8,i9,dihmat,internals)
!
  use iounit
  use sequen
  use math
!
  implicit none
!
  integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9
  integer rs1,rs2
  RTYPE dihmat(6),internals(12),getblen,getztor,getbang
!
  rs1 = crosslink(i)%rsnrs(1)
  rs2 = crosslink(i)%rsnrs(2)
!
  if (molofrs(rs1).ne.molofrs(rs2)) then
    write(ilog,*) 'Fatal. Encountered intermolecular crosslink in crosslink_refvals(...).&
 & Multiple intermolecular links for the same pair of molecules are not yet supported.'
    call fexit()
  end if
!
  if (crosslink(i)%itstype.eq.1) then
!   get/set old/required bond lengths
    internals(1) = getblen(i2,i3)
    internals(2) = getblen(i3,i4)
    internals(3) = getblen(i4,i5)
    internals(4) = 2.03 ! getblen(i5,i6) ! SG-SG
    internals(5) = getblen(i6,i7)
    internals(6) = getblen(i7,i8)
!
!   get/set old/required bond angles
    internals(7) = PI - getbang(i2,i3,i4)/RADIAN
    internals(8) = PI - getbang(i3,i4,i5)/RADIAN
    internals(9) = PI - 103.0/RADIAN ! PI - getbang(i4,i5,i6)/RADIAN ! CB-SG-SG
    internals(10) = PI - 103.0/RADIAN ! PI - getbang(i5,i6,i7)/RADIAN ! SG-SG-CB
    internals(11) = PI - getbang(i6,i7,i8)/RADIAN ! SG-CB-CA
    internals(12) = PI - getbang(i7,i8,i9)/RADIAN ! CB-CA-N
!
!   measure current (arbitrary) torsions
    dihmat(1) = getztor(i1,i2,i3,i4)/RADIAN ! phi rs2
    dihmat(2) = getztor(i2,i3,i4,i5)/RADIAN ! chi1 rs2
    dihmat(3) = getztor(i3,i4,i5,i6)/RADIAN ! "chi2" rs2
    dihmat(4) = getztor(i4,i5,i6,i7)/RADIAN ! CB-SG-SG-CB
    dihmat(5) = getztor(i5,i6,i7,i8)/RADIAN ! "chi2" rs1 
    dihmat(6) = getztor(i6,i7,i8,i9)/RADIAN ! chi1 rs1
  else
    write(ilog,*) 'Fatal. Encountered unsupported crosslink type in crosslink_refvals(...).&
 & This is most likely an omission bug.'
    call fexit()
  end if
!
end
!
!------------------------------------------------------------------------------------------------------
!
! find nine consecutive atoms defining 6 torsional degrees of freedom which we can use to close the chain
!
subroutine crosslink_picksolu(i,i1,i2,i3,i4,i5,i6,i7,i8,i9,solucount,solmat,ccmat,oldvals,vals,valsz)
!
  use iounit
  use sequen
  use atoms
  use ujglobals
  use fyoc
  use math
!
  implicit none
!
  integer i,solucount,i1,i2,i3,i4,i5,i6,i7,i8,i9,rs1,rs2,azero,picked,ptlst(5),valsz,nvals
  RTYPE solmat(MAXSOLU*4,6),fc,yc,cc(MAXCHI),ccmat(6),oldvals(12),ba,zt,getztor,vals(valsz)
!
  azero = 0
  rs1 = crosslink(i)%rsnrs(1)
  rs2 = crosslink(i)%rsnrs(2)
!
  if (crosslink(i)%itstype.eq.1) then
    do picked=1,solucount
!
      fc = phi(rs2) + RADIAN*solmat(picked,1)
      yc = psi(rs2)
      if (fc.gt.180.0) fc = fc - 360.0
      if (fc.lt.-180.0) fc = fc + 360.0
      call setfy(rs2,fc,yc,azero)
!
      cc(1) = chi(1,rs2) + RADIAN*solmat(picked,2)
      if (cc(1).gt.180.0) cc(i) = cc(1) - 360.0
      if (cc(1).lt.-180.0) cc(i) = cc(1) + 360.0
      call setchi(rs2,cc,azero)
!
      call makexyz_forbb(rs2)
!   
!     generate dummies for positions of SG and CB on rs1
      ptlst(1) = n+20-2
      ba = RADIAN*(PI-oldvals(9))
      zt = RADIAN*(ccmat(3) + solmat(picked,3))
      if (zt.gt.180.0) zt = zt - 360.0
      if (zt.lt.-180.0) zt = zt + 360.0
      call genxyz(ptlst(1),i5,oldvals(4),i4,ba,i3,zt,azero)
      x(i6) = x(ptlst(1))
      y(i6) = y(ptlst(1))
      z(i6) = z(ptlst(1))
      ptlst(2) = n+20-1
      ba = RADIAN*(PI-oldvals(10))
      zt = RADIAN*(ccmat(4) + solmat(picked,4))
      if (zt.gt.180.0) zt = zt - 360.0
      if (zt.lt.-180.0) zt = zt + 360.0
      call genxyz(ptlst(2),i6,oldvals(5),i5,ba,i4,zt,azero)
      x(i7) = x(ptlst(2))
      y(i7) = y(ptlst(2))
      z(i7) = z(ptlst(2))
!
!     now set actual sidechain dihedrals in rs1 and re-build
      cc(1) = getztor(i9,i8,i7,i6)
      call setchi(rs1,cc,azero) ! same as writing directly: chi(1,rs1)+RADIAN*solmat(picked,6)
      call makexyz_forsc(rs1)
      call eval_forstretch(rs1,rs2,vals,valsz,nvals,i)
!
!     undo changes
      fc = phi(rs2) - RADIAN*solmat(picked,1)
      yc = psi(rs2)
      if (fc.gt.180.0) fc = fc - 360.0
      if (fc.lt.-180.0) fc = fc + 360.0
      call setfy(rs2,fc,yc,azero)
      cc(1) = chi(1,rs2) - RADIAN*solmat(picked,2)
      if (cc(1).gt.180.0) cc(i) = cc(1) - 360.0
      if (cc(1).lt.-180.0) cc(i) = cc(1) + 360.0
      call setchi(rs2,cc,azero)
      call makexyz_forbb(rs2)
      chi(1,rs1) = chi(1,rs1) - RADIAN*solmat(picked,6)
      call makexyz_forsc(rs1)
    end do
!   
  else
    write(ilog,*) 'Fatal. Encountered unsupported crosslink type in crosslink_picksolu(...).&
 & This is most likely an omission bug.'
    call fexit()
  end if
!
end
!

