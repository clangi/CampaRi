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
! CONTRIBUTIONS: Xiaoling Wang                                             !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ##################################################
! #                                                # 
! #   the core file for concerted rotation moves:  #
! #   contains all relevant functions, mainly the  #
! #   two versions of the bulky deriv.-calculator  #
! #                                                #
! ##################################################
!
!--------------------------------------------------------------------------------
!
!
subroutine buildGMat(rsi,rsf,gmat)
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
  integer i,j,k,rsi,rsf
  RTYPE cosine,sine,eps,effbond,cphi,cpsi
  RTYPE gmat(MAXCRDOF,MAXCRDOF),jacmat(20,MAXCRDOF),nopl,nipl
  RTYPE alph,phival,phistar,psival,psistar,in_pl(3),or_pl(3)
  RTYPE pos1(3),pos2(3),pos3(3),pos4(3),con1(3),con2(3),c1,c2
!
#ifdef DISABLE_FLOAT
  eps = 0.0000001d0
#else
! the low-precision xyz-coordinates lead to some hysteresis error in converting back and forth
  eps = 0.02d0
#endif
!
! we loop over all the phi/psis involved in the move
  do i=rsi,rsf-1
!
!   we want to find the partials of three different atoms with the torsions above 
    do j=1,3
!
!     first get phi (C-N-CA-C)
      pos4(1) = x(ci(i-1))
      pos4(2) = y(ci(i-1))
      pos4(3) = z(ci(i-1))
      pos3(1) = x(ni(i))
      pos3(2) = y(ni(i))
      pos3(3) = z(ni(i))
      pos2(1) = x(cai(i))
      pos2(2) = y(cai(i))
      pos2(3) = z(cai(i))
      pos1(1) = x(ci(i))
      pos1(2) = y(ci(i))
      pos1(3) = z(ci(i))
      call dihed(pos1,pos2,pos3,pos4,phival)
      cphi = radian*phival
      if (abs(cphi-phi(i)).gt.eps) then
        sj_wrncnt(1) = sj_wrncnt(1) + 1
        if (sj_wrncnt(1).eq.sj_wrnlmt(1)) then
          write(ilog,*) 'Warning. Phi-angle of residue ',i,' does not seem to be computed with&
 & the expected accuracy (abs. deviation is ',abs(cphi-phi(i)),' degrees). This may be the &
 &result of precision-reducing floating-point optimizations by the compiler or a bug.'
          write(ilog,*) 'This was warning #',sj_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*sj_wrnlmt(1).gt.0.5*HUGE(sj_wrnlmt(1))) then
            sj_wrncnt(1) = 0 ! reset
          else
            sj_wrnlmt(1) = sj_wrnlmt(1)*10
          end if
        end if
      end if
!   
!     now get phi* for the current reference atom
      if (j.eq.1) then
        pos1(1) = x(ni(rsf))
        pos1(2) = y(ni(rsf))
        pos1(3) = z(ni(rsf))
      else if (j.eq.2) then
        pos1(1) = x(cai(rsf))
        pos1(2) = y(cai(rsf))
        pos1(3) = z(cai(rsf))
      else
        pos1(1) = x(ci(rsf))
        pos1(2) = y(ci(rsf))
        pos1(3) = z(ci(rsf))
      end if
      call dihed(pos1,pos2,pos3,pos4,phistar)
!
!     now the derivative is simply dref/dphistar * dphistar/dphi
!     - the latter is always unity
!     - the former can be computed from the analytical formula in genxyz
!
!     first get normalized connection vectors 23 and 34 to set up frame
      c1 = 0.0
      c2 = 0.0
      do k=1,3
        con1(k) = pos3(k) - pos2(k)
        c1 = c1 + con1(k)*con1(k)
        con2(k) = pos4(k) - pos3(k)
        c2 = c2 + con2(k)*con2(k)
      end do
      do k=1,3
        con1(k) = con1(k)/sqrt(c1)
        con2(k) = con2(k)/sqrt(c2)
      end do
!     a normal vector to the reference plane is given by 23x34
      call crossprod(con1,con2,or_pl,nopl)
!     normalize it
      cosine = con1(1)*con2(1)+con1(2)*con2(2)+con1(3)*con2(3)
      sine = sqrt(max(1.0-cosine**2,eps))
      if (abs(cosine) .ge. 1.0) then
        write(ilog,10)  i,phival,cosine
 10         format (/,' Bad phi for residue',i6,': ',2f12.4)
      end if
      do k=1,3
        or_pl(k) = or_pl(k)/sine
      end do
!     and the vector within the reference plane is given by (23x34)x34
      call crossprod(or_pl,con2,in_pl,nipl)
!  
!     the effective bond angles and bond lengths are readily calculated
      do k=1,3
        con1(k) = pos2(k) - pos1(k)
      end do
      effbond = sqrt(con1(1)**2 + con1(2)**2 + con1(3)**2)
      call bondang(pos1,pos2,pos3,alph)
!
!     finally, the formula for the derivatives
      do k=1,3
        jacmat(3*(j-1)+k,2*(i-rsi)+1) = &
 &          effbond*(-in_pl(k)*sin(alph)*sin(phistar) +&
 &                    or_pl(k)*sin(alph)*cos(phistar))
      end do
!
!     now the same procedure for psi (NCACN)
!
      pos4(1) = x(ni(i))
      pos4(2) = y(ni(i))
      pos4(3) = z(ni(i))
      pos3(1) = x(cai(i))
      pos3(2) = y(cai(i))
      pos3(3) = z(cai(i))
      pos2(1) = x(ci(i))
      pos2(2) = y(ci(i))
      pos2(3) = z(ci(i))
      pos1(1) = x(ni(i+1))
      pos1(2) = y(ni(i+1))
      pos1(3) = z(ni(i+1))
      call dihed(pos1,pos2,pos3,pos4,psival)
      cpsi = radian*psival
      cpsi = cpsi - psish(i)
      if (cpsi.gt.180.0) cpsi = cpsi - 360.0
      if (cpsi.lt.-180.0) cpsi = cpsi + 360.0
      if (abs(cpsi-psi(i)).gt.eps) then
        sj_wrncnt(1) = sj_wrncnt(1) + 1
        if (sj_wrncnt(1).eq.sj_wrnlmt(1)) then
          write(ilog,*) 'Warning. Psi-angle of residue ',i,' does not seem to be computed with&
 & the expected accuracy (abs. deviation is ',abs(cpsi-psi(i)),' degrees). This may be the &
 &result of precision-reducing floating-point optimizations by the compiler or a bug.'
          write(ilog,*) 'This was warning #',sj_wrncnt(1),' of this type not all of which may be displayed.'
          if (10.0*sj_wrnlmt(1).gt.0.5*HUGE(sj_wrnlmt(1))) then
            sj_wrncnt(1) = 0 ! reset
          else
            sj_wrnlmt(1) = sj_wrnlmt(1)*10
          end if
        end if
      end if
!   
!     now get psi* for the current reference atom
      if (j.eq.1) then
        pos1(1) = x(ni(rsf))
        pos1(2) = y(ni(rsf))
        pos1(3) = z(ni(rsf))
      else if (j.eq.2) then
        pos1(1) = x(cai(rsf))
        pos1(2) = y(cai(rsf))
        pos1(3) = z(cai(rsf))
      else
        pos1(1) = x(ci(rsf))
        pos1(2) = y(ci(rsf))
        pos1(3) = z(ci(rsf))
      end if
      call dihed(pos1,pos2,pos3,pos4,psistar)
!
!     now the derivative is simply dref/dpsistar * dpsistar/dpsi
!     - the latter is always unity
!     - the former can be computed from the analytical formula (compare genxyz(...)
!     -                                                         in topology.f)
!
!     first get normalized connection vectors 23 and 34 to set up frame
      c1 = 0.0
      c2 = 0.0
      do k=1,3
        con1(k) = pos3(k) - pos2(k)
        c1 = c1 + con1(k)*con1(k)
        con2(k) = pos4(k) - pos3(k)
        c2 = c2 + con2(k)*con2(k)
      end do
      do k=1,3
        con1(k) = con1(k)/sqrt(c1)
        con2(k) = con2(k)/sqrt(c2)
      end do
!     a normal vector to the reference plane is given by 23x34
      call crossprod(con1,con2,or_pl,nopl)
!     normalize it
      cosine = con1(1)*con2(1)+con1(2)*con2(2)+con1(3)*con2(3)
      sine = sqrt(max(1.0-cosine**2,eps))
      if (abs(cosine) .ge. 1.0) then
        write(ilog,10)  i,psival,cosine
 11         format (/,' Bad psi for residue',i6,': ',2f12.4)
      end if
      do k=1,3
        or_pl(k) = or_pl(k)/sine
      end do
!     and the vector within the reference plane is given by (23x34)x34
      call crossprod(or_pl,con2,in_pl,nipl)
!
!     the effective bond angles and bond lengths are readily calculated
      do k=1,3
        con1(k) = pos2(k) - pos1(k)
      end do
      effbond = sqrt(con1(1)**2 + con1(2)**2 + con1(3)**2)
      call bondang(pos1,pos2,pos3,alph)
!
!     finally, the formula for the derivatives
      do k=1,3
        jacmat(3*(j-1)+k,2*(i-rsi)+2) = &
 &          effbond*(-in_pl(k)*sin(alph)*sin(psistar) +&
 &                    or_pl(k)*sin(alph)*cos(psistar))
      end do
!
    end do !over the three reference atoms
!
  end do !over all relevant backbone angles
!
!    write(ilog,*) 'Asymmetric Jacobi matrix:'
!    do i=1,nr_crfit
!      write(ilog,27) (jacmat(i,j),j=1,nr_crdof)
!    end do
! 27     format(20f9.4)
!
  do i=1,nr_crdof
    do j=1,nr_crdof
      gmat(i,j) = 0.0
      do k=1,nr_crfit
        gmat(i,j) = gmat(i,j) + jacmat(k,i)*jacmat(k,j)
      end do
    end do
  end do
!
end
!
!----------------------------------------------------------------------------
!
!
!
!
! ###############################################
! ##                                           ##
! ##   buildGMat2 -- concerted rotation matrix ##
! ##                 via nested rot. matrices  ##
! ##                                           ##
! ###############################################
!
! compute matrix G defined in the ref.
! ref: G. Favrin, A.Irback, F.Sjunnesson, 'Monte Carlo update for chain 
! molecules: Biassed Gaussian steps in torsional space', J. Chem. Phys., 
! vol.14 (18), P8154-8158, 2001
! via nested rotation matrices
!
subroutine buildGmat2(rsi,rsf,g)
!
  use iounit
  use math
  use fyoc
  use movesets
!
  implicit none
!
  integer i,j,imod,rsi,rsf,ii,jj
  RTYPE g(MAXCRDOF,MAXCRDOF)
  RTYPE p(8),r1(3),r2(3),r3(3),dr1(8,3),dr2(8,3),dr3(8,3)
  RTYPE mtor(3,3),mang(3,3),mm(3,3),vtem(3),mv(3)
!
! p -- eight torsions subject to move (in degree)
! a -- returned 8-by-8 matrix defined in ref
!
! r1 - position vector of cai of fourth residue
! r2 - position vector of c of fourth residue
! r3 - position vector of O of fourth residue
!
! sanity check and get the torsions into vector p
  if (rsf-rsi.ne.4) then
    write(ilog,*) 'Fatal. Called mode 2 of concerted rotation algori&
 &thm with improper number of DOF: rsi: ',rsi,' rsf: ',rsf
    call fexit()
  end if
!
  do i=rsi,rsf-1
    p(2*(i-rsi)+1) = phi(i)
    p(2*(i-rsi)+2) = psi(i)
  end do
!
! cal. dr/d(phi(8))
  r3(1) = 0.0 
  r3(2) = bco*sin(cco/radian)*sin(p(8)/radian)
  r3(3) = bco*sin(cco/radian)*cos(p(8)/radian)

  do i=7,1,-1
     call rotmax(-p(i),1,mtor)
     imod = mod(i,2)
     if (imod .eq. 0) then
!        mang = rob
        do ii =1,3
          do jj = 1,3
            mang(ii,jj) = rob(ii,jj)
          end do
        end do
     else
!        mang = roa
        do ii =1,3
          do jj = 1,3
            mang(ii,jj) = roa(ii,jj)
          end do
        end do
     end if
     call max3max3(mtor,mang,mm)
     call max3vet(mm,r3,mv)
!     r3 = mv
     do ii=1,3
       r3(ii) = mv(ii)
     end do
     if ((imod .ne. 0) .and. (i .ne. 1)) then
        call max3vet(mcn,r3,mv)
!        r3 = mv
        do ii=1,3
          r3(ii) = mv(ii)
        end do
     end if
  end do      
!  dr1(8,1:3) = (/0.0, 0.0, 0.0/)
!  dr2(8,1:3) = (/0.0, 0.0, 0.0/)
!  dr3(8,1:3) = r3
  do ii=1,3
    dr1(8,ii) = 0.0
    dr2(8,ii) = 0.0
    dr3(8,ii) = r3(ii)
  end do
!  write(ilog,*) 8
!  write(ilog,*) dr1(8,1:3)
!  write(ilog,*) dr2(8,1:3)
!  write(ilog,*) dr3(8,1:3)
 
! cal dr/d(phi(i))
  do i=7,1,-1

! vectors in the local frame centered on Cai of fourth residue
!     r1 = (/0.0,0.0,0.0/)
!     r2 = (/bcc, 0.0, 0.0/)
     do ii=1,3
       r1(ii) = 0.0
       r2(ii) = 0.0
     end do
     r2(1) = bcc
     r3(1) = bcc - bco*cos(cco/radian) 
     r3(2) = -bco*sin(cco/radian)*cos(p(8)/radian)
     r3(3) = bco*sin(cco/radian)*sin(p(8)/radian)
   
     do j=7,i+1,-1
        call rotmax(-p(j),1,mtor)
        imod = mod(j,2)
        if (imod .eq. 0) then
!           mang = rob
           do ii =1,3
             do jj = 1,3
               mang(ii,jj) = rob(ii,jj)
             end do
           end do
!           vtem =  (/bcc, 0.0, 0.0/)
           vtem(1) = bcc
           vtem(2) = 0.0
           vtem(3) = 0.0
        else
!           mang = roa
           do ii =1,3
             do jj = 1,3
               mang(ii,jj) = roa(ii,jj)
             end do
           end do
!           vtem = (/bnc, 0.0, 0.0/)
           vtem(1) = bnc
           vtem(2) = 0.0
           vtem(3) = 0.0
        end if
        call max3max3(mtor,mang,mm)
        call max3vet(mm,r1,mv)            
!        r1 = mv + vtem
        do ii=1,3
          r1(ii) = mv(ii) + vtem(ii)
        end do
        call max3vet(mm,r2,mv)
!        r2 = mv + vtem
        do ii=1,3
          r2(ii) = mv(ii) + vtem(ii)
        end do
        call max3vet(mm,r3,mv)
!        r3 = mv + vtem
        do ii=1,3
          r3(ii) = mv(ii) + vtem(ii)
        end do
!
        if (imod .ne. 0) then
!            vtem = (/bcn, 0.0, 0.0/)
            vtem(1) = bcn
            vtem(2) = 0.0
            vtem(3) = 0.0
            call max3vet(mcn,r1,mv)
!            r1 = vtem + mv
            do ii=1,3
              r1(ii) = mv(ii) + vtem(ii)
            end do
            call max3vet(mcn,r2,mv)
!            r2 = vtem + mv
            do ii=1,3
              r2(ii) = mv(ii) + vtem(ii)
            end do
            call max3vet(mcn,r3,mv)
!            r3 = vtem + mv  
            do ii=1,3
              r3(ii) = mv(ii) + vtem(ii)
            end do
        end if
     end do
!
     do j=i,1,-1
        if (j .eq. i) then
           call drotmax(-p(j),1,mtor)
        else
           call rotmax(-p(j),1,mtor)
        end if
         
        imod = mod(j,2)
        if (imod .eq. 0) then
!           mang = rob
           do ii =1,3
             do jj = 1,3
               mang(ii,jj) = rob(ii,jj)
             end do
           end do
        else
!           mang = roa
           do ii =1,3
             do jj = 1,3
               mang(ii,jj) = roa(ii,jj)
             end do
           end do
        end if

        call max3max3(mtor,mang,mm)
        call max3vet(mm,r1,mv)
!        r1 = mv
        do ii=1,3
          r1(ii) = mv(ii)
        end do
        call max3vet(mm,r2,mv)
!        r2 = mv
        do ii=1,3
          r2(ii) = mv(ii)
        end do
        call max3vet(mm,r3,mv)
!        r3 = mv
        do ii=1,3
          r3(ii) = mv(ii)
        end do
!
        if ((imod .ne. 0) .and. (j .ne. 1)) then
           call max3vet(mcn,r1,mv)
!           r1 = mv
           do ii=1,3
             r1(ii) = mv(ii)
           end do
           call max3vet(mcn,r2,mv)
!           r2 = mv
           do ii=1,3
             r2(ii) = mv(ii)
           end do
           call max3vet(mcn,r3,mv)
!           r3 = mv
           do ii=1,3
             r3(ii) = mv(ii)
           end do
        end if
     end do
     do ii=1,3
       dr1(i,ii) = r1(ii)
       dr2(i,ii) = r2(ii)
       dr3(i,ii) = r3(ii)
     end do
!     do kkk=1,8
!       write(ilog,*) i
!       write(ilog,*) dr1(i,1:3)
!       write(ilog,*) dr2(i,1:3)
!       write(ilog,*) dr3(i,1:3)
!     end do 
  end do

! cal. matrix g
  do i=1,8
     do j=i,8
        do ii=1,3
          vtem(ii) = dr1(i,ii)*dr1(j,ii)
        end do
        g(i,j) = vtem(1) + vtem(2) + vtem(3)
        do ii=1,3
          vtem(ii) = dr2(i,ii)*dr2(j,ii)
        end do
        g(i,j) = g(i,j) + vtem(1) + vtem(2) + vtem(3)
        do ii=1,3
          vtem(ii) = dr3(i,ii)*dr3(j,ii)
        end do
        g(i,j) = g(i,j) + vtem(1) + vtem(2) + vtem(3)
     end do
  end do
! set up other elements of g according to symetry
  do i=1,8
     do j=1,(i-1)
        g(i,j) = g(j,i)
     end do
  end do
!
  end 
!
!-----------------------------------------------------------------------
!
! the following are helper routines which are non-generalizable, i.e.,
! specifically tailored toward being used with the concerted rotation
! approach.
! even though they partially perform basic LA tasks, they are not listed
! in math_utils.f
!
! rotmax -- set up the rotation matrix for rotation alpa angle
! about x or y or z axis
!
subroutine rotmax(alpa,axis,rom)
!
  use math
!
  implicit none
!
  integer axis
  RTYPE alpa,cc,ss,rom(3,3)

! axis == 1,2,3 coresponds to x,y,z axis, respectively
! alpa is the rotated angle, in unit of degree
  cc = cos(alpa/radian)
  ss = sin(alpa/radian)
 
  if (axis .eq. 1) then
     rom(1,1) = 1.0
     rom(1,2) = 0.0
     rom(1,3) = 0.0
     rom(2,1) = 0.0
     rom(2,2) = cc
     rom(2,3) = ss
     rom(3,1) = 0.0
     rom(3,2) = -ss
     rom(3,3) = cc
  else if (axis .eq. 2) then
     rom(1,1) = cc
     rom(1,2) = 0.0
     rom(1,3) = ss
     rom(2,1) = 0.0
     rom(2,2) = 1.0
     rom(2,3) = 0.0
     rom(3,1) = -ss
     rom(3,2) = 0.0
     rom(3,3) = cc
  else if (axis .eq. 3) then
     rom(1,1) = cc
     rom(1,2) = ss
     rom(1,3) = 0.0
     rom(2,1) = -ss
     rom(2,2) = cc
     rom(2,3) = 0.0
     rom(3,1) = 0.0
     rom(3,2) = 0.0
     rom(3,3) = 1.0
  end if
!
end
!  
!
!
! ###################################################
! ##                                               ##
! ##   drotmax -- derivative of rotation matrix    ##
! ##                                               ##
! ###################################################
!
! drotmax -- calculat the derivative of the 
! rotation matrix
! derivative respective to alpa

subroutine drotmax(alpa,axis,drom)
!
  use math
!
  implicit none
!
  integer axis
  RTYPE alpa,cc,ss,drom(3,3)

! axis == 1,2,3 coresponds to x,y,z axis, respectively
! alpa is the rotated angle
  cc = cos(alpa/radian)
  ss = sin(alpa/radian)
 
  if (axis .eq. 1) then
     drom(1,1) = 0.0
     drom(1,2) = 0.0
     drom(1,3) = 0.0
     drom(2,1) = 0.0
     drom(2,2) = -ss
     drom(2,3) = cc
     drom(3,1) = 0.0
     drom(3,2) = -cc
     drom(3,3) = -ss
  else if (axis .eq. 2) then
     drom(1,1) = -ss
     drom(1,2) = 0.0
     drom(1,3) = cc
     drom(2,1) = 0.0
     drom(2,2) = 0.0
     drom(2,3) = 0.0
     drom(3,1) = -cc
     drom(3,2) = 0.0
     drom(3,3) = -ss
  else if (axis .eq. 3) then
     drom(1,1) = -ss
     drom(1,2) = cc
     drom(1,3) = 0.0
     drom(2,1) = -cc
     drom(2,2) = -ss
     drom(2,3) = 0.0
     drom(3,1) = 0.0
     drom(3,2) = 0.0
     drom(3,3) = 0.0
  end if
!
end
!  
!-----------------------------------------------------------------------
!
!
! ############################################################
! ##                                                        ##
! ##   consmax -- constant matries for concerted rotation   ##
! ##                                                        ##
! ############################################################
!
! consmax -- compute constant matrices needed in cal. of matrix G 
!            defined in ref.
! Two adjustable parameters are also sep up
! ref: G. Favrin, A.Irback, F.Sjunnesson, 'Monte Carlo update for chain 
! molecules: Biassed Gaussian steps in torsional space', J. Chem. Phys., 
! vol.14 (18), P8154-8158, 2001
!
subroutine consmax
!
  use iounit
  use movesets
!
  implicit none
!
  RTYPE alph,beta,gama,rog(3,3),mtem(3,3),dummy
!
! alph  -- N-Cai-C bond angle
! beta  -- Cai-C-N bond angle
! gama  -- C-N-Cai bond angle
! roa   -- rotation matrix with alph
! rob   -- rotation matrix with beta
! rog   -- rotation matrix with gama
!
!
! initialize the commonly used bond lengths and angles
  bnc = 1.458d0
  bcc = 1.525d0
  bcn = 1.329d0
  bco = 1.231d0 
  cco = 120.5d0

!
  alph = 110.0
  beta = 116.2d0
  gama = 121.7d0
!
  alph = alph - 180.0
  beta = beta - 180.0
  gama = gama - 180.0
  dummy = -179.5d0
!
  call rotmax(alph,3,roa)
  call rotmax(beta,3,rob)
  call rotmax(gama,3,rog)
!
  call rotmax(dummy,1,mtem)
  call max3max3(mtem,rog,mcn)
!
  end 
!
!-----------------------------------------------------------------------
!
! cdecomp -- Cholesky decompose a symetric and positive matrix
! into a lower and upper triangular matrix
!
! note that this is specifically tailored for concerted rotation moves
!
subroutine cdecomp(M,nsz,nwk,L,trp)
!
  use iounit
  use movesets
!
  implicit none
!
  integer i,j,k,nsz,nwk
  logical trp
  RTYPE tem,M(nsz,nsz),L(nsz,nsz)
  RTYPE bu(nsz,nsz)
!
! M is the matrix subject to decomposition
! L is the lower triangular matrix
!
!  if (n.gt.MAXCRDOF) then
!    write(ilog,*) 'FATAL: Matrix dimensions exceeded in cdecomp.'
!    call fexit()
!  end if
!
! initialize L
  do i=1,nwk
    do j=1,nwk
      L(i,j) = 0.0
    end do
  end do
!
  do i=1,nwk
    tem = 0.0
    do k=1,(i-1)
      tem = tem + L(i,k)*L(i,k)
    end do
    tem = M(i,i) - tem
    if (tem .le.0.0) then 
      write(ilog,*) 'Fatal. Negative value in sqrt for i=',i,' during Cholesky decomposition.'
      call fexit()
    end if
    L(i,i) = sqrt(tem)
    if (abs(L(i,i)-0.0) .lt. 1.0d-8) then
      write(ilog,*) 'Warning. Zero entry in diagonal of Cholesky matrix during Cholesky decomposition.'
    end if
!     
    do j=(i+1),nwk
      tem = 0.0
      do k=1,(i-1)
        tem = tem +  L(i,k)*L(j,k)
      end do
      L(j,i) = (M(j,i) - tem)/L(i,i)
    end do
  end do
!
  if (trp.EQV..true.) then
    do i=1,nwk
      do j=1,nwk
        bu(i,j) = L(i,j)
      end do
    end do
    do i=1,nwk
      do j=1,nwk
        L(i,j) = bu(j,i)
      end do
    end do
  end if
!  write(ilog,*) 'Cholesky-decomposed matrix.'
!  do i=1,n
!      write(*,10) L(i,1:n)
!  end do
! 10  format(20f9.4)
  return
end
!
!--------------------------------------------------------------------------
!
! solves a triangular system of equations given by lmat*dphi = gvec
! lmat is upper triangle (nxn), dphi,gvec are (nx1)  
!
! note that this is specifically tailored for concerted rotation moves
!
subroutine trisolv_sj(lmat,dphi,gvec,nsz,nwk)
!
  use iounit
  use math
  use movesets
!
  implicit none
!
  integer nsz
  RTYPE lmat(nsz,nsz),gvec(nsz),random
  RTYPE dphi(nsz),eps
  integer nwk,i,j
!
! stupid sanity check
  eps = 0.0000001d0
  if (abs(lmat(nwk,1)).gt.eps) then
    write(ilog,*) 'Got the wrong triangular matrix for solver.'
    write(ilog,*) 'Please use transpose.'
    call fexit()
  end if
!
! first generate a vector of Gaussian-random numbers
  do i=1,nwk
    gvec(i) = sqrt(-log(random()))*cos(2*PI*random())
  end do
!  do i=1,n
!    gvec(i) = random()*(PI/50.0) - (PI/100.0)
!    write(ilog,*) gvec(i)**2
!  end do
!
! now solve the system of equations: lmat*dphi = gvec 
  do i=nwk,1,-1
    dphi(i) = gvec(i)
    do j=nwk,i+1,-1
      dphi(i) = dphi(i) - dphi(j)*lmat(i,j)
    end do
    dphi(i) = dphi(i)/lmat(i,i)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! solves a triangular system of equations given by lmat*dphi = gvec
! lmat is upper triangle (nxn), dphi,gvec are (nx1)
! nsz is allocation size for arrays, nwk is assumed size for algorithm 
!
subroutine trisolv(lmat,dphi,gvec,nsz,nwk)
!
  use iounit
  use math
  use movesets
!
  implicit none
!
  integer nsz
  RTYPE lmat(nsz,nsz),gvec(nsz)
  RTYPE dphi(nsz),eps
  integer nwk,i,j
!
! stupid sanity check
  eps = 0.0000001d0
  if ((nwk.gt.1).AND.(abs(lmat(nwk,1)).gt.eps)) then
    write(ilog,*) 'Got the wrong triangular matrix for solver.'
    write(ilog,*) 'Please use transpose.'
    call fexit()
  end if
!
! now solve the system of equations: lmat*dphi = gvec 
  do i=nwk,1,-1
    dphi(i) = gvec(i)
    do j=nwk,i+1,-1
      dphi(i) = dphi(i) - dphi(j)*lmat(i,j)
    end do
    dphi(i) = dphi(i)/lmat(i,i)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this creates lists with eligible final residues and minimal first residues for
! CR stretches of all types
!
subroutine setup_CR_eligible(mode,rslst,rslst2)
!
  use movesets
  use sequen
  use molecule
  use fyoc
  use ujglobals
  use iounit
  use polypep
  use zmatrix
!
  implicit none
!
  integer minsegl,rslst(nseq),rslst2(nseq),mode,rs
!
! assemble the list for eligible final residues in stretches for NUCCR
  if ((mode.eq.0).OR.(mode.eq.3)) then
    nuccrlst%nr = 0
    minsegl = 0
    do rs=1,nseq
      if (seqpolty(rs).ne.'N') then
        minsegl = 0
      else if ((seqflag(rs).eq.24).AND.(rs.eq.rsmol(molofrs(rs),1))) then
        minsegl = 0
      else if ((seqpolty(rs).eq.'N').AND.(rs.eq.rsmol(molofrs(rs),1))) then
        minsegl = 1
      else
        minsegl = minsegl + 1
      end if
      if (minsegl.ge.2) then
        nuccrlst%nr = nuccrlst%nr + 1
        rslst(nuccrlst%nr) = rs
        rslst2(nuccrlst%nr) = rs - minsegl + 1
        if ((rs-minsegl).eq.rsmol(molofrs(rs),1)) then
          if (seqflag(rs-minsegl).eq.24) then
            rslst2(nuccrlst%nr) = rs - minsegl
          end if
        end if
      end if
    end do
    if (nuccrlst%nr.gt.0) then
      if (allocated(nuccrlst%idx).EQV..true.) deallocate(nuccrlst%idx)
      if (allocated(nuccrlst%idx2).EQV..true.) deallocate(nuccrlst%idx2)
      allocate(nuccrlst%idx(nuccrlst%nr))
      allocate(nuccrlst%idx2(nuccrlst%nr))
      nuccrlst%idx(:) = rslst(1:nuccrlst%nr)
      nuccrlst%idx2(:) = rslst2(1:nuccrlst%nr)
    end if
  end if
  if ((mode.eq.0).OR.(mode.eq.4)) then
!   the same for UJCR
    ujcrlst%nr = 0
    minsegl = 0
    do rs=1,nseq
!     currently, these can only be intra-chain residues that are not cyclic
      if ((fline(rs).le.0).OR.(yline(rs).le.0).OR.(seqflag(rs).eq.5).OR.&
   &(rs.eq.rsmol(molofrs(rs),1)).OR.(rs.eq.rsmol(molofrs(rs),2)).OR.(seqpolty(rs).ne.'P')) then
        minsegl = 0
      else if ((izrot(fline(rs))%alsz.le.0).OR.(izrot(yline(rs))%alsz.le.0)) then
        minsegl = 0
      else
        minsegl = minsegl + 1
      end if
      if (minsegl.ge.(ujminsz)) then
        ujcrlst%nr = ujcrlst%nr + 1
        rslst(ujcrlst%nr) = rs
        rslst2(ujcrlst%nr) = rs - minsegl + 1
      end if
    end do
    if (ujcrlst%nr.gt.0) then
      if (allocated(ujcrlst%idx).EQV..true.) deallocate(ujcrlst%idx)
      if (allocated(ujcrlst%idx2).EQV..true.) deallocate(ujcrlst%idx2)
      allocate(ujcrlst%idx(ujcrlst%nr))
      allocate(ujcrlst%idx2(ujcrlst%nr))
      ujcrlst%idx(:) = rslst(1:ujcrlst%nr)
      ujcrlst%idx2(:) = rslst2(1:ujcrlst%nr)
    end if
  end if
  if ((mode.eq.0).OR.(mode.eq.1)) then
!   the same for DJCR
    djcrlst%nr = 0
    minsegl = 0
    do rs=1,nseq
      if (seqpolty(rs).ne.'P') then
        if (seqflag(rs).eq.10) then
          minsegl = 1
        else
          minsegl = 0
        end if
      else if ((fline(rs).le.0).OR.(yline(rs).le.0)) then
        if (seqflag(rs).eq.10) then
          minsegl = 1
        else
          minsegl = 0
        end if
      else if (((izrot(fline(rs))%alsz.le.0).AND.(pucline(rs).le.0)).OR.(izrot(yline(rs))%alsz.le.0)) then
        minsegl = 0
      else if ((seqflag(rs).eq.10).OR.((seqpolty(rs).eq.'P').AND.(rs.eq.rsmol(molofrs(rs),1)))) then
        minsegl = 1
      else
        minsegl = minsegl + 1
      end if
      if (minsegl.ge.4) then
        djcrlst%nr = djcrlst%nr + 1
        rslst(djcrlst%nr) = rs
        rslst2(djcrlst%nr) = rs - minsegl + 1
      end if
    end do
    if (djcrlst%nr.gt.0) then
      if (allocated(djcrlst%idx).EQV..true.) deallocate(djcrlst%idx)
      if (allocated(djcrlst%idx2).EQV..true.) deallocate(djcrlst%idx2)
      allocate(djcrlst%idx(djcrlst%nr))
      allocate(djcrlst%idx2(djcrlst%nr))
      djcrlst%idx(:) = rslst(1:djcrlst%nr)
      djcrlst%idx2(:) = rslst2(1:djcrlst%nr)
    end if
  end if
  if ((mode.eq.0).OR.(mode.eq.2)) then
!   the same for DOCR
    docrlst%nr = 0
    minsegl = 0
    do rs=1,nseq
      if (seqpolty(rs).ne.'P') then
        minsegl = 0
      else if ((fline(rs).le.0).OR.(yline(rs).le.0)) then
        minsegl = 0
      else if (((izrot(fline(rs))%alsz.le.0).AND.(pucline(rs).le.0)).OR.(izrot(yline(rs))%alsz.le.0)) then
        minsegl = 0
      else if ((seqpolty(rs).eq.'P').AND.(rs.eq.rsmol(molofrs(rs),1))) then
        minsegl = 1
      else
        minsegl = minsegl + 1
      end if
      if (minsegl.ge.3) then
        docrlst%nr = docrlst%nr + 1
        rslst(docrlst%nr) = rs
        rslst2(docrlst%nr) = rs - minsegl + 1
      end if
    end do
    if (docrlst%nr.gt.0) then
      if (allocated(docrlst%idx).EQV..true.) deallocate(docrlst%idx)
      if (allocated(docrlst%idx2).EQV..true.) deallocate(docrlst%idx2)
      allocate(docrlst%idx(docrlst%nr))
      allocate(docrlst%idx2(docrlst%nr))
      docrlst%idx(:) = rslst(1:docrlst%nr)
      docrlst%idx2(:) = rslst2(1:docrlst%nr)
    end if
  end if
  if ((mode.eq.0).OR.(mode.eq.5)) then
!   the same for SJCR
    if ((cr_mode.eq.2).AND.(nr_crres.ne.4)) then
      write(ilog,*) 'Warning. With FMCSC_CRMODE 2, stretch length for Sjunnesson-Favrin polypeptide concerted &
 &rotation moves is fixed to 4 residues. Overwriting ...' 
      nr_crdof = 8
      nr_crres = 4
    end if
    sjcrlst%nr = 0
    minsegl = 0
    do rs=1,nseq
!     currently, these can only be intra-chain residues that are not cyclic, in mode 1 extra buffer required
      if ((fline(rs).le.0).OR.(yline(rs).le.0).OR.(seqflag(rs).eq.5).OR.&
   &(rs.eq.rsmol(molofrs(rs),1)).OR.(rs.eq.rsmol(molofrs(rs),2)).OR.(seqpolty(rs).ne.'P')) then
        minsegl = 0
      else if ((izrot(fline(rs))%alsz.le.0).OR.(izrot(yline(rs))%alsz.le.0)) then
        minsegl = 0
      else if ((cr_mode.eq.1).AND.((rs.eq.(rsmol(molofrs(rs),1)+1)).OR.(rs.eq.(rsmol(molofrs(rs),2)-1)))) then
        minsegl = 0
      else
        minsegl = minsegl + 1
      end if
      if (minsegl.ge.nr_crres) then
        sjcrlst%nr = sjcrlst%nr + 1
        rslst(sjcrlst%nr) = rs
        rslst2(sjcrlst%nr) = rs - minsegl + 1
      end if
    end do
    if (sjcrlst%nr.gt.0) then
      if (allocated(sjcrlst%idx).EQV..true.) deallocate(sjcrlst%idx)
      if (allocated(sjcrlst%idx2).EQV..true.) deallocate(sjcrlst%idx2)
      allocate(sjcrlst%idx(sjcrlst%nr))
      allocate(sjcrlst%idx2(sjcrlst%nr))
      sjcrlst%idx(:) = rslst(1:sjcrlst%nr)
      sjcrlst%idx2(:) = rslst2(1:sjcrlst%nr)
    end if
  end if
! 
end
!
!-----------------------------------------------------------------------------------------------
!
