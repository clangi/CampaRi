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
! CONTRIBUTIONS: Hoang Tran, Rohit Pappu                                   !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!------------------------------------------------------------------------
!
!
!             ##################################
!             #                                #
!             # ZEROTH: GLOBAL ENERGY ROUTINES #
!             #                                #
!             ##################################
!
!------------------------------------------------------------------------
!
! routine for the smooth DSSP-based global potential energy term
! note that beta-bridges can very well stretch multiple molecules such
! that the DSSP-assignment is best done globally
!
subroutine en_dssp_gl(evec)
!
  use molecule
  use iounit
  use dssps
  use energies
!
  implicit none
!
  RTYPE esc,wesc,hsc,whsc,evec(MAXENERGYTERMS),esch,hsch
  RTYPE escm(npepmol),wescm(npepmol),hscm(npepmol),whscm(npepmol)
  logical atrue
!
  atrue = .true.
  call update_dssp(esc,wesc,esch,hsc,whsc,hsch,escm,wescm,hscm,whscm,atrue)
!
  evec(19) = evec(19) + scale_DSSP*&
 &            (par_DSSP(8)*(whsc - par_DSSP(7))**2 +&
 &             par_DSSP(10)*(wesc - par_DSSP(9))**2)
!
end
!
!------------------------------------------------------------------------
!
! routine for the smooth density restraint potential
!
subroutine en_emicro_gl(evec,mode,tpi)
!
  use ems
  use energies
  use units
  use atoms
  use iounit
  use mcsums
  use system
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,mode
!
  RTYPE evec(MAXENERGYTERMS)
  RTYPE accum,normer,unitvol,normer2
  integer ints(3),xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,lidx,intsm1(3)
  logical doline,doslice,doblocks,thelog
  integer sta,sto,sta2,sto2,sta3,sto3
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
#endif
!
  doslice = .false.
  if (allocated(emcoarsegrid%lslc).EQV..true.) then
    doslice = .true.
  end if
  doline = .false.
  if (allocated(emcoarsegrid%llin).EQV..true.) then
    doline = .true.
  end if
  doblocks = .false.
  if (allocated(emcoarsegrid%lblk).EQV..true.) then
    doblocks = .true.
  end if
!
#ifdef ENABLE_THREADS
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, force_emicro_gl(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    sta = thr_limits(51,tpi)
    sto = thr_limits(52,tpi)
    sta2 = thr_limits(49,tpi)
    sto2 = thr_limits(50,tpi)
    sta3 = thr_limits(57,tpi)
    sto3 = thr_limits(58,tpi)
  else
    if (doblocks.EQV..true.) then
      sta = 1
      sto = emcoarsegrid%dimc(3)
    else
      sta = 1
      sto = emcoarsegrid%dim(2)
    end if
    sta2 = 1
    sto2 = emgrid%dim(3)
    sta3 = 1
    sto3 = emcoarsegrid%dim(3)
  end if
#else
  if (doblocks.EQV..true.) then
    sta = 1
    sto = emcoarsegrid%dimc(3)
  else
    sta = 1
    sto = emcoarsegrid%dim(2)
  end if
  sta2 = 1
  sto2 = emgrid%dim(3)
  sta3 = 1
  sto3 = emcoarsegrid%dim(3)
#endif
!
  unitvol = emgrid%deltas(1)*emgrid%deltas(2)*emgrid%deltas(3)
  if (pdb_analyze.EQV..true.) curmassm = .false.
!
  thelog = .false.
  if ((mode.eq.2).OR.((mode.eq.0).AND.(curmassm.EQV..false.))) then
    if (mode.eq.2) thelog = .true.
#ifdef ENABLE_THREADS
    if (tpi.le.0) then
      call em_pop_mv(emgrid,thelog)
    else
!     ends with a barrier
      call em_pop_mv_threads(thelog,tpi)
    end if
#else
    call em_pop_mv(emgrid,thelog)
#endif
    if (mode.eq.0) curmassm = .true.
  else if ((mode.eq.1).OR.(mode.eq.0)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Called en_emicro_gl(...) with unsupported mode (',mode,').&
 & This is a bug.'
    call fexit()
  end if
!
  if (empotmode.eq.1) then
!
    normer = convdens/unitvol
    if (mode.le.1) then
      emgrid%dens(:,:,sta2:sto2) = normer*emgrid%mass(:,:,sta2:sto2)
    else if (mode.eq.2) then
      emgrid%dens(:,:,sta2:sto2) = normer*emgrid%mass2(:,:,sta2:sto2)
    end if
!
  else if (empotmode.eq.2) then
!
    if (mode.lt.0) then
      if (emgrid%cnt.gt.0) then
        normer = convdens/(unitvol*(emgrid%cnt))
        emgrid%dens(:,:,sta2:sto2) = normer*(emgrid%avgmass(:,:,sta2:sto2))
      else
        normer = convdens/unitvol
        emgrid%dens(:,:,sta2:sto2) = normer*emgrid%mass(:,:,sta2:sto2)
      end if
    else if (mode.le.1) then
      normer = convdens/unitvol
      if (emgrid%cnt.le.0) then
        emgrid%dens(:,:,sta2:sto2) = normer*emgrid%mass(:,:,sta2:sto2)
      else
        normer2 = (1.0-emiw)/(1.0*emgrid%cnt)
        emgrid%dens(:,:,sta2:sto2) = normer*(normer2*emgrid%avgmass(:,:,sta2:sto2)+emiw*emgrid%mass(:,:,sta2:sto2))
      end if
    else if (mode.eq.2) then
      normer = convdens/unitvol
      if (emgrid%cnt.le.0) then
        emgrid%dens(:,:,sta2:sto2) = normer*emgrid%mass2(:,:,sta2:sto2)
      else
        normer = convdens/unitvol
        normer2 = (1.0-emiw)/(1.0*emgrid%cnt)
        emgrid%dens(:,:,sta2:sto2) = normer*(normer2*emgrid%avgmass(:,:,sta2:sto2)+emiw*emgrid%mass2(:,:,sta2:sto2))
      end if
    end if
!
  else
    write(ilog,*) 'Fatal. Called en_emicro_gl(...) with unsupported potential flag (',empotmode,').&
 & This is a bug.'
    call fexit()
  end if
!
  ints(:) = nint(emingrid%deltas(:)/(emgrid%deltas(:))) 
  intsm1(:) = ints(:) - 1
  normer = 1.0/(ints(1)*ints(2)*ints(3))
  accum = 0.0
  if (mode.eq.2) then
    lidx = 2
  else
    lidx = 1
  end if
  emcoarsegrid%deltadens(:,:,sta3:sto3) = 0. ! may remain incomplete if heuristics are in use
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if (doblocks.EQV..true.) then
    do zk=sta,sto
      do yk=1,emcoarsegrid%dimc(2)
        do xk=1,emcoarsegrid%dimc(1)
          if (emcoarsegrid%lblk(xk,yk,zk,lidx).EQV..false.) cycle
          do yj=emcoarsegrid%blklms(yk,1,2),emcoarsegrid%blklms(yk,2,2)
            yi = yj*ints(2) - intsm1(2)
            do yl=0,ints(2)-1
              do zj=emcoarsegrid%blklms(zk,1,3),emcoarsegrid%blklms(zk,2,3)
                zi = zj*ints(3) - intsm1(3)
                do zl=0,ints(3)-1
                  do xj=emcoarsegrid%blklms(xk,1,1),emcoarsegrid%blklms(xk,2,1)
                    xi = xj*ints(1) - intsm1(1)
                    do xl=0,ints(1)-1
                      emcoarsegrid%deltadens(xj,yj,zj) = emcoarsegrid%deltadens(xj,yj,zj) + emgrid%dens(xi+xl,yi+yl,zi+zl)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    do zk=sta,sto
      do yk=1,emcoarsegrid%dimc(2)
        do xk=1,emcoarsegrid%dimc(1)
          if (emcoarsegrid%lblk(xk,yk,zk,lidx).EQV..false.) then
            accum = accum + emcoarsegrid%rblk(xk,yk,zk)
            cycle
          end if
          do yj=emcoarsegrid%blklms(yk,1,2),emcoarsegrid%blklms(yk,2,2)
            do zj=emcoarsegrid%blklms(zk,1,3),emcoarsegrid%blklms(zk,2,3)
              emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj) = &
 &                 normer*emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj) - &
 &                 emcoarsegrid%dens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj)
              accum = accum + sum(emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj)*&
 &                                emcoarsegrid%deltadens(emcoarsegrid%blklms(xk,1,1):emcoarsegrid%blklms(xk,2,1),yj,zj))
            end do
          end do
        end do
      end do
    end do
  else
    do yj=sta,sto ! 1,emcoarsegrid%dim(2)
      if (doslice.EQV..true.) then
        if (emcoarsegrid%lslc(yj,lidx).EQV..false.) cycle
      end if
      yi = yj*ints(2) - intsm1(2)
      do yk=0,ints(2)-1
        do zj=1,emcoarsegrid%dim(3)
          if (doline.EQV..true.) then
            if (emcoarsegrid%llin(yj,zj,lidx).EQV..false.) cycle
          end if
          zi = zj*ints(3) - intsm1(3)
          do zk=0,ints(3)-1
            do xj=1,emcoarsegrid%dim(1)
              xi = xj*ints(1) - intsm1(1)
              do xk=0,ints(1)-1
                emcoarsegrid%deltadens(xj,yj,zj) = emcoarsegrid%deltadens(xj,yj,zj) + emgrid%dens(xi+xk,yi+yk,zi+zk)
              end do
            end do
          end do
        end do
      end do
    end do
    do yj=sta,sto ! 1,emcoarsegrid%dim(2)
      if (doslice.EQV..true.) then
        if (emcoarsegrid%lslc(yj,lidx).EQV..false.) then
          accum = accum + emcoarsegrid%rslc(yj)
          cycle
        end if
      end if
      do zj=1,emcoarsegrid%dim(3)
        if (doline.EQV..true.) then
          if (emcoarsegrid%llin(yj,zj,lidx).EQV..false.) then
            accum = accum + emcoarsegrid%rlin(yj,zj)
            cycle
          end if
        end if
        emcoarsegrid%deltadens(:,yj,zj) = normer*emcoarsegrid%deltadens(:,yj,zj) - emcoarsegrid%dens(:,yj,zj)
        accum = accum + sum(emcoarsegrid%deltadens(:,yj,zj)*emcoarsegrid%deltadens(:,yj,zj))
      end do
    end do
  end if
!
  evec(21) = evec(21) + scale_EMICRO*accum
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!---------------------------------------------------------------------------------
!
!! routine for the potential based on linear combinations of system torsions
!!
!subroutine en_lc_tor(evec)
!!
!  use energies
!  use torsn
!  use aminos
!  use sequen
!  use fyoc
!  use math
!  use movesets
!!
!  implicit none
!!
!  integer i,j,dc,rs,bin
!  RTYPE evec(MAXENERGYTERMS),lcti
!!
!  dc = 0
!  do rs=1,nseq
!    if (wline(rs).gt.0) then
!      dc = dc + 1
!      curtvec(dc) = omega(rs)/RADIAN
!    end if
!    if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
!      dc = dc + 1
!      curtvec(dc) = phi(rs)/RADIAN
!    end if
!    if (yline(rs).gt.0) then
!      dc = dc + 1
!      curtvec(dc) = psi(rs)/RADIAN
!    end if
!    do i=1,nnucs(rs)
!      dc = dc + 1
!      curtvec(dc) = nucs(i,rs)/RADIAN
!    end do
!    do i=1,nchi(rs)
!      dc = dc + 1
!      curtvec(dc) = chi(i,rs)/RADIAN
!    end do
!  end do
!!
!  do i=1,ntorlcs
!    if (use_lctmoves.EQV..false.) then
!      lcti = 0.0
!      do j=1,ntorsn
!        if (torlcmode.eq.1) then
!          lcti = lcti + torlc_coeff(i,2*j-1)*cos(curtvec(j)) &
! &                  + torlc_coeff(i,2*j)*sin(curtvec(j))
!        else if (torlcmode.eq.2) then
!          lcti = lcti + torlc_coeff(i,j)*curtvec(j)
!        end if
!      end do
!    else
!      lcti = lct(i)
!    end if
!    bin = floor((lcti-par_LCTOR(2))/par_LCTOR(1)) + 1
!    if (bin.lt.1) then
!      evec(11) = evec(11) + 1.0*scale_LCTOR*lct_weight(i)*&
! &                                              lct_pot(i,1)
!    else if (bin.gt.par_LCTOR2(1)) then
!      evec(11) = evec(11) + 1.0*scale_LCTOR*lct_weight(i)*&
! &                                        lct_pot(i,par_LCTOR2(1))
!    else
!      evec(11) = evec(11) + 1.0*scale_LCTOR*lct_weight(i)*&
! & (1.0/par_LCTOR(1))*((lcti - lct_val(bin))*lct_pot(i,bin+1) + &
! &                     (lct_val(bin+1)-lcti)*lct_pot(i,bin))
!    end if
!  end do
!!
!end
!
!------------------------------------------------------------------------
!
!
!           #########################################
!           #                                       #
!           # FIRST: MOLECULE-BASED ENERGY ROUTINES #
!           #                                       #
!           #########################################
!
!
!------------------------------------------------------------------------
!
! routine for the order parameter-based global potential energy terms
!
subroutine en_zsec_gl(imol,evec)
!
  use energies
  use torsn
  use molecule
!
  implicit none
!
  integer imol
  RTYPE evec(MAXENERGYTERMS)
!
  if (do_tors(moltypid(imol)).EQV..false.) return
!
  call z_secondary(imol,z_alpha(imol),z_beta(imol))
!
  evec(8) = evec(8) + scale_ZSEC*&
 &            (par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))**2 +&
 &             par_ZSEC(4)*(z_beta(imol) - par_ZSEC(3))**2)
!
end
!
!-----------------------------------------------------------------------
!
! a simple derivative with respect to z_alpha_0 or z_beta_0
! assumes z_alpha and z_beta are set(!)
!
subroutine der_zsec_gl(imol,deriv,which)
!
  use energies
  use iounit
  use torsn
  use molecule
!
  implicit none
!
  integer imol,which
  RTYPE deriv
!
  if (do_tors(moltypid(imol)).EQV..false.) return
!
  if (which.eq.1) then
    deriv = deriv -2.0*scale_ZSEC*&
 &            par_ZSEC(2)*(z_alpha(imol) - par_ZSEC(1))
  else if (which.eq.2) then
    deriv = deriv -2.0*scale_ZSEC*&
 &            par_ZSEC(4)*(z_beta(imol) - par_ZSEC(3))
  else
    write(ilog,*) 'Fatal. Called der_zsec_gl(...) with unknown mode &
 &(offending mode is ',which,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------
!
! routine for the polymer property-based global potential energy terms
!
subroutine en_poly_gl(imol,evec,mode)
!
  use iounit
  use energies
  use polyavg
  use molecule
!
  implicit none
!
  integer imol,mode
  RTYPE evec(MAXENERGYTERMS),dlts,t1n,t2n,tps
!
  if (par_POLY2(imol).le.0) return
!
  if (atmol(imol,2).le.atmol(imol,1)) return
!
  if (mode.eq.0) then
    t1n = rgevs(imol,1)*rgevs(imol,2) + rgevs(imol,2)*rgevs(imol,3) &
 &                + rgevs(imol,1)*rgevs(imol,3)
    t2n = (rgevs(imol,1) + rgevs(imol,2) + rgevs(imol,3))**2
    dlts = 1.0 - 3.0*(t1n/t2n)
    tps = 2.5*(1.75*rgv(imol)/molcontlen(moltypid(imol)))&
 &          **(4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)
!
    evec(14) = evec(14) + scale_POLY*&
 &            (par_POLY(imol,3)*(tps - par_POLY(imol,1))**2 +&
 &             par_POLY(imol,4)*(dlts - par_POLY(imol,2))**2)
  else if (mode.eq.1) then
    t1n =           rgevsref(imol,1)*rgevsref(imol,2)&
 &                + rgevsref(imol,2)*rgevsref(imol,3)&
 &                + rgevsref(imol,1)*rgevsref(imol,3)
    t2n = (rgevsref(imol,1) + rgevsref(imol,2) +rgevsref(imol,3))**2
    dlts = 1.0 - 3.0*(t1n/t2n)
    tps = 2.5*(1.75*rgvref(imol)/molcontlen(moltypid(imol)))&
 &          **(4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)
!
    evec(14) = evec(14) + scale_POLY*&
 &            (par_POLY(imol,3)*(tps - par_POLY(imol,1))**2&
 &          +  par_POLY(imol,4)*(dlts - par_POLY(imol,2))**2)
  else
    write(ilog,*) 'Fatal. Called en_poly_gl(...) with unknown mode (&
 &',mode,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a simple derivative with respect to t_0 and dlts_0
! assumes rgevs and rgmol are set to current values
!
subroutine der_poly_gl(imol,deriv,which)
!
  use energies
  use iounit
  use polyavg
  use molecule
!
  implicit none
!
  integer imol,which
  RTYPE deriv,t1n,t2n,dlts,tps
!
  if (do_pol(moltypid(imol)).EQV..false.) return
!
  if (which.eq.1) then
    tps = 2.5*(1.75*rgv(imol)/molcontlen(moltypid(imol)))&
 &          **(4.0/(dble(rsmol(imol,2)-rsmol(imol,1)+1))**tp_exp)
    deriv = deriv - 2.0*scale_POLY*&
 &            par_POLY(imol,3)*(tps - par_POLY(imol,1))
  else if (which.eq.2) then
    t1n = rgevs(imol,1)*rgevs(imol,2) + rgevs(imol,2)*rgevs(imol,3) &
 &                + rgevs(imol,1)*rgevs(imol,3)
    t2n = (rgevs(imol,1) + rgevs(imol,2) + rgevs(imol,3))**2
    dlts = 1.0 - 3.0*(t1n/t2n)
    deriv = deriv - 2.0*scale_POLY*&
 &            par_POLY(imol,4)*(dlts - par_POLY(imol,2))
  else
    write(ilog,*) 'Fatal. Called der_poly_gl(...) with unknown mode &
 &(offending mode is ',which,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------
!
!
!
!             #########################################
!             #                                       #
!             # SECOND: RESIDUE-BASED ENERGY ROUTINES #
!             #                                       #
!             #########################################
!
!
!
!--------------------------------------------------------------------------
!
! this is a collection routine for all non-standard terms which are required
! to generate chemically reasonable structures, which are determined by effects
! outside of the scope of the standard NB terms
!
subroutine e_corrector(rs,evec)
!
  use energies
  use aminos
  use sequen
  use fyoc
  use movesets
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS),eomega,ephenol!,epucker
  integer rs
  character(3) resname
!
  resname = amino(seqtyp(rs))
!
! omega torsions (i.e., Z/E amide isomerization)
  if ((omegafreq.gt.0.0).AND.(wline(rs).gt.0)) then
   evec(2) = evec(2) + scale_CORR*eomega(rs,resname)
  end if
!
! proline puckering
!  if (resname.eq.'PRO') then
!    evec(2) = evec(2) + scale_CORR*epucker(rs)
!  end if
!
! phenolic polar hydrogen Z/E
  if ((resname.eq.'PCR').OR.(resname.eq.'TYR')) then
    evec(2) = evec(2) + scale_CORR*ephenol(rs,resname)
  end if
!  write(*,*) chi(3,rs),evec(2)
!
  end 
! 
!
!---------------------------------------------------------------------------
!
function fs_ipol(ati)
!
  use atoms
  use iounit
  use energies
  use math, ONLY: PI
!
  implicit none
!
  RTYPE fs_ipol,rat,hlp0,hlp60,hlp61,hlp1
  integer ati,stepi
!
  rat = atsav(ati)/atbvol(ati)
  if (rat.gt.atsavmaxfr(ati)) then
    fs_ipol = 1.0
  else if (rat.gt.par_IMPSOLV(5)) then
    if (par_IMPSOLV2(2).eq.1) then
      fs_ipol = 1.0/(1.0 + exp(-(rat-atsavprm(ati,1))/par_IMPSOLV(3)))
      fs_ipol = (fs_ipol-0.5)*atsavprm(ati,2) + atsavprm(ati,3)
    else if (par_IMPSOLV2(2).eq.2) then
      stepi = floor((rat-par_IMPSOLV(5))/atsavprm(ati,7)) + 1 ! which increment
      hlp0 = par_IMPSOLV(5)+(stepi-1)*atsavprm(ati,7) ! target lower plateau value in SAV fraction
      hlp60 = 1.0/(1.0 + exp(-(hlp0-atsavprm(ati,1))/par_IMPSOLV(3)))
      hlp60 = (hlp60-0.5)*atsavprm(ati,2) + atsavprm(ati,3) ! target lower plateau value in upsilon
      hlp61 = 1.0/(1.0 + exp(-(hlp0+atsavprm(ati,7)-atsavprm(ati,1))/par_IMPSOLV(3)))
      hlp61 = (hlp61-0.5)*atsavprm(ati,2) + atsavprm(ati,3) ! target next higher plateau value in upsilon
      hlp1 = hlp0 + par_IMPSOLV(13)*par_IMPSOLV(14)*atsavprm(ati,7)
      if (rat.lt.hlp1) then
        fs_ipol = hlp60 
      else if (rat.gt.(hlp0+(1.0-(1.0-par_IMPSOLV(14))*par_IMPSOLV(13))*atsavprm(ati,7))) then
        fs_ipol = hlp61
      else ! par_IMPSOLV(13) must not be exactly 1.0
        fs_ipol = hlp60 + 0.5*(hlp61-hlp60)*(1.0 - cos(((1./(1.-par_IMPSOLV(13)))*PI/atsavprm(ati,7))*(rat-hlp1)))
      end if
    end if
  else
    fs_ipol = 0.0
  end if
!
  return
!
end
!
function fos_Tdep(tdvs)
!
  use system
  use fos
!
  implicit none
!
  RTYPE tdvs(3),fos_Tdep,Trat
!
  if (use_Tdepfosvals.EQV..true.) then
    Trat = kelvin/fos_Tdepref
    fos_Tdep = Trat*(tdvs(1)-tdvs(2)) + tdvs(2) + 0.001*tdvs(3)*(kelvin*(1.0-log(Trat))-fos_Tdepref)
  else
    fos_Tdep = tdvs(1)
  end if
!
end
!
!
! a routine "computing" the hydration free energy for an individual residue based
! on its accessible volume (individual solvation groups separate)
!
subroutine en_freesolv(evec,rs)
!
  use energies
  use polypep
!
  implicit none
!
  integer i,rs,k
  RTYPE ssav,evec(MAXENERGYTERMS),fs_ipol,fos_Tdep
!
! all we need to do is sum up all solvation groups for that particular residue
! (remember that part of a solvation group might be formed by atoms on residue
! rs-1, though)
  do i=1,at(rs)%nfosgrps
    ssav = 0.0
    do k=1,at(rs)%fosgrp(i)%nats
      ssav = ssav + at(rs)%fosgrp(i)%wts(k)*&
 &                               fs_ipol(at(rs)%fosgrp(i)%ats(k))
    end do
    evec(4) = evec(4) + scale_IMPSOLV*ssav*fos_Tdep(at(rs)%fosgrp(i)%val(1:3))
  end do
!
end
!
!---------------------------------------------------------------------------
!
! "epucker" sets the value of a restraint energy to penalize excursions
! from canonical values for the phi angles; these restraints are
! set based on the pucker geometry at the gamma carbon which is realized
! by variation in the chi angle; for exo  (chi1 < 0, Cg far from Ci)
! the width of the potential attempts to reprdouce the distribution of
! phi angles seen in the PDB; similarly for endo (chi1 > 0, Cg near to Ci)
! the width of the potential attempts to reprdouce the distribution of
! phi angles seen in the PDB
!
! WARNING: epucker has to be scaled from the outside!
!
function epucker(rs)
!
  use fyoc
  use zmatrix
  use system
!
  implicit none
!
  integer rs,ff
  RTYPE delta,ee,epucker
  RTYPE kphi,phival
!
! first, for each pucker, P(phi) is obtained from the database
! this is followed by obtaining -log(P(phi)) which in turn is fit
! to a parabolic restraint potential of the form K*(phi-phio)**2
! an effective spring constant K' = RT*K = 0.006 (T=298K) is used
!
  epucker = 0.0d0
! in kcal/(mol*deg^2) as a function of temperature i believe
  kphi = 0.01/invtemp
  ff = fline(rs)
  phival = ztor(ff)
!
  if (chiral(rs).eq.1) then
    if(chi(1,rs) .lt. 0.0d0) then
      delta = phival + 60.8d0
      ee = kphi*delta**2
      epucker = epucker + ee
    else
      delta = phival + 71.0d0
      ee = kphi*delta**2
      epucker = epucker + ee
    end if
  else
    if(chi(1,rs) .gt. 0.0d0) then
      delta = phival - 60.8d0
      ee = kphi*delta**2
      epucker = epucker + ee
    else
      delta = phival - 71.0d0
      ee = kphi*delta**2
      epucker = epucker + ee
    end if
  end if
!
end
!
!--------------------------------------------------------------------------------
!
! the energy term needed to keep the peptide unit more or less planar, while
! allowing for cis/trans isomerization (unlike other torsions, this one is not
! guided by sterics, but by non-NB-capturable electronic effects
! scaled from the outside
!
function eomega(rs,resname)
!
  use fyoc
  use math
  use energies
  use sequen
  use polypep
  use system
!
  implicit none
!
  integer rs,shf
  character(3) resname
  RTYPE eomega,t1,t2,t3,t4,getztor
!
! note that we'll do this only for secondary amides, primary/tertiary ones are a bit too redundant
!
!  write(*,*) t1-t2,t3-t4,t2-t4,t1-t3
  if (resname.eq.'NMF') then
    t1=omega(rs)
    if (ua_model.gt.0) then
      t2=getztor(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(4))
      eomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) + &
 &             (2.300/2.0)*(1.0 + cos(1.0*t2/RADIAN)) + &
 &             (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN))
      return
    end if
    t2=getztor(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(5))
    t3=getztor(at(rs)%bb(4),at(rs)%bb(1),at(rs)%bb(3),at(rs)%sc(1))
    t4=getztor(at(rs)%bb(4),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(5))
!    write(*,*) t1-t2,t3-t4,t2-t4,t1-t3
  else if (resname.eq.'NMA') then
    t1=omega(rs)
    t2=getztor(at(rs)%bb(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(4))
    t3=getztor(at(rs)%sc(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%sc(1))
    t4=getztor(at(rs)%sc(2),at(rs)%bb(1),at(rs)%bb(3),at(rs)%bb(4))
  else if (resname.eq.'PRO') then
    shf = 0
    if (ua_model.gt.0) shf = 1
    t1 = getztor(oi(rs-1),ci(rs-1),ni(rs),cai(rs))
    t2 = getztor(oi(rs-1),ci(rs-1),ni(rs),at(rs)%sc(4-shf))
    if ((ua_model.gt.0).AND.(seqtyp(rs-1).eq.28)) then
      eomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) + &
 &             (2.300/2.0)*(1.0 + cos(1.0*t2/RADIAN)) + &
 &             (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN))
      return
    end if
    t3 = omega(rs)
    if (seqtyp(rs-1).eq.28) then
      t4 = getztor(at(rs-1)%bb(3),ci(rs-1),ni(rs),at(rs)%sc(4-shf))
    else
      t4 = getztor(cai(rs-1),ci(rs-1),ni(rs),at(rs)%sc(4-shf))
    end if
    eomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) + &
 &         (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t3/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t3/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t4/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t4/RADIAN))
!   proline is a tertiary amide -> different torsions
    return
  else
    t1 = getztor(oi(rs-1),ci(rs-1),ni(rs),cai(rs))
    t2 = getztor(oi(rs-1),ci(rs-1),ni(rs),hni(rs))
    if ((ua_model.gt.0).AND.(seqtyp(rs-1).eq.28)) then
      eomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) + &
 &             (2.300/2.0)*(1.0 + cos(1.0*t2/RADIAN)) + &
 &             (6.089/2.0)*(1.0 - cos(2.0*t2/RADIAN))
      return
    end if
    t3 = omega(rs)
    if (seqtyp(rs-1).eq.28) then
      t4 = getztor(at(rs-1)%bb(3),ci(rs-1),ni(rs),hni(rs))
    else
      t4 = getztor(cai(rs-1),ci(rs-1),ni(rs),hni(rs))
    end if
  end if
!
  eomega = (6.089/2.0)*(1.0 - cos(2.0*t1/RADIAN)) + &
 &         (4.900/2.0)*(1.0 - cos(2.0*t2/RADIAN)) +&
 &         (2.300/2.0)*(1.0 + cos(1.0*t3/RADIAN)) +&
 &         (6.089/2.0)*(1.0 - cos(2.0*t3/RADIAN)) +&
 &         (4.900/2.0)*(1.0 - cos(2.0*t4/RADIAN))
!
end
!
!--------------------------------------------------------------------------------
!
! the energy term needed to keep phenolic polar hydrogens more or less in the
! plane of the aromatic ring, while allowing for "cis/trans" isomerization.
! again, this is based on electronic effects and actually counterintuitive based
! on sterics alone
! scaled from the outside
!
function ephenol(rs,resname)
!
  use iounit
  use fyoc
  use math
  use system
  use polypep
!
  implicit none
!
  integer rs,shf,shf2,shf3
  RTYPE ephenol,t1,t2,getztor
  character(3) resname
!
  ephenol = 0.0
  shf = 0
  shf2 = 0
  shf3 = 0
  if (ua_model.gt.0) then
    shf = 1
    if (ua_model.eq.1) then
      shf2 = 3
      shf3 = 3
    else
      shf2 = 3
      shf3 = 7
    end if
  end if
!
  if (resname.eq.'TYR') then
    t1 = chi(3,rs)
    t2=getztor(at(rs)%sc(6-shf),at(rs)%sc(8-shf),at(rs)%sc(9-shf),at(rs)%sc(16-shf3))
    ephenol = 0.5*1.682*((1.0 - cos(2.0*t1/RADIAN)) + &
 &                   (1.0 - cos(2.0*t2/RADIAN)))
  else if (resname.eq.'PCR') then
    t1 = chi(2-shf,rs)
    t2=getztor(at(rs)%bb(6),at(rs)%bb(7),at(rs)%bb(8),at(rs)%sc(4-shf2))
    ephenol = 0.5*1.682*((1.0 - cos(2.0*t1/RADIAN)) + &
 &                   (1.0 - cos(2.0*t2/RADIAN)))
  else
    write(ilog,*) 'Fatal. Called ephenol(...) with an inapplicable r&
 &esidue (',rs,': ',resname,').'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------
!
! these are routines meant to provide energies from bonded terms
!
subroutine en_bonds(rs,evec)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use sequen
  use fyoc
!
  implicit none
!
  integer ibl1,ibl2,i,rs
  RTYPE dvec(nrsbleff(rs),3),d2(nrsbleff(rs)),d1(nrsbleff(rs))
  RTYPE term0,term1,evec(MAXENERGYTERMS),globscal,svec(3),svecref(3)
!
  if (nrsbleff(rs).eq.0) return
!
  globscal = scale_BOND(1)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(15).EQV..false.) then
        return
      else
        globscal = scale_FEGS(15)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsbleff(rs)
    ibl1 = iaa(rs)%bl(i,1)
    ibl2 = iaa(rs)%bl(i,2)
    svec(:) = 0.0
    if (atmres(ibl2).eq.disulf(rs)) svec(:) = svecref(:)
    dvec(i,1) = x(ibl2) - x(ibl1) + svec(1)
    dvec(i,2) = y(ibl2) - y(ibl1) + svec(2)
    dvec(i,3) = z(ibl2) - z(ibl1) + svec(3)
  end do
  d2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
  d1(:) = sqrt(d2(:))
!
  do i=1,nrsbleff(rs)
    ibl1 = iaa(rs)%bl(i,1)
    ibl2 = iaa(rs)%bl(i,2)
!   harmonic
    if (iaa(rs)%typ_bl(i).eq.1) then
      term0 = globscal*iaa(rs)%par_bl(i,1)
      evec(15) = evec(15) + term0*((d1(i)-iaa(rs)%par_bl(i,2))**2)
!   Morse (anharmonic)
    else if (iaa(rs)%typ_bl(i).eq.2) then
      term0 = globscal*iaa(rs)%par_bl(i,3)
      term1 = exp(-iaa(rs)%par_bl(i,1)*(d1(i)-iaa(rs)%par_bl(i,2)))
      evec(15) = evec(15) + term0*((1.0 - term1)**2)
!   GROMOS quartic
    else if (iaa(rs)%typ_bl(i).eq.3) then
      term0 = 0.25*globscal*iaa(rs)%par_bl(i,1)
      evec(15) = evec(15) + term0*((d2(i)-iaa(rs)%par_bl(i,3))**2)
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of bond length&
 & potential in en_bonds(...) (offending type is ',&
 &iaa(rs)%typ_bl(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the same only that it allows disregarding of crosslink bridges 
! 
subroutine en_bonds_alt(rs,evec,cycflag)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use sequen
  use fyoc
!
  implicit none
!
  logical, INTENT(IN):: cycflag
  integer, INTENT(IN):: rs
!
  integer ibl1,ibl2,i
  RTYPE dvec(nrsbleff(rs),3),d2(nrsbleff(rs)),d1(nrsbleff(rs))
  RTYPE term0,term1,evec(MAXENERGYTERMS),globscal,svec(3),svecref(3),rdum
!
  if (nrsbleff(rs).eq.0) return
!
  globscal = scale_BOND(1)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(15).EQV..false.) then
        return
      else
        globscal = scale_FEGS(15)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsbleff(rs)
    ibl1 = iaa(rs)%bl(i,1)
    ibl2 = iaa(rs)%bl(i,2)
    svec(:) = 0.0
    if (atmres(ibl2).eq.disulf(rs)) svec(:) = svecref(:)
    dvec(i,1) = x(ibl2) - x(ibl1) + svec(1)
    dvec(i,2) = y(ibl2) - y(ibl1) + svec(2)
    dvec(i,3) = z(ibl2) - z(ibl1) + svec(3)
  end do
  d2(:) = dvec(:,1)**2 + dvec(:,2)**2 + dvec(:,3)**2
  d1(:) = sqrt(d2(:))
!
  rdum = globscal
  do i=1,nrsbleff(rs)
    ibl1 = iaa(rs)%bl(i,1)
    ibl2 = iaa(rs)%bl(i,2)
    if (cycflag.EQV..true.) then
      if (disulf(atmres(ibl1)).eq.atmres(ibl2)) then
        rdum = 0.
      else
        rdum = globscal
      end if
    end if
!   harmonic
    if (iaa(rs)%typ_bl(i).eq.1) then
      term0 = rdum*iaa(rs)%par_bl(i,1)
      evec(15) = evec(15) + term0*((d1(i)-iaa(rs)%par_bl(i,2))**2)
!   Morse (anharmonic)
    else if (iaa(rs)%typ_bl(i).eq.2) then
      term0 = rdum*iaa(rs)%par_bl(i,3)
      term1 = exp(-iaa(rs)%par_bl(i,1)*(d1(i)-iaa(rs)%par_bl(i,2)))
      evec(15) = evec(15) + term0*((1.0 - term1)**2)
!   GROMOS quartic
    else if (iaa(rs)%typ_bl(i).eq.3) then
      term0 = 0.25*rdum*iaa(rs)%par_bl(i,1)
      evec(15) = evec(15) + term0*((d2(i)-iaa(rs)%par_bl(i,3))**2)
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of bond length&
 & potential in en_bonds_alt(...) (offending type is ',&
 &iaa(rs)%typ_bl(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! same for as en_bonds for bond angles
!
subroutine en_angles(rs,evec)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use math
  use sequen
  use fyoc
!
  implicit none
!
  integer iba1,iba2,iba3,i,rs
  RTYPE dv12(nrsbaeff(rs),3),id12(nrsbaeff(rs)),dv31(nrsbaeff(rs),3)
  RTYPE dv32(nrsbaeff(rs),3),id32(nrsbaeff(rs)),dotp(nrsbaeff(rs))
  RTYPE rang(nrsbaeff(rs)),vecs(nrsbaeff(rs),3),nvecs(nrsbaeff(rs)),globscal
  RTYPE term0,term2,evec(MAXENERGYTERMS),d2,d1,id1,si(nrsbaeff(rs))
  RTYPE svec(3),svecref(3)
!
  if (nrsbaeff(rs).eq.0) return
!
  globscal = scale_BOND(2)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(16).EQV..false.) then
        return
      else
        globscal = scale_FEGS(16)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsbaeff(rs)
    iba1 = iaa(rs)%ba(i,1)
    iba2 = iaa(rs)%ba(i,2)
    iba3 = iaa(rs)%ba(i,3)
    svec(:) = 0.0
    if (atmres(iba2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(iba1) - x(iba2) + svec(1)
    dv12(i,2) = y(iba1) - y(iba2) + svec(2)
    dv12(i,3) = z(iba1) - z(iba2) + svec(3)
    svec(:) = 0.0
    if ((atmres(iba3).eq.disulf(rs)).AND.(atmres(iba2).ne.disulf(rs))) svec(:) = svecref(:)
    dv32(i,1) = x(iba3) - x(iba2) + svec(1)
    dv32(i,2) = y(iba3) - y(iba2) + svec(2)
    dv32(i,3) = z(iba3) - z(iba2) + svec(3)
!   add the 13-vector for Urey-Bradley
    if (iaa(rs)%typ_ba(i).eq.2) then
      svec(:) = 0.0
      if (atmres(iba3).eq.disulf(rs)) svec(:) = svecref(:)
      dv31(i,1) = x(iba3) - x(iba1) + svec(1)
      dv31(i,2) = y(iba3) - y(iba1) + svec(2)
      dv31(i,3) = z(iba3) - z(iba1) + svec(3)
    end if
  end do
  id12(:) = 1.0/sqrt(dv12(:,1)**2 + dv12(:,2)**2 + dv12(:,3)**2)
  id32(:) = 1.0/sqrt(dv32(:,1)**2 + dv32(:,2)**2 + dv32(:,3)**2)
  dotp(:) = dv12(:,1)*dv32(:,1) + dv12(:,2)*dv32(:,2) + &
 &          dv12(:,3)*dv32(:,3)
  do i=1,nrsbaeff(rs)
    if (dotp(i).lt.0.0) then
      si(i) = -1.0
    else
      si(i) = 1.0
    end if
  end do
! this half-angle formula has the advantage of being safe to use in all circumstances
! as 0.5*nvecs never gets near the -1/1 limits
  vecs(:,1) = id32(:)*dv32(:,1) - si(:)*id12(:)*dv12(:,1)
  vecs(:,2) = id32(:)*dv32(:,2) - si(:)*id12(:)*dv12(:,2)
  vecs(:,3) = id32(:)*dv32(:,3) - si(:)*id12(:)*dv12(:,3)
  nvecs(:) = sqrt(vecs(:,1)**2 + vecs(:,2)**2 + vecs(:,3)**2)
  rang(:) = si(:)*2.0*RADIAN*asin(0.5*nvecs(:))
  do i=1,nrsbaeff(rs)
    if (dotp(i).lt.0.0) then
      rang(i) = rang(i) + RADIAN*PI
    end if
  end do
!
  do i=1,nrsbaeff(rs)
    iba1 = iaa(rs)%ba(i,1)
    iba2 = iaa(rs)%ba(i,2)
    iba3 = iaa(rs)%ba(i,3)
!   harmonic
    if (iaa(rs)%typ_ba(i).eq.1) then
      term0 = globscal*iaa(rs)%par_ba(i,1)
      evec(16) = evec(16) + term0*&
 &             ((rang(i)-iaa(rs)%par_ba(i,2))**2)
!   UB
    else if (iaa(rs)%typ_ba(i).eq.2) then
!     get the extra pseudo-bond term for Urey-Bradley
      d2 = dv31(i,1)**2 + dv31(i,2)**2 + dv31(i,3)**2
      d1 = sqrt(d2)
      id1 = 1.0/d1
      term0 = globscal*iaa(rs)%par_ba(i,1)
      term2 = globscal*iaa(rs)%par_ba(i,3)
      evec(16) = evec(16) + term0*&
 &             ((rang(i)-iaa(rs)%par_ba(i,2))**2) + &
 &         term2*((d1-iaa(rs)%par_ba(i,4))**2)
!   GROMOS harmonic cosine (note we're wasting the computational advantage of the cos-form here)
    else if (iaa(rs)%typ_ba(i).eq.3) then
       term0 = 0.5*globscal*iaa(rs)%par_ba(i,1)
       term2 = dotp(i)*(id12(i)*id32(i))
       evec(16) = evec(16) + term0*((term2 - iaa(rs)%par_ba(i,3))**2)
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of bond angle &
 &potential in en_angles(...) (offending type is ',&
 &iaa(rs)%typ_ba(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine en_angles_alt(rs,evec,cycflag)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use math
  use sequen
  use fyoc
!
  implicit none
!
  logical, INTENT(IN):: cycflag
  integer, INTENT(IN):: rs
!
  integer iba1,iba2,iba3,i
  RTYPE dv12(nrsbaeff(rs),3),id12(nrsbaeff(rs)),dv31(nrsbaeff(rs),3)
  RTYPE dv32(nrsbaeff(rs),3),id32(nrsbaeff(rs)),dotp(nrsbaeff(rs))
  RTYPE rang(nrsbaeff(rs)),vecs(nrsbaeff(rs),3),nvecs(nrsbaeff(rs)),globscal
  RTYPE term0,term2,evec(MAXENERGYTERMS),d2,d1,id1,si(nrsbaeff(rs))
  RTYPE svec(3),svecref(3),rdum
!
  if (nrsbaeff(rs).eq.0) return
!
  globscal = scale_BOND(2)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(16).EQV..false.) then
        return
      else
        globscal = scale_FEGS(16)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsbaeff(rs)
    iba1 = iaa(rs)%ba(i,1)
    iba2 = iaa(rs)%ba(i,2)
    iba3 = iaa(rs)%ba(i,3)
    svec(:) = 0.0
    if (atmres(iba2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(iba1) - x(iba2) + svec(1)
    dv12(i,2) = y(iba1) - y(iba2) + svec(2)
    dv12(i,3) = z(iba1) - z(iba2) + svec(3)
    svec(:) = 0.0
    if ((atmres(iba3).eq.disulf(rs)).AND.(atmres(iba2).ne.disulf(rs))) svec(:) = svecref(:)
    dv32(i,1) = x(iba3) - x(iba2) + svec(1)
    dv32(i,2) = y(iba3) - y(iba2) + svec(2)
    dv32(i,3) = z(iba3) - z(iba2) + svec(3)
!   add the 13-vector for Urey-Bradley
    if (iaa(rs)%typ_ba(i).eq.2) then
      svec(:) = 0.0
      if (atmres(iba3).eq.disulf(rs)) svec(:) = svecref(:)
      dv31(i,1) = x(iba3) - x(iba1) + svec(1)
      dv31(i,2) = y(iba3) - y(iba1) + svec(2)
      dv31(i,3) = z(iba3) - z(iba1) + svec(3)
    end if
  end do
  id12(:) = 1.0/sqrt(dv12(:,1)**2 + dv12(:,2)**2 + dv12(:,3)**2)
  id32(:) = 1.0/sqrt(dv32(:,1)**2 + dv32(:,2)**2 + dv32(:,3)**2)
  dotp(:) = dv12(:,1)*dv32(:,1) + dv12(:,2)*dv32(:,2) + &
 &          dv12(:,3)*dv32(:,3)
  do i=1,nrsbaeff(rs)
    if (dotp(i).lt.0.0) then
      si(i) = -1.0
    else
      si(i) = 1.0
    end if
  end do
! this half-angle formula has the advantage of being safe to use in all circumstances
! as 0.5*nvecs never gets near the -1/1 limits
  vecs(:,1) = id32(:)*dv32(:,1) - si(:)*id12(:)*dv12(:,1)
  vecs(:,2) = id32(:)*dv32(:,2) - si(:)*id12(:)*dv12(:,2)
  vecs(:,3) = id32(:)*dv32(:,3) - si(:)*id12(:)*dv12(:,3)
  nvecs(:) = sqrt(vecs(:,1)**2 + vecs(:,2)**2 + vecs(:,3)**2)
  rang(:) = si(:)*2.0*RADIAN*asin(0.5*nvecs(:))
  do i=1,nrsbaeff(rs)
    if (dotp(i).lt.0.0) then
      rang(i) = rang(i) + RADIAN*PI
    end if
  end do
!
  rdum = globscal
  do i=1,nrsbaeff(rs)
    iba1 = iaa(rs)%ba(i,1)
    iba2 = iaa(rs)%ba(i,2)
    iba3 = iaa(rs)%ba(i,3)
    if (cycflag.EQV..true.) then
      if ((disulf(atmres(iba1)).eq.atmres(iba2)).OR.(disulf(atmres(iba1)).eq.atmres(iba3)).OR.&
 &        (disulf(atmres(iba2)).eq.atmres(iba3))) then
        rdum = 0.
      else
        rdum = globscal
      end if
    end if
!   harmonic
    if (iaa(rs)%typ_ba(i).eq.1) then
      term0 = rdum*iaa(rs)%par_ba(i,1)
      evec(16) = evec(16) + term0*&
 &             ((rang(i)-iaa(rs)%par_ba(i,2))**2)
!   UB
    else if (iaa(rs)%typ_ba(i).eq.2) then
!     get the extra pseudo-bond term for Urey-Bradley
      d2 = dv31(i,1)**2 + dv31(i,2)**2 + dv31(i,3)**2
      d1 = sqrt(d2)
      id1 = 1.0/d1
      term0 = rdum*iaa(rs)%par_ba(i,1)
      term2 = rdum*iaa(rs)%par_ba(i,3)
      evec(16) = evec(16) + term0*&
 &             ((rang(i)-iaa(rs)%par_ba(i,2))**2) + &
 &         term2*((d1-iaa(rs)%par_ba(i,4))**2)
!   GROMOS harmonic cosine (note we're wasting the computational advantage of the cos-form here)
    else if (iaa(rs)%typ_ba(i).eq.3) then
       term0 = 0.5*rdum*iaa(rs)%par_ba(i,1)
       term2 = dotp(i)*(id12(i)*id32(i))
       evec(16) = evec(16) + term0*((term2 - iaa(rs)%par_ba(i,3))**2)
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of bond angle &
 &potential in en_angles_alt(...) (offending type is ',&
 &iaa(rs)%typ_ba(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine en_torsions(rs,evec)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use math
  use sequen
  use fyoc
  use forces
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,i,rs
  RTYPE dv12(nrsdieff(rs),3),dv23(nrsdieff(rs),3)
  RTYPE dv43(nrsdieff(rs),3),ndv23(nrsdieff(rs))
  RTYPE incp123(nrsdieff(rs)),incp234(nrsdieff(rs))
  RTYPE ncp123(nrsdieff(rs)),ncp234(nrsdieff(rs))
  RTYPE dotp(nrsdieff(rs))
  RTYPE cp123(nrsdieff(rs),3),cp234(nrsdieff(rs),3)
  RTYPE dotpah(nrsdieff(rs)),si,vv,svec(3),svecref(3)
  RTYPE termc(nrsdieff(rs)),terms(nrsdieff(rs))
  RTYPE inpm(nrsdieff(rs)),st2,ct2,cmin,smin,globscal
  RTYPE term0,ctp2,ctp3,ctp4,ctp5,ctp6,evec(MAXENERGYTERMS)
!
  if (nrsdieff(rs).eq.0) return
!
  globscal = scale_BOND(4)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(18).EQV..false.) then
        return
      else
        globscal = scale_FEGS(18)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsdieff(rs)
    ito1 = iaa(rs)%di(i,1)
    ito2 = iaa(rs)%di(i,2)
    ito3 = iaa(rs)%di(i,3)
    ito4 = iaa(rs)%di(i,4)
    svec(:) = 0.0
    if (atmres(ito2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(ito1) - x(ito2) + svec(1)
    dv12(i,2) = y(ito1) - y(ito2) + svec(2)
    dv12(i,3) = z(ito1) - z(ito2) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito3).eq.disulf(rs)).AND.(atmres(ito2).ne.disulf(rs))) svec(:) = -svecref(:)
    if ((atmres(ito3).ne.disulf(rs)).AND.(atmres(ito2).eq.disulf(rs))) svec(:) = svecref(:)
    dv23(i,1) = x(ito2) - x(ito3) + svec(1)
    dv23(i,2) = y(ito2) - y(ito3) + svec(2)
    dv23(i,3) = z(ito2) - z(ito3) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito4).eq.disulf(rs)).AND.(atmres(ito3).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito4).ne.disulf(rs)).AND.(atmres(ito3).eq.disulf(rs))) svec(:) = -svecref(:)
    dv43(i,1) = x(ito4) - x(ito3) + svec(1)
    dv43(i,2) = y(ito4) - y(ito3) + svec(2)
    dv43(i,3) = z(ito4) - z(ito3) + svec(3)
  end do
  cp123(:,1) = dv12(:,2)*dv23(:,3) - dv12(:,3)*dv23(:,2)
  cp123(:,2) = dv12(:,3)*dv23(:,1) - dv12(:,1)*dv23(:,3)
  cp123(:,3) = dv12(:,1)*dv23(:,2) - dv12(:,2)*dv23(:,1)
  cp234(:,1) = dv43(:,2)*dv23(:,3) - dv43(:,3)*dv23(:,2)
  cp234(:,2) = dv43(:,3)*dv23(:,1) - dv43(:,1)*dv23(:,3)
  cp234(:,3) = dv43(:,1)*dv23(:,2) - dv43(:,2)*dv23(:,1)
! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
  dotpah(:) = cp123(:,1)*dv43(:,1) + cp123(:,2)*dv43(:,2) +&
 &           cp123(:,3)*dv43(:,3)
  ndv23(:) = sqrt(dv23(:,1)**2 + dv23(:,2)**2 + dv23(:,3)**2)
!
! filter for exceptions which would otherwise NaN
  do i=1,nrsdieff(rs)
!   colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(8) = fo_wrncnt(8) + 1
      if (fo_wrncnt(8).eq.fo_wrnlmt(8)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f proper dihedral energies. This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(8),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(8).gt.0.5*HUGE(fo_wrnlmt(8))) then
          fo_wrncnt(8) = 0
        else
          fo_wrnlmt(8) = fo_wrnlmt(8)*10
        end if
      end if
      incp123(i) = 0.0
      incp234(i) = 0.0
    else
      incp123(i) = 1.0/ncp123(i)
      incp234(i) = 1.0/ncp234(i)
    end if
  end do
!
! now get the cosine and sine term (sine term from Blondel equ. 49)
  inpm(:) = sqrt(incp123(:)*incp234(:))
  termc(:) = dotp(:)*inpm(:)
  terms(:) = inpm(:)*ndv23(:)*dotpah(:)
!
  do i=1,nrsdieff(rs)
    if (iaa(rs)%typ_di(i).eq.0) cycle
!   Ryckaert-Bellemans
    if (iaa(rs)%typ_di(i).eq.1) then
      ctp2 = termc(i)*termc(i)
      ctp3 = ctp2*termc(i)
      ctp4 = ctp3*termc(i)
      ctp5 = ctp4*termc(i)
      ctp6 = ctp5*termc(i)
      evec(18) = evec(18) + globscal*(iaa(rs)%par_di(i,1) + &
 &                   iaa(rs)%par_di(i,2)*termc(i) +&
 &                   iaa(rs)%par_di(i,3)*ctp2 +&
 &                   iaa(rs)%par_di(i,4)*ctp3 +&
 &                   iaa(rs)%par_di(i,5)*ctp4 +&
 &                   iaa(rs)%par_di(i,6)*ctp5 +&
 &                   iaa(rs)%par_di(i,7)*ctp6)
!      write(*,*) bio_botyp(b_type(iaa(rs)%di(i,1:4)))
!      write(*,*) evec(18), globscal*(iaa(rs)%par_di(i,1) + &
! &                   iaa(rs)%par_di(i,2)*termc(i) +&
! &                   iaa(rs)%par_di(i,3)*ctp2 +&
! &                   iaa(rs)%par_di(i,4)*ctp3 +&
! &                   iaa(rs)%par_di(i,5)*ctp4 +&
! &                   iaa(rs)%par_di(i,6)*ctp5 +&
! &                   iaa(rs)%par_di(i,7)*ctp6),iaa(rs)%par_di(i,1)
!   harmonic
    else if (iaa(rs)%typ_di(i).eq.2) then

      cmin = cos(iaa(rs)%par_di(i,2)/RADIAN) 
      smin = sin(iaa(rs)%par_di(i,2)/RADIAN)
!     get cos(phi-phi_0) and sin(phi-phi_0)
      ct2 = cmin*termc(i) + smin*terms(i)
      st2 = cmin*terms(i) - smin*termc(i)
      if (st2.lt.0.0) then
        si = -1.0
      else
        si = 1.0
      end if
      if (abs(ct2).gt.0.1) then ! safe regime for asin
        vv = RADIAN*asin(st2)
        if (ct2.lt.0.0) vv = 180.0 - vv ! safe for asin but angle very far away
      else ! switch to cosine (means angle is far away from target)
        vv = si*RADIAN*acos(ct2)
      end if
      term0 = globscal*iaa(rs)%par_di(i,1)
!     make sure periodicity is accounted for (note that the shift is a linear
!     function with slope unity such that that there is no extra chain rule derivative)
      if (vv.lt.-180.0) vv = vv + 360.0
      if (vv.gt.180.0) vv = vv - 360.0
      evec(18) = evec(18) + 0.5*term0*(vv**2)
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in en_torsions(...) (offending type is ',&
 &iaa(rs)%typ_di(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine en_impropers(rs,evec)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use math
  use forces
  use sequen
  use fyoc
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,i,rs
  RTYPE dv12(nrsimpteff(rs),3),dv23(nrsimpteff(rs),3)
  RTYPE dv43(nrsimpteff(rs),3),ndv23(nrsimpteff(rs))
  RTYPE incp123(nrsimpteff(rs)),incp234(nrsimpteff(rs))
  RTYPE ncp123(nrsimpteff(rs)),ncp234(nrsimpteff(rs))
  RTYPE dotp(nrsimpteff(rs)),globscal
  RTYPE cp123(nrsimpteff(rs),3),cp234(nrsimpteff(rs),3)
  RTYPE dotpah(nrsimpteff(rs)),si,vv,svec(3),svecref(3)
  RTYPE termc(nrsimpteff(rs)),terms(nrsimpteff(rs))
  RTYPE inpm(nrsimpteff(rs)),st2,ct2,cmin,smin
  RTYPE term0,ctp2,ctp3,ctp4,ctp5,ctp6,evec(MAXENERGYTERMS)
!
  if (nrsimpteff(rs).eq.0) return
!
  globscal = scale_BOND(3)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(17).EQV..false.) then
        return
      else
        globscal = scale_FEGS(17)
      end if
    end if
  end if
!
  if (disulf(rs).gt.0) then
    if (molofrs(disulf(rs)).ne.molofrs(rs)) then
      call dis_bound_rs(rs,disulf(rs),svecref)
    else
      svecref(:) = 0.0
    end if
  else
    svecref(:) = 0.0
  end if
  do i=1,nrsimpteff(rs)
    ito1 = iaa(rs)%impt(i,improper_conv(1))
    ito2 = iaa(rs)%impt(i,2)
    ito3 = iaa(rs)%impt(i,improper_conv(2))
    ito4 = iaa(rs)%impt(i,4)
    svec(:) = 0.0
    if (atmres(ito2).eq.disulf(rs)) svec(:) = -svecref(:)
    dv12(i,1) = x(ito1) - x(ito2) + svec(1)
    dv12(i,2) = y(ito1) - y(ito2) + svec(2)
    dv12(i,3) = z(ito1) - z(ito2) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito3).eq.disulf(rs)).AND.(atmres(ito2).ne.disulf(rs))) svec(:) = -svecref(:)
    if ((atmres(ito3).ne.disulf(rs)).AND.(atmres(ito2).eq.disulf(rs))) svec(:) = svecref(:)
    dv23(i,1) = x(ito2) - x(ito3) + svec(1)
    dv23(i,2) = y(ito2) - y(ito3) + svec(2)
    dv23(i,3) = z(ito2) - z(ito3) + svec(3)
    svec(:) = 0.0
    if ((atmres(ito4).eq.disulf(rs)).AND.(atmres(ito3).ne.disulf(rs))) svec(:) = svecref(:)
    if ((atmres(ito4).ne.disulf(rs)).AND.(atmres(ito3).eq.disulf(rs))) svec(:) = -svecref(:)
    dv43(i,1) = x(ito4) - x(ito3) + svec(1)
    dv43(i,2) = y(ito4) - y(ito3) + svec(2)
    dv43(i,3) = z(ito4) - z(ito3) + svec(3)
  end do
  cp123(:,1) = dv12(:,2)*dv23(:,3) - dv12(:,3)*dv23(:,2)
  cp123(:,2) = dv12(:,3)*dv23(:,1) - dv12(:,1)*dv23(:,3)
  cp123(:,3) = dv12(:,1)*dv23(:,2) - dv12(:,2)*dv23(:,1)
  cp234(:,1) = dv43(:,2)*dv23(:,3) - dv43(:,3)*dv23(:,2)
  cp234(:,2) = dv43(:,3)*dv23(:,1) - dv43(:,1)*dv23(:,3)
  cp234(:,3) = dv43(:,1)*dv23(:,2) - dv43(:,2)*dv23(:,1)
! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
  dotpah(:) = cp123(:,1)*dv43(:,1) + cp123(:,2)*dv43(:,2) +&
 &           cp123(:,3)*dv43(:,3)
  ndv23(:) = sqrt(dv23(:,1)**2 + dv23(:,2)**2 + dv23(:,3)**2)
!
! filter for exceptions which would otherwise NaN
  do i=1,nrsimpteff(rs)
!   colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(7) = fo_wrncnt(7) + 1
      if (fo_wrncnt(7).eq.fo_wrnlmt(7)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f improper dihedral angle energies. This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(7),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(7).gt.0.5*HUGE(fo_wrnlmt(7))) then
          fo_wrncnt(7) = 0
        else
          fo_wrnlmt(7) = fo_wrnlmt(7)*10
        end if
      end if
      incp123(i) = 0.0
      incp234(i) = 0.0
    else
      incp123(i) = 1.0/ncp123(i)
      incp234(i) = 1.0/ncp234(i)
    end if
  end do
!
! now get the cosine and sine term (sine term from Blondel equ. 49)
  inpm(:) = sqrt(incp123(:)*incp234(:))
  termc(:) = dotp(:)*inpm(:)
  terms(:) = inpm(:)*ndv23(:)*dotpah(:)
!
  do i=1,nrsimpteff(rs)
    if (iaa(rs)%typ_impt(i).eq.0) cycle
!   Ryckaert-Bellemans
    if (iaa(rs)%typ_impt(i).eq.1) then
      ctp2 = termc(i)*termc(i)
      ctp3 = ctp2*termc(i)
      ctp4 = ctp3*termc(i)
      ctp5 = ctp4*termc(i)
      ctp6 = ctp5*termc(i)
      evec(17) = evec(17) + globscal*(iaa(rs)%par_impt(i,1) + &
 &                   iaa(rs)%par_impt(i,2)*termc(i) +&
 &                   iaa(rs)%par_impt(i,3)*ctp2 +&
 &                   iaa(rs)%par_impt(i,4)*ctp3 +&
 &                   iaa(rs)%par_impt(i,5)*ctp4 +&
 &                   iaa(rs)%par_impt(i,6)*ctp5 +&
 &                   iaa(rs)%par_impt(i,7)*ctp6 )
!   harmonic
    else if (iaa(rs)%typ_impt(i).eq.2) then

      cmin = cos(iaa(rs)%par_impt(i,2)/RADIAN) 
      smin = sin(iaa(rs)%par_impt(i,2)/RADIAN)
!     get cos(phi-phi_0) and sin(phi-phi_0)
      ct2 = cmin*termc(i) + smin*terms(i)
      st2 = cmin*terms(i) - smin*termc(i)
      if (st2.lt.0.0) then
        si = -1.0
      else
        si = 1.0
      end if
      if (abs(ct2).gt.0.1) then ! safe regime for asin
        vv = RADIAN*asin(st2)
        if (ct2.lt.0.0) vv = 180.0 - vv ! safe for asin but angle very far away
      else ! switch to cosine (means angle is far away from target)
        vv = si*RADIAN*acos(ct2)
      end if
      term0 = globscal*iaa(rs)%par_impt(i,1)
!     make sure periodicity is accounted for (note that the shift is a linear
!     function with slope unity such that that there is no extra chain rule derivative)
      if (vv.lt.-180.0) vv = vv + 360.0
      if (vv.gt.180.0) vv = vv - 360.0
      evec(17) = evec(17) + 0.5*term0*(vv**2)
!      write(*,*) evec(17),0.5*term0*(vv**2),iaa(rs)%par_impt(i,1)
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in en_impropers(...) (offending type is ',&
 &iaa(rs)%typ_impt(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! WARNING: This function assumes that the participating atoms are not
!          part of a crosslinked residue
!
subroutine en_CMAP(rs,evec)
!
  use iounit
  use params
  use energies
  use inter
  use atoms
  use math
  use forces
!
  implicit none
!
  integer ito1,ito2,ito3,ito4,ito5,i,rs
  RTYPE dv12(nrscmeff(rs),3),dv23(nrscmeff(rs),3),dv54(nrscmeff(rs),3)
  RTYPE dv43(nrscmeff(rs),3),ndv23(nrscmeff(rs)),ndv43(nrscmeff(rs))
  RTYPE incp123(nrscmeff(rs)),incp234(nrscmeff(rs)),incp345(nrscmeff(rs))
  RTYPE ncp123(nrscmeff(rs)),ncp234(nrscmeff(rs)),ncp345(nrscmeff(rs))
  RTYPE dotp(nrscmeff(rs)),dotp2(nrscmeff(rs))
  RTYPE cp123(nrscmeff(rs),3),cp234(nrscmeff(rs),3),cp345(nrscmeff(rs),3)
  RTYPE dotpah(nrscmeff(rs)),dotpah2(nrscmeff(rs))
  RTYPE termc(nrscmeff(rs)),terms(nrscmeff(rs))
  RTYPE termc2(nrscmeff(rs)),terms2(nrscmeff(rs))
  RTYPE inpm(nrscmeff(rs)),inpm2(nrscmeff(rs))
  RTYPE globscal,si,si2,vv,vv2,dvv,dvv2
  RTYPE term0,evec(MAXENERGYTERMS)
!
  if (nrscmeff(rs).eq.0) return
!
  globscal = scale_BOND(5)
  if ((use_FEG.EQV..true.).AND.(fegmode.eq.2)) then
    if (par_FEG(rs).EQV..true.) then
      if (use_FEGS(20).EQV..false.) then
        return
      else
        globscal = scale_FEGS(20)
      end if
    end if
  end if
!
  do i=1,nrscmeff(rs)
    ito1 = iaa(rs)%cm(i,1)
    ito2 = iaa(rs)%cm(i,2)
    ito3 = iaa(rs)%cm(i,3)
    ito4 = iaa(rs)%cm(i,4)
    ito5 = iaa(rs)%cm(i,5)
    dv12(i,1) = x(ito1) - x(ito2)
    dv12(i,2) = y(ito1) - y(ito2)
    dv12(i,3) = z(ito1) - z(ito2)
    dv23(i,1) = x(ito2) - x(ito3)
    dv23(i,2) = y(ito2) - y(ito3)
    dv23(i,3) = z(ito2) - z(ito3)
    dv43(i,1) = x(ito4) - x(ito3)
    dv43(i,2) = y(ito4) - y(ito3)
    dv43(i,3) = z(ito4) - z(ito3)
    dv54(i,1) = x(ito5) - x(ito4)
    dv54(i,2) = y(ito5) - y(ito4)
    dv54(i,3) = z(ito5) - z(ito4)
  end do
  cp123(:,1) = dv12(:,2)*dv23(:,3) - dv12(:,3)*dv23(:,2)
  cp123(:,2) = dv12(:,3)*dv23(:,1) - dv12(:,1)*dv23(:,3)
  cp123(:,3) = dv12(:,1)*dv23(:,2) - dv12(:,2)*dv23(:,1)
  cp234(:,1) = dv43(:,2)*dv23(:,3) - dv43(:,3)*dv23(:,2)
  cp234(:,2) = dv43(:,3)*dv23(:,1) - dv43(:,1)*dv23(:,3)
  cp234(:,3) = dv43(:,1)*dv23(:,2) - dv43(:,2)*dv23(:,1)
  cp345(:,1) = -dv54(:,2)*dv43(:,3) + dv54(:,3)*dv43(:,2)
  cp345(:,2) = -dv54(:,3)*dv43(:,1) + dv54(:,1)*dv43(:,3)
  cp345(:,3) = -dv54(:,1)*dv43(:,2) + dv54(:,2)*dv43(:,1)

! norms and dot-product
  ncp123(:) = cp123(:,1)**2 + cp123(:,2)**2 + cp123(:,3)**2
  ncp234(:) = cp234(:,1)**2 + cp234(:,2)**2 + cp234(:,3)**2
  ncp345(:) = cp345(:,1)**2 + cp345(:,2)**2 + cp345(:,3)**2

  dotp(:) = cp123(:,1)*cp234(:,1) + cp123(:,2)*cp234(:,2) + &
 &          cp123(:,3)*cp234(:,3)
  dotp2(:) = cp234(:,1)*cp345(:,1) + cp234(:,2)*cp345(:,2) + &
 &          cp234(:,3)*cp345(:,3)
  dotpah(:) = cp123(:,1)*dv43(:,1) + cp123(:,2)*dv43(:,2) +&
 &           cp123(:,3)*dv43(:,3)
  dotpah2(:) = cp234(:,1)*dv54(:,1) + cp234(:,2)*dv54(:,2) +&
 &           cp234(:,3)*dv54(:,3)

  ndv23(:) = sqrt(dv23(:,1)**2 + dv23(:,2)**2 + dv23(:,3)**2)
  ndv43(:) = sqrt(dv43(:,1)**2 + dv43(:,2)**2 + dv43(:,3)**2)

!
! filter for exceptions which would otherwise NaN
  do i=1,nrscmeff(rs)
!   colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
    if ((ncp123(i).le.0.0).OR.(ncp234(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(9) = fo_wrncnt(9) + 1
      if (fo_wrncnt(9).eq.fo_wrnlmt(9)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f CMAP energies. This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(9),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(9).gt.0.5*HUGE(fo_wrnlmt(9))) then
          fo_wrncnt(9) = 0
        else
          fo_wrnlmt(9) = fo_wrnlmt(9)*10
        end if
      end if
      incp123(i) = 0.0
      incp234(i) = 0.0
    else
      incp123(i) = 1.0/ncp123(i)
      incp234(i) = 1.0/ncp234(i)
    end if
    if ((ncp234(i).le.0.0).OR.(ncp345(i).le.0.0)) then
!     this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
      fo_wrncnt(9) = fo_wrncnt(9) + 1
      if (fo_wrncnt(9).eq.fo_wrnlmt(9)) then
        write(ilog,*) 'WARNING. Colinear reference atoms in computation o&
 &f CMAP energies. This indicates bad parameters, an unstable simulation, or a bug&
 &. Expect run to crash soon.'
        write(ilog,*) 'This was warning #',fo_wrncnt(9),' of this type not all of which may be displayed.'
        if (10.0*fo_wrnlmt(9).gt.0.5*HUGE(fo_wrnlmt(9))) then
          fo_wrncnt(9) = 0
        else
          fo_wrnlmt(9) = fo_wrnlmt(9)*10
        end if
      end if
      incp345(i) = 0.0
      incp234(i) = 0.0
    else
      incp345(i) = 1.0/ncp345(i)
      incp234(i) = 1.0/ncp234(i)
    end if

  end do
!
! now get the cosine and sine term (sine term from Blondel equ. 49)
  inpm(:) = sqrt(incp123(:)*incp234(:))
  inpm2(:) = sqrt(incp234(:)*incp345(:))
  termc(:) = dotp(:)*inpm(:)
  termc2(:) = dotp2(:)*inpm2(:)
  terms(:) = inpm(:)*ndv23(:)*dotpah(:)
  terms2(:) = inpm2(:)*ndv43(:)*dotpah2(:)
!
  do i=1,nrscmeff(rs)
    if (iaa(rs)%typ_cm(i).eq.0) cycle
!
    if (terms(i).lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    if (abs(termc(i)).gt.0.1) then ! safe regime for asin
      vv = RADIAN*asin(terms(i))
      if (termc(i).lt.0.0) vv = 180.0 - vv
    else 
      vv = si*RADIAN*acos(termc(i))
    end if
    if (terms2(i).lt.0.0) then
      si2 = -1.0
    else
      si2 = 1.0
    end if
    if (abs(termc2(i)).gt.0.1) then ! safe regime for asin
      vv2 = RADIAN*asin(terms2(i))
      if (termc2(i).lt.0.0) vv2 = 180.0 - vv2
    else ! switch to cosine (means angle is far away from target)
      vv2 = si2*RADIAN*acos(termc2(i))
    end if
    if (vv.lt.-180.0) vv = vv + 360.0
    if (vv.gt.180.0) vv = vv - 360.0
    if (vv2.lt.-180.0) vv2 = vv2 + 360.0
    if (vv2.gt.180.0) vv2 = vv2 - 360.0

!   vv is the first angle (in order of spec.), vv2 the second

    if ((cm_typ(iaa(rs)%typ_cm(i),1).eq.1).OR.(cm_typ(iaa(rs)%typ_cm(i),1).eq.2)) then
      call cardinal_bspline2D(iaa(rs)%typ_cm(i),vv,vv2,term0,dvv,dvv2)
      evec(20) = evec(20) + globscal*term0
    else if ((cm_typ(iaa(rs)%typ_cm(i),1).eq.3).OR.(cm_typ(iaa(rs)%typ_cm(i),1).eq.4)) then
      call bicubic_spline(iaa(rs)%typ_cm(i),vv,vv2,term0,dvv,dvv2)
      evec(20) = evec(20) + globscal*term0
    else
      write(ilog,*) 'Fatal. Encountered unsupported type of torsional p&
 &otential in en_torsions(...) (offending type is ',&
 &iaa(rs)%typ_di(i),').'
      call fexit()
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine cardinal_bspline2D(which,vv1,vv2,val,dvv1,dvv2)
!
  use params
!
  implicit none
!
  integer j,k,xi,yi,sh,flr(2)
  RTYPE vv1,vv2,dvv1,dvv2,ori(2),frac(2),val
  RTYPE dum1,cm_bspl(cm_splor,2),cm_bspld(cm_splor,2),temj,temk,temjd
  integer which,basel(2)
!
  sh = cm_splor - 1

  if (cm_typ(which,1).eq.2) then
    ori(1) = -180.0 - cm_par2(which)
    ori(2) = -180.0 - cm_par2(which)
  else
    ori(1) = -180.0 - 0.5*cm_par2(which)
    ori(2) = -180.0 - 0.5*cm_par2(which)
  end if
!
  dum1 = (vv1 - ori(1))/cm_par2(which)
  flr(1) = floor(dum1)
  frac(1) = dum1 - flr(1)
  dum1 = (vv2 - ori(2))/cm_par2(which)
  flr(2) = floor(dum1)
  frac(2) = dum1 - flr(2)
  call cardBspline(frac(1),cm_splor,cm_bspl(:,1),cm_bspld(:,1))
  call cardBspline(frac(2),cm_splor,cm_bspl(:,2),cm_bspld(:,2))
  cm_bspld(1:cm_splor,1) = cm_bspld(1:cm_splor,1)/cm_par2(which)
  cm_bspld(1:cm_splor,2) = cm_bspld(1:cm_splor,2)/cm_par2(which)
  basel(1) = flr(1)
  basel(2) = flr(2)
  if (basel(1).le.0) basel(1) = basel(1) + cm_typ(which,2)
  if (basel(2).le.0) basel(2) = basel(2) + cm_typ(which,2)
  val = 0.0
  dvv1 = 0.0
  dvv2 = 0.0
  do j=1,cm_splor
    temj = 1.0*cm_bspl(j,1)
    temjd = 1.0*cm_bspld(j,1)
    xi = basel(1) + j - sh
    if (xi.gt.cm_typ(which,2)) xi = xi - cm_typ(which,2)
    if (xi.lt.1) xi = xi + cm_typ(which,2)
    do k=1,cm_splor
      yi = basel(2) + k - sh
      if (yi.gt.cm_typ(which,2)) yi = yi - cm_typ(which,2)
      if (yi.lt.1) yi = yi + cm_typ(which,2)
      temk = temj*cm_bspl(k,2)*cm_par(which,xi,yi,1)
      val = val + temk
      dvv1 = dvv1 + temjd*cm_bspl(k,2)*cm_par(which,xi,yi,1)
      dvv2 = dvv2 + temj*cm_bspld(k,2)*cm_par(which,xi,yi,1)
!      write(*,*) vv1,vv2,val,dvv1,dvv2
    end do
  end do
!
 45 format(400(f12.5))
end
!
!----------------------------------------------------------------------------------
!
! this routine returns the bicubic spline-interpolated value of the CMAP grid
! in val and the (trivial) derivatives dECMAP/dphi(psi) in dvald1 and dvald2
!
subroutine bicubic_spline(which,vv1,vv2,val,dvaldv1,dvaldv2)
!
  use params
  use mcsums
!
  implicit none
!
  integer which,idx(2,2),i,j!,ll,mm
  RTYPE vv1,vv2,dvaldv1,dvaldv2,frac(2),ori(2),coeff
  RTYPE dum1,fracp(4,2),val!,vals(3610)
!
!  ki = ceiling((vv-leftlim)/dx)
  if (cm_typ(which,1).eq.4) then
    ori(1) = -180.0 - cm_par2(which)
    ori(2) = -180.0 - cm_par2(which)
  else
    ori(1) = -180.0 - 0.5*cm_par2(which)
    ori(2) = -180.0 - 0.5*cm_par2(which)
  end if
!
!  do ll=1,720
!  do mm=1,720
!  vv1 = mm*0.5-180.5
!  vv2 = ll*0.5-180.5
!
  dum1 = (vv1 - ori(1))/cm_par2(which)
  idx(1,1) = floor(dum1)
  frac(1) = dum1 - idx(1,1)
  dum1 = (vv2 - ori(2))/cm_par2(which)
  idx(2,1) = floor(dum1)
  frac(2) = dum1 - idx(2,1)
  do while (idx(1,1).gt.cm_typ(which,2))
    idx(1,1) = idx(1,1) - cm_typ(which,2)
  end do
  do while (idx(2,1).gt.cm_typ(which,2))
    idx(2,1) = idx(2,1) - cm_typ(which,2)
  end do
  do while (idx(1,1).lt.1)
    idx(1,1) = idx(1,1) + cm_typ(which,2)
  end do
  do while (idx(2,1).lt.1)
    idx(2,1) = idx(2,1) + cm_typ(which,2)
  end do
!
  fracp(1,:) = 1.0
  do i=2,4
    fracp(i,1) = frac(1)**(i-1)
    fracp(i,2) = frac(2)**(i-1)
  end do
  val = 0.0
  dvaldv1 = 0.0
  dvaldv2 = 0.0
  do i=1,4
    do j=1,4
      coeff = cm_par(which,idx(1,1),idx(2,1),4*(j-1)+i+1)
      val = val + fracp(i,1)*fracp(j,2)*coeff
      if (i.gt.1) then
        dvaldv1 = dvaldv1 + (i-1)*fracp(i-1,1)*fracp(j,2)*coeff/cm_par2(which)
      end if
      if (j.gt.1) then
        dvaldv2 = dvaldv2 + (j-1)*fracp(j-1,2)*fracp(i,1)*coeff/cm_par2(which)
      end if
    end do
  end do
!
!  vals(mm) = val !dvaldv2
!  end do
!  write(iene,56) (vals(mm),mm=1,720)
!  end do
! 56 format(3601(f12.5,1x))
!
end
!
!---------------------------------------------------------------------------------------
!
! wrapper around e_boundary_rs
! 
subroutine e_boundary(evec)
!
  use iounit
  use sequen
  use energies
!
  implicit none
!
  integer rs,aone
  RTYPE evec(MAXENERGYTERMS)
!
  aone = 1
!
  do rs=1,nseq
    call e_boundary_rs(rs,evec,aone)
  end do
!
end
!
!--------------------------------------------------------------------------------
!
! note that we could use e_boundary_rs for hard walls as well (just pick a large threshold
! value), but this is a) not 100% clean and b) slower
!
subroutine e_boundary_rs(rs,evec,mode)
!
  use iounit
  use atoms
  use system
  use aminos
  use sequen
  use polypep
  use energies
  use cutoffs, ONLY: rsinfo
!
  implicit none
!
  integer rs,mode,i,ii,j
  RTYPE dis2,dis,evec(MAXENERGYTERMS),dvec(3),hvec(3),hvec2(3),svec(3)
!
  if ((bnd_type.eq.1).OR.(bnd_type.eq.2)) then
    if (bnd_shape.eq.3) then ! periodic cylinder w/ axis along z (continuous tube) using ASWBC
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dvec(1) = (x(ii)-bnd_params(1))
        dvec(2) = (y(ii)-bnd_params(2))
        dis2 = dvec(1)*dvec(1) + dvec(2)*dvec(2)
        dis = sqrt(dis2) - bnd_params(4)
        if (dis.gt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dis*dis
        end if
      end do
    else
!     do nothing
    end if
  else if (bnd_type.eq.3) then
!   we'll use the reference atom for the particular residue to determine the boundary violation
!
    ii = refat(rs)
    if (bnd_shape.eq.1) then
      dvec(1) = (x(ii)-bnd_params(4))
      dvec(2) = (y(ii)-bnd_params(5))
      dvec(3) = (z(ii)-bnd_params(6))
      do j=1,3
        if (dvec(j).lt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dvec(j)*dvec(j)
        else if (dvec(j).gt.bnd_params(j)) then
          evec(12) = evec(12) + bnd_params(7)*(dvec(j)-bnd_params(j))**2
        end if
      end do
    else if (bnd_shape.eq.2) then
      dis2 = (x(ii)-bnd_params(1))*(x(ii)-bnd_params(1)) + &
 &           (y(ii)-bnd_params(2))*(y(ii)-bnd_params(2)) + &
 &           (z(ii)-bnd_params(3))*(z(ii)-bnd_params(3))
      dis = sqrt(dis2) - bnd_params(4)
      if (dis.gt.0.0) then
        evec(12) = evec(12) + bnd_params(7)*dis*dis
      end if
    else if (bnd_shape.eq.3) then
      dis2 = (x(ii)-bnd_params(1))*(x(ii)-bnd_params(1)) + &
 &           (y(ii)-bnd_params(2))*(y(ii)-bnd_params(2))
      dis = sqrt(dis2) - bnd_params(4)
      if (dis.gt.0.0) then
        evec(12) = evec(12) + bnd_params(7)*dis*dis
      end if
      dvec(3) = z(ii)-bnd_params(3)+0.5*bnd_params(6)
      if (dvec(3).lt.0.0) then
        evec(12) = evec(12) + bnd_params(7)*dvec(3)*dvec(3)
      else if (dvec(3).gt.bnd_params(6)) then
        evec(12) = evec(12) + bnd_params(7)*(dvec(3)-bnd_params(6))**2
      end if  
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in e_b&
 &oundary_rs() (code # ',bnd_shape,').'
      call fexit()
    end if
  else if (bnd_type.eq.4) then
    if (bnd_shape.eq.1) then
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dvec(1) = (x(ii)-bnd_params(4))
        dvec(2) = (y(ii)-bnd_params(5))
        dvec(3) = (z(ii)-bnd_params(6))
        do j=1,3
          if (dvec(j).lt.0.0) then
            evec(12) = evec(12) + bnd_params(7)*dvec(j)*dvec(j)
          else if (dvec(j).gt.bnd_params(j)) then
            evec(12) = evec(12) + bnd_params(7)*(dvec(j)-bnd_params(j))**2
          end if
        end do
      end do
    else if (bnd_shape.eq.2) then
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dis2 = 0.0
        dis2 = (x(ii)-bnd_params(1))*(x(ii)-bnd_params(1)) + &
 &             (y(ii)-bnd_params(2))*(y(ii)-bnd_params(2)) + &
 &             (z(ii)-bnd_params(3))*(z(ii)-bnd_params(3))
        dis = sqrt(dis2) - bnd_params(4)
        if (dis.gt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dis*dis
        end if
      end do
    else if (bnd_shape.eq.3) then
      do i=1,at(rs)%nbb+at(rs)%nsc
        if (i.le.at(rs)%nbb) then
          ii = at(rs)%bb(i)
        else
          ii = at(rs)%sc(i-at(rs)%nbb)
        end if
        dis2 = (x(ii)-bnd_params(1))*(x(ii)-bnd_params(1)) + &
 &             (y(ii)-bnd_params(2))*(y(ii)-bnd_params(2))
        dis = sqrt(dis2) - bnd_params(4)
        if (dis.gt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dis*dis
        end if
        dvec(3) = z(ii)-bnd_params(3)+0.5*bnd_params(6)
        if (dvec(3).lt.0.0) then
          evec(12) = evec(12) + bnd_params(7)*dvec(3)*dvec(3)
        else if (dvec(3).gt.bnd_params(6)) then
          evec(12) = evec(12) + bnd_params(7)*(dvec(3)-bnd_params(6))**2
        end if         
      end do
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in e_b&
 &oundary_rs() (code # ',bnd_shape,').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & in e_boundary_rs() (code # ',bnd_type,').'
    call fexit()
  end if
!
  if (use_OSMO.EQV..true.) then
    if (par_OSMO2(rs).gt.0) then
!     for a 3D periodic boundary box, all boundaries are new
!     note that any (planar) boundary orthogonal to a periodic dimension becomes infinite itself
!     due to image ambiguities no boundaries should be added that would effectively be just an outer boundary
      if ((bnd_shape.eq.1).AND.(bnd_type.eq.1)) then
        hvec(1:3) = bnd_params(1:3)
        hvec2(1:3) = bnd_params(4:6)
        if (par_OSMO(1).eq.1) then
          hvec(3) = hvec(3)/2.
          if (par_OSMO2(rs).eq.1) hvec2(3) = hvec2(3) + hvec(3)
        else if (par_OSMO(1).eq.2) then
          hvec(2:3) = hvec(2:3)/2.
          if ((par_OSMO2(rs).eq.1).OR.(par_OSMO2(rs).eq.3)) hvec2(3) = hvec2(3) + hvec(3)
          if (par_OSMO2(rs).le.2) hvec2(2) = hvec2(2) + hvec(2)
        else if (par_OSMO(1).eq.3) then
          hvec(1:3) = hvec(1:3)/2.
          if (mod(par_OSMO2(rs),2).eq.1) hvec2(3) = hvec2(3) + hvec(3)
          if (mod(par_OSMO2(rs)-1,4).le.1) hvec2(2) = hvec2(2) + hvec(2)
          if (par_OSMO2(rs).le.4) hvec2(1) = hvec2(1) + hvec(1)
        end if
        if (par_OSMO(1).eq.3) then
          svec(1) = bnd_params(1)*anint((x(refat(rs))-hvec2(1)-0.5*hvec(1))/bnd_params(1))
          do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
            dvec(1) = (x(ii)-hvec2(1)) - svec(1)
            if (dvec(1).lt.0.0) then
              evec(11) = evec(11) + bnd_params(7)*dvec(1)*dvec(1)
            else if (dvec(1).gt.hvec(1)) then
              evec(11) = evec(11) + bnd_params(7)*(dvec(1)-hvec(1))**2
            end if
          end do
        end if
        if ((par_OSMO(1).eq.2).OR.(par_OSMO(1).eq.3)) then
          svec(2) = bnd_params(2)*anint((y(refat(rs))-hvec2(2)-0.5*hvec(2))/bnd_params(2))
          do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
            dvec(2) = (y(ii)-hvec2(2)) - svec(2)
            if (dvec(2).lt.0.0) then
              evec(11) = evec(11) + bnd_params(7)*dvec(2)*dvec(2)
            else if (dvec(2).gt.hvec(2)) then
              evec(11) = evec(11) + bnd_params(7)*(dvec(2)-hvec(2))**2
            end if
          end do
        end if
        if ((par_OSMO(1).ge.1).AND.(par_OSMO(1).le.3)) then
          svec(3) = bnd_params(3)*anint((z(refat(rs))-hvec2(3)-0.5*hvec(3))/bnd_params(3))
          do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
            dvec(3) = (z(ii)-hvec2(3)) - svec(3)
            if (dvec(3).lt.0.0) then
              evec(11) = evec(11) + bnd_params(7)*dvec(3)*dvec(3)
            else if (dvec(3).gt.hvec(3)) then
              evec(11) = evec(11) + bnd_params(7)*(dvec(3)-hvec(3))**2
            end if
          end do
        end if
!       for any aperiodic boundary, treatment is the same (only the inner boundaries are added and need to be considered)
      else if (((bnd_type.ne.1).AND.(bnd_shape.ge.1).AND.(bnd_shape.le.3)).OR.((bnd_type.eq.1).AND.(bnd_shape.eq.3))) then
        if ((bnd_shape.eq.2).OR.(bnd_shape.eq.3)) then
          hvec(1:3) = bnd_params(1:3)
        else if (bnd_shape.eq.1) then
          hvec(1:3) = 0.5*bnd_params(1:3) + bnd_params(4:6)
        end if
!       xy-plane
        if ((mod(par_OSMO2(rs),2).eq.0).AND.(bnd_type.ne.1)) then
          do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
            dvec(3) = z(ii)-hvec(3)
            if (dvec(3).gt.0.0) then
              evec(11) = evec(11) + bnd_params(7)*dvec(3)*dvec(3)
            end if
          end do
        else if (bnd_type.ne.1) then
          do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
            dvec(3) = z(ii)-hvec(3)
            if (dvec(3).lt.0.0) then
              evec(11) = evec(11) + bnd_params(7)*dvec(3)*dvec(3)
            end if
          end do
        end if
!       xz-plane
        if (par_OSMO(1).ge.2) then
          if (mod(par_OSMO2(rs)-1,4).ge.2) then
            do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
              dvec(2) = y(ii)-hvec(2)
              if (dvec(2).gt.0.0) then
                evec(11) = evec(11) + bnd_params(7)*dvec(2)*dvec(2)
              end if
            end do
          else
            do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
              dvec(2) = y(ii)-hvec(2)
              if (dvec(2).lt.0.0) then
                evec(11) = evec(11) + bnd_params(7)*dvec(2)*dvec(2)
              end if
            end do
          end if
        end if
!       yz-plane
        if (par_OSMO(1).ge.3) then
          if (par_OSMO2(rs).ge.5) then
            do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
              dvec(1) = x(ii)-hvec(1)
              if (dvec(1).gt.0.0) then
                evec(11) = evec(11) + bnd_params(7)*dvec(1)*dvec(1)
              end if
            end do
          else
            do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
              dvec(1) = x(ii)-hvec(1)
              if (dvec(1).lt.0.0) then
                evec(11) = evec(11) + bnd_params(7)*dvec(1)*dvec(1)
              end if
            end do
          end if
        end if
!       add cylinder caps for periodic cylinder
        if ((bnd_type.eq.1).AND.(bnd_shape.eq.3)) then
          hvec(3) = bnd_params(6)/2.
          if (mod(par_OSMO2(rs),2).eq.0) then
            hvec2(3) = bnd_params(3) - hvec(3)
          else
            hvec2(3) = bnd_params(3)
          end if
          svec(3) = bnd_params(6)*anint((z(refat(rs))-hvec2(3)-0.5*hvec(3))/bnd_params(6))
          do ii=rsinfo(rs,1),rsinfo(rs,1)+at(rs)%na-1
            dvec(3) = (z(ii)-hvec2(3)) - svec(3)
            if (dvec(3).lt.0.0) then
              evec(11) = evec(11) + bnd_params(7)*dvec(3)*dvec(3)
            else if (dvec(3).gt.hvec(3)) then
              evec(11) = evec(11) + bnd_params(7)*(dvec(3)-hvec(3))**2
            end if
          end do
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported boundary condition / container in e_boundary_rs() (code # ',bnd_shape,').'
        call fexit()
      end if
    end if
  end if
!
end
!
!-----------------------------------------------------------------------------
!
! atom-atom distance restraints
!
subroutine edrest(evec,tpi)
!
  use energies
  use distrest
  use atoms
  use sequen
#ifdef ENABLE_THREADS
  use threads, ONLY: thr_limits
#endif
  use system, ONLY: bnd_shape,bnd_type,bnd_params
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  RTYPE edr,dis,evec(MAXENERGYTERMS),dvec(3),svec(3),dis2
  integer i,ii,kk,sta,sto
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    sta = thr_limits(47,tpi)
    sto = thr_limits(48,tpi)
  else
    sta = 1
    sto = ndrest
  end if
#else
  sta = 1
  sto = ndrest
#endif
!
  edr = 0.0
!
  do i=sta,sto
!
    ii = dresat(i,1)
    kk = dresat(i,2)
    if (kk.gt.0) then
      if (molofrs(atmres(ii)).ne.molofrs(atmres(kk))) then
        call dis_bound_rs(atmres(ii),atmres(kk),svec)
        call dis_bound5(ii,kk,svec,dis2,dvec)
        dis = sqrt(dis2)
      else
        dvec(1) = x(kk)-x(ii)
        dvec(2) = y(kk)-y(ii)
        dvec(3) = z(kk)-z(ii)
        dis = sqrt(dvec(1)*dvec(1)+dvec(2)*dvec(2)+dvec(3)*dvec(3))
      end if
    else if (kk.eq.-1) then
      if ((bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        dis = x(ii) - bnd_params(1)*anint((x(ii)-dresprm(i,1))/bnd_params(1))
      else
        dis = x(ii)
      end if
    else if (kk.eq.-2) then
      if ((bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        dis = y(ii) - bnd_params(2)*anint((y(ii)-dresprm(i,1))/bnd_params(2))
      else
        dis = y(ii)
      end if
    else ! -3
      if (bnd_type.eq.1) then
        if (bnd_shape.eq.1) then
          dis = z(ii) - bnd_params(3)*anint((z(ii)-dresprm(i,1))/bnd_params(3))
        else if (bnd_shape.eq.3) then
          dis = z(ii) - bnd_params(6)*anint((z(ii)-dresprm(i,1))/bnd_params(6))
        end if
      else
        dis = z(ii)
      end if
    end if
    if (drestyp(i).eq.1) then
      edr = edr + dresprm(i,2)*(dis-dresprm(i,1))*(dis-dresprm(i,1))
    else if (drestyp(i).eq.2) then
      if (dis.gt.dresprm(i,1)) edr = edr + dresprm(i,2)*(dis-dresprm(i,1))*(dis-dresprm(i,1))
    else if (drestyp(i).eq.3) then
      if (dis.lt.dresprm(i,1)) edr = edr + dresprm(i,2)*(dis-dresprm(i,1))*(dis-dresprm(i,1))
    end if
!
  end do
!
  evec(10) = evec(10) + scale_DREST*edr
!
end
!
!------------------------------------------------------------------------------
!
! a function to compute specific residue-based torsional (biasing) potentials
!
subroutine en_torrs(evec,rs)
!
  use iounit
  use energies
  use fyoc
  use aminos
  use torsn
  use sequen
  use polypep
  use interfaces
  use system, ONLY: ua_model
  use zmatrix, ONLY: izrot,ztor
  use movesets, ONLY: natlst
!
  implicit none
!
  RTYPE evec(MAXENERGYTERMS),vv(12),getpuckertor,torincr,getztor_inter
  integer rs,j,k,zl,shf
!
  if (par_TOR2(rs).eq.0) then
    write(ilog,*) 'Fatal. Called etorrs(rs) with residue, for which &
 &no torsional potential was set up.'
    call fexit()
  end if
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
  vv(:) = 0.0
!
  if (seqpolty(rs).eq.'N') then
    do k=1,nnucs(rs)
      if (par_TOR(rs,MAXTORRES+k).gt.0) vv(k) = nucs(k,rs) - par_TOR(rs,k)
    end do
    vv(nnucs(rs)+1) = getpuckertor(rs,zl) - par_TOR(rs,nnucs(rs)+1)
  else if ((fline(rs).gt.0).OR.(yline(rs).gt.0).OR.(wline(rs).gt.0)) then
    if (par_TOR(rs,MAXTORRES+1).gt.0.0) vv(1) = omega(rs) - par_TOR(rs,1)
    if (par_TOR(rs,MAXTORRES+2).gt.0.0) vv(2) = phi(rs) - par_TOR(rs,2) 
    if (par_TOR(rs,MAXTORRES+3).gt.0.0) vv(3) = psi(rs) - par_TOR(rs,3)
  end if
  do k=1,nchi(rs)
    if (par_TOR(rs,MAXTORRES+6+k).gt.0) vv(6+k) = chi(k,rs) - par_TOR(rs,6+k)
  end do
  if (disulf(rs).gt.0) then
    if (crosslink(crlk_idx(rs))%itstype.le.2) then
      vv(6+nchi(rs)+1) = getztor_inter(cai(rs),at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf))-par_TOR(rs,6+nchi(rs)+1)
      if (disulf(rs).gt.rs) then
        vv(6+nchi(rs)+2) = getztor_inter(at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf),at(disulf(rs))%sc(2-shf)) - &
 &                         par_TOR(rs,6+nchi(rs)+2)
      end if
    end if
  end if
  if (seqtyp(rs).eq.26) then
    if ((seqpolty(rs).ne.'P').AND.(seqpolty(rs).ne.'N')) then
      j = 0
    else
      j = 6
    end if
    do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%na-1
      if (izrot(k)%alsz.gt.0) then
        if (natlst%nr.gt.0) then
          call binary_search(natlst%nr,natlst%idx(1:natlst%nr),k,zl)
          if (natlst%idx(max(1,min(natlst%nr,zl))).eq.k) cycle
        end if
        if ((k.eq.fline(rs)).OR.(k.eq.yline(rs)).OR.(k.eq.fline2(rs)).OR.(k.eq.yline2(rs))) cycle ! these are so annoying
        j = j + 1
        if (j.gt.12) exit
        vv(j) = ztor(k) - par_TOR(rs,j)
      end if
    end do
  end if 
!
! periodic boundary conditions (so to speak)
  do j=1,MAXTORRES
    if (vv(j).lt.-180.0) vv(j) = vv(j) + 360.0
    if (vv(j).gt.180.0) vv(j) = vv(j) - 360.0
  end do
  torincr = 0.0
  do j=1,MAXTORRES
    if (par_TOR(rs,MAXTORRES+j).gt.0.0) torincr = torincr + par_TOR(rs,MAXTORRES+j)*vv(j)*vv(j)
  end do
! the potential is harmonic  
  if ((par_TOR2(rs).eq.1).OR.(par_TOR2(rs).eq.3)) then
    evec(7) = evec(7) + scale_TOR*torincr
! the potential is a Gaussian well
  else if ((par_TOR2(rs).eq.2).OR.(par_TOR2(rs).eq.4)) then
    evec(7) = evec(7) - scale_TOR*exp(-torincr)
  else
    write(ilog,*) 'Fatal. Called etorrs(rs) with unknown mode. Offen&
 &ding value is ',par_TOR2(rs),' (residue ',rs,').'
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------------
!
