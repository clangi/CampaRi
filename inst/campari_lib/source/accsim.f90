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
! 
subroutine els_init()
!
  use wl
  use atoms
  use iounit
  use energies
  use system
  use units
  use cutoffs
  use params
  use threads, ONLY: thrdat
!
  implicit none
!
  integer i,k
  integer curmaxets
  RTYPE dum
  RTYPE, ALLOCATABLE:: tmpa(:,:)
!
! always allocate
  allocate(hmjam%isin(MAXENERGYTERMS))
  hmjam%isin(:) = .false.
!
  curmaxets = 23
  if (hmjam%nlst.le.0) then
    do_accelsim = .false.
    hmjam%prtfrmwts = HUGE(hmjam%prtfrmwts)
    if (allocated(hmjam%lst).EQV..true.) deallocate(hmjam%lst)
    return
  end if
  allocate(tmpa(hmjam%nlst,4))
  tmpa(:,1:2) = hmjam%thresh(1:hmjam%nlst,1:2)
  tmpa(:,3:4) = hmjam%alpha(1:hmjam%nlst,1:2)
!
  do i=1,hmjam%nlst
    if ((hmjam%lst(i).le.1).OR.(hmjam%lst(i).gt.curmaxets).OR.(hmjam%lst(i).eq.13)) then
      write(ilog,*) 'Fatal. Specified illegal dimension for energy landscape sculpting.'
      call fexit()
    end if
    if ((hmjam%lst(i).eq.3).AND.(use_IPP.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.5).AND.(use_attLJ.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.4).AND.(use_CORR.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.6).AND.(use_IMPSOLV.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.7).AND.(use_WCA.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.8).AND.(use_POLAR.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.9).AND.(use_TOR.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.10).AND.(use_ZSEC.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.11).AND.(use_TABUL.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.12).AND.(use_DREST.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.14).AND.(bnd_type.ne.3).AND.(bnd_type.ne.4)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.15).AND.(ideal_run.EQV..true.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.16).AND.(use_POLY.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.17).AND.(use_BOND(1).EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.18).AND.(use_BOND(2).EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.19).AND.(use_BOND(3).EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.20).AND.(use_BOND(4).EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.22).AND.(use_BOND(5).EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.21).AND.(use_DSSP.EQV..false.)) hmjam%lst(i) = 0
    if ((hmjam%lst(i).eq.23).AND.(use_EMICRO.EQV..false.)) hmjam%lst(i) = 0
    if (hmjam%lst(i).eq.2) hmjam%lst(i) = 13 ! hijack code 13 for total energy 
  end do
!
  if (allocated(hmjam%thresh).EQV..false.) allocate(hmjam%thresh(MAXENERGYTERMS,2))
  if (allocated(hmjam%alpha).EQV..false.) allocate(hmjam%alpha(MAXENERGYTERMS,2))
  hmjam%thresh(:,1) = -1.0*(HUGE(dum)/1000.0)
  hmjam%thresh(:,2) = HUGE(dum)
  hmjam%alpha(:,:) = gasconst*kelvin*10.0 ! 10kT as default
!
  k = 0
  do i=1,hmjam%nlst
    if (hmjam%lst(i).le.1) cycle
    k = k + 1
    if (tmpa(i,1).gt.tmpa(i,2)) then
      write(ilog,*) 'Fatal. For energy landscape sculpting, an inconsistent pair of energy thresholds was provided for &
 &dimension with code ',hmjam%lst(i),'. Please check key-file (FMCSC_ELS_FILLS and FMCSC_ELS_SHAVES).'
      call fexit()
    end if
    if (hmjam%lst(i).eq.2) hmjam%lst(i) = 13 ! hijack code 13 (becomes 11) for total energy 
    hmjam%lst(i) = hmjam%lst(i)-2
    if (hmjam%isin(hmjam%lst(i)).EQV..true.) cycle ! idiotic
    hmjam%isin(hmjam%lst(i)) = .true.
    hmjam%thresh(hmjam%lst(i),1:2) = tmpa(i,1:2)
    if (tmpa(i,3).ge.0.0) hmjam%alpha(hmjam%lst(i),1) = tmpa(i,3)
    if (tmpa(i,4).ge.0.0) hmjam%alpha(hmjam%lst(i),2) = tmpa(i,4)
  end do
!
  if ((k.gt.1).AND.(hmjam%isin(11).EQV..true.)) then
    write(ilog,*) 'Fatal. For energy landscape sculpting, total energy is only available if it is the sole dimension.'
    call fexit()
  end if
!
  if ((use_hardsphere.EQV..true.).AND.((hmjam%isin(3).EQV..true.).OR.(hmjam%isin(11).EQV..true.))) then
    write(ilog,*) 'Warning. The simultaneous use of hard sphere potentials and energy landscape sculpting may produce &
 &nonsensical results.'
  end if
!
  if ((use_stericscreen.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
    write(ilog,*) 'Fatal. In Monte Carlo simulations, it is generally unsafe to not use an energetic pre-screen if the &
 &energy landscape sculpting method is in use. Any term that could produce very large positive energies (without sculpting) &
 &may cause the simulation to become meaningless due to numerical accuracy loss. This termination can be circumvented &
 &with keyword FMCSC_UNSAFE.'
    if (be_unsafe.EQV..false.) call fexit()
  end if
!
  if (use_dyn.EQV..true.) then
    if ((hmjam%isin(9).EQV..true.).OR.(hmjam%isin(6).EQV..true.)) then
      write(ilog,*) 'Fatal. For energy landscape sculpting in gradient-based simulations, all nonbonded interactions must be &
 &treated as a joint (sum) term and cannot be sculpted independently (only use code 15).'
      call fexit()
    end if
    if ((hmjam%isin(1).EQV..true.).OR.(hmjam%isin(3).EQV..true.).OR.(hmjam%isin(5).EQV..true.)) then
      write(ilog,*) 'Fatal. For energy landscape sculpting in gradient-based simulations, all nonbonded interactions must be &
 &treated as a joint (sum) term and cannot be sculpted independently (only use code 15).'
      call fexit()
    end if
  end if
!
  if (k.le.0) then
    write(ilog,*) 'Warning. Disabling energy landscape sculpting method due to lack of any eligible dimensions.'
    do_accelsim = .false.
    hmjam%prtfrmwts = HUGE(hmjam%prtfrmwts)
    if (allocated(hmjam%lst).EQV..true.) deallocate(hmjam%lst)
    return
  end if
!
  if (use_CORR.EQV..true.) then
    write(ilog,*) 'Fatal. The energy landscape sculpting method is incompatible with the quasi-obsolete electronic &
 &correction potential (FMCSC_SC_EXTRA). Use FMCSC_SC_BONDED_T instead.'
    call fexit()
  end if
!
  allocate(hmjam%boosts(MAXENERGYTERMS,2))
  hmjam%boosts(:,:) = 0.0
  if (use_dyn.EQV..true.) then
    allocate(hmjam%ca_f(n,3))
    allocate(hmjam%t_ca_f(3,n))
    hmjam%ca_f(:,:) = 0.0
    hmjam%t_ca_f(:,:) = 0.0
    if (hmjam%isin(6).EQV..true.) then
      if (use_cutoffs.EQV..true.) then
        allocate(hmjam%ca_f_tr(n,3))
        hmjam%ca_f_tr(:,:) = 0.0
      end if
      allocate(hmjam%ca_f_bu(n,3))
      hmjam%ca_f_bu(:,:) = 0.0
    end if
    allocate(hmjam%evec_thr(MAXENERGYTERMS,max(thrdat%maxn,2)))
    hmjam%evec_thr(:,:) = 0.0
  end if
!
  deallocate(tmpa)
  if (allocated(hmjam%lst).EQV..true.) deallocate(hmjam%lst)
!
 77 format('Term ',a,' with thresholds of ',g11.4,' and ',g11.4,' and buffer values (alpha) of ',g11.4,' and ',g11.4,'.')
  write(ilog,*) 
  write(ilog,*) '--- Summary of Sculpted Energy Terms (FMCSC_SCULPT) ---'
  write(ilog,*)
  if (hmjam%isin(11).EQV..true.) write(ilog,77) 'Total Energy ',hmjam%thresh(11,1),hmjam%thresh(11,2),hmjam%alpha(11,1),&
 &                                                hmjam%alpha(11,2)
  if (use_dyn.EQV..true.) then
    if (hmjam%isin(1).EQV..true.) write(ilog,77) 'Joint1 (IPP/ATTLJ/WCA) ',hmjam%thresh(1,1),hmjam%thresh(1,2),&
 &                                                hmjam%alpha(1,1),hmjam%alpha(1,2)
  else
    if (hmjam%isin(1).EQV..true.) write(ilog,77) 'IPP ',hmjam%thresh(1,1),hmjam%thresh(1,2),hmjam%alpha(1,1),hmjam%alpha(1,2)
    if (hmjam%isin(3).EQV..true.) write(ilog,77) 'ATTLJ ',hmjam%thresh(3,1),hmjam%thresh(3,2),hmjam%alpha(3,1),hmjam%alpha(3,2)
    if (hmjam%isin(5).EQV..true.) write(ilog,77) 'WCA ',hmjam%thresh(5,1),hmjam%thresh(5,2),hmjam%alpha(5,1),hmjam%alpha(5,2)
  end if
  if (use_dyn.EQV..true.) then
    if (hmjam%isin(6).EQV..true.) write(ilog,77) 'Joint2 (TABUL/POLAR) ',hmjam%thresh(6,1),hmjam%thresh(6,2),&
 &                                                hmjam%alpha(6,1),hmjam%alpha(6,2)
  else
    if (hmjam%isin(9).EQV..true.) write(ilog,77) 'TABUL ',hmjam%thresh(9,1),hmjam%thresh(9,2),hmjam%alpha(9,1),hmjam%alpha(9,2)
    if (hmjam%isin(6).EQV..true.) write(ilog,77) 'POLAR ',hmjam%thresh(6,1),hmjam%thresh(6,2),hmjam%alpha(6,1),hmjam%alpha(6,2)
  end if
  if (hmjam%isin(2).EQV..true.) write(ilog,77) 'EXTRA ',hmjam%thresh(2,1),hmjam%thresh(2,2),hmjam%alpha(2,1),hmjam%alpha(2,2)
  if (hmjam%isin(4).EQV..true.) write(ilog,77) 'IMPSOLV ',hmjam%thresh(4,1),hmjam%thresh(4,2),hmjam%alpha(4,1),hmjam%alpha(4,2)
  if (hmjam%isin(7).EQV..true.) write(ilog,77) 'TOR ',hmjam%thresh(7,1),hmjam%thresh(7,2),hmjam%alpha(7,1),hmjam%alpha(7,2)
  if (hmjam%isin(8).EQV..true.) write(ilog,77) 'ZSEC ',hmjam%thresh(8,1),hmjam%thresh(8,2),hmjam%alpha(8,1),hmjam%alpha(8,2)
  if (hmjam%isin(10).EQV..true.) write(ilog,77) 'DREST ',hmjam%thresh(10,1),hmjam%thresh(10,2),hmjam%alpha(10,1),hmjam%alpha(10,2)
  if (hmjam%isin(12).EQV..true.)write(ilog,77) 'SOFTWALL ',hmjam%thresh(12,1),hmjam%thresh(12,2),hmjam%alpha(12,1),hmjam%alpha(12,2)
  if (hmjam%isin(13).EQV..true.) write(ilog,77) 'Joint3 (nonbonded) ',hmjam%thresh(13,1),hmjam%thresh(13,2),&
 &                                                                    hmjam%alpha(13,1),hmjam%alpha(13,2)
  if (hmjam%isin(14).EQV..true.) write(ilog,77) 'POLY ',hmjam%thresh(14,1),hmjam%thresh(14,2),hmjam%alpha(14,1),hmjam%alpha(14,2)
  if (hmjam%isin(15).EQV..true.)write(ilog,77) 'BONDED_B ',hmjam%thresh(15,1),hmjam%thresh(15,2),hmjam%alpha(15,1),hmjam%alpha(15,2)
  if (hmjam%isin(16).EQV..true.)write(ilog,77) 'BONDED_A ',hmjam%thresh(16,1),hmjam%thresh(16,2),hmjam%alpha(16,1),hmjam%alpha(16,2)
  if (hmjam%isin(17).EQV..true.)write(ilog,77) 'BONDED_I ',hmjam%thresh(17,1),hmjam%thresh(17,2),hmjam%alpha(17,1),hmjam%alpha(17,2)
  if (hmjam%isin(18).EQV..true.)write(ilog,77) 'BONDED_T ',hmjam%thresh(18,1),hmjam%thresh(18,2),hmjam%alpha(18,1),hmjam%alpha(18,2)
  if (hmjam%isin(19).EQV..true.) write(ilog,77) 'DSSP ',hmjam%thresh(19,1),hmjam%thresh(19,2),hmjam%alpha(19,1),hmjam%alpha(19,2)
  if (hmjam%isin(20).EQV..true.)write(ilog,77) 'BONDED_M ',hmjam%thresh(20,1),hmjam%thresh(20,2),hmjam%alpha(20,1),hmjam%alpha(20,2)
  if (hmjam%isin(21).EQV..true.) write(ilog,77) 'EMICRO ',hmjam%thresh(21,1),hmjam%thresh(21,2),hmjam%alpha(21,1),hmjam%alpha(21,2)
  write(ilog,*)
  write(ilog,*) '--- End of Summary of Sculpted Energy Terms ---'
  write(ilog,*)  
!
end
!
!----------------------------------------------------------------------------------
!
! compute bias terms to unbiased energy vector
!
subroutine els_manage_justE(eunbias,mode)
!
  use energies
  use wl
  use iounit
  use system
!
  implicit none
!
  integer dimi,mode
  RTYPE eunbias(MAXENERGYTERMS),dum,eref
!
  if ((mode.lt.1).OR.(mode.gt.2)) then
    write(ilog,*) 'Called els_manage_justE(...) with unsupported mode. This is a bug.'
    call fexit()
  end if
  hmjam%boosts(:,mode) = 0.0
!
  do dimi=1,MAXENERGYTERMS
    if (hmjam%isin(dimi).EQV..false.) cycle
    eref = eunbias(dimi)
    if (use_dyn.EQV..true.) then! hybrid run
!     generate the required composite terms (otherwise inconsistent Hamiltonian relative to MD)
      if (dimi.eq.13) then
        eref = eunbias(1) + eunbias(3) + eunbias(5) + eunbias(6) + eunbias(9)
      end if
    end if
    if (hmjam%isin(11).EQV..true.) eref = sum(eunbias(:))
    if (eref.lt.hmjam%thresh(dimi,1)) then
      dum = hmjam%thresh(dimi,1) - eref
      if (hmjam%alpha(dimi,1).le.0.0) then
        hmjam%boosts(dimi,mode) = dum
      else
        hmjam%boosts(dimi,mode) = dum*dum/(hmjam%alpha(dimi,1) + dum)
      end if
    else if (eref.gt.hmjam%thresh(dimi,2)) then
      dum = hmjam%thresh(dimi,2) - eref
      if (hmjam%alpha(dimi,2).le.0.0) then
        hmjam%boosts(dimi,mode) = dum
      else
        hmjam%boosts(dimi,mode) = dum*dum/(dum - hmjam%alpha(dimi,2))
      end if
    end if
  end do
!
end
!
!----------------------------------------------------------------------------------
!
! compute bias terms to unbiased energy for an individual component and manage forces
! the critical (and difficult) bit is that hmajm%ca_f must contain the partial force
! for the boosted energy term
! unlike els_manage_justE, this fxn also increments biased energy right away
!
subroutine els_manage_one(dimi,ebias,ca_f,mode)
!
  use energies
  use wl
  use iounit
  use atoms
!
  implicit none
!
  integer dimi,mode
  RTYPE ebias,dum,ca_f(n,3),dbdU
!
  if ((mode.lt.1).OR.(mode.gt.3)) then
    write(ilog,*) 'Called els_manage_one(...) with unsupported mode. This is a bug.'
    call fexit()
  end if
!
  hmjam%boosts(dimi,1) = 0.0
!
  if (hmjam%isin(dimi).EQV..false.) return
!
  if (ebias.lt.hmjam%thresh(dimi,1)) then
    dum = hmjam%thresh(dimi,1) - ebias
    if (hmjam%alpha(dimi,1).le.0.0) then
      hmjam%boosts(dimi,1) = dum
    else
      hmjam%boosts(dimi,1) = dum*dum/(hmjam%alpha(dimi,1) + dum)
    end if
    dbdU = 1.0 - 2.0*dum/(hmjam%alpha(dimi,1) + dum) + dum*dum/((hmjam%alpha(dimi,1) + dum)**2)
  else if (ebias.gt.hmjam%thresh(dimi,2)) then
    dum = hmjam%thresh(dimi,2) - ebias
    if (hmjam%alpha(dimi,2).le.0.0) then
      hmjam%boosts(dimi,1) = dum
    else
      hmjam%boosts(dimi,1) = dum*dum/(dum - hmjam%alpha(dimi,2))
    end if
    dbdU = 1.0 - 2.0*dum/(dum - hmjam%alpha(dimi,2)) + dum*dum/((dum - hmjam%alpha(dimi,2))**2)
  else
    dbdU = 1.0
  end if
!
  if (mode.eq.1) then ! increment force passed as argument and zero hmjam%ca_f
    ca_f(:,:) = ca_f(:,:) + hmjam%ca_f(:,:)*dbdU
    hmjam%ca_f(:,:) = 0.0
  else if (mode.eq.2) then ! increment force passed as argument with bias force only
    ca_f(:,:) = ca_f(:,:) + hmjam%ca_f(:,:)*(dbdU - 1.0)
    hmjam%ca_f(:,:) = 0.0
  else if (mode.eq.3) then ! overwrite force passed as argument with net bias derived from itself
    ca_f(:,:) = ca_f(:,:)*dbdU
  end if
  ebias = ebias + hmjam%boosts(dimi,1)
!
end

!
!----------------------------------------------------------------------------------
!
! the same for transpose
!
subroutine els_manage_one_t(dimi,ebias,ca_f,mode)
!
  use energies
  use wl
  use iounit
  use atoms
  use cutoffs, ONLY: t_ca_f,t_ca_f_tr
  use forces, ONLY: cart_f,cart_f_tr
!
  implicit none
!
  integer dimi,mode
  RTYPE ebias,dum,ca_f(3,n),dbdU
!
  if ((mode.lt.1).OR.(mode.gt.4)) then
    write(ilog,*) 'Called els_manage_one_t(...) with unsupported mode. This is a bug.'
    call fexit()
  end if
!
  hmjam%boosts(dimi,1) = 0.0
!
  if (hmjam%isin(dimi).EQV..false.) return
!
  if (ebias.lt.hmjam%thresh(dimi,1)) then
    dum = hmjam%thresh(dimi,1) - ebias
    if (hmjam%alpha(dimi,1).le.0.0) then
      hmjam%boosts(dimi,1) = dum
    else
      hmjam%boosts(dimi,1) = dum*dum/(hmjam%alpha(dimi,1) + dum)
    end if
    dbdU = 1.0 - 2.0*dum/(hmjam%alpha(dimi,1) + dum) + dum*dum/((hmjam%alpha(dimi,1) + dum)**2)
  else if (ebias.gt.hmjam%thresh(dimi,2)) then
    dum = hmjam%thresh(dimi,2) - ebias
    if (hmjam%alpha(dimi,2).le.0.0) then
      hmjam%boosts(dimi,1) = dum
    else
      hmjam%boosts(dimi,1) = dum*dum/(dum - hmjam%alpha(dimi,2))
    end if
    dbdU = 1.0 - 2.0*dum/(dum - hmjam%alpha(dimi,2)) + dum*dum/((dum - hmjam%alpha(dimi,2))**2)
  else
    dbdU = 1.0
  end if
!
  if (mode.eq.1) then ! increment force passed as argument and zero hmjam%t_ca_f
    ca_f(:,:) = ca_f(:,:) + hmjam%t_ca_f(:,:)*dbdU
    hmjam%t_ca_f(:,:) = 0.0
    ebias = ebias + hmjam%boosts(dimi,1)
  else if (mode.eq.2) then ! increment force passed as argument with bias force only
    ca_f(:,:) = ca_f(:,:) + hmjam%t_ca_f(:,:)*(dbdU - 1.0)
    hmjam%t_ca_f(:,:) = 0.0
    ebias = ebias + hmjam%boosts(dimi,1)
  else if (mode.eq.3) then ! overwrite force passed as argument with net bias derived from itself
    ca_f(:,:) = ca_f(:,:)*dbdU
    ebias = ebias + hmjam%boosts(dimi,1)
  else if (mode.eq.4) then
    ca_f(:,:) = ca_f(:,:) + (t_ca_f(:,:)+t_ca_f_tr(:,:))*(dbdU-1.0) + transpose(cart_f(:,:)+cart_f_tr(:,:))*(dbdU-1.0)
    t_ca_f(:,:) = t_ca_f(:,:) + ca_f(:,:)
    ca_f(:,:) = 0.0
    ebias = hmjam%boosts(dimi,1)
  end if
!
end

!
!----------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
! for threads, we cannot really afford to make that many copies of cart_f -> slightly
! different strategy consists of merging information from thr_ca_f into t_ca_f, which 
! is otherwise unused in threaded force evaluation
!
subroutine els_manage_one_t_threads(dimi,mode,eref,tpi)
!
  use energies
  use wl
  use iounit
  use atoms
  use cutoffs, ONLY: t_ca_f
  use threads
  use forces, ONLY: cart_f,cart_f_tr
!
  implicit none
!
  integer, INTENT(IN):: tpi,mode,dimi
!
  integer i,j
  RTYPE eraw,eref,dum,dbdU
!
  logical OMP_IN_PARALLEL
!  integer(KIND=8) ttt(2)
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, els_manage_one_t_threads(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if

  if ((mode.lt.0).OR.(mode.gt.3)) then
    write(ilog,*) 'Called els_manage_one_t_threads(...) with unsupported mode. This is a bug.'
    call fexit()
  end if
!$OMP BARRIER
!
! simply store force in thr_ca_f so that we can use it for incremental forces
  if (mode.eq.0) then
    do i=0,thrdat%maxn-1
      j = mod(i+tpi,thrdat%maxn) + 1
      t_ca_f(:,thr_limits(1,j):thr_limits(2,j)) = t_ca_f(:,thr_limits(1,j):thr_limits(2,j)) + &
 &                                        thr_ca_f(:,thr_limits(1,j):thr_limits(2,j),tpi)
  !$OMP BARRIER
    end do
    thr_ca_f(:,:,tpi) = 0.0
    return
  end if
!
  eraw = sum(hmjam%evec_thr(dimi,:))
!
  if (tpi.le.1) hmjam%boosts(dimi,1) = 0.0
!
  if (hmjam%isin(dimi).EQV..false.) return
!
! energy bias is computed by every thread, but increment of %boosts is restricted
  if (eraw.lt.hmjam%thresh(dimi,1)) then
    dum = hmjam%thresh(dimi,1) - eraw
    if (hmjam%alpha(dimi,1).le.0.0) then
      if (tpi.le.1) hmjam%boosts(dimi,1) = dum
    else
      if (tpi.le.1) hmjam%boosts(dimi,1) = dum*dum/(hmjam%alpha(dimi,1) + dum)
    end if
    dbdU = 1.0 - 2.0*dum/(hmjam%alpha(dimi,1) + dum) + dum*dum/((hmjam%alpha(dimi,1) + dum)**2)
  else if (eraw.gt.hmjam%thresh(dimi,2)) then
    dum = hmjam%thresh(dimi,2) - eraw
    if (hmjam%alpha(dimi,2).le.0.0) then
      if (tpi.le.1) hmjam%boosts(dimi,1) = dum
    else
      if (tpi.le.1) hmjam%boosts(dimi,1) = dum*dum/(dum - hmjam%alpha(dimi,2))
    end if
    dbdU = 1.0 - 2.0*dum/(dum - hmjam%alpha(dimi,2)) + dum*dum/((dum - hmjam%alpha(dimi,2))**2)
  else
    dbdU = 1.0
  end if
!
! the update loop involves N_threads barriers
  if (mode.eq.1) then
    do i=0,thrdat%maxn-1
      j = mod(i+tpi,thrdat%maxn) + 1
      t_ca_f(:,thr_limits(1,j):thr_limits(2,j)) = t_ca_f(:,thr_limits(1,j):thr_limits(2,j)) + &
 &                                        thr_ca_f(:,thr_limits(1,j):thr_limits(2,j),tpi)*dbdU
!$OMP BARRIER
    end do
    thr_ca_f(:,:,tpi) = 0.0
  else if (mode.eq.2) then ! tally up everything but subtract the bias contribution only to avoid having to zero the force
    do i=0,thrdat%maxn-1
      j = mod(i+tpi,thrdat%maxn) + 1
      t_ca_f(:,thr_limits(1,j):thr_limits(2,j)) = t_ca_f(:,thr_limits(1,j):thr_limits(2,j)) - &
 &     (thr_ca_f_tr(:,thr_limits(1,j):thr_limits(2,j),tpi)+thr_ca_f(:,thr_limits(1,j):thr_limits(2,j),tpi))*(1.0-dbdU)
!$OMP BARRIER
    end do
    t_ca_f(:,thr_limits(1,tpi):thr_limits(2,tpi)) = t_ca_f(:,thr_limits(1,tpi):thr_limits(2,tpi)) - &
 &    transpose(cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:)+cart_f_tr(thr_limits(1,tpi):thr_limits(2,tpi),:))*(1.0-dbdU)
!$OMP BARRIER
    hmjam%evec_thr(:,tpi) = 0.0
  else if (mode.eq.3) then ! transform in place
    cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:) = cart_f(thr_limits(1,tpi):thr_limits(2,tpi),:)*dbdU
!$OMP BARRIER
  end if
!
  if (tpi.le.1) then
    if (mode.eq.1) then
      eref = eraw + hmjam%boosts(dimi,1)
    else if (mode.eq.2) then
      eref = hmjam%boosts(dimi,1)
    else if (mode.eq.3) then
      eref = eraw + hmjam%boosts(dimi,1)
    end if
  end if
!
end
!
#endif
!
!----------------------------------------------------------------------------------
!
subroutine els_prt_wts()
!
  use wl
  use mcsums
  use system
!
  implicit none
!
  RTYPE netboost
!
  netboost = exp(invtemp*sum(hmjam%boosts(:,1)))
  if (netboost.le.hmjam%threshwt) return
! 
 333 format(i14,1x,g19.11E3)
  write(hmjam%iwtsfile,333) nstep,netboost
!
end
!
!------------------------------------------------------------------------------------
!
