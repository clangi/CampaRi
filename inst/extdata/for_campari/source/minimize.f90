!-------------------------------------------------------------------------!
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
! MAIN AUTHOR:   Adam Steffen, Andreas Vitalis                             !
! CONTRIBUTIONS: Nicholas Lyle                                             !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! ##############################################################
! ##                                                          ##
! ## subroutine minimize -- Minimization over entire system   ##
! ##                                                          ##
! ##############################################################
!    
!
subroutine minimize(maxiter, grmsmin, istepsize, minimode, tpi)

  use iounit
  use molecule
  use forces
  use energies
  use movesets
  use system
  use mcsums
  use fyoc
  use cutoffs
  use mini
  use zmatrix
  use math
  use atoms
#ifdef ENABLE_THREADS
  use threads
#endif
!      
  implicit none
!      
  integer, INTENT(IN):: minimode,maxiter,tpi
!
  integer azero,iter,nrdof,imol,j,i,j1,j2,t,pmode,iter_bfgs,wrnlmt(2),wrncnt(2),iter_avg,stx(2)
  integer ndump,mstore,incr,bound,evals,redux_count,redux_max,upaccept,badaccept,avgsz
  RTYPE grmsmin,grms,stepsize,eprevious,istepsize
  RTYPE force3,scalefactor,sigma,b_bfgs,line_grad,hnot_bfgs
  logical atrue,bfgsfail
  integer(KIND=8) t1,t2
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
  integer stx2(2)
  RTYPE dumm,dumm2
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, minimize(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
! "vectorize" all degrees of freedom into a single vector with dimension "nrdof"
  nrdof = 0
  if (fycxyz.eq.1) then
    do imol=1,nmol
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet supported in minimization runs.'
        call fexit()
      end if
      if (tpi.le.1) molinfo(imol,3) = nrdof + 1
      do j=1,3
        if (dc_di(imol)%frz(j).EQV..true.) cycle
        nrdof = nrdof + 1
      end do
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            j = i+6 ! RB rotation
          else
            j = dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (j.le.0) cycle
          if (dc_di(imol)%frz(j).EQV..true.) cycle
          nrdof = nrdof + 1
        end do
      end if
    end do 
  else if (fycxyz.eq.2) then
! SHAKE-type constraints and drift removal must be handled separately (not mappable)
    do imol=1,nmol
      if (tpi.le.1) molinfo(imol,3) = nrdof + 1
      do i=atmol(imol,1),atmol(imol,2)
        do j=1,3
          if (cart_frz(i,j).EQV..true.) cycle
          nrdof = nrdof + 1
        end do
      end do
    end do
  else
    write(ilog,*) 'Fatal. Called minimize(...) with unsupported system representation. This is an omission bug.'
    call fexit()
  end if
!
  if (nrdof.le.0) then
    write(ilog,*) 'Fatal. No remaining degrees of freedom to minimize. Relieve constraints.'
    call fexit()
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    call threads_bounds(1,nrdof,tpi,thrdat%maxn,stx)
  else
    stx(1) = 1
    stx(2) = nrdof
  end if
#else
  stx(1) = 1
  stx(2) = nrdof
#endif
!     
10   format(i8,1x,i8,1x,20000(f18.10,1x))
!
! possible globals
  mstore = min(mini_mem,10) ! threshold for accepting uphill steps in BFGS (max 10, min memory)
  redux_max = 50  ! number of stepsize divisions tried before exiting (usually indicates uphill attempt)
  avgsz = 100
!
! initializations
  azero = 0
  evals = 0
  ndump = 0
  grms = 0.0
  redux_count = 0
  iter = 0
  iter_bfgs = 0
  iter_avg = 0
  upaccept = 0
  badaccept = 0
  wrncnt(:) = 0
  wrnlmt(:) = 1
  atrue = .true.
  eprevious = esave
  stepsize = 1.0
  bfgsfail = .false.
!
#ifdef ENABLE_THREADS
!$OMP SINGLE 
#endif
  do imol=1,nmol
    call update_rigidm(imol)
  end do 
!
! allocations used by all methods
  allocate(minivs%cgrv(nrdof))
  allocate(minivs%pgrv(nrdof))
  allocate(minivs%dposv(nrdof))
  allocate(minivs%pcgn(nrdof))
  allocate(minivs%peakg(nrdof))
  allocate(minivs%aheadg(nrdof))
  allocate(minivs%mshfs(3,nmol))
  allocate(minivs%cgn(nrdof))
  allocate(minivs%avstp(avgsz))
  minivs%dposv(:) = 0.0
  minivs%cgn(:) = 0.0
  minivs%pcgn(:) = 0.0
  minivs%cgrv(:) = 0.0
  minivs%avstp(:) = 0.0
!  
! quasi-newton BFGS allocations/inits
  if(minimode.eq.3) then
    allocate(minivs%b_s(nrdof,mini_mem))
    allocate(minivs%b_y(nrdof,mini_mem))
    allocate(minivs%b_a(mini_mem))
    allocate(minivs%b_q(nrdof))
    allocate(minivs%b_r(nrdof))
    minivs%b_s(:,:) = 0.0
    minivs%b_y(:,:) = 0.0
    minivs%b_a(:) = 0.0
    minivs%b_q(:) = 0.0
    minivs%b_r(:) = 0.0
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE 
#endif
!     
! ---------- MASTER LOOP -------------
!
  do while (evals.le.maxiter)
!
!   count the number of tried energies (peak and reduce method)
    evals = evals + 1
!
!   zero out internal coordinate forces
    if (fycxyz.eq.1) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        stx2(1) = thr_limits(61,tpi)
        stx2(2) = thr_limits(62,tpi)
      else
        stx2(1) = 1
        stx2(2) = nmol
      end if
      do imol=stx2(1),stx2(2)
#else
      do imol=1,nmol
#endif
        dc_di(imol)%f(:) = 0.0
      end do
    end if ! cart_f is zeroed elsewhere 
!
!   calculate energies and gradient at test point
    if (tpi.le.1) call System_Clock(t1)
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
!$OMP BARRIER
      call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
    else
#endif
    esave = force3(esterms,esterms_tr,esterms_lr,atrue)
#ifdef ENABLE_THREADS
    end if
#endif
    if (tpi.le.1) then
      ens%insR(10) = ens%insU
      ens%insU = esave
      call System_Clock(t2)
      time_energy = time_energy + 1.0*(t2 - t1)
    end if

    if (fycxyz.eq.1) then
      call cart2int_f(skip_frz,tpi)
    end if
    minivs%peakg(stx(1):stx(2)) = 0.0
    minivs%aheadg(stx(1):stx(2)) = 0.0
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
    call get_gradient(minivs%aheadg,minivs%grms,nrdof,tpi,stx)
!
    if (minimode.gt.1) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
!$OMP SINGLE
        minivs%hlp(3) = 0.0
!$OMP END SINGLE
        dumm = sum(minivs%aheadg(stx(1):stx(2))*minivs%cgrv(stx(1):stx(2)))
!$OMP CRITICAL(MINI_LINEG_UP)
        minivs%hlp(3) = minivs%hlp(3) + dumm
!$OMP END CRITICAL(MINI_LINEG_UP)
!$OMP BARRIER
        line_grad = minivs%hlp(3)
      else
#endif 
      line_grad = DOT_PRODUCT(minivs%aheadg,minivs%cgrv)
#ifdef ENABLE_THREADS
      end if
#endif
    end if
!
!   if the energy is lower (or nearly lower) then move on to a new point
    if ((evals.eq.1).OR.&
 &     ((minimode.eq.1).AND.(esave.lt.eprevious)).OR.&
 &     ((minimode.eq.2).AND.(esave.lt.eprevious).AND.(line_grad.gt.0.0)).OR.&
 &     ((minimode.eq.3).AND.(esave.lt.(eprevious + mini_uphill)).AND.((line_grad.gt.0.0).OR.(iter_bfgs.eq.1)))) then

      iter = iter + 1
      iter_bfgs = iter_bfgs + 1
      if (iter_avg.gt.0) then
        if (tpi.le.1) minivs%avstp(mod(iter_avg,avgsz)+1) = stepsize
      end if
      iter_avg = iter_avg + 1
!
      if (tpi.le.1) nstep = nstep + 1
!
      pmode = 0
      minivs%peakg(stx(1):stx(2)) = minivs%aheadg(stx(1):stx(2))
      grms = minivs%grms
!      write(ilog,10) evals,iter,stepsize,esave,grms
      if (tpi.le.1) call System_Clock(t1)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      call mcstat(iter,ndump,tpi)
      if (tpi.le.1) then
        call System_Clock(t2)
        time_analysis = time_analysis + 1.0*(t2 - t1)
      end if
!
!     store gradient information
      minivs%pgrv(stx(1):stx(2)) = minivs%cgrv(stx(1):stx(2))
      minivs%pcgn(stx(1):stx(2)) = minivs%cgn(stx(1):stx(2))
      minivs%cgrv(stx(1):stx(2)) = minivs%peakg(stx(1):stx(2))
      eprevious = esave
!
!     BFGS: update gradient stored value for l-bfgs
      if (minimode.eq.3) then
        j = MOD(iter_bfgs,mini_mem)
        if (j.eq.0) then
          j = mini_mem
        end if
        minivs%b_y(stx(1):stx(2),j) = minivs%cgrv(stx(1):stx(2))
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      end if
!  
!     gradient RMS convergence check
      if (grms.lt.grmsmin) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        write(ilog,*) 'Termination: Gradient RMS convergence reached.'
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
        exit
      end if
      if(mod(iter,nsancheck).eq.0) then
 2  format('Step: ',i8,' Current Energy: ',g14.6,' G-RMS: ',g14.6,' Avg Step Scale: ',g14.6)
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        write(ilog,2) iter,esave,grms,(1.0/(1.0*min(iter_avg,avgsz)))*sum(minivs%avstp(1:min(iter_avg,avgsz)))
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
      end if
!
!     increase the stepsize if sequential minimization moves are accepted, reset if last step was rejected
!     note that due to step-size reductions it's rather likely for the energies in subsequent steps to be floating-point
!     identical (.ge.)
      if (esave.ge.eprevious) then
        badaccept = badaccept + 1
        if (esave.gt.eprevious) upaccept = upaccept + 1
!        
!       warn against real uphill problems
        if ((upaccept.gt.mstore)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
          wrncnt(2) = wrncnt(2) + 1
          if (wrncnt(2).eq.wrnlmt(2)) then
            write(ilog,*) 'WARNING: Minimization has exceeded the permissible number of uphill steps.&
 & This means that the current direction(s) was (were) erroneous.'
            write(ilog,*) 'This was warning #',wrncnt(2),' of this type not all of which may be displayed.'
            wrnlmt(2) = wrnlmt(2)*10
          end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
        end if
!
!       reset
        if (badaccept.gt.mstore) then
          upaccept = 0
          badaccept = 0
          minivs%cgrv(stx(1):stx(2)) = minivs%aheadg(stx(1):stx(2))
          if (minimode.eq.2) then
            minivs%pgrv(stx(1):stx(2)) = minivs%cgrv(stx(1):stx(2))
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(minivs%avstp(1:min(iter_avg,avgsz)))
          end if
          if (minimode.eq.3) then
            minivs%b_y(stx(1):stx(2),:) = 0.0
            minivs%b_y(stx(1):stx(2),1) = minivs%cgrv(stx(1):stx(2))
            minivs%b_s(stx(1):stx(2),:) = 0.0
            iter_bfgs = 1
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(minivs%avstp(1:min(iter_avg,avgsz)))
          end if
        end if
      else
        upaccept = 0
        badaccept = 0
      end if
      if ((redux_count.eq.0).OR.(minimode.eq.2)) then
        stepsize = stepsize*1.618034
      else
        stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(minivs%avstp(1:min(iter_avg,avgsz)))
      end if
      redux_count = 0
!
!   If peak energy is not lower, reset structure and try a smaller step
    else
!
!     reduce step size
      stepsize = stepsize*0.618034
      redux_count = redux_count + 1
      pmode = 1
!
!     explicitly undo last move
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call set_position_threads(pmode,nrdof,tpi,stx)
      else
#endif
      call set_position(pmode,nrdof)
#ifdef ENABLE_THREADS
      end if
#endif
      pmode = 0
!
!     give up if max reduction in stepsize is reached (usually indicates an uphill attempt for CG and BFGS)
      if (redux_count.gt.redux_max) then
        if ((minimode.eq.2).AND.(grms.gt.(1.0*grmsmin))) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
          wrncnt(1) = wrncnt(1) + 1
          if (wrncnt(1).eq.wrnlmt(1)) then
            write(ilog,*) 'WARNING: The CG method has exceeded the permissible number of stepsize reductions.&
 & Resetting direction to steepest-descent and continuing minimization since GRMS is far away from requested value.'
            write(ilog,*) 'This was warning #',wrncnt(1),' of this type not all of which may be displayed.'
            wrnlmt(1) = wrnlmt(1)*10
          end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
          minivs%cgrv(stx(1):stx(2)) = minivs%aheadg(stx(1):stx(2))
          minivs%pgrv(stx(1):stx(2)) = minivs%cgrv(stx(1):stx(2))
          minivs%pcgn(stx(1):stx(2)) = 0.0
          if (iter_avg.gt.0) then
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(minivs%avstp(1:min(iter_avg,avgsz)))
          else
            stepsize = 1.0
          end if
          redux_count = 0
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
        else if ((minimode.eq.3).AND.(grms.gt.(1.0*grmsmin))) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
          wrncnt(1) = wrncnt(1) + 1
          if (wrncnt(1).eq.wrnlmt(1)) then
            write(ilog,*) 'WARNING: The BFGS method has exceeded the permissible number of stepsize reductions.&
 & This means either that a local minimum was found or that the current direction was erroneous. Resetting&
 & Hessian and continuing minimization since GRMS is far away from requested value.'
            write(ilog,*) 'This was warning #',wrncnt(1),' of this type not all of which may be displayed.'
            wrnlmt(1) = wrnlmt(1)*10
          end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
          minivs%b_y(stx(1):stx(2),:) = 0.0
          minivs%cgrv(stx(1):stx(2)) = minivs%aheadg(stx(1):stx(2))
          minivs%b_y(stx(1):stx(2),1) = minivs%cgrv(stx(1):stx(2))
          minivs%b_s(stx(1):stx(2),:) = 0.0
          iter_bfgs = 1
          if (iter_avg.gt.0) then
            stepsize = (1.0/(1.0*min(iter_avg,avgsz)))*sum(minivs%avstp(1:min(iter_avg,avgsz)))
          else
            stepsize = 1.0
          end if
          redux_count = 0
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
        else
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
          esave = eprevious
          write(ilog,*) 'Termination: Maximum stepsize reductions reached during minimization. This indicates a&
 &n unreasonably small GRMS request, a nearly singular function in the energy landscape, alternative termination&
 & of the BFGS method (GRMS near target) or possibly a bug.'
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
          exit
        end if
      end if
!
    end if
!
!  
!   Make minimizing move, mode based on minimode:
!
!   1) steepest decent
    if (minimode.eq.1) then
      minivs%dposv(stx(1):stx(2)) = stepsize*istepsize*minivs%cgrv(stx(1):stx(2))
!
!    2) conj. grad. (Polak and Ribiere form)
    else if(minimode.eq.2) then
      if(iter.eq.1) then
        minivs%dposv(stx(1):stx(2)) = stepsize*istepsize*minivs%cgrv(stx(1):stx(2))
      else
#ifdef ENABLE_THREADS
        if (tpi.gt.0) then
!$OMP SINGLE
          minivs%hlp(1:2) = 0.0
!$OMP END SINGLE
          dumm = sum(minivs%cgrv(stx(1):stx(2))*minivs%cgrv(stx(1):stx(2)))
          dumm2 = sum(minivs%pgrv(stx(1):stx(2))*minivs%pgrv(stx(1):stx(2)))
!$OMP CRITICAL(MINI_CG_UP)
          minivs%hlp(1) = minivs%hlp(1) + dumm
          minivs%hlp(2) = minivs%hlp(2) + dumm2
!$OMP END CRITICAL(MINI_CG_UP)
!$OMP BARRIER
          scalefactor = minivs%hlp(1)/minivs%hlp(2)
        else
#endif
        scalefactor = DOT_PRODUCT(minivs%cgrv(stx(1):stx(2)),minivs%cgrv(stx(1):stx(2))) / &
 &                    DOT_PRODUCT(minivs%pgrv(stx(1):stx(2)),minivs%pgrv(stx(1):stx(2)))
#ifdef ENABLE_THREADS
        end if
#endif
        minivs%cgn(stx(1):stx(2)) = minivs%cgrv(stx(1):stx(2)) + scalefactor*minivs%pcgn(stx(1):stx(2))
        minivs%dposv(stx(1):stx(2)) = istepsize*stepsize*minivs%cgn(stx(1):stx(2))
      end if
!
!   3) quasi-newton L-BGFS method (Nocedal method)
    else if(minimode.eq.3) then
      if ((iter_bfgs.eq.1).OR.(iter.eq.1)) then 
        minivs%dposv(stx(1):stx(2)) = stepsize*istepsize*minivs%cgrv(stx(1):stx(2))
      else
!
!       Two Loop Recursion (Nocedal)
        if(iter_bfgs.le.(mini_mem+1)) then 
          incr = 0
          bound = iter_bfgs
        else
          incr = iter_bfgs - mini_mem
          bound = mini_mem
        end if
        minivs%b_q(stx(1):stx(2)) = minivs%cgrv(stx(1):stx(2))
!
        j1 = MOD(iter_bfgs,mini_mem)
        if (j.eq.0) then
          j1 = mini_mem
        end if
        j2 = j1 - 1
        if (j2.eq.0) j2 = mini_mem
!
        do i=bound-1,1,-1
!
          j1 = MOD(i + incr,mini_mem)+1
          j2 = j1-1
          if (j1.eq.1) then
            j2 = mini_mem
          end if
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
!$OMP SINGLE
            minivs%hlp(6) = 0.0
            minivs%b_a(i) = 0.0
!$OMP END SINGLE
            dumm = sum((minivs%b_y(stx(1):stx(2),j1) - (minivs%b_y(stx(1):stx(2),j2)))*minivs%b_s(stx(1):stx(2),j1))
            dumm2 = sum(minivs%b_s(stx(1):stx(2),j1)*minivs%b_q(stx(1):stx(2)))
!$OMP CRITICAL(MINI_BFGS1_UP)
            minivs%hlp(6) = minivs%hlp(6) + dumm
            minivs%b_a(i) = minivs%b_a(i) + dumm2
!$OMP END CRITICAL(MINI_BFGS1_UP)
!$OMP BARRIER
            sigma = minivs%hlp(6)
            if (sigma.eq.0.0) then
              bfgsfail = .true.
              sigma = 1.0e-9
            end if
            minivs%b_q(stx(1):stx(2)) = minivs%b_q(stx(1):stx(2)) - (minivs%b_a(i)/sigma)*&
 &                                     (minivs%b_y(stx(1):stx(2),j1)-minivs%b_y(stx(1):stx(2),j2))
!$OMP BARRIER
!$OMP SINGLE
            minivs%b_a(i) = minivs%b_a(i)/sigma
!$OMP END SINGLE NOWAIT
          else
#endif
          sigma = 0.0
          do t=1,nrdof
            sigma = sigma + (minivs%b_y(t,j1)-minivs%b_y(t,j2))*minivs%b_s(t,j1)
          end do
          if (sigma.eq.0.0) then
            bfgsfail = .true.
            sigma = 1.0e-9
          end if
          minivs%b_a(i) = 0.0
          do t=1,nrdof
            minivs%b_a(i) = minivs%b_a(i) + (1.0/sigma)*minivs%b_s(t,j1)*minivs%b_q(t)
          end do
          minivs%b_q(:) = minivs%b_q(:) - minivs%b_a(i)*(minivs%b_y(:,j1)-minivs%b_y(:,j2))
#ifdef ENABLE_THREADS
          end if
#endif
        end do
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
        minivs%hlp(7) = 0.0
!$OMP END SINGLE
#endif
        hnot_bfgs = sum((minivs%b_y(stx(1):stx(2),j1)-minivs%b_y(stx(1):stx(2),j2))*&
 &                      (minivs%b_y(stx(1):stx(2),j1)-minivs%b_y(stx(1):stx(2),j2)))
#ifdef ENABLE_THREADS
!$OMP CRITICAL(MINI_BFGS2_UP)
        minivs%hlp(7) = minivs%hlp(7) + hnot_bfgs
!$OMP END CRITICAL(MINI_BFGS2_UP)
!$OMP BARRIER
        hnot_bfgs = minivs%hlp(7)
        if (hnot_bfgs.le.0.0) then
          bfgsfail = .true.
          hnot_bfgs = 1.0e-9
        end if
#endif
        minivs%b_r(stx(1):stx(2)) = minivs%b_q(stx(1):stx(2))*sigma/hnot_bfgs
!
!       second loop
        do i=1,bound-1
!
          j1 = MOD(i + incr,mini_mem)+1
          j2 = j1-1
          if (j1.eq.1) then
            j2 = mini_mem
          end if
!
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
!$OMP BARRIER
!$OMP SINGLE
            minivs%hlp(4:5) = 0.0
!$OMP END SINGLE
            dumm = sum((minivs%b_y(stx(1):stx(2),j1) - (minivs%b_y(stx(1):stx(2),j2)))*minivs%b_s(stx(1):stx(2),j1))
            dumm2 = sum((minivs%b_y(stx(1):stx(2),j1) - (minivs%b_y(stx(1):stx(2),j2)))*minivs%b_r(stx(1):stx(2)))
!$OMP CRITICAL(MINI_BFGS3_UP)
            minivs%hlp(4) = minivs%hlp(4) + dumm
            minivs%hlp(5) = minivs%hlp(5) + dumm2
!$OMP END CRITICAL(MINI_BFGS3_UP)
!$OMP BARRIER
            sigma = minivs%hlp(4)
            if (sigma.eq.0.0) then
              bfgsfail = .true.
              sigma = 1.0e-9
            end if
            b_bfgs = minivs%hlp(5)/sigma
          else
#endif
          sigma = 0.0
          do t=1,nrdof
            sigma = sigma + (minivs%b_y(t,j1)-minivs%b_y(t,j2))*(minivs%b_s(t,j1))
          end do
          if (sigma.eq.0.0) then
            bfgsfail = .true.
            sigma = 1.0e-9
          end if
!          write(*,*) 'S2',sigma
          b_bfgs = 0.0
          do t=1,nrdof
            b_bfgs = b_bfgs + (1.0/sigma)*(minivs%b_y(t,j1)-minivs%b_y(t,j2))*minivs%b_r(t)
          end do
#ifdef ENABLE_THREADS
          end if
#endif
          minivs%b_r(stx(1):stx(2)) = minivs%b_r(stx(1):stx(2)) + (minivs%b_a(i)-b_bfgs)*(minivs%b_s(stx(1):stx(2),j1))
        end do
!
!       Normalize to the unit direction
        sigma = -stepsize
#ifdef ENABLE_THREADS
!$OMP SINGLE
        minivs%hlp(8) = 0.0
!$OMP END SINGLE
#endif
        hnot_bfgs = sum(minivs%b_r(stx(1):stx(2))*minivs%b_r(stx(1):stx(2)))
#ifdef ENABLE_THREADS
!$OMP CRITICAL(MINI_BFGS4_UP)
        minivs%hlp(8) = minivs%hlp(8) + hnot_bfgs
!$OMP END CRITICAL(MINI_BFGS4_UP)
!$OMP BARRIER
        hnot_bfgs = sigma/sqrt(minivs%hlp(8))
#else
        hnot_bfgs = sigma/sqrt(hnot_bfgs)
#endif
        minivs%dposv(stx(1):stx(2)) = hnot_bfgs*minivs%b_r(stx(1):stx(2))
      end if
!
    end if ! which minimization algorithm
!
!   BFGS: exit if needed, otherwise update stored position values
    if (minimode.eq.3) then
      if (bfgsfail.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        write(ilog,*) 'WARNING: The BFGS method is trapped irrevocably and will now be terminated. This is not &
 &necessarily a problem and could simply indicate an overly stringent convergence criterion.'
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
        exit
      end if
      j = MOD(iter_bfgs,mini_mem) + 1
      minivs%b_s(stx(1):stx(2),j) = minivs%dposv(stx(1):stx(2))
    end if
!
!   set peak position, backing up current position
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      call set_position_threads(pmode,nrdof,tpi,stx) ! ends with a barrier
    else
#endif
    call set_position(pmode,nrdof)
#ifdef ENABLE_THREADS
    end if
#endif
!
  end do ! iterations of MASTER LOOP
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
! analysis if converged
  if (tpi.le.1) call System_Clock(t1)
  j = nsim  ! ensure printout
  call mcstat(j,ndump,tpi)
  if (tpi.le.1) then
    call System_Clock(t2)
    time_analysis = time_analysis + 1.0*(t2 - t1)
  end if
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
!    
! clean up
  if(minimode.eq.3) then
    deallocate(minivs%b_s)
    deallocate(minivs%b_y)
    deallocate(minivs%b_a)
    deallocate(minivs%b_q)
    deallocate(minivs%b_r)
  end if
  deallocate(minivs%cgrv)
  deallocate(minivs%pgrv)
  deallocate(minivs%dposv)
  deallocate(minivs%pcgn)
  deallocate(minivs%peakg)
  deallocate(minivs%mshfs)
  deallocate(minivs%cgn)
  deallocate(minivs%avstp)
!
! reached max iterations without convergence
  if (evals.gt.maxiter) then
    esave = eprevious
    do imol=1,nmol
      call getref_formol(imol)
    end do
    write(ilog,*) 'Termination: Max number of energy evaluations reached'
  end if
!
! Write summary
!      
  write(ilog,*) 'Minimization Summary:'
  write(ilog,*) 'Final Energy: ', esave
  write(ilog,*) 'Total iterations:   ', min(iter,maxiter)
  write(ilog,*) 'Energy evaluations: ', min(evals,maxiter)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!  
end
!
!----------------------------------------------------------------------
!
! The two following subroutines are helpers for the minimization functions
!
!----------------------------------------------------------------------
!
! populate the gradient vector (notice that the order of degrees of freedom in
! vec_tmp in internal coordinate space is nontrivial (follows recursive structure) 
!
subroutine get_gradient(vec_tmp,grms_tmp,size_vec,tpi,stx)
!   
  use mini
  use molecule
  use forces
  use system
  use iounit
  use atoms
  use zmatrix
#ifdef ENABLE_THREADS
  use cutoffs, ONLY: molinfo
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,size_vec,stx(2)
!
  integer apos,imol,i,j,ttc
  RTYPE vec_tmp(size_vec),grms_tmp
#ifdef ENABLE_THREADS
  integer stx2(2)
  RTYPE dumm,dumm2
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, get_gradient(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif

  apos = 0
!
  if (fycxyz.eq.1) then
!
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      stx2(1) = thr_limits(61,tpi)
      stx2(2) = thr_limits(62,tpi)
    else
      stx2(1) = 1
      stx2(2) = nmol
    end if
    if (stx2(2).ge.stx2(1)) then
      apos = molinfo(stx2(1),3) - 1
    else
      apos = 0
    end if
    do imol=stx2(1),stx2(2)
#else
    do imol=1,nmol
#endif
      do ttc=1,3
        if (dc_di(imol)%frz(ttc).EQV..true.) cycle
        apos = apos + 1
        vec_tmp(apos) = dc_di(imol)%f(ttc)*mini_xyzstep
      end do
!
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) cycle
          apos = apos + 1
          if (i.le.0) then
            vec_tmp(apos) = dc_di(imol)%f(ttc)*mini_rotstep
          else
            vec_tmp(apos) = dc_di(imol)%f(ttc)*mini_intstep
          end if 
        end do
      end if
    end do
!
  else if (fycxyz.eq.2) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      stx2(1) = thr_limits(1,tpi)
      stx2(2) = thr_limits(2,tpi)
      if (stx2(2).ge.stx2(1)) then
        apos = thr_limits(101,tpi)
      else
        apos = 0
      end if
    else
      stx2(1) = 1
      stx2(2) = n
      apos = 0
    end if
    do i=stx2(1),stx2(2)
#else
    do i=1,n
#endif
      do j=1,3
        if (cart_frz(i,j).EQV..true.) cycle
        apos = apos + 1
        vec_tmp(apos) = cart_f(i,j)*mini_xyzstep
      end do
    end do
!
  end if
!
! bug check
  if ((tpi.le.0).AND.(apos.ne.size_vec)) then
    write(ilog,*) 'Fatal. Number of degrees of freedom does not ma&
 &tch up in get_gradient(...). This is most certainly a bug.'
    call fexit()
  end if
!
! calculate GRMS for convergence
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
!$OMP SINGLE
    grms_tmp = 0.0
!$OMP END SINGLE
!
    dumm2 = 1.0/(1.0*size_vec)
    dumm = dumm2*sum(vec_tmp(stx(1):stx(2))*vec_tmp(stx(1):stx(2)))
!$OMP BARRIER
!$OMP CRITICAL(MINI_GRMS_UP)
    grms_tmp = grms_tmp + dumm
!$OMP END CRITICAL(MINI_GRMS_UP)
!$OMP BARRIER
!$OMP SINGLE
    grms_tmp = sqrt(grms_tmp)
!$OMP END SINGLE
  else
#endif
  grms_tmp = sqrt(DOT_PRODUCT(vec_tmp(stx(1):stx(2)),vec_tmp(stx(1):stx(2)))/(1.0*size_vec))
#ifdef ENABLE_THREADS
  end if
#endif
!
end subroutine get_gradient
!
!----------------------------------------------------------------------
!
! propagate coordinates based on increments in minivs%dposv
!
subroutine set_position(mode,size_vec)
!
  use iounit
  use molecule
  use forces
  use system
  use zmatrix
  use mini
  use cutoffs
  use movesets
  use math
  use atoms
!
  implicit none
! 
  integer, INTENT(IN):: mode,size_vec
!
  RTYPE bhlp,rotn
  integer imol,i,ttc,apos,aone,rs
!
  aone = 1
!
  if (mode.eq.1) then
    minivs%dposv(:) = -1.0*minivs%dposv(:) ! for undoing a move
  else if (mode.eq.0) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Called set_position(...) with unsupported mode (related to minimize(...)). &
 &This is most certainly a bug.'
    call fexit()
  end if
!
  apos = 0
!
  if (fycxyz.eq.1) then
!
    bhlp = 1.0/360.0
    do imol=1,nmol
      rotn = 0.0
      do ttc=1,3
        if (dc_di(imol)%frz(ttc).EQV..true.) then
          cur_trans(ttc) = 0.0
          cycle
        end if
        apos = apos + 1
        cur_trans(ttc) = minivs%dposv(apos)
      end do
!
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
          else
            apos = apos + 1
            dc_di(imol)%incr(ttc) = minivs%dposv(apos)
!           reset to correct period: the 2nd level minimizers sometimes create bad trials
            dc_di(imol)%incr(ttc) = dc_di(imol)%incr(ttc) - 360.0*anint(dc_di(imol)%incr(ttc)*bhlp)
          end if
          if (i.le.0) then
            rotn = rotn + sin(dc_di(imol)%incr(ttc)/(2.0*RADIAN))**2
          end if
        end do
        if (rotn.gt.1.0) then
          fo_wrncnt(11) = fo_wrncnt(11) + 1
          if (fo_wrncnt(11).eq.fo_wrnlmt(11)) then
            write(ilog,*) 'Warning. Extreme increment for rigid rotation of molecule ',imol,' during minimization. If this happens&
 & often, consider trying a different minimizer.'
            write(ilog,*) 'This was warning #',fo_wrncnt(11),' of this type not all of which may be displayed.'
            if (10.0*fo_wrnlmt(11).gt.0.5*HUGE(fo_wrnlmt(11))) then
              fo_wrncnt(11) = 0
            else
              fo_wrnlmt(11) = fo_wrnlmt(11)*10
            end if
          end if
          do while (sum(sin(dc_di(imol)%incr(4:6)/(2.0*RADIAN))**2).gt.1.0)
            dc_di(imol)%incr(4:6) = dc_di(imol)%incr(4:6)/2.0
          end do
        end if
!       transfer RB rotation increment
        cur_rot(1:3) = dc_di(imol)%incr(4:6)/RADIAN
!       note that this will edit the Z-matrix
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,mode)
      end if
!
      if (mode.eq.0) then
        call makeref_formol(imol)
!        call makeref_polym(imol)
        call IMD_prealign(imol,skip_frz)
        call makexyz_formol(imol)
!
!       finally: increment coordinates
!              note rotation is done first as it relies on comm
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rotxyzm(imol,aone)
        end if
        call transxyzm(imol,aone,minivs%mshfs(:,imol))
!       due to torsional moves we need to recompute rigid-body coordinates
        call update_rigidm(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      else
        call getref_formol(imol)
        call update_rigidm(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      end if
!
    end do
!
  else if (fycxyz.eq.2) then
!
    if (mode.eq.0) then
      do imol=1,nmol
        call makeref_formol(imol)
      end do
      do i=1,n
        do ttc=1,3
          if (cart_frz(i,ttc).EQV..true.) cycle
          apos = apos + 1
          if (ttc.eq.1) x(i) = x(i) + minivs%dposv(apos)
          if (ttc.eq.2) y(i) = y(i) + minivs%dposv(apos)
          if (ttc.eq.3) z(i) = z(i) + minivs%dposv(apos)
        end do
      end do
!
!     loop over all molecules
      do imol=1,nmol 
        call update_rigidm(imol)
        call update_comv(imol)
        call genzmat(imol)
!       and shift molecule into central cell if necessary
        call update_image(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      end do
!     update pointer arrays such that analysis routines can work properly 
      call zmatfyc2()
!
    else if (mode.eq.1) then
      apos = size_vec
      do imol=1,nmol
        call getref_formol(imol)
        call update_rigidm(imol)
        call update_comv(imol)
        call genzmat(imol)
!       and shift molecule into central cell if necessary
        call update_image(imol)
!       update grid association if necessary
        if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
          do rs=rsmol(imol,1),rsmol(imol,2)
            call updateresgp(rs)
          end do
        end if
      end do
!     update pointer arrays such that analysis routines can work properly 
      call zmatfyc2()
    end if
!
  end if
!
!  bug check
  if (apos.ne.size_vec) then
    write(ilog,*) 'Fatal. Number of degrees of freedom does not ma&
 &tch up in set_position(...). This is most certainly a bug.'
    call fexit()
  end if
!    
end subroutine set_position
!
!-------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine set_position_threads(mode,size_vec,tpi,stx)
!
  use threads
  use forces
  use system
  use mini
  use atoms
  use iounit
  use zmatrix
  use molecule
  use cutoffs, ONLY: molinfo,use_cutoffs,use_mcgrid
  use movesets, ONLY: skip_frz
  use math, ONLY: RADIAN
!
  implicit none
!
  integer, INTENT(IN):: tpi,stx(2),size_vec,mode
!
  logical atrue,jobflags(4)
  integer apos,imol,i,j,ttc,azero,atwo,ixx,rs
  RTYPE bhlp,rotn
  integer(KIND=8) ttimer
!
!$OMP BARRIER
  atrue = .true.
  azero = 0
  atwo = 2
!
  if (mode.eq.1) then
    minivs%dposv(stx(1):stx(2)) = -1.0*minivs%dposv(stx(1):stx(2)) ! for undoing a move
  else if (mode.eq.0) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Called set_position_threads(...) with unsupported mode (related to minimize(...)). &
 &This is most certainly a bug.'
    call fexit()
  end if
!
  if (fycxyz.eq.1) then
!
    bhlp = 1.0/360.0
    if (mode.eq.0) then
!     coordinate backup
      xref(thr_limits(1,tpi):thr_limits(2,tpi)) = x(thr_limits(1,tpi):thr_limits(2,tpi))
      yref(thr_limits(1,tpi):thr_limits(2,tpi)) = y(thr_limits(1,tpi):thr_limits(2,tpi))
      zref(thr_limits(1,tpi):thr_limits(2,tpi)) = z(thr_limits(1,tpi):thr_limits(2,tpi))
    else if (mode.eq.1) then
!     coordinate recovery
      x(thr_limits(1,tpi):thr_limits(2,tpi)) = xref(thr_limits(1,tpi):thr_limits(2,tpi))
      y(thr_limits(1,tpi):thr_limits(2,tpi)) = yref(thr_limits(1,tpi):thr_limits(2,tpi))
      z(thr_limits(1,tpi):thr_limits(2,tpi)) = zref(thr_limits(1,tpi):thr_limits(2,tpi))
    end if
!$OMP BARRIER
!
! first loop over all molecules that are treated by one thread only
    if (thr_dlb(12,1).gt.0) then
      if (tpi.eq.1) thr_dlb(12,2) = thr_dlb(12,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(27,tpi) = thr_timings(27,tpi) + ttimer
    end if
!
    if (thr_limits(62,tpi).ge.thr_limits(61,tpi)) then
      apos = molinfo(thr_limits(61,tpi),3) - 1
    else
      apos = 0
    end if
    do imol=thr_limits(61,tpi),thr_limits(62,tpi)
      rotn = 0.0
!     first: center-of-mass translation (linear motion)
      do j=1,3
        if (dc_di(imol)%frz(j).EQV..true.) then
          dc_di(imol)%incr(j) = 0.0
          cycle
        end if
        apos = apos + 1
        dc_di(imol)%incr(j) = minivs%dposv(apos)
      end do
!     now all rotational motion (including RB) 
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            ttc = i+6 ! RB rotation
          else
            ttc =dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (ttc.le.0) cycle
          if (dc_di(imol)%frz(ttc).EQV..true.) then
            dc_di(imol)%incr(ttc) = 0.0
            cycle
          end if
          apos = apos + 1
          dc_di(imol)%incr(ttc) = minivs%dposv(apos)
!         reset to correct period: the 2nd level minimizers sometimes create bad trials
          dc_di(imol)%incr(ttc) = dc_di(imol)%incr(ttc) - 360.0*anint(dc_di(imol)%incr(ttc)*bhlp)
          if (i.le.0) then
            rotn = rotn + sin(dc_di(imol)%incr(ttc)/(2.0*RADIAN))**2
          end if
        end do
        if (rotn.gt.1.0) then
          fo_wrncnt(11) = fo_wrncnt(11) + 1
          if (fo_wrncnt(11).eq.fo_wrnlmt(11)) then
            write(ilog,*) 'Warning. Extreme increment for rigid rotation of molecule ',imol,' during minimization. If this happens&
 & often, consider trying a different minimizer.'
            write(ilog,*) 'This was warning #',fo_wrncnt(11),' of this type not all of which may be displayed.'
            if (10.0*fo_wrnlmt(11).gt.0.5*HUGE(fo_wrnlmt(11))) then
              fo_wrncnt(11) = 0
            else
              fo_wrnlmt(11) = fo_wrnlmt(11)*10
            end if
          end if
          do while (sum(sin(dc_di(imol)%incr(4:6)/(2.0*RADIAN))**2).gt.1.0)
            dc_di(imol)%incr(4:6) = dc_di(imol)%incr(4:6)/2.0
          end do
        end if
!       now transfer
        if (dc_di(imol)%maxntor.gt.0) call tmdmove(imol,azero)
      end if
      if (thr_mlgix(imol).gt.0) cycle 
!
      if (mode.eq.0) then
        if (dc_di(imol)%maxntor.gt.0) then
          call IMD_prealign(imol,skip_frz)
          call makexyz_formol(imol)
        end if
        call rbcxyzm(imol,azero,atrue)
      end if
!     due to torsional moves we need to recompute rigid-body coordinates
      call update_rigidm(imol)
    end do
    if (thr_dlb(12,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(28,tpi) = thr_timings(28,tpi) + ttimer
    end if
!
!   now the same for any internally parallelized molecules
    if (nmlgs.gt.0) then
      jobflags(:) = .false.
      jobflags(1) = .true.
    end if
!$OMP BARRIER
    do ixx=1,nmlgs
      imol = mlg_limits(ixx,5,tpi)
      if (mode.eq.0) then
        if (dc_di(imol)%maxntor.gt.0) then
!$OMP SINGLE
          call IMD_prealign(imol,skip_frz)
!$OMP END SINGLE
          call makexyz_threads(ixx,tpi) ! ends with a barrier
        end if
!
        call rbcxyzm(ixx,tpi,atrue) ! ends with an implied barrier
      end if
!     due to torsional moves we need to recompute rigid-body coordinates
      call molops_threads(ixx,tpi,jobflags)
    end do
!
    if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
!$OMP BARRIER
!     upgrade grid association if needed 
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call updateresgp(rs)
      end do
    end if
!
  else if (fycxyz.eq.2) then
!
    if (mode.eq.0) then
!     coordinate backup
      xref(thr_limits(1,tpi):thr_limits(2,tpi)) = x(thr_limits(1,tpi):thr_limits(2,tpi))
      yref(thr_limits(1,tpi):thr_limits(2,tpi)) = y(thr_limits(1,tpi):thr_limits(2,tpi))
      zref(thr_limits(1,tpi):thr_limits(2,tpi)) = z(thr_limits(1,tpi):thr_limits(2,tpi))
!$OMP BARRIER
      if (thr_limits(2,tpi).ge.thr_limits(1,tpi)) then
        apos = thr_limits(101,tpi)
      end if
!     now loop over all molecules
      do i=thr_limits(1,tpi),thr_limits(2,tpi)
        do j=1,3
          if (cart_frz(i,j).EQV..true.) cycle
          apos = apos + 1
          if (j.eq.1) x(i) = x(i) + minivs%dposv(apos)
          if (j.eq.2) y(i) = y(i) + minivs%dposv(apos)
          if (j.eq.3) z(i) = z(i) + minivs%dposv(apos)
        end do
      end do
    else if (mode.eq.1) then
!     coordinate recovery
      x(thr_limits(1,tpi):thr_limits(2,tpi)) = xref(thr_limits(1,tpi):thr_limits(2,tpi))
      y(thr_limits(1,tpi):thr_limits(2,tpi)) = yref(thr_limits(1,tpi):thr_limits(2,tpi))
      z(thr_limits(1,tpi):thr_limits(2,tpi)) = zref(thr_limits(1,tpi):thr_limits(2,tpi))
    end if
!$OMP BARRIER
    if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
!     upgrade grid association if needed 
      do rs=thr_limits(3,tpi),thr_limits(4,tpi)
        call updateresgp(rs)
      end do
    end if
!   loop over all molecules again
    if (thr_dlb(10,1).gt.0) then
      if (tpi.eq.1) thr_dlb(10,2) = thr_dlb(10,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(23,tpi) = thr_timings(23,tpi) + ttimer
    end if
    do imol=thr_limits(35,tpi),thr_limits(36,tpi)
      if (thr_mlgix(imol).gt.0) cycle
      call update_rigidm(imol)
      call update_comv(imol)
      call genzmat(imol)
      call update_image(imol)
    end do
    if (thr_dlb(10,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(24,tpi) = thr_timings(24,tpi) + ttimer
    end if
    jobflags(:) = .true.
    do i=1,nmlgs
!     covers update_rigidm, update_comv, genzmat, update_image
      call molops_threads(i,tpi,jobflags)
    end do
!   update pointer arrays such that analysis routines can work properly
!$OMP BARRIER
    call zmatfyc_threads(tpi,atwo)
!
  end if
!
!$OMP BARRIER
!
end subroutine set_position_threads
!
#endif
!
!----------------------------------------------------------------------
!
! this routine takes advantage of the built-in dynamics routines while externally modifying 
! the target bath temperature 
!
subroutine stochastic_min(maxiter,tpi)
!  
  use iounit
  use molecule
  use forces
  use system
  use mcsums
  use mini
  use energies
  use cutoffs
  use atoms
  use shakeetal
  use mini
  use zmatrix
#ifdef ENABLE_THREADS
  use threads
#endif
!    
  implicit none
!
  integer, INTENT(IN):: tpi,maxiter
!
  integer i,istep,collectend,ndump,nstepi,nrdof,schkbuf,imol,j,dumm,stx(2)
  RTYPE tmax,targetrate,instrate,grms,buke,butau
  integer(KIND=8) t1,t2
  logical postanalyze,afalse,atrue
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, stochastic_min(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  dumm = 100
  atrue = .true.
  afalse = .false.
  grms = 0.0
!    
! take a few steps of steepest descent to relax major energy barriers
  nstepi = maxiter
  if (mini_sc_sdsteps.gt.0) then
    if ((no_shake.EQV..false.).AND.(fycxyz.eq.2)) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      write(ilog,*) 'Warning. Skipping requested steepest-descent steps due to presence of holonomic constraints for &
 &stochastic minimizer.'
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
    else
      j = 1
      call minimize(mini_sc_sdsteps,minivs%grms,mini_stepsize,j,tpi)
    end if
  end if
!
! "vectorize" all degrees of freedom into a single vector with dimension "nrdof"
  nrdof = 0
  if (fycxyz.eq.1) then
    do imol=1,nmol
      if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        write(ilog,*) 'Fatal. Two-atom molecules are not yet supported in minimization runs.'
        call fexit()
      end if
      if (tpi.le.1) molinfo(imol,3) = nrdof + 1
      do j=1,3
        if (dc_di(imol)%frz(j).EQV..true.) cycle
        nrdof = nrdof + 1
      end do
      if ((atmol(imol,2)-atmol(imol,1)).gt.1) then ! support for diatomic or linear mol.s missing
        do i=-2,dc_di(imol)%maxntor
          if (i.le.0) then
            j = i+6 ! RB rotation
          else
            j = dc_di(imol)%recurs(izrot(dc_di(imol)%recurs(i,1))%treevs(4),3)
          end if
          if (j.le.0) cycle
          if (dc_di(imol)%frz(j).EQV..true.) cycle
          nrdof = nrdof + 1
        end do
      end if
    end do 
  else if (fycxyz.eq.2) then
! SHAKE-type constraints and drift removal must be handled separately (not mappable)
    do imol=1,nmol
      if (tpi.le.1) molinfo(imol,3) = nrdof + 1
      do i=atmol(imol,1),atmol(imol,2)
        do j=1,3
          if (cart_frz(i,j).EQV..true.) cycle
          nrdof = nrdof + 1
        end do
      end do
    end do
  else
    write(ilog,*) 'Fatal. Called minimize(...) with unsupported system representation. This is an omission bug.'
    call fexit()
  end if
!
  if (nrdof.le.0) then
    write(ilog,*) 'Fatal. No remaining degrees of freedom to minimize. Relieve constraints.'
    call fexit()
  end if
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    call threads_bounds(1,nrdof,tpi,thrdat%maxn,stx)
  else
    stx(1) = 1
    stx(2) = nrdof
  end if
#else
  stx(1) = 1
  stx(2) = nrdof
#endif
!
  if (tpi.le.1) allocate(minivs%cgrv(nrdof))
!
  if (mini_sc_tbath.gt.kelvin) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
    write(ilog,*) 'Warning. For the stochastic minimizer, the supplied target temperature is larger than the &
 &initial one. This run will most likely terminate immediately after the number of steps defined by FMCSC_MINI_SC_HEAT.'
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
  end if
  buke = kelvin
  butau = tstat%params(1)
  if (tpi.le.1) kelvin = max(kelvin,mini_sc_tbath) ! target temperature (T of bath)
  collectend=floor(nstepi*mini_sc_heat)
  tmax=0.0
  ndump = 0
!    
  istep = 1
  if (tpi.le.1) nstep = 0
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if (fycxyz.eq.2) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      call cart_mdmove_threads(istep,dumm,atrue,tpi)
    else
#endif
    call cart_mdmove(istep,dumm,atrue)
#ifdef ENABLE_THREADS
    end if
#endif
  else if (fycxyz.eq.1) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      call int_mdmove_threads(istep,dumm,atrue,tpi)
    else
#endif
    call int_mdmove(istep,dumm,atrue)
#ifdef ENABLE_THREADS
    end if
#endif
  else
    write(ilog,*) 'Fatal. Called stochastic_min(...) with unsupported system representation. This is an omission bug.'
    call fexit()
  end if
  do istep=2,nstepi
    schkbuf = nsancheck 
    if (tpi.le.1) nsancheck = nstepi+1 ! suppress ensemble output from MD
    if (fycxyz.eq.2) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call cart_mdmove_threads(istep,dumm,afalse,tpi) ! ends with barrier
      else
#endif
      call cart_mdmove(istep,dumm,afalse)
#ifdef ENABLE_THREADS
      end if
#endif
    else if (fycxyz.eq.1) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call int_mdmove_threads(istep,dumm,afalse,tpi) ! ends with barrier
      else
#endif
      call int_mdmove(istep,dumm,afalse)
#ifdef ENABLE_THREADS
      end if
#endif
    end if
    if (tpi.le.1) nsancheck = schkbuf
    call get_gradient(minivs%cgrv,minivs%grms,nrdof,tpi,stx) ! has barrier
    if (tpi.le.1) then
      if(mod(istep,nsancheck).eq.0) then
 12  format('Step: ',i8,' Current Energy: ',g14.6,' Current Temperature: ',g14.6, ' GRMS: ',g14.6)
        write(ilog,12) istep,ens%insU,ens%insT,minivs%grms
      end if
    end if
!
    if(istep.lt.collectend) then
      if(tmax.lt.ens%insT) then  !T is initialized to FMCSC_TEMP, but it will rise as
        tmax = ens%insT          !PE is converted to KE.
      end if
    end if 
!      
! Monitor gradients here and terminate when criteria reached.
    if (istep.gt.collectend) then
      if (minivs%grms.lt.mini_econv) then
 14 format('Normalized RMS gradient is ',g10.4,' and has dropped below convergence criterion&
 & of ',g10.4,'.')
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        write(ilog,14) minivs%grms,mini_econv
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
        postanalyze = .true.
        exit
      else if(ens%insT.lt.mini_sc_tbath) then
 13 format('System temperature is ',g10.4,' and has dropped below ',g10.4,'. Terminating.')
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        write(ilog,13) ens%insT,mini_sc_tbath
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
        postanalyze = .true.
        exit
      else if(istep .eq. nstepi) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
        write(ilog,*) 'Maximum number of stochastic minimization steps exceeded. Forced termination.'
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
        postanalyze = .true.
        exit
      end if
    end if
!      
! Set target rate based on number of steps given
    if(istep.eq.collectend) then
      targetrate=1.25*(tmax-mini_sc_tbath)/(nstepi-collectend)
    end if
!      
! Dynamically update target temperature and rate for faster end convergence (should do after heating MD steps)
    if ((tpi.le.1).AND.(istep.gt.collectend)) then
      kelvin = mini_sc_tbath + buke*(exp(-((istep-collectend)/(0.2*(nstepi-collectend)))**2)) - buke*1.389e-11
      instrate=(tmax - ens%insT)/(istep-collectend)
      if ((instrate.gt.0.0).AND.(targetrate.gt.0.0)) then
        tstat%params(1) = tstat%params(1)*(instrate/targetrate)
        if (tstat%params(1).lt.(10.*dyn_dt)) then
          tstat%params(1) = 10.0*dyn_dt
        end if
        if (tstat%params(1).gt.butau) then
          tstat%params(1) = butau
        end if
      end if
    end if
!
    if (tpi.le.1) call System_Clock(t1)
    call mcstat(istep,ndump,tpi)
    if (tpi.le.1) then 
      call System_Clock(t2)
      time_analysis = time_analysis + 1.0*(t2 - t1)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!     
  end do
!    
! last-step analysis if converged
  if (postanalyze.EQV..true.) then
    if (tpi.le.1) call System_Clock(t1)
    j = nsim
    call mcstat(j,ndump,tpi)
    if (tpi.le.1) then
      call System_Clock(t2)
      time_analysis = time_analysis + 1.0*(t2 - t1)
    end if
  end if
!   
! Write summary
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  deallocate(minivs%cgrv)  
  write(ilog,*) 'Minimization Summary:'
  write(ilog,*) 'Final Energy: ', esave
  write(ilog,*) 'Total iterations: ', min(istep,nstepi)
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
!    
end
!
!---------------------------------------------------------------------------------------------------------
!
