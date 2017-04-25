c
#include "macros.i"
c
c     #######################################################
c     ##                                                   ##
c     ## This subroutine propagates internal degrees of    ##
c     ## freedom by Brownian (over-damped) Dynamics        ##
c     ##                                                   ##
c     #######################################################
c     
c
      subroutine int_bdmove(istep,ndump)
c
      use iounit
      use math
      use atoms
      use molecule
      use forces
      use energies
      use movesets
      use system
      use units
      use mcsums
      use movesets
      use fyoc
      use torsn
      use cutoffs
c
      implicit none
c
      integer imol,ii,j,istep,ndump,aone,nt,ttc,rs,i,boxdof,jj,atm,ints
      integer atmax,intmax,freeunit,k
      RTYPE energy3,ttmp,t1,t2,fr1,fr2,rand,rand1,rand2,random,t3,t4
      RTYPE tsc,force3,normal,boxovol,radius,vel(3),speed,nmm
      RTYPE, ALLOCATABLE:: pn(:,:),g_rand(:,:),fn(:),qn(:)
      RTYPE, ALLOCATABLE:: bma_ref(:,:,:),pos_ref(:,:,:)
      logical afalse
c
      afalse = .false.
      aone = 1
      nstep = nstep + 1
      intmax = 0
      atmax = 0
c
      do imol=1,nmol
c       sanity check in first step against (partially) unsupported 2-atom molecules
        if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          write(ilog,*) 'Fatal. Two-atom molecules are not yet support
     &ed in dynamics. Check back later.'
          call fexit()
        end if
        call update_rigidm(imol) 
c
c       find indices to allocate matrices to the largest molecule
        atm = 3*(atmol(imol,2)-atmol(imol,1)+1)
        if (atmol(imol,2).eq.atmol(imol,1)) then
          ints = ntormol(moltypid(imol))+3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ints = ntormol(moltypid(imol))+5
        else
          ints = ntormol(moltypid(imol))+6
        end if
        if (atm > atmax) then
          atmax = atm
        end if
        if (ints > intmax) then
          intmax = ints
        end if
      end do
c
c     allocate cross molecule matrices
      allocate(g_rand(atmax,nmol))
      allocate(pn(intmax,nmol))
      allocate(bma_ref(atmax,intmax,nmol))
      allocate(pos_ref(3,atmax/3,nmol))
c
c     get energies
      call CPU_time(t1)
      esave = force3(esterms,esterms_tr,esterms_lr,afalse)
      ens%insU = esave
      boxovol = ens%insV
      call CPU_time(t2)
      time_energy = time_energy + t2 - t1      
c
c     calculate internal forces
      call cart2int_f()
c
c     get the temperature re-scaling factor from thermostat
      if ((ens%flag.eq.1).OR.(ens%flag.eq.3)) then
        call thermostat(tsc)
      else
        tsc = 1.0
      end if
c
c     loop over each molecule
      do imol=1,nmol
c
        atm = 3*(atmol(imol,2)-atmol(imol,1)+1)
        if (atmol(imol,2).eq.atmol(imol,1)) then
          ints = ntormol(moltypid(imol))+3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ints = ntormol(moltypid(imol))+5
        else
          ints = ntormol(moltypid(imol))+6
        end if
c
        allocate(fn(ints))
        allocate(bma(atm,ints))
c
c       used to calculate velocities
        j = 0
        do i=atmol(imol,1),atmol(imol,2)
          j = j+1 
          pos_ref(1,j,imol) = x(i)
          pos_ref(2,j,imol) = y(i) 
          pos_ref(3,j,imol) = z(i)
        end do
c
c       calculate the base vectors and forces
        call CPU_time(t3)
        call metric_tensors_b(imol)
        call CPU_time(t4)
        time_b = time_b + t4 - t3
c
c       back up vector b(n) at step q(n) and generate random gaussian numbers
        do j=1,atm/3
          write(*,*) 'radius atom',j,atr(j)
          do k=1,3
            nmm = normal()
            g_rand(3*(j-1)+k,imol) =nmm*sqrt(12.0*u_dyn_kb*kelvin*PI
     &*fric_ga*atr(j)*dyn_dt)
            do i=1,ints
              bma_ref(3*(j-1)+k,i,imol) = bma(3*(j-1)+k,i)
            end do
          end do
        end do
c
c       increment position t to intermediate step: t + dt*
        do j=1,ints
          rand = 0.0
          do jj=1,atm
            rand = rand + bma(jj,j)*g_rand(jj,imol)
          end do
          pn(j,imol) =fn(j)*dyn_dt + rand!*sqrt(2*u_dyn_kb*kelvin*dyn_dt)
          write(bds,*) j,rand
        end do
c
c       gather the displacements in internal coordinate space, update velocities
c       first: center-of-mass translation
        do j=1,3
          if (dc_di(imol)%frz(j).EQV..true.) then
            cur_trans(j) = 0.0
            cycle
          else
            cur_trans(j) = pn(j,imol)
          end if
        end do
c
c       second: rigid-body rotation
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
c         do nothing
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        else
          do j=4,6
            if (dc_di(imol)%frz(j).EQV..true.) then
              cur_rot(j-3) = 0.0
              cycle
            else
              cur_rot(j-3) = pn(j,imol)/RADIAN
            end if
          end do
        end if
c    
 57     format('Mol. ',i6,' Res. ',i8,' #',i2,' :',g14.7)
 58     format('Mol. ',i6,' Res. ',i8,' :',g14.7)
c
c       third: torsional evolution
        do rs=rsmol(imol,1),rsmol(imol,2)
          if (wline(rs).gt.0) then
            ttc = wnr(rs)
            ttmp = omega(rs) + pn(ttc,imol)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = omega(rs)
            end if
            cur_omega = ttmp
            if ((cur_omega.gt.180.0).OR.(cur_omega.lt.-180.0)) then
              write(ilog,*) 'Velocity exception in Omega-Tors.:'
              write(ilog,58) imol,rs,pn(ttc,imol)
            end if
            call setw(rs,cur_omega,aone)
          end if
          cur_phi = phi(rs)
          if (fline(rs).gt.0) then
            ttc = fnr(rs)
            ttmp = phi(rs) + pn(ttc,imol)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = phi(rs)
            end if
            cur_phi = ttmp
          end if
          cur_psi = psi(rs)
          if (yline(rs).gt.0) then
            ttc = ynr(rs)
            ttmp = psi(rs) + pn(ttc,imol)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = psi(rs)
            end if
            cur_psi = ttmp
          end if
          if ((cur_phi.gt.180.0).OR.(cur_phi.lt.-180.0)) then
            write(ilog,*) 'Velocity exception in Psi-Tors.:'
            write(ilog,58) imol,rs,pn(ttc,imol)
          end if
          if ((cur_psi.gt.180.0).OR.(cur_psi.lt.-180.0)) then
            write(ilog,*) 'Velocity exception in Psi-Tors.:'
            write(ilog,58) imol,rs,pn(ttc,imol)
          end if
          if ((fline(rs).gt.0).OR.(yline(rs).gt.0)) then
            call setfy(rs,cur_phi,cur_psi,0)
          end if
          do nt=1,nnucs(rs)
            ttc = nucsnr(nt,rs)
            ttmp = nucs(nt,rs) + pn(ttc,imol)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = nucs(nt,rs)
            end if
            cur_nucs(nt) = ttmp
            if ((cur_nucs(nt).gt.180).OR.(cur_nucs(nt).lt.-180)) then
              write(ilog,*) 'Velocity exception in Nuc. Tors.:'
              write(ilog,57) imol,rs,nt,pn(ttc,imol)
            end if
            cur_nucflag(nt) = .true.
          end do
          if (nnucs(rs).gt.0) then
            call setnucs(rs,cur_nucs,aone)
          end if
          do nt=1,nchi(rs)
            ttc = chinr(nt,rs)
            ttmp = chi(nt,rs) + pn(ttc,imol)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = chi(nt,rs)
            end if
            cur_chis(nt) = ttmp
            if ((cur_chis(nt).gt.180).OR.(cur_chis(nt).lt.-180)) then
              write(ilog,*) 'Velocity exception in Chi-Tors.:'
              write(ilog,57) imol,rs,nt,pn(ttc,imol)
            end if
            cur_chiflag(nt) = .true.
          end do
          if (nchi(rs).gt.0) then
            call setchi(rs,cur_chis,aone)
          end if
        end do
        call makeref_formol(imol)
        call makexyz_formol(imol)
c    
c       finally: increment coordinates
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rotxyzm(imol,aone)
        end if
        call transxyzm(imol,aone)
c       due to torsional moves we need to recompute rigid-body coordinates
        call update_rigidm(imol)
c
        deallocate(fn)
        deallocate(bma)
c
        end do
c
c       re-calculate internal forces
        call cart2int_f()
c
c       re-loop over each molcule
        do imol=1,nmol
c
c       indices
        atm = 3*(atmol(imol,2)-atmol(imol,1)+1)
        if (atmol(imol,2).eq.atmol(imol,1)) then
          ints = ntormol(moltypid(imol))+3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ints = ntormol(moltypid(imol))+5
        else
          ints = ntormol(moltypid(imol))+6
        end if
c
        allocate(qn(ints))
        allocate(fn(ints))
        allocate(bma(atm,ints))
c
c       re-calculate the base vector
        call metric_tensors_b(imol)
c
c       for the second step: increment position t to final position: t + dt
        do j=1,ints
          do jj=1,atm
            rand1 = rand1 + bma(jj,j)*g_rand(jj,imol)
            rand2 = rand2 + bma_ref(jj,j,imol)*g_rand(jj,imol)
          end do
          qn(j) = .5*(rand2-rand1)
        end do
        write(*,*) 'qn',qn(:)
c
c       re-gather the displacements in internal coordinate space
c       first: center-of-mass translation
        do j=1,3
          if (dc_di(imol)%frz(j).EQV..true.) then
            cur_trans(j) = 0.0
            cycle
          else
            cur_trans(j) = qn(j)
          end if
        end do
c
c       second: rigid-body rotation
        if ((atmol(imol,2)-atmol(imol,1)).eq.0) then
c         do nothing
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
        else
          do j=4,6
            if (dc_di(imol)%frz(j).EQV..true.) then
              cur_rot(j-3) = 0.0
              cycle
            else
              cur_rot(j-3) = qn(j)/RADIAN
            end if
          end do
        end if
c    
c       third: torsional evolution
        do rs=rsmol(imol,1),rsmol(imol,2)
          if (wline(rs).gt.0) then
            ttc = wnr(rs)
            ttmp = omega(rs) + qn(ttc)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = omega(rs)
            end if
            cur_omega = ttmp
            if ((cur_omega.gt.180.0).OR.(cur_omega.lt.-180.0)) then
              write(ilog,*) 'Velocity exception in Omega-Tors.:'
              write(ilog,58) imol,rs,qn(ttc)
            end if
            call setw(rs,cur_omega,aone)
          end if
          cur_phi = phi(rs)
          if (fline(rs).gt.0) then
            ttc = fnr(rs)
            ttmp = phi(rs) + qn(ttc)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = phi(rs)
            end if
            cur_phi = ttmp
          end if
          cur_psi = psi(rs)
          if (yline(rs).gt.0) then
            ttc = ynr(rs)
            ttmp = psi(rs) + qn(ttc)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = psi(rs)
            end if
            cur_psi = ttmp
          end if
          if ((cur_phi.gt.180.0).OR.(cur_phi.lt.-180.0)) then
            write(ilog,*) 'Velocity exception in Phi-Tors.:'
            write(ilog,58) imol,rs,qn(ttc)
          end if
          if ((cur_psi.gt.180.0).OR.(cur_psi.lt.-180.0)) then
            write(ilog,*) 'Velocity exception in Psi-Tors.:'
            write(ilog,58) imol,rs,qn(ttc)
          end if
          if ((fline(rs).gt.0).OR.(yline(rs).gt.0)) then
            call setfy(rs,cur_phi,cur_psi,0)
          end if
          do nt=1,nnucs(rs)
            ttc = nucsnr(nt,rs)
            ttmp = nucs(nt,rs) + qn(ttc)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = nucs(nt,rs)
            end if
            cur_nucs(nt) = ttmp
            if ((cur_nucs(nt).gt.180).OR.(cur_nucs(nt).lt.-180)) then
              write(ilog,*) 'Velocity exception in Nuc. Tors.:'
              write(ilog,57) imol,rs,nt,qn(ttc)
            end if
            cur_nucflag(nt) = .true.
          end do
          if (nnucs(rs).gt.0) then
            call setnucs(rs,cur_nucs,aone)
          end if
          do nt=1,nchi(rs)
            ttc = chinr(nt,rs)
            ttmp = chi(nt,rs) + qn(ttc)
            if (ttmp.gt.180.0) ttmp = ttmp - 360.0
            if (ttmp.lt.-180.0) ttmp = ttmp + 360.0
            if (dc_di(imol)%frz(ttc).EQV..true.) then
              ttmp = chi(nt,rs)
            end if
            cur_chis(nt) = ttmp
            if ((cur_chis(nt).gt.180).OR.(cur_chis(nt).lt.-180)) then
              write(ilog,*) 'Velocity exception in Chi-Tors.:'
              write(ilog,57) imol,rs,nt,qn(ttc)
            end if
            cur_chiflag(nt) = .true.
          end do
          if (nchi(rs).gt.0) then
            call setchi(rs,cur_chis,aone)
          end if
        end do
        call makeref_formol(imol)
        call makexyz_formol(imol)
c    
c       finally: increment coordinates
        if ((atmol(imol,2)-atmol(imol,1)).gt.0) then
          call rotxyzm(imol,aone)
        end if
        call transxyzm(imol,aone)
c       due to torsional moves we need to recompute rigid-body coordinates
        call update_rigidm(imol)
c
c     output velocities
      speed = 0.0
      j=0
      do i=atmol(imol,1),atmol(imol,2)
        write(*,*) 'i',i
        j=j+1
        write(*,*) 'j',j
        write(*,*) dyn_dt
        write(*,*) pos_ref(:,j,imol)
        write(*,*) x(i),y(i),z(i)
        vel(1) = (pos_ref(1,j,imol)-x(i))/dyn_dt
        vel(2) = (pos_ref(2,j,imol)-y(i))/dyn_dt
        vel(3) = (pos_ref(3,j,imol)-z(i))/dyn_dt
        write(*,*) vel(:)
        speed = speed + sqrt(vel(1)**2+vel(2)**2+vel(3)**2)
        write(*,*) speed
      end do
      write(vbd,*) imol,speed
c      write(tif,*) nstep*dyn_dt,time_energy,time_b
c
      deallocate(qn)
      deallocate(fn)
      deallocate(bma)
c
      end do
c
      call get_ensv(boxovol)
      call drift_removal(fr1,fr2)
c
 997  format(' Kinetic E   | Potential E | Total E     | Te
     &mperature | Drift-xyz   | Drift-Euler |')
 998  format(' Kinetic E   | Potential E | Enthalpy    | Te
     &mperature | Pressure    | Box volume  |')
 999  format('---------------------------------------------
     &---------------------------------------')
      if (istep.eq.1) then
c       write header line
        if (nsim.ge.nsancheck) then
          if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
            write(ilog,997)
c
            write(ilog,999)
          else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
            write(ilog,998)
            write(ilog,999)
          end if
        end if
      end if
      if (mod(nstep,nsancheck).eq.0) then
        if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
          write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU,ens%insT,
     &fr1,fr2
        else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
          write(ilog,56) ens%insK,ens%insU,ens%insK + ens%insU + bnd_pV,
     &ens%insT,ens%insP,boxovol/1.0e6
         end if
      end if
c
   56 format(50(g13.6,1x))
c
      deallocate(bma_ref)
      deallocate(pn)
      deallocate(g_rand)
      deallocate(pos_ref)
c
      CONTAINS
c
c ------------------------------------------------------------------------------
c
        subroutine metric_tensors_b(imol)
c
        use iounit
        use forces
        use atoms
        use molecule
        use sequen
        use math
c 
        implicit none
c 
        integer imol,i,j,k,atm,freeunit,ints
        RTYPE, ALLOCATABLE:: pep(:,:),pei(:,:),iii(:,:),ama(:,:)
        RTYPE tem,t2,t1
c 
        atm = 3*(atmol(imol,2)-atmol(imol,1)+1)
        if (atmol(imol,2).eq.atmol(imol,1)) then
          ints = ntormol(moltypid(imol))+3
        else if ((atmol(imol,2)-atmol(imol,1)).eq.1) then
          ints = ntormol(moltypid(imol))+5
        else
          ints = ntormol(moltypid(imol))+6
        end if
c 
        allocate(ama(ints,atm))
        allocate(pep(ints,ints))
        allocate(pei(ints,ints))
        allocate(iii(ints,ints))
c 
        ama = dc_di(imol)%it
c
c       multiply the viscosity term into the aij matrix (dR_i*epsilon_i/dq)
        do i=1,ints
          do j=1,atm/3
            do k=1,3
              ama(i,3*(j-1)+k)=ama(i,3*(j-1)+k)*(6.0*PI*fric_ga*atr(j))
            end do
          end do
        end do
c
c       generate the covariant metric tensor
        do i=1,ints
          do j=1,ints
            pep(i,j) = 0.0
            pei(i,j) = 0.0
            iii(i,j) = 0.0
            do k=1,atm
              pep(i,j) = pep(i,j)+ama(i,k)*ama(j,k)
            end do
          end do
        end do
c 
c       get the lower triangular matrix (L) through Cholesky decomposition
        do i=1,ints
          tem = 0.0d0
          do k=1,(i-1)
            tem = tem + pei(k,i)*pei(k,i)
          end do
          tem = pep(i,i) - tem
          if (tem .le.0.0d0) then
            write(ilog,*) 'negative value in sqrt for i=',i,tem
            call fexit()
          end if
          pei(i,i) = sqrt(tem)
          if (abs(pei(i,i)-0.0d0) .lt. 1.0d-8) then
            write(ilog,*) 'zero entry in diagonal of Cholesky matrix',
     &           abs(pei(i,i)-0.0d0)
          end if
c 
          do j=(i+1),ints
            tem = 0.0d0
            do k=1,(i-1)
              tem = tem +  pei(k,i)*pei(k,j)
            end do
            pei(i,j) = (pep(j,i) - tem)/pei(i,i)
          end do
        end do
c 
c       get the inverse of L:
c       note we assume the lower triangular nature of the inverse
c       otherwise straightforward recursive back calculation
        do i=1,ints
          iii(i,i) = 1.0/pei(i,i)
          do j=i-1,1,-1
            tem = 0.0
            do k=1,i
              tem = tem + pei(k,i)*iii(k,j)
            end do
            iii(i,j) = -tem/pei(i,i)
          end do
        end do
c 
c       multiply L with its transpose to get the inverse of the covariant metric tensor
c       (i.e., the contravariant metric tensor, i.e., the inverse of pep)
        do i=1,ints
          do j=1,ints
            pei(j,i) = 0.0
            do k=1,ints
              pei(j,i) = pei(j,i)+iii(k,i)*iii(k,j)
            end do
          end do
        end do
c 
c       finally multiply the contravariant tensor with the transposed covariant base vectors
c       to get the contravariant base vectors
        do j=1,atm
          do i=1,ints
            bma(j,i) = 0.0
            do k=1,ints
              bma(j,i) = bma(j,i) + pei(k,i)*ama(k,j)!dc_di(imol)%it(k,j)
            end do
          end do
        end do
c 
c       now use the contravariant base vectors to get effective forces
        do i=1,ints
          fn(i) = 0.0
          do j=1,ints
            fn(i) = fn(i) + dc_di(imol)%f(j)*pei(i,j)
          end do
        end do
c 
        deallocate(ama)
        deallocate(pep)
        deallocate(pei)
        deallocate(iii)
c
        end subroutine metric_tensors_b
c
      end subroutine int_bdmove
