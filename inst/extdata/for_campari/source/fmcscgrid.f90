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
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!--------------------------------------------------------------------
!
subroutine setupmcgrid()
!
  use iounit
  use sequen
  use atoms
  use cutoffs
  use system
  use mcgrid
!
  implicit none
!
  integer i,j,which,gridnavi
  RTYPE pvd2,cutplusmargin2
  RTYPE refp(3)
!
  if (use_waterloops.EQV..true.) then
!    call map_molecules()
  end if
!
! first use boundary condition to determine grid origin and deltas
! periodic BC (PBC)
  if (n_topgrps.gt.0) then
    cutplusmargin2 = (mcel_cutoff + 2.0*topgrprad)**2.0
  else
    cutplusmargin2 = (mcel_cutoff + 2.0*mrrd)**2.0
  end if
  if (bnd_type.eq.1) then
!   rectangular box
    if (bnd_shape.eq.1) then
      do i=1,3
        grid%origin(i) = bnd_params(3+i)
        grid%deltas(i) = bnd_params(i)/(1.0*grid%dim(i))
      end do
!   cylinder
    else if (bnd_shape.eq.3) then
      grid%origin(1:2) = bnd_params(1:2) - 1.5*bnd_params(4)
      grid%deltas(1:2) = 3.0*bnd_params(4)/(1.0*grid%dim(1:2))
      grid%origin(3) = bnd_params(3) - 0.5*bnd_params(6)
      grid%deltas(3) = bnd_params(6)/(1.0*grid%dim(3))
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in set&
 &upmcgrid() (code # ',bnd_shape,').'
      call fexit()
    end if
! all-wall BCs
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   rectangular box
    if (bnd_shape.eq.1) then
      do i=1,3
        grid%origin(i) = bnd_params(3+i) - 0.25*bnd_params(i)
        grid%deltas(i) = (1.5*bnd_params(i))/(1.0*grid%dim(i))
      end do
!   spherical system
    else if (bnd_shape.eq.2) then
      do i=1,3
        grid%origin(i) = bnd_params(i)-1.5*bnd_params(4)
        grid%deltas(i) = 3.0*bnd_params(4)/(1.0*grid%dim(i))
      end do
!   cylinder
    else if (bnd_shape.eq.3) then
      grid%origin(1:2) = bnd_params(1:2) - 1.5*bnd_params(4)
      grid%deltas(1:2) = 3.0*bnd_params(4)/(1.0*grid%dim(1:2))
      grid%origin(3) = bnd_params(3) - 0.75*bnd_params(6)
      grid%deltas(3) = 1.5*bnd_params(6)/(1.0*grid%dim(3))
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in set&
 &upmcgrid() (code # ',bnd_shape,').'
      call fexit()
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & in setupmcgrid() (code # ',bnd_type,').'
    call fexit()
  end if
  grid%diags(1,1,1) = 0.0
  grid%diags(2,2,2) = sqrt(grid%deltas(1)*grid%deltas(1)&
 &                     + grid%deltas(2)*grid%deltas(2)&
 &                     + grid%deltas(3)*grid%deltas(3))
  grid%diags(2,2,1) = sqrt(grid%deltas(1)*grid%deltas(1)&
 &                     + grid%deltas(2)*grid%deltas(2))
  grid%diags(2,1,2) = sqrt(grid%deltas(1)*grid%deltas(1)&
 &                     + grid%deltas(3)*grid%deltas(3))
  grid%diags(1,2,2) = sqrt(grid%deltas(2)*grid%deltas(2)&
 &                     + grid%deltas(3)*grid%deltas(3))
  grid%diags(2,1,1) = grid%deltas(1)
  grid%diags(1,2,1) = grid%deltas(2)
  grid%diags(1,1,2) = grid%deltas(3)
  grid%pnts = grid%dim(1)*grid%dim(2)*grid%dim(3)
!
! now major allocation for NB lists (these arrays easily get very large)
  allocate(grid%nblist(grid%maxnbs,grid%pnts+1))
  allocate(grid%gpnblist(grid%pnts+1,grid%maxgpnbs))
  allocate(grid%gpnbmd(grid%pnts+1,grid%maxgpnbs))
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    allocate(grid%b_gpnblist(grid%pnts+1,grid%maxgpnbs))
    allocate(grid%b_gpnbmd(grid%pnts+1,grid%maxgpnbs))
  end if
!
  do i=1,grid%pnts
    do j=i+1,grid%pnts
      call gridmindiff(i,j,pvd2)
      if (pvd2.lt.cutplusmargin2) then
        grid%gpnbnr(i) = grid%gpnbnr(i) + 1
        if (grid%gpnbnr(i).gt.grid%maxgpnbs) then
          write(ilog,*) 'Fatal. Exceeded static neighbor number for grid points in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXGPNB or choose a coarser grid.'
          call fexit()
        end if
        grid%gpnblist(i,grid%gpnbnr(i)) = j
        grid%gpnbmd(i,grid%gpnbnr(i)) = floor(sqrt(pvd2))
!       for MC, we need the complete double-mapped list
        if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
          grid%b_gpnbnr(j) = grid%b_gpnbnr(j) + 1
          if (grid%b_gpnbnr(j).gt.grid%maxgpnbs) then
            write(ilog,*) 'Fatal. Exceeded static neighbor number for grid points in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXGPNB or choose a coarser grid.'
            call fexit()
          end if
          grid%b_gpnblist(j,grid%b_gpnbnr(j)) = i
          grid%b_gpnbmd(j,grid%b_gpnbnr(j)) = floor(sqrt(pvd2))
        end if    
      end if
    end do
  end do
!
! assign the various residues to 'their' grid points (via C-alpha or other reference atom)
!
  grid%maxr(:) = 0.0
!
  if (n_topgrps.gt.0) then
    grid%maxr(:) = 1.2
    do i=1,n_topgrps
      refp(1) = topgrpxyz(1,i)
      refp(2) = topgrpxyz(2,i)
      refp(3) = topgrpxyz(3,i)
      which = gridnavi(refp)
      grid%nbnr(which) = grid%nbnr(which) + 1
      if (grid%nbnr(which).gt.grid%maxnbs) then
        write(ilog,*) 'Fatal. Exceeded static maximum for groups associated to a single grid point in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXRSNB or choose a finer grid.'
        call fexit()
      end if
      grid%nblist(grid%nbnr(which),which) = i
      grid%tgrpgp(i) = which
    end do
  else
    do i=1,nseq
      refp(1) = x(refat(i))
      refp(2) = y(refat(i))
      refp(3) = z(refat(i))
      which = gridnavi(refp)
      grid%nbnr(which) = grid%nbnr(which) + 1
      grid%maxr(which) = max(grid%maxr(which),resrad(i))
      if (grid%nbnr(which).gt.grid%maxnbs) then
        write(ilog,*) 'Fatal. Exceeded static maximum for residues associated to a single grid point in setupmcgrid().&
 & Increase the value for FMCSC_GRIDMAXRSNB or choose a finer grid.'
        call fexit()
      end if
      grid%nblist(grid%nbnr(which),which) = i
      grid%resgp(i) = which
    end do
  end if
!
end
!
!-----------------------------------------------------------------------------
!
! an experimental routine to map identical particles spatially to optimize memory access
!
subroutine map_molecules()
!
  use atoms
  use cutoffs
  use interfaces
  use polypep
  use sequen
  use system
!
  implicit none
!
  integer i,j,k,ii,i0,ik,nwa,yi,zi,crosses,cccs(2),cc
  logical atrue
  RTYPE ccvs(2)
  integer bla(nseq,2)
!
  atrue = .true.
!
  if (use_waterloops.EQV..true.) then
!
    i0 = 1
    ii = nseq - rsw1 + 1
    nwa = at(rsw1)%na - 1
    ik = ii
    if (ii.le.10) return
!
    atinfo(1,1:n) = x(1:n)
    atinfo(2,1:n) = y(1:n)
    atinfo(3,1:n) = z(1:n)
!
    cccs(1:2) = nint(bnd_params(2:3)/(sqrt(2.0)*mcel_cutoff))
    ccvs(1:2) = bnd_params(2:3)/(1.0*cccs(1:2))
    crosses = cccs(1)*cccs(2)
    vhlper(1:ii,1) = x(refat(rsw1:nseq)) 
    do i=1,ii
      bla(i,1) = i
    end do
    call merge_sort(ldim=ik,up=atrue,list=vhlper(1:ii,1),olist=vhlper(1:ii,2),ilo=i0,ihi=ii,&
 &                  idxmap=bla(1:ii,1),olist2=bla(1:ii,2))
    ii = rsw1-1
    do cc=1,crosses
      do i=rsw1,nseq
        k = bla(i-rsw1+1,2) + rsw1 - 1
        yi = floor(min(bnd_params(2),max(0.0,(atinfo(2,refat(k)) - bnd_params(5))))/ccvs(1)) + 1
        zi = floor(min(bnd_params(3),max(0.0,(atinfo(3,refat(k)) - bnd_params(6))))/ccvs(2)) + 1
        if (yi.gt.cccs(1)) yi = 1
        if (yi.le.0) yi = cccs(1)
        if (zi.gt.cccs(2)) zi = 1
        if (zi.le.0) zi = cccs(2)
        if (((yi-1)*cccs(2)+zi).eq.cc) then
          ii = ii + 1
          do j=0,nwa
            x(at(ii)%bb(1)+j) = atinfo(1,at(k)%bb(1)+j)
            y(at(ii)%bb(1)+j) = atinfo(2,at(k)%bb(1)+j)
            z(at(ii)%bb(1)+j) = atinfo(3,at(k)%bb(1)+j)
          end do
        end if
      end do
    end do
    write(*,*) ii
    do i=rsw1,nseq
      write(*,*) x(at(i)%bb(1))
    end do    
!
  end if
!
end
!
! 
!
! this function must be called when the forces (incl. TR/LR) are all 0.0 or a going to be set to zero before 
! their next use
! this includes the solvation model derivatives/temp. arrays sav_dr, sum_scrcbs, sum_scrcbs_tr, svte, sisa%nix
!
subroutine tryit()
!
  use cutoffs
  use atoms, ONLY: x,y,z,n
  use sequen, ONLY: nseq
  use polypep, ONLY: at
  use energies, ONLY: is_tab
  use movesets, ONLY: skip_frz
  use mcgrid
!
  implicit none
!
  integer i,j,k,nwa,ii,jj,ki(3),ks(3),mapm(3),m1,m2,m3,k1,k2,k3
  integer map(nseq,2)
!
  if (use_waterloops.EQV..true.) then
    ii = at(rsw1)%bb(1)
    atinfo(1,ii:n) = x(ii:n)
    atinfo(2,ii:n) = y(ii:n)
    atinfo(3,ii:n) = z(ii:n)
    x(ii:n) = 0.
 
    map(:,:) = 0
    nwa = at(rsw1)%na - 1
!
    if (use_mcgrid.EQV..true.) then
      k = rsw1-1
      ks(:) = nint(0.5*grid%dim(1:3))
      do i=1,3
        mapm(i) = i
      end do
      m1 = 1
      m2 = 2
      m3 = 3
      j = grid%dim(2)*grid%dim(3)*(ks(1)-1) + grid%dim(3)*(ks(2)-1) + ks(3)
      do i=1,grid%nbnr(j)
        if (grid%nblist(i,j).ge.rsw1) then
          k = k + 1
          map(k,1) = grid%nblist(i,j)
          map(grid%nblist(i,j),2) = k
        end if
      end do
      write(*,*) ks,j
      do jj=1,maxval(ks)+1
        do k1=0,2*jj
          if (mod(k1,2).eq.0) then
            ki(1) = k1/2
          else
            ki(1) = -(k1+1)/2
          end if
          if ((ks(m1)+ki(1).le.0).OR.(ks(m1)+ki(1)).gt.grid%dim(m1)) cycle
          do k2=0,2*jj
            if (mod(k2,2).eq.0) then
              ki(2) = k2/2
            else
              ki(2) = -(k2+1)/2
            end if
            if ((ks(m2)+ki(2).le.0).OR.(ks(m2)+ki(2)).gt.grid%dim(m2)) cycle
            do k3=0,2*jj
              if (mod(k3,2).eq.0) then
                ki(3) = k3/2
              else
                ki(3) = -(k3+1)/2
              end if
              if (maxval(abs(ki)).lt.jj) cycle
              if ((ks(m3)+ki(3).le.0).OR.(ks(m3)+ki(3)).gt.grid%dim(m3)) cycle
              j = grid%dim(2)*grid%dim(3)*(ks(1)+ki(mapm(1))-1) + grid%dim(3)*(ks(2)+ki(mapm(2))-1) + ks(3)+ki(mapm(3))
!              write(*,*) j,grid%pnts,jj
              do i=1,grid%nbnr(j)
                if (grid%nblist(i,j).ge.rsw1) then
                  k = k + 1
                  map(k,1) = grid%nblist(i,j)
                  map(grid%nblist(i,j),2) = k
                end if
              end do
            end do
          end do
        end do
      end do         
!      do j=1,grid%pnts
!        do i=1,grid%nbnr(j)
!          if (grid%nblist(i,j).ge.rsw1) then
!            k = k + 1
!            map(k,1) = grid%nblist(i,j)
!            map(grid%nblist(i,j),2) = k
!          end if
!        end do
!      end do
      do i=rsw1,nseq
        do j=0,nwa
          x(at(i)%bb(1)+j) = atinfo(1,at(map(i,1))%bb(1)+j)
          y(at(i)%bb(1)+j) = atinfo(2,at(map(i,1))%bb(1)+j)
          z(at(i)%bb(1)+j) = atinfo(3,at(map(i,1))%bb(1)+j)
        end do
      end do
    else

      do i=rsw1,nseq
        if (rs_nbl(i)%nwnbs.gt.0) map(i,1) = maxval(rs_nbl(i)%wnb(1:rs_nbl(i)%nwnbs)) - minval(rs_nbl(i)%wnb(1:rs_nbl(i)%nwnbs))
      end do
      write(*,*) sum(map(rsw1:nseq,1)),maxval(map(rsw1:nseq,1))
      write(*,55) rs_nbl(1)%wnb(1:rs_nbl(1)%nwnbs)

      map(:,:) = 0
      k = rsw1-1
      do i=rsw1,nseq
        if (i.eq.rsw1) then
          k = k + 1
          map(i,2) = k
          map(k,1) = i
        end if
        do j=1,rs_nbl(i)%nwnbs
          if (map(rs_nbl(i)%wnb(j),2).le.0) then
            k = k + 1
            map(k,1) = rs_nbl(i)%wnb(j)
            map(rs_nbl(i)%wnb(j),2) = k
          end if
        end do
      end do
      write(*,*) k,rsw1,nseq
      do i=rsw1,nseq
        if (map(i,2).le.0) then
          k = k + 1
          write(*,*) i,k
          map(i,2) = k
          map(k,1) = i
        end if
        do j=0,nwa
          x(at(i)%bb(1)+j) = atinfo(1,at(map(i,1))%bb(1)+j)
          y(at(i)%bb(1)+j) = atinfo(2,at(map(i,1))%bb(1)+j)
          z(at(i)%bb(1)+j) = atinfo(3,at(map(i,1))%bb(1)+j)
        end do
      end do
    end if

!
    atinfo(1,ii:n) = x(ii:n)
    atinfo(2,ii:n) = y(ii:n)
    atinfo(3,ii:n) = z(ii:n)
!
    if ((use_mcgrid.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
      rs_nbl(:)%ntmpanb = 0
      do i=rsw1,nseq
        call updateresgp(i)
      end do
    end if
    call all_respairs_nbl(skip_frz,is_tab,0)
    map(:,:) = 0
    do i=rsw1,nseq
      if (rs_nbl(i)%nwnbs.gt.0) map(i,1) = maxval(rs_nbl(i)%wnb(1:rs_nbl(i)%nwnbs)) - minval(rs_nbl(i)%wnb(1:rs_nbl(i)%nwnbs))
    end do
    write(*,*) sum(map(rsw1:nseq,1)),maxval(map(rsw1:nseq,1))
    write(*,55) rs_nbl(1)%wnb(1:rs_nbl(1)%nwnbs)
 55 format(10000(i4,1x))
  end if
!
end
!
!
!------------------------------------------------------------------------
!
function gridnavi(refp)
!
  use iounit
  use mcgrid
  use system
  use movesets
!
  implicit none
!
  integer xi,yi,zi,gridnavi,k(3),i
  RTYPE refp(3)
! see above
!
! get the grid indexing
!
  xi = floor(((refp(1)-grid%origin(1))/grid%deltas(1))) + 1
  yi = floor(((refp(2)-grid%origin(2))/grid%deltas(2))) + 1
  zi = floor(((refp(3)-grid%origin(3))/grid%deltas(3))) + 1
!
! periodic BC (PBC) (note that global shifts are molecule-wise and that refat(rs) may consequently
! hang "off" the grid at all times) -> shift back
  if (bnd_type.eq.1) then
!   cubic box
    if (bnd_shape.eq.1) then
      do i=1,3
        k(i) = 0
        if (refp(i).lt.bnd_params(3+i)) then
          k(i) = floor((bnd_params(3+i)-refp(i))/bnd_params(i)) + 1
        else if(refp(i).gt.(bnd_params(3+i)+bnd_params(i))) then
          k(i) = floor((refp(i)-(bnd_params(3+i)+bnd_params(i)))/&
 &                                             bnd_params(i)) + 1
          k(i) = -k(i)
        end if
      end do
    else if (bnd_shape.eq.3) then
      k(:) = 0
      if (refp(3).lt.(bnd_params(3)-0.5*bnd_params(6))) then
        k(3) = floor((bnd_params(3)-0.5*bnd_params(6)-refp(3))/bnd_params(6)) + 1
      else if(refp(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
        k(3) = floor((refp(3)-(bnd_params(3)+0.5*bnd_params(6)))/bnd_params(6)) + 1
        k(3) = -k(3)
      end if
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in gri&
 &dnavi() (code # ',bnd_shape,').'
      call fexit()
    end if
    xi = xi + k(1)*grid%dim(1)
    yi = yi + k(2)*grid%dim(2)
    zi = zi + k(3)*grid%dim(3)
  end if
!
 45 format(' x: ',g12.4,' y: ',g12.4,' z: ',g12.4)
  if ((xi.le.0).OR.(xi.gt.grid%dim(1))) then
    xi = max(xi,1)
    xi = min(xi,grid%dim(1))
!   in MC this happens -> lump into outermost grid-point
    if ((use_dyn.EQV..true.).OR.((bnd_type.eq.1).AND.(bnd_shape.eq.1))) then
      if ((bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
      else if (.NOT.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..false.))) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
        write(ilog,*) 'In non-periodic, soft boundaries, this will usually indicate that&
 & the system is poorly equilibrated and/or that the simulation has become unstable.'
      end if
    end if
  end if
  if ((yi.le.0).OR.(yi.gt.grid%dim(2))) then
    yi = max(yi,1)
    yi = min(yi,grid%dim(2))
!   in MC this happens -> lump into outermost grid-point
    if ((use_dyn.EQV..true.).OR.((bnd_type.eq.1).AND.(bnd_shape.eq.1))) then
      if ((bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
      else if (.NOT.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..false.))) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
        write(ilog,*) 'In non-periodic, soft boundaries, this will usually indicate that&
 & the system is poorly equilibrated and/or that the simulation has become unstable.'
      end if
    end if
  end if
  if ((zi.le.0).OR.(zi.gt.grid%dim(3))) then
    zi = max(zi,1)
    zi = min(zi,grid%dim(3))
!   in MC this happens -> lump into outermost grid-point
    if ((use_dyn.EQV..true.).OR.(bnd_type.eq.1)) then
      if (bnd_type.eq.1) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
      else if (.NOT.(((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(in_dyncyc.EQV..false.))) then
        write(ilog,*) 'Warning. Position is off grid:'
        write(ilog,45) refp(1),refp(2),refp(3)
        write(ilog,*) 'In non-periodic, soft boundaries, this will usually indicate that&
 & the system is poorly equilibrated and/or that the simulation has become unstable.'
      end if
    end if
  end if
!
  gridnavi = grid%dim(2)*grid%dim(3)*(xi-1) + grid%dim(3)*(yi-1)+ zi
  return
!
end
!
!-----------------------------------------------------------------------
!
subroutine gridpvec(w,pvec)
!
  use mcgrid
!
  implicit none
!
  integer w,xi,yi,zi
  RTYPE pvec(3)
!
  call getgridtriple(w,xi,yi,zi)
  pvec(1) = grid%origin(1) + (xi-0.5)*grid%deltas(1)
  pvec(2) = grid%origin(2) + (yi-0.5)*grid%deltas(2)
  pvec(3) = grid%origin(3) + (zi-0.5)*grid%deltas(3)
!
end
!
!-----------------------------------------------------------------------
!
subroutine gridmindiff(w1,w2,dis2)
!
  use iounit
  use mcgrid
  use system
!
  implicit none
!
  integer w1,w2,xi1,xi2,yi1,yi2,zi1,zi2,dxi,dyi,dzi
  RTYPE dvec(3),dis2
!
  call getgridtriple(w1,xi1,yi1,zi1)
  call getgridtriple(w2,xi2,yi2,zi2)
!
  dxi = abs(xi2 - xi1)
  dyi = abs(yi2 - yi1)
  dzi = abs(zi2 - zi1)
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
      if (dxi.gt.(grid%dim(1)/2)) dxi = -dxi + grid%dim(1)
      if (dyi.gt.(grid%dim(2)/2)) dyi = -dyi + grid%dim(2)
      if (dzi.gt.(grid%dim(3)/2)) dzi = -dzi + grid%dim(3)
    else if (bnd_shape.eq.3) then
      if (dzi.gt.(grid%dim(3)/2)) dzi = -dzi + grid%dim(3)
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in grid&
 &diff(...) (code # ',bnd_shape,').'
      call fexit()
    end if
  else if ((bnd_type.ge.2).AND.(bnd_type.le.4)) then
!   do nothing
  else
    write(ilog,*) 'Fatal. Encountered unsupported boundary condition&
 & in griddiff(...) (code # ',bnd_type,').'
    call fexit()
  end if
  dvec(1) = max(0,dxi-1)*grid%deltas(1)
  dvec(2) = max(0,dyi-1)*grid%deltas(2)
  dvec(3) = max(0,dzi-1)*grid%deltas(3)
  dis2 = sum(dvec(:)*dvec(:))
!
end
!
!------------------------------------------------------------------------
!
subroutine getgridtriple(w,xi,yi,zi)
!
  use iounit
  use mcgrid
!
  implicit none
!
  integer buf,w,xi,yi,zi
!
  buf = w
  xi = int((w-1)/(grid%dim(2)*grid%dim(3))) + 1
  buf = buf - (xi-1)*grid%dim(2)*grid%dim(3)
!
  yi = int((buf-1)/grid%dim(3)) + 1
  zi = buf - (yi-1)*grid%dim(3)
!
  buf = grid%dim(2)*grid%dim(3)*(xi-1) + grid%dim(3)*(yi-1) + zi
  if (buf.ne.w) write(*,*) 'NOOOOOOOOOOOO'
!
end
!
!------------------------------------------------------------------------
!
subroutine updateresgp(rs)
!
  use iounit
  use sequen
  use atoms
  use mcgrid
!
  implicit none
!
  integer, INTENT(IN):: rs
!
  RTYPE refp(3)
  integer gridnavi,wo,wn,i,j
!
  wo = grid%resgp(rs)
  !write(ilog,*) 'GRP# ',grpnr, 'OLD# ',wo, 'W/ ',grid%nbnr(wo)
  refp(1) = x(refat(rs))
  refp(2) = y(refat(rs))
  refp(3) = z(refat(rs))
  wn = gridnavi(refp)
  if ((wo.le.0).OR.(wo.gt.grid%pnts)) then
    write(ilog,*) 'Fatal. Grid-point was misassigned. This is either an exploded simulation or a bug.'
    call fexit()
  end if
!
  if (wn.eq.wo) then !nothing to do
    return 
  else !gotta switch group from old gp to new one
!    write(*,*) 'updating ',wo,' to ',wn,' for ',rs
#ifdef ENABLE_THREADS
!$OMP CRITICAL(GPRS_UPDATE)
#endif
    do i=1,grid%nbnr(wo)
      if (rs.eq.grid%nblist(i,wo)) then
        do j=i+1,grid%nbnr(wo)
          grid%nblist(j-1,wo) = grid%nblist(j,wo)
        end do
        exit
      end if
    end do
    grid%nbnr(wo) = grid%nbnr(wo) - 1
    grid%nbnr(wn) = grid%nbnr(wn) + 1
    grid%maxr(wn) = max(grid%maxr(wn),resrad(rs))
    if (grid%nbnr(wn).gt.grid%maxnbs) then
      write(ilog,*) 'Warning. Exceeded static maximum for groups associated to a single grid point in updateresgp().&
 & Increasing the value for FMCSC_GRIDMAXRSNB dynamically which may cause memory exceptions.'
      call gridnblrsz()
    end if
    grid%nblist(grid%nbnr(wn),wn) = rs
    grid%resgp(rs) = wn
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(GPRS_UPDATE)
#endif
  end if
!
end
!
!------------------------------------------------------------------------
!
subroutine updatetgrpgp(tgi)
!
  use iounit
  use sequen
  use cutoffs
  use mcgrid
!
  implicit none
!
  integer, INTENT(IN):: tgi
!
  RTYPE refp(3)
  integer gridnavi,wo,wn,i,j
!
  wo = grid%tgrpgp(tgi)
  !write(ilog,*) 'GRP# ',grpnr, 'OLD# ',wo, 'W/ ',grid%nbnr(wo)
  refp(1) = topgrpxyz(1,tgi)
  refp(2) = topgrpxyz(2,tgi)
  refp(3) = topgrpxyz(3,tgi)
  wn = gridnavi(refp)
  if ((wo.le.0).OR.(wo.gt.grid%pnts)) then
    write(ilog,*) 'Fatal. Grid-point was misassigned. This is either an exploded simulation or a bug.'
    call fexit()
  end if
!
  if (wn.eq.wo) then !nothing to do
    return 
  else !gotta switch group from old gp to new one
!    write(*,*) 'updating ',wo,' to ',wn,' for ',tgi
#ifdef ENABLE_THREADS
!$OMP CRITICAL(GPRS_UPDATE)
#endif
    do i=1,grid%nbnr(wo)
      if (tgi.eq.grid%nblist(i,wo)) then
        do j=i+1,grid%nbnr(wo)
          grid%nblist(j-1,wo) = grid%nblist(j,wo)
        end do
        exit
      end if
    end do
    grid%nbnr(wo) = grid%nbnr(wo) - 1
    grid%nbnr(wn) = grid%nbnr(wn) + 1
    if (grid%nbnr(wn).gt.grid%maxnbs) then
      write(ilog,*) 'Warning. Exceeded static maximum for groups associated to a single grid point in updatetgrpgp().&
 & Increasing the value for FMCSC_GRIDMAXRSNB dynamically which may cause memory exceptions.'
      call gridnblrsz()
    end if
    grid%nblist(grid%nbnr(wn),wn) = tgi
    grid%tgrpgp(tgi) = wn
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(GPRS_UPDATE)
#endif
  end if
!
end
!
!-----------------------------------------------------------------
!
subroutine gridnblrsz()
!
  use mcgrid
  use iounit
!
  implicit none
!
  integer, ALLOCATABLE:: onelist(:,:)
  integer i
!
  allocate(onelist(grid%maxnbs,grid%pnts))
!
  do i=1,grid%pnts
    onelist(1:grid%maxnbs,i) = grid%nblist(1:grid%maxnbs,i)
  end do
  deallocate(grid%nblist)
  grid%maxnbs = ceiling(1.2*grid%maxnbs)
  allocate(grid%nblist(grid%maxnbs,grid%pnts+1))
  do i=1,grid%pnts
    grid%nblist(1:grid%maxnbs,i) = onelist(1:grid%maxnbs,i)
  end do
!
  deallocate(onelist)
!
  end
!
!-----------------------------------------------------------------
!
subroutine gridreport()
!
  use atoms
  use sequen
  use mcgrid
  use iounit
  use aminos
  use cutoffs
!
  implicit none
!
  integer i,j,k,kk,tg
  character(3) resname
!
 20   format('# ',i6,' (',a3,'): Ref. at ',f7.3,1x,f7.3,1x,f7.3)
 21   format('Gridpoint ',i6,' has ',i4,' residues.')
 22   format('There are ',i8,' residues on the grid.')
 23   format('Gridpoint ',i6,' has ',i4,' groups.')
 24   format('There are ',i8,' groups on the grid.')

  write(ilog,*)
  write(ilog,*) '--- Summary of grid occupation ---'
  write(ilog,*)
!
  if (n_topgrps.gt.0) then
    tg = 0
    do i=1,grid%dim(1)*grid%dim(2)*grid%dim(3)
      if (grid%nbnr(i).gt.0) then
        write(ilog,23) i,grid%nbnr(i)
      end if
      tg = tg + grid%nbnr(i)
    end do
    write(ilog,24) tg
  else
    tg = 0
    do i=1,grid%dim(1)*grid%dim(2)*grid%dim(3)
      if (grid%nbnr(i).gt.0) then
        write(ilog,21) i,grid%nbnr(i)
      end if
      tg = tg + grid%nbnr(i)
      do j=1,grid%nbnr(i)
        k = grid%nblist(j,i)
        resname = amino(seqtyp(k))
        kk = refat(k)
        write(ilog,20) k,resname,x(kk),y(kk),z(kk)
      end do
    end do
    write(ilog,22) tg
  end if
  write(ilog,*)
!
end
!
!-----------------------------------------------------------------------------
!
! in both grd_respairs() and respairs() we're going to use one diagonal
! half-matrix (first index larger) for the standard IPP/LJ cutoff
! and the other one (first index smaller) for the EL cutoff
! additional terms are turned on on the EL side for charged species
!
subroutine grd_respairs(rs)
!
  use iounit
  use cutoffs
  use sequen
  use mcgrid
  use energies
  use system
  use atoms
!
  implicit none
!
  integer rs,rs2,nbi,ii,i,w
  RTYPE dis,dis2
!
  w = grid%resgp(rs)
!
! add groups associated with same grid-point (incl. self-term)
  do nbi=1,grid%nbnr(w)
    rs2 = grid%nblist(nbi,w)
!   with extra check
    call dis_bound(refat(rs),refat(rs2),dis2)
    dis = sqrt(dis2)
!   now check for the short-range cutoff
    if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(rs2))) then
      if (rs2.ge.rs) then
        rsp_mat(rs2,rs) = 1
      else
        rsp_mat(rs,rs2) = 1
      end if
    end if
    if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!   within standard "twin"-range second cutoff (TABUL, POLAR)
    if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(rs2))) then
      if (rs2.ge.rs) then
        rsp_mat(rs,rs2) = 1
      else
        rsp_mat(rs2,rs) = 1
      end if 
    end if
!   this would be the version without extra check
!    rsp_mat(rs,rs2) = 1
!    rsp_mat(rs2,rs) = 1
  end do
! scan gp neighbors
  do i=1,grid%gpnbnr(w)
    ii = grid%gpnblist(w,i)
    do nbi=1,grid%nbnr(ii)
      rs2 = grid%nblist(nbi,ii)
!     with extra explicit check
      call dis_bound(refat(rs),refat(rs2),dis2)
      dis = sqrt(dis2)
!     now check for the short-range cutoff
      if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs2,rs) = 1
        else
          rsp_mat(rs,rs2) = 1
        end if
      end if
      if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!     within standard "twin"-range second cutoff (TABUL, POLAR)
      if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs,rs2) = 1
        else
          rsp_mat(rs2,rs) = 1
        end if 
      end if
!     this would be the version without extra explicit check
!     this additional screen is expensive but crucial if the resrads in the simulation are very different
!     such as in the case of a macromolecule in a solvent box
!      chkdis = 1.0*grid%gpnbmd(w,i)-resrad(rs)-resrad(rs2)
!      if (chkdis.gt.mcel_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs2,rs) = 1
!      else
!        rsp_mat(rs,rs2) = 1
!      end if
!      if (chkdis.gt.mcnb_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs,rs2) = 1
!      else
!        rsp_mat(rs2,rs) = 1
!      end if
    end do
  end do
  do i=1,grid%b_gpnbnr(w)
    ii = grid%b_gpnblist(w,i)
    do nbi=1,grid%nbnr(ii)
      rs2 = grid%nblist(nbi,ii)
!     with extra explicit check
      call dis_bound(refat(rs),refat(rs2),dis2)
      dis = sqrt(dis2)
!     now check for the short-range cutoff
      if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs2,rs) = 1
        else
          rsp_mat(rs,rs2) = 1
        end if
      end if
      if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!     within standard "twin"-range second cutoff (TABUL, POLAR)
      if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(rs2))) then
        if (rs2.ge.rs) then
          rsp_mat(rs,rs2) = 1
        else
          rsp_mat(rs2,rs) = 1
        end if 
      end if
!     this would be the version without extra explicit check
!     this additional screen is expensive but crucial if the resrads in the simulation are very different
!     such as in the case of a macromolecule in a solvent box
!      chkdis = 1.0*grid%b_gpnbmd(w,i)-resrad(rs)-resrad(rs2)
!      if (chkdis.gt.mcel_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs2,rs) = 1
!      else
!        rsp_mat(rs,rs2) = 1
!      end if
!      if (chkdis.gt.mcnb_cutoff) cycle
!      if (rs.ge.rs2) then
!        rsp_mat(rs,rs2) = 1
!      else
!        rsp_mat(rs2,rs) = 1
!      end if
    end do
  end do
!
  if (use_POLAR.EQV..false.) return
!
  if (lrel_mc.eq.1) then
    if (chgflag(rs).EQV..true.) then
      do rs2=1,rs
        if (rsp_mat(rs2,rs).eq.0) rsp_mat(rs2,rs) = 2
      end do
      do rs2=rs+1,nseq
        if (rsp_mat(rs,rs2).eq.0) rsp_mat(rs,rs2) = 2
      end do
    else
      do i=1,cglst%ncs
        ii = cglst%it(i)
        rs2 = atmres(cglst%it(i))
        if (chgflag(rs2).EQV..false.) cycle ! may be patched
        if (rs2.gt.rs) then
          if (rsp_mat(rs,rs2).eq.0) rsp_mat(rs,rs2) = 2
        else
          if (rsp_mat(rs2,rs).eq.0) rsp_mat(rs2,rs) = 2
        end if
      end do
    end if
  else if ((lrel_mc.eq.2).OR.(lrel_mc.eq.3)) then
    if (chgflag(rs).EQV..true.) then
      do i=1,cglst%ncs
        ii = cglst%it(i)
        rs2 = atmres(cglst%it(i))
        if (chgflag(rs2).EQV..false.) cycle ! may be patched
        if (rs2.gt.rs) then
          if (rsp_mat(rs,rs2).eq.0) rsp_mat(rs,rs2) = 2
        else
          if (rsp_mat(rs2,rs).eq.0) rsp_mat(rs2,rs) = 2
        end if
      end do
    end if
! do nothing for lrel_mc = 4
  end if
!
end
!
!-----------------------------------------------------------------------------
!
  function gridwrap(v1,d1)
!
  implicit none
!
  integer v1,d1,gridwrap
!
  if (v1.gt.d1) gridwrap = v1 - d1
  if (v1.lt.1) gridwrap = v1 + d1
!
  end 
!
!-----------------------------------------------------------------------------
!
! this is the dynamics-relevant neighbor list routine
! which creates (somewhat smart) neighbor list candidate index lists, neighbor lists themselves,
! and does all the memory allocation
!
subroutine all_respairs_nbl(cycle_frz,cycle_tab,tpi)
!
  use iounit
  use cutoffs
  use sequen
  use mcgrid
  use energies
  use system
  use tabpot
  use mcsums
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  logical, INTENT(IN):: cycle_frz,cycle_tab
  integer, INTENT(IN):: tpi
!
  integer i,j,rs1,rs2,nbi,nbj,jj,rs,sta,sto,aone,atwo,fixiv,sta2,sto2
  RTYPE frs1(grid%maxnbs)
  integer lrs1(grid%maxnbs),lrs2(grid%maxnbs)
#ifdef ENABLE_THREADS
  integer(KIND=8) ttimer
  logical OMP_IN_PARALLEL
#endif
!
  fixiv = 10
!
#ifdef ENABLE_THREADS
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, all_respairs_nbl(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
  if (tpi.gt.0) then
    if (thr_dlb(1,1).gt.0) then
      if (tpi.eq.1) thr_dlb(1,2) = thr_dlb(1,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(1,tpi) = thr_timings(1,tpi) + ttimer
    end if
    sta = thr_limits(5,tpi)
    sto = thr_limits(6,tpi)
    sta2 = thr_limits(7,tpi)
    sto2 = thr_limits(8,tpi)
  else
    sta = 1
    sto = grid%pnts
    sta2 = 1
    sto2 = nseq
  end if
#else
  sta = 1
  sto = grid%pnts
  sta2 = 1
  sto2 = nseq
#endif
!
  aone = 1
  atwo = 2
!
  if (use_cutoffs.EQV..true.) then
!
    if (use_mcgrid.EQV..true.) then
!
      do i=sta,sto
!       when using threads, this update creates a tolerable race condition for read access later
        grid%maxr(i) = maxval(resrad(grid%nblist(1:grid%nbnr(i),i)))
      end do
      if (cycle_tab.EQV..true.) then
        do i=sta,sto
!         add all pairs within two pre-lists
          do nbi=1,grid%nbnr(i)-1
            rs1 = grid%nblist(nbi,i)
            rs2 = 0
            do nbj=nbi+1,grid%nbnr(i)
              if (tbp%rsmat(rs1,grid%nblist(nbj,i)).gt.0) then
                rs2 = rs2 + 1
                lrs1(rs2) = grid%nblist(nbj,i)
              end if
            end do
            rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + rs2
            if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
              call nbl_resz(rs1,7)
            end if
            rs_nbl(rs1)%tmpanb((rs_nbl(rs1)%ntmpanb-rs2+1):rs_nbl(rs1)%ntmpanb) = lrs1(1:rs2)
          end do
!         for neighboring grid-points always add everything to residue associated with point i
          do j=1,grid%gpnbnr(i)
            jj = grid%gpnblist(i,j)
!
            if (mcel_cutoff.lt.(1.0*grid%gpnbmd(i,j) - grid%maxr(i) - grid%maxr(jj))) cycle
            lrs1(1:grid%nbnr(jj)) = grid%nblist(1:grid%nbnr(jj),jj)
!
            if (mcel_cutoff.ge.1.0*grid%gpnbmd(i,j)) then ! add all without checks
              do nbi=1,grid%nbnr(i)
                rs1 = grid%nblist(nbi,i)
                rs2 = 0
                do nbj=1,grid%nbnr(jj)
                  if (tbp%rsmat(rs1,lrs1(nbj)).gt.0) then
                    rs2 = rs2 + 1
                    lrs2(rs2) = lrs1(nbj)
                  end if
                end do
                rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + rs2
                if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
                  call nbl_resz(rs1,7)
                end if
                rs_nbl(rs1)%tmpanb((rs_nbl(rs1)%ntmpanb-rs2+1):rs_nbl(rs1)%ntmpanb) = lrs2(1:rs2)
              end do
            else
              frs1(1:grid%nbnr(jj)) = 1.0*grid%gpnbmd(i,j) - resrad(lrs1(1:grid%nbnr(jj)))
              do nbi=1,grid%nbnr(i)
                rs1 = grid%nblist(nbi,i)
                frs1(1:grid%nbnr(jj)) = frs1(1:grid%nbnr(jj)) - resrad(rs1)
                rs2 = 0
                do nbj=1,grid%nbnr(jj)
                  if ((frs1(nbj).le.mcel_cutoff).AND.(tbp%rsmat(rs1,lrs1(nbj)).gt.0)) then
                    rs2 = rs2 + 1
                    lrs2(rs2) = lrs1(nbj)
                  end if
                end do
                rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + rs2
                if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
                  call nbl_resz(rs1,7)
                end if
                rs_nbl(rs1)%tmpanb((rs_nbl(rs1)%ntmpanb-rs2+1):rs_nbl(rs1)%ntmpanb) = lrs2(1:rs2)
              end do
            end if
          end do
        end do
      else
        do i=sta,sto
          do nbi=1,grid%nbnr(i)-1
            rs1 = grid%nblist(nbi,i)
            lrs1(1:(grid%nbnr(i)-nbi)) = grid%nblist((nbi+1):grid%nbnr(i),i)
            rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + grid%nbnr(i) - nbi
            if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
              call nbl_resz(rs1,7)
            end if
            rs_nbl(rs1)%tmpanb((rs_nbl(rs1)%ntmpanb+nbi-grid%nbnr(i)+1):rs_nbl(rs1)%ntmpanb) = lrs1(1:(grid%nbnr(i)-nbi))
          end do
!         for neighboring grid-points always add everything to residue associated with point i
          do j=1,grid%gpnbnr(i)
            jj = grid%gpnblist(i,j)
!
            if (mcel_cutoff.lt.(1.0*grid%gpnbmd(i,j) - grid%maxr(i) - grid%maxr(jj))) cycle
            lrs1(1:grid%nbnr(jj)) = grid%nblist(1:grid%nbnr(jj),jj)
!
            if (mcel_cutoff.ge.1.0*grid%gpnbmd(i,j)) then ! add all without checks
              do nbi=1,grid%nbnr(i)
                rs1 = grid%nblist(nbi,i)
                rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + grid%nbnr(jj)
                if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
                  call nbl_resz(rs1,7)
                end if
                rs_nbl(rs1)%tmpanb((rs_nbl(rs1)%ntmpanb-grid%nbnr(jj)+1):rs_nbl(rs1)%ntmpanb) = lrs1(1:grid%nbnr(jj))
              end do
            else
              frs1(1:grid%nbnr(jj)) = 1.0*grid%gpnbmd(i,j) - resrad(lrs1(1:grid%nbnr(jj)))
              do nbi=1,grid%nbnr(i)
                rs1 = grid%nblist(nbi,i)
                frs1(1:grid%nbnr(jj)) = frs1(1:grid%nbnr(jj)) - resrad(rs1)
                rs2 = 0
                do nbj=1,grid%nbnr(jj)
                  if (frs1(nbj).le.mcel_cutoff) then
                    rs2 = rs2 + 1
                    lrs2(rs2) = lrs1(nbj)
                  end if
                end do
                rs_nbl(rs1)%ntmpanb = rs_nbl(rs1)%ntmpanb + rs2
                if (rs_nbl(rs1)%ntmpanb.gt.rs_nbl(rs1)%tmpalsz) then
                  call nbl_resz(rs1,7)
                end if
                rs_nbl(rs1)%tmpanb((rs_nbl(rs1)%ntmpanb-rs2+1):rs_nbl(rs1)%ntmpanb) = lrs2(1:rs2)
              end do
            end if
          end do
        end do

      end if ! if cycle_tab is true
!
    else if (use_rescrit.EQV..true.) then
!
!     do nothing
!
    else
!
      write(ilog,*) 'Fatal. Encountered unsupported cutoff treatment in all_respairs_nbl(...).&
 & Please report this bug.'
      call fexit()
    end if
!
  else
!
    write(ilog,*) 'Fatal. Called all_respairs_nbl(...) with cutoffs turned off.'
    call fexit()
!
  end if
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (thr_dlb(1,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(2,tpi) = thr_timings(2,tpi) + ttimer
    end if
  end if
#endif
  do rs=sta2,sto2
!   initialize
    rsp_vec(rs) = 0
    rs_nbl(rs)%nnbs = 0
    rs_nbl(rs)%nnbtrs = 0
    rs_nbl(rs)%nnblrs = 0
    rs_nbl(rs)%nnbats = 0
    rs_nbl(rs)%nnbtrats = 0
    rs_nbl(rs)%nnblrats = 0
  end do
  if (use_FEG.EQV..true.) then
    do rs=sta2,sto2
      rs_nbl(rs)%ngnbs = 0
      rs_nbl(rs)%ngnbtrs = 0
      rs_nbl(rs)%ngnbats = 0
      rs_nbl(rs)%ngnbtrats = 0
    end do
  end if
  if (use_waterloops.EQV..true.) then
    do rs=sta2,sto2
      rs_nbl(rs)%nwnbs = 0
      rs_nbl(rs)%nwnbtrs = 0
      rs_nbl(rs)%nwnbats = 0
      rs_nbl(rs)%nwnbtrats = 0
!      rs_nbl(rs)%nwnblrs = 0
!      rs_nbl(rs)%nwnblrats = 0
    end do
  end if
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
!$OMP BARRIER
    if (thr_dlb(2,1).gt.0) then
      if (tpi.eq.1) thr_dlb(2,2) = thr_dlb(2,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(3,tpi) = thr_timings(3,tpi) + ttimer
    end if
  end if
#endif
  if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
    do rs=sta2,sto2
      rs2 = rs_nbl(rs)%ntmpanb
      if (rs2.gt.0) call setup_nblidx(rs,cycle_frz,cycle_tab,atwo,tpi,rs2)
    end do
    if ((lrel_md.eq.4).OR.(lrel_md.eq.5)) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call complete_nbl_threads(cycle_frz,tpi)
      else
#endif
      do rs=sta2,sto2
        call complete_nbl(rs,cycle_frz)
      end do
#ifdef ENABLE_THREADS
      end if 
#endif
    end if
  else if (cycle_tab.EQV..true.) then
    do rs=sta2,sto2
      rs2 = rs_nbl(rs)%ntabnbs
      if (rs2.gt.0) call setup_nblidx(rs,cycle_frz,cycle_tab,aone,tpi,rs2)
    end do
  else
    do rs=sta2,sto2
      rs2 = nseq-rs
      if (rs2.gt.0) call setup_nblidx(rs,cycle_frz,cycle_tab,aone,tpi,rs2)
    end do
  end if
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    if (thr_dlb(2,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(4,tpi) = thr_timings(4,tpi) + ttimer
    end if
!$OMP BARRIER
  end if
#endif
!
end
!
!-----------------------------------------------------------------------------
!
! this is the MC-relevant "neighbor-list" routine, which simply populates
! the residue-by-residue matrix rsp_mat(i,j) for i > j SR terms, and for j > i
! LR terms
!
subroutine respairs(rs)
!
  use iounit
  use sequen
  use energies
  use cutoffs
#if ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer i,ii,kk,rs,sta,sto
  RTYPE dis,dis2
#ifdef ENABLE_THREADS
  integer tpi,tpn,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
!
  kk = refat(rs)
!
#if ENABLE_THREADS
  call omp_set_num_threads(min(nseq,thrdat%maxn))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tpn,tpi,sta,sto,i,ii,dis,dis2) FIRSTPRIVATE(kk)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  sto = tpi*nseq/tpn
  sta = (tpi-1)*nseq/tpn + 1
#else
  sta=1
  sto=nseq
#endif
  do i=sta,sto
    ii = refat(i)
    call dis_bound(ii,kk,dis2)
    dis = sqrt(dis2)
!   now check for the short-range cutoff
    if (dis.lt.(mcnb_cutoff+resrad(rs)+resrad(i))) then
      if (i.ge.rs) then
        rsp_mat(i,rs) = 1
      else
        rsp_mat(rs,i) = 1
      end if
    end if
    if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) cycle
!   within standard "twin"-range second cutoff (TABUL, POLAR)
    if (dis.lt.(mcel_cutoff+resrad(rs)+resrad(i))) then
      if (i.ge.rs) then
        rsp_mat(rs,i) = 1
      else
        rsp_mat(i,rs) = 1
      end if 
    else if (use_POLAR.EQV..true.) then
!     depending on LR electrostatics, turn on extra terms beyond that
      if (lrel_mc.eq.1) then
        if ((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.)) then
          if (i.ge.rs) then
            rsp_mat(rs,i) = 2
          else
            rsp_mat(i,rs) = 2
          end if
        end if
      else if ((lrel_mc.eq.2).OR.(lrel_mc.eq.3)) then
        if ((chgflag(rs).EQV..true.).AND.(chgflag(i).EQV..true.)) then
          if (i.ge.rs) then
            rsp_mat(rs,i) = 2
          else
            rsp_mat(i,rs) = 2
          end if
        end if
!     do nothing for lrel_mc = 4
      end if
    end if
  end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP END PARALLEL
#endif
!
end
!
!-----------------------------------------------------------------------------
!
! this neighbor list routine is only relevant for MC and populates an explicit list
! of pair interactions (per thread if so desired).
! this differs from the MD approach in that there are only single lists containing pair
! information rather than N lists containing neighbor information.
!
! the stretch irs-frs is the one supposed to have moved, with xrs the residue undergoing internal
! rearrangement: if frzflg is 0 this is marginal as all unique interactions are considered;
! with frzflg 1, only interactions between moving and nonmoving parts are considered + interactions
! involving xrs
!
! the stretch rsl-rsh is the one to be scanned.
!
! there are always two alternative lists selected by ix.
!
! in case of grid-based cutoffs, we pull the candidate list per residue from the lattice assignment.
! this causes complications with settings for lrel_mc not being 4, which do require augmentation.
! this is solved either by explicit accounting/adding via cglst (complicated due to having to track whether
! a pair has been added on account of distance already) or by override (for lrel_mc being 1, a charged
! residue requires all interactions anyway leading to simply forgetting the lattice information (bujcg is 2))
!
subroutine respairs_new(rsl,rsh,ix,irs,frs,xrs,initflag,frzflg,tpi)
!
  use iounit
  use sequen
  use energies
  use cutoffs
  use mcsums
  use polypep
  use inter, ONLY: nrsnb,nrsintra,nrpolnb,nrpolintra
  use mcgrid
  use tabpot
#if ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: rsl,rsh,tpi,ix,irs,frs,xrs,frzflg
  logical, INTENT(IN):: initflag
!
  integer aone,atwo,athree
  integer i,ii,iii,kk,sta,sto,rs,gpi,gpil,gpnx,bujcg
  integer(KIND=8):: ttt1,ttt2
  RTYPE dis2,hlp1,hlp2,hlp3,hlp4
  logical isinmov,isinset,rsisinmov,need_corrs,have_lr,mode2
#ifdef ENABLE_THREADS
  integer incr,ihlp,tpx
!
  if (tpi.gt.0) then
    sta = thr_limits(3,tpi)
    sto = thr_limits(4,tpi)
  else
    sta = 1
    sto = nseq
  end if
  if (tpi.le.1) call System_Clock(ttt1)
#else
!
  call System_Clock(ttt1)
  sta = 1
  sto = nseq
#endif
!
  have_lr = .true.
  if ((use_POLAR.EQV..false.).AND.(use_TABUL.EQV..false.)) have_lr = .false.
!
  if (ix.eq.2) then
#ifdef ENABLE_THREADS
    tpx = thrdat%maxn + tpi
#endif
  else
#ifdef ENABLE_THREADS
    tpx = tpi
#endif
  end if
  aone = 1
  atwo = 2
  athree = 3
  if (initflag.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      thr_rsp_nbl(tpx)%srnrs = 0
      thr_rsp_nbl(tpx)%trnrs = 0
      thr_rsp_nbl(tpx)%lrnrs = 0
      thr_rsp_nbl(tpx)%sr(1,3) = 0
      thr_rsp_nbl(tpx)%tr(1,3) = 0
      thr_rsp_nbl(tpx)%lr(1,3) = 0
    else
#endif
    rsp_nbl(ix)%srnrs = 0
    rsp_nbl(ix)%trnrs = 0
    rsp_nbl(ix)%lrnrs = 0
#ifdef ENABLE_THREADS
    end if
#endif
  end if
!
! with topology-assisted cutoffs, we have to consider all relevant pairs regardless -> simple 
!
  if (use_mcgrid.EQV..false.) then
!
    do rs=rsl,rsh
      kk = refat(rs)
      hlp1 = mcnb_cutoff+resrad(rs)
      hlp2 = mcel_cutoff+resrad(rs)
      rsisinmov = .false.
      if ((rs.ge.irs).AND.(rs.le.frs)) rsisinmov = .true.
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        incr = thrdat%maxn
        sta = tpi
        sto = nseq
      else
        incr = 1
        sta = 1
        sto = nseq
      end if
      do iii=sta,sto,incr
#else
      sta = 1
      sto = nseq
      do iii=sta,sto
#endif
        i = iii
        isinset = .false.
        isinmov = .false.
        if ((i.ge.rsl).AND.(i.le.rsh)) isinset = .true.
        if ((i.ge.irs).AND.(i.le.frs)) isinmov = .true.
!       don't double-count intra-set terms
        if ((isinset.EQV..true.).AND.(rs.gt.i)) cycle
!       don't count terms within a rigid moving part (frzflag 1 or 2) or at least don't double count (0)
        if ((isinmov.EQV..true.).AND.(rsisinmov.EQV..true.)) then
          if ((frzflg.eq.0).AND.(rs.gt.i)) cycle
          if ((frzflg.eq.1).AND.(rs.ne.xrs).AND.(i.ne.xrs)) cycle
        end if
!
        ii = refat(i)
        call dis_bound(ii,kk,dis2)
!
        hlp3 = (hlp1 + resrad(i))**2
        if (dis2.lt.hlp3) then
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
            thr_rsp_nbl(tpx)%srnrs = thr_rsp_nbl(tpx)%srnrs + 1
            if (thr_rsp_nbl(tpx)%srnrs.gt.thr_rsp_nbl(tpx)%sralsz) call thr_rsp_nbl_resz(tpi,aone,ix)
            thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,1) = min(rs,i)
            thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,2) = max(rs,i)
            if (abs(rs-i).gt.1) then
              thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = at(rs)%na*at(i)%na
            else
              if (rs.eq.i) then
                thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = max(1,nrsintra(rs))
              else
                thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = max(1,nrsnb(min(rs,i)))
              end if
            end if
            if (thr_rsp_nbl(tpx)%srnrs.gt.1) thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = &
 &         thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) + thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs-1,3)
          else
#endif
          rsp_nbl(ix)%srnrs = rsp_nbl(ix)%srnrs + 1
          if (rsp_nbl(ix)%srnrs.gt.rsp_nbl(ix)%sralsz) call rsp_nbl_resz(aone,ix)
          rsp_nbl(ix)%sr(rsp_nbl(ix)%srnrs,1) = min(rs,i)
          rsp_nbl(ix)%sr(rsp_nbl(ix)%srnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
          end if
#endif
        end if
!
        if (have_lr.EQV..false.) cycle
!        if ((use_ionloops.EQV..true.).AND.((rs.ge.rsmion1).OR.(i.ge.rsmion1))) cycle
!
        hlp4 = (hlp2 + resrad(i))**2
        if (dis2.lt.hlp4) then
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
            if (use_POLAR.EQV..true.) then
              thr_rsp_nbl(tpx)%trnrs = thr_rsp_nbl(tpx)%trnrs + 1
              if (thr_rsp_nbl(tpx)%trnrs.gt.thr_rsp_nbl(tpx)%tralsz) call thr_rsp_nbl_resz(tpi,atwo,ix)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,1) = min(rs,i)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,2) = max(rs,i)
              if (abs(rs-i).gt.1) then
                thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,at(rs)%npol*at(i)%npol)
              else
                if (rs.eq.i) then
                  thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,nrpolintra(rs))
                else
                  thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,nrpolnb(min(rs,i)))
                end if
              end if
              if (use_TABUL.EQV..true.) then
                ihlp = tbp%rsmat(max(i,rs),min(i,rs)) - tbp%rsmat(min(i,rs),max(i,rs)) + 1
                thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) + max(1,ihlp)
              end if
            else if (tbp%rsmat(max(i,rs),min(i,rs)).gt.0) then ! TABUL must be true in this branch
              thr_rsp_nbl(tpx)%trnrs = thr_rsp_nbl(tpx)%trnrs + 1
              if (thr_rsp_nbl(tpx)%trnrs.gt.thr_rsp_nbl(tpx)%tralsz) call thr_rsp_nbl_resz(tpi,atwo,ix)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,1) = min(rs,i)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,2) = max(rs,i)
              ihlp = tbp%rsmat(max(i,rs),min(i,rs)) - tbp%rsmat(min(i,rs),max(i,rs)) + 1
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,ihlp)
            end if
            if (thr_rsp_nbl(tpx)%trnrs.gt.1) thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = &
 &         thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) + thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs-1,3)
          else
#endif
          if (use_POLAR.EQV..true.) then
            rsp_nbl(ix)%trnrs = rsp_nbl(ix)%trnrs + 1
            if (rsp_nbl(ix)%trnrs.gt.rsp_nbl(ix)%tralsz) call rsp_nbl_resz(atwo,ix)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,1) = min(rs,i)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,2) = max(rs,i)
          else if (tbp%rsmat(max(i,rs),min(i,rs)).gt.0) then
            rsp_nbl(ix)%trnrs = rsp_nbl(ix)%trnrs + 1
            if (rsp_nbl(ix)%trnrs.gt.rsp_nbl(ix)%tralsz) call rsp_nbl_resz(atwo,ix)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,1) = min(rs,i)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,2) = max(rs,i)
          end if
#ifdef ENABLE_THREADS
          end if
#endif
!
        else if ((use_POLAR.EQV..true.).AND.(lrel_mc.ne.4)) then
          if (((lrel_mc.eq.1).AND.((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.))).OR.&
 &             ((chgflag(rs).EQV..true.).AND.(chgflag(i).EQV..true.))) then
#ifdef ENABLE_THREADS
            if (tpi.gt.0) then
              thr_rsp_nbl(tpx)%lrnrs = thr_rsp_nbl(tpx)%lrnrs + 1
              if (thr_rsp_nbl(tpx)%lrnrs.gt.thr_rsp_nbl(tpx)%lralsz) call thr_rsp_nbl_resz(tpi,athree,ix)
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,1) = min(rs,i)
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,2) = max(rs,i)
              if (lrel_mc.le.2) then
                thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = max(1,at(rs)%npol*at(i)%npol)
              else if (lrel_mc.eq.3) then
                thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = 1
              end if
              if (thr_rsp_nbl(tpx)%lrnrs.gt.1) thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = &
 &           thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) + thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs-1,3)
            else
#endif
            rsp_nbl(ix)%lrnrs = rsp_nbl(ix)%lrnrs + 1
            if (rsp_nbl(ix)%lrnrs.gt.rsp_nbl(ix)%lralsz) call rsp_nbl_resz(athree,ix)
            rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,1) = min(rs,i)
            rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
            end if
#endif
          end if
        end if
      end do
    end do
!
! with grid-assisted cutoffs, we pull a grid-based candidate list
! the problem is that additional interactions necessitated by lrel_mc 1-3 still have to be accounted for without
! distance evaluation
!
  else
!
    need_corrs = .false.
    if ((use_POLAR.EQV..true.).AND.(lrel_mc.ne.4).AND.(cglst%ncs.gt.0)) need_corrs = .true.
!
    if (need_corrs.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
     rsp_lmat(:,:) = .false.
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    end if
    gpil = 0
    do rs=rsl,rsh
      kk = refat(rs)
      hlp1 = mcnb_cutoff+resrad(rs)
      hlp2 = mcel_cutoff+resrad(rs)
      rsisinmov = .false.
      if ((rs.ge.irs).AND.(rs.le.frs)) rsisinmov = .true.
      gpi = grid%resgp(rs)
      if (gpi.ne.gpil) then
        sta = 1
        sto = 0
#ifdef ENABLE_THREADS
        if (tpi.gt.0) then
          if (grid%nbnr(gpi).gt.0) thr_rsp_nbl(tpi)%tmpl(1:grid%nbnr(gpi)) = grid%nblist(1:grid%nbnr(gpi),gpi)
          sto = grid%nbnr(gpi)
          do i=1,grid%gpnbnr(gpi)
            gpnx = grid%gpnblist(gpi,i)
            if (grid%nbnr(gpnx).gt.0) then
              thr_rsp_nbl(tpi)%tmpl((sto+1):(sto+grid%nbnr(gpnx))) = grid%nblist(1:grid%nbnr(gpnx),gpnx)
              sto = sto + grid%nbnr(gpnx)
            end if
          end do
          do i=1,grid%b_gpnbnr(gpi)
            gpnx = grid%b_gpnblist(gpi,i)
            if (grid%nbnr(gpnx).gt.0) then
              thr_rsp_nbl(tpi)%tmpl((sto+1):(sto+grid%nbnr(gpnx))) = grid%nblist(1:grid%nbnr(gpnx),gpnx)
              sto = sto + grid%nbnr(gpnx)
            end if
          end do
        else
#endif
        if (grid%nbnr(gpi).gt.0) grid%tmpl(1:grid%nbnr(gpi),1) = grid%nblist(1:grid%nbnr(gpi),gpi)
        sto = grid%nbnr(gpi)
        do i=1,grid%gpnbnr(gpi)
          gpnx = grid%gpnblist(gpi,i)
          if (grid%nbnr(gpnx).gt.0) then
            grid%tmpl((sto+1):(sto+grid%nbnr(gpnx)),1) = grid%nblist(1:grid%nbnr(gpnx),gpnx)
            sto = sto + grid%nbnr(gpnx)
          end if
        end do
        do i=1,grid%b_gpnbnr(gpi)
          gpnx = grid%b_gpnblist(gpi,i)
          if (grid%nbnr(gpnx).gt.0) then
            grid%tmpl((sto+1):(sto+grid%nbnr(gpnx)),1) = grid%nblist(1:grid%nbnr(gpnx),gpnx)
            sto = sto + grid%nbnr(gpnx)
          end if
        end do
        grid%ntmpl = sto
        if (grid%ntmpl.gt.nseq) call fexit()
#ifdef ENABLE_THREADS
        end if
#endif
        gpil = gpi
      end if
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        incr = thrdat%maxn
        sta = tpi
      else
        sta = 1
        sto = grid%ntmpl
        incr = 1
      end if
      do iii=sta,sto,incr
        if (tpi.gt.0) then
          i = thr_rsp_nbl(tpi)%tmpl(iii)
        else
          i = grid%tmpl(iii,1)
        end if
#else
      sta = 1
      sto = grid%ntmpl
      do iii=sta,sto
        i = grid%tmpl(iii,1)
#endif
        isinset = .false.
        isinmov = .false.
        if ((i.ge.rsl).AND.(i.le.rsh)) isinset = .true.
        if ((i.ge.irs).AND.(i.le.frs)) isinmov = .true.
!       don't double-count intra-set terms
        if ((isinset.EQV..true.).AND.(rs.gt.i)) cycle
!       don't count terms within a rigid moving part (frzflag 1 or 2) or at least don't double count (0)
        if ((isinmov.EQV..true.).AND.(rsisinmov.EQV..true.)) then
          if ((frzflg.eq.0).AND.(rs.gt.i)) cycle
          if ((frzflg.eq.1).AND.(rs.ne.xrs).AND.(i.ne.xrs)) cycle
        end if
        if (need_corrs.EQV..true.) then
          if (lrel_mc.eq.1) then
            if (chgflag(i).EQV..true.) rsp_lmat(cglst%irsl(i),rs) = .true.
            if (chgflag(rs).EQV..true.) rsp_lmat(cglst%irsl(rs),i) = .true.
          else
            if ((chgflag(i).EQV..true.).AND.(chgflag(rs).EQV..true.)) then
              rsp_lmat(cglst%irsl(rs),cglst%irsl(i)) = .true.
              rsp_lmat(cglst%irsl(i),cglst%irsl(rs)) = .true.
            end if
          end if
        end if
!
        ii = refat(i)
        call dis_bound(ii,kk,dis2)
!
        hlp3 = (hlp1 + resrad(i))**2
        if (dis2.lt.hlp3) then
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
            thr_rsp_nbl(tpx)%srnrs = thr_rsp_nbl(tpx)%srnrs + 1
            if (thr_rsp_nbl(tpx)%srnrs.gt.thr_rsp_nbl(tpx)%sralsz) call thr_rsp_nbl_resz(tpi,aone,ix)
            thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,1) = min(rs,i)
            thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,2) = max(rs,i)
            if (abs(rs-i).gt.1) then
              thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = at(rs)%na*at(i)%na
            else
              if (rs.eq.i) then
                thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = max(1,nrsintra(rs))
              else
                thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = max(1,nrsnb(min(rs,i)))
              end if
            end if
            if (thr_rsp_nbl(tpx)%srnrs.gt.1) thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) = &
 &         thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs,3) + thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%srnrs-1,3)
          else
#endif
          rsp_nbl(ix)%srnrs = rsp_nbl(ix)%srnrs + 1
          if (rsp_nbl(ix)%srnrs.gt.rsp_nbl(ix)%sralsz) call rsp_nbl_resz(aone,ix)
          rsp_nbl(ix)%sr(rsp_nbl(ix)%srnrs,1) = min(rs,i)
          rsp_nbl(ix)%sr(rsp_nbl(ix)%srnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
          end if
#endif
        end if
!
        if (have_lr.EQV..false.) cycle
!        if ((use_ionloops.EQV..true.).AND.((rs.ge.rsmion1).OR.(i.ge.rsmion1))) cycle
!
        hlp4 = (hlp2 + resrad(i))**2
        if (dis2.lt.hlp4) then
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
            if (use_POLAR.EQV..true.) then
              thr_rsp_nbl(tpx)%trnrs = thr_rsp_nbl(tpx)%trnrs + 1
              if (thr_rsp_nbl(tpx)%trnrs.gt.thr_rsp_nbl(tpx)%tralsz) call thr_rsp_nbl_resz(tpi,atwo,ix)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,1) = min(rs,i)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,2) = max(rs,i)
              if (abs(rs-i).gt.1) then
                thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,at(rs)%npol*at(i)%npol)
              else
                if (rs.eq.i) then
                  thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,nrpolintra(rs))
                else
                  thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,nrpolnb(min(rs,i)))
                end if
              end if
              if (use_TABUL.EQV..true.) then
                ihlp = tbp%rsmat(max(i,rs),min(i,rs)) - tbp%rsmat(min(i,rs),max(i,rs)) + 1
                thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) + max(1,ihlp)
              end if
            else if (tbp%rsmat(max(i,rs),min(i,rs)).gt.0) then ! TABUL must be true in this branch
              thr_rsp_nbl(tpx)%trnrs = thr_rsp_nbl(tpx)%trnrs + 1
              if (thr_rsp_nbl(tpx)%trnrs.gt.thr_rsp_nbl(tpx)%tralsz) call thr_rsp_nbl_resz(tpi,atwo,ix)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,1) = min(rs,i)
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,2) = max(rs,i)
              ihlp = tbp%rsmat(max(i,rs),min(i,rs)) - tbp%rsmat(min(i,rs),max(i,rs)) + 1
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,ihlp)
            end if
            if (thr_rsp_nbl(tpx)%trnrs.gt.1) thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = &
 &         thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) + thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs-1,3)
          else
#endif
          if (use_POLAR.EQV..true.) then
            rsp_nbl(ix)%trnrs = rsp_nbl(ix)%trnrs + 1
            if (rsp_nbl(ix)%trnrs.gt.rsp_nbl(ix)%tralsz) call rsp_nbl_resz(atwo,ix)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,1) = min(rs,i)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,2) = max(rs,i)
          else if (tbp%rsmat(max(i,rs),min(i,rs)).gt.0) then
            rsp_nbl(ix)%trnrs = rsp_nbl(ix)%trnrs + 1
            if (rsp_nbl(ix)%trnrs.gt.rsp_nbl(ix)%tralsz) call rsp_nbl_resz(atwo,ix)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,1) = min(rs,i)
            rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,2) = max(rs,i)
          end if
#ifdef ENABLE_THREADS
          end if
#endif
!
        else if ((use_POLAR.EQV..true.).AND.(lrel_mc.ne.4)) then
          if (((lrel_mc.eq.1).AND.((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.))).OR.&
 &             ((chgflag(rs).EQV..true.).AND.(chgflag(i).EQV..true.))) then
#ifdef ENABLE_THREADS
            if (tpi.gt.0) then
              thr_rsp_nbl(tpx)%lrnrs = thr_rsp_nbl(tpx)%lrnrs + 1
              if (thr_rsp_nbl(tpx)%lrnrs.gt.thr_rsp_nbl(tpx)%lralsz) call thr_rsp_nbl_resz(tpi,athree,ix)
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,1) = min(rs,i)
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,2) = max(rs,i)
              if (lrel_mc.le.2) then
                thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = max(1,at(rs)%npol*at(i)%npol)
              else if (lrel_mc.eq.3) then
                thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = 1
              end if
              if (thr_rsp_nbl(tpx)%lrnrs.gt.1) thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = &
 &           thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) + thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs-1,3)
            else
#endif
            rsp_nbl(ix)%lrnrs = rsp_nbl(ix)%lrnrs + 1
            if (rsp_nbl(ix)%lrnrs.gt.rsp_nbl(ix)%lralsz) call rsp_nbl_resz(athree,ix)
            rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,1) = min(rs,i)
            rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
            end if
#endif
          end if
        end if
      end do
    end do
!
    if (need_corrs.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      mode2 = .false.
      do rs=rsl,rsh
        if (lrel_mc.ne.1) then
          if (chgflag(rs).EQV..false.) cycle
        else
          if (chgflag(rs).EQV..false.) then
            mode2 = .false.
          else
            mode2 = .true.
          end if
        end if
        rsisinmov = .false.
        if ((rs.ge.irs).AND.(rs.le.frs)) rsisinmov = .true.
        bujcg = cglst%irsl(rs) ! 0 for neutral residue
#ifdef ENABLE_THREADS
        if (mode2.EQV..true.) then
          if (tpi.gt.0) then
            incr = thrdat%maxn
            sta = tpi
            sto = nseq
          else
            incr = 1
            sta = 1
            sto = nseq
          end if
        else
          if (tpi.gt.0) then
            incr = thrdat%maxn
            sta = tpi
            sto = cglst%ncrs
          else
            incr = 1
            sta = 1
            sto = cglst%ncrs
          end if
        end if
        do iii=sta,sto,incr
#else
        if (mode2.EQV..true.) then
          sta = 1
          sto = nseq
        else
          sta = 1
          sto = cglst%ncrs
        end if
        do iii=sta,sto
#endif 
          if (lrel_mc.eq.1) then
            if (mode2.EQV..true.) then
              if (rsp_lmat(bujcg,iii).EQV..true.) cycle
              i = iii
            else
              if (rsp_lmat(iii,rs).EQV..true.) cycle
              i = cglst%rsl(iii)
            end if
          else
            if (rsp_lmat(bujcg,iii).EQV..true.) cycle
            i = cglst%rsl(iii)
          end if
          if (lrel_mc.eq.1) then
            if ((chgflag(i).EQV..false.).AND.(chgflag(rs).EQV..false.)) cycle ! can happen with patches
          else
            if (chgflag(i).EQV..false.) cycle
          end if
          isinset = .false.
          isinmov = .false.
          if ((i.ge.rsl).AND.(i.le.rsh)) isinset = .true.
          if ((i.ge.irs).AND.(i.le.frs)) isinmov = .true.
!         don't double-count intra-set terms
          if ((isinset.EQV..true.).AND.(i.lt.rs)) cycle
!         don't count terms within a rigid moving part (frzflag 1 or 2) or at least don't double count (0)
          if ((isinmov.EQV..true.).AND.(rsisinmov.EQV..true.)) then
            if ((frzflg.eq.0).AND.(i.lt.rs)) cycle
            if ((frzflg.eq.1).AND.(rs.ne.xrs).AND.(i.ne.xrs)) cycle
          end if
!
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
            thr_rsp_nbl(tpx)%lrnrs = thr_rsp_nbl(tpx)%lrnrs + 1
            if (thr_rsp_nbl(tpx)%lrnrs.gt.thr_rsp_nbl(tpx)%lralsz) call thr_rsp_nbl_resz(tpi,athree,ix)
            thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,1) = min(rs,i)
            thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,2) = max(rs,i)
            if (lrel_mc.le.2) then
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = max(1,at(rs)%npol*at(i)%npol)
            else if (lrel_mc.eq.3) then
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = 1
            end if
            if (thr_rsp_nbl(tpx)%lrnrs.gt.1) thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = &
 &         thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) + thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs-1,3)
          else
#endif
          rsp_nbl(ix)%lrnrs = rsp_nbl(ix)%lrnrs + 1
          if (rsp_nbl(ix)%lrnrs.gt.rsp_nbl(ix)%lralsz) call rsp_nbl_resz(athree,ix)
          rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,1) = min(rs,i)
          rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
          end if
#endif
        end do
      end do
    end if
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if (tpi.le.1) then
!    if (tpi.eq.1) then
!      write(*,*) sum(thr_rsp_nbl(tpx:(tpx+thrdat%maxn-1))%srnrs),sum(thr_rsp_nbl(tpx:(tpx+thrdat%maxn-1))%trnrs),&
! &                 sum(thr_rsp_nbl(tpx:(tpx+thrdat%maxn-1))%lrnrs)
!    else
!      write(*,*) rsp_nbl(ix)%srnrs,rsp_nbl(ix)%trnrs,rsp_nbl(ix)%lrnrs
!    end if

!    if ((nstep.ge.1).AND.(tpi.le.1)) then
!      write(*,*) sum(thr_rsp_nbl(tpx:(tpx+thrdat%maxn-1))%srnrs),sum(thr_rsp_nbl(tpx:(tpx+thrdat%maxn-1))%trnrs),&
! &                 sum(thr_rsp_nbl(tpx:(tpx+thrdat%maxn-1))%lrnrs)
!    end if
    call System_Clock(ttt2)
    time_nbl = time_nbl + 1.0*(ttt2 - ttt1)
    time_energy = time_energy - 1.0*(ttt2 - ttt1)
  end if
!
end


!
!-----------------------------------------------------------------------------
! 
! a closely related routine that will augment existing neighbor lists with terms
! identified through rs_vec
!
! rs_vec-flagged residues are guaranteed to be in close proximity to the moving residues.
! this means that their putative interaction neighbors tend to have large overlap
! we do not, however, currently exploit this for grid-based cutoffs, which is due to the
! overhead of merging two pruned lists (the one with rs_vec nonzero is often even smaller)
!
! this means grid-based and topology-based cutoffs are not distinguishable in respairs_append
!
subroutine respairs_append(ix,irs,frs,xrs,stretchdone,tpi)
!
  use iounit
  use sequen
  use energies
  use cutoffs
  use mcsums
  use polypep
  use mcgrid
  use inter, ONLY: nrpolnb,nrpolintra
#if ENABLE_THREADS
  use threads
  use tabpot
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,ix,irs,frs,xrs,stretchdone
!
  integer aone,atwo,athree
  integer i,ii,kk,sta,sto,incr,rs
  integer(KIND=8):: ttt1,ttt2
  RTYPE dis2,hlp1,hlp2,hlp4
  logical isinmov,rsisinmov
#ifdef ENABLE_THREADS
  integer ihlp,tpx
!
  if ((use_IMPSOLV.EQV..false.).OR.(use_POLAR.EQV..false.)) return
!
  if (tpi.gt.0) then
    sta = thr_limits(3,tpi)
    sto = thr_limits(4,tpi)
  else
    sta = 1
    sto = nseq
  end if
  if (tpi.le.1) call System_Clock(ttt1)
#else
!
  call System_Clock(ttt1)
#endif
!
  if (ix.eq.2) then
#ifdef ENABLE_THREADS
    tpx = thrdat%maxn + tpi
#endif
  else
#ifdef ENABLE_THREADS
    tpx = tpi
#endif
  end if
  aone = 1
  atwo = 2
  athree = 3
!
  do rs=1,nseq
    if (rs.eq.xrs) cycle
    if (rs_vec(rs).le.0) cycle
    kk = refat(rs)
    hlp1 = mcnb_cutoff+resrad(rs)
    hlp2 = mcel_cutoff+resrad(rs)
    rsisinmov = .false.
    if ((rs.ge.irs).AND.(rs.le.frs)) rsisinmov = .true.
    if (rsisinmov.EQV..true.) then ! interactions between moving and non-moving parts are already covered
      sta = irs-1
      sto = frs
    else
      sta = 0
      sto = nseq
    end if
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      incr = thrdat%maxn
      sta = sta + tpi - incr
    else
      incr = 1
    end if
#else
    incr = 1
#endif
    i = sta
    do while ((i+incr).le.sto)
      i = i + incr
      isinmov = .false.
      if (i.eq.xrs) cycle
      if ((i.ge.irs).AND.(i.le.frs)) isinmov = .true.
      if ((isinmov.EQV..true.).AND.(rsisinmov.EQV..true.).AND.(stretchdone.eq.1)) cycle
      if ((rs.gt.i).AND.(rs_vec(i).gt.0)) cycle
      if (rsisinmov.EQV..false.) then
        if (isinmov.EQV..true.) then
!       MAKE MORE EFFICIENT
!          i = frs
#ifdef ENABLE_THREADS
!          if (tpi.gt.0) i = i + tpi - 1
#endif
          cycle
        end if
      end if
!     
      ii = refat(i)
      call dis_bound(ii,kk,dis2)
!
      hlp4 = (hlp2 + resrad(i))**2
      if (dis2.lt.hlp4) then
#ifdef ENABLE_THREADS
        if (tpi.gt.0) then
          thr_rsp_nbl(tpx)%trnrs = thr_rsp_nbl(tpx)%trnrs + 1
          if (thr_rsp_nbl(tpx)%trnrs.gt.thr_rsp_nbl(tpx)%tralsz) call thr_rsp_nbl_resz(tpi,atwo,ix)
          thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,1) = min(rs,i)
          thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,2) = max(rs,i)
          if (use_POLAR.EQV..true.) then
            if (abs(rs-i).gt.1) then
              thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,at(rs)%npol*at(i)%npol)
            else
              if (rs.eq.i) then
                thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,nrpolintra(rs))
              else
                thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = max(1,nrpolnb(min(rs,i)))
              end if
            end if
          end if
          if (use_TABUL.EQV..true.) then
            ihlp = tbp%rsmat(max(i,rs),min(i,rs)) - tbp%rsmat(min(i,rs),max(i,rs)) + 1
            thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) + max(1,ihlp)
          end if
          if (thr_rsp_nbl(tpx)%trnrs.gt.1) thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) = &
 &       thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs,3) + thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%trnrs-1,3)
        else
#endif
        rsp_nbl(ix)%trnrs = rsp_nbl(ix)%trnrs + 1
        if (rsp_nbl(ix)%trnrs.gt.rsp_nbl(ix)%tralsz) call rsp_nbl_resz(atwo,ix)
        rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,1) = min(rs,i)
        rsp_nbl(ix)%tr(rsp_nbl(ix)%trnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
        end if
#endif
      else if (lrel_mc.ne.4) then
        if (((lrel_mc.eq.1).AND.((chgflag(rs).EQV..true.).OR.(chgflag(i).EQV..true.))).OR.&
 &           ((chgflag(rs).EQV..true.).AND.(chgflag(i).EQV..true.))) then
#ifdef ENABLE_THREADS
          if (tpi.gt.0) then
            thr_rsp_nbl(tpx)%lrnrs = thr_rsp_nbl(tpx)%lrnrs + 1
            if (thr_rsp_nbl(tpx)%lrnrs.gt.thr_rsp_nbl(tpx)%lralsz) call thr_rsp_nbl_resz(tpi,athree,ix)
            thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,1) = min(rs,i)
            thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,2) = max(rs,i)
            if (lrel_mc.le.2) then
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = max(1,at(rs)%npol*at(i)%npol)
            else if (lrel_mc.eq.3) then
              thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = 1
            end if
            if (thr_rsp_nbl(tpx)%lrnrs.gt.1) thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) = &
 &         thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs,3) + thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lrnrs-1,3)
          else
#endif
          rsp_nbl(ix)%lrnrs = rsp_nbl(ix)%lrnrs + 1
          if (rsp_nbl(ix)%lrnrs.gt.rsp_nbl(ix)%lralsz) call rsp_nbl_resz(athree,ix)
          rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,1) = min(rs,i)
          rsp_nbl(ix)%lr(rsp_nbl(ix)%lrnrs,2) = max(rs,i)
#ifdef ENABLE_THREADS
          end if
#endif
        end if
      end if
    end do
  end do
!  write(*,*) irs,frs,rsp_nbl(ix)%lrnrs+rsp_nbl(ix)%trnrs,sum(rs_vec(1:nseq))
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if (tpi.le.1) then
    call System_Clock(ttt2)
    time_nbl = time_nbl + 1.0*(ttt2 - ttt1)
    time_energy = time_energy - 1.0*(ttt2 - ttt1)
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
end
!
!-----------------------------------------------------------------------------
!
! this extremely simple routine creates a later on static rs-rs neighbor list based
! on the presence of tabulated interactions (which may be sparse)
! note that this does not respect constraints (i.e, might be inefficient still)
!
  subroutine tab_respairs_nbl_pre(rs)
!
  use iounit
  use molecule
  use sequen
  use cutoffs
  use atoms
  use tabpot
!
  implicit none
!
  integer i,ii,rs
  integer, ALLOCATABLE:: inb(:)
!
  allocate(inb(nseq-rs))
!
  rs_nbl(rs)%ntabnbs = 0
  rs_nbl(rs)%ntabias = 0
!
  ii = refat(rs)
  do i=rs+1,nseq
    if (tbp%rsmat(rs,i).gt.0) then
      rs_nbl(rs)%ntabnbs = rs_nbl(rs)%ntabnbs + 1
      rs_nbl(rs)%ntabias = rs_nbl(rs)%ntabias + (tbp%rsmat(i,rs)-tbp%rsmat(rs,i)) + 1
      inb(rs_nbl(rs)%ntabnbs) = i
    end if
  end do
!  write(*,*) rs,rs_nbl(rs)%ntabnbs,rs_nbl(rs)%ntabias
!
  if (allocated(rs_nbl(rs)%tabnb).EQV..true.) then
    deallocate(rs_nbl(rs)%tabnb)
  end if
  if (rs_nbl(rs)%ntabnbs.gt.0) then
    allocate(rs_nbl(rs)%tabnb(rs_nbl(rs)%ntabnbs))
    rs_nbl(rs)%tabnb(1:rs_nbl(rs)%ntabnbs) = inb(1:rs_nbl(rs)%ntabnbs)
  end if
!
  deallocate(inb)
!
  end
!
!-----------------------------------------------------------------------------
!
! mode: 1 use all rsi > rs (assumes tmpanb not in use)
!       2 use those in rs_nbl(rs)%tmpanb
! cycle_frz: true exclude frozen ia.s from either
!            false do not exclude
! cycle_tab: true exclude non-tab-interacting ia.s from either
!            false do not exclude
!
subroutine setup_nblidx(rs,cycle_frz,cycle_tab,mode,tpi,talsz)
!
  use iounit
  use molecule
  use sequen
  use energies
  use cutoffs
  use polypep
  use system
  use atoms
  use tabpot
  use fyoc
!
  implicit none
!
  integer, INTENT(IN):: tpi,talsz,rs,mode
  logical, INTENT(IN):: cycle_frz,cycle_tab
!
  integer i,ii,j,hi,imol,jmol,tmpnb,sameresmol,hix(talsz)
  RTYPE lop1,lop3,thlper(talsz,6)
!
  imol = molofrs(rs)
  ii = refat(rs)
!
  if (mode.eq.1) then
!
    tmpnb = 0
!
    if (cycle_frz.EQV..true.) then
      if (cycle_tab.EQV..true.) then
        do j=1,rs_nbl(rs)%ntabnbs
          i = rs_nbl(rs)%tabnb(j)
          jmol = molofrs(i)
          if ((jmol.eq.imol).AND.(molfrzidx(imol).lt.3)) then
            if (molfrzidx(imol).eq.2) then
              if ((nchi(rs).eq.0).AND.(nchi(i).eq.0)) then
!               do nothing (backbone is frozen completely)
              else
                tmpnb = tmpnb + 1
                thlper(tmpnb,1) = x(refat(i)) - x(ii)
                thlper(tmpnb,2) = y(refat(i)) - y(ii)
                thlper(tmpnb,3) = z(refat(i)) - z(ii)
                hix(tmpnb) = i
              end if
            else
              tmpnb = tmpnb + 1
              thlper(tmpnb,1) = x(refat(i)) - x(ii)
              thlper(tmpnb,2) = y(refat(i)) - y(ii)
              thlper(tmpnb,3) = z(refat(i)) - z(ii)
              hix(tmpnb) = i
            end if
          else
            if ((molfrzidx(jmol).eq.4).AND.(molfrzidx(imol).eq.4)) then
!             do nothing (note there's no other way for the two molecules' interactions to never change)
            else
              tmpnb = tmpnb + 1
              thlper(tmpnb,1) = x(refat(i)) - x(ii)
              thlper(tmpnb,2) = y(refat(i)) - y(ii)
              thlper(tmpnb,3) = z(refat(i)) - z(ii)
              hix(tmpnb) = i
            end if
          end if
        end do
      else
        if (molfrzidx(imol).lt.3) then
          do i=rs+1,rsmol(imol,2)
            if (cycle_tab.EQV..true.) then
              if (tbp%rsmat(rs,i).le.0) cycle
            end if
            if (molfrzidx(imol).eq.2) then
              if ((nchi(rs).eq.0).AND.(nchi(i).eq.0)) then
!               do nothing (backbone is frozen completely)
              else
                tmpnb = tmpnb + 1
                thlper(tmpnb,1) = x(refat(i)) - x(ii)
                thlper(tmpnb,2) = y(refat(i)) - y(ii)
                thlper(tmpnb,3) = z(refat(i)) - z(ii)
                hix(tmpnb) = i
              end if
            else
              tmpnb = tmpnb + 1
              thlper(tmpnb,1) = x(refat(i)) - x(ii)
              thlper(tmpnb,2) = y(refat(i)) - y(ii)
              thlper(tmpnb,3) = z(refat(i)) - z(ii)
              hix(tmpnb) = i
            end if
          end do
        end if
        do jmol=imol+1,nmol
          if ((molfrzidx(jmol).eq.4).AND.(molfrzidx(imol).eq.4)) then
!           do nothing (note there's no other way for the two molecules' interactions to never change)
          else
            do i=rsmol(jmol,1),rsmol(jmol,2)
              if (cycle_tab.EQV..true.) then
                if (tbp%rsmat(rs,i).le.0) cycle
              end if
              tmpnb = tmpnb + 1
              thlper(tmpnb,1) = x(refat(i)) - x(ii)
              thlper(tmpnb,2) = y(refat(i)) - y(ii)
              thlper(tmpnb,3) = z(refat(i)) - z(ii)
              hix(tmpnb) = i
            end do
          end if
        end do
      end if
!
      hi = tmpnb
!
    else if (cycle_tab.EQV..true.) then ! but cycle_frz false
!
      tmpnb = 0
      do j=1,rs_nbl(rs)%ntabnbs
        i = rs_nbl(rs)%tabnb(j)
        if (tbp%rsmat(rs,i).le.0) cycle ! should never cycle
        tmpnb = tmpnb + 1
        thlper(tmpnb,1) = x(refat(i)) - x(ii)
        thlper(tmpnb,2) = y(refat(i)) - y(ii)
        thlper(tmpnb,3) = z(refat(i)) - z(ii)
        hix(tmpnb) = i
      end do
!
      hi = tmpnb
!
    else ! regular case
!
      hi = nseq-rs
      thlper(1:hi,1) = x(refat(rs+1:nseq)) - x(ii)
      thlper(1:hi,2) = y(refat(rs+1:nseq)) - y(ii)
      thlper(1:hi,3) = z(refat(rs+1:nseq)) - z(ii)
      do i=rs+1,nseq
        hix(i-rs) = i
      end do
!
    end if
!
  else if (mode.eq.2) then
!
    if (rs_nbl(rs)%ntmpanb.eq.0) return
    tmpnb = 0
!
    if (cycle_frz.EQV..true.) then
      do j=1,rs_nbl(rs)%ntmpanb
        i = rs_nbl(rs)%tmpanb(j)
        jmol = molofrs(i)
        if (cycle_tab.EQV..true.) then
          if (tbp%rsmat(rs,i).le.0) cycle
        end if
        if (imol.eq.jmol) then
          if (molfrzidx(imol).lt.3) then
            if (molfrzidx(imol).eq.2) then
!             the better approach here would be to scan only the segment between the two residues
!             for rigidity instead of relying in the unnecessarily stringent criterion in molfrzidx
              if ((nchi(rs).eq.0).AND.(nchi(i).eq.0)) then
!               do nothing (backbone is frozen completely)
              else
                tmpnb = tmpnb + 1
                thlper(tmpnb,1) = x(refat(i)) - x(ii)
                thlper(tmpnb,2) = y(refat(i)) - y(ii)
                thlper(tmpnb,3) = z(refat(i)) - z(ii)
                hix(tmpnb) = i
              end if
            else
              tmpnb = tmpnb + 1
              thlper(tmpnb,1) = x(refat(i)) - x(ii)
              thlper(tmpnb,2) = y(refat(i)) - y(ii)
              thlper(tmpnb,3) = z(refat(i)) - z(ii)
              hix(tmpnb) = i
            end if
          end if
        else      
          if ((molfrzidx(jmol).eq.4).AND.(molfrzidx(imol).eq.4)) then
!           do nothing (note there's no other way for the two molecules' interactions to never change)
          else
            tmpnb = tmpnb + 1
            thlper(tmpnb,1) = x(refat(i)) - x(ii)
            thlper(tmpnb,2) = y(refat(i)) - y(ii)
            thlper(tmpnb,3) = z(refat(i)) - z(ii)
            hix(tmpnb) = i
          end if
        end if
      end do
!
      hi = tmpnb
!
    else if (cycle_tab.EQV..true.) then ! but cycle_frz false (this may be inefficient)
!
      tmpnb = 0
!     note that the cutoff-aware pre-nb-list is also aware of tbp%rsmat
      do i=1,rs_nbl(rs)%ntmpanb
        if (tbp%rsmat(rs,rs_nbl(rs)%tmpanb(i)).le.0) cycle ! should never cycle
        tmpnb = tmpnb + 1
        thlper(tmpnb,1) = x(refat(rs_nbl(rs)%tmpanb(i))) - x(ii)
        thlper(tmpnb,2) = y(refat(rs_nbl(rs)%tmpanb(i))) - y(ii)
        thlper(tmpnb,3) = z(refat(rs_nbl(rs)%tmpanb(i))) - z(ii)
        hix(tmpnb) = rs_nbl(rs)%tmpanb(i)
      end do
      hi = tmpnb
!
    else ! regular case
!
      hi = rs_nbl(rs)%ntmpanb
      hix(1:hi) = refat(rs_nbl(rs)%tmpanb(1:hi))
      thlper(1:hi,1) = x(hix(1:hi)) - atinfo(1,ii)
      thlper(1:hi,2) = y(hix(1:hi)) - atinfo(2,ii)
      thlper(1:hi,3) = z(hix(1:hi)) - atinfo(3,ii)
      do i=1,hi
        hix(i) = rs_nbl(rs)%tmpanb(i)
      end do
!
    end if
!
  end if
!
! now take care of the shift vectors
!
! PBC
  if ((bnd_type.eq.1).AND.(hi.ge.1)) then
!   cubic box
    if (bnd_shape.eq.1) then
      do j=1,3
        lop1 = 1.0/bnd_params(j)
        lop3 = bnd_params(j)
        thlper(1:hi,j+3) = -lop3*anint(thlper(1:hi,j)*lop1)
      end do
    else if (bnd_shape.eq.3) then
      thlper(1:hi,4) = 0.0
      thlper(1:hi,5) = 0.0
      lop1 = 1.0/bnd_params(6)
      lop3 = bnd_params(6)
      thlper(1:hi,6) = -lop3*anint(thlper(1:hi,3)*lop1)
    else
      write(ilog,*) 'Fatal. Encountered unsupported box shape in res&
 &pairs_nbl() (code # ',bnd_shape,').'
      call fexit()
    end if
    sameresmol = rsmol(molofrs(rs),2)-rsmol(molofrs(rs),1)
    j = 0
    do i=1,hi
      if (j.ge.sameresmol) exit
      if (molofrs(hix(i)).eq.molofrs(rs)) then
        thlper(i,4:6) = 0.0
        j = j + 1
      end if
    end do
    thlper(1:hi,1) = sqrt((thlper(1:hi,1)+thlper(1:hi,4))**2 + (thlper(1:hi,2)+thlper(1:hi,5))**2 + &
 &                        (thlper(1:hi,3)+thlper(1:hi,6))**2 )
  else
    thlper(1:hi,4) = 0.0
    thlper(1:hi,5) = 0.0
    thlper(1:hi,6) = 0.0
    thlper(1:hi,1) = sqrt(thlper(1:hi,1)**2 + thlper(1:hi,2)**2 + thlper(1:hi,3)**2 )
  end if
!
!  hi2 = hi
!  if ((mode.eq.2).AND.((lrel_md.eq.4).OR.(lrel_md.eq.5))) then
!    ihix(:) = 0
!    do i=1,hi
!      ihix(hix(i)) = 1
!    end do
!    if (((lrel_md.eq.3).AND.(chargeflag(rs).EQV..false.)).OR.(lrel_md.eq.4)) then
!    else
!      do i=rs+2,nseq
!        jmol = molofrs(rs2)
!        if (ihix(i).eq.1) cycle
!        if (cycle_frz.EQV..true.) then
!          if (imol.eq.jmol) then
!            if (molfrzidx(imol).ge.3) cycle
!            if ((molfrzidx(imol).eq.2).AND.(nchi(rs).eq.0).AND.(nchi(i).eq.0)) cycle
!          else
!            if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
!          end if
!        end if
!        hi2 = hi2 + 1
!      end do
!  end if  
!
! now call the fxn which parses those arrays
  call populate_nbl(rs,hi,thlper(1:hi,4),thlper(1:hi,5),thlper(1:hi,6),thlper(1:hi,1),hix(1:hi),tpi)
!
  end
!
!-----------------------------------------------------------------------------
!
! the entries in the nb-list always have higher index
! advantages: list is safe to use with water loops
! disadvantage: lists are of uneven size
!
  subroutine populate_nbl(rs,hi,svecx,svecy,svecz,d1,idx,tpi)
!
  use iounit
  use molecule
  use sequen
  use energies
  use cutoffs
  use polypep
  use system
  use atoms
  use fyoc
!
  implicit none
!
  integer, INTENT(IN):: hi,rs,idx(hi),tpi
  RTYPE, INTENT(IN):: svecx(hi),svecy(hi),svecz(hi),d1(hi)
!
  integer i,ii,rs1,rs2,proxna
  RTYPE hlp1,hlp2,svsgn
  logical do_lr
!
  hlp1 = mcnb_cutoff+resrad(rs)
  hlp2 = mcel_cutoff+resrad(rs)
!
  if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.5).OR.(lrel_md.eq.4))) then
    do_lr = .true.
  else
    do_lr = .false.
  end if
  svsgn = 1.0
!
  if (use_FEG.EQV..false.) then
    do ii=1,hi
      i = idx(ii)
      if (tpi.eq.0) then
        rs1 = min(rs,i)
        rs2 = max(rs,i)
      else
        rs1 = rs ! min(rs,i)
        rs2 = i ! max(rs,i)
      end if
      if (do_lr.EQV..true.) then
        if (lrel_md.eq.4) then
          proxna = at(rs2)%nncgrps
        else
          proxna = at(rs2)%npol
        end if
      end if
      if ((use_waterloops.EQV..true.).AND.(rs2.ge.rsw1)) then
        if (d1(ii).lt.(hlp1+resrad(i))) then
          rs_nbl(rs1)%nwnbs = rs_nbl(rs1)%nwnbs + 1
          if (rs_nbl(rs1)%nwnbs.gt.rs_nbl(rs1)%wnbalsz) call nbl_resz(rs1,8)
          rs_nbl(rs1)%nwnbats = rs_nbl(rs1)%nwnbats + at(rs2)%na
          rs_nbl(rs1)%wnb(rs_nbl(rs1)%nwnbs) = rs2
        else if (d1(ii).lt.(hlp2+resrad(i))) then
          if (rs1.eq.rs) then
            svsgn = 1.0
          else
            svsgn = -1.0
          end if
          rs_nbl(rs1)%nwnbtrs = rs_nbl(rs1)%nwnbtrs + 1
          if (rs_nbl(rs1)%nwnbtrs.gt.rs_nbl(rs1)%wtralsz) call nbl_resz(rs1,9)
          rs_nbl(rs1)%nwnbtrats = rs_nbl(rs1)%nwnbtrats + at(rs2)%na
          rs_nbl(rs1)%wnbtr(rs_nbl(rs1)%nwnbtrs) = rs2
          if (molofrs(rs2).eq.molofrs(rs1)) then
            rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,:) = 0.0
          else
            rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,1) = svsgn*svecx(ii)
            rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,2) = svsgn*svecy(ii)
            rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,3) = svsgn*svecz(ii)
          end if
        else
          if (do_lr.EQV..true.) then
            if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!              if (rs2.eq.(rs1+1)) then
!              else if (rs1.eq.(rs2+1)) then
!               do nothing (detected as rsp_vec(i) being 0)
!              else
                if ((molofrs(rs1).eq.molofrs(rs2)).AND.(abs(rs2-rs1).eq.1)) call fexit()
                if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
                rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
                if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
                rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
                rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
!              end if
            end if
          end if
        end if
      else
        if ((rs1.eq.disulf(rs2)).OR.(rs2.eq.disulf(rs1))) cycle ! exclude crosslinked pairs
        if (d1(ii).lt.(hlp1+resrad(i))) then
          if (rs2.eq.(rs1+1)) then
            rsp_vec(rs1) = 1
          else if (rs1.eq.(rs2+1)) then
            rsp_vec(rs2) = 1
          else
            rs_nbl(rs1)%nnbs = rs_nbl(rs1)%nnbs + 1
            if (rs_nbl(rs1)%nnbs.gt.rs_nbl(rs1)%nbalsz) call nbl_resz(rs1,1)
            rs_nbl(rs1)%nnbats = rs_nbl(rs1)%nnbats + at(rs2)%na
            rs_nbl(rs1)%nb(rs_nbl(rs1)%nnbs) = rs2
          end if
        else if (d1(ii).lt.(hlp2+resrad(i))) then
          if (rs2.eq.(rs1+1)) then
            rsp_vec(rs1) = 2
           else if (rs1.eq.(rs2+1)) then
            rsp_vec(rs2) = 2
          else
            if (rs1.eq.rs) then
              svsgn = 1.0
            else
              svsgn = -1.0
            end if
            rs_nbl(rs1)%nnbtrs = rs_nbl(rs1)%nnbtrs + 1
            if (rs_nbl(rs1)%nnbtrs.gt.rs_nbl(rs1)%tralsz) call nbl_resz(rs1,2)
            rs_nbl(rs1)%nnbtrats = rs_nbl(rs1)%nnbtrats + at(rs2)%na
            rs_nbl(rs1)%nbtr(rs_nbl(rs1)%nnbtrs) = rs2
            if (molofrs(rs2).eq.molofrs(rs1)) then
              rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,:) = 0.0
            else
              rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,1) = svsgn*svecx(ii)
              rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,2) = svsgn*svecy(ii)
              rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,3) = svsgn*svecz(ii)
            end if
          end if
        else
          if (do_lr.EQV..true.) then
            if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!              if (rs2.eq.(rs1+1)) then
!              else if (rs1.eq.(rs2+1)) then
!               do nothing (detected as rsp_vec(i) being 0)
!              else
                if ((molofrs(rs1).eq.molofrs(rs2)).AND.(abs(rs2-rs1).eq.1)) call fexit()
                if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
                rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
                if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
                rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
                rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
!              end if
            end if
          end if
        end if
      end if
    end do
!
  else
!
!   for calculations with ghosts we have (annoyingly) separate nb-lists
    do ii=1,hi
      i = idx(ii)
      if (tpi.eq.0) then
        rs1 = min(rs,i)
        rs2 = max(rs,i)
      else
        rs1 = rs ! min(rs,i)
        rs2 = i ! max(rs,i)
      end if
      if (do_lr.EQV..true.) then
        if (lrel_md.eq.4) then
          proxna = at(rs2)%nncgrps
        else
          proxna = at(rs2)%npol
        end if
      end if
      if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!       do nothing (at least one of the residues is fully de-coupled)
      else if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
 &    ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
!       this is the case using the full Hamiltonian, same code as above
        if ((use_waterloops.EQV..true.).AND.(rs2.ge.rsw1)) then
          if (d1(ii).lt.(hlp1+resrad(i))) then
            rs_nbl(rs1)%nwnbs = rs_nbl(rs1)%nwnbs + 1
            if (rs_nbl(rs1)%nwnbs.gt.rs_nbl(rs1)%wnbalsz) call nbl_resz(rs1,8)
            rs_nbl(rs1)%nwnbats = rs_nbl(rs1)%nwnbats + at(rs2)%na
            rs_nbl(rs1)%wnb(rs_nbl(rs1)%nwnbs) = rs2
          else if (d1(ii).lt.(hlp2+resrad(i))) then
            if (rs1.eq.rs) then
              svsgn = 1.0
            else
              svsgn = -1.0
            end if
            rs_nbl(rs1)%nwnbtrs = rs_nbl(rs1)%nwnbtrs + 1
            if (rs_nbl(rs1)%nwnbtrs.gt.rs_nbl(rs1)%wtralsz) call nbl_resz(rs1,9)
            rs_nbl(rs1)%nwnbtrats = rs_nbl(rs1)%nwnbtrats + at(rs2)%na
            rs_nbl(rs1)%wnbtr(rs_nbl(rs1)%nwnbtrs) = rs2
            if (molofrs(rs2).eq.molofrs(rs1)) then
              rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,:) = 0.0
            else
              rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,1) = svsgn*svecx(ii)
              rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,2) = svsgn*svecy(ii)
              rs_nbl(rs1)%wtrsvec(rs_nbl(rs1)%nwnbtrs,3) = svsgn*svecz(ii)
            end if
          else
            if (do_lr.EQV..true.) then
              if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!                if (rs2.eq.(rs1+1)) then
!                else if (rs1.eq.(rs2+1)) then
!                 do nothing (detected as rsp_vec(i) being 0)
!                else
                  if ((molofrs(rs1).eq.molofrs(rs2)).AND.(abs(rs2-rs1).eq.1)) call fexit()
                  if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
                  rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
                  if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
                  rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
                  rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
!                end if
              end if
            end if
          end if
        else
          if ((rs1.eq.disulf(rs2)).OR.(rs2.eq.disulf(rs1))) cycle ! exclude crosslinked pairs
          if (d1(ii).lt.(hlp1+resrad(i))) then
            if (rs2.eq.(rs1+1)) then
              rsp_vec(rs1) = 1
            else if (rs1.eq.(rs2+1)) then
              rsp_vec(rs2) = 1
            else
              rs_nbl(rs1)%nnbs = rs_nbl(rs1)%nnbs + 1
              if (rs_nbl(rs1)%nnbs.gt.rs_nbl(rs1)%nbalsz) call nbl_resz(rs1,1)
              rs_nbl(rs1)%nnbats = rs_nbl(rs1)%nnbats + at(rs2)%na
              rs_nbl(rs1)%nb(rs_nbl(rs1)%nnbs) = rs2
            end if
          else if (d1(ii).lt.(hlp2+resrad(i))) then
            if (rs2.eq.(rs1+1)) then
              rsp_vec(rs1) = 2
             else if (rs1.eq.(rs2+1)) then
              rsp_vec(rs2) = 2
            else
              if (rs1.eq.rs) then
                svsgn = 1.0
              else
                svsgn = -1.0
              end if
              rs_nbl(rs1)%nnbtrs = rs_nbl(rs1)%nnbtrs + 1
              if (rs_nbl(rs1)%nnbtrs.gt.rs_nbl(rs1)%tralsz) call nbl_resz(rs1,2)
              rs_nbl(rs1)%nnbtrats = rs_nbl(rs1)%nnbtrats + at(rs2)%na
              rs_nbl(rs1)%nbtr(rs_nbl(rs1)%nnbtrs) = rs2
              if (molofrs(rs2).eq.molofrs(rs1)) then
                rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,:) = 0.0
              else
                rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,1) = svsgn*svecx(ii)
                rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,2) = svsgn*svecy(ii)
                rs_nbl(rs1)%trsvec(rs_nbl(rs1)%nnbtrs,3) = svsgn*svecz(ii)
              end if
            end if
          else
            if (do_lr.EQV..true.) then
              if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!                if (rs2.eq.(rs1+1)) then
!                else if (rs1.eq.(rs2+1)) then
!                 do nothing (detected as rsp_vec(i) being 0)
!                else
                  if ((molofrs(rs1).eq.molofrs(rs2)).AND.(abs(rs2-rs1).eq.1)) call fexit()
                  if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
                  rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
                  if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
                  rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
                  rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
!                end if
              end if
            end if
          end if
        end if
      else
!       and finally, this is the case for interactions using the ghosted Hamiltonian
        if ((rs1.eq.disulf(rs2)).OR.(rs2.eq.disulf(rs1))) cycle ! exclude crosslinked pairs
        if (d1(ii).lt.(hlp1+resrad(i))) then
          if ((use_waterloops.EQV..false.).OR.((use_waterloops.EQV..true.).AND.(rs2.lt.rsw1))) then
            if (rs2.eq.(rs1+1)) then
              rsp_vec(rs1) = 1
              cycle
            else if (rs1.eq.(rs2+1)) then
              rsp_vec(rs2) = 1
              cycle
            end if
          end if
          rs_nbl(rs1)%ngnbs = rs_nbl(rs1)%ngnbs + 1
          if (rs_nbl(rs1)%ngnbs.gt.rs_nbl(rs1)%gnbalsz) call nbl_resz(rs1,4)
          rs_nbl(rs1)%ngnbats = rs_nbl(rs1)%ngnbats + at(rs2)%na
          rs_nbl(rs1)%gnb(rs_nbl(rs1)%ngnbs) = rs2
        else if (d1(ii).lt.(hlp2+resrad(i))) then
          if ((use_waterloops.EQV..false.).OR.((use_waterloops.EQV..true.).AND.(rs2.lt.rsw1))) then
            if (rs2.eq.(rs1+1)) then
              rsp_vec(rs1) = 2
              cycle
            else if (rs1.eq.(rs2+1)) then
              rsp_vec(rs2) = 2
              cycle
            end if
          end if
          if (rs1.eq.rs) then
            svsgn = 1.0
          else
            svsgn = -1.0
          end if
          rs_nbl(rs1)%ngnbtrs = rs_nbl(rs1)%ngnbtrs + 1
          if (rs_nbl(rs1)%ngnbtrs.gt.rs_nbl(rs1)%gtralsz) call nbl_resz(rs1,5)
          rs_nbl(rs1)%ngnbtrats = rs_nbl(rs1)%ngnbtrats + at(rs2)%na
          rs_nbl(rs1)%gnbtr(rs_nbl(rs1)%ngnbtrs) = rs2
          if (molofrs(rs2).eq.molofrs(rs1)) then
            rs_nbl(rs1)%gtrsvec(rs_nbl(rs1)%ngnbtrs,:) = 0.0
          else
            rs_nbl(rs1)%gtrsvec(rs_nbl(rs1)%ngnbtrs,1) = svsgn*svecx(ii)
            rs_nbl(rs1)%gtrsvec(rs_nbl(rs1)%ngnbtrs,2) = svsgn*svecy(ii)
            rs_nbl(rs1)%gtrsvec(rs_nbl(rs1)%ngnbtrs,3) = svsgn*svecz(ii)
          end if
        else
          if (do_lr.EQV..true.) then
            if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!              if (rs2.eq.(rs1+1)) then
!              else if (rs1.eq.(rs2+1)) then
!               do nothing (detected as rsp_vec(i) being 0)
!              else
                if ((molofrs(rs1).eq.molofrs(rs2)).AND.(abs(rs2-rs1).eq.1)) call fexit()
                if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
                rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
                if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
                rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
                rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
!              end if
            end if
          end if
        end if
      end if
    end do
  end if
!
end
!
!-----------------------------------------------------------------------------
!
!! the same at the charge group level
!
!subroutine populate_cgnbl(rs,cgi,hi,svecx,svecy,svecz,d1,idx,tpi)
!
!  use iounit
!  use molecule
!  use sequen
!  use energies
!  use cutoffs
!  use polypep
!  use system
!  use atoms
!  use fyoc
!
!  implicit none
!
!  integer, INTENT(IN):: hi,cgi,idx(hi),tpi
!  RTYPE, INTENT(IN):: svecx(hi),svecy(hi),svecz(hi),d1(hi)
!
!  integer i,ii,rs1,rs2,rs
!  RTYPE hlp1,hlp2,svsgn
!  logical do_lr
!
!  hlp1 = mcnb_cutoff+cgrad(rs)
!  hlp2 = mcel_cutoff+cgrad(rs)
!
!  if ((use_POLAR.EQV..true.).AND.((lrel_md.eq.5).OR.(lrel_md.eq.4))) then
!    do_lr = .true.
!  else
!    do_lr = .false.
!  end if
!  svsgn = 1.0
!
!  if (use_FEG.EQV..false.) then
!    do ii=1,hi
!      i = idx(ii)
!      if (tpi.eq.0) then
!        rs1 = min(cgi,i)
!        rs2 = max(cgi,i)
!      else
!        rs1 = cgi ! min(rs,i)
!        rs2 = i ! max(rs,i)
!      end if
!      if ((use_waterloops.EQV..true.).AND.(rs2.ge.rsw1)) then
!        if (d1(ii).lt.(hlp1+resrad(i))) then
!          cg_nbl(rs1)%nwnbs = cg_nbl(rs1)%nwnbs + 1
!          if (cg_nbl(rs1)%nwnbs.gt.cg_nbl(rs1)%wnbalsz) call nbl_resz(rs1,8)
!          cg_nbl(rs1)%nwnbats = cg_nbl(rs1)%nwnbats + at(rs2)%na
!          cg_nbl(rs1)%wnb(cg_nbl(rs1)%nwnbs) = rs2
!        else if (d1(ii).lt.(hlp2+resrad(i))) then
!          if (rs1.eq.rs) then
!            svsgn = 1.0
!          else
!            svsgn = -1.0
!          end if
!          cg_nbl(rs1)%nwnbtrs = cg_nbl(rs1)%nwnbtrs + 1
!          if (cg_nbl(rs1)%nwnbtrs.gt.cg_nbl(rs1)%wtralsz) call nbl_resz(rs1,9)
!          cg_nbl(rs1)%nwnbtrats = cg_nbl(rs1)%nwnbtrats + at(rs2)%na
!          cg_nbl(rs1)%wnbtr(cg_nbl(rs1)%nwnbtrs) = rs2
!          if (molofrs(rs2).eq.molofrs(rs1)) then
!            cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,:) = 0.0
!          else
!            cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,1) = svsgn*svecx(ii)
!            cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,2) = svsgn*svecy(ii)
!            cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,3) = svsgn*svecz(ii)
!          end if
!        else
!          if (do_lr.EQV..true.) then
!            if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!              if (rs2.eq.(rs1+1)) then
!              else if (rs1.eq.(rs2+1)) then
!!               do nothing (detected as rsp_vec(i) being 0)
!              else
!                if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
!                cg_nbl(rs1)%nnblrs = cg_nbl(rs1)%nnblrs + 1
!                if (cg_nbl(rs1)%nnblrs.gt.cg_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
!                cg_nbl(rs1)%nnblrats = cg_nbl(rs1)%nnblrats + at(rs2)%na
!                cg_nbl(rs1)%nblr(cg_nbl(rs1)%nnblrs) = rs2
!              end if
!            end if
!          end if
!        end if
!      else
!        if ((rs1.eq.disulf(rs2)).OR.(rs2.eq.disulf(rs1))) cycle ! exclude crosslinked pairs
!        if (d1(ii).lt.(hlp1+resrad(i))) then
!          if (rs2.eq.(rs1+1)) then
!            rsp_vec(rs1) = 1
!          else if (rs1.eq.(rs2+1)) then
!            rsp_vec(rs2) = 1
!          else
!            cg_nbl(rs1)%nnbs = cg_nbl(rs1)%nnbs + 1
!            if (cg_nbl(rs1)%nnbs.gt.cg_nbl(rs1)%nbalsz) call nbl_resz(rs1,1)
!            cg_nbl(rs1)%nnbats = cg_nbl(rs1)%nnbats + at(rs2)%na
!            cg_nbl(rs1)%nb(cg_nbl(rs1)%nnbs) = rs2
!          end if
!        else if (d1(ii).lt.(hlp2+resrad(i))) then
!          if (rs2.eq.(rs1+1)) then
!            rsp_vec(rs1) = 2
!           else if (rs1.eq.(rs2+1)) then
!            rsp_vec(rs2) = 2
!          else
!            if (rs1.eq.rs) then
!              svsgn = 1.0
!            else
!              svsgn = -1.0
!            end if
!            cg_nbl(rs1)%nnbtrs = cg_nbl(rs1)%nnbtrs + 1
!            if (cg_nbl(rs1)%nnbtrs.gt.cg_nbl(rs1)%tralsz) call nbl_resz(rs1,2)
!            cg_nbl(rs1)%nnbtrats = cg_nbl(rs1)%nnbtrats + at(rs2)%na
!            cg_nbl(rs1)%nbtr(cg_nbl(rs1)%nnbtrs) = rs2
!            if (molofrs(rs2).eq.molofrs(rs1)) then
!              cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,:) = 0.0
!            else
!              cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,1) = svsgn*svecx(ii)
!              cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,2) = svsgn*svecy(ii)
!              cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,3) = svsgn*svecz(ii)
!            end if
!          end if
!        else
!          if (do_lr.EQV..true.) then
!            if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!              if (rs2.eq.(rs1+1)) then
!              else if (rs1.eq.(rs2+1)) then
!!               do nothing (detected as rsp_vec(i) being 0)
!              else
!                if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
!                cg_nbl(rs1)%nnblrs = cg_nbl(rs1)%nnblrs + 1
!                if (cg_nbl(rs1)%nnblrs.gt.cg_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
!                cg_nbl(rs1)%nnblrats = cg_nbl(rs1)%nnblrats + at(rs2)%na
!                cg_nbl(rs1)%nblr(cg_nbl(rs1)%nnblrs) = rs2
!              end if
!            end if
!          end if
!        end if
!      end if
!    end do
!
!  else
!
!!   for calculations with ghosts we have (annoyingly) separate nb-lists
!    do ii=1,hi
!      i = idx(ii)
!      if (tpi.eq.0) then
!        rs1 = min(rs,i)
!        rs2 = max(rs,i)
!      else
!        rs1 = rs ! min(rs,i)
!        rs2 = i ! max(rs,i)
!      end if
!      if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!!       do nothing (at least one of the residues is fully de-coupled)
!      else if (((par_FEG(rs1).EQV..false.).AND.(par_FEG(rs2).EQV..false.)).OR.&
! &    ((par_FEG(rs1).EQV..true.).AND.(par_FEG(rs2).EQV..true.).AND.(fegmode.eq.1))) then
!!       this is the case using the full Hamiltonian, same code as above
!        if ((use_waterloops.EQV..true.).AND.(rs2.ge.rsw1)) then
!          if (d1(ii).lt.(hlp1+resrad(i))) then
!            cg_nbl(rs1)%nwnbs = cg_nbl(rs1)%nwnbs + 1
!            if (cg_nbl(rs1)%nwnbs.gt.cg_nbl(rs1)%wnbalsz) call nbl_resz(rs1,8)
!            cg_nbl(rs1)%nwnbats = cg_nbl(rs1)%nwnbats + at(rs2)%na
!            cg_nbl(rs1)%wnb(cg_nbl(rs1)%nwnbs) = rs2
!          else if (d1(ii).lt.(hlp2+resrad(i))) then
!            if (rs1.eq.rs) then
!              svsgn = 1.0
!            else
!              svsgn = -1.0
!            end if
!            cg_nbl(rs1)%nwnbtrs = cg_nbl(rs1)%nwnbtrs + 1
!            if (cg_nbl(rs1)%nwnbtrs.gt.cg_nbl(rs1)%wtralsz) call nbl_resz(rs1,9)
!            cg_nbl(rs1)%nwnbtrats = cg_nbl(rs1)%nwnbtrats + at(rs2)%na
!            cg_nbl(rs1)%wnbtr(cg_nbl(rs1)%nwnbtrs) = rs2
!            if (molofrs(rs2).eq.molofrs(rs1)) then
!              cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,:) = 0.0
!            else
!              cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,1) = svsgn*svecx(ii)
!              cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,2) = svsgn*svecy(ii)
!              cg_nbl(rs1)%wtrsvec(cg_nbl(rs1)%nwnbtrs,3) = svsgn*svecz(ii)
!            end if
!          else
!            if (do_lr.EQV..true.) then
!              if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!                if (rs2.eq.(rs1+1)) then
!                else if (rs1.eq.(rs2+1)) then
!!                 do nothing (detected as rsp_vec(i) being 0)
!                else
!                  if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
!                  cg_nbl(rs1)%nnblrs = cg_nbl(rs1)%nnblrs + 1
!                  if (cg_nbl(rs1)%nnblrs.gt.cg_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
!                  cg_nbl(rs1)%nnblrats = cg_nbl(rs1)%nnblrats + at(rs2)%na
!                  cg_nbl(rs1)%nblr(cg_nbl(rs1)%nnblrs) = rs2
!                end if
!              end if
!            end if
!          end if
!        else
!          if ((rs1.eq.disulf(rs2)).OR.(rs2.eq.disulf(rs1))) cycle ! exclude crosslinked pairs
!          if (d1(ii).lt.(hlp1+resrad(i))) then
!            if (rs2.eq.(rs1+1)) then
!              rsp_vec(rs1) = 1
!            else if (rs1.eq.(rs2+1)) then
!              rsp_vec(rs2) = 1
!            else
!              cg_nbl(rs1)%nnbs = cg_nbl(rs1)%nnbs + 1
!              if (cg_nbl(rs1)%nnbs.gt.cg_nbl(rs1)%nbalsz) call nbl_resz(rs1,1)
!              cg_nbl(rs1)%nnbats = cg_nbl(rs1)%nnbats + at(rs2)%na
!              cg_nbl(rs1)%nb(cg_nbl(rs1)%nnbs) = rs2
!            end if
!          else if (d1(ii).lt.(hlp2+resrad(i))) then
!            if (rs2.eq.(rs1+1)) then
!              rsp_vec(rs1) = 2
!             else if (rs1.eq.(rs2+1)) then
!              rsp_vec(rs2) = 2
!            else
!              if (rs1.eq.rs) then
!                svsgn = 1.0
!              else
!                svsgn = -1.0
!              end if
!              cg_nbl(rs1)%nnbtrs = cg_nbl(rs1)%nnbtrs + 1
!              if (cg_nbl(rs1)%nnbtrs.gt.cg_nbl(rs1)%tralsz) call nbl_resz(rs1,2)
!              cg_nbl(rs1)%nnbtrats = cg_nbl(rs1)%nnbtrats + at(rs2)%na
!              cg_nbl(rs1)%nbtr(cg_nbl(rs1)%nnbtrs) = rs2
!              if (molofrs(rs2).eq.molofrs(rs1)) then
!                cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,:) = 0.0
!              else
!                cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,1) = svsgn*svecx(ii)
!                cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,2) = svsgn*svecy(ii)
!                cg_nbl(rs1)%trsvec(cg_nbl(rs1)%nnbtrs,3) = svsgn*svecz(ii)
!              end if
!            end if
!          else
!            if (do_lr.EQV..true.) then
!              if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!                if (rs2.eq.(rs1+1)) then
!                else if (rs1.eq.(rs2+1)) then
!!                 do nothing (detected as rsp_vec(i) being 0)
!                else
!                  if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
!                  cg_nbl(rs1)%nnblrs = cg_nbl(rs1)%nnblrs + 1
!                  if (cg_nbl(rs1)%nnblrs.gt.cg_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
!                  cg_nbl(rs1)%nnblrats = cg_nbl(rs1)%nnblrats + at(rs2)%na
!                  cg_nbl(rs1)%nblr(cg_nbl(rs1)%nnblrs) = rs2
!                end if
!              end if
!            end if
!          end if
!        end if
!      else
!!       and finally, this is the case for interactions using the ghosted Hamiltonian
!        if ((rs1.eq.disulf(rs2)).OR.(rs2.eq.disulf(rs1))) cycle ! exclude crosslinked pairs
!        if (d1(ii).lt.(hlp1+resrad(i))) then
!          if ((use_waterloops.EQV..false.).OR.((use_waterloops.EQV..true.).AND.(rs2.lt.rsw1))) then
!            if (rs2.eq.(rs1+1)) then
!              rsp_vec(rs1) = 1
!              cycle
!            else if (rs1.eq.(rs2+1)) then
!              rsp_vec(rs2) = 1
!              cycle
!            end if
!          end if
!          cg_nbl(rs1)%ngnbs = cg_nbl(rs1)%ngnbs + 1
!          if (cg_nbl(rs1)%ngnbs.gt.cg_nbl(rs1)%gnbalsz) call nbl_resz(rs1,4)
!          cg_nbl(rs1)%ngnbats = cg_nbl(rs1)%ngnbats + at(rs2)%na
!          cg_nbl(rs1)%gnb(cg_nbl(rs1)%ngnbs) = rs2
!        else if (d1(ii).lt.(hlp2+resrad(i))) then
!          if ((use_waterloops.EQV..false.).OR.((use_waterloops.EQV..true.).AND.(rs2.lt.rsw1))) then
!            if (rs2.eq.(rs1+1)) then
!              rsp_vec(rs1) = 2
!              cycle
!            else if (rs1.eq.(rs2+1)) then
!              rsp_vec(rs2) = 2
!              cycle
!            end if
!          end if
!          if (rs1.eq.rs) then
!            svsgn = 1.0
!          else
!            svsgn = -1.0
!          end if
!          cg_nbl(rs1)%ngnbtrs = cg_nbl(rs1)%ngnbtrs + 1
!          if (cg_nbl(rs1)%ngnbtrs.gt.cg_nbl(rs1)%gtralsz) call nbl_resz(rs1,5)
!          cg_nbl(rs1)%ngnbtrats = cg_nbl(rs1)%ngnbtrats + at(rs2)%na
!          cg_nbl(rs1)%gnbtr(cg_nbl(rs1)%ngnbtrs) = rs2
!          if (molofrs(rs2).eq.molofrs(rs1)) then
!            cg_nbl(rs1)%gtrsvec(cg_nbl(rs1)%ngnbtrs,:) = 0.0
!          else
!            cg_nbl(rs1)%gtrsvec(cg_nbl(rs1)%ngnbtrs,1) = svsgn*svecx(ii)
!            cg_nbl(rs1)%gtrsvec(cg_nbl(rs1)%ngnbtrs,2) = svsgn*svecy(ii)
!            cg_nbl(rs1)%gtrsvec(cg_nbl(rs1)%ngnbtrs,3) = svsgn*svecz(ii)
!          end if
!        else
!          if (do_lr.EQV..true.) then
!            if ((chgflag(rs1).EQV..true.).OR.(chgflag(rs2).EQV..true.)) then
!              if (rs2.eq.(rs1+1)) then
!              else if (rs1.eq.(rs2+1)) then
!!               do nothing (detected as rsp_vec(i) being 0)
!              else
!                if ((lrel_md.eq.4).AND.((chgflag(rs1).EQV..false.).OR.(chgflag(rs2).EQV..false.))) cycle
!                cg_nbl(rs1)%nnblrs = cg_nbl(rs1)%nnblrs + 1
!                if (cg_nbl(rs1)%nnblrs.gt.cg_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
!                cg_nbl(rs1)%nnblrats = cg_nbl(rs1)%nnblrats + at(rs2)%na
!                cg_nbl(rs1)%nblr(cg_nbl(rs1)%nnblrs) = rs2
!              end if
!            end if
!          end if
!        end if
!      end if
!    end do
!  end if
!!
!end
!
!-----------------------------------------------------------------------------
!
! a subroutine to resize a neighbor list by a fixed increment
!
subroutine nbl_resz(rs,which)
!
  use iounit
  use cutoffs
  use sequen
!
  implicit none
!
  integer rs,which
!
#ifdef ENABLE_THREADS
!$OMP CRITICAL(NBRSZ7)
#endif
! normal close-range
  if (which.eq.1) then
!   backup, deallocate, increment, reallocate, copy back
    rsinfo(1:rs_nbl(rs)%nbalsz,3) = rs_nbl(rs)%nb(1:rs_nbl(rs)%nbalsz)
    deallocate(rs_nbl(rs)%nb)
    rs_nbl(rs)%nbalsz = min(nseq,rs_nbl(rs)%nbalsz + 10)
    allocate(rs_nbl(rs)%nb(rs_nbl(rs)%nbalsz))
    rs_nbl(rs)%nb(1:rs_nbl(rs)%nbalsz) = rsinfo(1:rs_nbl(rs)%nbalsz,3)
! normal twin-range
  else if (which.eq.2) then
    rsinfo(1:rs_nbl(rs)%tralsz,3) = rs_nbl(rs)%nbtr(1:rs_nbl(rs)%tralsz)
    vhlper(1:rs_nbl(rs)%tralsz,1) = rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,1)
    vhlper(1:rs_nbl(rs)%tralsz,2) = rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,2)
    vhlper(1:rs_nbl(rs)%tralsz,3) = rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,3)
    deallocate(rs_nbl(rs)%nbtr)
    deallocate(rs_nbl(rs)%trsvec)
    rs_nbl(rs)%tralsz = min(nseq,rs_nbl(rs)%tralsz + 10)
    allocate(rs_nbl(rs)%nbtr(rs_nbl(rs)%tralsz))
    allocate(rs_nbl(rs)%trsvec(rs_nbl(rs)%tralsz,3))
    rs_nbl(rs)%nbtr(1:rs_nbl(rs)%tralsz) = rsinfo(1:rs_nbl(rs)%tralsz,3)
    rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,1) = vhlper(1:rs_nbl(rs)%tralsz,1)
    rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,2) = vhlper(1:rs_nbl(rs)%tralsz,2)
    rs_nbl(rs)%trsvec(1:rs_nbl(rs)%tralsz,3) = vhlper(1:rs_nbl(rs)%tralsz,3)
! normal long-range
  else if (which.eq.3) then
    rsinfo(1:rs_nbl(rs)%lralsz,3) = rs_nbl(rs)%nblr(1:rs_nbl(rs)%lralsz)
    deallocate(rs_nbl(rs)%nblr)
    rs_nbl(rs)%lralsz = min(nseq,rs_nbl(rs)%lralsz + 10)
    allocate(rs_nbl(rs)%nblr(rs_nbl(rs)%lralsz))
    rs_nbl(rs)%nblr(1:rs_nbl(rs)%lralsz) = rsinfo(1:rs_nbl(rs)%lralsz,3)
! ghosted close-range
  else if (which.eq.4) then
    rsinfo(1:rs_nbl(rs)%gnbalsz,3) = rs_nbl(rs)%gnb(1:rs_nbl(rs)%gnbalsz)
    deallocate(rs_nbl(rs)%gnb)
    rs_nbl(rs)%gnbalsz = min(nseq,rs_nbl(rs)%gnbalsz + 10)
    allocate(rs_nbl(rs)%gnb(rs_nbl(rs)%gnbalsz))
    rs_nbl(rs)%gnb(1:rs_nbl(rs)%gnbalsz) = rsinfo(1:rs_nbl(rs)%gnbalsz,3)
! ghosted twin-range
  else if (which.eq.5) then
    rsinfo(1:rs_nbl(rs)%gtralsz,3) = rs_nbl(rs)%gnbtr(1:rs_nbl(rs)%gtralsz)
    vhlper(1:rs_nbl(rs)%gtralsz,1) = rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,1)
    vhlper(1:rs_nbl(rs)%gtralsz,2) = rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,2)
    vhlper(1:rs_nbl(rs)%gtralsz,3) = rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,3)
    deallocate(rs_nbl(rs)%gnbtr)
    deallocate(rs_nbl(rs)%gtrsvec)
    rs_nbl(rs)%gtralsz = min(nseq,rs_nbl(rs)%gtralsz + 10)
    allocate(rs_nbl(rs)%gnbtr(rs_nbl(rs)%gtralsz))
    allocate(rs_nbl(rs)%gtrsvec(rs_nbl(rs)%gtralsz,3))
    rs_nbl(rs)%gnbtr(1:rs_nbl(rs)%gtralsz) = rsinfo(1:rs_nbl(rs)%gtralsz,3)
    rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,1) = vhlper(1:rs_nbl(rs)%gtralsz,1)
    rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,2) = vhlper(1:rs_nbl(rs)%gtralsz,2)
    rs_nbl(rs)%gtrsvec(1:rs_nbl(rs)%gtralsz,3) = vhlper(1:rs_nbl(rs)%gtralsz,3)
! temporary index list
  else if (which.eq.7) then
    rsinfo(1:rs_nbl(rs)%tmpalsz,3) = rs_nbl(rs)%tmpanb(1:rs_nbl(rs)%tmpalsz)
    deallocate(rs_nbl(rs)%tmpanb)
    rs_nbl(rs)%tmpalsz = min(nseq,max(rs_nbl(rs)%tmpalsz + 100,rs_nbl(rs)%ntmpanb))
    allocate(rs_nbl(rs)%tmpanb(rs_nbl(rs)%tmpalsz))
    rs_nbl(rs)%tmpanb(1:rs_nbl(rs)%tmpalsz) = rsinfo(1:rs_nbl(rs)%tmpalsz,3)
! water loops close-range
  else if (which.eq.8) then
    rsinfo(1:rs_nbl(rs)%wnbalsz,3) = rs_nbl(rs)%wnb(1:rs_nbl(rs)%wnbalsz)
    if (allocated(rs_nbl(rs)%wnb).EQV..true.) deallocate(rs_nbl(rs)%wnb)
    rs_nbl(rs)%wnbalsz = min(nseq,rs_nbl(rs)%wnbalsz + 10)
    allocate(rs_nbl(rs)%wnb(rs_nbl(rs)%wnbalsz))
    rs_nbl(rs)%wnb(1:rs_nbl(rs)%wnbalsz) = rsinfo(1:rs_nbl(rs)%wnbalsz,3)
! water loops twin-range
  else if (which.eq.9) then
    rsinfo(1:rs_nbl(rs)%wtralsz,3) = rs_nbl(rs)%wnbtr(1:rs_nbl(rs)%wtralsz)
    vhlper(1:rs_nbl(rs)%wtralsz,1) = rs_nbl(rs)%wtrsvec(1:rs_nbl(rs)%wtralsz,1)
    vhlper(1:rs_nbl(rs)%wtralsz,2) = rs_nbl(rs)%wtrsvec(1:rs_nbl(rs)%wtralsz,2)
    vhlper(1:rs_nbl(rs)%wtralsz,3) = rs_nbl(rs)%wtrsvec(1:rs_nbl(rs)%wtralsz,3)
    deallocate(rs_nbl(rs)%wnbtr)
    deallocate(rs_nbl(rs)%wtrsvec)
    rs_nbl(rs)%wtralsz = min(nseq,rs_nbl(rs)%wtralsz + 10)
    allocate(rs_nbl(rs)%wnbtr(rs_nbl(rs)%wtralsz))
    allocate(rs_nbl(rs)%wtrsvec(rs_nbl(rs)%wtralsz,3))
    rs_nbl(rs)%wnbtr(1:rs_nbl(rs)%wtralsz) = rsinfo(1:rs_nbl(rs)%wtralsz,3)
    rs_nbl(rs)%wtrsvec(1:rs_nbl(rs)%wtralsz,1) = vhlper(1:rs_nbl(rs)%wtralsz,1)
    rs_nbl(rs)%wtrsvec(1:rs_nbl(rs)%wtralsz,2) = vhlper(1:rs_nbl(rs)%wtralsz,2)
    rs_nbl(rs)%wtrsvec(1:rs_nbl(rs)%wtralsz,3) = vhlper(1:rs_nbl(rs)%wtralsz,3)
! water loops long-range
  else if (which.eq.10) then
    rsinfo(1:rs_nbl(rs)%wlralsz,3) = rs_nbl(rs)%wnblr(1:rs_nbl(rs)%wlralsz)
    deallocate(rs_nbl(rs)%wnblr)
    rs_nbl(rs)%wlralsz = min(nseq,rs_nbl(rs)%wlralsz + 10)
    allocate(rs_nbl(rs)%wnblr(rs_nbl(rs)%wlralsz))
    rs_nbl(rs)%wnblr(1:rs_nbl(rs)%wlralsz) = rsinfo(1:rs_nbl(rs)%wlralsz,3)
  else
    write(ilog,*) 'Fatal. Called nbl_resz(...) with unkown list identifier (offending &
 &code is ',which,'). Please report this bug.'
    call fexit()
  end if
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(NBRSZ7)
#endif
! 
end
!
!-----------------------------------------------------------------------------
!
! a subroutine to resize a pair neighbor list by a factor of 2
! a fixed increment would fail in cases where the list gets very large, which is easy
! to realize by certain choices of cutoffs, LREL_MC, system sizes, ...
!
subroutine rsp_nbl_resz(which,ix)
!
  use iounit
  use cutoffs
  use sequen
!
  implicit none
!
  integer, INTENT(IN):: which,ix
!
  integer, ALLOCATABLE:: hlper(:,:)
!
#ifdef ENABLE_THREADS
!$OMP CRITICAL(NBRSZ0)
#endif
! normal close-range
  if (which.eq.1) then
!   backup, deallocate, increment, reallocate, copy back
    allocate(hlper(rsp_nbl(ix)%sralsz,2))
    hlper(1:rsp_nbl(ix)%sralsz,1:2) = rsp_nbl(ix)%sr(1:rsp_nbl(ix)%sralsz,1:2)
    deallocate(rsp_nbl(ix)%sr)
    rsp_nbl(ix)%sralsz = rsp_nbl(ix)%sralsz*2
    allocate(rsp_nbl(ix)%sr(rsp_nbl(ix)%sralsz,2))
    rsp_nbl(ix)%sr(1:(rsp_nbl(ix)%sralsz/2),1:2) = hlper(1:(rsp_nbl(ix)%sralsz/2),1:2)
! normal twin-range
  else if (which.eq.2) then
    allocate(hlper(rsp_nbl(ix)%tralsz,2))
    hlper(1:rsp_nbl(ix)%tralsz,1:2) = rsp_nbl(ix)%tr(1:rsp_nbl(ix)%tralsz,1:2)
    deallocate(rsp_nbl(ix)%tr)
    rsp_nbl(ix)%tralsz = rsp_nbl(ix)%tralsz*2
    allocate(rsp_nbl(ix)%tr(rsp_nbl(ix)%tralsz,2))
    rsp_nbl(ix)%tr(1:(rsp_nbl(ix)%tralsz/2),1:2) = hlper(1:(rsp_nbl(ix)%tralsz/2),1:2)
! normal long-range
  else if (which.eq.3) then
    allocate(hlper(rsp_nbl(ix)%lralsz,2))
    hlper(1:rsp_nbl(ix)%lralsz,1:2) = rsp_nbl(ix)%lr(1:rsp_nbl(ix)%lralsz,1:2)
    deallocate(rsp_nbl(ix)%lr)
    rsp_nbl(ix)%lralsz = rsp_nbl(ix)%lralsz*2
    allocate(rsp_nbl(ix)%lr(rsp_nbl(ix)%lralsz,2))
    rsp_nbl(ix)%lr(1:(rsp_nbl(ix)%lralsz/2),1:2) = hlper(1:(rsp_nbl(ix)%lralsz/2),1:2)
  else
    write(ilog,*) 'Fatal. Called rsp_nbl_resz(...) with unkown list identifier (offending &
 &code is ',which,'). Please report this bug.'
    call fexit()
  end if
!
  if (allocated(hlper).EQV..true.) deallocate(hlper)
!
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(NBRSZ0)
#endif
! 
end
!
!-----------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
! a subroutine to resize a threads-only pair neighbor list by a factor of 2
! a fixed increment would fail in cases where the list gets very large, which is easy
! to realize by certain choices of cutoffs, LREL_MC, system sizes, ...
!
subroutine thr_rsp_nbl_resz(tpi,which,ix)
!
  use iounit
  use cutoffs
  use sequen
  use threads
!
  implicit none
!
  integer, INTENT(IN):: which,tpi,ix
!
  integer tpx
  integer, ALLOCATABLE:: hlper(:,:)
!
  if (ix.eq.2) then
    tpx = thrdat%maxn + tpi
  else 
    tpx = tpi
  end if
#ifdef ENABLE_THREADS
!$OMP CRITICAL(NBRSZ3)
#endif
! normal close-range
  if (which.eq.1) then
!   backup, deallocate, increment, reallocate, copy back
    allocate(hlper(thr_rsp_nbl(tpx)%sralsz,3))
    hlper(1:thr_rsp_nbl(tpx)%sralsz,1:3) = thr_rsp_nbl(tpx)%sr(1:thr_rsp_nbl(tpx)%sralsz,1:3)
    deallocate(thr_rsp_nbl(tpx)%sr)
    thr_rsp_nbl(tpx)%sralsz = thr_rsp_nbl(tpx)%sralsz*2
    allocate(thr_rsp_nbl(tpx)%sr(thr_rsp_nbl(tpx)%sralsz,3))
    thr_rsp_nbl(tpx)%sr(1:(thr_rsp_nbl(tpx)%sralsz/2),1:3) = hlper(1:(thr_rsp_nbl(tpx)%sralsz/2),1:3)
! normal twin-range
  else if (which.eq.2) then
    allocate(hlper(thr_rsp_nbl(tpx)%tralsz,3))
    hlper(1:thr_rsp_nbl(tpx)%tralsz,1:3) = thr_rsp_nbl(tpx)%tr(1:thr_rsp_nbl(tpx)%tralsz,1:3)
    deallocate(thr_rsp_nbl(tpx)%tr)
    thr_rsp_nbl(tpx)%tralsz = thr_rsp_nbl(tpx)%tralsz*2
    allocate(thr_rsp_nbl(tpx)%tr(thr_rsp_nbl(tpx)%tralsz,3))
    thr_rsp_nbl(tpx)%tr(1:(thr_rsp_nbl(tpx)%tralsz/2),1:3) = hlper(1:(thr_rsp_nbl(tpx)%tralsz/2),1:3)
! normal long-range
  else if (which.eq.3) then
    allocate(hlper(thr_rsp_nbl(tpx)%lralsz,3))
    hlper(1:thr_rsp_nbl(tpx)%lralsz,1:3) = thr_rsp_nbl(tpx)%lr(1:thr_rsp_nbl(tpx)%lralsz,1:3)
    deallocate(thr_rsp_nbl(tpx)%lr)
    thr_rsp_nbl(tpx)%lralsz = thr_rsp_nbl(tpx)%lralsz*2
    allocate(thr_rsp_nbl(tpx)%lr(thr_rsp_nbl(tpx)%lralsz,3))
    thr_rsp_nbl(tpx)%lr(1:(thr_rsp_nbl(tpx)%lralsz/2),1:3) = hlper(1:(thr_rsp_nbl(tpx)%lralsz/2),1:3)
  else
    write(ilog,*) 'Fatal. Called thr_rsp_nbl_resz(...) with unkown list identifier (offending &
 &code is ',which,'). Please report this bug.'
    call fexit()
  end if
!
  if (allocated(hlper).EQV..true.) deallocate(hlper)
!
#ifdef ENABLE_THREADS
!$OMP END CRITICAL(NBRSZ3)
#endif
! 
end
!
#endif
!
!-----------------------------------------------------------------------------
!
! with lrel_md set to 4 or 5, the entire local concept breaks down
! hence, we have to do a significant amount of clean-up work which may be prohibitively expensive
! this residue-based routine works only with rs_nbl-structures that have strict index
! compliance in that entries always exceed the residue number of the index into into rs_nbl itself
!
! this is why it cannot be used in threaded code where the index swapping is not supported (see populate_nbl)
!
subroutine complete_nbl(rs1,cycle_frz)
!
  use cutoffs
  use mcgrid
  use energies
  use sequen
  use polypep
  use fyoc
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: rs1
  logical, INTENT(IN):: cycle_frz
!
  integer imol,jmol,rs2,j,jj,proxna
!
  integer isthere(nseq)
!
  if (use_POLAR.EQV..false.) return
  if ((lrel_md.eq.4).AND.(chgflag(rs1).EQV..false.)) return
!
  isthere(:) = 0
! go through all neighbor list to tag ones we have in a list already
  do jj=1,rs_nbl(rs1)%nnbs
    isthere(rs_nbl(rs1)%nb(jj)) = 1
  end do
  do jj=1,rs_nbl(rs1)%nnbtrs
    isthere(rs_nbl(rs1)%nbtr(jj)) = 1
  end do
  do jj=1,rs_nbl(rs1)%nnblrs
    isthere(rs_nbl(rs1)%nblr(jj)) = 1
  end do
  do jj=1,rs_nbl(rs1)%nwnbs
    isthere(rs_nbl(rs1)%wnb(jj)) = 1
  end do
  do jj=1,rs_nbl(rs1)%nwnbtrs
    isthere(rs_nbl(rs1)%wnbtr(jj)) = 1
  end do
  do jj=1,rs_nbl(rs1)%ngnbs
    isthere(rs_nbl(rs1)%gnb(jj)) = 1
  end do
  do jj=1,rs_nbl(rs1)%ngnbtrs
    isthere(rs_nbl(rs1)%gnbtr(jj)) = 1
  end do
  if (rsp_vec(rs1).gt.0) then
    isthere(rs1+1) = 1
  end if
  imol = molofrs(rs1)
!
  if ((chgflag(rs1).EQV..false.).OR.(lrel_md.eq.4)) then ! just search other monopoles to complete list
    do j=1,cglst%ncrs
      rs2 = cglst%rsl(j)
      if (isthere(rs2).eq.1) cycle  ! exclude short- and mid-range interactions 
      if (chgflag(rs2).EQV..false.) cycle ! otherwise a patched net charge flag is not respected w/ grid-based cutoffs
      if (rs2.eq.disulf(rs1)) cycle ! exclude crosslinked pairs
      jmol = molofrs(rs2)
      if (rs2.le.rs1) cycle
      if (cycle_frz.EQV..true.) then! exclude frozen interactions 
        if (imol.eq.jmol) then
          if (molfrzidx(imol).ge.3) cycle
          if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
        else
          if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
        end if
      end if
      if (lrel_md.eq.4) then
        proxna = at(rs2)%nncgrps
      else
        proxna = at(rs2)%npol
      end if
      isthere(rs2) = 1 ! a residue may contain more than one monopole -> avoid double counting
      if (use_FEG.EQV..false.) then
        rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
        if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
        rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
        rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
      else
        if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!         do nothing (at least one of the residues is fully de-coupled)
        else
          rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
          if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
          rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + proxna
          rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
        end if
      end if
    end do 
  else ! search all other residues to complete list
    do rs2=rs1+1,nseq
      if (isthere(rs2).eq.1) cycle
      if (rs2.eq.disulf(rs1)) cycle ! exclude crosslinked pairs
      jmol = molofrs(rs2)
      if (cycle_frz.EQV..true.) then
        if (imol.eq.jmol) then
          if (molfrzidx(imol).ge.3) cycle
          if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
        else
          if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
        end if
      end if
      if (use_FEG.EQV..false.) then
        rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
        if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
        rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%npol
        rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
      else
        if ((par_FEG3(rs1).EQV..true.).OR.(par_FEG3(rs2).EQV..true.)) then
!         do nothing (at least one of the residues is fully de-coupled)
        else
          rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
          if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
          rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%npol
          rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
        end if
      end if
    end do 
  end if
!
end
!
!------------------------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS
!
subroutine complete_nbl_threads(cycle_frz,tpi)
!
  use cutoffs
  use threads
  use iounit
  use sequen, ONLY: nseq,chgflag,molofrs
  use polypep, ONLY: at
  use fyoc, ONLY: nchi
  use molecule, ONLY: molfrzidx
!
  implicit none
!
  integer, INTENT(IN):: tpi
  logical, INTENT(IN):: cycle_frz
!
  integer rs1,rs2,rsx,rsy,jj,imol,jmol
  integer(KIND=8):: ttimer
!
  if (cglst%ncrs.eq.0) return
!
  if (tpi.le.0) then
    write(ilog,*) 'Fatal. Subroutine complete_nbl_threads(...) must only be called from a parallel environment. This is a bug.'
    call fexit()
  end if
!
!$OMP SINGLE
  thr_lhlper(:,:) = .false.
!$OMP END SINGLE
  if (lrel_md.eq.5) then
!   interactions are unique in rs_nbl, therefore simultaneous write to thr_lhlper is valid
    do rs1=thr_limits(31,tpi),thr_limits(32,tpi)
      if (chgflag(rs1).EQV..true.) then
        rsx = cglst%irsl(rs1)
        if (rsp_vec(rs1).gt.0) thr_lhlper(rsx,rs1+1) = .true.
        do jj=1,rs_nbl(rs1)%nnbs
          thr_lhlper(rsx,rs_nbl(rs1)%nb(jj)) = .true.
        end do
        do jj=1,rs_nbl(rs1)%nnbtrs
          thr_lhlper(rsx,rs_nbl(rs1)%nbtr(jj)) = .true.
        end do
        do jj=1,rs_nbl(rs1)%nnblrs
          thr_lhlper(rsx,rs_nbl(rs1)%nblr(jj)) = .true.
        end do
        do jj=1,rs_nbl(rs1)%nwnbs
          thr_lhlper(rsx,rs_nbl(rs1)%wnb(jj)) = .true.
        end do
        do jj=1,rs_nbl(rs1)%nwnbtrs
          thr_lhlper(rsx,rs_nbl(rs1)%wnbtr(jj)) = .true.
        end do
        do jj=1,rs_nbl(rs1)%ngnbs
          thr_lhlper(rsx,rs_nbl(rs1)%gnb(jj)) = .true.
        end do
        do jj=1,rs_nbl(rs1)%ngnbtrs
          thr_lhlper(rsx,rs_nbl(rs1)%gnbtr(jj)) = .true.
        end do
      else
        if (rs1.lt.nseq) then
          if ((rsp_vec(rs1).gt.0).AND.(chgflag(rs1+1).EQV..true.)) thr_lhlper(cglst%irsl(rs1+1),rs1) = .true.
        end if
        do jj=1,rs_nbl(rs1)%nnbs
          if (chgflag(rs_nbl(rs1)%nb(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%nb(jj)),rs1) = .true.
          end if
        end do
        do jj=1,rs_nbl(rs1)%nnbtrs
          if (chgflag(rs_nbl(rs1)%nbtr(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%nbtr(jj)),rs1) = .true.
          end if
        end do
        do jj=1,rs_nbl(rs1)%nnblrs
          if (chgflag(rs_nbl(rs1)%nblr(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%nblr(jj)),rs1) = .true.
          end if
        end do
        do jj=1,rs_nbl(rs1)%ngnbs
          if (chgflag(rs_nbl(rs1)%gnb(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%gnb(jj)),rs1) = .true.
          end if
        end do
        do jj=1,rs_nbl(rs1)%ngnbtrs
          if (chgflag(rs_nbl(rs1)%gnbtr(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%gnbtr(jj)),rs1) = .true.
          end if
        end do
        do jj=1,rs_nbl(rs1)%nwnbs
          if (chgflag(rs_nbl(rs1)%wnb(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%wnb(jj)),rs1) = .true.
          end if
        end do
        do jj=1,rs_nbl(rs1)%nwnbtrs
          if (chgflag(rs_nbl(rs1)%wnbtr(jj)).EQV..true.) then
            thr_lhlper(cglst%irsl(rs_nbl(rs1)%wnbtr(jj)),rs1) = .true.
          end if
        end do
      end if
    end do
!$OMP BARRIER
    if (thr_dlb(5,1).gt.0) then
      if (tpi.eq.1) thr_dlb(5,2) = thr_dlb(5,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
    end if
    do rsx=thr_limits(27,tpi),thr_limits(28,tpi) ! 1,cglst%ncrs
      rs1 = cglst%rsl(rsx)
      imol = molofrs(rs1)
      do rs2=1,nseq
        if (rs1.eq.rs2) cycle
        if (cycle_frz.EQV..true.) then
          jmol = molofrs(rs2)
          if (imol.eq.jmol) then
            if (molfrzidx(imol).ge.3) cycle
            if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
          else
            if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
          end if
        end if
        if (chgflag(rs2).EQV..true.) then
          if (cglst%irsl(rs2).lt.rsx) cycle
          if ((thr_lhlper(rsx,rs2).EQV..false.).AND.(thr_lhlper(cglst%irsl(rs2),rs1).EQV..false.)) then
            rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
            if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
            rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%npol
            rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
          end if
        else
          if (thr_lhlper(rsx,rs2).EQV..false.) then
            rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
            if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
            rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%npol
            rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
          end if
        end if
      end do
    end do
    if (thr_dlb(5,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
    end if
!
  else if (lrel_md.eq.4) then
!   interactions are unique in rs_nbl, therefore simultaneous write to thr_lhlper is valid
    do rsx=tpi,cglst%ncrs,thrdat%maxn ! 1,cglst%ncrs
      rs1 = cglst%rsl(rsx)
      if ((rsp_vec(rs1).gt.0).AND.(chgflag(rs1+1).EQV..true.)) thr_lhlper(rsx,cglst%irsl(rs1+1)) = .true.
      do jj=1,rs_nbl(rs1)%nnbs
        if (chgflag(rs_nbl(rs1)%nb(jj)).EQV..true.) then
          thr_lhlper(rsx,cglst%irsl(rs_nbl(rs1)%nb(jj))) = .true.
        end if
      end do
      do jj=1,rs_nbl(rs1)%nnbtrs
        if (chgflag(rs_nbl(rs1)%nbtr(jj)).EQV..true.) then
          thr_lhlper(rsx,cglst%irsl(rs_nbl(rs1)%nbtr(jj))) = .true.
        end if
      end do
      do jj=1,rs_nbl(rs1)%nnblrs
        if (chgflag(rs_nbl(rs1)%nblr(jj)).EQV..true.) then
          thr_lhlper(rsx,cglst%irsl(rs_nbl(rs1)%nblr(jj))) = .true.
        end if
      end do
      do jj=1,rs_nbl(rs1)%nwnbs
        if (chgflag(rs_nbl(rs1)%wnb(jj)).EQV..true.) then
          thr_lhlper(rsx,cglst%irsl(rs_nbl(rs1)%wnb(jj))) = .true.
        end if
      end do
      do jj=1,rs_nbl(rs1)%nwnbtrs
        if (chgflag(rs_nbl(rs1)%wnbtr(jj)).EQV..true.) then
          thr_lhlper(rsx,cglst%irsl(rs_nbl(rs1)%wnbtr(jj))) = .true.
        end if
      end do
      do jj=1,rs_nbl(rs1)%ngnbs
        if (chgflag(rs_nbl(rs1)%gnb(jj)).EQV..true.) then
          thr_lhlper(rsx,cglst%irsl(rs_nbl(rs1)%gnb(jj))) = .true.
        end if
      end do
      do jj=1,rs_nbl(rs1)%ngnbtrs
        if (chgflag(rs_nbl(rs1)%gnbtr(jj)).EQV..true.) then
          thr_lhlper(rsx,rs_nbl(rs1)%gnbtr(jj)) = .true.
        end if
      end do
    end do
!$OMP BARRIER
    if (thr_dlb(5,1).gt.0) then
      if (tpi.eq.1) thr_dlb(5,2) = thr_dlb(5,2) + 1
      call System_Clock(count=ttimer)
      thr_timings(11,tpi) = thr_timings(11,tpi) + ttimer
    end if
    do rsx=thr_limits(27,tpi),thr_limits(28,tpi)
      rs1 = cglst%rsl(rsx)
      imol = molofrs(rs1)
      do rsy=rsx+1,cglst%ncrs
        rs2 = cglst%rsl(rsy)
        if (cycle_frz.EQV..true.) then
          jmol = molofrs(rs2)
          if (imol.eq.jmol) then
            if (molfrzidx(imol).ge.3) cycle
            if ((molfrzidx(imol).eq.2).AND.(nchi(rs1).eq.0).AND.(nchi(rs2).eq.0)) cycle
          else
            if ((molfrzidx(imol).eq.4).AND.(molfrzidx(jmol).eq.4)) cycle
          end if
        end if
        if ((thr_lhlper(rsx,rsy).EQV..false.).AND.(thr_lhlper(rsy,rsx).EQV..false.)) then
          rs_nbl(rs1)%nnblrs = rs_nbl(rs1)%nnblrs + 1
          if (rs_nbl(rs1)%nnblrs.gt.rs_nbl(rs1)%lralsz) call nbl_resz(rs1,3)
          rs_nbl(rs1)%nnblrats = rs_nbl(rs1)%nnblrats + at(rs2)%nncgrps
          rs_nbl(rs1)%nblr(rs_nbl(rs1)%nnblrs) = rs2
          thr_lhlper(rsx,rsy) = .true.
        end if
      end do
    end do
    if (thr_dlb(5,1).gt.0) then
      call System_Clock(count=ttimer)
      thr_timings(12,tpi) = thr_timings(12,tpi) + ttimer
    end if
  end if
!
end
!
#endif
!
!
!-----------------------------------------------------------------------------
!
subroutine contactpairs(rs,cutdis)
!
  use iounit
  use sequen
  use cutoffs
!
  implicit none
!
  integer i,ii,kk,rs
  RTYPE dis,dis2,cutdis
!
  kk = refat(rs)
!
  do i=1,nseq
    ii = refat(i)
    call dis_bound(ii,kk,dis2)
    dis = sqrt(dis2)
    if (dis.lt.(cutdis+resrad(rs)+resrad(i))) then
      if (i.ge.rs) then
        rsp_mat(i,rs) = 1
      else
        rsp_mat(rs,i) = 1
      end if
    end if
  end do
!
end
!
!--------------------------------------------------------------------------------
!
subroutine clear_rsp()
!
  use cutoffs
  use sequen
!
  implicit none
!
  rsp_mat(:,:) = 0
!
end
!
!--------------------------------------------------------------------------------
!
subroutine clear_rsp2(rs)
!
  use cutoffs
  use sequen
!
  implicit none
!
  integer rs
!
  rsp_mat(rs,:) = 0
  rsp_mat(:,rs) = 0
!
end
!
!---------------------------------------------------------------------------------
!
subroutine cutoff_check(tpi)
!
  use iounit
  use atoms
  use mcgrid
  use polypep
  use sequen
  use cutoffs
  use energies, ONLY: ideal_run
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer i,j,ij,ik,jj,kk,stas(2),stos(2),sta,sto
  RTYPE dis2,svec(3)
  integer(KIND=8):: cnt1,cnt2,cnt3
  logical atrue
!
  if (ideal_run.EQV..true.) return ! no nonbonded potentials
!
  cnt1 = 0
  cnt2 = 0
  cnt3 = 0
  atrue = .true.
!
  ij = 1
  ik = nseq
  j = nseq
  i = 1
  kk = 1
  jj = 0
  sta=  0
  call respairs_new(ij,ik,kk,i,j,jj,atrue,sta,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
  if (tpi.gt.0) then
    kk = 1
    jj = 2
    call get_thread_loop_bounds(kk,jj,stas,stos,tpi)
    do jj=stas(1),stos(1)
      sta = 1
      sto = thr_rsp_nbl(jj)%trnrs
      if (jj.eq.stas(1)) sta = stas(2)
      if (jj.eq.stos(1)) sto = stos(2)
      do j=sta,sto
        call dis_bound_rs(thr_rsp_nbl(jj)%tr(j,1),thr_rsp_nbl(jj)%tr(j,2),svec)
        do ij=rsinfo(thr_rsp_nbl(jj)%tr(j,1),1),rsinfo(thr_rsp_nbl(jj)%tr(j,1),1)+rsinfo(thr_rsp_nbl(jj)%tr(j,1),2)
          do ik=rsinfo(thr_rsp_nbl(jj)%tr(j,2),1),rsinfo(thr_rsp_nbl(jj)%tr(j,2),1)+rsinfo(thr_rsp_nbl(jj)%tr(j,2),2)
            dis2 = (x(ik) - x(ij) + svec(1))**2 + (y(ik) - y(ij) + svec(2))**2 + (z(ik) - z(ij) + svec(3))**2
            cnt1 = cnt1 + 1
            if (dis2.le.mcnb_cutoff2) cnt2 = cnt2 + 1
            if (dis2.le.mcel_cutoff2) cnt3 = cnt3 + 1
          end do
        end do
      end do
    end do
  else
    kk = 1
    do j=1,rsp_nbl(kk)%trnrs
    end do
  end if
#else
  kk = 1
  do j=1,rsp_nbl(kk)%trnrs
  end do
#endif
!
end
!
!---------------------------------------------------------------------------------
!
