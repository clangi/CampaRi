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
! CONTRIBUTIONS: Rohit Pappu                                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
!
!             ##################################
!             #                                #
!             #  A SET OF SUBROUTINES WHICH    #
!             #  DEAL WITH INTERNAL COORDINATE #
!             #  SPACE (Z-MATRIX) AND MOLEC.   #
!             #  TOPOLOGY IN GENERAL           #
!             #                                #
!             ##################################
!
!
!-----------------------------------------------------------------------------
!
! this subroutine appends the Z-matrix by adding the atom specified
! according to the input coordinate values as well as the input ref.
! atoms
!
subroutine z_add(bti,z1,z2,z3,iz1,iz2,iz3,iz4)
!
  use atoms
  use iounit
  use params
  use zmatrix
!
  implicit none
!
  integer bti,iz1,iz2,iz3,iz4
  RTYPE z1,z2,z3
!
!
  if (bti.gt.n_biotyp) then
!
    write(ilog,*) 'Fatal. Requested unknown biotype in z_add(...)&
 &or sidechain(...). Offending type is ',bti,' (atom #',n,').'
    call fexit()
!
  else if (bti.gt.0) then
!
    attyp(n) = bio_ljtyp(bti)
    b_type(n) = bti
    if ((attyp(n).lt.0).OR.(attyp(n).gt.n_ljtyp)) then
      write(ilog,*) 'Fatal. Encountered unknown LJ-type in z_add(...&
 &). Offending type is ',attyp(n),' for atom #',n,'.'
      call fexit()
    else if (attyp(n).eq.0) then
      write(ilog,*) 'Fatal. Parameter file lacks support for biotype&
 & ',b_type(n),' (atom #',n,').'
      call fexit()
    else
      atnam(n) = lj_symbol(attyp(n))
      mass(n) = lj_weight(attyp(n))
    end if
!
    blen(n) = z1
    bang(n) = z2
    ztor(n) = z3
    if (ztor(n).lt.-180.0) then
      ztor(n) = ztor(n) + 360.0
    else if (ztor(n).gt.180.0) then
      ztor(n) = ztor(n) - 360.0
    end if
    iz(1,n) = iz1
    iz(2,n) = iz2
    iz(3,n) = iz3
    iz(4,n) = iz4
!   we have to increase the number of atoms by 1
    n = n + 1
!
! we use the -1-signal to indicate a ring closure
!
  else if (bti.eq.-1) then
     nadd = nadd + 1
     iadd(1,nadd) = iz1
     iadd(2,nadd) = iz2
!
  else
!
    write(ilog,*) 'Fatal. Received unsupported (less or equal to zero excepting -1) biotype in&
 & z_add(...). This is a bug.'
    call fexit()
!
  end if
!
end
!
!------------------------------------------------------------------------------------
!
! the following two subroutines populate the Z-matrix differently: by extracting
! the internals from the Cartesian coordinates
!
subroutine genzmat(imol)
!
  use atoms
  use iounit
  use zmatrix
  use molecule
  use pdb, ONLY: pdb_force_box
!
  implicit none
!
  integer i,i1,i2,i3,i4,i5,imol
  integer iz0(0:n),iz1(n)
  RTYPE getbang,getztor,dasign,imshfvec(3),blen2
!
!
! zero out the values of the local defining atoms
!
  i1 = 0
  i2 = 0
  i3 = 0
  i4 = 0
  i5 = 0
  iz0(0) = 0
  do i = atmol(imol,1),atmol(imol,2)
    iz0(i) = 0
    iz1(i) = 0
  end do
!
! zero out the final internal coordinate values
!
  do i = atmol(imol,1),atmol(imol,2)
    blen(i) = 0.0d0
    bang(i) = 0.0d0
    ztor(i) = 0.0d0
  end do
!
! first, decide which of the atoms to define next
!
  do i = atmol(imol,1),atmol(imol,2)
    i1 = i
!
! define the bond length for the current atom
!
    if (i.ge.(atmol(imol,1)+1)) then
      i2 = iz(1,i1)
      if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
!       for input that was forced into the box, guess the correct image shift vector per atom
!       later rebuild and shift globally by COG/M
        imshfvec(:) = 0.0
        call dis_bound3t(i2,i1,imshfvec,blen2)
        blen(i1) = sqrt(blen2)
        x(i1) = x(i1) + imshfvec(1)
        y(i1) = y(i1) + imshfvec(2)
        z(i1) = z(i1) + imshfvec(3)
      else
        blen(i1) = sqrt((x(i1)-x(i2))**2 + (y(i1)-y(i2))**2 + (z(i1)-z(i2))**2)
      end if
    end if
!
! define the bond angle for the current atom
!
    if (i.ge.(atmol(imol,1)+2)) then
      i3 = iz(2,i1)
      bang(i1) = getbang(i1,i2,i3)
    end if
!
! decide whether to use a dihedral or second bond angle;
! then find the value of the angle
!
    if (i.ge.(atmol(imol,1)+3)) then
      if (iz(4,i1) .eq. 0) then
        i4 = iz(3,i1)
        i5 = 0
        ztor(i1) = getztor(i1,i2,i3,i4)
      else
        i4 = iz(3,i1)
        ztor(i1) = getbang(i1,i2,i4)
        i5 = 1
        dasign = getztor(i1,i2,i3,i4)
        if (dasign .gt. 0.0d0)  i5 = -1
      end if
    end if
!
! transfer defining atoms to the permanent array iz
!
    iz(1,i1) = iz0(i2)
    iz(2,i1) = iz0(i3)
    iz(3,i1) = iz0(i4)
    iz(4,i1) = i5
    iz0(i1) = i
    iz1(i1) = i2
  end do
!
end
!
!------------------------------------------------------------------------------------
!
! this modified routine is only used in readpdb.f when it 
! has to be able to correctly reflect varying geometries ...
! it also employs a bunch of sanity checks which allow it to probe a system for
! which bond lengths and angles are ill-defined initially, and to replace the latter
! with CAMPARI default values (useful for partial re-building)
! if the input is somehow ill-defined, the routine is quite fickle because of 
! dpendencies in the Z matrix that are difficult to circumvent with a simple
! and universal solution
!
subroutine genzmat_f(imol,steppi)
!
  use atoms
  use iounit
  use zmatrix
  use molecule
  use pdb
  use system
  use fyoc
  use params
  use sequen
  use polypep, ONLY: oi,ni,hni,cai,ci,at
!
  implicit none
!
  integer i,i1,i2,i3,i4,i5,imol,steppi,rs,i5bu,nrebuild,k
  integer iz0(0:n),iz1(n)
  RTYPE getbang_safe,getztor,dasign,blbu,babu,dhbu,anv,bll,blh,imshfvec(3),blen2,getblen,getbang,bamax
  logical flagi1,flagi1f,frameprob,flagi2,overelig,fndem(10)
  integer, ALLOCATABLE:: rebuildix(:)
!
! some initialization
  if (steppi.eq.1) then
    anv = pdb_tolerance(3)
    bll = pdb_tolerance(1)
    blh = pdb_tolerance(2)
! if we're referencing back to the last structures, the backup-
! values won't necessarily hold values close to the minimum -> bump up tolerance
  else if (pdb_analyze.EQV..true.) then
    anv = 2.0*pdb_tolerance(3)
    bll = max(1.0 - 2.0*(1.0-pdb_tolerance(1)),0.0d0)
    blh = max(1.0 - 2.0*(1.0-pdb_tolerance(2)),0.0d0)
  end if
!
! zero out the values of the local defining atoms
!
  i1 = 0
  i2 = 0
  i3 = 0
  i4 = 0
  i5 = 0
  iz0(0) = 0
  allocate(rebuildix(atmol(imol,2)-atmol(imol,1)+1))
  nrebuild = 0
!
  do i = atmol(imol,1),atmol(imol,2)
    iz0(i) = 0
    iz1(i) = 0
  end do
  frameprob = .false.
  overelig = .false.
  if (rsmol(imol,1).eq.rsmol(imol,2)) then
    if ((seqtyp(rsmol(imol,1)).eq.39).OR.(seqtyp(rsmol(imol,1)).eq.40).OR.&
 &      (seqtyp(rsmol(imol,1)).eq.45).OR.(seqtyp(rsmol(imol,1)).eq.46).OR.&
 &      (seqtyp(rsmol(imol,1)).eq.56).OR.(seqtyp(rsmol(imol,1)).eq.102).OR.&
 &      (seqtyp(rsmol(imol,1)).eq.103)) overelig = .true.
  end if
!
! first, decide which of the atoms to define next
!
  do i = atmol(imol,1),atmol(imol,2)
    blbu = blen(i)
    blen(i) = 0.0
    babu = bang(i)
    bang(i) = 0.0
    dhbu = ztor(i)
    ztor(i) = 0.0
    i1 = i
    flagi1 = .false.
    flagi1f = .false.
    flagi2 = .false.
    i5bu = iz(4,i1)
    i5 = 0
    rs = atmres(i1)
!
!   always flag the first PDB-derived HN in a polypeptide residue if an N-terminal stretch was fused to it 
!   (does not reset automatically due to Z matrix structure HN(rs) <- N(rs) <- CA(rs) <- C(rs)) 
    if (allocated(pdb_didread).EQV..true.) then
      if (size(pdb_didread,dim=1).ge.imol) then
        if ((pdb_didread(imol,1).gt.rsmol(imol,1)).AND.(atmres(i1).eq.pdb_didread(imol,1)).AND.(i1.eq.hni(atmres(i1)))) then
          flagi1 = .true.
        end if
      end if
    end if
!   define the bond length for the current atom
    if (i.ge.(atmol(imol,1)+1)) then
      i2 = iz(1,i1)
      if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
!       for input that was forced into the box, guess the correct image shift vector per atom
!       later rebuild and shift globally by COG/M
        imshfvec(:) = 0.0
        call dis_bound3t(i2,i1,imshfvec,blen2)
        blen(i1) = sqrt(blen2)
        x(i1) = x(i1) + imshfvec(1)
        y(i1) = y(i1) + imshfvec(2)
        z(i1) = z(i1) + imshfvec(3)
      else
        blen(i1) = sqrt((x(i1)-x(i2))**2 + (y(i1)-y(i2))**2 + (z(i1)-z(i2))**2)
      end if
      if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
        if(((blen(i1)/blbu).gt.blh).OR.((blen(i1)/blbu).lt.bll))then
          if (steppi.eq.1) write(ilog,*) 'Bond length exception at atoms ',i1,'and ',i2,'.'
          if (i.ge.(atmol(imol,1)+3)) then
            flagi1 = .true.
          else
            flagi2 = .true.
!            write(ilog,*) 'WARNING: Atom ',i1,' will NOT BE ADJUSTED!'
          end if
        end if
        if ((pdb_hmode.eq.2).AND.(mass(i).lt.3.5)) then
          flagi1 = .true.
        end if
      end if
    end if
!
! define the bond angle for the current atom
!
    if (i.ge.(atmol(imol,1)+2)) then
      i3 = iz(2,i1)
      bang(i1) = getbang_safe(i1,i2,i3)
      if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
        if ((abs(bang(i1)-babu).gt.anv).OR.(bang(i1).lt.0.0)) then
          bang(i1) = max(bang(i1),0.0)
          if (steppi.eq.1) write(ilog,*) 'Bond angle exception at atoms ',i1,'and ',i2,'and ',i3,'.'
          if (i.ge.(atmol(imol,1)+3)) then
            flagi1 = .true.
          else
            flagi2 = .true.
!            write(ilog,*) 'WARNING: Atom ',i1,' will NOT BE ADJUSTED!'
          end if
        end if
        if ((pdb_hmode.eq.2).AND.(mass(i).lt.3.5)) then
          flagi1 = .true.
        end if
      else if (bang(i1).lt.0.0) then
        bang(i1) = max(bang(i1),0.0)
      end if
    end if
!
! decide whether to use a dihedral or second bond angle;
! then find the value of the angle
!
    if (i.ge.(atmol(imol,1)+3)) then
      if (iz(4,i1).eq.0) then
        i4 = iz(3,i1)
        i5 = 0
        ztor(i1) = getztor(i1,i2,i3,i4)
!        if (i1.eq.hni(atmres(i1))) write(*,*) i1,ztor(i1),ztor(fline(atmres(i1))),flagi1
        if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
          if (((pdb_hmode.eq.2).OR.(flagi1.EQV..true.)).AND.(mass(i).lt.3.5)) then
            if (fline2(rs).eq.i1) then
              flagi1f = .true.
              if (ztor(fline(rs)).gt.0.0) then
                ztor(i1) = ztor(fline(rs)) - 180.0
              else
                ztor(i1) = ztor(fline(rs)) + 180.0
              end if
            else
              flagi1 = .true.
            end if
          end if
!
!         a problem with the y(f)line/y(f)line2 setup is that an actual dihedral
!         which is NOT a degree of freedom will be used
!         if it is ever requested to correct oxygen positions, use construct below (untested)
!
!          if (whatever) then
!            else if ((rs.gt.rsmol(imol,1)).AND.
! &                          (yline2(rs-1).eq.i1)) then
!              write(*,*) 'in with ',i1
!              if (ztor(yline(rs-1)).gt.0.0) then
!                ztor(i1) = ztor(yline(rs-1)) - 180.0
!              else
!                ztor(i1) = ztor(yline(rs-1)) + 180.0
!              end if
!            end if
!          end if
        end if
      else
        i4 = iz(3,i1)
        ztor(i1) = getbang_safe(i1,i2,i4)
        i5 = 1
        dasign = getztor(i1,i2,i3,i4)
        if (dasign .gt. 0.0d0)  i5 = -1
        if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
          if ((abs(ztor(i1)-dhbu).gt.anv).OR.(ztor(i1).lt.0.0)) then
            ztor(i1) = max(ztor(i1),0.0)
            if (steppi.eq.1) write(ilog,*) 'Bond angle exception at atoms ',i1,'and ',i2,'and ',i4,'.'
            if (i.ge.(atmol(imol,1)+3)) then
              flagi1 = .true.
            else
              flagi2 = .true.
!              write(ilog,*)'WARNING: Atom ',i1,' will NOT BE ADJUSTED!'
            end if
          end if
          if ((pdb_hmode.eq.2).AND.(mass(i).lt.3.5)) then
            flagi1 = .true.
          end if
        else if (ztor(i1).lt.0.0) then
          ztor(i1) = max(ztor(i1),0.0)
        end if
      end if
    end if
!
!   for a re-builder, restore every(!)thing except fline2 torsion (which is already set)
    if (flagi1.EQV..true.) then
      blen(i1) = blbu
      bang(i1) = babu
      if (flagi1f.EQV..false.) ztor(i1) = dhbu
      i5 = i5bu
    end if
!
!   atom within reference frame that has exception is only handled for special cases (otherwise everything shifts)
    if ((flagi2.EQV..true.).AND.(overelig.EQV..true.)) then
      frameprob = .true.
      blen(i1) = blbu
      bang(i1) = babu
    else
      flagi2 = .false.
    end if
!
!   force re-building message for phi-case
    if (flagi1f.EQV..true.) flagi1 = .true.
!
! transfer defining atoms to the permanent array iz
!
    iz(1,i1) = iz0(i2)
    iz(2,i1) = iz0(i3)
    iz(3,i1) = iz0(i4)
    iz(4,i1) = i5
    iz0(i1) = i
    iz1(i1) = i2
    if (((flagi1.EQV..true.).OR.(flagi2.EQV..true.)).AND.(steppi.eq.1)) then
      write(ilog,*) 'Re-building atom: ',i1,' (ty: ',b_type(i),' = '&
 &,bio_code(b_type(i)),').'
    end if
    if ((flagi1.EQV..true.).OR.(flagi2.EQV..true.)) then
!     if part of our reference frame is missing or we miss an atom occurring more than once in the Z matrix,
!     we will most likely have a problem which we deal with below
      nrebuild = nrebuild + 1
      rebuildix(nrebuild) = i1
    end if
  end do
!
! in this loop we check whether atoms being rebuilt rely on other such atoms and we provide the necessary information
  do rs=1,nrebuild
    i1 = rebuildix(rs)
    do i=i1+1,atmol(imol,2)
      if ((iz(1,i).eq.i1).OR.(iz(2,i).eq.i1).OR.(iz(3,i).eq.i1)) then
        do k=rs+1,nrebuild
          if (i.eq.rebuildix(k)) then
            exit
          else if (rebuildix(k).gt.i) then
!           in order to break the downstream behavior, we need to get the rebuilt coordinates, and reset the relevant
!           Z matrix entry
!           however, we can do this only if this does not create new problems (which would become a recursive problem)
            call genxyz(i1,iz(1,i1),blen(i1),iz(2,i1),bang(i1),iz(3,i1),ztor(i1),iz(4,i1))
            if (iz(3,i).eq.i1) then
              if (iz(4,i).ne.0) then
                babu = ztor(i)
                ztor(i) = getbang(i1,iz(1,i),i)
                if ((abs(ztor(i)-babu).gt.anv).OR.(ztor(i).lt.0.0)) then
                  if (steppi.eq.1) write(ilog,*) 'Warning. The rebuilding of a poorly placed or missing atom (',i1,') has caused &
 &a secondary problem (bond angle). Inspect input structure for position of atom ',i,'.'
                  ztor(i) = babu
                  exit
                end if
              else
                ztor(i) = getztor(i1,iz(2,i),iz(1,i),i)
              end if
            else if (iz(2,i).eq.i1) then
              babu = bang(i)
              bang(i) = getbang(i1,iz(1,i),i)
              if ((abs(bang(i)-babu).gt.anv).OR.(bang(i).lt.0.0)) then
                if (steppi.eq.1) write(ilog,*) 'Warning. The rebuilding of a poorly placed or missing atom (',i1,') has caused a &
 &secondary problem (bond angle). Inspect input structure for position of atom ',i,'.'
                bang(i) = babu
                exit
              end if
              if (iz(4,i).ne.0) then
!               do nothing (position of i1 is irrelevant for second bond angle)
              else
                ztor(i) = getztor(iz(3,i),i1,iz(1,i),i)
              end if
            else
              blbu = blen(i)
              blen(i) = getblen(i1,i)
              if (((blen(i)/blbu).gt.blh).OR.((blen(i)/blbu).lt.bll)) then
                if (steppi.eq.1) write(ilog,*) 'Warning. The rebuilding of a poorly placed or missing atom (',i1,') has caused a &
 &secondary problem (bond length). Inspect input structure for position of atom ',i,'.'
                blen(i) = blbu
                exit
              end if
              babu = bang(i)
              bang(i) = getbang(iz(2,i),i1,i)
              if ((abs(bang(i)-babu).gt.anv).OR.(bang(i).lt.0.0)) then
                if (steppi.eq.1) write(ilog,*) 'Warning. The rebuilding of a poorly placed or missing atom (',i1,') has caused a &
 &secondary problem (bond angle). Inspect input structure for position of atom ',i,'.'
                bang(i) = babu
                exit
              end if
              if (iz(4,i).ne.0) then
                babu = ztor(i)
                ztor(i) = getbang(iz(3,i),i1,i)
                if ((abs(ztor(i)-babu).gt.anv).OR.(ztor(i).lt.0.0)) then
                  if (steppi.eq.1) write(ilog,*) 'Warning. The rebuilding of a poorly placed or missing atom (',i1,') has caused a &
 &secondary problem (bond angle). Inspect input structure for position of atom ',i,'.'
                  ztor(i) = babu
                  exit
                end if
              else
                ztor(i) = getztor(iz(3,i),i1,iz(1,i),i)
              end if
            end if
            exit
          end if
        end do
      end if
    end do
  end do
! an indirect way of backward reliance are the (aaarghh) fline2/yline2 constructs and through wanting
! certain dihedral angles to be planar
  do rs=1,nrebuild
    i1 = rebuildix(rs)
    fndem(:) = .false.
    if ((i1.eq.oi(atmres(i1))).AND.(atmres(i1).lt.rsmol(molofrs(atmres(i1)),2))) then
      do k=1,nrebuild
        if (rs.eq.k) cycle
        if (rebuildix(k).eq.ni(atmres(i1)+1)) then
          fndem(1) = .true.
          ztor(oi(atmres(i1))) = -0.5
          ztor(ni(atmres(i1)+1)) = 179.5
        end if
      end do
      if (fndem(1).EQV..false.) then
        ztor(oi(atmres(i1))) = ztor(ni(atmres(i1)+1)) - 180.
      end if
    else if (i1.eq.ni(atmres(i1))) then
      do k=1,nrebuild
        if (rs.eq.k) cycle
        if (rebuildix(k).eq.hni(atmres(i1))) then
          fndem(6) = .true.
        else if (rebuildix(k).eq.cai(atmres(i1))) then
          fndem(7) = .true.
        end if
        if (atmres(i1).gt.rsmol(molofrs(atmres(i1)),1)) then
          if (rebuildix(k).eq.oi(atmres(i1)-1)) then
            fndem(8) = .true.
          end if
        end if
      end do
      if ((fndem(6).EQV..false.).AND.(hni(atmres(i1)).gt.0).AND.(ci(atmres(i1)).gt.0)) then
        ztor(hni(atmres(i1))) = ztor(ci(atmres(i1))) - 180.
      end if
      if ((fndem(7).EQV..false.).AND.(cai(atmres(i1)).gt.0)) then
        ztor(cai(atmres(i1))) = 180. ! restore planar amide group
      end if
      if ((fndem(8).EQV..false.).AND.(atmres(i1).gt.rsmol(molofrs(atmres(i1)),1))) then
        if ((oi(atmres(i1)-1).gt.0).AND.(ni(atmres(i1)).gt.0)) then
          ztor(oi(atmres(i1)-1)) = ztor(ni(atmres(i1))) - 180.
        end if
      end if
    else if (i1.eq.cai(atmres(i1))) then
      ztor(cai(atmres(i1))) = 180. ! restore planar amide group
    else if ((seqtyp(atmres(i1)).eq.30).AND.(i1.eq.at(atmres(i1))%bb(1)+1)) then
      ztor(i1) = 180.
    end if
  end do
! secondly, we check whether pairwise hydrogens rebuilt this way will cause a bond angle violation leading to
! overlapping atoms
  do i=atmol(imol,1),atmol(imol,2)-1
    if (((abs(iz(4,i)).eq.1).AND.(iz(4,i)*iz(4,i+1).eq.-1)).AND.(iz(3,i).eq.iz(3,i+1)).AND.(iz(2,i).eq.iz(2,i+1)).AND.&
 &      (iz(1,i).eq.iz(1,i+1))) then
      flagi1f = .false.
      do k=1,nrebuild
        if ((rebuildix(k).eq.iz(3,i)).OR.(rebuildix(k).eq.iz(2,i)).OR.(rebuildix(k).eq.iz(1,i))) then
          flagi1f = .true.
          exit
        end if
      end do
      if (flagi1f.EQV..true.) cycle ! we must not do this if a reference atom is missing
      bamax = (360.0 - getbang(iz(2,i),iz(1,i),iz(3,i)))
      if (bang(i)+bang(i+1).gt.bamax) then
        bamax = bamax/(bang(i)+bang(i+1))/0.9 ! the 0.9 is an empirical target ratio for (a_1 + a_2) / (360.0 - a_ref)
        bang(i) = bang(i)/bamax
        bang(i+1) = bang(i+1)/bamax
      end if
    end if
  end do
!
! the rather special case of proteins missing the second oxygen atom to generate a C-term. can lead to
! problems if the "wrong" right name (1OXT or 2OXT) is assigned to the lone O -> fix
  do k=1,nrebuild
    i1 = rebuildix(k)
    rs = atmres(i1)
    if (((i1.eq.yline(rs)).OR.(i1.eq.yline2(rs))).AND.(rs.eq.rsmol(imol,2)).AND.(seqpolty(rs).eq.'P')) then
      if ((yline(rs).gt.0).AND.(yline2(rs).gt.0)) then
        if ((b_type(yline(rs))).eq.(b_type(yline2(rs)))) then
          if (i1.eq.yline(rs)) then
            ztor(i1) = ztor(yline2(rs))
            if (ztor(i1).lt.0.0d0) then
              ztor(i1) = ztor(i1) + 180.0d0
            else
              ztor(i1) = ztor(i1) - 180.0d0
            end if
          else
            ztor(i1) = ztor(yline(rs))
            if (ztor(i1).lt.0.0d0) then
              ztor(i1) = ztor(i1) + 180.0d0
            else
              ztor(i1) = ztor(i1) - 180.0d0
            end if
          end if
        end if
      end if
    end if
  end do
!
  if ((steppi.eq.1).OR.(pdb_analyze.EQV..true.)) then
    do rs=rsmol(imol,1),rsmol(imol,2)
!     because of the way, phi- and psi-angles are coded up, there are two associated
!     ztor()-entries for each
!     in case the geometry comes from pdb we cannot assume what their (static) relationship
!     is and need to determine it explicitly
      if ((fline(rs).gt.0).AND.(fline2(rs).gt.0)) then
        phish(rs) = ztor(fline2(rs)) - ztor(fline(rs))
        if (phish(rs).gt.180.0) phish(rs) = phish(rs) - 360.0
        if (phish(rs).le.-180.0) phish(rs) = phish(rs) + 360.0
      end if
      if ((yline(rs).gt.0).AND.(yline2(rs).gt.0)) then 
        if (ztor(yline(rs)).lt.0.0) then
          psish(rs) = ztor(yline2(rs)) - (ztor(yline(rs)) + 180.0)
        else
          psish(rs) = ztor(yline2(rs)) - (ztor(yline(rs)) - 180.0)
        end if
        if (psish(rs).gt.180.0) psish(rs) = psish(rs) - 360.0
        if (psish(rs).le.-180.0) psish(rs) = psish(rs) + 360.0
      end if
    end do
  end if
  if ((pdb_ihlp.le.0).AND.(nrebuild.gt.0).AND.(steppi.gt.1)) then
    write(ilog,'(1x,a,i10,a)') 'Warning. One or more atoms were unexpectedly rebuilt for PDB analysis snapshot #',steppi,'.'
    write(ilog,'(1x,a)') 'This is most likely related to overly stringent settings for FMCSC_PDB_TOLERANCE_A and &
 &FMCSC_PDB_TOLERANCE_B.'
    write(ilog,*) 'All further warnings reporting on rebuilt atoms are omitted.'
  end if
  if ((pdb_ihlp.lt.(HUGE(pdb_ihlp)-nrebuild)).AND.(steppi.ge.1)) then
    pdb_ihlp = pdb_ihlp + nrebuild
  end if
  if ((nrebuild.gt.0).AND.(steppi.eq.1).AND.(pdb_analyze.EQV..true.)) then
    write(ilog,*) 'All further warnings reporting on rebuilt atoms are omitted.'
  end if
!
  call zmatfyc()
  if (frameprob.EQV..true.) then
    call makexyz_formol3(imol)
  else
    call makexyz_formol(imol)
  end if
!
  deallocate(rebuildix)
!
end
!
!----------------------------------------------------------------------------------
!
! a finally a routine which prints the Z-matrix 
!
subroutine prtzmat(iu)
!
  use atoms
  use molecule
  use zmatrix
  use iounit
!
  implicit none
!
  integer i,k,iu,nnn,kk,imol
  logical isopen
!
  inquire (unit=iu,opened=isopen)
  if (isopen.EQV..false.) then
    write(ilog,*) 'Fatal. Got a stale filehandle in prtzmat(...).'
    call fexit()
  end if
!
! these formattings are chosen to be compatible with TINKER (as much as possible
! since molecule handling is different) 
!
 51   format (i6)
 52   format (i6,2x,a3,i5)
 53   format (i6,2x,a3,i5,i6,f10.5)
 54   format (i6,2x,a3,i5,i6,f10.5,i6,f10.4)
 55   format (i6,2x,a3,i5,i6,f10.5,i6,f10.4,i6,f10.4,i6)
 56   format (2i6)
!
  write(iu,51) n
  do imol=1,nmol
!
!   first three are special cases
!
    nnn = atmol(imol,2)-atmol(imol,1)+1
    if (nnn.ge.1) then
      kk = atmol(imol,1)
      write(iu,52)  kk,atnam(kk),attyp(kk)
    end if
    if (nnn.ge.2) then
      kk = atmol(imol,1) + 1
      write(iu,53)  kk,atnam(kk),attyp(kk),iz(1,kk),blen(kk)
    end if
    if (nnn.ge.3) then
      kk = atmol(imol,1) + 2
      write(iu,54)  kk,atnam(kk),attyp(kk),iz(1,kk),blen(kk),&
 &                    iz(2,kk),bang(kk)
    end if
!
!   now standard rest
!
    do i=atmol(imol,1)+3,atmol(imol,2)
      write(iu,55)  i,atnam(i),attyp(i),iz(1,i),blen(i),&
 &                    iz(2,i),bang(i),iz(3,i),ztor(i),iz(4,i)
    end do
!
  end do
!
! ring closures
!
  if (nadd.ne.0)  write(iu,*)
  do i=1,nadd
    write(iu,56) (iadd(k,i),k=1,2)
  end do
!
end
!
!------------------------------------------------------------------------------
!
! a subroutine to build a working (and - if possible - chemically reasonable) Z-matrix from PDB for an
! unknown entity; also sets valency, atom types (through guess_types) and detects ring bonds
! works for intramolecular linkages as well
!
subroutine infer_topology(i,nl,unkbnd,xyzdum,d2max,pdbnm,rsnam)
!
  use pdb
  use zmatrix
  use sequen
  use molecule
  use math
  use params
  use atoms
  use polypep
  use iounit
  use aminos
  use system, ONLY: ua_model
!
  implicit none
!
  integer, INTENT(IN) :: i,nl
  integer, INTENT(INOUT) :: unkbnd(n_pdbunk,4)
  RTYPE xyzdum(3,nl),d2max(nl)
  character(4), INTENT(IN) :: pdbnm(nl)
  character(3), INTENT(IN) :: rsnam(nl)
  integer j,k,ii,kk,rs,dihk,knty,rcnt,startj,shf,ji,ki,imol,jmol,mshf,sta,sto
  integer, ALLOCATABLE:: invpdbmap(:),btyh(:)
  RTYPE d2,d2t
  logical dihthere,mustcopy
!
 77 format('Map extracted from PDB is incomplete.'/'Coordinates for atom ',i10,' of biotype ',a3,' in residue ',i8,' (',a3,') &
 &could not be read.'/'Some atoms in supported residues in immediate sequence context must be present and readable.')
  allocate(invpdbmap(max(n,nl)))
  invpdbmap(:) = 0
  rs = unkbnd(i,3)
  imol = molofrs(unkbnd(i,3))  
  sta = at(max(rsmol(imol,1),rs-1))%bb(1)
  sto = at(min(rsmol(imol,2),rs+1))%bb(1) + at(min(rsmol(imol,2),rs+1))%nbb + at(min(rsmol(imol,2),rs+1))%nsc - 1
!
  knty = 0
  mustcopy = .false.
  do j=1,i-1
    if (rsnam(unkbnd(i,1))(1:3).eq.rsnam(unkbnd(j,1))(1:3)) then
      if ((unkbnd(i,1).eq.unkbnd(j,1)).AND.(unkbnd(i,2).eq.unkbnd(j,2))) mustcopy = .true.
      knty = j
      unkbnd(i,4) = unkbnd(knty,3)
      exit
    end if
  end do
!
! if template does not actually provide explicit coordinates, no extra checks needed
  if (mustcopy.EQV..true.) then
!   implies single residue molecule
    jmol = molofrs(unkbnd(knty,3))
    mshf = atmol(imol,1)-atmol(jmol,1)
!   Z matrix
    iz(1:3,atmol(imol,1)) = 0
    do j=atmol(jmol,1)+1,atmol(jmol,2)
      iz(1,j+mshf) = iz(1,j) + mshf
      if (j.gt.(atmol(jmol,1)+1)) iz(2,j+mshf) = iz(2,j) + mshf
      if (j.gt.(atmol(jmol,1)+2)) iz(3,j+mshf) = iz(3,j) + mshf
    end do
    iz(4,atmol(imol,1):atmol(imol,2)) = iz(4,atmol(jmol,1):atmol(jmol,2))
    blen(atmol(imol,1):atmol(imol,2)) = blen(atmol(jmol,1):atmol(jmol,2))
    bang(atmol(imol,1):atmol(imol,2)) = bang(atmol(jmol,1):atmol(jmol,2))
    ztor(atmol(imol,1):atmol(imol,2)) = ztor(atmol(jmol,1):atmol(jmol,2))
    b_type(atmol(imol,1):atmol(imol,2)) = b_type(atmol(jmol,1):atmol(jmol,2))
    atnam(atmol(imol,1):atmol(imol,2)) = atnam(atmol(jmol,1):atmol(jmol,2))
    mass(atmol(imol,1):atmol(imol,2)) = mass(atmol(jmol,1):atmol(jmol,2))
    attyp(atmol(imol,1):atmol(imol,2)) = attyp(atmol(jmol,1):atmol(jmol,2))
    do j=1,nadd
      if ((iadd(1,j).ge.atmol(jmol,1)).AND.(iadd(1,j).le.atmol(jmol,2)).AND.&
 &        (iadd(2,j).ge.atmol(jmol,1)).AND.(iadd(2,j).le.atmol(jmol,2))) then
        nadd = nadd + 1
        iadd(1,nadd) = iadd(1,j) + mshf
        iadd(2,nadd) = iadd(2,j) + mshf
      end if
    end do
    return
  end if
!
  do j=sta,sto
    if (pdbmap(j).gt.0) then
      if (invpdbmap(pdbmap(j)).eq.0) invpdbmap(pdbmap(j)) = j
    else
      if ((iz(1,j).le.0).OR.(iz(2,j).le.0).OR.(iz(3,j).le.0)) then
        write(ilog,77) j,bio_code(b_type(j)),atmres(j),amino(seqtyp(atmres(j)))
        call fexit()
      end if
    end if
  end do
 !
  shf = 0
  if (ua_model.gt.0) shf = 1
!
! sanity
  if (knty.gt.0) then
    k = 0
    if ((unkbnd(i,2)-unkbnd(i,1)).ne.(unkbnd(knty,2)-unkbnd(knty,1))) then
      write(ilog,*) 'Fatal. When using multiple unsupported residues with the same residue name, &
 &it is mandatory that all those residues possess identical numbers and types of atoms in identical &
 &order.'
      call fexit()
    end if
    if (((unkbnd(knty,3).eq.rsmol(molofrs(unkbnd(knty,3)),1)).AND.(unkbnd(i,3).ne.rsmol(molofrs(unkbnd(i,3)),1))).OR.&
 &      ((unkbnd(knty,3).eq.rsmol(molofrs(unkbnd(knty,3)),2)).AND.(unkbnd(i,3).ne.rsmol(molofrs(unkbnd(i,3)),2))).OR.&
 &      ((unkbnd(knty,3).ne.rsmol(molofrs(unkbnd(knty,3)),2)).AND.(unkbnd(i,3).eq.rsmol(molofrs(unkbnd(i,3)),2))).OR.&
 &      ((unkbnd(knty,3).ne.rsmol(molofrs(unkbnd(knty,3)),2)).AND.(unkbnd(i,3).eq.rsmol(molofrs(unkbnd(i,3)),2)))) then
      write(ilog,*) 'Fatal. When using multiple unsupported residues with the same residue name, &
 &it is mandatory that these residues do not occur in both terminal and non-terminal positions.'
      call fexit()
    end if
    do j=unkbnd(i,1),unkbnd(i,2)
      if (pdbnm(j).ne.pdbnm(unkbnd(knty,1)+k)) then
        write(ilog,*) 'Fatal. When using multiple unsupported residues with the same residue name, &
 &it is mandatory that all those residues possess identical numbers and types of atoms in identical &
 &order.'
        call fexit()
      end if
      k = k + 1
    end do
  end if
!
  allocate(btyh(n_biotyp+unkbnd(i,2)-unkbnd(i,1)+1))
  btyh(:) = 0
  do j=unkbnd(i,1),unkbnd(i,2)
    d2 = HUGE(d2)
!   find the nearest atom listed prior that is part of the same molecule (partner for Z-matrix bond)
    startj = 0
    do k=j-1,1,-1
      if (invpdbmap(k).le.0) cycle
      if (atmres(invpdbmap(k)).eq.rsmol(molofrs(unkbnd(i,3)),1)) then
        startj = k
      else if (startj.gt.0) then
        exit
      end if
    end do
    if (startj.eq.0) startj = 1
    ii = 0
    do k=startj,j-1
      if (invpdbmap(k).le.0) cycle
      call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2t)
      if (d2t.lt.d2) then 
        d2 = d2t
        ii = k
      end if
    end do
!   get an estimate of the net valence (this is highly approximate of course and makes very many assumptions
!   about atomistic systems)
    kk = 0
    do k=1,nl
      if (k.eq.j) cycle
      if (invpdbmap(k).le.0) cycle
      call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2t)
      if ((k.ge.unkbnd(i,1)).AND.(k.le.unkbnd(i,2))) then
        if (d2t.lt.d2max(j)*d2max(k)) then
          kk = kk + 1
        end if
      else if (attyp(invpdbmap(k)).le.0) then
        if (d2t.lt.d2max(j)*d2max(k)) then
          kk = kk + 1
        end if
      else
        if (d2t.lt.d2max(j)*0.5*sqrt(lj_sig(attyp(invpdbmap(k)),attyp(invpdbmap(k))))) then
          kk = kk + 1
        end if
      end if
    end do
!
    if (knty.gt.0) then
      rcnt = at(unkbnd(i,4))%bb(1)+j-unkbnd(i,1)
      b_type(invpdbmap(j)) = b_type(rcnt)
      attyp(invpdbmap(j)) = attyp(rcnt)
      mass(invpdbmap(j)) = mass(rcnt)
    else
      dihthere = .false.
      do k=unkbnd(i,1),j-1
!       try to infer (stringently) equivalent atom types
        if ((pdbnm(j)(2:4).eq.pdbnm(k)(2:4)).AND.(iz(1,invpdbmap(k)).eq.invpdbmap(ii)).AND.&
 &((pdbnm(j)(1:1).eq.'1').OR.(pdbnm(j)(1:1).eq.'2').OR.(pdbnm(j)(1:1).eq.'3').OR.(pdbnm(j)(1:1).eq.'4')).AND.&
 &((pdbnm(k)(1:1).eq.'1').OR.(pdbnm(k)(1:1).eq.'2').OR.(pdbnm(k)(1:1).eq.'3').OR.(pdbnm(k)(1:1).eq.'4')).AND.&
 &(pdbnm(j)(1:1).ne.pdbnm(k)(1:1))) then
          b_type(invpdbmap(j)) = b_type(invpdbmap(k))
          btyh(b_type(invpdbmap(j))) = btyh(b_type(invpdbmap(j))) + 1
          attyp(invpdbmap(j)) = attyp(invpdbmap(k))
          mass(invpdbmap(j)) = mass(invpdbmap(k))
          dihthere = .true.
          exit
        end if
      end do
      if (dihthere.EQV..false.) then
        call guess_types(kk,mass(invpdbmap(j)),pdbnm(j),attyp(invpdbmap(j)),b_type(invpdbmap(j)))
        btyh(b_type(invpdbmap(j))) = btyh(b_type(invpdbmap(j))) + 1
!       mass should be matched to atom type (info about guessed mass would be lost -> undesirable)
        mass(invpdbmap(j)) = lj_weight(attyp(invpdbmap(j)))
      end if
    end if
!
    ji = invpdbmap(j)
!    write(*,*) attyp(j),': ',lj_symbol(attyp(j)),' and ',b_type(j),': ',bio_code(b_type(j)),' | ',kk
    iz(1:4,ji) = 0
    if ((j.eq.unkbnd(i,1)).AND.(rs.eq.rsmol(molofrs(rs),1))) then
      iz(1:4,ji) = 0
    else if ((j.eq.unkbnd(i,1)+1).AND.(rs.eq.rsmol(molofrs(rs),1))) then
      iz(1,ji) = invpdbmap(ii) ! i = i
      blen(ji) = sqrt(d2)
      iz(2:4,ji) = 0
    else
!     the angle reference is set up without any particular considerations (input atom order critical)
      iz(1,ji) = invpdbmap(ii)
      blen(ji) = sqrt(d2)
      if (iz(1,iz(1,ji)).gt.0) then
        iz(2,ji) = iz(1,iz(1,ji))
        call bondang2(xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),bang(ji))
        bang(ji) = bang(ji)*RADIAN
      else
        do k=1,j-1
          if (invpdbmap(k).le.0) cycle
          if (iz(1,invpdbmap(k)).eq.iz(1,ji)) exit
        end do
        iz(2,ji) = invpdbmap(k)
        call bondang2(xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),bang(ji))
        bang(ji) = bang(ji)*RADIAN
      end if
      if ((j.eq.unkbnd(i,1)+2).AND.(rs.eq.rsmol(molofrs(rs),1))) then
        iz(3:4,ji) = 0
      else
!       we should scan whether we already have an atom utilizing the same bond iz(1,ji),iz(2,ji)
!       if not, this will become a putative d.o.f-type dihedral, otherwise create indirect reference to fix
!       covalent geometry
        dihthere = .false.
        do k=1,j-1
          if (invpdbmap(k).le.0) cycle
          ki = invpdbmap(k)
          if (((iz(1,ki).eq.iz(2,ji)).AND.(iz(2,ki).eq.iz(1,ji))).OR.&
 &            ((iz(1,ki).eq.iz(1,ji)).AND.(iz(2,ki).eq.iz(2,ji)))) then
            dihk = ki
            dihthere = .true.
            exit
          end if
        end do
        if (dihthere.EQV..false.) then
          if ((iz(1,iz(2,ji)).gt.0).AND.(iz(1,iz(2,ji)).ne.iz(1,ji))) then
            iz(3,ji) = iz(1,iz(2,ji))
            call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
            ztor(ji) = ztor(ji)*RADIAN
            iz(4,ji) = 0
          else
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              ki = invpdbmap(k)
              if ((iz(1,ki).eq.iz(2,ji)).AND.(ki.ne.iz(1,ji))) exit
            end do
            iz(3,ji) = ki
            call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
            ztor(ji) = ztor(ji)*RADIAN
            iz(4,ji) = 0
          end if
          if (iz(3,ji).eq.0) then
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              ki = invpdbmap(k)
              if (((iz(1,ki).eq.iz(1,ji)).OR.(iz(1,iz(1,ji)).eq.ki)).AND.(ki.ne.iz(2,ji))) then
                iz(3,ji) = ki
                call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),&
 &                         xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
                iz(4,ji) = -1
                if (ztor(ji).lt.0.0) iz(4,ji) = 1
                call bondang2(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
                ztor(ji) = ztor(ji)*RADIAN
                exit
              end if
            end do
          end if
        else
          if ((iz(1,ji).eq.iz(1,dihk)).AND.(iz(2,ji).eq.(iz(2,dihk))).AND.(iz(3,dihk).gt.0)) then
            iz(3,ji) = dihk
            iz(4,ji) = 0
!           basically an improper (more stable)
            call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
            ztor(ji) = ztor(ji)*RADIAN
          end if
          if (iz(3,ji).eq.0) then
            do k=1,j-1
              if (invpdbmap(k).le.0) cycle
              ki = invpdbmap(k)
              if (((iz(1,ki).eq.iz(1,ji)).OR.(iz(1,iz(1,ji)).eq.ki)).AND.(ki.ne.iz(2,ji))) then
                iz(3,ji) = ki
                call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
                if ((iz(1,ki).eq.iz(1,ji))) then
                  iz(4,ji) = 0
                else            
                  iz(4,ji) = -1
                  if (ztor(ji).lt.0.0) iz(4,ji) = 1
                  call bondang2(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
                end if
                ztor(ji) = ztor(ji)*RADIAN
                exit
              end if
            end do
          end if
          if (iz(3,ji).eq.0) then
            if ((iz(1,iz(2,ji)).gt.0).AND.(iz(1,iz(2,ji)).ne.iz(1,ji))) then
              iz(3,ji) = iz(1,iz(2,ji))
              call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
              ztor(ji) = ztor(ji)*RADIAN
              iz(4,ji) = 0
            else
              do k=1,j-1
                if (invpdbmap(k).le.0) cycle
                ki = invpdbmap(k)
                if ((iz(1,ki).eq.iz(2,ji)).AND.(ki.ne.iz(1,ji))) exit
              end do
              iz(3,ji) = ki
              call dihed(xyzdum(:,pdbmap(iz(3,ji))),xyzdum(:,pdbmap(iz(2,ji))),xyzdum(:,pdbmap(iz(1,ji))),xyzdum(:,j),ztor(ji))
              ztor(ji) = ztor(ji)*RADIAN
              iz(4,ji) = 0
            end if
          end if
        end if
      end if
    end if
!    write(*,*) iz(1,ji),blen(ji),iz(2,ji),bang(ji),iz(3,ji),ztor(ji)    
  end do
!
  if (knty.gt.0) then
    do j=unkbnd(i,1),unkbnd(i,2)
      rcnt = at(unkbnd(i,4))%bb(1)+j-unkbnd(i,1)
      atnam(invpdbmap(j))(3:3) = atnam(rcnt)(3:3)
    end do
  else
    do j=1,n_biotyp
      if (btyh(j).le.1) then
        btyh(j) = -100
      else
        btyh(j) = 0
      end if
    end do
    k = 1
    do j=unkbnd(i,1),unkbnd(i,2)
      ji = invpdbmap(j)
      if (btyh(b_type(ji)).ge.0) then
        btyh(b_type(ji)) = btyh(b_type(ji)) + 1
        call int2str(btyh(b_type(ji)),atnam(ji)(3:3),k)
      end if
    end do
  end if
  deallocate(btyh)
!
! the linkage to the next residue also needs to be fixed;
! only do this if next residue is natively supported and not in a new molecule
  if (rs.lt.rsmol(molofrs(rs),2)) then
    if (seqtyp(rs+1).gt.0) then
!     these must be processed in their native order
      do j=at(rs+1)%bb(1),at(rs+1)%bb(1)+at(rs+1)%nbb+at(rs+1)%nsc-1
        if (((iz(1,j).le.0).OR.(iz(2,j).le.0).OR.(iz(3,j).le.0)).AND.(pdbmap(j).le.0)) then
          call fexit()
        end if
        if ((iz(1,j).gt.0).AND.(iz(2,j).gt.0).AND.(iz(3,j).gt.0)) cycle
!       bond lengths and angles across the linkage are handled pretty rigidly again
        if (iz(1,j).eq.0) then
          d2 = HUGE(d2)
          do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
            if (invpdbmap(k).le.0) cycle
            if (invpdbmap(k).ge.j) cycle
            call dis_bound2(xyzdum(:,pdbmap(j)),xyzdum(:,k),d2t)
            if (d2t.lt.d2) then 
              d2 = d2t
              ii = k
            end if
          end do
          iz(1,j) = invpdbmap(ii)
          blen(j) = sqrt(d2)
        end if
        if (iz(2,j).eq.0) then
          if (iz(1,iz(1,j)).gt.0) then
            iz(2,j) = iz(1,iz(1,j))
            call bondang2(xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),bang(j))
            bang(j) = bang(j)*RADIAN
          else
            do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
              if (invpdbmap(k).le.0) cycle
              if (invpdbmap(k).ge.j) cycle
              ki = invpdbmap(k)
              if (iz(1,ki).eq.iz(1,j)) then
                iz(2,j) = ki
                call bondang2(xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),bang(j))
                bang(j) = bang(j)*RADIAN
                exit
              end if
            end do
          end if
        end if
        if (iz(3,j).eq.0) then
          dihthere = .false.
          do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
            if (invpdbmap(k).le.0) cycle
            if (invpdbmap(k).ge.j) cycle
            ki = invpdbmap(k)
            if (((iz(1,ki).eq.iz(2,j)).AND.(iz(2,ki).eq.iz(1,j))).OR.&
 &              ((iz(1,ki).eq.iz(1,j)).AND.(iz(2,ki).eq.iz(2,j)))) then
              dihthere = .true.
              dihk = ki
              exit
            end if
          end do
          if (dihthere.EQV..false.) then
            if ((iz(1,iz(2,j)).gt.0).AND.(iz(1,iz(2,j)).ne.iz(1,j))) then
              iz(3,j) = iz(1,iz(2,j))
              call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                       xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
              ztor(j) = ztor(j)*RADIAN
              iz(4,j) = 0
            else
              do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
                if (invpdbmap(k).le.0) cycle
                if (invpdbmap(k).ge.j) cycle
                ki = invpdbmap(k)
                if ((iz(1,ki).eq.iz(2,j)).AND.(ki.ne.iz(1,j))) then
                  iz(3,j) = ki
                  call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                           xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                  ztor(j) = ztor(j)*RADIAN
                  iz(4,j) = 0
                  exit
                end if
              end do
            end if
            if (iz(3,j).eq.0) then
              do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
                if (invpdbmap(k).le.0) cycle
                if (invpdbmap(k).ge.j) cycle
                ki = invpdbmap(k)
                if (((iz(1,ki).eq.iz(1,j)).OR.(iz(1,iz(1,j)).eq.ki)).AND.(ki.ne.iz(2,j))) then
                  iz(3,j) = ki
                  call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                           xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                  iz(4,j) = -1
                  if (ztor(j).lt.0.0) iz(4,j) = 1
                  call bondang2(xyzdum(:,k),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                  ztor(j) = ztor(j)*RADIAN
                  exit
                end if
              end do
            end if
          else
            if ((iz(1,j).eq.iz(1,dihk)).AND.(iz(2,j).eq.(iz(2,dihk))).AND.(iz(3,dihk).gt.0)) then
              iz(3,j) = dihk
              iz(4,j) = 0
!             basically an improper (more stable)
              call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),xyzdum(:,pdbmap(iz(1,j))),&
 &                      xyzdum(:,pdbmap(j)),ztor(j))
              ztor(j) = ztor(j)*RADIAN
            end if
            do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
              if (invpdbmap(k).le.0) cycle
              if (invpdbmap(k).ge.j) cycle
              ki = invpdbmap(k)
              if (((iz(1,ki).eq.iz(1,j)).OR.(iz(1,iz(1,j)).eq.ki)).AND.(ki.ne.iz(2,j))) then
                iz(3,j) = ki
                call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                         xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                if ((iz(1,ki).eq.iz(1,j))) then
                  iz(4,i) = 0
                else
                  ztor(j) = -1
                  if (ztor(j).lt.0.0) iz(4,j) = 1
                  call bondang2(xyzdum(:,k),xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                end if
                ztor(j) = ztor(j)*RADIAN
                exit
              end if
            end do
            if (iz(3,j).eq.0) then
              if ((iz(1,iz(2,j)).gt.0).AND.(iz(1,iz(2,j)).ne.iz(1,j))) then
                iz(3,j) = iz(1,iz(2,j))
                call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                         xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                ztor(j) = ztor(j)*RADIAN
                iz(4,j) = 0
              else
                do k=unkbnd(i,1),min(nl,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc)
                  if (invpdbmap(k).le.0) cycle
                  if (invpdbmap(k).ge.j) cycle
                  ki = invpdbmap(k)
                  if ((iz(1,ki).eq.iz(2,j)).AND.(ki.ne.iz(1,j))) then
                    iz(3,j) = ki
                    call dihed(xyzdum(:,pdbmap(iz(3,j))),xyzdum(:,pdbmap(iz(2,j))),&
 &                             xyzdum(:,pdbmap(iz(1,j))),xyzdum(:,pdbmap(j)),ztor(j))
                    ztor(j) = ztor(j)*RADIAN
                    iz(4,j) = 0
                    exit
                  end if
                end do
              end if
            end if
          end if
        end if
      end do
!     correct for backwards topology reference, but watch out for proline sidechain atoms, which must be
!     defined independently of the previous residue
      do j=unkbnd(i,2)+1,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc
        if (invpdbmap(j).le.0) cycle
        ji = invpdbmap(j)
        if (iz(4,ji).eq.0) then
          if ((pdbmap(iz(1,ji)).gt.unkbnd(i,2)).AND.(pdbmap(iz(2,ji)).gt.unkbnd(i,2)).AND.(pdbmap(iz(3,ji)).le.unkbnd(i,2))) then
            do k=j+1,unkbnd(i,2)+at(rs+1)%nbb+at(rs+1)%nsc
              if (invpdbmap(k).le.0) cycle
              ki = invpdbmap(k)
              if ((seqtyp(rs+1).eq.9).OR.(seqtyp(rs+1).eq.25).OR.(seqtyp(rs+1).eq.32)) then
                if (ki.ge.at(rs+1)%sc(2-shf)) cycle
              end if
              if (((iz(1,ki).eq.iz(2,ji)).AND.(iz(2,ki).eq.iz(1,ji))).OR.&
 &                ((iz(1,ki).eq.iz(1,ji)).AND.(iz(2,ki).eq.iz(2,ji)).AND.(iz(4,ki).eq.0).AND.(iz(3,ki).ne.ji))) then
                if (iz(1,ki).eq.iz(2,ji)) then
                  iz(3,ki) = iz(3,ji)
                else if (iz(1,ki).eq.iz(1,ji)) then
                  iz(3,ki) = ji      
                end if
                call dihed(xyzdum(:,pdbmap(iz(3,ki))),xyzdum(:,pdbmap(iz(2,ki))),&
 &                         xyzdum(:,pdbmap(iz(1,ki))),xyzdum(:,k),ztor(ki))
                iz(4,ki) = 0
                ztor(ki) = ztor(ki)*RADIAN         
                exit
              end if
            end do
            exit
          end if
        end if
      end do
    end if
  end if
!
  do j=unkbnd(i,1),unkbnd(i,2)
    ji = invpdbmap(j)
    do k=j+1,unkbnd(i,2)
      ki = invpdbmap(k)
      call dis_bound2(xyzdum(:,j),xyzdum(:,k),d2t)
      if ((iz(1,ji).ne.ki).AND.(iz(1,ki).ne.ji).AND.(d2t.lt.d2max(j)*d2max(k))) then
        nadd = nadd + 1
        iadd(1,nadd) = ji
        iadd(2,nadd) = ki
      end if
    end do
  end do
!
  deallocate(invpdbmap)
!
end
!
!-----------------------------------------------------------------------------
!
! a subroutine to correct linkages if an unknown residue was identified as a supported polymer type
! also appends eligible residue lists for move types
!
subroutine correct_topology(rs,nl,xyzdum)
!
  use pdb
  use zmatrix
  use sequen
  use molecule
  use math
  use polypep
  use fyoc
  use aminos
  use params
  use movesets
  use atoms
  use system
  use iounit
!
  implicit none
!
  integer imol,shf,shf2,k,cbi,cgi,cdi,c2i,c1i,o4i,c3i,c4i
  integer, INTENT(IN) :: nl,rs
  RTYPE xyzdum(3,nl),blenp
  logical pckyes,renter
!
  pckyes = .false.
  imol = molofrs(rs)
!
  if (seqpolty(rs).eq.'P') then
    if (rs.gt.rsmol(imol,1)) then
      if ((ci(rs-1).gt.0).AND.(cai(rs-1).gt.0)) then
        iz(2,cai(rs)) = ci(rs-1)
        iz(3,cai(rs)) = cai(rs-1)
        wline(rs) = cai(rs)
        call bondang2(xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(ci(rs-1))),bang(cai(rs)))
        bang(cai(rs)) = RADIAN*bang(cai(rs))
        call dihed(xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(ci(rs-1))),&
 &                                                                        xyzdum(:,pdbmap(cai(rs-1))),ztor(cai(rs)))
        ztor(cai(rs)) = RADIAN*ztor(cai(rs))
      end if
      if (((seqtyp(rs-1).eq.26).AND.(seqpolty(rs-1).eq.'P')).OR.(aminopolty(seqtyp(rs-1)).eq.'P')) then
        if (yline(rs-1).gt.0) then
          if (iz(1,yline(rs-1)).eq.iz(1,ni(rs))) then
            yline2(rs-1) = ni(rs)
          end if
        end if
      end if
      if (ci(rs-1).gt.0) then ! should be redundant check
        iz(2,ci(rs)) = ni(rs)
        iz(3,ci(rs)) = ci(rs-1)
        fline(rs) = ci(rs)
        call bondang2(xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),bang(ci(rs)))
        bang(ci(rs)) = RADIAN*bang(ci(rs))
        call dihed(xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),&
 &                                                                        xyzdum(:,pdbmap(ci(rs-1))),ztor(ci(rs)))
        ztor(ci(rs)) = RADIAN*ztor(ci(rs))
      end if
      if (oi(rs).gt.0) then ! should also be redundant
        iz(2,oi(rs)) = cai(rs)
        iz(3,oi(rs)) = ni(rs)
        yline(rs) = oi(rs)
        call bondang2(xyzdum(:,pdbmap(oi(rs))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),bang(oi(rs)))
        bang(oi(rs)) = RADIAN*bang(oi(rs))
        call dihed(xyzdum(:,pdbmap(oi(rs))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                        xyzdum(:,pdbmap(ni(rs))),ztor(oi(rs)))
        ztor(oi(rs)) = RADIAN*ztor(oi(rs))
      end if
      if (hni(rs).gt.0) then
        iz(1,hni(rs)) = ni(rs)
        iz(2,hni(rs)) = cai(rs)
        iz(3,hni(rs)) = ci(rs)
        if ((fline(rs).ne.hni(rs)).AND.(fline(rs).gt.0)) fline2(rs) = hni(rs)
        call dis_bound2(xyzdum(:,pdbmap(hni(rs))),xyzdum(:,pdbmap(ni(rs))),blen(hni(rs)))
        blen(hni(rs)) = sqrt(blen(hni(rs)))
        call bondang2(xyzdum(:,pdbmap(hni(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),bang(hni(rs)))
        bang(hni(rs)) = RADIAN*bang(hni(rs))
        call dihed(xyzdum(:,pdbmap(hni(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                        xyzdum(:,pdbmap(ci(rs))),ztor(hni(rs)))
        ztor(hni(rs)) = RADIAN*ztor(hni(rs))
      end if
    end if
    shf = 0
    if (ua_model.ne.0) shf = 1
!   detect if there is a pucker to consider for sampling (exclude N-terminal residues for this)
    if ((at(rs)%nsc.ge.(4-shf)).AND.(rs.gt.rsmol(imol,1))) then
      if ((iz(1,at(rs)%sc(3-shf)).eq.at(rs)%sc(2-shf)).AND.(iz(1,at(rs)%sc(2-shf)).eq.cai(rs)).AND.&
 &        ((iz(1,at(rs)%sc(4-shf)).eq.at(rs)%sc(3-shf)).OR.(iz(1,at(rs)%sc(4-shf)).eq.ni(rs))).AND.(ci(rs-1).gt.0)) then
        call dis_bound2(xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(at(rs)%sc(4-shf))),blenp)
        if ((bio_code(b_type(at(rs)%sc(2-shf)))(1:2).eq.'CB').AND.(bio_code(b_type(at(rs)%sc(3-shf)))(1:2).eq.'CG').AND.&
 &          (bio_code(b_type(at(rs)%sc(4-shf)))(1:2).eq.'CD').AND.(blenp.le.1.7)) then
!         set Z-matrix to necessary format
          cbi = at(rs)%sc(2-shf)
          cgi = at(rs)%sc(3-shf)
          cdi = at(rs)%sc(4-shf)
          iz(1,cbi) = cai(rs)
          iz(2,cbi) = ni(rs)
          iz(3,cbi) = ci(rs-1)
          iz(4,cbi) = 0
          iz(1,cgi) = cbi
          iz(2,cgi) = cai(rs)
          iz(3,cgi) = ni(rs)
          iz(4,cgi) = 0
          iz(1,cdi) = ni(rs)
          iz(2,cdi) = cai(rs)
          iz(3,cdi) = cbi
          iz(4,cdi) = 0
          call dis_bound2(xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),blen(cbi))
          blen(cbi) = sqrt(blen(cbi))
          call dis_bound2(xyzdum(:,pdbmap(cgi)),xyzdum(:,pdbmap(cbi)),blen(cgi))
          blen(cgi) = sqrt(blen(cgi))
          call dis_bound2(xyzdum(:,pdbmap(cdi)),xyzdum(:,pdbmap(ni(rs))),blen(cdi))
          blen(cdi) = sqrt(blen(cdi))
          call bondang2(xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),bang(cbi))
          bang(cbi) = RADIAN*bang(cbi)
          call bondang2(xyzdum(:,pdbmap(cgi)),xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),bang(cgi))
          bang(cgi) = RADIAN*bang(cgi)
          call bondang2(xyzdum(:,pdbmap(cdi)),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),bang(cdi))
          bang(cdi) = RADIAN*bang(cdi)
          call dihed(xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(ci(rs-1))),ztor(cbi))
          ztor(cbi) = RADIAN*ztor(cbi)
          call dihed(xyzdum(:,pdbmap(cgi)),xyzdum(:,pdbmap(cbi)),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(ni(rs))),ztor(cgi))
          ztor(cgi) = RADIAN*ztor(cgi)
          call dihed(xyzdum(:,pdbmap(cdi)),xyzdum(:,pdbmap(ni(rs))),xyzdum(:,pdbmap(cai(rs))),xyzdum(:,pdbmap(cbi)),ztor(cdi))
          ztor(cdi) = RADIAN*ztor(cdi)
          pucline(rs) = ci(rs) 
          do k=1,nadd
            if ((iadd(1,k).eq.cbi).OR.(iadd(1,k).eq.cdi).OR.(iadd(2,k).eq.cbi).OR.(iadd(2,k).eq.cdi)) then
              iadd(1,k) = cgi
              iadd(2,k) = cdi
              exit
            end if
          end do
!         improve setup for atoms directly attached (primarily to make auxiliary bond angle variations consistent between systems)
          do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nsc+at(rs)%nbb
!           find atoms directly attached to the ring with a dihedral Z-matrix setup
            if (((iz(1,k).eq.cbi).OR.(iz(1,k).eq.cgi).OR.(iz(1,k).eq.cdi).OR.&
 &                (iz(1,k).eq.cai(rs)).OR.(iz(1,k).eq.ni(rs))).AND.(iz(4,k).eq.0).AND.(k.ne.ni(rs)).AND.(k.ne.cai(rs)).AND.&
 &               (k.ne.cbi).AND.(k.ne.cgi).AND.(k.ne.cdi).AND.(k.ne.ci(rs))) then
              if ((iz(1,k).eq.ni(rs)).AND.(iz(2,k).eq.cai(rs)).AND.(iz(3,k).ne.cdi).AND.(k.gt.cdi)) iz(3,k) = cdi
              if ((iz(1,k).eq.ni(rs)).AND.(iz(3,k).ne.cai(rs)).AND.(iz(2,k).eq.cdi).AND.(k.gt.cai(rs))) iz(3,k) = cai(rs)
              if ((iz(1,k).eq.cai(rs)).AND.(iz(2,k).eq.ni(rs)).AND.(iz(3,k).ne.cbi).AND.(k.gt.cbi)) iz(3,k) = cbi
              if ((iz(1,k).eq.cai(rs)).AND.(iz(3,k).ne.ni(rs)).AND.(iz(2,k).eq.cbi).AND.(k.gt.ni(rs))) iz(3,k) = ni(rs)
              if ((iz(1,k).eq.cbi).AND.(iz(2,k).eq.cai(rs)).AND.(iz(3,k).ne.cgi).AND.(k.gt.cgi)) iz(3,k) = cgi
              if ((iz(1,k).eq.cbi).AND.(iz(3,k).ne.cai(rs)).AND.(iz(2,k).eq.cgi).AND.(k.gt.cai(rs))) iz(3,k) = cai(rs)
              if ((iz(1,k).eq.cgi).AND.(iz(2,k).eq.cbi).AND.(iz(3,k).ne.cdi).AND.(k.gt.cdi)) iz(3,k) = cdi
              if ((iz(1,k).eq.cgi).AND.(iz(3,k).ne.cbi).AND.(iz(2,k).eq.cdi).AND.(k.gt.cbi)) iz(3,k) = cbi
              if ((iz(1,k).eq.cdi).AND.(iz(2,k).eq.cgi).AND.(iz(3,k).ne.ni(rs)).AND.(k.gt.ni(rs))) iz(3,k) = ni(rs)
              if ((iz(1,k).eq.cdi).AND.(iz(3,k).ne.cgi).AND.(iz(2,k).eq.ni(rs)).AND.(k.gt.cgi)) iz(3,k) = cgi
              call dihed(xyzdum(:,pdbmap(k)),xyzdum(:,pdbmap(iz(1,k))),xyzdum(:,pdbmap(iz(2,k))),xyzdum(:,pdbmap(iz(3,k))),&
 &                       ztor(k))
              ztor(k) = RADIAN*ztor(k)
            end if
          end do
        end if
      end if
    end if
    if (rs.lt.rsmol(imol,2)) then
      if (((seqtyp(rs+1).eq.26).AND.(seqpolty(rs+1).eq.'P')).OR.(aminopolty(seqtyp(rs+1)).eq.'P')) then
        iz(1,ni(rs+1)) = ci(rs)
        iz(2,ni(rs+1)) = cai(rs)
        iz(3,ni(rs+1)) = ni(rs)
        iz(4,ni(rs+1)) = 0
        iz(2,cai(rs+1)) = ci(rs)
        iz(3,cai(rs+1)) = cai(rs)
        iz(4,cai(rs+1)) = 0
        iz(3,ci(rs+1)) = ci(rs)
        iz(4,ci(rs+1)) = 0
        call dis_bound2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),blen(ni(rs+1)))
        blen(ni(rs+1)) = sqrt(blen(ni(rs+1)))
        call bondang2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),bang(ni(rs+1)))
        bang(ni(rs+1)) = RADIAN*bang(ni(rs+1))
        call bondang2(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),bang(cai(rs+1)))
        bang(cai(rs+1)) = RADIAN*bang(cai(rs+1))
        call dihed(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                     xyzdum(:,pdbmap(ni(rs))),ztor(ni(rs+1)))
        ztor(ni(rs+1)) = RADIAN*ztor(ni(rs+1))
        if (yline(rs).gt.0) yline2(rs) = ni(rs+1)
        call dihed(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),&
 &                                                                        xyzdum(:,pdbmap(cai(rs))),ztor(cai(rs+1)))
        ztor(cai(rs+1)) = RADIAN*ztor(cai(rs+1))
        wline(rs+1) = cai(rs+1)
        call dihed(xyzdum(:,pdbmap(ci(rs+1))),xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),&
 &                                                                        xyzdum(:,pdbmap(ci(rs))),ztor(ci(rs+1)))
        ztor(ci(rs+1)) = RADIAN*ztor(ci(rs+1))
        fline(rs+1) = ci(rs+1)
        if (hni(rs+1).gt.0) then
          iz(1,hni(rs+1)) = ni(rs+1)
          iz(2,hni(rs+1)) = cai(rs+1)
          iz(3,hni(rs+1)) = ci(rs+1)
          iz(4,hni(rs+1)) = 0
          fline2(rs+1) = hni(rs+1)
          call dihed(xyzdum(:,pdbmap(hni(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(cai(rs+1))),&
 &                                                                          xyzdum(:,pdbmap(ci(rs+1))),ztor(hni(rs+1)))
          ztor(hni(rs+1)) = RADIAN*ztor(hni(rs+1))
        end if
      else if ((cai(rs+1).le.0).OR.(ni(rs+1).le.0)) then
!       do nothing
      else if (iz(1,cai(rs+1)).eq.ni(rs+1)) then
        iz(1,ni(rs+1)) = ci(rs)
        iz(2,ni(rs+1)) = cai(rs)
        iz(3,ni(rs+1)) = ni(rs)
        iz(4,ni(rs+1)) = 0
        iz(2,cai(rs+1)) = ci(rs)
        iz(3,cai(rs+1)) = cai(rs)
        iz(4,cai(rs+1)) = 0
        call dis_bound2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),blen(ni(rs+1)))
        blen(ni(rs+1)) = sqrt(blen(ni(rs+1)))
        call bondang2(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),bang(ni(rs+1)))
        bang(ni(rs+1)) = RADIAN*bang(ni(rs+1))
        call bondang2(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),bang(cai(rs+1)))
        bang(cai(rs+1)) = RADIAN*bang(cai(rs+1))
        call dihed(xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),xyzdum(:,pdbmap(cai(rs))),&
 &                                                                     xyzdum(:,pdbmap(ni(rs))),ztor(ni(rs+1)))
        ztor(ni(rs+1)) = RADIAN*ztor(ni(rs+1))
        if (yline(rs).gt.0) yline2(rs) = ni(rs+1)
        call dihed(xyzdum(:,pdbmap(cai(rs+1))),xyzdum(:,pdbmap(ni(rs+1))),xyzdum(:,pdbmap(ci(rs))),&
 &                                                                        xyzdum(:,pdbmap(cai(rs))),ztor(cai(rs+1)))
        ztor(cai(rs+1)) = RADIAN*ztor(cai(rs+1))
        wline(rs+1) = cai(rs+1)
      end if
    end if
  else if (seqpolty(rs).eq.'N') then
    shf = 0
    if ((nuci(rs,5).le.0).AND.(nuci(rs,6).le.0).AND.(nuci(rs,4).gt.0)) shf = 2
!   make sure own linkage to supported residue prior in sequence is ok (note that d.o.f setup is not handled)
    if (rs.gt.rsmol(imol,1)) then
      if ((seqpolty(rs-1).eq.'N').OR.(aminopolty(seqtyp(rs-1)).eq.'N')) nucsline(6,rs-1) = nuci(rs,1)
      if ((seqtyp(rs-1).ne.26).AND.(aminopolty(seqtyp(rs-1)).eq.'N')) then
        shf2 = 0
        if ((seqtyp(rs-1).ge.76).AND.(seqtyp(rs-1).le.87)) shf2 = 2
        iz(2,nuci(rs,1)) = nuci(rs-1,5-shf2)
        iz(3,nuci(rs,1)) = nuci(rs-1,4-shf2)
        iz(4,nuci(rs,1)) = 0
        iz(2,nuci(rs,2)) = nuci(rs-1,6-shf2)
        iz(3,nuci(rs,2)) = nuci(rs-1,5-shf2)
        iz(4,nuci(rs,2)) = 0
        call bondang2(xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),xyzdum(:,pdbmap(nuci(rs-1,5-shf2))),&
 &                    bang(nuci(rs,1)))
        bang(nuci(rs,1)) = RADIAN*bang(nuci(rs,1))
        call dihed(xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),xyzdum(:,pdbmap(nuci(rs-1,5-shf2))),&
 &                                                                     xyzdum(:,pdbmap(nuci(rs-1,4-shf2))),ztor(nuci(rs,1)))
        ztor(nuci(rs,1)) = RADIAN*ztor(nuci(rs,1))
        call bondang2(xyzdum(:,pdbmap(nuci(rs,2))),xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),&
 &                    bang(nuci(rs,2)))
        bang(nuci(rs,2)) = RADIAN*bang(nuci(rs,2))
        call dihed(xyzdum(:,pdbmap(nuci(rs,2))),xyzdum(:,pdbmap(nuci(rs,1))),xyzdum(:,pdbmap(nuci(rs-1,6-shf2))),&
 &                                                                     xyzdum(:,pdbmap(nuci(rs-1,5-shf2))),ztor(nuci(rs,2)))
        ztor(nuci(rs,2)) = RADIAN*ztor(nuci(rs,2))
        nucsline(5-shf2,rs-1) = nuci(rs,2)
      end if
    end if
!   detect if the sugar is sampleable
    if (at(rs)%nsc.ge.3) then
      if ((iz(1,at(rs)%sc(2)).eq.at(rs)%sc(1)).AND.(iz(1,at(rs)%sc(1)).eq.nuci(rs,6-shf)).AND.&
 &        ((iz(1,at(rs)%sc(3)).eq.at(rs)%sc(2)).OR.(iz(1,at(rs)%sc(3)).eq.nuci(rs,5-shf)))) then
        call dis_bound2(xyzdum(:,pdbmap(nuci(rs,5-shf))),xyzdum(:,pdbmap(at(rs)%sc(3))),blenp)
        if ((bio_code(b_type(at(rs)%sc(1)))(1:2).eq.'C2').AND.(bio_code(b_type(at(rs)%sc(2)))(1:2).eq.'C1').AND.&
 &          (bio_code(b_type(at(rs)%sc(3)))(1:2).eq.'O4').AND.(blenp.le.1.7)) then
!         set Z-matrix to necessary format
          c2i = at(rs)%sc(1)
          c1i = at(rs)%sc(2)
          o4i = at(rs)%sc(3)
          c4i = nuci(rs,5-shf)
          c3i = nuci(rs,6-shf)
          iz(1,c2i) = c3i
          iz(2,c2i) = c4i
          iz(3,c2i) = nuci(rs,4-shf) ! C5*
          iz(4,c2i) = 0
          iz(1,c1i) = c2i
          iz(2,c1i) = c3i
          iz(3,c1i) = c4i
          iz(4,c1i) = 0
          iz(1,o4i) = c4i
          iz(2,o4i) = c3i
          iz(3,o4i) = c2i
          iz(4,o4i) = 0
          call dis_bound2(xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),blen(c2i))
          blen(c2i) = sqrt(blen(c2i))
          call dis_bound2(xyzdum(:,pdbmap(c1i)),xyzdum(:,pdbmap(c2i)),blen(c1i))
          blen(c1i) = sqrt(blen(c1i))
          call dis_bound2(xyzdum(:,pdbmap(o4i)),xyzdum(:,pdbmap(c4i)),blen(o4i))
          blen(o4i) = sqrt(blen(o4i))
          call bondang2(xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c4i)),bang(c2i))
          bang(c2i) = RADIAN*bang(c2i)
          call bondang2(xyzdum(:,pdbmap(c1i)),xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),bang(c1i))
          bang(c1i) = RADIAN*bang(c1i)
          call bondang2(xyzdum(:,pdbmap(o4i)),xyzdum(:,pdbmap(c4i)),xyzdum(:,pdbmap(c3i)),bang(o4i))
          bang(o4i) = RADIAN*bang(o4i)
          call dihed(xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c4i)),xyzdum(:,pdbmap(nuci(rs,4-shf))),ztor(c2i))
          ztor(c2i) = RADIAN*ztor(c2i)
          call dihed(xyzdum(:,pdbmap(c1i)),xyzdum(:,pdbmap(c2i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c4i)),ztor(c1i))
          ztor(c1i) = RADIAN*ztor(c1i)
          call dihed(xyzdum(:,pdbmap(o4i)),xyzdum(:,pdbmap(c4i)),xyzdum(:,pdbmap(c3i)),xyzdum(:,pdbmap(c2i)),ztor(o4i))
          ztor(o4i) = RADIAN*ztor(o4i)
          pckyes = .true.
          do k=1,nadd
            if ((iadd(1,k).eq.c2i).OR.(iadd(1,k).eq.o4i).OR.(iadd(2,k).eq.c2i).OR.(iadd(2,k).eq.o4i)) then
              iadd(1,k) = c1i
              iadd(2,k) = o4i
              exit
            end if
          end do
          if (rs.eq.rsmol(imol,2)) then
            renter = .false.
            do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1
              if ((iz(1,k).eq.c3i.AND.(iz(2,k).eq.c4i)).AND.(iz(4,k).eq.0).AND.(iz(1,iz(3,k)).ne.c3i)) then
                nucsline(6,rs) = k
                renter = .true.
                exit
              end if
            end do
            if (renter.EQV..false.) pckyes = .false.
          end if
          if (pckyes.EQV..true.) then
!           improve setup for atoms directly attached (primarily to make auxiliary bond angle variations consistent between systems)
            do k=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nsc+at(rs)%nbb
!             find atoms directly attached to the ring with a dihedral Z-matrix setup
              if (((iz(1,k).eq.c2i).OR.(iz(1,k).eq.c1i).OR.(iz(1,k).eq.o4i).OR.&
 &                  (iz(1,k).eq.c3i).OR.(iz(1,k).eq.c4i)).AND.(iz(4,k).eq.0).AND.(k.ne.c4i).AND.(k.ne.c3i).AND.&
 &                 (k.ne.c2i).AND.(k.ne.c1i).AND.(k.ne.o4i).AND.(k.ne.nuci(rs,4-shf))) then
                if ((iz(1,k).eq.c4i).AND.(iz(2,k).eq.c3i).AND.(iz(3,k).ne.o4i).AND.(k.gt.o4i)) iz(3,k) = o4i
                if ((iz(1,k).eq.c4i).AND.(iz(3,k).ne.c3i).AND.(iz(2,k).eq.o4i).AND.(k.gt.c3i)) iz(3,k) = c3i
                if ((iz(1,k).eq.c3i).AND.(iz(2,k).eq.c4i).AND.(iz(3,k).ne.c2i).AND.(k.gt.c2i)) iz(3,k) = c2i
                if ((iz(1,k).eq.c3i).AND.(iz(3,k).ne.c4i).AND.(iz(2,k).eq.c2i).AND.(k.gt.c4i)) iz(3,k) = c4i
                if ((iz(1,k).eq.c2i).AND.(iz(2,k).eq.c3i).AND.(iz(3,k).ne.c1i).AND.(k.gt.c1i)) iz(3,k) = c1i
                if ((iz(1,k).eq.c2i).AND.(iz(3,k).ne.c3i).AND.(iz(2,k).eq.c1i).AND.(k.gt.c3i)) iz(3,k) = c3i
                if ((iz(1,k).eq.c1i).AND.(iz(2,k).eq.c2i).AND.(iz(3,k).ne.o4i).AND.(k.gt.o4i)) iz(3,k) = o4i
                if ((iz(1,k).eq.c1i).AND.(iz(3,k).ne.c2i).AND.(iz(2,k).eq.o4i).AND.(k.gt.c2i)) iz(3,k) = c2i
                if ((iz(1,k).eq.o4i).AND.(iz(2,k).eq.c1i).AND.(iz(3,k).ne.c4i).AND.(k.gt.c4i)) iz(3,k) = c4i
                if ((iz(1,k).eq.o4i).AND.(iz(3,k).ne.c1i).AND.(iz(2,k).eq.c4i).AND.(k.gt.c1i)) iz(3,k) = c1i
                call dihed(xyzdum(:,pdbmap(k)),xyzdum(:,pdbmap(iz(1,k))),xyzdum(:,pdbmap(iz(2,k))),xyzdum(:,pdbmap(iz(3,k))),&
 &                         ztor(k))
                ztor(k) = RADIAN*ztor(k)
              end if
            end do
          end if
        end if
      end if
    end if
    if (rs.lt.rsmol(imol,2)) then
      if ((nuci(rs+1,1).gt.0).AND.(nuci(rs+1,2).gt.0)) then
        iz(1,nuci(rs+1,1)) = nuci(rs,6-shf)
        iz(2,nuci(rs+1,1)) = nuci(rs,5-shf)
        iz(3,nuci(rs+1,1)) = nuci(rs,4-shf)
        iz(4,nuci(rs+1,1)) = 0
        iz(2,nuci(rs+1,2)) = nuci(rs,6-shf)
        iz(3,nuci(rs+1,2)) = nuci(rs,5-shf)
        iz(4,nuci(rs+1,2)) = 0
        if (nuci(rs+1,3).gt.0) then
          iz(3,nuci(rs+1,3)) = nuci(rs,6-shf)
          iz(4,nuci(rs+1,3)) = 0
        end if
        call dis_bound2(xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),blen(nuci(rs+1,1)))
        blen(nuci(rs+1,1)) = sqrt(blen(nuci(rs+1,1)))
        call bondang2(xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),xyzdum(:,pdbmap(nuci(rs,5-shf))),&
 &                    bang(nuci(rs+1,1)))
        bang(nuci(rs+1,1)) = RADIAN*bang(nuci(rs+1,1))
        call bondang2(xyzdum(:,pdbmap(nuci(rs+1,2))),xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),&
 &                    bang(nuci(rs+1,2)))
        bang(nuci(rs+1,2)) = RADIAN*bang(nuci(rs+1,2))
        call dihed(xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),xyzdum(:,pdbmap(nuci(rs,5-shf))),&
 &                                                                     xyzdum(:,pdbmap(nuci(rs,4-shf))),ztor(nuci(rs+1,1)))
        ztor(nuci(rs+1,1)) = RADIAN*ztor(nuci(rs+1,1))
        call dihed(xyzdum(:,pdbmap(nuci(rs+1,2))),xyzdum(:,pdbmap(nuci(rs+1,1))),xyzdum(:,pdbmap(nuci(rs,6-shf))),&
 &                                                                        xyzdum(:,pdbmap(nuci(rs,5-shf))),ztor(nuci(rs+1,2)))
        ztor(nuci(rs+1,2)) = RADIAN*ztor(nuci(rs+1,2))
        nnucs(rs) = nnucs(rs) + 1
        nucsline(nnucs(rs),rs) = nuci(rs+1,2)
        if (nuci(rs+1,3).gt.0) then
          call dihed(xyzdum(:,pdbmap(nuci(rs+1,3))),xyzdum(:,pdbmap(nuci(rs+1,2))),xyzdum(:,pdbmap(nuci(rs+1,1))),&
 &                                                                        xyzdum(:,pdbmap(nuci(rs,6-shf))),ztor(nuci(rs+1,3)))
          ztor(nuci(rs+1,3)) = RADIAN*ztor(nuci(rs+1,3))
          if (nnucs(rs+1).gt.0) nucsline(1,rs+1) = nuci(rs+1,3)
        end if
      end if
    end if
  else if (rs.lt.rsmol(imol,2)) then
    if ((seqtyp(rs+1).ne.26).AND.(aminopolty(seqtyp(rs+1)).eq.'P')) then
      if ((ci(rs).le.0).OR.(iz(1,at(rs+1)%bb(1)).ne.ci(rs))) then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the phi angle in polypeptide residue ',rs+1,' may be &
 &altered.'
      end if
      if ((ci(rs).le.0).OR.(iz(1,at(rs+1)%bb(1)).ne.ci(rs)).OR.&
 &        (oi(rs).le.0)) then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the omega angle in polypeptide residue ',rs+1,' may be &
 &altered.'
      else if (iz(1,oi(rs)).ne.ci(rs)) then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the omega angle in polypeptide residue ',rs+1,' may be &
 &altered.'
      end if
    else if ((seqtyp(rs+1).ne.26).AND.(aminopolty(seqtyp(rs+1)).eq.'N')) then
      if (bio_code(b_type(iz(1,at(rs+1)%bb(1))))(1:2).ne.'C3') then
        write(ilog,*) 'Warning. Due to unsupported residue of unsupported polymer type being &
 &immediately N-terminal, the meaning of the first two nucleic acid backbone angles in polynucleotide &
 &residue ',rs+1,' may be altered.'
      end if
    end if
  end if
!
end
!
!-------------------------------------------------------------------------------------
!
! this routine assembles a list of small groups in simple data structures
! generally, these will simply lump hydrogens
!
subroutine setup_topgrps()
!
  use atoms
  use zmatrix
  use mcgrid
  use cutoffs
!
  implicit none
!
  integer i
!
  allocate(which_topgrp(n))
  which_topgrp(:) = 0
!
  n_topgrps = 0
  do i=1,n
    if (n12(i).eq.0) then
      n_topgrps = n_topgrps + 1
      which_topgrp(i) = n_topgrps
    else if ((n12(i).eq.1).AND.(blen(i).lt.topgrprad)) then
      if (which_topgrp(iz(1,i)).gt.0) then
        which_topgrp(i) = which_topgrp(iz(1,i))
      else
        n_topgrps = n_topgrps + 1
        which_topgrp(i) = n_topgrps
      end if
    else
      n_topgrps = n_topgrps + 1
      which_topgrp(i) = n_topgrps
    end if
  end do
!
  allocate(topgrplst(n_topgrps))
  topgrplst(:) = 0
  allocate(topgrpxyz(3,n_topgrps))
  if (use_mcgrid.EQV..true.) then
    allocate(grid%tgrpgp(n_topgrps))
  end if
  do i=1,n
    if (topgrplst(which_topgrp(i)).gt.0) then
      if ((n12(i).eq.1).AND.(blen(i).lt.topgrprad)) then
!       do nothing        
      else
        topgrplst(which_topgrp(i)) = i ! may overwrite for eligible X-X molecule
      end if
    else
      topgrplst(which_topgrp(i)) = i
    end if
  end do
!
end
!
!-----------------------------------------------------------------------------
!
! this routine generates the position of atom at based on the Cartesian
! coordinates of the reference atoms (lower index number) and the passed
! internals
! obviously has to be called in order (atom-wise)
! chirality indicates the chirality with respect to the ordering of ref.
! atoms (1 or -1), if set to 0 two bond angles are used instead (generally
! discouraged due to numerical in"stability")
! arguments are passed in distance units and degrees, respectively
! for the first few atoms, not all reference atoms might be provided
! -> build in a standard frame
!
subroutine genxyz(at,i2,bl,i3,ba,i4,baodi,chirality)
!
  use atoms
  use iounit
  use math
  use zmatrix
!
  implicit none
!
  integer at,i2,i3,i4,chirality
  RTYPE bl,ba,bl2,baodi,dv3(3),dv2(3),stm1,ctm1,dv4(3)
  RTYPE trx,trz,stm2,ctm2,s1,c1,s2,c2,bl3,cp1(3),np1,cp2(3)
!
!
 67   format(a,i6,1x,i6,1x,i6,a,i7,a)
 68   format('Fatal. Ill-defined dihedral angle for atom ',i7,' (1-2:',&
 &i7,': ',f9.4,' A, 1-3:',i7,': ',f9.4,' deg., 1-4:',i7,': ',f9.4,&
 &' deg.)')
 69   format('Warning. Ill-defined bond angles for atom ',i7,' (1-2:',&
 &i7,': ',f9.4,' A, 1-3:',i7,': ',f9.4,' deg., 1-3:',i7,': ',f9.4,&
 &' deg.)')
 70   format('Fatal. Colinear reference atoms for atom ',i7,' (1-2:',&
 &i7,': ',f9.4,' A, 1-3:',i7,': ',f9.4,' deg., 1-3:',i7,': ',f9.4,&
 &' deg.)')
! sanity check
  if ((i2.lt.0).OR.(i3.lt.0).OR.(i4.lt.0)) then
     write(*,*) i2,i3,i4
    write(ilog,67) 'Fatal. Bad atom index in genxyz (',&
 & i2,i3,i4,') for atom ',at,'. This is a bug.'
    call fexit()
  end if
!
! first atom in a molecule always at the origin (note that
! for re-building later these are skipped!), second on +z-axis,
! third in (+)xz-plane
!
  if (i2.eq.0) then
    x(at) = 0.0d0
    y(at) = 0.0d0
    z(at) = 0.0d0
  else if (i3.eq.0) then
    x(at) = x(i2)
    y(at) = y(i2)
    z(at) = z(i2) + bl
  else if (i4.eq.0) then
    stm1 = sin(ba/RADIAN)
    ctm1 = cos(ba/RADIAN)
    dv2(1) = x(i2) - x(i3)
    dv2(2) = y(i2) - y(i3)
    dv2(3) = z(i2) - z(i3)
    bl2 = sqrt(sum(dv2(:)**2))
    dv2(:) = dv2(:)/bl2
!   rel. xz in the assumed frame
    trx = bl*stm1
    trz = bl2 - bl*ctm1
!   shift/rotate into real frame 
    s1 = sqrt(sum(dv2(1:2)**2))
    if (s1.gt.0.0) then
      c1 = dv2(2)/s1
      s2 = dv2(1)/s1
    else
      c1 = 1.0
      s2 = 0.0
    end if
    x(at) = x(i3) + trx*c1 + trz*s2*s1
    y(at) = y(i3) - trx*s2 + trz*c1*s1
    z(at) = z(i3) + trz*dv2(3)
!
! now the general (and most often used) cases
!
  else
!
    if (chirality.eq.0) then
!
!     with a torsional angle
!
      stm1 = ba/RADIAN
      stm2 = baodi/RADIAN
      ctm1 = cos(stm1)
      ctm2 = cos(stm2)
      stm2 = sin(stm2)
      stm1 = sin(stm1)
      dv2(1) = x(i2) - x(i3)
      dv2(2) = y(i2) - y(i3)
      dv2(3) = z(i2) - z(i3)
      bl2 = sqrt(sum(dv2(:)**2))
      dv2(:) = dv2(:)/bl2
      if (stm1.le.0.0) then  ! at,i2,i3 colinear, torsion is irrelevant
        dv4(:) = -dv2(:)*ctm1
        x(at) = x(i2) + bl*dv4(1)
        y(at) = y(i2) + bl*dv4(2)
        z(at) = z(i2) + bl*dv4(3)
        return
      end if
      dv3(1) = x(i3) - x(i4)
      dv3(2) = y(i3) - y(i4)
      dv3(3) = z(i3) - z(i4)
      bl3 = sqrt(sum(dv3(:)**2))
      dv3(:) = dv3(:)/bl3
!     get the cross product of the normed bond vectors (cp1)
      call crossprod3(dv3,dv2,cp1)
!     get the dot product (acos) of the normed bond vectors
      np1 = sum(dv2(:)*dv3(:))
      if (abs(np1).ge.1.0) then
        write(ilog,*) 'Fatal. This indicates an instable simulation,&
 & a bug, or bad input in genxyz(...). Please report if necessary.'
!        write(ilog,68) at,i2,bl,i3,ba,i4,baodi
        call fexit()
      end if
      s1 = sqrt(1.0 - np1**2)
      cp1(:) = cp1(:)/s1
!     get the cross product of the first bond vector with cp1 
      call crossprod3(cp1,dv2,cp2)
!     build the displacement vector and generate coordinates
      dv4(:) = cp2(:)*stm1*ctm2 + cp1(:)*stm1*stm2 - dv2(:)*ctm1
      x(at) = x(i2) + bl*dv4(1)
      y(at) = y(i2) + bl*dv4(2)
      z(at) = z(i2) + bl*dv4(3)
! 
    else if ((chirality.eq.1).OR.(chirality.eq.-1)) then
!
!     with two bond angles
!
      stm1 = ba/RADIAN
      stm2 = baodi/RADIAN
      ctm1 = cos(stm1)
      ctm2 = cos(stm2)
      dv2(1) = x(i3) - x(i2)
      dv2(2) = y(i3) - y(i2)
      dv2(3) = z(i3) - z(i2)
      bl2 = sqrt(sum(dv2(:)**2))
      dv2(:) = dv2(:)/bl2
      dv3(1) = x(i2) - x(i4)
      dv3(2) = y(i2) - y(i4)
      dv3(3) = z(i2) - z(i4)
      bl3 = sqrt(sum(dv3(:)**2))
      dv3(:) = dv3(:)/bl3
      call crossprod3(dv3,dv2,cp1)
      np1 = sum(dv2(:)*dv3(:))
      if (np1.ge.1.0) then
!       i2, i3, and i4 are colinear, which is fatal as orientation around
!       axis becomes ill-defined (in other words, this case HAS to be expressed
!       as a torsional problem)
        write(ilog,70) at,i2,bl,i3,ba,i4,baodi
        call fexit()
      end if
      s2 = 1.0/(1.0 - np1**2)
      c1 = (-ctm2 - np1*ctm1)*s2
      c2 = (ctm1 + np1*ctm2)*s2
      s1 = (c1*ctm2 - c2*ctm1 + 1.0)*s2
      if (s1.gt.0.0) then
        s1 = chirality*sqrt(s1)
      else if (s1.lt.0.0) then
!       the two bond angle setup can potentially be incompatible at planar centers:
!       warn eventually ...
        if (s1.lt.(-sqrt(1.0*(10.0**(-precision(s1)))))) then
          dblba_wrncnt = dblba_wrncnt + 1
          if (dblba_wrncnt.eq.dblba_wrnlmt) then
            write(ilog,69) at,i2,bl,i3,ba,i4,baodi
            write(ilog,*) 'This is warning number #',dblba_wrncnt,' of this type not all of&
 & which may be displayed.'
            if (10.0*dblba_wrnlmt.gt.0.5*HUGE(dblba_wrnlmt)) then
              dblba_wrncnt = 0 ! reset
            else
              dblba_wrnlmt = dblba_wrnlmt*10
            end if
          end if
        end if
        s1 = sqrt(sum((c1*dv3(:) + c2*dv2(:))**2))
        c1 = c1/s1
        c2 = c2/s1
        s1 = 0.0
      end if
      dv4(:) = dv3(:)*c1 + dv2(:)*c2 + cp1(:)*s1
      x(at) = x(i2) + bl*dv4(1)
      y(at) = y(i2) + bl*dv4(2)
      z(at) = z(i2) + bl*dv4(3)
!
    end if
!
  end if
!
end
!
!----------------------------------------------------------------
!
subroutine buildxyz_from2(i2,i3,bl,ba,retv)
!
  use atoms, ONLY: x,y,z
  use math
!
  implicit none
!
  RTYPE, INTENT(IN):: ba,bl
  integer, INTENT(IN):: i2,i3
!
  RTYPE retv(3)
  RTYPE stm1,ctm1,bl2,dv2(3),c1,s2,s1,trx,trz
!
  stm1 = sin(ba/RADIAN)
  ctm1 = cos(ba/RADIAN)
  dv2(1) = x(i2) - x(i3)
  dv2(2) = y(i2) - y(i3)
  dv2(3) = z(i2) - z(i3)
  bl2 = sqrt(sum(dv2(:)**2))
  dv2(:) = dv2(:)/bl2
! rel. xz in the assumed frame
  trx = bl*stm1
  trz = bl2 - bl*ctm1
! shift/rotate into real frame 
  s1 = sqrt(sum(dv2(1:2)**2))
  if (s1.gt.0.0) then
    c1 = dv2(2)/s1
    s2 = dv2(1)/s1
  else
    c1 = 1.0
    s2 = 0.0
  end if
  retv(1) = x(i3) + trx*c1 + trz*s2*s1
  retv(2) = y(i3) - trx*s2 + trz*c1*s1
  retv(3) = z(i3) + trz*dv2(3)
!
end
!
!----------------------------------------------------------------
!
! next are a few wrappers for genxyz(...)
!
!----------------------------------------------------------------
!
subroutine makexyz_forsc(rs)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,firstat,rs,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! even though this is slightly redundant, we'll rebuild the whole residue
! regardless
! note that the definition of chi-angles and "sidechain" implies that a part
! of the molecule is branched off and that this is what is rebuilt here
! (i.e., this would NOT work for generic branched polymers)
! it is very important NOT to have parts in backbone arrays that move as a
! result of chi-perturbation
!
  if (at(rs)%nsc.eq.0) then
    write(ilog,*) 'Fatal. Called makexyz_forsc() with ineligible residue (zero sidechain&
 & atoms).'
    call fexit()
  end if
  imol = molofrs(rs)
  firstat = at(rs)%bb(1)
  if (firstat.lt.atmol(imol,1)+3) firstat = atmol(imol,1)+3
  do i=firstat,at(rs)%sc(at(rs)%nsc)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
  end do
!
end
!
!---------------------------------------------------------------------
!
subroutine makexyz_forbb(rs)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,rs,i2,i3,i4,chirali,firstat,imol
  RTYPE bl,ba,baodi
!
! rebuild all atoms from first one of pivot residue
!
  imol = molofrs(rs)
  firstat = at(rs)%bb(1)
  if (firstat.lt.atmol(imol,1)+3) firstat = atmol(imol,1)+3
  do i=firstat,atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chirali = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chirali)
  end do
!
end
!
!---------------------------------------------------------------------
!
! just rebuilds all atoms in a single residue
!
subroutine makexyz_forset(ilo,ihi)
!
  use iounit
  use atoms, ONLY: atmres
  use sequen, ONLY: molofrs
  use zmatrix
  use molecule
!
  implicit none
!
  integer, INTENT(IN):: ilo,ihi
!
  integer i,i2,i3,i4,chirali,firstat,imol
  RTYPE bl,ba,baodi
!
  imol = molofrs(atmres(ilo))
  firstat = ilo
  if (ilo.lt.atmol(imol,1)+3) firstat = atmol(imol,1)+3
  do i=firstat,ihi
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chirali = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chirali)
  end do
!
end
!
!---------------------------------------------------------------------
! 
! build all atoms in the molecule except reference frame atoms (first 3) 
! 
subroutine makexyz_formol(imol)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
  use atoms
!
  implicit none
!
  integer i,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! never build for first three atoms of a molecule (since those are only affected
! by rigid body moves ("feature" of the N->C building paradigm)) 
! alternative solution if this not feasible: maintaining a ref. array/frame
! for each molecule
! 
  do i=atmol(imol,1)+3,atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
  end do
!
end
!
!---------------------------------------------------------------------
!
! build all atoms in molecule except first one
!
subroutine makexyz_formol2(imol)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! initially, the first three atoms have to be built of course, which is done by this fxn
! 
  do i=atmol(imol,1),atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
  end do
!
end
!
!---------------------------------------------------------------------
!
! build atoms 2 and/or 3 in reference frame
!
subroutine makexyz_formol4(imol,doi1,doi2,doi3)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
  use atoms, ONLY: x,y,z
!
  implicit none
!
  integer i2,i3,imol
  RTYPE bl,ba,dv(3)
  logical doi1,doi2,doi3
!
! allow some exceptions for molecules that require hydrogens in their first few atoms (water, methane, ...)
! heavy atom is read-in -> orientation of rest is "by default"
  if ((doi1.EQV..true.).AND.(doi2.EQV..false.).AND.(doi3.EQV..false.)) then
    i2 = atmol(imol,1)+1
    i3 = atmol(imol,1)+2
    bl = blen(i2)
    ba = bang(i3)
    call buildxyz_from2(i2,i3,bl,ba,dv) 
    x(atmol(imol,1)) = dv(1)
    y(atmol(imol,1)) = dv(2)
    z(atmol(imol,1)) = dv(3)
  else if ((doi1.EQV..false.).AND.(doi2.EQV..true.).AND.(doi3.EQV..false.)) then
    i2 = atmol(imol,1)
    i3 = atmol(imol,1)+2
    dv(1) = x(i3) - x(i2)
    dv(2) = y(i3) - y(i2)
    dv(3) = z(i3) - z(i2)
    x(atmol(imol,1)+1) = x(i2) + 0.5*dv(1) - 0.8*dv(2)/sqrt(dv(2)**2+dv(1)**2)
    y(atmol(imol,1)+1) = y(i2) + 0.5*dv(2) + 0.8*dv(1)/sqrt(dv(2)**2+dv(1)**2)
    z(atmol(imol,1)+1) = z(i2) + 0.5*dv(3)
  else if ((doi1.EQV..false.).AND.(doi2.EQV..false.).AND.(doi3.EQV..true.)) then
    i2 = atmol(imol,1)+1
    i3 = atmol(imol,1)
    bl = blen(i2)
    ba = bang(atmol(imol,1)+2)
    call buildxyz_from2(i2,i3,bl,ba,dv) 
    x(atmol(imol,1)+2) = dv(1)
    y(atmol(imol,1)+2) = dv(2)
    z(atmol(imol,1)+2) = dv(3)
  else if ((doi1.EQV..true.).AND.(doi2.EQV..true.).AND.(doi3.EQV..false.)) then
    i2 = atmol(imol,1)+1
    i3 = atmol(imol,1)+2
    x(i2) = x(i3)
    y(i2) = y(i3)
    z(i2) = z(i3) + blen(i3)
    bl = blen(i2)
    ba = bang(i3)
    call buildxyz_from2(i2,i3,bl,ba,dv) 
    x(atmol(imol,1)) = dv(1)
    y(atmol(imol,1)) = dv(2)
    z(atmol(imol,1)) = dv(3)
  else if ((doi1.EQV..true.).AND.(doi2.EQV..false.).AND.(doi3.EQV..true.)) then
    i2 = atmol(imol,1)+1
    i3 = atmol(imol,1)+2
    x(i3) = x(i2)
    y(i3) = y(i2)
    z(i3) = z(i2) + blen(i3)
    bl = blen(i2)
    ba = bang(i3)
    call buildxyz_from2(i2,i3,bl,ba,dv) 
    x(atmol(imol,1)) = dv(1)
    y(atmol(imol,1)) = dv(2)
    z(atmol(imol,1)) = dv(3)
  else if ((doi1.EQV..false.).AND.(doi2.EQV..true.).AND.(doi3.EQV..true.)) then
    i2 = atmol(imol,1)+1
    i3 = atmol(imol,1)
    x(i2) = x(i3)
    y(i2) = y(i3)
    z(i2) = z(i3) + blen(i2)
    bl = blen(atmol(imol,1)+2)
    ba = bang(atmol(imol,1)+2)
    call buildxyz_from2(i2,i3,bl,ba,dv) 
    x(atmol(imol,1)+2) = dv(1)
    y(atmol(imol,1)+2) = dv(2)
    z(atmol(imol,1)+2) = dv(3)
  else 
    call fexit()
  end if
!
end
!
!---------------------------------------------------------------------
!  
subroutine makexyz_formol3(imol)
!
  use iounit
  use sequen
  use polypep
  use zmatrix
  use molecule
!
  implicit none
!
  integer i,i2,i3,i4,chiral,imol
  RTYPE bl,ba,baodi
!
! allow some exceptions for molecules that require hydrogens in their first few atoms (water, methane, ...)
! heavy atom is read-in -> orientation of rest is "by default"
! 
  do i=atmol(imol,1)+1,atmol(imol,2)
    i2 = iz(1,i)
    i3 = iz(2,i)
    i4 = iz(3,i)
    chiral = iz(4,i)
    bl = blen(i)
    ba = bang(i)
    baodi = ztor(i)
    call genxyz(i,i2,bl,i3,ba,i4,baodi,chiral)
  end do
!
end
!
!---------------------------------------------------------------------
!
! simplest of all wrappers
!   
subroutine makexyz2()
!
  use iounit
  use molecule
!
  implicit none
!
  integer imol
!
  do imol=1,nmol
    call makexyz_formol2(imol)
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this concludes the wrappers for genxyz(...)
!
!-----------------------------------------------------------------------
!
! the next routine is related in that it also builds coordinates
! rather than following the Z-matrix hierarchy, however, it assumes only a single
! rotation takes place and applies this by means of a quaternion
!
subroutine quatxyz_forrotlst(ati,alC,dw,tpi)
!
  use iounit
  use sequen
  use zmatrix
  use molecule
  use atoms
  use system
  use forces
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat
#endif
!
  implicit none
!
  integer, INTENT(IN):: ati,tpi
  logical, INTENT(IN):: alC
  RTYPE, INTENT(IN):: dw
!
  integer i,k,imol,ii,ik,shf,shf2
  RTYPE axl(3),mmm,q(4),qinv(4),refv(3),hlp1
  logical setalC
#ifdef ENABLE_THREADS
  integer sta,incr
!
  if (tpi.gt.0) then
    sta = tpi-1
    incr = thrdat%maxn
  else
    sta = 0
    incr = 1
  end if
#endif
!
  if ((izrot(ati)%alsz.le.0).OR.(allocated(izrot(ati)%rotis).EQV..false.)) return
!
  refv(1) = x(iz(1,ati))
  refv(2) = y(iz(1,ati))
  refv(3) = z(iz(1,ati))
  axl(1) = x(iz(1,ati)) - x(iz(2,ati))
  axl(2) = y(iz(1,ati)) - y(iz(2,ati))
  axl(3) = z(iz(1,ati)) - z(iz(2,ati))
  if (alC.EQV..true.) axl(:) = -1.0*axl(:)
!
  mmm = 1.0/sqrt(sum(axl(:)*axl(:)))
  axl(:) = axl(:)*mmm
!
! construct the quaternion
  q(1) = cos(0.5*dw)
  hlp1 = sin(0.5*dw)
  q(2:4) = axl(1:3)*hlp1
!
! the rotation set may already be flipped, detect this with the help of iz(3,ati)
  setalC = .false.
  do k=izrot(ati)%alsz,1,-1
    ii = izrot(ati)%rotis(k,1)
    ik = izrot(ati)%rotis(k,2)
    if ((iz(3,ati).ge.ii).AND.(iz(3,ati).le.ik)) then
      setalC = .true.
      exit
    end if
  end do
  if (alC.EQV..true.) then
    if (setalC.EQV..true.) then ! setalC true means a C-terminal rot list
      setalC = .false.
    else
      setalC = .true.
    end if
  end if
!  
! construct the multiplicative inverse
  mmm = 1.0/sum(q(:)*q(:))
  if (mmm.le.0.0) then
    write(ilog,*) 'Fatal error. Got a 0-quat. in quatxyz_forrotlst(...).'
    call fexit()
  end if
  qinv(1) = q(1)*mmm
  qinv(2:4) = -q(2:4)*mmm

  if (setalC.EQV..true.) then
    shf = 1
    shf2 = 1
    imol = molofrs(atmres(ati))
    if (izrot(ati)%rotis(1,1).eq.atmol(imol,1)) shf = 2
    if (izrot(ati)%rotis(izrot(ati)%alsz,2).eq.atmol(imol,2)) shf2 = 0
    do k=shf,izrot(ati)%alsz+shf2
      if (k.eq.1) then 
        ii = atmol(imol,1)
      else
        ii = izrot(ati)%rotis(k-1,2)+1
      end if
      if (k.gt.izrot(ati)%alsz) then 
        ik = atmol(imol,2)
      else
        ik = izrot(ati)%rotis(k,1)-1
      end if
#ifdef ENABLE_THREADS
      do i=ii+sta,ik,incr
#else
      do i=ii,ik
#endif
        call quat_conjugate3(q,qinv,i,refv)
      end do
    end do
  else
    do k=1,izrot(ati)%alsz
      ii = izrot(ati)%rotis(k,1)
      ik = izrot(ati)%rotis(k,2)
#ifdef ENABLE_THREADS
      do i=ii+sta,ik,incr
#else
      do i=ii,ik
#endif
        call quat_conjugate3(q,qinv,i,refv)
      end do
    end do
  end if
!
end
!
!---------------------------------------------------------------------
!
! for this next routine, see description below
! note this routine relies on (rather unfortunate) explicit listing
! of residue types via name -> needs to be appended if new residues
! are created
!
subroutine find_ntorsn()
!
  use fyoc
  use iounit
  use sequen
  use torsn
  use zmatrix
  use aminos
  use molecule
  use movesets
  use atoms
  use cutoffs, ONLY: rsinfo
!
  implicit none
!
  integer rs
  integer kk,mty,nmt,curmol,buk
  integer, ALLOCATABLE:: molnrs(:)
!
  ntorsn = 0
  ntorpuck = 0
  ndyntorsn = 0
  allocate(molnrs(nmol))
!
! go through the backbone and identify rotatable bonds
! now this was originally meant for a more general torsional implementation,
! i.e., to have arbitrary defined torsions
! since, however, the MC sampler needs to distinguish (efficiencies, movesets,
! etc., all depend on type of torsion) them anyway, we use more descript pointer
! arrays (see fyoc.i) directly into the relevant z-matrix array
! -> this function only serves to establish "ntorsn", "ndyntorsn", and "ntorpuck" and to populate rsinfo(:,4)
!
  curmol = -1
!
  do rs=1,nseq
!
    rsinfo(rs,4) = ntorpuck + 1
    if (molofrs(rs).ne.curmol) then
      curmol = molofrs(rs)
      if (atmol(molofrs(rs),1).eq.atmol(molofrs(rs),2)) then
        nmt = 3 ! will cycle o'TOR',ut anyway
      else if ((atmol(molofrs(rs),2)-atmol(molofrs(rs),2)).eq.1)&
 &                                                             then
        nmt = 5 ! will cycle out anyway (just for clarity)
      else
        nmt = 6
      end if
    end if
    mty = moltypid(molofrs(rs))
    buk = ntorpuck
!
    if (wline(rs).gt.0) then
      if (izrot(wline(rs))%alsz.gt.0) then
        ntorsn = ntorsn + 1
        ntorpuck = ntorpuck + 1
        nmt = nmt + 1
        wnr(rs) = nmt
        if (moltyp(mty,1).eq.molofrs(rs)) then
          ntormol(mty) = ntormol(mty) + 1
        end if
      end if
    end if
!
    if (fline(rs).gt.0) then
      if (izrot(fline(rs))%alsz.gt.0) then
        ntorsn = ntorsn + 1
        ntorpuck = ntorpuck + 1
        nmt = nmt + 1
        fnr(rs) = nmt
        if (moltyp(mty,1).eq.molofrs(rs)) then
          ntormol(mty) = ntormol(mty) + 1
        end if
!      else
!        ndyntorsn = ndyntorsn + 1
      end if
    end if
!
    if ((pucline(rs).gt.0).OR.(nucsline(6,rs).gt.0)) then
      ntorpuck = ntorpuck + 7
    end if
!
    if (yline(rs).gt.0) then
      if (izrot(yline(rs))%alsz.gt.0) then
        ntorsn = ntorsn + 1
        ntorpuck = ntorpuck + 1
        nmt = nmt + 1
        ynr(rs) = nmt
        if (moltyp(mty,1).eq.molofrs(rs)) then
          ntormol(mty) = ntormol(mty) + 1
        end if
      end if
!
    end if
!
    do kk=1,nnucs(rs)
      if (izrot(nucsline(kk,rs))%alsz.le.0) cycle
      ntorsn = ntorsn + 1
      ntorpuck = ntorpuck + 1
      nmt = nmt + 1
      nucsnr(kk,rs) = nmt
      if (moltyp(mty,1).eq.molofrs(rs)) then
        ntormol(mty) = ntormol(mty) + 1
      end if
    end do
!
    do kk=1,nchi(rs)
      if (izrot(chiline(kk,rs))%alsz.le.0) cycle
      ntorsn = ntorsn + 1
      ntorpuck = ntorpuck + 1
      nmt = nmt + 1
      chinr(kk,rs) = nmt
      if (moltyp(mty,1).eq.molofrs(rs)) then
        ntormol(mty) = ntormol(mty) + 1
      end if
    end do
!
    if (buk.eq.ntorpuck) notors(rs) = .true.
    molnrs(molofrs(rs)) = nmt
!
  end do
!
  do kk=1,unslst%nr
    rs = atmres(iz(1,unslst%idx(kk)))
    notors(rs) = .false.
    ntorsn = ntorsn + 1
    ntorpuck = ntorpuck + 1
    molnrs(molofrs(rs)) = molnrs(molofrs(rs)) + 1
    if (moltyp(moltypid(molofrs(rs)),1).eq.molofrs(rs)) then
      if (othidxmol(moltypid(molofrs(rs))).eq.0) othidxmol(moltypid(molofrs(rs))) = molnrs(molofrs(rs))
      ntormol(moltypid(molofrs(rs))) = ntormol(moltypid(molofrs(rs))) + 1
    end if
  end do
!
  do kk=1,unklst%nr
    rs = atmres(iz(1,unklst%idx(kk)))
    notors(rs) = .false.
    ntorsn = ntorsn + 1
    ntorpuck = ntorpuck + 1
    molnrs(molofrs(rs)) = molnrs(molofrs(rs)) + 1
    if (moltyp(moltypid(molofrs(rs)),1).eq.molofrs(rs)) then
      if (othidxmol(moltypid(molofrs(rs))).eq.0) othidxmol(moltypid(molofrs(rs))) = molnrs(molofrs(rs))
      ntormol(moltypid(molofrs(rs))) = ntormol(moltypid(molofrs(rs))) + 1
    end if
  end do
!
  ndyntorsn = ndyntorsn + ntorsn ! includes proline phi
!
  deallocate(molnrs)
!
end
!
!-----------------------------------------------------------------------
!
! transfer 
!
subroutine transfer_torsn(mode,tvec,tpi)
!
  use fyoc
  use iounit
  use sequen
  use zmatrix
  use system
  use polypep
  use torsn
  use movesets
#ifdef ENABLE_THREADS
  use threads
  use cutoffs, ONLY: rsinfo
#endif
!
  implicit none
!
  integer, INTENT(IN):: mode,tpi
!
  integer rs,ttc,shf,shfbu
  integer kk
  RTYPE tvec(ntorpuck)
#ifdef ENABLE_THREADS
  integer sta,sto
!
  if (tpi.gt.0) then
    sta = thr_limits(45,tpi)
    sto = thr_limits(46,tpi)
  else
    sta =  1
    sto = nseq
  end if
#endif
!
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
#ifdef ENABLE_THREADS
  if ((sta.ge.1).AND.(sta.le.nseq)) then
    ttc = rsinfo(sta,4) - 1
  else if (sto.ge.sta) then
    write(ilog,*) 'Fatal. Inconsistent thread bounds in transfer_torsn. This is a bug.'
    call fexit()
  end if
  do rs=sta,sto
#else
  ttc = 0
  do rs=1,nseq
#endif
!
    if (wline(rs).gt.0) then
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = omega(rs)
      if (mode.eq.2) omega(rs) = tvec(ttc)
    end if
!
    if ((fline(rs).gt.0).AND.(seqflag(rs).ne.5)) then
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = phi(rs)
      if (mode.eq.2) phi(rs) = tvec(ttc)
    end if
!
    if (yline(rs).gt.0) then
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = psi(rs)
      if (mode.eq.2) psi(rs) = tvec(ttc)
    end if
!
    do kk=1,nnucs(rs)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = nucs(kk,rs)
      if (mode.eq.2) nucs(kk,rs) = tvec(ttc)
    end do
!
    do kk=1,nchi(rs)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = chi(kk,rs)
      if (mode.eq.2) chi(kk,rs) = tvec(ttc)
    end do
!
    if (((seqflag(rs).eq.5).AND.(fline(rs).gt.0)).OR.((seqpolty(rs).eq.'N').AND.(nucsline(6,rs).gt.0))) then
      if (seqpolty(rs).eq.'N') then
        shfbu = shf
        shf = 1
      end if
!
      if (seqpolty(rs).eq.'N') then
        ttc = ttc + 1
        if (mode.eq.1) tvec(ttc) = ztor(nucsline(6,rs))
        if (mode.eq.2) ztor(nucsline(6,rs)) = tvec(ttc)
      else
        ttc = ttc + 1
        if (mode.eq.1) tvec(ttc) = ztor(fline(rs))
        if (mode.eq.2) then
          ztor(fline(rs)) = tvec(ttc)
          phi(rs) = tvec(ttc)
        end if
      end if
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = ztor(at(rs)%sc(2-shf))
      if (mode.eq.2) ztor(at(rs)%sc(2-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = ztor(at(rs)%sc(3-shf))
      if (mode.eq.2) ztor(at(rs)%sc(3-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = ztor(at(rs)%sc(4-shf))
      if (mode.eq.2) ztor(at(rs)%sc(4-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = bang(at(rs)%sc(2-shf))
      if (mode.eq.2) bang(at(rs)%sc(2-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = bang(at(rs)%sc(3-shf))
      if (mode.eq.2) bang(at(rs)%sc(3-shf)) = tvec(ttc)
      ttc = ttc + 1
      if (mode.eq.1) tvec(ttc) = bang(at(rs)%sc(4-shf))
      if (mode.eq.2) bang(at(rs)%sc(4-shf)) = tvec(ttc)
!
      if (seqpolty(rs).eq.'N') shf = shfbu
    end if
!
  end do
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    ttc = ntorpuck - unslst%nr - unklst%nr + thr_limits(63,tpi) - 1
    sta = thr_limits(63,tpi)
    sto = thr_limits(64,tpi)
  else
    sta = 1
    sto = unslst%nr
  end if
  do rs=sta,sto
#else
  do rs=1,unslst%nr
#endif
    ttc = ttc + 1
    if (mode.eq.1) tvec(ttc) = ztor(unslst%idx(rs))
    if (mode.eq.2) ztor(unslst%idx(rs)) = tvec(ttc)
  end do
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    ttc = ntorpuck - unklst%nr + thr_limits(65,tpi) - 1
    sta = thr_limits(65,tpi)
    sto = thr_limits(66,tpi)
  else
    sta = 1
    sto = unklst%nr
  end if
  do rs=sta,sto
#else
  do rs=1,unklst%nr
#endif
    ttc = ttc + 1
    if (mode.eq.1) tvec(ttc) = ztor(unklst%idx(rs))
    if (mode.eq.2) ztor(unklst%idx(rs)) = tvec(ttc)
  end do
!
end
!
!-----------------------------------------------------------------------------
!
! this subroutine populates the 1-2 bonded arrays which form the basis
! for the 1-3 and 1-4 ones (see below)
!
! strategy outline: every molecule effectively maintains an independent Z-matrix,
!                   i.e., would one write out the Z-matrix, the first three atoms
!                   of each molecule would indeed not have a well-defined reference
!                   this is more meaningful than the alternative solution to simply
!                   chain the Z-matrix through molecular boundaries leading to
!                   quite meaningless internals which are difficult to maintain
!
!
subroutine setup_connectivity_2()
!
  use zmatrix
  use iounit
  use molecule
  use atoms
  use aminos
  use sequen, ONLY: seqtyp
  use params, ONLY: bio_code
!
  implicit none
!
  integer imol,at,at2,i,j,k
!
  i12(:,:) = 0
  n12(:) = 0
!
 33 format ('(encountered for atom #',i10,' of biotype ',a3,' in residue #',i8,' of type ',a3,')')
  do imol=1,nmol
    do at=atmol(imol,1)+1,atmol(imol,2)
      at2 = iz(1,at)
!     sanity, should be handled by the +1-offset
      if (at2.le.0) then
        write(ilog,*) 'Warning! Missing reference atom for second atom&
 & in molecule ',imol,' (atom #',at,' returned ',at2,'). This is mos&
 &t likely a bug. Please report.'
        cycle
      end if
      n12(at) = n12(at) + 1
      n12(at2) = n12(at2) + 1
      if ((n12(at).gt.MAXVLNC).OR.(n12(at2).gt.MAXVLNC)) then
        write(ilog,*) 'Fatal. Maximum valency was exceeded. This is either a bug or caused by &
 &a Z matrix inferred from an input file with unsuitable atom order (unsupported residue).'
        if (n12(at).gt.MAXVLNC) write(ilog,33) at,bio_code(b_type(at)),atmres(at),amino(seqtyp(atmres(at)))
        if (n12(at2).gt.MAXVLNC) write(ilog,33) at2,bio_code(b_type(at2)),atmres(at2),amino(seqtyp(atmres(at2)))
        call fexit()
      end if
      i12(n12(at),at) = at2
      i12(n12(at2),at2) = at
    end do
  end do
!
! cyclic systems require the inclusion of extra 1-2 pairs not apparent
! from the simple Z-matrix term
! these were set up previously in the z-matrix generating routines
!
  do i=1,nadd
    do j=1,2
      k = iadd(j,i)
      n12(k) = n12(k) + 1
      if (n12(k).gt.MAXVLNC) then
        write(ilog,*) 'Fatal. Maximum valency was exceeded. This is either a bug or caused by &
 &a Z matrix inferred from an input file with unsuitable atom order (unsupported residue).'
        write(ilog,33) k,bio_code(b_type(k)),atmres(k),amino(seqtyp(atmres(k)))
        call fexit()
      end if
      i12(n12(k),k) = iadd(3-j,i)
    end do
  end do
!
! finally sort the list
!
  do at=1,n
    call sort(n12(at),i12(1,at))
  end do
!
end
!
!-----------------------------------------------------------------------
!
! this subroutine populates the 1-3 and 1-4 bonded arrays
!
subroutine setup_connectivity_34()
!
  use atoms
!
  implicit none
!
  integer at,i,ii,j,jj,k,kk,l
  logical docyc
!
! first 1-3
!
  i13(:,:) = 0
  do at=1,n
    n13(at) = 0
    do i=1,n12(at)
      ii = i12(i,at)
      do j=1,n12(ii)
        jj = i12(j,ii)
        docyc = .false.
        if (jj.eq.at) cycle
        do k=1,n12(at)
          if (jj.eq.i12(k,at)) then
            docyc = .true.
            exit
          end if
        end do
        if (docyc.EQV..true.) cycle
        n13(at) = n13(at) + 1
        i13(n13(at),at) = jj
      end do
    end do
    call sort(n13(at),i13(1,at))
  end do
!
! now 1-4
!
  i14(:,:) = 0
  do at=1,n
    n14(at) = 0
    do i=1,n12(at)
      ii = i12(i,at)
      do j=1,n12(ii)
        jj = i12(j,ii)
        do k=1,n12(jj)
          kk = i12(k,jj)
          docyc = .false.
          if (kk.eq.at) cycle
          do l=1,n12(at)
            if (kk.eq.i12(l,at)) then
              docyc = .true.
              exit
            end if
          end do
          do l=1,n13(at)
            if (kk.eq.i13(l,at)) then
              docyc = .true.
              exit
            end if
          end do
          if (docyc.EQV..true.) cycle
          n14(at) = n14(at) + 1
          i14(n14(at),at) = kk
        end do
      end do
    end do
    call sort(n14(at),i14(1,at))
  end do
!
end
!
!-----------------------------------------------------------------------
!
subroutine valence_check()
!
  use atoms
  use iounit
  use params
!
  implicit none
!
  integer at,k,kk,t1,t2
  logical, ALLOCATABLE:: chkd(:)
!
! since the parameter file has valence information, check whether
! any exceptions occur: good bug-check against builder+parameters
!    
 75   format('Fatal. Atom ',i7,' (biotype ',i5,'), which was assigned LJ&
 &-type ',i4,', has ',i2,' bound partner atoms instead of the ',&
 &i2,' defined for it in the parameter file (',a,'). This is the first error of this type&
 & for LJ-type ',i4,'. Further warnings omitted.')
!
  allocate(chkd(n_ljtyp))
  chkd(:) = .false.
  call strlims(paramfile,t1,t2)
!
  do at=1,n
    k = attyp(at)
    kk = bio_ljtyp(b_type(at))
    if (n12(at).ne.lj_val(k)) then
      if ((kk.ne.k).AND.(lj_patched.EQV..true.).AND.(lj_val(kk).eq.n12(at))) cycle
      if (chkd(k).EQV..false.) then
        write(ilog,75) at,b_type(at),k,n12(at),lj_val(k),&
 &paramfile(t1:t2),k
        chkd(k) = .true.
      end if
      if (be_unsafe.EQV..false.) call fexit()
    end if
  end do
!
  deallocate(chkd)
!
end
!
!-------------------------------------------------------------------------
!
! subroutine to determine whether an atom is part of a ring backbone
! it expects an integer array of size 10 to hold up to 10 rings
!
  subroutine determine_upto6ring(ati,nrgs,rnglst)
!
  use atoms
!
  implicit none
!
  integer ati,nrgs,rnglst(10),i,j,k,l,m
  logical isok
!
  nrgs = 0
!
! 3 
  do i=1,n12(ati)
    do j=i+1,n12(ati)
      do k=1,n12(i12(i,ati))
        if (i12(k,i12(i,ati)).eq.i12(j,ati)) then
          nrgs = nrgs + 1
          rnglst(nrgs) = 3
        end if
      end do
    end do
  end do
! 4 - fails for bicyclo [1.1.0]butane
  do i=1,n13(ati)
    do j=i+1,n13(ati)
      if (i13(i,ati).eq.i13(j,ati)) then
        nrgs = nrgs + 1
        rnglst(nrgs) = 4
      end if
    end do
  end do
! 5 - filters out atoms connected to 3-rings
  do i=1,n13(ati)
    do j=i+1,n13(ati)
      isok = .false.
      do k=1,n12(i13(i,ati))
        if (i12(k,i13(i,ati)).eq.i13(j,ati)) then
          isok = .true.
          exit
        end if
      end do
      if (isok.EQV..false.) cycle
      do k=1,n12(i13(i,ati))
        do l=1,n12(i13(j,ati))
          if (i12(k,i13(i,ati)).eq.i12(l,i13(j,ati))) isok = .false.
        end do
      end do
      if (isok.EQV..false.) cycle
      nrgs = nrgs + 1
      rnglst(nrgs) = 5
    end do
  end do
! 6 - filters out atoms connected to 4-rings
  do i=1,n13(ati)
    do j=i+1,n13(ati)
      isok = .false.
      do k=1,n13(i13(i,ati))
        if (i13(k,i13(i,ati)).eq.i13(j,ati)) then
          isok = .true.
          exit
        end if
      end do
      if (isok.EQV..false.) cycle
      do k=1,n12(ati)
        do l=1,n12(i13(i,ati))
          do m=1,n12(i13(j,ati))
            if ((i12(k,ati).eq.i12(l,i13(i,ati))).AND.(i12(k,ati).eq.i12(m,i13(j,ati)))) isok = .false.
          end do
        end do
      end do
      if (isok.EQV..false.) cycle
      nrgs = nrgs + 1
      rnglst(nrgs) = 6
    end do
  end do
!
  end 
!
!---------------------------------------------------------------------------
!
! this routine scans all bonds in a molecule, identifies those that are freely rotatable, picks a single
! atom to describe rotation around this bond (f,yline/f,yline2 exceptions apply), and allocates memory and
! stores lists of atoms rotating (default building direction) upon changes of the corresponding dihedral Z-matrix entry
! lists are stored as stretches of consecutive atoms
! this information can also be useful in identifying angles that should be eligible, but
! fail due to Z-matrix restrictions (usually caused by setting up multiple rotating atoms with iz(4,i) being 0)
! for this to work with hard crosslinks (crlk_mode is 2), crosslinks would have to first be established, intermolecular
! ones coperturbed (crosslink_follow), and checks as per genints(...) may have to be altered
!
subroutine find_rotlsts(imol)
!
  use zmatrix
  use sequen
  use molecule
  use atoms
  use polypep 
  use fyoc
  use inter
  use aminos
  use params
  use iounit
!
  implicit none
!
  integer imol,rs,tnat,ati,k,i,j,l,zid,zid2,critrs,rslo,rshi,atlo,athi,atcnt,sthlp(4),stsz(4),stlo(4),ll
  RTYPE incr,epsa,epsl,zbu,zbu2
  integer, ALLOCATABLE:: blst(:,:),startrs(:,:)
  RTYPE, ALLOCATABLE:: coords(:,:),lvec(:,:),ivec(:,:),tvec(:,:),avec(:,:)
  logical elig,isnew,proflag
!
  izrot(atmol(imol,1):atmol(imol,2))%alsz = 0
  izrot(atmol(imol,1):atmol(imol,2))%alsz2 = 0
!
 95 format('Bond in question: ',a4,'(',i10,')--',a4,'(',i10,')',/,'Terminal atoms in question: ',a4,'(',i10,') / ',a4,'(',i10,')')
! first make a list of all non-terminal bonds
  incr = 27.
  epsl = 1.0e-6
  epsa = 1.0e-4
  tnat = atmol(imol,2) - atmol(imol,1) + 1
  allocate(blst(tnat,MAXVLNC+1))
  blst(:,1) = 0
!
  do ati=atmol(imol,1),atmol(imol,2)
    if (n12(ati).le.1) cycle
    do k=1,n12(ati)
      if ((ati.lt.i12(k,ati)).AND.(n12(i12(k,ati)).gt.1)) then
        blst(ati-atmol(imol,1)+1,1) = blst(ati-atmol(imol,1)+1,1) + 1
        blst(ati-atmol(imol,1)+1,1+blst(ati-atmol(imol,1)+1,1)) = i12(k,ati)
      end if
    end do
  end do
!
! check whether the Z matrix contains any unexpected references to residues within this molecule
  critrs = rsmol(imol,2)+1
  do ati=atmol(imol,1),atmol(imol,2)
    if (iz(1,ati).le.0) cycle
    ll = 3
    if (iz(3,ati).le.0) ll = 2
    if (iz(2,ati).le.0) ll = 1
    j = max(maxval(atmres(iz(1:ll,ati))),atmres(ati))
    k = min(minval(atmres(iz(1:ll,ati))),atmres(ati))
    if ((j.gt.atmres(ati)).OR.(abs(j-k).gt.1)) then
      critrs = max(rsmol(imol,1),k-1)
      exit
    end if
  end do
!
! for each bond, check whether it is rotatable and whether it is already in a pointer array
  allocate(coords(tnat,3))
  coords(1:tnat,1) = x(atmol(imol,1):atmol(imol,2))
  coords(1:tnat,2) = y(atmol(imol,1):atmol(imol,2))
  coords(1:tnat,3) = z(atmol(imol,1):atmol(imol,2))
  allocate(startrs(rsmol(imol,2)-rsmol(imol,1)+1,4))
  startrs(:,:) = 0
  do rs=rsmol(imol,1)+1,rsmol(imol,2)
    startrs(rs-rsmol(imol,1)+1,1) = startrs(rs-rsmol(imol,1),1) + nrsbl(rs-1)
    startrs(rs-rsmol(imol,1)+1,2) = startrs(rs-rsmol(imol,1),2) + nrsba(rs-1)
    startrs(rs-rsmol(imol,1)+1,3) = startrs(rs-rsmol(imol,1),3) + nrsimpt(rs-1)
    startrs(rs-rsmol(imol,1)+1,4) = startrs(rs-rsmol(imol,1),4) + nrsdi(rs-1)
  end do
  allocate(lvec(max(1,sum(nrsbl(rsmol(imol,1):rsmol(imol,2)))),2))
  allocate(avec(max(1,sum(nrsba(rsmol(imol,1):rsmol(imol,2)))),2))
  allocate(ivec(max(1,sum(nrsimpt(rsmol(imol,1):rsmol(imol,2)))),2))
  allocate(tvec(max(1,sum(nrsdi(rsmol(imol,1):rsmol(imol,2)))),2))
  lvec(:,:) = 0.0
  avec(:,:) = 0.0
  ivec(:,:) = 0.0
  tvec(:,:) = 0.0
  sthlp(:) = startrs(1,:)
  stsz(:) = startrs(rsmol(imol,2)-rsmol(imol,1)+1,:)
  stsz(1) = max(1,stsz(1) + nrsbl(rsmol(imol,2)))
  stsz(2) = max(1,stsz(2) + nrsba(rsmol(imol,2)))
  stsz(3) = max(1,stsz(3) + nrsimpt(rsmol(imol,2)))
  stsz(4) = max(1,stsz(4) + nrsdi(rsmol(imol,2)))
  call genints(rsmol(imol,1),rsmol(imol,2),lvec(1:stsz(1),2),avec(1:stsz(2),2),ivec(1:stsz(3),2),tvec(1:stsz(4),2),sthlp,stsz)
!
  do i=1,tnat
    do k=1,blst(i,1)
      zid = 0
      zid2 = 0
      isnew = .false.
      proflag = .false.
      do rs=max(rsmol(imol,1),atmres(blst(i,1+k))-2),min(rsmol(imol,2),atmres(blst(i,1+k))+2)
        if (wline(rs).gt.0) then
          if (((iz(1,wline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,wline(rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,wline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,wline(rs)).eq.blst(i,1+k)))) then
            zid = wline(rs)
          end if
        end if
        if (fline(rs).gt.0) then
          if (((iz(1,fline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,fline(rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,fline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,fline(rs)).eq.blst(i,1+k)))) then
            zid = fline(rs)
            if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.25).OR.(seqtyp(rs).eq.32)) then
              proflag = .true.
            end if
            zid2 = fline2(rs)
          end if
        end if
        if (yline(rs).gt.0) then
          if (((iz(1,yline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,yline(rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,yline(rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,yline(rs)).eq.blst(i,1+k)))) then
            zid = yline(rs)
            zid2 = yline2(rs)
          end if
        end if
        do j=1,nchi(rs)
          if (((iz(1,chiline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,chiline(j,rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,chiline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,chiline(j,rs)).eq.blst(i,1+k)))) then
            zid = chiline(j,rs)
          end if
        end do
        do j=1,nnucs(rs)
          if (((iz(1,nucsline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(2,nucsline(j,rs)).eq.blst(i,1+k))).OR.&
 &            ((iz(2,nucsline(j,rs)).eq.(i+atmol(imol,1)-1)).AND.(iz(1,nucsline(j,rs)).eq.blst(i,1+k)))) then
            zid = nucsline(j,rs)
          end if
        end do
      end do
      if (zid.eq.0) then
        do j=1,n12(blst(i,1+k))
          if ((iz(1,i12(j,blst(i,1+k))).eq.blst(i,1+k)).AND.(iz(2,i12(j,blst(i,1+k))).eq.(i+atmol(imol,1)-1))) then
            elig = .true.
            do l=1,n12(blst(i,1+k))
              if (i12(l,blst(i,1+k)).eq.iz(3,i12(j,blst(i,1+k)))) then
                elig = .false.
                exit ! improper setup
              end if
            end do
!           further conditions are dihedral setup and having all reference atoms defined
            if ((elig.EQV..true.).AND.(iz(4,i12(j,blst(i,1+k))).eq.0).AND.(minval(iz(1:3,i12(j,blst(i,1+k)))).gt.0)) then
              zid = i12(j,blst(i,1+k))
              isnew = .true.
              exit
            end if
          end if
        end do
        do j=1,n12(i+atmol(imol,1)-1)
          if ((iz(1,i12(j,(i+atmol(imol,1)-1))).eq.(i+atmol(imol,1)-1)).AND.(iz(2,i12(j,(i+atmol(imol,1)-1))).eq.blst(i,1+k))) then
            elig = .true.
            do l=1,n12((i+atmol(imol,1)-1))
              if (i12(l,(i+atmol(imol,1)-1)).eq.iz(3,i12(j,(i+atmol(imol,1)-1)))) then
                elig = .false.
                exit ! improper setup
              end if
            end do
!           further conditions are dihedral setup and having all reference atoms defined
            if ((elig.EQV..true.).AND.(iz(4,i12(j,(i+atmol(imol,1)-1))).eq.0).AND.&
 &              (minval(iz(1:3,i12(j,(i+atmol(imol,1)-1)))).gt.0)) then
              zid2 = i12(j,(i+atmol(imol,1)-1))
              isnew = .true.
              exit
            end if
          end if
        end do
      end if
      if ((zid.le.0).AND.(zid2.le.0)) cycle
!
!     determine ranges we need to check
      ll = 3
      if (iz(3,blst(i,1+k)).le.0) ll = 2
      if (iz(2,blst(i,1+k)).le.0) ll = 1
      j = max(maxval(atmres(iz(1:ll,blst(i,1+k)))),atmres(blst(i,1+k)))
      l = min(minval(atmres(iz(1:ll,blst(i,1+k)))),atmres(blst(i,1+k)))

      rslo = min(critrs,l)
      if (l.ge.critrs) then
        rshi = rsmol(imol,2)
      else
        rshi = min(rsmol(imol,2),j+1)
      end if
      atlo = at(rslo)%bb(1)
      l = min(rshi+2,rsmol(imol,2))
      athi = at(l)%bb(1)+at(l)%nbb+at(l)%nsc-1
      atcnt = athi-atlo+1
      sthlp(:) = startrs(rslo-rsmol(imol,1)+1,:)
      stlo(:) = sthlp(:)+1
      stsz(:) = startrs(rshi-rsmol(imol,1)+1,:)
      stsz(1) = max(1,stsz(1) + nrsbl(rshi))
      stsz(2) = max(1,stsz(2) + nrsba(rshi))
      stsz(3) = max(1,stsz(3) + nrsimpt(rshi))
      stsz(4) = max(1,stsz(4) + nrsdi(rshi))
!
!      write(*,*) i,k,i+atmol(imol,1)-1,blst(i,1+k),isnew
!      if (zid.gt.0) write(*,*) 'trying ',zid,bio_code(b_type(zid)),atmres(zid),amino(seqtyp(atmres(zid))),' for rot ',isnew
!      if (zid2.gt.0) write(*,*) 'also trying ',zid2,bio_code(b_type(zid2)),atmres(zid2),amino(seqtyp(atmres(zid2))),' for rot '
      if (zid.gt.0) then
        zbu = ztor(zid)
        ztor(zid) = ztor(zid) + incr
        if (ztor(zid).gt.180.) ztor(zid) = ztor(zid) - 360.
        if (ztor(zid).le.-180.) ztor(zid) = ztor(zid) + 360.
        call makexyz_forset(atlo,athi)
        call genints(rslo,rshi,lvec(1:stsz(1),1),avec(1:stsz(2),1),ivec(1:stsz(3),1),tvec(1:stsz(4),1),sthlp,stsz)
        elig = .true.
        if ((maxval(abs(lvec(stlo(1):stsz(1),1)-lvec(stlo(1):stsz(1),2))).gt.epsl).OR.&
 &          (maxval(abs(avec(stlo(2):stsz(2),1)-avec(stlo(2):stsz(2),2))).gt.epsa).OR.&
 &          (maxval(abs(ivec(stlo(3):stsz(3),1)-ivec(stlo(3):stsz(3),2))).gt.epsa)) then
          elig = .false.
        end if
        if ((isnew.EQV..false.).AND.(zid2.gt.0).AND.(elig.EQV..true.)) then
          do rs=max(rsmol(imol,1),atmres(blst(i,1+k))-2),min(rsmol(imol,2),atmres(blst(i,1+k))+2)
            if (zid2.eq.fline2(rs)) then
              fline2(rs) = 0
              exit
            else if (zid2.eq.yline2(rs)) then
              yline2(rs) = 0
              exit
            end if
          end do
          if (seqtyp(rs).ne.26) then
            if ((rs.lt.rsmol(molofrs(rs),2)).AND.(seqtyp(min(rs+1,nseq)).eq.26)) then
              write(ilog,*) 'Warning. Existing and rotatable dihedral angle with double setup does not &
 &require second Z-matrix entry for proper rotation. This is likely to be caused by linkage to an unsupported &
 &residue and may change the meaning of phi/psi angles in the affected polypeptide residue.'
              write(ilog,95) bio_code(b_type(i+atmol(imol,1)-1)),i+atmol(imol,1)-1,bio_code(b_type(blst(i,1+k))),blst(i,1+k),&
 &                         bio_code(b_type(zid)),zid,bio_code(b_type(zid2)),zid2
            else
              write(ilog,*) 'Warning. Existing and rotatable dihedral angle with double setup does not &
 &require second Z-matrix entry for proper rotation. This is likely an inconsistency in setting up the Z-matrix or &
 &the result of linkage to a residue of unsupported polymer type.'
              write(ilog,95) bio_code(b_type(i+atmol(imol,1)-1)),i+atmol(imol,1)-1,bio_code(b_type(blst(i,1+k))),blst(i,1+k),&
 &                       bio_code(b_type(zid)),zid,bio_code(b_type(zid2)),zid2
            end if
          end if
          zid2 = 0
        else if ((isnew.EQV..true.).AND.(elig.EQV..true.)) then
          zid2 = 0
        end if
        if (elig.EQV..true.) call gen_rotlst(zid,zid2,imol,atlo,athi,coords(:,:),isnew)
!       restore
        x(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),1)
        y(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),2)
        z(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),3)
        ztor(zid) = zbu
        if ((elig.EQV..false.).AND.((zid2.le.0).OR.(isnew.EQV..true.))) zid = 0
      end if
!
      if ((zid2.gt.0).AND.(isnew.EQV..true.)) then
        zbu = ztor(zid2)
        ztor(zid2) = ztor(zid2) + incr
        if (ztor(zid2).gt.180.) ztor(zid2) = ztor(zid2) - 360.
        if (ztor(zid2).le.-180.) ztor(zid2) = ztor(zid2) + 360.
        call makexyz_forset(atlo,athi)
        call genints(rslo,rshi,lvec(1:stsz(1),1),avec(1:stsz(2),1),ivec(1:stsz(3),1),tvec(1:stsz(4),1),sthlp,stsz)
        elig = .true.
        if ((maxval(abs(lvec(stlo(1):stsz(1),1)-lvec(stlo(1):stsz(1),2))).gt.epsl).OR.&
 &          (maxval(abs(avec(stlo(2):stsz(2),1)-avec(stlo(2):stsz(2),2))).gt.epsa).OR.&
 &          (maxval(abs(ivec(stlo(3):stsz(3),1)-ivec(stlo(3):stsz(3),2))).gt.epsa)) then
          elig = .false.
        end if
        if (elig.EQV..true.) then
          zid = zid2
        end if
        if (elig.EQV..true.) call gen_rotlst(zid,zid2,imol,atlo,athi,coords(:,:),isnew)
        x(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),1)
        y(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),2)
        z(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),3)
        ztor(zid2) = zbu
        if ((elig.EQV..false.).AND.((zid.le.0).OR.(isnew.EQV..true.))) zid2 = 0 
      end if
      if ((zid2.gt.0).AND.(zid.gt.0).AND.(isnew.EQV..false.)) then
        zbu = ztor(zid)
        zbu2 = ztor(zid2)
        ztor(zid) = ztor(zid) + incr
        ztor(zid2) = ztor(zid2) + incr
        if (ztor(zid).gt.180.) ztor(zid) = ztor(zid) - 360.
        if (ztor(zid2).gt.180.) ztor(zid2) = ztor(zid2) - 360.
        if (ztor(zid).le.-180.) ztor(zid) = ztor(zid) + 360.
        if (ztor(zid2).le.-180.) ztor(zid2) = ztor(zid2) + 360.
        call makexyz_forset(atlo,athi)
        call genints(rslo,rshi,lvec(1:stsz(1),1),avec(1:stsz(2),1),ivec(1:stsz(3),1),tvec(1:stsz(4),1),sthlp,stsz)
        elig = .true.
        if ((maxval(abs(lvec(stlo(1):stsz(1),1)-lvec(stlo(1):stsz(1),2))).gt.epsl).OR.&
 &          (maxval(abs(avec(stlo(2):stsz(2),1)-avec(stlo(2):stsz(2),2))).gt.epsa).OR.&
 &          (maxval(abs(ivec(stlo(3):stsz(3),1)-ivec(stlo(3):stsz(3),2))).gt.epsa)) then
          elig = .false.
        end if
        if (elig.EQV..true.) call gen_rotlst(zid,zid2,imol,atlo,athi,coords(:,:),isnew)
        x(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),1)
        y(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),2)
        z(atlo:athi) = coords((atlo-atmol(imol,1)+1):(atlo-atmol(imol,1)+atcnt),3)
        ztor(zid) = zbu
        ztor(zid2) = zbu2
        if (elig.EQV..false.) then
          zid = 0
          zid2 = 0
        end if
      end if
      if ((isnew.EQV..false.).AND.(zid.le.0).AND.(proflag.EQV..false.)) then
        write(ilog,*) 'Fatal. A native dihedral angle appears to not be rotatable freely. This is probably &
 &caused by an unusual linkage between an unsupported and either a supported or an unsupported residue. &
 &The problem may be solved by masking the supported polymer type of unsupported residues, or by renaming &
 &additional residues to be considered as unsupported. It could also indicate a bug in setting up &
 &the Z-matrix.'
        call fexit()
      end if
    end do
  end do
!
  deallocate(tvec)
  deallocate(ivec)
  deallocate(avec)
  deallocate(lvec)
  deallocate(coords)
  deallocate(blst)
!
end 
!
!------------------------------------------------------------------------
!
! compute all internal coordinates in the complete lists (not Z-matrix!)
! this is used to test preservation of covalent geometry 
!
subroutine genints(rlo,rhi,vl,vb,vi,vt,strss,stsz)
!
  use inter
  use molecule
  use sequen
  use fyoc
  use atoms
  use iounit
  use polypep
  use system
!
  implicit none
!
  integer, INTENT(IN):: rlo,rhi,strss(4),stsz(4)
!
  integer refi
  RTYPE vl(stsz(1))
  RTYPE vb(stsz(2))
  RTYPE vi(stsz(3))
  RTYPE vt(stsz(4))
  integer rs,cn(4),i,rs1,rs2
  RTYPE getblen,getbang,getztor
!
  cn(:) = strss(:)
  do rs=rlo,rhi
    if ((disulf(rs).gt.0).AND.(crlk_mode.eq.1)) then ! crosslink treated as restraints has to be skipped
      rs1 = min(rs,disulf(rs))
      rs2 = max(rs,disulf(rs))
      do i=1,nrsbl(rs)
        cn(1) = cn(1) + 1
        if ((maxval(atmres(iaa(rs)%bl(i,1:2))).eq.rs2).AND.(minval(atmres(iaa(rs)%bl(i,1:2))).eq.rs1)) then
          vl(cn(1)) = 0.0
        else
          vl(cn(1)) = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
        end if
      end do
      do i=1,nrsba(rs)
        cn(2) = cn(2) + 1
        if ((maxval(atmres(iaa(rs)%ba(i,1:3))).eq.rs2).AND.(minval(atmres(iaa(rs)%ba(i,1:3))).eq.rs1)) then
          vb(cn(2)) = 0.0
        else
          vb(cn(2)) = getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
        end if
      end do
      do i=1,nrsimpt(rs)
        cn(3) = cn(3) + 1
        if ((maxval(atmres(iaa(rs)%impt(i,1:4))).eq.rs2).AND.(minval(atmres(iaa(rs)%impt(i,1:4))).eq.rs1)) then
          vi(cn(3)) = 0.0
        else
          vi(cn(3)) = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
        end if
      end do
      do i=1,nrsdi(rs)
        cn(4) = cn(4) + 1
        if ((maxval(atmres(iaa(rs)%di(i,1:4))).eq.rs2).AND.(minval(atmres(iaa(rs)%di(i,1:4))).eq.rs1)) then
          vt(cn(4)) = 0.0
        else
          vt(cn(4)) = getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4))
        end if
      end do
!     may need to add back terms
      if (((rs2-rs1).eq.1).AND.(molofrs(rs1).eq.molofrs(rs2))) then
        if (crosslink(crlk_idx(rs1))%itstype.eq.1) then ! CysCys S-S
          if (ua_model.eq.0) then
            refi = at(rs1)%sc(3) ! SG 
          else
            refi = at(rs1)%sc(2) ! SG
          end if
          do i=1,nrsbl(rs)
            if ((maxval(atmres(iaa(rs)%bl(i,1:2))).eq.rs2).AND.(minval(atmres(iaa(rs)%bl(i,1:2))).eq.rs1)) then
              if (minval(abs(iaa(rs)%bl(i,1:2)-refi)).gt.0) then
                vl(cn(1)-nrsbl(rs)+i) = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
              end if 
            end if
          end do
          do i=1,nrsba(rs)
            if ((maxval(atmres(iaa(rs)%ba(i,1:3))).eq.rs2).AND.(minval(atmres(iaa(rs)%ba(i,1:3))).eq.rs1)) then
              if (minval(abs(iaa(rs)%ba(i,1:3)-refi)).gt.0) then
                vb(cn(2)-nrsba(rs)+i) = getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
              end if 
            end if
          end do
          do i=1,nrsimpt(rs)
            if ((maxval(atmres(iaa(rs)%impt(i,1:4))).eq.rs2).AND.(minval(atmres(iaa(rs)%impt(i,1:4))).eq.rs1)) then
              if (minval(abs(iaa(rs)%impt(i,1:4)-refi)).gt.0) then
                vi(cn(3)-nrsimpt(rs)+i) = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
              end if 
            end if
          end do
          do i=1,nrsdi(rs)
            if ((maxval(atmres(iaa(rs)%di(i,1:4))).eq.rs2).AND.(minval(atmres(iaa(rs)%di(i,1:4))).eq.rs1)) then
              if (minval(abs(iaa(rs)%di(i,1:4)-refi)).gt.0) then
                vt(cn(4)-nrsdi(rs)+i) = getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4))
              end if 
            end if
          end do
        end if
      end if
    else if ((crlk_mode.eq.2).AND.(disulf(rs).gt.0)) then
      write(ilog,*) 'Fatal. Hard crosslink constraints are not yet supported in genints(...). This is an &
 &omission bug.'
      call fexit()
    else
      do i=1,nrsbl(rs)
        cn(1) = cn(1) + 1
        vl(cn(1)) = getblen(iaa(rs)%bl(i,1),iaa(rs)%bl(i,2))
      end do
      do i=1,nrsba(rs)
        cn(2) = cn(2) + 1
        vb(cn(2)) = getbang(iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3))
      end do
      do i=1,nrsimpt(rs)
        cn(3) = cn(3) + 1
        vi(cn(3)) = getztor(iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4))
      end do
      do i=1,nrsdi(rs)
        cn(4) = cn(4) + 1
        vt(cn(4)) = getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4))
      end do
    end if
  end do
!
end
!
!----------------------------------------------------------------------
!
! populate rotation list from coordinate set differences
! this assumes the natural building direction (N->C, branches outward)
! the complement set (sans axis atoms) is obtained in straightforward manner
!
subroutine gen_rotlst(z1,z2,imol,atlo,athi,cds,isnew)
!
  use molecule
  use zmatrix
  use atoms
  use iounit
!  use params
!
  implicit none
!
  integer, INTENT(IN):: z1,z2,imol,atlo,athi
  RTYPE, INTENT(IN):: cds(atmol(imol,2)-atmol(imol,1)+1,3)
!
  integer i,zf,nsgs
  RTYPE epsx
  integer, ALLOCATABLE:: tmplst(:,:)
  logical inseg,isnew
!
  allocate(tmplst(atmol(imol,2)-atmol(imol,1)+1,2))
  epsx = 1.0e-7
  nsgs = 0
  inseg = .false.
!
  do i=atlo,athi
    if ((abs(cds(i-atmol(imol,1)+1,1)-x(i)).gt.epsx).OR.(abs(cds(i-atmol(imol,1)+1,2)-y(i)).gt.epsx).OR.&
 &      (abs(cds(i-atmol(imol,1)+1,3)-z(i)).gt.epsx)) then
      if (inseg.EQV..false.) then
        nsgs = nsgs + 1
        tmplst(nsgs,1) = i
        inseg = .true.
      end if
    else
      if (inseg.EQV..true.) then
        tmplst(nsgs,2) = i-1
        inseg = .false.
      end if
    end if
  end do
!
  if (atmres(athi).lt.rsmol(imol,2)) then
!    we did not check all downstream atoms explicitly -> assume that, because i to i-1 connectivity was verified
!    explicitly outside, an open segment implies everything moves downstream whereas a lack of an open
!    segment means that nothing moves downstream (categorically)
    if (inseg.EQV..true.) then ! terminate last stretch if needed
      tmplst(nsgs,2) = atmol(imol,2)
    end if
  else
    if (inseg.EQV..true.) then ! terminate last stretch if needed
      tmplst(nsgs,2) = atmol(imol,2)
    end if
  end if
  if ((z1.gt.0).AND.(z2.gt.0)) then
    zf = z1
  else
    zf = z1
    if (z2.gt.0) zf = z2
  end if
!
  if (nsgs.gt.0) then
    allocate(izrot(zf)%rotis(nsgs,2))
    izrot(zf)%alsz = nsgs
    izrot(zf)%rotis(:,:) = tmplst(1:nsgs,:)
!    write(*,'(i10,a,a3,a,i4)') zf,' of type ',bio_code(b_type(zf)),': ',nsgs
!    do i=1,nsgs
!      write(*,'(i10,1x,i10)') izrot(zf)%rotis(i,1:2)
!    end do
  else if (isnew.EQV..false.) then
    write(ilog,*) 'Fatal. A dihedral angle explicitly defined as a degree of freedom does not rotate any atoms &
 &(Z-matrix line ',zf,'). This is usually the result of introducing colinear geometries, or could be indicative of a bug.'
!    call fexit()
  end if
!
  if (allocated(tmplst).EQV..true.) deallocate(tmplst)
!  
end
!
!----------------------------------------------------------------------------
!
subroutine trans_rotlstdof()
!
  use atoms
  use zmatrix
  use movesets
  use fyoc
  use molecule
  use sequen
!
  implicit none
!
  integer i,rs,rsi,j,imol
  logical skip
!
  unklst%nr = 0
  unslst%nr = 0
  natlst%nr = 0
  do i=1,n
    if (izrot(i)%alsz.gt.0) then
      rs = atmres(i)
      imol = molofrs(rs)
      skip = .false.
      do rsi=max(rs-1,rsmol(imol,1)),min(rs+1,rsmol(imol,2))
        if ((i.eq.wline(rsi)).OR.((i.eq.fline(rsi)).AND.(fline2(rsi).le.0)).OR.&
 &          ((i.eq.yline(rsi)).AND.(yline2(rsi).le.0))) then
          natlst%nr = natlst%nr + 1
          skip = .true.
          exit
        else if ((i.eq.fline(rsi)).OR.(i.eq.yline(rsi))) then
          skip=.true.
          exit
        end if
        do j=1,nchi(rsi)
          if (i.eq.chiline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit
        do j=1,nnucs(rsi)
          if (i.eq.nucsline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit 
      end do
      if (skip.EQV..false.) then
        if (seqtyp(atmres(i)).eq.26) then
          unklst%nr = unklst%nr + 1
        else
          unslst%nr = unslst%nr + 1
        end if
      end if
    end if
  end do
  if (natlst%nr.gt.0) then
    allocate(natlst%idx(natlst%nr))
    allocate(natlst%wt(natlst%nr))
  end if
  if (unklst%nr.gt.0) then
    allocate(unklst%idx(unklst%nr))
    allocate(unklst%wt(unklst%nr))
  end if
  if (unslst%nr.gt.0) then
    allocate(unslst%idx(unslst%nr))
    allocate(unslst%wt(unslst%nr))
  end if
!
  unklst%nr = 0
  unslst%nr = 0
  natlst%nr = 0
  do i=1,n
    if (izrot(i)%alsz.gt.0) then
      rs = atmres(i)
      imol = molofrs(rs)
      if (i.eq.atmol(imol,1)) cycle
      skip = .false.
      do rsi=max(rs-1,rsmol(imol,1)),min(rs+1,rsmol(imol,2))
        if ((i.eq.wline(rsi)).OR.((i.eq.fline(rsi)).AND.(fline2(rsi).le.0)).OR.&
 &          ((i.eq.yline(rsi)).AND.(yline2(rsi).le.0))) then
          natlst%nr = natlst%nr + 1
          natlst%idx(natlst%nr) = i
          skip = .true.
          exit
        else if ((i.eq.fline(rsi)).OR.(i.eq.yline(rsi))) then
          skip=.true.
          exit
        end if
        do j=1,nchi(rsi)
          if (i.eq.chiline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            natlst%idx(natlst%nr) = i
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit
        do j=1,nnucs(rsi)
          if (i.eq.nucsline(j,rsi)) then
            natlst%nr = natlst%nr + 1
            natlst%idx(natlst%nr) = i
            skip = .true.
            exit
          end if
        end do
        if (skip.EQV..true.) exit 
      end do
      if (skip.EQV..false.) then
        if (seqtyp(atmres(i)).eq.26) then
          unklst%nr = unklst%nr + 1
          unklst%idx(unklst%nr) = i
        else
          unslst%nr = unslst%nr + 1
          unslst%idx(unslst%nr) = i
        end if
      end if
    end if
  end do
!  write(*,*) 'UNS'
!  write(*,*)  unslst%idx(1:unslst%nr)
!  write(*,*) 'NAT'
!  write(*,*)  natlst%idx(1:natlst%nr)
!  write(*,*) 'UNK'
!  write(*,*)  unklst%idx(1:unklst%nr)
!
end
!
!------------------------------------------------------------------------------
!
! this routine parses the rotation set lists for a molecule to identify putative
! overlap relationships between them that are crucial for recursive computation
! of properties
!
subroutine parse_rotlsts(imol)
!
  use molecule
  use zmatrix
  use sequen
  use fyoc
  use atoms
  use interfaces
  use forces
  use iounit
!
  implicit none
!
  integer imol,i,k,l,nelg,kk,ll,rs,tmp(8),kid
  integer, ALLOCATABLE:: elgis(:,:),hlp(:,:),tmpl(:,:)
  logical fndit
!
  k = 0
  do i=atmol(imol,1)+3,atmol(imol,2)
    if (allocated(izrot(i)%rotis).EQV..true.) k = k + 1
  end do
!
  if (k.le.0) then
    izrot(atmol(imol,1))%alsz = 1
    if (allocated(izrot(atmol(imol,1))%rotis).EQV..true.) deallocate(izrot(atmol(imol,1))%rotis)
    allocate(izrot(atmol(imol,1))%rotis(1,2))
    izrot(atmol(imol,1))%rotis(1,1) = atmol(imol,1)
    izrot(atmol(imol,1))%rotis(1,2) = atmol(imol,2)
    if (allocated(izrot(atmol(imol,1))%diffis).EQV..true.) deallocate(izrot(atmol(imol,1))%diffis)
    izrot(atmol(imol,1))%alsz2 = 1
    allocate(izrot(atmol(imol,1))%diffis(izrot(atmol(imol,1))%alsz2,2))
    izrot(atmol(imol,1))%diffis(:,:) = izrot(atmol(imol,1))%rotis(:,:)
    if (allocated(izrot(atmol(imol,1))%treevs).EQV..true.) deallocate(izrot(atmol(imol,1))%treevs)
    allocate(izrot(atmol(imol,1))%treevs(5))
    izrot(atmol(imol,1))%treevs(:) = 0
    izrot(atmol(imol,1))%treevs(2) = 1
    izrot(atmol(imol,1))%treevs(4) = 1
    return
  end if
!
  allocate(elgis(k,9))
  elgis(:,:) = 0
  k = 0
  do i=atmol(imol,2),atmol(imol,1)+3,-1
    if (allocated(izrot(i)%rotis).EQV..true.) then
      k = k + 1
      elgis(k,1) = i
!     set total number of atoms for each group
      do kk=1,izrot(elgis(k,1))%alsz
        elgis(k,2) = elgis(k,2) + izrot(elgis(k,1))%rotis(kk,2) - izrot(elgis(k,1))%rotis(kk,1) + 1
      end do
    end if
  end do
!
! for each group, find the smallest larger group that completely contains it
  nelg = k
  do k=1,nelg
    do l=1,nelg
      if (k.eq.l) cycle
      if (elgis(k,2).gt.elgis(l,2)) cycle ! discard all smaller groups
      if (elgis(k,4).gt.0) then
        if (elgis(l,2).gt.elgis(elgis(k,4),2)) cycle ! discard all groups larger than the one currently found
      end if
      do kk=1,izrot(elgis(k,1))%alsz
        fndit = .false.
        do ll=1,izrot(elgis(l,1))%alsz
          if ((izrot(elgis(l,1))%rotis(ll,2).ge.izrot(elgis(k,1))%rotis(kk,2)).AND.&
 &            (izrot(elgis(l,1))%rotis(ll,1).le.izrot(elgis(k,1))%rotis(kk,1))) then
            fndit = .true.
            exit
          end if
        end do
        if (fndit.EQV..false.) exit
      end do
      if (fndit.EQV..false.) cycle
      elgis(k,4) = l
    end do
!    write(*,*) elgis(k,1),elgis(elgis(k,4),1)
    if (elgis(k,4).gt.0) elgis(elgis(k,4),3) = 1 ! flag the minimum container as a collection point
  end do
!
  l = 0
  do k=1,nelg
    if (elgis(k,3).gt.0) elgis(k,6) = -1
  end do
  do while (1.eq.1)
    fndit = .false.
    do k=1,nelg
      if (elgis(k,6).eq.l) then
        if (elgis(k,4).gt.0) then
          elgis(elgis(k,4),6) = l + 1 ! depth level for access
          fndit = .true.
        end if
      end if
    end do
    l = l + 1
    if (fndit.EQV..false.) exit
  end do
!
! sort according to depth
  do i=1,nelg
    elgis(i,5) = i
  end do
  kk = 1
  ll = nelg
  k = nelg
  fndit = .true.
  call merge_sort(k,fndit,elgis(:,6),elgis(:,7),kk,ll,elgis(:,5),elgis(:,8))
  do i=1,nelg
    elgis(elgis(i,8),5) = i
  end do
  elgis(:,6) = elgis(:,7)
  do i=1,nelg
    if (elgis(i,4).gt.0) elgis(i,4) = elgis(elgis(i,4),5)
  end do
  do i=1,nelg
    if (i.eq.elgis(i,8)) cycle
    tmp(1:4) = elgis(i,1:4)
    elgis(i,1:4) = elgis(elgis(i,8),1:4)
    elgis(elgis(i,8),1:4) = tmp(1:4)
    kk = elgis(i,5)
    ll = elgis(i,8)
    elgis(i,5) = elgis(ll,5)
    elgis(ll,5) = kk
    elgis(kk,8) = elgis(i,8)
  end do
!
! now index branches
  elgis(:,7:9) = 0
  kk = 1
  ll = 1
  do while (elgis(kk,6).eq.0)
    elgis(kk,7) = ll
    ll = ll + 1
    kk = kk + 1
    if (kk.gt.nelg) exit
  end do
  do k=1,nelg
    if (elgis(k,4).eq.0) cycle
    if (elgis(elgis(k,4),7).eq.0) then
      elgis(elgis(k,4),7) = elgis(k,7)
    else
!     this indicates branch merging
      if (elgis(elgis(k,4),7).lt.elgis(k,7)) then
        elgis(k,8) = elgis(elgis(k,4),7) ! the target container already has lower ID, end current branch
      else
        elgis(elgis(k,4),7) = elgis(k,7) ! the target container has higher ID, override and end relevant branches
        do i=1,k-1
          if (elgis(i,4).eq.elgis(k,4)) then
           elgis(i,8) = elgis(k,7)
          end if
        end do
      end if
    end if
  end do
!
! for all mergers, we need to choose the right place in hierarchy to perform merging
! first, populate a counter for how many branches are merging in at a given dihedral
  do k=1,nelg
    if (elgis(k,8).gt.0) then
      elgis(elgis(k,4),9) = elgis(elgis(k,4),9) + 1 ! elgis(k,7)
    end if
  end do  
!
! determine difference sets for rotation lists
  allocate(hlp(nelg,2))
  hlp(:,1) = 0
  do k=1,nelg
    if (elgis(k,4).gt.0) hlp(elgis(k,4),1) = hlp(elgis(k,4),1) + 1
  end do
  ll = maxval(hlp(:,1))
  deallocate(hlp)
  allocate(hlp(nelg,ll+1))
  hlp(:,1) = 0
  do k=1,nelg
    if (elgis(k,4).gt.0) then
      hlp(elgis(k,4),1) = hlp(elgis(k,4),1) + 1
      hlp(elgis(k,4),1+hlp(elgis(k,4),1)) = k
    end if
  end do
  allocate(tmpl(sum(izrot(atmol(imol,1):atmol(imol,2))%alsz),2))
  do k=1,nelg
    if (allocated(izrot(elgis(k,1))%diffis).EQV..true.) deallocate(izrot(elgis(k,1))%diffis)
    if (hlp(k,1).gt.0) then
      tmpl(1:izrot(elgis(k,1))%alsz,:) = izrot(elgis(k,1))%rotis(1:izrot(elgis(k,1))%alsz,:)
      izrot(elgis(k,1))%alsz2 = izrot(elgis(k,1))%alsz
      do rs=2,hlp(k,1)+1
        l = hlp(k,rs)
        do ll=1,izrot(elgis(l,1))%alsz
          kk = 1
          do while (kk.le.izrot(elgis(k,1))%alsz2)
!           note it can never happen that sets need to be joined because this would imply adding atoms
            if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
              izrot(elgis(k,1))%alsz2 = izrot(elgis(k,1))%alsz2 + 1
              tmpl(izrot(elgis(k,1))%alsz2,1) = izrot(elgis(l,1))%rotis(ll,2) + 1
              tmpl(izrot(elgis(k,1))%alsz2,2) = tmpl(kk,2)
              tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
              exit
            else if ((izrot(elgis(l,1))%rotis(ll,2).eq.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
              tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
              exit
            else if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).eq.tmpl(kk,1))) then
              tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
              exit
            else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,1))) then
              ! delete segment
              tmpl(kk:(izrot(elgis(k,1))%alsz2-1),:) = tmpl((kk+1):izrot(elgis(k,1))%alsz2,:)
              izrot(elgis(k,1))%alsz2 = izrot(elgis(k,1))%alsz2 - 1
            else
              kk = kk + 1
            end if
          end do
          do kk=1,izrot(elgis(k,1))%alsz2
            if ((izrot(elgis(l,1))%rotis(ll,2).gt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,2))) then
              tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
            else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,1)).AND.(izrot(elgis(l,1))%rotis(ll,1).lt.tmpl(kk,1))) then
              tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
            end if
          end do
        end do
      end do
      if (izrot(elgis(k,1))%alsz2.le.0) then
        write(ilog,*) 'Fatal. In computing difference lists of rotating atoms, list associated with atom ',elgis(k,1),' turns &
 &out empty. This is a bug.'
        call fexit()
      else
        allocate(izrot(elgis(k,1))%diffis(izrot(elgis(k,1))%alsz2,2))
        izrot(elgis(k,1))%diffis(:,:) = tmpl(1:izrot(elgis(k,1))%alsz2,:)
      end if
    else
      izrot(elgis(k,1))%alsz2 = 0 ! this is needed because a torsion is now terminal, but may not have been previously 
    end if
  end do
! sanity check
  deallocate(tmpl)
  deallocate(hlp)
  allocate(tmpl(atmol(imol,2)-atmol(imol,1)+1,2))
  tmpl(:,:) = 0
  do l=1,nelg
    if (elgis(l,6).eq.0) then
      do ll=1,izrot(elgis(l,1))%alsz
        do k=izrot(elgis(l,1))%rotis(ll,1),izrot(elgis(l,1))%rotis(ll,2)
          kk = k - atmol(imol,1) + 1
          if (tmpl(kk,1).eq.0) then
            tmpl(kk,1) = elgis(l,1)
          else
            write(ilog,*) 'Fatal. In setting up difference lists of rotating atoms, atom ',k,' is part of multiple lists.&
 & This is a bug.'
            call fexit()
          end if
        end do
      end do
    else
      do ll=1,izrot(elgis(l,1))%alsz2
        do k=izrot(elgis(l,1))%diffis(ll,1),izrot(elgis(l,1))%diffis(ll,2)
          kk = k - atmol(imol,1) + 1
          if (tmpl(kk,1).eq.0) then
            tmpl(kk,1) = elgis(l,1)
          else
            write(ilog,*) 'Fatal. In setting up difference lists of rotating atoms, atom ',k,' is part of multiple lists.&
 & This is a bug.'
            call fexit()
          end if
        end do
      end do
    end if
  end do
  deallocate(tmpl)
!
! allocate and transfer
  do k=1,nelg
     if (allocated(izrot(elgis(k,1))%treevs).EQV..true.) deallocate(izrot(elgis(k,1))%treevs)
     allocate(izrot(elgis(k,1))%treevs(5+elgis(k,9)))
  end do
  elgis(:,9) = 0
  do k=1,nelg
    if (elgis(k,4).gt.0) then
      izrot(elgis(k,1))%treevs(1) = elgis(elgis(k,4),1) ! minimal container
    else
      izrot(elgis(k,1))%treevs(1) = 0
    end if
    izrot(elgis(k,1))%treevs(2:3) = elgis(k,7:8) ! branch and merger indicator
    izrot(elgis(k,1))%treevs(3) = elgis(k,9)
    izrot(elgis(k,1))%treevs(4) = k ! rank
    izrot(elgis(k,1))%treevs(5) = elgis(k,6) ! depth
    if (elgis(k,8).gt.0) then
      elgis(elgis(k,4),9) = elgis(elgis(k,4),9) + 1
      izrot(elgis(elgis(k,4),1))%treevs(5+elgis(elgis(k,4),9)) = elgis(k,7)
    end if
  end do
!
! set up a rigid-body motion list pinned to the first atom
  izrot(atmol(imol,1))%alsz = 1
  if (allocated(izrot(atmol(imol,1))%rotis).EQV..true.) deallocate(izrot(atmol(imol,1))%rotis)
  allocate(izrot(atmol(imol,1))%rotis(1,2))
  izrot(atmol(imol,1))%rotis(1,1) = atmol(imol,1)
  izrot(atmol(imol,1))%rotis(1,2) = atmol(imol,2)
! difference list also
  allocate(tmpl(atmol(imol,2)-atmol(imol,1)+1,2))
  kid = 1
  tmpl(1,:) =  izrot(atmol(imol,1))%rotis(1,:)
  do l=1,nelg
    do ll=1,izrot(elgis(l,1))%alsz
!     subsets, supersets
      kk = 1
      do while (kk.le.kid)
        if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
          kid = kid + 1
          tmpl(kid,1) = izrot(elgis(l,1))%rotis(ll,2) + 1
          tmpl(kid,2) = tmpl(kk,2)
          tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
          exit
        else if ((izrot(elgis(l,1))%rotis(ll,2).eq.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).gt.tmpl(kk,1))) then
          tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
          exit
        else if ((izrot(elgis(l,1))%rotis(ll,2).lt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).eq.tmpl(kk,1))) then
          tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
          exit
        else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,1))) then
          ! delete segment
          tmpl(kk:(kid-1),:) = tmpl((kk+1):kid,:)
          kid = kid - 1
        else
          kk = kk + 1
        end if
      end do
!     overlaps
      do kk=1,kid
        if ((izrot(elgis(l,1))%rotis(ll,2).gt.tmpl(kk,2)).AND.(izrot(elgis(l,1))%rotis(ll,1).le.tmpl(kk,2))) then
          tmpl(kk,2) = izrot(elgis(l,1))%rotis(ll,1) - 1 ! existing segment shortened
        else if ((izrot(elgis(l,1))%rotis(ll,2).ge.tmpl(kk,1)).AND.(izrot(elgis(l,1))%rotis(ll,1).lt.tmpl(kk,1))) then
          tmpl(kk,1) = izrot(elgis(l,1))%rotis(ll,2) + 1 ! existing segment shortened
        end if
      end do
    end do
  end do
  if (allocated(izrot(atmol(imol,1))%diffis).EQV..true.) deallocate(izrot(atmol(imol,1))%diffis)
  izrot(atmol(imol,1))%alsz2 = kid
  allocate(izrot(atmol(imol,1))%diffis(kid,2))
  izrot(atmol(imol,1))%diffis(:,:) = tmpl(1:kid,:)
! finally, set treevs for base atom
  kk = 0
  do k=1,nelg
    if (elgis(k,4).le.0) kk = kk + 1
  end do
  if (allocated(izrot(atmol(imol,1))%treevs).EQV..true.) deallocate(izrot(atmol(imol,1))%treevs)
  allocate(izrot(atmol(imol,1))%treevs(5+kk))
  izrot(atmol(imol,1))%treevs(1) = 0
  izrot(atmol(imol,1))%treevs(2) = 1
  izrot(atmol(imol,1))%treevs(3) = kk
  izrot(atmol(imol,1))%treevs(4) = nelg + 1
  izrot(atmol(imol,1))%treevs(5) = maxval(elgis(1:nelg,6)) + 1
  kk = 0
  do k=1,nelg
    if (elgis(k,4).le.0) then
      kk = kk + 1
      izrot(atmol(imol,1))%treevs(5+kk) = elgis(k,7)
    end if
  end do
!
  deallocate(tmpl)
  deallocate(elgis)
!
end
!
!--------------------------------------------------------------------------------------------------------------
!
! a subroutine to determine the interaction status of two atoms via standard topology and rotation lists (including effective 14)
! this function has to deal with all specifics for both intra-residue and next-residue terms
!
subroutine ia_rotlsts(ii,kk,doia,are14)
!
  use atoms
  use sequen
  use molecule
  use zmatrix
  use inter
  use polypep
  use fyoc
  use system
!
  implicit none
!
  logical doia,are14
  integer j,jj,rs,rs1,rs2,ii,kk,shf,iij,kkj,rb1,rb2
  logical are14bu
!
  doia = .true.
  are14 = .false.
!
  if (molofrs(atmres(ii)).ne.(molofrs(atmres(kk)))) return
  if (abs(atmres(ii)-atmres(kk)).gt.1) return
!
  do j=1,n12(ii)
    if (i12(j,ii).eq.kk) doia = .false.
  end do
  do j=1,n13(ii)
    if (i13(j,ii).eq.kk) doia = .false.
  end do
  do j=1,n14(ii)
    if (use_14.EQV..false.) then
      if (i14(j,ii).eq.kk) doia = .false.
    else
      if (i14(j,ii).eq.kk) are14 = .true.
    end if
  end do
  are14bu = are14
!
! rotation list-based checks
  if (doia.EQV..true.) then
    iij = 0
    kkj = 0
    rs1 = min(atmres(ii),atmres(kk))
    rs2 = max(atmres(ii),atmres(kk))
    rs = rs1
    if (rs.gt.rsmol(molofrs(atmres(ii)),1)) rs = rs-1
    do j=at(rs)%bb(1),at(rs2)%bb(1)+at(rs2)%nbb+at(rs2)%nsc-1
      if (izrot(j)%alsz.le.0) cycle
      if (izrot(j)%treevs(5).eq.0) then
        do jj=1,izrot(j)%alsz
          if ((izrot(j)%rotis(jj,1).le.ii).AND.(izrot(j)%rotis(jj,2).ge.ii)) iij = j
          if ((izrot(j)%rotis(jj,1).le.kk).AND.(izrot(j)%rotis(jj,2).ge.kk)) kkj = j
        end do
      else
        do jj=1,izrot(j)%alsz2
          if ((izrot(j)%diffis(jj,1).le.ii).AND.(izrot(j)%diffis(jj,2).ge.ii)) iij = j
          if ((izrot(j)%diffis(jj,1).le.kk).AND.(izrot(j)%diffis(jj,2).ge.kk)) kkj = j
        end do
      end if
    end do
!
    if (iij.eq.kkj) then
      doia = .false.
    end if
!   need to check for corrections caused by axis atoms and set extended 14
    if (iij.gt.0) then
      if ((izrot(iij)%treevs(1).eq.kkj).AND.(mode_14.eq.2)) are14 = .true.
      if ((izrot(iij)%treevs(1).eq.0).AND.(kkj.eq.at(rs)%bb(1)).AND.(mode_14.eq.2)) are14 = .true.
      if ((iz(1,iij).eq.kk).OR.(iz(2,iij).eq.kk)) doia = .false.
      if (izrot(iij)%treevs(1).gt.0) then
        if (((iz(1,izrot(iij)%treevs(1)).eq.kk).OR.(iz(2,izrot(iij)%treevs(1)).eq.kk)).AND.(mode_14.eq.2)) are14 = .true.
      end if
    end if
    if (kkj.gt.0) then
      if ((izrot(kkj)%treevs(1).eq.0).AND.(iij.eq.at(rs)%bb(1)).AND.(mode_14.eq.2)) are14 = .true.
      if ((izrot(kkj)%treevs(1).eq.iij).AND.(mode_14.eq.2)) are14 = .true.
      if ((iz(1,kkj).eq.ii).OR.(iz(2,kkj).eq.ii)) doia = .false.
      if (izrot(kkj)%treevs(1).gt.0) then
        if (((iz(1,izrot(kkj)%treevs(1)).eq.ii).OR.(iz(2,izrot(kkj)%treevs(1)).eq.ii)).AND.(mode_14.eq.2)) are14 = .true.
      end if
    end if
!
!   we need an override if rotation lists are not appropriate to determine rigidity: this is the case
!   for flexible rings
    if ((seqpolty(rs1).eq.'N').AND.(nucsline(6,rs1).gt.0).AND.(at(rs1)%nsc.ge.3)) then
      rb1 = 0
      rb2 = 0 
      shf = 0
      if (seqflag(rs1).eq.24) shf = 2
      do j=1,n12(ii)
        if (((i12(j,ii).ge.at(rs1)%sc(1)).AND.(i12(j,ii).le.at(rs1)%sc(3))).OR.(i12(j,ii).eq.nuci(rs1,5-shf)).OR.&
 &          (i12(j,ii).eq.nuci(rs1,6-shf))) then
          rb1 = i12(j,ii)
          exit
        end if
      end do
      do j=1,n12(kk)
        if (((i12(j,kk).ge.at(rs1)%sc(1)).AND.(i12(j,kk).le.at(rs1)%sc(3))).OR.(i12(j,kk).eq.nuci(rs1,5-shf)).OR.&
 &          (i12(j,kk).eq.nuci(rs1,6-shf))) then
          rb2 = i12(j,kk)
          exit
        end if
      end do
      if ((rb1.gt.0).AND.(rb2.gt.0)) then
        doia = .true.
        are14 = are14bu
      end if
    else if ((seqpolty(rs1).eq.'P').OR.(seqpolty(rs2).eq.'P')) then
      shf = 0
      if (ua_model.ne.0) shf = 1
      do jj=rs1,rs2
        rb1 = 0
        rb2 = 0
        if ((at(jj)%nsc.lt.(4-shf)).OR.(seqpolty(jj).ne.'P')) cycle
        do j=1,n12(ii)
          if (((i12(j,ii).ge.at(jj)%sc(2-shf)).AND.(i12(j,ii).le.at(jj)%sc(4-shf))).OR.(i12(j,ii).eq.ni(jj)).OR.&
   &          (i12(j,ii).eq.cai(jj))) then
            rb1 = i12(j,ii)
            exit
          end if
        end do
        do j=1,n12(kk)
          if (((i12(j,kk).ge.at(jj)%sc(2-shf)).AND.(i12(j,kk).le.at(jj)%sc(4-shf))).OR.(i12(j,kk).eq.ni(jj)).OR.&
   &          (i12(j,kk).eq.cai(jj))) then
            rb2 = i12(j,kk)
            exit
          end if
        end do
        if ((rb1.gt.0).AND.(rb2.gt.0)) then
          doia = .true.
          are14 = are14bu
          exit
        end if
      end do
    end if
  end if
!
end
!
!--------------------------------------------------------------------------------------------------------------
!
! a subroutine to populate helper arrays for rotation lists
!
subroutine helper_rotlsts()
!
  use atoms
  use zmatrix
  use sequen
  use molecule
!
  implicit none
!
  integer ati,i,j,k
!
  do ati=1,n
    if (izrot(ati)%alsz.gt.0) then
      izrot(ati)%rsbnds(1,:) = nseq + 1
      izrot(ati)%rsbnds(2,:) = 0
      izrot(ati)%atmss(:) = 0.0
      do j=1,izrot(ati)%alsz
        izrot(ati)%rsbnds(1,1) = min(minval(atmres(izrot(ati)%rotis(j,1:2))),izrot(ati)%rsbnds(1,1))
        izrot(ati)%rsbnds(2,1) = max(maxval(atmres(izrot(ati)%rotis(j,1:2))),izrot(ati)%rsbnds(2,1))
        izrot(ati)%atmss(1) = izrot(ati)%atmss(1) + sum(mass(izrot(ati)%rotis(j,1):izrot(ati)%rotis(j,2)))
      end do
      if (ati.ne.atmol(molofrs(atmres(ati)),1)) then
        izrot(ati)%atmss(2) = molmass(moltypid(molofrs(atmres(ati)))) - izrot(ati)%atmss(1) - mass(iz(1,ati)) - mass(iz(2,ati))
      end if
      k = atmol(molofrs(atmres(ati)),1)
      i = 1
      do while (k.le.atmol(molofrs(atmres(ati)),2))
        if (i.le.izrot(ati)%alsz) then
          if ((k.ge.izrot(ati)%rotis(i,1)).AND.(k.le.izrot(ati)%rotis(i,2))) then
            k = izrot(ati)%rotis(i,2)+1
            i = i + 1
          else
            if ((k.ne.iz(1,ati)).AND.(k.ne.(iz(2,ati)))) then
              if (atmres(k).lt.izrot(ati)%rsbnds(1,2)) izrot(ati)%rsbnds(1,2) = atmres(k)
              if (atmres(k).gt.izrot(ati)%rsbnds(2,2)) izrot(ati)%rsbnds(2,2) = atmres(k)
            end if
            k = k + 1
          end if
        else
          if ((k.ne.iz(1,ati)).AND.(k.ne.(iz(2,ati)))) then
            if (atmres(k).lt.izrot(ati)%rsbnds(1,2)) izrot(ati)%rsbnds(1,2) = atmres(k)
            if (atmres(k).gt.izrot(ati)%rsbnds(2,2)) izrot(ati)%rsbnds(2,2) = atmres(k)
          end if
          k = k + 1
        end if
      end do
    end if
  end do
!
end
!
!-----------------------------------------------------------------------
!
! the following are a few subroutines/fxns to compute individual
! z-matrix terms (lengths, angles, dihedrals) from atomic xyzs
!
!-----------------------------------------------------------------------
!
!
! a function which gives the bond length
! returned as output for two atomic indices (i1,i2)
!
function getblen(i1,i2)
!
  use atoms
!
  implicit none
!
  integer i1,i2
  RTYPE c1(3),getblen
!
  c1(1) = x(i1) - x(i2)
  c1(2) = y(i1) - y(i2)
  c1(3) = z(i1) - z(i2)
  getblen = sqrt(sum(c1(:)**2.0))
!
end
!
!------------------------------------------------------------------
!
! the same operating on alternative coordinates
!
function getblen_ref(i1,i2)
!
  use atoms
!
  implicit none
!
  integer i1,i2
  RTYPE c1(3),getblen_ref
!
  c1(1) = xref(i1) - xref(i2)
  c1(2) = yref(i1) - yref(i2)
  c1(3) = zref(i1) - zref(i2)
  getblen_ref = sqrt(sum(c1(:)**2.0))
!
end
!
!-----------------------------------------------------------------------
!
! a function which gives bond angle mapped onto interval [0,180.0]
! returned as output for three atom indices (i1,i2,i3 in that order)
!
function getbang(i1,i2,i3)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,getbang,dotp,si
!
  c1(1) = x(i1) - x(i2)
  c1(2) = y(i1) - y(i2)
  c1(3) = z(i1) - z(i2)
  c2(1) = x(i3) - x(i2)
  c2(2) = y(i3) - y(i2)
  c2(3) = z(i3) - z(i2)
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  getbang = si*2.0*RADIAN*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    getbang = getbang + RADIAN*PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a function which gives bond angle mapped onto interval [0,180.0]
! returned as output for three atom indices (i1,i2,i3 in that order)
! differently from getbang, it explicitly checks for superimposed atoms and returns -1.0
! to indicate this condition
!
function getbang_safe(i1,i2,i3)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,getbang_safe,dotp,si
!
  c1(1) = x(i1) - x(i2)
  c1(2) = y(i1) - y(i2)
  c1(3) = z(i1) - z(i2)
  c2(1) = x(i3) - x(i2)
  c2(2) = y(i3) - y(i2)
  c2(3) = z(i3) - z(i2)
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  if ((cc1.le.0.0).OR.(cc2.le.0.0)) then
    getbang_safe = -1.0
    return
  end if
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  getbang_safe = si*2.0*RADIAN*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    getbang_safe = getbang_safe + RADIAN*PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same operating on alternative coordinates
!
function getbang_ref(i1,i2,i3)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,getbang_ref,dotp,si
!
  c1(1) = xref(i1) - xref(i2)
  c1(2) = yref(i1) - yref(i2)
  c1(3) = zref(i1) - zref(i2)
  c2(1) = xref(i3) - xref(i2)
  c2(2) = yref(i3) - yref(i2)
  c2(3) = zref(i3) - zref(i2)
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  getbang_ref = si*2.0*RADIAN*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    getbang_ref = getbang_ref + RADIAN*PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a function which gives torsional angle in degrees [-180:180] returned
! as output for four atom indices (i1,i2,i3,i4 in that order)
!
function getztor(i1,i2,i3,i4)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3,i4
  RTYPE getztor,c1(3),c2(3),c3(3),c4(3)
  RTYPE cp1(3),np1,cp2(3),np2,np4,si,dotp
!
  c1(1) = x(i2) - x(i1)
  c1(2) = y(i2) - y(i1)
  c1(3) = z(i2) - z(i1)
  c2(1) = x(i3) - x(i2)
  c2(2) = y(i3) - y(i2)
  c2(3) = z(i3) - z(i2)
  c3(1) = x(i4) - x(i3)
  c3(2) = y(i4) - y(i3)
  c3(3) = z(i4) - z(i3)
  call crossprod(c1,c2,cp1,np1)
  call crossprod(c2,c3,cp2,np2)
!
  if (np1*np2 .gt. 0.0d0) then
    dotp = sum(cp1(:)*cp2(:))
    if (dotp.lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    c4(:) = cp2(:)/np2 - si*cp1(:)/np1
    np4 = sqrt(sum(c4(:)**2.0))
    getztor = si*2.0*RADIAN*asin(0.5*np4)
    if (dotp.lt.0.0) then
      getztor = getztor + RADIAN*PI
    end if
! 
    si = sum(c1(:)*cp2(:))
    if (si.lt.0.0) getztor = -getztor
  else
!   torsion is ill-defined, left at 0.0 (not necessarily fatal, though)
    getztor = 0.0
  end if
!
end
!
!-------------------------------------------------------------------------
!
! the same as getztor but acting on xref, yref, zref
!
function getztor_ref(i1,i2,i3,i4)
!
  use atoms
  use math
!
  implicit none
!
  integer i1,i2,i3,i4
  RTYPE getztor_ref,c1(3),c2(3),c3(3),c4(3)
  RTYPE cp1(3),np1,cp2(3),np2,np4,si,dotp
!
  c1(1) = xref(i2) - xref(i1)
  c1(2) = yref(i2) - yref(i1)
  c1(3) = zref(i2) - zref(i1)
  c2(1) = xref(i3) - xref(i2)
  c2(2) = yref(i3) - yref(i2)
  c2(3) = zref(i3) - zref(i2)
  c3(1) = xref(i4) - xref(i3)
  c3(2) = yref(i4) - yref(i3)
  c3(3) = zref(i4) - zref(i3)
  call crossprod(c1,c2,cp1,np1)
  call crossprod(c2,c3,cp2,np2)
!
  if (np1*np2 .gt. 0.0d0) then
    dotp = sum(cp1(:)*cp2(:))
    if (dotp.lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    c4(:) = cp2(:)/np2 - si*cp1(:)/np1
    np4 = sqrt(sum(c4(:)**2.0))
    getztor_ref = si*2.0*RADIAN*asin(0.5*np4)
    if (dotp.lt.0.0) then
      getztor_ref = getztor_ref + RADIAN*PI
    end if
! 
    si = sum(c1(:)*cp2(:))
    if (si.lt.0.0) getztor_ref = -getztor_ref
  else
!   torsion is ill-defined, left at 0.0 (not necessarily fatal, though)
    getztor_ref = 0.0
  end if
!
end
!
!-----------------------------------------------------------------------
!
! the same as getztor while accounting for the possibility that atoms are 
! from different molecules and need to be image-shifted
! a function which gives torsional angle in degrees [-180:180] returned
! as output for four atom indices (i1,i2,i3,i4 in that order)
!
function getztor_inter(i1,i2,i3,i4)
!
  use atoms
  use math
  use sequen, ONLY: molofrs
!
  implicit none
!
  integer i1,i2,i3,i4
  RTYPE getztor_inter,c1(3),c2(3),c3(3),c4(3),svr(3),svr2(3)
  RTYPE cp1(3),np1,cp2(3),np2,np4,si,dotp
!
  svr(:) = 0.0
  if (molofrs(atmres(i1)).ne.molofrs(atmres(i2))) call dis_bound_rs(atmres(i1),atmres(i2),svr)
  c1(1) = x(i2) - x(i1) - svr(1)
  c1(2) = y(i2) - y(i1) - svr(2)
  c1(3) = z(i2) - z(i1) - svr(3)
  svr(:) = 0.0
  if (molofrs(atmres(i3)).ne.molofrs(atmres(i2))) call dis_bound_rs(atmres(i2),atmres(i3),svr)
  c2(1) = x(i3) - x(i2) + svr(1)
  c2(2) = y(i3) - y(i2) + svr(2)
  c2(3) = z(i3) - z(i2) + svr(3)
  svr2(:) = 0.0
  if (molofrs(atmres(i4)).ne.molofrs(atmres(i2))) call dis_bound_rs(atmres(i2),atmres(i4),svr2)
  c3(1) = x(i4) - x(i3) + svr2(1) -svr(1)
  c3(2) = y(i4) - y(i3) + svr2(2) -svr(2)
  c3(3) = z(i4) - z(i3) + svr2(3) -svr(3)
  call crossprod(c1,c2,cp1,np1)
  call crossprod(c2,c3,cp2,np2)
!
  if (np1*np2 .gt. 0.0d0) then
    dotp = sum(cp1(:)*cp2(:))
    if (dotp.lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    c4(:) = cp2(:)/np2 - si*cp1(:)/np1
    np4 = sqrt(sum(c4(:)**2.0))
    getztor_inter = si*2.0*RADIAN*asin(0.5*np4)
    if (dotp.lt.0.0) then
      getztor_inter = getztor_inter + RADIAN*PI
    end if
! 
    si = sum(c1(:)*cp2(:))
    if (si.lt.0.0) getztor_inter = -getztor_inter
  else
!   torsion is ill-defined, left at 0.0 (not necessarily fatal, though)
    getztor_inter = 0.0
  end if
!
end
!
!-------------------------------------------------------------------------
!
! a simple wrapper to substitute for a (missing) pointer array
!
function getpuckertor(rs,zline)
!
  use polypep
  use sequen
  use molecule
  use iounit
  use aminos
  use fyoc
!
  implicit none
!
  integer rs,zline,shf
  RTYPE getztor,getpuckertor
!
  if (seqpolty(rs).eq.'N') then
    shf = 0
    if (((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)).OR.((seqtyp(rs).eq.26).AND.(nuci(rs,5).le.0))) shf = 2
    getpuckertor = getztor(nuci(rs,4-shf),nuci(rs,5-shf),nuci(rs,6-shf),nucsline(6,rs))
    zline = nucsline(6,rs)
  else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.25).OR.(seqtyp(rs).eq.32)) then
    if (rs.eq.rsmol(molofrs(rs),1)) then
      getpuckertor = getztor(ci(rs),cai(rs),ni(rs),fline(rs))
      zline = fline(rs)
    else
      getpuckertor = getztor(ci(rs),cai(rs),ni(rs),ci(rs-1))
      zline = ci(rs)
    end if
  else
     write(ilog,*) 'Fatal. Called getpuckertor(...) with unsupported residue.'
     call fexit()
  end if
!
  return
!
end
!
!-----------------------------------------------------------------------
!
! gives bond angle in radian
! angle is always mapped onto interval [0,PI/2] (smaller inside angle only)
!
subroutine bondang(p1,p2,p3,ang)
!
  implicit none
!
  integer k
  RTYPE p1(3),p2(3),p3(3),c1(3),c2(3),c3(3),cc1,cc2,cc3,dotp,ang,si
!
  cc1 = 0.0
  cc2 = 0.0
  do k=1,3
    c1(k) = p1(k) - p2(k)
    cc1 = cc1 + c1(k)**2
    c2(k) = p2(k) - p3(k)
    cc2 = cc2 + c2(k)**2
  end do
  cc1 = sqrt(cc1)
  cc2 = sqrt(cc2)
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  ang = si*2.0*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    ang = -ang
  end if
!
end
!
!------------------------------------------------------------
!
! gives bond angle in radian, assumes atoms p1,p2,p3 connected
! in this order
! angle is always mapped onto interval [0,PI] (inside angles only)
!
subroutine bondang2(p1,p2,p3,ang)
!
  use math
!
  implicit none
!
  integer k
  RTYPE p1(3),p2(3),p3(3),c1(3),c2(3),c3(3),cc1,cc2,cc3,si,dotp,ang
!
  cc1 = 0.0
  cc2 = 0.0
  do k=1,3
    c1(k) = p1(k) - p2(k)
    cc1 = cc1 + c1(k)**2
    c2(k) = p3(k) - p2(k)
    cc2 = cc2 + c2(k)**2
  end do
  cc1 = sqrt(cc1)
  cc2 = sqrt(cc2)
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  ang = si*2.0*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    ang = ang + PI
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a similar subroutine which also gives the derivatives with respect to p1,p2,p3
!
subroutine bondang2_dcart(p1,p2,p3,ang,fvec1,fvec2,fvec3)
!
  use iounit
  use polypep
  use atoms
  use math
!
  implicit none
!
  RTYPE dv12(3),id12,p1(3),p2(3),p3(3),ang
  RTYPE dv32(3),id32,dotp
  RTYPE vecs(3),nvecs,frang
  RTYPE fvec1(3),fvec3(3),fvec2(3),si
!
  dv12(:) = p1(:) - p2(:)
  dv32(:) = p3(:) - p2(:)
  id12 = 1.0/sqrt(dv12(1)**2 + dv12(2)**2 + dv12(3)**2)
  id32 = 1.0/sqrt(dv32(1)**2 + dv32(2)**2 + dv32(3)**2)
  dotp = dv12(1)*dv32(1) + dv12(2)*dv32(2) + &
 &          dv12(3)*dv32(3)
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
! this half-angle formula has the advantage of being safe to use in all circumstances
! as 0.5*nvecs never gets near the -1/1 limits
! in addition there is considerable precision loss with the (simpler) acos-form
  vecs(1) = id32*dv32(1) - si*id12*dv12(1)
  vecs(2) = id32*dv32(2) - si*id12*dv12(2)
  vecs(3) = id32*dv32(3) - si*id12*dv12(3)
  nvecs = sqrt(vecs(1)**2 + vecs(2)**2 + vecs(3)**2)
  ang = si*2.0*asin(0.5*nvecs)
  frang = si*2.0/sqrt(1.0 - (0.5*nvecs)**2)
  if (dotp.lt.0.0) then
    ang = ang + PI
  end if
! there is a price, however, primarily in the derivatives
  fvec1(1) = 0.5*frang*(1.0/nvecs)*(vecs(1)*&
 &  (-si*id12 + si*dv12(1)*(id12**3)*dv12(1))&
 & + vecs(2)*(si*dv12(2)*(id12**3)*dv12(1))&
 & + vecs(3)*(si*dv12(3)*(id12**3)*dv12(1)) )
  fvec1(2) = 0.5*frang*(1.0/nvecs)*(vecs(2)*&
 &  (-si*id12 + si*dv12(2)*(id12**3)*dv12(2))&
 & + vecs(1)*(si*dv12(1)*(id12**3)*dv12(2))&
 & + vecs(3)*(si*dv12(3)*(id12**3)*dv12(2)) )
  fvec1(3) = 0.5*frang*(1.0/nvecs)*(vecs(3)*&
 &  (-si*id12 + si*dv12(3)*(id12**3)*dv12(3))&
 & + vecs(1)*(si*dv12(1)*(id12**3)*dv12(3))&
 & + vecs(2)*(si*dv12(2)*(id12**3)*dv12(3)) )
  fvec3(1) = 0.5*frang*(1.0/nvecs)*(vecs(1)*&
 &  (id32 - dv32(1)*(id32**3)*dv32(1))&
 & + vecs(2)*(-dv32(2)*(id32**3)*dv32(1))&
 & + vecs(3)*(-dv32(3)*(id32**3)*dv32(1)) )
  fvec3(2) = 0.5*frang*(1.0/nvecs)*(vecs(2)*&
 &  (id32 - dv32(2)*(id32**3)*dv32(2))&
 & + vecs(1)*(-dv32(1)*(id32**3)*dv32(2))&
 & + vecs(3)*(-dv32(3)*(id32**3)*dv32(2)) )
  fvec3(3) = 0.5*frang*(1.0/nvecs)*(vecs(3)*&
 &  (id32 - dv32(3)*(id32**3)*dv32(3))&
 & + vecs(1)*(-dv32(1)*(id32**3)*dv32(3))&
 & + vecs(2)*(-dv32(2)*(id32**3)*dv32(3)) )
  fvec2(:) = - (fvec1(:) + fvec3(:))
!
end
!
!------------------------------------------------------------
!
! gives angle between bond two vectors
! angle is always mapped onto interval [0,PI] (inside angle only)
!
subroutine bondang3(c1,c2,ang)
!
  use math
!
  implicit none
!
  RTYPE c1(3),c2(3),c3(3),cc1,cc2,cc3,dotp,ang,si
!
  cc1 = sqrt(sum(c1(:)**2.0))
  cc2 = sqrt(sum(c2(:)**2.0))
  dotp = sum(c1(:)*c2(:))
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c3(:) = c2(:)/cc2 - si*c1(:)/cc1
  cc3 = sqrt(sum(c3(:)**2.0))
!
  ang = si*2.0*asin(0.5*cc3)
  if (dotp.lt.0.0) then
    ang = ang + PI
  end if
!
end
!
!----------------------------------------------------------------------------
!
! gives torsion in radian, note sign convention
!
subroutine dihed(p1,p2,p3,p4,tor)
!
  use iounit
  use math
  use zmatrix
!
  implicit none
!
  integer k
  RTYPE p1(3),p2(3),p3(3),p4(3),c1(3),c2(3),c3(3),c4(3),np4
  RTYPE pl_1(3),pl_2(3),np1,np2,si,sign_conv,tor,dotp
!
  tor = 0.0
!
  do k=1,3
    c1(k) = p2(k) - p1(k)
    c2(k) = p3(k) - p2(k)
    c3(k) = p4(k) - p3(k)
  end do
  call crossprod(c1,c2,pl_1,np1)
  call crossprod(c2,c3,pl_2,np2)
!
  if (np1*np2 .gt. 0.0d0) then
    dotp = sum(pl_1(:)*pl_2(:))
    if (dotp.lt.0.0) then
      si = -1.0
    else
      si = 1.0
    end if
    c4(:) = pl_2(:)/np2 - si*pl_1(:)/np1
    np4 = sqrt(sum(c4(:)**2.0))
    tor = si*2.0*asin(0.5*np4)
    if (dotp.lt.0.0) then
      tor = tor + PI
    end if
! 
    sign_conv = c1(1)*pl_2(1) + c1(2)*pl_2(2) + c1(3)*pl_2(3)
    if (sign_conv .lt. 0.0d0)  tor = -tor
  else
!   torsion is ill-defined, left at 0.0 (not necessarily fatal, though)
    dihed_wrncnt = dihed_wrncnt + 1
    if (dihed_wrncnt.eq.dihed_wrnlmt) then
      write(ilog,*) 'Warning: Call to dihed(...) with bad vectors.'
      write(ilog,*) 'This is warning number #',dihed_wrncnt,' of this type not all of&
 & which may be displayed.'
      if (10.0*dihed_wrnlmt.gt.0.5*HUGE(dihed_wrnlmt)) then
        dihed_wrncnt = 0 ! reset
      else
        dihed_wrnlmt = dihed_wrnlmt*10
      end if
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a similar subroutine which lso gives derivatives with respect to
! p1-p4
!
subroutine dihed_dcart(p1,p2,p3,p4,tor,dtmd1,dtmd2,dtmd3,dtmd4)
!
  use iounit
  use polypep
  use atoms
  use math
!
  implicit none
!
  integer j
  RTYPE dv12(3),dv23(3),p1(3),p2(3),p3(3),p4(3)
  RTYPE dv43(3),ndv23,incp123,incp234,ncp123,ncp234
  RTYPE dotp,dotpfg,cp123(3),cp234(3),dtmd1(3),dtmd2(3)
  RTYPE dtmd3(3),dotphg,dtmd4(3),dotpah,si
  RTYPE tor,c4(3),np4
!
  dv12(:) = p1(:) - p2(:) 
  dv23(:) = p2(:) - p3(:)
  dv43(:) = p4(:) - p3(:)

  cp123(1) = dv12(2)*dv23(3) - dv12(3)*dv23(2)
  cp123(2) = dv12(3)*dv23(1) - dv12(1)*dv23(3)
  cp123(3) = dv12(1)*dv23(2) - dv12(2)*dv23(1)
  cp234(1) = dv43(2)*dv23(3) - dv43(3)*dv23(2)
  cp234(2) = dv43(3)*dv23(1) - dv43(1)*dv23(3)
  cp234(3) = dv43(1)*dv23(2) - dv43(2)*dv23(1)
! norms and dot-product
  ncp123 = cp123(1)**2 + cp123(2)**2 + cp123(3)**2
  ncp234 = cp234(1)**2 + cp234(2)**2 + cp234(3)**2
  dotp = cp123(1)*cp234(1) + cp123(2)*cp234(2) + &
 &          cp123(3)*cp234(3)
  dotpfg = dv12(1)*dv23(1) + dv12(2)*dv23(2) +&
 &           dv12(3)*dv23(3)
  dotphg = dv43(1)*dv23(1) + dv43(2)*dv23(2) +&
 &           dv43(3)*dv23(3)
  dotpah = cp123(1)*dv43(1) + cp123(2)*dv43(2) +&
 &           cp123(3)*dv43(3)
  ndv23 = sqrt(dv23(1)**2 + dv23(2)**2 + dv23(3)**2)
!
! filter for exceptions which would otherwise NaN
! colinear reference atoms (note that this also detects the equally fatal |dv23| = 0)
  if ((ncp123.le.0.0).OR.(ncp234.le.0.0)) then
!   this will set all derivatives to zero and the (now arbitrary) dihedral to zero deg.
    write(ilog,*) 'WARNING. Colinear reference atoms in computation of &
 &torsional derivative (dihed_dcart(...)). This indicates an unstabl&
 &e simulation or a bug. Expect run to crash soon.'
    incp123 = 0.0
    incp234 = 0.0
  else
    incp123 = 1.0/ncp123
    incp234 = 1.0/ncp234
  end if
!
! equations 27 in Blondel and Karplus
  do j=1,3
    dtmd1(j) = -cp123(j)*incp123*ndv23

    dtmd2(j) = cp123(j)*incp123*&
 &                  (ndv23 + dotpfg/ndv23) - &
 &               cp234(j)*dotphg*incp234/ndv23

    dtmd3(j) = cp234(j)*incp234*&
 &                  (dotphg/ndv23 - ndv23) - &
 &               cp123(j)*dotpfg*incp123/ndv23
    dtmd4(j) = cp234(j)*incp234*ndv23
  end do
!
  if (dotp.lt.0.0) then
    si = -1.0
  else
    si = 1.0
  end if
  c4(:) = cp234(:)*sqrt(incp234) - si*cp123(:)*sqrt(incp123)
  np4 = sqrt(sum(c4(:)**2.0))
  tor = si*2.0*asin(0.5*np4)
  if (dotp.lt.0.0) then
    tor = tor + PI
  end if
! 
  si = sum(-dv12(:)*cp234(:))
  if (si.lt.0.0) tor = -tor
!
end
!
!-----------------------------------------------------------------------
!
! a routine to check for colinear reference atoms for torsional tuples
!
function check_colinear(i1,i2,i3,i4)
!
  implicit none
!
  integer i1,i2,i3,i4
  logical check_colinear
  RTYPE ang,getbang
!
  check_colinear = .false.
!
  ang = getbang(i1,i2,i3)
  if ((abs(ang).le.0.01).OR.(abs(180.0-ang).le.0.01)) then
    check_colinear = .true.
    return
  end if
!
  ang = getbang(i2,i3,i4)
  if ((abs(ang).le.0.01).OR.(abs(180.0-ang).le.0.01)) check_colinear = .true.
!
  return
!
end
!
!--------------------------------------------------------------------------
!
! the same for the reference coordinates
!
function check_colinear_ref(i1,i2,i3,i4)
!
  implicit none
!
  integer i1,i2,i3,i4
  logical check_colinear_ref
  RTYPE ang,getbang_ref
!
  check_colinear_ref = .false.
!
  ang = getbang_ref(i1,i2,i3)
  if ((abs(ang).le.0.01).OR.(abs(180.0-ang).le.0.01)) then
    check_colinear_ref = .true.
    return
  end if
!
  ang = getbang_ref(i2,i3,i4)
  if ((abs(ang).le.0.01).OR.(abs(180.0-ang).le.0.01)) check_colinear_ref = .true.
!
  return
!
end
!
!--------------------------------------------------------------------------
!

