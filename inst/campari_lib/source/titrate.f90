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
! MAIN AUTHOR:   Hoang Tran                                                !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! attempts to (de)ionize residues
!
!
! pKa values from Li, Robertson, et al. (2005) "Very fast empirical
! prediction and rationalization of protein pKa values"
! Proteins 61(4):704-21
!
! -- assigns atom numbers for titration sites
!
subroutine ionize_setup()
!
  use iounit
  use aminos
  use system
  use molecule
  use ionize
  use polypep
  use sequen
  use movesets
  use mcsums
!
  implicit none
!
  integer imol,rs,aone,shf
  character(3) resname

!cc      integer i
!
  aone = 1
  nionres = 0
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
! do a dry run first
!
  do imol=1,nmol
     do rs=rsmol(imol,1),rsmol(imol,2)
        resname = amino(seqtyp(rs))
        if (resname.eq.'ASP') then
           nionres = nionres + 1
        else if (resname.eq.'GLU') then
           nionres = nionres + 1
        else if (resname.eq.'LYS') then
           nionres = nionres + 1
        else if (resname.eq.'ARG') then
           nionres = nionres + 1
        else if ((resname.eq.'HIE').or.(resname.eq.'HID')) then
           nionres = nionres + 1
        end if
     end do
  end do
!
  if (nionres.eq.0) then
    write(ilog,*) 'WARNING: No ionizable residues in sequence. Turni&
 &ng off titration moves.'
    write(ilog,*)
    phfreq = 0.0
    have_ph = .false.
    phout = nsim + 1
    return
  end if
!
  call allocate_ionize(aone)
!
  nionres = 0

  do imol=1,nmol
     do rs=rsmol(imol,1),rsmol(imol,2)
        resname = amino(seqtyp(rs))
        if (resname.eq.'ASP') then
           nionres = nionres + 1
           ionres(nionres) = rs
           ionstate(nionres) = -1
           pkai(nionres) = 3.8d0
           titsite(nionres) = at(rs)%sc(3-shf)
!           write(*,*)'ASP titsite type = ',b_type(titsite(nionres))
        else if (resname.eq.'GLU') then
           nionres = nionres + 1
           ionres(nionres) = rs
           ionstate(nionres) = -1
           pkai(nionres) = 4.5d0
           titsite(nionres) = at(rs)%sc(4-shf)
!           write(*,*)'GLU titsite type = ',b_type(titsite(nionres))
        else if (resname.eq.'LYS') then
           nionres = nionres + 1
           ionres(nionres) = rs
           ionstate(nionres) = 1
           pkai(nionres) = 10.5d0
           titsite(nionres) = at(rs)%sc(6-shf)
!           write(*,*)'LYS titsite type = ',b_type(titsite(nionres))
        else if (resname.eq.'ARG') then
           nionres = nionres + 1
           ionres(nionres) = rs
           ionstate(nionres) = 1
           pkai(nionres) = 12.5d0
           titsite(nionres) = at(rs)%sc(6-shf)
!           write(*,*)'ARG titsite type = ',b_type(titsite(nionres))
        else if ((resname.eq.'HIE').or.(resname.eq.'HID')) then
           nionres = nionres + 1
           ionres(nionres) = rs
           ionstate(nionres) = 0
           pkai(nionres) = 6.5d0
           titsite(nionres) = at(rs)%sc(7-shf)
!           write(*,*)'HIS titsite type = ',b_type(titsite(nionres))
        end if
     end do
  end do
!
!
end subroutine
!
!----------------------------------------------------------------------------
!
! (de)protonation move
! attempt to change protonation state of residues
!
subroutine protmov()
!
  use atoms
  use aminos
  use ionize
  use iounit
  use mcsums
  use system
  use sequen
  use units
  use accept
!
  implicit none
!
  integer ir,oldst,newst,i
  RTYPE random,dgprotint,diff,expterm,rcheck
  RTYPE xi,yi,zi,xj,yj,zj,dist2,dist
  RTYPE eep,dielec,dhfac
  character(3) resname

!
! dielectric constant.  water = 78.5
  dielec = 78.5d0
  newst = 0
!
! choose an ionizable residue randomly
  ir = int(random()*nionres) + 1
  diff = 0.0d0
  resname = amino(seqtyp(ionres(ir)))
  oldst = ionstate(ir)
!
! attempt to change protonation state

!
! "intrinsic" deltaG of protonation
  dgprotint = 2.303d0*(1.0d0/invtemp)*(phs - pkai(ir))
  if ((resname.eq.'ASP') .or. (resname.eq.'GLU')) then
     if (oldst .eq. -1) then
        diff = diff + dgprotint
        newst = 0
     else if (oldst .eq. 0) then
        diff = diff - dgprotint
        newst = -1
     else
        write(ilog,*)'Unknown protonation state'
        call fexit()
     end if
  else if ((resname.eq.'LYS') .or. (resname.eq.'ARG') .or.&
 &        (resname.eq.'HIE') .or. (resname.eq.'HID')) then
     if (oldst .eq. 0) then
        diff = diff + dgprotint
        newst = 1
     else if (oldst .eq. 1) then
        diff = diff - dgprotint
        newst = 0
     else
        write(ilog,*)'Unknown protonation state'
        call fexit()
     end if
  else
     write(ilog,*)'Unsupported ionizable residue: ',resname
     call fexit()
  end if
!
! electrostatic interactions with other titratable groups
  eep = 0.0d0
  xi = x(titsite(ir))
  yi = y(titsite(ir))
  zi = z(titsite(ir))
!  write(*,'(i3,"xi,yi,zi:",3f10.3)'),ionres(ir),xi,yi,zi
!  write(*,*)'---ionstate,newst---',ionstate(ir),newst
  do i=1,nionres
     if (i .ne. ir) then
        xj = x(titsite(i))
        yj = y(titsite(i))
        zj = z(titsite(i))
        dist2 = (xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2
        dist = sqrt(dist2)
! energy for placing a unit charge at site ir
        dhfac = exp(-invdhlen*dist)
        eep = eep + electric*ionstate(i)*dhfac/(dielec*dist)
!        write(*,*)'---dionst,energy--',ionstate(i),
! &           electric*ionstate(i)/(dielec*dist)
!        write(*,'(i3,"xj,yj,zj,dist to i:",4f10.3)')ionres(i),
! &           xj,yj,zj,dist
     end if
  end do
  if ((oldst .eq. 0) .and. (newst .eq. 1)) then
!   do nothing
  else if ((oldst .eq. 1) .and. (newst .eq. 0)) then
     eep = -eep
  else if ((oldst .eq. 0) .and. (newst .eq. -1)) then
     eep = -eep
  else if ((oldst .eq. -1) .and. (newst .eq. 0)) then
!   do nothing
  else
     write(ilog,*)' Unknown proton transition:'
     call fexit()
  end if
      
  diff = diff + eep

!  write(*,*)'-----',diff,eep,diff - eep

  if (diff .le. 0.0d0) then
     ionstate(ir) = newst
     acc%nph = acc%nph + 1
  else
     expterm = exp(-diff*invtemp)
     rcheck = random()
     if (expterm .gt. rcheck) then
        ionstate(ir) = newst
        acc%nph = acc%nph + 1
     end if
  end if
!
end subroutine
!
!---------------------------------------------------------------------------
!
! prints out the current ionization state
!
subroutine prt_phstate()
!
  use ionize
  use mcsums
  use molecule
  use sequen
  use grandensembles
  use system
!
  implicit none
!
  integer i,imol
  logical ismember
!
  prtionst(:) = ionstate(:)
  do i=1,nionres
    imol = molofrs(ionres(i))
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) prtionst(i) = ionstate(i) + 10
      end if
    end if
  end do
!
  write(ipht,10)(prtionst(i),i=1,nionres)
 10   format(2000000(i10))
!
end subroutine
!
!-----------------------------------------------------------------------------
!
! setup
!
! -- prints out the names of ionizable residues
! -- calculates Debye-Huckel screening length
!
! D-H screening term:  exp(-k r), k is invdhlen
! k^2 = [(8 pi Na ec^2 ) / (1000 ew kb T)] I
!
subroutine ph_setup()
!
  use aminos
  use ionize
  use math
  use mcsums
  use system
  use sequen
  use units
!
  implicit none
!
  integer ir
  RTYPE elemc,dielw,boltzcgs,bscr
! elementary charge, in statcoulombs, or esu's
  parameter (elemc=4.8032068d-10)
  parameter (boltzcgs=1.3803d-16)

! water dielectric
  dielw = 78.5

  write(ipht,20,advance="no")(ionres(ir),amino(seqtyp(ionres(ir))),ir=1,nionres)
 20   format(2000000(i6,1x,a3))
  write(ipht,*)''

! compare with table for parameter B on p211 of Modern Electrochemistry 1, Bockris and Reddy
  bscr = 8*pi*AVOGADRO*elemc*elemc/(1000*dielw*boltzcgs*kelvin)
  bscr = sqrt(bscr * 1d-16)
!  write(*,*)'-------bscr = ' ,bscr
! bscr is ~1/3.04 1/Angst. at 298 K
  invdhlen = bscr*sqrt(ionicstr)
!  write(*,*)'-------invdhlen = ' ,invdhlen
end
!
