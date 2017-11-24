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
!-----------------------------------------------------------------
!
! fill some fmsmc-specific atom arrays
!
subroutine absinth_atom()
!
  use atoms
  use params
  use sequen
  use polypep
  use energies
  use aminos
  use zmatrix
  use math
  use system
  use iounit
  use fyoc
!
  implicit none
!
  integer i,j,k,rs,ati,shf,shf2,shf3,atwo,afour
  logical sayyes
  RTYPE dum,dis,d1,d2
  logical saidit,saidit2
!
  atwo = 2
  afour = 4
  sayyes = .true.
  saidit = .false.
 55 format('Here, atom ',i7,' of type ',i4,' has radius ',g10.3,'.')
 56 format('Warning. Atom ',i7,' of biotype ',i5,' and LJ-type ',i5,' has ill-defined mass. Further warnings omitted.')
!
! default radius assignment
  atr(1:n) = lj_rad(attyp(1:n))
! potential patch overrides for radii
  call read_atpatchfile(afour)
!
! now checks
  do i=1,n
    if ((atr(i).le.0.0).AND.((use_IMPSOLV.EQV..true.).OR.(savcalc.le.nsim))) then
      write(ilog,*) 'Fatal. Atomic radii have to be finite for ABSINTH implicit solvent model&
 & or solvent-accessible volume analysis to work.'
      write(ilog,55) i,attyp(i),atr(i)
      write(ilog,*) 'Use the RADIUS keyword in the parameter file or a patch (FMCSC_RPATCHFILE) to solve this problem.'
      call fexit()
    end if
    if (mass(i).gt.0.0) then
      ivms(i) = 1.0/mass(i)
    else
      ivms(i) = 0.0
    end if
    atvol(i) = (4./3.)*PI*atr(i)**3
    atbvol(i) = (4./3.)*PI*(atr(i)+par_IMPSOLV(1))**3 - atvol(i)
    if ((mass(i).le.0.0).AND.(saidit.EQV..false.)) then
      saidit = .true.
      write(ilog,56) i,b_type(i),bio_ljtyp(b_type(i))
    end if
  end do
  do i=1,n
    atsavred(i) = atvol(i)
  end do
  do i=1,n
    k = iz(1,i)
!   disconnected atoms have to be excluded
!   note that even though the first atom in a multi-atom molecule
!   will always have k=0, the connectivity will still be handled
!   through the atom(s) actually connected to the first one
    if (k.le.0) cycle
    dis = sqrt((x(k) - x(i))**2&
 &           + (y(k) - y(i))**2&
 &           + (z(k) - z(i))**2)
    if (dis.lt.(atr(i)+atr(k))) then
      if (atr(i).ge.atr(k)) then
        if (dis.lt.(atr(i)-atr(k))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(k)
          atsavred(k) = atsavred(k) - 0.5*atvol(k)
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
          atsavred(k) = atsavred(k) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
        end if
      else
        if (dis.lt.(atr(k)-atr(i))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(i)
          atsavred(k) = atsavred(k) - 0.5*atvol(i)
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
          atsavred(k) = atsavred(k) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
        end if
      end if
    end if
  end do
!
! rings have additional bonds not showing up through iz - recover through nadd/iadd
  shf = 0
  shf2 = 0
  shf3 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
  end if
  do rs=1,nadd
    i = iadd(1,rs)
    k = iadd(2,rs)
    dis = sqrt((x(k) - x(i))**2&
 &           + (y(k) - y(i))**2&
 &           + (z(k) - z(i))**2)
!   distance override for crosslinks may be necessary
    do j=1,n_crosslinks
      if ((crosslink(j)%itstype.eq.1).OR.(crosslink(j)%itstype.eq.2)) then
        if (((i.eq.at(crosslink(j)%rsnrs(1))%sc(3-shf)).AND.(k.eq.at(crosslink(j)%rsnrs(2))%sc(3-shf))).OR.&
 &          ((k.eq.at(crosslink(j)%rsnrs(1))%sc(3-shf)).AND.(i.eq.at(crosslink(j)%rsnrs(2))%sc(3-shf)))) then
          dis = 2.03 ! make consistent with elsewhere
!              -> note we cannot infer the distance at this point since structure does not satisfy crosslink yet
          exit
        end if
      else
        write(ilog,*) 'Fatal. Encountered unsupported crosslink type in absinth_atom().&
 & This is most likely an omission bug.'
        call fexit()
      end if
    end do
    if (dis.lt.(atr(i)+atr(k))) then
      if (atr(i).ge.atr(k)) then
        if (dis.lt.(atr(i)-atr(k))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(k)
          atsavred(k) = atsavred(k) - 0.5*atvol(k)
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
          atsavred(k) = atsavred(k) - 0.5*atvol(k)*&
 &                       ((atr(i)+atr(k)-dis)/(2.0*atr(k)))
        end if
      else
        if (dis.lt.(atr(k)-atr(i))) then
          atsavred(i) = atsavred(i) - 0.5*atvol(i)
          atsavred(k) = atsavred(k) - 0.5*atvol(i)
        else
          atsavred(i) = atsavred(i) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
          atsavred(k) = atsavred(k) - 0.5*atvol(i)*&
 &                       ((atr(k)+atr(i)-dis)/(2.0*atr(i)))
        end if
      end if
    end if
  end do
!
  do i=1,n
    atsavred(i) = atsavred(i)/atvol(i)
  end do
  call read_atpatchfile(atwo)
  call setup_savol(sayyes,dum)
!
  if ((use_IMPSOLV.EQV..true.).OR.(savcalc.le.nsim)) then
    saidit = .false.
    saidit2 = .false.
    do ati=1,n
      if (atsavmaxfr(ati).gt.1.0) then
         write(ilog,*) 'Fatal. Setup of solvation parameters is incorrect for&
 & atom ',ati,' (type: ',b_type(ati),'). This is a bug.'
         call fexit()
      else if (atsavmaxfr(ati).le.0.0) then
         if (saidit.EQV..false.) then
           saidit = .true.
           write(ilog,*) 'Fatal. Maximum solvent-accessible volume fraction for&
 & atom ',ati,' (type: ',b_type(ati),') is zero. This could be a bug or indicates &
 &radii and bond lengths that are incompatible with underlying assumptions. Check &
 &FMCSC_SAVPATCHFILE, FMCSC_ASRPATCHFILE, and FMCSC_RPATCHFILE for workarounds or use FMCSC_UNSAFE to skip this error.&
 & Further warnings omitted.'
         end if
         if (be_unsafe.EQV..false.) call fexit()
      end if
      if (atsavred(ati).le.0.0) then
         if (saidit2.EQV..false.) then
           saidit2 = .true.
           write(ilog,*) 'Fatal. Atomic volume reduction factor for atom ',ati,' (type: ',b_type(ati),') &
 &is zero or negative. This indicates assumed radii that are incompatible with underlying assumptions. Check &
 &FMCSC_ASRPATCHFILE and FMCSC_RPATCHFILE for workarounds or use FMCSC_UNSAFE to skip this error. Further warnings omitted.'
         end if
         if (be_unsafe.EQV..false.) call fexit()
      end if
      atsavprm(ati,1) = par_IMPSOLV(6)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(6))*par_IMPSOLV(5)
      d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
      d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
      atsavprm(ati,2) = 1.0/(d2-d1)
      atsavprm(ati,3) = 1.0 - atsavprm(ati,2)*(d2-0.5)
      atsavprm(ati,4) = par_IMPSOLV(7)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(7))*par_IMPSOLV(5)
      d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
      d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
      atsavprm(ati,5) = 1.0/(d2-d1)
      atsavprm(ati,6) = 1.0 - atsavprm(ati,5)*(d2-0.5)
      atsavprm(ati,7) = (atsavmaxfr(ati)-par_IMPSOLV(5))/(1.0*max(1,nint(atbvol(ati)/par_IMPSOLV(12))))
      atsavprm(ati,8) = (atsavmaxfr(ati)-par_IMPSOLV(5))/(1.0*max(1,nint(atbvol(ati)/par_IMPSOLV(15))))
    end do
  end if
!
end
!
!-------------------------------------------------------------------
!
subroutine absinth_savprm()
!
  use energies
  use atoms
!
  implicit none
!
  integer ati
  RTYPE d1,d2
!
  do ati=1,n
    atsavprm(ati,1) = par_IMPSOLV(6)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(6))*par_IMPSOLV(5)
    d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
    d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,1))/&
 &par_IMPSOLV(3)))
    atsavprm(ati,2) = 1.0/(d2-d1)
    atsavprm(ati,3) = 1.0 - atsavprm(ati,2)*(d2-0.5)
    atsavprm(ati,4) = par_IMPSOLV(7)*atsavmaxfr(ati) + &
 &                    (1.0-par_IMPSOLV(7))*par_IMPSOLV(5)
    d1 = 1.0/(1.0 + exp(-(par_IMPSOLV(5)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
    d2 = 1.0/(1.0 + exp(-(atsavmaxfr(ati)-atsavprm(ati,4))/&
 &par_IMPSOLV(4)))
    atsavprm(ati,5) = 1.0/(d2-d1)
    atsavprm(ati,6) = 1.0 - atsavprm(ati,5)*(d2-0.5)
    atsavprm(ati,7) = (atsavmaxfr(ati)-par_IMPSOLV(5))/(1.0*nint(atbvol(ati)/par_IMPSOLV(12)))
    atsavprm(ati,8) = (atsavmaxfr(ati)-par_IMPSOLV(5))/(1.0*nint(atbvol(ati)/par_IMPSOLV(15)))
  end do
!
end
!
!-----------------------------------------------------------------------------------
!
! quantities: molvol(mt,1): net vdW-volume occupied by molecules of type mt
!             molvol(mt,2): net SAV (not overlap-, only topology-corrected)
!             resvol(rs,1): net vdW-volume occupied by residue rs
!             resvol(rs,2): net heavy atom SAV (not overlap- only topology-corrected) for res. rs
!             resvol(rs,3): number of heavy atoms in residue rs
!             sysvol(3)   : net vdW-volume
! 
subroutine absinth_molecule()
!
  use iounit
  use system
  use energies
  use molecule
  use atoms
  use polypep
  use sequen
  use math
  use aminos
  use torsn
!
  implicit none
!
  integer mt,i,rs,ii
  RTYPE d1,getblen
!
  longestmol = 1
  sysmass(1:2) = 0.
  do mt=1,nmoltyp
    do i=atmol(moltyp(mt,1),1),atmol(moltyp(mt,1),2)
      molvol(mt,1) = molvol(mt,1) + atsavred(i)*atvol(i)
      molvol(mt,2) = molvol(mt,2) + atsavmaxfr(i)*atbvol(i)
      molmass(mt) = molmass(mt) + mass(i)
    end do
    sysmass(1) = sysmass(1) + moltyp(mt,2)*molmass(mt)
    do rs=rsmol(moltyp(mt,1),1),rsmol(moltyp(mt,1),2)
      if (seqtyp(rs).ne.26) then
        molcontlen(mt) = molcontlen(mt) + rescontlen(seqtyp(rs))
      else if (seqpolty(rs).eq.'P') then
        molcontlen(mt) = molcontlen(mt) + rescontlen(1)
      else if ((seqpolty(rs).eq.'N').AND.((nuci(rs,6).gt.0).OR.(rs.gt.rsmol(moltyp(mt,1),1)))) then
        molcontlen(mt) = molcontlen(mt) + rescontlen(64)
      else if (seqpolty(rs).eq.'N') then
        molcontlen(mt) = molcontlen(mt) + rescontlen(76)
      else if ((rs.gt.rsmol(moltyp(mt,1),1)).AND.(rs.lt.rsmol(moltyp(mt,1),2))) then
        molcontlen(mt) = molcontlen(mt) + getblen(at(rs)%bb(1),at(rs+1)%bb(1))
      else
        molcontlen(mt) = molcontlen(mt) + resrad(rs)
      end if
    end do
    if ((rsmol(moltyp(mt,1),2)-rsmol(moltyp(mt,1),1)+1).gt.longestmol) then
      longestmol = rsmol(moltyp(mt,1),2)-rsmol(moltyp(mt,1),1)+1
    end if
  end do
  maxseglen = longestmol + 1
  sysvol(3) = 0.0
!
  do rs=1,nseq
    resvol(rs,:) = 0.0
    do i=1,at(rs)%nbb+at(rs)%nsc
      if (i.le.at(rs)%nbb) then
        ii = at(rs)%bb(i)
      else
        ii = at(rs)%sc(i-at(rs)%nbb)
      end if
      resvol(rs,1) = resvol(rs,1) + atsavred(ii)*atvol(ii)
      sysvol(1) = sysvol(1) + atsavred(ii)*atvol(ii)
      if (mass(ii).gt.5.0) then
        resvol(rs,2) = resvol(rs,2) + atbvol(ii)
        resvol(rs,3) = resvol(rs,3) + 1.0
      end if
    end do
    sysvol(3) = sysvol(3) + resvol(rs,2)
  end do
!
  d1 = ((3./(4.*PI))*(sysvol(1)))**(1./3.)
  sysvol(2) = 1.3504*(4./3.)*PI*((d1 + par_IMPSOLV(1))**3 - d1**3)
!  d1 =  sysvol(1)/(1.5*1.5*PI)
!  sysvol(3) = 1.5*PI*(d1+par_IMPSOLV(1))*((1.5+par_IMPSOLV(1))**2.0)
!
end
!
!----------------------------------------------------------------------------
!
! set a couple of reference arrays for residues
! ditto for expected location of fxn
!
subroutine absinth_residue()
!
  use iounit
  use polypep
  use sequen
  use atoms
  use molecule
  use aminos
!
  implicit none
!
  integer rs
  character(3) resname
!
  do rs=1,nseq
!
    at(rs)%na = at(rs)%nbb + at(rs)%nsc
    resname = amino(seqtyp(rs))
    if (resname.ne.'UNK') seqpolty(rs) = aminopolty(seqtyp(rs))
!
!   first set the reference atom for each residue
    if ((resname.eq.'NA+').OR.(resname.eq.'CL-').OR.&
 &      (resname.eq.'K+ ').OR.(resname.eq.'BR-').OR.&
 &      (resname.eq.'CS+').OR.(resname.eq.'I- ').OR.&
 &      (resname.eq.'NH4').OR.(resname.eq.'AC-').OR.&
 &      (resname.eq.'T3P').OR.(resname.eq.'SPC').OR.&
 &      (resname.eq.'URE').OR.(resname.eq.'FOR').OR.&
 &      (resname.eq.'NH2').OR.(resname.eq.'ACA').OR.&
 &      (resname.eq.'PPA').OR.(resname.eq.'NMA').OR.&
 &      (resname.eq.'CH4').OR.(resname.eq.'FOA').OR.&
 &      (resname.eq.'T4P').OR.(resname.eq.'DMA').OR.&
 &      (resname.eq.'MOH').OR.(resname.eq.'GDN').OR.&
 &      (resname.eq.'EOH').OR.(resname.eq.'MSH').OR.&
 &      (resname.eq.'1MN').OR.(resname.eq.'LCP').OR.&
 &      (resname.eq.'NO3').OR.(resname.eq.'T5P').OR.&
 &      (resname.eq.'T4E').OR.(resname.eq.'BEN').OR.&
 &      (resname.eq.'NAP').OR.(resname.eq.'O2 ')) then
       refat(rs) = at(rs)%bb(1)
    else if (resname.eq.'NMF') then
       refat(rs) = at(rs)%bb(3)
    else if ((resname.eq.'PCR').OR.(resname.eq.'TOL').OR.&
 &           (resname.eq.'IMD').OR.(resname.eq.'IME').OR.&
 &           (resname.eq.'PRP').OR.(resname.eq.'IBU').OR.&
 &           (resname.eq.'EMT').OR.(resname.eq.'NBU').OR.&
 &           (resname.eq.'2MN')) then
       refat(rs) = at(rs)%bb(2)
    else if (resname.eq.'THY') then
       refat(rs) = at(rs)%bb(8)
    else if ((resname.eq.'URA').OR.(resname.eq.'MIN')) then
       refat(rs) = at(rs)%bb(4)
    else if (resname.eq.'CYT') then
       refat(rs) = at(rs)%bb(7)
    else if (resname.eq.'ADE') then
       refat(rs) = at(rs)%bb(4)
    else if ((resname.eq.'GUA').OR.(resname.eq.'PUR')) then
       refat(rs) = at(rs)%bb(2)
    else if ((resname.eq.'R5P').OR.(resname.eq.'D5P')) then
!      the C5*
       refat(rs) = nuci(rs,4)
    else if ((resname.eq.'RIB').OR.(resname.eq.'DIB')) then
!      the C4*
       refat(rs) = nuci(rs,3)
    else if ((resname.eq.'RPU').OR.(resname.eq.'RPC').OR.&
 &           (resname.eq.'RPT').OR.(resname.eq.'RPA').OR.&
 &           (resname.eq.'RPG').OR.(resname.eq.'DPU').OR.&
 &           (resname.eq.'DPC').OR.(resname.eq.'DPT').OR.&
 &           (resname.eq.'DPA').OR.(resname.eq.'DPG').OR.&
 &           (resname.eq.'RIU').OR.(resname.eq.'RIC').OR.&
 &           (resname.eq.'RIT').OR.(resname.eq.'RIA').OR.&
 &           (resname.eq.'RIG').OR.(resname.eq.'DIU').OR.&
 &           (resname.eq.'DIC').OR.(resname.eq.'DIT').OR.&
 &           (resname.eq.'DIA').OR.(resname.eq.'DIG')) then
       refat(rs) = at(rs)%sc(2)
    else if ((resname.eq.'GLY').OR.(resname.eq.'PRO').OR.&
 &      (resname.eq.'ALA').OR.(resname.eq.'ABA').OR.&
 &      (resname.eq.'NVA').OR.(resname.eq.'VAL').OR.&
 &      (resname.eq.'LEU').OR.(resname.eq.'ILE').OR.&
 &      (resname.eq.'MET').OR.(resname.eq.'PHE').OR.&
 &      (resname.eq.'NLE').OR.(resname.eq.'SER').OR.&
 &      (resname.eq.'THR').OR.(resname.eq.'CYS').OR.&
 &      (resname.eq.'HIE').OR.(resname.eq.'HID').OR.&
 &      (resname.eq.'TYR').OR.(resname.eq.'TRP').OR.&
 &      (resname.eq.'ASP').OR.(resname.eq.'GLU').OR.&
 &      (resname.eq.'LYS').OR.(resname.eq.'ARG').OR.&
 &      (resname.eq.'GLN').OR.(resname.eq.'ASN').OR.&
 &      (resname.eq.'HYP').OR.(resname.eq.'ACE').OR.&
 &      (resname.eq.'NME').OR.(resname.eq.'AIB').OR.&
 &      (resname.eq.'DAB').OR.(resname.eq.'ORN').OR.&
 &      (resname.eq.'GAM').OR.(resname.eq.'PCA').OR.&
 &      (resname.eq.'GLH').OR.(resname.eq.'ASH').OR.&
 &      (resname.eq.'TYO').OR.(resname.eq.'CYX').OR.&
 &      (resname.eq.'LYD').OR.(resname.eq.'PTR').OR.&
 &      (resname.eq.'SEP').OR.(resname.eq.'TPO').OR.&
 &      (resname.eq.'KAC').OR.(resname.eq.'KM1').OR.&
 &      (resname.eq.'KM2').OR.(resname.eq.'KM3').OR.&
 &      (resname.eq.'HIP')) then
      refat(rs) = cai(rs)
    else if (resname.eq.'UNK') then
      cycle
    else
      write(ilog,*) 'Fatal. Reference atom undefined for residue ',&
 &resname,'. Check back later.'
      call fexit()
    end if
!
    if (refat(rs).le.0) then
      write(ilog,*) 'Fatal. Reference atom is not set for residue ',&
 &rs,' (',resname,'). This is an omission bug.'
      call fexit()
    end if
!
!   seqflag is a way to replace checks with simpler ones
!   other types are: 12 other N-cap for peptides, 13 other C-cap for peptides;
!                    26 other 5'-cap for nucleotides; 28 3'-cap for nucleotides; 50 - unknown polymer; 
    if (((seqtyp(rs).ge.1).AND.(seqtyp(rs).le.23).AND.(seqtyp(rs).ne.9)).OR.((seqtyp(rs).ge.33).AND.(seqtyp(rs).le.35)).OR.&
 &      (seqtyp(rs).eq.51).OR.(seqtyp(rs).eq.31).OR.((seqtyp(rs).ge.108).AND.(seqtyp(rs).le.119))) then
      seqflag(rs) = 2
    else if ((seqtyp(rs).eq.9).OR.(seqtyp(rs).eq.25).OR.(seqtyp(rs).eq.32)) then
      seqflag(rs) = 5
    else if (seqtyp(rs).eq.24) then
      seqflag(rs) = 8
    else if ((seqtyp(rs).eq.27).OR.(seqtyp(rs).eq.28)) then
      seqflag(rs) = 10 ! polypeptide N-cap supporting omega in rs+1 only
    else if (seqtyp(rs).eq.29) then
      seqflag(rs) = 11 ! polypeptide C-cap supporting omega on itself only
    else if (seqtyp(rs).eq.30) then
      seqflag(rs) = 13 ! polypeptide C-cap not supporting omega
    else if ((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)) then 
      seqflag(rs) = 24 ! nucleoside residue with or without sugar sampling (5'-terminal, the latter via nucsline(6,rs))
    else if ((seqtyp(rs).ge.64).AND.(seqtyp(rs).le.75)) then 
      seqflag(rs) = 22 ! nucleotide residue with or without sugar sampling
    else if ((seqtyp(rs).eq.36).OR.(seqtyp(rs).eq.37).OR.(seqtyp(rs).eq.39).OR.(seqtyp(rs).eq.40).OR.&
 &           (seqtyp(rs).eq.45).OR.(seqtyp(rs).eq.46).OR.((seqtyp(rs).ge.52).AND.(seqtyp(rs).le.56)).OR.&
 &           ((seqtyp(rs).ge.100).AND.(seqtyp(rs).le.103))) then
      seqflag(rs) = 101 ! mandatory small molecules with no possibility or intent for any torsional degrees of freedom
    else if ((seqtyp(rs).eq.60).OR.((seqtyp(rs).ge.104).AND.(seqtyp(rs).le.107))) then
      seqflag(rs) = 102 ! mandatory small molecules that are topologically completely rigidified in a torsional sense
    else if ((seqtyp(rs).eq.38).OR.(seqtyp(rs).eq.57).OR.(seqtyp(rs).eq.58).OR.(seqtyp(rs).eq.59).OR.&
 &           ((seqtyp(rs).ge.61).AND.(seqtyp(rs).le.63)).OR.((seqtyp(rs).ge.41).AND.(seqtyp(rs).le.44)).OR.&
 &           ((seqtyp(rs).ge.47).AND.(seqtyp(rs).le.50)).OR.((seqtyp(rs).ge.88).AND.(seqtyp(rs).le.99))) then
      seqflag(rs) = 103 ! mandatory small molecules with torsional degrees of freedom (whether active by default or not)
    end if
  end do
!
end
!
!--------------------------------------------------------------------
!
! a subroutine to assign parameters for bonded interactions from
! the parsed parameter file
!
subroutine assign_bndtprms()
!
  use params
  use iounit
  use polypep
  use sequen
  use aminos
  use inter
  use atoms
  use math
  use system
  use energies
  use zmatrix
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer rs,i,j,jj,k,kk,i1,i2,i3,i4,k1,k2,k3,k4,i5,k5,t1,t2
  integer, ALLOCATABLE:: dum(:,:),dum2(:,:)
  RTYPE, ALLOCATABLE:: pams(:,:)
  RTYPE testimp,getztor,getztor_ref
  logical missing,reverse,allmissing(5),check_colinear,check_colinear_ref,skiptest
!
 777 format('BLXX ',a4,a4,a4,a4,' ',6(f10.5,1x))
 778  format('BLBA ',a4,a4,a4,' ',2(f10.5,1x))
 76   format('Matching bond length potential (#',i4,') for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 86   format('Guessed bond length potential for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 88   format('Cannot guess bond length potential for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 87   format('Type: ',i3,' | Length: ',f7.2,' | Parameter 1: ',f7.2)
 77   format('Missing bond length potential specification for atoms ',&
 &i7,' and ',i7,' in residue # ',i6)
 78   format('--> Biotypes ',i4,' (',a4,') and ',i4,' (',a4,') in residu&
 &e type ',a3,'.')
 67   format('Missing bond angle potential specification for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 66   format('Matching bond angle potential (#',i4,') for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 96   format('Guessed bond angle potential for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 98   format('Cannot guess bond angle potential for atoms ',&
 &i7,',',i7,', and ',i7,' in residue # ',i6)
 97   format('Type: ',i3,' | Angle : ',f7.2,' | Parameter 1: ',f7.2)
 68   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), and ',i4,' &
 &(',a4,') in residue type ',a3,'.')
 56   format('Matching torsional potential (#',i4,') for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 57   format('No matching torsional potential specification for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)!,f10.2)
 58   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), and ',i4,' (',a4,') in residue type ',a3,'.')
 36   format('Matching CMAP potential (#',i4,') for atoms ',&
 &i7,',',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 37   format('No matching CMAP potential specification for atoms ',&
 &i7,',',i7,',',i7,',',i7,', and ',i7,' with start residue # ',i6)
 38   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), ',i4,' (',a4,'), and ',i4,' (',a4,') -> residue type of middle atom: ',a3,'.')
 46   format('Matching improper dihedral potential (#',i4,') for &
 &atoms ',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 47   format('No matching improper dihedral potential specification for &
 &atoms ',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
 48   format('--> Biotypes ',i4,' (',a4,'), ',i4,' (',a4,'), ',i4,&
 &' (',a4,'), and ',i4,' (',a4,') in residue type ',a3,'.')
116   format('Guessed improper dihedral potential for atoms ',&
 &i7,',',i7,',',i7,',and ',i7,' in residue # ',i6)
126   format('Guessed torsional potential for atoms ',&
 &i7,',',i7,',',i7,',and ',i7,' in residue # ',i6)
118   format('Cannot guess improper dihedral potential for atoms ',&
 &i7,',',i7,',',i7,', and ',i7,' in residue # ',i6)
117   format('Type: ',i3,' | Parameters : ',100(f9.3,1x))
 27   format('Proper or improper dihedral potentials disabled for &
 &atoms ',i7,',',i7,',',i7,', and ',i7,' in residue # ',i6,' due to colinear reference atoms.')
 111  format(5(a6,1x),7(g12.5))
 112  format(4(a6,1x),3(g12.5))
!
  call strlims(paramfile,t1,t2)
  allmissing(:) = .true.
  if ((bonded_report.EQV..true.).AND.(use_BOND(1).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------    Bond Length Terms    --------'
  end if
  do rs=1,nseq
    if (use_BOND(1).EQV..false.) exit
!   bond length terms (mandatory)
    do i=1,nrsbl(rs)
      missing = .true.
      i1 = bio_botyp(b_type(iaa(rs)%bl(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%bl(i,2)))
      do k=1,bo_lstsz
        k1 = bo_lst(k,1)
        k2 = bo_lst(k,2)
        if (((i1.eq.k1).AND.(i2.eq.k2)).OR.&
 &          ((i2.eq.k1).AND.(i1.eq.k2))) then
          allmissing(1) = .false.
          missing = .false.
          kk = bo_lst(k,3)
          iaa(rs)%typ_bl(i) = bo_typ(kk)
          if (bonded_report.EQV..true.) then
            write(ilog,76) kk,iaa(rs)%bl(i,1),iaa(rs)%bl(i,2),rs
            write(ilog,78) b_type(iaa(rs)%bl(i,1)),&
 &bio_code(b_type(iaa(rs)%bl(i,1))),b_type(iaa(rs)%bl(i,2)),&
 &bio_code(b_type(iaa(rs)%bl(i,2))),amino(seqtyp(rs))
          end if
          do j=1,MAXBOPAR
            iaa(rs)%par_bl(i,j) = bo_par(kk,j)
          end do
          if (iaa(rs)%typ_bl(i).eq.3) then ! GROMOS quartic: store squared equ length
            iaa(rs)%par_bl(i,3) = iaa(rs)%par_bl(i,2)*iaa(rs)%par_bl(i,2)
          end if
        end if
      end do
      if (missing.EQV..true.) then
        if (iaa(rs)%typ_bl(i).le.0) then
          write(ilog,77) iaa(rs)%bl(i,1),iaa(rs)%bl(i,2),rs
          write(ilog,78) b_type(iaa(rs)%bl(i,1)),&
 &bio_code(b_type(iaa(rs)%bl(i,1))),b_type(iaa(rs)%bl(i,2)),&
 &bio_code(b_type(iaa(rs)%bl(i,2))),amino(seqtyp(rs))
        else
          allmissing(1) = .false.
          write(ilog,86) iaa(rs)%bl(i,1),iaa(rs)%bl(i,2),rs
          write(ilog,78) b_type(iaa(rs)%bl(i,1)),&
 &bio_code(b_type(iaa(rs)%bl(i,1))),b_type(iaa(rs)%bl(i,2)),&
 &bio_code(b_type(iaa(rs)%bl(i,2))),amino(seqtyp(rs))
          write(ilog,87) iaa(rs)%typ_bl(i),iaa(rs)%par_bl(i,2),iaa(rs)%par_bl(i,1)
        end if
! &amino(seqtyp(rs))
      end if
    end do
  end do
  if ((allmissing(1).EQV..true.).AND.(use_BOND(1).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    write(ilog,*) 'No applicable assignments for bond length potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! bond angle terms (mandatory)
  if ((bonded_report.EQV..true.).AND.(use_BOND(2).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------    Bond Angle Terms     --------'
  end if
  do rs=1,nseq
    if (use_BOND(2).EQV..false.) exit
    do i=1,nrsba(rs)
      missing = .true.
      i1 = bio_botyp(b_type(iaa(rs)%ba(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%ba(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%ba(i,3)))
      do k=1,ba_lstsz
        k1 = ba_lst(k,1)
        k2 = ba_lst(k,2)
        k3 = ba_lst(k,3)
        if (i2.ne.k2) cycle
        if (((i1.eq.k1).AND.(i3.eq.k3)).OR.&
 &          ((i3.eq.k1).AND.(i1.eq.k3))) then
          missing = .false.
          allmissing(2) = .false.
          kk = ba_lst(k,4)
          iaa(rs)%typ_ba(i) = ba_typ(kk)
          if (bonded_report.EQV..true.) then
!          write(*,*) abs(iaa(rs)%par_ba(i,2)-ba_par(kk,2))
          write(ilog,66) kk,iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),&
 &iaa(rs)%ba(i,3),rs
          write(ilog,68) b_type(iaa(rs)%ba(i,1)),&
 &bio_code(b_type(iaa(rs)%ba(i,1))),b_type(iaa(rs)%ba(i,2)),&
 &bio_code(b_type(iaa(rs)%ba(i,2))),b_type(iaa(rs)%ba(i,3)),&
 &bio_code(b_type(iaa(rs)%ba(i,3))),amino(seqtyp(rs))
!              write(ilog,778) bio_code(b_type(iaa(rs)%ba(i,1))),&
! &bio_code(b_type(iaa(rs)%ba(i,2))),bio_code(b_type(iaa(rs)%ba(i,3))),ba_par(kk,1:2)
          end if
          do j=1,MAXBAPAR
            iaa(rs)%par_ba(i,j) = ba_par(kk,j)
          end do
! DEBUG          write(ilog,112) amino(seqtyp(rs)),bio_code(b_type(iaa(rs)%ba(i,1:3))),iaa(rs)%par_ba(i,1:3)
          if ((ba_typ(kk).ge.1).AND.(ba_typ(kk).le.2)) then ! note that unit conversion is not necessary for GROMOS cos-harmonic
            iaa(rs)%par_ba(i,1)=iaa(rs)%par_ba(i,1)/(RADIAN*RADIAN)
          end if
          if (iaa(rs)%typ_ba(i).eq.3) then ! GROMOS cos-harmonic: store cosine of equ angle
            iaa(rs)%par_ba(i,3) = cos(iaa(rs)%par_ba(i,2)/RADIAN)
          end if
        end if
      end do
      if (missing.EQV..true.) then
        if (iaa(rs)%typ_ba(i).le.0) then
          write(ilog,67) iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),&
 &iaa(rs)%ba(i,3),rs
          write(ilog,68) b_type(iaa(rs)%ba(i,1)),&
 &bio_code(b_type(iaa(rs)%ba(i,1))),b_type(iaa(rs)%ba(i,2)),&
 &bio_code(b_type(iaa(rs)%ba(i,2))),b_type(iaa(rs)%ba(i,3)),&
 &bio_code(b_type(iaa(rs)%ba(i,3))),amino(seqtyp(rs))
        else
          allmissing(2) = .false.
          write(ilog,96) iaa(rs)%ba(i,1),iaa(rs)%ba(i,2),iaa(rs)%ba(i,3),rs
          write(ilog,68) b_type(iaa(rs)%ba(i,1)),&
 &bio_code(b_type(iaa(rs)%ba(i,1))),b_type(iaa(rs)%ba(i,2)),&
 &bio_code(b_type(iaa(rs)%ba(i,2))),b_type(iaa(rs)%ba(i,3)),&
 &bio_code(b_type(iaa(rs)%ba(i,3))),amino(seqtyp(rs))
          write(ilog,97) iaa(rs)%typ_ba(i),iaa(rs)%par_ba(i,2),iaa(rs)%par_ba(i,1)
        end if
      end if
    end do
  end do
  if ((allmissing(2).EQV..true.).AND.(use_BOND(2).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    write(ilog,*) 'No applicable assignments for bond angle potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! improper dihedral terms (optional)
  if ((bonded_report.EQV..true.).AND.(use_BOND(3).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '-------- Improper Dihedral Terms --------'
  end if
  do rs=1,nseq
    if (use_BOND(3).EQV..false.) exit
    do i=1,nrsimpt(rs),3
      missing = .true.
      skiptest = .false.
      if (disulf(rs).gt.0) then
!       crosslinks are not yet satisfied (no input structure read or randomization performed)
!       note that this may skip valid tests for crosslinked adjacent residues also connected regularly -> sorry
        if ((atmres(iaa(rs)%impt(i,1)).eq.disulf(rs)).OR.(atmres(iaa(rs)%impt(i,2)).eq.disulf(rs)).OR.&
 &          (atmres(iaa(rs)%impt(i,3)).eq.disulf(rs)).OR.(atmres(iaa(rs)%impt(i,4)).eq.disulf(rs))) then
          skiptest = .true.
        end if
      end if
      if (skiptest.EQV..false.) then
        if ((seqtyp(atmres(iaa(rs)%impt(i,improper_conv(1)))).eq.26).OR.(seqtyp(atmres(iaa(rs)%impt(i,2))).eq.26).OR.&
 &        (seqtyp(atmres(iaa(rs)%impt(i,improper_conv(2)))).eq.26).OR.(seqtyp(atmres(iaa(rs)%impt(i,4))).eq.26)) then
          if (check_colinear_ref(iaa(rs)%impt(i,improper_conv(1)),iaa(rs)%impt(i,2),iaa(rs)%impt(i,improper_conv(2)),&
 &            iaa(rs)%impt(i,4)).EQV..true.) then
            write(ilog,27) iaa(rs)%impt(i,1:4),rs
            cycle
          end if
        else
          if (check_colinear(iaa(rs)%impt(i,improper_conv(1)),iaa(rs)%impt(i,2),iaa(rs)%impt(i,improper_conv(2)),&
 &            iaa(rs)%impt(i,4)).EQV..true.) then
            write(ilog,27) iaa(rs)%impt(i,1:4),rs
            cycle
          end if
        end if
      end if
      i1 = bio_botyp(b_type(iaa(rs)%impt(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%impt(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%impt(i,3)))
      i4 = bio_botyp(b_type(iaa(rs)%impt(i,4)))
!     all three same (triple degeneracy)
      if ((i2.eq.i3).AND.(i2.eq.i4)) then
        jj = -1
        do k=1,impt_lstsz
          k1 = impt_lst(k,1)
          k2 = impt_lst(k,2)
          k3 = impt_lst(k,3)
          k4 = impt_lst(k,4)
          if (.NOT.((k4.eq.k3).AND.(k4.eq.k2).AND.(k4.eq.i2))) cycle
          jj = jj + 1
          missing = .false.
          allmissing(3) = .false.
          kk = impt_lst(k,5)
          iaa(rs)%typ_impt(i+jj) = di_typ(kk)
          if (bonded_report.EQV..true.) then
            write(ilog,46) kk,iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4),rs
            write(ilog,48) b_type(iaa(rs)%impt(i+jj,1)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,1))),b_type(iaa(rs)%impt(i+jj,2)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,2))),b_type(iaa(rs)%impt(i+jj,3)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,3))),b_type(iaa(rs)%impt(i+jj,4)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,4))),amino(seqtyp(rs))
          end if
          do j=1,MAXDIPAR
            iaa(rs)%par_impt(i+jj,j) = di_par(kk,j)
          end do
!         convert assumed peptide convention RB to internal polymer convention RB
          if (di_typ(kk).eq.3) then
            iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
            iaa(rs)%par_impt(i+jj,4) = -iaa(rs)%par_impt(i+jj,4)
            iaa(rs)%par_impt(i+jj,6) = -iaa(rs)%par_impt(i+jj,6)
            iaa(rs)%par_impt(i+jj,8) = -iaa(rs)%par_impt(i+jj,8)
!         convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
          else if (di_typ(kk).eq.2) then
            iaa(rs)%par_impt(i+jj,1)=iaa(rs)%par_impt(i+jj,1)/(RADIAN*RADIAN)
          end if
!         re-flag polymer convention Ryckaert-Bellemans
          if (di_typ(kk).eq.3) then
            iaa(rs)%typ_impt(i+jj) = 1
          end if
        end do
!     single degeneracy
      else if (((i2.ne.i3).OR.(i2.ne.i4).OR.(i3.ne.i4)).AND.&
 &        ((i2.eq.i3).OR.(i2.eq.i4).OR.(i3.eq.i4))) then
!       in this case, the first entry will be the unique one with no relevant degeneracy (see unbond.f90)
!       note that here the biotype permutation XYZY vs. XZYY is in fact relevant since the Y's are degenerate
!       hence the redudancy check is simply omitted for a quick+dirty solution
        do jj=0,2
          i2 = bio_botyp(b_type(iaa(rs)%impt(i+jj,2)))
          i3 = bio_botyp(b_type(iaa(rs)%impt(i+jj,3)))
          i4 = bio_botyp(b_type(iaa(rs)%impt(i+jj,4)))
          do k=1,impt_lstsz
            k1 = impt_lst(k,1)
            k2 = impt_lst(k,2)
            k3 = impt_lst(k,3)
            k4 = impt_lst(k,4)
            if (.NOT.(((i1.eq.k1).AND.(i2.eq.k2)).AND.&
 &                 ((i3.eq.k3).AND.(i4.eq.k4)))) cycle
            missing = .false.
            allmissing(3) = .false.
            kk = impt_lst(k,5)
            iaa(rs)%typ_impt(i+jj) = di_typ(kk)
            if (bonded_report.EQV..true.) then
              write(ilog,46) kk,iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4),rs
              write(ilog,48) b_type(iaa(rs)%impt(i+jj,1)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,1))),b_type(iaa(rs)%impt(i+jj,2)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,2))),b_type(iaa(rs)%impt(i+jj,3)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,3))),b_type(iaa(rs)%impt(i+jj,4)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,4))),amino(seqtyp(rs))
            end if
            do j=1,MAXDIPAR
              iaa(rs)%par_impt(i+jj,j) = di_par(kk,j)
            end do
!           convert assumed peptide convention RB to internal polymer convention RB
            if (di_typ(kk).eq.3) then
              iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
              iaa(rs)%par_impt(i+jj,4) = -iaa(rs)%par_impt(i+jj,4)
              iaa(rs)%par_impt(i+jj,6) = -iaa(rs)%par_impt(i+jj,6)
              iaa(rs)%par_impt(i+jj,8) = -iaa(rs)%par_impt(i+jj,8)
!           convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
            else if (di_typ(kk).eq.2) then
              iaa(rs)%par_impt(i+jj,1)=iaa(rs)%par_impt(i+jj,1)/(RADIAN*RADIAN)
            end if
!           re-flag polymer convention Ryckaert-Bellemans
            if (di_typ(kk).eq.3) then
              iaa(rs)%typ_impt(i+jj) = 1
            end if
          end do
        end do
!     no (relevant) degeneracy
      else
        do jj=0,2
          i2 = bio_botyp(b_type(iaa(rs)%impt(i+jj,2)))
          i3 = bio_botyp(b_type(iaa(rs)%impt(i+jj,3)))
          i4 = bio_botyp(b_type(iaa(rs)%impt(i+jj,4)))
          do k=1,impt_lstsz
            k1 = impt_lst(k,1)
            k2 = impt_lst(k,2)
            k3 = impt_lst(k,3)
            k4 = impt_lst(k,4)
            if (.NOT.((((i1.eq.k1).AND.(i2.eq.k2)).AND.&
 &                 ((i3.eq.k3).AND.(i4.eq.k4))).OR.&
 &                (((i1.eq.k1).AND.(i2.eq.k3)).AND.&
 &                 ((i3.eq.k2).AND.(i4.eq.k4)))) ) cycle
            reverse = .true.
            if ((i2.eq.k2).AND.(i3.eq.k3)) reverse = .false.
            if ((improper_conv(1).ne.1).AND.(reverse.EQV..true.)) then
              j = iaa(rs)%impt(i+jj,2)
              iaa(rs)%impt(i+jj,2) = iaa(rs)%impt(i+jj,3)
              iaa(rs)%impt(i+jj,3) = j
            end if
            missing = .false.
            allmissing(3) = .false.
            kk = impt_lst(k,5)
            iaa(rs)%typ_impt(i+jj) = di_typ(kk)
            if (bonded_report.EQV..true.) then
              write(ilog,46) kk,iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4),rs
              write(ilog,48) b_type(iaa(rs)%impt(i+jj,1)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,1))),b_type(iaa(rs)%impt(i+jj,2)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,2))),b_type(iaa(rs)%impt(i+jj,3)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,3))),b_type(iaa(rs)%impt(i+jj,4)),&
 &bio_code(b_type(iaa(rs)%impt(i+jj,4))),amino(seqtyp(rs))
            end if
            do j=1,MAXDIPAR
              iaa(rs)%par_impt(i+jj,j) = di_par(kk,j)
            end do
!           convert assumed peptide convention RB to internal polymer convention RB
            if (di_typ(kk).eq.3) then
              iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
              iaa(rs)%par_impt(i+jj,4) = -iaa(rs)%par_impt(i+jj,4)
              iaa(rs)%par_impt(i+jj,6) = -iaa(rs)%par_impt(i+jj,6)
              iaa(rs)%par_impt(i+jj,8) = -iaa(rs)%par_impt(i+jj,8)
!           convert units from kcal/(mol*rad^2) to kcal/(mol*deg^2) and flip sign of 
!           equilibrium position based on crude heuristic (otherwise, we would have to
!           maintain six independent impropers for each site)
            else if (di_typ(kk).eq.2) then
              iaa(rs)%par_impt(i+jj,1)=iaa(rs)%par_impt(i+jj,1)/(RADIAN*RADIAN)
              if (iaa(rs)%par_impt(i+jj,2).ne.0.0) then

                if ((seqtyp(atmres(iaa(rs)%impt(i,1))).eq.26).OR.(seqtyp(atmres(iaa(rs)%impt(i,2))).eq.26).OR.&
 &                  (seqtyp(atmres(iaa(rs)%impt(i,3))).eq.26).OR.(seqtyp(atmres(iaa(rs)%impt(i,4))).eq.26)) then
                  testimp = getztor_ref(iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4))
                else
                  testimp = getztor(iaa(rs)%impt(i+jj,1),iaa(rs)%impt(i+jj,2),&
 &iaa(rs)%impt(i+jj,3),iaa(rs)%impt(i+jj,4))
                end if
                if (testimp*iaa(rs)%par_impt(i+jj,2).lt.0.0) then
                  iaa(rs)%par_impt(i+jj,2) = -iaa(rs)%par_impt(i+jj,2)
                end if
              end if
            end if
!           re-flag polymer convention Ryckaert-Bellemans
            if (di_typ(kk).eq.3) then
              iaa(rs)%typ_impt(i+jj) = 1
            end if
          end do
        end do
      end if
      if (missing.EQV..true.) then
        if (iaa(rs)%typ_impt(i).le.0) then
          if (bonded_report.EQV..true.) then
            write(ilog,47) iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),&
 &iaa(rs)%impt(i,3),iaa(rs)%impt(i,4),rs
            write(ilog,48) b_type(iaa(rs)%impt(i,1)),&
 &bio_code(b_type(iaa(rs)%impt(i,1))),b_type(iaa(rs)%impt(i,2)),&
 &bio_code(b_type(iaa(rs)%impt(i,2))),b_type(iaa(rs)%impt(i,3)),&
 &bio_code(b_type(iaa(rs)%impt(i,3))),b_type(iaa(rs)%impt(i,4)),&
 &bio_code(b_type(iaa(rs)%impt(i,4))),amino(seqtyp(rs))
          end if
        else
          allmissing(3) = .false.
          write(ilog,116) iaa(rs)%impt(i,1),iaa(rs)%impt(i,2),iaa(rs)%impt(i,3),iaa(rs)%impt(i,4),rs
          write(ilog,48) b_type(iaa(rs)%impt(i,1)),&
 &bio_code(b_type(iaa(rs)%impt(i,1))),b_type(iaa(rs)%impt(i,2)),&
 &bio_code(b_type(iaa(rs)%impt(i,2))),b_type(iaa(rs)%impt(i,3)),&
 &bio_code(b_type(iaa(rs)%impt(i,3))),b_type(iaa(rs)%impt(i,4)),&
 &bio_code(b_type(iaa(rs)%impt(i,4))),amino(seqtyp(rs))
          write(ilog,117) iaa(rs)%typ_impt(i),(iaa(rs)%par_impt(i,k1),k1=1,MAXDIPAR)
        end if
      end if
    end do
  end do
  if ((allmissing(3).EQV..true.).AND.(use_BOND(3).EQV..true.)) then
    write(ilog,*) 'No applicable assignments for improper dihedral potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! torsional terms (optional)
  if ((bonded_report.EQV..true.).AND.(use_BOND(4).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------     Torsional Terms     --------'
  end if
  do rs=1,nseq
    if (use_BOND(4).EQV..false.) exit
    do i=1,nrsdi(rs)
      missing = .true.
      skiptest = .false.
      if (disulf(rs).gt.0) then
!       crosslinks are not yet satisfied (no input structure read or randomization performed)
!       note that this may skip valid tests for crosslinked adjacent residues also connected regularly -> sorry
        if ((atmres(iaa(rs)%di(i,1)).eq.disulf(rs)).OR.(atmres(iaa(rs)%di(i,2)).eq.disulf(rs)).OR.&
 &          (atmres(iaa(rs)%di(i,3)).eq.disulf(rs)).OR.(atmres(iaa(rs)%di(i,4)).eq.disulf(rs))) then
          skiptest = .true.
        end if
      end if
      if (skiptest.EQV..false.) then
        if ((seqtyp(atmres(iaa(rs)%di(i,1))).eq.26).OR.(seqtyp(atmres(iaa(rs)%di(i,2))).eq.26).OR.&
 &        (seqtyp(atmres(iaa(rs)%di(i,3))).eq.26).OR.(seqtyp(atmres(iaa(rs)%di(i,4))).eq.26)) then
          if (check_colinear_ref(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4)).EQV..true.) then
            write(ilog,27) iaa(rs)%di(i,1:4),rs
            cycle
          end if
        else
          if (check_colinear(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4)).EQV..true.) then
            write(ilog,27) iaa(rs)%di(i,1:4),rs
            cycle
          end if
        end if
      end if
      i1 = bio_botyp(b_type(iaa(rs)%di(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%di(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%di(i,3)))
      i4 = bio_botyp(b_type(iaa(rs)%di(i,4)))
      do k=1,di_lstsz
        k1 = di_lst(k,1)
        k2 = di_lst(k,2)
        k3 = di_lst(k,3)
        k4 = di_lst(k,4)
        if (.NOT.((((i2.eq.k2).AND.(i3.eq.k3)).AND.&
 &                 ((i1.eq.k1).AND.(i4.eq.k4))).OR.&
 &                (((i2.eq.k3).AND.(i3.eq.k2)).AND.&
 &                 ((i1.eq.k4).AND.(i4.eq.k1)))) ) cycle
        reverse = .true.
        if ((i2.eq.k2).AND.(i3.eq.k3)) reverse = .false.
        missing = .false.
        allmissing(4) = .false.
        kk = di_lst(k,5)
        iaa(rs)%typ_di(i) = di_typ(kk)
        if (bonded_report.EQV..true.) then
          write(ilog,56) kk,iaa(rs)%di(i,1),iaa(rs)%di(i,2),&
 &iaa(rs)%di(i,3),iaa(rs)%di(i,4),rs
          write(ilog,58) b_type(iaa(rs)%di(i,1)),&
 &bio_code(b_type(iaa(rs)%di(i,1))),b_type(iaa(rs)%di(i,2)),&
 &bio_code(b_type(iaa(rs)%di(i,2))),b_type(iaa(rs)%di(i,3)),&
 &bio_code(b_type(iaa(rs)%di(i,3))),b_type(iaa(rs)%di(i,4)),&
 &bio_code(b_type(iaa(rs)%di(i,4))),amino(seqtyp(rs))
!              write(ilog,777) bio_code(b_type(iaa(rs)%di(i,1))),&
! &bio_code(b_type(iaa(rs)%di(i,2))),bio_code(b_type(iaa(rs)%di(i,3))),&
! &bio_code(b_type(iaa(rs)%di(i,4))),di_par(kk,1:4)

        end if
!        write(*,*) di_par(kk,1:6)
        do j=1,MAXDIPAR
          iaa(rs)%par_di(i,j) = di_par(kk,j)
        end do
!       convert Ryckaert-Bellemans from assumed polymer convention to internal peptide convention
        if (di_typ(kk).eq.3) then
          iaa(rs)%par_di(i,2) = -iaa(rs)%par_di(i,2)
          iaa(rs)%par_di(i,4) = -iaa(rs)%par_di(i,4)
          iaa(rs)%par_di(i,6) = -iaa(rs)%par_di(i,6)
          iaa(rs)%par_di(i,8) = -iaa(rs)%par_di(i,8)
!       convert units from kcal/(mol*deg^2) to kcal/(mol*rad^2)
        else if (di_typ(kk).eq.2) then
          iaa(rs)%par_di(i,1)=iaa(rs)%par_di(i,1)/(RADIAN*RADIAN)
        end if
!       re-flag polymer convention Ryckaert-Bellemans
        if (di_typ(kk).eq.3) then
          iaa(rs)%typ_di(i) = 1
        end if
! DEBUG        write(ilog,111) amino(seqtyp(rs)),bio_code(b_type(iaa(rs)%di(i,1:4))),iaa(rs)%par_di(i,1:7)
      end do
      if (missing.EQV..true.) then
        if (bonded_report.EQV..true.) then
          write(ilog,57) iaa(rs)%di(i,1),iaa(rs)%di(i,2),&
 &iaa(rs)%di(i,3),iaa(rs)%di(i,4),rs!,getztor(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4))
          write(ilog,58) b_type(iaa(rs)%di(i,1)),&
 &bio_code(b_type(iaa(rs)%di(i,1))),b_type(iaa(rs)%di(i,2)),&
 &bio_code(b_type(iaa(rs)%di(i,2))),b_type(iaa(rs)%di(i,3)),&
 &bio_code(b_type(iaa(rs)%di(i,3))),b_type(iaa(rs)%di(i,4)),&
 &bio_code(b_type(iaa(rs)%di(i,4))),amino(seqtyp(rs))
        end if
      end if
    end do
  end do
  if ((use_BOND(4).EQV..true.).AND.(guess_bonded.ge.3)) then
    allocate(dum(max(n,sum(nrsdi(1:nseq))),MAXVLNC+1))
    dum(:,:) = 0
    k = 0
    do rs=1,nseq
      do i=1,nrsdi(rs)
        if (check_colinear(iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4)).EQV..true.) then
          cycle
        end if
        if (iaa(rs)%typ_di(i).gt.0) cycle
        k = k + 1
        i1 = iaa(rs)%di(i,1)
        i2 = iaa(rs)%di(i,2)
        i3 = iaa(rs)%di(i,3)
        i4 = iaa(rs)%di(i,4)
        do kk=1,n12(i2)
          if (i12(kk,i2).eq.i3) cycle
          if ((izrot(i12(kk,i2))%alsz.gt.0).AND.(iz(1,i12(kk,i2)).eq.i2).AND.(iz(2,i12(kk,i2)).eq.i3)) then
            dum(k,1) = i12(kk,i2)
            exit
          end if
        end do
        if (dum(k,1).le.0) then
          do kk=1,n12(i3)
            if (i12(kk,i3).eq.i2) cycle
            if ((izrot(i12(kk,i3))%alsz.gt.0).AND.(iz(1,i12(kk,i3)).eq.i3).AND.(iz(2,i12(kk,i3)).eq.i2)) then
              dum(k,1) = i12(kk,i3)
              exit
            end if
          end do
        end if
        if (i3.gt.i2) then
          do kk=1,n12(i2)
            if (i12(kk,i2).eq.i3) then
              k1 = i2
              k2 = kk+1
              exit
            end if
          end do
        else
          do kk=1,n12(i3)
            if (i12(kk,i3).eq.i2) then
              k1 = i3
              k2 = kk+1
              exit
            end if
          end do
        end if
!        write(*,*) i1,i2,i3,i4,' maps to ',k1,k2,dum(k,1),dum(k1,k2)
        if (dum(k1,k2).eq.0) then
          call guess_bondtypes(i2,i3,dum(k,1),dum(k1,k2))
          if (dum(k1,k2).eq.-1) then
            write(ilog,*) 'Warning. UNCLASSIFIABLEBONDERROR.'
          end if
        end if
!        write(*,*) k,i2,i3,' is ',dum(k,2)
        call guess_torpot(i1,i2,i3,i4,rs,i,dum(k1,k2))
        if (iaa(rs)%typ_di(i).gt.0) then
          allmissing(4) = .false.
          write(ilog,126) iaa(rs)%di(i,1),iaa(rs)%di(i,2),iaa(rs)%di(i,3),iaa(rs)%di(i,4),rs
          write(ilog,48) b_type(iaa(rs)%di(i,1)),&
 &bio_code(b_type(iaa(rs)%di(i,1))),b_type(iaa(rs)%di(i,2)),&
 &bio_code(b_type(iaa(rs)%di(i,2))),b_type(iaa(rs)%di(i,3)),&
 &bio_code(b_type(iaa(rs)%di(i,3))),b_type(iaa(rs)%di(i,4)),&
 &bio_code(b_type(iaa(rs)%di(i,4))),amino(seqtyp(rs))
          write(ilog,117) iaa(rs)%typ_di(i),(iaa(rs)%par_di(i,k1),k1=1,MAXDIPAR)
        end if
      end do
    end do
    deallocate(dum)
  end if
  if ((allmissing(4).EQV..true.).AND.(use_BOND(4).EQV..true.)) then
    write(ilog,*) 'No applicable assignments for torsional angle potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! CMAP terms (utterly optional)
  if ((bonded_report.EQV..true.).AND.(use_BOND(5).EQV..true.)) then
    write(ilog,*)
    write(ilog,*) '--------       CMAP terms         --------'
  end if
  do rs=1,nseq
    if (use_BOND(5).EQV..false.) exit
    do i=1,nrscm(rs)
      missing = .true.
      i1 = bio_botyp(b_type(iaa(rs)%cm(i,1)))
      i2 = bio_botyp(b_type(iaa(rs)%cm(i,2)))
      i3 = bio_botyp(b_type(iaa(rs)%cm(i,3)))
      i4 = bio_botyp(b_type(iaa(rs)%cm(i,4)))
      i5 = bio_botyp(b_type(iaa(rs)%cm(i,5)))
      do k=1,cm_lstsz
        k1 = cm_lst(k,1)
        k2 = cm_lst(k,2)
        k3 = cm_lst(k,3)
        k4 = cm_lst(k,4)
        k5 = cm_lst(k,5)
!       match needs to be exact
        if (.NOT.(((i2.eq.k2).AND.(i3.eq.k3)).AND.&
 &                ((i1.eq.k1).AND.(i4.eq.k4)).AND.(i5.eq.k5)) ) cycle
        missing = .false.
        allmissing(5) = .false.
        kk = cm_lst(k,6)
        iaa(rs)%typ_cm(i) = kk
        if (bonded_report.EQV..true.) then
          write(ilog,36) iaa(rs)%typ_cm(i),iaa(rs)%cm(i,1),iaa(rs)%cm(i,2),&
 &iaa(rs)%cm(i,3),iaa(rs)%cm(i,4),iaa(rs)%cm(i,5),atmres(iaa(rs)%cm(i,3))
          write(ilog,38) b_type(iaa(rs)%cm(i,1)),&
 &bio_code(b_type(iaa(rs)%cm(i,1))),b_type(iaa(rs)%cm(i,2)),&
 &bio_code(b_type(iaa(rs)%cm(i,2))),b_type(iaa(rs)%cm(i,3)),&
 &bio_code(b_type(iaa(rs)%cm(i,3))),b_type(iaa(rs)%cm(i,4)),&
 &bio_code(b_type(iaa(rs)%cm(i,4))),b_type(iaa(rs)%cm(i,5)),&
 &bio_code(b_type(iaa(rs)%cm(i,5))),amino(seqtyp(atmres(iaa(rs)%cm(i,3))))
        end if
      end do
!     never report missing CMAPs 
    end do
  end do
  if ((allmissing(5).EQV..true.).AND.(use_BOND(5).EQV..true.)) then
    write(ilog,*) 'No applicable assignments for CMAP correction potentials found in parameter file ',paramfile(t1:t2),'.'
  end if
!
! eventually patch the default assignments from parameter file
  call read_bondedpatchfile(allmissing)
!
  do rs=1,nseq
!   sort terms such that ones with energy terms assigned are first
!   use second counter variable to mimick smaller array size
    if (use_BOND(1).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsbl(rs),MAXBOPAR+1))
      allocate(dum(nrsbl(rs),3))
      allocate(dum2(nrsbl(rs),2))
      do i=1,nrsbl(rs)
        if (iaa(rs)%typ_bl(i).gt.0) then
          k1 = k1 + 1
          do j=1,2
            dum(k1,j) = iaa(rs)%bl(i,j)
          end do
          dum(k1,3) = iaa(rs)%typ_bl(i)
          do j=1,MAXBOPAR+1
            pams(k1,j) = iaa(rs)%par_bl(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,2
            dum2(k2,j) = iaa(rs)%bl(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,2
          iaa(rs)%bl(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_bl(i) = dum(i,3)
        do j=1,MAXBOPAR+1
          iaa(rs)%par_bl(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,2
          iaa(rs)%bl(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_bl(i) = 0
      end do
      nrsbleff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(2).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsba(rs),MAXBAPAR+1))
      allocate(dum(nrsba(rs),4))
      allocate(dum2(nrsba(rs),3))
      do i=1,nrsba(rs)
        if (iaa(rs)%typ_ba(i).gt.0) then
          k1 = k1 + 1
          do j=1,3
            dum(k1,j) = iaa(rs)%ba(i,j)
          end do
          dum(k1,4) = iaa(rs)%typ_ba(i)
          do j=1,MAXBAPAR+1
            pams(k1,j) = iaa(rs)%par_ba(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,3
            dum2(k2,j) = iaa(rs)%ba(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,3
          iaa(rs)%ba(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_ba(i) = dum(i,4)
        do j=1,MAXBAPAR+1
          iaa(rs)%par_ba(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,3
          iaa(rs)%ba(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_ba(i) = 0
      end do
      nrsbaeff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(3).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsimpt(rs),MAXDIPAR))
      allocate(dum(nrsimpt(rs),5))
      allocate(dum2(nrsimpt(rs),4))
      do i=1,nrsimpt(rs)
        if (iaa(rs)%typ_impt(i).gt.0) then
          k1 = k1 + 1
          do j=1,4
            dum(k1,j) = iaa(rs)%impt(i,j)
          end do
          dum(k1,5) = iaa(rs)%typ_impt(i)
          do j=1,MAXDIPAR
            pams(k1,j) = iaa(rs)%par_impt(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,4
            dum2(k2,j) = iaa(rs)%impt(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,4
          iaa(rs)%impt(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_impt(i) = dum(i,5)
        do j=1,MAXDIPAR
          iaa(rs)%par_impt(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,4
          iaa(rs)%impt(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_impt(i) = 0
      end do
      nrsimpteff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(4).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(pams(nrsdi(rs),MAXDIPAR))
      allocate(dum(nrsdi(rs),5))
      allocate(dum2(nrsdi(rs),4))
      do i=1,nrsdi(rs)
        if (iaa(rs)%typ_di(i).gt.0) then
          k1 = k1 + 1
          do j=1,4
            dum(k1,j) = iaa(rs)%di(i,j)
          end do
          dum(k1,5) = iaa(rs)%typ_di(i)
          do j=1,MAXDIPAR
            pams(k1,j) = iaa(rs)%par_di(i,j)
          end do
        else
          k2 = k2 + 1
          do j=1,4
            dum2(k2,j) = iaa(rs)%di(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,4
          iaa(rs)%di(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_di(i) = dum(i,5)
        do j=1,MAXDIPAR
          iaa(rs)%par_di(i,j) = pams(i,j)
        end do
      end do
      do i=k1+1,k1+k2
        do j=1,4
          iaa(rs)%di(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_di(i) = 0
      end do
      nrsdieff(rs) = k1
      deallocate(pams)
      deallocate(dum)
      deallocate(dum2)
    end if
    if (use_BOND(5).EQV..true.) then
      k1 = 0
      k2 = 0
      allocate(dum(nrscm(rs),6))
      allocate(dum2(nrscm(rs),5))
      do i=1,nrscm(rs)
        if (iaa(rs)%typ_cm(i).gt.0) then
          k1 = k1 + 1
          do j=1,5
            dum(k1,j) = iaa(rs)%cm(i,j)
          end do
          dum(k1,6) = iaa(rs)%typ_cm(i)
        else
          k2 = k2 + 1
          do j=1,5
            dum2(k2,j) = iaa(rs)%cm(i,j)
          end do
        end if
      end do
      do i=1,k1
        do j=1,5
          iaa(rs)%cm(i,j) = dum(i,j)
        end do
        iaa(rs)%typ_cm(i) = dum(i,6)
      end do
      do i=k1+1,k1+k2
        do j=1,5
          iaa(rs)%cm(i,j) = dum2(i-k1,j)
        end do
        iaa(rs)%typ_cm(i) = 0
      end do
      nrscmeff(rs) = k1
      deallocate(dum)
      deallocate(dum2)
    end if
  end do
!
  if ((allmissing(1).EQV..true.).AND.(use_BOND(1).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    use_BOND(1) = .false.
    scale_BOND(1) = 0.0
  end if
  if ((allmissing(2).EQV..true.).AND.(use_BOND(2).EQV..true.).AND.&
 &    (fycxyz.ne.2)) then
    use_BOND(2) = .false.
    scale_BOND(2) = 0.0
  end if
  if ((allmissing(3).EQV..true.).AND.(use_BOND(3).EQV..true.)) then
    use_BOND(3) = .false.
    scale_BOND(3) = 0.0
  end if
  if ((allmissing(4).EQV..true.).AND.(use_BOND(4).EQV..true.)) then
    use_BOND(4) = .false.
    scale_BOND(4) = 0.0
  end if
  if ((allmissing(5).EQV..true.).AND.(use_BOND(5).EQV..true.)) then
    use_BOND(5) = .false.
    scale_BOND(5) = 0.0
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine guess_types(inval,inmass,innam,ljty,boty)
!
  use params
!
  implicit none
!
  character(4) innam
  integer inval,ljty,boty,k
  RTYPE inmass,thresh
!
  ljty = -1
  thresh = 0.01
  do while (thresh.lt.HUGE(thresh)/100.0) 
    do k=1,n_ljtyp
      if ((abs(inmass-lj_weight(k)).le.thresh).AND.(lj_val(k).eq.inval)) then
        ljty = k
        exit
      end if
    end do
    if (ljty.eq.-1) then ! has to trigger eventually
      do k=1,n_ljtyp
        if (abs(inmass-lj_weight(k)).le.thresh) then
          ljty = k
          exit
        end if
      end do
    end if
    if (ljty.eq.-1) then
      do k=1,n_ljtyp
        if (inval.eq.lj_val(k)) then
          ljty = k
          exit
        end if
      end do
    end if
    if (ljty.gt.0) exit
    thresh = thresh*10.0
  end do
!
  n_biotyp = n_biotyp + 1
  boty = n_biotyp
  if ((innam(1:1).eq.' ').OR.(innam(1:1).eq.'1').OR.(innam(1:1).eq.'2').OR.(innam(1:1).eq.'3').OR.(innam(1:1).eq.'4').OR.&
 &(innam(1:1).eq.'5').OR.(innam(1:1).eq.'6').OR.(innam(1:1).eq.'7').OR.(innam(1:1).eq.'8').OR.(innam(1:1).eq.'9').OR.&
 &(innam(1:1).eq.'0')) then
    bio_code(boty) = innam(2:4)
  else 
    bio_code(boty) = innam(1:3)
  end if
  bio_ljtyp(boty) = ljty
  bio_ctyp(boty) = 0
  bio_botyp(boty) = 0
!
end
!
!-----------------------------------------------------------------------------------------------------
!
subroutine guess_bondtypes(ia,ib,roti,bndty)
!
  use zmatrix
  use atoms
  use iounit
!
  implicit none
!
  integer ia,ib,roti,i,bndty,ii1,ii2,i2,i1,iii,isubs(2,20),noops,i3,i4,k
  RTYPE getztor,getbang,dih(MAXVLNC*MAXVLNC),impts(2),dtmp,threshs(3)
  character(3) attc
!
 344 format(i5,i5,a)
!
  threshs(1) = 5.0 ! for planar rings
  threshs(2) = 10.0 ! for linear planar groups
!
  if ((n12(ia).le.1).OR.(n12(ib).le.1)) then
    bndty = 0 ! terminal bond of any type
    call fexit()
    return 
  end if
!
  bndty = -1 ! general flag for unclassified bond
!
! compute some generic properties of the bond
  attc = '   '
  isubs(:,:) = 0
  i1 = ia ! reorder to C < N < O < P < S
  i2 = ib
  if ((atnam(ia).eq.'S  ').OR.&
 &    ((atnam(ia).eq.'P  ').AND.(atnam(ib).ne.'S  ')).OR.&
 &    ((atnam(ia).eq.'O  ').AND.(atnam(ib).ne.'S  ').AND.(atnam(ib).ne.'P  ')).OR.&
 &    ((atnam(ia).eq.'N  ').AND.(atnam(ib).eq.'C  '))) then
    i1 = ib
    i2 = ia
  end if
  if (atnam(i1).eq.atnam(i2)) then ! reorder by valence
    if (n12(i1).gt.n12(i2)) then
      ii1 = i1
      i1 = i2
      i2 = ii1
    end if 
  end if
  isubs(1:2,1) = 1
  do ii1=1,n12(i1)
    if (i12(ii1,i1).eq.i2) cycle
    if (attc.eq.'   ') attc = atnam(i12(ii1,i1))
    if ((isubs(1,1).eq.1).AND.(attc.eq.atnam(i12(ii1,i1)))) then
      isubs(1,1) = 1
    else
      isubs(1,1) = 0
    end if
    attc = atnam(i12(ii1,i1))
    if (attc.eq.'H  ') isubs(1,2) = isubs(1,2) + 1
    if (attc.eq.'H  ') isubs(1,20) = isubs(1,20) + 1
    if ((attc.eq.'CL ').AND.(n12(i12(ii1,i1)).eq.1)) isubs(1,2) = isubs(1,2) + 1
    if ((attc.eq.'F  ').AND.(n12(i12(ii1,i1)).eq.1)) isubs(1,2) = isubs(1,2) + 1
    if ((attc.eq.'BR ').AND.(n12(i12(ii1,i1)).eq.1)) isubs(1,2) = isubs(1,2) + 1
    if ((attc.eq.'I  ').AND.(n12(i12(ii1,i1)).eq.1)) isubs(1,2) = isubs(1,2) + 1
    if (attc.eq.'C  ') isubs(1,3) = isubs(1,3) + 1
    if (attc.eq.'N  ') isubs(1,4) = isubs(1,4) + 1
    if (attc.eq.'O  ') isubs(1,5) = isubs(1,5) + 1
    if (attc.eq.'P  ') isubs(1,8) = isubs(1,8) + 1
    if (attc.eq.'S  ') isubs(1,9) = isubs(1,9) + 1
    if (n12(i12(ii1,i1)).eq.1) isubs(1,11) = isubs(1,11) + 1
    if (((n12(i12(ii1,i1)).le.2).AND.(attc.eq.'N  ')).OR.((n12(i12(ii1,i1)).eq.3).AND.(attc.eq.'C  ')).OR.&
 &    ((n12(i12(ii1,i1)).le.2).AND.(attc.eq.'P  ')).OR.((n12(i12(ii1,i1)).eq.1).AND.((attc.eq.'S  ').OR.(attc.eq.'O  ')))) then
      isubs(1,12) = isubs(1,12) + 1
    end if
  end do
  do ii1=1,n12(i2)
    if (i12(ii1,i2).eq.i1) cycle
    if (attc.eq.'   ') attc = atnam(i12(ii1,i2))
    if ((isubs(2,1).eq.1).AND.(attc.eq.atnam(i12(ii1,i2)))) then
      isubs(2,1) = 1
    else
      isubs(2,1) = 0
    end if
    attc = atnam(i12(ii1,i2))
    if (attc.eq.'H  ') isubs(2,2) = isubs(2,2) + 1
    if (attc.eq.'H  ') isubs(2,20) = isubs(2,20) + 1
    if ((attc.eq.'CL ').AND.(n12(i12(ii1,i2)).eq.1)) isubs(2,2) = isubs(2,2) + 1
    if ((attc.eq.'F  ').AND.(n12(i12(ii1,i2)).eq.1)) isubs(2,2) = isubs(2,2) + 1
    if ((attc.eq.'BR ').AND.(n12(i12(ii1,i2)).eq.1)) isubs(2,2) = isubs(2,2) + 1
    if ((attc.eq.'I  ').AND.(n12(i12(ii1,i2)).eq.1)) isubs(2,2) = isubs(2,2) + 1
    if (attc.eq.'C  ') isubs(2,3) = isubs(2,3) + 1
    if (attc.eq.'N  ') isubs(2,4) = isubs(2,4) + 1
    if (attc.eq.'O  ') isubs(2,5) = isubs(2,5) + 1
    if (attc.eq.'P  ') isubs(2,8) = isubs(2,8) + 1
    if (attc.eq.'S  ') isubs(2,9) = isubs(2,9) + 1
    if (n12(i12(ii1,i2)).eq.1) isubs(2,11) = isubs(2,11) + 1
    if (((n12(i12(ii1,i2)).le.2).AND.(attc.eq.'N  ')).OR.((n12(i12(ii1,i2)).eq.3).AND.(attc.eq.'C  ')).OR.&
 &    ((n12(i12(ii1,i2)).le.2).AND.(attc.eq.'P  ')).OR.((n12(i12(ii1,i2)).eq.1).AND.((attc.eq.'S  ').OR.(attc.eq.'O  ')))) then
      isubs(2,12) = isubs(2,12) + 1
    end if
  end do
  isubs(1,6) = n12(i1) - sum(isubs(1,2:5)) - sum(isubs(1,8:9)) - 1 
  isubs(2,6) = n12(i2) - sum(isubs(2,2:5)) - sum(isubs(2,8:9)) - 1
  impts(:) = 0.0
  if (n12(i1).eq.3) then     
    impts(1) = getztor(i1,i12(1,i1),i12(2,i1),i12(3,i1))
  end if
  if (n12(i2).eq.3) then     
    impts(2) = getztor(i2,i12(1,i2),i12(2,i2),i12(3,i2))
  end if
  do i=1,2
    if (abs(impts(i)).gt.threshs(min(roti,1)+1)) then
      dtmp = impts(i) - 180.
      if (dtmp.lt.-180.) dtmp = dtmp + 360.
      if (abs(dtmp).gt. threshs(min(roti,1)+1)) then
        isubs(i,7) = isubs(i,7) + 1
      end if
    end if
  end do
  iii = 0
  do ii1=1,n12(i1)
    if (i12(ii1,i1).eq.i2) cycle
    do ii2=1,n12(i2)
      if (i12(ii2,i2).eq.i1) cycle
      iii = iii + 1
      dih(iii) = getztor(i12(ii1,i1),i1,i2,i12(ii2,i2))
    end do
  end do
  noops = 0
  do i=1,iii
    if (abs(dih(i)).gt.threshs(min(roti,1)+1)) then
      dtmp = dih(i) - 180.
      if (dtmp.lt.-180.) dtmp = dtmp + 360.
      if (abs(dtmp).gt. threshs(min(roti,1)+1)) then
        noops = noops + 1
      end if
    end if
  end do
  if (n12(i2).gt.3) isubs(2,7) = 1
  if (n12(i1).gt.3) isubs(1,7) = 1
  if (n12(i1).eq.2) then
    if (getbang(i12(1,i1),i1,i12(2,i1)).ge.160.0) isubs(1,10) = 1
  end if
  if (n12(i2).eq.2) then
    if (getbang(i12(1,i2),i2,i12(2,i2)).ge.160.0) isubs(2,10) = 1
  end if
  k = 15
  do ii1=1,n12(i1)
    iii = 0
    if (i12(ii1,i1).eq.i2) cycle
    do ii2=1,n12(i12(ii1,i1))
      if (i12(ii2,i12(ii1,i1)).eq.i1) cycle
      iii = iii + 1
      if (getbang(i1,i12(ii1,i1),i12(ii2,i12(ii1,i1))).ge.175.0) then
        dih(iii) = 0.0
      else
        dih(iii) = getztor(i2,i1,i12(ii1,i1),i12(ii2,i12(ii1,i1)))
      end if
    end do
    if (iii.eq.0) isubs(1,k) = -1
    do i=1,iii
      if (abs(dih(i)).gt.threshs(min(roti,1)+1)) then
        dtmp = dih(i) - 180.
        if (dtmp.lt.-180.) dtmp = dtmp + 360.
        if (abs(dtmp).gt. threshs(min(roti,1)+1)) then
          isubs(1,k) = isubs(1,k) + 1
        end if
      end if
    end do
    if (isubs(1,k).gt.0) isubs(1,13) = isubs(1,13) + 1
    if (isubs(1,k).eq.0) isubs(1,14) = isubs(1,14) + 1
    k = k + 1
  end do
  k = 15
  do ii1=1,n12(i2)
    iii = 0
    if (i12(ii1,i2).eq.i1) cycle
    do ii2=1,n12(i12(ii1,i2))
      if (i12(ii2,i12(ii1,i2)).eq.i2) cycle
      iii = iii + 1
      if (getbang(i2,i12(ii1,i2),i12(ii2,i12(ii1,i2))).ge.175.0) then
        dih(iii) = 0.0
      else
        dih(iii) = getztor(i1,i2,i12(ii1,i2),i12(ii2,i12(ii1,i2)))
      end if
    end do
    if (iii.eq.0) isubs(2,k) = -1
    do i=1,iii
      if (abs(dih(i)).gt.threshs(min(roti,1)+1)) then
        dtmp = dih(i) - 180.
        if (dtmp.lt.-180.) dtmp = dtmp + 360.
        if (abs(dtmp).gt. threshs(min(roti,1)+1)) then
          isubs(2,k) = isubs(2,k) + 1
        end if
      end if
    end do
    if (isubs(2,k).gt.0) isubs(2,13) = isubs(2,13) + 1
    if (isubs(2,k).eq.0) isubs(2,14) = isubs(2,14) + 1
    k = k + 1
  end do
!
 52 format(i4,1x,a3,' H: ',i1,' C: ',i1,' N: ',i1,' O: ',i1,' P: ',i1,' S: ',i1,' X: ',i1,' Sym: ',i1,' Ter: ',i1,&
 &' Hyp: ',i1,' Lin: ',i1,' OOP: ',i1,' OOPSubs: ',i1,' IPSubs: ',i1) ! 5(i2,1x))
!  write(*,*) noops
  write(*,52) i1,atnam(i1),(isubs(1,i),i=2,5),isubs(1,8),isubs(1,9),isubs(1,6),isubs(1,1),isubs(1,11),isubs(1,12),isubs(1,10),&
 &isubs(1,7),isubs(1,13),isubs(1,14) !(isubs(1,i),i=15,17)
  write(*,52) i2,atnam(i2),(isubs(2,i),i=2,5),isubs(2,8),isubs(2,9),isubs(2,6),isubs(2,1),isubs(2,11),isubs(2,12),isubs(2,10),&
 &isubs(2,7),isubs(2,13),isubs(2,14) !(isubs(2,i),i=15,17)
!
! first deal with generic cases where we do not need to distinguish further
! rotatable, aliphatic sp3-sp3 C/N
  if ((n12(ia).eq.4).AND.(n12(ib).eq.4).AND.(roti.gt.0).AND.((atnam(ia).eq.'N  ').OR.(atnam(ia).eq.'C  ')).AND.&
 &((atnam(ib).eq.'N  ').OR.(atnam(ib).eq.'C  '))) then
    bndty = 101
    write(ilog,344) ia,ib,' is a rotatable bond between sp3 N/C atoms'
    return
  end if
!
  if ((n12(i1).le.3).AND.(n12(i2).le.3).AND.(roti.le.0)) then
    if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(n12(i1).eq.3).AND.(n12(i2).eq.3)) then
      write(ilog,344) ia,ib,' is a planar ring bond between trigonal sp2 atoms.'
      bndty = 1
    else if ((noops.gt.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(n12(i1).eq.3).AND.(n12(i2).eq.3)) then
      write(ilog,344) ia,ib,' is an ring bond between two trigonal and separately planar but misaligned centers.'
      bndty = 201
    else if ((noops.le.0).AND.(n12(i1).eq.2).AND.(n12(i2).eq.2)) then
      write(ilog,344) ia,ib,' is most likely a ring bond between two heteroatoms (probably -N=N-, =N-S-, =N-O-, =N-N=).'
      bndty = 3
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0)) then
      write(ilog,344) ia,ib,' is a planar ring bond between a trigonal sp2 atom and a lower-valence planar atom (N/C=N, =C-S/O/N).'
      bndty = 2
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2)) then
      write(ilog,344) ia,ib,' is most likely a bond between heteroatoms in a flexible ring (S-S, O-S, O-O, =N-S, =N-O).'
      bndty = 102
    else if ((isubs(2,3).ge.1).AND.(isubs(1,5).eq.1).AND.(isubs(1,8).eq.0).AND.(atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'N  ')) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in a distorted lactam (e.g., beta-lactam).'
      bndty = 280
    else
      write(ilog,344) ia,ib,' is an unclassified bond (most likely between heteroatoms) in a flexible ring.'
    end if
    return
  end if
!
  if ((n12(i1).ge.4).AND.(n12(i2).ge.4).AND.(roti.le.0)) then
    if ((noops.gt.0).AND.((atnam(i1).eq.'C  ').OR.(atnam(i1).eq.'N  ')).AND.((atnam(i2).eq.'C  ').OR.(atnam(i2).eq.'N  ')).AND.&
 &      (sum(isubs(1,2:3)).eq.(n12(i1)-1)).AND.(sum(isubs(2,2:3)).eq.(n12(i2)-1))) then
      write(ilog,344) ia,ib,' is a single bond between two aliphatic sp3 N/C atoms in a saturated ring.'
      bndty = 206
    else if (((atnam(i1).eq.'C  ').OR.(atnam(i1).eq.'N  ')).AND.((atnam(i2).eq.'C  ').OR.(atnam(i2).eq.'N  '))) then
      write(ilog,344) ia,ib,' is a single bond between two sp3 N/C atoms in a saturated ring.'
      bndty = 207
    else if ((atnam(i2).eq.'S  ').AND.((atnam(i1).eq.'C  ').OR.(atnam(i1).eq.'N  '))) then
      write(ilog,344) ia,ib,' is most likely an S-C bond in a cyclic sulfonamide or sulfonyl compound.'
      bndty = 208
    else if ((atnam(i2).eq.'P  ').AND.((atnam(i1).eq.'C  ').OR.(atnam(i1).eq.'N  '))) then
      write(ilog,344) ia,ib,' is most likely a P-C bond in a cyclic phosphonamide, phosphonyl group, phosphonate, phosphinoxide, &
 &phosphinamide, iminophosphorane, phosphonium ion, etc.'
      bndty = 277
    else
      write(ilog,344) ia,ib,' is an unclassified bond in a flexible ring most likely involving heteroatoms.'
    end if
    return
  end if
!
  if ((n12(i1).ge.4).AND.(n12(i2).eq.3).AND.(roti.le.0).AND.(isubs(2,7).ge.1)) then
    if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'N  ')) then
      write(ilog,344) ia,ib,' is most likely a single N-C bond between sp3 atoms in a cyclic amine, hydrazone, hydrazide, &
 &N-oxide, or to a distorted lactam (e.g., beta-lactam), etc.'
      bndty = 278
    else if ((atnam(i1).eq.'C  ').AND.((atnam(i2).eq.'S  ').OR.(atnam(i2).eq.'P  '))) then
      write(ilog,344) ia,ib,' is most likely a single C-P/S bond between sp3 atoms in a cyclic sulfoxide (or phosphane).'
      bndty = 279
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'C  ').AND.((n12(i1)+n12(i2)).eq.7)) then
    if ((roti.gt.0).AND.(sum(isubs(1:2,2:3)).eq.5)) then
      write(ilog,344) ia,ib,' is a single C-C bond between sp3 and sp2 carbons in an unsaturated hydro/halocarbon moiety (e.g., &
 &alkylated aromatics).'
      bndty = 245
    else if ((roti.le.0).AND.(sum(isubs(1:2,2:3)).eq.5)) then
      write(ilog,344) ia,ib,' is a single C-C bond between sp3 and sp2 carbons in a partially unsaturated ring (e.g., &
 &cyclopentene).'
      bndty = 246
    else if (((isubs(1,5)+isubs(1,9)).eq.2).AND.(sum(isubs(2,2:3)).eq.3)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between sp3 and sp2 carbons in a carboxylate, carboxylic ester, or &
 &lactone (incl. sulfur derivatives).'
      bndty = 247
    else if (((isubs(1,5)+isubs(1,9)).eq.1).AND.(isubs(1,14).ge.1).AND.(sum(isubs(2,2:3)).eq.3)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between sp3 and sp2 carbons in an aliphatic amide, ketone, aldehyde&
 &, enol, lactam, etc. (including sulfur derivatives).'
      bndty = 248
    else if ((isubs(1,4).eq.2).AND.(isubs(1,14).ge.2).AND.(sum(isubs(2,2:3)).eq.3).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between sp3 and sp2 carbons in an aliphatic amidine or alkyl-&
 &substituted (aromatic) heterocycle.'
      bndty = 249
    else if ((isubs(1,4).eq.1).AND.(isubs(1,14).ge.1).AND.(sum(isubs(2,2:3)).eq.3).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between sp3 and sp2 carbons in an aliphatic oxime, enamine, &
 &hydrazone, or alkyl-substituted (aromatic) heterocycle.'
      bndty = 250
    else if (((isubs(1,5)+isubs(1,9)).ge.1).AND.(isubs(1,14).ge.1).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between sp3 and sp2 carbons in an alpha-functionalized amide, &
 &ketone, ester, aldehyde, carboxylate, etc. (including sulfur derivatives).'
      bndty = 251
    else if ((isubs(1,4).ge.1).AND.(isubs(1,14).ge.1).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between sp3 and sp2 carbons in an alpha-functionalized &
 &amidine, oxime, enamine, hydrazone, or alkyl-substituted (aromatic) heterocycle.'
      bndty = 256 
    else if ((roti.gt.0).AND.(sum(isubs(1,2:3)).eq.2).AND.((isubs(2,5)+isubs(2,9)).eq.1).AND.(sum(isubs(2,2:3)).eq.2)) then
      write(ilog,344) ia,ib,' is a single C-C bond between sp3 and sp2 carbons in an arylated or alkenylated alcohol, &
 &(thio)ether, ester, etc.'
      bndty = 254
    else if ((roti.gt.0).AND.(sum(isubs(1,2:3)).eq.2).AND.(isubs(2,4).eq.1).AND.(sum(isubs(2,2:3)).eq.2)) then
      write(ilog,344) ia,ib,' is a single C-C bond between sp3 and sp2 carbons in an arylated or alkenylated amine, amide, etc.'
      bndty = 255
    else if (roti.le.0) then
      write(ilog,344) ia,ib,' is a single C-C bond between sp3 and sp2 carbons in a heterocycle or polyfunctional ring.'
      bndty = 252
    else
      write(ilog,344) ia,ib,' is an unclassified sp2-sp3 C-C bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'C  ').AND.(n12(i1).eq.2).AND.(n12(i2).eq.4)) then
    if ((isubs(1,10).eq.1).AND.(isubs(1,11).eq.1)) then
      write(ilog,344) ia,ib,' is most likely an sp-sp3 C-C bond to a nitrile (cyanide) or terminal alkyne.'
      bndty = 257
    else if ((isubs(1,10).eq.1).AND.(isubs(1,3).eq.1)) then
      write(ilog,344) ia,ib,' is most likely an sp-sp3 C-C bond in an alkyne.'
      bndty = 257
    else
      write(ilog,344) ia,ib,' is an unclassified sp-sp3 C-C bond.'
      bndty = 257
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'C  ').AND.(n12(i1).eq.2).AND.(n12(i2).eq.3)) then
    if ((isubs(1,10).eq.1).AND.(isubs(1,11).eq.1)) then
      write(ilog,344) ia,ib,' is most likely an sp-sp2 C-C bond in a conjugated nitrile (cyanide) or terminal alkyne.'
      bndty = 257
    else if ((isubs(1,10).eq.1).AND.(isubs(1,3).eq.1).AND.(isubs(1,11).eq.0)) then
      write(ilog,344) ia,ib,' is most likely an sp-sp2 C-C bond in a conjugated alkyne or an allene derivative (conformational &
 &preferences not representable).'
      bndty = 401 
    else
      write(ilog,344) ia,ib,' is an unclassified sp-sp2 C-C bond.'
      bndty = 257
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'C  ').AND.(n12(i1).eq.2).AND.(n12(i2).eq.2)) then
    if ((isubs(1,10).eq.1).AND.(isubs(2,10).eq.1).AND.(minval(isubs(1:2,12)).le.0).AND.(sum(isubs(1:2,2:3)).eq.2)) then
      write(ilog,344) ia,ib,' is most likely an sp-sp triple bond in an alkyne.'
      bndty = 258
    else if ((isubs(1,10).eq.1).AND.(isubs(2,10).eq.1).AND.(sum(isubs(1:2,2:3)).eq.2)) then
      write(ilog,344) ia,ib,' is most likely an sp-sp single or triple bond in a polyyne or a triple bond in a fully &
 &conjugated alkyne (conformational preferences not representable).'
      bndty = 402
    else
      write(ilog,344) ia,ib,' is an unclassified sp-sp C-C bond.'
      bndty = 258
    end if
    return
  end if
!
! unlike C=N and N=N bonds, C=C bonds are typically so isomer-stable that "clean" compounds are purchasable;
! this is unfortunately weakened in conjugated systems: doi: 10.1021/cr0104375
  if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(roti.gt.0).AND.(atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'C  ')) then
    if ((sum(isubs(1,4:5))+isubs(1,9)).gt.(sum(isubs(2,4:5))+isubs(2,9))) then
      i3 = 2
      i4 = 1
    else
      i3 = 1
      i4 = 2
    end if
    if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.((isubs(i3,2).eq.2).OR.(isubs(i4,2).eq.2))) then
      write(ilog,344) ia,ib,' is most likely a C=C bond in a terminal alkene.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.&
 &           (((isubs(i3,13).eq.0).AND.(isubs(i3,14).eq.0).AND.(isubs(i3,5).eq.0).AND.(isubs(i3,9).eq.0)).OR.&
 &            ((isubs(i4,13).eq.0).AND.(isubs(i4,14).eq.0).AND.(isubs(i4,5).eq.0).AND.(isubs(i4,9).eq.0)))) then
      write(ilog,344) ia,ib,' is most likely a C=C bond in a terminal, halogenated alkene.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(sum(isubs(1:2,2:3)).eq.4).AND.(sum(isubs(1:2,14)).eq.0)) then
      write(ilog,344) ia,ib,' is most likely a nonisomerizable C=C bond with only saturated hydro/halocarbon substituents.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.((sum(isubs(1:2,4:5))+sum(isubs(1:2,9))).eq.4)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between trigonal, coplanar sp2 carbons in a conjugated &
 &polyfunctional vicinal moiety (e.g., amidinopyrimidine, oxalic acid ester, thiazole carboxylic acid, etc.)'
      if ((isubs(1,5).eq.1).AND.(isubs(2,5).eq.1)) then
        bndty = 312
      else if ((isubs(1,9).eq.1).AND.(isubs(2,9).eq.1)) then
        bndty = 313
      else
        bndty = 13
      end if
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.((sum(isubs(i4,4:5))+isubs(i4,9)).eq.2).AND.&
 &           (isubs(i4,13).eq.0).AND.(isubs(i4,11).gt.0)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between trigonal, coplanar sp2 carbons in a conjugated &
 &amide, ester, carboxylic acid, carbamate, etc. (e.g., benzamide).'
      bndty = 310
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(isubs(i4,4).eq.2)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between trigonal, coplanar sp2 carbons in an arylated/alkenylated &
 &heteroaromatic cycle or conjugated amidine.'
      bndty = 311
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(sum(isubs(i3,2:3)).eq.2).AND.(isubs(i4,4).eq.1)) then
      write(ilog,344) ia,ib,' is most likely a C-C or C=C bond to an alkenylated heteroaromatic cycle, in an enamine, or in a &
 &conjugated amide, oxime, etc.' 
      bndty = 301
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.((sum(isubs(i4,4:5))+isubs(i4,9)).eq.1).AND.&
 &           (sum(isubs(i4,13:14)).le.1).AND.(sum(isubs(i3,2:3)).eq.2).AND.(isubs(i3,13).eq.0).AND.(isubs(i4,11).gt.0)) then
      write(ilog,344) ia,ib,' is most likely a single C-C bond between trigonal, coplanar sp2 carbons in a conjugated &
 &ketone or the C=C bond in the corresponding enol(ate).'
      bndty = 314
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(sum(isubs(1:2,14)).eq.0).AND.(maxval(isubs(1:2,13)).gt.0)) then
      write(ilog,344) ia,ib,' is most likely a nonisomerizable, unconjugated C=C bond.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(sum(isubs(1:2,2:3)).eq.4).AND.(sum(isubs(1:2,11)).eq.4)) then
      write(ilog,344) ia,ib,' is most likely ethene or a halogenated derivative.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(sum(isubs(1:2,2:3)).eq.4).AND.(sum(isubs(1:2,13)).eq.0)) then
      write(ilog,344) ia,ib,' is most likely a conjugated C=C or C-C bond with only hydro/halocarbon substituents.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(sum(isubs(1:2,2:3)).eq.4)) then
      write(ilog,344) ia,ib,' is most likely a C=C with mixed saturated/unsaturated hydro/halocarbon substituents.' 
      bndty = 16
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(isubs(i3,4).eq.1).AND.(isubs(i4,4).eq.1).AND.&
 &           (sum(isubs(1:2,13)).eq.0)) then
      write(ilog,344) ia,ib,' is most likely a conjugated C-C bond between heteroaromatic cycles or in a diimine.' 
      bndty = 315
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0)) then
      write(ilog,344) ia,ib,' is an unclassified but planar C-C or C=C bond between two trigonal sp2 carbons.'
      bndty = 16
    else if ((noops.gt.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(maxval(isubs(1:2,13)).le.0).AND.(maxval(isubs(1:2,11)).le.0)) then
      write(ilog,344) ia,ib,' is a rotatable bond between two trigonal and separately planar but misaligned carbon centers &
 &(most likely between separate ring systems).'
      bndty = 202
    else if ((noops.gt.0).AND.(maxval(isubs(1:2,7)).le.0)) then
      write(ilog,344) ia,ib,' is a rotatable bond between two trigonal and separately planar but misaligned carbon centers &
 &(diones, derivatives of benzoic acid, biphenyl?).'
      bndty = 202
    else
      write(ilog,344) ia,ib,' is an unclassified bond between sp2 carbons.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.((atnam(i2).eq.'O  ').OR.(atnam(i2).eq.'S  '))) then
    if ((n12(i1).eq.3).AND.(n12(i2).eq.2).AND.(roti.gt.0)) then
      if ((isubs(1,5)+isubs(1,9)).eq.2) then
        bndty = 302
        write(ilog,344) ia,ib,' is most likely a C-O bond in a carbonate or sulfur derivative.'
      else if (((isubs(1,5)+isubs(1,9)).eq.1).AND.(isubs(1,11).ge.1).AND.(maxval(isubs(1:2,7)).le.0)) then
        write(ilog,344) ia,ib,' is most likely a C-O bond in an ester, carbamate, or carboxylic acid (incl. sulfur &
 &derivatives).'
        bndty = 302
      else 
        write(ilog,344) ia,ib,' is most likely a C-O bond in an enol or enol ether, most likely a (thio)phenol derivative.'
        bndty = 303
      end if
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is a C-O/S bond in an alkyine alcohol/ether or (thio)cyanate ester.'
      bndty = 403
    else if ((n12(i1).ge.3).AND.(n12(i2).eq.3).AND.(sum(isubs(2,4:5)).ge.2).AND.(roti.gt.0).AND.(isubs(2,7).ge.1)) then
      write(ilog,344) ia,ib,' is most likely a C-S bond in a sulfin-yl/ate or sulfinamide.'
      bndty = 203
    else if ((n12(i1).ge.3).AND.(n12(i2).eq.3).AND.(sum(isubs(2,4:5)).eq.1).AND.(roti.gt.0).AND.(isubs(2,7).ge.1)) then
      write(ilog,344) ia,ib,' is most likely a C-S bond in a sulfoxide, sulfilimine, or sulfinimide.'
      bndty = 204
    else if ((n12(i1).ge.3).AND.(n12(i2).eq.4).AND.(isubs(2,7).gt.0).AND.(isubs(2,5).eq.3).AND.(isubs(2,4).eq.0).AND.&
 &(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a C-S bond in a sulfonate.'
      bndty = 271
    else if ((n12(i1).ge.3).AND.(n12(i2).eq.4).AND.(isubs(2,7).gt.0).AND.(isubs(2,5).eq.2).AND.(isubs(2,4).eq.1).AND.&
 &(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a C-S bond in a sulfonamide or sulfonimidate.'
      bndty = 272
    else if ((n12(i1).ge.3).AND.(n12(i2).eq.4).AND.(isubs(2,7).gt.0).AND.(isubs(2,5).eq.1).AND.(isubs(2,4).eq.1).AND.&
 &(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a C-S bond in a sulfoximine.'
      bndty = 273
    else if ((n12(i1).ge.3).AND.(n12(i2).eq.4).AND.(isubs(2,7).gt.0).AND.(isubs(2,5).eq.1).AND.(isubs(2,4).eq.1).AND.&
 &(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a C-S bond in a sulfone.'
      bndty = 274
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.2).AND.(roti.gt.0)) then
      bndty = 205
      write(ilog,344) ia,ib,' is any linear C-O bond, e.g., to an ester or carbamate, or in an alcohol, ether, thioether, thiol,&
 & sulfenylhalide, etc.'
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.2).AND.(roti.le.0)) then
      bndty = 217
      write(ilog,344) ia,ib,' is a C(al)-O/S bond in a ring (lactam, hemiacetal, epoxide, azolidines, etc).'
    else if ((noops.le.0).AND.(maxval(isubs(1:2,7)).le.0).AND.(n12(i2).gt.2)) then
      bndty = 5
      write(ilog,344) ia,ib,' is an unclassified planar C-S/O bond.'
    else
      write(ilog,344) ia,ib,' is an unclassified C-S/O bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'N  ').AND.(roti.gt.0)) then
    if ((n12(i1).eq.4).AND.(n12(i2).eq.3).AND.(isubs(2,7).gt.0)) then
      write(ilog,344) ia,ib,' is most likely an aliphatic amine carbon-nitrogen linkage.'
      bndty = 209
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.3).AND.(isubs(2,7).le.0).AND.(sum(isubs(2,2:3)).eq.(n12(i2)-1))) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in an N-alkylsubstituted amide, to a heteroaromatic cycle, &
 &in amidine, guanidine, carbamates, urea derivatives, etc.'
      bndty = 210
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.3).AND.(isubs(2,7).le.0).AND.(sum(isubs(2,3:4)).eq.(n12(i2)-1))) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in an N-alkylsubstituted di- or triazole, or in an &
 &N,N-disubstituted hydrazone or hydrazide derivative.'
      bndty = 275
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.3).AND.(isubs(2,7).le.0).AND.(sum(isubs(2,2:4)).eq.(n12(i2)-1))) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in an N-alkylsubstituted hydrazone or hydrazide derivative.'
      bndty = 276
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.2).AND.(sum(isubs(2,2:3)).eq.(n12(i2)-1))) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in an N-alkylsubstituted imine, guanidine, amidine, etc.'
      bndty = 213
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.3).AND.(isubs(2,8).eq.1).AND.(sum(isubs(2,2:3)).eq.(n12(i2)-2))) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in an N-alkylsubstituted phosph(ate/onic/inic) amide, etc.'
      bndty = 236
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.2).AND.(isubs(2,8).eq.1)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in an N-alkylsubstituted iminophosph(or)ane.'
      bndty = 238
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2)) then
      write(ilog,344) ia,ib,' is most likely the C-N in an isonitril or cyanate.'
      bndty = 211
    else if ((n12(i1).eq.2).AND.(n12(i2).ge.3)) then
      write(ilog,344) ia,ib,' is most likely the C-N bond of an a cyanamide.'
      bndty = 212
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.2).AND.(isubs(2,9).eq.1).AND.(isubs(2,13).eq.1).AND.(isubs(2,12).eq.0)) then
      write(ilog,344) ia,ib,' is most likely a C-N bond to a deprotonated sulfonamide or sulfoximine.'
      bndty = 405
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.2).AND.(isubs(1,1).eq.1)) then
      write(ilog,344) ia,ib,' is most likely a C=N or conjugated C-N bond in a an imine, azo compound, etc.'
      bndty = 10
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.2).AND.(isubs(1,1).eq.0)) then
      write(ilog,344) ia,ib,' is most likely a C=N or conjugated C-N bond in a an asymmetric imine, azo compound, etc.'
      bndty = 301
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.&
 &           (isubs(2,20).eq.2)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in a primary amide, aniline, amidine, guanidine, etc.'
      bndty = 6
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.&
 &           ((isubs(1,5)+isubs(1,9)).eq.2)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in a carbamate or a sulfur analog.'
      bndty = 304
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.&
 &           ((isubs(1,5)+isubs(1,9)).eq.1).AND.(isubs(1,4).eq.1)) then
      write(ilog,344) ia,ib,' is most likely one of the single C-N bonds in a (thio)urea derivative.'
      bndty = 8
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.&
 &           (maxval(isubs(1:2,1)).eq.1)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in aniline or a symmetrically substituted amide, amidine, etc.'
      bndty = 7
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.(isubs(2,20).eq.1).AND.&
 &           (isubs(2,3).eq.1).AND.(isubs(1,5).eq.1).AND.(isubs(1,8).eq.0)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in a secondary (singly N-substituted) amide.'
      bndty = 305
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.(isubs(2,20).eq.1).AND.&
 &           (isubs(2,3).eq.1).AND.((isubs(1,3)+isubs(1,2)).eq.2)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in a secondary (singly N-substituted) enamine.'
      bndty = 308
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.(isubs(2,20).eq.1).AND.&
 &           (isubs(2,3).eq.1).AND.(isubs(1,4).eq.0).AND.(isubs(1,8).eq.0)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond in a secondary (singly N-substituted) amidine or guanidine.'
      bndty = 9 ! this one should sort itself out sterically
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(isubs(2,5).eq.2).AND.(maxval(isubs(1:2,7)).le.0)) then
      write(ilog,344) ia,ib,' is most likely the C-N bond in an aromatic nitro compound.'
      bndty = 13
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0).AND.&
 &(isubs(2,4).eq.1).AND.(isubs(1,5).eq.1)) then
      write(ilog,344) ia,ib,' is most likely the C-N bond in a hydrazide.'
      bndty = 311 
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(noops.le.0)) then
      write(ilog,*) ia,ib,' is an unclassified but planar C-N bond between trigonal sp2 centers.'
      bndty = 5
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(isubs(1,3).eq.2).AND.(isubs(1,13).le.0)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond between two trigonal and separately planar but misaligned centers &
 &in an aniline derivative.'
      bndty = 14
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(maxval(isubs(1:2,7)).le.0).AND.(isubs(2,3).eq.2).AND.(isubs(2,13).le.0)) then
      write(ilog,344) ia,ib,' is most likely the single C-N bond between two trigonal and separately planar but misaligned centers &
 &in an N-substituted heteroaromatic ring.'
      bndty = 14
    else if ((noops.gt.0).AND.(maxval(isubs(1:2,7)).le.0)) then
      write(ilog,344) ia,ib,' is a rotatable C-N bond between two trigonal and separately planar but misaligned centers &
 &(aryl amide?).'
      bndty = 216
    else
      write(ilog,344) ia,ib,' is an unclassified, rotatable C-N bond.'
    end if
    return
  else if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'N  ')) then ! many ring bonds taken care of above
    if ((n12(i1).eq.4).AND.(n12(i2).eq.2)) then
      write(ilog,344) ia,ib,' is most likely a single =N-C bond in a mixed ring.'
    else
      write(ilog,344) ia,ib,' is an unclassified, C-N or C=N bond in a ring.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'S  ').AND.(atnam(i2).eq.'S  ')) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(roti.gt.0)) then
      bndty = 218
      write(ilog,344) ia,ib,' is most likely an S-S bond in a linear disulfide.' ! risky as steric information is limited
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(roti.le.0)) then
      bndty = 219
      write(ilog,344) ia,ib,' is most likely an S-S bond in a cyclic disulfide.'
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable S-S bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'O  ').AND.(atnam(i2).eq.'O  ')) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(roti.gt.0)) then
      bndty = 218
      write(ilog,344) ia,ib,' is most likely an O-O bond in a linear peroxide.' ! risky as steric information is limited
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(roti.le.0)) then
      bndty = 219
      write(ilog,344) ia,ib,' is most likely an O-O bond in a cyclic disulfide.'
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable O-O bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'O  ').AND.(atnam(i2).eq.'S  ')) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.(isubs(2,5).ge.2)) then
      if (roti.gt.0) then
        bndty = 220
        write(ilog,344) ia,ib,' is most likely an S-O bond in a linear sulfonic ester, sulfuric amide, or sulfate.'
      else
        bndty = 221
        write(ilog,344) ia,ib,' is most likely an S-O bond in a cyclic sulfonic ester (sultone), sulfuric amide, or sulfate.'
      end if
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.3).AND.(isubs(2,5).eq.1)) then
      if (roti.gt.0) then
        bndty = 223
        write(ilog,344) ia,ib,' is most likely an S-O bond in a linear sulfinic ester.'
      else
        bndty = 224
        write(ilog,344) ia,ib,' is most likely an S-O bond in a cyclic sulfinic ester (sultine).'
      end if
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable O-S bond.' 
      bndty = 5
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable O-S bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'O  ').AND.(atnam(i2).eq.'P  ')) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.(isubs(2,11).le.0).AND.(isubs(2,7).ge.1)) then
      bndty = 240
      write(ilog,344) ia,ib,' is most likely a P-O bond in a phosphonium ion.'
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.((isubs(2,5)+isubs(2,9)).eq.3)) then
      if (roti.gt.0) then
        bndty = 225
        write(ilog,344) ia,ib,' is most likely a P-O bond in a linear (thio)phosphate (ester).'
      else
        bndty = 226
        write(ilog,344) ia,ib,' is most likely a P-O bond in a cyclic (thio)phosphodiester.'
      end if
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.1).AND.((isubs(2,2)+isubs(2,3)).eq.2)) then
      if (roti.gt.0) then
        bndty = 227
        write(ilog,344) ia,ib,' is most likely a P-O bond in a linear phosphinic ester.'
      else
        bndty = 228
        write(ilog,344) ia,ib,' is most likely a P-O bond in a cyclic phosphinic ester.'
      end if
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.2)) then
      if (roti.gt.0) then
        bndty = 229
        write(ilog,344) ia,ib,' is most likely a P-O bond in a linear phosphonate or phosphoric (mono)amide.'
      else
        bndty = 230
        write(ilog,344) ia,ib,' is most likely a P-O bond in a cyclic phosphonate or phosphoric (mono)amide.'
      end if
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable P-O bond.' 
      bndty = 16
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable P-O bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'P  ').AND.(atnam(i2).eq.'S  ')) then
    if ((n12(i2).eq.2).AND.(n12(i1).eq.4).AND.(isubs(1,11).le.0).AND.(isubs(1,7).ge.1)) then
      bndty = 259
      write(ilog,344) ia,ib,' is most likely a P-S bond in a phosphonium ion.'
    else if ((n12(i2).eq.2).AND.(n12(i1).eq.4).AND.((isubs(1,5)+isubs(1,9)).eq.3)) then
      if (roti.gt.0) then
        bndty = 260
        write(ilog,344) ia,ib,' is most likely a P-S bond in a linear (thio)phosphate (ester).'
      else
        bndty = 261
        write(ilog,344) ia,ib,' is most likely a P-S bond in a cyclic (thio)phosphodiester.'
      end if
    else if ((n12(i2).eq.2).AND.(n12(i1).eq.4).AND.((isubs(1,5)+isubs(1,9)).eq.1).AND.((isubs(1,2)+isubs(1,3)).eq.2)) then
      if (roti.gt.0) then
        bndty = 262
        write(ilog,344) ia,ib,' is most likely a P-S bond in a linear (thio)phosphinic ester.'
      else
        bndty = 263
        write(ilog,344) ia,ib,' is most likely a P-S bond in a cyclic (thio)phosphinic ester.'
      end if
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.2)) then
      if (roti.gt.0) then
        bndty = 264
        write(ilog,344) ia,ib,' is most likely a P-S bond in a linear (thio)phosphonate or (thio)phosphoric (mono)amide.'
      else
        bndty = 265
        write(ilog,344) ia,ib,' is most likely a P-S bond in a cyclic (thio)phosphonate or (thio)phosphoric (mono)amide.'
      end if
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable P-S bond.' 
      bndty = 16
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable P-S bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'N  ').AND.(atnam(i2).eq.'P  ')) then
    if ((n12(i1).eq.3).AND.(n12(i2).eq.4).AND.(isubs(2,11).le.0).AND.(isubs(2,7).ge.1)) then
      bndty = 239
      write(ilog,344) ia,ib,' is most likely a P-N bond in a phosphonium ion.'
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.4).AND.(isubs(2,4).eq.3)) then
      if (roti.gt.0) then
        bndty = 231
        write(ilog,344) ia,ib,' is most likely a P-N bond in a linear phosphazene.'
      else
        bndty = 232
        write(ilog,344) ia,ib,' is most likely a P-N bond in a cyclic phosphazene.'
      end if
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.4).AND.(isubs(2,4).eq.3).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a P=N bond in an iminophosphorane (phosphazene).'
      bndty = 5
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.3).AND.(isubs(2,4).eq.3).AND.(isubs(2,7).eq.0).AND.(noops.le.0).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a P=N bond at a trigonal-planar P atom.'
      bndty = 5
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(isubs(2,4).eq.3).AND.(isubs(2,7).ge.1).AND.(noops.gt.0).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a P-N bond in an aminophosphane.'
      bndty = 237
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.3).AND.(roti.gt.0)) then
      if (isubs(1,7).eq.0) then
        bndty = 310
      else
        bndty = 233
      end if
      write(ilog,344) ia,ib,' is most likely a P-N bond in a phosphoric (mono)amide.'
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.2).AND.(roti.gt.0)) then
      if (isubs(1,7).eq.0) then
        bndty = 310
      else
        bndty = 234
      end if
      write(ilog,344) ia,ib,' is most likely a P-N bond in a phosphonic (mono)amide or phosphodiamide.'
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.1).AND.((isubs(2,2)+isubs(2,3)).eq.2).AND.(roti.gt.0)) then
      if (isubs(1,7).eq.0) then
        bndty = 309
      else
        bndty = 235
      end if
      write(ilog,344) ia,ib,' is most likely a P-N bond in a phosphinamide.'
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(noops.le.0)) then
      write(ilog,344) ia,ib,' is most likely a linear P=N or conjugated P-N bond (conj. amino- or iminophosph(a/i)ne, etc).'
      bndty = 306
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable P-N bond.' 
      bndty = 5
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable P-N bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'N  ').AND.(atnam(i2).eq.'S  ')) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.(roti.gt.0)) then
      bndty = 266
      write(ilog,344) ia,ib,' is most likely an S-N bond in a sulfenylimine or S-nitrosothiol.'
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.4)) then
      if (roti.le.0) then
        bndty = 267
        write(ilog,344) ia,ib,' is most likely an S-N= bond in a cyclic sulfonamide, and S-N- bond in a cyclic sulfonamide anion, &
 &or an S=N bond in a cyclic sulfoximine or related compound.'
      else if ((roti.gt.0).AND.(isubs(2,11).eq.1).AND.(isubs(2,5).eq.1)) then
        bndty = 16
        write(ilog,344) ia,ib,' is most likely an S=N bond in a sulfoximine.'
      else if ((roti.gt.0).AND.(isubs(2,3).eq.1).AND.(isubs(2,5).eq.2).AND.(isubs(1,12).ge.1)) then
        bndty = 404
        write(ilog,344) ia,ib,' is most likely an S-N bond in a deprotonated sulfon(i/a)mide or an S=N bond in a &
 &sulfonylimidamide.'
      else if (roti.gt.0) then
        bndty = 5
        write(ilog,344) ia,ib,' is an unclassified S-N or S=N bond.'
      else
        write(ilog,344) ia,ib,' is an unclassified S-N or S=N bond in a ring.'
      end if
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.2).AND.(isubs(1,7).gt.0).AND.(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely an S-N bond in a sulfenamide.'
      bndty = 270
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.3).AND.(isubs(1,7).eq.0).AND.(isubs(2,7).ge.1).AND.(roti.gt.0).AND.&
 &(isubs(2,5).eq.1)) then
      write(ilog,344) ia,ib,' is most likely an S-N bond in a sulfinamide.'
      bndty = 268
    else if ((n12(i1).eq.3).AND.(n12(i2).eq.4).AND.(isubs(1,7).eq.0).AND.(isubs(2,7).ge.1).AND.(roti.gt.0).AND.&
 &(isubs(2,5).eq.2)) then
      write(ilog,344) ia,ib,' is most likely an S-N bond in a sulfonamide.'
      bndty = 269
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.3).AND.(sum(isubs(2,4:5)).eq.3).AND.(isubs(2,7).eq.0).AND.(noops.le.0).AND.& 
 &(roti.gt.0)) then
      write(ilog,344) ia,ib,' is most likely a S=N bond at a trigonal-planar S atom (e.g., S,S-dioxide).'
      bndty = 5
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable S-N bond.' 
      bndty = 5
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable S-N bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'C  ').AND.(atnam(i2).eq.'P  ')) then
    if ((n12(i1).eq.4).AND.(n12(i2).eq.4).AND.(isubs(2,11).le.0).AND.(isubs(2,7).ge.1)) then
      bndty = 241
      write(ilog,344) ia,ib,' is most likely a P-C bond in a phosphonium ion.'
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.3)) then
      bndty = 242
      write(ilog,344) ia,ib,' is most likely a P-C bond in an aliphatic phosphonate.'
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.2)) then
      bndty = 243
      write(ilog,344) ia,ib,' is most likely a P-C bond in an aliphatic phosphinate or phosphonic (mono)amide.'
    else if ((n12(i1).eq.4).AND.(n12(i2).eq.4).AND.(isubs(2,5).eq.1).AND.((isubs(2,2)+isubs(2,3)).eq.2)) then
      bndty = 244
      write(ilog,344) ia,ib,' is most likely a P-C bond in an aliphatic phosphine oxide.'
    else if (roti.gt.0) then
      write(ilog,344) ia,ib,' is an unclassified, rotatable C-P bond.'
    end if
    return
  end if
!
  if ((atnam(i1).eq.'N  ').AND.(atnam(i2).eq.'O  ').AND.(roti.gt.0)) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.2)) then
      write(ilog,344) ia,ib,' is most likely the N-O bond in an oxime or oxime derivative or a nitrite ester.'
      bndty = 306
    else if (isubs(1,5).eq.2) then
      write(ilog,344) ia,ib,' is most likely a nitrate ester.' ! somewhat oddly, these are planar; doi: 10.1139/v72-273
      bndty = 11
    else if ((isubs(1,7).eq.0).AND.(n12(i2).eq.2).AND.(noops.le.0)) then
      write(ilog,344) ia,ib,' is a rotatable, planar bond between an sp2 nitrogen and oxygen.'
      bndty = 307
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable N-O bond.' 
      bndty = 5
    else
      write(ilog,344) ia,ib,' is an unclassified and rotatable N-O bond.' 
      bndty = -1
    end if
    return
  end if
!
  if ((atnam(i1).eq.'N  ').AND.(atnam(i2).eq.'N  ').AND.(roti.gt.0)) then
    if ((n12(i1).eq.2).AND.(n12(i2).eq.2).AND.&
 &      (((isubs(1,4).eq.1).AND.(isubs(1,10).eq.1)).OR.((isubs(2,4).eq.1).AND.(isubs(2,10).eq.1)))) then
      write(ilog,344) ia,ib,' is most likely an azide N=N bond.'
      bndty = 214
    else if ((n12(i1).eq.2).AND.(n12(i2).eq.2)) then
      write(ilog,344) ia,ib,' is most likely a linear N=N or conjugated N-N bond (azo compound, etc).'
      bndty = 306
    else if ((n12(i1).ge.3).AND.(n12(i2).ge.3).AND.(minval(isubs(1:2,7)).ge.1)) then
      write(ilog,344) ia,ib,' is most likely an N-N bond in a hydrazine derivative.'
      bndty = 215 ! could lead to problems because cis is strongly disfavored electronically (lone pairs)
    else if (((n12(i1)+n12(i2)).eq.5).AND.(maxval(isubs(1:2,7)).eq.0).AND.(noops.le.0)) then
      write(ilog,344) ia,ib,' is most likely an N-N bond in a hydrazone or nitrosamine derivative.'
      bndty = 307 ! asymmetric only if with N-NHR
    else if (noops.le.0) then
      write(ilog,344) ia,ib,' is an unclassified but planar and rotatable N-N bond.' 
      bndty = 5
    else if ((n12(i1).ge.3).AND.(n12(i2).ge.3).AND.(minval(isubs(1:2,7)).eq.0)) then
      write(ilog,344) ia,ib,' is most likely an N-N bond in a hydrazide.'
      bndty = 253 ! could lead to problems because of lone pairs (see above)
    else
      write(ilog,344) ia,ib,' is an unclassified and rotatable N-N bond.' 
      bndty = -1
    end if
    return
  end if
!
  write(ilog,344) ia,ib,' is an unclassified bond CAMPARI knows absolutely nothing about.'
!
end
!
!-------------------------------------------------------------------------------------------------------
!
subroutine guess_torpot(i1,i2,i3,i4,rs,i,bndty)
!
  use inter
  use atoms
  use math
!
  implicit none
!
  integer rs,i,i1,i2,i3,i4,bndty
  RTYPE getztor
! 
  iaa(rs)%par_di(i,:) = 0.0
  if ((bndty.eq.1).OR.(bndty.eq.4).OR.((bndty.ge.6).AND.(bndty.le.12))) then ! twofold term 0,180
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 5.0
    iaa(rs)%par_di(i,3) = -5.0
  else if (bndty.eq.2) then ! twofold term 0,180
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 10.0
    iaa(rs)%par_di(i,3) = -10.0
  else if (bndty.eq.3) then ! twofold term 0,180
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 20.0
    iaa(rs)%par_di(i,3) = -20.0
  else if (bndty.eq.5) then ! unclassified planar -> maintain
    iaa(rs)%typ_di(i) = 2
    iaa(rs)%par_di(i,2) = getztor(i1,i2,i3,i4)
    iaa(rs)%par_di(i,1) = 10.0/(RADIAN*RADIAN)
  else if (bndty.eq.13) then ! twofold term 0,180
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 2.5
    iaa(rs)%par_di(i,3) = -2.5
  else if (bndty.eq.14) then ! twofold term 0,180
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 1.25
    iaa(rs)%par_di(i,3) = -1.25
  else if (bndty.eq.15) then ! twofold term 0,180
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 0.5
    iaa(rs)%par_di(i,3) = -0.5
  else if (bndty.eq.16) then ! unclassified planar -> maintain
    iaa(rs)%typ_di(i) = 2
    iaa(rs)%par_di(i,2) = getztor(i1,i2,i3,i4)
    if (abs(iaa(rs)%par_di(i,2)).le.20.0) iaa(rs)%par_di(i,2) = 0.0
    if (abs(iaa(rs)%par_di(i,2)).ge.160.0) iaa(rs)%par_di(i,2) = 180.0
    iaa(rs)%par_di(i,1) = 10.0/(RADIAN*RADIAN)
  else if (bndty.eq.301) then ! heavy atoms on opposite sides
    iaa(rs)%typ_di(i) = 1
    if ((atnam(i1).ne.'H  ').AND.(atnam(i4).ne.'H  ')) then
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = 0.5
      iaa(rs)%par_di(i,3) = -4.5 ! favors E-isomer by 1kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.302) then ! substituent on side of double-bonded O/S
    iaa(rs)%typ_di(i) = 1
    if (((atnam(i1).eq.'O  ').AND.(n12(i1).eq.1).AND.(n12(i2).eq.3)).OR.&
 &      ((atnam(i4).eq.'O  ').AND.(n12(i4).eq.1).AND.(n12(i3).eq.3))) then
      iaa(rs)%par_di(i,1) = 8.0
      iaa(rs)%par_di(i,2) = -1.5
      iaa(rs)%par_di(i,3) = -6.5 ! favors Z-isomer by 3kcal/mol
    else if (((atnam(i1).eq.'S  ').AND.(n12(i1).eq.1).AND.(n12(i2).eq.3)).OR.&
 &           ((atnam(i4).eq.'S  ').AND.(n12(i4).eq.1).AND.(n12(i3).eq.3))) then
      iaa(rs)%par_di(i,1) = 8.0
      iaa(rs)%par_di(i,2) = -1.0
      iaa(rs)%par_di(i,3) = -7.0 ! favors Z-isomer by 2kcal/mol
    else
      iaa(rs)%par_di(i,1) = 8.0
      iaa(rs)%par_di(i,3) = -8.0
    end if
  else if (bndty.eq.303) then ! double bond pointing away from lone pairs; doi: 10.1039/TF9646000634
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if system is conjugated
    if ((((atnam(i1).eq.'H  ').OR.(n12(i1).eq.4)).AND.(n12(i2).eq.3)).OR.&
 &      (((atnam(i4).eq.'H  ').OR.(n12(i4).eq.4)).AND.(n12(i3).eq.3))) then 
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = 0.5
      iaa(rs)%par_di(i,3) = -4.5 ! favors Z-isomer by 1kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.304) then ! hydrogen toward singly bonded O/S; doi:10.1021/ct1004017
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if system is conjugated
    if (((atnam(i1).eq.'H  ').AND.((atnam(i4).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &      ((atnam(i4).eq.'H  ').AND.((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then 
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = 0.75
      iaa(rs)%par_di(i,3) = -4.25 ! favors E-isomer by 1.5kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.305) then
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if trigonal atoms is not bound to H
    if ((atnam(i1).eq.'C  ').AND.(atnam(i4).eq.'C  ')) then 
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = 0.75
      iaa(rs)%par_di(i,3) = -4.25 ! favors Z-isomer by 1.5kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.306) then ! always trans (only a single term acting)
    iaa(rs)%typ_di(i) = 1
    iaa(rs)%par_di(i,1) = 20.0
    iaa(rs)%par_di(i,2) = 2.0
    iaa(rs)%par_di(i,3) = -18.0 ! favors E-isomer by 4.0kcal/mol
  else if (bndty.eq.307) then
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if trigonal atoms is not bound to H
    if ((atnam(i1).eq.'H  ').OR.(atnam(i4).eq.'H  ')) then 
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = -1.0
      iaa(rs)%par_di(i,3) = -4.0 ! favors E-isomer by 2kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.308) then ! hydrogen toward doubly bonded C; doi:10.1021/ct1004017
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if system is conjugated
    if (((atnam(i1).eq.'H  ').AND.(atnam(i4).eq.'C  ').AND.(n12(i4).eq.3)).OR.&
 &      ((atnam(i4).eq.'H  ').AND.(atnam(i1).eq.'C  ').AND.(n12(i1).eq.3))) then 
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = -0.75
      iaa(rs)%par_di(i,3) = -4.25 ! favors E-isomer by 1.5kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.309) then ! hydrogen toward doubly bonded O/S
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if the system is doubly substituted
    if (((atnam(i1).eq.'H  ').AND.((atnam(i4).eq.'O  ').OR.(atnam(i4).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &      ((atnam(i4).eq.'H  ').AND.((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then 
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,2) = -0.75
      iaa(rs)%par_di(i,3) = -4.25 ! favors E-isomer by 1.5kcal/mol
    else if ((((atnam(i4).eq.'O  ').OR.(atnam(i4).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &           (((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    else
      iaa(rs)%typ_di(i) = 0
    end if
  else if (bndty.eq.310) then ! hydrogen toward doubly bonded O/S
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if the system is doubly substituted
    if (((atnam(i1).eq.'H  ').AND.((atnam(i4).eq.'O  ').OR.(atnam(i4).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &      ((atnam(i4).eq.'H  ').AND.((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then 
      iaa(rs)%par_di(i,1) = 1.5
      iaa(rs)%par_di(i,2) = -0.75
      iaa(rs)%par_di(i,3) = -0.75 ! favors E-isomer by 1.5kcal/mol
    else if ((((atnam(i4).eq.'O  ').OR.(atnam(i4).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &           (((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then
      iaa(rs)%par_di(i,1) = 1.5
      iaa(rs)%par_di(i,3) = -1.5
    else
      iaa(rs)%typ_di(i) = 0
    end if
  else if (bndty.eq.311) then ! O and N on the same side 
    iaa(rs)%typ_di(i) = 1
    if (((atnam(i1).eq.'O  ').AND.(atnam(i4).eq.'N  ')).OR.((atnam(i1).eq.'N  ').AND.(atnam(i4).eq.'O  '))) then 
      iaa(rs)%par_di(i,1) = 5.
      iaa(rs)%par_di(i,2) = -0.5
      iaa(rs)%par_di(i,3) = -4.5 ! favors that isomer by 1kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  else if (bndty.eq.312) then ! Os opposite 
    iaa(rs)%typ_di(i) = 1
    if ((atnam(i1).eq.'O  ').AND.(atnam(i4).eq.'O  ')) then 
      iaa(rs)%par_di(i,1) = 2.5
      iaa(rs)%par_di(i,2) = 1.
      iaa(rs)%par_di(i,3) = -1.5 ! favors that isomer by 3kcal/mol
    else
      iaa(rs)%par_di(i,1) = 2.5
      iaa(rs)%par_di(i,3) = -2.5
    end if
  else if (bndty.eq.313) then ! Os opposite 
    iaa(rs)%typ_di(i) = 1
    if ((atnam(i1).eq.'S  ').AND.(atnam(i4).eq.'S  ').AND.(n12(i1).eq.1).AND.(n12(i2).eq.1)) then 
      iaa(rs)%par_di(i,1) = 2.5
      iaa(rs)%par_di(i,2) = 1.
      iaa(rs)%par_di(i,3) = -1.5 ! favors that isomer by 3kcal/mol
    else
      iaa(rs)%par_di(i,1) = 2.5
      iaa(rs)%par_di(i,3) = -2.5
    end if
  else if (bndty.eq.314) then ! hydrogen toward doubly bonded O/S
    iaa(rs)%typ_di(i) = 1
!   this will not trigger if the system is doubly substituted
    if (((atnam(i1).eq.'H  ').AND.((atnam(i4).eq.'O  ').OR.(atnam(i4).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &      ((atnam(i4).eq.'H  ').AND.((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then 
      iaa(rs)%par_di(i,1) = 3.5
      iaa(rs)%par_di(i,2) = -1.0
      iaa(rs)%par_di(i,3) = -2.5 ! favors E-isomer by 2kcal/mol
    else if ((((atnam(i4).eq.'O  ').OR.(atnam(i4).eq.'S  ')).AND.(n12(i4).eq.1)).OR.&
 &           (((atnam(i1).eq.'O  ').OR.(atnam(i1).eq.'S  ')).AND.(n12(i1).eq.1))) then
      iaa(rs)%par_di(i,1) = 3.5
      iaa(rs)%par_di(i,3) = -3.5
    else
      iaa(rs)%typ_di(i) = 0
    end if
  else if (bndty.eq.315) then ! Ns opposite 
    iaa(rs)%typ_di(i) = 1
    if ((atnam(i1).eq.'N  ').AND.(atnam(i4).eq.'N  ')) then 
      iaa(rs)%par_di(i,1) = 50
      iaa(rs)%par_di(i,2) = 3.5
      iaa(rs)%par_di(i,3) = -1.5 ! favors that isomer by 3kcal/mol
    else
      iaa(rs)%par_di(i,1) = 5.0
      iaa(rs)%par_di(i,3) = -5.0
    end if
  end if
!
end
!
!----------------------------------------------------------------------------------------------------------
!

