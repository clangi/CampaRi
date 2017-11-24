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
! CONTRIBUTIONS: Adam Steffen                                              !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!------------------------------------------------------------------------
!
subroutine FMSMC_pdb(ipdb,ndump)
!
  use iounit
  use polypep
  use sequen
  use system
  use atoms
  use params
  use aminos
  use molecule
  use pdb
  use mcsums
  use grandensembles
  use cutoffs, ONLY: rsw1
!
  implicit none
!
  integer i,ii,rs,k,ipdb,atti,nrch,imol,ndump,kk,kkk,mt,rsshft,rsshft2
  integer attiwr,rswr,nunk,backupval,ncyc,nins,lnins
  character(3) resname,resname2,resnameprt
  character(4) ana
  character(1) chch,dum
  character(1) extra,extra2
  character(6) tstr
  logical ismember,fixrnprt
  RTYPE occup,shf,zero,ninety,outv(3),shf3(3)
  character(MAXSTRLEN+50) fpdbstr
!
 30    format (a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,1x,f5.2,'  0.00') ! default PDB format, overridden
 31    format ('CRYST1',1x,3(f8.3,1x),' 90.00  90.00  90.00 P 1&
 &           1')
 29    format ('CRYST1',1x,3(f8.3,1x),3(f6.2,1x),' P 1           1')
 32    format ('MODEL ',i8)
 33    format ('ENDMDL')
 34    format ('END')
 35    format ('TITLE',4x,a19,': step = ',i10)
!
  backupval = rsw1
  rsw1 = -1
  extra2 = ' ' ! do not reproduce insertion code
!
  call strlims(pdb_formstr(1),i,ii)
  fpdbstr='('//pdb_formstr(1)(i:ii)//",1x,f5.2,'  0.00')"
  zero = 0.0
  ninety = 90.0
  if (((pdb_force_box.gt.1).OR.(pdb_rmol.gt.0)).AND.(bnd_type.eq.1).AND.(bnd_shape.ne.1).AND.(bnd_shape.ne.3)) then
    write(ilog,*) 'Fatal. Encountered unsupported box shape for image adjustment &
 &in FMCSC_pdb(...). This is an omission bug (please report).'
    call fexit()
  end if
  if (bnd_shape.eq.1) then
    write(ipdb,31) bnd_params(1),bnd_params(2),bnd_params(3)
  else if (bnd_shape.eq.2) then
    write(ipdb,31) 2.0*bnd_params(4),2.0*bnd_params(4),2.0*bnd_params(4)!,&
! &         bnd_params(1)-bnd_params(4),bnd_params(2)-bnd_params(4),bnd_params(3)-bnd_params(4)
  else if (bnd_shape.eq.3) then
    write(ipdb,31) 2.0*bnd_params(4),2.0*bnd_params(4),bnd_params(6)!,&
! &         bnd_params(1)-bnd_params(4),bnd_params(2)-bnd_params(4),bnd_params(3)-0.5*bnd_params(6)
  end if
  if (pdb_writemode.eq.2) then
    write(ipdb,35) basename(1:19),nstep
  end if
  if ((pdb_writemode.eq.2).AND.(ndump.ne.-1)) then
    write(ipdb,32) ndump
  end if
  nrch = 0
  atti = 0
  ana = '    '
  rsshft = 0
  shf = 0.0
  nunk = 0
  nins = 0
  lnins = 0
  do imol=1,nmol
   if ((use_trajidx.EQV..false.).AND.(just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
   if ((rsmol(imol,2)-rsmol(imol,1)).gt.0) then
     nrch = nrch + 1
     if (nrch.eq.1) chch='A'
     if (nrch.eq.2) chch='B'
     if (nrch.eq.3) chch='C'
     if (nrch.eq.4) chch='D'
     if (nrch.eq.5) chch='E'
     if (nrch.eq.6) chch='F'
     if (nrch.eq.7) chch='G'
     if (nrch.eq.8) chch='H'
     if (nrch.eq.9) chch='I'
     if (nrch.eq.10) chch='J'
     if (nrch.eq.11) chch='K'
     if (nrch.eq.12) chch='L'
     if (nrch.eq.13) chch='M'
     if (nrch.eq.14) chch='N'
     if (nrch.eq.15) chch='O'
     if (nrch.eq.16) chch='P'
     if (nrch.eq.17) chch='Q'
     if (nrch.eq.18) chch='R'
     if (nrch.eq.19) chch='S'
     if (nrch.eq.20) chch='T'
     if (nrch.eq.21) chch='U'
     if (nrch.eq.22) chch='V'
     if (nrch.eq.23) chch='W'
     if (nrch.eq.24) chch='X'
     if (nrch.eq.25) chch='Y'
     if (nrch.eq.26) chch='Z'
   else
     chch=' '
   end if
   mt = moltypid(imol)
   if (pdb_rmol.gt.0) then
     call shift_bound3(pdb_rmol,imol,shf3)
   else
     shf3(:) = 0.0
   end if
   if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
     if (ismember(fluctypes,mt).EQV..false.) then
       occup = 1.00
       shf = 0.0
     else
       if (ismember(ispresent,imol).EQV..true.) then
         occup = 1.00
         shf = 0.0
       else
         occup = 0.00
         if (bnd_shape.eq.1) then
           shf = 2.0*bnd_params(3)
         else if (bnd_shape.eq.2) then
           shf = 4.0*bnd_params(4)
         else if (bnd_shape.eq.3) then
           shf = 2.0*bnd_params(6)
         else
           call fexit()
         end if
       end if
     end if
     shf3(3) = shf3(3) + shf
   else
     occup = 1.00
   end if
   do rs=rsmol(imol,1),rsmol(imol,2)
    ncyc = 0
    lnins = nins
    nins = 0
    if (seqtyp(rs).ne.26) then ! not UNK
      resname = amino(seqtyp(rs))
      resnameprt = resname
    else
      nunk = nunk + 1
      resname = pdb_unknowns(nunk)
      resnameprt = resname
    end if
    fixrnprt = .false.
    rsshft2 = 0
    if (at(rs)%nbb.le.0) then
      write(ilog,*) 'Fatal error in prtpdb(...). Builder generated r&
 &esidue without any backbone atom. This is a bug.'
      call fexit()
    end if
    i = at(rs)%bb(1)
    do ii=i,i+at(rs)%nbb+at(rs)%nsc-1
      if (use_trajidx.EQV..true.) then
        if (pdboutlog(ii).EQV..false.) cycle
      end if
!      if (i.le.at(rs)%nbb) then
!        ii = at(rs)%bb(i)
!      else
!        ii = at(rs)%sc(i-at(rs)%nbb)
!      end if
      if (b_type(ii).le.n_biotyp) then
!       if we write out nucleotides in the exp. format, we need to omit the O3P at the beginning
!       of this non-N-terminal polynucleotide residue
        if (pdb_nucmode.eq.2) then
          if (seqpolty(rs).eq.'N') then
            if ((rs.gt.rsmol(molofrs(rs),1)).AND.(lnins.eq.1)) then
              if ((b_type(ii).eq.726).OR.(b_type(ii).eq.673).OR.&
 &                (b_type(ii).eq.780).OR.(b_type(ii).eq.836).OR.&
 &                (b_type(ii).eq.894).OR.(b_type(ii).eq.700).OR.&
 &                (b_type(ii).eq.647).OR.(b_type(ii).eq.753).OR.&
 &                (b_type(ii).eq.808).OR.(b_type(ii).eq.865).OR.&
 &                (b_type(ii).eq.626).OR.(b_type(ii).eq.938)) then
                if ((ii.ne.at(rs)%bb(1)).OR.(ii.ne.nuci(rs,1))) then
                  write(ilog,*) 'Warning. Re-assignment of oxygen at&
 &om to residue i-1 to comply with PDB-standard might be corrupted.'
                end if
                ncyc = ncyc + 1
                if (ncyc.eq.1) cycle
              else if ((seqtyp(rs).eq.26).AND.((bio_code(b_type(ii)).eq.'O3*').OR.(bio_code(b_type(ii)).eq."O3'"))) then
                ncyc = ncyc + 1
                if (ncyc.eq.1) cycle      
              end if
            end if
          end if
        end if
!
        atti = atti + 1
!
!       manually add numbering to non-distuingishable heavy atom types (same biotype)
!       overwrite some internal names for full pdb-compatibility
        ana(1:1) = ' '
        ana(2:4) = bio_code(b_type(ii))
        tstr(1:6) = 'ATOM  '
        extra = ' '
        call name_heavies(ii,rs,ana,ncyc)
!
!       now the same for hydrogens
        call name_hydrogens(ii,rs,imol,ana)
!
!       now the same for hydrogens
        if (seqtyp(rs).eq.26) ana(1:1) = atnam(ii)(3:3)
!
!       change gross formatting
        if (pdb_convention(1).eq.3) then ! CHARMM
          if ((ana(1:1).eq.'1').OR.(ana(1:1).eq.'2').OR.(ana(1:1).eq.'3')) then
            if (ana(4:4).ne.' ') then
              dum = ana(1:1)
              ana(1:3) = ana(2:4)
              ana(4:4) = dum
            else
              if (ana(3:3).eq.' ') then
                ana(3:3) = ana(1:1)
                ana(1:1) = ' '
              else
                ana(4:4) = ana(1:1)
                ana(1:1) = ' '
              end if
            end if
          end if
          if ((seqtyp(rs).eq.27).OR.(seqtyp(rs).eq.28)) then ! N-terminal cap residue is merged into next one
            if (seqpolty(rs+1).eq.'P') then ! only merge into valid peptide residues
              if (ii.eq.i) rsshft = rsshft + 1
              rsshft2 = 1
              resnameprt = amino(seqtyp(rs+1))
            end if
          else if ((seqtyp(rs).eq.29).OR.(seqtyp(rs).eq.30)) then ! C-terminal cap res. is merged into one before
            if (seqpolty(rs-1).eq.'P') then ! only merge into valid peptide residues
              if (ii.eq.i) rsshft = rsshft + 1
              resnameprt = amino(seqtyp(rs-1))
            end if
          end if
          if (resname.eq.'HIE') resnameprt = 'HSE'
          if (resname.eq.'HID') resnameprt = 'HSD'
          if (resname.eq.'HIP') resnameprt = 'HSP'
          if ((resname.eq.'DPU').OR.(resname.eq.'DIU').OR.(resname.eq.'RPU').OR.(resname.eq.'RIU')) resnameprt = 'URA'
          if ((resname.eq.'DPG').OR.(resname.eq.'DIG').OR.(resname.eq.'RPG').OR.(resname.eq.'RIG')) resnameprt = 'GUA'
          if ((resname.eq.'DPA').OR.(resname.eq.'DIA').OR.(resname.eq.'RPA').OR.(resname.eq.'RIA')) resnameprt = 'ADE'
          if ((resname.eq.'DPC').OR.(resname.eq.'DIC').OR.(resname.eq.'RPC').OR.(resname.eq.'RIC')) resnameprt = 'CYT'
          if ((resname.eq.'DPT').OR.(resname.eq.'DIT').OR.(resname.eq.'RPT').OR.(resname.eq.'RIT')) resnameprt = 'THY'
        else if (pdb_convention(1).eq.2) then ! GROMOS
          if (resname.eq.'NME') resnameprt = 'NAC'
          if ((ana(1:1).eq.'1').OR.(ana(1:1).eq.'2').OR.(ana(1:1).eq.'3')) then
            if (ana(4:4).ne.' ') then
              dum = ana(1:1)
              ana(1:3) = ana(2:4)
              ana(4:4) = dum
            else
              if (ana(3:3).eq.' ') then
                ana(3:3) = ana(1:1)
                ana(1:1) = ' '
              else
                ana(4:4) = ana(1:1)
                ana(1:1) = ' '
              end if
            end if
          end if
          if ((ii.eq.hni(rs)).AND.(seqpolty(rs).eq.'P')) ana(3:3) = ' '
        else if (pdb_convention(1).eq.4) then ! AMBER
          if (resname.eq.'T4P') resnameprt = 'TP4'
          if (resname.eq.'T3P') resnameprt = 'TP3'
          if (resname.eq.'NA+') resnameprt = 'Na+'
          if (resname.eq.'CL-') resnameprt = 'Cl-'
          if (resname.eq.'T5P') resnameprt = 'TP5'
          if ((ana(1:1).eq.'1').OR.(ana(1:1).eq.'2').OR.(ana(1:1).eq.'3')) then
            if (ana(4:4).ne.' ') then
              dum = ana(1:1)
              ana(1:3) = ana(2:4)
              ana(4:4) = dum
            else
              if (ana(3:3).eq.' ') then
                ana(3:3) = ana(1:1)
                ana(1:1) = ' '
              else
                ana(4:4) = ana(1:1)
                ana(1:1) = ' '
              end if
            end if
          end if
          if ((seqpolty(rs).eq.'N').AND.(fixrnprt.EQV..false.)) then
            if ((seqtyp(rs).ne.26).AND.(seqflag(rs).eq.24)) then ! DIB, RIB, RIX, or DIX
              if (rs.eq.rsmol(imol,2)) then
                resnameprt(2:2) = resnameprt(3:3)
                resnameprt(3:3) = 'N'
              else
                resnameprt(2:2) = resnameprt(3:3)
                resnameprt(3:3) = '5'
              end if
            else if (rs.eq.rsmol(imol,2)) then
              resnameprt(2:2) = resnameprt(3:3)
              resnameprt(3:3) = '3'
            else if (rs.gt.rsmol(imol,1)) then ! keep CAMPARI nomenclature for 5'-P residues
              resnameprt(2:2) = resnameprt(1:1)
              resnameprt(1:1) = ' '
            end if
            fixrnprt = .true.
          end if
        end if
!
!       override the ATOM-string for hetero atoms
        if (seqtyp(rs).eq.26) then
          if (rsmol(imol,1).eq.rsmol(imol,2)) tstr(1:6) = 'HETATM'
        else
          if ((b_type(ii).gt.443).AND.(.NOT.&
 &          ((b_type(ii).ge.519).AND.(b_type(ii).le.530))).AND.&
 &   (.NOT.((b_type(ii).ge.544).AND.(b_type(ii).le.557))).AND.&
 &   (.NOT.((b_type(ii).ge.624).AND.(b_type(ii).le.962))).AND.&
 &   (.NOT.((b_type(ii).ge.1074).AND.(b_type(ii).le.1259)))) then
            tstr(1:6) = 'HETATM'
         end if
        end if
!
!       write to file but watch for formatting errors (correct crudely)
        attiwr = atti
        if (atti.gt.99999) attiwr = mod(atti,100000)
        rswr = rs-rsshft+rsshft2
        if (rswr.gt.9999) rswr = mod(rswr,10000)
        if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
          outv(1) = x(ii)
          outv(2) = y(ii)
          outv(3) = z(ii)
          if (bnd_shape.eq.1) then
            do k=1,3  
              if (outv(k).lt.bnd_params(k+3)) then
                outv(k) = outv(k) + ceiling((bnd_params(3+k)-outv(k))/bnd_params(k))*bnd_params(k)
              else if (outv(k).gt.(bnd_params(k+3)+bnd_params(k))) then
                outv(k) = outv(k) - ceiling((outv(k)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
              end if
            end do
          else if (bnd_shape.eq.3) then
            if (outv(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
              outv(3) = outv(3) + ceiling((bnd_params(3)-0.5*bnd_params(6)-outv(3))/bnd_params(6))*bnd_params(6)
            else if (outv(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
              outv(3) = outv(3) - ceiling((outv(3)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
            end if
          end if
          write (ipdb,fpdbstr) tstr,attiwr,ana,extra,resnameprt,chch,rswr,extra2,outv(1),outv(2),outv(3)+shf,occup
        else       
          write (ipdb,fpdbstr) tstr,attiwr,ana,extra,resnameprt,chch,rswr,extra2,x(ii)+shf3(1),y(ii)+shf3(2),z(ii)+shf3(3),occup
        end if
!
!       if we write out nucleotides in the exp. format, we need to manually add
!       the 3'-O, which in CAMPARI notation is on the next residue (i+1)
        if (pdb_nucmode.eq.2) then
          if ((seqpolty(rs).eq.'N').AND.(nins.eq.0)) then
            if (rs.lt.rsmol(molofrs(rs),2)) then
!             identify the biotype for C3*
              if(((b_type(ii).eq.733).OR.(b_type(ii).eq.680).OR.&
 &                (b_type(ii).eq.787).OR.(b_type(ii).eq.843).OR.&
 &                (b_type(ii).eq.901).OR.(b_type(ii).eq.707).OR.&
 &                (b_type(ii).eq.654).OR.(b_type(ii).eq.760).OR.&
 &                (b_type(ii).eq.815).OR.(b_type(ii).eq.872).OR.&
 &                (b_type(ii).eq.633).OR.(b_type(ii).eq.945)).OR.&
 &    ((seqtyp(rs).eq.26).AND.((bio_code(b_type(ii)).eq.'C3*').OR.(bio_code(b_type(ii)).eq."C3'").OR.&
 &                             (ii.eq.(i+at(rs)%nbb+at(rs)%nsc-1))))) then
                if (seqflag(rs).eq.22) then
                  if (ii.ne.nuci(rs,6)) then
                    write(ilog,*) 'Warning. Re-assignment of oxygen at&
 &om to residue i-1 to comply with PDB-standard might be corrupted.'
                  end if
                else if (seqflag(rs).eq.24) then ! DIB, RIB, RIX, or DIX
                  if (ii.ne.nuci(rs,4)) then
                    write(ilog,*) 'Warning. Re-assignment of oxygen at&
 &om to residue i-1 to comply with PDB-standard might be corrupted.'
                  end if
                end if
                atti = atti + 1
                ana(1:1) = ' '
                if ((pdb_convention(1).eq.3).OR.(pdb_convention(1).eq.4)) then
                  ana(2:4) = "O3'"
                else
                  ana(2:4) = 'O3*'
                end if
                tstr(1:6) = 'ATOM  '
                extra = ' '
                resname2 = amino(seqtyp(rs+1))
                if (seqtyp(rs+1).ne.26) then
                  kk = nuci(rs+1,1)
                  if ((b_type(kk).eq.726).OR.(b_type(kk).eq.673).OR.&
   &                (b_type(kk).eq.780).OR.(b_type(kk).eq.836).OR.&
   &                (b_type(kk).eq.894).OR.(b_type(kk).eq.700).OR.&
   &                (b_type(kk).eq.647).OR.(b_type(kk).eq.753).OR.&
   &                (b_type(kk).eq.808).OR.(b_type(kk).eq.865).OR.&
   &                (b_type(kk).eq.626).OR.(b_type(kk).eq.938)) then
                    nins = nins + 1
                    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
                      outv(1) = x(kk)
                      outv(2) = y(kk)
                      outv(3) = z(kk)
                      if (bnd_shape.eq.1) then
                        do k=1,3  
                          if (outv(k).lt.bnd_params(k+3)) then
                            outv(k) = outv(k) + ceiling((bnd_params(3+k)-outv(k))/bnd_params(k))*bnd_params(k)
                          else if (outv(k).gt.(bnd_params(k+3)+bnd_params(k))) then
                            outv(k) = outv(k) - ceiling((outv(k)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
                          end if
                        end do
                      else if (bnd_shape.eq.3) then
                        if (outv(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
                          outv(3) = outv(3) + ceiling((bnd_params(3)-0.5*bnd_params(6)-outv(3))/bnd_params(6))*bnd_params(6)
                        else if (outv(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
                          outv(3) = outv(3) - ceiling((outv(3)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
                        end if
                      end if
                      write (ipdb,fpdbstr) tstr,atti,ana,extra,resnameprt,chch,rs,extra2,outv(1),outv(2),outv(3)+shf,occup
                    else       
                      write (ipdb,fpdbstr) tstr,atti,ana,extra,resnameprt,chch,rs,extra2,x(kk),y(kk),z(kk)+shf,occup
                    end if
                  else
                    write(ilog,*) 'Fatal. Re-assignment of oxygen atom&
   & to residue i-1 to comply with PDB-standard is corrupted.'
                    call fexit()
                  end if
                else
                  do kkk=at(rs+1)%bb(1),at(rs+1)%bb(1)+at(rs+1)%nbb+at(rs+1)%nsc-1
                    if ((bio_code(b_type(kkk)).eq.'O3*').OR.(bio_code(b_type(kkk)).eq."O3'")) then
                      kk = kkk
                      exit
                    end if
                  end do
                  nins = nins + 1
                  if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
                    outv(1) = x(kk)
                    outv(2) = y(kk)
                    outv(3) = z(kk)
                    if (bnd_shape.eq.1) then
                      do k=1,3  
                        if (outv(k).lt.bnd_params(k+3)) then
                          outv(k) = outv(k) + ceiling((bnd_params(3+k)-outv(k))/bnd_params(k))*bnd_params(k)
                        else if (outv(k).gt.(bnd_params(k+3)+bnd_params(k))) then
                          outv(k) = outv(k) - ceiling((outv(k)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
                        end if
                      end do
                    else if (bnd_shape.eq.3) then
                      if (outv(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
                        outv(3) = outv(3) + ceiling((bnd_params(3)-0.5*bnd_params(6)-outv(3))/bnd_params(6))*bnd_params(6)
                      else if (outv(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
                        outv(3) = outv(3) - ceiling((outv(3)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
                      end if
                    end if
                    write (ipdb,fpdbstr) tstr,atti,ana,extra,resnameprt,chch,rs,extra2,outv(1),outv(2),outv(3)+shf,occup
                  else       
                    write (ipdb,fpdbstr) tstr,atti,ana,extra,resnameprt,chch,rs,extra2,x(kk),y(kk),z(kk)+shf,occup
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end do
   end do
  end do
  if ((pdb_writemode.eq.1).OR.(ndump.eq.-1)) then
    write(ipdb,34)
  else
    write(ipdb,33)
  end if
  rsw1 = backupval
!
end
!
!-----------------------------------------------------------------------
!
! depending on biotype this routine gives individual names to chemically
! indistinguishable atom types
!
subroutine name_heavies(ii,rs,ana,ncyc)
!
  use iounit
  use atoms
  use pdb
  use polypep
  use system
  use aminos
  use sequen
!
  implicit none
!
  integer ii,rs,shf,ncyc
  character(4) ana
!
  if (seqtyp(rs).eq.26) return ! refuse to name atoms for unkown residues --> array violations
  shf = 0
  if (ua_model.gt.0) then
    shf = 1
  end if
!
  if ((ii.eq.hni(rs)).AND.(pdb_convention(1).eq.4)) ana=' H  '
!
  if (b_type(ii).eq.23) then
    if (ii.eq.at(rs)%sc(3-shf)) ana=' CG1'
    if (ii.eq.at(rs)%sc(4-shf)) ana=' CG2'
  else if (b_type(ii).eq.35) then
    if (ii.eq.at(rs)%sc(4-shf)) ana=' CD1'
    if (ii.eq.at(rs)%sc(5-shf)) ana=' CD2'
  else if (b_type(ii).eq.49) then
    if ((pdb_convention(1).ge.2).AND.(pdb_convention(1).le.3)) then
      ana = ' CD '
    end if
  else if (b_type(ii).eq.112) then
    if (ii.eq.at(rs)%sc(4-shf)) ana=' CD1'
    if (ii.eq.at(rs)%sc(5-shf)) ana=' CD2'
  else if (b_type(ii).eq.114) then
    if (ii.eq.at(rs)%sc(6-shf)) ana=' CE1'
    if (ii.eq.at(rs)%sc(7-shf)) ana=' CE2'
  else if ((b_type(ii).eq.127).OR.(b_type(ii).eq.1134).OR.(b_type(ii).eq.1149)) then
    if (ii.eq.at(rs)%sc(4-shf)) ana=' CD1'
    if (ii.eq.at(rs)%sc(5-shf)) ana=' CD2'
  else if ((b_type(ii).eq.129).OR.(b_type(ii).eq.1136).OR.(b_type(ii).eq.1151)) then
    if (ii.eq.at(rs)%sc(6-shf)) ana=' CE1'
    if (ii.eq.at(rs)%sc(7-shf)) ana=' CE2'
  else if (b_type(ii).eq.215) then
    if (ii.eq.at(rs)%sc(4-shf)) ana=' OD1'
    if (ii.eq.at(rs)%sc(5-shf)) ana=' OD2'
  else if (b_type(ii).eq.239) then
    if (ii.eq.at(rs)%sc(5-shf)) ana=' OE1'
    if (ii.eq.at(rs)%sc(6-shf)) ana=' OE2'
  else if (b_type(ii).eq.298) then
    if (ii.eq.at(rs)%sc(7-shf)) ana=' NH1'
    if (ii.eq.at(rs)%sc(8-shf)) ana=' NH2'
  else if (b_type(ii).eq.319) then
    if (ii.eq.at(rs)%sc(1)) ana=' CB1'
    if (ii.eq.at(rs)%sc(2)) ana=' CB2'
  else if ((b_type(ii).eq.357).OR.(b_type(ii).eq.363).OR.&
 &         (b_type(ii).eq.368).OR.(b_type(ii).eq.553)) then
    if (pdb_convention(1).eq.1) then
      if (ii.eq.at(rs)%bb(4)) ana='1OXT'
      if (ii.eq.at(rs)%bb(5)) ana='2OXT'
    else if (pdb_convention(1).eq.2) then
      if (ii.eq.at(rs)%bb(4)) ana=' O1 '
      if (ii.eq.at(rs)%bb(5)) ana=' O2 '
    else if (pdb_convention(1).eq.3) then
      if (ii.eq.at(rs)%bb(4)) ana=' OT1'
      if (ii.eq.at(rs)%bb(5)) ana=' OT2'
    else if (pdb_convention(1).eq.4) then
      if (ii.eq.at(rs)%bb(4)) ana=' O  '
      if (ii.eq.at(rs)%bb(5)) ana=' OXT'
    end if
  else if (b_type(ii).eq.372) then
    if (pdb_convention(1).eq.3) then
      ana(2:4) = 'CAY'
    else if (pdb_convention(1).eq.2) then
      ana(2:4) = 'CA '
    end if
  else if ((b_type(ii).eq.374).OR.(b_type(ii).eq.375)) then
    if (pdb_convention(1).eq.3) ana(3:3) = 'Y'
  else if (b_type(ii).eq.379) then
    if (pdb_convention(1).eq.3) ana(3:3) = 'T'
  else if (b_type(ii).eq.381) then
    if (pdb_convention(1).eq.3) then
      ana(2:4) = 'CAT'
    else if (pdb_convention(1).eq.2) then
      ana(2:4) = 'CA '
    end if
  else if (b_type(ii).eq.383) then
    if (pdb_convention(1).eq.3) ana(2:4) = 'NT '
  else if (b_type(ii).eq.444) then
    if (pdb_convention(1).eq.4) ana(1:4) = ' Na+'
  else if (b_type(ii).eq.445) then
    if (pdb_convention(1).eq.4) ana(1:4) = ' Cl-'
  else if (b_type(ii).eq.446) then
    if (pdb_convention(1).eq.2) ana(1:4) = ' CU '
  else if (b_type(ii).eq.447) then
    if (pdb_convention(1).eq.2) ana(1:4) = ' OU '
  else if (b_type(ii).eq.448) then
    if (pdb_convention(1).eq.2) then
      if (ii.eq.at(rs)%bb(3)) ana=' N1U'
      if (ii.eq.at(rs)%bb(4)) ana=' N2U'
    else
      if (ii.eq.at(rs)%bb(3)) ana=' N1 '
      if (ii.eq.at(rs)%bb(4)) ana=' N2 '
    end if
  else if ((b_type(ii).eq.450).OR.(b_type(ii).eq.452).OR.(b_type(ii).eq.1050).OR.(b_type(ii).eq.483).OR.(b_type(ii).eq.486)) then
     if (pdb_convention(1).eq.4) then
       if (ii.eq.at(rs)%bb(1)) ana =' O  '
     end if
  else if (b_type(ii).eq.505) then
    if (pdb_convention(1).eq.2) ana(1:4) = 'CMet'
  else if (b_type(ii).eq.506) then
    if (pdb_convention(1).eq.2) ana(1:4) = 'OMet'
  else if (b_type(ii).eq.512) then
    if (ii.eq.at(rs)%bb(3)) ana=' C21'
    if (ii.eq.at(rs)%bb(4)) ana=' C22' 
  else if (b_type(ii).eq.514) then
    if (ii.eq.at(rs)%bb(5)) ana=' C31'
    if (ii.eq.at(rs)%bb(6)) ana=' C32'
  else if (b_type(ii).eq.538) then
    if (ii.eq.at(rs)%bb(2)) ana='1OXT'
    if (ii.eq.at(rs)%bb(3)) ana='2OXT'
  else if (b_type(ii).eq.542) then
    if (ii.eq.at(rs)%bb(2)) ana=' N1 '
    if (ii.eq.at(rs)%bb(3)) ana=' N2 '
    if (ii.eq.at(rs)%bb(4)) ana=' N3 '
  else if ((b_type(ii).eq.648).OR.(b_type(ii).eq.674).OR.&
 &         (b_type(ii).eq.701).OR.(b_type(ii).eq.727).OR.&
 &         (b_type(ii).eq.754).OR.(b_type(ii).eq.781).OR.&
 &         (b_type(ii).eq.809).OR.(b_type(ii).eq.837).OR.&
 &         (b_type(ii).eq.866).OR.(b_type(ii).eq.895).OR.&
 &         (b_type(ii).eq.627).OR.(b_type(ii).eq.939)) then
    if (ii.eq.at(rs)%bb(4)) ana=' O1P'
    if (ii.eq.at(rs)%bb(5)) ana=' O2P'
  else if ((b_type(ii).eq.625).OR.(b_type(ii).eq.646).OR.&
 &         (b_type(ii).eq.672).OR.(b_type(ii).eq.699).OR.&
 &         (b_type(ii).eq.725).OR.(b_type(ii).eq.752).OR.&
 &         (b_type(ii).eq.779).OR.(b_type(ii).eq.807).OR.&
 &         (b_type(ii).eq.835).OR.(b_type(ii).eq.864).OR.&
 &         (b_type(ii).eq.893).OR.(b_type(ii).eq.937)) then
!   the vmd ribbon display doesn't like O3P/O5P, so we'll sugar em up
    ana=' O5*'
  else if ((b_type(ii).eq.626).OR.(b_type(ii).eq.647).OR.&
 &         (b_type(ii).eq.673).OR.(b_type(ii).eq.700).OR.&
 &         (b_type(ii).eq.726).OR.(b_type(ii).eq.753).OR.&
 &         (b_type(ii).eq.780).OR.(b_type(ii).eq.808).OR.&
 &         (b_type(ii).eq.836).OR.(b_type(ii).eq.865).OR.&
 &         (b_type(ii).eq.894).OR.(b_type(ii).eq.938)) then
!   the vmd ribbon display doesn't like O3P/O5P, so we'll sugar em up
    ana=' O3*'
!   in CAMPARI-convention there's two O3* in 3'-terminal residues, so
!   we'll modify the terminal one
  else if (b_type(ii).eq.774) then
    if (pdb_convention(1).eq.4) ana = ' C7 '
  else if (b_type(ii).eq.802) then
    if (pdb_convention(1).eq.4) ana = ' C7 '
  else if ((b_type(ii).eq.929).OR.(b_type(ii).eq.926).OR.&
 &         (b_type(ii).eq.634).OR.(b_type(ii).eq.946)) then
    if (ncyc.eq.0) then
      ana(1:1) = '2'
    end if
!   in PDB-convention there's two O3* in 5'-P-terminal residues, so
!   we'll modify the terminal one
  else if (b_type(ii).eq.923) then
    if (pdb_convention(1).eq.3) then
      ana=" O5T"
    else if ((pdb_nucmode.eq.2).AND.(pdb_convention(1).eq.1)) then
      ana(1:4) = '2O3*'
    else
      ana(1:4) = ' O3*'
    end if
  else if ((b_type(ii).eq.924).OR.(b_type(ii).eq.955).OR.(b_type(ii).eq.956)) then
    if ((pdb_convention(1).eq.3).OR.(pdb_convention(1).eq.4)) then
      ana=" H5T"
    end if
  else if ((b_type(ii).eq.644).OR.(b_type(ii).eq.954)) then
    if (pdb_convention(1).eq.3) then
      ana=" H1'"
    else
      ana='1HO*'
!      extra = '*'
    end if
  else if (((b_type(ii).ge.931).AND.(b_type(ii).le.935)).OR.&
 &         (b_type(ii).eq.640)) then
    if (pdb_convention(1).eq.3) then
      ana=" H2'"
    else if (pdb_convention(1).eq.2) then
      ana=" H2*"
    else if (pdb_convention(1).eq.4) then
      ana="HO'2"
    else
      ana='2HO*'
!      extra = '*'
    end if
  else if ((b_type(ii).eq.927).OR.(b_type(ii).eq.930).OR.&
 &         (b_type(ii).eq.948).OR.(b_type(ii).eq.636)) then
    if ((pdb_convention(1).eq.3).OR.(pdb_convention(1).eq.4)) then
      ana=' H3T'
    else if (pdb_convention(1).eq.2) then
      ana=' H3*' ! doesn't exist in GROMOS
    else
      ana='3HO*'
!      extra = '*'
    end if
  else if ((b_type(ii).eq.959).OR.(b_type(ii).eq.962)) then
    if ((pdb_convention(1).eq.3).OR.(pdb_convention(1).eq.4)) then
      ana=' H5T'
    else
      ana='5HO*'
!      extra = '*' 
    end if
  else if (b_type(ii).eq.963) then
    if (ii.eq.at(rs)%bb(1)) ana=' CT1'
    if (ii.eq.at(rs)%bb(3)) ana=' CT2' 
  else if (b_type(ii).eq.967) then
    if (ii.eq.at(rs)%bb(1)) ana=' CT1'
    if (ii.eq.at(rs)%sc(1)) ana=' CT2'
  else if (b_type(ii).eq.969) then
    if (ii.eq.at(rs)%bb(2)) ana=' CB1'
    if (ii.eq.at(rs)%bb(3)) ana=' CB2'
  else if (b_type(ii).eq.971) then
    if (ii.eq.at(rs)%bb(1)) ana=' CT1'
    if (ii.eq.at(rs)%bb(3)) ana=' CT2'
    if (ii.eq.at(rs)%bb(4)) ana=' CT3'
  else if (b_type(ii).eq.978) then
    if (ii.eq.at(rs)%bb(3)) ana=' C21'
    if (ii.eq.at(rs)%bb(4)) ana=' C22'
  else if (b_type(ii).eq.980) then
    if (ii.eq.at(rs)%bb(5)) ana=' C31'
    if (ii.eq.at(rs)%bb(6)) ana=' C32'
  else if (b_type(ii).eq.1044) then
    if (ii.eq.at(rs)%bb(1)) ana=' CT1'
    if (ii.eq.at(rs)%bb(3)) ana=' CT2'
  else if (b_type(ii).eq.1047) then
    if (ii.eq.at(rs)%bb(2)) ana=' O1 '
    if (ii.eq.at(rs)%bb(3)) ana=' O2 '
    if (ii.eq.at(rs)%bb(4)) ana=' O3 '
    if (ii.eq.at(rs)%bb(5)) ana=' O4 '
  else if (b_type(ii).eq.1049) then
    if (ii.eq.at(rs)%bb(2)) ana=' O1 '
    if (ii.eq.at(rs)%bb(3)) ana=' O2 '
    if (ii.eq.at(rs)%bb(4)) ana=' O3 '
  else if (b_type(ii).eq.1053) then
    if (ii.eq.at(rs)%bb(1)) ana=' C1 '
    if (ii.eq.at(rs)%bb(2)) ana=' C2 '
    if (ii.eq.at(rs)%bb(3)) ana=' C3 '
    if (ii.eq.at(rs)%bb(4)) ana=' C4 '
    if (ii.eq.at(rs)%bb(5)) ana=' C5 '
    if (ii.eq.at(rs)%bb(6)) ana=' C6 '
  else if (b_type(ii).eq.1055) then
    if (ii.eq.at(rs)%bb(1)) ana=' C11'
    if (ii.eq.at(rs)%bb(2)) ana=' C12'
  else if (b_type(ii).eq.1056) then
    if (ii.eq.at(rs)%bb(3)) ana=' C21'
    if (ii.eq.at(rs)%bb(4)) ana=' C22'
    if (ii.eq.at(rs)%bb(7)) ana=' C23'
    if (ii.eq.at(rs)%bb(8)) ana=' C24'
  else if (b_type(ii).eq.1057) then
    if (ii.eq.at(rs)%bb(5)) ana=' C31'
    if (ii.eq.at(rs)%bb(6)) ana=' C32'
    if (ii.eq.at(rs)%bb(9)) ana=' C33'
    if (ii.eq.at(rs)%bb(10)) ana=' C34'
  else if (b_type(ii).eq.1073) then
    if (ii.eq.at(rs)%bb(1)) ana=' O1 '
    if (ii.eq.at(rs)%bb(2)) ana=' O2 '
  else if (b_type(ii).eq.1156) then
    if (ii.eq.at(rs)%sc(12-shf)) ana=' O1P'
    if (ii.eq.at(rs)%sc(13-shf)) ana=' O2P'
  else if (b_type(ii).eq.1169) then
    if (ii.eq.at(rs)%sc(6-shf)) ana=' O1P'
    if (ii.eq.at(rs)%sc(7-shf)) ana=' O2P'
  else if (b_type(ii).eq.1184) then
    if (ii.eq.at(rs)%sc(7-shf)) ana=' O1P'
    if (ii.eq.at(rs)%sc(8-shf)) ana=' O2P'
  else if (b_type(ii).eq.1241) then
    if (ii.eq.at(rs)%sc(7-shf)) ana=' CH1'
    if (ii.eq.at(rs)%sc(8-shf)) ana=' CH2'
  else if (b_type(ii).eq.1259) then
    if (ii.eq.at(rs)%sc(7-shf)) ana=' CH1'
    if (ii.eq.at(rs)%sc(8-shf)) ana=' CH2'
    if (ii.eq.at(rs)%sc(9-shf)) ana=' CH3'
  end if
!
  if (seqpolty(atmres(ii)).eq.'N') then
    if ((pdb_convention(1).eq.3).OR.(pdb_convention(1).eq.4)) then
      if (ana(4:4).eq.'*') ana(4:4) = "'"
    end if
  end if
!
end
!
!-----------------------------------------------------------------------
!
! depending on biotype this routine gives individual names to chemically
! indistinguishable atom types
!
subroutine name_hydrogens(ii,rs,imol,ana)
!
  use iounit
  use atoms
  use pdb
  use polypep
  use molecule
  use sequen
  use fyoc
  use system
  use aminos
!
  implicit none
!
  integer ii,rs,imol,shf,shf2,shf3,shf5,shf7,shf9,shfx2,shfx1
  character(4) ana
!
  if (seqtyp(rs).eq.26) return ! refuse to name atoms for unkown residues --> array violations
  shf = 0
  shf2 = 0
  shf3 = 0
  shf5 = 0
  shf7 = 0
  shf9 = 0
  shfx2 = 0
  shfx1 = 0
  if (ua_model.gt.0) then
    shf = 1
    shf2 = 2
    shf3 = 3
    shf5 = 5
    shf7 = 7
    shf9 = 9
    if (ua_model.eq.2) then
      shfx2 = 2
      shfx1 = 1
    end if
  end if
!
  if (b_type(ii).eq.6) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(1)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(2)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(1)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(2)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.14) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.24) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(6)) ana='1HG1'
      if (ii.eq.at(rs)%sc(7)) ana='2HG1'
      if (ii.eq.at(rs)%sc(8)) ana='3HG1'
      if (ii.eq.at(rs)%sc(9)) ana='1HG2'
      if (ii.eq.at(rs)%sc(10)) ana='2HG2'
      if (ii.eq.at(rs)%sc(11)) ana='3HG2'
    end if
  else if (b_type(ii).eq.32) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.36) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(9)) ana='1HD1'
      if (ii.eq.at(rs)%sc(10)) ana='2HD1'
      if (ii.eq.at(rs)%sc(11)) ana='3HD1'
      if (ii.eq.at(rs)%sc(12)) ana='1HD2'
      if (ii.eq.at(rs)%sc(13)) ana='2HD2'
      if (ii.eq.at(rs)%sc(14)) ana='3HD2'
    end if
  else if ((b_type(ii).eq.46).OR.(b_type(ii).eq.100)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.48) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.50) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(12)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(13)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(14)) ana(1:1)='3'
      if (pdb_convention(1).eq.3) ana(4:4) = ' '
    end if
  else if ((b_type(ii).eq.58).OR.(b_type(ii).eq.80).OR.(b_type(ii).eq.90).OR.(b_type(ii).eq.1123)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(4)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(5)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(4)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(5)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.60).OR.(b_type(ii).eq.82)) then
    if (pdb_convention(1).eq.3) then
      ana(4:4) = '1'
    end if
  else if (b_type(ii).eq.72) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.98) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(5)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(5)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.102).OR.(b_type(ii).eq.110).OR.(b_type(ii).eq.1232)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.113) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%sc(11-shf3)) ana=' HD1'
      if (ii.eq.at(rs)%sc(12-shf3)) ana=' HD2'
    end if
  else if (b_type(ii).eq.115) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%sc(13-shf3)) ana=' HE1'
      if (ii.eq.at(rs)%sc(14-shf3)) ana=' HE2'
    end if
  else if ((b_type(ii).eq.125).OR.(b_type(ii).eq.1132).OR.(b_type(ii).eq.1216)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(11)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.128).OR.(b_type(ii).eq.1135)) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%sc(12-shf3)) ana=' HD1'
      if (ii.eq.at(rs)%sc(13-shf3)) ana=' HD2'
    end if
  else if ((b_type(ii).eq.130).OR.(b_type(ii).eq.1137)) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%sc(14-shf3)) ana=' HE1'
      if (ii.eq.at(rs)%sc(15-shf3)) ana=' HE2'
    end if
  else if ((b_type(ii).eq.141).OR.(b_type(ii).eq.1218)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(12)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(13)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(12)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(13)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.164).OR.(b_type(ii).eq.181).OR.&
 &               (b_type(ii).eq.197).OR.(b_type(ii).eq.1166).OR.(b_type(ii).eq.1214)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.213).OR.(b_type(ii).eq.223).OR.(b_type(ii).eq.1095)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.227) then
    if (ii.eq.at(rs)%sc(8-shf3)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(9-shf3)) ana(1:1)='2'
  else if ((b_type(ii).eq.235).OR.(b_type(ii).eq.247).OR.(b_type(ii).eq.274).OR.(b_type(ii).eq.1107).OR.(b_type(ii).eq.1081)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.237).OR.(b_type(ii).eq.249).OR.(b_type(ii).eq.276).OR.(b_type(ii).eq.290).OR.(b_type(ii).eq.1109).OR.&
 &         (b_type(ii).eq.1083)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.253) then
    if (ii.eq.at(rs)%sc(11-shf5)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(12-shf5)) ana(1:1)='2'
  else if (b_type(ii).eq.261) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.263) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(8)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(9)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.266) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(12)) ana(1:1)='3'
    end if
  else if ((b_type(ii).eq.278).OR.(b_type(ii).eq.292).OR.(b_type(ii).eq.1111).OR.(b_type(ii).eq.1234)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(12)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(11)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(12)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.280).OR.(b_type(ii).eq.294).OR.(b_type(ii).eq.1113).OR.(b_type(ii).eq.1236)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(13)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(14)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(13)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(14)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.282).OR.(b_type(ii).eq.1115)) then
    if (ii.eq.at(rs)%sc(15-shf9)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(16-shf9)) ana(1:1)='2'
    if (ii.eq.at(rs)%sc(17-shf9)) ana(1:1)='3'
  else if (b_type(ii).eq.299) then
    if (ii.eq.at(rs)%sc(16-shf7)) ana='1HH1'
    if (ii.eq.at(rs)%sc(17-shf7)) ana='2HH1'
    if (ii.eq.at(rs)%sc(18-shf7)) ana='1HH2'
    if (ii.eq.at(rs)%sc(19-shf7)) ana='2HH2'
  else if (b_type(ii).eq.307) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.309) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.311) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.313) then
    if (ii.eq.at(rs)%sc(12-shf7)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(13-shf7)) ana(1:1)='2'
    if (ii.eq.at(rs)%sc(14-shf7)) ana(1:1)='3'
  else if (b_type(ii).eq.320) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(3)) ana='1HB1'
      if (ii.eq.at(rs)%sc(4)) ana='2HB1'
      if (ii.eq.at(rs)%sc(5)) ana='3HB1'
      if (ii.eq.at(rs)%sc(6)) ana='1HB2'
      if (ii.eq.at(rs)%sc(7)) ana='2HB2'
      if (ii.eq.at(rs)%sc(8)) ana='3HB2'
    end if
  else if ((b_type(ii).eq.336).OR.(b_type(ii).eq.342).OR.(b_type(ii).eq.547)) then
    if (rs.eq.rsmol(molofrs(rs),2)) then
      if (ii.eq.at(rs)%bb(6)) ana='1H  '
      if (ii.eq.at(rs)%bb(7)) ana='2H  '
      if (moltermid(molofrs(rs),1).eq.1) then
        if (ii.eq.at(rs)%bb(8)) ana='3H  '
      end if
    else
      if (ii.eq.at(rs)%bb(5)) ana='1H  '
      if (ii.eq.at(rs)%bb(6)) ana='2H  '
      if (moltermid(molofrs(rs),1).eq.1) then
        if (ii.eq.at(rs)%bb(7)) ana='3H  '
      end if
    end if
    if (pdb_convention(1).eq.3) ana(3:4) = 'T'
  else if ((b_type(ii).eq.344).OR.(b_type(ii).eq.364)) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana='1HA '
      if (ii.eq.at(rs)%sc(2)) ana='2HA '
    end if
  else if (b_type(ii).eq.348) then
    if (rs.eq.rsmol(molofrs(rs),2)) then
      if (ii.eq.at(rs)%bb(6)) ana=' H  '
      if (moltermid(molofrs(rs),1).eq.1) then
        if (ii.eq.at(rs)%bb(6)) ana='1H  '
        if (ii.eq.at(rs)%bb(7)) ana='2H  '
      end if
    else
      if (ii.eq.at(rs)%bb(5)) ana=' H  '
      if (moltermid(molofrs(rs),1).eq.1) then
        if (ii.eq.at(rs)%bb(5)) ana='1H  '
        if (ii.eq.at(rs)%bb(6)) ana='2H  '
      end if
    end if
    if (pdb_convention(1).eq.3) ana(3:3) = 'N'
  else if (b_type(ii).eq.352) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.371) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.373) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(2)) ana='1H  '
      if (ii.eq.at(rs)%sc(3)) ana='2H  '
      if (ii.eq.at(rs)%sc(4)) ana='3H  '
      if (pdb_convention(1).eq.3) ana(3:3) = 'Y'
    end if
  else if (b_type(ii).eq.380) then
    if (pdb_convention(1).eq.3) ana(4:4) = 'T'
  else if (b_type(ii).eq.382) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(4)) ana='1H  '
      if (ii.eq.at(rs)%bb(5)) ana='2H  '
      if (ii.eq.at(rs)%bb(6)) ana='3H  '
      if (pdb_convention(1).eq.3) ana(3:4) = 'AT'
    end if
  else if (b_type(ii).eq.384) then
    if (ii.eq.at(rs)%bb(2)) ana='1H  '
    if (ii.eq.at(rs)%bb(3)) ana='2H  '
    if (pdb_convention(1).eq.3) ana(3:4) = 'T '
  else if (b_type(ii).eq.415) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.417) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.425) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.427) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.429) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.437) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.439) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.441) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.443) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(12)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(13)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(14)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.449) then
    if (pdb_convention(1).eq.2) then
      if (ii.eq.at(rs)%bb(5)) ana='H11U'
      if (ii.eq.at(rs)%bb(6)) ana='H12U'
      if (ii.eq.at(rs)%bb(7)) ana='H21U'
      if (ii.eq.at(rs)%bb(8)) ana='H22U'
    else if (pdb_convention(1).eq.4) then
      if (ii.eq.at(rs)%bb(5)) ana=' H11'
      if (ii.eq.at(rs)%bb(6)) ana=' H12'
      if (ii.eq.at(rs)%bb(7)) ana=' H21'
      if (ii.eq.at(rs)%bb(8)) ana=' H22'
    else
      if (ii.eq.at(rs)%bb(5)) ana='1HN1'
      if (ii.eq.at(rs)%bb(6)) ana='2HN1'
      if (ii.eq.at(rs)%bb(7)) ana='1HN2'
      if (ii.eq.at(rs)%bb(8)) ana='2HN2'
    end if
  else if (b_type(ii).eq.459) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.466) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.468) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='3'
    end if
  else if ((b_type(ii).eq.472).OR.(b_type(ii).eq.478)) then
    if (ii.eq.at(rs)%bb(4)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(5)) ana(1:1)='2'
  else if (b_type(ii).eq.474) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.480) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.482) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='3'
    end if
  else if ((b_type(ii).eq.484).OR.(b_type(ii).eq.487).OR.(b_type(ii).eq.1051).OR.(b_type(ii).eq.451).OR.(b_type(ii).eq.453)) then
    if (ii.eq.at(rs)%bb(2)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(3)) ana(1:1)='2'
    if (pdb_convention(1).eq.4) then
      if (ii.eq.at(rs)%bb(2)) ana(2:4) = 'H  '
      if (ii.eq.at(rs)%bb(3)) ana(2:4) = 'H  '
    end if
  else if ((b_type(ii).eq.485).OR.(b_type(ii).eq.1052)) then
    if (pdb_convention(1).eq.4) then
      if (ii.eq.at(rs)%bb(4)) ana(1:4)=" EPW"
    end if
  else if (b_type(ii).eq.488) then
    if (ii.eq.at(rs)%bb(4)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(5)) ana(1:1)='2'
    if (pdb_convention(1).eq.4) then
      if (ii.eq.at(rs)%bb(4)) ana(2:4)='EPW'
      if (ii.eq.at(rs)%bb(5)) ana(2:4)='EPW'
    end if
  else if (b_type(ii).eq.490) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(2)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(3)) ana(1:1)='2'
      if (ii.eq.at(rs)%bb(4)) ana(1:1)='3'
      if (ii.eq.at(rs)%bb(5)) ana(1:1)='4'
    end if
  else if (b_type(ii).eq.494) then
    if (ii.eq.at(rs)%bb(5-shf)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(6-shf)) ana(1:1)='2'
  else if (b_type(ii).eq.500) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(12)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.502) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.504) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='3'
    end if
  else if ((b_type(ii).eq.507).OR.(b_type(ii).eq.992)) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(2)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(3)) ana(1:1)='2'
      if (ii.eq.at(rs)%bb(4)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.508) then
     if (pdb_convention(1).eq.2) then
       ana = 'HMet'
     end if
  else if (b_type(ii).eq.510) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.513) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(9)) ana=' H21'
      if (ii.eq.at(rs)%bb(10)) ana=' H22'
    end if
  else if (b_type(ii).eq.515) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(11)) ana=' H31'
      if (ii.eq.at(rs)%bb(12)) ana=' H32'
    end if
  else if (b_type(ii).eq.526) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.528) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.530) then
    if (ii.eq.at(rs)%sc(9-shf5)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(10-shf5)) ana(1:1)='2'
    if (ii.eq.at(rs)%sc(11-shf5)) ana(1:1)='3'
  else if (b_type(ii).eq.536) then
    if (ii.eq.at(rs)%bb(2)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(3)) ana(1:1)='2'
    if (ii.eq.at(rs)%bb(4)) ana(1:1)='3'
    if (ii.eq.at(rs)%bb(5)) ana(1:1)='4'
  else if ((b_type(ii).eq.540).OR.(b_type(ii).eq.995)) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(4)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.543) then
    if (ii.eq.at(rs)%bb(5)) ana='1HN1'
    if (ii.eq.at(rs)%bb(6)) ana='2HN1'
    if (ii.eq.at(rs)%bb(7)) ana='1HN2'
    if (ii.eq.at(rs)%bb(8)) ana='2HN2'
    if (ii.eq.at(rs)%bb(9)) ana='1HN3'
    if (ii.eq.at(rs)%bb(10)) ana='2HN3'
  else if (b_type(ii).eq.565) then
    if (ii.eq.at(rs)%bb(12-shfx2)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(13-shfx2)) ana(1:1)='2'
  else if ((b_type(ii).eq.592).OR.(b_type(ii).eq.984).OR.&
 &         (b_type(ii).eq.976)) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.603) then
    if (ii.eq.at(rs)%bb(14-shfx2)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(15-shfx2)) ana(1:1)='2'
  else if (b_type(ii).eq.613) then
    if (ii.eq.at(rs)%bb(15-shfx1)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(16-shfx1)) ana(1:1)='2'
  else if ((b_type(ii).eq.650).OR.(b_type(ii).eq.676).OR.&
 &         (b_type(ii).eq.703).OR.(b_type(ii).eq.729).OR.&
 &         (b_type(ii).eq.756).OR.(b_type(ii).eq.783).OR.&
 &         (b_type(ii).eq.811).OR.(b_type(ii).eq.839).OR.&
 &         (b_type(ii).eq.868).OR.(b_type(ii).eq.897).OR.&
 &         (b_type(ii).eq.629).OR.(b_type(ii).eq.941)) then
    if (ua_model.eq.0) then
      if ((seqtyp(rs).ge.76).AND.(seqtyp(rs).le.87)) then ! DIB, RIB, RIX, or DIX
        if (pdb_convention(1).eq.3) then
          if (ii.eq.at(rs)%bb(5)) ana(1:4)=" H5'"
          if (ii.eq.at(rs)%bb(6)) ana(1:4)="H5''"
        else if (pdb_convention(1).eq.4) then
          if (ii.eq.at(rs)%bb(5)) ana(1:4)="H5'1"
          if (ii.eq.at(rs)%bb(6)) ana(1:4)="H5'2"
        else
          if (ii.eq.at(rs)%bb(5)) ana(1:1)='1'
          if (ii.eq.at(rs)%bb(6)) ana(1:1)='2'
        end if
      else if (rs.eq.rsmol(imol,2)) then
        if (pdb_convention(1).eq.3) then
          if (ii.eq.at(rs)%bb(10)) ana(1:4)=" H5'"
          if (ii.eq.at(rs)%bb(11)) ana(1:4)="H5''"
        else if (pdb_convention(1).eq.4) then
          if (ii.eq.at(rs)%bb(10)) ana(1:4)="H5'1"
          if (ii.eq.at(rs)%bb(11)) ana(1:4)="H5'2"
        else
          if (ii.eq.at(rs)%bb(10)) ana(1:1)='1'
          if (ii.eq.at(rs)%bb(11)) ana(1:1)='2'
        end if
      else
        if (pdb_convention(1).eq.3) then
          if (ii.eq.at(rs)%bb(9)) ana(1:4)=" H5'"
          if (ii.eq.at(rs)%bb(10)) ana(1:4)="H5''"
        else if (pdb_convention(1).eq.4) then
          if (ii.eq.at(rs)%bb(9)) ana(1:4)="H5'1"
          if (ii.eq.at(rs)%bb(10)) ana(1:4)="H5'2"
        else
          if (ii.eq.at(rs)%bb(9)) ana(1:1)='1'
          if (ii.eq.at(rs)%bb(10)) ana(1:1)='2' 
        end if
      end if
    end if
  else if ((b_type(ii).eq.657).OR.(b_type(ii).eq.710).OR.&
 &         (b_type(ii).eq.763).OR.(b_type(ii).eq.818).OR.&
 &         (b_type(ii).eq.875).OR.(b_type(ii).eq.950)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.3) then
        if (ii.eq.at(rs)%sc(4)) ana(1:4)=" H2'"
        if (ii.eq.at(rs)%sc(5)) ana(1:4)="H2''"
      else if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(4)) ana(1:4)="H2'1"
        if (ii.eq.at(rs)%sc(5)) ana(1:4)="H2'2"
      else
        if (ii.eq.at(rs)%sc(4)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(5)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.639).OR.(b_type(ii).eq.684).OR.&
 &         (b_type(ii).eq.737).OR.(b_type(ii).eq.791).OR.&
 &         (b_type(ii).eq.847).OR.(b_type(ii).eq.905)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.3) then
        ana(1:4)="H2''"
      else if (pdb_convention(1).eq.4) then
        ana(1:4)="H2'1"
      end if
    end if 
  else if (b_type(ii).eq.666) then
    if (ii.eq.at(rs)%sc(17-shf3-shfx2)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(18-shf3-shfx2)) ana(1:1)='2'
  else if (b_type(ii).eq.693) then
    if (ii.eq.at(rs)%sc(18-shf2-shfx2)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(19-shf2-shfx2)) ana(1:1)='2'
  else if (b_type(ii).eq.775) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(18)) ana=' H71'
        if (ii.eq.at(rs)%sc(19)) ana=' H72'
        if (ii.eq.at(rs)%sc(20)) ana=' H73'
      else if (pdb_convention(1).eq.3) then
        if (ii.eq.at(rs)%sc(18)) ana=' H51'
        if (ii.eq.at(rs)%sc(19)) ana=' H52'
        if (ii.eq.at(rs)%sc(20)) ana=' H53'
      else
        if (ii.eq.at(rs)%sc(18)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(19)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(20)) ana(1:1)='3'
      end if
    end if
  else if (b_type(ii).eq.803) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(19)) ana=' H71'
        if (ii.eq.at(rs)%sc(20)) ana=' H72'
        if (ii.eq.at(rs)%sc(21)) ana=' H73'
      else if (pdb_convention(1).eq.3) then
        if (ii.eq.at(rs)%sc(19)) ana=' H51'
        if (ii.eq.at(rs)%sc(20)) ana=' H52'
        if (ii.eq.at(rs)%sc(21)) ana=' H53'
      else
        if (ii.eq.at(rs)%sc(19)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(20)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(21)) ana(1:1)='3'
      end if
    end if
  else if (b_type(ii).eq.829) then
    if (ii.eq.at(rs)%sc(19-shf3-shfx2)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(20-shf3-shfx2)) ana(1:1)='2'
  else if (b_type(ii).eq.858) then
    if (ii.eq.at(rs)%sc(20-shf2-shfx2)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(21-shf2-shfx2)) ana(1:1)='2'
  else if (b_type(ii).eq.882) then
    if (ii.eq.at(rs)%sc(20-shf3-shfx1)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(21-shf3-shfx1)) ana(1:1)='2'
  else if (b_type(ii).eq.912) then
    if (ii.eq.at(rs)%sc(21-shf2-shfx1)) ana(1:1)='1'
    if (ii.eq.at(rs)%sc(22-shf2-shfx1)) ana(1:1)='2'
  else if (b_type(ii).eq.964) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana='1HT1'
      if (ii.eq.at(rs)%sc(2)) ana='2HT1'
      if (ii.eq.at(rs)%sc(3)) ana='3HT1'
      if (ii.eq.at(rs)%sc(4)) ana='1HT2'
      if (ii.eq.at(rs)%sc(5)) ana='2HT2'
      if (ii.eq.at(rs)%sc(6)) ana='3HT2'
    end if
  else if (b_type(ii).eq.966) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(4)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(5)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.968) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(2)) ana='1HT1'
      if (ii.eq.at(rs)%sc(3)) ana='2HT1'
      if (ii.eq.at(rs)%sc(4)) ana='3HT1'
      if (ii.eq.at(rs)%sc(7)) ana='1HT2'
      if (ii.eq.at(rs)%sc(8)) ana='2HT2'
      if (ii.eq.at(rs)%sc(9)) ana='3HT2'
    end if
  else if (b_type(ii).eq.970) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(4)) ana='1HB1'
      if (ii.eq.at(rs)%bb(5)) ana='2HB1'
      if (ii.eq.at(rs)%sc(5)) ana='1HB2'
      if (ii.eq.at(rs)%sc(6)) ana='2HB2'
    end if
  else if (b_type(ii).eq.972) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana='1HT1'
      if (ii.eq.at(rs)%sc(2)) ana='2HT1'
      if (ii.eq.at(rs)%sc(3)) ana='3HT1'
      if (ii.eq.at(rs)%sc(4)) ana='1HT2'
      if (ii.eq.at(rs)%sc(5)) ana='2HT2'
      if (ii.eq.at(rs)%sc(6)) ana='3HT2'
      if (ii.eq.at(rs)%sc(7)) ana='1HT3'
      if (ii.eq.at(rs)%sc(8)) ana='2HT3'
      if (ii.eq.at(rs)%sc(9)) ana='3HT3'
    end if
  else if (b_type(ii).eq.979) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(8)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(9)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.981) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(10)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(11)) ana(1:1)='2'
    end if
  else if ((b_type(ii).eq.985).OR.(b_type(ii).eq.1002).OR.&
 &         (b_type(ii).eq.1012).OR.(b_type(ii).eq.1022)) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.987) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(4)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(5)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.992) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%bb(3)) ana(1:1)='1'
      if (ii.eq.at(rs)%bb(4)) ana(1:1)='2'
      if (ii.eq.at(rs)%bb(5)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.998) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(5)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(6)) ana(1:1)='2'
    end if
  else if (b_type(ii).eq.1000) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(7)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(8)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(9)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.1039) then
    if (ii.eq.at(rs)%bb(3)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(4)) ana(1:1)='2'
    if (ii.eq.at(rs)%bb(5)) ana(1:1)='3'
  else if (b_type(ii).eq.1041) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(2)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(3)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.1043) then
    if (ii.eq.at(rs)%bb(4)) ana(1:1)='1'
    if (ii.eq.at(rs)%bb(5)) ana(1:1)='2'
  else if (b_type(ii).eq.1045) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(1)) ana(1:4)='1HT1'
      if (ii.eq.at(rs)%sc(2)) ana(1:4)='2HT1'
      if (ii.eq.at(rs)%sc(3)) ana(1:4)='3HT1'
      if (ii.eq.at(rs)%sc(4)) ana(1:4)='1HT2'
      if (ii.eq.at(rs)%sc(5)) ana(1:4)='2HT2'
      if (ii.eq.at(rs)%sc(6)) ana(1:4)='3HT2'
    end if
  else if (b_type(ii).eq.1054) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(7)) ana=' H1 '
      if (ii.eq.at(rs)%bb(8)) ana=' H2 '
      if (ii.eq.at(rs)%bb(9)) ana=' H3 '
      if (ii.eq.at(rs)%bb(10)) ana=' H4 '
      if (ii.eq.at(rs)%bb(11)) ana=' H5 '
      if (ii.eq.at(rs)%bb(12)) ana=' H6 '
    end if
  else if (b_type(ii).eq.1058) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(11)) ana=' H21'
      if (ii.eq.at(rs)%bb(12)) ana=' H22'
      if (ii.eq.at(rs)%bb(15)) ana=' H23'
      if (ii.eq.at(rs)%bb(16)) ana=' H24'
    end if
  else if (b_type(ii).eq.1059) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%bb(13)) ana=' H31'
      if (ii.eq.at(rs)%bb(14)) ana=' H32'
      if (ii.eq.at(rs)%bb(17)) ana=' H33'
      if (ii.eq.at(rs)%bb(18)) ana=' H34'
    end if
  else if (b_type(ii).eq.1182) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(12)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.1147) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(14)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(15)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(14)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(15)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.1150) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%sc(16-shf3)) ana=' HD1'
      if (ii.eq.at(rs)%sc(17-shf3)) ana=' HD2'
    end if
  else if (b_type(ii).eq.1152) then
    if (ua_model.lt.2) then
      if (ii.eq.at(rs)%sc(18-shf3)) ana=' HE1'
      if (ii.eq.at(rs)%sc(19-shf3)) ana=' HE2'
    end if
  else if ((b_type(ii).eq.1194).OR.(b_type(ii).eq.1250)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(11)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(10)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(11)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.1196).OR.(b_type(ii).eq.1252)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(12)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(13)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(12)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(13)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.1198).OR.(b_type(ii).eq.1220).OR.(b_type(ii).eq.1254)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(14)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(15)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(14)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(15)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.1200).OR.(b_type(ii).eq.1256)) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(16)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(17)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(16)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(17)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.1206) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(19)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(20)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(21)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.1222) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(16)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(17)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(16)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(17)) ana(1:1)='2'
      end if
    end if
  else if (b_type(ii).eq.1224) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(18)) ana(1:1)='1'
      if (ii.eq.at(rs)%sc(19)) ana(1:1)='2'
      if (ii.eq.at(rs)%sc(20)) ana(1:1)='3'
    end if
  else if (b_type(ii).eq.1238) then
    if (ua_model.eq.0) then
      if (pdb_convention(1).eq.4) then
        if (ii.eq.at(rs)%sc(15)) ana(1:1)='2'
        if (ii.eq.at(rs)%sc(16)) ana(1:1)='3'
      else
        if (ii.eq.at(rs)%sc(15)) ana(1:1)='1'
        if (ii.eq.at(rs)%sc(16)) ana(1:1)='2'
      end if
    end if
  else if ((b_type(ii).eq.1242).OR.(b_type(ii).eq.1259)) then
    if (ua_model.eq.0) then
      if (ii.eq.at(rs)%sc(18)) ana(1:4)='1HH1'
      if (ii.eq.at(rs)%sc(19)) ana(1:4)='2HH1'
      if (ii.eq.at(rs)%sc(20)) ana(1:4)='3HH1'
      if (ii.eq.at(rs)%sc(21)) ana(1:4)='1HH2'
      if (ii.eq.at(rs)%sc(22)) ana(1:4)='2HH2'
      if (ii.eq.at(rs)%sc(23)) ana(1:4)='3HH2'
      if (ii.eq.at(rs)%sc(24)) ana(1:4)='1HH3'
      if (ii.eq.at(rs)%sc(25)) ana(1:4)='2HH3'
      if (ii.eq.at(rs)%sc(26)) ana(1:4)='3HH3'
    end if
  end if
  if (seqpolty(atmres(ii)).eq.'N') then
    if ((pdb_convention(1).eq.3).OR.(pdb_convention(1).eq.4)) then
      if (ana(4:4).eq.'*') ana(4:4) = "'"
    end if
  end if
!
end
!
!------------------------------------------------------------------------
!
#ifdef LINK_XDR
!
subroutine FMSMC_prtxtc(xdr,ndump)
!
  use iounit
  use system
  use atoms
  use mcsums
  use pdb
  use molecule
  use sequen
  use grandensembles
!
  implicit none
!
  integer xdr,ndump,magic,ret,natsel,istp,i,j,k,kk,mt,imol,lmol
! note that the XDR-code only works with single precision variables (float in C)
  real(KIND=4) box(3,3),xtcprec,shf,shf3(3)
  RTYPE tvec(3)
  real(KIND=4), ALLOCATABLE:: coords(:)
  logical ismember
!
! magic number for compression routine
  magic = 1995
! box matrix
  do i=1,3
    do j=1,3
      box(i,j) = 0.0
    end do
  end do
  shf = 0.0
  if (((pdb_rmol.gt.0).OR.(pdb_force_box.gt.1)).AND.(bnd_type.eq.1).AND.(bnd_shape.ne.1).AND.(bnd_shape.ne.3)) then
    write(ilog,*) 'Fatal. Encountered unsupported box shape for image adjustment &
 &in FMCSC_prtxtc(...). This is an omission bug (please report).'
    call fexit()
  end if
  if (bnd_shape.eq.1) then
    do i=1,3
      box(i,i) = bnd_params(i)
    end do
  else if (bnd_shape.eq.2) then
    do i=1,3
      box(i,i) = 2.0*bnd_params(4)
    end do
  else if (bnd_shape.eq.3) then
    do i=1,2
      box(i,i) = 2.0*bnd_params(4)
    end do
    box(3,3) = bnd_params(6)
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in FMS&
 &MC_prtxtc(...) (code # ',bnd_shape,').'
    call fexit()
  end if
!
  if (ndump.lt.1) return ! should not happen
  allocate(coords(3*pdbeffn))
! coordinates are converted to nm to be consistent with GROMACS
! otherwise VMD plugin can't read it properly
  if (use_trajidx.EQV..true.) then
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      do i=1,pdbeffn
        mt = moltypid(molofrs(atmres(pdboutlst(i))))
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtxtc(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
        end if
        coords(3*i-2) = x(pdboutlst(i))
        coords(3*i-1) = y(pdboutlst(i))
        coords(3*i) = z(pdboutlst(i))
        if (bnd_shape.eq.1) then
          do k=1,3
            kk = k + 3*i - 3
            if (coords(kk).lt.bnd_params(k+3)) then
              coords(kk) = coords(kk) + ceiling((bnd_params(3+k)-coords(kk))/bnd_params(k))*bnd_params(k)
            else if (coords(kk).gt.(bnd_params(k+3)+bnd_params(k))) then
              coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
            end if
          end do
        else if (bnd_shape.eq.3) then
          kk = 3*i
          if (coords(kk).lt.(bnd_params(3)-0.5*bnd_params(6))) then
            coords(kk) = coords(kk) + ceiling((bnd_params(3)-0.5*bnd_params(6)-coords(kk))/bnd_params(6))*bnd_params(6)
          else if (coords(kk).gt.(bnd_params(3)+0.5*bnd_params(6))) then
            coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
          end if
        end if
        coords(3*i) = coords(3*i) + shf
        coords(3*i-2:3*i) = coords(3*i-2:3*i)/10.0
      end do
    else
      lmol = 0
      do i=1,pdbeffn
        imol = molofrs(atmres(pdboutlst(i)))
        mt = moltypid(imol)
        if (pdb_rmol.gt.0) then
          if (lmol.ne.imol) then
            call shift_bound3(pdb_rmol,imol,tvec)
            shf3(:) = tvec(:) ! conversion
          else
            shf3(:) = tvec(:)
          end if
          lmol = imol
        else
          shf3(:) = 0.0
        end if
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtxtc(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
          shf3(3) = shf3(3) + shf
        end if
        coords(3*i-2) = (x(pdboutlst(i))+shf3(1))/10.0
        coords(3*i-1) = (y(pdboutlst(i))+shf3(2))/10.0
        coords(3*i)   = (z(pdboutlst(i))+shf3(3))/10.0
      end do
    end if
  else
    j = 0
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      do imol=1,nmol
        if ((just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
        mt = moltypid(imol)
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtxtc(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
        end if
        do i=atmol(imol,1),atmol(imol,2)
          j = j + 1
          coords(3*j-2) = x(i)
          coords(3*j-1) = y(i)
          coords(3*j) = z(i)
          if (bnd_shape.eq.1) then
            do k=1,3
              kk = k + 3*j - 3
              if (coords(kk).lt.bnd_params(k+3)) then
                coords(kk) = coords(kk) + ceiling((bnd_params(3+k)-coords(kk))/bnd_params(k))*bnd_params(k)
              else if (coords(kk).gt.(bnd_params(k+3)+bnd_params(k))) then
                coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
              end if
            end do
          else if (bnd_shape.eq.3) then
            kk = 3*j
            if (coords(kk).lt.(bnd_params(3)-0.5*bnd_params(6))) then
              coords(kk) = coords(kk) + ceiling((bnd_params(3)-0.5*bnd_params(6)-coords(kk))/bnd_params(6))*bnd_params(6)
            else if (coords(kk).gt.(bnd_params(3)+0.5*bnd_params(6))) then
              coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
            end if
          end if
          coords(3*j) = coords(3*j) + shf
          coords(3*j-2:3*j) = coords(3*j-2:3*j)/10.0
        end do
      end do
    else
      do imol=1,nmol
        if ((just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
        mt = moltypid(imol)
        if (pdb_rmol.gt.0) then
          call shift_bound3(pdb_rmol,imol,tvec)
          shf3(:) = tvec(:) ! conversion
        else
          shf3(:) = 0.0
        end if
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtxtc(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
          shf3(3) = shf3(3) + shf
        end if
        do i=atmol(imol,1),atmol(imol,2)
          j = j + 1
          mt = moltypid(imol)
          coords(3*j-2) = (x(i)+shf3(1))/10.0
          coords(3*j-1) = (y(i)+shf3(2))/10.0
          coords(3*j)   = (z(i)+shf3(3))/10.0
        end do
      end do
    end if
  end if
! if we have nucleotides we might have to swap the O3*s place
!  if (pdb_nucmode.eq.2) then
!    kk = 0
!    lasti = 0
!    do rs=1,nseq
!      do ii=at(rs)%bb(1),at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1
!        kk = kk + 1
!        if (seqpolty(rs).eq.'N') then
!          if (rs.gt.rsmol(molofrs(rs),1)) then
!            if ((b_type(ii).eq.726).OR.(b_type(ii).eq.673).OR.
! &              (b_type(ii).eq.780).OR.(b_type(ii).eq.836).OR.
! &              (b_type(ii).eq.894).OR.(b_type(ii).eq.700).OR.
! &              (b_type(ii).eq.647).OR.(b_type(ii).eq.753).OR.
! &              (b_type(ii).eq.808).OR.(b_type(ii).eq.865).OR.
! &              (b_type(ii).eq.626).OR.(b_type(ii).eq.938)) then
!              lasti = ii
!              bufv(:) = coords(3*kk-2:3*kk)
!            end if
!          end if
!        end if
!        if ((lasti.ne.0).AND.
! &          (ii.lt.(at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1))) then
!          coords(3*kk-2:3*kk) = coords(3*kk+1:3*kk+3)
!        end if
!      end do
!      if (lasti.ne.0) then
!        coords(3*kk-2:3*kk) = bufv(:)
!        lasti = 0
!      end if
!
!          if (rs.lt.rsmol(molofrs(rs),2)) then
!            if ((b_type(ii).eq.733).OR.(b_type(ii).eq.680).OR.
! &              (b_type(ii).eq.787).OR.(b_type(ii).eq.843).OR.
! &              (b_type(ii).eq.901).OR.(b_type(ii).eq.707).OR.
! &              (b_type(ii).eq.654).OR.(b_type(ii).eq.760).OR.
! &              (b_type(ii).eq.815).OR.(b_type(ii).eq.872).OR.
! &              (b_type(ii).eq.633).OR.(b_type(ii).eq.945)) then
!                atti = atti + 1
!                ana(1:1) = ' '
!                ana(2:4) = 'O3*'
!                tstr(1:6) = 'ATOM  '
!                extra = ' '
!                resname2 = amino(seqtyp(rs+1))
!                kk = nuci(rs+1,1)
!                if ((b_type(kk).eq.726).OR.(b_type(kk).eq.673).OR.
! &                (b_type(kk).eq.780).OR.(b_type(kk).eq.836).OR.
! &                (b_type(kk).eq.894).OR.(b_type(kk).eq.700).OR.
! &                (b_type(kk).eq.647).OR.(b_type(kk).eq.753).OR.
! &                (b_type(kk).eq.808).OR.(b_type(kk).eq.865).OR.
! &                (b_type(kk).eq.626).OR.(b_type(kk).eq.938)) then
!                  write (ipdb,30)  tstr,atti,ana,extra,resname
! &                   ,chch,rs,x(kk),y(kk),z(kk)+shf,occup
!                else
!                  write(ilog,*) 'Fatal. Re-assignment of oxygen atom
! & to residue i-1 to comply with PDB-standard is corrupted.'
!                  call fexit()
!                end if
!              end if
!            end if
!          end if
!
!
!
!
!  end if
!
  natsel = pdbeffn
!
  xtcprec = xtc_prec
  istp = nstep
! the header line holds 4 4-bit elements (magic number, number of atoms, number of snapshot, and step number)
  call xdrfint(xdr,magic,ret)
  call xdrfint(xdr,natsel,ret)
  call xdrfint(xdr,ndump,ret)
  call xdrfint(xdr,istp,ret)
! write box matrix
  box(:,:) = box(:,:)/10.0
  do i=1,3
     do j=1,3
       call xdrffloat(xdr,box(i,j),ret)
     end do
  end do
! write coordinates
  call xdrf3dfcoord(xdr,coords,natsel,xtcprec,ret)
! de-allocate
  deallocate(coords)
!
end
!
#endif
!
!------------------------------------------------------------------------
!
!
subroutine FMSMC_prtdcd(dcdf,ndump,dohead)
!
  use iounit
  use system
  use atoms
  use mcsums
  use pdb
  use molecule
  use sequen
  use grandensembles
!
  implicit none
!
  integer ndump,i,mt,imol,lmol,dcdf,j,k,kk
  integer dt(8)
  integer dcdhs(20)
! note that the DCD standard implies single precision variables
  real(KIND=8) box(6)
  real(KIND=4), ALLOCATABLE:: coords(:)
  real(KIND=4) shf,tstp,shf3(3)
  RTYPE tvec(3)
  logical ismember,dohead
  character(4) shortti
  character(80) longti
  character(80) timest
  character(12) rc(3)
!
! box matrix
  shf = 0.0
  if (((pdb_force_box.gt.1).OR.(pdb_rmol.gt.0)).AND.(bnd_type.eq.1).AND.(bnd_shape.ne.1).AND.(bnd_shape.ne.3)) then
    write(ilog,*) 'Fatal. Encountered unsupported box shape for image adjustment &
 &in FMCSC_prtdcd(...). This is an omission bug (please report).'
    call fexit()
  end if
  if (bnd_shape.eq.1) then
    do i=1,3
      box(i) = bnd_params(i)
    end do
    do i=4,6
      box(i) = 90.0  ! WARNING: fix if needed
    end do
  else if (bnd_shape.eq.2) then
    do i=1,3
      box(i) = 2.0*bnd_params(4) ! WARNING: crude
    end do
    do i=4,6
      box(i) = 90.0  ! WARNING: crude
    end do
  else if (bnd_shape.eq.3) then
    do i=1,2
      box(i) = 2.0*bnd_params(4) ! WARNING: crude
    end do
    box(3) = bnd_params(6)
    do i=4,6
      box(i) = 90.0  ! WARNING: crude
    end do
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in FMS&
 &MC_prtdcd(...) (code # ',bnd_shape,').'
    call fexit()
  end if
!
  if (ndump.lt.1) return ! should not happen
!
  if (dohead.EQV..true.) then
!   write header
    call Date_and_Time(rc(1),rc(2),rc(3),dt)
    timest = 'DCD-file written at '//rc(2)(1:2)//':'//rc(2)(3:4)//&
 &' on '//rc(1)(5:6)//'/'//rc(1)(7:8)//'/'//rc(1)(1:4)//'.'
    shortti(1:4) = 'CORD'
    longti = 'CAMPARI Version XYZ run of job "'//basename(1:bleng)//'": '
    k = 0
    do i=nequil+1,nsim
      if (mod(i,xyzout).eq.0) k=k+1
    end do
    dcdhs(1:2) = 0
    do i=nequil+1,nsim
      if (mod(i,xyzout).eq.0) then
        dcdhs(1) = dcdhs(1) + 1
        if (dcdhs(2).le.0) dcdhs(2) = i
      end if
    end do
    dcdhs(3) = xyzout ! skip interval
    dcdhs(4) = nsim ! number of steps
    dcdhs(5:8) = 0 ! unused
    dcdhs(9) = 0 ! this is for "fixed atoms" support (typically: n_atoms - n_ref)
    tstp = 1.0e-12*dyn_dt/4.888821E-14 ! down-conversion from double and conversion to CHARMM ATNA units
    dcdhs(11) = 1 ! whether we have box information with each snap
    dcdhs(12:19) = 0 ! unused
    dcdhs(20) = 24 ! Version number -> CHARMM 24-like format
    write(dcdf) shortti,dcdhs(1:9),tstp,dcdhs(11:20)
    write(dcdf) 2,longti,timest
    write(dcdf) pdbeffn
  end if
!
  allocate(coords(3*pdbeffn))
!
  if (use_trajidx.EQV..true.) then
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      do i=1,pdbeffn
        mt = moltypid(molofrs(atmres(pdboutlst(i))))
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtdcd(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
        end if
        coords(i) = x(pdboutlst(i))
        coords(pdbeffn+i) = y(pdboutlst(i))
        coords(2*pdbeffn+i) = z(pdboutlst(i))
        if (bnd_shape.eq.1) then
          do k=1,3
            kk = (k-1)*pdbeffn + i
            if (coords(kk).lt.bnd_params(k+3)) then
              coords(kk) = coords(kk) + ceiling((bnd_params(3+k)-coords(kk))/bnd_params(k))*bnd_params(k)
            else if (coords(kk).gt.(bnd_params(k+3)+bnd_params(k))) then
              coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
            end if
          end do
        else if (bnd_shape.eq.3) then
          kk = 2*pdbeffn + i
          if (coords(kk).lt.(bnd_params(3)-0.5*bnd_params(6))) then
            coords(kk) = coords(kk) + ceiling((bnd_params(3)-0.5*bnd_params(6)-coords(kk))/bnd_params(6))*bnd_params(6)
          else if (coords(kk).gt.(bnd_params(3)+0.5*bnd_params(6))) then
            coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
          end if
        end if
        coords(2*pdbeffn+i) = coords(2*pdbeffn+i) + shf
      end do
    else
      lmol = 0
      do i=1,pdbeffn
        imol = molofrs(atmres(pdboutlst(i)))
        mt = moltypid(imol)
        if (pdb_rmol.gt.0) then
          if (lmol.ne.imol) then
            call shift_bound3(pdb_rmol,imol,tvec)
            shf3(:) = tvec(:) ! conversion
          else
            shf3(:) = tvec(:)
          end if
          lmol = imol
        else
          shf3(:) = 0.0
        end if
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtdcd(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
          shf3(3) = shf3(3) + shf
        end if
        coords(i) = x(pdboutlst(i))+shf3(1)
        coords(pdbeffn+i) = y(pdboutlst(i))+shf3(2)
        coords(2*pdbeffn+i) = z(pdboutlst(i))+shf3(3)
      end do
    end if
  else
    j = 0
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      do imol=1,nmol
        if ((just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
        mt = moltypid(imol)
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtdcd(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
        end if
        do i=atmol(imol,1),atmol(imol,2)
          j = j + 1
          coords(j) = x(i)
          coords(pdbeffn+j) = y(i)
          coords(2*pdbeffn+j) = z(i)
          if (bnd_shape.eq.1) then
            do k=1,3
              kk = (k-1)*pdbeffn + j
              if (coords(kk).lt.bnd_params(k+3)) then
                coords(kk) = coords(kk) + ceiling((bnd_params(3+k)-coords(kk))/bnd_params(k))*bnd_params(k)
              else if (coords(kk).gt.(bnd_params(k+3)+bnd_params(k))) then
                coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
              end if
            end do
          else if (bnd_shape.eq.3) then
            kk = 2*pdbeffn + j
            if (coords(kk).lt.(bnd_params(3)-0.5*bnd_params(6))) then
              coords(kk) = coords(kk) + ceiling((bnd_params(3)-0.5*bnd_params(6)-coords(kk))/bnd_params(6))*bnd_params(6)
            else if (coords(kk).gt.(bnd_params(3)+0.5*bnd_params(6))) then
              coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
            end if
          end if
          coords(2*pdbeffn+j) = coords(2*pdbeffn+j) + shf
        end do
      end do
    else
      do imol=1,nmol
        if ((just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
        mt = moltypid(imol)
        if (pdb_rmol.gt.0) then
          call shift_bound3(pdb_rmol,imol,tvec)
          shf3(:) = tvec(:) ! conversion
        else
          shf3(:) = 0.0
        end if
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtdcd(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
          shf3(3) = shf3(3) + shf
        end if
        do i=atmol(imol,1),atmol(imol,2)
          j = j + 1
          coords(j) = x(i)+shf3(1)
          coords(pdbeffn+j) = y(i)+shf3(2)
          coords(2*pdbeffn+j) = (z(i)+shf3(3))
        end do
      end do
    end if
  end if
!
! write box vector
  write(dcdf) box(1:6)
! write coordinates
  write(dcdf) coords(1:pdbeffn)
  write(dcdf) coords(pdbeffn+1:2*pdbeffn)
  write(dcdf) coords(2*pdbeffn+1:3*pdbeffn)
! de-allocate
  deallocate(coords)
!
end
!
!------------------------------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
subroutine FMSMC_prtnetcdf(ncid,ndump,dohead)
!
  use iounit
  use system
  use atoms
  use mcsums
  use pdb
  use molecule
  use sequen
  use grandensembles
  use netcdf
!
  implicit none
!
  integer ndump,i,j,mt,imol,lmol,ncid,intin,dims(3),k,kk,framei
!  integer nf90_redef,nf90_put_att,nf90_def_dim
! note that the XDR-code only works with single precision variables (float in C)
  real(KIND=4) shf,timepos,shf3(3)
  real(KIND=4), ALLOCATABLE:: coords(:)
  RTYPE box(3,3),tvec(3)
  logical ismember,dohead
  character(MAXSTRLEN) attstring
  character inp(3)
!
  if (dohead.EQV..true.) then
!   enable definition
!   set AMBER-required global attributes
    do i=1,MAXSTRLEN
      attstring(i:i) = " "
    end do
    attstring(1:5) = "AMBER"
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "Conventions", attstring(1:5)) )
    attstring(1:5) = "     "
    attstring(1:3) = "1.0"
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "ConventionVersion", attstring(1:3)) )
    attstring(1:3) = "   "
    attstring(1:7) = "CAMPARI"
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "program", attstring(1:7)) )
    attstring(1:7) = "       "
    attstring(1:3) = "XXX"
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "programVersion", attstring(1:3)) )
    attstring(1:3) = "   "
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "title", basename(1:bleng)) )
!   define AMBER-required dimensions
    attstring(1:5) = "frame"
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:5), NF90_UNLIMITED, netcdf_ids(1)) )
    attstring(1:5) = "     " 
    attstring(1:7) = "spatial"
    intin = 3
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:7), intin, netcdf_ids(2)) )
    attstring(1:7) = "       "
    attstring(1:4) = "atom"
    intin = pdbeffn
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:4), intin, netcdf_ids(3)) )
    attstring(1:4) = "    "
    attstring(1:12) = "cell_spatial"
    intin = 3
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:12), intin, netcdf_ids(4)) )
    attstring(1:12) = "            "
    attstring(1:12) = "cell_angular"
    intin = 3
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:12), intin, netcdf_ids(5)) )
    attstring(1:12) = "            "
    attstring(1:5) = "label"
    intin = 7
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:5), intin, netcdf_ids(6)) )
    attstring(1:5) = "     "
!   define (not set) AMBER-style label variables
    attstring(1:7) = "spatial"
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:7), NF90_CHAR, netcdf_ids(2), netcdf_ids(16)) )
    attstring(1:7) = "       "
    attstring(1:12) = "cell_spatial"
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:12), NF90_CHAR, netcdf_ids(4), netcdf_ids(17)) )
    attstring(1:12) = "            "
    attstring(1:12) = "cell_angular"
    dims(1:2) = (/ netcdf_ids(6), netcdf_ids(5) /)
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:12), NF90_CHAR, dims(1:2), netcdf_ids(18)) )
    attstring(1:12) = "            "
!   define (not set) AMBER-style physical variables
    attstring(1:4) = "time"
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:4), NF90_FLOAT, netcdf_ids(1), netcdf_ids(11)) )
    attstring(1:5) = "units"
    attstring(6:15) = "picosecond"
    call check_fileionetcdf( nf90_put_att(ncid, netcdf_ids(11), attstring(1:5), attstring(6:15)) )
    attstring(1:15) = "               "
    attstring(1:11) = "coordinates"
    dims(1:3) = (/ netcdf_ids(2), netcdf_ids(3), netcdf_ids(1) /)
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:11), NF90_FLOAT, dims(1:3), netcdf_ids(12)) )
    attstring(1:5) = "units"
    attstring(6:13) = "angstrom"
    call check_fileionetcdf( nf90_put_att(ncid, netcdf_ids(12), attstring(1:5), attstring(6:13)) )
    attstring(1:13) = "             "
    attstring(1:12) = "cell_lengths"
    dims(1:2) = (/ netcdf_ids(4), netcdf_ids(1) /)
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:12), NF90_DOUBLE, dims(1:2), netcdf_ids(13)) )
    attstring(1:5) = "units"
    attstring(6:13) = "angstrom"
    call check_fileionetcdf( nf90_put_att(ncid, netcdf_ids(13), attstring(1:5), attstring(6:13)) )
    attstring(1:13) = "             "
    attstring(1:11) = "cell_angles"
    dims(1:2) = (/ netcdf_ids(5), netcdf_ids(1) /)
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:11), NF90_DOUBLE, dims(1:2), netcdf_ids(14)) )
    attstring(1:5) = "units"
    attstring(6:11) = "degree"
    call check_fileionetcdf( nf90_put_att(ncid, netcdf_ids(14), attstring(1:5), attstring(6:11)) )
    attstring(1:11) = "           "
!   quit define mode
    call check_fileionetcdf( nf90_enddef(ncid) )
!   now write data to label variables
    attstring(1:3) = "xyz"
    inp(1) = 'a'
    call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(16), attstring(1:3)) )
    attstring(1:3) = "abc"
    call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(17), attstring(1:3)) )
    attstring(1:3) = "   "
    attstring(1:7) = "alpha  "
    call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(18),&
 &      attstring(1:7), start=(/ 1, 1 /), count=(/ netcdf_ids(6), 1 /)) )
    attstring(1:7) = "beta   "
    call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(18),&
 &      attstring(1:7), start=(/ 1, 2 /), count=(/ netcdf_ids(6), 1 /)) )
    attstring(1:7) = "gamma  "
    call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(18),&
 &      attstring(1:7), start=(/ 1, 3 /), count=(/ netcdf_ids(6), 1 /)) )
    attstring(1:15) = "               "
    call check_fileionetcdf( nf90_sync(ncid) )
  end if
!
  if (ndump.lt.1) return ! should not happen
!
! set frame offset when appending
  framei = netcdf_fr1 + ndump
!
  if (netcdf_ids(11).le.0) call get_netcdfvarids(ncid)
  if (dyn_mode.eq.1) then
    timepos = 1.0*nstep
  else
    timepos = dyn_dt*nstep
  end if
  dims(1:2) = (/ 1,framei /)
  call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(11), (/ timepos /), count = dims(1:1), start = dims(2:2)) )
!  box matrix
  box(:,:) = 0.0
  shf = 0.0
  if (((pdb_rmol.gt.0).OR.(pdb_force_box.gt.1)).AND.(bnd_type.eq.1).AND.(bnd_shape.ne.1).AND.(bnd_shape.ne.3)) then
    write(ilog,*) 'Fatal. Encountered unsupported box shape for image adjustment &
 &in FMCSC_prtnetcdf(...). This is an omission bug (please report).'
    call fexit()
  end if
  if (bnd_shape.eq.1) then
    do i=1,3
      box(i,1) = bnd_params(i)
    end do
    do i=1,3
!     WARNING: change this to actual angles should non-rectangular boxes ever be supported
      box(i,2) = 90.0
    end do
  else if (bnd_shape.eq.2) then
!   proxy settings
    box(:,1) = 2.0*bnd_params(4)
    box(:,2) = 90.0
  else if (bnd_shape.eq.3) then
!   proxy settings
    do i=1,2
      box(i,1) = 2.0*bnd_params(4)
    end do
    box(3,1) = bnd_params(6)
    box(:,2) = 90.0
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in FMS&
 &MC_prtnetcdf(...) (code # ',bnd_shape,').'
    call fexit()
  end if
  call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(13), box(1:3,1), start = (/ 1,framei /) , count = (/ 3,1 /)) )
  call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(14), box(1:3,2), start = (/ 1,framei /) , count = (/ 3,1 /)) )
!
! now coords
  allocate(coords(3*pdbeffn))
!
  if (use_trajidx.EQV..true.) then
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      do i=1,pdbeffn
        mt = moltypid(molofrs(atmres(pdboutlst(i))))
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtnetcdf(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
        end if
        coords(i) = x(pdboutlst(i))
        coords(pdbeffn+i) = y(pdboutlst(i))
        coords(2*pdbeffn+i) = z(pdboutlst(i))
        if (bnd_shape.eq.1) then
          do k=1,3
            kk = (k-1)*pdbeffn + i
            if (coords(kk).lt.bnd_params(k+3)) then
              coords(kk) = coords(kk) + ceiling((bnd_params(3+k)-coords(kk))/bnd_params(k))*bnd_params(k)
            else if (coords(kk).gt.(bnd_params(k+3)+bnd_params(k))) then
              coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
            end if
          end do
        else if (bnd_shape.eq.3) then
          kk = 2*pdbeffn + i
          if (coords(kk).lt.(bnd_params(3)-0.5*bnd_params(6))) then
            coords(kk) = coords(kk) + ceiling((bnd_params(3)-0.5*bnd_params(6)-coords(kk))/bnd_params(6))*bnd_params(6)
          else if (coords(kk).gt.(bnd_params(3)+0.5*bnd_params(6))) then
            coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
          end if
        end if
        coords(2*pdbeffn+i) = coords(2*pdbeffn+i) + shf
      end do
    else
      lmol = 0
      do i=1,pdbeffn
        imol = molofrs(atmres(pdboutlst(i)))
        mt = moltypid(imol)
        if (pdb_rmol.gt.0) then
          if (lmol.ne.imol) then
            call shift_bound3(pdb_rmol,imol,tvec)
            shf3(:) = tvec(:) ! conversion
          else
            shf3(:) = tvec(:)
          end if
          lmol = imol
        else
          shf3(:) = 0.0
        end if
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtnetcdf(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
          shf3(3) = shf3(3) + shf
        end if
        coords(i) = x(pdboutlst(i))+shf3(1)
        coords(pdbeffn+i) = y(pdboutlst(i))+shf3(2)
        coords(2*pdbeffn+i) = z(pdboutlst(i))+shf3(3)
      end do
    end if
  else
    j = 0
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      do imol=1,nmol
        if ((just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
        mt = moltypid(imol)
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtnetcdf(...). This is an omission bug.'
                call fexit()
              end if
            end if
          end if
        end if
        do i=atmol(imol,1),atmol(imol,2)
          j = j + 1
          coords(j) = x(i)
          coords(pdbeffn+j) = y(i)
          coords(2*pdbeffn+j) = z(i)
          if (bnd_shape.eq.1) then
            do k=1,3
              kk = (k-1)*pdbeffn + j
              if (coords(kk).lt.bnd_params(k+3)) then
                coords(kk) = coords(kk) + ceiling((bnd_params(3+k)-coords(kk))/bnd_params(k))*bnd_params(k)
              else if (coords(kk).gt.(bnd_params(k+3)+bnd_params(k))) then
                coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
              end if
            end do
          else if (bnd_shape.eq.3) then
            kk = 2*pdbeffn + j
            if (coords(kk).lt.(bnd_params(3)-0.5*bnd_params(6))) then
              coords(kk) = coords(kk) + ceiling((bnd_params(3)-0.5*bnd_params(6)-coords(kk))/bnd_params(6))*bnd_params(6)
            else if (coords(kk).gt.(bnd_params(3)+0.5*bnd_params(6))) then
              coords(kk) = coords(kk) - ceiling((coords(kk)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
            end if
          end if
          coords(2*pdbeffn+j) = coords(2*pdbeffn+j) + shf
        end do
      end do
    else
      do imol=1,nmol
        if ((just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
        mt = moltypid(imol)
        if (pdb_rmol.gt.0) then
          call shift_bound3(pdb_rmol,imol,tvec)
          shf3(:) = tvec(:) ! conversion
        else
          shf3(:) = 0.0
        end if
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,mt).EQV..false.) then
            shf = 0.0
          else
            if (ismember(ispresent,imol).EQV..true.) then
              shf = 0.0
            else
              if (bnd_shape.eq.1) then
                shf = 2.0*bnd_params(3)
              else if (bnd_shape.eq.2) then
                shf = 4.0*bnd_params(4)
              else if (bnd_shape.eq.3) then
                shf = 2.0*bnd_params(6)
              else
                write(ilog,*) 'Fatal. Encountered unsupported box shape &
   &in FMSMC_prtnetcdf(...). This is an omission bug.'
                call fexit()
              end if
            end if
            shf3(3) = shf3(3) + shf
          end if
        end if
        do i=atmol(imol,1),atmol(imol,2)
          j = j + 1
          coords(j) = x(i)+shf3(1)
          coords(pdbeffn+j) = y(i)+shf3(2)
          coords(2*pdbeffn+j) = z(i)+shf3(3)
        end do
      end do
    end if
  end if
!
  call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(12), coords(1:pdbeffn), &
 &                                       start = (/ 1,1,framei /) , count = (/ 1,pdbeffn,1 /)) )
  call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(12), coords(pdbeffn+1:2*pdbeffn),&
 &                                        start = (/ 2,1,framei /) , count = (/ 1,pdbeffn,1 /)) )
  call check_fileionetcdf( nf90_put_var(ncid, netcdf_ids(12), coords(2*pdbeffn+1:3*pdbeffn),&
 &                                        start = (/ 3,1,framei /) , count = (/ 1,pdbeffn,1 /)) )
!
  deallocate(coords)
!
!
end
!
!-----------------------------------------------------------------------
!
subroutine get_netcdfvarids(ncid)
!
  use pdb
  use netcdf
!
  implicit none
!
  integer ncid
!
  call check_fileionetcdf( nf90_inq_varid(ncid,"time",netcdf_ids(11)) )
  call check_fileionetcdf( nf90_inq_varid(ncid,"coordinates",netcdf_ids(12)) )
  call check_fileionetcdf( nf90_inq_varid(ncid,"cell_lengths",netcdf_ids(13)) )
  call check_fileionetcdf( nf90_inq_varid(ncid,"cell_angles",netcdf_ids(14)) )
!
end
!
!---------------------------------------------------------------------------
!
subroutine check_fileionetcdf(ret)
!
  use iounit
  use netcdf
!
  implicit none
!
  integer ret
!
  if (ret.ne.NF90_NOERR) then
    write(ilog,*) 'Fatal. File I/O error in NetCDF operation. Returned error status is&
 &',ret,', which is explained as:'
    write(ilog,*) nf90_strerror(ret)
    write(ilog,*) 'By the following NetCDF library linked to CAMPARI:'
    write(ilog,*) nf90_inq_libvers()
! ' Reading from xtc-file (',netcdfinfile(t1:t2),')&
! &exited with an error (got ',ret,'). Sequence mismatch? Number of s&
! &napshots incorrect?'
!    ret2 = nf90_close(ncid)
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------------
!
subroutine check_netcdfappend(ncid,head,dumpfile)
!
  use iounit
  use netcdf
  use pdb
!
  implicit none
!
  integer ncid,i,framei,natomi,freeunit,dimlen,ret
  character(MAXSTRLEN) trystr
  character(*) dumpfile
  logical head
!
  natomi = -1
  framei = -1
  do i=1,NF90_MAX_DIMS
    ret = nf90_inquire_dimension(ncid,i,trystr,dimlen)
    if (ret.eq.NF90_NOERR) then
      call toupper(trystr(1:5))
      if (trystr(1:5).eq.'FRAME') then
        framei = dimlen
      else if (trystr(1:4).eq.'ATOM') then
        natomi = dimlen
      end if
    end if
  end do
!
  if ((natomi.lt.0).OR.(natomi.ne.pdbeffn)) then
    write(ilog,*) 'Warning. Existing NetCDF file is incompatible with new simulation &
 &and will be overwritten.'
    call check_fileionetcdf( nf90_close(ncid) )
    ncid = freeunit()
    call check_fileionetcdf( nf90_create(dumpfile, NF90_CLOBBER, ncid) )
    head  = .true.
    return
  end if
  if (framei.lt.0) then
    write(ilog,*) 'Warning. Existing NetCDF file is likely corrupt &
 &and will be overwritten.'
    call check_fileionetcdf( nf90_close(ncid) )
    ncid = freeunit()
    call check_fileionetcdf( nf90_create(dumpfile, NF90_CLOBBER, ncid) )
    head  = .true.
    return
  end if
!
  netcdf_fr1 = framei
!
end
!
#endif
!
!--------------------------------------------------------------------------------
!
! a wrapper to write a file with name [N_*_]{basename}_{stringy}.pdb to disk  
!
subroutine FMCSC_dump(stringy)
!
  use pdb
  use system
#ifdef ENABLE_MPI
  use mpistuff
#endif
  use cutoffs, ONLY: rsw1
!
  implicit none
!
  character(MAXSTRLEN), INTENT(IN):: stringy
!
  integer iint,freeunit,t1,t2,t3,t4,m1,buml
  character(MAXSTRLEN) intfile
  logical exists
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
!
  if (rsw1.le.0) then
!   this is a dirty monitor of whether we are far enough along in the calculation to safely dump a structure file
    return
  end if
!
  call strlims(stringy,t1,t2)
  m1 = -1
!
#ifdef ENABLE_MPI
  call int2str(myrank,nod,re_aux(10))
  intfile = 'N_'//nod//'_'//basename(1:bleng)//'_'//stringy(t1:t2)//'.int'
#else
  intfile = basename(1:bleng)//'_'//stringy(t1:t2)//'.int'
#endif
  call strlims(intfile,t3,t4)
  inquire(file=intfile(t3:t4),exist=exists)
  if (exists.EQV..true.) then
     iint = freeunit()
     open (unit=iint,file=intfile(t3:t4),status='old')
     close(unit=iint,status='delete')
  end if
  iint = freeunit()
  open (unit=iint,file=intfile(t3:t4),status='new')
  call prtzmat(iint)
  close(unit=iint)
!
! dump out a final .pdb-file (in case someone forgot to set xyzout-param)
#ifdef ENABLE_MPI
  call int2str(myrank,nod,re_aux(10))
  intfile = 'N_'//nod//'_'//basename(1:bleng)//'_'//stringy(t1:t2)//'.pdb'
#else
  intfile = basename(1:bleng)//'_'//stringy(t1:t2)//'.pdb'
#endif
  call strlims(intfile,t3,t4)
  inquire(file=intfile(t3:t4),exist=exists)
  if (exists.EQV..true.) then
     iint = freeunit()
     open (unit=iint,file=intfile(t3:t4),status='old')
     close(unit=iint,status='delete')
  end if
  iint = freeunit()
  open (unit=iint,file=intfile(t3:t4),status='new')
  buml = pdb_rmol
  pdb_rmol = 0
  call FMSMC_pdb(iint,m1)
  pdb_rmol = buml
  close(unit=iint)
!
end
!
!-----------------------------------------------------------------------------------
!
! the scripts are based more or less directly on work by VMD-afficionados:
! Axel Kohlmeyer (general), Andrew Dalke (sscache-stuff), Olaf Lenz (box stuff)
!
subroutine FMSMC_pdb_helper(ipdbh,ndump)
!
  use iounit
  use polypep
  use sequen
  use system
  use atoms
  use params
  use aminos
  use molecule
  use pdb
  use mcsums
  use grandensembles
  use mpistuff
  use math
!
  implicit none
!
  integer i,ipdbh,ndump,t2,t3,k,imol
  character(3) resname
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  character(MAXSTRLEN) itn,itn2
  RTYPE midc,r1,r2,r3
  integer haveprotein,havenucleic,cup(2),cun(2)
!
  haveprotein = 0
  havenucleic = 0
  cup(2) = 0
  cun(2) = 0
  do imol=1,nmol
    cun(1) = 0
    cup(1) = 0
    do i=rsmol(imol,1),rsmol(imol,2)
      if (seqpolty(i).eq.'P') then
        haveprotein = haveprotein + 1
        cup(1) = cup(1) + 1
      end if
      if (seqpolty(i).eq.'N') then
        havenucleic = havenucleic + 1
        cun(1) = cun(1) + 1
      end if
    end do
    if (cun(1).gt.cun(2)) cun(2) = cun(1)
    if (cup(1).gt.cup(2)) cup(2) = cup(1)
  end do
!
#ifdef ENABLE_MPI
  call int2str(myrank,nod,re_aux(10))
  if ((use_REMC.EQV..true.).OR.((use_MPIAVG.EQV..true.).AND.(force_singlexyz.EQV..true.))) then
    t2 = bleng + 3 + re_aux(10)
    itn(1:t2) = 'N_'//nod(1:re_aux(10))//'_'//basename(1:bleng)
  else
    itn(1:bleng) = basename(1:bleng)
    t2 = bleng
  end if
  t3 = bleng + 3 + re_aux(10)
  itn2(1:t3) = 'N_'//nod(1:re_aux(10))//'_'//basename(1:bleng)
#else
  itn(1:bleng) = basename(1:bleng)
  itn2(1:bleng) = basename(1:bleng)
  t2 = bleng
  t3 = bleng
#endif
!
  if (ndump.eq.1) then
    write(ipdbh,*) '## Credits for these scripts go to:'
    write(ipdbh,*) '## Axel Kohlmeyer'
    write(ipdbh,*) '## Andrew Dalke'
    write(ipdbh,*) '## Olaf Lenz'
    write(ipdbh,*)
    write(ipdbh,*) 'color Display Background gray'
    write(ipdbh,*) 'material change opacity Transparent 0.7'
    write(ipdbh,*) 'display culling off'
    write(ipdbh,*) 'display depthcue on'
    write(ipdbh,*) 'display cuedensity 0.3'
    write(ipdbh,*) 'display rendermode GLSL'
    write(ipdbh,*) 'axes location off'
    do i=1,MAXAMINO
      resname = amino(i)
      if ((aminopolty(i).eq.'P').OR.(aminopolty(i).eq.'N').OR.&
 & (resname.eq.'ACE').OR.(resname.eq.'FOR').OR.&
 &(resname.eq.'NME').OR.(resname.eq.'NH2')) then
        write(ipdbh,*) 'lappend rsnms ',resname
      end if
    end do
    do i=1,n_pdbunk
      k = pdb_unkbnd(i,3)
      resname = pdb_unknowns(i)
      if ((seqpolty(pdb_unkbnd(i,3)).eq.'P').OR.(seqpolty(pdb_unkbnd(i,3)).eq.'N').OR.&
 &        (rsmol(molofrs(pdb_unkbnd(i,3)),2).ne.rsmol(molofrs(pdb_unkbnd(i,3)),1))) then
        write(ipdbh,*) 'lappend rsnms ',resname
      end if
    end do

    write(ipdbh,*)
    write(ipdbh,*) 'set rep1 Licorice'
    write(ipdbh,*) 'set rep2 NewCartoon'
    write(ipdbh,*) 'set rep3 VdW'
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (bnd_shape.eq.1) then
        midc = bnd_params(6) + 1.5*bnd_params(3)
      else if (bnd_shape.eq.2) then
        midc = bnd_params(3) + 2.0*bnd_params(4)
      else if (bnd_shape.eq.3) then
        midc = bnd_params(3) + 1.0*bnd_params(6)
      else
        write(ilog,*) 'Fatal. Encountered unsupported box shape in F&
 &MSMC_pdb_helper(...). This is an omission bug.'
        call fexit()
      end if
      write(ipdbh,*) 'set sel1 "resname $rsnms and z < ',midc,'"'
      write(ipdbh,*) 'set sel2 "resname $rsnms and z < ',midc,&
 &'"'
      write(ipdbh,*) 'set sel3 "(not resname $rsnms) and (not hydrogen)\'
      write(ipdbh,*) ' and (z < ',midc,')"'
      write(ipdbh,*) 
      write(ipdbh,*) 'set selg "not hydrogen and z > ',midc,'"'
      write(ipdbh,*) 'set repg "VdW"'
      write(ipdbh,*)
    else        
      write(ipdbh,*) 'set sel1 "resname $rsnms"'
      write(ipdbh,*) 'set sel2 "resname $rsnms"'
      write(ipdbh,*) 'set sel3 "(not resname $rsnms) and (not hydrogen)"'
    end if
    write(ipdbh,*) 'set trjsm 0'
    write(ipdbh,*)
    write(ipdbh,*) 'if { [glob -nocomplain ',itn2(1:t3),'_END.pdb ] \'
    write(ipdbh,*) '== "',itn2(1:t3),'_END.pdb" } then \'
    write(ipdbh,*) '{set ftemplate "',itn2(1:t3),'_END.pdb"} \'
    write(ipdbh,*) 'else {set ftemplate "',itn2(1:t3),'_START.pdb"}'
    write(ipdbh,*) 'mol load pdb $ftemplate'
    if (xyzmode.eq.1) then
      if (mod(nsim - nequil,xyzout).eq.0) then
        k = (nsim-nequil)/xyzout
      else
        k = floor(1.0*(nsim - nequil)/(1.0*xyzout)) + 1
      end if
      write(ipdbh,*) 'for {set i ',ndump,'} {$i<=',k,'} {incr i 1} {'
      write(ipdbh,*) '  set fileName [format ',itn(1:t2),'_%05d.arc $i]'
      write(ipdbh,*) '  mol addfile $fileName waitfor all'
      write(ipdbh,*) '}'
    else if (xyzmode.eq.2) then
      if (pdb_writemode.eq.2) then
        write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.pdb waitfor all'
      else if (pdb_writemode.eq.1) then
        if (mod(nsim - nequil,xyzout).eq.0) then
          k = (nsim-nequil)/xyzout
        else
          k = floor(1.0*(nsim - nequil)/(1.0*xyzout)) + 1
        end if
        write(ipdbh,*) 'for {set i ',ndump,'} {$i<=',k,'} {incr i 1} {'
        write(ipdbh,*) '  set fileName [format ',itn(1:t2),'_%05d.pdb $i]'
        write(ipdbh,*) '  mol addfile $fileName waitfor all'
        write(ipdbh,*) '}'
      end if
    else if (xyzmode.eq.3) then
      write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.dcd waitfor all'
    else if (xyzmode.eq.4) then
      write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.xtc waitfor all'
    else if (xyzmode.eq.5) then
      write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.nc waitfor all'
    end if
    write(ipdbh,*) 'mol modcolor 0 0 Type'
    write(ipdbh,*) 'mol modstyle 0 0 $rep1'
    write(ipdbh,*) 'mol modmaterial 0 0 Glossy'
    write(ipdbh,*) 'mol modselect 0 0 $sel1'
    write(ipdbh,*) 'mol addrep 0'
    write(ipdbh,*) 'mol modstyle 1 0 $rep2'
    write(ipdbh,*) 'mol modcolor 1 0 Index'
    write(ipdbh,*) 'mol modselect 1 0 $sel2'
    write(ipdbh,*) 'mol modmaterial 1 0 Transparent'
    write(ipdbh,*) 'mol addrep 0'
    write(ipdbh,*) 'mol modstyle 2 0 $rep3'
    write(ipdbh,*) 'mol modcolor 2 0 Type'
    write(ipdbh,*) 'mol modselect 2 0 $sel3'
    write(ipdbh,*) 'mol modmaterial 2 0 Glossy'
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      write(ipdbh,*) 'mol addrep 0'
      write(ipdbh,*) 'mol modstyle 3 0 $repg'
      write(ipdbh,*) 'mol modcolor 3 0 Type'
      write(ipdbh,*) 'mol modmaterial 3 0 BrushedMetal'
      write(ipdbh,*) 'mol modselect 3 0 $selg'
    end if
    write(ipdbh,*)
    write(ipdbh,*) 'set nreps [molinfo 0 get numreps]'
    write(ipdbh,*) 'set nfrms [molinfo 0 get numframes]'
    write(ipdbh,*) 
    write(ipdbh,*) 'for {set i 0} {$i < $nreps} {incr i} {'
    write(ipdbh,*) '  mol selupdate $i 0 on'
    write(ipdbh,*) '  mol smoothrep 0 $i $trjsm'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc start_sscache {{molid 0}} {'
    write(ipdbh,*) '  global sscache_data vmd_frame'
    write(ipdbh,*) '  trace variable vmd_frame($molid) w sscache'
    write(ipdbh,*) '  return'
    write(ipdbh,*) '}'
    write(ipdbh,*) 'proc stop_sscache {{molid 0}} {'
    write(ipdbh,*) '  global vmd_frame'
    write(ipdbh,*) '  trace vdelete vmd_frame($molid) w sscache'
    write(ipdbh,*) '  return'
    write(ipdbh,*) '}'
    write(ipdbh,*) 'proc reset_sscache {} {'
    write(ipdbh,*) '  global sscache_data'
    write(ipdbh,*) '  if [info exists sscache_data] {'
    write(ipdbh,*) '    unset sscache_data'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  return'
    write(ipdbh,*) '}'
    write(ipdbh,*) 'proc sscache {name index op} {'
    write(ipdbh,*) '  global sscache_data'
    write(ipdbh,*) '  set sel [atomselect $index "protein name CA"]'
    write(ipdbh,*) '  set frame [molinfo $index get frame]'
    write(ipdbh,*) '  if [info exists sscache_data($index,$frame)] {'
    write(ipdbh,*) '    $sel set structure $sscache_data($index,$frame)'
    write(ipdbh,*) '  return'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  vmd_calculate_structure $index'
    write(ipdbh,*) '  set sscache_data($index,$frame) [$sel get structure]'
    write(ipdbh,*) '  return'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    if ((havenucleic.eq.0).AND.(haveprotein.eq.0)) then
      write(ipdbh,*) 'proc align_it {{sel "all"} {rfmol 0} {reffr 0} {molid 0}} {'
    else if (havenucleic.eq.0) then
      write(ipdbh,*) 'proc align_it {{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
    else if (haveprotein.eq.0) then
      write(ipdbh,*) 'proc align_it {{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
    else if (ceiling(1.4*cun(2)).ge.cup(2)) then
      write(ipdbh,*) 'proc align_it {{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
    else
      write(ipdbh,*) 'proc align_it {{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
    end if
    write(ipdbh,*) '  set alignref [atomselect $rfmol $sel frame $reffr]'
    write(ipdbh,*) '  set nf [molinfo $molid get numframes]'
    write(ipdbh,*) '  for {set i 0} {$i < $nf} {incr i} {'
    write(ipdbh,*) '    set match [atomselect $molid $sel frame $i]'
    write(ipdbh,*) '    set tm [measure fit $match $alignref]'
    write(ipdbh,*) '    set alls [atomselect $molid "all" frame $i]'
    write(ipdbh,*) '    $alls move $tm'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc rmsd_it {{filen \'
    write(ipdbh,*) itn(1:t2),'_RMSD.dat} \'
    if ((havenucleic.eq.0).AND.(haveprotein.eq.0)) then
      write(ipdbh,*) '{sel "all"} {rfmol 0} {reffr 0} {molid 0}} {'
    else if (havenucleic.eq.0) then
      write(ipdbh,*) '{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
    else if (haveprotein.eq.0) then
      write(ipdbh,*) '{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
    else if (ceiling(1.4*cun(2)).ge.cup(2)) then
      write(ipdbh,*) '{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
    else
      write(ipdbh,*) '{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
    end if
    write(ipdbh,*) '  set fi [open $filen w]'
    write(ipdbh,*) '  set alignref [atomselect $rfmol $sel frame $reffr]'
    write(ipdbh,*) '  set refall [atomselect $rfmol "all" frame $reffr]'
    write(ipdbh,*) '  set nf [molinfo $molid get numframes]'
    write(ipdbh,*) '  for {set i 1} {$i < $nf} {incr i} {'
    write(ipdbh,*) '    set match [atomselect $molid $sel frame $i]'
    write(ipdbh,*) '    set tm [measure fit $match $alignref]'
    write(ipdbh,*) '    set alls [atomselect $molid "all" frame $i]'
    write(ipdbh,*) '    $alls move $tm'
    write(ipdbh,*) '    set match [atomselect $molid $sel frame $i]'
    write(ipdbh,*) '    set irmsd [measure rmsd $match $alignref]'
    write(ipdbh,*) '    set alls [atomselect $molid "all" frame $i]'
    write(ipdbh,*) '    set irmsd2 [measure rmsd $alls $refall]'
    write(ipdbh,*) '    puts $fi "$i  $irmsd  $irmsd2"'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  close $fi'
    write(ipdbh,*) '}'
    write(ipdbh,*)
   72 format('  lappend gids [graphics 0 line \')
   68 format('  "',a4,1x,a4,1x,a4,'"\')
   69 format('  "',a4,1x,a4,1x,a4,'" width 2]')
   73 format('  lappend gids [graphics 0 sphere \')
   79 format('  lappend gids [graphics 0 cylinder \')
   70 format('  "',a4,1x,a4,1x,a4,'"\')
   71 format('  radius "',a4,'" resolution "',a4,'"]')
   80 format('  radius "',a4,'" resolution "',a4,'" filled "yes"]')
   74 format('  return $gids')
    if (bnd_shape.eq.1) then
      write(ipdbh,*) 'proc see_box {{la 50} {lb 50} {lc 50} \'
      write(ipdbh,*) ' {aa 90.0} {ab 90.0} {ac 90.0} \'
      write(ipdbh,*) ' {ooa 0.0} {oob 0.0} {ooc 0.0}} {'
!      write(ipdbh,*) '  graphics 0 delete all'
!      write(ipdbh,*) '  graphics 0 color yellow'
      write(ipdbh,*) '  set radian [expr "180.0/3.14159265"];'
      write(ipdbh,*) '  set cax [expr "$ooa+$la"];'
      write(ipdbh,*) '  set cbx [expr "$ooa+$lb*cos($ac/$radian)"];'
      write(ipdbh,*) '  set cby [expr "$oob+$lb*sin($ac/$radian)"];'
      write(ipdbh,*) '  set ccx [expr "$ooa+$lc*cos($ab/$radian)"];'
      write(ipdbh,*) '  set ccy [expr "$oob+$lc*cos($aa/$radian)*sin($ab/$radian)"];'
      write(ipdbh,*) '  set ccz [expr "$ooc+$lc*sin($aa/$radian)*sin($ab/$radian)"];'
      write(ipdbh,72)
      write(ipdbh,68) '$ooa','$oob','$ooc'
      write(ipdbh,69) '$cax','$oob','$ooc'
      write(ipdbh,72)
      write(ipdbh,68) '$ooa','$oob','$ooc'
      write(ipdbh,69) '$cbx','$cby','$ooc'
      write(ipdbh,72)
      write(ipdbh,68) '$ooa','$oob','$ooc'
      write(ipdbh,69) '$ccx','$ccy','$ccz'
      write(ipdbh,*) '  set dbx [expr "$cax+$lb*cos($ac/$radian)"];'
      write(ipdbh,*) '  set dby [expr "$oob+$lb*sin($ac/$radian)"];'
      write(ipdbh,*) '  set dcx [expr "$cax+$lc*cos($ab/$radian)"];'
      write(ipdbh,*) '  set dcy [expr "$oob+$lc*cos($aa/$radian)*sin($ab/$radian)"];'
      write(ipdbh,*) '  set dcz [expr "$ooc+$lc*sin($aa/$radian)*sin($ab/$radian)"];'
      write(ipdbh,72)
      write(ipdbh,68) '$cax','$oob','$ooc'
      write(ipdbh,69) '$dbx','$dby','$ooc'
      write(ipdbh,72)
      write(ipdbh,68) '$cax','$oob','$ooc'
      write(ipdbh,69) '$dcx','$dcy','$dcz'
      write(ipdbh,*) '  set eax [expr "$cbx+$la"];'
      write(ipdbh,*) '  set ecx [expr "$cbx+$lc*cos($ab/$radian)"];'
      write(ipdbh,*) '  set ecy [expr "$cby+$lc*cos($aa/$radian)*sin($ab/$radian)"];'
      write(ipdbh,*) '  set ecz [expr "$ooc+$lc*sin($aa/$radian)*sin($ab/$radian)"];'
      write(ipdbh,72)
      write(ipdbh,68) '$cbx','$cby','$ooc'
      write(ipdbh,69) '$eax','$cby','$ooc'
      write(ipdbh,72)
      write(ipdbh,68) '$cbx','$cby','$ooc'
      write(ipdbh,69) '$ecx','$ecy','$ecz'
      write(ipdbh,*) '  set fax [expr "$ccx+$la"];'
      write(ipdbh,*) '  set fbx [expr "$ccx+$lb*cos($ac/$radian)"];'
      write(ipdbh,*) '  set fby [expr "$ccy+$lb*sin($ac/$radian)"];'
      write(ipdbh,72)
      write(ipdbh,68) '$ccx','$ccy','$ccz'
      write(ipdbh,69) '$fax','$ccy','$ccz'
      write(ipdbh,72)
      write(ipdbh,68) '$ccx','$ccy','$ccz'
      write(ipdbh,69) '$fbx','$fby','$ccz'
      write(ipdbh,*) '  set tax [expr "$dbx+$dcx-$cax"];'
      write(ipdbh,*) '  set tay [expr "$cby+$ccy-$oob"];'
      write(ipdbh,*) '  set taz [expr "$ccz"];'
      write(ipdbh,72)
      write(ipdbh,68) '$fbx','$fby','$ccz'
      write(ipdbh,69) '$tax','$tay','$taz'
      write(ipdbh,72)
      write(ipdbh,68) '$dcx','$dcy','$dcz'
      write(ipdbh,69) '$tax','$tay','$taz'
      write(ipdbh,72)
      write(ipdbh,68) '$dbx','$dby','$ooc'
      write(ipdbh,69) '$tax','$tay','$taz'
      write(ipdbh,74) 
      write(ipdbh,*) '}'
    else if (bnd_shape.eq.2) then
      write(ipdbh,*) 'proc see_box {{rad 50.0} {ooa 0.0} \'
      write(ipdbh,*) '{oob 0.0} {ooc 0.0} {res 200}} {'
!      write(ipdbh,*) '  graphics 0 delete all'
!      write(ipdbh,*) '  graphics 0 color ochre'
!      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set radian [expr "180.0/3.14159265"];'
      write(ipdbh,73)
      write(ipdbh,70) '$ooa','$oob','$ooc'
      write(ipdbh,71) '$rad','$res'
      write(ipdbh,74) 
      write(ipdbh,*) '}'
    else if (bnd_shape.eq.3) then
      write(ipdbh,*) 'proc see_box {{rad 50.0} {hght 50.0} {ooa 0.0} \'
      write(ipdbh,*) '{oob 0.0} {ooc 0.0} {res 200}} {'
      write(ipdbh,*) '  set radian [expr "180.0/3.14159265"];'
      write(ipdbh,*) '  set cax [expr "$ooc+$hght"];'
      write(ipdbh,79)
      write(ipdbh,70) '$ooa','$oob','$ooc'
      write(ipdbh,70) '$ooa','$oob','$cax'
      write(ipdbh,80) '$rad','$res'
      write(ipdbh,74) 
      write(ipdbh,*) '}'
    end if
  end if
  write(ipdbh,*)
!
   63 format('  see_box ',3(f9.3,1x),' \')
   64 format(3(f9.3,1x),' \')
   65 format(3(f9.3,1x))
   60 format(1(f9.3,1x),' \')
   681 format(2(f9.3,1x),' \')
   66 format('  see_box ',f9.3,' \')
   67 format(3(f9.3,1x),' 200')
   62 format('  get a frame $frn" ] ',3(f9.3,1x),'"')
   61 format('  ',3(f9.3,1x),'"')
  if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
    write(ipdbh,*) 'proc see_system {{molid 0}} {'
    write(ipdbh,*) '  global box_gids'
    write(ipdbh,*) '  if {[info exists box_gids]} then {'
    write(ipdbh,*) '    foreach g $box_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    if (bnd_shape.eq.1) then
      write(ipdbh,*) '  graphics 0 color yellow'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)
      write(ipdbh,64) 90.0,90.0,90.0 
      write(ipdbh,64) bnd_params(4),bnd_params(5),bnd_params(6)
      write(ipdbh,*) '"'
    else if (bnd_shape.eq.2) then
      write(ipdbh,*) '  graphics 0 color ochre'
      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,60) bnd_params(4)
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)
      write(ipdbh,*) '"'
    else if (bnd_shape.eq.3) then
      write(ipdbh,*) '  graphics 0 color ochre'
      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,681) bnd_params(4),bnd_params(6)
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)-0.5*bnd_params(6)
      write(ipdbh,*) '"'
    end if
    write(ipdbh,*) '  set box_gids [eval "see_box $box_prms" ]'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc hide_system {{molid 0}} {'
    write(ipdbh,*) '  global box_gids'
    write(ipdbh,*) '  if {[info exists box_gids]} then {'
    write(ipdbh,*) '    foreach g $box_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '}'
    write(ipdbh,*)
  else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
    write(ipdbh,*) 'set molid 0'
    write(ipdbh,*)
    write(ipdbh,*) 'proc box_draw {args} {'
    if (bnd_shape.eq.1) then
      write(ipdbh,*) '  global vmd_frame box_gids box_prms molid'
      write(ipdbh,*) '  set frn [molinfo $molid get frame]'
      write(ipdbh,*) '  set box_prms($frn) "[eval "molinfo $molid \'
      write(ipdbh,*) 'get \"a b c alpha beta gamma\" frame $frn" ]\'
      write(ipdbh,61) bnd_params(4),bnd_params(5),bnd_params(6)
      write(ipdbh,*) '  if {[info exists box_prms($frn)]} then {'
      write(ipdbh,*) '    if {[info exists box_gids]} then {'
      write(ipdbh,*) '      foreach g $box_gids {'
      write(ipdbh,*) '        graphics $molid delete $g'
      write(ipdbh,*) '      }'
      write(ipdbh,*) '    }'
      write(ipdbh,*) '    graphics 0 color ochre'
      write(ipdbh,*) '    graphics 0 material Glass1'
      write(ipdbh,*) '    set box_gids [eval "see_box $box_prms($frn)" ]'
      write(ipdbh,*) '  }'
    else if (bnd_shape.eq.2) then
      write(ipdbh,*) '  global vmd_frame box_gids box_prms molid'
      write(ipdbh,*) '  set frn [molinfo $molid get frame]'
      write(ipdbh,*) '  set box_prms($frn) "[eval "molinfo $molid \'
      write(ipdbh,62) bnd_params(1),bnd_params(2),bnd_params(3)
      write(ipdbh,*) '  if {[info exists box_prms($frn)]} then {'
      write(ipdbh,*) '    if {[info exists box_gids]} then {'
      write(ipdbh,*) '      foreach g $box_gids {'
      write(ipdbh,*) '        graphics $molid delete $g'
      write(ipdbh,*) '      }'
      write(ipdbh,*) '    }'
      write(ipdbh,*) '    graphics 0 color ochre'
      write(ipdbh,*) '    graphics 0 material Glass1'
      write(ipdbh,*) '    set box_gids [eval "see_box $box_prms($frn)" ]'
      write(ipdbh,*) '  }'
    end if
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc see_system {args} {'
    write(ipdbh,*) '  global vmd_frame box_gids box_prms molid'
    write(ipdbh,*) '  trace variable vmd_frame($molid) w box_draw'
    write(ipdbh,*) '  animate goto 0'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc hide_system {args} {'
    write(ipdbh,*) '  global vmd_frame box_gids box_prms molid'
    write(ipdbh,*) '  if {[info exists box_gids]} then {'
    write(ipdbh,*) '    foreach g $box_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  trace vdelete vmd_frame($molid) w box_draw'
    write(ipdbh,*) '  animate goto 0'
    write(ipdbh,*) '}'
    write(ipdbh,*)
  else if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    write(ipdbh,*) 'proc see_system {{molid 0}} {'
    write(ipdbh,*) '  global box_gids'
    write(ipdbh,*) '  if {[info exists box_gids]} then {'
    write(ipdbh,*) '    foreach g $box_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  set nreps [molinfo 0 get numreps]'
    write(ipdbh,*) '  set butone [expr $nreps-1]'
    write(ipdbh,*) '  for {set i 0} {$i < $butone} {incr i} {'
    write(ipdbh,*) '    mol showrep 0 $i on'
    write(ipdbh,*) '  }'
    if (bnd_shape.eq.1) then
      write(ipdbh,*) '  graphics 0 color yellow'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)
      write(ipdbh,64) 90.0,90.0,90.0 
      write(ipdbh,64) bnd_params(4),bnd_params(5),bnd_params(6)
      write(ipdbh,*) '"'
    else if (bnd_shape.eq.2) then
      write(ipdbh,*) '  graphics 0 color ochre'
      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,60) bnd_params(4)
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)
      write(ipdbh,*) '"'
    else if (bnd_shape.eq.3) then
      write(ipdbh,*) '  graphics 0 color ochre'
      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,681) bnd_params(4),bnd_params(6)
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)-0.5*bnd_params(6)
      write(ipdbh,*) '"'
    end if
    write(ipdbh,*) '  set box_gids [eval "see_box $box_prms" ]'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc hide_system {{molid 0}} {'
    write(ipdbh,*) '  global box_gids'
    write(ipdbh,*) '  if {[info exists box_gids]} then {'
    write(ipdbh,*) '    foreach g $box_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  set nreps [molinfo 0 get numreps]'
    write(ipdbh,*) '  set butone [expr $nreps-1]'
    write(ipdbh,*) '  for {set i 0} {$i < $butone} {incr i} {'
    write(ipdbh,*) '    mol showrep 0 $i off'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc see_bath {{molid 0}} {'
    write(ipdbh,*) '  global bth_gids'
    write(ipdbh,*) '  if {[info exists bth_gids]} then {'
    write(ipdbh,*) '    foreach g $bth_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  set nreps [molinfo 0 get numreps]'
    write(ipdbh,*) '  set butone [expr $nreps-1]'
    write(ipdbh,*) '  mol showrep 0 $butone on'
    if (bnd_shape.eq.1) then
      write(ipdbh,*) '  graphics 0 color yellow'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,64) bnd_params(1),bnd_params(2),bnd_params(3)
      write(ipdbh,64) 90.0,90.0,90.0 
      r1 = bnd_params(4)
      r2 = bnd_params(5)
      r3 = bnd_params(6) + 2.0*bnd_params(3)
      write(ipdbh,64) r1,r2,r3
      write(ipdbh,*) '"'
    else if (bnd_shape.eq.2) then
      write(ipdbh,*) '  graphics 0 color ochre'
      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,60) bnd_params(4)
      r1 = bnd_params(1)
      r2 = bnd_params(2)
      r3 = bnd_params(3) + 4.0*bnd_params(4)
      write(ipdbh,64) r1,r2,r3
      write(ipdbh,*) '"'
    else if (bnd_shape.eq.3) then
      write(ipdbh,*) '  graphics 0 color ochre'
      write(ipdbh,*) '  graphics 0 material Glass1'
      write(ipdbh,*) '  set box_prms "\'
      write(ipdbh,681) bnd_params(4),bnd_params(6)
      r1 = bnd_params(1)
      r2 = bnd_params(2)
      r3 = bnd_params(3) + 1.5*bnd_params(6)
      write(ipdbh,64) r1,r2,r3
      write(ipdbh,*) '"'
    end if
    write(ipdbh,*) '  set bth_gids [eval "see_box $box_prms" ]'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc hide_bath {{molid 0}} {'
    write(ipdbh,*) '  global bth_gids'
    write(ipdbh,*) '  if {[info exists bth_gids]} then {'
    write(ipdbh,*) '    foreach g $bth_gids {'
    write(ipdbh,*) '      graphics $molid delete $g'
    write(ipdbh,*) '    }'
    write(ipdbh,*) '  }'
    write(ipdbh,*) '  set nreps [molinfo 0 get numreps]'
    write(ipdbh,*) '  set butone [expr $nreps-1]'
    write(ipdbh,*) '  mol showrep 0 $butone off'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'proc see_all {{molid 0}} {'
    write(ipdbh,*) '  global box_gids bth_gids'
    write(ipdbh,*) '  see_system $molid'
    write(ipdbh,*) '  see_bath $molid'
    write(ipdbh,*) '}'
    write(ipdbh,*)
    write(ipdbh,*) 'see_all 0'
    write(ipdbh,*)
  end if
  write(ipdbh,*) 'start_sscache 0'
  write(ipdbh,*)
!
end
!
!-----------------------------------------------------------------------
!
! this subroutine writes atomic coordinates in TINKER's xyz-format
!
subroutine tinkerxyz(iu)
!
  use atoms
  use molecule
  use grandensembles
  use sequen
  use iounit
  use system
  use pdb
!
  implicit none
!
  integer i,j,k,iu,mt,imol,offset
  logical isopen,ismember
  RTYPE shf,outv(3),shf3(3),tvec(3)
!
!
 10   format (i6)
 30   format (i6,2x,a3,3f12.6,9i6)
!
! open output unit if not already done
!
  inquire (unit=iu,opened=isopen)
  if (isopen.EQV..false.) then
    write(ilog,*) 'Fatal. Got stale filehandle in tinkerxyz(...).'
    call fexit()
  end if
!
  shf = 0.0
  if (((pdb_rmol.gt.0).OR.(pdb_force_box.gt.1)).AND.(bnd_type.eq.1).AND.(bnd_shape.ne.1).AND.(bnd_shape.ne.3)) then
    write(ilog,*) 'Fatal. Encountered unsupported box shape for image adjustment &
 &in tinkerxyz(...). This is an omission bug (please report).'
    call fexit()
  end if
  write (iu,10) pdbeffn
  j = 0
  do imol=1,nmol
    if ((use_trajidx.EQV..false.).AND.(just_solutes.EQV..true.).AND.(is_solvent(imol).EQV..true.)) cycle
    offset = atmol(imol,1) - (j+1)
    if (pdb_rmol.gt.0) then
      call shift_bound3(pdb_rmol,imol,tvec)
    else
      tvec(:) = 0.0
    end if
    do i=atmol(imol,1),atmol(imol,2)
      shf3(:) = tvec(:)
      if (use_trajidx.EQV..true.) then
        if (pdboutlog(i).EQV..false.) cycle
      end if
      j = j + 1
      mt = moltypid(imol)
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,mt).EQV..false.) then
          shf = 0.0
        else
          if (ismember(ispresent,imol).EQV..true.) then
            shf = 0.0
          else
            if (bnd_shape.eq.1) then
              shf = 2.0*bnd_params(3)
            else if (bnd_shape.eq.2) then
              shf = 4.0*bnd_params(4)
            else if (bnd_shape.eq.3) then
              shf = 2.0*bnd_params(6)
            else
              write(ilog,*) 'Fatal. Encountered unsupported box shape &
 &in tinkerxyz(...). This is an omission bug.'
              call fexit()
            end if
          end if
        end if
        shf3(3) = shf3(3) + shf
      end if
      if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
        outv(1) = x(i)
        outv(2) = y(i)
        outv(3) = z(i)
        if (bnd_shape.eq.1) then
          do k=1,3  
            if (outv(k).lt.bnd_params(k+3)) then
              outv(k) = outv(k) + ceiling((bnd_params(3+k)-outv(k))/bnd_params(k))*bnd_params(k)
            else if (outv(k).gt.(bnd_params(k+3)+bnd_params(k))) then
              outv(k) = outv(k) - ceiling((outv(k)-bnd_params(3+k)-bnd_params(k))/bnd_params(k))*bnd_params(k)
            end if
          end do
        else if (bnd_shape.eq.3) then
          if (outv(3).lt.bnd_params(3)-0.5*bnd_params(6)) then
            outv(3) = outv(3) + ceiling((bnd_params(3)-0.5*bnd_params(6)-outv(3))/bnd_params(6))*bnd_params(6)
          else if (outv(3).gt.(bnd_params(3)+0.5*bnd_params(6))) then
            outv(3) = outv(3) - ceiling((outv(3)-bnd_params(3)-0.5*bnd_params(6))/bnd_params(6))*bnd_params(6)
          end if
        end if
        write (iu,30) j,atnam(i),outv(1),outv(2),outv(3)+shf,attyp(i),(i12(k,i)-offset,k=1,n12(i))
      else
        write (iu,30) j,atnam(i),x(i)+shf3(1),y(i)+shf3(2),z(i)+shf3(3),attyp(i),(i12(k,i)-offset,k=1,n12(i))
      end if
    end do
  end do
!
end
!
!-------------------------------------------------------------------
!
