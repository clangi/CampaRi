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
! CONTRIBUTIONS: Rohit Pappu, Hoang Tran                                   !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------------
!
! this routine parses the sequence and then sets up the simulated molecules
!
!------------------------------------------------------------------------------
!
subroutine makepept()
!
  use fyoc
  use iounit
  use sequen
  use molecule
  use aminos
  use zmatrix
  use atoms
  use polypep
  use system
  use pdb
  use cutoffs
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer rs,imol,jj,iomess
  integer freeunit,iseqfile,aone,atwo,sickle,athree
  integer j,k,next,length,i,t1,t2,p1,p2
  RTYPE rr,random!,getbang,getztor,getblen
  character(60) resname
  character(MAXSTRLEN) seqlin
  logical done,startmol,foundN,foundC,foundM,foundit,laststartmol,exists
!
! a little bit of init
  aone = 1
  atwo = 2
  athree = 3
!
! read the parameter file
  call readprm()
!
! open sequence input file
  call strlims(seqfile,t1,t2)
  inquire(file=seqfile(t1:t2),exist=exists)  !input ascii file
  if (exists.EQV..false.) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Sequence file ',seqfile(t1:t2),' not found.'
    call fexit()
  end if
  iseqfile = freeunit()
  open(unit=iseqfile,file=seqfile(t1:t2),status='old')
!
! first do a dry run to get the limits for the calculation
! also check for errors
!
! get the peptide sequence
  i = 0
  n_pdbunk = 0
  done = .false.
  startmol = .false.
  sickle = -1
  n_crosslinks = 0
  do while (done.EQV..false.)
!
     foundN = .false.
     foundC = .false.
     foundM = .false.
     laststartmol = startmol
     if (sickle.ne.-1) then ! last res. was CYS
       jj = 0
       read(seqlin(next:MAXSTRLEN),49,iostat=iomess) jj
       if ((iomess.eq.0).AND.(jj.gt.0).AND.(jj.ne.i)) then
         if (i.lt.jj) n_crosslinks = n_crosslinks + 1
       end if
       do k=1,MAXSTRLEN
         seqlin(k:k) = ' '
       end do
     end if

     sickle = -1
     foundit = .false.
!
 49      format(i100)
 50      format(A)
!
     read(iseqfile,50,iostat=iomess) seqlin
     if (iomess.eq.IOSTAT_END) then
       close(unit=iseqfile)
       exit
     else if (iomess.eq.2) then
       write(ilog,*) 'Fatal. File I/O error while processing sequenc&
 &e input file (got: ',seqfile(t1:t2),').'
       call fexit()
     end if
     call toupper(seqlin)
     next = 1
     call extract_str(seqlin,resname,next)
     call strlims(resname,p1,p2)
     length = p2-p1+1
!
!    first let's check for 3-letter-or-less codes
     if ((length.eq.3).OR.(length.eq.2)) then
        if (resname.eq.'END') then
          if (startmol.EQV..true.) then
            write(ilog,*) 'Fatal. Encountered end of sequence input in s&
 &equence file without end of previous molecule.'
            call fexit()
          end if
          done = .true.
          exit
        end if
        if (length.eq.2) then
          resname(3:3) = ' '
        end if
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            sickle = j
            foundit = .true.
          end if
        end do
!       N-terminal caps
        if ((foundit.EQV..true.).AND.((sickle.eq.27).OR.(sickle.eq.28).OR.&
 &          ((sickle.ge.76).AND.(sickle.le.87)))) then
          foundN = .true.
!       C-terminal caps
        else if ((foundit.EQV..true.).AND.((sickle.eq.29).OR.(sickle.eq.30))) then
          foundC = .true.
!       standard in-chain residues
        else if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &               ((sickle.ge.31).AND.(sickle.le.35)).OR.&
 &               (sickle.eq.51).OR.&
 &               ((sickle.ge.64).AND.(sickle.le.75)).OR.((sickle.ge.108).AND.(sickle.le.119)))) then
          foundM = .true.
!       small molecules and atoms
        else if ((foundit.EQV..true.).AND.(((sickle.ge.36).AND.(sickle.le.50)).OR.&
 &               ((sickle.ge.52).AND.(sickle.le.63)).OR.&
 &               ((sickle.ge.88).AND.(sickle.le.107)))) then
          foundN = .true.  
          foundC = .true.
!       unknown residues (have to always be explicitly N/C labeled)
        else ! if (pdb_analyze.EQV..true.) then
          foundM = .true.
          i = i + 1
        end if
!
!       now see what we've got
        if ((foundN.EQV..true.).AND.(foundC.EQV..true.)) then
          if (startmol.EQV..true.) then
            write(ilog,*) 'Fatal. Found start of new molecule in seq&
 &uence file without end of previous molecule (',resname(1:3),').'
            call fexit()
          end if
          nmol = nmol + 1
        else if (foundN.EQV..true.) then
          if (startmol.EQV..true.) then
            write(ilog,*) 'Fatal. Found start of new molecule in seq&
 &uence file without end of previous molecule (',resname(1:3),').'
            call fexit()
          end if
          nmol = nmol + 1
          startmol = .true.
        else if (foundC.EQV..true.) then
          if (startmol.EQV..false.) then
            write(ilog,*) 'Fatal. Found end of current molecule in s&
 &equence file without ever starting it (',resname(1:3),').'
            call fexit()
          end if
          startmol = .false.
        else if (foundM.EQV..true.) then
          if (startmol.EQV..false.) then
            write(ilog,*) 'Fatal. Found intra-chain residue without &
 &ever starting a chain in seq. file (',resname(1:3),').'
            call fexit()
          end if
        else
          write(ilog,*) 'Fatal. Found an unknown residue in seq. fil&
 &e, which is currently not supported (',resname(1:3),').'
          call fexit()
        end if
!
!
!    the 5-letter codes are really 3-letter codes plus an appendix (_D, _C, _N)
     else if (length.eq.5) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            sickle = j
            foundit = .true.
          end if
        end do
!
!       parse for appendix
        if ((resname(4:4).ne.'_').OR.((resname(5:5).ne.'D').AND.&
 &           (resname(5:5).ne.'N').AND.(resname(5:5).ne.'C'))) then
          write(ilog,*) 'Fatal. Found an unknown residue in seq. fil&
 &e, which is currently not supported (',resname(1:5),').'
          call fexit()
        end if
!
!       chirality request
        if (resname(5:5).eq.'D') then
!         D-aa in-chain residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.2).AND.(sickle.le.23)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.119)).OR.&
 &            (sickle.eq.25).OR.(sickle.eq.51))) then
            foundM = .true.              
          end if
!
          if (foundM.EQV..true.) then
            if (startmol.EQV..false.) then
              write(ilog,*) 'Fatal. Found intra-chain residue withou&
 &t ever starting a chain in seq. file (',resname(1:5),').'
              call fexit()
            end if
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible chirality request in seq. file (',resname(1:5),').'
            call fexit()
          end if
!
!       N-terminal residue request
        else if (resname(5:5).eq.'N') then
!         N-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.&
 &            (sickle.eq.51).OR.&
 &            ((sickle.ge.64).AND.(sickle.le.75)).OR.((sickle.ge.108).AND.(sickle.le.119)))) then
            foundN = .true.
!         unknown residues (have to always be explicitly N/C labeled)
          else if (foundit.EQV..false.) then
            foundN = .true.
            i = i + 1
!         allow redundant N-terminus spec.
          else if ((sickle.eq.27).OR.(sickle.eq.28).OR.((sickle.ge.76).AND.(sickle.le.87))) then
            foundN = .true.          
          end if
!
          if (foundN.EQV..true.) then
            if (startmol.EQV..true.) then
              write(ilog,*) 'Fatal. Found start of new molecule in s&
 &equence file without end of previous molecule (',resname(1:5),').'
              call fexit()
            end if
            nmol = nmol + 1
            startmol = .true.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible request for N-term. in seq. file (',resname(1:5),').'
            call fexit()
          end if
!
!       C-terminal residue request
        else if (resname(5:5).eq.'C') then
!         C-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.&
 &            (sickle.eq.51).OR.&
 &            ((sickle.ge.64).AND.(sickle.le.75)).OR.((sickle.ge.108).AND.(sickle.le.119)))) then
            foundC = .true.
!         unknown residues (have to always be explicitly N/C labeled)
          else if (foundit.EQV..false.) then
            foundC = .true.
            i = i + 1
          else if ((sickle.eq.29).OR.(sickle.eq.30)) then ! allow redundant C-terminus spec.
            foundC = .true.
          end if
!
          if (foundC.EQV..true.) then
            if (startmol.EQV..false.) then
              write(ilog,*) 'Fatal. Found end of current molecule in&
 & sequence file without ever starting it (',resname(1:5),').'
              call fexit()
            end if
            startmol = .false.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible request for C-term. in seq. file (',resname(1:5),').'
            call fexit()
          end if
!
        end if
!
!    the 7-letter codes are really 3-letter codes plus 2 appendices (_D, _C, _N)
!
     else if (length.eq.7) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            sickle = j
            foundit = .true.
          end if
        end do
!
!       parse for appendix
        if ((resname(4:4).ne.'_').OR.(resname(6:6).ne.'_').OR.&
 &          ((resname(5:5).ne.'D').AND.&
 &           (resname(5:5).ne.'N').AND.(resname(5:5).ne.'C')).OR.&
 &          ((resname(7:7).ne.'D').AND.&
 &           (resname(7:7).ne.'N').AND.(resname(7:7).ne.'C'))) then
          write(ilog,*) 'Fatal. Found an unknown residue in seq. fil&
 &e, which is currently not supported (',resname(1:7),').'
          call fexit()
        end if
!
!       chirality request and N-term. request
        if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'N')).OR.&
 &          ((resname(5:5).eq.'N').AND.(resname(7:7).eq.'D'))) then
!         D-aa in-chain residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.2).AND.(sickle.le.23)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.119)).OR.&
 &            (sickle.eq.25).OR.(sickle.eq.51))) then
            foundN = .true.              
          end if
!
          if (foundN.EQV..true.) then
            if (startmol.EQV..true.) then
              write(ilog,*) 'Fatal. Found start of new molecule in s&
 &equence file without end of previous molecule (',resname(1:7),').'
              call fexit()
            end if
            nmol = nmol + 1
            startmol = .true.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible N/D-request in seq. file (',resname(1:7),').'
            call fexit()
          end if
!
!       chirality request and C-term. request
        else if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'D'))) then
!         C-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.23)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.119)).OR.&
 &            (sickle.eq.51).OR.(sickle.eq.25))) then
            foundC = .true.              
          end if
!
          if (foundC.EQV..true.) then
            if (startmol.EQV..false.) then
              write(ilog,*) 'Fatal. Found end of current molecule in&
 & sequence file without ever starting it (',resname(1:7),').'
              call fexit()
            end if
            startmol = .false.
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible C/D-request in seq. file (',resname(1:7),').'
            call fexit()
          end if
!
!       N/C-terminal (free amino acid) request
        else if (((resname(5:5).eq.'N').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'N'))) then
!         N/C-terminal standard residues
          if ((foundit.EQV..true.).AND.(((sickle.ge.1).AND.(sickle.le.25)).OR.&
 &            ((sickle.ge.31).AND.(sickle.le.35)).OR.((sickle.ge.108).AND.(sickle.le.119)).OR.&
 &            (sickle.eq.51).OR.&
 &            ((sickle.ge.64).AND.(sickle.le.75)))) then
            foundC = .true.
            foundN = .true.
!         unknown residues (have to always be explicitly N/C labeled)
          else if (foundit.EQV..false.) then
            foundN = .true.
            foundC = .true.
            i = i + 1
          end if
!
          if ((foundC.EQV..true.).AND.(foundN.EQV..true.)) then
            if (startmol.EQV..true.) then
              write(ilog,*) 'Fatal. Found start of new molecule in s&
 &equence file without end of previous molecule (',resname(1:7),').'
              call fexit()
            end if
            nmol = nmol + 1
          else
            write(ilog,*) 'Fatal. Found an unknown res. or incompreh&
 &ensible request for N/C-term. in seq. file (',resname(1:7),').'
            call fexit()
          end if
!
        end if
!
!    the 9-letter codes are really 3-letter codes plus all 3 appendices (_D, _C, _N)
!
!     else if (length.eq.9) then
!        do j = 1,maxamino
!          if(resname(1:3) .eq. amino(j)) then
!            seqtyp(i) = j
!            chiral(i) = 1
!          end if
!        end do
!
     end if
!
     if (sickle.eq.-1) then
       if (resname(1:3).eq.'HIS') then
         write(ilog,*) 'Fatal. Residue name HIS is not allowed in sequence input (must be explicit about HIP/HID/HIE). If &
 &this is meant to be an unsupported residue, please rename it in both sequence and structure input.'
         call fexit()
       end if
       if ((laststartmol.EQV..false.).AND.(((resname(4:5).ne.'_N').AND.(length.eq.5)).OR.&
 &                                     ((resname(4:7).ne.'_N_C').AND.(length.eq.7)).OR.((length.ne.5).AND.(length.ne.7)))) then
         write(ilog,*) 'Fatal. An unsupported residue starting a new molecule must be a three-letter code &
 &followed by "_N" or "_N_C" (instead got: ',resname(p1:p2),').'
         call fexit()
       else if ((laststartmol.EQV..true.).AND.(((resname(4:5).ne.'_C').AND.(length.eq.5)).OR.&
 &                                         ((length.ne.3).AND.(length.ne.5)))) then
         write(ilog,*) 'Fatal. An unsupported residue continuing an existing molecule must be a three-letter code &
 &possibly followed by "_C" (instead got: ',resname(p1:p2),').'
         call fexit()
       end if
       n_pdbunk = n_pdbunk + 1
     end if
!
  end do
!
  if (nmol.eq.0) then
    write(ilog,*) 'Fatal error. Requested calculation without any mo&
 &lecules. Sane, but empty sequence file?'
    call fexit()
  end if
!
  if (done.EQV..false.) then
    write(ilog,*) 'Fatal error. Please terminate sequence file with &
 &END line.'
    call fexit()
  end if
!
! we now have the full residue and molecule information
  nseq = i
!
!  write(*,*) 'allocating pdbunk'
  if (n_pdbunk.gt.0) allocate(pdb_unknowns(n_pdbunk))
!  write(*,*) 'allocating seq'
  call allocate_sequen(aone)
!  write(*,*) 'allocating mol.'
  call allocate_molecule(aone)
!  write(*,*) 'allocating fyoc'
  call allocate_fyoc(aone)
!  write(*,*) 'allocating done'
!
  rewind(unit=iseqfile)
!  
! get the peptide sequence
  i = 0
  n_pdbunk = 0
  nmol = 0
  done = .false.
  startmol = .false.
  n_crosslinks = 0
  seqtyp(:) = -1
  chiral(:) = 1
  do while (done.EQV..false.)
!
     if (i.gt.0) then
       jj = 0
       read(seqlin(next:MAXSTRLEN),49,iostat=iomess) jj
       if ((iomess.eq.0).AND.(jj.gt.0).AND.(jj.ne.i).AND.(jj.le.nseq)) then
         if (i.lt.jj) then
           n_crosslinks = n_crosslinks + 1
           crosslink(n_crosslinks)%rsnrs(1) = i
           crosslink(n_crosslinks)%rsnrs(2) = jj
         end if
       end if
       do k=1,MAXSTRLEN
         seqlin(k:k) = ' '
       end do
     end if

     foundN = .false.
     foundC = .false.
     foundM = .false.
     foundit = .false.
!
     read(iseqfile,50,iostat=iomess) seqlin
     if (iomess.eq.IOSTAT_END) then
       close(unit=iseqfile)
       exit
     else if (iomess.eq.2) then
       write(ilog,*) 'Fatal. File I/O error while processing sequenc&
 &e input file (got: ',seqfile(t1:t2),').'
       call fexit()
     end if
     call toupper(seqlin)
     next = 1
     call extract_str(seqlin,resname,next)
     call strlims(resname,p1,p2)
     length = p2-p1+1
!
!    first let's check for N-terminal residues
     if ((length.eq.3).OR.(length.eq.2)) then
        if (resname.eq.'END') then
          done = .true.
          exit
        end if
        if (length.eq.2) then
          resname(3:3) = ' '
        end if
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            seqtyp(i) = j
            chiral(i) = 1
            foundit = .true.
          end if
        end do
        if (foundit.EQV..true.) then
!         N-terminal caps
          if ((seqtyp(i).eq.27).OR.(seqtyp(i).eq.28).OR.&
 &                 ((seqtyp(i).ge.76).AND.(seqtyp(i).le.87))) then
            foundN = .true.
!         C-terminal caps
          else if ((seqtyp(i).eq.29).OR.(seqtyp(i).eq.30)) then
            foundC = .true.
!         standard in-chain residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &                 ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &                 (seqtyp(i).eq.51).OR.&
 &                 ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &                 ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
            foundM = .true.
!         small molecules and atoms
          else if (((seqtyp(i).ge.36).AND.(seqtyp(i).le.50)).OR.&
 &                 ((seqtyp(i).ge.52).AND.(seqtyp(i).le.63)).OR.&
 &                 ((seqtyp(i).ge.88).AND.(seqtyp(i).le.107))) then
            foundN = .true.  
            foundC = .true.
!          else
!            write(ilog,*) 'Fatal. Unable to parse residue class in makepept(). This is a bug.'
!            call fexit()
          end if
!       unknown residues (have to always be explicitly N/C labeled)
        else ! if (pdb_analyze.EQV..true.) then
          foundM = .true.
          i = i + 1
        end if
!
!       now see what we've got
        if ((foundN.EQV..true.).AND.(foundC.EQV..true.)) then
          nmol = nmol + 1
          rsmol(nmol,1) = i
          rsmol(nmol,2) = i
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
        else if (foundN.EQV..true.) then
          nmol = nmol + 1
          rsmol(nmol,1) = i
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
          startmol = .true.
        else if (foundC.EQV..true.) then
          rsmol(nmol,2) = i
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
          startmol = .false.
        else if (foundM.EQV..true.) then
          phi(i) = 0.0d0
          psi(i) = 0.0d0
          omega(i) = 0.0d0
          do j = 1, 4
            chi(j,i) = 0.0d0
          end do
        end if
!
!    the 5-letter codes are really 3-letter codes plus an appendix (_D, _C, _N)
     else if (length.eq.5) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            seqtyp(i) = j
            chiral(i) = 1
            foundit = .true.
          end if
        end do
!
!       chirality request
        if (resname(5:5).eq.'D') then
!         D-aa in-chain residues
          if (foundit.EQV..true.) then
            if (((seqtyp(i).ge.2).AND.(seqtyp(i).le.23)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.25).OR.(seqtyp(i).eq.51).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
              foundM = .true.              
            end if
          end if
!
          if (foundM.EQV..true.) then
            chiral(i) = -1
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
          end if
!
!       N-terminal residue request
        else if (resname(5:5).eq.'N') then
!         unknown residues (have to always be explicitly N/C labeled)
          if (foundit.EQV..false.) then
            foundN = .true.
            i = i + 1
!         N-terminal standard residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.51).OR.&
 &              ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
            foundN = .true.
!         allow redundant N-terminus spec.
          else if ((seqtyp(i).eq.27).OR.(seqtyp(i).eq.28).OR.((seqtyp(i).ge.76).AND.(seqtyp(i).le.87))) then
            foundN = .true.
          end if
!
          if (foundN.EQV..true.) then
            nmol = nmol + 1
            rsmol(nmol,1) = i
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,1) = 1
            startmol = .true.
          end if
!
!       C-terminal residue request
        else if (resname(5:5).eq.'C') then
!         unknown residues (have to always be explicitly N/C labeled)
          if (foundit.EQV..false.) then
            foundC = .true.
            i = i + 1
!         C-terminal standard residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &             ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &             (seqtyp(i).eq.51).OR.&
 &             ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &             ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
            foundC = .true.
          else if ((seqtyp(i).eq.29).OR.(seqtyp(i).eq.30)) then ! allow redundant C-terminus spec.
            foundC = .true.
          end if
!
          if (foundC.EQV..true.) then
            rsmol(nmol,2) = i
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,2) = 1
            startmol = .false.
          end if
!
        end if
!
!    the 7-letter codes are really 3-letter codes plus 2 appendices (_D, _C, _N)
!
     else if (length.eq.7) then
        do j = 1,maxamino
          if(resname(1:3) .eq. amino(j)) then
            i = i + 1
            seqtyp(i) = j
            chiral(i) = 1
            foundit = .true.
          end if
        end do
!
!       chirality request and N-term. request
        if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'N')).OR.&
 &          ((resname(5:5).eq.'N').AND.(resname(7:7).eq.'D'))) then
!         D-aa in-chain residues
          if (foundit.EQV..true.) then
            if (((seqtyp(i).ge.2).AND.(seqtyp(i).le.23)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.25).OR.(seqtyp(i).eq.51).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
              foundN = .true.              
            end if
          end if
!
          if (foundN.EQV..true.) then
            nmol = nmol + 1
            rsmol(nmol,1) = i
            chiral(i) = -1
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,1) = 1
            startmol = .true.
          end if
!
!       chirality request and C-term. request
        else if (((resname(5:5).eq.'D').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'D'))) then
!         C-terminal standard residues
          if (foundit.EQV..true.) then
            if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.23)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.51).OR.(seqtyp(i).eq.25).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
              foundC = .true.              
            end if
          end if
!
          if (foundC.EQV..true.) then
            rsmol(nmol,2) = i
            chiral(i) = -1
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,2) = 1
            startmol = .false.
          end if
!
!       N/C-terminal (free amino acid) request
        else if (((resname(5:5).eq.'N').AND.(resname(7:7).eq.'C'))&
 &       .OR.((resname(5:5).eq.'C').AND.(resname(7:7).eq.'N'))) then
!         unknown residues (have to always be explicitly N/C labeled)
          if (foundit.EQV..false.) then
            foundN = .true.
            foundC = .true.
            i = i + 1 
!         N/C-terminal standard residues
          else if (((seqtyp(i).ge.1).AND.(seqtyp(i).le.25)).OR.&
 &              ((seqtyp(i).ge.31).AND.(seqtyp(i).le.35)).OR.&
 &              (seqtyp(i).eq.51).OR.&
 &              ((seqtyp(i).ge.64).AND.(seqtyp(i).le.75)).OR.&
 &              ((seqtyp(i).ge.108).AND.(seqtyp(i).le.119))) then
            foundC = .true.
            foundN = .true.
          end if
!
          if ((foundC.EQV..true.).AND.(foundN.EQV..true.)) then
            nmol = nmol + 1
            rsmol(nmol,1) = i
            rsmol(nmol,2) = i
            phi(i) = 0.0d0
            psi(i) = 0.0d0
            omega(i) = 0.0d0
            do j = 1, 4
              chi(j,i) = 0.0d0
            end do
            moltermid(nmol,1) = 1
            moltermid(nmol,2) = 1
          end if
!
        end if
!
!    the 9-letter codes are really 3-letter codes plus all 3 appendices (_D, _C, _N)
!
!     else if (length.eq.9) then
!        do j = 1,maxamino
!          if(resname(1:3) .eq. amino(j)) then
!            seqtyp(i) = j
!            chiral(i) = 1
!          end if
!        end do
!
     end if
!
     if (seqtyp(i).eq.-1) then
       n_pdbunk = n_pdbunk + 1
       pdb_unknowns(n_pdbunk)(1:3) = resname(1:3)
     end if
!
  end do
!
  if (done.EQV..false.) then
    write(ilog,*) 'Fatal error. Please terminate sequence file with &
 &END line.'
    call fexit()
  end if
!
  nseq = i
!
! set molofrs
  do imol=1,nmol
    do i=rsmol(imol,1),rsmol(imol,2)
      molofrs(i) = imol
    end do
  end do
!
! for any unsupported residues, grab as much info as possible from PDB
  if (n_pdbunk.gt.0) then
    call infer_from_pdb(1)
  else if ((use_pdb_template.EQV..true.).AND.(pdb_analyze.EQV..false.)) then
    write(ilog,*) 'Warning. Ignoring superfluous request for use of pdb-template (no unsupported residues &
 &and not in trajectory analysis mode).'
    use_pdb_template = .false.
  end if
!
! now let's initialize some more arrays
  call setup_resrad()
!  write(*,*) 'allocating polypep.'
  call allocate_polypep(aone)
!  write(*,*) 'allocating atomestls.'
  call allocate_atomestls(aone)
!  write(*,*) 'allocating done'
  if (n_crosslinks.gt.0) then
    call setup_crosslinks()
  end if
!
  do imol=1,nmol
    do rs = rsmol(imol,1),rsmol(imol,2)
      if (((seqtyp(rs).ge.1).AND.(seqtyp(rs).le.25)).OR.&
 &        ((seqtyp(rs).ge.31).AND.(seqtyp(rs).le.35)).OR.&
 &        (seqtyp(rs).eq.51).OR.((seqtyp(rs).ge.108).AND.(seqtyp(rs).le.119))) then
        phi(rs) = 179.5d0!-57.0
        psi(rs) = 179.5d0!-47.0
        omega(rs) = 180.0d0
        if (ua_model.gt.0) then
          if (rs.eq.rsmol(imol,1)+1) then
            if (seqtyp(rs-1).eq.28) then
              omega(rs) = 0.0
            end if
          end if
        end if
        do jj = 1,MAXCHI
          chi(jj,rs) = 180.0d0
!         override for P
          if (seqtyp(rs).eq.9) chi(jj,rs) = 0.0d0
        end do
!       chi3-override for Q/E
        if (seqtyp(rs).eq.19) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.18) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.108) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.113) chi(3,rs) = 90.0d0
        if (seqtyp(rs).eq.115) chi(3,rs) = 120.0d0
!       the hack of initializing pucker states randomly is made defunct (rr is always <= 1.0)
!       the chi-hijack still works to communicate the choice of pucker state to
!       sidechain(...)
!       do not change unless you understand the consequences
        if(seqtyp(rs).eq.9) then
          if (rs.gt.rsmol(imol,1)) then
            rr = 1.0 ! random()
            if (rr.le.1.5) then
! exo (C-gamma far from Ci; phi +/-60.8)
              if (chiral(rs).eq.1) then
                chi(1,rs) = -28.0d0
                chi(2,rs) =  4.3d0  ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -60.8d0
              else
                chi(1,rs) = 30.5d0  ! 28.0d0
                chi(2,rs) = -10.0d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 60.8d0
              end if
            else
! endo (C-gamma closer to Ci; phi -/+71.0)
              if (chiral(rs).eq.1) then
                chi(1,rs) = 31.0d0   ! 32.0d0
                chi(2,rs) =  -11.9d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -71.0d0
              else
                chi(1,rs) = -29.0d0  ! -31.0d0
                chi(2,rs) = 6.4d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 71.0d0
              end if
            end if
          else
            rr = 1.0 ! random()
            if (rr.le.1.5) then
              if (chiral(rs).eq.1) then
                chi(1,rs) = -28.0d0
                chi(2,rs) =  4.3d0  ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -0.8d0
              else
                chi(1,rs) = 30.5d0  ! 28.0d0
                chi(2,rs) = -10.0d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 120.8d0
              end if
            else
              if (chiral(rs).eq.1) then
                chi(1,rs) = 31.0d0   ! 32.0d0
                chi(2,rs) =  -11.9d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = -11.0d0
              else
                chi(1,rs) = -29.0d0  ! -31.0d0
                chi(2,rs) = 6.4d0 ! -1.3503d0*chi(1,rs) + 3.2392d0
                phi(rs) = 131.0d0
              end if
            end if
          end if
!       dealing with hydroxyproline
        else if (seqtyp(rs).eq.32) then
          phi(rs) = -60.0d0
          psi(rs) = 110.0d0
          chi(1,rs) = -29.4d0
          chi(2,rs) =  39.6d0 
        end if
!       some small molecules have residual degrees of freedom
      else if (((seqtyp(rs).ge.41).AND.(seqtyp(rs).le.44)).OR.&
 &             (seqtyp(rs).eq.48).OR.(seqtyp(rs).eq.50).OR.&
 &             ((seqtyp(rs).ge.88).AND.(seqtyp(rs).le.99))) then
        do jj = 1,MAXCHI
          chi(jj,rs) = 180.0d0
        end do
!       override for PPA
        if (seqtyp(rs).eq.44) then
          if (random().gt.0.5) then
            chi(1,rs) = 120.0d0
          else
            chi(1,rs) = -120.0d0
          end if
        end if
      else if ((seqtyp(rs).ge.29).AND.(seqtyp(rs).le.30)) then
        omega(rs) = 180.0d0
      else if ((seqtyp(rs).ge.41).AND.(seqtyp(rs).le.42)) then
        omega(rs) = 0.0d0
!      nucleic acids have all kinds of degrees of freedom
      else if ((seqtyp(rs).ge.64).AND.(seqtyp(rs).le.87)) then
        do jj = 1,6
          nucs(jj,rs) = 180.0d0
        end do
        do jj = 1,MAXCHI
          chi(jj,rs) = 180.0d0
        end do
!       sugar pucker missing
      end if
    end do
  end do
!     
! set up polypeptide chain geometry
  call proteus_init(aone)
  do imol=1,nmol
   call proteus(imol)
  end do
  call proteus_init(atwo)
!
  if (n_pdbunk.gt.0) call infer_from_pdb(2)
!
! top level patches first -> biotypes (these potentially reset everything beyond topology)
  call read_biotpatchfile()
!
! potential mass patches a.s.a.p.
  call read_atpatchfile(athree)
!
! from the fully assembled Z-matrix, build backbone coordinates
! and populate 1-2-topology arrays
  call makexyz2()
  call setup_connectivity_2()
!
! now regenerate internals from the newly found Cartesian coordinates
! this serves to remove intrinsic incompatibilities in the Z-matrix
! (like ill-defined bond angles) which were adjusted by genxyz(...) called
! from within makexyz2()
! also set up 1-3 and 1-4 topology arrays
  do imol=1,nmol
    call genzmat(imol)
  end do
  call makexyz2()
  call setup_connectivity_34()
!
! potential changes to LJ parameters must occur early due to reliance of other settings on them
  call read_ljpatchfile()
!
! now allocate memory for atomic parameters and assign their values
!  write(*,*) 'allocating atomprms.'
  call allocate_atomprms(aone)
!  write(*,*) 'done'
!
! check atomic valencies
  call valence_check()
!
! some high-level residue-wise parameters
  call absinth_residue()
!
! populate short-range interaction  and more atomic arrays: note that
! setup_srinter has to be called before setup_savol (within absinth_atom)
! can be called!
!  write(*,*) 'allocating inter.'
  call allocate_inter(aone)
!  write(*,*) 'done'
  call setup_srinter()
!
  call absinth_atom()
  call makexyz2()
!
! assign residue- and molecule-wise parameters
  call parse_sequence()
  call absinth_molecule()
  if (n_crosslinks.gt.0) then
    call setup_crosslinks2()
  end if
!
  do imol=1,nmol
    call find_rotlsts(imol)
    call parse_rotlsts(imol)
  end do
  call trans_rotlstdof()
  call find_ntorsn()
! forces requires find_ntorsn() to be done
!  write(*,*) 'allocating forces.'
  call allocate_forces(aone)
!  write(*,*) 'done'
  call correct_srinter()
  call correct_sampler()
!
! populate molecular RBC-arrays
  do imol=1,nmol
    call update_rigid(imol)
    call update_rigidm(imol)
  end do
!
! set some analysis flags and check water loops (note that BT/LJ patches have already been applied)
  call parse_moltyp()
!
! populate helper arrays
  do rs=1,nseq
    rsinfo(rs,1) = at(rs)%bb(1)
    rsinfo(rs,2) = at(rs)%na - 1
  end do
  do imol=1,nmol
    molinfo(imol,1) = atmol(imol,2)-atmol(imol,1)+1
  end do
!
! make a back up of Z matrix
  bangpr(1:n) = bang(1:n)
  blenpr(1:n) = blen(1:n)
!
end
!
!----------------------------------------------------------------------
!
subroutine setup_crosslinks()
!
  use sequen
  use aminos
  use iounit
  use fyoc
  use molecule
!
  implicit none
!
  integer i,j,rs1,rs2,imol
  logical samemol
!
  do i=1,n_crosslinks
    do j=i+1,n_crosslinks
      if ((crosslink(j)%rsnrs(1).eq.crosslink(i)%rsnrs(1)).OR.&
 &        (crosslink(j)%rsnrs(1).eq.crosslink(i)%rsnrs(2)).OR.&
 &        (crosslink(j)%rsnrs(2).eq.crosslink(i)%rsnrs(1)).OR.&
 &        (crosslink(j)%rsnrs(2).eq.crosslink(i)%rsnrs(2))) then
        write(ilog,*) 'Fatal. Only one type of chemical cross-link&
 & is possible for a given residue. Check sequence input.'
        call fexit()
      end if
    end do
    rs1 = crosslink(i)%rsnrs(1)
    rs2 = crosslink(i)%rsnrs(2)
    if ((seqtyp(rs1).eq.8).AND.(seqtyp(rs2).eq.8)) then
      samemol = .true.
      do imol=1,nmol
        if ((rsmol(imol,1).le.rs1).AND.(rsmol(imol,2).ge.rs1)) then
          if ((rsmol(imol,1).gt.rs2).OR.(rsmol(imol,2).lt.rs2)) then
            samemol = .false.
            exit
          end if
        end if
      end do
      if (samemol.EQV..true.) then
        crosslink(i)%itstype = 1 ! intramol. disulfide
      else
        crosslink(i)%itstype = 2 ! intermol. disulfide
      end if
      disulf(rs1) = rs2
      disulf(rs2) = rs1
    else
      write(ilog,*) 'Fatal. Chemical cross-links are currently not supported&
 & between residues ',amino(seqtyp(rs1)),' and ',amino(seqtyp(rs2)),'. Check back later.'
      call fexit()
    end if
  end do
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine setup_crosslinks2()
!
  use sequen
  use iounit
  use molecule
!
  implicit none
!
  integer i,j,k,l,kk,ends,cnt,nlks,imol,jmol,last,first,diffo
  integer, ALLOCATABLE:: lklst(:,:), counts(:), idx(:), counts2(:)
  logical, ALLOCATABLE:: donewith(:)
  logical notdone,override
!
  allocate(lklst(nseq,2))
  allocate(counts(nmol))
  allocate(counts2(nmol))
  allocate(donewith(nseq))
  allocate(idx(n_crosslinks+3))
!
! try to arrange intermolecular crosslinks such that they can be observed sequentially (inside-out)
  donewith(:) = .false.
  counts(:) = 0
!
  do imol=1,nmol
    do j=1,n_crosslinks
      if (donewith(j).EQV..true.) cycle
      if (molofrs(crosslink(j)%rsnrs(1)).ne.molofrs(crosslink(j)%rsnrs(2))) then
        if ((molofrs(crosslink(j)%rsnrs(1)).eq.imol).OR.(molofrs(crosslink(j)%rsnrs(2)).eq.imol)) then
          counts(imol) = counts(imol) + 1
        end if
      end if
    end do
  end do
!
  idx(:) = 0
  first = 1
  last = n_crosslinks
  override = .false.
  do while (first.le.last)
!
    counts2(:) = counts(:)
    diffo = last-first
    do j=1,n_crosslinks
      imol = molofrs(crosslink(j)%rsnrs(1))
      jmol = molofrs(crosslink(j)%rsnrs(2))
      if (donewith(j).EQV..true.) cycle
      if ((imol.ne.jmol).AND.((counts(imol).eq.1).OR.(counts(jmol).eq.1))) then
        idx(j) = last
        last = last - 1
        counts2(imol) = counts(imol) - 1 
        counts2(jmol) = counts(jmol) - 1
        if (counts(jmol).gt.1) then
          kk = crosslink(j)%rsnrs(1)
          crosslink(j)%rsnrs(1) = crosslink(j)%rsnrs(2)
          crosslink(j)%rsnrs(2) = kk
        end if
        donewith(j) = .true.
      else if (imol.eq.jmol) then
        idx(j) = first
        first = first + 1
        donewith(j) = .true.
      else
        if (override.EQV..true.) then
!         this means that are at least 2 intermolecular links per linked mol. (cycles), which will prevent randomization later
!         (as in: may crash)
          idx(j) = first
          first = first + 1
          donewith(j) = .true.
        end if
      end if
    end do
    if ((override.EQV..true.).AND.((last-first).eq.diffo)) then
      write(ilog,*) 'Interminable loop in setup_crosslink2(...). This is a bug.'
      call fexit()
    else if (override.EQV..true.) then
      override = .false.
    else if ((last-first).eq.diffo) then
      override = .true.
    end if
 !
    counts(:) = counts2(:)
  end do
!
  do j=1,n_crosslinks
    if (j.eq.idx(j)) cycle
    do i=1,n_crosslinks
      if (idx(i).eq.j) exit
    end do
    do k=1,2
      kk = crosslink(j)%rsnrs(k)
      crosslink(j)%rsnrs(k) = crosslink(i)%rsnrs(k)
      crosslink(i)%rsnrs(k) = kk
    end do
    kk = crosslink(j)%itstype
    crosslink(j)%itstype = crosslink(i)%itstype
    crosslink(i)%itstype = kk
    kk = idx(j)
    idx(j) = idx(i)
    idx(i) = kk
  end do
!
! set coupling sets for crosslinks (currently unused)
  nlks = 0
  donewith(:) = .false.
  override = .false.
  do i=1,n_crosslinks
    if (donewith(i).EQV..true.) cycle
    if (molofrs(crosslink(i)%rsnrs(2)).ne.molofrs(crosslink(i)%rsnrs(1))) then
      nlks = 1
      lklst(1,1) = i
      notdone = .true.
      imol = molofrs(crosslink(i)%rsnrs(1))
      jmol = molofrs(crosslink(i)%rsnrs(2))
      counts(1) = imol
      counts(2) = jmol
      donewith(i) = .true.
      cnt = 2
      ends = 2
      do while (notdone.EQV..true.)
        notdone = .false.
        do j=1,n_crosslinks
          if (donewith(j).EQV..true.) cycle
          if (molofrs(crosslink(j)%rsnrs(2)).ne.molofrs(crosslink(j)%rsnrs(1))) then
            do k=1,cnt
              if (molofrs(crosslink(j)%rsnrs(2)).eq.counts(k)) then
                donewith(j) = .true.
                notdone = .true.
                kk = 0
                do l=1,cnt
                  if (l.eq.k) cycle
                  if (counts(l).eq.molofrs(crosslink(j)%rsnrs(1))) then
                    ends = ends - 2
                    kk = 1
                    if (override.EQV..false.) then
                      write(ilog,*) 'Warning. Intermolecular crosslinks that produce&
 & ring topologies may cause CAMPARI to crash later.'
                      override = .true.
                    end if
                    exit
                  end if
                end do
                if (kk.eq.0) then
                  cnt = cnt + 1
                  counts(cnt) = molofrs(crosslink(j)%rsnrs(1)) 
                  nlks = nlks + 1
                  lklst(nlks,1) = j
                  exit
                end if
              else if (molofrs(crosslink(j)%rsnrs(1)).eq.counts(k)) then
                donewith(j) = .true.
                notdone = .true.
                kk = 0
                do l=1,cnt
                  if (l.eq.k) cycle
                  if (counts(l).eq.molofrs(crosslink(j)%rsnrs(2))) then
                    ends = ends - 2
                    kk = 1
                    if (override.EQV..false.) then
                      write(ilog,*) 'Warning. Intermolecular crosslinks that produce&
 & ring topologies may cause CAMPARI to crash later.'
                      override = .true.
                    end if
                    exit
                  end if
                end do
                if (kk.eq.0) then
                  cnt = cnt + 1
                  counts(cnt) = molofrs(crosslink(j)%rsnrs(2))
                  nlks = nlks + 1
                  lklst(nlks,1) = j
                  exit
                end if
              end if
            end do
          end if
        end do
        if (ends.le.0) notdone = .false.
      end do
      do k=1,nlks
        crosslink(lklst(k,1))%nolks = nlks-1
        kk = 0
        do l=1,nlks
          if (l.eq.k) cycle
          kk = kk + 1
          crosslink(lklst(k,1))%olks(kk) = lklst(l,1)
        end do
!        write(ilog,*) 'Found ',crosslink(k)%olks(1:crosslink(k)%nolks),' for ',lklst(k,1)
      end do
    end if
  end do
!
  deallocate(lklst)
  deallocate(donewith)
  deallocate(counts)
  deallocate(counts2)
  deallocate(idx)
!
  do i=1,n_crosslinks
    crlk_idx(crosslink(i)%rsnrs(1)) = i
    crlk_idx(crosslink(i)%rsnrs(2)) = i
  end do
  call crosslink_excludes()
!
end
!
!--------------------------------------------------------------------------
!
! a subroutine to parse the sequence (molecule-wise) to set some analysis
! flags and to allow pooling of the data for identical molecules in some
! cases
!
subroutine parse_sequence()
!
  use iounit
  use sequen
  use molecule
  use polypep
  use pdb
  use atoms
  use fyoc
!
  implicit none
!
  integer imol,rs,jmol,k,kkk,aone,cl
  integer, ALLOCATABLE:: dummy(:)
  logical newmt
!
  allocate(dummy(nmol))
!
! first a dry run to establish nmoltyp
  nmoltyp = 1
  dummy(1) = 1
!
  do imol=2,nmol
    newmt = .true.
    do jmol=1,nmoltyp
      if ((rsmol(dummy(jmol),2)-rsmol(dummy(jmol),1)).eq.&
 &        (rsmol(imol,2)-rsmol(imol,1))) then
        newmt = .false.
        do rs=1,rsmol(imol,2)-rsmol(imol,1)+1
          if (seqtyp(rs+rsmol(imol,1)-1).ne.&
 &            seqtyp(rs+rsmol(dummy(jmol),1)-1)) then
            newmt = .true.
            exit
          end if
          if (chiral(rs+rsmol(imol,1)-1).ne.&
 &            chiral(rs+rsmol(dummy(jmol),1)-1)) then
            newmt = .true.
            exit
          end if
          if (seqtyp(rs+rsmol(imol,1)-1).eq.26) then ! always declare new for unknowns unless biotype seq is identical
            cl = at(rs+rsmol(dummy(jmol),1)-1)%bb(1)
            if ((at(rs+rsmol(dummy(jmol),1)-1)%nbb.ne.at(rs+rsmol(imol,1)-1)%nbb).OR.&
 &              (at(rs+rsmol(dummy(jmol),1)-1)%nsc.ne.at(rs+rsmol(imol,1)-1)%nsc)) then
              newmt = .true.
              exit
            end if
            do k=at(rs+rsmol(imol,1)-1)%bb(1),at(rs+rsmol(imol,1)-1)%bb(1)+at(rs+rsmol(imol,1)-1)%nbb+at(rs+rsmol(imol,1)-1)%nsc-1
              if (b_type(k).ne.b_type(cl)) then
                newmt = .true.
                exit
              else
                cl = cl + 1
              end if
            end do
            if (newmt.EQV..true.) exit
          end if
        end do
        if (newmt.EQV..false.) then
!         always declare a new molecule type for crosslinked molecules (this may be overkill but so what)
          do cl=1,n_crosslinks
            if (((crosslink(cl)%rsnrs(1).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(1).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(dummy(jmol),2)))) then
              newmt = .true.
            end if
            if (((crosslink(cl)%rsnrs(2).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(2).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(dummy(jmol),2)))) then
              newmt = .true.
            end if
          end do
        end if
        if (newmt.EQV..false.) exit
      end if
    end do
    if (newmt.EQV..true.) then
      nmoltyp = nmoltyp + 1
      dummy(nmoltyp) = imol
    end if
  end do
!
! now allocate memory
  aone = 1
  call allocate_molecule_type(aone)
!
! the first molecule defines the first type and adds one count to that
  kkk = 1
  moltyp(1,1) = 1
  moltyp(1,2) = 1
  moltypid(1) = 1
!
  do imol=2,nmol
!
    newmt = .true.
!
    do jmol=1,kkk
      if ((rsmol(moltyp(jmol,1),2)-rsmol(moltyp(jmol,1),1)).eq.&
 &        (rsmol(imol,2)-rsmol(imol,1))) then
        newmt = .false.
        do rs=1,rsmol(imol,2)-rsmol(imol,1)+1
          if (seqtyp(rs+rsmol(imol,1)-1).ne.&
 &            seqtyp(rs+rsmol(moltyp(jmol,1),1)-1)) then
            newmt = .true.
            exit
          end if
          if (chiral(rs+rsmol(imol,1)-1).ne.&
 &            chiral(rs+rsmol(moltyp(jmol,1),1)-1)) then
            newmt = .true.
            exit
          end if
          if (seqtyp(rs+rsmol(imol,1)-1).eq.26) then ! always declare new for unknowns
            cl = at(rs+rsmol(moltyp(jmol,1),1)-1)%bb(1)
            if ((at(rs+rsmol(moltyp(jmol,1),1)-1)%nbb.ne.at(rs+rsmol(imol,1)-1)%nbb).OR.&
 &              (at(rs+rsmol(moltyp(jmol,1),1)-1)%nsc.ne.at(rs+rsmol(imol,1)-1)%nsc)) then
              newmt = .true.
              exit
            end if
            do k=at(rs+rsmol(imol,1)-1)%bb(1),at(rs+rsmol(imol,1)-1)%bb(1)+at(rs+rsmol(imol,1)-1)%nbb+at(rs+rsmol(imol,1)-1)%nsc-1
              if (b_type(k).ne.b_type(cl)) then
                newmt = .true.
                exit
              else
                cl = cl + 1
              end if
            end do
            if (newmt.EQV..true.) exit
          end if
        end do
      end if
!
      if (newmt.EQV..false.) then
!       always declare a new molecule type for crosslinked molecules (this may be overkill but so what)
        do cl=1,n_crosslinks
          if (((crosslink(cl)%rsnrs(1).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(1).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(1).le.rsmol(dummy(jmol),2)))) then
            newmt = .true.
          end if
          if (((crosslink(cl)%rsnrs(2).ge.rsmol(imol,1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(imol,2))).OR.&
 &((crosslink(cl)%rsnrs(2).ge.rsmol(dummy(jmol),1)).AND.(crosslink(cl)%rsnrs(2).le.rsmol(dummy(jmol),2)))) then
            newmt = .true.
          end if
        end do
      end if
!
      if (newmt.EQV..false.) then
        moltyp(jmol,2) = moltyp(jmol,2) + 1
        moltypid(imol) = jmol
        exit
      end if
    end do
! 
    if (newmt.EQV..true.) then
      kkk = kkk + 1
      moltyp(kkk,1) = imol
      moltyp(kkk,2) = 1
      moltypid(imol) = kkk
    end if
!
  end do
!
! set total number of rigid-body degrees of freedom, assign analysis groups (default), and set solvent/solute
  totrbd = 0
  nressolute = 0
  nsolutes = 0
  nangrps = nmoltyp
  do imol=1,nmol
    an_grp_mol(imol) = moltypid(imol)
    if ((rsmol(imol,2)-rsmol(imol,1)).eq.0) is_solvent(imol) = .true.
    if (atmol(imol,2)-atmol(imol,1).eq.0) then 
      totrbd = totrbd + 3
    else if (atmol(imol,2)-atmol(imol,1).eq.1) then
      totrbd = totrbd + 5
    else
      totrbd = totrbd + 6
    end if
  end do
!
! read (and eventually overwrite analysis group spec.s)
  call read_analysisgrpfile()
  call alloc_angrps(an_grp_mol,is_solvent)
  if (just_solutes.EQV..true.) then
    pdbeffn = 0
    do imol=1,nmol
      if (is_solvent(imol).EQV..false.) pdbeffn = pdbeffn + (atmol(imol,2)-atmol(imol,1)+1)
    end do
  else
    pdbeffn = n
  end if
!
  deallocate(dummy)
!
end
!
!---------------------------------------------------------------------------
!
! this routine sets analysis flags and checks for water loop support
! also provides sequence report (formerly in parse_sequence)
!
subroutine parse_moltyp()
!
  use iounit
  use sequen
  use molecule
  use atoms
  use polypep
  use zmatrix
  use cutoffs
  use params
  use energies
  use inter
  use fyoc
  use aminos
  use units
  use system
!
  implicit none
!
  integer jmol,ati,rs,k,j,rs2,ii,ixx
  RTYPE, ALLOCATABLE:: epsik(:)
  logical dropout,foundit
!
! fix resrad for polynucleotides with 5' phosphate
  do rs=1,nseq
    if ((seqflag(rs).eq.22).AND.(rs.eq.rsmol(molofrs(rs),1)).AND.(amino(seqtyp(rs)).ne.'R5P').AND.&
 &      (amino(seqtyp(rs)).ne.'D5P')) then
      resrad(rs) = resrad(rs) + 0.8
    end if
  end do
!
! go through the molecule types and analyze what can be done with them
  do jmol=1,nmoltyp
!
    do ati=atmol(moltyp(jmol,1),1)+1,atmol(moltyp(jmol,1),2)
      if (izrot(ati)%alsz.gt.0) then
        do_pol(jmol) = .true.
        exit
      end if
    end do
!
    if (rsmol(moltyp(jmol,1),2).gt.rsmol(moltyp(jmol,1),1)) then
!     do_pers requires two specific reference atoms per residues (caps are allowed to misbehave) and
!     at least two eligible residues
      do_pers(jmol) = .true.
      k = 0
      do rs=rsmol(moltyp(jmol,1),1),rsmol(moltyp(jmol,1),2)
        if ((ci(rs).gt.0).AND.(ni(rs).gt.0)) then
          k = k + 1
        else if ((nuci(rs,2).gt.0).AND.(nuci(rs,6).gt.0)) then
          k = k + 1
        else 
          if ((rs.gt.rsmol(moltyp(jmol,1),1)).AND.(rs.lt.rsmol(moltyp(jmol,1),2))) then
            do_pers(jmol) = .false.
            exit
          end if
        end if
      end do
      if (k.lt.2) do_pers(jmol) = .false.
!     terminal residues are never suited for traditional polypeptide torsional
!     analysis as they have at most one informative backbone angle
!     note that this has to change should one do generic torsional analysis
      do rs=rsmol(moltyp(jmol,1),1)+1,rsmol(moltyp(jmol,1),2)-1
        if (seqpolty(rs).eq.'P') then
          do_tors(jmol) = .true.
          exit
        end if
      end do
    end if
  end do
!
! finally, find if there is a way to speed up a prototypical calculation
! through dedicated loops
! the assumed architecture of the sequence file is:
! arbitrary residues / molecules first, terminated by water solvent of
! single type
  foundit = .false.
  use_waterloops = .false.
  is_3site = .false.
  is_4site = .false.
  is_5site = .false.
  rsw1 = nseq + 1
  k = 0
  do rs=1,nseq
    if (((seqtyp(rs).eq.39).OR.(seqtyp(rs).eq.40).OR.(seqtyp(rs).eq.102).OR.&
 &       (seqtyp(rs).eq.45).OR.(seqtyp(rs).eq.103)).AND.(foundit.EQV..false.)) then
      foundit = .true.
      k = seqtyp(rs)
      rsw1 = rs
      do rs2=rs+1,nseq
        if (seqtyp(rs2).ne.k) then
          rsw1 = nseq + 1
          foundit = .false.
          if ((seqtyp(rs2).ne.39).AND.(seqtyp(rs2).ne.40).AND.(seqtyp(rs2).ne.45).AND.&
 &            (seqtyp(rs2).ne.103).AND.(seqtyp(rs2).ne.102)) then
            write(ilog,*) 'Warning. This calculation might be able to run &
 &faster if water molecules are placed at the end of sequence input.'
          else
            write(ilog,*) 'Warning. Concurrent use of different water models &
 &means that specialized water loops are no longer usable.'
          end if
          write(ilog,*)
          exit
        end if
      end do
      if (rsw1.eq.nseq) then
        foundit = .false.
      end if
      if (foundit.EQV..false.) then
        use_waterloops = .false.
        exit
      end if
    end if
    if (foundit.EQV..true.) then
      use_waterloops = .true.
      if ((seqtyp(rsw1).eq.39).OR.(seqtyp(rsw1).eq.40)) is_3site = .true.
      if ((seqtyp(rsw1).eq.45).OR.(seqtyp(rsw1).eq.103)) is_4site = .true.
      if (seqtyp(rsw1).eq.102) is_5site = .true.
      exit
    end if
  end do
! check for general Hamiltonian settings
  if (use_waterloops.EQV..true.) then
    if ((use_IPP.EQV..true.).AND.(use_attLJ.EQV..true.).AND.&
 &      (use_WCA.EQV..false.).AND.(use_TABUL.EQV..false.).AND.&
 &      (use_IMPSOLV.EQV..false.).AND.(use_POLAR.EQV..true.)) then
!     do nothing
    else
      write(ilog,*) 'Warning. Specialized water loops are disabled since energy function contains &
 &nonstandard terms in nonbonded interactions or is missing one or more standard terms.'
      use_waterloops = .false.
    end if    
  end if
! check for parameter conformity and homogeneity
  if (use_waterloops.EQV..true.) then
    dropout = .false.
    allocate(epsik(n-at(rsw1)%bb(1)+1))
    ii = b_type(at(rsw1)%bb(1))
    ixx = at(rsw1)%na*at(rsw1)%na
    do rs=1,rsw1-1
      if (disulf(rs).ge.rsw1) then
        write(ilog,*) 'Warning. Crosslinks involving residues with optimized interaction loops disable these loops.'
        use_waterloops = .false.
        exit
      end if
    end do
    do rs=rsw1,nseq
      if (use_waterloops.EQV..false.) exit
      if ((attyp(at(rs)%bb(1)).ne.bio_ljtyp(b_type(at(rs)%bb(1)))).OR.(b_type(at(rs)%bb(1)).ne.ii)) then
        write(ilog,*) 'Warning. Applied biotype or LJ patches imply that for chosen rigid water model optimized &
 &loops are not available.'
        use_waterloops = .false.
        exit
      end if
      if ((nrsintra(rs).gt.0).OR.(nrpolintra(rs).gt.0).OR.((nrsnb(rs).ne.ixx).AND.(rs.ne.nseq)).OR.&
 &        ((nrexpolnb(rs).gt.0).AND.((lrel_md.eq.2).OR.(lrel_md.eq.3))) ) then
        write(ilog,*) 'Warning. Interaction model is incompatible with the use of optimized water loops. Disabled.'
        use_waterloops = .false.
        exit
      end if
      if (disulf(rs).gt.0) then
        write(ilog,*) 'Warning. Crosslinks involving residues with optimized interaction loops disable these loops.'
        use_waterloops = .false.
        exit
      end if
      j = n-at(rs)%bb(1)+1
      if (is_3site.EQV..true.) then  
        do k=at(rs)%bb(2),at(rs)%bb(3)
          epsik(:) = lj_eps(attyp(k),attyp(at(rsw1)%bb(1):n))
          if (maxval(epsik).gt.0.0) then
            write(ilog,*) 'Warning. LJ parameters imply that for chosen 3-site water model optimized &
 &loops are not available.'
            dropout = .true.
            use_waterloops = .false.
            exit
          end if
        end do
        if (dropout.EQV..true.) exit
      else if (is_4site.EQV..true.) then
        do k=at(rs)%bb(2),at(rs)%bb(4)
          epsik(:) = lj_eps(attyp(k),attyp(at(rsw1)%bb(1):n))
          if (maxval(epsik).gt.0.0) then
            write(ilog,*) 'Warning. LJ parameters imply that for chosen 4-site water model optimized &
 &loops are not available.'
            dropout = .true.
            use_waterloops = .false.
            exit
          end if
        end do
        if (dropout.EQV..true.) exit
      else if (is_5site.EQV..true.) then
        do k=at(rs)%bb(2),at(rs)%bb(5)
          epsik(:) = lj_eps(attyp(k),attyp(at(rsw1)%bb(1):n))
          if (maxval(epsik).gt.0.0) then
            write(ilog,*) 'Warning. LJ parameters imply that for chosen 4-site water model optimized &
 &loops are not available.'
            dropout = .true.
            use_waterloops = .false.
            exit
          end if
        end do
        if (dropout.EQV..true.) exit
      else
        write(ilog,*) 'Fatal. Error in setting up specialized water loops. This &
 &is most likely an omission bug.'
        call fexit()
      end if
    end do
    deallocate(epsik)
  end if
  if ((use_cutoffs.EQV..false.).AND.((is_3site.EQV..true.).OR.(is_4site.EQV..true.).OR.(is_5site.EQV..true.))) then
    write(ilog,*) 'Warning. Simulation performance would benefit twofold from cutoffs because specialized &
 &water loops are disabled by lack of use of cutoffs.'
    use_waterloops = .false.
  end if
  if (use_POLAR.EQV..false.) then
    use_waterloops = .false.
  end if
  if (use_waterloops.EQV..true.) then
    do rs=1,nseq
      rs_nbl(rs)%wnbalsz = min(10,nseq-rsw1+1)
      rs_nbl(rs)%wtralsz = min(10,nseq-rsw1+1)
!      rs_nbl(rs)%wlralsz = min(10,nseq-rsw1+1)
      allocate(rs_nbl(rs)%wnb(rs_nbl(rs)%wnbalsz))
      allocate(rs_nbl(rs)%wnbtr(rs_nbl(rs)%wtralsz))
      allocate(rs_nbl(rs)%wtrsvec(rs_nbl(rs)%wtralsz,3))
!      allocate(rs_nbl(rs)%wnblr(rs_nbl(rs)%wlralsz))
    end do
  else 
    rsw1 = nseq + 1 ! restore (important!)
  end if
!
! provide summary
!
 23   format('Molecule type ',i3,' has ',i6,' residues and a mass of ',g16.8,' D. ')
 24   format('   The first example molecule is # ',i6,' (residues ',i6,'-',&
 &i6,').')
 25   format('   There are ',i6,' molecules of this type.')
 251  format('   There are ',i6,' molecules of this type, which gives a bulk concentration of ',g12.5,'mM.')
 252  format('   There are ',i6,' molecules of this type, which gives a bulk concentration of atmost ',g12.5,'mM.')
 26   format('Analysis group # ',i6,' is equivalent to molecule type # ',i6,'.')
 27   format('   These molecule(s) are (is) considered solvent.')
 28   format('   These molecule(s) are (is) considered solute(s).')
 29   format('Analysis group # ',i6,' is a subset of molecule type # ',i6,' and constitutes:')
 30   format('   Molecule # ',i6)
!
  if (seq_report.EQV..true.) then
    write(ilog,*)
    write(ilog,*) '--- Overview of sequence ---'
    do k=1,nmoltyp
      write(ilog,*)
      write(ilog,23) k,rsmol(moltyp(k,1),2)-rsmol(moltyp(k,1),1)+1,sum(mass(atmol(moltyp(k,1),1):atmol(moltyp(k,1),2)))
      write(ilog,24) moltyp(k,1),rsmol(moltyp(k,1),1),&
 &                                 rsmol(moltyp(k,1),2)
      if ((moltyp(k,2).gt.1).AND.((ens%flag.eq.1).OR.(ens%flag.eq.2)).AND.(bnd_type.eq.1).AND.(bnd_shape.eq.1)) then
        write(ilog,251) moltyp(k,2),moltyp(k,2)/(ens%insV*avogadro*1.0e-30) ! mmol/l
      else if ((moltyp(k,2).gt.1).AND.((ens%flag.eq.1).OR.(ens%flag.eq.2))) then
        write(ilog,252) moltyp(k,2),moltyp(k,2)/(ens%insV*avogadro*1.0e-30) ! mmol/l
      else
        write(ilog,25) moltyp(k,2)
      end if
      if (do_pol(k).EQV..false.) then
        write(ilog,*) '  Polymer analysis is turned off for this molec&
 &ule type.'
      else
        write(ilog,*) '  Polymer analysis is turned on for this molecu&
 &le type.'
      end if
      if (do_pers(k).EQV..false.) then
        write(ilog,*) '  Turns/angular statistics are turned off for t&
 &his molecule type.'
      else
        write(ilog,*) '  Turns/angular statistics are turned on for th&
 &is molecule type.'
      end if
      if (do_tors(k).EQV..false.) then
        write(ilog,*) '  Torsional analysis is turned off for this mol&
 &ecule type.'
      else
        write(ilog,*) '  Torsional analysis is turned on for this mole&
 &cule type.'
      end if
    end do
    write(ilog,*)
    write(ilog,*)
    write(ilog,*) '--- Overview of analysis groups ---'
    write(ilog,*)
    do k=1,nangrps
      if (moltyp(moltypid(molangr(k,1)),2).eq.molangr(k,2)) then
        write(ilog,26) k,moltypid(molangr(k,1))
        if (is_solvent(molangr(k,1)).EQV..true.) then
          write(ilog,27)
        else
          write(ilog,28)
        end if
      else
        write(ilog,29) k,moltypid(molangr(k,1))
        do j=1,nmol
          if (an_grp_mol(j).eq.k) then
            write(ilog,30) j
          end if
        end do
        if (is_solvent(molangr(k,1)).EQV..true.) then
          write(ilog,27)
        else
          write(ilog,28)
        end if
      end if
    end do
    write(ilog,*)
  end if
!
end
!
!-------------------------------------------------------------------------------
!
! initial structure randomization (may be silent --> globrandomize)
!
subroutine randomize_bb()
!
  use iounit
  use sequen
  use energies
  use atoms
  use molecule
  use cutoffs
  use aminos
  use ujglobals
  use polypep
  use pdb
  use system, ONLY: pdb_analyze,dyn_mode
  use params, ONLY: be_unsafe
  use fyoc, ONLY: disulf
!
  implicit none
!
  RTYPE et(MAXENERGYTERMS),scale_1,cut,dis,dis2,random,etval
  RTYPE dum3(3),ccmat(6),oldvals(12),solmat(MAXSOLU*4,6),intsvec(9),fixvals(11)
  integer i,j,cnt,azero,imol,jmol,aone,nlks,solucount,lastrs,firstrs,rs,rsi,k,kk,rs1,rs2,lkcnt,imol2,avsz,iacnt,pdbrsize
  integer i1,i2,i3,i4,i5,i6,i7,i8,i9,ii,fixed_scs,maxmol
  integer, ALLOCATABLE:: lklst(:,:),molmap(:)
  logical log_0,log_1,log_2,log_3,log_4,log_5,log_6,log_7,log_hs,log_8
  logical notdone,dorand,doclnk(3),badflag,haveone,do_ids(4),atrue,afalse,ldum
  logical, ALLOCATABLE:: donewith(:)
  RTYPE, ALLOCATABLE:: allvals(:)
!
! a bit of initialization
!
  if (n_pdbunk.gt.0) call zmatfyc2() ! pointer arrays may not be up2date
!
  azero = 0
  aone = 1
  atrue = .true.
  afalse = .false.
  scale_1 = scale_IPP
  log_hs = use_hardsphere
  scale_IPP = 1.0 ! note that this is only relevant if use_IPP is true, i.e., it's completely random if FMCSC_SC_IPP = 0.0
  use_hardsphere = .false.
  cut = mcnb_cutoff 
  mcnb_cutoff = 5.0
  log_0 = use_IPP
  log_1 = use_attLJ
  log_2 = use_CORR
  log_3 = use_IMPSOLV
  log_4 = use_WCA
  log_5 = use_POLAR
  log_6 = use_ZSEC
  log_7 = use_TOR
  log_8 = use_FEG
  if (use_WCA.EQV..true.) then
    use_IPP = .true. ! pointless to randomize a structure without EV if later EV is added
  end if
  use_attLJ = .false.
  use_CORR = .false.  
  use_IMPSOLV = .false.
  use_WCA = .false.
  use_POLAR = .false.
  use_ZSEC = .false.
  use_TOR = .false.
  use_FEG = .false.
  et(:) = 0.0
 66 format(' Warning. Failed to build a low energy conformation for segment from residues ',i7,' (',a3,') to ',i7,' (',a3,').')
 67 format(' Warning. Failed to place molecule ',i7,' in a low energy manner (res. ',i7,' (',a3,') to ',i7,' (',a3,')).')
 77 format(' A starting structure with high forces is fatal in gradient-based calculations without prior relaxation (use &
 &FMCSC_UNSAFE to avoid this termination). Increasing the number of attempts (FMCSC_RANDOMATTS) might improve chances.')
!
  ldum = .false.
  if (globrandomize.eq.2) ldum = .true.
  allocate(molmap(nmol))
  molmap(:) = 0
  maxmol = 0
!
  pdbrsize = 0
  if (allocated(pdb_didread).EQV..true.) then
    pdbrsize = size(pdb_didread,dim=1)
    pdb_didread(:,3) = 0
  end if
  write(ilog,*)
  write(ilog,'(a)') '--- Now entering initial structure generation/augmentation routine ---'
!
! first remove clashes of randomly rebuilt, non CRLK-sidechains (if any) with read-in structures
  if (allocated(pdbmap2).EQV..true.) then
    do imol=1,pdbrsize
      do rs=pdb_didread(imol,1),pdb_didread(imol,2) ! rsmol(imol,1),rsmol(imol,2)
        if (disulf(rs).ne.0) cycle
        notdone = .false.
        do ii=1,at(rs)%nsc
          i = at(rs)%sc(ii)
          if ((pdbmap2(i).le.0).AND.((n12(i).gt.1).OR.(mass(i).gt.6.0))) then
            notdone = .true.
          end if
        end do
        if (notdone.EQV..true.) then
          pdb_didread(imol,3) = pdb_didread(imol,3) + 1
          if ((globrandomize.eq.1).OR.(globrandomize.eq.2)) call MC_randomize_sc(imol,rs,ldum,badflag)
        end if
      end do
    end do
  end if
!
! next built all missing tails onto existing partial internal conformation read from PDB
! restrict interaction checking to intramoelcular if RBC are later randomized anyway
  do imol=1,pdbrsize
!   check/warn intramolecular crosslinks
    do i=1,n_crosslinks
      if ((molofrs(crosslink(i)%rsnrs(1)).eq.imol).AND.(molofrs(crosslink(i)%rsnrs(2)).eq.imol)) then
        if ((minval(crosslink(i)%rsnrs(1:2)).lt.pdb_didread(imol,1)).AND.&
 &          (maxval(crosslink(i)%rsnrs(1:2)).gt.pdb_didread(imol,2))) then
          if ((globrandomize.eq.0).OR.(globrandomize.eq.3)) then
           write(ilog,*) 'Warning. Intramolecular crosslink spanning two extended tails in the same molecule will not &
 &be satisfied (residues ',crosslink(i)%rsnrs(1),' and ',crosslink(i)%rsnrs(2),').'
          else
            write(ilog,*) 'Warning. Intramolecular crosslink spanning two randomized tails in the same molecule will at best &
 &be obtained approximately (residues ',crosslink(i)%rsnrs(1),' and ',crosslink(i)%rsnrs(2),').'
          end if
        else if ((minval(crosslink(i)%rsnrs(1:2)).ge.pdb_didread(imol,1)).AND.&
 &               (maxval(crosslink(i)%rsnrs(1:2)).le.pdb_didread(imol,2))) then
          write(ilog,*) 'Warning. Existing intramolecular crosslink will be retained from structural input &
 &(residues ',crosslink(i)%rsnrs(1),' and ',crosslink(i)%rsnrs(2),').'
        else
          if ((globrandomize.eq.0).OR.(globrandomize.eq.3)) then
            write(ilog,*) 'Warning. Intramolecular crosslink involving an extended (default) tail will not be &
 &satisfied (residues ',crosslink(i)%rsnrs(1),' and ',crosslink(i)%rsnrs(2),').'
          else
            write(ilog,*) 'Warning. Intramolecular crosslink involving a randomized tail will at best &
 &be obtained approximately (residues ',crosslink(i)%rsnrs(1),' and ',crosslink(i)%rsnrs(2),').'
          end if
        end if
      end if
    end do
    if ((globrandomize.eq.0).OR.(globrandomize.eq.3)) exit
    if (pdb_didread(imol,1).gt.rsmol(imol,1)) then
      call MC_randomize(imol,rsmol(imol,1),pdb_didread(imol,1)-1,atrue,atrue,ldum,badflag)
      if (badflag.EQV..true.) then
        write(ilog,66) rsmol(imol,1),amino(seqtyp(rsmol(imol,1))),pdb_didread(imol,1)-1,amino(seqtyp(pdb_didread(imol,1)-1))
        if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
          write(ilog,77)
          if (be_unsafe.EQV..false.) call fexit()
        end if
      end if
    end if
    if (pdb_didread(imol,2).lt.rsmol(imol,2)) then
      call MC_randomize(imol,pdb_didread(imol,2)+1,rsmol(imol,2),afalse,atrue,ldum,badflag)
      if (badflag.EQV..true.) then
        write(ilog,66) pdb_didread(imol,2)+1,amino(seqtyp(pdb_didread(imol,2)+1)),rsmol(imol,2),amino(seqtyp(rsmol(imol,2)))
        if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
          write(ilog,77)
          if (be_unsafe.EQV..false.) call fexit()
        end if
      end if
    end if
  end do
  if ((globrandomize.ne.2).AND.(globrandomize.ne.3)) then
    molmap(1:pdbrsize) = 1
    maxmol = pdbrsize
  end if
!
! now loop over all molecules 
  do imol=1,nmol
    do_ids(:) = .true.
    if (pdbrsize.ge.imol) then
      if (allocated(pdbmap2).EQV..true.) then
!       molecules with completely missing reference frames that had sequence entries, may cause FPE
!       without displacement -> force RB randomization
        if ((pdbmap2(atmol(imol,1)).le.0).AND.(pdbmap2(min(atmol(imol,2),atmol(imol,1)+1)).le.0).AND.&
 &        (pdbmap2(min(atmol(imol,2),atmol(imol,1)+2)).le.0).AND.(rsmol(imol,1).eq.rsmol(imol,2))) then
          do_ids(1) = .false.
        else if (maxval(pdbmap2(atmol(imol,1):atmol(imol,2))).le.0) then
          do_ids(1) = .false.
        else
          do_ids(1) = .false.
          if ((globrandomize.ne.2).AND.(globrandomize.ne.3)) do_ids(2) = .false.
        end if
      else       
        do_ids(1) = .false.
        if ((globrandomize.ne.2).AND.(globrandomize.ne.3)) do_ids(2) = .false.
      end if
    else
      if ((globrandomize.eq.0).OR.(globrandomize.eq.3)) do_ids(1) = .false.
    end if
    if (ntormol(moltypid(imol)).le.0) then
      do_ids(1) = .false.
    end if
!
    nlks = 0
    if (rsmol(imol,2).gt.rsmol(imol,1)) then
      allocate(lklst(rsmol(imol,2)-rsmol(imol,1)+1,3))
      firstrs = rsmol(imol,2)
      lastrs = rsmol(imol,1)
      do i=1,n_crosslinks
        if ((molofrs(crosslink(i)%rsnrs(1)).eq.imol).AND.&
 &          (molofrs(crosslink(i)%rsnrs(2)).eq.imol)) then
          nlks = nlks + 1
          lklst(nlks,1) = crosslink(i)%rsnrs(1)
          lklst(nlks,2) = crosslink(i)%rsnrs(2)
          lklst(nlks,3) = i
          if (lklst(nlks,1).lt.firstrs) firstrs = lklst(nlks,1)
          if (lklst(nlks,1).gt.lastrs) lastrs = lklst(nlks,1)
          if (lklst(nlks,2).lt.firstrs) firstrs = lklst(nlks,2)
          if (lklst(nlks,2).gt.lastrs) lastrs = lklst(nlks,2)
          do j=1,nlks-1
            if (lklst(nlks,1).lt.lklst(j,1)) then
              write(ilog,*) 'Fatal. List of crosslinks is not in expected order. This is a bug.'
              call fexit()
            end if
          end do
        end if
      end do
    end if
    rs1 = rsmol(imol,1)
     
    if (do_ids(1).EQV..true.) then
!   building randomized conformations gets extremely messy for crosslinked peptides
!   due to the implied loop closures; without them, it's trivial
      if (nlks.eq.0) then
        call stretch_randomize(imol,rs1,rsmol(imol,2),badflag)
        if (badflag.EQV..true.) then
          write(ilog,66) rs1,amino(seqtyp(rs1)),rsmol(imol,2),amino(seqtyp(rsmol(imol,2)))
          if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
            write(ilog,77)
            if (be_unsafe.EQV..false.) call fexit()
          end if
        end if
      else
        call stretch_randomize(imol,rs1,firstrs-1,badflag)
        if (badflag.EQV..true.) then
           write(ilog,66) rs1,amino(seqtyp(rs1)),firstrs-1,amino(seqtyp(firstrs-1))
          if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
            write(ilog,77)
            if (be_unsafe.EQV..false.) call fexit()
          end if
        end if
        k = 0
        allocate(donewith(nlks))
        donewith(:) = .false.
        do while (k.lt.nlks)
!         pick a crosslink whose 2nd residue is within another crosslink stretch
          notdone = .true.
          do i=1,nlks
            if (donewith(i).EQV..true.) cycle
            do j=1,nlks
              if ((lklst(i,2).lt.lklst(j,2)).AND.(lklst(i,2).gt.lklst(j,1))) then
                notdone = .false.
                if (donewith(j).EQV..true.) then
                  write(ilog,*) 'Fatal. Crosslink topology is too complex for CAMPARI&
   & to close all loops consistently. Please check back later.'
                  call fexit()
                end if
                exit
              end if
            end do
            if (notdone.EQV..false.) exit
          end do
          if (notdone.EQV..true.) then
            do i=1,nlks
              if (donewith(i).EQV..false.) exit
            end do
          end if
!
!         now sample whatever we can randomly, then close crosslink via concerted rotation
          notdone = .true.
          kk = 0
          if (allocated(allvals).EQV..true.) deallocate(allvals)
          avsz = at(lklst(i,2))%bb(1)-at(lklst(i,1))%bb(1)+at(lklst(i,2))%nbb+at(lklst(i,2))%nsc+3
          allocate(allvals(avsz))
          allvals(1) = HUGE(allvals(1))
          allvals(avsz) = allvals(1)
          haveone = .false.
          do while (notdone.EQV..true.)
            kk = kk + 1
            if (kk.gt.globrdmatts) then
              write(ilog,66) lklst(i,1),amino(seqtyp(lklst(i,1))),lklst(i,2),amino(seqtyp(lklst(i,2)))
              write(ilog,*) 'This is caused by intramolecular crosslinks that could not be closed clash-free or even at all &
   &(molecule ',imol,'). It may be necessary to simplify crosslinks.'
              if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
                write(ilog,77)
                if (be_unsafe.EQV..false.) call fexit()
              end if
              if (haveone.EQV..true.) call apply_forstretch(lklst(i,1),lklst(i,2),allvals,avsz)
              exit
            end if
            do rs=lklst(i,1),lklst(i,2)
              dorand = .true.
              do j=1,nlks
                if (donewith(j).EQV..false.) cycle
                if ((rs.ge.lklst(j,1)).AND.(rs.le.lklst(j,2))) then
                  dorand = .false.
                  exit
                end if
              end do
              if (dorand.EQV..false.) cycle
              call stretch_randomize(imol,rs,rs,badflag)
!              if (badflag.EQV..true.) write(ilog,66) rs,amino(seqtyp(rs)),rs,amino(seqtyp(rs))
            end do
            call crosslink_indices(lklst(i,3),i1,i2,i3,i4,i5,i6,i7,i8,i9)
            call crosslink_refvals(lklst(i,3),i1,i2,i3,i4,i5,i6,i7,i8,i9,ccmat,oldvals)
            call torchainclose(i1,i2,i3,i7,i8,i9,ccmat,solmat,oldvals,solucount)
            if (solucount.ge.1) then
              haveone = .true.
              call crosslink_picksolu(lklst(i,3),i1,i2,i3,i4,i5,i6,i7,i8,i9,solucount,solmat,ccmat,oldvals,allvals,avsz)
              if (allvals(avsz).le.globrdmthresh) then
                call apply_forstretch(lklst(i,1),lklst(i,2),allvals,avsz)
                notdone = .false.
              end if
            end if
          end do
          donewith(i) = .true.
          k = k + 1
        end do
        if (allocated(allvals).EQV..true.) deallocate(allvals)
        deallocate(donewith)
        rs = firstrs+1
        do while (rs.lt.lastrs-1)
          dorand = .true.
          do i=1,nlks
            if ((rs.ge.lklst(i,1)).AND.(rs.le.lklst(i,2))) then
              dorand = .false.
              exit
            end if
          end do
          if (dorand.EQV..true.) then
            rsi = rs
            notdone = .true.
            do while (notdone.EQV..true.)
              rs = rs + 1
              if (rs.eq.lastrs-1) notdone = .false.
              do i=1,nlks
                if ((rs.ge.lklst(i,1)).AND.(rs.le.lklst(i,2))) then
                  notdone = .false.
                  rs = rs - 1
                end if
              end do
            end do
            dorand = .false.
            call stretch_randomize(imol,rsi,rs,badflag)
            if (badflag.EQV..true.) then
              write(ilog,66) rsi,amino(seqtyp(rsi)),rs,amino(seqtyp(rs))
              write(ilog,*) 'This is likely related to intramolecular crosslinks in molecule ',imol,'.'
              if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
                write(ilog,77)
                if (be_unsafe.EQV..false.) call fexit()
              end if
            end if
            rs = rs + 1
          else
            rs = rs + 1
          end if
        end do
        call stretch_randomize(imol,lastrs+1,rsmol(imol,2),badflag)
        if (badflag.EQV..true.) then
          write(ilog,66) lastrs+1,amino(seqtyp(lastrs+1)),rsmol(imol,2),amino(seqtyp(rsmol(imol,2)))
          write(ilog,*) 'This is likely related to intramolecular crosslinks in molecule ',imol,'.'
          if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
            write(ilog,77)
            if (be_unsafe.EQV..false.) call fexit()
          end if
        end if
      end if
    end if
    if (allocated(allvals).EQV..true.) deallocate(allvals)
    if (allocated(lklst).EQV..true.) deallocate(lklst)
!
    call update_rigid(imol)
    call update_rigidm(imol)
    imol2 = imol
!
    doclnk(:) = .false.
    do i=1,n_crosslinks
      rs1 = crosslink(i)%rsnrs(1)
      imol2 = molofrs(rs1)
      rs2 = crosslink(i)%rsnrs(2)
      jmol = molofrs(rs2) ! should always be larger than imol2
      if ((jmol.eq.imol).AND.(jmol.ne.imol2)) then
        if (doclnk(1).EQV..true.) then
          write(ilog,*) 'Warning. Intermolecular crosslinks having the same target molecule may cause &
 &CAMPARI to crash. If possible, it is better to put the molecule with multiple connections first in the &
 &sequence. Intermolecular crosslink between residues ',crosslink(lkcnt)%rsnrs(1),' and ',crosslink(lkcnt)%rsnrs(2),&
 &' will be ignored in any case.' 
        end if
        doclnk(1) = .true.
        lkcnt = i
      else if ((imol2.eq.imol).AND.(jmol.ne.imol2)) then
        doclnk(2) = .true.
      end if
    end do
!   if RBC are to remain fixed for this molecule (which only happens if it is structurally resolved (at least partially)
!   and RANDOMIZE was not set to 2 or 3)
    if (doclnk(1).EQV..true.) then
      do i=1,n_crosslinks
        rs1 = crosslink(i)%rsnrs(1)
        imol2 = molofrs(rs1)
        rs2 = crosslink(i)%rsnrs(2)
        jmol = molofrs(rs2)
        if (imol2.eq.jmol) cycle ! intra-crosslink
        if ((molmap(imol2).gt.0).AND.(molmap(jmol).gt.0)) cycle ! handled already
        if ((molmap(imol2).le.0).AND.(molmap(jmol).le.0)) cycle ! nothing to do yet
!
        if (((molmap(max(imol2,jmol)).gt.0).AND.(molmap(min(imol2,jmol)).le.0).AND.(min(imol2,jmol).lt.imol)).OR.&
 &          ((max(imol2,jmol).eq.imol).AND.(molmap(max(imol2,jmol)).le.0))) then
          if ((max(imol2,jmol).eq.imol).AND.(molmap(max(imol2,jmol)).le.0)) then
            imol2 = imol
          else
            imol2 = min(imol2,jmol)
          end if
          if (do_ids(2).EQV..true.) then
            if (crosslink(i)%itstype.eq.2) then
              intsvec(3) = 2.03 ! h.c. -S-S
              intsvec(6) = 1.822 ! h.c. - match with sidechain
              intsvec(9) = 1.53  ! h.c. - match with sidechain
              intsvec(2) = 103.0 ! h.c. -CB-S-S
              intsvec(5) = 103.0 ! h.c. -S-S-CB
              intsvec(8) = 114.4 ! h.c. - match with sidechain
            else
              write(ilog,*) 'Fatal. Encountered unsupported crosslink type in randomize_bb(...).&
         & This is most likely an omission bug.'
              call fexit()
            end if
            fixvals(:) = 0.0
            fixvals(10:11) = HUGE(fixvals(11))
            cnt = 0
!
            do while (fixvals(11).gt.globrdmthresh)
              cnt = cnt + 1
              if (cnt.gt.globrdmatts) then
                call crosslink_follow(i,fixvals(1:9)) ! re-apply best solution
                write(ilog,67) imol2,rsmol(imol2,1),amino(seqtyp(rsmol(imol2,1))),rsmol(imol2,2),amino(seqtyp(rsmol(imol2,2)))
                write(ilog,*) 'This is related to a problematic intermolecular crosslink.'
                if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
                  write(ilog,77)
                  if (be_unsafe.EQV..false.) call fexit()
                end if
                exit
              end if
              if (crosslink(i)%itstype.eq.2) then
                intsvec(1) = random()*360.0 - 180.0
                intsvec(4) = random()*360.0 - 180.0
                intsvec(7) = random()*360.0 - 180.0
              else
                write(ilog,*) 'Fatal. Encountered unsupported crosslink type in randomize_bb(...).&
         & This is most likely an omission bug.'
                call fexit()
              end if
              if (cnt.eq.1) fixvals(1:9) = intsvec(1:9) ! be safe
              call crosslink_follow(i,intsvec)
!             make sure intermolecular energies (only backward) are sane
              call evaluate_formol(imol2,maxmol,intsvec,fixvals,molmap,atrue)
            end do
          end if
          molmap(imol2) = 1
          maxmol = max(maxmol,imol2)
          call update_rigid(imol2)
          call update_rigidm(imol2)
        end if
      end do
!
    else
!
      imol2 = imol
!
      if (do_ids(2).EQV..true.) then
!
        cnt = 0
        fixvals(:) = 0.0
        fixvals(10:11) = HUGE(fixvals(11))
        do while (fixvals(11).gt.globrdmthresh)
          cnt = cnt + 1
!         obviously we got into a mess: clear-up
          if (cnt.gt.globrdmatts) then
            call apply_formol(imol2,fixvals(1:9))
            write(ilog,67) imol2,rsmol(imol2,1),amino(seqtyp(rsmol(imol2,1))),rsmol(imol2,2),amino(seqtyp(rsmol(imol2,2)))
            if (doclnk(2).EQV..false.) then
              write(ilog,*) 'This is a density or size problem.'
            else
              write(ilog,*) 'This is likely the result of an intermolecular crosslink that can not be satisfied exactly due &
 &to partial PDB input or complicated topology (check other warnings). It could also be a density or size problem.'
            end if
            if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
              write(ilog,77)
              if (be_unsafe.EQV..false.) call fexit()
            end if
            write(ilog,*) 'WARNING: Could not arrange molecules &
     &with non-dramatic overlap (for molecule ',imol2,').'
            exit
          end if
!         randomize pos for molecule imol2
          call randomize_trans(imol2,dum3)
          do i=atmol(imol2,1),atmol(imol2,2)
            x(i) = x(i) + dum3(1)
            y(i) = y(i) + dum3(2)
            z(i) = z(i) + dum3(3)
          end do
          call randomize_rot(imol2)
!         make sure intermolecular energies (only backward) are sane
          call evaluate_formol(imol2,maxmol,intsvec,fixvals,molmap,afalse)
        end do
      end if
      molmap(imol2) = 1
      maxmol = max(maxmol,imol2)
      call update_rigid(imol2)
      call update_rigidm(imol2)
    end if
  end do
!
  if (allocated(molmap).EQV..true.) deallocate(molmap)
  write(ilog,'(a)') '--- Now exiting initial structure generation/augmentation routine ---'
  write(ilog,*)
! 
!
! finalize a bit
  use_IPP = log_0
  use_attLJ = log_1
  use_CORR = log_2
  use_IMPSOLV = log_3
  use_WCA = log_4
  use_POLAR = log_5
  use_ZSEC = log_6
  use_TOR = log_7
  use_FEG = log_8
  scale_IPP = scale_1
  use_hardsphere = log_hs
  mcnb_cutoff = cut
!
end
!!
!!-------------------------------------------------------------------------------
!!
!! below is the old function (commented), which is obsolete
!!
!subroutine randomize2_bb()
!!
!  use iounit
!  use sequen
!  use energies
!  use atoms
!  use molecule
!  use cutoffs
!  use aminos
!  use ujglobals
!  use polypep
!  use pdb
!  use system, ONLY: pdb_analyze,dyn_mode
!  use params, ONLY: be_unsafe
!!
!  implicit none
!!
!  RTYPE et(MAXENERGYTERMS),scale_1,cut,dis,dis2,random,etval
!  RTYPE dum3(3),ccmat(6),oldvals(12),solmat(MAXSOLU*4,6),intsvec(9)
!  integer i,j,cnt,azero,imol,jmol,aone,nlks,solucount,lastrs,firstrs,rs,rsi,k,kk,rs1,rs2,frs,lkcnt,imol2,avsz,iacnt
!  integer i1,i2,i3,i4,i5,i6,i7,i8,i9
!  integer, ALLOCATABLE:: lklst(:,:)
!  logical log_0,log_1,log_2,log_3,log_4,atrue,log_5,log_6,log_7,log_hs,log_8,notdone,dorand,doclnk(3),badflag,haveone
!  logical, ALLOCATABLE:: donewith(:)
!  RTYPE, ALLOCATABLE:: allvals(:)
!!
!! a bit of initialization
!!
!  if (n_pdbunk.gt.0) call zmatfyc2() ! pointer arrays may not be up2date
!!
!  azero = 0
!  aone = 1
!  atrue = .true.
!  scale_1 = scale_IPP
!  log_hs = use_hardsphere
!  scale_IPP = 1.0 ! note that this is only relevant if use_IPP is true, i.e., it's completely random if FMCSC_SC_IPP = 0.0
!  use_hardsphere = .false.
!  cut = mcnb_cutoff 
!  mcnb_cutoff = 5.0
!  log_0 = use_IPP
!  log_1 = use_attLJ
!  log_2 = use_CORR
!  log_3 = use_IMPSOLV
!  log_4 = use_WCA
!  log_5 = use_POLAR
!  log_6 = use_ZSEC
!  log_7 = use_TOR
!  log_8 = use_FEG
!  if (use_WCA.EQV..true.) then
!    use_IPP = .true. ! pointless to randomize a structure without EV if later EV is added
!  end if
!  use_attLJ = .false.
!  use_CORR = .false.  
!  use_IMPSOLV = .false.
!  use_WCA = .false.
!  use_POLAR = .false.
!  use_ZSEC = .false.
!  use_TOR = .false.
!  use_FEG = .false.
!  et(:) = 0.0
!  if (pdb_ihlp.le.0) then
!    frs = nseq
!  else
!    frs = pdb_ihlp
!  end if
!  do imol=1,size(pdb_didread,dim=1)
!    write(*,*) imol,pdb_didread(imol,1:2)
!  end do
! 66 format(' Warning. Failed to build a clash-free segment for residues ',i7,' (',a3,') to ',i7,' (',a3,').')
! 67 format(' Warning. Failed to place molecule ',i7,' in a clash-free manner (res. ',i7,' (',a3,') to ',i7,' (',a3,')).')
! 77 format(' A starting structure with clashes is fatal in gradient-based calculations without prior relaxation (use FMCSC_UNSAFE &
! &to avoid this termination). Increase tolerance and attempts (FMCSC_RANDOMTHRESH and FMCSC_RANDOMATTS) to improve chances.')
!!
!  do imol=1,nmol
!!   for no internal randomization exit right away
!    if ((globrandomize.eq.0).OR.(globrandomize.eq.2)) exit
!    if ((globrandomize.ge.4).AND.(rsmol(imol,2).le.frs)) cycle
!!
!    allocate(lklst(rsmol(imol,2)-rsmol(imol,1)+1,3))
!    nlks = 0
!    firstrs = rsmol(imol,2)
!    lastrs = rsmol(imol,1)
!    do i=1,n_crosslinks
!      if ((molofrs(crosslink(i)%rsnrs(1)).eq.imol).AND.&
! &        (molofrs(crosslink(i)%rsnrs(2)).eq.imol)) then
!        nlks = nlks + 1
!        lklst(nlks,1) = crosslink(i)%rsnrs(1)
!        lklst(nlks,2) = crosslink(i)%rsnrs(2)
!        lklst(nlks,3) = i
!        if ((minval(lklst(nlks,1:2)).le.frs).AND.(maxval(lklst(nlks,1:2)).gt.frs)) then
!          write(ilog,*) 'Warning. Intramolecular crosslink is ignored during partial &
! &structure randomization (residues ',lklst(nlks,1),' and ',lklst(nlks,2),').'
!          nlks = nlks - 1
!          cycle
!        end if
!        if ((minval(lklst(nlks,1:2)).le.frs).AND.(maxval(lklst(nlks,1:2)).le.frs).AND.(globrandomize.ne.1)) then
!          write(ilog,*) 'Warning. Existing intramolecular crosslink is retained during partial &
! &structure randomization (residues ',lklst(nlks,1),' and ',lklst(nlks,2),').'
!          nlks = nlks - 1
!          cycle
!        end if
!        if (lklst(nlks,1).lt.firstrs) firstrs = lklst(nlks,1)
!        if (lklst(nlks,1).gt.lastrs) lastrs = lklst(nlks,1)
!        if (lklst(nlks,2).lt.firstrs) firstrs = lklst(nlks,2)
!        if (lklst(nlks,2).gt.lastrs) lastrs = lklst(nlks,2)
!        do j=1,nlks-1
!          if (lklst(nlks,1).lt.lklst(j,1)) then
!            write(ilog,*) 'Fatal. List of crosslinks is not in expected order. This is a bug.'
!            call fexit()
!          end if
!        end do
!      end if
!    end do
!    if ((rsmol(imol,1).le.frs).AND.(rsmol(imol,2).gt.frs)) then
!      rs1 = frs + 1
!    else
!      rs1 = rsmol(imol,1)
!    end if
!!   building randomized conformations gets extremely messy for crosslinked peptides
!!   due to the implied loop closures; without them, it's trivial
!    if (nlks.eq.0) then
!      call stretch_randomize(imol,rs1,rsmol(imol,2),badflag)
!      if (badflag.EQV..true.) write(ilog,66) rs1,amino(seqtyp(rs1)),rsmol(imol,2),amino(seqtyp(rsmol(imol,2)))
!      if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!        write(ilog,77)
!        if (be_unsafe.EQV..false.) call fexit()
!      end if
!    else
!      call stretch_randomize(imol,rs1,firstrs-1,badflag)
!      if (badflag.EQV..true.) write(ilog,66) rsi,amino(seqtyp(rsi)),firstrs-1,amino(seqtyp(firstrs-1))
!      if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!        write(ilog,77)
!        if (be_unsafe.EQV..false.) call fexit()
!      end if
!      k = 0
!      allocate(donewith(nlks))
!      donewith(:) = .false.
!      do while (k.lt.nlks)
!!       pick a crosslink whose 2nd residue is within another crosslink stretch
!        notdone = .true.
!        do i=1,nlks
!          if (donewith(i).EQV..true.) cycle
!          do j=1,nlks
!            if ((lklst(i,2).lt.lklst(j,2)).AND.(lklst(i,2).gt.lklst(j,1))) then
!              notdone = .false.
!              if (donewith(j).EQV..true.) then
!                write(ilog,*) 'Fatal. Crosslink topology is too complex for CAMPARI&
! & to close all loops consistently. Please check back later.'
!                call fexit()
!              end if
!              exit
!            end if
!          end do
!          if (notdone.EQV..false.) exit
!        end do
!        if (notdone.EQV..true.) then
!          do i=1,nlks
!            if (donewith(i).EQV..false.) exit
!          end do
!        end if
!!
!!       now sample whatever we can randomly, then close crosslink via concerted rotation
!        notdone = .true.
!        kk = 0
!        if (allocated(allvals).EQV..true.) deallocate(allvals)
!        avsz = at(lklst(i,2))%bb(1)-at(lklst(i,1))%bb(1)+at(lklst(i,2))%nbb+at(lklst(i,2))%nsc+2
!        allocate(allvals(avsz))
!        allvals(1) = HUGE(allvals(1))
!        haveone = .false.
!        do while (notdone.EQV..true.)
!          kk = kk + 1
!          if (kk.gt.globrdmatts) then
!            write(ilog,66) lklst(i,1),amino(seqtyp(lklst(i,1))),lklst(i,2),amino(seqtyp(lklst(i,2)))
!            write(ilog,*) 'This is caused by intramolecular crosslinks that could not be closed clash-free or even at all &
! &(molecule ',imol,'). It may be necessary to simplify crosslinks.'
!            if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!              write(ilog,77)
!              if (be_unsafe.EQV..false.) call fexit()
!            end if
!            if (haveone.EQV..true.) call apply_forstretch(lklst(i,1),lklst(i,2),allvals,avsz)
!            exit
!          end if
!          do rs=lklst(i,1),lklst(i,2)
!            dorand = .true.
!            do j=1,nlks
!              if (donewith(j).EQV..false.) cycle
!              if ((rs.ge.lklst(j,1)).AND.(rs.le.lklst(j,2))) then
!                dorand = .false.
!                exit
!              end if
!            end do
!            if (dorand.EQV..false.) cycle
!            call stretch_randomize(imol,rs,rs,badflag)
!!            if (badflag.EQV..true.) write(ilog,66) rs,amino(seqtyp(rs)),rs,amino(seqtyp(rs))
!          end do
!          call crosslink_indices(lklst(i,3),i1,i2,i3,i4,i5,i6,i7,i8,i9)
!          call crosslink_refvals(lklst(i,3),i1,i2,i3,i4,i5,i6,i7,i8,i9,ccmat,oldvals)
!          call torchainclose(i1,i2,i3,i7,i8,i9,ccmat,solmat,oldvals,solucount)
!          if (solucount.ge.1) then
!            haveone = .true.
!            call crosslink_picksolu(lklst(i,3),i1,i2,i3,i4,i5,i6,i7,i8,i9,solucount,solmat,ccmat,oldvals,allvals,avsz)
!            if (allvals(1).le.globrdmthresh) then
!              call apply_forstretch(lklst(i,1),lklst(i,2),allvals,avsz)
!              notdone = .false.
!            end if
!          end if
!        end do
!        donewith(i) = .true.
!        k = k + 1
!      end do
!      deallocate(allvals)
!      deallocate(donewith)
!      rs = firstrs+1
!      do while (rs.lt.lastrs-1)
!        dorand = .true.
!        do i=1,nlks
!          if ((rs.ge.lklst(i,1)).AND.(rs.le.lklst(i,2))) then
!            dorand = .false.
!            exit
!          end if
!        end do
!        if (dorand.EQV..true.) then
!          rsi = rs
!          notdone = .true.
!          do while (notdone.EQV..true.)
!            rs = rs + 1
!            if (rs.eq.lastrs-1) notdone = .false.
!            do i=1,nlks
!              if ((rs.ge.lklst(i,1)).AND.(rs.le.lklst(i,2))) then
!                notdone = .false.
!                rs = rs - 1
!              end if
!            end do
!          end do
!          dorand = .false.
!          call stretch_randomize(imol,rsi,rs,badflag)
!          if (badflag.EQV..true.) write(ilog,66) rsi,amino(seqtyp(rsi)),rs,amino(seqtyp(rs))
!          if (badflag.EQV..true.) write(ilog,*) 'This is likely related to intramolecular crosslinks in molecule ',imol,'.'
!          if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!            write(ilog,77)
!            if (be_unsafe.EQV..false.) call fexit()
!          end if
!          rs = rs + 1
!        else
!          rs = rs + 1
!        end if
!      end do
!      call stretch_randomize(imol,lastrs+1,rsmol(imol,2),badflag)
!      if (badflag.EQV..true.) write(ilog,66) lastrs+1,amino(seqtyp(lastrs+1)),rsmol(imol,2),amino(seqtyp(rsmol(imol,2)))
!      if (badflag.EQV..true.) write(ilog,*) 'This is likely related to intramolecular crosslinks in molecule ',imol,'.'
!      if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!        write(ilog,77)
!        if (be_unsafe.EQV..false.) call fexit()
!      end if
!    end if
!    deallocate(lklst)
!  end do
!  if (allocated(allvals).EQV..true.) deallocate(allvals)
!!
!  do imol=1,nmol
!    call update_rigid(imol)
!  end do
!!
!! now change the rigid-body degrees of freedom to eliminate intermolecular clashes
!  do imol=1,nmol
!    doclnk(:) = .false.
!    do i=1,n_crosslinks
!      rs1 = crosslink(i)%rsnrs(1)
!      imol2 = molofrs(rs1)
!      rs2 = crosslink(i)%rsnrs(2)
!      jmol = molofrs(rs2) ! should always be larger than imol2
!      if ((jmol.eq.imol).AND.(jmol.ne.imol2)) then
!        if (doclnk(1).EQV..true.) then
!          write(ilog,*) 'Warning. Intermolecular crosslinks having the same target molecule may cause &
! &CAMPARI to crash. At most the molecule latest molecule/crosslink in the sequence will behave as expected during &
! &randomization.' 
!        end if
!        doclnk(1) = .true.
!        lkcnt = i
!      end if
!    end do
!    if (doclnk(1).EQV..true.) then
!!
!      if (imol2.le.size(pdb_didread,dim=1)) then
!        if ((rs1.lt.pdb_didread(imol2,1)).OR.(rs1.gt.pdb_didread(imol2,2))) then
!          doclnk(2) = .true. ! missing from input but some of the molecule there
!        end if
!      end if
!      if (jmol.le.size(pdb_didread,dim=1)) then
!        if ((rs2.lt.pdb_didread(jmol,1)).OR.(rs2.gt.pdb_didread(jmol,2))) then
!          doclnk(3) = .true. ! missing from input but some of the molecule there
!        end if
!      end if
!      if ((doclnk(2).EQV..true.).OR.(doclnk(3).EQV..true.)) then
!        write(ilog,*) 'Fatal. The presence of partial coordinates for multiple chains connected by an intermolecular &
! &crosslink without the crosslinked residues resolved is a loop closure problem not currently supported by CAMPARI. &
! &Disable randomization (FCMSC_RANDOMIZE) if a conformation with an unsatisfied intermolecular crosslink is desired.'
!        call fexit()
!      end if
!      if (((globrandomize.eq.0).OR.(globrandomize.eq.3).OR.(globrandomize.eq.4)).AND.(imol2.le.size(pdb_didread,dim=1)).AND.&
! (jmol.le.size(pdb_didread,dim=1))) cycle ! both must have been read -> do not randomize
!      if (crosslink(lkcnt)%itstype.eq.2) then
!        intsvec(3) = 2.03 ! h.c. -S-S
!        intsvec(6) = 1.822 ! h.c. - match with sidechain
!        intsvec(9) = 1.53  ! h.c. - match with sidechain
!        intsvec(2) = 103.0 ! h.c. -CB-S-S
!        intsvec(5) = 103.0 ! h.c. -S-S-CB
!        intsvec(8) = 114.4 ! h.c. - match with sidechain
!      else
!        write(ilog,*) 'Fatal. Encountered unsupported crosslink type in randomize_bb(...).&
!   & This is most likely an omission bug.'
!        call fexit()
!      end if
!      etval = HUGE(etval)
!      cnt = 0
!!
!      do while (etval.gt.globrdmthresh)
!        cnt = cnt + 1
!        if (cnt.gt.globrdmatts) then
!          write(ilog,67) imol2,rsmol(imol2,1),amino(seqtyp(rsmol(imol2,1))),rsmol(imol2,2),amino(seqtyp(rsmol(imol2,2)))
!          write(ilog,*) 'This is related to a problematic intermolecular crosslink.'
!          if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!            write(ilog,77)
!            if (be_unsafe.EQV..false.) call fexit()
!          end if
!          exit
!        end if
!        if (crosslink(lkcnt)%itstype.eq.2) then
!          intsvec(1) = random()*360.0 - 180.0
!          intsvec(4) = random()*360.0 - 180.0
!          intsvec(7) = random()*360.0 - 180.0
!        else
!          write(ilog,*) 'Fatal. Encountered unsupported crosslink type in randomize_bb(...).&
!   & This is most likely an omission bug.'
!          call fexit()
!        end if
!        call crosslink_follow(lkcnt,imol2,jmol,intsvec)
!!       make sure intermolecular energies (only backward) are sane
!        et(1) = 0.0
!        et(12) = 0.0
!        iacnt = 0
!        do k=rsmol(jmol,1),rsmol(jmol,2)
!          call e_boundary_rs(k,et,azero)
!          do j=1,rsmol(imol,2)
!            call dis_bound(refat(j),refat(k),dis2) 
!            dis = sqrt(dis2)
!            if (dis.lt.(mcnb_cutoff+resrad(k)+resrad(j))) then
!              call Ven_rsp(et,j,k,svte,atrue)
!              iacnt = iacnt + 1
!            end if
!          end do
!          call Ven_rsp(et,min(k,rs2),max(k,rs2),svte,atrue)
!        end do
!        etval = (et(1)+et(12))/(1.0*iacnt)
!      end do
!!
!    else ! this is the normal part for just randomizing a molecule
!!
!      if (((globrandomize.eq.0).OR.(globrandomize.eq.3).OR.(globrandomize.eq.4)).AND.(imol.le.size(pdb_didread,dim=1))) cycle
!      cnt = 0
!      etval = HUGE(etval)
!      do while (etval.gt.globrdmthresh)
!        cnt = cnt + 1
!!       obviously we got into a mess: clear-up
!        if (cnt.gt.globrdmatts) then
!          write(ilog,66) imol,rsmol(imol,1),amino(seqtyp(rsmol(imol,1))),rsmol(imol,2),amino(seqtyp(rsmol(imol,2)))
!          write(ilog,*) 'This is a density or size problem.'
!          if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.2).OR.(dyn_mode.eq.3).OR.(dyn_mode.eq.6))) then
!            write(ilog,77)
!            if (be_unsafe.EQV..false.) call fexit()
!          end if

!          write(ilog,*) 'WARNING: Could not arrange molecules &
!   &with non-dramatic overlap (for molecule ',imol,').'
!          exit
!        end if
!!       randomize pos for molecule imol
!        call randomize_trans(imol,dum3)
!        do i=atmol(imol,1),atmol(imol,2)
!          x(i) = x(i) + dum3(1)
!          y(i) = y(i) + dum3(2)
!          z(i) = z(i) + dum3(3)
!        end do
!        call randomize_rot(imol)
!!       make sure intermolecular energies (only backward) are sane
!        et(1) = 0.0
!        et(12) = 0.0
!        iacnt = 0
!        do i=rsmol(imol,1),rsmol(imol,2)
!          call e_boundary_rs(i,et,azero)
!          do j=1,rsmol(imol,1)-1
!            if (j.lt.rsmol(imol,1)) then
!              call dis_bound(refat(i),refat(j),dis2) 
!              dis = sqrt(dis2)
!              if (dis.lt.(mcnb_cutoff+resrad(i)+resrad(j))) then
!                iacnt = iacnt + 1
!                call Ven_rsp(et,j,i,svte,atrue)
!              end if
!            end if
!          end do
!        end do
!        etval = (et(1)+et(12))/(1.0*iacnt)
!      end do
!    end if
!    call update_rigid(imol)
!  end do
!!
!  do imol=1,size(pdb_didread,dim=1)
!    if ((globrandomize.eq.0).OR.(globrandomize.eq.2).OR.(globrandomize.eq.3)) exit ! should never trigger with globrandomize 1
!    write(*,*) imol,pdb_didread(imol,1:2)
!    if (pdb_didread(imol,1).gt.rsmol(imol,1)) then
!      call MC_randomize(imol,rsmol(imol,1),pdb_didread(imol,1)-1,.true.,.true.,.false.)
!      write(ilog,*) 'Was doing N-terminal segment ',rsmol(imol,1),' to ',pdb_didread(imol,1)-1,' now.'
!    end if
!    if (pdb_didread(imol,2).lt.rsmol(imol,2)) then
!      call MC_randomize(imol,pdb_didread(imol,2)+1,rsmol(imol,2),.false.,.true.,.false.)
!      write(ilog,*) 'Was doing C-terminal segment ',pdb_didread(imol,2)+1,' to ',rsmol(imol,2),' now.'
!    end if
!  end do
!!  call fexit()

!! finalize a bit
!  use_IPP = log_0
!  use_attLJ = log_1
!  use_CORR = log_2
!  use_IMPSOLV = log_3
!  use_WCA = log_4
!  use_POLAR = log_5
!  use_ZSEC = log_6
!  use_TOR = log_7
!  use_FEG = log_8
!  scale_IPP = scale_1
!  use_hardsphere = log_hs
!  mcnb_cutoff = cut
!!
!end
!
!-----------------------------------------------------------------
!
subroutine stretch_randomize(imol,rsi,rsf,isbad)
!
  use aminos
  use iounit
  use energies
  use sequen
  use molecule
  use polypep
  use fyoc
  use cutoffs
  use zmatrix
  use atoms, ONLY: svte
!
  implicit none
!
  RTYPE et(MAXENERGYTERMS),dis,dis2,etbu,etpen,eto,etval,etbu2
  RTYPE, ALLOCATABLE:: curvals(:,:),potstr(:,:)
  integer imol,rsi,rsf,i,j,npots,cnt,azero,aone,atwo,athree,dofcnt,dofcntbu,dofcntbu2,alcsz,phs1,phs2,iacnt
  logical afalse,atrue,isbad
!
  isbad = .false.
  if (rsi.gt.rsf) return
!
  alcsz = at(rsf)%bb(1)-at(rsi)%bb(1)+at(rsf)%nbb+at(rsf)%nsc
  allocate(curvals(alcsz,2))
  afalse = .false.
  atrue = .true.
  azero = 0
  aone = 1
  atwo = 2
  athree = 3
  dofcnt = 0
  dofcntbu2 = 0
  phs1 = 1*globrdmatts/3 ! at least 1
  phs2 = 2*globrdmatts/3 ! at least 2
!
! setup up an extra potential if this is a randomization within an intramolecular crosslink loop
  if (n_crosslinks.gt.0) allocate(potstr(n_crosslinks,4))
  npots = 0
  do i=1,n_crosslinks
    if ((molofrs(crosslink(i)%rsnrs(1)).eq.imol).AND.&
 &      (molofrs(crosslink(i)%rsnrs(2)).eq.imol)) then
      if ((rsi.ge.minval(crosslink(i)%rsnrs(1:2))).AND.(rsf.le.maxval(crosslink(i)%rsnrs(1:2))).AND.&
 &        (crosslink(i)%itstype.eq.1)) then
        npots = npots + 1
        potstr(npots,1) = 1.0*cai(rsf) ! cai(crosslink(i)%rsnrs(1))
        potstr(npots,2) = 1.0*cai(crosslink(i)%rsnrs(2))
        potstr(npots,3) = 3.0 + 2.0*(maxval(crosslink(i)%rsnrs(1:2))-rsf) ! max dis before penalty sets in
        potstr(npots,4) = 1000 ! 0.1*globrdmthresh ! kcal/molA^2
      end if
    end if
  end do
! 
  do i=rsi,rsf
    et(18) = 0.
    if (use_BOND(4).EQV..true.) call en_torsions(i,et)
    etpen = et(18)
    et(1) = HUGE(etbu)
    etbu2 = et(1)
    etbu = et(1)
    cnt = 0
    if (i.gt.rsi) dofcntbu2 = dofcntbu
    dofcntbu = dofcnt
    call setrandom_forres(i,alcsz,curvals,dofcnt,azero,afalse)
    dofcnt = dofcntbu
    do while (etbu2.gt.globrdmthresh)
      cnt = cnt + 1
      if ((cnt.eq.(phs2+1)).AND.(i.gt.rsi)) then
        et(18) = 0.
        if (use_BOND(4).EQV..true.) call en_torsions(i-1,et)
        etpen = etpen + et(18)
      end if
!     obviously we got into a mess: clear-up
      if (cnt.gt.globrdmatts) exit
!
      iacnt = 0
      dofcnt = dofcntbu
      et(18) = 0.
      if ((cnt.gt.phs2).AND.(i.gt.rsi)) then
!       randomize bb-angles for residue i-1 (this becomes expensive quickly)
        dofcnt = dofcntbu2
        call setrandom_forres(i-1,alcsz,curvals,dofcnt,athree,atrue)
        if (dofcnt.ne.dofcntbu) then
          write(ilog,*) 'Fatal. D.o.f counter is corrupt in stretch_randomize(...). This is a bug.'
          call fexit()
        end if
        if (use_BOND(4).EQV..true.) call en_torsions(i-1,et)
      else if (cnt.gt.phs1) then
!       add side chains if still failing
        call setrandom_forres(i,alcsz,curvals,dofcnt,atwo,atrue)
        dofcnt = dofcntbu
      end if
!     randomize bb-angles for residue i
      call setrandom_forres(i,alcsz,curvals,dofcnt,aone,atrue)
      if (use_BOND(4).EQV..true.) call en_torsions(i,et)
      eto = et(18)
!     make sure backward-in-chain-energies are sane
      et(1) = 0.0
      if (use_IPP.EQV..true.) then
        do j=rsmol(imol,1),min(i+1,rsf)
          call dis_bound(refat(i),refat(j),dis2)
          dis = sqrt(dis2)
          if (dis.lt.(mcnb_cutoff+resrad(i)+resrad(j))) then
            call Ven_rsp(et,min(j,i),max(j,i),svte,atrue)
            iacnt = iacnt + 1
          end if
          if ((cnt.gt.phs2).AND.(j.ne.i).AND.(i.gt.rsi)) then
            call dis_bound(refat(i-1),refat(j),dis2) 
            dis = sqrt(dis2)
            if (dis.lt.(mcnb_cutoff+resrad(i-1)+resrad(j))) then
              call Ven_rsp(et,min(j,i-1),max(j,i-1),svte,atrue)
              iacnt = iacnt + 1
            end if
          end if
        end do
      else
        iacnt = 1
      end if
      do j=1,npots
        iacnt = iacnt + 1
        call dis_bound(nint(potstr(j,1)),nint(potstr(j,2)),dis2)
        dis = sqrt(dis2)
        if (dis.gt.potstr(j,3)) then
          et(1) = et(1) + potstr(j,4)*(dis-potstr(j,3))*(dis-potstr(j,3))
        end if
      end do
      etval = (et(1) + eto - etpen) ! /(1.0*iacnt)
      if ((etval.lt.etbu).AND.((dofcnt.le.0).OR.(dofcnt.gt.dofcntbu).OR.((cnt.gt.phs2).AND.(i.gt.rsi).AND.(dofcnt.gt.dofcntbu2))))&
 & then
        if ((cnt.gt.phs2).AND.(i.gt.rsi)) then
          curvals((dofcntbu2+1):dofcnt,2) = curvals((dofcntbu2+1):dofcnt,1)
        else
          curvals((dofcntbu+1):dofcnt,2) = curvals((dofcntbu+1):dofcnt,1)
        end if
        etbu = etval
        etbu2 = etbu/(1.0*iacnt)
      end if
      if ((dofcnt.eq.dofcntbu).AND.(cnt.lt.phs1)) then
!       turn on side chain addition right away
        cnt = phs1
      else if ((dofcnt.eq.dofcntbu).AND.(cnt.lt.phs2)) then
!       turn on panic mode right away
        cnt = phs2
      end if
      if (cnt.eq.globrdmatts) then ! still here means trouble - restore best guess
        isbad = .true.
        if (i.gt.rsi) then
          curvals((dofcntbu2+1):dofcnt,1) = curvals((dofcntbu2+1):dofcnt,2)
          dofcnt = dofcntbu2
          call setrandom_forres(i-1,alcsz,curvals,dofcnt,athree,afalse)
          call setrandom_forres(i,alcsz,curvals,dofcnt,athree,afalse)
        else
          curvals((dofcntbu+1):dofcnt,1) = curvals((dofcntbu+1):dofcnt,2)
          dofcnt = dofcntbu
          call setrandom_forres(i,alcsz,curvals,dofcnt,athree,afalse)
        end if
      end if
    end do
  end do
  deallocate(curvals)
  if (n_crosslinks.gt.0) deallocate(potstr)
!
end
!
!----------------------------------------------------------------------------------
!
! a small helper for stretch_randomize
!
subroutine setrandom_forres(i,alcsz,curvals,curcur,mode,dosam)
!
  use fyoc
  use polypep
  use molecule
  use atoms
  use zmatrix
  use sequen
!
  implicit none
!
  integer alcsz,curcur,mode,j,ati,i,azero
  RTYPE curvals(alcsz,2),random
  logical isnat,dosam
!
  azero = 0
!
  if ((mode.eq.1).OR.(mode.eq.3)) then
    if (fline(i).gt.0) then
      if (izrot(fline(i))%alsz.gt.0) then
        if (seqflag(i).ne.5) then
          if (dosam.EQV..true.) curvals(curcur+1,1) = random()*360.0 - 180.0
        else
          if (dosam.EQV..true.) curvals(curcur+1,1) = phi(i)
        end if
      else
        if (dosam.EQV..true.) curvals(curcur+1,1) = phi(i)
      end if
    end if
    if (yline(i).gt.0) then
      if (izrot(yline(i))%alsz.gt.0) then
        if (dosam.EQV..true.) curvals(curcur+2,1) = random()*360.0 - 180.0
      else
        if (dosam.EQV..true.) curvals(curcur+2,1) = psi(i)
      end if
    end if
  else if (mode.eq.0) then
    if (fline(i).gt.0) then
      curvals(curcur+1,1) = phi(i)
    end if
    if (yline(i).gt.0) then
      curvals(curcur+2,1) = psi(i)
    end if
  end if
  if ((fline(i).gt.0).OR.(yline(i).gt.0)) then
    if ((mode.eq.1).OR.(mode.eq.3)) call setfy(i,curvals(curcur+1,1),curvals(curcur+2,1),azero)
    curcur = curcur + 2
  end if
  do ati=at(i)%bb(1),at(i)%bb(1)+at(i)%nbb+at(i)%nsc-1
    if (izrot(ati)%alsz.gt.0) then
      if ((ati.eq.fline(i)).OR.(ati.eq.yline(i)).OR.(ati.eq.wline(i))) cycle
      isnat = .false.
      if ((ati.eq.nucsline(1,i)).OR.(ati.eq.nucsline(2,i)).OR.(ati.eq.nucsline(3,i)).OR.&
 &        (ati.eq.nucsline(4,i))) isnat = .true.
      if (i.gt.1) then
        if ((ati.eq.nucsline(5,i-1)).OR.(ati.eq.nucsline(3,i-1))) isnat = .true.
      end if
      if ((izrot(ati)%rotis(izrot(ati)%alsz,2).gt.(at(i)%bb(1)+at(i)%nbb+at(i)%nsc-1)).OR.(isnat.EQV..true.)) then
        curcur = curcur + 1
        if ((mode.eq.1).OR.(mode.eq.3)) then
          if (dosam.EQV..true.) curvals(curcur,1) = random()*360.0 - 180.0
          call setother(ati,curvals(curcur,1),azero,isnat)
        else if (mode.eq.0) then
          curvals(curcur,1) = ztor(ati)
        end if
      end if
    end if
  end do
  do ati=at(i)%bb(1),at(i)%bb(1)+at(i)%nbb+at(i)%nsc-1
    if (izrot(ati)%alsz.gt.0) then
      if ((ati.eq.fline(i)).OR.(ati.eq.yline(i)).OR.(ati.eq.wline(i))) cycle
      if ((ati.eq.nucsline(1,i)).OR.(ati.eq.nucsline(2,i)).OR.(ati.eq.nucsline(3,i)).OR.&
 &        (ati.eq.nucsline(4,i))) cycle
      if (i.gt.1) then
        if ((ati.eq.nucsline(5,i-1)).OR.(ati.eq.nucsline(3,i-1))) cycle
      end if
      isnat = .false.
      do j=1,nchi(i)
        if (ati.eq.chiline(j,i)) isnat = .true.
      end do
      if ((seqtyp(i).ne.26).AND.(isnat.EQV..false.)) cycle  ! don't perturb those in unslst
      if (atmres(izrot(ati)%rotis(izrot(ati)%alsz,2)).eq.atmres(ati)) then
        curcur = curcur + 1
        if ((mode.eq.2).OR.(mode.eq.3)) then
          if (dosam.EQV..true.) curvals(curcur,1) = random()*360.0 - 180.0
          call setother(ati,curvals(curcur,1),azero,isnat)
        else if (mode.eq.0) then
          curvals(curcur,1) = ztor(ati)
        end if
      end if
    end if
  end do
  if (mode.gt.0) call makexyz_forbb(i)
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! evaluate and possibly store current configuration for entire stretch at FYOC level
! 
subroutine eval_forstretch(rsi,rsf,vals,valsz,nvals,clki)
!
  use cutoffs
  use energies
  use sequen
  use atoms, ONLY: svte,atmres
  use fyoc
  use zmatrix
  use molecule
  use movesets, ONLY: unklst
!
  implicit none
!
  integer, INTENT(IN):: rsi,rsf,clki,valsz
  integer, INTENT(OUT):: nvals
  RTYPE, INTENT(INOUT):: vals(valsz)
!
  RTYPE et(MAXENERGYTERMS),dis,dis2,etval
  integer i,j,imol,iacnt
  logical atrue
!
  atrue = .true.
!
  et(:) = 0.0
  iacnt = 0
  imol = molofrs(rsi)
  do i=rsi,rsf
    if (use_BOND(4).EQV..true.) call en_torsions(i,et)
    if (use_IPP.EQV..true.) then
      do j=rsmol(imol,1),i
        call dis_bound(refat(i),refat(j),dis2)
        dis = sqrt(dis2)
        if (dis.lt.(mcnb_cutoff+resrad(i)+resrad(j))) then
          iacnt = iacnt + 1
          call Ven_rsp(et,min(j,i),max(j,i),svte,atrue)
        end if
      end do
    end if
  end do
  if (use_IPP.EQV..false.) iacnt = rsf-rsi+1
  etval = sum(et) ! /(1.0*iacnt)
  if (etval.lt.vals(1)) then
    vals(1) = etval
!
    nvals = 1
    do i=rsi,rsf
      if (wline(i).gt.0) then
        nvals = nvals + 1
        vals(nvals) = omega(i)
      end if
      if ((fline(i).gt.0).OR.(fline2(i).gt.0)) then
        nvals = nvals + 1
        vals(nvals) = phi(i)
      end if
      if ((yline(i).gt.0).OR.(yline2(i).gt.0)) then
        nvals = nvals + 1
        vals(nvals) = psi(i)
      end if
      do j=1,nchi(i)
        nvals = nvals + 1
        vals(nvals) = chi(j,i)
      end do
      do j=1,nnucs(i)
        nvals = nvals + 1
        vals(nvals) = nucs(j,i)
      end do
    end do
    do j=1,unklst%nr
      i = atmres(unklst%idx(j))
      if ((i.lt.rsi).OR.(i.gt.rsf)) cycle
      nvals = nvals + 1
      vals(nvals) = ztor(unklst%idx(j))
    end do
    vals(valsz) = vals(1)/(1.0*iacnt)
  end if
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! restore the best configuration for entire stretch at FYOC level
! 
subroutine apply_forstretch(rsi,rsf,vals,valsz)
!
  use sequen
  use atoms, ONLY: svte,atmres
  use fyoc
  use zmatrix
  use molecule
  use movesets, ONLY: unklst
!
  implicit none
!
  integer, INTENT(IN):: rsi,rsf,valsz
  RTYPE, INTENT(IN):: vals(valsz)
!
  integer i,j,nvals
!
  nvals = 1
  do i=rsi,rsf
    if (wline(i).gt.0) then
      nvals = nvals + 1
      omega(i) = vals(nvals)
      ztor(wline(i)) = vals(nvals)
    end if
    if ((fline(i).gt.0).OR.(fline2(i).gt.0)) then
      nvals = nvals + 1
      phi(i) = vals(nvals)
      if (fline(i).gt.0) ztor(fline(i)) = vals(nvals)
      if (fline2(i).gt.0) then
        ztor(fline2(i)) = ztor(fline(i)) + phish(i)
        if (ztor(fline2(i)).gt.180.0) ztor(fline2(i)) = ztor(fline2(i)) - 360.0
        if (ztor(fline2(i)).lt.-180.0) ztor(fline2(i)) = ztor(fline2(i)) + 360.0
      end if
    end if
    if ((yline(i).gt.0).OR.(yline2(i).gt.0)) then
      nvals = nvals + 1
      psi(i) = vals(nvals)
      if (yline(i).gt.0) ztor(yline(i)) = vals(nvals)
      if (yline2(i).gt.0) then
        ztor(yline2(i)) = ztor(yline(i)) + psish(i)
        if (ztor(yline2(i)).gt.180.0) ztor(yline2(i)) = ztor(yline2(i)) - 360.0
        if (ztor(yline2(i)).lt.-180.0) ztor(yline2(i)) = ztor(yline2(i)) + 360.0
      end if
      if (yline(i).gt.0) then
        if (ztor(yline(i)).lt.0.0) then
          ztor(yline(i)) = ztor(yline(i)) + 180.0
        else
          ztor(yline(i)) = ztor(yline(i)) - 180.0
        end if
      end if
    end if
    do j=1,nchi(i)
      nvals = nvals + 1
      chi(j,i) = vals(nvals)
      ztor(chiline(j,i)) = vals(nvals)
    end do
    do j=1,nnucs(i)
      nvals = nvals + 1
      nucs(j,i) = vals(nvals)
      ztor(nucsline(j,i)) = vals(nvals)
    end do
  end do
  do j=1,unklst%nr
    i = atmres(unklst%idx(j))
    if ((i.lt.rsi).OR.(i.gt.rsf)) cycle
    nvals = nvals + 1
    ztor(unklst%idx(j)) = vals(nvals)
  end do
  call makexyz_forbb(rsi)
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! evaluate and possibly store current configuration for rigid-body displacements
!
subroutine evaluate_formol(imol,mmol,invals,vals,molmap,is_clk)
!
  use cutoffs
  use energies
  use sequen
  use atoms, ONLY: svte,atmres,x,y,z
  use fyoc
  use zmatrix
  use molecule
  use movesets, ONLY: unklst
  use polypep
!
  implicit none
!
  integer, INTENT(IN):: imol,mmol,molmap(nmol)
  RTYPE, INTENT(IN):: invals(9)
  RTYPE, INTENT(INOUT):: vals(11)
  logical, INTENT(IN):: is_clk
!
  RTYPE et(MAXENERGYTERMS),dis,dis2,etval
  integer i,j,iacnt,lo,hi,hx,azero
  logical atrue
!
  atrue = .true.
  azero = 0
!
  et(:) = 0.0
  iacnt = 1 ! minimum
  do i=rsmol(imol,1),rsmol(imol,2)
    if (disulf(i).gt.0) then
      if (molmap(molofrs(disulf(i))).gt.0) then
        if (use_BOND(4).EQV..true.)  call en_torsions(i,et)
        if (use_BOND(1).EQV..true.)  call en_bonds(i,et)
        if (use_BOND(2).EQV..true.)  call en_angles(i,et)
        if (use_BOND(3).EQV..true.)  call en_impropers(i,et)
        if (use_BOND(5).EQV..true.)  call en_cmap(i,et)
      end if
    end if
    call e_boundary_rs(i,et,azero)
    if (mmol.gt.0) then
      do j=1,rsmol(mmol,2)
        if (molmap(molofrs(j)).le.0) cycle
        if (disulf(j).gt.0) then
          if ((disulf(j).eq.i).AND.(molofrs(j).ne.imol)) then
            if (use_BOND(4).EQV..true.) call en_torsions(j,et)
            if (use_BOND(1).EQV..true.)  call en_bonds(j,et)
            if (use_BOND(2).EQV..true.)  call en_angles(j,et)
            if (use_BOND(3).EQV..true.)  call en_impropers(j,et)
            if (use_BOND(5).EQV..true.)  call en_cmap(j,et)
            if (use_IPP.EQV..true.) then
              call Ven_rsp(et,min(j,disulf(j)),max(j,disulf(j)),svte,atrue)
              iacnt = iacnt + 1
            end if
          end if
        end if
        if (use_IPP.EQV..true.) then
          call dis_bound(refat(j),refat(i),dis2)
          dis = sqrt(dis2)
          if (dis.lt.(mcnb_cutoff+resrad(i)+resrad(j))) then
            call Ven_rsp(et,min(i,j),max(i,j),svte,atrue)
            iacnt = iacnt + 1
          end if
        end if
      end do
    end if
  end do
  if (use_IPP.EQV..false.) iacnt = rsmol(imol,2)-rsmol(imol,1)+1
  etval = sum(et)
!
  if (etval.lt.vals(10)) then
    vals(10) = etval
    if (is_clk.EQV..false.) then
      lo = atmol(imol,1)
      hi = min(atmol(imol,2),atmol(imol,1)+2)
      hx = hi-lo+1
      vals(1:hx) = x(lo:hi)
      vals(4:(hx+3)) = y(lo:hi)
      vals(7:(hx+6)) = z(lo:hi)
    else
      vals(1:9) = invals(1:9)
    end if
    vals(11) = vals(10)/(1.0*iacnt)
  end if
!
end
!
!--------------------------------------------------------------------------------------------------------
!
subroutine apply_formol(imol,vals)
!
  use molecule, ONLY: atmol
  use atoms, ONLY: x,y,z
!
  implicit none
!
  integer, INTENT(IN):: imol
  RTYPE, INTENT(IN):: vals(9)
!
  integer lo,hi,hx
!
  lo = atmol(imol,1)
  hi = min(atmol(imol,2),atmol(imol,1)+2)
  hx = hi-lo+1
  x(lo:hi) = vals(1:hx)
  y(lo:hi) = vals(4:(hx+3))
  z(lo:hi) =  vals(7:(hx+6))
!
  call makexyz_formol(imol)
!
end
!
!--------------------------------------------------------------------------------------------------------
!
! a routine to use something resembling a 0K Monte Carlo to relax a missing tail
! the energy criterion for convergence here is tighter (due to normalization) than in
! the other functions, and the whole tail is sampled
!
subroutine MC_randomize(imol,rsi,rsf,align_C,fone,singlemol,isbad)
!
  use atoms
  use zmatrix
  use polypep
  use fyoc
  use sequen
  use cutoffs
  use energies
  use math
  use molecule
  use pdb, ONLY: pdb_didread
!
  implicit none
!
  integer, INTENT(IN):: rsi,rsf,imol
  logical, INTENT(IN):: align_C,fone,singlemol
!
  integer rs,i,sta,sto,inc,azero,aone,afour,dofs,atl,ath,cnt,curl,pdbrsz
  RTYPE evec(MAXENERGYTERMS),random,picky,etbu,ettmp,oth,doth,etvs(2)
  logical needfy,needone,needC,noopt,isbad
!
  azero = 0
  afour = 4
  aone = 1
  isbad = .false.
!
  pdbrsz = 0
  if (allocated(pdb_didread).EQV..true.) then
    pdbrsz = size(pdb_didread,dim=1)
  end if
  if (align_C.EQV..true.) then
    sta = rsf
    sto = rsi
    inc = -1
  else
    sta = rsi
    sto = rsf 
    inc = 1
  end if
!
  rs = sta
  do while ((rs.ge.rsi).AND.(rs.le.rsf))
    cnt = 0
    atl = at(rsi)%bb(1)
    ath = at(rsf)%bb(1)+at(rsf)%nsc+at(rsf)%nbb-1
    
    dofs = 0
    needfy = .false.
    noopt = .false.
    do i=at(rs)%bb(1)+1,at(rs)%bb(1)+at(rs)%nsc+at(rs)%nbb-1
      if (izrot(i)%alsz.gt.0) then
        if ((fline(rs).eq.i).OR.(fline2(rs).eq.i)) then
          dofs = dofs + 1
          needfy = .true.
        else if ((yline(rs).eq.i).OR.(yline2(rs).eq.i)) then
          dofs = dofs + 1
          needfy = .true.
        else if (wline(rs).ne.i) then ! exclude omega
          dofs = dofs + 1
        end if
      end if
    end do
    if (dofs.le.0) then
      noopt = .true.
      picky = 0.
    else
      picky = min(1.0,2.0/(1.0*dofs))
    end if
    ztorpr(atl:ath) = ztor(atl:ath)
    xref(atl:ath) = x(atl:ath)
    yref(atl:ath) = y(atl:ath)
    zref(atl:ath) = z(atl:ath)
    evec(:) = 0.
    call randomize_energy_short(rs,evec,use_cutoffs,align_C,singlemol,pdbrsz,pdb_didread,curl,aone) 
    etbu = sum(evec)
    etvs(1) = etbu/curl
!
    needone = fone
    do while ((etvs(1).gt.globrdmthresh).OR.(needone.EQV..true.))
      cnt = cnt + 1
      if (cnt.gt.globrdmatts) exit
      dofs = 0
      if (needfy.EQV..true.) then  
        if (random().le.picky) then
          cur_phi = random()*360.0 - 180.0
          cur_psi = random()*360.0 - 180.0
          call sample_bb(rs,aone) ! may not use both
          call quatrot_pivot(aone,rs,azero)
          if (align_c.eqv..true.) call quatrot_pivot(azero,rs,azero)
          dofs = dofs + 2
        end if
      end if
      do i=at(rs)%bb(1)+1,at(rs)%bb(1)+at(rs)%nbb+at(rs)%nsc-1
        if ((i.eq.fline(rs)).OR.(i.eq.yline(rs)).OR.(i.eq.fline2(rs)).OR.(i.eq.yline2(rs))) cycle
        if (rs.gt.rsmol(imol,1)) then
          if (i.eq.yline2(rs-1)) cycle
        end if
        if (izrot(i)%alsz.gt.0) then
          if (i.eq.atmol(imol,1)) cycle
          if (i.eq.wline(rs)) cycle
          if (random().le.picky) then
            if (maxval(izrot(i)%rotis(:,2)).gt.ath) then
              needC = .true.
            else
              needC = .false.
            end if
            oth = random()*360.0 - 180.0
            doth = oth - ztor(i)
            if (doth.gt.180.0) doth = doth - 360.0
            if (doth.lt.-180.0) doth = doth + 360.0
            ztor(i) = oth
            doth = doth/RADIAN
            call quatxyz_forrotlst(i,needC,doth,azero)
            dofs = dofs + 1
          end if
        end if
      end do
      if ((dofs.le.0).AND.(noopt.EQV..false.)) cycle 
      evec(:) = 0.
      call randomize_energy_short(rs,evec,use_cutoffs,align_C,singlemol,pdbrsz,pdb_didread,curl,aone) 
      ettmp = sum(evec)
!
      if (ettmp.lt.etbu) then
        etbu = ettmp
        etvs(1) = etbu/curl
        ztorpr(atl:ath) = ztor(atl:ath)
        xref(atl:ath) = x(atl:ath)
        yref(atl:ath) = y(atl:ath)
        zref(atl:ath) = z(atl:ath)
        needone = .false.
      else
        ztor(atl:ath) = ztorpr(atl:ath)
        x(atl:ath) = xref(atl:ath)
        y(atl:ath) = yref(atl:ath)
        z(atl:ath) = zref(atl:ath)
      end if
      call zmatfyc_rs(rs)
      if (noopt.EQV..true.) exit ! single attempt
    end do
    rs = rs + inc
  end do
!
  evec(:) = 0.
  call randomize_energy_short(sta,evec,use_cutoffs,align_C,singlemol,pdbrsz,pdb_didread,curl,aone) 
  etbu = sum(evec)/curl
  if (etbu.gt.globrdmthresh) isbad = .true.
!
end
!
!-------------------------------------------------------------------------------------------------------------
!
! a similar routine to do just side chain sampling on a single residue
!
subroutine MC_randomize_sc(imol,rs,singlemol,isbad)
!
  use atoms
  use zmatrix
  use polypep
  use fyoc
  use sequen
  use cutoffs
  use energies
  use math
  use molecule
  use pdb, ONLY: pdb_didread,pdbmap2
!
  implicit none
!
  integer, INTENT(IN):: rs,imol
  logical, INTENT(IN):: singlemol
!
  integer i,azero,atwo,dofs,atl,ath,cnt,curl,pdbrsz,cdf,doflst(at(rs)%nbb+at(rs)%nsc),k,ki,kj,dref
  RTYPE evec(MAXENERGYTERMS),random,picky,etbu,ettmp,oth,doth,etvs(2),tm
  logical isbad,afalse,isvalid
!
  azero = 0
  atwo = 2
  isbad = .false.
  afalse = .false.
!
  pdbrsz = 0
  if (allocated(pdb_didread).EQV..true.) then
    pdbrsz = size(pdb_didread,dim=1)
  end if
!
  cnt = 0
  atl = at(rs)%bb(1)
  ath = at(rs)%bb(1)+at(rs)%nsc+at(rs)%nbb-1
  
  dofs = 0
  do i=atl,ath
    if (izrot(i)%alsz.gt.0) then
      if ((maxval(izrot(i)%rotis).le.ath).AND.(minval(izrot(i)%rotis).ge.atl)) then
        isvalid = .true.
        tm = 0.0
        dref = 0
        do k=1,izrot(i)%alsz
          ki = izrot(i)%rotis(k,1)
          kj = izrot(i)%rotis(k,2)
          if (minval(pdbmap2(ki:kj)).gt.0) then
            isvalid = .false.
            exit
          end if
          tm = tm + sum(mass(ki:kj))
          dref = max(dref,maxval(n12(ki:kj)))
        end do
        if ((dref.le.1).AND.(tm.le.15.0)) isvalid = .false.
        if (isvalid.EQV..true.) then
          dofs = dofs + 1
          doflst(dofs) = i
        end if
      end if
    end if
  end do
  if (dofs.le.0) then
    picky = 0.
  else
    picky = min(1.0,2.0/(1.0*dofs))
  end if
!
  ztorpr(atl:ath) = ztor(atl:ath)
  xref(atl:ath) = x(atl:ath)
  yref(atl:ath) = y(atl:ath)
  zref(atl:ath) = z(atl:ath)
  evec(:) = 0.
  call randomize_energy_short(rs,evec,use_cutoffs,afalse,singlemol,pdbrsz,pdb_didread,curl,atwo)
  etbu = sum(evec)
  etvs(1) = etbu/curl
!
  do while ((etvs(1).gt.globrdmthresh).AND.(dofs.gt.0))
    cnt = cnt + 1
    if (cnt.gt.globrdmatts) exit
    cdf = 0
    do ki=1,dofs
      i = doflst(ki)
      if (random().le.picky) then
        oth = random()*360.0 - 180.0
        doth = oth - ztor(i)
        if (doth.gt.180.0) doth = doth - 360.0
        if (doth.lt.-180.0) doth = doth + 360.0
        ztor(i) = oth
        doth = doth/RADIAN
        call quatxyz_forrotlst(i,afalse,doth,azero)
        cdf = cdf + 1
      end if
    end do
    if (cdf.le.0) cycle
    evec(:) = 0.
    call randomize_energy_short(rs,evec,use_cutoffs,afalse,singlemol,pdbrsz,pdb_didread,curl,atwo) 
    ettmp = sum(evec)
!
    if (ettmp.lt.etbu) then
      etbu = ettmp
      etvs(1) = etbu/curl
      ztorpr(atl:ath) = ztor(atl:ath)
      xref(atl:ath) = x(atl:ath)
      yref(atl:ath) = y(atl:ath)
      zref(atl:ath) = z(atl:ath)
    else
      ztor(atl:ath) = ztorpr(atl:ath)
      x(atl:ath) = xref(atl:ath)
      y(atl:ath) = yref(atl:ath)
      z(atl:ath) = zref(atl:ath)
    end if
    call zmatfyc_rs(rs)
  end do
!
  evec(:) = 0.
  call randomize_energy_short(rs,evec,use_cutoffs,afalse,singlemol,pdbrsz,pdb_didread,curl,atwo) 
  etbu = sum(evec)/curl
  if (etbu.gt.globrdmthresh) isbad = .true.
!
end
!
!-------------------------------------------------------------------------------------------------------------
!
