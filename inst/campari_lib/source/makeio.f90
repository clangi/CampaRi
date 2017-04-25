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
! MAIN AUTHOR:   Rohit Pappu                                               !
! CONTRIBUTIONS: Andreas Vitalis, Hoang Tran                               !
!                                                                          !
!--------------------------------------------------------------------------!
!
!
! #############################################################
! ##                                                         ##
! ## subroutine makeio -- open or close desired output files ##
! ##                                                         ##
! #############################################################
!
! "makeio" opens or closes output files for Monte Carlo simulations
! mode = 1: open
! mode = 2: or anything else close
!
! for parallel computing filehandles are a mess, here we open all files in the
! same (running) directory (REMC) or we suppress output (averaging)
! WARNING:
! this means that ALL GENERAL CHANGES HAVE TO BE MADE AT LEAST TWICE 
!
subroutine makeio(mode)
!
  use iounit
  use mcsums
  use polyavg
  use system
  use molecule
  use torsn
  use movesets
  use mpistuff
  use energies
  use dssps
  use fos
  use wl
  use paircorr
  use clusters
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer freeunit,lext,mode,mt,imol
  logical exists,doperser
  character(3) mtstr
  character(60) dihedfile,polmrfile
  character(60) holesfile,accfile,trcvfile
  character(60) persfile,sasafile,enfile,dsspfile
  character(60) phtfile,ensfile,elsfile,pcfile,rmsdfile
!  character(60) bprtnrfile
  logical fycopen,polymopen,ensopen,dsspopen,&
 &accopen,enopen,persopen,holesopen,sasaopen,phtopen,elsopen,pcopen,rmsdopen!,bpopen
#ifdef ENABLE_MPI
  integer masterrank,t1,t2,iomess
  character(3) xpont
  character(60) remcfile,retrfile
  logical remcopen,retropen
#endif
#ifdef ENABLE_THREADS
  character(60) threadsfile
  logical threadsopen
#endif
!
! the save/data construct ensures that the logicals are treated like subroutine-globals
! including initialization (once(!)) 
  data fycopen/.false./,&
 &polymopen/.false./,&
 &accopen/.false./,dsspopen/.false./,&
 &enopen/.false./,persopen/.false./,ensopen/.false./,&
 &holesopen/.false./,sasaopen/.false./,phtopen/.false./,elsopen/.false./,pcopen/.false./,rmsdopen/.false./!,bpopen/.false./
#ifdef ENABLE_MPI
  data remcopen/.false./,retropen/.false./
#endif
#ifdef ENABLE_THREADS
  data threadsopen/.false./
#endif
!
  save fycopen,polymopen,ensopen,&
 &accopen,enopen,persopen,holesopen,sasaopen,&
 &phtopen,dsspopen,elsopen,pcopen,rmsdopen!,bpopen
#ifdef ENABLE_MPI
  save remcopen,retropen
#endif
#ifdef ENABLE_THREADS
  save threadsopen
#endif
!
  doperser = .false.
  do imol=1,nmol
    if (do_pers(moltypid(imol)).EQV..true.) then
      doperser = .true.
      exit
    end if
  end do
!
! for trajectory and helper file we block two units permanently (this is respected in freeunit in the same way as the default
! I/O channels for Fortran)
  itraj = freeunit()
  itrajhlp = freeunit()  
!
#ifdef ENABLE_MPI
!
  masterrank = 0
  if (use_REMC.EQV..true.) then
!   in this case we have to open the output files for each node individually
    if (mode .eq. 1) then
      lext = 3
      call int2str(myrank,xpont,lext)
      dihedfile  = 'N_'//xpont(1:lext)//'_FYC.dat'
      polmrfile  = 'N_'//xpont(1:lext)//'_POLYMER.dat'
      holesfile  = 'N_'//xpont(1:lext)//'_HOLES.dat'
      enfile  = 'N_'//xpont(1:lext)//'_ENERGY.dat'
      if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        accfile = 'N_'//xpont(1:lext)//'_ACCEPTANCE.dat'
      end if
      if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        ensfile = 'N_'//xpont(1:lext)//'_ENSEMBLE.dat'
      end if
      persfile = 'N_'//xpont(1:lext)//'_PERSISTENCE.dat'
      sasafile = 'N_'//xpont(1:lext)//'_SAV.dat'
      dsspfile = 'N_'//xpont(1:lext)//'_DSSP_RUNNING.dat'
      phtfile = 'N_'//xpont(1:lext)//'_PHTIT.dat'
      remcfile = 'N_'//xpont(1:lext)//'_EVEC.dat'
      retrfile = 'N_'//xpont(1:lext)//'_REXTRACE.dat'
      elsfile = 'N_'//xpont(1:lext)//'_ELS_WFRAMES.dat'
      pcfile = 'N_'//xpont(1:lext)//'_GENERAL_DIS.dat'
      rmsdfile = 'N_'//xpont(1:lext)//'_RMSD.dat'
!      bprtnrfile = 'N_'//xpont(1:lext)//'_BOUND_PARTNERS.dat'
#ifdef ENABLE_THREADS
      threadsfile = 'N_'//xpont(1:lext)//'_THREADS.log'
#endif
!  
      inquire(file=dihedfile,exist=exists)
      if (exists.EQV..true.) then
        idihed = freeunit()
        open (unit=idihed,file=dihedfile,status='old')
        close(unit=idihed,status='delete')
      end if
      if (torout.le.nsim) then
        idihed = freeunit()
        open (unit=idihed,file=dihedfile,status='new')
        fycopen = .true.
      end if
!   
! setup filehandle for polymeric outputs
      inquire(file=polmrfile,exist=exists)
      if (exists.EQV..true.) then
        ipolmr = freeunit()
        open (unit=ipolmr,file=polmrfile,status='old')
        close(unit=ipolmr,status='delete')
      end if
      if (polout .le. nsim) then
        ipolmr = freeunit()
        open (unit=ipolmr,file=polmrfile,status='new') 
        polymopen = .true.
      end if 
!   
! setup filehandle for holes outputs
      inquire(file=holesfile,exist=exists)
      if (exists.EQV..true.) then
        iholes = freeunit()
        open (unit=iholes,file=holesfile,status='old')
        close(unit=iholes,status='delete')
      end if
      if (holescalc.le.nsim) then
        iholes = freeunit()
        open (unit=iholes,file=holesfile,status='new') 
        holesopen = .true.
      end if
! 
! setup filehandle for pH titration outputs
      inquire(file=phtfile,exist=exists)
      if (exists.EQV..true.) then
        ipht = freeunit()
        open (unit=ipht,file=phtfile,status='old')
        close(unit=ipht,status='delete')
      end if
      if ((phfreq.gt.0.0).AND.(phout.le.nsim)) then
        ipht = freeunit()
        open (unit=ipht,file=phtfile,status='new') 
        phtopen = .true.
      end if 
!   
! setup filehandle for energy outputs
      inquire(file=enfile,exist=exists)
      if (exists.EQV..true.) then
        iene = freeunit()
        open (unit=iene,file=enfile,status='old')
        close(unit=iene,status='delete')
      end if
      if (enout.le.nsim) then
        iene = freeunit()
        open (unit=iene,file=enfile,status='new')     
        enopen = .true.
      end if
!
! setup filehandle for acceptance outputs
      if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
        inquire(file=accfile,exist=exists)
        if (exists.EQV..true.) then
          iacc = freeunit()
          open (unit=iacc,file=accfile,status='old')
          close(unit=iacc,status='delete')
        end if
        if (accout.le.nsim) then
          iacc = freeunit()
          open (unit=iacc,file=accfile,status='new')
          accopen = .true.
        end if
      end if
!
! setup filehandle for ensemble outputs
      if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
        inquire(file=ensfile,exist=exists)
        if (exists.EQV..true.) then
          iens = freeunit()
          open (unit=iens,file=ensfile,status='old')
          close(unit=iens,status='delete')
        end if
        if (ensout.le.nsim) then
          iens = freeunit()
          open (unit=iens,file=ensfile,status='new')
          ensopen = .true.
        end if
      end if
!
! setup filehandle for persistence length calculation output file
      inquire(file=persfile,exist=exists)
      if (exists.EQV..true.) then
        ipers = freeunit()
        open (unit=ipers,file=persfile,status='old')
        close(unit=ipers,status='delete')
      end if
      if ((polcalc.le.nsim).AND.(doperser.EQV..true.)) then
        ipers = freeunit()
        open (unit=ipers,file=persfile,status='new')
        persopen = .true.
      end if
!
! setup filehandle for on-the-fly SAV output
      inquire(file=sasafile,exist=exists)
      if (exists.EQV..true.) then
        isasa = freeunit()
        open (unit=isasa,file=sasafile,status='old')
        close(unit=isasa,status='delete')
      end if
      if (savreq%instfreq.le.nsim) then
        isasa = freeunit()
        open (unit=isasa,file=sasafile,status='new')
        sasaopen = .true.
      end if
!
! setup filehandle for on-the-fly DSSP output
      inquire(file=dsspfile,exist=exists)
      if (exists.EQV..true.) then
        idssp = freeunit()
        open (unit=idssp,file=dsspfile,status='old')
        close(unit=idssp,status='delete')
      end if
      if ((dsspcalc.le.nsim).AND.(inst_dssp.EQV..true.)) then
        idssp = freeunit()
        open (unit=idssp,file=dsspfile,status='new')
        dsspopen = .true.
      end if
!
! setup filehandle for on-the-fly distance output
      inquire(file=pcfile,exist=exists)
      if (exists.EQV..true.) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='old')
        close(unit=igpc,status='delete')
      end if
      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='new')
        pcopen = .true.
      end if
!
!! setup filehandle for on-the-fly bound partner output FOR THE MOMENT CONTROLLED AS GENERAL_DIST.DAT
!      inquire(file=bprtnrfile,exist=exists)
!      if (exists.EQV..true.) then
!        ibprtnr = freeunit()
!        open (unit=ibprtnr,file=bprtnrfile,status='old')
!        close(unit=ibprtnr,status='delete')
!      end if
!      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
!        ibprtnr = freeunit()
!        open (unit=ibprtnr,file=bprtnrfile,status='new')
!        bpopen = .true.
!      end if
!!
! setup filehandle for on-the-fly RMSD output
      inquire(file=rmsdfile,exist=exists)
      if (exists.EQV..true.) then
        irmsd = freeunit()
        open (unit=irmsd,file=rmsdfile,status='old')
        close(unit=irmsd,status='delete')
      end if
      if ((align%calc.le.nsim).AND.(align%yes.EQV..true.).AND.(align%instrmsd.EQV..true.)) then
        irmsd = freeunit()
        open (unit=irmsd,file=rmsdfile,status='new')
        rmsdopen = .true.
      end if
!
#ifdef ENABLE_THREADS
! setup filehandle for threads reporting
      inquire(file=threadsfile,exist=exists)
      if (exists.EQV..true.) then
        ithread = freeunit()
        open (unit=ithread,file=threadsfile,status='old')
        close(unit=ithread,status='delete')
      end if
      ithread = freeunit()
      open (unit=ithread,file=threadsfile,status='new')
      threadsopen = .true.
#endif
!
! setup filehandle for on-the-fly cross-node energy output
      inquire(file=remcfile,exist=exists)
      if (exists.EQV..true.) then
        iremc = freeunit()
        open (unit=iremc,file=remcfile,status='old')
        close(unit=iremc,status='delete')
      end if
      if ((re_olcalc.le.nsim).AND.(inst_reol.gt.0)) then
        iremc = freeunit()
        open (unit=iremc,file=remcfile,status='new')
        remcopen = .true.
      end if
!
! setup filehandle for on-the-fly RE trace (original structure vs. condition)
      if ((myrank.eq.masterrank).AND.(pdb_analyze.EQV..false.)) then
        inquire(file=retrfile,exist=exists)
        if (exists.EQV..true.) then
          iretr = freeunit()
          open(unit=iretr,file=retrfile,status='old',position='append')
          if (do_restart.EQV..false.) close(unit=iretr,status='delete')
        end if
        inquire(file=retrfile,exist=exists)
        if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.((do_restart.EQV..false.).OR.(exists.EQV..false.))) then
          iretr = freeunit()
          open (unit=iretr,file=retrfile,status='new')
          retropen = .true.
        else if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.(do_restart.EQV..true.).AND.(exists.EQV..true.)) then
          retropen = .true.
          write(ilog,*) 'Warning. Trace file for replica exchange run exists and will be appended (this may be inappropriate).'
          close(unit=iretr)
          open(unit=iretr,file=retrfile,status='old',action='read')
          do while (.true.)
            read(iretr,*,iostat=iomess) mt,mpi_lmap(1:mpi_nodes,1)
            if (iomess.ne.0) exit
            mpi_lmap(1:mpi_nodes,2) = mpi_lmap(1:mpi_nodes,1)
          end do
          close(unit=iretr)
          open(unit=iretr,file=retrfile,status='old',position='append')
        else if ((do_restart.EQV..true.).AND.(exists.EQV..true.)) then
          close(unit=iretr)
        end if
      else if ((myrank.eq.masterrank).AND.(pdb_analyze.EQV..true.).AND.(re_aux(3).eq.1)) then
        call strlims(re_traceinfile,t1,t2)
        inquire(file=re_traceinfile(t1:t2),exist=exists)
        if (exists.EQV..true.) then
          iretr = freeunit()
          open(unit=iretr,file=re_traceinfile(t1:t2),status='old',action='read')
        else
          write(ilog,*) 'Warning. In MPI trajectory analysis mode, requested un"scrambling" of replica exchange trajectories &
 &cannot be performed because trace file cannot be opened (got: ',re_traceinfile(t1:t2),').'
          re_aux(3) = 0
        end if
      else if ((pdb_analyze.EQV..true.).AND.(re_aux(3).eq.1)) then
        call strlims(re_traceinfile,t1,t2)
        inquire(file=re_traceinfile(t1:t2),exist=exists)
        if (exists.EQV..false.) then
          write(ilog,*) 'Warning. In MPI trajectory analysis mode, requested un"scrambling" of replica exchange trajectories &
 &cannot be performed because trace file cannot be opened (got: ',re_traceinfile(t1:t2),').'
          re_aux(3) = 0
        end if
      end if
!
! setup filehandle for temporary TRCV (covariance of torsions) file
      if (covcalc.le.nsim) then
        do mt=1,nangrps
          call int2str(mt,mtstr,lext)
          trcvfile = 'N_'//xpont(1:lext)//'_TRCV_'//mtstr(1:lext)&
 &                                                //'.tmp'
          inquire(file=trcvfile,exist=exists)
          if (exists.EQV..true.) then
            itrcv(mt) = freeunit()
            open (unit=itrcv(mt),file=trcvfile,status='old')
            close(unit=itrcv(mt),status='delete')
          end if
!
          if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
          itrcv(mt) = freeunit()
          open (unit=itrcv(mt),file=trcvfile,status='new')
          trcvopen(mt) = .true.
        end do
      end if
!
! setup filehandle for frame weights from ELS (appended for restarts)
      inquire(file=elsfile,exist=exists)
      if (exists.EQV..true.) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        if ((do_restart.EQV..false.).AND.(hmjam%prtfrmwts.le.nsim)) then
          close(unit=hmjam%iwtsfile,status='delete')
        else
          close(unit=hmjam%iwtsfile)
        end if
      end if
      inquire(file=elsfile,exist=exists)
      if ((hmjam%prtfrmwts.le.nsim).AND.(exists.EQV..false.)) then
        hmjam%iwtsfile = freeunit()
        open (unit=hmjam%iwtsfile,file=elsfile,status='new')
        elsopen = .true.
      else if (hmjam%prtfrmwts.le.nsim) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        elsopen = .true.
      end if
!
!   close files and free up filehandles
    else
      if (fycopen) close(unit=idihed)
      if (polymopen) close(unit=ipolmr)
      if (enopen) close(unit=iene)
      if (accopen) close(unit=iacc)
      if (ensopen) close(unit=iens)       
      if (persopen) close(unit=ipers)
      if (sasaopen) close(unit=isasa)
      if (dsspopen) close(unit=idssp)
      if (phtopen) close(unit=ipht)
      if (remcopen) close(unit=iremc)
      if (retropen) close(unit=iretr)
      if (covcalc.le.nsim) then
        do mt=1,nangrps
          if (trcvopen(mt)) close(unit=itrcv(mt))
        end do
      end if
      if (elsopen) close(unit=hmjam%iwtsfile)
      if (pcopen) close(unit=igpc)
      if (rmsdopen) close(unit=irmsd)
!      if (bpopen) close(unit=ibprtnr)
#ifdef ENABLE_THREADS
      if (threadsopen) close(unit=ithread)
#endif
    end if
!
  else if (use_MPIAVG.EQV..true.) then
!
    if (mode .eq. 1) then
      lext = 3
      call int2str(myrank,xpont,lext)
      dihedfile  = 'FYC.dat'
      polmrfile  = 'POLYMER.dat'
      enfile  = 'ENERGY.dat'
      if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        accfile  = 'ACCEPTANCE.dat'
      end if
      if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
        ensfile = 'ENSEMBLE.dat'
      end if
      persfile = 'PERSISTENCE.dat'
      sasafile = 'SAV.dat'
      phtfile = 'PHTIT.dat'
      retrfile = 'N_'//xpont(1:lext)//'_PIGSTRACE.dat'
      elsfile = 'N_'//xpont(1:lext)//'_ELS_WFRAMES.dat'
      pcfile = 'N_'//xpont(1:lext)//'_GENERAL_DIS.dat'
!      bprtnrfile = 'N_'//xpont(1:lext)//'_BOUND_PARTNERS.dat'
#ifdef ENABLE_THREADS
      threadsfile = 'N_'//xpont(1:lext)//'_THREADS.log'
#endif
!
      if (myrank.eq.masterrank) then
!   setup filehandle for dihedral angle output
        inquire(file=dihedfile,exist=exists)
        if (exists.EQV..true.) then
          idihed = freeunit()
          open (unit=idihed,file=dihedfile,status='old')
          close(unit=idihed,status='delete')
        end if
        if (torout .le. nsim) then
          idihed = freeunit()
          open (unit=idihed,file=dihedfile,status='new')
          fycopen = .true.
        end if
!
!   setup filehandle for polymeric outputs
        inquire(file=polmrfile,exist=exists)
        if (exists.EQV..true.) then
          ipolmr = freeunit()
          open (unit=ipolmr,file=polmrfile,status='old')
          close(unit=ipolmr,status='delete')
        end if
        if (polout .le. nsim) then
          ipolmr = freeunit()
          open (unit=ipolmr,file=polmrfile,status='new')
          polymopen = .true.
        end if
!
!   setup filehandle for energy outputs
        inquire(file=enfile,exist=exists)
        if (exists.EQV..true.) then
          iene = freeunit()
          open (unit=iene,file=enfile,status='old')
          close(unit=iene,status='delete')
        end if
        if (enout.le.nsim) then
          iene = freeunit()
          open (unit=iene,file=enfile,status='new')
          enopen = .true.
        end if
!
!   setup filehandle for acceptance outputs
        if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
          inquire(file=accfile,exist=exists)
          if (exists.EQV..true.) then
            iacc = freeunit()
            open (unit=iacc,file=accfile,status='old')
            close(unit=iacc,status='delete')
          end if
          if (accout.le.nsim) then
            iacc = freeunit()
            open (unit=iacc,file=accfile,status='new')
            accopen = .true.
          end if
        end if
!
!   setup filehandle for ensemble outputs
        if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
          inquire(file=ensfile,exist=exists)
          if (exists.EQV..true.) then
            iens = freeunit()
            open (unit=iens,file=ensfile,status='old')
            close(unit=iens,status='delete')
          end if
          if (ensout.le.nsim) then
            iens = freeunit()
            open (unit=iens,file=ensfile,status='new')
            ensopen = .true.
          end if
        end if
!
!   setup filehandle for persistence length calculation output file
        inquire(file=persfile,exist=exists)
        if (exists.EQV..true.) then
          ipers = freeunit()
          open (unit=ipers,file=persfile,status='old')
          close(unit=ipers,status='delete')
        end if
        if ((polcalc.le.nsim).AND.(doperser.EQV..true.)) then
         ipers = freeunit()
         open (unit=ipers,file=persfile,status='new')
         persopen = .true.
        end if
!
!   setup filehandle for pH titration outputs
        inquire(file=phtfile,exist=exists)
        if (exists.EQV..true.) then
           ipht = freeunit()
           open (unit=ipht,file=phtfile,status='old')
           close(unit=ipht,status='delete')
        end if
        if ((phfreq.gt.0.0d0).AND.(phout.le.nsim)) then
           ipht = freeunit()
           open (unit=ipht,file=phtfile,status='new') 
           phtopen = .true.
        end if 
!
!   setup filehandle for on-the-fly SAV output
        inquire(file=sasafile,exist=exists)
        if (exists.EQV..true.) then
          isasa = freeunit()
          open (unit=isasa,file=sasafile,status='old')
          close(unit=isasa,status='delete')
        end if
        if (savreq%instfreq.le.nsim) then
          isasa = freeunit()
          open (unit=isasa,file=sasafile,status='new')
          sasaopen = .true.
        end if
!
      end if ! if master
!
! setup filehandle for on-the-fly distance output
      inquire(file=pcfile,exist=exists)
      if (exists.EQV..true.) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='old')
        close(unit=igpc,status='delete')
      end if
      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
        igpc = freeunit()
        open (unit=igpc,file=pcfile,status='new')
        pcopen = .true.
      end if
!
!! setup filehandle for on-the-fly bound partner output FOR THE MOMENT CONTROLLED AS GENERAL_DIST.DAT
!      inquire(file=bprtnrfile,exist=exists)
!      if (exists.EQV..true.) then
!        ibprtnr = freeunit()
!        open (unit=ibprtnr,file=bprtnrfile,status='old')
!        close(unit=ibprtnr,status='delete')
!      end if
!      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
!        ibprtnr = freeunit()
!        open (unit=ibprtnr,file=bprtnrfile,status='new')
!        bpopen = .true.
!      end if
!!
#ifdef ENABLE_THREADS
! setup filehandle for threads reporting
      inquire(file=threadsfile,exist=exists)
      if (exists.EQV..true.) then
        ithread = freeunit()
        open (unit=ithread,file=threadsfile,status='old')
        close(unit=ithread,status='delete')
      end if
      ithread = freeunit()
      open (unit=ithread,file=threadsfile,status='new')
      threadsopen = .true.
#endif
!
! setup filehandle for temporary TRCV (covariance of torsions) file
      lext = 3
      if ((covcalc.le.nsim).AND.(myrank.eq.masterrank)) then
        do mt=1,nangrps
          call int2str(mt,mtstr,lext)
          trcvfile = 'TRCV_'//mtstr(1:lext)//'.tmp'
          inquire(file=trcvfile,exist=exists)
          if (exists.EQV..true.) then
            itrcv(mt) = freeunit()
            open (unit=itrcv(mt),file=trcvfile,status='old')
            close(unit=itrcv(mt),status='delete')
          end if
!
          if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
          itrcv(mt) = freeunit()
          open (unit=itrcv(mt),file=trcvfile,status='new')
          trcvopen(mt) = .true.
        end do
      end if
!
! setup filehandle for on-the-fly RE trace (original structure vs. condition)
      if ((myrank.eq.masterrank).AND.(use_MPIMultiSeed.EQV..true.)) then
        inquire(file=retrfile,exist=exists)
        if (exists.EQV..true.) then
          iretr = freeunit()
          open(unit=iretr,file=retrfile,status='old',position='append')
          if ((pdb_analyze.EQV..false.).AND.(do_restart.EQV..false.)) close(unit=iretr,status='delete')
        end if
        inquire(file=retrfile,exist=exists)
        if (pdb_analyze.EQV..false.) then
          if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.((do_restart.EQV..false.).OR.(exists.EQV..false.))) then
            iretr = freeunit()
            open (unit=iretr,file=retrfile,status='new')
            retropen = .true.
          else if ((re_freq.le.nsim).AND.(inst_retr.EQV..true.).AND.(do_restart.EQV..true.).AND.(exists.EQV..true.)) then
            retropen = .true.
            write(ilog,*) 'Warning. Trace file for PIGS run exists and will be appended (this may be inappropriate).'
            close(unit=iretr)
            open(unit=iretr,file=retrfile,status='old',position='append')
          else if ((do_restart.EQV..true.).AND.(exists.EQV..true.)) then
            close(unit=iretr)
          end if
        else
          if ((exists.EQV..true.).AND.(inst_retr.EQV..true.)) then
            retropen = .true.
            write(ilog,*) 'Warning. Existing trace file will be appended by PIGS analysis run (this may be inappropriate).'
            close(unit=iretr)
            open(unit=iretr,file=retrfile,status='old',position='append')
          else if (inst_retr.EQV..false.) then
            close(unit=iretr)
          else
            iretr = freeunit()
            open (unit=iretr,file=retrfile,status='new')
            retropen = .true.
          end if
        end if
      end if
!
! setup filehandle for frame weights from ELS (appended)
      inquire(file=elsfile,exist=exists)
      if (exists.EQV..true.) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        if ((do_restart.EQV..false.).AND.(hmjam%prtfrmwts.le.nsim)) then
          close(unit=hmjam%iwtsfile,status='delete')
        else
          close(unit=hmjam%iwtsfile)
        end if
      end if
      inquire(file=elsfile,exist=exists)
      if ((hmjam%prtfrmwts.le.nsim).AND.(exists.EQV..false.)) then
        hmjam%iwtsfile = freeunit()
        open (unit=hmjam%iwtsfile,file=elsfile,status='new')
        elsopen = .true.
      else if (hmjam%prtfrmwts.le.nsim) then
        hmjam%iwtsfile = freeunit()
        open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
        elsopen = .true.
      end if
!
    else
      if (fycopen) close(unit=idihed)
      if (polymopen) close(unit=ipolmr)
      if (enopen) close(unit=iene)
      if (accopen) close(unit=iacc)
      if (ensopen) close(unit=iens)
      if (persopen) close(unit=ipers)
      if (sasaopen) close(unit=isasa)
      if (phtopen) close(unit=ipht)
      if (covcalc.le.nsim) then
        do mt=1,nangrps
          if (trcvopen(mt)) close(unit=itrcv(mt))
        end do
      end if
      if (elsopen) close(unit=hmjam%iwtsfile)
      if (pcopen) close(unit=igpc)
!      if (bpopen) close(unit=ibprtnr)
#ifdef ENABLE_THREADS
      if (threadsopen) close(unit=ithread)
#endif
    end if
!
  end if
#else
  if(mode .eq. 1) then
    dihedfile  = 'FYC.dat'
    polmrfile  = 'POLYMER.dat'
    holesfile = 'HOLES.dat'
    enfile  = 'ENERGY.dat'
    if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      accfile  = 'ACCEPTANCE.dat'
    end if
    if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      ensfile = 'ENSEMBLE.dat'
    end if
    persfile = 'PERSISTENCE.dat'
    sasafile = 'SAV.dat'
    dsspfile = 'DSSP_RUNNING.dat'
    pcfile = 'GENERAL_DIS.dat'
    trcvfile = 'TRCV.tmp'
    phtfile = 'PHTIT.dat'
    elsfile = 'ELS_WFRAMES.dat'
    rmsdfile = 'RMSD.dat'
!    bprtnrfile = 'BOUND_PARTNERS.dat'
#ifdef ENABLE_THREADS
    threadsfile = 'THREADS.log'
#endif
!   
! setup filehandle for dihedral angle output
    inquire(file=dihedfile,exist=exists)
    if (exists.EQV..true.) then
       idihed = freeunit()
       open (unit=idihed,file=dihedfile,status='old')
       close(unit=idihed,status='delete')
    end if
    if (torout .le. nsim) then
       idihed = freeunit()
       open (unit=idihed,file=dihedfile,status='new')
       fycopen = .true.
    end if
!
! setup filehandle for polymeric outputs
    inquire(file=polmrfile,exist=exists)
    if (exists.EQV..true.) then
       ipolmr = freeunit()
       open (unit=ipolmr,file=polmrfile,status='old')
       close(unit=ipolmr,status='delete')
    end if
    if (polout .le. nsim) then
       ipolmr = freeunit()
       open (unit=ipolmr,file=polmrfile,status='new') 
       polymopen = .true.
    end if  
!   
! setup filehandle for holes outputs
    inquire(file=holesfile,exist=exists)
    if (exists.EQV..true.) then
      iholes = freeunit()
      open (unit=iholes,file=holesfile,status='old')
      close(unit=iholes,status='delete')
    end if
    if (holescalc.le.nsim) then
      iholes = freeunit()
      open (unit=iholes,file=holesfile,status='new') 
      holesopen = .true.
    end if 
!   
! setup filehandle for energy outputs
    inquire(file=enfile,exist=exists)
    if (exists.EQV..true.) then
       iene = freeunit()
       open (unit=iene,file=enfile,status='old')
       close(unit=iene,status='delete')
    end if
    if (enout.le.nsim) then
      iene = freeunit()
      open (unit=iene,file=enfile,status='new')     
      enopen = .true.
    end if
!
! setup filehandle for acceptance outputs
    if ((pdb_analyze.EQV..false.).AND.((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
      inquire(file=accfile,exist=exists)
      if (exists.EQV..true.) then
        iacc = freeunit()
        open (unit=iacc,file=accfile,status='old')
        close(unit=iacc,status='delete')
      end if
      if (accout.le.nsim) then
        iacc = freeunit()
        open (unit=iacc,file=accfile,status='new')
        accopen = .true.
      end if
    end if
!
! setup filehandle for ensemble outputs
    if ((pdb_analyze.EQV..false.).AND.(((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))) then
      inquire(file=ensfile,exist=exists)
      if (exists.EQV..true.) then
        iens = freeunit()
        open (unit=iens,file=ensfile,status='old')
        close(unit=iens,status='delete')
      end if
      if (ensout.le.nsim) then
        iens = freeunit()
        open (unit=iens,file=ensfile,status='new')
        ensopen = .true.
      end if
    end if
!
! setup filehandle for persistence length calculation output file
    inquire(file=persfile,exist=exists)
    if (exists.EQV..true.) then
       ipers = freeunit()
       open (unit=ipers,file=persfile,status='old')
       close(unit=ipers,status='delete')
    end if
    if ((polcalc.le.nsim).AND.(doperser.EQV..true.)) then
      ipers = freeunit()
      open (unit=ipers,file=persfile,status='new')
      persopen = .true.
    end if
!
! setup filehandle for pH titration outputs
    
    if ((phfreq.gt.0.0).AND.(phout.le.nsim)) then
       inquire(file=phtfile,exist=exists)
       if (exists.EQV..true.) then
          ipht = freeunit()
          open (unit=ipht,file=phtfile,status='old')
          close(unit=ipht,status='delete')
       end if
       ipht = freeunit()
       open (unit=ipht,file=phtfile,status='new') 
       phtopen = .true.
    end if
!
! setup filehandle for on-the-fly SAV output
    inquire(file=sasafile,exist=exists)
    if (exists.EQV..true.) then
       isasa = freeunit()
       open (unit=isasa,file=sasafile,status='old')
       close(unit=isasa,status='delete')
    end if
    if (savreq%instfreq.le.nsim) then
      isasa = freeunit()
      open (unit=isasa,file=sasafile,status='new')
      sasaopen = .true.
    end if
!
! setup filehandle for on-the-fly DSSP output
    inquire(file=dsspfile,exist=exists)
    if (exists.EQV..true.) then
      idssp = freeunit()
      open (unit=idssp,file=dsspfile,status='old')
      close(unit=idssp,status='delete')
    end if
    if ((dsspcalc.le.nsim).AND.(inst_dssp.EQV..true.)) then
      idssp = freeunit()
      open (unit=idssp,file=dsspfile,status='new')
      dsspopen = .true.
    end if
!
! setup filehandle for on-the-fly distance output
    inquire(file=pcfile,exist=exists)
    if (exists.EQV..true.) then
      igpc = freeunit()
      open (unit=igpc,file=pcfile,status='old')
      close(unit=igpc,status='delete')
    end if
    if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
      igpc = freeunit()
      open (unit=igpc,file=pcfile,status='new')
      pcopen = .true.
    end if
!
! setup filehandle for on-the-fly RMSD output
    inquire(file=rmsdfile,exist=exists)
    if (exists.EQV..true.) then
      irmsd = freeunit()
      open (unit=irmsd,file=rmsdfile,status='old')
      close(unit=irmsd,status='delete')
    end if
    if ((align%calc.le.nsim).AND.(align%yes.EQV..true.).AND.(align%instrmsd.EQV..true.)) then
      irmsd = freeunit()
      open (unit=irmsd,file=rmsdfile,status='new')
      rmsdopen = .true.
    end if
!
!! setup filehandle for on-the-fly bound partner output FOR THE MOMENT CONTROLLED AS GENERAL_DIST.DAT
!      inquire(file=bprtnrfile,exist=exists)
!      if (exists.EQV..true.) then
!        ibprtnr = freeunit()
!        open (unit=ibprtnr,file=bprtnrfile,status='old')
!        close(unit=ibprtnr,status='delete')
!      end if
!      if ((pccalc.le.nsim).AND.(inst_gpc.gt.0).AND.(inst_gpc.le.(nsim/pccalc)).AND.(gpc%nos.gt.0)) then
!        ibprtnr = freeunit()
!        open (unit=ibprtnr,file=bprtnrfile,status='new')
!        bpopen = .true.
!      end if
!!
#ifdef ENABLE_THREADS
! setup filehandle for threads reporting
    inquire(file=threadsfile,exist=exists)
    if (exists.EQV..true.) then
      ithread = freeunit()
      open (unit=ithread,file=threadsfile,status='old')
      close(unit=ithread,status='delete')
    end if
    ithread = freeunit()
    open (unit=ithread,file=threadsfile,status='new')
    threadsopen = .true.
#endif
!
! setup filehandle for temporary TRCV (covariance of torsions) file
    lext = 3
    if (covcalc.le.nsim) then 
      do mt=1,nangrps
        call int2str(mt,mtstr,lext)
        trcvfile = 'TRCV_'//mtstr(1:lext)//'.tmp'
        inquire(file=trcvfile,exist=exists)
        if (exists.EQV..true.) then
          itrcv(mt) = freeunit()
          open (unit=itrcv(mt),file=trcvfile,status='old')
          close(unit=itrcv(mt),status='delete')
        end if
!
        if (do_pol(moltypid(molangr(mt,1))).EQV..false.) cycle
        itrcv(mt) = freeunit()
        open (unit=itrcv(mt),file=trcvfile,status='new')
        trcvopen(mt) = .true.
      end do
    end if
!
! setup filehandle for frame weights from ELS (appended)
    inquire(file=elsfile,exist=exists)
    if (exists.EQV..true.) then
      hmjam%iwtsfile = freeunit()
      open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
      if ((do_restart.EQV..false.).AND.(hmjam%prtfrmwts.le.nsim)) then
        close(unit=hmjam%iwtsfile,status='delete')
      else
        close(unit=hmjam%iwtsfile)
      end if
    end if
    inquire(file=elsfile,exist=exists)
    if ((hmjam%prtfrmwts.le.nsim).AND.(exists.EQV..false.)) then
      hmjam%iwtsfile = freeunit()
      open (unit=hmjam%iwtsfile,file=elsfile,status='new')
      elsopen = .true.
    else if (hmjam%prtfrmwts.le.nsim) then 
      hmjam%iwtsfile = freeunit()
      open(unit=hmjam%iwtsfile,file=elsfile,status='old',position='append')
      elsopen = .true.
    end if
!
! close files and free up filehandles
  else
    if (fycopen) close(unit=idihed)
    if (polymopen) close(unit=ipolmr)
    if (enopen) close(unit=iene)
    if (accopen) close(unit=iacc)
    if (persopen) close(unit=ipers)
    if (ensopen) close(unit=iens)
    if (sasaopen) close(unit=isasa)
    if (dsspopen) close(unit=idssp)
    if (phtopen) close (unit=ipht)
    if (covcalc.le.nsim) then
      do mt=1,nangrps
        if (trcvopen(mt).EQV..true.) close(unit=itrcv(mt))
      end do
    end if
    if (elsopen) close(unit=hmjam%iwtsfile)
    if (pcopen) close(unit=igpc)
    if (rmsdopen) close(unit=irmsd)
!    if (bpopen) close(unit=ibprtnr)
#ifdef ENABLE_THREADS
    if (threadsopen) close(unit=ithread)
#endif
  end if
!
#endif
!
end
!
!------------------------------------------------------------------
!
#ifdef ENABLE_MPI
!
subroutine makelogio(mode)
!
  use iounit
  use mcsums
  use mpistuff
!
  implicit none
!
  logical exists
  integer freeunit,lext,mode
  character(60) mclogfile
  character(3) xpont
!
  if (mode.eq.1) then
!
    lext = 3
    call int2str(myrank,xpont,lext) 
    mclogfile ='N_'//xpont(1:lext)//'.log' 
!
!   setup filehandle for log-file
    inquire(file=mclogfile,exist=exists)
    if (exists.EQV..true.) then
      ilog = freeunit()
      open (unit=ilog,file=mclogfile,status='old')
      close(unit=ilog,status='delete')
    end if
    ilog = freeunit()
    open (unit=ilog,file=mclogfile,status='new')
!
  else
!
    inquire(unit=ilog,opened=exists)
    if (exists.EQV..true.) close(unit=ilog)
!
  end if
end
!
#endif
!
!--------------------------------------------------------------------
!
subroutine flushopen()
!
  use mcsums
  use iounit
  use threads, ONLY: ithread
  use molecule, ONLY: nangrps
  use system, ONLY: nsim
  use molecule, ONLY: nangrps
  use torsn, ONLY: trcvopen,itrcv,covcalc
  use wl, ONLY: hmjam
!
  implicit none
!
  logical isop
  integer mt
!
  inquire(unit=ilog,opened=isop)
  if (isop.EQV..true.) flush(ilog)
#ifdef ENABLE_MPI
  inquire(unit=iremc,opened=isop)
  if (isop.EQV..true.) flush(iremc)
  inquire(unit=iretr,opened=isop)
  if (isop.EQV..true.) flush(iretr)
#endif
#ifdef ENABLE_THREADS
  inquire(unit=ithread,opened=isop)
  if (isop.EQV..true.) flush(ithread)
#endif
!
  inquire(unit=idihed,opened=isop)
  if (isop.EQV..true.) flush(idihed)
  inquire(unit=idihed_ring,opened=isop)
  if (isop.EQV..true.) flush(idihed_ring)
  inquire(unit=ipolmr,opened=isop)
  if (isop.EQV..true.) flush(ipolmr)
  inquire(unit=iholes,opened=isop)
  if (isop.EQV..true.) flush(iholes)
  inquire(unit=ipht,opened=isop)
  if (isop.EQV..true.) flush(ipht)
  inquire(unit=idssp,opened=isop)
  if (isop.EQV..true.) flush(idssp)
  inquire(unit=iacc,opened=isop)
  if (isop.EQV..true.) flush(iacc)
  inquire(unit=ipers,opened=isop)
  if (isop.EQV..true.) flush(ipers)
  inquire(unit=isasa,opened=isop)
  if (isop.EQV..true.) flush(isasa)
  inquire(unit=iene,opened=isop)
  if (isop.EQV..true.) flush(iene)
!  inquire(unit=ipdbtraj,opened=isop)
!  if (isop.EQV..true.) flush(ipdbtraj)
  inquire(unit=iens,opened=isop)
  if (isop.EQV..true.) flush(iens)
  inquire(unit=igpc,opened=isop)
  if (isop.EQV..true.) flush(igpc)
  inquire(unit=irmsd,opened=isop)
  if (isop.EQV..true.) flush(irmsd)
  if (covcalc.le.nsim) then
    do mt=1,nangrps
      if (trcvopen(mt).EQV..true.) flush(unit=itrcv(mt))
    end do
  end if
  inquire(unit=hmjam%iwtsfile,opened=isop)
  if (isop.EQV..true.) flush(hmjam%iwtsfile)
!
end
!
!--------------------------------------------------------------------
!
!
! ####################################################
! ##                                                ##
! ## subroutine contours -- file I/O for dipeptides ##
! ##                                                ##
! ####################################################
!
! "contours" opens or closes output files generating 
! dipeptide maps
! mode = 1: open
! mode = 2: or anything else close
!
!  subroutine contours(mode)
!  implicit none
!  include 'files.i"
!  integer mode,freeunit
!  integer ien,ignorm,is,iphi,ipsi
!  character*60 enfile,gnormfile
!  character*60 phifile,psifile
!  character*60 basinfile
!  logical exists
!c     
!c     open output files for explicit or averaged energies
!  if(mode .eq. 1) then
!     enfile = filename(1:leng)//'En.dat'
!     inquire(file=enfile,exist=exists)
!     if (exists.EQV..true.) then
!        ien = freeunit()
!        open (unit=ien,file=enfile,status='old')
!        close(unit=ien,status='delete')
!     end if
!     ien = freeunit()
!     open (unit=ien,file=enfile,status='new')
!c     
!c     open output files to store gradient norms
!     gnormfile = filename(1:leng)//'Gnorm.dat'
!     inquire(file=gnormfile,exist=exists)
!     if (exists.EQV..true.) then
!        ignorm = freeunit()
!        open (unit=ignorm,file=gnormfile,status='old')
!        close(unit=ignorm,status='delete')
!     end if
!     ignorm = freeunit()
!     open (unit=ignorm,file=gnormfile,status='new')
!c     
!c     open file to store phi angles
!     phifile = filename(1:leng)//'Phi.dat'
!     inquire(file=phifile,exist=exists)
!     if (exists.EQV..true.) then
!        iphi = freeunit()
!        open (unit=iphi,file=phifile,status='old')
!        close(unit=iphi,status='delete')
!     end if
!     iphi = freeunit()
!     open (unit=iphi,file=phifile,status='new')
!c         
!c     open file to store psi angles
!     psifile = filename(1:leng)//'Psi.dat'
!     inquire(file=psifile,exist=exists)
!     if (exists.EQV..true.) then
!c           ipsi = freeunit()
!        open (unit=ipsi,file=psifile,status='old')
!        close(unit=ipsi,status='delete')
!     end if
!     ipsi = freeunit()
!     open (unit=ipsi,file=psifile,status='new')
!c         
!c     open file to store inherent structures
!     basinfile = filename(1:leng)//'IS.dat'
!     inquire(file=basinfile,exist=exists)
!     if (exists.EQV..true.) then
!        is = freeunit()
!        open (unit=is,file=basinfile,status='old')
!        close(unit=is,status='delete')
!     end if
!     is = freeunit()
!     open (unit=is,file=basinfile,status='new')
!  else
!     close(unit=iphi)
!     close(unit=ipsi)
!     close(unit=ien)
!     close(unit=ignorm)
!     close(unit=is)
!  end if
!  return
!  end
!
