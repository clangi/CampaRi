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
!
#include "macros.i"
!
!
! ############################################
! ##                                        ##
! ## THIS IS THE CAMPARI MASTER SOURCE FILE ##
! ##                                        ##
! ############################################
!
!
program chainsaw
!
  use iounit
  use torsn
  use mcsums
  use energies
  use atoms
  use sequen
  use mcgrid
  use pdb
  use molecule
  use polyavg
  use movesets
  use paircorr
  use mpistuff
  use cutoffs
  use system
  use contacts
  use dipolavg
  use forces
  use diffrac
  use mini
  use dssps
  use clusters
  use shakeetal
  use grandensembles
  use ems
  use ewalds
  use interfaces
  use, INTRINSIC:: ISO_C_BINDING
  use threads
  use zmatrix
  use ncdm
  use fos, ONLY: savreq
  use, INTRINSIC:: ISO_FORTRAN_ENV
!  use math
!
  implicit none
!
  integer dt(8),modstep
  integer(KIND=8) ee2,ee1,t1,t2
  character (len=12) rc(3)
  integer rs,azero,aone,imol,tpi
  integer istep,lstep,ndump,m1
  integer firststep
  RTYPE ee,random,eit
  RTYPE edum(MAXENERGYTERMS),force1,force3
  RTYPE, ALLOCATABLE:: blas(:)
  RTYPE ed1(MAXENERGYTERMS),ed2(MAXENERGYTERMS)
  RTYPE tt!,fs_ipol,sq_ipol,sqdr_ipol,fsdr_ipol
#ifdef ENABLE_THREADS
  integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,tpn
#endif
  character(MAXSTRLEN) intfile
  logical logdummy,atrue,afalse,logd2,logbuff
#ifdef ENABLE_MPI
  integer masterrank,ierr
#else
  integer st1,st2
  character(len=1,kind=c_char), dimension(MAXSTRLEN+1):: honam(0:MAXSTRLEN)
  character(MAXSTRLEN) hnam
  integer(c_int) hnerr
  integer(c_size_t) hnml
#endif
!
#ifdef ENABLE_MPI
  masterrank = 0
#endif
!
  allocate(blas(3))
  atrue = .true.
  afalse = .false.
  azero = 0
  aone = 1
  m1 = -1
#ifdef ENABLE_MPI
!
! MPI: obviously some things need to be different ...
! this routine will have the master distribute the simulation information
! and run the setup procedures for all CPUs
  time_comm = 0.0
  call System_Clock(ee1)
  call MPI_STARTMC()
  call System_Clock(t2)
  time_comm = time_comm + 1.0*(t2 - ee1)
  call System_Clock(ee2)
#else
!
! when we start
  hnml = MAXSTRLEN
  honam(:) = ' '
  hnerr = handwrapped_gethostname(honam(0:(MAXSTRLEN-1)),hnml)
  do rs=1,MAXSTRLEN
    hnam(rs:rs) = honam(rs-1)
  end do
  call strlims(hnam,st1,st2)
  if (st2.gt.st1) then
    if (hnam(st2:st2).eq.C_NULL_CHAR) st2 = st2 - 1
  end if
  call Date_and_Time(rc(1),rc(2),rc(3),dt)
  write(*,*)
  call write_license()
  write(*,*) 'Execution started on host ',hnam(st1:st2),' at ',rc(2)(1:2),':',rc(2)(3:4),&
 &' on ',rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
  call System_Clock(ee1)
  ee2 = ee1
  write(*,*)
!     
! all initializations for start up
  call initial()
  logdummy = .false.
!
! we need to read (not yet parse) the key (input) file
  call getkey()
  call init_campprng(aone)
!
! level 1 input reader
  write(ilog,*)
  write(ilog,*) '---   Now parsing keywords ...            ---'
  write(ilog,*)
  call parsekey(1)
! level 2 input reader
  call parsekey(2)
! level 3 input reader
  call parsekey(3)
  write(ilog,*)
  write(ilog,*) '---   ... finished parsing keywords.      ---'
  write(ilog,*)
#endif
!
  logbuff = use_trajidx
  use_trajidx = .false.
!
! set the derived box variables
  call update_bound(atrue)
!
! build the system 
  call makepept()
!
! resolve problems with requests for structural input
  call strucinp_sanitychecks()
!
  if (use_cutoffs.EQV..false.) then
!   overwrite checkfreq
    if (dyn_mode.eq.1) then
!      nsancheck = nsim + 1
    end if
  end if 
!
  if (use_POLAR.EQV..true.) then
    call polar_groups()
  end if
!
  if (use_IMPSOLV.EQV..true.) then
    call setup_freesolv()
    call solvation_groups()
!    do rs=1,1000
!      atsav(1) = (rs/1000.)*atbvol(1)
!      write(0,*) rs/1000.,fs_ipol(1),sq_ipol(1)
!    end do
  end if
!
  if (use_TABUL.EQV..true.) then
    call read_tabfiles()
  end if
!
!  if (use_POLAR.EQV..true.) then
!    call setup_ionloops() ! requires TABUL setup to be complete
!  end if
!
  if (use_OSMO.EQV..true.) then
    call read_osmofile()
  end if
!
  if (use_POLY.EQV..true.) then
    call read_polfile()
  end if
!
  if (use_DREST.EQV..true.) then
    call read_drestfile()
  end if
!
  call assign_bndtprms()
!
  if (phfreq.gt.0.0) then
     call ionize_setup()
  end if
!
! Particle fluctuation setup also requires setup to be complete
  if ((ens%flag.ge.5).AND.(ens%flag.le.6)) then
    call read_particleflucfile()
  end if
!
  if (do_restart.EQV..true.) then
!   do nothing (see below)
  else if (pdbinput.EQV..true.) then
    if (pdb_readmode.eq.1) then
      call FMSMC_readpdb()
    else
      call FMSMC_readpdb3()
    end if
  end if

  call randomize_bb() ! this fxn may not change anything
  if (n_crosslinks.gt.0) then
    call correct_crosslinks()
  end if
!
! torsional setup requires all structure manipulation to be complete
  if (use_TOR.EQV..true.) then
    call read_torfile()
  end if
!
! FEG requires all other energy setup to be complete
  if (use_FEG.EQV..true.) then
    call read_fegfile()
  end if
!
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call remc_sanitychecks()
  end if
#endif
!
! if equilibration run manually turn off analysis, so we don't do
! unnecessary and potentially ill-defined setup work
  if (nequil.ge.nsim) then
    pccalc = nsim+1
    covcalc = nsim+1
    segcalc = nsim+1
    torlccalc = nsim+1
    polcalc = nsim+1
    rhcalc = nsim+1
    sctcalc = nsim+1
    holescalc = nsim+1
    savcalc = nsim+1
    contactcalc = nsim+1
    particlenumcalc = nsim+1
    clucalc = nsim+1
    angcalc = nsim+1
    dipcalc = nsim+1
    intcalc = nsim+1
    dsspcalc = nsim+1
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
!     do nothing - needed for PIGS
    else
      cstorecalc = nsim+1
    end if
#else
    cstorecalc = nsim+1
#endif
    xyzout = nsim+1
    phout = nsim+1
    diffrcalc = nsim+1
    emcalc = nsim+1
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      re_olcalc = nsim+1
    end if
#endif
    if (pdb_analyze.EQV..true.) then
      write(ilog,*) 'Warning. No analysis performed in PDB analysis &
 &mode due to equilibration period exceeding snapshot number.'
      write(ilog,*)
    end if
  end if
!
  if (pdb_analyze.EQV..true.) then
    if (use_dyn.EQV..true.) then
      write(ilog,*) 'Warning. Reconstruction of dynamical variables in PDB an&
 &alysis mode not yet supported. Nonetheless using force routines to compute energies (if &
 &applicable).'
      write(ilog,*)
    end if
    if (ens%flag.ne.1) then
      write(ilog,*) 'Fatal. Ensembles other than NVT are currently n&
 &ot supported by PDB analysis mode.'
      write(ilog,*)
      call fexit()
    end if
  end if
!
  if (just_solutes.EQV..true.) then
    if (nsolutes.le.0) then
      write(ilog,*) 'Warning. Suppressing structural output for a system &
 &composed exclusively of solvent is not supported. Turning off suppression&
 & flag.'
      just_solutes = .false.
      pdbeffn = n
    end if
  end if
!
! temporary shutdown of LCT stuff
  if (use_LCTOR.EQV..true.) then
    write(ilog,*) 'Warning. LCT potential term is currently not supp&
 &orted. Turning off.'
    use_LCTOR = .false.
    scale_LCTOR = 0.0
  end if
  if (use_lctmoves.EQV..true.) then
    write(ilog,*) 'Warning. LCT moves are currently not supported. T&
 &urning off.'
    use_lctmoves = .false.
  end if
  if (torlccalc.le.nsim) then
    write(ilog,*) 'Warning. LCT analysis is currently not supported.&
 & Turning off.'
    torlccalc = nsim + 1
  end if
!
! pdb analysis mode requires more setup work
  if (pdb_analyze.EQV..true.) then
!    call randomize_bb()
    if (align%yes.EQV..true.) then
      call read_alignfile(align)
      if (pdb_fileformat.le.2) then
        if (use_pdb_template.EQV..true.) then
          call read_pdb_template()
        end if
      end if
    end if
    if ((cmode.eq.7).AND.(cstorecalc.le.nsim)) then
!     do nothing
    else
      if (pdb_fileformat.eq.1) call setup_pdbtraj()
      if (pdb_fileformat.eq.2) call setup_pdbsnaps()
#ifdef LINK_XDR
      if (pdb_fileformat.eq.3) call setup_xtctraj()
#endif
      if (pdb_fileformat.eq.4) call setup_dcdtraj()
#ifdef LINK_NETCDF
      if (pdb_fileformat.eq.5) call setup_netcdftraj()
#endif
#ifdef ENABLE_MPI
      if ((select_frames.EQV..true.).AND.(use_REMC.EQV..false.)) then
        write(ilog,*) 'Warning. Trajectory subset analysis of user-selected framesfile is currently not &
 &supported in MPI averaging-type trajectory analysis mode. Disabled.'
        select_frames = .false.
      end if
#else
      if (select_frames.EQV..true.) then
        curframe = 0
        call read_framesfile()
      end if
#endif
    end if
  end if
!
  if (do_restart.EQV..false.) then
!   dump out a pdb- and Z matrix-file for sanity check
    if ((pdb_analyze.EQV..false.).OR.(n_pdbunk.gt.0)) then
      intfile = 'START '
      call FMCSC_dump(intfile)
    end if
  end if
!
! read in file containing information about regions in
! phi/psi space to sample (center of boxes, boxsizes)
  if (use_stericgrids.EQV..true.) then
    call allocate_stericgrid(aone)
    call readgrid()
  end if
!
! allocate torsional analysis arrays (up here because of ZSEC-
! coupling)
  call allocate_torsn(aone)
!
! setup for LC torsional analysis
  if ((covcalc.le.nsim).AND.(covmode.eq.3)) then
    call setup_lc_tor()
  end if
!
! setup work for DSSP analysis (here because of DSSP-restraint coupling)
  if ((dsspcalc.le.nsim).OR.(use_DSSP.EQV..true.)) then
    call setup_dssp()
  end if
!
#ifdef LINK_NETCDF
!
! read in file-based EM restraint map
  if (use_EMICRO.EQV..true.) then
    call read_netcdf3Dmap()
  end if
!
#endif
!
!
! now we proceed through some sanity checks for the completed settings
!
! setup for domain grids
  if ((do_restart.EQV..false.).AND.(use_mcgrid.EQV..true.)) then
    call allocate_mcgrid(aone)
    call setupmcgrid()
    if (grid%report.EQV..true.) then
      call gridreport()
    end if
  end if
!
! setup work for rigid-body moves
  rblst%nr = nmol
  allocate(rblst%idx(rblst%nr))
  do imol=1,nmol
    rblst%idx(imol) = imol
  end do
!
! Hamiltonian sanity checks
  call hamiltonian_sanitychecks()
!
! a sanity check against a meaningless calculation
  if ((fycxyz.eq.1).AND.(ntorpuck.eq.0).AND.(nmol.eq.1)) then
    write(ilog,*) 'Warning. Attempting a calculation with no relevan&
 &t degrees of freedom (1 rigid molecule in a box). Input error?'
  end if
!
! big sanity check routine for move set
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    call moveset_sanitychecks()
  end if
!
! setup work for torsional dynamics (must precede read_frzfile())
  if (((dyn_mode.ne.1).AND.(fycxyz.eq.1)).OR.((cstorecalc.le.nsim).AND.&
 &((cdis_crit.eq.2).OR.(cdis_crit.eq.4)))) then
    call set_IMDalign()
  end if
!
! handling of constraints for internal coordinate spaces and adjustment of MC sampling weights
  if (do_frz.EQV..true.) then
    call read_frzfile()
  end if
  call preferential_sampling_setup()
!
! setup work for system-level constraints
  if (dyn_mode.ne.1) then
    call set_movingmass()
    call dimensionality_check()
  end if
!
  call dynamics_sanitychecks()
  call loop_checks()
!
! setup work for constraints in Cartesian dynamics
  if ((dyn_mode.ne.1).AND.(fycxyz.eq.2)) then
    call shake_setup()
  end if
!
! more IMD setup/reporting work
  if ((dyn_mode.ne.1).AND.(fycxyz.eq.1)) then
    call set_mass_fudges()
  end if
  if (izrot_report.EQV..true.) then
    call rotlst_report()
  end if
!
! T-coupling groups for thermostatting (if nothing else, at least for analysis)
  if (dyn_mode.ne.1) then
    call read_tgrpfile()
    call init_thermostat()
  end if
!
! allocation for movesets
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    call allocate_moveset(aone)
  end if
!
! setup work for Ewald sums (might need loop_checks to check on Hamiltonian)
! and (G)RF electrostatics
  if (dyn_mode.ne.1) then
    if (lrel_md.eq.2) then
      call setup_ewald()
    else if (lrel_md.eq.3) then
      call grf_setup()
      if (use_FEG.EQV..true.) then
        call setup_rfcnst()
      end if
    end if
  end if
  call gentables_erfcoverr(ewpm,ewerfctol)
!
  if (dyn_mode.ne.1) then 
    if (use_TABUL.EQV..true.) then
      do rs=1,nseq
        call tab_respairs_nbl_pre(rs)
      end do
    end if
  end if
!
! setup work for WL sampling
  if (mc_acc_crit.eq.3) call wl_init()
!
! setup work for energy landscape sculpting
  call els_init()
!
  call allocate_threads(aone)
! 
! monitor initial setup time
  call System_Clock(t2)
  eit = 1.0*(t2 - ee2)
!
! monitor time for initial energy/force calculation
  if ((do_restart.EQV..false.).AND.(pdb_analyze.EQV..false.)) then
!
    call System_Clock(t1)
!
!    evtl.y perform a gradient check
    if ((use_dyn.EQV..true.).AND.(grad_check.EQV..true.)) then
!     for now, just use these hard-coded parameters
      blas(1) = 0.0000001
      blas(2) = 0.0001
      blas(3) = 0.0001
      call gradient_test(blas(1),blas(2),blas(3))
    end if
!
!   use force routines unless in straight Monte Carlo
    if (dyn_mode.eq.1) then
      if (do_n2loop.EQV..true.) then
        call energy(esterms,esave,azero)
      end if
      if (use_cutoffs.EQV..true.) then
        call energy3(edum,atrue,esavec,azero)
        if (do_n2loop.EQV..false.) then
          write(ilog,*) 'Warning. Using cutoff-based initial energy also as "reference" energy (due to FMCSC_N2LOOP).'
          esave = esavec
          esterms(:) = edum(:)
        end if
      end if
    else 
      if (do_n2loop.EQV..true.) then
        esave = force1(esterms)
      end if
      if (use_cutoffs.EQV..true.) then
        esavec = force3(edum,ed1,ed2,atrue)
        if (do_n2loop.EQV..false.) then
          write(ilog,*) 'Warning. Using cutoff-based initial energy also as "reference" energy (due to FMCSC_N2LOOP).'
          esave = esavec
          esterms(:) = edum(:)
        end if
      end if
    end if
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
    call System_Clock(t1) 
!
!   setup work for rigid body moves (not sure whether or why this would be needed)
    if (rigidfreq.gt.0.0) then
      do imol=1,nmol
        call update_rigid(imol)
      end do
    end if
!
!   cycle initialization for hybrid methods (this is overwritten in restart cases)
    if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      in_dyncyc = .false.
      curcyc_end = first_mccyclen
      curcyc_start = 1
      curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_MPI
      if (use_REMC.EQV..true.) then
        if (mod(curcyc_end,re_freq).eq.0) then
          curcyc_end = curcyc_end + 1
        end if
      end if
#endif
    end if
  else
    call System_Clock(t1)
  end if
!
! setup for pc analysis
  if (pccalc.le.nsim) then
    call read_gpcfile()
  end if
!
! sanity checks and setup for structural clustering
#ifdef ENABLE_MPI
  if ((cstorecalc.le.nsim).OR.((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.))) then
#else
  if (cstorecalc.le.nsim) then
#endif
    call read_clusteringfile()
  end if
!
! setup work for internal coordinate histograms
  if (intcalc.le.nsim) then
    call setup_inthists()
  end if
!
! additional SAV requests
  if (savcalc.le.nsim) then
    call read_savreqfile()
  end if
!
! various sanity checks for requested analyses
  call analysis_checks()
!
! set up file I/O (get handles/open files)
  call makeio(aone)
!
! setup for indexed trajectory output
  use_trajidx = logbuff
  if (use_trajidx.EQV..true.) call read_trajidxfile()
!
! setup for covariance analysis (needs I/O to be done)
  if (covcalc.le.nsim) then
    call setup_covar_int()
  end if
!
! setup for pH calculations (mostly informational output) (needs I/O to be done)
  if (phfreq.gt.0.0) then
     call ph_setup()
  end if
!
! setup for diffraction calculations
  if (diffrcalc.le.nsim) then
    call bessel_read()
  end if
!
! setup for torsional output (needs I/O to be done)
  if (torout.le.nsim) then
    call torsion_header()
  end if
!
! setup for distance set for alignment when there is no PDB ref 
  if ((align%yes.EQV..true.).AND.(align%instrmsd.EQV..true.).AND.(align%refset.EQV..false.)) then
    call init_aligndisset()
  end if
!
#ifdef ENABLE_MPI
  if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.)) then
    call PIGS_sanitychecks()
  end if
#endif
!
! before we proceed to the actual calculation we take care of
! all memory allocation which hasn't happened yet
  call allocate_rest()
! with memory allocated, get volume element for PC
  if (pccalc.le.nsim)  call get_voli_pc()
!
#ifdef ENABLE_THREADS
! setup work for threads
  call threads_init()
  call threads_setup_molecules()
  call threads_sanity_checks()
  if ((thrdat%test.EQV..true.).AND.(do_restart.EQV..false.)) call threads_test() ! this will terminate execution
#endif
!
! set local counter to zero
  ndump = 0
!
! switch back for restart purposes
  if (do_restart.EQV..true.) then
    call System_Clock(t2)
    eit = eit + 1.0*(t2 - t1)
    call System_Clock(t1)
!   grid can only be initialized properly after reading restart
!   so use residue-based cut in first force/energy calculation
    if (use_cutoffs.EQV..true.) then
      logdummy = use_mcgrid
      logd2 = use_rescrit
      use_mcgrid = .false.
      use_rescrit = .true.   
    end if
    call read_restart()
    if (use_cutoffs.EQV..true.) then
      use_mcgrid = logdummy
      use_rescrit = logd2
    end if
    if (use_mcgrid.EQV..true.) then
      call allocate_mcgrid(aone)
      call setupmcgrid()
      if (grid%report.EQV..true.) then
        call gridreport()
      end if
!     recalculate neighbor lists with grid setup
      call all_respairs_nbl(skip_frz,is_tab,azero)
    end if
#ifdef ENABLE_THREADS
!   may have to redo setup work for threads
    if (use_mcgrid.EQV..true.) call threads_init()
#endif
    call System_Clock(t2)
    time_energy = time_energy + 1.0*(t2 - t1)
    call System_Clock(t1)
  end if
  if (pccalc.le.nsim) call get_voli_pc()
!
  firststep = nstep + 1
!
! deal with hybrid method cycle counters
  if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    if (nstep.lt.nsim) then
      if (in_dyncyc.EQV..false.) then
        mvcnt%nmcseg = mvcnt%nmcseg + 1
        mvcnt%avgmcseglen = mvcnt%avgmcseglen + 1.0*(curcyc_end-max(nstep+1,curcyc_start)+1)
      else
        mvcnt%ndynseg = mvcnt%ndynseg + 1
        mvcnt%avgdynseglen = mvcnt%avgdynseglen + 1.0*(curcyc_end-max(nstep+1,curcyc_start)+1)
      end if
    end if
  end if
!
! print summary before proceeding
  if ((pdb_analyze.EQV..true.).AND.(select_frames.EQV..true.).AND.(curframe.gt.0)) then
    nsim = curframe
    call summary()
    nsim = framecnt
    curframe = 0
  else
    call summary()
    curframe = 0
  end if
!
! finalize time counter for initial setup
!
  call System_Clock(t2)
  eit = eit + 1.0*(t2 - t1)
!
#ifdef ENABLE_THREADS
  call omp_set_num_threads(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(istep,tpi,tpn,rs,logdummy,modstep)
  tpn = omp_get_num_threads()
  tpi = omp_get_thread_num() + 1
  thr_timings(:,tpi) = 0
  thr_timings(15,tpi) = firststep-1
!$OMP SINGLE
  curframe = 0
!$OMP END SINGLE
#else
  thr_timings(:,1) = 0
  thr_timings(15,1) = firststep-1
  curframe = 0
  tpi = 0
#endif
! thread-based grid NBL have different distribution -> recalculate
  if ((use_mcgrid.EQV..true.).AND.(do_restart.EQV..true.).AND.(use_dyn.EQV..true.).AND.(use_cutoffs.EQV..true.).AND.&
 &    (ideal_run.EQV..false.)) then
#ifdef ENABLE_THREADS
    rs_nbl(thr_limits(3,tpi):thr_limits(4,tpi))%ntmpanb = 0
!$OMP BARRIER     
#else
    rs_nbl(:)%ntmpanb = 0
#endif
    call all_respairs_nbl(skip_frz,is_tab,tpi)
  end if
#ifdef ENABLE_THREADS
  if (tpi.le.1) call System_Clock(thr_timings(7,tpi))
#else
  call System_Clock(thr_timings(7,1))
#endif
! a special case is given by analyses "jumping" through trajectories with a frames file
  if ((select_frames.EQV..true.).AND.(pdb_analyze.EQV..true.).AND.((pdb_fileformat.eq.5).OR.(pdb_fileformat.eq.2))) then
!   adjust simulation length for mcstat
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    nsim = framecnt
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
    thr_timings(15,tpi) = 0
    do rs=1,framecnt
#ifdef ENABLE_THREADS
      if (tpi.eq.1) then
#else
      tpi = 1
#endif
      call System_Clock(thr_timings(8,tpi))
      if (1.0*(thr_timings(8,tpi)-thr_timings(7,tpi))/thrdat%rate.gt.60.0*time_flush) then
        thr_timings(16,tpi) = rs
        tt = thrdat%rate*86.4e3*1.0*(thr_timings(16,tpi)-thr_timings(15,tpi))/&
 &      (1.0*(thr_timings(8,tpi)-thr_timings(7,tpi))) ! in steps per day
        dt(1) = floor(1.0*(nsim-rs+1)/tt)
        dt(2) = floor(24.0*(nsim-rs+1)/tt)
        dt(3) = floor(1440.0*(nsim-rs+1)/tt)
        call flushopen()
        thr_timings(15,tpi) = rs
        thr_timings(7,tpi) = thr_timings(8,tpi)
#ifdef ENABLE_THREADS
        if (thrdat%verbosity.gt.0) then
          write(ithread,455) tt
          write(ithread,456) dt(1),dt(2)-dt(1)*24,dt(3)-dt(2)*60
        end if
#else
        write(ilog,455) tt
        write(ilog,456) dt(1),dt(2)-dt(1)*24,dt(3)-dt(2)*60
#endif
      end if
#ifdef ENABLE_THREADS
      end if
#else
      tpi = 0
#endif
      istep = framelst(rs)
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      nstep = nstep + 1
      curframe = rs
      if (pdb_fileformat.eq.2) call FMSMC_readpdb2(istep)
#ifdef LINK_NETCDF
      if (pdb_fileformat.eq.5) call FMSMC_readnetcdf(istep)
#endif
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      if (tpi.le.1) call System_Clock(t1)
      if (use_dyn.EQV..true.) then
#ifdef ENABLE_THREADS
        if (tpi.gt.0) then
          call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
        else
          esave = force3(esterms,esterms_tr,esterms_lr,atrue)
        end if
#else
        esave = force3(esterms,esterms_tr,esterms_lr,atrue)
#endif
      else
        call energy3(esterms,atrue,esave,tpi)
      end if
      if (tpi.le.1) then
        call System_Clock(t2)
        time_energy = time_energy + 1.0*(t2 - t1)
      end if
      if ((align%yes.EQV..true.).AND.(mod(rs,align%calc).eq.0)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call struct_align()
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      end if
      if (tpi.le.1) call System_Clock(t1)
      call mcstat(rs,ndump,tpi)
      if (tpi.le.1) then
        call System_Clock(t2)
        time_analysis = time_analysis + 1.0*(t2 - t1)
      end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
    end do
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
    firststep = nsim+1 ! jump main loop
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
  end if
!
 455 format('--- Current step progression rate: ',g12.4,' steps per day ---')
 456 format('--- Expected time to completion  :',i4,' days : ',i2,' hrs : ',i2,' min ---')
! master loop through individual simulation or trajectory steps
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  curframe = 0
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
  do istep = firststep,nsim
#ifdef ENABLE_THREADS
    if (tpi.eq.1) then
#else
    tpi = 1
#endif
    call System_Clock(thr_timings(8,tpi))
    if (1.0*(thr_timings(8,tpi)-thr_timings(7,tpi))/thrdat%rate.gt.60.0*time_flush) then
      thr_timings(16,tpi) = istep
      tt = thrdat%rate*86.4e3*1.0*(thr_timings(16,tpi)-thr_timings(15,tpi))/&
 &    (1.0*(thr_timings(8,tpi)-thr_timings(7,tpi))) ! in steps per day
      dt(1) = floor(1.0*(nsim-istep+1)/tt)
      dt(2) = floor(24.0*(nsim-istep+1)/tt)
      dt(3) = floor(1440.0*(nsim-istep+1)/tt)
      call flushopen()
      thr_timings(15,tpi) = istep
      thr_timings(7,tpi) = thr_timings(8,tpi)
#ifdef ENABLE_THREADS
      if (thrdat%verbosity.gt.0) then
        write(ithread,455) tt
        write(ithread,456) dt(1),dt(2)-dt(1)*24,dt(3)-dt(2)*60
      end if
#else
      write(ilog,455) tt
      write(ilog,456) dt(1),dt(2)-dt(1)*24,dt(3)-dt(2)*60
#endif
      end if
#ifdef ENABLE_THREADS
    end if
#else
    tpi = 0
#endif
!
!   in the pdb-analysis mode a "move" is the next structure: this is always sequential for all file formats
!   except NetCDF trajectories and individual PDB files with frames files handled above (random access frames file)
    if (pdb_analyze.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
      nstep = nstep + 1
      lstep = istep
      if (pdb_fileformat.eq.1) call FMSMC_readpdb2(lstep)
      if (pdb_fileformat.eq.2) call FMSMC_readpdb2(lstep)
#ifdef LINK_XDR
      if (pdb_fileformat.eq.3) call FMSMC_readxtc(lstep)
#endif
      if (pdb_fileformat.eq.4) call FMSMC_readdcd(lstep)
#ifdef LINK_NETCDF
      if (pdb_fileformat.eq.5) call FMSMC_readnetcdf(lstep)
#endif
      if (lstep.le.0) then
#ifdef ENABLE_MPI
        write(ilog,*) 'Fatal. In MPI-parallel analysis tasks, no failed read-in from a trajectory is tolerated. Check &
 &input files and setting for FMCSC_NRSTEPS.'
        call fexit()
#endif
        write(ilog,'(1x,a,i10,a)') 'Warning. Because of failed structure read-in, the run stops at step #',nstep-1,'.'
        if ((istep-1).ge.nequil) then
          write(ilog,'(1x,a)') 'Built-in analysis routines with a calculation/collection frequency larger than 1 may return &
 &unreasonable results or even cause CAMPARI to crash.'
        end if
        xyzout = nsim + 1
        enout = nsim + 1
        ensout = nsim + 1
        accout = nsim + 1
        torout = nsim + 1
        polout = nsim + 1
        phout = nsim + 1
        savreq%instfreq = nsim + 1
        inst_dssp = .false.
        inst_gpc = 0
        nstep = nsim + 1
        nsim = istep - 1
      end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      if (nstep.ne.istep) then
        lstep = nsim
        if (tpi.le.1) call System_Clock(t1)
        call mcstat(lstep,ndump,tpi)
        if (tpi.le.1) then
          call System_Clock(t2)
          time_analysis = time_analysis + 1.0*(t2 - t1)
        end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        nsim = nstep - 1
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        exit
      end if
      if (select_frames.EQV..true.) then
        if (istep.eq.framelst(curframe+1)) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
          curframe = curframe + 1
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
        else
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
          cycle
        end if
      end if
#ifdef ENABLE_MPI
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      if (re_aux(3).eq.1) call MPI_RETrajMode(istep) ! un-REX trajectories
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#endif
      if (tpi.le.1) call System_Clock(t1)
      if (use_dyn.EQV..true.) then
#ifdef ENABLE_THREADS
        if (tpi.gt.0) then
          call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
        else
          esave = force3(esterms,esterms_tr,esterms_lr,atrue)
        end if
#else
        esave = force3(esterms,esterms_tr,esterms_lr,atrue)
#endif
      else
        call energy3(esterms,atrue,esave,tpi)
      end if
      if (tpi.le.1) then
        call System_Clock(t2)
        time_energy = time_energy + 1.0*(t2 - t1)
      end if
      if ((align%yes.EQV..true.).AND.(mod(istep,align%calc).eq.0)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call struct_align()
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
      end if
      if (tpi.le.1) call System_Clock(t1)
      call mcstat(istep,ndump,tpi)
      if (tpi.le.1) then
        call System_Clock(t2)
        time_analysis = time_analysis + 1.0*(t2 - t1)
      end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      cycle
    end if
!
    if (use_dyn.EQV..true.) then
!
      if (dyn_mode.eq.2) then
#ifdef ENABLE_MPI
        if (use_REMC.EQV..true.) then
          modstep = mod(istep,re_freq)
          if (modstep.eq.0) then
            call MPI_REMaster(istep,ndump,tpi)
            cycle
          end if
        end if
#endif
        if (fycxyz.eq.1) then
#ifdef ENABLE_THREADS
          call int_mdmove_threads(istep,ndump,afalse,tpi)
#else
          call int_mdmove(istep,ndump,afalse)
#endif
        else
#ifdef ENABLE_THREADS
          call cart_mdmove_threads(istep,ndump,afalse,tpi)
#else
          call cart_mdmove(istep,ndump,afalse)
#endif
        end if
!
      else if (dyn_mode.eq.3) then
#ifdef ENABLE_MPI
        if (use_REMC.EQV..true.) then
          modstep = mod(istep,re_freq)
          if (modstep.eq.0) then
            call MPI_REMaster(istep,ndump,tpi)
            cycle
          end if
        end if
#endif
        if (fycxyz.eq.1) then
          call int_ldmove(istep,ndump,afalse)
        else
#ifdef ENABLE_THREADS
          call cart_ldmove_threads(istep,ndump,afalse,tpi)
#else
          call cart_ldmove(istep,ndump,afalse)
#endif
        end if
!
      else if (dyn_mode.eq.6) then
!       this is a bit ugly, as the minimizers internalize the loop
!       hence, exit immediately when finished 
        if ((mini_mode.gt.0).AND.(mini_mode.lt.4)) then
          call minimize(nsim,mini_econv,mini_stepsize,mini_mode,tpi)
          exit
        else if (mini_mode.eq.4) then
          call stochastic_min(nsim,tpi)
          exit
        end if
!
!     hybrid method requires special adjustment
      else if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7)) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!       if new cycle, switch
!       note that RE-moves should never be placed on curcyc_start or curcyc_end
        if ((istep.eq.curcyc_start).AND.(istep.ne.firststep)) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
          if (in_dyncyc.EQV..true.) then
            in_dyncyc = .false.
          else
            in_dyncyc = .true.
          end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
        end if
#ifdef ENABLE_MPI
        if (use_REMC.EQV..true.) then
          modstep = mod(istep,re_freq)
          if (modstep.eq.0) then
            call MPI_REMaster(istep,ndump,tpi)
            if (in_dyncyc.EQV..false.) then
              if (mod(istep,nsancheck).eq.0) then
                if (tpi.le.1) call System_Clock(t1)
                call energy_check(tpi)
                if (tpi.le.1) then
                  call System_Clock(t2)
                  time_energy = time_energy + 1.0*(t2 - t1)
                end if
              end if
            end if
            cycle
          end if
        end if
#endif
        if (in_dyncyc.EQV..true.) then
          logdummy = .false.
          if (istep.eq.curcyc_start) then
            logdummy = .true.
            if ((fycxyz.eq.2).AND.(cart_cons_mode.gt.1).AND.(istep.gt.(first_mccyclen+1))) logdummy = .false.
          end if
          if (dyn_mode.eq.5) then
            if (fycxyz.eq.1) then
#ifdef ENABLE_THREADS
              call int_mdmove_threads(istep,ndump,logdummy,tpi)
#else
              call int_mdmove(istep,ndump,logdummy)
#endif
            else
#ifdef ENABLE_THREADS
              call cart_mdmove_threads(istep,ndump,logdummy,tpi)
#else
              call cart_mdmove(istep,ndump,logdummy)
#endif
            end if
          else if (dyn_mode.eq.7) then
            if (fycxyz.eq.1) then
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
              call int_ldmove(istep,ndump,logdummy)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
            else
#ifdef ENABLE_THREADS
              call cart_ldmove_threads(istep,ndump,logdummy,tpi)
#else
              call cart_ldmove(istep,ndump,logdummy)
#endif
            end if
          end if
          if (istep.eq.curcyc_end) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
            curcyc_start = istep + 1
            curcyc_end = istep + min_mccyclen + floor(random()*(max_mccyclen-min_mccyclen) + 0.5)
            curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
#ifdef ENABLE_MPI
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
            if (use_REMC.EQV..true.) then
!             important: re_freq is never less than 2
              if (mod(curcyc_start,re_freq).eq.0) then
                curcyc_start = curcyc_start + 1
              end if
              if (mod(curcyc_end,re_freq).eq.0) then
                curcyc_end = curcyc_end + 1
              end if
!             publish the new cycle to all nodes and receive on slaves (note that this will override
!             the cycle settings just obtained on the slave node)
              call MPI_SyncHybridCycle(istep)
            else if (use_MPIAVG.EQV..true.) then
              call MPI_SyncHybridCycle(istep)
            end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#endif
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
            if (istep.lt.nsim) then
              mvcnt%nmcseg = mvcnt%nmcseg + 1
              mvcnt%avgmcseglen = mvcnt%avgmcseglen + 1.0*(curcyc_end-curcyc_start+1)
            end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
          end if
        else
!         re-set energy to MC
          if (istep.eq.curcyc_start) then
            logdummy = .true.
!           this corrects for image mismatches due to comm / com - differences in PBC
            call update_rigid_mc_all(tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
            call energy3(esterms,logdummy,esave,tpi)
          end if
!         pick an elementary, non-RE move type we want
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
          call select_mcmove_tree()
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!         displace, evaluate and deal with it
          call mcmove(istep,tpi)
!
          if (mod(istep,nsancheck).eq.0) then
            if (tpi.le.1) call System_Clock(t1)
            call energy_check(tpi)
            if (tpi.le.1) then
              call System_Clock(t2)
              time_energy = time_energy + 1.0*(t2 - t1)
            end if
          end if
          if (istep.eq.curcyc_end) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
            curcyc_start = istep + 1
            curcyc_end = istep + min_dyncyclen + floor(random()*(max_dyncyclen-min_dyncyclen) + 0.5)
            curcyc_end = min(curcyc_end,nsim)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
#ifdef ENABLE_MPI
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
            if (use_REMC.EQV..true.) then
!             important: re_freq is never less than 2
              if (mod(curcyc_start,re_freq).eq.0) then
                curcyc_start = curcyc_start + 1
              end if
              if (mod(curcyc_end,re_freq).eq.0) then
                curcyc_end = curcyc_end + 1
              end if
!             publish the new cycle to all nodes and receive on slaves (note that this will override
!             the cycle settings just obtained on the slave node)
              call MPI_SyncHybridCycle(istep)
            else if (use_MPIAVG.EQV..true.) then
              call MPI_SyncHybridCycle(istep)
            end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#endif
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
            if (istep.lt.nsim) then
              mvcnt%ndynseg = mvcnt%ndynseg + 1
              mvcnt%avgdynseglen = mvcnt%avgdynseglen + 1.0*(curcyc_end-curcyc_start+1)
            end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
          end if
        end if
!
      end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      call threads_dlb(istep)
!$OMP END SINGLE
#endif
      if (tpi.le.1) call System_Clock(t1)
      call mcstat(istep,ndump,tpi)
      if (tpi.le.1) then
        call System_Clock(t2)
        time_analysis = time_analysis + 1.0*(t2 - t1)
      end if
!
#ifdef ENABLE_MPI
      if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.).AND.(istep.gt.1)) then
        modstep = mod(istep,re_freq)
        if ((modstep.eq.0).AND.(istep.lt.nsim)) then
          call MPI_ASMaster(istep,tpi) ! does not increment step -> do not cycle
        end if
      end if
#endif
!
    else
!

!   determine what movetype we want
#ifdef ENABLE_MPI
    if (use_REMC.EQV..true.) then
      modstep = mod(istep,re_freq)
      if (modstep.eq.0) then
        call MPI_REMaster(istep,ndump,tpi)
        if (mod(istep,nsancheck).eq.0) then
          if (tpi.le.1) call System_Clock(t1)
          call energy_check(tpi)
          if (tpi.le.1) then
            call System_Clock(t2)
            time_energy = time_energy + 1.0*(t2 - t1)
          end if
        end if
        cycle
      end if
    end if
#endif
!
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
!   pick an elementary, non-RE move type we want
    call select_mcmove_tree()
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!   displace, evaluate and deal with it
    call mcmove(istep,tpi)
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
      call threads_dlb(istep)
!$OMP END SINGLE
#endif
!
!   write out current statistics
    if (tpi.le.1) call System_Clock(t1)
    call mcstat(istep,ndump,tpi)
    if (tpi.le.1) then
      call System_Clock(t2)
      time_analysis = time_analysis + 1.0*(t2 - t1)
    end if
!
#ifdef ENABLE_MPI
    if ((use_MPIAVG.EQV..true.).AND.(use_MPIMultiSeed.EQV..true.).AND.(istep.gt.1)) then
      modstep = mod(istep,re_freq)
      if ((modstep.eq.0).AND.(istep.lt.nsim)) then
        call MPI_ASMaster(istep,tpi) ! does not increment step -> do not cycle
      end if
    end if
#endif
!
    if (mod(istep,nsancheck).eq.0) then
      if (tpi.le.1) call System_Clock(t1)
      call energy_check(tpi)
      if (tpi.le.1) then
        call System_Clock(t2)
        time_energy = time_energy + 1.0*(t2 - t1)
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
    end if
!
  end do
!
! in this mode of operation, we have skipped both prior loops
  if ((nsim.le.0).AND.(cstored.gt.0).AND.(cmode.eq.7)) then
    call file_clustering(cmode,nstruccls,tpi)
  end if
!
#ifdef ENABLE_THREADS
!$OMP END PARALLEL
#endif
!
  intfile = 'END '
  call FMCSC_dump(intfile)
!
! if in pdb-analysis mode, we might have to close the trajectory file
  if (pdb_analyze.EQV..true.) then
    if ((nsim.le.0).AND.(cstored.gt.0).AND.(cmode.eq.7)) then
!     do nothing
    else
      if (pdb_fileformat.eq.1) call close_pdbtraj()
#ifdef LINK_XDR
      if (pdb_fileformat.eq.3) call close_xtctraj()
#endif
      if (pdb_fileformat.eq.4) call close_dcdtraj()
#ifdef LINK_NETCDF
      if (pdb_fileformat.eq.5) call close_netcdftraj()
#endif
    end if
  end if
!
! finalize I/O (close files)
  call makeio(2)
  
!
! summarize performance information
  call System_Clock(ee2)
  ee = 1.0*(ee2 - ee1)
  write(ilog,*)
  write(ilog,*) 'Total CPU time elapsed [s]     : ',ee/thrdat%rate
  write(ilog,*) 'Fraction for energy functions  : ',100.0*time_energy/ee,'%'
  write(ilog,*) 'Fraction for analysis routines : ',100.0*time_analysis/ee,'%'
  if (time_struc.gt.0.0) write(ilog,*) 'Fraction for chain closure     : ',100.0*time_struc/ee,'%'
  if (time_holo.gt.0.0) write(ilog,*) 'Fraction for constraint solvers: ',100.0*time_holo/ee,'%'
  if (time_ph.gt.0.0) write(ilog,*) 'Fraction for pH routines       : ',100.0*time_ph/ee,'%'
  if (time_pme.gt.0.0) write(ilog,*) 'Fraction for PME reciprocal sum: ',100.0*time_pme/ee,'%'
  if (time_nbl.gt.0.0) write(ilog,*) 'Fraction for neighbor lists    : ',100.0*time_nbl/ee,'%'
  write(ilog,*) 'Fraction for initial setup     : ',100.0*eit/ee,'%'
#ifdef ENABLE_MPI
  write(ilog,*) 'Fraction for communication     : ',100.0*time_comm/ee,'%'
  write(ilog,*) 'Fraction of remainder          : ',100.0*(1.0 - time_energy/ee - time_analysis/ee - time_struc/ee -&
 &time_comm/ee - time_ph/ee - time_holo/ee - time_nbl/ee - time_pme/ee - eit/ee),'%'
#else
  write(ilog,*) 'Fraction of remainder          : ',100.0*(1.0 - time_energy/ee - time_analysis/ee - time_struc/ee -&
 & time_ph/ee - time_holo/ee - time_nbl/ee - time_pme/ee - eit/ee),'%'
#endif
  write(ilog,*)
!
! and clean-up of course ....
  call deallocate_all()
  deallocate(blas)
!
  call Date_and_Time(rc(1),rc(2),rc(3),dt)
  write(ilog,*) 'Execution ended at ',rc(2)(1:2),':',rc(2)(3:4),&
 &' on ',rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
!
#ifdef ENABLE_MPI
  call makelogio(2)
#endif
!
! in parallel mode, need to shutdown MPI universe
#ifdef ENABLE_MPI
  call MPI_StopMC()
  call MPI_FINALIZE(ierr)
#endif
!
end
!
!------------------------------------------------------------------------------
!
subroutine write_license()
!
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
  write(*,*) '                  ---        Welcome to CAMPARI       ---'
  write(*,*)
  write(*,*) '                  ---       Released Version 3.0b     ---'
  write(*,*)
  write(*,*) '    Copyright (C) 2017, CAMPARI Development Team'
  write(*,*)
  write(*,*) '    Websites: http://campari.sourceforge.net'
  write(*,*) '              http://sourceforge.net/projects/campari/'
  write(*,*) '              http://pappulab.wustl.edu/campari'
  write(*,*)
  write(*,*) '    CAMPARI is free software: you can redistribute it and/or modify'
  write(*,*) '    it under the terms of the GNU General Public License as published by'
  write(*,*) '    the Free Software Foundation, either version 3 of the License, or'
  write(*,*) '    (at your option) any later version.'
  write(*,*)
  write(*,*) '    CAMPARI is distributed in the hope that it will be useful,'
  write(*,*) '    but WITHOUT ANY WARRANTY; without even the implied warranty of'
  write(*,*) '    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
  write(*,*) '    GNU General Public License for more details. '
  write(*,*)
  write(*,*) '    You should have received a copy of the GNU General Public License' 
  write(*,*) '    along with CAMPARI. If not, see <http://www.gnu.org/licenses/>.'
  write(*,*)
  write(*,*)
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
  write(*,*)
!
end
!
!------------------------------------------------------------------------------
!

