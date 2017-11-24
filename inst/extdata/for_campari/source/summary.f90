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
! CONTRIBUTIONS: Marco Bacci, Rohit Pappu, Hoang Tran, Adam Steffen        !
!                Albert Mao                                                !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!------------------------------------------------------------------------------
!
! summarize parameters picked and simulation to be performed 
!
subroutine summary()
!
  use iounit
  use sequen
  use polypep
  use energies
  use system
  use cutoffs
  use contacts
  use inter
  use torsn
  use params
  use grids
  use mcgrid
  use pdb
  use ionize
  use polyavg
  use movesets
  use paircorr
  use mpistuff
  use tabpot
  use distrest
  use mcsums
  use molecule
  use dipolavg
  use diffrac
  use atoms
  use forces
  use ewalds
  use mini
  use grandensembles
  use dssps
  use ujglobals
  use math
  use threads
  use clusters
  use shakeetal
  use ems
  use fos
  use wl
!
  implicit none
!
  integer ti,ti2,ti3,i,j,s3,s4,s5,s6,dummy,t1,t2,ft1,ft2,ff1,ff2,tailscnt,sccnt,atomscnt,pdbrsize
  character(MAXSTRLEN) fn,fn2
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
! 
  5   format(a,1x,a)
  8   format(a,1x,i8,1x,i8,1x,i8)
  9   format(a,1x,f8.4,1x,f8.4,1x,f8.4)
  10  format(a,8x,f9.4)
  11  format(a,5x,f9.1)
  20  format(a,2x,f8.3)
  22  format(a,2x,f8.3,1x,f8.3)
  21  format(a,2x,i8)
  23  format(a,f9.3,a,f9.3,a,f9.3)
  24  format(a,i4,a,i4,a,i4)
  25  format(a,9x,f7.3)
  26  format(a,9x,f7.3)
  27  format(a,2x,f8.3,a)
  28  format(a,1x,g9.3)
  49  format(1x,a,i10)
  50  format(a,7x,f8.4)
  51  format(2x,a,i3,a,f9.3)
  52  format(2x,a,i3,a,i5)
  53  format(1x,a,i5,a,i5,a)
  54  format(1x,a,i6)
  55  format(a,10x,f6.3)
  56  format(1x,a,i8)
  70  format(a,5x,f13.3)
  71  format(a,1x,g14.7)
  72  format(a,5x,g13.3)
  73  format(a,5x,g13.6)
  74  format(a,3x,g13.3)
  101 format(a,8x,g9.4)
  102 format(a,8x,g9.3)
!
  write(ilog,*)
  if (do_restart.EQV..true.) then
    write(ilog,*) '-- THIS IS A RESTART OF A PREVIOUSLY ABORTED --'
    write(ilog,*) '-- RUN AND ON-THE-FLY ANALYSIS WILL ONLY INC --'
    write(ilog,*) '-- LUDE INFORMATION OF THE RESTARTED PORTION --'
  else
    write(ilog,*) '--- SUMMARY OF CALCULATION ---'
  end if
  write(ilog,*)
  write(ilog,*) '**IMPORTANT FILES (IF EXISTENT):'
  call strlims(paramfile,s3,s4)
  write(ilog,5) ' Parameter file used   : ',paramfile(s3:s4)
  call strlims(seqfile,s3,s4)
  write(ilog,5) ' Sequence file used    : ',seqfile(s3:s4)
  call strlims(angrpfile,s3,s4)
  write(ilog,5) ' Analy. grp. file used : ',angrpfile(s3:s4)
  if (segcalc.le.nsim) then
    call strlims(bbsegfile,s3,s4)
    write(ilog,5) ' BB-segment file used  : ',bbsegfile(s3:s4)
  end if
  if (mass_patched.EQV..true.) then
    call strlims(masspatchfile,s3,s4)
    write(ilog,5) ' Masses patch file used: ',masspatchfile(s3:s4)
  end if
  if (rad_patched.EQV..true.) then
    call strlims(radpatchfile,s3,s4)
    write(ilog,5) ' Radii patch file used : ',radpatchfile(s3:s4)
  end if
  if (bonded_patched.EQV..true.) then
    call strlims(bpatchfile,s3,s4)
    write(ilog,5) ' Bonded patch file used: ',bpatchfile(s3:s4)
  end if
  if (lj_patched.EQV..true.) then
    call strlims(ljpatchfile,s3,s4)
    write(ilog,5) ' LJ patch file used    : ',ljpatchfile(s3:s4)
  end if
  if (charge_patched.EQV..true.) then
    call strlims(cpatchfile,s3,s4)
    write(ilog,5) ' Charge patch file used: ',cpatchfile(s3:s4)
  end if
  if (nc_patched.EQV..true.) then
    call strlims(ncpatchfile,s3,s4)
    write(ilog,5) ' q-Flag patch file used: ',ncpatchfile(s3:s4)
  end if
  if (fos_patched.EQV..true.) then
    call strlims(fospatchfile,s3,s4)
    write(ilog,5) ' FOS patch file used   : ',fospatchfile(s3:s4)
  end if
  if (ard_patched.EQV..true.) then
    call strlims(ardpatchfile,s3,s4)
    write(ilog,5) ' Vred patch file used  : ',ardpatchfile(s3:s4)
  end if
  if (asm_patched.EQV..true.) then
    call strlims(asmpatchfile,s3,s4)
    write(ilog,5) ' AtSAV patch file used : ',asmpatchfile(s3:s4)
  end if
  call strlims(basename,s3,s4)
  write(ilog,5) ' Basename for output   : ',basename(s3:s4)
  call strlims(pdb_formstr(1),ft1,ft2)
  write(ilog,*) 'Output string for PDB :  ',pdb_formstr(1)(ft1:ft2)
  if (((use_pdb_template.EQV..true.).OR.(pdbinput.EQV..true.)).AND.(pdb_analyze.EQV..false.)) then
    call strlims(pdbinfile,s3,s4)
    call strlims(pdbtmplfile,s5,s6)
#ifdef ENABLE_MPI
    if (pdb_mpimany.EQV..true.) then
      write(ilog,5) ' Base filename for PDBs: ',pdbinfile(s3:s4)
    else if (pdbinput.EQV..true.) then
      write(ilog,5) ' Structural input file : ',pdbinfile(s3:s4)
    end if
    if (use_pdb_template.EQV..true.) write(ilog,5) ' PDB template file     : ',pdbtmplfile(s5:s6)
#else
    if (pdbinput.EQV..true.) write(ilog,5) ' Structural input file : ',pdbinfile(s3:s4)
    if (use_pdb_template.EQV..true.) write(ilog,5) ' PDB template file     : ',pdbtmplfile(s5:s6)
#endif
    write(ilog,*) 'Reading mode          : ',pdb_readmode
    write(ilog,*) 'Read-in format conv.  : ',pdb_convention(2)
    if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
      write(ilog,*) 'Assuming that molecules in input files are broken to fit into box.'
    else
      write(ilog,*) 'Assuming that molecules in input files are always intact.'
    end if
    write(ilog,*) 
    call strlims(pdb_formstr(2),ft1,ft2)
    write(ilog,*) 'Input string for PDB  :  ',pdb_formstr(2)(ft1:ft2)
    if (pdbinput.EQV..true.) then
      write(ilog,*) 'Hydrogen code         : ',pdb_hmode
      write(ilog,27) ' Bond length tol. low  : ',100.0*pdb_tolerance(1),' %'
      write(ilog,27) ' Bond length tol. high : ',100.0*pdb_tolerance(2),' %'
      write(ilog,27) ' Bond angle tolerance  : ',pdb_tolerance(3),' deg.'
    end if
  end if
#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call strlims(re_infile,s3,s4)
    write(ilog,5) ' Input file for MPI-RE : ',re_infile(s3:s4)
  end if
#endif
  if (use_TOR.EQV..true.) then
    call strlims(torfile,s3,s4)
    write(ilog,5) ' Torsion def. file used: ',torfile(s3:s4)
  end if
  if (use_POLY.EQV..true.) then
    call strlims(polfile,s3,s4)
    write(ilog,5) ' Polym. def. file used : ',polfile(s3:s4)
  end if
  if (use_TABUL.EQV..true.) then
    call strlims(tabcodefile,s3,s4)
    write(ilog,5) ' Atom code matrix file : ',tabcodefile(s3:s4)
    call strlims(tabpotfile,s3,s4)
    write(ilog,5) ' Tabulated pot. file   : ',tabpotfile(s3:s4)
    if (tbp%tangmode.gt.0) then
      call strlims(tabtangfile,s3,s4)
      write(ilog,5) ' Tabulated tang. file  : ',tabtangfile(s3:s4)
    end if
  end if
  if (use_DREST.EQV..true.) then
    call strlims(drestfile,s3,s4)
    write(ilog,5) ' Dist. rest. file used : ',drestfile(s3:s4)
  end if
  if (use_OSMO.EQV..true.) then
    call strlims(osmofile,s3,s4)
    write(ilog,5) ' Compartment file used : ',osmofile(s3:s4)
  end if
  if (do_frz.EQV..true.) then
    call strlims(frzfile,s3,s4)
    write(ilog,5) ' Constraint file used  : ',frzfile(s3:s4)
    if (use_dyn.EQV..true.) then
      if (skip_frz.EQV..true.) then
        write(ilog,*) 'Constrained interactions are ignored in force calculations.'
      else
        write(ilog,*) 'Constrained interactions are included in force calculations.'
      end if
    end if
  end if
  if (pdb_analyze.EQV..false.) then
    call strlims(prefsamplingfile,s3,s4)
    write(ilog,5) ' PSW file used (if any): ',prefsamplingfile(s3:s4)
  end if
  if (use_FEG.EQV..true.) then
    call strlims(fegfile,s3,s4)
    write(ilog,5) ' Ghosting file used    : ',fegfile(s3:s4)
  end if
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    call strlims(particleflucfile,s3,s4)
    write(ilog,5) ' Input file for (S)GCE : ',particleflucfile(s3:s4)
  end if
!  if (use_LCTOR.EQV..true.) then
!    call strlims(torlcfile,s3,s4)
!    write(ilog,5) ' LCT coefficients file : ',torlcfile(s3:s4)
!    if (use_lctmoves.EQV..true.) then
!      write(ilog,5) ' File for inv. trans.  : ',torlcfile2
!    end if
!    write(ilog,5) ' Tab. LCT pot. file    : ',lctpotfile
!  else if (torlccalc.le.nsim) then
!    write(ilog,5) ' LCT coefficients file : ',torlcfile
!    if (use_lctmoves.EQV..true.) then
!      write(ilog,5) ' File for inv. trans.  : ',torlcfile2
!    end if
!  end if
  write(ilog,*)
  write(ilog,*) '**GENERAL PARAMETERS:'
  if (be_unsafe.EQV..true.) then
    write(ilog,*) 'Warning. This run has some fatal exits in sanity checks disabled.&
 & Double-check "fatal" warnings above.'
  end if
  if (pdb_analyze.EQV..true.) then
    write(ilog,*) 'This is an analysis run for an input trajectory.'
    write(ilog,*) 'Not all warnings and settings reported above and below may be relevant.'
    if ((pdb_fileformat.eq.1).OR.(pdb_fileformat.eq.2)) then
      call strlims(pdbinfile,s3,s4)
      fn2(s3:s4) = pdbinfile(s3:s4)
    else if (pdb_fileformat.eq.3) then
      call strlims(xtcinfile,s3,s4)
      fn2(s3:s4) = xtcinfile(s3:s4)
    else if (pdb_fileformat.eq.4) then
      call strlims(dcdinfile,s3,s4)
      fn2(s3:s4) = dcdinfile(s3:s4)
    else if (pdb_fileformat.eq.5) then
      call strlims(netcdfinfile,s3,s4)
      fn2(s3:s4) = netcdfinfile(s3:s4)
    end if
#ifdef ENABLE_MPI
    do j=s4,s3,-1
      if (fn2(j:j).eq.SLASHCHAR) exit
    end do
    ti = j
    call int2str(myrank,nod,re_aux(10))
    if ((ti.ge.s3).AND.(ti.lt.s4)) then
      fn =  fn2(s3:ti)//'N_'//nod(1:re_aux(10))//'_'//fn2((ti+1):s4)
    else
      fn =  'N_'//nod(1:re_aux(10))//'_'//fn2(s3:s4)
    end if
    call strlims(fn,s3,s4)
#else
    fn = fn2
#endif
    write(ilog,5) ' Base file           : ',fn(s3:s4)
#ifdef LINK_XDR
    if (pdb_fileformat.eq.3) then
      write(ilog,11) ' Precision multiplier: ',xtc_prec
    end if
#endif
    if ((pdb_force_box.eq.1).OR.(pdb_force_box.eq.3)) then
      write(ilog,*) 'Assuming that molecules in input files are broken to fit into box.'
    else
      write(ilog,*) 'Assuming that molecules in input files are always intact.'
    end if
    if (use_pdb_template.EQV..true.) then
      call strlims(pdbtmplfile,s3,s4)
      write(ilog,5) ' Template file       : ',pdbtmplfile(s3:s4)
    end if
    if ((pdb_fileformat.eq.1).OR.(pdb_fileformat.eq.2).OR.(use_pdb_template.EQV..true.)) then
      write(ilog,*) 'Format code         : ',pdb_fileformat
      write(ilog,*) 'Hydrogen code       : ',pdb_hmode
      write(ilog,27) ' Bond length tol. low: ',100.0*pdb_tolerance(1),' %'
      write(ilog,27) ' Bond length tol. hi.: ',100.0*pdb_tolerance(2),' %'
      write(ilog,27) ' Bond angle tolerance: ',pdb_tolerance(3),' deg.'
      write(ilog,*) 'Read-in format conv.: ',pdb_convention(2)
      call strlims(pdb_formstr(2),ft1,ft2)
      write(ilog,*) 'Input string for PDB:  ',pdb_formstr(2)(ft1:ft2)
    end if
    if (select_frames.EQV..true.) then
      write(ilog,*) 'This run will only analyze user-selected frames.'
      if (use_frame_weights.EQV..true.) then
        write(ilog,*) 'All frames receive specific weights provided through the input file.'
      end if
      write(ilog,*) 'Number of frames    : ',framecnt
      call strlims(frameidxfile,s3,s4)
      write(ilog,*) 'Filename used       : ',frameidxfile(s3:s4)
    end if
    write(ilog,*)
  end if
  if (dyn_mode.ne.6) then
    if (pdb_analyze.EQV..false.) then
      write(ilog,*) 'Number of steps     : ',nsim
      if (do_restart.EQV..true.) then
        write(ilog,*) 'Restart from step   : ',nstep+1
      end if
    else  
      write(ilog,*) 'Min. traj. length   : ',nsim
    end if
    write(ilog,*) 'Equilibration steps : ',nequil
  end if
  if (phfreq.gt.0.0) then
     write(ilog,10)' pH                  : ',phs
  end if
  if (bnd_type.eq.1) then
    if (bnd_shape.eq.1) then
      write(ilog,*) '3D periodic boundary conditions (PBC) are used.'
    else if (bnd_shape.eq.3) then
      write(ilog,*) '1D periodic boundary conditions (PBC) are used along cylinder axis.'
      write(ilog,*) 'A soft-wall atom-based boundary condition (SWBC) is used elsewhere.'
      write(ilog,10) ' Force constant      : ',bnd_params(7)
    end if
  else if (bnd_type.eq.2) then
    write(ilog,*) 'A hard-wall boundary condition (HWBC) is used.'
  else if (bnd_type.eq.3) then
    write(ilog,*) 'A soft-wall residue-based boundary condition (SWBC) is used.'
    write(ilog,10) ' Force constant      : ',bnd_params(7)
  else if (bnd_type.eq.4) then
    write(ilog,*) 'A soft-wall atom-based boundary condition (SWBC) is used.'
    write(ilog,10) ' Force constant      : ',bnd_params(7)
  end if
  if (bnd_shape.eq.1) then
    write(ilog,*) 'The simulation system is a rectangular cuboid box.'
    write(ilog,9) ' Box lengths (A)      : ',bnd_params(1),bnd_params(2),bnd_params(3)
    write(ilog,9) ' Box origin (A)       : ',bnd_params(4),bnd_params(5),bnd_params(6)
  else if (bnd_shape.eq.2) then
    write(ilog,*) 'The simulation system is a sphere.'
    write(ilog,9) ' Sphere center (A)   : ',bnd_params(1),bnd_params(2),bnd_params(3)
    write(ilog,9) ' Sphere radius (A)   : ',bnd_params(4)
  else if (bnd_shape.eq.3) then
    write(ilog,*) 'The simulation system is a cylinder (z-axis aligned).'
    write(ilog,9) ' Cylinder center (A) : ',bnd_params(1),bnd_params(2),bnd_params(3)
    write(ilog,9) ' Cylinder radius (A) : ',bnd_params(4)
    write(ilog,9) ' Cylinder height (A) : ',bnd_params(6)
  end if
  if (pdb_analyze.EQV..false.) then
  if (dyn_mode.eq.1) then
    write(ilog,*) 'Simulation uses pure Monte Carlo sampling (see TE&
 &CHNICAL SETTINGS below) in rigid body/torsional space.'
  else if ((dyn_mode.eq.2).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Simulation uses ballistic, Newtonian molecular dy&
 &namics sampling in rigid body/torsional space.'
  else if ((dyn_mode.eq.2).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'Simulation uses ballistic, Newtonian molecular dy&
 &namics sampling in Cartesian space.'
  else if ((dyn_mode.eq.3).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Simulation uses stochastic (Langevin) dynamics sa&
 &mpling in rigid body/torsional space.'
  else if ((dyn_mode.eq.3).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'Simulation uses stochastic (Langevin) dynamics sa&
 &mpling in Cartesian space.'
  else if ((dyn_mode.eq.4).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Simulation uses stochastic (Brownian) dynamics sa&
 &mpling in rigid body/torsional space.'
  else if ((dyn_mode.eq.4).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'Simulation uses stochastic (Brownian) dynamics sa&
 &mpling in Cartesian space.'
  else if ((dyn_mode.eq.5).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Simulation uses mixed Newtonian dynam&
 &ics and Monte Carlo sampling in rigid body/torsional space.'
  else if ((dyn_mode.eq.5).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'Simulation uses mixed Newtonian dynam&
 &ics and Monte Carlo sampling in Cartesian space (MD only).'
  else if ((dyn_mode.eq.7).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Simulation uses mixed stochastic (Langevin) dynam&
 &ics and Monte Carlo sampling in rigid body/torsional space.'
  else if ((dyn_mode.eq.7).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'Simulation uses mixed stochastic (Langevin) dynam&
 &ics and Monte Carlo sampling in Cartesian space (LD only).'
  else if ((dyn_mode.eq.8).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'Simulation uses mixed stochastic (Brownian) dynam&
 &ics and Monte Carlo sampling in rigid body/torsional space.'
  else if ((dyn_mode.eq.8).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'Simulation uses mixed stochastic (Brownian) dynam&
 &ics and Monte Carlo sampling in Cartesian space (BD only).'
  else if ((dyn_mode.eq.6).AND.(fycxyz.eq.1)) then
    write(ilog,*) 'This run is a minimization run in rigid body/tors&
 &ional space.'
 else if ((dyn_mode.eq.6).AND.(fycxyz.eq.2)) then
    write(ilog,*) 'This run is a minimization run in Cartesian space.'
  end if
  if ((dyn_mode.ne.1).AND.(fycxyz.eq.1)) then
    if (align_NC.eq.1) then
      write(ilog,*) 'Torsional dynamics will use N-termini as centers of movement and kinetics will&
      & be slowed down artificially at those N-termini.'
    else if (align_NC.eq.2) then
      write(ilog,*) 'Torsional dynamics will use C-termini as centers of movement and kinetics will&
      & be slowed down artificially at those C-termini.'
    else if ((align_NC.eq.3).OR.(align_NC.eq.4)) then
      write(ilog,*) 'Torsional dynamics will use central portion of molecules as center of movement&
      & and kinetics will be comparable at both termini.'
    end if
    if ((unklst%nr.gt.0).OR.(unslst%nr.gt.0)) then
      if ((dyn_integrator_ops(10).eq.0).OR.((dyn_integrator_ops(10).eq.1).AND.(unklst%nr.eq.0)).OR.&
 &             ((dyn_integrator_ops(10).eq.2).AND.(unslst%nr.eq.0))) then
        write(ilog,*) 'Torsional dynamics will automatically constrain all degrees of freedom that are not native &
 &to CAMPARI (see documentation on FMCSC_SEQFILE).'
      else if ((dyn_integrator_ops(10).eq.3).OR.((dyn_integrator_ops(10).eq.1).AND.(unslst%nr.eq.0)).OR.&
 &             ((dyn_integrator_ops(10).eq.2).AND.(unklst%nr.eq.0))) then
        write(ilog,*) 'Torsional dynamics will consider all available degrees of freedom irrespective of whether they are native &
 &to CAMPARI or not (see documentation on FMCSC_SEQFILE). This is before application of any custom constraints (FMCSC_FRZFILE).'
      else if ((dyn_integrator_ops(10).eq.1).AND.(unslst%nr.gt.0).AND.(unklst%nr.gt.0)) then
        write(ilog,*) 'Torsional dynamics will consider available degrees of freedom in unsupported residues, but not any other &
 &nonnative degrees of freedom (see documentation on FMCSC_SEQFILE). This is before application of any custom constraints &
 &(FMCSC_FRZFILE).' 
      else if ((dyn_integrator_ops(10).eq.2).AND.(unslst%nr.gt.0).AND.(unklst%nr.gt.0)) then
        write(ilog,*) 'Torsional dynamics will consider nonnative degrees of freedom in supported residues, but not any available&
 & degrees of freedom in unsupported residues (see documentation on FMCSC_SEQFILE). This is before application of any custom &
 &constraints (FMCSC_FRZFILE).' 
      end if
    end if
  end if
  if ((dyn_mode.ge.2).AND.(dyn_mode.ne.6)) then
    if (((dyn_mode.eq.2).OR.(dyn_mode.eq.5)).AND.(fycxyz.eq.1)) then
      write(ilog,*) 'TMD integrator type : ',dyn_integrator
      write(ilog,*) 'Velocity updates    : ',dyn_integrator_ops(1)
    end if
    write(ilog,102) ' Integration timestep: ',dyn_dt
  end if
  if ((dyn_mode.ge.3).AND.(dyn_mode.ne.6).AND.(dyn_mode.ne.5)) then
    write(ilog,102) ' Friction Parameter  : ',fric_ga
  end if
  if ((ens%flag.eq.1).AND.((dyn_mode.ne.6).OR.((dyn_mode.eq.6).AND.(mini_mode.eq.4)))) then
    write(ilog,*) 'Simulation is performed in the isochoric-isothermal ensemble (NVT).'
    write(ilog,10) ' Desired temperature : ',kelvin
    if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5).OR.((dyn_mode.eq.6).AND.(mini_mode.eq.4))) then
      if (tstat%flag.eq.1) then
        write(ilog,*) 'Berendsen thermostat used to maintain temperature in Newtonian MD.'
        write(ilog,10) ' Coupling time in ps : ',tstat%params(1)
        write(ilog,*) 'Number of T-groups  : ',tstat%n_tgrps
      else if (tstat%flag.eq.2) then
        write(ilog,*) 'Andersen thermostat used to maintain temperature in Newtonian MD.'
        write(ilog,10) ' Decay time in ps    : ',tstat%params(1)
        write(ilog,*) 'Number of T-groups  : ',tstat%n_tgrps
      else if (tstat%flag.eq.4) then
        write(ilog,*) 'Bussi-Parrinello thermostat used to maintain temperature in Newtonian MD.'
        write(ilog,10) ' Decay time in ps    : ',tstat%params(1)
        write(ilog,*) 'Number of T-groups  : ',tstat%n_tgrps
      end if
    end if
  else if ((ens%flag.eq.2).AND.(dyn_mode.ne.6)) then
    write(ilog,*) 'Simulation is performed in the isochoric-isoenergetic ensemble (NVE).'
    write(ilog,10) ' Initial temperature : ',kelvin
  else if (ens%flag.eq.3) then
    write(ilog,*) 'Simulation is performed in the isobaric-isothermal ensemble (NPT).'
    write(ilog,10) ' Desired temperature : ',kelvin
    if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
      if (tstat%flag.eq.1) then
        write(ilog,*) 'Berendsen thermostat used to maintain temperature in Newtonian MD.'
        write(ilog,10) ' Coupling time in ps : ',tstat%params(1)
        write(ilog,*) 'Number of T-groups  : ',tstat%n_tgrps
      else if (tstat%flag.eq.2) then
        write(ilog,*) 'Andersen thermostat used to maintain temperature in Newtonian MD.'
        write(ilog,10) ' Decay time in ps    : ',tstat%params(1)
        write(ilog,*) 'Number of T-groups  : ',tstat%n_tgrps
      else if (tstat%flag.eq.4) then
        write(ilog,*) 'Bussi-Parrinello thermostat used to maintain temperature in Newtonian MD.'
        write(ilog,10) ' Decay time in ps    : ',tstat%params(1)
        write(ilog,*) 'Number of T-groups  : ',tstat%n_tgrps
      end if
    end if
    write(ilog,10) ' Desired pressure    : ',extpress
    if ((bnd_type.ge.3).AND.(bnd_type.le.4)) then
        write(ilog,*) 'Reflective droplet manostat used to maintain pressure in Newtonian MD.'
      write(ilog,10) ' Mass of boundary    : ',pstat%params(1)
    end if
  else if (ens%flag.eq.4) then
    write(ilog,*) 'Simulation is performed in the isobaric-isenthalpic ensemble (NPE).'
    write(ilog,10) ' Initial temperature : ',kelvin
    write(ilog,10) ' Desired pressure    : ',extpress
    if ((bnd_type.ge.3).AND.(bnd_type.le.4)) then
        write(ilog,*) 'Reflective droplet manostat used to maintain pressure in Newtonian MD.'
      write(ilog,10) ' Mass of boundary    : ',pstat%params(1)
    end if
  else if (ens%flag.eq.5) then
    write(ilog,*) 'Simulation is performed in the grand isochoric-isothermal ensemble (muVT).'
    write(ilog,10) ' Desired temperature : ',kelvin
  else if (ens%flag.eq.6) then
    write(ilog,*) 'Simulation is performed in the semigrand isochoric-isothermal ensemble (DeltamuVT).'
    write(ilog,10) ' Desired temperature : ',kelvin
  end if
  if ((use_dyn.EQV..true.).AND.(fycxyz.eq.2)) then
    if (cart_cons_grps.le.0) then
      write(ilog,*) 'There are no topology-derived or user-requested bond lengths or angles eligible for &
 &holonomic constraints.'
    else if ((add_settle.EQV..false.).OR.((settle_tip3ps.eq.0).AND.(settle_tip4ps.eq.0).AND.(settle_tip5ps.eq.0).AND.&
 &(settle_tip4pes.eq.0).AND.(settle_spcs.eq.0))) then
      if (cart_cons_mode.eq.1) then
        write(ilog,*) 'No constraints to bond lengths or other geometrical parameters are enforced.'
      else if (cart_cons_mode.eq.2) then
        write(ilog,*) 'Constraining all topology-derived bonds to terminal atoms interpreted to be hydrogen.'
      else if (cart_cons_mode.eq.3) then
        write(ilog,*) 'Constraining all topology-derived bond lengths'
      else if (cart_cons_mode.eq.4) then
        write(ilog,*) 'Constraining all topology-derived bond lengths and bond angles (explicit angle constraints).'
      else if (cart_cons_mode.eq.5) then
        write(ilog,*) 'Constraining all topology-derived bond lengths and bond angles (via additional distance constraints).'
      else if (cart_cons_mode.eq.6) then
        write(ilog,*) 'Constraining all bond lengths as defined by user selection.'
        call strlims(cart_cons_file,t1,t2)
        write(ilog,*) 'Constraint set file : ',cart_cons_file(t1:t2)
      end if
    else
      if (cart_cons_mode.eq.1) then
        write(ilog,*) 'Constraining eligible water models (TipXp, SPC) to be rigid.'
      else if (cart_cons_mode.eq.2) then
        write(ilog,*) 'Constraining all topology-derived bonds to terminal atoms interpreted to be hydrogen and &
 &constraining eligible water models (TipXp, SPC) to be rigid.'
      else if (cart_cons_mode.eq.3) then
        write(ilog,*) 'Constraining all topology-derived bond lengths and &
 &constraining eligible water models (TipXp, SPC) to be rigid.'
      else if (cart_cons_mode.eq.4) then
        write(ilog,*) 'Constraining all topology-derived bond lengths and bond angles (explicit angle constraints) and &
 &constraining eligible water models (TipXp, SPC) to be rigid.'
      else if (cart_cons_mode.eq.5) then
        write(ilog,*) 'Constraining all topology-derived bond lengths and bond angles (via additional distance constraints) and &
 &constraining eligible water models (TipXp, SPC) to be rigid.'
      else if (cart_cons_mode.eq.6) then
        write(ilog,*) 'Constraining all bond lengths as defined by user selection and &
 &constraining eligible water models (TipXp, SPC) to be rigid.'
        call strlims(cart_cons_file,t1,t2)
        write(ilog,*) 'Constraint set file : ',cart_cons_file(t1:t2)
      end if
    end if
    if (cart_cons_grps.gt.0) then
      if (shake_cnt.ne.cart_cons_grps) then
        write(ilog,*) 'Appropriate constraint groups are constrained by analytical SETTLE algorithm.'
      end if
      if (shake_cnt.gt.0) then
        if (cart_cons_method.eq.1) then
          write(ilog,*) 'Holonomic constraints are enforced via standard SHAKE algorithm.'
        else if (cart_cons_method.eq.2) then
          write(ilog,*) 'Holonomic constraints (distance only) are enforced via mixed SHAKE/P-SHAKE algorithm.'
        else if (cart_cons_method.eq.3) then
          write(ilog,*) 'Holonomic constraints (distance only) are enforced via P-SHAKE algorithm.'
        else if (cart_cons_method.eq.4) then
          write(ilog,*) 'Holonomic constraints (distance only) are enforced via LINCS algorithm.'
        end if
        if ((cart_cons_method.ge.1).AND.(cart_cons_method.le.3)) then
          write(ilog,'(1x,a,g12.4)') '10^4 x relative toler.: ',shake_tol*10000
          write(ilog,'(1x,a,g12.4)') '10^4 x absolute toler.: ',shake_atol*10000
          write(ilog,*) 'Maximum iterations  : ',shake_maxiter
        else if (cart_cons_method.eq.4) then
          write(ilog,'(1x,a,g12.4)') '10^4 x relative toler.: ',shake_tol*10000
          write(ilog,*) 'Start. expansion order: ',lincs_order
          write(ilog,*) 'Corr. iteration count : ',lincs_iter
        end if
      end if
      if (cart_cons_source.eq.1) then
        write(ilog,*) 'Values for constraints are taken from database information.'
      else if (cart_cons_source.eq.2) then
        write(ilog,*) 'Values for constraints are taken from force field in use (if possible).'
      else if (cart_cons_source.eq.3) then
        write(ilog,*) 'Values for constraints are taken from input structure (if possible).'
      end if
    end if
  end if
  if (((dyn_mode.ge.2).AND.(dyn_mode.le.5)).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    if (ens%sysfrz.eq.1) then
      write(ilog,*) 'No drift-removal for the system center-of-mass requested.'
    else if (ens%sysfrz.eq.2) then
      write(ilog,*) 'Removal of only the translational drift of the system center-of-mass requested.'
    else if (ens%sysfrz.eq.3) then
      write(ilog,*) 'Removal of both translational and rotational drift of the system center-of-mass requested.'
    end if
  end if
  if (dyn_mode.eq.5) then
    if (do_restart.EQV..true.) then
      if (in_dyncyc.EQV..true.) then
        write(ilog,*) 'Restarting in MD segment (ends at step ',curcyc_end,').'
      else
        write(ilog,*) 'Restarting in MC segment (ends at step ',curcyc_end,').'
      end if
    else
      write(ilog,21) ' Cycle len. for hyb. MD/MC init. MC: ',first_mccyclen
    end if
    write(ilog,21) ' Min. MC cyc. length for hyb. MD/MC: ',min_mccyclen
    write(ilog,21) ' Max. MC cyc. length for hyb. MD/MC: ',max_mccyclen
    write(ilog,21) ' Min. dyn. cyc. len. for hyb. MD/MC: ',min_dyncyclen
    write(ilog,21) ' Max. dyn. cyc. len. for hyb. MD/MC: ',max_dyncyclen
  else if (dyn_mode.eq.7) then
    if (do_restart.EQV..true.) then
      if (in_dyncyc.EQV..true.) then
        write(ilog,*) 'Restarting in LD segment (ends at step ',curcyc_end,').'
      else
        write(ilog,*) 'Restarting in MC segment (ends at step ',curcyc_end,').'
      end if
    else
      write(ilog,21) ' Cycle len. for hyb. LD/MC init. MC: ',first_mccyclen
    end if
    write(ilog,21) ' Min. MC cyc. length for hyb. LD/MC: ',min_mccyclen
    write(ilog,21) ' Max. MC cyc. length for hyb. LD/MC: ',max_mccyclen
    write(ilog,21) ' Min. dyn. cyc. len. for hyb. LD/MC: ',min_dyncyclen
    write(ilog,21) ' Max. dyn. cyc. len. for hyb. LD/MC: ',max_dyncyclen
  else if (dyn_mode.eq.8) then
    if (do_restart.EQV..true.) then
      if (in_dyncyc.EQV..true.) then
        write(ilog,*) 'Restarting in BD segment (ends at step ',curcyc_end,').'
      else
        write(ilog,*) 'Restarting in MC segment (ends at step ',curcyc_end,').'
      end if
    else
      write(ilog,21) ' Cycle len. for hyb. BD/MC init. MC: ',first_mccyclen
    end if
    write(ilog,21) ' Min. MC cyc. length for hyb. BD/MC: ',min_mccyclen
    write(ilog,21) ' Max. MC cyc. length for hyb. BD/MC: ',max_mccyclen
    write(ilog,21) ' Min. dyn. cyc. len. for hyb. BD/MC: ',min_dyncyclen
    write(ilog,21) ' Max. dyn. cyc. len. for hyb. BD/MC: ',max_dyncyclen
  end if
  else ! if not pdb_analyze 
  if ((ens%flag.eq.2).OR.(ens%flag.eq.1)) then
    write(ilog,*) 'Analysis assumes isochoric conditions (NVT/E).'
  else if ((ens%flag.eq.2).OR.(ens%flag.eq.1)) then
    write(ilog,*) 'Analysis assumes isobaric conditions (NPT/E).'
  else if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    write(ilog,*) 'Analysis assumes grand ensemble with isochoric conditions (muVT/E).'
  end if ! whether pdb_analyze
  end if
  if (dyn_mode.eq.6) then
    if (mini_mode.eq.1) then
      write(ilog,*) 'Using steepest-descent method as minimizer.'
    else if (mini_mode.eq.2) then
      write(ilog,*) 'Using conjugate-gradient method as minimizer.'
    else if (mini_mode.eq.3) then
      write(ilog,*) 'Using quasi-Newton (BFGS) method as minimizer.'
    else if (mini_mode.eq.4) then
      write(ilog,*) 'Using stochastic minimizer.'
    end if
    write(ilog,21) ' Maximum number of iterations      : ',nsim
    write(ilog,21) ' First steps discarded for analysis: ',nequil
    if (mini_mode.ne.4) then
      write(ilog,20) ' Global Iteration step size        : ',mini_stepsize
      if (fycxyz.eq.1) then
        write(ilog,20) ' Translational base step size      : ',mini_xyzstep
        write(ilog,20) ' Rotational (RB) base step size    : ',mini_rotstep
        write(ilog,20) ' Torsional base step size          : ',mini_intstep
      end if
    end if
    if (mini_mode.eq.3) then
      write(ilog,21) ' Memory span for BFGS Hessian upd. : ',mini_mem
      write(ilog,20) ' Uphill energy tolerance criterion : ',mini_uphill
    end if
    write(ilog,20) ' Convergence crit. (RMS gradient)  : ',mini_econv    
    if (mini_mode.eq.4) then
      if ((fycxyz.eq.1).AND.(mini_sc_sdsteps.gt.0)) then
        write(ilog,20) ' Translational base step (for GRMS): ',mini_xyzstep
        write(ilog,20) ' Rotational base step (for GRMS)   : ',mini_rotstep
        write(ilog,20) ' Torsional base step (for GRMS)    : ',mini_intstep
      end if
      write(ilog,21) ' Steepest-descent steps initially  : ',mini_sc_sdsteps
      write(ilog,20) ' Starting temperature              : ',mini_sc_tbath
      write(ilog,21) ' Steps before starting cooling     : ',floor(mini_sc_heat*nsim)
    end if
  end if
  if ((mc_acc_crit.eq.3).AND.(((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.&
 &(dyn_mode.eq.8)).AND.(pdb_analyze.EQV..false.))) then
    write(ilog,*) 'This simulation uses Wang-Landau (WL) sampling for MC acceptance.'
    write(ilog,*) 'Since the goal of such a calculation is conceptually different, the produced&
 & trajectory and analyses relying on it DO NOT CORRESPOND to well-defined averages or &
 &distributions. WEIGHTED REANALYSIS IS NECESSARY TO RECOVER BOLTZMANN-DISTRIBUTED &
 &ENSEMBLES.'
    if (dyn_mode.ne.1) then
      write(ilog,*) 'In case of a hybrid Monte Carlo / dynamics run, only the dynamics &
 &portion samples from the specified ensemble, and does not contribute per se to the &
 &WL procedure.'
    end if
    if (wld%wl_mode.eq.1) then
      write(ilog,*) 'WL sampling occurs in energy space (yields g(E)).'
    else if (wld%wl_mode.eq.2) then
      write(ilog,*) 'WL sampling occurs in the reaction coordinate space(s) of &
 &selected molecule(s).'
    else if ((wld%wl_mode.eq.3).AND.(wld%dimensionality.eq.1)) then
      write(ilog,*) 'WL sampling occurs in Boltzmann-weighted energy space.'
    else if ((wld%wl_mode.eq.3).AND.(wld%dimensionality.eq.2)) then
      write(ilog,*) 'WL sampling occurs in Boltzmann-weighted energy space and the &
 &reaction coordinate space of a selected molecule.'
    end if
    write(ilog,21) ' Dimensionality                    : ',wld%dimensionality
    write(ilog,21) ' Visitation criterion for f-updates: ',wld%hvmode
    write(ilog,21) ' Frequency for updating histograms : ',wld%hufreq
    write(ilog,21) ' Frequency for f-update checks     : ',wld%gh_flatcheck_freq
    write(ilog,21) ' Buffer steps for f-updates        : ',wld%buffer
    if (wld%freeze.EQV..true.) then
      write(ilog,*) 'Convergence analysis restricted to bins encountered during buffer steps.'
    else
      write(ilog,*) 'Convergence analysis can extend beyond bins encountered during buffer steps.'
    end if
    if (do_restart.EQV..false.) then
      write(ilog,20) ' Starting value for f-factor (log) : ',wld%fval
      if (wld%dimensionality.eq.1) then
        write(ilog,20) ' Bin size (units -> doc.)          : ',wld%g_binsz
        write(ilog,20) ' Center of largest (units -> doc.) : ',wld%g_max
      else if (wld%dimensionality.eq.2) then
        write(ilog,22) ' Bin sizes (units -> doc.)         : ',wld%g_binsz2d(:)
        write(ilog,22) ' Centers of largest (units -> doc.): ',wld%g_max2d(:)
      end if
      if (wld%exrule.eq.1) then
        write(ilog,*) 'WL histograms are fixed in size.'
      else if (wld%exrule.eq.2) then
        write(ilog,*) 'WL histograms will grow dynamically toward lower values.'
      else if (wld%exrule.eq.3) then
        write(ilog,*) 'WL histograms will grow dynamically in both directions.'
      end if
      if (finit_wanglandau.EQV..true.) then
        call strlims(wld%ginitfile,t1,t2)
        write(ilog,*) 'File for providing init. histogram: ',wld%ginitfile(t1:t2)
      end if
    end if
    if (wld%wl_mode.eq.2) then
      write(ilog,21) ' Mol. 1 for reaction-coordinate WL : ',wld%wl_mol(1)
      write(ilog,21) ' RC 1 for reaction-coordinate   WL : ',wld%wl_rc(1)
    end if
    if (wld%dimensionality.eq.2) then
      write(ilog,21) ' Mol. 2 for two-dimensional WL     : ',wld%wl_mol(2)
      write(ilog,21) ' RC 2 for two-dimensional WL       : ',wld%wl_rc(2)
    end if
    if (debug_wanglandau.EQV..true.) then
      write(ilog,*) 'Temporary histograms and debugging information for WL statistics will be written.'
    end if
  else if (((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.&
 &(dyn_mode.eq.8)).AND.(pdb_analyze.EQV..false.)) then
    if (mc_acc_crit.eq.1) then
      write(ilog,*) 'Using Metropolis-Hastings acceptance criterion for Monte Carlo moves.'
    else if (mc_acc_crit.eq.2) then
      write(ilog,*) 'Using Fermi acceptance probability for Monte Carlo moves.'
    end if
  end if
  if (((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)).AND.(pdb_analyze.EQV..false.)) then
    if (do_n2loop.EQV..true.) then
      write(ilog,*) 'Monte Carlo steps use N^2 sum (ignoring FMCSC_CUTOFFMODE) as reference energy &
 &(FMCSC_CHECKFREQ+FMCSC_N2LOOP).'
    else
      write(ilog,*) 'Monte Carlo steps use truncation scheme (FMCSC_CUTOFFMODE) as reference energy&
 & (FMCSC_CHECKFREQ+FMCSC_N2LOOP).'
    end if
  end if
  write(ilog,*)
  write(ilog,*) '**ENERGY TERMS:'
  if (use_hardsphere.EQV..true.) then
    if (scale_IPP.gt.0.0) then
      write(ilog,20) ' Hard sphere potential used instead of IPP.'
      write(ilog,20) ' Penalty for nuclear overlap       : ',screenbarrier
    else
      write(ilog,20) ' Scale for inverse power potential : ',scale_IPP
    end if
  else
    write(ilog,20) ' Scale for inverse power potential : ',scale_IPP
    if (scale_IPP.gt.0.0) then
      write(ilog,21) ' Exponent for IPP                  : ',nindx
    end if
  end if
  if ((scale_IPP.gt.0.0).OR.(scale_attLJ.gt.0.0)) then
    if (reduce.lt.0.0) then
      write(ilog,20) ' Scale for non-14 sigma parameters : ',abs(reduce)
    else
      write(ilog,20) ' Scale for all sigma parameters    : ',reduce
    end if
    write(ilog,21) ' Sigma combination rule code       : ',sigrule
    write(ilog,21) ' Epsilon combination rule code     : ',epsrule
  end if
  write(ilog,20) ' Scale for attractive Lennard-Jones: ',scale_attLJ
  write(ilog,20) ' Scale for WCA potential           : ',scale_WCA
  if (use_WCA.EQV..true.) then
    write(ilog,20) ' Attract. factor for WCA pot.      : ',par_WCA(2)
    write(ilog,20) ' Cutoff for WCA potential          : ',par_WCA(1)
    write(ilog,20) ' WCA attraction coefficent 1       : ',par_WCA(3)
    write(ilog,20) ' WCA attraction coefficent 2       : ',par_WCA(4)
  end if
  if ((use_attLJ.EQV..true.).OR.(use_IPP.EQV..true.).OR.(use_WCA.EQV..true.)) then
    write(ilog,20) ' Fudge factor for steric 1-4 terms : ',fudge_st_14
  end if
  write(ilog,20) ' Scale for struct. correction terms: ',scale_CORR
  write(ilog,20) ' Scale for bond length potentials  : ',scale_BOND(1)
  write(ilog,20) ' Scale for bond angle potentials   : ',scale_BOND(2)
  if ((use_BOND(1).EQV..true.).OR.(use_BOND(2).EQV..true.).OR.(use_BOND(3).EQV..true.).OR.(use_BOND(4).EQV..true.)) then
    if (guess_bonded.eq.1) then
      write(ilog,*) 'Missing bonded interaction parameters (bonds and angles)&
 & will be constructed from crystallographic or input geometries.'
    else if (guess_bonded.eq.2) then
      write(ilog,*) 'Missing bonded interaction parameters (bonds, angles, and impropers)&
 & will be constructed from crystallographic or input geometries.'
    else if (guess_bonded.eq.3) then
     write(ilog,*) 'Missing bonded interaction parameters (bonds, angles, impropers, and)&
 &torsional potentials describing electronic barriers) will be constructed from crystallographic or input geometries and &
 &atom type and connectivity information (for torsional potentials).'
    else
      write(ilog,*) 'Missing bonded interaction parameters (bonds, angles, etc.)&
 & will not be constructed. This may prevent simulations in Cartesian space.'
    end if
  end if
  write(ilog,20) ' Scale for improper dihedral terms : ',scale_BOND(3)
  if ((scale_BOND(3).gt.0.0).AND.(improper_conv(1).ne.1)) then
    write(ilog,*) 'Using alternate (AMBER/OPLS) convention for improper dihedral potentials.'
  end if
  write(ilog,20) ' Scale for torsional potentials    : ',scale_BOND(4)
  write(ilog,20) ' Scale for CMAP potentials         : ',scale_BOND(5)
  if (use_BOND(5).EQV..true.) then
    call strlims(cm_dir,t1,t2)
    write(ilog,*) 'Directory for CMAP grid corr.s    : ',cm_dir(t1:t2)
    write(ilog,21) ' B-Spline order (if applicable)    : ',cm_splor
  end if
  write(ilog,20) ' Scale for implicit solvent poten. : ',scale_IMPSOLV
  if (use_IMPSOLV.EQV..true.) then
    write(ilog,20) ' Thickness of solvation layer in A : ',par_IMPSOLV(1)
    write(ilog,21) ' Mode for sphere overlap calculat. : ',par_IMPSOLV2(1)
    write(ilog,20) ' Dielectric of implicit solvent    : ',par_IMPSOLV(2)
    if (par_IMPSOLV2(2).eq.1) then
      write(ilog,*) 'Using smooth sigmoidal function for SAV-to-solvation state mapping (F.o.S.).'
    else if (par_IMPSOLV2(2).eq.2) then
      write(ilog,*) 'Using stair-stepped sigmoidal function for SAV-to-solvation state mapping (F.o.S.).'
    end if
    if ((par_IMPSOLV2(2).eq.1).OR.(par_IMPSOLV2(2).eq.2)) then
      write(ilog,20) ' Sigmoidal decay constant (F.o.S.) : ',par_IMPSOLV(3)
      write(ilog,20) ' Rel. center of sigm. fxn (F.o.S.) : ',par_IMPSOLV(6)
    end if
    if (par_IMPSOLV2(2).eq.2) then
      write(ilog,20) ' Stair step granularity (F.o.S.)   : ',par_IMPSOLV(12)
      write(ilog,20) ' Stair step steepness (F.o.S.)     : ',par_IMPSOLV(13)
      write(ilog,20) ' Stair step position (F.o.S.)      : ',par_IMPSOLV(14)
    end if
    if (use_Tdepfosvals.EQV..true.) then
      write(ilog,*) 'Using T-dependent reference values for model compound free energies of solvation.'
    else
      write(ilog,*) 'Using fixed reference values for model compound free energies of solvation.'
    end if
    write(ilog,20) ' Assumed ref. temp. in K for F.o.S.: ',fos_Tdepref
    if (scrq_model.eq.1) then
      write(ilog,*) 'Charged-weighted, group-consistent dielectric screening is used.'
    else if (scrq_model.eq.2) then
      write(ilog,*) 'Purely atom-based dielectric screening is used.'
    else if (scrq_model.eq.3) then
      write(ilog,*) 'Mixed group-based environmental and distance-dependent screening is used (base model 1).'
    else if (scrq_model.eq.4) then
      write(ilog,*) 'Pure distance-dependent dielectric is used.'
    else if (scrq_model.eq.6) then
      write(ilog,*) 'Generalized atom-based dielectric screening is used.'
    else if (scrq_model.eq.5) then
      write(ilog,*) 'Generalized charge-weighted, group-consistent dielectric screening is used.'
    else if (scrq_model.eq.9) then
      write(ilog,*) 'Mixed atom-based environmental and distance-dependent screening is used (base model 2).'
    else if (scrq_model.eq.7) then
      write(ilog,*) 'Mixed generalized group-based environmental and distance-dependent screening is used (base model 5).'
    else if (scrq_model.eq.8) then
      write(ilog,*) 'Mixed generalized atom-based environmental and distance-dependent screening is used (base model 6).'
    end if
    if (scrq_model.ne.4) then
      if (par_IMPSOLV2(3).eq.1) then
        write(ilog,*) 'Using smooth sigmoidal function for SAV-to-solvation state mapping (q-Scr.).'
      else if (par_IMPSOLV2(3).eq.2) then
        write(ilog,*) 'Using stair-stepped sigmoidal function for SAV-to-solvation state mapping (q-Scr.).'
      end if
      if ((par_IMPSOLV2(3).eq.1).OR.(par_IMPSOLV2(3).eq.2)) then
        write(ilog,20) ' Sigmoidal decay constant (q-Scr.) : ',par_IMPSOLV(4)
        write(ilog,20) ' Rel. center of sigm. fxn (q-Scr.) : ',par_IMPSOLV(7)
      end if
      if (par_IMPSOLV2(3).eq.2) then
        write(ilog,20) ' Stair step granularity (q-Scr.)   : ',par_IMPSOLV(15)
        write(ilog,20) ' Stair step steepness (q-Scr.)     : ',par_IMPSOLV(16)
        write(ilog,20) ' Stair step position (q-Scr.)      : ',par_IMPSOLV(17)
      end if
    end if
    if ((scrq_model.eq.3).OR.(scrq_model.eq.4).OR.(scrq_model.gt.6)) then
      write(ilog,20) ' Value for contact dielectric      : ',1./par_IMPSOLV(8)
    end if
    if ((scrq_model.eq.3).OR.(scrq_model.gt.6)) then
      write(ilog,20) ' Mixing contr. for distance-dep.   : ',par_IMPSOLV(9)
    end if
    if ((scrq_model.ge.5).AND.(scrq_model.le.8)) then
      write(ilog,21) ' Order of generalized mean (q-scr) : ',i_sqm
    end if
!    if (par_IMPSOLV(10).ne.0.0) then
!      write(ilog,20) ' SAV-response compressibility      : ',&
! &par_IMPSOLV(10)
!    end if
!    if (par_IMPSOLV(11).ne.0.0) then
!      write(ilog,20) ' SEV-response compressibility      : ',&
! &par_IMPSOLV(11)
!    end if
  end if
  write(ilog,20) ' Scale for polar interaction term  : ',scale_POLAR
  if (use_POLAR.EQV..true.) then
    write(ilog,20) ' Fudge factor for polar 1-4 terms  : ',fudge_el_14
    write(ilog,20) ' Polarization shift for 1-amides   : ',primamide_cis_chgshft
    write(ilog,20) ' Fractional charge tolerance sett. : ',dpgrp_neut_tol
  end if
  write(ilog,20) ' Scale for torsional terms (bias)  : ',scale_TOR
  write(ilog,20) ' Scale for polymeric terms (bias)  : ',scale_POLY
  write(ilog,20) ' Scale for global 2-structure term : ',scale_ZSEC
  if (use_ZSEC.EQV..true.) then
    write(ilog,20) ' Target alpha-content              : ',par_ZSEC(1)
    write(ilog,20) ' Harmonic alpha-bias spring const. : ',par_ZSEC(2)
    write(ilog,22) ' Center of alpha-basin             : ',par_ZSEC2(1),par_ZSEC2(2)
    write(ilog,20) ' Radius of alpha-basin             : ',par_ZSEC2(3)
    write(ilog,20) ' Steepness of alpha-basin          : ',par_ZSEC2(4)
    write(ilog,20) ' Target beta-content               : ',par_ZSEC(3)
    write(ilog,20) ' Harmonic beta-bias spring const.  : ',par_ZSEC(4)
    write(ilog,22) ' Center of beta-basin              : ',par_ZSEC2(5),par_ZSEC2(6)
    write(ilog,20) ' Radius of beta-basin              : ',par_ZSEC2(7)
    write(ilog,20) ' Steepness of beta-basin           : ',par_ZSEC2(8)
  end if
  write(ilog,20) ' Scale for DSSP restraint term     : ',scale_DSSP
  if (use_DSSP.EQV..true.) then
    write(ilog,20) ' Target H-score                    : ',par_DSSP(7)
    write(ilog,20) ' Harmonic H-bias spring constant   : ',par_DSSP(8)
    write(ilog,20) ' Target E-score                    : ',par_DSSP(9)
    write(ilog,20) ' Harmonic E-bias spring constant   : ',par_DSSP(10)
  end if
  write(ilog,20) ' Scale for tabulated pot. term     : ',scale_TABUL
  if (use_TABUL.EQV..true.) then
    write(ilog,21) ' Highest potential no. referred to : ',tbp%nos
    write(ilog,21) ' Number of distance bins           : ',tbp%bins
    if (tbp%tangmode.eq.0) then
      write(ilog,*) 'Tangents are all constructed numerically (Kochanek-Bartels spline).'
    else if (tbp%tangmode.eq.1) then
      write(ilog,*) 'Tangents are all read from file (see above).'
    else if (tbp%tangmode.eq.2) then
      write(ilog,*) 'Tangents are partially constructed numerically (Kochanek-Bartels spline) and &
 &partially read from file.'
    end if
    if (tbp%tangmode.ne.1) then
      write(ilog,20) ' Tightness for numerical tangents  : ',tbp%ipprm(1)
      write(ilog,20) ' L/R bias for numerical tangents   : ',tbp%ipprm(2)
    end if
  end if
  write(ilog,20) ' Scale for distance restraint term : ',scale_DREST
  write(ilog,20) ' Scale for compartmentaliz. potent.: ',scale_OSMO
  if (use_OSMO.EQV..true.) then
    write(ilog,21) ' Number of symmetric compartments  : ',2**par_OSMO(1)
  end if
  write(ilog,20) ' Scale for EM restraint term       : ',scale_EMICRO
  if (use_EMICRO.EQV..true.) then
    write(ilog,21) ' Mode (instantaneous vs. averaged) : ',empotmode
    if (empotmode.eq.2) then
      write(ilog,28) ' Inst. weight in avg. :',emiw
    end if
    if (empotprop.eq.1) then
      write(ilog,*) 'Using mass as atomic property for density restraint.'
    else if (empotprop.eq.2) then
      write(ilog,*) 'Using atomic number as property for density restraint.'
    else if (empotprop.eq.3) then
      write(ilog,*) 'Using charge as atomic property for density restraint.'
    end if
    call strlims(emmapfile,t1,t2)
    if (allocated(emcoarsegrid%rslc).EQV..true.) then
      write(ilog,*) 'Using slice heuristic for energy and force evaluations.'
    else if (allocated(emcoarsegrid%rlin).EQV..true.) then
      write(ilog,*) 'Using line heuristic for energy and force evaluations.'
    else if (allocated(emcoarsegrid%rblk).EQV..true.) then
      write(ilog,*) 'Using block heuristic for energy and force evaluations.'
    else
      write(ilog,*) 'No heuristic for energy and force evaluations in use.'
    end if
    write(ilog,*) 'Reference map (filename)          : ',emmapfile(t1:t2)
    write(ilog,28) ' Threshold level for (bio)molecules: ',emthreshdensity
    write(ilog,28) ' Assumed masss of (bio)molecules   : ',emtotalmass
    write(ilog,28) ' Computed average mol. density     : ',emavgdensity
    write(ilog,28) ' Input background dens. (set/hist.): ',emoptbg
    if (emtruncate.gt.-1.0*(HUGE(emtruncate)/1000.0)) then
      write(ilog,28) ' Input dens. may be truncated below: ',emtruncate
    end if
    if (emflatval.lt.HUGE(emflatval)) then
      write(ilog,28) ' Input density may be flattened at : ',emflatval
    end if
    if ((fdlts(1).gt.0.0).AND.(fdlts(2).gt.0.0).AND.(fdlts(3).gt.0.0)) then
      write(ilog,*) 'Input map was adjusted in resolution before use.'
      write(ilog,*) 'Values given below are for the adjusted map.'
    end if
    write(ilog,9) ' Map increments (A)   : ',emingrid%deltas(1:3)
    write(ilog,9) ' Map origin (A)       : ',emingrid%origin(1:3)
    write(ilog,8) ' Map grid dimensions  : ',emingrid%dim(1:3)
    if (emcalc.gt.nsim) then
      if (empotprop.eq.1) then
        write(ilog,*) 'Computing spatial density based on atomic masses.'
      else if (empotprop.eq.2) then
        write(ilog,*) 'Computing spatial density based on atomic number.'
      else if (empotprop.eq.3) then
        write(ilog,*) 'Computing spatial density based on atomic charge.'
      else if (empotprop.eq.4) then
        call strlims(empropfile,t1,t2)
        write(ilog,*) 'Computing spatial density based on user-defined property.'
        write(ilog,*) 'Input file with propert.: ',empropfile(t1:t2)
      end if
      write(ilog,55) ' Background density      : ',embgdensity
      write(ilog,21) ' B-spline order          : ',emsplor
      if (.NOT.((bnd_shape.eq.1).AND.(bnd_type.eq.1))) write(ilog,55) ' Buffer scale factor     : ',embuf
      write(ilog,9) ' Spatial resol.  (x,y,z) : ',emgrid%deltas(1:3)
      write(ilog,9) ' Map origin      (x,y,z) : ',emgrid%origin(1:3)
      write(ilog,8) ' Grid dimensions (x,y,z) : ',emgrid%dim(1:3)
    else
      write(ilog,*) 'See below for further settings.'
    end if
  end if
!  write(ilog,20) ' Scale for LCT pot. term           : ',scale_LCTOR
!  if (use_LCTOR.EQV..true.) then
!    write(ilog,21) ' Number of LCs to use              : ',&
! &ntorlcs
!    write(ilog,20) ' LCT distribution bin size         : ',&
! &par_LCTOR(1)
!    write(ilog,21) ' Number of distribution bins       : ',&
! &par_LCTOR2(1)
!  end if
  if (use_FEG.EQV..true.) then
    dummy = 0
    do i=1,nseq
      if (par_FEG(i).EQV..true.) then
        dummy = dummy + 1
      end if
    end do
    write(ilog,21) ' Number of ghost residues          : ',dummy
    if (fegmode.eq.1) then
      write(ilog,*) 'Ghost-ghost interactions use the default (background) Hamiltonian.'
    else
      write(ilog,*) 'Ghost-ghost interactions use the scaled (ghosted) Hamiltonian.'
    end if
    write(ilog,20) ' Scale for ghost IPP               : ',scale_FEGS(1)
    write(ilog,20) ' Scale for ghost attractive LJ     : ',scale_FEGS(3)
    if ((scale_FEGS(1).gt.0.0).OR.(scale_FEGS(3).gt.0.0)) then
      if (fegljmode.eq.1) then
        write(ilog,*) 'FEG-LJ-scaling is linear.' 
      else if (fegljmode.eq.2) then
        write(ilog,*) 'Ghost LJ-term is a polynomially scaled soft-core potential.'
        write(ilog,20) ' Soft-core radius for ghost-LJ     : ',par_FEG2(3)
        write(ilog,20) ' Polynomial exponent for ghost-LJ  : ',par_FEG2(8)
        write(ilog,20) ' Soft-core exponent for ghost-LJ   : ',par_FEG2(7)
      else 
        write(ilog,*) 'Ghost LJ-term is an exponentially scaled soft-core potential.'
        write(ilog,20) ' Soft-core radius for ghost-LJ     : ',par_FEG2(3)
        write(ilog,20) ' Exponential factor for ghost-LJ   : ',par_FEG2(8)
        write(ilog,20) ' Soft-core exponent for ghost-LJ   : ',par_FEG2(7)
      end if
    end if
    write(ilog,20) ' Scale for ghost polar term        : ',scale_FEGS(6)
    if (scale_FEGS(6).gt.0.0) then
      if (fegcbmode.eq.1) then
        write(ilog,*) 'Ghost polar scaling is linear.' 
      else if (fegcbmode.eq.2) then
        write(ilog,*) 'Ghost polar term is a polynomially scaled soft-core potential.'
        write(ilog,20) ' Soft-core radius for ghost-POLAR  : ',par_FEG2(4)
        write(ilog,20) ' Polynomial exp. for ghost-POLAR   : ',par_FEG2(12)
        write(ilog,20) ' Soft-core exponent for ghost-POLAR: ',par_FEG2(11)
      end if
    end if
    if (fegmode.eq.2) then
      write(ilog,20) ' Scale for ghost bond length term  : ',scale_FEGS(15)
      write(ilog,20) ' Scale for ghost bond angle term   : ',scale_FEGS(16)
      write(ilog,20) ' Scale for ghost improper dih. term: ',scale_FEGS(17)
      write(ilog,20) ' Scale for ghost torsional term    : ',scale_FEGS(18)
    end if
  end if
  if (do_accelsim.EQV..true.) then
    write(ilog,*) 'Some energy terms may be sculpted (see summary above).'
    write(ilog,*) 'Frequency of instantaneous output of simulation step '
    write(ilog,*) 'numbers and corresponding weights : ',hmjam%prtfrmwts
    write(ilog,71) ' Threshold weight for print-out    : ',hmjam%threshwt
  else
    write(ilog,*) 'Using unperturbed energy landscape (sculpting is off -> FMCSC_SCULPT).'
  end if
!
  if (((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8))&
 &.AND.(pdb_analyze.EQV..false.)) then
  write(ilog,*)
  write(ilog,*) '**TECHNICAL SETTINGS:'
  write(ilog,25) ' Frequency of (S)GC moves      : ',particleflucfreq
  write(ilog,25) ' Frequency of rigid-body moves : ',rigidfreq
  if (rigidfreq.gt.0.0) then
    write(ilog,25) ' Freq. of single mol. RB moves : ',1.0-clurb_freq
    if ((1.0-clurb_freq).gt.0.0) then
      if (use_coupledrigid.EQV..true.) then
        write(ilog,*) 'RB rotation and translation are coupled.'
      else
        write(ilog,*) 'RB rotation and translation are decoupled.'
        write(ilog,25) ' Freq. of rot.-moves (in RB)   : ',rotfreq
      end if
    end if
    write(ilog,25) ' Freq. of cluster RB moves     : ',clurb_freq
    if (clurb_freq.gt.0.0) then
      write(ilog,25) ' Freq. of random cluster moves : ',clurb_rdfreq
      if (clurb_rdfreq.gt.0.0) then
        write(ilog,21) ' Maximum cluster size          : ',clurb_maxsz
      end if
      write(ilog,25) ' Freq. of struc. cluster moves : ',1.0-clurb_rdfreq
      if (clurb_rdfreq.lt.1.0) then
        write(ilog,25) ' Distance criterion (in A)     : ',clurb_discrit
        write(ilog,21) ' Cluster reset interval        : ',clurb_reset
        write(ilog,25) ' Bail-out frequency            : ',clurb_strbail
      end if
    end if
    write(ilog,*) 'Fraction of full randomization'
    write(ilog,25) ' attempts amongst all RB moves : ',rigid_randfreq
    if (rigid_randfreq.gt.0.0) then
      write(ilog,25) ' Max. box size ratio for pos.  : ',trans_xyzbuf
    end if
    if (rigid_randfreq.lt.1.0) then
     write(ilog,*) 'Fraction of local moves amongst'
     write(ilog,25) ' all RB moves                  : ',1.0-rigid_randfreq
     write(ilog,25) ' Translational step size (in A): ',trans_stepsz
     write(ilog,25) ' Rotational step size (in deg.): ',rot_stepsz 
    end if
  end if
  if (rigidfreq.lt.1.0) then
  write(ilog,*) 'Fraction of torsional moves'
  write(ilog,25) ' amongst total moves           : ',1.0-rigidfreq
  if (use_lctmoves.EQV..true.) then
!    write(ilog,*) 'All torsional sampling moves are LCT moves.'
!    write(ilog,25) ' Wall height for range determ. : ',par_LCTOR(3)
  else
    if (use_globmoves.EQV..true.) then
      write(ilog,25) ' Frequency of pure s.c. moves  : ',chifreq
    else
      write(ilog,25) ' Frequency of sidechain moves  : ',chifreq
    end if
    if (chifreq.gt.0.0) then
      if (phfreq.gt.0.0) then
        write(ilog,*) 'Frequency of titration moves '
        write(ilog,25) ' amongst chi-cycles            : ',phfreq
        write(ilog,21) ' Number of eligible residues   : ',nionres
        write(ilog,25) ' Continuum pH                  : ',phs
        write(ilog,25) ' Continuum ionic strength      : ',ionicstr
      end if
      if (phfreq.lt.1.0) then
        write(ilog,*) 'Frequency of dihedral moves '
        write(ilog,25)' amongst sidechain moves       : ',1.0-phfreq
        write(ilog,*) 'Number of eligible residues   : ',chilst%nr
        write(ilog,*) 'Number of chi-angles in move  : ',nrchis
        write(ilog,*) 'Fraction of full randomization'
        write(ilog,25) ' attempts amongst chi moves    : ',chi_randfreq 
        if (chi_randfreq.lt.1.0) then
          write(ilog,25) ' Step size for local chi moves : ',chi_stepsz
        end if
      end if
    end if
    if (chifreq.lt.1.0) then
      write(ilog,*) 'Chain alignment mode          : ',align_NC
      write(ilog,*) 'Fraction of concerted rotation moves amongst'
      write(ilog,25) ' remaining dihedral moves      : ',crfreq
      if (crfreq.gt.0.0) then
        if (angcrfreq.lt.1.0) then
          write(ilog,*) 'Fraction of purely torsional amongst '
          write(ilog,25) ' all concerted rotation moves  : ',1.0-angcrfreq
          if (torcrfreq.lt.1.0) then
            write(ilog,*) 'Fraction of Sjunnesson moves amongst'
            write(ilog,25) ' torsional CR moves            : ',1.0-torcrfreq
            write(ilog,*) 'Implementation mode           : ',cr_mode
            write(ilog,25) ' Distribution width (rad)      : ',1.0/cr_a
            write(ilog,*) 'Number of torsions used       : ',nr_crdof
            write(ilog,*) 'Number of biasing restraints  : ',nr_crfit
            write(ilog,25) ' Bias scaling factor           : ',cr_b
          end if
          if (torcrfreq.gt.0.0) then
            write(ilog,*) 'Fraction of Dinner-Ulmschneider moves'
            write(ilog,25) ' amongst torsional CR moves    : ',torcrfreq
            write(ilog,25) ' Fraction to use omega bonds   : ',torcrfreq_omega
            write(ilog,25) ' Pre-rotation scaling for omega: ',1.0/UJ_params(6)
            if (torcrfreq_omega.gt.0.0) then
              write(ilog,*) ' Minimum pre-rot. d.o.f.s     : ',torcrminsz
              write(ilog,*) ' Maximum pre-rot. d.o.f.s     : ',torcrmaxsz
            end if
            write(ilog,25) ' Frac. not to use omega bonds  : ',1.0-torcrfreq_omega
            if (torcrfreq_omega.lt.1.0) then
              write(ilog,*) ' Minimum pre-rot. d.o.f.s     : ',torcrminsz2
              write(ilog,*) ' Maximum pre-rot. d.o.f.s     : ',torcrmaxsz2
            end if
          end if
        end if
        if (angcrfreq.gt.0.0) then
          write(ilog,*) 'Fraction of Ulmschneider moves amongst'
          write(ilog,25) ' all concerted rotation moves  : ',angcrfreq
          write(ilog,*) 'Minimum pre-rotation length   : ',ujminsz-2
          write(ilog,*) 'Maximum pre-rotation length   : ',ujmaxsz-2
          write(ilog,25) '  Interval for 1D-root search  : ',UJ_params(5)
          write(ilog,26) '  Angular scaling factor       : ',1.0/UJ_params(3)
        end if
        if ((angcrfreq.gt.0.0).OR.(torcrfreq.gt.0.0)) then
          write(ilog,*) 'Parameters for all exact CR methods as follows:'
          write(ilog,*) ' Mode setting (if applicable)   : ',torcrmode
          write(ilog,25) '  Distribution width in deg.   : ',RADIAN*sqrt(1.0/UJ_params(1))
          write(ilog,26) '  Bias strength factor         : ',UJ_params(2)
          write(ilog,25) '  Stepsize for 1D-root search  : ',UJ_params(4)
          if (torcrmode.eq.2) then
            write(ilog,*) ' Max. nr. of closure attempts : ',maxcrtries
          end if
        end if
      end if
      if (crfreq.lt.1.0) then
        write(ilog,*) 'Fraction of omega moves amongst'
        write(ilog,25) ' remaining dihedral moves      : ',omegafreq
        if (omegafreq.gt.0.0) then
          write(ilog,*) 'Fraction of full randomization'
          write(ilog,25) ' attempts amongst omega moves  : ',omega_randfreq
          if (omega_randfreq.lt.1.0) then
            write(ilog,25) ' Step sz. for local omega moves: ',omega_stepsz
          end if
        end if
        if (omegafreq.lt.1.0) then
          write(ilog,*) 'Fraction of nucleic acid moves amongst '
          write(ilog,25) ' remaining dihedral moves      : ',nucfreq
          if ((nucfreq.gt.0.0).AND.(nuccrfreq.gt.0.0)) then
            write(ilog,*) 'Fraction of Dinner-Ulmschneider nucleic acid'
            write(ilog,25) ' CR moves amongst n.a. moves   : ',nuccrfreq
            write(ilog,*) ' Minimum pre-rot. d.o.f.s     : ',nuccrminsz
            write(ilog,*) ' Maximum pre-rot. d.o.f.s     : ',nuccrmaxsz
            if ((angcrfreq.le.0.0).AND.((torcrfreq.le.0.0).OR.(crfreq.le.0.0))) then
              write(ilog,*) 'Parameters for all exact CR methods as follows:'
              write(ilog,*) ' Mode setting                 : ',torcrmode
              write(ilog,25) '  Distribution width in deg.   : ',RADIAN*sqrt(1.0/UJ_params(1))
              write(ilog,26) '  Bias strength factor         : ',UJ_params(2)
              write(ilog,25) '  Stepsize for 1D-root search  : ',UJ_params(4)
              if (torcrmode.eq.2) then
                write(ilog,*) ' Max. nr. of closure attempts : ',maxcrtries
              end if
            end if
          end if
          if ((nucfreq.gt.0.0).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.gt.0.0)) then
            write(ilog,*) 'Fraction of sugar puckering moves amongst'
            write(ilog,25) ' nucleic acid moves            : ',nucpuckfreq
            write(ilog,*) 'Fraction of reflection (pucker inversion) moves'
            write(ilog,25) ' amongst sugar puckering moves : ',nucpuckrdfreq
            if (nucpuckrdfreq.lt.1.0) then
              write(ilog,25) '  Stepsize for dihedral pert.  : ',pucker_distp
              write(ilog,25) '  Stepsize for angular pert.   : ',pucker_anstp 
            end if
          end if
          if ((nucfreq.gt.0.0).AND.(nuccrfreq.lt.1.0).AND.(nucpuckfreq.lt.1.0)) then
            write(ilog,*) 'Number of dihedrals in a move : ',nrnucim
            write(ilog,*) 'Fraction of full randomization'
            write(ilog,25) ' attempts amongst n.a. moves   : ',nuc_randfreq
            if (nuc_randfreq.lt.1.0) then
              write(ilog,25) ' Step sz. for local n.a. moves : ',nuc_stepsz
            end if
          end if
          if (nucfreq.lt.1.0) then
            write(ilog,*) 'Fraction of polypeptide puckering moves amongst'
            write(ilog,25) ' remaining dihedral moves      : ',puckerfreq
            if (puckerfreq.gt.0.0) then
              write(ilog,*) 'Fraction of reflection (pucker inversion) moves'
              write(ilog,25) ' amongst those puckering moves : ',puckerrdfreq
            end if
            if ((puckerfreq.gt.0.0).AND.(puckerrdfreq.lt.1.0)) then
              write(ilog,25) '  Stepsize for dihedral pert.  : ',pucker_distp
              write(ilog,25) '  Stepsize for angular pert.   : ',pucker_anstp 
            end if
            if (puckerfreq.lt.1.0) then
              write(ilog,*) 'Fraction of single dihedral (OTHER) moves amongst'
              write(ilog,25) ' remaining dihedral moves      : ',otherfreq
              if (otherfreq.gt.0.0) then
                write(ilog,*) 'Fraction of OTHER moves operating on residues not '
                write(ilog,25) ' natively supported by CAMPARI :',other_unkfreq
                if (other_unkfreq.lt.1.0) then
                  write(ilog,*) 'Fraction of all OTHER moves operating on d.o.f.s '
                  write(ilog,25) ' natively supported by CAMPARI :',other_natfreq*(1.0-other_unkfreq)
                  write(ilog,*) 'Fraction of all OTHER moves operating on unsupported '
                  write(ilog,25) ' d.o.f.s in native CAMPARI res.:',(1.0-other_natfreq)*(1.0-other_unkfreq)
                end if
                write(ilog,*) 'Fraction of full randomization'
                write(ilog,25) ' attempts amongst OTHER moves  : ',other_randfreq 
                if (other_randfreq.lt.1.0) then
                  write(ilog,25) ' Step size for local OTHER m.  : ',other_stepsz
                end if
              end if
              if (otherfreq.lt.1.0) then
                if (use_globmoves.EQV..true.) then
                  write(ilog,*) 'Fraction of chi-coupled phi/psi-'
                  write(ilog,25) ' moves amongst bb pivot moves  : ',(1.0-puckerfreq)*(1.0-nucfreq)
                else
                  write(ilog,*) 'Fraction of regular (uncoupled) phi/psi-'
                  write(ilog,25) ' moves amongst bb pivot moves  : ',(1.0-puckerfreq)*(1.0-nucfreq)
                end if
                write(ilog,*) 'Fraction of full randomization'
                write(ilog,25) ' attempts amongst pivot moves  : ',pivot_randfreq 
                if (pivot_randfreq.lt.1.0) then
                  write(ilog,25) ' Step size for local pivot m.  : ',pivot_stepsz
                end if
              end if
            end if
          end if
        end if
      end if
    end if
  end if
  if (use_stericgrids.EQV..true.) then
    call strlims(griddir,s3,s4)
    write(ilog,*) 'Grids read from directory     : ',griddir(s3:s4)
    write(ilog,*) 'Resolution of grids           : '&
 &                  ,2*stgr%halfwindow
  end if
  end if
  end if
!
  write(ilog,*)
  write(ilog,*) '**OUTPUT AND ANALYSIS PARAMETERS:'
  write(ilog,*) 'Restart file output interval       : ',rstout
  write(ilog,*) 'Trajectory output interval         : ',xyzout
  if (xyzout.le.nsim) then
    if (xyzmode.eq.1) then
      write(ilog,*) 'Structures only written in Tinker arc(xyz)-format.'
    else if (xyzmode.eq.2) then
      write(ilog,*) 'Structures only written in pdb-format.'
    else if (xyzmode.eq.3) then
      write(ilog,*) 'Structures only written in (binary) dcd-format.'
    else if (xyzmode.eq.4) then
      write(ilog,*) 'Structures only written in (binary) xtc-format.'
    else if (xyzmode.eq.5) then
      write(ilog,*) 'Structures only written in (binary) nc-format (NetCDF).'
    end if
    write(ilog,*) 'PDB file write-out format conv.    : ',pdb_convention(1)
    if (use_trajidx.EQV..true.) then
      call strlims(trajidxfile,s3,s4)
      write(ilog,*) 'Only a user-selected subset of atoms is written to structural output.'
      write(ilog,*) '(Selected ',pdbeffn,' atoms via input file: ',trajidxfile(s3:s4),').'
    else if (just_solutes.EQV..true.) then
      write(ilog,*) 'Structural output suppressed for molecules tagged as solvent.'
    end if
    if ((pdb_force_box.gt.1).AND.(bnd_type.eq.1)) then
      write(ilog,*) 'Molecules are split on a per-atom basis when crossing periodic boundaries.'
    else if (bnd_type.eq.1) then
      write(ilog,*) 'Molecules are always left intact when crossing periodic boundaries.'
      if (pdb_rmol.gt.0) then
        write(ilog,*) 'Molecules are written as nearest images with respect to molecule #',pdb_rmol,'.'
      else
        write(ilog,*) 'The image whose centroid is in the box will be written for every molecule.'
      end if
    end if
  end if
  if (pdb_analyze.EQV..true.) then
    if (align%yes.EQV..true.) then
      write(ilog,*) 'Alignment interval                 : ',align%calc
    else
      write(ilog,*) 'Alignment interval                 : ',nsim+1
    end if
    if ((align%yes.EQV..true.).AND.(align%calc.le.nsim)) then
      call strlims(align%filen,t1,t2)
      write(ilog,*) 'File with atom index list          : ',align%filen(t1:t2)
      if (align%refset.EQV..true.) then
        call strlims(pdbtmplfile,t1,t2)
        write(ilog,*) 'Reference structure                : ',pdbtmplfile(t1:t2)
      else
        write(ilog,*) 'Alignment procedure always uses prior aligned frame as reference.'
      end if
      if (align%instrmsd.EQV..true.) then
        write(ilog,*) 'Requested instantaneous RMSD output.'
        if (align%diffnr.eq.0) then
          write(ilog,*) 'Difference set for RMSD is the same as alignment set.'
        else
          write(ilog,*) '  Number of atoms in difference set: ',align%diffnr
          call strlims(cfilen,t1,t2)
          write(ilog,*) '  Source file                      : ',cfilen(t1:t2)
        end if
      else
        write(ilog,*) 'No instantaneous RMSD output was requested.'
      end if
    end if
  end if
  write(ilog,*) 'Energy output interval             : ',enout
  if (dyn_mode.ne.1) then
    write(ilog,*) 'Ensemble output interval           : ',ensout
  end if
  write(ilog,*) 'Acceptance output interval         : ',accout
  write(ilog,*) 'Torsions output interval           : ',torout
  write(ilog,*) 'Sys. polymer stat. output interval : ',polout
  write(ilog,*) 'Angular distribution calc. interval: ',angcalc
  if (angcalc.le.nsim) then
    write(ilog,50) ' Resolution in deg.      : ',fyres
    if (nrspecrama.gt.0) then
      write(ilog,*) 'Nr. of specif. residues requested  : ',nrspecrama
    end if
    if (nrmolrama.gt.0) then
      write(ilog,*) 'Nr. of spec. analysis grp.s reques.: ',nrmolrama
    end if
  end if 
  write(ilog,*) 'Intern. coord. stat. calc. interval: ',intcalc
  if ((do_ints(1).EQV..true.).AND.(intcalc.le.nsim)) then
    write(ilog,*) 'Bond length analysis requested.'
    write(ilog,*) 'Number of terms         : ',sum(nrsbl(1:nseq))
    write(ilog,50) ' Resolution in Ang.      : ',intres(1)
  end if
  if ((do_ints(2).EQV..true.).AND.(intcalc.le.nsim)) then
    write(ilog,*) 'Bond angle analysis requested.'
    write(ilog,*) 'Number of terms         : ',sum(nrsba(1:nseq))
    write(ilog,50) ' Resolution in deg.      : ',intres(2)
  end if
  if ((do_ints(3).EQV..true.).AND.(intcalc.le.nsim)) then
    write(ilog,*) 'Improper torsion angle analysis requested.'
    write(ilog,*) 'Number of terms         : ',sum(nrsimpt(1:nseq))
    write(ilog,50) ' Resolution in deg.      : ',intres(3)
  end if
  if ((do_ints(4).EQV..true.).AND.(intcalc.le.nsim)) then
    write(ilog,*) 'Proper dihedral angle analysis requested.'
    write(ilog,*) 'Number of terms         : ',sum(nrsdi(1:nseq))
    write(ilog,50) ' Resolution in deg.      : ',intres(4)
  end if
  write(ilog,*) 'BB-segment distrib. calc. interval : ',segcalc
  if ((use_ZSEC.EQV..false.).AND.(segcalc.le.nsim)) then
    write(ilog,22) ' Center of alpha-basin    : ',&
 &                           par_ZSEC2(1),par_ZSEC2(2)
    write(ilog,20) ' Radius of alpha-basin    : ',&
 &par_ZSEC2(3)
    write(ilog,20) ' Steepness of alpha-basin : ',&
 &par_ZSEC2(4)
    write(ilog,22) ' Center of beta-basin     : ',&
 &                           par_ZSEC2(5),par_ZSEC2(6)
    write(ilog,20) ' Radius of beta-basin     : ',&
 &par_ZSEC2(7)
    write(ilog,20) ' Steepness of beta-basin  : ',&
 &par_ZSEC2(8)
  end if
  write(ilog,*) 'DSSP calculation interval          : ',dsspcalc
  if ((dsspcalc.le.nsim).AND.(inst_dssp.EQV..true.))  then
    write(ilog,*) 'Running DSSP output requested.'
  end if
  if (dsspcalc.le.nsim) then
    write(ilog,50) ' Energy for good H-bond   : ',&
 &par_DSSP(1)
    write(ilog,50) ' Worst H-bond energy      : ',&
 &par_DSSP(2)
    write(ilog,50) ' Best H-bond energy       : ',&
 &par_DSSP(3)
    write(ilog,50) ' Cutoff for H-bond search : ',&
 &par_DSSP(4)
    write(ilog,*) 'DSSP E/H-score mode      : ',&
 &par_DSSP2(2)
    write(ilog,*) 'DSSP E/H-score exponent  : ',&
 &par_DSSP2(1)
  end if
  write(ilog,*) 'Polymer stat. calculation interval : ',polcalc
  if (polcalc.le.nsim) then
    write(ilog,50) ' Bin-size for Rg-axis of histogram  : ',rg_binsz
  end if
  write(ilog,*) 'Hydrodynamic radius calc. interval : ',rhcalc
  write(ilog,*) 'Pair correlation calc. interval    : ',pccalc
  if (pccalc.le.nsim) then
    write(ilog,50) ' Bin-size in Angstrom    : ',pc_binsz
    write(ilog,*) 'Number of unique PCs    : ',gpc%nos 
    if (gpc%nos.gt.0) then
      if (inst_gpc.gt.0) then
        write(ilog,*) 'Relative output interval for instantaneous values for '
        write(ilog,*) 'selected distances      : ',inst_gpc
      end if
    end if
  end if
  write(ilog,*) 'Contact analysis calc. interval    : ',contactcalc
  if (contactcalc.lt.nsim) then
    write(ilog,*) 'Contact analysis offset : ',contact_off
    write(ilog,55) ' Contact c.o.m cutoff    : ',sqrt(contact_cuts(2))
    write(ilog,55) ' Contact minimum cutoff  : ',sqrt(contact_cuts(1))
    write(ilog,*) 'Cluster calc. interval  : ',clucalc
  end if
  write(ilog,*) 'Dipole analysis calc. interval     : ',dipcalc
!  write(ilog,*) 'LC of torsions calc. interval      : ',torlccalc
!  if (torlccalc.lt.nsim) then
!    write(ilog,*) 'Number of torsional LCs : ',ntorlcs
!    write(ilog,55) ' Resolution for tor. LCs : ',torlc_params(2)
!  end if
  write(ilog,*) 'Covariance analysis calc. interval : ',covcalc
  if (covcalc.lt.nsim) then
    write(ilog,*) 'Mode for cov. analysis  : ',covmode
  end if
  write(ilog,*) 'Solvent-acc. volume calc. interval : ',savcalc
  if (savcalc.le.nsim) then
    write(ilog,*) 'Unique atom reqests     : ',savreq%nats
    if (use_IMPSOLV.EQV..false.) then
        write(ilog,55) ' SAV probe radius        : ',savprobe
        write(ilog,21) ' Sphere overlap calc mode: ',par_IMPSOLV2(1)
      if (par_IMPSOLV2(2).eq.1) then
        write(ilog,*) 'Using smooth sigmoidal function for SAV-to-solvation state mapping (analysis).'
      else if (par_IMPSOLV2(2).eq.2) then
        write(ilog,*) 'Using stair-stepped sigmoidal function for SAV-to-solvation state mapping (analysis).'
      end if
      if ((par_IMPSOLV2(2).eq.1).OR.(par_IMPSOLV2(2).eq.2)) then
        write(ilog,20) ' Sigmoidal decay constant: ',par_IMPSOLV(3)
        write(ilog,20) ' Rel. center of sigm. fxn: ',par_IMPSOLV(6)
      end if
      if (par_IMPSOLV2(2).eq.2) then
        write(ilog,20) ' Stair step granularity  : ',par_IMPSOLV(12)
        write(ilog,20) ' Stair step steepness    : ',par_IMPSOLV(13)
        write(ilog,20) ' Stair step position     : ',par_IMPSOLV(14)
      end if
    end if
    write(ilog,*) 'Rel. inst. output freq. : ',savreq%instfreq
  end if
  write(ilog,*) 'Scattering calculation interval    : ',sctcalc
  if ((rhcalc.lt.nsim).AND.(sctcalc.lt.(nsim/rhcalc))) then
    write(ilog,*) 'Number of wave vectors  : ',nqv
    write(ilog,55) ' Res. for wave vectors   : ',qv_res
  end if
  write(ilog,*) 'Hole analysis calculation interval : ',holescalc
  write(ilog,*) 'Diffraction calculation interval   : ',diffrcalc
  if (diffrcalc.le.nsim) then
    write(ilog,*) 'Number of radial bins   : ',rr_max
    write(ilog,*) 'Number of axial bins    : ',zz_max
    write(ilog,55) ' Res. in radial dim.     : ',rr_res
    write(ilog,55) ' Res. in axial dim.      : ',zz_res
    write(ilog,*) 'Max. order of Bessel fxn: ',bes_max
    write(ilog,55) ' Res. for Bessel fxn.s   : ',bes_res
    if (diffr_uax.EQV..true.) then
      write(ilog,9) ' Cylindrical axis (x,y,z): ',diffr_axis(1),&
 &diffr_axis(2),diffr_axis(3)
      write(ilog,9) ' Axis origin (x,y,z)     : ',diffr_axon(1),&
 &diffr_axon(2),diffr_axon(3)
    else
      write(ilog,*) 'Using longest cylindrical axis as defined by sy&
 &stem.'
    end if
  end if
  write(ilog,*) 'Spatial density map calc. interval : ',emcalc
  if (emcalc.le.nsim) then
    if (empotprop.eq.1) then
      write(ilog,*) 'Computing spatial density based on atomic masses.'
    else if (empotprop.eq.2) then
      write(ilog,*) 'Computing spatial density based on atomic number.'
    else if (empotprop.eq.3) then
      write(ilog,*) 'Computing spatial density based on atomic charge.'
    else if (empotprop.eq.4) then
      call strlims(empropfile,t1,t2)
      write(ilog,*) 'Computing spatial density based on user-defined property.'
      write(ilog,*) 'Input file with propert.: ',empropfile(t1:t2)
    end if
    write(ilog,55) ' Background density      : ',embgdensity
    write(ilog,21) ' B-spline order          : ',emsplor
    if (.NOT.((bnd_shape.eq.1).AND.(bnd_type.eq.1))) write(ilog,55) ' Buffer scale factor     : ',embuf
    write(ilog,9) ' Spatial resol.  (x,y,z) : ',emgrid%deltas(1:3)
    write(ilog,9) ' Map origin      (x,y,z) : ',emgrid%origin(1:3)
    write(ilog,8) ' Grid dimensions (x,y,z) : ',emgrid%dim(1:3)
  end if
  if (particleflucfreq.gt.0.0) then
    write(ilog,*) 'GC number hist. calc. interval     : ',particlenumcalc
  else
    write(ilog,*) 'Subvol. number hist. calc. interval: ',particlenumcalc
  end if
  if (particlenumcalc.le.nsim) then
    write(ilog,*) 'Inst. number output (rel. interval): ',particlenuminst
  end if
  write(ilog,*) 'Struct. clust. collection interval : ',cstorecalc
  if (((cstorecalc.le.nsim).OR.(cmode.eq.7)).AND.(use_MPIMultiSeed.EQV..false.)) then
    if (cmode.ne.7) then
      write(ilog,*) 'Distance criterion index  : ',cdis_crit
      if (cprepmode.eq.0) then
        write(ilog,*) 'No signal preprocessing requested for structural clustering.'
      else if (cprepmode.eq.1) then
        write(ilog,*) 'Requested data for structural clustering to be centered.'
      else if (cprepmode.eq.2) then
        write(ilog,*) 'Requested data for structural clustering to be centered and normalized by standard deviation.'
      else if (cprepmode.eq.3) then
        write(ilog,*) 'Requested data for structural clustering to be smoothed.'
      else if (cprepmode.eq.4) then
        write(ilog,*) 'Requested data for structural clustering to be centered and smoothed.'
      else if (cprepmode.eq.5) then
        write(ilog,*) 'Requested data for structural clustering to be centered, normalized by standard deviation, and smoothed.'
      end if
      if (cprepmode.ge.3) write(ilog,*) 'Smoothing window (snaps.) : ',csmoothie
      if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.((cdis_crit.ge.8).AND.(cdis_crit.le.10))).AND.(cchangeweights.gt.0)) then
        write(ilog,*) 'Mode for altered weights  : ',cchangeweights
        write(ilog,*) 'Combin. rule for weights  : ',cwcombination
        if ((cchangeweights.ne.2).AND.(cdis_crit.ne.8)) write(ilog,*) 'Window size (snapshots)   : ',cwwindowsz
        if (((cchangeweights.ge.4).AND.(cchangeweights.le.9)).AND.(cdis_crit.ne.8)) then
          write(ilog,50) ' Buffer parameter for inv. : ',cdynwbuf
        end if
        if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
          write(ilog,*) 'Lag time for ACF (snaps.) : ',cwacftau
        end if
        if ((cprepmode.lt.3).AND.((cchangeweights.eq.5).OR.(cchangeweights.eq.7).OR.(cchangeweights.eq.9))) then
          write(ilog,*) 'Smoothing window (snaps.) : ',csmoothie
        end if
      end if
    end if
    if (cfeature_dump.EQV..true.) then
      write(ilog,*) 'Requested writing of the features extracted for clustering to a dedicated output file.'
    else
#ifdef LINK_NETCDF
      write(ilog,*) 'No writing of the features extracted for clustering requested.'
#endif
    end if
    if (cmode.eq.1) then
      write(ilog,*) 'Using simple leader-based clustering algorithm (does not require neighbor list).'
    else if (cmode.eq.2) then
      write(ilog,*) 'Using modified leader-based clustering algorithm (does not require neighbor list).'
    else if (cmode.eq.3) then
      write(ilog,*) 'Using hierarchical clustering algorithm (requires neighbor list).'
    else if ((cmode.eq.4).AND.(cprogindex.eq.1)) then
      write(ilog,*) 'Using cut-based one-shot clustering algorithm with exact progress sequence (requires neighbor list).'
    else if ((cmode.eq.4).AND.(cprogindex.eq.2)) then
      write(ilog,*) 'Using cut-based one-shot clustering algorithm with approximate progress sequence (does not require neighbor &
 &list).'
    else if (cmode.eq.5) then
      write(ilog,*) 'Using BIRCH-like clustering algorithm (does not require neighbor list).'
    else if (cmode.eq.6) then
      write(ilog,*) 'Reading back clustering from file (does not require neighbor list) with data read-in.'
    else if (cmode.eq.7) then
      write(ilog,*) 'Reading back clustering from file (does not require neighbor list) without data read-in (limited options).'
    end if
    if ((cmode.eq.3).OR.((cmode.eq.4).AND.(cprogindex.eq.1))) then
      write(ilog,50) ' Hard cutoff (var. units)  : ',chardcut
      write(ilog,50) ' Threshold for init. screen: ',cmaxrad
    end if
    if ((read_nbl_from_nc.EQV..false.).OR.(cmode.le.2)) then
      write(ilog,*) 'Leader directional flag   : ',cleadermode
    end if
    write(ilog,50) ' Thres. radius (var. units): ',cradius
    if ((cmode.eq.5).OR.((cmode.eq.4).AND.(cprogindex.eq.2))) then
      write(ilog,50) ' Coarsest radius crit.     : ',cmaxrad
#ifdef ENABLE_THREADS
      write(ilog,50) ' Merging parameter (rel.)  : ',cmergecrit
      write(ilog,*) 'Factor for parallel exec. : ',cbirchbatch
#endif
    end if
    if ((read_nbl_from_nc.EQV..true.).AND.((cmode.eq.3).OR.((cmode.eq.4).AND.(cprogindex.eq.1)))) then
      call strlims(nblfilen,t1,t2)
      write(ilog,*) 'Neighbor list NetCDF file : ',nblfilen(t1:t2)
    end if
    if (cmode.ne.7) then
      call strlims(cfilen,t1,t2)
      write(ilog,*) 'Input file (if existing)  : ',cfilen(t1:t2)
      if (cdis_crit.eq.1) then
        write(ilog,*) 'Number of torsions in set : ',calcsz
      else if (cdis_crit.eq.2) then
        write(ilog,*) 'Number of torsions in set : ',calcsz/2
      else if (cdis_crit.eq.3) then
        write(ilog,*) '# of Fourier terms in set : ',calcsz
      else if (cdis_crit.eq.4) then
        write(ilog,*) '# of Fourier terms in set : ',2*calcsz/3
      else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
        write(ilog,*) 'Number of atoms in set    : ',(cdofsbnds(2)-cdofsbnds(1)+1)/3
      else if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
        write(ilog,*) 'Number of distances in set: ',clstsz
        if (cdistransform.eq.0) then
          write(ilog,*) 'Distances are used as is.'
        else if (cdistransform.eq.1) then
          write(ilog,*) 'Distances are transformed (sigmoidal).'
          write(ilog,50) ' Mid-point parameter : ',cdistrans_params(1)
          write(ilog,50) ' Width parameter     : ',cdistrans_params(2)
        else if (cdistransform.eq.2) then
          write(ilog,*) 'Distances are transformed (hyperbolic).'
          write(ilog,50) ' Buffer parameter    : ',cdistrans_params(1)
          write(ilog,50) ' Decay parameter     : ',cdistrans_params(2)
        else if (cdistransform.eq.3) then
          write(ilog,*) 'Distances are transformed (sine).'
          write(ilog,50) ' Tapering threshold  : ',cdistrans_params(1)
        end if
      end if
      if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
        if (align_for_clustering.EQV..true.) then
          write(ilog,*) 'Structures will be aligned before computing distance.'
          if ((cdis_crit.eq.6).AND.((cdofsbnds(4).ne.cdofsbnds(2)).OR.(cdofsbnds(3).ne.cdofsbnds(1)))) then
            call strlims(align%filen,t1,t2)
            write(ilog,*) 'File with alignment set   : ',align%filen(t1:t2)
            write(ilog,*) 'Number of atoms in set    : ',(cdofsbnds(4)-cdofsbnds(3)+1)/3
            write(ilog,*) 'Number of at. in joint set: ',calcsz/3
          end if
        else
          write(ilog,*) 'Structures will NOT be aligned before computing distance.'
        end if
        if (pdb_rmol.gt.0) then
          write(ilog,*) 'Coordinates are taken as nearest images with respect to molecule #',pdb_rmol,'.'
        else
          write(ilog,*) 'Coordinates are taken such that centroid of every molecule is in the box.'
        end if
      end if
    end if
    if ((cmode.eq.2).OR.(cmode.eq.5).OR.((cmode.eq.4).AND.(cprogindex.eq.2))) then
      if (refine_clustering.EQV..true.) then
        write(ilog,'(1x,a)') 'Requested refinement of clustering results. This is not recommended and may compromise performance.'
      end if
    end if
    if ((cmode.ne.6).AND.(cmode.ne.7).AND.(cprogindex.ne.1)) then
      if (resort_clustering.EQV..true.) then
        write(ilog,'(1x,a)') 'Clusters with identical sizes will be resorted by the index of their centroid/origin &
 &representatives.'
      else
        write(ilog,'(1x,a)') 'Ties in cluster size will not be reprocessed for sorting.'
#ifdef ENABLE_THREADS
        if ((cmode.eq.5).OR.((cmode.eq.4).AND.(cprogindex.eq.2)).AND.(cbirchbatch.ge.1.0).AND.(thrdat%maxn.gt.1)) then
          write(ilog,'(1x,a)') 'Results from stable, threaded tree-based clustering will look different because of this &
 &even though they are not.' 
        end if
#endif
      end if
    end if
    if (cmode.eq.3) then
      write(ilog,*) 'Linkage criterion         : ',clinkage
    else if (cmode.eq.4) then
      if (cprogindstart.eq.0) then
        write(ilog,*) 'Multiple profiles generated from guessed maximum cut (basin) snapshots.'
        write(ilog,*) 'Min minimum search width  : ',csivmin
        write(ilog,*) 'Max minimum search width  : ',csivmax
      else
        write(ilog,*) 'Ref. snapshot for profile : ',cprogindstart
      end if
      if (cprogindex.eq.2) then
        write(ilog,*) 'Number of levels in tree  : ',c_nhier
        write(ilog,*) 'Multi-resolution depth    : ',c_multires
        write(ilog,*) 'Maximum search attempts   : ',cprogindrmax
        write(ilog,*) 'Auxiliary search depth    : ',cprogrdepth
        write(ilog,*) 'Batch size (random bran.) : ',cprogbatchsz
      end if
      write(ilog,*) 'Max. local partition width: ',cprogpwidth
      write(ilog,*) 'No. of MST edge folds     : ',cprogfold
    else if (cmode.eq.5) then
      write(ilog,*) 'Number of levels in tree  : ',c_nhier
      write(ilog,*) 'Multi-resolution depth    : ',c_multires
    else if ((cmode.eq.6).OR.(cmode.eq.7)) then
      call strlims(clufilen,t1,t2)
      write(ilog,*) 'File with cluster traject.: ',clufilen(t1:t2)
      write(ilog,*) 'Column to read            : ',clufcol 
    end if
    if (cmode.ne.7) then
      if (pcamode.eq.1) then
        write(ilog,*) 'No computation of principal components requested.'
      else if (pcamode.eq.2) then
        write(ilog,*) 'Computation of principal components requested (eigenvectors only).'
      else if (pcamode.eq.3) then
        write(ilog,*) 'Computation of principal components requested (including transformed data).'
      else if (pcamode.eq.4) then
        write(ilog,*) 'Computation of time structure-based independent components requested (eigenvectors only).'
      else if (pcamode.eq.5) then
        write(ilog,*) 'Computation of time structure-based independent components requested (including transformed data).'
      end if
      if ((pcamode.eq.4).OR.(pcamode.eq.5)) then
        write(ilog,*) 'Lag time for ACF (snaps.) : ',cwacftau
      end if
      if ((pcamode.eq.3).OR.(pcamode.eq.5)) then
        if (reduced_dim_clustering.gt.0) then
          write(ilog,*) 'Clustering algorithm will be run on a subset of transformed components.'
          if (cdis_crit.eq.4) then
            write(ilog,*) 'No. of components to use  : ',reduced_dim_clustering-mod(reduced_dim_clustering,2)
          else
            write(ilog,*) 'No. of components to use  : ',reduced_dim_clustering
          end if
        end if
      else
        if (reduced_dim_clustering.gt.0) then
          write(ilog,*) 'Data dimensions will be discarded at the end per user request.'
          write(ilog,*) 'No. of dimensions to use  : ',reduced_dim_clustering        
        end if
      end if
    end if
    if (ntbrks.gt.0) then
      call strlims(tbrkfilen,t1,t2)
      write(ilog,*) 'File with traject. breaks : ',tbrkfilen(t1:t2)
    end if
    if (ntlnks.gt.0) then
      call strlims(tlnkfilen,t1,t2)
      write(ilog,*) 'File with add. traj. links: ',tlnkfilen(t1:t2)
    end if
    if (re_aux(3).eq.2) then
      call strlims(re_traceinfile,s3,s4)
      write(ilog,*) 'Trace file in use               : ',re_traceinfile(s3:s4)
      write(ilog,54) ' Assumed equilibration period   : ',re_aux(5)
      write(ilog,54) ' Assumed trajectory saving freq.: ',re_aux(4)
      write(ilog,54) ' Assumed snapshots per replica  : ',re_aux(8)
    end if
    if (clagt_msm.eq.1) then
      write(ilog,*) 'Snapshots in memory are assumed to be connected in sequence by default.'
    else
      write(ilog,*) 'Snapshots in memory are assumed to be connected with a time lag by default.'
      write(ilog,54) 'Lag time in stored snapshots    : ',clagt_msm
    end if
    if (caddlkmd.ne.0) then
      write(ilog,*) 'Mode for add./remov. links: ',caddlkmd
!      write(ilog,50) ' Thresh. for added links   : ',caddlinks
      if ((abs(caddlkmd).le.3).OR.(abs(caddlkmd).eq.5)) then
        write(ilog,71) ' Weight of added links     : ',caddlinkwt
      end if
    else
      write(ilog,*) 'No cluster-based links (edges) are added to network (graph).'
    end if
    if ((cequil.gt.0).OR.(ccfepmode.ne.0).OR.(synmode.ne.0)) then
      write(ilog,71) ' Time cutoff for iterat. algorithms (if any) [h]: ',maxtime_eqmsm/(3600.*thrdat%rate)
    end if
    if (cequil.eq.1) then
      write(ilog,*) 'Requested network-based reweighting via steady state of MLE Markov model.'
    else if (cequil.eq.2) then
      write(ilog,*) 'Requested network-based reweighting via steady state of diffusion-based Markov model.'
    else if (cequil.eq.3) then
      write(ilog,*) 'Requested network-based reweighting via flat propagator method.'
    else if (cequil.eq.4) then
      write(ilog,*) 'Requested network-based reweighting via flat propagator method and diffusion information.'
    end if
    if (cequil.gt.0) then
      write(ilog,*) 'This implies independent treatment of strongly connected components.'
      write(ilog,50) ' Buffer for snapshot wts.  : ',cequilbuf
    else
      write(ilog,*) 'Complete graph is analyzed independently of connectedness. This may cause errors.'
    end if
    if (ccfepmode.eq.0) then
      write(ilog,*) 'No computation of cut-based pseudo free energy profiles performed.'
    else if (ccfepmode.eq.1) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using iterative MFPT values as order parameter.'
      if (inissnap_usrslct.gt.0) then
        write(ilog,*) 'Ref. snapshot for profile : ',inissnap_usrslct
      else if (inissnap_usrslct.eq.0) then
        write(ilog,*) 'Reference snapshot is centroid of largest cluster.'
      else if (inissnap_usrslct.eq.-1) then
        write(ilog,*) 'Reference snapshot are centroids of largest clusters per component.'
      else if (inissnap_usrslct.eq.-2) then
        write(ilog,*) 'Reference snapshot is centroid of largest cluster in largest strongly connected component.'
      end if
    else if (ccfepmode.eq.8) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using (+) committors as order parameter.'
    else if (ccfepmode.eq.9) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using (-) committors as order parameter.'
    else if (ccfepmode.eq.10) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using both (+) and (-) committors as order &
 &parameter.'
    else if (ccfepmode.eq.11) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using edge outward flux as order parameter.'
    else if (ccfepmode.eq.12) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using edge inward flux as order parameter.'
    end if
#ifdef LINK_HSL
    if (eigvalmd.ne.0) then
      write(ilog,'(A58)') ' Spectral analysis of transition matrix(ces) is requested.'
      write(ilog,56) 'Eigenvalues selection mode                              : ',eigvalmd
      write(ilog,56) 'Number of eigenvalues requested                         : ',numeig
      if (arnblks.ne.-1) then
        write(ilog,56) 'Number of blocks in Arnoldi method                      : ',arnblks
      else
        write(ilog,56) 'Number of blocks in Arnoldi method                      : ',numeig + 2  !set in sanity
      end if
      if (arnblks.ne.-1) then
        write(ilog,56) 'Number of steps in Arnoldi method                       : ',arnstps  
      else 
        write(ilog,56) 'Number of steps in Arnoldi method                       : ',ceiling((8.*numeig)/numeig + 2)  !set in sanity 
      end if
      write(ilog,56) 'Number of max restarts in Arnoldi method                : ',arnmaxitr
      write(ilog,56) 'Mode for transition matrix(ces) build-up                : ',which_timeflow
      write(ilog,56) 'Reference snapshot for spectral analysis                : ',inissnap_usrslct
      if (arntol.gt.0) then
        write(ilog,74) ' Tolerance for Arnoldi method                            : ',arntol
      else
        write(ilog,'(A82)') ' Tolerance for Arnoldi method                          : "machine precision"*10^3.'
      end if
      if (doeigvec_msm.EQV..true.) then
        write(ilog,'(A42)') ' Computation of eigenvectors is requested.'
      else
        write(ilog,'(A42)') ' No computation of eigenvectors requested.'
      end if
    else
      write(ilog,'(A53)') ' Computation of spectral properties is not requested.'
    end if
    if (dopfold.EQV..true.) then
      write(ilog,'(A59)') ' Computation of committor probability (pfold) is requested.'
      if (pfold_report.EQV..true.) then
        write(ilog,'(A59)') ' Report files for pfold-specific properties will be written.'
      end if
      call strlims(clufoldfile,t1,t2)
      write(ilog,5) ' File with folded set (if existent)  : ',clufoldfile(t1:t2)  
      call strlims(cluunfoldfile,t1,t2)
      write(ilog,5) ' File with unfolded set (if existent): ',cluunfoldfile(t1:t2)
      if (dopfoldminus.EQV..true.) then
        write(ilog,'(A60)') ' Computation of backward committor probability is requested.'
        write(ilog,'(A44)') ' Steady state will be computed if necessary.'
        if (eigvalmd.eq.0) then
          if (arnblks.ne.-1) then
            write(ilog,*) 'Number of blocks in Arnoldi method       : ',arnblks
          else
            write(ilog,*) 'Number of blocks in Arnoldi method       : ',3  !set in graph_alg
          end if
          if (arnblks.ne.-1) then
            write(ilog,*) 'Number of steps in Arnoldi method        : ',arnstps  
          else 
            write(ilog,*) 'Number of steps in Arnoldi method        : ',3  !set in graph_alg
          end if
            write(ilog,*) 'Number of max restarts in Arnoldi method : ',max(10,arnmaxitr)
        end if
      else
        write(ilog,*) 'No computation of backward committor probability requested.'
      end if
    else
      write(ilog,*) 'Computation of committor probability (pfold) is not requested.'
    end if
#endif
    if (synmode.ne.0) then
      write(ilog,'(A57)') ' Random walks (synthetic trajectories) will be generated.'
      write(ilog,49) 'Synthetic trajectory mode                              : ',synmode
      write(ilog,49) 'Target numb. of synthetic trajectories                 : ',nstraj
      write(ilog,49) 'Output frequency for synthetic trajectories            : ',nsskip
      if (inissnap_usrslct.gt.0) then
        write(ilog,49) 'Component for random walker hosts snapshot             : ',inissnap_usrslct
      else if (inissnap_usrslct.eq.-1) then
        write(ilog,'(A110)') ' Component for random walker                         : The one that hosts the centroid &
 &of the largest cluster.'
      else if (inissnap_usrslct.eq.-2) then
        write(ilog,'(A67)') ' Component for random walker                         : largest SCC.'
      end if 
      if (synmode.eq.1) then
        write(ilog,'(A95)') ' Start snapshot of synth. traj.                      : centroid of largest cluster in ref. SCC.'
        if (endssnap_usrslct.ne.0) then
          write(ilog,49) 'End (target) snaphot of synthetic trajectories         : ',endssnap_usrslct
        else
          write(ilog,'(A95)') ' End (target) snaphot of synth. traj.                   : the last one processed by clustering.'
        end if
        write(ilog,49) 'Max snaps. per synthetic trajectory                    : ',nssnap
      else if ((synmode.eq.2).or.(synmode.eq.3)) then
        if (synmode.eq.2) then
          write(ilog,'(A95)') ' Start snapshot of synth. traj.                      : centroid of largest cluster in ref. SCC.'
        else if (synmode.eq.3) then
          write(ilog,'(A99)') ' Start snapshot of synth. traj.                     : randomly varied according to cluster weights.'
        end if
        write(ilog,'(A76,I8,A7)') ' End (target) cluster of synth. traj.                   : The one hit after ', &
 &nssnap,' steps.'
        write(ilog,49) 'Snapshots per synthetic trajectory                     : ',nssnap
      else
        write(ilog,'(A48)') ' Unrecognized synth. traj. mode: This is a bug.'
        call fexit()
      end if
    else
      write(ilog,'(A70)') ' Generation of random walks (synthetic trajectories) is not requested.'
    end if

  end if
  if (phfreq.gt.0.0) then
    write(ilog,*) 'Titration state writeout interval  : ',phout
  end if
  write(ilog,*)
!  if (use_inherent.EQV..true.) then
!    write(ilog,*) 'Inherent structure calc. interval  : ',ncmap
!    write(ilog,*) 'Inherent structure xyz-interval    : ',ncmap_xyz
!    write(ilog,*) 'Convergence criterion for minimizer: ',convcrit
!    write(ilog,*) 'Output interval for minimizer      : ',minreport
!  else
!    write(ilog,*) 'No inherent structure calculations requested.'
!  end if
!  write(ilog,*)
  write(ilog,*) '**PARALLELISM SETTINGS:'
#ifdef ENABLE_MPI
#ifdef ENABLE_THREADS
  write(ilog,*) 'This is a shared memory (OpenMP-based) calculation beneath an MPI layer.'
  write(ilog,21) ' Max. number of threads per MPI process: ',thrdat%maxn
  write(ilog,21) ' Verbosity level for threads diagnost. : ',thrdat%verbosity
  write(ilog,21) ' Frequency for reenabling DLB          : ',thrdat%dlbfreq
  write(ilog,21) ' Maximum length of DLB segment         : ',thrdat%dlblen
  write(ilog,21) ' Scale for DLB collection              : ',thrdat%dlbscale
#else
  write(ilog,*) 'This is an MPI calculation with every MPI process using exactly one thread.'
#endif
  write(ilog,53) 'I am node ',myrank+1,' in a universe of size ',&
 &mpi_nodes,'.'
  if (use_MPIcolls.EQV..true.) then
    write(ilog,*) 'Higher level standard MPI collective communication routines may be in use.'
  else
    write(ilog,*) 'Higher level collective communication (if any) is always handled by custom routines.'
    write(ilog,49)'Assumed comm. granule size      : ',mpi_granularity
  end if
  if (use_REMC.EQV..true.) then
    if (pdb_analyze.EQV..false.) then
      write(ilog,*) 'This is an MPI Replica Exchange calculation.'
    else
      write(ilog,*) 'This is a parallel trajectory analysis run using the MPI Replica Exchange framework.'
    end if
    write(ilog,54)'Number of replicas (conditions) : ',re_conditions
    write(ilog,54)'Number of dimensions of those   : ',re_conddim
    if (pdb_analyze.EQV..false.) then
      write(ilog,49)'Deterministic frequency         : ',re_freq
      write(ilog,54)'Number of swap attempts in cycle: ',re_tryswap
      if (inst_retr.EQV..true.) write(ilog,*)' Trace file with swap history is being written.'
      if (inst_retr.EQV..false.) write(ilog,*)' Trace file with swap history is NOT being written.'
      if (re_nbmode.eq.2) then
        write(ilog,*) 'Only neighboring conditions are swapped with.'
      else
        write(ilog,*) 'All conditions are swapped with.'
      end if
      if ((fycxyz.eq.2).OR.(force_rexyz.EQV..true.)) then
        write(ilog,*) 'Structures are swapped via Cartesian coordinates.'
      else
        write(ilog,*) 'Structures are swapped via internal degrees of freedom (rebuilt).'
      end if
    end if
!
    write(ilog,*) 'My specifications (some may be relevant only for other replicas):'
    do i=1,re_conddim
      if (re_types(i).eq.1) write(ilog,51) 'Dimension ',i,'&
 & (temperature)                     : ',re_mat(myrank+1,i)
      if (re_types(i).eq.2) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for IPP)          : ',re_mat(myrank+1,i)
      if (re_types(i).eq.3) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for attLJ)        : ',re_mat(myrank+1,i)
      if (re_types(i).eq.4) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for WCA)          : ',re_mat(myrank+1,i)
      if (re_types(i).eq.5) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for POLAR)        : ',re_mat(myrank+1,i)
      if (re_types(i).eq.6) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for IMPSOLV)      : ',re_mat(myrank+1,i)
      if (re_types(i).eq.7) write(ilog,51) 'Dimension ',i,'&
 & (implicit solvent dielectric)     : ',re_mat(myrank+1,i)
      if (re_types(i).eq.8) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for TOR)          : ',re_mat(myrank+1,i)
      if (re_types(i).eq.9) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for ZSEC)         : ',re_mat(myrank+1,i)
      if (re_types(i).eq.10) write(ilog,51) 'Dimension ',i,'&
 & (target alpha-content for ZSEC)   : ',re_mat(myrank+1,i)
      if (re_types(i).eq.11) write(ilog,51) 'Dimension ',i,'&
 & (target beta-content for ZSEC)    : ',re_mat(myrank+1,i)
      if (re_types(i).eq.12) write(ilog,52) 'Dimension ',i,'&
 & (screening model for impl. sol.)  : ',nint(re_mat(myrank+1,i))
      if (re_types(i).eq.13) write(ilog,51) 'Dimension ',i,'&
 & (FOS-tau for implicit solvent)    : ',re_mat(myrank+1,i)
      if (re_types(i).eq.14) write(ilog,51) 'Dimension ',i,'&
 & (SCR-tau for implicit solvent)    : ',re_mat(myrank+1,i)
      if (re_types(i).eq.15) write(ilog,51) 'Dimension ',i,'&
 & (FOS-midpoint for implicit solv.) : ',re_mat(myrank+1,i)
      if (re_types(i).eq.16) write(ilog,51) 'Dimension ',i,'&
 & (SCR-midpoint for implicit solv.) : ',re_mat(myrank+1,i)
      if (re_types(i).eq.17) write(ilog,51) 'Dimension ',i,'&
 & (contact dielectric of impl. sol.): ',re_mat(myrank+1,i)
      if (re_types(i).eq.18) write(ilog,52) 'Dimension ',i,'&
 & (generalized mean of impl. sol.)  : ',nint(re_mat(myrank+1,i))
      if (re_types(i).eq.19) write(ilog,51) 'Dimension ',i,'&
 & (mix. contr. for dist. dep. diel.): ',re_mat(myrank+1,i)
      if (re_types(i).eq.20) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for ghost-IPP)    : ',re_mat(myrank+1,i)
      if (re_types(i).eq.21) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for ghost-attLJ)  : ',re_mat(myrank+1,i)
      if (re_types(i).eq.22) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for ghost-POLAR)  : ',re_mat(myrank+1,i)
      if (re_types(i).eq.23) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for TABUL)        : ',re_mat(myrank+1,i)
      if (re_types(i).eq.24) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for POLY)         : ',re_mat(myrank+1,i)
      if (re_types(i).eq.25) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for DREST)        : ',re_mat(myrank+1,i)
      if (re_types(i).eq.26) write(ilog,51) 'Dimension ',i,'&
 & (scaling fac. for ghost-BONDED_B) : ',re_mat(myrank+1,i)
      if (re_types(i).eq.27) write(ilog,51) 'Dimension ',i,'&
 & (scaling fac. for ghost-BONDED_A) : ',re_mat(myrank+1,i)
      if (re_types(i).eq.28) write(ilog,51) 'Dimension ',i,'&
 & (scaling fac. for ghost-BONDED_I) : ',re_mat(myrank+1,i)
      if (re_types(i).eq.29) write(ilog,51) 'Dimension ',i,'&
 & (scaling fac. for ghost-BONDED_T) : ',re_mat(myrank+1,i)
      if (re_types(i).eq.30) write(ilog,51) 'Dimension ',i,'&
 & (target beta-content for DSSP)    : ',re_mat(myrank+1,i)
      if (re_types(i).eq.31) write(ilog,51) 'Dimension ',i,'&
 & (target alpha-content for DSSP)   : ',re_mat(myrank+1,i)
      if (dyn_mode.ne.1) then
        if (re_types(i).eq.32) write(ilog,51) 'Dimension ',i,'&
 & (timestep)                        : ',re_mat(myrank+1,i)
      else
        if (re_types(i).eq.32) write(ilog,51) 'Dimension ',i,'&
 & (dummy)                           :  N/A'
      end if
      if (re_types(i).eq.33) write(ilog,51) 'Dimension ',i,'&
 & (scaling factor for EMICRO)       : ',re_mat(myrank+1,i)
      if (re_types(i).eq.34) write(ilog,51) 'Dimension ',i,'&
 & (threshold density for EMICRO)    : ',re_mat(myrank+1,i)
    end do
    write(ilog,54)'Overlap calculation frequency   : ',re_olcalc
    if ((inst_reol.gt.0).AND.(inst_reol*re_olcalc.le.nsim)) then
      write(ilog,54)'Instantaneous write-out frequ.  : ',inst_reol*re_olcalc
    else
      write(ilog,*) 'No instantaneous values of foreign energies are&
 & reported.'
    end if
    if (reol_all.EQV..true.) then
      write(ilog,*) 'Overlap metrics are computed between all possible replica&
 & pairs (full matrix generated).'
    else
      write(ilog,*) 'Overlap metrics are computed only with the (at most) two&
 & neighboring replicas (only tridiagonal band matrix generated).'
    end if
    if ((dyn_mode.ne.1).AND.(pdb_analyze.EQV..false.)) then
      if (re_velmode.eq.1) then
        write(ilog,*) 'For REMD swaps, velocities are randomly re-assigned to the temperature for&
 & the current node (just like at the beginning of a dynamics run).'
      else if ((re_velmode.eq.2).AND.(Tthere.EQV..true.)) then
        write(ilog,*) 'For REMD swaps, velocities are re-scaled by the temperature increment&
 & defined by those velocities` origin condition (temperature).'
      else if ((re_velmode.eq.3).OR.((re_velmode.eq.2).AND.(Tthere.EQV..false.))) then
        write(ilog,*) 'For REMD swaps, velocities remain associated with structures at all times.'
      end if
    end if
    if (pdb_analyze.EQV..true.) then
      if (re_aux(3).eq.1) then
        write(ilog,*) 'This run is meant to operate on an unscrambled trajectory.'
        call strlims(re_traceinfile,s3,s4)
        write(ilog,*) 'Trace file in use               : ',re_traceinfile(s3:s4)
        write(ilog,54) ' Assumed equilibration period   : ',re_aux(5)
        write(ilog,54) ' Assumed trajectory saving freq.: ',re_aux(4)
      else
        write(ilog,*) 'This run operates exactly on the supplied trajectory.'
      end if
    end if
  else if (use_MPIAVG.EQV..true.) then
    if (force_singlexyz.EQV..false.) then
      write(ilog,*) 'Enforcing the generation of a composite (striped) trajectory file containing information from &
 &all replicas.'
    end if
    if ((do_wanglandau.EQV..true.).AND.(dyn_mode.eq.1)) then
      write(ilog,*) 'This is a parallel Wang-Landau calculation.'
    else if (use_MPIMultiSeed.EQV..true.) then
      if (pdb_analyze.EQV..false.) then
        write(ilog,*) 'This is an MPI PIGS calculation.'
        write(ilog,54) ' Number of protected replicas    : ',re_aux(7)
        write(ilog,54) ' Reseeding interval              : ',re_freq
        write(ilog,70) '  Samples per replica and interval: ',1.0*re_freq/(1.0*cstorecalc)
        if (inst_retr.EQV..true.) write(ilog,*)' Trace file with reseeding history is being written.'
        if (inst_retr.EQV..false.) write(ilog,*)' Trace file with reseeding history is NOT being written.'
        if (use_dyn.EQV..true.) then
          if (re_velmode.eq.1) then
            write(ilog,*) ' For reseedings, velocities are randomly reassigned (just like at the beginning of a dynamics run).'
          else
            write(ilog,*) ' For reseedings, velocities remain associated with structures at all times.'
          end if
        end if
      else
        write(ilog,*) 'This is an MPI analysis run emulating a PIGS stretch.'
        write(ilog,54) ' Number of protected replicas   : ',re_aux(7)
        write(ilog,54) ' Samples per replica            : ',cmaxsnaps
      end if
      if ((pdb_analyze.EQV..false.).OR.((pdb_analyze.EQV..true.).AND.(mpi_nodes*cmaxsnaps.gt.4))) then
        if (cprepmode.eq.0) then
          write(ilog,*) ' No signal preprocessing requested for data from PIGS stretches.'
        else if (cprepmode.eq.1) then
          write(ilog,*) ' Requested data from PIGS stretches to be centered.'
        else if (cprepmode.eq.2) then
          write(ilog,*) ' Requested data from PIGS stretches to be centered and normalized by standard deviation.'
        else if (cprepmode.eq.3) then
          write(ilog,*) ' Requested data from PIGS stretches to be smoothed.'
        else if (cprepmode.eq.4) then
          write(ilog,*) ' Requested data from PIGS stretches to be centered and smoothed.'
        else if (cprepmode.eq.5) then
          write(ilog,*) ' Requested data from PIGS stretches to be centered, normalized by standard deviation, and smoothed.'
        end if
        if (cprepmode.ge.3) write(ilog,*) ' Smoothing window (snaps.) : ',csmoothie
        if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.((cdis_crit.ge.8).AND.(cdis_crit.le.10))).AND.(cchangeweights.gt.0)) then
          write(ilog,*) ' Mode for altered weights  : ',cchangeweights
          write(ilog,*) ' Combin. rule for weights  : ',cwcombination
          if ((cchangeweights.ne.2).AND.(cdis_crit.ne.8)) write(ilog,*) ' Window size (snapshots)   : ',cwwindowsz
          if (((cchangeweights.ge.4).AND.(cchangeweights.le.9)).AND.(cdis_crit.ne.8)) then
            write(ilog,50) '  Buffer parameter for inv. : ',cdynwbuf
          end if
          if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
            write(ilog,*) ' Lag time for ACF (snaps.) : ',cwacftau
          end if
          if ((cprepmode.lt.3).AND.((cchangeweights.eq.5).OR.(cchangeweights.eq.7).OR.(cchangeweights.eq.9))) then
            write(ilog,*) ' Smoothing window (snaps.) : ',csmoothie
          end if
        end if
        write(ilog,*) ' PCA/tICA-based dimensionality reduction is currently not supported for the PIGS heuristic.'
        write(ilog,*) ' Distance criterion index  : ',cdis_crit
        write(ilog,*) ' Using cut-based one-shot clustering algorithm with approximate progress sequence (does not &
   &require neighbor list).'
        write(ilog,*) ' Leader directional flag   : ',cleadermode
        write(ilog,50) '  Thres. radius (var. units): ',cradius
        write(ilog,50) '  Coarsest radius crit.     : ',cmaxrad
#ifdef   ENABLE_THREADS
        write(ilog,50) '  Merging parameter (rel.)  : ',cmergecrit
        write(ilog,*) ' Factor for parallel exec. : ',cbirchbatch
#endif
        write(ilog,*) ' Number of levels in tree  : ',c_nhier
        write(ilog,*) ' Multi-resolution depth    : ',c_multires
        call strlims(cfilen,t1,t2)
        write(ilog,*) ' Input file (if existing)  : ',cfilen(t1:t2)
        if (cdis_crit.eq.1) then
          write(ilog,*) ' Number of torsions in set : ',calcsz
        else if (cdis_crit.eq.2) then
          write(ilog,*) ' Number of torsions in set : ',calcsz/2
        else if (cdis_crit.eq.3) then
          write(ilog,*) ' # of Fourier terms in set : ',calcsz
        else if (cdis_crit.eq.4) then
          write(ilog,*) ' # of Fourier terms in set : ',2*calcsz/3
        else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
          write(ilog,*) ' Number of atoms in set    : ',(cdofsbnds(2)-cdofsbnds(1)+1)/3
        else if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
          write(ilog,*) ' Number of distances in set: ',clstsz
        end if
        if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
          if (align_for_clustering.EQV..true.) then
            write(ilog,*) ' Structures will be aligned before computing distance.'
            if ((cdis_crit.eq.6).AND.((cdofsbnds(4).ne.cdofsbnds(2)).OR.(cdofsbnds(3).ne.cdofsbnds(1)))) then
              call strlims(align%filen,t1,t2)
              write(ilog,*) ' File with alignment set   : ',align%filen(t1:t2)
              write(ilog,*) ' Number of atoms in set    : ',(cdofsbnds(4)-cdofsbnds(3)+1)/3
              write(ilog,*) ' Number of at. in joint set: ',calcsz/3
            end if
          else
            write(ilog,*) ' Structures will NOT be aligned before computing distance.'
          end if
        end if
        write(ilog,*) ' Representative snapshot of largest cluster serves are reference snapshot for progress index.'
        write(ilog,*) ' Maximum search attempts   : ',cprogindrmax
        write(ilog,*) ' Auxiliary search depth    : ',cprogrdepth
        write(ilog,*) ' Batch size (random bran.) : ',cprogbatchsz
      end if
    else
      if (pdb_analyze.EQV..false.) then
        write(ilog,*) 'This is an MPI averaging calculation.'
      else
        write(ilog,*) 'This is an MPI averaging calculation in trajectory analysis mode.'
      end if
    end if
  end if
#else
#ifdef ENABLE_THREADS
  write(ilog,*) 'This is a shared memory (OpenMP-based) calculation without MPI layer.'
  write(ilog,21) ' Maximum number of threads to use      : ',thrdat%maxn
  write(ilog,21) ' Verbosity level for threads diagnost. : ',thrdat%verbosity
  write(ilog,21) ' Frequency for reenabling DLB          : ',thrdat%dlbfreq
  write(ilog,21) ' Maximum length of DLB segment         : ',thrdat%dlblen
  write(ilog,21) ' Scale for DLB collection              : ',thrdat%dlbscale
#else
  write(ilog,*) 'This is a single-CPU calculation (no parallelism of any kind).'
#endif
#endif
!
  write(ilog,*)
  if ((dyn_mode.eq.1).OR.(dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
    write(ilog,*) '**MONTE CARLO CUTOFF SETTINGS:'
    if (use_cutoffs.EQV..true.) then
      if (mccutmode.eq.1) then
        write(ilog,*) 'Short-range cutoff is applied at interatomic resolution.'
      else if (mccutmode.eq.2) then
        write(ilog,*) 'Short-range cutoff is applied at inter-residue resolution.'
      end if
    end if
    if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
      write(ilog,70) ' Grid-based cutoff (IPP/LJ)    : ',mcnb_cutoff
      if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
        write(ilog,70) ' Grid-based cutoff (EL/TAB)    : ',mcel_cutoff
      end if
      write(ilog,23) ' Grid origin     : ',grid%origin(1),' , ',grid%origin(2),' , ',grid%origin(3)
      write(ilog,24) ' Grid dimensions : ',grid%dim(1),' x ',grid%dim(2),' x ',grid%dim(3)
      write(ilog,23) ' Grid spacings   : ',grid%deltas(1),' , ',grid%deltas(2),' , ',grid%deltas(3)
      write(ilog,21) ' Max. grid-point neighbors  : ',grid%maxgpnbs
      write(ilog,21) ' Max. groups per grid-point : ',grid%maxnbs
    else if ((use_cutoffs.EQV..true.).AND.(use_rescrit.EQV..true.)) then
      write(ilog,70) ' Residue-based cutoff (IPP/LJ) : ',mcnb_cutoff
      if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
        write(ilog,70) ' Residue-based cutoff (EL/TAB) : ',mcel_cutoff
      end if
    else if (use_cutoffs.EQV..true.) then
      write(ilog,70) ' Simple NB cutoff (IPP/LJ)     : ',mcnb_cutoff
      if ((use_POLAR.EQV..true.).OR.(use_TABUL.EQV..true.)) then
        write(ilog,70) ' Simple cutoff (EL/TAB)        : ',mcel_cutoff
      end if
    else
      write(ilog,*) 'No cutoffs used.'
    end if
    if ((use_cutoffs.EQV..true.).AND.(use_POLAR.EQV..true.)) then
      if (lrel_mc.eq.1) then
        write(ilog,*)'Computing all monopole-monopole and monopole-dipole &
 &terms beyond (EL/TAB) cutoff explicitly in minimum image convention. This applies to all residue pairs, for which &
 &at least one of the two is flagged as carrying a net charge.'
      else if (lrel_mc.eq.2) then
        write(ilog,*) 'Computing all monopole-monopole &
 &terms beyond (EL/TAB) cutoff explicitly in minimum image convention. This applies to all residue pairs, for which &
 &both are flagged as carrying a net charge.'
      else if (lrel_mc.eq.3) then
        write(ilog,*) 'Computing all monopole-monopole terms beyond (EL/TAB) cutoff in minimum image &
 &convention at reduced resolution (point approximation). This applies to all residue pairs, for which &
 &both are flagged as carrying a net charge.'
      else if (lrel_mc.eq.4) then
        write(ilog,*) 'Using truncation for polar (and all other) interactions beyond (EL/TAB) cutoff.'
      end if
    end if
    if ((pdb_analyze.EQV..false.).AND.(use_stericscreen.EQV..true.)) then
      write(ilog,73) ' SR/bias MC screen (kcal/mol)  : ',screenbarrier
    end if
  end if
  write(ilog,*) 
  if (dyn_mode.ne.1) then
    write(ilog,*) '**DYNAMICS CUTOFF AND LOOP SETTINGS:'
    if (use_cutoffs.EQV..true.) then
      write(ilog,*) 'Neighbor list update frequency      : ',nbl_up
!      write(ilog,*) 'Maximum neighbor number             : ',nbl_maxnbs
    end if
    if ((use_cutoffs.EQV..true.).AND.(use_mcgrid.EQV..true.)) then
      write(ilog,70) ' Grid-based short-range cut.: ',mcnb_cutoff
      write(ilog,70) ' Grid-based mid-range cut.  : ',mcel_cutoff
      write(ilog,23) ' Grid origin     : ',grid%origin(1),' , ',grid%origin(2),' , ',grid%origin(3)
      write(ilog,24) ' Grid dimensions : ',grid%dim(1),' x ',grid%dim(2),' x ',grid%dim(3)
      write(ilog,23) ' Grid spacings   : ',grid%deltas(1),' , ',grid%deltas(2),' , ',grid%deltas(3)
      write(ilog,21) ' Max. grid-point neighbors  : ',grid%maxgpnbs
      write(ilog,21) ' Max. groups per grid-point : ',grid%maxnbs
    else if ((use_cutoffs.EQV..true.).AND.(use_rescrit.EQV..true.)) then
      write(ilog,70) ' Residue-based short-range cut.: ',mcnb_cutoff
      write(ilog,70) ' Residue-based mid-range cut.  : ',mcel_cutoff
    else if (use_cutoffs.EQV..true.) then
      write(ilog,70) ' Simple short-range cutoff  : ',mcnb_cutoff
      write(ilog,70) ' Simple mid-range cutoff    : ',mcel_cutoff
    else
      write(ilog,*) 'No cutoffs used.'
    end if
    if (lrel_md.eq.1) then
      if ((use_cutoffs.EQV..true.).AND.(use_POLAR.EQV..true.)) then
        write(ilog,*) 'Using truncation for polar (and all other) interactions beyond mid-range cutoff.'
      end if
    else if (lrel_md.eq.4) then
      if ((use_cutoffs.EQV..true.).AND.(use_POLAR.EQV..true.)) then
        write(ilog,*)'Computing all monopole-monopole terms beyond mid-range cutoff in minimum image &
 &convention at reduced resolution (point approximation). This applies to all residue pairs, for which &
 &at least one of the two is flagged as carrying a net charge, and happens at every neighbor list update.'
      end if
    else if (lrel_md.eq.5) then
      if ((use_cutoffs.EQV..true.).AND.(use_POLAR.EQV..true.)) then
        write(ilog,*)'Computing all monopole-monopole and monopole-dipole &
 &terms beyond mid-range cutoff explicitly in minimum image convention. This applies to all residue pairs, for which &
 &both are flagged as carrying a net charge, and happens at every neighbor list update.'
      end if
    else if ((lrel_md.eq.2).AND.(use_POLAR.EQV..true.)) then
      if (ewald_mode.eq.2) then
        write(ilog,*) 'Using standard Ewald summation for long-range electrostatics treatment.'
        write(ilog,70) ' Ewald parameter in 1/A   : ',ewpm
        write(ilog,72) ' Error estimate (a.u.)    : ',ewetol
        write(ilog,70) ' Ewald constant (kcal/mol): ',ewcnst
        write(ilog,*) 'Real-space cutoff is same as short-range cutoff.'
        write(ilog,24) ' Reciprocal space cutoffs : ',&
 &   kdims(1,2)-kdims(1,1)+1,' x ',&
 &   kdims(2,2)-kdims(2,1)+1,' x ',&
 &   kdims(3,2)-kdims(3,1)+1
      else if (ewald_mode.eq.1) then
        write(ilog,*) 'Using Particle-Mesh Ewald (PME) summation for long-range electrostatics treatment.'
        write(ilog,70) ' Ewald parameter in 1/A   : ',ewpm
        write(ilog,72) ' Error estimate (a.u.)    : ',ewetol
        write(ilog,72) ' Tolerance for tab. erfc(): ',ewerfctol
        write(ilog,70) ' Ewald constant (kcal/mol): ',ewcnst
        write(ilog,21) ' Order of B-splines       : ',splor
        call strlims(fftfname,ff1,ff2)
        if (fftlevel.eq.5) then
          write(ilog,*) 'FFTW/DFFT wisdom is read from file ',fftfname(ff1:ff2),'.'
        else
          write(ilog,21) ' FFTW/DFFT planning level : ',fftlevel
          if (fftprint.eq.1) then
            write(ilog,*) 'FFTW/DFFT wisdom is written to file ',fftfname(ff1:ff2),'.'
          end if
        end if
        write(ilog,*) 'Real-space cutoff is same as short-range cutoff.'
        write(ilog,24) ' Mesh dimensions          : ',&
 &   kdims(1,2),' x ',&
 &   kdims(2,2),' x ',&
 &   kdims(3,2)
      end if
    else if ((lrel_md.eq.3).AND.(use_POLAR.EQV..true.)) then
      if (rf_mode.eq.1) then
        write(ilog,*) 'Using generalized reaction field correction for electrostatic interactions.'
        write(ilog,70) ' Continuum dielectric     : ',par_IMPSOLV(2)
        write(ilog,70) ' Ionic strength in #/nm^3 : ',1000.0*par_POLAR(3)
        write(ilog,70) ' RF cutoff in A           : ',mcel_cutoff
        write(ilog,70) ' RF constant (kcal/mol)   : ',rfcnst
      else if (rf_mode.eq.2) then
        write(ilog,*) 'Using standard reaction field correction for electrostatic interactions.'
        write(ilog,70) ' Continuum dielectric     : ',par_IMPSOLV(2)
        write(ilog,70) ' RF cutoff in A           : ',mcel_cutoff
        write(ilog,70) ' RF constant (kcal/mol)   : ',rfcnst
      end if
    end if
    if (is_lj.EQV..true.) then
      write(ilog,*) 'Using dedicated LJ loops for speed reasons.'
    else if (is_ev.EQV..true.) then
      write(ilog,*) 'Using dedicated EV loops for speed reasons.'
    else if (((is_plj.EQV..true.).OR.(is_pewlj.EQV..true.).OR.&
 &  (is_prflj.EQV..true.)).AND.(use_waterloops.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
      write(ilog,*) 'Using dedicated gas-phase loops including optim&
 &ized water loops for speed reasons.'
    else if ((is_plj.EQV..true.).OR.(is_pewlj.EQV..true.).OR.&
 &           (is_prflj.EQV..true.)) then
      write(ilog,*) 'Using dedicated gas-phase loops for speed reaso&
 &ns.'
!    else if (is_fegev.EQV..true.) then
!      write(ilog,*) 'Using dedicated EV loops with support for ghost&
! &ed particles for speed reasons.'
!    else if (is_feglj.EQV..true.) then
!      write(ilog,*) 'Using dedicated LJ loops with support for ghost&
! &ed particles for speed reasons.'
    else if ((is_fegplj.EQV..true.).AND.&
 &                   (use_waterloops.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
      write(ilog,*) 'Using dedicated gas-phase loops with support fo&
 &r ghosted particles and optimized water loops for speed reasons.'
    else if ((is_fegprflj.EQV..true.).AND.&
 &                   (use_waterloops.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
      write(ilog,*) 'Using dedicated gas-phase loops with support fo&
 &r ghosted particles, the RF method and optimized water loops for s&
 &peed reasons.'
    else if (is_fegplj.EQV..true.) then
      write(ilog,*) 'Using dedicated gas-phase loops with support fo&
 &r ghosted particles for speed reasons.'
    else if ((is_implj.EQV..true.).OR.(is_impljp.EQV..true.).AND.(use_cutoffs.EQV..true.)) then
      write(ilog,*) 'Using dedicated ABSINTH-loops for speed reasons.'
    else if (is_fegprflj.EQV..true.) then
      write(ilog,*) 'Using dedicated gas-phase loops with support fo&
 &r ghosted particles and the RF method for speed reasons.'
    else if (is_tab.EQV..true.) then
      write(ilog,*) 'Using dedicated loops with support for just tabulated&
 & potentials for speed reasons.'
    else
      write(ilog,*) 'Using standard (general and slightly slower) loops.'
    end if
  end if
!
  write(ilog,*)
  if (ua_model.eq.1) then
    write(ilog,*) 'Polymers are built without aliphatic hydrogen atoms (GROMOS-compatible United-Atom Model).'
  else if (ua_model.eq.2) then
    write(ilog,*) 'Polymers are built without aliphatic or aromatic hydrogen atoms (CHARMM19-compatible United-Atom Model).'
  else
    write(ilog,*) 'Polymers are built with all hydrogen atoms (All-Atom Model).'
  end if
  write(ilog,*)
  if (use_14.EQV..true.) then
    write(ilog,*) '14-Interactions included in energy calculation.'
  else
    write(ilog,*) '14-Interactions excluded in energy calculation.'
  end if 
  write(ilog,*)
  if (nbsr_model.eq.2) then
    write(ilog,*) 'Standard MMFF NB model used (all atoms separated &
 &by at least three bonds interact using 14-fudge factors, irrespect&
 &ive of constraints on the rotation around these bonds).'
  else if (nbsr_model.eq.3) then
    write(ilog,*) 'GROMOS MMFF NB model used (most atoms separated &
 &by at least three bonds interact using 14-fudge factors) with special&
 & exclusion rules (most rigid 14- and some rigid 15-interactions within&
 & aromatic ring systems are excluded).'
  else if (nbsr_model.eq.1) then
    write(ilog,*) 'Modified NB model used (all atoms separated by ri&
 &gid or quasi-rigid bonds only are excluded from NB energy evaluati&
 &ons).'
    if (mode_14.eq.2) then
      write(ilog,*) 'Fudge factors for 14-interactions are assigned &
 &whenever exactly one relevant, rotatable bond separates two atoms.&
 &'
    else if (mode_14.eq.1) then
      write(ilog,*) 'Fudge factors for 14-interactions are assigned &
 &whenever exactly three bonds (the middle one rotatable) separate t&
 &wo atoms.'
    end if
  end if
  write(ilog,*)
  if ((elec_model.eq.1).AND.(use_POLAR.EQV..true.)) then
    write(ilog,*) 'Standard MMFF electrostatics model used (all poin&
 &t charges separated by at least one rotatable bond interact). Fudg&
 &ing of 14-interactions applies.'
  else if ((elec_model.eq.2).AND.(use_POLAR.EQV..true.)) then
    write(ilog,*) 'Reduced (sane) group-based electrostatic model us&
 &ed (only full charge groups interact (fully, fully fudged or not a&
 &t all) as determined by automatic charge group finder).'
  end if
  write(ilog,*)
  if (globcyclic.EQV..true.) then
    write(ilog,*) 'Polymer chain(s) requested to be cyclic but your request was blatantly ignored by the software.'
  else
    write(ilog,*) 'Polymer chains (if any) are linear and non-cyclic.'
  end if
  write(ilog,*)
  if ((pdb_analyze.EQV..true.).OR.(do_restart.EQV..true.)) then
!   write nothing
  else
    tailscnt = 0
    pdbrsize = 0
    sccnt = 0
    if (allocated(pdb_didread).EQV..true.) then
      pdbrsize = size(pdb_didread,dim=1)
      sccnt = sum(pdb_didread(:,3))
      atomscnt = 0
      if (allocated(pdbmap2).EQV..true.) then
        do i=1,atmol(pdbrsize,2)
          if (pdbmap2(i).le.0) atomscnt = atomscnt + 1
        end do
      end if
    end if
    do i=1,pdbrsize
      tailscnt = tailscnt + pdb_didread(i,1)-rsmol(i,1) + rsmol(i,2)-pdb_didread(i,2)
    end do
    if (pdbrsize.eq.nmol) then
      if ((tailscnt.eq.0).AND.(globrandomize.ne.2).AND.(globrandomize.ne.3)) then
        write(ilog,'(1x,a)') 'System geometry is read entirely from PDB-file (see above).'
      else if (tailscnt.eq.0) then
        write(ilog,'(1x,a)') 'System geometry is read from PDB-file (see above), but only the intramolecular arrangement is kept.'
      else if ((globrandomize.ne.2).AND.(globrandomize.ne.3)) then
        write(ilog,'(1x,a,i5,a)') 'With the exception of ',tailscnt,' residues in missing tails, system geometry is read entirely&
 & from PDB-file (see above).'
      else
        write(ilog,'(1x,a,i5,a)') 'With the exception of ',tailscnt,' residues in missing tails, system geometry is read entirely&
 & from PDB-file (see above), but only the intramolecular arrangement is kept.'
      end if
      if (atomscnt.gt.0) then
        write(ilog,'(1x,a,i7,a)') 'A total of ',atomscnt,' atoms was missing in these molecules, which were built in default &
 &covalent geometries.'
      end if
      if ((sccnt.gt.0).AND.((globrandomize.eq.1).OR.(globrandomize.eq.2))) then
        write(ilog,'(1x,a,i7,a)') 'Missing heavy atoms led to the randomization of ',sccnt,' side chains.'
      else if (sccnt.gt.0) then
        write(ilog,'(1x,a,i7,a)') 'Missing, nonterminal heavy atoms were identified in ',sccnt,' side chains.'
      end if
    else if (pdbrsize.eq.0) then
      if ((globrandomize.eq.1).OR.(globrandomize.eq.2)) then
        write(ilog,'(1x,a)') 'System was initially randomized (using natively supported degrees of freedom and OTHER torsions &
 &in unsupported residues).'
        if (n_crosslinks.gt.0) then
          write(ilog,'(1x,a)') 'Any crosslinks were attempted to be satisfied exactly.'
        end if
      else
        write(ilog,'(1x,a)') 'System was initially randomized (rigid-body degrees of freedom). Any polymer chains are &
 &built in (mostly extended) default conformations.'
        if (n_crosslinks.gt.0) then
          write(ilog,'(1x,a)') 'Only intermolecular (if any) crosslinks were attempted to be satisfied exactly, and &
 &intramolecular ones (if any) were ignored.'
        end if
      end if
    else if (pdbrsize.gt.0) then
      if (tailscnt.gt.0) then
        write(ilog,'(1x,a,i5,a,i6,a)') 'With the exception of ',tailscnt,' residues in missing tails, the system geometry &
 &for the first ',pdbrsize,' molecule(s) was read from PDB-file (see above).'
      else
        write(ilog,'(1x,a,i6,a)') 'The system geometry for the first ',pdbrsize,' molecule(s) was read from &
 &PDB-file (see above).'
      end if
      if ((globrandomize.eq.2).OR.(globrandomize.eq.3)) then
        write(ilog,'(1x,a)') 'Only the intramolecular arrangement was kept.'
      end if
      if (atomscnt.gt.0) then
        write(ilog,'(1x,a,i7,a)') 'A total of ',atomscnt,' atoms was missing in these molecules, which were built in default &
 &covalent geometries.'
      end if
      if ((sccnt.gt.0).AND.((globrandomize.eq.1).OR.(globrandomize.eq.2))) then
        write(ilog,'(1x,a,i7,a)') 'Missing heavy atoms led to the randomization of ',sccnt,' side chains.'
      else if (sccnt.gt.0) then
        write(ilog,'(1x,a,i7,a)') 'Missing, nonterminal heavy atoms were identified in ',sccnt,' side chains.'
      end if
      if ((globrandomize.eq.2).OR.(globrandomize.eq.1)) then
        write(ilog,'(1x,a)') 'Any missing molecules were randomized internally (if necessary and possible) and then placed &
 &randomly (inevitable).'
      else
        write(ilog,'(1x,a)') 'Any missing molecules were left in default conformations (extended for polymers) and then placed &
 &randomly (inevitable).'
      end if
    end if
    if (((tailscnt.gt.0).OR.(sccnt.gt.0)).AND.(globrandomize.ne.0).AND.(globrandomize.ne.3)) then
      write(ilog,'(1x,a)') 'Tails and side chain randomizations used a hierarchical Monte Carlo procedure based on excluded &
 &volume, bonded, and boundary interactions (if any).' 
      write(ilog,70) ' En. threshold (kcal/mol) per res. pair: ',globrdmthresh
      write(ilog,21) ' Max. attempts per res. in every tail  : ',globrdmatts
      write(ilog,'(1x,a)') 'See documentation on FMCSC_RANDOMIZE for further information.'
    else if (tailscnt.gt.0) then
      write(ilog,'(1x,a)') 'Tails were simply added in default polymer conformations (usually extended).'
    end if
    if (size(pdb_didread,dim=1).lt.nmol) then
      write(ilog,'(1x,a)') 'All randomizations of completely missing molecules used a hierarchical Monte Carlo procedure based &
 &on excluded volume, torsional, and boundary interactions (if any).'
      write(ilog,70) ' En. threshold (kcal/mol) per res. pair: ',globrdmthresh
      write(ilog,21) ' Max. att. per unit (mol./loop/res.)   : ',globrdmatts
      write(ilog,'(1x,a)') 'See documentation on FMCSC_RANDOMIZE for further information.'
    end if
  end if
  write(ilog,*)
  ti = 0
  ti2 = 0
  ti3 = 0
  do i=1,nseq
    ti = ti + nrsintra(i)
    if (i.lt.nseq) ti2 = ti2 + nrsnb(i)
    do j=i+2,nseq
      ti3 = ti3 + (at(i)%nsc+at(i)%nbb)*(at(j)%nsc+at(j)%nbb)
    end do
  end do
  write(ilog,*)
  if (pdb_analyze.EQV..false.) then
    write(ilog,*) 'SUMMARY OF DEGREES OF FREEDOM IN SYSTEM:'
    if (fycxyz.eq.1) then
      write(ilog,49) 'Torsional degrees of freedom    : ',ntorsn
      write(ilog,49) 'Pucker degrees of freedom       : ',5*(ntorpuck-ntorsn)/7
      write(ilog,49) 'Rigid-body degrees of freedom   : ',totrbd
    else if (fycxyz.eq.2) then
      write(ilog,49) 'Cartesian degrees of freedom    : ',3*n
    end if
    if (dyn_mode.ne.1) then
      write(ilog,49) 'Number of constrained d.o.f.s   : ',n_constraints
    end if
  end if
  write(ilog,*) 'SUMMARY OF ATOM-ATOM INTERACTIONS IN SYSTEM:'
  write(ilog,49) 'Intra-residue                   : ',ti
  write(ilog,49) 'Neighboring residues            : ',ti2
  write(ilog,49) 'Remaining residue pairs         : ',ti3
  write(ilog,49) 'Total number                    : ',ti+ti2+ti3
!  write(ilog,*) 'Initial energy (old fxn) : ',esaver
  if ((pdb_analyze.EQV..false.).AND.(do_restart.EQV..false.)) then
    if (do_n2loop.EQV..true.) then
      write(ilog,71) ' Initial energy (min. image, N^2): ',esave
    end if
    if (use_cutoffs.EQV..true.) then
      write(ilog,71) ' Initial energy (cutoffs, LREL)  : ',esavec
    end if
  end if
  if ((dyn_mode.gt.1).AND.(dyn_mode.le.3)) then
    if (pdb_analyze.EQV..false.) then
      write(ilog,49) 'Dyn. progress output interval   : ',nsancheck
    end if
  else if (dyn_mode.eq.6) then
    write(ilog,49) 'Min. progress output interval   : ',nsancheck
  else if (dyn_mode.eq.1) then
    write(ilog,49) 'Energy reset interval           : ',nsancheck
  else 
    write(ilog,49) 'Reporting/reset interval        : ',nsancheck
  end if
!
  write(ilog,*)
  write(ilog,*) '--- END OF SUMMARY ---'
  write(ilog,*)
  flush(ilog)
!
end
!
!-----------------------------------------------------------------------
