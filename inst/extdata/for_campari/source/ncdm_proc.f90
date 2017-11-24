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
! MAIN AUTHOR:   Marco Bacci                                               !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!---------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
!---------------------------------------------------------------------------
!
subroutine ncdm_initial !just set what can be set immediately
!
  use system, ONLY: nequil
  use clusters, ONLY: cdis_crit,cdofset,sumssz
  use pdb, ONLY: select_frames,pdb_fileformat
!  use ncdm
 
!
  implicit none
!
! This we can set from the very beginning
  nequil = 0
  pdb_fileformat = 5      !This has to set free in case we will want to support sequential access frames files
  select_frames = .false.
!
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)) then
    sumssz = 5
  else
    sumssz = 2
  end if
!
! degrees of freedom for DRMS
  if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
    if (allocated(cdofset).EQV..true.) deallocate(cdofset)
  end if
!
end subroutine ncdm_initial
!
!---------------------------------------------------------------------------
!
subroutine ncdm_sanity_checks_1
!
  use iounit
  use clusters 
  use mpistuff, ONLY: re_traceinfile
  use ncdm
!
  implicit none
!
  integer t1,t2                       !to determine string limits
!
  logical exists
!
! If nothing is requested, chrash very bad
  if ((ncdm_doas.EQV..false.).AND.(ncdm_donc.EQV..false.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. It is required to specify EITHER an input ASCII file (FMCSC_NCDM_ASFILE) &
 &or a NetCDF file (FMCSC_NCDM_NCFILE).'
    call fexit()
  end if
! Should not be possible
  if ((ncdm_isdmonas.EQV..true.).AND.(ncdm_doas.EQV..false.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got ncdm_doas false with ncdm_isdmonas TRUE. This is a bug.'
    call fexit()
  end if
!
! We do not even welcome you if you misbehave so bad 
  if ((ncdm_doas.EQV..true.).AND.(ncdm_donc.EQV..true.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Requesting concurrent processing of ascii-file and net-cdf file is &
 &explicitly disallowed for many reasons. If the goal is to convert and/or analyse the ASCII file, just set &
 &NCDM_WRTINPUT to 1 and/or NCDM_ANONAS to 1. &
 &If the goal is to mine the net-cdf input file, just comment or remove keyword NCDM_ASFILE.'
    call fexit()
  end if
  if ((ncdm_doas.EQV..true.).AND.(ncdm_isdmonas.EQV..false.).AND.(ncdm_wrtinp.EQV..false.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. It is meaningless to specify an input ascii file without requiring neither conversion &
 &nor data analysis.'
    call fexit()
  end if
!
  write(ilog,*)
  write(ilog,*) '-- CAMPARI IS USED IN NETCDF-BASED DATA MINING ONLY    --'
  write(ilog,*) '-- ALL SETTINGS THAT ARE NOT RELEVANT ARE IGNORED      --'
  write(ilog,*)
!
  if (ncdm_doas.EQV..true.) then
    write(ilog,*) '-- AN EXTERNAL DATA FILE WILL BE PROCESSED   --'    
    if (ncdm_wrtinp.EQV..true.) then
      write(ilog,*) '-- AND CONVERTED TO NET-CDF FORMAT           --'    
    end if
    if (ncdm_isdmonas.EQV..true.) then
      write(ilog,*) '-- AND ANALYZED WITH SPECIFIC ROUTINES       --'
    end if
    write(ilog,*)
    if (ncdm_checkas.EQV..false.) then
      write(ilog,*) 'Warning. No check on input ascii file will be performed (default). If there is not the &
 &exact number of features expected per line, wrong results may be generated without a crash.'
    else 
      write(ilog,*) 'Warning. Checking of input ascii file requested. This may be very slow.'
    end if
  end if
  if (ncdm_donc.EQV..true.) then
    write(ilog,*) '-- AN EXTERNAL NETCDF FILE WILL BE READ AND  --'    
    write(ilog,*) '-- POSSIBLY MINED WITH SPECIFIC ROUTINES     --'    
    write(ilog,*)
  end if
!
  if ((ncdm_isdmonas.EQV..true.).AND.(ncdm_doas.EQV..false.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Requested analysis on converted ascii without supplying an input ascii file. This is a bug.'
    call fexit()
  end if
  if ((ncdm_isdmonas.EQV..true.).AND.(ncdm_donc.EQV..true.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Requested analysis on both converted ascii file and input net-cdf. This is a bug.'
    call fexit()
  end if
!
  if (cstorecalc.le.0) then
    write(ilog,*)
    write(ilog,*) 'Fatal. CCOLLECT is set to zero or negative number.'
    call fexit()
  end if
  if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then
    if (cstorecalc.gt.1) then
      write(ilog,*) 'Warning. CCOLLECT is greater than 1. This will entail subsampling the input data to the &
   &specified frequency, regardless their input origin (ascii or net-cdf) for analysis (but not for the possible &
   &conversion of the input ASCII, if any).'
      if (ncdm_isthere_framesfl.EQV..true.) then
        write(ilog,*)
        write(ilog,*) 'Fatal. The concurrent support of CCOLLECT and NCDM_FRAMESFILE is explicitly not offered. &
   &Please, either set CCOLLECT to 1 or do not use a FRAMESFILE.' 
        call fexit()
      end if
    end if
    if (ncdm_isdmonas.EQV..true.) then
      if (3.ge.int(floor(1.*ncdm_nframes_ini/cstorecalc))) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Due to the settings of CCOLLECT and NCDM_NRFRMS there will be no more than three &
 &frames to analyze and no data mining can be performed on the input ASCII file despite the request.'
        call fexit()
      end if
      if ((ncdm_is_periodic.EQV..true.).AND.(cdis_crit.gt.4)) then 
        write(ilog,*)
        if (ncdm_wrtinp.EQV..true.) then
          write(ilog,*) 'Warning. The specified CDISTANCE does not support periodic data and periodicity will not be accounted &
 &for in analysis in any way. However, the converted ascii file will contain the periodicity attribute.' !not touch ncdm_is_period
        else 
          write(ilog,*) 'Warning. The specified CDISTANCE does not support periodic data and periodicity will not be accounted &
          &for in analysis.'
          ncdm_is_periodic = .false. !okay as no conversion is requested, we are not going to miss the attribute anyhow
        end if 
      else if ((ncdm_is_periodic.EQV..false.).AND.(cdis_crit.le.4)) then
        write(ilog,*) 'The specified CDISTANCE requires that the data are periodic. Given the lack of specification in the input &
 &key-file, a periodic range from -180 to 180 is assumed in analysis. Input data will be scaled and centered within this interval.'
        if (ncdm_wrtinp.EQV..true.) then
          write(ilog,*) 'However, no periodicity attribute/rescaling will be present in the converted Net-cdf file.'
        end if
      end if
    end if
  else !only ascii conversion requested
    if (ncdm_is_periodic.EQV..true.) then
      write(ilog,*) 'Periodicity attribute will be present in the converted Net-cdf file but data are not scaled and centered &
 &within this interval in the converted output file.'
    end if
    if (cstorecalc.gt.1) then
      write(ilog,*) 'Warning. CCOLLECT is greater than 1 but this will be ignored in the requested conversion of the &
 &input ASCII file.'
      if (ncdm_isthere_framesfl.EQV..true.) then
        write(ilog,*)
        write(ilog,*) 'Fatal. The concurrent support of CCOLLECT and NCDM_FRAMESFILE is explicitly not offered. &
   &Not even when just converting an input ASCII file. Please, either set CCOLLECT to 1 or do not use a FRAMESFILE.' 
        call fexit()
      end if
    end if
  end if
!
  if (ncdm_doas.EQV..true.) then 
    call strlims(ncdm_fn_as,t1,t2)
    inquire(file=ncdm_fn_as(t1:t2),exist=exists)  !input ascii file
    if (exists.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. ASCII file ',ncdm_fn_as(t1:t2),' not found.' 
      call fexit()
    end if
    if (ncdm_nframes_ini.le.0) then
      write(ilog,*)
      write(ilog,*) 'Fatal. When converting an input ascii file to net-cdf it is mandatory to specify the &
 &number of frames (NCDM_NRFRMS).'
      call fexit()
    end if
    if (ncdm_nfeats_ini.le.0) then
      write(ilog,*)
      write(ilog,*) 'Fatal. When converting an input ascii file to net-cdf it is mandatory to specify the &
 &number of features (NCDM_NRFEATS).'
      call fexit()
    end if
  end if !if (ncdm_doas.EQV..true.) then
!
  if (ncdm_donc.EQV..true.) then
    call strlims(ncdm_fn_r,t1,t2)
    inquire(file=ncdm_fn_r(t1:t2),exist=exists)  !input net-cdf file
    if (exists.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Net-CDF file ',ncdm_fn_r(t1:t2),' not found.' 
      call fexit()
    end if
  end if
!
! cfile 
  if (ncdm_isthere_cfile.EQV..true.) then
    call strlims(ncdm_nm_cfile,t1,t2)
    inquire(file=ncdm_nm_cfile(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal, provided input cfile (NCDM_CFILE) ',ncdm_nm_cfile(t1:t2),' not found or corrupted.'
      ncdm_isthere_cfile = .false.
      call fexit()
    end if
    if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then
      write(ilog,*) 'Warning. A CFILE has been provided. This will entail using only a subset of the features regardless &
 &the input origin of the data (ascii or net-cdf) for analysis (but not for the possible conversion of the input ASCII, if any).'
    else !conversion only
      write(ilog,*) 'Warning. The specified CFILE is ignored in the requested conversion of the input ASCII file.'
    end if
  end if
!
! frames file
  if (ncdm_isthere_framesfl.EQV..true.) then
    call strlims(ncdm_nm_framesfl,t1,t2)
    inquire(file=ncdm_nm_framesfl(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal, provided input frames file (NCDM_FRAMESFILE) ',ncdm_nm_framesfl(t1:t2),' not found or corrupted.'
      ncdm_isthere_framesfl = .false.
      call fexit()
    end if
    if (cstorecalc.ne.1) then 
      write(ilog,*)
      write(ilog,*) 'Fatal. The concurrent support of NCDM_FRAMESFILE and CCOLLECT is explicitly not offered. &
 &Plese, either set CCOLLECT to 1 or do not use a FRAMESFILE.' 
      call fexit()
    end if
    if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then
      write(ilog,*) 'Warning. A FRAMESFILE has been provided. This will entail using only a subset of the frames regardless &
 &the input origin of the data (ascii or net-cdf) for analysis (but not for the possible conversion of the input ASCII, if any).'
    else !conversion only
      write(ilog,*) 'Warning. The specified FRAMESFILE is ignored in the requested conversion of the input ASCII file.'
    end if
  end if
!
! the rest matters only if analysis is actually performed
  if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then
!
! clustering and cut profiles
    if ( ((cdis_crit.gt.4).AND.(cdis_crit.lt.7)).OR.(cdis_crit.gt.9) ) then 
      write(ilog,*)
      write(ilog,*) 'Fatal. Unsupported distance criterion for Net-CDF data mining (see CDISTANCE). Gotten: ',cdis_crit
      call fexit()
    end if
    if ((cmode.eq.6).OR.(cmode.eq.7)) then
      write(ilog,*) 'Fatal. CMODE options 6 and 7 are not available in Net-CDF data mining mode.'
      call fexit()
    end if 
    if ((cmode.eq.5).OR.((cmode.eq.4).AND.(cprogindex.eq.2))) then
      if ((c_nhier-c_multires.lt.2).AND.(c_nhier.ge.2)) then
        c_multires = max(0,c_nhier - 2)
        write(ilog,*) 'Warning. The number of requested levels for multi-resolution clustering is incompatible with &
   &tree height (request for FMCSC_BIRCHMULTI changed to ',c_multires,').'
      end if
      if ((cmode.eq.4).AND.(cprogrdepth.gt.c_nhier)) then
        cprogrdepth = c_nhier
        write(ilog,*) 'Warning. The number of requested levels for additional search for approximate progress index &
   &is incompatible with tree height (request for FMCSC_CPROGRDEPTH changed to ',cprogrdepth,').'
      end if
    end if
    if ((ccfepmode.gt.0).AND.(cequil.le.0)) then
      write(ilog,*) 'Warning. It can produce misleading results to compute cut-based pseudo free energy profiles (setting for &
   &CMSMCFEP) without equilibrating the underlying Markov state model first (choice for CREWEIGHT).'
    end if
    if (ncdm_isdmonas.EQV..true.) then
      if ((cdis_crit.eq.2.OR.cdis_crit.eq.4.OR.cdis_crit.eq.9)) then 
        write(ilog,*) 'Warning. The current setting of CDISTANCE (to 2, 4 or 9) requires the specification of dynamic weigths. &
   &Flat (uniformative) weights will be used if no specification is found in the input ascii file. They might become informative &
   &in analysis depending on keyword CMODWEIGHTS and related ones.'
      else if (cdis_crit.eq.8) then
        write(ilog,*) 'Warning. The current setting of CDISTANCE to 8 requires the specification of static weigths. &
   &Flat (uniformative) weights will be used if no specification is found in the input ascii file. They might become informative &
   &in analysis depending on keyword CMODWEIGHTS and related ones.' 
      end if
    end if
    if ((ccfepmode.ge.8).AND.(ccfepmode.le.10)) then 
#ifdef LINK_HSL
      if (dopfold.EQV..false.) then
        write(ilog,*) 'Warning. For committor probability-based cut profiles (FMCSC_CMSMCFEP), it is required to &
   &perform at least the computation of (+) committors (FMCSC_DOPFOLD). Disabling.'
        ccfepmode = 0
      end if
      if ((ccfepmode.ge.9).AND.(ccfepmode.le.10)) then 
        if (dopfoldminus.EQV..false.) then
          if (ccfepmode.eq.9) then
            write(ilog,*) 'Warning. For the specified committor probability-based cut profile (FMCSC_CMSMCFEP), it is required &
   &to perform the computation of (-) committors (FMCSC_DOPFOLD_MINUS). Disabling.'
            ccfepmode = 0
          else
            write(ilog,*) 'Warning. Only part of specified committor probability-based cut profiles (FMCSC_CMSMCFEP) will be &
   &computed because (-) committors (FMCSC_DOPFOLD_MINUS) have not been requested.'
            ccfepmode = 8
          end if
        end if
      end if
#else
        write(ilog,*) 'Fatal. Attempting to use an unsupported feature. For committor probability-based cut profiles, it &
   &is required to link the code to the HSL library (see installation instructions).'
        call fexit()
#endif
    end if  
!
! data preprocessing
    if ((cdis_crit.ne.1).AND.(cdis_crit.ne.2).AND.((cprepmode.eq.2).OR.(cprepmode.eq.5))) then
      write(ilog,*) 'Warning. The normalization by sample variance (setting for CPREPMODE) can lead to nonsensical &
   &results or program crashes if one or more dimensions have negligible variance.'
    end if
!
! weigths for clustering
    if ((cchangeweights.gt.0).AND.(cdis_crit.ne.2).AND.(cdis_crit.ne.4).AND.(cdis_crit.ne.8).AND.(cdis_crit.ne.9)) then
      write(ilog,*) 'Warning. Requesting a replacement of weights for structural clustering (setting for CMODWEIGHTS) &
   &is only possible if an eligible distance function is used (CDISTANCE has to be 2, 4, or 8-9). Disabled.'
      cchangeweights = 0
    end if
!
! trace file is disabled
    t1 = 1
    t2 = 0
    call strlims(re_traceinfile,t1,t2)
    if (t1.gt.t2) then
      exists = .false.
    else
      inquire(file=re_traceinfile(t1:t2),exist=exists)
    end if
    if (exists.EQV..true.) then
      write(ilog,*) 'Fatal. Found that the TRACEFILE specified in the key file exists. This is not allowed.'
      call fexit()
    end if
!
! Sanity check against synthetic trajs., ergodicity, spectral analysis, pfold
    if ((synmode.ne.0).AND.((cmode.eq.4).AND.(cprogindex.eq.1))) then
      write(ilog,*) 'Warning. Disabling the generation of synthetic trajectories because no clusters are generated & 
   &when the exact progress index method is in use (CPROGINDMODE).'
      synmode = 0
    end if
    if ((synmode.eq.1).AND.(inissnap_usrslct.eq.endssnap_usrslct)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Initial and target snapshots for the generation of synth. trajs. &
    &are the same.'
      call fexit()
    end if
    if ((synmode.gt.1).AND.(nsskip.le.0)) then
      write(ilog,*) 'Fatal. It is nonsensical to request the generation of synthetic trajectories while suppressing all &
   &output (FMCSC_SYNTRAJOUT is 0) unless FMCSC_SYNTRAJ_MD is 1.'
      call fexit()
    end if
    if (eigvalmd.ne.0) then
#ifdef LINK_HSL
      if (arnblks.eq.-1) then  !the user has not specified it and it is still as in initial
        arnblks = numeig + 2
      end if
      if (arnstps.eq.-1) then  !the user has not specified it and it is still as in initial
         arnstps = ceiling((8.*numeig)/arnblks)
      end if
      if (eigvalmd.ne.2) then
        write(ilog,*) 'Warning. The chosen mode for the selection of the eigenvalues in the blocked Arnoldi method with &
   &Chebychev acceleration of the starting vectors differs from the default choice.'
      end if
      if (arnmaxitr.eq.0) then
        write(ilog,*) 'Warning. No restart of the Arnoldi method with Chebychev acceleration of the starting vectors is allowed.'
      end if
#else
      write(ilog,*)
      write(ilog,*) 'Fatal. Attempting to use an unsupported feature. For spectral analysis of transition matrix it is required &
   &to link the code to the HSL library (see installation instruction).'
      call fexit()
#endif
    end if
    if (dopfold.EQV..true.) then
#ifdef LINK_HSL
      inquire(file=clufoldfile,exist=exists)
      if (exists.EQV..false.) then
        write(ilog,*) 'Warning. Cannot open spec.d input file (',trim(clufoldfile),') for CLUFOLDFILE. &
   &Defaulting to largest cluster.'
        clufoldfile = 'Largest cluster will be used.'
      end if
      inquire(file=cluunfoldfile,exist=exists)
      if (exists.EQV..false.) then
        write(ilog,*) 'Warning. Cannot open spec.d input file (',trim(cluunfoldfile),') for CLUUNFOLDFILE. &
   &Defaulting to last (smallest) cluster.'
        cluunfoldfile = 'Last (smallest) cluster will be used.'
      end if
      if ((dopfoldminus.EQV..true.).AND.((caddlkmd.ge.4).AND.(caddlkmd.le.6))) then
        write(ilog,*) 'Warning. It is redundant to compute (-) committor probabilities (FMCSC_DOPFOLD_MINUS) if &
   &detailed balance is imposed (simply 1.0 - (+) committors).'
      end if
#else
      write(ilog,*)
      write(ilog,*) 'Fatal. Attempting to use an unsupported feature. For committor probabilities it is required &
   &to link the code to the HSL library (see installation instruction).'
      call fexit()
#endif
    end if
!
  end if !if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then
!
end subroutine ncdm_sanity_checks_1
!
!---------------------------------------------------------------------------
!
subroutine ncdm_summary_1
!
  use iounit
  use system, ONLY: basename
  use clusters
  use ncdm
  use threads, ONLY: thrdat !With no threads we need thrdat%rate (set up in initial), otherwise also thrdat%maxn.
!
  implicit none
!
  integer t1,t2
!
  logical exists
!
  5   format(a,1x,a)
  21  format(a,2x,i8)
  49  format(1x,a,i10)
  50  format(a,7x,f8.4)
  54  format(1x,a,i6)
  56  format(1x,a,i8)
  71  format(a,1x,g14.7)
  74  format(a,3x,g13.3)
!
  write(ilog,*)
  write(ilog,*) '-- SETTINGS RELEVANT TO NET-CDF DATA MINING (1st PART) --'
  if (ncdm_doas.EQV..true.) then
    write(ilog,*) 'Input ascii file (NCDM_ASFILE)         : ',trim(ncdm_fn_as)
    if (ncdm_wrtinp.EQV..true.) then
      write(ilog,*) 'Output net-cdf file will be written (NCDM_WRTINPUT).'
    end if
    if (ncdm_isdmonas.EQV..true.) then
      write(ilog,*) 'All analyses will refer to the above-specified ascii file (NCDM_ANONAS).'
    end if 
    if (ncdm_is_periodic.EQV..true.) then !user has specified keyword NCDM_PRDCRNG and only conv requested or cdist supports it
      if (ncdm_isdmonas.EQV..true.) then
        if (ncdm_wrtinp.EQV..true.) then 
          write(ilog,*) 'Periodicity range (NCDM_PRDCRNG) for conversion and analysis:',ncdm_periodic_rng(1),ncdm_periodic_rng(2)
        else 
          write(ilog,*) 'Periodicity range (NCDM_PRDCRNG) for analysis:',ncdm_periodic_rng(1),ncdm_periodic_rng(2)
          if (cdis_crit.gt.4) then
            write(ilog,*)
            write(ilog,*) 'Fatal. Gotten ncdm_is_periodic.EQV..true. with wrong cdistance and with no conversion. This is a bug.'
            call fexit()
          end if
        end if
      else if (ncdm_wrtinp.EQV..true.) then
        write(ilog,*) 'Periodicity range (NCDM_PRDCRNG) for attribute in conversion:',ncdm_periodic_rng(1),ncdm_periodic_rng(2)
      else 
        write(ilog,*)
        write(ilog,*) 'Fatal. Wrong combination of flags in summary. This is a bug. (1)'
        call fexit()    
      end if
    end if
    write(ilog,*) 'Number of features (NCDM_NRFEATS)        : ',ncdm_nfeats_ini
    write(ilog,*) 'Name of netcdf variable (NCDM_NMFV)      : ',trim(ncdm_nm_fvar)
    write(ilog,*) 'Number of frames (NCDM_NRFRAMES)         :',ncdm_nframes_ini
    if (ncdm_wrtinp.EQV..true.) then
      write(ilog,*) 'For conversion, all frames and features will be retained (with no adjustments/rescaling).'
    end if
  end if
!
! fast return for conversion only
  if ((ncdm_isdmonas.EQV..false.).AND.(ncdm_donc.EQV..false.)) then
    write(ilog,*) 'No analysis requested, only conversion of input ascii.'
    write(ilog,*) '------------------------------------------------------'
    write(ilog,*)
    if ((ncdm_doas.EQV..false.).AND.(ncdm_wrtinp.EQV..false.)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Wrong combination of flags in summary. This is a bug. (2)'
      call fexit()    
    end if
    return
  end if
!
  if (ncdm_donc.EQV..true.) then
    write(ilog,*) 'Input netcdf file (NCDM_NCFILE)        : ',trim(ncdm_fn_r)
    write(ilog,*) 'Assumed slowest dimension in input net-cdf file: frames'
    if (ncdm_wrtinp.EQV..true.) then
      write(ilog,*) 'Output net-cdf file will be written (NCDM_WRTINPUT). It will reflect many analysis settings.'
    end if
  end if
!
! framesfile, this matters regardless
  if (ncdm_isthere_framesfl.EQV..true.) then
    call strlims(ncdm_nm_framesfl,t1,t2)
    write(ilog,*) 'Input frames file (NCDM_FRAMESFILE) : ',ncdm_nm_framesfl(t1:t2)
    write(ilog,*) 'Frames will be subsampled accordingly.' 
  else 
    write(ilog,*) 'No input frames file provided (NCDM_FRAMESFILE).' 
  end if
!
! cfile, this matters regardless
  if (ncdm_isthere_cfile.EQV..true.) then
    call strlims(ncdm_nm_cfile,t1,t2)
    write(ilog,*) 'Input cfile (NCDM_CFILE) : ',ncdm_nm_cfile(t1:t2)
    write(ilog,*) 'Features will be subsampled accordingly.' 
  else 
    write(ilog,*) 'No input cfile provided (NCDM_CFILE). All features will be used.'
  end if
!
! data mining settings relevant regardless
  if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then
    call strlims(basename,t1,t2)
    write(ilog,5) ' Basename for output   (BASENAME) : ',basename(t1:t2)
    write(ilog,*) 'Distance criterion index  (CDISTANCE) : ',cdis_crit
!
    if (cdis_crit.le.4) then 
      write(ilog,*) 'The specified distance criterion implies that the data are periodic.'
      if (cdis_crit.ge.3) then
        write(ilog,*) 'The specified distance criterion implies that sine and cosine functions are meaningful &
 &to be computed for the input data.'
      end if
      if (ncdm_isdmonas.EQV..true.) then 
        if (ncdm_is_periodic.EQV..false.) then !the user has not specified the periodic range in the key-file but cdis requires it
          write(ilog,*) 'No further information have been specified with keyword NCDM_PRDCRNG. &
 &The assumed periodicity range is -180 to 180. Data exceeding these values will be rescaled within &
 &the reference interval in analysis.' 
        end if
      else
        write(ilog,*) 'Periodicity will be inferred directly from input net-cdf file from attribute periodic_range.'
        write(ilog,*) 'In case the attribute periodic_range is missing, the assumed range is -180 to 180. Data outside &
 &this range will be rescaled within the reference interval in analysis.' 
      end if
    end if
    write(ilog,*) 'Struct. clust. collection interval (CCOLLECT) : ',cstorecalc
!
    if (cstorecalc.eq.1) then
      write(ilog,*) 'All frames will be retained for analysis (CCOLLECT).'
    else if (cstorecalc.gt.1) then
        write(ilog,*) 'Due to the settings for CCOLLECT input frames will be subsampled at an interval of:',cstorecalc, &
 &' for analysis.'
    end if
!  
!   preprocessing
    if (cprepmode.eq.0) then
        write(ilog,*) 'No signal preprocessing requested for structural clustering (CPREPMODE).'
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
!
!   clustering 
    if ( (((cdis_crit.eq.2).OR.(cdis_crit.eq.4)).OR.((cdis_crit.ge.8).AND.(cdis_crit.le.9))).AND.(cchangeweights.gt.0)) then
      write(ilog,*) 'Mode for altered weights  : ',cchangeweights
      write(ilog,*) 'Combin. rule for weights  : ',cwcombination
      if ((cchangeweights.ne.2).AND.(cdis_crit.ne.8)) write(ilog,*) 'Window size (snapshots)   : ',cwwindowsz
      if (((cchangeweights.ge.4).AND.(cchangeweights.le.9)).AND.(cdis_crit.ne.8)) then
        write(ilog,*) ' Buffer parameter for inv. : ',cdynwbuf
      end if
      if ((cchangeweights.eq.2).OR.(cchangeweights.eq.3).OR.(cchangeweights.eq.6).OR.(cchangeweights.eq.7)) then
        write(ilog,*) 'Lag time for ACF (snaps.) : ',cwacftau
      end if
      if ((cprepmode.lt.3).AND.((cchangeweights.eq.5).OR.(cchangeweights.eq.7).OR.(cchangeweights.eq.9))) then
        write(ilog,*) 'Smoothing window (snaps.) : ',csmoothie
      end if
    end if
    if (cmode.eq.1) then
      write(ilog,*) 'Using simple leader-based clustering algorithm (does not require neighbor list) (CMODE).'
    else if (cmode.eq.2) then
      write(ilog,*) 'Using modified leader-based clustering algorithm (does not require neighbor list) (CMODE).'
    else if (cmode.eq.3) then
      write(ilog,*) 'Using hierarchical clustering algorithm (requires neighbor list) (CMODE).'
    else if ((cmode.eq.4).AND.(cprogindex.eq.1)) then
      write(ilog,*) 'Using cut-based one-shot clustering algorithm with exact progress sequence (requires neighbor list) (CMODE).'
    else if ((cmode.eq.4).AND.(cprogindex.eq.2)) then
      write(ilog,*) 'Using cut-based one-shot clustering algorithm with approximate progress sequence (does not require neighbor &
   &list) (CMODE).'
    else if (cmode.eq.5) then
      write(ilog,*) 'Using BIRCH-like clustering algorithm (does not require neighbor list) (CMODE).'
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
    if (cprogindex.ne.1) then
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
    end if
!
!   pca
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
        write(ilog,*) 'No. of components to use  : ',reduced_dim_clustering
      end if
    else
      if (reduced_dim_clustering.gt.0) then
        write(ilog,*) 'Data dimensions will be discarded at the end per user request.'
        write(ilog,*) 'No. of dimensions to use  : ',reduced_dim_clustering        
      end if
    end if
!
!   equilibration / iterations / profiles / networks
    call strlims(tbrkfilen,t1,t2)
    inquire(file=tbrkfilen(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &tbrkfilen(t1:t2),') for FMCSC_TRAJBREAKSFILE. All transitions will be kept (this may include replica-exchange swaps,&
 & transitions caused by trajectory concatenation, etc).'
    else 
      write(ilog,*) 'File with traject. breaks : ',tbrkfilen(t1:t2)
    end if
    call strlims(tlnkfilen,t1,t2)
    inquire(file=tlnkfilen(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*) 'WARNING: Cannot open spec.d input file (',&
 &tlnkfilen(t1:t2),') for FMCSC_TRAJLINKSFILE. No transitions will be added.'
    else
      write(ilog,*) 'File with add. traj. links: ',tlnkfilen(t1:t2)
    end if
    if (brklnk_report.gt.0) then
      write(ilog,*) 'In case a break and/or link file(s) have been provided, additional information on connections &
   &manipulation will be reported to standard output.'
    end if
    if (cequil.gt.0) then
      write(ilog,*) 'This implies independent treatment of strongly connected components.'
      write(ilog,50) ' Buffer for snapshot wts.  : ',cequilbuf
    else
      write(ilog,*) 'Complete graph is analyzed independently of connectedness. This may cause errors.'
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
    if (caddlkmd.gt.0) then
      write(ilog,*) 'Mode for adding links     : ',caddlkmd
      write(ilog,50) ' Thresh. for added links   : ',caddlinks
      write(ilog,50) ' Weight of added links     : ',caddlinkwt
    else
      write(ilog,*) 'No links (edges) are added to network (graph).'
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
    if ((tmat_report.EQV..true.).AND.((eigvalmd.ne.0).OR.(dopfold.EQV..true.).OR.(synmode.ne.0))) then
      write(ilog,*) 'Report file for transition matrix(ces) will be written.'
    end if
    if (clagt_msm.eq.1) then
      write(ilog,*) 'Snapshots in memory are assumed to be connected in sequence by default.'
    else
      write(ilog,*) 'Snapshots in memory are assumed to be connected with a time lag by default.'
      write(ilog,54) 'Lag time in stored snapshots    : ',clagt_msm
    end if
!
!   synthetic trajectories
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
!
#ifdef ENABLE_THREADS
    write(ilog,*) 'This is a shared memory (OpenMP-based) calculation without MPI layer.'
    write(ilog,21) ' Maximum number of threads to use      : ',thrdat%maxn
#else
    write(ilog,*) 'This is a single-CPU calculation (no parallelism of any kind).'
#endif
    write(ilog,*) '------------------------------------------------------'
    write(ilog,*)
!
  end if !if data mining is on, not just a conversion
!
end subroutine ncdm_summary_1
!
!---------------------------------------------------------------------------
!
subroutine ncdm_makeio(mode)
!
  use iounit
#ifdef ENABLE_THREADS
  use threads, ONLY: ithread
#endif
!
  implicit none 
!
  integer freeunit,mode
!
  logical exists
!
#ifdef ENABLE_THREADS
  character(60) threadsfile
  logical threadsopen
#endif
!
#ifdef ENABLE_THREADS
  data threadsopen/.false./
#endif
!
#ifdef ENABLE_THREADS
  save threadsopen
#endif
!
  if (mode.eq.1) then
#ifdef ENABLE_THREADS
    threadsfile = 'THREADS.log'
#endif
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
! close files and free up filehandles
  else if (mode.eq.2) then
#ifdef ENABLE_THREADS
    if (threadsopen) close(unit=ithread)
#endif
  else
    write(ilog,*)
    write(ilog,*) 'Fatal. Calling ncdm_makeio with unknown mode. This is a bug.'
    call fexit()
  end if ! if (mode.eq.1)
!   
end subroutine ncdm_makeio
!
!---------------------------------------------------------------------------
!
subroutine ncdm_sanity_checks_2 !refers to what is left for analysis
!
  use iounit
  use clusters, ONLY: cdis_crit,pcamode,clstsz,reduced_dim_clustering,cstorecalc
  use ncdm
!
  implicit none
!
  if (3.ge.int(floor(1.*ncdm_nframes/cstorecalc))) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Due to the settings for CCOLLECT, there are no more than 3 snapshots left for analysis.'
      call fexit()
  end if
  if (ncdm_nframes.le.0) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. The final number of frames is <= 0. This might be due to an input error or to a bug.'
    call fexit()
  end if
  if (ncdm_nfeats.le.0) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. The final number of features is <= 0. This might be due to an input error or to a bug.'
    call fexit()
  end if
  if (3.ge.ncdm_nframes) then
    write(ilog,*)
    write(ilog,*) 'Fatal. There are no more than three frames to analyze and no data mining can be performed despite the request.'
    call fexit()
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(ncdm_arethere_dywghts.EQV..true.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Both static and dynamic weights in use. This is a bug.'
    call fexit()
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(cdis_crit.ne.8)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Static weights in use and cdistance not 8. This is a bug.'
    call fexit()
  end if
  if ((ncdm_arethere_dywghts.EQV..true.).AND.((cdis_crit.ne.2).AND.(cdis_crit.ne.4).AND.(cdis_crit.ne.9))) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Static weights in use and cdistance does not support them. This is a bug.'
    call fexit()
  end if
  if ((ncdm_is_periodic.EQV..true.).AND.(cdis_crit.gt.4)) then !possible in case ascii conv is also requested but cdis > 4 
    write(ilog,*)
    write(ilog,*) 'Warning. Periodic indicator is set to .true. but the selected distance criterion (CDISTANCE) &
 &does not account for periodicity. Disabling periodicity correction.'
    ncdm_is_periodic = .false.
  end if
!
! Sanity check against bad clustering requests
  if (reduced_dim_clustering.gt.ncdm_nfeats) then
    write(ilog,*) 'Warning. Reduced dimensionality clustering implies &
 & that the chosen number of dimensions ( 7-9: distances; &
 &1-2: dihedral angles; 3-4: Fourier terms) is smaller than the requested one (adjust&
 & clustering subset and/or setting for FMCSC_CREDUCEDIM). Turning off dimensionality reduction.'
    reduced_dim_clustering = 0
  else if ((mod(reduced_dim_clustering,2).ne.0).AND.(cdis_crit.eq.4)) then
    write(ilog,*) 'Warning. Uneven numbers for FMCSC_CREDUCEDIM are not supported if FMCSC_CDISTANCE is 4. Adjusted.'
    if (reduced_dim_clustering.lt.clstsz) then
      reduced_dim_clustering = reduced_dim_clustering + 1
    else
      reduced_dim_clustering = reduced_dim_clustering - 1 ! should never trigger
    end if
  end if
!
#ifdef LINK_LAPACK
! Sanity check against PCA
  if (pcamode.gt.1) then
!
    if (cdis_crit.le.2) then
      write(ilog,*) 'Fatal. Principal components are currently not supported &
 &for CDISTANCE 1 or 2 (supposedly periodic data). Please, set PCAMODE to 1.'
      call fexit()
    end if
!
    if (((pcamode.eq.3).OR.(pcamode.eq.5)).AND.(reduced_dim_clustering.gt.0)) then
      if (reduced_dim_clustering.gt.clstsz) then
        write(ilog,*) 'Fatal: Reduced dimensionality clustering after PCA implies &
 &that the chosen number of dimensions is smaller than the requested one (adjust &
 &clustering subset and/or setting for CREDUCEDIM).'
        call fexit()
      end if
    else if (reduced_dim_clustering.gt.0) then
      write(ilog,*) 'Fatal. Reduced dimensionality clustering is requires option 3 or 5 for PCAMODE.'
      call fexit()
    end if 
!
    if (cdis_crit.eq.9) then
      write(ilog,*) 'Warning. For PCA, locally adaptive weights are converted to static ones, which &
 &are then used to scale the data.'
      if (reduced_dim_clustering.gt.0) then
        write(ilog,*) 'This may lead to undesired effects for subsequent processing of data in reduced dimensionality.'
      end if
    end if
  end if
#endif
!
end subroutine ncdm_sanity_checks_2
!
!---------------------------------------------------------------------------
!
subroutine ncdm_summary_2 !refers to what is left for analysis
!
  use iounit
  use clusters, ONLY: cdis_crit,inissnap_usrslct,endssnap_usrslct,synmode
  use ncdm
  use pdb, ONLY: framelst,pdb_fileformat
!
  implicit none
!
  write(ilog,*)
  write(ilog,*) '-- SETTINGS RELEVANT TO NET-CDF DATA MINING (2nd PART) --'
  write(ilog,*) 'Number of frames present in the input file     : ',ncdm_nframes_ini
  write(ilog,*) 'Number of features present in the input file   : ',ncdm_nfeats_ini
  write(ilog,*) 'Number of features actually used               : ',ncdm_nfeats
  write(ilog,*) 'Final number of frames actually used           : ',ncdm_nframes
  if (cdis_crit.le.4) then
      write(ilog,'(A50,2x,F12.6,2x,F12.6)') ' Assumed periodicity range                      : ',&
 &ncdm_periodic_rng(1),ncdm_periodic_rng(2)
  else if (ncdm_is_periodic.EQV..true.) then
      write(ilog,*) 'Warning. If data are periodic, CDISTANCE should be set to either 1, 2, 3 or 4.'
      write(ilog,*) 'With current value of CDISTANCE, no periodicity correction will be performed.'
      write(ilog,*) 
      write(ilog,*) 'Fatal. Bad set of flags in summary_2. This is a bug.' !fixed in sanity_checks_2
      call fexit()
  else
    write(ilog,*) 'No periodicity correction in use.'
  end if
  if (ncdm_arethere_stwghts.EQV..true.) then
    write(ilog,*) 'Static weights are in use.'
  end if
  if (ncdm_arethere_dywghts.EQV..true.) then
    write(ilog,*) 'Dynamic weights are in use.'
  end if
  !if ((pdb_fileformat.eq.5).AND.(ncdm_isthere_framesfl.EQV..true.)) then  !aka random access frames file
    if (inissnap_usrslct.gt.0) then
     write(ilog,*) 'INISYNSNAP is absolute snapshot : ',framelst(inissnap_usrslct)
    end if 
    if ((synmode.eq.1).AND.(endssnap_usrslct.ne.0)) then
      write(ilog,*) 'ENDSYNSNAP is absolute snapshot : ',framelst(endssnap_usrslct)
    end if
  !end if
  write(ilog,*) '---------------------------------------------------------'
  write(ilog,*)
!
end subroutine ncdm_summary_2
!
!---------------------------------------------------------------------------
!
subroutine ncdm_read_framesfl
!
  use iounit
  use ncdm
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none 
!
  integer iomess                     !error management
  integer iu,freeunit,t1,t2,nframes_framesfl
  integer dummyi,iframe,maxframe,counter
!
  character(MAXSTRLEN) str2
!
  logical arethere_wrongs            !to manage nulls and negatives. Duplicates are fine in framesfile
!
  nframes_framesfl = 0  !number of frames that are resulting from frames file
  maxframe = 0
  arethere_wrongs = .false.
  dummyi = 0 !never know
!
! find out actual number of frames and relevant maximum index (see ncdm_whichframes)
  iu = freeunit()
  call strlims(ncdm_nm_framesfl,t1,t2)
  open(unit=iu,file=ncdm_nm_framesfl(t1:t2),status='old')  !to count the amount of eligible frames (set in nframes_framesfl)
  do while(.true.)
    read(iu,*,iostat=iomess) dummyi !the other numbers possibly present on the line are discarded, no snapshot weights are read
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.ne.0) then 
      write(ilog,*) 'Fatal. File I/O error while reading frames file (1). Not an integer number? Scientific notation? Or a bug.'
      call fexit() 
    end if
    if (dummyi.gt.0) then
      maxframe = max(maxframe,dummyi)            !highest frame index that is in the frames file
      nframes_framesfl = nframes_framesfl + 1  !counting how many eligible frames are in the frames file
    else
      arethere_wrongs = .true.
    end if
  end do 
  close(unit=iu)
!
! checks 
  if (arethere_wrongs.EQV..true.) then
    write(ilog,*) 'Warning. It appears that there were negative and/or nulls in the specified frames file. &
 &Those will be ignored.'
  end if
  if (maxframe.eq.0) then
    write(ilog,*)
    write(ilog,*) 'Fatal. No eligible frame indicator in NCDM_FRAMESFILE. All nulls or negatives?'
    call fexit()
  end if
  if (nframes_framesfl.lt.1) then
    write(ilog,*)
    write(ilog,*) 'Fatal. No eligible frame indicator in NCDM_FRAMESFILE. Empty file?'
    call fexit()
  end if
!
! allocate and initialize the mask
  if (allocated(ncdm_whichframes).EQV..true.) deallocate(ncdm_whichframes) 
  allocate(ncdm_whichframes(nframes_framesfl)) !as many as there are eligible rows in the framesfile
  ncdm_whichframes(:) = 0
!
! read and assign eligible ones
  counter = 0
  iu = freeunit()
  call strlims(ncdm_nm_framesfl,t1,t2)
  open(unit=iu,file=ncdm_nm_framesfl(t1:t2),status='old')  !open frames file to store only eligible snapshot 
  do while(.true.)
    read(iu,*,iostat=iomess) dummyi
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.ne.0) then 
      write(ilog,*) 'Fatal. File I/O error while reading frames file (2). Not an integer number? Scientific notation?'
      call fexit() 
    end if
    if (dummyi.gt.0) then
      counter = counter + 1
      ncdm_whichframes(counter) = dummyi !store the sequence of eligible snapshot
    end if
  end do 
  close(unit=iu)
!
!checks 
  if (maxval(ncdm_whichframes).ne.maxframe) then !important to avoid possible bug on mask allocation 
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten maxval(ncdm_whichframes).ne.maxframe at the end of ncdm_read_framesfl. This is a bug.'
    call fexit()
  end if
  if (ncdm_doas.EQV..true.) then
    if (maxframe.gt.ncdm_nframes_ini) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Detected a frame indicator greater than the assumed number of frames (NCDM_NFRAMES) &
 &in input ascii file (NCDM_ASFILE).'
      call fexit()
    end if
  end if
!
! the mask
  if (ncdm_doas.EQV..true.) then            !this is important
    if (allocated(ncdm_framesmask).EQV..true.) deallocate(ncdm_framesmask)
    allocate(ncdm_framesmask(ncdm_nframes_ini)) !ncdm_nframes_ini is the number of rows of input ASCII
    ncdm_framesmask = .false.                   !adjusted below on eligible frames
    ncdm_framesmask(ncdm_whichframes) = .true.
  else                                          !this is useless, I do it so that I do not have to distinguish the two cases below
    if (allocated(ncdm_framesmask).EQV..true.) deallocate(ncdm_framesmask)
    allocate(ncdm_framesmask(maxframe))
    ncdm_framesmask = .false.  
    ncdm_framesmask(ncdm_whichframes) = .true.
  end if
!
! final checks 
  if (counter.ne.nframes_framesfl) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten a different number of eligible frames in the second reading. This is a bug.'
    call fexit()
  end if
  do iframe=1,nframes_framesfl
    if (ncdm_whichframes(iframe).le.0) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Gotten negative or null frame at the end of cleaning procedure. This is a bug.'
      call fexit()
    end if
  end do 
  if (nframes_framesfl.ne.size(ncdm_whichframes)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten nframes_framesfl.ne.size(ncdm_whichframes) at the end of ncdm_read_framesfl. This is a bug.'
    call fexit()
  end if
!
! set already the right value to ncdm_nframes
  ncdm_nframes = nframes_framesfl
!
end subroutine ncdm_read_framesfl
!
!---------------------------------------------------------------------------
!
subroutine ncdm_manage_whichframes(from_read_ncfl)
!
  use iounit
  use ncdm
!
  implicit none
!
  integer iframe,jframe
!
  logical, INTENT(IN)  :: from_read_ncfl
!
  if ((ncdm_doas.EQV..true.).AND.((from_read_ncfl.EQV..false.))) then !we cannot be here if we had a framesfile (see datasaw)
    if (ncdm_isthere_framesfl.EQV..true.) then
      write(ilog,*)   
      write(ilog,*) 'Fatal. Inside ncdm_manage_whichframes with the wrong flags combination (1). This is a bug.'
      call fexit()
    end if
    ncdm_nframes = ncdm_nframes_ini  !ccollect could still change it
    if (allocated(ncdm_whichframes).EQV..true.) deallocate(ncdm_whichframes)
    allocate(ncdm_whichframes(ncdm_nframes_ini))
    do iframe=1,ncdm_nframes_ini
      ncdm_whichframes(iframe) = iframe         !just use all of them
    end do
    if (allocated(ncdm_framesmask).EQV..true.) deallocate(ncdm_framesmask)
    allocate(ncdm_framesmask(ncdm_nframes_ini))
    ncdm_framesmask = .false.
    ncdm_framesmask(ncdm_whichframes) = .true.  !set all of them to eligible
  else if ((ncdm_donc.EQV..true.).AND.(from_read_ncfl.EQV..true.)) then !it already set many things (e.g. ncdm_nframes_ini)
    if (ncdm_isthere_framesfl.EQV..true.) then  !it already set many things in ncdm_read_framesfl (e.g. ncdm_nframes)
      if (allocated(ncdm_whichframes).EQV..false.) then
        write(ilog,*)   
        write(ilog,*) 'Fatal. In ncdm_manage_whichframes with framesfile and with ncdm_whichframes not allocated. This is a bug.'
        call fexit()
      end if 
      if (ncdm_nframes_ini.lt.maxval(ncdm_whichframes)) then  !first meaningful thing to do here is to check consistency
        write(ilog,*)
        write(ilog,*) 'Fatal. There are less snapshots in the input net-cdf file that what stated in the framesfile. (1)'
        call fexit()
      end if
      if (allocated(ncdm_framesmask)) deallocate(ncdm_framesmask)  !we have to do this as the previous assignment was incomplete
      allocate(ncdm_framesmask(ncdm_nframes_ini))  !now allocated to ncdm_nframes_ini rather than to maxframe
      ncdm_framesmask = .false.  
      ncdm_framesmask(ncdm_whichframes) = .true.   !set to true the ones corresponding to the eligible frames
    else  !so we are here for the first time with no framesfile but from the call in read_ncfl
      ncdm_nframes = ncdm_nframes_ini  
      if (allocated(ncdm_whichframes).EQV..true.) deallocate(ncdm_whichframes)  !no framesfile, so I can deallocate with no risk
      allocate(ncdm_whichframes(ncdm_nframes_ini))
      do iframe=1,ncdm_nframes_ini
        ncdm_whichframes(iframe) = iframe         !just use all of them
      end do
      if (allocated(ncdm_framesmask).EQV..true.) deallocate(ncdm_framesmask)
      allocate(ncdm_framesmask(ncdm_nframes_ini))
      ncdm_framesmask = .false.
      ncdm_framesmask(ncdm_whichframes) = .true.  !set all of them to eligible
    end if
  else if (ncdm_donc.EQV..false.) then
    write(ilog,*)   
    write(ilog,*) 'Fatal. Inside ncdm_manage_whichframes with the wrong flags combination (2). This is a bug.'
    call fexit()
  else     !if the call is from datasaw (which can happen only without framesfile) and ncdm_donc is true,
    return !it does nothing as there is yet no information on the number of frames in the net-cdf file
  end if
!
end subroutine ncdm_manage_whichframes
!
!---------------------------------------------------------------------------
!
subroutine ncdm_read_cfile 
!
  use iounit
  use interfaces
  use ncdm
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none 
!
  integer iomess                     !error management
  integer iu,freeunit,t1,t2,nfeats_cfile,counter,startfeat
  integer dummyi,ifeat,maxfeat,jfeat
!
  character(MAXSTRLEN) str2
!
  logical arethere_wrongs,atrue,notpositive
!
  integer, ALLOCATABLE :: tmparr(:)
!
  atrue = .true.
  nfeats_cfile = 0  !number of features that are resulting from cfile
  maxfeat = 0
  dummyi = 0        !never know
!
! find out actual number of features and relevant maximum index (see ncdm_whichfeats)
  iu = freeunit()
  call strlims(ncdm_nm_cfile,t1,t2)
  open(unit=iu,file=ncdm_nm_cfile(t1:t2),status='old')  !open cfile to get to know how many features will be used (if all ok)
  do while(.true.)
    read(iu,*,iostat=iomess) dummyi !all rest of colums excluded
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.ne.0) then 
      write(ilog,*) 'Fatal. File I/O error while reading cfile file (3). Not an integer number? Scientific notation? Or a bug.'
      call fexit() 
    end if
      maxfeat = max(maxfeat,dummyi)    !highest feature index that is in the cfile
      nfeats_cfile = nfeats_cfile + 1  !counting how many features are in the cfile
  end do 
  close(unit=iu)
!
! checks and possible fast return
  if (maxfeat.eq.0) then
    write(ilog,*)
    write(ilog,*) 'Fatal. No eligible feature indicator in NCDM_CFILE. All nulls or negatives?'
    call fexit()
  end if
  if (nfeats_cfile.lt.1) then
    write(ilog,*)
    write(ilog,*) 'Fatal. No eligible feature indicator in NCDM_CFILE. Empty file?'
    call fexit()
  end if
  if (nfeats_cfile.eq.1) then !only one feature 
    if (dummyi.ne.maxfeat) then !we should have had a fatal before, on maxfeat.eq.0 control if the only dummyi was <= 0
      write(ilog,*)
      write(ilog,*) 'Fatal. After maxfeat.eq.0 we have dummyi.ne.maxfeat and nfeats_cfile.eq.1. This is a bug.'
      call fexit()
    end if  
    if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
    allocate(ncdm_whichfeats(nfeats_cfile))
    ncdm_whichfeats(nfeats_cfile) = maxfeat
    if (ncdm_doas.EQV..true.) then  !specific checks for this case are required
      if (maxval(ncdm_whichfeats).gt.ncdm_nfeats_ini) then !featur index is out of range
        write(ilog,*)
        write(ilog,*) 'Fatal. Detected a feature indicator greater than the assumed number of features (NCDM_NFEATS) &
 &in input ascii file (NCDM_ASFILE). (1)'
        call fexit()
      end if
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(ncdm_nfeats_ini))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.
    else 
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(maxfeat))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.
    end if
    ncdm_nfeats = nfeats_cfile
    return !with one feature we are done, i.e. there was no need for sorting etc., as done below otherwise.
  end if
!
! allocate
  if (allocated(tmparr).EQV..true.) deallocate(tmparr)
  allocate(tmparr(nfeats_cfile))  !helper array for merge_sort
  tmparr(:) = 0
  if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
  allocate(ncdm_whichfeats(nfeats_cfile))  !will be old feature indexing of those that are part of the cfile
  ncdm_whichfeats(:) = 0
!
! read it first time
  iu = freeunit()
  call strlims(ncdm_nm_cfile,t1,t2)
  open(unit=iu,file=ncdm_nm_cfile(t1:t2),status='old')
  do ifeat=1,nfeats_cfile
    read(iu,*) tmparr(ifeat)  !original numbering of features
    ncdm_whichfeats(ifeat) = tmparr(ifeat)
  end do
  close(unit=iu)
!
! sort it
  call merge_isort(ldim=nfeats_cfile,up=atrue,list=ncdm_whichfeats,olist=tmparr,ilo=1,ihi=nfeats_cfile) !tmparr is gettin' sorted
!
! remove duplicates and nulls and negatives
  arethere_wrongs = .false.
  counter = 1
  if (tmparr(counter).gt.0) then                  !first check for negative and zeroes
    ncdm_whichfeats(counter) = tmparr(counter)    !get first value from sorted array tmparr
    startfeat = 2                                 !prepare to next
  else
    arethere_wrongs = .true.
    notpositive = .true.
    startfeat = 2
    do while (notpositive.EQV..true.) 
       if (startfeat.gt.nfeats_cfile) then
         write(ilog,*)
         write(ilog,*) 'Fatal. I do not know where to start in your cfile, there are no eligible positive integers...' !impossible
         call fexit()
       end if
       if (tmparr(startfeat).gt.0) then
         notpositive = .false.
         ncdm_whichfeats(counter) = tmparr(startfeat) !we start at the first that is gt 0 (tmparr is sorted already here)
         startfeat = startfeat + 1                    !prepare to next
       else
         startfeat = startfeat + 1
       end if
    end do
  end if
  if (nfeats_cfile.ge.startfeat) then
    do ifeat=startfeat,nfeats_cfile
      if (tmparr(ifeat).ne.tmparr(ifeat-1)) then  !not duplicated. At first cycle tmparr(ifeat-1) is already in ncdm_whichfeats(1)
        counter = counter + 1
        ncdm_whichfeats(counter) = tmparr(ifeat)  !getting not duplicated val in first free position -> ncdm_whichfeats gets sorted
      else                                        !just ignoring duplicated value
        arethere_wrongs = .true.
      end if
    end do
  end if
!
! some checks and possible fast return
  if ((arethere_wrongs.EQV..true.).AND.(nfeats_cfile.eq.counter)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten equal number of features before and after removing ascertained duplicates. This is a bug.'
    call fexit()
  end if 
  if ((arethere_wrongs.EQV..false.).AND.(nfeats_cfile.ne.counter)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten different number of features after removing unexisting duplicates. This is a bug.'
    call fexit()
  end if
  if (counter.eq.1) then !only one good value in ncdm_whichfeats(1), but more than one entry there (already out if nfeats_cfile=1)
    tmparr(1) = ncdm_whichfeats(1)
    if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
    allocate(ncdm_whichfeats(counter))
    ncdm_whichfeats(counter) = tmparr(1)
    ncdm_nfeats = counter
    nfeats_cfile = counter
    if (arethere_wrongs.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Gotten unexisting duplicates with counter 1. This is a bug as we should have returned earlier.'
      call fexit()
    else
      write(ilog,*) 'Warning. The input cfile (NCDM_CFILE) contained wrong feature indexes, which have been removed. &
 &It remains only a feature which is not wrong.'
    end if
    if (ncdm_doas.EQV..true.) then  !specific checks for this case are required
      if (maxval(ncdm_whichfeats).gt.ncdm_nfeats_ini) then !feature index is out of range
        write(ilog,*)
        write(ilog,*) 'Fatal. Detected a feature indicator greater than the assumed number of features (NCDM_NFEATS) &
 &in input ascii file (NCDM_ASFILE). (2)'
        call fexit()
      end if
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(ncdm_nfeats_ini))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.
    else 
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(maxfeat))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.
    end if
    return
  end if
  nfeats_cfile = counter      !counter is number of eligible features here
!
! checks before final assignment
  if (maxval(ncdm_whichfeats).ne.maxfeat) then !important to avoid possible bug on mask allocation 
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten maxval(ncdm_whichfeats).ne.maxfeat at the end of ncdm_read_featsfl. This is a bug.'
    call fexit()
  end if
  if (ncdm_doas.EQV..true.) then
    if (maxval(ncdm_whichfeats).gt.ncdm_nfeats_ini) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Detected a feature indicator greater than the assumed number of features (NCDM_NFEATS) &
 &in input ascii file (NCDM_ASFILE). (3)'
      call fexit()
    end if
  end if
!
! final assignment
  if (arethere_wrongs.EQV..true.) then    !then I have to reduce the dimension of ncdm_whichfeats, which was sorted before
    do ifeat=1,nfeats_cfile
      tmparr(ifeat) = ncdm_whichfeats(ifeat)  !storing only the right values (sorted and not duplicated)
    end do
    if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
    allocate(ncdm_whichfeats(nfeats_cfile))   !will be old feature indexing of those that are part of the cfile
    do ifeat=1,nfeats_cfile
      ncdm_whichfeats(ifeat) = tmparr(ifeat)  !storing the right value
    end do
    write(ilog,*) 'Warning. The input cfile (NCDM_CFILE) contained wrong features indexes, which have been removed.'
  end if
!
  if (allocated(tmparr).EQV..true.) deallocate(tmparr) 
!
! the mask
  if (ncdm_doas.EQV..true.) then  !this is important
    if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
    allocate(ncdm_featsmask(ncdm_nfeats_ini))
    ncdm_featsmask = .false.  
    ncdm_featsmask(ncdm_whichfeats) = .true.
  else  !this is useless, but I do it so that I do not have to distinguish the two cases below
    if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
    allocate(ncdm_featsmask(maxfeat))
    ncdm_featsmask = .false.  
    ncdm_featsmask(ncdm_whichfeats) = .true.
  end if
!
! last checks (mainly consistency within the routine)
  if (nfeats_cfile.ge.2) then
    do ifeat=2,nfeats_cfile
      if (ncdm_whichfeats(ifeat).eq.ncdm_whichfeats(ifeat-1)) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Gotten duplicated feature value at the end of the cleaning procedure. This is a bug or a &
  &mismatched data type in input cfile.'
        call fexit()
      end if
      if (ncdm_whichfeats(ifeat).le.0) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Gotten negative feature value at the end of the cleaning procedure. This is a bug or a &
  &mismatched data type in input cfile. (1)'
        call fexit()
      end if
    end do
!  else if (ncdm_whichfeats(1).le.0) then
!    write(ilog,*)
!    write(ilog,*) 'Fatal. Gotten negative feature value at the end of the cleaning procedure. This is a bug or a &
!  &mismatched data type in input cfile. (2)'
!    call fexit()
!  end if
  else
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten nfeats_cfile.ge.2 false. We should have exited earlier. This is a bug.'
    call fexit()
  end if
  if (nfeats_cfile.ne.size(ncdm_whichfeats)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten nfeats_cfile.ne.size(ncdm_whichfeats) at the end of ncdm_read_cfile. This is a bug. (2)'
    call fexit()
  end if
  if (nfeats_cfile.gt.maxval(ncdm_whichfeats)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten more features than maxval of indexes. This is a bug.'
    call fexit()
  end if
  if (ncdm_doas.EQV..true.) then
    if (nfeats_cfile.gt.ncdm_nfeats_ini) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Detected more features in CFILE than what specified with (NCDM_NFEATS). This is a bug.'
      call fexit()
    end if
  end if
!
! set already the right value to ncdm_nfeats
  ncdm_nfeats = nfeats_cfile 
!
end subroutine ncdm_read_cfile
!
!---------------------------------------------------------------------------
!
subroutine ncdm_manage_whichfeats(from_read_ncfl)
!
  use iounit
  use ncdm
!
  implicit none
!
  integer ifeat
!
  logical, INTENT(IN)  :: from_read_ncfl
!
  if ((ncdm_doas.EQV..true.).AND.((from_read_ncfl.EQV..false.))) then
    if (ncdm_isthere_cfile.EQV..true.) then
      write(ilog,*) 'Fatal. Inside ncdm_manage_whichfeats with the wrong flags combination. This is a bug. (1)'
      call fexit()
    end if
    ncdm_nfeats = ncdm_nfeats_ini  !will stay, cannot be changed by other keyword
    if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
    allocate(ncdm_whichfeats(ncdm_nfeats_ini))
    if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
    allocate(ncdm_featsmask(ncdm_nfeats_ini))
    do ifeat=1,ncdm_nfeats_ini
      ncdm_whichfeats(ifeat) = ifeat          !just use all of them
    end do
    ncdm_featsmask = .false.
    ncdm_featsmask(ncdm_whichfeats) = .true.  !set all of them to eligible
  else if ((ncdm_donc.EQV..true.).AND.(from_read_ncfl.EQV..true.)) then !it already set many things (e.g. ncdm_nfeats_ini)
    if (ncdm_isthere_cfile.EQV..true.) then !it already set many things (e.g. ncdm_nfeats)
      if (allocated(ncdm_whichfeats).EQV..false.) then
        write(ilog,*)   
        write(ilog,*) 'Fatal. Inside ncdm_manage_whichfeats with framesfile and with ncdm_whichfeats not allocated. This is a bug.'
        call fexit()
      end if 
      if (ncdm_nfeats_ini.lt.maxval(ncdm_whichfeats)) then  !only meaningful thing to do here is to check consistency
        write(ilog,*)
        write(ilog,*) 'Fatal. In the input net-cdf file, there are less features that what stated in the cfile. (1)'
        call fexit()
      end if
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(ncdm_nfeats_ini))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.  
    else  !so we are here for the first time without having been previously in ncdm_read_cfile and with no cfile
      ncdm_nfeats = ncdm_nfeats_ini  
      if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)  !no cfile, so I can do deallocate with no risk
      allocate(ncdm_whichfeats(ncdm_nfeats_ini))
      do ifeat=1,ncdm_nfeats_ini
        ncdm_whichfeats(ifeat) = ifeat          !just use all of them
      end do
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(ncdm_nfeats_ini))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.  !set all of them to eligible
    end if
  else if (ncdm_donc.EQV..false.) then
    write(ilog,*)   
    write(ilog,*) 'Fatal. Inside ncdm_manage_whichfeats with the wrong flags combination. This is a bug. (2)'
    call fexit()
  else     !if the call is from datasaw (which can happen only without cfile) and ncdm_donc is true,
    return !it does nothing as there is yet no information on the number of features in the net-cdf file
  end if 
!
end subroutine ncdm_manage_whichfeats
!
!---------------------------------------------------------------------------
!
subroutine ncdm_manage_ccollect(from_read_ncfl) !it should work also with frames file and ccollect
!
  use iounit
  use clusters, ONLY: cstorecalc
  use ncdm
!
  implicit none
!
  integer counter,nframes,iframe
!
  integer, ALLOCATABLE :: tmpvec(:)
  logical, INTENT(IN)  :: from_read_ncfl
!
  if ((ncdm_doas.EQV..true.).AND.(from_read_ncfl.EQV..false.)) then
    if (allocated(ncdm_whichframes).EQV..false.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got ncdm_whichframes not allocated in ncdm_manage_ccollect. This is a bug. (1)'
      call fexit()
    end if
    nframes = int(floor(1.*ncdm_nframes/cstorecalc))  !ncdm_nframes set by ncdm_manage_whichframes or framesfile (chaindata)
    allocate(tmpvec(nframes))                         !will be new ncdm_whichframes
    counter = 0
    do iframe=1,ncdm_nframes  !ncdm_nframes: either ncdm_nframes_ini or the number of frames in framesfile
      if (mod(iframe,cstorecalc).eq.0) then           !eligible
        counter = counter + 1                         !effective snapshot
        tmpvec(counter) = ncdm_whichframes(iframe)    !original numbering, possibly ccollecting a framesfile 
      end if
    end do
    if (nframes.ne.counter) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got a mismatch between counter and nframes in ncdm_manage_ccollect. This is a bug (1).'
      call fexit()
    end if
    deallocate(ncdm_whichframes)
    allocate(ncdm_whichframes(nframes))        !allocated to the used number of snapshot. Contains original (absl.) numbers
    do iframe=1,nframes
      ncdm_whichframes(iframe) = tmpvec(iframe)
    end do
    deallocate(tmpvec)
    ncdm_nframes = nframes !finally accounting for possible ccollect subsampling 
    ncdm_framesmask = .false.
    ncdm_framesmask(ncdm_whichframes) = .true. !here by now allocated to original number of frames in ascii, ncdm_nframes_ini
  else if ((ncdm_donc.EQV..true.).AND.(from_read_ncfl.EQV..true.)) then 
    if (allocated(ncdm_whichframes).EQV..false.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got ncdm_whichframes not allocated in ncdm_manage_ccollect. This is a bug. (2)'
      call fexit()
    end if
    nframes = int(floor(1.*ncdm_nframes/cstorecalc)) !ncdm_nframes by framesfl (chaindata) or by ncdm_manage_whichframes (read_ncfl)
    !if (nframes.le.3) then
    !  write(ilog,*) 
    !  write(ilog,*) 'Fatal. Due to the settings for CCOLLECT, there are no more than 3 snapshots left for analysis of net-cdf file.'
    !  call fexit()
    !end if
    allocate(tmpvec(nframes))                        !will be new ncdm_whichframes
    counter = 0
    do iframe=1,ncdm_nframes
      if (mod(iframe,cstorecalc).eq.0) then
        counter = counter + 1                         !effective snapshot
        tmpvec(counter) = ncdm_whichframes(iframe)    !original numbering, possibly ccollecting a framesfile 
      end if
    end do
    if (nframes.ne.counter) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got a mismatch between counter and nframes in ncdm_manage_ccollect. This is a bug (2)'
      call fexit()
    end if
    deallocate(ncdm_whichframes)
    allocate(ncdm_whichframes(nframes)) !allocated to the used number of snapshot. Contains original (absl.) numbers
    do iframe=1,nframes
      ncdm_whichframes(iframe) = tmpvec(iframe)
    end do
    deallocate(tmpvec)
    ncdm_nframes = nframes  !finally accounting for possible ccollect subsampling
    ncdm_framesmask = .false.
    ncdm_framesmask(ncdm_whichframes) = .true. !here by now allocated to original number of frames in netcdf, ncdm_nframes_ini
  else if ((ncdm_donc.EQV..false.)) then
    write(ilog,*)   
    write(ilog,*) 'Fatal. Inside ncdm_manage_ccollect with the wrong flags combination.'
    call fexit()
  else     !if the call is from datasaw and ncdm_donc is true,
    return !it does nothing as there is yet no information on the number of frames in the net-cdf file
  end if
!
end subroutine ncdm_manage_ccollect
!
!---------------------------------------------------------------------------
!
! I should split this routine in two. First read the ascii as it is and then rearrange cludata according to cdis layout
! I might want to delete all the checkas part...
!
subroutine ncdm_convert_ascii !from input ascii to cludata and possibly to net cdf by calling ncdm_write_ncfl
!
  use iounit
  use clusters, ONLY: cdis_crit,cludata,cstorecalc,cl_imvec
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
  use ncdm
!
  implicit none
!
  integer, parameter :: aone = 1        !to manage writing of check file according to cdistance
!
  integer freeunit,iu,iframe,counter,ifeat,featcounter
  integer jframe,endframe
  integer iomess,ipoint,howmany,ipointx3         !managing error in reading and string to reals
  integer dummyi                                 !helper integer
!
  integer, ALLOCATABLE :: featvals_inarr(:)
  integer, ALLOCATABLE :: dynwghts_inarr(:)
!
  integer, parameter :: charlen = 999999
!
  RTYPE dummyr
! 
  logical goon,afalse           !goon controls what to do while reading
!
  character(MAXSTRLEN) frmt
  character (len=:), ALLOCATABLE :: lnstrng
!
  RTYPE, ALLOCATABLE :: data_on_row(:)  !temporary vector to deal with reading only specific rows of input file
  RTYPE, ALLOCATABLE :: tmparr(:,:)     !temporary helper array to read data and store them when needed 
  RTYPE, ALLOCATABLE :: tmparrw(:,:)    !helper array for dynamic weight management
!
  iomess = 0
  goon = .false.
  afalse = .false.
!
! allocate the helper arrays 
  if (allocated(tmparr).EQV..true.) deallocate(tmparr)
  allocate(tmparr(ncdm_nfeats_ini,ncdm_nframes_ini))
  if (allocated(data_on_row).EQV..true.) deallocate(data_on_row) !to deal with possible feature reduction due to cfile
  allocate(data_on_row(ncdm_nfeats_ini))                         !the user-provided full rank
!
  iu = freeunit()
  open(unit=iu,file=ncdm_fn_as,status='old')  !opening ascii file
!
  if (ncdm_checkas.EQV..true.) then  !checking input string
    write(frmt,'("(a",I0,")")') charlen
    if (allocated(lnstrng).EQV..true.) deallocate(lnstrng)
    allocate(character(len=charlen) :: lnstrng)
    iframe = 1
    read(iu,frmt,iostat=iomess) lnstrng
    ipoint = 1
    call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
    if (howmany.ne.ncdm_nfeats_ini) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file or wrong format (e.g. scientific) (-1).'
      call fexit()
    else 
      call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse) !checking there is nothing else beyond
      if (howmany.ne.0) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file or wrong format (e.g. scientific) (0).'
        call fexit()
      end if 
    end if 
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Input ascii file appears to contain less lines than what specified with keyword NCDM_NFRAMES. (0)'
      call fexit()
    end if
    tmparr(1:ncdm_nfeats_ini,iframe) = data_on_row !values of features in array tmparr
    deallocate(lnstrng)
    ipointx3 = 3*ipoint                            !heuristic
    allocate(character(len=ipointx3) :: lnstrng)   !should be in the standard of most recent complilers
    write(frmt,'("(a",I0,")")') ipointx3
    do iframe=2,ncdm_nframes_ini                   !reading the values of the features in temporary array
      read(iu,frmt,iostat=iomess) lnstrng
      ipoint = 1
      call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
      if (howmany.ne.ncdm_nfeats_ini) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file or wrong format (e.g. scientific) (1).'
        write(ilog,*) lnstrng
        call fexit()
      else 
        call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse) !checking there is nothing else beyond
        if (howmany.ne.0) then
          write(ilog,*) 
          write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file or wrong format (e.g. scientific) (2).'
          call fexit()
        end if 
      end if 
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Input ascii file appears to contain less lines than what specified with keyword NCDM_NFRAMES. (1)'
        call fexit()
      end if
      tmparr(1:ncdm_nfeats_ini,iframe) = data_on_row  !values of features in upper half of array tmparr
    end do
!
    read(iu,frmt,iostat=iomess) lnstrng  !read one extra line, if not there it means there are no weights and do not go on reading
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Warning. Input ascii file does not contain any feature weight.'
      ncdm_arethere_dywghts = .false.
      ncdm_arethere_stwghts = .false.
      goon = .false.                  !finished, just copy temparr to cludata
    else  
      ipoint = 1
      call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse) !storing the extra line in data_on_row
      if (howmany.ne.ncdm_nfeats_ini) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (3).'
        call fexit()
      else 
        call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse)
        if (howmany.ne.0) then
          write(ilog,*) 
          write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (4).'
          call fexit()
        end if 
      end if 
      goon = .true.
    end if
!
    if (goon.EQV..true.) then  !there is at least one extra line 
      goon = .false.
      read(iu,frmt,iostat=iomess) lnstrng 
      if (iomess.eq.IOSTAT_END) then  !now we know there is only one extra line, i.e. there are only static weights
        ncdm_arethere_stwghts = .true.
        write(ilog,*) 'Warning. Assuming static weights from input file.'
        if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
        allocate(cl_imvec(1:ncdm_nfeats_ini))
        cl_imvec(:) = data_on_row 
      else
        ncdm_arethere_stwghts = .false.  !two lines detected after the last of frames, we assume there are dynamic weights
        ncdm_arethere_dywghts = .true.
        write(ilog,*) 'Warning. Assuming dynamic weights from input file.'
        if (allocated(tmparrw).EQV..true.) deallocate(tmparrw) !helper vector for storing weights
        allocate(tmparrw(ncdm_nfeats_ini,ncdm_nframes_ini))  !as many dyn weights as features for all snapshots
        iframe = 1
        tmparrw(1:ncdm_nfeats_ini,iframe) = data_on_row
        ipoint = 1
        call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
        if (howmany.ne.ncdm_nfeats_ini) then
          write(ilog,*) 
          write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (5).'
          call fexit()
        else 
          call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse)
          if (howmany.ne.0) then
            write(ilog,*) 
            write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (6).'
            call fexit()
          end if 
        end if 
        iframe = 2
        tmparrw(1:ncdm_nfeats_ini,iframe) = data_on_row
        do iframe=3,ncdm_nframes_ini  !already checked that are at least three frames. First two dynamic weights already read
          read(iu,frmt,iostat=iomess) lnstrng 
          if (iomess.eq.IOSTAT_END) then
            write(ilog,*)
            write(ilog,*) 'Fatal. Input ascii file appears to contain less weights than frames. (1)'
            call fexit()
          end if
          ipoint = 1
          call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
          if (howmany.ne.ncdm_nfeats_ini) then
            write(ilog,*) 
            write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (7).'
            call fexit()
          else 
            call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse)
            if (howmany.ne.0) then
              write(ilog,*) 
              write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (8).'
              call fexit()
            end if 
          end if 
          tmparrw(1:ncdm_nfeats_ini,iframe) = data_on_row
        end do
        read(iu,frmt,iostat=iomess) lnstrng  !looking for static weights at last line 
        if (iomess.ne.IOSTAT_END) then
          ipoint = 1
          call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
          if (howmany.ne.ncdm_nfeats_ini) then
            write(ilog,*) 
            write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (9).'
            call fexit()
          else 
            call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse)
            if (howmany.ne.0) then
              write(ilog,*) 
              write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file (10).'
              call fexit()
            end if 
          end if
          write(ilog,*) 'Warning. Assuming also static weights from input file.'
          ncdm_arethere_stwghts = .true.
          if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
          allocate(cl_imvec(ncdm_nfeats_ini))
          cl_imvec(:) = data_on_row
          read(iu,*,iostat=iomess) !just checking there is no an extra line...
          if (iomess.ne.IOSTAT_END) then
            write(ilog,*) 
            write(ilog,*) 'Fatal. Unrecognized total number of lines in input ascii file. It appears that there are more lines &
     &than assuming both dynamic and static weights. (1)'
            call fexit()
          end if
        else 
          write(ilog,*) 'Warning. No static weights detected in input file.'
        end if
      end if !dyn weights part
    end if !search for weights ends here
!
  else  !ncdm_checkas is false
!
    do iframe=1,ncdm_nframes_ini                !reading the values of the features in temporary array
      read(iu,*,iostat=iomess) data_on_row
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Input ascii file appears to contain less lines than what specified with keyword NCDM_NFRAMES. (2)'
        call fexit()
      end if
      tmparr(1:ncdm_nfeats_ini,iframe) = data_on_row  !values of features in upper half of array 
    end do
!
    read(iu,*,iostat=iomess) data_on_row  !read one extra line, if not there it means there are no weights and do not go on reading
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Warning. Input ascii file does not contain any feature weight.'
      ncdm_arethere_dywghts = .false.
      ncdm_arethere_stwghts = .false.
      goon = .false.                  !finished, just copy temparr to cludata
    else  
      goon = .true.
    end if
!
    if (goon.EQV..true.) then
      goon = .false.
      if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec) !just helper here
      allocate(cl_imvec(1:ncdm_nfeats_ini))
      cl_imvec(:) = data_on_row
      read(iu,*,iostat=iomess) data_on_row
      if (iomess.eq.IOSTAT_END) then  !now we know there is only one extra line, i.e. there are only static weights
        ncdm_arethere_stwghts = .true.
        write(ilog,*) 'Warning. Assuming static weights from input file.'
        goon = .false.                   
      else
        ncdm_arethere_stwghts = .false.  !two lines detected after the last of frames, we assume there are dynamic weights
        ncdm_arethere_dywghts = .true.
        write(ilog,*) 'Warning. Assuming dynamic weights from input file.'
        if (allocated(tmparrw).EQV..true.) deallocate(tmparrw) !helper vector for storing weights
        allocate(tmparrw(ncdm_nfeats_ini,ncdm_nframes_ini))  !as many dyn weights as features for all snapshots
        iframe = 1
        tmparrw(1:ncdm_nfeats_ini,iframe) = cl_imvec(:)
        deallocate(cl_imvec)
        iframe = 2
        tmparrw(1:ncdm_nfeats_ini,iframe) = data_on_row
        do iframe=3,ncdm_nframes_ini  !already checked that are at least three frames. First two dynamic weights already read
          read(iu,*,iostat=iomess) data_on_row
          if (iomess.eq.IOSTAT_END) then
            write(ilog,*)
            write(ilog,*) 'Fatal. Input ascii file appears to contain less weights than frames. (2)'
            call fexit()
          end if
          tmparrw(1:ncdm_nfeats_ini,iframe) = data_on_row
        end do
        read(iu,*,iostat=iomess) data_on_row  !looking for static weights at last line 
        if (iomess.ne.IOSTAT_END) then
          write(ilog,*) 'Warning. Assuming also static weights from input file.'
          ncdm_arethere_stwghts = .true.
          if (allocated(cl_imvec).EQV..true.) then
            write(ilog,*)
            write(ilog,*) 'Fatal. Found cl_imvec allocated in the wrong moment.'
            call fexit()
          end if
          allocate(cl_imvec(ncdm_nfeats_ini))
          cl_imvec(:) = data_on_row
          read(iu,*,iostat=iomess)  !just checkin extra line
          if (iomess.ne.IOSTAT_END) then
            write(ilog,*) 
            write(ilog,*) 'Fatal. Unrecognized total number of lines in input ascii file. It appears that there are more lines &
     &than assuming both dynamic and static weights. (2)'
            call fexit()
          end if
        else 
          write(ilog,*) 'Warning. No static weights detected in input file.'
        end if
      end if
    end if !search for weights ends here
!
  end if  !ncdm_checkas
!
  close(unit=iu)
!
! if write binary populate cludata for writing the input ascii in plain format
  if (ncdm_wrtinp.EQV..true.) then
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    if (ncdm_arethere_dywghts.EQV..true.) then
      dummyi = 9 !faking cdistance 9 to get style 9 since this is how we read the ascii for feats and weights
      call ncdm_set_cludatamasks(afalse,dummyi,ncdm_featvals_includata,ncdm_nfeats_ini,ncdm_dynwghts_includata)
      if (allocated(cludata).EQV..true.) deallocate(cludata)
      allocate(cludata(2*ncdm_nfeats_ini,ncdm_nframes_ini))
      if ((size(cludata,dim=1).ne.(2*size(tmparr,dim=1))).OR.(size(cludata,dim=1).ne.(2*size(tmparrw,dim=1)))) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Got mismatch on size cludata and helper arrays. This is a bug. (1)'
        call fexit()
      end if
      do iframe=1,ncdm_nframes_ini
        cludata(ncdm_featvals_includata,iframe) = tmparr(1:ncdm_nfeats_ini,iframe)
        cludata(ncdm_dynwghts_includata,iframe) = tmparrw(1:ncdm_nfeats_ini,iframe)
      end do
    else
      dummyi = 7 !faking cdistance 7 since this is how we read the ascii for feats
      call ncdm_set_cludatamasks(afalse,dummyi,ncdm_featvals_includata,ncdm_nfeats_ini) !no dyn weights
      if (allocated(cludata).EQV..true.) deallocate(cludata)
      allocate(cludata(ncdm_nfeats_ini,ncdm_nframes_ini))
      if ((size(cludata,dim=1).ne.(size(tmparr,dim=1))).OR.(allocated(tmparrw).EQV..true.)) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Got error on cludata or helper arrays. This is a bug. (1)'
        call fexit()
      end if
      do iframe=1,ncdm_nframes_ini
        cludata(ncdm_featvals_includata,iframe) = tmparr(1:ncdm_nfeats_ini,iframe)
      end do
    end if
!
!   a few checks
    if (ncdm_arethere_dywghts.EQV..true.) then !we have faked 9
      if ((size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1))) then !half values in cludata must be features
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1)) with ncdm_arethere_dywghts &
   &.true. in ncdm_convert_ascii. This is a bug.'
        call fexit()
      end if
      if ((size(cludata,dim=1)).ne.(2*size(ncdm_dynwghts_includata,dim=1))) then !half values in cludata must be weights
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1)) with cdis 2 or 9 in &
   &ncdm_convert_ascii. This is a bug.'
        call fexit()
      end if
    else !we have faked 7
      if ((size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1))) then !only features in cludata
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1)) ncdm_convert_ascii. This is a bug.'
        call fexit()
      end if
    end if
    call ncdm_write_ncfl(aone,afalse) !we always fake a cdistance different from 3 or 4 -> aone
  end if !(ncdm_wrtinp.EQV..true.)
!
! possible fast return if no analysis is to be performed
  if (ncdm_isdmonas.EQV..false.) then
    if (allocated(cludata).EQV..true.) deallocate(cludata)
    if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec) 
    if (allocated(data_on_row).EQV..true.) deallocate(data_on_row)
    if (allocated(tmparr).EQV..true.) deallocate(tmparr)
    if (allocated(tmparrw).EQV..true.) deallocate(tmparrw)
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    return
  end if
!
! finished reading and writing and we know what is in the ascii file. Allocate cludata, set masks 
  if (allocated(cludata).EQV..true.) deallocate(cludata)
  if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
  if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (cdis_crit.eq.2.OR.cdis_crit.eq.9) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)  
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))
    else if (cdis_crit.eq.4) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)  
      allocate(cludata(3*ncdm_nfeats,ncdm_nframes))
    else !since the cdistance does not support dyn weights, we discard them here
      ncdm_arethere_dywghts = .false.
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats)  
      if (cdis_crit.ne.3) then 
        allocate(cludata(ncdm_nfeats,ncdm_nframes)) 
      else
        allocate(cludata(2*ncdm_nfeats,ncdm_nframes)) 
      end if
    end if
  else if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)) then !no dyn weigths but cdist requires them
    if (cdis_crit.eq.2.OR.cdis_crit.eq.9) then
!     ncdm_arethere_dywghts F
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))
    else if (cdis_crit.eq.4) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)  
      allocate(cludata(3*ncdm_nfeats,ncdm_nframes))
    end if
  else 
    call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats) !no dyn weights
    if (cdis_crit.ne.3) then 
      allocate(cludata(ncdm_nfeats,ncdm_nframes))
    else 
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes)) 
    end if
  end if
  cludata(:,:) = 1 !to take uniformative flat weights by default 
!
! a few checks
  if (cdis_crit.eq.2.OR.cdis_crit.eq.9) then
    if (size(cludata,dim=1).ne.(2*ncdm_nfeats)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(2*ncdm_nfeats) with cdis 2 or 9 in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1))) then !half values in cludata must be features
      write(ilog,*)
      write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((size(cludata,dim=1)).ne.(2*size(ncdm_dynwghts_includata,dim=1))) then !half values in cludata must be weights
      write(ilog,*)
      write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_dynwghts_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  else if (cdis_crit.eq.4) then
    if (size(cludata,dim=1).ne.(3*ncdm_nfeats)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(3*ncdm_nfeats) with cdis 4 in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((size(cludata,dim=1)).ne.(3*size(ncdm_featvals_includata,dim=1)/2)) then !two thirds of values in cludata must be features
      write(ilog,*)
      write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(3*size(ncdm_featvals_includata,dim=1)/2) with cdis 4 in &
 &ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((size(cludata,dim=1)).ne.(3*size(ncdm_dynwghts_includata,dim=1))) then !one third of values in cludata must be weights
      write(ilog,*)
      write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(3*size(ncdm_dynwghts_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  else if (cdis_crit.ne.3) then
    if (size(cludata,dim=1).ne.(ncdm_nfeats)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(ncdm_nfeats) with cdis 1,7 or 8 in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((cdis_crit.ne.1).AND.(cdis_crit.ne.7).AND.(cdis_crit.ne.8)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got wrong cdistance in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1))) then !only features in cludata
      write(ilog,*)
      write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1)) with cdis 1,7 or 8 in &
 &ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  else !cdis_crit is 3
    if (cdis_crit.ne.3) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got cdis_crit not 3 in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if (size(cludata,dim=1).ne.(2*ncdm_nfeats)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(2*ncdm_nfeats) with cdis 3 in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if ((size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1))) then !only features in cludata
      write(ilog,*)
      write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1)) with cdis 3 in &
 &ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  end if
  if (size(ncdm_whichframes,dim=1).ne.ncdm_nframes) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(ncdm_whichframes,dim=1).ne.ncdm_nframes in ncdm_convert_ascii. This is a bug.'
    call fexit()
  end if
  if (size(cludata,dim=2).ne.ncdm_nframes) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(cludata,dim=2).ne.ncdm_nframes in ncdm_convert_ascii. This is a bug.'
    call fexit()
  end if
  if (size(ncdm_whichfeats,dim=1).ne.ncdm_nfeats) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(ncdm_whichfeats,dim=1).ne.ncdm_nfeats in ncdm_convert_ascii. This is a bug.'
    call fexit()
  end if
  if (maxval(ncdm_whichframes).gt.size(tmparr,dim=2)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got maxval(ncdm_whichframes).gt.size(tmparr,dim=2) in ncdm_convert_ascii. This is a bug.'
    call fexit()
  end if 
  if (maxval(ncdm_whichfeats).gt.size(tmparr,dim=1)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got maxval(ncdm_whichfeat).gt.size(tmparr,dim=2) in ncdm_convert_ascii. This is a bug.'
    call fexit()
  end if
  if (ncdm_arethere_dywghts.EQV..true.) then  !really there in the ascii
    if (maxval(ncdm_whichframes).gt.size(tmparrw,dim=2)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got maxval(ncdm_whichframes).gt.size(tmparrw,dim=2) in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if 
    if (maxval(ncdm_whichfeats).gt.size(tmparrw,dim=1)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got maxval(ncdm_whichfeat).gt.size(tmparrw,dim=2) in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    if (size(ncdm_dynwghts_includata,dim=1).ne.ncdm_nfeats) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(ncdm_dynwghts_includata,dim=1).ne.ncdm_nfeats in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  end if
  if (cdis_crit.ne.3.AND.cdis_crit.ne.4) then
    if (size(ncdm_featvals_includata,dim=1).ne.ncdm_nfeats) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(ncdm_featvals_includata,dim=1).ne.ncdm_nfeats in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  else if (cdis_crit.eq.3.OR.cdis_crit.eq.4) then
    if (size(ncdm_featvals_includata,dim=1).ne.(2*ncdm_nfeats)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(ncdm_featvals_includata,dim=1).ne.(2*ncdm_nfeats) in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
  end if
!
! populate cludata for good
  do iframe=1,ncdm_nframes 
    featcounter = 0
    do ifeat=1,ncdm_nfeats 
      featcounter = featcounter + 1
      cludata(ncdm_featvals_includata(featcounter),iframe) = tmparr(ncdm_whichfeats(ifeat),ncdm_whichframes(iframe))
      if (ncdm_arethere_dywghts.EQV..true.) then
        cludata(ncdm_dynwghts_includata(ifeat),iframe) = tmparrw(ncdm_whichfeats(ifeat),ncdm_whichframes(iframe))
      end if
      if (cdis_crit.eq.3.OR.cdis_crit.eq.4) then !duplicate feature val
        featcounter = featcounter + 1
        cludata(ncdm_featvals_includata(featcounter),iframe) = tmparr(ncdm_whichfeats(ifeat),ncdm_whichframes(iframe))
      end if
    end do
  end do
!
! manage static weights
  if (size(ncdm_whichfeats,dim=1).ne.ncdm_nfeats) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(ncdm_whichfeats,dim=1).ne.ncdm_nfeats in ncdm_convert_ascii at end. This is a bug.'
    call fexit()
  end if
  if ((cdis_crit.eq.8).AND.(ncdm_arethere_stwghts.EQV..true.)) then
    tmparr(1:ncdm_nfeats,1) = cl_imvec(ncdm_whichfeats) 
    deallocate(cl_imvec)
    allocate(cl_imvec(ncdm_nfeats))
    cl_imvec(:) = tmparr(1:ncdm_nfeats,1)
  else if ((cdis_crit.eq.8).AND.(ncdm_arethere_stwghts.EQV..false.)) then
    if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
    allocate(cl_imvec(ncdm_nfeats))
    cl_imvec(:) = 1 !flat
    ncdm_arethere_stwghts = .true.
  else if ((cdis_crit.ne.8).AND.(ncdm_arethere_stwghts.EQV..true.)) then !if cdistance is not 8, forget about them
    ncdm_arethere_stwghts = .false.
    deallocate(cl_imvec)
  end if
!
! manage dynamic weights
  if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)).AND.(ncdm_arethere_dywghts.EQV..false.)) then !flat but there
    ncdm_arethere_dywghts = .true.
  end if
!
  if (allocated(data_on_row).EQV..true.) deallocate(data_on_row)
  if (allocated(tmparr).EQV..true.) deallocate(tmparr)
  if (allocated(tmparrw).EQV..true.) deallocate(tmparrw)
!
end subroutine ncdm_convert_ascii
!
!---------------------------------------------------------------------------
!
! index of feats and wghts in cludata according to nfeats and cdis
subroutine ncdm_set_cludatamasks(frommain,cdis,featsvec,nfeats,dywghtsvec)
!
  use iounit
  use ncdm, ONLY: ncdm_donc,ncdm_doas
!
  implicit none
!
  logical, INTENT(IN) :: frommain
!
  integer, INTENT(IN)  :: cdis
  integer, INTENT(IN)  :: nfeats
!
  integer, ALLOCATABLE, INTENT(OUT) :: featsvec(:)
  integer, ALLOCATABLE, INTENT(OUT), OPTIONAL:: dywghtsvec(:)
!
  integer featcounter,ifeat
!
  if (allocated(featsvec).EQV..true.) then
    write(ilog,*)
    write(ilog,*) 'Fatal. allocated(featsvec).EQV..true. This is a bug.' 
    call fexit()
  end if
  if (present(dywghtsvec).EQV..true.) then
    if (allocated(dywghtsvec).EQV..true.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. allocated(dywghtsvec).EQV..true. This is a bug.' 
      call fexit()
    end if
  end if
!
  if (frommain.EQV..false.) then
!   masks depending on cdistance
    if ((cdis.eq.2).OR.(cdis.eq.4).OR.(cdis.eq.9)) then !the ones with dyn weights
      if (present(dywghtsvec).EQV..false.) then
        write(ilog,*)
        write(ilog,*) 'Fatal. present(dywghtsvec).EQV..false. with cdis that requires dyn weights (1). This is a bug.' 
        call fexit()
      end if
      if (cdis.eq.2) then
        allocate(featsvec(nfeats))
        allocate(dywghtsvec(nfeats))
        do ifeat=1,nfeats
          featsvec(ifeat) = 2*ifeat - 1
          dywghtsvec(ifeat) = 2*ifeat
        end do
      else if (cdis.eq.4) then
        allocate(featsvec(2*nfeats))
        allocate(dywghtsvec(nfeats))
        featcounter = 0
        do ifeat=1,nfeats
          featcounter = featcounter + 1
          featsvec(featcounter) = 3*ifeat - 2
          featcounter = featcounter + 1
          featsvec(featcounter) = 3*ifeat - 1
          dywghtsvec(ifeat) = 3*ifeat
        end do
      else if (cdis.eq.9) then
        allocate(featsvec(nfeats))
        allocate(dywghtsvec(nfeats))
        do ifeat=1,nfeats
          featsvec(ifeat) = ifeat
          dywghtsvec(ifeat) = ifeat + nfeats
        end do
      end if
    else
      if (present(dywghtsvec).EQV..true.) then
        write(ilog,*)
        write(ilog,*) 'Fatal. present(dywghtsvec).EQV..true. with cdis that does not require dyn weights (1). This is a bug.'
        call fexit()
      end if
      if (cdis.ne.3) then
        allocate(featsvec(nfeats))
        do ifeat=1,nfeats
          featsvec(ifeat) = ifeat
        end do
      else
        allocate(featsvec(2*nfeats))
        do ifeat=1,2*nfeats
          featsvec(ifeat) = ifeat
        end do
      end if
    end if
  else !we are calling from clustering.f90 
!   masks depending on cdistance
    if ((cdis.eq.2).OR.(cdis.eq.4).OR.(cdis.eq.9)) then !the ones with dyn weights
      if (present(dywghtsvec).EQV..false.) then
        write(ilog,*)
        write(ilog,*) 'Fatal. present(dywghtsvec).EQV..false. with cdis that requires dyn weights (2). This is a bug.' 
        call fexit()
      end if
      if (cdis.eq.2) then
        allocate(featsvec(nfeats))
        allocate(dywghtsvec(nfeats))
        do ifeat=1,nfeats
          featsvec(ifeat) = 2*ifeat - 1
          dywghtsvec(ifeat) = 2*ifeat
        end do
      else if (cdis.eq.4) then
        if (mod(nfeats,2).ne.0) then 
          write(ilog,*)
          write(ilog,*) 'Fatal. Gotten nfeats not divisible by 2 with CDIST 4. This is a bug.' 
          call fexit()
        end if
        allocate(featsvec(nfeats))   !nfeats already counts the duplicated features
        allocate(dywghtsvec(nfeats/2))
        featcounter = 0
        do ifeat=1,nfeats/2
          featcounter = featcounter + 1
          featsvec(featcounter) = 3*ifeat - 2
          featcounter = featcounter + 1
          featsvec(featcounter) = 3*ifeat - 1
          dywghtsvec(ifeat) = 3*ifeat
        end do
      else if (cdis.eq.9) then
        allocate(featsvec(nfeats))
        allocate(dywghtsvec(nfeats))
        do ifeat=1,nfeats
          featsvec(ifeat) = ifeat
          dywghtsvec(ifeat) = ifeat + nfeats
        end do
      end if
    else
      if (present(dywghtsvec).EQV..true.) then
        write(ilog,*)
        write(ilog,*) 'Fatal. present(dywghtsvec).EQV..true. with cdis that does not require dyn weights (2). This is a bug.'
        call fexit()
      end if
      if ((cdis.eq.3).AND.(mod(nfeats,2).ne.0)) then 
        write(ilog,*)
        write(ilog,*) 'Fatal. Gotten nfeats not divisible by 2 with CDIST 3. This is a bug.' 
        call fexit()
      end if
      allocate(featsvec(nfeats))
      do ifeat=1,nfeats
        featsvec(ifeat) = ifeat
      end do
    end if
  end if
!
end subroutine ncdm_set_cludatamasks
!
!---------------------------------------------------------------------------
!
subroutine ncdm_write_ncfl(featstep,dumpname)
!
  use iounit
  use system, ONLY: basename,bleng
  use clusters, ONLY: cdis_crit,cl_imvec,cludata
  use netcdf
  use ncdm
!
  implicit none
!
  integer, parameter :: nvar = 3      !how many variable do we need, 1 is the featureval cludataay, 2 and 3 are weights
  integer, parameter :: ndim = 2      !number of dimensions that we need to index our variables
!
  integer, INTENT(IN) :: featstep     !to define the step in do write to deal with duplicated features for cdist 3 and 4
  logical, INTENT(IN) :: dumpname     !to allow writing to a different non-NCDM filename
  integer t1,t2                       !to determine string limits
  integer i,freeunit,counter,ifeat,iframe
  integer nframes,nfeatsnc,nfeatsclu  !local values 
  integer ncid_db                     !net-cdf databese id
  integer ncid_dim(ndim)              !net-cdf dimensions ids
  integer ncid_var(nvar)              !net-cdf variables ids
  character(MAXSTRLEN) attstring
  character(len=(max(MAXSTRLEN,bleng+11))) ncfn_w    !net-cdf database (file) name to be written
  logical exists
!
  integer, ALLOCATABLE :: featvals_includata(:) !local to deal with duplicated features for cdist 3 and 4
!
! check
  if (featstep.lt.1) then
    write(ilog,*)
    write(ilog,*) 'Fatal. featstep.lt.1. This is a bug.'
    call fexit()
  end if
  if (allocated(ncdm_featvals_includata).EQV..false.) then
    write(ilog,*)
    write(ilog,*) 'Fatal. In ncdm_write_ncfl got ncdm_featvals_includata not allocated. This is a bug.'
  end if
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (allocated(ncdm_dynwghts_includata).EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. In ncdm_write_ncfl got ncdm_dynwghts_includata not allocated with dyn weights present. This is a bug.'
    end if
  end if
!
! initialize 
  nfeatsclu = size(ncdm_featvals_includata,dim=1) !number of features in cludata, including duplicates
  if (featstep.ne.1) then
    counter = 0
    do i=1,nfeatsclu,featstep
      counter = counter + 1
    end do
    nfeatsnc = counter !effective number of features that will be written down
    allocate(featvals_includata(nfeatsnc))
    counter = 0
    do i=1,nfeatsclu,featstep
      counter = counter + 1
      featvals_includata(counter) = ncdm_featvals_includata(i) !subsampling to discard duplicated dihedrals only
    end do
  else 
    nfeatsnc = nfeatsclu
    allocate(featvals_includata(nfeatsnc))
    do i=1,nfeatsclu
      featvals_includata(i) = ncdm_featvals_includata(i) !local copy
    end do
  end if
  nframes = size(cludata,dim=2)
!
  if (dumpname.EQV..true.) then
    ncfn_w = "CLUSTERING_FEATURES.nc"
    call strlims(ncfn_w,t1,t2)
  else
    if (ncdm_doas.EQV..true.) then                 
      t1 = 1
      t2 = bleng + 11
      ncfn_w(t1:t2) = basename(1:bleng) // "_convert.nc"
    else
      t1 = 1
      t2 = bleng + 11
      ncfn_w(t1:t2) = basename(1:bleng) // "_checked.nc"
    end if
  end if
!
  ncid_db = 0
  ncid_dim = 0
!
  inquire(file=ncfn_w(t1:t2),exist=exists)
  if (exists.EQV..true.) then
    ncid_db = freeunit()
    open(unit=ncid_db,file=ncfn_w(t1:t2),status='old')
    close(unit=ncid_db,status='delete')
  end if
!
! open the (potentially large) data file
  call check_fileionetcdf( nf90_create(path=ncfn_w(t1:t2), cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid_db) )
!
! global attributes
  attstring(1:19) = "CAMPARI DATA MINING"
  call check_fileionetcdf( nf90_put_att(ncid_db, NF90_GLOBAL, "program", attstring(1:19)) )
!
! dimensions
  attstring(1:6) = "nfeats" 
  call check_fileionetcdf( nf90_def_dim(ncid_db, attstring(1:6), nfeatsnc, ncid_dim(1)) )  
  attstring(1:7) = "nframes"
  call check_fileionetcdf( nf90_def_dim(ncid_db, attstring(1:7), NF90_UNLIMITED, ncid_dim(2)) )
!
! define variables and put relevat attributes
! featurevals variable
  call strlims(ncdm_nm_fvar,t1,t2)
  attstring(t1:t2) = ncdm_nm_fvar
  call check_fileionetcdf( nf90_def_var(ncid_db, attstring(t1:t2), NF90_FLOAT, (/ ncid_dim(1), ncid_dim(2) /), ncid_var(1)) )
  if ( (ncdm_is_periodic.EQV..true.).OR.( (dumpname.EQV..true.).AND.(cdis_crit.le.2) ) ) then 
    attstring(1:14) = "periodic_range"
    call check_fileionetcdf( nf90_put_att(ncid_db, ncid_var(1), attstring(1:14), ncdm_periodic_rng) ) 
    if (((ncdm_donc.EQV..true.).OR.(ncdm_isdmonas.EQV..true.)).AND.((cdis_crit.eq.1).OR.(cdis_crit.eq.2))) then
    if ( ( minval(cludata(featvals_includata,:)).lt.ncdm_periodic_rng(1) ).OR.( &
 &         maxval(cludata(featvals_includata,:)).gt.ncdm_periodic_rng(2)      ) ) then
        write(ilog,*) 'Warning. The provided input data exceed the reference periodic interval. This is not corrected &
 &in output files, but it is taken care of during analysis with periodic CDISTANCE, as specified.'
      end if
    end if
  end if
! possible weights
  if ( (ncdm_arethere_stwghts.EQV..true.).OR.( (dumpname.EQV..true.).AND.(cdis_crit.eq.8) ) ) then !output them
    attstring(1:14) = "static_weights"
    call check_fileionetcdf( nf90_def_var(ncid_db, attstring(1:14), NF90_FLOAT, ncid_dim(1), ncid_var(2)) )
  end if
  if ( (ncdm_arethere_dywghts.EQV..true.).OR.&
&      ( (dumpname.EQV..true.).AND.((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)) ) ) then 
    attstring(1:15) = "dynamic_weights"
    call check_fileionetcdf( nf90_def_var(ncid_db, attstring(1:15), NF90_FLOAT, (/ ncid_dim(1), ncid_dim(2) /), ncid_var(3)) )
  end if
!
! quit define mode
  call check_fileionetcdf( nf90_enddef(ncid_db) )
  call check_fileionetcdf( nf90_sync(ncid_db) ) 
!
! populate the variables
! cludata
  call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(1), values = cludata(featvals_includata,1:nframes), &
 &                                       start = (/ 1, 1 /), count = (/ nfeatsnc, nframes /) ) )
! static weights
  if ( (ncdm_arethere_stwghts.EQV..true.).OR.( (dumpname.EQV..true.).AND.(cdis_crit.eq.8) ) ) then
    if (size(cl_imvec,dim=1).ne.nfeatsnc) then
      write(ilog,*)
      write(ilog,*) 'Fatal. size(cl_imvec,dim=1).ne.nfeatsnc. This is a bug.'
      call fexit()
    end if
    if (dumpname.EQV..true.) then
      call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(2), values = cl_imvec(1:nfeatsnc)**2, &
 &                                           start = (/ 1 /) , count =  (/ nfeatsnc /)  ) ) 
    else 
      call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(2), values = cl_imvec(1:nfeatsnc), &
 &                                           start = (/ 1 /) , count =  (/ nfeatsnc /)  ) ) 
    end if
  end if 
! dynamic weights
  if ( (ncdm_arethere_dywghts.EQV..true.).OR.&
 &     ( (dumpname.EQV..true.).AND.((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)) ) ) then
    if ((dumpname.EQV..true.).AND.(cdis_crit.eq.4)) then !I do not want to make a new dimension for the weights. I duplicate them
      if ((2*size(ncdm_dynwghts_includata,dim=1)).ne.nfeatsnc) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. (2*size(ncdm_dynwghts_includata,dim=1)).ne.nfeatsnc). Wrong call to ncdm_set_cludatamasks before? &
 &This is a bug.'
        call fexit()
      end if
      do iframe=1,nframes
        counter = 0
        do ifeat=1,size(ncdm_dynwghts_includata,dim=1) !less weights than features
          counter = counter + 1
          call check_fileionetcdf(nf90_put_var(ncid_db, ncid_var(3),values=cludata(ncdm_dynwghts_includata(ifeat),iframe:iframe),&
 &                                start = (/ counter, iframe /), count = (/ 1, 1 /) ))
          counter = counter + 1
          call check_fileionetcdf(nf90_put_var(ncid_db, ncid_var(3),values=cludata(ncdm_dynwghts_includata(ifeat),iframe:iframe),&
 &                                start = (/ counter, iframe /), count = (/ 1, 1 /) ))
        end do
      end do
    else
      if (size(ncdm_dynwghts_includata,dim=1).ne.nfeatsnc) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. size(ncdm_dynwghts_includata,dim=1).ne.nfeatsnc. This is a bug.'
        call fexit()
      end if
      call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(3), values = cludata(ncdm_dynwghts_includata,1:nframes),&
 &                                           start = (/ 1, 1 /), count = (/ nfeatsnc, nframes /) ) )
    end if
  end if 
!
  call check_fileionetcdf( nf90_close(ncid_db) )
!
  if (allocated(featvals_includata).EQV..true.) deallocate(featvals_includata)
!
end subroutine ncdm_write_ncfl
!
!---------------------------------------------------------------------------
!
subroutine ncdm_read_ncfl
!
  use iounit
  use clusters, ONLY: cdis_crit,cstorecalc,cludata,cl_imvec
  use netcdf
  use ncdm
!
  implicit none
!
  integer, parameter :: nvar = 3        !variable ids, 1 is featurevals, 2 and 3 are possible weights
  integer, parameter :: ndim = 2        !number of dimensions that we need to index our variables
  integer, parameter :: aone = 1        !to manage writing of check file according to cdistance
  integer, parameter :: atwo = 2        !if cdistance is 3 or 4
!
  integer t1,t2,featcounter,ifeat,iframe,iatt,ivar,dummyi
  integer ncid_db                       !net-cdf databese id
  integer ncid_dim(ndim)                !net-cdf dimensions ids
  integer ncid_var(nvar)                !net-cdf variables ids
  integer nvar_db                       !how many variables are in the databased (used to look for weights)
  integer natts_vardb                   !number of attributes of db variable
!
  logical atrue,afalse
!
  character(MAXSTRLEN) attstring
!
  RTYPE, ALLOCATABLE :: tmparr(:,:)  !array that helps dealing with cfile and ccollect while keeping memory requirement low
!
  ncid_db = 0
  ncid_dim = 0
  ncid_var = 0
  attstring(:) = " "
  atrue = .true.
  afalse = .false.
!
! checks 
  if (ncdm_donc.EQV..false.)  then
    write(ilog,*)
    write(ilog,*) 'Fatal. This is a bug. ncdm_read_ncfl is called with ncdm_donc as .false.' 
    call fexit()
  end if
!
! open
  call check_fileionetcdf( nf90_open(ncdm_fn_r, NF90_NOWRITE, ncid_db) ) !read access only
!
! find all the necessary infos from the variable name
  call strlims(ncdm_nm_fvar,t1,t2)
  call check_fileionetcdf( nf90_inq_varid(ncid_db, ncdm_nm_fvar(t1:t2), ncid_var(1)) )
  call check_fileionetcdf( nf90_Inquire_Variable(ncid_db, ncid_var(1), dimids = ncid_dim , nAtts = natts_vardb ) ) 
  if (natts_vardb.gt.0) then  !number of attributes associated to variable 1
   do iatt=1,natts_vardb
     call check_fileionetcdf(nf90_inq_attname( ncid_db, ncid_var(1), iatt, attstring ) )
     call strlims(attstring,t1,t2)
     if (attstring(t1:t2).eq.'periodic_range') then
       call check_fileionetcdf( nf90_get_att(ncid_db, ncid_var(1), 'periodic_range', ncdm_periodic_rng) )
       ncdm_is_periodic = .true.
     end if
   end do 
  end if
!
! frames
  ncdm_nframes_ini = 0  
  call check_fileionetcdf( nf90_Inquire_Dimension(ncid_db, ncid_dim(2), len=ncdm_nframes_ini) )
  call ncdm_manage_whichframes(atrue)
!
! features
  ncdm_nfeats_ini = 0  
  call check_fileionetcdf( nf90_Inquire_Dimension(ncid_db, ncid_dim(1), len=ncdm_nfeats_ini) )
  call ncdm_manage_whichfeats(atrue)
!
! ccollect
  call ncdm_manage_ccollect(atrue)
!
! weights (if present)
  ncdm_arethere_stwghts = .false.
  ncdm_arethere_dywghts = .false.
  call check_fileionetcdf( nf90_Inquire(ncid_db, nVariables = nvar_db) )
  do ivar=1,nvar_db
    call check_fileionetcdf( nf90_Inquire_Variable(ncid_db, ivar, attstring) )
    if (attstring(1:14).eq.'static_weights') then
      ncdm_arethere_stwghts = .true.
      ncid_var(2) = ivar
    else if (attstring(1:15).eq.'dynamic_weights') then
      ncdm_arethere_dywghts = .true.
      ncid_var(3) = ivar
    end if
  end do
!
! some checks and fixes
  if ((ncdm_is_periodic.EQV..false.).AND.(cdis_crit.le.4)) then
    ncdm_is_periodic = .true.
    write(ilog,*) 'Warning. The current setting of CDISTANCE requires that the data are periodic. A periodic range between &
 &-180 and 180 will be used instead to circumvent the lack of specification or identification in input net-cdf file.'
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(ncdm_arethere_dywghts.EQV..true.)) then 
    write(ilog,*)
    write(ilog,*) 'Warning. Gotten both weights present in input net-cdf file. Only one type (if any) will be used in analysis, &
 &depending on the value of CDISTANCE.'
  end if
  if ((cdis_crit.eq.8).AND.(ncdm_arethere_stwghts.EQV..false.)) then
    write(ilog,*)
    write(ilog,*) 'Warning. The current setting of CDISTANCE to 8 requires the specification of static weigths. &
 &Flat (uniformative) weights will be used instead to circumvent the lack of specification or identification. &
 &They might become informative in analysis depending on keyword CMODWEIGHTS and related ones.'
    if (allocated(cl_imvec).EQV..true.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got cl_imvec allocated with arethere_stwghts .false. in ncdm_read_ncfl. This is a bug.'
      call fexit()
    end if
    allocate(cl_imvec(ncdm_nfeats))
    cl_imvec(:) = 1 !take flat ones
  end if
  if ( ( (cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9) ).AND.(ncdm_arethere_dywghts.EQV..false.) ) then  
    write(ilog,*) 'Warning. The current setting of CDISTANCE (to 2, 4 or 9) requires the specification of dynamic weigths. &
 &Flat (uniformative) weights will be used instead to circumvent the lack of specification or identification. &
 &They might become informative in analysis depending on keyword CMODWEIGHTS and related ones.'
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(cdis_crit.ne.8)) then  
    write(ilog,*) 'Warning. The detected static weights in the input net-cdf file will be discarded. This entails &
  &both the possible check-file and the analysis and is due to the value of CDISTANCE.'
    ncdm_arethere_stwghts = .false.
  end if
  if ((ncdm_arethere_dywghts.EQV..true.).AND.((cdis_crit.ne.2).AND.(cdis_crit.ne.4).AND.(cdis_crit.ne.9))) then 
    write(ilog,*) 'Warning. The detected dynamic weights in the input net-cdf file will be discarded. This entails &
  &both the possible check-file and the analysis and is due to the value of CDISTANCE.'
    ncdm_arethere_dywghts = .false.
  end if
!
! allocate cludata, set masks 
  if (allocated(cludata).EQV..true.) deallocate(cludata)
  if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
  if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (cdis_crit.eq.2.OR.cdis_crit.eq.9) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)  
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))
    else if (cdis_crit.eq.4) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)  
      allocate(cludata(3*ncdm_nfeats,ncdm_nframes))
    else 
      write(ilog,*)
      write(ilog,*) 'Fatal. Got ncdm_arethere_dywghts true with wrong cdistance. This is a bug.'
      call fexit()
    end if
  else if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)) then !no dyn weigths but cdist requires them
    if (cdis_crit.eq.2.OR.cdis_crit.eq.9) then
!     ncdm_arethere_dywghts F
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))
    else if (cdis_crit.eq.4) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)  
      allocate(cludata(3*ncdm_nfeats,ncdm_nframes))
    end if
  else 
    if (cdis_crit.ne.3) then
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats) !no dyn weights
      allocate(cludata(ncdm_nfeats,ncdm_nframes))
    else 
      call ncdm_set_cludatamasks(afalse,cdis_crit,ncdm_featvals_includata,ncdm_nfeats) !no dyn weights
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))
    end if
  end if
  cludata(:,:) = 1 !to take uniformative flat weights by default 
!
! populate cludata
  do iframe=1,ncdm_nframes
    featcounter = 0
    do ifeat=1,ncdm_nfeats
      featcounter = featcounter + 1
      call check_fileionetcdf( nf90_get_var(ncid_db, ncid_var(1), cludata(ncdm_featvals_includata(featcounter),iframe:iframe), &
 &                             start = (/ ncdm_whichfeats(ifeat), ncdm_whichframes(iframe) /), count = (/ 1, 1 /) )      )
      if (cdis_crit.eq.3.OR.cdis_crit.eq.4) then !duplicate feature val
        featcounter = featcounter + 1
        cludata(ncdm_featvals_includata(featcounter),iframe) = cludata(ncdm_featvals_includata(featcounter-1),iframe)
      end if
    end do
  end do 
!
! in case, account for weights
  if (ncdm_arethere_stwghts.EQV..true.) then !meaning there were in the netcdf file and cdistance is 8, otherwise no way to enter 
    allocate(cl_imvec(ncdm_nfeats))
    do ifeat=1,ncdm_nfeats
      call check_fileionetcdf( nf90_get_var(ncid_db, ncid_var(2), cl_imvec(ifeat:ifeat), &
 &                             start = (/ ncdm_whichfeats(ifeat) /), count = (/ 1 /) ))
    end do
  end if
  if ((ncdm_arethere_dywghts.EQV..true.)) then !meaning there were in the netcdf file and cdistance is 2, 4 or 9
    do iframe=1,ncdm_nframes
      do ifeat=1,ncdm_nfeats
        call check_fileionetcdf( nf90_get_var(ncid_db, ncid_var(3), cludata(ncdm_dynwghts_includata(ifeat),iframe:iframe), &
 &                               start = (/ ncdm_whichfeats(ifeat), ncdm_whichframes(iframe) /), count = (/ 1, 1 /) ))
      end do
    end do 
  end if
!
  if (ncdm_wrtinp.EQV..true.) then 
    if (cdis_crit.ne.3.AND.cdis_crit.ne.4) then
      call ncdm_write_ncfl(aone,afalse)
    else 
      call ncdm_write_ncfl(atwo,afalse)
    end if
  end if 
!
! final set (we do not turn on weights before to not write flat weights in net-cdf output for debug)
  if ((cdis_crit.eq.8).AND.(ncdm_arethere_stwghts.EQV..false.)) then
    ncdm_arethere_stwghts = .true.
  end if
  if (((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)).AND.(ncdm_arethere_dywghts.EQV..false.)) then
    ncdm_arethere_dywghts = .true.
  end if
!
! a few checks
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (cdis_crit.eq.2.OR.cdis_crit.eq.9) then
      if (size(cludata,dim=1).ne.(2*ncdm_nfeats)) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(2*ncdm_nfeats) with cdis 2 or 9 in ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
      if ((size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1))) then !half values in cludata must be features
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
      if ((size(cludata,dim=1)).ne.(2*size(ncdm_dynwghts_includata,dim=1))) then !half values in cludata must be weights
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_dynwghts_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
    else if (cdis_crit.eq.4) then
      if (size(cludata,dim=1).ne.(3*ncdm_nfeats)) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(3*ncdm_nfeats) with cdis 4 in ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
      if ((size(cludata,dim=1)).ne.(3*size(ncdm_featvals_includata,dim=1)/2)) then !two third of values in cludata must be features
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(2*size(ncdm_featvals_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
      if ((size(cludata,dim=1)).ne.(3*size(ncdm_dynwghts_includata,dim=1))) then !one third values in cludata must be weights
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(3*size(ncdm_dynwghts_includata,dim=1)) with cdis 2 or 9 in &
 &ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
    else 
      write(ilog,*)
      write(ilog,*) 'Fatal. Got ncdm_arethere_dywghts true with wrong cdistance in ncdm_read_ncfl. This is a bug.'
      call fexit()
    end if
  else 
    if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got ncdm_arethere_dywghts false with cdis 2, 4 or 9 in ncdm_read_ncfl. This is a bug.'
      call fexit()
    else
      if (cdis_crit.ne.3) then
        if (size(cludata,dim=1).ne.(ncdm_nfeats)) then
          write(ilog,*)
          write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(ncdm_nfeats) with cdis 1,7 or 8 in ncdm_read_ncfl. This is a bug.'
          call fexit()
        end if
        if ((cdis_crit.ne.1).AND.(cdis_crit.ne.7).AND.(cdis_crit.ne.8)) then
          write(ilog,*)
          write(ilog,*) 'Fatal. Got wrong cdistance in ncdm_read_ncfl. This is a bug.'
          call fexit()
        end if
      else
        if (size(cludata,dim=1).ne.(2*ncdm_nfeats)) then
          write(ilog,*)
          write(ilog,*) 'Fatal. Got size(cludata,dim=1).ne.(2*ncdm_nfeats) with cdis 3 in ncdm_read_ncfl. This is a bug.'
          call fexit()
        end if
      end if
      if ((size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1))) then !only features in cludata
        write(ilog,*)
        write(ilog,*) 'Fatal. (size(cludata,dim=1)).ne.(size(ncdm_featvals_includata,dim=1)) with no dyn. weigths in &
 &ncdm_read_ncfl. This is a bug.'
        call fexit()
      end if
    end if
  end if
  if (size(ncdm_whichframes,dim=1).ne.ncdm_nframes) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(ncdm_whichframes,dim=1).ne.ncdm_nframes in ncdm_read_ncfl. This is a bug.'
    call fexit()
  end if
  if (size(cludata,dim=2).ne.ncdm_nframes) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(cludata,dim=2).ne.ncdm_nframes in ncdm_read_ncfl. This is a bug.'
    call fexit()
  end if
  if (size(ncdm_whichfeats,dim=1).ne.ncdm_nfeats) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Got size(ncdm_whichfeats,dim=1).ne.ncdm_nfeats in ncdm_read_ncfl. This is a bug.'
    call fexit()
  end if
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (size(ncdm_dynwghts_includata,dim=1).ne.ncdm_nfeats) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(ncdm_dynwghts_includata,dim=1).ne.ncdm_nfeats in ncdm_read_ncfl. This is a bug.'
      call fexit()
    end if
  end if
  if (cdis_crit.ne.3.AND.cdis_crit.ne.4) then
    if (size(ncdm_featvals_includata,dim=1).ne.ncdm_nfeats) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(ncdm_featvals_includata,dim=1).ne.ncdm_nfeats in ncdm_read_ncfl. This is a bug.'
      call fexit()
    end if
  else if (cdis_crit.eq.3.OR.cdis_crit.eq.4) then
    if (size(ncdm_featvals_includata,dim=1).ne.(2*ncdm_nfeats)) then !afalse
      write(ilog,*)
      write(ilog,*) 'Fatal. Got size(ncdm_featvals_includata,dim=1).ne.(2*ncdm_nfeats) in ncdm_read_ncfl. This is a bug.'
      call fexit()
    end if
  end if
!
end subroutine ncdm_read_ncfl
!
!---------------------------------------------------------------------------
!
subroutine ncdm_fill_auxcludata 
!
  use iounit
  use system, ONLY: nsim
  use clusters, ONLY: cdis_crit,calcsz,clstsz,cstored,cl_imvec,cludata,cstorecalc
  use pdb, ONLY: framecnt,framelst
  use ncdm
 
!
  implicit none
!
  integer, parameter :: aone = 1
!
  integer iu,freeunit
  integer iframe,ifeat,iweight,counter
!
  RTYPE, ALLOCATABLE:: tmparr(:,:)
!
! from distance criterion to allocation of features and weights. calcsz is dim 1 of cludata (include dyn weights)
  if (cdis_crit.eq.1) calcsz = ncdm_nfeats   ! dihedrals
  if (cdis_crit.eq.2) calcsz = 2*ncdm_nfeats ! dihedrals with time-dep. weights
  if (cdis_crit.eq.3) calcsz = 2*ncdm_nfeats 
  if (cdis_crit.eq.4) calcsz = 3*ncdm_nfeats 
  if (cdis_crit.eq.7) calcsz = ncdm_nfeats 
  if (cdis_crit.eq.8) calcsz = ncdm_nfeats    ! DRMS with static weights -----> in cl_imvec
  if (cdis_crit.eq.9) calcsz = 2*ncdm_nfeats  ! DRMS with time-dep. weights
!
! from distance criterion to allocation of features only. clstsz is how many features are in cludata (excludes dyn weights)
  if (cdis_crit.ne.3.AND.cdis_crit.ne.4) then
    clstsz = ncdm_nfeats    
  else 
    clstsz = 2*ncdm_nfeats 
  end if
!
! checks
  if ( (cdis_crit.eq.8).AND.(allocated(cl_imvec).EQV..false.) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten cl_imvec not allocated with cdistance 8 in ncdm_fill_auxcludata. This is a bug.'
    call fexit()
  end if
  if ( (cdis_crit.eq.8).AND.(ncdm_arethere_stwghts.EQV..false.) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten ncdm_arethere_stwghts false with cdistance 8 in ncdm_fill_auxcludata. This is a bug.'
    call fexit()
  end if
  if ( (cdis_crit.ne.8).AND.(allocated(cl_imvec).EQV..true.) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten cl_imvec allocated with cdistance not 8 in ncdm_fill_auxcludata. This is a bug.'
    call fexit()
  end if
  if ( ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9)).AND.(ncdm_arethere_dywghts.EQV..false.) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten ncdm_arethere_dywghts false with cdistance 2, 4 or 9 in ncdm_fill_auxcludata. This is a bug.'
    call fexit() 
  end if
  if ( ((cdis_crit.eq.2).OR.(cdis_crit.eq.9)).AND.(size(cludata,dim=1).ne.(2*ncdm_nfeats)) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten mismatch on size of dim 1 of cludata with cdistance 2 or 9 in ncdm_fill_auxcludata. &
 &This is a bug.'
    call fexit()
  else if ( ((cdis_crit.eq.1).OR.(cdis_crit.eq.7).OR.(cdis_crit.eq.8)).AND.(size(cludata,dim=1).ne.(ncdm_nfeats)) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten mismatch on size of dim 1 of cludata with cdistance 1, 7 or 8 in ncdm_fill_auxcludata. &
 &This is a bug.',ncdm_nfeats,size(cludata,dim=1)
    call fexit()
  end if
  if ( size(cludata,dim=2).ne.(ncdm_nframes) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten mismatch on size of dim 2 of cludata in ncdm_fill_auxcludata. This is a bug.'
    call fexit()
  end if
!
! other allocations and initializations
  cstored = ncdm_nframes
  nsim = ncdm_nframes
  framecnt = ncdm_nframes
  cstorecalc = 1
!
  if (allocated(framelst).EQV..true.) then
    deallocate(framelst)
  end if
  allocate(framelst(size(ncdm_whichframes)))
  framelst = ncdm_whichframes
!
end subroutine ncdm_fill_auxcludata
!
!---------------------------------------------------------------------------
!
!
subroutine ncdm_scale_periodicity_andsincos
!
  use iounit
  use clusters, ONLY: cludata,cdis_crit
  use math, ONLY: PI
  use ncdm
!
  implicit none
!
  integer iframe,ifeat
  integer iu,freeunit
!
  logical dousrscale
!
  RTYPE offset,tmpr,rng
!
  RTYPE, parameter:: leftbound = -180.0
  RTYPE, parameter:: rightbound = 180.0
!
! checks
  if (ncdm_arethere_dywghts.EQV..false.) then
    if (cdis_crit.ne.1.AND.cdis_crit.ne.3) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Wrong combination of parameters in ncdm_scale_periodicity (1). This is a bug.'
      call fexit()
    end if
  else
    if (cdis_crit.ne.2.AND.cdis_crit.ne.4) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Wrong combination of parameters in ncdm_scale_periodicity (2). This is a bug.'
      call fexit()
    end if
  end if  
!
! scale data to the user interval: only the features that are out of range
  dousrscale = .false. 
  if ( ( minval(cludata(ncdm_featvals_includata,:)).lt.ncdm_periodic_rng(1) ).OR.( &
 &       maxval(cludata(ncdm_featvals_includata,:)).gt.ncdm_periodic_rng(2)      ) ) then
    write(ilog,*) 'Warning. The provided input data exceed the reference periodic interval. This is corrected for analysis.'
    dousrscale = .true.
  end if
  if (dousrscale.EQV..true.) then
    iu = freeunit()
    rng = ncdm_periodic_rng(2) - ncdm_periodic_rng(1)
    do iframe=1,ncdm_nframes
      do ifeat=1,size(ncdm_featvals_includata,dim=1) !due to cdist 3 or 4
        if ((cludata(ncdm_featvals_includata(ifeat),iframe).lt.ncdm_periodic_rng(1)).OR. &
 &          (cludata(ncdm_featvals_includata(ifeat),iframe).gt.ncdm_periodic_rng(2))) then
          tmpr = cludata(ncdm_featvals_includata(ifeat),iframe)
          cludata(ncdm_featvals_includata(ifeat),iframe) = tmpr - (ceiling( (tmpr - ncdm_periodic_rng(2)) / rng ) * rng)
        end if  
      end do
    end do
  end if
!
! procedure to map the data to [-180:180]
  offset = (ncdm_periodic_rng(2) + ncdm_periodic_rng(1))/2. !center of periodicity interval
  ncdm_periodic_rng(2) = ncdm_periodic_rng(2) - offset      !re-centered margins of periodic interval   
  ncdm_periodic_rng(1) = ncdm_periodic_rng(1) - offset
  do iframe=1,ncdm_nframes  !mapping to [-180:180]: all features of all snapshosts
    do ifeat=1,size(ncdm_featvals_includata,dim=1) !due to cdist 3 or 4
      tmpr = cludata(ncdm_featvals_includata(ifeat),iframe)
      cludata(ncdm_featvals_includata(ifeat),iframe) = ( (tmpr - offset) * rightbound )/ncdm_periodic_rng(2)
    end do
  end do
!
! final check
  do iframe=1,ncdm_nframes
    do ifeat=1,size(ncdm_featvals_includata,dim=1) !due to cdist 3 or 4
      if ((cludata(ncdm_featvals_includata(ifeat),iframe).lt.leftbound).OR. &
 &        (cludata(ncdm_featvals_includata(ifeat),iframe).gt.rightbound)) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Periodicity bound exceeded. This is a bug.'
        call fexit()
      end if  
    end do
  end do
!
! deal with sin and cosine
  if (cdis_crit.eq.3.OR.cdis_crit.eq.4) then
    do iframe=1,ncdm_nframes !,2
      do ifeat=1,size(ncdm_featvals_includata,dim=1),2 
        cludata(ncdm_featvals_includata(ifeat),iframe) = dsin(cludata(ncdm_featvals_includata(ifeat),iframe)*PI/rightbound)
      end do
      do ifeat=2,size(ncdm_featvals_includata,dim=1),2 
        cludata(ncdm_featvals_includata(ifeat),iframe) = dcos(cludata(ncdm_featvals_includata(ifeat),iframe)*PI/rightbound)
      end do
    end do
  end if
!
end subroutine ncdm_scale_periodicity_andsincos
!---------------------------------------------------------------------------
!
#endif
!
!---------------------------------------------------------------------------
