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
! MAIN AUTHOR:   Andreas Vitalis                                           !
! CONTRIBUTIONS: Pippo Pazzo                                               !
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
subroutine ncdm_sanity_checks_1
!
  use iounit
  use clusters 
  use params, ONLY: be_unsafe
  use ncdm
!
  implicit none
!
  integer t1,t2                       !to determine string limits
!
  logical exists
!
! I do not even welcome you if you misbehave so bad
  if ((ncdm_doasconv.EQV..true.).AND.(ncdm_donc.EQV..true.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Requesting conversion of an ascii-file and analysis on an input net-cdf file is &
 &explicitly disallowed for many reasons. If the goal is to convert and analyse the same file, just set &
 &NCDM_ANONAS to 1. If the goal is to mine the net-cdf input file, just comment keyword NCDM_ASFILE.'
    call fexit()
  end if
!
  write(ilog,*)
  write(ilog,*) '-- CAMPARI IS USED IN NETCDF-BASED DATA MINING ONLY    --'
  write(ilog,*) '-- ALL SETTINGS THAT ARE NOT RELEVANT ARE IGNORED      --'
  write(ilog,*)
!
  if (ncdm_doasconv.EQV..true.) then
    write(ilog,*) '-- AN EXTERNAL DATA FILE WILL BE CONVERTED   --'    
    write(ilog,*) '-- TO NET-CDF FORMAT                         --'    
    if (ncdm_isdmonas.EQV..true.) then
      write(ilog,*) '-- AND ANALYZED WITH SPECIFIC ROUTINES       --'
    end if
    write(ilog,*)
  end if
  if (ncdm_donc.EQV..true.) then
    write(ilog,*) '-- AN EXTERNAL NETCDF FILE WILL BE READ AND  --'    
    write(ilog,*) '-- POSSIBLY MINED WITH SPECIFIC ROUTINES     --'    
    write(ilog,*)
  end if
!
  if (be_unsafe.EQV..true.) then
    write(ilog,*) 'Warning. No check on input ascii file will be performed. If there is not the &
 &exact number of features expected per line, wrong results may be generated without a crash.'
  end if
!
  if ((ncdm_isdmonas.EQV..true.).AND.(ncdm_doasconv.EQV..false.)) then
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
  if (cstorecalc.gt.1) then
    write(ilog,*) 'Warning. CCOLLECT is greater than 1. This will entail subsampling the input data to the &
 &specified frequency, regardless their input origin (ascii or net-cdf) for analysis.'
    if (ncdm_isthere_framesfl.EQV..true.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. The concurrent support of CCOLLECT and NCDM_FRAMESFILE is explicitly not offered. &
 &Plese, either set CCOLLECT to 1 or do not use a FRAMESFILE.' 
      call fexit()
    end if
  end if
!
  if (ncdm_doasconv.EQV..true.) then 
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
    if (ncdm_isdmonas.EQV..true.) then
      if (3.gt.int(floor(1.*ncdm_nframes_ini/cstorecalc))) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Due to the settings of CCOLLECT and NCDM_NRFRMS there will be less than three &
 &frames to analyze and no data mining can be performed despite the request.'
        call fexit()
      end if
    end if
  end if !if (ncdm_doasconv.EQV..true.) then
!
! ncdm_donc
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
! clustering
!
  if ((cdis_crit.ne.1).AND.(cdis_crit.ne.2).AND.((cdis_crit.lt.7).OR.(cdis_crit.gt.9))) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Unsupported distance criterion for Net-CDF data mining (see CDISTANCE). Gotten: ',cdis_crit
    call fexit()
  end if
!
! cut profiles
!
  if ((ccfepmode.gt.0).AND.(cequil.le.0)) then
    write(ilog,*) 'Warning. It can produce misleading results to compute cut-based pseudo free energy profiles (setting for &
 &CMSMCFEP) without equilibrating the underlying Markov state model first (choice for CREWEIGHT).'
  end if
!
! data preprocessing
!
  if ((cdis_crit.ne.1).AND.(cdis_crit.ne.2).AND.((cprepmode.eq.2).OR.(cprepmode.eq.5))) then
    write(ilog,*) 'Warning. The normalization by sample variance (setting for CPREPMODE) can lead to nonsensical &
 &results or program crashes if one or more dimensions have negligible variance.'
  end if
!
! weigths for clustering
!
  if ((cchangeweights.gt.0).AND.(cdis_crit.ne.2).AND.(cdis_crit.ne.8).AND.(cdis_crit.ne.9)) then
    write(ilog,*) 'Warning. Requesting a replacement of weights for structural clustering (setting for CMODWEIGHTS) &
 &is only possible if an eligible distance function is used (CDISTANCE has to be 2, 4, or 8-10). Disabled.'
    cchangeweights = 0
  end if
!
! cfile 
!
  if (ncdm_isthere_cfile.EQV..true.) then
    call strlims(ncdm_nm_cfile,t1,t2)
    inquire(file=ncdm_nm_cfile(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal, provided input cfile (NCDM_CFILE) ',ncdm_nm_cfile(t1:t2),' not found or corrupted.'
      ncdm_isthere_cfile = .false.
      call fexit()
    end if
    write(ilog,*) 'Warning. A CFILE has been provided. This will entail using only a subset of the features regardless &
 &the input origin of the data (ascii or net-cdf) for analysis.'
  end if
!
! frames file
!
  if (ncdm_isthere_framesfl.EQV..true.) then
    call strlims(ncdm_nm_framesfl,t1,t2)
    inquire(file=ncdm_nm_framesfl(t1:t2),exist=exists)
    if (exists.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal, provided input cfile (NCDM_FRAMESFILE) ',ncdm_nm_framesfl(t1:t2),' not found or corrupted.'
      ncdm_isthere_framesfl = .false.
      call fexit()
    end if
    if (cstorecalc.ne.1) then 
      write(ilog,*)
      write(ilog,*) 'Fatal. The concurrent support of NCDM_FRAMESFILE and CCOLLECT is explicitly not offered. &
 &Plese, either set CCOLLECT to 1 or do not use a FRAMESFILE.' 
      call fexit()
    end if
    write(ilog,*) 'Warning. A FRAMESFILE has been provided. This will entail using only a subset of the frames regardless &
 &the input origin of the data (ascii or net-cdf) for analysis.'
  end if
!
! Sanity check against cfep, synthetic trajs., ergodicity, spectral analysis, pfold
!
  if ((synmode.ne.0).AND.((cmode.eq.4).AND.(cprogindex.eq.1))) then
    write(ilog,*) 'Warning. Disabling the generation of synthetic trajectories because no clusters are generated & 
 &when the exact progress index method is in use (CPROGINDMODE).'
    synmode = 0
  end if
  if ((synmode.ne.0).AND.((caddlkmd.ne.1).AND.(caddlkmd.ne.3).AND.(caddlkmd.ne.4))) then
    write(ilog,*) 'Warning. Generation of synthetic trajectories is performed without imposition of reversibility.'
  end if
  if(((synmode.eq.1).or.(synmode.eq.2)).AND.(inissnap_usrslct.eq.endssnap_usrslct)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Initial and target snapshots for the generation of synth. trajs. &
  &are the same.'
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
!    if (arnblks.ne.(numeig + 2)) then
!      write(ilog,*) 'Warning. The number of blocks used in the blocked Arnoldi method with Chebychev acceleration of &
! &the starting vectors differs from the default choice.'
!    end if
!    if (arnstps.ne.ceiling((8.*numeig)/arnblks)) then
!      write(ilog,*) 'Warning. The number of steps used in the blocked Arnoldi method with Chebychev acceleration of &
! &the starting vectors differs from the default choice.'
!    end if
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
!
  if (dopfold.EQV..true.) then
#ifdef LINK_HSL
    if (clagt_msm.gt.1) then
      write(ilog,*) 'Warning. Specifying a non-unit MSM lagtime is somewhat untested. Carefully inspect your results.'
    end if
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
#else
    write(ilog,*)
    write(ilog,*) 'Fatal. Attempting to use an unsupported feature. For committor probabilities it is required &
 &to link the code to the HSL library (see installation instruction).'
    call fexit()
#endif
  end if
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
#ifdef ENABLE_THREADS
  use threads, ONLY: thrdat
#endif

!
  implicit none
!
  integer t1,t2
  integer, parameter:: sectenyrs=315576000
!
  5   format(a,1x,a)
  21  format(a,2x,i8)
  50  format(a,7x,f8.4)
  
!
  write(ilog,*)
  write(ilog,*) '-- SETTINGS RELEVANT TO NET-CDF DATA MINING (1st PART) --'
  if (ncdm_doasconv.EQV..true.) then
    write(ilog,*) 'Input ascii file (NCDM_ASFILE)         : ',trim(ncdm_fn_as)
    if (ncdm_isdmonas.EQV..true.) then
      write(ilog,*) 'All analyses will refer to the above-specified ascii file (NCDM_ANONAS).'
    end if 
    if (ncdm_is_periodic.EQV..true.) then
      write(ilog,*) 'Periodicity range (NCDM_PRDCRNG):',ncdm_periodic_rng(1),ncdm_periodic_rng(2)
    end if
  end if
  if (ncdm_donc.EQV..true.) then
    write(ilog,*) 'Input netcdf file (NCDM_NCFILE)        : ',trim(ncdm_fn_r)
    write(ilog,*) 'Assumed slowest dimension in input net-cdf file: frames'
  end if
  if ((ncdm_donc.EQV..true.).AND.(ncdm_docheck_r.EQV..true.)) then
    write(ilog,*) 'Output net-cdf file will be written (NCDM_DOINPTCHECK).'
  end if
  if (ncdm_doasconv.EQV..true.) then
    write(ilog,*) 'Number of features (NCDM_NRFEATS)        : ',ncdm_nfeats_ini
    write(ilog,*) 'Name of netcdf variable (NCDM_NMFV)      : ',trim(ncdm_nm_feats)
    write(ilog,*) 'Number of frames (NCDM_NRFRAMES)         :',ncdm_nframes_ini
    write(ilog,*) 'For conversion, all frames will be retained.'
  end if
!
! framesfile, this matters regardless
  if (ncdm_isthere_framesfl.EQV..true.) then
    call strlims(ncdm_nm_framesfl,t1,t2)
    write(ilog,*) 'Input framesfile (NCDM_FRAMESFILE) : ',ncdm_nm_framesfl(t1:t2)
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
  if ((ncdm_isdmonas.EQV..true.).OR.(ncdm_donc.EQV..true.)) then !data mining settings relevant regardless
    call strlims(basename,t1,t2)
    write(ilog,5) ' Basename for output   (BASENAME) : ',basename(t1:t2)
    write(ilog,*) 'Distance criterion index  (CDISTANCE) : ',cdis_crit
!
    if ((cdis_crit.eq.1).OR.(cdis_crit.eq.2)) then 
      write(ilog,*) 'The specified distance criterion implies that the data are periodic.'
      if (ncdm_isdmonas.EQV..true.) then 
        if (ncdm_is_periodic.EQV..false.) then
          write(ilog,*) 'No further information have been specified with keyword NCDM_PRDCRNG. &
   &The assumed periodicity range is -180 to 180. Data exceeding these values will be rescaled within &
 &the reference interval.' !'provoke a failure.'
        end if
      else
        write(ilog,*) 'Periodicity will be inferred directly from input net-cdf file from attribute periodic_range.'
        write(ilog,*) 'In case the attribute periodic_range is missing, the assumed range is -180 to 180. Data outside &
 &this range will be rescaled within the reference interval.' !'provoke a failure.'
      end if
    end if
    if (ncdm_donc.EQV..true.) then !just written in the other case
      write(ilog,*) 'Struct. clust. collection interval (CCOLLECT) : ',cstorecalc
    end if
    if (cstorecalc.eq.1) then
      write(ilog,*) 'All frames will be retained for analysis.'
    else if (cstorecalc.gt.1) then
      if (ncdm_donc.EQV..true.) then
        write(ilog,*) 'Due to the settings for CCOLLECT input frames will be subsampled at an interval of:',cstorecalc, &
 &' for analysis.'
      else if (ncdm_isdmonas.EQV..true.) then
        write(ilog,*) 'The analyzed frames are a subset of the converted ones according to the value for CCOLLECT.'
      else
        write(ilog,*) 'No analysis is requested. Conversion discards the value for CCOLLECT, all frames retained.'
      end if
    end if
!  
!   preprocessing
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
!
!   clustering 
    if (((cdis_crit.eq.2).OR.((cdis_crit.ge.8).AND.(cdis_crit.lt.10))).AND.(cchangeweights.gt.0)) then
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
    if ((pcamode.eq.3).OR.(pcamode.eq.5)) then
      if (reduced_dim_clustering.gt.0) then
        write(ilog,*) 'Clustering algorithm will be run on a subset of transformed components.'
        write(ilog,*) 'No. of components to use  : ',reduced_dim_clustering
      end if
    end if
    if ((pcamode.eq.4).OR.(pcamode.eq.5)) then
      write(ilog,*) 'Lag time for ACF (snaps.) : ',cwacftau
    end if
    if (ntbrks.gt.0) then
      call strlims(tbrkfilen,t1,t2)
      write(ilog,*) 'File with traject. breaks : ',tbrkfilen(t1:t2)
    end if
    if (ntlnks.gt.0) then
      call strlims(tlnkfilen,t1,t2)
      write(ilog,*) 'File with add. traj. links: ',tlnkfilen(t1:t2)
    end if
    if (brklnk_report.gt.0) then
      write(ilog,*) 'In case a break, link and/or trace file(s) have been provided, additional information on connections &
   &manipulation will be reported to standard output.'
    end if
    if (caddlkmd.gt.0) then
      write(ilog,*) 'Mode for adding links     : ',caddlkmd
      write(ilog,50) ' Thresh. for added links   : ',caddlinks
      write(ilog,50) ' Weight of added links     : ',caddlinkwt
    else
      write(ilog,*) 'No links (edges) are added to network (graph).'
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
!
!   equilibration
    if (cequil.gt.0) then
      write(ilog,*) 'This implies independent treatment of strongly connected components.'
      write(ilog,50) ' Buffer for snapshot wts.  : ',cequilbuf
! ------------------------ TO BE RETESTED -----------------------------------------------------
      if (maxtime_eqmsm.lt.sectenyrs) then
        write(ilog,*) 'Time threshold to exit MSM equilibration :',maxtime_eqmsm,' [s].'
      else 
        write(ilog,*) 'MSM equilibration can potentially continue for more than ten years.'
      end if
! ---------------------------------------------------------------------------------------------
    else
      write(ilog,*) 'Complete graph is analyzed independently of connectedness. This may cause errors.'
    end if
    if (ccfepmode.eq.0) then
      write(ilog,*) 'No computation of cut-based pseudo free energy profiles performed.'
    else if (ccfepmode.eq.1) then
      write(ilog,*) 'Requested computation of cut-based pseudo free energy profile using iterative MFPT values as order parameter.'
      write(ilog,*) 'Ref. node for profile     : ',cprogindstart
! ------------------------ TO BE RETESTED -----------------------------------------------------
      if (maxtime_eqmsm.lt.sectenyrs) then
        write(ilog,*) 'Time threshold to exit MSM equilibration :',maxtime_eqmsm,' [s].'
      else
        write(ilog,*) 'MSM equilibration can potentially continue for more than ten years.'
      end if
! ---------------------------------------------------------------------------------------------
    end if
!
!   transition matrix(ces)
    if ((tmat_report.EQV..true.).AND.((eigvalmd.ne.0).OR.(dopfold.EQV..true.).OR.(synmode.ne.0))) then
      write(ilog,*) 'Report file for transition matrix(ces) will be written.'
    end if
!
!   synthetic trajectories
    if (synmode.ne.0) then
      write(ilog,*) 'Synthetic trajectory mode           : ',synmode
      if (synmode.eq.1) then
        write(ilog,*) 'Target numb. of synth. traj.        : ',nstraj
        if (inissnap_usrslct.gt.0) then
          write(ilog,*) 'Start snaphot of synth. traj.       : ',inissnap_usrslct,' of whole-length trajectory (may be remapped).'
          write(ilog,*) 'Component for random walker         : The one that hosts snapshot ',inissnap_usrslct, &
 &' of whole-length trajectory (may be remapped).'
        else
          write(ilog,*) 'Start snaphot of synth. traj.       : centroid of largest cluster.'
          write(ilog,*) 'Component for random walker         : The one that hosts the centroid &
 &of the largest cluster'
        end if
      else if (synmode.eq.2) then 
        if (inissnap_usrslct.gt.0) then
          write(ilog,*) 'Start snaphot of synth. traj.       : ',inissnap_usrslct,' of whole-length trajectory (may be remapped).'
          write(ilog,*) 'Component for random walker         : The one that hosts snapshot ',inissnap_usrslct, &
 &' of whole-length trajectory (may be remapped).'
        else
          write(ilog,*) 'Start snaphot of synth. traj.       : centroid of largest cluster.'
          write(ilog,*) 'Component for random walker         : The one that hosts snapshot ',inissnap_usrslct, &
 &' of the whole-length trajectory (may be remapped).'
        end if
      else
        if (inissnap_usrslct.gt.0) then
          write(ilog,*) 'Component for random walker         : The one that hosts snapshot ',inissnap_usrslct, &
 &' of the whole-length trajectory (may be remapped).'
        else
          write(ilog,*) 'Component for random walker         : The one that will host the centroid &
 &of the largest cluster'
        end if
      end if 
      if (synmode.eq.1) then
        write(ilog,*) 'End (target) snaphot of synth. traj.: ',endssnap_usrslct,' of whole-length trajectory (may be remapped).'
        write(ilog,*) 'Max snaps. per synthetic trajectory : ',nssnap
        write(ilog,*) 'Max number of syn. traj.            : ',nstraj
      else if ((synmode.eq.2).or.(synmode.eq.3)) then
        write(ilog,*) 'End (target) cluster of synth. traj.: The one hit after',nssnap,' steps.'
        write(ilog,*) 'Snaps. per synthetic trajectory : ',nssnap
        write(ilog,*) 'Number of syn. traj.  : ',nstraj
      else
        write(ilog,*) 'Unrecognized synth. traj. mode: This is a bug.'
        call fexit()
      end if
      write(ilog,*) 'Lag time for sliding-window network-based analyses : ',clagt_msm  !sliding window only?
    else
      write(ilog,*) 'Generation of synthetic trajectories is not requested.'
    end if
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
subroutine ncdm_sanity_checks_2
!
  use iounit
  use clusters, ONLY: cdis_crit,pcamode,clstsz,reduced_dim_clustering
  use ncdm
!
  implicit none
!
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
  if (3.gt.ncdm_nframes) then
    write(ilog,*)
    write(ilog,*) 'Fatal. There are less than three frames to analyze and no data mining can be performed despite the request.'
    call fexit()
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(ncdm_arethere_dywghts.EQV..true.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Both static and dynamic weights in use. This is a bug.'
    call fexit()
  end if
  if ((ncdm_is_periodic.EQV..true.).AND.((cdis_crit.ne.1).AND.(cdis_crit.ne.2))) then
    write(ilog,*)
    write(ilog,*) 'Warning. Periodic indicator is set to .true. but the selected distance criterion (CDISTANCE) &
 &does not account for periodicity.'
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
subroutine ncdm_summary_2
!
  use iounit
  use clusters, ONLY: cdis_crit
  use ncdm
!
  implicit none
!
  write(ilog,*)
  write(ilog,*) '-- SETTINGS RELEVANT TO NET-CDF DATA MINING (2nd PART) --'
  write(ilog,*) 'Number of frames present in the input file     : ',ncdm_nframes_ini
  write(ilog,*) 'Number of features present in the input file   : ',ncdm_nfeats_ini
  write(ilog,*) 'Number of features actually used               : ',ncdm_nfeats
  write(ilog,*) 'Final number of frames actually used           : ',ncdm_nframes
  if (ncdm_is_periodic.EQV..true.) then
    if ((cdis_crit.eq.1).OR.(cdis_crit.eq.2)) then
      write(ilog,*) 'Assumed periodicity range: ',ncdm_periodic_rng(1),ncdm_periodic_rng(2)
    else
      write(ilog,*) 'Warning. If data are periodic, CDISTANCE should be set to either 1 or 2.'
      write(ilog,*) 'With current value of CDISTANCE, no periodicity correction is performed.'
    end if
  else
    write(ilog,*) 'No periodicity correction in use.'
  end if
  if (ncdm_arethere_stwghts.EQV..true.) then
    write(ilog,*) 'Static weights are in use.'
  end if
  if (ncdm_arethere_dywghts.EQV..true.) then
    write(ilog,*) 'Dynamic weights are in use.'
  end if
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
!
! find out actual number of frames and relevant maximum index (see ncdm_whichframes)
  iu = freeunit()
  call strlims(ncdm_nm_framesfl,t1,t2)
  open(unit=iu,file=ncdm_nm_framesfl(t1:t2),status='old')  !open frames file to get to know how many frames will be used
  do while(.true.)
    read(iu,*,iostat=iomess) dummyi
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.ne.0) then 
      write(ilog,*) 'Fatal. File I/O error while reading frames file.'
      call fexit() 
    end if
    maxframe = max(maxframe,dummyi)              !highest frame index that is in the frames file
    if (dummyi.gt.0) then
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
  allocate(ncdm_whichframes(nframes_framesfl))
  ncdm_whichframes(:) = 0
!
! read and assign eligible ones
  counter = 0
  iu = freeunit()
  call strlims(ncdm_nm_framesfl,t1,t2)
  open(unit=iu,file=ncdm_nm_framesfl(t1:t2),status='old')  !open frames file to get to know how many frames will be used
  do while(.true.)
    read(iu,*,iostat=iomess) dummyi
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.ne.0) then 
      write(ilog,*) 'Fatal. File I/O error while reading frames file.'
      call fexit() 
    end if
    if (dummyi.gt.0) then
      counter = counter + 1
      ncdm_whichframes(counter) = dummyi
    end if
  end do 
  close(unit=iu)
!
! the mask
  if (ncdm_doasconv.EQV..true.) then  !this is important
    if (allocated(ncdm_framesmask).EQV..true.) deallocate(ncdm_framesmask)
    allocate(ncdm_framesmask(ncdm_nframes_ini))
    ncdm_framesmask = .false.  !adjusted below on eligible frames
    ncdm_framesmask(ncdm_whichframes) = .true.
  else  !this is useless, but I do it so that I do not have to distinguish the two cases below
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
  if (ncdm_doasconv.EQV..true.) then
    if (maxval(ncdm_whichframes).gt.ncdm_nframes_ini) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Detected a frame indicator greater than the assumed number of frames (NCDM_NFRAMES) &
 &in input ascii file (NCDM_ASFILE).'
      call fexit()
    end if
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
  if ((ncdm_doasconv.EQV..true.).AND.((from_read_ncfl.EQV..false.))) then
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
  else if ((ncdm_donc.EQV..true.).AND.(from_read_ncfl.EQV..true.)) then
    if (ncdm_isthere_framesfl.EQV..true.) then  !it already set many things in ncdm_read_framesfl
      if (ncdm_nframes_ini.lt.maxval(ncdm_whichframes)) then  !first meaningful thing to do here is to check consistency
        write(ilog,*)
        write(ilog,*) 'Fatal. There are less snapshots in the input net-cdf file that what stated in the framesfile. (1)'
        call fexit()
      end if
      if (allocated(ncdm_framesmask)) deallocate(ncdm_framesmask)  !we have to do this as the previous assignment was incomplete
      allocate(ncdm_framesmask(ncdm_nframes_ini))  !now allocated to ncdm_nframes_ini rather than to maxframe
      ncdm_framesmask = .false.  
      ncdm_framesmask(ncdm_whichframes) = .true.   !set to true the ones corresponding to the eligible frames
    else  !so we are here for the first time with no framesfile
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
    write(ilog,*) 'Fatal. Inside ncdm_manage_whichframes with the wrong flags combination.'
    call fexit()
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
  integer, allocatable :: tmparr(:)
!
  atrue = .true.
  nfeats_cfile = 0  !number of features that are resulting from cfile
  maxfeat = 0
!
! find out actual number of features and relevant maximum index (see ncdm_whichfeats)
  iu = freeunit()
  call strlims(ncdm_nm_cfile,t1,t2)
  open(unit=iu,file=ncdm_nm_cfile(t1:t2),status='old')  !open cfile to get to know how many features will be used (if all ok)
  do while(.true.)
    read(iu,*,iostat=iomess) dummyi
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.ne.0) then 
      write(ilog,*) 'Fatal. File I/O error while reading cfile file.'
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
    if (dummyi.ne.maxfeat) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Gotten Juventus.'
      call fexit()
    end if  
    if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
    allocate(ncdm_whichfeats(nfeats_cfile))
    ncdm_whichfeats(nfeats_cfile) = maxfeat
    if (ncdm_doasconv.EQV..true.) then  !specific checks for this case are required
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
    return
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
  if (counter.eq.1) then !only one feature -> only one good value in ncdm_whichfeats(1)
    tmparr(1) = ncdm_whichfeats(1)
    if (allocated(ncdm_whichfeats).EQV..true.) deallocate(ncdm_whichfeats)
    allocate(ncdm_whichfeats(counter))
    ncdm_whichfeats(counter) = tmparr(1)
    ncdm_nfeats = counter
    if (arethere_wrongs.EQV..false.) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Gotten unexisting duplicates with counter 1. This is a bug as we should have returned earlier.'
      call fexit()
    else
      write(ilog,*) 'Warning. The input cfile (NCDM_CFILE) contained wrong feature indexes, which have been removed. &
 &It remains only a feature which is not wrong.'
    end if
    if (ncdm_doasconv.EQV..true.) then  !specific checks for this case are required
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
  nfeats_cfile = counter      !number of not duplicated features
!
! checks before final assignment
  if (ncdm_doasconv.EQV..true.) then
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
  if (ncdm_doasconv.EQV..true.) then  !this is important
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
      else if (ncdm_whichfeats(ifeat).le.0) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Gotten negative feature value at the end of the cleaning procedure. This is a bug or a &
  &mismatched data type in input cfile. (1)'
        call fexit()
      end if
    end do
  else if (ncdm_whichfeats(1).le.0) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten negative feature value at the end of the cleaning procedure. This is a bug or a &
  &mismatched data type in input cfile. (2)'
    call fexit()
  end if
  if (nfeats_cfile.ne.size(ncdm_whichfeats)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten nfeats_cfile.ne.size(ncdm_whichfeats) at the end of ncdm_read_cfile. This is a bug.'
    call fexit()
  end if
  if (nfeats_cfile.gt.maxval(ncdm_whichfeats)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Gotten more features than maxval of indexes. This is a bug.'
    call fexit()
  end if
  if (ncdm_doasconv.EQV..true.) then
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
  if ((ncdm_doasconv.EQV..true.).AND.((from_read_ncfl.EQV..false.))) then
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
  else if ((ncdm_donc.EQV..true.).AND.(from_read_ncfl.EQV..true.)) then
    if (ncdm_isthere_cfile.EQV..true.) then !it already set many things in ncdm_read_cfile
      if (ncdm_nfeats_ini.lt.maxval(ncdm_whichfeats)) then  !only meaningful thing to do here is to check consistency
        write(ilog,*)
        write(ilog,*) 'Fatal. In the input net-cdf file, there are less features that what stated in the cfile. (1)'
        call fexit()
      end if
      if (allocated(ncdm_featsmask).EQV..true.) deallocate(ncdm_featsmask)
      allocate(ncdm_featsmask(ncdm_nfeats_ini))
      ncdm_featsmask = .false.
      ncdm_featsmask(ncdm_whichfeats) = .true.  
    else  !so we are here for the first time without having been previously in ncdm_read_cfile
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
    write(ilog,*) 'Fatal. Inside ncdm_manage_whichfeats with the wrong flags combination.'
    call fexit()
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
  integer, allocatable :: tmpvec(:)
  logical, INTENT(IN)  :: from_read_ncfl
!
  if ((ncdm_doasconv.EQV..true.).AND.(from_read_ncfl.EQV..false.)) then
    nframes = int(floor(1.*ncdm_nframes/cstorecalc))  !ncdm_nframes already set right by ncdm_manage_whichframes or framesfile
    allocate(tmpvec(nframes))                         !will be new ncdm_whichframes
    counter = 0
    do iframe=1,ncdm_nframes  !ncdm_nframes: either ncdm_nframes_ini or the number of frames in framesfile
      if (mod(iframe,cstorecalc).eq.0) then           !eligible
        counter = counter + 1                         !effective snapshot
        tmpvec(counter) = ncdm_whichframes(iframe)    !original numbering 
      end if
    end do
    if (nframes.ne.counter) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got a mismatch between counter and nframes in ncdm_manage_ccollect. This is a bug (1).'
      call fexit()
    end if
    deallocate(ncdm_whichframes)
    allocate(ncdm_whichframes(nframes))
    do iframe=1,nframes
      ncdm_whichframes(iframe) = tmpvec(iframe)
    end do
    deallocate(tmpvec)
    ncdm_nframes = nframes  
    ncdm_framesmask = .false.
    ncdm_framesmask(ncdm_whichframes) = .true.
  else if ((ncdm_donc.EQV..true.).AND.(from_read_ncfl.EQV..true.)) then !ncdm_nframes is set by framesfl or by ncdm_manage_nframes
    nframes = int(floor(1.*ncdm_nframes/cstorecalc))
    allocate(tmpvec(nframes))                         !will be new ncdm_whichframes
    counter = 0
    do iframe=1,ncdm_nframes
      if (mod(iframe,cstorecalc).eq.0) then
        counter = counter + 1                         !effective snapshot
        tmpvec(counter) = ncdm_whichframes(iframe)    
      end if
    end do
    if (nframes.ne.counter) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got a mismatch between counter and nframes in ncdm_manage_ccollect. This is a bug (2)'
      call fexit()
    end if
    deallocate(ncdm_whichframes)
    allocate(ncdm_whichframes(nframes))
    do iframe=1,nframes
      ncdm_whichframes(iframe) = tmpvec(iframe)
    end do
    deallocate(tmpvec)
    ncdm_nframes = nframes  !finally accounting for possible ccollect subsampling
    ncdm_framesmask = .false.
    ncdm_framesmask(ncdm_whichframes) = .true.
  else if ((ncdm_donc.EQV..false.)) then
    write(ilog,*)   
    write(ilog,*) 'Fatal. Inside ncdm_manage_ccollect with the wrong flags combination.'
    call fexit()
  end if
!
end subroutine ncdm_manage_ccollect
!
!---------------------------------------------------------------------------
!
subroutine ncdm_convert_ascii !from input ascii to net cdf by calling ncdm_write_ncfl
!
  use iounit
  use clusters, ONLY: cdis_crit,cludata,cstorecalc,cl_imvec
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
  use params, ONLY: be_unsafe
  use ncdm
!
  implicit none
!
  integer freeunit,iu,iframe,counter,ifeat,featcounter
  integer jframe,endframe
  integer iomess,ipoint,howmany,ipointx3         !managing error in reading and string to reals
!
  integer, allocatable :: featvals_inarr(:)
  integer, allocatable :: dynwghts_inarr(:)
!
  integer, parameter :: charlen = 999999
!
  RTYPE dummyr
! 
  logical changed,goon,afalse           !goon controls what to do while reading
  logical write_weights(2)              !whether or not should be written to file. (1) is static, (2) dynamic
!
  character(MAXSTRLEN) frmt
  character (len=:), allocatable :: lnstrng
!
  RTYPE, allocatable :: data_on_row(:)  !temporary vector to deal with reading only specific rows of input file
  RTYPE, allocatable :: tmpvec(:)       !helper vector for weight management
  RTYPE, allocatable :: tmparr(:,:)     !temporary helper array to read data and store them when needed 
!
  iomess = 0
  goon = .false.
  afalse = .false.
  write_weights(:) = .true.  !if there are weights, we assume they are good
!
  changed = .false.
  if (ncdm_docheck_r.EQV..true.) then
    ncdm_docheck_r = .false.
    changed = .true.
  end if
!
! allocate the helper arrays 
  if (allocated(data_on_row).EQV..true.) deallocate(data_on_row) !to deal with possible feature reduction due to cfile
  allocate(data_on_row(ncdm_nfeats_ini))                         !the user-provided full rank
  if (allocated(tmparr).EQV..true.) deallocate(tmparr)
  allocate(tmparr(ncdm_nfeats_ini,ncdm_nframes_ini)) 
!
  iu = freeunit()
  open(unit=iu,file=ncdm_fn_as,status='old')  !opening ascii file
!
  if (be_unsafe.EQV..false.) then  !checking input string
    write(frmt,'("(a",I0,")")') charlen
    if (allocated(lnstrng).EQV..true.) deallocate(lnstrng)
    allocate(character(len=charlen) :: lnstrng)
    iframe = 1
    read(iu,frmt,iostat=iomess) lnstrng
    ipoint = 1
    call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
    if (howmany.ne.ncdm_nfeats_ini) then
      !write(ilog,*) lnstrng,data_on_row
      write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file or wrong format (e.g. scientific) (-1).',&
 &howmany,ncdm_nfeats_ini
      call fexit()
    else 
      call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse)
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
    tmparr(1:ncdm_nfeats_ini,iframe) = data_on_row  !values of features in upper half of array cludata
    deallocate(lnstrng)
    ipointx3 = 3*ipoint 
    allocate(character(len=ipointx3) :: lnstrng) !should be in the standard of most recent complilers
    write(frmt,'("(a",I0,")")') ipointx3 
    do iframe=2,ncdm_nframes_ini                 !reading the values of the features in temporary array
      read(iu,frmt,iostat=iomess) lnstrng
      ipoint = 1
      call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
      if (howmany.ne.ncdm_nfeats_ini) then
        write(ilog,*) 'Fatal. Wrong number of entries per line in input ascii file or wrong format (e.g. scientific) (1).'
        call fexit()
      else 
        call get_reals_from_str(lnstrng,1,dummyr,ipoint,howmany,afalse)
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
      tmparr(1:ncdm_nfeats_ini,iframe) = data_on_row  !values of features in upper half of array cludata
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
      call get_reals_from_str(lnstrng,ncdm_nfeats_ini,data_on_row,ipoint,howmany,afalse)
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
      if (allocated(tmpvec).EQV..true.) deallocate(tmpvec)  !we do not know yet whether static or dynamic or both
      allocate(tmpvec(1:ncdm_nfeats_ini))
      goon = .true.
    end if
!
    if (goon.EQV..true.) then  !there is at least one extra line 
      read(iu,frmt,iostat=iomess) lnstrng 
      if (iomess.eq.IOSTAT_END) then  !now we know there is only one extra line, i.e. there are only static weights
        ncdm_arethere_stwghts = .true.
        write(ilog,*) 'Warning. Assuming static weights from input file.'
        if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
        allocate(cl_imvec(1:ncdm_nfeats_ini))
        cl_imvec(:) = data_on_row
        goon = .false.                !this discrimiates below the copying of tmparr to cludata
      else
        ipoint = 1
        call get_reals_from_str(lnstrng,ncdm_nfeats_ini,tmpvec,ipoint,howmany,afalse)
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
        ncdm_arethere_stwghts = .false.  !two lines detected after the last of frames, we assume there are dynamic weights
        ncdm_arethere_dywghts = .true.
        write(ilog,*) 'Warning. Assuming dynamic weights from input file.'
        if (allocated(cludata).EQV..true.) deallocate(cludata)
        allocate(cludata(2*ncdm_nfeats_ini,ncdm_nframes_ini))
        if ((cdis_crit.ne.2).OR.(ncdm_isdmonas.EQV..false.)) then      
          cludata(1:ncdm_nfeats_ini,:) = tmparr(1:ncdm_nfeats_ini,:)   !features in upper half of array cludata
          deallocate(tmparr)
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,1) = data_on_row !dynamic weights of first snapshot
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,2) = tmpvec      !dynamic weights of second snapshot
        else
          do iframe=1,ncdm_nframes_ini
            featcounter = 0
            do ifeat=1,ncdm_nfeats_ini
              featcounter = featcounter + 1
              cludata(2*featcounter-1,iframe) = tmparr(ifeat,iframe)  !populating the odd lines of the columns with the features
            end do
          end do
          deallocate(tmparr)
          featcounter = 0
          do ifeat=1,ncdm_nfeats_ini
            featcounter = featcounter + 1
            cludata(2*featcounter,1) = data_on_row(ifeat)  !dynamic weights of first snapshot on even lines of first column
            cludata(2*featcounter,2) = tmpvec(ifeat)       !dynamic weights of second snapshot on even lines of second column
          end do
        end if
      end if
    end if
!
    if (goon.EQV..false.) then  !no dynamic weights present, just copy temparr to cludata
      if (ncdm_arethere_dywghts.EQV..true.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Gotten wrong combination of ncdm_arethere_dywghts and goon. This is a bug (1).'
        call fexit()
      end if
      if (((cdis_crit.ne.2).AND.(cdis_crit.ne.9)).OR.(ncdm_isdmonas.EQV..false.)) then !static weights cared not here if present
        if (allocated(cludata).EQV..true.) deallocate(cludata)
        allocate(cludata(ncdm_nfeats_ini,ncdm_nframes_ini))
        do iframe=1,ncdm_nframes_ini
          cludata(:,iframe) = tmparr(:,iframe)
        end do
        deallocate(tmparr)
      else if (ncdm_arethere_dywghts.EQV..false.) then
        ncdm_arethere_dywghts = .true.  !flat, but still...
        write_weights(2) = .false.      !do not write them in converted file though
        write(ilog,*) 'Warning. The current setting of CDISTANCE (to 2 or 9) requires the specification of dynamic weigths. &
 &Flat (uniformative) weights will be used instead to circumvent the lack of specification. They might become informative in & 
 &analysis depending on keyword CMODWEIGHTS and related ones.'
        if (allocated(cludata).EQV..true.) then
          write(ilog,*)
          write(ilog,*) 'Fatal. Gotten allocated cludata in the wrong moment. This is a bug. (2)'
          call fexit()
        end if
        allocate(cludata(2*ncdm_nfeats_ini,ncdm_nframes_ini))
        if (cdis_crit.eq.9) then      
          cludata(1:ncdm_nfeats_ini,:) = tmparr(1:ncdm_nfeats_ini,:)   !features in upper half of array cludata
          deallocate(tmparr)
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,:) = 1 !flat weights in lower half
        else if (cdis_crit.eq.2) then
          do iframe=1,ncdm_nframes_ini
            featcounter = 0
            do ifeat=1,ncdm_nfeats_ini
              featcounter = featcounter + 1
              cludata(2*featcounter-1,iframe) = tmparr(ifeat,iframe)  !populating the odd lines of the columns with the features
              cludata(2*featcounter,iframe) = 1                       !flat weights on even rows
            end do
          end do
        else
          write(ilog,*)
          write(ilog,*) 'Fatal. Wrong CDISTANCE. This is a bug.'
          call fexit()
        end if
      else
        write(ilog,*)
        write(ilog,*) 'Fatal. Undetermined combo. This is a bug. (1)'
        call fexit()
      end if
    else if (goon.EQV..true.) then  !we are then dealing with dynamic weights and cludata has been already allocated
      if (ncdm_arethere_dywghts.EQV..false.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong combination between goon and ncdm_arethere_dywghts. This is a bug. (1)'
        call fexit()
      end if
      if (ncdm_arethere_stwghts.EQV..true.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong combination between goon and ncdm_arethere_stwghts. This is a bug. (1)'
        call fexit()
      end if
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
        if ((cdis_crit.ne.2).OR.(ncdm_isdmonas.EQV..false.)) then  
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,iframe) = data_on_row  !weights in lower half of array cludata
        else
          featcounter = 0
          do ifeat=1,ncdm_nfeats_ini
            featcounter = featcounter + 1
            cludata(2*featcounter,iframe) = data_on_row(ifeat)  !weights on even lines of column iframe
          end do 
        end if
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
      else 
        goon = .false.
        write(ilog,*) 'Warning. No static weights detected in input file.'
      end if
      if (goon.EQV..true.) then  !just checking there is no an extra line...
        read(iu,*,iostat=iomess)
        if (iomess.ne.IOSTAT_END) then
          write(ilog,*) 
          write(ilog,*) 'Fatal. Unrecognized total number of lines in input ascii file. It appears that there are more lines &
   &than assuming both dynamic and static weights. (1)'
          call fexit()
        end if
      end if
    end if
!
  else  !be_unsafe
!
    do iframe=1,ncdm_nframes_ini                !reading the values of the features in temporary array
      read(iu,*,iostat=iomess) data_on_row
      if (iomess.eq.IOSTAT_END) then
        write(ilog,*)
        write(ilog,*) 'Fatal. Input ascii file appears to contain less lines than what specified with keyword NCDM_NFRAMES. (2)'
        call fexit()
      end if
      tmparr(1:ncdm_nfeats_ini,iframe) = data_on_row  !values of features in upper half of array cludata
    end do
!
    read(iu,*,iostat=iomess) data_on_row  !read one extra line, if not there it means there are no weights and do not go on reading
    if (iomess.eq.IOSTAT_END) then
      write(ilog,*) 'Warning. Input ascii file does not contain any feature weight.'
      ncdm_arethere_dywghts = .false.
      ncdm_arethere_stwghts = .false.
      goon = .false.                  !finished, just copy temparr to cludata
    else  
      if (allocated(tmpvec).EQV..true.) deallocate(tmpvec)  !we do not know yet whether static or dynamic or both
      allocate(tmpvec(1:ncdm_nfeats_ini))
      goon = .true.
    end if
!
    if (goon.EQV..true.) then
      read(iu,*,iostat=iomess) tmpvec
      if (iomess.eq.IOSTAT_END) then  !now we know there is only one extra line, i.e. there are only static weights
        ncdm_arethere_stwghts = .true.
        write(ilog,*) 'Warning. Assuming static weights from input file.'
        if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
        allocate(cl_imvec(1:ncdm_nfeats_ini))
        cl_imvec(:) = data_on_row
        goon = .false.                !this discrimiates below the copying of tmparr to cludata
      else
        ncdm_arethere_stwghts = .false.  !two lines detected after the last of frames, we assume there are dynamic weights
        ncdm_arethere_dywghts = .true.
        write(ilog,*) 'Warning. Assuming dynamic weights from input file.'
        if (allocated(cludata).EQV..true.) deallocate(cludata)
        allocate(cludata(2*ncdm_nfeats_ini,ncdm_nframes_ini))
        if ((cdis_crit.ne.2).OR.(ncdm_isdmonas.EQV..false.)) then      
          cludata(1:ncdm_nfeats_ini,:) = tmparr(1:ncdm_nfeats_ini,:)   !features in upper half of array cludata
          deallocate(tmparr)
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,1) = data_on_row !dynamic weights of first snapshot
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,2) = tmpvec      !dynamic weights of second snapshot
        else
          do iframe=1,ncdm_nframes_ini
            featcounter = 0
            do ifeat=1,ncdm_nfeats_ini
              featcounter = featcounter + 1
              cludata(2*featcounter-1,iframe) = tmparr(ifeat,iframe)  !populating the odd lines of the columns with the features
            end do
          end do
          deallocate(tmparr)
          featcounter = 0
          do ifeat=1,ncdm_nfeats_ini
            featcounter = featcounter + 1
            cludata(2*featcounter,1) = data_on_row(ifeat)  !dynamic weights of first snapshot on even lines of first column
            cludata(2*featcounter,2) = tmpvec(ifeat)       !dynamic weights of second snapshot on even lines of first column
          end do
        end if
      end if
    end if
!
    if (goon.EQV..false.) then  !no dynamic weights present, just copy temparr to cludata
      if (ncdm_arethere_dywghts.EQV..true.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Gotten wrong combination of ncdm_arethere_dywghts and goon. This is a bug (2).'
        call fexit()
      end if
      if (((cdis_crit.ne.2).AND.(cdis_crit.ne.9)).OR.(ncdm_isdmonas.EQV..false.)) then !static weights cared not here if present
        if (allocated(cludata).EQV..true.) deallocate(cludata)
        allocate(cludata(ncdm_nfeats_ini,ncdm_nframes_ini))
        do iframe=1,ncdm_nframes_ini
          cludata(:,iframe) = tmparr(:,iframe)
        end do
        deallocate(tmparr)
      else if (ncdm_arethere_dywghts.EQV..false.) then
        ncdm_arethere_dywghts = .true.  !flat, but still...
        write_weights(2) = .false.      !do not write them in converted file though
        write(ilog,*) 'Warning. The current setting of CDISTANCE (to 2 or 9) requires the specification of dynamic weigths. &
 &Flat (uniformative) weights will be used instead to circumvent the lack of specification. They might become informative in & 
 &analysis depending on keyword CMODWEIGHTS and related ones.'
        if (allocated(cludata).EQV..true.) then
          write(ilog,*)
          write(ilog,*) 'Fatal. Gotten allocated cludata in the wrong moment. This is a bug. (2)'
          call fexit()
        end if
        allocate(cludata(2*ncdm_nfeats_ini,ncdm_nframes_ini))
        if (cdis_crit.eq.9) then      
          cludata(1:ncdm_nfeats_ini,:) = tmparr(1:ncdm_nfeats_ini,:)   !features in upper half of array cludata
          deallocate(tmparr)
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,:) = 1 !flat weights in lower half
        else if (cdis_crit.eq.2) then
          do iframe=1,ncdm_nframes_ini
            featcounter = 0
            do ifeat=1,ncdm_nfeats_ini
              featcounter = featcounter + 1
              cludata(2*featcounter-1,iframe) = tmparr(ifeat,iframe)  !populating the odd lines of the columns with the features
              cludata(2*featcounter,iframe) = 1                       !flat weights on even rows
            end do
          end do
        else
          write(ilog,*)
          write(ilog,*) 'Fatal. Wrong CDISTANCE. This is a bug.'
          call fexit()
        end if
      else
        write(ilog,*)
        write(ilog,*) 'Fatal. Undetermined combo. This is a bug. (2)'
        call fexit()
      end if
    else if (goon.EQV..true.) then  !we are then dealing with dynamic weights and cludata has been already allocated
      if (ncdm_arethere_dywghts.EQV..false.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong combination between goon and ncdm_arethere_dywghts. This is a bug. (2)'
        call fexit()
      end if
      if (ncdm_arethere_stwghts.EQV..true.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. Wrong combination between goon and ncdm_arethere_stwghts. This is a bug. (2)'
        call fexit()
      end if
      do iframe=3,ncdm_nframes_ini  !already checked that are at least three frames. First two dynamic weights already read
        read(iu,*,iostat=iomess) data_on_row
        if (iomess.eq.IOSTAT_END) then
          write(ilog,*)
          write(ilog,*) 'Fatal. Input ascii file appears to contain less weights than frames. (2)'
          call fexit()
        end if
        if ((cdis_crit.ne.2).OR.(ncdm_isdmonas.EQV..false.)) then  
          cludata(ncdm_nfeats_ini+1:2*ncdm_nfeats_ini,iframe) = data_on_row  !weights in lower half of array cludata
        else
          featcounter = 0
          do ifeat=1,ncdm_nfeats_ini
            featcounter = featcounter + 1
            cludata(2*featcounter,iframe) = data_on_row(ifeat)  !weights on even lines of column iframe
          end do 
        end if
      end do
      read(iu,*,iostat=iomess) data_on_row  !looking for static weights at last line 
      if (iomess.ne.IOSTAT_END) then
        write(ilog,*) 'Warning. Assuming also static weights from input file.'
        ncdm_arethere_stwghts = .true.
        if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
        allocate(cl_imvec(ncdm_nfeats_ini))
        cl_imvec(:) = data_on_row
      else 
        goon= .false.
        write(ilog,*) 'Warning. No static weights detected in input file.'
      end if
      if (goon.EQV..true.) then  !just checking there is no an extra line...
        read(iu,*,iostat=iomess)
        if (iomess.ne.IOSTAT_END) then
          write(ilog,*) 
          write(ilog,*) 'Fatal. Unrecognized total number of lines in input ascii file. It appears that there are more lines &
   &than assuming both dynamic and static weights. (2)'
        end if
      end if
    end if
!
  end if  !be_unsafe
!
  close(unit=iu)
  if (allocated(tmparr).EQV..true.) deallocate(tmparr)
!
! checks
  if (size(cludata,dim=2).ne.ncdm_nframes_ini) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Mismatch on size of dim 2 of cludata in ncdm_convert_ascii. This is a bug.'
    call fexit()
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(ncdm_arethere_dywghts.EQV..true.)) then 
    if (size(cludata,dim=1).ne.(2*ncdm_nfeats_ini)) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Mismatch on size of dim 1 of cludata in ncdm_convert_ascii (1). This is a bug.'
      call fexit()
    end if
    write(ilog,*)
    write(ilog,*) 'Warning. Gotten both weights present in input ascii file. Only one type will be used (if any) in analysis, &
 &depending on the value of CDISTANCE.'
    if (allocated(cl_imvec).EQV..false.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got cl_imvec not allocated with arethere_stwghts .true. in ncdm_convert_ascii (1). This is a bug.'
      call fexit()
    else if(size(cl_imvec).ne.ncdm_nfeats_ini) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got cl_imvec with wrong size with arethere_stwghts .true. in ncdm_convert_ascii (1). This is a bug.'
      call fexit()
    end if
  else if ((ncdm_arethere_dywghts.EQV..true.).AND.(size(cludata,dim=1).ne.(2*ncdm_nfeats_ini))) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Mismatch on size of dim 1 of cludata in ncdm_convert_ascii (2). This is a bug.'
    call fexit()
  else if ((ncdm_arethere_dywghts.EQV..false.).AND.(size(cludata,dim=1).ne.(ncdm_nfeats_ini))) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Mismatch on size of dim 1 of cludata in ncdm_convert_ascii (3). This is a bug.'
    call fexit()
  else if (ncdm_arethere_stwghts.EQV..true.) then
    if (allocated(cl_imvec).EQV..false.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got cl_imvec not allocated with arethere_stwghts .true. in ncdm_convert_ascii (2). This is a bug.'
      call fexit()
    else if(size(cl_imvec).ne.ncdm_nfeats_ini) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got cl_imvec with wrong size with arethere_stwghts .true. in ncdm_convert_ascii (2). This is a bug.'
      call fexit()
    end if
  end if
  if ((cdis_crit.eq.8).AND.(ncdm_arethere_stwghts.EQV..false.).AND.(ncdm_isdmonas.EQV..true.)) then !done on the fly for dyn. wghts
    write(ilog,*)
    write(ilog,*) 'Warning. The current setting of CDISTANCE to 8 requires the specification of static weigths. &
 &Flat (uniformative) weights will be used instead to circumvent the lack of specification. They might become informative in & 
 &analysis depending on keyword CMODWEIGHTS and related ones.'
    if (allocated(cl_imvec).EQV..true.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Got cl_imvec allocated with arethere_stwghts .false. in ncdm_convert_ascii. This is a bug.'
      call fexit()
    end if
    allocate(cl_imvec(ncdm_nfeats_ini))
    cl_imvec(:) = 1.
    ncdm_arethere_stwghts = .true.  !well flat, but are there now
    write_weights(1) = .false.      !do not write them in converted file though
  end if
  if (((cdis_crit.eq.2).OR.(cdis_crit.eq.9)).AND.(ncdm_arethere_dywghts.EQV..false.).AND.(ncdm_isdmonas.EQV..true.)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. With this setting of CDISTANCE and analysis we should not be here anymore. This is a bug.'
    call fexit()
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(cdis_crit.ne.8).AND.(ncdm_isdmonas.EQV..true.)) then
    write(ilog,*) 'Warning. The detected static weights in the input net-cdf file will be discarded in analysis &
 &due to the value of CDISTANCE.'
  end if
  if ((ncdm_arethere_dywghts.EQV..true.).AND.((cdis_crit.ne.2).AND.(cdis_crit.ne.9)).AND.(ncdm_isdmonas.EQV..true.)) then
    write(ilog,*) 'Warning. The detected dynamic weights in the input ascii file will be discarded in analysis &
 &due to the value of CDISTANCE.'
  end if
!
! here we have finished reading and we know what is in the ascii file (all flags are set). Set the masks to the lines of cludata
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    allocate(ncdm_featvals_includata(ncdm_nfeats_ini))
    allocate(ncdm_dynwghts_includata(ncdm_nfeats_ini))
    call ncdm_set_cludatamasks(ncdm_featvals_includata,ncdm_nfeats_ini,ncdm_dynwghts_includata)  !populatin globals for write_ncfl
  else
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    allocate(ncdm_featvals_includata(ncdm_nfeats_ini))
    call ncdm_set_cludatamasks(ncdm_featvals_includata,ncdm_nfeats_ini)
  end if
!
! write binary
  call ncdm_write_ncfl(write_weights)
  if (changed.EQV..true.) then
    ncdm_docheck_r = .true.
    changed = .false.
  end if
!
! possible fast return if no analysis is to be performed
  if (ncdm_isdmonas.EQV..false.) then
    if (allocated(cludata).EQV..true.) deallocate(cludata)
    if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec) 
    if (allocated(data_on_row).EQV..true.) deallocate(data_on_row)
    if (allocated(tmpvec).EQV..true.) deallocate(tmpvec)
    if (allocated(tmparr).EQV..true.) deallocate(tmparr)
    return
  end if
!
! adjust cludata if necessary
 if ((cstorecalc.ne.1).OR.(ncdm_isthere_framesfl.EQV..true.).OR.(ncdm_isthere_cfile.EQV..true.)) then
    if (allocated(tmparr).EQV..true.) deallocate(tmparr)
    if ((ncdm_arethere_dywghts.EQV..true.).AND.((cdis_crit.eq.2).OR.(cdis_crit.eq.9))) then 
      allocate(tmparr(2*ncdm_nfeats,ncdm_nframes)) !pruned data set 
        if (allocated(featvals_inarr).EQV..true.) deallocate(featvals_inarr)
        if (allocated(dynwghts_inarr).EQV..true.) deallocate(dynwghts_inarr)
        allocate(featvals_inarr(ncdm_nfeats))
        allocate(dynwghts_inarr(ncdm_nfeats))
        call ncdm_set_cludatamasks(featvals_inarr,ncdm_nfeats,dynwghts_inarr)  !indexes to tmparr  
        tmparr(featvals_inarr,:) = cludata(ncdm_featvals_includata(ncdm_whichfeats),ncdm_whichframes)
        tmparr(dynwghts_inarr,:) = cludata(ncdm_dynwghts_includata(ncdm_whichfeats),ncdm_whichframes)
        if (allocated(featvals_inarr).EQV..true.) deallocate(featvals_inarr)
        if (allocated(dynwghts_inarr).EQV..true.) deallocate(dynwghts_inarr)
    else 
      if (ncdm_arethere_dywghts.EQV..true.) then
        write(ilog,*) 'Warning. Discarding dynamic weights due to value of CDISTANCE.'
        ncdm_arethere_dywghts = .false.
      end if
      allocate(tmparr(ncdm_nfeats,ncdm_nframes))
      if (allocated(featvals_inarr).EQV..true.) deallocate(featvals_inarr)
      allocate(featvals_inarr(ncdm_nfeats))
      call ncdm_set_cludatamasks(featvals_inarr,ncdm_nfeats)  !indexes to tmparr  
      tmparr(featvals_inarr,:) = cludata(ncdm_featvals_includata(ncdm_whichfeats),ncdm_whichframes)
      if (allocated(featvals_inarr).EQV..true.) deallocate(featvals_inarr)
      if (allocated(dynwghts_inarr).EQV..true.) deallocate(dynwghts_inarr)
    end if
!
    if (allocated(tmpvec).EQV..true.) deallocate(tmpvec)
    deallocate(cludata)
    if ((ncdm_arethere_dywghts.EQV..true.)) then
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))
      if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
      if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
      allocate(ncdm_featvals_includata(ncdm_nfeats))
      allocate(ncdm_dynwghts_includata(ncdm_nfeats))
      call ncdm_set_cludatamasks(ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata) !re-populating globals on final dim 
    else 
      allocate(cludata(ncdm_nfeats,ncdm_nframes))
      if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
      if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
      allocate(ncdm_featvals_includata(ncdm_nfeats))
      call ncdm_set_cludatamasks(ncdm_featvals_includata,ncdm_nfeats) !re-populating globals on final dim
    end if
    do iframe=1,ncdm_nframes
      cludata(:,iframe) = tmparr(:,iframe)
    end do
!
    if ((ncdm_arethere_stwghts.EQV..true.).AND.(cdis_crit.eq.8)) then
      tmparr(1:ncdm_nfeats,1) = cl_imvec(ncdm_whichfeats)
      deallocate(cl_imvec)
      allocate(cl_imvec(ncdm_nfeats))
      cl_imvec(:) = tmparr(1:ncdm_nfeats,1)
    else if (ncdm_arethere_stwghts.EQV..true.) then
      write(ilog,*) 'Warning. Discarding static weights due to value of CDISTANCE.'
      ncdm_arethere_stwghts = .false.
      deallocate(cl_imvec)
    end if
!
    if (allocated(tmparr).EQV..true.) deallocate(tmparr)
!
  else  !if ((cstorecalc.ne.1).OR.(ncdm_isthere_framesfl.EQV..true.).OR.(ncdm_isthere_cfile.EQV..true.))
!
! just reduce allocation if possible due to cdistance (we have done this already in case framesfile, ccollect or cfile are present)
    if (ncdm_nframes.ne.ncdm_nframes_ini) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Gotten ncdm_nframes.ne.ncdm_nframes_ini with no frames file, ccollect or cfile. This is a bug.'
      call fexit()
    end if
    if (ncdm_nfeats.ne.ncdm_nfeats_ini) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Gotten ncdm_nfeats.ne.ncdm_nfeats_ini with no frames file, ccollect or cfile. This is a bug.'
      call fexit()
    end if
    if ( (ncdm_arethere_dywghts.EQV..true.) ) then
      if ((cdis_crit.ne.2).AND.(cdis_crit.ne.9)) then  !it means read in 9 style
        write(ilog,*) 'Warning. Discarding dynamic weights due to value of CDISTANCE.'
        ncdm_arethere_dywghts = .false.
        if (allocated(tmparr).EQV..true.) deallocate(tmparr)
        allocate(tmparr(ncdm_nfeats,ncdm_nframes)) 
        tmparr(:,:) = cludata(ncdm_whichfeats,ncdm_whichframes)  !read 9 style
        deallocate(cludata)
        allocate(cludata(ncdm_nfeats,ncdm_nframes))
        cludata(:,:) = tmparr(:,:)
      end if
    end if
    if ((ncdm_arethere_stwghts.EQV..true.).AND.(cdis_crit.eq.8)) then
      if (allocated(cl_imvec).EQV..false.) then
        write(ilog,*) 
        write(ilog,*) 'Fatal. cl_imvec not allocated with static weights and CDISTANCE 8. This is a bug.' 
        call fexit()
      end if
      if (allocated(tmpvec).EQV..true.) deallocate(tmpvec)
      allocate(tmpvec(ncdm_nfeats))
      tmpvec(:) = cl_imvec(ncdm_whichfeats)
      deallocate(cl_imvec)
      allocate(cl_imvec(ncdm_nfeats))
      cl_imvec(:) = tmpvec(:)
    else if (ncdm_arethere_stwghts.EQV..true.) then
      write(ilog,*) 'Warning. Discarding static weights due to value of CDISTANCE.'
      ncdm_arethere_stwghts = .false.
      deallocate(cl_imvec)  
    end if
  end if  !if ((cstorecalc.ne.1).OR.(ncdm_isthere_framesfl.EQV..true.).OR.(ncdm_isthere_cfile.EQV..true.))
!
  if (allocated(data_on_row).EQV..true.) deallocate(data_on_row)
  if (allocated(tmpvec).EQV..true.) deallocate(tmpvec)
  if (allocated(tmparr).EQV..true.) deallocate(tmparr)
  if (allocated(featvals_inarr).EQV..true.) deallocate(featvals_inarr)
  if (allocated(dynwghts_inarr).EQV..true.) deallocate(dynwghts_inarr)
!
end subroutine ncdm_convert_ascii
!
!---------------------------------------------------------------------------
!
subroutine ncdm_set_cludatamasks(featsvec,nfeats,dywghtsvec) !return index of feats and weights when dim is nfeats according to cdis
!
  use iounit
  use clusters, ONLY: cdis_crit
  use ncdm, ONLY: ncdm_arethere_dywghts,ncdm_donc,ncdm_isdmonas
!
  implicit none
!
  integer, INTENT(IN)  :: nfeats
!
  integer, INTENT(OUT) :: featsvec(nfeats)
  integer, INTENT(OUT), OPTIONAL:: dywghtsvec(nfeats)
!
  integer featcounter,ifeat
!
! initialize dimensions
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (present(dywghtsvec).EQV..false.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Calling ncdm_set_cludatamasks with wrong input combinations (1). This is a bug.'
      call fexit()
    end if
  else
    if (present(dywghtsvec).EQV..true.) then
      write(ilog,*) 
      write(ilog,*) 'Fatal. Calling ncdm_set_cludatamasks with wrong input combinations (2). This is a bug.'
      call fexit()
    end if
  end if
!
! masks depending on whether we do the analysis or not and possibly on cdistance
  if ((ncdm_donc.EQV..true.).OR.(ncdm_isdmonas.EQV..true.)) then
    if (ncdm_arethere_dywghts.EQV..true.) then
      if (cdis_crit.eq.2) then
        featcounter = 0
        do ifeat=1,nfeats
          featcounter = featcounter + 1
          featsvec(ifeat) = 2*featcounter - 1
          dywghtsvec(ifeat) = 2*featcounter
        end do
      else  !9 style
        do ifeat=1,nfeats
          featsvec(ifeat) = ifeat
          dywghtsvec(ifeat) = ifeat + nfeats
        end do
      end if
    else
      do ifeat=1,nfeats
        featsvec(ifeat) = ifeat
      end do
    end if
  else
    if (ncdm_arethere_dywghts.EQV..true.) then
      do ifeat=1,nfeats
        featsvec(ifeat) = ifeat
        dywghtsvec(ifeat) = ifeat + nfeats
      end do
    else
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
subroutine ncdm_write_ncfl(write_weights)
!
  use iounit
  use system, ONLY: basename,bleng
  use clusters, ONLY: cdis_crit,cludata,cl_imvec
  use netcdf
  use ncdm
!
  implicit none
!
  integer, parameter :: nvar = 3      !how many variable do we need, 1 is the featureval array, 2 and 3 are weights
  integer, parameter :: ndim = 2      !number of dimensions that we need to index our variables
!
  integer t1,t2                       !to determine string limits
  integer i,freeunit,counter,ifeat
  integer nframes,nfeats              !local values from size(cludata)
  integer ncid_db                     !net-cdf databese id
  integer ncid_dim(ndim)              !net-cdf dimensions ids
  integer ncid_var(nvar)              !net-cdf variables ids
  integer featcounter
  character(MAXSTRLEN) attstring
  character(len=(bleng+11)) ncfn_w    !net-cdf database (file) name to be written
  logical exists
  logical, INTENT(IN):: write_weights(2)
!
  if (ncdm_docheck_r.EQV..false.) then                 !not here cause checking what has been read as input net-cdf
    t1 = 1
    t2 = bleng + 11
    ncfn_w(t1:t2) = basename(1:bleng) // "_convert.nc"
  else
    t1 = 1
    t2 = bleng + 11
    ncfn_w(t1:t2) = basename(1:bleng) // "_checked.nc"
  end if
!
  ncid_db = 0
  ncid_dim = 0
!
! initialize dimensions
  if (ncdm_arethere_dywghts.EQV..false.) then
    nfeats = size(cludata,dim=1)
  else
    nfeats = size(cludata,dim=1)/2
  end if
  nframes = size(cludata,dim=2)
!
  inquire(file=ncfn_w,exist=exists)
  if (exists.EQV..true.) then
    ncid_db = freeunit()
    open(unit=ncid_db,file=ncfn_w,status='old')
    close(unit=ncid_db,status='delete')
  end if
!
  call check_fileionetcdf( nf90_create(path=ncfn_w, cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid_db) ) !large data file
!
! global attributes
  attstring(1:19) = "CAMPARI DATA MINING"
  call check_fileionetcdf( nf90_put_att(ncid_db, NF90_GLOBAL, "program", attstring(1:19)) )
!
! dimensions
  attstring(1:6) = "nfeats" 
  call check_fileionetcdf( nf90_def_dim(ncid_db, attstring(1:6), nfeats, ncid_dim(1)) )  
  attstring(1:7) = "nframes"
  call check_fileionetcdf( nf90_def_dim(ncid_db, attstring(1:7), NF90_UNLIMITED, ncid_dim(2)) )
!
! define variables and put relevat attributes
  !featurevals variable
  call strlims(ncdm_nm_feats,t1,t2)
  attstring(t1:t2) = ncdm_nm_feats
  call check_fileionetcdf( nf90_def_var(ncid_db, attstring(t1:t2), NF90_FLOAT, (/ ncid_dim(1), ncid_dim(2) /), ncid_var(1)) )
  if (ncdm_is_periodic.EQV..true.) then 
    attstring(1:14) = "periodic_range"
    call check_fileionetcdf( nf90_put_att(ncid_db, ncid_var(1), attstring(1:14), ncdm_periodic_rng) ) 
    if (((ncdm_donc.EQV..true.).OR.(ncdm_isdmonas.EQV..true.)).AND.((cdis_crit.eq.1).OR.(cdis_crit.eq.2))) then
      if ( ( minval(cludata(ncdm_featvals_includata,:)).lt.ncdm_periodic_rng(1) ).OR.( &
 &           maxval(cludata(ncdm_featvals_includata,:)).gt.ncdm_periodic_rng(2)      ) ) then
        write(ilog,*) 'Warning. The provided input data exceed the reference periodic interval. This is not corrected &
 &in output files, but it is taken care of during analysis with periodic CDISTANCE, as specified.'
      end if
    end if
  end if
  !possible weights
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(write_weights(1).EQV..true.)) then  !then we want to output them
    attstring(1:14) = "static_weights"
    call check_fileionetcdf( nf90_def_var(ncid_db, attstring(1:14), NF90_FLOAT, ncid_dim(1), ncid_var(2)) )
  end if
  if ((ncdm_arethere_dywghts.EQV..true.).AND.(write_weights(2).EQV..true.)) then  !then we want to output them
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
  call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(1), values = cludata(ncdm_featvals_includata,1:nframes), &
 &                                       start = (/ 1, 1 /), count = (/ nfeats, nframes /) ) )
!  static weights
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(write_weights(1).EQV..true.)) then
    call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(2), values = cl_imvec(1:nfeats), &
 &                                       start = (/ 1 /) , count =  (/ nfeats /)  ) ) 
  end if 
!  dynamic weights
  if ((ncdm_arethere_dywghts.EQV..true.).AND.(write_weights(2).EQV..true.)) then
  call check_fileionetcdf( nf90_put_var( ncid_db, ncid_var(3), values = cludata(ncdm_dynwghts_includata,1:nframes), &
 &                                           start = (/ 1, 1 /), count = (/ nfeats, nframes /) ) )
  end if 
!
  call check_fileionetcdf( nf90_close(ncid_db) )
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
  integer, parameter :: nvar = 3               !variable ids, 1 is featurevals, 2 and 3 are possible weights
  integer, parameter :: ndim = 2               !number of dimensions that we need to index our variables
!
  integer t1,t2,counter,ifeat,iframe,iatt,ivar
  integer ncid_db                       !net-cdf databese id
  integer ncid_dim(ndim)                !net-cdf dimensions ids
  integer ncid_var(nvar)                !net-cdf variables ids
  integer nvar_db                       !how many variables are in the databased (used to look for weights)
  integer natts_vardb                   !number of attributes of db variable
!
  logical atrue
  logical write_weights(2)              !whether or not should be written to file. (1) is static, (2) dynamic
!
  character(MAXSTRLEN) attstring
!
  RTYPE, allocatable :: tmparr(:,:)  !array that helps dealing with cfile and ccollect while keeping memory requirement low
!
  ncid_db = 0
  ncid_dim = 0
  ncid_var = 0
  attstring(:) = " "
  atrue = .true.
  write_weights(:) = .true.  !if there are weights, we assume they are good
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
  call strlims(ncdm_nm_feats,t1,t2)
  call check_fileionetcdf( nf90_inq_varid(ncid_db, ncdm_nm_feats(t1:t2), ncid_var(1)) )
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
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(ncdm_arethere_dywghts.EQV..true.)) then 
    write(ilog,*)
    write(ilog,*) 'Warning. Gotten both weights present in input ascii file. Only one type will be used (if any) in analysis, &
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
    ncdm_arethere_stwghts = .true.  !pretending there are
    write_weights(1) = .false.      !very important also here not only in writing check output to discern whether pretending of not
  end if
  if (((cdis_crit.eq.2).OR.(cdis_crit.eq.9)).AND.(ncdm_arethere_dywghts.EQV..false.)) then  !locally adaptive weights
    ncdm_arethere_dywghts = .true.  !pretending there are
    write_weights(2) = .false.      !very important also here not only in writing check output to discern whether pretending of not
    write(ilog,*) 'Warning. The current setting of CDISTANCE (to 2 or 9) requires the specification of dynamic weigths. &
 &Flat (uniformative) weights will be used instead to circumvent the lack of specification or identification. &
 &Their might become informative in analysis depending on keyword CMODWEIGHTS and related ones.'
  end if
  if ((ncdm_arethere_stwghts.EQV..true.).AND.(cdis_crit.ne.8)) then  !very important
    write(ilog,*) 'Warning. The detected static weights in the input net-cdf file will be discarded. This entails &
  &both the possible check-file and the analysis and is due to the value of CDISTANCE.'
    ncdm_arethere_stwghts = .false.
  end if
  if ((ncdm_arethere_dywghts.EQV..true.).AND.((cdis_crit.ne.2).AND.(cdis_crit.ne.9))) then !very important
    write(ilog,*) 'Warning. The detected dyatic weights in the input net-cdf file will be discarded. This entails &
  &both the possible check-file and the analysis and is due to the value of CDISTANCE.'
    ncdm_arethere_dywghts = .false.
  end if
!
! allocate cludata and possibly cl_imvec
  if (allocated(cludata).EQV..true.) deallocate(cludata)
  if (ncdm_arethere_stwghts.EQV..true.) then  
    if (allocated(cl_imvec).EQV..true.) deallocate(cl_imvec)
    allocate(cl_imvec(ncdm_nfeats))
    if (ncdm_arethere_dywghts.EQV..false.) then
      allocate(cludata(ncdm_nfeats,ncdm_nframes))    !already pruned
    else
      allocate(cludata(2*ncdm_nfeats,ncdm_nframes))    !already pruned
    end if
  else if (ncdm_arethere_dywghts.EQV..true.) then  
    allocate(cludata(2*ncdm_nfeats,ncdm_nframes))  !space also for dynamic weights
  else 
    allocate(cludata(ncdm_nfeats,ncdm_nframes))
  end if
!
! establishes the masks
  if (ncdm_arethere_dywghts.EQV..true.) then
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    allocate(ncdm_featvals_includata(ncdm_nfeats))
    allocate(ncdm_dynwghts_includata(ncdm_nfeats))
    call ncdm_set_cludatamasks(ncdm_featvals_includata,ncdm_nfeats,ncdm_dynwghts_includata)
  else 
    if (allocated(ncdm_featvals_includata).EQV..true.) deallocate(ncdm_featvals_includata)
    if (allocated(ncdm_dynwghts_includata).EQV..true.) deallocate(ncdm_dynwghts_includata)
    allocate(ncdm_featvals_includata(ncdm_nfeats))
    call ncdm_set_cludatamasks(ncdm_featvals_includata,ncdm_nfeats)
  end if
!
! populate cludata
!  if (cdis_crit.ne.2) then
    do iframe=1,ncdm_nframes
      do ifeat=1,ncdm_nfeats
        call check_fileionetcdf( nf90_get_var(ncid_db, ncid_var(1), cludata(ncdm_featvals_includata(ifeat),iframe:iframe), &
 &           start = (/ ncdm_whichfeats(ifeat), ncdm_whichframes(iframe) /), count = (/ 1, 1 /) ))
      end do
    end do 
!
! in case, account for weights, ncdm_docheck_r is to write proper check file, but it slows down if cdis_crit neq 8 (see fill_)
  if (ncdm_arethere_stwghts.EQV..true.) then 
    if (write_weights(1).EQV..true.) then !no fake
      do ifeat=1,ncdm_nfeats
        call check_fileionetcdf( nf90_get_var(ncid_db, ncid_var(2), cl_imvec(ifeat:ifeat), &
 &           start = (/ ncdm_whichfeats(ifeat) /), count = (/ 1 /) ))
      end do
    else
      cl_imvec(:) = 1.
    end if
  end if
  if ((ncdm_arethere_dywghts.EQV..true.).AND.((cdis_crit.eq.9).OR.(cdis_crit.eq.2))) then
    if (write_weights(2).EQV..true.) then !no fake
      do iframe=1,ncdm_nframes
        do ifeat=1,ncdm_nfeats
          call check_fileionetcdf( nf90_get_var(ncid_db, ncid_var(3), cludata(ncdm_dynwghts_includata(ifeat),iframe:iframe), &
 &             start = (/ ncdm_whichfeats(ifeat), ncdm_whichframes(iframe) /), count = (/ 1, 1 /) ))
        end do
      end do 
    else
      do iframe=1,ncdm_nframes
        do ifeat=1,ncdm_nfeats
          cludata(ncdm_dynwghts_includata(ifeat),iframe) = 1.
        end do
      end do
    end if
  end if
!
  if (ncdm_docheck_r.EQV..true.) then 
    call ncdm_write_ncfl(write_weights)
  end if 
!
end subroutine ncdm_read_ncfl
!
!---------------------------------------------------------------------------
!
subroutine ncdm_fill_cludata 
!
  use iounit
  use system, ONLY: nequil,nsim
  use clusters, ONLY: cdis_crit,calcsz,clstsz,cstored,cdofset,cl_imvec,sumssz,cludata,clagt_msm
  use pdb, ONLY: framecnt,select_frames,framelst
  use ncdm
 
!
  implicit none
!
  integer iu,freeunit,aone
  integer iframe,ifeat,iweight,counter
!
  RTYPE, allocatable:: tmparr(:,:)
!
  aone = 1
!
! from distance criterion to allocation of features and weights
  if (cdis_crit.eq.1) calcsz = ncdm_nfeats   ! dihedrals
  if (cdis_crit.eq.2) calcsz = 2*ncdm_nfeats ! dihedrals with time-dep. weights
  if (cdis_crit.eq.7) calcsz = ncdm_nfeats 
  if (cdis_crit.eq.8) calcsz = ncdm_nfeats    ! DRMS with static weights -----> in cl_imvec
  if (cdis_crit.eq.9) calcsz = 2*ncdm_nfeats  ! DRMS with time-dep. weights
!
! from distance criterion to allocation of features only
  clstsz = ncdm_nfeats    ! always this way in here
!
! degrees of freedom for DRMS
  if ((cdis_crit.ge.7).AND.(cdis_crit.le.9)) then
    if (allocated(cdofset).EQV..true.) deallocate(cdofset)
  end if
!
! checks
  if ( (cdis_crit.eq.8).AND.(allocated(cl_imvec).EQV..false.) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten cl_imvec not allocated with cdistance 8 in ncdm_fill_cludata. This is a bug.'
    call fexit()
  end if
  if ( (cdis_crit.ne.8).AND.(allocated(cl_imvec).EQV..true.) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten cl_imvec allocated with cdistance not 8 in ncdm_fill_cludata. This is a bug.'
    call fexit()
  end if
  if ( ((cdis_crit.eq.2).OR.(cdis_crit.eq.9)).AND.(size(cludata,dim=1).ne.(2*ncdm_nfeats)) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten mismatch on size of dim 1 of cludata with cdistance 2 or 9 in ncdm_fill_cludata. &
 &This is a bug.'
    call fexit()
  else if ( ((cdis_crit.eq.1).OR.(cdis_crit.eq.7).OR.(cdis_crit.eq.8)).AND.(size(cludata,dim=1).ne.(ncdm_nfeats)) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten mismatch on size of dim 1 of cludata with cdistance 1, 7 or 8 in ncdm_fill_cludata. &
 &This is a bug.',ncdm_nfeats,size(cludata,dim=1)
    call fexit()
  end if
  if ( size(cludata,dim=2).ne.(ncdm_nframes) ) then
    write(ilog,*) 
    write(ilog,*) 'Fatal. Gotten mismatch on size of dim 2 of cludata in ncdm_fill_cludata. This is a bug.'
    call fexit()
  end if
!
! other allocations and initializations
  cstored = ncdm_nframes
  nequil = 0
  nsim = ncdm_nframes_ini !cstored
  framecnt = ncdm_nframes
  framelst = ncdm_whichframes
  select_frames = .false.
  !write(ilog,*) 'framelst',framelst
  !if (ncdm_isthere_framesfl.EQV..false.) then
  !  nsim = ncdm_nframes_ini  !correct for both input ascii file or input net-cdf file
  !  select_frames = .false.
  !else 
  !  nsim = ncdm_nframes 
  !  framecnt = ncdm_nframes
  !  framelst = ncdm_whichframes
  !  select_frames = .true.
  !end if
  if ((cdis_crit.eq.2).OR.(cdis_crit.eq.9)) then
    sumssz = 5
  else
    sumssz = 2
  end if
!
  call build_sconnect(clagt_msm,aone)
!
end subroutine ncdm_fill_cludata
!
!---------------------------------------------------------------------------
!
!
subroutine ncdm_scale_periodicity
!
  use iounit
  use clusters, ONLY: cludata,cdis_crit
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
  if (ncdm_periodic_rng(1).ge.ncdm_periodic_rng(2)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. The values detected periodicity for the periodicity range are ineligible, &
 &as they are the same value. Please, either specify eligible values in the net-cdf file or remove the attribute.'
    call fexit()
  end if
  if (ncdm_arethere_dywghts.EQV..false.) then
    if (cdis_crit.ne.1) then
      write(ilog,*)
      write(ilog,*) 'Fatal. Wrong combination of parameters in ncdm_scale_periodicity (1). This is a bug.'
      call fexit()
    end if
  else
    if (cdis_crit.ne.2) then
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
      do ifeat=1,ncdm_nfeats
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
  offset = (ncdm_periodic_rng(2) + ncdm_periodic_rng(1))/2  !center of periodicity interval
  ncdm_periodic_rng(2) = ncdm_periodic_rng(2) - offset      !re-centered margins of periodic interval   
  ncdm_periodic_rng(1) = ncdm_periodic_rng(1) - offset
  do iframe=1,ncdm_nframes  !mapping to [-180:180]: all features of all snapshosts
    do ifeat=1,ncdm_nfeats  
      tmpr = cludata(ncdm_featvals_includata(ifeat),iframe)
      cludata(ncdm_featvals_includata(ifeat),iframe) = ( (tmpr - offset) * rightbound )/ncdm_periodic_rng(2)
    end do
  end do
!
! final check
  if ((minval(cludata).lt.leftbound).OR.(maxval(cludata).gt.rightbound)) then
    write(ilog,*)
    write(ilog,*) 'Fatal. Periodicity bound exceeded. This is a bug.'
    call fexit()
  end if
!
end subroutine ncdm_scale_periodicity
!---------------------------------------------------------------------------
!
#endif
!
!---------------------------------------------------------------------------
