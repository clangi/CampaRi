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
!--------------------------------------------------------------------------!
!
!
! ncdm_nframes_ini        : preliminary number of frames, i.e. before ccollect or framesfile
! ncdm_nfeats_ini         : preliminary number of features, i.e. before possible cfile
! ncdm_nframes            : eff. number of frames after all changes, the only ones that is seen by cludata and mining routines
! ncdm_nfeats             : eff. number of features after all changes, the only ones that is seen by cludata and mining routines
!
! ncdm_fn_r               : input net-cdf file for reading
! ncdm_fn_as              : file-name for the possible feature ascii file
! ncdm_nm_feats           : name of the variable that contains featurevals in the input net-cdf file
! ncdm_docheck_r          : whether or not to print out what campari has read from input net-cdf (mainly for debug)
! ncdm_donc               : whether or not do the reading of the net-df file and all the associated analysis
! ncdm_doasconv           : whether or not to convert the ascii input file  
! ncdm_isdmonas           : in case we specify both an ascii (to be converted) and a net-cdf, is the ascii to be data-mined?
!
! ncdm_periodic_rng       : in case we have, during conversion, to write down the periodic attribute
! ncdm_is_periodic        : just a flag that tells whether the data have to be assumed periodic
!
! ncdm_nm_cfile           : name of possible cfile
! ncdm_isthere_cfile      : logical flag that answer: was a cfile specified in the key file? 
! ncdm_whichfeats         : index array that contains the original numbering of the features that will be used 
! ncdm_featsmask          : logical array the entry of which is true if the corresponding feature is eligible for analysis
!
! ncdm_nm_framesfl        : name of possible frames file 
! ncdm_isthere_framesfl   : logical flag that answer: was a cfile specified in the key file? 
! ncdm_whichframes        : index array that contains the original numbering of the frames that will be used
! ncdm_framesmask         : logical array the entry of which is true if the corresponding frame is eligible for analysis
!
! ncdm_arethere_stwgths   : are there static weights or not?
! ncdm_arethere_dywgths   : are there dynamics weights or not?
!
! ncdm_featvals_includata : what lines of cludata store the values of the features?
! ncdm_dynwghts_includata : what lines of cludata store the values of the dynamic weights?
! 
!--------------------------------------------------------------------------!
module ncdm          ! net-cdf data mining
!
  integer ncdm_nframes_ini,ncdm_nfeats_ini
  integer ncdm_nframes,ncdm_nfeats
!
  character(MAXSTRLEN) ncdm_fn_r,ncdm_fn_as
  character(MAXSTRLEN) ncdm_nm_feats
  character(MAXSTRLEN) ncdm_nm_cfile,ncdm_nm_framesfl  
!
  RTYPE ncdm_periodic_rng(2)
!
  logical ncdm_docheck_r,ncdm_donc,ncdm_doasconv,ncdm_isdmonas, &
 &ncdm_is_periodic,ncdm_isthere_cfile,ncdm_isthere_framesfl, &
 &ncdm_arethere_stwghts,ncdm_arethere_dywghts
!
  integer, allocatable:: ncdm_whichfeats(:),ncdm_whichframes(:)
  integer, allocatable:: ncdm_featvals_includata(:),ncdm_dynwghts_includata(:)
!
  logical, allocatable:: ncdm_framesmask(:),ncdm_featsmask(:)
!
! interfaces
  interface ncdm_set_cludatamasks
    subroutine ncdm_set_cludatamasks(featsvec,nfeats,dywghtsvec)
      integer, INTENT(IN)  :: nfeats
      integer, INTENT(OUT) :: featsvec(nfeats)
      integer, INTENT(OUT), OPTIONAL:: dywghtsvec(nfeats)
      integer featcounter,ifeat
    end subroutine ncdm_set_cludatamasks
  end interface ncdm_set_cludatamasks
!
end module ncdm
!
