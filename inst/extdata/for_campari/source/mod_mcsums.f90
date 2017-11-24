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
! MAIN AUTHOR:   Rohit Pappu                                               !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! file handles for running output:          idihed,ipolmr,....
! write-out frequencies for running output: xyzout,enout,accout,phout,ensout,rstout
! timing variables:                         time_energy,...
! flush interval for running output:        time_flush
! firsthole:                                true if no holes calculations have been done
!
module mcsums
!
  logical firsthole
!
! time variables
  RTYPE time_energy,time_analysis,time_struc,time_comm
  RTYPE time_ph,time_holo,time_nbl,time_pme,time_flush
!
! file handles
  integer idihed,ipolmr,iholes,ipht,iremc,idssp,idihed_ring
  integer iacc,ipers,isasa,iene,ipdbtraj,iens,iretr,igpc,irmsd,itraj,itrajhlp,ibprtnr
!
! control variables for running output
  integer xyzout,xyzmode,enout,accout,phout,nstep,rstout,ensout
!
end module mcsums
!
