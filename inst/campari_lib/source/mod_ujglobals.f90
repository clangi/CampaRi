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
! MAIN AUTHOR:   Adam Steffen                                              !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module ujglobals
!
   integer MAXUJDOF,MAXSOLU
   parameter (MAXUJDOF=90)
   parameter (MAXSOLU=40)
!
   integer cur_ujsz,uj_deadcnt,dj_deadcnt,do_deadcnt,nc_deadcnt,ujmaxsz,ujminsz,npucmoves,torcrmode
   integer dj_totcnt,do_totcnt,nc_totcnt,uj_totcnt,dj_wrncnt(5),do_wrncnt(5),nc_wrncnt(5)
   integer dj_wrnlmt(5),do_wrnlmt(5),nc_wrnlmt(5),maxcrtries,pc_wrncnt(6),pc_wrnlmt(6)
   integer cur_torcrsz,torcrmaxsz,torcrminsz,nuccrmaxsz,nuccrminsz,torcrmaxsz2,torcrminsz2
   RTYPE UJ_params(6)
   RTYPE, ALLOCATABLE:: obang(:)
!
end module ujglobals
!
