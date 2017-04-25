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
! CONTRIBUTIONS: Nicholas Lyle                                             !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
module interfaces
!
  interface
    function handwrapped_getpid() bind(C, name='getpid')
      use, INTRINSIC:: ISO_C_BINDING
      implicit none
      integer(c_int) :: handwrapped_getpid
    end function handwrapped_getpid
  end interface
!
  interface
    function handwrapped_gethostname(hnam,hnaml) bind(C, name='gethostname')
      use, INTRINSIC:: ISO_C_BINDING
      implicit none
      character(len=1,kind=c_char), dimension(*), intent(out) :: hnam
      integer(c_size_t):: hnaml
!      type(C_PTR):: hnam
      integer(c_int) :: handwrapped_gethostname
    end function handwrapped_gethostname
  end interface
!
  interface custom_verfcoverr
    function custom_verfcoverr(d1v,pfac)
      implicit none
      RTYPE, INTENT(IN):: pfac
      RTYPE, INTENT(IN):: d1v(:)
      RTYPE  custom_verfcoverr(size(d1v),2)
    end function custom_verfcoverr
  end interface custom_verfcoverr
!
  interface fsort
    subroutine fsort(ldim,up,flist,ilist,ilist2,ilist3,flist2,flist3)
      integer ldim
      logical up
      RTYPE flist(*)
      integer, OPTIONAL:: ilist(*),ilist2(*),ilist3(*)
      RTYPE, OPTIONAL:: flist2(*),flist3(*)
    end subroutine fsort
  end interface fsort
!
  interface isort
    subroutine isort(ldim,up,ilist,flist,ilist2,ilist3,flist2,flist3)
      integer ldim
      logical up
      integer ilist(*)
      integer, OPTIONAL:: ilist2(*),ilist3(*)
      RTYPE, OPTIONAL:: flist(*),flist2(*),flist3(*)
    end subroutine isort
  end interface isort
!
  interface merge_sort
    recursive subroutine merge_isort(ldim,up,list,olist,ilo,ihi,idxmap,olist2)
      integer ldim,ilo,ihi
      logical up
      integer list(*),olist(*)
      integer, OPTIONAL:: idxmap(*),olist2(*)
    end subroutine merge_isort
    recursive subroutine merge_fsort(ldim,up,list,olist,ilo,ihi,idxmap,olist2)
      integer ldim,ilo,ihi
      logical up
      RTYPE list(*),olist(*)
      integer, OPTIONAL:: idxmap(*),olist2(*)
    end subroutine merge_fsort
    recursive subroutine merge_ffsort(ldim,up,list,olist,ilo,ihi,idxmap,olist2)
      integer ldim,ilo,ihi
      logical up
      real(KIND=4) list(*),olist(*)
      integer, OPTIONAL:: idxmap(*),olist2(*)
    end subroutine merge_ffsort
  end interface merge_sort
!
  interface resize_vec
    subroutine resize_ivec(oldsz,iv,newsz,preserve)
      integer, INTENT(IN):: oldsz,newsz
      logical, INTENT(IN):: preserve
      integer, ALLOCATABLE, INTENT(INOUT):: iv(:)
    end subroutine resize_ivec
    subroutine resize_fvec(oldsz,iv,newsz,preserve)
      integer, INTENT(IN):: oldsz,newsz
      logical, INTENT(IN):: preserve
      RTYPE, ALLOCATABLE, INTENT(INOUT):: iv(:)
    end subroutine resize_fvec
    subroutine resize_ffvec(oldsz,iv,newsz,preserve)
      integer, INTENT(IN):: oldsz,newsz
      logical, INTENT(IN):: preserve
      real(KIND=4), ALLOCATABLE, INTENT(INOUT):: iv(:)
    end subroutine resize_ffvec
  end interface resize_vec
!
  interface scluster_resizelst
    subroutine scluster_resizelst(currentalsz,it)
      use clusters
      integer currentalsz
      type(t_scluster), ALLOCATABLE, INTENT(IN OUT):: it(:)
    end subroutine scluster_resizelst
  end interface scluster_resizelst
!
  interface binary_search
    subroutine fbinary_search(arsz,arvec,keyv,ifound)
      integer arsz,ifound
      RTYPE keyv
      RTYPE arvec(*)
    end subroutine fbinary_search
    subroutine ffbinary_search(arsz,arvec,keyv,ifound)
      integer arsz,ifound
      real(KIND=4) keyv
      real(KIND=4) arvec(*)
    end subroutine ffbinary_search
    subroutine ibinary_search(arsz,arvec,keyv,ifound)
      integer arsz,keyv,ifound
      integer arvec(*)
    end subroutine ibinary_search
  end interface binary_search
!
#ifdef ENABLE_MPI
  interface wl_parallel_combine
    subroutine wl_parallel_combine(gtmp,gtmp2d)
      integer, OPTIONAL, ALLOCATABLE, INTENT(IN OUT):: gtmp2d(:,:,:),gtmp(:,:)
    end subroutine wl_parallel_combine
  end interface wl_parallel_combine
#endif
!
end module interfaces 
