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
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
! WRAPPER: Davide Garolini                                                 !
!                                                                          !
!--------------------------------------------------------------------------!

module m_gen_nbls

  real hardcut
  integer maxalcsz ! default variable for the total max allocation size

  ! cnblst,cnc_ids : if cmode = 3(4): snapshot NB-list variables
  type t_cnblst
    integer nbs, alsz ! actual number of neighbors and current allocation size
    real, ALLOCATABLE :: dis(:) ! list of distances
    integer, ALLOCATABLE :: idx(:) ! and associated indices
  end type t_cnblst
  ! type t_adjlist
  !   integer deg ! degree of vertex
  !   integer alsz ! allocation size
  !   integer, ALLOCATABLE:: adj(:) ! list of adjacent vertices
  !   real, ALLOCATABLE:: dist(:) ! distance to the adjacent vertices
  ! end type t_adjlist

  type(t_cnblst), ALLOCATABLE :: cnblst(:)
  ! type(t_adjlist), ALLOCATABLE:: approxmst(:)

  contains
    subroutine gen_nb(trj)
      use gutenberg
      use m_clustering
      use m_variables_gen
      implicit none

      real, intent(in) :: trj(n_snaps,n_xyz)
      ! integer oldsz
      real tmp_d ! temporary variable for radiuses and distances
      integer i, ii, k, kk, j, l, ll, mi, mj, u, h1
      integer testcnt, testcnt2, overbounds
      real vecti(n_xyz),vectj(n_xyz)
      real maxx !for dist_methods balancing (maximum values)
      ! call spr("DEBUGGING"
      ! call spr("trj_input",trj(1:10,1:10)
      ! call spr("scluster 1", scluster(2)%snaps
      allocate(cnblst(n_snaps))
      cnblst(:)%nbs = 0 ! number of snapshots that are connected to one snap
      cnblst(:)%alsz = 4 ! allocation size
      ! maxalcsz = 4 ! default variable for the total max allocation size
      do i=1,n_snaps
        allocate(cnblst(i)%idx(cnblst(i)%alsz))
        allocate(cnblst(i)%dis(cnblst(i)%alsz))
      end do
      maxx = 0
      ii = 0
      testcnt = 0
      testcnt2 = 0
      overbounds = 0
      call sl()
      call sipr("Distance method: ", dis_method)
      call sl()
      do i=1,nclu ! for each cluster (intra cluster distances)
        ii = ii + 1
        do kk=1,scluster(i)%nmbrs ! for each snapshot in a cluster
          k = scluster(i)%snaps(kk)
          do ll=kk+1,scluster(i)%nmbrs ! for each l != k snapshot
            l = scluster(i)%snaps(ll)
            if(dis_method .eq. 5 .or. dis_method .eq. 1 &
            & .or. dis_method .eq. 12) then
              vecti = trj(k,1:n_xyz)
              vectj = trj(l,1:n_xyz)
            else if(dis_method.eq.11) then
              if(k.eq.1) then !to avoid error I can double the first value
                 vecti = trj(k+1,n_xyz) - trj(k,n_xyz)
              else
                vecti = trj(k,n_xyz) - trj(k-1,n_xyz)
              end if
              if(l.eq.1) then !to avoid error I can double the first value
                vectj = trj(l+1,n_xyz) - trj(l,n_xyz)
              else
                vectj = trj(l,n_xyz) - trj(l-1,n_xyz)
              end if
            end if
            call distance(tmp_d,vecti,vectj)
            if(abs(tmp_d).gt.maxx) maxx = abs(tmp_d)
            ! if(cnblst(k)%nbs.ge.1499) call spr(tmp_d,l,k
            ! compute the distance between all cluster-internal snapshots
            testcnt = testcnt + 1
            if (tmp_d.lt.hardcut) then ! hardcut ext-var CCUTOFF
              testcnt2 = testcnt2 + 1
              ! a_mat(l,k) = tmp_d
              ! a_mat(k,l) = tmp_d
              cnblst(k)%nbs = cnblst(k)%nbs + 1

              if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(cnblst(k))
              ! if (cnblst(k)%alsz.gt.maxalcsz) maxalcsz = cnblst(k)%alsz
              cnblst(k)%idx(cnblst(k)%nbs) = l !indexes
              cnblst(k)%dis(cnblst(k)%nbs) = tmp_d !distance val
              cnblst(l)%nbs = cnblst(l)%nbs + 1
              if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(cnblst(l))
              cnblst(l)%idx(cnblst(l)%nbs) = k
              cnblst(l)%dis(cnblst(l)%nbs) = tmp_d
            end if
          end do
        end do
        do j=i+1,nclu
          if (scluster(i)%nmbrs.le.scluster(j)%nmbrs) then
            mj = j !mj is major
            mi = i !mi is minor
          else
            mj = i
            mi = j
          end if
          do kk=1,scluster(mi)%nmbrs
            k = scluster(mi)%snaps(kk)
            if(dis_method .eq. 5 .or. dis_method .eq. 1 &
            & .or. dis_method .eq. 12) then
              vecti = trj(k,1:n_xyz)
            else if(dis_method.eq.11) then
              if(k.eq.1) then !to avoid error I can double the first value
                 vecti = trj(k+1,n_xyz) - trj(k,n_xyz)
              else
                vecti = trj(k,n_xyz) - trj(k-1,n_xyz)
              end if
            end if
            call snap_to_cluster_d(tmp_d,scluster(mj),vecti)
            !all the minorcluster snaps vs the mjcluster center
            do ll=1,scluster(mj)%nmbrs
              l = scluster(mj)%snaps(ll)
              if(dis_method .eq. 5 .or. dis_method .eq. 1 &
              & .or. dis_method .eq. 12) then
                vectj = trj(l,1:n_xyz)
              else if(dis_method.eq.11) then
                if(l.eq.1) then !to avoid error I can double the first value
                  vectj = trj(l+1,n_xyz) - trj(l,n_xyz)
                else
                  vectj = trj(l,n_xyz) - trj(l-1,n_xyz)
                end if
              end if
              if (scluster(mj)%nmbrs.eq.1) then
                call cluster_to_cluster_d(tmp_d,scluster(mi),scluster(mj))
                testcnt = testcnt + 1
              else
                call distance(tmp_d,vecti,vectj)
                if(abs(tmp_d).gt.maxx) maxx = abs(tmp_d)
                !then you do a complete graph doing euristic nothing
                testcnt = testcnt + 1
              end if
              if (tmp_d.le.hardcut) then
                testcnt2 = testcnt2 + 1
                cnblst(k)%nbs = cnblst(k)%nbs + 1
                if (cnblst(k)%nbs.gt.cnblst(k)%alsz) call cnbl_resz(cnblst(k))
                cnblst(k)%idx(cnblst(k)%nbs) = l !indexes
                cnblst(k)%dis(cnblst(k)%nbs) = tmp_d !distance val
                cnblst(l)%nbs = cnblst(l)%nbs + 1
                if (cnblst(l)%nbs.gt.cnblst(l)%alsz) call cnbl_resz(cnblst(l))
                cnblst(l)%idx(cnblst(l)%nbs) = k
                cnblst(l)%dis(cnblst(l)%nbs) = tmp_d
              end if
            end do
          end do
        end do
      end do

      ! Normalizing
      if(normalize_dis) then
        call srpr("Normalization mode active. Max value that will be &
        & used to normalize (in abs):", maxx)
        call sl()
        do i=1,n_snaps
          do u=1,n_snaps
            if(dis_method.eq.11) then
              !TO DO: understand this problem
              if(cnblst(i)%dis(u).gt.20.or.cnblst(i)%dis(u).lt.(-20)) then
                cnblst(i)%nbs = cnblst(i)%nbs + 1
                cnblst(i)%dis(u) = maxx
                cnblst(i)%idx(u) = n_snaps
                if(i.eq.n_snaps) cnblst(i)%idx(u) = n_snaps - 1
                overbounds = overbounds + 1
              end if
            end if
          end do
          cnblst(i)%dis = cnblst(i)%dis / abs(maxx)
        end do
      end if
      if(overbounds.gt.0) call sipr('Out of bounds variables:', overbounds)
      ! Adding previous distance if more than one is selected
      call spr('...done')
      call srpr( 'Distances computed (%): ',(100.0*testcnt)/(0.5*n_snaps*(n_snaps-1)))
      call srpr( 'Distances kept after cutoff (%): ',(100.0*testcnt2)/(1.0*testcnt))
    end subroutine gen_nb
end module m_gen_nbls
