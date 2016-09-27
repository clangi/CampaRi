module m_gen_nbls

  real(KIND=4) hardcut
  integer maxalcsz ! default variable for the total max allocation size

  ! cnblst,cnc_ids : if cmode = 3(4): snapshot NB-list variables
    type t_cnblst
      integer nbs,alsz ! actual number of neighbors and current allocation size
      real(KIND=4), ALLOCATABLE:: dis(:) ! list of distances
      integer, ALLOCATABLE:: idx(:) ! and associated indices
    end type t_cnblst
    ! type t_adjlist
    !   integer deg ! degree of vertex
    !   integer alsz ! allocation size
    !   integer, ALLOCATABLE:: adj(:) ! list of adjacent vertices
    !   real(KIND=4), ALLOCATABLE:: dist(:) ! distance to the adjacent vertices
    ! end type t_adjlist

    type(t_cnblst), ALLOCATABLE :: cnblst(:)
    ! type(t_adjlist), ALLOCATABLE:: approxmst(:)

  contains
    subroutine gen_nb(trj)
      use m_clustering
      use m_var_nbls_clu
      implicit none

      real(KIND=4), intent(in) :: trj(n_snaps,n_xyz)
      ! integer oldsz
      real(KIND=4) tmp_d ! temporary variable for radiuses and distances
      integer i, ii, k, kk, j, l, ll, mi, mj, u
      integer testcnt, testcnt2
      real vecti(n_xyz),vectj(n_xyz)
      ! write(*,*) "DEBUGGING"
      ! write(*,*) "trj_input",trj(1:10,1:10)
      ! write(*,*) "scluster 1", scluster(2)%snaps
      write(*,*)
      write(*,*)

      ii = 0
      testcnt = 0
      testcnt2 = 0
      do i=1,nclu ! for each cluster
        ii = ii + 1
        do kk=1,scluster(i)%nmbrs ! for each snapshot in a cluster
          k = scluster(i)%snaps(kk)
          do ll=kk+1,scluster(i)%nmbrs ! for each l != k snapshot
            l = scluster(i)%snaps(ll)
            if(dis_method.eq.5) then
              vecti = trj(k,1:n_xyz)
              vectj = trj(l,1:n_xyz)
            else if(dis_method.eq.11) then
              ! veci(1:n_xyz) = 1
              if(k.eq.1) k = 2 !to avoid error I can double the first value
              if(l.eq.1) l = 2 !to avoid error I can double the first value
              vecti = trj(k,n_xyz) - trj(k-1,n_xyz)
              vectj = trj(l,n_xyz) - trj(l-1,n_xyz)
            end if
            call distance(tmp_d,vecti,vectj)
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
            if(dis_method.eq.5) then
              vecti = trj(k,1:n_xyz)
            else if(dis_method.eq.11) then
              if(k.eq.1) k = 2 !to avoid error I can double the first value
              vecti = trj(k,n_xyz) - trj(k-1,n_xyz)
            end if
            call snap_to_cluster_d(tmp_d,scluster(mj),vecti)
            !all the minorcluster snaps vs the mjcluster center
            do ll=1,scluster(mj)%nmbrs
              l = scluster(mj)%snaps(ll)
              if(dis_method.eq.5) then
                vectj = trj(l,1:n_xyz)
              else if(dis_method.eq.11) then
                if(l.eq.1) l = 2 !to avoid error I can double the first value
                vectj = trj(l,n_xyz) - trj(l-1,n_xyz)
              end if
              if (scluster(mj)%nmbrs.eq.1) then
                call cluster_to_cluster_d(tmp_d,scluster(mi),scluster(mj))
                testcnt = testcnt + 1
              else
                call distance(tmp_d,vecti,vectj)
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
      write(*,*) '... done after computing ',(100.0*testcnt)/(0.5*n_snaps*(n_snaps-1)),'% of &
      &possible terms with ',(100.0*testcnt2)/(1.0*testcnt),'% successful.'
    end subroutine gen_nb
end module m_gen_nbls
