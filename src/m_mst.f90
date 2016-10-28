module m_mst
  logical mst


  contains
    !---------------------------------------------------------------------------------
    !
    ! this subroutine generates an exact MST assuming it is provided with a nb-list
    ! object that holds all the necessary edges
    ! this routine is very memory-intensive due to the duplication of the already large nb-list object
    ! it is highly related to hierarchical clustering with minimum linkage and max threshold
    !
    subroutine gen_MST_from_nbl(adjl_deg2,adjl_ix2,adjl_dis2,max_degree)
    !
      use m_clustering
      use m_gen_nbls
      use m_variables_gen
      ! cstored = n_snaps
      ! cnblst
      implicit none

      integer, intent(inout) :: adjl_deg2(n_snaps)
      integer, intent(inout) :: adjl_ix2(n_snaps,n_snaps)
      real, intent(inout) :: adjl_dis2(n_snaps,n_snaps)

      integer max_degree !made for successive cooperation with R wrapper
      integer allnbs !total number of unique connections
      integer i,j,k
      integer globi !global index of numbered links
      integer i_of_lnk,j_of_lnk !snapshot indices for a specific link
      integer ntrees, nlnks !identificatio of the tree number and number of links
      real tmp_dist
      !They are the transpose one of the other. alllinks is also sorted
      integer, ALLOCATABLE :: tmp_all_lnks(:,:), alllnks(:,:)

      !distances and indexes
      real(KIND=4), ALLOCATABLE :: alldiss(:) !vector of all distances
      integer, ALLOCATABLE :: ix(:)

      !for sorting
      real, ALLOCATABLE :: Talldiss(:)
      integer, ALLOCATABLE :: Tix(:)

      integer, ALLOCATABLE :: iv1(:)
      type(t_scluster), ALLOCATABLE :: it(:)
      logical notdone
    !
      write(ilog,*)
      write(ilog,*) 'Now creating global sorted list of neighbor pairs ...'
    !
      notdone = .true. !
      allnbs = sum(cnblst(1:n_snaps)%nbs)/2 ! Total number of unique connections
      allocate(tmp_all_lnks(allnbs,2)) !
      allocate(alldiss(allnbs)) ! NOT
      allocate(Talldiss((allnbs+1)/2))
      allocate(ix(allnbs))!
      allocate(Tix((allnbs+1)/2))!
      allocate(alllnks(2,allnbs))!
        ! allocate(iv2(allnbs,2)) tmp_all_lnks
        ! allocate(alldiss(allnbs)) same
        ! allocate(tmpv(allnbs)) Talldiss
        ! allocate(iv1(allnbs)) ix
        ! allocate(iv3(allnbs)) Tix

      j = 0
      do i=1,n_snaps ! n_snaps
        do k=1,cnblst(i)%nbs
          if (cnblst(i)%idx(k).gt.i) then
            j = j + 1
            tmp_all_lnks(j,1) = i
            tmp_all_lnks(j,2) = cnblst(i)%idx(k)
            alldiss(j) = cnblst(i)%dis(k)
            !i is connected to k if k > i.
            !tmp_all_lnks(j,1) has the i for every connection
            !tmp_all_lnks(j,2) has the k for every connection
            !tmpv(j) has the distance of this connection
            !THEN j is the number of connections
          end if
        end do
      end do

      !deallocate cnblst
      do i=1,n_snaps
        if (allocated(cnblst(i)%dis).EQV..true.) deallocate(cnblst(i)%dis)
        if (allocated(cnblst(i)%idx).EQV..true.) deallocate(cnblst(i)%idx)
      end do
      deallocate(cnblst)

      do j=1,allnbs !indexes to sort
        ix(j) = j
      end do

      ! call merge_sort(ldim=allnbs,up=atrue,list=tmpv(1:allnbs),olist=alldiss(1:allnbs),&
     ! &                ilo=aone,ihi=allnbs,idxmap=iv1(1:allnbs),olist2=ix(1:allnbs))
     !the output of this ordering is iv3 and alldiss <---
     !iv1 is only a seq of 1:allnbs
     !atrue -> mystery
     !ilo (input lower) ihi (input higher) I guess
     !
     !The ordered distances are made on tmpv and output in alldiss while the ordered
     !indexes are output in iv3

     !SIMPLER MERGE_SORT:

        ! integer, parameter :: N = 8 !== allnbs
        ! real, dimension(N) :: A = alldiss
        ! integer, dimension(N) :: ix iv3

      call MergeSort(alldiss,ix,allnbs,Talldiss,Tix)

  !     call merge_sort(ldim=allnbs,up=atrue,list=tmpv(1:allnbs),olist=alldiss(1:allnbs),&
  !  &                ilo=aone,ihi=allnbs,idxmap=iv1(1:allnbs),olist2=iv3(1:allnbs))
  ! allocate(iv2(allnbs,2)) tmp_all_lnks
  ! allocate(alldiss(allnbs)) same
  ! allocate(tmpv(allnbs)) Talldiss
  ! allocate(iv1(allnbs)) ix
  ! allocate(iv3(allnbs)) Tix
      ! 1->allnbs (3)
      ! 3->Talldiss (4)
      ! 4->alldiss (1)
      ! 6->allnbs (3)
      ! 7->ix (2)
      ! 8->Tix (5)
      ! allocate(alldiss(allnbs)) ! NOT
      ! allocate(Talldiss((allnbs+1)/2))
      ! allocate(ix(allnbs))!
      ! allocate(Tix((allnbs+1)/2))!
      ! allocate(alllnks(2,allnbs))!
      do i=1,allnbs
        alllnks(:,i) = tmp_all_lnks(ix(i),:) !using the new order iv3 to order it
    !    write(ilog,*) alllnks(1:2,i),alldiss(i)
      end do
      deallocate(tmp_all_lnks)
      deallocate(ix)
      deallocate(Tix)
      deallocate(Talldiss)
      write(ilog,*) '... done.'
      write(ilog,*)
    !
      write(ilog,*) 'Now generating MST by considering shortest remaining link and merging ...'
      allocate(iv1(n_snaps))
      ! variables initialization
      iv1(:) = 0 !empty index vector
      globi = 1 !global index
      nlnks = 0
      ntrees = 0
      nclalcsz = 10 !number of clusters
      allocate(it(nclalcsz))
      do i=1,nclalcsz
        it(i)%alsz = 0
        it(i)%nmbrs = 0
      end do
      !
      do while (notdone.EQV..true.)
        tmp_dist = alldiss(globi) !dist of connection globi (alldiss is ordered)
        i_of_lnk = alllnks(1,globi) !i
        j_of_lnk = alllnks(2,globi) !k
        globi = globi + 1
        if (globi.gt.allnbs) exit
        if ((iv1(i_of_lnk).le.0).AND.(iv1(j_of_lnk).le.0)) then
          !if i and j are not in any tree
          nlnks = nlnks + 1
          ntrees = ntrees + 1
          iv1(i_of_lnk) = ntrees !add the tree!!
          iv1(j_of_lnk) = ntrees
          if (ntrees.gt.nclalcsz) call cluster_lst_resize(it)
          !adding vertices i and j to tree iv1(i_of_lnk)
          call cluster_addsnap_ix(it(iv1(i_of_lnk)),i_of_lnk)
          call cluster_addsnap_ix(it(iv1(i_of_lnk)),j_of_lnk)
        else if ((iv1(i_of_lnk).gt.0).AND.(iv1(j_of_lnk).gt.0)) then !i,j are in trees
          if (iv1(i_of_lnk).eq.iv1(j_of_lnk)) cycle !they belong to the same tree!!
          nlnks = nlnks + 1
          if (it(iv1(i_of_lnk))%nmbrs.gt.it(iv1(j_of_lnk))%nmbrs) then
            !if they are not equal but both are != 0 -> join the cluster (smaller)
            j = iv1(j_of_lnk)
            do i=1,it(j)%nmbrs
              !make all the elements of the smaller clu to belong to the big tree
              iv1(it(j)%snaps(i)) = iv1(i_of_lnk)
            end do
            call join_clusters(it(iv1(i_of_lnk)),it(j)) !'physically' join them
          else
            j = iv1(i_of_lnk)
            do i=1,it(j)%nmbrs
              iv1(it(j)%snaps(i)) = iv1(j_of_lnk)
            end do
            call join_clusters(it(iv1(j_of_lnk)),it(j))
          end if
        else if (iv1(i_of_lnk).gt.0) then !so the other is 0 not belonging to any tree
          nlnks = nlnks + 1
          iv1(j_of_lnk) = iv1(i_of_lnk) !they belong to the same tree
          call cluster_addsnap_ix(it(iv1(i_of_lnk)),j_of_lnk)
        else
          nlnks = nlnks + 1
          iv1(i_of_lnk) = iv1(j_of_lnk)
          call cluster_addsnap_ix(it(iv1(j_of_lnk)),i_of_lnk)
        end if
        alllnks(1,nlnks) = i_of_lnk
        alllnks(2,nlnks) = j_of_lnk
        alldiss(nlnks) = tmp_dist
        if (nlnks.eq.(n_snaps-1)) exit
      end do
    !
      if (nlnks.ne.(n_snaps-1)) then
        write(ilog,*) 'Fatal. Neighbor list is insufficient to create minimum spanning tree. &
     &Increase relevant thresholds.'
      end if
    !
      deallocate(iv1)
      do i=1,nclalcsz
        if (allocated(it(i)%snaps).EQV..true.) deallocate(it(i)%snaps)
        if (allocated(it(i)%sums).EQV..true.) deallocate(it(i)%sums)
      end do
      deallocate(it)
      call gen_MST(adjl_deg2,adjl_ix2,adjl_dis2,&
      &alllnks(:,1:(n_snaps-1)),alldiss(1:(n_snaps-1)),max_degree)
      deallocate(alldiss)
      deallocate(alllnks)

      write(ilog,*) '... done.'
      write(ilog,*)
    !
  end subroutine gen_MST_from_nbl
    !--------------------------------------------------------------------------------
    !
    ! this routines transcribes a list of edges with their distances, that effectively
    ! describe the MST, into an array of adjacency list objects
    !
    subroutine gen_MST(adjl_deg3,adjl_ix3,adjl_dis3,mstedges,lmstedges, maxdeg)
    !
      use m_gen_nbls
      use m_variables_gen
    !
      implicit none
    !
      integer, intent(inout) :: adjl_deg3(n_snaps)
      integer, intent(inout) :: adjl_ix3(n_snaps,n_snaps)
      real, intent(inout) :: adjl_dis3(n_snaps,n_snaps)
      integer maxdeg
      integer e,v
      integer mstedges(2,n_snaps-1) !mstedges(v,e) are the indexes for an edge
      REAL(kind=4) lmstedges(n_snaps-1)


      !transform edgelist of MST (mstedges) to adjacencylist:
      do e=1,n_snaps-1
        do v=1,2 !for each vertex(v) belonging to each edges (e)
          adjl_deg3(mstedges(v,e)) = adjl_deg3(mstedges(v,e)) + 1
          adjl_ix3(mstedges(v,e),adjl_deg3(mstedges(v,e))) = mstedges(3-v,e)
          adjl_dis3(mstedges(v,e),adjl_deg3(mstedges(v,e))) = lmstedges(e)
          if(e.eq.1.and.v.eq.1) then
            maxdeg = adjl_deg3(mstedges(v,e))
          else if(maxdeg.lt.adjl_deg3(mstedges(v,e))) then
            maxdeg = adjl_deg3(mstedges(v,e))
          end if
        end do
      end do

  end subroutine gen_MST
    !------------------------------------------------------------------------------------------------------
    !
    ! a simple helper to grow as conservatively (and slowly) as possible an adjacency list object
    !
    ! subroutine extend_adjlst_byone(mstnode)
    ! !
    !   use m_gen_nbls
    ! !
    !   implicit none
    !   integer size
    !   integer, ALLOCATABLE:: itmp1(:)
    !   real(KIND=4), ALLOCATABLE:: rtmp1(:)
    !   type(t_adjlist) mstnode
    ! !
    !   allocate(itmp1(mstnode%deg))
    !   allocate(rtmp1(mstnode%deg))
    ! !
    !   itmp1(1:mstnode%deg) = mstnode%adj(1:mstnode%deg)
    !   rtmp1(1:mstnode%deg) = mstnode%dist(1:mstnode%deg)
    !   deallocate(mstnode%adj)
    !   deallocate(mstnode%dist)
    !   size = mstnode%deg*2
    !
    !   ! size = mstnode%deg*2
    !   ! if(size .gt. n_snaps-1) then
    !   !   size = n_snaps-1
    !   !   if(verbose) write(ilog,*) "Size limit reached in deg"
    !   ! end if
    !   allocate(mstnode%adj(size))
    !   allocate(mstnode%dist(size))
    !   mstnode%adj(1:mstnode%deg) = itmp1(1:mstnode%deg)
    !   mstnode%dist(1:mstnode%deg) = rtmp1(1:mstnode%deg)
    ! !
    !   deallocate(itmp1)
    !   deallocate(rtmp1)
    ! !
    ! end
end module m_mst
