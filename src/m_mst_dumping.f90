module m_mst_dumping
  type t_progindextree
    integer nsnaps ! number of snapshots in tree
    integer nsiblings ! number of trees to merge with
    integer nsibalsz ! alloc size for that
    integer, ALLOCATABLE:: snaps(:) ! indices of snapshots in tree
    integer, ALLOCATABLE:: siblings(:) ! tree indices of tree to merge with
    integer mine(2) ! shortest edge leaving the tree
    real mind ! distance to nearest snapshot of tree, equals length(mine)
    integer ptr ! simple pointer
  end type t_progindextree
  type(t_progindextree), ALLOCATABLE:: tmptree(:)
  integer cdevalcnt !faith mystery (some kind of count)
  integer cprogbatchsz ! cprogbatchsz       : if cmode = 4 and cprogindex = 2: batch size for random stretches aka dim of random branches
  integer cprogrdepth ! cprogrdepth        : if cmode = 4 and cprogindex = 2: auxiliary search depth
  integer cprogindrmax ! this must be an INPUT MAIN: it is the number of guesses

  !this will be used in the dumping
  type t_adjlist
    integer deg ! degree of vertex
    integer alsz ! allocation size
    integer, ALLOCATABLE :: adj(:) ! list of adjacent vertices
    real, ALLOCATABLE :: dist(:) ! distance to the adjacent vertices
    ! logical, ALLOCATABLE:: tagged(:) ! helper flag for reading netcdf
  end type t_adjlist
  type(t_adjlist), allocatable :: approxmst(:) !output of everything



  contains
    !---------------------------------------------------------------------------------
    !
    ! this subroutine generates an exact MST assuming it is provided with a nb-list
    ! object that holds all the necessary edges
    ! this routine is very memory-intensive due to the duplication of the already large nb-list object
    ! it is highly related to hierarchical clustering with minimum linkage and max threshold
    !
    subroutine gen_MST_from_nbl_w()
    !
      use m_clustering
      use m_gen_nbls
      use m_variables_gen

      implicit none

      integer allnbs !total number of unique connections
      integer i,j,k
      integer globi !global index of numbered links
      integer i_of_lnk, j_of_lnk !snapshot indices for a specific link
      integer ntrees, nlnks !identificatio of the tree number and number of links
      real tmp_dist
      !They are the transpose one of the other. alllinks is also sorted
      integer, ALLOCATABLE :: tmp_all_lnks(:,:), alllnks(:,:)

      !distances and indexes
      real, ALLOCATABLE :: alldiss(:) !vector of all distances
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

      call MergeSort(alldiss,ix,allnbs,Talldiss,Tix,.true.)

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
      n_clu_alc_sz_gen = 10 !number of clusters
      allocate(it(n_clu_alc_sz_gen))
      do i=1,n_clu_alc_sz_gen
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
          if (ntrees.gt.n_clu_alc_sz_gen) call cluster_lst_resize(it)
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
      do i=1,n_clu_alc_sz_gen
        if (allocated(it(i)%snaps).EQV..true.) deallocate(it(i)%snaps)
        if (allocated(it(i)%sums).EQV..true.) deallocate(it(i)%sums)
        if (allocated(it(i)%tmpsnaps).EQV..true.) deallocate(it(i)%tmpsnaps)
        if (allocated(it(i)%children).EQV..true.) deallocate(it(i)%children)
      end do
      deallocate(it)
      allocate(approxmst(n_snaps))
      approxmst(1:n_snaps)%deg = 0
      call gen_MST_w(alllnks(:,1:(n_snaps-1)),alldiss(1:(n_snaps-1)),approxmst)
      deallocate(alldiss)
      deallocate(alllnks)

      write(ilog,*) '... done.'
      write(ilog,*)
    !
    end subroutine gen_MST_from_nbl_w
    !--------------------------------------------------------------------------------
    !
    ! this routines transcribes a list of edges with their distances, that effectively
    ! describe the MST, into an array of adjacency list objects
    !
    subroutine gen_MST_w(mstedges,lmstedges,mst)
    !
      use m_variables_gen
    !
      implicit none
    !
      integer e,v
      type(t_adjlist) mst(n_snaps)
      integer mstedges(2,n_snaps-1)
      real lmstedges(n_snaps-1)
    !
    ! transform edgelist of MST (mstedges) to adjacencylist:
      do e=1,n_snaps
        allocate(mst(e)%adj(1))
        allocate(mst(e)%dist(1))
      end do
      do e=1,n_snaps-1
        do v=1,2
          if (mst(mstedges(v,e))%deg.gt.0) call extend_adjlst_byone(mst(mstedges(v,e)))
          mst(mstedges(v,e))%deg = mst(mstedges(v,e))%deg + 1
          mst(mstedges(v,e))%adj(mst(mstedges(v,e))%deg) = mstedges(3-v,e)
          mst(mstedges(v,e))%dist(mst(mstedges(v,e))%deg) = lmstedges(e)
        end do
      end do
    !
    end
    !------------------------------------------------------------------------------------------------------
    !
    ! a simple helper to grow as conservatively (and slowly) as possible an adjacency list object
    !
    subroutine extend_adjlst_byone(mstnode)
    !
      implicit none
    !
      integer, ALLOCATABLE:: itmp1(:)
      real, ALLOCATABLE:: rtmp1(:)
      type(t_adjlist) mstnode
    !
      allocate(itmp1(mstnode%deg))
      allocate(rtmp1(mstnode%deg))
    !
      itmp1(1:mstnode%deg) = mstnode%adj(1:mstnode%deg)
      rtmp1(1:mstnode%deg) = mstnode%dist(1:mstnode%deg)
      deallocate(mstnode%adj)
      deallocate(mstnode%dist)
    !
      allocate(mstnode%adj(mstnode%deg+1))
      allocate(mstnode%dist(mstnode%deg+1))
      mstnode%adj(1:mstnode%deg) = itmp1(1:mstnode%deg)
      mstnode%dist(1:mstnode%deg) = rtmp1(1:mstnode%deg)
    !
      deallocate(itmp1)
      deallocate(rtmp1)
    !
    end
    !-----------------------------------------------------------------------------------------
    !
    ! this subroutine generates a set of links and their lengths that constitute an approximate
    ! minimum spanning tree (SST) based on results from tree-based clustering in the birchtree object
    ! it then transcribes this into an adjacency list object
    ! the accuracy of the SST depends on the number of guesses (cprogindrmax), the properties
    ! of the clustering, and the auxiliary search depth being utilized (cprogrdepth)
    !
    subroutine gen_MST_from_treeclustering_w(trj2)
    !
      use m_clustering
      use m_variables_gen
      ! use iounit
      ! use interfaces
    !
      implicit none
    !
      integer N_NEARS
      parameter(N_NEARS=5)  ! number of guesses per snapshot to keep in memory

      !my adds
      real, intent(in) :: trj2(n_snaps,n_xyz)
      integer max_degree !made for successive cooperation with R wrapper

      integer i,j,k,ki,e,i1,i2,l,m,mm,ixx,ixx2,cursnp,thismode,b1,b2,mycl,cbnds(2)
      integer tlstsz,ishfx,kshfx
      integer(8) testcnt,trcnt,tdcnt
      integer ntrees,oldntrees ! number of active trees
      integer, ALLOCATABLE:: tmpix(:,:)
      real, ALLOCATABLE:: tmpdis(:,:)
      real tmvars(2),tmpcld(n_xyz) !calcsz doesn't exist nomore -> n_xyz
      ! real random ! random number
      real random_or
      real jkdist ! distance between two snapshots
      integer boruvkasteps! counts how many boruvka steps are needed
      integer nmstedges   ! current number of edges in the growing trees
      integer, ALLOCATABLE :: testlst(:,:) ! temporary list of eligible candidates for random search with annotations
      integer, ALLOCATABLE :: to_order(:,:)
      integer, ALLOCATABLE :: snplst(:)    ! a temp array to hold lists of snapshots to scan for eligible distances
      integer a           ! numbers of snapshots, used for sophisticated nearest neighbor guessing
      integer kk,kix,kixi ! local snapshot index within cluster, used for sophisticated nearest neighbor guessing
      integer, ALLOCATABLE :: satisfied(:) ! array to indicate whether enough distances have been probed
      integer, ALLOCATABLE :: mstedges(:,:)! contains all the edges in the SST
                                          ! Indices: edgenumber (1-nedges), endpoint (1-2)
      real, ALLOCATABLE:: lmstedges(:) ! contains the lengths of all the edges in the SST


      !myvars
      integer csnap2tree_ik_eq_cnt !counts of the "should-be" exits
      real t4,t3

      !from mod_clusters
      integer, ALLOCATABLE:: csnap2tree(:),csnap2clus(:,:)
      real vecti(n_xyz)


      !myvars defaults
      csnap2tree_ik_eq_cnt = 0


      ntrees = n_snaps
      boruvkasteps = 0
      tlstsz = n_snaps
      testcnt = 0
      trcnt = 0
      tdcnt = 0
      nmstedges = 0
    !
      allocate(tmptree(ntrees))
      allocate(csnap2tree(n_snaps))
      allocate(csnap2clus(n_snaps,c_nhier+1))
      allocate(mstedges(2,n_snaps-1))
      allocate(lmstedges(n_snaps-1))
    !
      write(ilog,*)
      write(ilog,*) 'Now generating approximate MST (SST) based on tree-based clustering ...'
      do i=1,ntrees !ntrees is initialized with number of snapshots (cstored)
        allocate(tmptree(i)%snaps(1)) !there are n_snaps trees allocated. For each the first snap is allocated
        tmptree(i)%snaps(1) = i
        csnap2tree(i) = i
      end do
      tmptree(:)%nsnaps = 1
      tmptree(:)%mine(1) = -1
      tmptree(:)%mine(2) = -1
      tmptree(:)%mind = HUGE(jkdist)
      tmptree(:)%nsibalsz = 0
    !
    ! Generate map (snapshot, tree level) -> cluster based on birchtree;
    ! in the process remove double entries for non-terminal levels and populate persistent tree index vector per cluster
      do l=c_nhier+1,1,-1
        csnap2clus(:,l) = 0
        do k=1,birchtree(l)%ncls
          if (allocated(birchtree(l)%cls(k)%tmpsnaps).EQV..true.) deallocate(birchtree(l)%cls(k)%tmpsnaps)
          allocate(birchtree(l)%cls(k)%tmpsnaps(birchtree(l)%cls(k)%nmbrs))
          do j=1,birchtree(l)%cls(k)%nmbrs
            birchtree(l)%cls(k)%tmpsnaps(j) = csnap2tree(birchtree(l)%cls(k)%snaps(j))
            if (csnap2clus(birchtree(l)%cls(k)%snaps(j),l).gt.0) then
              if (birchtree(l)%cls(csnap2clus(birchtree(l)%cls(k)%snaps(j),l))%nmbrs.gt.birchtree(l)%cls(k)%nmbrs) then
                call preprocess_snapshot(trj2,birchtree(l)%cls(k)%snaps(j),vecti)
                call cluster_removesnap(birchtree(l)%cls(k),birchtree(l)%cls(k)%snaps(j),vecti)
              else
                call preprocess_snapshot(trj2,birchtree(l)%cls(k)%snaps(j),vecti)
                call cluster_removesnap(birchtree(l)%cls(csnap2clus(birchtree(l)%cls(k)%snaps(j),l)),&
                  birchtree(l)%cls(k)%snaps(j),vecti)
                csnap2clus(birchtree(l)%cls(k)%snaps(j),l) = k !because the k cls has more elements than the previous stored a new maximum is assigned
              end if
            else
              csnap2clus(birchtree(l)%cls(k)%snaps(j),l) = k !k has always the maximum number of snaps
            end if
          end do
        end do
      end do
    !
    ! temporary working variables
      allocate(satisfied(tlstsz))
      satisfied(:) = 0
      allocate(snplst(tlstsz))
      snplst(:) = c_nhier + 1
      allocate(tmpix(N_NEARS,tlstsz))
      allocate(tmpdis(N_NEARS,tlstsz)) ! N_NEARS+1))
      tmpdis(:,:) = HUGE(tmpdis(1,1)) ! initialize
      tmpix(:,:) = 0                  ! initialize
      allocate(testlst(n_snaps,4))
      allocate(to_order(n_snaps,2))
      testlst(:,4) = 0
      ! cbnds(:) = 0 !myadd
      ! thismode = 1 !myadd
    !
      call CPU_time(tmvars(1))
    !
    ! allocation and initialization for temporary working variables
    !
      do while (ntrees.ge.2)
    !
        if (boruvkasteps.gt.0) then
    !     sort snapshots in clusters according to tree membership (this may destroy the integrity of some aspects of the tree structure
    !      -> don't use these aspects hereafter) and update the persistent sorted list itself
          do l=c_nhier+1,1,-1
            do i=1,birchtree(l)%ncls
              if (birchtree(l)%cls(i)%nmbrs.le.1) cycle
              do j=1,birchtree(l)%cls(i)%nmbrs
                testlst(j,4) = csnap2tree(birchtree(l)%cls(i)%snaps(j))
                testlst(j,1) = j
              end do
              ixx = 1
              ixx2 = birchtree(l)%cls(i)%nmbrs
              mm = birchtree(l)%cls(i)%nmbrs
              to_order(1:mm,1) = testlst(1:mm,4)
              to_order(1:mm,2) = testlst(1:mm,1)
              call MergeSort(to_order(1:mm,1),to_order(1:mm,2),mm,testlst(1:mm,3),testlst(1:mm,2),.true.)
    !           call merge_sort(ldim=mm,up=atrue,list=testlst(1:mm,4),olist=testlst(1:mm,3),ilo=ixx,ihi=ixx2,&
    !  &                        idxmap=testlst(1:mm,1),olist2=testlst(1:mm,2))
    !          write(*,*) l,i,':'
    !          write(*,555) testlst(1:mm,3)
              do j=1,mm
                testlst(j,4) = birchtree(l)%cls(i)%snaps(to_order(j,2))
              end do
              birchtree(l)%cls(i)%snaps(1:mm) = testlst(1:mm,4)
              birchtree(l)%cls(i)%tmpsnaps(1:mm) = to_order(1:mm,1)
    !          write(*,555) birchtree(l)%cls(i)%snaps(1:mm) ! testlst(1:mm,3)
    !          write(*,555) snap2tree(birchtree(l)%cls(i)%snaps(1:mm))
            end do
          end do
        end if
    !
    !   first get a set of guesses based on fixed trees, clustering data structure, and prior guesses
        do i=1,n_snaps

          ishfx = i
          l = snplst(ishfx) ! this is remembered across Boruvka stages
          ! write(ilog,*) "asdadsdads", snplst(ishfx)
          do while (satisfied(ishfx).lt.cprogindrmax)
    !
            mycl = csnap2clus(i,l)
            tmpcld(1:n_xyz) = trj2(i,1:n_xyz)
            a = birchtree(l)%cls(mycl)%nmbrs
            ixx = a
    !       we can quickly test whether search space is large enough to trigger random method
            ixx2 = 0
            i1 = -1
            i2 = -1
            if ((cprogindrmax-satisfied(ishfx)).gt.a) then
              thismode = 1 ! determ required
              cbnds(:) = 0
            else if (boruvkasteps.gt.0) then
    !          if (csnap2tree(birchtree(l)%cls(mycl)%snaps(ixx)).lt.csnap2tree(i)) then
    !            i1 = ixx
    !            if (ixx.lt.a) then
    !              if (csnap2tree(birchtree(l)%cls(mycl)%snaps(ixx+1)).lt.csnap2tree(i)) i1 = ixx + 1
    !            end if
    !          else if (csnap2tree(birchtree(l)%cls(mycl)%snaps(1)).lt.csnap2tree(i)) then
              if (csnap2tree(birchtree(l)%cls(mycl)%snaps(1)).lt.csnap2tree(i)) then
    !            testlst(1:ixx,4) = csnap2tree(birchtree(l)%cls(mycl)%snaps(1:ixx))
                call ibinary_search(ixx,birchtree(l)%cls(mycl)%tmpsnaps(1:ixx),csnap2tree(i)-1,ixx2)
                cbnds(1) = ixx2
                i1 = ixx2
              else
                i1 = 0
                cbnds(1) = a
              end if
    !          if (csnap2tree(birchtree(l)%cls(mycl)%snaps(a-ixx+1)).gt.csnap2tree(i)) then
    !            i2 = ixx
    !            if ((a-ixx+1).gt.1) then
    !              if (csnap2tree(birchtree(l)%cls(mycl)%snaps(a-ixx)).gt.csnap2tree(i)) i2 = ixx + 1
    !            end if
    !          else if (csnap2tree(birchtree(l)%cls(mycl)%snaps(a)).gt.csnap2tree(i)) then
              if (csnap2tree(birchtree(l)%cls(mycl)%snaps(ixx)).gt.csnap2tree(i)) then
    !            testlst(1:ixx,4) = csnap2tree(birchtree(l)%cls(mycl)%snaps(1:ixx))
                call ibinary_search(ixx,birchtree(l)%cls(mycl)%tmpsnaps(1:ixx),csnap2tree(i),ixx2)
                cbnds(2) = ixx2+1
                i2 = ixx - ixx2
              else
                i2 = 0
                cbnds(2) = 1
              end if
              if ((i1+i2).gt.(cprogindrmax-satisfied(ishfx))) then
                thismode = 0
              else
                thismode = 1
              end if
            else if (boruvkasteps.eq.0) then
              cbnds(:) = 0
              if ((a-1).gt.(cprogindrmax-satisfied(ishfx))) then
                thismode = 0
              else
                thismode = 1
              end if
            end if
            if (thismode.eq.1) then ! deterministic
              do j=1,a
                k = birchtree(l)%cls(mycl)%snaps(j)
                if (csnap2tree(k).eq.csnap2tree(i)) exit
                testcnt = testcnt + 1
                tdcnt = tdcnt + 1
                call distance(jkdist,tmpcld(1:n_xyz),trj2(k,1:n_xyz))
    !   update shortest eligible edges per snap
                if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,ishfx).eq.k) exit
                    if (jkdist.lt.tmpdis(ixx,ishfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                        tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                      end do
                      tmpix(ixx,ishfx) = k
                      tmpdis(ixx,ishfx) = jkdist
                      exit
                    end if
                  end do
                end if
                kshfx = k
                if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,kshfx).eq.i) exit
                    if (jkdist.lt.tmpdis(ixx,kshfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                        tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                      end do
                      tmpix(ixx,kshfx) = i
                      tmpdis(ixx,kshfx) = jkdist
                      exit
                    end if
                  end do
                end if
                satisfied(ishfx) = satisfied(ishfx) + 1
              end do
              do j=a,1,-1
                k = birchtree(l)%cls(mycl)%snaps(j)
                if (csnap2tree(k).eq.csnap2tree(i)) exit
                tdcnt = tdcnt + 1
                testcnt = testcnt + 1
                call distance(jkdist,tmpcld(1:n_xyz),trj2(k,1:n_xyz))
    !   update shortest eligible edges per snap
                if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,ishfx).eq.k) exit
                    if (jkdist.lt.tmpdis(ixx,ishfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                        tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                      end do
                      tmpix(ixx,ishfx) = k
                      tmpdis(ixx,ishfx) = jkdist
                      exit
                    end if
                  end do
                end if
                kshfx = k
                if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
                  do ixx=1,N_NEARS
                    if (tmpix(ixx,kshfx).eq.i) exit
                    if (jkdist.lt.tmpdis(ixx,kshfx)) then
                      do ixx2=N_NEARS,ixx+1,-1 ! shift
                        tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                        tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                      end do
                      tmpix(ixx,kshfx) = i
                      tmpdis(ixx,kshfx) = jkdist
                      exit
                    end if
                  end do
                end if
                satisfied(ishfx) = satisfied(ishfx) + 1
              end do
            else ! random
              m = a !a is the numbers of elements
              b1 = 0 ! max(1,a/2)
              b2 = 0 ! b1
              if (cbnds(1).gt.0) then
                if (cbnds(2).eq.1) then
                  b2 = 1
                  m = cbnds(1)
                else if (cbnds(1).gt.cbnds(2)) then
                  b2 = cbnds(2)
                  m = cbnds(1) - cbnds(2) + 1
                else if (cbnds(2).gt.cbnds(1)) then
                  b2 = cbnds(2)
                  m = a - (cbnds(2) - cbnds(1) - 1)
                end if
              end if
              if (boruvkasteps.gt.0) then
                do while (satisfied(ishfx).lt.cprogindrmax)
                  kk = int(random_or()*m) + b2
                  do ki=1,cprogbatchsz
                    kix = kk
                    if (kix.gt.a) kix = kix - a
                    if (kix.le.0) kix = kix + a
                    k = birchtree(l)%cls(mycl)%snaps(kix)
                    if (csnap2tree(i).eq.csnap2tree(k)) then
                      csnap2tree_ik_eq_cnt = csnap2tree_ik_eq_cnt + 1
                      !TODO
                      if(csnap2tree_ik_eq_cnt.gt.100) then
                        write(ilog,*) "FATAL : Something went wrong. &
                        &csnap2tree has similar i-k. &
                        &Please consider changing tree height"
                        ! cycle
                        call fexit()
                      end if
                    end if
                    trcnt = trcnt + 1
                    testcnt = testcnt + 1
                    call distance(jkdist,tmpcld(1:n_xyz),trj2(k,1:n_xyz))
    !                update shortest eligible edges per snap
                    if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
                      do ixx=1,N_NEARS
                        if (tmpix(ixx,ishfx).eq.k) exit
                        if (jkdist.lt.tmpdis(ixx,ishfx)) then
                          do ixx2=N_NEARS,ixx+1,-1 ! shift
                            tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                            tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                          end do
                          tmpix(ixx,ishfx) = k
                          tmpdis(ixx,ishfx) = jkdist
                          exit
                        end if
                      end do
                    end if
                    kshfx = k
                    if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
                      do ixx=1,N_NEARS
                        if (tmpix(ixx,kshfx).eq.i) exit
                        if (jkdist.lt.tmpdis(ixx,kshfx)) then
                          do ixx2=N_NEARS,ixx+1,-1 ! shift
                            tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                            tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                          end do
                          tmpix(ixx,kshfx) = i
                          tmpdis(ixx,kshfx) = jkdist
                          exit
                        end if
                      end do
                    end if
                    kk = kk + 1
                    if (kk.gt.(m+b2-1)) kk = b2
                    satisfied(ishfx) = satisfied(ishfx) + 1
                    if (satisfied(ishfx).eq.cprogindrmax) exit
                  end do
                end do
              else !boruvkasteps == 0
                kixi = max(1,a/cprogbatchsz)
                do while (satisfied(ishfx).lt.cprogindrmax)
                  kk = int(random_or()*a)
                  do ki=1,cprogbatchsz
                    kix = kk
                    if (kix.gt.a) kix = kix - a
                    if (kix.le.0) kix = kix + a
                    k = birchtree(l)%cls(mycl)%snaps(kix)
                    if (csnap2tree(i).eq.csnap2tree(k)) then
                      cycle ! in first stage, only a single snapshot per cluster needs to be cycled
                    end if
                    call distance(jkdist,tmpcld(1:n_xyz),trj2(k,1:n_xyz))
                    trcnt = trcnt + 1
                    testcnt = testcnt + 1
    !   update shortest eligible edges per snap
                    if (jkdist.lt.tmpdis(N_NEARS,ishfx)) then
                      do ixx=1,N_NEARS
                        if (tmpix(ixx,ishfx).eq.k) exit
                        if (jkdist.lt.tmpdis(ixx,ishfx)) then
                          do ixx2=N_NEARS,ixx+1,-1 ! shift
                            tmpix(ixx2,ishfx) = tmpix(ixx2-1,ishfx)
                            tmpdis(ixx2,ishfx) = tmpdis(ixx2-1,ishfx)
                          end do
                          tmpix(ixx,ishfx) = k
                          tmpdis(ixx,ishfx) = jkdist
                          exit
                        end if
                      end do
                    end if
                    kshfx = k
                    if (jkdist.lt.tmpdis(N_NEARS,kshfx)) then
                      do ixx=1,N_NEARS
                        if (tmpix(ixx,kshfx).eq.i) exit
                        if (jkdist.lt.tmpdis(ixx,kshfx)) then
                          do ixx2=N_NEARS,ixx+1,-1 ! shift
                            tmpix(ixx2,kshfx) = tmpix(ixx2-1,kshfx)
                            tmpdis(ixx2,kshfx) = tmpdis(ixx2-1,kshfx)
                          end do
                          tmpix(ixx,kshfx) = i
                          tmpdis(ixx,kshfx) = jkdist
                          exit
                        end if
                      end do
                    end if
                    kk = kk + kixi
                    if (kk.gt.a) kk = kk - a
                    satisfied(ishfx) = satisfied(ishfx) + 1
                    if (satisfied(ishfx).eq.cprogindrmax) exit
                  end do
                end do
              end if
            end if
            l = l - 1
            if (satisfied(ishfx).eq.0) snplst(ishfx) = l
            if ((snplst(ishfx)-(l+1)).eq.cprogrdepth) exit
            if (l.eq.0) exit
          end do
        end do
    !
        testlst(1:ntrees,3) = 0
    !
        do cursnp=1,n_snaps
          ishfx = cursnp
          j = csnap2tree(cursnp)
          if (testlst(j,3).eq.0) then ! first edge for that tree
            testlst(j,1) = cursnp
            testlst(j,2) = tmpix(1,ishfx)
            testlst(j,3) = ishfx
            testlst(j,4) = 1
          else if (tmpdis(1,ishfx).lt.tmpdis(testlst(j,4),testlst(j,3))) then
            testlst(j,1) = cursnp
            testlst(j,2) = tmpix(1,ishfx)
            testlst(j,3) = ishfx
            testlst(j,4) = 1
          end if
          do ixx=1,N_NEARS
            if (tmpix(ixx,ishfx).le.0) cycle
            k = csnap2tree(tmpix(ixx,ishfx))
            if (testlst(k,3).eq.0) then ! first edge for that tree
              testlst(k,1) = tmpix(ixx,ishfx)
              testlst(k,2) = cursnp
              testlst(k,3) = ishfx
              testlst(k,4) = ixx
            else if (tmpdis(ixx,ishfx).lt.tmpdis(testlst(k,4),testlst(k,3))) then
              testlst(k,1) = tmpix(ixx,ishfx)
              testlst(k,2) = cursnp
              testlst(k,3) = ishfx
              testlst(k,4) = ixx
            end if
          end do
        end do
        do k=1,ntrees ! could be made faster for first stage at least
          if (testlst(k,3).gt.0) then
            if (tmpdis(testlst(k,4),testlst(k,3)).lt.tmptree(k)%mind) then
              tmptree(k)%mind = tmpdis(testlst(k,4),testlst(k,3))
              tmptree(k)%mine(1) = testlst(k,1)
              tmptree(k)%mine(2) = testlst(k,2)
            end if
          end if
        end do
    !
    ! merge trees (& update approximate minimum spanning tree):
        tmptree(1:ntrees)%nsiblings = 0
        do k=1,ntrees
          tmptree(k)%ptr = k
        end do
    !
        do e=1,ntrees ! nedges
          i1 = tmptree(e)%ptr                               ! tree index of first edge endpoint
          i2 = tmptree(csnap2tree(tmptree(e)%mine(2)))%ptr  ! tree index of second edge endpoint
          if (i1.ne.i2) then ! edge does not introduce cycle
    ! add corresponding edge to approximate minimum spanning tree:
            nmstedges = nmstedges+1
            if (nmstedges.ge.n_snaps) then
              write(ilog,*) 'Fatal. The number of edges in the approximate minimum spanning tree is not correct. Please report &
     &this bug.'
              call fexit()
            end if
            mstedges(1,nmstedges) = tmptree(e)%mine(1) ! edgelst(e,1)
            mstedges(2,nmstedges) = tmptree(e)%mine(2) ! edgelst(e,2)
            lmstedges(nmstedges) = tmptree(e)%mind     ! ledgelst(e)
    ! update pointer for tree and all its siblings obtained through prior merge operations
            tmptree(i1)%ptr = i2
            tmptree(i2)%nsiblings = tmptree(i2)%nsiblings + 1
            if (tmptree(i2)%nsiblings.gt.tmptree(i2)%nsibalsz) call pidxtree_growsiblings(tmptree(i2))
            tmptree(i2)%siblings(tmptree(i2)%nsiblings) = i1
            do j=1,tmptree(i1)%nsiblings
              mm = tmptree(i1)%siblings(j)
              tmptree(i2)%nsiblings = tmptree(i2)%nsiblings + 1
              if (tmptree(i2)%nsiblings.gt.tmptree(i2)%nsibalsz) call pidxtree_growsiblings(tmptree(i2))
              tmptree(i2)%siblings(tmptree(i2)%nsiblings) = mm
              tmptree(mm)%ptr = i2
            end do
    ! decrease forest size
          end if
        end do
    ! reassign csnap2tree based on pointer
        do j=1,n_snaps
          csnap2tree(j) = tmptree(csnap2tree(j))%ptr
        end do
    ! now construct a new set of sorted from csnap2tree and readjust csnap2tree (this does not scale)
        mm = 0
        oldntrees = ntrees
        ntrees = 0
        do j=1,n_snaps
          k = csnap2tree(j)
          if (testlst(k,3).ge.0) then
            ntrees = ntrees + 1
            testlst(k,3) = -ntrees
          end if
        end do
    !
        do j=1,n_snaps
          csnap2tree(j) = -testlst(csnap2tree(j),3)
        end do
        do j=ntrees+1,oldntrees
          if (allocated(tmptree(j)%siblings).EQV..true.) deallocate(tmptree(j)%siblings)
        end do
        tmptree(1:ntrees)%mind = HUGE(jkdist)
    !
        boruvkasteps = boruvkasteps+1
    ! manage stored guesses
        do j=1,n_snaps
          ishfx = j
          satisfied(ishfx) = 0
          ixx2 = 0
          do ixx=1,N_NEARS
            if (tmpix(ixx,ishfx).gt.0) then
              if (csnap2tree(j).ne.csnap2tree(tmpix(ixx,ishfx))) then
                ixx2 = ixx2 + 1
                tmpix(ixx2,ishfx) = tmpix(ixx,ishfx)
                tmpdis(ixx2,ishfx) = tmpdis(ixx,ishfx)
              end if
            end if
          end do
          do ixx=max(1,ixx2),N_NEARS
            tmpix(ixx,ishfx) = 0
            tmpdis(ixx,ishfx) = HUGE(tmpdis(ixx,ishfx))
          end do
        end do
    !
     567 format('... time for Boruvka stage ',i4,': ',g11.3,'s ...')
        ! write(ilog,*) "ntrees:", ntrees
    !  456 format(i4,8(g12.4,1x),i10,i10,i6,i6,1x,g12.4,g12.4,g12.4)
        call CPU_time(tmvars(2))
        write(ilog,567) boruvkasteps, tmvars(2)-tmvars(1)
        tmvars(1) = tmvars(2)
        if(csnap2tree_ik_eq_cnt.gt.0) then
          write(ilog,*) "Warning: we found exit flags for csnap2tree identity. Amount:", csnap2tree_ik_eq_cnt
        end if
    !
      end do
    !
      cdevalcnt = testcnt
    !
      deallocate(testlst)
      deallocate(snplst)
      deallocate(satisfied)
      deallocate(tmpdis)
      deallocate(tmpix)
    !
      if (n_snaps-1.ne.nmstedges) then
        write(ilog,*) 'Fatal. The number of edges in the approximate minimum spanning tree is not correct. Please report &
     &this bug.'
        call fexit()
      end if
     587 format(' Weight of Short Spanning Tree: ',1x,g12.5,a)
      write(ilog,*)
      if (tmp_dis_method.le.2) then
        write(ilog,587) sum(lmstedges(1:n_snaps-1)),' degrees'
      else if (tmp_dis_method.le.4) then
        write(ilog,587) sum(lmstedges(1:n_snaps-1)),' '
      else if (tmp_dis_method.le.10) then
        write(ilog,587) sum(lmstedges(1:n_snaps-1)),' Angstrom'
      else if (tmp_dis_method.eq.11) then
        write(ilog,587) sum(lmstedges(1:n_snaps-1)),' degrees'
      end if

      write(ilog,*)
      allocate(approxmst(n_snaps))
      approxmst(1:n_snaps)%deg = 0
      write(ilog,*) '...ordering mst into an adjlist...'
      call CPU_time(t3)
      call gen_MST_w(mstedges,lmstedges,approxmst)
      call CPU_time(t4)
      write(ilog,*) '...done using ',t4-t3, ' [s]'
      do i=1,ntrees
        if (allocated(tmptree(i)%snaps).EQV..true.) deallocate(tmptree(i)%snaps)
        if (allocated(tmptree(i)%siblings).EQV..true.) deallocate(tmptree(i)%siblings)
      end do
      deallocate(tmptree)
      deallocate(csnap2tree)
      deallocate(csnap2clus)
      deallocate(mstedges)
      deallocate(lmstedges)
      if (allocated(birchtree).EQV..true.) then
        do i=1,c_nhier+1
          do j=1,size(birchtree(i)%cls) ! birchtree(i)%ncls
            if (allocated(birchtree(i)%cls(j)%snaps).EQV..true.) deallocate(birchtree(i)%cls(j)%snaps)
            if (allocated(birchtree(i)%cls(j)%tmpsnaps).EQV..true.) deallocate(birchtree(i)%cls(j)%tmpsnaps)
            if (allocated(birchtree(i)%cls(j)%sums).EQV..true.) deallocate(birchtree(i)%cls(j)%sums)
            if (allocated(birchtree(i)%cls(j)%children).EQV..true.) deallocate(birchtree(i)%cls(j)%children)
          end do
          deallocate(birchtree(i)%cls)
        end do
        deallocate(birchtree)
      end if
    !
    !  77 format(a,20(i18,1x))
      write(ilog,*) '... and after ',cdevalcnt,' additional distance evaluations.'
      write(ilog,*)
    !
    end
    !
    !-----------------------------------------------------------------------------------
    !
    subroutine pidxtree_growsiblings(t1)
    !
      type(t_progindextree) t1
      integer tmp(max(t1%nsibalsz,1)),k
    !
      if (t1%nsibalsz.eq.0) then
        t1%nsibalsz = 5
        allocate(t1%siblings(t1%nsibalsz))
      else
        tmp(:) = t1%siblings(1:t1%nsibalsz)
        deallocate(t1%siblings)
        k = t1%nsibalsz
        t1%nsibalsz = t1%nsibalsz*2
        allocate(t1%siblings(t1%nsibalsz))
        t1%siblings(1:k) = tmp(:)
      end if
    !
    end
    !
    !-----------------------------------------------------------------------------------------

end module m_mst_dumping
