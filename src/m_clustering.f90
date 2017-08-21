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
!                        Nicholas Lyle, Nicolas Bloechliger,               !
!                        Davide Garolini                                   !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN:      Andreas Vitalis                                               !
! WRAPPER:   Davide Garolini                                               !
!                                                                          !
!--------------------------------------------------------------------------!
! trj_data            : data extracted for clustering (all read into memory -> large)
! n_snaps            : number of snapshots stored for clustering (-> trj_data)
! cmode              : choice of algorithm for clustering or similar tasks (1-5)
! dis_method          : the distance criterion underlying the clustering (1-8)
! cstorecalc         : frequency for collecting data for clustering
! cmaxsnaps          : allocation size for trj_data in numbers of snapshots dimension
! calcsz,clstsz      : dimensionality size variables for stored data (-> trj_data)
! cfilen             : input file for defining clustering subset
! cdofset/cl_imvec   : auxiliary variables for clustering data
! cradius            : general size threshold for clustering
! cprepmode          : options for signal (data) preprocessing
! cchangeweights     : for eligible dis_method, replace weights for distance evaluations
! cl_mwvec           : keep track of mean weight per dimension
! cwcombination      : for locally adaptive weights, how to combine (L_? norm)
! cwwindowsz         : for locally adaptive weights, how big a window to use
! cdynwbuf           : for LAWs dependent on transition counts, buffer parameter for taking inverse
! cwacftau           : for ACF-based weights, fixed lag time
! csmoothie          : order for smoothing fxn (sets window size)
! refine_clustering  : general flag to trigger refinement in clustering
! birchtree,scluster : arrays to hold clustering result
! sumssz             : global size variable for (t_scluster)%sums(:,:)
! cleadermode        : if cmode = 1/2/5(4): processing direction flags
! clinkage           : if cmode = 3: linkage criterion (1-3)
! chardcut           : if cmode = 3(4): truncation cut for snapshot NB-list
! nblfilen,read_nbl_from_nc,cnblst,cnc_ids : if cmode = 3(4): snapshot NB-list variables
! caddlkmd           : whether and how to augment network (add missing links)
! caddlinkwt         : what type of FP weight to use for added links
! cprogindex         : if cmode = 4: exact or approximate scheme (1-2)
! cprogindrmax       : if cmode = 4: number of nearest neighbor guesses for approximate scheme
! cprogindstart      : if cmode = 4: snapshot to start from
! cprogpwidth        : if cmode = 4: the (maximum) size of A/B partitions in localized cut
! cprogfold          : if cmode = 4: how often to collapse terminal vertices of the MST
! cprogrdepth        : if cmode = 4 and cprogindex = 2: auxiliary search depth
! cprogbatchsz       : if cmode = 4 and cprogindex = 2: batch size for random stretches
! csivmax,csivmin    : if cmode = 4: used to identify starting points automatically
! c_nhier            : if cmode = 5(4): tree height
! cmaxrad            : if cmode = 5(4): root level size threshold
! c_multires         : if cmode = 5(4): how many levels to recluster in multi-pass and print
! ccfepmode          : if any clustering is performed, type of cFEP to produce
! align_for_clustering: if dis_method = 5/6 whether to do alignment
! cdofsetbnds        : if dis_method = 6, selectors for separating sets for alignment and distance compu.
! pcamode            : whether to do PCA and what output to produce (1-3)
! reduced_dim_clustering: whether to proceed with clustering in reduced dimensional space after PCA
! align              : structure used for trajectory alignment (ALIGNFILE)
! csnap2tree         : temporary array for cmode = 4 and cprogindex = 2
! csnap2clus         : temporary array for cmode = 4 and cprogindex = 2
! tmptree            : temporary array for cmode = 4 and cprogindex = 2
! ntbrks,ntbrks2,itbrklst : management of user-requested trajectory breaks (counts and counter)
! ntlnks             : management of user-requested link additions (count)
!
module m_clustering
  
  ! distance matrix must be allocatable
  real, allocatable :: distance_mat(:,:)

  ! cluster object
  type t_scluster
    integer nmbrs, alsz !alsz is the allocation size of the snaps principally
    real sqsum
    real, ALLOCATABLE :: sums(:)
    integer, ALLOCATABLE :: snaps(:) ! snapshot lists
    integer, ALLOCATABLE :: children(:) !for tree-based clustering
    ! for tree-based clustering
    integer nchildren
    integer childr_alsz
    integer parent
    !for quality
    real diam
    real radius
    real quality
    integer center
    !for SST
    integer, ALLOCATABLE :: tmpsnaps(:)
    ! nbalsz - resizenb - expand graph nb lists - we assume unallocated member arrays to be dealt with
    ! elsewhere
    ! like scluster_resize must never be called with allocation size currently zero
  end type t_scluster

  type t_ctree
    type(t_scluster), ALLOCATABLE:: cls(:) ! nothing but an array of clusters
    integer ncls, n_clu_alsz_tree ! real and allocation sizes
  end type t_ctree


  type(t_scluster), ALLOCATABLE :: scluster(:)
  type(t_ctree), ALLOCATABLE :: birchtree(:)
  real radius !radius for clusters belonging limit
  integer n_clu_alc_sz_gen ! number of clusters allocatable in total
  integer nclu ! number of clusters allocated
  integer c_nhier !height of the tree
  integer c_multires !FMCSC_BIRCHMULTI initial default = 0
  integer ordering
  real cmaxrad
  logical clu_summary !helper var for precise_clu_descr
  integer precise_clu_descr !the clusters are shortened and presented properly
  ! This are the sizes for sums variable that keeps the centroid measures.
  ! In the rmsd (dis_method = 5) only needs sumssz = 1 and clstsz = calcsz
  ! (calcsz is the number of atoms * coordinates)
  contains
    !--------------------------------------------------------------------
    ! LEADER CLUSTERING
    !--------------------------------------------------------------------
    ! each snapshot is added to a cluster on the base of a defined threshold.
    !The original structure is based on:
    !leader_clustering, gen_MST_from_nbl,do_prog_ind
    subroutine leader_clustering(trj)
      use m_variables_gen
      use gutenberg
      implicit none

      real, intent(in) :: trj(n_snaps,n_xyz)
      real tmp_d ! temporary variable for radiuses and distances
      integer i, j, ii, u
      real vecti(n_xyz)

      n_clu_alc_sz_gen = 2 ! a priori initial cluster numbers
      allocate(scluster(n_clu_alc_sz_gen))
      scluster(:)%nmbrs = 0 ! number of elements in each cluster
      scluster(:)%alsz = 0 ! allocation size of each cluster


      call sl()
      call spr('-----------------------------------')
      call spr('Now using truncated leader algorithm for pre-clustering...')
      call sl()

      nclu = 0 !cluster index
      do i=1,n_snaps
        ii = -1
        do j=nclu,max(nclu-500,1),-1 ! when nclu = 0 it skips this do  [j=0,1,-1]
          ! if(superver) write(ilog,*) "now checking dist(cluster,snap) : (",j,",",i,")"

          ! j is going down of -1, not nclu. nclu is growing if (ii.eq.-1)(new cluster is
          ! added). The number of clusters that can be considered in order to put i
          ! (snapshot) in one of those is 500. Once exceeded this upper bound
          ! it starts to shift on that limit forwards.
          call preprocess_snapshot(trj,i, vecti)
          call snap_to_cluster_d(tmp_d,scluster(j),vecti)
          if (tmp_d.lt.radius) then
            call cluster_addsnap(scluster(j),vecti,i)
            ! add snapshot i in the cluster j
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          ! a new cluster is added to the cluster list because tmp_d .gt. radius
          nclu = nclu + 1
          ! enlarge the size of the clusters array by the double
          if (nclu.gt.n_clu_alc_sz_gen) then
            call cluster_lst_resize(scluster)
          end if

          call cluster_addsnap(scluster(nclu),vecti,i) ! nclu new cluster
          ! if(superver) write(ilog,*) 'Cluster number ',nclu, ' created for snap ', i
        end if
      end do

      call sipr('... done with total number of clusters equal: ', nclu)
    end subroutine

! helper function to resize the cluster
    subroutine cluster_lst_resize(clu)
      use m_variables_gen
      implicit none
      type(t_scluster), ALLOCATABLE :: tmpcls(:)
      type(t_scluster), ALLOCATABLE :: clu(:)
      integer oldsz
      allocate(tmpcls(n_clu_alc_sz_gen))
      oldsz = n_clu_alc_sz_gen
      tmpcls = clu
      deallocate(clu)
      n_clu_alc_sz_gen = n_clu_alc_sz_gen * 2
      allocate(clu(n_clu_alc_sz_gen))
      clu(1:oldsz) = tmpcls(1:oldsz)
      clu(oldsz+1:n_clu_alc_sz_gen)%nmbrs = 0
      clu(oldsz+1:n_clu_alc_sz_gen)%alsz = 0
      deallocate(tmpcls)
    end subroutine cluster_lst_resize
    !--------------------------------------------------------------------
    ! BIRCH CLUSTERING
    !--------------------------------------------------------------------
    !it is based on a kind of unbalanced BIRCH algorithm
    !The original structure is based on:
    !birch_clustering,gen_MST_from_treeclustering, do_prog_ind
    !original INPUT
    !cmode = 4 (prog index) 5 (only clustering) (this var is not used here)
    !nstruccls = nnodes  This variable is not used outside this function.
    ! it will be reasonably forced in it
    !OLD vars mapping
    !c_nhier = height of the tree
    !cleadermode = if cmode = 1/2/5(4): processing direction flags
    !scrcts = treshold array.
    !scrcts(1) = root level. The thresh is not imporant. All the data-set
    !scrcts(2) = cmaxrads
    !scrcts(c_nhier+1) = cradius
    !scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - cradius)
    subroutine birch_clustering(trj)
    !
      use gutenberg
      use m_variables_gen
    !
      implicit none
    !
      integer i,j,k,jj,kk,ll,mm
      integer ii !level of the tree
      integer fail,kkf,thekk,nlst1,nlst2
      integer snapstart,snapend,snapinc !sequence of clustering(in,en,step)
      real, intent(in) :: trj(n_snaps,n_xyz)
      real vecti(n_xyz)

      ! integer nnodes ! number of nodes ? it is assigned here and used here only
      ! integer modei ! not used here
      integer, ALLOCATABLE :: kklst(:,:),kkhistory(:)
      integer cnt1,cnt2 !counter variables
      real, ALLOCATABLE:: scrcts(:) !threshold array
      real tmp_d !temporary variable for the distance
      real min_d,maxd(2),qualmet(4)
      ! real helper, normer(4) !no more used
      logical notdone, cycit
      integer cl_one_elem !count of clusters with one single element
      integer cnhmax ! Toask: levels of the tree that are recalculated?
      integer max_show_clu ! maximum number of cl to show.

    ! NOT USED it seems
    !  33 format('ERROR B: ',i6,' could be part of ',i5,' at ',g14.7,' (last d',g14.6,')')
    !  34 format('NEXT HIGHER: ',i6,' could have been part of ',i5,' at ',g14.7,'.')
    !  35 format('ERROR C: ',i6,' is meant to a child of ',i5,' at level ',i5,' but has distance ',g14.7,'.')
    !
      if (ordering.le.2) then ! cleadermode = processing direction flags (default = 1)
        snapstart = 1
        snapend = n_snaps
        snapinc = 1
      else
        snapstart = n_snaps
        snapend = 1
        snapinc = -1
      end if
      ! Variables initialization
      allocate(kklst(n_snaps,2)) !This is two layers of snapshots (integer)
      allocate(kkhistory(c_nhier+1))
      allocate(birchtree(c_nhier+1))
      do ii=1,c_nhier+1
        allocate(birchtree(ii)%cls(10))
        birchtree(ii)%ncls = 0
        birchtree(ii)%n_clu_alsz_tree = 10
      end do
      allocate(scrcts(c_nhier+1))
      birchtree(1)%ncls = 1 !number of clusters
      scrcts(c_nhier+1) = radius !leaves level
      if (c_nhier.gt.1) then
        scrcts(2) = cmaxrad
        do i=3,c_nhier
          scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - radius) ! linear so far
        end do
      end if

      do i=1,c_nhier+1
        birchtree(i)%cls(:)%nmbrs = 0
        ! birchtree(i)%cls(:)%nb = 0
        birchtree(i)%cls(:)%nchildren = 0
        birchtree(i)%cls(:)%alsz = 0
        birchtree(i)%cls(:)%childr_alsz = 0
        birchtree(i)%cls(:)%parent = 0
      end do

      allocate(birchtree(1)%cls(1)%snaps(2))
      birchtree(1)%cls(1)%alsz = 2
      allocate(birchtree(1)%cls(1)%children(2))
      birchtree(1)%cls(1)%childr_alsz = 2
    !
      call sl()
      call spr('---------------------------------------------------------------------')
      call spr('Now performing tree-based clustering ...')
      call sl()
      cycit = .false.
      cnt1 = 0 !number of calculated distances
      cnt2 = 0
      do i=snapstart,snapend,snapinc
        kk = 1
        notdone = .true.
        fail = -1
        nlst1 = 1
        kklst(1,1) = kk
    !    cnt1 = 0
        do ii=2,c_nhier !for each level (from 2)
          jj = -1
          min_d = HUGE(min_d) !the biggest number available for that (real) type
          nlst2 = 0
          do mm=1,nlst1
            kk = kklst(mm,1)
            do j=1,birchtree(ii-1)%cls(kk)%nchildren !for each children (nothing done if no children)
              ll = birchtree(ii-1)%cls(kk)%children(j)
              cnt1 = cnt1 + 1 !number of calculated distances
              if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.n_snaps)) then
                call spr("BUG: the center is not allocated. Please contact the admin.")
                call fexit()
              end if
              !distance between the one element i with all the children clusters ll
              call preprocess_snapshot(trj, i, vecti)
              call snap_to_cluster_d(tmp_d, birchtree(ii)%cls(ll), vecti)
              if (tmp_d.lt.min_d) then !I care only about the minimum distance for storing
                min_d = tmp_d !minimum distance is in the j-th child
                jj = j !the children cluster (ll) number
                thekk = kk
              end if
              if (tmp_d.lt.scrcts(ii)) then ! snap i belong to the cluster
                nlst2 = nlst2 + 1 !number of adds
                kklst(nlst2,2) = ll !cluster that has it added
              end if
            end do
          end do
    !     store the path
          if (jj.eq.-1) then
            kkhistory(ii-1) = 1 !e.g. no children
          else
            kkhistory(ii-1) = thekk !the minimum!
          end if
          !STANDARD ADD - there is a nearest cluster and it is not a leaf
          if ((ii.le.(c_nhier-1)).AND.(jj.gt.0)) then !jj gt 0 it means only if it i is a minimum of some cluster ll
            if (nlst2.le.0) then ! absolutely nothing nearby (using the scrcts threshold)
              nlst1 = 1
              kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj) !it stores the nearest cluster that is NOT near enough
              if (fail.eq.-1) then
    !            write(ilog,*) i,' failed at ',ii,' w/ ',min_d
                fail = ii
                kkf = thekk
              end if
            else !  if NOT (nlst2.le.0) (it is not alone)
              if (fail.gt.0) then
                do j=fail,ii-1 !from the point it found a value alone far away from the others until the last level
                  birchtree(j)%ncls = birchtree(j)%ncls + 1 !generate a new cluster
                  if (birchtree(j)%ncls.gt.birchtree(j)%n_clu_alsz_tree) &
                    call scluster_resizelst(birchtree(j)%n_clu_alsz_tree,birchtree(j)%cls)
                  kkhistory(j) = birchtree(j)%ncls
                  if (j.gt.fail) then
                    call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
     &                                    birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                    cnt2 = cnt2 + 1
                  end if
                end do
                call cluster_addchild(birchtree(fail-1)%cls(kkf),kkf,&
                                      birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
                call cluster_addchild(birchtree(ii-1)%cls(birchtree(ii-1)%ncls),&
     &                                birchtree(ii)%cls(jj)%parent,birchtree(ii)%cls(jj),jj)
                fail = -1
              end if
              nlst1 = 1
              kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
            end if
            cycit = .true.
          ! LEAF ADD - to access this you must have found min_d (a minimum distance INSIDE a cluster (<scrcts(lev)))
          else if ((ii.eq.c_nhier).AND.(min_d.lt.scrcts(ii))) then
            kk = birchtree(ii-1)%cls(thekk)%children(jj)
            call preprocess_snapshot(trj, i, vecti)
            call cluster_addsnap(birchtree(ii)%cls(kk),vecti,i)
            do j=2,c_nhier
              call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),vecti,i)
            end do
            notdone = .false.
          ! NEW CLUSTER
          else  !case in which you must add one cluster
            if (fail.eq.-1) then !e.g. first add
              fail = ii
              kkf = kk !kk = kklst(mm,1) mm=1,nlst1 (kk=1 first add)
            end if
            kk = kkf
            do j=fail,c_nhier !create a new cluster and addsnap all the way to the leaf
              birchtree(j)%ncls = birchtree(j)%ncls + 1 !for each level add 1 cluster
              if (birchtree(j)%ncls.gt.birchtree(j)%n_clu_alsz_tree) &
                call scluster_resizelst(birchtree(j)%n_clu_alsz_tree, birchtree(j)%cls)
              call preprocess_snapshot(trj, i, vecti) !working
              call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),vecti,i) !working
              if (j.gt.fail) then
                call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
     &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                cnt2 = cnt2 + 1
              end if
            end do
            call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)!last level up
            cnt2 = cnt2 + 1
            do j=2,fail
              call preprocess_snapshot(trj, i, vecti)
              call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),vecti,i)
            end do
            call preprocess_snapshot(trj, i, vecti)
            call snap_to_cluster_d(maxd(1),birchtree(fail-1)%cls(kk),vecti)
            notdone = .false.
          end if
          if (cycit.EQV..true.) then
            cycit = .false.
            cycle
          end if
          if (notdone.EQV..false.) exit
        end do
      end do
      ! do i=1,c_nhier+1
      !   do j=1,birchtree(i)%ncls
      !     write(ilog,*) "lev/ncls: ", i,j
      !     write(ilog,*) "numbers: ", birchtree(i)%cls(j)%nmbrs
      !     write(ilog,*) "center: ", birchtree(i)%cls(j)%center
      !   end do
      ! end do
      call spr("multi-passing...")
      call sl()
    !  multi-pass
      do cnhmax=c_nhier+1,min(c_nhier+1,max(3,c_nhier+1-c_multires)),-1
        if (cnhmax.lt.(c_nhier+1)) then !not done if == (always the case with c_multires = 0)
          ! write(ilog,*) "and the center?"
          birchtree(cnhmax)%ncls = 0
          birchtree(cnhmax)%cls(:)%nmbrs = 0
          birchtree(cnhmax)%cls(:)%parent = 0
          birchtree(cnhmax-1)%cls(:)%nchildren = 0
        end if
        ! write(ilog,*) "cnhmax: ",cnhmax
    !
        cycit = .false.
        do i=snapstart,snapend,snapinc
          kk = 1
          notdone = .true.
          fail = -1
          nlst1 = 1
          kklst(1,1) = kk
          do ii=2,cnhmax
            jj = -1
            min_d = HUGE(min_d)
            nlst2 = 0
            do mm=1,nlst1
              kk = kklst(mm,1)
              do j=1,birchtree(ii-1)%cls(kk)%nchildren
                ll = birchtree(ii-1)%cls(kk)%children(j)
                cnt1 = cnt1 + 1
                if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.n_snaps))then
                  ! if(i.lt.10) then
                  !    write(ilog,*) "center: ",birchtree(ii)%cls(ll)%center
                  !    write(ilog,*) "ii-lev/ll-childrenID(upperlev):",ii,ll
                  !    write(ilog,*) "number of elem:", birchtree(ii)%cls(ll)%nmbrs
                  !  end if
                  call spr("BUG: the center is not allocated. Please contact the admin.")
                  call fexit()
                end if
                call preprocess_snapshot(trj, i, vecti) !working
                call snap_to_cluster_d(tmp_d,birchtree(ii)%cls(ll),vecti)
                if (tmp_d.lt.min_d) then
                  min_d = tmp_d
                  jj = j
                  thekk = kk
                end if
                if (tmp_d.lt.scrcts(ii)) then
                  nlst2 = nlst2 + 1
                  kklst(nlst2,2) = ll
                end if
              end do
            end do
    !       store the path
            if (jj.eq.-1) then
              kkhistory(ii-1) = 1
            else
              kkhistory(ii-1) = thekk
            end if
            if ((ii.lt.cnhmax).AND.(jj.gt.0)) then
              if (nlst2.le.0) then ! absolutely nothing nearby
                nlst1 = 1
                kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
                if (fail.eq.-1) then
    !              write(*,*) i,' failed at ',ii,' w/ ',min_d
                  fail = ii
                  kkf = thekk
                end if
              else
                nlst1 = 1
                kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
              end if
              cycit = .true.
    !       leaf
            else if ((ii.eq.cnhmax).AND.(min_d.lt.scrcts(ii))) then
              kk = birchtree(ii-1)%cls(thekk)%children(jj)
              call preprocess_snapshot(trj, i, vecti)
              call cluster_addsnap(birchtree(ii)%cls(kk),vecti,i)
              notdone = .false.
            else if (fail.gt.0) then
              kk = kkf
              do j=fail,cnhmax
                birchtree(j)%ncls = birchtree(j)%ncls + 1
                if (birchtree(j)%ncls.gt.birchtree(j)%n_clu_alsz_tree) &
                  call scluster_resizelst(birchtree(j)%n_clu_alsz_tree,birchtree(j)%cls)
                call preprocess_snapshot(trj, i, vecti)
                call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),vecti,i)
                if (j.gt.fail) then
                  call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
       &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
                           cnt2 = cnt2 + 1
                end if
              end do
              call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
              cnt2 = cnt2 + 1
              notdone = .false.
            else
              fail = ii
              j = cnhmax
              birchtree(j)%ncls = birchtree(j)%ncls + 1
              if (birchtree(j)%ncls.gt.birchtree(j)%n_clu_alsz_tree) &
                call scluster_resizelst(birchtree(j)%n_clu_alsz_tree,birchtree(j)%cls)
              call preprocess_snapshot(trj, i, vecti) !working
              call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),vecti,i)
              call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
              notdone = .false.
            end if
            if (cycit.EQV..true.) then
              cycit = .false.
              cycle
            end if
            if (notdone.EQV..false.) exit
          end do
        end do
      end do

      ! do i=1,c_nhier+1
      !   do j=1,birchtree(i)%ncls
      !     write(ilog,*) "lev/ncls: ", i,j
      !     write(ilog,*) "numbers: ", birchtree(i)%cls(j)%nmbrs
      !     write(ilog,*) "center: ", birchtree(i)%cls(j)%center
      !   end do
      ! end do

      !PRINTING RESULTS
      call spr("--------------------------- TREE SUMMARY ----------------------------")
      ! call spr('Thresholds: ')
      ! call rvpr(scrcts)
      call spr('Level    #Clusters     Threshold      TotalSnaps    Tot Children')
      ! call spr('Level    #Clusters      TotalSnaps    Tot Children')
    !  67 format(i5,3x,i10,4x,g14.4,1x,i11,4x,i12)
    !  68 format(i5,3x,i10,7x,a7,5x,i11,4x,i12)
    !   write(ilog,66)
    !   write(ilog,68) 1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
    !   sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
      do i=1,c_nhier+1  ! i was 2 with MAXIMAL level
        if(i.eq.1) then
          call clu_summary_maximal(i,birchtree(i)%ncls,"MAXIMAL",&
          sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
          sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren))
        else
          call clu_summary_linepr(i,birchtree(i)%ncls,scrcts(i),&
          sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
          sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren))
        end if
      end do
      call spr('---------------------------------------------------------------------')
      call sl()
      call sipr('... done after total distance evaluations of: ',cnt1)
      call sl()

    !
      ! do j=1,birchtree(c_nhier+1)%ncls
      !   call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
      ! end do
    !
    !   if (refine_clustering.EQV..true.) then !this is the refinement !TODO
    !     call quality_of_clustering(birchtree(c_nhier+1)%ncls,birchtree(c_nhier+1)%cls(1:c_nhier+1),scrcts(c_nhier+1),qualmet)
    !     write(ilog,*) 'Now merging clusters that yield joint reduced average intracluster distance ...'
    !  77 format('Would join ',i5,' (',i6,') and ',i5,'(',i6,') from: ',/,'Diam: ',g10.4,' / ',g10.4,'; Rad.: ',&
    !  &g10.4,' / ',g10.4,' to ',g10.4,' / ',g10.4,'.')
    !     cnt1 = 0
    !     cnt2 = 0
    !     do i=1,birchtree(c_nhier)%ncls
    !       nlst1 = 0
    !       do kkf=1,birchtree(c_nhier)%ncls
    !         if (i.eq.kkf) cycle
    !         call cluster_to_cluster_d(tmp_d,birchtree(c_nhier)%cls(i),birchtree(c_nhier)%cls(kkf))
    !         cnt1 = cnt1 + 1
    !         if (tmp_d.lt.scrcts(c_nhier)) then
    !           nlst1 = nlst1 + 1
    !           kklst(nlst1,1) = kkf
    !         end if
    !       end do
    !       do j=1,birchtree(c_nhier)%cls(i)%nchildren
    !         jj = birchtree(c_nhier)%cls(i)%children(j)
    !         if (birchtree(c_nhier+1)%cls(jj)%nmbrs.le.0) cycle
    !         do kkf=1,nlst1
    !           do k=1,birchtree(c_nhier)%cls(kklst(kkf,1))%nchildren
    !             kk = birchtree(c_nhier)%cls(kklst(kkf,1))%children(k)
    !             if (birchtree(c_nhier+1)%cls(kk)%nmbrs.le.0) cycle
    !             if (jj.eq.kk) cycle
    !             call clusters_joint_diam(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk),tmp_d,helper)
    !             cnt1 = cnt1 + 1
    !             normer(1) = 0.5*birchtree(c_nhier+1)%cls(jj)%nmbrs*(birchtree(c_nhier+1)%cls(jj)%nmbrs-1.0)
    !             normer(2) = 0.5*birchtree(c_nhier+1)%cls(kk)%nmbrs*(birchtree(c_nhier+1)%cls(kk)%nmbrs-1.0)
    !             normer(3) = 1.0/(1.0*(birchtree(c_nhier+1)%cls(jj)%nmbrs + birchtree(c_nhier+1)%cls(kk)%nmbrs))
    !             normer(4) = 0.0
    !             if (sum(normer(1:2)).gt.0) normer(4) = 1.0/sum(normer(1:2))
    !             maxd(1) =  normer(4)*(normer(1)*birchtree(c_nhier+1)%cls(jj)%diam + &
    !               &                   normer(2)*birchtree(c_nhier+1)%cls(kk)%diam)
    !             maxd(2) =  normer(3)*(birchtree(c_nhier+1)%cls(jj)%nmbrs*birchtree(c_nhier+1)%cls(jj)%radius + &
    !    &                              birchtree(c_nhier+1)%cls(kk)%nmbrs*birchtree(c_nhier+1)%cls(kk)%radius)
    !             if ((helper.le.maxd(2)).OR.(tmp_d.le.maxd(1))) then
    !               if (birchtree(c_nhier+1)%cls(jj)%nmbrs.gt.birchtree(c_nhier+1)%cls(kk)%nmbrs) then
    !                 call join_clusters(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk))
    !                 birchtree(c_nhier+1)%cls(jj)%diam = tmp_d
    !                 birchtree(c_nhier+1)%cls(jj)%radius = helper
    !                 cnt2 = cnt2 + 1
    !               else
    !                 call join_clusters(birchtree(c_nhier+1)%cls(kk),birchtree(c_nhier+1)%cls(jj))
    !                 birchtree(c_nhier+1)%cls(kk)%diam = tmp_d
    !                 birchtree(c_nhier+1)%cls(kk)%radius = helper
    !                 cnt2 = cnt2 + 1
    !                 exit
    !               end if
    !             end if
    !           end do
    !           if (birchtree(c_nhier+1)%cls(jj)%nmbrs.eq.0) exit
    !         end do
    !       end do
    !     end do
    !     write(ilog,*) '... done after a total of ',cnt2,' merges requiring ',cnt1,' additional &
    !  &distance or joint size evaluations.'
    !     write(ilog,*)
    !   end if
    !
      ! do j=1,birchtree(c_nhier+1)%ncls
      !   call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
      ! end do

    ! now shorten list and resorti


    ! ! lastly, copy into global cluster array !NOT USED IF NOT IN PRINTING THE GRAPH
    !   allocate(scluster(birchtree(c_nhier+1)%ncls))
    !   !initializing the vars
    !   scluster(:)%diam = 0
    !   scluster(:)%radius = 0
    !   scluster(:)%quality = 0
    !   do i=1,birchtree(c_nhier+1)%ncls
    !     call copy_cluster(birchtree(c_nhier+1)%cls(i),scluster(i))
    !   end do
      ! nnodes = birchtree(c_nhier+1)%ncls
      ! do k=2,c_nhier+1

      do k=1,c_nhier+1
        call clusters_shorten(birchtree(k)%cls, birchtree(k)%ncls)
        call clusters_sort(birchtree(k)%cls, birchtree(k)%ncls)
        do i=1,birchtree(k)%ncls
          call cluster_getcenter(birchtree(k)%cls(i),trj)
        end do
        do i=1,birchtree(k)%ncls
          call cluster_calc_params(birchtree(k)%cls(i),scrcts(k))
        end do
      end do

      ! call cluster_getcenter(birchtree(1)%cls(1),trj)
      ! call cluster_calc_params(birchtree(1)%cls(1),huge(scrcts(k)))
      if(clu_summary) then
        call spr('------------------------- CLUSTERS SUMMARY -----------------------------')
        call spr('level    threshold         #     No.    Center      Diameter      Radius      ')
        do k=1,c_nhier+1
          cl_one_elem = 0
          j = 0
    !    do i=1,birchtree(c_nhier+1)%ncls
    !      write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
    !    end do
          ! do i=1,birchtree(k)%ncls
          if(precise_clu_descr.eq.0.or.precise_clu_descr.ge.birchtree(k)%ncls) then
            max_show_clu = birchtree(k)%ncls
          else
            max_show_clu = precise_clu_descr
          end if

          do i=1,birchtree(k)%ncls
            if(birchtree(k)%cls(i)%nmbrs.ne.1) then
              if(j.le.max_show_clu) then
                ! 63 format(i6,3x,g12.5,1x,i7,1x,i7,2x,i8,3x,2(1x,g11.5))
                if(k.eq.1) then
                  ! write(ilog,63) k,"",i,birchtree(k)%cls(i)%nmbrs,birchtree(k)%cls(i)%center,&
                  ! birchtree(k)%cls(i)%diam,birchtree(k)%cls(i)%radius
                  call clu_summary_details_pr(k,0.0,i,birchtree(k)%cls(i)%nmbrs,&
                  birchtree(k)%cls(i)%center,&
                  birchtree(k)%cls(i)%diam, birchtree(k)%cls(i)%radius)
                else
                  call clu_summary_details_pr(k,scrcts(k),i,birchtree(k)%cls(i)%nmbrs,&
                  birchtree(k)%cls(i)%center,&
                  birchtree(k)%cls(i)%diam,birchtree(k)%cls(i)%radius)
                end if
                j = j + 1
              end if
            else
              cl_one_elem = cl_one_elem + 1
            end if
          end do
          if(cl_one_elem.gt.0) then
            ! 70 format(' -> level ',i3,' has ',i7,' clusters with one element')
            call sivpr(' -> level, number of 1-elem clusters: ',(/k, cl_one_elem/))
            ! write(ilog,70) k, cl_one_elem
          end if
        end do

        call sl()
        call sl()
        do k=2,c_nhier+1
          if(birchtree(k)%ncls.ne.n_snaps) then
            call quality_of_clustering(birchtree(k)%ncls,&
              birchtree(k)%cls(1:birchtree(k)%ncls),scrcts(k),qualmet)
          end if
        end do
      end if
      ! call gen_graph_from_clusters(scluster,nnodes,atrue)
      ! call graphml_helper_for_clustering(scluster,nnodes)
      ! call vmd_helper_for_clustering(scluster,nnodes)
    !
      deallocate(kklst)
      deallocate(kkhistory)
      deallocate(scrcts)
    !
    end
    !
    !-------------------------------------------------------------------------------
    !HELPER FUNCTION FOR BIRCH
    subroutine scluster_resizelst(currentalsz,it)
    !
      implicit none
    !
      integer i,oldsz,currentalsz
      type(t_scluster), ALLOCATABLE, INTENT(IN OUT):: it(:)
      type(t_scluster), ALLOCATABLE:: tmpcls(:)
    !
      oldsz = currentalsz
      allocate(tmpcls(oldsz))
      do i=1,oldsz
        call copy_cluster(it(i),tmpcls(i))
      end do
    !
      deallocate(it)
    !
      currentalsz = currentalsz*2
      allocate(it(currentalsz))
    !
      do i=1,oldsz
        call copy_cluster(tmpcls(i),it(i))
      end do
      it(oldsz+1:currentalsz)%nmbrs = 0
      it(oldsz+1:currentalsz)%nchildren = 0
      it(oldsz+1:currentalsz)%alsz = 0
      it(oldsz+1:currentalsz)%childr_alsz = 0
      it(oldsz+1:currentalsz)%parent = 0
    !
      deallocate(tmpcls)
    !
    end
    !
    !-------------------------------------------------------------------------------
    !


end module
