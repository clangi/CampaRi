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

  type t_scluster
    integer nmbrs, alsz
    real(KIND=4) sqsum
    real(KIND=4), ALLOCATABLE :: sums(:)
    integer, ALLOCATABLE:: snaps(:) ! snapshot lists
  end type t_scluster

  type(t_scluster), ALLOCATABLE :: scluster(:)
  real(KIND=4) radius
  integer nclalcsz ! number of clusters allocatable
  integer nclu ! number of clusters allocated

  ! This are the sizes for sums variable that keeps the centroid measures.
  ! In the rmsd (dis_method = 5) only needs sumssz = 1 and clstsz = calcsz
  ! (calcsz is the number of atoms * coordinates)
  contains
    subroutine leader_clustering(trj)
      use m_var_nbls_clu
      implicit none

      real(KIND=4), intent(in) :: trj(n_snaps,n_xyz)
      real(KIND=4) tmp_d ! temporary variable for radiuses and distances
      integer i, i2, j, ii, u
      real vecti(n_xyz)

      tmp_dis_method = dis_method(1) !This is hard. Can be eventually taylored

      nclu = 0 !cluster index
      do i=1,n_snaps
        ii = -1
        do j=nclu,max(nclu-500,1),-1 ! when nclu = 0 it skips this do  [j=0,1,-1]
          ! if(superver) write(*,*) "now checking dist(cluster,snap) : (",j,",",i,")"

          ! j is going down of -1, not nclu. nclu is growing if (ii.eq.-1)(new cluster is
          ! added). The number of clusters that can be considered in order to put i
          ! (snapshot) in one of those is 500. Once exceeded this upper bound
          ! it starts to shift on that limit forwards.
          if(tmp_dis_method.eq.5) then
            vecti = trj(i,1:n_xyz)
          else if(tmp_dis_method.eq.11) then
            ! veci(1:n_xyz) = 1
            i2 = i
            if(i.eq.1) i2 = 2 !to avoid error I can double the first value
            vecti = trj(i2,n_xyz) - trj(i2-1,n_xyz)
          end if
          call snap_to_cluster_d(tmp_d,scluster(j),vecti)
          if(superver .and. (i.lt.10)) write(*,*) "The distance is : ",tmp_d
          if (tmp_d.lt.radius) then
            call cluster_addsnap(scluster(j),vecti,i)
            ! add snapshot i in the cluster j
            ii = j
            exit
          end if
        end do
        if (ii.eq.-1) then
          ! a new cluster is added to the cluster list because tmp_d .gt. radius
          if(superver .and. (i.lt.10)) then
            write(*,*) "new cluster created: n", nclu + 1
          end if
          nclu = nclu + 1
          ! enlarge the size of the clusters array by the double
          if (nclu.gt.nclalcsz) then
            call cluster_lst_resize(scluster)
            ! allocate(tmpcls(nclalcsz))
            ! oldsz = nclalcsz
            ! tmpcls = scluster
            ! deallocate(scluster)
            ! nclalcsz = nclalcsz * 2
            ! allocate(scluster(nclalcsz))
            ! scluster(1:oldsz) = tmpcls(1:oldsz)
            ! scluster(oldsz+1:nclalcsz)%nmbrs = 0
            ! scluster(oldsz+1:nclalcsz)%alsz = 0
            ! deallocate(tmpcls)
          end if

          call cluster_addsnap(scluster(nclu),vecti,i) ! nclu new cluster
          ! if(superver) write(*,*) 'Cluster number ',nclu, ' created for snap ', i
        end if
      end do

      write(*,*) '... done with initial cluster generation for distances list with &
      &a total number of ', nclu , ' clusters'
      if(superver) then
        write(*,*) 'Number of elements per cluster: '
        do u=1,nclu
          write(*,*) 'Cl',u, ': nmbr of elements = ', scluster(u)%nmbrs
          ! write(*,*) 'Snaps contained: ', scluster(u)%snaps
        end do
      end if
    end subroutine
    subroutine cluster_lst_resize(clu)
      use m_var_nbls_clu
      implicit none
      type(t_scluster), ALLOCATABLE :: tmpcls(:)
      type(t_scluster), ALLOCATABLE :: clu(:)
      integer oldsz
      allocate(tmpcls(nclalcsz))
      oldsz = nclalcsz
      tmpcls = clu
      deallocate(clu)
      nclalcsz = nclalcsz * 2
      allocate(clu(nclalcsz))
      clu(1:oldsz) = tmpcls(1:oldsz)
      clu(oldsz+1:nclalcsz)%nmbrs = 0
      clu(oldsz+1:nclalcsz)%alsz = 0
      deallocate(tmpcls)

    end subroutine cluster_lst_resize

end module
