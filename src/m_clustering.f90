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
    !--------------------------------------------------------------------
    ! LEADER CLUSTERING
    !--------------------------------------------------------------------
    ! each snapshot is added to a cluster on the base of a defined threshold.
    subroutine leader_clustering(trj)
      use m_variables_gen
      implicit none

      real(KIND=4), intent(in) :: trj(n_snaps,n_xyz)
      real(KIND=4) tmp_d ! temporary variable for radiuses and distances
      integer i, j, ii, u
      real vecti(n_xyz)

      tmp_dis_method = dis_method(1) !This is hard. Can be eventually taylored

      nclu = 0 !cluster index
      do i=1,n_snaps
        ii = -1
        do j=nclu,max(nclu-500,1),-1 ! when nclu = 0 it skips this do  [j=0,1,-1]
          ! if(superver) write(ilog,*) "now checking dist(cluster,snap) : (",j,",",i,")"

          ! j is going down of -1, not nclu. nclu is growing if (ii.eq.-1)(new cluster is
          ! added). The number of clusters that can be considered in order to put i
          ! (snapshot) in one of those is 500. Once exceeded this upper bound
          ! it starts to shift on that limit forwards.
          if(tmp_dis_method.eq.5) then
            vecti = trj(i,1:n_xyz)
          else if(tmp_dis_method.eq.11) then
            if(i.eq.1) then !to avoid error I can double the first value
               vecti = trj(i+1,n_xyz) - trj(i,n_xyz)
            else
              vecti = trj(i,n_xyz) - trj(i-1,n_xyz)
            end if
          end if
          call snap_to_cluster_d(tmp_d,scluster(j),vecti)
          if(superver .and. (i.lt.10)) write(ilog,*) "The distance is : ",tmp_d
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
            write(ilog,*) "new cluster created: n", nclu + 1
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
          ! if(superver) write(ilog,*) 'Cluster number ',nclu, ' created for snap ', i
        end if
      end do

      write(ilog,*) '... done with initial cluster generation for distances list with &
      &a total number of ', nclu , ' clusters'
      if(superver) then
        write(ilog,*) 'Number of elements per cluster: '
        do u=1,nclu
          write(ilog,*) 'Cl',u, ': nmbr of elements = ', scluster(u)%nmbrs
          ! write(ilog,*) 'Snaps contained: ', scluster(u)%snaps
        end do
      end if
    end subroutine
    subroutine cluster_lst_resize(clu)
      use m_variables_gen
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
    !--------------------------------------------------------------------
    ! BIRCH CLUSTERING
    !--------------------------------------------------------------------
    ! it is based on a kind of unbalanced BIRCH algorithm
    ! subroutine birch_clustering(modei,nnodes)
    ! !
    !   use m_variables_gen
    ! !
    !   implicit none
    ! !
    !   integer i,j,k,ii,jj,kk,ll,mm,atwo,fail,kkf,thekk,nlst1,nlst2,snapstart,snapend,snapinc,nnodes,modei
    !   integer, ALLOCATABLE:: kklst(:,:),errcnt(:),kkhistory(:)
    !   integer(KIND=8) cnt1,cnt2
    !   RTYPE, ALLOCATABLE:: scrcts(:)
    !   RTYPE rdv,mind,helper,maxd(2),normer(4),qualmet(4)
    !   logical atrue,afalse,notdone
    ! !
    !  33 format('ERROR B: ',i6,' could be part of ',i5,' at ',g14.7,' (last d',g14.6,')')
    !  34 format('NEXT HIGHER: ',i6,' could have been part of ',i5,' at ',g14.7,'.')
    !  35 format('ERROR C: ',i6,' is meant to a child of ',i5,' at level ',i5,' but has distance ',g14.7,'.')
    ! !
    !   atrue = .true.
    !   afalse = .false.
    !   atwo = 2
    !   if (cmaxrad.le.cradius) cmaxrad = 2.0*cradius
    ! !
    !   if (cleadermode.le.2) then
    !     snapstart = 1
    !     snapend = cstored
    !     snapinc = 1
    !   else
    !     snapstart = cstored
    !     snapend = 1
    !     snapinc = -1
    !   end if
    ! !
    !   allocate(kklst(cstored,2))
    !   allocate(kkhistory(c_nhier+1))
    !   allocate(errcnt(c_nhier+1))
    !   allocate(birchtree(c_nhier+1))
    !   do ii=1,c_nhier+1
    !     allocate(birchtree(ii)%cls(10))
    !     birchtree(ii)%ncls = 0
    !     birchtree(ii)%nclsalsz = 10
    !   end do
    !   allocate(scrcts(c_nhier+1))
    !   errcnt(:) = 0
    !   birchtree(1)%ncls = 1
    !   scrcts(c_nhier+1) = cradius
    !   if (c_nhier.gt.1) then
    !     scrcts(2) = cmaxrad
    !     do i=3,c_nhier
    !       scrcts(i) = cmaxrad - ((i-2.0)/(c_nhier-1.0))*(cmaxrad - cradius) ! linear so far
    !     end do
    !   end if
    !   do i=1,c_nhier+1
    !     birchtree(i)%cls(:)%nmbrs = 0
    !     birchtree(i)%cls(:)%nb = 0
    !     birchtree(i)%cls(:)%nchildren = 0
    !     birchtree(i)%cls(:)%alsz = 0
    !     birchtree(i)%cls(:)%nbalsz = 0
    !     birchtree(i)%cls(:)%chalsz = 0
    !     birchtree(i)%cls(:)%parent = 0
    !   end do
    !
    !   allocate(birchtree(1)%cls(1)%snaps(2))
    !   birchtree(1)%cls(1)%alsz = 2
    !   allocate(birchtree(1)%cls(1)%children(2))
    !   birchtree(1)%cls(1)%chalsz = 2
    ! !
    !   write(ilog,*)
    !   write(ilog,*) 'Now performing tree-based clustering ...'
    !   cnt1 = 0
    !   cnt2 = 0
    !   do i=snapstart,snapend,snapinc
    !     kk = 1
    !     notdone = .true.
    !     fail = -1
    !     nlst1 = 1
    !     kklst(1,1) = kk
    ! !    cnt1 = 0
    !     do ii=2,c_nhier
    !       jj = -1
    !       mind = HUGE(mind)
    !       nlst2 = 0
    !       do mm=1,nlst1
    !         kk = kklst(mm,1)
    !         do j=1,birchtree(ii-1)%cls(kk)%nchildren
    !           ll = birchtree(ii-1)%cls(kk)%children(j)
    !           cnt1 = cnt1 + 1
    !           if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
    !           call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
    !           if (rdv.lt.mind) then
    !             mind = rdv
    !             jj = j
    !             thekk = kk
    !           end if
    !           if (rdv.lt.scrcts(ii)) then
    !             nlst2 = nlst2 + 1
    !             kklst(nlst2,2) = ll
    !           end if
    !         end do
    !       end do
    ! !     store the path
    !       if (jj.eq.-1) then
    !         kkhistory(ii-1) = 1
    !       else
    !         kkhistory(ii-1) = thekk
    !       end if
    !       if ((ii.le.(c_nhier-1)).AND.(jj.gt.0)) then
    !         if (nlst2.le.0) then ! absolutely nothing nearby
    !           nlst1 = 1
    !           kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
    !           if (fail.eq.-1) then
    ! !            write(ilog,*) i,' failed at ',ii,' w/ ',mind
    !             fail = ii
    !             kkf = thekk
    !           end if
    !         else
    !           if (fail.gt.0) then
    !             do j=fail,ii-1
    !               birchtree(j)%ncls = birchtree(j)%ncls + 1
    !               if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
    !               kkhistory(j) = birchtree(j)%ncls
    !               if (j.gt.fail) then
    !                 call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
    !  &                                    birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
    !                 cnt2 = cnt2 + 1
    !               end if
    !             end do
    !             call cluster_addchild(birchtree(fail-1)%cls(kkf),kkf,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
    !             call cluster_addchild(birchtree(ii-1)%cls(birchtree(ii-1)%ncls),&
    !  &                                birchtree(ii)%cls(jj)%parent,birchtree(ii)%cls(jj),jj)
    !             fail = -1
    !           end if
    !           nlst1 = 1
    !           kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
    !         end if
    !         cycle
    ! !     leaf
    !       else if ((ii.eq.c_nhier).AND.(mind.lt.scrcts(ii))) then
    !         kk = birchtree(ii-1)%cls(thekk)%children(jj)
    !         call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
    !         do j=2,c_nhier
    !           call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),i,rdv)
    !         end do
    !         notdone = .false.
    !       else
    !         if (fail.eq.-1) then
    !           fail = ii
    !           kkf = kk
    !         end if
    !         kk = kkf
    !         do j=fail,c_nhier
    !           birchtree(j)%ncls = birchtree(j)%ncls + 1
    !           if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
    !           call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
    !           if (j.gt.fail) then
    !             call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
    !  &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
    !             cnt2 = cnt2 + 1
    !           end if
    !         end do
    !         call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
    !         cnt2 = cnt2 + 1
    !         do j=2,fail
    !           call cluster_addsnap(birchtree(j-1)%cls(kkhistory(j-1)),i,rdv)
    !         end do
    !         call snap_to_cluster_d(maxd(1),birchtree(fail-1)%cls(kk),i)
    !         notdone = .false.
    !       end if
    !       if (notdone.EQV..false.) exit
    !     end do
    !   end do
    !
    !   do i=snapstart,snapend,snapinc
    !     kk = 1
    !     notdone = .true.
    !     fail = -1
    !     nlst1 = 1
    !     kklst(1,1) = kk
    !     do ii=2,c_nhier+1
    !       jj = -1
    !       mind = HUGE(mind)
    !       nlst2 = 0
    !       do mm=1,nlst1
    !         kk = kklst(mm,1)
    !         do j=1,birchtree(ii-1)%cls(kk)%nchildren
    !           ll = birchtree(ii-1)%cls(kk)%children(j)
    !           cnt1 = cnt1 + 1
    !           if ((birchtree(ii)%cls(ll)%center.le.0).OR.(birchtree(ii)%cls(ll)%center.gt.cstored)) call fexit()
    !           call snap_to_cluster_d(rdv,birchtree(ii)%cls(ll),i)
    !           if (rdv.lt.mind) then
    !             mind = rdv
    !             jj = j
    !             thekk = kk
    !           end if
    !           if (rdv.lt.scrcts(ii)) then
    !             nlst2 = nlst2 + 1
    !             kklst(nlst2,2) = ll
    !           end if
    !         end do
    !       end do
    ! !     store the path
    !       if (jj.eq.-1) then
    !         kkhistory(ii-1) = 1
    !       else
    !         kkhistory(ii-1) = thekk
    !       end if
    !       if ((ii.le.c_nhier).AND.(jj.gt.0)) then
    !         if (nlst2.le.0) then ! absolutely nothing nearby
    !           nlst1 = 1
    !           kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
    !           if (fail.eq.-1) then
    ! !            write(ilog,*) i,' failed at ',ii,' w/ ',mind
    !             fail = ii
    !             kkf = thekk
    !           end if
    !         else
    !           nlst1 = 1
    !           kklst(1:nlst1,1) = birchtree(ii-1)%cls(thekk)%children(jj)
    !         end if
    !         cycle
    ! !     leaf
    !       else if ((ii.eq.(c_nhier+1)).AND.(mind.lt.scrcts(ii))) then
    !         kk = birchtree(ii-1)%cls(thekk)%children(jj)
    !         call cluster_addsnap(birchtree(ii)%cls(kk),i,rdv)
    !         notdone = .false.
    !       else if (fail.gt.0) then
    !         kk = kkf
    !         do j=fail,c_nhier+1
    !           birchtree(j)%ncls = birchtree(j)%ncls + 1
    !           if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
    !           call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
    !           if (j.gt.fail) then
    !             call cluster_addchild(birchtree(j-1)%cls(birchtree(j-1)%ncls),birchtree(j-1)%ncls,&
    !  &                                birchtree(j)%cls(birchtree(j)%ncls),birchtree(j)%ncls)
    !                      cnt2 = cnt2 + 1
    !           end if
    !         end do
    !         call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
    !         cnt2 = cnt2 + 1
    !         notdone = .false.
    !       else
    !         fail = ii
    !         j = c_nhier+1
    !         birchtree(j)%ncls = birchtree(j)%ncls + 1
    !         if (birchtree(j)%ncls.gt.birchtree(j)%nclsalsz) call scluster_resizelst(birchtree(j)%nclsalsz,birchtree(j)%cls)
    !         call cluster_addsnap(birchtree(j)%cls(birchtree(j)%ncls),i,rdv)
    !         call cluster_addchild(birchtree(fail-1)%cls(kk),kk,birchtree(fail)%cls(birchtree(fail)%ncls),birchtree(fail)%ncls)
    !         notdone = .false.
    !       end if
    !       if (notdone.EQV..false.) exit
    !     end do
    !   end do
    !
    !  66 format('Level    # Clusters     Threshold     Total Snaps    Total Children')
    !  67 format(i9,i10,1x,g14.4,4x,i12,4x,i12)
    !  68 format(i9,i10,5x,a7,7x,i12,4x,i12)
    !   write(ilog,66)
    !   write(ilog,68) c_nhier+1,birchtree(1)%ncls,'MAXIMAL',sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nmbrs),&
    !  &               sum(birchtree(1)%cls(1:birchtree(1)%ncls)%nchildren)
    !   do i=2,c_nhier+1
    !     write(ilog,67) c_nhier+2-i,birchtree(i)%ncls,scrcts(i),sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nmbrs),&
    !  &               sum(birchtree(i)%cls(1:birchtree(i)%ncls)%nchildren)
    !   end do
    !   write(ilog,*) '---------------------------------------------------------------------'
    !   write(ilog,*)
    !   write(ilog,*) '... done after a total of ',cnt1,' distance evaluations.'
    !   write(ilog,*)
    ! !
    !   do j=1,birchtree(c_nhier+1)%ncls
    !     call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
    !   end do
    ! !
    !   if (refine_clustering.EQV..true.) then
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
    !         call cluster_to_cluster_d(rdv,birchtree(c_nhier)%cls(i),birchtree(c_nhier)%cls(kkf))
    !         cnt1 = cnt1 + 1
    !         if (rdv.lt.scrcts(c_nhier)) then
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
    !             call clusters_joint_diam(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk),rdv,helper)
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
    !             if ((helper.le.maxd(2)).OR.(rdv.le.maxd(1))) then
    !               if (birchtree(c_nhier+1)%cls(jj)%nmbrs.gt.birchtree(c_nhier+1)%cls(kk)%nmbrs) then
    !                 call join_clusters(birchtree(c_nhier+1)%cls(jj),birchtree(c_nhier+1)%cls(kk))
    !                 birchtree(c_nhier+1)%cls(jj)%diam = rdv
    !                 birchtree(c_nhier+1)%cls(jj)%radius = helper
    !                 cnt2 = cnt2 + 1
    !               else
    !                 call join_clusters(birchtree(c_nhier+1)%cls(kk),birchtree(c_nhier+1)%cls(jj))
    !                 birchtree(c_nhier+1)%cls(kk)%diam = rdv
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
    ! !
    !   do j=1,birchtree(c_nhier+1)%ncls
    !     call cluster_calc_params(birchtree(c_nhier+1)%cls(j),scrcts(c_nhier+1))
    !   end do
    !
    ! ! now shorten list and resorti
    !   call clusters_shorten(birchtree(c_nhier+1)%cls,birchtree(c_nhier+1)%ncls)
    !   call clusters_sort(birchtree(c_nhier+1)%cls,birchtree(c_nhier+1)%ncls,afalse)
    !   do i=1,birchtree(c_nhier+1)%ncls
    !     call cluster_getcenter(birchtree(c_nhier+1)%cls(i))
    !   end do
    ! ! lastly, copy into global cluster array
    !   allocate(scluster(birchtree(c_nhier+1)%ncls))
    !   do i=1,birchtree(c_nhier+1)%ncls
    !     call copy_cluster(birchtree(c_nhier+1)%cls(i),scluster(i))
    !   end do
    !   nnodes = birchtree(c_nhier+1)%ncls
    !   do k=2,c_nhier+1
    !     do i=1,birchtree(k)%ncls
    !       call cluster_calc_params(birchtree(k)%cls(i),scrcts(k))
    !     end do
    !   end do
    ! !
    !  63 format(i7,1x,i7,1x,i8,1000(1x,g12.5))
    !  64 format(1000(g12.5,1x))
    ! !
    !   write(ilog,*) '------------- CLUSTER SUMMARY ------------------'
    !   write(ilog,*) ' #       No.     "Center"  Diameter     Radius      '
    !   do i=1,birchtree(c_nhier+1)%ncls
    !     write(ilog,63) i,scluster(i)%nmbrs,scluster(i)%center,scluster(i)%diam,scluster(i)%radius
    !   end do
    !   write(ilog,*) '------------------------------------------------'
    !   write(ilog,*)
    ! !
    !   atrue = .true.
    !   call quality_of_clustering(nnodes,scluster,scrcts(c_nhier+1),qualmet)
    !   call gen_graph_from_clusters(scluster,nnodes,atrue)
    !   call graphml_helper_for_clustering(scluster,nnodes)
    !   call vmd_helper_for_clustering(scluster,nnodes)
    ! !
    !   deallocate(kklst)
    !   deallocate(kkhistory)
    !   deallocate(scrcts)
    !   deallocate(errcnt)
    ! !
    ! end
    !
    !
end module
