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

! CHANGES IN VARIABLES
! dis_method = cdis_crit
! n_xyz = cdofsbnds !principally the length and starting point of the coordinate
! n_snaps = cstored
! trj_data(n_snaps,n_xyz) = cludata(n_xyz,n_snaps)
! tmp_d = rmsdtmp
! n_xyz = calcsz
! childr_alsz = chalsz
! I eliminated maxrads and the |clust(1|2)| = 1 distance for redundancy
!
! OLDERIES
! real, ALLOCATABLE :: maxrads(:) !vector of radius for each cluster
! This is a different variable from cluster radius because it is needed that
! cluster is complete to calculate it.
! I decided not to use the maxrads because it is not a useful calculation
! if the data-set is big enough to have cluster that has maxrads ~ radius
! dis_method = 'FMCSC_CDISTANCE'
!
!
subroutine generate_neighbour_list( &
  trj_data, n_xyz_in, n_snaps_in, dfffo, clu_radius_in, clu_hardcut_in, & !input
  adjl_deg, adjl_ix, adjl_dis, max_degr, & !output
  dis_method_in, birch_in, mst_in, & !algorithm details
  rootmax_rad_in, tree_height_in, n_search_attempts_in,& !sst details
  normalize_dis_in, return_tree_in_r, mute_in) !modes

  use gutenberg
  use m_variables_gen
  use m_clustering
  use m_gen_nbls
#ifdef LINK_NETCDF
  use netcdf
  use m_mst_dumping
#else
  use m_mst !   use m_mst_dumping
#endif
  implicit none

  ! INPUT VARIABLES
  integer, intent(in) :: n_xyz_in !numbers of xyz (atoms*3)
  integer, intent(in) :: n_snaps_in !number of snapshots in input
  integer, intent(in) :: dfffo !dimensional_flag_for_fixed_out (netcdf workaround)
  integer, intent(in) :: dis_method_in !distance method
  integer, intent(in) :: tree_height_in !number of levels in the tree
  integer, intent(in) :: n_search_attempts_in !number of search attempts in sst
  real, intent(in) :: trj_data(n_snaps_in, n_xyz_in) !trajectory input
  real, intent(in) :: clu_radius_in !defines the max cluster sizes
  real, intent(in) :: clu_hardcut_in !threshold between cluster snaps
  real, intent(in) :: rootmax_rad_in !first level (non-root) rad-thresh
  logical, intent(in) :: birch_in !flag for birch clustering
  logical, intent(in) :: mst_in !make already the ordered mst(true) or not?
  logical, intent(in) :: normalize_dis_in !flag for normalize the distance matrix
  logical, intent(in) :: return_tree_in_r !if I want to use R even if with netcdf..
  logical, intent(in) :: mute_in !verbose terminal output
  ! the mst option could generate an adjlist which was not an mst - only NO BIRCH

  ! HELPING AND DEBUGGING VARIABLES
  integer i ! helper variables
  integer nzeros ! nzeros = number of not connected components (dis .le. 0)
  ! character(len=1024) :: format_var !format for above mentioned dumping
  real t2,t1 !timing variables


  ! OUTPUT VARIABLES !strict output collaboration with R ! TODO: dummy var to lower memory need
  integer, intent(inout) :: max_degr !maximum degree of the adjlist
  integer, intent(inout) :: adjl_deg(dfffo)
  integer, intent(inout) :: adjl_ix(dfffo,dfffo)
  real, intent(inout) :: adjl_dis(dfffo,dfffo)

  ! CHEKS IF THE USER HAS WHAT HE WANTS - NETCDF
  ! ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---
#ifdef LINK_NETCDF
  if(return_tree_in_r) then
    call sl()
    call spr('------------------------------------------------------------------')
    call spr('ATTENTION: Even if CampaRi was installed using netcdf support, you selected &
    to use the R data management system. Both options will be followed at the same time.')
    call spr('------------------------------------------------------------------')
  end if
#else
  if(.not.return_tree_in_r) then
    call sl()
    call spr('------------------------------------------------------------------')
    call spr('ATTENTION: Even if CampaRi was installed without the netcdf support, &
    the user tried to use the netcdf dumping functionality. This run will follow the &
    usual flow without netcdf. If you want to use it install the full version of the package.')
    call spr('------------------------------------------------------------------')
  end if
#endif
  ! ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---

  ! N B :
  ! clu_hardcut is used for distances between different clusters snapshos as
  ! a threshold in leader clusting (MST).
  ! If the intent is in this cannot be an ALLOCATABLE variable

  n_xyz = n_xyz_in
  n_snaps = n_snaps_in
  mute = mute_in
  radius = clu_radius_in
  dis_method = dis_method_in
  birch = birch_in
  do_mst = mst_in
  normalize_dis = normalize_dis_in

  ! Internal defaults
  if (hardcut.lt.radius) hardcut = 2.0*radius
  firstcall = .true.

  !debugging flags
  rand_seed = 10 !it is used to have fixed seed
  !if it is 0 it uses the standard random_seed
  clu_summary = .true. !showing or not the clu summary
  precise_clu_descr = 10 !dev var that shows the first 10 clusters

  !SST
  cmaxrad = rootmax_rad_in !root level size threshold
  c_nhier = tree_height_in !height of the tree
  ! default vars that should be set from main function TODO
  ! cprogindrmax = max(floor((n_snaps*7.0)/100),min(n_snaps,5)) !def 7/100
  cprogindrmax = n_search_attempts_in
  cprogbatchsz = 1  !  batch size for random stretches aka dim of random branches def = 1 TODO
  cprogrdepth = 0   !  auxiliary search depth def = 0 TODO
  if(cprogrdepth.gt.c_nhier) cprogrdepth = c_nhier
  if(cmaxrad.le.radius) cmaxrad = 2.0*radius
  c_multires = 0 !inital value of FMCSC_BIRCHMULTI TODO
  ordering = 1 !TODO


  !LOGGING id
  ! defining the logging id. 0 is sterror, 5 stinput, 6 stoutput
  ! ilog = 6

  ! Logging function - LOGGING is R-sided now
  ! if(log_print) then
  !   inquire(file="campari.log", exist=exist)
  !   if (exist) then
  !     open(ilog, file="campari.log", status="replace", action="write")
  !     write(ilog,*) 'NB: campari.log already exists. It has been overwritten'
  !   else
  !     open(ilog, file="campari.log", status="new", action="write")
  !   end if
  ! end if


  call sl()
  call sl()
#ifdef LINK_NETCDF
  call spr("Starting analysis in Fortran [with netcdf dumping]...")
#else
  call spr("Starting analysis in Fortran [no netcdf dumping]...")
#endif
  call spr('------------------------------------------------------------')
  call sipr("Selected distance: ",dis_method)
  call sl()
  call sipr("Number of snapshots: ",n_snaps)
  call sipr("Number of variables: ",n_xyz)
  call sl()

  ! if(superver) then
  !    write(ilog,*) "Input example", trj_data(1:10,1:10)
  !    write(ilog,*)
  !  end if


  if(.not.birch) then
    call CPU_time(t1)
    call leader_clustering(trj_data)
    ! now compare all blocks to each other (the slowest part) taking advantage
    ! of information generated previously (otherwise intractable)
    call spr('------------------------------------------------------------')
    call spr( 'Now computing cutoff-assisted neighbor list...')
    call sl()
    call gen_nb(trj_data)
    call CPU_time(t2)
    call srpr( 'Time elapsed for neighbor list creation/reading (s): ', t2-t1)
    call sl()
    call spr( 'Neighbor list generated.')
    call sl()
      nzeros = 0
      do i=1,n_snaps
        if (cnblst(i)%nbs.le.0) then
          nzeros = nzeros + 1
          ! call sipr('Warning. Snapshot # ',i,' is with 0 degree'
       end if
      end do
      if (nzeros.gt.0) then
        call sipr( 'Warning. The following snapshots are without a neighbor &
        &(similar) structure. This may in some cases cause the clustering &
        &algorithm to misbehave.',nzeros)
      end if
      call sl()
      do i=1,n_clu_alc_sz_gen
        if (allocated(scluster(i)%snaps).EQV..true.) deallocate(scluster(i)%snaps)
        if (allocated(scluster(i)%sums).EQV..true.) deallocate(scluster(i)%sums)
        if (allocated(scluster(i)%tmpsnaps).EQV..true.) deallocate(scluster(i)%tmpsnaps)
        if (allocated(scluster(i)%children).EQV..true.) deallocate(scluster(i)%children)
      end do
      deallocate(scluster)

#ifdef LINK_NETCDF
      call CPU_time(t1)
      call gen_MST_from_nbl_w()
      call CPU_time(t2)
      if(.not.mute) call intpr( 'TIME elapsed for MST building (s): ',-1,t2-t1,1)
#else
      if(do_mst) then
        call CPU_time(t1)
        call gen_MST_from_nbl(adjl_deg,adjl_ix,adjl_dis,max_degr)
        call CPU_time(t2)
        call srpr( 'TIME elapsed for MST building (s): ',t2-t1)
      else
        do i=1,n_snaps
          adjl_deg(i) = cnblst(i)%nbs
          adjl_ix(i,:) = cnblst(i)%idx
          adjl_dis(i,:) = cnblst(i)%dis
          ! the max search is again needed if not MST is built
          if(i .eq. 1 .OR. max_degr .lt. cnblst(i)%nbs) max_degr = cnblst(i)%nbs
          if (allocated(cnblst(i)%dis).EQV..true.) deallocate(cnblst(i)%dis)
          if (allocated(cnblst(i)%idx).EQV..true.) deallocate(cnblst(i)%idx)
        end do
        deallocate(cnblst)
      end if
#endif
  else
    call CPU_time(t1)
    call birch_clustering(trj_data)
    call CPU_time(t2)
    call srpr('Time elapsed for birch_clustering (s): ',t2-t1)
    call spr('Birch clustering completed.')
    call spr('------------------------------------------------------------')
    call CPU_time(t1)
#ifdef LINK_NETCDF
    call gen_MST_from_treeclustering_w(trj_data)
#else
    call gen_MST_from_treeclustering(adjl_deg, adjl_ix, adjl_dis, max_degr, trj_data)
#endif
    call CPU_time(t2)
    call srpr('Time elapsed for SST building (s): ',t2-t1)
    call spr('SST generated from tree-based clustering successfully.')
    call spr('------------------------------------------------------------')
    call sl()
  end if

  ! Dumping to netcdf if necessary
#ifdef LINK_NETCDF
  call spr('DUMPING THE MST/SST to netcdf...')
  call sl()
  ! #ifdef LINK_NETCDF
  call CPU_time(t1)
  call dump_nbl_nc()
  call CPU_time(t2)
  call sl()
  call srpr('Time elapsed for mst dumping (s): ',t2-t1)
  call spr('...file-mst dumping done.')
  if(return_tree_in_r) then
    call sl()
    call spr('As you have the R-based tree handling active, now the tree &
    will return to the glorious hands of R...')
    do i=1,n_snaps
      adjl_deg(i) = approxmst(i)%deg
      adjl_ix(i,:) = approxmst(i)%adj(1:approxmst(i)%deg)
      adjl_dis(i,:) = approxmst(i)%dist(1:approxmst(i)%deg)
      if(i .eq. 1 .OR. max_degr .lt. approxmst(i)%deg) max_degr = approxmst(i)%deg
    end do
    call spr('...done')
  end if
  ! deallcoating the sst used for dumping to netcdf
  if (allocated(approxmst).EQV..true.) then
    do i=1,size(approxmst)
      if (allocated(approxmst(i)%adj).EQV..true.) deallocate(approxmst(i)%adj)
      if (allocated(approxmst(i)%dist).EQV..true.) deallocate(approxmst(i)%dist)
    end do
    deallocate(approxmst)
  end if
  call spr('------------------------------------------------------------')
#endif
  ! if(log_print) close(ilog)
  firstcall = .false.

end
