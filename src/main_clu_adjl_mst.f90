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
! tmp_dis_method = cdis_crit
! n_xyz = cdofsbnds !principally the length and starting point of the coordinate
! n_snaps = cstored
! trj_data(n_snaps,n_xyz) = cludata(n_xyz,n_snaps)
! tmp_d = rmsdtmp
! n_xyz = calcsz
! childr_alsz = chalsz
! I eliminated maxrads and the |clust(1|2)| = 1 distance for redundancy
!
! OLDERIES
! real(KIND=4), ALLOCATABLE :: maxrads(:) !vector of radius for each cluster
! This is a different variable from cluster radius because it is needed that
! cluster is complete to calculate it.
! I decided not to use the maxrads because it is not a useful calculation
! if the data-set is big enough to have cluster that has maxrads ~ radius
! dis_method = 'FMCSC_CDISTANCE'
!
!
subroutine generate_neighbour_list( &
  trj_data, n_xyz_in, n_snaps_in, clu_radius_in, clu_hardcut_in, & !input
  adjl_deg, adjl_ix, adjl_dis, max_degr, & !output
  dis_method_in, dis_weight_in, birch_in, mst_in, & !algorithm details
  rootmax_rad_in, tree_height_in, n_search_attempts_in,& !sst details
  data_meth_in, normalize_dis_in, log_print_in, verbose_in) !modes

  use m_variables_gen
  use m_clustering
  use m_gen_nbls
  use m_mst
  implicit none

  integer, INTENT(IN) :: n_xyz_in !numbers of xyz (atoms*3)
  integer, INTENT(IN) :: n_snaps_in !number of snapshots in input
  integer, intent(in) :: dis_method_in(11) !distance method
  integer, intent(in) :: data_meth_in !data managing method TODO netcdf
  integer, intent(in) :: tree_height_in !number of levels in the tree
  integer, intent(in) :: n_search_attempts_in !number of search attempts in sst
  real(KIND=4), INTENT(IN) :: trj_data(n_snaps_in, n_xyz_in) !trajectory input
  real(KIND=4), intent(in) :: dis_weight_in(11) !distance weight when averaged
  real(KIND=4), intent(in) :: clu_radius_in !defines the max cluster sizes
  real(KIND=4), intent(in) :: clu_hardcut_in !threshold between cluster snaps
  real(KIND=4), intent(in) :: rootmax_rad_in !first level (non-root) rad-thresh
  logical, intent(in) :: verbose_in !verbose terminal output
  logical, intent(in) :: mst_in !make already the ordered mst(true) or not?
  logical, intent(in) :: normalize_dis_in !flag for normalize the distance matrix
  logical, intent(in) :: birch_in !flag for birch clustering
  logical, intent(in) :: log_print_in !flag for printing log file or on terminal
  !if the intent is in this cannot be an ALLOCATABLE variable


  ! N B :
  ! clu_hardcut is used for distances between different clusters snapshos as
  ! a threshold.

  integer i,j,ii,kk,ll,k,l,mi,mj,u,u2 ! helper variables
  integer nzeros ! nzeros = number of not connected components (dis .le. 0)
  integer mst_file_unit, freeunit ! must have for using the same-name-function
  logical exist, mst_print, log_print !file dumping for mst tree for debugging
  character(len=1024) :: format_var !format for above mentioned dumping
  real t2,t1 !timing variables

  integer, intent(inout) :: max_degr !maximum degree of the adjlist
  integer, intent(inout) :: adjl_deg(n_snaps_in)
  integer, intent(inout) :: adjl_ix(n_snaps_in,n_snaps_in)
  real(KIND=4), intent(inout) :: adjl_dis(n_snaps_in,n_snaps_in)
  ! real(KIND=4), intent(inout) :: adj_mat(:,:) !adjac matrix
  !to have an inout intent the thing cant be allocatable...



  ! Intents are not allocatable
  n_xyz = n_xyz_in
  n_snaps = n_snaps_in
  ver = verbose_in
  log_print = log_print_in !logging flag
  radius = clu_radius_in
  dis_method = dis_method_in
  birch = birch_in
  mst = mst_in
  normalize_dis = normalize_dis_in

  !defaults
  dis_weight = 1
  n_dis_method = 0
  if (hardcut.lt.radius) hardcut = 2.0*radius

  !debugging flag
  superver = .false.  !dev flag for superverbose output
  mst_print = .false. !dev flag for adjlist (mst) checking
  rand_seed_cnt = 1 !it is used to have fixed seed (different for each call)
  !if it is 0 it uses the standard random_seed

  !SST
  cmaxrad = rootmax_rad_in !root level size threshold
  c_nhier = tree_height_in
  ! default vars that should be set from main function
  ! cprogindrmax = max(floor((n_snaps*7.0)/100),min(n_snaps,5)) !def 7/100
  cprogindrmax = n_search_attempts_in
  cprogbatchsz = 1  !  batch size for random stretches aka dim of random branches def = 1 TODO
  cprogrdepth = 0   !  auxiliary search depth def = 0 TODO
  if(cprogrdepth.gt.c_nhier) cprogrdepth = c_nhier


! defining the logging id. 0 is sterror, 5 stinput, 6 stoutput
  if(log_print) then
    ilog = freeunit()
  else
    ilog = 6
  end if


  ! Logging function
  if(log_print) then
    inquire(file="campari.log", exist=exist)
    if (exist) then
      open(ilog, file="campari.log", status="replace", action="write")
      write(ilog,*) 'NB: campari.log already exists. It has been overwritten'
    else
      open(ilog, file="campari.log", status="new", action="write")
    end if
  end if

  write(ilog,*)
  write(ilog,*) '-----------------------------------'
  write(ilog,*) 'WELCOME TO CAMPARI ANALYSIS TOOL'
  write(ilog,*) '-----------------------------------'
  write(ilog,*)

  do i=1,11
    if(dis_method(i).ge.1.and.dis_method(i).le.11) n_dis_method = n_dis_method + 1
    if(dis_weight_in(i).ge.0.and.dis_weight_in(i).le.1) dis_weight(i) = dis_weight_in(i)
  end do
  write(ilog,*)
  if(n_dis_method.gt.1) write(ilog,*) 'Multiple distance functions selected. They will be &
  & balanced using clustering valuse, then summed toghether in a old fashion way.\n&
  &Remember that the leader clustering algorithm will consider only the first value &
  &of the distances'
  write(ilog,*)
  write(ilog,*) "Selected distances: ", dis_method(1:n_dis_method)
  write(ilog,*) "Distance weights:", dis_weight(1:n_dis_method)
  write(ilog,*)

  allocate(cnblst(n_snaps))
  if(n_dis_method.gt.1) allocate(cnblst_dis(n_snaps))


  cnblst(:)%nbs = 0 ! number of snapshots that are connected to one snap
  cnblst(:)%alsz = 4 ! allocation size
  ! maxalcsz = 4 ! default variable for the total max allocation size
  do i=1,n_snaps
    allocate(cnblst(i)%idx(cnblst(i)%alsz))
    allocate(cnblst(i)%dis(cnblst(i)%alsz))
    if(n_dis_method.gt.1) allocate(cnblst_dis(i)%dis(cnblst(i)%nbs))
  end do

  write(ilog,*) "Input dimensions:", n_snaps, " row and ", n_xyz," col"
  write(ilog,*)
  if(superver) then
     write(ilog,*) "clu_radius_in", clu_radius_in
     write(ilog,*) "clu_hardcut_in", clu_hardcut_in
     write(ilog,*)
     write(ilog,*) "Input example", trj_data(1:10,1:10)
     write(ilog,*)
   end if

  write(ilog,*)


  if(.not.birch) then
    call CPU_time(t1)
    call leader_clustering(trj_data)

    ! now compare all blocks to each other (the slowest part) taking advantage
    ! of information generated previously (otherwise intractable)
    write(ilog,*) '-----------------------------------'
    write(ilog,*) 'Now computing cutoff-assisted neighbor list...'
    write(ilog,*)

    call gen_nb(trj_data)
    call CPU_time(t2)
    write(ilog,*) 'TIME elapsed for neighbor list creation/reading: ',t2-t1, ' [s]'

    write(ilog,*)
    write(ilog,*) 'Neighbor list generated.'
    write(ilog,*)

      nzeros = 0
      do i=1,n_snaps
        if (cnblst(i)%nbs.le.0) then
          nzeros = nzeros + 1
         if(superver) write(ilog,*) 'Warning. Snapshot # ',i,' is without a &
         &neighbor (similar) structure. This may cause the clustering algorithm &
         &to crash or misbehave otherwise.'
       end if
      end do
      if (nzeros.gt.0) then
        write(ilog,*) 'Warning. ',nzeros,' snapshots are without a neighbor &
        &(similar) structure. This may in some cases cause the clustering &
        &algorithm to misbehave.'
      end if
      write(ilog,*)

      do i=1,n_clu_alc_sz_gen
        if (allocated(scluster(i)%snaps).EQV..true.) deallocate(scluster(i)%snaps)
        if (allocated(scluster(i)%sums).EQV..true.) deallocate(scluster(i)%sums)
      end do
      deallocate(scluster)

      if(mst) then
        call gen_MST_from_nbl(adjl_deg,adjl_ix,adjl_dis,max_degr)
        call CPU_time(t2)
        write(ilog,*) 'TIME elapsed for MST: ',t2-t1, ' [s]'
      else
        do i=1,n_snaps
          adjl_deg(i) = cnblst(i)%nbs
          adjl_ix(i,:) = cnblst(i)%idx
          adjl_dis(i,:) = cnblst(i)%dis
          if(i .eq. 1 .OR. max_degr .lt. cnblst(i)%nbs) max_degr = cnblst(i)%nbs
          if (allocated(cnblst(i)%dis).EQV..true.) deallocate(cnblst(i)%dis)
          if (allocated(cnblst(i)%idx).EQV..true.) deallocate(cnblst(i)%idx)
        end do
        deallocate(cnblst)
      end if
  else
    call CPU_time(t1)
    call birch_clustering(trj_data)
    call CPU_time(t2)
    write(ilog,*) '>>TIME<< elapsed for birch_clustering: ',t2-t1, ' [s]'
    write(ilog,*) 'BIRCH DONE'
    call gen_MST_from_treeclustering(adjl_deg,adjl_ix,adjl_dis,max_degr,trj_data)
    call CPU_time(t2)
    write(ilog,*) '>>TIME<< elapsed for SST building: ',t2-t1, ' [s]'
  end if


  ! Eventual file-dumping for mst (debugging)
  if(mst_print) then
    write(ilog,*) '-----------------------------------'
    write(ilog,*) 'DUMPING OF THE MST FOR DEBUGGING REASONS'
    inquire(file="mst_new.txt", exist=exist)
    if (exist) then
      mst_file_unit = freeunit()
      write(ilog,*) 'mst_new.txt already exists. It will be overwritten'
      open(mst_file_unit, file="mst_new.txt", status="replace", action="write")
    else
      open(mst_file_unit, file="mst_new.txt", status="new", action="write")
    end if
    do i=1,n_snaps
      write(format_var,*) adjl_deg(i)
      if(adjl_deg(i).gt.max_degr) write(ilog,*) 'During dumping the degree exceeded maximum possible for element',i
      write (format_var, "(A1,I"//adjustl(trim(format_var))//",A7)") "(", adjl_deg(i), "f15.22)"
      write(mst_file_unit,format_var)  adjl_dis(i,:)
    end do

    close(mst_file_unit)
    write(ilog,*) '...file-mst dumping done.'
  end if

  if(log_print) close(ilog)

end
