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
! trj_data = cludata
! tmp_d = rmsdtmp
!
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
  dis_method_in, dis_weight_in, mst_log_in, data_meth_in, normalize_dis_in, verbose_in) !modes

  use m_var_nbls_clu
  use m_clustering
  use m_gen_nbls
  use m_mst
  implicit none

  integer, INTENT(IN) :: n_xyz_in !numbers of xyz (atoms*3)
  integer, INTENT(IN) :: n_snaps_in !number of snapshots in input
  real(KIND=4), INTENT(IN) :: trj_data(n_snaps_in, n_xyz_in) !trajectory input
  integer, intent(in) :: dis_method_in(11) !distance method
  real(KIND=4), intent(in) :: dis_weight_in(11) !distance weight when averaged
  real(KIND=4), intent(in) :: clu_radius_in !defines the max cluster sizes
  real(KIND=4), intent(in) :: clu_hardcut_in !threshold between cluster snaps
  logical, intent(in) :: verbose_in !verbose terminal output
  integer, intent(in) :: data_meth_in !data managing method TODO netcdf
  logical, intent(in) :: mst_log_in !make already the ordered mst(true) or not?
  logical, intent(in) :: normalize_dis_in !flag for normalize the distance matrix
  !if the intent is in this cannot be an ALLOCATABLE variable


  ! N B :
  ! clu_hardcut is used for distances between different clusters snapshos as
  ! a threshold.

  integer i,j,ii,kk,ll,k,l,mi,mj,u,u2 ! helper variables
  integer nzeros ! nzeros = number of not connected components (dis .le. 0)
  logical exist, mst_print !file dumping for mst tree for debugging
  character(len=1024) :: format_var !format for above mentioned dumping

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
  superver = .false.  !dev flag for superverbose output
  mst_print = .false. !dev flag for adjlist (mst) checking
  radius = clu_radius_in
  hardcut = clu_hardcut_in
  dis_method = dis_method_in
  dis_weight = 1
  mst_log = mst_log_in
  normalize_dis = normalize_dis_in


  n_dis_method = 0
  do i=1,11
    if(dis_method(i).ge.1.and.dis_method(i).le.11) n_dis_method = n_dis_method + 1
    if(dis_weight_in(i).ge.0.and.dis_weight_in(i).le.1) dis_weight(i) = dis_weight_in(i)
  end do
  write(*,*)
  if(n_dis_method.gt.1) write(*,*) 'Multiple distance functions selected. They will be &
  & balanced using clustering valuse, then summed toghether in a old fashion way.\n&
  &Remember that the leader clustering algorithm will consider only the first value &
  &of the distances'
  write(*,*)
  write(*,*) "Selected distances: ", dis_method(1:n_dis_method)
  write(*,*) "Distance weights:", dis_weight(1:n_dis_method)
  write(*,*)

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

  write(*,*) "Input dimensions:", n_snaps, " row and ", n_xyz," col"
  write(*,*)
  if(superver) write(*,*) "clu_radius_in", clu_radius_in
  if(superver) write(*,*) "clu_hardcut_in", clu_hardcut_in
  write(*,*)
  if(superver) write(*,*) "Input example", trj_data(1:10,1:10)
  write(*,*)
  write(*,*)
  write(*,*) 'Now using truncated leader algorithm for pre-clustering&
   &in neighbor list generation ...'

  nclalcsz = 2 ! a priori cluster numbers
  allocate(scluster(nclalcsz))
  scluster(:)%nmbrs = 0 ! number of elements in each cluster
  scluster(:)%alsz = 0 ! allocation size of each cluster

  call leader_clustering(trj_data)

  ! now compare all blocks to each other (the slowest part) taking advantage
  ! of information generated previously (otherwise intractable)
  write(*,*)
  write(*,*) 'Now computing cutoff-assisted neighbor list...'
  write(*,*)

  call gen_nb(trj_data)

  write(*,*)
  write(*,*) 'Neighbor list generated.'
  write(*,*)

  nzeros = 0
  do i=1,n_snaps
    if (cnblst(i)%nbs.le.0) then
      nzeros = nzeros + 1
     write(*,*) 'Warning. Snapshot # ',i,' is without a neighbor (similar) &
     &structure. This may cause the clustering algorithm to crash or &
     &misbehave otherwise.'
   end if
  end do
  if (nzeros.gt.0) then
    write(*,*) 'Warning. ',nzeros,' snapshots are without a neighbor &
    &(similar) structure. This may in some cases cause the clustering &
    &algorithm to misbehave.'
  end if
  write(*,*)

  do i=1,nclalcsz
    if (allocated(scluster(i)%snaps).EQV..true.) deallocate(scluster(i)%snaps)
    if (allocated(scluster(i)%sums).EQV..true.) deallocate(scluster(i)%sums)
  end do
  deallocate(scluster)

  if(mst_log) then
    call gen_MST_from_nbl(adjl_deg,adjl_ix,adjl_dis,max_degr)
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

  ! Eventual file-dumping for mst (debugging)
  if(mst_print) then
    inquire(file="mst_original.txt", exist=exist)
    if (exist) then
      write(*,*) 'mst_original.txt already exists. It will be overwritten'
      open(378, file="mst_original.txt", status="replace", action="write")
    else
      open(378, file="mst_original.txt", status="new", action="write")
    end if
    write(*,*) 'the problem is not here'
    do i=1,n_snaps
      write (format_var, "(A1,I1,A7)") "(", adjl_deg(i), "f15.10)"
      write(378, format_var)  adjl_dis(i,:)
    end do
    close(378)
    write(*,*) '...file-mst dumping done.'
  end if


end
