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

! Notes on n_xyz: this variable could be found in the chainsaw program (main)
! There is called this function:
! ! sanity checks and setup for structural clustering
!   if (cstorecalc.le.nsim) then
!     call read_clusteringfile()
!   end if
! In this function the n_xyz variable is set to keep the coordinates bounds
! (e.g. xyz = 1 and 3)
!
!-----------------------------------------------------------------------------
!
subroutine snap_to_cluster_d(val, it, sn2) !i={1:n_snaps}
  use m_clustering
  use m_var_nbls_clu
  implicit none

  real , INTENT(IN) :: sn2(n_xyz) ! Array of coo of a single snap (i)
  type(t_scluster), INTENT(IN) :: it ! reference cluster
  real vec1(n_xyz), hlp ! n_xyz is the coordinates allocation size
  real, INTENT(OUT) :: val ! this is the distance from the clu (tmp_d)
  if(it%nmbrs.eq.0) then
    val = 0
  else

    if (dis_method.eq.5.or.dis_method.eq.11) then
      vec1 = it%sums(1:n_xyz)/(1.0 * it%nmbrs) ! cluster center
    end if
    call distance(val,vec1,sn2) ! val=tmp_d, vec1=centroid, sn=snapshot
  end if
end
!
!------------------------------------------------------------------------------------
!
subroutine cluster_to_cluster_d(val, it1, it2)
  use m_var_nbls_clu
  use m_clustering
  implicit none

  type(t_scluster), INTENT(IN):: it1, it2
  real vec1(n_xyz), hlp, vec2(n_xyz)
  integer k
  real, INTENT(OUT) :: val

  if (dis_method.eq.5.or.dis_method.eq.11) then
    vec1 = it1%sums(:)/(1.0 * it1%nmbrs) ! center it1
    vec2 = it2%sums(:)/(1.0 * it2%nmbrs) ! center it2
  end if
  call distance(val, vec1, vec2)
end
!
!------------------------------------------------------------------------------------------
!
subroutine distance(val, veci, vecj)
  use m_var_nbls_clu
! CALLED from:
! vec1 = it%sums(1:n_xyz,1)/(1.0*it%nmbrs)
! call distance(val,vec1,trj_data(1:n_xyz,i)) !val=tmp_d,it=scluster(j),i=i
! As i is the snapshot index could be that trj_data is taking the 1:n_xyz sized array from
! the big data-set trj_data            : data extracted for clustering (all read into memory -> large)


  implicit none

  real, INTENT(IN) :: veci(n_xyz), vecj(n_xyz)
  real vec_ref(n_xyz)
  real  hlp, hlp2
  real, INTENT(OUT) :: val

  vec_ref = 1
  val = 0.0
  hlp = 0.0
  hlp2 = 0.0

  if (dis_method.eq.5) then
    hlp = sum((veci(1:n_xyz) - vecj(1:n_xyz))**2)
    val = sqrt((3.0 * hlp)/(1.0 * (n_xyz))) ! why *3?????
  else if (dis_method.eq.11) then
    hlp = ACOS(dot_product(vec_ref,veci)/&
    (dot_product(vec_ref,vec_ref) + &
    dot_product(veci,veci)))
    hlp2 = ACOS(dot_product(vec_ref,vecj)/&
    (dot_product(vec_ref,vec_ref) + &
    dot_product(vecj,vecj)))
    val = hlp2 - hlp
  end if
end

!
!-------------------------------------------------------------------------------
!
subroutine cluster_addsnap_ix(it,i)
  use m_clustering
  use m_var_nbls_clu
  implicit none

  type (t_scluster), INTENT(INOUT) :: it !cluster it is modified by adding i
  integer, intent(in) :: i

  if (it%nmbrs.eq.0) then !case in which it is the first element of the cluster
    if (allocated(it%sums).EQV..false.) then
      allocate(it%sums(n_xyz))
      it%sums(:) = 0.0
    else
      it%sums(:) = 0.0
    end if
    it%sqsum = 0.0
    it%nmbrs = 1
    if (it%nmbrs.gt.it%alsz) then
      it%alsz = 1
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
      allocate(it%snaps(it%alsz))
    end if
  else ! regular add
    it%nmbrs = it%nmbrs + 1 !number of elements in the cluster updating
    if (it%nmbrs.gt.it%alsz) then
      call cluster_resize(it)
    end if
  end if

  it%snaps(it%nmbrs) = i !add the snap

end subroutine cluster_addsnap_ix
!-------------------------------------------------------------------------------
subroutine cluster_addsnap(it,sn1,i)
  use m_clustering
  use m_var_nbls_clu
  implicit none

  type (t_scluster), INTENT(INOUT) :: it !cluster it is modified by adding i
  real , INTENT(IN) :: sn1(n_xyz)
  integer, intent(in) :: i
  ! real , intent(in) :: r ! distance from the center of the cluster

  if (it%nmbrs.eq.0) then !case in which it is the first element of the cluster
    if (allocated(it%sums).EQV..false.) then
      allocate(it%sums(n_xyz))
      it%sums(:) = 0.0
    else
      it%sums(:) = 0.0
    end if
    it%sqsum = 0.0
    it%nmbrs = 1
    if (it%nmbrs.gt.it%alsz) then
      it%alsz = 1
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
      allocate(it%snaps(it%alsz))
    end if
    it%sums(1:n_xyz) = it%sums(1:n_xyz) + sn1 !
  else ! regular add
    it%nmbrs = it%nmbrs + 1 !number of elements in the cluster updating
    if (it%nmbrs.gt.it%alsz) then
      call cluster_resize(it)
    end if
  end if

  it%snaps(it%nmbrs) = i !add the snap

  if((dis_method.eq.5.or.dis_method.eq.11)) then
    it%sums(1:n_xyz) = it%sums(1:n_xyz) + sn1 !
    it%sqsum = it%sqsum + dot_product(sn1,sn1)

  end if
end
!
!-------------------------------------------------------------------------------
!
subroutine cluster_resize(it)
  use m_clustering
  use m_var_nbls_clu
  implicit none

  type(t_scluster) it
  integer tmpa(it%alsz), oldsz

  ! if (it%alsz.ge.it%nmbrs) return ! This should never happen because the
  ! function is called only in the opposite case

  oldsz = it%alsz
  tmpa(:) = it%snaps(:)
  deallocate(it%snaps)
  it%alsz = it%alsz * 2
  if(it%alsz.gt.n_snaps) it%alsz = n_snaps ! dont exagerate with space
  allocate(it%snaps(it%alsz)) ! reallocate snaps using allocation size of *2
  it%snaps(:) = 0 ! dont make fortran invent variables when not initialized
  it%snaps(1:oldsz) = tmpa(:) ! copy the old array into the new
end
!
!-----------------------------------------------------------------------------------------
!
subroutine join_clusters(itl,its)

  use m_clustering
  use m_var_nbls_clu
!
  implicit none
!

  type (t_scluster) its, itl
!
  if ((dis_method.eq.5).OR.(dis_method.eq.6)) then
    itl%sums(1:n_xyz) = itl%sums(1:n_xyz) + its%sums(1:n_xyz)
    itl%sqsum = itl%sqsum + its%sqsum
  end if
!
  itl%nmbrs = itl%nmbrs + its%nmbrs
  do while (itl%nmbrs.gt.itl%alsz)
    call cluster_resize(itl)
  end do
  itl%snaps(itl%nmbrs-its%nmbrs+1:itl%nmbrs) = its%snaps(1:its%nmbrs)
  its%nmbrs = 0
  its%alsz = 0
  deallocate(its%snaps)
  deallocate(its%sums)

end
!
!-----------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!
subroutine cnbl_resz(i2)
  use m_gen_nbls
  implicit none

  type(t_cnblst), INTENT(INOUT) :: i2 !it will be the i-element of the cnblst
  ! type(t_cnblst) cnblst
  integer oldalsz
  real distmp(i2%alsz)
  integer idxtmp(i2%alsz)
!
  distmp(:) = i2%dis(:)
  idxtmp(:) = i2%idx(:)
  deallocate(i2%dis)
  deallocate(i2%idx)
  oldalsz = i2%alsz
  i2%alsz = i2%alsz*2
  allocate(i2%dis(i2%alsz))
  allocate(i2%idx(i2%alsz))
  i2%dis(1:oldalsz) = distmp(:)
  i2%idx(1:oldalsz) = idxtmp(:)
!
end
!
!-------------------------------------------------------------------------------
!
