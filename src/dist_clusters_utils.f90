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
  use m_variables_gen
  implicit none

  real , INTENT(IN) :: sn2(n_xyz) ! Array of coo of a single snap (i)
  type(t_scluster), INTENT(IN) :: it ! reference cluster
  real, INTENT(OUT) :: val ! this is the distance from the clu (tmp_d)
  real vec1(n_xyz) ! n_xyz is the coordinates allocation size
  if(it%nmbrs.eq.0) then
    val = 0
  else

    if (tmp_dis_method.eq.5.or.tmp_dis_method.eq.11.or.tmp_dis_method.eq.1) then
      vec1 = it%sums(1:n_xyz)/(1.0 * it%nmbrs) ! cluster center
    end if
    call distance(val,vec1,sn2) ! val=tmp_d, vec1=centroid, sn=snapshot
  end if
end
!
!------------------------------------------------------------------------------------
!
subroutine cluster_to_cluster_d(val, it1, it2)
  use m_variables_gen
  use m_clustering
  implicit none

  type(t_scluster), INTENT(IN):: it1, it2
  real vec1(n_xyz), vec2(n_xyz)
  real, INTENT(OUT) :: val

  if (tmp_dis_method.eq.5.or.tmp_dis_method.eq.11.or.tmp_dis_method.eq.1) then
    vec1 = it1%sums(:)/(1.0 * it1%nmbrs) ! center it1
    vec2 = it2%sums(:)/(1.0 * it2%nmbrs) ! center it2
  end if
  call distance(val, vec1, vec2)
end
!
!------------------------------------------------------------------------------------------
!
subroutine distance(val2, veci, vecj)
  use m_variables_gen
! CALLED from:
! vec1 = it%sums(1:n_xyz,1)/(1.0*it%nmbrs)
! call distance(val,vec1,trj_data(1:n_xyz,i)) !val=tmp_d,it=scluster(j),i=i
! As i is the snapshot index could be that trj_data is taking the 1:n_xyz sized array from
! the big data-set trj_data            : data extracted for clustering (all read into memory -> large)


  implicit none

  real, INTENT(IN) :: veci(n_xyz), vecj(n_xyz)
  real vec_ref(n_xyz)
  real  hlp, hlp2
  real, INTENT(OUT) :: val2
  integer k

  vec_ref = 1
  val2 = 0.0
  hlp = 0.0
  hlp2 = 0.0

  if (tmp_dis_method.eq.1) then
    do k=1,n_xyz
      hlp = abs(veci(k) - vecj(k))
      if (hlp.gt.180.0) hlp = abs(hlp - 360.0)
      val2 = val2 + hlp*hlp
    end do
    val2 = sqrt(val2/(1.0*n_xyz))
  else if(tmp_dis_method.eq.5) then
    hlp = sum((veci(1:n_xyz) - vecj(1:n_xyz))**2)
    val2 = sqrt((3.0 * hlp)/(1.0 * (n_xyz))) ! why *3????? number of coo
  else if (tmp_dis_method.eq.11) then
    hlp = ACOS(dot_product(vec_ref,veci)/&
    (dot_product(vec_ref,vec_ref) + &
    dot_product(veci,veci)))
    hlp2 = ACOS(dot_product(vec_ref,vecj)/&
    (dot_product(vec_ref,vec_ref) + &
    dot_product(vecj,vecj)))
    val2 = hlp2 - hlp
    if(val2.gt.20.or.val2.lt.(-20)) write(ilog,*) "ATTENTION: an angle out of &
    &control = ", val2
  end if
end

!
!-------------------------------------------------------------------------------
!
subroutine cluster_addsnap_ix(it,i)
  use m_clustering
  use m_variables_gen
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
      it%alsz = 2
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
  use m_variables_gen
  implicit none

  type (t_scluster), INTENT(INOUT) :: it !cluster it is modified by adding i
  real, INTENT(IN) :: sn1(n_xyz) !snap that must be added (vector of positions)
  integer, intent(in) :: i !id of the snap
  ! for dist = 1 dihedral
  integer clstsz, k
  real normer, normer2, incr
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
      it%alsz = 2
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
      allocate(it%snaps(it%alsz))
    end if
    it%center = i
  else ! regular add
    it%nmbrs = it%nmbrs + 1 !number of elements in the cluster updating
    if (it%nmbrs.gt.it%alsz) then
      call cluster_resize(it)
    end if
  end if

  it%snaps(it%nmbrs) = i !add the snap
  clstsz = n_xyz
  !cludata = sn1
  if(tmp_dis_method.eq.1) then
    if (it%nmbrs.gt.1) then
      normer2 = 1.0/(1.0*it%nmbrs-1.0)
    else
      normer2 = 0.0
    end if
    normer = 1.0/(1.0*it%nmbrs)
  !   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (sn1(k)-normer2*it%sums(k).lt.-180.0) incr = 360.0
      if (sn1(k)-normer2*it%sums(k).gt.180.0) incr = -360.0
      incr = incr + sn1(k)
      it%sums(k) = it%sums(k) + incr
      it%sqsum = it%sqsum + incr*incr
      if (normer*it%sums(k).lt.-180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 + 720.0*it%sums(k)
        it%sums(k) = it%sums(k) + it%nmbrs*360.0
      end if
      if (normer*it%sums(k).gt.180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 - 720.0*it%sums(k)
        it%sums(k) = it%sums(k) - it%nmbrs*360.0
      end if
    end do
  else if((tmp_dis_method.eq.5.or.tmp_dis_method.eq.11)) then
    it%sums(1:n_xyz) = it%sums(1:n_xyz) + sn1 !
    it%sqsum = it%sqsum + dot_product(sn1,sn1)
  end if
end
!
!-------------------------------------------------------------------------------
!
subroutine cluster_resize(it)
  use m_clustering
  use m_variables_gen
  implicit none

  type(t_scluster) it
  integer tmpa(it%alsz), oldsz

  if (allocated(it%snaps).EQV..false.) then
    allocate(it%snaps(2))
    it%alsz = 2
    if (it%alsz.ge.it%nmbrs) return
  end if
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
  use m_variables_gen
!
  implicit none
!

  type (t_scluster) its, itl

  !for dihedral angles dist = 1
  real normer, normer2, normer3
  integer clstsz
  integer k
  real incr

  clstsz = n_xyz

  if (tmp_dis_method.eq.1) then
    normer2 = 1.0/(1.0*itl%nmbrs)
    normer3 = 1.0/(1.0*its%nmbrs)
    normer = 1.0/(1.0*itl%nmbrs+its%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    itl%sqsum = itl%sqsum + its%sqsum
    do k=1,clstsz
      incr = 0.0
      if (normer3*its%sums(k)-normer2*itl%sums(k).lt.-180.0) then
        incr = its%nmbrs*360.0
        itl%sqsum = itl%sqsum + its%nmbrs*360.0*360.0 + 720.0*its%sums(k)
      end if
      if (normer3*its%sums(k)-normer2*itl%sums(k).gt.180.0) then
        incr = -its%nmbrs*360.0
        itl%sqsum = itl%sqsum + its%nmbrs*360.0*360.0 - 720.0*its%sums(k)
      end if
      incr = incr + its%sums(k)
      itl%sums(k) = itl%sums(k) + incr
      if (normer*itl%sums(k).lt.-180.0) then
        itl%sqsum = itl%sqsum + (itl%nmbrs+its%nmbrs)*360.0*360.0 + 720.0*itl%sums(k)
        itl%sums(k) = itl%sums(k) + (itl%nmbrs+its%nmbrs)*360.0
      end if
      if (normer*itl%sums(k).gt.180.0) then
        itl%sqsum = itl%sqsum + (itl%nmbrs+its%nmbrs)*360.0*360.0 - 720.0*itl%sums(k)
        itl%sums(k) = itl%sums(k) - (itl%nmbrs+its%nmbrs)*360.0
      end if
    end do
  else if ((tmp_dis_method.eq.5).OR.(tmp_dis_method.eq.6)) then
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
subroutine preprocess_snapshot(trj3, i2, vecti2)
  use m_variables_gen
  implicit none

  real, intent(in) :: trj3(n_snaps,n_xyz)
  real, intent(inout) :: vecti2(n_xyz)
  integer,intent(in) :: i2
  if(tmp_dis_method.eq.5) then
    vecti2 = trj3(i2,1:n_xyz)
  else if(tmp_dis_method.eq.11) then
    if(i2.eq.1) then !to avoid error I can double the first value
       vecti2 = trj3(i2+1,1:n_xyz) - trj3(i2,1:n_xyz)
    else
      vecti2 = trj3(i2,1:n_xyz) - trj3(i2-1,1:n_xyz)
    end if
  end if
end
!
!------------------------------------------------------------------------------------------
!
subroutine cluster_addchild(it,itsidx,itschild,itschildidx)
!
  use m_clustering
!
  implicit none
!
  integer itschildidx,itsidx
  type(t_scluster) it,itschild
!
  it%nchildren = it%nchildren + 1
  if (it%childr_alsz.eq.0) then
    it%childr_alsz = 2
    if (allocated(it%children).EQV..true.) deallocate(it%children)
    allocate(it%children(it%childr_alsz))
  end if
  if (it%nchildren.gt.it%childr_alsz) then
    call scluster_resizech(it)
  end if
  it%children(it%nchildren) = itschildidx
  itschild%parent = itsidx
!
end
!
!-------------------------------------------------------------------------------
!
  subroutine scluster_resizech(it)
!
  use m_clustering
!
  implicit none
!
  type(t_scluster) it
  integer tmpa(it%childr_alsz),oldsz
!
  if (allocated(it%children).EQV..false.) then
    allocate(it%children(2))
    it%childr_alsz = 2
    if (it%childr_alsz.ge.it%nchildren) return
  end if
  oldsz = it%childr_alsz
  tmpa(:) = it%children(:)
  deallocate(it%children)
  it%childr_alsz = it%childr_alsz*4
  allocate(it%children(it%childr_alsz))
  it%children(1:oldsz) = tmpa(:)
!
end
!
!-------------------------------------------------------------------
!
! shorten list of basins based on allocation status of snaps
!
subroutine clusters_shorten(it,nba) !this should be done on the spot for memory reason
!
  use m_clustering
!
  implicit none
!
  integer i,nba,j
  type(t_scluster) it(nba)
!
  i = 1
  do while (i.le.nba)
    if (it(i)%nmbrs.eq.0) then !sort of fill of the list (it is a collapsing)
      do j=i,nba-1
        call copy_cluster(it(j+1),it(j))
      end do
      nba = nba - 1
    else
      i = i + 1
    end if
  end do
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine clusters_sort(it, nba)
!nba is the number of clusters
  use m_clustering
!
  implicit none
!
  integer nba,i,ii,jj
  integer, ALLOCATABLE:: iv(:,:)
    integer, ALLOCATABLE:: tmp_iv(:,:)
  type(t_scluster) it(nba),tmpit
!
  if (nba.le.1) return
!
! sort according to size
  allocate(iv(nba,3))
  allocate(tmp_iv((nba+1)/2,2))
  do i=1,nba
    iv(i,2) = i !(/ 1:nba /) !indexes
  end do
  iv(:,1) = it(1:nba)%nmbrs

  call MergeSort(iv(:,1),iv(:,2),nba,tmp_iv(:,1),tmp_iv(:,2),.false.)
  !call merge_sort(nba,upornot,iv(:,1),iv(:,2),ii,jj,iv(:,3),iv(:,4))
  ! 1->allnbs (3) !number of elements [nba]
  ! 3->Talldiss (4) !temporary balue for sorting iv(:,1)=> iv(:,2)
  ! 4->alldiss (1) WHAT will be sorted with [iv(:,2)] => iv(:,1)
  ! 6->allnbs (3)
  ! 7->ix (2) The indexes of WHAT will be sorted [iv(:,3)]
  ! 8->Tix (5) the temporary indexes [iv(:,4)]

  !as I believe that the order is based on the clusters' number of elements I will
  !sort that one
  ! allocate(alldiss(allnbs)) ! NOT
  ! allocate(Talldiss((allnbs+1)/2))
  ! allocate(ix(allnbs))!
  ! allocate(Tix((allnbs+1)/2))!
  do i=1,nba
    iv(iv(i,2),3) = i
  end do
!
! now use the maps - we need to move stuff around continuously to preserve all necessary information
  do i=1,nba
    if (i.eq.iv(i,2)) cycle
    call copy_cluster(it(i),tmpit)
    call copy_cluster(it(iv(i,2)),it(i))
    call copy_cluster(tmpit,it(iv(i,2)))
    ii = iv(i,3)
    jj = iv(i,2)
    iv(i,3) = iv(jj,3)
    iv(jj,3) = ii
    iv(ii,2) = iv(i,2)
  end do
!
  deallocate(iv)
!
end
!
!----------------------------------------------------------------------------------
!
subroutine cluster_getcenter(it,trj4)
!
  use m_clustering
  use m_variables_gen
!
  implicit none
!
  type(t_scluster) it
  real, intent(in) :: trj4(n_snaps,n_xyz)
  integer i
  real maxdis, tmp_di
  real vecti(n_xyz)
!
  maxdis = HUGE(maxdis)
  do i=1,it%nmbrs
    call preprocess_snapshot(trj4, it%snaps(i), vecti)
    call snap_to_cluster_d(tmp_di,it,vecti)
    if (tmp_di.lt.maxdis) then
      it%center = it%snaps(i)
      maxdis = tmp_di
    end if
  end do
!
end
!
!-------------------------------------------------------------------------------
!
subroutine copy_cluster(cfrom,cto)
!
  use m_clustering
  use m_variables_gen
!
  implicit none
!
  type(t_scluster) cfrom,cto
!
  if (allocated(cfrom%tmpsnaps).EQV..true.) then
    if (allocated(cto%tmpsnaps).EQV..true.) deallocate(cto%tmpsnaps)
    allocate(cto%tmpsnaps(cfrom%alsz))
    cto%tmpsnaps(1:cfrom%alsz) = cfrom%tmpsnaps(1:cfrom%alsz)
  end if
  if (allocated(cfrom%snaps).EQV..true.) then
    if (allocated(cto%snaps).EQV..true.) deallocate(cto%snaps)
    allocate(cto%snaps(cfrom%alsz))
    cto%snaps(1:cfrom%alsz) = cfrom%snaps(1:cfrom%alsz)
  end if
  if (allocated(cfrom%children).EQV..true.) then
    if (allocated(cto%children).EQV..true.) deallocate(cto%children)
    allocate(cto%children(cfrom%childr_alsz))
    cto%children(1:cfrom%childr_alsz) = cfrom%children(1:cfrom%childr_alsz)
  end if
  if (allocated(cfrom%sums).EQV..true.) then
    if (allocated(cto%sums).EQV..true.) deallocate(cto%sums)
    allocate(cto%sums(n_xyz))
    cto%sums(1:n_xyz) = cfrom%sums(1:n_xyz)
  end if
  cto%nmbrs = cfrom%nmbrs
  cto%sqsum = cfrom%sqsum
  cto%alsz = cfrom%alsz
! for tree-based clustering
  cto%nchildren = cfrom%nchildren
  cto%childr_alsz = cfrom%childr_alsz
  cto%parent = cfrom%parent
  !for quality
  cto%diam = cfrom%diam
  cto%radius = cfrom%radius
  cto%quality = cfrom%quality
  cto%center = cfrom%center
!
end
!
!--------------------------------------------------------------------------------------
!
subroutine cluster_calc_params(it,targetsz)
!
  use m_clustering
  use m_variables_gen
!
  implicit none
!
  type (t_scluster) it
  real helper, targetsz
  ! this is occasionally inaccurate for mode 1 due to the approximate treatment of wraparound
  if (tmp_dis_method.eq.1) then
    helper = (it%nmbrs*it%sqsum - dot_product(it%sums(:),it%sums(:)))/(1.0*n_xyz)
  !  this is approximate for all clusters exceeding size 2 if alignment is used
  else if(tmp_dis_method.eq.5.or.tmp_dis_method.eq.11) then
  helper = 3.0*(it%nmbrs*it%sqsum - dot_product(it%sums(1:n_xyz),it%sums(1:n_xyz)))&
  &              /(1.0*(n_xyz)) !sasdaafsjghdsviubwFB TODO
  end if

  if ((it%nmbrs.gt.1).AND.(helper.ge.0.0)) then
    it%radius = sqrt(helper/(1.0*it%nmbrs*it%nmbrs))
    it%diam = sqrt(2.0*helper/(1.0*it%nmbrs*(it%nmbrs-1.0)))
    it%quality = 1.0 - it%radius/targetsz
    ! write(ilog, *) "-----"
    ! write(ilog,*) helper,it%nmbrs
    ! write(ilog, * ) it%sqsum,it%sums(1)
    ! write(ilog,*)it%radius,it%diam
  else if (it%nmbrs.eq.1) then
    it%diam = 0.0
    it%quality = 1.0
    it%radius = 0.0
! because of FPE/vectorization, this number can occasionally be small and negative if the real variance is zero
  else if (abs(helper)/radius.le.1.0e-8) then
    write(ilog,*) "this happened"
    it%diam = 0.0
    it%quality = 1.0
    it%radius = 0.0
  else
    write(ilog,*) 'Warning. Encountered bad numbers in simplified computation of cluster &
 &parameters in cluster_calc_params(...).'
    it%diam = HUGE(it%diam)
    it%radius = HUGE(it%radius)
    it%quality = 0.0
  end if
!
end
!
!-------------------------------------------------------------------------------
!
subroutine quality_of_clustering(ncls,it,radcrit,quals)
!
  use m_clustering
  use m_variables_gen
!
  implicit none
!
  integer ncls,normer,alln,i
  type (t_scluster) it(ncls)
  real val1,val2,quals(3),unif,radcrit
!
  alln = sum(it(1:ncls)%nmbrs)
  unif = ncls*(1.0/(1.0*ncls))*log(1.0/(1.0*ncls))
  val1 = 0.0
  val2 = 0.0
  normer = 0
  do i=1,ncls
    if (it(i)%nmbrs.le.0) cycle
    if (it(i)%nmbrs.gt.1) normer = normer + it(i)%nmbrs
    val1 = val1 + it(i)%nmbrs*it(i)%radius
    val2 = val2 + it(i)%nmbrs/(1.0*alln)*log(it(i)%nmbrs/(1.0*alln))
  end do
!
 44 format('Quality of clustering at threshold of ',g12.5,':',/,'Tightness: ',g12.5,/,'Singles  : '&
 &,g12.5,/,'Entropy  : ',g12.5,/,'Total    : ',g12.5,/)
  if (normer.gt.0) then
    quals(1) = 1.0 - val1/(normer*radcrit)
  else
    quals(1) = 1.0
  end if
  quals(2) = 1.0 - (1.0*(alln-normer))/(1.0*alln)
  if (unif.ge.0.0) then
    quals(3) = 0.0
  else
    quals(3) = 1.0 - val2/unif
  end if
  write(ilog,44) radcrit,quals(:),sum(quals)/3.0
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine cluster_removesnap(it,i,vecti2)
!
  use m_clustering
  use m_variables_gen
!
  implicit none
!
  real, intent(in) :: vecti2(n_xyz)
  type (t_scluster) it
  integer i,k
  !for dist = 1
  integer clstsz
  real normer, normer2, incr

!
  do k=1,it%nmbrs
    if (it%snaps(k).eq.i) exit
  end do
  if (k.gt.it%nmbrs) then
    write(ilog,*) 'Fatal. Attempting to remove snap from cluster that is not a member of that cluster.&
 & This is a bug.'
    call fexit()
  end if
  it%snaps(k) = it%snaps(it%nmbrs) !is it removed like this?
  it%nmbrs = it%nmbrs - 1
!
! be safe (this condition should only occur in specialized applications)
  if (it%nmbrs.eq.0) then
    it%sums(:) = 0.0
    it%sqsum = 0.0
    if (it%alsz.gt.0) then
      it%alsz = 0
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
    end if
    return
  end if
!
  clstsz = n_xyz
  if (tmp_dis_method.eq.1) then
    if (it%nmbrs.gt.1) then
      normer2 = 1.0/(1.0*it%nmbrs+1.0)
    else
      normer2 = 0.0
    end if
    normer = 1.0/(1.0*it%nmbrs)
  !   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (vecti2(k)-normer2*it%sums(k).lt.-180.0) incr = 360.0
      if (vecti2(k)-normer2*it%sums(k).gt.180.0) incr = -360.0
      incr = -incr - vecti2(k)
      it%sums(k) = it%sums(k) + incr
      it%sqsum = it%sqsum + incr*incr
      if (normer*it%sums(k).lt.-180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 + 720.0*it%sums(k)
        it%sums(k) = it%sums(k) + it%nmbrs*360.0
      end if
      if (normer*it%sums(k).gt.180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 - 720.0*it%sums(k)
        it%sums(k) = it%sums(k) - it%nmbrs*360.0
      end if
    end do
  else if ((tmp_dis_method.eq.5).OR.(tmp_dis_method.eq.11)) then
    it%sums(1:n_xyz) = it%sums(1:n_xyz) - vecti2(1:n_xyz)
    it%sqsum = it%sqsum - dot_product(vecti2(1:n_xyz),vecti2(1:n_xyz))
  end if
!
end
!
!------------------------------------------------------------------------------------------
!
