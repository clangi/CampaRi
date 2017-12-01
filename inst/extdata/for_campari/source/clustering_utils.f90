!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 3.0                                                           !
!                                                                          !
!    Copyright (C) 2017, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  !
!                        Davide Garolini, Jiri Vymetal                     !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    CAMPARI is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-------------------------------------------------------------------------
!
subroutine snap_to_snap_d(val,i,j)
!
  use clusters
!
  implicit none
!
  integer i,j
  RTYPE val
!
  call clustering_distance(val,cludata(1:calcsz,i),cludata(1:calcsz,j))
!
end
!
!------------------------------------------------------------------------------------
!
subroutine snap_to_cluster_d(val,it,i)
!
  use clusters
  use atoms
!
  implicit none
!
  integer, INTENT(IN):: i
  type(t_scluster), INTENT(IN):: it
!
  integer k
  RTYPE val,vec1(calcsz),hlp
!
  if (cdis_crit.eq.2) then
    do k=1,clstsz
      vec1(2*k-1) = it%sums(k,1)/(1.0*it%nmbrs)
      vec1(2*k) = it%sums(k,5)/(1.0*it%nmbrs)
    end do
  else if (cdis_crit.eq.4) then
    do k=1,calcsz/3
      vec1(3*k-2) = it%sums(2*k-1,1)/(1.0*it%nmbrs)
      vec1(3*k-1) = it%sums(2*k,1)/(1.0*it%nmbrs)
      vec1(3*k) = it%sums(2*k,5)/(1.0*it%nmbrs)
    end do
  else if (cdis_crit.eq.8) then
    hlp = sum((it%sums(1:calcsz,1)/(1.0*it%nmbrs) - cl_imvec(1:calcsz)*cludata(1:calcsz,i))**2)
    val = sqrt(hlp)
    return
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    vec1(1:clstsz) = it%sums(:,1)/(1.0*it%nmbrs)
    vec1((clstsz+1):calcsz) = it%sums(:,5)/(1.0*it%nmbrs)
  else
    vec1 = it%sums(1:calcsz,1)/(1.0*it%nmbrs)
  end if
  call clustering_distance(val,vec1,cludata(1:calcsz,i))
!
end
!
!------------------------------------------------------------------------------------
!
subroutine cluster_to_cluster_d(val,it1,it2)
!
  use clusters
  use atoms
!
  implicit none
!
  integer k
  type(t_scluster) it1,it2
  RTYPE val,vec1(calcsz),hlp,vec2(calcsz)
!
  if (cdis_crit.eq.2) then
    do k=1,clstsz
      vec1(2*k-1) = it1%sums(k,1)/(1.0*it1%nmbrs)
      vec1(2*k) = it1%sums(k,5)/(1.0*it1%nmbrs)
      vec2(2*k-1) = it2%sums(k,1)/(1.0*it2%nmbrs)
      vec2(2*k) = it2%sums(k,5)/(1.0*it2%nmbrs)
    end do
  else if (cdis_crit.eq.4) then
    do k=1,calcsz/3
      vec1(3*k-2) = it1%sums(2*k-1,1)/(1.0*it1%nmbrs)
      vec1(3*k-1) = it1%sums(2*k,1)/(1.0*it1%nmbrs)
      vec1(3*k) = it1%sums(2*k,5)/(1.0*it1%nmbrs)
      vec2(3*k-2) = it2%sums(2*k-1,1)/(1.0*it2%nmbrs)
      vec2(3*k-1) = it2%sums(2*k,1)/(1.0*it2%nmbrs)
      vec2(3*k) = it2%sums(2*k,5)/(1.0*it2%nmbrs)
    end do
  else if (cdis_crit.eq.8) then
    hlp = sum((it1%sums(:,1)/(1.0*it1%nmbrs) - it2%sums(:,1)/(1.0*it2%nmbrs))**2)
    val = sqrt(hlp)
    return
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    vec1(1:clstsz) = it1%sums(:,1)/(1.0*it1%nmbrs)
    vec1((clstsz+1):calcsz) = it1%sums(:,5)/(1.0*it1%nmbrs)
    vec2(1:clstsz) = it2%sums(:,1)/(1.0*it2%nmbrs)
    vec2((clstsz+1):calcsz) = it2%sums(:,5)/(1.0*it2%nmbrs)
  else
    vec1 = it1%sums(:,1)/(1.0*it1%nmbrs)
    vec2 = it2%sums(:,1)/(1.0*it2%nmbrs)
  end if
  call clustering_distance(val,vec1,vec2)
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine clustering_distance(val,veci,vecj)
!
  use clusters
  use atoms
!
  implicit none
!
  RTYPE, INTENT(IN):: veci(calcsz),vecj(calcsz)
!
  integer k,catnrs
  RTYPE vec1(calcsz),vec2(calcsz),vec3(calcsz),val,tvec(3),qrot(4),cto(3),ctn(3),vectmp(3)
  RTYPE hlp,hlp2,whlp,genmu
!
  val = 0.0
  hlp2 = 0.0
  hlp = 0.0
!
  if (cdis_crit.eq.1) then
    do k=1,calcsz
      hlp = abs(veci(k) - vecj(k))
      if (hlp.gt.180.0) hlp = abs(hlp - 360.0)
      val = val + hlp*hlp
    end do
    val = sqrt(val/(1.0*calcsz))
  else if (cdis_crit.eq.2) then
    do k=1,calcsz/2
      hlp = abs(veci(2*k-1) - vecj(2*k-1))
      if (hlp.gt.180.0) hlp = abs(hlp - 360.0)
      whlp = genmu(veci(2*k),vecj(2*k),cwcombination)
      val = val + whlp*hlp*hlp
      hlp2 = hlp2 + whlp
    end do
    val = sqrt(val/hlp2)
  else if (cdis_crit.eq.3) then
    hlp = sum((veci(:) - vecj(:))**2)
    val = sqrt(hlp/(1.0*calcsz))
  else if (cdis_crit.eq.4) then
    do k=1,calcsz/3
      hlp = (veci(3*k-2) - vecj(3*k-2))**2 + (veci(3*k-1) - vecj(3*k-1))**2
      whlp = genmu(veci(3*k),vecj(3*k),cwcombination)
      val = val + whlp*hlp
      hlp2 = hlp2 + 2.0*whlp
    end do
    val = sqrt(val/hlp2)
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    if (align_for_clustering.EQV..true.) then
      catnrs = (cdofsbnds(4)-cdofsbnds(3)+1)/3
      vec1(:) = vecj(:)
      vec2(:) = veci(:)
      call align_3D(catnrs,vec1(cdofsbnds(3):cdofsbnds(4)),vec2(cdofsbnds(3):cdofsbnds(4)),tvec,qrot,cto,ctn)
      do k=1,calcsz/3
        call quat_conjugatei(qrot,vec1(3*k-2:3*k),vectmp,cto)
        vec3(3*k-2:3*k) = vectmp(:)
        vec3(3*k-2:3*k) = vec3(3*k-2:3*k) + tvec(:)
      end do
      hlp = sum((vec2(cdofsbnds(1):cdofsbnds(2)) - vec3(cdofsbnds(1):cdofsbnds(2)))**2)
      val = sqrt(3.0*hlp/(1.0*(cdofsbnds(2)-cdofsbnds(1)+1)))
    else
      hlp = sum((veci(cdofsbnds(1):cdofsbnds(2)) - vecj(cdofsbnds(1):cdofsbnds(2)))**2)
      val = sqrt(3.0*hlp/(1.0*(cdofsbnds(2)-cdofsbnds(1)+1)))
    end if
  else if (cdis_crit.eq.7) then
    hlp = sum((veci(:) - vecj(:))**2)
    val = sqrt(hlp/(1.0*calcsz))
  else if (cdis_crit.eq.8) then
    vec1(:) = vecj(:)*cl_imvec(1:calcsz)
    vec2(:) = veci(:)*cl_imvec(1:calcsz)
    hlp =  sum((vec2(:) - vec1(:))**2)
    val = sqrt(hlp)
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    do k=clstsz+1,calcsz
      vec3(k-clstsz) = genmu(veci(k),vecj(k),cwcombination)
    end do
    hlp2 = 1.0/sum(vec3(1:clstsz))
    if (cdis_crit.eq.10) hlp2 = 3.0*hlp2
    vec1(1:clstsz) = veci(1:clstsz)-vecj(1:clstsz)
    hlp = sum(vec3(1:clstsz)*vec1(1:clstsz)*vec1(1:clstsz))
    val = sqrt(hlp*hlp2)
  end if
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine cnbl_resz(i)
!
  use clusters
!
  implicit none
!
  integer i,oldalsz
  RTYPE distmp(cnblst(i)%alsz)
  integer idxtmp(cnblst(i)%alsz)
!
  distmp(:) = cnblst(i)%dis(:)
  idxtmp(:) = cnblst(i)%idx(:)
  deallocate(cnblst(i)%dis)
  deallocate(cnblst(i)%idx)
  oldalsz = cnblst(i)%alsz
  cnblst(i)%alsz = cnblst(i)%alsz*2
  allocate(cnblst(i)%dis(cnblst(i)%alsz))
  allocate(cnblst(i)%idx(cnblst(i)%alsz))
  cnblst(i)%dis(1:oldalsz) = distmp(:)
  cnblst(i)%idx(1:oldalsz) = idxtmp(:)
!
end
!
!-----------------------------------------------------------------------------------------------
!
subroutine reduce_set(effsetsz,effsetis,tryat)
!
  use clusters
!
  implicit none
!
  integer i,ii,k
  integer effsetsz,effsetis(cstored),tryat(cstored),bla(effsetsz)
!
  bla(:) = effsetis(1:effsetsz)
!
  do ii=1,effsetsz
    i = bla(ii)
    if (tryat(i).gt.cnblst(i)%nbs) then
      bla(ii) = 0
    end if
  end do
!
  k = 0
  do ii=1,effsetsz
    if (bla(ii).gt.0) then
      k = k + 1
      effsetis(k) = bla(ii)
    end if
  end do
  effsetsz = k
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine cluster_addsnap(it,i,rdv)
!
  use clusters
  use math
  use atoms
!
  implicit none
!
  type (t_scluster) it
  integer i,k,catnrs
  RTYPE incr,vec1(calcsz),vec2(calcsz),vectmp(calcsz),qrot(4),cto(3),ctn(3),tvec(3),normer,normer2,rdv
!
  if (it%nmbrs.eq.0) then
    if (allocated(it%sums).EQV..false.) then
      allocate(it%sums(clstsz,sumssz))
      it%sums(:,:) = 0.0
    else
      it%sums(:,:) = 0.0
    end if
    it%sqsum = 0.0
    it%nmbrs = 1
    if (it%nmbrs.gt.it%alsz) then
      it%alsz = 2
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
      allocate(it%snaps(it%alsz))
    end if
    it%center = i 
  else
!   regular add
    it%nmbrs = it%nmbrs + 1
    if (it%nmbrs.gt.it%alsz) then
      call scluster_resize(it)
    end if
  end if
  it%snaps(it%nmbrs) = i
  if (cdis_crit.eq.1) then
    if (it%nmbrs.gt.1) then
      normer2 = 1.0/(1.0*it%nmbrs-1.0)
    else
      normer2 = 0.0
    end if
    normer = 1.0/(1.0*it%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (cludata(k,i)-normer2*it%sums(k,1).lt.-180.0) incr = 360.0
      if (cludata(k,i)-normer2*it%sums(k,1).gt.180.0) incr = -360.0
      incr = incr + cludata(k,i)
      it%sums(k,1) = it%sums(k,1) + incr
      it%sqsum = it%sqsum + incr*incr
      if (normer*it%sums(k,1).lt.-180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 + 720.0*it%sums(k,1)
        it%sums(k,1) = it%sums(k,1) + it%nmbrs*360.0
      end if
      if (normer*it%sums(k,1).gt.180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 - 720.0*it%sums(k,1)
        it%sums(k,1) = it%sums(k,1) - it%nmbrs*360.0
      end if
    end do
  else if (cdis_crit.eq.2) then
    if (it%nmbrs.gt.1) then
      normer2 = 1.0/(1.0*it%nmbrs-1.0)
    else
      normer2 = 0.0
    end if
    normer = 1.0/(1.0*it%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (cludata(2*k-1,i)-normer2*it%sums(k,1).lt.-180.0) incr = 360.0
      if (cludata(2*k-1,i)-normer2*it%sums(k,1).gt.180.0) incr = -360.0
      incr = incr + cludata(2*k-1,i)
      it%sums(k,1) = it%sums(k,1) + incr
      it%sums(k,2) = it%sums(k,2) + incr*incr
      it%sums(k,3) = it%sums(k,3) + incr*cludata(2*k,i)
      it%sums(k,4) = it%sums(k,4) + incr*incr*cludata(2*k,i)
      it%sums(k,5) = it%sums(k,5) + cludata(2*k,i)
      if (normer*it%sums(k,1).lt.-180.0) then
        it%sums(k,2) = it%sums(k,2) + it%nmbrs*360.0*360.0 + 720.0*it%sums(k,1)
        it%sums(k,4) = it%sums(k,4) + it%sums(k,5)*360.0*360.0 + 720.0*it%sums(k,3)
        it%sums(k,1) = it%sums(k,1) + it%nmbrs*360.0
        it%sums(k,3) = it%sums(k,3) + it%sums(k,5)*360.0
      end if
      if (normer*it%sums(k,1).gt.180.0) then
        it%sums(k,2) = it%sums(k,2) + it%nmbrs*360.0*360.0 - 720.0*it%sums(k,1)
        it%sums(k,4) = it%sums(k,4) + it%sums(k,5)*360.0*360.0 - 720.0*it%sums(k,3)
        it%sums(k,1) = it%sums(k,1) - it%nmbrs*360.0
        it%sums(k,3) = it%sums(k,3) - it%sums(k,5)*360.0
      end if
    end do
  else if ((cdis_crit.eq.3).OR.(cdis_crit.eq.7)) then
    it%sums(1:calcsz,1) = it%sums(1:calcsz,1) + cludata(1:calcsz,i)
    it%sqsum = it%sqsum + dot_product(cludata(1:calcsz,i),cludata(1:calcsz,i))
  else if (cdis_crit.eq.4) then
    do k=1,calcsz/3
      it%sums(2*k-1:2*k,1) = it%sums(2*k-1:2*k,1) + cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,2) = it%sums(2*k-1:2*k,2) + cludata(3*k-2:3*k-1,i)*cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,3) = it%sums(2*k-1:2*k,3) + cludata(3*k,i)*cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,4) = it%sums(2*k-1:2*k,4) + cludata(3*k,i)*cludata(3*k-2:3*k-1,i)*cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,5) = it%sums(2*k-1:2*k,5) + cludata(3*k,i)
    end do
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    tvec(:) = 0.0
    if (align_for_clustering.EQV..true.) then
      catnrs = (cdofsbnds(4)-cdofsbnds(3)+1)/3
      vec1(:) = cludata(1:calcsz,i)
      vec2(:) = it%sums(1:calcsz,1)/(1.0*it%nmbrs)
      call align_3D(catnrs,vec1(cdofsbnds(3):cdofsbnds(4)),vec2(cdofsbnds(3):cdofsbnds(4)),tvec,qrot,cto,ctn)
      do k=1,calcsz/3
        call quat_conjugatei(qrot,vec1(3*k-2:3*k),vectmp(1:3),cto)
        vec2(3*k-2:3*k) = vectmp(1:3)
        vec2(3*k-2:3*k) = vec2(3*k-2:3*k) + tvec(:)
      end do
      it%sums(1:calcsz,1) = it%sums(1:calcsz,1) + vec2(:)
      it%sqsum = it%sqsum + dot_product(vec2(cdofsbnds(1):cdofsbnds(2)),vec2(cdofsbnds(1):cdofsbnds(2)))
    else
      it%sums(1:calcsz,1) = it%sums(1:calcsz,1) + cludata(1:calcsz,i)
      it%sqsum = it%sqsum + dot_product(cludata(cdofsbnds(1):cdofsbnds(2),i),cludata(cdofsbnds(1):cdofsbnds(2),i))
    end if
  else if (cdis_crit.eq.8) then
    vectmp(:) = cl_imvec(1:calcsz)*cludata(1:calcsz,i)
    it%sums(1:calcsz,1) = it%sums(1:calcsz,1) + vectmp(:)
    it%sqsum = it%sqsum + dot_product(vectmp(:),vectmp(:))
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    it%sums(:,1) = it%sums(:,1) + cludata(1:clstsz,i)
    it%sums(:,2) = it%sums(:,2) + cludata(1:clstsz,i)*cludata(1:clstsz,i)
    it%sums(:,3) = it%sums(:,3) + cludata((clstsz+1):calcsz,i)*cludata(1:clstsz,i)
    it%sums(:,4) = it%sums(:,4) + cludata((clstsz+1):calcsz,i)*cludata(1:clstsz,i)*cludata(1:clstsz,i)
    it%sums(:,5) = it%sums(:,5) + cludata((clstsz+1):calcsz,i)
  end if
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine cluster_addjustsnap(it,i)
!
  use clusters
!
  implicit none
!
  type (t_scluster) it
  integer i
!
  if (it%nmbrs.eq.0) then
    it%nmbrs = 1
    if (it%nmbrs.gt.it%alsz) then
      it%alsz = 2
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
      allocate(it%snaps(it%alsz))
    end if
    it%center = i 
  else
!   regular add
    it%nmbrs = it%nmbrs + 1
    if (it%nmbrs.gt.it%alsz) then
      call scluster_resize(it)
    end if
  end if
  it%snaps(it%nmbrs) = i
!
  end
!
!-------------------------------------------------------------------------------------------
!
subroutine cluster_removesnap(it,i)
!
  use clusters
  use math
  use atoms
  use iounit
!
  implicit none
!
  type (t_scluster) it
  integer i,k,catnrs
  RTYPE incr,vec1(calcsz),vec2(calcsz),vectmp(calcsz),qrot(4),cto(3),ctn(3),tvec(3),normer,normer2
!
  do k=1,it%nmbrs
    if (it%snaps(k).eq.i) exit
  end do
  if (k.gt.it%nmbrs) then
    write(ilog,*) 'Fatal. Attempting to remove snap from cluster that is not a member of that cluster.&
 & This is a bug.'
    call fexit()
  end if
  it%snaps(k) = it%snaps(it%nmbrs)
  it%nmbrs = it%nmbrs - 1
!
! be safe (this condition should only occur in specialized applications)
  if (it%nmbrs.eq.0) then
    it%sums(:,:) = 0.0
    it%sqsum = 0.0
    if (it%alsz.gt.0) then
      it%alsz = 0
      if (allocated(it%snaps).EQV..true.) deallocate(it%snaps)
    end if
    return
  end if
!
  if (cdis_crit.eq.1) then
    if (it%nmbrs.gt.1) then
      normer2 = 1.0/(1.0*it%nmbrs+1.0)
    else
      normer2 = 0.0
    end if
    normer = 1.0/(1.0*it%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (cludata(k,i)-normer2*it%sums(k,1).lt.-180.0) incr = 360.0
      if (cludata(k,i)-normer2*it%sums(k,1).gt.180.0) incr = -360.0
      incr = -incr - cludata(k,i)
      it%sums(k,1) = it%sums(k,1) + incr
      it%sqsum = it%sqsum + incr*incr
      if (normer*it%sums(k,1).lt.-180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 + 720.0*it%sums(k,1)
        it%sums(k,1) = it%sums(k,1) + it%nmbrs*360.0
      end if
      if (normer*it%sums(k,1).gt.180.0) then
        it%sqsum = it%sqsum + it%nmbrs*360.0*360.0 - 720.0*it%sums(k,1)
        it%sums(k,1) = it%sums(k,1) - it%nmbrs*360.0
      end if
    end do
  else if (cdis_crit.eq.2) then
    if (it%nmbrs.gt.1) then
      normer2 = 1.0/(1.0*it%nmbrs+1.0)
    else
      normer2 = 0.0
    end if
    normer = 1.0/(1.0*it%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (cludata(2*k-1,i)-normer2*it%sums(k,1).lt.-180.0) incr = 360.0
      if (cludata(2*k-1,i)-normer2*it%sums(k,1).gt.180.0) incr = -360.0
      incr = -incr - cludata(2*k-1,i)
      it%sums(k,1) = it%sums(k,1) + incr
      it%sums(k,2) = it%sums(k,2) - incr*incr
      it%sums(k,3) = it%sums(k,3) + incr*cludata(2*k,i)
      it%sums(k,4) = it%sums(k,4) - incr*incr*cludata(2*k,i)
      it%sums(k,5) = it%sums(k,5) - cludata(2*k,i)
      if (normer*it%sums(k,1).lt.-180.0) then
        it%sums(k,2) = it%sums(k,2) + it%nmbrs*360.0*360.0 + 720.0*it%sums(k,1)
        it%sums(k,4) = it%sums(k,4) + it%sums(k,5)*360.0*360.0 + 720.0*it%sums(k,3)
        it%sums(k,1) = it%sums(k,1) + it%nmbrs*360.0
        it%sums(k,3) = it%sums(k,3) + it%sums(k,5)*360.0
      end if
      if (normer*it%sums(k,1).gt.180.0) then
        it%sums(k,2) = it%sums(k,2) + it%nmbrs*360.0*360.0 - 720.0*it%sums(k,1)
        it%sums(k,4) = it%sums(k,4) + it%sums(k,5)*360.0*360.0 - 720.0*it%sums(k,3)
        it%sums(k,1) = it%sums(k,1) - it%nmbrs*360.0
        it%sums(k,3) = it%sums(k,3) - it%sums(k,5)*360.0
      end if
    end do
  else if ((cdis_crit.eq.3).OR.(cdis_crit.eq.7)) then
    it%sums(1:calcsz,1) = it%sums(1:calcsz,1) - cludata(1:calcsz,i)
    it%sqsum = it%sqsum - dot_product(cludata(1:calcsz,i),cludata(1:calcsz,i))
  else if (cdis_crit.eq.4) then
    do k=1,calcsz/3
      it%sums(2*k-1:2*k,1) = it%sums(2*k-1:2*k,1) - cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,2) = it%sums(2*k-1:2*k,2) - cludata(3*k-2:3*k-1,i)*cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,3) = it%sums(2*k-1:2*k,3) - cludata(3*k,i)*cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,4) = it%sums(2*k-1:2*k,4) - cludata(3*k,i)*cludata(3*k-2:3*k-1,i)*cludata(3*k-2:3*k-1,i)
      it%sums(2*k-1:2*k,5) = it%sums(2*k-1:2*k,5) - cludata(3*k,i)
    end do
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    if (align_for_clustering.EQV..true.) then
      catnrs = (cdofsbnds(4)-cdofsbnds(3)+1)/3
      vec1(:) = cludata(1:calcsz,i)
      vec2(:) = it%sums(1:calcsz,1)/(1.0*it%nmbrs)
      call align_3D(catnrs,vec1(cdofsbnds(3):cdofsbnds(4)),vec2(cdofsbnds(3):cdofsbnds(4)),tvec,qrot,cto,ctn)
      do k=1,calcsz/3
        call quat_conjugatei(qrot,vec1(3*k-2:3*k),vectmp(1:3),cto)
        vec2(3*k-2:3*k) = vectmp(1:3)
        vec2(3*k-2:3*k) = vec2(3*k-2:3*k) + tvec(:)
      end do
      it%sums(1:calcsz,1) = it%sums(1:calcsz,1) - vec2(1:calcsz)
      it%sqsum = it%sqsum - dot_product(vec2(cdofsbnds(1):cdofsbnds(2)),vec2(cdofsbnds(1):cdofsbnds(2)))
    else
      it%sums(1:calcsz,1) = it%sums(1:calcsz,1) - cludata(1:calcsz,i)
      it%sqsum = it%sqsum - dot_product(cludata(cdofsbnds(1):cdofsbnds(2),i),cludata(cdofsbnds(1):cdofsbnds(2),i))
    end if
  else if (cdis_crit.eq.8) then
    vectmp(:) = cl_imvec(1:calcsz)*cludata(1:calcsz,i)
    it%sums(1:calcsz,1) = it%sums(1:calcsz,1) - vectmp(:)
    it%sqsum = it%sqsum - dot_product(vectmp(:),vectmp(:))
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    it%sums(:,1) = it%sums(:,1) - cludata(1:clstsz,i)
    it%sums(:,2) = it%sums(:,2) - cludata(1:clstsz,i)*cludata(1:clstsz,i)
    it%sums(:,3) = it%sums(:,3) - cludata((clstsz+1):calcsz,i)*cludata(1:clstsz,i)
    it%sums(:,4) = it%sums(:,4) - cludata((clstsz+1):calcsz,i)*cludata(1:clstsz,i)*cludata(1:clstsz,i)
    it%sums(:,5) = it%sums(:,5) - cludata((clstsz+1):calcsz,i)
  end if
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine cluster_addchild(it,itsidx,itschild,itschildidx)
!
  use clusters
!
  implicit none
!
  integer itschildidx,itsidx
  type(t_scluster) it,itschild
!
  it%nchildren = it%nchildren + 1
  if (it%chalsz.eq.0) then
    it%chalsz = 2
    if (allocated(it%children).EQV..true.) deallocate(it%children)
    allocate(it%children(it%chalsz))
  end if
  if (it%nchildren.gt.it%chalsz) then
    call scluster_resizech(it)
  end if
  it%children(it%nchildren) = itschildidx
  itschild%parent = itsidx
!
end
!
!--------------------------------------------------------------------------------------
!
subroutine cluster_calc_params(it,targetsz)
!
  use clusters
  use atoms 
  use iounit
!
  implicit none
!
  type (t_scluster) it
  RTYPE helper,targetsz
!
! if we do not have data to calculate cluster properties, default
  if (cmode.eq.7) then
    helper = 0.
! this is occasionally inaccurate for mode 1 due to the approximate treatment of wraparound
  else if ((cdis_crit.eq.7).OR.(cdis_crit.eq.3).OR.(cdis_crit.eq.1)) then
    helper = (it%nmbrs*it%sqsum - dot_product(it%sums(:,1),it%sums(:,1)))/(1.0*calcsz)
!  this is approximate for all clusters exceeding size 2 if alignment is used
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    helper = 3.0*(it%nmbrs*it%sqsum - dot_product(it%sums(cdofsbnds(1):cdofsbnds(2),1),it%sums(cdofsbnds(1):cdofsbnds(2),1)))&
 &              /(1.0*(cdofsbnds(2)-cdofsbnds(1)+1))
! this is approximate due to the normalizer being illegally averaged over
! and also due to periodic wraparound for mode 2
  else if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    helper = 0.5*it%nmbrs*(dot_product(it%sums(:,2),it%sums(:,5)) + it%nmbrs*sum(it%sums(:,4)) &
 &        - 2.0*dot_product(it%sums(:,1),it%sums(:,3)))/sum(it%sums(:,5))
    if (cdis_crit.eq.10) helper = 3.0*helper
  else if (cdis_crit.eq.8) then
    helper = (it%nmbrs*it%sqsum - dot_product(it%sums(:,1),it%sums(:,1)))
  end if
!
  if ((it%nmbrs.gt.1).AND.(helper.ge.0.0)) then
    it%diam = sqrt(2.0*helper/(1.0*it%nmbrs*(it%nmbrs-1.0)))
    it%radius = sqrt(helper/(1.0*it%nmbrs*it%nmbrs))
    it%quality = 1.0 - it%radius/targetsz
  else if (it%nmbrs.eq.1) then
    it%diam = 0.0
    it%quality = 1.0
    it%radius = 0.0
! because of FPE/vectorization, this number can occasionally be small and negative if the real variance is zero
  else if (abs(helper)/cradius.le.1.0e-8) then
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
!--------------------------------------------------------------------------------------
!
subroutine clusters_joint_diam(it1,it2,jointd,jointr)
!
  use clusters
!
  implicit none
!
  type (t_scluster) it1,it2,it3,it4
  RTYPE dummy,jointr,jointd
!
  if (it1%nmbrs.ge.it2%nmbrs) then
    call copy_cluster(it1,it3)
    call copy_cluster(it2,it4)
  else
    call copy_cluster(it1,it4)
    call copy_cluster(it2,it3)
  end if
  call join_clusters(it3,it4)
  dummy = cradius
  call cluster_calc_params(it3,dummy)
  jointd = it3%diam
  jointr = it3%radius
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine clusters_joint_radius(itl,its,jointr,jointd)
!
  use clusters
!
  implicit none
!
  type (t_scluster), INTENT(IN):: itl,its
  RTYPE, INTENT(OUT):: jointr,jointd
!
  integer netnum,k,catnrs
  RTYPE jsums(clstsz,sumssz),jsqsum,normer,normer2,normer3,vec1(calcsz),vec2(calcsz),vec3(calcsz),incr,helper
  RTYPE qrot(4),cto(3),ctn(3),tvec(3),vectmp(calcsz)
!
  netnum = itl%nmbrs + its%nmbrs
!
  if (cdis_crit.eq.1) then
    normer2 = 1.0/(1.0*itl%nmbrs)
    normer3 = 1.0/(1.0*its%nmbrs)
    normer = 1.0/(1.0*netnum)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    jsqsum = itl%sqsum + its%sqsum
    do k=1,clstsz
      incr = 0.0
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).lt.-180.0) then
        incr = its%nmbrs*360.0
        jsqsum = jsqsum + its%nmbrs*360.0*360.0 + 720.0*its%sums(k,1)
      end if
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).gt.180.0) then
        incr = -its%nmbrs*360.0
        jsqsum = jsqsum + its%nmbrs*360.0*360.0 - 720.0*its%sums(k,1)
      end if
      incr = incr + its%sums(k,1)
      jsums(k,1) = itl%sums(k,1) + incr
      if (normer*jsums(k,1).lt.-180.0) then
        jsqsum = jsqsum + netnum*360.0*360.0 + 720.0*jsums(k,1)
        jsums(k,1) = jsums(k,1) + netnum*360.0
      end if
      if (normer*jsums(k,1).gt.180.0) then
        jsqsum = jsqsum + netnum*360.0*360.0 - 720.0*jsums(k,1)
        jsums(k,1) = jsums(k,1) - netnum*360.0
      end if
    end do
  else if (cdis_crit.eq.2) then
    normer2 = 1.0/(1.0*itl%nmbrs)
    normer3 = 1.0/(1.0*its%nmbrs)
    normer = 1.0/(1.0*itl%nmbrs+its%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).lt.-180.0) then
        jsums(k,2) = itl%sums(k,2) + its%nmbrs*360.0*360.0 + 720.0*its%sums(k,1)
        jsums(k,4) = itl%sums(k,4) + its%sums(k,5)*360.0*360.0 + 720.0*its%sums(k,3)
        jsums(k,1) = itl%sums(k,1) + its%nmbrs*360.0
        jsums(k,3) = itl%sums(k,3) + its%sums(k,5)*360.0
      end if
      if (normer3*its%sums(k,1)-normer2*jsums(k,1).gt.180.0) then
        jsums(k,2) = itl%sums(k,2) + its%nmbrs*360.0*360.0 - 720.0*its%sums(k,1)
        jsums(k,4) = itl%sums(k,4) + its%sums(k,5)*360.0*360.0 - 720.0*its%sums(k,3)
        jsums(k,1) = itl%sums(k,1) - its%nmbrs*360.0
        jsums(k,3) = itl%sums(k,3) - its%sums(k,5)*360.0
      end if
    end do
    jsums(:,1:4) = jsums(:,1:4) + its%sums(:,1:4)
    jsums(:,5) = itl%sums(:,5) + its%sums(:,5)
    do k=1,clstsz
      if (normer*jsums(k,1).lt.-180.0) then
        jsums(k,2) = jsums(k,2) + netnum*360.0*360.0 + 720.0*jsums(k,1)
        jsums(k,4) = jsums(k,4) + jsums(k,5)*360.0*360.0 + 720.0*jsums(k,3)
        jsums(k,1) = jsums(k,1) + netnum*360.0
        jsums(k,3) = jsums(k,3) + jsums(k,5)*360.0
      end if
      if (normer*jsums(k,1).gt.180.0) then
        jsums(k,2) = jsums(k,2) + netnum*360.0*360.0 - 720.0*jsums(k,1)
        jsums(k,4) = jsums(k,4) + jsums(k,5)*360.0*360.0 - 720.0*jsums(k,3)
        jsums(k,1) = jsums(k,1) - netnum*360.0
        jsums(k,3) = jsums(k,3) - jsums(k,5)*360.0
      end if
    end do
  else if ((cdis_crit.eq.3).OR.(cdis_crit.eq.7).OR.(cdis_crit.eq.8)) then
    jsums(1:calcsz,1) = itl%sums(1:calcsz,1) + its%sums(1:calcsz,1)
    jsqsum = itl%sqsum + its%sqsum
  else if ((cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    jsums(1:clstsz,1:sumssz) = itl%sums(1:clstsz,1:sumssz) + its%sums(1:clstsz,1:sumssz)
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    if (align_for_clustering.EQV..true.) then
      catnrs = (cdofsbnds(4)-cdofsbnds(3)+1)/3
      vec1(:) = its%sums(1:calcsz,1)/(1.0*its%nmbrs)
      vec2(:) = itl%sums(1:calcsz,1)/(1.0*itl%nmbrs)
      call align_3D(catnrs,vec1(cdofsbnds(3):cdofsbnds(4)),vec2(cdofsbnds(3):cdofsbnds(4)),tvec,qrot,cto,ctn)
      do k=1,calcsz/3
        call quat_conjugatei(qrot,vec1(3*k-2:3*k),vectmp(1:3),cto)
        vec3(3*k-2:3*k) = vectmp(1:3)
        vec3(3*k-2:3*k) = vec3(3*k-2:3*k) + tvec(:)
      end do
      jsums(1:calcsz,1) = itl%sums(1:calcsz,1) + its%nmbrs*vec3(:)
      jsqsum = itl%sqsum + its%nmbrs*dot_product(vec3(cdofsbnds(1):cdofsbnds(2)),vec3(cdofsbnds(1):cdofsbnds(2)))
    else
      jsums(1:calcsz,1) = itl%sums(1:calcsz,1) + its%sums(1:calcsz,1)
      jsqsum = itl%sqsum + its%sqsum
    end if
  end if

  if ((cdis_crit.eq.7).OR.(cdis_crit.eq.3).OR.(cdis_crit.eq.1)) then
    helper = (netnum*jsqsum - dot_product(jsums(:,1),jsums(:,1)))/(1.0*calcsz)
!  this is approximate for all clusters exceeding size 2 if alignment is used
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    helper = 3.0*(netnum*jsqsum - dot_product(jsums(cdofsbnds(1):cdofsbnds(2),1),jsums(cdofsbnds(1):cdofsbnds(2),1)))&
 &              /(1.0*(cdofsbnds(2)-cdofsbnds(1)+1))
! this is approximate due to the normalizer being illegally averaged over
! and also due to periodic wraparound for mode 2
  else if ((cdis_crit.eq.2).OR.(cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    helper = 0.5*netnum*(dot_product(jsums(:,2),jsums(:,5)) + netnum*sum(jsums(:,4)) &
 &        - 2.0*dot_product(jsums(:,1),jsums(:,3)))/sum(jsums(:,5))
    if (cdis_crit.eq.10) helper = 3.0*helper
  else if (cdis_crit.eq.8) then
    helper = (netnum*jsqsum - dot_product(jsums(:,1),jsums(:,1)))
  end if
!
  if ((netnum.gt.1).AND.(helper.ge.0.0)) then
    jointd = sqrt(2.0*helper/(1.0*netnum*(netnum-1.0)))
    jointr = sqrt(helper/(1.0*netnum*netnum))
  end if
!
! now test explicitly
  normer = 0.0
  normer2 = 0.0
!
  if (cdis_crit.eq.2) then
    do k=1,clstsz
      vec1(2*k-1) = jsums(k,1)/(1.0*netnum)
      vec1(2*k) = jsums(k,5)/(1.0*netnum)
    end do
  else if (cdis_crit.eq.4) then
    do k=1,calcsz/3
      vec1(3*k-2) = jsums(2*k-1,1)/(1.0*netnum)
      vec1(3*k-1) = jsums(2*k,1)/(1.0*netnum)
      vec1(3*k) = jsums(2*k,5)/(1.0*netnum)
    end do
  else if (cdis_crit.eq.8) then
    do k=1,calcsz
      if (cl_imvec(k).gt.0.0) then
        vec1(k) = jsums(k,1)/(cl_imvec(k)*netnum)
      else
        vec1(k) = 0.0
      end if
    end do
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    vec1(1:clstsz) = jsums(:,1)/(1.0*netnum)
    vec1((clstsz+1):calcsz) = jsums(:,5)/(1.0*netnum)
  else
    vec1 = jsums(1:calcsz,1)/(1.0*netnum)
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
subroutine clusters_mean_diff(it1,it2,val)
!
  use clusters
!
  implicit none
!
  type (t_scluster) it1,it2
  integer catnrs,k
  RTYPE helper,val,vec1(calcsz),vec2(calcsz),tvec(3),qrot(4),cto(3),ctn(3),vectmp(3),hsqsum,shf
!
  shf = 0.0
  if (cdis_crit.eq.1) then
    vec2(:) = it1%sums(1:calcsz,1)/(1.0*it1%nmbrs) - it2%sums(1:calcsz,1)/(1.0*it2%nmbrs)
    vec1(:) = 0.0
    do k=1,calcsz
      if (vec2(k).gt.180.0) then
        vec1(k) = -360.0
      end if
      if (vec2(k).lt.-180.0) then
        vec1(k) = 360.0
      end if
    end do
    shf = 2.0*(dot_product(it1%sums(1:calcsz,1),vec1)*it2%nmbrs - &
 &             dot_product(it2%sums(1:calcsz,1),vec1)*it1%nmbrs) + &
 &        it1%nmbrs*it2%nmbrs*dot_product(vec1,vec1)
  else if (cdis_crit.eq.2) then
    vec2(1:clstsz) = it1%sums(1:clstsz,1)/(1.0*it1%nmbrs) - it2%sums(1:clstsz,1)/(1.0*it2%nmbrs)
    vec1(:) = 0.0
    do k=1,clstsz
      if (vec2(k).gt.180.0) then
        vec1(k) = -360.0
      end if
      if (vec2(k).lt.-180.0) then
        vec1(k) = 360.0
      end if
    end do
    shf = 2.0*(it2%nmbrs*dot_product(it1%sums(1:clstsz,3),vec1(1:clstsz)) - &
 &             it1%nmbrs*dot_product(it2%sums(1:clstsz,3),vec1(1:clstsz)) -  &
 &             dot_product(it2%sums(1:clstsz,1)*it1%sums(1:clstsz,5),vec1(1:clstsz)) + &
 &             dot_product(it1%sums(1:clstsz,1)*it2%sums(1:clstsz,5),vec1(1:clstsz))) +&
 &             it2%nmbrs*dot_product(vec1(1:clstsz)*vec1(1:clstsz),it1%sums(1:clstsz,5)) + &
 &             it1%nmbrs*dot_product(vec1(1:clstsz)*vec1(1:clstsz),it2%sums(1:clstsz,5))
  end if
  if ((cdis_crit.eq.3).OR.(cdis_crit.eq.1).OR.(cdis_crit.eq.7)) then
    helper = (shf + it1%nmbrs*it2%sqsum+it2%nmbrs*it1%sqsum - &
 &            2.0*dot_product(it1%sums(1:calcsz,1),it2%sums(1:calcsz,1)))/(1.0*calcsz*it1%nmbrs*it2%nmbrs)
    val = sqrt(helper)
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    if (align_for_clustering.EQV..true.) then
      catnrs = (cdofsbnds(4)-cdofsbnds(3)+1)/3
      vec1(:) = it1%sums(1:calcsz,1)/(1.0*it1%nmbrs)
      vec2(:) = it2%sums(1:calcsz,1)/(1.0*it2%nmbrs)
      call align_3D(catnrs,vec1(cdofsbnds(3):cdofsbnds(4)),vec2(cdofsbnds(3):cdofsbnds(4)),tvec,qrot,cto,ctn)
      do k=1,calcsz/3
        call quat_conjugatei(qrot,vec1(3*k-2:3*k),vectmp(1:3),cto)
        vec2(3*k-2:3*k) = vectmp(1:3)
        vec2(3*k-2:3*k) = vec2(3*k-2:3*k) + tvec(:)
      end do
      vec2(:) = vec2(:)*it1%nmbrs
      hsqsum = dot_product(vec2(cdofsbnds(1):cdofsbnds(2)),vec2(cdofsbnds(1):cdofsbnds(2)))/(1.0*it1%nmbrs)
      helper = 3.0*(it2%nmbrs*hsqsum+it1%nmbrs*it2%sqsum - &
 &   2.0*dot_product(it2%sums(cdofsbnds(1):cdofsbnds(2),1),vec2(cdofsbnds(1):cdofsbnds(2))))/(1.0*it1%nmbrs*it2%nmbrs)
      val = sqrt(helper/(1.0*(cdofsbnds(2)-cdofsbnds(1)+1)))
    else
      helper = 3.0*(it1%nmbrs*it2%sqsum+it2%nmbrs*it1%sqsum - &
 &   2.0*dot_product(it1%sums(cdofsbnds(1):cdofsbnds(2),1),it2%sums(cdofsbnds(1):cdofsbnds(2),1)))/(1.0*it1%nmbrs*it2%nmbrs)
      val = sqrt(helper/(1.0*(cdofsbnds(2)-cdofsbnds(1)+1)))
    end if
  else if (cdis_crit.eq.8) then
    helper = (it1%nmbrs*it2%sqsum+it2%nmbrs*it1%sqsum - 2.0*dot_product(it1%sums(1:calcsz,1),it2%sums(1:calcsz,1)))/&
 &           (1.0*it1%nmbrs*it2%nmbrs)
    val = sqrt(helper)
  else if ((cdis_crit.eq.4).OR.(cdis_crit.eq.2).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    helper = (shf + dot_product(it1%sums(:,2),it2%sums(:,5)) + dot_product(it2%sums(:,2),it1%sums(:,5)) + &
 &            it2%nmbrs*sum(it1%sums(:,4)) + it1%nmbrs*sum(it2%sums(:,4)) - &
 &            2.0*(dot_product(it1%sums(:,1),it2%sums(:,3)) + dot_product(it2%sums(:,1),it1%sums(:,3))))/&
 &            (it2%nmbrs*sum(it1%sums(:,5)) + it1%nmbrs*sum(it2%sums(:,5)))
    if (cdis_crit.eq.10) helper = 3.0*helper
    val = sqrt(helper)
  end if
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine clusters_sort(it,nba,upornot)
!
  use clusters
  use interfaces
!
  implicit none
!
  integer, INTENT(IN):: nba
  logical, INTENT(IN):: upornot
  type(t_scluster) it(nba)
!
  integer i,j,ii,jj,k,kk,klo,khi
  integer, ALLOCATABLE:: iv(:,:)
  logical lhlp
  type(t_scluster) tmpit
!
  if (nba.le.1) return
!
! sort according to size
  allocate(iv(nba,5))
  do i=1,nba
    iv(i,3) = i !(/ 1:nba /)
  end do
  iv(:,1) = it(1:nba)%nmbrs
  ii = 1
  jj = nba
  call merge_sort(nba,upornot,iv(:,1),iv(:,2),ii,jj,iv(:,3),iv(:,4))
  do i=1,nba
    iv(iv(i,4),5) = i
  end do
!
! now use the maps - we need to move stuff around continuously to preserve all necessary information
  do i=1,nba
    if (i.eq.iv(i,4)) cycle
    call copy_cluster(it(i),tmpit)
    call copy_cluster(it(iv(i,4)),it(i))
    call copy_cluster(tmpit,it(iv(i,4)))
    ii = iv(i,5)
    jj = iv(i,4)
    iv(i,5) = iv(jj,5)
    iv(jj,5) = ii
    iv(ii,4) = iv(i,4)
  end do
!
  if (resort_clustering.EQV..true.) then
    i = nba
    do while (i.ge.1)
      k = i-1
      khi = i
      do while (it(i)%nmbrs.eq.it(k)%nmbrs)
        k = k - 1
        if (k.le.0) exit
      end do
      klo = k+1
      if (khi.gt.(klo+1)) then
        do j=klo,khi
          iv(j,3) = j-klo+1
          iv(j,1) = it(j)%center
        end do
        ii = 1
        jj = khi-klo+1
        kk = jj
        lhlp = .true.
        call merge_sort(kk,lhlp,iv(klo:khi,1),iv(klo:khi,2),ii,jj,iv(klo:khi,3),iv(klo:khi,4))
        iv(klo:khi,4) = iv(klo:khi,4) + klo - 1
        do j=klo,khi
          iv(iv(j,4),5) = j
        end do
        do j=klo,khi
          if (j.eq.iv(j,4)) cycle
          call copy_cluster(it(j),tmpit)
          call copy_cluster(it(iv(j,4)),it(j))
          call copy_cluster(tmpit,it(iv(j,4)))
          ii = iv(j,5)
          jj = iv(j,4)
          iv(j,5) = iv(jj,5)
          iv(jj,5) = ii
          iv(ii,4) = iv(j,4)
        end do
      else if (khi.eq.(klo+1)) then
        if (it(khi)%center.lt.it(klo)%center) then
          call copy_cluster(it(klo),tmpit)
          call copy_cluster(it(khi),it(klo))
          call copy_cluster(tmpit,it(khi))
        end if
      end if
      i = k
      if (i.le.1) exit
    end do
  end if
!
  deallocate(iv)
!
end
!
!-------------------------------------------------------------------
!
! shorten list of basins based on allocation status of snaps
!
subroutine clusters_shorten(it,nba)
!
  use clusters
!
  implicit none
!
  integer i,nba,j
  type(t_scluster) it(nba)
!
  i = 1
  do while (i.le.nba)
    if (it(i)%nmbrs.eq.0) then
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
!----------------------------------------------------------------------------------
!
subroutine cluster_getcenter(it)
!
  use clusters
  use iounit
!
  implicit none
!
  type(t_scluster) it
  integer i
  RTYPE maxdis,rdv
!
  if ((it%nmbrs.eq.1).OR.(it%nmbrs.eq.2)) then
!   metric implies symmetry, so for n=2 it is safe to set the first as centroid
    it%center = minval(it%snaps(1:it%nmbrs))
    return
  end if
!
  maxdis = HUGE(maxdis)
  do i=1,it%nmbrs
    call snap_to_cluster_d(rdv,it,it%snaps(i))
    if (rdv.lt.maxdis) then
      it%center = it%snaps(i)
      maxdis = rdv
    end if
  end do
!
end
!
!----------------------------------------------------------------------------------
!
subroutine cluster_getcenter_heur(which,helper,alsz,mode)
!
  use clusters
  use iounit
!
  implicit none
!
  integer alsz,which,mode,j,k,jj,totnum,trystruc,maxtnum
  integer helper(alsz)
  RTYPE totdis,coreval
  
  helper(:) = 0
  do j=1,scluster(which)%nmbrs
    helper(scluster(which)%snaps(j)) = 1
  end do
  if (mode.eq.1) then
    maxtnum = 0
    coreval = cstored*chardcut
    do j=1,scluster(which)%nmbrs
      totnum = 0
      totdis = 0.0
      k = scluster(which)%snaps(j)
      do jj=1,cnblst(k)%nbs
        if (helper(cnblst(k)%idx(jj)).eq.1) then
          totdis = totdis + cnblst(k)%dis(jj)
          totnum = totnum + 1
        end if
      end do
      if ((totnum.gt.maxtnum).OR.((totnum.eq.maxtnum).AND.(totdis.lt.coreval))) then
        maxtnum = totnum
        coreval = totdis
        trystruc = k
      end if
    end do
    scluster(which)%center = trystruc
  else
    write(ilog,*) 'Fatal. Encountered unsupported mode in cluster_getcenter_heur(...). This is an &
 &omission bug.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
subroutine join_clusters(itl,its)
!
  use clusters
!
  implicit none
!
  integer k,catnrs
  RTYPE normer,normer2,normer3,vec1(calcsz),vec2(calcsz),vec3(calcsz),incr
  RTYPE vectmp(3),cto(3),tvec(3),ctn(3),qrot(4)
  type (t_scluster) its,itl
!
  if (cdis_crit.eq.1) then
    normer2 = 1.0/(1.0*itl%nmbrs)
    normer3 = 1.0/(1.0*its%nmbrs)
    normer = 1.0/(1.0*itl%nmbrs+its%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    itl%sqsum = itl%sqsum + its%sqsum
    do k=1,clstsz
      incr = 0.0
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).lt.-180.0) then
        incr = its%nmbrs*360.0
        itl%sqsum = itl%sqsum + its%nmbrs*360.0*360.0 + 720.0*its%sums(k,1)
      end if
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).gt.180.0) then
        incr = -its%nmbrs*360.0
        itl%sqsum = itl%sqsum + its%nmbrs*360.0*360.0 - 720.0*its%sums(k,1)
      end if
      incr = incr + its%sums(k,1)
      itl%sums(k,1) = itl%sums(k,1) + incr
      if (normer*itl%sums(k,1).lt.-180.0) then
        itl%sqsum = itl%sqsum + (itl%nmbrs+its%nmbrs)*360.0*360.0 + 720.0*itl%sums(k,1)
        itl%sums(k,1) = itl%sums(k,1) + (itl%nmbrs+its%nmbrs)*360.0
      end if
      if (normer*itl%sums(k,1).gt.180.0) then
        itl%sqsum = itl%sqsum + (itl%nmbrs+its%nmbrs)*360.0*360.0 - 720.0*itl%sums(k,1)
        itl%sums(k,1) = itl%sums(k,1) - (itl%nmbrs+its%nmbrs)*360.0
      end if
    end do
  else if (cdis_crit.eq.2) then
    normer2 = 1.0/(1.0*itl%nmbrs)
    normer3 = 1.0/(1.0*its%nmbrs)
    normer = 1.0/(1.0*itl%nmbrs+its%nmbrs)
!   this is approximate: the average mustn't shift very much and/or the points need to be well-clustered
    do k=1,clstsz
      incr = 0.0
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).lt.-180.0) then
        itl%sums(k,2) = itl%sums(k,2) + its%nmbrs*360.0*360.0 + 720.0*its%sums(k,1)
        itl%sums(k,4) = itl%sums(k,4) + its%sums(k,5)*360.0*360.0 + 720.0*its%sums(k,3)
        itl%sums(k,1) = itl%sums(k,1) + its%nmbrs*360.0
        itl%sums(k,3) = itl%sums(k,3) + its%sums(k,5)*360.0
      end if
      if (normer3*its%sums(k,1)-normer2*itl%sums(k,1).gt.180.0) then
        itl%sums(k,2) = itl%sums(k,2) + its%nmbrs*360.0*360.0 - 720.0*its%sums(k,1)
        itl%sums(k,4) = itl%sums(k,4) + its%sums(k,5)*360.0*360.0 - 720.0*its%sums(k,3)
        itl%sums(k,1) = itl%sums(k,1) - its%nmbrs*360.0
        itl%sums(k,3) = itl%sums(k,3) - its%sums(k,5)*360.0
      end if
    end do
    itl%sums(:,:) = itl%sums(:,:) + its%sums(:,:)
    do k=1,clstsz
      if (normer*itl%sums(k,1).lt.-180.0) then
        itl%sums(k,2) = itl%sums(k,2) + (itl%nmbrs+its%nmbrs)*360.0*360.0 + 720.0*itl%sums(k,1)
        itl%sums(k,4) = itl%sums(k,4) + itl%sums(k,5)*360.0*360.0 + 720.0*itl%sums(k,3)
        itl%sums(k,1) = itl%sums(k,1) + (itl%nmbrs+its%nmbrs)*360.0
        itl%sums(k,3) = itl%sums(k,3) + itl%sums(k,5)*360.0
      end if
      if (normer*itl%sums(k,1).gt.180.0) then
        itl%sums(k,2) = itl%sums(k,2) + (itl%nmbrs+its%nmbrs)*360.0*360.0 - 720.0*itl%sums(k,1)
        itl%sums(k,4) = itl%sums(k,4) + itl%sums(k,5)*360.0*360.0 - 720.0*itl%sums(k,3)
        itl%sums(k,1) = itl%sums(k,1) - (itl%nmbrs+its%nmbrs)*360.0
        itl%sums(k,3) = itl%sums(k,3) - itl%sums(k,5)*360.0
      end if
    end do
  else if ((cdis_crit.eq.3).OR.(cdis_crit.eq.7).OR.(cdis_crit.eq.8)) then
    itl%sums(1:calcsz,1) = itl%sums(1:calcsz,1) + its%sums(1:calcsz,1)
    itl%sqsum = itl%sqsum + its%sqsum
  else if ((cdis_crit.eq.4).OR.(cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    itl%sums(1:clstsz,1:sumssz) = itl%sums(1:clstsz,1:sumssz) + its%sums(1:clstsz,1:sumssz)
  else if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
    if (align_for_clustering.EQV..true.) then
      catnrs = (cdofsbnds(4)-cdofsbnds(3)+1)/3
      vec1(:) = its%sums(1:calcsz,1)/(1.0*its%nmbrs)
      vec2(:) = itl%sums(1:calcsz,1)/(1.0*itl%nmbrs)
      call align_3D(catnrs,vec1(cdofsbnds(3):cdofsbnds(4)),vec2(cdofsbnds(3):cdofsbnds(4)),tvec,qrot,cto,ctn)
      do k=1,calcsz/3
        call quat_conjugatei(qrot,vec1(3*k-2:3*k),vectmp(1:3),cto)
        vec3(3*k-2:3*k) = vectmp(1:3)
        vec3(3*k-2:3*k) = vec3(3*k-2:3*k) + tvec(:)
      end do
      itl%sums(1:calcsz,1) = itl%sums(1:calcsz,1) + its%nmbrs*vec3(:)
      itl%sqsum = itl%sqsum + its%nmbrs*dot_product(vec3(cdofsbnds(1):cdofsbnds(2)),vec3(cdofsbnds(1):cdofsbnds(2)))
    else
      itl%sums(1:calcsz,1) = itl%sums(1:calcsz,1) + its%sums(1:calcsz,1)
      itl%sqsum = itl%sqsum + its%sqsum
    end if
  end if
!
  itl%nmbrs = itl%nmbrs + its%nmbrs
  do while (itl%nmbrs.gt.itl%alsz)
    call scluster_resize(itl)
  end do
  itl%snaps(itl%nmbrs-its%nmbrs+1:itl%nmbrs) = its%snaps(1:its%nmbrs)
  its%nmbrs = 0
  its%alsz = 0
  deallocate(its%snaps)
  deallocate(its%sums)
!
end
!
!-----------------------------------------------------------------------------------------
!
subroutine merge_clusters(maxj,i,helper,alsz,isup,isadd)
!
  use clusters
  use interfaces
!
  implicit none
!
  integer maxj,i,ii,jj,kk,compnum,alsz
  integer helper(alsz)
  integer, ALLOCATABLE:: iv2(:)
  logical candid,isup,isadd
!
  kk = 1
  compnum = 0
  do ii=1,scluster(i)%nmbrs
    candid = .false.
    do while (scluster(maxj)%snaps(kk).lt.scluster(i)%snaps(ii))
      kk = kk + 1
      if (kk.eq.scluster(maxj)%nmbrs) exit
    end do
    if (scluster(maxj)%snaps(kk).eq.scluster(i)%snaps(ii)) then
      if (isadd.EQV..false.) then
        compnum = compnum + 1
        helper(compnum) = scluster(i)%snaps(ii)
      end if
    else if ((kk.eq.scluster(maxj)%nmbrs).AND.(scluster(maxj)%snaps(kk).lt.scluster(i)%snaps(ii))) then
      if (isadd.EQV..true.) then
        compnum = compnum + 1
        helper(compnum) = scluster(i)%snaps(ii)
      end if
      exit
    else
      if (isadd.EQV..true.) then
        compnum = compnum + 1
        helper(compnum) = scluster(i)%snaps(ii)
      end if
    end if
    kk = kk - 1
  end do
  if (isadd.EQV..true.) then
    allocate(iv2(scluster(maxj)%nmbrs))
    iv2(:) = scluster(maxj)%snaps(:)
    deallocate(scluster(maxj)%snaps)
    allocate(scluster(maxj)%snaps(scluster(maxj)%nmbrs+compnum))
    scluster(maxj)%snaps(1:scluster(maxj)%nmbrs) = iv2(:)
    scluster(maxj)%snaps(scluster(maxj)%nmbrs+1:scluster(maxj)%nmbrs+compnum) = helper(1:compnum)
    scluster(maxj)%nmbrs = scluster(maxj)%nmbrs + compnum
    deallocate(iv2)
  else
    deallocate(scluster(maxj)%snaps)
    allocate(scluster(maxj)%snaps(compnum))
    scluster(maxj)%snaps(1:compnum) = helper(1:compnum)
    scluster(maxj)%nmbrs = compnum
  end if
  allocate(iv2(scluster(maxj)%nmbrs))
  iv2(:) = scluster(maxj)%snaps(1:scluster(maxj)%nmbrs)
  ii = 1
  jj = scluster(maxj)%nmbrs
  call merge_sort(ldim=scluster(maxj)%nmbrs,up=isup,list=iv2,&
 &                olist=scluster(maxj)%snaps(1:scluster(maxj)%nmbrs),ilo=ii,ihi=jj)
  deallocate(iv2)
  deallocate(scluster(i)%snaps)
!
end
!
!----------------------------------------------------------------------------------
!
subroutine cluster_sortsnaps(it)
!
  use clusters
  use interfaces
!
  implicit none
!
  type (t_scluster) it
  integer ii,jj
  integer, ALLOCATABLE:: iv1(:) ! iv1(it%nmbrs)
  logical atrue
!
  atrue = .true.
!
  allocate(iv1(it%nmbrs))
  iv1(:) = it%snaps(1:it%nmbrs)
  ii = 1
  jj = it%nmbrs
  call merge_sort(ldim=it%nmbrs,up=atrue,list=iv1,olist=it%snaps(1:it%nmbrs),ilo=ii,ihi=jj)
  deallocate(iv1)
!
end
!
!---------------------------------------------------------------------------------
!
subroutine cluster_removeframes(it,frlist,lilen)
!
  use clusters
!
  implicit none
!
  type(t_scluster) it
  integer lilen,icnt,jcnt,ncnt
  integer frlist(lilen)
  integer, ALLOCATABLE:: tmplist(:)
!
  if (it%nmbrs.le.0) return
  allocate(tmplist(it%nmbrs))
  jcnt = 1
  icnt = 1
  ncnt = 0
  do while (icnt.le.it%nmbrs)
    do while (frlist(jcnt).lt.it%snaps(icnt))
      jcnt = jcnt + 1
      if (jcnt.gt.lilen) exit
    end do
    if (jcnt.gt.lilen) exit
    if (it%snaps(icnt).eq.frlist(jcnt)) then
      tmplist(icnt+1:it%nmbrs) = it%snaps(icnt+1:it%nmbrs)
      it%nmbrs = it%nmbrs-1
      it%snaps(icnt:it%nmbrs) = tmplist(icnt+1:it%nmbrs+1)
      icnt = icnt - 1
      ncnt = ncnt + 1
    end if
    if (ncnt.eq.lilen) exit
    if (icnt.ge.it%nmbrs) exit
    icnt = icnt + 1
  end do 
  if (it%nmbrs.eq.0) deallocate(it%snaps)
  deallocate(tmplist)
!
end
!
!---------------------------------------------------------------------------------
!
subroutine cluster_transferframes(it1,it2,frlist,lilen)
!
  use clusters
!
  implicit none
!
  type(t_scluster) it1,it2,it3
  integer lilen,i,j
  integer frlist(lilen)
  RTYPE rdv
!
  rdv = cradius
  if (lilen.eq.it1%nmbrs) then
    call join_clusters(it2,it1)
    call cluster_sortsnaps(it2)
    call cluster_calc_params(it2,rdv)
    return
  end if
!
  it3%nmbrs = 0
  it3%alsz = 0
  do i=1,lilen
    call cluster_addsnap(it3,frlist(i),rdv)
  end do
  call join_clusters(it2,it3)
!
  it3%nmbrs = 0 ! it3 is wiped by call to join_clusters
  j = 1
  
  do i=1,it1%nmbrs
    if (j.le.lilen) then
      do while (frlist(j).lt.it1%snaps(i))
        if (frlist(j).eq.it1%snaps(i)) exit
        j = j + 1
        if (j.gt.lilen) exit
      end do
      if (j.le.lilen) then
        if (frlist(j).ne.it1%snaps(i)) call cluster_addsnap(it3,it1%snaps(i),rdv)
      else
        call cluster_addsnap(it3,it1%snaps(i),rdv)
      end if
    else
      call cluster_addsnap(it3,it1%snaps(i),rdv)
    end if
  end do
!
  call copy_cluster(it3,it1)
  call cluster_calc_params(it1,rdv)
!
  call cluster_sortsnaps(it2)
!
end
!
!---------------------------------------------------------------------------------
!
function profile_min(cutv,setsz,i,mins,maxs,ints,siv)
!
  implicit none
!
  integer, INTENT(INOUT):: siv
  integer mins,maxs,ints,setsz,sumle,sumri,sum1,hsiv,i,j,jj
  logical profile_min
  integer cutv(setsz)
!
  profile_min = .false.
  do siv=mins,min(min(maxs,i-1),setsz-i),ints
    hsiv = floor(0.5*siv)
    sumle = sum(cutv(i-siv:i-1))
    sum1 = sum(cutv(i-hsiv:i+hsiv))
    sumri = sum(cutv(i+1:i+siv)) 
    if ((sumle.gt.sum1).AND.(sumri.gt.sum1).AND.((sumri+sumle).gt.2*(sum1+4*siv))) then
      profile_min = .true.
      if (sum(cutv(i-siv:i-hsiv-2)).le.sum(cutv(i-hsiv:i-1))) profile_min = .false.
      if (sum(cutv(i+hsiv+2:i+siv)).le.sum(cutv(i+1:i+hsiv))) profile_min = .false.
      jj = cutv(i)
      do j=i-siv,i-1
        if (cutv(j).lt.jj) profile_min = .false.
      end do
      do j=i+1,i+siv
        if (cutv(j).le.jj) profile_min = .false.
      end do
    end if
    if (profile_min.EQV..true.) return
  end do
!
  return
!
end
!
!-------------------------------------------------------------------------------
!
  subroutine scluster_resize(it)
!
  use clusters
!
  implicit none
!
  type(t_scluster) it
  integer oldsz
  integer, ALLOCATABLE:: tmpa(:) ! tmpa(it%alsz),oldsz
!
  if (allocated(it%snaps).EQV..false.) then
    allocate(it%snaps(2))
    it%alsz = 2
    if (it%alsz.ge.it%nmbrs) return
  end if
  oldsz = it%alsz
  allocate(tmpa(oldsz))
  tmpa(:) = it%snaps(:)
  deallocate(it%snaps)
  if (it%alsz.lt.100) then
    it%alsz = it%alsz*4
  else
    it%alsz = min(cstored+1,it%alsz*2)
  end if
  allocate(it%snaps(it%alsz))
  it%snaps(1:oldsz) = tmpa(:)
  deallocate(tmpa)
!
end
!
!-------------------------------------------------------------------------------
!
  subroutine scluster_resizech(it)
!
  use clusters
!
  implicit none
!
  type(t_scluster) it
  integer oldsz
  integer, ALLOCATABLE:: tmpa(:) ! tmpa(it%chalsz),oldsz
!
  if (allocated(it%children).EQV..false.) then
    allocate(it%children(2))
    it%chalsz = 2
    if (it%chalsz.ge.it%nchildren) return
  end if
  oldsz = it%chalsz
  allocate(tmpa(oldsz))
  tmpa(:) = it%children(:)
  deallocate(it%children)
  if (it%chalsz.lt.100) then
    it%chalsz = it%chalsz*4
  else
    it%chalsz = min(cstored+1,it%chalsz*2)
  end if
  allocate(it%children(it%chalsz))
  it%children(1:oldsz) = tmpa(:)
  deallocate(tmpa)
!
end
!
!-------------------------------------------------------------------------------
!
! expand graph nb lists - we assume unallocated member arrays to be dealt with
! elsewhere
! like scluster_resize must never be called with allocation size currently zero
!
  subroutine scluster_resizenb(it)
!
  use clusters
!
  implicit none
!
  type(t_scluster) it
  integer tmpa(it%nbalsz,2),oldsz,newsz
  RTYPE tmpb(it%nbalsz,4)
!
  oldsz = it%nbalsz
  newsz = oldsz*2
  if (allocated(it%lstnb).EQV..true.) then
    tmpa(:,1) = it%lstnb(1:oldsz)
    deallocate(it%lstnb)
    allocate(it%lstnb(newsz))
    it%lstnb(1:oldsz) = tmpa(:,1)
  end if
!
  if (allocated(it%map).EQV..true.) then
    tmpa(:,1) = it%map(1:oldsz)
    deallocate(it%map)
    allocate(it%map(newsz))
    it%map(1:oldsz) = tmpa(:,1)
  end if
!
  if (allocated(it%flwnb).EQV..true.) then
    tmpa(:,:) = it%flwnb(1:oldsz,:)
    deallocate(it%flwnb)
    allocate(it%flwnb(newsz,2))
    it%flwnb(1:oldsz,:) = tmpa(:,:)
  end if
!
  if (allocated(it%wghtsnb).EQV..true.) then
    tmpa(:,:) = it%wghtsnb(1:oldsz,:)
    deallocate(it%wghtsnb)
    allocate(it%wghtsnb(newsz,2))
    it%wghtsnb(1:oldsz,:) = tmpa(:,:)
  end if
!
  if (allocated(it%lensnb).EQV..true.) then
    tmpb(:,1:2) = it%lensnb(1:oldsz,:)
    deallocate(it%lensnb)
    allocate(it%lensnb(newsz,2))
    it%lensnb(1:oldsz,:) = tmpb(:,1:2)
  end if
!
  if (allocated(it%fewtsnb).EQV..true.) then
    tmpb(:,:) = it%fewtsnb(1:oldsz,:)
    deallocate(it%fewtsnb)
    allocate(it%fewtsnb(newsz,4))
    it%fewtsnb(1:oldsz,:) = tmpb(:,:)
  end if
!
  it%nbalsz = newsz
!
end
!
!-------------------------------------------------------------------------------
!
!subroutine scope_estimate_for_clustering()
!!
!  use clusters
!  use iounit
!  use, INTRINSIC:: ISO_C_BINDING
!!
!  implicit none
!!
!  integer i
!  integer(KIND=C_INT) C_SIZEOF
!  real(KIND=4) flt
!  RTYPE sz_cl(3),szcd,sznb
!!
!! estimates for final scluster structure 
! ! sz_cl(1) = ((clstsz*sumssz+5 + 2*log(1.0*cmaxsnaps))*C_SIZEOF(sz_cl(1)) + (14+4*log(1.0*cmaxsnaps))*C_SIZEOF(i))*0.1*cmaxsnaps
!  sz_cl(2) = sz_cl(1)/10.
!  sz_cl(3) = sz_cl(1)/100.
!! cludata
!!  szcd = C_SIZEOF(flt)*calcsz*cmaxsnaps
!! neighbor list
!  if ((cmode.eq.3).OR.((cmode.eq.4).AND.(cprogindex.eq.1))) then
!    sznb = 0.0 !(C_SIZEOF(flt)+C_SIZEOF(i))*cmaxsnaps*cmaxsnaps/2. + 5*C_SIZEOF(i)*cmaxsnaps
!  else
!    sznb = 0.0
!  end if
!! 
!!
! 65 format('Postprocessing could see an estimated memory usage of ',g12.5,' GB.')
!!  write(ilog,65) (szcd+sznb+sz_cl(1))/1.0e9
!!  call fexit()
!!
!end
!
!-------------------------------------------------------------------------------
!
subroutine quality_of_clustering(ncls,it,radcrit,quals)
!
  use clusters
  use iounit
!
  implicit none
!
  integer ncls,normer,alln,i
  type (t_scluster) it(ncls)
  RTYPE val1,val2,quals(3),unif,radcrit
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
!-------------------------------------------------------------------------------
!
subroutine copy_cluster(cfrom,cto)
!
  use clusters
!
  implicit none
!
  integer tmp1
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
    allocate(cto%children(cfrom%chalsz))
    cto%children(1:cfrom%chalsz) = cfrom%children(1:cfrom%chalsz)
  end if
  if (allocated(cfrom%wghtsnb).EQV..true.) then
    if (allocated(cto%wghtsnb).EQV..true.) deallocate(cto%wghtsnb)
    allocate(cto%wghtsnb(cfrom%nbalsz,2))
    cto%wghtsnb(1:cfrom%nbalsz,:) = cfrom%wghtsnb(1:cfrom%nbalsz,:)
  end if
  if (allocated(cfrom%fewtsnb).EQV..true.) then
    if (allocated(cto%fewtsnb).EQV..true.) deallocate(cto%fewtsnb)
    allocate(cto%fewtsnb(cfrom%nbalsz,4))
    cto%fewtsnb(1:cfrom%nbalsz,:) = cfrom%fewtsnb(1:cfrom%nbalsz,:)
  end if
  if (allocated(cfrom%lensnb).EQV..true.) then
    if (allocated(cto%lensnb).EQV..true.) deallocate(cto%lensnb)
    allocate(cto%lensnb(cfrom%nbalsz,2))
    cto%lensnb(1:cfrom%nbalsz,:) = cfrom%lensnb(1:cfrom%nbalsz,:)
  end if
  if (allocated(cfrom%map).EQV..true.) then
    tmp1 = size(cfrom%map)
    if (allocated(cto%map).EQV..true.) deallocate(cto%map)
    allocate(cto%map(tmp1))
    cto%map(1:tmp1) = cfrom%map(1:tmp1)
  end if
  if (allocated(cfrom%lstnb).EQV..true.) then
    tmp1 = size(cfrom%lstnb)
    if (allocated(cto%lstnb).EQV..true.) deallocate(cto%lstnb)
    allocate(cto%lstnb(tmp1))
    cto%lstnb(1:tmp1) = cfrom%lstnb(1:tmp1)
  end if
  if (allocated(cfrom%flwnb).EQV..true.) then
    if (allocated(cto%flwnb).EQV..true.) deallocate(cto%flwnb)
    allocate(cto%flwnb(cfrom%nbalsz,2))
    cto%flwnb(1:cfrom%nbalsz,:) = cfrom%flwnb(1:cfrom%nbalsz,:)
  end if
  if (allocated(cfrom%sums).EQV..true.) then
    if (allocated(cto%sums).EQV..true.) deallocate(cto%sums)
    allocate(cto%sums(clstsz,sumssz))
    cto%sums(1:clstsz,1:sumssz) = cfrom%sums(1:clstsz,1:sumssz)
  end if
  cto%nmbrs = cfrom%nmbrs
  cto%center = cfrom%center
  cto%geni = cfrom%geni
  cto%genidx = cfrom%genidx
  cto%alsz = cfrom%alsz
  cto%nb = cfrom%nb
  cto%ldis = cfrom%ldis
  cto%active = cfrom%active
  cto%rflw = cfrom%rflw
  cto%inscc = cfrom%inscc
  cto%diam = cfrom%diam
  cto%radius = cfrom%radius
  cto%nodewt(:) = cfrom%nodewt(:)
  cto%quality = cfrom%quality
  cto%sqsum = cfrom%sqsum
  cto%nbalsz = cfrom%nbalsz
  cto%nchildren = cfrom%nchildren
  cto%chalsz = cfrom%chalsz
  cto%parent = cfrom%parent
!
end
!
!-------------------------------------------------------------------------------
!
subroutine scluster_resizelst(currentalsz,it)
!
  use clusters
!
  implicit none
!
  integer i,oldsz,currentalsz
  type(t_scluster), ALLOCATABLE, INTENT(IN OUT):: it(:)
  type(t_scluster), ALLOCATABLE:: tmpcls(:)
!
  oldsz = currentalsz
  if (oldsz.gt.0) then
    allocate(tmpcls(oldsz))
  end if
  do i=1,oldsz
    call copy_cluster(it(i),tmpcls(i))
  end do
!
  if (currentalsz.le.0) then
    currentalsz = 10
  else
    deallocate(it)
    currentalsz = currentalsz*2
  end if
!
  allocate(it(currentalsz))
!
  do i=1,oldsz
    call copy_cluster(tmpcls(i),it(i))
  end do
  it(oldsz+1:currentalsz)%nmbrs = 0
  it(oldsz+1:currentalsz)%nchildren = 0
  it(oldsz+1:currentalsz)%nb = 0
  it(oldsz+1:currentalsz)%alsz = 0  
  it(oldsz+1:currentalsz)%nbalsz = 0
  it(oldsz+1:currentalsz)%chalsz = 0
  it(oldsz+1:currentalsz)%parent = 0
  it(oldsz+1:currentalsz)%geni = 0
!
  if (allocated(tmpcls).EQV..true.) deallocate(tmpcls)
!
end
!
!-------------------------------------------------------------------------------
!
#ifdef LINK_LAPACK
!
subroutine get_principal_components(sz1,sz2,sz3)
!
  use clusters
  use atoms
  use iounit
#ifdef ENABLE_MPI
  use mpistuff
#endif
!
  implicit none
!
  integer i,j,sizen,ipca,sz1,sz2,sz3,effsz1,lda,ldvt,ldu,lwork,info,freeunit,ii,jj
  RTYPE normer,qrot(4),tvec(3),cto(3),ctn(3)
  RTYPE, ALLOCATABLE:: bla1(:,:),bla2(:),bla3(:,:),bla4(:),bla5(:,:),bla6(:),bla7(:,:),bla8(:),bla9(:,:)
  character(1) cu,cvt
  character(MAXSTRLEN) fn
  logical exists
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
!
  sizen = sz2
!
  if (sz2.le.sz3) then
    write(ilog,*) 'Warning. Skipping computation of principal components due to &
 &limited amount of data (not more snapshots than dimensions).'
    return
  end if
!
! for options 4, 9, and 10, create an average static weight from the local weights and scale data
  effsz1 = sz3
  if (cdis_crit.eq.4) then
    allocate(bla8(sz1/3))
    allocate(bla7(effsz1,sz2)) 
    do i=1,sz1/3
      bla8(i) = sum(cludata(3*i,1:sz2))
    end do
    normer = (1.0*(sz1/3))/sum(bla8(:))
    bla8(:) = sqrt(normer*bla8(:))
    do i=1,sz1/3
      bla7(2*i-1,1:sz2) = cludata(3*i-2,1:sz2)*bla8(i)
      bla7(2*i,1:sz2) = cludata(3*i-1,1:sz2)*bla8(i)
    end do
  else if ((cdis_crit.eq.9).OR.(cdis_crit.eq.10)) then
    allocate(bla8(sz3))
    allocate(bla7(effsz1,sz2)) 
    do i=1,sz3
      bla8(i) = sum(cludata(sz3+i,1:sz2))
    end do
    normer = (1.0*sz3)/sum(bla8(:))
    bla8(:) = sqrt(normer*bla8(:))
    effsz1 = 0
    do i=1,sz3
      if (bla8(i).gt.0.0) then
        effsz1 = effsz1 + 1
        bla7(effsz1,1:sz2) = cludata(i,1:sz2)*bla8(i)
      end if
    end do
    if (effsz1.ne.sz3) then
      write(ilog,*) 'Warning. Dimensions with zero weight have been removed for PCA/tICA computation.'
      reduced_dim_clustering = min(reduced_dim_clustering,effsz1)
    end if
! for static weights, skip the first step
  else if (cdis_crit.eq.8) then
    allocate(bla7(sz1,sz2))
    do i=1,sz2
      bla7(1:sz1,i) = sqrt(1.0*sz3)*cludata(1:sz1,i)*cl_imvec(1:sz1)
    end do
! in all other allowed cases, simply copy
  else
    allocate(bla7(sz1,sz2))
    bla7(:,:) = cludata(1:sz1,1:sz2)
  end if
  if (allocated(bla8).EQV..true.) deallocate(bla8)
!
! for RMSD with alignment, prealign trajectory data to last snapshot (may be a poor choice of course)
  if (((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)).AND.(align_for_clustering.EQV..true.)) then
    allocate(bla8(effsz1))
    if (effsz1.ne.calcsz) then
      write(ilog,*) 'Fatal. When using RMSD with alignment, principal components were passed incompatible data. This &
 &is a bug.'
      call fexit()
    end if
    ii = (cdofsbnds(4)-cdofsbnds(3)+1)/3
    do jj=1,sz2-1
      call align_3D(ii,bla7(cdofsbnds(3):cdofsbnds(4),jj),bla7(cdofsbnds(3):cdofsbnds(4),sz2),tvec,qrot,cto,ctn)
      do j=1,effsz1/3
        call quat_conjugatei(qrot,bla7(3*j-2:3*j,jj),bla8(3*j-2:3*j),cto)
        bla8(3*j-2:3*j) = bla8(3*j-2:3*j) + tvec(:)
      end do
      bla7(1:effsz1,jj) = bla8(:)
    end do
    deallocate(bla8)
    if (cdis_crit.eq.6) then
      effsz1 = cdofsbnds(2)-cdofsbnds(1)+1 ! cdofsbnds(1) must always be 1
    end if
  end if
!
! subtract sample mean (this may be redundant)
  allocate(bla5(effsz1,sz2))
  allocate(bla6(effsz1))
  bla6(1:effsz1) = (1.0/(1.0*sz2))*sum(bla7(1:effsz1,1:sz2),dim=2)
  do i=1,sz2
    bla5(1:effsz1,i) = bla7(1:effsz1,i) - bla6(1:effsz1)
  end do
!
! standard PCA
  if (pcamode.le.3) then
    ldu = effsz1
    allocate(bla1(ldu,ldu))
    lda = min(sz2,effsz1)
    allocate(bla2(lda))
    ldvt = 1
    allocate(bla3(ldvt,ldvt))
    lwork = 5*min(sz2,effsz1)
    allocate(bla4(lwork))
!   request only left EVs, PCA via SVD
    cvt = 'N'
    cu = 'A'
    jj = effsz1
    call dgesvd(cu,cvt,jj,sizen,bla5,lda,bla2,bla1,ldu,bla3,ldvt,bla4,lwork,info)
! tICA
  else if ((pcamode.eq.4).OR.(pcamode.eq.5)) then
!   compute covariance matrix
    allocate(bla1(effsz1,effsz1))
    bla1(:,:) = 0.
    do i=1,sz2
      do j=1,effsz1
        bla1(:,j) = bla1(:,j) + bla5(1:effsz1,i)*bla5(j,i)
      end do
    end do
    bla1(:,:) = bla1(:,:)/(1.0*(sz2)) ! N or N-1 estimator?
!   time-lagged covariance matrix
    allocate(bla3(effsz1,effsz1))
    allocate(bla9(effsz1,effsz1))
    bla3(:,:) = 0.
    do i=1+cwacftau,sz2
      do j=1,effsz1
        bla9(:,j) = bla5(1:effsz1,i)*bla5(j,(i-cwacftau))
      end do
      bla3(:,:) = bla3(:,:) + 0.5*(bla9(:,:) + transpose(bla9(:,:))) ! symmetrization
    end do
    bla3(:,:) = bla3(:,:)/(1.0*(sz2-cwacftau)) ! N or N-1 estimator?
    lda = effsz1 ! min(sz2,effsz1)
    allocate(bla2(lda))
    lwork = 5*min(sz2,effsz1)
    allocate(bla4(lwork))
    cvt = 'V'
    cu = 'U'
    i = 1
    ii = effsz1
    jj = effsz1
    call dsygv(i,cvt,cu,lda,bla3(:,:),ii,bla1(:,:),jj,bla2(:),bla4(:),lwork,info) ! ascending order of eigenvalues -> switch
    bla4(1:effsz1) = bla2(1:effsz1)
    do i=1,effsz1
      bla2(i) = bla4(effsz1-i+1) ! eigenvalues switched
    end do
    bla9(:,:) = bla3(:,:)
    do i=1,effsz1
      bla3(:,i) = bla9(:,effsz1-i+1) ! eigenvectors switched
    end do
    deallocate(bla9)
  end if
  if (info.ne.0) then
    write(ilog,*) 'Warning. The linear algebra routine required for the selected option for keyword FMCSC_PCAMODE returned &
 &a nonzero error code (',info,'). Results in PRINCIPAL_COMPONENTS.dat and PRINCIPAL_COMPONENTS.evs may be meaningless.'
  end if
!
! now use eigenvectors as per user request
  if ((pcamode.eq.3).OR.(pcamode.eq.5)) then ! transform, write both, potentially reduce dimensionality for subsequent clustering
    deallocate(bla5)
    allocate(bla5(sz2,effsz1))
    bla5(:,:) = 0.0
!   for reasons unknown, all attempts to use matmul or dgemm failed here for large matrices whether with intrinsic, LAPACK, or MKL
!   speculation has it that this due to a size violation on a hidden temporary array or simply a stack problem
    jj = effsz1
    if (pcamode.eq.3) then
      call matmul2(sz2,effsz1,jj,bla7(1:effsz1,1:sz2),bla1,bla5)
!      bla5(1:sz2,1:sz1) = matmul(transpose(cludata(1:sz1,1:sz2)),bla1)
!      call dgemm(cvt,cvt,sizen,sizem,sizem,aone,transpose(indata(1:sz1,1:sz2)),sizen,bla1,sizem,azero,bla5,sizen)
!
!     compute the resultant variances for the transformed space (two-pass)
      normer = (1.0/(sz2*1.0))
      bla4(1:effsz1) = normer*sum(bla5(1:sz2,1:effsz1),dim=1)
      do i=1,effsz1
        bla6(i) = normer*sum((bla5(1:sz2,i)-bla4(i))**2)
      end do
    else if (pcamode.eq.5) then
     call matmul2(sz2,effsz1,jj,bla7(1:effsz1,1:sz2),bla3(:,:),bla5(:,:))
!      bla5(1:sz2,1:sz1) = matmul(transpose(cludata(1:sz1,1:sz2)),bla1)
!      call dgemm(cvt,cvt,sizen,sizem,sizem,aone,transpose(indata(1:sz1,1:sz2)),sizen,bla1,sizem,azero,bla5,sizen)
!
!     compute the resultant autocorrelation values at fixed lag for the transformed space
      ii = sz2
      do i=1,effsz1
        call vacf_fixtau(bla5(1:sz2,i),ii,cwacftau,bla6(i))
      end do
    end if
!
#ifdef ENABLE_MPI
    call int2str(myrank,nod,re_aux(10))
    if (use_REMC.EQV..true.) then
      fn = 'N_'//nod(1:re_aux(10))//'_PRINCIPAL_COMPONENTS.dat'
    else
      fn = 'PRINCIPAL_COMPONENTS.dat'
    end if
#else
    fn = 'PRINCIPAL_COMPONENTS.dat'
#endif
!
    call strlims(fn,ii,jj)
    inquire(file=fn(ii:jj),exist=exists)
    if(exists) then
      ipca = freeunit()
      open(unit=ipca,file=fn(ii:jj),status='old',position='append')
      close(unit=ipca,status='delete')
    end if
    ipca=freeunit()
    open(unit=ipca,file=fn(ii:jj),status='new')

 67 format('#',g10.4,99999(1x,g10.4))
 66 format(100000000(1x,g10.4))
!
    write(ipca,67) bla6(1:min(effsz1,sz2)) ! maximized measure (variance or ACF)
    do i=1,effsz1
      bla4(i) = sum(bla5(:,i))/(1.0*sz2)
    end do
    do i=1,sz2
      write(ipca,66) (bla5(i,j) - bla4(j),j=1,effsz1)      ! transformed and centered data
    end do
    close(unit=ipca)
!
!   if the user requested to keep only part of the data for clustering, data need to be transferred and settings adjusted
    if (reduced_dim_clustering.gt.0) then
      deallocate(cludata)
      allocate(cludata(reduced_dim_clustering,cstored))
      do i=1,reduced_dim_clustering
        cludata(i,1:sz2) = bla5(1:sz2,i)
      end do
      if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6).OR.(cdis_crit.eq.10)) then
        if (align_for_clustering.EQV..true.) then
          align_for_clustering = .false.
        end if
      end if
      sz1 = reduced_dim_clustering
      sz3 = reduced_dim_clustering
      cdis_crit = 7 ! simplest type of distance 
    end if
!
  end if
!
  deallocate(bla5)
!
#ifdef ENABLE_MPI
  call int2str(myrank,nod,re_aux(10))
  if (use_REMC.EQV..true.) then
    fn = 'N_'//nod(1:re_aux(10))//'_PRINCIPAL_COMPONENTS.evs'
  else
    fn = 'PRINCIPAL_COMPONENTS.evs'
  end if
#else
  fn = 'PRINCIPAL_COMPONENTS.evs'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ipca = freeunit()
    open(unit=ipca,file=fn(ii:jj),status='old',position='append')
    close(unit=ipca,status='delete')
  end if
  ipca=freeunit()
  open(unit=ipca,file=fn(ii:jj),status='new')
!
  write(ipca,67) bla2(1:min(effsz1,sz2))
  if (pcamode.le.3) then
    do i=1,effsz1
      write(ipca,66) (bla1(i,j),j=1,effsz1)
    end do
  else
    do i=1,effsz1
      write(ipca,66) (bla3(i,j),j=1,effsz1)
    end do
  end if
  close(unit=ipca)
!
  deallocate(bla7)
  deallocate(bla6)
  deallocate(bla1)
  deallocate(bla2)
  deallocate(bla4)
  deallocate(bla3)
!
end
!
#endif
!
!--------------------------------------------------------------------------------------------------------
!
! compatibility fxn to remove dimensions from the end
! this is almost exclusively for the case of mimicking arbitrary matrix input as an RMSD (no alignment) problem
!
subroutine reduce_cludimensionality()
!
  use clusters
  use iounit
!
  implicit none
!
  integer k
  RTYPE, ALLOCATABLE:: bla5(:,:)
!
  if (reduced_dim_clustering.le.0) return
!
! if the user requested to keep only part of the data for clustering, data need to be transferred
  allocate(bla5(calcsz,cstored))
  bla5(1:calcsz,1:cstored) = cludata(1:calcsz,1:cstored)
!
  if ((cdis_crit.eq.9).OR.((cdis_crit.eq.10).AND.(align_for_clustering.EQV..false.))) then
    if ((reduced_dim_clustering.gt.0).AND.(reduced_dim_clustering.lt.clstsz)) then
      deallocate(cludata)
      k = clstsz
      calcsz = reduced_dim_clustering*2
      clstsz = reduced_dim_clustering
      allocate(cludata(calcsz,cstored))
      cludata(1:clstsz,1:cstored) =  bla5(1:clstsz,1:cstored)
      cludata((clstsz+1):calcsz,1:cstored) = bla5((k+1):(k+clstsz),1:cstored)
      if (cdis_crit.eq.10) then
        cdofsbnds(2) = clstsz
        cdofsbnds(4) = clstsz
      end if
    else
      write(ilog,*) 'Warning. Request to reduce dimensionality to ',reduced_dim_clustering,' was ignored (out of range).'
    end if
  else if ((((cdis_crit.eq.5).OR.(cdis_crit.eq.6)).AND.(align_for_clustering.EQV..false.)).OR.&
 &         (cdis_crit.eq.7).OR.(cdis_crit.eq.1).OR.(cdis_crit.eq.3)) then
    if ((reduced_dim_clustering.gt.0).AND.(reduced_dim_clustering.lt.calcsz)) then
      deallocate(cludata)
      calcsz = reduced_dim_clustering
      clstsz = reduced_dim_clustering
      allocate(cludata(calcsz,cstored))
      cludata(1:calcsz,1:cstored) =  bla5(1:calcsz,1:cstored)
      if ((cdis_crit.eq.5).OR.(cdis_crit.eq.6)) then
        cdofsbnds(2) = calcsz
        cdofsbnds(4) = calcsz
      end if
    else
      write(ilog,*) 'Warning. Request to reduce dimensionality to ',reduced_dim_clustering,' was ignored (out of range).'
    end if
  else if (cdis_crit.eq.2) then
    if ((reduced_dim_clustering.gt.0).AND.(reduced_dim_clustering.lt.clstsz)) then
      deallocate(cludata)
      calcsz = 2*reduced_dim_clustering
      clstsz = reduced_dim_clustering
      allocate(cludata(calcsz,cstored))
      cludata(1:calcsz,1:cstored) =  bla5(1:calcsz,1:cstored)
    else
      write(ilog,*) 'Warning. Request to reduce dimensionality to ',reduced_dim_clustering,' was ignored (out of range).'
    end if
  else if (cdis_crit.eq.4) then
    if ((reduced_dim_clustering.gt.0).AND.(reduced_dim_clustering.lt.clstsz)) then
      deallocate(cludata)
      calcsz = 3*reduced_dim_clustering/2
      clstsz = reduced_dim_clustering
      allocate(cludata(calcsz,cstored))
      cludata(1:calcsz,1:cstored) =  bla5(1:calcsz,1:cstored)
    else
      write(ilog,*) 'Warning. Request to reduce dimensionality to ',reduced_dim_clustering,' was ignored (out of range).'
    end if
  else if (cdis_crit.eq.8) then
    if ((reduced_dim_clustering.gt.0).AND.(reduced_dim_clustering.lt.calcsz)) then
      deallocate(cludata)
      k = calcsz
      calcsz = reduced_dim_clustering
      clstsz = reduced_dim_clustering
      allocate(cludata(calcsz,cstored))
      cludata(1:calcsz,1:cstored) =  bla5(1:calcsz,1:cstored)
      bla5(1:k,1) = cl_imvec(1:k)
      deallocate(cl_imvec)
      allocate(cl_imvec(calcsz))
      cl_imvec(1:calcsz) = bla5(1:calcsz,1)
      bla5(1,1) = sum(cl_imvec(1:calcsz)**2)
      if (bla5(1,1).le.0.0) then ! no variance in data
        cl_imvec(1:calcsz) = 1.0
      else
        cl_imvec(1:calcsz) = sqrt((cl_imvec(1:calcsz)**2)/bla5(1,1))
      end if
    else
      write(ilog,*) 'Warning. Request to reduce dimensionality to ',reduced_dim_clustering,' was ignored (out of range).'
    end if

  end if
!
  deallocate(bla5)
!
end
!
!--------------------------------------------------------------------------------------------------
!
! the scripts are based more or less directly on work by VMD-afficionados:
! Axel Kohlmeyer (general), Andrew Dalke (sscache-stuff), Olaf Lenz (box stuff)
!
subroutine vmd_helper_for_clustering(it,nbasins)
!
  use iounit
  use mpistuff
  use clusters
  use molecule
  use sequen
  use aminos
  use system
  use mcsums
  use pdb
!
  implicit none
!
  integer i,j,ii,jj,ipdbh,t2,t3,t4,t5,k,imol,repcnt,freeunit
  character(MAXSTRLEN) fn,fn2
#ifdef ENABLE_MPI
  integer tslsh
  character(re_aux(10)) nod
#endif
  character(MAXSTRLEN) itn,itn2
  integer haveprotein,havenucleic,cup(2),cun(2),nbasins
  logical exists
  type(t_scluster) it(nbasins)
!
  haveprotein = 0
  havenucleic = 0
  cup(2) = 0
  cun(2) = 0
  do imol=1,nmol
    cun(1) = 0
    cup(1) = 0
    do i=rsmol(imol,1),rsmol(imol,2)
      if (seqpolty(i).eq.'P') then
        haveprotein = haveprotein + 1
        cup(1) = cup(1) + 1
      end if
      if (seqpolty(i).eq.'N') then
        havenucleic = havenucleic + 1
        cun(1) = cun(1) + 1
      end if
    end do
    if (cun(1).gt.cun(2)) cun(2) = cun(1)
    if (cup(1).gt.cup(2)) cup(2) = cup(1)
  end do
!
#ifdef ENABLE_MPI
  call int2str(myrank,nod,re_aux(10))
  if (use_REMC.EQV..true.) then
    t2 = bleng + 3 + re_aux(10)
    itn(1:t2) = 'N_'//nod(1:re_aux(10))//'_'//basename(1:bleng)
    fn = 'N_'//nod(1:re_aux(10))//'_STRUCT_CLUSTERING.vmd'
  else
    itn(1:bleng) = basename(1:bleng)
    t2 = bleng
    fn = 'STRUCT_CLUSTERING.vmd'
  end if
  t3 = bleng + 3 + re_aux(10)
  itn2(1:t3) = 'N_'//nod(1:re_aux(10))//'_'//basename(1:bleng)
  if (use_REMC.EQV..true.) then
    if ((re_aux(3).eq.1).AND.(pdb_analyze.EQV..true.)) then
      call strlims(fn,ii,jj)
      write(ilog,*) 'Warning. In structural clustering analysis on trajectories unscrambled via keyword &
 &FMCSC_TRACEFILE, it is not possible to obtain the auxiliary VMD script (',fn(ii:jj),'). It is &
 &recommended to perform this operation in two steps (a) transcribe trajectories, b) run clustering).'
      return
    end if
  end if
#else
  itn(1:bleng) = basename(1:bleng)
  itn2(1:bleng) = basename(1:bleng)
  t2 = bleng
  t3 = bleng
  fn = 'STRUCT_CLUSTERING.vmd'
#endif
!
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj),exist=exists)
  if(exists) then
    ipdbh = freeunit()
    open(unit=ipdbh,file=fn(ii:jj),status='old',position='append')
    close(unit=ipdbh,status='delete')
  end if
#ifdef ENABLE_MPI
  if (use_MPIAVG.EQV..true.) return ! never write vmd file for MPI averaging run
#endif
  ipdbh=freeunit()
  open(unit=ipdbh,file=fn(ii:jj),status='new')
!
  write(ipdbh,*) '## Credits for these scripts go to:'
  write(ipdbh,*) '## Axel Kohlmeyer'
  write(ipdbh,*) '## Andrew Dalke'
  write(ipdbh,*) '## Olaf Lenz'
  write(ipdbh,*)
  write(ipdbh,*) 'color Display Background white'
  write(ipdbh,*) 'material change opacity Transparent 0.7'
  write(ipdbh,*) 'display culling off'
  write(ipdbh,*) 'display depthcue on'
  write(ipdbh,*) 'display cuedensity 0.3'
  write(ipdbh,*) 'display rendermode GLSL'
  write(ipdbh,*) 'axes location off'
  write(ipdbh,*)
  write(ipdbh,*) 'set rep1 Licorice'
  write(ipdbh,*) 'set rep2 NewCartoon'
  write(ipdbh,*) 'set rep3 NewRibbons'  !Lines'
  write(ipdbh,*) 'set sel1 "all"'
  write(ipdbh,*) 'set sel2 "all"'
  write(ipdbh,*) 'set trjsm 0'
  write(ipdbh,*)
  if ((use_pdb_template.EQV..false.).OR.(pdb_analyze.EQV..false.)) then
    write(ipdbh,*) 'if { [glob -nocomplain ',itn2(1:t3),'_END.pdb ] \'
    write(ipdbh,*) '== "',itn2(1:t3),'_END.pdb" } then \'
    write(ipdbh,*) '{set ftemplate "',itn2(1:t3),'_END.pdb"} \'
    write(ipdbh,*) 'else {set ftemplate "',itn2(1:t3),'_START.pdb"}'
  else
    call strlims(pdbtmplfile,t4,t5)
    write(ipdbh,*) 'set ftemplate "',pdbtmplfile(t4:t5),'"'
  end if
  write(ipdbh,*) 'mol load pdb $ftemplate'
  if (pdb_analyze.EQV..false.) then
    if (xyzmode.eq.1) then
      k = floor(1.0*(nsim - nequil)/(1.0*xyzout)) + 1
      write(ipdbh,*) 'for {set i 1} {$i<=',k,'} {incr i 1} {'
      write(ipdbh,*) '  set fileName [format ',itn(1:t2),'_%05d.arc $i]'
      write(ipdbh,*) '  mol addfile $fileName waitfor all'
      write(ipdbh,*) '}'
    else if (xyzmode.eq.2) then
      if (pdb_writemode.eq.2) then
        write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.pdb waitfor all'
      else if (pdb_writemode.eq.1) then
        k = floor(1.0*(nsim - nequil)/(1.0*xyzout)) + 1
        write(ipdbh,*) 'for {set i 1} {$i<=',k,'} {incr i 1} {'
        write(ipdbh,*) '  set fileName [format ',itn(1:t2),'_%05d.pdb $i]'
        write(ipdbh,*) '  mol addfile $fileName waitfor all'
        write(ipdbh,*) '}'
      end if
    else if (xyzmode.eq.3) then
      write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.dcd waitfor all'
    else if (xyzmode.eq.4) then
      write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.xtc waitfor all'
    else if (xyzmode.eq.5) then
      write(ipdbh,*) 'mol addfile ',itn(1:t2),'_traj.nc waitfor all'
    end if
  else
    do i=nequil+1,nsim
      if (mod(i,cstorecalc).eq.0) exit
    end do
    i = min(i,nsim)
    if (pdb_fileformat.eq.2) then
      write(ipdbh,*) 'for {set i ',i+pdb_format(4),'} {$i<=',nsim+pdb_format(4),'} \'
      write(ipdbh,*) '{incr i ',cstorecalc,'} {'
      if ((pdb_format(1).gt.0).AND.(pdb_format(2).gt.0)) then
        write(ipdbh,*) '  set fileName [format ',pdb_pref(1:pdb_format(1)),'%05d',pdb_suff(1:pdb_format(2)),' $i]'
      else if (pdb_format(1).gt.0) then
        write(ipdbh,*) '  set fileName [format ',pdb_pref(1:pdb_format(1)),'%05d $i]'
      else if (pdb_format(2).gt.0) then
        write(ipdbh,*) '  set fileName [format %05d',pdb_suff(1:pdb_format(2)),' $i]'
      else
        write(ipdbh,*) '  set fileName [format %05d $i]'
      end if
      write(ipdbh,*) '  mol addfile $fileName waitfor all'
      write(ipdbh,*) '}'
    else if (((pdb_fileformat.ge.3).AND.(pdb_fileformat.le.5)).OR.(pdb_fileformat.eq.1)) then
      if (pdb_fileformat.eq.1) then
        call strlims(pdbinfile,t4,t5)
        fn2(t4:t5) = pdbinfile(t4:t5)
      else if (pdb_fileformat.eq.3) then
        call strlims(xtcinfile,t4,t5)
        fn2(t4:t5) = xtcinfile(t4:t5)
      else if (pdb_fileformat.eq.4) then
        call strlims(dcdinfile,t4,t5)
        fn2(t4:t5) = dcdinfile(t4:t5)
      else if (pdb_fileformat.eq.5) then
        call strlims(netcdfinfile,t4,t5)
        fn2(t4:t5) = netcdfinfile(t4:t5)
      end if
#ifdef ENABLE_MPI
      do j=t5,t4,-1
        if (fn2(j:j).eq.SLASHCHAR) exit
      end do
      tslsh = j
      call int2str(myrank,nod,re_aux(10))
      if ((tslsh.ge.t4).AND.(tslsh.lt.t5)) then
        fn =  fn2(t4:tslsh)//'N_'//nod(1:re_aux(10))//'_'//fn2((tslsh+1):t5)
      else
        fn =  'N_'//nod(1:re_aux(10))//'_'//fn2(t4:t5)
      end if
      call strlims(fn,t4,t5)
#else
      fn = fn2
#endif
      write(ipdbh,*) 'mol addfile ',fn(t4:t5),' \'
      write(ipdbh,*) 'step ',cstorecalc,' first ',i-1,' \'
      write(ipdbh,*) 'last ',nsim-1,' waitfor all'
    end if
  end if
  write(ipdbh,*) 'mol modcolor 0 0 Type'
  write(ipdbh,*) 'mol modstyle 0 0 $rep1'
  write(ipdbh,*) 'mol modmaterial 0 0 Glossy'
  write(ipdbh,*) 'mol modselect 0 0 $sel1'
  write(ipdbh,*) 'mol addrep 0'
  write(ipdbh,*) 'mol modstyle 1 0 $rep2'
  write(ipdbh,*) 'mol modcolor 1 0 Index'
  write(ipdbh,*) 'mol modselect 1 0 $sel2'
  write(ipdbh,*) 'mol modmaterial 1 0 Transparent'
  write(ipdbh,*) 'mol showrep 0 0 off'
  write(ipdbh,*) 'mol showrep 0 1 off'
  repcnt = 2
  do i=1,nbasins
    if (it(i)%nmbrs.lt.50) cycle
    write(ipdbh,*) 'mol addrep 0'
    write(ipdbh,*) 'mol modstyle ',repcnt,' 0 $rep3'
    write(ipdbh,*) 'mol modmaterial ',repcnt,'  0 Glossy'
    write(ipdbh,*) 'mol modselect ',repcnt,'  0 $sel1'
    write(ipdbh,*) 'mol modcolor ',repcnt,' 0 ColorID ',mod((repcnt-2)/2,33)
    write(ipdbh,*) 'mol drawframes 0 ',repcnt,' " \'
    if (pdb_analyze.EQV..false.) then
      do j=1,it(i)%nmbrs
        if (mod(cstorecalc*it(i)%snaps(j),xyzout).eq.0) then
          write(ipdbh,*) (cstorecalc*it(i)%snaps(j))/xyzout,'\'
        end if
      end do
    else
      if (select_frames.EQV..true.) then
        do j=1,it(i)%nmbrs
          write(ipdbh,*) framelst(it(i)%snaps(j)),'\'
        end do
      else
        do j=1,it(i)%nmbrs
          write(ipdbh,*) it(i)%snaps(j),'\'
        end do
      end if
    end if
    write(ipdbh,*) '"'
    write(ipdbh,*)
    write(ipdbh,*) 'mol addrep 0'
    write(ipdbh,*) 'mol modstyle ',repcnt+1,' 0 $rep1'
    write(ipdbh,*) 'mol modmaterial ',repcnt+1,'  0 Glossy'
    write(ipdbh,*) 'mol modselect ',repcnt+1,'  0 $sel1'
    write(ipdbh,*) 'mol modcolor ',repcnt+1,' 0 ColorID ',mod((repcnt-2)/2,33)
    if (pdb_analyze.EQV..false.) then
      if (mod(cstorecalc*it(i)%center,xyzout).eq.0) then
        write(ipdbh,*) 'mol drawframes 0 ',repcnt+1,' "',(cstorecalc*it(i)%center)/xyzout,'"'
      end if
    else
      if (select_frames.EQV..true.) then
        write(ipdbh,*) 'mol drawframes 0 ',repcnt+1,' "',framelst(it(i)%center),'"'
      else
        write(ipdbh,*) 'mol drawframes 0 ',repcnt+1,' "',it(i)%center,'"'
      end if
    end if
    write(ipdbh,*)
    write(ipdbh,*) 'mol showrep 0 ',repcnt,' off'
    write(ipdbh,*) 'mol showrep 0 ',repcnt+1,' off'
    repcnt = repcnt + 2
  end do
!
  write(ipdbh,*)
  write(ipdbh,*) 'set nreps [molinfo 0 get numreps]'
  write(ipdbh,*) 'set nfrms [molinfo 0 get numframes]'
  write(ipdbh,*) 
  write(ipdbh,*) 'for {set i 0} {$i < $nreps} {incr i} {'
  write(ipdbh,*) '  mol selupdate $i 0 on'
  write(ipdbh,*) '  mol smoothrep 0 $i $trjsm'
  write(ipdbh,*) '}'
  write(ipdbh,*)
  write(ipdbh,*) 'proc start_sscache {{molid 0}} {'
  write(ipdbh,*) '  global sscache_data vmd_frame'
  write(ipdbh,*) '  trace variable vmd_frame($molid) w sscache'
  write(ipdbh,*) '  return'
  write(ipdbh,*) '}'
  write(ipdbh,*) 'proc stop_sscache {{molid 0}} {'
  write(ipdbh,*) '  global vmd_frame'
  write(ipdbh,*) '  trace vdelete vmd_frame($molid) w sscache'
  write(ipdbh,*) '  return'
  write(ipdbh,*) '}'
  write(ipdbh,*) 'proc reset_sscache {} {'
  write(ipdbh,*) '  global sscache_data'
  write(ipdbh,*) '  if [info exists sscache_data] {'
  write(ipdbh,*) '    unset sscache_data'
  write(ipdbh,*) '  }'
  write(ipdbh,*) '  return'
  write(ipdbh,*) '}'
  write(ipdbh,*) 'proc sscache {name index op} {'
  write(ipdbh,*) '  global sscache_data'
  write(ipdbh,*) '  set sel [atomselect $index "protein name CA"]'
  write(ipdbh,*) '  set frame [molinfo $index get frame]'
  write(ipdbh,*) '  if [info exists sscache_data($index,$frame)] {'
  write(ipdbh,*) '    $sel set structure $sscache_data($index,$frame)'
  write(ipdbh,*) '  return'
  write(ipdbh,*) '  }'
  write(ipdbh,*) '  vmd_calculate_structure $index'
  write(ipdbh,*) '  set sscache_data($index,$frame) [$sel get structure]'
  write(ipdbh,*) '  return'
  write(ipdbh,*) '}'
  write(ipdbh,*)
  if ((havenucleic.eq.0).AND.(haveprotein.eq.0)) then
    write(ipdbh,*) 'proc align_it {{sel "all"} {rfmol 0} {reffr 0} {molid 0}} {'
  else if (havenucleic.eq.0) then
    write(ipdbh,*) 'proc align_it {{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
  else if (haveprotein.eq.0) then
    write(ipdbh,*) 'proc align_it {{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
  else if (ceiling(1.4*cun(2)).ge.cup(2)) then
    write(ipdbh,*) 'proc align_it {{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
  else
    write(ipdbh,*) 'proc align_it {{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
  end if
  write(ipdbh,*) '  set alignref [atomselect $rfmol $sel frame $reffr]'
  write(ipdbh,*) '  set nf [molinfo $molid get numframes]'
  write(ipdbh,*) '  for {set i 0} {$i < $nf} {incr i} {'
  write(ipdbh,*) '    set match [atomselect $molid $sel frame $i]'
  write(ipdbh,*) '    set tm [measure fit $match $alignref]'
  write(ipdbh,*) '    set alls [atomselect $molid "all" frame $i]'
  write(ipdbh,*) '    $alls move $tm'
  write(ipdbh,*) '  }'
  write(ipdbh,*) '}'
  write(ipdbh,*)
  write(ipdbh,*) 'proc rmsd_it {{filen \'
  write(ipdbh,*) itn(1:t2),'_RMSD.dat} \'
  if ((havenucleic.eq.0).AND.(haveprotein.eq.0)) then
    write(ipdbh,*) '{sel "all"} {rfmol 0} {reffr 0} {molid 0}} {'
  else if (havenucleic.eq.0) then
    write(ipdbh,*) '{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
  else if (haveprotein.eq.0) then
    write(ipdbh,*) '{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
  else if (ceiling(1.4*cun(2)).ge.cup(2)) then
    write(ipdbh,*) '{sel "nucleic"} {rfmol 0} {reffr 0} {molid 0}} {'
  else
    write(ipdbh,*) '{sel "protein"} {rfmol 0} {reffr 0} {molid 0}} {'
  end if
  write(ipdbh,*) '  set fi [open $filen w]'
  write(ipdbh,*) '  set alignref [atomselect $rfmol $sel frame $reffr]'
  write(ipdbh,*) '  set refall [atomselect $rfmol "all" frame $reffr]'
  write(ipdbh,*) '  set nf [molinfo $molid get numframes]'
  write(ipdbh,*) '  for {set i 1} {$i < $nf} {incr i} {'
  write(ipdbh,*) '    set match [atomselect $molid $sel frame $i]'
  write(ipdbh,*) '    set tm [measure fit $match $alignref]'
  write(ipdbh,*) '    set alls [atomselect $molid "all" frame $i]'
  write(ipdbh,*) '    $alls move $tm'
  write(ipdbh,*) '    set match [atomselect $molid $sel frame $i]'
  write(ipdbh,*) '    set irmsd [measure rmsd $match $alignref]'
  write(ipdbh,*) '    set alls [atomselect $molid "all" frame $i]'
  write(ipdbh,*) '    set irmsd2 [measure rmsd $alls $refall]'
  write(ipdbh,*) '    puts $fi "$i  $irmsd  $irmsd2"'
  write(ipdbh,*) '  }'
  write(ipdbh,*) '  close $fi'
  write(ipdbh,*) '}'
  write(ipdbh,*)
  write(ipdbh,*)
  write(ipdbh,*) 'start_sscache 0'
  write(ipdbh,*) 'align_it'
  write(ipdbh,*)
!
end
!
!-----------------------------------------------------------------------------------
!
