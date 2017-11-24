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
!
subroutine gen2Demproj()
!
  use atoms
  use iounit
!
  implicit none
!
end
!
!----------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
subroutine read_netcdf3Dmap()
!
  use ems
  use iounit
  use netcdf
  use system
  use units
  use energies
  use atoms
!
  implicit none
!
  integer iemmap,t1,t2,dimlen,i,nspats(3),ret,fndds(5),lenner,xtyper,k
  integer iomess,aone,atwo,freeunit
  character(MAXSTRLEN) trystr,ucstr
  RTYPE scalefac,bgd
  logical exists
!
  aone = 1
  atwo = 2
  call strlims(emmapfile,t1,t2) 
  inquire (file=emmapfile(t1:t2),exist=exists)
! 
  if (exists.EQV..false.) then
    write(ilog,*) 'Fatal. Cannot open spec.d input file for reading &
 &from NetCDF (',emmapfile(t1:t2),') in read_netcdf3Dmap().'
    call fexit()
  end if
!
  if ((empotprop.gt.2).AND.(use_EMICRO.EQV..true.)) then
    write(ilog,*) 'Warning. The chosen property (',empotprop,') is currently not available &
 &in runs with the density restraint potential. Disabling term.'
    use_EMICRO = .false.
    scale_EMICRO = 0.0
    return
  end if
!
 44 format('Warning. Ambigous dimensions in NetCDF file for input density map. Encountered ',a,' twice but keeping&
 & only first.')
!
! open
  iemmap = freeunit()
  call strlims(emmapfile,t1,t2)
  call check_fileionetcdf( nf90_open(emmapfile(t1:t2), NF90_NOWRITE, iemmap) )
!
! find the necessary dimensions: X Y Z required
  fndds(:) = 0
  nspats(3) = 0
  emmapsz(:) = 0
  do i=1,NF90_MAX_DIMS
    ret = nf90_inquire_dimension(iemmap,i,trystr,dimlen)
    if (ret.eq.NF90_NOERR) then
      ucstr(1:15) = trystr(1:15)
      call toupper(ucstr(1:15))
      if (ucstr(1:1).eq.'Z') then
        if (fndds(1).eq.1) then 
          write(ilog,44) ucstr(1:1)
        else
          emmapsz(3) = dimlen
          fndds(1) = 1
          emnetcdf_ids(1) = i
        end if
      else if (ucstr(1:1).eq.'X') then
        if (fndds(2).eq.1) then 
          write(ilog,44) ucstr(1:1)
        else
          emmapsz(1) = dimlen
          fndds(2) = 1
          emnetcdf_ids(2) = i
        end if
      else if (ucstr(1:1).eq.'Y') then
        if (fndds(3).eq.1) then 
          write(ilog,44) ucstr(1:1)
        else
          emmapsz(2) = dimlen
          fndds(3) = 1
          emnetcdf_ids(3) = i
        end if
      else if (ucstr(1:7).eq.'spatial') then
        if (fndds(4).eq.1) then 
          write(ilog,44) ucstr(1:7)
        else
          fndds(4) = 1
          emnetcdf_ids(4) = i
        end if
      end if

    else if (ret.eq.NF90_EBADDIM) then
!     do nothing
    else ! get us out of here
     call check_fileionetcdf( nf90_inquire_dimension(iemmap,i,trystr,dimlen) )
    end if
  end do
!
  if ((emmapsz(3).lt.1).OR.(emmapsz(2).lt.1).OR.(emmapsz(1).lt.1)) then
    write(ilog,*) 'Fatal. NetCDF-file (',emmapfile(t1:t2),') has at least one &
 &ill-defined spatial dimension.'
    call fexit()
  end if
  emingrid%dim(:) = emmapsz(1:3)
!

! now find the necessary variables: grid spacings and densities are required, origin is optional
  fndds(:) = 0
  ret = nf90_inquire_attribute(iemmap,NF90_GLOBAL,"xyz_step",xtyper,lenner)
  if (ret.eq.NF90_NOERR) then
    if (((xtyper.eq.NF90_FLOAT).OR.(xtyper.eq.NF90_DOUBLE)).AND.(lenner.eq.3)) then
      call check_fileionetcdf(  nf90_get_att(iemmap,NF90_GLOBAL,"xyz_step",emingrid%deltas(:)) )
      fndds(1) = 1
    else if (xtyper.eq.NF90_CHAR) then
      call check_fileionetcdf(  nf90_get_att(iemmap,NF90_GLOBAL,"xyz_step",trystr(1:lenner)) )
      read(trystr(1:lenner),*,iostat=iomess) emingrid%deltas(1),emingrid%deltas(2),emingrid%deltas(3)
      if (iomess.eq.0) fndds(1) = 1
    end if
    if (fndds(1).ne.1) then
      write(ilog,*) 'Found global attribute "xyz_step", but could not parse it. Will &
 &look for variable "deltas" in input density NetCDF-file (',emmapfile(t1:t2),') instead.'
    else
      write(ilog,*) 'Assuming spatial units of Angstroms for global attribute "xyz_step"&
 & in input density NetCDF-file (',emmapfile(t1:t2),').'
      scalefac = 1.0
    end if
  end if
!
  if (fndds(1).ne.1) then
    ret = nf90_inq_varid(iemmap,"deltas",emnetcdf_ids(5))
    if (ret.eq.NF90_NOERR) then
      fndds(1) = 1
    else if (ret.eq.NF90_ENOTVAR) then
      write(ilog,*) 'Fatal. Variable "deltas" not found in input density NetCDF-file (',&
 &emmapfile(t1:t2),'). Use NetCDFs ncdump utility to check file.'
      call fexit()
    else
      call check_fileionetcdf( nf90_inq_varid(iemmap,"deltas",emnetcdf_ids(5)) )
    end if
!
    ret = nf90_inquire_attribute(iemmap,emnetcdf_ids(5),"units",xtyper,lenner)
    if ((ret.eq.NF90_NOERR).AND.(xtyper.eq.NF90_CHAR)) then
      call check_fileionetcdf(  nf90_get_att(iemmap,emnetcdf_ids(5),"units",trystr(1:lenner)) )
      emspatunitstr(1:lenner) = trystr(1:lenner)
    else 
      write(ilog,*) 'Fatal. Variable "deltas" in input density NetCDF-file (',&
 &emmapfile(t1:t2),') has no interpretable, associated unit.'
      call fexit()
    end if
!
    call tolower(emspatunitstr(1:lenner))
    if ((emspatunitstr(1:lenner).eq.'angstrom').OR.(emspatunitstr(1:lenner).eq.'a').OR.&
 &      (emspatunitstr(1:lenner).eq.'ang')) then
      scalefac = 1.0
    else if ((emspatunitstr(1:lenner).eq.'nanometer').OR.(emspatunitstr(1:lenner).eq.'nm')) then
      scalefac = 0.1
    else if ((emspatunitstr(1:lenner).eq.'picometer').OR.(emspatunitstr(1:lenner).eq.'pm')) then
      scalefac = 100
    else if ((emspatunitstr(1:lenner).eq.'micrometer').OR.(emspatunitstr(1:lenner).eq.'um').OR.&
 &           (emspatunitstr(1:lenner).eq.'microns')) then
      scalefac = 0.0001
    else if ((emspatunitstr(1:lenner).eq.'meter').OR.(emspatunitstr(1:lenner).eq.'m')) then
      scalefac = 1.0e-10
    else 
      write(ilog,*) 'Fatal. Variable "deltas" in input density NetCDF-file (',&
 &emmapfile(t1:t2),') has an unparsable unit (use Angstrom, nm, um, m, or pm).'
      call fexit()
    end if
    call check_fileionetcdf( nf90_get_var(iemmap, emnetcdf_ids(5), emingrid%deltas(:)) )
  end if
!
 47 format('Warning. In dimension ',i2,', the provided input density map (',a,') does not line up with &
 &the system dimensions (',g18.12,' vs. ',g18.12,' Angstrom).')
  if (bnd_shape.eq.1) then
    do k=1,3
      if (abs(emingrid%dim(k)*scalefac*emingrid%deltas(k) - bnd_params(k)).gt.1.0e-7) then
        write(ilog,47) k,emmapfile(t1:t2),bnd_params(k),emingrid%dim(k)*scalefac*emingrid%deltas(k)
      end if
    end do
  end if
!
  scalefac = -1.0
  ret = nf90_inquire_attribute(iemmap,NF90_GLOBAL,"xyz_origin",xtyper,lenner)
  if (ret.eq.NF90_NOERR) then
    if (((xtyper.eq.NF90_FLOAT).OR.(xtyper.eq.NF90_DOUBLE)).AND.(lenner.eq.3)) then
      call check_fileionetcdf(  nf90_get_att(iemmap,NF90_GLOBAL,"xyz_origin",emingrid%origin(:)) )
      fndds(3) = 1
    else if (xtyper.eq.NF90_CHAR) then
      call check_fileionetcdf(  nf90_get_att(iemmap,NF90_GLOBAL,"xyz_origin",trystr(1:lenner)) )
      read(trystr(1:lenner),*,iostat=iomess) emingrid%origin(1),emingrid%origin(2),emingrid%origin(3)
      if (iomess.eq.0) fndds(3) = 1
    end if
    if (fndds(3).ne.1) then
      write(ilog,*) 'Found global attribute "xyz_origin", but could not parse it. Will &
 &look for variable "origin" in input density NetCDF-file (',emmapfile(t1:t2),') instead.'
    else
      write(ilog,*) 'Assuming spatial units of Angstroms for global attribute "xyz_origin"&
 & in input density NetCDF-file (',emmapfile(t1:t2),').'
      scalefac = 1.0
    end if
  end if
!
  if (fndds(3).ne.1) then
    ret = nf90_inq_varid(iemmap,"origin",emnetcdf_ids(5))
    if (ret.eq.NF90_NOERR) then
      fndds(3) = 1
    else if (ret.eq.NF90_ENOTVAR) then
      write(ilog,*) 'Warning. Variable "origin" not found in input density NetCDF-file (',&
 &emmapfile(t1:t2),'). Assuming system origin.'
      if (bnd_shape.eq.1) then
        emingrid%origin(:) = bnd_params(4:6)
      else if (bnd_shape.eq.2) then
        emingrid%origin(:) = bnd_params(1:3)
      end if
    else
      call check_fileionetcdf( nf90_inq_varid(iemmap,"origin",emnetcdf_ids(5)) )
    end if
!
    ret = nf90_inquire_attribute(iemmap,emnetcdf_ids(5),"units",xtyper,lenner)
    if ((ret.eq.NF90_NOERR).AND.(xtyper.eq.NF90_CHAR)) then
      call check_fileionetcdf(  nf90_get_att(iemmap,emnetcdf_ids(5),"units",trystr(1:lenner)) )
      emoriunitstr(1:lenner) = trystr(1:lenner)
    else
      if (fndds(3).eq.1) then 
        write(ilog,*) 'Warning. Variable "origin" in input density NetCDF-file (',&
 &emmapfile(t1:t2),') has no interpretable, associated unit. Assuming Angstrom.'
      end if
      scalefac = 1.0
    end if
!
    if (scalefac.lt.0.0) then
      call tolower(emoriunitstr(1:lenner))
      if ((emoriunitstr(1:lenner).eq.'angstrom').OR.(emoriunitstr(1:lenner).eq.'a').OR.&
 &        (emoriunitstr(1:lenner).eq.'ang')) then
        scalefac = 1.0
      else if ((emoriunitstr(1:lenner).eq.'nanometer').OR.(emoriunitstr(1:lenner).eq.'nm')) then
          scalefac = 0.1
      else if ((emoriunitstr(1:lenner).eq.'picometer').OR.(emoriunitstr(1:lenner).eq.'pm')) then
        scalefac = 100
      else if ((emoriunitstr(1:lenner).eq.'micrometer').OR.(emoriunitstr(1:lenner).eq.'um').OR.&
 &             (emoriunitstr(1:lenner).eq.'microns')) then
        scalefac = 0.0001
      else if ((emoriunitstr(1:lenner).eq.'meter').OR.(emoriunitstr(1:lenner).eq.'m')) then
        scalefac = 1.0e-10
      else 
        write(ilog,*) 'Warning. Variable "origin" in input density NetCDF-file (',&
 &emmapfile(t1:t2),') has an unparsable unit (use Angstrom, nm, um, m, or pm). Assuming Angstrom.'
        scalefac = 1.0
      end if
      call check_fileionetcdf( nf90_get_var(iemmap, emnetcdf_ids(5), emingrid%origin(:)) )
    end if
  end if
  if (scalefac.gt.0.0) emingrid%origin(:) = scalefac*emingrid%origin
!
 57 format('Warning. In dimension ',i2,', the provided input density map (',a,') does not line up with &
 &the system origin (',g18.12,' vs. ',g18.12,' Angstrom).')
  if (bnd_shape.eq.1) then
    do k=1,3
      if (abs(emingrid%origin(k) - bnd_params(k+3)).gt.1.0e-7) then
        write(ilog,57) k,emmapfile(t1:t2),bnd_params(k+3),scalefac*emingrid%origin(k)
      end if
    end do
  end if
!
  ret = nf90_inq_varid(iemmap,"densities",emnetcdf_ids(11))
  if (ret.eq.NF90_NOERR) then
    fndds(2) = 1
  else if (ret.eq.NF90_ENOTVAR) then
    write(ilog,*) 'Fatal. Variable "densities" not found in input density NetCDF-file (',&
 &emmapfile(t1:t2),'). Use NetCDFs ncdump utility to check file.'
    call fexit()
  else
    call check_fileionetcdf( nf90_inq_varid(iemmap,"densities",emnetcdf_ids(11)) )
  end if
!
! see if we can find units
  ret = nf90_inquire_attribute(iemmap,emnetcdf_ids(11),"units",xtyper,lenner)
  if ((ret.eq.NF90_NOERR).AND.(xtyper.eq.NF90_CHAR)) then
    call check_fileionetcdf(  nf90_get_att(iemmap,emnetcdf_ids(11),"units",trystr(1:lenner)) )
    emdensunitstr(1:lenner) = trystr(1:lenner)
  else 
    write(ilog,*) 'Warning. Variable "densities" in input density NetCDF-file (',&
 &emmapfile(t1:t2),') has no interpretable, associated unit. Assuming CAMPARI convention.'
  end if
!
! this block is useful for increasing the size of a map to specific input (note initialization)
!
!    if (bnd_shape.eq.1) then
!      emgrid%dim(:) = ceiling(bnd_params(1:3)/emingrid%deltas(:))
!      emgrid%origin(:) = bnd_params(4:6)
!    else if (bnd_shape.eq.2) then
!      emgrid%dim(:) = ceiling(3*bnd_params(4)/emingrid%deltas(:))
!      emgrid%origin(:) = bnd_params(1:3) - embuf*bnd_params(4)
!    end if
!  
!  allocate(emingrid%optdens(max(emgrid%dim(1),emingrid%dim(1)),max(emgrid%dim(2),emingrid%dim(2)),&
! &                          max(emgrid%dim(3),emingrid%dim(3))))
!  emingrid%optdens(:,:,:) = 0.0
!  call check_fileionetcdf( nf90_get_var(iemmap, emnetcdf_ids(11), &
! &emingrid%optdens(1:emingrid%dim(1),1:emingrid%dim(2),1:emingrid%dim(3))) )
!    emingrid%dim(:) = emgrid%dim(:)
!
  allocate(emingrid%optdens(emingrid%dim(1),emingrid%dim(2),emingrid%dim(3)))
  call check_fileionetcdf( nf90_get_var(iemmap, emnetcdf_ids(11), emingrid%optdens(:,:,:)) )
!
  if ((fdlts(1).gt.0.0).AND.(fdlts(2).gt.0.0).AND.(fdlts(3).gt.0.0)) then
    call filter_emmap(emingrid,fdlts)
  else if ((fdlts(1).eq.0.0).AND.(fdlts(2).eq.0.0).AND.(fdlts(3).eq.0.0)) then
!   do nothing
  else
    write(ilog,*) 'Warning. Received uninterpretable input for reducing resolution of input density NetCDF-file (',&
 &emmapfile(t1:t2),'). Request in FMCSC_EMREDUCE is ignored.'
  end if
!
  allocate(emingrid%dens(emingrid%dim(1),emingrid%dim(2),emingrid%dim(3)))
!
  call check_fileionetcdf( nf90_close(iemmap) )
!
  call allocate_ems(aone)
!
  call scale_emmap(emthreshdensity,bgd)
  call choose_simplification()
  call precompute_diff_emmap()
!
  call prt_emmap(emingrid,atwo)
!
end
!
#endif
!
!------------------------------------------------------------------------------------------------
!
subroutine  filter_emmap(thegrid,targetres)
!
  use ems
  use iounit
  use atoms
  use units
!
  implicit none
!
  type(t_emgrid) thegrid
  RTYPE targetres(3),fracx,fracy,fracz,frac
  RTYPE cur(2,3),vxszcorr
  integer xi,yi,zi,xj,yj,zj,tgsz(3),afive
!
  afive = 5
  tgsz = nint(thegrid%dim(:)*thegrid%deltas(:)/targetres(:))
  targetres(:) = thegrid%dim(:)*thegrid%deltas(:)/(1.0*tgsz(:))
  allocate(thegrid%dens(tgsz(1),tgsz(2),tgsz(3)))
  thegrid%dens(:,:,:) = 0.0
!
  cur(:,:) = 0.0
  xj = 1
!
  vxszcorr = thegrid%deltas(1)*thegrid%deltas(2)*thegrid%deltas(3)/(targetres(1)*targetres(2)*targetres(3))
  do xi=1,thegrid%dim(1)
    cur(1,1) = thegrid%origin(1)+xi*thegrid%deltas(1)
    cur(2,1) = thegrid%origin(1)+xj*targetres(1)
    if (cur(1,1).gt.cur(2,1)) then
      xj = min(xj+1,tgsz(1))
      cur(2,1) = thegrid%origin(1)+xj*targetres(1)
    end if
    fracx = min(1.0,(cur(1,1) - cur(2,1) + targetres(1))/thegrid%deltas(1))
    yj = 1
    do yi=1,thegrid%dim(2)
      cur(1,2) = thegrid%origin(2)+yi*thegrid%deltas(2)
      cur(2,2) = thegrid%origin(2)+yj*targetres(2)
      if (cur(1,2).gt.cur(2,2)) then
        yj = min(yj+1,tgsz(2))
        cur(2,2) = thegrid%origin(2)+yj*targetres(2)
      end if
      fracy = min(1.0,(cur(1,2) - cur(2,2) + targetres(2))/thegrid%deltas(2))
      zj = 1
      do zi=1,thegrid%dim(3)
        cur(1,3) = thegrid%origin(3)+zi*thegrid%deltas(3)
        cur(2,3) = thegrid%origin(3)+zj*targetres(3)
        if (cur(1,3).gt.cur(2,3)) then
          zj = min(zj+1,tgsz(3))
          cur(2,3) = thegrid%origin(3)+zj*targetres(3)
        end if
        fracz = min(1.0,(cur(1,3) - cur(2,3) + targetres(3))/thegrid%deltas(3))
        frac = fracx*fracy*fracz
        thegrid%dens(xj,yj,zj) =  thegrid%dens(xj,yj,zj) + frac*thegrid%optdens(xi,yi,zi)
        if (xj.gt.1) then
          frac = (1.0-fracx)*fracy*fracz
          thegrid%dens(xj-1,yj,zj) =  thegrid%dens(xj-1,yj,zj) + frac*thegrid%optdens(xi,yi,zi)
        end if
        if (yj.gt.1) then
          frac = (1.0-fracy)*fracx*fracz
          thegrid%dens(xj,yj-1,zj) =  thegrid%dens(xj,yj-1,zj) + frac*thegrid%optdens(xi,yi,zi)
        end if
        if (zj.gt.1) then
          frac = (1.0-fracz)*fracy*fracx
          thegrid%dens(xj,yj,zj-1) =  thegrid%dens(xj,yj,zj-1) + frac*thegrid%optdens(xi,yi,zi)
        end if
        if ((xj.gt.1).AND.(yj.gt.1)) then
          frac = (1.0-fracx)*(1.0-fracy)*fracz
          thegrid%dens(xj-1,yj-1,zj) =  thegrid%dens(xj-1,yj-1,zj) + frac*thegrid%optdens(xi,yi,zi)
        end if
        if ((xj.gt.1).AND.(zj.gt.1)) then
          frac = (1.0-fracx)*(1.0-fracz)*fracy
          thegrid%dens(xj-1,yj,zj-1) =  thegrid%dens(xj-1,yj,zj-1) + frac*thegrid%optdens(xi,yi,zi)
        end if
        if ((zj.gt.1).AND.(yj.gt.1)) then
          frac = (1.0-fracz)*(1.0-fracy)*fracx
          thegrid%dens(xj,yj-1,zj-1) =  thegrid%dens(xj,yj-1,zj-1) + frac*thegrid%optdens(xi,yi,zi)
        end if
        if ((xj.gt.1).AND.(yj.gt.1).AND.(zj.gt.1)) then
          frac = (1.0-fracz)*(1.0-fracy)*(1.0-fracx)
          thegrid%dens(xj-1,yj-1,zj-1) =  thegrid%dens(xj-1,yj-1,zj-1) + frac*thegrid%optdens(xi,yi,zi)
        end if
      end do
    end do
  end do
!
  thegrid%dens(:,:,:) = vxszcorr*thegrid%dens(:,:,:) 
  deallocate(thegrid%optdens)
  allocate(thegrid%optdens(tgsz(1),tgsz(2),tgsz(3)))
  thegrid%dim(:) = tgsz(:)
  thegrid%deltas(:) = targetres(:)
  thegrid%optdens(1:tgsz(1),1:tgsz(2),1:tgsz(3)) = thegrid%dens(1:tgsz(1),1:tgsz(2),1:tgsz(3))
  deallocate(thegrid%dens)
!
end
!
!------------------------------------------------------------------------------------------------
!
! modifies only emingrd%dens and emcoarsegrid%dens, does not require emgrid at all
! 
subroutine scale_emmap(threshold,bgd)
!
  use ems
  use iounit
  use atoms
  use units
  use params
!
  implicit none
!
  integer nrbins,xi,yi,zi,k,t1,t2,xj,yj,zj
  integer, ALLOCATABLE:: histd(:)
  RTYPE shfv,minde,maxde,emsf2,totm,tots,totv,voli,bgd,threshold,thresher,emfl,truncer
!
! get range of input data
  minde = minval(emingrid%optdens)
  maxde = maxval(emingrid%optdens)
!  
! shift such that minimum value equals physical bg density
  shfv = minde - embgdensity
  emingrid%dens(:,:,:) = emingrid%optdens(:,:,:) - shfv
  threshold = threshold - shfv
  maxde = maxde - shfv
  minde = minde - shfv
!
! bin the densities unless input provided
  if (emoptbg.eq.HUGE(emoptbg)) then
    nrbins = max(10,emingrid%dim(1)*emingrid%dim(2)*emingrid%dim(3)/20)
    allocate(histd(nrbins))
    histd(:) = 0
    do xi=1,emingrid%dim(1)
      do yi=1,emingrid%dim(2)
        do zi=1,emingrid%dim(3)
          k = min(nrbins,max(1,floor(nrbins*((emingrid%dens(xi,yi,zi)-minde)/(maxde-minde))) + 1))
          histd(k) = histd(k) + 1
        end do
      end do
    end do
    xi = 0
    do k=1,nrbins
      if (histd(k).gt.xi) then
        xi = histd(k)
        yi = k
      end if
    end do 
    bgd = minde + (yi-0.5)*(maxde-minde)/nrbins
    emoptbg = bgd + shfv
    deallocate(histd)
  else
    bgd = emoptbg - shfv
  end if
!
! flatten eventually (after histogram construction so as to not obscure solvent peak identification
  emfl = emflatval - shfv
  if ((emfl.gt.minde).AND.(emfl.le.maxde)) then
    if (emfl.le.threshold) then
      write(ilog,*) 'Warning. Requested flattening of input density map is inconsistent with chosen &
 &threshold value (',threshold+shfv,'). No flattening occurred.'
    else
      do xi=1,emingrid%dim(1)
        do yi=1,emingrid%dim(2)
          do zi=1,emingrid%dim(3)
            if (emingrid%dens(xi,yi,zi).gt.emfl) then
              emingrid%dens(xi,yi,zi) = emfl
            end if
          end do
        end do
      end do
    end if
  else if (emflatval.lt.HUGE(emflatval)) then
    write(ilog,*) 'Warning. Requested flattening of input density map is inconsistent with range &
 &of values encountered (',minde+shfv,' to ',maxde+shfv,'). No flattening occurred.'
  end if
!
! now take the shifted values and count the excess mass with the determined peak position
! over the range of values that exceed the threshold (this happens on the shifted, unscaled values
! with the shifted, unscaled threshold)
  totm = 0.0
  tots = 0.0
  totv = 0.0
  voli = emingrid%deltas(1)*emingrid%deltas(2)*emingrid%deltas(3)
  do xi=1,emingrid%dim(1)
    do yi=1,emingrid%dim(2)
      do zi=1,emingrid%dim(3)
        if (emingrid%dens(xi,yi,zi).gt.threshold) then
          totm = totm + (emingrid%dens(xi,yi,zi)-bgd)*voli
          tots = tots + embgdensity*voli
          totv = totv + voli
        end if
      end do
    end do
  end do
  if (totm.eq.0.0) then
    write(ilog,*) 'Fatal. Integrated excess signal is exactly zero for spatial density restraint. Alter threshold level &
 &and/or verify that background density is reasonable.'
    call fexit()
  end if
  if (emtotalmass.lt.1.0) then
    if (empotprop.eq.1) then
      emtotalmass = sum(mass(1:n))
    else if (empotprop.eq.2) then
      emtotalmass = 1.0*sum(lj_atnum(attyp(1:n)))
    end if
  end if
!
! from the specified mass and the threshold, we can already compute a mean physical density
  emavgdensity = convdens*emtotalmass/totv
!  write(*,*) totv,sum(atsavred(1:n)*atvol(1:n)),emtotalmass
! assuming a proper unit conversion factor, we can compare the measured excess mass from the input
! (totm*voli) to the theoretical excess
  emsf2 = (convdens*emtotalmass - tots)/(totm)
!  write(*,*) emsf2
!
! if the excess is negative, throw exception
  if (emsf2.le.0.0) then
    write(ilog,*) 'Fatal. Density conversion relying on provided mass and threshold values &
 &cannot establish a suitable excess over determined background density (',bgd+shfv,' input units). &
 &Increase threshold level and/or assumed mass and/or verify that background density is reasonable.'
    call fexit()
  end if
!
! transform data according to scaling factors (note that emsf2 is irrelevant for all density values
! equivalent to bgd -> solvent peak does not shift on account of emsf2)
  do yi=1,emingrid%dim(2)
    do zi=1,emingrid%dim(3)
      emingrid%dens(:,yi,zi) = embgdensity + emsf2*(emingrid%dens(:,yi,zi)-bgd)
    end do
  end do
!  write(*,*) maxval(emingrid%dens),maxval(emingrid%optdens),minval(emingrid%dens)
!
  thresher = embgdensity + emsf2*(threshold-bgd)
  call strlims(emmapfile,t1,t2)
  if (embgdensity.gt.thresher) then
    write(ilog,*) 'Warning. Assumed threshold signal from input density map (',&
 &emmapfile(t1:t2),') is lower than determined solvent (peak) signal.'
!    write(*,*) bgd+shfv,thresher,threshold,shfv
  end if
!
  truncer = embgdensity + emsf2*((emtruncate-shfv)-bgd)
  if (truncer.ge.thresher) then
    write(ilog,*) 'Warning. Requested truncation level (EMTRUNCATE) that exceeds threshold level (EMTHRESHOLD)&
 & for input density map. Request ignored.'
    emtruncate = -1.0*(HUGE(emtruncate)/1000.0)
  else if (truncer.ge.(embgdensity + emsf2*(minde-bgd))) then
    do xi=1,emingrid%dim(1)
      do yi=1,emingrid%dim(2)
        do zi=1,emingrid%dim(3)
          if (emingrid%dens(xi,yi,zi).lt.truncer) then
            emingrid%dens(xi,yi,zi) = embgdensity
          end if
        end do
      end do
    end do
  end if
!
  threshold = threshold + shfv
!
! transfer to coarse map - note that cells align (otherwise allocate_ems exits fatally)
  do xi=1,emcoarsegrid%dim(1)
    xj = xi + nint((emcoarsegrid%origin(1)-emingrid%origin(1))/emingrid%deltas(1))
    do yi=1,emcoarsegrid%dim(2)
      yj = yi + nint((emcoarsegrid%origin(2)-emingrid%origin(2))/emingrid%deltas(2))
      do zi=1,emcoarsegrid%dim(3)
        zj = zi + nint((emcoarsegrid%origin(3)-emingrid%origin(3))/emingrid%deltas(3))
        if ((xj.ge.1).AND.(yj.ge.1).AND.(zj.ge.1).AND.(xj.le.emingrid%dim(1)).AND.&
 &(yj.le.emingrid%dim(2)).AND.(zj.le.emingrid%dim(3))) then
          emcoarsegrid%dens(xi,yi,zi) = emingrid%dens(xj,yj,zj)
        else
          emcoarsegrid%dens(xi,yi,zi) = embgdensity
        end if
      end do
    end do
  end do
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine choose_simplification()
!
  use ems
  use iounit
  use system
  use atoms
  use sequen
  use polypep
  use params
!
  implicit none
!
  integer i,npts,rfac,ndims,yk,xi,xk,incri
  RTYPE fracvol,fracvol2
!
  fracvol = sysvol(1)/(emgrid%deltas(1)*emgrid%deltas(2)*emgrid%deltas(3)*emgrid%dim(1)*emgrid%dim(2)*emgrid%dim(3))
  fracvol2 = sysvol(1)/(emingrid%deltas(1)*emingrid%deltas(2)*emingrid%deltas(3)*emingrid%dim(1)*emingrid%dim(2)*emingrid%dim(3))
  write(ilog,*) 'With current settings, fractional volume from overlap-corrected atomic volumes is expected &
 &to be around ',fracvol*100,'% in evaluation grid coverage and ',fracvol2*100,'% for input density map coverage&
 & for density restraint potential.'
!
!! WARNING HACK!!!!!!!!!!!!
!  do i=1,n
!    if ((i.eq.at(atmres(i))%bb(5)).AND.(mass(i).gt.15.0).AND.(seqpolty(atmres(i)).eq.'P')) then
!      mass(i) = 1.0
!!      lj_atnum(attyp(i)) = 1
!      write(*,*) i,' is 2OXT'
!    end if
!  end do
!  write(*,*) sum(mass(1:n)),sum(lj_atnum(attyp(1:n)))
!
  if (emheuristic.eq.1) then
    if (emcoarsegrid%dim(1).le.5) then
      write(ilog,*) 'Warning. First dimension should have a significant number of bins when using &
 &slice heuristic in density restraint potential.'
    end if
    allocate(emcoarsegrid%lslc(emcoarsegrid%dim(2),2))
    allocate(emcoarsegrid%rslc(emcoarsegrid%dim(2)))
  else if (emheuristic.eq.2) then
    allocate(emcoarsegrid%llin(emcoarsegrid%dim(2),emcoarsegrid%dim(3),2))
    allocate(emcoarsegrid%rlin(emcoarsegrid%dim(2),emcoarsegrid%dim(3)))
  else if (emheuristic.eq.3) then
    npts = emcoarsegrid%dim(1)*emcoarsegrid%dim(2)*emcoarsegrid%dim(3)
    emcoarsegrid%dimc(:) = 0
    rfac = 1
    ndims = 3
    do i=1,3
      if (emcoarsegrid%dim(i).le.3) then
        emcoarsegrid%dimc(i) = emcoarsegrid%dim(i)
        ndims = ndims - 1
      else if (emcoarsegrid%dim(i).le.10) then
        if (mod(emcoarsegrid%dim(i),2).eq.0) then! 4,6,8,10
          emcoarsegrid%dimc(i) = emcoarsegrid%dim(i)/2
          rfac = rfac*2
          ndims = ndims - 1
        else
          emcoarsegrid%dimc(i) = ceiling(emcoarsegrid%dim(i)/2.0)
          rfac = rfac*2
          ndims = ndims - 1
        end if
      end if
    end do
    if (ndims.le.0) then
!     do nothing
    else
      if (npts.ge.8000) then ! 20x20x20
!       aim for about 1000 blocks
        fracvol = (dble(npts)/(rfac*1000.0))**(1.0/dble(ndims))
        do i=1,3
          if (emcoarsegrid%dimc(i).le.0) then
            emcoarsegrid%dimc(i) = ceiling(emcoarsegrid%dim(i)/max(fracvol,2.0))
          end if
        end do
      else if (npts.gt.100) then
        fracvol = (dble(npts)/(rfac*100.0))**(1.0/dble(ndims))
        if (fracvol.gt.1.5) then
          do i=1,3
            if (emcoarsegrid%dimc(i).le.0) then
              emcoarsegrid%dimc(i) = ceiling(emcoarsegrid%dim(i)/max(fracvol,2.0))
            end if
          end do
        else
          do i=1,3
            if (emcoarsegrid%dimc(i).le.0) then
              emcoarsegrid%dimc(i) = emcoarsegrid%dim(i)
            end if
          end do
        end if
      end if
    end if
    allocate(emcoarsegrid%blklms(maxval(emcoarsegrid%dimc(:)),2,3))
    do yk=1,3
      xk = 0
      if (mod(emcoarsegrid%dim(yk),emcoarsegrid%dimc(yk)).eq.0) then
        incri = emcoarsegrid%dim(yk)/emcoarsegrid%dimc(yk)
      else
        incri = ceiling(dble(emcoarsegrid%dim(yk))/dble(emcoarsegrid%dimc(yk)))
      end if
      do xi=1,emcoarsegrid%dim(yk),incri
        xk = xk + 1
        emcoarsegrid%blklms(xk,1,yk) = xi
        emcoarsegrid%blklms(xk,2,yk) = min(xi+incri-1,emcoarsegrid%dim(yk))
      end do
      emcoarsegrid%dimc(yk) = xk
    end do
    allocate(emcoarsegrid%lblk(emcoarsegrid%dimc(1),emcoarsegrid%dimc(2),emcoarsegrid%dimc(3),2))
    allocate(emcoarsegrid%rblk(emcoarsegrid%dimc(1),emcoarsegrid%dimc(2),emcoarsegrid%dimc(3)))
  else
!   do nothing
  end if
!
end
!
!------------------------------------------------------------------------------------------
!
! this routine precomputes the difference contributions for empty yz-slices and empty z-lines
! these values are used in energy computation to simplify cases for large grids
!
subroutine precompute_diff_emmap()
!
  use ems
!
  implicit none
!
  RTYPE cubesum2,xsc,ysc,zsc
  integer xj,yj,zj,xk,yk,zk,ints(3)
  logical doslice,doline,doblocks
!
  doslice = .false.
  doline = .false.
  doblocks = .false.
!
  ints(:) = nint(emcoarsegrid%deltas(:)/(emgrid%deltas(:)))
! 
  if (allocated(emcoarsegrid%rslc).EQV..true.) then
    emcoarsegrid%rslc(:) = 0.0
    doslice = .true.
  end if
  if (allocated(emcoarsegrid%rlin).EQV..true.) then
    emcoarsegrid%rlin(:,:) = 0.0
    doline = .true.
  end if
  if (allocated(emcoarsegrid%rblk).EQV..true.) then
    emcoarsegrid%rblk(:,:,:) = 0.0
    xsc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(1))/dble(emcoarsegrid%dimc(1))))
    ysc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(2))/dble(emcoarsegrid%dimc(2))))
    zsc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(3))/dble(emcoarsegrid%dimc(3))))
    doblocks = .true.
  end if
  if ((doslice.EQV..false.).AND.(doline.EQV..false.).AND.(doblocks.EQV..false.)) return
!
! this loop triggers (most likely) a vectorization bug in ORACLE Studio 12.3 with -xO3 or higher
  if (doblocks.EQV..true.) then
    do zj=1,emcoarsegrid%dim(3)
      zk = ceiling(zsc*(zj-0.5))
      do yj=1,emcoarsegrid%dim(2)
        yk = ceiling(ysc*(yj-0.5))
        do xj=1,emcoarsegrid%dim(1)
          xk = ceiling(xsc*(xj-0.5))
          cubesum2 = (embgdensity - emcoarsegrid%dens(xj,yj,zj))*(embgdensity - emcoarsegrid%dens(xj,yj,zj))
          emcoarsegrid%rblk(xk,yk,zk) = emcoarsegrid%rblk(xk,yk,zk) + cubesum2
        end do
      end do
    end do
  else if (doline.EQV..true.) then
    do zj=1,emcoarsegrid%dim(3)
      do yj=1,emcoarsegrid%dim(2)
        do xj=1,emcoarsegrid%dim(1)
          cubesum2 = (embgdensity - emcoarsegrid%dens(xj,yj,zj))*(embgdensity - emcoarsegrid%dens(xj,yj,zj))
          emcoarsegrid%rlin(yj,zj) = emcoarsegrid%rlin(yj,zj) + cubesum2
        end do
      end do
    end do
  else if (doslice.EQV..true.) then
    do zj=1,emcoarsegrid%dim(3)
      do yj=1,emcoarsegrid%dim(2)
        do xj=1,emcoarsegrid%dim(1)
          cubesum2 = (embgdensity - emcoarsegrid%dens(xj,yj,zj))*(embgdensity - emcoarsegrid%dens(xj,yj,zj))
          emcoarsegrid%rslc(yj) = emcoarsegrid%rslc(yj) + cubesum2
        end do
      end do
    end do
  end if
!
end
!
!-----------------------------------------------------------------------------------------------
!
! a routine to prepopulate some maps to avoid expensive mapping in inner loops
! should be called after all other setup of EM is complete
!
subroutine em_getintmaps()
!
  use ems
  use energies
!
  implicit none
!
  integer i,xi,yi,zi,sh,ints(3)
!
  if (mod(emsplor,2).eq.0) then
    sh = emsplor/2 
!    oriv2(:) = it%origin(:) - it%deltas(:)
  else
    sh = floor(emsplor/2.0)
!    oriv2(:) = it%origin(:) - 0.5*it%deltas(:)
  end if
!
  allocate(emgrid%xmap(3*emgrid%dim(1)))
  allocate(emgrid%ymap(3*emgrid%dim(2)))
  allocate(emgrid%zmap(3*emgrid%dim(3)))
!
  do i=1,3*emgrid%dim(1)
    xi = i - sh
    if (xi.gt.emgrid%dim(1)) xi = xi - emgrid%dim(1)
    if (xi.gt.emgrid%dim(1)) xi = xi - emgrid%dim(1)
    if (xi.gt.emgrid%dim(1)) xi = xi - emgrid%dim(1)
    if (xi.lt.1) xi = xi + emgrid%dim(1)
    if (xi.lt.1) xi = xi + emgrid%dim(1)
    if (xi.lt.1) xi = xi + emgrid%dim(1)
    emgrid%xmap(i) = xi
  end do
  do i=1,3*emgrid%dim(2)
    yi = i - sh
    if (yi.gt.emgrid%dim(2)) yi = yi - emgrid%dim(2)
    if (yi.gt.emgrid%dim(2)) yi = yi - emgrid%dim(2)
    if (yi.gt.emgrid%dim(2)) yi = yi - emgrid%dim(2)
    if (yi.lt.1) yi = yi + emgrid%dim(2)
    if (yi.lt.1) yi = yi + emgrid%dim(2)
    if (yi.lt.1) yi = yi + emgrid%dim(2)
    emgrid%ymap(i) = yi
  end do
  do i=1,3*emgrid%dim(3)
    zi = i - sh
    if (zi.gt.emgrid%dim(3)) zi = zi - emgrid%dim(3)
    if (zi.gt.emgrid%dim(3)) zi = zi - emgrid%dim(3)
    if (zi.gt.emgrid%dim(3)) zi = zi - emgrid%dim(3)
    if (zi.lt.1) zi = zi + emgrid%dim(3)
    if (zi.lt.1) zi = zi + emgrid%dim(3)
    if (zi.lt.1) zi = zi + emgrid%dim(3)
    emgrid%zmap(i) = zi
  end do
!
  if (use_EMICRO.EQV..true.) then
!
    ints(:) = nint(emcoarsegrid%deltas(:)/(emgrid%deltas(:)))
    allocate(emgrid%xmap2(emgrid%dim(1)))
    allocate(emgrid%ymap2(emgrid%dim(2)))
    allocate(emgrid%zmap2(emgrid%dim(3)))
!
    do i=1,emgrid%dim(1)
      emgrid%xmap2(i) = ceiling((i*1.0-0.5)/(1.0*ints(1)))
    end do
    do i=1,emgrid%dim(2)
      emgrid%ymap2(i) = ceiling((i*1.0-0.5)/(1.0*ints(2)))
    end do
    do i=1,emgrid%dim(3)
      emgrid%zmap2(i) = ceiling((i*1.0-0.5)/(1.0*ints(3)))
    end do
!
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
subroutine do_em_density_map(tpi)
!
  use ems
  use energies
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
#ifdef ENABLE_THREADS
  integer(KIND=8) ttt(2)
#endif
!
  logical afalse
!
  afalse = .false.
  if (tpi.le.1) emgrid%cnt = emgrid%cnt + 1
! we assume that we do have a current version available in %mass
  if (use_EMICRO.EQV..false.) then
#ifdef ENABLE_THREADS
    if (tpi.eq.1) call System_Clock(count=ttt(1))
    call em_pop_mv_threads(afalse,tpi)
    if (tpi.eq.1) call System_Clock(count=ttt(2))
!    if (tpi.eq.1) write(*,*) 1.0e-3*(ttt(2)-ttt(1))
#else
    call em_pop_mv(emgrid,afalse)
#endif
  end if
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    emgrid%avgmass(:,:,thr_limits(49,tpi):thr_limits(50,tpi)) = emgrid%avgmass(:,:,thr_limits(49,tpi):thr_limits(50,tpi)) + &
 &                                                            emgrid%mass(:,:,thr_limits(49,tpi):thr_limits(50,tpi))
  else
#endif
  emgrid%avgmass(:,:,:) = emgrid%avgmass(:,:,:) + emgrid%mass(:,:,:)
#ifdef ENABLE_THREADS
  end if
#endif
!
end
!
!------------------------------------------------------------------------------------------
!
subroutine em_transfer_m(tpi)
!
  use ems
  use energies
#ifdef ENABLE_THREADS
  use threads, ONLY: thr_limits
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
#ifdef ENABLE_THREADS
  if (tpi.gt.0) then
    emgrid%mass(:,:,thr_limits(49,tpi):thr_limits(50,tpi)) = emgrid%mass2(:,:,thr_limits(49,tpi):thr_limits(50,tpi))
    if (allocated(emcoarsegrid%lslc).EQV..true.) &
 &emcoarsegrid%lslc(thr_limits(51,tpi):thr_limits(52,tpi),1) = emcoarsegrid%lslc(thr_limits(51,tpi):thr_limits(52,tpi),2)
    if (allocated(emcoarsegrid%llin).EQV..true.) &
 &emcoarsegrid%llin(:,thr_limits(57,tpi):thr_limits(58,tpi),1) = emcoarsegrid%llin(:,thr_limits(57,tpi):thr_limits(58,tpi),2)
    if (allocated(emcoarsegrid%lblk).EQV..true.) &
 &emcoarsegrid%lblk(:,:,thr_limits(51,tpi):thr_limits(52,tpi),1) = emcoarsegrid%lblk(:,:,thr_limits(51,tpi):thr_limits(52,tpi),2)
  else
#endif
  emgrid%mass(:,:,:) = emgrid%mass2(:,:,:)
  if (allocated(emcoarsegrid%lslc).EQV..true.) emcoarsegrid%lslc(:,1) = emcoarsegrid%lslc(:,2)
  if (allocated(emcoarsegrid%llin).EQV..true.) emcoarsegrid%llin(:,:,1) = emcoarsegrid%llin(:,:,2)
  if (allocated(emcoarsegrid%lblk).EQV..true.) emcoarsegrid%lblk(:,:,:,1) = emcoarsegrid%lblk(:,:,:,2)
#ifdef ENABLE_THREADS
  end if
#endif
!
end

!
!----------------------------------------------------------------------------------
!
! note that this routine silently sets the global, atom-based B-spline parameters and
! the global heuristic arrays in emcoarsegrid despite taking an argument for the emgrid structure 
! 
subroutine em_pop_mv(it,usemass2)
!
  use atoms
  use ems
  use units
  use iounit
  use params
  use system
  use grandensembles
  use molecule
  use energies
!
  implicit none
!
  type(t_emgrid) it
  integer xi,yi,zi,imol,ii,pv(3),j,k,l,zk,yk,idx,ints(3)
  logical doblocks,doslice,doline,usemass2,ismember
  RTYPE unitvol,soldens,gdzv(3),mvec(3),temj,temk,teml,tem0,term0,oriv2(3),shf,xsc,ysc,zsc
!
  unitvol = it%deltas(1)*it%deltas(2)*it%deltas(3)
  if (use_EMICRO.EQV..true.) then
    ints(:) = nint(emcoarsegrid%deltas(:)/(emgrid%deltas(:)))
  else
    ints(:) = 1
  end if
  if (empotprop.le.2) then
    soldens = embgdensity/convdens ! in a.m.u per A^3
  else if (empotprop.ge.3) then
    soldens = embgdensity ! assumed in whatever/A^3
  else
    write(ilog,*) 'Fatal. Encountered unsupported property for density analysis in em_pop_mv(...). This is a bug.'
    call fexit()
  end if
  term0 = soldens*unitvol
  if (usemass2.EQV..true.) then
    it%mass2(:,:,:) = term0
    idx = 2
  else
    it%mass(:,:,:) = term0
    idx = 1
  end if
  doslice = .false.
  if (allocated(emcoarsegrid%lslc).EQV..true.) then
    emcoarsegrid%lslc(:,idx) = .false.
    doslice = .true.
  end if
  doline = .false.
  if (allocated(emcoarsegrid%llin).EQV..true.) then
    emcoarsegrid%llin(:,:,idx) = .false.
    doline = .true.
  end if
  doblocks = .false.
  if (allocated(emcoarsegrid%lblk).EQV..true.) then
    doblocks = .true.
    emcoarsegrid%lblk(:,:,:,idx) = .false.
    xsc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(1))/dble(emcoarsegrid%dimc(1))))
    ysc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(2))/dble(emcoarsegrid%dimc(2))))
    zsc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(3))/dble(emcoarsegrid%dimc(3))))
    shf = -0.5
  end if
!
! spread scaled coordinates on grid by incrementing the spline array
! note that we use PBC(!)
  gdzv(:) = 1.0/it%deltas(:)
  if (mod(emsplor,2).eq.0) then
    oriv2(:) = it%origin(:) - it%deltas(:)
  else
    oriv2(:) = it%origin(:) - 0.5*it%deltas(:)
  end if
  if (soldens.eq.0.0) then
    do imol=1,nmol
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
          if (ismember(ispresent,imol).EQV..false.) cycle
        end if
      end if
      do ii=atmol(imol,1),atmol(imol,2)
        if (empotprop.eq.1) then
          tem0 = mass(ii)
        else if (empotprop.eq.2) then
          tem0 = 1.0*lj_atnum(attyp(ii))
        else if (empotprop.eq.3) then
          tem0 = atq(ii)
        else if (empotprop.eq.4) then
          tem0 = emcustom(ii,1)
        end if
        if (tem0.eq.0.0) cycle
        mvec(1) = gdzv(1)*(x(ii)-oriv2(1))
        mvec(2) = gdzv(2)*(y(ii)-oriv2(2))
        mvec(3) = gdzv(3)*(z(ii)-oriv2(3))
        pv(1:3) = floor(mvec(1:3))
!       first get the Mn(xi), i.e., the relevant B-spline values for the grid-points 
!       in the vicinity of atom ii: this can be done separately for all axes
!       get the d/dxi (Mn(xi)): note that the grid axes are aligned with the coordinate axes,
!       so only diagonal terms contribute, i.e., d/dyi (Mn(xi)) is obviously zero
        call cardBspline(mvec(1)-1.0*pv(1),emsplor,embspl(:,1,ii),embspld(:,1,ii))
        call cardBspline(mvec(2)-1.0*pv(2),emsplor,embspl(:,2,ii),embspld(:,2,ii))
        call cardBspline(mvec(3)-1.0*pv(3),emsplor,embspl(:,3,ii),embspld(:,3,ii))
        embspld(1:emsplor,1,ii) = embspld(1:emsplor,1,ii)*gdzv(1)
        embspld(1:emsplor,2,ii) = embspld(1:emsplor,2,ii)*gdzv(2)
        embspld(1:emsplor,3,ii) = embspld(1:emsplor,3,ii)*gdzv(3)
        emgfls(ii,:) = pv(:) + emgrid%dim(1:3)
        do k=1,emsplor
          yi = emgrid%ymap(emgfls(ii,2) + k)
          temk = embspl(k,2,ii)
          if (doslice.EQV..true.) emcoarsegrid%lslc(emgrid%ymap2(yi),idx) = .true.
          if (doblocks.EQV..true.) yk = ceiling(ysc*(emgrid%ymap2(yi)+shf))
          do l=1,emsplor
            zi = emgrid%zmap(emgfls(ii,3) + l) 
            teml = temk*embspl(l,3,ii)
            if (doline.EQV..true.) emcoarsegrid%llin(emgrid%ymap2(yi),emgrid%zmap2(zi),idx) = .true.
            if (doblocks.EQV..true.) zk = ceiling(zsc*(emgrid%zmap2(zi)+shf))
            do j=1,emsplor
              temj = teml*embspl(j,1,ii)
              xi = emgrid%xmap(emgfls(ii,1) + j)
              emgrid%xmap2(xi) = emgrid%xmap2(xi) 
              if (usemass2.EQV..true.) then
                it%mass2(xi,yi,zi) = it%mass2(xi,yi,zi) + temj*tem0
              else
                it%mass(xi,yi,zi) = it%mass(xi,yi,zi) + temj*tem0
              end if
              if (doblocks.EQV..true.) emcoarsegrid%lblk(ceiling(xsc*(emgrid%xmap2(xi)+shf)),yk,zk,idx) = .true.
            end do
          end do
        end do
      end do
    end do
  else
    do imol=1,nmol
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
          if (ismember(ispresent,imol).EQV..false.) cycle
        end if
      end if
      do ii=atmol(imol,1),atmol(imol,2)
        term0 = soldens*atsavred(ii)*atvol(ii)
        if (empotprop.eq.1) then
          tem0 = mass(ii) - term0
        else if (empotprop.eq.2) then
          tem0 = 1.0*lj_atnum(attyp(ii)) - term0
        else if (empotprop.eq.3) then
          tem0 = atq(ii) - term0
        else if (empotprop.eq.4) then
          tem0 = emcustom(ii,1) - term0
        end if
        mvec(1) = gdzv(1)*(x(ii)-oriv2(1))
        mvec(2) = gdzv(2)*(y(ii)-oriv2(2))
        mvec(3) = gdzv(3)*(z(ii)-oriv2(3))
        pv(1:3) = floor(mvec(1:3))
!       first get the Mn(xi), i.e., the relevant B-spline values for the grid-points 
!       in the vicinity of atom ii: this can be done separately for all axes
!       get the d/dxi (Mn(xi)): note that the grid axes are aligned with the coordinate axes,
!       so only diagonal terms contribute, i.e., d/dyi (Mn(xi)) is obviously zero
        call cardBspline(mvec(1)-1.0*pv(1),emsplor,embspl(:,1,ii),embspld(:,1,ii))
        call cardBspline(mvec(2)-1.0*pv(2),emsplor,embspl(:,2,ii),embspld(:,2,ii))
        call cardBspline(mvec(3)-1.0*pv(3),emsplor,embspl(:,3,ii),embspld(:,3,ii))
        embspld(1:emsplor,1,ii) = embspld(1:emsplor,1,ii)*gdzv(1)
        embspld(1:emsplor,2,ii) = embspld(1:emsplor,2,ii)*gdzv(2)
        embspld(1:emsplor,3,ii) = embspld(1:emsplor,3,ii)*gdzv(3)
        emgfls(ii,:) = pv(:) + emgrid%dim(1:3)
        do k=1,emsplor
          yi = emgrid%ymap(emgfls(ii,2) + k)
          temk = embspl(k,2,ii)
          if (doslice.EQV..true.) emcoarsegrid%lslc(emgrid%ymap2(yi),idx) = .true.
          if (doblocks.EQV..true.) yk = ceiling(ysc*(emgrid%ymap2(yi)+shf))
          do l=1,emsplor
            zi = emgrid%zmap(emgfls(ii,3) + l)
            teml = temk*embspl(l,3,ii)
            if (doline.EQV..true.) emcoarsegrid%llin(emgrid%ymap2(yi),emgrid%zmap2(zi),idx) = .true.
            if (doblocks.EQV..true.) zk = ceiling(zsc*(emgrid%zmap2(zi)+shf))
            do j=1,emsplor
              temj = teml*embspl(j,1,ii)
              xi = emgrid%xmap(emgfls(ii,1) + j)
              emgrid%xmap2(xi) = emgrid%xmap2(xi)
              if (usemass2.EQV..true.) then
                it%mass2(xi,yi,zi) = it%mass2(xi,yi,zi) + temj*tem0
              else
                it%mass(xi,yi,zi) = it%mass(xi,yi,zi) + temj*tem0
              end if
              if (doblocks.EQV..true.) emcoarsegrid%lblk(ceiling(xsc*(emgrid%xmap2(xi)+shf)),yk,zk,idx) = .true.
            end do
          end do
        end do
      end do
    end do
  end if
!
end
!
!----------------------------------------------------------------------------------
!
#ifdef ENABLE_THREADS

! a threaded alternative for which the grid argument is not passed 
! 
subroutine em_pop_mv_threads(usemass2,tpi)
!
  use atoms
  use ems
  use units
  use iounit
  use params
  use system
  use grandensembles
  use molecule
  use energies
  use threads
  use sequen
!
  implicit none
!
  integer, INTENT(IN):: tpi
  logical, INTENT(IN):: usemass2
!
  integer xi,yi,zi,imol,ii,pv(3),i,j,k,l,idx,ints(3),mlo,mhi
  logical doblocks,doslice,doline,ismember,slcl(emcoarsegrid%dim(2)),lncl(emcoarsegrid%dim(2),emcoarsegrid%dim(3))
  logical blkl(emcoarsegrid%dimc(1),emcoarsegrid%dimc(2),emcoarsegrid%dimc(3))
  RTYPE unitvol,soldens,gdzv(3),mvec(3),ts(3),temj,temk,teml,tem0,term0,oriv2(3),shf,xsc,ysc,zsc
  integer xmap2(emgrid%dim(1)),ymap2(emgrid%dim(2)),zmap2(emgrid%dim(3))
  integer(KIND=8) ttimer
!
  if (tpi.le.0) return
!
  unitvol = emgrid%deltas(1)*emgrid%deltas(2)*emgrid%deltas(3)
  if (use_EMICRO.EQV..true.) then
    ints(:) = nint(emcoarsegrid%deltas(:)/(emgrid%deltas(:)))
  else
    ints(:) = 1
  end if
  if (empotprop.le.2) then
    soldens = embgdensity/convdens ! in a.m.u per A^3
  else if (empotprop.ge.3) then
    soldens = embgdensity ! assumed in whatever/A^3
  else
    write(ilog,*) 'Fatal. Encountered unsupported property for density analysis in em_pop_mv_threads(...). This is a bug.'
    call fexit()
  end if
  term0 = soldens*unitvol
!
  if (usemass2.EQV..true.) then
    emgrid%mass2(:,:,thr_limits(49,tpi):thr_limits(50,tpi)) = term0
    idx = 2
  else
    emgrid%mass(:,:,thr_limits(49,tpi):thr_limits(50,tpi)) = term0
    idx = 1
  end if
  doslice = .false.
  if (allocated(emcoarsegrid%lslc).EQV..true.) then
    emcoarsegrid%lslc(thr_limits(51,tpi):thr_limits(52,tpi),idx) = .false.
    doslice = .true.
    slcl(:) = .false.
  end if
  doline = .false.
  if (allocated(emcoarsegrid%llin).EQV..true.) then
    emcoarsegrid%llin(:,thr_limits(57,tpi):thr_limits(58,tpi),idx) = .false.
    doline = .true.
    lncl(:,:) = .false.
  end if
  doblocks = .false.
  if (allocated(emcoarsegrid%lblk).EQV..true.) then
    doblocks = .true.
    emcoarsegrid%lblk(:,:,thr_limits(51,tpi):thr_limits(52,tpi),idx) = .false.
    xsc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(1))/dble(emcoarsegrid%dimc(1))))
    ysc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(2))/dble(emcoarsegrid%dimc(2))))
    zsc = 1.0/dble(ceiling(dble(emcoarsegrid%dim(3))/dble(emcoarsegrid%dimc(3))))
    shf = -0.5
    blkl(:,:,:) = .false.
  end if
!
! spread scaled coordinates on grid by incrementing the spline array
! note that we use PBC(!)
  gdzv(:) = 1.0/emgrid%deltas(:)
  if (mod(emsplor,2).eq.0) then
    oriv2(:) = emgrid%origin(:) - emgrid%deltas(:)
  else
    oriv2(:) = emgrid%origin(:) - 0.5*emgrid%deltas(:)
  end if
!
  if (doblocks.EQV..true.) then
    do i=1,emgrid%dim(1)
      xmap2(i) = ceiling(xsc*(ceiling((1.0*i-0.5)/(1.0*ints(1)))+shf))
    end do
    do i=1,emgrid%dim(2)
      ymap2(i) = ceiling(ysc*(ceiling((1.0*i-0.5)/(1.0*ints(2)))+shf))
    end do
    do i=1,emgrid%dim(3)
      zmap2(i) = ceiling(zsc*(ceiling((1.0*i-0.5)/(1.0*ints(3)))+shf))
    end do
  end if
!
!$OMP BARRIER
!
  if (thr_limits(54,tpi).ge.thr_limits(53,tpi)) then
    mlo = molofrs(atmres(thr_limits(53,tpi)))
    mhi = molofrs(atmres(thr_limits(54,tpi)))
  else
    mlo = 1
    mhi = 0
  end if
  do imol=mlo,mhi
    if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
      if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
        if (ismember(ispresent,imol).EQV..false.) cycle
      end if
    end if
    do ii=max(thr_limits(53,tpi),atmol(imol,1)),min(thr_limits(54,tpi),atmol(imol,2))
      if (empotprop.eq.1) then
        tem0 = mass(ii)
      else if (empotprop.eq.2) then
        tem0 = 1.0*lj_atnum(attyp(ii))
      else if (empotprop.eq.3) then
        tem0 = atq(ii)
      else if (empotprop.eq.4) then
        tem0 = emcustom(ii,1)
      end if
      if ((tem0.eq.0.0).AND.(soldens.eq.0.0)) cycle
      mvec(1) = gdzv(1)*(x(ii)-oriv2(1))
      mvec(2) = gdzv(2)*(y(ii)-oriv2(2))
      mvec(3) = gdzv(3)*(z(ii)-oriv2(3))
      pv(1:3) = floor(mvec(1:3))
      ts(1:3) = mvec(1:3) - pv(1:3) 
      call cardBspline(ts(1),emsplor,embspl(:,1,ii),embspld(:,1,ii))
      call cardBspline(ts(2),emsplor,embspl(:,2,ii),embspld(:,2,ii))
      call cardBspline(ts(3),emsplor,embspl(:,3,ii),embspld(:,3,ii))
      embspld(1:emsplor,1,ii) = embspld(1:emsplor,1,ii)*gdzv(1)
      embspld(1:emsplor,2,ii) = embspld(1:emsplor,2,ii)*gdzv(2)
      embspld(1:emsplor,3,ii) = embspld(1:emsplor,3,ii)*gdzv(3)
      emgfls(ii,:) = pv(:) + emgrid%dim(1:3)
    end do
  end do
!$OMP BARRIER
  if (thr_dlb(11,1).gt.0) then
    if (tpi.eq.1) thr_dlb(11,2) = thr_dlb(11,2) + 1
    call System_Clock(count=ttimer)
    thr_timings(25,tpi) = thr_timings(25,tpi) + ttimer
  end if
  if (soldens.eq.0.0) then
    do imol=1,nmol
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
          if (ismember(ispresent,imol).EQV..false.) cycle
        end if
      end if
      do ii=atmol(imol,1),atmol(imol,2)
        if (empotprop.eq.1) then
          tem0 = mass(ii)
        else if (empotprop.eq.2) then
          tem0 = 1.0*lj_atnum(attyp(ii))
        else if (empotprop.eq.3) then
          tem0 = atq(ii)
        else if (empotprop.eq.4) then
          tem0 = emcustom(ii,1)
        end if
        if (tem0.eq.0.0) cycle
        do l=1,emsplor
          zi = emgrid%zmap(emgfls(ii,3) + l)
          if (zi.lt.thr_limits(55,tpi)) cycle
          if (zi.gt.thr_limits(56,tpi)) cycle
          teml = embspl(l,3,ii)
          do k=1,emsplor
            yi = emgrid%ymap(emgfls(ii,2) + k)
            temk = teml*embspl(k,2,ii)
            if (doslice.EQV..true.) emcoarsegrid%lslc(emgrid%ymap2(yi),idx) = .true.
            if (doline.EQV..true.) emcoarsegrid%llin(emgrid%ymap2(yi),emgrid%zmap2(zi),idx) = .true.
            do j=1,emsplor
              temj = temk*embspl(j,1,ii)
              xi = emgrid%xmap(emgfls(ii,1) + j)
              if (usemass2.EQV..true.) then
                emgrid%mass2(xi,yi,zi) = emgrid%mass2(xi,yi,zi) + temj*tem0
              else
                emgrid%mass(xi,yi,zi) = emgrid%mass(xi,yi,zi) + temj*tem0
              end if
              if (doblocks.EQV..true.) emcoarsegrid%lblk(xmap2(xi),ymap2(yi),zmap2(zi),idx) = .true.
            end do
          end do
        end do
      end do
    end do
  else
    do imol=1,nmol
      if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
        if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
          if (ismember(ispresent,imol).EQV..false.) cycle
        end if
      end if
      do ii=atmol(imol,1),atmol(imol,2)
        term0 = soldens*atsavred(ii)*atvol(ii)
        if (empotprop.eq.1) then
          tem0 = mass(ii) - term0
        else if (empotprop.eq.2) then
          tem0 = 1.0*lj_atnum(attyp(ii)) - term0
        else if (empotprop.eq.3) then
          tem0 = atq(ii) - term0
        else if (empotprop.eq.4) then
          tem0 = emcustom(ii,1) - term0
        end if
        do l=1,emsplor
          zi = emgrid%zmap(emgfls(ii,3) + l)
          if (zi.lt.thr_limits(55,tpi)) cycle
          if (zi.gt.thr_limits(56,tpi)) cycle
          teml = embspl(l,3,ii)
          do k=1,emsplor
            yi = emgrid%ymap(emgfls(ii,2) + k)
            temk = teml*embspl(k,2,ii)
            if (doslice.EQV..true.) slcl(emgrid%ymap2(yi)) = .true.
            if (doline.EQV..true.) lncl(emgrid%ymap2(yi),emgrid%zmap2(zi)) = .true.
            do j=1,emsplor
              temj = temk*embspl(j,1,ii)
              xi = emgrid%xmap(emgfls(ii,1) + j)
              if (usemass2.EQV..true.) then
                emgrid%mass2(xi,yi,zi) = emgrid%mass2(xi,yi,zi) + temj*tem0
              else
                emgrid%mass(xi,yi,zi) = emgrid%mass(xi,yi,zi) + temj*tem0
              end if
              if (doblocks.EQV..true.) blkl(xmap2(xi),ymap2(yi),zmap2(zi)) = .true.
            end do
          end do
        end do
      end do
    end do
  end if

!$OMP CRITICAL(EMHEURISTIC_UP)
  if (doslice.EQV..true.) then
    do i=1,emcoarsegrid%dim(2)
      if (slcl(i).EQV..true.) emcoarsegrid%lslc(i,idx) = .true.
    end do
  else if (doline.EQV..true.) then
    do j=1,emcoarsegrid%dim(3)
      do i=1,emcoarsegrid%dim(2)
        if (lncl(i,j).EQV..true.) emcoarsegrid%llin(i,j,idx) = .true.
      end do
    end do
  else if (doblocks.EQV..true.) then
    do j=1,emcoarsegrid%dimc(3)
      do i=1,emcoarsegrid%dimc(2)
        do k=1,emcoarsegrid%dimc(1)
          if (blkl(k,i,j).EQV..true.) emcoarsegrid%lblk(k,i,j,idx) = .true.
        end do
      end do
    end do
  end if
!$OMP END CRITICAL(EMHEURISTIC_UP)

  if (thr_dlb(11,1).gt.0) then
    call System_Clock(count=ttimer)
    thr_timings(26,tpi) = thr_timings(26,tpi) + ttimer
  end if
!
!$OMP BARRIER
!
end
!
#endif
!
!------------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
subroutine prt_emmap(thatgrid,mode)
!
  use iounit
  use ems
  use netcdf
  use system
  use mpistuff
  use units
!
  implicit none
!
  type(t_emgrid) thatgrid
  integer intin,dims(3),ii,jj,ncid,freeunit,mode
  character(MAXSTRLEN) attstring,dumpfile
  RTYPE normer
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  logical dohead,exists
!
  if (mode.eq.1) then ! analysis print out
!
    if (thatgrid%cnt.le.0) return
    normer = convdens/(thatgrid%cnt*(thatgrid%deltas(1)*thatgrid%deltas(2)*thatgrid%deltas(3)))
    thatgrid%dens(:,:,:) = normer*thatgrid%avgmass(:,:,:)
!
#ifdef ENABLE_MPI
    call int2str(myrank,nod,re_aux(10))
    if (use_MPIAVG.EQV..true.) then
      dumpfile = 'DENSITY.nc'
    else if (use_REMC.EQV..true.) then
      dumpfile = 'N_'//nod(1:re_aux(10))//'_DENSITY.nc'
    end if
#else
    dumpfile = 'DENSITY.nc'
#endif
!
  else if (mode.eq.2) then ! read-in print out
!
#ifdef ENABLE_MPI
    if (myrank.ne.0) return
    dumpfile = 'DENSITY_INPUT_PHYS.nc'
#else
    dumpfile = 'DENSITY_INPUT_PHYS.nc'
#endif
!
  else if (mode.eq.3) then ! read-in print out
!
#ifdef ENABLE_MPI
    if (myrank.ne.0) return
    dumpfile = 'DENSITY_INPUT_PHYS2.nc'
#else
    dumpfile = 'DENSITY_INPUT_PHYS2.nc'
#endif
!
  else 
!
    write(ilog,*) 'Fatal. Called prt_emmap(...) with unsupported mode. This is a bug.'
    call fexit()
!
  end if
!
  call strlims(dumpfile,ii,jj)
  inquire(file=dumpfile(ii:jj),exist=exists)
  ncid = freeunit()
  call check_fileionetcdf( nf90_create(dumpfile(ii:jj), NF90_CLOBBER, ncid) )
!
  dohead = .TRUE.
!
  if (dohead.EQV..true.) then
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "title", basename(1:bleng)) )
    attstring(1:3) = "   "
    attstring(1:7) = "CAMPARI"
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "program", attstring(1:7)) )
    attstring(1:7) = "       "
    attstring(1:3) = "1.1"
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "programVersion", attstring(1:3)) )
    attstring(1:3) = "   "
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "xyz_origin", thatgrid%origin(1:3)) )
    call check_fileionetcdf( nf90_put_att(ncid, NF90_GLOBAL, "xyz_step", thatgrid%deltas(1:3)) )
!   grid extent dimensions
    attstring(1:1) = "X"
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:1), thatgrid%dim(1), emnetcdf_ids(1)) )
    attstring(1:1) = "Y"
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:1), thatgrid%dim(2), emnetcdf_ids(2)) )
    attstring(1:1) = "Z"
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:1), thatgrid%dim(3), emnetcdf_ids(3)) )
!   3D
    attstring(1:7) = "spatial"
    intin = 3
    call check_fileionetcdf( nf90_def_dim(ncid, attstring(1:7), intin, emnetcdf_ids(4)) )
    attstring(1:7) = "       "
!   grid spacings and units
    attstring(1:6) = "deltas"
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:6), NF90_FLOAT, emnetcdf_ids(4), emnetcdf_ids(5)) )
    attstring(1:6) = "      "
    attstring(1:5) = "units"
    attstring(6:13) = "angstrom"
    call check_fileionetcdf( nf90_put_att(ncid, emnetcdf_ids(5), attstring(1:5), attstring(6:13)) )
!   densities
    dims(1:3) = (/ emnetcdf_ids(1), emnetcdf_ids(2), emnetcdf_ids(3) /)
    attstring(1:9) = "densities"
    call check_fileionetcdf( nf90_def_var(ncid, attstring(1:9), NF90_FLOAT, dims(1:3), emnetcdf_ids(11)) )
    attstring(1:19) = "                   "
    attstring(1:5) = "units"
    attstring(6:11) = "g/cm^3"
    call check_fileionetcdf( nf90_put_att(ncid, emnetcdf_ids(11), attstring(1:5), attstring(6:11)) )
    call check_fileionetcdf( nf90_enddef(ncid) )

!   write spacing data
    call check_fileionetcdf( nf90_put_var(ncid, emnetcdf_ids(5), thatgrid%deltas(1:3)) )
  end if
!
  call check_fileionetcdf( nf90_put_var(ncid, emnetcdf_ids(11), thatgrid%dens) )
!
  call check_fileionetcdf( nf90_close(ncid) )
!
end
!
#else
!
subroutine prt_emmap(thatgrid)
!
  use iounit
  use ems
!
  implicit none
!
  type(t_emgrid) thatgrid
!
  write(ilog,*) 'Fatal. Alternative format writing not yet implemented for density maps. This is an omission bug.'
  call fexit()
!
end
!
#endif
